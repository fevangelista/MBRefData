"""Coordinate whole-geometry ownership and four-thread HBCI sector workers on the Sparks."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
import time

from grid import target, uncertainty_multiplier

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent
STATE = ROOT / ".cluster"
STATE.mkdir(exist_ok=True)
PRIMARY, SECONDARY = "spark-66ea", "spark-6ad0"
HOST = socket.gethostname().split(".")[0]
CONFIG = json.loads((BASE / "all-electron-FCI/input.json").read_text())
GRID = [f"{r:.2f}" for r in CONFIG["distances_angstrom"]]
EXTRA = [f"{r:.2f}" for r in CONFIG["additional_hbci_distances_angstrom"]]
SECTORS = [s["label"] for s in CONFIG["states"]]
PYTHON = Path.home() / "miniconda3/envs/forte2/bin/python"


def owners():
    file = STATE / "owners.json"
    return json.loads(file.read_text()) if file.exists() else {}


def save(data):
    tmp = STATE / "owners.tmp"
    tmp.write_text(json.dumps(data, indent=2) + "\n")
    tmp.replace(STATE / "owners.json")


def done(r, sector):
    fit = ROOT / "full" / f"R-{r}" / f"{sector}-fit.json"
    if not fit.exists():
        return False
    data = json.loads(fit.read_text())
    return (data["qualified"] and data["target_Eh"] == target
            and data["uncertainty_multiplier"] == uncertainty_multiplier)


def processes():
    active = set()
    endpoints = set()
    legacy = set()
    for proc in Path("/proc").glob("[0-9]*/cmdline"):
        try:
            argv = proc.read_bytes().decode(errors="replace").strip("\0").split("\0")
            if not argv or "python" not in Path(argv[0]).name:
                continue
            for i, arg in enumerate(argv):
                name = Path(arg).name
                tail = argv[i + 1:]
                if name == "threshold.py" and tail[:1] == ["full"]:
                    if "--buffer-roots" in tail:
                        endpoints.add(tail[2])
                    else:
                        active.add((f"{float(tail[1]):.2f}", tail[2]))
                if name == "cluster_queue.py" and tail[:1] == ["sector"]:
                    active.add((tail[1], tail[2]))
                if name == "check_endpoint.py":
                    endpoints.add(SECTORS[int("--triplet" in tail)])
                if name == "strict-grid.py":
                    legacy.update(tail)
        except (OSError, ProcessLookupError):
            continue
    return active, endpoints, legacy


def files(r):
    point = ROOT / "full" / f"R-{r}"
    raw, fits = [], []
    for s in SECTORS:
        for path in sorted(point.glob(s + "-*.json")):
            if path.name == s + "-fit.json":
                fits.append(str(path.relative_to(BASE)))
                continue
            log = path.with_suffix(".out")
            if log.is_file():
                with log.open("rb") as stream:
                    stream.seek(max(0, log.stat().st_size - 1000))
                    if b"Finished " + path.name.encode() in stream.read():
                        raw.extend(str(p.relative_to(BASE)) for p in (log, path))
    return raw, fits


def integral(r):
    return BASE / "all-electron-FCI" / f"R-{r}" / "hamiltonian.npz"


def step(action, args):
    owned = owners()
    active, endpoints, legacy = processes()
    # Only the primary writes authoritative reservations. The secondary activates
    # a copy after inputs and the immutable threshold history have been verified.
    if HOST == PRIMARY and "original" in legacy:
        raise RuntimeError("Original static dispatcher still present")
    if HOST == PRIMARY:
        for r, _ in active:
            if r in GRID:
                if owned.get(r, PRIMARY) != PRIMARY:
                    raise RuntimeError(f"Two hosts claim {r}")
                owned[r] = PRIMARY
        if action == "tick" and owned != owners():
            save(owned)
    if action == "reserve":
        assert HOST == PRIMARY
        # Idempotent across a lost response. Deliver all previous reservations first.
        ready = set(args)
        outstanding = [r for r in GRID if owned.get(r) == SECONDARY and r not in ready]
        candidates = [r for r in reversed(GRID) if r not in owned and any(not done(r,s) for s in SECTORS)]
        r = next(iter(outstanding or candidates), None)
        if r:
            owned[r] = SECONDARY
        save(owned)
        return r
    if action == "seed":
        r, = args
        assert HOST == PRIMARY and owned[r] == SECONDARY
        raw, fits = files(r)
        geom = integral(r).parent
        inputs = [str((geom / name).relative_to(BASE)) for name in
                  ("hamiltonian.npz","geometry.xyz","prepare.out","results.json","README.md")]
        return dict(paths=inputs+raw+fits, sha=hashlib.sha256(integral(r).read_bytes()).hexdigest())
    if action == "ready":
        r, sha = args
        assert HOST == SECONDARY and r in GRID
        if hashlib.sha256(integral(r).read_bytes()).hexdigest() != sha:
            raise ValueError(f"Hamiltonian checksum mismatch at {r}")
        raw, fits = files(r)
        available = {Path(p).name for p in raw}
        for p in fits:
            if not set(json.loads((BASE/p).read_text())["threshold_files"]) <= available:
                raise ValueError(f"Incomplete seed history at {p}")
        for p in raw:
            if p.endswith(".json") and json.loads((BASE/p).read_text())["sha256"] != sha:
                raise ValueError(f"Threshold Hamiltonian mismatch at {p}")
        owned[r] = SECONDARY
        save(owned)
        return dict(distance=r, threshold_pairs=len(raw)//2, sha=sha)
    if action == "export":
        assert HOST == SECONDARY
        raw, fits = [], []
        for r in owned:
            rr, ff = files(r)
            raw.extend(rr)
            fits.extend(ff)
        return dict(raw=raw, fits=fits)
    assert action in ("tick","status")
    dedicated = {(r,s) for r in EXTRA for s in SECTORS if not done(r,s)} if (
        HOST == SECONDARY and "extension" in legacy) else set()
    slots = (16 if HOST == PRIMARY else 20) - 2*len(endpoints) - 4*len(active | dedicated)
    launched = []
    # Prefer to finish the other sector at a geometry we already own.
    eligible = [r for r in GRID if owned.get(r) == HOST]
    if HOST == PRIMARY:
        eligible += [r for r in GRID if r not in owned]
    pending = [(r,s) for r in eligible for s in SECTORS if (r,s) not in active and not done(r,s)]
    for r,s in pending:
        if action != "tick" or slots < 4:
            break
        if (STATE / f"{r}-{s}.failed").exists():
            continue
        owned[r] = HOST
        save(owned)
        with (STATE / f"{r}-{s}.out").open("a") as log:
            subprocess.Popen([str(PYTHON), "-u", str(Path(__file__).resolve()), "sector", r,s],
                             env=dict(os.environ, OMP_NUM_THREADS="4", OPENBLAS_NUM_THREADS="1"),
                             cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                             start_new_session=True)
        slots -= 4
        launched.append([r,s])
    return dict(host=HOST, owners=owned, active=sorted(active), endpoint_sectors=sorted(endpoints),
                free_slots=max(0,slots//4), launched=launched,
                failures=sorted(p.name for p in STATE.glob("*.failed")),
                complete=all(done(r,s) for r in GRID for s in SECTORS) if HOST == PRIMARY else None)


def locked(action, args):
    with (STATE / "queue.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        return step(action,args)


if __name__ == "__main__":
    action, *args = sys.argv[1:]
    if HOST not in (PRIMARY,SECONDARY):
        raise RuntimeError("Run on a campaign Spark")
    if action == "drive":
        while True:
            result = locked("tick", [])
            if result["launched"]:
                print(time.ctime(),result["launched"],flush=True)
            if result["complete"]:
                break
            time.sleep(15)
    elif action == "sector":
        r,s = args
        try:
            from campaign import block
            protocol=json.loads((ROOT/"calibration/results.json").read_text())
            assert protocol["target_Eh"] <= target and protocol["calibrated_fit_multiplier"] == uncertainty_multiplier
            block((float(r),s),"full",protocol["minimum_refinement_to_epsilon"],uncertainty_multiplier)
        except Exception:
            (STATE / f"{r}-{s}.failed").touch()
            raise
    else:
        print(json.dumps(locked(action,args)))
