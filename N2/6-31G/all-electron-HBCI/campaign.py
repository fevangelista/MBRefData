"""Run a halving sequence, validate it on exact controls and apply it to the curve."""
import concurrent.futures
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np

from grid import settings, target, screening_target, uncertainty_multiplier

root = Path(__file__).resolve().parent
parent = root.parent / "all-electron-FCI"
sectors = [s["label"] for s in settings["states"]]
anchors = [1.10, 2.10]
epsilons = [0.002 * 0.5**i for i in range(23)]

def write(file, obj):
    tmp = file.with_suffix(f".{os.getpid()}.tmp")
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False)+"\n")
    tmp.replace(file)

def solve(space, r, sector, epsilon, strict=False, buffer_roots=None):
    point = root / space / f"R-{r:.2f}"
    if buffer_roots:
        point = root/"diagnostics"/point.name if space=="full" else point/"buffered"
    point.mkdir(parents=True, exist_ok=True)
    stem = f"{sector}-{epsilon:.8g}"
    base = point / (stem+".json")
    tighter = point / (stem+"-pt2.json")

    def read(file):
        data = json.loads(file.read_text())
        data["_source"] = file.name
        return data

    def run(file, flags):
        cmd = [sys.executable, "-u", str(root/"threshold.py"), space, str(r), sector, str(epsilon), *flags]
        if buffer_roots:
            cmd.extend(["--buffer-roots",str(buffer_roots)])
        with file.with_suffix(".out").open("a") as out:
            subprocess.run(cmd, cwd=root, stdout=out, stderr=subprocess.STDOUT, check=True)

    if tighter.exists():
        return read(tighter)
    if not base.exists():
        run(base, [])
    record = read(base)
    if strict and record["pt2_screening_change_Eh"] > screening_target:
        run(tighter, ["--tight-pt2-copy"])
        return read(tighter)
    return record

def match_bounds(a, b):
    common = a.keys() & b.keys()
    dot = sum(a[k]*b[k] for k in common)
    rem_a = max(0, 1-sum(a[k]**2 for k in common))
    rem_b = max(0, 1-sum(b[k]**2 for k in common))
    bound = (rem_a*rem_b)**0.5
    return max(0,abs(dot)-bound), min(1,abs(dot)+bound)

def fits(records, sector, level):
    rows = records[-5:]
    e = np.array([q["E_var_Eh"] for q in rows])
    p = np.array([q["pt2_history"][-1]["pt2_Eh"] for q in rows])
    s2 = 0 if sector==sectors[0] else 2
    estimates = []
    for i in range(e.shape[1]):
        def intercept(x, y):
            if np.ptp(x) < 1e-14:
                return float(y[-1]+x[-1])
            return float(np.polynomial.polynomial.polyfit(x,y,1)[0])
        now = intercept(p[-4:,i], e[-4:,i])
        prior = intercept(p[:-1,i], e[:-1,i])
        three = intercept(p[-3:,i], e[-3:,i])
        screening = max(q["pt2_screening_change_Eh"] for q in rows)
        stable = max(abs(now-prior),abs(now-three))
        err = level*(stable + screening)
        matched = True
        for a,b in zip(rows,rows[1:]):
            c = a["dominant_coefficients"][i]
            cs = b["dominant_coefficients"]
            low = match_bounds(c,cs[i])[0]
            competing=max(match_bounds(c,other)[1] for j,other in enumerate(cs) if j!=i)
            matched &= low>0.5 and low>competing
        estimates.append(dict(rank=i+1, extrapolated_total_Eh=now,
            last_variational_Eh=float(e[-1,i]),last_PT2_Eh=float(p[-1,i]),
            fit_stability_Eh=stable, pt2_screening_change_Eh=screening, estimated_uncertainty_Eh=err,
            root_match=bool(matched), spin_squared=rows[-1]["spin_squared"][i],
            diagnostics_pass=bool(matched and abs(rows[-1]["spin_squared"][i]-s2)<1e-5
                                  and screening<=screening_target and err<=target)))
    return estimates

def block(job, space, minimum, factor):
    r,sector=job
    point=root/space/f"R-{r:.2f}"
    records=[]
    truth=None
    if space=="calibration":
        exact_file=point/"exact-tight.json"
        if not exact_file.exists():
            exact_file=point/"exact.json"
        truth=json.loads(exact_file.read_text())[sector]
    for eps in epsilons:
        q=solve(space,r,sector,eps)
        records.append(q)
        print(space,r,sector,eps,q["ndet"],flush=True)
        if len(records)<5:
            continue
        estimates=fits(records,sector,factor)
        # Refresh a loose historical PT2 sum once the variational fit is near the target.
        if eps<=minimum and (space!="full" or max(s["fit_stability_Eh"] for s in estimates)<=5*target):
            for j in range(len(records)-5,len(records)):
                if records[j]["pt2_screening_change_Eh"]>screening_target:
                    records[j]=solve(space,r,sector,records[j]["epsilon_var"],strict=True)
            estimates=fits(records,sector,factor)
        if truth:
            for s,exact,residual in zip(estimates,truth["E_FCI_Eh"],truth["residual_Eh"]):
                s["exact_total_Eh"]=exact
                s["error_vs_FCI_Eh"]=s["extrapolated_total_Eh"]-exact
                s["exact_residual_Eh"]=residual
        ok=all(s["diagnostics_pass"] for s in estimates) and eps<=minimum
        if truth:
            ok &= all(abs(s["error_vs_FCI_Eh"])+s["exact_residual_Eh"]<=target for s in estimates)
        out=dict(space=space, distance_angstrom=r, sector=sector, qualified=bool(ok),
                 epsilon_var=eps, target_Eh=target, pt2_screening_target_Eh=screening_target,
                 uncertainty_multiplier=factor, levels=[v["epsilon_var"] for v in records],
                 threshold_files=[v["_source"] for v in records], fits=estimates)
        write(point/(sector+"-fit.json"),out)
        if ok: return out
    raise RuntimeError(f"Target unresolved for {job}. See archived history.")

def phase(space, jobs, minimum, factor, workers):
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        tasks=[pool.submit(block,j,space,minimum,factor) for j in jobs]
        return [f.result() for f in concurrent.futures.as_completed(tasks)]

def calibrate(workers):
    controls=phase("calibration",[(r,s) for r in anchors for s in sectors],0.002,uncertainty_multiplier,workers)
    protocol=dict(minimum_refinement_to_epsilon=min(b["epsilon_var"] for b in controls),
                  calibrated_fit_multiplier=uncertainty_multiplier,target_Eh=target,
                  pt2_screening_target_Eh=screening_target,controls=controls)
    write(root/"calibration"/"results.json",protocol)
    return protocol

if __name__=="__main__":
    workers=int(os.environ.get("HBCI_WORKERS","2"))
    protocol=calibrate(workers)
    distances=anchors+[r for r in settings["distances_angstrom"] if r not in anchors]
    phase("full",[(r,s) for r in distances for s in sectors],
          protocol["minimum_refinement_to_epsilon"],uncertainty_multiplier,workers)
