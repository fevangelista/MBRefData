"""Run the grid serially; completed blocks are skipped when resuming."""
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
settings = json.loads((root / "input.json").read_text())
os.chdir(root)
grid = settings["distances_angstrom"]
order = [1.10, 2.10] + [r for r in grid if r not in (1.10, 2.10)]
# Only one large spin-adapted block fits at a time on the Spark.
with (root / ".grid.lock").open("a") as lock:
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    for r in order:
        point = root / f"R-{r:.2f}"
        for spec in settings["states"]:
            sector = spec["label"]
            data = json.loads((point / "results.json").read_text())
            if any(x["sector"] == sector for x in data["fci_states"]):
                continue
            log_path = point / f"{sector}.out"
            print(f"Starting {point.name} {sector}", flush=True)
            with log_path.open("a") as out:
                subprocess.run([sys.executable, "-u", "solve.py", str(r), sector],
                               stdout=out, stderr=subprocess.STDOUT, check=True)
            print(f"Finished {point.name} {sector}", flush=True)
