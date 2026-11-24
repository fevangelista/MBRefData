"""Assemble per-distance and curve results from completed threshold runs and fit files."""
import json
from pathlib import Path

from grid import settings, target, screening_target, uncertainty_multiplier

root = Path(__file__).resolve().parent
parent = root.parent / "all-electron-FCI"
sectors = [s["label"] for s in settings["states"]]

def load(path):
    return json.loads(path.read_text())

def write(path, value):
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")
    tmp.replace(path)

for space in ["calibration","full"]:
    base = root / space
    base.mkdir(exist_ok=True)
    distances = [1.1,2.1] if space=="calibration" else settings["distances_angstrom"]
    curve = []
    for r in distances:
        point = base/f"R-{r:.2f}"
        point.mkdir(exist_ok=True)
        states, trials = [], {}
        rank_checks = {}
        if space=="full":
            diagnostic=root/"diagnostics"/point.name
            for name in ["results.json","triplet-rank-check.json"]:
                file=diagnostic/name
                if file.exists():
                    check=load(file)
                    if "rank_confirmed" in check:
                        rank_checks[(check["sector"],check["rank"])]=(check,f"../../diagnostics/{point.name}/{name}")
        for sector in sectors:
            measurements = {}
            for file in sorted(point.glob(sector+"-*.json"),key=lambda f:f.stem.endswith("-pt2")):
                if file.name.endswith("-fit.json"):
                    continue
                data=load(file)
                data["_source"]=file.name
                measurements[data["epsilon_var"]]=data
            levels=list(measurements.values())
            levels.sort(key=lambda q:-q["epsilon_var"])
            trials[sector] = [dict(epsilon_var=q["epsilon_var"], ndet=q["ndet"],
                 E_var_plus_PT2_Eh=q["E_var_plus_PT2_Eh"],
                 pt2_threshold=q["pt2_history"][-1]["threshold"],file=q["_source"]) for q in levels]
            fit_file = point/(sector+"-fit.json")
            if fit_file.exists():
                fit=load(fit_file)
                if fit["qualified"] and fit.get("target_Eh")==target and fit.get("uncertainty_multiplier")==uncertainty_multiplier:
                    assert len(fit["fits"]) == 2
                    for s in fit["fits"]:
                        assert s["diagnostics_pass"] and s["estimated_uncertainty_Eh"]<=target and s["pt2_screening_change_Eh"]<=screening_target
                        if (sector,s["rank"]) in rank_checks:
                            continue  # The four-root comparison supplies this reference after resolving its rank.
                        states.append(dict(sector=sector,epsilon_var=fit["epsilon_var"],**s))
        for (sector,rank),(check,source) in rank_checks.items():
            if check["rank_confirmed"] and check.get("target_Eh")==target:
                state=check["fits"][rank-1]
                assert state["diagnostics_pass"] and state["estimated_uncertainty_Eh"]<=target
                states.append(dict(sector=sector,epsilon_var=check["epsilon_var"],source=source,**state))
        states.sort(key=lambda state:(state["sector"],state["rank"]))
        status = "qualified" if len(states)==4 else "partial" if states else "in_progress" if any(trials.values()) else "pending"
        result=dict(distance_angstrom=r,space=space,status=status,target_per_root_Eh=target,
                    exact_control=load(point/"exact-tight.json") if (point/"exact-tight.json").exists() else load(point/"exact.json") if (point/"exact.json").exists() else None,
                    qualified_states=states, threshold_measurements=trials)
        if rank_checks:
            result["independent_rank_checks"]=[dict(sector=sector,rank=rank,file=source,
                rank_confirmed=bool(check["rank_confirmed"] and check.get("target_Eh")==target))
                for (sector,rank),(check,source) in rank_checks.items()]
        write(point/"results.json",result)
        curve.append(dict(distance_angstrom=r,status=status,file=f"{point.name}/results.json",
                          qualified_states=states))
        if not (point/"README.md").exists():
            (point/"README.md").write_text(f"# N₂ at {r:.2f} Å\n\nThis {space} calculation uses all 14 electrons. See the [dataset methods](../../README.md) and [results](results.json) for state identities, selection thresholds, and qualification status. Each threshold JSON has the corresponding raw program output in the file of the same stem with extension .out.\n")
    if space=="full":
        summary=dict(method="HBCI with deterministic EN-PT2 and zero-PT2 extrapolation",
                     basis="6-31G",electrons=14,orbitals=18,cholesky_tolerance=1e-12,
                     target_per_root_Eh=target,pt2_screening_target_Eh=screening_target,
                     uncertainty_multiplier=uncertainty_multiplier,
                     qualified_roots=sum(len(p["qualified_states"]) for p in curve),
                     requested_roots=len(distances)*4,points=curve)
        write(base/"results.json",summary)
        write(root/"results.json",dict(summary,points=[dict(p,file="full/"+p["file"]) for p in curve]))
