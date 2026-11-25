"""Refine an endpoint four-root HBCI calculation and check its spectral ordering."""
import argparse
import json
from pathlib import Path
import sys

root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root))
from campaign import control_target, epsilons, fits, screening_target, solve, target, write

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--triplet",action="store_true")
parser.add_argument("--validation",action="store_true",help="Use the 12-orbital endpoint control")
args=parser.parse_args()
triplet=args.triplet
sector="triplet-B1u" if triplet else "singlet-Ag"
report_file="triplet-rank-check.json" if triplet else "results.json"
space="validation" if args.validation else "full"
point=root/"validation/R-4.00/buffered" if args.validation else root/"diagnostics/R-4.00"
point.mkdir(exist_ok=True,parents=True)
protocol=json.loads((root/"calibration/results.json").read_text())
assert protocol["target_Eh"]==control_target
limit=control_target if args.validation else target
minimum=protocol["minimum_refinement_to_epsilon"]
factor=protocol["calibrated_fit_multiplier"]
measurements=[]
for eps in epsilons:
    measurements.append(solve(space,4.0,sector,eps,buffer_roots=4))
    if len(measurements)<5:
        continue
    estimates=fits(measurements,sector,factor,limit)
    if eps<=minimum and estimates[1]["fit_stability_Eh"]<=5*limit:
        for i in range(len(measurements)-5,len(measurements)):
            if measurements[i]["pt2_screening_change_Eh"]>screening_target:
                measurements[i]=solve(space,4.0,sector,measurements[i]["epsilon_var"],strict=True,buffer_roots=4)
        estimates=fits(measurements,sector,factor,limit)
    first,second=estimates[:2]
    lower=lambda s: s["extrapolated_total_Eh"]-s["estimated_uncertainty_Eh"]
    upper=lambda s: s["extrapolated_total_Eh"]+s["estimated_uncertainty_Eh"]
    separated=upper(first)<lower(second) and upper(second)<min(lower(s) for s in estimates[2:])
    confirm=bool(eps<=minimum and second["diagnostics_pass"] and separated and first["root_match"] and
                 all(abs(s["spin_squared"]-(2 if triplet else 0))<1e-5 for s in estimates))
    report=dict(sector=sector,rank=2,rank_confirmed=confirm,target_Eh=limit,
                epsilon_var=eps,threshold_files=[q["_source"] for q in measurements],
                spectral_intervals_separated=separated,uncertainty_multiplier=factor,
                fits=estimates,coarse_comparison=("../first-strict-check.json" if args.validation else "../../validation/R-4.00/first-strict-check.json") if triplet else "initial-comparison.json")
    write(point/report_file,report)
    if args.validation and triplet:
        write(point/"results.json",report)
    print(eps,[round(s["estimated_uncertainty_Eh"]*1e6,4) for s in estimates],separated,confirm,flush=True)
    if confirm:
        break
else:
    raise RuntimeError("The rank or target precision of the second endpoint singlet remains unresolved")
