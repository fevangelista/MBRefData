"""Refine the four-root endpoint singlet calculation and check its spectral ordering."""
import argparse
import json
from pathlib import Path
import sys

root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root))
from campaign import epsilons, fits, screening_target, solve, target, write

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--triplet",action="store_true")
triplet=parser.parse_args().triplet
sector="triplet-B1u" if triplet else "singlet-Ag"
report_file="triplet-rank-check.json" if triplet else "results.json"
point=root/"diagnostics/R-4.00"
protocol=json.loads((root/"calibration/results.json").read_text())
assert protocol["target_Eh"]==target
minimum=protocol["minimum_refinement_to_epsilon"]
factor=protocol["calibrated_fit_multiplier"]
measurements=[]
for eps in epsilons:
    measurements.append(solve("full",4.0,sector,eps,buffer_roots=4))
    if len(measurements)<5:
        continue
    estimates=fits(measurements,sector,factor)
    if eps<=minimum and estimates[1]["fit_stability_Eh"]<=5*target:
        for i in range(len(measurements)-5,len(measurements)):
            if measurements[i]["pt2_screening_change_Eh"]>screening_target:
                measurements[i]=solve("full",4.0,sector,measurements[i]["epsilon_var"],strict=True,buffer_roots=4)
        estimates=fits(measurements,sector,factor)
    first,second=estimates[:2]
    lower=lambda s: s["extrapolated_total_Eh"]-s["estimated_uncertainty_Eh"]
    upper=lambda s: s["extrapolated_total_Eh"]+s["estimated_uncertainty_Eh"]
    separated=upper(first)<lower(second) and upper(second)<min(lower(s) for s in estimates[2:])
    confirm=bool(eps<=minimum and second["diagnostics_pass"] and separated and
                 all(s["root_match"] and abs(s["spin_squared"]-(2 if triplet else 0))<1e-5 for s in estimates))
    report=dict(sector=sector,rank=2,rank_confirmed=confirm,target_Eh=target,
                epsilon_var=eps,threshold_files=[q["_source"] for q in measurements],
                spectral_intervals_separated=separated,uncertainty_multiplier=factor,
                fits=estimates,coarse_comparison="../../validation/R-4.00/first-strict-check.json" if triplet else "initial-comparison.json")
    write(point/report_file,report)
    print(eps,[round(s["estimated_uncertainty_Eh"]*1e6,4) for s in estimates],separated,confirm,flush=True)
    if confirm:
        break
else:
    raise RuntimeError("The rank or target precision of the second endpoint singlet remains unresolved")
