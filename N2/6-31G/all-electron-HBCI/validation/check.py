"""Apply the calibrated production rule before comparing to an independent exact calculation."""
import json
from pathlib import Path
import sys

root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root))
from campaign import phase, sectors, target, write

distance = 4.0
point = root / "validation" / f"R-{distance:.2f}"
protocol = json.loads((root / "calibration/results.json").read_text())
if protocol["target_Eh"] != target:
    raise ValueError("Run the 0.1 µEh calibration first")
blocks = phase("validation", [(distance, label) for label in sectors],
               protocol["minimum_refinement_to_epsilon"], protocol["calibrated_fit_multiplier"], 2)
exact_file = point / "exact-tight.json"
exact = json.loads((exact_file if exact_file.exists() else point/"exact.json").read_text())
states = []
for block in sorted(blocks, key=lambda block: block["sector"]):
    sector = block["sector"]
    for fit in block["fits"]:
        i = fit["rank"] - 1
        energy = exact[sector]["E_FCI_Eh"][i]
        deviation = fit["extrapolated_total_Eh"] - energy
        states.append(dict(sector=sector, rank=i+1, epsilon_var=block["epsilon_var"],
                           E_FCI_Eh=energy, FCI_residual_Eh=exact[sector]["residual_Eh"][i],
                           extrapolated_total_Eh=fit["extrapolated_total_Eh"],
                           estimated_uncertainty_Eh=fit["estimated_uncertainty_Eh"],
                           error_vs_FCI_Eh=deviation, within_target=abs(deviation)+exact[sector]["residual_Eh"][i]<=protocol["target_Eh"]))
result = dict(distance_angstrom=distance, target_Eh=protocol["target_Eh"], electrons=14,
              spatial_orbitals=12, states=states, all_within_target=all(s["within_target"] for s in states))
write(point / "results.json", result)
write(point.parent / "results.json", result)
