"""Rebuild the curve summary from individual, verified calculation outputs."""
import json
from pathlib import Path

root = Path(__file__).resolve().parent
settings = json.loads((root / "input.json").read_text())
records = []
for r in settings["distances_angstrom"]:
    point = root / f"R-{r:.2f}"
    record = json.loads((point / "results.json").read_text())
    assert record["distance_angstrom"] == r
    records.append(dict(distance_angstrom=r, file=f"{point.name}/results.json",
                        fci_status=record["fci_status"], fci_states=record["fci_states"]))
result = dict(molecule="N2", method="FCI", basis=settings["basis"], electrons=14,
              frozen_core=0, spatial_orbitals=18, cholesky_tolerance=settings["cholesky_tol"],
              energy_unit="hartree", energy_includes_nuclei=True,
              n_geometries=len(records), requested_states_per_geometry=settings["states"],
              converged_state_count=sum(len(x["fci_states"]) for x in records), points=records)
temp = root / "results.json.tmp"
temp.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
temp.replace(root / "results.json")
