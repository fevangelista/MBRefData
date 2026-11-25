"""Read the HBCI grid and common all-electron Hamiltonian options."""
import json
from pathlib import Path

input_file = Path(__file__).resolve().parent.parent / "all-electron-FCI" / "input.json"
settings = json.loads(input_file.read_text())
settings["distances_angstrom"] = sorted(
    settings["distances_angstrom"] + settings.get("additional_hbci_distances_angstrom", [])
)

# The original 4 Å control underestimated a 0.3173 µEh error as 0.0932 µEh.
# Round that observed ratio upward for the tighter campaign.
target = 1e-6
control_target = 1e-7
screening_target = 2e-8
uncertainty_multiplier = 4.0
