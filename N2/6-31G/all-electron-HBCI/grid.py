"""Read the HBCI grid and common all-electron Hamiltonian options."""
import json
from pathlib import Path

input_file = Path(__file__).resolve().parent.parent / "all-electron-FCI" / "input.json"
settings = json.loads(input_file.read_text())
settings["distances_angstrom"] = sorted(
    settings["distances_angstrom"] + settings.get("additional_hbci_distances_angstrom", [])
)
