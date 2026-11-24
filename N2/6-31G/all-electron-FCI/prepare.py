"""Prepare an all-electron Hamiltonian and verify it against RHF for one grid point."""
import json
from pathlib import Path
import os
import sys

import forte2
import numpy as np
from forte2.jkbuilder.mointegrals import RestrictedMOIntegrals

here = Path(__file__).resolve().parent
os.chdir(here)  # Forte2 normalizes basis names to lowercase, including any path prefix.
options = json.loads((here / "input.json").read_text())
distance = float(sys.argv[1])
additional = options.get("additional_hbci_distances_angstrom", [])
if distance not in options["distances_angstrom"] + additional:
    raise ValueError("Distance outside the requested grid")
point = here / f"R-{distance:.2f}"
if (point / "results.json").exists():
    previous = json.loads((point / "results.json").read_text())
    if previous.get("fci_states"):
        raise RuntimeError("An FCI result exists. Regenerate in a separate checkout.")
lines = (point / "geometry.xyz").read_text().splitlines()
system = forte2.System(xyz="\n".join(lines[2:]), unit="angstrom",
                       basis_set=options["basis_file"], symmetry=True,
                       cholesky_tei=True, cholesky_tol=options["cholesky_tol"])
rhf = forte2.RHF(charge=0, e_tol=options["rhf_energy_tol"], d_tol=options["rhf_density_tol"],
                 maxiter=200, die_if_not_converged=True)(system)
rhf.run()
norb = options["spatial_orbitals"]
assert system.nbf == norb and rhf.nel == 14 and rhf.converged
assert system.point_group.upper() == "D2H"
C = rhf.mos.C[0]
irrep = np.asarray(rhf.mos.irrep_indices[0], dtype=int)
labels = list(rhf.mos.irrep_labels[0])
assert all(("ag", "b1g", "b2g", "b3g", "au", "b1u", "b2u", "b3u")[i] == str(l).lower()
           for i, l in zip(irrep, labels))
ints = RestrictedMOIntegrals(system, C, range(norb))
nocc = rhf.na
vocc = ints.V[:nocc, :nocc, :nocc, :nocc]
energy_check = (ints.E + 2 * np.trace(ints.H[:nocc, :nocc])
                + 2 * np.einsum("ijij->", vocc) - np.einsum("ijji->", vocc))
if abs(energy_check - rhf.E) > 1e-9:
    raise ValueError("Integral reconstruction of RHF failed")
temp = point / "hamiltonian.tmp.npz"
np.savez_compressed(temp, scalar_Eh=ints.E, h=ints.H, V=ints.V, C=C, orbital_irreps=irrep)
temp.replace(point / "hamiltonian.npz")
record = {"distance_angstrom":distance,"rhf_total_energy_Eh":float(rhf.E),
          "nuclear_repulsion_Eh":float(ints.E),"rhf_reconstruction_error_Eh":float(energy_check-rhf.E),
          "coordinate_file":"geometry.xyz","integral_file":"hamiltonian.npz","orbitals":norb,
          "electrons":rhf.nel,"cholesky_tolerance":options["cholesky_tol"],
          "fci_states":[],"fci_status":"not_requested" if distance in additional else "pending"}
temp = point / "results.json.tmp"
temp.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
temp.replace(point / "results.json")
print(f"Verified Hamiltonian, RHF reconstruction error {energy_check-rhf.E:.3e} Eh.")
