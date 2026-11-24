"""Compute one spin- and symmetry-resolved FCI block from an archived Forte2 Hamiltonian."""
import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import time

import numpy as np
from forte2 import MOSpace, State
from forte2.base_classes import CIParams, DavidsonLiuParams
from forte2.ci.ci import _CISingleStateSolver
from forte2.helpers import logger

root = Path(__file__).resolve().parent
settings = json.loads((root / "input.json").read_text())
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("distance", type=float, choices=settings["distances_angstrom"])
parser.add_argument("sector", choices=[s["label"] for s in settings["states"]])
args = parser.parse_args()
point = root / f"R-{args.distance:.2f}"
target = next(s for s in settings["states"] if s["label"] == args.sector)
result_path = point / "results.json"
result = json.loads(result_path.read_text())
if any(s["sector"] == args.sector for s in result["fci_states"]):
    raise SystemExit("Sector already archived")

logger.set_verbosity_level(logger.VERBOSITY_INFO1)
start = time.perf_counter()
archive = point / result["integral_file"]
digest = hashlib.sha256(archive.read_bytes()).hexdigest()
with np.load(archive, allow_pickle=False) as data:
    ints = SimpleNamespace(E=float(data["scalar_Eh"]), H=data["h"], V=data["V"])
    orb_irreps = data["orbital_irreps"].tolist()
norb = settings["spatial_orbitals"]
assert norb == len(orb_irreps) == ints.H.shape[0] == 18
state = State(nel=settings["electrons"], multiplicity=target["multiplicity"],
              ms=target["ms"], symmetry=target["irrep_index"], symmetry_label=target["irrep"])
# This is the same spin-adapted worker that the Forte2 CISolver uses. Passing the
# saved integrals makes the Hamiltonian independent of repeated SCF diagonalizations.
solver = _CISingleStateSolver(
    mo_space=MOSpace(nmo=norb, active_orbitals=list(range(norb))),
    state=state, ints=ints, nroot=target["nroots"], active_orbsym=[orb_irreps],
    ci_params=CIParams(ci_algorithm="hz", ci_builder_memory=2048),
    davidson_liu_params=DavidsonLiuParams(
        basis_per_root=6, collapse_per_root=2, maxiter=200,
        e_tol=settings["ci_energy_tol"], r_tol=settings["ci_residual_tol"]),
    log_level=logger.VERBOSITY_INFO1, die_if_not_converged=True,
)
print(f"FCI {args.distance:.2f} angstrom {args.sector}; Hamiltonian SHA256 {digest}", flush=True)
solver.run()
assert solver.eigensolver.converged
apply_h = solver._make_sigma_builder()
states = []
for index, energy in enumerate(solver.evals):
    ci_vector = solver.evecs[:, index:index+1]
    sigma = np.empty_like(ci_vector)
    apply_h(ci_vector, sigma)
    residual = np.linalg.norm(sigma - energy * ci_vector)
    if not np.isfinite(residual) or residual > settings["ci_residual_tol"]:
        raise RuntimeError(f"Root {index+1} residual is {residual}")
    states.append(dict(sector=args.sector, rank_in_sector=index+1,
                       multiplicity=target["multiplicity"], ms=target["ms"], irrep=target["irrep"],
                       total_energy_Eh=float(energy), residual_norm_Eh=float(residual),
                       n_determinants=int(solver.ndet), n_csfs=int(solver.basis_size),
                       hamiltonian_sha256=digest))
    print(json.dumps(states[-1]), flush=True)
result["fci_states"].extend(states)
expected = sum(s["nroots"] for s in settings["states"])
result["fci_status"] = "complete" if len(result["fci_states"]) == expected else "partial"
temp = point / "results.json.tmp"
temp.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
temp.replace(result_path)
print(f"Verified all {len(states)} roots for this block in {time.perf_counter()-start:.1f} s.")
