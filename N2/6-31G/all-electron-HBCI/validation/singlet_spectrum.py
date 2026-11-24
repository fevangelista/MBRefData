"""Resolve nearby singlet or triplet levels in the 4 Å, 12-orbital validation Hamiltonian."""
import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from forte2 import MOSpace, State
from forte2.base_classes import CIParams, DavidsonLiuParams
from forte2.ci.ci import _CISingleStateSolver
from forte2.helpers import logger

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--triplet",action="store_true")
triplet=parser.parse_args().triplet
label="triplet-B1u" if triplet else "singlet-Ag"
state=State(nel=14,multiplicity=3 if triplet else 1,ms=int(triplet),
            symmetry=5 if triplet else 0,symmetry_label="B1u" if triplet else "Ag")
point = Path(__file__).resolve().parent / "R-4.00"
archive = point / "hamiltonian.npz"
with np.load(archive) as data:
    integrals = SimpleNamespace(E=float(data["scalar_Eh"]), H=data["h"], V=data["V"])
    irreps = data["orbital_irreps"].tolist()
assert len(irreps) == 12
solver = _CISingleStateSolver(
    mo_space=MOSpace(nmo=12, active_orbitals=list(range(12))),
    state=state,
    ints=integrals, nroot=4, active_orbsym=[irreps], ci_params=CIParams(ci_algorithm="hz"),
    davidson_liu_params=DavidsonLiuParams(e_tol=1e-13,r_tol=1e-9,maxiter=400,
                                         basis_per_root=32,collapse_per_root=6),
    log_level=logger.VERBOSITY_INFO1, die_if_not_converged=True)
solver.run()
sigma = np.empty_like(solver.evecs)
solver._make_sigma_builder()(solver.evecs,sigma)
residuals = np.linalg.norm(sigma - solver.evecs*solver.evals,axis=0)
assert solver.eigensolver.converged and max(residuals) <= 1e-9
output = dict(hamiltonian_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
              sector=label,E_FCI_Eh=solver.evals.tolist(),residual_Eh=residuals.tolist(),
              adjacent_gaps_Eh=np.diff(solver.evals).tolist())
(point/("triplet-spectrum.json" if triplet else "singlet-spectrum.json")).write_text(json.dumps(output,indent=2)+"\n")
print(output)
