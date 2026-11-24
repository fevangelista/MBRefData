"""Tighten reduced-basis exact-control residuals on the archived Hamiltonian."""
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
from campaign import root, settings, write

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("space",choices=["calibration","validation"])
parser.add_argument("distance",type=float)
args=parser.parse_args()
point=root/args.space/f"R-{args.distance:.2f}"
archive=point/"hamiltonian.npz"
digest=hashlib.sha256(archive.read_bytes()).hexdigest()
with np.load(archive) as data:
    ints=SimpleNamespace(E=float(data["scalar_Eh"]), H=data["h"], V=data["V"])
    irreps=data["orbital_irreps"].tolist()
assert len(irreps)==12
result={}
for sector in settings["states"]:
    label=sector["label"]
    spectrum=point/"singlet-spectrum.json"
    if spectrum.exists() and label=="singlet-Ag":
        old=json.loads(spectrum.read_text())
        assert old["hamiltonian_sha256"]==digest and max(old["residual_Eh"])<1e-9
        result[label]={name:old[name][:2] for name in ("E_FCI_Eh","residual_Eh")}
        continue
    solver=_CISingleStateSolver(
        mo_space=MOSpace(nmo=12,active_orbitals=list(range(12))),
        state=State(nel=14,multiplicity=sector["multiplicity"],ms=sector["ms"],
                    symmetry=sector["irrep_index"],symmetry_label=sector["irrep"]),
        ints=ints,nroot=2,active_orbsym=[irreps],ci_params=CIParams(ci_algorithm="hz"),
        davidson_liu_params=DavidsonLiuParams(e_tol=1e-13,r_tol=1e-9,maxiter=400,
                                             basis_per_root=32,collapse_per_root=6),
        log_level=logger.VERBOSITY_INFO1,die_if_not_converged=True)
    solver.run()
    sigma=np.empty_like(solver.evecs)
    solver._make_sigma_builder()(solver.evecs,sigma)
    residual=np.linalg.norm(sigma-solver.evecs*solver.evals,axis=0)
    assert solver.eigensolver.converged and max(residual)<=1e-9
    result[label]=dict(E_FCI_Eh=solver.evals.tolist(),residual_Eh=residual.tolist())
result["hamiltonian_sha256"]=digest
write(point/"exact-tight.json",result)
print(result,flush=True)
