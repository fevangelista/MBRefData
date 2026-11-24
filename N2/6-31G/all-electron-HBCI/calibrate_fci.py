"""Build the [3s1p] nitrogen subset and diagonalize both spin sectors exactly."""
import json
import os
import sys
from pathlib import Path
import numpy as np
from forte2 import CISolver, RHF, State, System
from forte2.base_classes import CIParams, DavidsonLiuParams
from forte2.jkbuilder.mointegrals import RestrictedMOIntegrals

root = Path(__file__).resolve().parent
parent = root.parent / "all-electron-FCI"
options = json.loads((parent / "input.json").read_text())
r = float(sys.argv[1])
assert r in (1.1, 2.1)
point = root / "calibration" / f"R-{r:.2f}"
point.mkdir(exist_ok=True, parents=True)
os.chdir(point)
basis = json.loads((parent / options["basis_file"]).read_text())
shell = basis["elements"]["7"]["electron_shells"][-1]
assert shell["angular_momentum"] == [0, 1]
shell["angular_momentum"] = [0]
shell["coefficients"] = shell["coefficients"][:1]
Path("basis.json").write_text(json.dumps(basis, indent=2) + "\n")
xyz = (parent / point.name / "geometry.xyz").read_text().splitlines()[2:]
system = System(xyz="\n".join(xyz), unit="angstrom", basis_set="basis.json", symmetry=True,
                cholesky_tei=True, cholesky_tol=options["cholesky_tol"])
hf = RHF(charge=0, e_tol=1e-12, d_tol=1e-9, maxiter=200)(system)
hf.run()
assert system.nbf == 12 and hf.nel == 14 and hf.converged
C = hf.mos.C[0]
ints = RestrictedMOIntegrals(system, C, range(12))
np.savez_compressed("hamiltonian.npz", scalar_Eh=ints.E, h=ints.H, V=ints.V,
                    C=C, orbital_irreps=hf.mos.irrep_indices[0])
results = {}
for s in options["states"]:
    state = State(nel=14, multiplicity=s["multiplicity"], ms=s["ms"],
                  symmetry=s["irrep_index"], symmetry_label=s["irrep"])
    solver = CISolver(states=state, nroots=2, active_orbitals=12, core_orbitals=0,
                      die_if_not_converged=True, ci_params=CIParams(ci_algorithm="hz"),
                      davidson_liu_params=DavidsonLiuParams(e_tol=1e-12,r_tol=5e-8,
                          basis_per_root=32, collapse_per_root=6, maxiter=200))(hf)
    solver.run()
    worker = solver.sub_solvers[0]
    assert worker.eigensolver.converged
    sigma = np.empty_like(worker.evecs)
    worker._make_sigma_builder()(worker.evecs, sigma)
    residuals = np.linalg.norm(sigma-worker.evecs*worker.evals, axis=0)
    if max(residuals) > 5e-8:
        raise RuntimeError(f"FCI residuals {residuals}")
    results[s["label"]] = dict(E_FCI_Eh=worker.evals.tolist(), residual_Eh=residuals.tolist())
    print(s["label"],results[s["label"]],flush=True)
Path("exact.json.tmp").write_text(json.dumps(results, indent=2) + "\n")
Path("exact.json.tmp").replace("exact.json")
