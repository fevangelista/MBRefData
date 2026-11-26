"""Solve one variational threshold and converge deterministic EN-PT2 at fixed wavefunction."""
import argparse
import hashlib
import fcntl
import json
import socket
from math import comb
from pathlib import Path
import time

import numpy as np
from forte2 import MOSpace, State
from forte2.base_classes import DavidsonLiuParams, SelectedCIParams
from forte2.helpers import logger
from forte2.lib.det import Determinant, hilbert_space
from forte2.sci.sci import _SelectedCISingleStateSolver
from types import SimpleNamespace

from grid import settings as options, screening_target

root = Path(__file__).resolve().parent
parent = root.parent / "all-electron-FCI"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("space", choices=["calibration", "validation", "full"])
parser.add_argument("r", type=float, choices=options["distances_angstrom"])
parser.add_argument("sector", choices=[s["label"] for s in options["states"]])
parser.add_argument("epsilon", type=float)
parser.add_argument("--tight-pt2-copy",action="store_true",help="Archive a stricter recomputation alongside the original threshold")
parser.add_argument("--buffer-roots", type=int, choices=[3, 4, 5],
                    help="Run a separate wide-guess diagnostic with extra roots")
args = parser.parse_args()
point = root / args.space / f"R-{args.r:.2f}"
out_dir = point
if args.buffer_roots:
    out_dir = root/"diagnostics"/point.name if args.space=="full" else point/"buffered"
out_dir.mkdir(exist_ok=True, parents=True)
suffix = "-pt2" if args.tight_pt2_copy else ""
file = out_dir / f"{args.sector}-{args.epsilon:.8g}{suffix}.json"
# A coarse grid survey can meet a refinement worker at the same threshold.
lock = file.with_suffix(".lock").open("a")
fcntl.flock(lock, fcntl.LOCK_EX)
if file.exists():
    raise SystemExit(0)
archive = (parent / point.name if args.space == "full" else point) / "hamiltonian.npz"
target = next(s for s in options["states"] if s["label"] == args.sector)
state = State(nel=14, multiplicity=target["multiplicity"], ms=target["ms"],
              symmetry=target["irrep_index"], symmetry_label=target["irrep"])
with np.load(archive) as a:
    irreps = a["orbital_irreps"].astype(int)
    ints = SimpleNamespace(E=float(a["scalar_Eh"]), H=a["h"], V=a["V"])
pair = irreps[:, None] ^ irreps[None, :]
# Forbidden operators have zero matrix elements inside a single target irrep.
# The Frobenius bound also limits the effect if one compares the entire operator.
forbidden_bound = 0.0
for integrals, forbidden, count in [(ints.H, pair != 0, 14),
        (ints.V, pair[:, :, None, None] != pair[None, None, :, :], comb(14, 2))]:
    forbidden_bound += count * np.linalg.norm(integrals[forbidden])
    integrals[forbidden] = 0.0
if forbidden_bound >= 1e-6:
    raise ValueError(f"Inspect symmetry labels: discarded-operator bound {forbidden_bound} Eh")
norb = len(irreps)
# Start from determinants of the requested irrep, including single and double excitations.
# The guess restriction is lifted during heat-bath selection.
seeds = hilbert_space(10, state.na, state.nb, 8, irreps[:10].tolist(), state.symmetry)
nroots = args.buffer_roots or target["nroots"]
if args.buffer_roots:
    ref = Determinant("2"*state.nb+"a"*(state.na-state.nb))
    seeds = hilbert_space(norb,state.na,state.nb,ref,2,8,irreps.tolist(),state.symmetry)
dav = DavidsonLiuParams(e_tol=1e-11, r_tol=5e-8, basis_per_root=18, collapse_per_root=4,
                       maxiter=250, ndets_per_guess=20)
params = SelectedCIParams(guess_dets=seeds, var_threshold=args.epsilon,
                          pt2_threshold=args.epsilon, e_tol=1e-13, maxcycle=120,
                          num_batches_per_thread=2)
logger.set_verbosity_level(logger.VERBOSITY_INFO1)
t = time.monotonic()
solver = _SelectedCISingleStateSolver(mo_space=MOSpace(nmo=norb, active_orbitals=list(range(norb))),
    state=state, ints=ints, nroot=nroots, active_orbsym=[irreps.tolist()], sci_params=params,
    davidson_liu_params=dav, die_if_not_converged=True, log_level=logger.VERBOSITY_INFO1)
solver.run()
assert solver.eigensolver.converged and solver.sci_helper.num_new_dets_var() == 0
energies = np.asarray(solver.e_var)
assert np.all(np.diff(energies)>=0)
helper = solver.sci_helper
assert helper.ndets()==solver.evecs.shape[0]
# No promotions during these calls: PT2 thresholds are compared in the identical space.
pt2_history = []
for pt_threshold in [1e-5, 2e-6, 4e-7, 8e-8, 1.6e-8, 3.2e-9, 0.0]:
    helper.select_hbci(float("inf"), pt_threshold)
    assert helper.ndets()==len(solver.dets)
    correction = np.array(helper.ept2_pt())
    assert np.isfinite(correction).all()
    pt2_history.append(dict(threshold=pt_threshold, pt2_Eh=correction.tolist()))
    if len(pt2_history)>=2:
        change = float(np.max(abs(correction-pt2_history[-2]["pt2_Eh"])))
        if change <= screening_target:
            break
if change > screening_target and pt_threshold != 0.0:
    raise RuntimeError("PT2 screening did not converge")
if pt_threshold == 0.0:
    change = 0.0  # The external-space sum is unscreened.
dominant = []
for c in solver.evecs.T:
    idx = np.argsort(abs(c))[-128:]
    dominant.append({solver.dets[int(i)].str(norb):float(c[i]) for i in idx})
record = dict(execution_host=socket.gethostname(), space=args.space, distance_angstrom=args.r, sector=args.sector,
              forbidden_operator_bound_Eh=float(forbidden_bound), nroots=nroots, wide_guess=bool(args.buffer_roots), epsilon_var=args.epsilon, ndet=helper.ndets(), E_var_Eh=energies.tolist(),
              E_var_plus_PT2_Eh=(energies+correction).tolist(), pt2_history=pt2_history,
              pt2_screening_change_Eh=change, pt2_screening_target_Eh=screening_target, spin_squared=list(helper.compute_spin2()),
              sha256=hashlib.sha256(archive.read_bytes()).hexdigest(), seconds=time.monotonic()-t,
              dominant_coefficients=dominant)
temp=file.with_suffix(".tmp")
temp.write_text(json.dumps(record,allow_nan=False)+"\n")
temp.replace(file)
print("Finished",file.name,record["ndet"],record["E_var_plus_PT2_Eh"],record["spin_squared"],flush=True)
