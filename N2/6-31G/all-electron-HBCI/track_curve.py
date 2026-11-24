"""Compare neighboring wavefunctions, accounting for orbital changes and omitted CI amplitudes."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment
from forte2 import System
from forte2.orbitals.orbital_overlap import mo_overlap

root = Path(__file__).resolve().parent
parent = root.parent / "all-electron-FCI"
settings = json.loads((parent/"input.json").read_text())
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--epsilon", type=float, help="Use a common survey threshold instead of qualified roots")
parser.add_argument("--buffer", action="store_true", help="Inspect the additional four-root singlet calculations")
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
os.chdir(parent)  # use the archived nitrogen basis with a relative, lower-case name


def basis_at(r):
    point = parent / f"R-{r:.2f}"
    xyz = "\n".join((point/"geometry.xyz").read_text().splitlines()[2:])
    system = System(xyz=xyz, unit="angstrom", basis_set=settings["basis_file"], symmetry=True)
    with np.load(point/"hamiltonian.npz") as data:
        C = data["C"]
    if not np.allclose(mo_overlap(C,system,C),np.eye(18),rtol=0,atol=1e-9):
        raise ValueError("Different AO convention from archived Hamiltonian")
    return system,C


def components(state, sector):
    dets = list(state)
    occ = [v.strip("|>") for v in dets]
    assert all(len(v)==18 and set(v)<=set("2ab0") for v in occ)
    alpha = np.array([[i for i,t in enumerate(v) if t in "a2"] for v in occ])
    beta = np.array([[i for i,t in enumerate(v) if t in "b2"] for v in occ])
    assert alpha.shape[1] == 7 + sector["ms"] and beta.shape[1] == 7 - sector["ms"]
    return alpha,beta,np.array([state[d] for d in dets])


def overlap(left,right,M):
    determinants=np.ones((len(left[2]),len(right[2])))
    for x,y in zip(left[:2],right[:2]):
        minors=M[x[:,None,:,None], y[None,:,None,:]]
        determinants *= np.linalg.det(minors)
    v=float(left[2] @ determinants @ right[2])
    nl,nr=np.linalg.norm(left[2]),np.linalg.norm(right[2])
    tl,tr=np.sqrt(max(0.,1-nl**2)),np.sqrt(max(0.,1-nr**2))
    # Subtract the known projection into the opposite retained determinant space.
    left_out=np.sqrt(max(0.,nl**2-np.linalg.norm(left[2] @ determinants)**2))
    right_out=np.sqrt(max(0.,nr**2-np.linalg.norm(determinants @ right[2])**2))
    bound=min(tl+nl*tr,tr+nr*tl,left_out*tr+right_out*tl+tl*tr)
    return abs(v), max(0.,abs(v)-bound), min(1.,abs(v)+bound)


def roots_at(r,sector):
    point=root/("diagnostics" if args.buffer else "full")/f"R-{r:.2f}"
    eps=args.epsilon
    if eps is None:
        fit=point/(sector["label"]+"-fit.json")
        if not fit.exists(): return None
        block=json.loads(fit.read_text())
        if not block["qualified"]: return None
        eps=block["epsilon_var"]
    f=point/f"{sector['label']}-{eps:.8g}.json"
    if not f.exists(): return None
    coeff=json.loads(f.read_text())["dominant_coefficients"]
    return eps,[components(c,sector) for c in coeff]


distances=settings["distances_angstrom"]
systems=[basis_at(r) for r in distances]
rows=[]
for i,(r0,r1) in enumerate(zip(distances,distances[1:])):
    if args.buffer and (r0,r1) not in [(0.8,0.9),(1.05,1.1)]: continue
    a,ca=systems[i]; b,cb=systems[i+1]
    M=mo_overlap(ca,a,cb,b)
    for sector in settings["states"]:
        left,right=roots_at(r0,sector),roots_at(r1,sector)
        if left is None or right is None: continue
        x=np.array([[overlap(l,t,M) for t in right[1]] for l in left[1]])
        _,best=linear_sum_assignment(-x[:,:,0])
        certified=[bool(x[k,j,1]>max(np.delete(x[k,:,2],j))) for k,j in enumerate(best)]
        unique=all(certified)
        rows.append(dict(left_angstrom=r0,right_angstrom=r1,sector=sector["label"],
                         epsilon_left=left[0],epsilon_right=right[0],
                         overlap_magnitude_truncated=x[:,:,0].tolist(),
                         lower=x[:,:,1].tolist(),upper=x[:,:,2].tolist(),
                         best_matching_right_ranks=[int(j)+1 for j in best],
                         row_match_certified=certified,
                         certified_within_root_pair=bool(unique)))
result=dict(threshold_mode="qualified" if args.epsilon is None else "common survey threshold",
            neighboring_sector_pairs=len(rows),pairs=rows)
outfile=args.output if args.output.is_absolute() else root/args.output
outfile.write_text(json.dumps(result,indent=2)+"\n")
print("Neighbor pairs",len(rows),"certified within the two roots",
      sum(row["certified_within_root_pair"] for row in rows))
