# All-electron N₂ FCI in 6-31G

[Results](results.json) · [Calculation options and bond lengths](input.json) · [Nitrogen basis definition](nitrogen-6-31g.json)

The curve covers 20 N–N distances between 0.80 and 2.50 Å. At each distance the target is the two lowest singlet roots in Ag and the two lowest triplet roots in B1u in D2h, with the bond on the z axis. The triplets use Ms = 1. Root numbers are one-based energy ranks within a spin and spatial symmetry sector at a single geometry. The dataset does not assign spectroscopic labels or establish continuity of wavefunction character at avoided crossings.

## Hamiltonian and solver

The nonrelativistic clamped-nuclei Hamiltonian correlates all 14 electrons in all 18 spatial orbitals. The nitrogen 6-31G basis contains only s and p functions. The integral preparation uses Forte2 Cholesky decomposition with a threshold of 1e-12, followed by transformation into symmetry-labeled canonical RHF orbitals. The RHF energy and density convergence tolerances are 1e-12 Eh and 1e-9, respectively.

The spin-adapted Forte2 FCI driver applies the Hamiltonian with the Harrison–Zarrabian algorithm and a 2048 MB sigma-builder setting. The Davidson–Liu solver requests two roots per sector with an energy tolerance of 1e-10 Eh, a residual tolerance of 1e-7 Eh, and up to 200 iterations. Each stored FCI eigenpair must also pass a separate application of the Hamiltonian to check its residual. Energies are in hartree and include nuclear repulsion.

## Files and provenance

Each R-* directory contains the Cartesian input in angstrom (geometry.xyz), the integral preparation output (prepare.out), the Hamiltonian (hamiltonian.npz), and results.json. The files singlet-Ag.out and triplet-B1u.out contain stdout and stderr for their respective calculation blocks as they become available. A pending block contributes zero energy entries to the index. Only converged and residual-checked FCI energies appear in fci_states. The RHF energy is stored separately so the integral preparation can be independently verified.

The compressed NumPy archive stores scalar_Eh (nuclear repulsion), h (one-electron integrals), V (spatial two-electron integrals in physicists’ order), C (AO-by-MO coefficients), and orbital_irreps (zero-based Cotton order). All integrals are in hartree. For a closed shell with occupied indices i,j the reference energy reconstructed from the archive is scalar_Eh + 2 sum_i h[i,i] + sum_ij (2 V[i,j,i,j] − V[i,j,j,i]).

The preparation and production machine was spark-66ea, Linux aarch64. Forte2 2026.11.2 was rebuilt in its main worktree at [8378e933](https://github.com/evangelistalab/forte2/commit/8378e93380939a6c93ad1d647f6087845a8ff850) in the forte2 conda environment. The N₂ symmetry CI regression passed after the rebuild. Python was 3.12.12, NumPy 2.3.5 and SciPy 1.16.3. The archived Basis Set Exchange nitrogen JSON fixes the basis independently of the installed BSE version. Integral preparation used four OpenMP threads and one OpenBLAS thread. Production FCI uses twelve OpenMP threads and four OpenBLAS threads, with OMP_MAX_ACTIVE_LEVELS=1.

## Reproducing the calculation

From this directory, with the specified Forte2 revision installed, prepare and solve a single block as follows.

```sh
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 python -u prepare.py 1.10 > R-1.10/prepare.out 2>&1
OMP_NUM_THREADS=12 OPENBLAS_NUM_THREADS=4 OMP_MAX_ACTIVE_LEVELS=1 python -u solve.py 1.10 singlet-Ag > R-1.10/singlet-Ag.out 2>&1
python summarize.py
```

For a new run from scratch use fresh output directories containing copies of the archived geometry.xyz files. With preparation complete for all distances, python -u run_grid.py runs one block at a time, prioritizing 1.10 and 2.10 Å and skipping sectors already saved. Rebuild the aggregate results.json with python summarize.py after transferring completed blocks.
