# N₂ all-electron heat-bath CI in 6-31G

[Full-basis results](results.json) · [Exact controls and threshold calibration](calibration/) · [Geometry grid, state definitions and integral options](../all-electron-FCI/input.json)

The targets are the two lowest singlet Ag roots and two lowest triplet B1u roots in D2h at each of the 23 bond lengths from 0.80 to 4.00 Å. Root labels denote energy order within a sector at each geometry. The production runs read the same 18-orbital, 14-electron Hamiltonians as the [FCI dataset](../all-electron-FCI/README.md). All total energies include nuclear repulsion. The integral Cholesky threshold is 1e-12.

## Live job console

On the Mac, double-click [monitor.command](monitor.command) or run `python3 monitor.py` from this directory. The console reads both Sparks over SSH every 30 seconds. It shows all 23 geometries, active selection thresholds, per-root uncertainties in µEh, and qualification against the current target. A question mark marks an excited-state assignment awaiting an independent rank check. Arrow keys select a geometry, Enter expands its solver and log details, **a** filters to active geometries, **r** refreshes, and **q** quits. The monitor reads job data without sending signals to the solvers.

Use `python3 monitor.py --once` for a plain-text snapshot or `python3 monitor.py --published --once` for the local checkout without SSH. If a host goes offline, its row identifies the data being shown as cached. Live mode uses the existing `spark-66ea` and `spark-6ad0` SSH aliases.

## Calibration and acceptance

At 1.10 and 2.10 Å we remove only the outer p contraction of nitrogen 6-31G, retaining all three s contractions and the inner p shell on each atom. This gives 12 spatial orbitals. Both N 1s orbitals remain in the correlation space. The original exact controls have directly verified residuals below 5e-8 Eh. For the tighter campaign, precision_fci.py requests 1e-9 Eh residuals on the same archived integrals.

Each HBCI calculation jointly selects determinants for two roots of a single spin and irrep. The initial determinant pool is restricted to the target irrep among the first ten orbitals. Selection can excite electrons throughout the entire orbital basis. The calculation projects the archived Hamiltonian onto D2h symmetry. Forbidden elements have zero matrix elements within the target sector. At 0.80 Å the largest forbidden two-electron integral is 1.47e-9 Eh. The conservative bound on the discarded operator over the full 14-electron space is below 1e-6 Eh at every geometry (14 times the one-body Frobenius norm plus 91 times the two-body Frobenius norm). The selected-CI spin projection is enabled. Each accepted result includes the calculated spin expectation.

Selection starts at 0.002 and the variational threshold is halved on each refinement. After variational selection stops adding determinants, the perturbative correction is recomputed at fixed wavefunction. Successive deterministic Epstein–Nesbet PT2 screening thresholds decrease by a factor of five from 1e-5 until both energies change by at most 0.02 microhartree, or the external sum is completely unscreened. Historical measurements with larger screening changes are retained. A precision rerun under the same selection threshold uses the suffix -pt2 and retains its own raw log. The raw JSON retains both roots, every PT2 evaluation, spin expectations, dominant determinants, and the SHA-256 of the input Hamiltonian.

The reported estimate is the intercept of a linear fit of four variational energies against their signed PT2 corrections. Its stability is checked against the newest three points and the preceding four, which exclude the newest calculation. A bound from the stored determinant amplitudes checks that the matched root has greater overlap than its partner across successive thresholds. The current target is **0.1 µEh (1e-7 Eh) per root**. A sector passes once both roots meet the spin, overlap, PT2, and empirical uncertainty checks. The error estimate is four times the sum of the fit spread and the largest fixed-wavefunction screening change over the five thresholds used for validation. The multiplier rounds up the largest observed ratio of error to the original fit estimate in the earlier 4.00 Å test (0.3173 / 0.0932). This precision change uses that endpoint as calibration evidence. Such estimates remain empirical in the full basis. The prior 1 µEh qualifications are available in [the previous dataset](https://github.com/fevangelista/MBRefData/tree/a181bdd/N2/6-31G/all-electron-HBCI).

In the reduced-basis calibration the absolute error plus the exact residual must also fit within the 0.1 µEh target. The tightest required selection threshold from the original two calibration geometries sets the first eligible full-basis threshold. The calculation continues to halve the selection threshold for each sector individually. A fit records its target, multiplier, and the exact filenames of its threshold inputs. The index includes only fits using the current protocol. The [4.00 Å comparison](validation/) and its four-root exact spectrum give an additional check in the dissociation region.

## Program and reproduction

The original 20 geometries use Forte2 2026.11.2 at upstream commit 8378e93380939a6c93ad1d647f6087845a8ff850 on spark-66ea. The 3.00, 3.50, and 4.00 Å extension uses a fresh build of origin/main at 692ae30049bebc71d5f1c4dc01f8bbbce596d751 on spark-6ad0, with Python 3.12.13, NumPy 2.4.4, and SciPy 1.17.1. The single intervening commit changes the Knowles–Handy FCI sigma builder and its symmetry test. Both N₂ FCI symmetry test cases and the multiroot N₂ selected-CI regression passed on the second machine. Repeating the [2.10 Å two-singlet calibration](diagnostics/cross-host/) on the same Hamiltonian selected 43,302 determinants on both hosts and reproduced the two variational-plus-PT2 energies within 7e-13 Eh. The exact-control inputs are in calibrate_fci.py. The HBCI and PT2 options are in threshold.py. The threshold sequence and acceptance rules are in campaign.py. Four calibration sectors can run independently using four Forte2 threads and one OpenBLAS thread per sector.

From this directory run the following in the Forte2 environment. Raw calculation logs are preserved alongside the corresponding machine-readable data.

```sh
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 HBCI_WORKERS=4
for r in 1.10 2.10; do
    mkdir -p "calibration/R-$r"
    python -u calibrate_fci.py "$r" > "calibration/R-$r/fci.out" 2>&1
done
python -u campaign.py
python summarize.py
```

The results index counts only roots passing the protocol. Per-distance results also include intermediate energies at completed thresholds so the convergence can be inspected.

## Following states along the curve

The [neighboring-geometry overlap survey](survey-tracking.json) compares the two roots in each sector between successive geometries using their wavefunctions at the common variational threshold of 0.001. The overlap includes the cross-geometry AO overlap and both MO coefficient matrices. The retained contribution uses up to 128 leading determinants per root, with separate alpha- and beta-spin determinant overlaps. The lower and upper bounds include the omitted coefficient norms and the known projections onto the opposite retained determinant space. The bounds concern these HBCI wavefunctions at the stated threshold.

On the original 20-point grid, the survey resolves 36 of 38 neighboring two-root assignments within the computed pairs. The singlet Ag assignments at 0.80–0.90 Å and 1.05–1.10 Å are ambiguous, driven by very small rank-2 overlaps with either state at the next geometry. Preserve energy rank at each distance and inspect higher roots if a continuous character assignment across those intervals is needed. For a comparison using only qualified calculations, run `python track_curve.py --output qualified-tracking.json` in the Forte2 environment. The survey can be reproduced by adding `--epsilon 0.001`.

A [four-root diagnostic](diagnostics/README.md) resolves the rank-2 to rank-3 character change between 1.05 and 1.10 Å and identifies a rank-4 to rank-2 match between 0.80 and 0.90 Å.

The extended survey contains 44 adjacent geometry/sector comparisons. On the 0.50 Å steps from 2.50 to 4.00 Å, the direct many-electron overlaps are small and the bounds from the leading determinants do not resolve either sector. A [second survey](survey-transported-tracking.json) compares CI character in Löwdin-orthogonalized AOs attached to the moving atoms. The AO labels are matched by atom, shell, and component. The MO coefficients in each Löwdin frame are S^(1/2) C. Identifying the two AO frames gives a unitary orbital map and the same determinant-tail bounds apply. Use `--epsilon 0.001 --transport` to reproduce this diagnostic.

In the atom-centered comparison, all six extended geometry/sector pairs favor rank 1 → 1 and rank 2 → 2 with separated bounds. The weakest retained same-rank overlap is 0.472 for the second triplet between 3.50 and 4.00 Å. These assignments compare the two calculated roots within each sector. The direct physical overlaps remain in the original survey file.

The [full-basis endpoint check](diagnostics/R-4.00/) tests the ordering of the second singlet at 4.00 Å with four independently selected roots. Independent four-root rank confirmations are required for both the second singlet and second triplet at 4.00 Å. In the reduced-basis validation, four-root selection resolved a 0.314804 µEh triplet splitting and recovered the second root to 0.00116 µEh.
