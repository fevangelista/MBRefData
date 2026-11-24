# Independent 4.00 Å control

This stretched geometry tests the threshold protocol calibrated at 1.10 and 2.10 Å. The control uses all 14 electrons in the same 12-orbital 6-31G subset, with the outer p contraction removed on each nitrogen. The calculation ran on spark-6ad0 with Forte2 commit 692ae30049bebc71d5f1c4dc01f8bbbce596d751, two OpenMP threads per solver, and one OpenBLAS thread.

The HBCI stopping rule uses only energies, PT2 screening convergence, spin expectations, and root-overlap bounds from the threshold sequence. The [comparison script](check.py) reads the independent exact energies after applying this rule. The [results](results.json) report signed errors (HBCI extrapolation minus exact FCI) and estimated uncertainties in hartree.

### Current 0.1 µEh check

The first application of the tighter protocol passed for three of four states. The second triplet's estimated uncertainty was 0.0469 µEh at a selection threshold of 3.90625e-6. Its error against the tightened exact calculation was 0.3159 µEh. Three extra threshold halvings reduced the estimate to 0.0145 µEh with a measured error of 0.3149 µEh. The [first comparison](R-4.00/first-strict-check.json) preserves this failure. The [four-triplet exact spectrum](R-4.00/triplet-spectrum.json) resolves ranks 2 and 3 with a gap of 0.314804 µEh. The two-root HBCI calculation followed rank 3. A [four-root HBCI calculation](R-4.00/buffered/) resolves rank 2 at a selection threshold of 3.90625e-6. Its fit estimates 0.0555 µEh uncertainty and differs from exact rank 2 by 0.00116 µEh. With this assignment all four endpoint controls pass the 0.1 µEh target. The [results](results.json) identify the fit used for each state.

### Previous 1 µEh check

Under the previous 1 µEh stopping rule, all four errors were below 1 µEh. The second triplet has the largest observed error, 0.317 µEh, against a fit-based estimate of 0.093 µEh. Thus the fit uncertainty can underestimate the error even when the reference meets the 1 µEh target. In that run the triplets qualified at a selection threshold of 1.5625e-5. Consecutive overlap checks for the second singlet required refinement to 1.220703125e-7.

A [four-singlet exact calculation](R-4.00/singlet-spectrum.json) places roots 2 and 3 only 0.929 µEh apart and roots 3 and 4 a further 1.050 µEh apart. The largest verified residual among these four eigenpairs is below 0.000026 µEh. Several coarser HBCI threshold steps switch excited-state character in this narrow spacing. The reported singlet energies pass the consecutive-overlap criterion at the final threshold.

From the HBCI dataset directory, use the Forte2 environment and run

```sh
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
mkdir -p validation/R-4.00
python -u calibrate_fci.py 4.00 --validation > validation/R-4.00/fci.out 2>&1
python -u precision_fci.py validation 4.00 > validation/R-4.00/fci-tight.out 2>&1
python -u diagnostics/check_endpoint.py --triplet --validation
python -u validation/check.py > validation/R-4.00/comparison-buffered.out 2>&1
python -u validation/singlet_spectrum.py > validation/R-4.00/singlet-spectrum.out 2>&1
python -u validation/singlet_spectrum.py --triplet > validation/R-4.00/triplet-spectrum.out 2>&1
```
