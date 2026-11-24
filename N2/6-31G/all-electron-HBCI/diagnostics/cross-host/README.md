# Forte2 cross-host calibration check

We repeated the [2.10 Å reduced-basis singlet calculation](../../calibration/R-2.10/) on spark-6ad0 with the original Hamiltonian, two roots, four OpenMP threads, one OpenBLAS thread, and a variational threshold of 1.5625e-5. Both machines used [threshold.py](../../threshold.py) and its deterministic PT2 stopping rule.

The [comparison](results.json) records the two source revisions, the integral digest, determinant counts, and energy differences. The [raw program output](singlet-Ag-1.5625e-05.out) and [machine-readable calculation](singlet-Ag-1.5625e-05.json) document the repeat on spark-6ad0.
