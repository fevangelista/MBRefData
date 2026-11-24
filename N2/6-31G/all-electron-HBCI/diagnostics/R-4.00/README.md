# Full-basis excited singlets at 4.00 Å

The independent [12-orbital spectrum](../../validation/R-4.00/singlet-spectrum.json) identified three closely spaced excited singlets. This four-root diagnostic checks the ordering in the full 18-orbital 6-31G reference.

At the initial comparison threshold of 0.0005, the two-root ground state matches the four-root ground state with overlap at least 0.919. The second two-root wavefunction overlaps both buffered ranks 2 and 4, with overlapping bounds. The four-root corrected energy of variational rank 3 lies 9.21 µEh below rank 2. These coarse calculations leave the identity of the second singlet open. The [comparison](results.json) records the overlap intervals and variational and corrected energies. The production index holds this excited state pending confirmation from tighter four-root calculations.

The Hamiltonian and basis are the same as for the [two-root production calculations](../../full/R-4.00/). Run the following in the dataset directory, changing the threshold for subsequent refinements.

```sh
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 python -u threshold.py full 4.0 singlet-Ag 0.0005 --buffer-roots 4
```
