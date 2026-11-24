# Full-basis excited states at 4.00 Å

The independent [12-orbital spectrum](../../validation/R-4.00/singlet-spectrum.json) identified three closely spaced excited singlets. This four-root diagnostic checks the ordering in the full 18-orbital 6-31G reference.

At the initial comparison threshold of 0.0005, the two-root ground state matches the four-root ground state with overlap at least 0.919. The second two-root wavefunction overlaps both buffered ranks 2 and 4, with overlapping bounds. The four-root corrected energy of variational rank 3 lies 9.21 µEh below rank 2. These coarse calculations leave the identity of the second singlet open. The [initial comparison](initial-comparison.json) records the overlap intervals and variational and corrected energies. The production index holds this excited state pending confirmation from tighter four-root calculations. The [precision driver](../check_endpoint.py) applies the 0.1 µEh uncertainty checks, tests each root's overlap against all competitors, and checks that the extrapolated intervals for rank 2 separate from ranks 1, 3, and 4. Once these checks pass, the index takes the excited-singlet energy directly from that four-root calculation. The [running result](results.json) records its status.

The Hamiltonian and basis are the same as for the [two-root production calculations](../../full/R-4.00/). Run the following in the dataset directory, changing the threshold for subsequent refinements.

```sh
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 python -u threshold.py full 4.0 singlet-Ag 0.0005 --buffer-roots 4
```

For the second triplet, the [reduced-basis four-root control](../../validation/R-4.00/buffered/) detected a 0.314804 µEh spacing to rank 3. The [full-basis triplet verification](triplet-rank-check.json) applies the same rank and precision checks. Reproduce it with `python diagnostics/check_endpoint.py --triplet`.
