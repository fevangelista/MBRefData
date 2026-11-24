# Buffered singlet overlaps

At each endpoint of the two ambiguous intervals we targeted the four lowest singlet Ag roots at variational threshold 0.001, with the normal fixed-wavefunction PT2 checks. The initial determinant pool includes symmetry-compatible single and double excitations from RHF into all 18 orbitals. These calculations test continuity of character. The reference energy index uses the separately refined two-root sectors.

- At 0.80–0.90 Å the fourth root on the left matches the second on the right (truncated overlap magnitude 0.72). The leading candidate for the second root on the left is the fourth on the right (0.64), but the omitted-norm bounds leave that match ambiguous.
- At 1.05–1.10 Å the second root on the left matches the third on the right (0.90), with a unique assignment within the four computed roots under the omitted-norm bounds.

[Overlap matrices, bounds and individual row confidence](results.json) · [Methods](../README.md)

To reproduce a calculation, run `python threshold.py full 1.05 singlet-Ag 0.001 --buffer-roots 4` from the dataset directory. After all four points are available, run `python track_curve.py --buffer --epsilon 0.001 --output diagnostics/results.json`.

A separate [four-root check at 4.00 Å](R-4.00/README.md) tests the rank of the second full-basis singlet at the dissociation end of the curve.
