# Reduced-basis controls

The two distance folders contain the [3s1p] nitrogen basis, MO Hamiltonian, exact FCI output, raw HBCI threshold calculations and extrapolation errors. All 14 electrons are active. Removing the outer p shell from each atom leaves 12 spatial orbitals. The [protocol JSON](results.json) records the criteria transferred to production.

| R (Å) | Sector | Tightest selection threshold used | Root 1 error (µEh) | Root 2 error (µEh) |
|---:|---|---:|---:|---:|
| 1.10 | singlet-Ag | 3.125e-05 | +0.0089 | +0.0394 |
| 1.10 | triplet-B1u | 3.125e-05 | +0.0204 | +0.0300 |
| 2.10 | singlet-Ag | 1.5625e-05 | -0.0004 | +0.0027 |
| 2.10 | triplet-B1u | 3.125e-05 | -0.0065 | -0.0499 |

Errors are extrapolated HBCI total energies minus exact FCI total energies for the same Hamiltonian. The full-basis protocol uses 1.5625e-5 as its first eligible selection threshold. Every full-basis sector refines further until both roots pass their own stability, spin, PT2-screening and root-matching checks. [Methods and reproduction commands](../README.md)
