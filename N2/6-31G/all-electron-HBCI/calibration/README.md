# Reduced-basis controls

The two distance folders contain the [3s1p] nitrogen basis, MO Hamiltonian, exact FCI output, raw HBCI threshold calculations and extrapolation errors. All 14 electrons are active. Removing the outer p shell from each atom leaves 12 spatial orbitals. The [protocol JSON](results.json) records the criteria transferred to production.

### Current 0.1 µEh calibration

| R (Å) | Sector | Final selection threshold | Root 1 error (µEh) | Root 2 error (µEh) |
|---:|---|---:|---:|---:|
| 1.10 | singlet-Ag | 3.90625e-06 | +0.0046 | +0.0033 |
| 1.10 | triplet-B1u | 7.8125e-06 | +0.0025 | +0.0035 |
| 2.10 | singlet-Ag | 7.8125e-06 | -0.0018 | -0.0011 |
| 2.10 | triplet-B1u | 7.8125e-06 | +0.0011 | +0.0017 |

All eight exact residuals are below 0.001 µEh. The largest absolute extrapolation error plus its exact residual is 0.00463 µEh. The first eligible full-basis variational threshold is 3.90625e-6. Each fit records the individual threshold files, PT2 screening changes, uncertainties, and directly verified exact residuals. The files ending in -pt2 retain the stricter perturbative recalculations alongside the original datasets.

### Previous 1 µEh calibration

| R (Å) | Sector | Tightest selection threshold used | Root 1 error (µEh) | Root 2 error (µEh) |
|---:|---|---:|---:|---:|
| 1.10 | singlet-Ag | 3.125e-05 | +0.0089 | +0.0394 |
| 1.10 | triplet-B1u | 3.125e-05 | +0.0204 | +0.0300 |
| 2.10 | singlet-Ag | 1.5625e-05 | -0.0004 | +0.0027 |
| 2.10 | triplet-B1u | 3.125e-05 | -0.0065 | -0.0499 |

Errors are extrapolated HBCI total energies minus exact FCI total energies for the same Hamiltonian. The preceding 1 µEh protocol used 1.5625e-5 as its first eligible selection threshold. Every full-basis sector refines further until both roots pass their own stability, spin, PT2-screening and root-matching checks. [Methods and reproduction commands](../README.md)
