# N₂ at 4.00 Å in the 12-orbital validation basis

The [results](results.json) compare both singlets and both triplets to [residual-checked exact FCI](exact.json). The geometry is in the [shared input](../../../all-electron-FCI/R-4.00/geometry.xyz). This folder contains the reduced basis, the all-electron Hamiltonian, the exact calculation stdout, and the HBCI raw output and JSON for each selection threshold. See the [validation procedure](../README.md) for the independent stopping rule and reproduction commands.

The [four-root singlet spectrum](singlet-spectrum.json) gives exact adjacent gaps and eigenpair residuals for the near-degenerate excited roots. Its [raw output](singlet-spectrum.out) and [input script](../singlet_spectrum.py) are archived.
