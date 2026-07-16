---
source_file: "DVR.f"
type: "code"
community: "Eigenvalue Output & Spectroscopy"
location: "DVR.f:190"
tags:
  - graphify/code
  - graphify/EXTRACTED
  - community/Eigenvalue_Output__Spectroscopy
---

# EIGCLL (STLM wrapper and results writer)

## Connections
- [[BUILD (assemble sinc-DVR Hamiltonian)]] - `shares_data_with` [EXTRACTED]
- [[EXCHG (read eigenvectors from unit 19)]] - `calls` [AMBIGUOUS]
- [[OPTMZD (build X coordinate matrix)]] - `shares_data_with` [INFERRED]
- [[PROGRAM DVR (main driver)]] - `calls` [EXTRACTED]
- [[STLM (shift-invert Lanczos main driver)]] - `shares_data_with` [EXTRACTED]
- [[Shift-invert Lanczos eigensolver (STLM library)]] - `conceptually_related_to` [EXTRACTED]
- [[Spectroscopic constants we, wexe, weye, alfae, gamae (finite differences of first 4 levels)]] - `implements` [EXTRACTED]
- [[Unit conversions 219474.631 hartree to cm-1; 1822.88839 amu to electron mass]] - `references` [EXTRACTED]
- [[fort.18  fort.19 (STLM scratch KFILE=18 sequential, DAFILE=19 direct-access)]] - `references` [EXTRACTED]
- [[fort.4 (output energies, spectrum, level differences, spectroscopic constants in cm-1)]] - `references` [EXTRACTED]
- [[fort.97 (output we in cm-1)]] - `references` [EXTRACTED]

#graphify/code #graphify/EXTRACTED #community/Eigenvalue_Output__Spectroscopy