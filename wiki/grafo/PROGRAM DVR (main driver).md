---
source_file: "DVR.f"
type: "code"
community: "DVR Pipeline & Plots"
location: "DVR.f:5"
tags:
  - graphify/code
  - graphify/EXTRACTED
  - community/DVR_Pipeline__Plots
---

# PROGRAM DVR (main driver)

## Connections
- [[BUILD (assemble sinc-DVR Hamiltonian)]] - `calls` [EXTRACTED]
- [[DVR_OP (optimized DVR output to fort.11)]] - `calls` [EXTRACTED]
- [[DVR_PR (plot primitive DVR basis function)]] - `calls` [EXTRACTED]
- [[EIGCLL (STLM wrapper and results writer)]] - `calls` [EXTRACTED]
- [[FBR (finite basis representation plot, stub)]] - `calls` [EXTRACTED]
- [[OPTMZD (build X coordinate matrix)]] - `calls` [EXTRACTED]
- [[READ_I (read fort.3 input)]] - `calls` [EXTRACTED]
- [[Sinc-DVR method (Colbert-Miller, J. Chem. Phys. 96, 1982 (1992))]] - `implements` [EXTRACTED]

#graphify/code #graphify/EXTRACTED #community/DVR_Pipeline__Plots