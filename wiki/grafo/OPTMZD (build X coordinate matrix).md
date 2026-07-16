---
source_file: "DVR.f"
type: "code"
community: "DVR Pipeline & Plots"
location: "DVR.f:437"
tags:
  - graphify/code
  - graphify/EXTRACTED
  - community/DVR_Pipeline__Plots
---

# OPTMZD (build X coordinate matrix)

## Connections
- [[DVR_OP (optimized DVR output to fort.11)]] - `shares_data_with` [INFERRED]
- [[EIGCLL (STLM wrapper and results writer)]] - `shares_data_with` [INFERRED]
- [[Optimized DVR  X-matrix pipeline (partially disabled)]] - `implements` [EXTRACTED]
- [[PROGRAM DVR (main driver)]] - `calls` [EXTRACTED]
- [[Packed lower-triangular storage in array W (H followed by overlap block)]] - `implements` [EXTRACTED]

#graphify/code #graphify/EXTRACTED #community/DVR_Pipeline__Plots