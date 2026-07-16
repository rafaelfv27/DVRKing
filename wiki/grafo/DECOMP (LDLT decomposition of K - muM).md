---
source_file: "DVR.f"
type: "code"
community: "Shift-Invert Solver Driver"
location: "DVR.f:1380"
tags:
  - graphify/code
  - graphify/EXTRACTED
  - community/Shift-Invert_Solver_Driver
---

# DECOMP (LDL^T decomposition of K - mu*M)

## Connections
- [[LDLT (LDLT factorization dispatcher)]] - `calls` [EXTRACTED]
- [[STLM (shift-invert Lanczos main driver)]] - `calls` [EXTRACTED]
- [[WRINFO (STLM diagnostics writer)]] - `calls` [EXTRACTED]

#graphify/code #graphify/EXTRACTED #community/Shift-Invert_Solver_Driver