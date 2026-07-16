---
source_file: "DVR.f"
type: "code"
community: "Shift-Invert Solver Driver"
location: "DVR.f:5375"
tags:
  - graphify/code
  - graphify/EXTRACTED
  - community/Shift-Invert_Solver_Driver
---

# STLM (shift-invert Lanczos main driver)

## Connections
- [[DECOMP (LDLT decomposition of K - muM)]] - `calls` [EXTRACTED]
- [[EIGCLL (STLM wrapper and results writer)]] - `shares_data_with` [EXTRACTED]
- [[ERROR (STLM error reporting)]] - `calls` [EXTRACTED]
- [[FRSTIT (first-iteration setup and tolerance estimation)]] - `calls` [EXTRACTED]
- [[INITD (initialize STLM data and COMMON state)]] - `calls` [EXTRACTED]
- [[INITEW (initialize eigenvector workspace addresses)]] - `calls` [EXTRACTED]
- [[Packed lower-triangular storage in array W (H followed by overlap block)]] - `conceptually_related_to` [INFERRED]
- [[Python conversion plan (port physics core, replace STLM with SciPy)]] - `references` [EXTRACTED]
- [[REPL (update control variables between shifts)]] - `calls` [EXTRACTED]
- [[SELDOM (recover missing eigenvalues)]] - `calls` [EXTRACTED]
- [[SciPy eigensolver replacement (eigh subset_by_value  eigsh)]] - `semantically_similar_to` [INFERRED]
- [[Shift-invert Lanczos eigensolver (STLM library)]] - `implements` [EXTRACTED]
- [[TERMIN (termination decision)]] - `calls` [EXTRACTED]
- [[WRINFO (STLM diagnostics writer)]] - `calls` [EXTRACTED]
- [[always()]] - `calls` [EXTRACTED]
- [[second()]] - `calls` [EXTRACTED]

#graphify/code #graphify/EXTRACTED #community/Shift-Invert_Solver_Driver