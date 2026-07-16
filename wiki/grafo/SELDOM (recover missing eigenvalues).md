---
source_file: "DVR.f"
type: "code"
community: "STLM Lanczos Core"
location: "DVR.f:5020"
tags:
  - graphify/code
  - graphify/EXTRACTED
  - community/STLM_Lanczos_Core
---

# SELDOM (recover missing eigenvalues)

## Connections
- [[INITLA (initialize a Lanczos run)]] - `calls` [EXTRACTED]
- [[LANCZO (single Lanczos recurrence step)]] - `calls` [EXTRACTED]
- [[PREPSV (prepare saved vectors, orthogonalization)]] - `calls` [EXTRACTED]
- [[RANDVC (generate random start vector)]] - `calls` [EXTRACTED]
- [[SAVEXL (save converged eigenpairs)]] - `calls` [EXTRACTED]
- [[STLM (shift-invert Lanczos main driver)]] - `calls` [EXTRACTED]
- [[TRIDIG (tridiagonal analysis Ritz values, convergence, next shift)]] - `calls` [EXTRACTED]
- [[WRINFO (STLM diagnostics writer)]] - `calls` [EXTRACTED]

#graphify/code #graphify/EXTRACTED #community/STLM_Lanczos_Core