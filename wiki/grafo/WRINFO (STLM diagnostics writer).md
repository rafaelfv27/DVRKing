---
source_file: "DVR.f"
type: "code"
community: "STLM Utilities & Diagnostics"
location: "DVR.f:6776"
tags:
  - graphify/code
  - graphify/EXTRACTED
  - community/STLM_Utilities__Diagnostics
---

# WRINFO (STLM diagnostics writer)

## Connections
- [[CHECK (validate STLM input parameters)]] - `calls` [EXTRACTED]
- [[DECOMP (LDLT decomposition of K - muM)]] - `calls` [EXTRACTED]
- [[FRSTIT (first-iteration setup and tolerance estimation)]] - `calls` [EXTRACTED]
- [[NEWPC (update Lanczos stepscopies optimization)]] - `calls` [EXTRACTED]
- [[SELDOM (recover missing eigenvalues)]] - `calls` [EXTRACTED]
- [[STLM (shift-invert Lanczos main driver)]] - `calls` [EXTRACTED]
- [[STLM minor helper routines (grouped LDL2LDL3LDLSUB, MULVMUL3, SUBVSUB3, SCP3, SOL2SOL3SOLVE, OPM2OPM3OPMSUB, HALF, MAXNRM, MCHECK, FREEID, INITAD, RAN3, ROTATE, SORTP, NEWMU, UI, WRINF3WRINF4)]] - `calls` [EXTRACTED]
- [[TRIDIG (tridiagonal analysis Ritz values, convergence, next shift)]] - `calls` [EXTRACTED]
- [[always()]] - `calls` [EXTRACTED]

#graphify/code #graphify/EXTRACTED #community/STLM_Utilities__Diagnostics