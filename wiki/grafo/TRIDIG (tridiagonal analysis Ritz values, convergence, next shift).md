---
source_file: "DVR.f"
type: "code"
community: "STLM Utilities & Diagnostics"
location: "DVR.f:5912"
tags:
  - graphify/code
  - graphify/EXTRACTED
  - community/STLM_Utilities__Diagnostics
---

# TRIDIG (tridiagonal analysis: Ritz values, convergence, next shift)

## Connections
- [[CONVER (convergence test for Ritz values)]] - `calls` [EXTRACTED]
- [[IMTQL2 (implicit QL tridiagonal eigensolver, EISPACK-style)]] - `calls` [EXTRACTED]
- [[SELDOM (recover missing eigenvalues)]] - `calls` [EXTRACTED]
- [[STLM minor helper routines (grouped LDL2LDL3LDLSUB, MULVMUL3, SUBVSUB3, SCP3, SOL2SOL3SOLVE, OPM2OPM3OPMSUB, HALF, MAXNRM, MCHECK, FREEID, INITAD, RAN3, ROTATE, SORTP, NEWMU, UI, WRINF3WRINF4)]] - `calls` [EXTRACTED]
- [[WRINFO (STLM diagnostics writer)]] - `calls` [EXTRACTED]
- [[always()]] - `calls` [EXTRACTED]
- [[compl()]] - `calls` [EXTRACTED]
- [[deldup()]] - `calls` [EXTRACTED]

#graphify/code #graphify/EXTRACTED #community/STLM_Utilities__Diagnostics