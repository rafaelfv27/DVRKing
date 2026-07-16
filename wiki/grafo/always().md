---
source_file: "DVR.f"
type: "code"
community: "STLM Utilities & Diagnostics"
location: "L680"
tags:
  - graphify/code
  - graphify/EXTRACTED
  - community/STLM_Utilities__Diagnostics
---

# always()

## Connections
- [[DVR.f]] - `defines` [EXTRACTED]
- [[INITLA (initialize a Lanczos run)]] - `calls` [EXTRACTED]
- [[LANCZO (single Lanczos recurrence step)]] - `calls` [EXTRACTED]
- [[NEWPC (update Lanczos stepscopies optimization)]] - `calls` [EXTRACTED]
- [[PREPSV (prepare saved vectors, orthogonalization)]] - `calls` [EXTRACTED]
- [[RANDVC (generate random start vector)]] - `calls` [EXTRACTED]
- [[SAVEXL (save converged eigenpairs)]] - `calls` [EXTRACTED]
- [[STLM (shift-invert Lanczos main driver)]] - `calls` [EXTRACTED]
- [[STLM minor helper routines (grouped LDL2LDL3LDLSUB, MULVMUL3, SUBVSUB3, SCP3, SOL2SOL3SOLVE, OPM2OPM3OPMSUB, HALF, MAXNRM, MCHECK, FREEID, INITAD, RAN3, ROTATE, SORTP, NEWMU, UI, WRINF3WRINF4)]] - `calls` [EXTRACTED]
- [[TRIDIG (tridiagonal analysis Ritz values, convergence, next shift)]] - `calls` [EXTRACTED]
- [[WRINFO (STLM diagnostics writer)]] - `calls` [EXTRACTED]
- [[second()]] - `calls` [EXTRACTED]

#graphify/code #graphify/EXTRACTED #community/STLM_Utilities__Diagnostics