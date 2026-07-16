---
source_file: "DVR.f"
type: "code"
community: "Shift-Invert Solver Driver"
location: "DVR.f:3662"
tags:
  - graphify/code
  - graphify/EXTRACTED
  - community/Shift-Invert_Solver_Driver
---

# LDLT (LDL^T factorization dispatcher)

## Connections
- [[DECOMP (LDLT decomposition of K - muM)]] - `calls` [EXTRACTED]
- [[IO (direct-access scratch file IO handler)]] - `calls` [EXTRACTED]
- [[STLM minor helper routines (grouped LDL2LDL3LDLSUB, MULVMUL3, SUBVSUB3, SCP3, SOL2SOL3SOLVE, OPM2OPM3OPMSUB, HALF, MAXNRM, MCHECK, FREEID, INITAD, RAN3, ROTATE, SORTP, NEWMU, UI, WRINF3WRINF4)]] - `calls` [EXTRACTED]
- [[Shift-invert Lanczos eigensolver (STLM library)]] - `conceptually_related_to` [INFERRED]

#graphify/code #graphify/EXTRACTED #community/Shift-Invert_Solver_Driver