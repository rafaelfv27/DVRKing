---
source_file: "DVR.f"
type: "code"
community: "STLM Lanczos Core"
location: "DVR.f:4526"
tags:
  - graphify/code
  - graphify/EXTRACTED
  - community/STLM_Lanczos_Core
---

# RANDVC (generate random start vector)

## Connections
- [[ALLOC (identifier-to-workspace storage manager)]] - `calls` [EXTRACTED]
- [[FRSTIT (first-iteration setup and tolerance estimation)]] - `calls` [EXTRACTED]
- [[IO (direct-access scratch file IO handler)]] - `calls` [EXTRACTED]
- [[OPINV (apply (K - muM)-1 via triangular solves)]] - `calls` [EXTRACTED]
- [[OPM (apply massoverlap operator M)]] - `calls` [EXTRACTED]
- [[SELDOM (recover missing eigenvalues)]] - `calls` [EXTRACTED]
- [[STLM minor helper routines (grouped LDL2LDL3LDLSUB, MULVMUL3, SUBVSUB3, SCP3, SOL2SOL3SOLVE, OPM2OPM3OPMSUB, HALF, MAXNRM, MCHECK, FREEID, INITAD, RAN3, ROTATE, SORTP, NEWMU, UI, WRINF3WRINF4)]] - `calls` [EXTRACTED]
- [[always()]] - `calls` [EXTRACTED]

#graphify/code #graphify/EXTRACTED #community/STLM_Lanczos_Core