---
source_file: "DVR.f"
type: "code"
community: "STLM Lanczos Core"
location: "DVR.f:3294"
tags:
  - graphify/code
  - graphify/EXTRACTED
  - community/STLM_Lanczos_Core
---

# LANCZO (single Lanczos recurrence step)

## Connections
- [[MULVEC (managed vector scaling)]] - `calls` [EXTRACTED]
- [[OPINV (apply (K - muM)-1 via triangular solves)]] - `calls` [EXTRACTED]
- [[OPM (apply massoverlap operator M)]] - `calls` [EXTRACTED]
- [[SCPROD (managed scalar product)]] - `calls` [EXTRACTED]
- [[SELDOM (recover missing eigenvalues)]] - `calls` [EXTRACTED]
- [[SUBVEC (managed vector subtraction)]] - `calls` [EXTRACTED]
- [[Shift-invert Lanczos eigensolver (STLM library)]] - `implements` [EXTRACTED]
- [[always()]] - `calls` [EXTRACTED]

#graphify/code #graphify/EXTRACTED #community/STLM_Lanczos_Core