---
source_file: "DVR.f"
type: "code"
community: "STLM Lanczos Core"
location: "DVR.f:549"
tags:
  - graphify/code
  - graphify/EXTRACTED
  - community/STLM_Lanczos_Core
---

# ALLOC (identifier-to-workspace storage manager)

## Connections
- [[ERROR (STLM error reporting)]] - `calls` [EXTRACTED]
- [[IO (direct-access scratch file IO handler)]] - `calls` [EXTRACTED]
- [[MULVEC (managed vector scaling)]] - `calls` [EXTRACTED]
- [[OPINV (apply (K - muM)-1 via triangular solves)]] - `calls` [EXTRACTED]
- [[OPM (apply massoverlap operator M)]] - `calls` [EXTRACTED]
- [[RANDVC (generate random start vector)]] - `calls` [EXTRACTED]
- [[SCPROD (managed scalar product)]] - `calls` [EXTRACTED]
- [[SUBVEC (managed vector subtraction)]] - `calls` [EXTRACTED]

#graphify/code #graphify/EXTRACTED #community/STLM_Lanczos_Core