---
source_file: "DVR.f"
type: "code"
community: "STLM Lanczos Core"
location: "DVR.f:3070"
tags:
  - graphify/code
  - graphify/EXTRACTED
  - community/STLM_Lanczos_Core
---

# INITLA (initialize a Lanczos run)

## Connections
- [[MULVEC (managed vector scaling)]] - `calls` [EXTRACTED]
- [[OPM (apply massoverlap operator M)]] - `calls` [EXTRACTED]
- [[SCPROD (managed scalar product)]] - `calls` [EXTRACTED]
- [[SELDOM (recover missing eigenvalues)]] - `calls` [EXTRACTED]
- [[always()]] - `calls` [EXTRACTED]

#graphify/code #graphify/EXTRACTED #community/STLM_Lanczos_Core