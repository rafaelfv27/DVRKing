---
source_file: "DVR.f"
type: "code"
community: "STLM Lanczos Core"
location: "DVR.f:4397"
tags:
  - graphify/code
  - graphify/EXTRACTED
  - community/STLM_Lanczos_Core
---

# PREPSV (prepare saved vectors, orthogonalization)

## Connections
- [[OPM (apply massoverlap operator M)]] - `calls` [EXTRACTED]
- [[SCPROD (managed scalar product)]] - `calls` [EXTRACTED]
- [[SELDOM (recover missing eigenvalues)]] - `calls` [EXTRACTED]
- [[SUBVEC (managed vector subtraction)]] - `calls` [EXTRACTED]
- [[always()]] - `calls` [EXTRACTED]

#graphify/code #graphify/EXTRACTED #community/STLM_Lanczos_Core