---
source_file: "DVR.f"
type: "code"
community: "STLM Lanczos Core"
location: "DVR.f:1653"
tags:
  - graphify/code
  - graphify/EXTRACTED
  - community/STLM_Lanczos_Core
---

# ERROR (STLM error reporting)

## Connections
- [[ALLOC (identifier-to-workspace storage manager)]] - `calls` [EXTRACTED]
- [[CHECK (validate STLM input parameters)]] - `calls` [EXTRACTED]
- [[IMTQL2 (implicit QL tridiagonal eigensolver, EISPACK-style)]] - `calls` [EXTRACTED]
- [[IO (direct-access scratch file IO handler)]] - `calls` [EXTRACTED]
- [[STLM (shift-invert Lanczos main driver)]] - `calls` [EXTRACTED]

#graphify/code #graphify/EXTRACTED #community/STLM_Lanczos_Core