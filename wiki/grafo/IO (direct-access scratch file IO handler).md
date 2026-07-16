---
source_file: "DVR.f"
type: "code"
community: "STLM Lanczos Core"
location: "DVR.f:3175"
tags:
  - graphify/code
  - graphify/EXTRACTED
  - community/STLM_Lanczos_Core
---

# IO (direct-access scratch file I/O handler)

## Connections
- [[ALLOC (identifier-to-workspace storage manager)]] - `calls` [EXTRACTED]
- [[ERROR (STLM error reporting)]] - `calls` [EXTRACTED]
- [[EXCHG (read eigenvectors from unit 19)]] - `shares_data_with` [INFERRED]
- [[LDLT (LDLT factorization dispatcher)]] - `calls` [EXTRACTED]
- [[MULVEC (managed vector scaling)]] - `calls` [EXTRACTED]
- [[RANDVC (generate random start vector)]] - `calls` [EXTRACTED]
- [[SCPROD (managed scalar product)]] - `calls` [EXTRACTED]
- [[SUBVEC (managed vector subtraction)]] - `calls` [EXTRACTED]
- [[fort.18  fort.19 (STLM scratch KFILE=18 sequential, DAFILE=19 direct-access)]] - `references` [INFERRED]

#graphify/code #graphify/EXTRACTED #community/STLM_Lanczos_Core