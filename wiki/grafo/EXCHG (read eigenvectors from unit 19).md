---
source_file: "DVR.f"
type: "code"
community: "Eigenvalue Output & Spectroscopy"
location: "DVR.f:406"
tags:
  - graphify/code
  - graphify/AMBIGUOUS
  - community/Eigenvalue_Output__Spectroscopy
---

# EXCHG (read eigenvectors from unit 19)

## Connections
- [[EIGCLL (STLM wrapper and results writer)]] - `calls` [AMBIGUOUS]
- [[IO (direct-access scratch file IO handler)]] - `shares_data_with` [INFERRED]
- [[fort.18  fort.19 (STLM scratch KFILE=18 sequential, DAFILE=19 direct-access)]] - `references` [EXTRACTED]

#graphify/code #graphify/AMBIGUOUS #community/Eigenvalue_Output__Spectroscopy