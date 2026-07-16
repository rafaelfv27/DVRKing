---
type: community
cohesion: 0.67
members: 4
---

# Eigenvalue Output & Spectroscopy

**Cohesion:** 0.67 - moderately connected
**Members:** 4 nodes

## Members
- [[EIGCLL (STLM wrapper and results writer)]] - code - DVR.f
- [[EXCHG (read eigenvectors from unit 19)]] - code - DVR.f
- [[Spectroscopic constants we, wexe, weye, alfae, gamae (finite differences of first 4 levels)]] - concept - DVR.f
- [[fort.18  fort.19 (STLM scratch KFILE=18 sequential, DAFILE=19 direct-access)]] - concept - DVR.f

## Live Query (requires Dataview plugin)

```dataview
TABLE source_file, type FROM #community/Eigenvalue_Output__Spectroscopy
SORT file.name ASC
```

## Connections to other communities
- 3 edges to [[_COMMUNITY_Hamiltonian & Potential Physics]]
- 3 edges to [[_COMMUNITY_File IO & Verification]]
- 2 edges to [[_COMMUNITY_DVR Pipeline & Plots]]
- 2 edges to [[_COMMUNITY_Shift-Invert Solver Driver]]
- 2 edges to [[_COMMUNITY_STLM Lanczos Core]]

## Top bridge nodes
- [[EIGCLL (STLM wrapper and results writer)]] - degree 11, connects to 4 communities
- [[fort.18  fort.19 (STLM scratch KFILE=18 sequential, DAFILE=19 direct-access)]] - degree 4, connects to 2 communities
- [[EXCHG (read eigenvectors from unit 19)]] - degree 3, connects to 1 community
- [[Spectroscopic constants we, wexe, weye, alfae, gamae (finite differences of first 4 levels)]] - degree 2, connects to 1 community