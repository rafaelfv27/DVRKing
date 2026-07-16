---
type: community
cohesion: 0.28
members: 9
---

# Shift-Invert Solver Driver

**Cohesion:** 0.28 - loosely connected
**Members:** 9 nodes

## Members
- [[DECOMP (LDLT decomposition of K - muM)]] - code - DVR.f
- [[INITD (initialize STLM data and COMMON state)]] - code - DVR.f
- [[INITEW (initialize eigenvector workspace addresses)]] - code - DVR.f
- [[LDLT (LDLT factorization dispatcher)]] - code - DVR.f
- [[REPL (update control variables between shifts)]] - code - DVR.f
- [[STLM (shift-invert Lanczos main driver)]] - code - DVR.f
- [[SciPy eigensolver replacement (eigh subset_by_value  eigsh)]] - rationale - CLAUDE.md
- [[Shift-invert Lanczos eigensolver (STLM library)]] - concept - DVR.f
- [[TERMIN (termination decision)]] - code - DVR.f

## Live Query (requires Dataview plugin)

```dataview
TABLE source_file, type FROM #community/Shift-Invert_Solver_Driver
SORT file.name ASC
```

## Connections to other communities
- 8 edges to [[_COMMUNITY_STLM Utilities & Diagnostics]]
- 5 edges to [[_COMMUNITY_STLM Lanczos Core]]
- 2 edges to [[_COMMUNITY_Eigenvalue Output & Spectroscopy]]
- 2 edges to [[_COMMUNITY_Hamiltonian & Potential Physics]]
- 1 edge to [[_COMMUNITY_DVR Pipeline & Plots]]

## Top bridge nodes
- [[STLM (shift-invert Lanczos main driver)]] - degree 16, connects to 5 communities
- [[Shift-invert Lanczos eigensolver (STLM library)]] - degree 5, connects to 2 communities
- [[LDLT (LDLT factorization dispatcher)]] - degree 4, connects to 2 communities
- [[SciPy eigensolver replacement (eigh subset_by_value  eigsh)]] - degree 3, connects to 1 community
- [[DECOMP (LDLT decomposition of K - muM)]] - degree 3, connects to 1 community