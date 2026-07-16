---
type: community
cohesion: 0.27
members: 11
---

# DVR Pipeline & Plots

**Cohesion:** 0.27 - loosely connected
**Members:** 11 nodes

## Members
- [[DVR_OP (optimized DVR output to fort.11)]] - code - DVR.f
- [[DVR_PR (plot primitive DVR basis function)]] - code - DVR.f
- [[FBR (finite basis representation plot, stub)]] - code - DVR.f
- [[OPTMZD (build X coordinate matrix)]] - code - DVR.f
- [[Optimized DVR  X-matrix pipeline (partially disabled)]] - concept - DVR.f
- [[PROGRAM DVR (main driver)]] - code - DVR.f
- [[Packed lower-triangular storage in array W (H followed by overlap block)]] - concept - DVR.f
- [[READ_I (read fort.3 input)]] - code - DVR.f
- [[Sinc-DVR method (Colbert-Miller, J. Chem. Phys. 96, 1982 (1992))]] - concept - DVR.f
- [[fort.10 (output primitive DVR function plot)]] - concept - DVR.f
- [[fort.11 (output grid x vs V from DVR_OP)]] - concept - DVR.f

## Live Query (requires Dataview plugin)

```dataview
TABLE source_file, type FROM #community/DVR_Pipeline__Plots
SORT file.name ASC
```

## Connections to other communities
- 4 edges to [[_COMMUNITY_Hamiltonian & Potential Physics]]
- 3 edges to [[_COMMUNITY_File IO & Verification]]
- 2 edges to [[_COMMUNITY_Eigenvalue Output & Spectroscopy]]
- 1 edge to [[_COMMUNITY_Shift-Invert Solver Driver]]

## Top bridge nodes
- [[PROGRAM DVR (main driver)]] - degree 8, connects to 2 communities
- [[Packed lower-triangular storage in array W (H followed by overlap block)]] - degree 3, connects to 2 communities
- [[DVR_OP (optimized DVR output to fort.11)]] - degree 5, connects to 1 community
- [[OPTMZD (build X coordinate matrix)]] - degree 5, connects to 1 community
- [[Sinc-DVR method (Colbert-Miller, J. Chem. Phys. 96, 1982 (1992))]] - degree 4, connects to 1 community