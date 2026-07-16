---
type: community
cohesion: 0.47
members: 6
---

# File I/O & Verification

**Cohesion:** 0.47 - moderately connected
**Members:** 6 nodes

## Members
- [[File IO map (Fortran unit numbers = fort.N filenames)]] - document - CLAUDE.md
- [[Regression verification against committed fort.4  fort.97]] - rationale - CLAUDE.md
- [[fort.3 (input grid AB, NPOINT, AMASS; ELOWEHIGHNUMEIG)]] - concept - DVR.f
- [[fort.4 (output energies, spectrum, level differences, spectroscopic constants in cm-1)]] - concept - DVR.f
- [[fort.97 (output we in cm-1)]] - concept - DVR.f
- [[fort.98  fort.99 (V(x) evaluation log)]] - concept - DVR.f

## Live Query (requires Dataview plugin)

```dataview
TABLE source_file, type FROM #community/File_I/O__Verification
SORT file.name ASC
```

## Connections to other communities
- 3 edges to [[_COMMUNITY_DVR Pipeline & Plots]]
- 3 edges to [[_COMMUNITY_Eigenvalue Output & Spectroscopy]]
- 3 edges to [[_COMMUNITY_Hamiltonian & Potential Physics]]

## Top bridge nodes
- [[File IO map (Fortran unit numbers = fort.N filenames)]] - degree 7, connects to 2 communities
- [[fort.4 (output energies, spectrum, level differences, spectroscopic constants in cm-1)]] - degree 4, connects to 2 communities
- [[Regression verification against committed fort.4  fort.97]] - degree 4, connects to 1 community
- [[fort.3 (input grid AB, NPOINT, AMASS; ELOWEHIGHNUMEIG)]] - degree 3, connects to 1 community
- [[fort.97 (output we in cm-1)]] - degree 3, connects to 1 community