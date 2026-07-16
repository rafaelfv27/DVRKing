---
type: community
cohesion: 0.31
members: 17
---

# STLM Lanczos Core

**Cohesion:** 0.31 - loosely connected
**Members:** 17 nodes

## Members
- [[ALLOC (identifier-to-workspace storage manager)]] - code - DVR.f
- [[ERROR (STLM error reporting)]] - code - DVR.f
- [[FRSTIT (first-iteration setup and tolerance estimation)]] - code - DVR.f
- [[IMTQL2 (implicit QL tridiagonal eigensolver, EISPACK-style)]] - code - DVR.f
- [[INITLA (initialize a Lanczos run)]] - code - DVR.f
- [[IO (direct-access scratch file IO handler)]] - code - DVR.f
- [[LANCZO (single Lanczos recurrence step)]] - code - DVR.f
- [[MULVEC (managed vector scaling)]] - code - DVR.f
- [[OPINV (apply (K - muM)-1 via triangular solves)]] - code - DVR.f
- [[OPM (apply massoverlap operator M)]] - code - DVR.f
- [[PREPSV (prepare saved vectors, orthogonalization)]] - code - DVR.f
- [[RANDVC (generate random start vector)]] - code - DVR.f
- [[SAVEXL (save converged eigenpairs)]] - code - DVR.f
- [[SCPROD (managed scalar product)]] - code - DVR.f
- [[SELDOM (recover missing eigenvalues)]] - code - DVR.f
- [[SUBVEC (managed vector subtraction)]] - code - DVR.f
- [[TRANSF (Ritz vector transformation)]] - code - DVR.f

## Live Query (requires Dataview plugin)

```dataview
TABLE source_file, type FROM #community/STLM_Lanczos_Core
SORT file.name ASC
```

## Connections to other communities
- 16 edges to [[_COMMUNITY_STLM Utilities & Diagnostics]]
- 5 edges to [[_COMMUNITY_Shift-Invert Solver Driver]]
- 2 edges to [[_COMMUNITY_Eigenvalue Output & Spectroscopy]]

## Top bridge nodes
- [[FRSTIT (first-iteration setup and tolerance estimation)]] - degree 10, connects to 2 communities
- [[IO (direct-access scratch file IO handler)]] - degree 9, connects to 2 communities
- [[LANCZO (single Lanczos recurrence step)]] - degree 8, connects to 2 communities
- [[SELDOM (recover missing eigenvalues)]] - degree 8, connects to 2 communities
- [[ERROR (STLM error reporting)]] - degree 5, connects to 2 communities