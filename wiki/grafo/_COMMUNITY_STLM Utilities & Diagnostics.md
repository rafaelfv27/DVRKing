---
type: community
cohesion: 0.28
members: 15
---

# STLM Utilities & Diagnostics

**Cohesion:** 0.28 - loosely connected
**Members:** 15 nodes

## Members
- [[CHECK (validate STLM input parameters)]] - code - DVR.f
- [[CONVER (convergence test for Ritz values)]] - code - DVR.f
- [[DVR.f]] - code - DVR.f
- [[NEWPC (update Lanczos stepscopies optimization)]] - code - DVR.f
- [[RANNUM (Knuth subtractive RNG)]] - code - DVR.f
- [[STLM minor helper routines (grouped LDL2LDL3LDLSUB, MULVMUL3, SUBVSUB3, SCP3, SOL2SOL3SOLVE, OPM2OPM3OPMSUB, HALF, MAXNRM, MCHECK, FREEID, INITAD, RAN3, ROTATE, SORTP, NEWMU, UI, WRINF3WRINF4)]] - code - DVR.f
- [[TRIDIG (tridiagonal analysis Ritz values, convergence, next shift)]] - code - DVR.f
- [[WRINFO (STLM diagnostics writer)]] - code - DVR.f
- [[always()]] - code - DVR.f
- [[cmprss()]] - code - DVR.f
- [[compl()]] - code - DVR.f
- [[dcheck()]] - code - DVR.f
- [[deldup()]] - code - DVR.f
- [[ranfx (RNG wrapper over rand())]] - code - DVR.f
- [[second()]] - code - DVR.f

## Live Query (requires Dataview plugin)

```dataview
TABLE source_file, type FROM #community/STLM_Utilities__Diagnostics
SORT file.name ASC
```

## Connections to other communities
- 16 edges to [[_COMMUNITY_STLM Lanczos Core]]
- 8 edges to [[_COMMUNITY_Shift-Invert Solver Driver]]

## Top bridge nodes
- [[STLM minor helper routines (grouped LDL2LDL3LDLSUB, MULVMUL3, SUBVSUB3, SCP3, SOL2SOL3SOLVE, OPM2OPM3OPMSUB, HALF, MAXNRM, MCHECK, FREEID, INITAD, RAN3, ROTATE, SORTP, NEWMU, UI, WRINF3WRINF4)]] - degree 16, connects to 2 communities
- [[always()]] - degree 12, connects to 2 communities
- [[WRINFO (STLM diagnostics writer)]] - degree 9, connects to 2 communities
- [[CHECK (validate STLM input parameters)]] - degree 5, connects to 2 communities
- [[TRIDIG (tridiagonal analysis Ritz values, convergence, next shift)]] - degree 8, connects to 1 community