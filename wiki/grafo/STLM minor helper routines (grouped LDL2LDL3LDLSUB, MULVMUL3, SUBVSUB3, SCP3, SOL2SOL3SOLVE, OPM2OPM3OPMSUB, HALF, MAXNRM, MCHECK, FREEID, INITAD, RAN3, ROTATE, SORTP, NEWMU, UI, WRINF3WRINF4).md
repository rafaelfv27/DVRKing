---
source_file: "DVR.f"
type: "code"
community: "STLM Utilities & Diagnostics"
location: "DVR.f:1956-6775"
tags:
  - graphify/code
  - graphify/EXTRACTED
  - community/STLM_Utilities__Diagnostics
---

# STLM minor helper routines (grouped: LDL2/LDL3/LDLSUB, MULV/MUL3, SUBV/SUB3, SCP3, SOL2/SOL3/SOLVE, OPM2/OPM3/OPMSUB, HALF, MAXNRM, MCHECK, FREEID, INITAD, RAN3, ROTATE, SORTP, NEWMU, UI, WRINF3/WRINF4)

## Connections
- [[CHECK (validate STLM input parameters)]] - `calls` [EXTRACTED]
- [[FRSTIT (first-iteration setup and tolerance estimation)]] - `calls` [EXTRACTED]
- [[INITEW (initialize eigenvector workspace addresses)]] - `calls` [EXTRACTED]
- [[LDLT (LDLT factorization dispatcher)]] - `calls` [EXTRACTED]
- [[MULVEC (managed vector scaling)]] - `calls` [EXTRACTED]
- [[NEWPC (update Lanczos stepscopies optimization)]] - `calls` [EXTRACTED]
- [[OPINV (apply (K - muM)-1 via triangular solves)]] - `calls` [EXTRACTED]
- [[OPM (apply massoverlap operator M)]] - `calls` [EXTRACTED]
- [[RANDVC (generate random start vector)]] - `calls` [EXTRACTED]
- [[RANNUM (Knuth subtractive RNG)]] - `calls` [INFERRED]
- [[SUBVEC (managed vector subtraction)]] - `calls` [EXTRACTED]
- [[TERMIN (termination decision)]] - `calls` [EXTRACTED]
- [[TRIDIG (tridiagonal analysis Ritz values, convergence, next shift)]] - `calls` [EXTRACTED]
- [[WRINFO (STLM diagnostics writer)]] - `calls` [EXTRACTED]
- [[always()]] - `calls` [EXTRACTED]
- [[deldup()]] - `calls` [EXTRACTED]

#graphify/code #graphify/EXTRACTED #community/STLM_Utilities__Diagnostics