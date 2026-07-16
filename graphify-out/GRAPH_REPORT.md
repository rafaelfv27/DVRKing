# Graph Report - .  (2026-07-13)

## Corpus Check
- Corpus is ~28,156 words - fits in a single context window. You may not need a graph.

## Summary
- 69 nodes · 172 edges · 7 communities
- Extraction: 92% EXTRACTED · 7% INFERRED · 1% AMBIGUOUS · INFERRED: 12 edges (avg confidence: 0.87)
- Token cost: 122,340 input · 0 output

## Community Hubs (Navigation)
- STLM Lanczos Core
- STLM Utilities & Diagnostics
- DVR Pipeline & Plots
- Shift-Invert Solver Driver
- Hamiltonian & Potential Physics
- File I/O & Verification
- Eigenvalue Output & Spectroscopy

## God Nodes (most connected - your core abstractions)
1. `STLM (shift-invert Lanczos main driver)` - 16 edges
2. `STLM minor helper routines (grouped: LDL2/LDL3/LDLSUB, MULV/MUL3, SUBV/SUB3, SCP3, SOL2/SOL3/SOLVE, OPM2/OPM3/OPMSUB, HALF, MAXNRM, MCHECK, FREEID, INITAD, RAN3, ROTATE, SORTP, NEWMU, UI, WRINF3/WRINF4)` - 16 edges
3. `always()` - 12 edges
4. `EIGCLL (STLM wrapper and results writer)` - 11 edges
5. `FRSTIT (first-iteration setup and tolerance estimation)` - 10 edges
6. `Python conversion plan (port physics core, replace STLM with SciPy)` - 10 edges
7. `IO (direct-access scratch file I/O handler)` - 9 edges
8. `WRINFO (STLM diagnostics writer)` - 9 edges
9. `PROGRAM DVR (main driver)` - 8 edges
10. `ALLOC (identifier-to-workspace storage manager)` - 8 edges

## Surprising Connections (you probably didn't know these)
- `SciPy eigensolver replacement (eigh subset_by_value / eigsh)` --semantically_similar_to--> `STLM (shift-invert Lanczos main driver)`  [INFERRED] [semantically similar]
  CLAUDE.md → DVR.f
- `Python conversion plan (port physics core, replace STLM with SciPy)` --references--> `Spectroscopic constants we, wexe, weye, alfae, gamae (finite differences of first 4 levels)`  [EXTRACTED]
  CLAUDE.md → DVR.f
- `File I/O map (Fortran unit numbers = fort.N filenames)` --references--> `fort.10 (output: primitive DVR function plot)`  [EXTRACTED]
  CLAUDE.md → DVR.f
- `File I/O map (Fortran unit numbers = fort.N filenames)` --references--> `fort.11 (output: grid x vs V from DVR_OP)`  [EXTRACTED]
  CLAUDE.md → DVR.f
- `Python conversion plan (port physics core, replace STLM with SciPy)` --references--> `STLM (shift-invert Lanczos main driver)`  [EXTRACTED]
  CLAUDE.md → DVR.f

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **STLM shift-invert Lanczos main iteration loop** — dvr_stlm, dvr_decomp, dvr_frstit, dvr_initew, dvr_always, dvr_seldom, dvr_termin, dvr_repl [EXTRACTED 1.00]
- **STLM routines coupled through /STLM../ COMMON blocks** — dvr_stlm, dvr_initd, dvr_decomp, dvr_frstit, dvr_always, dvr_seldom, dvr_tridig, dvr_lanczo, dvr_ldlt, dvr_termin, dvr_wrinfo, dvr_io [EXTRACTED 1.00]
- **Optimized DVR pipeline (second diagonalization commented out)** — dvr_fbr, dvr_optmzd, dvr_dvr_op [EXTRACTED 1.00]

## Communities (7 total, 0 thin omitted)

### Community 0 - "STLM Lanczos Core"
Cohesion: 0.31
Nodes (17): ALLOC (identifier-to-workspace storage manager), ERROR (STLM error reporting), FRSTIT (first-iteration setup and tolerance estimation), IMTQL2 (implicit QL tridiagonal eigensolver, EISPACK-style), INITLA (initialize a Lanczos run), IO (direct-access scratch file I/O handler), LANCZO (single Lanczos recurrence step), MULVEC (managed vector scaling) (+9 more)

### Community 1 - "STLM Utilities & Diagnostics"
Cohesion: 0.28
Nodes (14): always(), CHECK (validate STLM input parameters), cmprss(), compl(), CONVER (convergence test for Ritz values), dcheck(), deldup(), NEWPC (update Lanczos steps/copies optimization) (+6 more)

### Community 2 - "DVR Pipeline & Plots"
Cohesion: 0.27
Nodes (11): Optimized DVR / X-matrix pipeline (partially disabled), Packed lower-triangular storage in array W (H followed by overlap block), Sinc-DVR method (Colbert-Miller, J. Chem. Phys. 96, 1982 (1992)), PROGRAM DVR (main driver), DVR_OP (optimized DVR output to fort.11), DVR_PR (plot primitive DVR basis function), FBR (finite basis representation plot, stub), OPTMZD (build X coordinate matrix) (+3 more)

### Community 3 - "Shift-Invert Solver Driver"
Cohesion: 0.28
Nodes (9): SciPy eigensolver replacement (eigh subset_by_value / eigsh), Shift-invert Lanczos eigensolver (STLM library), DECOMP (LDL^T decomposition of K - mu*M), INITD (initialize STLM data and COMMON state), INITEW (initialize eigenvector workspace addresses), LDLT (LDL^T factorization dispatcher), REPL (update control variables between shifts), STLM (shift-invert Lanczos main driver) (+1 more)

### Community 4 - "Hamiltonian & Potential Physics"
Cohesion: 0.52
Nodes (7): Python conversion plan (port physics core, replace STLM with SciPy), Colbert-Miller particle-in-a-box kinetic energy matrix elements, Constantes2.txt (missing INCLUDE file: De, c1..c6, R0, r_m, J), Extended-Rydberg potential with centrifugal J(J+1)/(2*mu*r^2) term, Unit conversions: 219474.631 hartree to cm-1; 1822.88839 amu to electron mass, BUILD (assemble sinc-DVR Hamiltonian), V (extended-Rydberg potential function)

### Community 5 - "File I/O & Verification"
Cohesion: 0.47
Nodes (6): File I/O map (Fortran unit numbers = fort.N filenames), Regression verification against committed fort.4 / fort.97, fort.3 (input: grid A/B, NPOINT, AMASS; ELOW/EHIGH/NUMEIG), fort.4 (output: energies, spectrum, level differences, spectroscopic constants in cm-1), fort.97 (output: we in cm-1), fort.98 / fort.99 (V(x) evaluation log)

### Community 6 - "Eigenvalue Output & Spectroscopy"
Cohesion: 0.67
Nodes (4): Spectroscopic constants we, wexe, weye, alfae, gamae (finite differences of first 4 levels), EIGCLL (STLM wrapper and results writer), EXCHG (read eigenvectors from unit 19), fort.18 / fort.19 (STLM scratch: KFILE=18 sequential, DAFILE=19 direct-access)

## Ambiguous Edges - Review These
- `EIGCLL (STLM wrapper and results writer)` → `EXCHG (read eigenvectors from unit 19)`  [AMBIGUOUS]
  DVR.f · relation: calls

## Knowledge Gaps
- **2 isolated node(s):** `ranfx (RNG wrapper over rand())`, `REPL (update control variables between shifts)`
  These have ≤1 connection - possible missing edges or undocumented components.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **What is the exact relationship between `EIGCLL (STLM wrapper and results writer)` and `EXCHG (read eigenvectors from unit 19)`?**
  _Edge tagged AMBIGUOUS (relation: calls) - confidence is low._
- **Why does `STLM (shift-invert Lanczos main driver)` connect `Shift-Invert Solver Driver` to `STLM Lanczos Core`, `STLM Utilities & Diagnostics`, `DVR Pipeline & Plots`, `Hamiltonian & Potential Physics`, `Eigenvalue Output & Spectroscopy`?**
  _High betweenness centrality (0.424) - this node is a cross-community bridge._
- **Why does `EIGCLL (STLM wrapper and results writer)` connect `Eigenvalue Output & Spectroscopy` to `DVR Pipeline & Plots`, `Shift-Invert Solver Driver`, `Hamiltonian & Potential Physics`, `File I/O & Verification`?**
  _High betweenness centrality (0.238) - this node is a cross-community bridge._
- **Why does `Python conversion plan (port physics core, replace STLM with SciPy)` connect `Hamiltonian & Potential Physics` to `Shift-Invert Solver Driver`, `File I/O & Verification`, `Eigenvalue Output & Spectroscopy`?**
  _High betweenness centrality (0.178) - this node is a cross-community bridge._
- **Are the 2 inferred relationships involving `STLM (shift-invert Lanczos main driver)` (e.g. with `SciPy eigensolver replacement (eigh subset_by_value / eigsh)` and `Packed lower-triangular storage in array W (H followed by overlap block)`) actually correct?**
  _`STLM (shift-invert Lanczos main driver)` has 2 INFERRED edges - model-reasoned connections that need verification._
- **What connects `ranfx (RNG wrapper over rand())`, `REPL (update control variables between shifts)` to the rest of the system?**
  _2 weakly-connected nodes found - possible documentation gaps or missing edges._