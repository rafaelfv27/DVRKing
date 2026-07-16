# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project goal

Convert the legacy Fortran 77 program `DVR.f` to Python. The physics core is small (~500 lines); the remaining ~6300 lines are the STLM Lanczos eigensolver library, which must NOT be ported — replace it entirely with NumPy/SciPy.

## What the program does

Sinc-DVR (Discrete Variable Representation) solver for the 1D vibrational Schrödinger equation of a diatomic molecule, per Colbert & Miller, J. Chem. Phys. 96, 1982 (1992). Written by J.J. Soares Neto (1995). Pipeline in `PROGRAM DVR` (DVR.f:5):

1. `READ_I` (DVR.f:538) — reads input from `fort.3`
2. `BUILD` (DVR.f:98) — builds the DVR Hamiltonian in packed lower-triangular storage (array `W`), kinetic term analytic (Colbert-Miller particle-in-a-box formulas), potential `V(x)` on diagonal
3. `DVR_PR` (DVR.f:159) — writes primitive DVR basis function to `fort.10` (plot)
4. `EIGCLL` (DVR.f:190) — wrapper around STLM shift-invert Lanczos; diagonalizes H, writes vibrational energies, spectrum, level differences, and spectroscopic constants (we, wexe, weye, alfae, gamae — finite-difference formulas from the first 4 levels, DVR.f:342-392) to `fort.4`, and we to `fort.97`
5. `FBR`, `OPTMZD`, `DVR_OP` — build X matrix / optimized DVR (partially commented out; `DVR_OP` writes `fort.11`)

Potential `V(XI)` (DVR.f:6826): extended-Rydberg-style form `-De*(1+c1*xx+...+c6*xx^6)*exp(-c1*xx) + De + J(J+1)/(2*mu*r^2)` (centrifugal term for rotation). Parameters (`De, c1..c6, R0, r_m, J`) come from `INCLUDE 'Constantes2.txt'` — **this file is missing from the repo** and is required to recompile or to port V(x). Ask the user for it or recover it before converting the potential. V also logs every evaluation to `fort.98`/`fort.99`.

Unit conventions: atomic units internally; energies converted with 219474.631 (hartree → cm⁻¹); mass converted with 1822.88839 (amu → electron masses, applied inside V to `r_m`).

## File I/O map (Fortran unit numbers = literal filenames `fort.N`)

| File | Direction | Content |
|------|-----------|---------|
| fort.3 | input | line 1: `A B NPOINT AMASS` (grid start/end in bohr, # points, reduced mass); line 2: `ELOW EHIGH NUMEIG` (eigenvalue window in hartree, # of eigenvalues) |
| fort.4 | output | energies, spectrum, differences, spectroscopic constants (cm⁻¹) |
| fort.97 | output | we (cm⁻¹) |
| fort.10 | output | primitive DVR function plot |
| fort.11 | output | grid x vs V (DVR_OP) |
| fort.98/99 | output | V(x) evaluation log |
| fort.18/19 | scratch | STLM internal (drop in Python) |

## Legacy Fortran commands

```
gfortran -o DVR.exe DVR.f     # requires Constantes2.txt present (INCLUDE)
./DVR.exe                      # reads fort.3 from cwd; ends with READ * (press Enter)
```

`fort.4`, `fort.97`, etc. from the last run are in the repo — use them as regression references for the Python port (energies should match to ~1e-6 cm⁻¹).

## Python conversion plan

- **NumPy** for the Hamiltonian: build dense symmetric matrix directly (no packed triangular storage — NPOINT ≤ 1000, dense is trivial).
- **SciPy** replaces all of STLM: `scipy.linalg.eigh(H, subset_by_value=(ELOW, EHIGH))` or `subset_by_index` reproduces the eigenvalue-window behavior; `scipy.sparse.linalg.eigsh` (Lanczos, same algorithm family as STLM) only if N grows large.
- Keep the Colbert-Miller kinetic matrix elements exactly as in `BUILD` (uniform grid, particle-in-a-box sinc-DVR).
- Potential parameters: load from a config file replacing `Constantes2.txt`.
- Preserve output values/units of `fort.4` (energies, spectrum, spectroscopic constants) for comparison; file naming can be modernized.
- Do not port: STLM library (everything from `ALLOC` at DVR.f:549 onward except `V` at the end), `RANNUM`, `ranfx`, `second`, `EXCHG`, scratch-file logic.

## Verification

Run Python port with the same `fort.3` input (A=13.04, B=30.20, NPOINT=500, AMASS=218947.07152) and diff eigenvalues against the committed `fort.4` (first level ≈ 9.7211221 cm⁻¹) and `fort.97` (we ≈ 19.5309196 cm⁻¹).
