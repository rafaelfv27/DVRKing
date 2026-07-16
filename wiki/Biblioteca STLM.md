---
title: Biblioteca STLM (Lanczos Shift-Invert)
tags:
  - dvr
  - fortran
  - stlm
  - nao-portar
source: DVR.f:549-6820
---

# Biblioteca STLM — **NÃO PORTAR**

`DVR.f:549` (`ALLOC`) até `DVR.f:6820` (`WRINFO`), ~6300 linhas. Solver de autovalores **S**pectral **T**ransformation **L**anczos **M**ethod (Ericsson & Ruhe) para problema generalizado $Kx = \lambda Mx$ com matrizes em perfil skyline, janela espectral $[a,b]$, shifts múltiplos e restart. Tudo isso vira **uma chamada** `scipy.linalg.eigh`.

## Entrada principal

`SUBROUTINE STLM` (`DVR.f:5375`): loop principal — decomposição $LDL^T$ de $K-\mu M$ (`DECOMP`), iteração Lanczos (`LANCZO`), tridiagonal QL (`IMTQL2`), seleção de shifts (`NEWMU`), término (`TERMIN`). Estado global em ~25 blocos `COMMON /STLM**/`.

## Inventário de rotinas (para auditoria de completude do porte)

| Grupo | Rotinas | Função |
|---|---|---|
| Gerência de memória/ID | `ALLOC`, `FREEID`, `INITAD`, `INITEW`, `IO`, `REPL`, `UI` | vetores nomeados dentro de W, paging para fort.19 |
| Fatoração | `DECOMP`, `LDLT`, `LDL2`, `LDL3`, `LDLSUB` | $LDL^T$ skyline de $K-\mu M$ |
| Lanczos | `LANCZO`, `FRSTIT`, `INITLA`, `ALWAYS`, `SELDOM`, `TRIDIG` | recorrência de três termos, reortogonalização |
| Autovalores tridiagonais | `IMTQL2`, `CONVER`, `COMPL` | QL implícito, conversão λ=μ+1/ν |
| Shifts/janela | `NEWMU`, `NEWPC`, `HALF`, `CHECK`, `SAVEXL`, `ROTATE`, `DELDUP`, `CMPRSS` | estratégia de shift, dedup de autovalores |
| Produtos matriz-vetor | `OPM`, `OPM2`, `OPM3`, `OPMSUB`, `OPINV`, `PREPSV`, `MULVEC`, `MUL3`, `MULV`, `SUBVEC`, `SUB3`, `SUBV`, `TRANSF`, `MAXNRM` | y=Mx, y=(K−μM)⁻¹x, álgebra vetorial |
| Aleatórios | `RANDVC`, `RAN3`, `RANNUM`, `ranfx` | vetor inicial aleatório |
| Diagnóstico | `DCHECK`, `MCHECK`, `ERROR`, `WRINFO`, `WRINF3`, `WRINF4`, `TERMIN` | validação, mensagens, término |
| Relógio | `second` (usa `etime`) | timing |
| Utilidade DVR | `EXCHG` (`DVR.f:406`) | lê autovetores de fort.19 — chamada comentada |

## Equivalência SciPy

| STLM | SciPy |
|---|---|
| janela $[a,b]$ + shifts | `eigh(H, subset_by_value=(ELOW, EHIGH))` |
| Lanczos shift-invert | `scipy.sparse.linalg.eigsh(H, k, sigma=...)` (só se N ≫ 1000) |
| skyline + scratch fort.18/19 | desnecessário — denso em RAM (N=500 ⇒ 2 MB) |
| `TCONV`, `ERRNO`, restarts | desnecessário — eigh é direto e determinístico |

> [!success] Justificativa
> N ≤ 1000 (NDIM). Matriz densa 1000×1000 = 8 MB. `eigh` (LAPACK dsyevr) resolve em ms com garantia de todos os autovalores na janela — elimina risco de nível perdido do Lanczos (o próprio STLM tem lógica `SELDOM`/`NLEFT` para autovalores faltantes).

Relacionado: [[EIGCLL - Diagonalizacao]], [[Plano de Porte Python]].
