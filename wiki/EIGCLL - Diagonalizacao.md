---
title: EIGCLL — Diagonalização e Saídas
tags:
  - dvr
  - fortran
  - lanczos
  - espectroscopia
source: DVR.f:190-403
---

# `SUBROUTINE EIGCLL(ELOW,EHIGH,NUMEIG,W,IW,NPOINT,NDIM,EIGVCT,NUMEV)`

Localização: `DVR.f:190-403`. Wrapper do STLM (escrito em 1991) + toda a escrita de resultados científicos. É a rotina que o porte substitui por `scipy.linalg.eigh`.

## Parte 1 — Setup STLM (`DVR.f:224-258`)

| Parâmetro | Valor | Significado |
|---|---|---|
| `A, B` | `ELOW, EHIGH` | janela de autovalores em hartree (0.0–1000.0 no input atual = "todos") |
| `DIAGM` | `.FALSE.` | overlap não é tratada como diagonal trivial |
| `MAXW/MAXIW` | 820000 / 16000 | tamanhos de workspace |
| `MAXL` | `2*NUMEIG` | máximo de passos Lanczos |
| `PROFIL=1, PMAX=50, MXREST=0, MSGLVL=2` | — | perfil skyline, shifts, restarts, verbosidade |
| `DAFILE=19, KFILE=18` | scratch | direct-access (RECL=`8*NPOINT`) e sequencial |

Abre `fort.19`/`fort.18` como `STATUS='SCRATCH'`, chama:

```fortran
CALL STLM(NPOINT,A,B,MAXL,PROFIL,PMAX,MXREST,MSGLVL,MAXW,
*         MAXIW,DAFILE,MAXREC,KFILE,X,BG,TCONV,NLEFT,ERRNO,W,IW)
```

Retorno: autovalores convergidos **nas primeiras posições de `W`** (sobrescreve o H!), `TCONV` = nº convergidos. `NUMEV = MIN(NUMEIG, TCONV)`.

> [!warning] W é reciclado
> Depois de STLM, `W(1..NUMEV)` são os autovalores em hartree. O Hamiltoniano foi destruído. `EIGVCT` **nunca é preenchido** — a chamada `EXCHG` que leria os vetores de `fort.19` está comentada (`DVR.f:314`).

## Parte 2 — Saídas em `fort.4` (cm⁻¹, fator 219474.631)

1. **ENERGIA VIBRACIONAL** (`DVR.f:277-285`): `219474.631*W(I)`, I=1..NUMEV
2. **ESPECTRO VIBRACIONAL** (`DVR.f:301-317`): `219474.631*(W(I+1)-W(1))` — transições a partir do fundamental. ⚠️ lê `W(NUMEV+1)` — ver [[Bugs Legados#Off-by-one no espectro]]
3. **DIFERENÇA ENTRE NÍVEIS** (`DVR.f:328-338`): `219474.631*(W(I+1)-W(I))` — mesmo off-by-one

## Parte 3 — Constantes espectroscópicas (`DVR.f:342-392`)

Duas fases, com pegadinha:

**Fase A** (`DVR.f:342-353`): usa valores **hardcoded** `we=29.4345888...`, `wexe=1.9986e-2`, `weye=2.1474e-5` para calcular:
- `alfae = ⅛[−12(E₂−E₁) + 4(E₃−E₁)]·219474.631 + 4we − 23weye` *(mistura hardcoded + calculado)*
- `gamae = ¼[−2(E₂−E₁) + (E₃−E₁)]·219474.631 + 2wexe − 9weye`

**Fase B** (`DVR.f:356-364`): **recalcula** we, wexe, weye por diferenças finitas dos 4 primeiros níveis (em hartree):
- `we   = 1/24·[141(E₂−E₁) − 93(E₃−E₁) + 23(E₄−E₁)]`
- `wexe = 1/4·[13(E₂−E₁) − 11(E₃−E₁) + 3(E₄−E₁)]`
- `weye = 1/6·[3(E₂−E₁) − 3(E₃−E₁) + (E₄−E₁)]`

`fort.97` recebe `we*219474.631` (Fase B). `fort.4` imprime we/wexe/weye da Fase B (convertidos) e alfae/gamae da Fase A (já em cm⁻¹, contaminados pelos hardcoded). Comentário no código: `'WE(cm-1)=' ... !! só imprimir essa`.

> [!bug] alfae/gamae não são confiáveis
> Dependem de constantes hardcoded de outro sistema molecular. No porte: reproduzir para regressão, mas sinalizar. Ver [[Bugs Legados]].

## Porte Python

```python
E = scipy.linalg.eigh(H, eigvals_only=True,
                      subset_by_value=(ELOW, EHIGH))[:NUMEIG]
cm = 219474.631
we   = (141*(E[1]-E[0]) - 93*(E[2]-E[0]) + 23*(E[3]-E[0])) / 24 * cm
wexe = ( 13*(E[1]-E[0]) - 11*(E[2]-E[0]) +  3*(E[3]-E[0])) / 4  * cm
weye = (  3*(E[1]-E[0]) -  3*(E[2]-E[0]) +     (E[3]-E[0])) / 6 * cm
```

Relacionado: [[Biblioteca STLM]], [[Mapa de IO]], [[Plano de Porte Python]].
