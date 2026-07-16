---
title: Pipeline Principal (PROGRAM DVR)
tags:
  - dvr
  - fortran
  - pipeline
source: DVR.f:5-91
---

# Pipeline Principal — `PROGRAM DVR`

Localização: `DVR.f:5-91`. Programa principal, sequência fixa de 7 chamadas.

## Declarações-chave

| Símbolo | Valor/Forma | Papel |
|---|---|---|
| `NDIM` | `PARAMETER (NDIM=1000)` | dimensão máxima da grade |
| `W` | `DIMENSION W(3*NDIM*(NDIM+1)/2)` | array de trabalho gigante: H empacotado + overlap + workspace STLM |
| `IW` | `DIMENSION IW(NDIM)` | índices diagonais do perfil (`IW(I)=I*(I+1)/2`) |
| `EIGVCT` | `(NDIM,NDIM)` | autovetores (nunca preenchido de fato — ver [[Bugs Legados]]) |
| `XVCT` | `(NDIM,NDIM)` | autovetores da matriz X (fluxo desativado) |
| `IMPLICIT REAL*8 (A-H,O-Z)` | — | tipagem implícita em todo o código |

## Sequência de execução

1. `CALL READ_I(A,B,NPOINT,AMASS,ELOW,EHIGH,NUMEIG)` — lê 2 linhas de `fort.3` ([[Rotinas Secundarias#READ_I]])
2. `CALL BUILD(A,B,NPOINT,AMASS,W,IW)` — monta Hamiltoniano ([[BUILD - Hamiltoniano]])
3. `CALL DVR_PR(A,B,NPOINT)` — plota função DVR primitiva em `fort.10`
4. `CALL EIGCLL(ELOW,EHIGH,NUMEIG,W,IW,NPOINT,NDIM,EIGVCT,NUMEV)` — diagonaliza ([[EIGCLL - Diagonalizacao]])
5. `CALL FBR(A,B,NPOINT,NDIM,EIGVCT)` — efetivamente no-op (print apenas)
6. `CALL OPTMZD(A,B,NDIM,EIGVCT,NPOINT,NUMEV,W,IW)` — monta matriz X ⟨i|x|j⟩
7. ~~`CALL EIGCLL(...XVCT...)`~~ — **comentado** (linha 76): a diagonalização da matriz X está desativada
8. `CALL DVR_OP(A,B,NPOINT,NUMEV,W,NDIM,EIGVCT,XVCT)` — escreve `fort.11`

Fecha com `cpu_time`, imprime `tempo`, e `READ *` (espera Enter — remover no porte).

> [!important] O que realmente importa
> Só os passos 1, 2 e 4 produzem os resultados científicos (`fort.4`, `fort.97`). Passos 3, 5, 6, 8 são plots/diagnósticos com o fluxo de DVR otimizado **parcialmente desativado** (passo 7 comentado ⇒ `XVCT` nunca é preenchido; `NUMEV` vem de EIGCLL).

## Porte Python

```python
def main():
    params = read_input()          # READ_I
    H = build_hamiltonian(params)  # BUILD (denso, simétrico)
    E, C = diagonalize(H, params)  # EIGCLL → scipy.linalg.eigh
    write_outputs(E, params)       # fort.4 / fort.97
```

Ver [[Plano de Porte Python]].
