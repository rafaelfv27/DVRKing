---
title: Rotinas Secundárias
tags:
  - dvr
  - fortran
---

# Rotinas Secundárias

## READ_I

`DVR.f:538-547`. Lê `fort.3`, formato livre:
```
linha 1: A B NPOINT AMASS     → 13.04 30.20 500 218947.07152
linha 2: ELOW EHIGH NUMEIG    → 0.0 1000.0 10
```
- `A, B` — limites da grade em bohr
- `NPOINT` — nº de pontos (≤ NDIM=1000)
- `AMASS` — massa reduzida **em massas de elétron** (não amu!)
- `ELOW, EHIGH` — janela de autovalores em hartree
- `NUMEIG` — nº de autovalores desejados

`fort.3` tem linhas extras de comentário após as 2 lidas — ignoradas pelo Fortran; o leitor Python deve ignorar também.

## DVR_PR

`DVR.f:159-186`. Escreve em `fort.10` a função de base DVR primitiva centrada no 1º ponto da grade ($x_1$):

$$\phi_1(x) = \frac{2}{N+1}\sum_{n=1}^{N} \sin\!\frac{n\pi(x-A)}{B-A}\,\sin\!\frac{n\pi(x_1-A)}{B-A}$$

Varre $x$ de A até B com passo $\Delta x$. Só visualização — sem efeito no cálculo. Loop `GOTO 100` = while.

## FBR

`DVR.f:420-433`. **No-op efetivo**: loop sobre a grade com o `PRINT` interno comentado. Só imprime `'FBR'` no console. Não portar (ou portar como stub).

## OPTMZD

`DVR.f:437-463`. Monta matriz de posição $X_{ij} = \sum_k C_{ik}\,x_k\,C_{jk}$ (elementos ⟨ψᵢ|x|ψⱼ⟩ na base de autovetores) em armazenamento empacotado + overlap identidade. **Inútil no fluxo atual**: `EIGVCT` está zerado (EXCHG comentado) ⇒ X = 0, e a diagonalização de X (passo 7 do [[Pipeline Principal]]) está comentada. Parte do fluxo "DVR otimizado" abandonado.

## DVR_OP

`DVR.f:131-152`. Deveria plotar o DVR otimizado; na prática escreve em `fort.11` pares `(XI, V(dble(I)))` — ver [[Bugs Legados#V avaliado no índice]]. O loop interno com `EIGVCT·XVCT` calcula `FXI` e **nunca usa**.

## EXCHG

`DVR.f:406-416`. Leria autovetor I de `fort.19` (direct access) para `EIGVCT(I,:)`. Chamada comentada em `DVR.f:314`. No porte: desnecessário — `eigh` retorna vetores diretamente.

> [!note] Resumo de porte
> Portar de verdade: `READ_I` (parser de input), `DVR_PR` (se quiser o plot). Stub/descartar: `FBR`, `OPTMZD`, `DVR_OP`, `EXCHG`. Se o usuário quiser o fluxo de DVR otimizado funcional, é *feature nova* — no legado está morto.

Relacionado: [[Pipeline Principal]], [[Mapa de IO]].
