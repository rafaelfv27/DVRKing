---
title: BUILD — Construção do Hamiltoniano
tags:
  - dvr
  - fortran
  - hamiltoniano
  - colbert-miller
source: DVR.f:98-127
---

# `SUBROUTINE BUILD(A,B,NPOINT,AMASS,W,IW)`

Localização: `DVR.f:98-127`. Monta o Hamiltoniano sinc-DVR em **armazenamento triangular inferior empacotado** no array `W`, com a matriz de overlap (identidade) na segunda metade.

## Grade

Grade uniforme de partícula-na-caixa em $[A, B]$, **excluindo as bordas**:

$$x_i = A + (B-A)\,\frac{i}{N+1}, \qquad i = 1,\dots,N$$

com $N = \text{NPOINT}$. Espaçamento $\Delta x = (B-A)/(N+1)$.

## Fatores pré-computados

| Fortran | Fórmula | Papel |
|---|---|---|
| `FCT1` | $\pi^2 / [4\,\mu\,(B-A)^2]$ | escala cinética ($\hbar=1$, u.a.) |
| `FCT2` | $\pi / [2(N+1)]$ | argumento dos senos off-diagonal |
| `FCT3` | $\pi / (N+1)$ | argumento do seno diagonal |

`AMASS` = massa reduzida $\mu$ **já em massas de elétron** (218947.07152 no `fort.3`).

## Elementos de matriz (Colbert-Miller, caixa $[a,b]$, Apêndice A)

Off-diagonal ($i \neq j$), `DVR.f:113-114`:

$$T_{ij} = (-1)^{i-j}\,\frac{\pi^2}{4\mu(B-A)^2}\left[\frac{1}{\sin^2\!\frac{(i-j)\pi}{2(N+1)}} - \frac{1}{\sin^2\!\frac{(i+j)\pi}{2(N+1)}}\right]$$

Diagonal ($i = j$), `DVR.f:119-120`:

$$H_{ii} = \frac{\pi^2}{4\mu(B-A)^2}\left[\frac{2(N+1)^2+1}{3} - \frac{1}{\sin^2\!\frac{i\pi}{N+1}}\right] + V(x_i)$$

O potencial [[Funcao V - Potencial]] entra **só na diagonal** — essência do DVR.

## Layout de memória (crítico para entender o legado)

- `W(1 .. N(N+1)/2)` — triangular inferior de H, ordem linha a linha: `(1,1),(2,1),(2,2),(3,1)...`
- `W(N(N+1)/2+1 .. N(N+1))` — matriz de overlap M empacotada: 0 fora da diagonal, 1 na diagonal (identidade ⇒ problema de autovalor padrão, mas STLM trata como generalizado com `DIAGM=.FALSE.`)
- `IW(I) = I*(I+1)/2` — ponteiros de perfil para a diagonal (formato skyline do STLM)
- Resto de `W` — workspace do STLM (`MAXW=820000`)

> [!tip] Porte Python
> Nada disso sobrevive. Matriz densa `H = np.zeros((N, N))` direta:
> ```python
> i, j = np.indices((N, N)) + 1
> fct1 = np.pi**2 / (4 * mu * (B - A)**2)
> T = (-1.0)**(i - j) * fct1 * (
>     1/np.sin((i - j)*np.pi/(2*(N + 1)))**2
>     - 1/np.sin((i + j)*np.pi/(2*(N + 1)))**2)   # dá inf na diagonal — sobrescrever
> d = np.arange(1, N + 1)
> T[np.diag_indices(N)] = fct1 * ((2*(N + 1)**2 + 1)/3
>     - 1/np.sin(d*np.pi/(N + 1))**2)
> x = A + (B - A)*d/(N + 1)
> H = T + np.diag(Vfunc(x))
> ```
> Overlap/empacotamento/IW: descartar.

Relacionado: [[Pipeline Principal]], [[EIGCLL - Diagonalizacao]], [[Unidades e Constantes]].
