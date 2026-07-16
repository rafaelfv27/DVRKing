---
title: Função V — Potencial Rydberg Estendido
tags:
  - dvr
  - fortran
  - potencial
  - rydberg
source: DVR.f:6826-6854
---

# `REAL*8 FUNCTION V(XI)`

Localização: `DVR.f:6826-6854` (fim do arquivo). Potencial de Rydberg estendido (ordem 6) + termo centrífugo rotacional.

## Forma funcional

$$V(r) = -D_e\left(1 + \sum_{k=1}^{6} c_k\,(r-R_0)^k\right)e^{-c_1 (r-R_0)} + D_e + \frac{J(J+1)}{2\,\mu\,r^2}$$

Código (`DVR.f:6842-6846`):
```fortran
V=-De*(1+c1*xx+c2*xx**2+c3*xx**3+c4*xx**4+c5*xx**5+c6*xx**6)*
 >exp(-c1*xx)+ De
 >+(J*(J+1.0d0))/(2*r_m*(XI**2))
```
com `xx = XI - XE`, `XE = R0`.

> [!danger] Bloqueador do porte
> `INCLUDE 'Constantes2.txt'` (`DVR.f:6837`) — **arquivo ausente do repositório**. Contém `De, c1..c6, R0, r_m, J`. Sem ele não dá para recompilar nem portar V. Pedir ao usuário ou recuperar.
> Parcialmente recuperável: `fort.11`/`fort.98`/`fort.99` do último run contêm pares $(r, V(r))$ — um ajuste dos 9 parâmetros à curva logada é possível se o arquivo não aparecer.

## Detalhes traiçoeiros

1. **Conversão de massa embutida**: `r_m = r_m*1822.88839` (`DVR.f:6839`) — `r_m` vem em amu do include e é convertido para massas de elétron **a cada chamada**. Como `r_m` é variável de include (não SAVE explícito), em Fortran 77 o comportamento entre chamadas depende do compilador — na prática gfortran re-executa o INCLUDE-init? Não: INCLUDE só injeta declarações/DATA; se `r_m` vier de `DATA`, a multiplicação acumula a cada chamada **se** a variável for estática. Ver [[Bugs Legados#Suspeita r_m cumulativa]].
2. **J é inteiro implícito**: `J` cai na regra `IMPLICIT` I-N inteiro. `J(J+1.0d0)` mistura int e real.
3. **Logging em toda avaliação**: escreve em `fort.99` (sempre) e `fort.98` (só quando `mod(XI,1.) ≠ 0`). N chamadas de BUILD ⇒ N linhas. No porte: eliminar ou tornar opcional.
4. **Termo centrífugo usa r absoluto** (`XI**2`), não deslocado.

## Porte Python

```python
def V(r, p):  # p = dataclass com De, c[1..6], R0, mu_amu, J
    xx = r - p.R0
    poly = 1 + sum(p.c[k] * xx**k for k in range(1, 7))
    mu = p.mu_amu * 1822.88839
    return -p.De * poly * np.exp(-p.c[1] * xx) + p.De \
           + p.J * (p.J + 1) / (2 * mu * r**2)
```
Parâmetros via arquivo de config (JSON/TOML) substituindo `Constantes2.txt`.

Relacionado: [[BUILD - Hamiltoniano]], [[Mapa de IO]], [[Unidades e Constantes]], [[Bugs Legados]].
