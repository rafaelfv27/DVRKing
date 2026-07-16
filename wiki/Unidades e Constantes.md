---
title: Unidades e Constantes
tags:
  - dvr
  - unidades
---

# Unidades e Constantes

Internamente tudo em **unidades atômicas** (ħ=1, mₑ=1, bohr, hartree).

| Constante | Valor no código | Conversão | Onde |
|---|---|---|---|
| hartree → cm⁻¹ | `219474.631` | energias de saída | `DVR.f:281,303,332,347-353,366,378-390` |
| amu → mₑ | `1822.88839` | massa reduzida dentro de V | `DVR.f:6839` |
| π | `3.1415926535898D0` | PARAMETER em BUILD e DVR_PR | `DVR.f:102,163` |

## Duas massas distintas — não confundir

1. `AMASS` (fort.3, = 218947.07152) — massa reduzida da **energia cinética**, já em mₑ. Equivale a ≈ 120.11 amu (218947.07152/1822.88839).
2. `r_m` (Constantes2.txt) — massa reduzida do **termo centrífugo** em V, em amu, convertida ×1822.88839 dentro de V.

Fisicamente deveriam ser a mesma molécula; no código são fontes independentes (risco de inconsistência no input legado).

## Grade do run de referência

- A=13.04, B=30.20 bohr; NPOINT=500 ⇒ Δx = 17.16/501 ≈ 0.034252 bohr
- Janela: ELOW=0, EHIGH=1000 hartree (pega tudo), NUMEIG=10

> [!tip] Porte
> Usar `scipy.constants`? Não — **manter os valores literais do legado** (219474.631, 1822.88839) para reproduzir `fort.4` a 1e-6 cm⁻¹. Valores CODATA modernos diferem no 4º decimal e quebram a regressão.

Relacionado: [[BUILD - Hamiltoniano]], [[Funcao V - Potencial]].
