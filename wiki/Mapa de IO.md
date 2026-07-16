---
title: Mapa de I/O (arquivos fort.N)
tags:
  - dvr
  - fortran
  - io
---

# Mapa de I/O

Fortran usa unidades numéricas ⇒ arquivos literais `fort.N` no diretório de execução.

| Arquivo | Direção | Escrito por | Conteúdo |
|---|---|---|---|
| `fort.3` | entrada | — (lido por READ_I) | linha 1: `A B NPOINT AMASS`; linha 2: `ELOW EHIGH NUMEIG` |
| `fort.4` | saída | EIGCLL | energias, espectro, diferenças, constantes espectroscópicas (cm⁻¹) — **arquivo de regressão principal** |
| `fort.97` | saída | EIGCLL | só `we` em cm⁻¹ (1 linha) |
| `fort.10` | saída | DVR_PR | plot da função DVR primitiva `(x, φ₁(x))` |
| `fort.11` | saída | DVR_OP | pares `(xᵢ, V(i))` — bugado, ver [[Bugs Legados#V avaliado no índice]] |
| `fort.98` | saída | V | log de V quando `mod(XI,1.)≠0` (500 linhas no último run) |
| `fort.99` | saída | V | log de **toda** avaliação de V (1000 linhas: BUILD 500 + DVR_OP 500) |
| `fort.18` | scratch | STLM (KFILE) | sequencial unformatted — descartar no porte |
| `fort.19` | scratch | STLM (DAFILE) | direct access RECL=8·NPOINT, autovetores — descartar no porte |

## Valores de regressão (run committado, NPOINT=500)

- `fort.4`: nível 1 = **9.7211221261708740** cm⁻¹; nível 2 = 28.870890234651771; nível 3 = 47.635691885033019
- `fort.97`: we = **19.530919626995409** cm⁻¹
- Tolerância alvo do porte: ~1e-6 cm⁻¹

> [!tip] Porte
> Manter valores/unidades de `fort.4` e `fort.97` byte-comparáveis numericamente; nomes de arquivo podem modernizar (`energies.txt`, config por CLI). Logs `fort.98/99`: eliminar ou flag `--debug`.

Relacionado: [[EIGCLL - Diagonalizacao]], [[Plano de Porte Python]].
