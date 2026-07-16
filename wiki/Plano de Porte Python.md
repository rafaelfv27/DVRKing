---
title: Plano de Porte Python
tags:
  - dvr
  - python-port
  - plano
status: pronto-para-implementar
---

# Plano de Porte Python

## Arquitetura alvo

Um arquivo `dvr.py` (~150 linhas) + config de potencial. Sem classes desnecessárias.

```
dvr.py
  read_input(path)        ← READ_I         (fort.3 ou argparse)
  V(r, params)            ← FUNCTION V     (params de constantes.toml)
  build_hamiltonian(...)  ← BUILD          (denso, NumPy vetorizado)
  diagonalize(H, window)  ← EIGCLL+STLM    (scipy.linalg.eigh)
  spectro_constants(E)    ← DVR.f:356-364  (diferenças finitas)
  write_outputs(...)      ← fort.4/fort.97
constantes.toml           ← Constantes2.txt (De, c1..c6, R0, r_m, J)
```

## Mapeamento 1:1

| Fortran | Python | Nota |
|---|---|---|
| `W` empacotado + `IW` skyline | `np.ndarray (N,N)` | [[BUILD - Hamiltoniano]] tem o snippet vetorizado |
| STLM (6300 linhas) | `scipy.linalg.eigh(H, subset_by_value=(ELOW,EHIGH))` | [[Biblioteca STLM]] |
| fort.18/19 scratch | — | eliminado |
| `EXCHG` | `eigh` já retorna autovetores | |
| `RANNUM`/`ranfx`/`second` | — | eliminado |
| logs fort.98/99 | flag `--debug` ou nada | |
| `READ *` final | — | eliminado |

## Bloqueador

> [!danger] Constantes2.txt ausente
> `De, c1..c6, R0, r_m, J` desconhecidos. Opções: (1) usuário fornece; (2) ajustar aos pares (r,V) de `fort.11`/`fort.98`/`fort.99` — atenção: `fort.11` tem V(índice) bugado, usar `fort.99` filtrado (ver [[Bugs Legados]]); (3) recompilar impossível sem ele.
> Tudo exceto V é portável **agora** — dá para validar cinética + eigh com potencial dummy lido de fort.99 interpolado.

## Verificação (regressão)

1. Input: A=13.04, B=30.20, NPOINT=500, AMASS=218947.07152, janela 0–1000 h, NUMEIG=10
2. Comparar com `fort.4`: E₁ = 9.7211221261708740 cm⁻¹ (tolerância 1e-6)
3. Comparar `fort.97`: we = 19.530919626995409 cm⁻¹
4. Espectro/diferenças: comparar só até NUMEV−1 (off-by-one legado, [[Bugs Legados#Off-by-one no espectro]])
5. alfae/gamae: reproduzir com os mesmos hardcoded se quiser diff exato

## Decisões registradas

- Constantes de conversão **literais do legado** (219474.631, 1822.88839) — não CODATA ([[Unidades e Constantes]])
- Fluxo DVR otimizado (FBR/OPTMZD/DVR_OP) morto no legado ⇒ não portar ([[Rotinas Secundarias]])
- Bugs legados corrigidos no porte, documentados em [[Bugs Legados]] para explicar diffs

Relacionado: [[Home]], [[Pipeline Principal]].
