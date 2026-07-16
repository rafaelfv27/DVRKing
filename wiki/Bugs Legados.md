---
title: Bugs Legados
tags:
  - dvr
  - bugs
  - fortran
---

# Bugs Legados

Defeitos no Fortran original. Decisão de porte: reproduzir saídas científicas, **não** reproduzir bugs — mas documentá-los para explicar qualquer diferença.

## Off-by-one no espectro

`DVR.f:303` e `DVR.f:332`: loops `I=1..NUMEV` leem `W(I+1)` ⇒ na última iteração leem `W(NUMEV+1)`, **além dos autovalores convergidos**. `W(NUMEV+1)` é lixo do workspace STLM. As últimas linhas de "ESPECTRO" e "DIFERENÇA" em `fort.4` são inválidas.
**Porte**: iterar até `NUMEV-1` ou calcular só transições válidas. ^off-by-one

## V avaliado no índice

`DVR.f:146-147` (DVR_OP): `TMPkcm = dble(I)` e `WRITE(11,*) XI, V(TMPkcm)` — escreve $V(i)$ (índice 1..500 como coordenada!) contra $x_i$. Deveria ser `V(XI)`. `fort.11` do repo é fisicamente sem sentido (V avaliado em r=1,2,...,500 bohr).
Efeito colateral: são essas 500 avaliações inteiras que alimentam o filtro `mod(XI,1.)≠0` de `fort.98`.
**Porte**: escrever `(x_i, V(x_i))` correto. ^v-indice

## alfae/gamae contaminados

`DVR.f:342-353`: `alfae` e `gamae` calculados com `we/wexe/weye` **hardcoded** (29.4345888..., de outro sistema/run antigo), misturados a diferenças de níveis do run atual. Fórmulas de Dunham exigiriam também dependência em J (rotacional) — o comentário `*WRITE(4,*) 'ENERGIA ROVIBRACIONAL J=1...'` sugere fluxo rovibracional desativado.
**Porte**: reproduzir para regressão com flag "legacy", sinalizar como não confiável. Comentário no código: `!! só imprimir essa` (a linha do WE). ^alfae

## EIGVCT nunca preenchido

`CALL EXCHG` comentado (`DVR.f:314`) ⇒ `EIGVCT` fica zerado. FBR, OPTMZD, DVR_OP operam sobre zeros. Todo o fluxo "DVR otimizado" (X matrix, quadratura otimizada) está morto desde que a 2ª chamada `EIGCLL` foi comentada (`DVR.f:76`). ^eigvct-zero

## Suspeita r_m cumulativa

`DVR.f:6839`: `r_m = r_m*1822.88839` dentro de V. Se `Constantes2.txt` define `r_m` via `DATA` (estática), a multiplicação **acumula a cada chamada** (V é chamada 1000×) e o termo centrífugo colapsa para ~0 depois de poucas chamadas. Se define via atribuição executável... INCLUDE não pode conter atribuição executável antes das declarações — precisa ver o arquivo para julgar. Com J=0 o termo é nulo e o bug fica invisível. **Verificar quando Constantes2.txt aparecer.** ^rm-cumulativo

## Menores

- `READ *` final (`DVR.f:90`) — trava esperando Enter; remover.
- `DIMENSION W(6)` em EIGCLL (`DVR.f:205`) — declaração mentirosa (array real é gigante); legal em F77, mas confunde.
- `IMPLICIT REAL*8` + variáveis I-N inteiras: `J` no potencial é inteiro implícito.
- Formato `200 FORMAT` órfão (`DVR.f:306`) — nunca referenciado.

Relacionado: [[EIGCLL - Diagonalizacao]], [[Funcao V - Potencial]], [[Rotinas Secundarias]].
