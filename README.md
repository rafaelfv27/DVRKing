# DVRpy

Solver sinc-DVR para a equação de Schrödinger vibracional 1D de uma molécula diatômica
(Colbert & Miller, *J. Chem. Phys.* **96**, 1982, 1992).

Porte em Python/NumPy/SciPy do programa Fortran 77 `DVR.f` (J.J. Soares Neto, 1995), mantido
no repositório como referência. A biblioteca Lanczos STLM do original (~6300 linhas) foi
substituída por `scipy.linalg.eigh`.

## O que o script faz

Entrada: um CSV com duas colunas — distância internuclear `r` e energia potencial `V`.
Saída: níveis vibracionais, espectro e constantes espectroscópicas, em cm⁻¹.

Pipeline (`dvr.py`):

1. **Lê o CSV** e converte para bohr/hartree. As unidades são detectadas automaticamente
   (nome da coluna primeiro, magnitude dos dados depois) e podem ser forçadas com
   `--r-unit` / `--e-unit`.
2. **Ajusta um potencial Rydberg estendido de ordem 6** aos pontos, por projeção variável
   (mínimos quadrados lineares em `De, xe, c2..c6` dentro de um ajuste 2D em `Re, c1`).
   Devolve `De, Re, xe, c1..c6` e o rms do ajuste.
3. **Monta o Hamiltoniano DVR** numa grade uniforme de 500 pontos cobrindo o intervalo de `r`
   do CSV. Energia cinética pelas fórmulas analíticas de Colbert-Miller (partícula na caixa);
   potencial na diagonal.
4. **Diagonaliza** para J=0 e J=1 (o termo centrífugo `J(J+1)/2μr²` entra no potencial).
   Todos os estados ligados abaixo de `De` são retornados — o número de níveis é automático,
   não é um parâmetro. Níveis espúrios de "caixa" perto do limiar são descartados.
5. **Calcula as constantes**: `we, wexe, weye` das diferenças dos 4 primeiros níveis;
   `alfae, gamae` pela metodologia dos dois runs do legado (constantes de J=0 + espaçamentos
   de J=1); `Be` e `Bv` de `Bv = (E_v(J=1) − E_v(J=0))/2`.
6. **Escreve** `<nome>_out.txt` e, por padrão, um relatório LaTeX/PDF `<nome>_report.pdf`
   com três tabelas e quatro figuras (curva de energia potencial + ajuste, resíduo do ajuste,
   escada de níveis, roll-off dos espaçamentos) — `dvrreport.py`.

Tudo além da massa reduzida vem do CSV: a grade é o intervalo de `r`, o número de níveis são
os estados ligados abaixo do `De` ajustado.

## Como rodar

Requisitos: Python 3, `numpy`, `scipy`, `matplotlib` (só para o PDF) e `pdflatex` no PATH
(só para o PDF).

```bash
pip install numpy scipy matplotlib
```

Uso:

```bash
python dvr.py <arquivo.csv> (--mass-amu M | --mass M) [--r-unit U] [--e-unit U] [--no-pdf]
```

A massa reduzida é obrigatória, nas unidades que preferir:

- `--mass-amu` — em uma (convertida internamente por 1822.88839)
- `--mass` — em massas de elétron (unidades atômicas)

Exemplos:

```bash
# Li–C70, massa reduzida 6.884168 amu -> Li_Omega_out.txt + Li_Omega_report.pdf
python dvr.py Li_Omega.csv --mass-amu 6.884168

# só o texto, sem LaTeX
python dvr.py potential.csv --mass 218947.07152 --no-pdf

# forçando as unidades do CSV
python dvr.py meudado.csv --mass-amu 6.884168 --r-unit angstrom --e-unit eV
```

Formato do CSV: duas colunas, separadas por vírgula ou espaço, com ou sem cabeçalho.
Linhas de cabeçalho nomeando a unidade (`radii_bohr,E_hartree`) evitam a heurística de
detecção. Unidades aceitas: `r` em `bohr`/`angstrom`; `V` em `hartree`, `cm-1`, `eV`,
`kcal/mol`, `kJ/mol`.

## Testes

```bash
python -m pytest test_dvr.py
```

Regressão contra o run de referência do `DVR.exe` legado. As saídas `fort.*` daquele run não
estão mais no repositório: os 500 pontos (r, V) que ele registrou viraram `potential.csv`, e os
autovalores/we de referência estão fixados no próprio `test_dvr.py`. O núcleo DVR os reproduz
em 1e-6 cm⁻¹; o caminho completo (CSV → ajuste → níveis) em 0.1 cm⁻¹ — limitado pelo ajuste, já
que os dados legados desviam de um Rydberg-6 perfeito por rms 1.2e-7 hartree.

## Arquivos

| Arquivo | Conteúdo |
|---|---|
| `dvr.py` | solver: leitura do CSV, ajuste, Hamiltoniano, diagonalização, constantes |
| `dvrreport.py` | figuras e relatório LaTeX/PDF (sem física) |
| `test_dvr.py` | regressão vs. o legado |
| `DVR.f`, `DVR.exe` | programa Fortran original (referência) |
| `Li_Omega.csv` | curva Li–C70 (μ = 6.884168 amu) |
| `potential.csv` | curva do run de referência do legado |
| `wiki/` | documentação do código Fortran e do plano de porte |

## Notas

- Unidades internas atômicas; energias convertidas com 219474.631 (hartree → cm⁻¹), o literal
  do legado, mantido para a regressão.
- O legado imprimia os níveis com esse fator em precisão simples (219474.625) — o teste desfaz
  isso antes de comparar com os valores de referência. Ver `wiki/Bugs Legados.md`.
- `alfae`/`gamae` só são definidos para J ≠ 0; ver `nota_alfae_gamae.pdf`.
