# DVRKing

Solver sinc-DVR para a equação de Schrödinger vibracional 1D de uma molécula diatômica
(Colbert & Miller, *J. Chem. Phys.* **96**, 1982, 1992), em Python/NumPy/SciPy.

Dá a ele uma curva de energia potencial e ele devolve os níveis vibracionais, o espectro e as
constantes espectroscópicas.

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
4. **Diagonaliza** (`scipy.linalg.eigh`) para J=0 e J=1 — o termo centrífugo `J(J+1)/2μr²`
   entra no potencial. Todos os estados ligados abaixo de `De` são retornados: o número de
   níveis é automático, não é um parâmetro. Níveis espúrios de "caixa" perto do limiar são
   descartados.
5. **Calcula as constantes**: `we, wexe, weye` das diferenças dos 4 primeiros níveis;
   `alfae, gamae` combinando as constantes de J=0 com os espaçamentos de J=1; `Be` e `Bv` de
   `Bv = (E_v(J=1) − E_v(J=0))/2`.
6. **Extrapola por Dunham** para J arbitrário (o DVR só diagonaliza em J=0 e J=1) e acha os
   dois limites de soma: `vmax` de `ε(v,0) < De` e `vmax` de `Bv > 0`.
7. **Soma a função de partição rovibracional** e daí a termodinâmica:

   ```
   Q(T) = Σ_v Σ_J (2J+1) exp(−[ε(v,J) − ε(0,0)]·hc/kT)
   ```

   Soma dupla explícita — sem Euler-Maclaurin, então nada assume T alta. Corta em
   `vmax` do `Bv > 0` (além dele ε cai com J e a soma diverge) e em J onde o termo já vale
   `e⁻⁴⁰`. Saem `Q, U, Cv, S, G` por mol, só rovibracional: sem translação, sem eletrônica,
   sem `pV` (logo `H = U` e `G = A`). Zero em `ε(0,0)`, então `Q → 1`, `U → 0`, `S → 0`
   quando `T → 0`.
8. **Escreve** `<nome>_out.txt` e, por padrão, um relatório LaTeX/PDF `<nome>_report.pdf`
   com cinco tabelas e quatro figuras (curva de energia potencial + ajuste, resíduo do ajuste,
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
python dvr.py <arquivo.csv> (--mass-amu M | --mass M) [--r-unit U] [--e-unit U]
                            [--jmax J] [--temps LISTA] [--no-pdf]
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

# temperaturas próprias na tabela termodinâmica (--temps '' pula a tabela)
python dvr.py Li_Omega.csv --mass-amu 6.884168 --temps 273.15,298.15,500,1000
```

Formato do CSV: duas colunas, separadas por vírgula ou espaço, com ou sem cabeçalho.
Linhas de cabeçalho nomeando a unidade (`radii_bohr,E_hartree`) evitam a heurística de
detecção. Unidades aceitas: `r` em `bohr`/`angstrom`; `V` em `hartree`, `cm-1`, `eV`,
`kcal/mol`, `kJ/mol`.

## Arquivos

| Arquivo | Conteúdo |
|---|---|
| `dvr.py` | solver: leitura do CSV, ajuste, Hamiltoniano, diagonalização, constantes, termodinâmica |
| `dvrreport.py` | figuras e relatório LaTeX/PDF (sem física) |
| `test_thermo.py` | checagens de `dvr.thermo` — `python3 test_thermo.py` |
## Notas

- Unidades internas atômicas; energias convertidas com 219474.631 (hartree → cm⁻¹) e massas
  com 1822.88839 (amu → massas de elétron).
