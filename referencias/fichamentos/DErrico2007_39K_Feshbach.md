# D'Errico et al. 2007 — Feshbach resonances in ultracold ³⁹K

**Referência:** New J. Phys. **9**, 223 (2007)
**DOI:** 10.1088/1367-2630/9/7/223 · **Acesso:** aberto

## Em uma frase
Primeiro levantamento sistemático das ressonâncias de Feshbach do ³⁹K
(8 ressonâncias, l = 0), com posições experimentais e cálculo de canais
acoplados.

## Convenções
a(B) = a_bg[1 − Δ/(B − B₀)] · C5. Larguras tabeladas como **−Δ**;
atenção ao sinal ao importar.

## Tabela 1 (B_th, −Δ, a_bg) — a base do nosso `feshbach.py`
| m_Fa,m_Fb | B_exp (G) | B_th (G) | −Δ (G) | a_bg (a₀) |
|---|---|---|---|---|
| 1,1 | 25,85(10) | 25,9 | 0,47 | −33 |
| 1,1 | 403,4(7) | 402,4 | 52 | −29 |
| 1,1 | — | 745,1 | 0,4 | −35 |
| 1,1 | 752,3(1) | 752,4 | 0,4 | −35 |
| 0,0 | 59,3(6) | 58,8 | 9,6 | −18 |
| 0,0 | 66,0(9) | 65,6 | 7,9 | −18 |
| 0,0 | — | 471 | 72 | −28 |
| 0,0 | — | 490 | 5 | −28 |
| 0,0 | — | 825 | 0,032 | −36 |
| 0,0 | — | 832 | 0,52 | −36 |
| −1,−1 | 32,6(1,5) | 33,6 | −55 | −19 |
| −1,−1 | 162,8(9) | 162,3 | 37 | −19 |
| −1,−1 | 562,2(1,5) | 560,7 | 56 | −29 |

⚠️ O rótulo de estado aqui é (m_Fa, m_Fb) e no Etrych é \|F,m_F⟩ —
correspondência: 1,1 ↔ \|1,1⟩; 0,0 ↔ \|1,0⟩; −1,−1 ↔ \|1,−1⟩.

## O que disso vira código
O `feshbach.py` já usa esta tabela. **Duas correções pendentes:** faltam
as duas ressonâncias de alto campo do estado \|1,0⟩ (825 e 832 G), e as
posições devem migrar para os valores de Etrych 2023 onde existirem.

## Status
Superado em precisão por Etrych 2023 para as ressonâncias medidas lá,
mas continua sendo a fonte mais **completa** (inclui as de alto campo que
o Etrych não cobre).
