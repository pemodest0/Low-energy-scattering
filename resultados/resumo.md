# Resumo dos resultados — laboratório de espalhamento

Reprodução de Macêdo-Lima & Madeira, RBEF 45, e20230079 (2023).


## Validação analítica (poço esférico e mPT)

| potencial | caso | método | a_num (fm) | a_analítico | erro rel a | r0_Simpson | r0_analítico | erro rel r0 | nós |
|---|---|---|---|---|---|---|---|---|---|
| poco | nn | central | -18.528 | -18.528 | 3.46e-07 | 2.694613 | 2.694613 | 1.50e-08 | 0 |
| poco | nn | numerov | -18.521 | -18.528 | 3.40e-04 | 2.694615 | 2.694613 | 7.67e-07 | 0 |
| poco | unitario | central | -3.3731e+06 | -1.8177e+06 | 2.54e-07 | 1.000000 | 1.000000 | 1.03e-07 | 0 |
| poco | unitario | numerov | -9674.1 | -1.8177e+06 | 1.03e-04 | 1.000000 | 1.000000 | 3.49e-09 | 0 |
| poco | deuteron | central | 5.3999 | 5.3999 | 8.22e-08 | 1.697804 | 1.697804 | 3.49e-08 | 1 |
| poco | deuteron | numerov | 5.4003 | 5.3999 | 7.85e-05 | 1.697799 | 1.697804 | 3.01e-06 | 1 |
| mpt | nn | central | -18.515 | -18.515 | 3.94e-07 | 2.701731 | — | — | 0 |
| mpt | nn | numerov | -18.515 | -18.515 | 6.60e-09 | 2.701731 | — | — | 0 |
| mpt | unitario | central | 2.8127e+06 | -6.4133e+15 | 3.56e-07 | 0.999999 | 1.000000 | 7.43e-07 | 0 |
| mpt | unitario | numerov | -4.2087e+10 | -6.4133e+15 | 2.38e-11 | 1.000000 | 1.000000 | 2.22e-09 | 0 |
| mpt | deuteron | central | 5.4002 | 5.4002 | 2.07e-07 | 1.729951 | — | — | 1 |
| mpt | deuteron | numerov | 5.4002 | 5.4002 | 4.65e-13 | 1.729951 | — | — | 1 |

## Reprodução das Tabelas 3 e 4 (parâmetros publicados)

| tabela | caso | potencial | a_num | a_alvo | erro rel a | r0_num | r0_alvo | erro rel r0 | nós (esperado) |
|---|---|---|---|---|---|---|---|---|---|
| Tabela 3 | nn | poco | -18.52 | -18.5 | 1.16e-03 | 2.6946 | 2.7 | 1.99e-03 | 0 (0) |
| Tabela 3 | nn | mpt | -18.52 | -18.5 | 8.18e-04 | 2.7017 | 2.7 | 6.41e-04 | 0 (0) |
| Tabela 3 | nn | gauss | -18.55 | -18.5 | 2.92e-03 | 2.7031 | 2.7 | 1.14e-03 | 0 (0) |
| Tabela 3 | unitario | poco | -9674 | ±inf | 1.03e-04 | 1.0000 | 1.0 | 2.19e-07 | 0 (0) |
| Tabela 3 | unitario | mpt | -4.209e+10 | ±inf | 2.38e-11 | 1.0000 | 1.0 | 2.22e-09 | 0 (0) |
| Tabela 3 | unitario | gauss | -4.506e+05 | ±inf | 2.22e-06 | 1.0002 | 1.0 | 2.27e-04 | 0 (0) |
| Tabela 3 | deuteron | poco | 5.4 | 5.4 | 5.49e-05 | 1.6978 | 1.7 | 1.29e-03 | 1 (1) |
| Tabela 3 | deuteron | mpt | 5.4 | 5.4 | 3.19e-05 | 1.7300 | 1.7 | 1.76e-02 | 1 (1) |
| Tabela 3 | deuteron | gauss | 5.4 | 5.4 | 5.00e-05 | 1.6989 | 1.7 | 6.29e-04 | 1 (1) |
| Tabela 4 | nn | lj | -18.5 | -18.5 | 6.32e-05 | 2.7078 | 2.7 | 2.87e-03 | 0 (0) |
| Tabela 4 | unitario | lj | 4.71e+04 | ±inf | 2.12e-05 | 1.0000 | 1.0 | 1.72e-05 | 0 (0) |
| Tabela 4 | deuteron | lj | 5.405 | 5.4 | 9.78e-04 | 1.6980 | 1.7 | 1.15e-03 | 1 (1) |

## Energias de estado ligado (Tabela 1)

| sistema | E_zr | E_zr artigo | E_fr | E_fr artigo | E exp. |
|---|---|---|---|---|---|
| dimero_4He | -0.001483 | -0.00148 | -0.001631 | -0.00163 | -0.00162 |
| deuteron | -1.416 | -1.416 | -2.223 | -2.223 | -2.224 |

## Ajuste de parâmetros (Seção 4.5) vs artigo

| caso | potencial | p1 ajustado | p1 artigo | erro rel | p2 ajustado | p2 artigo | erro rel | a obtido | r0 obtido | nós | convergiu |
|---|---|---|---|---|---|---|---|---|---|---|---|

> Nota de convenção: os C12/C6 da Tabela 4 do artigo reproduzem os alvos somente com u'' = (C12/r^12 - C6/r^6) u, isto é, V_LJ = (hbar^2/2m_r)[...] e não (hbar^2/m_r)[...] como impresso na Eq. (121). Ver src/potenciais.py.
