# -*- coding: utf-8 -*-
"""
Ressonâncias de Feshbach do 39K — a ponte entre o "botão v" dos nossos
potenciais-modelo e o BOTÃO REAL que se gira num laboratório de átomos
frios: o campo magnético B.

Perto de uma ressonância isolada (Eq. 1 de D'Errico et al., New J. Phys.
9, 223 (2007); arXiv:0705.3036):

    a(B) = a_bg * [ 1 - Delta / (B - B0) ]

  a_bg  : comprimento de espalhamento de fundo (longe da ressonância)
  B0    : centro da ressonância (a diverge: unitariedade!)
  Delta : largura; o zero de a(B) fica em B = B0 + Delta

Dados das ressonâncias em onda-s do 39K (Tabela I do artigo; colunas
B0_th, Delta, a_bg; a coluna do artigo lista -Delta).  Parâmetros do
modelo de colisão: a_S = 138.90 a0, a_T = -33.3 a0, C6 = 3921 u.a.

A ressonância mais usada (e a do grupo do 39K) é a LARGA em 402.4 G do
estado |F=1, mF=1>: Delta = -52 G, a_bg = -29 a0, zero de a em 350.4 G.
É nela que se estuda física de Efimov no 39K.
"""
import numpy as np

# ---------------------------------------------------------- constantes
A0_EM_ANGSTROM = 0.529177210903      # raio de Bohr em Angstrom
A0_EM_FM = 0.529177210903e5          # raio de Bohr em fm (1 A = 1e5 fm)

# ------------------------------------------------- dados do 39K (l = 0)
# (estado, B0 [G], Delta [G], a_bg [a0]); Tabela I de D'Errico 2007
RESSONANCIAS_39K = {
    "1,1": [
        {"B0": 25.9, "Delta": -0.47, "abg": -33.0},
        # A ressonância larga (valores refinados de Zaccanti et al. 2009,
        # Nat. Phys. 5, 586: B0 = 402.50(3) G, Delta = -52.1(1) G,
        # a_bg = -29.0(3) a0; D'Errico 2007 dava B0 = 402.4(2) G)
        {"B0": 402.50, "Delta": -52.1, "abg": -29.0},
        {"B0": 745.1, "Delta": -0.4, "abg": -35.0},
        {"B0": 752.4, "Delta": -0.4, "abg": -35.0},
    ],
    "1,0": [
        {"B0": 58.8, "Delta": -9.6, "abg": -18.0},
        {"B0": 65.6, "Delta": -7.9, "abg": -18.0},
        {"B0": 471.0, "Delta": -72.0, "abg": -28.0},
        {"B0": 490.0, "Delta": -5.0, "abg": -28.0},
    ],
    "1,-1": [
        {"B0": 33.6, "Delta": 55.0, "abg": -19.0},
        {"B0": 162.3, "Delta": -37.0, "abg": -19.0},
        {"B0": 560.7, "Delta": -56.0, "abg": -29.0},
    ],
}

C6_39K_AU = 3921.0        # coeficiente de van der Waals (u.a.)
M_39K_U = 38.9637064864   # massa do 39K em u


def a_de_B(B, B0, Delta, abg):
    """Eq. (1): a(B) em unidades de a0.  Diverge em B = B0."""
    B = np.asarray(B, dtype=float)
    return abg * (1.0 - Delta / (B - B0))


def B_de_a(a_alvo, B0, Delta, abg):
    """Inverte a(B): o campo que produz um dado a (mesmo ramo da ressonância)."""
    return B0 + Delta * abg / (abg - a_alvo)


def zero_de_a(B0, Delta):
    """Campo onde a(B) = 0 (ex.: 402.50 - 52.1 = 350.40 G na larga do |1,1>)."""
    return B0 + Delta


def comprimento_vdw_a0(C6_au=C6_39K_AU, m_u=M_39K_U):
    """l_vdW = (1/2) (2 m_r C6 / hbar^2)^(1/4) em unidades de a0
    (definição de Chin et al., Rev. Mod. Phys. 82, 1225 (2010)).

    Em unidades atômicas (hbar = m_e = a0 = 1): C6 já está em u.a. e
    m_r = (m_u/2) * 1822.888486 m_e.  Escala natural do potencial real:
    para o 39K dá ~64.6 a0."""
    m_r_me = 0.5 * m_u * 1822.888486209
    return 0.5 * (2.0 * m_r_me * C6_au)**0.25


def resumo_ressonancia_larga():
    """A ressonância de trabalho do 39K |1,1> em números."""
    r = RESSONANCIAS_39K["1,1"][1]
    lvdw = comprimento_vdw_a0()
    return {
        "estado": "|F=1, mF=1>",
        "B0_G": r["B0"], "Delta_G": r["Delta"], "abg_a0": r["abg"],
        "zero_de_a_G": zero_de_a(r["B0"], r["Delta"]),
        "l_vdW_a0": lvdw,
        # a1- do Efimov medido no 39K: -1500(90) a0 (Zaccanti et al. 2009;
        # o pico em ~-650 a0 é o TETRÂMERO associado, não o trímero)
        "B_do_primeiro_trimero_G": B_de_a(-1500.0, **r),
    }
