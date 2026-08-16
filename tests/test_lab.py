# -*- coding: utf-8 -*-
"""
Suite de testes do laboratório — rede de segurança para a reescrita v1.

Rode na raiz do projeto:   python -m pytest tests/ -v

Cada teste ancora o código num resultado ANALÍTICO ou num valor publicado
no artigo (RBEF 45, e20230079).  As tolerâncias vêm dos erros medidos no
estudo de convergência (resultados/convergencia.csv); se um teste quebrar
durante a reescrita, o módulo recém-alterado é o suspeito.
"""
import math

import pytest

from src import espalhamento, analitico
from src.potenciais import PocoEsferico, PoschlTeller, Gaussiano, LennardJones
from src.bench_referencias import EsferaMole, a_esfera_mole, a_fase_variavel
from src import constantes as cte


# ---------------------------------------------------------------- digamma
def test_digamma_valores_conhecidos():
    # psi(1) = -gamma ; psi(1/2) = -gamma - 2 ln 2 ; psi(2) = 1 - gamma
    g = cte.EULER_GAMMA
    assert analitico.digamma(1.0) == pytest.approx(-g, abs=1e-12)
    assert analitico.digamma(0.5) == pytest.approx(-g - 2 * math.log(2), abs=1e-12)
    assert analitico.digamma(2.0) == pytest.approx(1.0 - g, abs=1e-12)


# ------------------------------------------------- poço esférico (Eq. 80/92)
def test_poco_a_e_r0_analiticos():
    v, mu = 1.1096, 0.3918                      # caso n-n da Tabela 3
    R = 1.0 / mu
    res = espalhamento.calcular(PocoEsferico(v, mu), dr=1e-3, metodo="central")
    assert res.a == pytest.approx(analitico.a_poco(v, R), rel=5e-6)
    assert res.r0_simp == pytest.approx(analitico.r0_poco(v, R), rel=5e-7)
    assert res.nos == 0


def test_poco_limiar_estado_ligado():
    # logo abaixo do limiar v = pi^2/8: a << 0; logo acima: a >> 0 (1 nó)
    eps = 5e-3
    abaixo = espalhamento.calcular(PocoEsferico(analitico.V_LIMIAR_POCO - eps, 1.0), dr=5e-4)
    acima = espalhamento.calcular(PocoEsferico(analitico.V_LIMIAR_POCO + eps, 1.0), dr=5e-4)
    assert abaixo.a < -50 and abaixo.nos == 0
    assert acima.a > 50 and acima.nos == 1


# ------------------------------------------ Pöschl-Teller mod. (Eq. 117/118)
def test_mpt_a_analitico_numerov():
    v, mu = 0.9071, 0.7991                      # caso n-n da Tabela 3
    res = espalhamento.calcular(PoschlTeller(v, mu), dr=1e-3, metodo="numerov")
    assert res.a == pytest.approx(analitico.a_mpt(v, mu), rel=1e-7)


def test_mpt_unitario_r0():
    res = espalhamento.calcular(PoschlTeller(1.0, 2.0), dr=1e-3)
    assert abs(res.a) > 1e9                     # unitariedade
    assert res.r0_simp == pytest.approx(analitico.r0_mpt_unitario(2.0), rel=1e-7)
    assert res.nos == 0


# ------------------------------------------------ gaussiano (Tabela 3)
def test_gauss_deuteron_tabela3():
    res = espalhamento.calcular(Gaussiano(1.9102, 0.6754), dr=1e-3)
    assert res.a == pytest.approx(5.4, rel=2e-3)
    assert res.r0_simp == pytest.approx(1.70, rel=2e-3)
    assert res.nos == 1


# -------------------------------------------- Lennard-Jones (Tabela 4)
@pytest.mark.parametrize("C12,C6,a_alvo,r0_alvo,nos", [
    (3.08836698, 9.86668911, -18.5, 2.71, 0),      # n-n
    (0.90485319, 6.81472000, 5.4, 1.70, 1),        # dêuteron
])
def test_lj_tabela4(C12, C6, a_alvo, r0_alvo, nos):
    res = espalhamento.calcular(LennardJones(C12, C6), dr=2e-3)
    assert res.a == pytest.approx(a_alvo, rel=5e-3)
    assert res.r0_simp == pytest.approx(r0_alvo, rel=5e-3)
    assert res.nos == nos


def test_lj_unitario_tabela4():
    res = espalhamento.calcular(LennardJones(0.00034068, 0.26462461), dr=1e-3)
    assert abs(res.a) > 1e3                     # |a| enorme na unitariedade
    assert res.r0_simp == pytest.approx(1.00, rel=2e-3)
    assert res.nos == 0


# ------------------------------------- esfera mole (Pera & Boronat [17])
def test_esfera_mole_analitico_e_limite_duro():
    res = espalhamento.calcular(EsferaMole(10.0, 1.0), dr=5e-4)
    assert res.a == pytest.approx(a_esfera_mole(10.0, 1.0), rel=1e-4)
    dura = espalhamento.calcular(EsferaMole(1e8, 1.0), dr=5e-4)
    assert dura.a == pytest.approx(1.0, rel=1e-3)          # a -> R
    assert dura.r0_simp == pytest.approx(2.0 / 3.0, rel=1e-3)  # r0 -> 2R/3


# ------------------------------ fase variável (Calogero / Viterbo [19])
def test_fase_variavel_cruzada():
    pot = Gaussiano(1.2121, 0.5672)             # n-n, sem polo
    a_fv = a_fase_variavel(pot, dr=2e-4)
    a_num = espalhamento.calcular(pot, dr=1e-3).a
    assert a_fv == pytest.approx(a_num, rel=1e-6)


def test_fase_variavel_atravessa_polo():
    pot = Gaussiano(1.9102, 0.6754)             # dêuteron: a(r) tem polo
    a_fv = a_fase_variavel(pot, dr=2e-4)
    assert a_fv == pytest.approx(5.4, rel=2e-3)


# ------------------------------------------- energias (Tabela 1 / Eq. 93/95)
def test_energias_deuteron():
    h22m = cte.hbar2_over_2mr_MeV_fm2(cte.M_R_NP_MEV)
    assert analitico.energia_zr(5.4112, h22m) == pytest.approx(-1.416, rel=1e-3)
    assert analitico.energia_fr(5.4112, 1.7436, h22m) == pytest.approx(-2.223, rel=1e-3)


def test_energias_dimero_he():
    h22m = cte.hbar2_over_2mr_K_A2(cte.M_R_HE_U)
    assert analitico.energia_zr(90.4, h22m) * 1e3 == pytest.approx(-1.48, rel=5e-3)
    assert analitico.energia_fr(90.4, 8.0, h22m) * 1e3 == pytest.approx(-1.63, rel=5e-3)


# ------------------------------------------------ propriedades estruturais
def test_normalizacao_u_iguala_g0_em_R():
    """C u(R) deve valer g0(R) = 1 - R/a (Eq. 111)."""
    res = espalhamento.calcular(PoschlTeller(0.9071, 0.7991), dr=1e-3)
    N = round((res.R_match - res.r[0]) / res.dr)
    g0_R = 1.0 - res.R_match / res.a
    assert res.u[N] == pytest.approx(g0_R, rel=1e-12)


def test_metodos_concordam():
    for pot in (PoschlTeller(1.4388, 0.8631), Gaussiano(1.2121, 0.5672)):
        a1 = espalhamento.calcular(pot, dr=1e-3, metodo="central").a
        a2 = espalhamento.calcular(pot, dr=1e-3, metodo="numerov").a
        assert a1 == pytest.approx(a2, rel=1e-5)


# --------------------------- gaussiano vs Jeszenszki et al. [34], Eq. (25)
def test_gauss_vs_jeszenszki():
    from src.bench_referencias import a_gauss_jeszenszki, W_JESZ
    for v in (0.5, 1.0, 1.5, 3.0):
        a_num = espalhamento.calcular(Gaussiano(v, 1.0), dr=1e-3).a
        assert a_num == pytest.approx(a_gauss_jeszenszki(v), rel=2e-4)
    # limiar do 1º estado ligado: W1/2 = 1.3420023 (nosso ajuste: 1.3420)
    assert W_JESZ[0] / 2 == pytest.approx(1.3420, rel=1e-4)


# ------------------------------------ Feshbach 39K (D'Errico et al. 2007)
def test_feshbach_39K():
    from src import feshbach as fb
    r = fb.RESSONANCIAS_39K["1,1"][1]           # a ressonância larga
    # zero de a em B0 + Delta = 350.4 G (valor citado no artigo)
    assert fb.zero_de_a(r["B0"], r["Delta"]) == pytest.approx(350.4, abs=0.01)
    assert fb.a_de_B(350.4, **r) == pytest.approx(0.0, abs=1e-10)
    # longe da ressonância, a -> a_bg
    assert fb.a_de_B(200.0, **r) == pytest.approx(r["abg"], rel=0.3)
    # entre o zero e B0 o a é positivo e cresce ao se aproximar de B0
    assert fb.a_de_B(400.0, **r) > fb.a_de_B(380.0, **r) > 0
    # inversão B(a) consistente
    assert fb.a_de_B(fb.B_de_a(-1500.0, **r), **r) == pytest.approx(-1500.0)
    # l_vdW do potássio ~ 64.6 a0 (Chin et al. 2010)
    assert fb.comprimento_vdw_a0() == pytest.approx(64.6, rel=0.01)


# ------------------------------------------ Aziz HFD-B liga o dímero [25]
def test_aziz_liga_dimero():
    from src.bench_referencias import AzizHFDB
    res = espalhamento.calcular(AzizHFDB(), dr=4e-3)
    assert res.a == pytest.approx(88.4, rel=0.01)   # literatura: ~88.5 A
    assert res.nos == 1                              # o dímero EXISTE
    # profundidade do mínimo = eps = 10.948 K
    import numpy as np
    pot = AzizHFDB()
    rr = np.linspace(2.0, 5.0, 4000)
    assert float(pot._V_K(rr).min()) == pytest.approx(-10.948, rel=1e-3)


# ---------------------------------- singleto n-p: estado virtual, não ligado
def test_singleto_np_virtual():
    a_s, r_s = -23.74, 2.77
    kappa = (1 - math.sqrt(1 - 2 * r_s / a_s)) / r_s
    assert kappa < 0                                 # virtual (não ligado)
    h22m = cte.hbar2_over_2mr_MeV_fm2(cte.M_R_NP_MEV)
    assert abs(h22m * kappa**2) * 1e3 == pytest.approx(66.0, rel=0.01)  # keV


# --------------------------------- torre de Efimov (equação hiperradial)
def test_efimov_s0_e_torre():
    from src import efimov
    # s0 resolvido da equação transcendental, sem valor tabelado
    assert efimov.S0 == pytest.approx(1.00624, abs=2e-5)
    ks = efimov.torre(4)
    assert len(ks) == 4
    rz = efimov.razoes(ks)
    # níveis profundos sentem o corte; os rasos convergem ao universal
    esperado = math.exp(math.pi / efimov.S0)
    assert rz[-1] == pytest.approx(esperado, rel=1e-4)
    assert rz[-1]**2 == pytest.approx(515.03, rel=2e-4)   # razão de energias


# ------------------------- Schrödinger canônica (src/schrodinger.py)
def test_oscilador_harmonico():
    from src import schrodinger as sq
    for E, n in sq.autovalores(sq.V_oscilador, -8, 8, 0.0, 5.0,
                               n_estados=4, dx=4e-3):
        assert E == pytest.approx(sq.E_oscilador(n), abs=2e-6)


def test_hidrogenio_radial():
    from src import schrodinger as sq
    niveis = sq.autovalores(sq.V_hidrogenio(0), 1e-6, 60, -0.6, -0.03,
                            n_estados=2, dx=2e-3)
    for E, n in niveis:                      # n nós -> estado n+1: E=-1/2(n+1)^2
        assert E == pytest.approx(sq.E_hidrogenio(n + 1), rel=2e-4)


def test_morse_anarmonico():
    from src import schrodinger as sq
    assert sq.n_max_morse() == 3             # 4 níveis para D=10, a=1
    for E, n in sq.autovalores(sq.V_morse(), 0.05, 12, -9.99, -0.02,
                               n_estados=4, dx=4e-3):
        assert E == pytest.approx(sq.E_morse(n), rel=1e-4)


def test_tunelamento_barreira():
    from src import schrodinger as sq
    for E in (1.0, 3.0, 6.0, 8.0):
        Tn = sq.transmissao(sq.V_barreira(), E, -0.5, 1.5, dx=1e-4)
        assert Tn == pytest.approx(sq.T_barreira_analitico(E), abs=1e-4)
    # sanidade física: T cresce com E e satura em 1
    assert sq.transmissao(sq.V_barreira(), 30.0, -0.5, 1.5, dx=1e-4) > 0.99
