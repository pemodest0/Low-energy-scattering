# -*- coding: utf-8 -*-
"""
Testes de âncora do src/trimero.py.

Regra da casa: todo teste ancora em resultado ANALÍTICO ou PUBLICADO,
nunca em "o que o código deu ontem".

Os três últimos são testes de REGRESSÃO de bugs reais encontrados em
2026-08. Cada um passou despercebido por um tempo e produziu número errado
sem levantar exceção — o pior tipo.
"""
import math
import numpy as np
import pytest

from src import trimero as T


# =============================================================================
#  ÂNCORAS ANALÍTICAS
# =============================================================================
def test_s0_bate_com_braaten_hammer():
    """s0 = 1.0062378 — Braaten & Hammer, Phys. Rep. 428, 259 (2006)."""
    assert T.S0 == pytest.approx(1.0062378, abs=1e-7)


def test_razoes_universais():
    """e^(pi/s0) = 22.694 e e^(2pi/s0) = 515.03."""
    assert T.RAZAO_ESCALA == pytest.approx(22.694, abs=1e-3)
    assert T.RAZAO_ENERGIA == pytest.approx(515.03, abs=0.05)


def test_transcendental_e_de_fato_satisfeita():
    """A raiz resolve  s cosh(pi s/2) = (8/sqrt3) sinh(pi s/6)."""
    s = T.S0
    esq = s * math.cosh(math.pi * s / 2)
    dir_ = (8 / math.sqrt(3)) * math.sinh(math.pi * s / 6)
    assert esq == pytest.approx(dir_, rel=1e-12)


# =============================================================================
#  A TORRE
# =============================================================================
@pytest.fixture(scope="module")
def torre_4():
    return T.torre(4, npts=60000)


def test_torre_e_geometrica_nos_excitados(torre_4):
    """kappa_n/kappa_(n+1) -> e^(pi/s0) para os estados excitados.

    O FUNDAMENTAL fica de fora de propósito: ele é o mais próximo da parede
    e portanto o menos universal (desvio medido: +0.05%).
    """
    ks = torre_4
    razoes = [ks[i - 1] / ks[i] for i in range(2, len(ks))]
    for r in razoes:
        assert r == pytest.approx(T.RAZAO_ESCALA, rel=1e-5)


def test_razao_de_energias(torre_4):
    """E_n/E_(n+1) = (kappa_n/kappa_(n+1))^2 = 515.03."""
    ks = torre_4
    r = (ks[-2] / ks[-1]) ** 2
    assert r == pytest.approx(T.RAZAO_ENERGIA, rel=1e-4)


def test_fundamental_e_o_menos_universal(torre_4):
    """O desvio se CONCENTRA no fundamental — se estivesse espalhado,
    seria sinal de erro numérico e não de física."""
    ks = torre_4
    desvio_fund = abs(ks[0] / ks[1] / T.RAZAO_ESCALA - 1)
    desvio_exc = abs(ks[1] / ks[2] / T.RAZAO_ESCALA - 1)
    assert desvio_fund > 10 * desvio_exc
    assert desvio_fund < 0.01          # mas ainda pequeno: 0.05%


def test_espectro_escala_com_a_parede():
    """rho_0 é a ÚNICA escala: dobrar a parede divide todo kappa por 2.

    Consequência direta da invariância de escala do problema.
    """
    k_a = T.kappa_n(1, x0=0.0, npts=40000)
    k_b = T.kappa_n(1, x0=math.log(2.0), npts=40000)
    assert k_a / k_b == pytest.approx(2.0, rel=2e-3)


# =============================================================================
#  A ARMADILHA
# =============================================================================
@pytest.mark.parametrize("l_ho,n_esperado", [(1e4, 3), (1e3, 2), (1e2, 1), (10.0, 0)])
def test_armadilha_trunca_a_torre(l_ho, n_esperado):
    """N = (s0/pi) ln(l_ho/rho_0) estados de Efimov sobrevivem.

    A torre infinita fica finita: o degrau mais alto é o maior, e é o
    primeiro a não caber na armadilha.
    """
    omega = 1.0 / l_ho ** 2
    n = 0
    for nivel in range(6):
        E = T.energia_n(nivel, omega, npts=25000, iteracoes=55)
        if E < 0:
            n += 1
        else:
            break
    assert n == n_esperado


def test_previsao_de_contagem_bate_com_o_numerico():
    """A fórmula N ~ (s0/pi) ln(l_ho/rho0) arredonda para o valor achado."""
    for l_ho in (1e4, 1e3, 1e2):
        omega = 1.0 / l_ho ** 2
        prev = T.n_efimov_previsto(omega)
        n = 0
        for nivel in range(6):
            if T.energia_n(nivel, omega, npts=25000, iteracoes=55) < 0:
                n += 1
            else:
                break
        assert abs(prev - n) < 1.0


# =============================================================================
#  REGRESSÃO — bugs reais de 2026-08
# =============================================================================
def test_regressao_underflow_apaga_os_nos():
    """BUG 1 (o pior): o reescalonamento contra overflow divide o array todo
    por 1e250. Com armadilha forte isso dispara várias vezes e a região
    oscilatória do INÍCIO faz underflow para zero exato — os nós somem.

    Sintoma: a torre desaparecia com armadilha, sem erro nenhum.
    Conserto: contar os nós DURANTE a integração.
    """
    x, w, nos = T.integrar(0.0, omega=1e-8, npts=40000)
    assert nos == 3                                  # contados na hora: certo
    assert T.contar_nos(w) == 0                      # contados no fim: errado
    # ^ a segunda linha documenta o bug; contar_nos está marcada como depreciada


def test_regressao_overflow_no_produto_de_sinais():
    """BUG 2: detectar nó com w[i+1]*w[i] < 0 estoura quando ambos ~1e250.
    Conserto: comparar sinais em vez de multiplicar.
    """
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error")               # overflow vira exceção
        T.integrar(-0.5 * 0.05 ** 2, 0.0, npts=20000)


def test_regressao_colisao_de_nome_n():
    """BUG 3: 'n' era o índice do estado E o tamanho da grade.
    A grade agora se chama npts.
    """
    T.kappa_n(1, npts=20000)                          # não deve levantar
    with pytest.raises(TypeError):
        T.integrar(-1e-3, n=20000)                    # o nome antigo sumiu


# =============================================================================
#  CONSISTÊNCIA INTERNA
# =============================================================================
def test_sem_armadilha_e_o_limite_de_armadilha_fraca():
    """omega -> 0 tem que reproduzir o resultado sem armadilha."""
    k1 = T.kappa_n(1, npts=40000)
    E_sem = -0.5 * k1 ** 2
    E_com = T.energia_n(1, omega=1e-14, npts=40000, iteracoes=70)
    assert E_com == pytest.approx(E_sem, rel=0.02)
