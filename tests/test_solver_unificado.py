# -*- coding: utf-8 -*-
"""
Prova de que a equação radial é UMA só: u'' = [2V + (L²-¼)/r²] u, com

    L = l + (d-2)/2

Cobre 1D, 2D, 3D, ondas parciais, e verifica a consistência com o d = 6
do hiper-raio de Efimov.

Rode:  python -m pytest tests/test_solver_unificado.py -v
"""
import numpy as np
import pytest

from src import solvers
from src.solvers import L_efetivo, termo_centrifugo
from src.potenciais import PocoEsferico


# ============================================================ identidades
@pytest.mark.parametrize("l", [0, 1, 2, 3, 4, 5])
def test_3D_reproduz_l_vezes_l_mais_1(l):
    """d=3: (L²-¼) tem que ser EXATAMENTE l(l+1) — o termo centrífugo clássico."""
    L = L_efetivo(l, d=3)
    assert L == pytest.approx(l + 0.5, abs=1e-15)
    assert L**2 - 0.25 == pytest.approx(l * (l + 1), abs=1e-12)


@pytest.mark.parametrize("m", [0, 1, 2, 3])
def test_2D_da_m2_menos_um_quarto(m):
    """d=2: (L²-¼) = m² - ¼ — o termo centrífugo bidimensional."""
    L = L_efetivo(m, d=2)
    assert L == pytest.approx(m, abs=1e-15)
    assert L**2 - 0.25 == pytest.approx(m**2 - 0.25, abs=1e-12)


def test_1D_onda_s_sem_termo_centrifugo():
    """d=1, l=0: L = -½ e o termo se anula — partícula livre 1D."""
    assert L_efetivo(0, d=1)**2 - 0.25 == pytest.approx(0.0, abs=1e-15)


def test_d6_e_a_equacao_de_efimov():
    """d=6: L = s, logo (L²-¼) = s²-¼.

    Com s = i·s0 (imaginário), s² = -s0² e o termo vira -(s0²+¼)/R² —
    exatamente o potencial hiperradial usado em src/efimov.py.  O ¼ é
    geometria (jacobiano em 6D), não física de Efimov.
    """
    from src import efimov
    s0 = efimov.s0_universal()

    # A cadeia geométrica: em d=6 o jacobiano gera (d-1)(d-3)/4 = 15/4,
    # e o autovalor hiperangular de Efimov é lambda = s^2 - 4.
    assert (6 - 1) * (6 - 3) / 4 == pytest.approx(3.75, abs=1e-15)
    lam = -s0**2 - 4.0                       # com s = i*s0
    coef = lam + 3.75                        # (L^2 - 1/4)
    assert coef == pytest.approx(-(s0**2 + 0.25), abs=1e-12), (
        "d=6 tem que reproduzir o -(s0^2 + 1/4)/R^2 do efimov.py")

    # e L_efetivo NÃO se aplica aqui: o autovalor não é l(l+d-2)
    assert L_efetivo(s0, d=6) != pytest.approx(s0, abs=1e-6), (
        "L_efetivo(s0, 6) daria s0+2 — o hiper-raio usa L=s direto")


def test_compatibilidade_onda_s_3D():
    """O padrão (l=0, d=3) tem termo centrífugo EXATAMENTE nulo."""
    r = np.linspace(0.0, 5.0, 50)
    assert np.all(termo_centrifugo(r, L_efetivo(0, 3)) == 0.0)


def test_origem_nao_gera_nan():
    r = np.linspace(0.0, 5.0, 50)
    t = termo_centrifugo(r, L_efetivo(2, 3))
    assert np.isfinite(t).all()
    assert t[0] == 0.0


# ====================================== validação contra solução exata
class _Livre:
    """V = 0 em toda parte. Solução exata: u = r^{L+1/2}."""
    nome, rotulo = "livre", "Partícula livre"
    r_min, R = 0.0, 10.0

    def V(self, r):
        return np.zeros_like(np.asarray(r, dtype=float))


# L = 0 (o caso d=2, m=0) fica de fora de propósito: ali as duas soluções
# são DEGENERADAS (sqrt(r) e sqrt(r)·ln r) e não há potência pura.
@pytest.mark.parametrize("l,d", [(0, 3), (1, 3), (2, 3), (3, 3),
                                 (1, 2), (2, 2),
                                 (0, 1)])
def test_particula_livre_da_potencia_exata(l, d):
    """Com V=0, u'' = (L²-¼)/r² u tem solução exata u ∝ r^{L+½}.

    Este é o teste que prova a unificação: uma única linha de código
    reproduz a potência certa em três dimensionalidades diferentes.
    """
    L = L_efetivo(l, d)
    r, u, N, _, _ = solvers.integrar(_Livre(), dr=2e-4, l=l, d=d)

    # compara a razão u(r1)/u(r2) com (r1/r2)^{L+1/2}, longe da origem
    i, j = int(0.7 * N), int(0.35 * N)
    obtido = u[i] / u[j]
    exato = (r[i] / r[j]) ** (abs(L) + 0.5)     # solução regular
    assert obtido == pytest.approx(exato, rel=2e-4), (
        f"l={l}, d={d}, L={L}: esperava expoente {abs(L) + 0.5}")


# ================================================== física: barreira
def test_barreira_centrifuga_desliga_o_estado_ligado():
    """O poço do dêuteron liga em onda-s (1 nó). Com l=1 a barreira
    centrífuga tira a ligação (0 nós). É o efeito físico do termo novo."""
    # poço fundo o bastante para ter nó DENTRO do alcance:
    # K0*R = sqrt(2v) = 4 > pi (onda-s liga com nó interior), mas
    # 4 < 4.4934 (1o zero de j1) => onda-p ainda NÃO liga.
    poco = PocoEsferico(8.0, 1.0)

    def nos(l):
        _, u, N, _, _ = solvers.integrar(poco, dr=1e-3, l=l, d=3)
        s = np.sign(u[2:N + 1])
        return int(np.sum(np.diff(s) != 0))

    assert nos(0) == 1, "onda-s deve ter 1 nó interior"
    assert nos(1) == 0, "a barreira centrífuga deve desligar o estado"


def test_barreira_cresce_com_l():
    """Quanto maior l, mais a função é empurrada para fora: u(r) na região
    do poço cai monotonicamente com l (normalizando no mesmo ponto)."""
    poco = PocoEsferico(8.0, 1.0)
    pesos = []
    for l in (0, 1, 2, 3):
        r, u, N, _, _ = solvers.integrar(poco, dr=1e-3, l=l, d=3)
        u = u / np.max(np.abs(u[: N + 1]))
        dentro = r[: N + 1] < poco.R * 0.5
        pesos.append(np.trapezoid(np.abs(u[: N + 1][dentro]),
                                  r[: N + 1][dentro]))
    assert all(pesos[i] > pesos[i + 1] for i in range(len(pesos) - 1)), pesos
