# -*- coding: utf-8 -*-
"""
Extração do comprimento de espalhamento a e do alcance efetivo r0
a partir da solução numérica u(r) (Seção 3 do artigo).

  a  (Eq. 110):  a = R - 2 dr u(R) / [u(R+dr) - u(R-dr)]
                 (derivada logarítmica de u casada com g0(r) = 1 - r/a em R)

  r0 (Eq. 56):   r0 = 2 * int_0^R [ g0^2(r) - u0^2(r) ] dr
                 com u0 = C u normalizada por C = (1 - R/a)/u(R)  (Eq. 111),
                 de modo que C u(R) = g0(R).

A integral é feita com trapézio E Simpson, para comparação.
Para o LJ (r_min > 0) o trecho [0, r_min], onde u = 0, contribui com
int_0^{r_min} g0^2 dr = r_min - r_min^2/a + r_min^3/(3 a^2)  (analítico).
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import trapezoid, simpson

from . import solvers


@dataclass
class Resultado:
    a: float            # comprimento de espalhamento (fm)
    r0_trap: float      # alcance efetivo, trapézio (fm)
    r0_simp: float      # alcance efetivo, Simpson (fm)
    nos: int            # número de nós de u em (0, R)
    r: np.ndarray       # grade radial
    u: np.ndarray       # solução normalizada (C u)
    R_match: float      # ponto de casamento (>= alcance R)
    dr: float
    metodo: str

    @property
    def r0(self):       # valor "oficial": Simpson
        return self.r0_simp


def _trapezio(y, dr):
    """v1: scipy.integrate.trapezoid (a v0 tinha a soma explícita)."""
    return trapezoid(y, dx=dr)


def _simpson(y, dr):
    """v1: scipy.integrate.simpson, que trata sozinho o caso de número
    ímpar de intervalos (a v0 fechava o último intervalo com trapézio)."""
    return simpson(y, dx=dr)


def contar_nos(u, N, V, dr):
    """Número de trocas de sinal de u no intervalo aberto (0, R].

    Pontos do caroço repulsivo numericamente não resolvido
    (2 dr^2 V > 1, região classicamente proibida onde u é
    exponencialmente pequena e não tem nós físicos) são ignorados:
    lá o método de Numerov pode oscilar de sinal espuriamente.
    """
    ok = 2.0 * dr * dr * V[1:N + 1] <= 1.0
    ui = u[1:N + 1][ok]
    sinais = np.sign(ui[np.abs(ui) > 0.0])
    return int(np.count_nonzero(sinais[1:] * sinais[:-1] < 0))


def calcular(pot, dr=1e-3, metodo="numerov"):
    """Resolve o potencial e devolve a, r0 (trapézio e Simpson) e nós."""
    r, u, N, V, dr = solvers.integrar(pot, dr, metodo)   # dr = passo efetivo
    R = r[N]                                        # borda exata do alcance

    # --- comprimento de espalhamento (Eq. 110) ---------------------------
    # na unitariedade exata u(R+dr) = u(R-dr) e a diverge: +-inf é um
    # resultado legítimo aqui (1/a = 0), então silenciamos o aviso.
    with np.errstate(divide="ignore"):
        a = R - 2.0 * dr * u[N] / (u[N + 1] - u[N - 1])

    # --- normalização (Eq. 111): C u(R) = g0(R) = 1 - R/a -----------------
    C = (1.0 - R / a) / u[N]
    un = C * u

    # --- alcance efetivo (Eq. 56) ----------------------------------------
    g0 = 1.0 - r / a
    integrando = g0[:N + 1]**2 - un[:N + 1]**2
    # trecho analítico [0, r_min] (u = 0 dentro do caroço do LJ)
    rm = pot.r_min
    extra = 2.0 * (rm - rm**2 / a + rm**3 / (3.0 * a**2)) if rm > 0 else 0.0

    r0_trap = 2.0 * _trapezio(integrando, dr) + extra
    r0_simp = 2.0 * _simpson(integrando, dr) + extra

    # --- número de nós -----------------------------------------------------
    # nós dentro de (0, R] + o nó EXTERIOR da solução g0 = 1 - r/a em r = a,
    # que existe quando 0 < R < a (ex.: poço do dêuteron, a = 5.4 > R = 2 fm).
    # Na unitariedade (a -> +infinito) o nó está "no infinito" e não é
    # contado; usamos o corte pragmático a < 1e4 fm.
    nos = contar_nos(u, N, V, dr)
    if R < a < 1e4:
        nos += 1

    return Resultado(a=a, r0_trap=r0_trap, r0_simp=r0_simp,
                     nos=nos, r=r, u=un,
                     R_match=R, dr=dr, metodo=metodo)
