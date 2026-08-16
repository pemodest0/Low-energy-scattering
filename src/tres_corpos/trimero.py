# -*- coding: utf-8 -*-
"""
Trímero de Efimov no hiper-raio — com e sem armadilha.

===============================================================================
 A EQUAÇÃO
===============================================================================

Três bósons idênticos, alcance zero, na unitariedade. Depois de separar o
centro de massa e resolver o problema hiperangular, sobra UMA equação radial no
hiper-raio rho (ver notas_teoria/DoZeroAoTrimero_ParteII, Sua vez 6 e 9):

    u''(rho) = [ -2E + w^2 rho^2 - (s0^2 + 1/4)/rho^2 ] u(rho)

com s0 = 1.0062378 (raiz de  s cosh(pi s/2) = (8/sqrt3) sinh(pi s/6) ).

O termo de armadilha w^2 rho^2 é EXATO, não aproximado: para três partículas
numa armadilha harmônica ISOTRÓPICA,

    sum_i (1/2) m w^2 r_i^2 = (1/2) m w^2 (3 R_cm^2) + (1/2) m w^2 rho^2

porque rho^2 = xi_1^2 + xi_2^2 e sum_i r_i^2 = 3 R_cm^2 + xi_1^2 + xi_2^2.
A armadilha isotrópica preserva a simetria hiperesférica.

(A armadilha ANISOTRÓPICA — a do experimento — NÃO preserva, e é aí que o
problema fica difícil. Este módulo é a âncora antes daquele passo.)

===============================================================================
 A TRANSFORMAÇÃO QUE FAZ TUDO FUNCIONAR
===============================================================================

Com rho = e^x  e  u = e^(x/2) w(x):

    w''(x) = [ w^2 e^(4x) - 2E e^(2x) - s0^2 ] w(x)

Vantagens, e são grandes:

  * a log-periodicidade fica MANIFESTA — sem armadilha e para kappa rho << 1
    sobra w'' = -s0^2 w, ou seja w = sin(s0 x + phi);
  * a torre de Efimov vive em décadas de rho; em x ela vira uma malha uniforme;
  * a grade fica bem-condicionada onde o problema em rho seria impossível.

===============================================================================
 CONDIÇÕES DE CONTORNO
===============================================================================

  rho -> 0 : parede em rho_0 (o PARÂMETRO DE TRÊS CORPOS). Sem ele o problema
             não está definido — é a anomalia de escala. Aqui: w(x_0) = 0.
  rho -> oo: w -> 0.

===============================================================================
 ÂNCORAS (tests/)
===============================================================================

  * sem armadilha: kappa_n/kappa_(n+1) = e^(pi/s0) = 22.6944
  * sem armadilha: E_n/E_(n+1) = e^(2pi/s0) = 515.03
  * o desvio se concentra no estado FUNDAMENTAL (é o menos universal)
"""

from __future__ import annotations
import math
import numpy as np
from scipy.optimize import brentq

# ----------------------------------------------------------------- constante
def _resolve_s0():
    """Raiz da equação transcendental de Efimov, 3 bósons idênticos, 3D."""
    f = lambda s: s * math.cosh(math.pi * s / 2) - (8 / math.sqrt(3)) * math.sinh(math.pi * s / 6)
    return brentq(f, 0.5, 1.5, xtol=1e-15)

S0 = _resolve_s0()                      # 1.0062378...
RAZAO_ESCALA = math.exp(math.pi / S0)   # 22.6944
RAZAO_ENERGIA = math.exp(2 * math.pi / S0)  # 515.03


# ------------------------------------------------------------------ integrar
def integrar(E, omega=0.0, x0=0.0, xmax=None, npts=120000, s0=S0):
    """Integra w'' = [w^2 e^4x - 2E e^2x - s0^2] w a partir da parede em x0.

    E     : energia (unidades hbar = m = 1; para estado ligado use E < 0)
    omega : frequência da armadilha isotrópica (0 = sem armadilha)
    x0    : ln(rho_0), a parede — é o parâmetro de três corpos

    Devolve (x, w).
    """
    if xmax is None:
        if omega > 0:
            # vai até bem além do comprimento do oscilador
            xmax = x0 + 4.0 - 0.5 * math.log(omega)
        else:
            kap = math.sqrt(max(-2.0 * E, 1e-300))
            xmax = x0 + 6.0 - math.log(kap)
    x = np.linspace(x0, xmax, npts)
    h = x[1] - x[0]
    W = (omega ** 2) * np.exp(4 * x) - 2.0 * E * np.exp(2 * x) - s0 * s0
    w = np.zeros(npts)
    w[1] = 1e-8
    nos = 0
    for i in range(1, npts - 1):
        w[i + 1] = 2 * w[i] - w[i - 1] + h * h * W[i] * w[i]
        # CONTA OS NOS AQUI, e nao no fim: o reescalonamento abaixo faz a parte
        # inicial da solucao fazer underflow para zero exato, apagando os nos.
        # (bug real, achado em 2026-08: a torre sumia com armadilha forte)
        if (w[i + 1] < 0.0) != (w[i] < 0.0) and w[i] != 0.0:
            nos += 1
        if abs(w[i + 1]) > 1e250:            # a equação é linear: reescalona
            w[: i + 2] /= 1e250
    return x, w, nos


def contar_nos(w):
    """DEPRECIADO: use o terceiro valor devolvido por integrar().

    Contar no fim nao e confiavel — o reescalonamento contra overflow apaga
    (por underflow) a parte inicial da solucao.
    """
    s = np.sign(w[np.abs(w) > 0])
    return int(np.count_nonzero(s[1:] * s[:-1] < 0))


# ------------------------------------------------------ espectro sem armadilha
def kappa_n(nivel, x0=0.0, lo=1e-9, hi=5.0, iteracoes=70, **kw):
    """n-ésimo estado de Efimov SEM armadilha. Devolve kappa_n (E = -kappa^2/2).

    A bisseção é feita em LOG, porque o espectro é geométrico.
    """
    f = lambda k: integrar(-0.5 * k * k, 0.0, x0, **kw)[2]
    for _ in range(iteracoes):
        m = math.sqrt(lo * hi)
        if f(m) > nivel:
            lo = m
        else:
            hi = m
    return math.sqrt(lo * hi)


def torre(n_estados=5, x0=0.0, **kw):
    """Os primeiros kappa_n. Devolve lista."""
    ks, lo, hi = [], 1e-9, 5.0
    for n in range(n_estados):
        k = kappa_n(n, x0, lo, hi, **kw)
        ks.append(k)
        hi = k * 0.9
    return ks


# ------------------------------------------------------ espectro com armadilha
def energia_n(nivel, omega, x0=0.0, lo=None, hi=None, iteracoes=90, **kw):
    """n-ésimo nível COM armadilha isotrópica de frequência omega."""
    if lo is None:
        lo = -2.0 * math.exp(-2 * x0)      # bem abaixo: perto da parede
    if hi is None:
        hi = 6.0 * omega                    # bem acima: níveis de armadilha
    f = lambda E: integrar(E, omega, x0, **kw)[2]
    for _ in range(iteracoes):
        m = 0.5 * (lo + hi)
        if f(m) > nivel:
            hi = m
        else:
            lo = m
    return 0.5 * (lo + hi)


def n_efimov_previsto(omega, x0=0.0, s0=S0):
    """Quantos estados de Efimov cabem entre a parede e o comprimento do oscilador.

    A torre vive entre rho_0 e l_ho = 1/sqrt(omega). Cada degrau multiplica rho
    por e^(pi/s0), logo:

        N ~ (s0/pi) * ln(l_ho / rho_0)
    """
    l_ho = 1.0 / math.sqrt(omega)
    rho0 = math.exp(x0)
    if l_ho <= rho0:
        return 0.0
    return (s0 / math.pi) * math.log(l_ho / rho0)
