# -*- coding: utf-8 -*-
"""
Resultados analíticos usados para validar o código.

Poço esférico:
    a  (Eq. 80):  a = R [1 - tan(sqrt(2 v))/sqrt(2 v)]
    r0 (Eq. 92):  forma fechada obtida integrando a Eq. (56) com a solução
                  exata (Eq. 90); implementamos a integral em forma fechada,
                  equivalente à Eq. (92) do artigo.
    Limiar de estado ligado: v > pi^2/8.

Pöschl-Teller modificado (v = lambda(lambda-1)/2):
    a  (Eq. 117): a mu = (pi/2) cot(pi lambda/2) + gamma + Psi(lambda)
    Na unitariedade (lambda = 2, v = 1): r0 = 2/mu, e a solução radial é
    u0(r) = tanh(mu r)/tanh(mu R)  (Eq. 118).

Energias de estado ligado (Tabela 1):
    E_zr (Eq. 93): E = -hbar^2/(2 m_r a^2)
    E_fr (Eq. 95): 1/a = kappa - r0 kappa^2/2  ->  E = -hbar^2 kappa^2/(2 m_r)

digamma via scipy.special (v1); a v0 tinha implementação caseira.
"""
import math
import numpy as np

from .constantes import EULER_GAMMA


# ------------------------------------------------------------------ digamma
# v1: usamos scipy.special.digamma (validada e vetorizada) no lugar da
# implementação caseira da v0 (recorrência + série assintótica).  O wrapper
# preserva a assinatura usada no resto do código e nos testes.
from scipy.special import digamma as _digamma_scipy


def digamma(x):
    """Psi(x) via scipy.special.digamma (x > 0 no nosso uso)."""
    if x <= 0:
        raise ValueError("digamma implementada apenas para x > 0")
    return float(_digamma_scipy(x))


# ------------------------------------------------------- poço esférico
V_LIMIAR_POCO = math.pi**2 / 8.0     # limiar do 1º estado ligado (Eq. 86)


def a_poco(v, R=1.0):
    """Eq. (80): a = R [1 - tan(x)/x], x = sqrt(2v) = k0 R."""
    x = math.sqrt(2.0 * v)
    return R * (1.0 - math.tan(x) / x)


def r0_poco(v, R=1.0):
    """r0 = 2 int_0^R [g0^2 - u0^2] dr em forma fechada (equiv. à Eq. 92),
    com u0 da Eq. (90): u0 = (1 - R/a) sin(k0 r)/sin(k0 R), g0 = 1 - r/a."""
    a = a_poco(v, R)
    k0 = math.sqrt(2.0 * v) / R
    I1 = R - R**2 / a + R**3 / (3.0 * a**2)                       # int g0^2
    I2 = (1.0 - R / a)**2 / math.sin(k0 * R)**2 \
        * (R / 2.0 - math.sin(2.0 * k0 * R) / (4.0 * k0))         # int u0^2
    return 2.0 * (I1 - I2)


def u_poco(r, v, R=1.0):
    """Solução exata de energia zero (Eq. 90), normalizada a g0 fora."""
    a = a_poco(v, R)
    k0 = math.sqrt(2.0 * v) / R
    r = np.asarray(r, dtype=float)
    dentro = (1.0 - R / a) * np.sin(k0 * r) / math.sin(k0 * R)
    fora = 1.0 - r / a
    return np.where(r < R, dentro, fora)


# ------------------------------------------- Pöschl-Teller modificado
def lambda_mpt(v):
    """v = lambda(lambda-1)/2  ->  lambda = (1 + sqrt(1+8v))/2 (raiz > 1)."""
    return 0.5 * (1.0 + math.sqrt(1.0 + 8.0 * v))


def a_mpt(v, mu=1.0):
    """Eq. (117): a mu = (pi/2) cot(pi lambda/2) + gamma + Psi(lambda)."""
    lam = lambda_mpt(v)
    cot = math.cos(math.pi * lam / 2.0) / math.sin(math.pi * lam / 2.0)
    return ((math.pi / 2.0) * cot + EULER_GAMMA + digamma(lam)) / mu


def r0_mpt_unitario(mu=1.0):
    """Na unitariedade (v = 1): r0 = 2/mu (integrando a Eq. 118 na Eq. 56)."""
    return 2.0 / mu


def u_mpt_unitario(r, mu=1.0, R=None):
    """Eq. (118): u0(r) = tanh(mu r)/tanh(mu R); com R -> infinito, tanh(mu R) = 1."""
    r = np.asarray(r, dtype=float)
    denom = math.tanh(mu * R) if R is not None else 1.0
    return np.tanh(mu * r) / denom


# --------------------------------------------- energias de estado ligado
def energia_zr(a, h2_2mr):
    """Eq. (93): E_zr = -hbar^2/(2 m_r a^2).  h2_2mr = hbar^2/(2 m_r)."""
    return -h2_2mr / a**2


def energia_fr(a, r0, h2_2mr):
    """Eq. (95): resolve 1/a = kappa - r0 kappa^2/2 (raiz física, que
    tende a 1/a quando r0 -> 0) e devolve E_fr = -hbar^2 kappa^2/(2 m_r)."""
    kappa = (1.0 - math.sqrt(1.0 - 2.0 * r0 / a)) / r0
    return -h2_2mr * kappa**2
