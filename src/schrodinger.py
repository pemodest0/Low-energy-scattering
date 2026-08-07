# -*- coding: utf-8 -*-
"""
Equação de Schrödinger de livro-texto com a máquina do laboratório:
ESTADOS LIGADOS (autovalores/autofunções) e TUNELAMENTO, sempre
numérico × analítico.

Equação resolvida (unidades com hbar = m = 1):

    u''(x) = 2 [V(x) - E] u(x)

- Estados ligados: "shooting" — u(x0) = 0, marcha com Numerov até x1;
  E é autovalor quando u(x1) = 0.  Varremos E, localizamos as trocas de
  sinal de u(x1; E) e refinamos com brentq.  O número de nós rotula n.
- Tunelamento: matriz de transferência construída com DUAS soluções
  reais integradas através da barreira; T = 4/|M11+M22 + i(kM12 - M21/k)|².

Sistemas com resposta exata (validação nos testes):
  oscilador harmônico  V = x²/2        -> E_n = n + 1/2
  hidrogênio (radial)  V = -1/r        -> E_n = -1/(2n²)   (unid. atômicas)
  Morse                V = D(1-e^{-a(r-re)})² - D
                       -> E_n = -D + w(n+½) - w²(n+½)²/(4D),  w = a√(2D)
  barreira retangular  T analítico com sinh²(kappa L)

Conexões com o mestrado: o oscilador é a ARMADILHA dos átomos frios
(confinamento 3D->2D = oscilador anisotrópico); o hidrogênio usa a MESMA
equação radial u(r) do espalhamento; o Morse é o primo do He2/LJ.
"""
import math

import numpy as np
from scipy.optimize import brentq


# ------------------------------------------------ marcha de Numerov
def _integra(Veff, E, x0, x1, dx):
    """u'' = 2(Veff-E)u com u(x0)=0, u(x0+dx)=1e-8; devolve (x, u)."""
    n = int(round((x1 - x0) / dx)) + 1
    x = x0 + dx * np.arange(n)
    W = -2.0 * (Veff(x) - E)                 # u'' = -W u
    u = np.empty(n)
    u[0], u[1] = 0.0, 1e-8
    h2 = dx * dx / 12.0
    for i in range(1, n - 1):
        u[i + 1] = (2.0 * u[i] * (1.0 - 5.0 * h2 * W[i])
                    - u[i - 1] * (1.0 + h2 * W[i - 1])) / (1.0 + h2 * W[i + 1])
        if abs(u[i + 1]) > 1e250:
            u[: i + 2] /= 1e250
    return x, u


def _nos(u):
    s = np.sign(u[np.abs(u) > 0])
    return int(np.count_nonzero(s[1:] * s[:-1] < 0)) - 1   # ignora u(x1)~0


# ------------------------------------------------ autovalores (shooting)
def _conta_nos_E(V, E, x0, x1, dx):
    _, u = _integra(V, E, x0, x1, dx)
    s = np.sign(u[np.abs(u) > 0])
    return int(np.count_nonzero(s[1:] * s[:-1] < 0)), u[-1]


def autovalor_n(V, n, x0, x1, E_min, E_max, dx=2e-3):
    """E_n (estado com n nós) por bisseção na CONTAGEM DE NÓS + brentq.

    Teorema da oscilação: o nº de nós de u(x; E) em (x0, x1) é o nº de
    autovalores abaixo de E.  Bisseção até encurralar a transição
    n -> n+1 nós; dentro desse intervalo u(x1; E) troca de sinal
    exatamente uma vez -> brentq.  ~40 integrações por estado.
    """
    lo, hi = E_min, E_max
    for _ in range(60):                       # encurrala a transição de nós
        mid = 0.5 * (lo + hi)
        nos, _ = _conta_nos_E(V, mid, x0, x1, dx)
        if nos <= n:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-9 * max(1.0, abs(hi)):
            break
    f = lambda E: _integra(V, E, x0, x1, dx)[1][-1]
    # u(x1) muda de sinal ao cruzar o autovalor dentro de [lo-, hi+]
    eps = max(1e-8, 1e-6 * (E_max - E_min))
    a, b = lo - eps, hi + eps
    fa, fb = f(a), f(b)
    if fa * fb > 0:                           # salvaguarda: alarga um pouco
        a, b = lo - 100 * eps, hi + 100 * eps
        fa, fb = f(a), f(b)
    return brentq(f, a, b, xtol=1e-12, rtol=1e-12)


def autovalores(V, x0, x1, E_min, E_max, n_estados=4, dx=2e-3):
    """Lista [(E_n, n)] para n = 0 .. n_estados-1."""
    return [(autovalor_n(V, n, x0, x1, E_min, E_max, dx), n)
            for n in range(n_estados)]


def autofuncao(V, E, x0, x1, dx=1e-3):
    """Autofunção normalizada (norma L2 = 1) para energia E."""
    x, u = _integra(V, E, x0, x1, dx)
    u = u / math.sqrt(np.trapezoid(u * u, x))
    return x, u


# ----------------------------------------------------- sistemas canônicos
def V_oscilador(x):
    return 0.5 * np.asarray(x, dtype=float)**2


def E_oscilador(n):
    return n + 0.5


def V_hidrogenio(l=0):
    """V efetivo radial do H em unidades atômicas: -1/r + l(l+1)/2r²."""
    def V(r):
        r = np.asarray(r, dtype=float)
        rs = np.where(r > 0, r, 1e-12)
        return -1.0 / rs + l * (l + 1) / (2.0 * rs**2)
    return V


def E_hidrogenio(n):
    return -0.5 / n**2


def V_morse(D=10.0, a=1.0, re=2.0):
    def V(r):
        e = np.exp(-a * (np.asarray(r, dtype=float) - re))
        return D * (1.0 - e)**2 - D
    return V


def E_morse(n, D=10.0, a=1.0):
    w = a * math.sqrt(2.0 * D)
    return -D + w * (n + 0.5) - w**2 * (n + 0.5)**2 / (4.0 * D)


def n_max_morse(D=10.0, a=1.0):
    """Número de níveis ligados: n <= sqrt(2D)/a - 1/2."""
    return int(math.floor(math.sqrt(2.0 * D) / a - 0.5))


# ------------------------------------------------------- tunelamento
def transmissao(V, E, x0, x1, dx=1e-4):
    """T(E) por matriz de transferência: integra duas soluções reais
    (ICs (1,0) e (0,1)) através de [x0,x1], com ondas planas fora
    (V=0, k=sqrt(2E)).  T = 4/|M11+M22 + i(k M12 - M21/k)|²."""
    n = int(round((x1 - x0) / dx)) + 1
    x = x0 + dx * np.arange(n)
    W = -2.0 * (V(x) - E)
    h2 = dx * dx / 12.0

    def marcha(u0, u1):
        u = np.empty(n)
        u[0], u[1] = u0, u1
        for i in range(1, n - 1):
            u[i + 1] = (2.0 * u[i] * (1.0 - 5.0 * h2 * W[i])
                        - u[i - 1] * (1.0 + h2 * W[i - 1])) / (1.0 + h2 * W[i + 1])
        return u

    # solução A: y(x0)=1, y'(x0)=0 ; solução B: y(x0)=0, y'(x0)=1
    # (derivada inicial via passo de Taylor com y'' = -W y)
    yA = marcha(1.0, 1.0 + 0.5 * dx * dx * (-W[0]))
    yB = marcha(0.0, dx * (1.0 + dx * dx * (-W[1]) / 6.0))
    # derivada EM x1 por estêncil unilateral de 5 pontos, O(dx^4)
    # f'(x_N) = (25 f_N - 48 f_{N-1} + 36 f_{N-2} - 16 f_{N-3} + 3 f_{N-4})/(12 dx)
    def deriv(y):
        return (25 * y[-1] - 48 * y[-2] + 36 * y[-3]
                - 16 * y[-4] + 3 * y[-5]) / (12 * dx)
    M11, M21 = yA[-1], deriv(yA)
    M12, M22 = yB[-1], deriv(yB)
    k = math.sqrt(2.0 * E)
    den = complex(M11 + M22, k * M12 - M21 / k)
    return 4.0 / abs(den)**2


def V_barreira(V0=5.0, L=1.0):
    def V(x):
        x = np.asarray(x, dtype=float)
        return np.where((x >= 0) & (x <= L), V0, 0.0)
    return V


def T_barreira_analitico(E, V0=5.0, L=1.0):
    """Barreira retangular exata (E < V0 usa sinh; E > V0 usa sin)."""
    if E == V0:
        k2L = math.sqrt(2.0 * E) * L
        return 1.0 / (1.0 + (k2L / 2.0)**2 * 0 + (math.sqrt(2*E)*L/2)**2)  # limite
    if E < V0:
        kap = math.sqrt(2.0 * (V0 - E))
        s = math.sinh(kap * L)
        return 1.0 / (1.0 + V0**2 * s * s / (4.0 * E * (V0 - E)))
    kk = math.sqrt(2.0 * (E - V0))
    s = math.sin(kk * L)
    return 1.0 / (1.0 + V0**2 * s * s / (4.0 * E * (E - V0)))
