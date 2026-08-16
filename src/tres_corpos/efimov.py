# -*- coding: utf-8 -*-
"""
Efimov na unitariedade — a torre de trímeros saindo do NOSSO código.

Física (três bósons idênticos, |a| -> infinito):

1) O parâmetro universal s0 é a raiz da equação transcendental de
   Efimov (ver Braaten & Hammer, Phys. Rep. 428, 259 (2006)):

       s0 * cosh(pi s0 / 2)  =  (8/sqrt(3)) * sinh(pi s0 / 6)

   Aqui NÃO usamos o valor tabelado 1.00624: nós o RESOLVEMOS (brentq).

2) No canal hiperesférico mais baixo, a equação hiperradial de energia
   E = -hbar^2 kappa^2 / (2m) é

       f''(R) + [ (s0^2 + 1/4)/R^2 - kappa^2 ] f(R) = 0,

   com R o hiper-raio.  O potencial -1/R^2 é invariante de escala
   (colapso de Thomas): é preciso UMA condição de curto alcance — o
   "parâmetro de três corpos".  Usamos a mais simples: parede dura em
   R = R0, f(R0) = 0.

3) Com a substituição  x = ln(R/R0),  f = e^{x/2} g(x), a equação vira

       g''(x) + [ s0^2 - (kappa R0)^2 e^{2x} ] g(x) = 0,

   um problema de autovalor 1D numa grade UNIFORME em x (cobre décadas
   de R sem esforço).  Os autovalores kappa_n formam a torre geométrica:

       kappa_n / kappa_{n+1} -> e^{pi/s0} ~ 22.694
       E_n / E_{n+1}        -> e^{2 pi/s0} ~ 515.03

   O fator 22.7 tem que EMERGIR do cálculo — é o teste.
"""
import math

import numpy as np
from scipy.optimize import brentq


# ------------------------------------------------ 1) s0 do zero
def s0_universal():
    """Resolve s cosh(pi s/2) = (8/sqrt(3)) sinh(pi s/6) para s > 0."""
    F = lambda s: s * math.cosh(math.pi * s / 2) \
        - (8.0 / math.sqrt(3.0)) * math.sinh(math.pi * s / 6)
    return brentq(F, 0.5, 1.5, xtol=1e-14)


S0 = s0_universal()          # ~ 1.0062378...


# ------------------------------- 2-3) torre por shooting (Numerov em x)
def _g_final(kR0, s0=S0, dx=2e-3, margem=5.0):
    """Integra g'' = -[s0^2 - (kR0)^2 e^{2x}] g de x=0 (g=0) até depois
    do ponto de retorno; devolve o sinal/valor de g no fim (a solução
    ligada é a que decai: autovalor quando g_final cruza zero)."""
    x_ret = math.log(s0 / kR0) if kR0 < s0 else 0.0
    x_max = x_ret + margem
    n = int(x_max / dx) + 2
    x = dx * np.arange(n)
    W = s0**2 - (kR0 * np.exp(x))**2          # g'' = -W g  (xi = W)
    g = np.empty(n)
    g[0], g[1] = 0.0, 1e-8
    h2 = dx * dx / 12.0
    for i in range(1, n - 1):
        g[i + 1] = (2.0 * g[i] * (1.0 - 5.0 * h2 * W[i])
                    - g[i - 1] * (1.0 + h2 * W[i - 1])) / (1.0 + h2 * W[i + 1])
        if abs(g[i + 1]) > 1e250:
            g[: i + 2] /= 1e250
    return g[-1]


def torre(n_niveis=4, s0=S0, kR0_max=3.0, dx=2e-3):
    """Autovalores kappa_n R0 (do mais fundo ao mais raso) por varredura
    geométrica + bisseção nos cruzamentos de sinal de g_final(kappa)."""
    ks = np.geomspace(kR0_max, kR0_max * math.exp(-math.pi / s0 * (n_niveis + 1.5)),
                      60 * (n_niveis + 2))
    sinais = [math.copysign(1.0, _g_final(k, s0, dx)) for k in ks]
    raizes = []
    for i in range(len(ks) - 1):
        if sinais[i] * sinais[i + 1] < 0:
            r = brentq(lambda k: _g_final(k, s0, dx), ks[i + 1], ks[i],
                       xtol=1e-12, rtol=1e-10)
            raizes.append(r)
            if len(raizes) >= n_niveis:
                break
    return np.array(raizes)          # kappa_n R0, n = 0, 1, ...


def razoes(kappas):
    """Razões consecutivas kappa_n/kappa_{n+1} (esperado: e^{pi/s0})."""
    return kappas[:-1] / kappas[1:]


if __name__ == "__main__":
    print(f"s0 (resolvido do zero)  : {S0:.10f}   (literatura: 1.00624)")
    print(f"e^(pi/s0)  esperado     : {math.exp(math.pi / S0):.4f}")
    ks = torre(5)
    print(f"kappa_n R0              : {[f'{k:.6e}' for k in ks]}")
    rz = razoes(ks)
    print(f"razões kappa_n/kappa_n+1: {[f'{r:.4f}' for r in rz]}")
    print(f"razões de ENERGIA       : {[f'{r*r:.2f}' for r in rz]}  (esperado ~515.03)")
