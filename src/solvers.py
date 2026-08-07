# -*- coding: utf-8 -*-
"""
Integração da equação radial de energia zero, onda-s:

    u''(r) = 2 V(r) u(r)          (unidades adimensionais: hbar = m_r = 1)

Dois métodos (Seção 3.1 do artigo):

  Método A - diferença central de 2a ordem (Eq. 99):
      u_{i+1} = 2 u_i - u_{i-1} + 2 (dr)^2 V(r_i) u_i

  Método B - Numerov (Eq. 101), para y'' = -xi(x) y + s(x) com s = 0 e
  xi(r) = -2 V(r):
      u_{i+1} = [ 2 u_i (1 - 5 h2 xi_i) - u_{i-1} (1 + h2 xi_{i-1}) ]
                / (1 + h2 xi_{i+1}),    h2 = (dr)^2/12

Condições de contorno: u(r_min) = 0 e u(r_min + dr) = 1 (a normalização
é livre, a equação é linear).  Integramos até R_match + dr, onde
R_match = r_min + N dr >= R é o primeiro ponto de grade além do alcance
(lá V ~ 0, então casar em R_match em vez de R não introduz erro).

Para o LJ o caroço repulsivo faz u crescer exponencialmente; usamos
re-escalonamento periódico (linearidade) para evitar overflow.
"""
import math
import numpy as np

ESCALA_MAX = 1e250   # limiar de re-escalonamento contra overflow


def grade(pot, dr):
    """Grade radial: r_i = r_min + i dr_eff, i = 0..N+1, com r_N = R exato.

    O passo pedido dr é levemente ajustado para dr_eff = (R - r_min)/N,
    de modo que a borda do potencial caia exatamente sobre um ponto da
    grade.  Isso elimina o erro de desalinhamento O(dr) na borda do poço
    esférico (sawtooth em função de mu) e não altera nada nos potenciais
    suaves.  O casamento com g0 = 1 - r/a é feito em r_N = R.
    """
    N = int(math.ceil((pot.R - pot.r_min) / dr - 1e-12))
    if N < 20:
        raise ValueError("Grade grossa demais: diminua dr.")
    dr_eff = (pot.R - pot.r_min) / N
    r = pot.r_min + dr_eff * np.arange(N + 2)
    return r, N, dr_eff


def integrar(pot, dr, metodo="numerov"):
    """Integra u''=2Vu.  Devolve (r, u, N, V, dr_eff) com r[N] = R."""
    r, N, dr = grade(pot, dr)          # dr passa a ser o passo efetivo
    V = pot.V(r)
    u = np.empty(N + 2)
    u[0] = 0.0
    u[1] = 1.0

    if metodo == "central":
        c = 2.0 * dr * dr
        for i in range(1, N + 1):
            u[i + 1] = 2.0 * u[i] - u[i - 1] + c * V[i] * u[i]
            if abs(u[i + 1]) > ESCALA_MAX:            # re-escalona (linearidade)
                u[: i + 2] /= ESCALA_MAX
    elif metodo == "numerov":
        xi = -2.0 * V
        h2 = dr * dr / 12.0
        for i in range(1, N + 1):
            u[i + 1] = (2.0 * u[i] * (1.0 - 5.0 * h2 * xi[i])
                        - u[i - 1] * (1.0 + h2 * xi[i - 1])) \
                       / (1.0 + h2 * xi[i + 1])
            if abs(u[i + 1]) > ESCALA_MAX:
                u[: i + 2] /= ESCALA_MAX
    else:
        raise ValueError(f"método desconhecido: {metodo}")

    return r, u, N, V, dr
