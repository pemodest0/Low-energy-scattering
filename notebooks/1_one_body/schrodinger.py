# -*- coding: utf-8 -*-
"""
Textbook Schrodinger problems solved with the machinery of the laboratory:
BOUND STATES (eigenvalues and eigenfunctions) and TUNNELLING, each one
checked against a closed form.

Equation solved (units with hbar = m = 1):

    u''(x) = 2 [V(x) - E] u(x)

- Bound states: shooting. Set u(x0) = 0, march with Numerov up to x1;
  E is an eigenvalue when u(x1) = 0. Sweep E, find the sign changes of
  u(x1; E), refine with brentq. The node count labels n.
- Tunnelling: transfer matrix built from TWO real solutions integrated
  across the barrier; T = 4/|M11 + M22 + i(k M12 - M21/k)|^2.

Systems with an exact answer, which is what makes them worth solving:
  harmonic oscillator  V = x^2/2       -> E_n = n + 1/2
  hydrogen (radial)    V = -1/r        -> E_n = -1/(2n^2)   (atomic units)
  Morse                V = D(1-e^{-a(r-re)})^2 - D
                       -> E_n = -D + w(n+1/2) - w^2(n+1/2)^2/(4D),
                          w = a sqrt(2D)
  square barrier       T in closed form, sinh^2(kappa L)

Why these three and not others: the oscillator is the TRAP of a cold atomic
gas (3D->2D confinement is an anisotropic oscillator); hydrogen uses the SAME
radial equation u(r) as the scattering problem; the Morse potential is the
close relative of the He2 / Lennard-Jones well.
"""
import math

import numpy as np
from scipy.optimize import brentq


# ------------------------------------------------------- Numerov march
def _march(Veff, E, x0, x1, dx):
    """u'' = 2(Veff-E)u with u(x0)=0, u(x0+dx)=1e-8; returns (x, u)."""
    n = int(round((x1 - x0) / dx)) + 1
    x = x0 + dx * np.arange(n)
    W = -2.0 * (Veff(x) - E)                 # so that u'' = -W u
    u = np.empty(n)
    u[0], u[1] = 0.0, 1e-8
    h2 = dx * dx / 12.0
    for i in range(1, n - 1):
        u[i + 1] = (2.0 * u[i] * (1.0 - 5.0 * h2 * W[i])
                    - u[i - 1] * (1.0 + h2 * W[i - 1])) / (1.0 + h2 * W[i + 1])
        if abs(u[i + 1]) > 1e250:             # rescale before overflow
            u[: i + 2] /= 1e250
    return x, u


# --------------------------------------------------- eigenvalues (shooting)
def _nodes_and_end(V, E, x0, x1, dx):
    _, u = _march(V, E, x0, x1, dx)
    s = np.sign(u[np.abs(u) > 0])
    return int(np.count_nonzero(s[1:] * s[:-1] < 0)), u[-1]


def eigenvalue(V, n, x0, x1, E_min, E_max, dx=2e-3):
    """E_n, the state with n nodes, by bisection on the NODE COUNT + brentq.

    The oscillation theorem for Sturm-Liouville problems: the number of nodes
    of u(x; E) in (x0, x1) equals the number of eigenvalues below E. So bisect
    until the n -> n+1 transition is cornered; inside that bracket u(x1; E)
    changes sign exactly once, which is what brentq needs. About 40
    integrations per state.
    """
    lo, hi = E_min, E_max
    for _ in range(60):                       # corner the node transition
        mid = 0.5 * (lo + hi)
        nodes, _ = _nodes_and_end(V, mid, x0, x1, dx)
        if nodes <= n:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-9 * max(1.0, abs(hi)):
            break
    f = lambda E: _march(V, E, x0, x1, dx)[1][-1]
    # u(x1) changes sign on crossing the eigenvalue, somewhere in [lo-, hi+]
    eps = max(1e-8, 1e-6 * (E_max - E_min))
    a, b = lo - eps, hi + eps
    fa, fb = f(a), f(b)
    if fa * fb > 0:                           # safeguard: widen the bracket
        a, b = lo - 100 * eps, hi + 100 * eps
    return brentq(f, a, b, xtol=1e-12, rtol=1e-12)


def eigenvalues(V, x0, x1, E_min, E_max, n_states=4, dx=2e-3):
    """List of (E_n, n) for n = 0 .. n_states-1."""
    return [(eigenvalue(V, n, x0, x1, E_min, E_max, dx), n)
            for n in range(n_states)]


def eigenfunction(V, E, x0, x1, dx=1e-3):
    """Eigenfunction normalised to L2 norm 1, at energy E."""
    x, u = _march(V, E, x0, x1, dx)
    u = u / math.sqrt(np.trapezoid(u * u, x))
    return x, u


# ------------------------------------------------------- the four systems
def V_oscillator(x):
    return 0.5 * np.asarray(x, dtype=float)**2


def E_oscillator(n):
    return n + 0.5


def V_hydrogen(l=0):
    """Effective radial potential of H in atomic units: -1/r + l(l+1)/2r^2."""
    def V(r):
        r = np.asarray(r, dtype=float)
        rs = np.where(r > 0, r, 1e-12)
        return -1.0 / rs + l * (l + 1) / (2.0 * rs**2)
    return V


def E_hydrogen(n):
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
    """How many bound levels the well holds: n <= sqrt(2D)/a - 1/2."""
    return int(math.floor(math.sqrt(2.0 * D) / a - 0.5))


# ----------------------------------------------------------- tunnelling
def transmission(V, E, x0, x1, dx=1e-4):
    """T(E) by transfer matrix: integrate two real solutions, with initial
    conditions (1,0) and (0,1), across [x0,x1], matching plane waves outside
    (V=0, k=sqrt(2E)). T = 4/|M11 + M22 + i(k M12 - M21/k)|^2."""
    n = int(round((x1 - x0) / dx)) + 1
    x = x0 + dx * np.arange(n)
    W = -2.0 * (V(x) - E)
    h2 = dx * dx / 12.0

    def march(u0, u1):
        u = np.empty(n)
        u[0], u[1] = u0, u1
        for i in range(1, n - 1):
            u[i + 1] = (2.0 * u[i] * (1.0 - 5.0 * h2 * W[i])
                        - u[i - 1] * (1.0 + h2 * W[i - 1])) / (1.0 + h2 * W[i + 1])
        return u

    # solution A: y(x0)=1, y'(x0)=0 ; solution B: y(x0)=0, y'(x0)=1
    # (the second point comes from a Taylor step using y'' = -W y)
    yA = march(1.0, 1.0 + 0.5 * dx * dx * (-W[0]))
    yB = march(0.0, dx * (1.0 + dx * dx * (-W[1]) / 6.0))

    # derivative AT x1 by a one-sided 5-point stencil, order dx^4:
    # f'(x_N) = (25 f_N - 48 f_{N-1} + 36 f_{N-2} - 16 f_{N-3} + 3 f_{N-4})/(12 dx)
    def deriv(y):
        return (25 * y[-1] - 48 * y[-2] + 36 * y[-3]
                - 16 * y[-4] + 3 * y[-5]) / (12 * dx)

    M11, M21 = yA[-1], deriv(yA)
    M12, M22 = yB[-1], deriv(yB)
    k = math.sqrt(2.0 * E)
    den = complex(M11 + M22, k * M12 - M21 / k)
    return 4.0 / abs(den)**2


def V_barrier(V0=5.0, L=1.0):
    def V(x):
        x = np.asarray(x, dtype=float)
        return np.where((x >= 0) & (x <= L), V0, 0.0)
    return V


def T_barrier_exact(E, V0=5.0, L=1.0):
    """Square barrier, closed form. Below the top the answer carries sinh,
    above it sin, and exactly at E = V0 both degenerate to the same limit."""
    if E == V0:
        return 1.0 / (1.0 + V0 * L * L / 2.0)
    if E < V0:
        kap = math.sqrt(2.0 * (V0 - E))
        s = math.sinh(kap * L)
        return 1.0 / (1.0 + V0**2 * s * s / (4.0 * E * (V0 - E)))
    kk = math.sqrt(2.0 * (E - V0))
    s = math.sin(kk * L)
    return 1.0 / (1.0 + V0**2 * s * s / (4.0 * E * (E - V0)))
