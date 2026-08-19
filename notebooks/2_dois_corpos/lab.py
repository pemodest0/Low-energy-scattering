"""Two-body low-energy scattering: potentials, solver and tabulated data.

The notebook next to this file needs nothing else.

UNITS
    hbar = reduced mass = 1, lengths in fm, so the radial equation is
        u'' = 2 (V - E) u
    and physical energies come from  E_phys = 2 * (hbar^2 / 2 mu) * E_code.

THE ONE IDEA
    The potential only acts up to a radius R. Past R we have V = 0, and there
    the equation has a solution we already know:

        at E = 0:  u'' = 0        so u is a STRAIGHT LINE
        at E < 0:  u'' = -2 E u   so u is exp(-kappa r), kappa = sqrt(-2E)

    So we never integrate past R. We integrate from 0 to R, read off the value
    and the slope there, and glue the known formula on. All three quantities we
    want come out of that gluing:

        a  = where the straight line crosses zero
        r0 = integral of (line - true solution) between 0 and R
        E  = the energy that makes the solution glue onto the exponential

    The integration is Numerov, order h^4, on a grid uniform in x = ln r.
    See numerov() for why the logarithm, and why the grid ends exactly on R.

SOURCES
    [1] Macedo-Lima, M.; Madeira, L. Rev. Bras. Ensino Fis. 45, e20230079 (2023)
        doi:10.1590/1806-9126-RBEF-2023-0079
        The four potentials (Eqs. 74, 116, 118, 120), the targets and published
        parameters (Tables 2, 3, 4), and the square-well closed forms (80, 92).
    [2] Hackenburg, R. W. Phys. Rev. C 73, 044002 (2006)
        doi:10.1103/PhysRevC.73.044002
        Deuteron a, r0 and binding energy. Ref. [23] of [1].
    [3] Cencek, W. et al. J. Chem. Phys. 136, 224303 (2012)
        doi:10.1063/1.4712218
        Helium-4 dimer a, r0 and binding energy. Ref. [22] of [1].
"""
import math
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import brentq

# Below V_ZERO we treat the potential as gone. Above V_CORE the Lennard-Jones
# core is too steep to integrate through, and the wavefunction is zero there
# anyway, so we start outside it.
V_ZERO = 1e-12
V_CORE = 1e5


# =============================================================================
#  The four potentials
#
#  Each class holds two things: the formula for V(r), and where it stops (R).
#  R is written in closed form -- we solve |V(R)| = V_ZERO by hand for each
#  shape, so no search is needed.
# =============================================================================
class Well:
    """V = -v mu^2 inside r < R = 1/mu, zero outside.  Eq. (74) of [1].

    Discontinuous on purpose: it is the hard case for any numerical method.
    """
    name = "Square well"

    def __init__(self, v, mu):
        self.v = v
        self.mu = mu
        self.r_min = 0.0
        self.R = 1.0 / mu               # exact: the edge is where the well ends

    def V(self, r):
        depth = -self.v * self.mu ** 2
        return np.where(r <= self.R, depth, 0.0)


class PoschlTeller:
    """V = -v mu^2 / cosh^2(mu r).  Eq. (116) of [1].  Smooth and solvable."""
    name = "Poschl-Teller"

    def __init__(self, v, mu):
        self.v = v
        self.mu = mu
        self.r_min = 0.0
        # v mu^2 / cosh^2(mu R) = V_ZERO
        self.R = math.acosh(math.sqrt(v * mu ** 2 / V_ZERO)) / mu

    def V(self, r):
        return -self.v * self.mu ** 2 / np.cosh(self.mu * r) ** 2


class Gaussian:
    """V = -v mu^2 exp(-mu^2 r^2).  Eq. (118) of [1].  Tail dies very fast."""
    name = "Gaussian"

    def __init__(self, v, mu):
        self.v = v
        self.mu = mu
        self.r_min = 0.0
        # v mu^2 exp(-mu^2 R^2) = V_ZERO
        self.R = math.sqrt(math.log(v * mu ** 2 / V_ZERO)) / mu

    def V(self, r):
        return -self.v * self.mu ** 2 * np.exp(-(self.mu * r) ** 2)


class LennardJones:
    """V = (1/2) (C12/r^12 - C6/r^6).  Eq. (120) of [1].  Hard core + vdW tail.

    The 1/2 is a CONVENTION and it is not universal: much of the literature
    writes C12/r^12 - C6/r^6 without it, which shifts C6 by a factor of two
    while the result still looks correct. [1] uses the 1/2, so we do.
    """
    name = "Lennard-Jones"

    def __init__(self, C6, C12):
        self.C6 = C6
        self.C12 = C12
        self.r_min = (0.5 * C12 / V_CORE) ** (1.0 / 12.0)   # core: 0.5 C12/r^12 = V_CORE
        self.R = (0.5 * C6 / V_ZERO) ** (1.0 / 6.0)         # tail: 0.5 C6/r^6  = V_ZERO

    def V(self, r):
        return 0.5 * (self.C12 / r ** 12 - self.C6 / r ** 6)


# The strength parameter always comes first, the scale parameter second.
POTENTIALS = {"well": Well, "mpt": PoschlTeller,
              "gauss": Gaussian, "lj": LennardJones}


# One grid size for all four potentials -- see numerov() for why a logarithmic
# grid lets a single number work here. 8001 is not arbitrary: total error is
# truncation (falls with more points) plus roundoff (grows with them), and this
# sits at the minimum. The three smooth potentials are at machine precision,
# 1e-11, and are insensitive. The Lennard-Jones is the one that sets the value:
# past ~8000 points its r0 degrades again, by 3e-5 at 32001 and 3e-4 at 64001.
POINTS = 8001


# =============================================================================
#  Numerov on a logarithmic grid
# =============================================================================
def numerov(pot, E):
    """Integrate the radial equation and return r, u, x, h, slope at the edge.

    WHY THE LOGARITHM
        Numerov needs h * sqrt(|W|) << 1 to represent the solution. On a
        uniform grid the Lennard-Jones fails that badly: its core reaches
        V = 1e5, so sqrt(2V) = 447 and the natural length there is 0.002 fm,
        while the tail forces the grid out to 125 fm. One step cannot serve
        both. Measured consequence: r0 drifted 1.74257 -> 1.74215 -> 1.74165
        as the points doubled, instead of converging.

        With r = e^x the grid is dense at small r and sparse at large r, which
        is exactly what a hard core plus a power-law tail asks for. Substituting
        u = e^(x/2) w turns the equation into

            w'' = [ 2 (V - E) e^(2x) + 1/4 ] w

        which is again u'' = W u, so the Numerov step is unchanged: only what
        goes into it changes.

    WHY THE GRID ENDS ON R
        If a step straddled the square well's discontinuity, Numerov would drop
        from order 4 to order 1 -- measured, not assumed. Ending on the edge
        means no step ever crosses it. The slope there comes from a five-point
        one-sided difference, order h^4, built from interior points only.
    """
    # Where to start: outside the hard core if there is one, otherwise close
    # enough to the origin that the piece we skip is analytic (see scattering).
    r_start = pot.r_min if pot.r_min > 0.0 else pot.R * 1e-6

    x = np.linspace(math.log(r_start), math.log(pot.R), POINTS)
    h = x[1] - x[0]
    r = np.exp(x)
    W = 2.0 * (pot.V(r) - E) * r * r + 0.25
    f = 1.0 - h * h * W / 12.0

    w = np.zeros(POINTS)
    if pot.r_min > 0.0:
        w[1] = h * math.exp(x[0] / 2.0)      # hard core: u = 0, u' = 1 there
    else:
        w[0] = math.sqrt(r[0])               # regular origin: u ~ r, so w ~ sqrt(r)
        w[1] = math.sqrt(r[1])
    for i in range(1, POINTS - 1):
        w[i + 1] = ((12.0 - 10.0 * f[i]) * w[i] - f[i - 1] * w[i - 1]) / f[i + 1]

    u = np.exp(x / 2.0) * w

    # du/dr = e^(-x/2) (w/2 + dw/dx), with dw/dx from the one-sided formula.
    dwdx = (25.0 * w[-1] - 48.0 * w[-2] + 36.0 * w[-3]
            - 16.0 * w[-4] + 3.0 * w[-5]) / (12.0 * h)
    slope = math.exp(-x[-1] / 2.0) * (0.5 * w[-1] + dwdx)

    return r, u, x, h, slope, r_start


# =============================================================================
#  Zero energy: the scattering length and the effective range
# =============================================================================
def scattering(pot):
    """Return a, r0 and the number of bound states, all at E = 0."""
    r, u, x, h, slope, r_start = numerov(pot, 0.0)

    # Past the edge V = 0, so u is the straight line C (1 - r/a). Where that
    # line crosses zero is the scattering length -- exact, not a fit.
    if slope == 0.0:
        a = math.inf        # exact unitarity: the line is flat, it never crosses
    else:
        a = pot.R - u[-1] / slope

    # Rescale u so it joins the line at the edge, then compare the two.
    u = u * (1.0 - pot.R / a) / u[-1]
    line = 1.0 - r / a

    # r0 = 2 * integral of (line^2 - u^2) dr, Eq. (56) of [1]. On this grid
    # dr = r dx, hence the extra factor of r.
    r0 = 2.0 * simpson((line ** 2 - u ** 2) * r, x=x)

    # The grid starts at r_start, so add the piece below it by hand. There u is
    # either zero (inside a hard core) or proportional to r, and in both cases
    # its contribution is O(r_start^3) while the line's is O(r_start).
    s = r_start
    r0 = r0 + 2.0 * (s - s ** 2 / a + s ** 3 / (3.0 * a ** 2))

    # Every zero of the E = 0 solution is one bound state (Levinson's theorem).
    nodes = 0
    for i in range(1, POINTS - 1):
        if u[i] * u[i + 1] < 0.0:
            nodes = nodes + 1
    if pot.R < a < 1e4:
        # The line crosses zero at r = a, outside the potential. When a is
        # effectively infinite that zero sits at infinity and does not count.
        nodes = nodes + 1

    return a, r0, nodes


# =============================================================================
#  Negative energy: the bound state
# =============================================================================
def bound_energy(pot, guess):
    """Ground-state energy, in code units, refined around `guess`.

    Past the edge the solution has to be the decaying exponential exp(-kappa r)
    with kappa = sqrt(-2E), so the right energy is the one where

        slope(R) + kappa * u(R) = 0

    We add instead of writing slope/u + kappa because u(R) can be zero -- the
    deuteron has a node -- and dividing would invent a pole.

    `guess` is the finite-range formula, already within about 1%, so there is
    no search here: the physics gives us the bracket.
    """
    def gap(E):
        r, u, x, h, slope, r_start = numerov(pot, E)
        kappa = math.sqrt(-2.0 * E)
        return (slope + kappa * u[-1]) / max(abs(u[-1]), abs(slope), 1.0)

    return brentq(gap, 1.8 * guess, 0.3 * guess, xtol=1e-15)


# =============================================================================
#  The inverse problem: given (a, r0), find the parameters
# =============================================================================
def find_root(f, start):
    """Widen a bracket around `start` until f changes sign, then bisect."""
    low = start
    high = start
    for attempt in range(40):
        low = low / 1.3
        high = high * 1.3
        if f(low) * f(high) < 0.0:
            return brentq(f, low, high, xtol=1e-12)
    raise RuntimeError("no sign change found around %g" % start)


def match_a(name, scale, strength_guess, inv_a_target):
    """Inner loop: move the strength until 1/a hits the target, scale fixed.

    We chase 1/a and never a. `a` has poles -- one every time a new bound state
    appears -- and bisection dies on a pole but works fine on a zero. As a
    bonus, unitarity is simply 1/a = 0 and needs no special case.
    """
    def error(strength):
        pot = POTENTIALS[name](strength, scale)
        a, r0, nodes = scattering(pot)
        return 1.0 / a - inv_a_target

    return find_root(error, strength_guess)


def tune(name, a_target, r0_target, strength, scale, nodes_target=None):
    """Find the parameters of `name` that give (a_target, r0_target).

    Outer loop moves the scale so r0 matches; the inner loop above moves the
    strength so 1/a matches. One equation, one variable, at a time.

    A 2D Newton on both at once looks tidier and is worse: the level sets of
    1/a and r0 are nearly parallel over much of the plane, so its Jacobian is
    ill-conditioned exactly where you need it.

    Returns strength, scale, a, r0, nodes, ok.
    """
    if math.isinf(a_target):
        inv_a_target = 0.0
    else:
        inv_a_target = 1.0 / a_target

    for step in range(25):
        strength = match_a(name, scale, strength, inv_a_target)
        a, r0, nodes = scattering(POTENTIALS[name](strength, scale))
        if abs(r0 - r0_target) < 1e-7:
            break

        def r0_error(trial_scale):
            s = match_a(name, trial_scale, strength, inv_a_target)
            a2, r02, n2 = scattering(POTENTIALS[name](s, trial_scale))
            return r02 - r0_target

        scale = find_root(r0_error, scale)

    strength = match_a(name, scale, strength, inv_a_target)
    a, r0, nodes = scattering(POTENTIALS[name](strength, scale))

    ok = abs(r0 - r0_target) < 1e-6
    if nodes_target is not None and nodes != nodes_target:
        # Two solutions can share (a, r0) and differ in node count, and then
        # they are not the same physical state. This check is not optional.
        ok = False

    return strength, scale, a, r0, nodes, ok


# =============================================================================
#  Closed forms, used to check the solver
# =============================================================================
def a_well(v, R=1.0):
    """a/R = 1 - tan(x)/x with x = sqrt(2v).  Eq. (80) of [1]."""
    x = math.sqrt(2.0 * v)
    return R * (1.0 - math.tan(x) / x)


def r0_well(v, R=1.0):
    """r0/R for the square well.  Eq. (92) of [1].

    At every pole of a, that is x = pi/2 + n pi, this returns exactly 1: the
    effective range equals the range of the potential. Sharpest single check
    we have on the whole solver.
    """
    a = a_well(v, R)
    k = math.sqrt(2.0 * v) / R
    line_part = R - R ** 2 / a + R ** 3 / (3.0 * a ** 2)
    wave_part = (1.0 - R / a) ** 2 / math.sin(k * R) ** 2 \
        * (R / 2.0 - math.sin(2.0 * k * R) / (4.0 * k))
    return 2.0 * (line_part - wave_part)


def E_zero_range(a, h2_2mu):
    """E = -(hbar^2/2mu) / a^2. Throws away the size of the potential."""
    return -h2_2mu / a ** 2


def E_finite_range(a, r0, h2_2mu):
    """Solves kappa = 1/a + r0 kappa^2 / 2. Keeps r0, drops only O(k^4)."""
    kappa = (1.0 - math.sqrt(1.0 - 2.0 * r0 / a)) / r0
    return -h2_2mu * kappa ** 2


# =============================================================================
#  Tabulated data
# =============================================================================
# Table 2 of [1]: the three targets. "nodes" is the bound-state count required.
TARGETS = {
    "nn":        {"a": -18.5,     "r0": 2.7, "nodes": 0},   # neutron-neutron
    "unitarity": {"a": math.inf,  "r0": 1.0, "nodes": 0},   # |a| -> infinity
    "deuteron":  {"a": 5.4,       "r0": 1.7, "nodes": 1},
}

# Tables 3 and 4 of [1]: published parameters. Used as the starting guess and
# as the thing we have to reproduce. "r0_pub" is what those parameters actually
# give, which is not always the Table 2 target -- see SOURCES["D1"].
PUBLISHED = {
    ("nn", "well"):         {"p1": 1.1096,     "p2": 0.3918,     "r0_pub": 2.70},
    ("nn", "mpt"):          {"p1": 0.9071,     "p2": 0.7991,     "r0_pub": 2.70},
    ("nn", "gauss"):        {"p1": 1.2121,     "p2": 0.5672,     "r0_pub": 2.70},
    ("nn", "lj"):           {"p1": 9.86668911, "p2": 3.08836698, "r0_pub": 2.71},
    ("unitarity", "well"):  {"p1": 1.2337,     "p2": 1.0000,     "r0_pub": 1.00},
    ("unitarity", "mpt"):   {"p1": 1.0000,     "p2": 2.0000,     "r0_pub": 1.00},
    ("unitarity", "gauss"): {"p1": 1.3420,     "p2": 1.4349,     "r0_pub": 1.00},
    ("unitarity", "lj"):    {"p1": 0.26462461, "p2": 0.00034068, "r0_pub": 1.00},
    ("deuteron", "well"):   {"p1": 1.7575,     "p2": 0.5000,     "r0_pub": 1.70},
    ("deuteron", "mpt"):    {"p1": 1.4388,     "p2": 0.8631,     "r0_pub": 1.73},
    ("deuteron", "gauss"):  {"p1": 1.9102,     "p2": 0.6754,     "r0_pub": 1.70},
    ("deuteron", "lj"):     {"p1": 6.81472,    "p2": 0.90485319, "r0_pub": 1.70},
}

# Table 1 of [1]: two real systems. E is measured and is never an input.
# hbar^2/2mu is tabulated nowhere, so we get it by inverting the zero-range
# formula on the published pair (a, E_zr). Recovering the known constants
# 41.47 MeV.fm^2 and 12.12 K.A^2 is our check that this table is consistent.
SYSTEMS = {
    "deuteron":  {"a": 5.4112, "r0": 1.7436, "E": -2.224,
                  "unit": "MeV", "E_zr": -1.416, "source": "[2]"},
    "he4_dimer": {"a": 90.4,   "r0": 8.0,    "E": -1.62e-3,
                  "unit": "K",  "E_zr": -1.48e-3, "source": "[3]"},
}
for system in SYSTEMS.values():
    system["h2_2mu"] = -system["E_zr"] * system["a"] ** 2

# Starting points: the deuteron solution. The problem is scale invariant --
# multiplying all lengths by L takes (a, r0) to (L a, L r0) -- so for another
# system we divide the scale parameter by L, and C12 by L^6 because it
# multiplies r^-12.
GUESS = {"well": (1.7806, 0.48415), "mpt": (1.44397, 0.853766),
         "gauss": (1.93585, 0.652551), "lj": (7.8433, 1.27426)}

SOURCES = {
    "[1]": "Macedo-Lima & Madeira, Rev. Bras. Ensino Fis. 45, e20230079 (2023), "
           "doi:10.1590/1806-9126-RBEF-2023-0079",
    "[2]": "Hackenburg, Phys. Rev. C 73, 044002 (2006), "
           "doi:10.1103/PhysRevC.73.044002",
    "[3]": "Cencek et al., J. Chem. Phys. 136, 224303 (2012), "
           "doi:10.1063/1.4712218",
    "D1":  "For (nn, lj) and (deuteron, mpt) the parameters published in [1] "
           "give r0 = 2.71 and 1.73, while Table 2 of the same paper lists the "
           "targets as 2.70 and 1.70.",
}
