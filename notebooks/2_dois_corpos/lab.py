"""Two-body low-energy scattering: potentials, solver, tabulated data.

Reasoning, derivations and the measurements behind every constant are in
notas_teoria/Theory_and_Implementation.pdf. This file is the implementation.

    UNITS         The solver fixes no length unit. It carries whatever unit the
                  potential's parameters use: fm for the model potentials,
                  angstrom for the helium ones. Energies and potentials are
                  divided by hbar^2/m_r, so [V] = [E] = (that length)^-2.
                  Back to physics: V_phys = 2 (hbar^2/2m_r) V_code.

    CONSTANTS     EPS_TAIL, EPS_CORE   where the potential is cut, dimensionless
                  UNITARITY            when |r0/a| counts as zero
                  R_START_FRACTION     where integration begins, fraction of R
                  POINTS               grid size

    POTENTIALS    Well, PoschlTeller, Gaussian, LennardJones -- two parameters
                  each, tunable. HFDFamily (HFDHE2, HFD-B, LM2M2) and
                  TangToenniesYiu -- measured helium potentials, no free
                  parameters, so their (a, r0) are predictions.

    SOLVER        numerov      integrates u'' = 2(V-E)u on a grid uniform in ln r
                  scattering   reads a, r0 and the node count at E = 0
                  bound_energy finds E < 0 by matching onto exp(-kappa r)
                  tune         inverse problem: given (a, r0), find parameters

    CLOSED FORMS  a_well, r0_well, E_zero_range, E_finite_range -- used to check
                  the solver, never inside it.

    DATA          TARGETS, PUBLISHED, SYSTEMS, GUESS, SOURCES, all with citation.

Sources: [1] Macedo-Lima & Madeira, Rev. Bras. Ensino Fis. 45, e20230079 (2023).
[2] Hackenburg, Phys. Rev. C 73, 044002 (2006). [3] Cencek et al., J. Chem.
Phys. 136, 224303 (2012). [4] Aziz et al. 1979, 1987, 1991; Tang, Toennies &
Yiu 1995. [5] Motovilov et al., Eur. Phys. J. D 13, 33 (2001), Tables IX, X
and I -- potential parameters come from there, never from memory.
"""
import math
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import brentq
from scipy.special import gammainc

# Every value below was fixed by measurement, not chosen. The convergence
# study and the numbers behind each one are in Secs. 2 and 19 of the notes.

EPS_TAIL = 1e-6        # cut the tail at |V(R)| R^2. Dimensionless: [V] = L^-2
EPS_CORE = 1e4         # cut the core at |V(r_min)| r_min^2, same reason
UNITARITY = 1e-3       # below this |r0/a|, a counts as infinite
R_START_FRACTION = 1e-6  # where integration begins, as a fraction of R
BREAKDOWN_MARGIN = 1.15  # TTY only: how far outside its validity limit to start
H2_2MU_HE4 = 12.11932  # hbar^2/2m_r for the He-4 dimer, K.A^2, CODATA.
                       # SYSTEMS holds an inferred value 0.20% lower; that is a
                       # check on Table 1 and must not be used for a potential
                       # with no free parameters, where 0.20% becomes 6% in a.


def raio_de_corte(V, r_dentro, fator=1e5):
    """Smallest R beyond r_dentro with |V(R)| R^2 = EPS_TAIL.

    The cut is dimensionless because [V] = (length)^-2, so |V| r^2 is a pure
    number whatever unit the potential was written in. Four potentials need
    this and none of them can solve it on paper.
    """
    return brentq(lambda r: abs(V(r)) * r * r - EPS_TAIL,
                  r_dentro, r_dentro * fator, xtol=1e-13)


# =============================================================================
#  The four potentials
#
#  Each class holds two things: the formula for V(r), and where it stops (R).
#  R comes from the dimensionless cut |V(R)| R^2 = EPS_TAIL. Closed form for
#  the well and the Lennard-Jones; brentq for the other two.
# =============================================================================
class Well:
    """V = -v mu^2 inside r < R = 1/mu, zero outside.  Eq. (70) of [1].

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
    """V = -v mu^2 / cosh^2(mu r).  Eq. (116) of [1].  Smooth.

    [1] gives a closed form for the scattering length. It does NOT give one for
    r0 or for the zero-energy wavefunction except at unitarity, so this is not a
    potential we can check the solver against the way we check the well.
    """
    name = "Poschl-Teller"

    def __init__(self, v, mu):
        self.v = v
        self.mu = mu
        self.r_min = 0.0
        # v mu^2 R^2 / cosh^2(mu R) = EPS_TAIL. No closed form once the r^2 is
        # there, so brentq. Paid once, at construction.
        # cosh overflows past ~350/mu, so the bracket stops short of it
        self.R = raio_de_corte(self.V, 1.0 / mu, fator=300.0)

    def V(self, r):
        return -self.v * self.mu ** 2 / np.cosh(self.mu * r) ** 2


class Gaussian:
    """V = -v mu^2 exp(-mu^2 r^2).  Eq. (120) of [1].  Tail dies very fast."""
    name = "Gaussian"

    def __init__(self, v, mu):
        self.v = v
        self.mu = mu
        self.r_min = 0.0
        self.R = raio_de_corte(self.V, 1.0 / mu, fator=1e2)

    def V(self, r):
        return -self.v * self.mu ** 2 * np.exp(-(self.mu * r) ** 2)


class LennardJones:
    """V = (1/2)(C12/r^12 - C6/r^6). The constructor takes C6 FIRST.

    The factor 1/2 is not in Eq. (121) of [1] and is needed: without it the
    published deuteron parameters give a = 1.435 instead of 5.405.
    """
    name = "Lennard-Jones"

    def __init__(self, C6, C12):
        self.C6 = C6
        self.C12 = C12
        # Both cuts keep a closed form after the r^2, with the exponent shifted
        # by two: 0.5 C12/r^10 = EPS_CORE and 0.5 C6/R^4 = EPS_TAIL.
        self.r_min = (0.5 * C12 / EPS_CORE) ** (1.0 / 10.0)
        self.R = (0.5 * C6 / EPS_TAIL) ** (1.0 / 4.0)

    def V(self, r):
        return 0.5 * (self.C12 / r ** 12 - self.C6 / r ** 6)


class HFDFamily:
    """Aziz semi-empirical helium potentials: one form, three parameter sets.

    V = eps [A exp(-a x + b x^2) - F(x)(C6/x^6 + C8/x^8 + C10/x^10)] + Va,
    with x = r/rm and F(x) = exp(-(D/x - 1)^2) below x = D, 1 above.
    Va is the LM2M2 bump and is zero for the other two.
    Parameters: [5], Table IX. Soft core, so r_min = 0.
    """
    name = "HFD family"
    AA = 0.0                    # Va is absent unless a subclass says otherwise
    Z1 = Z2 = 0.0

    def __init__(self, h2_2mu):
        self.h2_2mu = h2_2mu
        self.r_min = 0.0                    # soft core: see the docstring

        self.R = raio_de_corte(self.V, self.RM)

    def V(self, r):
        x = np.asarray(r, dtype=float) / self.RM
        F = np.where(x <= self.D, np.exp(-(self.D / x - 1.0) ** 2), 1.0)
        disp = self.C6 / x ** 6 + self.C8 / x ** 8 + self.C10 / x ** 10
        Vb = self.A * np.exp(-self.ALPHA * x + self.BETA * x * x) - F * disp
        if self.AA != 0.0:
            inside = (x >= self.Z1) & (x <= self.Z2)
            arg = 2.0 * math.pi * (x - self.Z1) / (self.Z2 - self.Z1) - math.pi / 2.0
            Vb = Vb + np.where(inside, self.AA * (np.sin(arg) + 1.0), 0.0)
        return self.EPS * Vb / (2.0 * self.h2_2mu)


class AzizHFDHE2(HFDFamily):
    """HFDHE2, Aziz et al., J. Chem. Phys. 70, 4330 (1979). The 1979 original."""
    name = "Aziz HFDHE2"
    EPS, RM, D = 10.8, 2.9673, 1.241314
    A, ALPHA, BETA = 544850.4, 13.353384, 0.0
    C6, C8, C10 = 1.3732412, 0.4253785, 0.178100


class AzizHFDB(HFDFamily):
    """HFD-B, Aziz, McCourt & Wong, Molec. Phys. 61, 1487 (1987)."""
    name = "Aziz HFD-B"
    EPS, RM, D = 10.948, 2.963, 1.4826
    A, ALPHA, BETA = 184431.01, 10.43329537, -2.27965105
    C6, C8, C10 = 1.36745214, 0.42123807, 0.17473318


class AzizLM2M2(HFDFamily):
    """LM2M2, Aziz & Slaman, J. Chem. Phys. 94, 8047 (1991).

    HFD-B plus the Va bump. That small addition moves the dimer energy from
    -1.685 mK to -1.303 mK, a 23% change from a term whose peak is 0.0052 in
    units of eps -- another reading of the halo amplification.
    """
    name = "Aziz LM2M2"
    EPS, RM, D = 10.97, 2.9695, 1.4088
    A, ALPHA, BETA = 189635.353, 10.70203539, -1.90740649
    C6, C8, C10 = 1.34687065, 0.41308398, 0.17060159
    AA, Z1, Z2 = 0.0026, 1.003535949, 1.454790369


class TangToenniesYiu:
    """TTY helium potential, purely theoretical. Parameters: [5], Table X.

    V = A [D x^p exp(-2 b x) - sum_n C_2n f_2n(x) x^-2n], x in bohr, V in K.
    p = 7/(2b) - 1 is fixed by atomic asymptotics, not fitted. C_2n follow the
    published recurrence. f_2n is gammainc(2n+1, bx): written out as
    1 - exp(-bx) sum(...) it cancels catastrophically and V oscillates in sign
    between 0.36 and 0.44 A.

    The form breaks below x = p/(2b), where b(x) changes sign and V diverges to
    -infinity. Integration starts BREAKDOWN_MARGIN times outside that radius.
    """
    name = "TTY"

    A_K, BETA, D = 315766.2067, 1.3443, 7.449
    C6, C8, C10 = 1.461, 14.11, 183.5
    N = 12
    BOHR = 0.52917                      # angstrom per bohr, the value used in [5]
    RM = 2.97                           # only a bracket start, not a parameter

    def __init__(self, h2_2mu):
        self.h2_2mu = h2_2mu
        # build the dispersion coefficients once, by the published recurrence
        c = [self.C6, self.C8, self.C10]
        while len(c) < self.N - 2:
            c.append((c[-1] / c[-2]) ** 3 * c[-3])
        self.C = c
        self.p = 7.0 / (2.0 * self.BETA) - 1.0

        # r_min is NOT set by EPS_CORE here: the core is soft, so |V| r^2 goes
        # to zero at both ends and never reaches it. What sets
        # r_min is the validity of the published form: b(x) changes sign at
        # x = p/(2 beta), and the expression is meaningless below it. We start
        # a margin outside that radius. See BREAKDOWN_MARGIN.
        self.r_breakdown = self.p / (2.0 * self.BETA) * self.BOHR
        self.r_min = BREAKDOWN_MARGIN * self.r_breakdown
        self.R = raio_de_corte(self.V, self.RM)

    def V(self, r):
        x = np.atleast_1d(np.asarray(r, dtype=float)) / self.BOHR   # angstrom -> bohr
        bx = 2.0 * self.BETA * x - self.p                   # b(x) * x, algebraically
        Vex = self.D * x ** self.p * np.exp(-2.0 * self.BETA * x)

        # The damping function IS the regularised lower incomplete gamma:
        #
        #     f_2n(x) = 1 - exp(-bx) sum_{k=0}^{2n} (bx)^k / k! = P(2n+1, bx)
        #
        # because Q(m, y) = exp(-y) sum_{k=0}^{m-1} y^k/k! for integer m. Using
        # gammainc instead of the written-out sum is not a convenience: the
        # literal form subtracts two nearly equal numbers, and C_2n/x^(2n)
        # reaches 1e12 by n = 12, so the surviving roundoff is amplified twelve
        # orders. Measured with the literal form, V oscillates in sign between
        # 0.36 and 0.44 A -- +3.5e7, -9.4e6, +3.0e6 K -- which is cancellation,
        # not physics. gammainc is stable over the whole range and handles the
        # large-bx limit, where f_2n -> 1, without overflowing.
        Vdisp = np.zeros_like(x)
        for j, C2n in enumerate(self.C):
            n = j + 3
            Vdisp = Vdisp - C2n * gammainc(2 * n + 1, bx) / x ** (2 * n)

        out = self.A_K * (Vex + Vdisp) / (2.0 * self.h2_2mu)
        return out if np.ndim(r) else out[0]


# The strength parameter always comes first, the scale parameter second.
# AzizHFDB is deliberately NOT in here: it takes no (strength, scale) pair, so
# it cannot be tuned, and tune() must never be handed it.
POTENTIALS = {"well": Well, "mpt": PoschlTeller,
              "gauss": Gaussian, "lj": LennardJones}


# One grid size for all four potentials -- see numerov() for why a logarithmic
# grid lets a single number work here. 8001 is not arbitrary: total error is
# truncation (falls with more points) plus roundoff (grows with them), and this
# sits at the minimum. The three purely attractive potentials are at machine
# precision, 1e-11, and are insensitive. The Lennard-Jones sets the value:
# past ~8000 points its r0 degrades again, by 3e-5 at 32001 and 3e-4 at 64001.
POINTS = 8001


# =============================================================================
#  Numerov on a logarithmic grid
# =============================================================================
def numerov(pot, E):
    """Integrate u'' = 2(V-E)u from the start radius out to R.

    Grid uniform in x = ln r, where u = e^(x/2) w turns the equation into
    w'' = [2(V-E)e^(2x) + 1/4] w. The last point is assigned, not computed:
    exp(log(R)) can overshoot R and the well then reads V = 0 there.
    Slope at the edge from a five-point one-sided stencil, interior points only.
    Returns r, u, x, h, slope, r_start.
    """
    # Where to start: outside the hard core if there is one, otherwise close
    # enough to the origin that the piece we skip is analytic (see scattering).
    r_start = pot.r_min if pot.r_min > 0.0 else pot.R * R_START_FRACTION

    x = np.linspace(math.log(r_start), math.log(pot.R), POINTS)
    h = x[1] - x[0]
    r = np.exp(x)

    # PIN THE LAST POINT. exp(log(R)) is not guaranteed to return R, and when
    # it overshoots by one ulp the square well evaluates V = 0 there instead of
    # the depth, because its V tests r <= R. Measured: that single point takes
    # the error in a from 1.6e-10 to 1.1e-3, seven orders, and it happens for
    # 14% of mu values -- one case in seven, silently. The published parameters
    # avoid it by luck, which is why it never showed until a rescaling test hit
    # mu = 1/3.
    r[-1] = pot.R
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
    """a, r0 and the number of bound states, all read off at E = 0.

    Beyond R the potential is gone and u is the straight line C(1 - r/a), so
    all three come from the value and slope at the edge. Derivations are in
    the notes, Secs. 9 and 10; the measurements behind each choice are there
    too and are not repeated here.
    """
    r, u, x, h, slope, r_start = numerov(pot, 0.0)

    # 1/a, not a: it passes through zero at a resonance while a runs to
    # -infinity and returns from +infinity. No branch, no sign to invent.
    inv_a = slope / (slope * pot.R - u[-1])
    a = 1.0 / inv_a

    # rescale u onto the line, then r0 = 2 int (line^2 - u^2) dr, Eq. (56).
    # The factor r is the Jacobian: on this grid dr = r dx.
    u = u * (1.0 - pot.R / a) / u[-1]
    line = 1.0 - r / a
    r0 = 2.0 * simpson((line ** 2 - u ** 2) * r, x=x)

    # the stretch below r_start, where u is zero or proportional to r, added
    # analytically rather than integrated
    r0 += 2.0 * (r_start - r_start ** 2 / a + r_start ** 3 / (3.0 * a ** 2))

    # zeros of u count the bound states. The line contributes one more at
    # r = a, which stops counting as unitarity is approached -- and "approach"
    # has to be dimensionless, hence r0/a.
    nodes = int(np.count_nonzero(u[1:-1] * u[2:] < 0.0))
    if pot.R < a and abs(r0 / a) > UNITARITY:
        nodes += 1

    return a, r0, nodes

# =============================================================================
#  Negative energy: the bound state
# =============================================================================
def bound_energy(pot, guess):
    """Ground-state energy in code units, bracketed around `guess`.

    Matches onto the decaying exponential: u'(R) + kappa u(R) = 0 with
    kappa = sqrt(-2E), divided by max(|u|, |u'|, 1) so the residual is order
    one for every potential. Raw it runs from 1e-2 to 1e21.
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
    """Parameters of `name` that reproduce (a_target, r0_target).

    Outer loop moves the scale until r0 matches, inner moves the strength
    until 1/a does. Chases 1/a rather than a, because a has poles.
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
    """E = -(hbar^2/2mu) / a^2. Keeps a and nothing else: r0 is dropped."""
    return -h2_2mu / a ** 2


def E_finite_range(a, r0, h2_2mu):
    """Solves kappa = 1/a + r0 kappa^2 / 2. Keeps r0, drops only O(k^4)."""
    kappa = (1.0 - math.sqrt(1.0 - 2.0 * r0 / a)) / r0
    return -h2_2mu * kappa ** 2


# =============================================================================
#  Tabulated data
# =============================================================================
# Table 2 of [1]: the three targets. The bound-state count is NOT a column of
# that table -- it comes from the text of Sec. 4.5, which asks for a nodeless
# u(r) when a < 0 and one node when a > 0.
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
# multiplying all lengths by L takes (a, r0) to (L a, L r0). Since [V_code] is
# fm^-2, the rescaled potential must satisfy V'(L r) = V(r) / L^2. For the three
# one-scale potentials that means dividing the scale parameter by L. For the
# Lennard-Jones it means C6 -> L^4 C6 and C12 -> L^10 C12, which follows from
# [C6] = fm^4 and [C12] = fm^10.
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
