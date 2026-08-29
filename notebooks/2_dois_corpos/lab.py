"""Two-body low-energy scattering: potentials, solver and tabulated data.

The notebook next to this file needs nothing else.

UNITS
    Distances are in fm and are NOT rescaled. Energies and potentials are,
    by m_r / hbar^2:

        V_code = (m_r / hbar^2) V_phys,   E_code = (m_r / hbar^2) E_phys

    so the radial equation becomes  u'' = 2 (V - E) u.

    Note what this does and does not say. It is NOT "hbar = m_r = 1 and
    everything is dimensionless": with r still in fm,

        [V_code] = [E_code] = fm^-2.

    That single line fixes the units of the Lennard-Jones coefficients,
    [C6] = fm^4 and [C12] = fm^10, and with them the exponents in the scaling
    rule under GUESS. Going back:

        E_phys = 2 (hbar^2 / 2 m_r) E_code,   V_phys = 2 (hbar^2 / 2 m_r) V_code

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
        The four potentials (Eqs. 70, 116, 120, 121), the targets (Table 2), the
        published parameters (Tables 3 and 4), and the square-well closed
        forms (Eqs. 80 and 92).
    [2] Hackenburg, R. W. Phys. Rev. C 73, 044002 (2006)
        doi:10.1103/PhysRevC.73.044002
        Deuteron a, r0 and binding energy. Ref. [23] of [1].
    [3] Cencek, W. et al. J. Chem. Phys. 136, 224303 (2012)
        doi:10.1063/1.4712218
        Helium-4 dimer a, r0 and binding energy. Ref. [22] of [1].
    [4] The helium-helium potentials, original references:
        HFDHE2  Aziz, R. A. et al. J. Chem. Phys. 70, 4330 (1979)
        HFD-B   Aziz, R. A.; McCourt, F. R. W.; Wong, C. C. K.
                Molec. Phys. 61, 1487 (1987)
        LM2M2   Aziz, R. A.; Slaman, M. J. J. Chem. Phys. 94, 8047 (1991)
        TTY     Tang, K. T.; Toennies, J. P.; Yiu, C. L.
                Phys. Rev. Lett. 74, 1546 (1995)
    [5] Motovilov, A. K.; Sandhas, W.; Sofianos, S. A.; Kolganova, E. A.
        Eur. Phys. J. D 13, 33 (2001), arXiv:physics/9910016
        Appendix, Tables IX and X: the parameters of all four potentials in [4],
        tabulated side by side. Table I: the dimer energies and scattering
        lengths they obtain, which is what test_helium_benchmark checks against.
        PARAMETERS ARE TAKEN FROM HERE, not from memory or from secondary
        summaries. A mistyped potential parameter fails nowhere -- it just
        produces a plausible wrong answer.
"""
import math
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import brentq
from scipy.special import gammainc

# Below V_ZERO we treat the potential as gone. Above V_CORE the Lennard-Jones
# core is too steep to integrate through, and the wavefunction is zero there
# anyway, so we start outside it.
#
# [1] prescribes V(r_min) ~ 1e10 for this cut. We use 1e5 because it was
# measured to be already converged: r0 changes by 2e-6 between 1e4 and 1e6,
# while 1e3 is wrong by 0.2%. Starting deeper only costs stiffness.
V_ZERO = 1e-12
V_CORE = 1e5

# How close to unitarity counts as unitarity, for the bound-state count only.
#
# The question "10 to which power is infinity?" has no answer while it is asked
# about `a` itself, because `a` carries a length and the answer would then
# depend on the system. Asked about a dimensionless ratio it does have one.
# The ratio to use is r0/a: it is invariant under rescaling, and it vanishes at
# unitarity by definition. R/a would not do -- R is a truncation we chose, not
# physics, and for the Lennard-Jones at unitarity R/a = 1.5e-3, which is not
# small at all.
#
# Measured separation, over the twelve published cases:
#     |r0/a| at unitarity      2.2e-11 ... 2.1e-05
#     |r0/a| where a node counts          3.14e-01
# a gap of about 15000. The threshold below sits 47x above the first group and
# 314x below the second.
UNITARITY = 1e-3

# Where integration starts when there is no hard core, as a fraction of R.
# Below this radius u is proportional to r and the piece is added analytically
# in scattering(), so the value is a truncation knob and belongs in the error
# budget with V_ZERO and V_CORE rather than buried in numerov().
R_START_FRACTION = 1e-6

# hbar^2 / 2 mu for the He-4 dimer, from CODATA, in K.A^2.
#
# SYSTEMS below carries an INFERRED value, 12.0948, obtained by inverting the
# zero-range formula on the published pair. That inference is a consistency
# check on Table 1 and it is 0.20% below this one. For the four tunable
# potentials the difference is invisible, because they are fitted to a target
# and absorb it. For a potential with no free parameters it is not:
#
#     d(ln a) / d(ln h2_2mu) = 48.7
#
# near unitarity, so 0.20% in the constant becomes 6% in a. The He-4 dimer is a
# halo state and amplifies everything. Parameter-free potentials use this value.
H2_2MU_HE4 = 12.11932

# How far outside its breakdown radius the TTY potential is allowed to start,
# as a multiple of that radius. Its damping argument b(x) changes sign at
# x = p/(2 beta) = 0.3156 A and the published expression means nothing below.
# This is a truncation knob like V_ZERO and V_CORE and belongs in the error
# budget: the check is that (a, r0) do not move when it is varied, which they
# do not, because V there is 2e5 K and the wavefunction is zero to machine
# precision long before the grid reaches it.
BREAKDOWN_MARGIN = 1.15


# =============================================================================
#  The four potentials
#
#  Each class holds two things: the formula for V(r), and where it stops (R).
#  R is written in closed form -- we solve |V(R)| = V_ZERO by hand for each
#  shape, so no search is needed.
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
        # v mu^2 / cosh^2(mu R) = V_ZERO
        self.R = math.acosh(math.sqrt(v * mu ** 2 / V_ZERO)) / mu

    def V(self, r):
        return -self.v * self.mu ** 2 / np.cosh(self.mu * r) ** 2


class Gaussian:
    """V = -v mu^2 exp(-mu^2 r^2).  Eq. (120) of [1].  Tail dies very fast."""
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
    """V = (1/2) (C12/r^12 - C6/r^6).  Hard core plus van der Waals tail.

    THE FACTOR 1/2 IS NOT IN THE PAPER, AND IT IS NEEDED.

    Eq. (121) of [1] reads V = (hbar^2 / m_r)(C12/r^12 - C6/r^6). In the units
    used here that is C12/r^12 - C6/r^6, with no 1/2. Taken literally it does
    not reproduce Table 4 of the same paper -- measured, for the deuteron row:

        published C12 = 0.90485319, C6 = 6.81472, and Table 4 reports
        a = 5.4 fm, r0 = 1.70 fm

        with the 1/2:  a = 5.405,  r0 = 1.699   <- reproduces Table 4
        without it:    a = 1.435,  r0 = 1.412   <- does not

    Eqs. (70) and (120), for the well and the Gaussian, carry the same
    hbar^2/m_r prefactor and DO reproduce Table 3 as printed. So the mismatch
    is specific to Eq. (121); the likely reading is 2 m_r in its denominator.
    """
    name = "Lennard-Jones"

    def __init__(self, C6, C12):
        self.C6 = C6
        self.C12 = C12
        self.r_min = (0.5 * C12 / V_CORE) ** (1.0 / 12.0)   # core: 0.5 C12/r^12 = V_CORE
        self.R = (0.5 * C6 / V_ZERO) ** (1.0 / 6.0)         # tail: 0.5 C6/r^6  = V_ZERO

    def V(self, r):
        return 0.5 * (self.C12 / r ** 12 - self.C6 / r ** 6)


class HFDFamily:
    """The Aziz semi-empirical helium potentials. HFDHE2, HFD-B, LM2M2.

    THESE HAVE NO FREE PARAMETERS. The four potentials above are shapes we tune
    until they sit at a chosen (a, r0). These are fixed by spectroscopy and
    ab initio theory, so feeding one in and reading (a, r0) out is a prediction,
    not a fit -- the only falsifiable objects in this file.

    All three share one functional form, differing only in the parameter table:

        V(r) = eps [ Vb(x) + Va(x) ],                    x = r / rm
        Vb(x) = A exp(-alpha x + beta x^2)
                - F(x) (C6/x^6 + C8/x^8 + C10/x^10)
        F(x)  = exp(-(D/x - 1)^2) for x <= D, else 1

    Va is zero for HFDHE2 and HFD-B. LM2M2 adds a bump over a finite interval:

        Va(x) = Aa [ sin(2 pi (x - z1)/(z2 - z1) - pi/2) + 1 ],  z1 <= x <= z2

    HFDHE2 has beta = 0, which is what "B-type" refers to: HFD-B added that
    quadratic term to the exponent.

    UNITS ARE THE POINT HERE
        The published parameters are in kelvin and angstrom, not in code units.
        By the definitions at the top of this file,

            V_code = V_phys / (2 * h2_2mu),   h2_2mu = 12.09 K.A^2 for He2

        so lengths come out in angstrom and V_code in A^-2. This is the first
        place in the laboratory where the unit convention has to be used rather
        than assumed, and it is a check on it: the potential must reproduce its
        own tabulated minimum, V(rm) = -eps, which it does to every digit.

    THE CORE IS SOFT, AND THAT CHANGES THE START
        The Lennard-Jones diverges as 1/r^12, so the code steps around it by
        starting at r_min where V = V_CORE. Aziz does not diverge: as r -> 0 the
        exponential saturates and F(x) kills the dispersion terms, so V climbs
        to a CEILING of eps*A/(2 h2_2mu) = 8.35e4 in code units -- below V_CORE.
        There is no radius where V = V_CORE, so r_min = 0 and the regular-origin
        start applies, the same one the well and the Gaussian use. Trying to
        solve V = V_CORE here raises "f(a) and f(b) must have different signs",
        which is the bracket telling the truth.

    R IS NOT CLOSED FORM
        Unlike the other four, |V(R)| = V_ZERO cannot be solved on paper here,
        so R comes from brentq. Honest cost of using a real potential.

    Parameters from Table IX of [5], which tabulates all three side by side.
    """
    name = "HFD family"
    AA = 0.0                    # Va is absent unless a subclass says otherwise
    Z1 = Z2 = 0.0

    def __init__(self, h2_2mu):
        self.h2_2mu = h2_2mu
        self.r_min = 0.0                    # soft core: see the docstring
        # the tail: |V| = V_ZERO, bracketed beyond the minimum
        self.R = brentq(lambda r: abs(self.V(r)) - V_ZERO, self.RM, 1e4, xtol=1e-12)

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
    """TTY, Tang, Toennies & Yiu, Phys. Rev. Lett. 74, 1546 (1995).

    Purely theoretical, from perturbation theory -- no fit to experiment at all.
    Structurally different from the Aziz family, and it lands within 0.2% of
    LM2M2 on both the dimer energy and the scattering length, which is a
    non-trivial agreement between a semi-empirical fit and a first-principles
    calculation.

        V(x) = A [ Vex(x) + Vdisp(x) ],        x in bohr, V in kelvin
        Vex(x)   = D x^p exp(-2 beta x),       p = 7/(2 beta) - 1
        Vdisp(x) = -sum_{n=3}^{N} C_2n f_2n(x) / x^(2n)
        f_2n(x)  = 1 - exp(-b x) sum_{k=0}^{2n} (b x)^k / k!
        b(x)     = 2 beta - (7/(2 beta) - 1) / x

    Only C6, C8 and C10 are given; the rest come from the recurrence

        C_2n = (C_{2n-2} / C_{2n-4})^3 C_{2n-6}

    up to n = N = 12, so the series runs C6 through C24.

    Parameters from Table X of [5]. Note the length unit is the bohr, not the
    angstrom: this is the one potential here whose native length is atomic, so
    r has to be converted on the way in.

    THE FORM BREAKS DOWN AT SMALL r, AND THAT IS NOT A BUG HERE
        The damping argument b(x) = 2 beta - p/x changes sign at
        x = p/(2 beta) = 0.596 bohr = 0.3156 A. Below that, b x is negative, the
        truncated series stops damping anything, and Vdisp runs off to -infinity
        -- V(0.2 A) evaluates to -7e13 K. The published form is simply not meant
        to be evaluated there. Ref. [5] never sees it: they impose a hard core at
        c = 1.0 A.

        We do the same thing the Lennard-Jones does: start at r_min where
        V = V_CORE, which lands at ~0.44 A, well outside the breakdown. The
        wavefunction is zero to machine precision through a wall that high, so
        nothing physical is lost -- but the potential must never be sampled
        below r_min, and that is why r_min is found by bracketing DOWN from the
        minimum rather than up from zero.
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

        # r_min is NOT set by V_CORE here. The core is soft -- V saturates at
        # about 2.0e5 K, which is 8.5e3 in code units, below V_CORE. What sets
        # r_min is the validity of the published form: b(x) changes sign at
        # x = p/(2 beta), and the expression is meaningless below it. We start
        # a margin outside that radius. See BREAKDOWN_MARGIN.
        self.r_breakdown = self.p / (2.0 * self.BETA) * self.BOHR
        self.r_min = BREAKDOWN_MARGIN * self.r_breakdown
        self.R = brentq(lambda s: abs(self.V(s)) - V_ZERO, self.RM, 1e4, xtol=1e-12)

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
    r_start = pot.r_min if pot.r_min > 0.0 else pot.R * R_START_FRACTION

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
    #
    # The branch below is mathematically right and numerically unreachable: a
    # float slope is never exactly 0.0. Swept across the well's resonance it
    # goes 7.9e-3, 7.9e-4, 3.5e-7, -7.8e-4, -7.8e-3 -- small, then small and
    # negative, never zero. So `a` is never inf in practice; it is large and
    # it CHANGES SIGN through the resonance, which is the physical behaviour.
    # Returning +inf here would be picking one side of that sign flip by fiat.
    # Kept only so the function is total.
    if slope == 0.0:
        a = math.inf
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

    # Count the zeros of the E = 0 solution: they give the number of bound
    # states. [1] uses this as a mandatory check -- "checking the number of
    # nodes of the radial function is necessary to guarantee that it is indeed
    # the situation we wanted to reproduce" -- expecting a nodeless u(r) for
    # a < 0 and one node for a > 0.
    nodes = 0
    for i in range(1, POINTS - 1):
        if u[i] * u[i + 1] < 0.0:
            nodes = nodes + 1
    if pot.R < a and abs(r0 / a) > UNITARITY:
        # The line crosses zero at r = a, outside the potential. That zero is a
        # node when a is finite, and marches off to infinity as we approach
        # unitarity, where by convention the new state does not count yet.
        # "Approach unitarity" has to be a dimensionless statement -- see
        # UNITARITY above for why r0/a and not a itself.
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

    WHY THE SUM AND NOT slope/u + kappa
        Measured: with the ratio, brentq fails to bracket for three of the four
        potentials -- it reports the same sign at both ends of the interval.
        The sum brackets for all four.

        An earlier version of this docstring said the reason was that u(R) can
        pass through zero, "the deuteron has a node". That is wrong twice. The
        deuteron's bound state is the ground state and is nodeless; what carries
        a node is the zero-energy scattering solution of a potential that
        supports a bound state, which is the function scattering() looks at, not
        this one. And u(R) never comes near zero here: swept across the bracket,
        its smallest value over the four potentials is 7e-3.

    WHY THE DENOMINATOR
        Normalisation, nothing else. It moves the root by 3e-16. What it fixes
        is scale: the raw residual is about 1e-2 for the square well and 1e21
        for the Lennard-Jones, whose u climbs out of the hard core before
        reaching the edge. One absolute tolerance cannot serve both, so we
        divide by the largest thing in play and every potential lands at order
        one.

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
    ill-conditioned exactly where it matters.

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
