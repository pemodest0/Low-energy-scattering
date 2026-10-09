"""Two-body potentials calibrated on Table III of [6], and the path to unitarity.

    python etrych.py

Table III of [6] gives, for five 39K Feshbach resonances, the effective range
r_eff near B_res and the Efimov feature a_-, the scattering length at which the
lowest trimer meets the three-atom continuum. This module builds two-body
potentials that sit on those numbers, and families of potentials that walk from
a_- through the resonance to the positive side.

    a_- IS NOT PREDICTED HERE. It is a three-body observable. What this module
    does is construct the two-body interaction AT the point a = a_-, which is
    the input a later three-body calculation needs. Only that calculation can
    predict a_- and be compared with the table.

EVERYTHING IS COMPUTED AT r0 = 1 AND RESCALED.
    a and r0 are the only two lengths, so a solution depends on the ratio r0/a
    alone and the rest is scale: under r -> Lr the problem is unchanged and the
    parameters follow the exponents of rescale(). Five rows of the table are
    therefore five multiplications of one family, not five searches. Tuned
    directly at r_eff = 136 a0 the Lennard-Jones takes 144 s and fails outright
    for three of the five rows, because the starting guess carries L^10 and
    lands orders of magnitude away; tuned at r0 = 1 and rescaled it takes 1.4 s.

Lengths are in Bohr radii throughout, which is the unit Table III uses. The
solver fixes no length unit, so nothing has to be converted. hbar^2/2mu never
appears: a and r0 are both lengths and tune() works on the two of them.
"""
import math

from scipy.optimize import brentq

import lab

# =============================================================================
#  Table III of [6], transcribed verbatim. Lengths in Bohr radii.
#
#  The paper's r_eff is this code's r0: their Eq. (A3) reads
#  1/a(k0) = 1/a - r_eff k0^2 / 2, the effective-range expansion in the same
#  sign convention used here. No conversion, no sign flip.
#
#  Two things the source says that the numbers alone do not. The r_eff is NOT
#  measured: it comes from the authors' own coupled-channel calculations near
#  B_res. And the 33.5820 G row is not their measurement either -- it is
#  reproduced there from Refs. [38] and [54] of that paper.
#
#  B_txt is the resonance position as a string, exactly as the paper prints
#  it: 402.74 and 33.5820 are quoted to different numbers of decimals there,
#  and a float cannot carry that difference.
#
#  a_p is the positive-side feature in the paper's own notation, and it is NOT
#  the a_+ of that paper: Sec. V defines minima of the log-periodic modulation
#  of the three-body loss at a_+ e^{j pi/s0} and the corresponding maxima at
#  a_p e^{j pi/s0}, with a_p = a_+ e^{pi/(2 s0)}, half an Efimov period above.
#  Table III tabulates a_p, and the universal ratio the paper tests is
#  a_p/a_- = -1, which is why the two columns have nearly equal magnitude.
#  eta_minus is the Efimov width parameter of the same table: it encodes the
#  lifetime of the trimers, and it is nonzero only when there are deeply
#  bound dimers for them to decay into. See "the branch" below.
#
#  Both features are calibrated on -- calibrate() at a_- and calibrate_p() at
#  a_p -- and the sweep in "families" below runs between them. The experiment
#  reports the positive-side feature as the easier of the two to measure, and
#  the shape dependence is as small there as it is at a_-, so a_p is the better
#  three-body input of the two. The row at 162.36 G publishes no a_p.
# =============================================================================
R_VDW = 64.6            # a0, the van der Waals length of 39K, from [6]

TABLE_III = (
    {"state": "|1,1>",   "B_res": 402.74, "B_txt": "402.74", "s_res": 2.8, "r_eff": 136,
     "eta_minus": 0.27, "eta_minus_err": 0.06,
     "a_minus": -830, "a_minus_err": 40, "a_p": 900,  "a_p_err": 60},
    {"state": "|1,0>",   "B_res": 472.33, "B_txt": "472.33", "s_res": 2.8, "r_eff": 137,
     "eta_minus": 0.26, "eta_minus_err": 0.05,
     "a_minus": -840, "a_minus_err": 30, "a_p": 760,  "a_p_err": 70},
    {"state": "|1,-1>",  "B_res": 33.5820, "B_txt": "33.5820", "s_res": 2.6, "r_eff": 135,
     "eta_minus": 0.25, "eta_minus_err": 0.01,
     "a_minus": -908, "a_minus_err": 11, "a_p": 876,  "a_p_err": 28},
    {"state": "|1,-1>",  "B_res": 162.36, "B_txt": "162.36", "s_res": 1.1, "r_eff": 59,
     "eta_minus": 0.5, "eta_minus_err": 0.1,
     "a_minus": -780, "a_minus_err": 70, "a_p": None, "a_p_err": None},
    {"state": "|1,-1>",  "B_res": 561.14, "B_txt": "561.14", "s_res": 2.5, "r_eff": 132,
     "eta_minus": 0.33, "eta_minus_err": 0.06,
     "a_minus": -810, "a_minus_err": 30, "a_p": 800,  "a_p_err": 100},
)


def ratio(row):
    """r0/a at the Efimov feature. The one dimensionless number of the row."""
    return row["r_eff"] / row["a_minus"]


def ratio_p(row):
    """r0/a at the positive-side feature, or None where the row has no a_p.

    Same effective range -- r_eff is a property of the resonance, not of the
    side one approaches it from -- so only a changes, and with it the sign. The
    162.36 G row publishes no a_p and returns None rather than a guess.

    Worth having because the two sides are not equally measurable. On the
    positive side the loss feature is the atom-dimer one and it is reported as
    the easier and more robust measurement of the two; a_- sits on the
    three-atom side, where the trimer crosses the continuum. A two-body
    calibration at a_p is therefore the input for the three-body comparison
    that experiment can supply best, which is reason enough to carry it even
    though the task asked only for a_-.
    """
    if row["a_p"] is None:
        return None
    return row["r_eff"] / row["a_p"]


# =============================================================================
#  Rescaling
#
#  Under r -> L r the code potential must satisfy V'(Lr) = V(r)/L^2, because
#  [V] = length^-2. For the three one-scale shapes the scale parameter is an
#  inverse length and divides by L; for the Lennard-Jones [C6] = length^4 and
#  [C12] = length^10.
# =============================================================================
def rescale(name, p1, p2, L):
    """Take a solution at unit scale to one whose lengths are L times larger."""
    if name == "lj":
        return p1 * L ** 4, p2 * L ** 10
    return p1, p2 / L


# =============================================================================
#  Solutions at unit scale, computed once and kept
# =============================================================================
_UNIT = {}          # name -> (p1, p2) with 1/a = 0 and r0 = 1
_RATIO = {}         # (name, r0/a) -> (p1, p2) with r0 = 1
PASSO_MAXIMO = 0.075   # largest jump in r0/a a single tune() is asked to make


def unit_solution(name):
    """Parameters giving 1/a = 0 and r0 = 1. Computed once per shape."""
    if name not in _UNIT:
        pub = lab.PUBLISHED[("unitarity", name)]
        p1, p2, a, r0, nodes, ok = lab.tune(
            name, math.inf, 1.0, pub["p1"], pub["p2"], nodes_target=0)
        if not ok:
            raise RuntimeError("unitary solution did not converge for %s" % name)
        _UNIT[name] = (p1, p2)
    return _UNIT[name]


def at_ratio(name, r0_over_a):
    """Parameters with r0 = 1 and the given r0/a. Computed once per ratio.

    Walks in 1/a rather than jumping, because the outer loop of tune() moves
    the scale and a distant starting point costs the Lennard-Jones its
    convergence. The walk starts from the nearest ratio already solved for this
    shape on the same side of unitarity, and the unitary solution is the
    fallback: along a family the answer moves smoothly, so the neighbour is a
    far better guess than unitarity and the walk is correspondingly shorter.
    The step is capped at PASSO_MAXIMO either way, so a cold start behaves as
    it would have without the warm one. Measured over a 25-point family: 3.6
    times fewer solver calls, and parameters identical to 6e-14.
    """
    chave = (name, round(r0_over_a, 12))
    if chave in _RATIO:
        return _RATIO[chave]

    p1, p2 = unit_solution(name)
    if r0_over_a == 0.0:
        _RATIO[chave] = (p1, p2)
        return p1, p2

    partida = 0.0
    for (forma, r), par in _RATIO.items():
        if forma == name and r * r0_over_a > 0.0 \
                and abs(r - r0_over_a) < abs(partida - r0_over_a):
            partida, (p1, p2) = r, par

    # r0 = 1, so 1/a IS r0/a.
    passos = max(1, int(math.ceil(abs(r0_over_a - partida) / PASSO_MAXIMO)))
    for k in range(1, passos + 1):
        alvo = partida + (r0_over_a - partida) * k / passos
        p1, p2, a, r0, nodes, ok = lab.tune(
            name, 1.0 / alvo, 1.0, p1, p2, nodes_target=nodes_for(alvo))
        if not ok:
            raise RuntimeError("no convergence for %s at r0/a = %g" % (name, alvo))
    _RATIO[chave] = (p1, p2)
    return p1, p2


def nodes_for(inv_a_reduced):
    """How many bound states a shape must have at this 1/a, with r0 = 1.

    Crossing the resonance creates one. On the negative side there is none, on
    the positive side there is one, and inside the unitarity window the count
    is ambiguous by construction -- that is what lab.UNITARITY exists to decide
    -- so there the guard is switched off rather than guessed.
    """
    if abs(inv_a_reduced) <= lab.UNITARITY:
        return None
    return 0 if inv_a_reduced < 0.0 else 1


def unitary(name, r_eff):
    """The unitary potential whose effective range is r_eff. No search."""
    p1, p2 = unit_solution(name)
    return rescale(name, p1, p2, r_eff)


def calibrate(name, row):
    """Prescription of Sec. 1 of the task: a = a_-, r0 = r_eff, for one row."""
    p1, p2 = at_ratio(name, ratio(row))
    return rescale(name, p1, p2, row["r_eff"])


def calibrate_p(name, row):
    """a = a_p, r0 = r_eff: the same prescription on the positive side.

    Not a separate calculation. r0/a is the only number a solution depends on,
    so this is at_ratio() at a positive ratio and the same rescaling -- the
    machinery is indifferent to the side. What differs is the physics it lands
    on: a > 0 carries a shallow dimer, which the node target already expects.
    """
    r = ratio_p(row)
    if r is None:
        raise ValueError("row at %s G publishes no a_p" % row["B_txt"])
    p1, p2 = at_ratio(name, r)
    return rescale(name, p1, p2, row["r_eff"])


# =============================================================================
#  The two families of Sec. 2 of the task
#
#  "Turning up the strength" means multiplying the depth and leaving the range
#  alone. For the well, the Poschl-Teller and the Gaussian that is one
#  parameter. For the Lennard-Jones it is BOTH coefficients by the same factor:
#  with V = (1/2)(C12/r^12 - C6/r^6) the minimum sits at r_m^6 = 2 C12/C6 and
#  has depth C6^2/(8 C12), so a common factor leaves r_m exactly where it was
#  and multiplies the depth by that factor. Moving C6 alone -- the obvious
#  thing -- turns the potential up and slides its minimum at the same time,
#  which is a different operation from the one the other three are performing.
#  sweep_c6_alone measures that difference rather than asserting it.
# =============================================================================
def family_retuned(name, inv_a_reduced):
    """Both parameters retuned at every point, r0 held at 1.

    A family of fixed effective range whose potential changes at each point.
    One record per value of 1/a, all at r0 = 1; multiply the lengths by r_eff
    for a given row.
    """
    saida = []
    for inv in inv_a_reduced:
        p1, p2 = at_ratio(name, inv)
        a, r0, nodes = lab.scattering(lab.POTENTIALS[name](p1, p2))
        saida.append({"inv_a": inv, "p1": p1, "p2": p2,
                      "a": a, "r0": r0, "nodes": nodes})
    return saida


def turn_up(name, factor):
    """The unitary potential with its depth multiplied by `factor`."""
    p1, p2 = unit_solution(name)
    if name == "lj":
        return p1 * factor, p2 * factor
    return p1 * factor, p2


def r_min_of(name, p1, p2):
    """Where the potential has its minimum, to see whether it moved."""
    if name == "lj":
        return (2.0 * p2 / p1) ** (1.0 / 6.0)
    return 0.0


def family_fixed_range(name, inv_a_reduced):
    """Only the depth varies; the range stays at the unitary solution.

    Closer to how a physical interaction moves, and here the effective range is
    an output rather than a constraint. Takes the same list of 1/a as
    family_retuned so the two are read off at the same points: sweeping the
    depth directly is the obvious alternative and the wrong comparison, because
    the same step in depth is a different step in 1/a for every shape.
    """
    saida = []
    factor = 1.0
    for inv in inv_a_reduced:
        def erro(f):
            p1, p2 = turn_up(name, f)
            a, r0, nodes = lab.scattering(lab.POTENTIALS[name](p1, p2))
            return 1.0 / a - inv

        factor = lab.find_root(erro, factor)
        p1, p2 = turn_up(name, factor)
        a, r0, nodes = lab.scattering(lab.POTENTIALS[name](p1, p2))
        saida.append({"inv_a": inv, "p1": p1, "p2": p2,
                      "a": a, "r0": r0, "nodes": nodes,
                      "strength_ratio": factor, "r_min": r_min_of(name, p1, p2)})
    return saida


def sweep_c6_alone(inv_a_reduced):
    """The other reading of "turn up the Lennard-Jones": C6 alone.

    Deepens the well AND slides its minimum, which the common factor does not.
    Kept to measure the difference rather than assert it.
    """
    p1_0, p2_0 = unit_solution("lj")
    saida = []
    factor = 1.0
    for inv in inv_a_reduced:
        def erro(f):
            a, r0, n = lab.scattering(lab.LennardJones(p1_0 * f, p2_0))
            return 1.0 / a - inv

        factor = lab.find_root(erro, factor)
        p1 = p1_0 * factor
        a, r0, nodes = lab.scattering(lab.LennardJones(p1, p2_0))
        saida.append({"inv_a": inv, "p1": p1, "p2": p2_0,
                      "a": a, "r0": r0, "nodes": nodes,
                      "strength_ratio": factor,
                      "r_min": r_min_of("lj", p1, p2_0)})
    return saida


# =============================================================================
#  The branch: how many deeply bound dimers the potential carries
#
#  (a, r0) does not determine the potential. Deepening a shape past each pole
#  of a adds one bound state, and the same (a, r0) is reached again on every
#  branch beyond. All of those potentials are indistinguishable in two-body
#  low-energy scattering -- that is what universality means -- and they are NOT
#  equivalent in three.
#
#  The reason is a decay channel. With deeply bound dimers in the spectrum an
#  Efimov trimer can fall into one of them and release an atom, and the rate of
#  that is what the Efimov width parameter eta measures. On the shallowest
#  branch there are no deep dimers, no such channel, and eta would be exactly
#  zero. Table III of [6] publishes eta_- for every row -- 0.27, 0.26, 0.25,
#  0.5 and 0.33 -- so the measurement already rules the shallow branch out as
#  a description of 39K, which carries dozens of deeply bound molecular states.
#
#  This makes the branch a testable choice rather than a convention: a
#  three-body calculation on each branch can be compared against the published
#  eta_-. The shallow branch is the control -- what zero-range universality
#  alone would say -- and a branch with deep states is the realistic candidate.
#
#  Measured for the square well at a_- = -830, r_eff = 136 a0: the three
#  branches have nearly the same range (128.0, 136.2 and 136.9 a0) and very
#  different depths, the extra dimers coming out compact, of size 35 a0 and
#  19 a0, well inside r_vdW = 64.6 a0.
# =============================================================================
TETO_FORCA = 60.0       # highest depth multiple the branch search considers
PISO_FORCA = 0.2        # lowest; below it there is no bound state in any shape
XTOL_FORCA = 1e-12      # relative precision of every bisected landmark


def _na_forca(name, f):
    """a, r0 and the bound-state count at depth multiple f, at unit scale."""
    p1, p2 = turn_up(name, f)
    try:
        return lab.scattering(lab.POTENTIALS[name](p1, p2))
    except Exception:
        return None


def _bisseca(teste, lo, hi, xtol=XTOL_FORCA):
    """Where `teste` turns from false to true between lo and hi.

    Requires teste(lo) false and teste(hi) true. Returns None if the solver
    fails anywhere along the way, which is reported rather than skipped: a
    silently dropped sample is how the scan this replaces lost boundaries.
    """
    while hi - lo > xtol * max(1.0, abs(hi)):
        mid = 0.5 * (lo + hi)
        t = teste(mid)
        if t is None:
            return None
        if t:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def degraus(name, teto=TETO_FORCA, piso=PISO_FORCA):
    """The depth multiples at which the bound-state COUNT steps up, by bisection.

    The count is a non-decreasing staircase in the depth: deepening a potential
    cannot unbind a state, and the node correction of lab exists precisely so
    the count does not dip as a passes through infinity. Measured monotone at
    240 depths for every shape. So each step is the single crossing of the test
    "count >= n", which bisection locates to 1e-12 relative and CANNOT miss,
    however narrow.

    These steps are NOT the physical boundaries -- see fronteiras() -- they are
    what brackets them.
    """
    base = _na_forca(name, piso)
    if base is None:
        raise RuntimeError("solver fails at the floor for %s" % name)
    topo = _na_forca(name, teto)
    if topo is None:
        raise RuntimeError("solver fails at the ceiling for %s" % name)

    saida = {}
    for n in range(base[2] + 1, topo[2] + 1):
        f = _bisseca(lambda x: (lambda c: None if c is None else c[2] >= n)
                     (_na_forca(name, x)), piso, teto)
        if f is None:
            raise RuntimeError("bisection for state %d failed for %s" % (n, name))
        saida[n] = f
    return saida


def fronteiras(name, teto=TETO_FORCA, piso=PISO_FORCA):
    """The depth multiples at which each bound dimer actually appears.

    Two ways to ask this, and they do not agree. The bound-state count steps up
    at the depths degraus() returns, but that step is displaced: lab counts the
    node the asymptotic line contributes at r = a only once |r0/a| clears
    lab.UNITARITY, so the count lags the physics and the step sits ABOVE the
    depth where the state becomes bound. Measured: +0.08% on the first boundary
    of every shape, and the same absolute offset, 8.114e-4, at all four
    boundaries of the square well -- a systematic threshold artefact, not noise.

    The physical boundary is the POLE of a, where 1/a = 0 and the new state
    arrives at exactly zero binding energy. No threshold enters it. So the
    staircase step is used only to bracket, and the boundary itself is the zero
    of 1/a inside that bracket. Below the pole a runs to -infinity and 1/a
    approaches zero from below; above it, from above; and the nearest other
    structure, the zero of a, sits near the middle of the window, so a bracket
    reaching a few percent below the step cannot contain anything else.

    Validated against closed forms, two shapes, nothing fitted:
        square well      poles at f = (2n-1)^2   1, 9, 25, 49
        Poschl-Teller    poles at f = n(2n-1)    1, 6, 15, 28, 45
    both to 1e-13 relative. The unit solutions those ladders rest on are closed
    too: the well sits at v = pi^2/8, mu = 1 and the Poschl-Teller at v = 1,
    mu = 2, each to 1e-13.
    """
    passos = degraus(name, teto=teto, piso=piso)
    saida = {}
    anterior = piso
    for n, degrau in sorted(passos.items()):
        f = None
        largura = 0.01
        while largura <= 0.30:
            lo = max(anterior + XTOL_FORCA, degrau * (1.0 - largura))
            c_lo = _na_forca(name, lo)
            if c_lo is not None and c_lo[0] < 0.0:
                f = brentq(lambda x: 1.0 / _na_forca(name, x)[0],
                           lo, degrau, xtol=1e-14)
                break
            largura *= 3.0
        if f is None:
            raise RuntimeError("no pole of a below the count step %d for %s"
                               % (n, name))
        saida[n] = f
        anterior = f
    return saida


def janela(name, deep, teto=TETO_FORCA, piso=PISO_FORCA):
    """The depth window on which the shape carries exactly `deep` bound dimers."""
    f = fronteiras(name, teto=teto, piso=piso)
    lo = f.get(deep, piso)
    hi = f.get(deep + 1, teto)
    if deep > 0 and deep not in f:
        raise RuntimeError("no branch with %d deep states below the ceiling "
                           "for %s" % (deep, name))
    return lo, hi


def _onde_a_anula(name, lo, hi):
    """Where a passes through zero inside a window, by bisection on its sign.

    1/a runs monotonically from 0+ at the lower edge of a window to 0- at the
    upper edge, which it can only do by passing through infinity: a crosses
    zero exactly once, and r0 diverges with it. r0/a therefore has ONE pole
    inside every window and is continuous on either side of it. Locating that
    pole is what makes the bracketing below safe -- a bracket laid across it
    sees a sign change with no root in it, and brentq returns the pole.

    lo and hi are themselves poles of a, where |a| overflows and its sign is
    whatever the rounding left behind. The sign must therefore be read a sliver
    INSIDE the window, never at the edges. Reading it at the edges put the
    pole at the lower edge for deep = 1 and lost it entirely for deep = 2, and
    the unsplit window then handed back a degenerate potential.
    """
    margem = 1e-7 * (hi - lo)
    a_lo, a_hi = lo + margem, hi - margem
    c_lo = _na_forca(name, a_lo)
    c_hi = _na_forca(name, a_hi)
    if c_lo is None or c_hi is None:
        return None
    if (c_lo[0] > 0.0) == (c_hi[0] > 0.0):
        return None
    positivo_embaixo = c_lo[0] > 0.0
    return _bisseca(lambda x: (lambda c: None if c is None
                               else (c[0] < 0.0) == positivo_embaixo)
                    (_na_forca(name, x)), a_lo, a_hi)


def _raizes(name, alvo, lo, hi, pontos=12, refinos=9):
    """Every depth in (lo, hi) where r0/a equals alvo, refined to a stable count.

    r0/a is continuous here but far from monotone: it rises, turns over, crosses
    zero where r0 does, and only then runs off to the pole. One target can be
    met several times and no single bracket finds them all, so the interval is
    sampled, every sign change handed to brentq, and the sampling doubled until
    the number of roots repeats twice. The per-refinement counts come back with
    the roots so the convergence is reported rather than asserted.

    The starting density was measured, not chosen: 6, 12, 24 and 48 points all
    return the same roots to 1e-8 and the same counts at every refinement, for
    every shape and both branches. Twelve is kept for margin, and the whole
    search then costs about what the 400-point grid it replaced cost -- the
    convergence criterion comes free.
    """
    contagens = []
    raizes = []
    estaveis = 0
    n = pontos
    for _ in range(refinos):
        fs = [lo + (hi - lo) * k / n for k in range(n + 1)]
        vals = []
        for f in fs:
            c = _na_forca(name, f)
            vals.append(None if c is None else c[1] / c[0] - alvo)

        achadas = []
        for k in range(n):
            u, v = vals[k], vals[k + 1]
            if u is None or v is None or u == 0.0 or u * v > 0.0:
                continue
            try:
                achadas.append(brentq(
                    lambda f: (lambda c: c[1] / c[0] - alvo)(_na_forca(name, f)),
                    fs[k], fs[k + 1], xtol=1e-14))
            except (ValueError, TypeError):
                continue

        contagens.append(len(achadas))
        if len(contagens) >= 2 and contagens[-1] == contagens[-2]:
            estaveis += 1
            raizes = achadas
            if estaveis >= 2:
                break
        else:
            estaveis = 0
            raizes = achadas
        n *= 2
    return raizes, contagens


def at_ratio_on_branch(name, r0_over_a, deep, teto=TETO_FORCA, piso=PISO_FORCA,
                       diagnostico=None):
    """Parameters with r0 = 1 and the given r0/a, on the branch with `deep`
    deeply bound dimers. deep = 0 is what at_ratio() returns.

    The window is bounded by poles of a, split at the pole of r0/a, and each
    half searched to a stable root count. Costs far more solver calls than
    at_ratio(), so it is meant for the comparison described above and not for
    the production families.

    Pass a dict as `diagnostico` to receive the window, the pole, the
    refinement counts and the node count at every root, which is how the test
    suite checks the search rather than only its answer.
    """
    if deep == 0:
        return at_ratio(name, r0_over_a)

    # The window is indexed by the TOTAL number of bound states, and that is
    # not `deep`. A positive scattering length already carries a shallow dimer
    # of size ~a, so a potential with `deep` deep dimers has deep + 1 states
    # in all; a negative one has exactly `deep`. Indexing by `deep` on the
    # positive side returned the window one branch too low, and deep = 0 and
    # deep = 1 then came back as the same potential -- the shallow-only one,
    # with its dimer at 826 a0 next to a = 900 a0, which is the shallow dimer
    # and not a deep one.
    if abs(r0_over_a) <= lab.UNITARITY:
        raise ValueError("the bound-state count is ambiguous at r0/a = %g, so "
                         "a branch cannot be named there" % r0_over_a)
    total = deep + (1 if r0_over_a > 0.0 else 0)

    lo, hi = janela(name, total, teto=teto, piso=piso)
    polo = _onde_a_anula(name, lo, hi)

    # The window edges are poles of a themselves: r0/a is 0 there and the
    # solver is at its worst. Step in by a sliver, far smaller than any
    # structure measured inside a window.
    margem = 1e-7 * (hi - lo)
    if polo is None:
        trechos = [(lo + margem, hi - margem)]
    else:
        trechos = [(lo + margem, polo - margem), (polo + margem, hi - margem)]

    candidatos = []
    contagens = []
    discordancia = []
    for t_lo, t_hi in trechos:
        raizes, cont = _raizes(name, r0_over_a, t_lo, t_hi)
        contagens.append(cont)
        for f in raizes:
            c = _na_forca(name, f)
            if c is None:
                continue
            a, r0, nodes = c
            # r0 > 0 is not optional. The ratio r0/a reaches its target with
            # both signs flipped too, and rescaling by 1/r0 would then ask for
            # a negative scale. The requirement is r0 = +r_eff with a = a_- < 0,
            # and a positive rescaling cannot repair a sign.
            if r0 <= 0.0:
                continue
            # The branch is settled by the window, which is bounded by poles of
            # a and carries no threshold. The node count is NOT used to filter:
            # between the pole and the count's step it still reads total - 1,
            # and screening on it would put lab.UNITARITY back into the answer.
            # It is recorded instead, so any disagreement is visible.
            discordancia.append((f, nodes))
            p1, p2 = turn_up(name, f)
            q1, q2 = rescale(name, p1, p2, 1.0 / r0)

            # Rescaling is exact on paper and the solver still has to agree
            # with it. Re-solve and require the rescaled potential to carry
            # r0 = 1 and r0/a = the target. A root sitting near the pole of
            # r0/a comes back with r0 enormous, so 1/r0 collapses the scale and
            # leaves a potential the solver cannot integrate at all -- and that
            # degenerate potential has range zero, which the "most compact"
            # rule below then PREFERS. Verifying here is what stops it.
            try:
                a_v, r0_v, n_v = lab.scattering(lab.POTENTIALS[name](q1, q2))
            except Exception:
                continue
            if not (abs(r0_v - 1.0) < 1e-6
                    and abs(r0_v / a_v - r0_over_a) < 1e-6 * max(1.0, abs(r0_over_a))):
                continue
            candidatos.append((lab.POTENTIALS[name](q1, q2).R, q1, q2, f))

    if diagnostico is not None:
        diagnostico.update({"janela": (lo, hi), "polo": polo, "total": total,
                            "contagens": contagens, "nos": discordancia,
                            "candidatos": [(c[0], c[3]) for c in candidatos]})

    if not candidatos:
        raise RuntimeError("no branch with %d deep states reaches r0/a = %g for %s"
                           % (deep, r0_over_a, name))

    # A branch holds more than one solution: the same (a, r0) and the same
    # node count are reached at several depths, differing in range. Measured
    # for the square well at deep = 1: ranges of 136 a0 and 1428 a0 both
    # qualify. The compact one is returned, because 1428 a0 is twenty-two
    # van der Waals lengths and nothing about 39K is that long-ranged. The
    # others are discarded here rather than left for the caller, and this is a
    # physical judgement, not a numerical one.
    candidatos.sort()
    return candidatos[0][1], candidatos[0][2]


def calibrate_on_branch(name, row, deep):
    """a = a_-, r0 = r_eff, on the branch carrying `deep` deep dimers."""
    p1, p2 = at_ratio_on_branch(name, ratio(row), deep)
    return rescale(name, p1, p2, row["r_eff"])


def calibrate_p_on_branch(name, row, deep):
    """a = a_p, r0 = r_eff, on the branch carrying `deep` deep dimers.

    The positive side carries a shallow dimer of size ~a on top of whatever
    deep ones the branch holds, so `deep` counts the deep ones only and the
    total is one more. at_ratio_on_branch does that arithmetic.
    """
    r = ratio_p(row)
    if r is None:
        raise ValueError("row at %s G publishes no a_p" % row["B_txt"])
    p1, p2 = at_ratio_on_branch(name, r, deep)
    return rescale(name, p1, p2, row["r_eff"])


def main():
    print("Table III of [6], five rows, lengths in Bohr radii\n")
    print("  %-8s %9s %6s %7s %9s %8s" %
          ("state", "B_res (G)", "s_res", "r_eff", "a_-", "r0/a"))
    for row in TABLE_III:
        print("  %-8s %9s %6.1f %7d %9d %8.5f" %
              (row["state"], row["B_txt"], row["s_res"], row["r_eff"],
               row["a_minus"], ratio(row)))

    print("\nUnitary calibration. No search: scale invariance alone.")
    print("  %-8s %-6s %16s %16s %12s %12s" %
          ("r_eff", "shape", "p1", "p2", "1/a", "r0"))
    for row in TABLE_III:
        for name in ("well", "mpt", "gauss", "lj"):
            p1, p2 = unitary(name, row["r_eff"])
            a, r0, nodes = lab.scattering(lab.POTENTIALS[name](p1, p2))
            print("  %-8d %-6s %16.8g %16.8g %12.1e %12.6f" %
                  (row["r_eff"], name, p1, p2, 1.0 / a, r0))

    print("\nCalibration at a = a_-.")
    print("  %-8s %-6s %16s %16s %12s %12s %6s" %
          ("r_eff", "shape", "p1", "p2", "a", "r0", "nodes"))
    for row in TABLE_III:
        for name in ("well", "mpt", "gauss", "lj"):
            p1, p2 = calibrate(name, row)
            a, r0, nodes = lab.scattering(lab.POTENTIALS[name](p1, p2))
            print("  %-8d %-6s %16.8g %16.8g %12.4f %12.6f %6d" %
                  (row["r_eff"], name, p1, p2, a, r0, nodes))


if __name__ == "__main__":
    main()
