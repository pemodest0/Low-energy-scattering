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
#  It is carried because it is the natural far end of the path in "families"
#  below; the row at 162.36 G has none published.
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
def _ramos(name, deep, passos=400, teto=60.0):
    """Walk the strength upward at unit scale, recording r0/a and the nodes.

    tune() cannot be pointed at a branch: its bracket widens until the sign
    changes and it finds its way back to the shallowest solution whatever it
    is started from. So the branch is located first, by scanning, and only
    then handed to a root find that stays inside it.
    """
    p1_0, p2_0 = unit_solution(name)
    anterior = None
    fora = []
    for k in range(1, passos + 1):
        f = teto * k / passos
        p1, p2 = turn_up(name, f)
        try:
            a, r0, nodes = lab.scattering(lab.POTENTIALS[name](p1, p2))
        except Exception:
            anterior = None
            continue
        if anterior is not None and anterior[1] == nodes == deep:
            fora.append((anterior[0], f, anterior[2], r0 / a))
        anterior = (f, nodes, r0 / a)
    return fora


def at_ratio_on_branch(name, r0_over_a, deep, **kwargs):
    """Parameters with r0 = 1 and the given r0/a, on the branch with `deep`
    deeply bound dimers. deep = 0 is what at_ratio() returns.

    Costs a scan of the strength, so it is far slower than at_ratio() and is
    meant for the comparison described above, not for the production families.
    """
    if deep == 0:
        return at_ratio(name, r0_over_a)

    def razao(f):
        p1, p2 = turn_up(name, f)
        a, r0, nodes = lab.scattering(lab.POTENTIALS[name](p1, p2))
        return r0 / a - r0_over_a

    candidatos = []
    for lo, hi, raz_lo, raz_hi in _ramos(name, deep, **kwargs):
        if (raz_lo - r0_over_a) * (raz_hi - r0_over_a) > 0.0:
            continue
        f = brentq(razao, lo, hi, xtol=1e-14)
        p1, p2 = turn_up(name, f)
        a, r0, nodes = lab.scattering(lab.POTENTIALS[name](p1, p2))
        # r0 > 0 is not optional. The ratio r0/a reaches its target with both
        # signs flipped too, and rescaling by 1/r0 would then ask for a
        # negative scale. The requirement is r0 = +r_eff with a = a_- < 0, and
        # a positive rescaling cannot repair a sign.
        if nodes != deep or r0 <= 0.0:
            continue
        q1, q2 = rescale(name, p1, p2, 1.0 / r0)
        candidatos.append((lab.POTENTIALS[name](q1, q2).R, q1, q2))

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
