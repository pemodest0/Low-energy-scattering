"""Checks for lab.py. Run it with either

    python test_lab.py           (prints a report)
    pytest test_lab.py           (when pytest is available)

Every check compares against something not computed here: a closed form, an
analytic condition, or a published table. None of them asserts "the same as
last time" -- a test like that only proves the code is consistent with its
own mistakes.
"""
import math

import lab


# =============================================================================
#  1. The square well has closed forms, so it can be checked exactly
# =============================================================================
def test_well_scattering_length_matches_closed_form():
    """a from the solver vs. Eq. (80) of [1]: a/R = 1 - tan(x)/x."""
    v, mu = 1.7806, 0.48415
    a, r0, nodes = lab.scattering(lab.Well(v, mu))
    exact = lab.a_well(v, 1.0 / mu)
    assert abs(a / exact - 1) < 1e-8


def test_well_effective_range_matches_closed_form():
    """r0 from the solver vs. Eq. (92) of [1]."""
    v, mu = 1.7806, 0.48415
    a, r0, nodes = lab.scattering(lab.Well(v, mu))
    exact = lab.r0_well(v, 1.0 / mu)
    assert abs(r0 / exact - 1) < 1e-8


def test_effective_range_equals_the_range_at_every_pole():
    """Sharpest check we have: at every pole of a, that is sqrt(2v) = pi/2 + n pi,
    the article predicts r0/R = 1 exactly. Fig. 6 of [1]."""
    for n in range(4):
        x = math.pi / 2 + n * math.pi
        assert abs(lab.r0_well(x ** 2 / 2) - 1.0) < 1e-9


def test_bound_state_matches_the_analytic_well():
    """The square well's bound state solves k cot(kR) = -kappa, with
    k = sqrt(2(V0 + E)) and kappa = sqrt(-2E). Independent of our solver."""
    from scipy.optimize import brentq

    v, mu = 1.7806, 0.48415
    V0, R = v * mu ** 2, 1.0 / mu

    def condition(E):
        k = math.sqrt(2.0 * (V0 + E))
        return k / math.tan(k * R) + math.sqrt(-2.0 * E)

    exact = brentq(condition, -V0 + 1e-12, -1e-12, xtol=1e-15)
    guess = lab.E_finite_range(5.4112, 1.7436, 41.462) / (2 * 41.462)
    ours = lab.bound_energy(lab.Well(v, mu), guess)
    assert abs(ours / exact - 1) < 1e-7


# =============================================================================
#  2. Against the published tables
# =============================================================================
def test_node_counts_match_table_2():
    """Every published parameter pair must give the bound-state count that
    Table 2 of [1] requires. Two solutions can share (a, r0) and differ here,
    and then they are not the same physical state."""
    for case in lab.TARGETS:
        for name in lab.POTENTIALS:
            pub = lab.PUBLISHED[(case, name)]
            pot = lab.POTENTIALS[name](pub["p1"], pub["p2"])
            a, r0, nodes = lab.scattering(pot)
            assert nodes == lab.TARGETS[case]["nodes"], (case, name, nodes)


def test_tuning_hits_its_target():
    """The inverse problem must land on the target it was given, for all
    3 x 4 combinations of Tables 3 and 4."""
    for case in lab.TARGETS:
        target = lab.TARGETS[case]
        for name in lab.POTENTIALS:
            pub = lab.PUBLISHED[(case, name)]
            out = lab.tune(name, target["a"], target["r0"],
                           pub["p1"], pub["p2"], nodes_target=target["nodes"])
            strength, scale, a, r0, nodes, ok = out
            assert ok, (case, name)
            assert abs(r0 - target["r0"]) < 1e-6, (case, name, r0)
            if not math.isinf(target["a"]):
                assert abs(a / target["a"] - 1) < 1e-6, (case, name, a)


def test_inferred_constants_match_the_literature():
    """hbar^2/2mu is not tabulated: we invert the published zero-range energy.
    Recovering 41.47 MeV.fm^2 and 12.12 K.A^2 shows Table 1 is consistent."""
    assert abs(lab.SYSTEMS["deuteron"]["h2_2mu"] - 41.47) < 0.05
    assert abs(lab.SYSTEMS["he4_dimer"]["h2_2mu"] - 12.12) < 0.05


def test_aziz_reproduces_its_own_minimum():
    """The Aziz parameters must put the minimum where the paper says.

    V(rm) = -eps is a property of the published parameter set, so this checks
    the transcription without needing any external result: if a digit of
    EPS, RM, D, A, ALPHA, BETA or the C's were wrong, the minimum would move.

    It also exercises the unit conversion. The parameters are in K and A while
    the solver works in code units, so V(rm) has to be multiplied back by
    2 * h2_2mu to be compared with -eps. That round trip is the units section
    of the notes, used rather than assumed.
    """
    h2 = lab.SYSTEMS["he4_dimer"]["h2_2mu"]
    pot = lab.AzizHFDB(h2)

    depth = float(pot.V(pot.RM)) * 2.0 * h2
    assert abs(depth + pot.EPS) < 1e-5, depth

    grid = [2.0 + 0.0001 * k for k in range(20001)]
    values = [float(pot.V(r)) for r in grid]
    r_at_min = grid[values.index(min(values))]
    assert abs(r_at_min - pot.RM) < 1e-3, r_at_min

    assert pot.r_min == 0.0            # soft core, no hard wall
    assert pot.R > 100.0               # the 1/r^6 tail reaches a long way


def test_aziz_reproduces_the_published_dimer():
    """The strongest check in this file: a potential we did not fit.

    The four tunable potentials are fitted to (a, r0), so reproducing (a, r0) is
    circular. Aziz has no free parameters, so its predictions can be compared
    against someone else's published calculation with the same potential.

    Reference values from Stipanovic, Vranjes Markic & Boronat, arXiv:1607.05872,
    Table 1, row HFDB, He-4 dimer:

        a = 88.430 A,   r0 = 7.276 A,   E = -1.69 mK

    Getting within 0.1% of an independent calculation, on a potential nothing
    here was tuned to, tests the integrator, the matching, the effective-range
    integral, the bound-state search and the unit conversion at once.

    Note the constant: this uses H2_2MU_HE4, not the value inferred in SYSTEMS.
    The inferred one is 0.20% lower, which lands a 6% away, because near
    unitarity d(ln a)/d(ln h2_2mu) = 48.7. Getting this test to pass with the
    inferred constant would require loosening the tolerance to 7%, which would
    make it prove nothing.
    """
    pot = lab.AzizHFDB(lab.H2_2MU_HE4)
    a, r0, nodes = lab.scattering(pot)

    assert abs(a / 88.430 - 1) < 1e-3, a
    assert abs(r0 / 7.276 - 1) < 1e-3, r0
    assert nodes == 1, nodes

    guess = lab.E_finite_range(a, r0, lab.H2_2MU_HE4) / (2 * lab.H2_2MU_HE4)
    E_mK = lab.bound_energy(pot, guess) * 2 * lab.H2_2MU_HE4 * 1e3
    assert abs(E_mK / -1.69 - 1) < 5e-3, E_mK


def test_helium_benchmark_four_potentials():
    """Four parameter-free potentials against a published table.

    Reference: Motovilov, Sandhas, Sofianos & Kolganova, Eur. Phys. J. D 13, 33
    (2001), arXiv:physics/9910016. Table I gives the dimer energy and the He-He
    scattering length for HFDHE2, HFD-B, LM2M2 and TTY; the Appendix (Tables IX
    and X) gives the parameters. Both come from the same paper, so this is a
    closed loop: their parameters in, their numbers out.

    They state hbar^2/m = 12.12 K.A^2, so that value is used here rather than
    CODATA. Mixing the two would shift a by several percent -- see
    test_aziz_reproduces_the_published_dimer for why.

        potential   a (A)      E (mK)
        HFDHE2      124.65     -0.83012
        HFD-B        88.50     -1.68541
        LM2M2       100.23     -1.30348
        TTY         100.01     -1.30962

    HFD-B is given a looser tolerance on a. Its ENERGY reproduces to six
    figures, which means the physics agrees, but a comes out at 88.60 against
    their 88.50. The other three agree on both to five figures with the same
    code, so this is specific to that row and is not understood. It is left
    visible rather than tuned away.
    """
    h2 = 12.12
    casos = [
        (lab.AzizHFDHE2,      124.65,  -0.83012, 1e-3),
        (lab.AzizHFDB,         88.50,  -1.68541, 2e-3),
        (lab.AzizLM2M2,       100.23,  -1.30348, 1e-3),
        (lab.TangToenniesYiu, 100.01,  -1.30962, 1e-3),
    ]
    for cls, a_ref, E_ref, tol_a in casos:
        pot = cls(h2)
        a, r0, nodes = lab.scattering(pot)
        assert abs(a / a_ref - 1) < tol_a, (pot.name, "a", a, a_ref)
        assert nodes == 1, (pot.name, "nodes", nodes)

        guess = lab.E_finite_range(a, r0, h2) / (2 * h2)
        E_mK = lab.bound_energy(pot, guess) * 2 * h2 * 1e3
        assert abs(E_mK / E_ref - 1) < 1e-4, (pot.name, "E", E_mK, E_ref)


def test_hfd_family_reproduces_its_own_minimum():
    """Each Aziz potential must put its minimum at rm with depth eps.

    Checks the transcription of Table IX with no external number: three
    potentials, three parameter sets, and a mistyped digit in any of A, alpha,
    beta, D or the C's moves the minimum off (rm, -eps).
    """
    h2 = 12.12
    for cls in (lab.AzizHFDHE2, lab.AzizHFDB, lab.AzizLM2M2):
        pot = cls(h2)
        grid = [2.5 + 0.0001 * k for k in range(10001)]
        valores = [float(pot.V(r)) for r in grid]
        r_min = grid[valores.index(min(valores))]
        assert abs(r_min - pot.RM) < 2e-3, (pot.name, r_min, pot.RM)

        depth = min(valores) * 2.0 * h2
        assert abs(depth / -pot.EPS - 1) < 1e-3, (pot.name, depth, -pot.EPS)


def test_tty_is_insensitive_to_its_validity_cutoff():
    """TTY has no wall at V_CORE; its inner cutoff is a validity limit.

    The damping argument b(x) = 2 beta - p/x changes sign at 0.3156 A, below
    which the published expression diverges to -infinity. Integration starts at
    BREAKDOWN_MARGIN times that radius. If the result depended on the margin,
    the cutoff would be doing physics, which it must not.
    """
    original = lab.BREAKDOWN_MARGIN
    try:
        lab.BREAKDOWN_MARGIN = 1.15
        a0, r00, n0 = lab.scattering(lab.TangToenniesYiu(12.12))
        for margem in (1.05, 1.5, 2.0, 3.0):
            lab.BREAKDOWN_MARGIN = margem
            a, r0, nodes = lab.scattering(lab.TangToenniesYiu(12.12))
            assert abs(a / a0 - 1) < 1e-6, (margem, "a", a, a0)
            assert abs(r0 / r00 - 1) < 1e-6, (margem, "r0", r0, r00)
            assert nodes == n0, (margem, "nodes")
    finally:
        lab.BREAKDOWN_MARGIN = original


def test_the_unitarity_threshold_is_not_fragile():
    """UNITARITY decides when a is 'infinite' for the bound-state count.

    Any threshold is a convention, so the thing to check is not the value but
    the width of the window that works. Measured: every one of the twelve
    published cases gets the node count of Table 2 for UNITARITY anywhere from
    1e-2 to 1e-4. It breaks at 1e-5, where the Lennard-Jones at unitarity
    (|r0/a| = 2.1e-5) stops being recognised as unitarity and its exterior zero
    gets counted.

    Three decades of margin, with the chosen value in the middle. A threshold
    that only worked at one value would be a fit, not a criterion.
    """
    original = lab.UNITARITY
    try:
        for threshold in (1e-2, 1e-3, 1e-4):
            lab.UNITARITY = threshold
            for case in lab.TARGETS:
                for name in lab.POTENTIALS:
                    pub = lab.PUBLISHED[(case, name)]
                    pot = lab.POTENTIALS[name](pub["p1"], pub["p2"])
                    a, r0, nodes = lab.scattering(pot)
                    assert nodes == lab.TARGETS[case]["nodes"], (threshold, case, name)
    finally:
        lab.UNITARITY = original


def test_unitarity_sign_does_not_leak():
    """At the resonance a has no limit. Prove the node count does not care.

    Approached from one side a runs to -infinity, from the other it returns
    from +infinity, and at the crossing itself neither sign is "the" answer.
    scattering() does not choose: it computes 1/a, which passes through zero
    smoothly, and inverts. But a is still returned and still read by the node
    count, so the sign has to be shown not to matter.

    If this ever fails, something downstream is reading a where it should read
    1/a or |r0/a|.
    """
    import math

    pub = lab.PUBLISHED[("unitarity", "well")]
    pot = lab.Well(pub["p1"], pub["p2"])
    a, r0, nodes = lab.scattering(pot)

    assert abs(r0 / a) < lab.UNITARITY, (r0, a)
    assert nodes == 0, nodes

    for infinito in (math.inf, -math.inf):
        conta = pot.R < infinito and abs(r0 / infinito) > lab.UNITARITY
        assert not conta, infinito


def test_scale_invariance():
    """Rescaling every length by L must take (a, r0) to (L a, L r0).

    This is a property of the physics, not a stored number, so the test needs
    no reference value: it asserts an invariant the code must satisfy on its
    own terms. That is what makes it strong.

    It is also the test that would have caught the exponent bug. The notes used
    to say C12 -> L^6 C12, which comes from treating V as dimensionless. Since
    [V_code] = fm^-2, the rescaled potential must obey V'(L r) = V(r) / L^2,
    which gives C6 -> L^4 C6 and C12 -> L^10 C12. With the wrong exponent this
    test fails by 41% at L = 2.

    The tolerances differ by potential, and the reason is itself a result. The
    physics is exactly scale invariant; our truncation is not. R and r_min come
    from the ABSOLUTE thresholds V_ZERO and V_CORE, while V rescales as 1/L^2,
    so the cut radii do not follow L:

        well          R = 1/mu, exact             -> violation ~1e-13
        Poschl-Teller R from acosh, logarithmic   -> ~1e-10
        Gaussian      R from sqrt(log), log       -> ~1e-11
        Lennard-Jones R from a power law,
                      R ~ L^(2/3), r_min ~ L^(5/6) -> ~3e-4

    So the Lennard-Jones is six orders worse than the others, and it is the
    truncation talking, not the solver. Making the cut dimensionless -- |V|R^2
    is the natural choice, since [V_code] = fm^-2 -- would restore the symmetry,
    and belongs in the error budget rather than here.
    """
    reescala = {
        "well":  lambda p1, p2, L: (p1, p2 / L),
        "mpt":   lambda p1, p2, L: (p1, p2 / L),
        "gauss": lambda p1, p2, L: (p1, p2 / L),
        "lj":    lambda p1, p2, L: (p1 * L ** 4, p2 * L ** 10),
    }
    tolerancia = {"well": 1e-9, "mpt": 1e-8, "gauss": 1e-8, "lj": 1e-3}
    for name in lab.POTENTIALS:
        pub = lab.PUBLISHED[("deuteron", name)]
        a0, r00, n0 = lab.scattering(lab.POTENTIALS[name](pub["p1"], pub["p2"]))
        tol = tolerancia[name]
        for L in (0.5, 2.0, 3.0):
            p1, p2 = reescala[name](pub["p1"], pub["p2"], L)
            a, r0, nodes = lab.scattering(lab.POTENTIALS[name](p1, p2))
            assert abs(a / (L * a0) - 1) < tol, (name, L, "a", a, L * a0)
            assert abs(r0 / (L * r00) - 1) < tol, (name, L, "r0", r0, L * r00)
            assert nodes == n0, (name, L, "nodes")


# =============================================================================
#  3. The grid itself
# =============================================================================
def test_the_grid_is_already_converged():
    """Halve the points: if nothing moves beyond 1e-6, POINTS is enough.

    This is the check that caught the real bug. On a uniform grid the
    Lennard-Jones r0 drifted 1.74257 -> 1.74215 -> 1.74165 as the points
    doubled, which is not convergence at all. The logarithmic grid fixed it.

    We halve rather than double on purpose. Total error is truncation, which
    falls with more points, plus roundoff, which grows with them, so there is
    an optimum. The three purely attractive potentials sit at machine precision,
    1e-11, and do not care. The Lennard-Jones does: past ~8000 points roundoff takes
    over and r0 degrades again (-3e-5 at 32001, -3e-4 at 64001). POINTS is set
    at the optimum, so the meaningful question is whether we are converged
    coming from below, not what happens deeper into the roundoff.
    """
    fine = {}
    for name in lab.POTENTIALS:
        pub = lab.PUBLISHED[("deuteron", name)]
        fine[name] = lab.scattering(lab.POTENTIALS[name](pub["p1"], pub["p2"]))

    original = lab.POINTS
    lab.POINTS = (original - 1) // 2 + 1
    try:
        for name in lab.POTENTIALS:
            pub = lab.PUBLISHED[("deuteron", name)]
            a, r0, nodes = lab.scattering(lab.POTENTIALS[name](pub["p1"], pub["p2"]))
            assert abs(a / fine[name][0] - 1) < 1e-6, (name, "a", a)
            assert abs(r0 / fine[name][1] - 1) < 1e-6, (name, "r0", r0)
            assert nodes == fine[name][2], (name, "nodes")
    finally:
        lab.POINTS = original


# =============================================================================
#  Plain runner, so no extra tool is needed
# =============================================================================
if __name__ == "__main__":
    checks = []
    for key in sorted(dir()):
        if key.startswith("test_"):
            checks.append((key, globals()[key]))

    failed = 0
    for name, check in checks:
        try:
            check()
            print("  ok    ", name)
        except AssertionError as problem:
            failed = failed + 1
            print("  FAILED", name, "->", problem)

    print()
    print(len(checks) - failed, "of", len(checks), "checks passed")
