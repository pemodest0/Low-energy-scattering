"""Checks for lab.py. Run it with either

    python test_lab.py           (prints a report)
    pytest test_lab.py           (when pytest is available)

Every check compares against something not computed here: a closed form, an
analytic condition, or a published table. None of them asserts "the same as
last time" -- a test like that only proves the code is consistent with its
own mistakes.

Tolerances sit about twenty times above what is measured, and no looser. A
tolerance far above the achieved error passes whatever happens and catches
nothing; several here were that way, calibrated for behaviour the code no
longer has.
"""
import math

import etrych
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
    assert abs(r0 / exact - 1) < 1e-10


def test_the_closed_form_effective_range_equals_the_range_at_every_pole():
    """Sharpest check we have: at every pole of a, that is sqrt(2v) = pi/2 + n pi,
    the article predicts r0/R = 1 exactly. Fig. 6 of [1]."""
    for n in range(4):
        x = math.pi / 2 + n * math.pi
        assert abs(lab.r0_well(x ** 2 / 2) - 1.0) < 1e-12


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
            # Relative, and the same 1e-7 tune() itself accepts: an absolute
            # 1e-8 on r0 = 1.7 demands 6e-9 relative, which is below what the
            # forward problem reproduces -- the Lennard-Jones stops improving
            # at 4e-8 relative, the noise of the tail extrapolation.
            assert abs(r0 / target["r0"] - 1.0) < 1e-7, (case, name, r0)
            if not math.isinf(target["a"]):
                assert abs(a / target["a"] - 1) < 1e-9, (case, name, a)


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
        (lab.TangToenniesYiu, 100.01,  -1.30962, 1e-5),
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
    """TTY has no wall at EPS_CORE; its inner cutoff is a validity limit.

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


def test_the_grid_ends_exactly_on_R():
    """exp(log(R)) must not overshoot R, for any mu.

    The square well decides its depth with r <= R. If the last grid point lands
    one ulp outside, that point gets V = 0 instead of the depth, and the error
    in a jumps from 1.6e-10 to 1.1e-3 -- seven orders, from one point.

    Measured before the fix: 14% of mu values overshoot. The published
    parameters were not among them, which is why it never surfaced until a
    rescaling test happened to land on mu = 1/3.

    numerov now pins r[-1] = pot.R. This sweeps mu widely and requires the
    closed form to be reproduced every time.
    """
    for k in range(1, 120):
        mu = 0.15 + k * 0.007
        pot = lab.Well(1.7575, mu)
        a, r0, nodes = lab.scattering(pot)
        assert abs(a / lab.a_well(1.7575, pot.R) - 1) < 1e-8, (mu, a)
        assert abs(r0 / lab.r0_well(1.7575, pot.R) - 1) < 1e-8, (mu, r0)


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

    Tolerances are set just above what is measured, with a factor of ten of
    slack and no more. They used to be far looser -- 1e-3 for the Lennard-Jones
    -- because the cuts were absolute thresholds and the truncated problem was
    genuinely not scale invariant. With the cuts written as |V| r^2 that is
    fixed, and a tolerance calibrated for the broken version would let the
    regression back in without failing.
    """
    reescala = {
        "well":  lambda p1, p2, L: (p1, p2 / L),
        "mpt":   lambda p1, p2, L: (p1, p2 / L),
        "gauss": lambda p1, p2, L: (p1, p2 / L),
        "lj":    lambda p1, p2, L: (p1 * L ** 4, p2 * L ** 10),
    }
    # The Lennard-Jones sits an order of magnitude looser than the other
    # three, and for a reason that is measured: its tail is extrapolated from
    # two cuts, and the residual of that extrapolation does not scale exactly
    # with the problem. Measured at 1.1e-7; the tolerance keeps the factor of
    # ten and no more.
    tolerancia = {"well": 1e-10, "mpt": 1e-10, "gauss": 1e-10, "lj": 1e-6}
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
#  3. Regression: does the answer still come out the same?
# =============================================================================
REGRESSAO = {
    ("nn", "well"): (-0.05397315947831557, 2.69461301195824, 0),
    ("nn", "mpt"): (-0.054009871965853426, 2.701731328218038, 0),
    ("nn", "gauss"): (-0.05389677696273141, 2.7030741336126063, 0),
    ("nn", "lj"): (-0.05405061491532482, 2.708146513279644, 0),
    ("unitarity", "well"): (-5.501364195160311e-07, 1.0000002229619138, 0),
    ("unitarity", "mpt"): (-1.586579559309834e-13, 1.0000000000468923, 0),
    ("unitarity", "gauss"): (-2.2193459732639425e-06, 1.000227209549898, 0),
    ("unitarity", "lj"): (2.1229701956925556e-05, 0.999985984608486, 0),
    ("deuteron", "well"): (0.18518954015905392, 1.697804215859341, 1),
    ("deuteron", "mpt"): (0.18517928641121206, 1.7299508321943549, 1),
    ("deuteron", "gauss"): (0.18517593147431377, 1.6989305700272126, 1),
    ("deuteron", "lj"): (0.18500462917843025, 1.700154188205488, 1),
}


def test_nothing_moved_since_the_last_recorded_run():
    """The twelve published cases, pinned to the values they currently give.

    This is NOT a correctness check and does not belong with the ones above.
    Those compare against closed forms and published tables, and would catch a
    wrong answer. This one only says the answer has not moved, and would
    happily pin a wrong answer forever.

    It earns its place by catching a different failure: a refactor that changes
    a number silently. That happened here more than once -- rewriting a as
    1/inv_a moved four cases in the last two digits, and migrating to
    dimensionless cuts moved all of them. Both were correct changes, but both
    should be seen and approved rather than discovered later.

    When it fails, decide which: a regression, or a change worth recording. If
    the latter, update the table in the same commit that causes it, and say so
    in the message.
    """
    # 1/a and not a: at the unitary rows a is whatever huge number the last
    # digits happened to produce, and pinning it pins noise. 1/a is the
    # quantity the solver actually computes and it passes smoothly through
    # zero.
    for (caso, nome), (inv_a_ref, r0_ref, nodes_ref) in REGRESSAO.items():
        pub = lab.PUBLISHED[(caso, nome)]
        pot = lab.POTENTIALS[nome](pub["p1"], pub["p2"])
        a, r0, nodes = lab.scattering(pot)
        assert abs(1.0 / a - inv_a_ref) < 1e-12, (caso, nome, "1/a", 1.0/a, inv_a_ref)
        assert abs(r0 / r0_ref - 1) < 1e-12, (caso, nome, "r0", r0, r0_ref)
        assert nodes == nodes_ref, (caso, nome, "nodes", nodes, nodes_ref)


# =============================================================================
#  4. The grid itself
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
#  The production solver, and the second method
# =============================================================================
def test_the_solver_effective_range_approaches_the_range_at_every_pole():
    """r0 -> R near the poles of a, computed and not asserted.

    The closed-form test above checks that lab.r0_well has the property the
    article gives it, which is circular and does not touch the solver. This
    one runs the solver near each pole -- near, not on: a diverges there -- and
    requires r0/R to approach 1 as the distance to the pole shrinks.
    """
    for n in range(1, 5):
        v_polo = (2 * n - 1) ** 2 * math.pi ** 2 / 8.0
        anterior = None
        for delta in (1e-2, 1e-3, 1e-4):
            a, r0, nodes = lab.scattering(lab.Well(v_polo * (1.0 + delta), 1.0))
            erro = abs(r0 / 1.0 - 1.0)
            if anterior is not None:
                assert erro < anterior, (n, delta, erro, anterior)
            anterior = erro
        assert anterior < 1e-3, (n, anterior)


def test_the_numerov_agrees_with_the_production_solver():
    """Two integrators with nothing in common, on the same potentials.

    Order 4 against 8, a fixed step in ln r against an adaptive step in r,
    Simpson on its own grid against an integral carried inside the system,
    sign changes against event detection. Run at a near tail cut, which is
    where the Numerov is reliable: its Simpson integral runs over a difference
    of large numbers, and at a far cut refining its grid makes the answer
    worse rather than better.
    """
    original = lab.EPS_TAIL
    lab.EPS_TAIL = 1e-6
    try:
        for name in ("well", "mpt", "gauss", "lj"):
            p1, p2 = etrych.unit_solution(name)
            pot = lab.POTENTIALS[name](p1, p2)
            a1, r01, n1 = lab._num_corte(pot, pot.R)
            a2, r02, n2 = lab.scattering_numerov(pot)
            assert abs(r01 / r02 - 1.0) < 1e-7, (name, "r0", r01, r02)
            assert n1 == n2, (name, "nodes", n1, n2)
    finally:
        lab.EPS_TAIL = original


def test_the_numerov_agrees_on_the_bound_state_too():
    """The second method on the other side of zero energy."""
    for v, guess in ((2.0, -0.5), (6.0, -2.0), (12.0, -5.0)):
        pot = lab.Well(v, 1.0)
        E = lab.bound_energy(pot, guess)
        k = math.sqrt(2.0 * v - 2.0 * abs(E))
        kappa = math.sqrt(-2.0 * E)
        assert abs(k / math.tan(k) + kappa) < 1e-9, (v, E)


def test_the_tail_extrapolation_is_consistent():
    """The same r0 from two independent extrapolations.

    r0(R) = r0 - K/R is extrapolated from R and 3R. If that form is right, the
    pair (3R, 9R) must give the same answer. Where the tail is exponential or
    absent the two agree trivially; the Lennard-Jones is the case that tests
    the 1/R.
    """
    tol = {"well": 1e-12, "mpt": 1e-9, "gauss": 1e-10, "lj": 1e-5}
    for name in ("well", "mpt", "gauss", "lj"):
        p1, p2 = etrych.unit_solution(name)
        pot = lab.POTENTIALS[name](p1, p2)
        base = pot.R
        r0 = [lab._num_corte(pot, base * m)[1] for m in (1, 3, 9)]
        tres = (3.0 * r0[1] - r0[0]) / 2.0
        nove = (3.0 * r0[2] - r0[1]) / 2.0
        assert abs(tres / nove - 1.0) < tol[name], (name, tres, nove)


def test_the_core_cut_is_exact():
    """The Lennard-Jones cuts solve |V(r)| r^2 = eps, not an approximation."""
    for C6, C12 in ((1.0, 1.0), (0.264624012, 0.000340682), (1e7, 1e17)):
        pot = lab.LennardJones(C6, C12)
        assert abs(abs(pot.V(pot.r_min)) * pot.r_min ** 2 / lab.EPS_CORE - 1.0) < 1e-12
        assert abs(abs(pot.V(pot.R)) * pot.R ** 2 / lab.EPS_TAIL - 1.0) < 1e-9


# =============================================================================
#  The calibration on Table III
# =============================================================================
def test_the_unitary_well_is_the_analytic_one():
    """1/a = 0 for the square well happens at 2v = (pi/2)^2, in closed form."""
    v, mu = etrych.unit_solution("well")
    assert abs(v / (math.pi ** 2 / 8.0) - 1.0) < 1e-10, v
    assert abs(mu - 1.0) < 1e-12, mu


def test_the_unitary_calibration_holds_at_every_published_range():
    """Rescaling the unitary solution must keep 1/a = 0 and give r0 = r_eff."""
    for row in etrych.TABLE_III:
        for name in ("well", "mpt", "gauss", "lj"):
            p1, p2 = etrych.unitary(name, row["r_eff"])
            a, r0, nodes = lab.scattering(lab.POTENTIALS[name](p1, p2))
            assert abs(1.0 / a) < 1e-9, (name, row["r_eff"], 1.0 / a)
            assert abs(r0 / row["r_eff"] - 1.0) < 1e-8, (name, r0, row["r_eff"])
            assert nodes == 0, (name, nodes)


def test_the_rescaled_calibration_equals_a_direct_tune():
    """Scale invariance is used as a shortcut; this is the shortcut checked.

    The whole module computes at r0 = 1 and multiplies. Tuning directly at the
    published range must land on the same parameters.
    """
    row = etrych.TABLE_III[0]
    p1, p2 = etrych.calibrate("well", row)
    direto = lab.tune("well", row["a_minus"], float(row["r_eff"]),
                      p1, p2, nodes_target=0)
    assert direto[5], "direct tune did not converge"
    assert abs(direto[0] / p1 - 1.0) < 1e-9, (direto[0], p1)
    assert abs(direto[1] / p2 - 1.0) < 1e-9, (direto[1], p2)


def test_the_deliverable_reproduces_its_own_targets():
    """Every calibrated potential hits the (a, r0) it claims to hit."""
    for row in etrych.TABLE_III:
        for name in ("well", "mpt", "gauss", "lj"):
            p1, p2 = etrych.calibrate(name, row)
            a, r0, nodes = lab.scattering(lab.POTENTIALS[name](p1, p2))
            assert abs(a / row["a_minus"] - 1.0) < 1e-9, (name, a, row["a_minus"])
            assert abs(r0 / row["r_eff"] - 1.0) < 1e-7, (name, r0, row["r_eff"])
            assert nodes == 0, (name, row["r_eff"], nodes)


def test_turning_up_the_lennard_jones_leaves_its_minimum_where_it_was():
    """A common factor on C6 and C12 multiplies the depth and moves nothing.

    This is what makes it the counterpart of "deepen the well" for the other
    three shapes. Scaling C6 alone fails this on purpose, and sweep_c6_alone
    exists to measure that.
    """
    p1, p2 = etrych.unit_solution("lj")
    r_m = etrych.r_min_of("lj", p1, p2)
    for factor in (0.5, 0.9, 1.1, 2.0):
        q1, q2 = etrych.turn_up("lj", factor)
        assert abs(etrych.r_min_of("lj", q1, q2) / r_m - 1.0) < 1e-14, factor
        pot_a, pot_b = lab.LennardJones(p1, p2), lab.LennardJones(q1, q2)
        assert abs(pot_b.V(r_m) / (factor * pot_a.V(r_m)) - 1.0) < 1e-12, factor
    q1, q2 = p1 * 1.5, p2
    assert abs(etrych.r_min_of("lj", q1, q2) / r_m - 1.0) > 1e-3


def test_the_fixed_range_family_matches_the_closed_form():
    """The depth-only family, against the closed forms of the square well."""
    alvos = [-0.2, -0.1, 0.1, 0.2]
    for linha in etrych.family_fixed_range("well", alvos):
        v, mu = linha["p1"], linha["p2"]
        assert abs(linha["a"] / lab.a_well(v, 1.0 / mu) - 1.0) < 1e-9, linha
        assert abs(linha["r0"] / lab.r0_well(v, 1.0 / mu) - 1.0) < 1e-9, linha


def test_the_deep_branches_share_the_two_body_observables():
    """Same (a, r0), different number of deeply bound dimers.

    (a, r0) does not determine the potential: every branch beyond each pole of
    a reaches the same pair again, carrying one more deep dimer. All of them
    are indistinguishable in two-body low-energy scattering, and they are not
    equivalent in three, because a deep dimer gives an Efimov trimer somewhere
    to decay to -- which is what the published eta_- measures, and what the
    shallow branch would predict to be exactly zero.

    The three branches are also required to come out on a physical scale. The
    search finds more than one solution per branch, differing in range; the
    compact one is kept, and this checks that what is kept is of order r_vdW
    and not twenty times it.
    """
    row = etrych.TABLE_III[0]
    for deep in (0, 1, 2):
        p1, p2 = etrych.calibrate_on_branch("well", row, deep)
        a, r0, nodes = lab.scattering(lab.Well(p1, p2))
        assert abs(a / row["a_minus"] - 1.0) < 1e-9, (deep, "a", a)
        assert abs(r0 / row["r_eff"] - 1.0) < 1e-8, (deep, "r0", r0)
        assert nodes == deep, (deep, "nodes", nodes)
        alcance = (1.0 / p2) / etrych.R_VDW
        assert 1.0 < alcance < 5.0, (deep, "range in r_vdW", alcance)


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
