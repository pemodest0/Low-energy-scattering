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
    an optimum. The three smooth potentials sit at machine precision (1e-11)
    and do not care. The Lennard-Jones does: past ~8000 points roundoff takes
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
