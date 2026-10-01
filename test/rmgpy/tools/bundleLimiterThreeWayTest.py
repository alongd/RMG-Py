"""Three-way cross-check of the near-exhaustion bundle-limiter law.

The law lives in three places on purpose (see ``bundle_limiter_three_way.py``
for the inventory).  ``TestConeMarginBandParity`` in
``polymerMomentsConsumerTest.py`` already pins copies 1 and 2 across the
stage-2 M band on NON-end-group rows.  This file extends the same idea to the
third copy -- ``_s_eff`` in ``test/rmgpy/solver/solverPolymerTest.py``, the
expected-value generator ~30 assertions there compare solver output against --
and to the stage-1 E band, the end-group surface, the early-return branches and
the mu3 closure.  It does not weaken or duplicate that pin; it brackets it.

THE BAR
-------
Exact equality is the right bar only where the three evaluate the same
operations in the same order, and a round number is no bar anywhere.  So:

* **bitwise** wherever the operation order is provably identical.  On the M
  axis with stage 1 inactive, that is the whole path: ``s_free`` is ``s_base``
  by early return in all three; ``q10 = mu1 - mu0`` is EXACT by Sterbenz's
  lemma (mu0 and mu1 within a factor of two throughout the M band, so their
  difference is representable and the cancellation costs nothing); ``V_poly ==
  1.0`` makes copy 1's ``V_poly*(b1c - 1)`` and copy 3's ``(b1c - 1)`` the same
  operation; and ``softmin_p([s_free, s_cone])`` is written the same way, in
  the same term order, in all three.  Nothing there can differ by rounding, so
  anything short of equality is a divergence.  Same for every branch whose
  answer is an exact ``0.0`` or an exact ``s_base`` passthrough.

* **derived** where stage 1 is active.  There copies 1/2 fold ``s_base`` into
  the bundle cap as a NESTED soft-min while copy 3 evaluates the FLAT soft-min
  over the same term set; the reciprocal p-norm is exactly associative, so the
  two agree algebraically and differ only in rounding.  Copy 1 also forms E as
  ``softmin_p(mu_k/f_k)`` (per-moment floors) where copies 2/3 form it as
  ``softmin_p(mu_k)/f`` (scalar floor) -- identical on a uniform-atol deck up
  to the division's rounding, and reaching the answer only through the
  smoothstep weight.  The bound is computed per sample point from the actual
  operation counts in ``bundle_limiter_three_way.derived_abs_bound``; it is
  never a literal here.

KNOWN DIVERGENCES
-----------------
The sweep found three places where the three copies genuinely disagree.  They
are NOT repaired here (deciding which of three implementations is right is not
this file's call, and a silent repair would destroy the measurement); they are
characterized in ``TestKnownThreeWayDivergences`` with the numbers, so that a
later change to any copy has to come past them.
"""

import numpy as np
import pytest

import rmgpy.solver.polymer as solver_mod

import bundle_limiter_three_way as tw
from bundle_limiter_three_way import (
    FLOOR,
    RTOL,
    U,
    V_POLY,
    Case,
    derived_abs_bound,
    evaluate,
    seff_cap_terms,
    solver_cap_terms,
)
from numpy_moments_consumer import softmin_p

E_LO = solver_mod.BUNDLE_LIMITER_E_LO
E_HI = solver_mod.BUNDLE_LIMITER_E_HI
M_LO = solver_mod.CONE_MARGIN_M_LO
M_HI = solver_mod.CONE_MARGIN_M_HI


# ---------------------------------------------------------------------------
# Sample-point families.  Each family is shaped so that ONE regime varies and
# everything else is pinned out of the way, and every point's placement is
# re-asserted from the same cancellation the three laws see.
# ---------------------------------------------------------------------------

def e_axis_case(e_target, end_group, regime):
    """Stage-1 E axis.  Shape (x, 1e9 x, 5e10 x): mu0 is the soft-min by nine
    orders, so E == mu0/FLOOR to the last bit, while q10 = mu1 - mu0 ~ 1e9*mu0
    puts m_dist ~ 1e9*E far above M_hi at every E sampled -- so stage 2 is
    inactive in all three copies (including copy 2, which unlike copies 1 and 3
    does NOT skip stage 2 on end-group rows; see the known divergences below).
    The bundle is non-degenerate: all three cap terms mu_k/b_k equal mu0, so
    the soft-min genuinely binds against s_base rather than passing it
    through."""
    x = e_target * FLOOR
    mu = (x, 1.0e9 * x, 5.0e10 * x)
    s_base = (mu[0] if end_group else mu[1]) / V_POLY
    tag = "eg" if end_group else "neg"
    return Case(label=f"{tag}/E={e_target:.6g}", mu=mu, end_group=end_group,
                s_base=s_base, regime=regime, bar="derived")


def m_axis_case(m_target, regime):
    """Stage-2 M axis: the existing two-way pin's own family (mu0 = 1e6 floors
    so stage 1 is bulk-inactive, b1 = mu2/mu1 = 50 > 1 so every point reaches
    the M-band decision), extended to the third copy."""
    mu0 = 1.0e6 * FLOOR
    mu1 = mu0 + m_target * FLOOR
    mu2 = 50.0 * mu1
    return Case(label=f"M={m_target:.6g}", mu=(mu0, mu1, mu2), end_group=False,
                s_base=mu1 / V_POLY, regime=regime, bar="bitwise")


def b1_surface_case(b1n_target, regime, m_target=40.0):
    """Inside the dead band (M = 40 floors), across the degenerate b1c == 1
    surface where the I-090 divergence-repairing completion runs.  b1_n is the
    position within the noise band b1_band = (rtol*(mu1 + mu2) + f1 + f2)/mu1;
    b1_n >= 1 hands back to the exact hard zero."""
    mu0 = 1.0e6 * FLOOR
    mu1 = mu0 + m_target * FLOOR
    band = max((RTOL * (mu1 + mu1) + 2.0 * FLOOR) / mu1,
               solver_mod.CONE_B1_NOISE_REL_FLOOR)
    mu2 = mu1 * (1.0 + b1n_target * band)
    return Case(label=f"b1_n={b1n_target:.6g}", mu=(mu0, mu1, mu2),
                end_group=False, s_base=mu1 / V_POLY, regime=regime,
                bar="bitwise")


_MU0 = 1.0e6 * FLOOR

E_AXIS = [
    (1.0e6, "E bulk"),
    (E_HI * (1.0 + 1.0e-6), "E bulk, just outside the E_hi edge"),
    (E_HI * (1.0 - 1.0e-6), "E blend, just inside the E_hi edge"),
    (5.0e3, "E blend, mid-band"),
    (1.0e3, "E blend"),
    (E_LO * (1.0 + 1.0e-6), "E blend, just outside the E_lo edge"),
    (E_LO * (1.0 - 1.0e-6), "E tail, just inside the E_lo edge"),
    (1.0e1, "E tail"),
    (1.0e-3, "E tail, deep"),
]

M_AXIS = [
    (1.0e-2, "M dead band, deep"),
    (4.0e1, "M dead band"),
    (M_LO * (1.0 - 1.0e-6), "M dead band, just inside the M_lo edge"),
    (M_LO * (1.0 + 1.0e-6), "M blend, just outside the M_lo edge"),
    (1.0e3, "M blend, mid-band"),
    (M_HI * (1.0 - 1.0e-6), "M blend, just inside the M_hi edge"),
    (M_HI * (1.0 + 1.0e-6), "M bulk, just outside the M_hi edge"),
    (1.0e6, "M bulk"),
]

B1_SURFACE = [1.0e-6, 0.25, 0.5, 0.75, 1.0 - 1.0e-4, 1.0 + 1.0e-4, 2.0, 1.0e1]

EARLY = [
    Case("b1c<1 (mu2 < mu1)", (_MU0, 2 * _MU0, 1.5 * _MU0), False, 2 * _MU0,
         "early return: b1c <= 1"),
    Case("b1c==1 exactly", (_MU0, 2 * _MU0, 2 * _MU0), False, 2 * _MU0,
         "early return: b1c <= 1"),
    Case("q10==0", (_MU0, _MU0, 50 * _MU0), False, _MU0,
         "early return: q10 <= 0"),
    Case("q10<0 (out of cone)", (_MU0, 0.5 * _MU0, 50 * _MU0), False,
         0.5 * _MU0, "early return: q10 <= 0"),
    Case("empty pool", (0.0, 0.0, 0.0), False, 1.0, "empty bundle"),
    Case("mu1==0", (_MU0, 0.0, 0.0), False, 1.0, "empty bundle"),
    Case("near-empty, 1e-31 (< SMALL_EPS)", (1e-31, 1e-31, 1e-31), False,
         1e-31, "empty bundle"),
    Case("eg empty pool", (0.0, 0.0, 0.0), True, 1.0, "empty bundle"),
    Case("eg mu0==0", (0.0, 1e-9, 1e-8), True, 1.0, "empty bundle"),
    Case("both stages in blend", (4e-11, 5e-11, 2.5e-9), False, 5e-11,
         "E blend AND M blend", bar="derived"),
    Case("mu3/mu2 <= SMALL_EPS", (1e-11, 1e-11, 1e-31), False, 1e-11,
         "mu3 closure: mu2 below SMALL_EPS", note="bundle-differs"),
    Case("mu3/Cauchy-Schwarz violation", (1e-11, 3e-11, 4e-11), False, 3e-11,
         "mu3 closure: mu0*mu2 < mu1^2"),
    Case("mu3/ln-overflow guard", (1.0e-12, 1.0e-9, 1.0e97), False, 1.0e-9,
         "mu3 closure: past LN_EXP_OVERFLOW_GUARD"),
]

AGREEING_CASES = (
    [e_axis_case(t, eg, r) for t, r in E_AXIS for eg in (True, False)]
    + [m_axis_case(t, r) for t, r in M_AXIS]
    + [b1_surface_case(b, "dead band, across the b1c == 1 surface")
       for b in B1_SURFACE]
    + EARLY
)


# ---------------------------------------------------------------------------


class TestThreeWayAdapter:
    """The adapter, not the law: are the three copies being asked the same
    question?  An adapter that silently feeds different bundles to different
    copies is exactly how a cross-check comes out green for the wrong reason."""

    def test_preconditions_hold_on_the_live_objects(self):
        tw.assert_preconditions()

    def test_the_third_copy_now_reads_the_solver_band_edges(self):
        """Before this ticket ``_s_eff`` spelled every band edge out as a
        literal, so the existing two-way constants pin could not reach it and a
        moved edge would have left ~30 expected-value assertions silently
        evaluating the old band.  Identity, not equality: the same objects."""
        seff = tw.seff_mod
        assert seff.BUNDLE_LIMITER_E_LO is solver_mod.BUNDLE_LIMITER_E_LO
        assert seff.BUNDLE_LIMITER_E_HI is solver_mod.BUNDLE_LIMITER_E_HI
        assert seff.CONE_MARGIN_M_LO is solver_mod.CONE_MARGIN_M_LO
        assert seff.CONE_MARGIN_M_HI is solver_mod.CONE_MARGIN_M_HI
        assert (seff.BUNDLE_LIMITER_SOFTMIN_P
                is solver_mod.BUNDLE_LIMITER_SOFTMIN_P)
        assert (seff.CONE_B1_NOISE_REL_FLOOR
                is solver_mod.CONE_B1_NOISE_REL_FLOOR)
        assert seff.SMALL_EPS is solver_mod.SMALL_EPS
        assert seff.EXHAUSTION_FLOOR_K is solver_mod.EXHAUSTION_FLOOR_K
        # and the soft-min exponent is the solver's, not a second copy of 8.0
        assert seff._softmin_p.__defaults__[0] is (
            solver_mod.BUNDLE_LIMITER_SOFTMIN_P)

    def test_no_band_edge_literal_survives_in_the_third_copy(self):
        """A regression guard on the edit above: if someone re-introduces a
        literal edge in ``_s_eff``, the constants pin goes back to being
        unable to reach it."""
        import inspect

        src = inspect.getsource(tw.seff_mod._s_eff)
        body = src[src.index('"""', src.index('"""') + 3) + 3:]
        for literal in ("1.0e2", "1.0e4", "1e-30", "100.0 * atol", "8.0",
                        "np.finfo"):
            assert literal not in body, (
                f"{literal!r} is back in _s_eff's body; the third copy has "
                "stopped following the solver's constants")

    @pytest.mark.parametrize(
        "case", AGREEING_CASES, ids=[c.label for c in AGREEING_CASES])
    def test_all_three_are_fed_the_same_per_chain_bundle(self, case):
        """Copies 1 and 2 expose ``_chain_bundle``; copy 3 recomputes its
        bundle inline and exposes nothing, so it is read back out of the term
        list it hands to its own ``_softmin_p``.  Asserted, never assumed."""
        res = evaluate(case)
        assert res.bundle_solver == res.bundle_oracle, (
            f"{case.label}: copies 1 and 2 built different bundles "
            f"{res.bundle_solver} vs {res.bundle_oracle}")
        s_terms, e_terms = solver_cap_terms(res), seff_cap_terms(res)
        if e_terms is None:
            # no stage-1 cap was formed: either stage 1 took the bulk early
            # return, or the bundle is empty and copy 3 bailed at the top
            assert res.e_dist >= E_HI or res.bundle_solver[0] <= 0.0
            assert (res.bundle_solver[0] > 0.0) == (s_terms is not None)
            return
        if case.note == "bundle-differs":
            # documented structural difference, characterized below in
            # TestKnownThreeWayDivergences -- assert it POSITIVELY so that
            # closing it has to come past this file
            assert len(s_terms) != len(e_terms)
            return
        assert s_terms is not None and len(s_terms) == len(e_terms), (
            f"{case.label}: cap term counts differ, "
            f"{s_terms} vs {e_terms}")
        for k, (a, b) in enumerate(zip(s_terms, e_terms)):
            if k == 3:
                # the mu2 cap term is built from the mu3 closure, which copies
                # 1/2 evaluate through log/exp and copy 3 through direct
                # powers.  Same closure, different rounding, and NOT a fixed
                # number of ULP -- see mu3_route_relerr's derivation.
                tol = tw.mu3_route_relerr(*case.mu) * max(abs(a), abs(b))
                assert abs(a - b) <= tol, (
                    f"{case.label}: cap term 3: {a!r} vs {b!r}, "
                    f"|d|={abs(a - b):.3e} > derived {tol:.3e}")
                continue
            assert a == b, f"{case.label}: cap term {k}: {a!r} != {b!r}"


class TestThreeWaySweep:
    """The law itself, across every regime the sweep identified."""

    @pytest.mark.parametrize(
        "case", AGREEING_CASES, ids=[c.label for c in AGREEING_CASES])
    def test_three_copies_agree(self, case):
        res = evaluate(case)
        print(res.line())

        if case.bar == "bitwise":
            assert res.solver == res.oracle == res.seff, (
                f"{case.label} ({case.regime}): the three copies evaluate the "
                "identical operations in the identical order here, so this is "
                f"a divergence, not rounding.  solver={res.solver!r} "
                f"oracle={res.oracle!r} seff={res.seff!r}")
            return

        bound = derived_abs_bound(res)
        for name, other in (("oracle", res.oracle), ("seff", res.seff)):
            assert abs(res.solver - other) <= bound, (
                f"{case.label} ({case.regime}): solver vs {name} differ by "
                f"{abs(res.solver - other):.6e}, derived bound "
                f"{bound:.6e} (= {bound / max(U, abs(res.solver) * U):.2f} "
                "ULP of the result)")

    @pytest.mark.parametrize("e_target,regime", E_AXIS,
                             ids=[f"E={c[0]:g}" for c in E_AXIS])
    @pytest.mark.parametrize("end_group", [True, False], ids=["eg", "neg"])
    def test_e_axis_points_land_where_claimed(self, e_target, regime,
                                              end_group):
        """A comparator that samples one regime certifies one regime.  Each
        point's placement is re-derived from the state, never from the nominal
        target -- and the E-bulk points must be an EXACT s_base passthrough,
        which proves the family brackets the band rather than sitting in one
        corner of it."""
        res = evaluate(e_axis_case(e_target, end_group, regime))
        if "bulk" in regime:
            assert res.e_dist >= E_HI
            assert res.solver == res.case.s_base
        elif "blend" in regime:
            assert E_LO < res.e_dist < E_HI
            assert 0.0 < res.solver < res.case.s_base
        else:
            assert res.e_dist <= E_LO
            assert 0.0 < res.solver < res.case.s_base
        # stage 2 must be inactive in ALL THREE, or this family is not
        # measuring stage 1
        assert res.m_dist >= M_HI

    @pytest.mark.parametrize("m_target,regime", M_AXIS,
                             ids=[f"M={c[0]:g}" for c in M_AXIS])
    def test_m_axis_points_land_where_claimed(self, m_target, regime):
        res = evaluate(m_axis_case(m_target, regime))
        assert res.e_dist >= E_HI                 # stage 1 bulk-inactive
        if "dead" in regime:
            assert res.m_dist <= M_LO
            assert res.solver == 0.0
            # ...and the zero is a CHOICE, not a degenerate state
            s_cone = res.q10 / (V_POLY * (res.b1c - 1.0))
            assert softmin_p([res.case.s_base, s_cone]) > 0.0
            assert res.b1_n >= 1.0                # off the b1c == 1 surface
        elif "blend" in regime:
            assert M_LO < res.m_dist < M_HI
            assert 0.0 < res.solver < res.case.s_base
        else:
            assert res.m_dist >= M_HI
            assert res.solver == res.case.s_base

    @pytest.mark.parametrize("b1n", B1_SURFACE, ids=[f"b1n={b:g}"
                                                    for b in B1_SURFACE])
    def test_b1_surface_points_land_where_claimed(self, b1n):
        res = evaluate(b1_surface_case(
            b1n, "dead band, across the b1c == 1 surface"))
        assert res.e_dist >= E_HI
        assert res.m_dist <= M_LO                 # inside the dead band
        # placement re-derived from the state, not taken from the target
        assert res.b1_n == pytest.approx(b1n, rel=1e-5, abs=1e-12)
        assert (res.b1_n >= 1.0) == (b1n >= 1.0)
        if b1n >= 1.0:
            assert res.solver == 0.0              # bit-for-bit, not merely tiny
        else:
            # the completion runs, and hands back CONTINUOUSLY to that zero
            assert res.solver > 0.0
            assert res.solver <= res.case.s_base * (1.0 + 1e-9)


class TestThreeWayFuzz:
    """The enumerated families certify the regimes someone thought of.  This
    one samples the domain at random (fixed seed, so a failure is
    reproducible) and reports the worst residual per regime, which is what
    catches a regime nobody enumerated."""

    def test_random_states_agree_outside_the_known_divergences(self):
        rng = np.random.default_rng(20260905)
        worst = {}
        n_kept = 0
        for _ in range(4000):
            mu = tuple(10.0 ** rng.uniform(-24.0, -4.0, 3))
            if rng.random() < 0.5:                # bias toward ordered states
                mu = tuple(sorted(mu))
            end_group = bool(rng.random() < 0.5)
            s_base = mu[0] if end_group else mu[1]
            case = Case("fuzz", mu, end_group, s_base, "fuzz", bar="derived")
            if _is_known_divergence(case):
                continue
            res = evaluate(case)
            if not (np.isfinite(res.solver) and np.isfinite(res.oracle)
                    and np.isfinite(res.seff)):
                continue
            n_kept += 1
            key = _classify(res)
            bound = derived_abs_bound(res)
            err = max(abs(res.solver - res.oracle),
                      abs(res.solver - res.seff))
            prev = worst.get(key)
            if prev is None or err / max(bound, 5e-324) > prev[0]:
                worst[key] = (err / max(bound, 5e-324), err, bound, res)

        print(f"\nfuzz: {n_kept} states kept, worst |residual| per regime "
              f"(as a fraction of that point's derived bound):")
        for key in sorted(worst):
            frac, err, bound, res = worst[key]
            print(f"  {key:<26s} err={err:.6e} bound={bound:.6e} "
                  f"frac={frac:.4f}   mu={res.case.mu}")
        assert n_kept > 500, f"fuzz kept too few states ({n_kept})"
        for key, (frac, err, bound, res) in worst.items():
            assert err <= bound, (
                f"fuzz regime {key}: solver/oracle/seff differ by {err:.6e} "
                f"against a derived bound of {bound:.6e} at mu={res.case.mu}, "
                f"end_group={res.case.end_group}")


def _is_known_divergence(case):
    """The three regimes where the copies genuinely disagree, characterized in
    TestKnownThreeWayDivergences below.  Excluded here so the fuzz measures
    rounding rather than re-finding what is already on the record."""
    mu0, mu1, mu2 = case.mu
    e_dist = softmin_p([mu0, mu1, mu2]) / FLOOR
    if case.end_group:
        # D-A: copy 2 does not skip stage 2 on end-group rows
        return mu0 > 0.0 and mu1 > mu0 and 0.0 < (mu1 - mu0) < M_HI * FLOOR
    # D-B/D-C: the mu3 realizability guard, which only copies 1 and 2 have,
    # and which they implement differently from each other on the noise band
    return mu1 < mu0 and e_dist < E_HI


def _classify(res):
    c = res.case
    stage1 = ("E-bulk" if res.e_dist >= E_HI else
              "E-blend" if res.e_dist > E_LO else "E-tail")
    if c.end_group:
        return f"eg/{stage1}"
    if not np.isfinite(res.b1c) or res.b1c <= 1.0:
        return f"neg/{stage1}/b1c<=1"
    if res.q10 <= 0.0:
        return f"neg/{stage1}/q10<=0"
    stage2 = ("M-bulk" if res.m_dist >= M_HI else
              "M-blend" if res.m_dist > M_LO else "M-dead")
    return f"neg/{stage1}/{stage2}"


class TestKnownThreeWayDivergences:
    """Characterization, NOT approval.

    The sweep found three points where the three copies do not implement one
    law.  Repairing them is out of this ticket's scope -- deciding which of
    three implementations is correct is a separate call, and making it
    silently would destroy the measurement.  They are pinned here with their
    numbers so that (a) they are on the record and (b) whoever resolves them
    has to come past this file and update it deliberately.
    """

    def test_D_A_oracle_does_not_skip_stage_2_on_end_group_rows(self):
        """Copy 1 (``if end_group: return s_free``, the round-62 N5b
        adjudicated fix) and copy 3 skip the stage-2 cone gate entirely for
        end-group rows.  Copy 2 does not: it forms ``b1c = mu1/mu0`` and falls
        through into the M-band decision, which is the PRE-round-62 behaviour.

        Reachability is not exotic -- any end-group row on a pool with mean
        chain length > 1 whose cone margin mu1 - mu0 sits below M_hi floors.
        The existing two-way pin cannot see it: it is non-end-group only.
        """
        a = 1.0e3 * FLOOR / softmin_p([1.0, 4.0, 20.0])
        blend = Case("D-A/blend", (a, 4 * a, 20 * a), True, a / V_POLY, "D-A")
        res = evaluate(blend)
        assert res.e_dist == pytest.approx(1.0e3, rel=1e-12)
        # copy 2's stage-2 gate is ACTIVE here (its M blend); copies 1 and 3
        # never reach it at all on an end-group row
        assert M_LO < res.m_dist < M_HI
        assert res.b1c == pytest.approx(4.0, rel=1e-12)   # = mu1/mu0 > 1
        assert res.solver == res.seff       # copies 1 and 3 are one law
        assert res.oracle < res.solver
        assert res.rel_12 == pytest.approx(2.2497e-2, rel=1e-3), res.rel_12

        b = 1.0e1 * FLOOR / softmin_p([1.0, 4.0, 20.0])
        tail = Case("D-A/tail", (b, 4 * b, 20 * b), True, b / V_POLY, "D-A")
        res = evaluate(tail)
        assert res.m_dist <= M_LO           # copy 2's dead band: hard zero
        assert res.solver == pytest.approx(8.408980191297415e-14, rel=1e-12)
        assert res.seff == pytest.approx(res.solver, rel=1e-15)
        assert res.oracle == 0.0, (
            "copy 2 zeroes an end-group drain that copies 1 and 3 keep")

    def test_D_B_mu3_noise_band_splits_all_three(self):
        """``_safe_mu3_from_mu012`` (copy 1) returns mu1 for a realizability
        violation inside MU3_CLOSURE_BOUNDARY_REL of the boundary; ``safe_mu3``
        (copy 2) returns 0.0 for ANY mu1 < mu0 -- it never got the N5b
        boundary-consistent band; ``_s_eff`` (copy 3) has no realizability
        guard at all and computes mu0*(mu2/mu1)**3 regardless.  Three laws."""
        mu0 = 1.0e-11
        mu1 = mu0 * (1.0 - 1.0e-9)          # violation of 1e-9 rel, inside the
        assert 0 < (mu0 - mu1) <= solver_mod.MU3_CLOSURE_BOUNDARY_REL * mu0
        case = Case("D-B", (mu0, mu1, 0.5e-11), False, mu0, "D-B")
        res = evaluate(case)
        assert E_LO < res.e_dist < E_HI     # stage 1 active, so mu3 matters
        assert res.b1c <= 1.0               # stage 2 returns s_free untouched
        # copy 1: mu3 -> mu1, so b2 = mu1/mu1 = 1
        assert res.bundle_solver[2] == pytest.approx(1.0, rel=1e-8)
        # copy 2: mu3 -> 0, so b2 == 0 and the mu2 cap term is dropped entirely
        assert res.bundle_oracle[2] == 0.0
        assert res.rel_12 == pytest.approx(0.4528, rel=1e-3), res.rel_12
        assert res.rel_23 == pytest.approx(9.4688e-7, rel=1e-3), res.rel_23

    def test_D_C_copy_3_has_no_realizability_guard_at_all(self):
        """A DEEP cone violation (mu1 = mu0/2): copies 1 and 2 agree exactly
        (both zero the closure), copy 3 extrapolates through it."""
        case = Case("D-C", (1e-11, 0.5e-11, 0.2e-11), False, 1e-11, "D-C")
        res = evaluate(case)
        assert res.solver == res.oracle
        assert res.rel_13 == pytest.approx(1.6108e-3, rel=1e-3), res.rel_13

    def test_the_ln_overflow_window_is_structurally_unreachable(self):
        """Copy 3 has no LN_EXP_OVERFLOW_GUARD, so on paper it disagrees with
        copies 1/2 for mu3 in (exp(700), DBL_MAX).  It cannot be reached
        through this law: the mu3 closure runs only on non-end-group rows,
        which need E < E_hi, i.e. min(mu) < E_hi*FLOOR = 1e-10 mol; with
        mu1 > mu0 that puts ln(mu0) < -23, so the window needs
        (mu2/mu1)**3 > exp(723) = 1.8e314, which overflows to inf in copy 3's
        own expression before mu0 ever multiplies it -- and an inf mu3 drops
        the b2 term in all three.  Demonstrated, not asserted from the
        algebra."""
        case = Case("ln-overflow", (1.0e-12, 1.0e-9, 1.0e97), False, 1.0e-9,
                    "ln overflow")
        res = evaluate(case)
        ln_mu3 = (3.0 * np.log(1.0e97) - 3.0 * np.log(1.0e-9)
                  + np.log(1.0e-12))
        assert ln_mu3 > solver_mod.LN_EXP_OVERFLOW_GUARD    # in the window
        assert solver_mod._safe_mu3_from_mu012(1.0e-12, 1.0e-9,
                                               1.0e97) == float("inf")
        with np.errstate(over="ignore"):
            assert 1.0e-12 * (np.float64(1.0e97) / 1.0e-9) ** 3 == float("inf")
        assert res.solver == res.oracle == res.seff

    def test_mu2_below_small_eps_is_structural_only(self):
        """Copies 1/2 drop the b2 cap term when mu2 <= SMALL_EPS; copy 3 keeps
        a b2 of ~1e-60, whose cap term mu2/b2 is ~1e29.  A term that large
        contributes ~1e-232 to the reciprocal 8-norm, so the soft-min is
        unchanged to the last bit: a real structural difference with no
        numerical consequence.  Recorded so it is not re-found as new."""
        case = Case("mu2<eps", (1e-11, 1e-11, 1e-31), False, 1e-11, "structural")
        res = evaluate(case)
        assert res.bundle_solver[2] == 0.0 and res.bundle_oracle[2] == 0.0
        terms = seff_cap_terms(res)
        assert terms is not None and len(terms) == 4      # copy 3 keeps it
        assert terms[3] > 1.0e20
        assert res.solver == res.oracle == res.seff       # ...and it vanishes
