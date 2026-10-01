"""One driver, three implementations of the near-exhaustion bundle-limiter law.

The law (`_bundle_limited_site`: stage-1 exhaustion tail limiter on the E band,
stage-2 cone-margin drain gate on the M band) is written out THREE times in
this repository, on purpose -- three independent implementations are what make
a cross-check mean anything:

  1. the compiled solver
     ``rmgpy/solver/polymer.pyx``  :: ``HybridPolymerSystem._bundle_limited_site``
     signature ``(pool_idx, y, V_poly, end_group, s_base)``
  2. the rmgpy-free numpy oracle
     ``test/rmgpy/tools/numpy_moments_consumer.py`` :: ``ArtifactConsumer._bundle_limited_site``
     signature ``(pool_label, y, end_group, s_base)``
  3. the test-side expected-value generator
     ``test/rmgpy/solver/solverPolymerTest.py`` :: ``_s_eff``
     signature ``(mu, end_group, s_base, v_poly, atol, rtol)``

Copy 3 is the dangerous one: ~30 assertions in solverPolymerTest.py compare
solver output against ``kf * _s_eff(...)``.  A divergence in copy 2 makes a test
fail, which announces itself.  A divergence in copy 3 rewrites the standard and
every test keeps passing against a number nobody checked.

WHY THIS MODULE EXISTS HERE
---------------------------
Neither ``test/rmgpy/solver/`` nor ``test/rmgpy/tools/`` is a Python package
(no ``__init__.py``), so neither copy 2 nor copy 3 is importable by a dotted
name.  Existing tests reach copy 2 with a ``sys.path.insert`` of their own
directory.  This module lives beside the oracle and the existing two-way pin
(``TestConeMarginBandParity`` in ``polymerMomentsConsumerTest.py``, whose deck
it reuses verbatim) and reaches copy 3 by inserting the SIBLING solver test
directory, computed from ``__file__``.  Nothing here depends on pytest
collection order or on any other test module having been imported first, so a
single node run standalone works exactly like a full-suite run.

The module deliberately does NOT re-implement the law.  Everything it computes
independently (``e_dist``, ``m_dist``, ``b1_band``) is placement arithmetic used
to assert that a sample point lands in the regime it claims, from the same
cancellation the three laws see -- never from the nominal target.
"""

import os
import sys
from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np

import rmgpy.solver.polymer as solver_mod

_TOOLS_DIR = os.path.dirname(os.path.abspath(__file__))
_SOLVER_TEST_DIR = os.path.join(os.path.dirname(_TOOLS_DIR), "solver")
for _p in (_TOOLS_DIR, _SOLVER_TEST_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy_moments_consumer as consumer_mod  # noqa: E402
import polymerMomentsConsumerTest as parity_mod  # noqa: E402
import solverPolymerTest as seff_mod  # noqa: E402

from numpy_moments_consumer import softmin_p  # noqa: E402

# ---------------------------------------------------------------------------
# Canonical case description.  ONE description -> three call forms.
# ---------------------------------------------------------------------------

V_POLY = parity_mod.V_POLY            # 1.0 m^3 of polymer phase
ATOL = 1.0e-16                        # the deck's scalar absolute tolerance
RTOL = 1.0e-8                         # the deck's scalar relative tolerance
FLOOR = max(solver_mod.SMALL_EPS, solver_mod.EXHAUSTION_FLOOR_K * ATOL)  # 1e-14

U = 2.0 ** -53                        # IEEE-754 binary64 unit roundoff


@dataclass(frozen=True)
class Case:
    """A point on the law's domain, in the units the SOLVER uses (moles)."""

    label: str
    mu: Tuple[float, float, float]    # (mu0, mu1, mu2) [mol]
    end_group: bool
    s_base: float                     # [mol/m^3_poly]
    regime: str                       # what the sweep claims about this point
    bar: str = "bitwise"              # "bitwise" | "derived"
    note: str = ""


@dataclass
class Result:
    case: Case
    solver: float = float("nan")
    oracle: float = float("nan")
    seff: float = float("nan")
    e_dist: float = float("nan")      # placement, computed independently
    m_dist: float = float("nan")
    b1c: float = float("nan")
    q10: float = float("nan")
    b1_band: float = float("nan")
    b1_n: float = float("nan")
    bundle_solver: tuple = ()
    bundle_oracle: tuple = ()
    seff_softmin_calls: list = field(default_factory=list)

    # -- residuals ---------------------------------------------------------
    def _rel(self, a, b):
        if a == b:
            return 0.0
        scale = max(abs(a), abs(b))
        if scale == 0.0 or not np.isfinite(scale):
            return float("inf")
        return abs(a - b) / scale

    @property
    def rel_12(self):
        return self._rel(self.solver, self.oracle)

    @property
    def rel_13(self):
        return self._rel(self.solver, self.seff)

    @property
    def rel_23(self):
        return self._rel(self.oracle, self.seff)

    @property
    def max_rel(self):
        return max(self.rel_12, self.rel_13, self.rel_23)

    @property
    def all_bitwise(self):
        return self.solver == self.oracle == self.seff

    def line(self):
        return (
            f"{self.case.label:<34s} "
            f"E={self.e_dist:11.4e} M={self.m_dist:11.4e} "
            f"b1c-1={self.b1c - 1.0:+.6e} b1_n={self.b1_n:9.3e} | "
            f"solver={self.solver!r:>26s} oracle={self.oracle!r:>26s} "
            f"seff={self.seff!r:>26s} | "
            f"r12={self.rel_12:.3e} r13={self.rel_13:.3e} r23={self.rel_23:.3e}"
        )


# ---------------------------------------------------------------------------
# The deck.  Built once; the SAME single-pool deck the existing two-way pin
# uses, so the third copy is being asked the question the first two already
# answer there.
# ---------------------------------------------------------------------------

_PAIR = None


def triplet():
    """(rs, consumer) for the canonical single-pool deck (cached)."""
    global _PAIR
    if _PAIR is None:
        _PAIR = parity_mod._cone_band_pair()
    return _PAIR


def assert_preconditions():
    """The adapter's load-bearing part: the three copies must be reading the
    SAME numbers out of the same state before any comparison of their outputs
    means anything.  Asserted from the live objects, never assumed."""
    rs, consumer = triplet()

    # 1. per-moment floors (copy 1) vs the single scalar floor (copies 2, 3).
    #    These coincide ONLY on a uniform-atol deck; assert the uniformity
    #    rather than relying on it.
    floors = np.asarray(rs._pool_mu_floors)
    assert floors.shape == (1, 3), floors.shape
    assert np.all(floors == FLOOR), floors
    assert consumer.mu_floor == FLOOR
    assert max(1e-30, 100.0 * ATOL) == FLOOR          # copy 3's own expression
    # copy 1 forms the stage-2 floor as max(f0, f1); copies 2/3 use the scalar
    assert max(floors[0, 0], floors[0, 1]) == FLOOR

    # 2. moment slots
    i012 = tuple(rs.polymer_pools[0].mu_indices)
    assert i012 == tuple(consumer.pools["poly"]["mu"]), (i012, consumer.pools)

    # 3. the b1-noise band's rtol half
    assert rs._cone_b1_rtol == RTOL
    assert consumer.cone_b1_rtol == RTOL

    # 4. volume basis
    assert consumer.V_poly == V_POLY
    assert V_POLY == 1.0        # copy 3 takes moles and divides by v_poly

    # 5. softclamp vs max(0, .): identical on non-negative states, which is
    #    what every sample point uses (asserted per case in evaluate()).
    assert rs._softclamp_lam > 0.0

    # 6. band constants: all three now read ONE set of module-level names
    assert consumer_mod.BUNDLE_LIMITER_E_LO == solver_mod.BUNDLE_LIMITER_E_LO
    assert consumer_mod.BUNDLE_LIMITER_E_HI == solver_mod.BUNDLE_LIMITER_E_HI
    assert consumer_mod.CONE_MARGIN_M_LO == solver_mod.CONE_MARGIN_M_LO
    assert consumer_mod.CONE_MARGIN_M_HI == solver_mod.CONE_MARGIN_M_HI
    assert (consumer_mod.BUNDLE_LIMITER_SOFTMIN_P
            == solver_mod.BUNDLE_LIMITER_SOFTMIN_P)
    assert (consumer_mod.CONE_B1_NOISE_REL_FLOOR
            == solver_mod.CONE_B1_NOISE_REL_FLOOR)
    return rs, consumer


class _SoftminRecorder:
    """Wraps copy 3's module-level ``_softmin_p`` so the term lists it builds
    are observable.  This is how the comparator proves copy 3 is fed the same
    per-chain bundle as copies 1 and 2 -- copy 3 recomputes its bundle inline
    and exposes nothing, so an adapter that merely *assumed* they matched would
    be the classic 'agrees for the wrong reason' failure."""

    def __init__(self):
        self.calls = []

    def __enter__(self):
        self._orig = seff_mod._softmin_p

        def spy(terms, p=self._orig.__defaults__[0]):
            terms = list(terms)
            out = self._orig(terms, p)
            self.calls.append((terms, p, out))
            return out

        seff_mod._softmin_p = spy
        return self

    def __exit__(self, *exc):
        seff_mod._softmin_p = self._orig
        return False


def evaluate(case: Case) -> Result:
    """Drive all three copies from ONE case description."""
    rs, consumer = triplet()
    assert all(m >= 0.0 for m in case.mu), (
        "sample points are non-negative so copy 1's softclamp and copies 2/3's "
        "max(0, .) are the identity on them -- otherwise the three would be "
        f"reading different numbers: {case.mu}")

    i0, i1, i2 = tuple(rs.polymer_pools[0].mu_indices)
    y = np.asarray(rs.y, dtype=np.float64).copy()
    y[i0], y[i1], y[i2] = case.mu
    mu0, mu1, mu2 = case.mu

    res = Result(case=case)

    # -- placement arithmetic (independent of all three laws) --------------
    res.e_dist = softmin_p([mu0, mu1, mu2]) / FLOOR
    res.q10 = mu1 - mu0
    res.m_dist = res.q10 / FLOOR
    # b1c as each stage-2 branch forms it.  NB copies 1 and 3 skip stage 2 for
    # end-group rows entirely (round-62 N5b); copy 2 does not -- the end-group
    # b1c below is the one copy 2 alone still forms.
    if case.end_group:
        res.b1c = (mu1 / mu0) if mu0 > 0.0 else float("nan")
    else:
        res.b1c = (mu2 / mu1) if mu1 > 0.0 else float("nan")
    if mu1 > 0.0:
        band = (RTOL * (mu1 + mu2) + FLOOR + FLOOR) / mu1
        res.b1_band = max(band, solver_mod.CONE_B1_NOISE_REL_FLOOR)
        if np.isfinite(res.b1c):
            res.b1_n = (res.b1c - 1.0) / res.b1_band

    # -- bundles: the adapter's "same question" assertion -------------------
    res.bundle_solver = tuple(
        rs._chain_bundle(0, y, V_POLY, case.end_group))
    res.bundle_oracle = tuple(
        consumer._chain_bundle("poly", y, case.end_group))

    # -- the three laws -----------------------------------------------------
    res.solver = rs._bundle_limited_site(0, y, V_POLY, case.end_group,
                                         case.s_base)
    res.oracle = consumer._bundle_limited_site("poly", y, case.end_group,
                                               case.s_base)
    with _SoftminRecorder() as rec:
        res.seff = seff_mod._s_eff(case.mu, end_group=case.end_group,
                                   s_base=case.s_base, v_poly=V_POLY,
                                   atol=ATOL, rtol=RTOL)
    res.seff_softmin_calls = rec.calls
    return res


def seff_cap_terms(res: Result):
    """Copy 3's stage-1 cap term list, or None if no stage-1 cap was ever
    formed (stage 1 took the bulk early return, or copy 3 bailed on an empty
    bundle).  Call ordering inside ``_s_eff`` is: #0 = the E-band floor
    distance, #1 = the stage-1 cap (only when E < E_hi AND the bundle is
    non-empty), #2 = the stage-2 softmin_p(s_free, s_cone) (only when the M
    band is entered) -- so #1 must be gated on E, not taken positionally."""
    if (res.e_dist >= solver_mod.BUNDLE_LIMITER_E_HI
            or len(res.seff_softmin_calls) < 2):
        return None
    return res.seff_softmin_calls[1][0]


def solver_cap_terms(res: Result):
    """Copy 1's stage-1 cap terms, reconstructed from its OWN bundle, in the
    order copy 3 lists them: [s_base, mu0, mu1/b1, mu2/b2] over positive b_k."""
    b0, b1, b2, _ok = res.bundle_solver
    if b0 <= 0.0:
        return None
    mu0, mu1, mu2 = res.case.mu
    terms = [res.case.s_base, mu0 / (V_POLY * b0)]
    if b1 > 0.0:
        terms.append(mu1 / (V_POLY * b1))
    if b2 > 0.0:
        terms.append(mu2 / (V_POLY * b2))
    return terms


# ---------------------------------------------------------------------------
# Agreement bar.  DERIVED from the floating-point structure, per regime.
# ---------------------------------------------------------------------------

def softmin_relerr(k, p=None):
    """Bound on the relative error of ``softmin_p`` over k terms, evaluated in
    the shared factored form ``m*(sum (m/x_i)**p)**(-1/p)``:

        m/x_i            1 rounding                        -> 1.0*U
        (.)**p           amplifies relative error by p,
                         plus the pow's own rounding       -> p*1.0*U + 0.5*U
        sum of k terms   k-1 additions                     -> +(k-1)*U
        (.)**(-1/p)      SHRINKS relative error by p,
                         plus the pow's own rounding       -> /p + 0.5*U
        m*(.)            1 rounding                        -> +0.5*U

    m itself is a copy of an input, exact."""
    if p is None:
        p = solver_mod.BUNDLE_LIMITER_SOFTMIN_P
    ratio = 1.0 * U
    powed = p * ratio + 1.0 * U
    summed = powed + (k - 1) * U
    rooted = summed / p + 1.0 * U
    return rooted + 1.0 * U


def mu3_route_relerr(mu0, mu1, mu2):
    """The two routes to the SAME closure mu3 = mu0*(mu2/mu1)**3.

    Copies 1 and 2 evaluate it as ``exp(3*log(mu2) - 3*log(mu1) + log(mu0))``
    (the solver's ``_safe_mu3_from_mu012`` / the oracle's ``safe_mu3``); copy 3
    evaluates the powers directly.  Algebraically identical -- numerically not,
    and not by a fixed number of ULP either: the relative error of ``exp(s)``
    is the ABSOLUTE error of ``s``, and ``s`` is built from logarithms whose
    magnitudes are ~46 at the moment scales this law runs at (1e-20 mol), so
    the log/exp route carries a few hundred U rather than a few.  This
    propagates into the stage-1 cap whenever the mu2 term binds the soft-min.

    Bound: |ds| <= 2U*(3|ln mu2| + 3|ln mu1| + |ln mu0|)   [the logs and their
    scalings] + 3U*|s| [the two adds and the final rounding], and exp adds its
    own U.  Copy 3's own route contributes ~6U, folded in."""
    if not (mu0 > 0.0 and mu1 > 0.0 and mu2 > 0.0):
        return 0.0
    l0, l1, l2 = np.log(mu0), np.log(mu1), np.log(mu2)
    s = 3.0 * l2 - 3.0 * l1 + l0
    mags = 3.0 * abs(l2) + 3.0 * abs(l1) + abs(l0)
    return float(2.0 * U * mags + 3.0 * U * abs(s) + 7.0 * U)


def stage1_cap_relerr(n_bundle_terms):
    """Copies 1/2 fold s_base into the bundle cap as a NESTED softmin --
    softmin_p(s_base, softmin_p(t_0..t_n)) -- while copy 3 evaluates the FLAT
    softmin over [s_base, t_0..t_n].  The reciprocal p-norm is exactly
    associative, so the two agree algebraically and differ only by rounding.
    Bound = nested error + flat error.

    In the nested form the inner result's relative error eps_in enters the
    outer ratio, is amplified by p in ``**p`` and shrunk by p in
    ``**(-1/p)``, i.e. it passes through with unit gain."""
    inner = softmin_relerr(n_bundle_terms)
    p = solver_mod.BUNDLE_LIMITER_SOFTMIN_P
    outer_ratio = 1.0 * U + inner
    outer_powed = p * outer_ratio + 0.5 * U
    outer_summed = outer_powed + 1.0 * U
    outer = outer_summed / p + 0.5 * U + 0.5 * U
    flat = softmin_relerr(n_bundle_terms + 1)
    return outer + flat


def smoothstep_abs_err(dist, lo, hi, dist_relerr):
    """Absolute error of w = n*n*(3 - 2n), n = (dist - lo)/(hi - lo), given a
    relative error ``dist_relerr`` on ``dist``.  |dw/dn| = |6n(1-n)| <= 1.5."""
    d_abs = abs(dist) * dist_relerr + 0.5 * U * abs(dist - lo)
    d_n = d_abs / (hi - lo) + 0.5 * U
    return 1.5 * d_n


def derived_abs_bound(res: Result):
    """Absolute agreement bound for a sample point whose stage 1 is ACTIVE
    (the only regime where the three do the same arithmetic in different
    orders).  Everywhere else the operation order is identical and the bar is
    bitwise -- see the module docstring of the test."""
    case = res.case
    if res.e_dist >= solver_mod.BUNDLE_LIMITER_E_HI:
        # stage 1 took the bulk early return in all three: s_free IS s_base,
        # bit for bit, and every stage-2 branch below it is order-identical.
        return 0.0
    b0, b1, b2, _ok = res.bundle_solver
    n_terms = 1 + (1 if b1 > 0.0 else 0) + (1 if b2 > 0.0 else 0)
    cap_rel = stage1_cap_relerr(n_terms)

    seff_terms = seff_cap_terms(res)
    cap_val = res.seff_softmin_calls[1][2] if seff_terms is not None else 0.0
    bound = cap_rel * abs(cap_val)

    # the mu2 cap term is built from the mu3 closure, whose two evaluation
    # routes differ by far more than the soft-min's own rounding.  softmin_p
    # is monotone with d(softmin)/dx_i = (softmin/x_i)**(p+1) <= 1, so the
    # term's relative error reaches the cap with gain at most 1.
    if not case.end_group and b2 > 0.0:
        bound += mu3_route_relerr(*case.mu) * abs(cap_val)

    if res.e_dist > solver_mod.BUNDLE_LIMITER_E_LO:
        # blend: s_free = w*s_base + (1 - w)*cap.  Copy 1 divides each moment
        # by its own floor then soft-mins; copies 2/3 soft-min then divide by
        # the scalar floor -- algebraically identical, so E can differ by a
        # couple of ULP, and that difference reaches the answer only through w.
        e_rel = softmin_relerr(3) + 1.0 * U
        d_w = smoothstep_abs_err(res.e_dist, solver_mod.BUNDLE_LIMITER_E_LO,
                                 solver_mod.BUNDLE_LIMITER_E_HI, e_rel)
        bound = (bound + d_w * abs(case.s_base - cap_val)
                 + 3.0 * U * abs(res.solver))
    # stage 2, when it is entered at all, is order-identical in all three and
    # propagates s_free's error with gain <= 1: softmin_p is monotone with
    # d(softmin)/dx_i = (softmin/x_i)**(p+1) <= 1, and the v-blend is convex.
    # Only the extra roundings of that path are added here.
    return bound + 4.0 * U * abs(res.solver)
