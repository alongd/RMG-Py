"""
A kernel whose flows are linear in its pool's own moments must keep one
slope through zero.

On PP, DASPK's step is pinned at ~1e-5 s by nonlinear-convergence failures,
with no error-test failures.  The culprits are near-empty mid-chain-radical
pools drained by moves at k ~ 1e7 s^-1.  The BDF predictor extrapolates their
decay below zero.  There the pool loop skipped every kernel (``mu0 > 0``, all
flows >= 0) and read the moments through the softclamp, so the iteration
matrix was built with the drain off, while at the solution the drain is on.
With h*k in the hundreds, the corrector cannot converge.

A move with a = 0, an interior scission, an unsited transfer and a partnered
recombination are each linear in the pool's own moments and need no closure.
Their exact flows are odd in the moments, so a negative excursion is drained
back, at the same rate constant as a positive one.  Kernels on the mu3
closure keep their guards.
"""
import numpy as np
import pytest

from rmgpy.solver.polymer import HybridPolymerSystem
from polymerSubAtolPoolTest import (PP_RA, PP_RM, _KickedPairSystem, _index, _moments, _simulate,
                                    pair_deck)

# a near-empty pool's moments, PP's noise sizes (mol)
NOISE = (1.0e-16, 1.5e-15, 1.3e-14)


def _rows(rs, labels):
    return [int(i) for label in labels for i in rs.polymer_pools[_index(rs, label)].mu_indices]


def _rates(rs, y):
    return np.asarray(rs.residual(0.0, y, np.zeros(rs.neq), rs.senpar)[0], dtype=float)


def _odd_part(rs, label, rows, scale=1.0, moments=NOISE):
    """The residual's change on ``rows`` when the pool's moments go from 0 to
    +scale*moments and to -scale*moments (every other slot at y0)."""
    y0 = np.array(rs.y0, dtype=float)
    base = _rates(rs, y0)
    out = []
    for sign in (1.0, -1.0):
        y = y0.copy()
        for i, m in zip(rs.polymer_pools[_index(rs, label)].mu_indices, moments):
            y[int(i)] = sign * scale * m
        out.append((_rates(rs, y) - base)[rows])
    return out


class TestOneSlopeThroughZero:

    @pytest.mark.parametrize("label, kernels", [("Rm", "move"), ("Ra", "move + interior scission"),
                                                ("Re", "unsited transfer")])
    def test_a_linear_kernel_is_odd_in_its_pools_moments(self, label, kernels):
        """Every row the pool's linear kernels write changes sign, not size,
        when the pool's moments do."""
        rs = pair_deck()
        rows = _rows(rs, [p.label for p in rs.polymer_pools])
        up, down = _odd_part(rs, label, rows)
        assert np.max(np.abs(up)) > 0.0
        np.testing.assert_allclose(down, -up, rtol=1e-12, atol=0.0)

    def test_the_drain_slope_is_the_rate_constant_on_both_sides_of_zero(self):
        """Rm's own rows: d(dmu_k/dt)/dmu_k = -k_move at +/- half an error weight."""
        rs = pair_deck()
        rows = _rows(rs, ["Rm"])
        for scale in (0.5, 1.0):
            up, down = _odd_part(rs, "Rm", rows, scale)
            np.testing.assert_allclose(up, -3.7e7 * scale * np.array(NOISE), rtol=1e-12)
            np.testing.assert_allclose(down, 3.7e7 * scale * np.array(NOISE), rtol=1e-12)

    def test_noise_without_a_distributions_shape_still_drains_linearly(self):
        """Inside the noise scale the fragment law is not judged on the
        moments' shape: a mean below one (Re's mu1 credit -m0 + m1/2 < 0)
        still breaks and drains, odd in the moments.  A filled pool outside
        the law's domain is still skipped (polymerInteriorScissionTest)."""
        rs = pair_deck()
        rows = _rows(rs, [p.label for p in rs.polymer_pools])
        up, down = _odd_part(rs, "Ra", rows, moments=(1.0e-16, 0.5e-16, 1.0e-17))
        assert abs(up[rows.index(_rows(rs, ["Re"])[1])]) > 0.0
        np.testing.assert_allclose(down, -up, rtol=1e-12, atol=0.0)


class TestNegativeExcursion:

    def test_a_negative_excursion_is_pulled_back(self):
        """The pair kicked to PP's refused (negative) state drains back toward
        zero at the moves' and the scission's rates (held at the kick before);
        the feed then refills the pair to a few thousandths of an error weight."""
        rs = pair_deck(cls=_KickedPairSystem)
        _simulate(rs, 1.0e-3)
        for label, kicked in (("Rm", PP_RM), ("Ra", PP_RA)):
            for m, m_kick in zip(_moments(rs, label), kicked):
                assert abs(m) < 1.0e-2 * abs(m_kick)


class _KickedPositivePairSystem(_KickedPairSystem):
    """The pair kicked to the size of PP's refused state, positive, then left
    to drain across zero (the fast drain's predictor overshoots)."""

    def residual(self, t, y, dydt, senpar=None):
        if senpar is None:
            delta, ires = HybridPolymerSystem.residual(self, t, y, dydt)
        else:
            delta, ires = HybridPolymerSystem.residual(self, t, y, dydt, senpar)
        if t < self.T_KICK:
            for label, target in (("Rm", PP_RM), ("Ra", PP_RA)):
                for i, m in zip(self.polymer_pools[_index(self, label)].mu_indices, target):
                    delta[i] += abs(m) / self.T_KICK
        return delta, ires

    def step(self, tout):
        self.accepted_steps = getattr(self, "accepted_steps", 0) + 1
        return super().step(tout)


class TestCorrectorAcrossZero:

    def test_a_fast_drain_crossing_zero_does_not_pin_the_step(self):
        """With the drain switched off below zero, this run took ~5900
        accepted steps to 1e3 s (2877 nonlinear-convergence failures); the
        same deck without the kick takes ~100."""
        rs = pair_deck(cls=_KickedPositivePairSystem)
        _simulate(rs, 1.0e3)
        assert rs.t >= 1.0e3
        assert rs.accepted_steps < 600


class _SinkSystem(HybridPolymerSystem):
    """A sign error: a constant sink of RATE mol/s in each moment of POOL."""

    POOL, RATE = "Ut", 1.0e-9

    def residual(self, t, y, dydt, senpar=None):
        if senpar is None:
            delta, ires = HybridPolymerSystem.residual(self, t, y, dydt)
        else:
            delta, ires = HybridPolymerSystem.residual(self, t, y, dydt, senpar)
        for i in self.polymer_pools[_index(self, self.POOL)].mu_indices:
            delta[i] -= self.RATE
        return delta, ires


class TestSignErrorStillRefused:

    def test_a_sink_on_a_pool_without_a_drain_is_refused(self):
        rs = pair_deck(cls=_SinkSystem)
        with pytest.raises(ValueError, match="beyond the exhaustion floor"):
            _simulate(rs, 1.0)
        assert rs.t < 1.0e-3

    def test_a_sink_the_drain_cannot_restore_is_refused(self):
        """On Ra the linear drains hold a sink at about -RATE/k; one past
        k times the floor (k ~ 2.5e6 s^-1 for the pair, floor 1e-14 mol) is
        refused."""
        rs = pair_deck(cls=type("_SinkRa", (_SinkSystem,), {"POOL": "Ra", "RATE": 1.0e-6}))
        with pytest.raises(ValueError, match="beyond the exhaustion floor"):
            _simulate(rs, 1.0)
        assert rs.t < 1.0e-3
