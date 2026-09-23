#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
Tests for :class:`TerminationSteadyState` -- "integrate until it stops changing".

The numbers embedded here are not invented. They are lifted from two trajectories
measured on 2026-08-30 and recorded in ``docs/i170-steady-state/FINDING.md``:

* the I-164 5 torr argon deck, whose composition is *exactly* constant for all 51
  solver steps; and
* a lithium-seeded argon deck at the same conditions, which does relax to an
  ionisation balance over 497 steps.

The second of those is why the arming latch exists, and
:func:`test_early_quiescence_does_not_fire_the_criterion` is the regression test for it.
"""

import logging

import numpy as np
import pytest

import rmgpy.constants as constants
from rmgpy.kinetics import Arrhenius
from rmgpy.kinetics.arrhenius import TwoTemperaturePlasma
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.solver.simple import SimpleReactor
from rmgpy.solver.termination import TerminationSteadyState, TerminationTime
from rmgpy.species import Species
from rmgpy.thermo import ThermoData

T_GAS = 1000.0    # K
T_E = 11604.5     # K (~1 eV)
P0 = 1.0e5        # Pa

ATOL = 1e-16      # the floor used throughout: the simulator's own absolute tolerance


class TerminationSteadyStateResidualTest:
    """The residual itself: R = max_i |d ln x_i / d ln t|, on hand-built steps."""

    def test_flat_step_has_zero_residual(self):
        """Nothing changed, so nothing is changing: R = 0."""
        y = np.array([1.0, 2.0, 3.0])
        r, _ = TerminationSteadyState.compute_residual(y, np.e, y, 1.0, ATOL)
        assert r == 0.0

    def test_residual_is_the_log_log_slope(self):
        """Over one e-fold of time, a species that changes by 2x gives R = ln 2."""
        y_prev = np.array([1.0, 1.0])
        y_now = np.array([1.5, 0.5])     # total moles deliberately unchanged
        r, worst = TerminationSteadyState.compute_residual(
            y_now, np.e, y_prev, 1.0, ATOL, labels=['a', 'b'])
        assert abs(r - np.log(2.0)) < 1e-12
        assert worst == 'b'              # the halving beats the 1.5x

    def test_residual_is_taken_on_mole_fractions_not_moles(self):
        """
        Doubling every species leaves the COMPOSITION untouched, so R must be 0. If the
        residual were taken on moles it would read ln 2 -- and a constant-pressure plasma
        whose total moles drift with the two-temperature equation of state would then look
        like it was still reacting when it was not.
        """
        y_prev = np.array([1.0, 1.0])
        r, _ = TerminationSteadyState.compute_residual(y_prev * 2.0, np.e, y_prev, 1.0, ATOL)
        assert r == 0.0

    def test_first_step_is_not_evaluable(self):
        """t_prev = 0 has no logarithm; the step carries no information either way."""
        y = np.array([1.0, 1.0])
        r, worst = TerminationSteadyState.compute_residual(y, 1e-15, y, 0.0, ATOL)
        assert np.isnan(r)
        assert worst is None

    def test_species_appearing_is_infinitely_far_from_steady(self):
        """A species crossing up through the floor is the opposite of stationary."""
        y_prev = np.array([1.0, 0.0])
        y_now = np.array([1.0, 1e-9])
        r, worst = TerminationSteadyState.compute_residual(
            y_now, np.e, y_prev, 1.0, ATOL, labels=['a', 'b'])
        assert np.isinf(r)
        assert worst == 'b'

    def test_species_below_the_floor_are_ignored(self):
        """
        Species beneath the integrator's absolute tolerance are noise, not chemistry. The
        argon deck's core carries three bath gases at exactly zero; without the floor their
        relative change is 0/0.
        """
        y_prev = np.array([1.0, 0.0, 1e-30])
        y_now = np.array([1.0, 0.0, 1e-25])       # a 1e5x change, entirely below the floor
        r, _ = TerminationSteadyState.compute_residual(y_now, np.e, y_prev, 1.0, ATOL)
        assert r == 0.0

    def test_a_large_negative_population_poisons_the_residual(self):
        """
        Finding 3. A negative population is a solver failure, not a settled composition. The
        live mask tests ``y > floor``, so a large-negative species would be dropped from the
        residual entirely -- neither counted nor reported -- and the criterion could declare
        steady state while it sat negative. Here two species are large and flat (so on their
        own they read as perfectly steady) while a third has gone negative beyond the floor;
        the residual must be inf, which resets the streak and can never arm, and must name the
        offending species rather than ignore it.
        """
        y_prev = np.array([1.0e6, 1.0e6, 0.3])
        y_now = np.array([1.0e6, 1.0e6, -0.3])    # species 'C' went negative beyond the floor
        r, worst = TerminationSteadyState.compute_residual(
            y_now, np.e, y_prev, 1.0, ATOL, labels=['A', 'B', 'C'])
        assert np.isinf(r), "a negative population must poison the residual, not be dropped from it"
        assert worst == 'C', "the residual must name the negative species, got {0!r}".format(worst)
        # And the latch consequence: fed to update(), an inf residual cannot arm or advance.
        term = TerminationSteadyState(tolerance=1e-6, window=2)
        term.armed = True                          # even a system that had armed earlier
        assert term.update(y_now, np.e, y_prev, 1.0, ATOL, labels=['A', 'B', 'C']) is False
        assert term.streak == 0

    def test_a_negative_below_the_floor_is_still_treated_as_noise(self):
        """
        The poison is for a population negative *beyond* the resolution floor. A tiny negative
        within ``[-floor, 0)`` is round-off around zero, not a failure, and must not trip the
        guard -- otherwise ordinary integrator noise would forbid steady state.
        """
        y_prev = np.array([1.0, 1.0])
        y_now = np.array([1.0, -ATOL / 2.0])       # negative, but within the floor: noise
        r, _ = TerminationSteadyState.compute_residual(y_now, np.e, y_prev, 1.0, ATOL)
        assert np.isfinite(r), "a sub-floor negative is noise and must not poison the residual"

    def test_no_live_species_is_not_evaluable(self):
        y = np.array([0.0, 0.0])
        r, _ = TerminationSteadyState.compute_residual(y, np.e, y, 1.0, ATOL)
        assert np.isnan(r)


class TerminationSteadyStateLatchTest:
    """The arming latch and the window -- the part that keeps the criterion honest."""

    def _feed(self, term, residual_sequence):
        """
        Drive `term` through a sequence of residuals by synthesising steps that produce
        them. One species at 1.0 and one at exp(+/-R) over an e-fold of time; the constant
        species dominates the total, so the mole fraction of the second tracks it.

        Physical time ADVANCES across steps -- ``t`` goes e^j -> e^(j+1), one e-fold per
        step -- because a criterion whose persistence is measured in physical time (so it
        does not depend on the integrator's step size) cannot be exercised by a clock that
        never moves. The clock belongs to THIS HARNESS, not to `term`: it is kept in
        ``self._feed_clock`` per term, and `term` is driven only through its real
        :meth:`update` API. An earlier version stamped the clock onto ``term._feed_step`` --
        an attribute production never sets -- which proved nothing about how the criterion is
        actually called; the harness owning the schedule is the honest arrangement.
        """
        clocks = self.__dict__.setdefault('_feed_clock', {})
        j = clocks.get(id(term), 0)
        fired_at = None
        for k, r in enumerate(residual_sequence):
            big = 1e12
            y_prev = np.array([big, 1.0])
            y_now = np.array([big, float(np.exp(r))])
            t_prev = float(np.exp(j))
            t_now = float(np.exp(j + 1))
            j += 1
            clocks[id(term)] = j
            if term.update(y_now, t_now, y_prev, t_prev, ATOL):
                fired_at = k
                break
        return fired_at

    def test_early_quiescence_does_not_fire_the_criterion(self):
        """
        THE regression test. On the measured lithium trajectory the residual is 2.83e-10
        at t = 3.0e-15 s -- step 2 of 497 -- because the integration has not yet reached
        the fastest chemistry, not because anything has converged. The true endpoint is
        1.84e-4 s. Without the latch the run stops there, 6e10 too early, and reports
        success. The tolerance cannot rescue it: excluding 2.83e-10 needs a tolerance
        below the deck's own rtol of 1e-8, i.e. below the integrator's noise.
        """
        term = TerminationSteadyState(tolerance=1e-6, window=3)
        # the quiescent head of the measured trajectory, verbatim
        assert self._feed(term, [2.116e-10, 2.830e-10, 9.988e-10] * 10) is None
        assert term.armed is False
        assert term.streak >= 3          # flat by the tolerance...
        assert term.residual < 1e-6      # ...and yet correctly refusing to fire

    def test_criterion_fires_once_armed_and_flat(self):
        """Transient (R rises past 1), then the flat tail: this is a real steady state.

        At the default window of two, and with each synthesised step spanning one e-fold of
        time, two flat samples fill the window AND span the e-fold the physical persistence
        requires, so firing lands on the second flat step (index 1)."""
        term = TerminationSteadyState(tolerance=1e-6, window=2)
        # rise through the measured peak of 15.5, then the measured collapse
        rise = [1e-9, 1e-4, 0.44, 2.13, 15.48, 3.45, 0.573, 3.328e-2, 1.355e-4]
        assert self._feed(term, rise) is None
        assert term.armed is True
        assert self._feed(term, [1.036e-7, 2.889e-8, 4.845e-9]) == 1   # 2 samples span an e-fold

    def test_two_flat_samples_are_the_minimum_and_not_sufficient(self):
        """Round 100 HIGH 2: the step count is no longer the criterion, only a fluke guard.
        A single flat sample is not an interval and cannot terminate however long its step;
        two samples are the minimum. But two samples are NOT sufficient on their own -- the
        physical span must also hold (tested in test_persistence_is_physical_time...)."""
        term = TerminationSteadyState(tolerance=1e-6, window=2)
        assert self._feed(term, [2.0]) is None            # arm it
        # one flat sample: even though this synthesised step spans a full e-fold, streak is
        # only 1, so the two-sample guard blocks it.
        assert self._feed(term, [1e-9]) is None
        assert term.streak == 1
        # the second flat sample makes an interval, and the e-fold span already holds
        assert self._feed(term, [1e-9]) == 0

    def test_a_single_moving_step_resets_the_window(self):
        term = TerminationSteadyState(tolerance=1e-6, window=3)
        self._feed(term, [2.0])
        self._feed(term, [1e-9, 1e-9])
        assert term.streak == 2
        self._feed(term, [0.5])                            # still moving
        assert term.streak == 0

    def test_never_converges_never_fires_and_keeps_its_residual(self):
        """
        A composition that keeps moving must never satisfy the criterion, and must leave
        behind the residual it reached so the run can say what it got to.
        """
        term = TerminationSteadyState(tolerance=1e-6, window=3)
        assert self._feed(term, [0.5] * 500) is None
        assert term.armed is False
        assert abs(term.residual - 0.5) < 1e-12

    def test_the_measured_argon_trajectory_never_fires(self):
        """
        The I-164 argon deck's composition is bit-identical at every step, so R is
        identically zero and the criterion never arms. That is correct: a model with no
        chemistry has not been shown to reach a stationary balance. The solver terminates
        such a run through the no-flux path instead, which reports NO steady state.
        """
        term = TerminationSteadyState(tolerance=1e-6, window=3)
        assert self._feed(term, [0.0] * 51) is None
        assert term.armed is False

    @staticmethod
    def _feed_r(term, target_r, t_prev, t_now, relax=float('nan')):
        """Produce a residual of exactly ``target_r`` at ANY step size: the residual is
        |d ln x / d ln t|, so the composition change must scale with the time step,
        d ln x = target_r * d ln t. (A fixed d ln x would make the residual blow up as the
        step shrinks -- the very step dependence under test.)"""
        big = 1e12
        dlnt = np.log(t_now) - np.log(t_prev)
        y_prev = np.array([big, 1.0])
        y_now = np.array([big, float(np.exp(target_r * dlnt))])
        return term.update(y_now, t_now, y_prev, t_prev, ATOL, relaxation_time=relax)

    def test_persistence_is_physical_time_not_a_bare_step_count(self):
        """Round 100 HIGH 2: the flat streak counted accepted solver STEPS, so a plateau
        spanning many steps in negligible physical time filled the window and fired while a
        plateau spanning few steps over a genuine relaxation did not -- the verdict depended
        on the step controller. The binding requirement is a physical span; the step count
        survives only as a two-sample fluke guard, never sufficient. And there is NO
        exact-zero waiver -- a residual of zero earns the same span confirmation as any other
        flat tail, because equal endpoints do not prove a frozen structure. (No relaxation
        time supplied here, so the span is one e-fold of absolute time.)"""
        term = TerminationSteadyState(tolerance=1e-6, window=3)
        self._feed_r(term, 2.0, 1.0, np.e)                  # arm
        # Flat steps in factor-1.1 increments. The two-sample guard is met at step 2 and the
        # bare step `window` of 3 at step 3, yet none fires until the flat run has spanned a
        # full e-fold of time -- which factor-1.1 steps reach only after ceil(1/ln 1.1) of
        # them. The verdict tracks the physical span, not the step count.
        t = np.e
        fired_step = None
        for k in range(30):
            t_next = t * 1.1
            got = self._feed_r(term, 1e-7, t, t_next)
            t = t_next
            if got:
                fired_step = k
                break
        assert term.streak >= term.window          # the step window filled long before...
        assert fired_step == int(np.ceil(1.0 / np.log(1.1))) - 1   # ...it fired at the e-fold

        # A residual of exactly zero gets NO waiver: it must span the e-fold like any tail.
        term0 = TerminationSteadyState(tolerance=1e-6, window=3)
        self._feed_r(term0, 2.0, 1.0, np.e)
        assert not self._feed_r(term0, 0.0, np.e, np.e * 1.1)      # frozen, but < 1 e-fold
        assert not self._feed_r(term0, 0.0, np.e * 1.1, np.e * 1.21)
        assert self._feed_r(term0, 0.0, np.e * 1.21, np.e * 3.0)   # now spans an e-fold

    def test_persistence_anchors_to_the_relaxation_time_when_supplied(self):
        """Round 100 HIGH 3: with a system relaxation time the flat run must persist for one
        tau of PHYSICAL time, not for a factor of e of absolute time. Anchored to tau, the
        confirmation is the same however large the absolute clock is -- it does not demand
        another 1.718*t0 just because the window opened late."""
        tau = 5.0
        # Window opens at a LARGE absolute time; only one tau of flatness is needed, not an
        # e-fold of the (large) absolute time.
        term = TerminationSteadyState(tolerance=1e-6, window=2)
        self._feed_r(term, 2.0, 1.0e6, 1.0e6 * np.e, relax=tau)          # arm at large t
        assert not self._feed_r(term, 1e-7, 1.0e6 * np.e, 1.0e6 * np.e + 2.0, relax=tau)  # 2 s < tau
        # cross one tau of elapsed flat time -> fires, though the absolute span is a tiny
        # fraction of an e-fold at t ~ 1e6.
        assert self._feed_r(term, 1e-7, 1.0e6 * np.e + 2.0, 1.0e6 * np.e + 6.0, relax=tau)

    def test_the_verdict_is_independent_of_when_the_flat_window_opens(self):
        """Round 100 HIGH 3 acceptance: the same physics -- arm, then a flat run that lasts a
        fixed number of relaxation times -- must reach the same verdict whether the flat
        window opens early or late in absolute time. Anchored to tau (not to an absolute
        e-fold, which would demand ever more elapsed time the later the window opens), it
        does."""
        tau = 3.0

        def fires_when_opening_at(t0):
            term = TerminationSteadyState(tolerance=1e-6, window=3)
            self._feed_r(term, 2.0, t0, t0 + tau, relax=tau)          # arm near t0
            t = t0 + tau
            for _ in range(8):
                t_next = t + 0.6 * tau
                if self._feed_r(term, 1e-7, t, t_next, relax=tau):
                    return True
                t = t_next
            return False

        assert fires_when_opening_at(1.0) is True
        assert fires_when_opening_at(1.0e3) is True     # far later in absolute time: same verdict
        assert fires_when_opening_at(1.0e9) is True

    def test_the_same_physics_gives_the_same_verdict_under_coarse_and_fine_stepping(self):
        """Round 100 HIGH 2 acceptance: a flat plateau one relaxation time long must reach
        the SAME verdict whether the integrator crossed it in few coarse steps or many fine
        ones. Driven here as two schedules over the identical physical trajectory."""
        tau = 4.0

        def run(step):
            term = TerminationSteadyState(tolerance=1e-6, window=3)
            self._feed_r(term, 2.0, 1.0, np.e, relax=tau)   # arm
            t = np.e
            for _ in range(int(round(3.0 * tau / step)) + 2):
                t_next = t + step
                if self._feed_r(term, 1e-8, t, t_next, relax=tau):
                    return True
                t = t_next
            return False

        assert run(0.5) is True         # fine: many small steps across the tau plateau
        assert run(3.9) is True         # coarse: a couple of big steps across the same plateau

    def test_the_external_arm_does_not_vouch_for_a_generic_channel_still_rising(self):
        """Round 100 HIGH 1: the external (electron) arm may license the criterion only while
        the generic channel is not still departing, judged over a WINDOW of samples, not the
        two most recent. A neutral residual rising throughout but carrying one flat/noisy
        sample must NOT arm the criterion; a genuinely settled channel does."""
        # rising sub-tolerance neutral with one flat sample, electron past relaxation
        term = TerminationSteadyState(tolerance=1e-6, window=3)
        rising = [1e-9, 2e-9, 3e-9, 3e-9, 4e-9, 5e-9, 6e-9, 7e-9, 8e-9, 9e-9]
        for k, rg in enumerate(rising):
            j = float(k)
            y_prev = np.array([1e12, 1.0])
            y_now = np.array([1e12, float(np.exp(rg))])
            term.update(y_now, np.exp(j + 1), y_prev, np.exp(j), ATOL, external_residual=1e-9, external_armed=True)
        assert term.armed_external is False and term.armed is False

        # a genuinely settled neutral, electron past relaxation: DOES arm
        term2 = TerminationSteadyState(tolerance=1e-6, window=3)
        settled = [0.5, 0.9, 0.95, 1e-9, 1e-9, 1e-9, 1e-9]
        for k, rg in enumerate(settled):
            j = float(k)
            y_prev = np.array([1e12, 1.0])
            y_now = np.array([1e12, float(np.exp(rg))])
            term2.update(y_now, np.exp(j + 1), y_prev, np.exp(j), ATOL, external_residual=1e-9, external_armed=True)
        assert term2.armed_external is True

    def test_the_external_arm_is_not_a_permanent_latch(self):
        """Round 100 HIGH 1: the external arm is a conjunction with a LIVE condition (the
        generic channel not departing), so it must be re-evaluated, not latched. A channel
        that settles (arming it) and then resumes moving must un-arm."""
        term = TerminationSteadyState(tolerance=1e-6, window=3)
        seq = [1e-9, 1e-9, 1e-9, 1e-9, 2e-9, 3e-9, 4e-9, 5e-9]   # settles, then departs
        for k, rg in enumerate(seq):
            j = float(k)
            y_prev = np.array([1e12, 1.0])
            y_now = np.array([1e12, float(np.exp(rg))])
            term.update(y_now, np.exp(j + 1), y_prev, np.exp(j), ATOL, external_residual=1e-9, external_armed=True)
        assert term.armed_external is False   # the resumed rise dropped the arm

    def test_a_species_rising_under_the_aggregate_maximum_still_blocks_the_arm(self):
        """Round 106 HIGH 1: the departing test must be PER SPECIES, not on the aggregate
        maximum residual. Here two neutrals sit below tolerance: one ('a') is RISING toward
        its own transient (log-log slope 5e-8 -> 1e-7 -> 1.5e-7 -> ...), the other ('b') is
        larger and FALLING (6.77e-7 -> 6.01e-7 -> 5.26e-7 -> ...). The aggregate maximum is
        'b' at every step, and it is falling -- so a departing test that reads only the
        maximum sees "nothing rising", stops treating the system as departing, and lets the
        external (electron) arm through. But 'a' is a live species climbing toward its own
        timescale: the system is NOT settled, and must not arm. The per-species trend catches
        the species the aggregate hides."""
        BIG = 1e12
        a_slopes = [5e-8, 1e-7, 1.5e-7, 2e-7, 2.5e-7, 3e-7, 3.5e-7, 4e-7]
        b_slopes = [6.77e-7, 6.01e-7, 5.26e-7, 5.0e-7, 4.8e-7, 4.7e-7, 4.6e-7, 4.55e-7]
        term = TerminationSteadyState(tolerance=1e-6, window=3)
        a = b = 1.0
        for k in range(len(a_slopes)):
            t_prev, t_now = float(np.exp(k)), float(np.exp(k + 1))     # dlnt = 1
            a_next, b_next = a * np.exp(a_slopes[k]), b * np.exp(b_slopes[k])
            term.update(np.array([BIG, a_next, b_next]), t_now,
                        np.array([BIG, a, b]), t_prev, ATOL,
                        labels=['big', 'a', 'b'], external_residual=1e-9, external_armed=True)
            a, b = a_next, b_next
        assert term.worst_label == 'b'          # the maximum is the falling species, as claimed
        assert term.armed_external is False      # ...yet the rising 'a' keeps the arm shut
        assert term.armed is False

        # Control: when BOTH species genuinely settle (every live slope flat for a full
        # window), the external arm IS granted -- the fix does not over-block a real steady
        # state.
        settled_a = [3e-7, 2e-7, 1e-7, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9]
        settled_b = [6e-7, 4e-7, 2e-7, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9]
        term2 = TerminationSteadyState(tolerance=1e-6, window=3)
        a = b = 1.0
        for k in range(len(settled_a)):
            t_prev, t_now = float(np.exp(k)), float(np.exp(k + 1))
            a_next, b_next = a * np.exp(settled_a[k]), b * np.exp(settled_b[k])
            term2.update(np.array([BIG, a_next, b_next]), t_now,
                         np.array([BIG, a, b]), t_prev, ATOL,
                         labels=['big', 'a', 'b'], external_residual=1e-9, external_armed=True)
            a, b = a_next, b_next
        assert term2.armed_external is True

    def test_the_window_parameter_is_honoured_as_the_minimum_sample_span(self):
        """Round 106 MEDIUM: ``window`` is documented as the minimum number of accepted
        samples a flat interval must span, but the production code hard-coded a floor of two
        regardless of the value passed, so window=4 terminated after two flat samples like
        window=2. Honoured, a larger window demands that many flat samples before firing --
        even when the physical span is already satisfied. Each synthesised step here spans a
        full e-fold, so the physical persistence (one e-fold) holds from the second flat
        sample on; only the sample-count floor separates window=2 from window=4."""
        def flat_sample_that_fires(window):
            term = TerminationSteadyState(tolerance=1e-6, window=window)
            assert self._feed(term, [2.0]) is None                      # arm
            return self._feed(term, [1e-9] * 10)                        # first flat step to fire
        assert flat_sample_that_fires(2) == 1     # two samples: fires on the 2nd flat step
        assert flat_sample_that_fires(4) == 3     # four samples: not until the 4th

    def test_a_species_appearing_while_the_external_channel_is_flat_does_not_return_steady(self):
        """Round 109 HIGH: a species crossing UP through the floor makes _slope_analysis return
        inf -- "emphatically not steady". The fold that adds the external (electron) channel
        tested ``not np.isfinite(r)``, which is true for inf as well as nan, so a small flat
        external residual REPLACED the inf poison instead of losing a MAX to it; and the
        per-species counters, updated only on an 'ok' step, kept their stale 'settled' values,
        so armed_external stayed True and the streak advanced. The appearance step then returned
        steady -- byte-identically to a control where nothing appears. A species crossing the
        floor is ordinary under model enlargement, so this is a live path. The appearance must
        NOT return steady; the control, genuinely settled, still must."""
        atol = 1e-16
        big, ext, tau = 1e12, 1e-9, 1e-30
        labels = ['big', 'settler', 'newcomer']

        def primed():
            # arm generic (settler climbs two decades), then one flat sample: armed and settled,
            # streak 1, one more flat step from terminating.
            t = TerminationSteadyState(tolerance=1e-6, window=2)
            t.update(np.array([big, 7.389, 0.0]), np.e ** 1, np.array([big, 1.0, 0.0]), np.e ** 0,
                     atol, labels=labels, external_residual=ext, external_armed=True, relaxation_time=tau)
            t.update(np.array([big, 7.389, 0.0]), np.e ** 2, np.array([big, 7.389, 0.0]), np.e ** 1,
                     atol, labels=labels, external_residual=ext, external_armed=True, relaxation_time=tau)
            return t

        # Control: nothing appears -> genuinely settled -> terminates.
        control = primed()
        assert control.update(np.array([big, 7.389, 0.0]), np.e ** 3, np.array([big, 7.389, 0.0]),
                              np.e ** 2, atol, labels=labels, external_residual=ext,
                              external_armed=True, relaxation_time=tau) is True

        # Appearance: the newcomer crosses up through the floor on this step.
        appear = primed()
        r_alone, lbl = appear.compute_residual(np.array([big, 7.389, atol * 10]), np.e ** 3,
                                               np.array([big, 7.389, 0.0]), np.e ** 2, atol, labels=labels)
        assert not np.isfinite(r_alone) and r_alone > 0.0 and lbl == 'newcomer'   # inf, the newcomer
        got = appear.update(np.array([big, 7.389, atol * 10]), np.e ** 3, np.array([big, 7.389, 0.0]),
                            np.e ** 2, atol, labels=labels, external_residual=ext,
                            external_armed=True, relaxation_time=tau)
        assert got is False                          # the inf survives the fold: not steady
        assert not np.isfinite(appear.residual)      # ...and is reported, not the electron's 1e-9
        assert appear.armed_external is False         # the appearance dropped the external arm

    def test_nan_and_inf_steps_treat_the_departing_counters_differently(self):
        """Round 109: a step the slope analysis cannot evaluate is either nan (NO information --
        degenerate interval, no live species) or inf (STRONG NEGATIVE information -- a species
        appeared or went negative). They are different answers and the counters do different
        things: a nan step leaves the per-species trend memory untouched (nothing was observed),
        an inf step DISCARDS it (the composition changed structurally, so a stale 'settled' count
        must not vouch for the new state)."""
        atol, big = 1e-16, 1e12
        labels = ['big', 's', 'new']

        def build():
            # two 'ok' steps build live counters for big and s (newcomer stays below the floor)
            t = TerminationSteadyState(tolerance=1e-6, window=2)
            t.update(np.array([big, 7.389, 0.0]), np.e ** 1, np.array([big, 1.0, 0.0]), np.e ** 0, atol, labels=labels)
            t.update(np.array([big, 7.389, 0.0]), np.e ** 2, np.array([big, 7.389, 0.0]), np.e ** 1, atol, labels=labels)
            return t

        # nan step: a degenerate interval (t_now == t_prev) -> no information -> counters UNCHANGED.
        t_nan = build()
        before = dict(t_nan._steps_since_rise)
        assert before != {}
        t_nan.update(np.array([big, 7.389, 0.0]), np.e ** 2, np.array([big, 7.389, 0.0]), np.e ** 2, atol, labels=labels)
        assert dict(t_nan._steps_since_rise) == before      # left exactly as they were

        # inf step: the newcomer appears -> structural change -> counters DISCARDED.
        t_inf = build()
        assert t_inf._steps_since_rise != {}
        t_inf.update(np.array([big, 7.389, atol * 10]), np.e ** 3, np.array([big, 7.389, 0.0]), np.e ** 2, atol, labels=labels)
        assert t_inf._steps_since_rise == {} and t_inf._slope_prev == {}

    def test_raising_window_makes_the_criterion_no_less_conservative_in_both_roles(self):
        """Round 109 MEDIUM 1: `window` gates two quantities -- the flat-streak length AND the
        number of consecutive non-rises before a species stops counting as 'departing'. Raising
        it makes BOTH more conservative (a longer flat streak is required, and a species is held
        as departing longer), so a larger window can only push termination LATER, never earlier;
        the single knob coherently means "how much evidence before trusting a 'no longer
        changing' verdict". (The finding's premise that the two move in opposite directions is
        inverted -- see the docstring on TerminationSteadyState.) This pins it: across a fixed
        trajectory driven through the external-arm gate, so both roles are active, the firing
        step is monotonically non-decreasing in `window`."""
        atol = 1e-16
        big, ext, tau = 1e12, 1e-9, 1e-30
        labels = ['b', 's']

        def fire_index(window):
            t = TerminationSteadyState(tolerance=1e-6, window=window)
            t.update(np.array([big, 7.389]), np.e ** 1, np.array([big, 1.0]), np.e ** 0, atol,
                     labels=labels, external_residual=ext, external_armed=True, relaxation_time=tau)
            for k in range(20):
                j = k + 1
                if t.update(np.array([big, 7.389]), np.e ** (j + 1), np.array([big, 7.389]), np.e ** j,
                            atol, labels=labels, external_residual=ext, external_armed=True,
                            relaxation_time=tau):
                    return k
            return None

        idx = [fire_index(w) for w in (2, 3, 4, 5)]
        assert all(a is not None for a in idx), idx
        assert all(idx[i] <= idx[i + 1] for i in range(len(idx) - 1)), idx   # never earlier
        assert idx[0] < idx[-1], idx                                          # and the knob is not inert

    def test_departing_counters_are_keyed_by_index_and_cleared_between_runs(self):
        """Round 109 MEDIUM 2: `_slope_prev`/`_steps_since_rise` are keyed by core-species
        integer index. Within one simulate() the core is fixed, so an index means the same
        species for the whole integration. The core only GROWS across simulate() boundaries
        (model enlargement), and each simulate() restarts from t=0 having called reset() first
        (base.pyx). This grows the core between steps to show the hazard concretely -- a stale
        counter at an index is silently re-read against whatever species now holds that index --
        and the guard that removes it: reset() clears both dicts, and a grown core only ever
        arrives on a fresh simulate()."""
        atol, big = 1e-16, 1e12

        # Run 1: a two-species core builds a counter at index 1 (species A).
        t = TerminationSteadyState(tolerance=1e-6, window=2)
        t.update(np.array([big, 7.389]), np.e ** 2, np.array([big, 7.389]), np.e ** 1, atol, labels=['b', 'A'])
        assert 1 in t._steps_since_rise

        # WITHOUT reset, feeding an ENLARGED core re-reads index 1 -- now species B -- and
        # inherits A's stale count. This is exactly the hazard the finding names.
        t.update(np.array([big, 3.0, 7.389]), np.e ** 3, np.array([big, 3.0, 7.389]), np.e ** 2, atol, labels=['b', 'B', 'A'])
        assert t._steps_since_rise.get(1, 0) >= 2      # B inherited the count that belonged to A

        # The production contract removes it: base.pyx calls reset() before every simulate(),
        # and a grown core only ever arrives on a fresh simulate(). After reset the dicts are
        # empty, so no index can carry a meaning from a previous, smaller core.
        t.reset()
        assert t._steps_since_rise == {} and t._slope_prev == {}
        t.update(np.array([big, 3.0, 7.389]), np.e ** 2, np.array([big, 3.0, 7.389]), np.e ** 1, atol, labels=['b', 'B', 'A'])
        assert set(t._steps_since_rise) <= {0, 1, 2}   # only the new core's indices, nothing stale

    def test_reset_clears_the_latch(self):
        """
        Each simulate() restarts from t=0. A latch carried over from the previous
        enlargement's simulation would arm the criterion for a trajectory that never
        earned it.
        """
        term = TerminationSteadyState(tolerance=1e-6, window=3)
        self._feed(term, [2.0])
        assert term.armed is True
        term.reset()
        assert term.armed is False
        assert term.streak == 0
        assert np.isnan(term.residual)

    def test_rejects_nonsense_configuration(self):
        with pytest.raises(ValueError):
            TerminationSteadyState(tolerance=0.0)
        with pytest.raises(ValueError):
            TerminationSteadyState(tolerance=-1e-6)
        with pytest.raises(ValueError):
            TerminationSteadyState(tolerance=float('nan'))
        with pytest.raises(ValueError):
            TerminationSteadyState(window=0)


# ---------------------------------------------------------------------------------
# End-to-end, through the solver
# ---------------------------------------------------------------------------------

def _thermo(h298_kj, s298):
    return ThermoData(
        Tdata=([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0], 'K'),
        Cpdata=([20.8] * 7, 'J/(mol*K)'),
        H298=(h298_kj, 'kJ/mol'),
        S298=(s298, 'J/(mol*K)'),
    )


def _relaxing_system(termination):
    """
    A plasma reactor that genuinely relaxes: an ionisation/recombination pair on the
    charged side and a reversible thermal pair on the neutral side, so both the electron
    chemistry and ordinary chemistry have to settle before the criterion may fire.
    """
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    ar_ion = Species(label='Ar+').from_adjacency_list('1 Ar u1 p3 c+1')
    spc_a = Species(label='A').from_adjacency_list('1 Ar u0 p4 c0')
    spc_a.thermo = _thermo(0.0, 150.0)
    spc_b = Species(label='B').from_adjacency_list('1 Ar u0 p4 c0')
    spc_b.thermo = _thermo(-15.0, 145.0)

    core_species = [ar, electron, ar_ion, spc_a, spc_b]
    core_reactions = [
        Reaction(reactants=[electron, ar], products=[ar_ion, electron, electron],
                 reversible=False,
                 kinetics=TwoTemperaturePlasma(A=(1.0e-3, 'm^3/(mol*s)'), n=0.5,
                                               Ea_g=(4000.0, 'J/mol'), Ea_e=(60000.0, 'J/mol'))),
        Reaction(reactants=[electron, ar_ion], products=[ar], reversible=False,
                 kinetics=TwoTemperaturePlasma(A=(5.0e5, 'm^3/(mol*s)'), n=-0.5,
                                               Ea_g=(0.0, 'J/mol'), Ea_e=(0.0, 'J/mol'))),
        Reaction(reactants=[spc_a], products=[spc_b], reversible=True,
                 kinetics=Arrhenius(A=(100.0, 's^-1'), n=0.0, Ea=(10.0, 'kJ/mol'))),
    ]
    imf = {electron: 1.0e-4, ar: 1.0, ar_ion: 1.0e-4, spc_a: 0.1, spc_b: 0.05}
    reactor = PlasmaReactor(T_GAS, P0, imf, (T_E, 'K'), n_sims=1, termination=termination)
    return reactor, core_species, core_reactions


def _inert_system(termination):
    """The shape of the argon deck: a core with species but not one reaction."""
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    ar_ion = Species(label='Ar+').from_adjacency_list('1 Ar u1 p3 c+1')
    core_species = [ar, electron, ar_ion]
    imf = {electron: 1.0e-4, ar: 1.0, ar_ion: 1.0e-4}
    reactor = PlasmaReactor(T_GAS, P0, imf, (T_E, 'K'), n_sims=1, termination=termination)
    return reactor, core_species, []


def _simulate(reactor, core_species, core_reactions):
    return reactor.simulate(
        core_species, core_reactions, [], [], [], [],
        model_settings=ModelSettings(tol_keep_in_edge=0, tol_move_to_core=1e5,
                                     tol_interrupt_simulation=1e8),
        simulator_settings=SimulatorSettings(),
    )


def _hydrocarbon(smiles, cp, h298, s298):
    tdata = ([300, 400, 500, 600, 800, 1000, 1500], "K")
    return Species(molecule=[Molecule().from_smiles(smiles)],
                   thermo=ThermoData(Tdata=tdata, Cpdata=(cp, "cal/(mol*K)"),
                                     H298=(h298, "kcal/mol"), S298=(s298, "cal/(mol*K)")))


def _inert_core_with_nonempty_edge(termination):
    """
    The shape of a real inert plasma deck: empty chemistry (no flux anywhere) but a NON-EMPTY
    edge. The core has no reactions, and the single edge reaction carries no flux (A = 0), so
    ``char_rate`` is zero and every edge rate is zero, while the edge species list is not
    empty. This is the case the empty-edge :func:`_inert_system` cannot exercise: with a
    non-empty edge the zero-flux promotion block runs first and breaks (reporting ZERO FLUX),
    so control never reached the steady-state termination logic. Returns the full (core, edge)
    argument set for ``simulate``.
    """
    ch4 = _hydrocarbon("C", [8.615, 9.687, 10.963, 12.301, 14.841, 16.976, 20.528], -17.714, 44.472)
    c2h6 = _hydrocarbon("CC", [12.684, 15.506, 18.326, 20.971, 25.500, 29.016, 34.595], -19.521, 54.799)
    ch3 = _hydrocarbon("[CH3]", [9.397, 10.123, 10.856, 11.571, 12.899, 14.055, 16.195], 9.357, 45.174)
    # A structurally-present edge reaction that carries no flux (A = 0): the edge is non-empty
    # but the chemistry is empty, so char_rate and every edge rate are exactly zero.
    edge_rxn = Reaction(reactants=[c2h6], products=[ch3, ch3],
                        kinetics=Arrhenius(A=(0.0, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K")))
    core_species = [ch4, c2h6]
    edge_species = [ch3]
    reactor = SimpleReactor(T=1000.0, P=1.0e5, initial_mole_fractions={ch4: 0.5, c2h6: 0.5},
                            n_sims=1, termination=termination)
    reactor.initialize_model(core_species, [], edge_species, [edge_rxn])
    return reactor, core_species, [], edge_species, [edge_rxn]


class SteadyStateTerminationInSolverTest:
    """The criterion driven by the real integrator, not by synthetic steps."""

    def test_relaxing_system_terminates_on_steady_state(self):
        """It stops before the backstop, and it stops with a residual to report."""
        steady = TerminationSteadyState(tolerance=1e-6, window=3)
        backstop = TerminationTime((1.0e4, 's'))
        reactor, core_species, core_reactions = _relaxing_system([steady, backstop])
        terminated = _simulate(reactor, core_species, core_reactions)[0]

        assert terminated
        assert reactor.steady_state_reached is True
        assert reactor.steady_state_residual < 1e-6
        assert reactor.t < 1.0e4                    # the backstop was NOT what stopped it
        assert steady.armed is True                 # it went through a real transient

    def test_never_converges_terminates_on_the_backstop_without_claiming_success(self):
        """
        Verifier 2. An unreachable tolerance must not hang, and must not report success:
        the run ends on its backstop time, `steady_state_reached` stays False, and the
        residual it actually reached survives for the log to name.
        """
        steady = TerminationSteadyState(tolerance=1e-30, window=3)
        backstop = TerminationTime((1.0e-6, 's'))
        reactor, core_species, core_reactions = _relaxing_system([steady, backstop])
        terminated = _simulate(reactor, core_species, core_reactions)[0]

        assert terminated                                   # it stopped
        assert reactor.steady_state_reached is False        # and did not claim victory
        assert reactor.t > 1.0e-6                           # on the backstop
        assert np.isfinite(reactor.steady_state_residual)   # naming the residual it got to

    def test_no_flux_system_terminates_but_demonstrates_no_steady_state(self):
        """
        A core with no reactions cannot change, so the integration stops -- but a system
        that never started has converged to nothing. `steady_state_reached` must stay
        False: reporting "trivially stationary" as a satisfied criterion is exactly the
        confusion this criterion exists to remove.
        """
        steady = TerminationSteadyState(tolerance=1e-6, window=3)
        backstop = TerminationTime((1.0, 's'))
        reactor, core_species, core_reactions = _inert_system([steady, backstop])
        terminated = _simulate(reactor, core_species, core_reactions)[0]

        assert terminated
        assert reactor.steady_state_reached is False

    def test_absent_criterion_leaves_the_reactor_untouched(self):
        """A deck that does not ask for steady state behaves exactly as it did before."""
        reactor, core_species, core_reactions = _relaxing_system([TerminationTime((1.0e-6, 's'))])
        _simulate(reactor, core_species, core_reactions)
        assert reactor.steady_state_reached is False
        assert np.isnan(reactor.steady_state_residual)

    def test_latch_does_not_leak_between_simulations(self):
        """
        The termination object outlives one simulation and is reused across the whole
        model-generation run, so simulate() must clear the latch on entry.
        """
        steady = TerminationSteadyState(tolerance=1e-6, window=3)
        reactor, core_species, core_reactions = _relaxing_system(
            [steady, TerminationTime((1.0e4, 's'))])
        _simulate(reactor, core_species, core_reactions)
        assert steady.armed is True

        reactor2, core_species2, core_reactions2 = _inert_system(
            [steady, TerminationTime((1.0, 's'))])
        _simulate(reactor2, core_species2, core_reactions2)
        assert steady.armed is False
        assert reactor2.steady_state_reached is False

    def test_inert_core_with_nonempty_edge_reports_no_steady_state_instead_of_promoting(self, caplog):
        """
        Finding 2. The steady-state "NO STEADY STATE WAS DEMONSTRATED" report was unreachable
        whenever the edge was non-empty: the zero-flux promotion block runs BEFORE the
        steady-state termination logic and breaks, so an inert plasma deck (empty chemistry
        but a non-empty edge) took the promotion path instead of reporting honestly. The
        existing no-flux coverage uses an empty edge, which never exercises that ordering.

        With a steady-state criterion declared, the zero-flux block must stand down and let the
        steady-state block own the no-flux exit: nothing is promoted, `steady_state_reached`
        stays False, and the honest report is emitted.
        """
        steady = TerminationSteadyState(tolerance=1e-6, window=3)
        backstop = TerminationTime((1.0, 's'))
        reactor, core_species, core_reactions, edge_species, edge_reactions = \
            _inert_core_with_nonempty_edge([steady, backstop])

        # Guard the premise: the edge is non-empty (so the zero-flux block's len()>0 gate is
        # met and it runs first), exactly the ordering the empty-edge case never reaches.
        assert len(edge_species) > 0

        with caplog.at_level(logging.INFO):
            terminated, _res, invalid_objects, _ss, _sr, _t, _x = reactor.simulate(
                core_species, core_reactions, edge_species, edge_reactions, [], [],
                model_settings=ModelSettings(tol_keep_in_edge=0, tol_move_to_core=1e5,
                                             tol_interrupt_simulation=1e8),
                simulator_settings=SimulatorSettings(),
            )

        assert np.all(reactor.core_species_rates == 0.0), "the core is meant to be inert here"
        assert np.all(reactor.edge_species_rates == 0.0), "empty chemistry: no edge flux either"
        assert terminated
        assert invalid_objects == [], (
            "an inert core with a non-empty edge promoted {0!r} instead of reporting that no "
            "steady state was demonstrated".format(invalid_objects)
        )
        assert reactor.steady_state_reached is False
        assert not any(
            "added to model core to avoid singularity" in r.getMessage() for r in caplog.records
        ), "the zero-flux promotion path fired even though a steady-state criterion was declared"
        assert not any(
            r.getMessage().startswith("ZERO FLUX:") for r in caplog.records
        ), "the zero-flux block reported instead of deferring to the steady-state criterion"
        assert any(
            "NO STEADY STATE WAS DEMONSTRATED" in r.getMessage() for r in caplog.records
        ), "the honest steady-state report was never reached: {0!r}".format(
            [r.getMessage() for r in caplog.records]
        )


# ===========================================================================
# Round 110 -- the value is ignored, the permission is not
# ===========================================================================

def test_an_unusable_external_residual_cannot_authorise_termination():
    """Round 110 HIGH 1. Round 109 made the residual FOLD ignore a non-finite
    ``external_residual``, but left the boolean ``external_armed`` -- which arrives on the
    same ``update()`` call and is independent -- to grant the arm on its own. So a channel
    that has just reported NO USABLE NUMBER (``nan``) or a value that blew up (``inf``) still
    authorises termination on a flat generic channel, byte-identically to a valid residual.

    Reproduces the byte-identical verdict list the reviewer measured, then pins the fix: the
    value and the permission must travel together, so an external arm is licensed only while
    the external residual is a usable finite number.
    """
    y = np.array([1.0, 0.5])                 # two live species, flat (identical) across steps
    t = [1.0e-6, 1.0e-3, 1.0, 1.0e3]         # log-spaced, so the e-fold persistence span holds
    floor = 1.0e-30

    def verdicts(external_residual):
        term = TerminationSteadyState(tolerance=1e-6, window=2)
        return [bool(term.update(y, t[k], y, t[k - 1], floor,
                                 external_residual=external_residual, external_armed=True,
                                 relaxation_time=float('nan')))
                for k in range(1, len(t))]

    finite = verdicts(1.0e-9)
    nan = verdicts(float('nan'))
    pos_inf = verdicts(float('inf'))
    neg_inf = verdicts(float('-inf'))

    # A usable finite external residual still arms and terminates -- the valid path is intact.
    assert any(finite), finite
    # ...but a nan (no usable number), +inf or -inf (invalid) must NOT authorise termination:
    # a non-finite residual can never grant convergence through any channel.
    assert nan == [False, False, False], nan
    assert pos_inf == [False, False, False], pos_inf
    assert neg_inf == [False, False, False], neg_inf
    # The finite verdict list is no longer byte-identical to the non-finite ones, which was
    # the defect (all four were [False, True, True] before the fix).
    assert finite != nan and finite != pos_inf and finite != neg_inf


def test_a_nonfinite_external_residual_poisons_even_an_armed_generic_channel():
    """Round 110 HIGH 1 + round 111 HIGH 1. A non-finite external residual must never authorise
    termination through ANY channel -- and in particular must not be dropped as 'criterion
    unavailable' while the remaining (generic) channel goes on to terminate the run, which is
    exactly the old behaviour. Arm the generic channel through a real transient and take it
    flat, then have the external channel report a non-finite value: the folded residual is
    poisoned, the step cannot terminate, and the invalid value is named in the diagnostic.

    This now covers the WHOLE non-finite class, not just the enumerated infinities (round 111):
      * ``+inf``/``-inf`` poison UNCONDITIONALLY -- only an active channel can compute an
        infinite slope, so an ordinary reactor never emits one and gating on the arm is
        unnecessary;
      * ``nan`` poisons when the channel is ARMED -- a live discharge whose electron has passed
        ``t*nu_wall >= 1`` yet produced no usable number (the subnormal-fraction underflow that
        round 110's own addendum rerouted from ``+inf`` onto the ``nan`` branch). Round 110
        poisoned the infinities and left this ``nan`` on the benign 'no information' path, so a
        flat generic channel still fired on it. A guard written as a comparison lets ``nan``
        through by construction; the fix dispatches on ``np.isfinite`` instead."""
    big = 1.0e12

    def feed(term, target_r, t_prev, t_now, external_residual, external_armed):
        dlnt = np.log(t_now) - np.log(t_prev)
        y_prev = np.array([big, 1.0])
        y_now = np.array([big, float(np.exp(target_r * dlnt))])
        return term.update(y_now, t_now, y_prev, t_prev, ATOL,
                           external_residual=external_residual, external_armed=external_armed,
                           relaxation_time=float('nan'))

    for poison, armed in ((float('inf'), False), (float('-inf'), False), (float('nan'), True)):
        term = TerminationSteadyState(tolerance=1e-6, window=2)
        # Arm the generic channel through a real transient, with a benign (unarmed) external.
        feed(term, 2.0, 1.0, np.e, float('nan'), False)
        assert term.armed_generic is True
        # Generic now flat, but the external channel reports a non-finite residual it cannot
        # stand behind: no termination, even though the armed generic channel alone would fire.
        v1 = feed(term, 0.0, np.e, np.e ** 2, poison, armed)
        v2 = feed(term, 0.0, np.e ** 2, np.e ** 3, poison, armed)
        assert v1 is False and v2 is False, (poison, armed, v1, v2)
        assert not np.isfinite(term.residual), (poison, armed)   # the fold is poisoned, not dropped
        assert 'external channel' in str(term.worst_label), (poison, term.worst_label)


def test_a_nan_external_residual_from_an_armed_channel_cannot_authorise_termination():
    """Round 111 HIGH 1, the reviewer's exact reproduction. The non-finite authorisation
    matrix that round 110 left behind read ``+inf -> 0, -inf -> 0, nan -> 1 unsafe path``: a
    ``nan`` external residual still authorised termination on one path. Arm the generic channel
    with ``R = 2`` over ``[1, e]``, then supply two flat generic intervals with
    ``external_residual = nan`` while the external channel is ARMED (``external_armed=True`` --
    the electron has passed its relaxation time but its fraction underflowed to an unresolvable
    zero). The second flat interval returned ``True`` before the fix. A channel that is vouching
    yet cannot produce a finite number must certify nothing."""
    big = 1.0e12

    def feed(term, target_r, t_prev, t_now):
        dlnt = np.log(t_now) - np.log(t_prev)
        y_prev = np.array([big, 1.0])
        y_now = np.array([big, float(np.exp(target_r * dlnt))])
        return term.update(y_now, t_now, y_prev, t_prev, ATOL,
                           external_residual=float('nan'), external_armed=True,
                           relaxation_time=float('nan'))

    term = TerminationSteadyState(tolerance=1e-6, window=2)
    feed(term, 2.0, 1.0, np.e)                     # R = 2 over [1, e] arms the generic channel
    assert term.armed_generic is True
    v1 = feed(term, 0.0, np.e, np.e ** 2)          # generic now flat; external channel = nan
    v2 = feed(term, 0.0, np.e ** 2, np.e ** 3)     # the interval that returned True before
    assert v1 is False and v2 is False, (v1, v2)
    # Both must be a Python bool, not a numpy.bool_ leaked from the relaxation-time fallback.
    assert type(v2) is bool, type(v2)
    assert not np.isfinite(term.residual)          # the fold is poisoned, not silently dropped
    assert 'external channel' in str(term.worst_label), term.worst_label


def test_update_returns_a_python_bool_on_the_relaxation_time_fallback():
    """Round 111 MEDIUM. The numpy.bool_ leak was fixed at its source for the ``armed_external``
    attribute (round 110), but the DEFAULT relaxation-time fallback still computes ``span_ok``
    as ``(np.log(t_now) - np.log(t_flat_start)) >= 1.0`` -- a numpy.float64 comparison yielding a
    numpy.bool_ -- and returned ``self.armed and span_ok`` unchanged. So a genuine termination on
    that path returned a numpy.bool_, not a Python bool, and an ``is True`` caller would fail.
    ``update()`` must return a Python bool on EVERY path: the not-evaluable early return, the
    not-yet-persisted return, and the fired verdict on both the supplied-relaxation-time branch
    and the e-fold fallback."""
    # An ordinary reactor: no external channel, no relaxation time -> the e-fold fallback path.
    term = TerminationSteadyState(tolerance=1e-6, window=2)
    flat = np.array([1.0, 0.5])
    dlnt = np.log(np.e) - np.log(1.0)
    moved = np.array([1.0, float(np.exp(2.0 * dlnt))])
    assert term.update(moved, np.e, flat, 1.0, 1e-30) is False        # a real transient arms generic
    assert term.armed_generic is True
    v1 = term.update(flat, np.e ** 2, flat, np.e, 1e-30)              # flat, streak 1: not yet
    v2 = term.update(flat, np.e ** 3, flat, np.e ** 2, 1e-30)        # flat, streak 2: fires on fallback
    assert v1 is False, v1
    assert v2 is True, v2
    assert type(v1) is bool and type(v2) is bool, (type(v1), type(v2))


def test_rate_ratio_criterion_evaluates_only_for_a_finite_positive_denominator():
    """Round 110 HIGH 2 (owner's ruling; the denominator matrix). The edge/core/network rate
    ratios are dimensionless only when divided by a finite, strictly positive characteristic
    rate. ``_rate_ratios_or_zero`` evaluates the ratio ONLY then; for a zero, negative, or
    non-finite denominator the ratio is undefined, so the relative criterion ABSTAINS and
    returns zeros -- promoting and terminating nothing -- rather than dividing by a floored
    1.0/epsilon that would compare a dimensional rate against a dimensionless tolerance."""
    from rmgpy.solver.base import ReactionSystem
    rates = np.array([2.0, -4.0, 0.0])

    # positive denominator -> the ordinary applicable ratio, unchanged.
    np.testing.assert_allclose(ReactionSystem._rate_ratios_or_zero(rates, 4.0), [0.5, 1.0, 0.0])
    # zero / negative / NaN / +inf / -inf denominator -> undefined -> abstain -> all zeros,
    # and finite (no NaN or inf laundered into argmax or the branching numbers).
    for denom in (0.0, -4.0, float('nan'), float('inf'), float('-inf')):
        out = ReactionSystem._rate_ratios_or_zero(rates, denom)
        assert np.all(out == 0.0), (denom, out)
        assert np.isfinite(out).all(), (denom, out)


def _zero_core_flux_with_large_edge_rate(termination):
    """Round 110 HIGH 2 fixture: a NON-PLASMA reactor whose core carries no reaction
    (``char_rate`` exactly 0, hence ``total_char_rate`` 0 -- no transport), but whose
    non-empty edge has a LARGE flux. Edge reactions do not feed core-species derivatives
    (simple.pyx), so the core rate stays exactly zero while the edge rate is real -- and here
    large enough (~1.2e7, A=1e6) that, divided by the old floored denominator of 1.0, it would
    clear a dimensionless ``tol_move_to_core``: the dimensional comparison the fix removes."""
    ch4 = _hydrocarbon("C", [8.615, 9.687, 10.963, 12.301, 14.841, 16.976, 20.528], -17.714, 44.472)
    c2h6 = _hydrocarbon("CC", [12.684, 15.506, 18.326, 20.971, 25.500, 29.016, 34.595], -19.521, 54.799)
    ch3 = _hydrocarbon("[CH3]", [9.397, 10.123, 10.856, 11.571, 12.899, 14.055, 16.195], 9.357, 45.174)
    edge_rxn = Reaction(reactants=[c2h6], products=[ch3, ch3],
                        kinetics=Arrhenius(A=(1.0e6, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K")))
    core_species = [ch4, c2h6]
    edge_species = [ch3]
    reactor = SimpleReactor(T=1000.0, P=1.0e5, initial_mole_fractions={ch4: 0.5, c2h6: 0.5},
                            n_sims=1, termination=termination)
    reactor.initialize_model(core_species, [], edge_species, [edge_rxn])
    return reactor, core_species, [], edge_species, [edge_rxn]


def test_a_zero_core_flux_reactor_does_not_promote_by_dimensional_comparison(caplog):
    """Round 110 HIGH 2 (owner's ruling; matrix case 'zero-core-flux reactor'). With
    ``char_rate == 0`` the ratio criterion is undefined and must ABSTAIN: it may neither
    promote an edge species nor terminate. The old ``char_rate ... else 1.0`` compared the raw
    dimensional edge rate against a dimensionless ``tol_move_to_core`` and would have promoted
    this large edge rate; after the fix nothing is promoted through the ratio, and the
    dimensioned ABSOLUTE criterion (the zero-flux / steady-state-inert band on
    ``total_char_rate``) governs the outcome deterministically instead."""
    steady = TerminationSteadyState(tolerance=1e-6, window=2)
    backstop = TerminationTime((1.0, 's'))
    reactor, core_species, core_reactions, edge_species, edge_reactions = \
        _zero_core_flux_with_large_edge_rate([steady, backstop])

    # Guard the premise: char_rate is exactly zero, and the edge rate is large enough that the
    # old dimensional comparison edge/1.0 WOULD have exceeded tol_move_to_core below.
    y0 = np.array(reactor.y0, float)
    reactor.residual(0.0, y0.copy(), np.zeros_like(y0))
    assert np.all(reactor.core_species_rates == 0.0), "core must be inert (char_rate == 0)"
    max_edge = float(np.max(np.abs(reactor.edge_species_rates)))
    assert max_edge > 1.0e5, "edge rate must exceed tol_move_to_core so the dimensional bug WOULD fire"

    with caplog.at_level(logging.INFO):
        _terminated, _res, invalid_objects, _ss, _sr, _t, _x = reactor.simulate(
            core_species, core_reactions, edge_species, edge_reactions, [], [],
            model_settings=ModelSettings(tol_keep_in_edge=0, tol_move_to_core=1e5,
                                         tol_interrupt_simulation=1e8),
            simulator_settings=SimulatorSettings(),
        )

    # The undefined ratio must not promote the edge species by its raw dimensional magnitude.
    assert edge_species[0] not in invalid_objects, (
        "a zero-core-flux reactor promoted {0!r} through a dimensional ratio comparison; the "
        "ratio criterion must abstain when its denominator is zero".format(invalid_objects))
    # The dimensioned absolute criterion governs deterministically instead.
    assert any('this model never started' in r.getMessage() for r in caplog.records), \
        "the absolute zero-flux criterion did not govern the zero-core-flux case"


def _reversible_surface_reaction(h298_b, termination, x_a=0.5, x_b=0.5):
    """A reversible reaction ``A <=> B`` whose two species carry IDENTICAL Cp/S and a settable
    B enthalpy, with B declared a SURFACE species on that reaction. At ``h298_b == 0`` the two
    thermo sets are identical (Keq == 1); started at equal moles the forward and reverse fluxes
    are then equal at every step, so the NET rate of each species is exactly zero -- hence
    ``char_rate == 0`` -- while the GROSS production and consumption of B are positive and equal.
    A small non-zero ``h298_b`` breaks the tie: ``char_rate`` becomes positive and the surface
    ratio ``max(production, consumption) / char_rate`` is a finite, large, applicable number."""
    tdata = ([300, 400, 500, 600, 800, 1000, 1500], "K")

    def spc(smiles, h):
        return Species(molecule=[Molecule().from_smiles(smiles)],
                       thermo=ThermoData(Tdata=tdata, Cpdata=([12.0] * 7, "cal/(mol*K)"),
                                         H298=(h, "kcal/mol"), S298=(50.0, "cal/(mol*K)")))

    a, b = spc("CC", 0.0), spc("CCC", h298_b)
    rxn = Reaction(reactants=[a], products=[b], reversible=True,
                   kinetics=Arrhenius(A=(1.0e3, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K")))
    reactor = SimpleReactor(T=1000.0, P=1.0e5, initial_mole_fractions={a: x_a, b: x_b},
                            n_sims=1, termination=termination)
    reactor.initialize_model([a, b], [rxn], [], [], [b], [rxn])
    return reactor, [a, b], [rxn], b


def _surface_model_settings():
    # use_dynamics on (a finite edge->core tolerance), and a surface->core species tolerance the
    # ratio is measured against; the reaction->core tolerance is set out of reach so only the
    # SPECIES promotion is exercised.
    return ModelSettings(tol_keep_in_edge=0, tol_move_to_core=1e5, tol_interrupt_simulation=1e8,
                         tol_move_edge_rxn_to_core=0.5, tol_move_surface_spc_to_core=0.2,
                         tol_move_surface_rxn_to_core=1e9)


def test_a_reversible_surface_reaction_at_equal_flux_does_not_divide_by_zero(caplog):
    """Round 111 HIGH 2. The surface-species-to-core promotion ratio in base.pyx divided
    ``max(|production|, |consumption|)`` by ``char_rate`` DIRECTLY, bypassing the round-110
    abstaining helper -- the second entry point the helper's install did not route through. On a
    reversible surface reaction with equal forward and reverse flux the net rates cancel
    (``char_rate == 0``) while the gross rates are positive, so the division is ``positive / 0``:
    before the fix the production entry raised ZeroDivisionError (a numpy divide would instead
    have promoted B through an undefined criterion). Routed through ``_rate_ratios_or_zero`` it
    abstains -- the ratio is zero, nothing is promoted -- and the run completes."""
    reactor, core, rxns, b = _reversible_surface_reaction(0.0, [TerminationTime((1.0e-3, 's'))])
    # Premise: char_rate is exactly zero while the surface species' gross rates are positive.
    y0 = np.array(reactor.y0, float)
    reactor.residual(0.0, y0.copy(), np.zeros_like(y0))
    assert float(np.sqrt(np.sum(np.asarray(reactor.core_species_rates, float) ** 2))) == 0.0, \
        "equal forward/reverse flux must give a zero characteristic chemistry rate"
    ib = reactor.get_species_index(b)
    assert np.asarray(reactor.core_species_production_rates, float)[ib] > 0.0
    assert np.asarray(reactor.core_species_consumption_rates, float)[ib] > 0.0

    with caplog.at_level(logging.INFO):
        result = reactor.simulate(core, rxns, [], [], [b], rxns,
                                  model_settings=_surface_model_settings(),
                                  simulator_settings=SimulatorSettings())
    # The undefined ratio abstains: the surface species is NOT moved to the core...
    assert not any('Moving species' in r.getMessage() and str(b) in r.getMessage()
                   for r in caplog.records), "a 0/0 surface ratio promoted the surface species"
    # ...and the run completes to its backstop instead of raising ZeroDivisionError.
    assert result[0] is True


def test_a_surface_species_is_promoted_when_the_characteristic_rate_is_positive(caplog):
    """Round 111 HIGH 2, the positive-denominator production-path arm. Abstaining on a zero
    denominator is only correct if the criterion still PROMOTES where its denominator is defined
    -- otherwise a guard broken to never promote is indistinguishable from a working one. Break
    the forward/reverse tie by a small B enthalpy so ``char_rate > 0`` and the surface ratio is a
    finite, large, applicable number: B is moved from surface to core, on the same production
    path, proving the fix removed only the undefined case."""
    reactor, core, rxns, b = _reversible_surface_reaction(-0.2, [TerminationTime((1.0e-3, 's'))])
    y0 = np.array(reactor.y0, float)
    reactor.residual(0.0, y0.copy(), np.zeros_like(y0))
    assert float(np.sqrt(np.sum(np.asarray(reactor.core_species_rates, float) ** 2))) > 0.0, \
        "a broken tie must give a positive characteristic chemistry rate (an applicable denominator)"

    with caplog.at_level(logging.INFO):
        reactor.simulate(core, rxns, [], [], [b], rxns,
                         model_settings=_surface_model_settings(),
                         simulator_settings=SimulatorSettings())
    assert any('Moving species' in r.getMessage() and str(b) in r.getMessage()
               for r in caplog.records), \
        "a surface species with a defined, large rate ratio was not promoted -- the guard is inert"


def test_steady_state_report_names_the_criterion_that_fired(caplog):
    """Round 110 MEDIUM. Any of several steady-state criteria can terminate, but the readback
    and the success log always used criterion zero (``steady_state_terms[0]``). With a first
    criterion whose tolerance (1e-30) it never meets and a second (1e-6) it does, the success
    line quoted criterion zero -- 'below tolerance 1.0000e-30 for 0 consecutive steps' -- a
    statement false on its face. Report the criterion that actually fired."""
    never = TerminationSteadyState(tolerance=1e-30, window=3)   # criterion 0: never satisfied
    fires = TerminationSteadyState(tolerance=1e-6, window=3)     # criterion 1: the one that fires
    backstop = TerminationTime((1.0e4, 's'))
    reactor, core_species, core_reactions = _relaxing_system([never, fires, backstop])

    with caplog.at_level(logging.INFO):
        terminated = _simulate(reactor, core_species, core_reactions)[0]

    assert terminated
    assert reactor.steady_state_reached is True
    msgs = [r.getMessage() for r in caplog.records if 'reached steady state' in r.getMessage()]
    assert msgs, "no steady-state success line was logged: {0!r}".format(
        [r.getMessage() for r in caplog.records])
    reached = msgs[-1]
    # The report must quote the tolerance and streak of the criterion that FIRED (1e-6, a full
    # window), not criterion zero's (1e-30, zero consecutive steps).
    assert '1.0000e-06' in reached, reached
    assert '1.0000e-30' not in reached, reached
    assert 'for 0 consecutive steps' not in reached, reached
    # steady_state_residual is the fired criterion's, and it cleared 1e-6.
    assert reactor.steady_state_residual < 1e-6


def test_an_armed_false_nan_external_residual_cannot_authorise_termination():
    """Round 112 HIGH 1, first reproduction (the steady-state decision site). Round 111 poisoned
    a ``nan`` external residual only when the channel was ARMED; with ``external_armed=False`` a
    supplied ``nan`` slipped onto the benign 'no information' path and a flat, armed generic
    channel fired on it. The value could not be told apart from 'no external channel', so the
    guard keyed on the arm FLAG. The fix makes ABSENCE a distinct sentinel (``external_residual
    is None``): a channel that is SUPPLIED (any non-None value) and non-finite authorises
    nothing, armed or not; an ABSENT channel (None, the default) leaves the generic channel to
    decide, byte-for-byte. Before the fix the second flat interval returned ``True``."""
    big = 1.0e12

    def feed(term, target_r, t_prev, t_now, **ext):
        dlnt = np.log(t_now) - np.log(t_prev)
        y_prev = np.array([big, 1.0])
        y_now = np.array([big, float(np.exp(target_r * dlnt))])
        return term.update(y_now, t_now, y_prev, t_prev, ATOL, relaxation_time=float('nan'), **ext)

    # A SUPPLIED nan with the channel UNARMED must refuse -- the round-111 gap.
    term = TerminationSteadyState(tolerance=1e-6, window=2)
    feed(term, 2.0, 1.0, np.e, external_residual=float('nan'), external_armed=False)
    assert term.armed_generic is True
    v1 = feed(term, 0.0, np.e, np.e ** 2, external_residual=float('nan'), external_armed=False)
    v2 = feed(term, 0.0, np.e ** 2, np.e ** 3, external_residual=float('nan'), external_armed=False)
    assert v1 is False and v2 is False, (v1, v2)
    assert type(v2) is bool, type(v2)
    assert not np.isfinite(term.residual), term.residual
    assert 'external channel' in str(term.worst_label), term.worst_label

    # Control: an ABSENT channel (external_residual omitted -> None) is the ordinary reactor and
    # MUST still terminate on the flat generic channel -- the byte-for-byte pre-plasma path.
    ctrl = TerminationSteadyState(tolerance=1e-6, window=2)
    feed(ctrl, 2.0, 1.0, np.e)
    assert ctrl.armed_generic is True
    feed(ctrl, 0.0, np.e, np.e ** 2)
    assert feed(ctrl, 0.0, np.e ** 2, np.e ** 3) is True


def test_a_nonfinite_per_species_slope_poisons_the_step_instead_of_dropping_it():
    """Round 112 HIGH 1, second reproduction (the same steady-state decision, via the slope
    analysis). ``_slope_analysis`` took each live species' log-log slope and DROPPED any that
    came out non-finite, letting a finite neighbour govern -- so a run could be called steady
    while one species blew up. With moles ``[1e308, 1e-20]`` and a floor below ``1e-20`` the
    second species is LIVE (its moles clear the floor) but its mole FRACTION ``1e-20 / 1e308``
    underflows to zero, ``log(0) = -inf``, and ``inf - inf`` makes its slope ``nan``. The finite
    heavy species then held the maximum and, once armed and flat, two intervals returned
    ``False, True``. A non-finite slope must poison the whole step, never be dropped."""
    floor = 1.0e-30
    big = 1.0e308

    def feed(term, y_prev, t_prev, y_now, t_now):
        return term.update(np.asarray(y_now, dtype=float), t_now,
                           np.asarray(y_prev, dtype=float), t_prev, floor,
                           labels=['heavy', 'trace'])

    term = TerminationSteadyState(tolerance=1e-6, window=2)
    # Arm the generic channel with a real transient while both fractions are O(1) and both live.
    feed(term, [1.0, 1.0], 1.0, [float(np.exp(2.0)), 1.0], np.e)
    assert term.armed_generic is True
    # Now the heavy species is flat while the trace fraction underflows: slope[trace] = nan.
    v1 = feed(term, [big, 1e-20], np.e, [big, 1e-20], np.e ** 2)
    v2 = feed(term, [big, 1e-20], np.e ** 2, [big, 1e-20], np.e ** 3)
    assert v1 is False and v2 is False, (v1, v2)
    assert not np.isfinite(term.residual), term.residual   # poisoned, not dropped to a finite max
    assert 'trace' in str(term.worst_label), term.worst_label


def test_rate_ratios_refuse_a_nonfinite_ratio_over_a_finite_denominator():
    """Round 112 HIGH 2 (the interrupt, network-promotion and surface-promotion decision sites).
    Every enlargement/interrupt/promotion ratio -- core, edge, network-leak and surface -- routes
    through ``_rate_ratios_or_zero``; it is the single decision-input gate for all four sites. A
    finite, positive denominator does NOT guarantee a finite ratio: a +/-inf rate makes the ratio
    non-finite (the network-leak ``lr > tol_interrupt`` interrupt), and a finite but huge gross
    rate OVERFLOWS to +inf even over a finite denominator (the surface ``ratio > tol`` promotion).
    A non-finite ratio then interrupts the run and promotes numerical garbage. The helper must
    refuse -- stop loudly -- rather than return it, while still ABSTAINING (zeros) on a
    zero/negative/non-finite denominator and passing finite ratios straight through."""
    from rmgpy.solver.base import ReactionSystem
    ratios = ReactionSystem._rate_ratios_or_zero

    # Interrupt / network-leak site: a +/-inf or nan rate over a finite, positive denominator.
    for bad in (np.inf, -np.inf, np.nan):
        with pytest.raises(ValueError, match='(?i)non-finite'):
            ratios(np.array([1.0, bad, 2.0]), 4.0)

    # Surface-promotion site: finite gross rate, finite denominator, the RATIO overflows to +inf.
    with pytest.raises(ValueError, match='(?i)non-finite'):
        ratios(np.array([1.0e308, 1.0e308]), 1.0e-300)

    # Round 113 BLOCKING 2: a non-finite NUMERATOR stops loudly REGARDLESS of the denominator.
    # The zero/non-finite-denominator abstention is NOT a licence to launder an inf rate to zero
    # (core_species_rates=[0], char_rate=0, network_leak_rates=[inf] must not quietly promote or
    # interrupt on a zeroed inf). This inverts the round-112 assertion, which blessed exactly that
    # laundering.
    for denom in (0.0, -1.0, float('nan'), float('inf')):
        with pytest.raises(ValueError, match='(?i)non-finite'):
            ratios(np.array([np.inf, 1.0]), denom)
        # ...but a FINITE numerator over a zero/non-finite denominator still ABSTAINS (zeros).
        np.testing.assert_array_equal(ratios(np.array([2.0, 1.0]), denom), np.zeros(2))
    # Preserved: finite rates over a finite positive denominator pass straight through.
    np.testing.assert_allclose(ratios(np.array([2.0, -4.0, 0.0]), 4.0), [0.5, 1.0, 0.0])
