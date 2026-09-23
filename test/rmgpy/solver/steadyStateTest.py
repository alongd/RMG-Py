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
            term.update(y_now, np.exp(j + 1), y_prev, np.exp(j), ATOL, external_armed=True)
        assert term.armed_external is False and term.armed is False

        # a genuinely settled neutral, electron past relaxation: DOES arm
        term2 = TerminationSteadyState(tolerance=1e-6, window=3)
        settled = [0.5, 0.9, 0.95, 1e-9, 1e-9, 1e-9, 1e-9]
        for k, rg in enumerate(settled):
            j = float(k)
            y_prev = np.array([1e12, 1.0])
            y_now = np.array([1e12, float(np.exp(rg))])
            term2.update(y_now, np.exp(j + 1), y_prev, np.exp(j), ATOL, external_armed=True)
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
            term.update(y_now, np.exp(j + 1), y_prev, np.exp(j), ATOL, external_armed=True)
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
                        labels=['big', 'a', 'b'], external_armed=True)
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
                         labels=['big', 'a', 'b'], external_armed=True)
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
