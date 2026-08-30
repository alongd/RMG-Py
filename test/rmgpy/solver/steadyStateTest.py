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

import numpy as np
import pytest

import rmgpy.constants as constants
from rmgpy.kinetics import Arrhenius
from rmgpy.kinetics.arrhenius import TwoTemperaturePlasma
from rmgpy.reaction import Reaction
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
from rmgpy.solver.plasma import PlasmaReactor
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

    def test_no_live_species_is_not_evaluable(self):
        y = np.array([0.0, 0.0])
        r, _ = TerminationSteadyState.compute_residual(y, np.e, y, 1.0, ATOL)
        assert np.isnan(r)


class TerminationSteadyStateLatchTest:
    """The arming latch and the window -- the part that keeps the criterion honest."""

    @staticmethod
    def _feed(term, residual_sequence):
        """
        Drive `term` through a sequence of residuals by synthesising steps that produce
        them. One species at 1.0 and one at exp(+/-R) over an e-fold of time; the constant
        species dominates the total, so the mole fraction of the second tracks it.
        """
        fired_at = None
        for k, r in enumerate(residual_sequence):
            big = 1e12
            y_prev = np.array([big, 1.0])
            y_now = np.array([big, float(np.exp(r))])
            if term.update(y_now, np.e, y_prev, 1.0, ATOL):
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
        """Transient (R rises past 1), then the flat tail: this is a real steady state."""
        term = TerminationSteadyState(tolerance=1e-6, window=3)
        # rise through the measured peak of 15.5, then the measured collapse
        rise = [1e-9, 1e-4, 0.44, 2.13, 15.48, 3.45, 0.573, 3.328e-2, 1.355e-4]
        assert self._feed(term, rise) is None
        assert term.armed is True
        assert self._feed(term, [1.036e-7, 2.889e-8, 4.845e-9]) == 2   # window of 3

    def test_window_is_enforced(self):
        """Two flat steps are not three."""
        term = TerminationSteadyState(tolerance=1e-6, window=3)
        assert self._feed(term, [2.0]) is None            # arm it
        assert self._feed(term, [1e-9, 1e-9]) is None     # only two flat steps
        assert term.streak == 2
        assert self._feed(term, [1e-9]) == 0              # the third fires

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
