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

import json
import multiprocessing
import os
import signal
import time
from pathlib import Path

import numpy as np
import pytest
from rmgpy.tools.eedf.validation import split_branches
from rmgpy.tools.eedf.branches import (
    CERTIFICATE_LIMITS,
    AdapterContractError,
    PlasmaReactorAdapter,
    ReactorRun,
    Seed,
    continue_branch,
    declared_branch,
    matching_branch,
    path_consistency,
    record_continuation,
    scan,
    stability_analysis,
)


def _complete(record):
    """Declare what every adapter must supply: half-step Jacobian and conservation."""
    record = dict(record)
    record.setdefault('_half_step_jacobian', record['_jacobian'])
    record.setdefault('_conservation_constraints', [])
    return record


class _BlockingAdapter:
    """A subprocess fixture whose worker cannot cooperate with shutdown."""

    def __init__(self, pid_directory):
        self.pid_directory = Path(pid_directory)

    def integrate(self, seed, power=None):
        (self.pid_directory / 'worker-{}.pid'.format(int(seed.u))).write_text(
            str(os.getpid()))
        if seed.u == 0.:
            signal.signal(signal.SIGTERM, signal.SIG_IGN)
        while True:
            time.sleep(1.)


class _PostResultAdapter:
    def __init__(self, pid_path):
        self.pid_path = Path(pid_path)


def _send_result_then_block(connection, adapter, seed, power):
    adapter.pid_path.write_text(str(os.getpid()))
    connection.send(('ok', _complete({
        'u': seed.u,
        'n_e': seed.n_e,
        '_state_vector': [seed.u],
        '_jacobian': [[-1.]],
    })))
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    while True:
        time.sleep(1.)


def _run_blocking_scan(output, pid_directory, queue):
    started = time.monotonic()
    payload = scan(
        _BlockingAdapter(pid_directory), (0., 1.), output, u_steps=2,
        electron_density_decades=1, processes=2, timeout=.05, reference_power=.5,
        recheck_timeouts=False)
    queue.put((time.monotonic() - started, payload['attempts']))


def _run_post_result_scan(output, pid_path, queue):
    import rmgpy.tools.eedf.branches as branches

    branches._integration_child = _send_result_then_block
    started = time.monotonic()
    payload = branches.scan(
        _PostResultAdapter(pid_path), (0., 0.), output, u_steps=1,
        electron_density_decades=1, processes=1, timeout=.05, reference_power=.5)
    queue.put((time.monotonic() - started, payload['attempts']))


def _watchdog_scan(target, args, pid_paths):
    """Run a potentially hanging baseline in an outer process and reap it."""
    context = multiprocessing.get_context('fork')
    queue = context.Queue()
    runner = context.Process(target=target, args=(*args, queue))
    runner.start()
    runner.join(1.)
    hung = runner.is_alive()
    if hung:
        runner.kill()
        runner.join(.2)
    if hung:
        for path in pid_paths:
            path = Path(path)
            if not path.exists():
                continue
            pid = int(path.read_text())
            try:
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
    assert not hung, 'scan exceeded its deadline because worker shutdown blocked'
    assert runner.exitcode == 0
    return queue.get(timeout=.2)


def test_two_different_scan_paths_are_not_averaged():
    def row(energy):
        return {'swarm': {'mean_energy_eV': energy}, 'k_ine': np.array([energy * 1e-15]),
                'k_sup': np.array([0.]), 'power_groups': {'field': energy * 1e-16},
                'f0': np.array([energy]), 'channel_power': np.array([energy * 1e-16]),
                'attachment_energy_eV': np.array([0.])}
    low = [row(1.), row(2.)]
    high = [row(3.), row(4.)]
    branches, verdict = split_branches({'up': low, 'down': high[::-1], 'cold': low},
                                      {'rtol': 1e-6, 'atol': 0.})
    assert len(branches) == 2
    assert not verdict['agreement']
    assert branches['branch_0'][1]['swarm']['mean_energy_eV'] == 2.
    assert branches['branch_1'][1]['swarm']['mean_energy_eV'] == 4.


def test_reactor_scan_finds_stable_and_unstable_states_and_fold(tmp_path):
    """The analytic saddle-node x**2=P has known roots, stability and fold."""
    class Synthetic:
        def __init__(self):
            self.calls = []

        def integrate(self, seed, power=None):
            self.calls.append((seed.u, power))
            power = 1.0 if power is None else power
            if power < 0.:
                raise RuntimeError('no real equilibrium below the fold')
            x = (-1. if seed.u < 0. else 1.) * np.sqrt(power)
            return {'u': x, 'n_e': 1.e16 * np.exp(x), '_state_vector': [x],
                    'self_sustained': True, 'discharge_class': 'self-sustained'}

        def jacobian(self, state):
            return np.array([[2.0 * state['u']]])

        half_step_jacobian = jacobian

        @staticmethod
        def conservation_constraints(state):
            return []

        @staticmethod
        def is_fold_failure(failure):
            return 'no real equilibrium' in failure

    synthetic = Synthetic()
    result = scan(synthetic, (-1., 1.), tmp_path / 'branches.json', u_steps=2,
                  electron_density_decades=2, powers=(1., 0.), reference_power=1.,
                  continuation_stability_resolution=.05 / 64.)
    assert len(result['branches']) == 2
    assert [item['stability']['stable'] for item in result['branches']] == [True, False]
    stable = result['branches'][0]
    assert synthetic.calls[4] == (-1., 1.)
    assert stable['continuation_path'][0]['converged']
    assert stable['continuation_path'][-1]['fold_detected'] is True
    lower, upper = stable['continuation_path'][-1]['fold_bracket']
    assert lower <= 0. <= upper
    assert upper - lower < 1.e-3
    assert declared_branch(tmp_path / 'branches.json', stable['id'])['id'] == stable['id']


def test_failed_seed_is_reported_not_dropped(tmp_path):
    class Broken:
        def integrate(self, seed, power=None):
            raise RuntimeError('synthetic nonconvergence')

        def jacobian(self, state):
            raise AssertionError('unreachable')

    result = scan(Broken(), (0., 0.), tmp_path / 'branches.json', u_steps=1,
                  electron_density_decades=1, reference_power=.5)
    assert result['branches'] == []
    assert result['attempts'][0]['converged'] is False
    assert 'synthetic nonconvergence' in result['attempts'][0]['failure']


def test_conservation_modes_are_removed_before_stability_is_decided():
    jacobian = np.array([[-1., 1., 0.], [1., -1., 0.], [0., 0., -3.]])
    constraints = [{'name': 'element:Ar', 'vector': [1., 1., 0.]}]

    result = stability_analysis(jacobian, constraints, half_step_jacobian=jacobian)

    assert result['classification'] == 'stable'
    assert result['stable'] is True
    assert result['conservation_modes'] == ['element:Ar']
    assert len(result['raw_eigenvalues']) == 3
    assert len(result['eigenvalues']) == 2
    assert result['step_halving_same_verdict'] is True


def test_nonconservation_zero_mode_is_marginal_not_stable():
    result = stability_analysis(np.diag([-2., 0.]), [])

    assert result['classification'] == 'marginal'
    assert result['stable'] is False


def test_step_halving_disagreement_is_marginal_not_admissible():
    result = stability_analysis(
        np.array([[-1.]]), [], half_step_jacobian=np.array([[1.]]))

    assert result['full_step_classification'] == 'stable'
    assert result['half_step_classification'] == 'unstable'
    assert result['classification'] == 'marginal'
    assert result['stable'] is False
    assert result['step_halving_same_verdict'] is False


def test_continuation_timeout_is_recorded():
    class Slow:
        def integrate(self, seed, power=None):
            time.sleep(.2)

    result = continue_branch(
        Slow(), {'u': 0., 'n_e': 1.e16, 'state': [0.]}, [.5], timeout=.02)

    assert result[0]['converged'] is False
    assert 'TimeoutError' in result[0]['failure']


def test_plasma_adapter_runs_engine_and_reads_engine_discharge_class():
    class Reactor:
        def __init__(self):
            self.t = 2.
            self.y = np.array([2., 3., 0.5])
            self.dydt = np.zeros(3)
            self.atol_array = np.full(3, 1.e-16)
            self.num_core_species = 2
            self.electron_index = 0
            self.te_index = 2
            self.steady_state_reached = False
            self.energy_budget = {'n_e': 4.e15}
            self.simulated = False
            self.initial_mole_fractions = {}
            self.electron_kinetics = {'initial_reduced_field': (1., 'Td')}

        def simulate(self, *args, **kwargs):
            assert self.initial_mole_fractions == {'e-': 2., 'Ar': 3.}
            assert self.electron_kinetics['initial_reduced_field'][0] == pytest.approx(
                np.exp(.5))
            assert self.electron_kinetics['initial_reduced_field'][1] == 'Td'
            self.simulated = True
            self.steady_state_reached = True

        def residual(self, t, y, dydt):
            return np.array([-2. * (y[0] - 2.), -(y[1] - 3.), -3. * (y[2] - .5)]), 0

        def _jacobian(self, t, y, dydt, cj):
            return np.diag([-2., -1., -3.])

        def discharge_state(self):
            return 'source-supported'

    reactor = Reactor()

    def factory(seed, power):
        assert seed == Seed(0.5, 4.e15, (2., 3., .5))
        assert power == .5
        return ReactorRun(reactor, ['e-', 'Ar'], [], object(), object())

    state = PlasmaReactorAdapter(factory).integrate(
        Seed(.5, 4.e15, (2., 3., .5)), power=.5)

    assert reactor.simulated
    assert state['u'] == .5
    assert state['n_e'] == 4.e15
    assert state['residual_max_abs'] == 0.
    assert state['discharge_class'] == 'source-supported'
    assert state['self_sustained'] is False
    assert np.asarray(state['_jacobian']) == pytest.approx(np.diag([-2., -1., -3.]))
    assert np.asarray(state['_half_step_jacobian']) == pytest.approx(np.diag([-2., -1., -3.]))
    assert state['_state_vector'] == pytest.approx([2., 3., .5])


def test_conservation_verification_uses_resolved_perturbation_ensemble():
    class Reactor:
        t = 0.
        y = np.array([1., 1.])
        atol_array = np.array([1.e-16, 1.e-16])
        num_core_species = 1

        def residual(self, t, state, dydt):
            delta = state[0] - self.y[0]
            if delta == 0.:
                # An almost-cancelled steady residual has a poor pointwise relative
                # ratio even though resolved perturbations conserve exactly.
                return np.array([1.e-20, -0.9999998e-20]), 0
            return np.array([delta, -delta]), 0

    constraint = {'name': 'element:Ar', 'vector': [1., 1.]}

    verified = PlasmaReactorAdapter(lambda *_: None)._verify_constraints(
        Reactor(), [constraint])

    assert verified[0]['verified'] is True
    assert verified[0]['verification_relative_residual'] < 1.e-7


def test_scan_refuses_more_than_six_worker_processes(tmp_path):
    with pytest.raises(ValueError, match='at most 6'):
        scan(object(), (0., 1.), tmp_path / 'branches.json', processes=7)


@pytest.mark.parametrize('offset', [0., 100.])
def test_continuation_follows_cubic_branch_and_stops_at_bracketed_fold(offset):
    """u'=P+u-u**3 must stop at, not jump across, its lower fold."""
    class Cubic:
        def integrate(self, seed, power=None):
            roots = np.roots([1., 0., -1., -float(power)])
            real = [offset + root.real for root in roots if abs(root.imag) < 1.e-10]
            predicted = seed.state[0]
            u = min(real, key=lambda value: abs(value - predicted))
            physical_u = u - offset
            return _complete({
                'u': u,
                'n_e': 1.e16,
                'self_sustained': True,
                'discharge_class': 'self-sustained',
                '_state_vector': [u],
                '_jacobian': [[1. - 3. * physical_u * physical_u]],
            })

    path = continue_branch(
        Cubic(), {'u': offset - 1., 'n_e': 1.e16, 'state': [offset - 1.]}, [0.5],
        initial_power=0., max_step=0.05, min_step=1.e-4,
        predictor_tolerance=0.2)

    assert path[-1]['fold_detected'] is True
    lower, upper = path[-1]['fold_bracket']
    exact_fold = 2. / np.sqrt(27.)
    assert lower < exact_fold < upper
    assert upper - lower <= 0.01
    assert all(entry['u'] < offset for entry in path if entry.get('converged'))
    assert not any(entry.get('power', 0.) > upper for entry in path)


def test_continuation_checks_first_interval_for_an_eigenvalue_crossing():
    class FirstStepCrossing:
        def integrate(self, seed, power=None):
            power = float(power)
            return _complete({
                'u': power,
                'n_e': 1.e16,
                '_state_vector': [power],
                '_jacobian': [[power - 0.05]],
            })

    path = continue_branch(
        FirstStepCrossing(), {'u': 0., 'n_e': 1.e16, 'state': [0.]}, [.1],
        initial_power=0., max_step=.1)

    assert path[-1]['event_type'] == 'stability'
    assert path[-1]['fold_detected'] is False
    assert path[-1]['failure'] == 'reduced eigenvalue crossed zero'
    lower, upper = path[-1]['event_bracket']
    assert lower <= .05 <= upper
    assert upper - lower <= .1 / 64.


@pytest.mark.parametrize('max_step', [.025, .0125])
@pytest.mark.parametrize('state_scale', [.02, .002])
def test_continuation_scale_aware_tangent_brackets_rescaled_cubic_fold(state_scale, max_step):
    """Changing state units cannot turn a branch jump into an accepted step."""
    offset = 2.

    class RescaledCubic:
        def integrate(self, seed, power=None):
            roots = np.roots([1., 0., -1., -float(power)])
            states = [offset + state_scale * root.real
                      for root in roots if abs(root.imag) < 1.e-10]
            u = min(states, key=lambda value: abs(value - seed.state[0]))
            physical = (u - offset) / state_scale
            return _complete({
                'u': u,
                'n_e': 1.e16,
                '_state_vector': [u],
                '_jacobian': [[1. - 3. * physical * physical]],
            })

    path = continue_branch(
        RescaledCubic(),
        {'u': offset - state_scale, 'n_e': 1.e16,
         'state': [offset - state_scale]},
        [.5], initial_power=0., max_step=max_step, min_step=1.e-4)

    event = path[-1]
    exact_fold = 2. / np.sqrt(27.)
    assert event['event_type'] == 'fold'
    assert event['fold_detected'] is True
    assert event['fold_bracket'][0] < exact_fold < event['fold_bracket'][1]
    assert event['fold_bracket'][1] - event['fold_bracket'][0] <= .001
    assert not any(entry.get('power', 0.) > event['fold_bracket'][1]
                   for entry in path)


def test_continuation_certification_brackets_paired_stability_crossing():
    """Endpoint-stable samples cannot hide the unstable interval (0.01, 0.015)."""
    class PairedCrossing:
        def integrate(self, seed, power=None):
            power = float(power)
            eigenvalue = -10000. * (power - .01) * (power - .015)
            return _complete({
                'u': 0.,
                'n_e': 1.e16,
                '_state_vector': [0.],
                '_jacobian': [[eigenvalue]],
            })

    path = continue_branch(
        PairedCrossing(), {'u': 0., 'n_e': 1.e16, 'state': [0.]},
        [.025], initial_power=0., max_step=.025, min_step=.00025)

    event = path[-1]
    assert event['event_type'] == 'stability'
    assert event['event_bracket'][0] < .01 < event['event_bracket'][1]
    assert event['event_bracket'][1] - event['event_bracket'][0] <= .00025
    assert event['fold_detected'] is False


def test_continuation_certification_finds_off_center_stability_excursion():
    """A narrow unstable interval away from the midpoint must be sampled."""
    class OffCenterCrossing:
        def integrate(self, seed, power=None):
            power = float(power)
            eigenvalue = -1. + 2. * np.exp(-((power - .006) / .0005) ** 2)
            return _complete({
                'u': power,
                'n_e': 1.e16,
                '_state_vector': [power],
                '_jacobian': [[eigenvalue]],
            })

    path = continue_branch(
        OffCenterCrossing(), {'u': 0., 'n_e': 1.e16, 'state': [0.]},
        [.025], initial_power=0., max_step=.025, min_step=.00025)

    event = path[-1]
    assert event['event_type'] == 'stability'
    crossing = .006 - np.sqrt(np.log(2.)) * .0005
    assert event['event_bracket'][0] <= crossing <= event['event_bracket'][1]
    assert event['event_bracket'][1] - event['event_bracket'][0] <= .00025


def test_default_certification_grid_matches_the_declared_minimum_step():
    class StableCounter:
        def __init__(self):
            self.calls = 0

        def integrate(self, seed, power=None):
            self.calls += 1
            return _complete({
                'u': 0.,
                'n_e': 1.e16,
                '_state_vector': [0.],
                '_jacobian': [[-1.]],
            })

    adapter = StableCounter()
    path = continue_branch(
        adapter, {'u': 0., 'n_e': 1.e16, 'state': [0.]},
        [.1], initial_power=0., max_step=.1)

    assert path[-1]['converged'] is True
    assert path[-1]['stability_resolution'] == pytest.approx(.1 / 64.)
    # Two initial solves (state and its noise re-solve), the predictor solve,
    # then 64 chain samples, each with one reverse solve.
    assert adapter.calls == 2 + 1 + 64 + 64
    assert path[-1]['certificate']['sub_steps'] == 64
    assert path[-1]['certificate_solves'] == 65


def test_fast_single_valued_branch_crossing_is_not_labeled_fold():
    """A fast state tangent alone is not evidence of a saddle-node fold."""
    class FastStablePath:
        def integrate(self, seed, power=None):
            power = float(power)
            u = np.exp(20. * power)
            return _complete({
                'u': u,
                'n_e': 1.e16,
                '_state_vector': [u],
                '_jacobian': [[power - .05]],
            })

    path = continue_branch(
        FastStablePath(), {'u': 1., 'n_e': 1.e16, 'state': [1.]},
        [.1], initial_power=0., max_step=.025, min_step=.001)

    event = path[-1]
    assert event['event_type'] == 'stability'
    assert event['fold_detected'] is False
    assert event['event_bracket'][0] <= .05 <= event['event_bracket'][1]


def test_steep_smooth_stable_branch_completes_without_a_false_fold():
    class SmoothStable:
        def integrate(self, seed, power=None):
            u = 2. + .5 * np.tanh(2000. * float(power))
            return _complete({
                'u': u,
                'n_e': 1.e16,
                '_state_vector': [u],
                '_jacobian': [[-1.]],
            })

    path = continue_branch(
        SmoothStable(), {'u': 2., 'n_e': 1.e16, 'state': [2.]},
        [.5], initial_power=0.)

    assert path[-1]['power'] == pytest.approx(.5)
    assert path[-1]['converged'] is True
    assert all(entry.get('event_type') is None for entry in path)


@pytest.mark.parametrize('steps', [{}, {'max_step': .025},
                                   {'max_step': .01, 'min_step': .001}])
def test_smooth_state_extremum_is_neither_a_fold_nor_unresolved(steps):
    """A stable branch whose state peaks keeps one orientation in (P, x)."""
    class Peak:
        def integrate(self, seed, power=None):
            u = 1. - 100. * (float(power) - .05) ** 2
            return _complete({
                'u': u,
                'n_e': 1.e16,
                '_state_vector': [u],
                '_jacobian': [[-1.]],
            })

    path = continue_branch(
        Peak(), {'u': .75, 'n_e': 1.e16, 'state': [.75]},
        [.1], initial_power=0., **steps)

    assert path[-1]['power'] == pytest.approx(.1)
    assert path[-1]['converged'] is True
    assert all(entry.get('event_type') is None for entry in path)


def test_continuation_snaps_roundoff_sized_remainder_to_exact_target():
    """Accumulated power roundoff cannot create a spurious terminal step."""
    class LinearStable:
        def integrate(self, seed, power=None):
            power = float(power)
            return _complete({
                'u': 1. + power,
                'n_e': 1.e16,
                '_state_vector': [1. + power],
                '_jacobian': [[-1.]],
            })

    path = continue_branch(
        LinearStable(), {'u': 1.5, 'n_e': 1.e16, 'state': [1.5]},
        [.25, .5, 1.], initial_power=.5, max_step=.0125,
        min_step=.001, stability_resolution=.0125)

    assert path[-1]['power'] == 1.
    assert path[-1]['converged'] is True
    assert len(path) == 81
    assert all(entry.get('event_type') is None for entry in path)


def test_continuation_refines_a_tangent_departure_without_an_eigenvalue_signal():
    """A branch hop is refused even when both branches linearise alike."""
    class SilentHop:
        def integrate(self, seed, power=None):
            power = float(power)
            u = 10. + power + (0. if power < .3 else .2)
            return _complete({
                'u': u,
                'n_e': 1.e16,
                '_state_vector': [u],
                '_jacobian': [[-1.]],
            })

    path = continue_branch(
        SilentHop(), {'u': 10., 'n_e': 1.e16, 'state': [10.]},
        [.5], initial_power=0.)

    event = path[-1]
    assert event['event_type'] == 'unresolved'
    assert event['fold_detected'] is False
    assert event['event_bracket'][0] < .3 <= event['event_bracket'][1]
    assert event['event_bracket'][1] - event['event_bracket'][0] <= 2. * .5 / 20. / 64.
    assert all(entry['u'] - 10. - entry['power'] < .1
               for entry in path if entry.get('converged'))


def test_continuation_refines_a_hop_seen_only_in_a_trace_component():
    """Departure is judged per component, so a 1e-9-sized species counts."""
    class TraceHop:
        def integrate(self, seed, power=None):
            power = float(power)
            trace = 1.e-9 * (1. + power + (0. if power < .3 else .5))
            u = 10. + power
            return _complete({
                'u': u,
                'n_e': 1.e16,
                '_state_vector': [trace, u],
                '_jacobian': np.diag([-1., -2.]).tolist(),
            })

    path = continue_branch(
        TraceHop(), {'u': 10., 'n_e': 1.e16, 'state': [1.e-9, 10.]},
        [.5], initial_power=0.)

    event = path[-1]
    assert event['event_type'] == 'unresolved'
    assert event['event_bracket'][0] < .3 <= event['event_bracket'][1]
    assert all(entry['state'][0] < 1.e-9 * (1.25 + entry['power'])
               for entry in path if entry.get('converged'))


def test_eigenvalue_motion_bound_refines_toward_a_narrow_excursion():
    """An endpoint in the tail of an unstable excursion forces refinement."""
    class TailCrossing:
        def integrate(self, seed, power=None):
            power = float(power)
            eigenvalue = -1. + 2. * np.exp(-((power - .0228) / .002) ** 2)
            return _complete({
                'u': 0.,
                'n_e': 1.e16,
                '_state_vector': [0.],
                '_jacobian': [[eigenvalue]],
            })

    path = continue_branch(
        TailCrossing(), {'u': 0., 'n_e': 1.e16, 'state': [0.]},
        [.025], initial_power=0., max_step=.025, min_step=.00025,
        stability_resolution=.025)

    event = path[-1]
    crossing = .0228 - np.sqrt(np.log(2.)) * .002
    assert event['event_type'] == 'stability'
    assert event['event_bracket'][0] <= crossing <= event['event_bracket'][1]


def test_fast_companion_component_cannot_dilute_a_branch_jump():
    """The reviewer's v = (u - 2) / 0.02 cubic beside a fast second species."""
    class CubicWithCompanion:
        def integrate(self, seed, power=None):
            power = float(power)
            roots = np.roots([1., 0., -1., -power])
            states = [2. + .02 * root.real
                      for root in roots if abs(root.imag) < 1.e-10]
            u = min(states, key=lambda value: abs(value - seed.state[0]))
            physical = (u - 2.) / .02
            return _complete({
                'u': u,
                'n_e': 1.e16,
                '_state_vector': [u, 1. + power],
                '_jacobian': [[1. - 3. * physical * physical, 0.], [0., -1.]],
            })

    start = 2. - .02
    path = continue_branch(
        CubicWithCompanion(), {'u': start, 'n_e': 1.e16, 'state': [start, 1.]},
        [.5], initial_power=0.)

    event = path[-1]
    exact_fold = 2. / np.sqrt(27.)
    assert event['event_type'] == 'fold'
    assert event['fold_bracket'][0] < exact_fold < event['fold_bracket'][1]
    assert not any(entry.get('power', 0.) > event['fold_bracket'][1]
                   for entry in path)


def test_tangent_reversal_is_refined_and_reported_unresolved_not_fold():
    """A backward hop reverses the secant; that alone is not fold evidence."""
    class BackwardHop:
        def integrate(self, seed, power=None):
            power = float(power)
            u = 10. + power - (0. if power < .3 else .2)
            return _complete({
                'u': u,
                'n_e': 1.e16,
                '_state_vector': [u],
                '_jacobian': [[-1.]],
            })

    path = continue_branch(
        BackwardHop(), {'u': 10., 'n_e': 1.e16, 'state': [10.]},
        [.5], initial_power=0.)

    event = path[-1]
    assert event['event_type'] == 'unresolved'
    assert event['failure'] == 'continuation tangent reversed orientation'
    assert event['fold_detected'] is False
    assert event['event_bracket'][0] < .3 <= event['event_bracket'][1]


def test_continuation_never_steps_below_the_declared_floor():
    """A tail shorter than the minimum step joins the step before it."""
    class LinearStable:
        def integrate(self, seed, power=None):
            power = float(power)
            return _complete({
                'u': 1. + power,
                'n_e': 1.e16,
                '_state_vector': [1. + power],
                '_jacobian': [[-1.]],
            })

    path = continue_branch(
        LinearStable(), {'u': 1.5, 'n_e': 1.e16, 'state': [1.5]},
        [.5125 + 1.e-9], initial_power=.5, max_step=.0125, min_step=.001)

    assert path[-1]['power'] == .5125 + 1.e-9
    assert path[-1]['converged'] is True
    assert min(entry['step_size'] for entry in path[1:]) >= .001


def test_corrector_noise_is_not_read_as_a_tangent_reversal():
    """Steps inside the measured corrector scatter carry no orientation."""
    class NoisyFlat:
        def __init__(self):
            self.calls = 0

        def integrate(self, seed, power=None):
            self.calls += 1
            u = 2. + 1.e-9 * np.sin(1.7 * self.calls)
            return _complete({
                'u': u,
                'n_e': 1.e16,
                '_state_vector': [u],
                '_jacobian': [[-1.]],
            })

    path = continue_branch(
        NoisyFlat(), {'u': 2., 'n_e': 1.e16, 'state': [2.]},
        [.1], initial_power=0., max_step=.025, min_step=.005,
        stability_resolution=.025)

    assert path[0]['corrector_reproducibility'] == pytest.approx(
        1.e-9 * abs(np.sin(3.4) - np.sin(1.7)) / 2., rel=1.e-6)
    assert path[-1]['power'] == pytest.approx(.1)
    assert all(entry.get('event_type') is None for entry in path)
    assert all(entry['tangent_cosine'] is None for entry in path)


def test_hopf_crossing_is_labeled_hopf_not_fold():
    class HopfCrossing:
        def integrate(self, seed, power=None):
            real = float(power) - .05
            return _complete({
                'u': 0.,
                'n_e': 1.e16,
                '_state_vector': [0., 0.],
                '_jacobian': [[real, -2.], [2., real]],
            })

    path = continue_branch(
        HopfCrossing(), {'u': 0., 'n_e': 1.e16, 'state': [0., 0.]},
        [.1], initial_power=0., max_step=.1, min_step=.001)

    event = path[-1]
    assert event['event_type'] == 'hopf'
    assert event['hopf_detected'] is True
    assert event['fold_detected'] is False
    assert event['event_bracket'][0] <= .05 <= event['event_bracket'][1]


def test_worker_timeout_and_post_result_shutdown_are_deadline_bounded(tmp_path):
    pid_directory = tmp_path / 'blocking-pids'
    pid_directory.mkdir()
    elapsed, attempts = _watchdog_scan(
        _run_blocking_scan,
        (tmp_path / 'blocking.json', pid_directory),
        [pid_directory / 'worker-0.pid', pid_directory / 'worker-1.pid'])

    assert elapsed < .7
    assert all('TimeoutError' in attempt['failure'] for attempt in attempts)

    post_result_pid = tmp_path / 'post-result.pid'
    elapsed, attempts = _watchdog_scan(
        _run_post_result_scan,
        (tmp_path / 'post-result.json', post_result_pid),
        [post_result_pid])

    assert elapsed < .7
    assert attempts[0]['converged'] is True


def test_continuation_carries_the_full_state_between_corrector_steps():
    """A metastable population, not only (u, n_e), selects this branch."""
    class TwoState:
        def __init__(self):
            self.corrector_seeds = []

        def integrate(self, seed, power=None):
            self.corrector_seeds.append(tuple(seed.state))
            if seed.state[0] > 0.5:
                u, metastable = -1. + float(power), 1. - float(power)
            else:
                u, metastable = 1. + float(power), 0.
            return _complete({
                'u': u,
                'n_e': 1.e16,
                'self_sustained': True,
                'discharge_class': 'self-sustained',
                '_state_vector': [metastable, u],
                '_jacobian': np.diag([-1., -2.]).tolist(),
            })

    adapter = TwoState()
    path = continue_branch(
        adapter, {'u': -1., 'n_e': 1.e16, 'state': [1., -1.]}, [0.2],
        initial_power=0., max_step=0.05)

    assert path[-1]['power'] == pytest.approx(0.2)
    assert all(entry['u'] < 0. for entry in path if entry.get('converged'))
    assert all(seed[0] > 0.5 for seed in adapter.corrector_seeds)
    assert path[-1]['state'] == pytest.approx([0.8, -0.8])


def test_crashed_worker_is_recorded_and_sleeping_sibling_is_reaped(tmp_path):
    """The crash is recorded however loaded the machine is.

    The crashing worker exits at once and the sibling sleeps until it has
    exited, so the only timing left is the 2 s deadline, far longer than a
    fork and an exit. The crash is never relabelled a timeout.
    """
    sibling_pid = tmp_path / 'sibling.pid'
    crashed = tmp_path / 'crashed'

    class CrashAndSleep:
        def integrate(self, seed, power=None):
            if seed.u == 0.:
                crashed.write_text('exiting')
                os._exit(7)
            sibling_pid.write_text(str(os.getpid()))
            while not crashed.exists():
                time.sleep(.01)
            time.sleep(60.)

    result = scan(
        CrashAndSleep(), (0., 1.), tmp_path / 'branches.json', u_steps=2,
        electron_density_decades=1, processes=2, timeout=2., reference_power=.5,
        recheck_timeouts=False)

    assert len(result['attempts']) == 2
    assert 'exit code 7' in result['attempts'][0]['failure']
    assert 'TimeoutError' not in result['attempts'][0]['failure']
    assert 'TimeoutError' in result['attempts'][1]['failure']
    pid = int(sibling_pid.read_text())
    with pytest.raises(ProcessLookupError):
        os.kill(pid, 0)


def test_matching_branch_uses_artifact_tolerances_and_refuses_nonunique_matches(tmp_path):
    path = tmp_path / 'branches.json'
    payload = {
        'u_tolerance': 1.e-6,
        'log_ne_tolerance': 1.e-6,
        'branches': [
            {'id': 'reactor-000', 'u': 1., 'n_e': 1.e16},
            {'id': 'reactor-001', 'u': 1.00005, 'n_e': 1.e16},
        ],
    }
    path.write_text(json.dumps(payload))

    assert matching_branch(path, {'u': 1.00005, 'n_e': 1.e16}) == 'reactor-001'
    with pytest.raises(ValueError, match='no recorded branch'):
        matching_branch(path, {'u': 2., 'n_e': 1.e16})

    payload['u_tolerance'] = 1.e-4
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match='more than one recorded branch'):
        matching_branch(path, {'u': 1.000025, 'n_e': 1.e16})


def test_matching_branch_uses_only_certified_path_at_run_power(tmp_path):
    path = tmp_path / 'branches.json'
    branch = {
        'id': 'reactor-000',
        'u': 1.5,
        'n_e': 2.e16,
        'continuation_path': [
            {'power': .25, 'u': 1., 'n_e': 1.e16, 'converged': True,
             'segment_certified': True},
            {'power': .5, 'u': 1.5, 'n_e': 2.e16, 'converged': True,
             'segment_certified': True,
             'interpolation_segments': [{
                 'certificate': 'declared-error-bound-v1',
                 'power': [.25, .5],
                 'u': [1., 1.5],
                 'n_e': [1.e16, 2.e16],
                 'u_error_bound': 0.,
                 'log_ne_error_bound': 0.,
                 'resolution': .125,
             }]},
        ],
    }
    # A matchable branch needs a recorded half-step path; here it holds the
    # same states, so it agrees everywhere.
    branch['continuation_path_half_step'] = branch['continuation_path']
    payload = {
        'reference_power': .5,
        'u_tolerance': 1.e-6,
        'log_ne_tolerance': 1.e-6,
        'branches': [branch],
    }
    path.write_text(json.dumps(payload))

    assert matching_branch(
        path, {'u': 1., 'n_e': 1.e16}, power=.25) == 'reactor-000'
    expected_ne = np.exp((np.log(1.e16) + np.log(2.e16)) / 2.)
    assert matching_branch(
        path, {'u': 1.25, 'n_e': expected_ne}, power=.375) == 'reactor-000'

    branch['continuation_path'][1]['segment_certified'] = False
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match='power not certified for branch reactor-000'):
        matching_branch(
            path, {'u': 1.25, 'n_e': expected_ne}, power=.375)


def test_matching_branch_uses_recorded_certification_points_on_curved_path(tmp_path):
    path = tmp_path / 'branches.json'
    branch = {
        'id': 'reactor-000',
        'u': .25,
        'n_e': 1.e16,
        'continuation_path': [
            {'power': .25, 'u': .25 ** 2, 'n_e': 1.e16,
             'converged': True, 'segment_certified': True},
            {'power': .5, 'u': .5 ** 2, 'n_e': 1.e16,
             'converged': True, 'segment_certified': True,
             'certification_points': [
                 {'power': .3995, 'u': .3995 ** 2, 'n_e': 1.e16},
                 {'power': .4005, 'u': .4005 ** 2, 'n_e': 1.e16},
             ],
             'interpolation_segments': [{
                 'certificate': 'declared-error-bound-v1',
                 'power': [.3995, .4005],
                 'u': [.3995 ** 2, .4005 ** 2],
                 'n_e': [1.e16, 1.e16],
                 'u_error_bound': 2.5e-7,
                 'log_ne_error_bound': 0.,
                 'resolution': .0005,
             }]},
        ],
    }
    branch['continuation_path_half_step'] = branch['continuation_path']
    payload = {
        'reference_power': .5,
        'u_tolerance': 1.e-6,
        'log_ne_tolerance': 1.e-6,
        'branches': [branch],
    }
    path.write_text(json.dumps(payload))

    assert matching_branch(
        path, {'u': .4 ** 2, 'n_e': 1.e16}, power=.4
    ) == 'reactor-000'


def test_scanner_emits_interpolation_only_from_adapter_declared_bound(tmp_path):
    class CertifiedLinear:
        def integrate(self, seed, power=None):
            power = float(power)
            return _complete({
                'u': power,
                'n_e': 1.e16,
                '_state_vector': [power],
                '_jacobian': [[-1.]],
            })

        @staticmethod
        def certify_interpolation(left, middle, right):
            return {'u_error_bound': 0., 'log_ne_error_bound': 0.}

    continuation = continue_branch(
        CertifiedLinear(), {'u': .25, 'n_e': 1.e16, 'state': [.25]},
        [.5], initial_power=.25, max_step=.25, min_step=.01,
        stability_resolution=.05, u_tolerance=1.e-6,
        log_ne_tolerance=1.e-6)
    half_step = continue_branch(
        CertifiedLinear(), {'u': .25, 'n_e': 1.e16, 'state': [.25]},
        [.5], initial_power=.25, max_step=.125, min_step=.005,
        stability_resolution=.025, u_tolerance=1.e-6,
        log_ne_tolerance=1.e-6)
    artifact = {
        'reference_power': .5,
        'u_tolerance': 1.e-6,
        'log_ne_tolerance': 1.e-6,
        'branches': [{
            'id': 'reactor-000', 'u': .5, 'n_e': 1.e16,
            'continuation_path': continuation,
            'continuation_path_half_step': half_step,
        }],
    }
    path = tmp_path / 'branches.json'
    path.write_text(json.dumps(artifact))

    assert matching_branch(
        path, {'u': .4, 'n_e': 1.e16}, power=.4
    ) == 'reactor-000'

    class LooseCertificate(CertifiedLinear):
        @staticmethod
        def certify_interpolation(left, middle, right):
            return {'u_error_bound': 1.e-3, 'log_ne_error_bound': 0.}

    loose = continue_branch(
        LooseCertificate(), {'u': .25, 'n_e': 1.e16, 'state': [.25]},
        [.5], initial_power=.25, max_step=.25, min_step=.01,
        stability_resolution=.05, u_tolerance=1.e-6,
        log_ne_tolerance=1.e-6)
    assert loose[-1]['converged'] is True
    assert all(not entry.get('interpolation_segments') for entry in loose)


@pytest.mark.parametrize('charge_row_scale', [1., 4.])
def test_algebraic_electron_stability_is_invariant_to_charge_row_scale(charge_row_scale):
    """The analytic DAE has physical eigenvalues [-1, -1] at either scale."""
    jacobian = np.array([
        [-charge_row_scale, charge_row_scale, 0.],
        [-3., 2., 0.],
        [0., 0., -1.],
    ])

    result = stability_analysis(
        jacobian, [], half_step_jacobian=jacobian, algebraic_indices=[0])

    assert result['classification'] == 'stable'
    assert sorted(value['real'] for value in result['eigenvalues']) == pytest.approx([-1., -1.])
    assert result['algebraic_indices'] == [0]


def test_branch_scanner_has_the_project_mit_header():
    source = Path(__file__).parents[4] / 'rmgpy/tools/eedf/branches.py'
    prefix = source.read_text().splitlines()[:28]

    assert "Permission is hereby granted, free of charge" in '\n'.join(prefix)


# ---------------------------------------------------------------------------
# The reversible-chain certificate (review r220) and its counterexamples.

REAL_SETTINGS = {'initial_power': .5, 'max_step': .025, 'min_step': .00625,
                 'stability_resolution': .00625}
REAL_TARGETS = (.25, .5, 1.)


def _accepted_states(path):
    """Every accepted endpoint and recorded certification point of a path."""
    states = []
    for entry in path:
        if not entry.get('converged'):
            continue
        states.extend(entry.get('certification_points', []))
        states.append(entry)
    return states


class _ParallelPair:
    """Stable states z = 0 (A) and z = gap (B) of y' = -k z (z - gap/2)(z - gap).

    z = y - a(P). A time-integrating corrector settles in its seed's basin;
    the saddle sits at z = gap / 2. Both states linearise to -1, so no
    eigenvalue tells them apart.
    """

    def __init__(self, a, gap):
        self.a = a
        self.gap = gap

    def branch(self, u, power):
        return 'A' if abs(u - self.a(power)) < abs(self.gap) / 2. else 'B'

    def integrate(self, seed, power=None):
        power = float(power)
        z = seed.state[0] - self.a(power)
        settled = 0. if z * np.sign(self.gap) < abs(self.gap) / 2. else self.gap
        u = self.a(power) + settled
        return _complete({'u': u, 'n_e': 1.e16, '_state_vector': [u],
                          '_jacobian': [[-1.]]})


class _TimeCubic:
    """v' = P + v - v**3 integrated in time; u = offset + shift(P) + scale * v.

    A seed below the middle root settles on the lower root, above it on the
    upper root. Past the fold only the upper root exists.
    """

    def __init__(self, scale=.02, offset=2., stretch=0., slow_window=None):
        self.scale = scale
        self.offset = offset
        self.stretch = stretch
        self.slow_window = slow_window

    def shift(self, power):
        return self.stretch * (1. + np.tanh((power - .1) / .002))

    def v(self, u, power):
        return (u - self.offset - self.shift(power)) / self.scale

    def integrate(self, seed, power=None):
        power = float(power)
        if self.slow_window and self.slow_window[0] < power < self.slow_window[1]:
            raise RuntimeError('corrector time budget exhausted')
        roots = sorted(root.real for root in np.roots([1., 0., -1., -power])
                       if abs(root.imag) < 1.e-10)
        start = self.v(seed.state[0], power)
        v = roots[0] if len(roots) == 1 or start < roots[1] else roots[-1]
        u = self.offset + self.shift(power) + self.scale * v
        return _complete({'u': u, 'n_e': 1.e16, '_state_vector': [u],
                          '_jacobian': [[1. - 3. * v * v]]})

    @staticmethod
    def is_fold_failure(failure):
        return failure == 'RuntimeError: reactor returned without reaching steady state'


class _SilentHop:
    """A corrector that ignores its seed and jumps by 0.2 at P = 0.3."""

    def __init__(self, fast_stretch=False, scatter=0.):
        self.fast_stretch = fast_stretch
        self.scatter = scatter
        self.calls = 0

    def u(self, power):
        u = 10. + power + (0. if power < .3 else .2)
        if self.fast_stretch:
            u += 5. * (1. + np.tanh((power - .1) / .002))
        return u

    def integrate(self, seed, power=None):
        self.calls += 1
        u = self.u(float(power))
        if self.calls <= 2:
            u *= 1. + self.scatter * (1. if self.calls == 1 else -1.)
        return _complete({'u': u, 'n_e': 1.e16, '_state_vector': [u],
                          '_jacobian': [[-1.]]})


@pytest.mark.parametrize('max_step', [None, .0125])
def test_r220_case_a_parallel_branch_hop_is_never_accepted(max_step):
    """Review r220 case A: branches 0.004 apart (40 u_tol), identical eigenvalues."""
    adapter = _ParallelPair(lambda power: 10. - 10. * power ** 3, .004)

    path = continue_branch(
        adapter, {'u': 10., 'n_e': 1.e16, 'state': [10.]}, [.5],
        initial_power=0., max_step=max_step)

    accepted = _accepted_states(path)
    assert len(accepted) > 10
    assert all(adapter.branch(entry['u'], entry['power']) == 'A' for entry in accepted)
    assert path[-1]['event_type'] is not None or path[-1]['power'] == pytest.approx(.5)


def test_r220_case_b_earlier_fast_stretch_cannot_loosen_a_later_step():
    """Review r220 case B: a smooth fast stretch at 0.1 must not hide the 0.3 hop."""
    adapter = _SilentHop(fast_stretch=True)
    start = adapter.u(0.)

    path = continue_branch(
        adapter, {'u': start, 'n_e': 1.e16, 'state': [start]}, [.5],
        initial_power=0.)

    event = path[-1]
    assert event['event_type'] == 'unresolved'
    assert event['event_bracket'][0] < .3 <= event['event_bracket'][1]
    assert not any(entry['power'] >= .3 for entry in _accepted_states(path))
    assert any(entry['power'] > .11 for entry in _accepted_states(path))


@pytest.mark.parametrize('scatter', [2.e-3, 4.e-6])
def test_r220_case_c_noisy_initial_resolve_cannot_widen_the_noise_floor(scatter):
    """Review r220 case C: initial scatter must reproduce within the branch tolerance."""
    adapter = _SilentHop(scatter=scatter)

    path = continue_branch(
        adapter, {'u': 10., 'n_e': 1.e16, 'state': [10.]}, [.5],
        initial_power=0.)

    event = path[-1]
    assert event['event_type'] == 'unresolved'
    if scatter > 1.e-5:
        # 2e-3 of u = 10 is 200 u_tol: the corrector cannot certify identity.
        assert len(path) == 1
        assert 'does not reproduce the initial state' in event['failure']
    else:
        assert event['event_bracket'][0] < .3 <= event['event_bracket'][1]
        assert not any(entry['power'] >= .3 for entry in _accepted_states(path))


@pytest.mark.parametrize('max_step', [None, .0125])
def test_r220_case_f_smooth_shift_cannot_carry_the_path_onto_the_upper_branch(max_step):
    """Review r220 case F: the rescaled cubic plus a smooth shift at P = 0.1."""
    adapter = _TimeCubic(stretch=.5)
    start = adapter.offset + adapter.shift(0.) - adapter.scale
    fold = 2. / np.sqrt(27.)

    path = continue_branch(
        adapter, {'u': start, 'n_e': 1.e16, 'state': [start]}, [.5],
        initial_power=0., max_step=max_step)

    accepted = _accepted_states(path)
    assert accepted
    assert all(adapter.v(entry['u'], entry['power']) < -1. / np.sqrt(3.)
               for entry in accepted)
    assert all(entry['power'] < fold for entry in accepted)
    assert path[-1]['event_type'] in ('fold', 'unresolved')


@pytest.mark.parametrize('gap', [-.01, -.004])
def test_r220_case_g_first_step_of_a_leg_is_certified_at_real_settings(gap):
    """Review r220 case G: a parallel branch beside the real deck's slope."""
    adapter = _ParallelPair(lambda power: 1.75 - .86 * (power - .5), gap)

    path = continue_branch(
        adapter, {'u': 1.75, 'n_e': 1.e16, 'state': [1.75]}, REAL_TARGETS,
        **REAL_SETTINGS)

    accepted = _accepted_states(path)
    assert all(adapter.branch(entry['u'], entry['power']) == 'A' for entry in accepted)
    assert path[-1]['event_type'] is not None or path[-1]['power'] == pytest.approx(1.)
    branch = {'id': 'reactor-000', 'u': 1.75, 'n_e': 1.e16, 'continuation_path': path}
    consistency = path_consistency(branch, 1.e-4, 1.e-3)
    assert consistency['consistent'] is False
    assert consistency['disagreeing_powers'] == []
    assert consistency['half_step_recorded'] is False


class _RotatingPair:
    """Stable states c +- (cos wP, sin wP); a seed settles on the one it faces."""

    centre = np.array([3., 3.])

    def __init__(self, omega):
        self.omega = omega

    def direction(self, power):
        angle = self.omega * float(power)
        return np.array([np.cos(angle), np.sin(angle)])

    def faces(self, state, power):
        return float(np.dot(np.asarray(state) - self.centre, self.direction(power))) > 0.

    def integrate(self, seed, power=None):
        axis = self.direction(power)
        state = self.centre + (axis if self.faces(seed.state, power) else -axis)
        rotation = np.array([[axis[0], -axis[1]], [axis[1], axis[0]]])
        jacobian = rotation @ np.diag([-2., -1.]) @ rotation.T
        return _complete({'u': float(state[0]), 'n_e': 1.e16 * float(np.exp(state[1])),
                          '_state_vector': state.tolist(),
                          '_jacobian': jacobian.tolist()})


def test_certificate_sees_branches_that_exchange_places_within_one_step():
    """Own counterexample 1: two stable branches swap basins between P0 and P1.

    Over one 0.1 W step the pair turns by 0.9 pi. A zero-order solve over the
    whole step lands on the other branch, and solving back returns to the
    start, so a certificate that checks only the two endpoints is fooled.
    The chain's 0.025 W gaps each turn by 0.225 pi and keep the branch.
    """
    adapter = _RotatingPair(9. * np.pi)
    start = adapter.centre + adapter.direction(0.)
    seed = {'u': float(start[0]), 'n_e': 1.e16 * float(np.exp(start[1])),
            'state': start.tolist()}

    forward = adapter.integrate(Seed(seed['u'], seed['n_e'], tuple(start)), power=.1)
    assert not adapter.faces(forward['_state_vector'], .1)
    back = adapter.integrate(
        Seed(forward['u'], forward['n_e'], tuple(forward['_state_vector'])), power=0.)
    assert back['_state_vector'] == pytest.approx(start.tolist())

    path = continue_branch(
        adapter, seed, [.5], initial_power=0., max_step=.1, min_step=.00625,
        stability_resolution=.025)

    accepted = _accepted_states(path)
    assert len(accepted) > 20
    assert all(adapter.faces(entry['state'], entry['power']) for entry in accepted)
    assert path[-1]['event_type'] is not None or path[-1]['power'] == pytest.approx(.5)


def test_aliasing_at_the_last_gap_is_exposed_by_step_halving_and_refused(tmp_path):
    """Own counterexample 2, the declared limit: a hop reversible at its own gap.

    Branch A is u = 1 - P**3 and B sits 0.109 above it, so the basin boundary
    is 0.0547 above A. The step path takes 0 -> 0.5 W in 0.125 W gaps; over
    the last gap A falls by 0.0723, past the boundary, and solving back from
    B lands on A again: that gap aliases. The half-step path's gaps move A by
    at most 0.0413 and keep it, so the two paths disagree at 0.5 W, and
    matching refuses there.
    """
    gap = 2. * 28. / 512.
    adapter = _ParallelPair(lambda power: 1. - power ** 3, gap)
    seed = {'u': 1., 'n_e': 1.e16, 'state': [1.]}

    main = continue_branch(adapter, seed, [.5], initial_power=0., max_step=.5,
                           min_step=.125, stability_resolution=.125)
    half = continue_branch(adapter, seed, [.5], initial_power=0., max_step=.25,
                           min_step=.015625, stability_resolution=.0625)

    assert [entry['power'] for entry in main] == [0., .5]
    assert adapter.branch(main[-1]['u'], .5) == 'B'
    assert all(adapter.branch(entry['u'], entry['power']) == 'A'
               for entry in _accepted_states(main)[:-1])
    assert all(adapter.branch(entry['u'], entry['power']) == 'A'
               for entry in _accepted_states(half))
    assert half[-1]['power'] == .5 and half[-1]['converged'] is True

    payload = {
        'reference_power': 0.,
        'u_tolerance': 1.e-4,
        'log_ne_tolerance': 1.e-3,
        'branches': [{'id': 'reactor-000', 'u': 1., 'n_e': 1.e16,
                      'continuation_path': main,
                      'continuation_path_half_step': half}],
    }
    record_continuation(payload)
    assert payload['path_consistency']['reactor-000']['consistent'] is False
    assert payload['path_consistency']['reactor-000']['disagreeing_powers'] == [.5]
    path = tmp_path / 'branches.json'
    path.write_text(json.dumps(payload))

    with pytest.raises(ValueError, match='power not certified for branch reactor-000'):
        matching_branch(path, {'u': main[-1]['u'], 'n_e': 1.e16}, power=.5)
    assert matching_branch(path, {'u': 1., 'n_e': 1.e16}, power=0.) == 'reactor-000'

    # A half-step path that stops early leaves the later step states unconfirmed.
    payload['branches'][0]['continuation_path_half_step'] = half[:2]
    record_continuation(payload)
    assert payload['path_consistency']['reactor-000']['unconfirmed_by_half_step'][-1] == .5
    assert payload['path_consistency']['reactor-000']['consistent'] is False


def test_marginal_sample_stops_the_path_with_an_event():
    """Review r220 P2-1: a sample whose own verdict is marginal is not continued."""
    class MarginalPatch:
        def integrate(self, seed, power=None):
            power = float(power)
            half = 1.e-3 if .2 < power < .3 else -1.
            return _complete({'u': 1. + power, 'n_e': 1.e16,
                              '_state_vector': [1. + power], '_jacobian': [[-1.]],
                              '_half_step_jacobian': [[half]]})

    path = continue_branch(
        MarginalPatch(), {'u': 1., 'n_e': 1.e16, 'state': [1.]}, [.5],
        initial_power=0.)

    event = path[-1]
    assert event['event_type'] == 'unresolved'
    assert "verdict 'marginal'" in event['failure']
    assert event['event_bracket'][0] <= .2 < event['event_bracket'][1]
    assert all(entry['stability']['stable'] for entry in path if entry.get('converged'))
    assert not any(entry['power'] > .2 for entry in _accepted_states(path))


def test_continuation_refuses_a_non_stable_initial_state():
    class MarginalStart:
        def integrate(self, seed, power=None):
            return _complete({'u': 1., 'n_e': 1.e16, '_state_vector': [1.],
                              '_jacobian': [[-1.]], '_half_step_jacobian': [[1.e-3]]})

    path = continue_branch(
        MarginalStart(), {'u': 1., 'n_e': 1.e16, 'state': [1.]}, [.5],
        initial_power=0.)

    assert len(path) == 1
    assert path[0]['event_type'] == 'unresolved'
    assert path[0]['failure'] == 'initial state is marginal; a non-stable state is not continued'


class _LinearStable:
    def integrate(self, seed, power=None):
        power = float(power)
        return _complete({'u': 1. + power, 'n_e': 1.e16 * (1. + power),
                          '_state_vector': [1. + power], '_jacobian': [[-1.]]})


def test_off_grid_power_is_refused_as_not_certified_with_certified_powers(tmp_path):
    """Review r220 P2-2: no certified state at the run power is named as such."""
    tolerances = {'u_tolerance': 1.e-6, 'log_ne_tolerance': 1.e-6}
    path = continue_branch(
        _LinearStable(), {'u': 1.5, 'n_e': 1.5e16, 'state': [1.5]}, [.25, .5],
        initial_power=.5, max_step=.125, min_step=.03125,
        stability_resolution=.0625, **tolerances)
    half_step = continue_branch(
        _LinearStable(), {'u': 1.5, 'n_e': 1.5e16, 'state': [1.5]}, [.25, .5],
        initial_power=.5, max_step=.0625, min_step=.015625,
        stability_resolution=.03125, **tolerances)
    payload = {'reference_power': .5, 'u_tolerance': 1.e-6, 'log_ne_tolerance': 1.e-6,
               'branches': [{'id': 'reactor-000', 'u': 1.5, 'n_e': 1.5e16,
                             'continuation_path': path,
                             'continuation_path_half_step': half_step}]}
    record_continuation(payload)
    coverage = payload['certified_coverage']['reactor-000']
    assert coverage['powers'][0] == pytest.approx(.25)
    assert .3125 in coverage['powers']
    assert coverage['segments'] == []
    artifact = tmp_path / 'branches.json'
    artifact.write_text(json.dumps(payload))

    assert matching_branch(
        artifact, {'u': 1.3125, 'n_e': 1.3125e16}, power=.3125) == 'reactor-000'
    with pytest.raises(ValueError) as refused:
        matching_branch(artifact, {'u': 1.3123, 'n_e': 1.3123e16}, power=.3123)
    message = str(refused.value)
    assert message.startswith('power not certified for branch reactor-000 at 0.3123 W')
    assert ("certified powers: ['0.25', '0.28125', '0.3125', '0.34375', '0.375', "
            "'0.40625', '0.4375', '0.46875', '0.5']") in message
    assert 'no recorded branch' not in message


def test_artifact_without_reference_power_refuses_a_power_aware_match(tmp_path):
    path = tmp_path / 'branches.json'
    path.write_text(json.dumps({
        'u_tolerance': 1.e-6, 'log_ne_tolerance': 1.e-6,
        'branches': [{'id': 'reactor-000', 'u': 1., 'n_e': 1.e16}]}))

    with pytest.raises(ValueError, match='records no reference power'):
        matching_branch(path, {'u': 1., 'n_e': 1.e16}, power=.7)


def test_scan_records_certificate_resolution_coverage_and_consistency(tmp_path):
    """Review r220 P3-2 and P2-2: the artifact states what it certifies."""
    class Linear:
        def integrate(self, seed, power=None):
            power = float(power)
            return _complete({'u': 1. + power, 'n_e': 1.e16, '_state_vector': [1. + power],
                              '_jacobian': [[-1.]], 'self_sustained': True,
                              'discharge_class': 'self-sustained'})

    payload = scan(Linear(), (1., 1.), tmp_path / 'branches.json', u_steps=1,
                   electron_density_decades=1, powers=(.5, .25), reference_power=.5,
                   continuation_max_step=.125, continuation_stability_resolution=.0625)

    written = json.loads((tmp_path / 'branches.json').read_text())
    assert written['stability_resolution'] == .0625
    assert written['continuation_certificate']['name'] == 'reversible-chain-v1'
    assert 'reversible within the branch tolerances' in written['continuation_certificate']['statement']
    assert written['continuation_certificate']['limits']
    assert written['continuation_settings']['reactor-000']['continuation_path'][
        'stability_resolution'] == .0625
    assert written['certified_coverage']['reactor-000']['powers'][0] == .25
    assert written['path_consistency']['reactor-000']['consistent'] is False
    assert written['path_consistency']['reactor-000']['half_step_recorded'] is False
    assert written['continuation_solves']['reactor-000']['continuation_path'][
        'certificate_solves'] > 0
    assert payload['timeout_rechecks'] == []


def test_fold_bracket_is_unresolved_when_a_non_fold_failure_interrupts_it():
    """Review r220 P3-1: a timeout near the fold cannot narrow the fold bracket."""
    fold = 2. / np.sqrt(27.)
    adapter = _TimeCubic(slow_window=(.379, fold))
    start = adapter.offset + adapter.scale * min(
        root.real for root in np.roots([1., 0., -1., -.35]) if abs(root.imag) < 1.e-10)

    path = continue_branch(
        adapter, {'u': start, 'n_e': 1.e16, 'state': [start]}, [.5],
        initial_power=.35, max_step=.025, min_step=1.e-4, stability_resolution=.025)

    event = path[-1]
    assert event['event_type'] == 'unresolved'
    assert event['failure'].startswith(
        'corrector failed while bracketing a fold: RuntimeError: corrector time budget')
    assert event['event_bracket'][0] < fold < event['event_bracket'][1]
    assert event['fold_bracket'] is None


def test_fast_mode_stiffening_far_from_zero_is_not_refused_at_real_settings():
    """Review r220 P3-3 case I: only modes near zero carry the motion bound."""
    class Stiffening:
        def integrate(self, seed, power=None):
            power = float(power)
            fast = -1.e5 * (1. + 5. * (1. + np.tanh((power - .7) / .002)))
            u = 1.75 - .86 * (power - .5)
            return _complete({'u': u, 'n_e': 1.e16, '_state_vector': [u, 1.],
                              '_jacobian': np.diag([-483., fast]).tolist()})

    path = continue_branch(
        Stiffening(), {'u': 1.75, 'n_e': 1.e16, 'state': [1.75, 1.]}, REAL_TARGETS,
        **REAL_SETTINGS)

    assert path[-1]['power'] == 1.
    assert path[-1]['converged'] is True
    assert all(entry.get('event_type') is None for entry in path)


@pytest.mark.parametrize('missing, message', [
    ('_state_vector', 'no full _state_vector'),
    ('_half_step_jacobian', 'no half-step Jacobian'),
    ('_conservation_constraints', 'no conservation constraints'),
])
def test_adapter_missing_data_is_refused_not_defaulted(missing, message, tmp_path):
    """Review r220 P3-4: an incomplete adapter stops both scan and continuation."""
    class Incomplete:
        def integrate(self, seed, power=None):
            state = _complete({'u': 1., 'n_e': 1.e16, '_state_vector': [1.],
                               '_jacobian': [[-1.]], 'self_sustained': True,
                               'discharge_class': 'self-sustained'})
            del state[missing]
            return state

    with pytest.raises(AdapterContractError, match=message):
        continue_branch(Incomplete(), {'u': 1., 'n_e': 1.e16, 'state': [1.]}, [.5],
                        initial_power=0.)
    with pytest.raises(AdapterContractError, match=message):
        scan(Incomplete(), (1., 1.), tmp_path / 'branches.json', u_steps=1,
             electron_density_decades=1, powers=(.5, .25), reference_power=.5,
             continuation_max_step=.125, continuation_stability_resolution=.0625)


def test_scan_and_continuation_refuse_missing_reference_power_or_seed_state(tmp_path):
    with pytest.raises(ValueError, match='finite reference_power'):
        scan(_LinearStable(), (1., 1.), tmp_path / 'branches.json', u_steps=1,
             electron_density_decades=1)
    with pytest.raises(AdapterContractError, match='no full state'):
        continue_branch(_LinearStable(), {'u': 1., 'n_e': 1.e16}, [.5], initial_power=0.)


def test_predictor_density_is_extrapolated_in_log_space_not_substituted():
    """Review r220 P3-4: a falling n_e is predicted, never replaced by the last value."""
    class FallingDensity:
        def __init__(self):
            self.seeds = []

        def integrate(self, seed, power=None):
            self.seeds.append((float(power), seed.n_e))
            power = float(power)
            return _complete({'u': 1. + power, 'n_e': 1.e16 * np.exp(-40. * power),
                              '_state_vector': [1. + power], '_jacobian': [[-1.]]})

    adapter = FallingDensity()
    path = continue_branch(
        adapter, {'u': 1., 'n_e': 1.e16, 'state': [1.]}, [.5], initial_power=0.,
        max_step=.05, min_step=.025, stability_resolution=.05)

    assert path[-1]['power'] == pytest.approx(.5)
    predictor_seeds = [entry['seed'] for entry in path[2:]]
    assert predictor_seeds
    for entry, seed in zip(path[2:], predictor_seeds):
        assert seed['n_e'] == pytest.approx(1.e16 * np.exp(-40. * entry['power']), rel=1.e-9)


def test_timed_out_seed_is_rechecked_and_recorded_in_the_artifact(tmp_path):
    """Review r220 P3-5: the recheck of a timed-out seed lives in branches.json."""
    marker = tmp_path / 'first-attempt'

    class SlowOnce:
        def integrate(self, seed, power=None):
            if not marker.exists():
                marker.write_text('slow')
                time.sleep(30.)
            return _complete({'u': 1., 'n_e': 1.e16, '_state_vector': [1.],
                              '_jacobian': [[-1.]]})

    payload = scan(SlowOnce(), (1., 1.), tmp_path / 'branches.json', u_steps=1,
                   electron_density_decades=1, timeout=2., reference_power=.5)

    attempt = json.loads((tmp_path / 'branches.json').read_text())['attempts'][0]
    assert 'TimeoutError' in attempt['first_attempt_failure']
    assert attempt['recheck'] == {'converged': True, 'failure': None}
    assert attempt['converged'] is True
    assert payload['timeout_rechecks'] == [0]
    assert len(payload['branches']) == 1


def test_steps_never_exceed_max_step_and_a_step_plus_floor_tail_is_split():
    """A remainder of step + floor is split, so a refined step is never snapped back."""
    path = continue_branch(
        _LinearStable(), {'u': 1., 'n_e': 1.e16, 'state': [1.]}, [.5],
        initial_power=0., max_step=.2, min_step=.1, stability_resolution=.1)

    sizes = [entry['step_size'] for entry in path[1:]]
    assert sizes == pytest.approx([.2, .2, .1])
    assert max(sizes) <= .2 * (1. + 1.e-9)


def test_eigen_event_bracket_is_unresolved_when_a_probe_fails():
    """Review r220 P3-1: a failed bisection probe cannot narrow a crossing bracket.

    The crossing is at 0.045 W and the probe at 0.04375 W times out. Treating
    the failure as the crossing's side would report [., 0.04375], which
    excludes the crossing.
    """
    class Crossing:
        def integrate(self, seed, power=None):
            power = float(power)
            if .041 < power < .044:
                raise RuntimeError('corrector time budget exhausted')
            return _complete({'u': power, 'n_e': 1.e16, '_state_vector': [power],
                              '_jacobian': [[power - .045]]})

    path = continue_branch(
        Crossing(), {'u': 0., 'n_e': 1.e16, 'state': [0.]}, [.1], initial_power=0.,
        max_step=.1, min_step=.001, stability_resolution=.05)

    event = path[-1]
    assert event['event_type'] == 'unresolved'
    assert event['failure'].startswith(
        'corrector failed while bracketing a stability event: RuntimeError')
    assert event['event_bracket'][0] <= .045 <= event['event_bracket'][1]


def test_per_component_departure_sees_a_jump_a_fast_companion_would_dilute():
    """A 10 % jump in one species beside a fast one stays a departure.

    The RMS corrector distance stays under the 0.25 predictor bound and the
    seed-independent jump is reversible, so only the per-component secant test
    sees it; an RMS departure test would divide it by the companion's motion.
    """
    class DilutedJump:
        def integrate(self, seed, power=None):
            power = float(power)
            species = 1.e-9 * (1. + power) * (1. if power < .3 else 1.1)
            fast = 1. + 50. * power
            return _complete({'u': 10. + power, 'n_e': 1.e16,
                              '_state_vector': [species, fast],
                              '_jacobian': np.diag([-1., -2.]).tolist()})

    path = continue_branch(
        DilutedJump(), {'u': 10., 'n_e': 1.e16, 'state': [1.e-9, 1.]}, [.5],
        initial_power=0.)

    event = path[-1]
    assert event['event_type'] == 'unresolved'
    assert event['event_bracket'][0] < .3 <= event['event_bracket'][1]
    assert all(entry['state'][0] < 1.e-9 * (1.05 + entry['power'])
               for entry in _accepted_states(path))


HALF_SETTINGS = {'initial_power': .5, 'max_step': .0125, 'min_step': .003125,
                 'stability_resolution': .003125}


def _with_paths(main, half, reference_power=.5, u=1.75, n_e=1.e16):
    payload = {'reference_power': reference_power, 'u_tolerance': 1.e-4,
               'log_ne_tolerance': 1.e-3,
               'branches': [{'id': 'reactor-000', 'u': u, 'n_e': n_e,
                             'continuation_path': main,
                             'continuation_path_half_step': half}]}
    return record_continuation(payload)


def test_limit_3_interior_single_gap_alias_is_declared_not_caught(tmp_path):
    """Review r221 case J: a hop confined to one gap passes both paths.

    A has the real deck's slope; B sits 0.02 (200 u_tol) above it. Both shift
    down by 0.012 through a tanh 2e-4 W wide inside one half-step gap, so that
    gap alone moves A by more than the 0.01 to the basin boundary. Every gap
    is reversible, both paths hop, and matching accepts B and refuses A. The
    recorded aliasing radius is the motion between recorded states; across the
    hop gap that is A's motion less the 0.02 offset, so the radius stays below
    the boundary distance. This is
    CERTIFICATE_LIMITS[2] as declared, not a detection.
    """
    neighbour, shift, width, centre = .02, .012, 2.e-4, .6015
    adapter = _ParallelPair(
        lambda power: 1.75 - .86 * (power - .5)
        - shift * (1. + np.tanh((power - centre) / width)) / 2., neighbour)
    seed = {'u': adapter.a(.5), 'n_e': 1.e16, 'state': [adapter.a(.5)]}

    main = continue_branch(adapter, seed, REAL_TARGETS, **REAL_SETTINGS)
    half = continue_branch(adapter, seed, REAL_TARGETS, **HALF_SETTINGS)

    for path in (main, half):
        assert path[-1]['power'] == 1. and path[-1]['event_type'] is None
        assert all(adapter.branch(entry['u'], entry['power']) == 'A'
                   for entry in _accepted_states(path) if entry['power'] < centre)
        assert all(adapter.branch(entry['u'], entry['power']) == 'B'
                   for entry in _accepted_states(path) if entry['power'] > .61)
    payload = _with_paths(main, half)
    assert payload['path_consistency']['reactor-000']['consistent'] is True
    radius = payload['aliasing_radius']['reactor-000']
    assert radius['continuation_path']['max_u_motion'] == pytest.approx(.86 * .00625, rel=1.e-6)
    for key in ('continuation_path', 'continuation_path_half_step'):
        assert radius[key]['max_u_motion'] < neighbour / 2.
    artifact = tmp_path / 'branches.json'
    artifact.write_text(json.dumps(payload))
    on_a = adapter.a(.7)
    assert matching_branch(artifact, {'u': on_a + neighbour, 'n_e': 1.e16}, power=.7) == 'reactor-000'
    with pytest.raises(ValueError, match='no recorded branch'):
        matching_branch(artifact, {'u': on_a, 'n_e': 1.e16}, power=.7)
    limit = CERTIFICATE_LIMITS[2]
    assert 'confined to one gap, the hop is undetected' in limit
    assert 'step halving does not expose it when the motion is narrower than the half-step gap' in limit
    assert 'the radius does not bound such a hop' in limit
    assert 'A steady neighbour then fails the next gap' not in limit


def test_aliasing_radius_is_recorded_per_path_absolutely_and_in_tolerances(tmp_path):
    """Review r221 P2-1: each path's largest per-gap motion, in both units."""
    seed = {'u': 1.5, 'n_e': 1.5e16, 'state': [1.5]}
    main = continue_branch(_LinearStable(), seed, [.25, .5], initial_power=.5,
                           max_step=.125, min_step=.03125, stability_resolution=.0625)
    half = continue_branch(_LinearStable(), seed, [.25, .5], initial_power=.5,
                           max_step=.0625, min_step=.015625, stability_resolution=.03125)

    radius = _with_paths(main, half, u=1.5, n_e=1.5e16)['aliasing_radius']['reactor-000']

    for key, gap in (('continuation_path', .0625), ('continuation_path_half_step', .03125)):
        assert radius[key]['max_u_motion'] == pytest.approx(gap)
        assert radius[key]['max_u_motion_in_u_tolerances'] == pytest.approx(gap / 1.e-4)
        assert radius[key]['max_log_ne_motion'] == pytest.approx(np.log((1.25 + gap) / 1.25))
        assert radius[key]['max_log_ne_motion_in_log_ne_tolerances'] == pytest.approx(
            np.log((1.25 + gap) / 1.25) / 1.e-3)
    assert radius['continuation_path']['gaps'] == 8
    assert radius['continuation_path_half_step']['gaps'] == 16


class _ScanPair(_ParallelPair):
    """The last-gap alias of the r221 case L, reached through scan()."""

    def integrate(self, seed, power=None):
        if not seed.state:
            seed = Seed(seed.u, seed.n_e, (self.a(float(power)),))
        state = super().integrate(seed, power)
        state.update(self_sustained=True, discharge_class='self-sustained')
        return state


def test_scan_without_half_step_path_is_unconfirmed_and_refused_off_reference(tmp_path):
    """Review r221 case L: production scan() checks nothing against a half path.

    The step path aliases onto B at its last gap (see the last-gap test). With
    no half-step path the artifact must not read consistent, and matching
    refuses both branches' states at 0.5 W instead of accepting B.
    """
    adapter = _ScanPair(lambda power: 1. - power ** 3, 2. * 28. / 512.)
    output = tmp_path / 'branches.json'

    payload = scan(adapter, (1., 1.), output, u_steps=1, electron_density_decades=1,
                   powers=(.5,), reference_power=0., continuation_max_step=.5,
                   continuation_min_step=.125, continuation_stability_resolution=.125)

    path = payload['branches'][0]['continuation_path']
    assert adapter.branch(path[-1]['u'], .5) == 'B'
    for u in (path[-1]['u'], adapter.a(.5)):
        with pytest.raises(ValueError, match='no half-step path is recorded'):
            matching_branch(output, {'u': u, 'n_e': 1.e16}, power=.5)
    assert matching_branch(output, {'u': 1., 'n_e': 1.e16}, power=0.) == 'reactor-000'
    consistency = json.loads(output.read_text())['path_consistency']['reactor-000']
    assert consistency['consistent'] is False
    assert consistency['half_step_recorded'] is False
    assert consistency['unconfirmed_by_half_step'] == [0., .125, .25, .375, .5]
    # A branch whose step path certified nothing was cross-checked by nothing.
    assert path_consistency({'id': 'reactor-001', 'continuation_path': []},
                            1.e-4, 1.e-3)['consistent'] is False


def test_scan_refuses_continuation_without_a_declared_resolution(tmp_path):
    """Review r221 P3-6: no silent 64-gap default through scan(powers=...)."""
    with pytest.raises(ValueError, match='declared continuation_stability_resolution'):
        scan(_LinearStable(), (1., 1.), tmp_path / 'branches.json', u_steps=1,
             electron_density_decades=1, powers=(.25,), reference_power=.5)
    assert not (tmp_path / 'branches.json').exists()


@pytest.mark.parametrize('initial_power', [.29, .2999])
def test_limit_4_first_step_of_a_leg_has_no_departure_trigger(initial_power):
    """Review r221 case K: a seed-independent 0.2 jump inside a leg's first step.

    The jump is 2000 u_tol at 0.3 W. The first step has no secant, and the
    corrector ignores its seed, so every gap reverses exactly: the step and the
    20 states after it are accepted. CERTIFICATE_LIMITS[3] declares this.
    """
    adapter = _SilentHop()
    start = adapter.u(initial_power)

    path = continue_branch(adapter, {'u': start, 'n_e': 1.e16, 'state': [start]},
                           [.5], initial_power=initial_power)

    accepted = _accepted_states(path)
    assert path[-1]['power'] == .5 and path[-1]['event_type'] is None
    assert sum(entry['power'] >= .3 for entry in path if entry.get('converged')) == 20
    assert all(entry['u'] == adapter.u(entry['power']) for entry in accepted)
    certificate = path[1]['certificate']
    assert certificate['sub_steps'] == 64
    assert certificate['max_reverse_u_difference'] == 0.
    assert certificate['forward_u_difference'] == 0.
    assert path[1]['corrector_to_tangent_ratio'] is None
    assert ('The first step of each leg has no secant, so there no trigger '
            'remains and a jump inside that step is accepted.') in CERTIFICATE_LIMITS[3]


def test_initial_state_on_another_branch_than_its_seed_is_refused():
    """Review r221 probe M: the seed settles 0.005 (50 u_tol) away, on B."""
    class Pair:
        def integrate(self, seed, power=None):
            u = 1.76 if seed.state[0] >= 1.7549 else 1.75
            return _complete({'u': u, 'n_e': 1.e16, '_state_vector': [u],
                              '_jacobian': [[-1.]]})

    path = continue_branch(Pair(), {'u': 1.755, 'n_e': 1.e16, 'state': [1.755]},
                           [.6], initial_power=.5)

    assert len(path) == 1
    assert path[0]['event_type'] == 'unresolved'
    assert "is not the seed's branch" in path[0]['failure']
    on_branch = continue_branch(Pair(), {'u': 1.75, 'n_e': 1.e16, 'state': [1.75]},
                                [.6], initial_power=.5)
    assert on_branch[-1]['power'] == .6 and on_branch[-1]['converged'] is True


def test_chain_sample_verdicts_and_code_hash_are_recorded(tmp_path):
    """Review r221 P3-3: the artifact shows each sample's verdict and the code."""
    seed = {'u': 1.5, 'n_e': 1.5e16, 'state': [1.5]}
    main = continue_branch(_LinearStable(), seed, [.25], initial_power=.5,
                           max_step=.125, min_step=.03125, stability_resolution=.0625)
    half = continue_branch(_LinearStable(), seed, [.25], initial_power=.5,
                           max_step=.0625, min_step=.015625, stability_resolution=.03125)

    payload = _with_paths(main, half, u=1.5, n_e=1.5e16)

    for entry in main[1:]:
        assert entry['certificate']['sample_verdicts'] == ['stable', 'stable']
        assert [point['stability'] for point in entry['certification_points']] == ['stable']
    source = Path(__file__).parents[4] / 'rmgpy/tools/eedf/branches.py'
    import hashlib
    assert payload['continuation_certificate']['branches_py_sha256'] == hashlib.sha256(
        source.read_bytes()).hexdigest()


def test_last_chain_sample_verdict_is_enforced_although_the_endpoint_replaces_it():
    """The certificate statement says every sample is stable, the last included.

    The last chain sample is solved from the sample before it, the endpoint
    from the step's start, so the two can differ. Here only the last sample,
    seeded from 0.375 W, is unstable at 0.5 W; the endpoint is stable.
    """
    class SeedDependentVerdict:
        def integrate(self, seed, power=None):
            power = float(power)
            unstable = abs(power - .5) < 1.e-12 and seed.state[0] > 1.3
            return _complete({'u': 1. + power, 'n_e': 1.e16,
                              '_state_vector': [1. + power],
                              '_jacobian': [[1. if unstable else -1.]]})

    path = continue_branch(
        SeedDependentVerdict(), {'u': 1.25, 'n_e': 1.e16, 'state': [1.25]}, [.5],
        initial_power=.25, max_step=.25, min_step=.125, stability_resolution=.125)

    event = path[-1]
    assert event['event_type'] == 'stability'
    assert 'of the last chain sample at 0.5 W' in event['failure']
    assert not any(entry['power'] > .25 for entry in _accepted_states(path))


def test_coverage_spans_both_paths_once_and_matching_uses_half_step_states(tmp_path):
    """Review r221 P3-4: legs that meet within roundoff are one covered power."""
    class Linear:
        def integrate(self, seed, power=None):
            power = float(power)
            u = 1.75 - .86 * (power - .5)
            return _complete({'u': u, 'n_e': 1.e16, '_state_vector': [u],
                              '_jacobian': [[-1.]]})

    seed = {'u': 1.75, 'n_e': 1.e16, 'state': [1.75]}
    main = continue_branch(Linear(), seed, REAL_TARGETS, **REAL_SETTINGS)
    half = continue_branch(Linear(), seed, REAL_TARGETS, **HALF_SETTINGS)

    payload = _with_paths(main, half)

    main_powers = [entry['power'] for entry in main if entry.get('converged')]
    main_powers += [point['power'] for entry in main
                    for point in entry.get('certification_points', [])]
    half_only = next(point['power'] for entry in half
                     for point in entry.get('certification_points', [])
                     if min(abs(point['power'] - other) for other in main_powers) > 1.e-6)
    artifact = tmp_path / 'branches.json'
    artifact.write_text(json.dumps(payload))
    u = 1.75 - .86 * (half_only - .5)
    assert matching_branch(artifact, {'u': u, 'n_e': 1.e16}, power=half_only) == 'reactor-000'
    coverage = payload['certified_coverage']['reactor-000']
    powers = coverage['powers']
    assert half_only in powers
    assert all(right - left > 1.e-9 for left, right in zip(powers, powers[1:]))
    assert len(powers) == 241
    assert coverage['paths'] == ['continuation_path', 'continuation_path_half_step']


def test_record_continuation_refuses_paths_certified_at_other_tolerances():
    """Review r221 P3-5: a path judged at u_tol 0.1 says nothing at 1e-4."""
    seed = {'u': 1., 'n_e': 1.e16, 'state': [1.]}
    loose = continue_branch(_LinearStable(), seed, [.5], initial_power=0.,
                            max_step=.25, min_step=.0625, stability_resolution=.125,
                            u_tolerance=.1)
    payload = {'reference_power': 0., 'u_tolerance': 1.e-4, 'log_ne_tolerance': 1.e-3,
               'branches': [{'id': 'reactor-000', 'u': 1., 'n_e': 1.e16,
                             'continuation_path': loose}]}

    with pytest.raises(ValueError, match="continuation_path: settings u_tolerance 0.1 "
                                         "differs from the artifact's 0.0001"):
        record_continuation(payload)
    del loose[0]['continuation_settings']
    with pytest.raises(ValueError, match='certificate at 0.25 W u_tolerance 0.1'):
        record_continuation(payload)


class _BudgetStopCubic(_TimeCubic):
    """Critical slowing down: simulations short of the fold exhaust their budget."""

    def integrate(self, seed, power=None):
        if self.slow_window[0] < float(power) < self.slow_window[1]:
            raise RuntimeError('reactor returned without reaching steady state')
        return super().integrate(seed, power)


def test_limit_6_time_budget_stop_on_the_real_adapter_narrows_a_fold_bracket():
    """Review r221 P3-7: the real adapter cannot tell a budget stop from a fold.

    PlasmaReactorAdapter raises one failure whenever simulate returns without
    steady state, which includes exhausting the deck's time budget, and its
    is_fold_failure accepts it. Before the fold at 0.3849 W such stops narrow
    the bracket below the fold and label it fold.
    """
    class BudgetReactor:
        steady_state_reached = False
        initial_mole_fractions = {}
        electron_kinetics = {'initial_reduced_field': (1., 'Td')}

        def simulate(self, *args, **kwargs):
            pass

    adapter = PlasmaReactorAdapter(
        lambda seed, power: ReactorRun(BudgetReactor(), ['e-'], [], object(), object()))
    with pytest.raises(RuntimeError) as stopped:
        adapter.integrate(Seed(1., 1.e16, (1., 0.)), power=.5)
    failure = '{}: {}'.format(type(stopped.value).__name__, stopped.value)
    assert PlasmaReactorAdapter.is_fold_failure(failure)
    assert not PlasmaReactorAdapter.is_fold_failure('TimeoutError: worker exceeded 1 s')

    fold = 2. / np.sqrt(27.)
    cubic = _BudgetStopCubic(slow_window=(.379, fold))
    cubic.is_fold_failure = PlasmaReactorAdapter.is_fold_failure
    start = cubic.offset + cubic.scale * min(
        root.real for root in np.roots([1., 0., -1., -.35]) if abs(root.imag) < 1.e-10)
    path = continue_branch(
        cubic, {'u': start, 'n_e': 1.e16, 'state': [start]}, [.5],
        initial_power=.35, max_step=.025, min_step=1.e-4, stability_resolution=.025)

    event = path[-1]
    assert event['event_type'] == 'fold'
    assert event['fold_bracket'][1] < fold
    assert 'a time-budget stop narrows a fold bracket and can label it fold' in CERTIFICATE_LIMITS[5]


def test_limit_2_states_differing_only_off_identity_are_one_branch(tmp_path):
    """Review r221 P3-8: identity is (u, ln n_e); other components do not split.

    Two seeds settle on states with equal u and n_e but a hidden component of
    0 and 1. The scan clusters them as one branch with two members.
    """
    class HiddenComponent:
        def integrate(self, seed, power=None):
            hidden = 0. if seed.n_e < 1.e16 else 1.
            return _complete({'u': 1., 'n_e': 1.e16, '_state_vector': [1., hidden],
                              '_jacobian': [[-1., 0.], [0., -1.]]})

    payload = scan(HiddenComponent(), (1., 1.), tmp_path / 'branches.json', u_steps=1,
                   electron_density_decades=2, reference_power=.5)

    assert len(payload['branches']) == 1
    assert payload['branches'][0]['members'] == 2
    assert [attempt['seed']['n_e'] for attempt in payload['attempts']] == [1.e14, 1.e18]
    assert [HiddenComponent().integrate(Seed(1., n_e))['_state_vector'][1]
            for n_e in (1.e14, 1.e18)] == [0., 1.]
    assert "components are one branch by the scan's definition" in CERTIFICATE_LIMITS[1]


@pytest.mark.parametrize('width, caught', [(.004, True), (.001, False), (.0004, False)])
def test_limit_5_unstable_window_narrower_than_the_resolution_can_be_missed(width, caught):
    """Review r221 probe H at the real settings: resolution 0.00625 W.

    A Gaussian excursion pushes the slowest mode (-483) past zero over an
    interval of 0.93 * width. At 0.0037 W it is caught; at 0.00093 and
    0.00037 W, centred between samples, the path runs through it.
    """
    centre = .3021

    class Excursion:
        def integrate(self, seed, power=None):
            power = float(power)
            mode = -483. + 600. * np.exp(-((power - centre) / width) ** 2)
            u = 1.75 - .86 * (power - .5)
            return _complete({'u': u, 'n_e': 1.e16, '_state_vector': [u],
                              '_jacobian': [[mode]]})

    path = continue_branch(Excursion(), {'u': 1.75, 'n_e': 1.e16, 'state': [1.75]},
                           REAL_TARGETS, **REAL_SETTINGS)

    unstable = 2. * width * np.sqrt(np.log(600. / 483.))
    if caught:
        assert unstable < REAL_SETTINGS['stability_resolution']
        assert path[-1]['event_type'] == 'stability'
        low, high = path[-1]['event_bracket']
        assert low < centre + unstable / 2. and high > centre - unstable / 2.
    else:
        assert path[-1]['power'] == 1. and path[-1]['event_type'] is None
        assert any(entry['power'] < centre for entry in _accepted_states(path))
    assert 'Unstable windows narrower than the stability resolution can be missed' in \
        CERTIFICATE_LIMITS[4]
