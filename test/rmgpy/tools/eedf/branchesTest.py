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

import copy
import ctypes
import json
import logging
import multiprocessing
import os
import re
import select
import signal
import struct
import subprocess
import sys
import threading
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from rmgpy.tools.eedf.validation import split_branches
from rmgpy.tools.eedf.branches import (
    CERTIFICATE_LIMITS,
    CERTIFICATE_STATEMENT,
    CONTINUATION_CERTIFICATE,
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


REPOSITORY = Path(__file__).resolve().parents[4]


def _repository_env():
    """Environment for a Python subprocess that must import this checkout's rmgpy.

    The interpreter's own path may link another RMG-Py checkout; this
    repository goes first on PYTHONPATH whatever the caller exported.
    """
    env = dict(os.environ)
    inherited = env.get('PYTHONPATH')
    env['PYTHONPATH'] = os.pathsep.join([str(REPOSITORY)] + ([inherited] if inherited else []))
    return env


def _native_probe(mode):
    _run_native_script('_native_protocol_probe.py', mode, 'NATIVE_PROTOCOL')


def _native_interrupt(mode):
    _run_native_script('_native_interrupt_probe.py', mode, 'NATIVE_INTERRUPT')


def _run_native_script(script, mode, marker):
    environment = dict(os.environ, PATH='')
    environment.pop('PYTHONPATH', None)
    result = subprocess.run([sys.executable, str(Path(__file__).with_name(script)), mode],
                            env=environment, capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, result.stdout + result.stderr
    assert marker in result.stdout, result.stdout
    print(result.stdout)


class _RawSigaction(ctypes.Structure):
    # glibc's struct sigaction on Linux.
    _fields_ = [('handler', ctypes.c_void_p), ('mask', ctypes.c_ulong * 16),
                ('flags', ctypes.c_int), ('restorer', ctypes.c_void_p)]


@pytest.fixture(autouse=True)
def _python_owns_sigint():
    """Give each test the SIGINT handler it assumes: Python's own.

    Earlier code in the same pytest process can leave a C handler behind
    signal.getsignal()'s back: RDKit's embedding (run by the excited-export
    census) leaves one that swallows the next SIGINT and resets the signal
    to SIG_DFL. The exact prior disposition is put back afterwards.
    """
    linux = sys.platform.startswith('linux')
    saved = _RawSigaction()
    if linux:
        ctypes.CDLL(None).sigaction(signal.SIGINT, None, ctypes.byref(saved))
    previous = signal.signal(signal.SIGINT, signal.default_int_handler)
    yield
    signal.signal(signal.SIGINT, previous)
    if linux:
        ctypes.CDLL(None).sigaction(signal.SIGINT, ctypes.byref(saved), None)


STALE_V1_ARTIFACT = (Path(__file__).parents[2] / 'test_data' / 'operating_branches'
                     / 'reversible_chain_v1_marginal_final_sample.json')
class _Branches:
    """Stable states u = c + shape(P), n_e = 1e16: one branch per offset c.

    The packed state carries c, so a continuation step follows its own branch
    however close another lies; a scan seed takes the offset nearest
    u - shape(P).
    """

    def __init__(self, *offsets, shape=lambda power: power):
        self.offsets = offsets
        self.shape = shape

    def integrate(self, seed, power=None):
        power = float(power)
        guess = seed.state[1] if len(seed.state) == 2 else seed.u - self.shape(power)
        offset = min(self.offsets, key=lambda value: abs(value - guess))
        u = offset + self.shape(power)
        return {'u': u, 'n_e': 1.e16, '_state_vector': [u, offset],
                '_jacobian': [[-1., 0.], [0., -1.]],
                '_half_step_jacobian': [[-1., 0.], [0., -1.]],
                '_conservation_constraints': [], 'self_sustained': True,
                'discharge_class': 'self-sustained'}


def _scanned(tmp_path, adapter, u_limits, u_steps, powers, reference_power=.5,
             main=(.125, .03125, .0625), half=(.0625, .015625, .03125),
             u_tolerance=1.e-6, log_ne_tolerance=1.e-6, half_step=True):
    """Write ``branches.json`` as the real-deck driver does, and return it and its payload.

    scan() solves, clusters and judges the seeds and continues each stable
    branch on the step path; the half-step path is continued from the same
    state and the artifact recorded again. Malformed cases edit this output.
    """
    output = tmp_path / 'branches.json'
    payload = scan(adapter, u_limits, output, u_steps=u_steps, electron_density_decades=1,
                   electron_density_limits=(1.e16, 1.e16), reference_power=reference_power,
                   powers=powers, continuation_max_step=main[0], continuation_min_step=main[1],
                   continuation_stability_resolution=main[2], u_tolerance=u_tolerance,
                   log_ne_tolerance=log_ne_tolerance)
    if half_step:
        for branch in payload['branches']:
            start = branch['continuation_path'][0]
            branch['continuation_path_half_step'] = continue_branch(
                adapter, {key: start[key] for key in ('u', 'n_e', 'state')}, powers,
                initial_power=reference_power, max_step=half[0], min_step=half[1],
                stability_resolution=half[2], u_tolerance=u_tolerance,
                log_ne_tolerance=log_ne_tolerance)
        output.write_text(json.dumps(record_continuation(payload)))
    return output, json.loads(output.read_text())


def _write(path, payload):
    path.write_text(json.dumps(payload))
    return path


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
    try:
        branches.scan(_PostResultAdapter(pid_path), (0., 0.), output, u_steps=1,
                      electron_density_decades=1, processes=1, timeout=.05, reference_power=.5)
    except branches.NativeWorkerStatusError:
        assert not Path(output).exists()
        queue.put((time.monotonic() - started, 'refused'))
    else:
        raise AssertionError('pipe result accepted before normal completion')


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
    with pytest.raises(ValueError, match="'reactor-000' cannot be matched: no half-step path"):
        declared_branch(tmp_path / 'branches.json', stable['id'])
    with pytest.raises(ValueError, match="'reactor-001' cannot be matched: its reference "
                                         "verdict is 'unstable', so it is recorded for "
                                         "reporting only"):
        declared_branch(tmp_path / 'branches.json', 'reactor-001')


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
    assert attempts == 'refused'  # Native N1-b: send is not normal completion.


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
    # scan returns only after reaping every worker it started, so the
    # sibling is no longer a child of this process, dead or alive.
    with pytest.raises(ChildProcessError):
        os.waitpid(int(sibling_pid.read_text()), os.WNOHANG)


class _SlowToDie:
    """A worker the kernel takes longer than every grace period to reap, as under heavy load."""

    def __init__(self, ignores_sigterm=True, dies_after=.5):
        self.ignores_sigterm = ignores_sigterm
        self.dies_after = dies_after
        self.death = None
        self.reaped = False
        # A real running child stands in for the OS's view of the worker.
        self.stand_in = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)'],
                                         env=_repository_env())
        self.pid = self.stand_in.pid

    def release(self):
        self.stand_in.kill()
        self.stand_in.wait()

    def is_alive(self):
        return not self.reaped

    def terminate(self):
        if not self.ignores_sigterm and self.death is None:
            self.death = time.monotonic() + self.dies_after

    def kill(self):
        if self.death is None:
            self.death = time.monotonic() + self.dies_after

    def join(self, timeout=None):
        if self.death is None:
            time.sleep(timeout)
            return
        wait = self.death - time.monotonic()
        if timeout is not None and wait > timeout:
            time.sleep(timeout)
            return
        time.sleep(max(wait, 0.))
        self.reaped = True


@pytest.mark.parametrize('ignores_sigterm', [True, False])
def test_a_killed_worker_is_reaped_before_shutdown_returns(ignores_sigterm):
    import rmgpy.tools.eedf.branches as branches

    process = _SlowToDie(ignores_sigterm)
    try:
        branches._reap_worker(process, terminate=True)
    finally:
        process.release()
    assert process.reaped


def test_a_worker_that_outlives_sigkill_is_logged_not_waited_on(monkeypatch, caplog):
    import rmgpy.tools.eedf.branches as branches

    monkeypatch.setattr(branches, 'WORKER_KILL_TIMEOUT', .05)
    process = _SlowToDie(True, dies_after=3600.)
    started = time.monotonic()
    try:
        with pytest.raises(branches.NativeWorkerStatusError, match=str(process.pid)):
            branches._reap_worker(process, terminate=True)
    finally:
        process.release()
    assert time.monotonic() - started < 2.
    assert not process.reaped
    # Native bounded refusal replaces the old warning-only outcome. The
    # retained exec regression verifies scan-level pending-PID logging.


class _SleepsIgnoringSigterm:
    def integrate(self, seed, power=None):
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
        time.sleep(60.)


class _ReportsSignalMask:
    def integrate(self, seed, power=None):
        return sorted(int(s) for s in signal.pthread_sigmask(signal.SIG_BLOCK, set()))


class _InterruptingContext:
    """Delivers a real SIGINT to this process at a chosen boundary.

    ``where`` is 'after-start' (just after a worker has started, before the
    caller can register it) or 'first-join' (the first join, which only the
    cleanup in ``finally`` reaches because poll() fails on the first pass).
    """

    def __init__(self, where):
        self.where = where
        self.pids = []
        self.fired = False
        from rmgpy.tools.eedf import _native_context
        real = _native_context.NativeForkContext()
        outer = self

        class Process(real.Process):
            def start(self):
                super().start()
                outer.pids.append(self.pid)
                if outer.where == 'after-start' and not outer.fired:
                    outer.fired = True
                    os.kill(os.getpid(), signal.SIGINT)

            def join(self, timeout=None):
                if outer.where == 'first-join' and not outer.fired:
                    outer.fired = True
                    os.kill(os.getpid(), signal.SIGINT)
                super().join(timeout)

        class Connection:
            def __init__(self, connection):
                self._connection = connection

            def poll(self, *args):
                if outer.where == 'first-join':
                    raise OSError('parent poll failure')
                return self._connection.poll(*args)

            def __getattr__(self, name):
                return getattr(self._connection, name)

        self._real = real
        self.Process = lambda **kwargs: Process(**kwargs)
        self._connection = Connection

    def Pipe(self, **kwargs):
        parent, child = self._real.Pipe(**kwargs)
        return self._connection(parent), child


@pytest.mark.skipif(not hasattr(signal, 'pthread_sigmask'), reason='needs pthread_sigmask')
@pytest.mark.parametrize('where, seeds', [('after-start', [0]), ('first-join', [0, 1])])
def test_sigint_inside_worker_bookkeeping_still_reaps_every_started_worker(
        monkeypatch, where, seeds):
    import rmgpy.tools.eedf.branches as branches

    context = _InterruptingContext(where)
    native_context = branches._native_context
    monkeypatch.setattr(native_context, 'NativeForkContext', lambda: context)
    monkeypatch.setattr(native_context, 'NativeForkProcess', context.Process)
    original_cleanup = native_context.Owner.cleanup
    def cleanup(owner, seconds=5.):
        if where == 'first-join' and not context.fired:
            context.fired = True
            os.kill(os.getpid(), signal.SIGINT)
        return original_cleanup(owner, seconds)
    monkeypatch.setattr(native_context.Owner, 'cleanup', cleanup)
    try:
        with pytest.raises(KeyboardInterrupt):
            branches._integrate_many(_SleepsIgnoringSigterm(), seeds, len(seeds), 30.)
        assert context.fired
        assert len(context.pids) == len(seeds)
        for pid in context.pids:
            # waitpid-level: reaped means the kernel no longer has the child.
            with pytest.raises(ChildProcessError):
                os.waitpid(pid, os.WNOHANG)
    finally:
        for pid in context.pids:
            try:
                os.kill(pid, signal.SIGKILL)
                os.waitpid(pid, 0)
            except (ProcessLookupError, ChildProcessError):
                pass
    blocked = signal.pthread_sigmask(signal.SIG_BLOCK, set())
    assert signal.SIGINT not in blocked and signal.SIGTERM not in blocked


class _ReportsHandlers:
    def integrate(self, seed, power=None):
        return [repr(signal.getsignal(signal.SIGINT)), repr(signal.getsignal(signal.SIGTERM)),
                sorted(int(s) for s in signal.pthread_sigmask(signal.SIG_BLOCK, set()))]


def _raise_system_exit(signum, frame):
    raise SystemExit(128 + signum)


@pytest.mark.skipif(not hasattr(signal, 'pthread_sigmask'), reason='needs pthread_sigmask')
def test_workers_run_with_the_callers_handlers_and_mask_and_the_scan_restores_them():
    """The scan's deferral handlers never reach a worker or outlive the scan."""
    import rmgpy.tools.eedf.branches as branches

    previous = signal.signal(signal.SIGTERM, _raise_system_exit)
    try:
        before = signal.pthread_sigmask(signal.SIG_BLOCK, set())
        results = branches._integrate_many(_ReportsHandlers(), [0, 1], 2, 30.)
        assert signal.pthread_sigmask(signal.SIG_BLOCK, set()) == before
        assert signal.getsignal(signal.SIGINT) is signal.default_int_handler
        assert signal.getsignal(signal.SIGTERM) is _raise_system_exit
    finally:
        signal.signal(signal.SIGTERM, previous)
    for status, (sigint, sigterm, mask) in results:
        assert status == 'ok'
        assert sigint == repr(signal.default_int_handler)
        assert sigterm == repr(_raise_system_exit)
        assert mask == sorted(int(s) for s in before)


def test_a_worker_reaped_elsewhere_is_not_reported_alive(caplog):
    """A caller's SIGCHLD reaper leaves multiprocessing believing the worker runs."""
    import rmgpy.tools.eedf.branches as branches

    owner = branches._native_context.Owner(60.)
    process = owner.process(lambda: None, ())
    process.start()
    os.waitpid(process.pid, 0)
    started = time.monotonic()
    branches._reap_worker(process, terminate=True)
    assert time.monotonic() - started < 1.
    assert 'still alive' not in caplog.text
    with pytest.raises(branches.NativeWorkerStatusError):
        branches._native_context.accepted(process, ('ok', None))
    assert owner.cleanup()


# Cancellation probes. Each trial runs in a controller that is a Linux child
# subreaper: the scan's owner is a fork of the controller, and a worker the
# owner did not reap before it exited is adopted by the controller, alive or
# dead, so "reaped by its owner" is observed, not inferred.

_needs_subreaper = pytest.mark.skipif(
    not sys.platform.startswith('linux'), reason='needs a Linux child subreaper')
_ROLES = ('worker', 'inner')


def _emit(fd, **event):
    os.write(fd, (json.dumps(event) + '\n').encode())


class _Sleeper:
    """Reports its pid, ignores SIGTERM and sleeps: only SIGKILL or its lifeline ends it."""

    def __init__(self, fd, role='worker', ready=None):
        self.fd, self.role, self.ready = fd, role, ready

    def integrate(self, seed, power=None):
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
        _emit(self.fd, **{self.role: os.getpid()})
        if self.ready is not None:
            os.write(self.ready, b'r')
        time.sleep(30.)


class _NestedScan:
    """Each worker owns a scan of its own, with one sleeping inner worker."""

    def __init__(self, fd, output):
        self.fd, self.output = fd, Path(output)

    def integrate(self, seed, power=None):
        _emit(self.fd, worker=os.getpid())
        return scan(_Sleeper(self.fd, 'inner'), (0., 0.),
                    self.output / 'inner-{}.json'.format(os.getpid()), u_steps=1,
                    electron_density_decades=1, processes=1, timeout=30.,
                    reference_power=.5, recheck_timeouts=False)


def _scan_target(fd, output, adapter=None, *, seeds=2, sigterm=signal.SIG_DFL,
                 idle_thread=False, main_blocks_sigint=False, in_thread=False,
                 fork_helper=False):
    """The owner: one public scan() under the given signal arrangement."""
    signal.signal(signal.SIGINT, signal.default_int_handler)
    signal.signal(signal.SIGTERM, sigterm)
    if idle_thread:
        threading.Thread(target=time.sleep, args=(30.,), daemon=True).start()
    if main_blocks_sigint:
        signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGINT})
    if fork_helper:
        def helper():
            # An unrelated fork made while the workers run inherits every
            # descriptor the owner holds at that moment.
            time.sleep(.3)
            # 'native': a C library forking directly; Python's fork hooks do not run.
            pid = ctypes.CDLL(None).fork() if fork_helper == 'native' else os.fork()
            if pid == 0:
                time.sleep(30.)
                os._exit(0)
            _emit(fd, helper=pid)
        threading.Thread(target=helper, daemon=True).start()

    def run():
        try:
            scan(adapter or _Sleeper(fd), (0., 1.), output, u_steps=seeds,
                 electron_density_decades=1, processes=seeds, timeout=30.,
                 reference_power=.5, recheck_timeouts=False)
            _emit(fd, outcome='return')
        except BaseException as exc:
            _emit(fd, outcome=type(exc).__name__)
    if in_thread:
        thread = threading.Thread(target=run)
        thread.start()
        thread.join()
    else:
        run()


def _fate(pid):
    """Classify a worker once the owner has exited; reaps an adopted corpse."""
    try:
        waited, _ = os.waitpid(pid, os.WNOHANG)
    except ChildProcessError:
        return 'reaped-by-owner'
    return 'alive' if waited == 0 else 'outlived-owner'


def _gone_within(pid, timeout):
    """Seconds until an adopted pid exits, or None when it outlives timeout."""
    started = time.monotonic()
    while time.monotonic() - started < timeout:
        try:
            waited, _ = os.waitpid(pid, os.WNOHANG)
        except ChildProcessError:
            return time.monotonic() - started
        if waited:
            return time.monotonic() - started
        time.sleep(.005)
    return None


def _owner_trial(target, on_event, deadline=20.):
    """Fork target(fd) as the owner; feed its events to on_event until it exits."""
    read_end, write_end = os.pipe()
    owner = os.fork()
    if owner == 0:
        os.close(read_end)
        try:
            target(write_end)
        finally:
            os._exit(0)
    os.close(write_end)
    os.set_blocking(read_end, False)
    events, buffer, status = [], b'', None

    def drain():
        nonlocal buffer
        while True:
            try:
                chunk = os.read(read_end, 65536)
            except BlockingIOError:
                break
            if not chunk:
                break
            buffer += chunk
        while b'\n' in buffer:
            line, buffer = buffer.split(b'\n', 1)
            events.append(json.loads(line))
            on_event(events[-1], owner, events)

    end = time.monotonic() + deadline
    while status is None and time.monotonic() < end:
        select.select([read_end], [], [], .005)
        drain()
        waited, raw = os.waitpid(owner, os.WNOHANG)
        if waited:
            status = os.waitstatus_to_exitcode(raw)
    if status is None:
        os.kill(owner, signal.SIGKILL)
        os.waitpid(owner, 0)
        status = 'timeout'
    drain()
    os.close(read_end)
    pids = {}
    for event in events:
        for role in _ROLES:
            if role in event:
                pids[event[role]] = role
    return status, events, pids


def _as_subreaper(trial):
    """Run trial() in a forked child subreaper; return its JSON result."""
    read_end, write_end = os.pipe()
    controller = os.fork()
    if controller == 0:
        os.close(read_end)
        code = 0
        try:
            assert ctypes.CDLL(None).prctl(36, 1, 0, 0, 0) == 0  # PR_SET_CHILD_SUBREAPER
            result = {'result': trial()}
        except BaseException as exc:
            result = {'error': '{}: {}'.format(type(exc).__name__, exc)}
            code = 1
        try:
            os.write(write_end, (json.dumps(result) + '\n').encode())
        finally:
            # Remove everything adopted, then exit: no trial outlives its test.
            for _ in range(50):
                children = Path('/proc/self/task/{}/children'.format(os.getpid())).read_text().split()
                if not children:
                    break
                for child in children:
                    try:
                        os.kill(int(child), signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                while True:
                    try:
                        if os.waitpid(-1, os.WNOHANG)[0] == 0:
                            break
                    except ChildProcessError:
                        break
                time.sleep(.01)
            os._exit(code)
    os.close(write_end)
    buffer = b''
    end = time.monotonic() + 60.
    while b'\n' not in buffer and time.monotonic() < end:
        if select.select([read_end], [], [], .1)[0]:
            chunk = os.read(read_end, 65536)
            if not chunk:
                break
            buffer += chunk
    os.close(read_end)
    os.waitpid(controller, 0)
    result = json.loads(buffer.split(b'\n', 1)[0])
    assert 'error' not in result, result['error']
    return result['result']


def _on_workers(count, *signals, gap=.02, roles=_ROLES):
    """Send signals to the owner once count workers have reported."""
    def on_event(event, owner, events):
        if any(role in event for role in roles):
            seen = sum(any(role in e for role in roles) for e in events)
            if seen == count:
                for number, signum in enumerate(signals):
                    if number:
                        time.sleep(gap)
                    os.kill(owner, signum)
    return on_event


def _cancel_trial(tmp_path, signals, count=2, **arrangement):
    def trial():
        status, events, pids = _owner_trial(
            lambda fd: _scan_target(fd, tmp_path / 'branches.json', **arrangement),
            _on_workers(count, *signals))
        outcome = next((e['outcome'] for e in events if 'outcome' in e), None)
        return {'status': status, 'outcome': outcome,
                'fates': {str(pid): _fate(pid) for pid in pids}}
    return _as_subreaper(trial)


@_needs_subreaper
@pytest.mark.parametrize('arrangement', [
    {'idle_thread': True},
    {'idle_thread': True, 'main_blocks_sigint': True},
])
@pytest.mark.parametrize('trial', range(2))
def test_sigint_delivered_through_another_thread_still_reaps_every_worker(
        tmp_path, arrangement, trial):
    """r255 P2-1: two SIGINTs, the second inside cleanup, with an unmasked thread."""
    result = _cancel_trial(tmp_path, (signal.SIGINT, signal.SIGINT), **arrangement)
    assert result['outcome'] == 'KeyboardInterrupt', result
    assert result['status'] == 0, result
    assert len(result['fates']) == 2 and set(result['fates'].values()) == {'reaped-by-owner'}, result


@_needs_subreaper
@pytest.mark.parametrize('signals', [
    (signal.SIGTERM,),
    (signal.SIGINT, signal.SIGTERM),
], ids=['sigterm', 'sigint-then-sigterm'])
def test_default_sigterm_reaps_every_worker_then_terminates_with_its_status(tmp_path, signals):
    """r255 P2-2: SIGTERM keeps its -15 status, after the scan reaped its workers."""
    if len(signals) == 1:
        result = _cancel_trial(tmp_path, signals)
    else:
        # Native cleanup can finish before the inherited 20 ms sender gap.
        # Deliver TERM while cleanup is actually pending, not to a dead owner.
        from rmgpy.tools.eedf import _native_context as native_context
        from rmgpy.tools.eedf import branches

        def trial():
            def target(fd):
                cleanup = native_context.Owner.cleanup

                def at_cleanup(owner, seconds=5.):
                    _emit(fd, cleanup_ready=True)
                    until = time.monotonic() + 3.
                    while signal.SIGTERM not in branches._DEFERRAL.pending:
                        assert time.monotonic() < until, 'TERM was not delivered during cleanup'
                        time.sleep(.001)
                    return cleanup(owner, seconds)

                native_context.Owner.cleanup = at_cleanup
                _scan_target(fd, tmp_path / 'branches.json')

            on_workers = _on_workers(2, signal.SIGINT)

            def on_event(event, owner, events):
                on_workers(event, owner, events)
                if event.get('cleanup_ready'):
                    os.kill(owner, signal.SIGTERM)

            status, events, pids = _owner_trial(target, on_event)
            return {'status': status, 'cleanup_ready': any(e.get('cleanup_ready') for e in events),
                    'fates': {str(pid): _fate(pid) for pid in pids}}

        result = _as_subreaper(trial)
        assert result['cleanup_ready'], result
    assert result['status'] == -signal.SIGTERM, result
    assert len(result['fates']) == 2 and set(result['fates'].values()) == {'reaped-by-owner'}, result


@_needs_subreaper
def test_a_python_sigterm_handler_is_chained_after_every_worker_is_reaped(tmp_path):
    result = _cancel_trial(tmp_path, (signal.SIGTERM,), sigterm=_raise_system_exit)
    assert result['outcome'] == 'SystemExit', result
    assert set(result['fates'].values()) == {'reaped-by-owner'}, result


def _sigterm_inside_start(fd, output):
    """Default SIGTERM raised by the owner itself right after a fork, before registration."""
    from rmgpy.tools.eedf import _native_context
    real = _native_context.NativeForkContext()

    class Process(real.Process):
        def start(self):
            super().start()
            _emit(fd, worker=self.pid)
            os.kill(os.getpid(), signal.SIGTERM)

    class Context:
        Pipe = staticmethod(real.Pipe)

        @staticmethod
        def Process(**kwargs):
            return Process(**kwargs)

    _native_context.NativeForkProcess = Process
    _scan_target(fd, output, seeds=1)


@_needs_subreaper
def test_default_sigterm_during_worker_start_reaps_the_worker(tmp_path):
    """r255 P2-2 startup: the started worker is reaped before the owner dies of SIGTERM."""
    def trial():
        status, events, pids = _owner_trial(
            lambda fd: _sigterm_inside_start(fd, tmp_path / 'branches.json'),
            lambda event, owner, events: None)
        return {'status': status, 'fates': {str(pid): _fate(pid) for pid in pids}}
    result = _as_subreaper(trial)
    assert result['status'] == -signal.SIGTERM, result
    assert result['fates'] and set(result['fates'].values()) == {'reaped-by-owner'}, result


def _sigint_at_nth_line_after_poll_failure(fd, output, nth):
    """An ordinary poll failure, then a real SIGINT at the nth traced line after it.

    Line by line this walks the exception transfer into cleanup, the stretch
    r255 P2-3 measured as unprotected before the cleanup mask.
    """
    import rmgpy.tools.eedf.branches as branches

    signal.signal(signal.SIGINT, signal.default_int_handler)
    real = multiprocessing.get_context('fork')
    ready_r, ready_w = os.pipe()
    owner = os.getpid()
    lines = [None]

    class Connection:
        def __init__(self, connection):
            self.connection = connection

        def poll(self, *args):
            needed = 2
            while needed:
                needed -= len(os.read(ready_r, needed))
            lines[0] = 0
            raise OSError('parent poll failure')

        def __getattr__(self, name):
            return getattr(self.connection, name)

    class Context:
        Process = staticmethod(real.Process)

        @staticmethod
        def Pipe(**kwargs):
            parent, child = real.Pipe(**kwargs)
            return Connection(parent), child

    branches._native_context.NativeForkContext = Context

    def trace(frame, event, arg):
        if os.getpid() != owner:
            return None
        if event == 'line' and lines[0] is not None:
            if lines[0] == nth:
                lines[0] = None
                os.kill(owner, signal.SIGINT)
            else:
                lines[0] += 1
        return trace

    sys.settrace(trace)
    try:
        branches._integrate_many(_Sleeper(fd, ready=ready_w), [0, 1], 2, 30.)
        _emit(fd, outcome='return')
    except BaseException as exc:
        sys.settrace(None)
        _emit(fd, outcome=type(exc).__name__, fired=lines[0] is None)


@_needs_subreaper
@pytest.mark.parametrize('nth', [0, 1, 2, 3, 4, 6, 9, 14])
def test_sigint_entering_cleanup_after_an_exception_still_reaps_every_worker(tmp_path, nth):
    def trial():
        status, events, pids = _owner_trial(
            lambda fd: _sigint_at_nth_line_after_poll_failure(fd, tmp_path / 'b.json', nth),
            lambda event, owner, events: None)
        outcome = next((e for e in events if 'outcome' in e), {})
        return {'status': status, 'outcome': outcome,
                'fates': {str(pid): _fate(pid) for pid in pids}}
    result = _as_subreaper(trial)
    assert result['outcome'] == {'outcome': 'KeyboardInterrupt', 'fired': True}, result
    assert len(result['fates']) == 2 and set(result['fates'].values()) == {'reaped-by-owner'}, result


@_needs_subreaper
def test_cancelling_a_nested_scan_leaves_no_worker_at_either_level(tmp_path):
    """r255 nested cancel: the inner worker never outlives its own owner."""
    def trial():
        status, events, pids = _owner_trial(
            lambda fd: _scan_target(fd, tmp_path / 'outer.json',
                                    _NestedScan(fd, tmp_path), seeds=1),
            _on_workers(1, signal.SIGINT, roles=('inner',)))
        outcome = next((e['outcome'] for e in events if 'outcome' in e), None)
        return {'status': status, 'outcome': outcome,
                'gone_after': {pids[pid]: _gone_within(pid, 5.) for pid in pids}}
    result = _as_subreaper(trial)
    assert result['outcome'] == 'KeyboardInterrupt', result
    assert set(result['gone_after']) == {'worker', 'inner'}, result
    assert all(seconds is not None for seconds in result['gone_after'].values()), result


@_needs_subreaper
@pytest.mark.parametrize('arrangement', [
    {},
    {'in_thread': True},
    {'fork_helper': True},
], ids=['main-thread', 'scan-thread', 'unrelated-fork'])
def test_workers_exit_when_their_owner_is_killed(tmp_path, arrangement):
    """Layer 1: SIGKILL of the owner mid-scan leaves no worker at any level.

    One worker owns a nested scan, so the lifeline must not be held open by
    a sibling, a nested worker or an unrelated fork of the owner.
    """
    expected = 3 + bool(arrangement.get('fork_helper'))

    def on_event(event, owner, events):
        reported = [e for e in events if e.keys() & {'worker', 'inner', 'helper'}]
        if event in reported and len(reported) == expected:
            os.kill(owner, signal.SIGKILL)

    def trial():
        status, events, pids = _owner_trial(
            lambda fd: _scan_target(fd, tmp_path / 'b.json',
                                    _MixedNesting(fd, tmp_path), **arrangement),
            on_event)
        return {'status': status,
                'gone_after': {str(pid): _gone_within(pid, 5.) for pid in pids}}
    result = _as_subreaper(trial)
    assert result['status'] == -signal.SIGKILL, result
    assert len(result['gone_after']) == 3, result
    assert all(seconds is not None for seconds in result['gone_after'].values()), result


class _MixedNesting:
    """Seed 0 owns a nested scan; seed 1 just sleeps."""

    def __init__(self, fd, output):
        self.nested = _NestedScan(fd, output)
        self.sleeper = _Sleeper(fd)

    def integrate(self, seed, power=None):
        if seed.u == 0.:
            return self.nested.integrate(seed, power)
        return self.sleeper.integrate(seed, power)


class _HoldsTheGil:
    """A native call that keeps the GIL: the worker runs no Python meanwhile."""

    def __init__(self, fd):
        self.fd = fd

    def integrate(self, seed, power=None):
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
        _emit(self.fd, worker=os.getpid())
        ctypes.PyDLL(None).sleep(12)


class _ForksNatively:
    """An adapter whose native code forks a helper of its own (not a scan process)."""

    def __init__(self, fd):
        self.fd = fd

    def integrate(self, seed, power=None):
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
        _emit(self.fd, worker=os.getpid())
        pid = ctypes.CDLL(None).fork()
        if pid == 0:
            time.sleep(30.)
            os._exit(0)
        _emit(self.fd, helper=pid)
        time.sleep(30.)


def _kill_owner_when(expected):
    def on_event(event, owner, events):
        reported = [e for e in events if e.keys() & {'worker', 'inner', 'helper'}]
        if event in reported and len(reported) == expected:
            time.sleep(.05)  # let a native call begin
            os.kill(owner, signal.SIGKILL)
    return on_event


@_needs_subreaper
@pytest.mark.parametrize('case', ['gil-held-native-call', 'native-fork-in-owner',
                                  'native-fork-in-worker'])
def test_a_worker_dies_with_its_owner_whatever_runs_in_it(tmp_path, case):
    """r257 P2-1: the owner-death deadline holds without the worker's GIL or fd hygiene."""
    def trial():
        if case == 'gil-held-native-call':
            target = lambda fd: _scan_target(fd, tmp_path / 'b.json', _HoldsTheGil(fd), seeds=1)
            expected = 1
        elif case == 'native-fork-in-owner':
            target = lambda fd: _scan_target(fd, tmp_path / 'b.json', seeds=1, fork_helper='native')
            expected = 2
        else:
            target = lambda fd: _scan_target(fd, tmp_path / 'b.json', _ForksNatively(fd), seeds=1)
            expected = 2
        status, events, pids = _owner_trial(target, _kill_owner_when(expected))
        return {'status': status,
                'gone_after': {str(pid): _gone_within(pid, 5.) for pid in pids}}
    result = _as_subreaper(trial)
    assert result['status'] == -signal.SIGKILL, result
    assert len(result['gone_after']) == 1, result
    assert all(seconds is not None for seconds in result['gone_after'].values()), result


@_needs_subreaper
def test_a_worker_dies_with_an_owner_left_unreaped():
    """r257 bar 3: the owner is killed and stays a zombie; its worker still goes."""
    def trial():
        read_end, write_end = os.pipe()
        owner = os.fork()
        if owner == 0:
            os.close(read_end)
            try:
                _scan_target(write_end, Path(os.environ.get('TMPDIR', '/tmp')) / 'zombie.json',
                             seeds=1)
            finally:
                os._exit(0)
        os.close(write_end)
        worker = int(json.loads(os.read(read_end, 4096).split(b'\n')[0])['worker'])
        os.kill(owner, signal.SIGKILL)
        killed = time.monotonic()
        state = None
        while state != 'Z' and time.monotonic() - killed < 1.:
            # Once the owner is a zombie its children are this subreaper's.
            state = Path('/proc/{}/stat'.format(owner)).read_text().rsplit(')', 1)[1].split()[0]
        gone = _gone_within(worker, 5.)
        still_zombie = Path('/proc/{}/stat'.format(owner)).read_text().rsplit(')', 1)[1].split()[0]
        os.waitpid(owner, 0)
        os.close(read_end)
        return {'owner_state': state, 'owner_state_after': still_zombie, 'gone_after': gone}
    result = _as_subreaper(trial)
    assert result['owner_state'] == result['owner_state_after'] == 'Z', result
    assert result['gone_after'] is not None, result


class _Sigaction(ctypes.Structure):
    # glibc's struct sigaction on Linux.
    _fields_ = [('handler', ctypes.c_void_p), ('mask', ctypes.c_ulong * 16),
                ('flags', ctypes.c_int), ('restorer', ctypes.c_void_p)]


def _raw_disposition(signum):
    action = _Sigaction()
    assert ctypes.CDLL(None).sigaction(signum, None, ctypes.byref(action)) == 0
    # glibc fills only the kernel's 64-bit part of sa_mask; the rest is not state.
    return action.handler, action.flags, action.mask[0]


class _InterruptedBySigint:
    def __init__(self, fd):
        self.fd = fd

    def integrate(self, seed, power=None):
        _emit(self.fd, worker=os.getpid())
        time.sleep(30.)


def _foreign_dispositions_target(fd, output):
    """A C SIGTERM handler installed behind Python's back, and SIGINT set to restart."""
    libc = ctypes.CDLL(None)
    libc.signal.argtypes = (ctypes.c_int, ctypes.c_void_p)
    libc.signal.restype = ctypes.c_void_p
    native = ctypes.cast(libc.getpid, ctypes.c_void_p).value  # returns; harmless as a handler
    libc.signal(signal.SIGTERM, native)
    signal.signal(signal.SIGINT, signal.default_int_handler)
    signal.siginterrupt(signal.SIGINT, False)
    before = {s: _raw_disposition(s) for s in (signal.SIGINT, signal.SIGTERM)}

    class Returns:
        def integrate(self, seed, power=None):
            time.sleep(.2)
            return _complete({'u': 0., 'n_e': 1.e16, '_state_vector': [0.],
                              '_jacobian': [[-1.]]})

    rows = {}
    # 1. A real SIGTERM during the scan reaches the C handler; the scan returns.
    threading.Thread(target=lambda: (time.sleep(.1), os.kill(os.getpid(), signal.SIGTERM)),
                     daemon=True).start()
    scan(Returns(), (0., 0.), output, u_steps=1, electron_density_decades=1, processes=1,
         timeout=30., reference_power=.5, recheck_timeouts=False)
    rows['after_return'] = {s: _raw_disposition(s) == before[s] for s in before}
    # 2. The scan raises KeyboardInterrupt.
    threading.Thread(target=lambda: (time.sleep(.2), os.kill(os.getpid(), signal.SIGINT)),
                     daemon=True).start()
    try:
        scan(_InterruptedBySigint(fd), (0., 0.), output, u_steps=1,
             electron_density_decades=1, processes=1, timeout=30., reference_power=.5,
             recheck_timeouts=False)
        rows['raised'] = None
    except KeyboardInterrupt:
        rows['raised'] = 'KeyboardInterrupt'
    rows['after_raise'] = {s: _raw_disposition(s) == before[s] for s in before}
    _emit(fd, outcome=json.dumps(rows))


@_needs_subreaper
def test_foreign_and_restart_dispositions_are_byte_identical_after_a_scan(tmp_path):
    """r257 P2-2: a C handler and SA_RESTART survive a scan that returns and one that raises."""
    def trial():
        status, events, pids = _owner_trial(
            lambda fd: _foreign_dispositions_target(fd, tmp_path / 'b.json'),
            lambda event, owner, events: None)
        outcome = next((e['outcome'] for e in events if 'outcome' in e), None)
        return {'status': status, 'outcome': outcome and json.loads(outcome),
                'fates': {str(pid): _fate(pid) for pid in pids}}
    result = _as_subreaper(trial)
    assert result['status'] == 0, result
    rows = result['outcome']
    assert rows['raised'] == 'KeyboardInterrupt', result
    assert all(rows['after_return'].values()) and all(rows['after_raise'].values()), result
    assert set(result['fates'].values()) == {'reaped-by-owner'}, result


def _default_behind_pythons_back_target(fd, output):
    """SIGINT reset to SIG_DFL by C code while Python's table still says
    default_int_handler: the state RDKit's handler leaves once it fires."""
    signal.signal(signal.SIGINT, signal.default_int_handler)
    libc = ctypes.CDLL(None)
    libc.signal.argtypes = (ctypes.c_int, ctypes.c_void_p)
    libc.signal.restype = ctypes.c_void_p
    libc.signal(signal.SIGINT, None)  # SIG_DFL
    _scan_until_killed(fd, output)


def _scan_until_killed(fd, output):
    try:
        scan(_Sleeper(fd), (0., 1.), output, u_steps=2, electron_density_decades=1,
             processes=2, timeout=30., reference_power=.5, recheck_timeouts=False)
        _emit(fd, outcome='return')
    except BaseException as exc:
        _emit(fd, outcome=type(exc).__name__)


@_needs_subreaper
def test_a_default_sigint_set_behind_pythons_back_still_reaps_every_worker(tmp_path):
    """The kernel's SIG_DFL, not Python's stale table, decides: the scan defers it,
    reaps its workers, and the process then dies of SIGINT as it would have."""
    def trial():
        status, events, pids = _owner_trial(
            lambda fd: _default_behind_pythons_back_target(fd, tmp_path / 'b.json'),
            _on_workers(2, signal.SIGINT))
        return {'status': status, 'fates': {str(pid): _fate(pid) for pid in pids}}
    result = _as_subreaper(trial)
    assert result['status'] == -signal.SIGINT, result
    assert len(result['fates']) == 2 and set(result['fates'].values()) == {'reaped-by-owner'}, result


def _rdkit_embedding_target(fd, output, send_sigint):
    """A real RDKit embedding before the scan, as RMG's QM geometry code does."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    signal.signal(signal.SIGINT, signal.default_int_handler)
    AllChem.EmbedMolecule(Chem.AddHs(Chem.MolFromSmiles('CCO')), randomSeed=1)
    before = _raw_disposition(signal.SIGINT)

    class Returns:
        def integrate(self, seed, power=None):
            _emit(fd, worker=os.getpid())
            time.sleep(.4)
            return _complete({'u': seed.u, 'n_e': seed.n_e, '_state_vector': [seed.u],
                              '_jacobian': [[-1.]]})
    if send_sigint:
        threading.Thread(target=lambda: (time.sleep(.2), os.kill(os.getpid(), signal.SIGINT)),
                         daemon=True).start()
    try:
        scan(Returns(), (0., 1.), output, u_steps=2, electron_density_decades=1,
             processes=2, timeout=30., reference_power=.5, recheck_timeouts=False)
        outcome = 'return'
    except BaseException as exc:
        outcome = type(exc).__name__
    after = _raw_disposition(signal.SIGINT)
    _emit(fd, outcome=outcome, rdkit_handler=bool(before[0]) and before[0] != after[0] or None,
          identical=after == before, after_is_default=not after[0])


@_needs_subreaper
@pytest.mark.parametrize('send_sigint', [False, True], ids=['no-signal', 'sigint'])
def test_rdkits_leftover_sigint_handler_is_left_to_act_as_it_would(tmp_path, send_sigint):
    """Limit: a C handler is not deferred. Untouched by a scan; when SIGINT
    arrives, RDKit's handler swallows it and resets the signal to SIG_DFL itself."""
    pytest.importorskip('rdkit')

    def trial():
        status, events, pids = _owner_trial(
            lambda fd: _rdkit_embedding_target(fd, tmp_path / 'b.json', send_sigint),
            lambda event, owner, events: None)
        report = next((e for e in events if 'outcome' in e), {})
        return {'status': status, 'report': report,
                'fates': {str(pid): _fate(pid) for pid in pids}}
    result = _as_subreaper(trial)
    report = result['report']
    assert result['status'] == 0 and report['outcome'] == 'return', result
    assert set(result['fates'].values()) == {'reaped-by-owner'}, result
    if send_sigint:
        assert report['after_is_default'], result  # RDKit's own reset, not the scan's
    else:
        assert report['identical'], result


def _raising_handler_inside_start_target(fd, output):
    """An application's SIGUSR1 handler raises; the signal lands inside
    Process.start() right after the fork, before multiprocessing records the child."""
    class AppSignal(Exception):
        pass

    def raise_app_signal(signum, frame):
        raise AppSignal()
    signal.signal(signal.SIGUSR1, raise_app_signal)
    owner = os.getpid()
    armed = [True]

    def trace(frame, event, arg):
        # A real SIGUSR1, timed by tracing: the first line multiprocessing runs
        # in the parent after os.fork() returned, before Process.start() can
        # record the child.
        if os.getpid() != owner:
            return None
        if (armed[0] and event == 'line' and frame.f_code.co_name == '_launch'
                and getattr(frame.f_locals.get('self'), 'pid', None)):
            armed[0] = False
            os.kill(owner, signal.SIGUSR1)
        return trace
    sys.settrace(trace)
    try:
        scan(_Sleeper(fd), (0., 0.), output, u_steps=1, electron_density_decades=1,
             processes=1, timeout=30., reference_power=.5, recheck_timeouts=False)
        outcome = 'return'
    except BaseException as exc:
        outcome = type(exc).__name__
    finally:
        sys.settrace(None)
    children = Path('/proc/self/task/{}/children'.format(os.getpid())).read_text().split()
    _emit(fd, outcome=outcome, children_after_scan=children)


@_needs_subreaper
def test_a_start_interrupted_after_its_fork_still_reaps_the_child(tmp_path):
    """r257 P3 follow-up: a raising handler for an undeferred signal inside
    Process.start() must not leave the forked child unknown and unreaped."""
    def trial():
        status, events, pids = _owner_trial(
            lambda fd: _raising_handler_inside_start_target(fd, tmp_path / 'b.json'),
            lambda event, owner, events: None)
        return {'status': status, 'report': next((e for e in events if 'outcome' in e), {})}
    result = _as_subreaper(trial)
    assert result['report'].get('outcome') == 'AppSignal', result
    assert result['report']['children_after_scan'] == [], result


_APP_RESETS_FDS = r"""
import json, os, sys, threading, time
from pathlib import Path
local = threading.local(); opened = []; report_r, report_w = os.pipe()
def app_child_setup():
    # An application's own child setup, registered before the scanner:
    # drop inherited descriptors and open its own resources.
    if not getattr(local, 'app_fork', False):
        return
    for entry in list(Path('/proc/self/fd').iterdir()):
        fd = int(entry.name)
        if fd > 2 and fd != report_w:
            try:
                os.close(fd)
            except OSError:
                pass
    for _ in range(16):
        opened.append(os.open('/dev/null', os.O_WRONLY))
os.register_at_fork(after_in_child=app_child_setup)
from rmgpy.tools.eedf import branches as b
ready_r, ready_w = os.pipe()
class Slow:
    def integrate(self, s, power=None):
        os.write(ready_w, b'r'); time.sleep(2)
def run():
    b.scan(Slow(), (0, 0), sys.argv[1], u_steps=1, electron_density_decades=1, processes=1,
           timeout=.6, reference_power=.5, recheck_timeouts=False)
bad = []
for trial in range(int(sys.argv[2])):
    t = threading.Thread(target=run); t.start(); os.read(ready_r, 1)
    local.app_fork = True; pid = os.fork()
    if not pid:
        broken = []
        for fd in opened:
            try:
                os.write(fd, b'x')
            except OSError as e:
                broken.append(fd)
        os.write(report_w, (json.dumps(broken) + '\n').encode()); os._exit(0)
    local.app_fork = False; os.waitpid(pid, 0)
    bad.extend(json.loads(os.read(report_r, 65536).split(b'\n')[0])); t.join()
print(json.dumps({'closed_app_descriptors': bad}))
"""


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='reads /proc/self/fd')
def test_the_fork_hook_never_closes_a_descriptor_it_did_not_open(tmp_path):
    """r257 P2-3: an earlier application fork hook reuses fd numbers in the child."""
    script = tmp_path / 'app.py'
    script.write_text(_APP_RESETS_FDS)
    out = subprocess.run([sys.executable, str(script), str(tmp_path / 'b.json'), '3'],
                         capture_output=True, env=_repository_env(), text=True, timeout=120)
    assert out.returncode == 0, out.stderr
    assert json.loads(out.stdout.splitlines()[-1]) == {'closed_app_descriptors': []}


# r258: whose a signal is, decided from the kernel before anything is written.
_FIRST_SCAN_OWNERSHIP = r"""
import ctypes, json, os, signal, sys, threading, time
from rmgpy.tools.eedf import branches as b
mode, output, sigterm = sys.argv[1], sys.argv[2], sys.argv[3]
if sys.argv[4:] == ['stripped']:
    b._elf_function = lambda *args, **kwargs: None  # an interpreter without .symtab
libc = ctypes.CDLL(None, use_errno=True)


class Value(ctypes.Union):
    _fields_ = [('integer', ctypes.c_int), ('pointer', ctypes.c_void_p)]


class Info(ctypes.Structure):
    _fields_ = [('signo', ctypes.c_int), ('errno', ctypes.c_int), ('code', ctypes.c_int),
                ('pad', ctypes.c_int), ('pid', ctypes.c_int), ('uid', ctypes.c_uint),
                ('value', Value)]


libc.sigqueue.argtypes = (ctypes.c_int, ctypes.c_int, Value)
libc.read.argtypes = (ctypes.c_int, ctypes.c_void_p, ctypes.c_size_t)
libc.read.restype = ctypes.c_ssize_t
deliveries = []


@ctypes.CFUNCTYPE(None, ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p)
def native(signum, info, context):
    info = ctypes.cast(info, ctypes.POINTER(Info)).contents
    deliveries.append([info.code, info.value.integer])


def raw(signum):
    action = b._raw_disposition(signum)
    return [action.handler, action.flags, action.mask[0], action.restorer]


def python_sigterm(signum, frame):
    raise SystemExit(143)


# An application's policy, set before its first scan: SIGTERM ignored, left
# at its default, or given a Python handler; SIGINT owned by Python, or by
# native code behind Python's own callable.
signal.signal(signal.SIGINT, signal.default_int_handler)
signal.signal(signal.SIGTERM, {'ignored': signal.SIG_IGN, 'default': signal.SIG_DFL,
                               'python': python_sigterm}[sigterm])
if mode == 'faulthandler':
    import faulthandler
    faulthandler.register(signal.SIGINT)
elif mode == 'faulthandler-siginterrupt':
    # faulthandler's handler, in the interpreter's own image, with flags 0.
    import faulthandler
    faulthandler.register(signal.SIGINT, chain=False)
    signal.siginterrupt(signal.SIGINT, True)
elif mode == 'rdkit':
    from rdkit import Chem
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(Chem.AddHs(Chem.MolFromSmiles('CCO')), randomSeed=1)
elif mode != 'python':
    action = b._raw_disposition(signal.SIGINT)
    if mode == 'foreign-pointer':
        action.handler = ctypes.cast(libc.getpid, ctypes.c_void_p).value
    else:
        action.handler = ctypes.cast(native, ctypes.c_void_p).value
        action.flags |= 0x4 | 0x10000000  # SA_SIGINFO | SA_RESTART
    assert b._SIGACTION(signal.SIGINT, ctypes.byref(action), None) == 0
before = {signum: raw(signum) for signum in (signal.SIGINT, signal.SIGTERM)}
table = signal.getsignal(signal.SIGINT)
me = os.getpid()
row = {'changed': [], 'sent': 0, 'early': [], 'sigterm_deferred': False,
       'sigint_deferred': False, 'written_before_decision': [], 'installs': 0,
       'table_changed': []}
reader = {}
if mode == 'eintr':
    pipe_r, pipe_w = os.pipe()
    started, done = threading.Event(), threading.Event()

    def read():
        reader['tid'] = threading.get_native_id()
        started.set()
        buffer = ctypes.create_string_buffer(1)
        reader['result'] = libc.read(pipe_r, buffer, 1)
        reader['errno'] = ctypes.get_errno()
        done.set()
    thread = threading.Thread(target=read)
    thread.start()
    started.wait()
    time.sleep(.02)
setup = tuple(getattr(b, name).__code__ for name in (
    '_defer_signals', '_install_recorder', '_interpreter_handler', '_image_base',
    '_trampoline_offset', '_mapped_file', '_elf_function') if hasattr(b, name))

injected = set()


def trace(frame, event, arg):
    if os.getpid() != me:
        return None
    if event == 'call' and frame.f_code is b._install_recorder.__code__:
        row['installs'] += 1
    if event != 'line' or frame.f_code.co_filename != b.__file__:
        return trace
    if raw(signal.SIGINT) != before[signal.SIGINT]:
        row['changed'].append([frame.f_code.co_name, frame.f_lineno])
    if signal.getsignal(signal.SIGINT) is not table:
        row['table_changed'].append([frame.f_code.co_name, frame.f_lineno])
    if not row['installs'] and any(raw(s) != before[s] for s in before):
        row['written_before_decision'].append([frame.f_code.co_name, frame.f_lineno])
    if frame.f_code.co_name == '_integrate_many':
        if signal.getsignal(signal.SIGTERM) is b._record_signal:
            row['sigterm_deferred'] = True
        if signal.getsignal(signal.SIGINT) is b._record_signal:
            row['sigint_deferred'] = True
    if frame.f_code not in setup or (frame.f_code, frame.f_lineno) in injected:
        return trace
    injected.add((frame.f_code, frame.f_lineno))  # every line once, not every pass of a loop
    if mode == 'siginfo':
        assert libc.sigqueue(me, signal.SIGINT, Value(integer=12345)) == 0
        row['sent'] += 1
    elif mode == 'eintr':
        assert libc.syscall(234, me, reader['tid'], int(signal.SIGINT)) == 0  # tgkill
        row['sent'] += 1
        time.sleep(.02)
        if done.is_set():
            row['early'].append(frame.f_lineno)
    return trace


class Fails:
    def integrate(self, seed, power=None):
        raise ValueError('an ordinary failed seed')


sys.settrace(trace)
try:
    b.scan(Fails(), (0, 0), output, u_steps=1, electron_density_decades=1, processes=1,
           timeout=1, reference_power=.5, recheck_timeouts=False)
finally:
    sys.settrace(None)
if mode == 'eintr':
    os.write(pipe_w, b'x')
    thread.join(5)
    row['read'] = reader['result']
learnt = b._PYTHON_HANDLER
signal.signal(signal.SIGUSR1, python_sigterm)  # CPython's own handler, for comparison
row.update(deliveries=deliveries, learnt=learnt is not None,
           learnt_wrong=learnt not in (None, raw(signal.SIGUSR1)[0]),
           identical={str(int(s)): raw(s) == before[s] for s in before})
print(json.dumps(row))
"""


def _first_scan(tmp_path, mode, sigterm, *extra):
    script = tmp_path / 'first_scan.py'
    script.write_text(_FIRST_SCAN_OWNERSHIP)
    out = subprocess.run([sys.executable, str(script), mode, str(tmp_path / 'b.json'), sigterm,
                          *extra],
                         capture_output=True, env=_repository_env(), text=True, timeout=120)
    assert out.returncode == 0, out.stderr
    return json.loads(out.stdout.splitlines()[-1])


_BOTH_IDENTICAL = {str(int(signal.SIGINT)): True, str(int(signal.SIGTERM)): True}


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='glibc sigaction')
@pytest.mark.parametrize('sigterm', ['ignored', 'python'])
@pytest.mark.parametrize('mode', ['siginfo', 'eintr', 'foreign-pointer', 'faulthandler',
                                  'faulthandler-siginterrupt', 'rdkit'])
def test_a_native_sigint_handler_is_never_written_at_the_first_scan(tmp_path, mode, sigterm):
    """r258 P2-2 and R5: SIGINT is held by native code behind Python's callable
    (a ctypes callback, a libc function, faulthandler, RDKit's embedding) at
    the first scan of a process. It is left foreign: its kernel action and
    Python's table entry are identical at every line the scanner runs, a
    queued signal keeps its SI_QUEUE payload, and a native read under
    SA_RESTART is not interrupted. faulthandler made to look like CPython's
    handler (r260: same image, flags 0 after siginterrupt) is still told
    apart, by its symbol. With SIGTERM ignored nothing is identified or
    cached; a Python SIGTERM handler is identified, without a write, as
    CPython's own, and deferred."""
    if mode == 'rdkit':
        pytest.importorskip('rdkit')
    row = _first_scan(tmp_path, mode, sigterm)
    assert row['changed'] == [] and row['table_changed'] == [] and not row['sigint_deferred'], row
    assert not row['learnt_wrong'], row
    assert row['written_before_decision'] == [], row
    assert row['identical'] == _BOTH_IDENTICAL, row
    assert row['learnt'] == row['sigterm_deferred'] == (sigterm == 'python'), row
    if mode == 'siginfo':
        assert row['sent'] > 0 and row['deliveries'] == [[-1, 12345]] * row['sent'], row
    if mode == 'eintr':
        assert row['sent'] > 0 and len(row['deliveries']) == row['sent'], row
        assert row['early'] == [] and row['read'] == 1, row


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='glibc sigaction')
def test_python_handlers_installed_before_the_first_scan_are_deferred(tmp_path):
    """R5: no signal is at SIG_DFL (the r257 contexts probe's arrangement), yet
    both Python handlers are deferred at the first scan, with nothing written
    before the decision, and restored exactly."""
    row = _first_scan(tmp_path, 'python', 'python')
    assert row['learnt'] and not row['learnt_wrong'], row
    assert row['sigterm_deferred'] and row['sigint_deferred'], row
    assert row['written_before_decision'] == [], row
    assert row['identical'] == _BOTH_IDENTICAL, row


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='glibc sigaction')
def test_without_a_symbol_table_python_handlers_are_left_alone(tmp_path):
    """R6: an interpreter whose symbol table cannot be read (stripped) identifies
    nothing, so with Python handlers on both signals the first scan fails
    closed: nothing deferred, nothing cached, nothing written."""
    row = _first_scan(tmp_path, 'python', 'python', 'stripped')
    assert not row['learnt'] and not row['sigterm_deferred'] and not row['sigint_deferred'], row
    assert row['installs'] == 0 and row['changed'] == [] and row['table_changed'] == [], row
    assert row['identical'] == _BOTH_IDENTICAL, row


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='glibc sigaction')
def test_a_default_signal_teaches_the_python_handler_without_touching_a_native_one(tmp_path):
    """r258 R1: SIGTERM at SIG_DFL is deferred, and its install teaches CPython's
    handler address; SIGINT's native handler is still never written."""
    row = _first_scan(tmp_path, 'siginfo', 'default')
    assert row['sigterm_deferred'] and row['learnt'], row
    assert row['changed'] == [], row
    assert row['identical'] == _BOTH_IDENTICAL, row


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='glibc sigaction')
def test_the_symbol_table_locates_the_live_python_handler_and_only_it():
    """R6: the ELF reader finds signalmodule.c's signal_handler at the offset of
    the handler CPython actually installs, and not at faulthandler's."""
    import faulthandler
    from rmgpy.tools.eedf import branches
    saved = {signum: _RawSigaction() for signum in (signal.SIGUSR1, signal.SIGUSR2)}
    libc = ctypes.CDLL(None)
    for signum, action in saved.items():
        libc.sigaction(signum, None, ctypes.byref(action))
    previous = signal.signal(signal.SIGUSR1, lambda signum, frame: None)
    try:
        faulthandler.register(signal.SIGUSR2)
        python = branches._raw_disposition(signal.SIGUSR1).handler
        native = branches._raw_disposition(signal.SIGUSR2).handler
        base = branches._image_base(python)
        assert base == branches._image_base(native)
        assert branches._trampoline_offset(python) == python - base
        assert branches._trampoline_offset(native) != native - base
        path, inode = branches._mapped_file(python)
        assert os.stat(path).st_ino == inode
        assert branches._elf_function(path, inode + 1) is None  # another file is not read
    finally:
        faulthandler.unregister(signal.SIGUSR2)
        signal.signal(signal.SIGUSR1, previous)
        for signum, action in saved.items():
            libc.sigaction(signum, ctypes.byref(action), None)


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='glibc sigaction')
def test_handlers_that_disagree_are_not_taken_for_pythons(monkeypatch):
    """R6: one symbol has one address; two qualifying addresses identify nothing."""
    from rmgpy.tools.eedf import branches
    monkeypatch.setattr(branches, '_image_base', lambda address: 4096)
    monkeypatch.setattr(branches, '_trampoline_offset', lambda address: address - 4096)
    monkeypatch.setattr(branches.signal, 'getsignal', lambda signum: signal.default_int_handler)

    def action(handler):
        return SimpleNamespace(handler=handler, flags=0)
    assert branches._interpreter_handler({2: action(8192), 15: action(8192)}) == 8192
    assert branches._interpreter_handler({2: action(8192), 15: action(12288)}) is None
    monkeypatch.setattr(branches, '_trampoline_offset', lambda address: None)
    assert branches._interpreter_handler({2: action(8192)}) is None


def _stale_ignore_target(fd, output, signum):
    """Python's table says SIG_IGN; C code has set the kernel back to SIG_DFL."""
    signal.signal(signal.SIGINT, signal.default_int_handler)
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    signal.signal(signum, signal.SIG_IGN)
    libc = ctypes.CDLL(None)
    libc.signal.argtypes = (ctypes.c_int, ctypes.c_void_p)
    libc.signal.restype = ctypes.c_void_p
    libc.signal(signum, None)  # SIG_DFL
    _scan_until_killed(fd, output)


@_needs_subreaper
@pytest.mark.parametrize('signum', [signal.SIGTERM, signal.SIGINT], ids=['sigterm', 'sigint'])
def test_a_default_hidden_behind_a_cached_ignore_still_reaps_every_worker(tmp_path, signum):
    """r258 P2-3: the kernel's SIG_DFL decides, not Python's cached SIG_IGN; the
    owner reaps its workers and then dies of the signal."""
    def trial():
        status, events, pids = _owner_trial(
            lambda fd: _stale_ignore_target(fd, tmp_path / 'b.json', signum),
            _on_workers(2, signum))
        return {'status': status, 'fates': {str(pid): _fate(pid) for pid in pids}}
    result = _as_subreaper(trial)
    assert result['status'] == -signum, result
    assert len(result['fates']) == 2 and set(result['fates'].values()) == {'reaped-by-owner'}, result


def _rdkit_after_ignored_sigint_target(fd, output):
    """RMG's QM path with SIGINT ignored: RDKit's handler swallows one SIGINT and
    resets the kernel to SIG_DFL while Python's table still says SIG_IGN."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    AllChem.EmbedMolecule(Chem.AddHs(Chem.MolFromSmiles('CCO')), randomSeed=1)
    embedded = _raw_disposition(signal.SIGINT)[0]
    os.kill(os.getpid(), signal.SIGINT)
    time.sleep(.05)
    _emit(fd, primed={'rdkit_handler': bool(embedded) and embedded != 1,
                      'kernel_default': not _raw_disposition(signal.SIGINT)[0],
                      'table_ignores': signal.getsignal(signal.SIGINT) == signal.SIG_IGN})
    _scan_until_killed(fd, output)


@_needs_subreaper
def test_rdkits_reset_behind_an_ignored_sigint_still_reaps_every_worker(tmp_path):
    """r258 P2-3, real RDKit: the scan defers the kernel's SIG_DFL behind the
    stale SIG_IGN, so the owner reaps its workers before SIGINT ends it."""
    pytest.importorskip('rdkit')

    def trial():
        status, events, pids = _owner_trial(
            lambda fd: _rdkit_after_ignored_sigint_target(fd, tmp_path / 'b.json'),
            _on_workers(2, signal.SIGINT))
        primed = next((e['primed'] for e in events if 'primed' in e), None)
        return {'status': status, 'primed': primed,
                'fates': {str(pid): _fate(pid) for pid in pids}}
    result = _as_subreaper(trial)
    assert result['primed'] == {'rdkit_handler': True, 'kernel_default': True,
                                'table_ignores': True}, result
    assert result['status'] == -signal.SIGINT, result
    assert len(result['fates']) == 2 and set(result['fates'].values()) == {'reaped-by-owner'}, result


# r258 P2-1: an application at-fork child hook registered before the scanner's.
_INTERRUPTED_START = r"""
import ctypes, json, logging, os, signal, sys, threading, time
from pathlib import Path
variant, output = sys.argv[1], sys.argv[2]
local_owner = None
WAITS = ('handler-fork-in-window', 'unreadable')  # the worker is known only by its report
# An earlier hook holds the worker 12 s, past the 2 s reclaim deadline set here.
HELD = ('as-written', 'idle-thread', 'renamed', 'rename-close', 'app-fork-python',
        'app-fork-native', 'app-fork-thread-exits', 'new-thread-exits')
HOOK = 1.5 if variant in WAITS else 12


def app_child():
    # The application's own child initialization: slow, and first.
    if local_owner is not None and os.getpid() != local_owner:
        if variant in ('renamed', 'rename-close'):
            # It names the process and, in 'rename-close', drops every
            # inherited descriptor: neither name nor descriptors identify it.
            ctypes.CDLL(None).prctl(15, b'app-child', 0, 0, 0)
            if variant == 'rename-close':
                os.closerange(3, 1024)
        time.sleep(HOOK)


os.register_at_fork(after_in_child=app_child)
from rmgpy.tools.eedf import branches as b
messages = []


class Capture(logging.Handler):
    def emit(self, record):
        messages.append(record.getMessage())


logging.getLogger(b.__name__).addHandler(Capture())
from rmgpy.tools.eedf import _native_context as c
records = []
original_process = c.Owner.process
def tracked(owner, target, args):
    process = original_process(owner, target, args)
    records.append(process._native_record)
    return process
c.Owner.process = tracked
deferred, reaped = {}, []


class AppSignal(Exception):
    pass


app = {}


def fork_app(native=False):
    pid = ctypes.CDLL(None).fork() if native else os.fork()
    if pid == 0:
        time.sleep(60)
        os._exit(0)
    app['pid'] = pid


def interrupted(signum, frame):
    if variant.startswith('handler-fork') or variant == 'unreadable-before':
        # The handler forks on the starting thread, inside the start.
        fork_app(native=variant.endswith('native'))
    raise AppSignal()


class Slow:
    def integrate(self, seed, power=None):
        time.sleep(60)


def children():
    found = set()
    for task in os.listdir('/proc/self/task'):
        try:
            found.update(Path('/proc/self/task/%s/children' % task).read_text().split())
        except OSError:
            pass
    return sorted(found)


local_owner = os.getpid()
signal.signal(signal.SIGUSR1, interrupted)
if variant == 'idle-thread':
    # An idle thread that blocks nothing: any thread may take the signal.
    threading.Thread(target=time.sleep, args=(60,), daemon=True).start()
fork_now, forked = threading.Event(), threading.Event()


def app_thread():
    fork_now.wait()
    fork_app(native=variant == 'app-fork-native')
    forked.set()
    if variant not in ('app-fork-thread-exits', 'new-thread-exits'):
        time.sleep(60)


if variant.startswith('app-fork') or variant == 'unreadable':
    app['thread'] = threading.Thread(target=app_thread, daemon=True)
    app['thread'].start()
early = late = False


def trace(frame, event, arg):
    global early, late
    if os.getpid() != local_owner:
        return None
    if event != 'line' or frame.f_code.co_name != '_launch' or late:
        return trace
    if not early:
        # The first line of the start's fork step, before the worker's fork.
        early = True
        if variant == 'new-thread-exits':
            # A thread created during the start forks and exits; its child
            # moves to another thread's list.
            app['thread'] = threading.Thread(target=app_thread, daemon=True)
            app['thread'].start()
        if 'thread' in app:
            fork_now.set()
            forked.wait()
            if variant in ('app-fork-thread-exits', 'new-thread-exits'):
                app['thread'].join()
        if variant.startswith('handler-fork-before') or variant == 'unreadable-before':
            late = True
            os.kill(local_owner, signal.SIGUSR1)  # the worker is never forked
        return trace
    pid = getattr(frame.f_locals.get('self'), 'pid', None)
    if pid:
        late = True
        app['worker'] = pid
        if variant in ('renamed', 'rename-close'):
            while Path('/proc/%d/comm' % pid).read_text() != 'app-child\n':
                time.sleep(.001)
            while variant == 'rename-close' and len(os.listdir('/proc/%d/fd' % pid)) > 3:
                time.sleep(.001)
        os.kill(local_owner, signal.SIGUSR1)
    return trace


sys.settrace(trace)
started = time.monotonic()
try:
    b.scan(Slow(), (0, 0), output, u_steps=1, electron_density_decades=1, processes=1,
           timeout=60, reference_power=.5, recheck_timeouts=False)
    outcome = 'return'
except BaseException as exc:
    outcome = type(exc).__name__
finally:
    sys.settrace(None)
elapsed = time.monotonic() - started


def state(pid):
    try:
        return Path('/proc/%s/stat' % pid).read_text().rsplit(')', 1)[1].split()[0]
    except OSError:
        return 'gone'


left = children()
states = {pid: state(pid) for pid in left}
worker = str(app['worker']) if 'worker' in app else None
gone_at = time.monotonic() - started if worker is not None and state(worker) == 'gone' else None
facts = [c.native.info(record) for record in records]
assert all(not r['pid'] or r['reaped'] == 1 for r in facts), facts
final = children()
print(json.dumps({'outcome': outcome, 'fired': late, 'elapsed': elapsed,
                  'children_after_scan': left, 'states': states,
                  'app_child': str(app['pid']) if 'pid' in app else None, 'worker': worker,
                  'messages': messages, 'worker_gone_at': gone_at, 'children_final': final,
                  'deferred_reaps': reaped, 'native_records': facts,
                  'states_final': {pid: state(pid) for pid in final},
                  'reapers_after': sum(t.name == 'eedf-scan-reaper' for t in threading.enumerate())}))
for pid in final:
    os.kill(int(pid), signal.SIGKILL)
    os.waitpid(int(pid), 0)
"""
_SURVIVES = ('app-fork-python', 'app-fork-native', 'app-fork-thread-exits', 'new-thread-exits',
             'handler-fork-before-python', 'handler-fork-before-native',
             'handler-fork-in-window', 'unreadable', 'unreadable-before')
_HELD = ('as-written', 'idle-thread', 'renamed', 'rename-close', 'app-fork-python',
         'app-fork-native', 'app-fork-thread-exits', 'new-thread-exits')


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='reads /proc children')
@pytest.mark.parametrize('variant', ['as-written', 'idle-thread', 'renamed', 'rename-close',
                                     *_SURVIVES])
def test_an_interrupted_start_reaps_its_child_however_long_an_earlier_hook_runs(
        tmp_path, variant):
    """Native identity survives a real raising SIGUSR1 inside start, while an earlier application child hook sleeps, renames the worker, or closes descriptors. Ordinary children forked by Python/native code, other threads or handlers survive. The exact worker is reaped before scan raises without waiting for the hook."""
    script = tmp_path / 'interrupted.py'
    script.write_text(_INTERRUPTED_START)
    out = subprocess.run([sys.executable, str(script), variant, str(tmp_path / 'b.json')],
                         capture_output=True, env=_repository_env(), text=True, timeout=120)
    assert out.returncode == 0, out.stderr
    row = json.loads(out.stdout.splitlines()[-1])
    assert row['outcome'] == 'AppSignal' and row['fired'], row
    assert 'Exception ignored' not in out.stderr, out.stderr
    app = [row['app_child']] if variant in _SURVIVES else []
    if variant in _HELD:
        worker = row['worker']
        # N1: native identity and cancellation precede the blocked child hook.
        assert sorted(row['children_after_scan']) == sorted(app), row
        assert row['elapsed'] < 4., row
        assert row['worker_gone_at'] is not None and row['worker_gone_at'] < 4., row
        assert any(r['pid'] == int(worker) and r['reaped'] == 1 for r in row['native_records']), row
        assert row['children_final'] == app and row['reapers_after'] == 0, row
        if app:
            assert row['states_final'][row['app_child']] != 'Z', row
        return
    assert row['children_after_scan'] == app, row
    if app:
        assert row['states'][row['app_child']] != 'Z', row
    assert row['elapsed'] < 4., row  # No report or application hook is awaited.


# r267 / R8: identity is unique, a report is authoritative, cleanup is bounded.
_RECLAIM_EDGES = r"""
import ctypes, json, logging, os, resource, signal, sys, threading, time
from pathlib import Path
mode, output = sys.argv[1], sys.argv[2]
DEADLINE = 2.0
HOOK = 4.0  # 'sleeping-hook': longer than the deadline
libc = ctypes.CDLL(None)
local_owner = None
fork_now, app_ready, nproc = threading.Event(), threading.Event(), []


def app_child_hook():
    # The application's own child initialization, registered before the scanner's.
    if local_owner is None or os.getpid() == local_owner:
        return
    if mode == 'exec':
        os.execv('/bin/sleep', ['sleep', '30'])  # the worker becomes another program
    if mode in ('zombie', 'window-exited', 'window-app-unread'):
        os._exit(0)  # dies before the scanner's hook can report it
    if mode == 'two-holders':
        time.sleep(30)  # never reports within the deadline
    if mode == 'sleeping-hook':
        time.sleep(HOOK)  # reports after the deadline


def app_before():
    # 'native-fork-error' (r267): an application thread forks natively inside
    # the window, then the worker's own fork fails.
    if mode == 'native-fork-error' and local_owner == os.getpid() and 'pid' not in app:
        fork_now.set()
        app_ready.wait(5)
        nproc.append(resource.getrlimit(resource.RLIMIT_NPROC))
        resource.setrlimit(resource.RLIMIT_NPROC, (0, nproc[0][1]))


def app_after_in_parent():
    if nproc:  # the failed fork is over: a transient EAGAIN
        resource.setrlimit(resource.RLIMIT_NPROC, nproc.pop())


os.register_at_fork(before=app_before, after_in_parent=app_after_in_parent,
                    after_in_child=app_child_hook)
from rmgpy.tools.eedf import branches as b
from rmgpy.tools.eedf import _native_context as c
records = []
original_process = c.Owner.process
def tracked(owner, target, args):
    process = original_process(owner, target, args)
    records.append(process._native_record)
    return process
c.Owner.process = tracked
if mode in ('report-zombie', 'app-unread-reported'):
    # Registered after the scanner hook: native identity was published earlier.
    os.register_at_fork(after_in_child=lambda: os._exit(0) if os.getpid() != local_owner else None)
messages = []


class Capture(logging.Handler):
    def emit(self, record):
        messages.append(record.getMessage())


logging.getLogger(b.__name__).addHandler(Capture())


class AppSignal(Exception):
    pass


raised = [0]


def interrupted(signum, frame):
    raised[0] += 1
    if mode == 'window-no-worker' and 'pid' not in app:
        fork_app()  # the handler forks inside the window; the worker never is
    raise AppSignal()


class Slow:
    def integrate(self, seed, power=None):
        time.sleep(30)


app, worker, marks = {}, {}, {}


def fork_app(unreadable=False):
    pid = libc.fork()  # a native fork: no Python hooks run in it
    if pid == 0:
        if unreadable:
            libc.prctl(4, 0, 0, 0, 0)  # PR_SET_DUMPABLE 0: its memory cannot be read
        time.sleep(30)
        os._exit(0)
    app['pid'] = pid


def state(pid):
    try:
        return Path('/proc/%d/stat' % pid).read_text().rsplit(')', 1)[1].split()[0]
    except OSError:
        return 'gone'


def wait_until(condition, limit=5.):
    end = time.monotonic() + limit
    while not condition() and time.monotonic() < end:
        time.sleep(.001)


def trace(frame, event, arg):
    if os.getpid() != local_owner:
        return None
    if event != 'line':
        return trace
    code = frame.f_code
    if code is not c.NativePopen._launch.__code__ or 'interrupted' in marks:
        return trace
    line = Path(c.__file__).read_text().splitlines()[frame.f_lineno - 1].strip()
    if mode == 'app-unread-reported' and 'pid' not in app:
        fork_app(unreadable=True)  # before the worker's fork, so not a candidate
    if line == 'self.pid = native.fork(self.record)' and 'pid' not in app and mode in (
            'window-exited', 'window-app-unread', 'two-holders'):
        fork_app(unreadable=mode == 'window-app-unread')  # inside the nonce window
    if line == 'self.pid = native.fork(self.record)' and mode == 'window-no-worker':
        marks['interrupted'] = time.monotonic()
        os.kill(local_owner, signal.SIGUSR1)
        return trace
    pid = getattr(frame.f_locals.get('self'), 'pid', None)
    if not pid:
        return trace
    worker['pid'] = pid
    marks['forked'] = time.monotonic()
    if mode == 'exec':
        wait_until(lambda: os.readlink('/proc/%d/exe' % pid).endswith('/sleep'))
    if mode in ('zombie', 'report-zombie', 'app-unread-reported', 'window-exited',
                'window-app-unread'):
        wait_until(lambda: state(pid) == 'Z')
    marks['interrupted'] = time.monotonic()
    os.kill(local_owner, signal.SIGUSR1)
    return trace


if mode == 'interrupt-reclaim':
    cleanup = c.Owner.cleanup
    def interrupted_cleanup(owner, seconds=5.):
        if 'reclaim' not in marks:
            marks['reclaim'] = True
            os.kill(local_owner, signal.SIGUSR1)
        return cleanup(owner, seconds)
    c.Owner.cleanup = interrupted_cleanup
local_owner = os.getpid()
signal.signal(signal.SIGUSR1, interrupted)
if mode == 'native-fork-error':
    def native_fork_thread():
        fork_now.wait()
        fork_app()
        app_ready.set()
    threading.Thread(target=native_fork_thread, daemon=True).start()
sys.settrace(trace)
try:
    b.scan(Slow(), (0, 0), output, u_steps=1, electron_density_decades=1, processes=1,
           timeout=30, reference_power=.5, recheck_timeouts=False)
    outcome = 'return'
except BaseException as exc:
    outcome = type(exc).__name__
finally:
    sys.settrace(None)
returned = time.monotonic()
left = set()
for task in os.listdir('/proc/self/task'):
    try:
        left.update(int(p) for p in Path('/proc/self/task/%s/children' % task).read_text().split())
    except OSError:
        pass
left_states = {str(p): state(p) for p in sorted(left)}


def reapers():
    return sum(t.name == 'eedf-scan-reaper' for t in threading.enumerate())


later = {'reapers_at_return': reapers()}
if mode in ('window-no-worker', 'native-fork-error'):
    # The application child is absent from every native record and remains untouched.
    time.sleep(1.)
    later['app_state'] = state(app['pid'])
    try:
        later['app_wait'] = os.waitpid(app['pid'], os.WNOHANG)
    except ChildProcessError:
        later['app_wait'] = 'already-reaped-by-scan'
    later['reapers'] = reapers()
if mode == 'sleeping-hook':
    # Native cleanup has already reaped the worker without waiting for this hook.
    wait_until(lambda: state(worker['pid']) == 'gone', HOOK + 4)
    later['worker_gone_after_fork'] = (time.monotonic() - marks['forked']
                                       if state(worker['pid']) == 'gone' else None)
    wait_until(lambda: not reapers(), 1.)
    later['reapers'] = reapers()
print(json.dumps({'outcome': outcome, 'raised': raised[0],
                  'after_interrupt': returned - marks.get('interrupted', returned),
                  'worker': worker.get('pid'), 'app': app.get('pid'),
                  'left': left_states, 'messages': messages, 'later': later, 'native_records': [c.native.info(r) for r in records]}))
for p in left:
    try:
        os.kill(p, signal.SIGKILL)
        os.waitpid(p, 0)
    except (ProcessLookupError, ChildProcessError):
        pass
"""


def _reclaim_edge(tmp_path, mode, outcome='AppSignal'):
    script = tmp_path / 'edges.py'
    script.write_text(_RECLAIM_EDGES)
    out = subprocess.run([sys.executable, str(script), mode, str(tmp_path / 'b.json')],
                         capture_output=True, env=_repository_env(), text=True, timeout=90)
    assert out.returncode == 0, out.stderr
    assert 'Exception ignored' not in out.stderr, out.stderr  # a live daemon reaper at exit
    row = json.loads(out.stdout.splitlines()[-1])
    assert row['outcome'] == outcome, row
    return row




@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='reads /proc children')
@pytest.mark.parametrize('mode', ['report-zombie', 'app-unread-reported'])
def test_a_worker_that_dies_during_its_start_is_reaped_not_left_a_zombie(tmp_path, mode):
    """A worker exits in an early or late child hook. Its C-recorded PID remains authoritative before Python registration returns; exact cleanup leaves no worker zombie and preserves the unreadable application child."""
    row = _reclaim_edge(tmp_path, mode)
    assert str(row['worker']) not in row['left'], row
    assert 'Z' not in row['left'].values(), row
    if row['app']:
        assert row['left'] == {str(row['app']): 'S'}, row
    else:
        assert row['left'] == {}, row


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='reads /proc children')
@pytest.mark.parametrize('mode', ['exec', 'zombie', 'window-app-unread', 'window-exited',
                                  'two-holders'])
def test_an_undecided_child_is_released_and_logged_at_the_reclaim_deadline(tmp_path, mode):
    """Native-plan remapping of the inherited identity edge cases. Exec replacement ends the watchdog guarantee: bounded typed refusal and exact pending-PID logging preserve prior artifacts while the replacement remains alive. Other workers have C-recorded identity and are reaped even when they exit early or share ancestry with unreadable application children."""
    if mode == 'exec':
        # Native-plan amendment: exec removes the watchdog. Keep this inherited
        # row, checking bounded typed refusal, exact pending PID/log and both
        # absent and prior-artifact preservation. Never claim a successful reap.
        environment = dict(os.environ, PATH='')
        environment.pop('PYTHONPATH', None)
        probe = Path(__file__).with_name('_native_protocol_probe.py')
        result = subprocess.run([sys.executable, str(probe), mode], env=environment,
                                capture_output=True, text=True, timeout=20)
        assert result.returncode == 0, result.stdout + result.stderr
        assert 'NATIVE_PROTOCOL exec PASS' in result.stdout
        return
    row = _reclaim_edge(tmp_path, mode)
    assert row['after_interrupt'] < 2.0 + 1.0, row
    assert str(row['worker']) not in row['left'], row
    assert row['left'] == ({str(row['app']): 'S'} if row['app'] else {}), row
    assert any(r['pid'] == row['worker'] and r['reaped'] == 1 for r in row['native_records']), row


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='reads /proc children')
@pytest.mark.parametrize('mode', ['window-no-worker', 'native-fork-error'])
def test_an_application_child_forked_in_the_window_survives_when_no_worker_is_forked(
        tmp_path, mode):
    """An application handler/thread forks before the worker fork is interrupted or fails. The native record has PID zero, so cleanup never signals or reaps the application child. Fork failure is a typed native startup refusal."""
    row = _reclaim_edge(tmp_path, mode,
                        'NativeWorkerStartupError' if mode == 'native-fork-error' else 'AppSignal')
    assert row['worker'] is None and row['app'], row
    assert row['after_interrupt'] < 2.0 + 1.0, row
    assert row['left'] == {str(row['app']): 'S'}, row
    assert all(r['pid'] == 0 for r in row['native_records']), row
    later = row['later']
    assert later['app_state'] == 'S' and later['app_wait'] == [0, 0], row
    assert later['reapers'] == 0, row  # No report consumer or candidate discovery.


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='reads /proc children')
def test_a_worker_held_by_a_sleeping_hook_is_released_at_the_deadline_then_reaped(tmp_path):
    """The native watchdog is already live when an earlier child hook sleeps. Native cancellation and exact reap complete without waiting for a report or for the hook to return."""
    row = _reclaim_edge(tmp_path, 'sleeping-hook')
    worker = str(row['worker'])
    assert row['after_interrupt'] < 2.0 + 1.0, row
    assert row['left'] == {}, row  # Reaped before returning, despite the held hook.
    later = row['later']
    assert row['later']['reapers_at_return'] == 0, row
    assert later['worker_gone_after_fork'] is not None, row
    assert later['worker_gone_after_fork'] < 4.0 + 2.0, row
    assert later['reapers'] == 0, row


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='needs Linux prctl and /proc')
def test_an_unreadable_held_worker_is_reaped_by_its_late_report(tmp_path):
    """Native identity/reaping remains exact with unreadable child memory."""
    _native_interrupt('unreadable')


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='reads /proc children')
def test_a_later_interrupt_does_not_restart_the_cleanup_deadline(tmp_path):
    """A real raising signal after C cleanup cannot create a second budget."""
    _native_probe('exec-retry')


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='reads /proc children')
def test_an_interrupt_in_the_release_warning_preserves_it_without_extending_cleanup(tmp_path):
    """Exec-limit logging retries within the already spent native deadline."""
    _native_probe('exec-warning')


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='real Linux signal retry boundaries')
def test_release_survives_a_signal_in_the_retry_bookkeeping(tmp_path):
    """The C guard cancels and exact-reaps before public unwind when a third signal escapes Python retry bookkeeping."""
    _native_interrupt('bookkeeping')


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='fork lease controls')
def test_a_healthy_scan_renews_workers_across_several_run_leases(tmp_path, monkeypatch):
    from rmgpy.tools.eedf import branches
    class Slow(_Branches):
        def integrate(self, seed, power=None):
            time.sleep(.65)
            return super().integrate(seed, power)
    before = time.monotonic()
    results = branches._integrate_many(Slow(0.), [Seed(0., 1.e16), Seed(1., 1.e16)],
                                       processes=2, timeout=3, power=.5, lease_s=.15)
    assert time.monotonic() - before > 4 * .15
    assert len(results) == 2 and all(status == 'ok' for status, _ in results), results


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='real main-thread stall after fork')
@pytest.mark.parametrize('phase', ['after-fork', 'before-fork'])
def test_a_stalled_scan_owner_cannot_return_adopted_results(tmp_path, monkeypatch, phase):
    from rmgpy.tools.eedf import branches
    context = branches._native_context
    original = context.NativePopen._launch
    stalled, launches = [], []
    def launch(popen, process):
        original(popen, process)
        launches.append(popen.pid)
        if phase == 'after-fork':
            time.sleep(.5)
            stalled.append(True)
    monkeypatch.setattr(context.NativePopen, '_launch', launch)
    def before_fork():
        # Initial pre-fork preparation has no worker lease. Stall the second
        # launch while the first worker is already admitted and must expire.
        if phase == 'before-fork' and launches:
            time.sleep(.5)
            stalled.append(True)
    sys._eedf_application_before_fork = before_fork
    try:
        with pytest.raises(branches.ScanLeaseExpired):
            scan(_Branches(0.), (0., 1.), tmp_path / 'stalled.json', u_steps=2,
                 electron_density_decades=1, processes=2, timeout=3, lease_s=.15,
                 reference_power=.5, recheck_timeouts=False)
    finally:
        sys._eedf_application_before_fork = None
    assert stalled
    assert not (tmp_path / 'stalled.json').exists()


def test_owner_stall_at_final_native_retirement_refuses_results(tmp_path, monkeypatch):
    """The final entry must not erase an elapsed owner lease after a successful reap."""
    from rmgpy.tools.eedf import branches
    native = branches._native_context.native
    original = native.retire
    observed = []

    def stalled(record):
        facts = native.info(record)
        assert facts['phase'] == 'normal' and facts['reaped'] == 1
        observed.append(facts)
        time.sleep(.2)  # Honest application callback at the last Python/C boundary.
        return original(record)

    monkeypatch.setattr(native, 'retire', stalled)
    output = tmp_path / 'retirement.json'
    output.write_bytes(b'previous accepted artifact\n')
    before = (output.stat().st_mtime_ns, output.read_bytes())
    with pytest.raises(branches.ScanLeaseExpired):
        scan(_Branches(0.), (0., 0.), output, u_steps=1, electron_density_decades=1,
             processes=1, timeout=2., lease_s=.1, reference_power=.5)
    assert observed
    assert (output.stat().st_mtime_ns, output.read_bytes()) == before


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='real Linux signal storms')
def test_release_survives_twenty_fresh_signal_storms(tmp_path):
    """Twenty fresh owners under real 2 ms raising SIGUSR1, native exit and exact reap."""
    for _ in range(20):
        _native_interrupt('storm')


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='real raising SIGUSR1 storms')
def test_release_survives_twenty_storms_continuing_past_unwind(tmp_path):
    """The C guard cancels and exact-reaps before public unwind; each real storm continues for 200 ms afterward without another reap."""
    for _ in range(20):
        _native_interrupt('storm-unbounded')


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='Linux child subreaper')
def test_release_warning_is_written_when_the_owner_exits_immediately(tmp_path):
    """Native cancellation/reap precedes owner exit; no release-warning guarantee remains."""
    _native_interrupt('exit-immediately')


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='fork and CPython logging reinit')
def test_a_worker_can_log_while_the_prefork_consumer_holds_a_logging_lock(tmp_path):
    """An ordinary parent thread holds a handler lock across native worker fork."""
    from rmgpy.tools.eedf import branches
    entered, release = threading.Event(), threading.Event()
    owner = os.getpid()
    read_end, write_end = os.pipe()
    class HeldLog(logging.Handler):
        def emit(self, record):
            if os.getpid() == owner:
                entered.set()
                assert release.wait(10)
            else:
                os.write(write_end, b'logged')
    class Logs(_Branches):
        def integrate(self, seed, power=None):
            logging.getLogger(branches.__name__).warning('worker logging after fork')
            return super().integrate(seed, power)
    handler = HeldLog()
    logger = logging.getLogger(branches.__name__)
    logger.addHandler(handler)
    thread = threading.Thread(target=logger.warning, args=('parent holds its handler lock',))
    thread.start()
    try:
        assert entered.wait(2)
        result = scan(Logs(0.), (0., 0.), tmp_path / 'fork-log.json', u_steps=1,
                      electron_density_decades=1, processes=1, timeout=2,
                      reference_power=.5, recheck_timeouts=False)
        assert result['branches'], result
        assert select.select([read_end], [], [], 1)[0]
        assert os.read(read_end, 6) == b'logged'
        assert not release.is_set()
    finally:
        release.set()
        thread.join(2)
        logger.removeHandler(handler)
        os.close(read_end)
        os.close(write_end)




@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='reads /proc children')
def test_an_interrupt_while_reaping_the_reported_worker_still_reaps_it(tmp_path):
    """R8.4: a second raising signal arrives while cleanup kills the reported
    worker; cleanup resumes, and the worker is reaped before scan() lets the
    exception leave."""
    row = _reclaim_edge(tmp_path, 'interrupt-reclaim')
    assert row['raised'] == 2 and row['left'] == {}, row


_OWNER_DIES_IN_AN_EARLIER_HOOK = r"""
import ctypes, json, os, select, signal, sys, time
assert ctypes.CDLL(None).prctl(36, 1, 0, 0, 0) == 0  # PR_SET_CHILD_SUBREAPER
r, w = os.pipe()
owner = None


def app_hook():
    # The application's own child initialization: slow, and first.
    if owner is not None and os.getpid() != owner:
        os.write(w, b'%d\n' % os.getpid())
        time.sleep(12)
        os.write(w, b'done\n')


os.register_at_fork(after_in_child=app_hook)
from rmgpy.tools.eedf import branches as b


class Slow:
    def integrate(self, seed, power=None):
        time.sleep(60)


pid = os.fork()
if pid == 0:
    os.close(r)
    owner = os.getpid()
    try:
        b.scan(Slow(), (0, 0), sys.argv[1], u_steps=1, electron_density_decades=1,
               processes=1, timeout=60, reference_power=.5, recheck_timeouts=False)
    finally:
        os._exit(0)
os.close(w)
buffer = b''


def line(timeout):
    global buffer
    end = time.monotonic() + timeout
    while b'\n' not in buffer and time.monotonic() < end:
        if select.select([r], [], [], .05)[0]:
            buffer += os.read(r, 100)
    if b'\n' not in buffer:
        return None
    first, buffer = buffer.split(b'\n', 1)
    return first


worker = int(line(30))
os.kill(pid, signal.SIGKILL)
os.waitpid(pid, 0)
killed = time.monotonic()
# Native watchdog exits before the hook's 12-second sleep can finish.
returned = killed
gone = None
while time.monotonic() - returned < 10:
    if os.waitpid(worker, os.WNOHANG)[0]:
        gone = time.monotonic() - returned
        break
    time.sleep(.005)
if gone is None:
    os.kill(worker, signal.SIGKILL)
    os.waitpid(worker, 0)
assert line(.1) != b'done'
print(json.dumps({'hook_did_not_return': True, 'gone_after_owner_death': gone}))
"""


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='needs a Linux child subreaper')
def test_a_worker_whose_owner_died_in_an_earlier_hook_ends_within_5_s_of_it(tmp_path):
    """An earlier child hook sleeps for 12 seconds. Killing the real owner stops its worker within five seconds, before that hook returns; the controller verifies the exact reap."""
    script = tmp_path / 'owner_dies.py'
    script.write_text(_OWNER_DIES_IN_AN_EARLIER_HOOK)
    out = subprocess.run([sys.executable, str(script), str(tmp_path / 'b.json')],
                         capture_output=True, env=_repository_env(), text=True, timeout=120)
    assert out.returncode == 0, out.stderr
    row = json.loads(out.stdout.splitlines()[-1])
    assert row['hook_did_not_return'], row
    assert row['gone_after_owner_death'] is not None and row['gone_after_owner_death'] <= 5., row


def _nested_reentry_target(fd, output):
    """A SIGUSR1 handler re-enters scan(); SIGTERM arrives inside the inner scan."""
    signal.signal(signal.SIGTERM, _raise_system_exit)
    inner = []

    def monitor(signum, frame):
        inner.append(True)
        try:
            scan(_Sleeper(fd, 'inner'), (0., 0.), output, u_steps=1,
                 electron_density_decades=1, processes=1, timeout=30., reference_power=.5,
                 recheck_timeouts=False)
        finally:
            inner.pop()
    signal.signal(signal.SIGUSR1, monitor)
    try:
        scan(_Sleeper(fd), (0., 1.), output, u_steps=2, electron_density_decades=1,
             processes=2, timeout=30., reference_power=.5, recheck_timeouts=False)
        _emit(fd, outcome='return')
    except BaseException as exc:
        _emit(fd, outcome=type(exc).__name__)


@_needs_subreaper
def test_a_nested_scan_delivering_early_still_reaps_the_outer_workers_first(tmp_path):
    """r257 P3 limit: delivery at the inner checkpoint; outer workers reaped before it leaves."""
    def on_event(event, owner, events):
        if 'worker' in event and sum('worker' in e for e in events) == 2:
            os.kill(owner, signal.SIGUSR1)
        if 'inner' in event:
            time.sleep(.02)
            os.kill(owner, signal.SIGTERM)

    def trial():
        status, events, pids = _owner_trial(
            lambda fd: _nested_reentry_target(fd, tmp_path / 'b.json'), on_event)
        outcome = next((e['outcome'] for e in events if 'outcome' in e), None)
        return {'status': status, 'outcome': outcome,
                'fates': {pids[pid] + ':' + str(pid): _fate(pid) for pid in pids}}
    result = _as_subreaper(trial)
    assert result['outcome'] == 'SystemExit', result
    assert len(result['fates']) == 3, result
    assert set(result['fates'].values()) == {'reaped-by-owner'}, result


def test_matching_branch_uses_artifact_tolerances_and_refuses_nonunique_matches(tmp_path):
    path, payload = _scanned(tmp_path, _Branches(.5, .50005), (1., 1.00005), 2, (.25,))
    assert [(branch['id'], branch['u']) for branch in payload['branches']] == [
        ('reactor-000', 1.), ('reactor-001', pytest.approx(1.00005))]

    assert matching_branch(path, {'u': 1.00005, 'n_e': 1.e16}) == 'reactor-001'
    with pytest.raises(ValueError, match='no recorded branch'):
        matching_branch(path, {'u': 2., 'n_e': 1.e16})

    # Widening the tolerance after the scan is refused, not applied (review r234).
    _write(path, {**payload, 'u_tolerance': 1.e-4})
    with pytest.raises(ValueError, match='matches at tolerances its certificates were not '
                                         'judged with: certificate identity u_tolerance'):
        matching_branch(path, {'u': 1.000025, 'n_e': 1.e16})

    duplicate = copy.deepcopy(payload['branches'][0])
    duplicate['id'] = 'reactor-002'
    _write(path, {**payload, 'branches': payload['branches'] + [duplicate]})
    with pytest.raises(ValueError, match='more than one recorded branch'):
        matching_branch(path, {'u': 1., 'n_e': 1.e16}, power=.5)


def test_matching_branch_accepts_only_solved_powers_not_recorded_segments(tmp_path):
    """An older artifact's declared segment certifies nothing between solves."""
    path, payload = _scanned(tmp_path, _Branches(1.), (1.5, 1.5), 1, (.25,))
    branch = payload['branches'][0]
    powers = payload['certified_coverage']['reactor-000']['powers']
    assert powers[:2] == [.25, .28125]
    # The older producer's segment record, with a zero error bound, on the solve at .25 W.
    branch['continuation_path'][-1]['interpolation_segments'] = [{
        'certificate': 'declared-error-bound-v1', 'power': [.25, .28125],
        'u': [1.25, 1.28125], 'n_e': [1.e16, 1.e16],
        'u_error_bound': 0., 'log_ne_error_bound': 0., 'resolution': .03125}]
    _write(path, payload)

    assert matching_branch(path, {'u': 1.25, 'n_e': 1.e16}, power=.25) == 'reactor-000'
    with pytest.raises(ValueError) as refused:
        matching_branch(path, {'u': 1.265625, 'n_e': 1.e16}, power=.265625)
    assert str(refused.value) == (
        'power not certified for branch reactor-000 at 0.265625 W; certified powers: {}'
        .format(['{:.9g}'.format(value) for value in powers]))

    for key in ('continuation_path', 'continuation_path_half_step'):
        assert branch[key][-1]['power'] == .25
        branch[key][-1]['segment_certified'] = False
    _write(path, payload)
    with pytest.raises(ValueError, match='power not certified for branch reactor-000 at 0.25 W'):
        matching_branch(path, {'u': 1.25, 'n_e': 1.e16}, power=.25)


def test_matching_branch_uses_recorded_certification_points_on_curved_path(tmp_path):
    path, payload = _scanned(tmp_path, _Branches(0., shape=lambda power: power ** 2),
                             (.25, .25), 1, (.25,))
    step = payload['branches'][0]['continuation_path'][-1]
    point = step['certification_points'][0]
    assert point['power'] not in (step['power'], .5)

    assert matching_branch(path, point, power=point['power']) == 'reactor-000'
    unsolved = point['power'] - .015625 / 2.
    assert unsolved not in payload['certified_coverage']['reactor-000']['powers']
    with pytest.raises(ValueError, match='power not certified for branch reactor-000 at '
                                         '{!r} W'.format(unsolved)):
        matching_branch(path, {'u': unsolved ** 2, 'n_e': 1.e16}, power=unsolved)


def test_interpolation_estimates_are_recorded_as_uncertified_diagnostics(tmp_path):
    """The adapter's midpoint estimate is labelled, kept, and never matched on."""
    class EstimatedLinear(_Branches):
        @staticmethod
        def estimate_interpolation_error(left, middle, right):
            return {'u_error_estimate': 0., 'log_ne_error_estimate': 0.}

    path, payload = _scanned(tmp_path, EstimatedLinear(0.), (.25, .25), 1, (.5,),
                             reference_power=.25, main=(.25, .01, .05),
                             half=(.125, .005, .025))
    continuation = payload['branches'][0]['continuation_path']
    diagnostics = [item for entry in continuation
                   for item in entry.get('interpolation_diagnostics', [])]
    assert diagnostics
    assert all(item['certified'] is False and item['u_error_estimate'] == 0.
               for item in diagnostics)
    assert not any('interpolation_segments' in entry for entry in continuation)
    powers = payload['certified_coverage']['reactor-000']['powers']
    assert 'segments' not in payload['certified_coverage']['reactor-000']

    assert matching_branch(path, {'u': powers[1], 'n_e': 1.e16}, power=powers[1]) == 'reactor-000'
    assert .41 not in powers
    with pytest.raises(ValueError, match='power not certified for branch reactor-000 at 0.41 W'):
        matching_branch(path, {'u': .41, 'n_e': 1.e16}, power=.41)


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
                      'stability': main[0]['stability'],
                      'continuation_path': main,
                      'continuation_path_half_step': half}],
    }
    record_continuation(payload)
    assert payload['path_consistency']['reactor-000']['consistent'] is False
    assert payload['path_consistency']['reactor-000']['disagreeing_powers'] == [.5]
    path = tmp_path / 'branches.json'
    path.write_text(json.dumps(payload))

    for state, power in (({'u': main[-1]['u'], 'n_e': 1.e16}, .5),
                         ({'u': 1., 'n_e': 1.e16}, 0.)):
        with pytest.raises(ValueError, match=r'branch reactor-000 cannot be matched: its '
                                             r"recorded states disagree at \['0.5'\] W"):
            matching_branch(path, state, power=power)

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
                             'stability': path[0]['stability'], 'continuation_path': path,
                             'continuation_path_half_step': half_step}]}
    record_continuation(payload)
    coverage = payload['certified_coverage']['reactor-000']
    assert coverage['powers'][0] == pytest.approx(.25)
    assert .3125 in coverage['powers']
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


@pytest.mark.parametrize('reference, refusal', [
    (None, 'records no valid reference power: None'),
    ('half', "records no valid reference power: 'half'"),
    (float('nan'), 'records no valid reference power: nan'),
    (.5001, 'records reference power 0.5001 W, but its branch reactor-000 continuation_path '
            'settings starts at 0.5 W'),
])
def test_reference_power_is_validated_and_reconciled_with_the_paths_at_load(
        tmp_path, reference, refusal):
    """Review r238 P3-1: the reference power is read once, at the loader, for every route.

    Each path starts from the reference state, so its settings and initial
    state record the reference power; an edit or omission is refused before
    declaration or matching, not left to the end of a run.
    """
    path, payload = _scanned(tmp_path, _Branches(1.), (1.5, 1.5), 1, (.25,))
    payload.pop('reference_power')
    if reference is not None:
        payload['reference_power'] = reference
    _write(path, payload)
    for call in (lambda: matching_branch(path, {'u': 1.5, 'n_e': 1.e16}, power=.5),
                 lambda: matching_branch(path, {'u': 1.5, 'n_e': 1.e16}),
                 lambda: declared_branch(path, 'reactor-000')):
        with pytest.raises(ValueError, match=re.escape(refusal)):
            call()

    # Every top-level record edited together is still refused by the paths themselves.
    payload['reference_power'] = .5001
    for paths in payload['continuation_settings'].values():
        for settings in paths.values():
            settings['initial_power'] = .5001
    _write(path, payload)
    with pytest.raises(ValueError, match=re.escape(
            'records reference power 0.5001 W, but its branch reactor-000 continuation_path '
            'initial state starts at 0.5 W')):
        declared_branch(path, 'reactor-000')


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
    assert written['continuation_certificate']['name'] == 'reversible-chain-v2'
    assert written['schema'] == 2
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


class _PhysicalLinear:
    """Review r238 P2-1: u' = -(u - (1 + P)), integrated by LSODA."""

    def integrate(self, seed, power=None):
        from scipy.integrate import solve_ivp
        target = 1. + float(power)
        solution = solve_ivp(lambda _t, y: [-(y[0] - target)], (0., 100.),
                             list(seed.state[:1]) or [seed.u], method='LSODA',
                             rtol=1.e-11, atol=1.e-13)
        assert solution.success
        return _complete({'u': float(solution.y[0, -1]), 'n_e': 1.e16,
                          '_state_vector': solution.y[:, -1].tolist(), '_jacobian': [[-1.]],
                          'self_sustained': True, 'discharge_class': 'self-sustained'})


def _step_sizes(path):
    return [entry['step_size'] for entry in path[1:] if entry.get('converged')]


@pytest.mark.parametrize('target', [.4, .41, .45, .49, .5, .51, .55])
def test_a_short_tail_is_split_rather_than_taken_above_max_step(target):
    """Review r238 P2-1: snapping the tail to the target exceeded the recorded ceiling.

    Remainders strictly between max_step and max_step + min_step (0.41, 0.45
    and 0.49 W here) are taken in two equal steps; the rest are controls.
    """
    path = continue_branch(_PhysicalLinear(), {'u': 1., 'n_e': 1.e16, 'state': [1.]},
                           [target], initial_power=0., max_step=.2, min_step=.1,
                           stability_resolution=.1)

    sizes = _step_sizes(path)
    assert path[-1]['converged'] is True and path[-1]['power'] == pytest.approx(target)
    assert path[0]['continuation_settings']['max_step'] == .2
    assert max(sizes) <= .2 * (1. + 1.e-9)
    assert min(sizes) >= .1 * (1. - 1.e-9)


def test_the_first_step_of_a_leg_stays_within_max_step_at_real_settings(tmp_path):
    """Review r238 P2-1: 0.5 -> 0.53 W at the development settings took one 0.03-W step."""
    path = continue_branch(_PhysicalLinear(), {'u': 1.5, 'n_e': 1.e16, 'state': [1.5]},
                           [.53], **REAL_SETTINGS)
    assert _step_sizes(path) == pytest.approx([.015, .015])
    assert path[-1]['power'] == pytest.approx(.53)

    # Through scan(), whose half-step path the real-deck driver adds.
    output, payload = _scanned(tmp_path, _PhysicalLinear(), (1., 1.), 1, (.45,),
                               reference_power=0., main=(.2, .1, .1), half=(.1, .05, .05),
                               u_tolerance=1.e-4, log_ne_tolerance=1.e-3)
    branch = payload['branches'][0]
    assert _step_sizes(branch['continuation_path']) == pytest.approx([.2, .125, .125])
    assert max(_step_sizes(branch['continuation_path_half_step'])) <= .1 * (1. + 1.e-9)
    end = branch['continuation_path'][-1]
    assert matching_branch(output, end, power=end['power']) == 'reactor-000'


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
                             'stability': main[0]['stability'],
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


def test_scan_without_half_step_path_is_unconfirmed_and_refused(tmp_path):
    """Review r221 case L: production scan() checks nothing against a half path.

    The step path aliases onto B at its last gap (see the last-gap test). With
    no half-step path the artifact must not read consistent, and matching
    refuses both branches' states at 0.5 W instead of accepting B, and the
    reference state too (review r224 P1-2).
    """
    adapter = _ScanPair(lambda power: 1. - power ** 3, 2. * 28. / 512.)
    output = tmp_path / 'branches.json'

    payload = scan(adapter, (1., 1.), output, u_steps=1, electron_density_decades=1,
                   powers=(.5,), reference_power=0., continuation_max_step=.5,
                   continuation_min_step=.125, continuation_stability_resolution=.125)

    path = payload['branches'][0]['continuation_path']
    assert adapter.branch(path[-1]['u'], .5) == 'B'
    for u, power in ((path[-1]['u'], .5), (adapter.a(.5), .5), (1., 0.)):
        with pytest.raises(ValueError, match='no half-step path is recorded'):
            matching_branch(output, {'u': u, 'n_e': 1.e16}, power=power)
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


EDITED_TOLERANCES = {'u_tolerance': .01, 'log_ne_tolerance': .1}


@pytest.mark.parametrize('edited', ['u_tolerance', 'log_ne_tolerance', 'both',
                                    'every top-level record'])
def test_matching_tolerances_edited_after_the_scan_are_refused_not_applied(tmp_path, edited):
    """Review r234 P1-1: a settings-only edit must not widen what a certificate matched.

    The artifact is the producer's own; only top-level settings are edited,
    never a state, verdict or step certificate. The terminal state lies
    16x the scan's u tolerance and 9.5x its density tolerance from the
    solved 0.5-W state, as the reviewer's real-deck geometry variant did.
    Every top-level record edited together is still refused by the paths.
    """
    seed = {'u': 1.5, 'n_e': 1.5e16, 'state': [1.5]}
    main = continue_branch(_LinearStable(), seed, [.5, .55], **REAL_SETTINGS)
    half = continue_branch(_LinearStable(), seed, [.5, .55], **HALF_SETTINGS)
    payload = json.loads(json.dumps(_with_paths(main, half, u=1.5, n_e=1.5e16)))
    state = {'u': 1.5 + .0016, 'n_e': 1.5e16 * float(np.exp(.0095))}
    artifact = tmp_path / 'branches.json'
    artifact.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match='terminal state matches no recorded branch'):
        matching_branch(artifact, state, power=.5)

    edits = EDITED_TOLERANCES if edited in ('both', 'every top-level record') else {
        edited: EDITED_TOLERANCES[edited]}
    payload.update(edits)
    if edited == 'every top-level record':
        payload['continuation_certificate']['identity'].update(edits)
        for paths in payload['continuation_settings'].values():
            for settings in paths.values():
                settings.update(edits)
    artifact.write_text(json.dumps(payload))
    refusal = ("matches at tolerances its certificates were not judged with: .*differs from "
               "the artifact's .*; restore the scan's tolerances or re-run the scan")
    for call in (lambda: matching_branch(artifact, state, power=.5),
                 lambda: matching_branch(artifact, {'u': 1.5, 'n_e': 1.5e16}, power=.5),
                 lambda: declared_branch(artifact, 'reactor-000')):
        with pytest.raises(ValueError, match=refusal):
            call()


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


class _InvariantBasinPair(PlasmaReactorAdapter):
    """Review r222 P1-1: two stable branches 0.002 (20 u_tol) apart at every power.

    The seed's w component selects the basin (w < 0.5 settles to 0, else 1),
    and u settles to a + (b - a) w. Sampled on the 0.0625 W grid branch A is
    flat at 1.5, while between samples it dips by 0.00198 to almost B's value.
    The production midpoint estimate inherited from PlasmaReactorAdapter.
    """

    def __init__(self):
        pass

    @staticmethod
    def branches(power):
        a = 1.5 - .00198 * np.sin(16. * np.pi * power) ** 2
        return a, a + .002

    def integrate(self, seed, power=None):
        a, b = self.branches(float(power))
        w = 0. if seed.state[1] < .5 else 1.
        u = a + (b - a) * w
        jacobian = [[-2., 2. * (b - a)], [0., -(3. * w * w - 3. * w + .5)]]
        return _complete({'u': u, 'n_e': 1.e16, '_state_vector': [u, w],
                          '_jacobian': jacobian})


def test_matching_refuses_between_solved_powers_and_never_assigns_b_to_a(tmp_path):
    """Review r222 P1-1 off the solved grid: no interpolated candidate exists.

    At 0.03 W the true B state lies 3e-5 from A's straight line between its
    solves at 0 and 0.0625 W, and the adapter's midpoint estimate there is
    roundoff. Interpolating assigned B to A; matching now refuses, and at the
    solved 0.03125 W each true state still matches its own branch.
    """
    adapter = _InvariantBasinPair()
    branches = []
    for number, (u, w) in enumerate([(1.5, 0.), (1.502, 1.)]):
        seed = {'u': u, 'n_e': 1.e16, 'state': [u, w]}
        main = continue_branch(adapter, seed, [.5], initial_power=0., max_step=.125,
                               min_step=.0625, stability_resolution=.0625)
        half = continue_branch(adapter, seed, [.5], initial_power=0., max_step=.0625,
                               min_step=.03125, stability_resolution=.03125)
        assert main[-1]['power'] == .5 and half[-1]['power'] == .5
        branches.append({'id': 'reactor-{:03d}'.format(number), 'u': u, 'n_e': 1.e16,
                         'stability': main[0]['stability'],
                         'continuation_path': main, 'continuation_path_half_step': half})
    payload = record_continuation({'reference_power': 0., 'u_tolerance': 1.e-4,
                                   'log_ne_tolerance': 1.e-3, 'branches': branches})
    artifact = tmp_path / 'branches.json'
    artifact.write_text(json.dumps(payload))

    def true_state(power, w):
        a, b = adapter.branches(power)
        return {'u': a + (b - a) * w, 'n_e': 1.e16}

    with pytest.raises(ValueError) as refused:
        matching_branch(artifact, true_state(.03, 1.), power=.03)
    assert str(refused.value).startswith('power not certified for branch reactor-000 at 0.03 W')
    with pytest.raises(ValueError, match='power not certified'):
        matching_branch(artifact, true_state(.03, 0.), power=.03)
    assert matching_branch(artifact, true_state(.03125, 1.), power=.03125) == 'reactor-001'
    assert matching_branch(artifact, true_state(.03125, 0.), power=.03125) == 'reactor-000'
    estimate = branches[0]['continuation_path'][1]['interpolation_diagnostics'][0]
    assert estimate['certified'] is False and estimate['u_error_estimate'] < 1.e-13
    assert 'nothing is certified between solved powers' in CERTIFICATE_LIMITS[6]


def test_marginal_final_chain_sample_from_another_basin_stops_the_step():
    """Review r222 P1-2: the final zero-order sample lands in a marginal basin.

    x has stable A (eigenvalue -1) and B, 0.0026 above, with eigenvalue -1e-9
    (marginal); u = 1 + P on both, so the identity cannot tell them apart. On
    0.025 -> 0.05 W the predictor endpoint stays on A while the last chain
    sample, seeded from 0.0375 W, crosses into B; reversing returns to A.
    """
    gap = .0026

    class FinalSampleBasin:
        @staticmethod
        def a(power):
            return 1. - .025 * power - power * power

        def integrate(self, seed, power=None):
            power = float(power)
            a = self.a(power)
            on_b = seed.state[0] - a >= gap / 2.
            mode = -1.e-9 if on_b else -1.
            return _complete({'u': 1. + power, 'n_e': 1.e16,
                              '_state_vector': [a + (gap if on_b else 0.), 1. + power],
                              '_jacobian': [[mode, 0.], [0., -2.]]})

    path = continue_branch(FinalSampleBasin(), {'u': 1., 'n_e': 1.e16, 'state': [1., 1.]},
                           [.05], initial_power=0., max_step=.025, min_step=.0125,
                           stability_resolution=.0125)

    event = path[-1]
    assert event['converged'] is False and event['event_type'] == 'unresolved'
    assert "verdict 'marginal' of the last chain sample at 0.05 W" in event['failure']
    assert [entry['power'] for entry in path if entry.get('converged')] == [0., .025]


def test_artifact_from_the_producer_before_the_final_sample_gate_is_refused(tmp_path):
    """Review r224 P1-1: a genuine v1 artifact certified a marginal final sample.

    The file was written unedited by 7471e69cb's scan, continue_branch and
    record_continuation (review r224, spec/stale_bump.py): on both paths the
    last chain sample at 0.05 W is marginal, a verdict v1 never read, and both
    accept 0.05 W. Its paths agree, so only the certificate version tells it
    apart. The loader refuses the file, and writing this code's metadata over
    it is refused too: that cannot restore the verdict v1 discarded.
    """
    payload = json.loads(STALE_V1_ARTIFACT.read_text())
    branch = payload['branches'][0]
    assert payload['continuation_certificate']['name'] == 'reversible-chain-v1'
    assert branch['continuation_path'][-1]['power'] == .05
    assert branch['continuation_path'][-1]['converged'] is True
    assert path_consistency(branch, payload['u_tolerance'],
                            payload['log_ne_tolerance'])['consistent'] is True

    for call in (lambda: matching_branch(STALE_V1_ARTIFACT, {'u': 1.05, 'n_e': 1.e16}, power=.05),
                 lambda: matching_branch(STALE_V1_ARTIFACT, {'u': 1., 'n_e': 1.e16}, power=0.),
                 lambda: declared_branch(STALE_V1_ARTIFACT, 'reactor-000')):
        with pytest.raises(ValueError, match="predates the current continuation certificate: "
                                             "artifact certificate 'reversible-chain-v1' is not "
                                             "'reversible-chain-v2'; re-run the scan"):
            call()

    with pytest.raises(ValueError, match="branch reactor-000 continuation_path: the accepted "
                                         "state at 0.0 W was certified with 'reversible-chain-v1'"):
        record_continuation(payload)
    relabelled = tmp_path / 'branches.json'
    payload['continuation_certificate']['name'] = CONTINUATION_CERTIFICATE
    relabelled.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="branch reactor-000 continuation_path: the accepted "
                                         "state at 0.0 W was certified with 'reversible-chain-v1'"):
        matching_branch(relabelled, {'u': 1.05, 'n_e': 1.e16}, power=.05)


def test_fresh_paths_carry_the_current_certificate_on_every_accepted_state():
    """The version is written where it is checked: per state, not only at top level."""
    seed = {'u': 1.5, 'n_e': 1.5e16, 'state': [1.5]}
    path = continue_branch(_LinearStable(), seed, [.25], initial_power=.5,
                           max_step=.125, min_step=.03125, stability_resolution=.0625)

    assert path[0]['continuation_settings']['certificate'] == 'reversible-chain-v2'
    assert {entry['certificate']['name'] for entry in path[1:]} == {'reversible-chain-v2'}


class _SlowPhysicalMode:
    """Review r224 P1-2: u' = -1e-9 (u - (1 + P)), integrated by LSODA.

    Its one mode lies inside the marginal band, so scan() classifies the
    reference state marginal from its true Jacobian and continues nothing.
    """

    def integrate(self, seed, power=None):
        from scipy.integrate import solve_ivp
        target = 1. + float(power)
        solution = solve_ivp(lambda _t, y: [-1.e-9 * (y[0] - target)], (0., 1.e12),
                             [seed.u], method='LSODA', rtol=1.e-11, atol=1.e-13)
        assert solution.success
        return _complete({'u': float(solution.y[0, -1]), 'n_e': 1.e16,
                          '_state_vector': solution.y[:, -1].tolist(),
                          '_jacobian': [[-1.e-9]], 'self_sustained': True,
                          'discharge_class': 'self-sustained'})


def test_marginal_reference_branch_is_reported_but_never_matched(tmp_path):
    """Review r224 P1-2: scan() keeps a marginal branch; the reference power matched it.

    Matching at the reference power returned the raw branch whatever its
    verdict. A branch is now matchable only with a stable reference verdict
    and confirmed paths, at the reference power as anywhere else.
    """
    artifact = tmp_path / 'branches.json'
    payload = scan(_SlowPhysicalMode(), (1.5, 1.5), artifact, u_steps=1,
                   electron_density_decades=1, electron_density_limits=(1.e16, 1.e16),
                   reference_power=.5, powers=[.6], continuation_stability_resolution=.005)
    branch = payload['branches'][0]
    assert branch['stability']['classification'] == 'marginal'
    assert branch['continuation_path'] == []

    refusal = ("branch reactor-000 cannot be matched: its reference verdict is 'marginal', "
               "so it is recorded for reporting only")
    for power in (.5, None):
        with pytest.raises(ValueError, match=refusal):
            matching_branch(artifact, {'u': 1.5, 'n_e': 1.e16}, power=power)
    with pytest.raises(ValueError, match="operating branch 'reactor-000' cannot be matched: "
                                         "its reference verdict is 'marginal'"):
        declared_branch(artifact, 'reactor-000')


def test_stable_reference_branch_without_a_certified_path_is_never_matched(tmp_path):
    """Review r224 P1-2: a stable verdict alone certifies no state."""
    path, payload = _scanned(tmp_path, _Branches(1.), (1.5, 1.5), 1, (), half_step=False)
    assert payload['branches'][0]['stability']['stable'] is True
    assert payload['branches'][0]['continuation_path'] == []

    with pytest.raises(ValueError, match='branch reactor-000 cannot be matched: it has no '
                                         'certified continuation path'):
        matching_branch(path, {'u': 1.5, 'n_e': 1.e16}, power=.5)


def test_marginal_reference_verdict_withholds_a_branch_even_with_certified_paths(tmp_path):
    """Review r224 P1-2: the verdict gate holds where a caller attached paths.

    The real-deck driver continues branches itself after scan() and attaches
    both paths. Here a consistent, certified branch carries the verdict scan()
    gave the marginal reference state of the slow mode; it stays
    reporting-only even with certified states at the run power.
    """
    marginal = scan(_SlowPhysicalMode(), (1.5, 1.5), tmp_path / 'marginal.json', u_steps=1,
                    electron_density_decades=1, electron_density_limits=(1.e16, 1.e16),
                    reference_power=.5)['branches'][0]['stability']
    assert marginal['classification'] == 'marginal'
    path, payload = _scanned(tmp_path, _Branches(1.), (1.5, 1.5), 1, (.25,))
    payload['branches'][0]['stability'] = marginal
    _write(path, payload)

    with pytest.raises(ValueError, match="branch reactor-000 cannot be matched: its reference "
                                         "verdict is 'marginal'"):
        matching_branch(path, {'u': 1.5, 'n_e': 1.e16}, power=.5)


def test_branch_fixtures_carry_what_the_current_writer_records(tmp_path):
    """Review r238 P2-2: a fixture is the writer's output, not a stamped name."""
    _, payload = _scanned(tmp_path, _Branches(1.), (1.5, 1.5), 1, (.25,))
    assert payload['continuation_certificate']['identity']['u_tolerance'] == 1.e-6
    assert set(payload['continuation_settings']['reactor-000']) == {
        'continuation_path', 'continuation_path_half_step'}
    assert payload['path_consistency']['reactor-000']['consistent'] is True


def test_user_guide_states_the_certificate_and_its_limits_verbatim():
    """Review r224 P3-1: the documented limits are the recorded ones, not a paraphrase."""
    guide = Path(__file__).parents[4] / 'documentation/source/users/rmg/input.rst'
    text = ' '.join(guide.read_text().replace('\\|', '|').split())
    for claim in (CERTIFICATE_STATEMENT, *CERTIFICATE_LIMITS):
        assert ' '.join(claim.split()) in text, claim[:60]
    assert "'operatingBranch': {'id': 'reactor-000', 'path': 'branches/branches.json'}" in text
