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

"""Reactor operating-branch scans for prescribed-power EEDF reactors.

This module deliberately knows nothing about EEDF *table* branches. A caller
supplies a factory for one reactor configuration; the scanner records every
converged reactor steady state, including failed and timed-out seeds.
"""

from __future__ import annotations

import ctypes
import hashlib
import json
import math
import logging
import mmap
import multiprocessing
import os
import platform
import signal
import struct
import sys
import threading
import time
from contextlib import contextmanager, nullcontext
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from rmgpy.data.kinetics.family import install_complete_reducers
from . import _native_context
from ._native_context import (ScanLeaseExpired, NativeWorkerUnsupported,
                              NativeWorkerStartupError, NativeWorkerStatusError)


# The scan sends reactor factories and converged states through multiprocessing
# transports. Register RMG's complete molecule reducers before any child starts.
install_complete_reducers()


U_TOLERANCE = 1.0e-4
LOG_NE_TOLERANCE = 1.0e-3
ENGINE_FD_RELATIVE_STEP = 6.0e-6
MAX_PROCESSES = 6
WORKER_TERMINATE_GRACE = 0.1
# Fixed bound on owner cleanup/reaping; kernel scheduling or D state can delay
# exit observation. The worker watchdog never needs the owner's Python thread.
WORKER_KILL_TIMEOUT = 5.0
# Provisional default; R11.10's measured owner gaps decide the product value.
RUN_LEASE = 60.0
# Signals a scan on the main thread defers while it owns workers, highest
# rank first: SIGTERM outranks SIGINT.
_DEFERRED_SIGNALS = (signal.SIGTERM, signal.SIGINT)
# Reduced modes farther from zero than this multiple of the slowest mode's
# distance are assumed unable to reach zero inside one sample gap.
EIGENVALUE_BAND = 10.0
# Relative roundoff allowed when a step is split into gaps of the declared width.
SAMPLE_ROUNDOFF = 1.0e-9
CONTINUATION_PATHS = ('continuation_path', 'continuation_path_half_step')
# v2: the last chain sample's own stability verdict gates its step. Paths
# certified before that discarded the verdict, so they cannot be re-certified.
CONTINUATION_CERTIFICATE = 'reversible-chain-v2'
CERTIFICATE_STATEMENT = (
    'Each accepted step is reversible within the branch tolerances at the '
    'declared step and floor. The step is split into gaps no wider than the '
    'stability resolution. Each sample is corrected from the previous sample '
    '(zero-order); correcting it back to the previous power returns to the '
    'previous sample within (u_tolerance, log_ne_tolerance); and the last '
    'sample matches the predictor-corrected endpoint within the same '
    'tolerances. This applies to every accepted step, the first step of each '
    'leg included. Every sample and endpoint is classified stable.')
CERTIFICATE_LIMITS = (
    'The corrector integrates the reactor in time, so it selects a state by '
    'basin of attraction: the certificate never shows that no other branch '
    'exists.',
    'Branch identity is (u, ln n_e) within (u_tolerance, log_ne_tolerance), '
    'the criterion cluster_states uses; states that differ only in other '
    'components are one branch by the scan\'s definition.',
    'Aliasing: when a branch moves, within one gap, by at least its distance '
    'to the basin boundary of a neighbouring branch, a hop to that branch can '
    'be reversible at that gap. If that motion is confined to one gap, the hop '
    'is undetected: the following gaps are reversible on the neighbour, and '
    'step halving does not expose it when the motion is narrower than the '
    'half-step gap, because both paths then hop. A hop is exposed only when '
    'the step and half-step paths, or the two legs, reach different branches '
    'at a common power; matching then refuses the branch, and it refuses any '
    'branch without a recorded half-step path. '
    'aliasing_radius records each path\'s largest per-gap |delta u| and '
    '|delta ln n_e| between recorded states, absolutely and in branch '
    'tolerances; a neighbour whose basin boundary lies within it is not '
    'excluded, and the motion across a gap where a hop occurred is not '
    'recorded, so the radius does not bound such a hop.',
    'A corrector whose result does not depend on its seed makes the '
    'reversibility test vacuous. For it only the secant departure trigger '
    'remains, which detects a jump larger than departure_fraction times the '
    'local motion over the floor step. The first step of each leg has no '
    'secant, so there no trigger remains and a jump inside that step is '
    'accepted.',
    'Unstable windows narrower than the stability resolution can be missed; '
    'reduced modes farther from zero than eigenvalue_band times the slowest '
    'mode are assumed not to reach zero within one gap.',
    'A fold bracket is where the corrector stops following the branch: only '
    'failures the adapter declares as no-steady-state (is_fold_failure) '
    'narrow it, so it encloses a fold only if those failures are '
    'fold-caused; any other corrector failure ends unresolved. On the real '
    'adapter (PlasmaReactorAdapter) a simulation that exhausts its time '
    'budget without steady state, such as critical slowing down short of a '
    'fold, raises that same failure, so a time-budget stop narrows a fold '
    'bracket and can label it fold; only a wall-clock TimeoutError does not.',
    'Coverage is the solved powers only: a run is matched to a branch only at '
    'a power where the step or half-step path solved a certified state, and '
    'nothing is certified between solved powers. interpolation_diagnostics '
    'are adapter estimates of the midpoint error, recorded for inspection; '
    'they bound nothing and matching never reads them.',
)


class AdapterContractError(ValueError):
    """An adapter omitted data the scan needs; it is never replaced by a default."""


@dataclass(frozen=True)
class Seed:
    """A continuation predictor, including the complete packed reactor state."""

    u: float
    n_e: float
    state: tuple = ()


@dataclass
class ReactorRun:
    """Everything needed to execute one production ``PlasmaReactor.simulate`` call."""

    reactor: object
    core_species: list
    core_reactions: list
    model_settings: object
    simulator_settings: object
    edge_species: tuple = ()
    edge_reactions: tuple = ()
    surface_species: tuple = ()
    surface_reactions: tuple = ()
    development_unqualified: bool = False


def seed_lattice(u_limits, u_steps=12, electron_density_decades=5,
                 electron_density_limits=(1.0e14, 1.0e18)):
    """Return the design lattice (12 fields x 5 electron-density decades)."""
    lo, hi = map(float, u_limits)
    if not (np.isfinite(lo) and np.isfinite(hi) and hi >= lo and u_steps > 0):
        raise ValueError('u_limits must be finite ordered values and u_steps positive')
    nlo, nhi = map(float, electron_density_limits)
    if not (nlo > 0.0 and nhi >= nlo and electron_density_decades > 0):
        raise ValueError('electron-density limits must be ordered positive values')
    return [Seed(float(u), float(n_e))
            for u in np.linspace(lo, hi, int(u_steps))
            for n_e in np.geomspace(nlo, nhi, int(electron_density_decades))]


def _state_values(state):
    try:
        return float(state['u']), float(state['n_e'])
    except KeyError as exc:
        raise AdapterContractError(
            'adapter state has no {!r} value'.format(exc.args[0])) from exc


def _identity_gap(left, right):
    """Return |delta u| and |delta ln n_e| between two states."""
    lu, ln = _state_values(left)
    ru, rn = _state_values(right)
    log_gap = (abs(float(np.log(ln)) - float(np.log(rn)))
               if ln > 0.0 and rn > 0.0 else float('inf'))
    return abs(lu - ru), log_gap


def same_cluster(left, right, u_tolerance=U_TOLERANCE,
                 log_ne_tolerance=LOG_NE_TOLERANCE):
    """Whether two states meet the specified, non-averaging cluster criterion."""
    lu, ln = _state_values(left)
    ru, rn = _state_values(right)
    return (abs(lu - ru) < u_tolerance and ln > 0.0 and rn > 0.0
            and abs(np.log(ln) - np.log(rn)) < log_ne_tolerance)


def cluster_states(records, u_tolerance=U_TOLERANCE,
                   log_ne_tolerance=LOG_NE_TOLERANCE):
    """Group converged records without altering any member's state."""
    clusters = []
    for record in records:
        if not record.get('converged', False):
            continue
        for cluster in clusters:
            if same_cluster(cluster[0]['state'], record['state'], u_tolerance,
                            log_ne_tolerance):
                cluster.append(record)
                break
        else:
            clusters.append([record])
    return clusters


def _eigenvalues(jacobian):
    values = np.linalg.eigvals(np.asarray(jacobian, dtype=float))
    return [{'real': float(value.real), 'imag': float(value.imag)} for value in values]


def _classification(eigenvalues, zero_tolerance):
    largest = max(value['real'] for value in eigenvalues)
    if largest > zero_tolerance:
        return 'unstable'
    if largest >= -zero_tolerance:
        return 'marginal'
    return 'stable'


def _constraint_basis(jacobian, constraints, half_step_jacobian=None,
                      conservation_rtol=1.0e-7):
    """Return the tangent basis after retaining only measured left-null modes."""
    jacobian = np.asarray(jacobian, dtype=float)
    half = None if half_step_jacobian is None else np.asarray(half_step_jacobian, dtype=float)
    retained = []
    names = []
    residuals = {}
    for item in constraints:
        vector = np.asarray(item['vector'], dtype=float)
        if vector.shape != (jacobian.shape[0],):
            raise ValueError('conservation vector {!r} has shape {}, expected {}'.format(
                item['name'], vector.shape, (jacobian.shape[0],)))
        if 'verification_relative_residual' in item:
            relative = float(item['verification_relative_residual'])
        else:
            denominator = max(float(np.sum(
                np.abs(vector) * np.linalg.norm(jacobian, axis=1))),
                np.finfo(float).tiny)
            relative = float(np.linalg.norm(vector @ jacobian) / denominator)
            if half is not None:
                half_denominator = max(float(np.sum(
                    np.abs(vector) * np.linalg.norm(half, axis=1))),
                    np.finfo(float).tiny)
                relative = max(relative, float(np.linalg.norm(vector @ half) / half_denominator))
        residuals[item['name']] = relative
        verified = item.get('verified', relative <= conservation_rtol)
        if verified and relative <= conservation_rtol:
            candidate = np.asarray(retained + [vector])
            if np.linalg.matrix_rank(candidate) > len(retained):
                retained.append(vector)
                names.append(item['name'])
    if not retained:
        return np.identity(jacobian.shape[0]), names, residuals
    matrix = np.asarray(retained)
    _u, singular, right = np.linalg.svd(matrix, full_matrices=True)
    tolerance = max(matrix.shape) * singular[0] * np.finfo(float).eps
    rank = int(np.count_nonzero(singular > tolerance))
    basis = right[rank:].T
    if basis.shape[1] == 0:
        raise ValueError('conservation constraints remove the entire dynamic state')
    return basis, names, residuals


def _eliminate_algebraic_rows(jacobian, constraints, algebraic_indices):
    """Eliminate index-1 algebraic variables from a residual Jacobian."""
    jacobian = np.asarray(jacobian, dtype=float)
    indices = tuple(int(index) for index in algebraic_indices)
    if len(indices) > 1:
        raise ValueError('only one algebraic row is supported')
    if not indices:
        return jacobian, list(constraints), np.identity(jacobian.shape[0])
    algebraic = indices[0]
    if not 0 <= algebraic < jacobian.shape[0]:
        raise ValueError('algebraic row index is outside the Jacobian')
    differential = [index for index in range(jacobian.shape[0]) if index != algebraic]
    row = jacobian[algebraic]
    pivot = row[algebraic]
    scale = max(float(np.linalg.norm(row, ord=np.inf)), np.finfo(float).tiny)
    if not np.isfinite(pivot) or abs(pivot) <= np.finfo(float).eps * scale:
        raise ValueError('algebraic row cannot eliminate its state variable')
    tangent = np.zeros((jacobian.shape[0], len(differential)))
    tangent[differential, np.arange(len(differential))] = 1.0
    tangent[algebraic] = -row[differential] / pivot
    reduced = jacobian[differential] @ tangent
    reduced_constraints = []
    for item in constraints:
        vector = np.asarray(item['vector'], dtype=float)
        if vector.shape != (jacobian.shape[0],):
            raise ValueError('conservation vector {!r} has shape {}, expected {}'.format(
                item['name'], vector.shape, (jacobian.shape[0],)))
        transformed = dict(item)
        transformed['vector'] = (vector @ tangent).tolist()
        reduced_constraints.append(transformed)
    return reduced, reduced_constraints, tangent


def stability_analysis(jacobian, constraints, *, half_step_jacobian=None,
                       zero_tolerance=1.0e-8, conservation_rtol=1.0e-7,
                       algebraic_indices=()):
    """Classify a steady state after projecting out measured conservation modes.

    A near-zero reduced eigenvalue is marginal. It is never counted as stable.
    If a half-step Jacobian is supplied, both matrices must give the same verdict.
    """
    jacobian = np.asarray(jacobian, dtype=float)
    if jacobian.ndim != 2 or jacobian.shape[0] != jacobian.shape[1]:
        raise ValueError('steady-state Jacobian must be square')
    if not np.isfinite(jacobian).all():
        raise ValueError('steady-state Jacobian contains a non-finite value')
    half = None if half_step_jacobian is None else np.asarray(half_step_jacobian, dtype=float)
    if half is not None and half.shape != jacobian.shape:
        raise ValueError('half-step Jacobian shape differs from the engine Jacobian')
    dynamic, dynamic_constraints, _ = _eliminate_algebraic_rows(
        jacobian, constraints, algebraic_indices)
    half_dynamic = None
    if half is not None:
        half_dynamic, _, _ = _eliminate_algebraic_rows(
            half, constraints, algebraic_indices)
    basis, names, residuals = _constraint_basis(
        dynamic, dynamic_constraints, half_dynamic, conservation_rtol)
    reduced = basis.T @ dynamic @ basis
    values = _eigenvalues(reduced)
    full_step_verdict = _classification(values, zero_tolerance)
    half_values = None
    half_verdict = None
    if half is not None:
        half_values = _eigenvalues(basis.T @ half_dynamic @ basis)
        half_verdict = _classification(half_values, zero_tolerance)
    same_verdict = None if half is None else half_verdict == full_step_verdict
    # A step-sensitive sign is numerical uncertainty, not evidence of stability.
    verdict = ('marginal' if same_verdict is False else full_step_verdict)
    return {
        'classification': verdict,
        'full_step_classification': full_step_verdict,
        'stable': verdict == 'stable',
        'zero_tolerance': float(zero_tolerance),
        'raw_eigenvalues': _eigenvalues(jacobian),
        'algebraic_indices': list(algebraic_indices),
        'eigenvalues': values,
        'half_step_eigenvalues': half_values,
        'half_step_classification': half_verdict,
        'step_halving_same_verdict': same_verdict,
        'conservation_modes': names,
        'conservation_left_null_relative_residuals': residuals,
        'reduced_determinant': float(np.linalg.det(reduced)),
        'reduced_dimension': int(reduced.shape[0]),
        'near_zero_reduced_mode': any(
            abs(complex(value['real'], value['imag'])) <= zero_tolerance
            for value in values),
    }


def _species_constraints(core_species, reactor):
    """Build element, charge and candidate total-amount conservation covectors."""
    size = int(reactor.num_core_species) + 1
    constraints = []
    elements = set()
    counts = []
    species_index = getattr(reactor, 'species_index', None)
    for fallback_index, species in enumerate(core_species):
        molecule = getattr(species, 'molecule', None)
        value = molecule[0].get_element_count() if molecule else {}
        value = {name: number for name, number in value.items()
                 if name not in ('e', 'X')}
        index = fallback_index if species_index is None else species_index[species]
        counts.append((index, value))
        elements.update(value)
    for element in sorted(elements):
        vector = np.zeros(size)
        for index, row in counts:
            vector[index] = row.get(element, 0)
        if np.any(vector):
            constraints.append({'name': 'element:' + element, 'vector': vector.tolist()})
    charges = getattr(reactor, 'species_charges', None)
    if charges is not None and len(charges) == len(core_species):
        vector = np.zeros(size)
        for fallback_index, species in enumerate(core_species):
            index = fallback_index if species_index is None else species_index[species]
            vector[index] = charges[fallback_index]
        if np.any(vector):
            constraints.append({'name': 'charge', 'vector': vector.tolist()})
    total = np.zeros(size)
    for fallback_index, species in enumerate(core_species):
        index = fallback_index if species_index is None else species_index[species]
        total[index] = 1.0
    constraints.append({'name': 'total_amount', 'vector': total.tolist()})
    return constraints


class PlasmaReactorAdapter:
    """Drive real ``PlasmaReactor`` runs and return portable scan records.

    ``factory(seed, power)`` returns :class:`ReactorRun`. The factory owns deck
    construction at the requested seed and power; this adapter owns the solver
    call and reads the accepted state, residual, physical steady-state Jacobian,
    conservation vectors and engine discharge class.
    """

    def __init__(self, factory, finite_difference_relative_step=ENGINE_FD_RELATIVE_STEP):
        self.factory = factory
        self.finite_difference_relative_step = float(finite_difference_relative_step)
        if not self.finite_difference_relative_step > 0.0:
            raise ValueError('finite-difference relative step must be positive')

    @staticmethod
    def is_fold_failure(failure):
        """Recognize the engine's typed no-steady-state continuation outcome."""
        return failure == 'RuntimeError: reactor returned without reaching steady state'

    @staticmethod
    def estimate_interpolation_error(left, middle, right):
        """Estimate the interpolation error of a smooth branch from three solves.

        Four times the midpoint discrepancy, plus roundoff: a diagnostic for
        the artifact, not a bound. Matching never uses it.
        """
        fraction = ((float(middle['power']) - float(left['power']))
                    / (float(right['power']) - float(left['power'])))
        estimated_u = float(left['u']) + fraction * (
            float(right['u']) - float(left['u']))
        estimated_log_ne = np.log(float(left['n_e'])) + fraction * (
            np.log(float(right['n_e'])) - np.log(float(left['n_e'])))
        u_error = abs(estimated_u - float(middle['u']))
        log_ne_error = abs(estimated_log_ne - np.log(float(middle['n_e'])))
        u_roundoff = 32.0 * np.finfo(float).eps * max(
            abs(float(left['u'])), abs(float(middle['u'])),
            abs(float(right['u'])), 1.0)
        log_roundoff = 32.0 * np.finfo(float).eps * max(
            abs(np.log(float(left['n_e']))),
            abs(np.log(float(middle['n_e']))),
            abs(np.log(float(right['n_e']))), 1.0)
        return {
            'u_error_estimate': 4.0 * u_error + u_roundoff,
            'log_ne_error_estimate': 4.0 * log_ne_error + log_roundoff,
        }

    @staticmethod
    def _context(development_unqualified):
        if not development_unqualified:
            return nullcontext()
        from rmgpy.solver.eedf_provider import development_unqualified_table_route
        return development_unqualified_table_route()

    @staticmethod
    def _residual(reactor, state):
        value = reactor.residual(reactor.t, state, np.zeros_like(state))[0]
        return np.asarray(value, dtype=float).copy()

    def _finite_difference_jacobian(self, reactor, relative_step):
        state = np.asarray(reactor.y, dtype=float).copy()
        base = None
        matrix = np.zeros((len(state), len(state)))
        atols = np.asarray(reactor.atol_array, dtype=float)
        for column in range(len(state)):
            step = relative_step * max(abs(state[column]), atols[column])
            plus = state.copy()
            plus[column] += step
            if column < reactor.num_core_species and state[column] - step <= 0.0:
                if base is None:
                    base = self._residual(reactor, state)
                matrix[:, column] = (self._residual(reactor, plus) - base) / step
                continue
            minus = state.copy()
            minus[column] -= step
            matrix[:, column] = (
                self._residual(reactor, plus) - self._residual(reactor, minus)) / (2.0 * step)
        self._residual(reactor, state)
        return matrix

    def _verify_constraints(self, reactor, constraints):
        """Check conservation on contracted residuals before row cancellation."""
        state = np.asarray(reactor.y, dtype=float).copy()
        atols = np.asarray(reactor.atol_array, dtype=float)
        samples = [self._residual(reactor, state)]
        for relative_step in (self.finite_difference_relative_step,
                              self.finite_difference_relative_step / 2.0):
            for column in range(len(state)):
                step = relative_step * max(abs(state[column]), atols[column])
                plus = state.copy()
                plus[column] += step
                samples.append(self._residual(reactor, plus))
                if column >= reactor.num_core_species or state[column] - step > 0.0:
                    minus = state.copy()
                    minus[column] -= step
                    samples.append(self._residual(reactor, minus))
        self._residual(reactor, state)
        verified = []
        for item in constraints:
            vector = np.asarray(item['vector'], dtype=float)
            # Measure the conservation contraction over the perturbation ensemble.
            # A per-sample maximum is ill-conditioned at a steady state: both its
            # numerator and denominator are roundoff-sized, so that one otherwise
            # harmless base-state cancellation can reject an exact conservation
            # mode.  The norm ratio retains every sampled violation while scaling
            # against the ensemble's resolved residual signal.
            contractions = [float(vector @ residual) for residual in samples]
            scales = [float(np.sum(np.abs(vector * residual))) for residual in samples]
            relative = float(
                np.linalg.norm(contractions) /
                max(np.linalg.norm(scales), np.finfo(float).tiny))
            row = dict(item)
            row['verification_relative_residual'] = relative
            row['verified'] = relative <= 1.0e-7
            verified.append(row)
        return verified

    @staticmethod
    def _apply_continuation_state(run, seed):
        """Replace feed initialization with the predictor's complete packed state."""
        if not seed.state:
            return
        reactor = run.reactor
        state = np.asarray(seed.state, dtype=float)
        expected = len(run.core_species) + 1
        if state.shape != (expected,) or not np.isfinite(state).all():
            raise ValueError(
                'continuation state has shape {}, expected {}, with finite values'.format(
                    state.shape, (expected,)))
        if np.any(state[:len(run.core_species)] < 0.0):
            raise ValueError('continuation state has a negative species amount')
        reactor.initial_mole_fractions = {
            species: float(state[index])
            for index, species in enumerate(run.core_species)
        }
        reactor.electron_kinetics['initial_reduced_field'] = (
            float(np.exp(state[-1])), 'Td')

    def integrate(self, seed, power=None):
        run = self.factory(seed, power)
        if not isinstance(run, ReactorRun):
            raise TypeError('reactor factory must return ReactorRun')
        reactor = run.reactor
        self._apply_continuation_state(run, seed)
        if run.development_unqualified and hasattr(reactor, 'configure_development_run'):
            reactor.configure_development_run(progress_interval_seconds=60.0)
        with self._context(run.development_unqualified):
            reactor.simulate(
                list(run.core_species), list(run.core_reactions),
                list(run.edge_species), list(run.edge_reactions),
                list(run.surface_species), list(run.surface_reactions),
                model_settings=run.model_settings,
                simulator_settings=run.simulator_settings)
            if reactor.steady_state_reached is not True:
                raise RuntimeError('reactor returned without reaching steady state')
            state = np.asarray(reactor.y, dtype=float).copy()
            residual = self._residual(reactor, state)
            jacobian = np.asarray(reactor._jacobian(
                reactor.t, state, np.zeros_like(state), 0.0), dtype=float).copy()
            half = self._finite_difference_jacobian(
                reactor, self.finite_difference_relative_step / 2.0)
            constraints = self._verify_constraints(
                reactor, _species_constraints(run.core_species, reactor))
            discharge = reactor.discharge_state()
            n_e = float(reactor.energy_budget['n_e'])
        return {
            'u': float(state[reactor.te_index]),
            'n_e': n_e,
            'residual_max_abs': float(np.max(np.abs(residual))),
            'discharge_class': discharge,
            'self_sustained': discharge == 'self-sustained',
            'finite_difference_relative_step': self.finite_difference_relative_step,
            'half_finite_difference_relative_step': self.finite_difference_relative_step / 2.0,
            '_jacobian': jacobian.tolist(),
            '_half_step_jacobian': half.tolist(),
            '_conservation_constraints': constraints,
            '_algebraic_indices': (
                [int(reactor.electron_index)]
                if getattr(reactor, 'quasineutral_electron', False) else []),
            '_state_vector': state.tolist(),
        }

    @staticmethod
    def jacobian(state):
        return state['_jacobian']

    @staticmethod
    def half_step_jacobian(state):
        return state['_half_step_jacobian']

    @staticmethod
    def conservation_constraints(state):
        return state['_conservation_constraints']


def _integration_child(connection, adapter, seed, power):
    try:
        connection.send(('ok', adapter.integrate(seed, power=power)))
    except (ScanLeaseExpired, NativeWorkerUnsupported, NativeWorkerStartupError, NativeWorkerStatusError) as exc:
        connection.send(('lifecycle', exc))
    except BaseException as exc:
        connection.send(('error', '{}: {}'.format(type(exc).__name__, exc)))
    finally:
        connection.close()


_LIBC = ctypes.CDLL(None, use_errno=True) if sys.platform.startswith('linux') else None


def _worker(connection, adapter, seed, power):
    # The native watchdog is already running, including throughout Python's
    # child hooks. Restore only the calling application's signal semantics.
    _DEFERRAL.release()
    _integration_child(connection, adapter, seed, power)


# Exact signal dispositions (glibc on Linux). Deferral replaces a handler and
# must put back exactly what was there, flags (SA_RESTART) included; where
# the disposition cannot be read and restored this way, nothing is deferred.
class _Sigaction(ctypes.Structure):
    _fields_ = [('handler', ctypes.c_void_p), ('mask', ctypes.c_ulong * 16),
                ('flags', ctypes.c_int), ('restorer', ctypes.c_void_p)]


_SIGACTION = None
if _LIBC is not None and platform.libc_ver()[0] == 'glibc':
    _SIGACTION = _LIBC.sigaction
    _SIGACTION.argtypes = (ctypes.c_int, ctypes.POINTER(_Sigaction), ctypes.POINTER(_Sigaction))
    _SIGACTION.restype = ctypes.c_int
# Address of CPython's C-level signal handler, once learnt (_defer_signals).
_PYTHON_HANDLER = None
_SIG_IGN_HANDLER = 1  # SIG_IGN as a raw sa_handler value
_TRAMPOLINES = {}  # (path, inode) -> signal_handler's offset from the image base, or None


class _DlInfo(ctypes.Structure):
    _fields_ = [('fname', ctypes.c_char_p), ('fbase', ctypes.c_void_p),
                ('sname', ctypes.c_char_p), ('saddr', ctypes.c_void_p)]


def _image_base(address):
    """Load address of the object (executable or library) that contains ``address``."""
    dladdr = getattr(_LIBC, 'dladdr', None) if _SIGACTION is not None else None
    if dladdr is None:
        return None
    info = _DlInfo()
    if not dladdr(ctypes.c_void_p(address), ctypes.byref(info)):
        return None
    return info.fbase


def _mapped_file(address):
    """(path, inode) of the file this process has mapped at ``address``, or None.

    dladdr names the main program by argv[0], which may be relative or
    another file's name, so the file is found from /proc/self/maps.
    """
    try:
        with open('/proc/self/maps') as maps:
            for line in maps:
                fields = line.split(None, 5)
                low, high = (int(bound, 16) for bound in fields[0].split('-'))
                if low <= address < high:
                    return (fields[5].rstrip('\n'), int(fields[4])) if len(fields) == 6 else None
    except (OSError, ValueError, IndexError):
        pass
    return None


def _elf_function(path, inode, name=b'signal_handler', source=b'signalmodule.c'):
    """Offset from the image base of a local function, read from an ELF64 .symtab.

    The function must be the only local one called ``name``, or the only one
    listed under a ``source`` file symbol where the table names its sources
    (CPython's own build leaves those names empty). None when the file is not
    that mapping's (another inode), not ELF64, stripped, or holds no single
    such function.
    """
    try:
        with open(path, 'rb') as image:
            if os.fstat(image.fileno()).st_ino != inode:
                return None
            ident = image.read(64)
            if ident[:4] != b'\x7fELF' or ident[4] != 2 or ident[5] not in (1, 2):
                return None
            order = '<' if ident[5] == 1 else '>'
            phoff, shoff = struct.unpack_from(order + 'QQ', ident, 0x20)
            phentsize, phnum, shentsize, shnum = struct.unpack_from(order + 'HHHH', ident, 0x36)
            image.seek(phoff)
            headers = image.read(phentsize * phnum)
            loads = [struct.unpack_from(order + 'II8xQ', headers, i * phentsize) for i in range(phnum)]
            lowest = min(vaddr for kind, _, vaddr in loads if kind == 1)  # PT_LOAD
            image.seek(shoff)
            table = image.read(shentsize * shnum)
            sections = [struct.unpack_from(order + 'II16xQQI', table, i * shentsize)
                        for i in range(shnum)]
            symtab = next((section for section in sections if section[1] == 2), None)  # SHT_SYMTAB
            if symtab is None:
                return None
            image.seek(sections[symtab[4]][2])
            strings = image.read(sections[symtab[4]][3])
            image.seek(symtab[2])
            symbols = image.read(symtab[3])
    except (OSError, ValueError, IndexError, struct.error):
        return None
    found, file = {}, None
    for offset in range(0, len(symbols) - 23, 24):
        label, info, _, _, value, _ = struct.unpack_from(order + 'IBBHQQ', symbols, offset)
        symbol = strings[label:strings.find(b'\0', label)]
        if info & 0xf == 4:  # STT_FILE: the locals that follow come from it
            file = symbol
        elif info == 2 and symbol == name:  # STB_LOCAL, STT_FUNC
            found[value] = file
    named = {value for value, file in found.items() if file == source}
    values = named or set(found)
    if len(values) != 1:
        return None
    return values.pop() - (lowest & ~(mmap.PAGESIZE - 1))


def _trampoline_offset(address):
    mapped = _mapped_file(address)
    if mapped is None:
        return None
    if mapped not in _TRAMPOLINES:
        _TRAMPOLINES[mapped] = _elf_function(*mapped)
    return _TRAMPOLINES[mapped]


def _interpreter_handler(actions):
    """CPython's handler address, identified from the kernel's actions alone, or None.

    A handler function is CPython's own only if Python's table holds a
    callable for its signal, it lies in the object that exports
    PyOS_getsig (the interpreter, or libpython), and its offset there is
    that of signalmodule.c's ``signal_handler`` in the object's symbol
    table. Another function in the same object (faulthandler's, say) has
    another offset; without a symbol table nothing is identified.
    """
    try:
        interpreter = _image_base(ctypes.cast(ctypes.pythonapi.PyOS_getsig, ctypes.c_void_p).value)
    except (AttributeError, ValueError):
        return None
    if interpreter is None:
        return None
    found = set()
    for signum, action in actions.items():
        if (action is not None and action.handler and action.handler != _SIG_IGN_HANDLER
                and callable(signal.getsignal(signum))
                and _image_base(action.handler) == interpreter
                and action.handler - interpreter == _trampoline_offset(action.handler)):
            found.add(action.handler)
    # One symbol has one address: two would be a bug here, so nothing is identified.
    return found.pop() if len(found) == 1 else None


def _raw_disposition(signum):
    if _SIGACTION is None:
        return None
    action = _Sigaction()
    if _SIGACTION(signum, None, ctypes.byref(action)) != 0:
        return None
    return action


def _restore_disposition(signum, handler, action):
    # Python's table first (None, a handler set outside Python, is not
    # settable: SIG_DFL stands in), then the exact kernel action over it.
    signal.signal(signum, signal.SIG_DFL if handler is None else handler)
    _SIGACTION(signum, ctypes.byref(action), None)


class _ScanCancelled(BaseException):
    """A deferred signal whose action can run only once the scan is cleaned up.

    The outermost deferral re-delivers the signal after every worker is
    reaped. It escapes only when that signal's default action is held by a
    signal mask the caller set in every thread.
    """


class _Deferral:
    """SIGINT and SIGTERM recorded on the main thread while a scan owns workers.

    Python runs signal handlers on the main thread whichever OS thread the
    signal reached, so the recording handler defers a signal however the
    process routes it, and at any bytecode. A recorded signal is acted on
    only where every started worker is owned: a checkpoint in the scan loop
    calls the caller's own Python handler there (so KeyboardInterrupt and a
    raising SIGTERM handler unwind through the cleanup, which still defers).
    A signal whose previous disposition is SIG_DFL cancels the scan and is
    re-delivered after the cleanup, so a default SIGTERM still ends the
    process with status -15. SIGTERM outranks SIGINT: once a handler has
    raised, only a higher-ranked signal is delivered after it. Deferral
    nests: only the outermost scope installs and restores handlers; a scan
    re-entered on the same thread (from a signal handler, say) acts on a
    signal at its own checkpoint, and its exception unwinds through the
    outer scan's cleanup before it leaves ``scan()``.
    """

    def __init__(self):
        self.previous = {}
        self.pending = {}
        self.raised = None

    def checkpoint(self):
        for rank, signum in enumerate(_DEFERRED_SIGNALS):
            if signum not in self.pending:
                continue
            handler = self.previous[signum][2]
            if not callable(handler):
                raise _ScanCancelled(signum)
            self._call(rank, signum, handler, self.pending.pop(signum))

    def _call(self, rank, signum, handler, frame):
        try:
            handler(signum, frame)
        except BaseException:
            if self.raised is None or rank < self.raised:
                self.raised = rank
            raise

    def release(self):
        """Restore the caller's handlers, then deliver what is still pending."""
        taken, self.pending = self.pending, {}
        previous, self.previous = self.previous, {}
        for signum, (handler, action, _action_now) in previous.items():
            _restore_disposition(signum, handler, action)
        taken.update(self.pending)  # recorded between the swap and the restore
        self.pending = {}
        raised, self.raised = self.raised, None
        for rank, signum in enumerate(_DEFERRED_SIGNALS):
            if signum not in taken or (raised is not None and rank >= raised):
                continue
            handler = previous[signum][2]
            if callable(handler):
                self._call(rank, signum, handler, taken[signum])
            else:
                os.kill(os.getpid(), signum)


class _NoDeferral:
    """Off the main thread a scan cannot install handlers; the main thread's own govern."""

    @staticmethod
    def checkpoint():
        pass


_DEFERRAL = _Deferral()
def _record_signal(signum, frame):
    _DEFERRAL.pending[signum] = frame


@contextmanager
def _cancellation_deferred():
    """Defer SIGINT/SIGTERM from before the first worker until the last is reaped.

    Yields an object whose ``checkpoint()`` acts on deferred signals. A
    signal is left alone when it is ignored, when its disposition belongs to
    C code rather than to Python, or when it cannot be saved and restored
    exactly. Off the main thread nothing is deferred. The owner-death
    backstop covers every one of these cases.
    """
    if threading.current_thread() is not threading.main_thread():
        yield _NoDeferral
        return
    if any(signal.getsignal(signum) is _record_signal for signum in _DEFERRED_SIGNALS):
        yield _DEFERRAL  # nested: the outermost scope delivers
        return
    try:
        _DEFERRAL.pending = {}
        _DEFERRAL.raised = None
        _defer_signals()
        yield _DEFERRAL
    finally:
        _DEFERRAL.release()


def _defer_signals():
    """Install the recording handler on each deferred signal C code does not own.

    Ownership is read from the kernel before anything is written, never from
    signal.getsignal(), which reports Python's table and can be stale: a C
    library may install its own handler, or reset the signal to SIG_DFL,
    behind Python's back (RDKit's embedding leaves a SIGINT handler that does
    the latter when it fires). A kernel SIG_DFL is deferred as SIG_DFL,
    whatever the table says, so the scan reaps its workers before the signal
    ends the process. An ignored signal stays ignored. A handler function is
    Python's only if it is CPython's own C handler and Python's table holds a
    callable; any other is C code's and is never written, not even briefly.

    CPython's handler address is not public. It is learnt by reading back
    the first install on a signal whose kernel action is SIG_DFL, which the
    scan defers anyway; failing that, it is identified, without writing,
    from the actions read and the interpreter's symbol table
    (_interpreter_handler). Until one of the two succeeds no handler
    function can be told apart from C code's, so all are left alone: with
    SIGTERM ignored and SIGINT held by a ctypes callback at the first scan,
    say, or with Python handlers on both signals under a stripped
    interpreter.
    """
    global _PYTHON_HANDLER
    actions = {signum: _raw_disposition(signum) for signum in _DEFERRED_SIGNALS}
    for signum, action in actions.items():
        if action is not None and not action.handler:
            _install_recorder(signum, action, signal.SIG_DFL)
            if _PYTHON_HANDLER is None:
                learnt = _raw_disposition(signum)
                if learnt is not None and learnt.handler:
                    _PYTHON_HANDLER = learnt.handler
    if _PYTHON_HANDLER is None:
        _PYTHON_HANDLER = _interpreter_handler(actions)
    for signum, action in actions.items():
        if (action is None or not action.handler or _PYTHON_HANDLER is None
                or action.handler != _PYTHON_HANDLER):
            continue
        handler = signal.getsignal(signum)
        if callable(handler):
            _install_recorder(signum, action, handler)


def _install_recorder(signum, action, effective):
    # Each entry: Python's table value and the raw action, both restored
    # exactly, and the action the signal would have taken (delivered later).
    _DEFERRAL.previous[signum] = (signal.getsignal(signum), action, effective)
    signal.signal(signum, _record_signal)


def _after_fork_in_child():
    # Native identity is set before Python repair and cannot be inherited by an
    # unrelated fork: the native child PID must match this process's PID.
    worker = (_native_context.native is not None
              and _native_context.native.is_worker())
    _DEFERRAL.pending = {}
    _DEFERRAL.raised = None
    if worker:
        return
    for signum, (handler, action, _action_now) in _DEFERRAL.previous.items():
        if signal.getsignal(signum) is _record_signal:
            _restore_disposition(signum, handler, action)
    _DEFERRAL.previous = {}


if hasattr(os, 'register_at_fork'):
    os.register_at_fork(after_in_child=_after_fork_in_child)


def _running(process):
    return process.is_alive()


def _join(process, timeout):
    deadline = time.monotonic() + timeout
    while _running(process):
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return
        process.join(min(remaining, 0.05))


def _reap_worker(process, *, terminate=False):
    """Request native self-exit if necessary, then reap the exact recorded PID."""
    if terminate:
        process.terminate()
    _join(process, WORKER_TERMINATE_GRACE)
    if _running(process):
        process.kill()  # Private NativePopen: cancellation, never os.kill.
        _join(process, WORKER_KILL_TIMEOUT)
    if _running(process):
        raise NativeWorkerStatusError(
            'native worker pid {} has not exited within the cleanup bound'.format(process.pid))


def _integrate_many(adapter, seeds, processes, timeout, power=None, lease_s=RUN_LEASE):
    """Integrate in order with native worker lifetime and exact exit attribution.

    Worker mode requires qualified Linux x86_64/glibc 2.39/CPython 3.9, lock-free
    native records and a successful watchdog handshake. There is no fallback.
    Only the calling scan thread renews its finite ownership lease. A long
    honest callback can expire it; ScanLeaseExpired is a lifecycle cancellation,
    never a scientific failure. Raise lease_s for known long owner operations.

    The watchdog starts before Python child hooks and exits without the GIL or
    application locks. Reaping still requires the owner: while the owner is
    wedged the stopped worker can remain a zombie until the owner resumes and
    reaps it, or exits. SIGKILL/power loss may leave no owner diagnostic or
    durable artifact. Same-UID hostile tampering is out of scope. OS scheduling
    and uninterruptible waits can delay exit observation. A native startup hang
    is an error that stops further launches, not an admitted worker guarantee.

    An exec inside a worker replaces the image and ends its watchdog guarantee.
    The replacement program can run past lease expiry or owner death. RMG
    reports its exact pending PID and bounded cleanup failure, but does not
    kill that replacement or accept its results.

    A caller-thread native guard reaps or reports every started worker before
    the scan unwinds, even if signals escape the Python cleanup retry. Retaining
    the exception does not postpone that accounting. Cleanup uses one deadline;
    reporting can take longer. Application logging that blocks forever or an
    endless raising-signal storm can delay scan return indefinitely. The report
    obligation remains; interruption after a log append may produce duplicates.

    Signal deferral retains the caller's handlers/masks on the main thread;
    other signals, native handlers and off-main-thread handlers run normally.
    Helper processes created by adapters remain the application's. External
    reaping makes a result unverifiable. Acceptance requires native normal
    completion, exact reap, exit zero and a complete protocol envelope.

    One process with no timeout runs in the caller, without worker isolation.
    """
    if processes == 1 and timeout is None:
        results = []
        for seed in seeds:
            try:
                results.append(('ok', adapter.integrate(seed, power=power)))
            except (ScanLeaseExpired, NativeWorkerUnsupported, NativeWorkerStartupError,
                    NativeWorkerStatusError):
                raise
            except BaseException as exc:
                results.append(('error', '{}: {}'.format(type(exc).__name__, exc)))
        return results
    owner = _native_context.Owner(lease_s)
    with _cancellation_deferred() as deferral:
        return _native_context.native.guard(
            owner.handle,
            lambda: _integrate_owned(owner, deferral, adapter, seeds, processes, timeout, power),
            lambda wait_error: _native_cleanup_error(owner, wait_error),
            WORKER_KILL_TIMEOUT,
        )


def _native_cleanup_error(owner, wait_error=0):
    """Deliver the exact native accounting failure and return its typed error."""
    facts = owner.facts()
    pending = [r['pid'] for r in facts['records'] if r['pid'] and not r['reaped']]
    external = [r['pid'] for r in facts['records'] if r['reaped'] == 2]
    reason = 'native cleanup deadline reached' if pending else 'external reap prevents attribution'
    message = '{}; pending child PIDs {!r}; native facts: {!r}'.format(reason, pending, facts)
    if wait_error:
        message += '; waitpid errno {}'.format(wait_error)
    if external:
        message += '; externally reaped PIDs {!r}'.format(external)
    error = NativeWorkerStatusError(message)
    logging.getLogger(__name__).error(message)
    return error


def _integrate_owned(owner, deferral, adapter, seeds, processes, timeout, power):
    """Worker loop enclosed by the caller-thread native final-sweep guard."""
    native = _native_context.native
    context = _native_context.NativeForkContext()
    waiting = list(enumerate(seeds))
    active = {}
    results = [None] * len(seeds)
    try:
        while waiting or active:
            owner.checkpoint()
            deferral.checkpoint()
            while waiting and len(active) < processes:
                owner.checkpoint()
                deferral.checkpoint()
                index, seed = waiting.pop(0)
                parent, child = context.Pipe(duplex=False)
                try:
                    process = owner.process(_worker, (child, adapter, seed, power))
                    active[index] = (process, parent, time.monotonic())
                    owner.checkpoint()
                    process.start()
                    owner.checkpoint()
                except BaseException:
                    parent.close()
                    raise
                finally:
                    child.close()
            progressed = False
            for index, (process, connection, started) in list(active.items()):
                owner.checkpoint()
                if connection.poll():
                    owner.checkpoint()
                    try:
                        envelope = connection.recv()
                    except (EOFError, OSError) as exc:
                        _reap_worker(process, terminate=True)
                        owner.checkpoint()
                        facts = native.info(process._native_record)
                        if facts['reaped'] != 1:
                            raise NativeWorkerStatusError('external reap prevents attribution: {!r}'.format(facts))
                        results[index] = ('error', '{}: reactor worker receive failed after exit code {}: {}'.format(
                            type(exc).__name__, process.exitcode, exc))
                    else:
                        owner.checkpoint()
                        _reap_worker(process)
                        owner.checkpoint()
                        results[index] = _native_context.accepted(process, envelope)
                elif timeout is not None and time.monotonic() - started >= timeout:
                    _reap_worker(process, terminate=True)
                    owner.checkpoint()
                    results[index] = ('error', 'TimeoutError: exceeded {:.6g} s'.format(timeout))
                elif not _running(process):
                    _reap_worker(process)
                    owner.checkpoint()
                    facts = native.info(process._native_record)
                    if facts['reaped'] != 1:
                        raise NativeWorkerStatusError('external reap prevents attribution: {!r}'.format(facts))
                    results[index] = ('error', 'ChildProcessError: reactor worker exited without a result '
                                      '(exit code {})'.format(process.exitcode))
                else:
                    continue
                owner.checkpoint()
                if not native.retire(process._native_record):
                    owner.checkpoint()  # Typed refusal if retirement observed lease loss.
                connection.close()
                del active[index]
                progressed = True
            if active and not progressed:
                time.sleep(0.01)
    finally:
        # Fast path only: the C guard sweeps if any retry statement escapes.
        # Both paths use the same native cleanup deadline.
        interrupted = None
        while True:
            try:
                complete = owner.cleanup(WORKER_KILL_TIMEOUT)
                for process, connection, _started in active.values():
                    connection.close()
                if not complete:
                    raise _native_cleanup_error(owner)
                break
            except (NativeWorkerStatusError, ScanLeaseExpired):
                raise
            except BaseException as exc:
                interrupted = exc
            # Outside the retry catch: a copied process/thread must reach the
            # native guard's ownership refusal, not retry cleanup forever.
            owner.facts()
        if interrupted is not None:
            raise interrupted
    return results


def _public_state(state):
    return {key: value for key, value in state.items() if not key.startswith('_')}


def _required(adapter, state, key, method_name, label):
    """Read state[key] or adapter.method_name(state); refuse when neither exists."""
    if key in state:
        value = state[key]
    else:
        method = getattr(adapter, method_name, None)
        if method is None:
            raise AdapterContractError('adapter supplies no {}'.format(label))
        value = method(state)
    if value is None:
        raise AdapterContractError('adapter supplies no {}'.format(label))
    return value


def _analysis(adapter, state):
    jacobian = _required(adapter, state, '_jacobian', 'jacobian',
                         'steady-state Jacobian')
    half = _required(adapter, state, '_half_step_jacobian', 'half_step_jacobian',
                     'half-step Jacobian')
    constraints = _required(adapter, state, '_conservation_constraints',
                            'conservation_constraints', 'conservation constraints')
    algebraic = state.get('_algebraic_indices', ())
    return stability_analysis(
        jacobian, constraints, half_step_jacobian=half,
        algebraic_indices=algebraic)


def scan(adapter, u_limits, output, *, u_steps=12, electron_density_decades=5,
         electron_density_limits=(1.0e14, 1.0e18), powers=(), processes=1,
         timeout=None, u_tolerance=U_TOLERANCE,
         reference_power=None, continuation_max_step=None,
         continuation_min_step=None, continuation_stability_resolution=None,
         log_ne_tolerance=LOG_NE_TOLERANCE, recheck_timeouts=True,
         lease_s=RUN_LEASE):
    """Scan an adapter and write a portable ``branches.json`` report.

    Every seed is solved at ``reference_power``, which is required: an
    artifact without it cannot be matched at a run's power. A seed that times
    out is solved once more when ``recheck_timeouts`` is set, and the attempt
    records both outcomes. Continuing to ``powers`` needs a declared
    ``continuation_stability_resolution``. The scan continues each branch on
    the step path only; until a caller records a half-step path,
    ``path_consistency`` reports the branch unconfirmed and matching refuses
    it. Every converged cluster is recorded, but only a branch with a stable
    reference verdict and a consistent certified path can be matched
    (``matchable_issue``).

    With ``processes`` above one or a ``timeout``, workers require the qualified
    native fork launcher (Linux x86_64, glibc 2.39, CPython 3.9, lock-free records
    and a successful watchdog handshake). Unsupported modes raise
    ``NativeWorkerUnsupported``; failed admission raises ``NativeWorkerStartupError``.
    There is no legacy fallback. One process without a timeout runs in the caller.

    ``lease_s`` is finite and positive (default 60 seconds). Only the calling
    scan thread renews it at checkpoints; long worker integration is allowed
    while the caller polls. A long honest owner callback or unpickle can expire
    the lease; choose a larger value for known long owner operations.
    ``ScanLeaseExpired`` is a distinct cancellation, never a scientific failure.
    These errors write no new artifact and leave any earlier output file
    unchanged; that file is not a result of the failed call.

    A wedged owner may leave stopped workers as zombies until it resumes and
    reaps or exits. SIGKILL/power loss may leave no owner-side diagnostic or
    durable output. Same-UID hostile tampering is out of scope. The native
    watchdog does not use Python, the GIL, logging or application locks; OS
    scheduling and uninterruptible kernel waits can still delay exit observation.
    Native startup hangs stop further launches rather than promising cleanup.
    External reaping prevents result acceptance. On the main thread SIGINT and
    SIGTERM are deferred until cleanup, then retain the caller's semantics.
    An exec inside a worker replaces its image and ends the watchdog guarantee.
    The replacement program can run past lease expiry or owner death. RMG reports
    its exact pending PID and bounded cleanup failure in NativeWorkerStatusError
    and the log, but does not kill it or accept results from it. A terminal-state
    timestamp is not a reap. A native guard reaps or reports every started worker
    before returning to caller code, including when Python cleanup is interrupted.
    Permanently blocked application logging or an endless raising-signal storm
    can delay return indefinitely; the reporting obligation is never dropped.
    Retried reporting may emit duplicate lines. See ``_integrate_many`` for the
    complete limits.
    """
    if not isinstance(processes, int) or isinstance(processes, bool) or not 1 <= processes <= MAX_PROCESSES:
        raise ValueError('processes must be an integer from 1 to at most 6')
    if not math.isfinite(lease_s) or lease_s <= 0.0:
        raise ValueError('lease_s must be finite and positive')
    if timeout is not None and (not math.isfinite(timeout) or timeout <= 0.0):
        raise ValueError('timeout must be finite and positive')
    try:
        reference_power = float(reference_power)
    except (TypeError, ValueError) as exc:
        raise ValueError('scan needs a finite reference_power') from exc
    if not math.isfinite(reference_power):
        raise ValueError('scan needs a finite reference_power')
    if powers and continuation_stability_resolution is None:
        raise ValueError(
            'scan continuation needs a declared continuation_stability_resolution; '
            'the certificate solves twice per gap of that width')
    seeds = seed_lattice(u_limits, u_steps, electron_density_decades,
                         electron_density_limits)
    outcomes = _integrate_many(
        adapter, seeds, processes, timeout, power=reference_power, lease_s=lease_s)
    first_failures = {}
    if recheck_timeouts:
        timed_out = [index for index, outcome in enumerate(outcomes)
                     if outcome[0] != 'ok' and outcome[1].startswith('TimeoutError')]
        if timed_out:
            repeated = _integrate_many(
                adapter, [seeds[index] for index in timed_out], processes,
                timeout, power=reference_power, lease_s=lease_s)
            for index, outcome in zip(timed_out, repeated):
                first_failures[index] = outcomes[index][1]
                outcomes[index] = outcome
    records = []
    attempts = []
    for index, (seed, outcome) in enumerate(zip(seeds, outcomes)):
        if outcome[0] == 'ok':
            state = dict(outcome[1])
            state['u'], state['n_e'] = _state_values(state)
            record = {'seed': asdict(seed), 'converged': True, 'state': state}
            attempt = {'seed': asdict(seed), 'converged': True,
                       'state': _public_state(state)}
        else:
            record = {'seed': asdict(seed), 'converged': False, 'failure': outcome[1]}
            attempt = record.copy()
        if index in first_failures:
            attempt['first_attempt_failure'] = first_failures[index]
            attempt['recheck'] = {
                'converged': outcome[0] == 'ok',
                'failure': None if outcome[0] == 'ok' else outcome[1],
            }
        records.append(record)
        attempts.append(attempt)
    branches = []
    for number, cluster in enumerate(cluster_states(
            records, u_tolerance, log_ne_tolerance)):
        state = dict(cluster[0]['state'])
        stability = _analysis(adapter, state)
        discharge = state.get('discharge_class')
        if discharge is None:
            discharge = getattr(adapter, 'classify_discharge', lambda value: None)(state)
        sustained = bool(state.get('self_sustained', discharge == 'self-sustained'))
        branch = {
            'id': 'reactor-{:03d}'.format(number),
            'seed': cluster[0]['seed'],
            'u': state['u'],
            'n_e': state['n_e'],
            'stability': stability,
            'self_sustained': sustained,
            'discharge_class': discharge,
            'members': len(cluster),
            'continuation_path': [],
        }
        if stability['stable'] and sustained and powers:
            converged_seed = {
                'u': state['u'],
                'n_e': state['n_e'],
                'state': _continuation_state(state).tolist(),
            }
            branch['continuation_path'] = continue_branch(
                adapter, converged_seed, powers, timeout=timeout,
                initial_power=reference_power,
                max_step=continuation_max_step,
                min_step=continuation_min_step,
                stability_resolution=continuation_stability_resolution,
                u_tolerance=u_tolerance,
                log_ne_tolerance=log_ne_tolerance, lease_s=lease_s)
        branches.append(branch)
    payload = {
        'schema': 2,
        'reference_power': reference_power,
        'u_tolerance': float(u_tolerance),
        'log_ne_tolerance': float(log_ne_tolerance),
        'attempts': attempts,
        'timeout_rechecks': sorted(first_failures),
        'branches': branches,
    }
    record_continuation(payload)
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + '\n')
    return payload


def _certified_states(path):
    """Converged, certified (power, u, n_e) records of one path, samples included."""
    states = []
    for entry in path or []:
        if not (entry.get('converged') and entry.get('segment_certified', False)):
            continue
        for point in entry.get('certification_points', []):
            if point.get('converged', True) and point.get('segment_certified', True):
                states.append(point)
        states.append(entry)
    return [entry for entry in states
            if 'power' in entry and 'u' in entry and 'n_e' in entry]


def _branch_certified_states(branch):
    """Certified states of the step and half-step paths of one branch."""
    states = []
    for key in CONTINUATION_PATHS:
        states.extend(_certified_states(branch.get(key)))
    return states


def _distinct_powers(powers):
    """Sorted powers with values equal within the sample roundoff merged."""
    distinct = []
    for power in sorted(float(value) for value in powers):
        if not distinct or abs(power - distinct[-1]) > SAMPLE_ROUNDOFF * max(abs(power), 1.0):
            distinct.append(power)
    return distinct


def certified_coverage(branch):
    """Powers at which a branch has solved, certified states.

    Both paths contribute; powers equal within the sample roundoff, such as
    one power reached on the outgoing and the return leg, are listed once.
    Nothing between these powers is certified (``CERTIFICATE_LIMITS[6]``).
    """
    return {
        'paths': [key for key in CONTINUATION_PATHS if branch.get(key)],
        'powers': _distinct_powers(
            entry['power'] for entry in _branch_certified_states(branch)),
    }

def path_consistency(branch, u_tolerance, log_ne_tolerance):
    """Check that every recorded state of a branch at one power is one state.

    States from the step and half-step paths, and from the outgoing and return
    legs of each, are grouped by power; each group must be one cluster under
    the scan's own tolerances. A group that splits is a hop one path took and
    another did not. Every certified power of the step path must also lie
    within the half-step path's certified range: a step state the half-step
    path never reached is unconfirmed. Without a recorded half-step path
    nothing was cross-checked, so every step state is unconfirmed and the
    branch is never consistent.
    """
    records = _branch_certified_states(branch)
    half_recorded = bool(branch.get(CONTINUATION_PATHS[1]))
    half_powers = [float(entry['power'])
                   for entry in _certified_states(branch.get(CONTINUATION_PATHS[1]))]
    unconfirmed = []
    for entry in _certified_states(branch.get(CONTINUATION_PATHS[0])):
        power = float(entry['power'])
        margin = SAMPLE_ROUNDOFF * max(abs(power), 1.0)
        if not half_powers or not (min(half_powers) - margin <= power
                                   <= max(half_powers) + margin):
            unconfirmed.append(power)
    records.sort(key=lambda entry: float(entry['power']))
    groups = []
    for entry in records:
        power = float(entry['power'])
        if groups and abs(power - groups[-1][0]) <= SAMPLE_ROUNDOFF * max(abs(power), 1.0):
            groups[-1][1].append(entry)
        else:
            groups.append((power, [entry]))
    compared = 0
    largest_u = 0.0
    largest_log_ne = 0.0
    disagreements = []
    for power, group in groups:
        if len(group) < 2:
            continue
        compared += 1
        for left_index, left in enumerate(group):
            for right in group[left_index + 1:]:
                u_gap, log_gap = _identity_gap(left, right)
                largest_u = max(largest_u, u_gap)
                largest_log_ne = max(largest_log_ne, log_gap)
                if not same_cluster(left, right, u_tolerance, log_ne_tolerance):
                    disagreements.append(power)
    return {
        'compared_powers': compared,
        'max_u_difference': largest_u,
        'max_log_ne_difference': largest_log_ne,
        'disagreeing_powers': sorted(set(disagreements)),
        'unconfirmed_by_half_step': _distinct_powers(unconfirmed),
        'half_step_recorded': half_recorded,
        'consistent': half_recorded and not disagreements and not unconfirmed,
    }


def aliasing_radius(path, u_tolerance, log_ne_tolerance):
    """Largest per-gap motion between consecutive recorded states of a path.

    Consecutive states are the chain samples and endpoints in path order. A
    neighbour whose basin boundary lies closer than this motion can alias at
    one gap (``CERTIFICATE_LIMITS[2]``).
    """
    states = _certified_states(path)
    largest_u = 0.0
    largest_log_ne = 0.0
    for left, right in zip(states, states[1:]):
        u_gap, log_gap = _identity_gap(left, right)
        largest_u = max(largest_u, u_gap)
        largest_log_ne = max(largest_log_ne, log_gap)
    return {
        'gaps': max(len(states) - 1, 0),
        'max_u_motion': largest_u,
        'max_log_ne_motion': largest_log_ne,
        'max_u_motion_in_u_tolerances': largest_u / u_tolerance,
        'max_log_ne_motion_in_log_ne_tolerances': largest_log_ne / log_ne_tolerance,
    }


def _path_certificate_issue(path):
    """Name the first accepted entry of a path not certified by this code's certificate.

    A step records its certificate; the initial state records the settings
    of the run that certified the path. An accepted entry with neither, or
    with an older certificate, was not certified with the last-sample gate.
    """
    for entry in path:
        if not (entry.get('converged') and entry.get('segment_certified', False)):
            continue
        record = entry.get('certificate') or entry.get('continuation_settings') or {}
        name = record.get('name', record.get('certificate'))
        if name != CONTINUATION_CERTIFICATE:
            return 'the accepted state at {!r} W was certified with {!r}'.format(
                entry.get('power'), name)
    return None


def certificate_issue(payload):
    """Why an artifact's continuation was not certified by this code, or None.

    ``CONTINUATION_CERTIFICATE`` names the version. Artifacts written before
    v2 never judged the last chain sample of a step, and the verdict they
    discarded cannot be restored from the file, so they must be re-scanned.
    """
    recorded = (payload.get('continuation_certificate') or {}).get('name')
    if recorded != CONTINUATION_CERTIFICATE:
        return 'artifact certificate {!r} is not {!r}'.format(
            recorded, CONTINUATION_CERTIFICATE)
    for branch in payload.get('branches', []):
        for key in CONTINUATION_PATHS:
            issue = _path_certificate_issue(branch.get(key) or [])
            if issue is not None:
                return 'branch {} {}: {}'.format(branch.get('id'), key, issue)
    return None


def _tolerance_mismatch(path, u_tolerance, log_ne_tolerance):
    """Name the first certificate tolerance in a path that differs from the artifact's."""
    declared = path[0].get('continuation_settings')
    records = [('settings', declared)] if declared is not None else []
    records.extend(('certificate at {!r} W'.format(entry.get('power')), entry['certificate'])
                   for entry in path if entry.get('certificate'))
    for label, record in records:
        for key, expected in (('u_tolerance', u_tolerance),
                              ('log_ne_tolerance', log_ne_tolerance)):
            if key in record and float(record[key]) != expected:
                return '{} {} {!r} differs from the artifact\'s {!r}'.format(
                    label, key, record[key], expected)
    return None


def record_continuation(payload):
    """Write the certificate, settings, coverage and consistency at top level.

    A path whose certificate was judged with tolerances other than the
    artifact's is refused: its certificate says nothing about the artifact's
    branch identity. So is a path certified by an older certificate: writing
    this code's metadata over it would not restore what it never checked.
    """
    u_tolerance = float(payload['u_tolerance'])
    log_ne_tolerance = float(payload['log_ne_tolerance'])
    settings = {}
    coverage = {}
    consistency = {}
    radius = {}
    solves = {}
    resolutions = []
    for branch in payload.get('branches', []):
        branch_id = branch.get('id')
        for key in CONTINUATION_PATHS:
            path = branch.get(key) or []
            if not path:
                continue
            mismatch = (_tolerance_mismatch(path, u_tolerance, log_ne_tolerance)
                        or _path_certificate_issue(path))
            if mismatch is not None:
                raise ValueError('branch {} {}: {}'.format(branch_id, key, mismatch))
            radius.setdefault(branch_id, {})[key] = aliasing_radius(
                path, u_tolerance, log_ne_tolerance)
            declared = path[0].get('continuation_settings')
            if declared is not None:
                settings.setdefault(branch_id, {})[key] = declared
                resolutions.append(float(declared['stability_resolution']))
            solves.setdefault(branch_id, {})[key] = {
                'corrector_solves': sum(int(entry.get('corrector_solves', 0))
                                        for entry in path),
                'certificate_solves': sum(int(entry.get('certificate_solves', 0))
                                          for entry in path),
            }
        if branch.get('continuation_path'):
            coverage[branch_id] = certified_coverage(branch)
            consistency[branch_id] = path_consistency(
                branch, u_tolerance, log_ne_tolerance)
    payload['continuation_certificate'] = {
        'name': CONTINUATION_CERTIFICATE,
        'statement': CERTIFICATE_STATEMENT,
        'limits': list(CERTIFICATE_LIMITS),
        'identity': {'u_tolerance': u_tolerance,
                     'log_ne_tolerance': log_ne_tolerance,
                     'components': ['u', 'ln n_e']},
        'branches_py_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    payload['continuation_settings'] = settings
    payload['stability_resolution'] = max(resolutions) if resolutions else None
    payload['certified_coverage'] = coverage
    payload['path_consistency'] = consistency
    payload['aliasing_radius'] = radius
    payload['continuation_solves'] = solves
    return payload


def _continuation_state(state):
    values = state.get('_state_vector')
    if values is None:
        raise AdapterContractError('adapter state has no full _state_vector')
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or not len(array) or not np.isfinite(array).all():
        raise ValueError('corrector returned an invalid full state')
    return array


def _component_scale(*vectors):
    values = np.vstack([np.abs(np.asarray(vector, dtype=float)) for vector in vectors])
    return np.maximum(np.max(values, axis=0), np.finfo(float).tiny)


def _scaled_norm(delta, scale):
    return float(np.linalg.norm(np.asarray(delta, dtype=float) / scale) /
                 np.sqrt(len(scale)))


def _corrector_metrics(current, predicted, corrected):
    """Return the RMS state-scaled corrector departure and tangent step."""
    scale = _component_scale(current, predicted, corrected)
    correction = _scaled_norm(corrected - predicted, scale)
    tangent_step = _scaled_norm(predicted - current, scale)
    return correction, tangent_step


def _component_noise(first, second):
    """Per-component scatter of two solves of one state, in its own scale."""
    return np.abs(second - first) / _component_scale(first, second)


def _departure(current, predicted, corrected, power_step, speeds,
               departure_fraction, noise):
    """Test every component's secant departure against its arclength step.

    Component ``i`` departs when ``|x_i - p_i| > departure_fraction *
    hypot(power_step * speeds_i, |p_i - x0_i|) + noise_i``, all in units of the
    component's own magnitude. Numerator and denominator share those units, so
    no rescaling or translation of any component changes the verdict, and one
    fast component cannot dilute a jump in another. ``speeds_i`` is the largest
    scaled speed component ``i`` has shown on this path, so its arclength step
    does not vanish where that component passes through an extremum. On one
    smooth branch the departure is O(h^2) against an O(h) step and falls under
    refinement; a jump keeps an O(1) departure and does not.
    """
    scale = _component_scale(current, predicted, corrected)
    departure = np.abs(corrected - predicted) / scale
    arclength = np.hypot(power_step * speeds, np.abs(predicted - current) / scale)
    ratio = float(np.max(departure / np.maximum(arclength, np.finfo(float).tiny)))
    departed = bool(np.any(departure > departure_fraction * arclength + noise))
    return departed, ratio


def _tangent_cosine(previous, current, corrected, previous_power,
                    current_power, corrected_power, speeds, noise):
    """Smallest per-component secant cosine in weighted (power, x_i) planes.

    A component whose step is within its corrector noise has no measurable
    direction and is skipped; ``None`` means no component had one.
    """
    scale = _component_scale(previous, current, corrected)
    old = (current - previous) / scale
    new = (corrected - current) / scale
    resolvable = (np.abs(old) > noise) & (np.abs(new) > noise)
    if not resolvable.any():
        return None
    old_power = (current_power - previous_power) * speeds[resolvable]
    new_power = (corrected_power - current_power) * speeds[resolvable]
    old, new = old[resolvable], new[resolvable]
    cosines = (old * new + old_power * new_power) / (
        np.hypot(old, old_power) * np.hypot(new, new_power))
    return float(np.min(cosines))


def _matched_eigenvalues(previous, current):
    left = [complex(value['real'], value['imag'])
            for value in previous['eigenvalues']]
    remaining = [complex(value['real'], value['imag'])
                 for value in current['eigenvalues']]
    if len(left) != len(remaining):
        return None
    pairs = []
    for value in sorted(left, key=lambda item: (abs(item.real), abs(item.imag))):
        index = min(range(len(remaining)), key=lambda item: abs(remaining[item] - value))
        pairs.append((value, remaining.pop(index)))
    return pairs


def _eigenvalue_motion(previous, current, fraction, band=EIGENVALUE_BAND):
    """Certify one sample gap from matched reduced eigenvalues.

    Any matched pair whose real part changes sign is a crossing. The motion
    bound ``|delta Re| <= fraction * min(|Re|)`` applies only to the modes
    nearest zero, those within ``band`` times the slowest mode's distance: a
    fast mode stiffening far from zero is no evidence of a stability change.
    """
    pairs = _matched_eigenvalues(previous, current)
    if pairs is None:
        return {'certified': False, 'ratio': float('inf'), 'bounded_modes': 0,
                'crossing': None, 'reason': 'reduced eigenspace dimension changed'}
    zero_tolerance = max(float(previous['zero_tolerance']),
                         float(current['zero_tolerance']))
    distances = [max(min(abs(left.real), abs(right.real)), zero_tolerance)
                 for left, right in pairs]
    nearest = min(distances, default=zero_tolerance)
    ratios = []
    crossing = None
    for (left, right), distance in zip(pairs, distances):
        if distance <= band * nearest:
            ratios.append(abs(right.real - left.real) / distance)
        if left.real * right.real <= 0.0 and left.real != right.real:
            complex_pair = max(abs(left.imag), abs(right.imag)) > zero_tolerance
            crossing = 'hopf' if complex_pair else (crossing or 'stability')
    ratio = max(ratios, default=0.0)
    return {
        'certified': crossing is None and ratio <= fraction,
        'ratio': float(ratio),
        'bounded_modes': len(ratios),
        'crossing': crossing,
        'reason': ('reduced eigenvalue crossed zero' if crossing
                   else 'reduced eigenvalues moved too far to certify the interval'),
    }


def _sample_count(interval, resolution):
    """Even number of gaps no wider than ``resolution`` (to SAMPLE_ROUNDOFF)."""
    cells = interval / (2.0 * resolution)
    return max(2, 2 * int(np.ceil(cells * (1.0 - SAMPLE_ROUNDOFF))))


def _nearest_real_distance(stability):
    return min(abs(float(value['real'])) for value in stability['eigenvalues'])


def _event_entry(kind, reason, bracket, power, predictor, last_power,
                 correction, tangent_cosine=None):
    lower, upper = sorted(float(value) for value in bracket)
    boundary = max((lower, upper), key=lambda value: abs(value - float(last_power)))
    fold = kind == 'fold'
    return {
        'power': float(boundary),
        'seed': asdict(predictor),
        'converged': False,
        'failure': reason,
        'event_type': kind,
        'event_bracket': [lower, upper],
        'fold_detected': fold,
        'fold_reason': reason if fold else None,
        'fold_bracket': [lower, upper] if fold else None,
        'hopf_detected': kind == 'hopf',
        'unresolved': kind == 'unresolved',
        'last_accepted_power': float(last_power),
        'predictor_relative_distance': correction,
        'corrector_scaled_departure': correction,
        'tangent_cosine': tangent_cosine,
        'segment_certified': False,
    }


def continue_branch(adapter, seed, powers, timeout=None, *, initial_power=None,
                    max_step=None, min_step=None, predictor_tolerance=0.25,
                    departure_fraction=2.0,
                    eigenvalue_fraction=0.5, eigenvalue_band=EIGENVALUE_BAND,
                    stability_resolution=None,
                    u_tolerance=U_TOLERANCE,
                    log_ne_tolerance=LOG_NE_TOLERANCE, lease_s=RUN_LEASE):
    """Follow one branch with certified adaptive predictor-corrector steps.

    A step is accepted only with the reversible-chain certificate
    (``CERTIFICATE_STATEMENT``). The step is split into gaps no wider than
    ``stability_resolution``. Each sample is corrected from the previous one,
    and each sample corrected back to the previous power returns to the
    previous sample. The last sample matches the predictor-corrected endpoint.
    All of these hold within (``u_tolerance``, ``log_ne_tolerance``), the
    tolerances that define a distinct branch. The first step of a leg is
    certified the same way.

    The secant tests are extra refinement triggers, not the certificate. The
    corrector departure is scaled component by component and must satisfy
    ``departure <= departure_fraction * arclength_step + noise``. Each
    component's arclength weights power by its speed over the last two
    accepted steps. ``noise`` is the corrector scatter on re-solving the
    initial state, which must itself reproduce within the branch tolerances.
    A secant that reverses orientation is refined, not called a fold.

    Every sample is classified, and a sample that is not ``stable`` ends the
    path with an event. On every gap, the matched reduced eigenvalues near
    zero must satisfy
    ``|delta Re(lambda)| <= eigenvalue_fraction * min(|Re(lambda)|)``, so an
    unstable interval at least as wide as the resolution cannot be skipped.
    A failed certificate is refined; at the minimum step it is ``unresolved``
    with its bracket. A fold needs fold-test or corrector-failure evidence,
    and a complex crossing is Hopf. If a corrector fails while an event is
    being bracketed, the bracket is ``unresolved``, unless the adapter
    declares the failure a missing steady state.

    ``lease_s`` has the same finite, positive ownership-lease meaning as in
    ``scan()`` when a timeout requests isolated corrector workers. Lifecycle
    cancellation propagates; it is not classified as a corrector failure.
    """
    if not math.isfinite(lease_s) or lease_s <= 0.0:
        raise ValueError('lease_s must be finite and positive')
    targets = [float(power) for power in powers]
    if not targets:
        return []
    if not (np.isfinite(targets).all()
            and np.isfinite(predictor_tolerance) and predictor_tolerance > 0.0
            and np.isfinite(departure_fraction) and departure_fraction > 0.0
            and np.isfinite(eigenvalue_fraction) and 0.0 < eigenvalue_fraction < 1.0
            and np.isfinite(eigenvalue_band) and eigenvalue_band >= 1.0):
        raise ValueError('continuation powers and certification fractions must be finite')
    current_power = float(targets[0] if initial_power is None else initial_power)
    if not np.isfinite(current_power):
        raise ValueError('initial continuation power must be finite')
    if 'state' not in seed:
        raise AdapterContractError('continuation seed carries no full state')
    initial_vector = np.asarray(seed['state'], dtype=float)
    if initial_vector.ndim != 1 or not len(initial_vector) or not np.isfinite(initial_vector).all():
        raise ValueError('initial continuation state must be a finite vector')
    current = {'u': float(seed['u']), 'n_e': float(seed['n_e']),
               'vector': initial_vector}
    if not (np.isfinite(current['u']) and np.isfinite(current['n_e'])
            and current['n_e'] > 0.0):
        raise ValueError('initial continuation u and n_e must be finite, n_e positive')
    span = max([abs(target - current_power) for target in targets] or [0.0])
    step_ceiling = float(
        max_step if max_step is not None else max(span / 20.0, 1.e-6))
    step_floor = float(
        min_step if min_step is not None else max(step_ceiling / 64.0, 1.e-12))
    stability_resolution = float(
        step_floor if stability_resolution is None else stability_resolution)
    u_tolerance = float(u_tolerance)
    log_ne_tolerance = float(log_ne_tolerance)
    if not (np.isfinite(step_ceiling) and np.isfinite(step_floor)
            and np.isfinite(stability_resolution) and step_ceiling > 0.0
            and 0.0 < step_floor <= stability_resolution <= step_ceiling
            and np.isfinite(u_tolerance) and u_tolerance > 0.0
            and np.isfinite(log_ne_tolerance) and log_ne_tolerance > 0.0):
        raise ValueError(
            'continuation bounds, resolution and branch tolerances must '
            'be finite, positive and ordered')
    settings = {
        'certificate': CONTINUATION_CERTIFICATE,
        'targets': targets,
        'initial_power': current_power,
        'max_step': step_ceiling,
        'min_step': step_floor,
        'stability_resolution': stability_resolution,
        'sample_roundoff': SAMPLE_ROUNDOFF,
        'predictor_tolerance': float(predictor_tolerance),
        'departure_fraction': float(departure_fraction),
        'eigenvalue_fraction': float(eigenvalue_fraction),
        'eigenvalue_band': float(eigenvalue_band),
        'u_tolerance': u_tolerance,
        'log_ne_tolerance': log_ne_tolerance,
    }
    fold_failure = getattr(adapter, 'is_fold_failure', lambda _failure: False)
    solves = {'corrector': 0}

    def same_branch(left, right):
        return same_cluster(left, right, u_tolerance, log_ne_tolerance)

    def correct(power, predicted_u, predicted_ne, predicted_vector):
        predictor = Seed(float(predicted_u), float(predicted_ne), tuple(predicted_vector))
        solves['corrector'] += 1
        outcome = _integrate_many(adapter, [predictor], 1, timeout, power=power, lease_s=lease_s)[0]
        if outcome[0] != 'ok':
            return predictor, None, None, None, outcome[1]
        try:
            state = dict(outcome[1])
            state['u'], state['n_e'] = _state_values(state)
            corrected = _continuation_state(state)
            if corrected.shape != np.asarray(predicted_vector).shape:
                raise ValueError(
                    'corrector state shape {} differs from predictor {}'.format(
                        corrected.shape, np.asarray(predicted_vector).shape))
            stability = _analysis(adapter, state)
        except AdapterContractError:
            raise
        except Exception as exc:
            return predictor, None, None, None, '{}: {}'.format(type(exc).__name__, exc)
        return predictor, state, corrected, stability, None

    def secant_ne(current_ne, old_ne, ratio):
        # Extrapolate ln n_e, so the predictor density is always positive.
        return float(np.exp(np.log(current_ne) + ratio * (
            np.log(current_ne) - np.log(old_ne))))

    def refine_eigen_event(left, right, kind, reason):
        """Bisect a measured real/complex crossing without accepting probes."""
        while abs(right[0] - left[0]) > step_floor:
            power = (left[0] + right[0]) / 2.0
            vector = (left[2] + right[2]) / 2.0
            u = (left[1]['u'] + right[1]['u']) / 2.0
            n_e = float(np.sqrt(left[1]['n_e'] * right[1]['n_e']))
            _predictor, state, corrected, stability, mid_failure = correct(
                power, u, n_e, vector)
            if mid_failure is not None:
                return ('unresolved',
                        'corrector failed while bracketing a {} event: {}'.format(
                            kind, mid_failure),
                        (left[0], right[0]))
            middle = (power, state, corrected, stability)
            if _eigenvalue_motion(left[3], stability, eigenvalue_fraction,
                                  eigenvalue_band)['crossing']:
                right = middle
            else:
                left = middle
        return (kind, reason, (left[0], right[0]))

    def refine_predictor_fold(previous_power, previous_state, left_power,
                              left_state, right_power, kind, reason):
        """Bisect where the natural corrector can no longer follow its secant.

        Only a failure the adapter declares a missing steady state moves the
        upper end. Any other failure, such as a timeout, ends ``unresolved``.
        """
        while abs(right_power - left_power) > step_floor:
            power = (left_power + right_power) / 2.0
            ratio = (power - left_power) / (left_power - previous_power)
            vector = left_state['vector'] + ratio * (
                left_state['vector'] - previous_state['vector'])
            u = left_state['u'] + ratio * (left_state['u'] - previous_state['u'])
            n_e = secant_ne(left_state['n_e'], previous_state['n_e'], ratio)
            _predictor, state, corrected, _stability, mid_failure = correct(
                power, u, n_e, vector)
            if mid_failure is not None:
                if fold_failure(mid_failure):
                    right_power = power
                    continue
                return ('unresolved',
                        'corrector failed while bracketing a fold: ' + mid_failure,
                        (left_power, right_power))
            middle = {'u': state['u'], 'n_e': state['n_e'], 'vector': corrected}
            correction, _tangent_step = _corrector_metrics(
                left_state['vector'], vector, corrected)
            departed, _ratio = _departure(
                left_state['vector'], vector, corrected,
                abs(power - left_power), speeds, departure_fraction, noise)
            if departed or correction > predictor_tolerance:
                right_power = power
            else:
                previous_power, previous_state = left_power, left_state
                left_power, left_state = power, middle
        return (kind, reason, (left_power, right_power))

    path = []
    initial_predictor, initial_result, initial_corrected, initial_stability, failure = correct(
        current_power, current['u'], current['n_e'], current['vector'])
    initial_correction = float('inf')
    reproducibility = float('inf')
    noise = None
    reproduced = False
    if failure is None:
        initial_correction, _ = _corrector_metrics(
            current['vector'], current['vector'], initial_corrected)
        # Re-solve the corrected state: its scatter is the corrector's own
        # convergence noise, the floor below which a departure is no evidence.
        _repeat_predictor, repeat, repeat_corrected, _repeat_stability, failure = correct(
            current_power, initial_result['u'], initial_result['n_e'],
            initial_corrected)
        if failure is None:
            reproducibility, _ = _corrector_metrics(
                initial_corrected, initial_corrected, repeat_corrected)
            noise = np.maximum(
                64.0 * np.finfo(float).eps,
                4.0 * _component_noise(initial_corrected, repeat_corrected))
            reproduced = same_branch(initial_result, repeat)
    seed_branch = failure is None and same_branch(current, initial_result)
    if (failure is not None or initial_correction > predictor_tolerance
            or reproducibility > predictor_tolerance or not reproduced
            or not seed_branch):
        if failure is not None:
            reason = 'corrector failed: ' + failure
        elif not reproduced:
            reason = ('corrector does not reproduce the initial state within '
                      'the branch tolerances')
        elif not seed_branch:
            reason = ('initial state (u={!r}, n_e={!r}) is not the seed\'s branch '
                      '(u={!r}, n_e={!r}) within the branch tolerances'.format(
                          initial_result['u'], initial_result['n_e'],
                          current['u'], current['n_e']))
        else:
            reason = 'corrector left scaled predictor neighbourhood'
        path.append(_event_entry(
            'unresolved', reason, (current_power, current_power), current_power,
            initial_predictor, current_power,
            None if not np.isfinite(initial_correction) else initial_correction))
        return path
    if initial_stability['near_zero_reduced_mode']:
        path.append(_event_entry(
            'unresolved', 'initial reduced Jacobian is singular',
            (current_power, current_power), current_power, initial_predictor,
            current_power, initial_correction))
        return path
    if not initial_stability['stable']:
        path.append(_event_entry(
            'unresolved',
            'initial state is {}; a non-stable state is not continued'.format(
                initial_stability['classification']),
            (current_power, current_power), current_power, initial_predictor,
            current_power, initial_correction))
        return path
    path.append({
        'power': current_power,
        'seed': asdict(initial_predictor),
        'converged': True,
        'u': initial_result['u'],
        'n_e': initial_result['n_e'],
        'state': initial_corrected.tolist(),
        'determinant': initial_stability['reduced_determinant'],
        'stability': initial_stability,
        'discharge_class': initial_result.get('discharge_class'),
        'event_type': None,
        'fold_detected': False,
        'fold_reason': None,
        'hopf_detected': False,
        'unresolved': False,
        'predictor_relative_distance': initial_correction,
        'corrector_scaled_departure': initial_correction,
        'corrector_reproducibility': reproducibility,
        'corrector_noise_floor': noise.tolist(),
        'departure_fraction': departure_fraction,
        'tangent_cosine': None,
        'step_size': 0.0,
        'segment_certified': True,
        'stability_resolution': stability_resolution,
        'continuation_settings': settings,
        'corrector_solves': solves['corrector'],
        'certificate_solves': 0,
    })
    solves['corrector'] = 0
    current = {'u': initial_result['u'], 'n_e': initial_result['n_e'],
               'vector': initial_corrected}
    history = [(current_power, current.copy())]
    previous_stability = initial_stability
    # Each component's power axis is weighted by its scaled speed over the
    # last two accepted steps. The arclength step then does not vanish at an
    # extremum, and an earlier fast stretch cannot loosen a later step.
    recent_speeds = []
    speeds = np.zeros_like(initial_corrected)
    previous_direction = 0.0

    for target in targets:
        direction = float(np.sign(target - current_power))
        if direction == 0.0:
            continue
        if previous_direction and direction != previous_direction:
            history = [(current_power, current.copy())]
        previous_direction = direction
        step = step_ceiling
        while direction * (target - current_power) > 0.0:
            remaining = abs(target - current_power)
            # Never leave a remainder below the minimum step: a step shorter
            # than the declared floor (roundoff, or the tail after refinement)
            # resolves nothing and its corrector noise can mimic a fold. A
            # remainder of exactly step + floor is split, so a refined step
            # is taken rather than snapped back to the target. A short tail
            # longer than the ceiling is taken in two equal steps, so no step
            # exceeds max_step; each half is at least the floor whenever
            # max_step is at least twice min_step.
            if remaining >= step + step_floor * (1.0 - SAMPLE_ROUNDOFF):
                trial_power = current_power + direction * step
            elif remaining <= step_ceiling * (1.0 + SAMPLE_ROUNDOFF):
                trial_power = target
            else:
                trial_power = current_power + direction * remaining / 2.0
            trial_interval = abs(trial_power - current_power)
            has_tangent = len(history) >= 2 and history[-1][0] != history[-2][0]
            if has_tangent:
                old_power, old = history[-2]
                ratio = (trial_power - current_power) / (current_power - old_power)
                predicted_vector = current['vector'] + ratio * (current['vector'] - old['vector'])
                predicted_u = current['u'] + ratio * (current['u'] - old['u'])
                predicted_ne = secant_ne(current['n_e'], old['n_e'], ratio)
            else:
                predicted_vector = current['vector'].copy()
                predicted_u = current['u']
                predicted_ne = current['n_e']
            predictor, state, corrected, stability, failure = correct(
                trial_power, predicted_u, predicted_ne, predicted_vector)
            correction = float('inf')
            correction_ratio = None
            departed = False
            tangent_cosine = None
            issue = None
            issue_samples = None
            issue_source = None
            if failure is None:
                correction, _tangent_step = _corrector_metrics(
                    current['vector'], predicted_vector, corrected)
                if has_tangent:
                    departed, correction_ratio = _departure(
                        current['vector'], predicted_vector, corrected,
                        trial_interval, speeds, departure_fraction, noise)
                    tangent_cosine = _tangent_cosine(
                        history[-2][1]['vector'], current['vector'], corrected,
                        history[-2][0], current_power, trial_power, speeds,
                        noise)
                previous_test = (
                    _nearest_real_distance(path[-2]['stability'])
                    if len(path) >= 2 and path[-2].get('converged') else None)
                current_test = _nearest_real_distance(previous_stability)
                trial_test = _nearest_real_distance(stability)
                fold_test_rebounded = (
                    previous_test is not None
                    and current_test < previous_test
                    and trial_test > current_test)
                beyond_failures = []
                if has_tangent and stability['near_zero_reduced_mode']:
                    old_power, old = history[-2]
                    for distance in (step_floor, 2.0 * step_floor):
                        beyond_power = trial_power + direction * distance
                        beyond_ratio = ((beyond_power - current_power)
                                        / (current_power - old_power))
                        beyond_vector = current['vector'] + beyond_ratio * (
                            current['vector'] - old['vector'])
                        beyond_u = current['u'] + beyond_ratio * (
                            current['u'] - old['u'])
                        beyond_ne = secant_ne(current['n_e'], old['n_e'], beyond_ratio)
                        (_beyond_predictor, _beyond_state,
                         _beyond_corrected, _beyond_stability,
                         beyond_failure) = correct(
                            beyond_power, beyond_u, beyond_ne, beyond_vector)
                        beyond_failures.append(beyond_failure)
                repeated_beyond_failure = beyond_failures and all(beyond_failures)
                confirmed_fold_failure = repeated_beyond_failure and all(
                    fold_failure(item) for item in beyond_failures)
                if repeated_beyond_failure and not confirmed_fold_failure:
                    issue = (
                        'unresolved',
                        'corrector failure while probing beyond singular endpoint',
                        (current_power, trial_power),
                    )
                elif confirmed_fold_failure:
                    issue = (
                        'fold',
                        'corrector failed at two powers beyond singular reduced Jacobian',
                        (trial_power, trial_power + direction * step_floor),
                    )
                    issue_source = 'eigen'
                elif departed and fold_test_rebounded:
                    issue = ('fold',
                             'fold test function rebounded after tangent predictor departure',
                             (current_power, trial_power))
                    issue_source = 'predictor'
                elif tangent_cosine is not None and tangent_cosine <= 0.0:
                    issue = ('unresolved', 'continuation tangent reversed orientation',
                             (current_power, trial_power))
                elif has_tangent and correction > predictor_tolerance:
                    issue = ('unresolved', 'corrector left scaled predictor neighbourhood',
                             (current_power, trial_power))
                elif departed:
                    issue = ('unresolved',
                             'corrector departed from the secant tangent by more '
                             'than the declared fraction',
                             (current_power, trial_power))
            else:
                issue = ('unresolved', 'corrector failed: ' + failure,
                         (current_power, trial_power))

            certification_points = []
            interpolation_checks = []
            interpolation_diagnostics = []
            motions = []
            certificate = None
            sample_count = 0
            if issue is None:
                # The reversible-chain certificate: zero-order samples across
                # the step, each corrected back to the sample before it.
                sample_count = _sample_count(trial_interval, stability_resolution)
                chain = [(
                    current_power,
                    {'u': current['u'], 'n_e': current['n_e']},
                    current['vector'],
                    previous_stability,
                )]
                reverse_gaps = []
                for sample_number in range(1, sample_count + 1):
                    sample_power = (
                        trial_power if sample_number == sample_count
                        else current_power + float(sample_number) / sample_count * (
                            trial_power - current_power))
                    left = chain[-1]
                    (_sample_predictor, sample_state, sample_corrected,
                     sample_stability, sample_failure) = correct(
                        sample_power, left[1]['u'], left[1]['n_e'], left[2])
                    if sample_failure is not None:
                        issue = (
                            'unresolved',
                            'certificate corrector failed: ' + sample_failure,
                            (left[0], sample_power),
                        )
                        break
                    (_back_predictor, back_state, _back_corrected,
                     _back_stability, back_failure) = correct(
                        left[0], sample_state['u'], sample_state['n_e'],
                        sample_corrected)
                    if back_failure is not None:
                        issue = (
                            'unresolved',
                            'certificate reverse corrector failed: ' + back_failure,
                            (left[0], sample_power),
                        )
                        break
                    reverse_gaps.append(_identity_gap(back_state, left[1]))
                    if not same_branch(back_state, left[1]):
                        issue = (
                            'unresolved',
                            'certificate: correcting back from {!r} W does not '
                            'return to the branch at {!r} W'.format(
                                sample_power, left[0]),
                            (left[0], sample_power),
                        )
                        break
                    chain.append((sample_power, sample_state, sample_corrected,
                                  sample_stability))
                if issue is None:
                    forward_gap = _identity_gap(chain[-1][1], state)
                    if not same_branch(chain[-1][1], state):
                        issue = (
                            'unresolved',
                            'certificate: the zero-order chain and the predictor '
                            'corrector reach different branches at {!r} W'.format(
                                trial_power),
                            (chain[-2][0], trial_power),
                        )
                if issue is None:
                    certificate = {
                        'name': CONTINUATION_CERTIFICATE,
                        'sub_steps': sample_count,
                        'gap': float(trial_interval / sample_count),
                        'max_reverse_u_difference': max(gap[0] for gap in reverse_gaps),
                        'max_reverse_log_ne_difference': max(
                            gap[1] for gap in reverse_gaps),
                        'forward_u_difference': forward_gap[0],
                        'forward_log_ne_difference': forward_gap[1],
                        'u_tolerance': u_tolerance,
                        'log_ne_tolerance': log_ne_tolerance,
                    }
                    certificate['sample_verdicts'] = [
                        sample[3]['classification'] for sample in chain[1:]]
                    for sample in chain[1:-1]:
                        certification_points.append({
                            'power': float(sample[0]),
                            'u': float(sample[1]['u']),
                            'n_e': float(sample[1]['n_e']),
                            'state': sample[2].tolist(),
                            'stability': sample[3]['classification'],
                            'converged': True,
                            'segment_certified': True,
                        })
                    samples = chain[:-1] + [(trial_power, state, corrected, stability)]
                    uncertified = None
                    for left, right in zip(samples, samples[1:]):
                        motion = _eigenvalue_motion(
                            left[3], right[3], eigenvalue_fraction, eigenvalue_band)
                        motions.append(motion)
                        bracket = (left[0], right[0])
                        if motion['crossing'] is not None:
                            issue = (
                                motion['crossing'], motion['reason'], bracket)
                            issue_samples = (left, right)
                            issue_source = 'eigen'
                            break
                        if not right[3]['stable']:
                            verdict = right[3]['classification']
                            issue = (
                                'stability' if verdict == 'unstable' else 'unresolved',
                                'stability verdict {!r} at {!r} W; a non-stable '
                                'state is not continued'.format(verdict, right[0]),
                                bracket)
                            issue_source = 'verdict'
                            break
                        if not motion['certified'] and uncertified is None:
                            uncertified = (
                                'unresolved', motion['reason'], bracket)
                    if issue is None and not chain[-1][3]['stable']:
                        # The endpoint replaces the last sample in the
                        # recorded path; the sample's own verdict still counts.
                        verdict = chain[-1][3]['classification']
                        issue = (
                            'stability' if verdict == 'unstable' else 'unresolved',
                            'stability verdict {!r} of the last chain sample at '
                            '{!r} W; a non-stable state is not continued'.format(
                                verdict, trial_power),
                            (chain[-2][0], trial_power))
                        issue_source = 'verdict'
                    if issue is None:
                        issue = uncertified
                    if issue is None:
                        for left, middle, right in zip(
                                samples[::2], samples[1::2], samples[2::2]):
                            fraction = ((middle[0] - left[0])
                                        / (right[0] - left[0]))
                            estimated_u = float(left[1]['u']) + fraction * (
                                float(right[1]['u']) - float(left[1]['u']))
                            estimated_log_ne = np.log(float(left[1]['n_e'])) + fraction * (
                                np.log(float(right[1]['n_e']))
                                - np.log(float(left[1]['n_e'])))
                            u_error = abs(estimated_u - float(middle[1]['u']))
                            log_ne_error = abs(
                                estimated_log_ne - np.log(float(middle[1]['n_e'])))
                            interpolation_checks.append({
                                'power': [float(left[0]), float(right[0])],
                                'midpoint_u_error': float(u_error),
                                'midpoint_log_ne_error': float(log_ne_error),
                                'resolution': float(max(
                                    abs(middle[0] - left[0]),
                                    abs(right[0] - middle[0]))),
                            })
                            estimator = getattr(
                                adapter, 'estimate_interpolation_error', None)
                            if estimator is not None:
                                estimate = estimator(
                                    {'power': left[0], 'u': left[1]['u'],
                                     'n_e': left[1]['n_e']},
                                    {'power': middle[0], 'u': middle[1]['u'],
                                     'n_e': middle[1]['n_e']},
                                    {'power': right[0], 'u': right[1]['u'],
                                     'n_e': right[1]['n_e']})
                                interpolation_diagnostics.append({
                                    'certified': False,
                                    'kind': 'adapter-estimated-midpoint-error',
                                    'power': [float(left[0]), float(right[0])],
                                    'u': [float(left[1]['u']),
                                          float(right[1]['u'])],
                                    'n_e': [float(left[1]['n_e']),
                                            float(right[1]['n_e'])],
                                    'u_error_estimate': float(
                                        estimate['u_error_estimate']),
                                    'log_ne_error_estimate': float(
                                        estimate['log_ne_error_estimate']),
                                    'resolution': float(max(
                                        abs(middle[0] - left[0]),
                                        abs(right[0] - middle[0]))),
                                })

            can_refine = (trial_interval > step_floor * (1.0 + 1.0e-12)
                          and step > step_floor * (1.0 + 1.0e-12))
            if issue is not None and (issue[0] in ('fold', 'hopf', 'stability')
                                      or issue_source == 'verdict'):
                kind, reason, bracket = issue
                if issue_source == 'predictor':
                    kind, reason, bracket = refine_predictor_fold(
                        history[-2][0], history[-2][1], current_power, current,
                        trial_power, kind, reason)
                elif issue_samples is not None:
                    kind, reason, bracket = refine_eigen_event(
                        *issue_samples, kind, reason)
                path.append(_event_entry(
                    kind, reason, bracket, trial_power, predictor, current_power,
                    None if not np.isfinite(correction) else correction,
                    tangent_cosine=tangent_cosine))
                return path
            if issue is not None and can_refine:
                step = max(trial_interval / 2.0, step_floor)
                continue
            if issue is not None:
                kind, reason, bracket = issue
                path.append(_event_entry(
                    kind, reason, bracket, trial_power, predictor, current_power,
                    None if not np.isfinite(correction) else correction,
                    tangent_cosine=tangent_cosine))
                return path

            maximum_eigen_ratio = max(
                (motion['ratio'] for motion in motions), default=0.0)
            entry = {
                'power': float(trial_power),
                'seed': asdict(predictor),
                'converged': True,
                'u': state['u'],
                'n_e': state['n_e'],
                'state': corrected.tolist(),
                'determinant': stability['reduced_determinant'],
                'stability': stability,
                'discharge_class': state.get('discharge_class'),
                'event_type': None,
                'fold_detected': False,
                'fold_reason': None,
                'hopf_detected': False,
                'unresolved': False,
                'predictor_relative_distance': correction,
                'corrector_scaled_departure': correction,
                'corrector_to_tangent_ratio': correction_ratio,
                'tangent_cosine': tangent_cosine,
                'eigenvalue_change_ratio': maximum_eigen_ratio,
                'eigenvalue_fraction': eigenvalue_fraction,
                'stability_probe_power': (
                    certification_points[len(certification_points) // 2]['power']
                    if certification_points else None),
                'stability_resolution': float(trial_interval / sample_count),
                'certification_points': certification_points,
                'interpolation_checks': interpolation_checks,
                'interpolation_diagnostics': interpolation_diagnostics,
                'certificate': certificate,
                'corrector_solves': solves['corrector'],
                'certificate_solves': sample_count + 1,
                'step_size': float(trial_interval),
                'segment_certified': True,
            }
            path.append(entry)
            solves['corrector'] = 0
            recent_speeds = (recent_speeds + [
                np.abs(corrected - current['vector']) / (
                    _component_scale(current['vector'], corrected) * trial_interval)])[-2:]
            speeds = np.max(np.vstack(recent_speeds), axis=0)
            current_power = trial_power
            current = {'u': state['u'], 'n_e': state['n_e'], 'vector': corrected}
            history.append((current_power, current.copy()))
            history = history[-2:]
            previous_stability = stability
            if correction < predictor_tolerance / 4.0 and maximum_eigen_ratio < eigenvalue_fraction / 4.0:
                step = min(step * 1.5, step_ceiling)
            direction = float(np.sign(target - current_power))
    return path


def _load_artifact(path):
    """Read ``branches.json``, refusing one not certified by this code."""
    try:
        payload = json.loads(Path(path).read_text())
    except (OSError, ValueError) as exc:
        raise ValueError('cannot read operating branches {}: {}'.format(path, exc)) from exc
    if not isinstance(payload, dict):
        raise ValueError('operating branches {} is not a branch artifact'.format(path))
    issue = certificate_issue(payload)
    if issue is not None:
        raise ValueError(
            'operating branches {} predates the current continuation certificate: '
            '{}; re-run the scan'.format(path, issue))
    issue = tolerance_issue(payload)
    if issue is not None:
        raise ValueError(
            'operating branches {} matches at tolerances its certificates were not judged '
            'with: {}; restore the scan\'s tolerances or re-run the scan'.format(path, issue))
    issue = reference_issue(payload)
    if issue is not None:
        raise ValueError('operating branches {} {}'.format(path, issue))
    return payload


def _artifact_tolerances(payload):
    try:
        u_tolerance = float(payload['u_tolerance'])
        log_ne_tolerance = float(payload['log_ne_tolerance'])
    except (KeyError, TypeError, ValueError, OverflowError) as exc:
        raise ValueError(
            'operating branches artifact has no valid recorded tolerances') from exc
    if not (np.isfinite(u_tolerance) and u_tolerance > 0.0
            and np.isfinite(log_ne_tolerance) and log_ne_tolerance > 0.0):
        raise ValueError('operating branches artifact tolerances must be finite and positive')
    return u_tolerance, log_ne_tolerance


def tolerance_issue(payload):
    """Name the first recorded tolerance that differs from the matching tolerances, or None.

    Matching reads the top-level ``u_tolerance`` and ``log_ne_tolerance``, but
    the certificate identity, each path's settings and each step's certificate
    record the tolerances the scan judged. Edited settings are refused, not
    applied: the certificates say nothing about branch identity at other
    tolerances.
    """
    u_tolerance, log_ne_tolerance = _artifact_tolerances(payload)
    identity = (payload.get('continuation_certificate') or {}).get('identity') or {}
    records = [('certificate identity', identity)]
    records.extend(('branch {} {} settings'.format(branch_id, key), settings)
                   for branch_id, paths in (payload.get('continuation_settings') or {}).items()
                   for key, settings in paths.items())
    for label, record in records:
        for key, expected in (('u_tolerance', u_tolerance),
                              ('log_ne_tolerance', log_ne_tolerance)):
            if key in record and float(record[key]) != expected:
                return '{} {} {!r} differs from the artifact\'s {!r}'.format(
                    label, key, record[key], expected)
    for branch in payload.get('branches', []):
        for key in CONTINUATION_PATHS:
            path = branch.get(key) or []
            issue = _tolerance_mismatch(path, u_tolerance, log_ne_tolerance) if path else None
            if issue is not None:
                return 'branch {} {}: {}'.format(branch.get('id'), key, issue)
    return None


def reference_issue(payload):
    """Why the recorded reference power is unusable or disagrees with the paths, or None.

    Every seed was solved at the reference power and every continuation path
    starts from it, so each path's settings and initial state record it.
    """
    recorded = payload.get('reference_power')
    try:
        reference_power = float(recorded)
    except (TypeError, ValueError, OverflowError):
        return 'records no valid reference power: {!r}'.format(recorded)
    if isinstance(recorded, bool) or not np.isfinite(reference_power):
        return 'records no valid reference power: {!r}'.format(recorded)
    records = [('branch {} {} settings'.format(branch_id, key), settings.get('initial_power'))
               for branch_id, paths in (payload.get('continuation_settings') or {}).items()
               for key, settings in paths.items()]
    for branch in payload.get('branches', []):
        for key in CONTINUATION_PATHS:
            path = branch.get(key) or []
            if path:
                start = path[0]
                records.append(('branch {} {} initial state'.format(branch.get('id'), key),
                                start.get('power')))
                records.append(('branch {} {} initial settings'.format(branch.get('id'), key),
                                (start.get('continuation_settings') or {}).get('initial_power')))
    for label, power in records:
        if power is not None and float(power) != reference_power:
            return 'records reference power {!r} W, but its {} starts at {!r} W'.format(
                recorded, label, power)
    return None


def matchable_issue(branch, u_tolerance, log_ne_tolerance):
    """Why a recorded branch can never be matched, or None.

    The scan records every converged cluster, but a run is assigned only to a
    branch whose reference verdict is stable and whose step and half-step
    paths certify consistent states. The reference power gets no exception:
    there, too, a candidate is a state a path solved.
    """
    stability = branch.get('stability') or {}
    if stability.get('stable') is not True:
        return 'its reference verdict is {!r}, so it is recorded for reporting only'.format(
            stability.get('classification'))
    if not _branch_certified_states(branch):
        return 'it has no certified continuation path, so it is recorded for reporting only'
    consistency = path_consistency(branch, u_tolerance, log_ne_tolerance)
    if not consistency['half_step_recorded']:
        return 'no half-step path is recorded, so nothing confirms its step path'
    if not consistency['consistent']:
        return ('its recorded states disagree at {} W and the half-step path does '
                'not reach {} W'.format(
                    ['{:.9g}'.format(value) for value in consistency['disagreeing_powers']],
                    ['{:.9g}'.format(value)
                     for value in consistency['unconfirmed_by_half_step']]))
    return None


def declared_branch(path, branch_id):
    """Read one matchable reactor operating branch by id from ``branches.json``."""
    payload = _load_artifact(path)
    for branch in payload.get('branches', []):
        if branch.get('id') == branch_id:
            issue = matchable_issue(branch, *_artifact_tolerances(payload))
            if issue is not None:
                raise ValueError('operating branch {!r} cannot be matched: {}'.format(
                    branch_id, issue))
            return branch
    raise ValueError('operating branch {!r} is absent from {}'.format(branch_id, path))


def _branch_states_at_power(branch, power):
    """Return the certified states either path of a matchable branch solved at ``power``.

    Only solved states count, at the reference power as anywhere else:
    nothing between two corrector solves is certified.
    """
    power_tolerance = 16.0 * np.finfo(float).eps * max(abs(float(power)), 1.0)
    return [{'u': float(entry['u']), 'n_e': float(entry['n_e'])}
            for entry in _branch_certified_states(branch)
            if abs(float(entry['power']) - power) <= power_tolerance]


def _uncertified_power_message(branch, power, issue):
    if issue is not None:
        return 'branch {} cannot be matched: {}'.format(branch.get('id'), issue)
    return 'power not certified for branch {} at {!r} W; certified powers: {}'.format(
        branch.get('id'), power,
        ['{:.9g}'.format(value) for value in certified_coverage(branch)['powers']])


def matching_branch(path, state, power=None):
    """Return the unique branch matching state at the run's absorbed power.

    When the state matches no certified state, and some branch has no
    certified state at ``power``, the refusal names that branch, its
    certified powers, rather than reporting no match. Only powers a path
    solved are certified; nothing is interpolated between them. A branch
    that ``matchable_issue`` refuses has no candidate at any power, and an
    artifact certified by an older certificate is refused whole.
    """
    payload = _load_artifact(path)
    u_tolerance, log_ne_tolerance = _artifact_tolerances(payload)
    if power is not None:
        try:
            power = float(power)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError('terminal matching power must be finite') from exc
        if not np.isfinite(power):
            raise ValueError('terminal matching power must be finite')
    matches = []
    uncovered = []
    for branch in payload.get('branches', []):
        issue = matchable_issue(branch, u_tolerance, log_ne_tolerance)
        if issue is not None:
            candidates = []
        elif power is None:
            candidates = [branch]
        else:
            candidates = _branch_states_at_power(branch, power)
        if not candidates:
            uncovered.append((branch, issue))
        if any(same_cluster(candidate, state, u_tolerance, log_ne_tolerance)
               for candidate in candidates):
            matches.append(branch.get('id'))
    if not matches and uncovered:
        raise ValueError('; '.join(
            _uncertified_power_message(branch, power, issue)
            for branch, issue in uncovered))
    if not matches:
        raise ValueError('terminal state matches no recorded branch')
    if len(matches) > 1:
        raise ValueError(
            'terminal state matches more than one recorded branch: {}'.format(matches))
    return matches[0]
