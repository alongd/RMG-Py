#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the "Software"),  #
# to deal in the Software without restriction, including without limitation  #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""Reactor-level ownership and caching for fingerprinted EEDF tables."""

import copy
from collections.abc import Mapping
from contextlib import contextmanager
from contextvars import ContextVar
import json
import os
from pathlib import Path
import re
import threading
import time
from types import MappingProxyType
import weakref

import numpy as np

from rmgpy.solver.eedf import (
    AmbiguousBranch,
    EEDFError,
    EEDFRow,
    EEDFTable,
    FingerprintMismatch,
    IllConditionedCoordinate,
    OutOfDomain,
)
from rmgpy.tools.eedf.schema import (
    canonical_json,
    content_hash,
    qualification_rule_for_quantity,
)


DEVELOPMENT_UNQUALIFIED_STATUS = (
    'DEVELOPMENT ONLY \u2014 TABLE QUALIFICATION FAILED')
DEVELOPMENT_DIAGNOSTIC_STATUS = 'DEVELOPMENT \u2014 not a qualification'
_development_unqualified_route = ContextVar(
    'eedf_development_unqualified_route', default=False)


class DevelopmentWallBudgetExceeded(RuntimeError):
    """The private development runner exhausted its configured wall budget."""


@contextmanager
def development_unqualified_table_route():
    """Temporarily admit an unqualified table for an explicit development run.

    This process-local context is deliberately absent from the input DSL and
    reactor constructor state.  It therefore cannot be selected by ``rmg.py``,
    serialized into an input file, or restored through ``__reduce__``.
    """
    token = _development_unqualified_route.set(True)
    try:
        yield
    finally:
        _development_unqualified_route.reset(token)


def unqualified_artifact_blockers(manifest):
    """Return the artifact's qualification blockers without changing them."""
    blockers = []
    verdicts = manifest.get('held_out_verdicts', [])
    failed = [verdict for verdict in verdicts if verdict.get('passed') is not True]
    if manifest.get('accepted') is not True:
        blockers.append('artifact manifest is not accepted')
    if failed:
        blockers.append(
            '{0} of {1} held-out verdicts failed'.format(len(failed), len(verdicts)))
    for branch, certification in sorted(manifest.get('branch_certification', {}).items()):
        if not str(certification).lower().startswith('certified'):
            blockers.append('{0}: {1}'.format(branch, certification))
    return tuple(blockers)


def _freeze(value):
    """Return a recursively immutable defensive value."""
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, np.ndarray):
        source = np.ascontiguousarray(value)
        # A read-only flag on an owning ndarray is reversible with
        # ``setflags(write=True)``.  Back the public row with immutable bytes so
        # callers cannot turn mutation back on and corrupt a cached evaluation.
        result = np.frombuffer(source.tobytes(), dtype=source.dtype).reshape(source.shape)
        return result
    return copy.deepcopy(value)


def _freeze_row(row):
    return EEDFRow(**{name: _freeze(value) for name, value in vars(row).items()})


class EEDFProvider:
    """Immutable reactor-owned view of one accepted table branch.

    Rows returned by this provider are immutable handles issued by this exact
    provider. Consumers cannot pass a row from another table or provider into
    rate, transport, or energy accessors.
    """

    def __init__(self, path, current_model, *, artifact_sha256, reaction_map,
                 branch, empirical_laws=None):
        development = bool(_development_unqualified_route.get())
        table = EEDFTable.load(path, current_model,
                               artifact_sha256=artifact_sha256,
                               require_accepted=not development)
        if branch not in table.manifest['branches']:
            raise AmbiguousBranch('declared branch ' + str(branch))
        energy_index = table._layout[branch]['swarm/mean_energy_eV'][0]
        if np.any(np.diff(table._data[branch][..., energy_index], axis=0) <= 0.0):
            raise IllConditionedCoordinate(
                'mean energy must increase strictly with u=ln(E/N)')
        self._table = table
        # Reactor callbacks need repeated read-only access to artifact metadata.
        # Freeze one private snapshot here so those callbacks never pay for the
        # public manifest property's deliberately defensive deep copy.
        self._runtime_metadata = _freeze(table.manifest)
        self._branch = str(branch)
        self._reaction_map = MappingProxyType(self._validate_reaction_map(reaction_map))
        if empirical_laws is not None and not isinstance(empirical_laws, dict):
            raise TypeError('empirical_laws must be a mapping')
        self._empirical_laws = _freeze(empirical_laws or {})
        if development and table.manifest.get('accepted') is True:
            raise EEDFError(
                'development unqualified-table route refuses an already qualified artifact')
        self._development_blockers = (unqualified_artifact_blockers(table.manifest)
                                      if development else ())
        self._cache_identity = None
        self._cached_row = None
        self._issued_rows = weakref.WeakValueDictionary()
        self._cache_lock = threading.RLock()
        self._qualification_status = None
        self._qualification_record = None
        self._qualification_state_resolver = None

    def _validate_reaction_map(self, reaction_map):
        if not isinstance(reaction_map, dict):
            raise TypeError('reaction_map must be a mapping')
        channels = self._runtime_metadata['channel_map']
        validated = {}
        for reaction_index, location in reaction_map.items():
            if (isinstance(reaction_index, bool) or
                    not isinstance(reaction_index, (int, np.integer)) or
                    reaction_index < 0):
                raise ValueError('invalid reaction index ' + repr(reaction_index))
            if not isinstance(location, (tuple, list)) or len(location) != 2:
                raise ValueError('reaction map entry must be (column, side)')
            column, side = location
            if (isinstance(column, bool) or
                    not isinstance(column, (int, np.integer)) or
                    not 0 <= column < len(channels)):
                raise ValueError('invalid channel column ' + repr(column))
            channel = channels[int(column)]
            if channel.get('classification') != 'A' or channel.get('reaction') is None:
                raise ValueError('channel column is not reaction-owned A chemistry')
            if side not in ('ine', 'sup'):
                raise ValueError('invalid channel side ' + repr(side))
            validated[int(reaction_index)] = (int(column), side)
        return validated

    @property
    def branch(self):
        return self._branch

    @property
    def reaction_map(self):
        return self._reaction_map

    @property
    def empirical_laws(self):
        return self._empirical_laws

    @property
    def manifest(self):
        """A defensive copy of the checked artifact manifest."""
        manifest = copy.deepcopy(self._table.manifest)
        notice = self.development_notice
        if notice is not None:
            manifest.update(notice)
        return manifest

    @property
    def runtime_metadata(self):
        """Immutable metadata for reactor-internal callback use."""
        return self._runtime_metadata

    @property
    def development_notice(self):
        """Serializable warning carried by every development diagnostic."""
        if not self._development_blockers:
            return None
        return {
            'scientific_status': DEVELOPMENT_UNQUALIFIED_STATUS,
            'artifact_blockers': list(self._development_blockers),
            'export_allowed': False,
            'qualification_allowed': False,
        }

    def require_qualified(self, operation):
        """Refuse export/qualification consumers for a development provider."""
        if self._development_blockers:
            raise EEDFError(
                '{0} refuses {1}: {2}'.format(
                    DEVELOPMENT_UNQUALIFIED_STATUS, operation,
                    '; '.join(self._development_blockers)))
        if self._qualification_status != 'PASS':
            status = self._qualification_status or 'NOT QUALIFIED'
            raise EEDFError('EEDF qualification {0} refuses {1}'.format(status, operation))

    def bind_qualification_state(self, resolver):
        """Bind the reactor callback that describes an accepted terminal state."""
        if not callable(resolver):
            raise TypeError('qualification state resolver must be callable')
        self._qualification_state_resolver = resolver

    @staticmethod
    def _json_value(value):
        if isinstance(value, MappingProxyType):
            value = dict(value)
        if isinstance(value, dict):
            return {str(key): EEDFProvider._json_value(item)
                    for key, item in value.items()}
        if isinstance(value, (tuple, list)):
            return [EEDFProvider._json_value(item) for item in value]
        if isinstance(value, np.ndarray):
            return EEDFProvider._json_value(value.tolist())
        if isinstance(value, np.generic):
            return EEDFProvider._json_value(value.item())
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, float) and not np.isfinite(value):
            if np.isnan(value):
                return 'NaN'
            return 'Infinity' if value > 0. else '-Infinity'
        return value

    def _qualification_state(self, y, t):
        if isinstance(y, dict):
            state = copy.deepcopy(y)
        elif self._qualification_state_resolver is not None:
            state = self._qualification_state_resolver(y, t)
        else:
            raise EEDFError('terminal qualification state is unavailable')
        for name in ('u', 'composition'):
            if name not in state:
                raise EEDFError('terminal qualification state has no ' + name)
        state['t'] = float(t)
        state['u'] = float(state['u'])
        state['composition'] = {str(key): float(value)
                                for key, value in state['composition'].items()}
        return state

    def _write_qualification(self, record, state):
        path = Path(state.get('record_path', 'eedf_qualification.json'))
        path.parent.mkdir(parents=True, exist_ok=True)
        self._qualification_record = record
        path.write_text(json.dumps(self._json_value(record), indent=2,
                                   sort_keys=True, allow_nan=False) + '\n')
        return path

    @staticmethod
    def _tolerance_check(checks, rule, quantity, actual, reference, tolerance,
                         labels=None, relative_scale=None):
        actual = np.asarray(actual, dtype=float)
        reference = np.asarray(reference, dtype=float)
        if actual.shape != reference.shape:
            raise EEDFError('{0} shape mismatch for {1}'.format(rule, quantity))
        error = np.abs(actual - reference)
        scale = np.abs(reference) if relative_scale is None else np.asarray(
            relative_scale, dtype=float)
        allowed = (np.asarray(tolerance['atol'], dtype=float) +
                   float(tolerance['rtol']) * scale)
        allowed = np.broadcast_to(allowed, error.shape)
        passed = np.isfinite(error) & (error <= allowed)
        if error.shape:
            scores = np.full(error.shape, np.inf, dtype=float)
            finite = np.isfinite(error) & np.isfinite(allowed)
            scores[finite] = error[finite] - allowed[finite]
            flat = int(np.argmax(scores))
            index = list(np.unravel_index(flat, error.shape))
            label = labels[flat] if labels is not None and flat < len(labels) else index
            absolute_error = float(error.flat[flat])
            allowed_error = float(allowed.flat[flat])
            relative_denominator = float(np.broadcast_to(scale, error.shape).flat[flat])
        else:
            index = []
            label = labels[0] if labels else quantity
            absolute_error = float(error)
            allowed_error = float(allowed)
            relative_denominator = float(scale)
        checks.append({
            'rule': rule,
            'quantity': quantity,
            'argmax': label,
            'index': index,
            'absolute_error': absolute_error,
            'relative_error': (absolute_error / abs(relative_denominator)
                               if relative_denominator and
                               np.isfinite(absolute_error) and
                               np.isfinite(relative_denominator) else None),
            'allowed_error': allowed_error,
            'passed': bool(np.all(passed)),
        })

    @staticmethod
    def _quantity_drift(quantity, actual, reference, labels=None):
        actual = np.asarray(actual, dtype=float)
        reference = np.asarray(reference, dtype=float)
        if actual.shape != reference.shape:
            return {'quantity': quantity, 'shape_mismatch': [list(actual.shape),
                                                              list(reference.shape)]}
        error = np.abs(actual - reference)
        if error.shape:
            flat = int(np.argmax(error))
            index = list(np.unravel_index(flat, error.shape))
            argmax = labels[flat] if labels is not None and flat < len(labels) else index
            absolute = float(error.flat[flat])
            ref = float(reference.flat[flat])
        else:
            index, argmax = [], quantity
            absolute, ref = float(error), float(reference)
        return {'quantity': quantity, 'argmax': argmax, 'index': index,
                'absolute_error': absolute,
                'relative_error': (absolute / abs(ref)
                                   if ref and np.isfinite(absolute) and
                                   np.isfinite(ref) else None)}

    def _all_quantity_drifts(self, predicted, direct):
        channels = [entry.get('description', str(i))
                    for i, entry in enumerate(self._runtime_metadata['channel_map'])]
        drifts = [self._quantity_drift('EN_Td', predicted['EN_Td'], direct['EN_Td'])]
        for name in sorted(set(predicted['swarm']) | set(direct['swarm'])):
            if name in predicted['swarm'] and name in direct['swarm']:
                drifts.append(self._quantity_drift(
                    'swarm.' + name, predicted['swarm'][name], direct['swarm'][name]))
        for name in ('k_ine', 'k_sup', 'channel_power', 'attachment_energy_eV',
                     'target_fractions', 'product_fractions'):
            if name in predicted and name in direct:
                drifts.append(self._quantity_drift(name, predicted[name], direct[name], channels))
        for name in sorted(set(predicted['power_groups']) | set(direct['power_groups'])):
            if name in predicted['power_groups'] and name in direct['power_groups']:
                drifts.append(self._quantity_drift(
                    'power.' + name, predicted['power_groups'][name],
                    direct['power_groups'][name]))
        if 'f0' in predicted and 'f0' in direct:
            drifts.append(self._quantity_drift('f0', predicted['f0'], direct['f0']))
        return drifts

    def _verify_generation_spec(self, spec, coordinates=None, fields_Td=None):
        """Verify that a saved direct-solver request owns this fingerprint."""
        from rmgpy.tools.eedf.loki import qualification_setup_sha256

        manifest = self._json_value(self._runtime_metadata)
        row_inputs = manifest['row_inputs']
        if content_hash(row_inputs) != manifest['fingerprint']:
            raise EEDFError('generation spec mismatch: manifest fingerprint')
        stored_setup = row_inputs.get('qualification_setup_sha256')
        if not stored_setup:
            raise EEDFError(
                'generation spec mismatch: artifact has no physical setup fingerprint')
        if manifest.get('qualification_setup_sha256') != stored_setup:
            raise EEDFError(
                'generation spec mismatch: physical setup fingerprint record')
        solver = manifest.get('solver')
        if not isinstance(solver, dict):
            raise EEDFError('generation spec mismatch: solver identity')

        def require_equal(label, actual, expected):
            if canonical_json(actual) != canonical_json(expected):
                raise EEDFError('generation spec mismatch: ' + label)

        require_equal('solver commit', spec.get('loki_commit'),
                      solver.get('commit'))
        require_equal('binary', spec.get('binary', {}).get('sha256'),
                      solver.get('binary_sha256'))
        require_equal('compiler', spec.get('compiler'), solver.get('compiler'))
        require_equal('cmake cache', spec.get('cmake_cache', {}).get('sha256'),
                      solver.get('cmake_cache_sha256'))
        spec_options = spec.get('solver_options')
        saved_options = solver.get('options')
        if (isinstance(spec_options, dict) and isinstance(saved_options, dict) and
                spec_options.get('ionizationOperatorType') !=
                saved_options.get('ionizationOperatorType')):
            raise EEDFError('generation spec mismatch: operator')
        require_equal('solver options', spec_options, saved_options)
        require_equal('working conditions', spec.get('working_conditions'),
                      solver.get('working_conditions'))
        spec_inputs = {
            name: {'sha256': item.get('sha256'), 'kind': item.get('kind')}
            for name, item in spec.get('input_files', {}).items()
        }
        require_equal('input SHAs', spec_inputs, solver.get('input_files'))

        cross_sections = {name: item['sha256']
                          for name, item in spec_inputs.items()
                          if item['kind'] == 'cross_section'}
        auxiliary = {name: item['sha256']
                     for name, item in spec_inputs.items()
                     if item['kind'] != 'cross_section'}
        for label, key, value in (
                ('solver commit', 'loki_commit', spec.get('loki_commit')),
                ('binary', 'binary_sha256',
                 spec.get('binary', {}).get('sha256')),
                ('solver options', 'solver_options', spec_options),
                ('working conditions', 'working_conditions',
                 spec.get('working_conditions')),
                ('qualification quantity rules', 'qualification_quantity_rules',
                 spec.get('qualification_quantity_rules')),
                ('input SHAs', 'cross_sections', cross_sections),
                ('input SHAs', 'auxiliary_inputs', auxiliary)):
            if key not in row_inputs:
                raise EEDFError(
                    'generation spec mismatch: fingerprint has no ' + key)
            require_equal(label, value, row_inputs[key])
        if qualification_setup_sha256(spec, coordinates, fields_Td) != stored_setup:
            raise EEDFError(
                'generation spec mismatch: physical setup fingerprint')

    @staticmethod
    def _declared_values(declarations, label):
        """Parse setup declarations numerically while retaining their names."""
        if not isinstance(declarations, list):
            raise EEDFError(label + ' declarations are invalid')
        values = {}
        pattern = re.compile(
            r'^\s*(.+?)\s*=\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*$')
        for declaration in declarations:
            match = pattern.fullmatch(str(declaration))
            if match is None:
                raise EEDFError(label + ' declaration is invalid')
            name, value = match.groups()
            value = float(value)
            if name in values or not np.isfinite(value):
                raise EEDFError(label + ' declaration is invalid')
            values[name] = value
        return values

    @staticmethod
    def _terminal_values(state, key, expected):
        values = state.get(key)
        if not isinstance(values, Mapping) or set(values) != set(expected):
            raise EEDFError('terminal ' + key.replace('_', ' ') + ' are incomplete')
        result = {str(name): float(value) for name, value in values.items()}
        if not all(np.isfinite(value) for value in result.values()):
            raise EEDFError('terminal ' + key.replace('_', ' ') + ' are incomplete')
        return result

    @staticmethod
    def _validate_qualification_row(row, label):
        """Reject an incomplete, non-finite, or non-converged solver row."""
        required = {
            'u', 'EN_Td', 'composition', 'swarm', 'gas_fractions',
            'state_populations', 'k_ine', 'k_sup', 'channel_power',
            'power_groups', 'attachment_energy_eV', 'target_fractions',
            'product_fractions', 'rate_floors', 'below_floor', 'f0',
            'energy_eV', 'energy_edges_eV', 'converged', 'iteration_count',
        }
        if not isinstance(row, Mapping):
            raise EEDFError('invalid ' + label + ' row: not a mapping')
        missing = sorted(required - set(row))
        if missing:
            raise EEDFError(
                'invalid ' + label + ' row: missing ' + ', '.join(missing))
        if row['converged'] is not True:
            raise EEDFError('invalid ' + label + ' row: converged')
        mappings = ('composition', 'swarm', 'gas_fractions',
                    'state_populations', 'power_groups')
        for name in mappings:
            values = row[name]
            if not isinstance(values, Mapping) or not values:
                if name == 'composition' and isinstance(values, Mapping):
                    continue
                raise EEDFError('invalid ' + label + ' row: ' + name)
            try:
                finite = all(np.all(np.isfinite(np.asarray(value, dtype=float)))
                             for value in values.values())
            except (TypeError, ValueError):
                finite = False
            if not finite:
                raise EEDFError('invalid ' + label + ' row: ' + name)
        numeric = required - set(mappings) - {'converged'}
        for name in numeric:
            try:
                finite = np.all(np.isfinite(np.asarray(row[name], dtype=float)))
            except (TypeError, ValueError):
                finite = False
            if not finite:
                raise EEDFError('invalid ' + label + ' row: ' + name)

    def _check_all_recorded_quantities(self, checks, predicted, direct, tolerances):
        """Apply only a-posteriori A budgets to every supplementary row value."""
        channels = [entry.get('description', str(i))
                    for i, entry in enumerate(self._runtime_metadata['channel_map'])]
        rules = self._runtime_metadata['row_inputs'].get(
            'qualification_quantity_rules')
        if not isinstance(rules, Mapping):
            raise EEDFError('artifact has no qualification quantity rules')

        def check(quantity, actual, reference, labels=None):
            rule = rules.get(quantity)
            allowed = qualification_rule_for_quantity(quantity)
            if rule != allowed:
                if rule is not None and allowed is not None:
                    raise EEDFError(
                        'invalid qualification rule for ' + quantity +
                        ': expected ' + allowed + ', got ' + str(rule))
                raise EEDFError(
                    'artifact has no qualification rule for ' + quantity)
            self._tolerance_check(checks, rule, quantity, actual, reference,
                                  tolerances[rule], labels)

        check('EN_Td', predicted['EN_Td'], direct['EN_Td'])
        predicted_names = set(predicted['swarm'])
        direct_names = set(direct['swarm'])
        if predicted_names != direct_names:
            missing = sorted(predicted_names ^ direct_names)
            raise EEDFError(
                'A2 quantity set mismatch for swarm: ' + ', '.join(missing))
        for name in sorted(predicted_names - {
                'mean_energy_eV', 'mobility_N', 'diffusion_N'}):
            check('swarm.' + name, predicted['swarm'][name],
                  direct['swarm'][name])
        predicted_names = set(predicted['power_groups'])
        direct_names = set(direct['power_groups'])
        if predicted_names != direct_names:
            missing = sorted(predicted_names ^ direct_names)
            raise EEDFError(
                'A4 quantity set mismatch for power_groups: ' + ', '.join(missing))
        for name in sorted(predicted_names):
            check('power.' + name, predicted['power_groups'][name],
                  direct['power_groups'][name])
        if 'channel_power' not in predicted or 'channel_power' not in direct:
            raise EEDFError('A5 quantity set mismatch: channel_power')
        predicted_field = abs(float(predicted['power_groups']['field']))
        direct_field = abs(float(direct['power_groups']['field']))
        if not predicted_field or not direct_field:
            raise EEDFError('A5 cannot scale zero field power')
        check('channel_power_fraction',
              np.asarray(predicted['channel_power']) / predicted_field,
              np.asarray(direct['channel_power']) / direct_field, channels)
        for name in ('attachment_energy_eV', 'target_fractions',
                     'product_fractions'):
            if name not in predicted or name not in direct:
                raise EEDFError('qualification quantity set mismatch: ' + name)
            check(name, predicted[name], direct[name], channels)
        if any(name not in predicted or name not in direct
               for name in ('f0', 'energy_eV', 'energy_edges_eV')):
            raise EEDFError('A1 quantity set mismatch: f0')
        if (not np.array_equal(predicted['energy_eV'], direct['energy_eV']) or
                not np.array_equal(predicted['energy_edges_eV'],
                                   direct['energy_edges_eV'])):
            raise EEDFError('A5 f0 energy grid mismatch')
        weights = (np.sqrt(np.asarray(direct['energy_eV'], dtype=float)) *
                   np.diff(np.asarray(direct['energy_edges_eV'], dtype=float)))
        distance = np.sum(
            np.abs(np.asarray(predicted['f0'], dtype=float) -
                   np.asarray(direct['f0'], dtype=float)) * weights)
        check('f0.weighted_L1', distance, 0.)

    def _compare_qualification_rows(self, predicted, direct, state, tolerances):
        required = tuple('A' + str(i) for i in range(1, 7))
        missing = [name for name in required if name not in tolerances]
        if missing:
            raise EEDFError('artifact has no runtime tolerance ' + ', '.join(missing))
        checks = []
        self._check_all_recorded_quantities(
            checks, predicted, direct, tolerances)
        channels = [entry.get('description', str(i))
                    for i, entry in enumerate(self._runtime_metadata['channel_map'])]
        self._tolerance_check(checks, 'A1', 'mean_energy_eV',
                              predicted['swarm']['mean_energy_eV'],
                              direct['swarm']['mean_energy_eV'], tolerances['A1'])
        swarm_names = ('mobility_N', 'diffusion_N')
        for name in swarm_names:
            if name not in predicted['swarm'] or name not in direct['swarm']:
                raise EEDFError('A2 quantity set mismatch: ' + name)
            self._tolerance_check(checks, 'A2', 'swarm.' + name,
                                  predicted['swarm'][name], direct['swarm'][name],
                                  tolerances['A2'])
        for name in ('k_ine', 'k_sup'):
            identities = []
            for index, channel in enumerate(self._runtime_metadata['channel_map']):
                match = re.fullmatch(
                    r'e \+ (.+?) (?:<->|->) e \+ (.+), [A-Za-z]+',
                    channel['description'])
                if match:
                    source, product = match.groups()
                else:
                    source = product = channel.get(
                        'flux_group', channel['description'])
                if name == 'k_sup':
                    source, product = product, source
                identities.append((source, product))
            floors = self._runtime_metadata['floors']

            def materiality(row):
                fractions = np.asarray(
                    row['target_fractions'] if name == 'k_ine'
                    else row['product_fractions'], dtype=float)
                rates = np.asarray(row[name], dtype=float)
                fluxes = np.abs(fractions * rates)
                losses = {}
                productions = {}
                for index, (source, product) in enumerate(identities):
                    losses[source] = losses.get(source, 0.) + fluxes[index]
                    productions[product] = (
                        productions.get(product, 0.) + fluxes[index])
                dynamic = float(floors.get('eedf_dynamic_range', 0.))
                sigma_max = np.maximum(
                    np.abs(np.asarray(row['k_ine'], dtype=float)),
                    np.abs(np.asarray(row['k_sup'], dtype=float)))
                rate_floors = np.maximum(
                    float(floors['rate_absolute']), dynamic * sigma_max)
                result = []
                for index, (source, product) in enumerate(identities):
                    totals = [value for value in
                              (losses[source], productions[product]) if value > 0.]
                    scale = min(totals, default=0.)
                    flux_floor = (
                        float(floors['absolute_flux_fraction']) * scale /
                        abs(fractions[index]) if fractions[index] else 0.)
                    load_bearing = (
                        scale > 0. and
                        fluxes[index] >=
                        float(floors['rate_flux_fraction']) * scale and
                        abs(rates[index]) >= rate_floors[index])
                    result.append((load_bearing, rate_floors[index], flux_floor))
                return result

            predicted_materiality = materiality(predicted)
            direct_materiality = materiality(direct)
            for index, channel in enumerate(self._runtime_metadata['channel_map']):
                classification = channel.get('classification')
                predicted_load = predicted_materiality[index][0]
                direct_load = direct_materiality[index][0]
                if classification != 'A' and not (
                        classification == 'B' and
                        (predicted_load or direct_load)):
                    continue
                tolerance = dict(tolerances['A3'])
                tolerance['atol'] = max(
                    float(tolerance['atol']),
                    float(predicted_materiality[index][1]),
                    float(direct_materiality[index][1]),
                    float(predicted_materiality[index][2]),
                    float(direct_materiality[index][2]))
                if not predicted_load and not direct_load:
                    tolerance['rtol'] = 0.
                self._tolerance_check(
                    checks, 'A3', name + ':' + str(index),
                    predicted[name][index], direct[name][index], tolerance,
                    [channels[index]])
        loss_names = ('elastic_gain', 'elastic_loss', 'car_gain', 'car_loss',
                      'excitation_loss', 'excitation_gain', 'vibrational_loss',
                      'vibrational_gain', 'rotational_loss', 'rotational_gain',
                      'ionization', 'attachment')
        self._tolerance_check(
            checks, 'A4', 'power.total_loss',
            -sum(float(predicted['power_groups'].get(name, 0.)) for name in loss_names),
            -sum(float(direct['power_groups'].get(name, 0.)) for name in loss_names),
            tolerances['A4'])
        field = abs(float(direct['power_groups']['field']))
        a5 = tolerances['A5']
        significant = [i for i, value in enumerate(direct['channel_power'])
                       if abs(value) >= float(a5.get('share_min', 0.)) * field]
        for index in significant:
            predicted_fraction = (float(predicted['channel_power'][index]) /
                                  abs(float(predicted['power_groups']['field'])))
            direct_fraction = float(direct['channel_power'][index]) / field
            self._tolerance_check(checks, 'A5', 'channel_power_fraction:' + str(index),
                                  predicted_fraction, direct_fraction, a5,
                                  [channels[index]])
        budget = state.get('energy_budget', {})
        if budget:
            predicted_field = float(predicted['power_groups']['field'])
            if predicted_field == 0.:
                raise EEDFError('A6 cannot scale a zero interpolated field power')
            scale = float(budget['joule_power']) / predicted_field
            direct_joule = scale * float(direct['power_groups']['field'])
            engine_target = (float(budget['P_abs']) -
                             float(budget.get('Q_wall_electron', 0.)) -
                             float(budget.get('Q_wall_ion', 0.)) -
                             float(budget.get('Q_flow', 0.)) +
                             scale * float(direct['power_groups'].get('growth', 0.)))
            self._tolerance_check(checks, 'A6', 'joule_power',
                                  direct_joule, engine_target,
                                  tolerances['A6'],
                                  relative_scale=abs(float(budget['P_abs'])))
        else:
            raise EEDFError('terminal qualification state has no energy_budget')
        return checks

    def _run_qualification(self, y, t, runner, development):
        from rmgpy.exceptions import PlasmaStateError
        from rmgpy.tools.eedf.loki import LoKIDriver, enrich_row

        start = time.monotonic()
        state = {}
        if isinstance(y, dict) and 'record_path' in y:
            state['record_path'] = y['record_path']
        record = {
            'status': DEVELOPMENT_DIAGNOSTIC_STATUS if development else 'FAIL',
            'comparison_passed': False,
            'terminal_state': state,
            'artifact_sha256': self._runtime_metadata['artifact_sha256'],
            'generation_spec': str(self._table.path / 'generation_spec.json'),
            'branch': self._branch,
        }
        try:
            state['t'] = float(t)
            state = self._qualification_state(y, t)
            predicted = self.row(state['u'], state['composition']).as_dict()
            self._validate_qualification_row(predicted, 'interpolated')
            record['terminal_state'] = state
            record['interpolated_row'] = predicted
            spec_path = self._table.path / 'generation_spec.json'
            spec = json.loads(spec_path.read_text())
            self._verify_generation_spec(spec)
            if 'scratch_root' in state:
                spec['scratch_root'] = str(state['scratch_root'])
            for key, target in (('gas_fractions', 'gas_properties'),
                                ('state_populations', 'state_properties')):
                field = 'fraction' if key == 'gas_fractions' else 'population'
                saved = self._declared_values(spec[target].get(field), key)
                terminal = self._terminal_values(state, key, saved)
                spec[target][field] = [
                    '{0} = {1:.17g}'.format(name, terminal[name])
                    for name in sorted(terminal)]
            saved_weights = self._declared_values(
                spec['state_properties'].get('statisticalWeight'),
                'state statistical weights')
            terminal_weights = self._terminal_values(
                state, 'state_statistical_weights', saved_weights)
            if any(terminal_weights[name] != saved_weights[name]
                   for name in saved_weights):
                raise EEDFError(
                    'generation spec mismatch: physical setup fingerprint '
                    '(terminal statistical weights)')
            # Fractions, populations and E/N are the only legitimate terminal
            # substitutions.  Recheck the digest against the setup that will
            # actually be rendered so a changed weight or other template value
            # cannot enter between provenance verification and execution.
            terminal_fields = [float(np.exp(state['u']))]
            self._verify_generation_spec(
                spec, state['composition'], terminal_fields)
            if runner is None:
                binary = spec.get('binary', {}).get('path')
                if not binary or not Path(binary).is_file() or not os.access(binary, os.X_OK):
                    raise EEDFError('missing or unexecutable LoKI-B binary')
                driver = LoKIDriver(spec)
                raw = driver.run(state['composition'], terminal_fields,
                                 'qualification-{0}'.format(time.time_ns()))[0]
                channel_map = self._json_value(self._runtime_metadata['channel_map'])
                direct = enrich_row(raw, channel_map, spec)
            else:
                direct = runner(spec, state['composition'],
                                terminal_fields, 'qualification')[0]
            record['solver_row'] = direct
            self._validate_qualification_row(direct, 'direct')
            tolerances = self._runtime_metadata['tolerances']
            if development and not all('A' + str(i) in tolerances for i in range(1, 7)):
                # Historical development artifacts predate A1--A6. Reuse their
                # explicit held-out budgets for measurement only; this route can
                # never confer qualification.
                tolerances = dict(tolerances)
                for target, source in (('A1', 'H1'), ('A2', 'H2'), ('A3', 'H3'),
                                       ('A4', 'H4'), ('A5', 'H4')):
                    if target not in tolerances and source in tolerances:
                        tolerances[target] = dict(tolerances[source])
                if 'A5' in tolerances:
                    tolerances['A5'].setdefault(
                        'share_min', self._runtime_metadata.get('floors', {}).get(
                            'relative_power_share', 0.))
                if 'A6' not in tolerances and 'energy_budget' in state:
                    budget = state['energy_budget']
                    tolerances['A6'] = {'rtol': float(budget['A6b_tolerance']), 'atol': 0.}
                record['tolerance_source'] = 'development artifact held-out budgets'
            checks = self._compare_qualification_rows(predicted, direct, state, tolerances)
            record.update(checks=checks,
                          quantity_drifts=self._all_quantity_drifts(predicted, direct),
                          drift_maxima={rule: max(
                              (item for item in checks if item['rule'] == rule),
                              key=lambda item: (item['absolute_error']
                                                if np.isfinite(item['absolute_error'])
                                                else float('inf')),
                              default=None)
                              for rule in sorted({item['rule'] for item in checks})})
            budget = state.get('energy_budget', {})
            if 'A6b_relative' in budget:
                record['engine_A6b'] = {
                    'relative_error': float(budget['A6b_relative']),
                    'tolerance': float(budget['A6b_tolerance']),
                    'passed': bool(float(budget['A6b_relative']) <=
                                   float(budget['A6b_tolerance'])),
                }
            record['comparison_passed'] = all(item['passed'] for item in checks)
            record['wall_time_s'] = time.monotonic() - start
            if development:
                self.annotate_development_diagnostic(record)
                self._write_qualification(record, state)
                return record
            failed = next((item for item in checks if not item['passed']), None)
            if failed is not None:
                raise PlasmaStateError('EEDF qualification failed: ' + failed['rule'])
            record['status'] = 'PASS'
            self._qualification_status = 'PASS'
            self._write_qualification(record, state)
            return record
        except Exception as error:
            record['wall_time_s'] = time.monotonic() - start
            record['error'] = str(error)
            if development:
                self.annotate_development_diagnostic(record)
                self._write_qualification(record, state)
                return record
            self._qualification_status = 'FAILED'
            self._write_qualification(record, state)
            if isinstance(error, PlasmaStateError):
                raise
            raise PlasmaStateError('EEDF qualification failed: ' + str(error)) from error

    def record_terminal_failure(self, y, t, error):
        """Replace any prior PASS with one recorded, terminal FAILED latch."""
        if self._qualification_status == 'FAILED' and self._qualification_record:
            return self._qualification_record
        state = {}
        if isinstance(y, dict) and 'record_path' in y:
            state['record_path'] = y['record_path']
        try:
            state = self._qualification_state(y, t)
        except Exception:
            try:
                state['t'] = float(t)
            except Exception:
                pass
        record = {
            'status': 'FAIL',
            'comparison_passed': False,
            'terminal_state': state,
            'artifact_sha256': self._runtime_metadata['artifact_sha256'],
            'generation_spec': str(self._table.path / 'generation_spec.json'),
            'branch': self._branch,
            'error': str(error),
        }
        self._qualification_status = 'FAILED'
        self._write_qualification(record, state)
        return record

    def qualify(self, y, t, *, runner=None):
        """Qualify one accepted terminal state with an independent direct solve."""
        from rmgpy.exceptions import PlasmaStateError

        if self._development_blockers:
            raise EEDFError(DEVELOPMENT_UNQUALIFIED_STATUS + ' refuses qualification')
        if self._qualification_status == 'FAILED':
            raise PlasmaStateError(
                'EEDF qualification failed: provider is latched FAILED')
        return self._run_qualification(y, t, runner, development=False)

    def development_diagnostic(self, y, t, *, runner=None):
        """Measure direct/interpolated drift without changing development status."""
        if not self._development_blockers:
            raise EEDFError('development diagnostic requires a development provider')
        return self._run_qualification(y, t, runner, development=True)

    def annotate_development_diagnostic(self, diagnostic):
        """Add the mandatory development warning to one diagnostic mapping."""
        notice = self.development_notice
        if notice is not None:
            diagnostic.update(notice)
        return diagnostic

    @property
    def axis_names(self):
        return tuple(self._table.axis_names)

    def axis(self, name):
        try:
            index = self._table.axis_names.index(name)
        except ValueError as error:
            raise KeyError(name) from error
        result = self._table.axes[index].copy()
        result.setflags(write=False)
        return result

    @property
    def cache_identity(self):
        with self._cache_lock:
            return self._cache_identity

    @property
    def cached_row(self):
        with self._cache_lock:
            return self._cached_row

    def _identity(self, u, composition):
        if not isinstance(composition, dict):
            raise OutOfDomain('composition must be a named coordinate mapping')
        return (float(u), tuple(sorted((name, float(value))
                                       for name, value in composition.items())))

    def domain_check(self, u, composition, y=None, context=None):
        """Validate an accepted state; y/context are reactor-hook context."""
        try:
            return self._table.domain_check(u, composition)
        except OutOfDomain:
            values = [u]
            if isinstance(composition, dict):
                values.extend(composition.get(name) for name in self.axis_names[1:])
            if len(values) == len(self.axis_names):
                for name, axis, value in zip(self.axis_names, self._table.axes, values):
                    if (value is not None and
                            (not np.isfinite(value) or not axis[0] <= value <= axis[-1])):
                        raise OutOfDomain(
                            '{}={!r} outside [{!r}, {!r}]'.format(
                                name, value, axis[0], axis[-1]))
            raise

    def row(self, u, composition, y=None, context=None):
        """Return this provider's immutable row handle for one full identity."""
        with self._cache_lock:
            self.domain_check(u, composition, y=y, context=context)
            identity = self._identity(u, composition)
            cached = self._cached_row
            if identity == self._cache_identity and cached is not None:
                return cached
            local_row = _freeze_row(self._table.row(u, composition, self._branch))
            if local_row.branch_id != self._branch:
                raise AmbiguousBranch('row branch differs from declared branch')
            self._issued_rows[id(local_row)] = local_row
            # Publish the row and its full identity as one critical section.
            self._cached_row = local_row
            self._cache_identity = identity
            return local_row

    def _require_row(self, row):
        with self._cache_lock:
            if row is None:
                row = self._cached_row
            if row is None:
                raise RuntimeError('EEDF row has not been evaluated')
            if self._issued_rows.get(id(row)) is not row:
                raise FingerprintMismatch('row was not issued by this EEDF provider')
            if row.branch_id != self._branch:
                raise AmbiguousBranch('row branch differs from declared branch')
            if row.fingerprint != self._runtime_metadata['fingerprint']:
                raise FingerprintMismatch('row fingerprint')
            return row

    def reaction_rate(self, reaction_index, row):
        """Return a mapped SI rate from an explicit provider-issued row."""
        if reaction_index not in self._reaction_map:
            raise KeyError('unmapped reaction index ' + repr(reaction_index))
        row = self._require_row(row)
        column, side = self._reaction_map[reaction_index]
        return float(getattr(row, 'k_' + side)[column])

    def transport(self, row):
        """Return immutable reduced transport from an explicit row handle."""
        return self._require_row(row).swarm

    def mean_energy_eV(self, row):
        """Return mean electron energy from an explicit row handle."""
        row = self._require_row(row)
        try:
            return float(row.swarm['mean_energy_eV'])
        except KeyError as error:
            raise EEDFError('row has no mean_energy_eV') from error
