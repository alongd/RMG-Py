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

"""Shared EEDF provenance and generation-spec validation.

Tolerance values are supplied by the spec; this module has no numerical defaults.
"""
import copy
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import yaml


class EEDFError(ValueError):
    """A table or generation request is scientifically inadmissible."""


class FingerprintMismatch(EEDFError):
    """Row inputs differ from the qualified inputs."""


class SpecError(EEDFError):
    """An incomplete or inconsistent generation specification."""


def canonical_json(value):
    """Serialize manifest values independently of dictionary order."""
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def content_hash(value):
    """Hash content, excluding no fields implicitly."""
    return hashlib.sha256(canonical_json(value).encode()).hexdigest()


def qualification_setup_identity(spec):
    """Return the complete direct-solver setup identity with exact-state inputs blanked.

    Start from the whole generation spec so new physical fields cannot silently
    escape the digest.  Remove only non-acceptance policy and execution locations,
    then remove the three values that terminal qualification must replace:
    gas fractions, state populations, and the reduced-field scan.
    """
    identity = copy.deepcopy(spec)
    for name in ('output_root', 'scratch_root', 'omp_threads', 'max_parallel',
                 'timeout_s', 'repositories', 'screen', 'held_out',
                 'refinement'):
        identity.pop(name, None)
    for name in ('binary', 'cmake_cache', 'channel_map'):
        if name in identity:
            identity[name].pop('path', None)
    for item in identity.get('input_files', {}).values():
        item.pop('path', None)
    shared = identity.get('shared_objects')
    if isinstance(shared, dict):
        normalized = {}
        for path, digest in shared.items():
            name = Path(path).name
            if name in normalized and normalized[name] != digest:
                raise SpecError('duplicate shared-object basename in setup identity')
            normalized[name] = digest
        identity['shared_objects'] = normalized
    for properties, field in (('gas_properties', 'fraction'),
                              ('state_properties', 'population')):
        if field in identity.get(properties, {}):
            # The rendered setup below retains the resolved declaration names.
            # Blank the spec-side source wholesale so direct declarations and
            # hashed property-file declarations have the same terminal seam.
            identity[properties][field] = '<terminal>'
    if 'reducedField' in identity.get('working_conditions', {}):
        identity['working_conditions']['reducedField'] = '<terminal>'
    return _canonical_physical_value(identity)


_DECLARATION = re.compile(
    r'^\s*(.+?)\s*=\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*$')


def _canonical_physical_value(value):
    """Normalize numerical declaration spellings throughout a physical spec."""
    if isinstance(value, dict):
        return {key: _canonical_physical_value(item)
                for key, item in value.items()}
    if isinstance(value, list):
        return [_canonical_physical_value(item) for item in value]
    if isinstance(value, str):
        match = _DECLARATION.fullmatch(value)
        if match is not None:
            return {'name': match.group(1), 'value': float(match.group(2))}
    return value


def qualification_rule_for_quantity(quantity):
    """Return the dimensionally allowed a-posteriori rule for a row quantity."""
    if quantity == 'attachment_energy_eV':
        return 'A1'
    if quantity == 'EN_Td' or quantity.startswith('swarm.'):
        return 'A2'
    if quantity.startswith('power.'):
        return 'A4'
    if quantity in {'channel_power_fraction', 'target_fractions',
                    'product_fractions', 'f0.weighted_L1'}:
        return 'A5'
    return None


def file_hash(path):
    """SHA-256 of the bytes of an input or artifact."""
    with open(path, 'rb') as stream:
        digest = hashlib.sha256()
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def check_model(manifest, current, relative_tolerance):
    """Check all row inputs. Repository SHAs belong only to provenance.

    Tg and pressure allow the explicitly recorded relative comparison tolerance.
    All other fields, including arrays and dictionary key sets, compare exactly.
    """
    stored = manifest['row_inputs']
    if content_hash(stored) != manifest['fingerprint']:
        raise FingerprintMismatch('manifest fingerprint')
    current = copy.deepcopy(current)
    if set(current) != set(stored):
        raise FingerprintMismatch('row_inputs fields')
    for key, value in stored.items():
        if key in ('Tg_K', 'P_Pa'):
            if not np.isclose(current[key], value, rtol=relative_tolerance, atol=0):
                raise FingerprintMismatch(key)
        elif canonical_json(current[key]) != canonical_json(value):
            raise FingerprintMismatch(key)


def load_spec(path):
    """Load and validate a YAML/JSON spec; paths resolve beside that spec."""
    path = Path(path).resolve()
    with path.open() as stream:
        spec = json.load(stream) if path.suffix == '.json' else yaml.safe_load(stream)
    validate_spec(spec)
    spec = copy.deepcopy(spec)
    spec['shared_objects'] = {str((path.parent / name).resolve()): digest for name, digest in spec['shared_objects'].items()}
    for key in ('binary', 'cmake_cache', 'channel_map'):
        spec[key]['path'] = str((path.parent / spec[key]['path']).resolve())
    for item in spec['input_files'].values():
        item['path'] = str((path.parent / item['path']).resolve())
    for key in ('output_root', 'scratch_root'):
        spec[key] = str((path.parent / spec[key]).resolve())
    return spec


def validate_spec(spec):
    """Reject unspecified policy instead of freezing tolerances in code."""
    required = {'schema_version', 'arm', 'Tg_K', 'P_Pa', 'binary', 'shared_objects', 'cmake_cache',
                'compiler', 'loki_commit', 'input_files', 'channel_map', 'axes',
                'screen', 'envelopes', 'gas_properties', 'state_properties',
                'working_conditions', 'solver_options', 'tolerances', 'floors',
                'held_out', 'repositories', 'output_root', 'scratch_root',
                'omp_threads', 'max_parallel', 'timeout_s', 'refinement'}
    optional = {'engine_state_map', 'qualification_quantity_rules'}
    if (not isinstance(spec, dict) or not required <= set(spec) or
            set(spec) - required - optional):
        fields = ((required - set(spec or {})) |
                  (set(spec or {}) - required - optional))
        raise SpecError('spec fields: ' + str(fields))
    if 'engine_state_map' in spec:
        fields = {'formula', 'electronic_state', 'vibrational_level',
                  'multiplicity', 'loki_state', 'statistical_weight'}
        if (not isinstance(spec['engine_state_map'], list) or
                not spec['engine_state_map']):
            raise SpecError('engine_state_map')
        for entry in spec['engine_state_map']:
            if (not isinstance(entry, dict) or set(entry) != fields or
                    not isinstance(entry['formula'], str) or
                    not isinstance(entry['electronic_state'], str) or
                    not isinstance(entry['loki_state'], str) or
                    not isinstance(entry['multiplicity'], int) or
                    entry['multiplicity'] < 1 or
                    not isinstance(entry['statistical_weight'], (int, float)) or
                    not np.isfinite(entry['statistical_weight']) or
                    entry['statistical_weight'] <= 0):
                raise SpecError('engine_state_map')
    if not isinstance(spec['shared_objects'], dict) or not spec['shared_objects']:
        raise SpecError('shared_objects pins required')
    try:
        numeric = [spec['Tg_K'], spec['P_Pa'], spec['timeout_s']]
        numeric.extend(spec['floors'].values())
        for entry in spec['tolerances'].values():
            numeric.extend(v for k, v in entry.items() if k != 'norm') if isinstance(entry, dict) else numeric.append(entry)
        if any(not isinstance(v, (int, float)) or not np.isfinite(v) for v in numeric):
            raise SpecError('numeric policy must be finite numbers')
    except (TypeError, ValueError) as error:
        raise SpecError('numeric policy') from error
    if spec['schema_version'] != 1:
        raise SpecError('schema_version')
    for name in ('Tg_K', 'P_Pa', 'timeout_s'):
        if not np.isfinite(spec[name]) or spec[name] <= 0:
            raise SpecError(name)
    if (type(spec['max_parallel']) is not int or type(spec['omp_threads']) is not int
            or not 1 <= spec['max_parallel'] <= 5 or spec['omp_threads'] < 1):
        raise SpecError('parallelism')
    axes = spec['axes']
    if not axes or next(iter(axes)) != 'u':
        raise SpecError('u must be the first axis')
    for name, values in axes.items():
        a = np.asarray(values, dtype=float)
        if a.ndim != 1 or len(a) < 2 or not np.all(np.isfinite(a)) or np.any(np.diff(a) <= 0):
            raise SpecError('axis ' + name)
    screen = spec['screen']
    if set(screen) != {'rule', 'threshold_fraction', 'candidates'} or screen['rule'] != 'any_quantity_exceeds':
        raise SpecError('screen rule')
    if screen['threshold_fraction'] <= 0:
        raise SpecError('screen threshold_fraction')
    for candidate in screen['candidates']:
        if set(candidate) != {'name', 'reference', 'values'} or len(candidate['values']) < 2:
            raise SpecError('screen candidate')
        if candidate['name'] not in axes and candidate['name'] not in spec['envelopes']:
            raise SpecError('screen coordinate lacks axis/envelope')
    for name, bounds in spec['envelopes'].items():
        if set(bounds) != {'reference', 'min', 'max'} or not bounds['min'] <= bounds['reference'] <= bounds['max']:
            raise SpecError('envelope ' + name)
    options = spec['solver_options']
    if (options.get('eedfType') != 'boltzmann' or
            options.get('ionizationOperatorType') != 'usingSDCS' or
            options.get('growthModelType') != 'temporal'):
        raise SpecError('LoKI model/operator/growth convention')
    numerics = options['numerics']
    for value in (numerics['maxPowerBalanceRelError'],
                  numerics['nonLinearRoutines']['maxEedfRelError']):
        if value <= 0:
            raise SpecError('convergence tolerance')
    tolerances = spec['tolerances']
    base_tolerances = {'G1', 'H1', 'H2', 'H3', 'H4', 'u_condition',
                       'fingerprint_rtol', 'normalization_atol', 'F0',
                       'channel_power_sum'}
    qualification_tolerances = {'A1', 'A2', 'A3', 'A4', 'A5', 'A6'}
    extra_tolerances = set(tolerances) - base_tolerances
    if (not base_tolerances <= set(tolerances) or
            extra_tolerances not in (set(), qualification_tolerances)):
        raise SpecError('tolerances fields')
    for key in ('G1', 'H1', 'H2', 'H3', 'H4'):
        if set(tolerances[key]) != ({'rtol', 'atol', 'total_rtol'} if key == 'H4' else {'rtol', 'atol'}) or any(v < 0 or not np.isfinite(v) for v in tolerances[key].values()):
            raise SpecError(key)
    for key in qualification_tolerances & set(tolerances):
        fields = {'rtol', 'atol', 'share_min'} if key == 'A5' else {'rtol', 'atol'}
        if (set(tolerances[key]) != fields or
                any(v < 0 or not np.isfinite(v)
                    for v in tolerances[key].values())):
            raise SpecError(key)
    rules = spec.get('qualification_quantity_rules')
    if qualification_tolerances <= set(tolerances):
        if (not isinstance(rules, dict) or not rules or
                any(not isinstance(quantity, str) or not quantity or
                    rule not in qualification_tolerances
                    for quantity, rule in rules.items())):
            raise SpecError('qualification_quantity_rules')
        invalid = [quantity for quantity, rule in rules.items()
                   if qualification_rule_for_quantity(quantity) != rule]
        if invalid:
            raise SpecError(
                'invalid qualification rule for ' + ', '.join(sorted(invalid)))
    elif rules is not None:
        raise SpecError('qualification_quantity_rules require A1-A6')
    if set(tolerances['u_condition']) != {'min_abs_denergy_du', 'min_scaled_denergy_du'}:
        raise SpecError('u_condition')
    if any(v < 0 for v in tolerances['u_condition'].values()):
        raise SpecError('u_condition thresholds')
    if any(tolerances[k] < 0 for k in ('fingerprint_rtol', 'normalization_atol')):
        raise SpecError('tolerance')
    if set(tolerances['F0']) != {'norm', 'atol'} or tolerances['F0']['norm'] != 'weighted_L1' or tolerances['F0']['atol'] < 0:
        raise SpecError('F0 norm/tolerance')
    if set(tolerances['channel_power_sum']) != {'rtol', 'atol'} or any(v < 0 for v in tolerances['channel_power_sum'].values()):
        raise SpecError('channel_power_sum')
    frequency = spec['working_conditions'].get('excitationFrequency')
    if not isinstance(frequency, (int, float)) or frequency != 0:
        raise SpecError('DC only: excitationFrequency')
    grid = options['numerics']['energyGrid']
    if set(grid) != {'cellNumber', 'maxEnergy'} or grid['cellNumber'] < 2 or grid['maxEnergy'] <= 0:
        raise SpecError('fixed uniform energy grid required')
    floors = spec['floors']
    if set(floors) != {'rate_absolute', 'eedf_dynamic_range', 'rate_flux_fraction',
                       'absolute_flux_fraction', 'relative_power_share', 'absolute_power_share'}:
        raise SpecError('floors fields')
    if floors['rate_absolute'] <= 0 or any(v < 0 or not np.isfinite(v) for v in floors.values()):
        raise SpecError('floors')
    if set(spec['held_out']) != {'lhs_count', 'seed', 'cell_midpoints'} or not spec['held_out']['cell_midpoints']:
        raise SpecError('held_out scheme')
    if spec['held_out']['lhs_count'] < 0:
        raise SpecError('held_out lhs_count')
    if set(spec['refinement']) != {'max_rounds'} or spec['refinement']['max_rounds'] < 0:
        raise SpecError('refinement')


def interpolant_identity():
    """Version and source hash checked at runtime."""
    from rmgpy.tools.eedf import __version__
    return {'version': __version__, 'source_sha256': file_hash(Path(__file__).parents[2] / 'solver' / 'eedf.py')}


def row_inputs(spec):
    """Content-address the physical inputs and numerical row-generation choices."""
    result = {key: copy.deepcopy(spec[key]) for key in (
        'arm', 'Tg_K', 'P_Pa', 'loki_commit', 'axes', 'gas_properties',
        'state_properties', 'working_conditions', 'solver_options', 'screen', 'envelopes',
        'floors', 'tolerances', 'held_out', 'refinement')}
    if 'engine_state_map' in spec:
        result['engine_state_map'] = copy.deepcopy(spec['engine_state_map'])
    if 'qualification_quantity_rules' in spec:
        result['qualification_quantity_rules'] = copy.deepcopy(
            spec['qualification_quantity_rules'])
    return result
