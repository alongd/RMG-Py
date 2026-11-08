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
from string import Formatter

import numpy as np

from rmgpy.diagnostic_mode import MAXWELLIAN_MODE
import yaml


class EEDFError(ValueError):
    """A table or generation request is scientifically inadmissible."""


class FingerprintMismatch(EEDFError):
    """Row inputs differ from the qualified inputs."""


class SpecError(EEDFError):
    """An incomplete or inconsistent generation specification."""


def _provider_kind(value):
    """Validate the literal provider name without importing a row generator."""
    if not isinstance(value, str) or value not in ('loki-table', 'maxwellian-table'):
        raise FingerprintMismatch('unknown EEDF provider kind')
    return value


def _artifact_provider_kind(manifest):
    """Old untagged artifacts are LoKI only; comparator identity is explicit."""
    kind = _provider_kind(manifest.get('provider_kind', 'loki-table'))
    inputs = manifest.get('row_inputs', {})
    if 'provider_kind' in inputs and inputs['provider_kind'] != kind:
        raise FingerprintMismatch('provider kind differs from row inputs')
    if kind == 'maxwellian-table' and inputs.get('provider_kind') != kind:
        raise FingerprintMismatch('comparator provider kind missing from row inputs')
    return kind


COMPARATOR_NA = ('branch_scan', 'self_consistency', 'field', 'A6b-source',
                 'A6b-runtime', 'T10 interpolation qualification', 'terminal_qualification')


def _comparator_metadata(manifest):
    """Validate numerical usability separately from production acceptance.

    This runtime contract contains no generator imports. Its complete return
    value is bound to the HDF5 bytes by the writer and checked by the loader.
    """
    metadata = manifest.get('comparator')
    if not isinstance(metadata, dict) or set(metadata) != {'source', 'applicability', 'numerical_validation', 'label'}:
        raise FingerprintMismatch('comparator metadata')
    if manifest.get('accepted') is not False:
        raise FingerprintMismatch('comparator cannot carry production acceptance')
    source = metadata['source']
    if (not isinstance(source, dict) or source != manifest['row_inputs'].get('comparator_source')
            or not isinstance(source.get('artifact_sha256'), str)
            or re.fullmatch('[0-9a-f]{64}', source['artifact_sha256']) is None
            or not isinstance(source.get('branch_id'), str)
            or source.get('retained_u') != manifest['axes']['u']
            or not isinstance(source.get('excluded_u'), list)):
        raise FingerprintMismatch('comparator source provenance')
    applicability = metadata['applicability']
    if not isinstance(applicability, dict) or set(applicability) != set(COMPARATOR_NA):
        raise FingerprintMismatch('comparator applicability')
    for name, record in applicability.items():
        if (not isinstance(record, dict) or set(record) != {'status', 'reason'}
                or record['status'] != 'NOT_APPLICABLE'
                or not isinstance(record['reason'], str) or not record['reason']):
            raise FingerprintMismatch('comparator applicability: ' + name)
    for name in ('A6b-source', 'A6b-runtime', 'T10 interpolation qualification'):
        if manifest.get(name) != applicability[name]:
            raise FingerprintMismatch('comparator production status: ' + name)
    if (manifest.get('branches') != ['maxwellian']
            or manifest.get('branch_certification') != {
                'maxwellian': 'NOT_APPLICABLE: prescribed Maxwellian comparator'}
            or manifest.get('branch_detection') != applicability['branch_scan']
            or manifest.get('held_out_verdicts') != []):
        raise FingerprintMismatch('comparator production branch claims')
    validation = metadata['numerical_validation']
    if (not isinstance(validation, dict) or validation.get('passed') is not True
            or validation.get('training_count') != int(np.prod([len(a) for a in manifest['axes'].values()]))
            or validation.get('held_out_count') != len(manifest.get('held_out', []))
            or not manifest.get('held_out')):
        raise FingerprintMismatch('comparator numerical validation')
    verdicts = validation.get('verdicts')
    if not isinstance(verdicts, list) or len(verdicts) != validation['held_out_count']:
        raise FingerprintMismatch('comparator numerical verdicts')
    for point, verdict in zip(manifest['held_out'], verdicts):
        if (verdict.get('point') != point or verdict.get('passed') is not True
                or not verdict.get('checks')):
            raise FingerprintMismatch('comparator numerical verdict')
        quantities = {'mean_energy_eV', 'mobility_N', 'diffusion_N', 'energy_mobility_N',
                      'energy_diffusion_N', 'k_ine', 'k_sup', 'attachment_energy_eV',
                      'channel_power', 'collision_loss', 'f0_weighted_L1'}
        if (len(verdict['checks']) != len(quantities)
                or {check.get('quantity') for check in verdict['checks']} != quantities):
            raise FingerprintMismatch('comparator numerical quantities')
        for check in verdict['checks']:
            error = check.get('max_tolerance_fraction')
            if (check.get('passed') is not True or type(error) not in (int, float)
                    or not np.isfinite(error) or not 0 <= error <= 1):
                raise FingerprintMismatch('comparator numerical check')
    resolution = validation.get('max_relative_mean_error')
    if (type(resolution) not in (int, float) or not np.isfinite(resolution)
            or not 0 <= resolution <= manifest['tolerances']['H2']['rtol']):
        raise FingerprintMismatch('comparator mean resolution')
    if metadata['label'] != MAXWELLIAN_MODE:
        raise FingerprintMismatch('comparator label')
    return copy.deepcopy(metadata)


def _canonical_json(value):
    """Serialize manifest values independently of dictionary order."""
    value = copy.deepcopy(value)
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def _content_hash(value):
    """Hash content, excluding no fields implicitly."""
    return hashlib.sha256(_canonical_json(value).encode()).hexdigest()


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
    for name in terminal_working_condition_fields(identity):
        identity['working_conditions'][name] = '<terminal>'
    if 'reducedField' in identity.get('working_conditions', {}):
        identity['working_conditions']['reducedField'] = '<terminal>'
    return _canonical_physical_value(identity)


def terminal_working_condition_fields(spec):
    """Return working-condition fields rendered from declared envelope values."""
    envelopes = set(spec.get('envelopes', {}))
    result = []
    for name, value in spec.get('working_conditions', {}).items():
        if not isinstance(value, str):
            continue
        fields = {
            field_name.split('.', 1)[0].split('[', 1)[0]
            for _, field_name, _, _ in Formatter().parse(value)
            if field_name
        }
        if fields and fields <= envelopes:
            result.append(name)
    return tuple(result)


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


def _file_hash(path):
    """SHA-256 of the bytes of an input or artifact."""
    path = Path(copy.deepcopy(path))
    with open(path, 'rb') as stream:
        digest = hashlib.sha256()
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


ARGON_VALIDITY_FIELDS = ('excited_fraction_max',
                         'argon_elimination_electron_density_max_m3')


def _argon_validity_declarations(document, error=SpecError):
    """Validate optional limits for a ground-population argon EEDF table.

    The excited fraction is total neutral Ar(4s) / neutral Ar ground; the
    electron-density limit applies to radiatively eliminated upper levels.
    Absence preserves tables that do not use this elimination convention.
    """
    result = {}
    for key in ARGON_VALIDITY_FIELDS:
        if key in document:
            value = document[key]
            if type(value) not in (int, float) or not np.isfinite(value) or value <= 0:
                raise error(key + ' must be finite and positive')
            result[key] = value
    return result


def _check_model(manifest, current, relative_tolerance):
    """Check all row inputs. Repository SHAs belong only to provenance.

    Tg and pressure allow the explicitly recorded relative comparison tolerance.
    All other fields, including arrays and dictionary key sets, compare exactly.
    """
    stored = manifest['row_inputs']
    if _content_hash(stored) != manifest['fingerprint']:
        raise FingerprintMismatch('manifest fingerprint')
    declarations = _argon_validity_declarations(manifest, FingerprintMismatch)
    recorded = _argon_validity_declarations(stored, FingerprintMismatch)
    for key in ARGON_VALIDITY_FIELDS:
        if declarations.get(key) != recorded.get(key):
            raise FingerprintMismatch(key + ' differs from row_inputs fingerprint')
    from rmgpy.tools.eedf.ownership import is_branch_family, validate_family
    mapped = [channel['reaction'] for channel in manifest.get('channel_map', [])
              if channel.get('reaction') is not None]
    if any(is_branch_family(identity) for identity in mapped + stored.get('reactions', [])):
        for identity in mapped:
            if is_branch_family(identity):
                validate_family(identity)
        if _canonical_json(mapped) != _canonical_json(stored.get('reactions')):
            raise FingerprintMismatch('branch ownership differs from row_inputs fingerprint')
    current = copy.deepcopy(current)
    if set(current) != set(stored):
        raise FingerprintMismatch('row_inputs fields')
    for key, value in stored.items():
        if key in ('Tg_K', 'P_Pa'):
            if not np.isclose(current[key], value, rtol=relative_tolerance, atol=0):
                raise FingerprintMismatch(key)
        elif _canonical_json(current[key]) != _canonical_json(value):
            raise FingerprintMismatch(key)


def load_spec(path):
    """Load and validate a YAML/JSON spec; paths resolve beside that spec."""
    path = Path(path).resolve()
    with path.open() as stream:
        spec = json.load(stream) if path.suffix == '.json' else yaml.safe_load(stream)
    _validate_spec(spec)
    spec = copy.deepcopy(spec)
    spec['shared_objects'] = {str((path.parent / name).resolve()): digest for name, digest in spec['shared_objects'].items()}
    for key in ('binary', 'cmake_cache', 'channel_map'):
        spec[key]['path'] = str((path.parent / spec[key]['path']).resolve())
    for item in spec['input_files'].values():
        item['path'] = str((path.parent / item['path']).resolve())
    for key in ('output_root', 'scratch_root'):
        spec[key] = str((path.parent / spec[key]).resolve())
    return spec


def _validate_magnitude_tolerances(policy):
    """Require explicit, finite reactor-impact budgets and their time sources."""
    numeric = {'species_epsilon', 'electron_epsilon', 'energy_epsilon',
               'neutral_floor_fraction', 'electron_fraction_max', 'time_max_s',
               'target_fraction_floor', 'underflow_absolute'}
    if (not isinstance(policy, dict) or set(policy) != numeric | {'time_sources'}
            or any(type(policy[k]) not in (int, float) or
                   not np.isfinite(policy[k]) or policy[k] <= 0 for k in numeric)):
        raise SpecError('magnitude_tolerances requires finite positive parameters')
    sources = policy['time_sources']
    if (not isinstance(sources, list) or not sources
            or any(not isinstance(s, str) or not s.strip() for s in sources)):
        raise SpecError('magnitude_tolerances requires time_sources')


def _validate_spec(spec):
    """Reject unspecified policy instead of freezing tolerances in code."""
    required = {'schema_version', 'arm', 'Tg_K', 'P_Pa', 'binary', 'shared_objects', 'cmake_cache',
                'compiler', 'loki_commit', 'input_files', 'channel_map', 'axes',
                'screen', 'envelopes', 'gas_properties', 'state_properties',
                'working_conditions', 'solver_options', 'tolerances', 'floors',
                'held_out', 'repositories', 'output_root', 'scratch_root',
                'omp_threads', 'max_parallel', 'timeout_s', 'refinement'}
    optional = {'engine_state_map', 'qualification_quantity_rules', 'vibrational_ladders', 'argon_4s_population', 'interpolation_accuracy', 'magnitude_tolerances',
                *ARGON_VALIDITY_FIELDS}
    if (not isinstance(spec, dict) or not required <= set(spec) or
            set(spec) - required - optional):
        fields = ((required - set(spec or {})) |
                  (set(spec or {}) - required - optional))
        raise SpecError('spec fields: ' + str(fields))
    _argon_validity_declarations(spec)
    from rmgpy.tools.eedf.argon_population import validate_population_axis
    validate_population_axis(spec)
    if 'magnitude_tolerances' in spec:
        _validate_magnitude_tolerances(spec['magnitude_tolerances'])
    if 'interpolation_accuracy' in spec:
        from rmgpy.tools.eedf.validation import _freeze_accuracy_criteria
        _freeze_accuracy_criteria(spec['interpolation_accuracy'])
    if 'vibrational_ladders' in spec:
        from rmgpy.tools.eedf.vibrational import validate_ladders
        validate_ladders(spec['vibrational_ladders'])
        if any(item['axis'] not in spec['axes'] for item in spec['vibrational_ladders']):
            raise SpecError('vibrational ladder requires its Tv axis')
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
    held_out_fields = set(spec['held_out'])
    if (held_out_fields not in ({'lhs_count', 'seed', 'cell_midpoints'},
                                {'lhs_count', 'seed', 'cell_midpoints', 'off_centre'})
            or not spec['held_out']['cell_midpoints']):
        raise SpecError('held_out scheme')
    if spec['held_out']['lhs_count'] < 0:
        raise SpecError('held_out lhs_count')
    for point in spec['held_out'].get('off_centre', []):
        if (not isinstance(point, dict) or set(point) != set(axes)
                or any(not isinstance(point[name], (int, float)) or
                       not np.isfinite(point[name]) or
                       not axes[name][0] <= point[name] <= axes[name][-1]
                       for name in axes)):
            raise SpecError('held_out off_centre')
    if set(spec['refinement']) != {'max_rounds'} or spec['refinement']['max_rounds'] < 0:
        raise SpecError('refinement')


def _interpolant_identity():
    """Version and source hash checked at runtime."""
    from rmgpy.tools.eedf import __version__
    return {'version': __version__, 'source_sha256': _content_hash({
        'table': _file_hash(Path(__file__).parents[2] / 'solver' / 'eedf.py'),
        'argon_population': _file_hash(Path(__file__).with_name('argon_population.py'))})}


def _row_inputs(spec):
    """Content-address the physical inputs and numerical row-generation choices."""
    result = {key: copy.deepcopy(spec[key]) for key in (
        'arm', 'Tg_K', 'P_Pa', 'loki_commit', 'axes', 'gas_properties',
        'state_properties', 'working_conditions', 'solver_options', 'screen', 'envelopes',
        'floors', 'tolerances', 'held_out', 'refinement')}
    result.update(_argon_validity_declarations(spec))
    if 'magnitude_tolerances' in spec:
        _validate_magnitude_tolerances(spec['magnitude_tolerances'])
        result['magnitude_tolerances'] = copy.deepcopy(spec['magnitude_tolerances'])
    if 'interpolation_accuracy' in spec:
        result['interpolation_accuracy'] = copy.deepcopy(spec['interpolation_accuracy'])
    if 'argon_4s_population' in spec:
        result['argon_4s_population'] = copy.deepcopy(spec['argon_4s_population'])
    if 'engine_state_map' in spec:
        result['engine_state_map'] = copy.deepcopy(spec['engine_state_map'])
    if 'vibrational_ladders' in spec:
        result['vibrational_ladders'] = copy.deepcopy(spec['vibrational_ladders'])
    if 'qualification_quantity_rules' in spec:
        result['qualification_quantity_rules'] = copy.deepcopy(
            spec['qualification_quantity_rules'])
    return result
