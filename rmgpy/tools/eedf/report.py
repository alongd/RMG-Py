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

"""Electron-power balance reports for latched EEDF budgets.

The report deliberately consumes the solver's accepted-state budget rather than
reconstructing powers from rates or reaction enthalpies.
"""
import argparse
from dataclasses import InitVar, dataclass, field
import html
import json
import math
from numbers import Real
from pathlib import Path
import re
from types import MappingProxyType


ARITHMETIC_RTOL = 1e-10
ARITHMETIC_ATOL_W = 1e-15
P1_ATOL_W = 1e-15
_VALIDATED_REPORT = object()

_REQUIRED_NUMERIC_FIELDS = {
    'P_abs', 'Q_elastic', 'Q_flow', 'Q_inelastic', 'Q_wall_electron',
    'Q_wall_ion', 'closure', 'dU_dt',
}
_REQUIRED_VECTOR_FIELDS = {
    'Q_heavy_particle_electron_by_reaction', 'Q_inelastic_channels',
    'Q_superelastic_channels',
}
_OPTIONAL_NUMERIC_FIELDS = {
    'A6a_power_scale', 'A6a_relative', 'A6a_steady_relative',
    'A6b_denominator', 'A6b_numerator', 'A6b_relative', 'A6b_tolerance',
    'EN_Td', 'Te', 'Te_eff', 'V', 'composition_energy_derivative', 'dNe_dt',
    'epsilon_k_eV', 'joule_power', 'mean_energy_eV', 'n_e', 'nu_ionisation',
    'nu_loss', 'nu_source', 'power_mismatch_fraction', 't', 'u',
}
_OPTIONAL_VECTOR_FIELDS = {
    'Q_inelastic_by_reaction', 'Q_inelastic_by_reaction_gross',
    'Q_superelastic_by_reaction',
}
_OPTIONAL_BOOLEAN_FIELDS = {
    'A6a_passed', 'A6a_steady_passed', 'A6b_passed', 'export_allowed',
    'qualification_allowed',
}
_OPTIONAL_TEXT_FIELDS = {
    'Te_eff_basis', 'discharge_state', 'scientific_status', 'wall_sheath_basis',
}
_QUALIFICATION_RECORD_FIELDS = {
    'A6b-runtime', 'A6b-source', 'T10 interpolation qualification',
}
_BUDGET_FIELDS = (
    _REQUIRED_NUMERIC_FIELDS | _REQUIRED_VECTOR_FIELDS | _OPTIONAL_NUMERIC_FIELDS
    | _OPTIONAL_VECTOR_FIELDS | _OPTIONAL_BOOLEAN_FIELDS | _OPTIONAL_TEXT_FIELDS
    | _QUALIFICATION_RECORD_FIELDS
    | {'Q_energy_only', 'artifact_blockers', 'electron_energy_disclosures'}
)
_CHANNEL_MAP_FIELDS = {
    'classification', 'cross_section', 'description', 'flux_group', 'kind',
    'mass_ratio', 'opb_eV', 'product_fraction', 'reaction', 'sigma_max_m2',
    'statistical_weight_ratio', 'target_fraction', 'threshold_eV',
}
_CHANNEL_FIELDS = _CHANNEL_MAP_FIELDS | {'core_reaction_index'}


class ReportRefusal(ValueError):
    """The persisted budget does not contain enough declared information."""


@dataclass(frozen=True)
class PowerReport:
    rows: tuple
    excitation_channels: tuple
    closure: float
    arithmetic_residual: float
    qualification: str
    p1_tolerance: float
    _validated: InitVar[object] = None
    unit: str = 'W (reactor inventory)'
    energy_only_groups: tuple = field(default_factory=tuple)
    channel_assignments: tuple = field(default_factory=tuple)
    show_excitation_channels: bool = False
    electron_energy_disclosure: str = ''
    a6b_runtime: object = None
    a6b_source: object = None
    t10_interpolation_qualification: object = None

    def __post_init__(self, _validated):
        if _validated is not _VALIDATED_REPORT:
            raise ReportRefusal('power report construction requires a validated budget')
        for name in ('rows', 'excitation_channels', 'energy_only_groups',
                     'channel_assignments'):
            records = tuple(MappingProxyType(dict(record))
                            for record in getattr(self, name))
            object.__setattr__(self, name, records)
        for name in ('a6b_runtime', 'a6b_source',
                     't10_interpolation_qualification'):
            record = getattr(self, name)
            if record is not None:
                object.__setattr__(self, name, MappingProxyType(dict(record)))

    @property
    def p1(self):
        absorbed = next(
            row['power_W'] for row in self.rows
            if row['name'] == 'absorbed electron power')
        return _p1_verdict(self.closure, absorbed, self.p1_tolerance)

    def as_dict(self):
        result = {
            'rows': [dict(row) for row in self.rows],
            'excitation_channels': [dict(item) for item in self.excitation_channels],
            'closure': self.closure,
            'arithmetic_residual': self.arithmetic_residual,
            'p1': self.p1,
            'qualification': self.qualification,
            'p1_tolerance': self.p1_tolerance,
            'unit': self.unit,
            'energy_only_groups': [dict(item) for item in self.energy_only_groups],
            'channel_assignments': [dict(item) for item in self.channel_assignments],
            'show_excitation_channels': self.show_excitation_channels,
            'electron_energy_disclosure': self.electron_energy_disclosure,
        }
        for key, name in (
                ('A6b-runtime', 'a6b_runtime'),
                ('A6b-source', 'a6b_source'),
                ('T10 interpolation qualification',
                 't10_interpolation_qualification')):
            record = getattr(self, name)
            if record is not None:
                result[key] = dict(record)
        return result


@dataclass(frozen=True)
class _ValidatedChannel:
    channel_id: int
    description: str
    kind: str
    classification: str
    flux_group: str
    group_id: int
    reaction_index: object
    core_reaction_position: object
    row: str


@dataclass(frozen=True)
class _ValidatedInputs:
    numbers: dict
    powers: list
    superelastic_powers: list
    heavy_powers: list
    energy_only: dict
    channels: list
    disclosures: list
    a6b_runtime: object
    a6b_source: object
    t10_interpolation_qualification: object
    qualification: str
    p1_rtol: float


def _number(value, label='power'):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ReportRefusal('invalid numeric ' + label)
    number = float(value)
    if not math.isfinite(number):
        raise ReportRefusal('non-finite ' + label)
    return number


def _finite_sum(values, label):
    total = 0.
    for value in values:
        total = _number(total + value, label)
    return total


def _finite_add(left, right, label):
    return _number(left + right, label)


def _finite_subtract(left, right, label):
    return _number(left - right, label)


def _finite_product(left, right, label):
    return _number(left * right, label)


def _relative_tolerance(rtol, scale, atol, label):
    return max(_finite_product(rtol, scale, label), atol)


def _p1_verdict(closure, absorbed, tolerance):
    allowed = _relative_tolerance(tolerance, absorbed, P1_ATOL_W,
                                  'P1 tolerance scale')
    return 'PASS' if abs(closure) <= allowed else 'FAIL'


def _single_line(value, label):
    if not isinstance(value, str) or not value or '\n' in value or '\r' in value:
        raise ReportRefusal('invalid ' + label)
    return value


def _qualification_record(value, name):
    """Validate one persisted qualification record as a closed structure."""
    numeric_fields = {
        'A6b-runtime': {
            'numerator', 'denominator', 'relative_error', 'tolerance',
        },
        'A6b-source': {'tolerance', 'maximum_relative_error'},
        'T10 interpolation qualification': set(),
    }
    integer_fields = {
        'A6b-runtime': set(),
        'A6b-source': set(),
        'T10 interpolation qualification': {'verdict_count'},
    }
    nullable_text_fields = {
        'A6b-runtime': set(),
        'A6b-source': set(),
        'T10 interpolation qualification': {'blocker'},
    }
    required = (
        {'check', 'criterion', 'passed'} | numeric_fields[name]
        | integer_fields[name] | nullable_text_fields[name]
    )
    if type(value) is not dict:
        raise ReportRefusal(name + ' must be a structured record')
    unknown = set(value) - required
    missing = required - set(value)
    if unknown or missing:
        detail = []
        if unknown:
            detail.append('unknown fields ' + ', '.join(sorted(unknown)))
        if missing:
            detail.append('missing fields ' + ', '.join(sorted(missing)))
        raise ReportRefusal(name + ' has ' + '; '.join(detail))
    if value['check'] != name:
        raise ReportRefusal(name + ' has an invalid check name')
    normalized = {
        'check': name,
        'criterion': _single_line(value['criterion'], name + ' criterion'),
    }
    for field_name in numeric_fields[name]:
        normalized[field_name] = _number(
            value[field_name], name + ' ' + field_name)
    for field_name in integer_fields[name]:
        field_value = value[field_name]
        if type(field_value) is not int or field_value < 0:
            raise ReportRefusal(name + ' has an invalid ' + field_name)
        normalized[field_name] = field_value
    for field_name in nullable_text_fields[name]:
        field_value = value[field_name]
        if field_value is not None:
            field_value = _single_line(field_value, name + ' ' + field_name)
        normalized[field_name] = field_value
    if type(value['passed']) is not bool:
        raise ReportRefusal(name + ' has an invalid passed verdict')
    normalized['passed'] = value['passed']
    return normalized


def _escape_markdown(value):
    escaped = html.escape(value, quote=True).replace('\\', '\\\\')
    return re.sub(r'([`*_[\]{}|])', r'\\\1', escaped)


def _vector(value, label):
    if not isinstance(value, list) or any(
            type(item) is not float or not math.isfinite(item) for item in value):
        raise ReportRefusal(label + ' must be a list of finite floats')
    return list(value)


def _reaction_id(reaction):
    if reaction is None:
        return None
    if not isinstance(reaction, dict):
        raise ReportRefusal('invalid mapped reaction')
    unknown = set(reaction) - {'library', 'index', 'repr', 'core_reaction_index', 'core_index'}
    if unknown:
        raise ReportRefusal('unknown reaction fields: ' + ', '.join(sorted(unknown)))
    reaction_index = reaction.get('index')
    if isinstance(reaction_index, bool) or not isinstance(reaction_index, int) \
            or reaction_index < 0:
        raise ReportRefusal('invalid reaction index')
    return reaction_index


def _declared_process(channel):
    """Return the manifest's declared process without parsing reaction repr text."""
    return _single_line(channel.get('description'), 'channel description')


def _is_electron_species(species):
    return species.lower() in ('e', 'e-')


_SPECIES = (
    r'[A-Za-z][A-Za-z0-9]*(?:[+-])?'
    r'(?:\((?:[^()\[\]\s<>]|\[[^\[\]\s<>]+\])*\))*\*?'
)
_STOICHIOMETRIC_TERM = re.compile(
    r'(?:(?P<multiplicity>[1-9][0-9]*) )?(?P<species>' + _SPECIES + r')\Z')


def _stoichiometric_side(text):
    """Resolve one process side using the supported collision grammar."""
    if not isinstance(text, str) or not text:
        return None
    raw_terms = []
    start = 0
    depth = 0
    for index, character in enumerate(text):
        if character == '(':
            depth += 1
        elif character == ')':
            depth -= 1
            if depth < 0:
                return None
        elif character == '+' and depth == 0 and text[index + 1:].strip():
            raw_terms.append(text[start:index].strip())
            start = index + 1
    if depth != 0:
        return None
    raw_terms.append(text[start:].strip())
    terms = []
    for raw in raw_terms:
        match = _STOICHIOMETRIC_TERM.fullmatch(raw)
        if match is None:
            return None
        terms.append((int(match.group('multiplicity') or '1'), match.group('species')))
    return terms


def _reaction_class(channel):
    """Classify the two cases absent from channel.kind from declared identity."""
    identity = _declared_process(channel)
    if ', ' in identity:
        equation, _display_kind = identity.rsplit(', ', 1)
    else:
        equation = identity
    kind = channel.get('kind')
    equation_match = re.fullmatch(r'(.+?) (<->|->) (.+)', equation)
    if equation_match is None:
        return None
    left_species = _stoichiometric_side(equation_match.group(1))
    right_species = _stoichiometric_side(equation_match.group(3))
    if left_species is None or right_species is None:
        return None
    left_electrons = sum(count for count, species in left_species
                         if _is_electron_species(species))
    right_electrons = sum(count for count, species in right_species
                          if _is_electron_species(species))
    left_heavy = [(count, species) for count, species in left_species
                  if not _is_electron_species(species)]
    right_heavy = [(count, species) for count, species in right_species
                   if not _is_electron_species(species)]
    left_heavy_count = sum(count for count, _ in left_heavy)
    right_heavy_count = sum(count for count, _ in right_heavy)
    if left_electrons != 1 or left_heavy_count != 1 or right_heavy_count == 0:
        return None
    if left_electrons > right_electrons:
        return 'attachment' if left_electrons - right_electrons == 1 \
            and kind == 'attachment' else None
    if right_electrons > left_electrons:
        # The declared collision category is structured metadata; punctuation
        # inside a species/state label is never interpreted as charge.
        if right_electrons - left_electrons != 1:
            return None
        if kind == 'ionization':
            return 'ionization'
        if kind == 'attachment':
            return 'detachment'
        return None
    if right_heavy_count > left_heavy_count:
        return 'dissociation'
    if left_electrons == right_electrons and left_heavy_count == right_heavy_count \
            and kind in ('elastic', 'excitation', 'vibrational', 'rotational'):
        return kind
    return None


def _row_for_channel(channel):
    """Return the single §17 row owned by one declared channel."""
    kind = channel.get('kind')
    if kind not in ('elastic', 'excitation', 'vibrational', 'rotational',
                    'ionization', 'attachment'):
        return None
    resolved = _reaction_class(channel)
    if resolved is None:
        return None
    direct = {
        'elastic': 'elastic',
        'vibrational': 'vibrational',
        'rotational': 'rotational',
        'ionization': 'ionization',
        'attachment': 'attachment',
    }
    if kind in direct:
        return direct[kind]
    return {'elastic': 'elastic', 'excitation': 'excitation',
            'vibrational': 'vibrational', 'rotational': 'rotational',
            'ionization': 'ionization', 'attachment': 'attachment',
            'dissociation': 'dissociation', 'detachment': 'attachment'}.get(resolved)


def _validate_channel_metadata(channel):
    for name in ('threshold_eV', 'sigma_max_m2', 'mass_ratio', 'opb_eV',
                 'statistical_weight_ratio'):
        if name in channel:
            value = _number(channel[name], 'channel ' + name)
            if value < 0.:
                raise ReportRefusal('invalid channel ' + name)
    for name in ('target_fraction', 'product_fraction'):
        if name in channel:
            value = _number(channel[name], 'channel ' + name)
            if value < 0. or value > 1.:
                raise ReportRefusal('invalid channel ' + name)
    if 'cross_section' not in channel:
        return
    section = channel['cross_section']
    if not isinstance(section, dict) or set(section) != {'energy_eV', 'sigma_m2'}:
        raise ReportRefusal('invalid channel cross section fields')
    if not isinstance(section['energy_eV'], list) \
            or not isinstance(section['sigma_m2'], list):
        raise ReportRefusal('invalid channel cross section vectors')
    energy = [
        _number(value, 'channel cross section energy')
        for value in section['energy_eV']
    ]
    sigma = [
        _number(value, 'channel cross section sigma')
        for value in section['sigma_m2']
    ]
    if len(energy) < 2 or len(sigma) != len(energy) \
            or any(value < 0. for value in energy + sigma) \
            or any(right <= left for left, right in zip(energy, energy[1:])):
        raise ReportRefusal('invalid channel cross section content')
    if 'sigma_max_m2' in channel \
            and float(channel['sigma_max_m2']).hex() != max(sigma).hex():
        raise ReportRefusal('channel sigma_max_m2 differs from cross section')


def _core_reaction_position(channel):
    positions = []
    if 'core_reaction_index' in channel:
        positions.append(channel['core_reaction_index'])
    reaction = channel['reaction']
    if isinstance(reaction, dict) and 'core_reaction_index' in reaction:
        positions.append(reaction['core_reaction_index'])
    if isinstance(reaction, dict) and 'core_index' in reaction:
        positions.append(reaction['core_index'])
    if not positions:
        return None
    for position in positions:
        if isinstance(position, bool) or not isinstance(position, int) or position < 0:
            raise ReportRefusal('invalid core reaction position: ' + str(position))
    if any(position != positions[0] for position in positions[1:]):
        raise ReportRefusal('conflicting core reaction positions')
    return positions[0]


def _qualification_text(qualification):
    if qualification is None:
        return 'not run'
    return 'unverified record'


def _validate_boundary(budget, channels, qualification, p1_rtol):
    """Validate and normalize every caller-controlled report input once."""
    if not isinstance(budget, dict) or 'Q_inelastic_channels' not in budget:
        raise ReportRefusal('legacy (non-EEDF) budget: Q_inelastic_channels is absent')
    unknown_budget = set(budget) - _BUDGET_FIELDS
    if unknown_budget:
        raise ReportRefusal('unknown budget fields: ' + ', '.join(sorted(unknown_budget)))
    required = _REQUIRED_NUMERIC_FIELDS | _REQUIRED_VECTOR_FIELDS | {'Q_energy_only'}
    missing = [key for key in required if key not in budget]
    if missing:
        raise ReportRefusal('EEDF budget is missing declared fields: ' + ', '.join(missing))
    numbers = {key: _number(budget[key], key) for key in _REQUIRED_NUMERIC_FIELDS}
    for key in _OPTIONAL_NUMERIC_FIELDS:
        if key in budget:
            _number(budget[key], key)
    for key in _OPTIONAL_VECTOR_FIELDS:
        if key in budget:
            _vector(budget[key], key)
    for key in _OPTIONAL_BOOLEAN_FIELDS:
        if key in budget and type(budget[key]) is not bool:
            raise ReportRefusal('invalid boolean ' + key)
    for key in _OPTIONAL_TEXT_FIELDS:
        if key in budget:
            _single_line(budget[key], key)
    if 'artifact_blockers' in budget:
        if not isinstance(budget['artifact_blockers'], list):
            raise ReportRefusal('invalid artifact blockers')
        for blocker in budget['artifact_blockers']:
            _single_line(blocker, 'artifact blocker')
    qualification_records = {
        name: _qualification_record(budget[name], name)
        if name in budget else None
        for name in _QUALIFICATION_RECORD_FIELDS
    }
    if numbers['P_abs'] <= 0.:
        raise ReportRefusal('P_abs must be positive')
    powers = _vector(budget['Q_inelastic_channels'], 'Q_inelastic_channels')
    superelastic = _vector(
        budget['Q_superelastic_channels'], 'Q_superelastic_channels')
    heavy = _vector(
        budget['Q_heavy_particle_electron_by_reaction'],
        'Q_heavy_particle_electron_by_reaction')
    if not powers:
        raise ReportRefusal('legacy (non-EEDF) budget: Q_inelastic_channels is empty')
    if not isinstance(channels, list):
        raise ReportRefusal('declared channel inventory is not a list')
    if len(channels) != len(powers):
        raise ReportRefusal('channel count does not match Q_inelastic_channels')
    if len(superelastic) != len(channels):
        raise ReportRefusal('channel count does not match Q_superelastic_channels')
    if not isinstance(budget['Q_energy_only'], dict):
        raise ReportRefusal('Q_energy_only is not a declared inventory')
    energy_only = {}
    for group, power in budget['Q_energy_only'].items():
        name = _single_line(group, 'Q_energy_only group')
        energy_only[name] = _number(power, 'Q_energy_only[' + name + ']')
    disclosures = budget['electron_energy_disclosures'] \
        if 'electron_energy_disclosures' in budget else None
    if disclosures is not None and not isinstance(disclosures, list):
        raise ReportRefusal('invalid electron-energy disclosure')
    if disclosures is not None:
        disclosures = [
            _single_line(item, 'electron-energy disclosure') for item in disclosures
        ]
    tolerance = _number(p1_rtol, 'P1 tolerance')
    if tolerance < 0.:
        raise ReportRefusal('invalid P1 tolerance')

    normalized_channels = []
    group_ids = {}
    supported_kinds = {
        'elastic', 'excitation', 'vibrational', 'rotational', 'ionization',
        'attachment',
    }
    for channel_id, channel in enumerate(channels):
        if not isinstance(channel, dict):
            raise ReportRefusal('invalid declared channel at index ' + str(channel_id))
        unknown_channel = set(channel) - _CHANNEL_FIELDS
        if unknown_channel:
            raise ReportRefusal('unknown channel fields: ' + ', '.join(sorted(unknown_channel)))
        missing_channel = [
            key for key in ('description', 'kind', 'classification', 'reaction')
            if key not in channel
        ]
        if missing_channel:
            raise ReportRefusal('channel is missing declared fields: '
                                + ', '.join(missing_channel))
        _validate_channel_metadata(channel)
        description = _single_line(channel['description'], 'channel description')
        flux_group = _single_line(channel['flux_group'], 'channel flux_group') \
            if 'flux_group' in channel else ''
        kind = channel['kind']
        if kind not in supported_kinds:
            raise ReportRefusal('unclassifiable channel: ' + description)
        classification = channel['classification']
        if classification not in ('A', 'B'):
            raise ReportRefusal('invalid channel classification: ' + description)
        reaction = channel['reaction']
        if classification == 'A' and reaction is None:
            raise ReportRefusal('mapped reaction unavailable: ' + description)
        if classification == 'B' and reaction is not None:
            raise ReportRefusal('channel classification disagrees with mapped reaction: '
                                + description)
        reaction_index = _reaction_id(reaction)
        if isinstance(reaction, dict) and 'repr' in reaction:
            _single_line(reaction['repr'], 'reaction repr')
        if isinstance(reaction, dict) and 'library' in reaction:
            _single_line(reaction['library'], 'reaction library')
        position = _core_reaction_position(channel)
        if position is not None and (
                isinstance(position, bool) or not isinstance(position, int)
                or position < 0 or position >= len(heavy)):
            raise ReportRefusal('invalid core reaction position: ' + str(position))
        group_id = -1
        if classification == 'B':
            if 'flux_group' in channel:
                flux_group = _single_line(
                    channel['flux_group'], 'channel flux_group')
            elif (all(record is not None
                      for record in qualification_records.values())
                  and description in energy_only):
                # The current runtime inventory uses its exact declared
                # description as the independently persisted B-group key.
                flux_group = description
            else:
                raise ReportRefusal(
                    'channel is missing declared flux_group: ' + description)
            if flux_group not in group_ids:
                group_ids[flux_group] = len(group_ids)
            group_id = group_ids[flux_group]
        row = _row_for_channel(channel)
        if row is None:
            raise ReportRefusal(
                'unclassifiable channel: ' + description
                + ' (unresolved product syntax)')
        normalized_channels.append(_ValidatedChannel(
            channel_id, description, kind, classification, flux_group, group_id,
            reaction_index, position, row,
        ))
    return _ValidatedInputs(
        numbers, powers, superelastic, heavy, energy_only, normalized_channels,
        disclosures, qualification_records['A6b-runtime'],
        qualification_records['A6b-source'],
        qualification_records['T10 interpolation qualification'],
        _qualification_text(qualification), tolerance,
    )


def power_table(budget, channels, *, qualification=None, p1_rtol=1e-3):
    """Return the §17 table from one accepted-state EEDF budget."""
    validated = _validate_boundary(budget, channels, qualification, p1_rtol)
    numbers = validated.numbers
    powers = validated.powers
    superelastic_powers = validated.superelastic_powers
    heavy = validated.heavy_powers
    declared_channels = validated.channels
    claimed_inelastic = numbers['Q_inelastic']
    gross_inelastic = _finite_sum(powers, 'gross Q_inelastic')
    superelastic_total = _finite_sum(
        superelastic_powers, 'summed Q_superelastic')
    heavy_total = _finite_sum(heavy, 'summed heavy transfer')
    reconstructed_inelastic = _finite_add(
        _finite_subtract(gross_inelastic, superelastic_total,
                         'net Q_inelastic'),
        heavy_total, 'reconstructed Q_inelastic')
    inelastic_scale = max(
        abs(claimed_inelastic),
        _finite_sum(
            (abs(value) for value in powers + superelastic_powers + heavy),
            'Q_inelastic magnitude scale'),
    )
    inelastic_difference = _finite_subtract(
        claimed_inelastic, reconstructed_inelastic, 'Q_inelastic difference')
    if abs(inelastic_difference) > _relative_tolerance(
            ARITHMETIC_RTOL, inelastic_scale, ARITHMETIC_ATOL_W,
            'Q_inelastic tolerance scale'):
        raise ReportRefusal('Q_inelastic disagrees with net channels plus heavy transfer')
    mapped = {}
    for channel in declared_channels:
        if channel.reaction_index is not None:
            mapped.setdefault(channel.reaction_index, []).append(channel)
    for rid, mapped_channels in mapped.items():
        positions = {channel.core_reaction_position for channel in mapped_channels}
        if None in positions:
            if any(value != 0. for value in heavy):
                raise ReportRefusal('core reaction position unavailable: ' + str(rid))
            continue
        for position in positions:
            if heavy[position] != 0.:
                raise ReportRefusal('double-counted core reaction: ' + str(position))
    values = {'elastic': numbers['Q_elastic'], 'excitation': 0.,
              'ionization': 0., 'attachment': 0., 'dissociation': 0.,
              'vibrational': 0., 'rotational': 0., 'superelastic': 0.}
    excitation = []
    channel_assignments = []
    for channel in declared_channels:
        value = powers[channel.channel_id]
        values[channel.row] = _finite_add(
            values[channel.row], value, 'aggregated row ' + channel.row)
        if channel.row == 'excitation':
            excitation.append({'channel': channel.description, 'power_W': value})
        values['superelastic'] = _finite_add(
            values['superelastic'], superelastic_powers[channel.channel_id],
            'aggregated row superelastic')
        channel_assignments.append({
            'channel': channel.description,
            'classification': channel.classification,
            'row': channel.row,
            'power_W': value,
        })
    energy_only_groups = []
    by_group = {}
    power_by_group = {}
    group_names = {}
    for channel in declared_channels:
        if channel.classification == 'B':
            group_id = channel.group_id
            group_names[group_id] = channel.flux_group
            by_group.setdefault(group_id, set()).add(channel.row)
            net_power = _finite_subtract(
                powers[channel.channel_id],
                superelastic_powers[channel.channel_id],
                'energy-only channel power')
            power_by_group[group_id] = _finite_add(
                power_by_group.get(group_id, 0.), net_power,
                'energy-only group power')
    expected_groups = set(group_names.values())
    supplied_groups = set(validated.energy_only)
    missing_groups = expected_groups - supplied_groups
    unexpected_groups = supplied_groups - expected_groups
    if missing_groups or unexpected_groups:
        details = []
        if missing_groups:
            details.append('missing ' + ', '.join(sorted(missing_groups)))
        if unexpected_groups:
            details.append('unexpected ' + ', '.join(sorted(unexpected_groups)))
        raise ReportRefusal('Q_energy_only inventory does not match declared B flux_group: '
                            + '; '.join(details))
    for group_id, group in group_names.items():
        rows_for_group = by_group[group_id]
        if len(rows_for_group) != 1:
            raise ReportRefusal('ambiguous energy-only group: ' + str(group))
        value = validated.energy_only[group]
        expected = power_by_group[group_id]
        group_scale = max(abs(value), abs(expected))
        group_difference = _finite_subtract(
            value, expected, 'Q_energy_only difference')
        if abs(group_difference) > _relative_tolerance(
                ARITHMETIC_RTOL, group_scale, ARITHMETIC_ATOL_W,
                'Q_energy_only tolerance scale'):
            raise ReportRefusal('Q_energy_only power disagrees with declared B channels: '
                                + group)
        energy_only_groups.append({'group': group, 'row': next(iter(rows_for_group)),
                                   'power_W': value})
    heavy_power = heavy_total
    disclosure_text = ('unavailable: persisted result has no verified reaction inventory '
                       'and electron_energies assignment coverage')
    if validated.disclosures:
        disclosure_text += ('; unverified caller disclosure: '
                            + '; '.join(validated.disclosures))
    rows = [{'name': 'absorbed electron power', 'power_W': numbers['P_abs']},
            {'name': 'elastic transfer', 'power_W': values['elastic']},
            {'name': 'excitation', 'power_W': values['excitation']},
            {'name': 'ionisation', 'power_W': values['ionization']},
            {'name': 'attachment/detachment', 'power_W': values['attachment']},
            {'name': 'dissociation', 'power_W': values['dissociation']},
            {'name': 'superelastic return', 'power_W': -values['superelastic']},
            {'name': 'heavy-particle electron transfer', 'power_W': heavy_power},
            {'name': 'wall loss', 'power_W': _finite_add(
                numbers['Q_wall_electron'], numbers['Q_wall_ion'],
                'aggregated row wall loss')},
            {'name': 'flow loss', 'power_W': numbers['Q_flow']},
            {'name': 'stored-energy change', 'power_W': numbers['dU_dt']}]
    rows += [{'name': 'vibrational', 'power_W': values['vibrational']},
             {'name': 'rotational', 'power_W': values['rotational']}]
    rows += [{'name': 'energy-only: ' + item['group'], 'power_W': item['power_W'],
              'included_in_balance': False, 'level': 1, 'row': item['row']}
             for item in energy_only_groups]
    for row in rows:
        _number(row['power_W'], 'aggregated row ' + row['name'])
    # Exclude the source row and compare the engine equation exactly.
    balance_rows = [
        row['power_W'] for row in rows[1:]
        if row.get('included_in_balance', True)
    ]
    reproduced = _finite_subtract(
        rows[0]['power_W'],
        _finite_sum(balance_rows, 'summed report rows'),
        'reproduced closure')
    residual = _finite_subtract(
        reproduced, numbers['closure'], 'reproduced residual')
    arithmetic_scale = max(
        abs(rows[0]['power_W']),
        _finite_sum((abs(value) for value in balance_rows),
                    'report magnitude scale'))
    arithmetic_tolerance = _relative_tolerance(
        ARITHMETIC_RTOL, arithmetic_scale, ARITHMETIC_ATOL_W,
        'report tolerance scale')
    if abs(residual) > arithmetic_tolerance:
        raise ReportRefusal('report arithmetic does not reproduce engine closure')
    rows.append({'name': 'numerical residual', 'power_W': numbers['closure'],
                 'included_in_balance': False})
    report = PowerReport(
        rows, excitation, numbers['closure'], residual,
        validated.qualification, validated.p1_rtol,
        _validated=_VALIDATED_REPORT,
        energy_only_groups=energy_only_groups,
        channel_assignments=channel_assignments,
        show_excitation_channels=(
            values['excitation'] >= _finite_product(
                .22, numbers['P_abs'], 'excitation display threshold')),
        electron_energy_disclosure=disclosure_text,
        a6b_runtime=validated.a6b_runtime,
        a6b_source=validated.a6b_source,
        t10_interpolation_qualification=(
            validated.t10_interpolation_qualification))
    report.p1
    return report


def render_markdown(report):
    if not isinstance(report, PowerReport):
        raise ReportRefusal('renderer requires a validated power report')
    lines = [f'| Quantity | Power — {report.unit} |', '|---|---:|']
    lines += [
        f"| {'&nbsp;&nbsp;↳ ' if row.get('level') else ''}{_escape_markdown(row['name'])} "
        f"| {row['power_W']:.12g} |"
        for row in report.rows
    ]
    lines += ['', f'Engine closure: `{report.closure:.12g}` W',
              f'Arithmetic residual: `{report.arithmetic_residual:.12g}` W',
              f'P1: **{report.p1}** (proposed, not frozen; tolerance {report.p1_tolerance:g})',
              'a-posteriori qualification: ' + report.qualification,
              'Electron-energy disclosure: ' + _escape_markdown(report.electron_energy_disclosure)]
    if report.show_excitation_channels and report.excitation_channels:
        lines += ['', 'Excitation channels:']
        lines += [f"- {_escape_markdown(item['channel'])}: {item['power_W']:.12g} W"
                  for item in report.excitation_channels]
    return '\n'.join(lines)


def _reconcile_channel_labels(data, channels):
    try:
        persisted = data['run_manifest']['a_posteriori_LoKI_B'][
            'direct_result']['channels']
    except (KeyError, TypeError):
        raise ReportRefusal('persisted result has no ordered channel labels')
    if not isinstance(persisted, list) or not isinstance(channels, list) \
            or len(persisted) != len(channels):
        raise ReportRefusal('persisted channel labels disagree with --channels')
    persisted = [
        _single_line(label, 'persisted channel label') for label in persisted
    ]
    supplied = [
        _single_line(channel.get('description'), 'channel description')
        if isinstance(channel, dict) else None
        for channel in channels
    ]
    if supplied != persisted:
        raise ReportRefusal('persisted channel labels disagree with --channels')


def _unique_object(pairs):
    """Refuse duplicate persisted names before record validation."""
    result = {}
    for name, value in pairs:
        if name in result:
            raise ReportRefusal('duplicate JSON name: ' + name)
        result[name] = value
    return result


def _main():
    parser = argparse.ArgumentParser()
    parser.add_argument('budget', type=Path)
    parser.add_argument('--channels', type=Path)
    args = parser.parse_args()
    data = json.loads(args.budget.read_text(), object_pairs_hook=_unique_object)
    if not isinstance(data, dict) or 'energy_budget' not in data:
        raise ReportRefusal('persisted result has no declared energy budget')
    budget = data['energy_budget']
    if not isinstance(budget, dict):
        raise ReportRefusal('persisted result has no declared energy budget')
    if 'Q_energy_only' not in budget or not isinstance(budget['Q_energy_only'], dict):
        raise ReportRefusal('persisted result has no declared Q_energy_only inventory')
    if args.channels is not None:
        channels = json.loads(args.channels.read_text())
    else:
        raise ReportRefusal('persisted result has no independent channel classification inventory; '
                            'provide --channels')
    _reconcile_channel_labels(data, channels)
    print(render_markdown(power_table(budget, channels)))


if __name__ == '__main__':
    _main()
