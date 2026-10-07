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
from dataclasses import dataclass, field
import argparse
import html
import json
import math
from numbers import Real
from pathlib import Path
import re


ARITHMETIC_RTOL = 1e-10
ARITHMETIC_ATOL_W = 1e-15
P1_ATOL_W = 1e-15

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
_BUDGET_FIELDS = (
    _REQUIRED_NUMERIC_FIELDS | _REQUIRED_VECTOR_FIELDS | _OPTIONAL_NUMERIC_FIELDS
    | _OPTIONAL_VECTOR_FIELDS | _OPTIONAL_BOOLEAN_FIELDS | _OPTIONAL_TEXT_FIELDS
    | {'Q_energy_only', 'artifact_blockers', 'electron_energy_disclosures'}
)
_CHANNEL_FIELDS = {
    'classification', 'cross_section', 'description', 'flux_group', 'kind',
    'mass_ratio', 'opb_eV', 'product_fraction', 'reaction', 'sigma_max_m2',
    'target_fraction', 'threshold_eV', 'core_reaction_index',
}


class ReportRefusal(ValueError):
    """The persisted budget does not contain enough declared information."""


@dataclass
class PowerReport:
    rows: list
    excitation_channels: list
    closure: float
    arithmetic_residual: float
    p1: str
    qualification: str
    p1_tolerance: float
    unit: str = 'W (reactor inventory)'
    energy_only_groups: list = field(default_factory=list)
    channel_assignments: list = field(default_factory=list)
    show_excitation_channels: bool = False
    electron_energy_disclosure: str = ''

    def as_dict(self):
        return {
            'rows': self.rows,
            'excitation_channels': self.excitation_channels,
            'closure': self.closure,
            'arithmetic_residual': self.arithmetic_residual,
            'p1': self.p1,
            'qualification': self.qualification,
            'p1_tolerance': self.p1_tolerance,
            'unit': self.unit,
            'energy_only_groups': self.energy_only_groups,
            'channel_assignments': self.channel_assignments,
            'show_excitation_channels': self.show_excitation_channels,
            'electron_energy_disclosure': self.electron_energy_disclosure,
        }


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
    qualification: str
    p1_rtol: float


def _number(value, label='power'):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ReportRefusal('invalid numeric ' + label)
    number = float(value)
    if not math.isfinite(number):
        raise ReportRefusal('non-finite ' + label)
    return number


def _single_line(value, label):
    if not isinstance(value, str) or not value or '\n' in value or '\r' in value:
        raise ReportRefusal('invalid ' + label)
    return value


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


_STOICHIOMETRIC_TERM = re.compile(
    r'(?:(?P<multiplicity>[1-9][0-9]*) )?(?P<species>[A-Za-z][^\s<>]*)\Z')


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
    try:
        equation, _display_kind = identity.rsplit(', ', 1)
    except ValueError:
        return None
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
    direct = {
        'elastic': 'elastic',
        'vibrational': 'vibrational',
        'rotational': 'rotational',
        'ionization': 'ionization',
    }
    if kind in direct:
        return direct[kind]
    special = _reaction_class(channel)
    if special is None:
        return None
    return {'elastic': 'elastic', 'excitation': 'excitation',
            'vibrational': 'vibrational', 'rotational': 'rotational',
            'ionization': 'ionization', 'attachment': 'attachment',
            'dissociation': 'dissociation', 'detachment': 'attachment'}.get(special)


def _core_reaction_position(channel):
    if 'core_reaction_index' in channel:
        return channel['core_reaction_index']
    reaction = channel['reaction']
    if isinstance(reaction, dict) and 'core_reaction_index' in reaction:
        return reaction['core_reaction_index']
    if isinstance(reaction, dict) and 'core_index' in reaction:
        return reaction['core_index']
    return None


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
        description = _single_line(channel['description'], 'channel description')
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
        flux_group = ''
        group_id = -1
        if classification == 'B':
            if 'flux_group' not in channel:
                raise ReportRefusal('channel is missing declared flux_group: ' + description)
            flux_group = _single_line(channel['flux_group'], 'channel flux_group')
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
        disclosures, _qualification_text(qualification), tolerance,
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
    reconstructed_inelastic = _number(
        sum(powers) - sum(superelastic_powers) + sum(heavy),
        'reconstructed Q_inelastic')
    inelastic_scale = max(
        abs(claimed_inelastic),
        sum(abs(value) for value in powers + superelastic_powers + heavy),
    )
    if abs(claimed_inelastic - reconstructed_inelastic) > max(
            ARITHMETIC_RTOL * inelastic_scale, ARITHMETIC_ATOL_W):
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
        values[channel.row] += value
        if channel.row == 'excitation':
            excitation.append({'channel': channel.description, 'power_W': value})
        values['superelastic'] += superelastic_powers[channel.channel_id]
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
            net_power = (powers[channel.channel_id]
                         - superelastic_powers[channel.channel_id])
            power_by_group[group_id] = power_by_group.get(group_id, 0.) + net_power
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
        if abs(value - expected) > max(
                ARITHMETIC_RTOL * group_scale, ARITHMETIC_ATOL_W):
            raise ReportRefusal('Q_energy_only power disagrees with declared B channels: '
                                + group)
        energy_only_groups.append({'group': group, 'row': next(iter(rows_for_group)),
                                   'power_W': value})
    heavy_power = sum(heavy)
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
            {'name': 'wall loss', 'power_W': numbers['Q_wall_electron'] + numbers['Q_wall_ion']},
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
    reproduced = rows[0]['power_W'] - sum(row['power_W'] for row in rows[1:]
                                          if row.get('included_in_balance', True))
    reproduced = _number(reproduced, 'reproduced closure')
    residual = _number(reproduced - numbers['closure'],
                       'reproduced residual')
    arithmetic_scale = max(abs(rows[0]['power_W']), sum(
        abs(row['power_W']) for row in rows[1:]
        if row.get('included_in_balance', True)))
    arithmetic_tolerance = max(ARITHMETIC_RTOL * arithmetic_scale, ARITHMETIC_ATOL_W)
    if abs(residual) > arithmetic_tolerance:
        raise ReportRefusal('report arithmetic does not reproduce engine closure')
    p1_tolerance_W = max(
        validated.p1_rtol * numbers['P_abs'], P1_ATOL_W)
    p1 = 'PASS' if abs(numbers['closure']) <= p1_tolerance_W else 'FAIL'
    rows.append({'name': 'numerical residual', 'power_W': numbers['closure'],
                 'included_in_balance': False})
    return PowerReport(rows, excitation, numbers['closure'], residual, p1,
                       validated.qualification, validated.p1_rtol,
                       energy_only_groups=energy_only_groups,
                       channel_assignments=channel_assignments,
                       show_excitation_channels=(
                           values['excitation'] >= .22 * numbers['P_abs']),
                       electron_energy_disclosure=disclosure_text)


def render_markdown(report):
    if isinstance(report, dict):
        try:
            report = PowerReport(**report)
        except TypeError:
            raise ReportRefusal('invalid serialized power report')
    if not isinstance(report, PowerReport):
        raise ReportRefusal('invalid power report')
    if report.unit != 'W (reactor inventory)':
        raise ReportRefusal('invalid report unit')
    if report.p1 not in ('PASS', 'FAIL'):
        raise ReportRefusal('invalid P1 status')
    if report.qualification not in ('not run', 'unverified record'):
        raise ReportRefusal('invalid qualification status')
    _number(report.closure, 'report closure')
    _number(report.arithmetic_residual, 'report arithmetic residual')
    tolerance = _number(report.p1_tolerance, 'report P1 tolerance')
    if tolerance < 0:
        raise ReportRefusal('invalid report P1 tolerance')
    _single_line(report.electron_energy_disclosure, 'electron-energy disclosure')
    if not isinstance(report.rows, list):
        raise ReportRefusal('invalid report rows')
    required_rows = {
        'absorbed electron power', 'elastic transfer', 'excitation', 'ionisation',
        'attachment/detachment', 'dissociation', 'superelastic return',
        'heavy-particle electron transfer', 'wall loss', 'flow loss',
        'stored-energy change', 'vibrational', 'rotational', 'numerical residual',
    }
    for row in report.rows:
        if not isinstance(row, dict) or 'name' not in row or 'power_W' not in row:
            raise ReportRefusal('invalid report row')
        _single_line(row['name'], 'report row name')
        _number(row['power_W'], 'report row power')
    rows_by_name = {
        name: [row for row in report.rows if row['name'] == name]
        for name in required_rows
    }
    if any(len(rows) != 1 for rows in rows_by_name.values()):
        raise ReportRefusal('serialized report does not contain complete report rows')
    if not isinstance(report.energy_only_groups, list):
        raise ReportRefusal('invalid energy-only group inventory')
    energy_only_names = set()
    for item in report.energy_only_groups:
        if not isinstance(item, dict) or set(item) != {'group', 'row', 'power_W'}:
            raise ReportRefusal('invalid energy-only group')
        group = _single_line(item['group'], 'energy-only group')
        if group in energy_only_names:
            raise ReportRefusal('duplicate energy-only group: ' + group)
        energy_only_names.add(group)
        if item['row'] not in {
                'elastic', 'excitation', 'ionization', 'attachment', 'dissociation',
                'vibrational', 'rotational'}:
            raise ReportRefusal('invalid energy-only group row')
        _number(item['power_W'], 'energy-only group power')
    expected_row_names = required_rows | {
        'energy-only: ' + group for group in energy_only_names
    }
    actual_row_names = [row['name'] for row in report.rows]
    if len(actual_row_names) != len(expected_row_names) \
            or set(actual_row_names) != expected_row_names:
        raise ReportRefusal('serialized report contains undeclared report rows')
    for row in report.rows:
        if row['name'].startswith('energy-only: '):
            allowed = {'name', 'power_W', 'included_in_balance', 'level', 'row'}
            if set(row) != allowed or row['included_in_balance'] is not False \
                    or type(row['level']) is not int or row['level'] != 1:
                raise ReportRefusal('invalid serialized energy-only row')
            group = row['name'][len('energy-only: '):]
            matches = [item for item in report.energy_only_groups if item['group'] == group]
            if len(matches) != 1 or row['row'] != matches[0]['row'] \
                    or row['power_W'] != matches[0]['power_W']:
                raise ReportRefusal('serialized energy-only row disagrees with group inventory')
        elif row['name'] == 'numerical residual':
            if set(row) != {'name', 'power_W', 'included_in_balance'} \
                    or row['included_in_balance'] is not False:
                raise ReportRefusal('invalid serialized numerical-residual row')
        elif set(row) != {'name', 'power_W'}:
            raise ReportRefusal('unknown serialized report-row fields')
    absorbed = rows_by_name['absorbed electron power'][0]['power_W']
    if absorbed <= 0.:
        raise ReportRefusal('report P_abs must be positive')
    expected_p1 = ('PASS' if abs(report.closure)
                   <= max(tolerance * absorbed, P1_ATOL_W) else 'FAIL')
    if report.p1 != expected_p1:
        raise ReportRefusal('stored P1 verdict disagrees with validated report numbers')
    if rows_by_name['numerical residual'][0]['power_W'] != report.closure:
        raise ReportRefusal('stored numerical residual disagrees with report closure')
    reproduced = absorbed - sum(
        row['power_W'] for row in report.rows
        if row['name'] != 'absorbed electron power'
        and row.get('included_in_balance', True))
    expected_residual = reproduced - report.closure
    render_scale = max(
        abs(absorbed),
        sum(abs(row['power_W']) for row in report.rows
            if row['name'] != 'absorbed electron power'
            and row.get('included_in_balance', True)),
    )
    if abs(expected_residual - report.arithmetic_residual) > max(
            ARITHMETIC_RTOL * render_scale, ARITHMETIC_ATOL_W):
        raise ReportRefusal('stored arithmetic residual disagrees with report rows')
    if not isinstance(report.show_excitation_channels, bool):
        raise ReportRefusal('invalid excitation-channel switch')
    if not isinstance(report.excitation_channels, list):
        raise ReportRefusal('invalid excitation-channel inventory')
    for item in report.excitation_channels:
        if not isinstance(item, dict) or set(item) != {'channel', 'power_W'}:
            raise ReportRefusal('invalid excitation channel')
        _single_line(item['channel'], 'excitation channel')
        _number(item['power_W'], 'excitation channel power')
    if not isinstance(report.channel_assignments, list):
        raise ReportRefusal('invalid channel-assignment inventory')
    for item in report.channel_assignments:
        if not isinstance(item, dict) or set(item) != {
                'channel', 'classification', 'row', 'power_W'}:
            raise ReportRefusal('invalid channel assignment')
        _single_line(item['channel'], 'assigned channel')
        if item['classification'] not in ('A', 'B') or item['row'] not in {
                'elastic', 'excitation', 'ionization', 'attachment', 'dissociation',
                'vibrational', 'rotational'}:
            raise ReportRefusal('invalid channel assignment metadata')
        _number(item['power_W'], 'assigned channel power')
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


def _main():
    parser = argparse.ArgumentParser()
    parser.add_argument('budget', type=Path)
    parser.add_argument('--channels', type=Path)
    args = parser.parse_args()
    data = json.loads(args.budget.read_text())
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
    print(render_markdown(power_table(budget, channels)))


if __name__ == '__main__':
    _main()
