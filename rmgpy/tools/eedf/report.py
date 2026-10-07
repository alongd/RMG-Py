"""Electron-power balance reports for latched EEDF budgets.

The report deliberately consumes the solver's accepted-state budget rather than
reconstructing powers from rates or reaction enthalpies.
"""
from dataclasses import dataclass, field
import argparse
import ast
import json
import math
from numbers import Real
from pathlib import Path
import re


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


def _vector(value, label):
    if isinstance(value, (str, bytes)):
        raise ReportRefusal(label + ' is not a vector')
    try:
        return list(value)
    except TypeError:
        raise ReportRefusal(label + ' is not a vector')


def _reaction_id(reaction):
    if reaction is None:
        return None
    if isinstance(reaction, dict):
        return reaction.get('id', reaction.get('index', reaction.get('repr')))
    return reaction


def _declared_process(channel):
    """Return the channel-map process, checking any mapped marker agrees."""
    process = _single_line(channel.get('description'), 'channel description')
    reaction = channel.get('reaction')
    if reaction is None:
        return process
    raw = reaction.get('repr') if isinstance(reaction, dict) else reaction
    try:
        call = ast.parse(raw, mode='eval').body
        if not isinstance(call, ast.Call) or not isinstance(call.func, ast.Name) \
                or call.func.id != 'EEDFChannel':
            raise ValueError
        fields = {item.arg: ast.literal_eval(item.value)
                  for item in call.keywords}
        mapped_process = fields.get('process')
    except (AttributeError, KeyError, SyntaxError, TypeError, ValueError):
        raise ReportRefusal('unresolved EEDFChannel identity: ' + repr(raw))
    if not isinstance(mapped_process, str) or not mapped_process:
        raise ReportRefusal('unresolved EEDFChannel identity: ' + repr(raw))
    if mapped_process != process:
        raise ReportRefusal('mapped EEDFChannel process differs from declared channel: ' + process)
    return process


def _is_electron_species(species):
    return species.lower() in ('e', 'e-')


_STOICHIOMETRIC_TERM = re.compile(
    r'(?:(?P<multiplicity>[1-9][0-9]*) )?(?P<species>[A-Za-z][^\s<>]*)\Z')


def _stoichiometric_side(text):
    """Resolve one process side using the supported collision grammar."""
    if not isinstance(text, str) or not text:
        return None
    terms = []
    for raw in text.split(' + '):
        match = _STOICHIOMETRIC_TERM.fullmatch(raw)
        if match is None:
            return None
        terms.append((int(match.group('multiplicity') or '1'), match.group('species')))
    return terms


def _reaction_class(channel):
    """Classify the two cases absent from channel.kind from declared identity."""
    identity = _declared_process(channel)
    try:
        equation, declared_kind = identity.rsplit(', ', 1)
    except ValueError:
        return None
    kind = channel.get('kind')
    if declared_kind.lower() != kind:
        return None
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
    special = _reaction_class(channel)
    if special is None:
        return None
    return {'elastic': 'elastic', 'excitation': 'excitation',
            'vibrational': 'vibrational', 'rotational': 'rotational',
            'ionization': 'ionization', 'attachment': 'attachment',
            'dissociation': 'dissociation', 'detachment': 'attachment'}.get(special)


def _core_reaction_position(channel):
    reaction = channel.get('reaction') or {}
    if 'core_reaction_index' in channel:
        return channel['core_reaction_index']
    if isinstance(reaction, dict):
        return reaction.get('core_reaction_index', reaction.get('core_index'))
    return None


def _qualification_text(qualification):
    if qualification is None:
        return 'not run'
    return 'unverified record'


def power_table(budget, channels, *, qualification=None, p1_rtol=1e-3):
    """Return the §17 table from one accepted-state EEDF budget."""
    if not isinstance(budget, dict) or 'Q_inelastic_channels' not in budget:
        raise ReportRefusal('legacy (non-EEDF) budget: Q_inelastic_channels is absent')
    required = ('P_abs', 'Q_inelastic', 'Q_elastic', 'Q_wall_electron',
                'Q_wall_ion', 'Q_flow', 'dU_dt', 'closure', 'Q_inelastic_channels',
                'Q_superelastic_channels', 'Q_heavy_particle_electron_by_reaction',
                'Q_energy_only')
    missing = [key for key in required if key not in budget]
    if missing:
        raise ReportRefusal('EEDF budget is missing declared fields: ' + ', '.join(missing))
    if channels is None:
        raise ReportRefusal('declared channel inventory is absent')
    channels = _vector(channels, 'declared channel inventory')
    powers = _vector(budget['Q_inelastic_channels'], 'Q_inelastic_channels')
    superelastic_powers = _vector(
        budget['Q_superelastic_channels'], 'Q_superelastic_channels')
    heavy = _vector(
        budget['Q_heavy_particle_electron_by_reaction'],
        'Q_heavy_particle_electron_by_reaction')
    if not powers:
        raise ReportRefusal('legacy (non-EEDF) budget: Q_inelastic_channels is empty')
    if len(channels) != len(powers):
        raise ReportRefusal('channel count does not match Q_inelastic_channels')
    if len(superelastic_powers) != len(channels):
        raise ReportRefusal('channel count does not match Q_superelastic_channels')
    for key in ('P_abs', 'Q_inelastic', 'Q_elastic', 'Q_wall_electron',
                'Q_wall_ion', 'Q_flow', 'dU_dt', 'closure'):
        _number(budget[key], key)
    powers = [_number(value, 'Q_inelastic_channels') for value in powers]
    superelastic_powers = [_number(value, 'Q_superelastic_channels')
                          for value in superelastic_powers]
    heavy = [_number(value, 'Q_heavy_particle_electron_by_reaction')
             for value in heavy]
    mapped = {}
    supported_kinds = ('elastic', 'excitation', 'vibrational', 'rotational',
                       'ionization', 'attachment')
    for index, channel in enumerate(channels):
        if not isinstance(channel, dict):
            raise ReportRefusal('invalid declared channel at index ' + str(index))
        missing_channel = [key for key in ('description', 'kind', 'classification', 'reaction')
                           if key not in channel]
        if missing_channel:
            raise ReportRefusal('channel is missing declared fields: '
                                + ', '.join(missing_channel))
        description = _single_line(channel['description'], 'channel description')
        if 'kind' not in channel or channel['kind'] not in supported_kinds:
            raise ReportRefusal('unclassifiable channel: ' + description)
        if channel['classification'] not in ('A', 'B', 'C', 'D'):
            raise ReportRefusal('invalid channel classification: ' + description)
        if channel['classification'] == 'B':
            if 'flux_group' not in channel:
                raise ReportRefusal('channel is missing declared flux_group: ' + description)
            _single_line(channel['flux_group'], 'channel flux_group')
        reaction = channel.get('reaction')
        if reaction is not None and channel.get('classification') != 'A':
            raise ReportRefusal('channel classification disagrees with mapped reaction: '
                                + description)
        if (reaction is None and channel.get('classification') == 'A'
                and any(abs(value) > 0 for value in heavy)):
            raise ReportRefusal('mapped reaction unavailable for nonzero heavy power: '
                                + description)
        if reaction is not None:
            _declared_process(channel)
        rid = _reaction_id(reaction)
        if rid is not None:
            mapped.setdefault(rid, []).append(index)
    for rid, indices in mapped.items():
        positions = {_core_reaction_position(channels[index]) for index in indices}
        if None in positions:
            if any(abs(value) > 0 for value in heavy):
                raise ReportRefusal('core reaction position unavailable: ' + str(rid))
            continue
        for position in positions:
            if (not isinstance(position, int) or position < 0 or position >= len(heavy)):
                raise ReportRefusal('invalid core reaction position: ' + str(position))
            if abs(heavy[position]) > 0:
                raise ReportRefusal('double-counted core reaction: ' + str(position))
    values = {'elastic': _number(budget['Q_elastic']), 'excitation': 0.,
              'ionization': 0., 'attachment': 0., 'dissociation': 0.,
              'vibrational': 0., 'rotational': 0., 'superelastic': 0.}
    excitation = []
    channel_assignments = []
    for index, (channel, power) in enumerate(zip(channels, powers)):
        value = power
        row = _row_for_channel(channel)
        if row is None:
            raise ReportRefusal('unclassifiable channel: ' + channel['description'])
        values[row] += value
        if row == 'excitation':
            excitation.append({'channel': channel['description'], 'power_W': value})
        values['superelastic'] += superelastic_powers[index]
        channel_assignments.append({'channel': channel['description'],
                                    'classification': channel['classification'],
                                    'row': row, 'power_W': value})
    energy_only = budget['Q_energy_only']
    if not isinstance(energy_only, dict):
        raise ReportRefusal('Q_energy_only is not a declared inventory')
    energy_only_groups = []
    by_group = {}
    power_by_group = {}
    for index, channel in enumerate(channels):
        if channel.get('classification') == 'B':
            group = channel['flux_group']
            row = _row_for_channel(channel)
            if row is None:
                raise ReportRefusal('unclassifiable energy-only group: ' + str(group))
            by_group.setdefault(group, set()).add(row)
            power_by_group[group] = power_by_group.get(group, 0.) + powers[index]
    for group in energy_only:
        _single_line(group, 'Q_energy_only group')
    missing_groups = set(by_group) - set(energy_only)
    unexpected_groups = set(energy_only) - set(by_group)
    if missing_groups or unexpected_groups:
        details = []
        if missing_groups:
            details.append('missing ' + ', '.join(sorted(missing_groups)))
        if unexpected_groups:
            details.append('unexpected ' + ', '.join(sorted(unexpected_groups)))
        raise ReportRefusal('Q_energy_only inventory does not match declared B flux_group: '
                            + '; '.join(details))
    for group, power in energy_only.items():
        rows_for_group = by_group.get(group)
        if not rows_for_group:
            raise ReportRefusal('unclassifiable energy-only group: ' + str(group))
        if len(rows_for_group) != 1:
            raise ReportRefusal('ambiguous energy-only group: ' + str(group))
        value = _number(power, 'Q_energy_only[' + str(group) + ']')
        if not math.isclose(value, power_by_group[group], rel_tol=1e-12, abs_tol=1e-12):
            raise ReportRefusal('Q_energy_only power disagrees with declared B channels: '
                                + group)
        energy_only_groups.append({'group': group, 'row': next(iter(rows_for_group)),
                                   'power_W': value})
    heavy_power = sum(heavy)
    disclosure = budget.get('electron_energy_disclosures')
    disclosure_text = ('unavailable: persisted result has no verified reaction inventory '
                       'and electron_energies assignment coverage')
    if disclosure is not None and (not isinstance(disclosure, list)
                                   or not all(isinstance(item, str) for item in disclosure)):
        raise ReportRefusal('invalid electron-energy disclosure')
    if disclosure:
        for item in disclosure:
            _single_line(item, 'electron-energy disclosure')
        disclosure_text += '; unverified caller disclosure: ' + '; '.join(disclosure)
    rows = [{'name': 'absorbed electron power', 'power_W': _number(budget['P_abs'])},
            {'name': 'elastic transfer', 'power_W': values['elastic']},
            {'name': 'excitation', 'power_W': values['excitation']},
            {'name': 'ionisation', 'power_W': values['ionization']},
            {'name': 'attachment/detachment', 'power_W': values['attachment']},
            {'name': 'dissociation', 'power_W': values['dissociation']},
            {'name': 'superelastic return', 'power_W': -values['superelastic']},
            {'name': 'heavy-particle electron transfer', 'power_W': heavy_power},
            {'name': 'wall loss', 'power_W': _number(budget['Q_wall_electron']) + _number(budget['Q_wall_ion'])},
            {'name': 'flow loss', 'power_W': _number(budget['Q_flow'])},
            {'name': 'stored-energy change', 'power_W': _number(budget['dU_dt'])}]
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
    residual = _number(reproduced - _number(budget['closure'], 'closure'),
                       'reproduced residual')
    if abs(residual) > max(1e-9, abs(_number(budget['P_abs'])) * 1e-10):
        raise ReportRefusal('report arithmetic does not reproduce engine closure')
    p1_rtol = _number(p1_rtol, 'P1 tolerance')
    if p1_rtol < 0:
        raise ReportRefusal('invalid P1 tolerance')
    p1 = 'PASS' if abs(_number(budget['closure'])) <= p1_rtol * abs(_number(budget['P_abs'])) else 'FAIL'
    rows.append({'name': 'numerical residual', 'power_W': _number(budget['closure'], 'closure'),
                 'included_in_balance': False})
    return PowerReport(rows, excitation, _number(budget['closure'], 'closure'), residual, p1,
                       _qualification_text(qualification), p1_rtol,
                       energy_only_groups=energy_only_groups,
                       channel_assignments=channel_assignments,
                       show_excitation_channels=(
                           values['excitation'] >= .22 * _number(budget['P_abs'], 'P_abs')),
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
    for row in report.rows:
        if not isinstance(row, dict) or 'name' not in row or 'power_W' not in row:
            raise ReportRefusal('invalid report row')
        _single_line(row['name'], 'report row name')
        _number(row['power_W'], 'report row power')
    if not isinstance(report.show_excitation_channels, bool):
        raise ReportRefusal('invalid excitation-channel switch')
    if not isinstance(report.excitation_channels, list):
        raise ReportRefusal('invalid excitation-channel inventory')
    for item in report.excitation_channels:
        if not isinstance(item, dict) or 'channel' not in item or 'power_W' not in item:
            raise ReportRefusal('invalid excitation channel')
        _single_line(item['channel'], 'excitation channel')
        _number(item['power_W'], 'excitation channel power')
    lines = [f'| Quantity | Power — {report.unit} |', '|---|---:|']
    lines += [f"| {'&nbsp;&nbsp;↳ ' if row.get('level') else ''}{row['name']} | {row['power_W']:.12g} |"
              for row in report.rows]
    lines += ['', f'Engine closure: `{report.closure:.12g}` W',
              f'Arithmetic residual: `{report.arithmetic_residual:.12g}` W',
              f'P1: **{report.p1}** (proposed, not frozen; tolerance {report.p1_tolerance:g})',
              'a-posteriori qualification: ' + report.qualification,
              'Electron-energy disclosure: ' + report.electron_energy_disclosure]
    if report.show_excitation_channels and report.excitation_channels:
        lines += ['', 'Excitation channels:']
        lines += [f"- {item['channel']}: {item['power_W']:.12g} W" for item in report.excitation_channels]
    return '\n'.join(lines)


def _main():
    parser = argparse.ArgumentParser()
    parser.add_argument('budget', type=Path)
    parser.add_argument('--channels', type=Path)
    args = parser.parse_args()
    data = json.loads(args.budget.read_text())
    budget = data.get('energy_budget', data)
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
