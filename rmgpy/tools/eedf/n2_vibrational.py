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

"""Offline N2(X,v) ladder generation with the D6 top-leak closure.

Inputs are local, cited, SHA256-pinned records. Cross-section payloads are
never embedded in output: each collision references its pinned source block.
Qualified hierarchical population handling belongs to the coupling layer.
"""
import argparse
import hashlib
from io import BytesIO
import json
import math
from pathlib import Path
import re
import tempfile
import zipfile

N2_GRAPH = '1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}\n'


class LadderError(ValueError):
    """A source, coverage, or D6 requirement is unsatisfied."""


class MissingPairsError(LadderError):
    """A declared source domain has uncovered upward pairs."""

    def __init__(self, pairs):
        self.pairs = tuple(pairs)
        super().__init__('missing N2 source pairs: ' + ', '.join(
            f'({i},{j})' for i, j in self.pairs))


def _number(value, name, minimum=0):
    if type(value) not in (int, float) or not math.isfinite(value) or value < minimum:
        raise LadderError(f'{name} must be finite and >= {minimum}')
    return value


def _integer(value, name, minimum, maximum):
    if type(value) is not int or not minimum <= value <= maximum:
        raise LadderError(f'{name} must be an integer in [{minimum}, {maximum}]')
    return value


def _pin(pin):
    """Read and hash the same snapshot; do not re-open mutable source paths."""
    if pin.get('provenance') != 'sourced':
        raise LadderError('only sourced inputs are supported; estimated is reserved')
    if not isinstance(pin.get('citation'), str) or not pin['citation'].strip():
        raise LadderError('sourced input requires a citation')
    if not re.fullmatch('[a-f0-9]{64}', pin.get('sha256', '')):
        raise LadderError('sourced input requires a SHA256')
    path = Path(pin['path'])
    if not path.is_absolute():
        raise LadderError('source paths must be absolute')
    snapshot = path.read_bytes()
    if hashlib.sha256(snapshot).hexdigest() != pin['sha256']:
        raise LadderError('source hash mismatch: ' + str(path))
    member = pin.get('member')
    if member is None:
        return snapshot
    if not isinstance(member, str) or not member:
        raise LadderError('archive member must be a nonempty path')
    try:
        with zipfile.ZipFile(BytesIO(snapshot)) as archive:
            return archive.read(member)
    except (KeyError, zipfile.BadZipFile) as exc:
        raise LadderError(f'unreadable archive member {member!r}: {path}') from exc


def read_cross_sections(pin):
    """Parse upward/diagonal LXCat Laporta N2 blocks, returning references.

    The one-based block index counts EXCITATION blocks. Unsupported identities,
    downward records, duplicate pairs and malformed grids refuse. Statistical
    weights for detailed balance come from the explicit spectroscopy record.
    """
    lines = _pin(pin).decode('utf-8').splitlines()
    result = []
    position = 0
    while position < len(lines):
        if lines[position].strip() != 'EXCITATION':
            position += 1
            continue
        position += 1
        match = re.fullmatch(r'N2\(X,v=(\d+)\) -> N2\(X,v=(\d+)\)', lines[position].strip())
        if not match:
            raise LadderError('unsupported N2 source identity')
        lower, upper = map(int, match.groups())
        position += 1
        threshold_fields = lines[position].split()
        if len(threshold_fields) != 1:
            raise LadderError('expected one sourced excitation threshold')
        threshold = _number(float(threshold_fields[0]), 'source threshold')
        if upper < lower or (lower == upper and threshold != 0):
            raise LadderError('only upward and zero-energy diagonal source pairs are supported')
        position += 1
        while position < len(lines) and not re.fullmatch('-{5,}', lines[position].strip()):
            if lines[position].strip() == 'EXCITATION':
                raise LadderError('missing cross-section grid')
            position += 1
        position += 1
        count, last_energy = 0, -1.
        while position < len(lines) and not re.fullmatch('-{5,}', lines[position].strip()):
            fields = lines[position].split()
            if len(fields) != 2:
                raise LadderError('malformed cross-section row')
            energy, sigma = map(float, fields)
            _number(energy, 'cross-section energy')
            _number(sigma, 'cross section')
            if energy <= last_energy:
                raise LadderError('cross-section energies must strictly increase')
            last_energy = energy
            count += 1
            position += 1
        if position >= len(lines) or count < 2:
            raise LadderError('incomplete cross-section grid')
        result.append({'initial_v': lower, 'final_v': upper,
                       'threshold_eV': threshold, 'source_block': len(result) + 1,
                       'source': dict(pin)})
        position += 1
    if not result:
        raise LadderError('no N2 cross-section blocks')
    return result


def require_complete_pairs(sections, upper_v=None, *, initial_levels=None, final_v=None):
    """Require every upward pair for covered initial levels through ``final_v``.

    The legacy ``upper_v`` form means initial and final levels are both 0..upper_v.
    """
    if initial_levels is None and final_v is None:
        _integer(upper_v, 'source upper_v', 1, 1000)
        initial_levels = range(upper_v + 1)
        final_v = upper_v
    elif upper_v is not None or initial_levels is None or final_v is None:
        raise LadderError('provide either upper_v or both initial_levels and final_v')
    if (not isinstance(initial_levels, (list, tuple, range))
            or isinstance(initial_levels, (str, bytes)) or not initial_levels):
        raise LadderError('initial_levels must be a nonempty ordered sequence')
    checked_initials = []
    for initial_v in initial_levels:
        checked_initials.append(_integer(initial_v, 'initial level', 0, 1000))
    if checked_initials != sorted(set(checked_initials)):
        raise LadderError('initial_levels must be unique and increasing')
    _integer(final_v, 'final_v', 1, 1000)
    if checked_initials[-1] > final_v:
        raise LadderError('initial_levels cannot exceed final_v')
    pairs = set()
    for entry in sections:
        pair = entry['initial_v'], entry['final_v']
        if pair[0] not in checked_initials or not 0 <= pair[0] <= pair[1] <= final_v:
            raise LadderError(f'source pair outside declared domain: {pair}')
        if pair in pairs:
            raise LadderError(f'duplicate source pair: {pair}')
        pairs.add(pair)
    missing = [(i, j) for i in checked_initials
               for j in range(i + 1, final_v + 1) if (i, j) not in pairs]
    if missing:
        raise MissingPairsError(missing)


def _policy(policy, allow_unset=False):
    if policy.get('set_by') != 'owner':
        raise LadderError('D6 thresholds must be explicitly owner-set')
    for name in ('top_fraction_limit', 'top_leak_share_limit'):
        if allow_unset and policy.get(name) is None:
            continue
        value = _number(policy.get(name), name)
        if value > 1:
            raise LadderError(name + ' must be <= 1')
    _integer(policy.get('K_step'), 'K_step', 1, 30)


def assess_truncation(K, populations, forward_vibrational_power, top_leak_power, policy):
    """Pure D6 decision: return accept/regenerate, or refuse at K=30.

    Populations are nonnegative amounts for v=0..K, normalized here. Powers
    share one caller-chosen unit; the denominator is gross forward vibrational
    electron power, including above-K excitations. No gas heating is implied.
    Unset thresholds always refuse. The coupling layer owns execution/restart.
    """
    _integer(K, 'K', 1, 30)
    _policy(policy)
    if len(populations) != K + 1:
        raise LadderError('populations must cover exactly v=0..K')
    total = sum(_number(n, 'population') for n in populations)
    if not math.isfinite(total) or total <= 0:
        raise LadderError('total population must be positive and finite')
    forward = _number(forward_vibrational_power, 'forward vibrational power')
    leak = _number(top_leak_power, 'top-leak power')
    if leak > forward:
        raise LadderError('top-leak power exceeds gross forward vibrational power')
    fraction = populations[K] / total
    share = leak / forward if forward else 0.
    reasons = []
    if fraction > policy['top_fraction_limit']:
        reasons.append('top_fraction_limit')
    if share > policy['top_leak_share_limit']:
        reasons.append('top_leak_share_limit')
    if reasons and K == 30:
        raise LadderError('D6 ceiling K=30 exceeded: ' + ', '.join(reasons))
    return {'action': 'regenerate' if reasons else 'accept',
            'K': min(30, K + policy['K_step']) if reasons else K,
            'top_fraction': fraction, 'top_leak_share': share, 'exceeded': reasons}


def _adjacency(v):
    return f'vibrationallevel {v}\n' + N2_GRAPH


def _validate_thermo(data, levels, K):
    """Load the pinned library snapshot and require exact X-state fixed-v entries."""
    from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary
    from rmgpy.molecule.molecule import Molecule
    from rmgpy.species import Species
    from rmgpy.thermo import NASA, ThermoData, Wilhoit

    database = ThermoDatabase()
    library = ThermoLibrary(label='N2Vibrational')
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / 'thermo.py'
        path.write_bytes(data)
        library.load(str(path), database.local_context, {})
    for v in range(K + 1):
        label = levels[v]['thermo_entry']
        entry = library.entries.get(label)
        expected = Molecule().from_adjacency_list(_adjacency(v))
        if entry is None or entry.data is None or not entry.item.is_isomorphic(expected):
            raise LadderError(f'missing exact library thermo for N2(X,v={v}): {label}')
        if not isinstance(entry.data, (NASA, Wilhoit, ThermoData)):
            raise LadderError(f'unsupported library thermo for N2(X,v={v}): {label}')
        # Exercise the same exact lookup as the emitted library's consumer,
        # without entering get_thermo_data's estimation fallback.
        try:
            match = database.get_thermo_data_from_library(
                Species(molecule=[expected]), library)
        except Exception as exc:
            raise LadderError(f'unusable library thermo for N2(X,v={v}): {label}') from exc
        if (match is None or match[2] is not entry
                or not isinstance(match[0], (NASA, Wilhoit, ThermoData))):
            raise LadderError(f'missing exact library thermo for N2(X,v={v}): {label}')


def generate(manifest, K=10, production=True):
    """Return deterministic artifact bytes after validating all input snapshots.

    ``manifest`` is a decoded JSON object. ``production=False`` permits unset
    owner thresholds for a draft only; source coverage and thermo still refuse.
    Output consists of manifest.json, input.py and an unchanged thermo.py library
    snapshot. The deck is a loadable species/manifold fragment, not a reactor run.
    """
    _integer(K, 'K', 1, 30)
    if manifest.get('schema') != 'n2-vibrational-source-v1':
        raise LadderError('unsupported N2 source schema')
    if not isinstance(manifest.get('D6'), dict):
        raise LadderError('D6 policy must be an object')
    _policy(manifest['D6'], allow_unset=not production)
    cross_section_pins = manifest.get('cross_sections')
    if not isinstance(cross_section_pins, list) or not cross_section_pins:
        raise LadderError('cross_sections must be a nonempty list of source pins')
    if any(not isinstance(pin, dict) for pin in cross_section_pins):
        raise LadderError('each cross-section source pin must be an object')
    level_pin = manifest['level_source']
    spectroscopy = json.loads(_pin(level_pin))
    if (spectroscopy.get('molecule') != '14N2' or spectroscopy.get('electronic_state') != 'X1Sigma_g+'
            or spectroscopy.get('energy_unit') != 'eV' or spectroscopy.get('energy_zero') != 'v=0'):
        raise LadderError('require explicit 14N2 X1Sigma_g+ term energies in eV relative to v=0')
    levels = spectroscopy['levels']
    if not isinstance(levels, list) or len(levels) < K + 1:
        raise LadderError('level source does not cover K')
    source_pins = manifest.get('cross_sections')
    if not isinstance(source_pins, list) or not source_pins or any(not isinstance(pin, dict) for pin in source_pins):
        raise LadderError('cross_sections must be a nonempty list of source pins')
    has_domain = 'source_domain' in manifest
    has_assignments = 'source_assignments' in manifest
    has_pin_ids = [isinstance(pin.get('source_id'), str) and bool(pin['source_id']) for pin in source_pins]
    legacy = not has_domain and not has_assignments and not any('source_id' in pin for pin in source_pins)
    modern = has_domain and has_assignments and all(has_pin_ids)
    if not (legacy or modern):
        raise LadderError('source_domain, source_assignments, and source_id must use one complete manifest form')
    source_domain = manifest.get('source_domain')
    if legacy:
        source_domain = {'initial_levels': list(range(len(levels))), 'final_v': len(levels) - 1}
    if not isinstance(source_domain, dict):
        raise LadderError('source_domain must declare initial_levels and final_v')
    initial_levels = source_domain.get('initial_levels')
    final_v = source_domain.get('final_v')
    if (not isinstance(initial_levels, list) or not initial_levels
            or any(type(v) is not int for v in initial_levels)):
        raise LadderError('source_domain.initial_levels must be a nonempty integer list')
    _integer(final_v, 'source_domain.final_v', 1, 1000)
    if K not in initial_levels:
        raise LadderError(f'requested K={K} requires uncovered initial v={K}')
    if len(levels) != final_v + 1:
        raise LadderError('level source must cover the declared final_v domain')
    previous = -1.
    for v, level in enumerate(levels):
        if set(level) != {'v', 'energy_eV', 'statistical_weight', 'thermo_entry'}:
            raise LadderError('level records require exactly v, energy_eV, statistical_weight, thermo_entry')
        if type(level['v']) is not int or level['v'] != v:
            raise LadderError('levels must be contiguous and ordered from v=0')
        energy = _number(level['energy_eV'], 'level energy')
        weight = _number(level['statistical_weight'], 'statistical weight')
        if weight == 0 or energy <= previous or (v == 0 and energy != 0):
            raise LadderError('invalid term energies or statistical weights')
        if v <= K and (not isinstance(level.get('thermo_entry'), str) or not level['thermo_entry']):
            raise LadderError(f'level v={v} requires a thermo entry label for requested K={K}')
        if v > K and level.get('thermo_entry') is not None and not isinstance(level['thermo_entry'], str):
            raise LadderError(f'invalid thermo entry label for source-only level v={v}')
        previous = energy
    sections = [section for pin in cross_section_pins for section in read_cross_sections(pin)]
    assignments = manifest.get('source_assignments')
    if legacy:
        assignments = {str(v): None for v in initial_levels}
    else:
        if not isinstance(assignments, dict):
            raise LadderError('source_assignments must map initial levels to source ids')
        source_ids = {pin.get('source_id'): pin for pin in cross_section_pins}
        if (len(source_ids) != len(cross_section_pins)
                or any(not isinstance(source_id, str) or not source_id for source_id in source_ids)):
            raise LadderError('every cross-section pin requires a unique source_id')
        if (any(not isinstance(source_id, str) for source_id in assignments)
                or any(not isinstance(source_id, str) for source_id in assignments.values())):
            raise LadderError('source_assignments keys and values must be strings')
        if any(source_id not in source_ids for source_id in assignments.values()):
            raise LadderError('source assignment refers to an unknown source id')
        if set(assignments) != {str(v) for v in initial_levels}:
            raise LadderError('source_assignments must exactly match source_domain.initial_levels')
        sections = [section for section in sections
                    if assignments.get(str(section['initial_v'])) == section['source'].get('source_id')]
    required_initials = list(range(K + 1))
    uncovered = [v for v in required_initials if str(v) not in assignments]
    if uncovered:
        raise LadderError(f'requested K={K} requires uncovered initial v={uncovered[0]}')
    require_complete_pairs(sections, initial_levels=initial_levels, final_v=final_v)
    tolerance = _number(manifest['energy_tolerance_eV'], 'energy tolerance')
    for section in sections:
        i, j = section['initial_v'], section['final_v']
        gap = levels[j]['energy_eV'] - levels[i]['energy_eV']
        if abs(section['threshold_eV'] - gap) > tolerance:
            raise LadderError(f'sourced threshold disagrees with term energies for ({i},{j})')
    thermo_data = _pin(manifest['thermo_library'])
    _validate_thermo(thermo_data, levels, K)
    retained = [{'v': v, 'label': 'N2' if v == 0 else f'N2v{v}',
                 'adjacency_list': _adjacency(v), **levels[v],
                 'provenance': dict(level_pin),
                 'thermo': {**manifest['thermo_library'], 'entry': levels[v]['thermo_entry']}}
                for v in range(K + 1)]
    audit, channels, ledger = [], [], []
    for pair_index, section in enumerate(sorted(sections, key=lambda x: (x['initial_v'], x['final_v']))):
        i, j = section['initial_v'], section['final_v']
        pair_id = f'N2:v{i}-v{j}'
        record = dict(section, source_pair_index=pair_index, pair_id=pair_id)
        if i > K:
            status = 'unresolved_initial_above_K'
        elif j > K:
            status = 'redirected_to_top'
        elif i == j:
            status = 'retained_diagonal'
        else:
            status = 'retained'
        audit.append(dict(record, disposition=status))
        # No population or reverse is invented for an unresolved initial level.
        if i > K:
            continue
        product_v = min(j, K)
        leak = levels[j]['energy_eV'] - levels[product_v]['energy_eV']
        retained_gap = levels[product_v]['energy_eV'] - levels[i]['energy_eV']
        entry = dict(record, kind='vibrational', direction='forward', reversible=False,
                     initial_v=i, final_v=product_v, source_final_v=j,
                     reactants={'electron': 1, retained[i]['label']: 1},
                     products={'electron': 1, retained[product_v]['label']: 1},
                     electron_energy_loss_eV=section['threshold_eV'],
                     retained_energy_gain_eV=retained_gap, top_leak_eV=leak,
                     source_energy_residual_eV=section['threshold_eV'] - retained_gap - leak,
                     rate_binding={'source_pair_index': pair_index, 'column': 'k_ine'})
        channels.append(entry)
        if j > K:
            ledger.append({'name': 'top-leak:' + pair_id, 'pair_id': pair_id,
                           'source_pair_index': pair_index, 'energy_eV_per_event': leak,
                           'molecular_gain_eV': leak, 'electron_loss_eV': -leak,
                           'gas_heating_eV': 0., 'provenance': dict(section['source'])})
        elif i != j:
            channels.append(dict(entry, direction='superelastic', initial_v=j, final_v=i,
                                 reactants=entry['products'], products=entry['reactants'],
                                 electron_energy_loss_eV=-section['threshold_eV'],
                                 retained_energy_gain_eV=-retained_gap, top_leak_eV=0.,
                                 source_energy_residual_eV=-entry['source_energy_residual_eV'],
                                 rate_binding={'source_pair_index': pair_index, 'column': 'k_sup'},
                                 detailed_balance={'source_pair_index': pair_index,
                                                  'weight_ratio_initial_over_final':
                                                  levels[i]['statistical_weight'] / levels[j]['statistical_weight']}))
    ready = all(manifest['D6'].get(key) is not None
                for key in ('top_fraction_limit', 'top_leak_share_limit'))
    output = {'schema': 'n2-vibrational-ladder-v1', 'K': K, 'source_upper_v': final_v,
              'source_complete': True, 'production_manifest': production and ready,
              'runtime_ready': False,
              'runtime_requirements': (['qualified EEDF hierarchical N2(X,v) selectors',
                                        'coupling executes D6 regeneration before accepting a run'] if legacy else
                                       ['coupling automatically executes D6 regeneration before accepting a run',
                                        'physical Ar+N2 arm qualification requires complete channel maps and class-A owners in their configured per-gas libraries']),
              'D6': dict(manifest['D6'], maximum_K=30, top_v=K,
                         evaluator='rmgpy.tools.eedf.n2_vibrational.assess_truncation',
                         power_denominator='gross forward vibrational electron power'),
              'energy_tolerance_eV': tolerance, 'states': retained, 'channels': channels,
              'top_leak_ledger': ledger, 'transition_audit': audit,
              'provenance_classes': {'sourced': 'citation and SHA256 required',
                                     'estimated': 'reserved; not supported by this generator'}}
    if not legacy:
        output['source_initial_levels'] = initial_levels
        output['source_domain'] = {'initial_levels': initial_levels, 'final_v': final_v}
    deck = ["# Species/manifold fragment. Install the supplied thermo.py as N2Vibrational.",
            "# Electron mappings and D6 execution require the separate coupling layer.",
            "database(thermoLibraries=['N2Vibrational'], kineticsFamilies='none')",
            "species(label='N2', structure=SMILES('N#N'))"]
    for v in range(1, K + 1):
        deck.append(f'species(label="N2v{v}", structure=adjacencyList({_adjacency(v)!r}))')
    deck.extend(["vibrationalManifold(species='N2')", 'simulator(atol=1e-16, rtol=1e-8)',
                 'model(toleranceMoveToCore=0.1, toleranceInterruptSimulation=1.0)'])
    return {'manifest.json': (json.dumps(output, sort_keys=True, indent=2, allow_nan=False) + '\n').encode(),
            'input.py': ('\n'.join(deck) + '\n').encode(), 'thermo.py': thermo_data}


def main():
    """Write a new artifact directory; never overwrite an earlier generation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source_manifest', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--K', type=int, default=10)
    parser.add_argument('--draft', action='store_true', help='permit unset D6 limits; runtime still refuses')
    args = parser.parse_args()
    artifacts = generate(json.loads(args.source_manifest.read_bytes()), args.K, not args.draft)
    args.output.mkdir(parents=True, exist_ok=False)
    for name, data in artifacts.items():
        (args.output / name).write_bytes(data)


if __name__ == '__main__':
    main()
