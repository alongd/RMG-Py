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

"""Gas-independent vibrational coordinates and terminal D6 qualification.

The table owns a snapshot of the level energies, identities and collision
bindings. A source-pair index is never treated as a solver column index.
"""
import copy
import json
import math
import re

from scipy.constants import Boltzmann, elementary_charge

from rmgpy.tools.eedf.schema import EEDFError, SpecError, _content_hash
from rmgpy.tools.eedf.n2_vibrational import LadderError, _pin, _policy, assess_truncation


class VibrationalPopulationError(EEDFError):
    """A resolved vibrational population cannot define a positive finite Tv."""


class D6RegenerationRequired(EEDFError):
    """Qualification requires an explicit higher-K generation and restart."""

    def __init__(self, axis, decision):
        self.axis, self.decision = axis, decision
        super().__init__('D6_REGENERATE ' + axis + ' K=' + str(decision['K']) +
                         ': ' + ', '.join(decision['exceeded']))


class D6CeilingExceeded(EEDFError):
    """The owner limits cannot be met at the maximum supported K."""


def temperature(levels, populations):
    """Return Tv in K and (dTv/dn0, dTv/dn1) in K per population unit.

    Tv = (E1-E0)/kB / log(n0*g1/(n1*g0)). Both populations must be
    positive. Inversion is judged on weight-reduced populations n/g;
    n1/g1 >= n0/g0 refuses. The separate n1>n0 restriction also refuses
    even when unequal weights would permit a positive temperature.
    Missing levels and zero populations refuse; there is no finite clamp.
    """
    if (len(levels) < 2 or len(populations) < 2
            or levels[0].get('v') != 0 or levels[1].get('v') != 1):
        raise VibrationalPopulationError('Tv requires levels v=0 and v=1')
    n0, n1 = map(float, populations[:2])
    gap = (float(levels[1]['energy_eV']) - float(levels[0]['energy_eV'])) * elementary_charge / Boltzmann
    g0, g1 = (float(levels[i]['statistical_weight']) for i in (0, 1))
    if not all(math.isfinite(x) and x > 0 for x in (n0, n1, gap, g0, g1)):
        raise VibrationalPopulationError('Tv requires positive finite populations, gap and weights')
    ratio_log = math.log(n0) - math.log(n1) + math.log(g1) - math.log(g0)
    if n1 > n0 or ratio_log <= 0:
        raise VibrationalPopulationError('Tv population inversion or infinite-temperature limit')
    tv = gap / ratio_log
    gradient = (-tv / ratio_log / n0, tv / ratio_log / n1)
    if not all(math.isfinite(x) for x in (tv, *gradient)):
        raise VibrationalPopulationError('Tv or its derivative is nonfinite')
    return tv, gradient


def boltzmann_populations(ladder, tv):
    """Truncated Boltzmann table distribution; unresolved source levels are zero."""
    tv = float(tv)
    if not math.isfinite(tv) or tv <= 0:
        raise VibrationalPopulationError('Tv must be positive and finite')
    levels, K = ladder['levels'], ladder['K']
    logs = [math.log(s['statistical_weight']) - s['energy_eV'] * elementary_charge / (Boltzmann * tv)
            for s in levels[:K + 1]]
    peak = max(logs)
    weights = [math.exp(x - peak) for x in logs]
    total = sum(weights)
    result = [x / total for x in weights]
    # Table points must obey the same domain as reactor coordinates.
    temperature(levels, result)
    return result + [0.] * (len(levels) - K - 1)


def ladder_coupling(manifest, bindings, *, axis, parent, d6_policy=None):
    """Adapt an emitted ladder to explicitly bound LoKI collision descriptions.

    ``bindings`` contains {source_sha256, source_block, description} for every
    audited source pair, including unresolved initial states. The source hash
    and block, not the generator pair index, identify each solver collision.
    The returned snapshot goes in spec['vibrational_ladders']; no input is
    mutated. The parent names the solver's neutral electronic state.
    """
    from rmgpy.molecule.molecule import Molecule
    from rmgpy.tools.eedf.integrity import _neutral_state_key
    parent = _neutral_state_key(parent)
    if ',v=' in parent:
        raise SpecError('vibrational parent must be an electronic state')
    states = manifest['states']
    spectrum = json.loads(_pin(states[0]['provenance']))
    source_levels = spectrum['levels']
    indexed = {}
    for binding in bindings:
        key = (binding['source_sha256'], binding['source_block'])
        if key in indexed:
            raise SpecError('duplicate source collision binding')
        indexed[key] = binding['description']
    levels = []
    ground = Molecule().from_adjacency_list(states[0]['adjacency_list'])
    formula = parent.split('(', 1)[0]
    if ground.get_formula() != formula:
        raise SpecError('ladder parent differs from molecular formula')
    for v, level in enumerate(source_levels):
        if v <= manifest['K']:
            molecule = Molecule().from_adjacency_list(states[v]['adjacency_list'])
            if (molecule.get_formula() != formula or molecule.vibrational_level != v
                    or molecule.electronic_state != ground.electronic_state
                    or molecule.multiplicity != ground.multiplicity
                    or any(states[v][key] != level[key]
                           for key in ('v', 'energy_eV', 'statistical_weight'))):
                raise SpecError('ladder state differs from sourced level identity')
        levels.append({'v': v, 'energy_eV': level['energy_eV'],
                       'statistical_weight': level['statistical_weight'],
                       'formula': formula, 'electronic_state': ground.electronic_state,
                       'multiplicity': ground.multiplicity, 'vibrational_level': v,
                       'loki_state': parent[:-1] + ',v=' + str(v) + ')'})
    channels = []
    for pair in manifest['transition_audit']:
        key = (pair['source']['sha256'], pair['source_block'])
        if key not in indexed:
            raise SpecError('missing source collision binding: ' + str(key))
        i, j = pair['initial_v'], pair['final_v']
        channels.append({'description': indexed.pop(key), 'initial_v': i, 'source_final_v': j,
                         'final_v': min(j, manifest['K']), 'threshold_eV': pair['threshold_eV'],
                         'top_leak_eV': (source_levels[j]['energy_eV'] -
                                         source_levels[min(j, manifest['K'])]['energy_eV'])
                                        if i <= manifest['K'] else 0.,
                         'source_sha256': key[0], 'source_block': key[1]})
    if indexed:
        raise SpecError('unexpected source collision binding')
    result = {'axis': axis, 'parent': parent, 'K': manifest['K'],
              'levels': levels, 'channels': channels,
              'D6': copy.deepcopy(manifest['D6'] if d6_policy is None else d6_policy),
              'source_manifest_sha256': _content_hash(manifest),
              'level_source': copy.deepcopy(states[0]['provenance'])}
    if 'source_domain' in manifest:
        result['source_domain'] = copy.deepcopy(manifest['source_domain'])
    validate_ladders([result])
    return result


def validate_ladders(ladders):
    """Validate the frozen coordinate/identity/power schema without filesystem IO."""
    from rmgpy.tools.eedf.integrity import _neutral_state_key
    axes, parents, descriptions, identities = set(), set(), set(), set()
    for ladder in ladders:
        axis, parent, K = ladder['axis'], ladder['parent'], ladder['K']
        if (not isinstance(axis, str) or not re.fullmatch('[A-Za-z_][A-Za-z_0-9]*', axis)
                or axis == 'u' or axis in axes or parent in parents
                or _neutral_state_key(parent) != parent or ',v=' in parent
                or type(K) is not int or not 1 <= K <= 30):
            raise SpecError('invalid/duplicate vibrational axis or parent')
        axes.add(axis)
        parents.add(parent)
        _policy(ladder['D6'])
        levels = ladder['levels']
        if len(levels) <= K:
            raise SpecError('vibrational levels do not cover K')
        tail = ladder.get('zero_population_tail')
        if tail is not None:
            expected_tail = {level['loki_state'] for level in levels if level['v'] > K}
            if (not isinstance(tail, dict) or set(tail) != expected_tail
                    or any(type(value) not in (int, float) or value != 0. for value in tail.values())):
                raise SpecError('invalid explicit zero-population vibrational tail')
        previous = -1.
        for v, level in enumerate(levels):
            energy, weight = level['energy_eV'], level['statistical_weight']
            if (level['v'] != v or level['vibrational_level'] != v
                    or level['loki_state'] != parent[:-1] + ',v=' + str(v) + ')'
                    or level['formula'] != parent.split('(')[0]
                    or not isinstance(level['electronic_state'], str)
                    or type(level['multiplicity']) is not int or level['multiplicity'] < 1
                    or not math.isfinite(energy) or energy <= previous or (v == 0 and energy != 0)
                    or not math.isfinite(weight) or weight <= 0):
                raise SpecError('invalid vibrational level identity, energy or weight')
            previous = energy
            identity = (level['formula'], level['electronic_state'], v, level['multiplicity'])
            if identity in identities:
                raise SpecError('vibrational ladders overlap the same engine state')
            identities.add(identity)
        source_domain = ladder.get('source_domain')
        if 'source_domain' not in ladder:
            # Old and synthetic coupling snapshots imply complete coverage
            # only over their own level list.
            initial_levels = list(range(len(levels)))
            final_v = len(levels) - 1
        else:
            if (not isinstance(source_domain, dict)
                    or set(source_domain) != {'initial_levels', 'final_v'}):
                raise SpecError('invalid vibrational source-domain metadata')
            initial_levels = source_domain['initial_levels']
            final_v = source_domain['final_v']
            if (not isinstance(initial_levels, list) or not initial_levels
                    or any(type(v) is not int for v in initial_levels)
                    or initial_levels != sorted(set(initial_levels))
                    or type(final_v) is not int or not 1 <= final_v < len(levels)
                    or initial_levels[0] < 0 or initial_levels[-1] > final_v):
                raise SpecError('invalid vibrational source-domain coverage')
            if final_v < K or any(v not in initial_levels for v in range(K + 1)):
                raise SpecError('vibrational source-domain does not cover every retained initial level')
        pairs = set()
        for channel in ladder['channels']:
            i, j = channel['initial_v'], channel['source_final_v']
            description = channel['description']
            if (type(i) is not int or type(j) is not int or not 0 <= i <= j < len(levels)
                    or (i, j) in pairs or description in descriptions
                    or channel['final_v'] != min(j, K)):
                raise SpecError('invalid/duplicate vibrational collision binding')
            if i not in initial_levels or j > final_v:
                raise SpecError('vibrational collision lies outside declared source-domain coverage')
            pairs.add((i, j))
            descriptions.add(description)
            target, product = levels[i]['loki_state'], levels[j]['loki_state']
            if description not in (f'e + {target} <-> e + {product}, Vibrational',
                                    f'e + {target} -> e + {product}, Vibrational'):
                raise SpecError('vibrational collision binding has different state identities')
            if i < j <= K and ' <-> ' not in description:
                raise SpecError('retained vibrational pair requires a superelastic solver column')
            if (not math.isfinite(channel['threshold_eV']) or channel['threshold_eV'] < 0
                    or channel['top_leak_eV'] != (levels[j]['energy_eV'] - levels[min(j, K)]['energy_eV']
                                                  if i <= K else 0.)):
                raise SpecError('invalid vibrational collision energy ledger')
        if any((i, j) not in pairs for i in initial_levels for j in range(i + 1, final_v + 1)):
            raise SpecError('incomplete vibrational source-pair bindings')


def resolve_levels(ladder, species):
    """Resolve each retained level exactly once by identity, never display label."""
    indices = []
    ground = ladder['levels'][0]
    for spc in species:
        if spc.get_net_charge() == 0 and any(
                m.get_formula() == ground['formula']
                and m.electronic_state == ground['electronic_state']
                and m.multiplicity == ground['multiplicity']
                and m.vibrational_level > ladder['K'] for m in spc.molecule):
            raise VibrationalPopulationError('resolved species above ladder K cannot be assigned zero population')
    for level in ladder['levels'][:ladder['K'] + 1]:
        matches = []
        for index, spc in enumerate(species):
            if not spc.molecule or spc.get_net_charge() != 0:
                continue
            identities = {(m.get_formula(), m.electronic_state, m.vibrational_level, m.multiplicity)
                          for m in spc.molecule}
            expected = (level['formula'], level['electronic_state'], level['v'], level['multiplicity'])
            if identities == {expected}:
                matches.append(index)
        if len(matches) != 1:
            raise VibrationalPopulationError('missing or ambiguous exact level ' + level['loki_state'])
        indices.append(matches[0])
    return tuple(indices)


def validate_bound_channels(ladder, channels, reactions, reaction_map, indices, species):
    """Require exact inelastic/superelastic stoichiometry and all retained owners."""
    columns = {c['description']: j for j, c in enumerate(channels)}
    for binding in ladder['channels']:
        if binding['description'] not in columns:
            raise SpecError('missing vibrational solver column')
        column = columns[binding['description']]
        channel = channels[column]
        if channel['kind'] != 'vibrational' or channel['threshold_eV'] != binding['threshold_eV']:
            raise SpecError('vibrational column differs from source threshold/kind')
        i, j, K = binding['initial_v'], binding['source_final_v'], ladder['K']
        owners = {side: index for index, (col, side, *_) in reaction_map.items() if col == column}
        expected_sides = ({'ine', 'sup'} if i < j <= K else {'ine'} if i < K and j > K else set())
        if set(owners) != expected_sides:
            raise SpecError('vibrational collision requires exact retained ine/sup owners')
        if not owners and channel.get('classification') != 'B':
            raise SpecError('unresolved/diagonal/top-to-top vibrational collision must be energy-only B')
        for side, index in owners.items():
            reaction = reactions[index]
            initial, final = (i, min(j, K)) if side == 'ine' else (j, i)
            for actual, target in ((reaction.reactants, species[indices[initial]]),
                                   (reaction.products, species[indices[final]])):
                if len(actual) != 2 or sum(s.is_electron() for s in actual) != 1 or target not in actual:
                    raise SpecError('vibrational reaction does not bind the exact level population')
            if reaction.reversible:
                raise SpecError('vibrational EEDF reactions must be irreversible')


def energy_ledger(ladder, row, channel_map, populations, neutral, power_factor):
    """Gross forward power and leak diagnostics, never extra electron charges.

    Rates use the actual resolved target amounts and sourced thresholds,
    including the full threshold for a transition above K. The electron
    energy row already charges channel power; this ledger adds no charge.
    ``discretization_power`` is the table's net power minus its net flux
    times the sourced threshold, and never enters gross forward power/D6.
    """
    columns = {c['description']: j for j, c in enumerate(channel_map)}
    forward = leak = retained = residual = discretization = 0.
    per_leak = {}
    for binding in ladder['channels']:
        i = binding['initial_v']
        if i > ladder['K']:
            continue
        column = columns[binding['description']]
        target = float(row.target_fractions[column]) if populations is None else populations[i] / neutral
        flux = power_factor * target * float(row.k_ine[column])
        net = (float(row.target_fractions[column]) * float(row.k_ine[column]) -
               float(row.product_fractions[column]) * float(row.k_sup[column]))
        loss = flux * binding['threshold_eV']
        delta = ladder['levels'][binding['final_v']]['energy_eV'] - ladder['levels'][i]['energy_eV']
        entry_leak = flux * binding['top_leak_eV']
        forward += loss
        leak += entry_leak
        retained += flux * delta
        residual += flux * (binding['threshold_eV'] - delta - binding['top_leak_eV'])
        discretization += power_factor * (float(row.channel_power[column]) - net * binding['threshold_eV'])
        if binding['top_leak_eV']:
            per_leak['top-leak:' + binding['description']] = entry_leak
    return {'forward_power': forward, 'top_leak_power': leak, 'retained_power': retained,
            'source_residual_power': residual, 'discretization_power': discretization,
            'top_leak_ledger': per_leak}


def qualify_ladder(ladder, populations, ledger):
    """Apply owner-set D6 limits before a terminal state can qualify."""
    try:
        decision = assess_truncation(ladder['K'], populations, ledger['forward_power'],
                                     ledger['top_leak_power'], ladder['D6'])
    except LadderError as error:
        if ladder['K'] == 30 and 'D6 ceiling' in str(error):
            raise D6CeilingExceeded('D6_REFUSE_K30 ' + ladder['axis'] + ': ' + str(error)) from error
        raise
    if decision['action'] == 'regenerate':
        raise D6RegenerationRequired(ladder['axis'], decision)
    return decision
