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

"""Synthetic ladder coupling; no physical molecular data or production claims."""
import copy
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from rmgpy.species import Species
from rmgpy.tools.eedf.integrity import _neutral_state_key, _channel_fractions, _resolved_properties
from rmgpy.tools.eedf.schema import SpecError, _content_hash
from rmgpy.tools.eedf.vibrational import (
    temperature, boltzmann_populations, ladder_coupling, resolve_levels,
    validate_ladders, energy_ledger, qualify_ladder, VibrationalPopulationError,
    D6RegenerationRequired, D6CeilingExceeded,
)

# Reuse the generator's own fixture without depending on pytest collection order.
_spec = importlib.util.spec_from_file_location(
    'vibrational_generator_fixture', Path(__file__).parents[1] / 'n2VibrationalTest.py')
_generator = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_generator)
source = _generator.source


@pytest.fixture
def small_ladder(source):
    levels_path = Path(source['level_source']['path'])
    data = json.loads(levels_path.read_text())
    data['levels'] = data['levels'][:3]
    levels_path.write_text(json.dumps(data))
    source['level_source'] = _generator.pin(levels_path)
    section_path = Path(source['cross_sections'][0]['path'])
    section_path.write_text(_generator.sections_text(2))
    section_pin = _generator.pin(section_path)
    section_pin['source_id'] = 'synthetic'
    source['cross_sections'] = [section_pin]
    source['source_domain'] = {'initial_levels': [0, 1, 2], 'final_v': 2}
    source['source_assignments'] = {str(v): 'synthetic' for v in range(3)}
    source['D6'].update(top_fraction_limit=1., top_leak_share_limit=1.)
    manifest = _generator.decoded(source, K=1)
    bindings = [{'source_sha256': pair['source']['sha256'], 'source_block': pair['source_block'],
                 'description': f"e + N2(X,v={pair['initial_v']}) <-> e + N2(X,v={pair['final_v']}), Vibrational"}
                for pair in reversed(manifest['transition_audit'])]
    coupling = ladder_coupling(manifest, bindings, axis='Tv_N2', parent='N2(X)')
    return manifest, coupling


@pytest.mark.parametrize('gas', ['N2', 'O2', 'NO'])
def test_hierarchical_selectors(gas):
    assert _neutral_state_key(gas + '(X,v=1)') == gas + '(X,v=1)'
    assert _neutral_state_key(gas + '(,X,v=1,)') == gas + '(X,v=1)'
    for bad in ('X,v=-1', 'X,v=1.0', 'X,v=*', 'X,v=1,J=0', '+,X,v=1'):
        with pytest.raises(SpecError):
            _neutral_state_key(gas + '(' + bad + ')')


def test_exact_species_identity_not_label(small_ladder):
    manifest, ladder = small_ladder
    species = [Species(label='deliberately misleading').from_adjacency_list(s['adjacency_list'])
               for s in reversed(manifest['states'])]
    assert resolve_levels(ladder, species) == (1, 0)
    species[0].molecule[0].vibrational_level = -1
    with pytest.raises(VibrationalPopulationError, match='missing'):
        resolve_levels(ladder, species)
    species = [Species().from_adjacency_list(s['adjacency_list']) for s in manifest['states']]
    with pytest.raises(VibrationalPopulationError, match='ambiguous'):
        resolve_levels(ladder, species + [species[0]])


def population_spec(ladder):
    return {'envelopes': {}, 'input_files': {}, 'solver_options': {},
            'gas_properties': {'fraction': ['N2 = 1']},
            'state_properties': {
                'population': ['N2(X) = 1'] + [
                    level['loki_state'] + ' = boltzmann:Tv_N2:' + str(level['v']) for level in ladder['levels']],
                'statisticalWeight': [s['loki_state'] + ' = ' + str(s['statistical_weight']) for s in ladder['levels']],
                'energy': [s['loki_state'] + ' = ' + str(s['energy_eV']) for s in ladder['levels']]},
            'vibrational_ladders': [ladder]}


def test_population_generation_and_terminal_numeric_override(small_ladder):
    _, ladder = small_ladder
    spec = population_spec(ladder)
    _, states, gas, populations = _resolved_properties(spec, {'Tv_N2': 800.})
    expected = boltzmann_populations(ladder, 800.)
    assert [populations[s['loki_state']] for s in ladder['levels']] == expected
    spec['state_properties'] = states
    spec['state_properties']['population'] = ['N2(X) = 1', 'N2(X,v=0) = 0.7',
                                               'N2(X,v=1) = 0.3', 'N2(X,v=2) = 0']
    _, _, _, exact = _resolved_properties(spec, {'Tv_N2': 800.})
    assert exact['N2(X,v=1)'] == .3
    channel = {'description': 'e + N2(X,v=0) <-> e + N2(X,v=1), Vibrational',
               'kind': 'vibrational', 'target_fraction': 'state:N2(X,v=0)',
               'product_fraction': 'state:N2(X,v=1)'}
    assert _channel_fractions([channel], gas, exact, {}) == ([.7], [.3])
    # Parent population participates in both collision directions.
    exact['N2(X)'] = .5
    assert _channel_fractions([channel], {'N2': .2}, exact, {}) == ([.2 * .5 * .7], [.2 * .5 * .3])
    del exact['N2(X)']
    with pytest.raises(SpecError, match='missing hierarchical'):
        _channel_fractions([channel], gas, exact, {})


@pytest.mark.parametrize('tv', [400., 800., 1200.])
def test_tv_roundtrip_and_gradient(small_ladder, tv):
    _, ladder = small_ladder
    populations = np.array(boltzmann_populations(ladder, tv)[:2]) * 2.3
    actual, derivative = temperature(ladder['levels'], populations)
    assert actual == pytest.approx(tv, rel=1e-13)
    for j in (0, 1):
        step = populations[j] * 1e-5
        low, high = populations.copy(), populations.copy()
        low[j] -= step
        high[j] += step
        fd = (temperature(ladder['levels'], high)[0] - temperature(ladder['levels'], low)[0]) / (2 * step)
        assert derivative[j] == pytest.approx(fd, rel=1e-8)


@pytest.mark.parametrize('populations', [[0., 1.], [1., 0.], [-1., .1], [1., 2.], [1., float('nan')], [1.]])
def test_tv_refuses_bad_population_limits(small_ladder, populations):
    with pytest.raises(VibrationalPopulationError):
        temperature(small_ladder[1]['levels'], populations)


def test_tv_infinite_limit_and_missing_level(small_ladder):
    levels = copy.deepcopy(small_ladder[1]['levels'])
    levels[1]['statistical_weight'] = levels[0]['statistical_weight']
    with pytest.raises(VibrationalPopulationError, match='infinite'):
        temperature(levels, [1., 1.])
    with pytest.raises(VibrationalPopulationError, match='requires levels'):
        temperature(levels[:1], [1., .1])


def test_binding_source_identity_and_policy(small_ladder):
    manifest, ladder = small_ladder
    assert len(ladder['channels']) == len(manifest['transition_audit'])
    for record in ladder['channels']:
        pair = next(p for p in manifest['transition_audit'] if p['source_block'] == record['source_block'])
        assert record['initial_v'] == pair['initial_v']
    bad = copy.deepcopy(ladder)
    bad['channels'][0]['description'] = bad['channels'][1]['description']
    with pytest.raises(SpecError):
        validate_ladders([bad])
    bad = copy.deepcopy(ladder)
    bad['D6']['top_fraction_limit'] = None
    with pytest.raises(_generator.LadderError):
        validate_ladders([bad])
    with pytest.raises(SpecError, match='missing source collision'):
        ladder_coupling(manifest, [], axis='Tv', parent='N2(X)')


def test_legacy_ladder_metadata_is_complete_only_over_its_levels(small_ladder):
    _, ladder = small_ladder
    legacy = copy.deepcopy(ladder)
    legacy.pop('source_domain')
    validate_ladders([legacy])
    legacy['channels'] = [c for c in legacy['channels']
                          if not (c['initial_v'] == 1 and c['source_final_v'] == 2)]
    with pytest.raises(SpecError, match='incomplete vibrational source-pair bindings'):
        validate_ladders([legacy])


@pytest.mark.parametrize('K', [10, 29])
def test_partial_source_domain_adapter_witness(source, tmp_path, K):
    spectrum_path = Path(source['level_source']['path'])
    spectrum = json.loads(spectrum_path.read_text())
    for v in range(len(spectrum['levels']), 59):
        spectrum['levels'].append({'v': v, 'energy_eV': .2 * v, 'statistical_weight': v + 1,
                                   'thermo_entry': None})
    spectrum_path.write_text(json.dumps(spectrum))
    source['level_source'] = _generator.pin(spectrum_path)
    sections_path = tmp_path / 'partial.txt'
    sections_path.write_text(_generator.partial_sections_text(range(30), 58))
    pin = _generator.pin(sections_path)
    pin['source_id'] = 'partial'
    source['cross_sections'] = [pin]
    source['source_domain'] = {'initial_levels': list(range(30)), 'final_v': 58}
    source['source_assignments'] = {str(v): 'partial' for v in range(30)}
    manifest = _generator.decoded(source, K=K)
    bindings = [{'source_sha256': pair['source']['sha256'], 'source_block': pair['source_block'],
                 'description': f"e + N2(X,v={pair['initial_v']}) <-> e + N2(X,v={pair['final_v']}), Vibrational"}
                for pair in manifest['transition_audit']]
    ladder = ladder_coupling(manifest, bindings, axis='Tv_N2', parent='N2(X)')
    assert ladder['source_domain'] == {'initial_levels': list(range(30)), 'final_v': 58}
    validate_ladders([ladder])
    missing = copy.deepcopy(ladder)
    missing['channels'] = [c for c in missing['channels']
                           if not (c['initial_v'] == 29 and c['source_final_v'] == 58)]
    with pytest.raises(SpecError, match='incomplete vibrational source-pair bindings'):
        validate_ladders([missing])


def test_d6_named_requests_and_closed_energy_partition(small_ladder):
    _, ladder = small_ladder
    channels = ladder['channels']
    row = SimpleNamespace(k_ine=np.full(len(channels), 2.), target_fractions=np.ones(len(channels)),
                          k_sup=np.zeros(len(channels)), product_fractions=np.zeros(len(channels)),
                          channel_power=np.array([2. * c['threshold_eV'] for c in channels]))
    ledger = energy_ledger(ladder, row, channels, [.9, .1], 1., 1.)
    assert ledger['forward_power'] == pytest.approx(
        ledger['retained_power'] + ledger['top_leak_power'] + ledger['source_residual_power'])
    assert ledger['top_leak_power'] > 0
    assert sum(ledger['top_leak_ledger'].values()) == pytest.approx(ledger['top_leak_power'])
    assert qualify_ladder(ladder, [.9, .1], ledger)['action'] == 'accept'
    ladder['D6']['top_fraction_limit'] = .05
    with pytest.raises(D6RegenerationRequired, match='D6_REGENERATE Tv_N2 K=6') as error:
        qualify_ladder(ladder, [.9, .1], ledger)
    assert error.value.decision['exceeded'] == ['top_fraction_limit']
    ladder['D6'].update(top_fraction_limit=1., top_leak_share_limit=.01)
    with pytest.raises(D6RegenerationRequired, match='top_leak_share_limit'):
        qualify_ladder(ladder, [.9, .1], ledger)
    ladder['K'] = 30
    with pytest.raises(D6CeilingExceeded, match='D6_REFUSE_K30'):
        qualify_ladder(ladder, [1.] + [0.] * 30, ledger)
    ladder['D6']['top_leak_share_limit'] = None
    with pytest.raises(_generator.LadderError):
        qualify_ladder(ladder, [1.] + [0.] * 30, ledger)


def test_grid_transfer_residual_is_visible_not_an_extra_charge(small_ladder):
    ladder = small_ladder[1]
    channels = ladder['channels']
    row = SimpleNamespace(k_ine=np.ones(len(channels)), k_sup=np.zeros(len(channels)),
                          target_fractions=np.ones(len(channels)), product_fractions=np.zeros(len(channels)),
                          channel_power=np.array([.9 * c['threshold_eV'] for c in channels]))
    ledger = energy_ledger(ladder, row, channels, None, 1., 1.)
    assert ledger['discretization_power'] < 0
    assert ledger['forward_power'] == pytest.approx(sum(ledger[name] for name in (
        'retained_power', 'top_leak_power', 'source_residual_power')))
    assert ledger['forward_power'] + ledger['discretization_power'] == pytest.approx(sum(row.channel_power))


def test_table_power_residual_survives_zero_net_flux(small_ladder):
    ladder = small_ladder[1]
    channels = ladder['channels']
    row = SimpleNamespace(k_ine=np.ones(len(channels)), k_sup=np.ones(len(channels)),
                          target_fractions=np.ones(len(channels)), product_fractions=np.ones(len(channels)),
                          channel_power=np.full(len(channels), .01))
    ledger = energy_ledger(ladder, row, channels, None, 1., 2.)
    retained_channels = [c for c in channels if c['initial_v'] <= ladder['K']]
    assert ledger['forward_power'] == pytest.approx(2. * sum(c['threshold_eV'] for c in retained_channels))
    assert ledger['discretization_power'] == pytest.approx(.02 * len(retained_channels))


def test_report_accepts_hierarchical_population_metadata():
    from rmgpy.tools.eedf.report import _validate_channel_metadata, ReportRefusal
    _validate_channel_metadata({'target_fraction': 'state:N2(X,v=0)',
                                'product_fraction': 'state:N2(X,v=1)'})
    with pytest.raises(ReportRefusal):
        _validate_channel_metadata({'target_fraction': 'state:N2(X,v=-1)'})


def test_ladder_metadata_changes_fingerprint(small_ladder):
    from rmgpy.tools.eedf.schema import _check_model, FingerprintMismatch
    model = {'vibrational_ladders': [small_ladder[1]]}
    manifest = {'row_inputs': model, 'fingerprint': _content_hash(model)}
    current = copy.deepcopy(model)
    current['vibrational_ladders'][0]['levels'][1]['energy_eV'] *= 1.01
    with pytest.raises(FingerprintMismatch, match='vibrational_ladders'):
        _check_model(manifest, current, 1e-9)


@pytest.mark.parametrize('property_name', ['energy', 'statisticalWeight'])
def test_solver_level_properties_must_match_source(small_ladder, property_name):
    spec = population_spec(small_ladder[1])
    spec['state_properties'][property_name][1] = 'N2(X,v=1) = 100'
    with pytest.raises(SpecError, match='differs from source'):
        _resolved_properties(spec, {'Tv_N2': 800.})


def test_loki_prepared_terminal_populations_normalize_per_parent(small_ladder):
    from rmgpy.tools.eedf.loki import LoKIDriver, LoKIError
    spec = population_spec(small_ladder[1])
    spec.update(Tg_K=300., P_Pa=100., working_conditions={},
                binary={'sha256': 'a' * 64}, cmake_cache={'sha256': 'b' * 64},
                channel_map={'sha256': 'c' * 64})
    _, properties, gas, populations = _resolved_properties(spec, {'Tv_N2': 800.})
    state = {'gas_fractions': gas, 'state_populations': populations, 'Tg_K': 300., 'P_Pa': 100.,
             'state_statistical_weights': {s['loki_state']: s['statistical_weight'] for s in small_ladder[1]['levels']}}
    driver = object.__new__(LoKIDriver)
    driver.spec = spec
    bundle = driver.prepare({'Tv_N2': 800.}, [2.], 'synthetic', terminal_state=state, verify_dependencies=False)
    assert b'N2(X,v=1)' in bundle.setup_bytes
    # Gas-level normalization would incorrectly sum the parent and children to two.
    assert sum(populations.values()) == pytest.approx(2.)
    spec['state_properties'] = properties
    spec['state_properties']['population'] = ['N2(X) = 1', 'N2(X,v=0) = 0.5',
                                               'N2(X,v=1) = 0.1', 'N2(X,v=2) = 0']
    state['state_populations'] = {'N2(X)': 1., 'N2(X,v=0)': .5, 'N2(X,v=1)': .1, 'N2(X,v=2)': 0.}
    with pytest.raises(LoKIError, match='normalize'):
        driver.prepare({'Tv_N2': 800.}, [2.], 'synthetic', terminal_state=state, verify_dependencies=False)


def test_tv_cannot_substitute_other_levels_for_zero_and_one(small_ladder):
    with pytest.raises(VibrationalPopulationError, match='requires levels'):
        temperature(small_ladder[1]['levels'][1:], [1., .1])


def test_two_ladders_cannot_alias_one_engine_state(small_ladder):
    original = small_ladder[1]
    alias = copy.deepcopy(original)
    alias['axis'], alias['parent'] = 'Tv_other', 'N2(A)'
    for level in alias['levels']:
        level['loki_state'] = level['loki_state'].replace('N2(X', 'N2(A')
    for channel in alias['channels']:
        channel['description'] = channel['description'].replace('N2(X', 'N2(A')
    with pytest.raises(SpecError, match='overlap'):
        validate_ladders([original, alias])


def test_above_k_engine_species_cannot_be_silently_zeroed(small_ladder):
    manifest, ladder = small_ladder
    species = [Species().from_adjacency_list(s['adjacency_list']) for s in manifest['states']]
    extra = species[-1].copy(deep=True)
    extra.molecule[0].vibrational_level = ladder['K'] + 1
    with pytest.raises(VibrationalPopulationError, match='above ladder K'):
        resolve_levels(ladder, species + [extra])
