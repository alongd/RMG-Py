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

import copy
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from rmgpy.tools.eedf import channel_map as adapter
from rmgpy.tools.eedf.schema import FingerprintMismatch, SpecError, _file_hash, load_spec
from branchYieldTest import PROCESS, branch_fixture, branch_reactions, database_cascade_channels, source_files


@pytest.fixture
def conversion_inputs(tmp_path):
    # Small classic-format source, independent of the converter and database.
    process = 'e + Ar(1S0) -> e + Ar(3P2), Excitation'
    section = tmp_path / 'Ar_LXCat.txt'
    section.write_text('PARAM.: E = 11.548 eV\n'
                       '[Ar(1S0) + e -> Ar(3P2) + e, Excitation]\n'
                       '-----\n11.548 0\n13 3.237651e-22\n-----\n')
    physical = {'description': process, 'kind': 'excitation',
                'classification': 'B', 'threshold_eV': 11.548,
                'target_fraction': 1., 'product_fraction': 0.,
                'sigma_max_m2': 3.237651e-22, 'flux_group': process,
                'reaction': None,
                'cross_section': {'energy_eV': [11.548, 13.],
                                  'sigma_m2': [0., 3.237651e-22]}}
    mapping = tmp_path / 'channels.json'
    mapping.write_text(json.dumps([physical]))
    spec = {'channel_map': {'path': str(mapping), 'sha256': _file_hash(mapping)},
            'input_files': {'Ar_LXCat.txt': {'path': str(section),
                            'sha256': _file_hash(section), 'kind': 'cross_section'}},
            'axes': {'u': [0., 1.]}, 'envelopes': {},
            'gas_properties': {'fraction': ['Ar = 1']},
            'state_properties': {'population': ['Ar(1S0) = 1']},
            'solver_options': {}}
    frozen = {'schema': 'eedf-channel-map/2',
              'loki_files': {'Ar_LXCat.txt': _file_hash(section)},
              'class_a_gaps': {}, 'class_counts': {'A': 1, 'B': 0, 'C': 0, 'D': 0},
              'channels': [{'id': 2, 'loki_process': process, 'lxcat_file': 'Ar_LXCat.txt',
                            'threshold_eV': 11.548, 'class': 'A',
                            'rmg': {'library': 'PlasmaArgon', 'index': 87,
                                    'status': 'present', 'reaction': 'Ar + e- => Ars + e-'}}]}
    reaction = SimpleNamespace(reactants=[SimpleNamespace(label=s) for s in ['Ar', 'e-']],
                               products=[SimpleNamespace(label=s) for s in ['Ars', 'e-']],
                               reversible=False)
    library = SimpleNamespace(entries={87: SimpleNamespace(item=reaction)})
    return spec, frozen, library


def test_converts_complete_map_without_mutating_history(conversion_inputs, tmp_path):
    spec, frozen, library = conversion_inputs
    original = Path(spec['channel_map']['path']).read_bytes()
    converted, report = adapter.convert_channel_map(spec, frozen, library, 'test-Ar')
    assert report['class_counts'] == {'A': 1}
    assert report['missing_owners'] == []
    assert converted[0]['classification'] == 'A'
    assert converted[0]['reaction'] == {
        'library': 'PlasmaArgon', 'index': 87,
        'repr': "EEDFChannel(process='e + Ar(1S0) -> e + Ar(3P2), Excitation', "
                "collision_set='test-Ar', side='ine')"}
    assert converted[0]['cross_section'] == {
        'energy_eV': [11.548, 13.], 'sigma_m2': [0., 3.237651e-22]}
    assert Path(spec['channel_map']['path']).read_bytes() == original
    destination = tmp_path / 'converted'
    adapter._write_outputs(destination, report, spec, converted)
    output = json.loads((destination / 'pure-ar-spec.json').read_text())
    assert output['channel_map']['sha256'] == _file_hash(destination / 'channel-map.json')
    assert output['input_files'] == spec['input_files']
    assert Path(spec['channel_map']['path']).read_bytes() == original
    with pytest.raises(FileExistsError):
        adapter._write_outputs(destination, report, spec, converted)
    assert output == json.loads((destination / 'pure-ar-spec.json').read_text())


def test_owner_problem_keeps_single_library_call_compatibility(conversion_inputs):
    _, frozen, library = conversion_inputs
    assert adapter._owner_problem(frozen['channels'][0], library) is None


def _two_gas_fixture(conversion_inputs, tmp_path):
    spec, ar_map, ar_library = conversion_inputs
    spec = copy.deepcopy(spec)
    ar_map = copy.deepcopy(ar_map)
    n2_path = tmp_path / 'N2_LXCat.txt'
    n2_path.write_text('PARAM.: E = 0.01 eV\n'
                       '[N2(X,v=0) + e -> N2(X,v=1) + e, Vibrational]\n'
                       '-----\n0.01 0\n1 1e-22\n-----\n')
    process = 'e + N2(X,v=0) -> e + N2(X,v=1), Vibrational'
    n2_hash = _file_hash(n2_path)
    spec['input_files']['N2_LXCat.txt'] = {
        'path': str(n2_path), 'sha256': n2_hash, 'kind': 'cross_section'}
    spec['axes']['x_N2'] = [.5]
    spec['gas_properties']['fraction'] = ['Ar = 1-{x_N2}', 'N2 = {x_N2}']
    spec['state_properties']['population'] = [
        'Ar(1S0) = 1', 'N2(X) = 1', 'N2(X,v=0) = 1', 'N2(X,v=1) = 1']
    channel_map = json.loads(Path(spec['channel_map']['path']).read_text())
    channel_map[0]['target_fraction'] = .5
    channel_map[0]['product_fraction'] = 0.
    channel_map.append({
        'description': process, 'kind': 'vibrational', 'classification': 'B',
        'threshold_eV': .01, 'target_fraction': .5, 'product_fraction': .5,
        'sigma_max_m2': 1e-22, 'flux_group': process,
        'reaction': None,
        'cross_section': {'energy_eV': [.01, 1.], 'sigma_m2': [0., 1e-22]}})
    map_path = Path(spec['channel_map']['path'])
    map_path.write_text(json.dumps(channel_map))
    spec['channel_map']['sha256'] = _file_hash(map_path)
    n2_map = copy.deepcopy(ar_map)
    n2_map['loki_files'] = {'N2_LXCat.txt': n2_hash}
    n2_map['class_a_gaps'] = {}
    n2_map['class_counts'] = {'A': 1, 'B': 0, 'C': 0, 'D': 0}
    n2_map['channels'] = [{
        'id': 9, 'loki_process': process, 'lxcat_file': 'N2_LXCat.txt',
        'threshold_eV': .01, 'class': 'A',
        'rmg': {'library': 'PlasmaNitrogen', 'index': 12, 'status': 'present',
                'reaction': 'N2 + e- => N2v1 + e-'}}]
    reaction = SimpleNamespace(reactants=[SimpleNamespace(label=s) for s in ['N2', 'e-']],
                               products=[SimpleNamespace(label=s) for s in ['N2v1', 'e-']],
                               reversible=False)
    n2_library = SimpleNamespace(entries={12: SimpleNamespace(item=reaction)})
    return spec, ar_map, n2_map, {'PlasmaArgon': ar_library, 'PlasmaNitrogen': n2_library}


def test_multigas_conversion_uses_each_configured_owner_library(conversion_inputs, tmp_path, monkeypatch):
    spec, ar_map, n2_map, libraries = _two_gas_fixture(conversion_inputs, tmp_path)
    channels, report = adapter.convert_channel_maps(
        spec, {'Ar': ar_map, 'N2': n2_map}, libraries, 'mix-Ar-N2',
        library_by_gas={'Ar': 'PlasmaArgon', 'N2': 'PlasmaNitrogen'})
    assert len(channels) == 2
    assert [channel['reaction']['library'] for channel in channels] == [
        'PlasmaArgon', 'PlasmaNitrogen']
    assert report['class_counts'] == {'A': 2}
    import rmgpy.tools.eedf.vibrational as vibrational
    ladder = {'axis': 'Tv_K', 'parent': 'N2(X)', 'K': 1, 'levels': [
        {'v': 0, 'loki_state': 'N2(X,v=0)', 'energy_eV': 0., 'statistical_weight': 1.},
        {'v': 1, 'loki_state': 'N2(X,v=1)', 'energy_eV': .29, 'statistical_weight': 1.}],
        'source_domain': {'initial_levels': list(range(30)), 'final_v': 58}}
    monkeypatch.setattr(vibrational, 'ladder_coupling', lambda *args, **kwargs: ladder)
    mixture = adapter.build_mixture_spec(spec, {}, [])
    destination = tmp_path / 'mixture'
    adapter._write_outputs(destination, report, spec, channels, mixture)
    published = json.loads((destination / 'ar-n2-spec.json').read_text())
    assert json.loads((destination / 'channel-map.json').read_text()) == channels
    assert published['channel_map']['sha256'] == _file_hash(destination / 'channel-map.json')
    assert published['axes']['Tv_K'] == [1000., 2000.]
    assert not (destination / 'pure-ar-spec.json').exists()
    from rmgpy.tools.eedf.integrity import _resolved_properties
    _, _, fractions, populations = _resolved_properties(
        published, {'x_N2': .25, 'Tv_K': 1000.})
    assert fractions == {'Ar': .75, 'N2': .25}
    assert populations['N2(X,v=0)'] + populations['N2(X,v=1)'] == pytest.approx(1.)


def _three_gas_state_maps():
    def entry(g, rmg_species):
        return {'paschen': 'ground', 'g': g, 'excitation_eV': 0.0,
                'rmg_species': rmg_species}
    return {'Ar': {'states': {'Ar(1S0)': entry(1, 'Ar')}},
            'N2': {'states': {'N2(X)': entry(1, 'N2')}},
            'NO': {'states': {'NO(X)': entry(4, 'NO')}}}


def _fake_ladder(monkeypatch):
    import rmgpy.tools.eedf.vibrational as vibrational
    ladder = {'axis': 'Tv_K', 'parent': 'N2(X)', 'K': 1, 'levels': [
        {'v': 0, 'loki_state': 'N2(X,v=0)', 'energy_eV': 0., 'statistical_weight': 1.},
        {'v': 1, 'loki_state': 'N2(X,v=1)', 'energy_eV': .29, 'statistical_weight': 1.}],
        'source_domain': {'initial_levels': list(range(30)), 'final_v': 58}}
    monkeypatch.setattr(vibrational, 'ladder_coupling', lambda *args, **kwargs: ladder)


def test_three_gas_mixture_declares_named_bath_fractions_and_map_states(
        conversion_inputs, monkeypatch):
    spec, _, _ = conversion_inputs
    _fake_ladder(monkeypatch)
    maps = _three_gas_state_maps()
    maps = {gas: maps[gas] for gas in ('NO', 'N2', 'Ar')}
    composition = {'bath': 'Ar', 'fractions': {
        'N2': {'axis': 'x_N2', 'values': [.1]},
        'NO': {'axis': 'x_NO', 'values': [.001, .005]}}}
    mixture = adapter.build_mixture_spec(spec, {}, [], gas_maps=maps,
                                         composition=composition)
    assert 'N2 = 0.1' in mixture['gas_properties']['fraction']
    assert 'x_N2' not in mixture['axes']
    assert mixture['gas_properties']['fraction'] == [
        'NO = {x_NO}', 'N2 = 0.1', 'Ar = 1-{x_NO}-0.1']
    assert 'NO(X) = 1' in mixture['state_properties']['population']
    assert 'NO(X) = 0.0' in mixture['state_properties']['energy']
    assert 'NO(X) = 4' in mixture['state_properties']['statisticalWeight']
    assert mixture['axes']['x_NO'] == [.001, .005]
    from rmgpy.tools.eedf.integrity import _resolved_properties
    _, _, fractions, _ = _resolved_properties(
        mixture, {'x_N2': .1, 'x_NO': .001, 'Tv_K': 300.})
    assert fractions == {'NO': .001, 'N2': .1, 'Ar': .899}


def test_three_gas_mixture_passes_physical_map_state_validation(
        conversion_inputs, monkeypatch):
    """The declared path must survive the same state checks as conversion."""
    spec, _, _ = conversion_inputs
    spec['state_properties']['population'].extend(
        f'N2(X,v={v}) = {1 if v == 0 else 0}' for v in range(59))
    _fake_ladder(monkeypatch)
    maps = _three_gas_state_maps()
    mixture = adapter.build_mixture_spec(
        spec, {}, [], gas_maps=maps,
        composition={'bath': 'Ar', 'fractions': {
            'N2': {'axis': 'x_N2', 'values': [.1]},
            'NO': {'axis': 'x_NO', 'values': [.001, .005]}}})

    # Keep _physical_properties and its selector conflict detection real; the
    # fixture's physical sections are irrelevant to this state regression.
    from rmgpy.tools.eedf import channels
    monkeypatch.setattr(channels, '_validate_map', lambda *_args: None)
    monkeypatch.setattr(channels, '_classic_sections', lambda *args: [])
    monkeypatch.setattr(channels, '_resolved_properties',
                        lambda *_args: (None, {}, None, None))
    channels._validate_physical_map([], mixture, {})


@pytest.mark.parametrize('composition,maps,error', [
    ({'fractions': {}}, _three_gas_state_maps(), 'bath'),
    ({'bath': 'Ne', 'fractions': {}}, _three_gas_state_maps(), 'bath'),
    ({'bath': 'Ar', 'fractions': {
        'N2': {'axis': 'x_N2', 'values': [.1]},
        'NO': {'axis': 'x_NO', 'values': [.1]}}},
     {**_three_gas_state_maps(), 'NO': {'states': {}}}, 'ground state'),
    ({'bath': 'Ar', 'fractions': {'N2': {'axis': 'x_N2', 'values': [.6]},
                                  'NO': {'axis': 'x_NO', 'values': [.5]}}},
     _three_gas_state_maps(), 'fractions'),
    ({'bath': 'Ar', 'fractions': {'N2': {'axis': 'x_N2', 'values': [.1]},
                                  'NO': {'axis': 'x_NO', 'values': [.1]}}},
     {'Ar': _three_gas_state_maps()['Ar'], 'N2': _three_gas_state_maps()['N2']}, 'map'),
])
def test_general_mixture_refuses_invalid_composition_or_missing_state_map(
        conversion_inputs, monkeypatch, composition, maps, error):
    spec, _, _ = conversion_inputs
    _fake_ladder(monkeypatch)
    with pytest.raises(SpecError, match=error):
        adapter.build_mixture_spec(spec, {}, [], gas_maps=maps, composition=composition)


def test_prepare_refuses_duplicate_gas_identity(conversion_inputs, tmp_path, monkeypatch):
    spec, _, _ = conversion_inputs
    monkeypatch.setattr(adapter, 'load_spec', lambda _: copy.deepcopy(spec))
    monkeypatch.setattr(adapter, '_diagnostic_d6', lambda _: None)
    path = tmp_path / 'argon.yaml'
    path.write_text(yaml.safe_dump(_three_gas_state_maps()['Ar']))
    digest = _file_hash(path)
    with pytest.raises(SpecError, match='gas map'):
        adapter.prepare_mixture_spec(
            'unused', [
                {'gas': 'Ar', 'path': str(path), 'sha256': digest, 'library': 'PlasmaArgon'},
                {'gas': 'Ar', 'path': 'unused', 'sha256': 'unused', 'library': 'PlasmaArgon'}],
            'unused', 'mixture', {'D6': None})


def test_prepare_refuses_ungoverned_multigas_legacy_composition(
        conversion_inputs, monkeypatch):
    spec, _, _ = conversion_inputs
    monkeypatch.setattr(adapter, 'load_spec', lambda _: copy.deepcopy(spec))
    monkeypatch.setattr(adapter, '_diagnostic_d6', lambda _: None)
    with pytest.raises(SpecError, match='explicitly named bath'):
        adapter.prepare_mixture_spec(
            'unused', [{'gas': gas, 'path': 'unused', 'sha256': 'unused', 'library': gas}
                       for gas in ('Ar', 'N2', 'NO')],
            'unused', 'mixture', {'D6': None})


def test_prepare_three_gas_conversion_passes_declared_bath_and_states_to_map_converter(
        conversion_inputs, tmp_path, monkeypatch):
    spec, _, _ = conversion_inputs
    masses = tmp_path / 'masses.txt'
    masses.write_text('Ar 6.6e-26\nN2 4.6e-26\nNO 5.0e-26\n')
    spec['input_files']['masses.txt'] = {
        'path': str(masses), 'sha256': _file_hash(masses), 'kind': 'property'}
    spec['gas_properties']['mass'] = 'masses.txt'
    _fake_ladder(monkeypatch)
    monkeypatch.setattr(adapter, 'load_spec', lambda _: copy.deepcopy(spec))
    monkeypatch.setattr(adapter, '_diagnostic_d6', lambda _: None)
    monkeypatch.setattr(adapter, '_load_library', lambda database, scratch, label:
                        (SimpleNamespace(entries={}), {'fixture': label}))
    monkeypatch.setattr(adapter, '_cascade_snapshot', lambda *args, **kwargs: ({}, None))
    import rmgpy.tools.eedf.n2_vibrational as n2_vibrational
    monkeypatch.setattr(n2_vibrational, 'generate', lambda *args, **kwargs:
                        {'manifest.json': json.dumps({'transition_audit': []})})
    converted = [{'description': 'fixture collision'}]
    report = {'outcome': 'converted', 'missing_owners': []}
    observed = {}

    def convert(spec_snapshot, maps, libraries, collision_set, yields, library_by_gas):
        observed['fractions'] = list(spec_snapshot['gas_properties']['fraction'])
        observed['population'] = list(spec_snapshot['state_properties']['population'])
        observed['energy'] = list(spec_snapshot['state_properties']['energy'])
        observed['weight'] = list(spec_snapshot['state_properties']['statisticalWeight'])
        return copy.deepcopy(converted), copy.deepcopy(report)

    monkeypatch.setattr(adapter, 'convert_channel_maps', convert)
    gas_map_inputs = []
    for gas, frozen, library in (
            ('Ar', _three_gas_state_maps()['Ar'], 'PlasmaArgon'),
            ('N2', _three_gas_state_maps()['N2'], 'PlasmaAir'),
            ('NO', _three_gas_state_maps()['NO'], 'PlasmaNitricOxide')):
        frozen = dict(frozen, schema='eedf-channel-map/2', set=gas,
                      loki_files={'fixture.txt': 'fixture'},
                      class_a_gaps={}, class_counts={'A': 0, 'B': 0, 'C': 0, 'D': 0},
                      channels=[])
        path = tmp_path / f'{gas}.yaml'
        path.write_text(yaml.safe_dump(frozen))
        gas_map_inputs.append({'gas': gas, 'path': str(path), 'sha256': _file_hash(path),
                               'library': library})
    composition = {'bath': 'Ar', 'fractions': {
        'N2': {'axis': 'x_N2', 'values': [.1]},
        'NO': {'axis': 'x_NO', 'values': [.001, .005]}}}
    _, returned_channels, returned_report, mixture = adapter.prepare_mixture_spec(
        'unused', gas_map_inputs, 'unused', 'test-Ar-N2-NO', {'D6': None},
        composition=composition, tv_values=(300., 3000.))
    assert returned_channels == converted and returned_report['outcome'] == 'converted'
    assert observed['fractions'] == [
        'Ar = 1-0.1-{x_NO}', 'N2 = 0.1', 'NO = {x_NO}']
    assert 'NO(X) = 1' in observed['population']
    assert 'NO(X) = 0.0' in observed['energy']
    assert 'NO(X) = 4' in observed['weight']
    assert mixture['axes']['x_NO'] == [.001, .005]


def test_supplemented_gas_collision_passes_physical_map_validation(
        tmp_path, monkeypatch):
    from rmgpy.tools.eedf import channels

    process = 'e + NO(X) -> e + NO(V1), Excitation'
    section = {'description': process, 'kind': 'excitation', 'threshold_eV': .232,
               'cross_section': {'energy_eV': [.1, 1.], 'sigma_m2': [0., 2.e-20]}}
    frozen = {'id': 3, 'loki_process': process, 'class': 'B', 'rmg': None,
              'flux_group': 'NO_lumped_vibrational_heating'}
    channel = adapter._physical_channel_from_frozen(section, frozen, 'fixture')
    assert channel['target_fraction'] == 'state:NO(X)'
    assert channel['product_fraction'] == 0.
    assert channel['flux_group'] == 'NO_lumped_vibrational_heating'

    source = tmp_path / 'NO.txt'
    source.write_text('fixture bytes')
    spec = {
        'input_files': {'NO.txt': {'path': str(source), 'sha256': _file_hash(source),
                                   'kind': 'cross_section'}},
        'gas_properties': {'fraction': ['NO = 1']},
        'state_properties': {'population': ['NO(X) = 1'], 'energy': ['NO(X) = 0'],
                             'statisticalWeight': ['NO(X) = 4']},
        'solver_options': {},
    }
    monkeypatch.setattr(channels, '_classic_sections', lambda *_args: [section])
    monkeypatch.setattr(channels, '_resolved_properties',
                        lambda *_args: (None, {}, None, None))
    channels._validate_physical_map([channel], spec, {})


def test_gas_map_species_supplies_missing_loki_particle_mass():
    species = SimpleNamespace(label='NO', molecule=[SimpleNamespace(
        get_molecular_weight=lambda: .030006258748471737)])
    library = SimpleNamespace(entries={1: SimpleNamespace(item=SimpleNamespace(
        reactants=[species], products=[]))})
    frozen = {'states': {'NO(X)': {'paschen': 'ground', 'rmg_species': 'NO'}}}
    mass = adapter._gas_mass_from_library('NO', frozen, library)
    assert mass == pytest.approx(.030006258748471737 / 6.02214076e23)


def test_prepare_mixture_spec_uses_real_59_level_source_manifest(
        conversion_inputs, tmp_path, monkeypatch):
    import importlib.util
    from rmgpy.tools.eedf.integrity import _resolved_properties
    database = os.environ.get('EEDF_DATABASE_ROOT')
    if not database:
        pytest.skip('requires EEDF_DATABASE_ROOT for sourced N2 ladder')
    checker_path = Path(database) / 'scripts/check_n2_vibrational_data.py'
    source = importlib.util.spec_from_file_location('i369_n2_checker', checker_path)
    checker = importlib.util.module_from_spec(source)
    source.loader.exec_module(checker)
    data_root = os.environ.get('EEDF_DATA_ROOT')
    if not data_root:
        pytest.skip('requires EEDF_DATA_ROOT for sourced N2 cross sections')
    manifest = checker.resolve_generator_manifest(data_root)

    spec, ar_map, n2_map, libraries = _two_gas_fixture(conversion_inputs, tmp_path)
    monkeypatch.setattr(adapter, 'load_spec', lambda _: copy.deepcopy(spec))
    monkeypatch.setattr(adapter, '_load_library',
                        lambda _, __, label: (libraries[label], {'fixture': label}))
    channels = json.loads(Path(spec['channel_map']['path']).read_text())
    report = {'outcome': 'converted', 'selected_channel_count': len(channels),
              'missing_owners': [], 'shared_owners': [], 'branch_families': []}
    monkeypatch.setattr(adapter, 'convert_channel_maps',
                        lambda *args: (channels, copy.deepcopy(report)))
    map_paths = []
    for gas, frozen, label in (('Ar', ar_map, 'PlasmaArgon'),
                               ('N2', n2_map, 'PlasmaNitrogen')):
        path = tmp_path / f'{gas}-map.yaml'
        path.write_text(yaml.safe_dump(frozen))
        map_paths.append({'gas': gas, 'path': str(path), 'sha256': _file_hash(path),
                          'library': label})
    prepared, converted, conversion_report, mixture = adapter.prepare_mixture_spec(
        'unused', map_paths, database, 'mix-Ar-N2', manifest)
    assert len(mixture['vibrational_ladders'][0]['levels']) == 59
    assert mixture['vibrational_ladders'][0]['K'] == 10
    assert 'N2(X,v=58)' in mixture['vibrational_ladders'][0]['zero_population_tail']
    assert 'N2(X,v=58) = 0' not in mixture['state_properties']['population']
    mixture.update(working_conditions={}, solver_options={}, Tg_K=300., P_Pa=101325.)
    from rmgpy.tools.eedf.loki import _setup_text
    setup = yaml.safe_load(_setup_text(
        mixture, {'x_N2': .25, 'Tv_K': 1000.}, [5.], 'test-output'))
    for field in ('population', 'energy', 'statisticalWeight'):
        entries = setup['electronKinetics']['stateProperties'][field]
        names = [entry.split(' = ', 1)[0] for entry in entries]
        assert not any(name.startswith('N2(X,v=') and
                       int(name.rsplit('=', 1)[1][:-1]) > 10 for name in names)
    destination = tmp_path / 'real-ladder-mixture'
    adapter._write_outputs(destination, conversion_report, prepared, converted, mixture)
    published = json.loads((destination / 'ar-n2-spec.json').read_text())
    _, _, fractions, populations = _resolved_properties(
        published, {'x_N2': .25, 'Tv_K': 1000.})
    assert fractions == {'Ar': .75, 'N2': .25}
    assert len(populations) == 61
    assert all(populations[f'N2(X,v={v})'] == 0 for v in range(11, 59))


@pytest.mark.parametrize('failure', ['wrong_library', 'missing_owner', 'cross_gas_duplicate'])
def test_multigas_owner_failures_refuse(conversion_inputs, tmp_path, failure):
    spec, ar_map, n2_map, libraries = _two_gas_fixture(conversion_inputs, tmp_path)
    owner_library = 'PlasmaNitrogen'
    if failure == 'wrong_library':
        owner_library = 'PlasmaArgon'
    elif failure == 'missing_owner':
        n2_map['channels'][0]['rmg'] = None
    else:
        owner_library = 'PlasmaArgon'
        n2_map['channels'][0]['rmg'].update(library='PlasmaArgon', index=87,
                                             reaction='Ar + e- => Ars + e-')
    with pytest.raises(SpecError):
        adapter.convert_channel_maps(
            spec, {'Ar': ar_map, 'N2': n2_map}, libraries, 'mix-Ar-N2',
            library_by_gas={'Ar': 'PlasmaArgon', 'N2': owner_library})


@pytest.mark.parametrize('problem', ['normal', 'swapped_libraries', 'undeclared_map',
                                     'target_map_mismatch', 'missing_gas_map'])
def test_mixture_binds_physical_gases_to_maps_and_libraries(
        conversion_inputs, tmp_path, problem):
    spec, ar_map, n2_map, libraries = _two_gas_fixture(conversion_inputs, tmp_path)
    maps = {'Ar': ar_map, 'N2': n2_map}
    configured = {'Ar': 'PlasmaArgon', 'N2': 'PlasmaNitrogen'}
    if problem == 'swapped_libraries':
        configured = {'Ar': 'PlasmaNitrogen', 'N2': 'PlasmaArgon'}
    elif problem == 'undeclared_map':
        maps['He'] = maps.pop('N2')
        configured['He'] = configured.pop('N2')
    elif problem == 'target_map_mismatch':
        n2_map['channels'][0]['loki_process'] = ar_map['channels'][0]['loki_process']
    elif problem == 'missing_gas_map':
        maps.pop('N2')
        configured.pop('N2')
    if problem == 'normal':
        channels, _ = adapter.convert_channel_maps(
            spec, maps, libraries, 'mix-Ar-N2', library_by_gas=configured)
        assert len(channels) == 2
    else:
        expected = {
            'swapped_libraries': 'configured owner library mismatch for gas Ar',
            'undeclared_map': 'channel map gas is not declared: He',
            'target_map_mismatch': 'channel target gas differs from map gas N2',
            'missing_gas_map': 'missing channel map for gas N2'}[problem]
        with pytest.raises(SpecError, match=expected):
            adapter.convert_channel_maps(
                spec, maps, libraries, 'mix-Ar-N2', library_by_gas=configured)


def test_effective_source_maps_to_one_unowned_elastic_heating_row(
        conversion_inputs, tmp_path):
    """Fixture excerpt from LoKI input/Nitrogen/N2_LXCat.txt (SHA256 d7299f42250e6a8897f93f81acb9d2a4eaafc0bea684bddfa92d79ad43647e1a)."""
    from scipy.constants import electron_mass
    from rmgpy.tools.eedf.moments import _distribution_moments
    from rmgpy.tools.eedf.validation import _channel_power_checks
    from rmgpy.tools.eedf.loki import _partition_effective_elastic_power

    spec, ar_map, n2_map, libraries = _two_gas_fixture(conversion_inputs, tmp_path)
    fixture = Path(__file__).parent / 'fixtures/n2-effective-minimal.txt'
    effective_path = tmp_path / 'N2_effective_LXCat.txt'
    effective_path.write_bytes(fixture.read_bytes())
    mass_path = tmp_path / 'masses.txt'
    from scipy.constants import atomic_mass
    n2_mass = 28.0134 * atomic_mass
    mass_path.write_text(f'Ar {39.948 * atomic_mass}\nN2 {n2_mass}\n')
    spec['state_properties'].update(
        energy=['N2(X,v=0) = 0', 'N2(X,v=1) = .29'],
        statisticalWeight=['N2(X,v=0) = 1', 'N2(X,v=1) = 1'])
    spec['input_files']['N2_effective_LXCat.txt'] = {
        'path': str(effective_path), 'sha256': _file_hash(effective_path), 'kind': 'cross_section'}
    spec['input_files']['masses.txt'] = {
        'path': str(mass_path), 'sha256': _file_hash(mass_path), 'kind': 'gas_property'}
    spec['gas_properties']['mass'] = 'masses.txt'
    from rmgpy.tools.eedf.channels import _classic_sections
    effective = _classic_sections(str(effective_path), _file_hash(effective_path))[0]
    effective.update(classification='B', target_fraction=.5, product_fraction=.5,
                     sigma_max_m2=1.62e-20,
                     flux_group=effective['description'], reaction=None)
    physical = json.loads(Path(spec['channel_map']['path']).read_text())
    physical.append(effective)
    channel_path = Path(spec['channel_map']['path'])
    channel_path.write_text(json.dumps(physical))
    spec['channel_map']['sha256'] = _file_hash(channel_path)
    n2_map['loki_files']['N2_effective_LXCat.txt'] = _file_hash(effective_path)
    n2_map['class_counts'] = {'A': 1, 'B': 1, 'C': 0, 'D': 0}
    n2_map['channels'].append({
        'id': 10, 'loki_process': effective['description'],
        'lxcat_file': 'N2_effective_LXCat.txt', 'threshold_eV': 0., 'class': 'B', 'rmg': None})

    channels, _ = adapter.convert_channel_maps(
        spec, {'Ar': ar_map, 'N2': n2_map}, libraries, 'mix-Ar-N2',
        library_by_gas={'Ar': 'PlasmaArgon', 'N2': 'PlasmaNitrogen'})
    effective_row = next(row for row in channels if row['kind'] == 'effective')
    elastic_row = next(row for row in channels if row['kind'] == 'elastic')
    assert effective_row['description'] == effective['description']
    assert elastic_row['description'].endswith(', Elastic')
    assert elastic_row['classification'] == 'B'
    assert elastic_row['reaction'] is None
    assert elastic_row['mass_ratio'] == pytest.approx(electron_mass / n2_mass, rel=1e-12, abs=0.)

    energy = [.0005, .0015, .0025, .0035]
    edges = [0., .001, .002, .003, .004]
    row = {'energy_eV': energy, 'energy_edges_eV': edges,
           'f0': [0.2, 0.15, 0.08, 0.02], 'target_fractions': [.5] * len(channels),
           'product_fractions': [.5] * len(channels), 'gas_temperature_K': 300.}
    moments = _distribution_moments(row, channels)
    eff_index = next(i for i, item in enumerate(channels) if item['kind'] == 'effective')
    elastic_index = next(i for i, item in enumerate(channels) if item['kind'] == 'elastic')
    assert moments['channel_power'][eff_index] == 0.
    elastic_power = moments['channel_power'][elastic_index]
    assert elastic_power != 0.
    powers = list(moments['channel_power'])
    elastic_indices = [i for i, item in enumerate(channels) if item['kind'] == 'elastic']
    derived_index = next(i for i in elastic_indices
                         if channels[i]['description'] ==
                         effective_row['description'].replace(', Effective', ', Elastic'))
    solver_elastic_net = sum(powers[i] for i in elastic_indices) + abs(elastic_power) * .2
    _partition_effective_elastic_power(channels, powers,
                                      {'elastic_gain': 0., 'elastic_loss': -solver_elastic_net})
    assert powers[eff_index] == 0.
    assert powers[derived_index] == pytest.approx(
        solver_elastic_net - sum(powers[i] for i in elastic_indices if i != derived_index), rel=1e-12, abs=0.)
    assert sum(powers[i] for i in elastic_indices) == pytest.approx(solver_elastic_net, rel=1e-12, abs=0.)
    powers_by_kind = dict(zip((item['kind'] for item in channels), moments['channel_power']))
    assert powers_by_kind['effective'] == 0.
    power_row = {'channel_power': powers,
                 'power_groups': {'elastic_gain': 0., 'elastic_loss': -solver_elastic_net,
                                  'excitation_loss': 0., 'excitation_gain': 0.,
                                  'vibrational_loss': 0., 'vibrational_gain': 0.,
                                  'rotational_loss': 0., 'rotational_gain': 0.,
                                  'ionization': 0., 'attachment': 0., 'field': abs(solver_elastic_net)},}
    tolerances = {'channel_power_sum': {'atol': abs(solver_elastic_net) * 1e-12, 'rtol': 1e-10}}
    checks = _channel_power_checks(power_row, {
        'channel_map': channels, 'tolerances': tolerances,
        'floors': {'relative_power_share': .01, 'absolute_power_share': 1e-12}})
    elastic_check = next(check for check in checks if check['quantity'] == 'power.sum.elastic')
    assert elastic_check['passed']
    assert elastic_power != pytest.approx(-power_row['power_groups']['elastic_loss'], rel=1e-12, abs=0.)
    power_row['channel_power'][derived_index] += abs(solver_elastic_net) * .01
    checks = _channel_power_checks(power_row, {
        'channel_map': channels, 'tolerances': tolerances,
        'floors': {'relative_power_share': .01, 'absolute_power_share': 1e-12}})
    assert not next(check for check in checks if check['quantity'] == 'power.sum.elastic')['passed']


def test_multiple_effective_sources_share_the_solver_elastic_balance():
    from rmgpy.tools.eedf.loki import _partition_effective_elastic_power

    channels = [
        {'kind': 'effective', 'description': 'e + N2(X) -> e + N2(X), Effective'},
        {'kind': 'effective', 'description': 'e + NO(X) -> e + NO(X), Effective'},
        {'kind': 'elastic', 'description': 'e + Ar(1S0) -> e + Ar(1S0), Elastic'},
        {'kind': 'elastic', 'description': 'e + N2(X) -> e + N2(X), Elastic'},
        {'kind': 'elastic', 'description': 'e + NO(X) -> e + NO(X), Elastic'},
    ]
    powers = [0., 0., 3., 4., 5.]
    _partition_effective_elastic_power(
        channels, powers, {'elastic_gain': 0., 'elastic_loss': -15.})
    assert powers[3] == 4.
    assert powers[4] == 8.
    assert sum(powers[index] for index in (2, 3, 4)) == 15.


def test_irreversible_excitation_group_has_no_product_population_requirement():
    from rmgpy.tools.eedf.integrity import _channel_fractions
    channel = {'description': 'e + N2(X) -> e + N2(A3Su+,v=0-4), Excitation',
               'kind': 'excitation', 'target_fraction': 'state:N2(X)',
               'product_fraction': 0.}
    targets, products = _channel_fractions([channel], {'N2': 1.},
                                           {'N2(X)': 1.}, {})
    assert targets == [1.]
    assert products == [0.]


def test_mixture_spec_writer_emits_composition_and_tv_ladder(conversion_inputs, monkeypatch):
    import rmgpy.tools.eedf.vibrational as vibrational
    from rmgpy.tools.eedf.integrity import _resolved_properties
    spec, _, _ = conversion_inputs
    levels = [{'v': 0, 'loki_state': 'N2(X,v=0)'},
              {'v': 1, 'loki_state': 'N2(X,v=1)'}]
    levels[0].update(energy_eV=0., statistical_weight=1.)
    levels[1].update(energy_eV=.29, statistical_weight=1.)
    ladder = {'axis': 'Tv_K', 'parent': 'N2(X)', 'K': 1, 'levels': levels,
              'source_domain': {'initial_levels': list(range(30)), 'final_v': 58}}
    monkeypatch.setattr(vibrational, 'ladder_coupling', lambda *args, **kwargs: ladder)
    mixture = adapter.build_mixture_spec(spec, {}, [])
    assert mixture['axes']['Tv_K'] == [1000., 2000.]
    assert mixture['gas_properties']['fraction'] == ['Ar = 1-{x_N2}', 'N2 = {x_N2}']
    assert 'N2(X,v=1) = boltzmann:Tv_K:1' in mixture['state_properties']['population']
    assert mixture['vibrational_ladders'] == [ladder]
    _, _, fractions, populations = _resolved_properties(mixture, {'x_N2': .25, 'Tv_K': 1000.})
    assert fractions == {'Ar': .75, 'N2': .25}
    assert populations['N2(X,v=0)'] + populations['N2(X,v=1)'] == pytest.approx(1.)


@pytest.mark.parametrize('problem', ['absent', 'missing_index', 'wrong_products', 'wrong_direction'])
def test_refuses_unowned_a_without_demotion(conversion_inputs, problem):
    spec, frozen, library = conversion_inputs
    owner = frozen['channels'][0]['rmg']
    if problem == 'absent':
        owner.update(status='absent', index=None)
        frozen['class_a_gaps'] = {'[absent] Ar + e- => Ars + e-': [2]}
    elif problem == 'missing_index':
        owner['index'] = None
    elif problem == 'wrong_products':
        library.entries[87].item.products = [SimpleNamespace(label='Arp')]
    else:
        library.entries[87].item.reversible = True
    with pytest.raises(adapter.MissingOwners) as caught:
        adapter.convert_channel_map(spec, frozen, library, 'test-Ar')
    report = caught.value.report
    assert report['outcome'] == 'refused'
    assert report['class_counts'] == {'A': 1}
    assert report['missing_owners'] == [{
        'id': 2, 'process': 'e + Ar(1S0) -> e + Ar(3P2), Excitation',
        'classification': 'A', 'products': 'e + Ar(3P2)', 'threshold_eV': 11.548,
        'class_a_gaps_status': owner['status'], 'requested_owner': owner,
        'reason': {'absent': 'absent', 'missing_index': 'requested library entry is absent',
                   'wrong_products': 'requested reaction differs from loaded library entry',
                   'wrong_direction': 'requested reaction differs from loaded library entry'}[problem]}]


@pytest.mark.parametrize('problem', [
    'schema', 'class', 'reaction', 'gaps', 'duplicates', 'threshold', 'file_hash', 'coverage'])
def test_frozen_map_fails_closed(conversion_inputs, problem):
    spec, frozen, library = conversion_inputs
    row = frozen['channels'][0]
    if problem == 'schema':
        frozen['schema'] = 'unknown'
    elif problem == 'class':
        del row['class']
    elif problem == 'reaction':
        del row['rmg']
    elif problem == 'gaps':
        row['rmg'].update(status='absent', index=None)
    elif problem == 'duplicates':
        frozen['channels'].append(copy.deepcopy(row))
    elif problem == 'threshold':
        row['threshold_eV'] = 5.
    elif problem == 'file_hash':
        frozen['loki_files']['Ar_LXCat.txt'] = '0' * 64
    else:
        row['loki_process'] = 'e + Ar(1S0) -> e + Ar(3P0), Excitation'
    with pytest.raises((SpecError, FingerprintMismatch)):
        adapter.convert_channel_map(spec, frozen, library, 'test-Ar')


def test_old_all_b_map_remains_loadable_without_frozen_yaml(conversion_inputs):
    from rmgpy.tools.eedf.channels import _load_channel_map, _validate_physical_map
    spec, _, _ = conversion_inputs
    channels = _load_channel_map(spec['channel_map'])
    _validate_physical_map(channels, spec, {})
    assert channels[0]['classification'] == 'B'
    assert channels[0]['reaction'] is None


def test_refusal_cli_writes_only_report(conversion_inputs, tmp_path, monkeypatch):
    spec, frozen, library = conversion_inputs
    frozen['channels'][0]['rmg'].update(status='absent', index=None)
    frozen['class_a_gaps'] = {'[absent] Ar + e- => Ars + e-': [2]}

    def prepare(*args):
        return adapter.convert_channel_map(spec, frozen, library, 'test-Ar')

    monkeypatch.setattr(adapter, 'prepare_spec', prepare)
    output = tmp_path / 'refusal'
    argv = ['--base-spec', 'unused', '--frozen-map', 'unused',
            '--frozen-map-sha256', 'unused', '--database', 'unused',
            '--collision-set', 'test-Ar', '--output-directory', str(output)]
    assert adapter._main(argv) == 2
    assert [p.name for p in output.iterdir()] == ['conversion-report.json']
    report = json.loads((output / 'conversion-report.json').read_text())
    assert [row['id'] for row in report['missing_owners']] == [2]
    original = (output / 'conversion-report.json').read_bytes()
    assert adapter._main(argv) == 2
    assert (output / 'conversion-report.json').read_bytes() == original


@pytest.mark.parametrize('problem', ['missing', 'malformed', 'hash', 'duplicate_key'])
def test_prepare_rejects_yaml_before_library_load(conversion_inputs, tmp_path, monkeypatch, problem):
    spec, _, _ = conversion_inputs
    monkeypatch.setattr(adapter, 'load_spec', lambda path: spec)
    monkeypatch.setattr(adapter, '_load_library', lambda *args: pytest.fail('library loaded'))
    path = tmp_path / 'map.yaml'
    path.write_text('class: A\nclass: B\n' if problem == 'duplicate_key' else 'channels: [')
    pin = _file_hash(path)
    if problem == 'missing':
        path.unlink()
    elif problem == 'hash':
        pin = '0' * 64
    with pytest.raises((SpecError, FingerprintMismatch)):
        adapter.prepare_spec('unused', path, pin, 'unused', 'test-Ar')


def test_publication_failure_preserves_report_and_never_overwrites(conversion_inputs, tmp_path, monkeypatch):
    spec, frozen, library = conversion_inputs
    channels, report = adapter.convert_channel_map(spec, frozen, library, 'test-Ar')
    destination = tmp_path / 'interrupted'
    original_open = Path.open

    def fail_spec(path, *args, **kwargs):
        if path.name == 'pure-ar-spec.json':
            raise OSError('injected publication failure')
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, 'open', fail_spec)
    with pytest.raises(OSError, match='injected'):
        adapter._write_outputs(destination, report, spec, channels)
    assert sorted(path.name for path in destination.iterdir()) == [
        'channel-map.json', 'conversion-report.json']
    before = {path.name: path.read_bytes() for path in destination.iterdir()}
    with pytest.raises(FileExistsError):
        adapter._write_outputs(destination, report, spec, channels)
    assert before == {path.name: path.read_bytes() for path in destination.iterdir()}


@pytest.fixture
def real_inputs():
    # Do not opt the generator's expensive real-solver tests in just to inspect metadata.
    names = ['EEDF_CHANNEL_MAP_BASE_SPEC', 'EEDF_FROZEN_CHANNEL_MAP', 'EEDF_FROZEN_CHANNEL_MAP_SHA256',
             'EEDF_DATABASE_ROOT']
    missing = [name for name in names if name not in os.environ]
    if missing:
        pytest.skip('real conversion inputs require ' + ', '.join(missing))
    return [os.environ[name] for name in names]


def test_real_map_literal_missing_owners(real_inputs):
    spec, frozen, digest, database = real_inputs
    with pytest.raises(adapter.MissingOwners) as caught:
        adapter.prepare_spec(spec, frozen, digest, database, 'IST-Lisbon-Ar-2021')
    report = caught.value.report
    assert report['class_counts'] == {'A': 38, 'B': 1}
    assert report['selected_channel_count'] == 39
    assert [row['id'] for row in report['missing_owners']] == [
        3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
        21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38]
    assert all(row['class_a_gaps_status'] == 'absent' for row in report['missing_owners'])
    assert report['shared_owners'] == [{'library': 'PlasmaArgon', 'index': 87, 'ids': [2, 4]}]


@pytest.mark.parametrize('channel_id, classification, index, threshold, count, energy, sigma', [
    (1, 'B', None, 0., 60, 0., 7.5e-20),
    (2, 'A', 87, 11.548, 27, 13., 3.237651e-22),
    (39, 'A', 86, 15.76, 28, 16., 2.024e-22),
    (4, 'A', 87, 11.723, 27, 13., 2.383485e-23),
])
def test_real_named_channels_have_literal_owners_and_cross_sections(
        real_inputs, channel_id, classification, index, threshold, count, energy, sigma):
    spec_path, frozen_path, digest, _ = real_inputs
    frozen = yaml.safe_load(adapter._read_pinned_bytes(frozen_path, digest))
    spec = load_spec(spec_path)
    original = adapter._load_channel_map(spec['channel_map'])
    row = next(row for row in frozen['channels'] if row['id'] == channel_id)
    physical = next(row for row in original if row['description'] == frozen['channels'][channel_id-1]['loki_process'])
    result = adapter._convert_channel(physical, row, 'IST-Lisbon-Ar-2021')
    assert result['classification'] == classification
    assert result['threshold_eV'] == threshold
    assert len(result['cross_section']['energy_eV']) == count
    point = result['cross_section']['energy_eV'].index(energy)
    assert result['cross_section']['sigma_m2'][point] == sigma
    assert result['cross_section'] == physical['cross_section']
    if index is None:
        assert result['reaction'] is None
    else:
        assert result['reaction']['library'] == 'PlasmaArgon'
        assert result['reaction']['index'] == index
        processes = {
            2: 'e + Ar(1S0) -> e + Ar(3P2), Excitation',
            4: 'e + Ar(1S0) -> e + Ar(3P0), Excitation',
            39: 'e + Ar(1S0) -> e + e + Ar(+,gnd), Ionization'}
        assert result['reaction']['repr'] == (
            'EEDFChannel(process=' + repr(processes[channel_id]) +
            ", collision_set='IST-Lisbon-Ar-2021', side='ine')")


@pytest.fixture
def cascade_conversion(conversion_inputs, branch_fixture):
    spec, frozen, _ = conversion_inputs
    snapshot, reactions, family = branch_fixture
    path = Path(spec['input_files']['Ar_LXCat.txt']['path'])
    path.write_text('PARAM.: E = 12.907 eV\n'
                    '[Ar(1S0) + e -> Ar(4p[1/2]1) + e, Excitation]\n'
                    '-----\n12.907 0\n13 3.237651e-22\n-----\n')
    spec['input_files']['Ar_LXCat.txt']['sha256'] = _file_hash(path)
    frozen['loki_files']['Ar_LXCat.txt'] = _file_hash(path)
    path = Path(spec['channel_map']['path'])
    physical = json.loads(path.read_text())[0]
    physical.update(description=PROCESS, threshold_eV=12.907, flux_group=PROCESS)
    physical['cross_section']['energy_eV'][0] = 12.907
    path.write_text(json.dumps([physical]))
    spec['channel_map']['sha256'] = _file_hash(path)
    row = frozen['channels'][0]
    row.update(id=6, loki_process=PROCESS, threshold_eV=12.907)
    row['rmg'].update(status='absent', index=None)
    frozen['class_a_gaps'] = {'[absent] cascade source': [6]}
    library = SimpleNamespace(entries={r.entry.index: r.entry for r in reactions})
    return spec, frozen, library, snapshot, family


def test_complete_family_converts_from_absent_placeholder(cascade_conversion, tmp_path):
    spec, frozen, library, snapshot, family = cascade_conversion
    converted, report = adapter.convert_channel_map(spec, frozen, library, 'test-Ar', snapshot)
    assert len(converted) == 1
    assert converted[0]['reaction'] == family
    assert report['branch_families'] == [{'id': 6, 'family': family}]
    assert report['missing_owners'] == []
    adapter._write_outputs(tmp_path / 'complete', report, spec, converted)
    published = json.loads((tmp_path / 'complete/channel-map.json').read_text())
    assert published == converted


@pytest.mark.parametrize('problem', ['missing', 'extra', 'unweighted', 'source', 'unresolved',
                                     'wrong_product', 'stale_yield', 'ordinary_overlap'])
def test_converter_refuses_incomplete_or_unmatched_family(cascade_conversion, problem):
    from rmgpy.data.base import Entry
    from rmgpy.kinetics import EEDFChannel
    spec, frozen, library, snapshot, _ = cascade_conversion
    entry = next(iter(library.entries.values()))
    if problem == 'missing':
        library.entries.pop(entry.index)
    elif problem == 'extra':
        library.entries[9000] = copy.deepcopy(entry)
    elif problem == 'unweighted':
        library.entries[9000] = Entry(index=9000, item=entry.item,
                                      data=EEDFChannel(process=PROCESS, collision_set='test-Ar'))
    elif problem == 'source':
        snapshot = None
    elif problem == 'unresolved':
        snapshot['channels'][0].update(status='unresolved', product_yields=None)
    elif problem == 'wrong_product':
        entry.item.products[1].molecule[0].electronic_state = '3P1'
    elif problem == 'stale_yield':
        snapshot['channels'][0]['product_yields']['Ar(1s5)'] += .001
        snapshot['channels'][0]['product_yields']['Ar(1s4)'] -= .001
    else:
        frozen['channels'][0]['rmg'].update(status='present', index=87)
        frozen['class_a_gaps'] = {}
    with pytest.raises(adapter.MissingOwners) as caught:
        adapter.convert_channel_map(spec, frozen, library, 'test-Ar', snapshot)
    assert caught.value.report['branch_families'] == []
    assert [row['id'] for row in caught.value.report['missing_owners']] == [6]
    assert caught.value.report['missing_owners'][0]['reason'] not in ('absent', '')


def test_map_family_requires_matching_physics(cascade_conversion):
    spec, frozen, library, snapshot, _ = cascade_conversion
    channels, _ = adapter.convert_channel_map(spec, frozen, library, 'test-Ar', snapshot)
    channels[0]['threshold_eV'] = 11.
    with pytest.raises(SpecError, match='collision identity'):
        adapter._validate_map(channels)
    channels[0]['threshold_eV'] = 12.907
    channels[0]['product_fraction'] = .1
    with pytest.raises(SpecError, match='forward excitation'):
        adapter._validate_map(channels)


def test_converter_unweighted_duplicate_refuses(conversion_inputs):
    from rmgpy.kinetics import EEDFChannel
    spec, frozen, library = conversion_inputs
    entry = library.entries[87]
    entry.data = EEDFChannel(process=frozen['channels'][0]['loki_process'], collision_set='test-Ar')
    library.entries[88] = copy.deepcopy(entry)
    with pytest.raises(adapter.MissingOwners, match='lack owning') as caught:
        adapter.convert_channel_map(spec, frozen, library, 'test-Ar')
    assert 'double-count one tabulated rate' in caught.value.report['missing_owners'][0]['reason']


def test_prepare_family_reads_pinned_sources(cascade_conversion, source_files, tmp_path, monkeypatch):
    spec, frozen, library, _, family = cascade_conversion
    frozen_path = tmp_path / 'frozen.yaml'
    frozen_path.write_text(yaml.safe_dump(frozen))
    monkeypatch.setattr(adapter, 'load_spec', lambda path: spec)
    monkeypatch.setattr(adapter, '_load_library', lambda *args: (library, {'fixture': 'digest'}))
    yields, source, _ = source_files
    _, converted, report = adapter.prepare_spec(
        'unused', frozen_path, _file_hash(frozen_path), tmp_path, 'test-Ar',
        yield_table=yields, radiative_source=source,
        yield_table_sha256=_file_hash(yields), source_sha256=_file_hash(source))
    assert converted[0]['reaction'] == family
    assert report['provenance']['cascade_sources']['yield_table_sha256'] == _file_hash(yields)
    with pytest.raises(FingerprintMismatch, match='source hash'):
        adapter.prepare_spec('unused', frozen_path, _file_hash(frozen_path), tmp_path, 'test-Ar',
                             yield_table=yields, radiative_source=source, yield_table_sha256='0' * 64)


def test_mixture_cli_loads_pinned_cascade_sources_and_converts_family(
        conversion_inputs, cascade_conversion, source_files, tmp_path, monkeypatch):
    import importlib.util
    spec, ar_map, _, cascade_snapshot, _ = cascade_conversion
    ar_library = SimpleNamespace(entries={
        reaction.entry.index: reaction.entry
        for reaction in branch_reactions(cascade_snapshot, 6, 'mix-Ar-N2')})
    vib_spec, _, n2_map, libraries = _two_gas_fixture(conversion_inputs, tmp_path)
    spec = copy.deepcopy(spec)
    spec['input_files']['N2_LXCat.txt'] = copy.deepcopy(vib_spec['input_files']['N2_LXCat.txt'])
    spec['axes']['x_N2'] = [.5]
    spec['gas_properties']['fraction'] = ['Ar = 1-{x_N2}', 'N2 = {x_N2}']
    spec['state_properties'] = copy.deepcopy(vib_spec['state_properties'])
    physical = [row for row in json.loads(Path(spec['channel_map']['path']).read_text())
                if row['description'] == PROCESS]
    physical.extend(row for row in json.loads(Path(vib_spec['channel_map']['path']).read_text())
                    if row['description'].startswith('e + N2('))
    channel_path = Path(spec['channel_map']['path'])
    channel_path.write_text(json.dumps(physical))
    spec['channel_map']['sha256'] = _file_hash(channel_path)
    loaded_libraries = {'PlasmaArgon': ar_library,
                        'PlasmaNitrogen': libraries['PlasmaNitrogen']}
    monkeypatch.setattr(adapter, 'load_spec', lambda _: copy.deepcopy(spec))
    monkeypatch.setattr(adapter, '_load_library',
                        lambda _, __, label: (loaded_libraries[label], {'fixture': label}))

    database = os.environ.get('EEDF_DATABASE_ROOT')
    data_root = os.environ.get('EEDF_DATA_ROOT')
    if not database or not data_root:
        pytest.skip('requires EEDF_DATABASE_ROOT and EEDF_DATA_ROOT for the real ladder')
    checker_path = Path(database) / 'scripts/check_n2_vibrational_data.py'
    checker_spec = importlib.util.spec_from_file_location('i369_n2_checker_cli', checker_path)
    checker = importlib.util.module_from_spec(checker_spec)
    checker_spec.loader.exec_module(checker)
    manifest = checker.resolve_generator_manifest(data_root)
    manifest_path = tmp_path / 'ladder-manifest.json'
    manifest_path.write_text(json.dumps(manifest))
    map_options = []
    for gas, frozen, label in (('Ar', ar_map, 'PlasmaArgon'),
                               ('N2', n2_map, 'PlasmaNitrogen')):
        path = tmp_path / f'{gas}-cascade-map.yaml'
        path.write_text(yaml.safe_dump(frozen))
        map_options.extend(['--gas-map', gas, str(path), _file_hash(path), label])
    yields, radiative_source, _ = source_files
    output = tmp_path / 'mixture-cli'
    argv = ['--base-spec', 'unused', *map_options,
            '--ladder-manifest', str(manifest_path), '--database', str(tmp_path),
            '--collision-set', 'mix-Ar-N2', '--output-directory', str(output),
            '--yield-table', str(yields), '--radiative-source', str(radiative_source),
            '--yield-table-sha256', _file_hash(yields),
            '--source-sha256', _file_hash(radiative_source)]
    assert adapter._main(argv) == 0
    published = json.loads((output / 'conversion-report.json').read_text())
    assert [row['id'] for row in published['branch_families']] == [6]
    assert published['branch_families'][0]['family']['collision_set'] == 'mix-Ar-N2'
    assert published['provenance']['cascade_sources'] == {
        'yield_table_sha256': _file_hash(yields),
        'source_sha256': _file_hash(radiative_source),
        'radiative_model': published['provenance']['cascade_sources']['radiative_model']}
    assert (output / 'ar-n2-spec.json').is_file()


def test_i353_fixture_families_are_reported(tmp_path, monkeypatch):
    """Complete cascades convert; an in-memory missing owner is reported exactly."""
    from rmgpy.tools.eedf.ownership import load_cascade_yields
    database = os.environ.get('EEDF_DATABASE_ROOT')
    spec_path = os.environ.get('EEDF_CHANNEL_MAP_BASE_SPEC')
    if not database or not spec_path:
        pytest.skip('requires EEDF_DATABASE_ROOT and EEDF_CHANNEL_MAP_BASE_SPEC')
    root = Path(database) / 'input/eedf'
    snapshot = load_cascade_yields(root / 'argon_cascade_yields.yaml', root / 'argon_radiative_data.yaml')
    assert snapshot['yield_table_sha256'] == _file_hash(root / 'argon_cascade_yields.yaml')
    assert snapshot['channels'] == database_cascade_channels(Path(database))
    complete = [row['id'] for row in snapshot['channels'] if row['status'] == 'complete_for_E1']
    unresolved = [row['id'] for row in snapshot['channels'] if row['status'] == 'unresolved']
    assert not unresolved
    library, hashes = adapter._load_library(database, tmp_path)
    collision_set = 'IST-Lisbon-Ar-2021'
    for channel_id in complete:
        row = next(row for row in snapshot['channels'] if row['id'] == channel_id)
        if any(getattr(entry.data, 'process', None) == row['process']
               for entry in library.entries.values()):
            continue
        for reaction in branch_reactions(snapshot, channel_id, collision_set):
            assert reaction.entry.index not in library.entries
            library.entries[reaction.entry.index] = reaction.entry
    monkeypatch.setattr(adapter, '_load_library', lambda *args: (library, hashes))
    frozen_path = root / 'argon_channel_map.yaml'
    _, _, report = adapter.prepare_spec(
        spec_path, frozen_path, _file_hash(frozen_path), database, collision_set)
    assert [row['id'] for row in report['branch_families']] == list(range(6, 39))
    assert report['missing_owners'] == []
    assert report['shared_owners'] == []

    missing_id = 16
    missing_process = next(row['process'] for row in snapshot['channels'] if row['id'] == missing_id)
    incomplete_library = copy.deepcopy(library)
    removed = [index for index, entry in incomplete_library.entries.items()
               if getattr(entry.data, 'process', None) == missing_process]
    assert removed
    for index in removed:
        del incomplete_library.entries[index]
    monkeypatch.setattr(adapter, '_load_library', lambda *args: (incomplete_library, hashes))
    with pytest.raises(adapter.MissingOwners) as caught:
        adapter.prepare_spec(spec_path, frozen_path, _file_hash(frozen_path), database, collision_set)
    report = caught.value.report
    assert [row['id'] for row in report['branch_families']] == [i for i in range(6, 39) if i != missing_id]
    assert [row['id'] for row in report['missing_owners']] == [missing_id]
    assert report['shared_owners'] == []
    adapter._write_outputs(tmp_path / 'refusal', report)
    assert [path.name for path in (tmp_path / 'refusal').iterdir()] == ['conversion-report.json']


@pytest.mark.parametrize('location', ['configuration', 'frozen', 'lookup'])
def test_list_owner_label_refuses_by_name(conversion_inputs, tmp_path, location):
    spec, ar_map, n2_map, libraries = _two_gas_fixture(conversion_inputs, tmp_path)
    configured = {'Ar': 'PlasmaArgon', 'N2': 'PlasmaNitrogen'}
    if location == 'configuration':
        configured['Ar'] = ['PlasmaArgon']
    else:
        ar_map['channels'][0]['rmg']['library'] = ['PlasmaArgon']
    with pytest.raises(SpecError, match='owner library|reaction metadata'):
        if location == 'lookup':
            adapter._owner_problem(ar_map['channels'][0], libraries)
        else:
            adapter.convert_channel_maps(spec, {'Ar': ar_map, 'N2': n2_map},
                                         libraries, 'mix-Ar-N2', library_by_gas=configured)


@pytest.mark.parametrize('prepare', [False, True])
def test_mixture_builder_refuses_supplied_owner_d6(conversion_inputs, prepare):
    spec, _, _ = conversion_inputs
    owner = {'set_by': 'owner', 'top_fraction_limit': .01,
             'top_leak_share_limit': .02, 'K_step': 5}
    manifest = {'D6': owner}
    original = copy.deepcopy(manifest)
    with pytest.raises(SpecError, match='diagnostic D6 conflicts'):
        if prepare:
            adapter.prepare_mixture_spec('unused', [], 'unused', 'mixture', manifest)
        else:
            adapter.build_mixture_spec(spec, manifest, [])
    assert manifest == original
