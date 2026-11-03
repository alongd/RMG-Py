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
from branchYieldTest import PROCESS, branch_fixture, branch_reactions, source_files


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


def test_i353_fixture_families_are_reported(tmp_path, monkeypatch):
    """Use the real read-only library plus private branches, with no generation."""
    from rmgpy.tools.eedf.ownership import load_cascade_yields
    database = os.environ.get('EEDF_DATABASE_ROOT')
    spec_path = os.environ.get('EEDF_CHANNEL_MAP_BASE_SPEC')
    if not database or not spec_path:
        pytest.skip('requires EEDF_DATABASE_ROOT and EEDF_CHANNEL_MAP_BASE_SPEC')
    root = Path(database) / 'input/eedf'
    snapshot = load_cascade_yields(root / 'argon_cascade_yields.yaml', root / 'argon_radiative_data.yaml')
    # The live source contract has an independent literal byte pin.
    assert snapshot['yield_table_sha256'] == '200cd309f4915122117f37293bb74fcb785d0b5805d3d9868d3fe5eec2dac0f5'
    library, hashes = adapter._load_library(database, tmp_path)
    for channel_id in [6, 7, 8, 9, 10, 12, 13, 14, 15]:
        for reaction in branch_reactions(snapshot, channel_id):
            assert reaction.entry.index not in library.entries
            library.entries[reaction.entry.index] = reaction.entry
    monkeypatch.setattr(adapter, '_load_library', lambda *args: (library, hashes))
    frozen_path = root / 'argon_channel_map.yaml'
    with pytest.raises(adapter.MissingOwners) as caught:
        adapter.prepare_spec(spec_path, frozen_path, _file_hash(frozen_path), database, 'test-Ar')
    report = caught.value.report
    assert [row['id'] for row in report['branch_families']] == [6, 7, 8, 9, 10, 12, 13, 14, 15]
    assert [row['id'] for row in report['missing_owners']] == [
        11, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38]
    assert report['shared_owners'] == []
    adapter._write_outputs(tmp_path / 'refusal', report)
    assert [path.name for path in (tmp_path / 'refusal').iterdir()] == ['conversion-report.json']
