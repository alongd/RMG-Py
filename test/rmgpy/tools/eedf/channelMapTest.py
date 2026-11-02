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
