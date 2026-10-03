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
import runpy

import h5py
import numpy as np
import pytest

from rmgpy.solver.eedf import EEDFTable
from rmgpy.tools.eedf import artifact
from rmgpy.tools.eedf.generate import generate, model_inputs
from rmgpy.tools.eedf.integrity import physical_properties, resolved_properties
from rmgpy.tools.eedf.loki import LoKIDriver, enrich_row, setup_text
from rmgpy.tools.eedf.schema import FingerprintMismatch, SpecError, file_hash

helpers = runpy.run_path(str(Path(__file__).with_name('rework3Test.py')))
fixture_table = helpers['fixture_table']
request_and_channels = helpers['request_and_channels']
reversible_request = helpers['reversible_request']


def repin(path):
    manifest = json.loads((path / 'manifest.json').read_text())
    manifest['artifact_sha256'] = file_hash(path / 'table.h5')
    (path / 'manifest.json').write_text(json.dumps(manifest))
    return manifest['artifact_sha256']


@pytest.mark.parametrize('storage', ['external_raw', 'virtual', 'external_link',
                                    'soft_link', 'plugin', 'compact',
                                    'attribute_reference', 'dataset_reference',
                                    'hard_alias', 'hard_cycle'])
def test_r4_regression_storage_is_self_contained_before_reading(tmp_path, storage, monkeypatch):
    path, model = fixture_table(tmp_path / 'tables')
    with h5py.File(path / 'table.h5', 'r+') as h5:
        group = h5['branches/branch_0']
        values = group['k_ine'][...]
        if storage in ('external_raw', 'virtual', 'external_link', 'soft_link', 'plugin', 'compact'):
            del group['k_ine']
        if storage == 'external_raw':
            raw = tmp_path / 'rates.raw'
            group.create_dataset('k_ine', data=values, external=[(str(raw), 0, h5py.h5f.UNLIMITED)])
            pin_contents = raw.read_bytes()
            raw.write_bytes((values * 1e6).tobytes())
            assert raw.read_bytes() != pin_contents
        elif storage in ('virtual', 'external_link'):
            other = tmp_path / 'other.h5'
            with h5py.File(other, 'w') as source:
                source['rates'] = values
            if storage == 'virtual':
                layout = h5py.VirtualLayout(shape=values.shape, dtype=values.dtype)
                layout[:] = h5py.VirtualSource(str(other), 'rates', shape=values.shape)
                group.create_virtual_dataset('k_ine', layout)
            else:
                group['k_ine'] = h5py.ExternalLink(str(other), '/rates')
        elif storage == 'soft_link':
            h5['rates'] = values
            group['k_ine'] = h5py.SoftLink('/rates')
        elif storage in ('plugin', 'compact'):
            space = h5py.h5s.create_simple(values.shape)
            properties = h5py.h5p.create(h5py.h5p.DATASET_CREATE)
            if storage == 'plugin':
                properties.set_chunk(values.shape)
                # Optional unregistered filter: baseline silently reads unfiltered bytes.
                properties.set_filter(32001, h5py.h5z.FLAG_OPTIONAL)
            else:
                properties.set_layout(h5py.h5d.COMPACT)
            dataset = h5py.h5d.create(group.id, b'k_ine', h5py.h5t.py_create(values.dtype), space, dcpl=properties)
            dataset.write(h5py.h5s.ALL, h5py.h5s.ALL, values)
            dataset.close()
        elif storage == 'attribute_reference':
            h5.attrs['uninterpreted_reference'] = group['k_ine'].ref
        elif storage == 'dataset_reference':
            h5.create_dataset('uninterpreted_reference', data=[group['k_ine'].ref], dtype=h5py.ref_dtype)
        elif storage == 'hard_alias':
            h5['alias'] = group
        elif storage == 'hard_cycle':
            h5['cycle'] = h5['/']
    pin = repin(path)
    # All bytes in the HDF5 pin are unchanged when an external source changes.
    # Refuse structurally before interpreting any dataset or attribute value.
    reads = []
    for cls in (h5py.Dataset, h5py.AttributeManager):
        original = cls.__getitem__
        def read(self, key, original=original):
            reads.append(key)
            return original(self, key)
        monkeypatch.setattr(cls, '__getitem__', read)
    with pytest.raises(FingerprintMismatch, match='HDF5 storage'):
        EEDFTable.load(path, model, artifact_sha256=pin)
    assert not reads


def test_r4_regression_real_reversible_table_qualifies(tmp_path):
    spec, channels = reversible_request(tmp_path)
    mapping = tmp_path / 'reversible-map.json'
    mapping.write_text(json.dumps(channels))
    spec['channel_map'] = {'path': str(mapping), 'sha256': file_hash(mapping)}
    spec['axes']['u'] = np.linspace(np.log(17), np.log(17.1), 3).tolist()
    spec['held_out']['lhs_count'] = 0
    spec['refinement']['max_rounds'] = 0
    path = generate(spec)
    manifest = json.loads((path / 'manifest.json').read_text())
    reverse = [c for v in manifest['held_out_verdicts'] for c in v['checks']
               if c['rule'] == 'H3' and c['quantity'] == 'k_sup:1']
    assert reverse and all(c['passed'] for c in reverse), reverse
    assert manifest['accepted'], [v for v in manifest['held_out_verdicts'] if not v['passed']]
    table = EEDFTable.load(path, model_inputs(spec), artifact_sha256=file_hash(path / 'table.h5'))
    row = table.row(np.log(17.05), {})
    assert row.k_sup[1] > 0 and row.product_fractions[1] == 0
    assert (path / 'generation_spec.json').is_file() and (path / 'validation.md').is_file()
    from rmgpy.tools.eedf.validation import compare_row
    # A zero product population must not hide a bad reverse coefficient.
    with h5py.File(path / 'table.h5', 'r') as h5:
        direct = {name: h5['held_out'][name][0] for name in artifact.NUMERIC}
        for name in artifact.MAPPINGS:
            direct[name] = {key: dataset[0] for key, dataset in h5['held_out'][name].items()}
    predicted = table.row(float(direct['u']), {}).as_dict()
    predicted['k_sup'][1] *= 1.1
    checks = [c for c in compare_row(predicted, direct, manifest)
              if c['rule'] == 'H3' and c['quantity'] == 'k_sup:1']
    assert checks and not checks[0]['passed']


@pytest.mark.parametrize('property_name', ['statisticalWeight', 'energy'])
@pytest.mark.parametrize('alias', ['Ar(1S0,)', 'Ar(,1S0)', 'Ar(,1S0,)'])
@pytest.mark.parametrize('source', ['inline', 'legacy_file', 'json_file'])
def test_r4_regression_conflicting_selector_aliases_refuse(tmp_path, property_name, alias, source):
    spec, channels = reversible_request(tmp_path)
    entries = ['Ar(1S0) = 1', alias + ' = 100', 'Ar(3P2) = 5']
    if source == 'inline':
        spec['state_properties'][property_name] = entries
    else:
        leaf = tmp_path / ('selectors.json' if source == 'json_file' else 'selectors.txt')
        if source == 'json_file':
            leaf.write_text(json.dumps({'states': {s.split(' = ')[0]: {'type': 'constant', 'value': float(s.split(' = ')[1])} for s in entries}}))
        else:
            leaf.write_text('\n'.join(s.replace(' = ', ' ') for s in entries) + '\n')
        spec['input_files'][leaf.name] = {'path': str(leaf), 'sha256': file_hash(leaf), 'kind': 'property'}
        spec['state_properties'][property_name] = [leaf.name]
    with pytest.raises(SpecError, match='selector'):
        physical_properties(spec)
    with pytest.raises(SpecError, match='selector'):
        resolved_properties(spec, {})


@pytest.mark.parametrize('property_name', ['energy', 'statisticalWeight'])
@pytest.mark.parametrize('selector', ['Ar(*)', 'Ar()', 'Ar(+,1S0)', 'Ar(1S0,v=0)'])
def test_r4_regression_unsupported_state_tree_overrides_refuse(tmp_path, property_name, selector):
    spec, _ = request_and_channels(tmp_path)
    spec['state_properties'][property_name] = [selector + ' = 100']
    with pytest.raises(SpecError, match='selector'):
        resolved_properties(spec, {})


def test_r4_regression_artifact_is_complete_when_published(tmp_path, monkeypatch):
    original = Path.rename
    observed = []
    def rename(source, target):
        if source.name.startswith('.writing-'):
            observed.append(source)
            assert (source / 'manifest.json').is_file(), 'manifest absent at publication'
            manifest = json.loads((source / 'manifest.json').read_text())
            assert manifest['artifact_sha256'] == file_hash(source / 'table.h5')
            assert Path(target).name == manifest['artifact_sha256']
        return original(source, target)
    monkeypatch.setattr(Path, 'rename', rename)
    fixture_table(tmp_path)
    assert len(observed) == 1


def test_r4_regression_incomplete_artifact_has_named_refusal(tmp_path):
    path, model = fixture_table(tmp_path)
    (path / 'manifest.json').unlink()
    with pytest.raises(FingerprintMismatch, match='manifest'):
        EEDFTable.load(path, model, artifact_sha256=file_hash(path / 'table.h5'))


def test_r4_regression_manifest_failure_does_not_publish(tmp_path, monkeypatch):
    original = Path.write_text
    def write(path, text, *args, **kwargs):
        if path.name == 'manifest.json':
            raise OSError('interrupted manifest write')
        return original(path, text, *args, **kwargs)
    monkeypatch.setattr(Path, 'write_text', write)
    with pytest.raises(OSError, match='interrupted'):
        fixture_table(tmp_path)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('alias', ['Ar(1S0,)', 'Ar(,1S0)', 'Ar(,1S0,)'])
def test_resolved_single_alias_matches_real_loki(tmp_path, alias):
    spec, channels = reversible_request(tmp_path)
    spec['state_properties']['population'] = [alias + ' = 1']
    spec['state_properties']['statisticalWeight'] = [alias + ' = 1', 'Ar(3P2) = 5']
    spec['state_properties']['energy'] = [alias + ' = 0']
    row = enrich_row(LoKIDriver(spec).run({}, [17.], 'single-alias')[0], channels, spec)
    assert row['state_populations'] == {'Ar(1S0)': 1.}
    from rmgpy.tools.eedf.moments import distribution_moments
    moments = distribution_moments(dict(row, gas_temperature_K=spec['Tg_K']), channels)
    assert moments['k_sup'][1] == pytest.approx(row['k_sup'][1], rel=1e-6)


def test_writer_refuses_external_storage_before_publication(tmp_path, monkeypatch):
    original = artifact.write_rows
    def write(group, rows):
        original(group, rows)
        data = group['k_ine'][...]
        del group['k_ine']
        group.create_dataset('k_ine', data=data, external=[(str(tmp_path / 'rates.raw'), 0, h5py.h5f.UNLIMITED)])
    monkeypatch.setattr(artifact, 'write_rows', write)
    with pytest.raises(FingerprintMismatch, match='HDF5 storage'):
        fixture_table(tmp_path / 'tables')
    assert list((tmp_path / 'tables').iterdir()) == []


@pytest.fixture(scope='module')
def real_forward_table(tmp_path_factory):
    root = tmp_path_factory.mktemp('r4-qualified-forward')
    spec, _ = request_and_channels(root)
    spec['axes']['u'] = np.linspace(np.log(17), np.log(17.1), 3).tolist()
    spec['held_out']['lhs_count'] = 0
    spec['refinement']['max_rounds'] = 0
    path = generate(spec)
    pin = file_hash(path / 'table.h5')
    EEDFTable.load(path, model_inputs(spec), artifact_sha256=pin)
    return path, model_inputs(spec)


@pytest.mark.parametrize('storage', ['external_raw', 'virtual'])
def test_r4_regression_real_qualified_storage_refuses(tmp_path, real_forward_table, storage):
    import shutil
    source, model = real_forward_table
    path = tmp_path / 'table'
    shutil.copytree(source, path)
    raw = tmp_path / 'rates.raw'
    with h5py.File(path / 'table.h5', 'r+') as h5:
        group = h5['branches/branch_0']
        rates = group['k_ine'][...]
        del group['k_ine']
        if storage == 'external_raw':
            group.create_dataset('k_ine', data=rates, external=[(str(raw), 0, h5py.h5f.UNLIMITED)])
        else:
            other = tmp_path / 'source.h5'
            with h5py.File(other, 'w') as h5source:
                h5source['rates'] = rates
            layout = h5py.VirtualLayout(shape=rates.shape, dtype=rates.dtype)
            layout[:] = h5py.VirtualSource(str(other), 'rates', shape=rates.shape)
            group.create_virtual_dataset('k_ine', layout)
    pin = repin(path)
    if storage == 'external_raw':
        raw.write_bytes((rates * 1e6).tobytes())
    assert file_hash(path / 'table.h5') == pin
    with pytest.raises(FingerprintMismatch, match='HDF5 storage'):
        EEDFTable.load(path, model, artifact_sha256=pin)
