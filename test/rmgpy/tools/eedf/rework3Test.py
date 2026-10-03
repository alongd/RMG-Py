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

import h5py
import numpy as np
import pytest

from rmgpy.solver.eedf import EEDFTable
from rmgpy.tools.eedf import artifact
from rmgpy.tools.eedf.generate import generate, model_inputs, _validate_map
from rmgpy.tools.eedf.loki import LoKIDriver, enrich_row
from rmgpy.tools.eedf.schema import EEDFError, SpecError, content_hash, file_hash, load_spec, validate_spec
import runpy
fixture_table = runpy.run_path(str(Path(__file__).parents[2] / 'solver' / 'eedfTableTest.py'))['fixture_table']


def request_and_channels(tmp_path):
    spec = load_spec(os.environ['EEDF_REAL_SPEC'])
    spec['scratch_root'] = str(tmp_path / 'loki')
    spec['output_root'] = str(tmp_path / 'tables')
    channels = json.loads(Path(spec['channel_map']['path']).read_text())
    return spec, channels


def test_r3_01_fabricated_sigma_max_cannot_disable_qualification(tmp_path):
    spec, channels = request_and_channels(tmp_path)
    row = LoKIDriver(spec).run({}, [17.], 'sigma-bound')[0]
    channels[1]['sigma_max_m2'] = 1e100
    with pytest.raises(SpecError, match='sigma_max'):
        enrich_row(row, channels, spec)


@pytest.mark.parametrize('classification', ['C', 'D'])
def test_r3_02_nonproduction_classification_refuses_load(tmp_path, classification):
    path, model = fixture_table(tmp_path)
    manifest = json.loads((path / 'manifest.json').read_text())
    manifest['channel_map'][0]['classification'] = classification
    with h5py.File(path / 'table.h5', 'r+') as h5:
        h5.attrs['channel_map_content_sha256'] = content_hash(manifest['channel_map'])
    manifest['artifact_sha256'] = file_hash(path / 'table.h5')
    (path / 'manifest.json').write_text(json.dumps(manifest))
    # Offline inspection remains possible, but these classes are never production rows.
    EEDFTable.load(path, model, artifact_sha256=manifest['artifact_sha256'], require_accepted=False)
    with pytest.raises(EEDFError, match='classification'):
        EEDFTable.load(path, model, artifact_sha256=manifest['artifact_sha256'])


def reversible_request(tmp_path):
    spec, channels = request_and_channels(tmp_path)
    entry = next(v for v in spec['input_files'].values() if v['kind'] == 'cross_section')
    original = Path(entry['path']).read_text()
    old = '[Ar(1S0) + e -> Ar(3P2) + e, Excitation]'
    assert old in original
    section = tmp_path / 'reversible-ar.txt'
    section.write_text(original.replace(old, old.replace(' -> ', ' <-> '), 1))
    entry.update(path=str(section), sha256=file_hash(section))
    spec['state_properties']['statisticalWeight'] = ['Ar(1S0) = 1', 'Ar(3P2) = 5']
    channels[1]['description'] = channels[1]['description'].replace(' -> ', ' <-> ')
    channels[1]['flux_group'] = channels[1]['description']
    channels[1]['statistical_weight_ratio'] = .2
    return spec, channels


def test_r3_03_real_reversible_loki_solve_enriches(tmp_path):
    spec, channels = reversible_request(tmp_path)
    row = LoKIDriver(spec).run({}, [17.], 'reversible')[0]
    assert row['channels'][1] == channels[1]['description']
    assert '<->' in row['channels'][1] and row['k_sup'][1] > 0
    enriched = enrich_row(row, channels, spec)
    from rmgpy.tools.eedf.moments import distribution_moments
    moments = distribution_moments(dict(enriched, gas_temperature_K=spec['Tg_K']), channels)
    assert moments['k_sup'][1] == pytest.approx(row['k_sup'][1], rel=1e-6)
    assert enriched['product_fractions'][1] == 0
    mapping = tmp_path / 'reversible-map.json'
    mapping.write_text(json.dumps(channels))
    spec['channel_map'] = {'path': str(mapping), 'sha256': file_hash(mapping)}
    spec['axes']['u'] = np.linspace(np.log(17), np.log(17.1), 3).tolist()
    spec['held_out']['lhs_count'] = 0
    spec['refinement']['max_rounds'] = 0
    path = generate(spec)
    table = EEDFTable.load(path, model_inputs(spec), artifact_sha256=file_hash(path / 'table.h5'))
    assert table.manifest['accepted']
    assert table.row(np.log(17.05), {}).k_sup[1] > 0


def test_r3_04_all_storage_identifiers_round_trip(tmp_path, monkeypatch):
    names = ["Ar(4p[1/2]1)", "Ar([state])", "Ar(state')", "Ar(state+)",
             "Ar(state*)", "Ar(state name)", 'a/b', 'a%2Fb', '.', '', 'μ']
    original = artifact.write_rows
    def write(group, rows):
        rows = copy.deepcopy(np.asarray(rows, dtype=object))
        for row in rows.flat:
            row['state_populations'].update({key: float(i) for i, key in enumerate(names)})
        original(group, rows)
    monkeypatch.setattr(artifact, 'write_rows', write)
    path, model = fixture_table(tmp_path)
    table = EEDFTable.load(path, model, artifact_sha256=file_hash(path / 'table.h5'))
    state = table.row(1., {}).state_populations
    for i, key in enumerate(names):
        assert state[key] == i
    assert len(state) == len(names) + 1


def test_r3_04_real_generation_with_slash_state_name(tmp_path):
    spec, _ = request_and_channels(tmp_path)
    spec['state_properties']['population'].append('Ar(4p[1/2]1) = 0')
    spec['axes']['u'] = np.linspace(np.log(17), np.log(17.1), 3).tolist()
    spec['held_out']['lhs_count'] = 0
    spec['refinement']['max_rounds'] = 0
    path = generate(spec)
    table = EEDFTable.load(path, model_inputs(spec), artifact_sha256=file_hash(path / 'table.h5'))
    assert table.row(np.log(17.05), {}).state_populations['Ar(4p[1/2]1)'] == 0


@pytest.mark.parametrize('value', [1.5, True])
def test_r3_05_omp_threads_requires_positive_integer(tmp_path, value):
    spec, _ = request_and_channels(tmp_path)
    spec['omp_threads'] = value
    with pytest.raises(SpecError, match='parallelism'):
        validate_spec(spec)


@pytest.mark.parametrize('field,value', [('threshold_eV', 1.), ('kind', 'ionization'),
                                        ('mass_ratio', 1.), ('opb_eV', 100.),
                                        ('cross_section', None)])
def test_r3_audit_physical_map_metadata_checked_against_solver_inputs(tmp_path, field, value):
    spec, channels = request_and_channels(tmp_path)
    index = 0 if field == 'mass_ratio' else (38 if field == 'opb_eV' else 1)
    if field == 'cross_section':
        channels[index][field]['sigma_m2'][2] *= 2
        channels[index]['sigma_max_m2'] = max(channels[index][field]['sigma_m2'])
    else:
        channels[index][field] = value
    mapping = tmp_path / 'metadata-map.json'
    mapping.write_text(json.dumps(channels))
    spec['channel_map'] = {'path': str(mapping), 'sha256': file_hash(mapping)}
    with pytest.raises(SpecError, match=field):
        model_inputs(spec)


def test_r3_audit_statistical_weight_ratio_checked_against_setup(tmp_path, monkeypatch):
    spec, channels = reversible_request(tmp_path)
    from rmgpy.tools.eedf import loki
    row = LoKIDriver(spec).run({}, [17.], 'wrong-weights')[0]
    original = loki.channel_fractions
    # Population parsing has its own regression. Feed its known forward-arrow
    # grammar here so the old SHA reaches the missing weight validation.
    def fractions(mapping, gas, populations, coordinates):
        forward = [dict(channel, description=channel['description'].replace(' <-> ', ' -> '))
                   for channel in mapping]
        return original(forward, gas, populations, coordinates)
    monkeypatch.setattr(loki, 'channel_fractions', fractions)
    channels[1]['statistical_weight_ratio'] = 100.
    with pytest.raises(SpecError, match='statistical'):
        enrich_row(row, channels, spec)


def test_r3_audit_flux_group_cannot_increase_rate_allowance(tmp_path):
    from rmgpy.tools.eedf.validation import compare_row
    spec, channels = request_and_channels(tmp_path)
    direct = enrich_row(LoKIDriver(spec).run({}, [17.], 'flux-group')[0], channels, spec)
    predicted = copy.deepcopy(direct)
    predicted['k_ine'][5] *= 2
    manifest = {'channel_map': channels, 'tolerances': spec['tolerances'],
                'floors': spec['floors'], 'row_inputs': {'Tg_K': spec['Tg_K']}}
    reference = [c for c in compare_row(predicted, direct, manifest) if c['rule'] == 'H3']
    for channel in channels:
        channel['flux_group'] = 'fabricated-total-including-elastic'
    altered = [c for c in compare_row(predicted, direct, manifest) if c['rule'] == 'H3']
    assert altered == reference


@pytest.mark.parametrize('classification,reaction', [('A', None),
    ('B', {'library': 'fabricated', 'index': 1, 'repr': 'x'})])
def test_r3_audit_classification_agrees_with_reaction_identity(tmp_path, classification, reaction):
    _, channels = request_and_channels(tmp_path)
    channels[1].update(classification=classification, reaction=reaction)
    with pytest.raises(SpecError, match='reaction'):
        _validate_map(channels)


def test_r3_audit_axis_branch_and_json_identifiers_round_trip(tmp_path, monkeypatch):
    namespace = fixture_table.__globals__
    original = namespace['write_artifact']
    axis, branch = "x / [state]' + *", "branch / [0]' + *"
    def write(root, manifest, branches, held):
        manifest = copy.deepcopy(manifest)
        manifest['axes'][axis] = manifest['axes'].pop('x')
        manifest['row_inputs']['axes'] = copy.deepcopy(manifest['axes'])
        manifest['fingerprint'] = content_hash(manifest['row_inputs'])
        manifest['branch_certification'] = {branch: 'uncertified: unseeded scans'}
        for verdict in manifest['held_out_verdicts']:
            verdict['branch_id'] = branch
        return original(root, manifest, {branch: next(iter(branches.values()))}, held)
    monkeypatch.setitem(namespace, 'write_artifact', write)
    path, model = fixture_table(tmp_path, composition=True)
    model['axes'][axis] = model['axes'].pop('x')
    table = EEDFTable.load(path, model, artifact_sha256=file_hash(path / 'table.h5'))
    row = table.row(1., {axis: .4}, branch)
    assert row.branch_id == branch and row.composition == {axis: .4}
    assert table.axis_names == ['u', axis]
