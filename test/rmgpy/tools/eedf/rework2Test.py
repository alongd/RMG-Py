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

from pathlib import Path
import copy
import json
import os
from pathlib import Path

import numpy as np
import pytest

from rmgpy.tools.eedf.generate import model_inputs
from rmgpy.tools.eedf.integrity import physical_properties, solver_environment
from rmgpy.tools.eedf.loki import LoKIDriver, LoKIError, enrich_row, parse_output, setup_text
from rmgpy.tools.eedf.schema import SpecError, file_hash, load_spec
from rmgpy.tools.eedf.validation import channel_power_checks

FIXTURE = Path(__file__).parent / 'fixtures' / 'ar'


def spec():
    return load_spec(os.environ['EEDF_REAL_SPEC'])


def test_r2_01_ld_audit_and_all_parent_environment_are_dropped(monkeypatch):
    request = spec()
    monkeypatch.setattr(os, 'environ', {})
    for name in ('LD_AUDIT', 'LD_DEBUG', 'LD_BIND_NOW', 'LD_PROFILE',
                 'GLIBC_TUNABLES', 'MALLOC_PERTURB_', 'LC_NUMERIC', 'LANG',
                 'OMP_DYNAMIC', 'GOMP_CPU_AFFINITY', 'LOKI_UNLISTED_SWITCH'):
        monkeypatch.setenv(name, 'untrusted')
    env = solver_environment(request)
    assert env == {'PATH': '/usr/bin:/bin', 'LANG': 'C', 'LC_ALL': 'C', 'TZ': 'UTC',
                   'OMP_NUM_THREADS': '1', 'OMP_DYNAMIC': 'FALSE'}


def test_r2_02_population_filename_ending_equals_one_must_be_hashed(tmp_path):
    request = spec()
    filename = tmp_path / 'population=1'
    filename.write_text('Ar(1S0) 1.0\n')
    request['state_properties']['population'] = [str(filename)]
    with pytest.raises(SpecError, match='unhashed'):
        model_inputs(request)


def test_r2_03_wrong_map_target_population_is_refused(tmp_path):
    request = spec()
    channels = json.loads(Path(request['channel_map']['path']).read_text())
    # This channel is too small for the former fixed absolute power floor.
    channels[5]['target_fraction'] = .1
    request['scratch_root'] = str(tmp_path / 'loki')
    row = LoKIDriver(request).run({}, [3.], 'wrong-population')[0]
    with pytest.raises((SpecError, LoKIError), match='population'):
        enrich_row(row, channels, request)


def test_r2_03_group_budget_cannot_hide_ninety_percent_material_channel_error():
    request = spec()
    # The absolute floor previously swallowed the entire difference.
    row = {'channel_power': np.array([1e-25]),
           'power_groups': {'field': 1e-23, 'excitation_loss': -1e-24}}
    manifest = {'channel_map': [{'kind': 'excitation'}],
                'tolerances': request['tolerances'], 'floors': request['floors']}
    checks = channel_power_checks(row, manifest)
    assert not next(v for v in checks if v['quantity'] == 'power.sum.excitation')['passed']


@pytest.mark.parametrize('property_name', ['population', 'energy', 'statisticalWeight'])
@pytest.mark.parametrize('syntax', ['equal', 'quoted', 'comment', 'json_nested', 'legacy_nested'])
def test_state_file_grammar_closes_all_property_keys(tmp_path, property_name, syntax):
    request = spec()
    name = 'leaf=1'
    value = tmp_path / name
    value.write_text('Ar(1S0) 1.0\n')
    reference = name
    if syntax == 'quoted':
        reference = json.dumps(name)
    elif syntax == 'comment':
        reference += ' % annotation'
    elif syntax in ('json_nested', 'legacy_nested'):
        parent_name = 'parent.json' if syntax == 'json_nested' else 'parent.in'
        parent = tmp_path / parent_name
        parent.write_text(json.dumps({'files': [name]}) if syntax == 'json_nested' else name+'\n')
        request['input_files'][parent_name] = {'path': str(parent), 'sha256': file_hash(parent), 'kind': 'property'}
        reference = parent_name
    request['state_properties'][property_name] = [reference]
    with pytest.raises(SpecError, match='unhashed'):
        physical_properties(request)
    request['input_files'][name] = {'path': str(value), 'sha256': file_hash(value), 'kind': 'property'}
    _, states = physical_properties(request)
    assert states[property_name] == ['Ar(1S0) = 1.0']


@pytest.mark.parametrize('key', ['LXCatFiles', 'LXCatFilesExtra', 'effectiveCrossSectionPopulations'])
def test_every_cross_section_file_key_requires_a_hash(key):
    request = spec()
    request['solver_options'][key] = ['unlisted=1']
    with pytest.raises(SpecError, match='unhashed'):
        physical_properties(request)


@pytest.mark.parametrize('syntax', ['duplicate', 'wildcard', 'hierarchy', 'charge', 'prefix', 'symbolic'])
def test_ambiguous_population_selectors_are_refused(syntax):
    request = spec()
    entries = {'duplicate': ['Ar(1S0) = 1', 'Ar(1S0) = 0.1'],
               'wildcard': ['Ar(*) = 1'], 'hierarchy': ['Ar(1S0,v=0) = 1'],
               'charge': ['Ar(+) = 1'], 'prefix': ['ignored Ar(1S0) = 1'], 'symbolic': ['Ar(1S0) = constantValue@1']}
    request['state_properties']['population'] = entries[syntax]
    with pytest.raises(SpecError):
        setup_text(request, {}, [17.], 'solve')


@pytest.mark.parametrize('direction', ['target_fraction', 'product_fraction'])
def test_both_map_populations_compare_bits(direction):
    request = spec()
    channels = json.loads(Path(request['channel_map']['path']).read_text())
    original = channels[1][direction]
    channels[1][direction] = float(np.nextafter(original, np.inf))
    row = parse_output(FIXTURE, request['tolerances']['normalization_atol'])
    row['composition'] = {}
    with pytest.raises(SpecError, match='population'):
        enrich_row(row, channels, request)


def test_real_solve_resolves_hashed_population_file_and_ignores_parent_environment(tmp_path, monkeypatch):
    request = spec()
    parent = tmp_path / 'population.json'
    leaf = tmp_path / 'population=1'
    leaf.write_text('Ar(1S0) 1.0 % numeric leaf\n')
    parent.write_text(json.dumps({'files': ['population=1']}))
    for file in (parent, leaf):
        request['input_files'][file.name] = {'path': str(file), 'sha256': file_hash(file), 'kind': 'property'}
    request['state_properties']['population'] = [parent.name]
    request['scratch_root'] = str(tmp_path / 'loki')
    baseline = model_inputs(request)
    monkeypatch.setenv('LD_AUDIT', '/untrusted/audit.so')
    monkeypatch.setenv('GLIBC_TUNABLES', 'untrusted')
    monkeypatch.setenv('MALLOC_PERTURB_', '255')
    channels = json.loads(Path(request['channel_map']['path']).read_text())
    row = LoKIDriver(request).run({}, [17.], 'hashed-population')[0]
    enrich_row(row, channels, request)
    assert row['target_fractions'][1] == 1.
    assert 'audit' not in (tmp_path / 'loki/hashed-population/stderr.log').read_text().lower()
    leaf.write_text('Ar(1S0) 1.00\n')
    request['input_files'][leaf.name]['sha256'] = file_hash(leaf)
    assert model_inputs(request) != baseline


def test_gas_fraction_exponent_is_normalized_for_pinned_parser():
    request = spec()
    request['gas_properties']['fraction'] = ['Ar = 1e0']
    assert 'Ar = 1' in setup_text(request, {}, [17.], 'solve')


def test_loader_rejects_map_population_disagreement_even_with_current_artifact_pin(tmp_path):
    import h5py
    from rmgpy.tools.eedf.schema import content_hash
    from rmgpy.tools.eedf.generate import generate
    from rmgpy.solver.eedf import EEDFTable
    request = spec()
    request['axes']['u'] = np.linspace(np.log(16), np.log(18), 9).tolist()
    request['held_out']['lhs_count'] = 0
    request['screen']['candidates'] = []
    for bounds in request['envelopes'].values():
        bounds['min'] = bounds['reference']
        bounds['max'] = bounds['reference']
    request['refinement']['max_rounds'] = 0
    request['scratch_root'] = str(tmp_path / 'loki')
    request['output_root'] = str(tmp_path / 'tables')
    path = generate(request)
    manifest_path = path / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    manifest['channel_map'][1]['target_fraction'] = .1
    with h5py.File(path / 'table.h5', 'r+') as h5:
        h5.attrs['channel_map_content_sha256'] = content_hash(manifest['channel_map'])
    digest = file_hash(path / 'table.h5')
    manifest['artifact_sha256'] = digest
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='populations'):
        EEDFTable.load(path, model_inputs(request), artifact_sha256=digest, require_accepted=False)


def test_group_budget_catches_channel_undercounted_below_materiality_threshold():
    request = spec()
    # The smaller channel should contribute 1e-20 (0.1% of field), but
    # contributes only 1e-21. A 0.5%-of-total budget would hide its 90% error.
    row = {'channel_power': np.array([9.99e-18, 1e-21]),
           'power_groups': {'field': 1e-17, 'excitation_loss': -1e-17}}
    manifest = {'channel_map': [{'kind': 'excitation'}, {'kind': 'excitation'}],
                'tolerances': request['tolerances'], 'floors': request['floors']}
    checks = channel_power_checks(row, manifest)
    assert not next(v for v in checks if v['quantity'] == 'power.sum.excitation')['passed']


@pytest.mark.parametrize('name', ['mass', 'harmonicFrequency', 'anharmonicFrequency',
                                  'rotationalConstant', 'electricQuadrupoleMoment', 'OPBParameter'])
def test_every_gas_property_file_key_requires_a_hash(name):
    request = spec()
    request['gas_properties'][name] = 'unlisted=1'
    with pytest.raises(SpecError, match='unhashed'):
        physical_properties(request)


@pytest.mark.parametrize('entry', ['Ar(1S0)=1', 'Ar(1S0)= 1', 'Ar(1S0) =1'])
def test_equals_without_surrounding_whitespace_is_a_file_reference(entry):
    request = spec()
    request['state_properties']['population'] = [entry]
    with pytest.raises(SpecError, match='unhashed'):
        physical_properties(request)


@pytest.mark.parametrize('entry', ['Ar(1S0) 1.2.3', 'Ar(1S0) = 1', '# ignored', '// ignored'])
def test_ambiguous_legacy_property_file_lines_are_refused(tmp_path, entry):
    request = spec()
    file = tmp_path / 'population.in'
    file.write_text(entry+'\n')
    request['input_files'][file.name] = {'path': str(file), 'sha256': file_hash(file), 'kind': 'property'}
    request['state_properties']['population'] = [file.name]
    with pytest.raises(SpecError):
        physical_properties(request)


def test_setup_list_cannot_inject_another_file_key():
    request = spec()
    request['working_conditions']['extra'] = ['1\n  LXCatFilesExtra:\n    - /unlisted']
    with pytest.raises(LoKIError, match='ambiguous'):
        setup_text(request, {}, [17.], 'solve')
