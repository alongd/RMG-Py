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
import inspect
import json
import os
import shutil
from pathlib import Path

import h5py
import numpy as np
import pytest

from rmgpy.solver.eedf import EEDFTable
from rmgpy.tools.eedf import loki
from rmgpy.tools.eedf.generate import model_inputs
from rmgpy.tools.eedf.schema import SpecError, FingerprintMismatch, file_hash, content_hash
from rmgpy.tools.eedf.validation import compare_row
import runpy
fixture_table = runpy.run_path(str(Path(__file__).parents[2] / 'solver' / 'eedfTableTest.py'))['fixture_table']

FIXTURE = Path(__file__).parent / 'fixtures' / 'ar'


def raw_spec():
    return json.loads(Path(os.environ['EEDF_REAL_SPEC']).read_text())


def pinned(path, model, **kw):
    # Execute the same behavioral regression against both API generations.
    if 'artifact_sha256' in inspect.signature(EEDFTable.load).parameters:
        kw['artifact_sha256'] = kw.get('artifact_sha256', file_hash(path / 'table.h5'))
    return EEDFTable.load(path, model, **kw)


def mutate(path, edit):
    manifest_path = path / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    with h5py.File(path / 'table.h5', 'r+') as h5:
        edit(h5, manifest)
    manifest['artifact_sha256'] = file_hash(path / 'table.h5')
    manifest_path.write_text(json.dumps(manifest))


def copied(tmp_path):
    dest=tmp_path/'output'
    shutil.copytree(FIXTURE,dest)
    return dest


def test_mf01_shared_solver_hash_is_checked():
    spec=raw_spec()
    # Supplying a false pin must be noticed even though the executable is intact.
    build=Path(spec['cmake_cache']['path']).parent
    spec['shared_objects']={str(build/'lib/libloki-b.so'): '0'*64}
    with pytest.raises(FingerprintMismatch, match='shared'):
        loki.LoKIDriver(spec)


def test_mf02_absolute_unlisted_property_file_is_refused():
    spec=raw_spec()
    spec['gas_properties']['mass']=spec['input_files']['Databases/masses.txt']['path']
    del spec['input_files']['Databases/masses.txt']
    with pytest.raises(SpecError, match='unhashed'):
        model_inputs(spec)


def test_mf03_rehashed_hdf5_tamper_does_not_change_callers_pin(tmp_path):
    path,model=fixture_table(tmp_path)
    original=file_hash(path/'table.h5')
    mutate(path, lambda h,m: h['branches/branch_0/k_ine'].__setitem__(Ellipsis, h['branches/branch_0/k_ine'][...]*1e6))
    kw={'artifact_sha256': original} if 'artifact_sha256' in inspect.signature(EEDFTable.load).parameters else {}
    with pytest.raises(FingerprintMismatch, match='artifact'):
        pinned(path,model,**kw)


def test_mf03_stored_acceptance_cannot_override_failing_checks(tmp_path):
    path,model=fixture_table(tmp_path)
    def edit(h,m):
        m['accepted']=True
        m['held_out_verdicts']=[{'passed':True,'checks':[{'passed':False}]}]
        h.attrs['qualification_sha256']=content_hash({k:m.get(k) for k in ('accepted','held_out_verdicts','screen_results','branch_detection')})
    mutate(path,edit)
    with pytest.raises(FingerprintMismatch, match='qualification'):
        pinned(path,model)


def test_mf03_loader_requires_external_artifact_pin(tmp_path):
    path,model=fixture_table(tmp_path)
    with pytest.raises((TypeError,FingerprintMismatch)):
        EEDFTable.load(path,model)


@pytest.mark.parametrize('name', ['rate_absolute','eedf_dynamic_range'])
def test_mf04_rate_policy_changes_model_identity(name):
    spec=raw_spec()
    first=model_inputs(spec)
    spec['floors'][name] *= 10
    assert model_inputs(spec)!=first


def test_mf05_perturbed_setup_preserves_temperature_and_pressure():
    spec=raw_spec()
    text=loki.setup_text(spec,{'Tg_K':350.,'P_Pa':777.},[17.],'solve')
    assert 'gasTemperature: 350.0' in text
    assert 'gasPressure: 777.0' in text


def test_mf06_reversed_eedf_fails_held_out_qualification():
    spec=raw_spec()
    row=loki.parse_output(FIXTURE,1e-8)
    row['composition']={}
    channels=json.loads(Path(spec['channel_map']['path']).read_text())
    loki.enrich_row(row,channels,spec)
    row['product_fractions']=np.zeros(len(channels))
    predicted=copy.deepcopy(row)
    predicted['f0']=predicted['f0'][::-1]
    verdicts=compare_row(predicted,row,{'tolerances':spec['tolerances'],'floors':spec['floors'],'channel_map':channels,'row_inputs': {'Tg_K':spec['Tg_K']}})
    assert any(not v['passed'] and 'f0' in v['quantity'] for v in verdicts)


def test_mf07_row_grid_mismatch_refuses_storage(tmp_path):
    from rmgpy.tools.eedf.artifact import write_artifact
    path,model=fixture_table(tmp_path/'initial')
    table=pinned(path,model)
    manifest=copy.deepcopy(table.manifest)
    rows=np.empty((4,),dtype=object)
    for i in range(4):
        rows[i]=table.row(float(i),{}).as_dict()
    rows[2]['energy_eV']=np.array([.6])
    with pytest.raises((SpecError,FingerprintMismatch),match='grid'):
        write_artifact(tmp_path/'bad',manifest,{'branch_0':rows},[])


def test_mf08_constant_symbolic_population_refuses(tmp_path):
    path,model=fixture_table(tmp_path)
    def edit(h,m):
        h['branches/branch_0/state_populations'].attrs['unresolved']=json.dumps(['boltzmannPopulation(298.15)']*4)
    mutate(path,edit)
    with pytest.raises(FingerprintMismatch,match='symbolic'):
        pinned(path,model)


@pytest.mark.parametrize('filename,old,new',[
    ('rateCoefficients.txt','Ine.R.Coeff.(m3/s)   Sup.R.Coeff.(m3/s)','Sup.R.Coeff.(m3/s)   Ine.R.Coeff.(m3/s)'),
    ('eedf.txt','Energy(eV)','Unknown(eV)'),
])
def test_mf09_swapped_or_unknown_headers_refuse(tmp_path,filename,old,new):
    folder=copied(tmp_path)
    p=folder/filename
    p.write_text(p.read_text().replace(old,new,1))
    with pytest.raises(loki.LoKIError,match='header'):
        loki.parse_output(folder,1e-8)


def test_mf09_fortran_exponents_keep_exponent(tmp_path):
    folder=copied(tmp_path)
    for p in folder.iterdir():
        p.write_text(p.read_text().replace('e+','D+').replace('e-','D-'))
    row=loki.parse_output(folder,1e-8)
    assert row['swarm']['mobility_N']==pytest.approx(1.01724123397156e24)
    assert row['k_ine'][-1]==pytest.approx(7.44149321789454e-21,abs=0)


def test_mf09_negative_energy_refuses_even_with_loose_normalization(tmp_path):
    folder=copied(tmp_path)
    p=folder/'eedf.txt'
    p.write_text(p.read_text().replace('5.00000000000000e-03','-5.00000000000000e-03',1))
    with pytest.raises(loki.LoKIError,match='grid'):
        loki.parse_output(folder,1.)


def test_mf10_elastic_channel_power_is_retained():
    spec=raw_spec()
    row=loki.parse_output(FIXTURE,1e-8)
    row['composition']={}
    channels=json.loads(Path(spec['channel_map']['path']).read_text())
    loki.enrich_row(row,channels,spec)
    assert row['channel_power'][0]==pytest.approx(-row['power_groups']['elastic_loss']-row['power_groups']['elastic_gain'],rel=1e-12,abs=0)


def test_mf10_wrong_threshold_fails_group_power_check():
    spec=raw_spec()
    row=loki.parse_output(FIXTURE,1e-8)
    row['composition']={}
    channels=json.loads(Path(spec['channel_map']['path']).read_text())
    channels[1]['threshold_eV']*=100
    with pytest.raises((loki.LoKIError, SpecError),match='power|threshold'):
        loki.enrich_row(row,channels,spec)


def test_mf11_composition_fixture_rates_depend_on_composition(tmp_path):
    path,model=fixture_table(tmp_path,composition=True)
    table=pinned(path,model)
    assert table.row(1.,{'x':0.}).k_ine[0]!=table.row(1.,{'x':1.}).k_ine[0]


def test_mf11_numeric_policy_test_reaches_numeric_validation():
    from rmgpy.tools.eedf.schema import validate_spec
    spec=raw_spec()
    spec['tolerances']['H1']['rtol']='not a number'
    with pytest.raises(SpecError,match='numeric'):
        validate_spec(spec)


def test_shared_loader_environment_clears_overrides(monkeypatch):
    from rmgpy.tools.eedf.integrity import solver_environment
    monkeypatch.setenv('LD_LIBRARY_PATH', '/untrusted')
    monkeypatch.setenv('LD_PRELOAD', '/untrusted.so')
    env = solver_environment(raw_spec())
    assert 'LD_LIBRARY_PATH' not in env and 'LD_PRELOAD' not in env


def test_relative_nested_population_file_must_be_hashed(tmp_path):
    spec = raw_spec()
    population = tmp_path / 'population.txt'
    population.write_text('nested.txt\n')
    spec['input_files']['population.txt'] = {'path': str(population), 'sha256': file_hash(population), 'kind': 'property'}
    spec['state_properties']['population'] = ['population.txt']
    with pytest.raises(SpecError, match='unhashed'):
        model_inputs(spec)


def test_numeric_population_file_is_resolved_and_stored(tmp_path):
    spec = raw_spec()
    population = tmp_path / 'population.txt'
    population.write_text('Ar(1S0) 1.0\n')
    spec['input_files']['population.txt'] = {'path': str(population), 'sha256': file_hash(population), 'kind': 'property'}
    spec['state_properties']['population'] = ['population.txt']
    row = loki.parse_output(FIXTURE, 1e-8)
    row['composition'] = {}
    channels = json.loads(Path(spec['channel_map']['path']).read_text())
    loki.enrich_row(row, channels, spec)
    assert row['state_populations'] == {'Ar(1S0)': 1.}
    assert 'Ar(1S0) = 1' in loki.setup_text(spec, {}, [17.], 'solve')


def test_constant_symbolic_population_is_refused_before_execution():
    spec = raw_spec()
    spec['state_properties']['population'] = ['Ar(1S0) = boltzmannPopulation(298.15)']
    with pytest.raises(SpecError, match='symbolic'):
        model_inputs(spec)


@pytest.mark.parametrize('filename,old,new', [
    ('rateCoefficients.txt', '1.22259597013409e-13', 'NaN'),
    ('swarmParameters.txt', '1.01724123397156e+24', 'Inf'),
    ('eedf.txt', '6.98108305338813e-02', 'NaN'),
])
def test_nonfinite_solver_values_refuse(tmp_path, filename, old, new):
    folder = copied(tmp_path)
    path = folder / filename
    text = path.read_text()
    assert old in text
    path.write_text(text.replace(old, new, 1))
    with pytest.raises(loki.LoKIError):
        loki.parse_output(folder, 1e-8)


def test_loader_hashes_and_reads_one_snapshot(tmp_path, monkeypatch):
    path, model = fixture_table(tmp_path)
    expected = file_hash(path / 'table.h5')
    real_open = Path.open
    count = []
    def opened(self, *args, **kw):
        stream = real_open(self, *args, **kw)
        if self == path / 'table.h5':
            count.append(self)
            # Replace the path after its descriptor was opened. The consumed
            # inode must still be the bytes covered by the caller's pin.
            replacement = tmp_path / 'replacement'
            replacement.write_bytes(b'untrusted HDF5')
            replacement.replace(self)
        return stream
    monkeypatch.setattr(Path, 'open', opened)
    table = EEDFTable.load(path, model, artifact_sha256=expected)
    assert len(count) == 1
    assert table.row(1., {}).k_ine[0] == pytest.approx(2e-15, rel=1e-12, abs=0)


def test_normalization_overflow_refuses(tmp_path):
    path, model = fixture_table(tmp_path)
    def edit(h, m):
        h['branches/branch_0/f0'][...] = 1e308
        m['energy_eV'] = [1e308]
        m['energy_edges_eV'] = [0., 1e308]
        h['energy_eV'][...] = m['energy_eV']
        h['energy_edges_eV'][...] = m['energy_edges_eV']
        h['branches/branch_0/energy_eV'][...] = 1e308
        h['branches/branch_0/energy_edges_eV'][...] = [0., 1e308]
    mutate(path, edit)
    with pytest.raises(Exception, match='normalization'):
        pinned(path, model).row(1., {})


def test_hf_refuses_at_spec_validation():
    from rmgpy.tools.eedf.schema import validate_spec
    spec = raw_spec()
    spec['working_conditions']['excitationFrequency'] = 1.
    with pytest.raises(SpecError, match='excitationFrequency'):
        validate_spec(spec)


def test_refinement_runs_and_keeps_original_midpoint_held_out(tmp_path, monkeypatch):
    import importlib
    generation = importlib.import_module('rmgpy.tools.eedf.generate')
    from rmgpy.tools.eedf.schema import load_spec
    spec = load_spec(os.environ['EEDF_REAL_SPEC'])
    spec['axes'] = {'u': [0., 1.]}
    spec['refinement']['max_rounds'] = 1
    spec['held_out']['lhs_count'] = 0
    spec['scratch_root'] = str(tmp_path / 'scratch')
    spec['output_root'] = str(tmp_path / 'tables')
    direct = loki.parse_output(FIXTURE, 1e-8)
    direct['composition'] = {}
    channels = json.loads(Path(spec['channel_map']['path']).read_text())
    loki.enrich_row(direct, channels, spec)
    direct['product_fractions'] = np.zeros(len(channels))
    direct['power_absolute'] = 0.
    def solve(driver, channel_map, local_spec, coordinates, fields, label):
        rows = []
        for field in fields:
            row = copy.deepcopy(direct)
            row['EN_Td'], row['u'], row['composition'], row['setup'] = field, float(np.log(field)), coordinates, label
            row['swarm']['mean_energy_eV'] = 1 + row['u']
            rows.append(row)
        return rows
    monkeypatch.setattr(generation, '_solve', solve)
    monkeypatch.setattr(generation, '_screen', lambda *args: [])
    def compare(predicted, reference, manifest):
        passed = len(manifest['axes']['u']) > 2
        return [{'passed': passed, 'rule': 'H1', 'relative_error': 0., 'quantity': 'fixture', 'absolute_error': 0., 'allowed_error': 1.}]
    monkeypatch.setattr(generation, 'compare_row', compare)
    path = generation.generate(spec)
    manifest = json.loads((path / 'manifest.json').read_text())
    assert len(manifest['refinement_history']) == 2
    assert len(manifest['axes']['u']) == 4
    assert .5 not in manifest['axes']['u']
    with h5py.File(path / 'table.h5', 'r') as h5:
        assert .5 in h5['held_out/u'][...]
        assert .5 not in h5['branches/branch_0/u'][...]
    assert manifest['accepted']
