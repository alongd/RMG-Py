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

"""Checks that benchmark closeout cannot promote partial results to validation."""
import copy
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tarfile

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('benchmark_report', ROOT / 'scripts/plasmaBenchmarkReport.py')
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)


def test_direction_failure_does_not_borrow_downward_value():
    rows = [{'direction': 'ascending', 'pressure_torr': 5, 'power_W': .01, 'status': 'failed',
             'error': 'PlasmaStateError("negative population")'},
            {'direction': 'descending', 'pressure_torr': 5, 'power_W': .01, 'status': 'steady',
             'ne_m3': 7., 'Te_eV': 1.}]
    pair = report.compare_directions(rows)[0]
    assert pair['ascending_status'] == 'solver_refused'
    assert pair['ascending_ne_m3'] is None
    assert pair['descending_ne_m3'] == 7.
    assert pair['ne_m3_relative_difference'] is None
    with pytest.raises(ValueError, match='Missing continuation'):
        report.compare_directions(rows[:1])
    with pytest.raises(ValueError, match='duplicate'):
        report.compare_directions(rows + rows)


def test_nonfinite_steady_endpoint_is_rejected():
    row = {'direction': 'ascending', 'pressure_torr': 5, 'power_W': .01, 'status': 'steady',
           'ne_m3': float('nan'), 'Te_eV': 1.}
    with pytest.raises(ValueError, match='Nonfinite'):
        report.compare_directions([row, dict(row, direction='descending')])


def power_fixture():
    base = {'native_field_W_m3': 10., 'P_abs_W_m3': 9., 'loss_W_m3': 9.,
            'reaction_power_W_m3': [3., 1., -.5], 'consumed_electron_power_W_m3': [0., .2, 0.],
            'elastic_W_m3': 4., 'wall_electron_W_m3': .5, 'wall_ion_W_m3': .8}
    manifest = {'reactions': [{'id': name} for name in ('exc', 'ion', 'sup')]}
    groups = {'excitation': ['exc'], 'ionization': ['ion']}
    text = '\n'.join(['Field = +1.0e+01 (eVm^3s^-1)',
                      'Elastic collisions (net) = -6.0e+00 (eVm^3s^-1)',
                      'Excitation collisions (net) = -3.0e+00 (eVm^3s^-1)',
                      'Ionization collisions = -1.0e+00 (eVm^3s^-1)'])
    return base, text, manifest, groups


def test_power_identity_alignment_and_signed_terms():
    base, text, manifest, groups = power_fixture()
    original = report.power_accounting(base, text, manifest, groups)
    assert original['ablation_confirmed'] is False
    assert original['rows'][3]['matched_W_m3'] == -.5
    assert original['gap_to_matched_loss_W_m3'] == pytest.approx(1.)
    permuted = copy.deepcopy(base)
    permuted['reaction_power_W_m3'] = list(reversed(base['reaction_power_W_m3']))
    permuted['consumed_electron_power_W_m3'] = list(reversed(base['consumed_electron_power_W_m3']))
    manifest['reactions'].reverse()
    assert report.power_accounting(permuted, text, manifest, groups) == original


def test_incomplete_power_inventory_is_rejected():
    base, text, manifest, groups = power_fixture()
    with pytest.raises(ValueError, match='do not close'):
        report.power_accounting(base, text.replace('-6.0e+00', '-5.0e+00'), manifest, groups)
    groups['ionization'] = ['unknown']
    with pytest.raises(ValueError, match='identity'):
        report.power_accounting(base, text, manifest, groups)
    with pytest.raises(ValueError, match='Missing or invalid'):
        report.power_accounting(base, 'empty', manifest, groups)


def native_fixture(tmp_path):
    stage, fetch = tmp_path/'stage', tmp_path/'fetch'
    stage.mkdir();(fetch/'N0').mkdir(parents=True)
    ids = ['N0'] + ['arm' + str(i) for i in range(16)]
    (stage/'native-modes.json').write_text(json.dumps([{'id': name} for name in ids]))
    (stage/'source.txt').write_text('synthetic test source')
    for kind, name in [('upstream', 'source.txt'), ('adapter', 'native-modes.json')]:
        (fetch/(kind+'.sha256')).write_text(report.digest(stage/name)+'  '+name+'\n')
    for name, value in [('job', 1), ('matlab-stage', 1), ('pack', 0), ('upstream-check', 0), ('adapter-check', 0)]:
        (fetch/(name+'.exit')).write_text(str(value))
    (fetch/'N0/trial-01.exit').write_text('1')
    config = {'workingConditions': {'electronTemperature': 1., 'electronDensity': 10.},
              'chemistry': {'iterationSchemes': {'neutralityMaxIterations': 3, 'neutralityRelError': 1e-7}}}
    (fetch/'N0/trial-01-config.json').write_text(json.dumps(config))
    stdout = '\n'.join('New neutrality cycle iteration ({}): relative error = {}\n'
                       'Updated electron temperature (0.9 eV).'.format(i, error)
                       for i, error in enumerate(['1.0e-1', '-1.0e-2', '2.0e-3'], 1))
    (fetch/'N0/trial-01.stdout').write_text(stdout)
    (fetch/'N0/trial-01.stderr').write_text('Error using Chemistry/neutralityCycle\n'
        'Maximum number of iterations reached for the neutrality cycle without\nconvergence.')
    (fetch/'native-index.json').write_text(json.dumps([
        {'arm_id': name, 'status': 'failed' if i == 0 else 'not_run'} for i, name in enumerate(ids)]))
    payload = b'partial output'
    name = 'Output/log.txt'
    (fetch/'output.sha256').write_text(hashlib.sha256(payload).hexdigest()+'  '+name+'\n')
    with tarfile.open(fetch/'loki-output.tar.gz', 'w:gz') as archive:
        member = tarfile.TarInfo(name);member.size = len(payload)
        archive.addfile(member, io.BytesIO(payload))
    return fetch, stage


def test_native_failure_is_not_an_accepted_endpoint(tmp_path):
    fetch, stage = native_fixture(tmp_path)
    result = report.audit_native(fetch, stage)
    assert result['status'] == 'native-non-convergent'
    assert result['iteration_cap'] == 3
    assert result['other_arms_not_run'] == 16
    assert result['ablation_confirmed'] is False
    assert result['best_iteration']['iteration'] == 3
    assert result['sign_changes'] == 2


@pytest.mark.parametrize('corruption', ['source', 'archive', 'missing_iteration', 'false_success', 'unexpected_arm'])
def test_native_integrity_and_incomplete_history_fail_closed(tmp_path, corruption):
    fetch, stage = native_fixture(tmp_path)
    if corruption == 'source':
        (stage/'source.txt').write_text('modified')
    elif corruption == 'archive':
        (fetch/'output.sha256').write_text('0'*64+'  Output/log.txt\n')
    elif corruption == 'missing_iteration':
        file = fetch/'N0/trial-01.stdout'
        file.write_text(file.read_text().replace('iteration (2)', 'iteration (4)'))
    elif corruption == 'false_success':
        (fetch/'job.exit').write_text('0')
    else:
        rows = json.loads((fetch/'native-index.json').read_text())
        rows[-1]['status'] = 'steady'
        (fetch/'native-index.json').write_text(json.dumps(rows))
    with pytest.raises(ValueError):
        report.audit_native(fetch, stage)


def test_empty_hash_manifest_does_not_verify_vacuously(tmp_path):
    manifest = tmp_path/'manifest';manifest.write_text('')
    with pytest.raises(ValueError, match='Empty hash'):
        report.verify_manifest(manifest, tmp_path)
