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

"""Failure evidence and source-identity guards for the field-drive job."""
import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

RECIPE = Path(__file__).resolve().parents[4] / 'scripts/eedf_field_drive'


def load_job_module(name):
    spec = importlib.util.spec_from_file_location(name, RECIPE / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FieldDriveJobTest:
    def test_refusal_preserves_exception_and_exact_log_tails(self, tmp_path, monkeypatch):
        module = load_job_module('job_loki')
        work = tmp_path / 'loki-runs/refused'
        work.mkdir(parents=True)
        for name in ('stdout.log', 'stderr.log', 'setup.in'):
            (work / name).write_bytes((name.encode() + b'\n') * 20000)
        error = RuntimeError('strict solver refusal')

        def refuse(self, bundle):
            raise error

        monkeypatch.setattr(module.LoKIDriver, 'run_prepared', refuse)
        module.install_refusal_capture(tmp_path)
        bundle = SimpleNamespace(spec={'scratch_root': str(work.parent)},
                                 job_name=work.name)
        with pytest.raises(RuntimeError) as caught:
            module.LoKIDriver.run_prepared(None, bundle)
        assert caught.value is error
        destination = tmp_path / 'results/loki-refusals/refused'
        for name in ('stdout.log', 'stderr.log', 'setup.in'):
            assert (destination / name).read_bytes() == (work / name).read_bytes()[-65536:]
        assert 'LOKI_REFUSAL_CAPTURED' in (tmp_path / 'results/steps.log').read_text()

    def test_success_keeps_return_value(self, tmp_path, monkeypatch):
        module = load_job_module('job_loki')
        rows = object()
        monkeypatch.setattr(module.LoKIDriver, 'run_prepared', lambda self, bundle: rows)
        module.install_refusal_capture(tmp_path)
        assert module.LoKIDriver.run_prepared(None, object()) is rows
        assert not (tmp_path / 'results').exists()

    def test_capture_failure_does_not_replace_solver_error(self, tmp_path, monkeypatch):
        module = load_job_module('job_loki')
        error = RuntimeError('original refusal')

        def refuse(self, bundle):
            raise error

        def broken_capture(*args):
            raise OSError('evidence destination unavailable')

        monkeypatch.setattr(module.LoKIDriver, 'run_prepared', refuse)
        monkeypatch.setattr(module, 'capture_logs', broken_capture)
        module.install_refusal_capture(tmp_path)
        bundle = SimpleNamespace(spec={'scratch_root': str(tmp_path / 'loki-runs')},
                                 job_name='refused')
        with pytest.raises(RuntimeError) as caught:
            module.LoKIDriver.run_prepared(None, bundle)
        assert caught.value is error

    @pytest.mark.parametrize('tamper', ['patch', 'source', 'expected_output'])
    def test_patch_identity_failure_never_writes_source(self, tmp_path, tamper):
        module = load_job_module('apply_loki_patch')
        original = (b'    // CIEff is <nu_eff>/N\n'
                    b'CIEffNew = mixingParameter * CIEffNew + (1 - mixingParameter) * CIEffOld;\n')
        changed = (b'    // CIEff is <nu_eff>/N\n'
                   b'CIEffNew = CIEffOld + mixingParameter * (CIEffNew - CIEffOld);\n')
        source = tmp_path / 'loki-src/source/ElectronKinetics.cpp'
        source.parent.mkdir(parents=True)
        source.write_bytes(original)
        patch = tmp_path / 'growth-mixing.patch'
        patch.write_bytes(b'test patch identity')
        digest = lambda data: hashlib.sha256(data).hexdigest()
        pin = {'patch_sha256': digest(patch.read_bytes()),
               'original_source_sha256': digest(original),
               'patched_source_sha256': digest(changed)}
        if tamper == 'patch':
            patch.write_bytes(b'changed patch')
        elif tamper == 'source':
            source.write_bytes(original + b'changed source')
        else:
            pin['patched_source_sha256'] = '0' * 64
        (tmp_path / 'growth-mixing-provenance.json').write_text(json.dumps(pin))
        before = source.read_bytes()
        with pytest.raises(ValueError):
            module.apply_growth_patch(tmp_path)
        assert source.read_bytes() == before


def load_scaling(monkeypatch):
    import sys
    job = load_job_module('job_loki')
    monkeypatch.setitem(sys.modules, 'job_loki', job)
    return load_job_module('node_scaling')


def test_node_budget_gate_uses_actual_measurements(monkeypatch):
    m = load_scaling(monkeypatch)
    fast = m.projection(41.34543256700272, 18.00126225402346, 100)
    assert fast['continue'] and fast['projected_total_hours'] < 16.8
    slow = m.projection(120, 100, 100)
    assert not slow['continue'] and slow['projected_total_hours'] > 16.8
    # Boundary includes elapsed setup/probe time and the remaining build reserve.
    overhead = m.projection(40, 20, 0)['projected_total_seconds']
    assert m.projection(40, 20, m.MAX_PROJECTED_SECONDS - overhead)['continue']
    assert not m.projection(40, 20, m.MAX_PROJECTED_SECONDS - overhead + 1)['continue']
    with pytest.raises(ValueError):
        m.projection(0, 1, 0)


def test_execution_equivalence_bounds_every_group_and_residual(monkeypatch):
    import copy
    import numpy as np
    m = load_scaling(monkeypatch)
    reference = {'eedf': np.array([[1., .2, 1e-30]]),
                 'rates': np.array([[1e-13, 1e-30]]),
                 'swarm': {'mobility': 1e24},
                 'power': {'Field': 1e-15, 'Power Balance': 1e-28}}
    candidate = copy.deepcopy(reference)
    candidate['eedf'][0, 1] *= 1 + 9.5e-13
    candidate['power']['Power Balance'] = -5e-28
    assert m.compare_values(reference, candidate)
    for group in reference:
        bad = copy.deepcopy(reference)
        if isinstance(bad[group], dict):
            key = next(iter(bad[group]))
            bad[group][key] *= 1.01
        else:
            bad[group] *= 1.01
        with pytest.raises(ValueError, match='equivalence'):
            m.compare_values(reference, bad)
    bad = copy.deepcopy(reference)
    bad['power']['Power Balance'] = 2e-20
    with pytest.raises(ValueError, match='residual'):
        m.compare_values(reference, bad)
    bad = copy.deepcopy(reference)
    bad['eedf'][0, 0] = np.nan
    with pytest.raises(ValueError, match='invalid numeric'):
        m.compare_values(reference, bad)


@pytest.mark.parametrize('refuse', [False, True])
def test_row_timing_preserves_results_and_refusals(tmp_path, monkeypatch, refuse):
    module = load_job_module('job_loki')
    error = RuntimeError('strict refusal')
    rows = [{'value': 1}]
    def run(self, bundle):
        if refuse:
            raise error
        return rows
    monkeypatch.setattr(module.LoKIDriver, 'run_prepared', run)
    module.install_row_timing(tmp_path)
    bundle = SimpleNamespace(job_name='timed', coordinate_values={'f_Ar4s': 1e-4},
                             fields_Td=(1.,), spec={'timeout_s': 6200})
    if refuse:
        with pytest.raises(RuntimeError) as caught:
            module.LoKIDriver.run_prepared(None, bundle)
        assert caught.value is error
    else:
        assert module.LoKIDriver.run_prepared(None, bundle) is rows
    start, end = [json.loads(line) for line in
                  (tmp_path / 'results/row-timings.jsonl').read_text().splitlines()]
    assert start['event'] == 'start' and end['event'] == 'end'
    assert start['coordinate'] == end['coordinate'] == {'f_Ar4s': 1e-4}
    assert start['started_utc'].endswith('Z') and end['ended_utc'].endswith('Z')
    assert end['status'] == ('refused' if refuse else 'passed')
    assert all(end[key] >= 0 for key in
               ('wall_seconds', 'child_user_seconds', 'child_system_seconds'))


@pytest.mark.parametrize('proceed', [False, True])
def test_pbs_budget_refusal_exits_cleanly_before_native_work(tmp_path, proceed):
    import os
    import subprocess
    import sys
    root = tmp_path / 'field-drive-job'
    root.mkdir()
    for name in ('job.pbs', 'run_table.py', 'job_loki.py', 'node_scaling.py',
                 'loki_preflight.py', 'apply_loki_patch.py',
                 'growth-mixing.patch', 'growth-mixing-provenance.json'):
        (root / name).write_bytes((RECIPE / name).read_bytes())
    (root / 'node_scaling.py').write_text(
        'import json\nfrom pathlib import Path\n'
        'Path("results/node-scaling.json").write_text(json.dumps('
        + repr({'continue': proceed, 'status': 'passed' if proceed else 'budget_refused'})
        + '))\n')
    pbs = (RECIPE / 'job.pbs').read_text()
    header = pbs[:pbs.index('mark environment\n')]
    decision = pbs[pbs.index('# First scientific execution:'):
                   pbs.index('step loki-single-coordinate-preflight')]
    script = tmp_path / 'test.pbs'
    script.write_text(header + decision + '\nstep native-build touch native-started\n')
    result = subprocess.run(
        ['bash', str(script)], env=dict(os.environ, PBS_O_WORKDIR=str(tmp_path),
                                      PATH=str(Path(sys.executable).parent) + ':' + os.environ['PATH']),
        capture_output=True, text=True, timeout=10)
    assert result.returncode == 0, result.stderr
    assert (root / 'native-started').exists() is proceed
    assert (root / 'results/exit_code').read_text().strip() == '0'
    assert (root / 'results/current-step').read_text().strip() == (
        'native-build' if proceed else 'node-scaling-refused')
    assert json.loads((root / 'results/node-scaling.json').read_text())['continue'] is proceed


def test_probe_reads_anisotropy_and_extra_rate_block(tmp_path, monkeypatch):
    m = load_scaling(monkeypatch)
    folder = tmp_path / 'output/solve/reducedField_30'
    folder.mkdir(parents=True)
    (folder / 'eedf.txt').write_text('energy eedf anisotropy\n1 .2 .3\n2 .4 .5\n')
    (folder / 'rateCoefficients.txt').write_text(
        'ine sup description\n1e-14 2e-14 first\n\n---------------------------\n'
        '* Extra Rate Coefficients *\n---------------------------\n3e-14 4e-14 extra\n')
    (folder / 'swarmParameters.txt').write_text('mobility = +1e24 unit\n')
    (folder / 'powerBalance.txt').write_text(
        'Field = +1e-15 unit\nPower Balance = -1e-28 unit\n'
        'Relative Power Balance = -3e-11 (%)\n')
    result = m.numeric_outputs(tmp_path)
    assert result['eedf'].shape == (2, 3)
    assert result['eedf'][1, 2] == .5
    assert result['rates'].shape == (2, 2)
    assert result['rates'][1, 1] == 4e-14
    assert result['power']['Power Balance'] == -1e-28
    assert 'Relative Power Balance' not in result['power']
