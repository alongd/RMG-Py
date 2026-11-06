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
