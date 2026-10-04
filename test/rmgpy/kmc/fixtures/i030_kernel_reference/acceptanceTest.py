"""Fast regressions for binding acceptance and the collected build entry."""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parent
BUILD_ENTRY = ROOT.parent.parent / "kernelReferenceBuildTest.py"


def run_python(program, optimized=False):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[key] = "1"
    return subprocess.run(
        [sys.executable, "-B", *(["-O"] if optimized else []), "-c", program],
        cwd=ROOT, env=env, capture_output=True, text=True, timeout=60,
    )


def test_committed_failing_report_is_rejected_under_optimization():
    completed = run_python("""
import json
import checks
from pathlib import Path
report = json.loads(Path('candidate_results.json').read_text())
if report['adoption_pass']:
    raise RuntimeError('regression fixture must retain its failing report')
checks.assert_scientific_acceptance(report)
""", optimized=True)
    assert completed.returncode == 1, completed.stdout + completed.stderr
    assert "AssertionError: kernel exceeds numerical plus proposed model band" in completed.stderr


@pytest.mark.parametrize("optimized", [False, True], ids=["normal", "optimized"])
@pytest.mark.parametrize("gate", ["valid", "transport", "branch", "reference", "comparison"])
def test_acceptance_gates_remain_binding(gate, optimized):
    completed = run_python("""
import json
import checks
from pathlib import Path
report = json.loads(Path('candidate_results.json').read_text())
report['transport_pass'] = report['branch_pass'] = True
for row in report['comparison']:
    row['reference_qualified'] = row['within_error_budget'] = True
gate = %r
if gate in ('transport', 'branch'):
    report[gate + '_pass'] = False
elif gate == 'reference':
    report['comparison'][0]['reference_qualified'] = False
elif gate == 'comparison':
    report['comparison'][0]['within_error_budget'] = False
checks.assert_scientific_acceptance(report)
""" % gate, optimized=optimized)
    messages = {
        "transport": "target transport or Rg differs from literal inputs",
        "branch": "target min/floor/two-chain radius rule differs from literal contract",
        "reference": "Rouse oracle is not numerically qualified",
        "comparison": "kernel exceeds numerical plus proposed model band",
    }
    if gate == "valid":
        assert completed.returncode == 0, completed.stderr
    else:
        assert completed.returncode == 1, completed.stdout + completed.stderr
        assert "AssertionError: " + messages[gate] in completed.stderr


def test_live_bad_sampler_variance_is_rejected_under_optimization():
    completed = run_python("""
import numpy as np
import audit_sampler as audit
cfg, _ = audit.oracle.read_parameters()
cfg = dict(cfg, cases=[next(c for c in cfg['cases'] if c['name'] == 'N16')],
           pair_classes=['mid/mid'])
advance = audit.oracle.advance_modes
def double_noise(q, dt, rates, variances, rng):
    mean = q * np.exp(-dt[:, None] * rates)[..., None]
    return mean + 2 * (advance(q, dt, rates, variances, rng) - mean)
audit.oracle.advance_modes = double_noise
audit.stationary_variance_checks(cfg)
""", optimized=True)
    assert completed.returncode == 1, completed.stdout + completed.stderr
    assert "AssertionError: actual sampled stationary variance failed" in completed.stderr


def test_optimized_adoption_refuses_before_a_campaign(tmp_path):
    completed = run_python("""
import proposed_adoption_test
from pathlib import Path
proposed_adoption_test.test_met_kernel_reference(Path(%r))
""" % str(tmp_path), optimized=True)
    assert completed.returncode == 1, completed.stdout + completed.stderr
    assert "RuntimeError: MET kernel sampler qualification requires Python without optimization" in completed.stderr
    assert not list(tmp_path.iterdir())


def test_missing_build_fixture_fails_collection(tmp_path):
    copied = tmp_path / BUILD_ENTRY.name
    copied.write_bytes(BUILD_ENTRY.read_bytes())
    completed = subprocess.run(
        [sys.executable, "-B", "-m", "pytest", "-c", "/dev/null", "-p", "no:cacheprovider",
         "--collect-only", "-q", str(copied)],
        cwd=tmp_path, env=dict(os.environ, PYTEST_DISABLE_PLUGIN_AUTOLOAD="1"),
        capture_output=True, text=True, timeout=30,
    )
    assert completed.returncode == 2, completed.stdout + completed.stderr
    assert "FileNotFoundError: MET kernel reference build fixture is incomplete" in completed.stdout
    assert "1 error" in completed.stdout


def test_build_entry_uses_fixed_fixture(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("kernel_reference_build_entry", BUILD_ENTRY)
    entry = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(entry)
    observed = []
    class Adoption:
        @staticmethod
        def test_met_kernel_reference(output):
            observed.append((output, os.environ["MET_KERNEL_REFERENCE_PACK"]))
    monkeypatch.setenv("RMG_KMC_SLOW", "1")
    monkeypatch.setenv("MET_KERNEL_REFERENCE_PACK", str(tmp_path / "wrong"))
    monkeypatch.setattr(entry.importlib.util, "module_from_spec", lambda spec: Adoption)
    class Loader:
        @staticmethod
        def exec_module(module):
            pass
    class Spec:
        loader = Loader()
    monkeypatch.setattr(entry.importlib.util, "spec_from_file_location", lambda *args: Spec())
    entry.test_met_kernel_reference_build(tmp_path, monkeypatch)
    assert observed == [(tmp_path, str(ROOT))]
