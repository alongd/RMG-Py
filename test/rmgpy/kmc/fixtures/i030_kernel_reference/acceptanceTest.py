"""Fast regressions for binding acceptance and the collected build entry."""
import importlib.util
import json
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


@pytest.mark.parametrize("optimized", [False, True], ids=["normal", "optimized"])
@pytest.mark.parametrize("fault", [
    "none", "additional_failure", "N4_N16_ratio", "N16_N4_ratio",
    "missing_row", "duplicate_row", "unknown_row", "known_now_passes",
    "unqualified", "unapproved", "additional_exception", "nonzero_model",
    "duplicate_exception", "loosen_record", "record_date",
])
def test_owner_ruling_is_exact_and_binding(fault, optimized):
    completed = run_python("""
import json
import checks
from pathlib import Path
cfg = json.loads(Path('reference/parameters.json').read_text())
report = json.loads(Path('candidate_results.json').read_text())
policy, digest = checks.read_policy(cfg)
fault = %r
if fault == 'additional_failure':
    report['comparison'][0]['within_error_budget'] = False
    report['comparison'][0]['accepted'] = False
elif fault.endswith('_ratio'):
    for row in report['comparison']:
        if row['case'] == fault[:-6] and row['pair_class'] == 'end/end':
            row['candidate_over_reference'] *= 1.001
elif fault == 'missing_row':
    report['comparison'].pop()
elif fault == 'duplicate_row':
    report['comparison'].append(report['comparison'][0])
elif fault == 'unknown_row':
    report['comparison'][0]['case'] = 'unapproved_case'
elif fault == 'known_now_passes':
    for row in report['comparison']:
        if row['case'] == 'N4_N16' and row['pair_class'] == 'end/end':
            row['within_error_budget'] = row['accepted'] = True
elif fault == 'unqualified':
    report['comparison'][0]['reference_qualified'] = False
elif fault == 'unapproved':
    policy['model_tolerance']['status'] = 'PROPOSED'
elif fault == 'additional_exception':
    policy['owner_ruling']['known_deviations'].append(dict(case='N8', pair_class='mid/mid'))
elif fault == 'nonzero_model':
    policy['model_tolerance']['relative_fraction'] = .01
elif fault == 'duplicate_exception':
    policy['owner_ruling']['known_deviations'][1] = policy['owner_ruling']['known_deviations'][0]
elif fault == 'loosen_record':
    policy['owner_ruling']['known_deviations'][0]['record_relative_tolerance'] = .1
elif fault == 'record_date':
    policy['owner_ruling']['known_deviations'][0]['owner_approved_on'] = 'unapproved'
checks.read_policy = lambda cfg: (policy, digest)
result = checks.assert_owner_approved_adoption(report, cfg)
if len(result['passing_rows']) != 16 or len(result['known_deviations']) != 2:
    raise RuntimeError('binding result lost its complete reported coverage')
print(json.dumps(result))
""" % fault, optimized=optimized)
    if fault == "none":
        assert completed.returncode == 0, completed.stderr
        record = json.loads(completed.stdout)
        assert record["binding_adoption_pass"]
        assert {row["case"] for row in record["known_deviations"]} == {"N4_N16", "N16_N4"}
    else:
        assert completed.returncode == 1, completed.stdout + completed.stderr
        assert "AssertionError:" in completed.stderr


def test_adoption_emits_and_retains_known_deviations(tmp_path, monkeypatch, capsys):
    # Exercise the real binding body using the committed reference as test input;
    # the campaign runner and complete preflight have their own live regressions.
    import shutil
    sys.path.insert(0, str(ROOT))
    import proposed_adoption_test
    import audit_sampler
    monkeypatch.setenv("MET_KERNEL_REFERENCE_PACK", str(ROOT))
    monkeypatch.setattr(audit_sampler, "run_audits", lambda cfg: json.loads((ROOT / "sampler_audits.json").read_text()))
    def committed_reference(args, check):
        output = Path(args[-1])
        output.mkdir()
        shutil.copyfile(ROOT / "reference/results/results.json", output / "results.json")
    monkeypatch.setattr(subprocess, "run", committed_reference)
    with pytest.warns(RuntimeWarning, match="OWNER-APPROVED KNOWN DEVIATION") as warnings:
        proposed_adoption_test.test_met_kernel_reference(tmp_path)
    assert len(warnings) == 2
    stdout = capsys.readouterr().out
    assert "N4_N16" in stdout and "N16_N4" in stdout
    record = json.loads((tmp_path / "owner_adoption.json").read_text())
    assert len(record["passing_rows"]) == 16
    assert len(record["known_deviations"]) == 2
    assert all(row["status"] == "reported known deviation; scientific comparison fails"
               for row in record["known_deviations"])
