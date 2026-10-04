"""Owner-approved binding test, directly runnable from the repository.

This file is the exact test embedded in pack.md and executed by the Verifier.
"""


def test_met_kernel_reference(tmp_path):
    import importlib.util
    import json
    import os
    from pathlib import Path
    import subprocess
    import sys
    import warnings
    sys.dont_write_bytecode = True
    root = Path(os.environ.get("MET_KERNEL_REFERENCE_PACK", str(Path(__file__).resolve().parent)))
    if not (root / "reference/run.py").is_file():
        import pytest
        pytest.skip("MET kernel reference fixture is not provisioned at " + str(root))
    sys.path.insert(0, str(root))
    import checks
    import audit_sampler
    spec = importlib.util.spec_from_file_location("test_rouse_inputs", root / "reference/run.py")
    oracle = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(oracle)
    cfg, digest = oracle.read_parameters()
    # Complete preflight, including innovation-noise stationary variance.
    audits = audit_sampler.run_audits(cfg)
    (tmp_path / "sampler_audits.json").write_text(json.dumps(audits, indent=2) + "\n")
    output = tmp_path / "fresh-reference"
    subprocess.run([
        sys.executable, "-B", "-u", str(root / "reference/run.py"),
        "--workers", os.environ.get("MET_KERNEL_REFERENCE_WORKERS", "8"),
        "--output", str(output),
    ], check=True)
    reference = json.loads((output / "results.json").read_text())
    assert reference["parameters_sha256"] == digest
    report = checks.target_checks(checks.load_target(), cfg, reference)
    (tmp_path / "candidate_results.json").write_text(json.dumps(report, indent=2) + "\n")
    adoption = checks.assert_owner_approved_adoption(report, cfg)
    (tmp_path / "owner_adoption.json").write_text(json.dumps(adoption, indent=2) + "\n")
    for deviation in adoption["known_deviations"]:
        message = "OWNER-APPROVED KNOWN DEVIATION: " + json.dumps(deviation, sort_keys=True)
        print(message, flush=True)
        warnings.warn(message, RuntimeWarning)
