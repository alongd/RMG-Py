"""Collected, opt-in build check for the finite-chain kernel fixture."""
import importlib.util
import os
from pathlib import Path

import pytest


FIXTURE = Path(__file__).resolve().parent / "fixtures/i030_kernel_reference"


def require_fixture():
    required = (
        "proposed_adoption_test.py", "checks.py", "audit_sampler.py",
        "reference/run.py", "reference/parameters.json", "acceptance_policy.json",
        "candidate_results.json", "reference/results/results.json",
    )
    missing = [str(FIXTURE / name) for name in required if not (FIXTURE / name).is_file()]
    if missing:
        raise FileNotFoundError("MET kernel reference build fixture is incomplete: " + ", ".join(missing))
    return FIXTURE


# Missing fixtures fail collection even when the expensive test is disabled.
require_fixture()


@pytest.mark.slow
def test_met_kernel_reference_build(tmp_path, monkeypatch):
    root = require_fixture()
    if os.environ.get("RMG_KMC_SLOW") != "1":
        pytest.skip("set RMG_KMC_SLOW=1 for the full MET kernel reference build check")
    monkeypatch.setenv("MET_KERNEL_REFERENCE_PACK", str(root))
    spec = importlib.util.spec_from_file_location("met_kernel_build_adoption", root / "proposed_adoption_test.py")
    adoption = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adoption)
    adoption.test_met_kernel_reference(tmp_path)
