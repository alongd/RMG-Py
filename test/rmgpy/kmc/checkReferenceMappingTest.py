"""Policy tests for the declaration-generated reference mapping checker."""

import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).parents[3] / "scripts" / "check-kmc-reference-mapping.py"
SPEC = importlib.util.spec_from_file_location("check_kmc_reference_mapping", SCRIPT)
CHECKER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECKER)


def test_firing_statuses_are_not_named_refusals_even_with_a_reason():
    assert CHECKER.named_refusal({"status": "enabled", "status_reason": "diagnostic"}) is None
    assert CHECKER.named_refusal({"status": "irreversible",
                                  "status_reason": "RMG reaction is irreversible"}) is None
    assert CHECKER.named_refusal({"status": "refused", "status_reason": "policy"}) == \
        "refused:record_status:refused"


def test_numeric_proxy_depths_are_added_from_artifact_site_types():
    records = [{"site_type": "pristine"}, {"site_type": "end_radical@5"},
               {"site_type": "other@12"}, {"site_type": "not-a-depth@x"}]
    assert CHECKER.reference_dps(records, (1, 2, 3, 8)) == (1, 2, 3, 5, 8, 12)
