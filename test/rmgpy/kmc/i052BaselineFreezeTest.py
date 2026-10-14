"""Replay the immutable I-052 Step-1 legacy baseline.

The baseline pins absolute paths, historical commits and the exact
environment of the machine that froze it; where those inputs are absent the
tests skip with the reason, and where they are present they stay strict.
"""

import json
import os

import pytest

from i052_baseline_freeze import (
    FIXTURE,
    MANIFEST,
    baseline_skip_reason,
    compare_frozen,
    generate_observations,
    verify_manifest,
)


def _skip_unless_inputs_present(manifest):
    reason = baseline_skip_reason(manifest)
    if reason is not None:
        pytest.skip(f"frozen I-052 baseline not checkable here: {reason}")


def test_manifest_hashes_are_immutable():
    manifest = json.loads(MANIFEST.read_text())
    _skip_unless_inputs_present(manifest)
    verify_manifest(manifest)


def test_replay_frozen_observations():
    _skip_unless_inputs_present(json.loads(MANIFEST.read_text()))
    expected = json.loads(FIXTURE.read_text())
    actual = generate_observations(
        remove_attachment=os.environ.get("RMG_KMC_BASELINE_NEGATIVE_CONTROL")
        == "remove_attachment"
    )
    compare_frozen(expected, actual)
