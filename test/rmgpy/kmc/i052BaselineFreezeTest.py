"""Replay the immutable I-052 Step-1 legacy baseline."""

import json
import os

from i052_baseline_freeze import (
    FIXTURE,
    MANIFEST,
    compare_frozen,
    generate_observations,
    verify_manifest,
)


def test_manifest_hashes_are_immutable():
    verify_manifest(json.loads(MANIFEST.read_text()))


def test_replay_frozen_observations():
    expected = json.loads(FIXTURE.read_text())
    actual = generate_observations(
        remove_attachment=os.environ.get("RMG_KMC_BASELINE_NEGATIVE_CONTROL")
        == "remove_attachment"
    )
    compare_frozen(expected, actual)
