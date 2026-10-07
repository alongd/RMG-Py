"""Behavioral contract for the persistent-carbene applicability policy."""

import copy
import json
from pathlib import Path

import pytest

from rmgpy.kmc.compiler import (
    PERSISTENT_CARBENE_POLICY_VERSION,
    apply_persistent_carbene_policy,
    classify_persistent_carbene,
    duplicate_transition_groups,
)


FIXTURES = Path(__file__).with_name("fixtures")


def _h_abstraction_record(event_id="evt_forward", reverse_of="evt_reverse"):
    return {
        "event_id": event_id,
        "reverse_of": reverse_of,
        "family": "H_Abstraction",
        "template": "synthetic",
        "status": "enabled",
        "rate_source": {"kind": "RMG family estimate", "available": True},
        "reactant_graphs": [
            """multiplicity 2
1 *1 C u1 p0 c0 {2,S} {3,S} {4,S}
2    H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
4    H u0 p0 c0 {1,S}
""",
            """multiplicity 2
1 *2 N u1 p1 c0 {2,S} {3,S}
2    H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
""",
        ],
        "product_graphs": [
            """multiplicity 3
1 *1 C u2 p0 c0 {2,S} {3,S}
2    H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
""",
            """1 *2 N u0 p1 c0 {2,S} {3,S} {4,S}
2    H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
4    H u0 p0 c0 {1,S}
""",
        ],
        "atom_map": {0: 0, 1: 1},
        "bond_ops": [
            {"action": "break", "atoms": [0, 3], "order": "1.0"},
            {"action": "form", "atoms": [4, 3], "order": "1.0"},
            {"action": "set_radical", "atom": 0, "value": 2},
            {"action": "set_radical", "atom": 4, "value": 0},
        ],
    }


def _pair():
    forward = _h_abstraction_record()
    reverse = copy.deepcopy(forward)
    reverse.update(
        event_id="evt_reverse",
        reverse_of="evt_forward",
        reactant_graphs=copy.deepcopy(forward["product_graphs"]),
        product_graphs=copy.deepcopy(forward["reactant_graphs"]),
    )
    return [forward, reverse]


def _root():
    return {
        "reactant_atom_index": 0,
        "recipe_label": "*1",
        "family_forward_role": "product",
    }


def _source(label="*1", role="product", persistent=True, source="training 629"):
    return [{
        "source": source,
        "mapped_roots": [{
            "recipe_label": label,
            "family_forward_role": role,
            "persistent_neutral_divalent_carbon": persistent,
        }],
    }]


def test_unsupported_persistent_carbene_refuses_both_directions():
    published, refusal = apply_persistent_carbene_policy(
        _pair(), _root(), _source(persistent=False)
    )
    assert published == []
    assert refusal["policy_version"] == PERSISTENT_CARBENE_POLICY_VERSION
    assert refusal["record_ids"] == ["evt_forward", "evt_reverse"]
    assert refusal["reason"] == "no selected rate contributor supports persistent u2 in the mapped role"
    assert refusal["mapped_root"]["recipe_label"] == "*1"
    assert refusal["rate_source"] == _source(persistent=False)


def test_supported_same_role_is_retained_with_transfer_caveat():
    published, refusal = apply_persistent_carbene_policy(
        _pair(), _root(), _source()
    )
    assert refusal is None
    assert len(published) == 2
    assert {record["reverse_of"] for record in published} == {
        "evt_forward", "evt_reverse"
    }
    annotations = [record["rate_source"]["applicability"] for record in published]
    assert all(item["disposition"] == "retained-supported-transfer" for item in annotations)
    assert all(item["caveat"] == "same reacting role; donor/acceptor environment may differ" for item in annotations)


def test_unserialized_source_domain_is_retained_as_unresolved():
    published, refusal = apply_persistent_carbene_policy(
        _pair(), _root(), None
    )
    assert refusal is None
    assert len(published) == 2
    assert all(
        record["rate_source"]["applicability"]["disposition"]
        == "retained-unresolved-applicability"
        for record in published
    )


def test_u2_contributor_in_wrong_recipe_role_does_not_authorize_transfer():
    published, refusal = apply_persistent_carbene_policy(
        _pair(), _root(), _source(label="*2", role="reactant")
    )
    assert published == []
    assert refusal["disposition"] == "refused-unsupported-transfer"


def test_nine_form_delocalised_witness_stays_enabled_and_identity_bound():
    witness = json.loads(
        (FIXTURES / "i078_resonance_witness.json").read_text()
    )
    decision = classify_persistent_carbene(
        witness,
        {"reactant_atom_index": 74, "recipe_label": "*3", "family_forward_role": "product"},
        _source(persistent=False),
    )
    assert witness["event_id"] == "evt_453febdaa6658ad2fa22f6375f21b3c46d0eb8ac92c534f69fb7dfafb72665b5"
    assert decision["disposition"] == "not-applicable"
    assert decision["mapped_root"]["product_rewrite_verified"] is True
    assert decision["mapped_root"]["resonance_form_count"] == 9
    assert decision["mapped_root"]["persistent_neutral_divalent_carbon"] is False


def test_ordinary_radical_channel_is_byte_unchanged():
    pair = _pair()
    before = json.dumps(pair, sort_keys=True)
    ordinary_root = {
        "reactant_atom_index": 4,
        "recipe_label": "*2",
        "family_forward_role": "product",
    }
    published, refusal = apply_persistent_carbene_policy(
        pair, ordinary_root, _source()
    )
    assert refusal is None
    assert json.dumps(published, sort_keys=True) == before


def test_retained_pair_remains_reciprocal():
    published, _ = apply_persistent_carbene_policy(_pair(), _root(), _source())
    by_id = {record["event_id"]: record for record in published}
    assert all(by_id[record["reverse_of"]]["reverse_of"] == record["event_id"] for record in published)


def _transition_pair(prefix, junction_kind=None, rate=1.0):
    forward, reverse = _pair()
    forward["event_id"], reverse["event_id"] = f"{prefix}_f", f"{prefix}_r"
    forward["reverse_of"], reverse["reverse_of"] = reverse["event_id"], forward["event_id"]
    for record in (forward, reverse):
        record["k_table"] = {"T": [1000.0], "k": [rate]}
        record["junction_ops"] = None
    if junction_kind:
        label = junction_kind.rsplit("_", 1)[-1]
        forward["junction_ops"] = [{
            "action": "create", "junction_kind": junction_kind,
            "attacked_atom_label": label,
            "chosen_resonance_localisation": f"ortho:{label}",
            "reflection": None if label == "S7" else {"S6": "S7", "S7": "S6"},
            "reverse_event_handle": reverse["event_id"],
        }]
        reverse["junction_ops"] = [{
            **forward["junction_ops"][0], "action": "dissociate",
            "reverse_event_handle": forward["event_id"],
        }]
    return [forward, reverse]


def test_s6_s7_split_rates_and_reverse_histories_are_not_duplicates():
    records = _transition_pair("s6", "J_ortho_S6", 0.5) + _transition_pair(
        "s7", "J_ortho_S7", 0.5
    )
    assert sum(record["k_table"]["k"][0] for record in records[::2]) == pytest.approx(1.0)
    assert duplicate_transition_groups(records) == []
    by_id = {record["event_id"]: record for record in records}
    assert all(by_id[record["reverse_of"]]["reverse_of"] == record["event_id"] for record in records)


def test_synthetic_duplicated_reciprocal_pair_is_detected():
    records = _transition_pair("first") + _transition_pair("duplicate")
    assert duplicate_transition_groups(records) == [
        [("duplicate_f", "duplicate_r"), ("first_f", "first_r")]
    ]
