#!/usr/bin/env python3
"""
Tests for rmgpy.kmc.event_record module.
"""
from dataclasses import FrozenInstanceError

import pytest

from rmgpy.kmc.event_record import EventRecord


def test_record_roundtrip():
    """Record should serialize to JSON and deserialize identically."""
    rec = EventRecord(
        event_id="evt_test123",
        family="H_Abstraction",
        template="template_1",
        arity=2,
        participant_site_types=["radical", "stable"],
        reactant_multiplicities=[2, 1],
        raw_path_degeneracy=3.0,
        ssa_multiplier=0.5,
        rate_order=2,
        rate_units="cm^3/mol/s",
        k_table={"T": [300.0, 500.0], "k": [1.0e10, 2.0e10]},
        atom_map={1: 1, 2: 2, 3: 3},
        bond_ops=[
            {"action": "break", "atoms": [1, 2]},
            {"action": "form", "atoms": [2, 3]},
        ],
        element_delta={"C": 0, "H": 0},
        radical_delta=0,
        status="enabled",
        status_reason="",
        orientation="as_generated",
        thermo_provenance={"method": "group_additivity"},
        provenance={"rmgpy_sha": "abc123", "database_sha": "def456"},
    )

    json_str = rec.to_json()
    rec2 = EventRecord.from_json(json_str)

    assert rec2.to_dict() == rec.to_dict()
    assert rec2.event_id == rec.event_id


def test_unknown_version_rejection():
    """from_dict should reject unknown schema versions."""
    data = {
        "schema_version": "rewrite_record/99.99",
        "event_id": "evt_test",
        "family": "Test",
        "template": "t",
        "arity": 1,
        "participant_site_types": [],
        "reactant_multiplicities": [],
        "raw_path_degeneracy": 1.0,
        "ssa_multiplier": 1.0,
        "rate_order": 0,
        "rate_units": "",
        "k_table": None,
        "atom_map": {},
        "bond_ops": [],
        "element_delta": {},
        "radical_delta": 0,
        "status": "enabled",
        "status_reason": "",
        "orientation": "as_generated",
        "thermo_provenance": {},
        "provenance": {},
        "cut_offset": None,
        "feature_ops": None,
        "inheritance": None,
        "junction_ops": None,
    }

    with pytest.raises(ValueError, match="Unknown schema version"):
        EventRecord.from_dict(data)


def test_missing_field_rejection():
    """from_dict should reject missing required fields."""
    data = {
        "schema_version": "rewrite_record/0.1",
        "event_id": "evt_test",
        "family": "Test",
        # missing many fields
    }

    with pytest.raises(ValueError, match="Missing required fields"):
        EventRecord.from_dict(data)


def test_unknown_field_rejection():
    """from_dict should reject unknown fields."""
    data = {
        "schema_version": "rewrite_record/0.1",
        "event_id": "evt_test",
        "family": "Test",
        "template": "t",
        "arity": 1,
        "participant_site_types": [],
        "reactant_multiplicities": [],
        "raw_path_degeneracy": 1.0,
        "ssa_multiplier": 1.0,
        "rate_order": 0,
        "rate_units": "",
        "k_table": None,
        "atom_map": {},
        "bond_ops": [],
        "element_delta": {},
        "radical_delta": 0,
        "status": "enabled",
        "status_reason": "",
        "orientation": "as_generated",
        "thermo_provenance": {},
        "provenance": {},
        "cut_offset": None,
        "feature_ops": None,
        "inheritance": None,
        "junction_ops": None,
        "unknown_field": "should_fail",
    }

    with pytest.raises(ValueError, match="Unknown fields"):
        EventRecord.from_dict(data)


def test_invalid_status():
    """validate() should reject invalid status."""
    rec = EventRecord(
        event_id="evt_test",
        family="Test",
        template="t",
        arity=1,
        participant_site_types=[],
        reactant_multiplicities=[],
        raw_path_degeneracy=1.0,
        ssa_multiplier=1.0,
        rate_order=0,
        rate_units="",
        k_table=None,
        atom_map={},
        bond_ops=[],
        element_delta={},
        radical_delta=0,
        status="invalid_status",
        status_reason="",
        orientation="as_generated",
        thermo_provenance={},
        provenance={},
    )

    with pytest.raises(ValueError, match="Invalid status"):
        rec.validate()


def test_invalid_orientation():
    """validate() should reject invalid orientation."""
    rec = EventRecord(
        event_id="evt_test",
        family="Test",
        template="t",
        arity=1,
        participant_site_types=[],
        reactant_multiplicities=[],
        raw_path_degeneracy=1.0,
        ssa_multiplier=1.0,
        rate_order=0,
        rate_units="",
        k_table=None,
        atom_map={},
        bond_ops=[],
        element_delta={},
        radical_delta=0,
        status="enabled",
        status_reason="",
        orientation="invalid_orientation",
        thermo_provenance={},
        provenance={},
    )

    with pytest.raises(ValueError, match="Invalid orientation"):
        rec.validate()


def test_k_table_validation():
    """validate() should check k_table structure."""
    rec = EventRecord(
        event_id="evt_test",
        family="Test",
        template="t",
        arity=1,
        participant_site_types=[],
        reactant_multiplicities=[],
        raw_path_degeneracy=1.0,
        ssa_multiplier=1.0,
        rate_order=0,
        rate_units="",
        k_table={"T": [300.0], "k": [1.0, 2.0]},  # mismatched lengths
        atom_map={},
        bond_ops=[],
        element_delta={},
        radical_delta=0,
        status="enabled",
        status_reason="",
        orientation="as_generated",
        thermo_provenance={},
        provenance={},
    )

    with pytest.raises(ValueError, match="same length"):
        rec.validate()

    # Missing keys
    rec2 = EventRecord(
        event_id="evt_test",
        family="Test",
        template="t",
        arity=1,
        participant_site_types=[],
        reactant_multiplicities=[],
        raw_path_degeneracy=1.0,
        ssa_multiplier=1.0,
        rate_order=0,
        rate_units="",
        k_table={"temperature": [300.0]},  # wrong key
        atom_map={},
        bond_ops=[],
        element_delta={},
        radical_delta=0,
        status="enabled",
        status_reason="",
        orientation="as_generated",
        thermo_provenance={},
        provenance={},
    )

    with pytest.raises(ValueError, match="must have 'T' and 'k' keys"):
        rec2.validate()


def test_strand_placeholders_null():
    """Strand-specific fields should be None by default."""
    rec = EventRecord(
        event_id="evt_test",
        family="Test",
        template="t",
        arity=1,
        participant_site_types=[],
        reactant_multiplicities=[],
        raw_path_degeneracy=1.0,
        ssa_multiplier=1.0,
        rate_order=0,
        rate_units="",
        k_table=None,
        atom_map={},
        bond_ops=[],
        element_delta={},
        radical_delta=0,
        status="enabled",
        status_reason="",
        orientation="as_generated",
        thermo_provenance={},
        provenance={},
    )

    assert rec.cut_offset is None
    assert rec.feature_ops is None
    assert rec.inheritance is None
    assert rec.junction_ops is None

    # They should serialize as null
    data = rec.to_dict()
    assert data["cut_offset"] is None
    assert data["feature_ops"] is None
    assert data["inheritance"] is None
    assert data["junction_ops"] is None


def test_event_id_deterministic():
    """event_id should be deterministic for same inputs."""
    rec1 = EventRecord(
        family="H_Abstraction",
        template="template_1",
        arity=2,
        participant_site_types=["radical", "stable"],
        reactant_multiplicities=[2, 1],
        raw_path_degeneracy=3.0,
        ssa_multiplier=0.5,
        rate_order=2,
        rate_units="cm^3/mol/s",
        k_table={"T": [300.0], "k": [1.0e10]},
        atom_map={1: 1, 2: 2},
        bond_ops=[{"action": "break", "atoms": [1, 2]}],
        element_delta={"C": 0, "H": 0},
        radical_delta=0,
        status="enabled",
        status_reason="",
        orientation="as_generated",
        thermo_provenance={},
        provenance={},
    )

    rec2 = EventRecord(
        family="H_Abstraction",
        template="template_1",
        arity=2,
        participant_site_types=["radical", "stable"],
        reactant_multiplicities=[2, 1],
        raw_path_degeneracy=3.0,
        ssa_multiplier=0.5,
        rate_order=2,
        rate_units="cm^3/mol/s",
        k_table={"T": [300.0], "k": [1.0e10]},
        atom_map={1: 1, 2: 2},
        bond_ops=[{"action": "break", "atoms": [1, 2]}],
        element_delta={"C": 0, "H": 0},
        radical_delta=0,
        status="enabled",
        status_reason="",
        orientation="as_generated",
        thermo_provenance={},
        provenance={},
    )

    assert rec1.event_id == rec2.event_id


def _semantic_record():
    methyl = (
        "multiplicity 2\n"
        "1 C u1 p0 c0 {2,S} {3,S} {4,S}\n"
        "2 H u0 p0 c0 {1,S}\n"
        "3 H u0 p0 c0 {1,S}\n"
        "4 H u0 p0 c0 {1,S}"
    )
    ethane = (
        "1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}\n"
        "2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}\n"
        "3 H u0 p0 c0 {1,S}\n"
        "4 H u0 p0 c0 {1,S}\n"
        "5 H u0 p0 c0 {1,S}\n"
        "6 H u0 p0 c0 {2,S}\n"
        "7 H u0 p0 c0 {2,S}\n"
        "8 H u0 p0 c0 {2,S}"
    )
    return EventRecord(
        family="R_Recombination",
        template="J_para",
        site_type="junction_radical+end_radical",
        arity=2,
        participant_site_types=["junction_radical", "end_radical"],
        reactant_multiplicities=[1, 1],
        raw_path_degeneracy=1.0,
        degeneracy=1.0,
        ssa_multiplier=1.0,
        rate_order=2,
        rate_units="m^3/(mol*s)",
        k_table={
            "T": [600.0],
            "k": [2.0],
            "interpolation": "linear-ln-k",
            "extrapolation": "refuse",
        },
        atom_map={0: 0, 1: 1},
        bond_ops=[{"action": "form", "atoms": [0, 1], "order": "1.0"}],
        element_delta={"C": 0},
        formula_delta={"C": 0},
        radical_delta=-2,
        status="enabled",
        status_reason="",
        orientation="as_generated",
        thermo_provenance={"gas_phase": "RMG"},
        provenance={"rmgpy_sha": "a" * 40},
        rate_source={"kind": "RMG family estimate"},
        reactant_graphs=[methyl, methyl],
        product_graphs=[ethane],
        rate_witness_reactant_graphs=[methyl, methyl],
        rate_witness_product_graphs=[ethane],
        proxy_padding={
            "status": "padded",
            "minimum_heavy_bond_distance": 15,
            "root_validation": {"method": "synthetic test fixture"},
            "executable_to_witness_projection": {
                "reactants": [
                    {
                        "participant_index": 0,
                        "executable_atom_index": 0,
                        "witness_atom_index": 0,
                    },
                    {
                        "participant_index": 1,
                        "executable_atom_index": 0,
                        "witness_atom_index": 0,
                    },
                ],
                "products": [
                    {
                        "participant_index": 0,
                        "executable_atom_index": 0,
                        "witness_atom_index": 0,
                    },
                    {
                        "participant_index": 0,
                        "executable_atom_index": 1,
                        "witness_atom_index": 1,
                    },
                ],
            },
        },
        coproducts=[],
        inventory_class="R1:J_ring",
        feature_ops=[{"action": "set_radical", "atom": 0, "value": 0}],
        inheritance={"policy": "template-local"},
        cut_offset=0,
        junction_ops=[{"junction_kind": "J_para"}],
    )


def test_event_id_is_full_sha256_and_record_is_frozen():
    record = _semantic_record()

    assert record.event_id.startswith("evt_")
    assert len(record.event_id) == 4 + 64
    int(record.event_id[4:], 16)
    record.validate()
    with pytest.raises(FrozenInstanceError):
        record.status = "refused"


def test_padded_record_requires_explicit_executable_projection():
    data = _semantic_record().to_dict()
    del data["proxy_padding"]["executable_to_witness_projection"]
    record = EventRecord.from_dict({**data, "event_id": ""})

    with pytest.raises(ValueError, match="two-sided executable projection"):
        record.validate()


@pytest.mark.parametrize(
    ("field_name", "tampered"),
    [
        ("inventory_class", "R1:quinoid_disproportionation"),
        ("junction_ops", [{"junction_kind": "J_ortho_S7"}]),
        ("reverse_of", "evt_" + "b" * 64),
        ("status", "irreversible"),
        ("status_reason", "tampered"),
        ("rate_source", {"kind": "tampered"}),
        (
            "k_table",
            {
                "T": [600.0],
                "k": [3.0],
                "interpolation": "linear-ln-k",
                "extrapolation": "refuse",
            },
        ),
        ("coproducts", [{"formula": {"H": 2}}]),
        ("feature_ops", []),
        ("inheritance", {"policy": "tampered"}),
        ("cut_offset", 1),
        ("degeneracy", 2.0),
        ("ssa_multiplier", 2.0),
        ("rate_witness_reactant_graphs", ["tampered reactant witness"]),
        ("rate_witness_product_graphs", ["tampered product witness"]),
        ("proxy_padding", {"status": "tampered"}),
        ("provenance", {"rmgpy_sha": "b" * 40}),
    ],
)
def test_validate_rejects_semantic_field_tampering(field_name, tampered):
    record = _semantic_record()
    data = record.to_dict()
    data[field_name] = tampered

    corrupted = EventRecord.from_dict(data)
    rehashed = EventRecord.from_dict({**data, "event_id": ""})
    assert rehashed.event_id != record.event_id
    with pytest.raises(ValueError):
        corrupted.validate()


def test_absent_witness_fields_preserve_legacy_event_id():
    record = EventRecord(family="legacy")
    data = record.to_dict()
    for name in (
        "rate_witness_reactant_graphs",
        "rate_witness_product_graphs",
        "proxy_padding",
    ):
        data.pop(name)

    recovered = EventRecord.from_dict(data)

    assert recovered.event_id == record.event_id
    recovered.validate()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
