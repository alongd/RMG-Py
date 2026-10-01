#!/usr/bin/env python3
"""
Tests for rmgpy.kmc.event_record module.
"""
import json
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
        bond_ops=[{"action": "break", "atoms": [1, 2]}, {"action": "form", "atoms": [2, 3]}],
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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])