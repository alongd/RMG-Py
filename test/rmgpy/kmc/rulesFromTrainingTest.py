"""Fast preparation and cached real-artifact regressions (never compile)."""

import copy
import importlib.util
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from rmgpy.kmc.compiler import prepare_rate_rules
from rmgpy.rmg import input as rmg_input


@pytest.mark.parametrize("excluded", [False, True])
def test_rmg_preparation_conditions_constraints_and_reuse(monkeypatch, excluded):
    calls = []
    constraints = {"maximumCarbonAtoms": 1}
    polymer = {"maximumCarbonAtoms": 2}
    owner = SimpleNamespace(species_constraints=constraints, polymer_constraints=polymer)
    monkeypatch.setattr(rmg_input, "rmg", owner)

    def add(*, thermo_database):
        assert owner.species_constraints == {}
        assert owner.polymer_constraints is None
        assert thermo_database == "thermo"
        calls.append("training")

    def fill(*, verbose):
        assert owner.species_constraints is constraints
        assert owner.polymer_constraints is polymer
        assert verbose is True
        calls.append("average")

    family = SimpleNamespace(auto_generated=False, rules=SimpleNamespace(entries={}),
                             add_rules_from_training=add, fill_rules_by_averaging_up=fill)
    tree = SimpleNamespace(auto_generated=True)
    database = SimpleNamespace(families={"ordinary": family, "tree": tree})
    options = {"kinetics_depositories": ["!training"] if excluded else ["training"]}
    first = prepare_rate_rules(database, "thermo", **options)
    assert calls == (["average"] if excluded else ["training", "average"])
    assert prepare_rate_rules(database, "thermo", **options) == first
    assert owner.species_constraints is constraints and owner.polymer_constraints is polymer
    assert vars(tree) == {"auto_generated": True}
    assert first["families"]["ordinary"]["training_rules_added"] is not excluded


def test_training_failure_restores_global_constraints(monkeypatch):
    owner = SimpleNamespace(species_constraints={"limit": 1}, polymer_constraints={"limit": 2})
    original = vars(owner).copy()
    monkeypatch.setattr(rmg_input, "rmg", owner)

    def fail(**kwargs):
        raise RuntimeError("training failed")

    family = SimpleNamespace(auto_generated=False, rules=SimpleNamespace(entries={}),
                             add_rules_from_training=fail)
    with pytest.raises(RuntimeError, match="training failed"):
        prepare_rate_rules(SimpleNamespace(families={"ordinary": family}), "thermo")
    assert vars(owner) == original
    assert rmg_input.rmg is owner
    assert not hasattr(family, "_kmc_rate_rule_preparation")


def test_offline_preparation_restores_absent_global_input(monkeypatch):
    monkeypatch.setattr(rmg_input, "rmg", None)
    family = SimpleNamespace(
        auto_generated=False, rules=SimpleNamespace(entries={}),
        add_rules_from_training=lambda **kwargs: None,
        fill_rules_by_averaging_up=lambda **kwargs: None,
    )
    prepare_rate_rules(SimpleNamespace(families={"ordinary": family}), "thermo")
    assert rmg_input.rmg is None


def load_probe():
    path = Path(__file__).with_name("fixtures") / "i046_probe/compare_artifacts.py"
    spec = importlib.util.spec_from_file_location("i046_compare", path)
    probe = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(probe)
    return probe


def test_comparison_matches_chemistry_across_changed_nested_inverse_ids():
    probe = load_probe()
    old = {
        "family": "R_Recombination", "event_id": "evt_" + "1" * 64,
        "junction_ops": [{"action": "create", "junction_kind": "J_para",
                          "reverse_event_handle": "evt_" + "2" * 64}],
        "rate_source": {"source": "old"}, "k_table": {"k": [1.0]},
    }
    new = copy.deepcopy(old)
    new["event_id"] = "evt_" + "3" * 64
    new["junction_ops"][0]["reverse_event_handle"] = "evt_" + "4" * 64
    new["rate_source"] = {"source": "new"}
    new["k_table"] = {"k": [2.0]}
    assert probe.structural_key(old) == probe.structural_key(new)
    new["junction_ops"][0]["junction_kind"] = "J_ortho_S8"
    assert probe.structural_key(old) != probe.structural_key(new)


def test_tree_invariance_allows_added_channels_and_corrected_site_labels():
    probe = load_probe()
    old = {"records": [
        {"family": family, "site_type": "old", "participant_site_types": ["end_radical"],
         "k_table": {"k": [1]}, "rate_source": {"source": "pinned"}}
        for family in probe.TREE_FAMILIES
    ]}
    new = copy.deepcopy(old)
    new["provenance"] = {"rmg_database_sha": "pinned"}
    for record in new["records"]:
        record["site_type"] = "benzylic_context"
        record["participant_site_types"] = ["interior_radical"]
    new["records"].append({"family": "R_Recombination", "bond_ops": [{"action": "new"}],
                           "k_table": {"k": [2]}, "rate_source": {"source": "new"}})
    assert probe.verify_tree_invariance(old, new) == {family: 1 for family in probe.TREE_FAMILIES}
    new["records"][0]["k_table"]["k"] = [3]
    with pytest.raises(AssertionError):
        probe.verify_tree_invariance(old, new)


def test_tree_invariance_rejects_rate_and_non_provenance_source_mutations():
    probe = load_probe()
    old = {"records": [
        {"family": family, "k_table": {"T": [600.0], "k": [1.0]},
         "rate_source": {"kind": "RMG family estimate",
                          "reference_thermo": {"rmg_database_sha": "pinned"}}}
        for family in probe.TREE_FAMILIES
    ]}
    new = copy.deepcopy(old)
    new["provenance"] = {"rmg_database_sha": "pinned"}

    new["records"][0]["k_table"]["k"][0] *= 1.0 + 1.0e-9
    with pytest.raises(AssertionError):
        probe.verify_tree_invariance(old, new)

    new = copy.deepcopy(old)
    new["provenance"] = {"rmg_database_sha": "pinned"}
    new["records"][0]["rate_source"]["kind"] = "mutated"
    with pytest.raises(AssertionError):
        probe.verify_tree_invariance(old, new)


def test_tree_invariance_requires_database_provenance_pin():
    probe = load_probe()
    old = {"records": [
        {"family": family, "k_table": {"T": [600.0], "k": [1.0]},
         "rate_source": {"kind": "RMG family estimate"}}
        for family in probe.TREE_FAMILIES
    ]}
    new = copy.deepcopy(old)
    new["provenance"] = {}

    with pytest.raises(AssertionError, match="provenance.rmg_database_sha"):
        probe.verify_tree_invariance(old, new)


@pytest.fixture(scope="module")
def cached_artifacts():
    supplied = os.environ.get("RMG_KMC_ARTIFACT")
    if not supplied:
        pytest.skip("set RMG_KMC_ARTIFACT to the single compiled artifact")
    probe = load_probe()
    return probe, json.loads(probe.OLD_ARTIFACT.read_text()), json.loads(Path(supplied).read_text())


@pytest.mark.parametrize("family, root", [
    ("R_Addition_MultipleBond", "R_R;YJ"),
    ("H_Abstraction", "X_H_or_Xrad_H_Xbirad_H_Xtrirad_H;Y_rad_birad_trirad_quadrad"),
    ("intra_H_migration", "RnH;Y_rad_out;XH_out"),
])
def test_compiled_representative_selects_training_derived_rules(cached_artifacts, family, root):
    probe, old, new = cached_artifacts
    record = probe.representative(old, new, family)
    source = record["rate_source"]
    assert probe.has_training(source), (family, source)
    assert not probe.has_default_or_root(source, root), (family, source)
    assert source.get("comment"), (family, "missing RMG estimation comment")
    assert source["template"] == record["template"]


def test_tree_family_rate_tables_and_sources_are_bit_identical(cached_artifacts):
    probe, old, new = cached_artifacts
    assert probe.verify_tree_invariance(old, new) == {
        "Disproportionation": 12680, "R_Recombination": 594,
    }
