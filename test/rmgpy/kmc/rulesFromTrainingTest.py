"""Fast preparation and cached real-artifact regressions (never compile)."""

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


@pytest.fixture(scope="module")
def cached_artifacts():
    supplied = os.environ.get("RMG_KMC_ARTIFACT")
    if not supplied:
        pytest.skip("set RMG_KMC_ARTIFACT to the single compiled artifact")
    path = Path(__file__).with_name("fixtures") / "i046_probe/compare_artifacts.py"
    spec = importlib.util.spec_from_file_location("i046_compare", path)
    probe = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(probe)
    return probe, json.loads(probe.OLD_ARTIFACT.read_text()), json.loads(Path(supplied).read_text())


@pytest.mark.parametrize("family", ["R_Addition_MultipleBond", "H_Abstraction", "intra_H_migration"])
def test_compiled_representative_selects_training_derived_rules(cached_artifacts, family):
    probe, old, new = cached_artifacts
    record = probe.representative(old, new, family)
    source = record["rate_source"]
    assert source.get("comment"), (family, "missing RMG estimation comment")
    assert probe.has_training(source), (family, source)
    assert not probe.has_default(source), (family, source)
    assert source["template"] == record["template"]


def test_tree_family_rate_tables_and_sources_are_bit_identical(cached_artifacts):
    probe, old, new = cached_artifacts
    assert probe.verify_tree_invariance(old, new) == {
        "Disproportionation": 12680, "R_Recombination": 594,
    }
