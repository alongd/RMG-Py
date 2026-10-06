import copy
import json
from pathlib import Path

import pytest

from rmgpy.kmc.equivalence import compare_artifacts
from rmgpy.molecule.molecule import Molecule


REAL_ARTIFACT = Path(
    "/home/alon/runs/i046-rules-from-training/cache/artifact/"
    "7491ed418f3ce2f0633688a9487cd12da278704b8f57e4f595176ddf80fb0296.json"
)
GRID = [500.0, 600.0]


def adjacency(smiles):
    return Molecule(smiles=smiles).to_adjacency_list(remove_h=False)


def record(
    event_id,
    reactants=None,
    products=None,
    atom_map=None,
    rates=None,
    degeneracy=2.0,
    status="enabled",
    reverse_of=None,
    coproducts=None,
):
    reactants = reactants or [adjacency("[CH3]"), adjacency("[OH]")]
    products = products or [adjacency("CO")]
    return {
        "event_id": event_id,
        "reactant_graphs": reactants,
        "product_graphs": products,
        "atom_map": atom_map or {0: 0, 1: 1},
        "bond_ops": [],
        "coproducts": coproducts or [],
        "status": status,
        "reverse_of": reverse_of,
        "rate_order": len(reactants),
        "rate_units": "m^3/(mol*s)" if len(reactants) == 2 else "s^-1",
        "reactant_pair_convention": "distinct-pair N_A*N_B",
        "ssa_multiplier": 1.0,
        "degeneracy": degeneracy,
        "raw_path_degeneracy": degeneracy,
        "k_table": {
            "T": GRID,
            "k": list(rates or [10.0, 20.0]),
            "interpolation": "linear-ln-k",
            "extrapolation": "refuse",
        },
    }


def artifact(records):
    return {
        "schema_version": "kmc_event_set/0.1",
        "inputs": {"temperature_grid": GRID},
        "records": records,
    }


def reversible_artifact():
    forward = record("forward", reverse_of="reverse")
    reverse = record(
        "reverse",
        reactants=forward["product_graphs"],
        products=forward["reactant_graphs"],
        atom_map={0: 0, 1: 1},
        rates=[2.0, 4.0],
        reverse_of="forward",
    )
    return artifact([forward, reverse])


@pytest.fixture(scope="module")
def real_pair():
    if not REAL_ARTIFACT.is_file():
        pytest.skip(f"external compiled-event fixture is absent: {REAL_ARTIFACT}")
    with REAL_ARTIFACT.open(encoding="utf-8") as handle:
        source = json.load(handle)
    by_id = {item["event_id"]: item for item in source["records"]}
    forward = next(
        item for item in source["records"] if item.get("reverse_of") in by_id
    )
    reverse = by_id[forward["reverse_of"]]
    result = {
        "schema_version": source["schema_version"],
        "inputs": {"temperature_grid": source["inputs"]["temperature_grid"]},
        "records": copy.deepcopy([forward, reverse]),
    }
    del source
    return result


def test_artifact_equals_itself():
    sample = reversible_artifact()
    report = compare_artifacts(sample, sample)
    assert report["equal"]
    assert report["left"]["channels"] == 2


def test_real_records_with_rehashed_ids_and_shuffled_order_compare_equal(real_pair):
    mutated = copy.deepcopy(real_pair)
    replacements = {
        item["event_id"]: f"replacement-{index}"
        for index, item in enumerate(mutated["records"])
    }
    for item in mutated["records"]:
        item["event_id"] = replacements[item["event_id"]]
        item["reverse_of"] = replacements[item["reverse_of"]]
    mutated["records"].reverse()
    assert compare_artifacts(real_pair, mutated)["equal"]


def test_permuted_participants_compare_equal():
    original = artifact([record("original")])
    permuted = copy.deepcopy(original)
    item = permuted["records"][0]
    item["reactant_graphs"].reverse()
    item["atom_map"] = {0: 1, 1: 0}
    assert compare_artifacts(original, permuted)["equal"]


def test_deleted_channel_is_reported():
    original = reversible_artifact()
    mutated = artifact([copy.deepcopy(original["records"][0])])
    report = compare_artifacts(original, mutated)
    assert not report["equal"]
    assert len(report["missing"]) == 1


def test_changed_rate_is_reported():
    original = artifact([record("original")])
    mutated = copy.deepcopy(original)
    mutated["records"][0]["k_table"]["k"][0] *= 1.1
    reasons = compare_artifacts(original, mutated)["changed"][0]["reasons"]
    assert {item["reason"] for item in reasons} >= {
        "summed channel rate",
        "summed physical propensity",
    }


def test_changed_degeneracy_is_reported():
    original = artifact([record("original")])
    mutated = copy.deepcopy(original)
    mutated["records"][0]["degeneracy"] = 3.0
    reasons = compare_artifacts(original, mutated)["changed"][0]["reasons"]
    assert "degeneracy" in {item["reason"] for item in reasons}


def test_flipped_reversibility_is_reported():
    original = reversible_artifact()
    mutated = copy.deepcopy(original)
    mutated["records"][0]["status"] = "irreversible"
    mutated["records"][0]["reverse_of"] = None
    reasons = compare_artifacts(original, mutated)["changed"][0]["reasons"]
    assert {item["reason"] for item in reasons} >= {
        "reversibility/status",
        "reverse linkage",
    }


def test_distinct_broken_reverse_linkages_are_reported():
    original = reversible_artifact()
    original["records"][1]["reverse_of"] = None
    mutated = copy.deepcopy(original)
    linked_half = mutated["records"][0]
    for key in ("degeneracy", "raw_path_degeneracy"):
        linked_half[key] /= 2.0
    linked_half["k_table"]["k"] = [value / 2.0 for value in linked_half["k_table"]["k"]]
    dangling_half = copy.deepcopy(linked_half)
    linked_half["event_id"] = "linked-half"
    dangling_half["event_id"] = "dangling-half"
    dangling_half["reverse_of"] = "missing-event"
    mutated["records"].append(dangling_half)

    report = compare_artifacts(original, mutated)
    reasons = report["changed"][0]["reasons"]
    assert "reverse linkage validity" in {item["reason"] for item in reasons}


def test_equivalent_split_channel_sums_propensity():
    original = artifact([record("original")])
    split = copy.deepcopy(original["records"][0])
    for key in ("degeneracy", "raw_path_degeneracy"):
        split[key] /= 2.0
    split["k_table"]["k"] = [value / 2.0 for value in split["k_table"]["k"]]
    other_half = copy.deepcopy(split)
    split["event_id"] = "split-a"
    other_half["event_id"] = "split-b"
    report = compare_artifacts(original, artifact([split, other_half]))
    assert report["equal"]
    assert report["left"]["records"] == 1
    assert report["right"]["records"] == 2


def test_overcounted_split_channel_is_reported():
    original = artifact([record("original")])
    duplicate = copy.deepcopy(original["records"][0])
    duplicate["event_id"] = "duplicate"
    report = compare_artifacts(original, artifact([original["records"][0], duplicate]))
    reasons = report["changed"][0]["reasons"]
    assert "summed physical propensity" in {item["reason"] for item in reasons}


def test_changed_product_fate_is_reported():
    carbon = adjacency("C")
    oxygen = adjacency("O")
    original_record = record(
        "fate",
        reactants=[adjacency("CO")],
        products=[carbon, oxygen],
        atom_map={0: 0, 1: 1},
        coproducts=[{"formula": {"H": 2, "O": 1}}],
    )
    original = artifact([original_record])
    mutated = copy.deepcopy(original)
    mutated["records"][0]["coproducts"] = [{"formula": {"C": 1, "H": 4}}]
    reasons = compare_artifacts(original, mutated)["changed"][0]["reasons"]
    assert "product fate" in {item["reason"] for item in reasons}
