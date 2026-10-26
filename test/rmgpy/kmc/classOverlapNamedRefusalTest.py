"""A class overlap nested neither way is a named refusal, and so is what it takes with it (I-251).

Two loaded templates of one (family, direction) that make one transition at
the same weight, each also making transitions the other cannot (near the
chain ends of their witnesses), cannot be resolved: covering either side
loses transitions and loading both makes the shared ones twice.
``resolve_overlaps`` refuses the component ``refused:class_overlap``.  That
refusal is correct, so it is counted by name, not failed as a gate; a
component whose templates' weights disagree stays a gate.

Its reciprocal records then fail the cycle check for no reason of their own:
their reverse record loads no template, or the loaded reverse transition on
a witness was the refused template's.  On the PP catalog at 67eabc9fb these
were 11 "no loaded reverse transition" pairs and 3 one-way records.  A pair
whose every problem is gone once the overlap-refused templates are loaded
again is refused ``refused:class_overlap`` too, by name; any other problem
keeps it ``refused:cycle``.

The fixture is four H_Abstraction records of the 67eabc9fb PP shards, two
reciprocal pairs (evt_1a22/evt_d89a, evt_73f3/evt_c6b5), pruned to their base
templates: bee0 of evt_1a22 and 1fd6 of evt_73f3 overlap nested neither way
at weights 0 decades apart; evt_73f3 then loads nothing (one-way) and 38a2 of
evt_d89a finds no loaded reverse transition on its witness.
"""

import copy
import gzip
import json
import os
from pathlib import Path

import pytest

from rmgpy.kmc import template_emission as te

REPOSITORY = Path(__file__).resolve().parents[3]
FIXTURE = REPOSITORY / "test/rmgpy/kmc/fixtures/i251/neither_nested_cascade_payloads.json.gz"
DATABASE = Path(os.environ.get("RMG_DATABASE_PATH", "/nonexistent"))
needs_database = pytest.mark.skipif(
    not (DATABASE / "input/kinetics/families").is_dir(),
    reason="RMG_DATABASE_PATH is not a database checkout",
)
OVERLAP = ("1fd63df94565", "bee00028af2b")
PAIRS = {("evt_1a226b67a7a3", "evt_d89ab32993d5"), ("evt_73f3f571b4f5", "evt_c6b52c556d42")}


@pytest.fixture(scope="module")
def services():
    from rmgpy.kmc import emit
    return emit.load_services(str(DATABASE), ["H_Abstraction"])


@pytest.fixture(scope="module")
def payloads():
    with gzip.open(FIXTURE, "rt") as handle:
        return json.load(handle)


def _scaled(payloads, key, factor):
    out = copy.deepcopy(payloads)
    for p in out:
        for t in p["templates"]:
            if t.get("key", "").startswith(key):
                t["weight"]["w"] = [w * factor for w in t["weight"]["w"]]
    return out


def _overlaps(payloads):
    merged, _collisions = te.resolve_collisions(payloads)
    return te.resolve_overlaps(merged)


def _cycles(payloads, services):
    out, table = _overlaps(payloads)
    out = te.derive_weights(out, counter=services.counter)
    checked, cycles = te.check_cycles(out, counter=services.counter, kc=services.kc)
    return table, checked, {(r["forward"][:16], r["reverse"][:16]): r for r in cycles}


def test_a_neither_nested_overlap_at_one_weight_is_named_not_a_gate(payloads):
    _out, table = _overlaps(payloads)
    (row,) = [r for r in table if r["disposition"] != "covered"]
    assert sorted(t["key"][:12] for t in row["templates"]) == list(OVERLAP)
    assert row["disposition"] == "refused:class_overlap"
    assert "neither partition contains the other" in row["reason"]
    assert row["named"] is True
    assert te.overlap_gates(table) == []


def test_a_neither_nested_overlap_whose_weights_disagree_stays_a_gate(payloads):
    _out, table = _overlaps(_scaled(payloads, OVERLAP[1], 1.05))
    (row,) = [r for r in table if r["disposition"] != "covered"]
    assert row["disposition"] == "refused:class_overlap" and "decades" in row["reason"]
    assert row["named"] is False
    (gate,) = te.overlap_gates(table)
    assert "class overlap" in gate and "decades" in gate


@needs_database
def test_the_pairs_an_overlap_refusal_takes_with_it_are_named(payloads, services):
    _table, checked, cycles = _cycles(payloads, services)
    for pair in PAIRS:
        assert cycles[pair]["status"] == "refused:class_overlap", cycles[pair]["reason"]
    assert te.cycle_gates(list(cycles.values()), checked) == []
    statuses = {p["record_event_id"][:16]: (p.get("refusal") or {}).get("status")
                for p in checked if p.get("status") == "refused"}
    assert set(statuses) == {r for pair in PAIRS for r in pair}
    assert set(statuses.values()) == {"refused:class_overlap"}


@needs_database
def test_a_pair_with_a_problem_of_its_own_stays_a_cycle_gate(payloads, services):
    """32d4 (evt_d89a, in no overlap) is derived from d36a; its own kinetics
    moved 3 % off that derivation: pair evt_1a22/evt_d89a has a problem the
    overlap did not make, so it stays refused:cycle, a gate."""
    moved = copy.deepcopy(payloads)
    t = next(t for p in moved for t in p["templates"] if t["key"].startswith("32d4fbb08a52"))
    t["weight"]["w"] = [w * 1.03 for w in t["weight"]["w"]]
    out, table = _overlaps(moved)
    assert all("32d4fbb08a52" not in str(r["templates"]) for r in table)
    _table, checked, cycles = _cycles(moved, services)
    assert cycles[("evt_1a226b67a7a3", "evt_d89ab32993d5")]["status"] == "refused:cycle"
    assert te.cycle_gates(list(cycles.values()), checked)
