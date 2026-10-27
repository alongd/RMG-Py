"""A class overlap nested neither way at one weight loads both sides, each shared transition once (I-271).

Two loaded templates of one (family, direction) that make one transition at
the same weight, each also making transitions the other cannot (near the
chain ends of their witnesses), cannot be resolved by covering: covering
either side loses transitions, and loading both as they are makes the shared
ones twice.  ``resolve_overlaps`` deduplicates the component: both sides
load, the templates of the reversed side defer to the as_generated templates
that share their transitions (``defers_to``), and a deferring template's
transition that a template it defers to also makes is not its own (the
matcher drops it).  A component whose templates' weights disagree stays a
gate.  (I-251 refused the component by name, and with it every cycle pair
the refusal left without its templates.)

The fixture is four H_Abstraction records of the 67eabc9fb PP shards, two
reciprocal pairs (evt_1a22/evt_d89a, evt_73f3/evt_c6b5), pruned to their base
templates: bee0 of evt_1a22 (as_generated) and 1fd6 of evt_73f3 (reversed)
overlap nested neither way at weights 0 decades apart.
"""

import copy
import gzip
import json
import os
from collections import Counter
from pathlib import Path

import pytest

from rmgpy.kmc import template_emission as te
from rmgpy.kmc import window_matcher as wm

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


def deduplicated_once(out, row, temperature):
    """On every stored witness of the component's templates (both sides):
    (transitions the loaded catalog makes more than once, transitions it
    loses against the union of what each template makes alone)."""
    members = {(t["record"], t["key"]) for t in row["templates"]}
    loaded = te.loaded_templates(out, temperature)
    alone = [wm.Template(t.name, t.family, t.atoms, list(t.bonds), t.edits, weight=t.weight)
             for t in loaded]
    twice, lost, states = 0, 0, 0
    for p in out:
        for t in p["templates"]:
            if (p["record_event_id"], t["key"]) not in members:
                continue
            for label in te.CHECK_LABELS:
                wit = ((t.get("maps") or {}).get("witnesses") or {}).get(label)
                if wit is None:
                    continue
                for side in ("r", "p"):
                    state = te.melt_state(te.rows_from_json(wit[side]))
                    made = Counter(tr.identity[1:] for tr in te.loaded_transitions(loaded, state))
                    union = {tr.identity[1:] for tr in te.loaded_transitions(alone, state)}
                    assert Counter(tr.identity[1:] for tr in wm.match(loaded, state).transitions) == made
                    twice += sum(1 for n in made.values() if n > 1)
                    lost += len(union - set(made))
                    states += 1
    assert states
    return twice, lost


def test_a_neither_nested_overlap_at_one_weight_loads_both_sides_once(payloads):
    out, table = _overlaps(payloads)
    (row,) = [r for r in table if r["disposition"] != "covered"]
    assert sorted(t["key"][:12] for t in row["templates"]) == list(OVERLAP)
    assert row["disposition"] == "deduplicated", row.get("reason")
    assert te.overlap_gates(table) == []
    status = {t["key"][:12]: t for p in out for t in p["templates"]}
    assert all(status[k]["status"] == "emitted" for k in OVERLAP)
    # the reversed side (1fd6) leaves the shared transitions to the as_generated one (bee0)
    assert [d["key"][:12] for d in status["1fd63df94565"].get("defers_to", ())] == ["bee00028af2b"]
    assert "defers_to" not in status["bee00028af2b"]
    assert deduplicated_once(out, row, 700.0) == (0, 0)


def test_a_neither_nested_overlap_whose_weights_disagree_stays_a_gate(payloads):
    _out, table = _overlaps(_scaled(payloads, OVERLAP[1], 1.05))
    (row,) = [r for r in table if r["disposition"] != "covered"]
    assert row["disposition"] == "refused:class_overlap" and "decades" in row["reason"]
    assert row["named"] is False
    (gate,) = te.overlap_gates(table)
    assert "class overlap" in gate and "decades" in gate


@needs_database
def test_the_reciprocal_pairs_of_a_deduplicated_overlap_close(payloads, services):
    _table, checked, cycles = _cycles(payloads, services)
    for pair in PAIRS:
        assert cycles[pair]["status"] == "ok", cycles[pair].get("reason")
    assert te.cycle_gates(list(cycles.values()), checked) == []
    assert not [p["record_event_id"][:16] for p in checked if p.get("status") == "refused"]


@needs_database
def test_a_pair_with_a_problem_of_its_own_stays_a_cycle_gate(payloads, services):
    """32d4 (evt_d89a, in no overlap) is derived from d36a; its own kinetics
    moved 3 % off that derivation: pair evt_1a22/evt_d89a has a problem of its
    own, so it is refused:cycle, a gate."""
    moved = copy.deepcopy(payloads)
    t = next(t for p in moved for t in p["templates"] if t["key"].startswith("32d4fbb08a52"))
    t["weight"]["w"] = [w * 1.03 for w in t["weight"]["w"]]
    out, table = _overlaps(moved)
    assert all("32d4fbb08a52" not in str(r["templates"]) for r in table)
    _table, checked, cycles = _cycles(moved, services)
    assert cycles[("evt_1a226b67a7a3", "evt_d89ab32993d5")]["status"] == "refused:cycle"
    assert te.cycle_gates(list(cycles.values()), checked)
