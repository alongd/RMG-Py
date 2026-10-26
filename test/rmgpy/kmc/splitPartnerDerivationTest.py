"""A loaded template that undoes two primaries is derived when they derive it alike (I-251).

``derive_weights`` derives a template's weight from its partner, the one loaded
template that undoes its transitions, by Kc exactly at T; a template with
several partners stayed independent (ln w linear in 1/T) and broke detailed
balance between grid points by the curvature of ln Kc.  Class coverage makes
such templates: on the PP catalog at 67eabc9fb, template 971c of evt_9bf4
(as_generated) covers both 39dc and 9e8e of evt_1d88, which at 448dc6247 were
derived one from 2735 and one from 46f7 of evt_b6fe.  971c then undoes both,
was never derived, and both reciprocal pairs (evt_1d88/evt_b6fe on 2735,
evt_7627/evt_9bf4 on 971c) were refused off the grid: S_f/S_r =
17.290370078603811 vs ssa ratio * Kc = 17.286109706722698 at 312.5 K.

The two primaries derive 971c alike (equal weight tables and counts; their
witnesses' Kc agree), so one derivation balances both.  Partners whose
derivations differ cannot be balanced by one loaded weight: the pair is a
named refusal (``refused:split_partners``), counted by name, not a gate.

The fixture is the four records' payloads of the 67eabc9fb PP shards, each
pruned to its base template and the templates of the split.
"""

import copy
import gzip
import json
import os
from pathlib import Path

import pytest

from rmgpy.kmc import template_emission as te

REPOSITORY = Path(__file__).resolve().parents[3]
FIXTURE = REPOSITORY / "test/rmgpy/kmc/fixtures/i251/split_partner_payloads.json.gz"
DATABASE = Path(os.environ.get("RMG_DATABASE_PATH", "/nonexistent"))
needs_database = pytest.mark.skipif(
    not (DATABASE / "input/kinetics/families").is_dir(),
    reason="RMG_DATABASE_PATH is not a database checkout",
)
SPLIT, PRIMARIES = "971cc825dd7f", ("273551e841c3", "46f75130660e")
PAIRS = {("evt_1d888f8e7224", "evt_b6fe3341ac28"), ("evt_7627ca42557e", "evt_9bf40ccf95c4")}


@pytest.fixture(scope="module")
def services():
    from rmgpy.kmc import emit
    return emit.load_services(str(DATABASE), ["H_Abstraction"])


@pytest.fixture(scope="module")
def payloads():
    with gzip.open(FIXTURE, "rt") as handle:
        return json.load(handle)


def _resolved(payloads):
    out, _table = te.resolve_collisions(payloads)
    return te.resolve_overlaps(out)


def _templates(out, status="emitted"):
    return {t["key"][:12]: t for p in out for t in p.get("templates", [])
            if t.get("status") == status}


def _cycles(payloads, services):
    out, _overlaps = _resolved(payloads)
    out = te.derive_weights(out, counter=services.counter)
    checked, table = te.check_cycles(out, counter=services.counter, kc=services.kc)
    return out, checked, {(r["forward"][:16], r["reverse"][:16]): r for r in table}


def test_the_fixture_is_the_split_partner_case(payloads):
    """971c covers both of evt_1d88's templates that 2735 and 46f7 undo; 2735
    and 46f7 cover evt_7627's 1c5b; all three load."""
    out, _overlaps = _resolved(payloads)
    covered = _templates(out, "covered")
    assert {r["key"][:12] for r in covered["39dcbd6991a0"]["covered_by"]} == {SPLIT}
    assert {r["key"][:12] for r in covered["9e8e722c7e39"]["covered_by"]} == {SPLIT}
    assert {r["key"][:12] for r in covered["1c5b2bcf0080"]["covered_by"]} == set(PRIMARIES)
    assert {SPLIT, *PRIMARIES} <= set(_templates(out))


@needs_database
def test_a_template_undoing_two_primaries_that_derive_it_alike_is_derived(payloads, services):
    out, checked, table = _cycles(payloads, services)
    split = _templates(out)[SPLIT]
    assert split.get("derived"), "971c undoes 2735 and 46f7 and stays an independent table"
    assert split["derived"]["primary"]["key"][:12] == min(PRIMARIES)
    assert {r["key"][:12] for r in split["derived"]["partners"]} == set(PRIMARIES)
    assert not any(t.get("derived") for k, t in _templates(out).items() if k in PRIMARIES)
    for pair in PAIRS:
        assert table[pair]["status"] == "ok", table[pair]["reason"]
    assert te.cycle_gates(list(table.values()), checked) == []


@needs_database
def test_primaries_that_derive_it_differently_are_a_named_refusal(payloads, services):
    """46f7's stored thermo form is moved by 10 J/mol on one product, so the
    two derivations of 971c differ by ~4e-3 at 300 K: no one weight balances
    both.  Both pairs are refused by name, and no gate fails for them."""
    moved = copy.deepcopy(payloads)
    t = next(t for p in moved for t in p["templates"] if t["key"].startswith(PRIMARIES[1]))
    t["thermo"]["form"]["products"][0]["H298"] += 10.0
    out, checked, table = _cycles(moved, services)
    split = _templates(out)[SPLIT]
    assert not split.get("derived")
    assert split["split_partners"]["status"] == "refused:split_partners"
    assert {r["key"][:12] for r in split["split_partners"]["partners"]} == set(PRIMARIES)
    for pair in PAIRS:
        assert table[pair]["status"] == "refused:split_partners", table[pair]["reason"]
    assert te.cycle_gates(list(table.values()), checked) == []
    refused = {(p["record_event_id"][:16], p["refusal"]["status"]) for p in checked
               if p.get("status") == "refused"}
    assert {r for r, _ in refused} == {r for pair in PAIRS for r in pair}
    assert {s for _, s in refused} == {"refused:split_partners"}


ONE_SIDED = REPOSITORY / "test/rmgpy/kmc/fixtures/i251/one_sided_partner_payloads.json.gz"


@needs_database
def test_a_template_whose_partner_shows_only_on_the_partner_s_witnesses_is_derived():
    """The partner relation is mutual, though a witness shows it one way only.

    On the PP catalog at 67eabc9fb (intra_H_migration), 77c1 of evt_bd1c covers
    0b44 and 8046, the reverses of dbb8 and 738f.  738f's witnesses show 77c1
    undoing it; 77c1's own witnesses show only dbb8.  So 738f (the lesser key)
    stayed primary against 77c1, 77c1 stayed primary against dbb8 (which was
    derived from it), and pair evt_6d5b/evt_df41 was refused off the grid:
    S_f/S_r = 17.290370078603431 vs 17.286109706723089 at 312.5 K.  The
    fixture is the six records of that component, each pruned to its base
    template and the templates of the component."""
    from rmgpy.kmc import emit
    services = emit.load_services(str(DATABASE), ["intra_H_migration"])
    with gzip.open(ONE_SIDED, "rt") as handle:
        payloads = json.load(handle)
    out, checked, table = _cycles(payloads, services)
    loaded = _templates(out)
    assert {r["key"][:12] for r in loaded["77c1d362bce2"]["derived"]["partners"]} == \
        {"738f06527b5a", "dbb8ec0a5629"}
    assert loaded["77c1d362bce2"]["derived"]["primary"]["key"][:12] == "738f06527b5a"
    assert not loaded["738f06527b5a"].get("derived")
    assert not loaded["dbb8ec0a5629"].get("derived")
    assert table and all(r["status"] == "ok" for r in table.values()), \
        [r["reason"] for r in table.values()]
    assert te.cycle_gates(list(table.values()), checked) == []
