"""Loaded templates that make one transition twice (I-210).

The fixture is the six H_Abstraction payloads of the three PE reciprocal pairs
the I-199 PE catalog refused (``fixtures/i210_pe_exchange_payloads.zlib.b64``,
copied unchanged from ``pe.rows.jsonl``): a secondary radical and a secondary
C-H exchange their H, both sites at least five carbons from a chain end.  The
balance checks run with the hydrogen-leaf quotient (``loaded_transitions``),
which these tests hold equal to the matcher's full enumeration.
"""

import base64
import json
import zlib
from pathlib import Path

import pytest

from rmgpy.kmc import template_emission as te
from rmgpy.kmc import window_matcher as wm

HERE = Path(__file__).resolve().parent
DATABASE = Path("/home/alon/runs/phase2b/database")
EXCHANGE = ("3e484af0e1", "53b80ec7cc", "a7562b8ed1", "247abac61b")


def exchange_payloads():
    path = HERE / "fixtures" / "i210_pe_exchange_payloads.zlib.b64"
    return json.loads(zlib.decompress(base64.b64decode(path.read_bytes())))


def bucket(payloads, keys=EXCHANGE):
    """The loaded H_Abstraction templates whose key starts with one of ``keys``,
    as matcher templates named by their key prefix."""
    loaded = te._Loaded(payloads)
    entries = [t for t in loaded.by[("H_Abstraction", "forward")] if t["key"][:10] in keys]
    return [te.matcher_from_payload(dict(t, family="H_Abstraction"), 0.0, weight=1.0,
                                    name=t["key"][:10]) for t in entries]


def witness_state(payloads, key, side):
    t = next(t for p in payloads for t in p["templates"] if t["key"].startswith(key))
    return te.melt_state(te.rows_from_json(t["maps"]["witnesses"]["emitted"][side]))


def transition_fields(tr):
    """Everything a balance check reads of a transition (not the embedding
    count, nor the mapping of unedited hydrogen leaves)."""
    return (tr.template, tr.family, tr.class_index, tr.edit_key, tr.components, tr.edits,
            tr.form_bonds, tr.form_radicals, tr.identity, tr.delta, tr.weight,
            tr.state_revision, tr.state_hash)


@pytest.mark.parametrize("key,side", [("3e484af0e1", "r"), ("53b80ec7cc", "r")])
def test_loaded_transitions_equal_the_matchers_on_pe_witnesses(key, side):
    payloads = te.resolve_collisions(exchange_payloads())[0]
    templates = bucket(payloads)
    state = witness_state(payloads, key, side)
    full = wm.match(templates, state).transitions
    quotient = te.loaded_transitions(templates, state)
    assert len(full) > 0
    assert [transition_fields(x) for x in quotient] == [transition_fields(x) for x in full]
    assert sum(x.embeddings for x in quotient) < sum(x.embeddings for x in full)


# --------------------------------------------------------------------------
# the refusal and its cause
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def services():
    if not DATABASE.is_dir():
        pytest.skip(f"needs the RMG database at {DATABASE}")
    from rmgpy.kmc import emit
    return emit.load_services(str(DATABASE), ["H_Abstraction"])


def collided():
    return te.resolve_collisions(exchange_payloads())[0]


def keys_of(refs):
    return sorted(r["key"][:10] for r in refs)


def counts_in(reason):
    import re
    return [int(n) for group in re.findall(r"counts \[([0-9, ]+)\]", reason)
            for n in group.split(",")]


def test_unresolved_exchange_pairs_are_refused_two_undoing_transitions_each(services):
    """The PE catalog's refusal, reproduced from its rows: every forward
    transition of the three pairs is undone by two loaded reverse transitions."""
    payloads = te.derive_weights(collided(), counter=services.counter)
    _out, table = te.check_cycles(payloads, counter=services.counter, kc=services.kc,
                                  complete=True)
    assert len(table) == 3
    assert [r["status"] for r in table] == ["refused:cycle"] * 3
    for row in table:
        counts = counts_in(row["reason"])
        assert row["bijection_failures"] > 0 and counts and set(counts) == {2}, row["reason"]


def test_the_two_undoing_transitions_are_one_edit_made_by_two_templates():
    """Not two symmetry images, not two copies of one row: one physical
    transition (same components, same product rows) made by the reversed
    record's base template and by an as_generated record's context template."""
    payloads = collided()
    loaded = te._Loaded(payloads)
    entries = loaded.by[("H_Abstraction", "forward")]
    templates = [te.matcher_from_payload(dict(t, family="H_Abstraction"), 0.0, weight=1.0,
                                         name=f"L{k}") for k, t in enumerate(entries)]
    rows = te.rows_from_json(next(t for p in payloads for t in p["templates"]
                                  if t["key"].startswith("3e484af0e1")
                                  and t["status"] == "emitted")["maps"]["witnesses"]["emitted"]["r"])
    state = te.melt_state(rows)
    first = te.first_form_rows(rows)
    forward = [tr for tr in te.loaded_transitions(templates, state)
               if entries[int(tr.template[1:])]["key"].startswith("3e484af0e1")]
    assert forward
    pairs = set()
    for tr in forward:
        touched = wm.product_rows(state, tr)
        after = dict(rows)
        after.update(touched)
        check = te.melt_state(after)
        back = [rt for rt in te.loaded_transitions(templates, check)
                if set(wm.product_rows(check, rt)) == set(touched)
                and all(te._norm_row(wm.product_rows(check, rt)[a]) == te._norm_row(first[a])
                        for a in touched)]
        assert len(back) == 2
        assert len({rt.identity[1:] for rt in back}) == 1
        owners = sorted((loaded.owner[id(entries[int(rt.template[1:])])][1][:10],
                         next(p["orientation"] for p in payloads if p["record_event_id"] ==
                              loaded.owner[id(entries[int(rt.template[1:])])][0]))
                        for rt in back)
        assert sorted(o for _k, o in owners) == ["as_generated", "reversed"]
        pairs.add(tuple(k for k, _o in owners))
    assert ("53b80ec7cc", "a7562b8ed1") in pairs
    assert {k for pair in pairs for k in pair} <= set(EXCHANGE)


# --------------------------------------------------------------------------
# the fix: overlap resolution at assembly
# --------------------------------------------------------------------------
def test_overlap_census_is_two_partitions_of_one_region():
    """One overlap component, the as_generated side {3e48, a756} and the
    reversed side {53b8, 247a}; each co-made transition is made by one template
    of each side, the weights agree, the reversed side is covered."""
    out, table = te.resolve_overlaps(collided())
    assert len(table) == 1
    row = table[0]
    assert (row["family"], row["direction"]) == ("H_Abstraction", "forward")
    assert row["disposition"] == "covered", row
    assert keys_of(row["covered"]) == ["247abac61b", "53b80ec7cc"]
    assert keys_of(row["covered_by"]) == ["3e484af0e1", "a7562b8ed1"]
    assert row["transitions"] > 0 and row["log10_diameter"] <= 1e-12
    statuses = {}
    for p in out:
        for t in p["templates"]:
            statuses.setdefault(t["key"][:10], set()).add(t["status"])
    assert statuses["53b80ec7cc"] == statuses["247abac61b"] == {"covered"}
    assert statuses["3e484af0e1"] == statuses["a7562b8ed1"] == {"emitted", "merged"}
    emitted = {(p["record_event_id"], t["key"]) for p in out for t in p["templates"]
               if t["status"] == "emitted"}
    for p in out:
        for t in p["templates"]:
            if t["status"] == "covered":
                assert keys_of(t["covered_by"]) == ["3e484af0e1", "a7562b8ed1"]
                assert {(c["record"], c["key"]) for c in t["covered_by"]} <= emitted
                assert "merged_into" not in t
    assert {p["status"] for p in out} == {"emitted"}


def test_exchange_pairs_balance_after_overlap_resolution(services):
    """Every reciprocal pair checked, every pair ok, no gate; the covered
    templates' witnesses are still checked rows; neither covered template loads.
    a756 undoes its own transitions (self-inverse, primary) and 3e48's, which
    nothing else undoes: 3e48 is derived from it (I-251: the partners of
    ``derive_weights`` are mutual; before, 3e48 stayed an independent table)."""
    before = te.derive_weights(collided(), counter=services.counter)
    _o, unresolved = te.check_cycles(before, counter=services.counter, kc=services.kc,
                                     complete=True)
    resolved = te.resolve_overlaps(collided())[0]
    payloads = te.derive_weights(resolved, counter=services.counter)
    out, table = te.check_cycles(payloads, counter=services.counter, kc=services.kc,
                                 complete=True)
    assert len(table) == 3 and [r["status"] for r in table] == ["ok"] * 3, table
    assert te.cycle_gates(table, out) == []
    assert sum(r["rows"] for r in table) == sum(r["rows"] for r in unresolved)
    assert max(r["max_ratio_error"] for r in table) <= te.CYCLE_RTOL
    loaded = te.loaded_templates(out, 700.0)
    assert sorted(t.name.split(":")[1][:10] for t in loaded
                  if t.name.split(":")[1][:10] in EXCHANGE) == ["3e484af0e1", "a7562b8ed1"]
    exchange = {t["key"][:10]: t for p in out for t in p["templates"]
                if t["status"] == "emitted" and t["key"][:10] in EXCHANGE}
    assert exchange["a7562b8ed1"].get("self_inverse") is True
    assert not exchange["a7562b8ed1"].get("derived")
    assert exchange["3e484af0e1"]["derived"]["primary"]["key"][:10] == "a7562b8ed1"


def _scaled(payloads, key, factor):
    out = json.loads(json.dumps(payloads))
    for p in out:
        for t in p["templates"]:
            if t.get("key", "").startswith(key):
                t["weight"]["w"] = [w * factor for w in t["weight"]["w"]]
    return out


def _refused_by_overlap(out):
    return {t["key"][:10] for p in out for t in p["templates"]
            if (t.get("refusal") or {}).get("status") == "refused:class_overlap"}


def test_overlapping_templates_whose_weights_differ_are_refused_by_name():
    """Two templates making one transition at weights 0.021 decades apart (>
    the locality tolerance): no choice of one is right, the component is
    refused refused:class_overlap with the numbers."""
    out, table = te.resolve_overlaps(_scaled(collided(), "247abac61b", 1.05))
    assert [r["disposition"] for r in table] == ["refused:class_overlap"]
    assert "decades" in table[0]["reason"]
    assert _refused_by_overlap(out) == set(EXCHANGE)
    assert not any(t["status"] in ("emitted", "merged", "covered") for p in out
                   for t in p["templates"] if t["key"][:10] in EXCHANGE)


def test_an_overlap_inside_one_orientation_is_refused_by_name():
    """With every record as_generated, both partitions are on one side: no
    side may be dropped, the component is refused."""
    payloads = collided()
    for p in payloads:
        p["orientation"] = "as_generated"
    out, table = te.resolve_overlaps(payloads)
    assert [r["disposition"] for r in table] == ["refused:class_overlap"]
    assert _refused_by_overlap(out) == set(EXCHANGE)


def test_a_contained_as_generated_partition_is_covered_by_the_reversed_one():
    """Without a756, the reversed side {53b8, 247a} makes transitions no
    as_generated template makes, and every transition of the as_generated side
    {3e48} is made by a reversed template: the contained side is covered,
    whichever side generated it (I-249), and no transition is lost.  (Before
    I-249 only the reversed side could be covered, and this component was
    refused.)"""
    payloads = collided()
    for p in payloads:
        for t in p["templates"]:
            if t.get("key", "").startswith("a7562b8ed1"):
                t["status"] = "refused"
                t["refusal"] = {"status": "refused:test", "reason": "removed by the test"}
                t.pop("merged_into", None)
    out, table = te.resolve_overlaps(payloads)
    assert [r["disposition"] for r in table] == ["covered"]
    assert keys_of(table[0]["covered"]) == ["3e484af0e1"]
    assert keys_of(table[0]["covered_by"]) == ["247abac61b", "53b80ec7cc"]
    assert _refused_by_overlap(out) == set()
    emitted = {(p["record_event_id"], t["key"]) for p in out for t in p["templates"]
               if t["status"] == "emitted"}
    for p in out:
        for t in p["templates"]:
            if t["status"] == "covered":
                assert {(c["record"], c["key"]) for c in t["covered_by"]} <= emitted


def digest(value):
    """A cheap comparison of large JSON (a failing ``==`` on megabytes of
    nested data makes pytest's report, not the test, run for minutes)."""
    import hashlib
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def _relabelled(payloads, mapping):
    text = json.dumps(payloads, sort_keys=True)
    for old, new in mapping.items():
        text = text.replace(old, new)
    return json.loads(text)


def test_overlap_resolution_does_not_follow_the_record_ids():
    payloads = collided()
    ids = sorted(p["record_event_id"] for p in payloads)
    tmp = {old: "evt_" + f"{k:064x}" for k, old in enumerate(ids)}
    swap = dict(zip(ids, reversed(ids)))
    there = _relabelled(_relabelled(payloads, tmp), {tmp[k]: v for k, v in swap.items()})
    out_here, table_here = te.resolve_overlaps(payloads)
    out_there, table_there = te.resolve_overlaps(there)
    back = {v: k for k, v in swap.items()}
    tmp_back = {old: "evt_" + f"{k + 100:064x}" for k, old in enumerate(ids)}
    again = _relabelled(_relabelled(out_there, tmp_back),
                        {tmp_back[k]: v for k, v in back.items()})
    assert digest([{k: v for k, v in p.items() if k != "sha256"} for p in again]) == \
        digest([{k: v for k, v in p.items() if k != "sha256"} for p in out_here])
    assert [te.partition_free_identity(p) for p in out_there] == \
        [te.partition_free_identity(p) for p in out_here]
    assert digest(table_here) == digest(
        _relabelled(_relabelled(table_there, tmp_back), {tmp_back[k]: v for k, v in back.items()}))
    covered = [(p, t) for p in out_here for t in p["templates"] if t["status"] == "covered"]
    assert covered
    for p, t in covered:
        q = next(x for x in again if x["record_event_id"] == p["record_event_id"])
        u = next(x for x in q["templates"] if x["key"] == t["key"])
        assert te.partition_free_identity(p, t) == te.partition_free_identity(q, u)


# --------------------------------------------------------------------------
# a covered template in subsets and catalog sets
# --------------------------------------------------------------------------
def test_a_subset_with_a_covered_template_carries_its_coverers():
    out = te.resolve_overlaps(collided())[0]
    reversed_only = [p for p in out if p["orientation"] == "reversed"]
    with pytest.raises(te.EmissionError, match="covered"):
        te.loaded_templates(reversed_only, 700.0)
    carried = te.carry_owners(reversed_only, out)
    names = {t.name.split(":")[1][:10] for t in te.loaded_templates(carried, 700.0)}
    assert {"3e484af0e1", "a7562b8ed1"} <= names
    assert not names & {"53b80ec7cc", "247abac61b"}


def test_catalog_set_links_name_the_coverers():
    from rmgpy.kmc import catalog_set
    out = te.resolve_overlaps(collided())[0]
    links = [x for x in catalog_set.payload_links(out) if x[0] == "covered_by"]
    covered = sum(len(t["covered_by"]) for p in out for t in p["templates"]
                  if t["status"] == "covered")
    assert covered and len(links) == covered
    assert {x[3][:10] for x in links} == {"3e484af0e1", "a7562b8ed1"}
