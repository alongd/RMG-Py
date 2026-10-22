"""One status per template across the records that emit it (I-237).

A template is cut from one record's witness, but its key (its executable
content) is shared by every record that cuts the same window, and the
catalog loads it by key.  The checks after the cut (locality, orbit,
spread) run on each record's own witnesses, so two records can disagree:
on the full PE catalog, Disproportionation template e06bf2aa493546b3 was
refused ``refused:non_transferable`` on one record and emitted on three, and
its reverse likewise.  Assembly kept the status of the LAST copy it read, so
the catalog depended on the shard order (a reciprocal split in some orders,
a template both loaded and refused in others).

The rule: a template refused at emission on any record's witnesses is
refused on every record, and so is the reciprocal partner it is decided
with; the refusal names the records whose witnesses refused it.
"""

import gzip
import itertools
import json
from pathlib import Path

import pytest

from rmgpy.kmc import template_emission as te

import rmg_moves_harness as harness
import templateEmissionTest as tet

ROWS = Path(__file__).resolve().parent / "fixtures" / "reciprocal_split" / \
    "pe_disproportionation_rows.jsonl.gz"
FORWARD, REVERSE = "e06bf2aa493546b3", "a859581e53b180dd"
# per key, the record whose own witnesses refuse it (the two are each other's reverse)
REFUSING = {FORWARD: {"evt_56574907b51fda5ab0a4"}, REVERSE: {"evt_83be647047c53854d9a8"}}


def _rows():
    with gzip.open(ROWS, "rt") as handle:
        return [json.loads(line) for line in handle]


def _statuses(payloads):
    out = {}
    for p in payloads:
        for t in p["templates"]:
            out.setdefault(t["key"], set()).add("emitted" if t["status"] == "merged"
                                                else t["status"])
    return out


def _refuse(payload, reason="port 0:tail +1: not one orbit"):
    p = json.loads(json.dumps(payload))
    t = tet.base_template(p)
    t["status"] = "refused"
    t["refusal"] = {"status": "refused:non_transferable", "reason": reason}
    p["status"], p["refusal"] = te.record_status(p["templates"])
    return te.finish(p)


def test_a_template_refused_on_one_record_is_refused_on_every_record_in_any_order():
    """Three forward records share a key, one refuses it on its own
    witnesses; their reverses share another key, all emitted.  In every
    input order both keys are refused on every record, naming the refusing
    record, and the pair is never split."""
    kf, kr = "a" * 64, "b" * 64
    payloads = tet._reciprocal([kf, kf, kf], [kr, kr, kr])
    payloads[1] = _refuse(payloads[1])
    canon = set()
    for order in itertools.permutations(payloads):
        out, table = te.resolve_collisions(list(order))
        assert _statuses(out) == {kf: {"refused"}, kr: {"refused"}}
        assert not [r for r in table if r["disposition"] != "merged"]
        for p in out:
            t = tet.base_template(p)
            assert t["refusal"]["status"] == "refused:non_transferable"
            assert t["refusal"]["records"] == [tet._rid(2)]
            te.validate_payload(p)
        canon.add(tet._canon(out))
    assert len(canon) == 1


def test_the_full_pe_split_rows_give_one_status_per_template_in_any_order():
    """The eight Zeus rows of the PE N=2e4 catalog that split the
    Disproportionation pair: in shard order, reversed and rotated, each of
    the two keys has one status on every record (refused, naming the record
    whose witnesses refuse it); the other templates keep theirs."""
    payloads = [r["payload"] for r in _rows()]
    orders = [payloads, payloads[::-1]] + [payloads[k:] + payloads[:k] for k in range(1, 8)]
    for order in orders:
        out, _table = te.resolve_collisions(list(order))
        by = {k[:16]: v for k, v in _statuses(out).items()}
        assert by[FORWARD] == by[REVERSE] == {"refused"}
        assert all(v == {"emitted"} for k, v in by.items() if k not in (FORWARD, REVERSE))
        for p in out:
            for t in p["templates"]:
                if t["key"][:16] in (FORWARD, REVERSE):
                    assert {r[:24] for r in t["refusal"]["records"]} == REFUSING[t["key"][:16]]
                    assert t["refusal"]["status"] == "refused:non_transferable"
            te.validate_payload(p)


@pytest.mark.skipif(not harness.DATABASE.is_dir(), reason=f"needs the RMG database at {harness.DATABASE}")
def test_the_refusal_is_the_record_witness_not_the_run():
    """The disagreement is the witnesses, not nondeterminism: on the refusing
    record's port 0:tail +1 witness (hept-3-yl + hept-4-yl) the H transfer in
    either role gives heptane + hept-3-ene, so the channel holds the
    template's transitions of both roles, not one orbit; on an emitting
    record's witness of the same label the channel is one orbit."""
    from rmgpy.kmc import emit
    families = emit.load_database(["Disproportionation"], harness.DATABASE).kinetics.families
    family = families["Disproportionation"]
    counter = te.RmgCounter(families)
    verdict = {}
    for row in _rows():
        if not row["event_id"].startswith(("evt_56574907", "evt_20bdc63f")):
            continue
        payload = row["payload"]
        t = next(t for t in payload["templates"] if t["key"].startswith(FORWARD))
        w = t["maps"]["witnesses"]["port 0:tail +1"]
        r, p = te.rows_from_json(w["r"]), te.rows_from_json(w["p"])
        count, inside, _channel = te.count_on(
            [te.matcher_from_payload(t, 700.0, weight=1.0)], r, p, family=family,
            direction=payload["direction"], counter=counter,
            pathway=te.count_pathway(family, payload["pathway"], te.template_on_products(payload)))
        verdict[row["event_id"][:12]] = (count.matcher, count.equal,
                                         te.orbit_problem(r, [tr.edits for tr, _k in inside]) is None)
    assert verdict == {"evt_56574907": (6, True, False), "evt_20bdc63f": (2, True, True)}
