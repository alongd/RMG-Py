"""A loaded template is derived from its partner wherever its copy is held (I-263).

``resolve_collisions`` makes one copy of a template the loaded one (the
representative); which record holds it is a content order, not a chemical
choice.  ``derive_weights`` then looked for the template's partners only
through that holding record's own ``reverse_of`` record.  When the holder's
reverse record was absent from the catalog (an incomplete shard set), the
template was never derived although its partner was loaded: both directions
stayed independent tables, ln w linear in 1/T, and detailed balance broke
between grid points by the curvature of ln Kc (about 1e-4 at the midpoints).

Since partners are read as one mutual graph (0188f0cdf), a link found from
either side derives the pair, so the fix (``_reverse_side``) carries weight
only when the PARTNER is held the same way too: neither loaded template's own
record has a reverse record that loads.  The fixture is four
Disproportionation payloads of the PP emission at 448dc6247 (the rows of
``/home/alon/runs/pk-i251/zeus-i249pp/shards``, the first fixture's source),
two reciprocal pairs less one record each:

- evt_250e5c346757 (reversed) and its reverse evt_a10b5a9f81d8, both present;
- evt_1ea01095f857 (reversed), which holds the loaded copy of fcb908cd280c
  over evt_250e's; its reverse evt_adefdefc8fa2 is left out;
- evt_7a543fca74b6, which holds the loaded copy of eb0da01b4aeb over
  evt_a10b's; its reverse evt_893cc6f09d1b is left out.

fcb908cd280c and eb0da01b4aeb are partners.  With the fix they are found
through the reverse records of evt_250e and evt_a10b, and fcb908cd280c is
derived from eb0da01b4aeb.  With the fix reverted neither is derived, and the
evt_250e/evt_a10b pair misses detailed balance off the grid by 1.1e-4.

The first fixture (evt_1f45/evt_024f, and evt_2783 without evt_84b2) stopped
covering the case at f431feaf0.  With each template's certificate dropped from
the assembly identity, evt_024f's copy loads, and evt_024f's reverse is in the
catalog.  It had already stopped needing the fix at 0188f0cdf.  So the tests
find the held copies from the ``resolve_collisions`` output, never from record
ids that an identity change can move.
"""

import gzip
import json
import os
from pathlib import Path

import pytest

from rmgpy.kmc import template_emission as te

REPOSITORY = Path(__file__).resolve().parents[3]
FIXTURE = REPOSITORY / "test/rmgpy/kmc/fixtures/i263/held_copy_payloads.json.gz"
DATABASE = Path(os.environ.get("RMG_DATABASE_PATH", "/nonexistent"))
needs_database = pytest.mark.skipif(
    not (DATABASE / "input/kinetics/families").is_dir(),
    reason="RMG_DATABASE_PATH is not a database checkout",
)


@pytest.fixture(scope="module")
def services():
    from rmgpy.kmc import emit
    return emit.load_services(str(DATABASE), ["Disproportionation"])


@pytest.fixture(scope="module")
def payloads():
    with gzip.open(FIXTURE, "rt") as handle:
        return json.load(handle)


def _derived(payloads, services):
    out, _table = te.resolve_collisions(payloads)
    out, _overlaps = te.resolve_overlaps(out)
    return te.derive_weights(out, counter=services.counter)


def _loaded(out):
    return {(p["record_event_id"], t["key"]): t
            for p in out for t in p.get("templates", []) if t.get("status") == "emitted"}


def _held_copy_cases(out):
    """(holder, key, records) of every loaded template whose holder's
    ``reverse_of`` record is not in the catalog while another record holding a
    copy of it (``records``) has a reverse record that is in the catalog and
    loads: the case the I-263 fix exists for, found from the assembly itself
    and not from record ids that an identity change can move."""
    by_record = {p["record_event_id"]: p for p in out}
    cases = []
    for p in out:
        if p.get("reverse_of") in by_record:
            continue
        for t in p.get("templates", []):
            if t.get("status") != "emitted":
                continue
            others = [r for r in (t.get("collision") or {}).get("records", [])
                      if r != p["record_event_id"]
                      and by_record.get(by_record.get(r, {}).get("reverse_of")) is not None
                      and te._loads(by_record[by_record[r]["reverse_of"]])]
            if others:
                cases.append((p["record_event_id"], t["key"], others))
    return cases


def test_the_fixture_is_the_held_copy_case(payloads):
    """Some loaded template is held by a record whose own reverse record is
    not in the catalog, while another copy of it sits on a record whose
    reverse loads.  Without such a template the derivation tests below pass
    whatever ``derive_weights`` does with the holders."""
    out, _table = te.resolve_collisions(payloads)
    assert _held_copy_cases(out)


@needs_database
def test_a_template_held_by_a_record_without_its_reverse_is_derived_from_its_loaded_partner(
        payloads, services):
    """Every held copy is in a derivation: it is derived from a loaded partner,
    or it is the primary another loaded template is derived from.  At least
    one is derived, and its primary is a held copy too.  That primary is the
    partner reached only through the reverse records of the records that
    hold the other copies."""
    held = {(r, k) for r, k, _o in _held_copy_cases(te.resolve_collisions(payloads)[0])}
    loaded = _loaded(_derived(payloads, services))
    primaries = {(t["derived"]["primary"]["record"], t["derived"]["primary"]["key"])
                 for t in loaded.values() if t.get("derived")}
    for ref in sorted(held):
        assert loaded[ref].get("derived") is not None or ref in primaries, ref
    derived = sorted(ref for ref in held if loaded[ref].get("derived"))
    assert derived, held
    for ref in derived:
        primary = loaded[ref]["derived"]["primary"]
        assert (primary["record"], primary["key"]) in held, (ref, primary)


@needs_database
def test_the_pair_keeps_detailed_balance_off_the_grid(payloads, services):
    """Every pair with both records in the catalog checks to 1e-12, at the
    grid, at its midpoints and at its seeded random temperatures.  The
    holders' pairs are not evaluable in this incomplete catalog, and nothing
    is refused."""
    out, table = te.check_cycles(_derived(payloads, services), counter=services.counter,
                                 kc=services.kc, complete=False)
    present = {p["record_event_id"] for p in out}
    holders = {r for r, _k, _o in _held_copy_cases(te.resolve_collisions(payloads)[0])}
    evaluable = [r for r in table if r["forward"] in present and r["reverse"] in present]
    assert evaluable
    for row in evaluable:
        assert row["status"] == "ok", row["reason"]
        assert row["max_offgrid_error"] <= 1e-12 and row["max_grid_error"] <= 1e-12, row
    assert [r["status"] for r in table if {r["forward"], r["reverse"]} & holders] == \
        ["not_evaluable"] * len(holders)
    assert te.cycle_gates(table, out, complete=False) == []
