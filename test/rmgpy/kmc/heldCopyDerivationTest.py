"""A loaded template is derived from its partner wherever its copy is held (I-263).

``resolve_collisions`` makes one copy of a template the loaded one (the
representative); which record holds it is a content order, not a chemical
choice.  ``derive_weights`` then looked for the template's partners only
through that holding record's own ``reverse_of`` record.  When the holder's
reverse record was absent from the catalog (an incomplete shard set), the
template was never derived although its partner was loaded: both directions
stayed independent tables, ln w linear in 1/T, and detailed balance broke
between grid points by the curvature of ln Kc (a few 1e-4 at the midpoints).

The fixture is three Disproportionation payloads of the PP emission at
448dc6247: a reciprocal pair (evt_1f45 forward, evt_024f reverse) and a
third record, evt_2783, that holds the loaded copies of the pair's reverse
templates; evt_2783's own reverse record (evt_84b2) is left out, as the
ff8d76f12 shards left it out.  The pair's numbers at 312.5 K were the gate's:
S_f/S_r = 6.4258425945324543e+45 against ssa ratio * Kc = 6.4280317979011599e+45.
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
FORWARD, REVERSE, HOLDER = "evt_1f45397da116", "evt_024f51f5d491", "evt_2783dc444441"


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
    return {(p["record_event_id"][:16], t["key"][:12]): t
            for p in out for t in p.get("templates", []) if t.get("status") == "emitted"}


@needs_database
def test_the_fixture_is_the_held_copy_case(payloads, services):
    """The pair's reverse templates load from the holder, whose own reverse
    record is not in the catalog; the pair's forward templates load from the
    pair itself."""
    out, _table = te.resolve_collisions(payloads)
    loaded = _loaded(out)
    holder = next(p for p in out if p["record_event_id"].startswith(HOLDER))
    assert {k for r, k in loaded if r == HOLDER} == {"1fc6e3998bf0", "f27b25aca594"}
    assert {k for r, k in loaded if r == FORWARD} == {"3ab8efd07f72", "f33f5cd2eccc"}
    assert not any(p["record_event_id"] == holder["reverse_of"] for p in out)


@needs_database
def test_a_template_held_by_a_record_without_its_reverse_is_derived_from_its_loaded_partner(
        payloads, services):
    loaded = _loaded(_derived(payloads, services))
    for reverse_key, forward_key in (("f27b25aca594", "3ab8efd07f72"),
                                     ("1fc6e3998bf0", "f33f5cd2eccc")):
        derived = loaded[(HOLDER, reverse_key)].get("derived")
        assert derived is not None, reverse_key
        assert (derived["primary"]["record"][:16], derived["primary"]["key"][:12]) == \
            (FORWARD, forward_key)


@needs_database
def test_the_pair_keeps_detailed_balance_off_the_grid(payloads, services):
    """The pair's check (at the grid, its midpoints and its seeded random
    temperatures) holds to 1e-12; the holder's own pair is not evaluable in
    this incomplete catalog, and nothing is refused."""
    out, table = te.check_cycles(_derived(payloads, services), counter=services.counter,
                                 kc=services.kc, complete=False)
    rows = {(r["forward"][:16], r["reverse"][:16]): r for r in table}
    pair = rows[tuple(sorted((FORWARD, REVERSE)))]
    assert pair["status"] == "ok", pair["reason"]
    assert pair["max_offgrid_error"] <= 1e-12 and pair["max_grid_error"] <= 1e-12
    assert [r["status"] for k, r in rows.items() if HOLDER in k] == ["not_evaluable"]
    assert te.cycle_gates(table, out, complete=False) == []
