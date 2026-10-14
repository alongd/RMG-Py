"""One inventory choice (R0 or R1) must drive the whole SSA engine, both ways."""

import hashlib
import json
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import pytest

from rmgpy.kmc.met import (
    TRANSPORT_ARMS,
    CompiledChannel,
    CompiledTerminationTable,
    RateTable,
)
from rmgpy.kmc.ssa import IsothermalSSA, SiteIndex, eligible_candidates
from rmgpy.kmc.state import AtomRef, KMCState, Strand
from ssaTest import _leak_records, _one_atom_radical_record

ARM = TRANSPORT_ARMS["A0_REF_H_CROSS_NEc"]
J_RING = "R1:J_ring"
QUINOID = "R1:quinoid_disproportionation"
DP4_ARTIFACT = Path(
    "/home/alon/runs/i127/probe/fullcompile/dp4_closed+pe/"
    "8c7a95f3b46af9f4c402f67940f10d658ac9897b039b704168dfef8ae1f56d8d.json"
)


def _radical_strand(strand_id, atom_uuid, length):
    ref = AtomRef(atom_uuid, strand_id, 0, "chain_end")
    return Strand(
        strand_id,
        length,
        formula={"C": 1},
        radical_count=1,
        atom_refs={0: ref},
        atom_graph={
            atom_uuid: {
                "element": "C",
                "radical": 1,
                "charge": 0,
                "lone_pairs": 0,
                "implicit_hydrogens": 0,
                "edges": {},
            }
        },
    )


def _state():
    return KMCState((_radical_strand("left", "left-C", 5), _radical_strand("right", "right-C", 5)))


def _as_dict(record, event_id, inventory_class, reverse_of):
    data = asdict(record)
    data.pop("_KNOWN_FIELDS", None)
    data.update(
        event_id=event_id, inventory_class=inventory_class, reverse_of=reverse_of
    )
    return data


def _termination(event_id, inventory_class=None, reverse_of=None):
    base = replace(
        _one_atom_radical_record(),
        event_id="",
        family="R_Recombination",
        arity=2,
        reactant_multiplicities=[2],
        radical_delta=-2,
        rate_order=2,
        rate_units="m^3 mol^-1 s^-1",
        status="irreversible",
        status_reason="toy termination",
        k_table={
            "T": [600.0, 800.0],
            "k": [1e12, 1e12],
            "interpolation": "linear-ln-k",
            "extrapolation": "refuse",
        },
        reactant_graphs=["1 *1 C u1 p0 c0", "1 *1 C u1 p0 c0"],
        product_graphs=["1 *1 C u0 p0 c0 {2,S}\n2 *2 C u0 p0 c0 {1,S}"],
        bond_ops=[
            {"action": "set_radical", "atom": 0, "value": 0},
            {"action": "set_radical", "atom": 1, "value": 0},
            {"action": "form", "atoms": [0, 1], "order": "S"},
        ],
        reactant_pair_convention="unordered-identical-pair N(N-1)/2",
    )
    return _as_dict(base, event_id, inventory_class, reverse_of)


def _reverse(event_id, inventory_class=None, reverse_of=None):
    return _as_dict(_one_atom_radical_record(), event_id, inventory_class, reverse_of)


def _records(forward_tag=J_RING, reverse_tag=J_RING):
    """A kept termination K-f plus a reciprocal pair J-f / J-r."""
    return (
        _termination("K-f"),
        _termination("J-f", forward_tag, "J-r"),
        _reverse("J-r", reverse_tag, "J-f"),
    )


def _channel(record, inventory_class="base"):
    return CompiledChannel(
        record["event_id"],
        record["event_id"],
        record["family"],
        inventory_class,
        RateTable((600.0, 800.0), (1e12, 1e12)),
        1.0,
        False,
        None,
        None,
        None,
        None,
        (),
    )


def _table(inventory, records, ids):
    by_id = {record["event_id"]: record for record in records}
    channels = tuple(
        _channel(by_id[key], J_RING if key == "J-f" else "base") for key in ids
    )
    return CompiledTerminationTable("bulk", ARM, inventory, channels)


def _engine(
    records,
    *,
    inventory,
    table=None,
    incremental=True,
    state=None,
    seed=1,
    **extra,
):
    return IsothermalSSA(
        state if state is not None else _state(),
        records,
        temperature=700.0,
        volume=1.0,
        rng=np.random.default_rng(seed),
        termination_table=table,
        incremental=incremental,
        inventory=inventory,
        **extra,
    )


def _retained(engine):
    return {record["event_id"] for record in engine.records}


MODES = pytest.mark.parametrize("incremental", (True, False), ids=("incr", "full"))
TABLES = pytest.mark.parametrize("with_table", (True, False), ids=("table", "notable"))


def _setup(inventory, records, with_table):
    if not with_table:
        return None
    ids = ("K-f",) if inventory == "R0" else ("K-f", "J-f")
    return _table(inventory, records, ids)


def test_fixture_j_ring_forward_has_eligible_sites():
    records = _records()
    state = _state()
    index = SiteIndex(state, records)
    forward = next(record for record in records if record["event_id"] == "J-f")
    assert len(eligible_candidates(state, index, forward)) > 0
    reverse = next(record for record in records if record["event_id"] == "J-r")
    assert len(eligible_candidates(state, index, reverse)) > 0


@MODES
@TABLES
def test_r1_positive_control_keeps_j_ring_pair_with_nonzero_propensity(
    incremental, with_table
):
    records = _records()
    engine = _engine(
        records,
        inventory="R1",
        table=_setup("R1", records, with_table),
        incremental=incremental,
    )
    assert _retained(engine) == {"K-f", "J-f", "J-r"}
    report = engine.channel_propensities()
    assert report.propensities["J-f"] > 0.0
    assert report.propensities["J-r"] > 0.0


@MODES
@TABLES
def test_r0_removes_j_ring_pair_everywhere(incremental, with_table):
    records = _records()
    state = _state()
    # Eligibility is explicit, so zero propensity below is not vacuous.
    forward = records[1]
    assert len(eligible_candidates(state, SiteIndex(state, records), forward)) > 0
    engine = _engine(
        records,
        inventory="R0",
        table=_setup("R0", records, with_table),
        incremental=incremental,
        state=state,
    )
    report = engine.channel_propensities()
    assert report.propensities.get("J-f", 0.0) == 0.0
    assert _retained(engine) == {"K-f"}
    assert "J-f" not in report.propensities
    assert "J-r" not in report.propensities
    assert report.propensities["K-f"] > 0.0
    assert report.enabled_channel_ids == frozenset({"K-f"})
    assert engine.excluded_event_ids == frozenset({"J-f", "J-r"})


@MODES
@TABLES
@pytest.mark.parametrize(
    ("forward_tag", "reverse_tag"),
    ((J_RING, None), (None, J_RING), (QUINOID, None), (None, QUINOID)),
    ids=("fwd-ring", "rev-ring", "fwd-quinoid", "rev-quinoid"),
)
def test_reciprocal_closure_removes_untagged_partner(
    incremental, with_table, forward_tag, reverse_tag
):
    records = _records(forward_tag, reverse_tag)
    engine = _engine(
        records,
        inventory="R0",
        table=_setup("R0", records, with_table),
        incremental=incremental,
    )
    assert _retained(engine) == {"K-f"}
    report = engine.channel_propensities()
    assert "J-f" not in report.propensities
    assert "J-r" not in report.propensities


def test_untagged_pair_survives_r0():
    records = _records(None, None)
    engine = _engine(records, inventory="R0")
    assert _retained(engine) == {"K-f", "J-f", "J-r"}


def test_r0_met_layout_and_ordinary_pool_exclude_removed_records():
    records = _records()
    engine = _engine(
        records, inventory="R0", table=_setup("R0", records, True), incremental=True
    )
    maintenance = engine._incremental_propensities()
    maintenance.sync()
    removed = {"J-f", "J-r"}
    assert not removed & set(maintenance.layout.records)
    assert not removed & set(maintenance.layout.record_order)
    assert not removed & set(maintenance.records)
    assert not removed & set(engine._records_by_id)
    assert not removed & set(maintenance.enabled)
    # K-f is owned by its table channel and never by the ordinary pool.
    assert "K-f" in maintenance.layout.met_event_ids
    assert "K-f" not in maintenance.records


def test_partition_each_termination_record_has_exactly_one_owner():
    records = _records(None, None) + (_termination("K2-f"),)
    table = _table("R1", records, ("K-f",))
    engine = _engine(records, inventory="R1", table=table)
    maintenance = engine._incremental_propensities()
    maintenance.sync()
    owned_by_table = set(maintenance.layout.met_event_ids)
    pool = set(maintenance.records)
    assert not owned_by_table & pool
    assert owned_by_table | pool == _retained(engine)
    assert "K-f" in owned_by_table and "K2-f" in pool


def test_stale_table_row_for_removed_record_is_refused():
    records = _records()
    stale = _table("R0", records, ("K-f", "J-f"))
    with pytest.raises(ValueError, match="removed"):
        _engine(records, inventory="R0", table=stale)


def test_table_row_for_absent_record_is_refused():
    records = _records()
    ghost = _termination("ghost-f")
    table = CompiledTerminationTable("bulk", ARM, "R0", (_channel(ghost),))
    with pytest.raises(ValueError, match="no record|absent|missing"):
        _engine(records, inventory="R0", table=table)


@pytest.mark.parametrize("engine_inventory", ("R0", "R1"))
def test_inventory_mismatch_between_engine_and_table_is_refused(engine_inventory):
    records = _records(None, None)
    other = "R1" if engine_inventory == "R0" else "R0"
    table = _table(other, records, ("K-f",))
    with pytest.raises(ValueError, match="inventory"):
        _engine(records, inventory=engine_inventory, table=table)


def test_inventory_choice_is_required_and_validated():
    records = _records(None, None)
    with pytest.raises(TypeError):
        IsothermalSSA(
            _state(),
            records,
            temperature=700.0,
            volume=1.0,
            rng=np.random.default_rng(1),
        )
    with pytest.raises(ValueError, match="inventory"):
        _engine(records, inventory="R2")


@pytest.mark.parametrize("incremental", (True, False), ids=("incr", "full"))
def test_dangling_reverse_of_to_missing_record_is_refused(incremental):
    records = (_termination("K-f", None, "nowhere"), *_records()[1:])
    with pytest.raises(ValueError, match="reverse_of"):
        _engine(records, inventory="R0", incremental=incremental)


def test_dangling_reverse_of_to_removed_record_is_refused_not_cleared():
    # K-f points one-way at the removed J-f; J-f/J-r are the reciprocal pair.
    records = (_termination("K-f", None, "J-f"), *_records()[1:])
    original = {record["event_id"]: record["reverse_of"] for record in records}
    with pytest.raises(ValueError, match="reverse_of"):
        _engine(records, inventory="R0")
    assert {record["event_id"]: record["reverse_of"] for record in records} == original


def test_non_reciprocal_tag_does_not_pull_in_a_third_record():
    # A one-way link must be refused, never silently widen the removal.
    records = (_termination("K-f", None, "J-f"), *_records()[1:])
    with pytest.raises(ValueError):
        _engine(records, inventory="R0")


def test_absent_partner_is_refused_by_default_and_admitted_with_partial_flag():
    records = (_termination("K-f", None, "outside"), *_records(None, None)[1:])
    with pytest.raises(ValueError, match="missing"):
        _engine(records, inventory="R0")
    engine = _engine(records, inventory="R0", partial_record_set=True)
    assert engine.external_reverse_ids == frozenset({"outside"})
    assert "K-f" in _retained(engine)
    assert _engine(_records(None, None), inventory="R0").external_reverse_ids == frozenset()


def test_removed_partner_is_refused_even_with_partial_flag():
    records = (_termination("K-f", None, "J-f"), *_records()[1:])
    with pytest.raises(ValueError, match="removed"):
        _engine(records, inventory="R0", partial_record_set=True)


def test_excluded_record_with_absent_partner_is_removed_under_partial_flag():
    records = (_termination("K-f"), _termination("J-f", J_RING, "outside"))
    engine = _engine(records, inventory="R0", partial_record_set=True)
    assert _retained(engine) == {"K-f"}
    assert engine.external_reverse_ids == frozenset()


def _digest(result):
    return hashlib.sha256(
        json.dumps(
            [(repr(e.time), e.event_id, e.state_hash, e.irreversible) for e in result.events]
        ).encode()
    ).hexdigest()


def _leak_trace(inventory):
    state = KMCState((_radical_strand("strand", "strand-C", 5),))
    engine = IsothermalSSA(
        state,
        _leak_records(1.0, 0.2),
        temperature=700.0,
        volume=1.0,
        rng=np.random.default_rng(6_013_700),
        inventory=inventory,
    )
    result = engine.run(until=2.0)
    return result, _digest(result)


# Pinned on the base commit, before inventory filtering existed.
LEAK_TRACE = (
    1,
    "a31300d5317e759c2c1b7981f6a0e1996792616317626ad4d06df1c5345230d0",
    0.33333633032773147,
    0.6666546786890744,
)
MET_TRACE = (
    12,
    "b0512b8e958a2579d13d76e68be9155df59d12700e90b7739988015def783e99",
    0.0,
    12.514585783018111,
)


def test_r1_leak_scenario_without_excluded_classes_is_identical_to_base():
    result, digest = _leak_trace("R1")
    assert (len(result.events), digest, result.leak_num, result.leak_den) == LEAK_TRACE


@MODES
def test_r1_met_table_scenario_without_excluded_classes_is_identical_to_base(
    incremental,
):
    records = (_termination("K-f"), _reverse("A-1"))
    state = KMCState(
        tuple(_radical_strand(f"s{i}", f"s{i}-C", 5) for i in range(12))
    )
    engine = _engine(
        records,
        inventory="R1",
        table=_table("R1", records, ("K-f",)),
        state=state,
        incremental=incremental,
        seed=7,
    )
    result = engine.run(max_events=40)
    assert (
        len(result.events),
        _digest(result),
        result.leak_num,
        result.leak_den,
    ) == MET_TRACE


def _closure(records):
    by_id = {record["event_id"]: record for record in records}
    tagged = {
        key
        for key, record in by_id.items()
        if record.get("inventory_class") in {J_RING, QUINOID}
    }
    closure = set(tagged)
    for key in list(tagged):
        partner = by_id[key].get("reverse_of")
        if partner in by_id and by_id[partner].get("reverse_of") == key:
            closure.add(partner)
    return tagged, closure


@pytest.mark.skipif(not DP4_ARTIFACT.exists(), reason="DP4 artifact absent")
def test_dp4_artifact_removed_count_equals_reciprocal_closure():
    records = tuple(json.loads(DP4_ARTIFACT.read_text())["records"])
    tagged, closure = _closure(records)
    assert len(tagged) == 180 and len(closure) == 360  # partners are untagged
    state = KMCState((_radical_strand("a", "a-C", 5),))
    engine = _engine(records, inventory="R0", state=state)
    removed = {record["event_id"] for record in records} - _retained(engine)
    assert removed == closure
    r1 = _engine(records, inventory="R1", state=state)
    assert len(r1.records) == len(records)
