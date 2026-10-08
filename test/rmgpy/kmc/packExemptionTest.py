import copy

import pytest

import compilerRealTest as real


def _artifact(records):
    return {"records": records}


def _fixture_records():
    return copy.deepcopy(real._pre_change_artifact()["records"])


def _j_para_records(records):
    return [
        record
        for record in records
        if (record.get("junction_ops") or [{}])[0].get("junction_kind") == "J_para"
    ]


def _repinned_records():
    records = _fixture_records()
    for record in _j_para_records(records):
        record["rate_source"]["database_sha"] = next(
            sha
            for sha in real.ARCHIVED_PACK_DATABASE_SHAS
            if sha != "4a12d36fcdc193ede82c8d1ab5c1653495d445bc"
        )
    return records


def test_repinned_sha_is_accepted_and_four_other_pack_records_are_exact():
    records = _repinned_records()
    exemptions = real._pack_exemptions(_artifact(records))
    assert {entry["pack_event_id"] for entry in exemptions.values()} == set(
        real.ARCHIVED_PACK_EVENT_IDS
    )
    assert all(
        entry["reason"] == real.ARCHIVED_PACK_EXEMPTION_REASON
        and "I037_db_repin_equivalence" in entry["reason"]
        for entry in exemptions.values()
    )
    baseline = {
        record["event_id"]: record for record in real._pre_change_artifact()["records"]
    }
    for event_id in real.ARCHIVED_PACK_EVENT_IDS[2:]:
        actual = next(record for record in records if record["event_id"] == event_id)
        assert real._normalized_base_record(actual) == real._normalized_base_record(
            baseline[event_id]
        )


def test_unknown_sha_is_rejected():
    records = _repinned_records()
    _j_para_records(records)[0]["rate_source"]["database_sha"] = "unknown"
    with pytest.raises(AssertionError):
        real._pack_exemptions(_artifact(records))


def _mutate_arrhenius(record):
    record["rate_source"]["arrhenius"]["A_m3_mol_s"] += 1.0


def _mutate_atom_map(record):
    record["atom_map"]["0"] = 999


def _mutate_rewrite(record):
    record["junction_ops"][0]["exact_inverse_rewrite"].pop()


def _mutate_degeneracy(record):
    record["degeneracy"] += 1.0


def _mutate_k_table(record):
    record["k_table"]["k"][0] *= 2.0


@pytest.mark.parametrize(
    "mutation",
    [
        pytest.param(_mutate_arrhenius, id="arrhenius"),
        pytest.param(_mutate_atom_map, id="atom-map"),
        pytest.param(_mutate_rewrite, id="rewrite"),
        pytest.param(_mutate_degeneracy, id="degeneracy"),
        pytest.param(_mutate_k_table, id="k-table"),
    ],
)
def test_repinned_sha_does_not_hide_record_mutations(mutation):
    records = _repinned_records()
    candidate = _j_para_records(records)[0]
    mutation(candidate)
    with pytest.raises(AssertionError):
        real._pack_exemptions(_artifact(records))


@pytest.mark.parametrize(
    "mutation",
    [
        pytest.param(_mutate_arrhenius, id="arrhenius"),
        pytest.param(_mutate_atom_map, id="atom-map"),
        pytest.param(_mutate_rewrite, id="rewrite"),
        pytest.param(_mutate_degeneracy, id="degeneracy"),
        pytest.param(_mutate_k_table, id="k-table"),
    ],
)
def test_assert_pack_exemption_rejects_archived_record_mutations(mutation):
    records = _repinned_records()
    candidate = _j_para_records(records)[0]
    exemption = {
        "reason": real.ARCHIVED_PACK_EXEMPTION_REASON,
        "pack_event_id": candidate["event_id"],
    }
    mutation(candidate)
    with pytest.raises(AssertionError):
        real._assert_pack_exemption(candidate, exemption)
