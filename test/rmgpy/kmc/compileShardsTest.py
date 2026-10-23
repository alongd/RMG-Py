"""The catalog compile split into shards: census, recount, build, merge.

The merged artifact must be the single-process artifact byte for byte, for any
shard counts; a missing, failed, incomplete or foreign shard must stop the
next stage, which names it."""

import gzip
import pickle
import os
from pathlib import Path

import pytest

from rmgpy.kmc import compile as compile_module
from rmgpy.kmc.compile_shards import (
    ShardCounts,
    ShardInputError,
    blocks,
    intercept_recounts,
    read_stage,
    round_robin,
    run_merge,
    run_sharded_in_process,
    serving,
    stage_path,
    write_stage,
)

REPOSITORY = Path(__file__).resolve().parents[3]
DECKS = REPOSITORY / "test/rmgpy/kmc/fixtures/decks"
DATABASE = Path(os.environ.get("RMG_DATABASE_PATH", "/nonexistent"))
needs_database = pytest.mark.skipif(
    not (DATABASE / "input/kinetics/families").is_dir(),
    reason="RMG_DATABASE_PATH is not a database checkout",
)
FAMILIES = ("R_Recombination", "H_Abstraction", "R_Addition_MultipleBond",
            "intra_H_migration", "Disproportionation")


# --- partitions ----------------------------------------------------------------


@pytest.mark.parametrize("count", [1, 2, 7])
@pytest.mark.parametrize("n", [0, 1, 5, 39, 100])
def test_blocks_cover_every_index_once_in_order_and_balance(n, count):
    parts = blocks(n, count)
    assert len(parts) == count
    assert [i for start, stop in parts for i in range(start, stop)] == list(range(n))
    sizes = [stop - start for start, stop in parts]
    assert max(sizes) - min(sizes) <= 1


@pytest.mark.parametrize("count", [1, 2, 7])
def test_round_robin_assigns_every_index_to_exactly_one_shard(count):
    owned = [round_robin(23, shard, count) for shard in range(count)]
    assert sorted(i for indices in owned for i in indices) == list(range(23))
    assert all(indices == sorted(indices) for indices in owned)


# --- the recount interception ---------------------------------------------------


class _Family:
    """calculate_degeneracy of 'outer' itself recounts 'inner' (a nested call)."""

    def __init__(self):
        self.calls = []

    def calculate_degeneracy(self, reaction, resonance=True):
        self.calls.append(reaction)
        if reaction == "outer":
            return self.calculate_degeneracy("inner", resonance=resonance) + 1.0
        return 2.0


def test_interception_sees_top_level_recounts_by_ordinal_and_restores_the_family():
    family = _Family()
    seen = []

    def handler(ordinal, reaction, compute):
        seen.append((ordinal, reaction))
        return compute() if reaction == "outer" else -1.0

    with intercept_recounts(family, handler):
        assert family.calculate_degeneracy("outer") == 3.0
        assert family.calculate_degeneracy("skipped") == -1.0
    # the nested 'inner' recount ran inside the original, unseen by the handler
    assert seen == [(0, "outer"), (1, "skipped")]
    assert family.calls == ["outer", "inner"]
    assert "calculate_degeneracy" not in vars(family)
    assert family.calculate_degeneracy("x") == 2.0


def test_a_served_recount_is_checked_against_the_census_key_and_refuses_drift():
    family = _Family()
    handler = serving("unit 3 (site+site, Fam)", ["k0", "k1"], {0: 4.0, 1: 6.0}, key=str)
    with intercept_recounts(family, handler):
        with pytest.raises(ShardInputError, match=r"unit 3 \(site\+site, Fam\).*recount 0"):
            family.calculate_degeneracy("not-k0")
    with intercept_recounts(family, serving("u", ["k0"], {0: 4.0}, key=str)):
        assert family.calculate_degeneracy("k0") == 4.0
        with pytest.raises(ShardInputError, match="u.*recount 1.*census recorded 1"):
            family.calculate_degeneracy("k1")
    assert family.calls == []


# --- stage files: a missing, failed, incomplete or foreign shard is named ------------


HEADER = {"request": {"families": ["F"]}, "code": "abc"}


def test_read_stage_returns_the_payloads_in_shard_order(tmp_path):
    for shard in (1, 0, 2):
        write_stage(tmp_path, "build", shard, 3, HEADER, {"shard": shard})
    assert [p["shard"] for p in read_stage(tmp_path, "build", 3, HEADER)] == [0, 1, 2]


def test_a_missing_shard_is_named(tmp_path):
    write_stage(tmp_path, "build", 0, 3, HEADER, {})
    write_stage(tmp_path, "build", 2, 3, HEADER, {})
    with pytest.raises(ShardInputError, match="build shard 1 of 3.*missing"):
        read_stage(tmp_path, "build", 3, HEADER)


def test_a_shard_that_failed_while_writing_leaves_nothing_and_is_named(tmp_path):
    write_stage(tmp_path, "recount", 0, 2, HEADER, {})
    with pytest.raises(Exception):
        write_stage(tmp_path, "recount", 1, 2, HEADER, {"value": lambda: None})
    assert sorted(p.name for p in (tmp_path / "recount").iterdir()) == [
        stage_path(tmp_path, "recount", 0, 2).name]
    with pytest.raises(ShardInputError, match="recount shard 1 of 2.*missing"):
        read_stage(tmp_path, "recount", 2, HEADER)


def test_a_truncated_shard_is_named(tmp_path):
    for shard in range(2):
        write_stage(tmp_path, "census", shard, 2, HEADER, {"units": list(range(50))})
    path = stage_path(tmp_path, "census", 1, 2)
    path.write_bytes(path.read_bytes()[:-20])
    with pytest.raises(ShardInputError, match="census shard 1 of 2.*unreadable"):
        read_stage(tmp_path, "census", 2, HEADER)


def test_a_shard_marked_incomplete_is_named(tmp_path):
    write_stage(tmp_path, "census", 0, 1, HEADER, {})
    path = stage_path(tmp_path, "census", 0, 1)
    document = pickle.loads(gzip.decompress(path.read_bytes()))
    document["complete"] = False
    path.write_bytes(gzip.compress(pickle.dumps(document)))
    with pytest.raises(ShardInputError, match="census shard 0 of 1.*not complete"):
        read_stage(tmp_path, "census", 1, HEADER)


@pytest.mark.parametrize("field, value", [("request", {"families": ["G"]}), ("code", "abd")])
def test_a_shard_of_another_request_or_code_is_named(tmp_path, field, value):
    write_stage(tmp_path, "build", 0, 2, HEADER, {})
    write_stage(tmp_path, "build", 1, 2, {**HEADER, field: value}, {})
    with pytest.raises(ShardInputError, match=f"build shard 1 of 2.*{field}"):
        read_stage(tmp_path, "build", 2, HEADER)


def test_a_shard_written_for_another_shard_count_is_named(tmp_path):
    write_stage(tmp_path, "build", 0, 2, HEADER, {})
    write_stage(tmp_path, "build", 1, 2, HEADER, {})
    path = stage_path(tmp_path, "build", 1, 2)
    document = pickle.loads(gzip.decompress(path.read_bytes()))
    document["of"] = 3
    path.write_bytes(gzip.compress(pickle.dumps(document)))
    with pytest.raises(ShardInputError, match="build shard 1 of 2.*3 shards"):
        read_stage(tmp_path, "build", 2, HEADER)


# --- real compiles: the merged artifact is the single-process artifact -------------


# deck, proxy units, families, sites (a reduced set where the full one is too slow here)
CASES = {
    "pe": ("pe_toy_deck.py", None, FAMILIES, None),
    "pp": ("pp_deck.py", None, FAMILIES,
           ("interior_radical:1", "pristine", "tail_end_radical+pendant_radical:0",
            "head_end_radical_short+monomer")),
    # reverse-generated Disproportionation: every reaction's degeneracy is recounted
    "ps": ("ps_deck.py", 2, FAMILIES,
           ("pendant_radical:2_short+monomer", "interior_radical:0_short+monomer")),
}


def _subset(monkeypatch, sites):
    if sites is None:
        return
    original = compile_module._proxy_inputs

    def filtered(request):
        proxies, source, options, discovery = original(request)
        return [p for p in proxies if p.site_type in sites], source, options, discovery

    monkeypatch.setattr(compile_module, "_proxy_inputs", filtered)


def _request(case, output):
    deck, units, families, _ = CASES[case]
    return compile_module.CompileRequest(
        database_path=str(DATABASE), output_dir=str(output), families=families,
        proxy_units=units, deck=str(DECKS / deck), database_sha="sharding-identity")


@needs_database
@pytest.mark.parametrize("case", sorted(CASES))
def test_the_merged_artifact_is_the_single_process_artifact(case, tmp_path, monkeypatch):
    _subset(monkeypatch, CASES[case][3])
    single, single_id = compile_module.compile_artifact(_request(case, tmp_path / "single"))
    database = None
    for count in (1, 2, 7):
        work = tmp_path / f"work{count}"
        merged, merged_id, database = run_sharded_in_process(
            _request(case, tmp_path / f"merged{count}"), work,
            ShardCounts(census=count, recount=count, build=count), database=database)
        assert merged_id == single_id, f"K={count}"
        assert merged.read_bytes() == single.read_bytes(), f"K={count}"
        census = read_stage(work, "census", count, None)
        recounts = sum(len(unit["recount_keys"]) for part in census for unit in part["units"])
        if case == "ps":
            assert recounts > 0
    # the merge refuses partial input and names the shard
    stage_path(tmp_path / "work2", "build", 1, 2).unlink()
    with pytest.raises(ShardInputError, match="build shard 1 of 2.*missing"):
        run_merge(_request(case, tmp_path / "partial"), tmp_path / "work2",
                  ShardCounts(census=2, recount=2, build=2), database=database)
    assert not (tmp_path / "partial").exists() or not list((tmp_path / "partial").iterdir())
