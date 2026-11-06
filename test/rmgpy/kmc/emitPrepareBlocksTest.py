"""Emit prepare in J blocks gives the serial prepare's shard inputs (I-304b).

``emit prepare-plan`` fixes J contiguous ranges of artifact order,
``prepare-block`` certifies one range in a fresh process and computes its
records' reference-case keys, ``prepare-finish`` proves every block's exact
range from its bytes and runs the gate, the placement, the republication and
the writes over all records.  No group is presorted, so a block's certificate
may list a tree node's atoms in another order than the serial one
(``certificates._same_certificate``).  The criterion (``compare_prepared``):
every manifest field but a run's own timing and gzip bytes equal, every shard
input item equal, certificates equal up to node atom order.

- J = 1, 2, 4, 7 against serial: stub certification refusing chosen records,
  and real certification of the PE toy and PP short fixtures, each block in
  its own process;
- reciprocal pairs split across blocks where the estimate, the reverse or both
  fail: the reverse is certified through its partner read from the artifact,
  and the refusal crosses the block boundary in the finishing step;
- a missing, duplicate, foreign, short, truncated or altered block, and a block
  killed mid-write, each refused by name, with nothing published;
- two attempts of one block: one seal, the other superseded, a third never runs;
- a block that died of a resource error never seals and blocks publication;
- a published generation is never rewritten.
"""
import gzip
import hashlib
import json
import os
import pickle
import signal
import subprocess
import sys
from pathlib import Path

import pytest

from rmgpy.kmc import emit
from rmgpy.kmc import emit_prepare as ep
from emitShardMemoryTest import (DATABASE, StubCertifier, artifact, needs_database,  # noqa: F401
                                 records, roster_and_doc, stubbed)

HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[2]
PP_SHORT = HERE / "fixtures" / "i304b" / "pp_short_artifact.json.gz"
PP_ROSTER = HERE / "fixtures" / "i304b" / "pp_short.roster.json.gz"


def _serial(artifact, roster, tmp_path, shards=2):
    out = tmp_path / "serial"
    assert emit.main(["prepare", str(artifact), str(roster), "--database", str(DATABASE),
                      "--shards", str(shards), "--out", str(out)]) == 0
    return out


def _plan(artifact, roster, run, blocks):
    assert emit.main(["prepare-plan", str(artifact), str(roster), "--database", str(DATABASE),
                      "--blocks", str(blocks), "--run", str(run)]) == 0


def _blocks(artifact, roster, tmp_path, blocks, shards=2, name=None):
    run = tmp_path / (name or f"run{blocks}")
    _plan(artifact, roster, run, blocks)
    for j in range(blocks):
        assert emit.main(["prepare-block", str(run), "--block", str(j)]) == 0
    assert emit.main(["prepare-finish", str(run), "--shards", str(shards)]) == 0
    return run


def _identical(serial, prepared):
    result = ep.compare_prepared(serial, prepared)
    assert result["differences"] == []
    assert result["items"] > 0
    return result


def _items(prepared):
    manifest = json.loads((Path(prepared) / ep.MANIFEST).read_text())
    return [ep._items(Path(prepared) / f) for f in manifest["files"]]


def _ids(records, kind):
    return [r["event_id"] for r in records if r["rate_source"]["kind"] == kind]


@pytest.mark.parametrize("blocks", [1, 2, 4, 7])
def test_blocks_give_the_serial_shard_inputs(tmp_path, artifact, roster_and_doc, records,
                                             stubbed, blocks):
    roster = roster_and_doc[0]
    stubbed.failing.update(_ids(records, "RMG family estimate")[::3])
    serial = _serial(artifact, roster, tmp_path)
    run = _blocks(artifact, roster, tmp_path, blocks)
    result = _identical(serial, run / "current")
    assert result["node_order"] == []   # stub certificates carry no tree nodes
    manifest = json.loads((run / "current" / ep.MANIFEST).read_text())
    assert manifest["republished"] > 0
    assert [s["block"] for s in manifest["blocks"]["seals"]] == list(range(blocks))
    # each block certified exactly its range, in artifact order; serial certified all
    order = [r["event_id"] for r in records]
    seen = [c.seen for c in stubbed.certifiers]
    assert seen[0] == order
    assert sum(seen[1:], []) == order and len(seen) == 1 + blocks


def test_a_family_selection_in_blocks_is_the_serial_one(tmp_path, artifact, roster_and_doc,
                                                         records, stubbed):
    """With ``--family`` a reverse whose partner is outside the selection has
    no rate owner, in blocks as in serial."""
    roster = roster_and_doc[0]
    stubbed.failing.update(_ids(records, "RMG family estimate")[::2])
    flags = ["--family", "R_Recombination", "--family", "intra_H_migration"]
    serial = tmp_path / "serial"
    assert emit.main(["prepare", str(artifact), str(roster), "--database", str(DATABASE),
                      "--shards", "2", "--out", str(serial), *flags]) == 0
    run = tmp_path / "run"
    assert emit.main(["prepare-plan", str(artifact), str(roster), "--database", str(DATABASE),
                      "--blocks", "3", "--run", str(run), *flags]) == 0
    for j in range(3):
        assert emit.main(["prepare-block", str(run), "--block", str(j)]) == 0
    assert emit.main(["prepare-finish", str(run), "--shards", "2"]) == 0
    _identical(serial, run / "current")
    manifest = json.loads((run / "current" / ep.MANIFEST).read_text())
    assert manifest["only"] == ["R_Recombination", "intra_H_migration"]
    assert len(manifest["artifact_ids"]) < len(records)


def _pair_split(records, blocks):
    """(estimate, reverse) of a reciprocal pair whose records fall in two
    different blocks of ``blocks``."""
    order = [r["event_id"] for r in records]
    ranges = ep.block_ranges(len(order), blocks)
    block_of = {i: next(j for j, (lo, hi) in enumerate(ranges) if lo <= k < hi)
                for k, i in enumerate(order)}
    for r in records:
        partner = r.get("reverse_of")
        if r["rate_source"]["kind"] == "reference-thermo reverse" and partner in block_of and \
                block_of[partner] != block_of[r["event_id"]]:
            return partner, r["event_id"]
    raise AssertionError("no pair split across blocks")


@pytest.mark.parametrize("failing", ["estimate", "reverse", "both"])
def test_a_reciprocal_pair_split_across_blocks_refuses_as_serial(
        tmp_path, artifact, roster_and_doc, records, stubbed, failing):
    roster = roster_and_doc[0]
    estimate, reverse = _pair_split(records, 2)
    stubbed.failing.update({"estimate": [estimate], "reverse": [reverse],
                            "both": [estimate, reverse]}[failing])
    serial = _serial(artifact, roster, tmp_path)
    run = _blocks(artifact, roster, tmp_path, 2)
    _identical(serial, run / "current")
    published = {}
    for items in _items(run / "current"):
        for item in items:
            if item["own"]:
                published[item["source"]] = item
    # both records refused and republished under new ids, whichever failed
    new_estimate = published[estimate]["record"]["event_id"]
    new_reverse = published[reverse]["record"]["event_id"]
    assert new_estimate != estimate and new_reverse != reverse
    assert {published[i]["record"]["status"] for i in (estimate, reverse)} == {"refused"}
    if failing == "reverse":
        assert "its linked partner" in published[estimate]["record"]["status_reason"]
    # the reverse's link and its certificate's rate owner follow the partner's new id
    assert published[reverse]["record"]["reverse_of"] == new_estimate
    assert published[reverse]["certificate"]["owner"]["rate_owner_event_id"] == new_estimate
    # the block holding the reverse read its rate owner (the other block's raw record)
    blocks = stubbed.certifiers[1:]
    holder = next(c for c in blocks if reverse in c.seen)
    assert estimate not in holder.seen
    serial_items = {item["source"]: item for items in _items(serial) for item in items
                    if item["own"]}
    assert serial_items[reverse]["certificate"] == published[reverse]["certificate"]


def test_the_key_pass_from_any_start_is_the_whole_pass(monkeypatch):
    """The species cache decides whether a key collision raises, so a block's
    key pass replays the whole pass's cache history from its reset point."""
    from rmgpy.kmc import channels as chm
    monkeypatch.setattr(ep, "SPECIES_RESET", 4)
    universe = [ep.astream.Entry(f"e{k:02d}", "F", 0, 0, ()) for k in range(23)]

    class Memo:
        def __init__(self):
            self.seen = []

        def __call__(self, mol):
            self.seen.append(mol)
            return f"{mol}@{len(self.seen)}"

    def parts(certifier, d, species):
        key = species(d["event_id"])
        if len(species.seen) == 3:
            raise ValueError("a collision the cache saw")
        return ("F", (key,), "t", ()), [], []

    monkeypatch.setattr(chm, "SpeciesKeys", Memo)
    monkeypatch.setattr(ep, "_one", lambda artifact, e: type("R", (), {
        "to_dict": lambda self: {"event_id": e.event_id}})())
    monkeypatch.setattr(emit, "record_parts", parts)
    whole = [(e.event_id, name) for e, name, _m, _o in ep.record_keys(None, universe, None)]
    assert sum(name is None for _i, name in whole) > 1
    for lo in range(len(universe)):
        for hi in range(lo + 1, len(universe) + 1):
            assert [(e.event_id, name) for e, name, _m, _o in
                    ep.record_keys(None, universe, None, lo, hi)] == whole[lo:hi]
    assert [ep.window_start(k) for k in range(9)] == [0, 0, 0, 3, 3, 3, 3, 7, 7]


# --- refusals ---------------------------------------------------------------
def _sealed_run(artifact, roster, tmp_path, blocks=2):
    run = tmp_path / "run"
    _plan(artifact, roster, run, blocks)
    for j in range(blocks):
        assert ep.prepare_block(run, j)["state"] == "sealed"
    return run


def _finish_refused(run, capsys, *needles):
    capsys.readouterr()
    assert emit.main(["prepare-finish", str(run), "--shards", "2"]) == 1
    out = capsys.readouterr().out
    assert "refused:blocks" in out
    for needle in needles:
        assert needle in out, out
    assert not (run / "current").exists()
    assert not (run / "prepared").exists() or not any((run / "prepared").iterdir())
    return out


def _seal(run, j, blocks=2):
    path = run / "seals" / f"{ep.block_name(j, blocks)}.json"
    return path, json.loads(path.read_text())


def _rewrite_meta(run, j, change, blocks=2):
    """Rewrite block j's meta stream through ``change`` (items -> items) and
    re-seal it: a writer that sealed wrong bytes."""
    path, seal = _seal(run, j, blocks)
    meta = run / "attempts" / seal["meta"]["file"]
    items = [i for i in ep._stream(meta) if i is not None]
    os.chmod(meta, 0o644)   # sealed attempt files are read-only
    with gzip.open(meta, "wb") as handle:
        for item in change(items):
            pickle.dump(item, handle)
    seal["meta"]["sha256"] = ep.file_sha256(meta)
    path.write_bytes(ep.canonical_bytes(seal))


def test_a_missing_block_is_refused_by_name(tmp_path, artifact, roster_and_doc, stubbed,
                                            capsys):
    run = _sealed_run(artifact, roster_and_doc[0], tmp_path)
    _seal(run, 1)[0].unlink()
    _finish_refused(run, capsys, "block b1of2 is missing: no sealed attempt",
                    "sealed, but its seal is gone")


def test_a_foreign_block_is_refused_by_name(tmp_path, artifact, roster_and_doc, records,
                                            stubbed, capsys):
    run = _sealed_run(artifact, roster_and_doc[0], tmp_path)
    path, seal = _seal(run, 0)
    # a seal of another plan's run copied in
    path.write_bytes(ep.canonical_bytes(dict(seal, plan_sha256="0" * 64)))
    _finish_refused(run, capsys, "block b0of2: foreign seal")
    path.write_bytes(ep.canonical_bytes(seal))
    # a seal of another block count
    (run / "seals" / "b2of3.json").write_bytes(ep.canonical_bytes(seal))
    _finish_refused(run, capsys, "foreign seal b2of3.json")
    (run / "seals" / "b2of3.json").unlink()
    # a block holding a record no block of the plan holds
    _rewrite_meta(run, 0, lambda items: items[:2] + [dict(items[2], source="evt_" + "f" * 60)]
                  + items[3:])
    _finish_refused(run, capsys, "block b0of2: foreign: holds evt_ffff", "short: 1 of its")


def test_a_duplicate_block_is_refused_by_name(tmp_path, artifact, roster_and_doc, records,
                                              stubbed, capsys):
    run = _sealed_run(artifact, roster_and_doc[0], tmp_path)
    other = records[-1]["event_id"]   # block 1's last record
    _rewrite_meta(run, 0, lambda items: items[:-2] + [dict(items[1], source=other)]
                  + items[-2:])
    _finish_refused(run, capsys, f"block b0of2: duplicate: holds {other}, a record of block "
                                 "b1of2")
    # a seal of block 1 naming block 0's files
    run2 = _sealed_run(artifact, roster_and_doc[0], tmp_path / "second")
    path, seal = _seal(run2, 1)
    _p0, seal0 = _seal(run2, 0)
    path.write_bytes(ep.canonical_bytes(dict(seal, meta=seal0["meta"], certs=seal0["certs"])))
    _finish_refused(run2, capsys, "block b1of2: foreign block")


def test_a_short_truncated_or_altered_block_is_refused_by_name(
        tmp_path, artifact, roster_and_doc, stubbed, capsys):
    run = _sealed_run(artifact, roster_and_doc[0], tmp_path)
    path, seal = _seal(run, 1)
    meta = run / "attempts" / seal["meta"]["file"]
    certs = run / "attempts" / seal["certs"]["file"]
    good = {f: f.read_bytes() for f in (path, meta, certs)}
    assert not os.access(meta, os.W_OK) and not os.access(certs, os.W_OK)
    for f in (meta, certs):
        os.chmod(f, 0o644)

    def restore():
        for f, data in good.items():
            f.write_bytes(data)

    # altered after sealing: the sealed digest no longer holds
    meta.write_bytes(gzip.compress(gzip.decompress(good[meta])[:-30]))
    _finish_refused(run, capsys, "block b1of2: its meta file", "is not the sealed bytes")
    restore()
    # sealed without its last record (items: header, records, subsets, trailer)
    _rewrite_meta(run, 1, lambda items: items[:-3] + items[-2:])
    _finish_refused(run, capsys, "block b1of2: short: 1 of its")
    restore()
    # sealed without its trailer
    _rewrite_meta(run, 1, lambda items: items[:-1])
    _finish_refused(run, capsys, "block b1of2: truncated")
    restore()
    # a certificate stream that is not the meta stream's records
    items = [i for i in ep._stream(certs) if i is not None]
    with gzip.open(certs, "wb") as handle:
        for item in [items[1], items[0]] + items[2:]:
            pickle.dump(item, handle)
    seal["certs"]["sha256"] = ep.file_sha256(certs)
    path.write_bytes(ep.canonical_bytes(seal))
    _finish_refused(run, capsys, "block b1of2: its certificates are not its records' in order")
    restore()
    assert emit.main(["prepare-finish", str(run), "--shards", "2"]) == 0


def test_a_block_file_replaced_after_its_check_is_not_read(tmp_path, artifact, roster_and_doc,
                                                           records, stubbed, monkeypatch):
    """The finishing step hashes each block file and unpickles that same open
    file: one replaced on disk after the check (while the step runs) changes
    nothing it publishes."""
    roster = roster_and_doc[0]
    stubbed.failing.update(_ids(records, "RMG family estimate")[:3])
    run = _sealed_run(artifact, roster, tmp_path)
    _p, seal = _seal(run, 1)
    certs = run / "attempts" / seal["certs"]["file"]
    read = ep.read_blocks
    swapped = []

    def then_swap(*args, **kwargs):
        out = read(*args, **kwargs)
        items = [i for i in ep._stream(certs) if i is not None]
        forged = tmp_path / "forged.pkl.gz"
        with gzip.open(forged, "wb") as handle:
            for item in items[:-1]:
                rc = dict(item["certificate"], certified=False, gaps=["forged"])
                pickle.dump(dict(item, certificate=rc), handle)
            pickle.dump(items[-1], handle)
        os.replace(forged, certs)
        swapped.append(True)
        return out

    monkeypatch.setattr(ep, "read_blocks", then_swap)
    manifest, generation = ep.prepare_finish(run, 2)
    assert swapped
    _identical(_serial(artifact, roster, tmp_path), generation)
    monkeypatch.setattr(ep, "read_blocks", read)
    # read afresh, the forged file fails its sealed digest
    with pytest.raises(ep.BlockRefusal, match="is not the sealed bytes"):
        ep.prepare_finish(run, 2)


def test_a_block_killed_mid_write_never_seals(tmp_path, artifact, roster_and_doc, stubbed,
                                              monkeypatch, capsys):
    roster = roster_and_doc[0]
    run = tmp_path / "run"
    _plan(artifact, roster, run, 2)
    assert ep.prepare_block(run, 0)["state"] == "sealed"

    class Killed(StubCertifier):
        def certify_record(self, record, rate_owner=None):
            if len(self.seen) == 3:
                os.kill(os.getpid(), signal.SIGKILL)
            return super().certify_record(record, rate_owner)

    monkeypatch.setattr(emit, "make_certifier", lambda *a, **k: Killed())
    pid = os.fork()
    if pid == 0:   # the block process
        try:
            ep.prepare_block(run, 1)
        finally:
            os._exit(0)
    _pid, status = os.waitpid(pid, 0)
    assert os.WIFSIGNALED(status) and os.WTERMSIG(status) == signal.SIGKILL
    assert list((run / "attempts").glob("b1of2.*.certs.pkl.gz.part"))
    assert not list((run / "attempts").glob("b1of2.*.certs.pkl.gz"))
    _finish_refused(run, capsys, "block b1of2 is missing", "running, or killed (no seal)")
    # a retry seals, and the run publishes the serial inputs
    monkeypatch.setattr(emit, "make_certifier", lambda *a, **k: StubCertifier())
    assert ep.prepare_block(run, 1)["state"] == "sealed"
    assert emit.main(["prepare-finish", str(run), "--shards", "2"]) == 0
    _identical(_serial(artifact, roster, tmp_path), run / "current")


def test_two_attempts_of_one_block_seal_once(tmp_path, artifact, roster_and_doc, stubbed,
                                             monkeypatch):
    """Attempt B passes the start check; attempt A runs and seals while B
    certifies (the condition changes while B waits); B's seal is refused at
    the atomic link and it reports superseded.  A later run does no work."""
    roster = roster_and_doc[0]
    run = tmp_path / "run"
    _plan(artifact, roster, run, 2)
    assert ep.prepare_block(run, 0)["state"] == "sealed"
    certify = ep.certify_range
    started, nested = [], []

    def racing(*args, **kwargs):
        if not started:
            started.append(True)
            nested.append(ep.prepare_block(run, 1))   # attempt A, start to seal
        return certify(*args, **kwargs)

    monkeypatch.setattr(ep, "certify_range", racing)
    b = ep.prepare_block(run, 1)
    monkeypatch.setattr(ep, "certify_range", certify)
    a = nested[0]
    assert (a["state"], b["state"]) == ("sealed", "superseded")
    _p, seal = _seal(run, 1)
    assert seal["meta"]["file"] == a["meta"]["file"] != b["meta"]["file"]
    states = sorted(json.loads(p.read_text())["state"]
                    for p in (run / "attempts").glob("b1of2.*.status.json"))
    assert states == ["sealed", "superseded"]
    files = sorted(p.name for p in (run / "attempts").iterdir())
    assert ep.prepare_block(run, 1)["state"] == "already-sealed"
    assert sorted(p.name for p in (run / "attempts").iterdir()) == files
    assert emit.main(["prepare-finish", str(run), "--shards", "2"]) == 0
    _identical(_serial(artifact, roster, tmp_path), run / "current")


def test_a_block_that_died_of_a_resource_error_blocks_publication(
        tmp_path, artifact, roster_and_doc, records, stubbed, monkeypatch, capsys):
    roster = roster_and_doc[0]
    run = tmp_path / "run"
    _plan(artifact, roster, run, 2)
    assert ep.prepare_block(run, 0)["state"] == "sealed"

    class Exhausted(StubCertifier):
        def certify_record(self, record, rate_owner=None):
            if len(self.seen) == 2:
                raise MemoryError("the block ran out of memory")
            return super().certify_record(record, rate_owner)

    monkeypatch.setattr(emit, "make_certifier", lambda *a, **k: Exhausted())
    capsys.readouterr()
    assert emit.main(["prepare-block", str(run), "--block", "1"]) == ep.RESOURCE_EXIT
    assert "BLOCK FAILED resource MemoryError" in capsys.readouterr().out
    assert not (run / "seals" / "b1of2.json").exists()
    (status,) = [json.loads(p.read_text()) for p in (run / "attempts").glob("b1of2.*.status.json")]
    assert status["state"] == "failed" and status["resource"] is True
    _finish_refused(run, capsys, "block b1of2 is missing", "failed: resource MemoryError")


def test_a_plan_is_written_once_and_recomputed_by_every_stage(
        tmp_path, artifact, roster_and_doc, stubbed, capsys):
    roster = roster_and_doc[0]
    run = tmp_path / "run"
    _plan(artifact, roster, run, 2)
    data = (run / ep.PLAN).read_bytes()
    _plan(artifact, roster, run, 2)   # the same plan: a no-op
    assert (run / ep.PLAN).read_bytes() == data
    capsys.readouterr()
    assert emit.main(["prepare-plan", str(artifact), str(roster), "--database", str(DATABASE),
                      "--blocks", "3", "--run", str(run)]) == 1
    assert "holds another plan" in capsys.readouterr().out
    assert emit.main(["prepare-plan", str(artifact), str(roster), "--database", str(DATABASE),
                      "--blocks", "99", "--run", str(tmp_path / "x")]) == 1
    assert "at most one per record" in capsys.readouterr().out
    # an edited plan steers nothing
    plan = json.loads(data)
    plan["ranges"][0]["hi"] -= 1
    (run / ep.PLAN).write_bytes(ep.canonical_bytes(plan))
    assert emit.main(["prepare-block", str(run), "--block", "0"]) == 1
    assert "not the plan of its inputs" in capsys.readouterr().out
    (run / ep.PLAN).write_bytes(data)
    assert ep.prepare_block(run, 0)["state"] == "sealed"
    assert ep.prepare_block(run, 1)["state"] == "sealed"
    # the artifact's bytes changed after the blocks: the finishing step refuses
    good = artifact.read_bytes()
    artifact.write_bytes(good.replace(b'"schema_version":', b'"schema_version": ', 1))
    assert emit.main(["prepare-finish", str(run), "--shards", "2"]) == 1
    assert "not the plan of its inputs" in capsys.readouterr().out
    artifact.write_bytes(good)
    assert emit.main(["prepare-finish", str(run), "--shards", "2"]) == 0


def test_a_published_generation_is_never_rewritten(tmp_path, artifact, roster_and_doc, records,
                                                   stubbed):
    roster = roster_and_doc[0]
    run = _blocks(artifact, roster, tmp_path, 2, shards=2)
    first = (run / "current").resolve()
    before = {p.name: p.read_bytes() for p in first.iterdir()}
    assert all(not os.access(p, os.W_OK) for p in first.iterdir())
    manifest, second = ep.prepare_finish(run, 3)
    assert second != first and (run / "current").resolve() == second
    assert {p.name: p.read_bytes() for p in first.iterdir()} == before
    # an emission shard reads the generation it was given, unchanged
    rows = tmp_path / "rows.0of2.jsonl"
    assert emit.main([str(artifact), str(roster), "--database", str(DATABASE), "--out", str(rows),
                      "--chunk", "0/2", "--prepared", str(first)]) == 0
    serial = _serial(artifact, roster, tmp_path)
    serial_rows = tmp_path / "serial.0of2.jsonl"
    assert emit.main([str(artifact), str(roster), "--database", str(DATABASE), "--out",
                      str(serial_rows), "--chunk", "0/2", "--prepared", str(serial)]) == 0
    assert rows.read_bytes() == serial_rows.read_bytes()


# --- real certification, each block in its own process -------------------------
def _write_artifact(directory: Path, doc) -> Path:
    raw = json.dumps(doc, sort_keys=True, separators=(",", ":")).encode()
    path = directory / f"{hashlib.sha256(raw).hexdigest()}.json"
    path.write_bytes(raw)
    return path


def _python(*args):
    env = dict(os.environ, PYTHONPATH=str(REPOSITORY), PYTHONHASHSEED="0")
    done = subprocess.run([sys.executable, "-m", "rmgpy.kmc.emit", *map(str, args)], env=env,
                          cwd=str(REPOSITORY), capture_output=True, text=True, timeout=3000)
    assert done.returncode == 0, done.stdout[-2000:] + done.stderr[-4000:]
    return done.stdout


@pytest.fixture(scope="module")
def real_fixtures(tmp_path_factory, roster_and_doc):
    out = tmp_path_factory.mktemp("real")
    pe = out / "pe"
    pe.mkdir()
    pe_artifact = _write_artifact(pe, roster_and_doc[1])
    pp = out / "pp"
    pp.mkdir()
    # the pk-i261 PP short compile (24 R_Recombination records), its database
    # file list dropped from the provenance
    pp_artifact = _write_artifact(pp, json.loads(gzip.decompress(PP_SHORT.read_bytes())))
    pp_roster = pp / "roster.json"
    pp_roster.write_bytes(gzip.decompress(PP_ROSTER.read_bytes()))
    return {"pe_toy": (pe_artifact, roster_and_doc[0]), "pp_short": (pp_artifact, pp_roster)}


@needs_database
@pytest.mark.parametrize("name", ["pe_toy", "pp_short"])
def test_real_blocks_in_fresh_processes_give_the_serial_inputs(tmp_path, real_fixtures, name):
    """Real certification, serial and every block each in a fresh process:
    J = 1, 2, 4, 7 against serial.  Only certificate node atom order may
    differ (no presort); the premise is refuted by any other difference."""
    artifact, roster = real_fixtures[name]
    serial = tmp_path / "serial"
    _python("prepare", artifact, roster, "--database", DATABASE, "--shards", 2, "--out", serial)
    moved = {}
    for blocks in (1, 2, 4, 7):
        run = tmp_path / f"run{blocks}"
        _python("prepare-plan", artifact, roster, "--database", DATABASE, "--blocks", blocks,
                "--run", run)
        for j in range(blocks):
            assert "sealed" in _python("prepare-block", run, "--block", j)
        manifest, generation = ep.prepare_finish(run, 2)
        result = _identical(serial, generation)
        moved[blocks] = len(result["node_order"])
        if blocks == 1:
            assert result["node_order"] == []   # one block is the serial order
    print(f"{name}: certificates equal only up to node atom order, by J: {moved}")


def test_the_comparison_reports_a_difference_outside_node_order(tmp_path, artifact, roster_and_doc,
                                                                records, stubbed):
    """``compare_prepared`` is the identity instrument: a changed certificate
    field, a changed record field or a dropped item is a difference, never
    node order."""
    roster = roster_and_doc[0]
    stubbed.failing.update(_ids(records, "RMG family estimate")[:2])
    serial = _serial(artifact, roster, tmp_path)
    run = _blocks(artifact, roster, tmp_path, 2)
    generation = (run / "current").resolve()
    manifest = json.loads((generation / ep.MANIFEST).read_text())
    shard = generation / manifest["files"][0]

    def altered(change):
        items = _items(generation)[0]
        os.chmod(shard, 0o644)
        with gzip.open(shard, "wb") as handle:
            for item in change(items) + [{"end": len(items)}]:
                pickle.dump(item, handle)
        result = ep.compare_prepared(serial, generation)
        with gzip.open(shard, "wb") as handle:
            for item in items + [{"end": len(items)}]:
                pickle.dump(item, handle)
        return result

    own = next(k for k, item in enumerate(_items(generation)[0]) if item["own"])

    def certified(items):
        item = dict(items[own], certificate=dict(items[own]["certificate"], certified=False,
                                                 gaps=["altered"]))
        return items[:own] + [item] + items[own + 1:]

    def status(items):
        record = dict(items[own]["record"], status_reason="altered")
        return items[:own] + [dict(items[own], record=record)] + items[own + 1:]

    for change, path in ((certified, ".certificate.gaps"),
                         (status, ".record.status_reason"),
                         (lambda items: items[:-1], "shard0[len]")):
        result = altered(change)
        assert result["node_order"] == []
        assert any(path in d for d in result["differences"]), result["differences"]
    _identical(serial, generation)
