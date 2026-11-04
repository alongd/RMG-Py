"""An emission shard's memory is its own records', not the artifact's (I-295).

A shard used to parse the whole artifact and certify every record, then emit
its ``1/N`` slice: memory and time grew with the whole catalog, and a 5.3 GB
artifact failed every shard in the JSON parse.  Certification cannot be split
per shard (a certificate's bytes depend on every record certified before it:
RMG re-sorts tree-group atoms in place), so ``emit prepare`` certifies once,
streaming, in artifact order, and writes each shard's input.

- a prepared shard never parses the artifact or certifies; it holds only its
  own records and the records they link (for reverse pairs: their partners);
- its rows file (header and rows) is byte-identical to the whole-artifact
  shard's, for every shard count, with certification refusing records so
  that ids are republished and refusals propagate across shard boundaries;
  with ``--family`` and ``--limit`` too; and from an index written once;
- prepare certifies in the whole-artifact shard's order, the artifact's, also
  when that is not id order;
- a prepared input of other artifact bytes, another roster or shard count,
  or a truncated one, is refused by name.

The records are real: 9 reverse pairs and one orphan (its partner absent) of
a PE toy compile (``fixtures/i295/pe_toy_subset_artifact.json.gz``, the head
cut to ``provenance.proxy_source``).  Certification is a stub that refuses
chosen records, so the certification steps republish them and their partners
for real; emission is a stub row digesting what emission reads (the published
record, its certificate, its rate owner's published form, the upstream
refusal).  The last test certifies for real against the database.
"""
import gzip
import hashlib
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from rmgpy.kmc import artifact_stream as astream
from rmgpy.kmc import emit_prepare
from rmgpy.kmc import certificates as cm
from rmgpy.kmc import emit
from rmgpy.kmc import template_emission as te

HERE = Path(__file__).resolve().parent
FIXTURE = HERE / "fixtures" / "i295" / "pe_toy_subset_artifact.json.gz"
DECK = HERE / "fixtures" / "decks" / "pe_toy_deck.py"
DATABASE = Path(os.environ.get("RMG_DATABASE_PATH", "/nonexistent"))
needs_database = pytest.mark.skipif(
    not (DATABASE / "input/kinetics/families").is_dir(),
    reason="RMG_DATABASE_PATH is not a database checkout")


def _canon(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


@pytest.fixture(scope="module")
def roster_and_doc(tmp_path_factory):
    from rmgpy.kmc.run import read_deck
    out = tmp_path_factory.mktemp("roster")
    roster = emit.write_roster(emit.declared_reference_systems(read_deck(str(DECK)).spec, 3),
                               out / "roster.json")
    doc = json.loads(gzip.decompress(FIXTURE.read_bytes()))
    # bind the compiled records to this checkout's declaration of the same deck
    doc["provenance"]["proxy_source"].update(roster.digests())
    return out / "roster.json", doc


@pytest.fixture
def artifact(tmp_path, roster_and_doc):
    _roster, doc = roster_and_doc
    raw = json.dumps(doc, sort_keys=True, separators=(",", ":")).encode()
    path = tmp_path / "art" / f"{hashlib.sha256(raw).hexdigest()}.json"
    path.parent.mkdir()
    path.write_bytes(raw)
    return path


@pytest.fixture
def records(roster_and_doc):
    return roster_and_doc[1]["records"]


class StubCertifier:
    """Certifies every record but ``failing`` (and a reverse without its rate
    owner), with the owner block ``certify_record`` builds; records what it
    certified."""

    def __init__(self, failing=()):
        self.failing, self.seen = set(failing), []

    def certify_record(self, record, rate_owner=None):
        d = cm._as_dict(record)
        self.seen.append(d["event_id"])
        owner = cm.record_owner(d)
        kind = (d.get("rate_source") or {}).get("kind")
        ro = d if kind == "RMG family estimate" else (
            cm._as_dict(rate_owner) if rate_owner is not None else None)
        owner.update(rate_owner_event_id=ro and ro["event_id"],
                     rate_source_sha256=ro and cm._sha(ro.get("rate_source") or {}))
        refused = d["event_id"] in self.failing or ro is None
        return {"schema": cm.RECORD_CERTIFICATE_SCHEMA, "owner": owner,
                "event_id": d["event_id"], "certified": not refused,
                "gaps": [f"stub refused {d['event_id'][:16]}"] if refused else [],
                "owner_rate_source": ro and ro.get("rate_source")}


def digest_row(self, record):
    """What emission reads of a record, digested."""
    S = self.state
    rc = emit.record_certificate(record, S["store"])
    owner = cm._rate_owner(record, S["by_id"])
    seen = {"record": record, "certificate": rc, "rate_owner": owner,
            "partner": S["by_id"].get(record.get("reverse_of")),
            "upstream": te.upstream_refusal(record, rc, S["by_id"])}
    return {"kind": "row", "event_id": record["event_id"][:24], "family": record["family"],
            "orientation": record["orientation"], "status": "refused",
            "digest": hashlib.sha256(_canon(seen).encode()).hexdigest(),
            "upstream": seen["upstream"], "seconds": 0.0, "rss_mb": 0}


@pytest.fixture
def stubbed(monkeypatch):
    """No database, reference cases or emission: certification by ``StubCertifier``."""
    from rmgpy.kmc import reference_thermo
    state = SimpleNamespace(failing=set(), certifiers=[], read=[])
    db = SimpleNamespace(kinetics=SimpleNamespace(families={}), thermo=None)

    def certifier(*args, **kwargs):
        made = StubCertifier(state.failing)
        state.certifiers.append(made)
        return made

    read = astream.read_records

    def spy(path, entries):
        entries = list(entries)
        state.read.append([e.event_id for e in entries])
        return read(path, entries)

    monkeypatch.setattr(emit, "load_database", lambda families, path: db)
    monkeypatch.setattr(emit, "make_certifier", certifier)
    monkeypatch.setattr(emit, "reference_thermo", lambda db: None)
    monkeypatch.setattr(emit, "reference_cases", lambda *a, **k: ({}, {}))
    monkeypatch.setattr(emit, "build_reference_systems", lambda roster: ({}, []))
    monkeypatch.setattr(emit, "database_content_sha256", lambda path: "h" * 64)
    monkeypatch.setattr(emit, "RatesByGrid", lambda *a, **k: None)
    monkeypatch.setattr(te, "RmgCounter", lambda fam: None)
    monkeypatch.setattr(te, "ThermoKc", lambda thermo: None)
    monkeypatch.setattr(te, "ProductionRates", lambda *a, **k: None)
    monkeypatch.setattr(reference_thermo, "GasPhaseRMGReferenceThermo", lambda *a: None)
    monkeypatch.setattr(astream, "read_records", spy)
    monkeypatch.setattr(emit.Emitter, "run_record", digest_row)
    return state


def _emit(artifact, roster, out, *flags):
    assert emit.main([str(artifact), str(roster), "--database", str(DATABASE), "--out", str(out),
                      *flags]) == 0
    return out.read_bytes()


def _prepare(artifact, roster, tmp_path, n, *flags, name="prepared"):
    out = tmp_path / f"{name}{n}"
    assert emit.main(["prepare", str(artifact), str(roster), "--database", str(DATABASE),
                      "--shards", str(n), "--out", str(out), *flags]) == 0
    return out


def _shards(artifact, roster, tmp_path, n, tag, *flags):
    return [_emit(artifact, roster, tmp_path / f"{tag}.{k}of{n}.jsonl", "--chunk", f"{k}/{n}",
                  *flags) for k in range(n)]


def _prepared_shards(artifact, roster, tmp_path, n, *flags, prepare_flags=()):
    prepared = _prepare(artifact, roster, tmp_path, n, *flags, *prepare_flags)
    return _shards(artifact, roster, tmp_path, n, "prepared", "--prepared", str(prepared), *flags)


def _partners(records, ids):
    by_id = {r["event_id"]: r for r in records}
    return {r["reverse_of"] for i in ids for r in [by_id[i]] if r.get("reverse_of") in by_id}


@pytest.mark.parametrize("n", [2, 3])
def test_a_prepared_shard_holds_only_its_records_and_their_partners(
        tmp_path, artifact, roster_and_doc, records, stubbed, monkeypatch, n):
    roster = roster_and_doc[0]
    prepared = _prepare(artifact, roster, tmp_path, n)
    held = []
    init = emit.Emitter.init

    def spied(self, *args, **kwargs):
        init(self, *args, **kwargs)
        held.append({r["event_id"] for r in self.records})

    def never(*args, **kwargs):
        raise AssertionError("a prepared shard parsed or certified the artifact")

    read_text, read_bytes = Path.read_text, Path.read_bytes
    monkeypatch.setattr(Path, "read_text", lambda self, *a, **k: never() if self == artifact
                        else read_text(self, *a, **k))
    monkeypatch.setattr(Path, "read_bytes", lambda self: never() if self == artifact
                        else read_bytes(self))
    monkeypatch.setattr(emit, "certified_records", never)
    monkeypatch.setattr(astream, "read_records", never)
    monkeypatch.setattr(emit.Emitter, "init", spied)
    ids = sorted(r["event_id"] for r in records)
    for k in range(n):
        _emit(artifact, roster, tmp_path / f"rows.{k}of{n}.jsonl", "--chunk", f"{k}/{n}",
              "--prepared", str(prepared))
        mine = set(ids[k::n])
        assert held[-1] == mine | _partners(records, mine), \
            f"shard {k}/{n} held {len(held[-1] - mine - _partners(records, mine))} other records"


# (shard count, records certification refuses): a refused record republishes
# itself and its partner under new ids, often across a shard boundary
CASES = [(1, ()), (2, (0,)), (3, (0, 4, 8)), (5, (2, 7)), (19, (1, 3))]


@pytest.mark.parametrize("n, failing", CASES)
def test_prepared_rows_are_the_whole_artifact_rows(tmp_path, artifact, roster_and_doc, records,
                                                   stubbed, n, failing):
    roster = roster_and_doc[0]
    estimates = sorted(r["event_id"] for r in records
                       if r["rate_source"]["kind"] == "RMG family estimate")
    stubbed.failing.update(estimates[i] for i in failing)
    whole = _shards(artifact, roster, tmp_path, n, "whole")
    assert _prepared_shards(artifact, roster, tmp_path, n) == whole
    heads = [json.loads(text.splitlines()[0]) for text in whole]
    republished = sum(a != s for h in heads
                      for a, s in zip(h["assigned"], h["partition"]["sources"]))
    # every refused estimate takes its partner with it
    assert republished == 2 * len(failing)
    # the whole-artifact path certified every record in every shard; prepare
    # certified each once, and the prepared shards none
    everything = sorted(r["event_id"] for r in records)
    assert [sorted(c.seen) for c in stubbed.certifiers] == \
        [everything] * n + [everything] + [[]] * n


def test_prepare_certifies_in_the_whole_artifact_order(tmp_path, roster_and_doc, stubbed):
    """Certificate bytes depend on what was certified before (tree-group atom
    order), so prepare must certify in the order the whole-artifact shard does:
    the artifact's.  An artifact's records are not always in id order (the
    compiler appends ortho-junction records after the others), so the
    fixture's records are reversed here."""
    roster, doc = roster_and_doc
    doc = json.loads(json.dumps(doc))
    doc["records"].reverse()
    for index, record in enumerate(doc["records"]):
        record["canonical_index"] = index
    raw = json.dumps(doc, sort_keys=True, separators=(",", ":")).encode()
    artifact = tmp_path / "art" / f"{hashlib.sha256(raw).hexdigest()}.json"
    artifact.parent.mkdir()
    artifact.write_bytes(raw)
    order = [r["event_id"] for r in doc["records"]]
    assert order != sorted(order)
    whole = _shards(artifact, roster, tmp_path, 2, "whole")
    assert _prepared_shards(artifact, roster, tmp_path, 2) == whole
    whole_shard, prepare = stubbed.certifiers[0], stubbed.certifiers[2]
    assert whole_shard.seen == order
    assert prepare.seen == whole_shard.seen


@pytest.mark.parametrize("flags", [("--family", "R_Recombination"), ("--limit", "2"),
                                   ("--family", "H_Abstraction", "--family", "intra_H_migration",
                                    "--limit", "3")])
def test_family_and_limit_selections_are_the_whole_artifact_ones(
        tmp_path, artifact, roster_and_doc, records, stubbed, flags):
    roster = roster_and_doc[0]
    stubbed.failing.update(r["event_id"] for r in records[::4]
                           if r["rate_source"]["kind"] == "RMG family estimate")
    family = [f for k, f in enumerate(flags) if k and flags[k - 1] == "--family"]
    prepare_flags = [x for f in family for x in ("--family", f)]
    prepared = _prepare(artifact, roster, tmp_path, 2, *prepare_flags)
    assert _shards(artifact, roster, tmp_path, 2, "p", "--prepared", str(prepared), *flags) == \
        _shards(artifact, roster, tmp_path, 2, "w", *flags)


def test_an_index_written_once_prepares_the_same_inputs(tmp_path, artifact, roster_and_doc,
                                                        stubbed, capsys):
    roster = roster_and_doc[0]
    index = tmp_path / "artifact.index.json"
    assert emit.main(["index", str(artifact), "--out", str(index)]) == 0
    assert _prepared_shards(artifact, roster, tmp_path, 3, prepare_flags=("--index", str(index))) \
        == _shards(artifact, roster, tmp_path, 3, "w")
    other = tmp_path / "other.json"
    other.write_bytes(artifact.read_bytes().replace(b'"schema_version":', b'"schema_version": ', 1))
    capsys.readouterr()
    assert emit.main(["prepare", str(other), str(roster), "--database", str(DATABASE),
                      "--shards", "2", "--out", str(tmp_path / "x"), "--index", str(index)]) == 1
    assert "refused:index" in capsys.readouterr().out


def _refused(artifact, roster, tmp_path, capsys, *flags):
    capsys.readouterr()
    out = tmp_path / "refused.jsonl"
    assert emit.main([str(artifact), str(roster), "--database", str(DATABASE), "--out", str(out),
                      *flags]) == 1
    assert not out.exists()
    return capsys.readouterr().out


def test_a_prepared_input_that_is_not_this_shards_is_refused(
        tmp_path, artifact, roster_and_doc, stubbed, capsys):
    roster = roster_and_doc[0]
    prepared = _prepare(artifact, roster, tmp_path, 2)
    assert "prepared for 2 shards" in _refused(artifact, roster, tmp_path, capsys,
                                               "--chunk", "0/3", "--prepared", str(prepared))
    assert "prepared for families" in _refused(artifact, roster, tmp_path, capsys, "--chunk",
                                               "0/2", "--family", "R_Recombination",
                                               "--prepared", str(prepared))
    other = tmp_path / "other.json"
    other.write_bytes(artifact.read_bytes().replace(b'"schema_version":', b'"schema_version": ', 1))
    assert "another artifact" in _refused(other, roster, tmp_path, capsys, "--chunk", "0/2",
                                          "--prepared", str(prepared))
    # the same declaration, another roster file (its systems differ)
    doc = json.loads(roster.read_text())
    doc["systems"] = dict(list(doc["systems"].items())[:1])
    second = tmp_path / "roster2.json"
    second.write_text(json.dumps(doc))
    assert "roster_sha256" in _refused(artifact, second, tmp_path, capsys, "--chunk", "0/2",
                                       "--prepared", str(prepared))
    manifest = json.loads((prepared / emit_prepare.MANIFEST).read_text())
    shard = prepared / manifest["files"][1]
    whole = shard.read_bytes()
    shard.write_bytes(gzip.compress(gzip.decompress(whole)[:-40]))
    assert "truncated" in _refused(artifact, roster, tmp_path, capsys, "--chunk", "1/2",
                                   "--prepared", str(prepared))
    shard.write_bytes(whole[:len(whole) // 2] + b"\0" * 64)
    assert "refused:prepared" in _refused(artifact, roster, tmp_path, capsys, "--chunk", "1/2",
                                          "--prepared", str(prepared))
    shard.write_bytes(whole)
    # a prepare rerun removes the manifest before anything else: no shard pairs
    # the old manifest with new inputs
    from rmgpy.kmc import emit_prepare as ep
    seen = []
    replace = ep.os.replace

    def watch(src, dst):
        seen.append((Path(dst).name, (prepared / ep.MANIFEST).exists()))
        return replace(src, dst)

    ep.os.replace = watch
    try:
        _prepare(artifact, roster, tmp_path, 2)
    finally:
        ep.os.replace = replace
    assert [name for name, _ in seen][-1] == ep.MANIFEST
    assert not any(present for _, present in seen)
    (prepared / emit_prepare.MANIFEST).unlink()
    assert "has no manifest.json" in _refused(artifact, roster, tmp_path, capsys, "--chunk",
                                              "0/2", "--prepared", str(prepared))


@needs_database
def test_real_certification_prepared_once_is_every_whole_artifact_shards(
        tmp_path, artifact, roster_and_doc, monkeypatch):
    """A shard certifying fewer, other records first publishes other
    certificate bytes (tree-group atom order); prepared once, in artifact
    order, the shards' records and certificates are the whole path's."""
    roster = roster_and_doc[0]
    monkeypatch.setattr(emit, "reference_cases", lambda *a, **k: ({}, {}))
    monkeypatch.setattr(emit, "build_reference_systems", lambda roster: ({}, []))
    monkeypatch.setattr(emit, "database_content_sha256", lambda path: "h" * 64)
    monkeypatch.setattr(emit.Emitter, "run_record", digest_row)
    whole = _shards(artifact, roster, tmp_path, 2, "whole")
    assert _prepared_shards(artifact, roster, tmp_path, 2) == whole
    heads = [json.loads(text.splitlines()[0]) for text in whole]
    # the compile's real certification refuses records here: republishing is exercised
    assert any(a != s for h in heads for a, s in zip(h["assigned"], h["partition"]["sources"]))
