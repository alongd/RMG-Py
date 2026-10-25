#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Sharded emission partitions the artifact's own records (I-265).

Every shard certified the whole artifact and took ``sorted(certified ids)[k::n]``;
certification republishes a record it refuses under a new id, so a shard that
saw another certification outcome (one transient MemoryError swallowed as a
certification failure) sliced another list: records fell between shards or
landed in two.  On PP 64 shards assigned 2065 of 2386 records, 310 twice.

- the shard partition is the artifact's record ids, sorted, ``k::n``, whatever
  a shard's certification republishes; a republished record stays with the
  shard owning its source id;
- a resource error (MemoryError, RecursionError) in certification or emission
  fails the shard; it is never an error certificate, a refused record or an
  error row;
- assembly names every artifact record no shard assigns, and every one two
  shards assign.
"""
import json
from collections import Counter
from pathlib import Path

import pytest

from rmgpy.kmc import certificates as cm
from rmgpy.kmc import emit
from rmgpy.kmc import template_emission as te

from emitTest import Services, chain, fixture_rows, header, record_id, write_roster, write_rows

N = 4
IDS = [f"evt_{k:02d}{'0' * 60}" for k in range(12)]
# a record certification refuses under memory pressure, with its linked partner
FAILED = {IDS[3]: "evt_ff" + "1" * 60, IDS[7]: "evt_00" + "2" * 59 + "a"}


class _CertifyingEmitter:
    """An emitter whose certification republishes FAILED (outcome "B") or
    nothing (outcome "A"), as two shards under different memory saw it."""

    outcome = "A"

    def __init__(self, artifact, database, systems, **kwargs):
        self.artifact, self.database, self.systems = Path(artifact), database, systems
        self.records, self.state = [], {}

    def init(self):
        replacements = dict(FAILED) if self.outcome == "B" else {}
        self.records = [{"event_id": replacements.get(i, i), "family": "F",
                         "orientation": "as_generated"} for i in IDS]
        self.state.update(artifact_records=len(IDS), families=["F"],
                          replacements=replacements)

    def run_record(self, record):
        return {"kind": "row", "event_id": record["event_id"][:24], "family": "F",
                "orientation": "as_generated", "status": "refused", "refusal": "stub",
                "seconds": 0.0, "rss_mb": 0}


def _shard(tmp_path, monkeypatch, k, outcome):
    artifact = tmp_path / "artifact.json"
    artifact.write_text(json.dumps({"provenance": {}, "records": [
        {"event_id": i, "family": "F", "orientation": "as_generated"} for i in IDS]}))
    roster = write_roster(tmp_path, "pe", {"a": [chain(2)]})
    monkeypatch.setattr(emit, "build_reference_systems", lambda roster: ({}, []))
    monkeypatch.setattr(emit, "database_content_sha256", lambda path: "h" * 64)
    monkeypatch.setattr(_CertifyingEmitter, "outcome", outcome)
    monkeypatch.setattr(emit, "Emitter", _CertifyingEmitter)
    out = tmp_path / f"rows.{k}of{N}.jsonl"
    assert emit.main([str(artifact), str(roster), "--database", str(tmp_path), "--out", str(out),
                      "--chunk", f"{k}/{N}", "--legacy-artifact"]) == 0
    return json.loads(out.read_text().splitlines()[0])


def test_shards_that_certify_differently_partition_the_artifact(tmp_path, monkeypatch):
    back = {new: old for old, new in FAILED.items()}
    heads = {k: _shard(tmp_path, monkeypatch, k, "A" if k % 2 == 0 else "B") for k in range(N)}
    sources = Counter(back.get(i, i) for h in heads.values() for i in h["assigned"])
    assert sorted(sources) == IDS and set(sources.values()) == {1}, sources
    for k, head in heads.items():
        assert sorted(back.get(i, i) for i in head["assigned"]) == sorted(IDS)[k::N]
    # the republished record is emitted, under its new id, by the shard owning its source
    owner = sorted(IDS).index(IDS[3]) % N
    assert FAILED[IDS[3]] in heads[owner]["assigned"]


@pytest.mark.parametrize("error", [MemoryError, RecursionError])
def test_a_resource_error_in_certification_fails_the_shard(error):
    class Certifier:
        def certify_record(self, record, rate_owner=None):
            raise error()
    records = [{"event_id": IDS[0], "family": "F", "template": "T",
                "orientation": "as_generated", "status": "enabled"}]
    with pytest.raises(error):
        cm.certify_records(records, {}, Certifier(), required=False)


@pytest.mark.parametrize("error", [MemoryError, RecursionError])
def test_a_resource_error_in_emission_fails_the_shard(monkeypatch, error):
    def boom(*args, **kwargs):
        raise error()
    monkeypatch.setattr(te, "upstream_refusal", lambda record, rc, by_id: None)
    monkeypatch.setattr(emit, "record_key", boom)
    emitter = emit.Emitter("artifact.json", "db", {})
    emitter.state.update(certifier=None, fam={"F": None}, store={}, by_id={}, species=None,
                         cases={})
    record = {"event_id": IDS[0], "family": "F", "orientation": "as_generated",
              "template": "T", "status": "enabled"}
    with pytest.raises(error):
        emitter.run_record(record)


# -- assembly names the ids


@pytest.fixture
def stubbed(monkeypatch):
    monkeypatch.setattr(emit, "assembler_identity", lambda database: {
        "database_content_sha256": "h" * 64, "code_sources_sha256": "s" * 64,
        "code_extensions_sha256": "x" * 64})
    monkeypatch.setattr(te, "derive_weights", lambda payloads, counter: payloads)
    monkeypatch.setattr(te, "check_cycles", lambda payloads, counter, kc, complete=True, **kw:
                        (payloads, []))


def _partitioned(tmp_path, assignment):
    """Shards over the fixture rows: shard k assigns the rows ``assignment[k]``
    (indices), its header carrying the artifact's ids and each row's source."""
    rows = fixture_rows()
    ids = sorted(record_id(r) for r in rows)
    paths = []
    for k, picks in enumerate(assignment):
        mine = [rows[i] for i in picks]
        head = header(f"{k}/{len(assignment)}", mine, partition={
            "basis": "artifact_event_ids", "artifact_ids": [i[:24] for i in ids],
            "sources": [record_id(r)[:24] for r in mine]})
        paths.append(write_rows(tmp_path / f"rows.{k}of{len(assignment)}.jsonl", head, mine))
    return rows, paths


def test_assembly_names_every_record_no_shard_or_two_shards_assign(tmp_path, stubbed):
    rows, paths = _partitioned(tmp_path, [[0, 2, 4], [1, 2, 3]])
    catalog = emit.assemble([str(p) for p in paths], services=Services())
    failures = catalog["gates"]["failures"]
    lost, twice = record_id(rows[5])[:24], record_id(rows[2])[:24]
    assert any(lost in f and "no shard" in f for f in failures), failures
    assert any(twice in f and "2 shards" in f for f in failures), failures
    assert catalog["complete"] is False


def test_an_exact_partition_passes_the_assignment_gate(tmp_path, stubbed):
    _rows, paths = _partitioned(tmp_path, [[0, 2, 4], [1, 3, 5]])
    catalog = emit.assemble([str(p) for p in paths], services=Services())
    assert not [f for f in catalog["gates"]["failures"] if "shard" in f or "assign" in f]
