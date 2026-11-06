"""Emission shards own whole reference-case keys, and their inputs are bound
to the manifest (I-304).

A prepared shard builds only the reference cases its own records read
(``emit.case_filter``), which bounds its memory only when one key's records
sit on one shard: ``sorted(ids)[k::n]`` scattered every key over every shard.
``emit_prepare.assign_shards`` places each key's records on one shard,
largest estimated cost first on the least loaded shard, and splits a key
costing more than a shard's share explicitly.  The manifest
(``kmc-emit-prepared/2``) records the assignment, each shard input's sha256,
checked on the bytes read before anything is unpickled, and a digest of each
shard's (source id, published id) pairs; a ``kmc-emit-prepared/1`` manifest is
refused by name.
"""
import gzip
import json
import pickle
from collections import Counter

import pytest

from rmgpy.kmc import emit
from rmgpy.kmc import emit_prepare as ep
from rmgpy.kmc.emit import PreparedRefusal

from emitShardMemoryTest import (_prepare, _refused, artifact,  # noqa: F401
                                 roster_and_doc, stubbed)


def _units(sizes, subsets=None):
    subsets = subsets or {}
    out, n = [], 0
    for name, size in sizes.items():
        out.append(ep.Unit(name, [f"evt_{name}_{i:04d}" for i in range(size)],
                           subsets.get(name, 0)))
        n += size
    return out


def _check_placed(units, shard_of, plan, shards):
    ids = [i for u in units for i in u.ids]
    assert sorted(shard_of) == sorted(ids)
    assert all(0 <= k < shards for k in shard_of.values())
    split = {s["unit"]: s for s in plan["split"]}
    for u in units:
        on = {shard_of[i] for i in u.ids}
        if u.name in split:
            assert sorted(on) == split[u.name]["shards"]
            assert len(on) == split[u.name]["parts"] > 1
        else:
            assert len(on) == 1, f"unit {u.name} scattered over {sorted(on)}"


def test_equal_singleton_units_are_the_round_robin_of_their_sorted_ids():
    units = _units({f"u{i:02d}": 1 for i in range(11)})
    shard_of, plan = ep.assign_shards(units, 3)
    ids = sorted(i for u in units for i in u.ids)
    assert [shard_of[i] for i in ids] == [p % 3 for p in range(len(ids))]
    assert plan["split"] == [] and plan["case_bound"] == []


def test_a_key_stays_on_one_shard_and_shards_balance():
    sizes = {f"k{i:02d}": s for i, s in enumerate([9, 7, 7, 5, 4, 4, 3, 2, 2, 1, 1, 1])}
    units = _units(sizes, {"k03": 40, "k09": 20})
    shard_of, plan = ep.assign_shards(units, 4)
    _check_placed(units, shard_of, plan, 4)
    assert plan["split"] == []
    loads = Counter()
    for u in units:
        loads[shard_of[u.ids[0]]] += ep.unit_cost(len(u.ids), u.subsets)
    largest = max(ep.unit_cost(len(u.ids), u.subsets) for u in units)
    assert max(loads.values()) - min(loads.values()) <= largest
    assert [round(loads[k], 3) for k in range(4)] == plan["loads"]
    # deterministic in the units, whatever their order
    again, plan2 = ep.assign_shards(list(reversed(units)), 4)
    assert again == shard_of and plan2 == plan


def test_a_skewed_key_too_large_for_one_shard_is_split_and_named():
    """One key holds 70% of the records: it cannot sit on one shard without
    that shard carrying more than twice its share, so it is split over
    several shards, each run carrying the key's case cost, and named."""
    sizes = {"big": 140, **{f"s{i:02d}": 3 for i in range(20)}}
    units = _units(sizes, {"big": 10})
    shard_of, plan = ep.assign_shards(units, 4)
    _check_placed(units, shard_of, plan, 4)
    (entry,) = plan["split"]
    assert entry["unit"] == "big" and entry["records"] == 140 and entry["parts"] >= 2
    share = plan["share"]
    assert max(plan["loads"]) <= 2 * share
    # the runs are contiguous in the key's sorted ids
    big = sorted(next(u for u in units if u.name == "big").ids)
    runs = [shard_of[i] for i in big]
    assert sum(1 for a, b in zip(runs, runs[1:]) if a != b) == entry["parts"] - 1
    # the runs of a split key sit on distinct shards, even where the other
    # keys leave one shard far less loaded than the rest
    units = _units({"big": 400, "mid": 10, "a": 1}, {"big": 0})
    for shards in (2, 3, 5, 8):
        shard_of, plan = ep.assign_shards(units, shards)
        _check_placed(units, shard_of, plan, shards)
        (entry,) = plan["split"]
        assert len(set(entry["shards"])) == len(entry["shards"]) == entry["parts"]
    # a one-record key never splits, and a key whose cases alone exceed a
    # shard's share is named as case-bound
    units = _units({"one": 1, "x": 2, "y": 2}, {"one": 1000})
    shard_of, plan = ep.assign_shards(units, 3)
    _check_placed(units, shard_of, plan, 3)
    assert plan["split"] == [] and plan["case_bound"] == ["one"]


def test_more_shards_than_units_leaves_shards_empty_and_names_must_be_unique():
    units = _units({"a": 2, "b": 1})
    shard_of, plan = ep.assign_shards(units, 5)
    _check_placed(units, shard_of, plan, 5)
    with pytest.raises(ValueError, match="share a name"):
        ep.assign_shards(units + _units({"a": 1}), 2)


def _manifest(prepared):
    return json.loads((prepared / ep.MANIFEST).read_text())


def _write(prepared, manifest):
    (prepared / ep.MANIFEST).write_text(json.dumps(manifest))


def test_an_old_manifest_is_refused_by_name(
        tmp_path, artifact, roster_and_doc, stubbed, capsys):  # noqa: F811
    roster = roster_and_doc[0]
    prepared = _prepare(artifact, roster, tmp_path, 2)
    manifest = _manifest(prepared)
    assert manifest["schema"] == "kmc-emit-prepared/2"
    manifest["schema"] = "kmc-emit-prepared/1"
    _write(prepared, manifest)
    out = _refused(artifact, roster, tmp_path, capsys, "--chunk", "0/2", "--prepared",
                   str(prepared))
    assert "kmc-emit-prepared/1 manifest" in out and "prepare again" in out


def test_a_shard_input_is_checked_against_its_sha256_before_it_is_unpickled(
        tmp_path, artifact, roster_and_doc, stubbed, monkeypatch, capsys):  # noqa: F811
    roster = roster_and_doc[0]
    prepared = _prepare(artifact, roster, tmp_path, 2)
    manifest = _manifest(prepared)
    for k, name in enumerate(manifest["files"]):
        assert ep.file_sha256(prepared / name) == manifest["files_sha256"][k]
    shard = prepared / manifest["files"][1]
    whole = shard.read_bytes()
    # a well-formed input of other content: it unpickles fine, so only the
    # digest can refuse it, and nothing may be unpickled first
    items = []
    with gzip.open(shard, "rb") as handle:
        while True:
            try:
                items.append(pickle.load(handle))
            except EOFError:
                break
    items[0]["record"] = dict(items[0]["record"], template="altered")
    shard.write_bytes(gzip.compress(b"".join(pickle.dumps(i) for i in items)))
    loads = []
    real = ep.pickle.load
    monkeypatch.setattr(ep.pickle, "load", lambda handle: loads.append(1) or real(handle))
    out = _refused(artifact, roster, tmp_path, capsys, "--chunk", "1/2", "--prepared",
                   str(prepared))
    assert "sha256 is not the manifest's" in out
    assert loads == []
    shard.write_bytes(whole)
    assert emit.main([str(artifact), str(roster), "--database", "/x", "--out",
                      str(tmp_path / "ok.jsonl"), "--chunk", "1/2", "--prepared",
                      str(prepared)]) == 0
    assert loads


def test_the_reader_checks_the_writers_assignment_and_published_ids(
        tmp_path, artifact, roster_and_doc, stubbed, capsys):  # noqa: F811
    roster = roster_and_doc[0]
    stubbed.failing.update(r["event_id"] for r in roster_and_doc[1]["records"][:2])
    prepared = _prepare(artifact, roster, tmp_path, 2)
    manifest = _manifest(prepared)
    shard_of = manifest["assignment"]["shard_of"]
    assert sorted(Counter(shard_of)) == [0, 1]
    assert manifest["own"] == [shard_of.count(0), shard_of.count(1)]
    # every shard reads exactly what the manifest assigns it
    for k in range(2):
        got = ep.read_shard(prepared, manifest, f"{k}/2")
        assert sorted(emit.source_ids(got.picked, got.replacements)) == ep.assigned_ids(manifest, k)
        assert got.picked == sorted(got.picked, key=lambda r: emit.source_ids(
            [r], got.replacements)[0])
    assert any(ep.read_shard(prepared, manifest, f"{k}/2").replacements for k in range(2))
    # an id moved to the other shard in the manifest only
    moved = dict(manifest, assignment=dict(manifest["assignment"], shard_of=[
        1 - shard_of[0]] + shard_of[1:]))
    with pytest.raises(PreparedRefusal, match="does not hold the records the manifest assigns"):
        ep.read_shard(prepared, moved, f"{shard_of[0]}/2")
    # the published ids differ from the writer's
    for k in range(2):
        other = dict(manifest, own_sha256=list(manifest["own_sha256"]))
        other["own_sha256"][k] = "0" * 64
        with pytest.raises(PreparedRefusal, match="source and published ids"):
            ep.read_shard(prepared, other, f"{k}/2")
    # an assignment naming a shard that does not exist
    bad = dict(manifest, assignment=dict(manifest["assignment"], shard_of=[2] + shard_of[1:]))
    with pytest.raises(PreparedRefusal, match="does not place every artifact id"):
        ep.read_shard(prepared, bad, "0/2")


# -- the prefilter, on real records, families and roster ----------------------


def _case_content(case, specs):
    """Everything a kept case carries into emission, in comparable form."""
    w = case.witness
    return (case.label, case.family, case.direction, case.pathway, case.support, case.species,
            case.kind, repr(sorted(case.wanted, key=repr)), repr(case.moves), repr(w),
            repr(w.extension_words), [m.to_adjacency_list() for m in case.reactants],
            [[f.to_adjacency_list() for f in forms] for forms in case.forms], repr(specs))


@pytest.fixture(scope="module")
def pe_cases():
    import gzip as gz
    import os
    from pathlib import Path

    from rmgpy.kmc import channels as chm
    from rmgpy.kmc import run as kr
    from rmgpy.kmc import windows as win
    from rmgpy.kmc.event_record import EventRecord
    database = Path(os.environ.get("RMG_DATABASE_PATH", "/nonexistent"))
    if not (database / "input/kinetics/families").is_dir():
        pytest.skip("RMG_DATABASE_PATH is not a database checkout")
    here = Path(__file__).resolve().parent
    spec = kr.read_deck(here / "fixtures" / "decks" / "pe_toy_deck.py").spec
    systems, refusals = emit.build_reference_systems(emit.declared_reference_systems(spec, 3))
    assert not refusals
    doc = json.loads(gz.decompress(
        (here / "fixtures" / "i295" / "pe_toy_subset_artifact.json.gz").read_bytes()))
    records = [EventRecord.from_dict(r).to_dict() for r in doc["records"]]
    db = emit.load_database(sorted({r["family"] for r in records}), database)
    fam = db.kinetics.families
    declared = emit.declared_provider(spec)
    policy = win.declared_interim(declared.declared_elements())
    certifier = emit.make_certifier(fam, {"domain": win.certificate_domain(policy),
                                          "repeat": declared.unit_formula(),
                                          "growth": declared}, emit.reference_thermo(db))
    full, _stats = emit.reference_cases(certifier, fam, systems, chm.SpeciesKeys())
    return certifier, fam, systems, records, full


def test_the_prefilter_keeps_every_case_its_records_read_in_order(pe_cases):
    """The cases a filtered build keeps under each of the records' keys are the
    unfiltered build's, in order and in full content (participant metadata,
    construction ids, specs, mapping multiplicity and so ``found[0]``).  The
    records' keys read forward and reverse cases (a reverse case's
    participants are its other side), keys with several cases of one species
    and keys whose cases carry different participant metadata; most subsets
    are skipped before ``member_cases``."""
    from rmgpy.kmc import channels as chm
    certifier, fam, systems, records, full = pe_cases
    species = chm.SpeciesKeys()
    wanted = emit.case_filter(certifier, records, species)
    kept, stats = emit.reference_cases(certifier, fam, systems, species, wanted=wanted)
    assert set(kept) == wanted.keys & set(full)
    assert set(kept) == wanted.keys, "every fixture record's key has reference cases"
    for key in wanted.keys:
        assert [_case_content(*c) for c in kept[key]] == [_case_content(*c) for c in full[key]]
    directions = Counter(c.direction for key in kept for c, _s in kept[key])
    assert directions["forward"] and directions["reverse"]
    assert any(len(kept[k]) > 1 and len({c.species for c, _ in kept[k]}) == 1 for k in kept)
    assert any(len({repr(s) for _, s in kept[k]}) > 1 for k in kept), \
        "no key whose cases carry different participant metadata"
    assert stats["subsets_skipped"] > stats["subsets_enumerated"] > 0
    # the condition is direction-aware: a filter reading the template side for
    # reverse cases (the generated products) loses those cases
    swapped = emit.CaseFilter(wanted.keys, {
        (f, d): wanted.sides.get((f, "reverse" if d == "forward" else "forward"), frozenset())
        for f, d in wanted.sides})
    lost, _ = emit.reference_cases(certifier, fam, systems, species, wanted=swapped)

    def reverse(cases):
        return sum(c.direction == "reverse" for found in cases.values() for c, _s in found)
    assert reverse(lost) < reverse(kept)


def test_records_whose_key_cannot_be_read_add_nothing_to_the_filter(pe_cases):
    from rmgpy.kmc import channels as chm
    certifier, _fam, _systems, records, _full = pe_cases
    broken = dict(records[0], **{k: ["not a graph"] for k in records[0]
                                 if k.endswith("graphs") or k in ("reactants", "products")})
    species = chm.SpeciesKeys()
    assert emit.case_filter(certifier, [broken], species).keys == frozenset()
    assert emit.case_filter(certifier, [broken] + records[:1], species).keys == \
        emit.case_filter(certifier, records[:1], species).keys
