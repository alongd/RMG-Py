"""Tests for slice S-g: one family is the unit of work, end to end.

Database-free.  The stub catalogs are the toy PE templates of meltRunTest
split by family into two one-family catalogs; the emit selection is tested
with a fake emission session; ``compile --family`` with a spy in place of the
child process.
"""

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

import melt_engine_support as S
from meltRunTest import PE_LIVE, toy_payloads

from rmgpy.kmc import compile as kcompile
from rmgpy.kmc import emit
from rmgpy.kmc import run as kr

GRID = (300.0, 1000.0)
FAMILIES_A = ("toy_homolysis", "toy_recombination")
SHARED = {"database_path": "/db", "database_content_sha256": "h" * 64, "code_sources_sha256": "s" * 64,
          "code_extensions_sha256": "x" * 64, "roster_sha256": "r" * 64, "database_head": "d" * 40,
          "code_sha": "c" * 40}


def stub_payloads():
    """The toy PE payloads with a temperature-dependent weight, split into the
    families of catalog A and the rest (catalog B)."""
    payloads = toy_payloads([t for t in S.pe_templates() if t.name in PE_LIVE], grid=GRID)
    for k, p in enumerate(payloads):
        for t in p["templates"]:
            t["weight"]["w"] = [t["weight"]["w"][0] * (1 + k), t["weight"]["w"][0] * (3 + k)]
    a = [p for p in payloads if p["family"] in FAMILIES_A]
    b = [p for p in payloads if p["family"] not in FAMILIES_A]
    assert a and b
    return a, b


def write_stub(path, payloads, *, complete=True, entry_keys=None, drop=(), **shared):
    """A catalog of ``payloads``; ``entry_keys``: the keys that entered its
    collision resolution (default: its emitted keys); ``drop``: provenance
    fields left out."""
    provenance = dict(SHARED, **shared)
    provenance.update(artifact_id=f"art-{Path(path).stem}", families=sorted({p["family"] for p in payloads}),
                      temperatures=list(GRID))
    for field in drop:
        provenance.pop(field)
    if entry_keys is None:
        entry_keys = sorted({t["key"] for p in payloads for t in p["templates"] if t.get("status") == "emitted"})
    collisions = {} if entry_keys == "absent" else {"entry_keys": list(entry_keys)}
    doc = {"schema": "kmc-template-catalog/1", "provenance": provenance, "payloads": payloads,
           "refusals": [], "gates": {"failures": [] if complete else ["injected"], "collisions": collisions},
           "completeness": {}, "complete": complete}
    Path(path).write_text(json.dumps(doc, sort_keys=True, indent=1) + "\n")
    return Path(path)


def link_pair(a, b):
    """A reciprocal pair across the two catalogs: the first record of each."""
    a[0]["reverse_of"] = b[0]["record_event_id"]
    b[0]["reverse_of"] = a[0]["record_event_id"]


def build_set(paths, out):
    emit.build_catalog_set([str(p) for p in paths], str(out))
    return Path(out)


def loaded(paths, T):
    catalog, records = kr.load_catalogs([str(p) for p in paths], T)
    return sorted((t.name, t.weight) for t in catalog.templates(T)), records


def test_a_catalog_set_loads_as_the_two_family_catalog(tmp_path):
    a, b = stub_payloads()
    fa = write_stub(tmp_path / "a.json", a)
    fb = write_stub(tmp_path / "b.json", b)
    both = write_stub(tmp_path / "both.json", a + b)
    manifest = build_set([fa, fb], tmp_path / "set.json")
    doc = json.loads(manifest.read_text())
    assert doc["schema"] == "kmc-catalog-set/1"
    assert [m["families"] for m in doc["catalogs"]] == [sorted({p["family"] for p in a}),
                                                        sorted({p["family"] for p in b})]
    assert doc["shared"] == SHARED
    for T in (300.0, 650.0, 1000.0):
        want, _ = loaded([both], T)
        got, records = loaded([manifest], T)
        assert len(want) == len(a) + len(b)
        assert got == want
        assert len(records) == 2


def test_a_set_missing_one_side_of_a_cross_family_dependency_is_refused_by_name(tmp_path):
    a, b = stub_payloads()
    link_pair(a, b)
    fa = write_stub(tmp_path / "a.json", a)
    fb = write_stub(tmp_path / "b.json", b)
    ida, idb = a[0]["record_event_id"], b[0]["record_event_id"]
    manifest = build_set([fa, fb], tmp_path / "set.json")
    doc = json.loads(manifest.read_text())
    kinds = {(d["kind"], d["record_event_id"], d["requires_record_event_id"]) for d in doc["dependencies"]}
    assert kinds == {("reverse_of", ida, idb), ("reverse_of", idb, ida)}
    loaded([manifest], 650.0)                                  # whole: loads
    with pytest.raises(ValueError, match=idb):                 # the builder, one side only
        build_set([fa], tmp_path / "half.json")
    half = dict(doc, catalogs=doc["catalogs"][:1])             # a set with a member dropped
    (tmp_path / "half_set.json").write_text(json.dumps(half))
    with pytest.raises(kr.RunError) as caught:
        kr.load_catalogs([str(tmp_path / "half_set.json")], 650.0)
    assert ida in str(caught.value) and idb in str(caught.value) and "reverse_of" in str(caught.value)


def test_the_manifest_dependency_check_stands_for_incomplete_members(tmp_path):
    """For incomplete members (loaded with ``allow_incomplete``) a dangling
    reverse_of is not an error by itself, so the manifest's own dependency
    list is the only guard: a set with a member dropped is refused by name."""
    a, b = stub_payloads()
    link_pair(a, b)
    fa = write_stub(tmp_path / "a.json", a, complete=False)
    fb = write_stub(tmp_path / "b.json", b, complete=False)
    ida, idb = a[0]["record_event_id"], b[0]["record_event_id"]
    manifest = build_set([fa, fb], tmp_path / "set.json")
    doc = json.loads(manifest.read_text())
    half = dict(doc, catalogs=doc["catalogs"][:1])
    (tmp_path / "half_set.json").write_text(json.dumps(half))
    with pytest.raises(kr.RunError) as caught:
        kr.load_catalogs([str(tmp_path / "half_set.json")], 650.0, allow_incomplete=True)
    message = str(caught.value)
    assert "lacks" in message and ida in message and idb in message and "reverse_of" in message


def test_a_template_link_across_catalogs_is_refused(tmp_path):
    """merged_into / derived primaries load only inside their own catalog."""
    a, b = stub_payloads()
    donor = b[0]["templates"][0]
    copy = dict(donor, status="merged", merged_into={"record": b[0]["record_event_id"], "key": donor["key"]})
    a.append({"record_event_id": "evt_copy", "family": a[0]["family"], "reverse_of": None, "templates": [copy]})
    fa = write_stub(tmp_path / "a.json", a)
    fb = write_stub(tmp_path / "b.json", b)
    with pytest.raises(ValueError, match="merged_into"):
        build_set([fa, fb], tmp_path / "set.json")


def test_a_member_whose_sha_changed_is_refused(tmp_path):
    a, b = stub_payloads()
    fa = write_stub(tmp_path / "a.json", a)
    fb = write_stub(tmp_path / "b.json", b)
    manifest = build_set([fa, fb], tmp_path / "set.json")
    write_stub(tmp_path / "b.json", b + [], complete=False)    # same payloads, other bytes
    with pytest.raises(kr.RunError, match=r"b\.json.*sha256|sha256.*b\.json"):
        kr.load_catalogs([str(manifest)], 650.0)


def test_members_with_differing_shared_inputs_are_refused_by_name(tmp_path):
    a, b = stub_payloads()
    fa = write_stub(tmp_path / "a.json", a)
    fb = write_stub(tmp_path / "b.json", b, database_content_sha256="e" * 64)
    with pytest.raises(ValueError, match="database_content_sha256"):
        build_set([fa, fb], tmp_path / "set.json")
    fb = write_stub(tmp_path / "b.json", b, database_head="e" * 40)
    with pytest.raises(ValueError, match="database_head"):
        build_set([fa, fb], tmp_path / "set.json")
    fb = write_stub(tmp_path / "b.json", b)
    manifest = build_set([fa, fb], tmp_path / "set.json")
    fb = write_stub(tmp_path / "b.json", b, roster_sha256="s" * 64)   # changed after the build
    doc = json.loads(manifest.read_text())
    import hashlib
    doc["catalogs"][1]["sha256"] = hashlib.sha256(fb.read_bytes()).hexdigest()
    manifest.write_text(json.dumps(doc))
    with pytest.raises(kr.RunError, match="roster_sha256"):
        kr.load_catalogs([str(manifest)], 650.0)


def test_a_set_with_a_shared_family_or_record_is_refused(tmp_path):
    a, b = stub_payloads()
    fa = write_stub(tmp_path / "a.json", a)
    twin = write_stub(tmp_path / "twin.json", a)
    with pytest.raises(ValueError, match="toy_homolysis"):
        build_set([fa, twin], tmp_path / "set.json")


# --------------------------------------------------------------------------
# compile --family
# --------------------------------------------------------------------------
@pytest.fixture
def spy(monkeypatch, tmp_path):
    seen = []

    def fake(request):
        seen.append(request)
        return tmp_path / "artifact.json", "artifact"
    monkeypatch.setattr(kcompile, "compile_in_child", fake)
    return seen


def toy_database(tmp_path, families=("R_Recombination", "H_Abstraction")):
    root = tmp_path / "db"
    for family in families:
        (root / "input/kinetics/families" / family).mkdir(parents=True)
        (root / "input/kinetics/families" / family / "groups.py").write_text("# groups\n")
    return root


def test_compile_family_requests_exactly_that_family(spy, tmp_path):
    db = toy_database(tmp_path)
    assert kcompile.main([str(db), str(tmp_path / "out"), "--family", "R_Recombination"]) == 0
    assert [r.families for r in spy] == [("R_Recombination",)]


@pytest.mark.parametrize("name", ["all", "ALL", "none", "default", "!H_Abstraction", "Not_A_Family",
                                  "../db", "families", "R_Recombination,H_Abstraction"])
def test_compile_family_refuses_selectors_and_unknown_names_before_anything_loads(spy, tmp_path, name,
                                                                                   capsys):
    """``all``, a recommended set, ``!F`` and anything that is not one family
    directory of the database would expand into several families (or none)."""
    db = toy_database(tmp_path)
    with pytest.raises(SystemExit) as caught:
        kcompile.main([str(db), str(tmp_path / "out"), "--family", name])
    assert caught.value.code != 0
    assert spy == []
    assert name in capsys.readouterr().err


def test_compile_families_still_reads_the_json_list(spy, tmp_path):
    listing = tmp_path / "families.json"
    listing.write_text(json.dumps(["A", "B"]))
    assert kcompile.main([str(tmp_path), str(tmp_path / "out"), "--families", str(listing)]) == 0
    assert [r.families for r in spy] == [("A", "B")]


@pytest.mark.parametrize("extra", [
    [],                                          # no family at all
    ["--family", "A", "--family", "B"],          # two families
    ["--family", "A", "--families", "@list"],  # both options
])
def test_compile_family_zero_or_two_is_refused_before_anything_loads(spy, tmp_path, extra):
    (tmp_path / "list.json").write_text(json.dumps(["A"]))
    extra = [str(tmp_path / "list.json") if x == "@list" else x for x in extra]
    with pytest.raises(SystemExit) as caught:
        kcompile.main([str(tmp_path), str(tmp_path / "out"), *extra])
    assert caught.value.code != 0
    assert spy == []


# --------------------------------------------------------------------------
# a one-family artifact emits a complete (non-partial) selection
# --------------------------------------------------------------------------
class FakeEmitter:
    families_in_artifact = ["R_Recombination"]

    def __init__(self, artifact, database, systems, families=None, contexts=True):
        self.artifact, self.database, self.systems = Path(artifact), database, systems
        self.only = list(families) if families else None
        self.state = {}

    def init(self):
        records = [{"event_id": f"evt{k}", "family": f} for k, f in
                   enumerate(self.families_in_artifact * 2)]
        self.state.update(artifact_records=len(records),
                          families=sorted({r["family"] for r in records
                                           if not self.only or r["family"] in self.only}))
        self.records = [r for r in records if not self.only or r["family"] in self.only]

    def run_record(self, record):
        return {"event_id": record["event_id"], "family": record["family"], "orientation": "forward",
                "status": "emitted", "seconds": 0.0, "rss_mb": 0}


def selection_of(monkeypatch, tmp_path, artifact_families, **args):
    FakeEmitter.families_in_artifact = artifact_families
    artifact = tmp_path / "artifact.json"
    artifact.write_text("{}")
    monkeypatch.setattr(emit, "Emitter", FakeEmitter)
    monkeypatch.setattr(emit, "load_roster", lambda path: SimpleNamespace(sha256="r" * 64))
    monkeypatch.setattr(emit, "build_reference_systems", lambda roster: ({}, []))
    monkeypatch.setattr(emit, "git_head", lambda path: "d" * 40)
    out = tmp_path / "rows.jsonl"
    database = args.pop("database", "/db")
    ns = SimpleNamespace(roster="r", artifact=str(artifact), database=database, out=str(out), chunk=None,
                         family=None, limit=None, no_contexts=True)
    for k, v in args.items():
        if k != "whole_header":
            setattr(ns, k, v)
    emit.run_emit(ns)
    head = json.loads(out.read_text().splitlines()[0])
    return head if args.get("whole_header") else head["selection"]


def test_a_one_family_artifact_is_a_complete_selection(monkeypatch, tmp_path):
    assert selection_of(monkeypatch, tmp_path, ["R_Recombination"]) == "all"
    assert selection_of(monkeypatch, tmp_path, ["R_Recombination"], family=["R_Recombination"]) == "all"


def test_a_subset_of_families_or_a_limit_stays_partial(monkeypatch, tmp_path):
    assert selection_of(monkeypatch, tmp_path, ["A", "B"], family=["A"]) == "partial"
    assert selection_of(monkeypatch, tmp_path, ["R_Recombination"], limit=1) == "partial"


def test_the_shard_header_carries_content_identities_that_never_need_git(monkeypatch, tmp_path):
    """A database outside git (database_head null) and code outside git
    (code_sha null) are still identified: by the database's content digest
    and the sha256 of the kMC sources."""
    from rmgpy.kmc.database_provenance import database_content_digest
    db = toy_database(tmp_path)
    (db / "input/thermo").mkdir()
    (db / "input/thermo/x.py").write_text("x = 1\n")
    monkeypatch.setattr(emit, "code_sha", lambda: None)
    head = selection_of(monkeypatch, tmp_path, ["R_Recombination"], database=str(db), whole_header=True)
    assert head["database_content_sha256"] == database_content_digest(db)[0]
    assert head["code_sources_sha256"] == emit.code_sources_sha256()
    assert len(head["code_sources_sha256"]) == 64 and head["code_sha"] is None


# --------------------------------------------------------------------------
# rework 1: a set gives what one joint assemble would, or is refused
# --------------------------------------------------------------------------
def colliding_stubs(tmp_path, *, refused_in_a=False):
    """Catalog A gains a record of A's family whose template has B's first
    executable key and weights: one joint assemble would resolve the two
    copies together.  ``refused_in_a``: A's copy was refused inside A (its key
    entered A's collision resolution but A loads nothing under it)."""
    a, b = stub_payloads()
    donor = b[0]["templates"][0]
    copy = json.loads(json.dumps(donor))
    copy["family"] = a[0]["family"]
    host = {"record_event_id": "evt_" + "9" * 64, "family": a[0]["family"], "reverse_of": None,
            "templates": [copy]}
    keys = sorted({t["key"] for p in a for t in p["templates"]} | {donor["key"]})
    if not refused_in_a:
        a.append(host)
    fa = write_stub(tmp_path / "a.json", a, entry_keys=keys)
    fb = write_stub(tmp_path / "b.json", b)
    return fa, fb, donor["key"], host["record_event_id"], b[0]["record_event_id"]


@pytest.mark.parametrize("refused_in_a", [False, True])
def test_a_cross_member_executable_collision_is_refused_by_name_when_the_set_is_built(tmp_path, refused_in_a):
    fa, fb, key, rid_a, rid_b = colliding_stubs(tmp_path, refused_in_a=refused_in_a)
    with pytest.raises(ValueError) as caught:
        build_set([fa, fb], tmp_path / "set.json")
    message = str(caught.value)
    assert key in message and "a.json" in message and "b.json" in message
    assert rid_b in message
    assert not (tmp_path / "set.json").exists()


def test_a_hand_made_set_with_a_cross_member_collision_is_refused_by_the_set_check(tmp_path):
    """The loader repeats the build check: the refusal names the key and the
    members, not the runner's generic duplicate-executable error."""
    fa, fb, key, _ra, _rb = colliding_stubs(tmp_path)
    a, b = stub_payloads()
    ok_a = write_stub(tmp_path / "ok_a.json", a)
    manifest = build_set([ok_a, fb], tmp_path / "set.json")
    doc = json.loads(manifest.read_text())
    import hashlib
    doc["catalogs"][0] = dict(doc["catalogs"][0], path="a.json",
                              sha256=hashlib.sha256(fa.read_bytes()).hexdigest())
    manifest.write_text(json.dumps(doc))
    with pytest.raises(kr.RunError) as caught:
        kr.load_catalogs([str(manifest)], 650.0)
    message = str(caught.value)
    assert key in message and "a.json" in message and "b.json" in message
    assert "fire it twice" not in message


def test_a_member_that_does_not_list_its_entry_keys_is_refused(tmp_path):
    a, b = stub_payloads()
    fa = write_stub(tmp_path / "a.json", a, entry_keys="absent")
    fb = write_stub(tmp_path / "b.json", b)
    with pytest.raises(ValueError, match="entry_keys"):
        build_set([fa, fb], tmp_path / "set.json")


def test_a_set_without_collisions_is_the_joint_catalog(tmp_path):
    """Disjoint entry keys: the set loads what the joint catalog loads."""
    a, b = stub_payloads()
    fa = write_stub(tmp_path / "a.json", a)
    fb = write_stub(tmp_path / "b.json", b)
    both = write_stub(tmp_path / "both.json", a + b)
    manifest = build_set([fa, fb], tmp_path / "set.json")
    for T in (300.0, 650.0, 1000.0):
        assert loaded([manifest], T)[0] == loaded([both], T)[0]


# --------------------------------------------------------------------------
# rework 1: the runner loads the bytes the set verified, once
# --------------------------------------------------------------------------
def test_the_runner_loads_the_bytes_the_set_verified(monkeypatch, tmp_path):
    """A member replaced on disk right after its sha256 was checked: the
    runner still loads the verified content, never the replacement."""
    from rmgpy.kmc import catalog_set as cs
    a, b = stub_payloads()
    fa = write_stub(tmp_path / "a.json", a)
    fb = write_stub(tmp_path / "b.json", b)
    manifest = build_set([fa, fb], tmp_path / "set.json")
    want = loaded([manifest], 650.0)[0]
    tampered = json.loads(json.dumps(b))
    for p in tampered:
        for t in p["templates"]:
            t["weight"]["w"] = [999.0 for _ in t["weight"]["w"]]
    original = cs.read_bytes

    def swap_after_read(path):
        raw = original(path)
        if Path(path).name == "b.json":
            write_stub(tmp_path / "b.json", tampered)
        return raw
    monkeypatch.setattr(cs, "read_bytes", swap_after_read)
    got = loaded([manifest], 650.0)[0]
    assert got == want
    assert 999.0 not in [w for _n, w in got]


# --------------------------------------------------------------------------
# rework 1: unknown provenance is not shared provenance
# --------------------------------------------------------------------------
@pytest.mark.parametrize("field", ["database_path", "database_content_sha256", "code_sources_sha256",
                                   "code_extensions_sha256", "roster_sha256"])
@pytest.mark.parametrize("how", ["null", "missing"])
def test_a_member_with_an_unknown_shared_input_is_refused_by_name(tmp_path, field, how):
    a, b = stub_payloads()
    fa = write_stub(tmp_path / "a.json", a, **({field: None} if how == "null" else {}),
                    drop=(field,) if how == "missing" else ())
    fb = write_stub(tmp_path / "b.json", b, **({field: None} if how == "null" else {}),
                    drop=(field,) if how == "missing" else ())
    with pytest.raises(ValueError) as caught:
        build_set([fa, fb], tmp_path / "set.json")
    assert field in str(caught.value) and "a.json" in str(caught.value)


@pytest.mark.parametrize("field", ["database_head", "code_sha"])
def test_a_git_identity_known_for_one_member_only_is_refused(tmp_path, field):
    """database_head and code_sha are null outside git (the content digests
    then identify the inputs); known for one member and unknown for another
    is not the same input."""
    a, b = stub_payloads()
    fa = write_stub(tmp_path / "a.json", a)
    fb = write_stub(tmp_path / "b.json", b, **{field: None})
    with pytest.raises(ValueError, match=field):
        build_set([fa, fb], tmp_path / "set.json")
    fa = write_stub(tmp_path / "a.json", a, **{field: None})
    manifest = build_set([fa, fb], tmp_path / "set.json")          # unknown for both: the digests decide
    assert json.loads(manifest.read_text())["shared"][field] is None


# --------------------------------------------------------------------------
# rework 2: entry keys must cover the loaded keys; the runner checks the union
# --------------------------------------------------------------------------
def test_a_member_loading_a_key_its_entry_keys_omit_is_refused(tmp_path):
    a, b = stub_payloads()
    keys = sorted({t["key"] for p in a for t in p["templates"]})
    fa = write_stub(tmp_path / "a.json", a, entry_keys=keys[1:])
    fb = write_stub(tmp_path / "b.json", b)
    with pytest.raises(ValueError) as caught:
        build_set([fa, fb], tmp_path / "set.json")
    assert keys[0] in str(caught.value) and "entry_keys" in str(caught.value)


def plain_or_sets(tmp_path, fa, fb, how):
    """The two catalogs as the runner gets them: plain, or one set each."""
    if how == "plain":
        return [fa, fb]
    return [build_set([fa], tmp_path / "set_a.json"), build_set([fb], tmp_path / "set_b.json")]


@pytest.mark.parametrize("how", ["plain", "two_sets"])
def test_catalogs_loaded_together_with_different_inputs_are_refused(tmp_path, how):
    """The set checks hold for everything the runner loads together, not only
    inside one set: two complete catalogs emitted against different databases
    are refused by field name, given directly or through separate sets."""
    a, b = stub_payloads()
    fa = write_stub(tmp_path / "a.json", a)
    fb = write_stub(tmp_path / "b.json", b, database_content_sha256="e" * 64)
    with pytest.raises(kr.RunError, match="database_content_sha256"):
        kr.load_catalogs([str(p) for p in plain_or_sets(tmp_path, fa, fb, how)], 650.0)


@pytest.mark.parametrize("how", ["plain", "two_sets"])
def test_a_historical_key_collision_is_refused_even_with_allow_incomplete(tmp_path, how):
    """A key refused inside incomplete catalog A (entered its collision
    resolution, loads nothing) and loaded by B: --allow-incomplete relaxes
    completeness only, never the collision check."""
    a, b = stub_payloads()
    key = b[0]["templates"][0]["key"]
    keys = sorted({t["key"] for p in a for t in p["templates"]} | {key})
    fa = write_stub(tmp_path / "a.json", a, entry_keys=keys, complete=False)
    fb = write_stub(tmp_path / "b.json", b)
    with pytest.raises(kr.RunError) as caught:
        kr.load_catalogs([str(p) for p in plain_or_sets(tmp_path, fa, fb, how)], 650.0,
                         allow_incomplete=True)
    assert key in str(caught.value) and "collision" in str(caught.value)


def test_one_catalog_alone_needs_no_shared_inputs(tmp_path):
    """Nothing is shared by a catalog loaded alone: an old catalog without the
    content digests or entry keys still runs by itself."""
    a, _b = stub_payloads()
    fa = write_stub(tmp_path / "a.json", a, entry_keys="absent",
                    drop=("database_content_sha256", "code_sources_sha256"))
    assert loaded([fa], 650.0)[0]
