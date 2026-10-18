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
SHARED = {"database_path": "/db", "database_head": "d" * 40, "code_sha": "c" * 40,
          "roster_sha256": "r" * 64}


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


def write_stub(path, payloads, *, complete=True, **shared):
    provenance = dict(SHARED, **shared)
    provenance.update(artifact_id=f"art-{Path(path).stem}", families=sorted({p["family"] for p in payloads}),
                      temperatures=list(GRID))
    doc = {"schema": "kmc-template-catalog/1", "provenance": provenance, "payloads": payloads,
           "refusals": [], "gates": {"failures": [] if complete else ["injected"]},
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


def test_compile_family_requests_exactly_that_family(spy, tmp_path):
    assert kcompile.main([str(tmp_path), str(tmp_path / "out"), "--family", "R_Recombination"]) == 0
    assert [r.families for r in spy] == [("R_Recombination",)]


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
    ns = SimpleNamespace(roster="r", artifact=str(artifact), database="/db", out=str(out), chunk=None,
                         family=None, limit=None, no_contexts=True)
    for k, v in args.items():
        setattr(ns, k, v)
    emit.run_emit(ns)
    return json.loads(out.read_text().splitlines()[0])["selection"]


def test_a_one_family_artifact_is_a_complete_selection(monkeypatch, tmp_path):
    assert selection_of(monkeypatch, tmp_path, ["R_Recombination"]) == "all"
    assert selection_of(monkeypatch, tmp_path, ["R_Recombination"], family=["R_Recombination"]) == "all"


def test_a_subset_of_families_or_a_limit_stays_partial(monkeypatch, tmp_path):
    assert selection_of(monkeypatch, tmp_path, ["A", "B"], family=["A"]) == "partial"
    assert selection_of(monkeypatch, tmp_path, ["R_Recombination"], limit=1) == "partial"
