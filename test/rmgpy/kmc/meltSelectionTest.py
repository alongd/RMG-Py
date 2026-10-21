"""A run restricted to a declared subset of the loaded templates.

``run(..., templates=PATH[, template_classes=...])`` selects loaded templates
from a ``kmc-kernel-rates/1`` table (every template of the named kernel types
or kernel ids, every kernel when none is named) or from a list of template
ids.  The rest stay loaded with their weights but never fire, take no part in
the MET activation, and their hazard is integrated as ``excluded``, never as
leak.  A selection that splits a reciprocal group (``run.reciprocal_groups``)
or names an unknown template or class is refused by name; ``run.json``
records the selection."""

import json
from pathlib import Path

import numpy
import pytest

import melt_engine_support as S
from meltRunTest import lines, pe_deck, toy_payloads, write_catalog, write_deck, PE_LIVE

from rmgpy.kmc import kernel_rates as krt
from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import run as kr

T = 700.0


def paired_catalog(tmp_path, pair=("homolysis", "recombination"), name="paired.json"):
    """The toy PE catalog with the records of ``pair`` declared reciprocal."""
    payloads = toy_payloads([t for t in S.pe_templates() if t.name in PE_LIVE])
    by = {p["record_event_id"]: p for p in payloads}
    a, b = (f"evt_{x}" for x in pair)
    by[a]["reverse_of"], by[b]["reverse_of"] = b, a
    return write_catalog(tmp_path / name, payloads, tag="paired", bind=S.pe_mapping({1: 1}))


def names_of(catalog):
    _cat, records = kr.load_catalogs([catalog], T)
    return [n for r in records for n in r.template_names]


def named(names, *stems):
    return [n for n in names if n.split(":")[0] in {s[:12] for s in stems}]   # names keep 12 characters of the record


def id_list(path, ids):
    Path(path).write_text("".join(f"{x}\n" for x in ids))
    return Path(path)


# --------------------------------------------------------------------------
# reciprocal groups
# --------------------------------------------------------------------------
def _rid(c):
    return "evt_" + c * 12


def _tpl(c, status="emitted", **extra):
    return dict({"key": c * 64, "status": status}, **extra)


def _label(r, k):
    return f"{_rid(r)[4:16]}:{(k * 64)[:12]}"


def test_reciprocal_groups_follow_records_merges_and_covers():
    ref = lambda r, k: {"record": _rid(r), "key": k * 64}
    payloads = [
        {"record_event_id": _rid("a"), "reverse_of": _rid("b"), "templates": [_tpl("1"), _tpl("2")]},
        {"record_event_id": _rid("b"), "reverse_of": _rid("a"), "templates": [_tpl("3", derived={"primary": ref("a", "1")})]},
        {"record_event_id": _rid("c"), "reverse_of": None,
         "templates": [_tpl("4"), _tpl("5", "merged", merged_into=ref("d", "6"))]},
        {"record_event_id": _rid("d"), "reverse_of": None, "templates": [_tpl("6")]},
        {"record_event_id": _rid("e"), "reverse_of": None,
         "templates": [_tpl("7", "covered", covered_by=[ref("f", "8"), ref("g", "9")])]},
        {"record_event_id": _rid("f"), "reverse_of": None, "templates": [_tpl("8")]},
        {"record_event_id": _rid("g"), "reverse_of": None, "templates": [_tpl("9")]},
        {"record_event_id": _rid("h"), "reverse_of": None, "templates": [_tpl("0"), _tpl("x", "refused")]},
    ]
    names = [_label(r, k) for r, k in (("a", "1"), ("a", "2"), ("b", "3"), ("c", "4"), ("d", "6"),
                                       ("f", "8"), ("g", "9"), ("h", "0"))]
    groups = kr.reciprocal_groups(payloads, names)
    assert sorted(groups) == sorted([
        tuple(sorted([_label("a", "1"), _label("a", "2"), _label("b", "3")])),
        tuple(sorted([_label("c", "4"), _label("d", "6")])),
        tuple(sorted([_label("f", "8"), _label("g", "9")])),
        (_label("h", "0"),)])
    # a loaded name the payloads do not give is refused (another label rule)
    with pytest.raises(kr.RunError, match="zzzz"):
        kr.reciprocal_groups(payloads, names + ["zzzz:yyyy"])


# --------------------------------------------------------------------------
# selection
# --------------------------------------------------------------------------
def test_a_split_pair_or_an_unknown_template_is_refused_by_name(tmp_path):
    catalog = paired_catalog(tmp_path)
    names = names_of(catalog)
    deck = pe_deck(tmp_path)
    homolysis, = named(names, "homolysis")
    recombination, = named(names, "recombination")
    with pytest.raises(kr.RunError, match=recombination):
        kr.run(deck, [catalog], T, [0.0, 1.0], [1], tmp_path / "split",
               templates=id_list(tmp_path / "split.txt", [homolysis]))
    with pytest.raises(kr.RunError, match="nonexistent:000000000000"):
        kr.run(deck, [catalog], T, [0.0, 1.0], [1], tmp_path / "unknown",
               templates=id_list(tmp_path / "unknown.txt", [homolysis, recombination, "nonexistent:000000000000"]))
    with pytest.raises(kr.RunError, match="template list"):
        kr.run(deck, [catalog], T, [0.0, 1.0], [1], tmp_path / "classes",
               templates=id_list(tmp_path / "pair.txt", [homolysis, recombination]),
               template_classes=["homolysis"])
    assert not (tmp_path / "split").exists() and not (tmp_path / "unknown").exists()


def test_restricted_run_fires_only_the_selection_and_counts_excluded_hazard_apart(tmp_path):
    """Homolysis + recombination selected on a slow-transport deck: radicals
    build up, so the excluded H transfer, backbiting and beta scission carry
    hazard; it is ``excluded``, and none of it is leak (the selection has no
    branching transition)."""
    catalog = paired_catalog(tmp_path)
    names = names_of(catalog)
    selected = named(names, "homolysis", "recombination")
    ids = id_list(tmp_path / "ids.txt", selected)
    deck = write_deck(tmp_path / "pe.py", S.pe_mapping({10: 2, 20: 2}))
    out = tmp_path / "out"
    manifest = kr.run(deck, [catalog], T, [0.0, 1.0, 3.0], [1, 2], out, templates=ids)
    on_disk = json.loads((out / "run.json").read_text())
    assert on_disk == manifest and manifest["status"] == "complete"
    sel = manifest["template_selection"]
    assert sel["restricted"] is True
    assert sel["templates"] == sorted(selected)
    assert sel["excluded"] == sorted(set(names) - set(selected))
    assert sel["source"]["kind"] == "template_list" and sel["source"]["sha256"] == kr._sha256(ids.read_bytes())
    assert sel["classes"] is None
    fired = set()
    for r in manifest["runs"]:
        fired.update(r["report"]["template_counts"])
    assert fired and fired <= set(selected)
    per_seed = manifest["leak"]["per_seed"]
    assert [x["leak_num"] for x in per_seed] == [0.0, 0.0]
    assert all(x["excluded_num"] > 0.0 for x in per_seed)
    excluded = manifest["leak"]["excluded"]
    assert excluded["templates"] == sel["excluded"]
    assert excluded["excluded_num"] == sum(x["excluded_num"] for x in per_seed)
    assert manifest["leak"]["pooled"]["leak_num"] == 0.0


def test_restricted_engine_is_the_filtered_chemistry(tmp_path):
    """On every state of a restricted trajectory, the propensities and the MET
    pair activations equal those of an engine loading only the selection."""
    catalog_path = paired_catalog(tmp_path)
    catalog, _records = kr.load_catalogs([catalog_path], T)
    names = names_of(catalog_path)
    # every route an excluded hazard takes: unimolecular (beta scission,
    # backbiting), bimolecular non-MET (H abstraction) and MET rows sharing a
    # pair with the kept recombination (disproportionation)
    keep = named(names, "homolysis", "recombination")
    excluded = [n for n in names if n not in keep]

    class Filtered(me.PayloadCatalog):
        def templates(self, temperature):
            return [t for t in super().templates(temperature) if t.name in keep]

    spec = S.pe_spec({10: 2, 20: 2})
    engine = me.MeltEngine(spec, catalog, T, numpy.random.default_rng(3), excluded=excluded)
    filtered = Filtered(catalog.payloads)
    checked = shared = hazard = 0
    for _ in range(40):
        ref = me.MeltEngine(spec, filtered, T, numpy.random.default_rng(0), state=engine.state)
        mine, theirs = engine.propensities(), ref.propensities()
        assert (mine.enabled, mine.leak, mine.thinned, mine.non_thinned) == \
            (theirs.enabled, theirs.leak, theirs.thinned, theirs.non_thinned)
        assert theirs.excluded == 0.0
        a = {k: (v.k_act, v.phi, v.allowed, v.leak) for k, v in engine.met_pair_rates().items() if v.k_act > 0}
        b = {k: (v.k_act, v.phi, v.allowed, v.leak) for k, v in ref.met_pair_rates().items()}
        assert a == b
        checked += bool(a)
        shared += any(v.k_act > 0 and v.excluded > 0 for v in engine.met_pair_rates().values())
        hazard += mine.excluded > 0
        if engine.run(max_events=1).events == 0:
            break
    assert checked and shared and hazard
    assert not set(engine.template_counts) & set(excluded)


def test_excluded_names_must_be_loaded(tmp_path):
    catalog, _records = kr.load_catalogs([paired_catalog(tmp_path)], T)
    with pytest.raises(me.DeclarationError, match="nonexistent"):
        me.MeltEngine(S.pe_spec({10: 1}), catalog, T, numpy.random.default_rng(0), excluded=["nonexistent"])


def test_selection_of_every_template_is_byte_identical_to_the_full_run(tmp_path):
    catalog = paired_catalog(tmp_path)
    deck = write_deck(tmp_path / "pe.py", S.pe_mapping({10: 2, 20: 2}))
    full = kr.run(deck, [catalog], T, [0.0, 1.0, 3.0], [1, 2], tmp_path / "full")
    every = kr.run(deck, [catalog], T, [0.0, 1.0, 3.0], [1, 2], tmp_path / "every",
                   templates=id_list(tmp_path / "all.txt", names_of(catalog)))
    assert full["template_selection"] == {"restricted": False}
    assert every["template_selection"]["restricted"] is True and every["template_selection"]["excluded"] == []
    assert [r["observables_sha256"] for r in every["runs"]] == [r["observables_sha256"] for r in full["runs"]]
    assert all(r["events"] > 0 for r in full["runs"])


@pytest.fixture(scope="module")
def kernel_table(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("kernels")
    catalog = paired_catalog(tmp)
    path = tmp / "toy.kernel_rates.json"
    path.write_text(json.dumps(krt.kernel_rates(json.loads(catalog.read_text()), S.pe_spec({10: 1}), T)))
    return catalog, path


def test_selection_from_kernel_rates_classes(tmp_path, kernel_table):
    catalog, table = kernel_table
    doc = json.loads(table.read_text())
    deck = write_deck(tmp_path / "pe.py", S.pe_mapping({10: 2, 20: 2}))
    manifest = kr.run(deck, [catalog], T, [0.0, 1.0], [1], tmp_path / "out", templates=table,
                      template_classes=["homolysis", "recombination"])
    sel = manifest["template_selection"]
    want = sorted({x["template"] for k in doc["kernels"].values() if k["type"] in ("homolysis", "recombination")
                   for x in k["templates"]})
    assert sel["templates"] == want and len(want) == 2
    assert sel["source"] == {"kind": "kmc-kernel-rates/1", "path": str(table),
                             "sha256": kr._sha256(table.read_bytes()), "temperature_K": T}
    assert sel["classes"] == ["homolysis", "recombination"]
    assert sorted(sel["kernels"]) == sorted(k for k, v in doc["kernels"].items()
                                            if v["type"] in ("homolysis", "recombination"))
    # a kernel id selects that kernel; an unknown class is refused by name
    kid = next(k for k, v in doc["kernels"].items() if v["type"] == "h_transfer")
    m2 = kr.run(deck, [catalog], T, [0.0, 1.0], [1], tmp_path / "kid", templates=table,
                template_classes=["homolysis", "recombination", kid])
    assert set(m2["template_selection"]["templates"]) > set(want)
    with pytest.raises(kr.RunError, match="no_such_kernel"):
        kr.run(deck, [catalog], T, [0.0, 1.0], [1], tmp_path / "bad", templates=table,
               template_classes=["homolysis", "no_such_kernel"])
    # homolysis alone splits its pair with recombination
    with pytest.raises(kr.RunError, match="recombinatio:"):
        kr.run(deck, [catalog], T, [0.0, 1.0], [1], tmp_path / "split", templates=table,
               template_classes=["homolysis"])


def test_a_kernel_table_of_another_chemistry_is_refused(tmp_path, kernel_table):
    catalog, table = kernel_table
    doc = json.loads(table.read_text())
    doc["catalog"]["chemistry_sha256"] = "0" * 64
    other = tmp_path / "other.json"
    other.write_text(json.dumps(doc))
    with pytest.raises(kr.RunError, match="chemistry"):
        kr.run(write_deck(tmp_path / "pe.py", S.pe_mapping({10: 2})), [catalog], T, [0.0, 1.0], [1],
               tmp_path / "out", templates=other)


def test_cli_templates(tmp_path):
    catalog = paired_catalog(tmp_path)
    names = names_of(catalog)
    ids = id_list(tmp_path / "ids.txt", named(names, "homolysis", "recombination"))
    deck = write_deck(tmp_path / "pe.py", S.pe_mapping({10: 1, 20: 1}))
    out = tmp_path / "cli"
    argv = [str(deck), str(catalog), "--temperature", "700", "--end", "2", "--snapshots", "2",
            "--seeds", "3", "--out", str(out), "--templates", str(ids)]
    assert kr.main(argv) == 0
    manifest = json.loads((out / "run.json").read_text())
    assert manifest["template_selection"]["templates"] == sorted(named(names, "homolysis", "recombination"))
    assert lines(out / "seed-3" / "observables.jsonl")
    # classes without a kernel table to read them from
    bad = [str(deck), str(catalog), "--temperature", "700", "--end", "2", "--snapshots", "2",
           "--seeds", "3", "--out", str(tmp_path / "cli2"), "--template-classes", "homolysis"]
    assert kr.main(bad) == 2
