"""Red tests: the melt kMC run from one deck and catalog files (``rmgpy.kmc.run``).

Expected values are derived without the run module: the declaration from
``LinearHomopolymerSpec.from_mapping`` of the mapping the deck was written
from, the initial moments from ``macroscopic_moments_mol``, the pooled leak
from ``pooled_leak_fraction`` over the per-seed parts, the comparison from
``side_by_side`` over the per-seed observation files.

Two polymers: the sampled PS catalog fixtures (``melt_engine_support.PS_CATALOG``)
wrapped as ``kmc-template-catalog/1`` files, and the PE toy (``pe_spec``, the
hand-written toy templates of ``melt_engine_support.pe_templates`` written out
as emitted payloads, constant weights over 300-1000 K).
"""

import ast
import csv
import hashlib
import inspect
import json
import math
import re
import statistics
from pathlib import Path

import pytest

import melt_engine_support as S
from meltStateTest import ps_mapping

from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import melt_observables as mo
from rmgpy.kmc import run as kr
from rmgpy.polymer_spec import LinearHomopolymerSpec

PE_LIVE = ("homolysis", "recombination", "h_abstraction", "disproportionation", "beta_scission",
           "mid_end_recombination", "intra_mid_end")
REL = 1.0e-12
INPUT_PY = Path(kr.__file__).resolve().parents[1] / "rmg" / "input.py"


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


# --------------------------------------------------------------------------
# inputs: decks and catalog files
# --------------------------------------------------------------------------
def ps_deck_mapping(counts):
    """The mapping ``melt_engine_support.ps_spec(counts)`` declares (no styrene)."""
    mapping = ps_mapping(counts, mobile={})
    mapping["small_species"] = {}
    mapping["transport_model"] = {"met": S.ps_met()}
    return mapping


DECK_HEAD = """\
database(
    thermoLibraries=['primaryThermoLibrary'],
    reactionLibraries=[],
    seedMechanisms=[],
    kineticsDepositories='default',
    kineticsFamilies='default',
    kineticsEstimator='rate rules',
)
species(label='N2', reactive=False, structure=SMILES('N#N'))
"""
DECK_TAIL = """\
simpleReactor(
    temperature=(700, 'K'),
    pressure=(1.0, 'bar'),
    initialMoleFractions={'N2': 1.0},
    terminationTime=(1.0, 's'),
)
model(toleranceKeepInEdge=0.0, toleranceMoveToCore=0.1, toleranceInterruptSimulation=0.1)
options(units='si', generateOutputHTML=False)
"""


def write_deck(path, mapping, *, extra="", cutoff=None, via_variable=False, sgh=None):
    args = [repr(mapping["label"])]
    if via_variable:
        body = f"spec = {mapping!r}\n"
        args.append("linear_homopolymer_spec=spec")
    else:
        body = ""
        args.append(f"linear_homopolymer_spec={mapping!r}")
    if cutoff is not None:
        args.append(f"cutoff={cutoff}")
    if sgh is not None:
        args.append(f"side_group_homolysis={sgh!r}")
    text = DECK_HEAD + body + f"polymer({', '.join(args)})\n" + extra + DECK_TAIL
    Path(path).write_text(text)
    return Path(path)


def toy_payloads(templates, grid=(300.0, 1000.0)):
    """Toy ``window_matcher.Template`` objects as emitted payloads (constant weight)."""
    out = []
    for t in templates:
        final_bond = {}
        final_atom = {}
        for edit in t.edits:
            if hasattr(edit, "order"):
                final_bond[(min(edit.a, edit.b), max(edit.a, edit.b))] = edit.order
            else:
                final_atom[edit.atom] = {k: getattr(edit, k) for k in ("radical", "charge", "lone_pairs")
                                         if getattr(edit, k) is not None}
        keys = sorted(set(t.bond) | set(final_bond))
        bonds = [[i, j, t.bond.get((i, j), 0), final_bond.get((i, j), t.bond.get((i, j), 0))] for i, j in keys]
        atoms = [{"element": a.element, "participant": a.participant, "kind": a.kind,
                  "exterior": [list(x) for x in a.exterior], "isotope": a.isotope, "charge": a.charge,
                  "radical": a.radical, "lone_pairs": a.lone_pairs, "final": final_atom.get(i, {})}
                 for i, a in enumerate(t.atoms)]
        key = hashlib.sha256(t.name.encode()).hexdigest()
        out.append({"record_event_id": f"evt_{t.name}", "family": t.family, "reverse_of": None,
                    "templates": [{"key": key, "name": t.name, "family": t.family, "status": "emitted",
                                   "direction": "forward",
                                   "executable": {"atoms": atoms, "bonds": bonds},
                                   "weight": {"T": list(grid), "w": [t.weight] * len(grid)}}]})
    return out


# the chemistry a catalog is bound to: these fields of the normalized declaration
# the declaration fields that do NOT enter the chemistry digest (everything else does)
CHEMISTRY_EXCLUDED = ("label", "initial_ensemble", "initial_mobile_counts", "initial_mass_kg", "Mn_g_mol",
                      "Mw_g_mol", "density_kg_m3", "density_basis", "transport_model",
                      "observable_relative_tolerance", "diagnostics")


def chemistry_of_normalized(normalized):
    return hashlib.sha256(json.dumps({k: v for k, v in normalized.items() if k not in CHEMISTRY_EXCLUDED},
                                     sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def chemistry_sha256(mapping):
    return chemistry_of_normalized(json.loads(LinearHomopolymerSpec.from_mapping(mapping).normalized_json))


def declaration_sha256(mapping):
    return hashlib.sha256(LinearHomopolymerSpec.from_mapping(mapping).normalized_json).hexdigest()


def write_catalog(path, payloads, *, complete=True, failures=(), tag="test", bind=None, declaration=None):
    """``bind``: a declaration mapping whose chemistry digest the provenance
    carries; ``declaration``: one whose full declaration digest it carries."""
    failures = list(failures)
    provenance = {"artifact_id": f"artifact-{tag}", "families": sorted({p["family"] for p in payloads}),
                  "code_sha": None, "temperatures": [300.0, 1000.0]}
    if bind is not None:
        provenance["chemistry_sha256"] = chemistry_sha256(bind)
    if declaration is not None:
        provenance["declaration_sha256"] = declaration_sha256(declaration)
    doc = {"schema": "kmc-template-catalog/1",
           "provenance": provenance,
           "payloads": payloads, "refusals": [],
           "gates": {"failures": failures}, "completeness": {}, "complete": complete}
    Path(path).write_text(json.dumps(doc, sort_keys=True, indent=1) + "\n")
    return Path(path)


def pe_catalog_file(tmp_path, name="pe_toy_catalog.json", **kw):
    kw.setdefault("bind", S.pe_mapping({1: 1}))
    return write_catalog(tmp_path / name, toy_payloads(
        [t for t in S.pe_templates() if t.name in PE_LIVE]), tag="pe-toy", **kw)


def ps_catalog_files(tmp_path, bind=True):
    return [write_catalog(tmp_path / f"ps_{i}.json", S.payloads(name), tag=f"ps-{i}",
                          bind=ps_deck_mapping({1: 1}) if bind else None)
            for i, name in enumerate(S.PS_CATALOG)]


def pe_deck(tmp_path, counts=None, **kw):
    return write_deck(tmp_path / "pe_deck.py", S.pe_mapping(counts or {10: 1, 20: 1}), **kw)


def ps_deck(tmp_path, counts=None, **kw):
    return write_deck(tmp_path / "ps_deck.py", ps_deck_mapping(counts or {6: 3}), **kw)


def lines(path):
    """The observation lines of a time series (typed records, e.g. the stop record, skipped)."""
    return [r for r in (json.loads(x) for x in Path(path).read_text().splitlines()) if "record" not in r]


def records(path):
    return [r for r in (json.loads(x) for x in Path(path).read_text().splitlines()) if "record" in r]


def strip_wall(obj):
    if isinstance(obj, dict):
        return {k: strip_wall(v) for k, v in obj.items() if k != "wall_s"}
    if isinstance(obj, list):
        return [strip_wall(v) for v in obj]
    return obj


# --------------------------------------------------------------------------
# deck parsing
# --------------------------------------------------------------------------
@pytest.mark.parametrize("polymer", ["ps", "pe"])
def test_deck_declaration_equals_from_mapping(tmp_path, polymer):
    mapping = ps_deck_mapping({6: 3}) if polymer == "ps" else S.pe_mapping({10: 1, 20: 1})
    path = write_deck(tmp_path / "deck.py", mapping, cutoff=3)
    deck = kr.read_deck(path)
    want = LinearHomopolymerSpec.from_mapping(mapping)
    assert deck.spec == want
    assert deck.spec.normalized_json == want.normalized_json
    assert deck.xs == 3 and deck.label == mapping["label"]
    assert deck.sha256 == hashlib.sha256(path.read_bytes()).hexdigest()
    assert deck.side_group_homolysis is None


def test_deck_without_cutoff_uses_the_polymer_default(tmp_path):
    import rmgpy.rmg.input as rmg_input
    default = inspect.signature(rmg_input.polymer).parameters["cutoff"].default
    deck = kr.read_deck(pe_deck(tmp_path))
    assert deck.xs == default
    assert kr.read_deck(pe_deck(tmp_path, cutoff=2)).xs == 2


def test_deck_reader_binds_polymer_arguments_in_the_input_order():
    """Positional polymer() arguments bind in the input reader's order."""
    import rmgpy.rmg.input as rmg_input
    assert tuple(inspect.signature(rmg_input.polymer).parameters) == kr.POLYMER_PARAMETERS


def test_deck_with_no_or_several_declarations_is_refused(tmp_path):
    mapping = S.pe_mapping({10: 1})
    none = tmp_path / "none.py"
    none.write_text(DECK_HEAD + DECK_TAIL)
    with pytest.raises(kr.RunError, match="no polymer"):
        kr.read_deck(none)
    legacy = tmp_path / "legacy.py"
    legacy.write_text(DECK_HEAD + "polymer('pe', monomer='[*:1]CC[*:2]', end_groups=['[H][*:1]', '[H][*:2]'], "
                      "Mn=1000.0, Mw=2000.0)\n" + DECK_TAIL)
    with pytest.raises(kr.RunError, match="linear_homopolymer_spec"):
        kr.read_deck(legacy)
    second = f"polymer({mapping['label']!r}, linear_homopolymer_spec={mapping!r}, cutoff=2)\n"
    with pytest.raises(kr.RunError, match="2 polymer"):
        kr.read_deck(write_deck(tmp_path / "two.py", mapping, extra=second))
    other = S.pe_mapping({30: 1})
    other["label"] = "second pool"
    second = f"polymer('second pool', linear_homopolymer_spec={other!r})\n"
    with pytest.raises(kr.RunError, match="2 polymer"):
        kr.read_deck(write_deck(tmp_path / "two_labels.py", mapping, extra=second))
    second = "polymer('legacy', monomer='[*:1]CC[*:2]', end_groups=['[H][*:1]', '[H][*:2]'], Mn=1000.0, Mw=2000.0)\n"
    with pytest.raises(kr.RunError, match="2 polymer"):
        kr.read_deck(write_deck(tmp_path / "spec_and_legacy.py", mapping, extra=second))
    mislabelled = tmp_path / "mislabelled.py"
    mislabelled.write_text(DECK_HEAD + f"polymer('other', linear_homopolymer_spec={mapping!r})\n" + DECK_TAIL)
    with pytest.raises(kr.RunError, match="label"):
        kr.read_deck(mislabelled)
    broken = tmp_path / "broken.py"
    broken.write_text(DECK_HEAD + "polymer('pe', linear_homopolymer_spec=undefined_name)\n")
    with pytest.raises(kr.RunError):
        kr.read_deck(broken)


# --------------------------------------------------------------------------
# catalogs
# --------------------------------------------------------------------------
def test_catalog_schema_matches_the_emitter():
    from rmgpy.kmc import emit
    assert kr.CATALOG_SCHEMA == emit.SCHEMA


def test_incomplete_catalog_is_refused_and_labelled_with_the_flag(tmp_path):
    deck = pe_deck(tmp_path)
    failures = ["cycles: 2 of 7 pairs refused", "midpoint: 1 template above 0.01"]
    catalog = pe_catalog_file(tmp_path, complete=False, failures=failures)
    with pytest.raises(kr.RunError, match="incomplete"):
        kr.run(deck, [catalog], 700.0, [0.0, 1.0], [1], tmp_path / "refused")
    assert not (tmp_path / "refused" / "run.json").exists()
    out = tmp_path / "accepted"
    manifest = kr.run(deck, [catalog], 700.0, [0.0, 1.0], [1, 2], out, allow_incomplete=True)
    on_disk = json.loads((out / "run.json").read_text())
    for m in (manifest, on_disk):
        assert m["catalog_complete"] is False
        assert [f["failures"] for f in m["catalog_gate_failures"]] == [failures]
        assert m["catalogs"][0]["complete"] is False
    for seed in (1, 2):
        rows = lines(out / f"seed-{seed}" / "observables.jsonl")
        assert len(rows) == 2
        for row in rows:
            assert row["catalog_complete"] is False
            assert [f["failures"] for f in row["catalog_gate_failures"]] == [failures]
    head = (out / "comparison.csv").read_text().splitlines()[0]
    assert head.startswith("#") and "catalog_complete=false" in head and failures[0] in head
    # a complete catalog is labelled complete
    good = kr.run(deck, [pe_catalog_file(tmp_path, "complete.json")], 700.0, [0.0], [1], tmp_path / "good")
    assert good["catalog_complete"] is True and good["catalog_gate_failures"] == []
    assert lines(tmp_path / "good" / "seed-1" / "observables.jsonl")[0]["catalog_complete"] is True


def test_catalog_schema_flag_and_collisions_are_checked(tmp_path):
    deck = pe_deck(tmp_path)
    catalog = pe_catalog_file(tmp_path)
    doc = json.loads(catalog.read_text())
    wrong = tmp_path / "wrong_schema.json"
    wrong.write_text(json.dumps(dict(doc, schema="kmc-template-catalog/0")))
    with pytest.raises(kr.RunError, match="schema"):
        kr.load_catalogs([wrong], 700.0)
    flagless = tmp_path / "flagless.json"
    flagless.write_text(json.dumps({k: v for k, v in doc.items() if k != "complete"}))
    with pytest.raises(kr.RunError, match="complete"):
        kr.load_catalogs([flagless], 700.0)
    inconsistent = tmp_path / "inconsistent.json"
    inconsistent.write_text(json.dumps(dict(doc, gates={"failures": ["x"]})))
    with pytest.raises(kr.RunError, match="complete"):
        kr.load_catalogs([inconsistent], 700.0)
    copy = tmp_path / "copy.json"
    copy.write_text(catalog.read_text())
    with pytest.raises(kr.RunError, match="collide"):
        kr.load_catalogs([catalog, copy], 700.0)
    with pytest.raises(kr.RunError, match="collide"):
        kr.run(deck, [catalog, copy], 700.0, [0.0], [1], tmp_path / "collide")
    # two disjoint catalogs load together and each provenance is recorded
    loaded, records = kr.load_catalogs(ps_catalog_files(tmp_path), 700.0)
    assert [r.provenance["artifact_id"] for r in records] == ["artifact-ps-0", "artifact-ps-1"]
    names = [t.name for t in loaded.templates(700.0)]
    assert sorted(names) == sorted(t.name for t in S.ps_catalog().templates(700.0))


# --------------------------------------------------------------------------
# the run
# --------------------------------------------------------------------------
def test_t0_observation_equals_the_declared_macroscopic_moments(tmp_path):
    counts = {6: 3, 9: 1}
    deck = ps_deck(tmp_path, counts=counts, cutoff=3)
    out = tmp_path / "ps"
    kr.run(deck, ps_catalog_files(tmp_path), 700.0, [0.0, 2.0e6], [3], out)
    spec = LinearHomopolymerSpec.from_mapping(ps_deck_mapping(counts))
    first = lines(out / "seed-3" / "observables.jsonl")[0]
    assert first["time_s"] == 0.0 and first["xs"] == 3
    for got, want in zip(first["tier2"]["dp_moments_mol"], spec.macroscopic_moments_mol):
        assert rel(got, want) <= REL, (got, want)
    assert first["tier2"]["dp_moments_representative"] == [4, 27, 3 * 36 + 81]


def test_same_seed_is_byte_identical_and_a_different_seed_differs(tmp_path):
    deck, catalog = pe_deck(tmp_path), pe_catalog_file(tmp_path)
    times = [0.0, 1.0, 3.0]
    a = kr.run(deck, [catalog], 700.0, times, [1, 2], tmp_path / "a")
    b = kr.run(deck, [catalog], 700.0, times, [2], tmp_path / "b")
    c = kr.run(deck, [catalog], 700.0, times, [1, 2], tmp_path / "c")
    seed2 = (tmp_path / "a" / "seed-2" / "observables.jsonl").read_bytes()
    assert seed2 == (tmp_path / "b" / "seed-2" / "observables.jsonl").read_bytes()
    assert seed2 == (tmp_path / "c" / "seed-2" / "observables.jsonl").read_bytes()
    seed1 = (tmp_path / "a" / "seed-1" / "observables.jsonl").read_bytes()
    assert seed1 == (tmp_path / "c" / "seed-1" / "observables.jsonl").read_bytes()
    assert seed1 != seed2
    assert lines(tmp_path / "a" / "seed-1" / "observables.jsonl")[-1]["counts"]["events"] > 0
    # the manifest differs only in wall times
    on_a = json.loads((tmp_path / "a" / "run.json").read_text())
    on_c = json.loads((tmp_path / "c" / "run.json").read_text())
    assert strip_wall(on_a) == strip_wall(on_c) == strip_wall(a) == strip_wall(c)
    assert on_a["runs"][1]["observables_sha256"] == b["runs"][0]["observables_sha256"]
    assert on_a["runs"][1]["observables_sha256"] == hashlib.sha256(seed2).hexdigest()
    with pytest.raises(kr.RunError, match="seed"):
        kr.run(deck, [catalog], 700.0, times, [1, 1], tmp_path / "dup")


def test_manifest_leak_is_the_pooled_ensemble_estimate(tmp_path):
    deck, catalog = pe_deck(tmp_path, counts={10: 2, 20: 2}), pe_catalog_file(tmp_path)
    out = tmp_path / "leak"
    manifest = kr.run(deck, [catalog], 700.0, [0.0, 4.0], [1, 2, 3, 4], out)
    leak = manifest["leak"]
    per_seed = leak["per_seed"]
    assert [p["seed"] for p in per_seed] == [1, 2, 3, 4]
    for p in per_seed:
        last = lines(out / f"seed-{p['seed']}" / "observables.jsonl")[-1]["leak_estimate"]
        for part in ("leak_num", "leak_den", "non_thinned_integral", "phi_sum", "thinned_proposals"):
            assert p[part] == last[part]
        assert "fraction" not in p
    want = me.pooled_leak_fraction(per_seed)
    assert leak["pooled"] == want
    assert sum(p["leak_num"] for p in per_seed) > 0.0
    fractions = [p["leak_num"] / p["leak_den"] for p in per_seed]
    assert rel(statistics.fmean(fractions), want["fraction"]) > 1e-9   # the test tells the two apart
    assert "leak_fraction" not in json.dumps(manifest["runs"])
    for r in manifest["runs"]:
        assert "fraction" not in json.dumps(r)


def test_event_cap_stops_the_seed_and_records_the_stop(tmp_path):
    deck, catalog = pe_deck(tmp_path), pe_catalog_file(tmp_path)
    out = tmp_path / "cap"
    manifest = kr.run(deck, [catalog], 700.0, [0.0, 1.0, 3.0, 50.0], [1], out, max_events=5)
    run = manifest["runs"][0]
    assert run["stop_reason"] == "max_events" and run["events"] == 5
    rows = lines(out / "seed-1" / "observables.jsonl")
    assert run["observations"] == len(rows) < 4
    assert all(r["counts"]["events"] <= 5 for r in rows)
    assert run["engine_time_s"] < 50.0
    free = kr.run(deck, [catalog], 700.0, [0.0, 1.0], [1], tmp_path / "free")
    assert free["runs"][0]["stop_reason"] in ("horizon", "absorbing")


def test_jobs_gives_the_sequential_result(tmp_path, monkeypatch):
    deck, catalog = pe_deck(tmp_path), pe_catalog_file(tmp_path)
    times = [0.0, 2.0]
    seq = kr.run(deck, [catalog], 700.0, times, [1, 2, 3], tmp_path / "seq")

    class Inline:
        """An executor that runs each task at submission (no extra process)."""

        def __init__(self, jobs):
            self.jobs = jobs

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def map(self, fn, *iterables):
            return [fn(*args) for args in zip(*iterables)]

    monkeypatch.setattr(kr, "_executor", lambda jobs, context: Inline(jobs))
    par = kr.run(deck, [catalog], 700.0, times, [1, 2, 3], tmp_path / "par", jobs=2)
    for seed in (1, 2, 3):
        assert ((tmp_path / "seq" / f"seed-{seed}" / "observables.jsonl").read_bytes()
                == (tmp_path / "par" / f"seed-{seed}" / "observables.jsonl").read_bytes())
    assert par["jobs"] == 2 and seq["jobs"] == 1
    assert strip_wall(dict(par, jobs=1)) == strip_wall(seq)


def test_times_from_end_and_snapshot_count():
    assert kr.observation_times(end=10.0, snapshots=3) == [0.0, 5.0, 10.0]
    assert kr.observation_times(times=[0.0, 1.0]) == [0.0, 1.0]
    for bad in (dict(times=[1.0, 0.5]), dict(times=[-1.0]), dict(end=1.0, snapshots=1),
                dict(times=[0.0], end=1.0, snapshots=2), dict(), dict(times=[math.inf])):
        with pytest.raises(kr.RunError):
            kr.observation_times(**bad)


# --------------------------------------------------------------------------
# CLI end to end
# --------------------------------------------------------------------------
@pytest.mark.parametrize("polymer", ["ps", "pe"])
def test_cli_end_to_end(tmp_path, polymer):
    if polymer == "ps":
        deck, catalogs, end = ps_deck(tmp_path, cutoff=3), ps_catalog_files(tmp_path), 1.0e7
    else:
        deck, catalogs, end = pe_deck(tmp_path, cutoff=3), [pe_catalog_file(tmp_path)], 3.0
    out = tmp_path / "out"
    argv = [str(deck), *map(str, catalogs), "--temperature", "700", "--end", str(end), "--snapshots", "3",
            "--seeds", "5", "6", "--out", str(out)]
    assert kr.main(argv) == 0
    manifest = json.loads((out / "run.json").read_text())
    assert manifest["schema"] == kr.RUN_SCHEMA and manifest["status"] == "complete"
    assert manifest["times_s"] == [0.0, end / 2, end]
    assert manifest["temperature_K"] == 700.0 and manifest["seeds"] == [5, 6]
    assert manifest["deck"]["sha256"] == hashlib.sha256(deck.read_bytes()).hexdigest()
    assert manifest["deck"]["declaration"] == json.loads(
        LinearHomopolymerSpec.from_mapping(
            ps_deck_mapping({6: 3}) if polymer == "ps" else S.pe_mapping({10: 1, 20: 1})).normalized_json)
    assert len(manifest["catalogs"]) == len(catalogs)
    for record, path in zip(manifest["catalogs"], catalogs):
        assert record["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
        assert record["provenance"] == json.loads(path.read_text())["provenance"]
    code = manifest["code"]
    assert re.fullmatch(r"[0-9a-f]{40}", code["sha"]) and isinstance(code["dirty"], bool)
    for r in manifest["runs"]:
        assert r["stop_reason"] in ("horizon", "absorbing")
        assert r["wall_s"] >= 0.0 and r["observations"] == 3
        rows = lines(out / f"seed-{r['seed']}" / "observables.jsonl")
        assert [x["time_s"] for x in rows] == manifest["times_s"]
        assert all(x["schema"] == mo.OBSERVABLES_SCHEMA and x["invariant_checks"] is not None for x in rows)
        assert rows[-1]["invariant_checks"] == rows[-1]["counts"]["events"]
        assert "declared_mobile_participation" in r["report"]
    assert (out / "comparison.csv").exists()
    # a non-empty output directory is refused
    assert kr.main(argv) == 2


def test_run_module_names_no_polymer():
    source = Path(kr.__file__).read_text()
    code = re.sub(r'"""[\s\S]*?"""', "", source)
    code = re.sub(r"#.*", "", code)
    assert not re.search(r"\b(PS|PE|PP|styrene|polystyrene|polyethylene|phenyl|benzylic|aryl)\b", code, re.I)
    for name in ("PS_M0", "PS_C_R2", "rmgpy.kmc.met", "pe_", "ps_"):
        assert name not in code


# --------------------------------------------------------------------------
# side by side with the MoM solver
# --------------------------------------------------------------------------
def read_comparison(path):
    text = Path(path).read_text().splitlines()
    assert text[0].startswith("#")
    return list(csv.DictReader(text[1:]))


def test_mom_side_by_side_gives_the_expected_differences(tmp_path):
    deck, catalog = pe_deck(tmp_path, cutoff=3), pe_catalog_file(tmp_path)
    times = [0.0, 2.0, 5.0]
    seeds = [1, 2, 3]
    base = tmp_path / "base"
    kr.run(deck, [catalog], 700.0, times, seeds, base)
    obs = {s: lines(base / f"seed-{s}" / "observables.jsonl") for s in seeds}
    # synthetic MoM numbers: the seed-1 values, shifted
    entries = []
    for k, t in enumerate(times[:2]):
        q = mo.mom_reported_quantities(obs[1][k])
        q["mu0_mol"] = q["mu0_mol"] * 1.25
        q["polymer_condensed_mass_g"] = q["polymer_condensed_mass_g"] - 0.5
        q["explicit_dp_mol"] = dict(q["explicit_dp_mol"], **{"2": 7.0e-6})
        entries.append({"time_s": t, "quantities": q})
    mom = tmp_path / "mom.json"
    mom.write_text(json.dumps(entries))
    out = tmp_path / "cmp"
    kr.run(deck, [catalog], 700.0, times, seeds, out, mom=mom)
    for s in seeds:
        assert ((out / f"seed-{s}" / "observables.jsonl").read_bytes()
                == (base / f"seed-{s}" / "observables.jsonl").read_bytes())
    rows = read_comparison(out / "comparison.csv")
    by = {(float(r["time_s"]), r["quantity"]): r for r in rows}
    for k, t in enumerate(times):
        mom_q = entries[k]["quantities"] if k < len(entries) else None
        per_seed = {}
        for s in seeds:
            for r in mo.side_by_side(obs[s][k], mom_q):
                per_seed.setdefault(r["quantity"], []).append(r)
        checked = ("mu0_mol", "polymer_condensed_mass_g", "Mn_g_mol")
        for quantity in checked + (("explicit_dp_mol.2",) if mom_q is not None else ()):
            row = by[(t, quantity)]
            values = [r["kmc"] for r in per_seed[quantity]]
            values = [0.0 if v is None else v for v in values]
            mean = math.fsum(values) / len(values)
            assert int(row["seeds"]) == len(seeds)
            assert rel(float(row["kmc_mean"]), mean) <= REL or mean == float(row["kmc_mean"])
            se = statistics.stdev(values) / math.sqrt(len(values))
            assert abs(float(row["kmc_stderr"]) - se) <= 1e-12 * max(abs(se), abs(mean), 1e-300)
            if mom_q is None:
                assert row["mom"] == "" and row["difference"] == ""
                continue
            want = per_seed[quantity][0]["mom"]
            assert float(row["mom"]) == want
            assert rel(float(row["difference"]), mean - want) <= 1e-9
    # t = 0: every seed is the declaration, so the differences are the shifts exactly
    q0 = mo.mom_reported_quantities(obs[1][0])
    assert rel(float(by[(0.0, "mu0_mol")]["difference"]), q0["mu0_mol"] - q0["mu0_mol"] * 1.25) <= 1e-12
    assert rel(float(by[(0.0, "polymer_condensed_mass_g")]["difference"]), 0.5) <= 1e-12
    assert float(by[(0.0, "mu0_mol")]["kmc_stderr"]) == 0.0
    assert by[(0.0, "convention_gap_g")]["mom"] == ""
    manifest = json.loads((out / "run.json").read_text())
    assert manifest["mom"]["sha256"] == hashlib.sha256(mom.read_bytes()).hexdigest()
    # MoM times the run did not observe are refused
    bad = tmp_path / "bad_mom.json"
    bad.write_text(json.dumps([{"time_s": 1.5, "quantities": entries[0]["quantities"]}]))
    with pytest.raises(kr.RunError, match="time"):
        kr.run(deck, [catalog], 700.0, times, seeds, tmp_path / "bad", mom=bad)


def test_deck_side_group_channels_reach_the_observables(tmp_path):
    """The deck's side_group_homolysis (the MoM channels) partitions the
    observables; a channel the declaration contradicts is refused before any output."""
    mapping = S.para_bromo_ps_spec({6: 2}).as_mapping()
    mapping["transport_model"] = {"met": S.ps_met()}
    channels = [S.sgh_channel("aryl_C-Br", "aryl", 1.0)]
    deck = write_deck(tmp_path / "deck.py", mapping, cutoff=3, sgh=channels)
    assert kr.read_deck(deck).side_group_homolysis == channels
    out = tmp_path / "sgh"
    manifest = kr.run(deck, ps_catalog_files(tmp_path, bind=False), 700.0, [0.0, 1.0e6], [1], out,
                      allow_unbound_catalog=True, allow_unmatched_catalog=True)
    assert manifest["deck"]["side_group_homolysis"] == channels
    first = lines(out / "seed-1" / "observables.jsonl")[0]
    (row,) = first["side_group_channels"]
    assert row["label"] == "aryl_C-Br" and row["released_tracked"] is True
    unit = first["scale_factor"] / first["avogadro"]
    assert rel(row["inventory_tail_mol"] + row["inventory_explicit_mol"], 12 * unit) <= REL
    wrong = write_deck(tmp_path / "wrong.py", mapping, sgh=[S.sgh_channel("aryl_C-Br", "aryl", 2.0)])
    with pytest.raises(kr.RunError):
        kr.run(wrong, ps_catalog_files(tmp_path, bind=False), 700.0, [0.0], [1], tmp_path / "wrong_out",
               allow_unbound_catalog=True, allow_unmatched_catalog=True)
    assert not (tmp_path / "wrong_out").exists()


# --------------------------------------------------------------------------
# rework 1 (spar r108)
# --------------------------------------------------------------------------
def renamed(payloads, prefix):
    out = json.loads(json.dumps(payloads))
    for p in out:
        p["record_event_id"] = "evt_" + prefix + p["record_event_id"][4:]
    return out


def test_catalogs_sharing_an_executable_are_refused(tmp_path):
    """Item 1: the same executable under another record id (so another template
    name) in a second catalog would double its rate."""
    deck = pe_deck(tmp_path, counts={5: 1})
    payloads = toy_payloads([t for t in S.pe_templates() if t.name in PE_LIVE])
    first = write_catalog(tmp_path / "first.json", payloads, bind=S.pe_mapping({1: 1}))
    homolysis = [p for p in renamed(payloads, "x") if p["templates"][0]["name"] == "homolysis"]
    second = write_catalog(tmp_path / "second.json", homolysis, bind=S.pe_mapping({1: 1}))
    names_first = {t.name for t in me.PayloadCatalog(payloads).templates(700.0)}
    names_second = {t.name for t in me.PayloadCatalog(homolysis).templates(700.0)}
    assert not names_first & names_second            # the name check alone lets them through
    with pytest.raises(kr.RunError, match="executable"):
        kr.load_catalogs([first, second], 700.0)
    with pytest.raises(kr.RunError, match="executable"):
        kr.run(deck, [first, second], 700.0, [0.0], [1], tmp_path / "doubled")
    assert not (tmp_path / "doubled").exists()
    within = write_catalog(tmp_path / "within.json", payloads + homolysis, bind=S.pe_mapping({1: 1}))
    with pytest.raises(kr.RunError, match="executable"):
        kr.load_catalogs([within], 700.0)


def test_catalog_bound_to_another_declaration_is_refused(tmp_path):
    """Item 2a: a catalog whose provenance digest names another declaration."""
    pe_map = S.pe_mapping({10: 1, 20: 1})
    pe = pe_deck(tmp_path)
    with pytest.raises(kr.RunError, match="chemistry_sha256"):
        kr.run(pe, ps_catalog_files(tmp_path), 700.0, [0.0], [1], tmp_path / "a",
               allow_unbound_catalog=True, allow_unmatched_catalog=True)
    wrong = pe_catalog_file(tmp_path, "wrong_declaration.json", bind=None, declaration=S.pe_mapping({30: 1}))
    with pytest.raises(kr.RunError, match="declaration_sha256"):
        kr.run(pe, [wrong], 700.0, [0.0], [1], tmp_path / "b", allow_unbound_catalog=True)
    exact = pe_catalog_file(tmp_path, "exact.json", bind=None, declaration=pe_map)
    m = kr.run(pe, [exact], 700.0, [0.0], [1], tmp_path / "c")
    assert m["catalog_bound"] is True and m["catalogs"][0]["bound_by"] == ["declaration_sha256"]
    # the digests the manifest documents are the ones a catalog carries
    assert m["deck"]["chemistry_sha256"] == chemistry_sha256(pe_map)
    assert m["deck"]["declaration_sha256"] == declaration_sha256(pe_map)
    assert kr.CHEMISTRY_EXCLUDED == CHEMISTRY_EXCLUDED
    bound = kr.run(pe, [pe_catalog_file(tmp_path)], 700.0, [0.0], [1], tmp_path / "d")
    assert bound["catalogs"][0]["bound_by"] == ["chemistry_sha256"]


def test_unbound_catalog_needs_the_flag_and_is_labelled(tmp_path):
    """Item 2b."""
    deck = pe_deck(tmp_path)
    unbound = pe_catalog_file(tmp_path, "unbound.json", bind=None)
    with pytest.raises(kr.RunError, match="allow-unbound-catalog"):
        kr.run(deck, [unbound], 700.0, [0.0, 1.0], [1], tmp_path / "refused")
    assert not (tmp_path / "refused").exists()
    out = tmp_path / "unbound"
    m = kr.run(deck, [unbound], 700.0, [0.0, 1.0], [1, 2], out, allow_unbound_catalog=True)
    assert m["catalog_bound"] is False and m["catalogs"][0]["bound_by"] == []
    for seed in (1, 2):
        assert all(r["catalog_bound"] is False for r in lines(out / f"seed-{seed}" / "observables.jsonl"))
    assert "catalog_bound=false" in (out / "comparison.csv").read_text().splitlines()[0]
    good = tmp_path / "bound"
    m = kr.run(deck, [pe_catalog_file(tmp_path)], 700.0, [0.0], [1], good)
    assert m["catalog_bound"] is True
    assert lines(good / "seed-1" / "observables.jsonl")[0]["catalog_bound"] is True
    assert "catalog_bound=true" in (good / "comparison.csv").read_text().splitlines()[0]
    argv = [str(deck), str(unbound), "--temperature", "700", "--times", "0", "--seeds", "1"]
    assert kr.main(argv + ["--out", str(tmp_path / "cli_refused")]) == 2
    assert kr.main(argv + ["--out", str(tmp_path / "cli"), "--allow-unbound-catalog"]) == 0


def test_catalog_matching_no_initial_site_is_refused_or_labelled(tmp_path):
    """Item 2c: the wrong polymer's catalog matches nothing in the initial melt."""
    deck = pe_deck(tmp_path, counts={5: 1})
    ps = ps_catalog_files(tmp_path, bind=False)
    with pytest.raises(kr.RunError, match="no catalog template matches"):
        kr.run(deck, ps, 700.0, [0.0, 1.0], [1], tmp_path / "refused", allow_unbound_catalog=True)
    assert not (tmp_path / "refused").exists()
    out = tmp_path / "flagged"
    m = kr.run(deck, ps, 700.0, [0.0, 1.0], [1], out, allow_unbound_catalog=True, allow_unmatched_catalog=True)
    assert m["catalog_matches_initial_state"] is False and m["initial_match"]["templates"] == 0
    assert all(r["catalog_matches_initial_state"] is False for r in lines(out / "seed-1" / "observables.jsonl"))
    assert "catalog_matches_initial_state=false" in (out / "comparison.csv").read_text().splitlines()[0]
    ok = kr.run(deck, [pe_catalog_file(tmp_path)], 700.0, [0.0], [1], tmp_path / "ok")
    assert ok["catalog_matches_initial_state"] is True and ok["initial_match"]["templates"] > 0
    argv = [str(deck), *map(str, ps), "--temperature", "700", "--times", "0", "--seeds", "1",
            "--allow-unbound-catalog", "--allow-unmatched-catalog", "--out", str(tmp_path / "cli")]
    assert kr.main(argv) == 0


def test_deck_is_read_without_executing_it(tmp_path):
    """Item 3: nothing in the deck runs."""
    mapping = S.pe_mapping({10: 1})
    sentinel = tmp_path / "executed"
    attack = f"polymer.__globals__['__builtins__']['open']({str(sentinel)!r}, 'w').write('x')\n"
    with pytest.raises(kr.RunError):
        kr.read_deck(write_deck(tmp_path / "attack.py", mapping, extra=attack))
    assert not sentinel.exists()
    other = f"species(label=open({str(sentinel)!r}, 'w').write('x'))\n"
    deck = kr.read_deck(write_deck(tmp_path / "other.py", mapping, extra=other))
    assert deck.spec == LinearHomopolymerSpec.from_mapping(mapping)
    assert not sentinel.exists()


@pytest.mark.parametrize("variant", ["conditional", "function", "loop", "alias", "variable", "starstar",
                                     "call_argument", "cutoff_expression", "nested"])
def test_deck_polymer_call_must_be_one_literal_module_level_call(tmp_path, variant):
    """Item 3: one module-level polymer() call with literal arguments, or a refusal."""
    mapping = S.pe_mapping({10: 1})
    label, lit = repr(mapping["label"]), repr(mapping)
    call = f"polymer({label}, linear_homopolymer_spec={lit})"
    body = {
        "conditional": f"if True:\n    {call}\n",
        "function": f"def make():\n    {call}\nmake()\n",
        "loop": f"for _ in [1]:\n    {call}\n",
        "alias": f"p = polymer\np({label}, linear_homopolymer_spec={lit})\n",
        "variable": f"spec = {lit}\npolymer({label}, linear_homopolymer_spec=spec)\n",
        "starstar": f"kw = {{'linear_homopolymer_spec': {lit}}}\npolymer({label}, **kw)\n",
        "call_argument": f"polymer({label}, linear_homopolymer_spec=dict({lit}))\n",
        "cutoff_expression": f"polymer({label}, linear_homopolymer_spec={lit}, cutoff=int('3'))\n",
        "nested": f"x = [{call}]\n",
    }[variant]
    path = tmp_path / "deck.py"
    path.write_text(DECK_HEAD + body + DECK_TAIL)
    with pytest.raises(kr.RunError):
        kr.read_deck(path)


def test_capped_seeds_censor_the_comparison(tmp_path):
    """Item 4: a row averaged over fewer seeds than requested is censored."""
    deck, catalog = pe_deck(tmp_path), pe_catalog_file(tmp_path)
    times = [0.0, 1.0, 4.0]
    seeds = [1, 2, 3, 4]
    free = kr.run(deck, [catalog], 700.0, times, seeds, tmp_path / "free")
    events = {r["seed"]: r["events"] for r in free["runs"]}
    cap = max(events.values()) - 1
    capped = {s for s, n in events.items() if n >= cap}
    assert min(events.values()) < cap and capped, events
    out = tmp_path / "capped"
    m = kr.run(deck, [catalog], 700.0, times, seeds, out, max_events=cap)
    censoring = m["censoring"]
    assert censoring["seeds_requested"] == 4
    stopped = {s["seed"]: s for s in censoring["stopped_early"]}
    assert set(stopped) == capped
    assert all(s["stop_reason"] == "max_events" and s["observations"] < 3 for s in stopped.values())
    assert 4.0 in censoring["censored_times_s"] and 0.0 not in censoring["censored_times_s"]
    rows = read_comparison(out / "comparison.csv")
    first = [r for r in rows if float(r["time_s"]) == 0.0]
    assert first and all(r["censored"] == "false" and r["kmc_mean"] != "" for r in first)
    last = [r for r in rows if float(r["time_s"]) == 4.0]
    assert last and all(r["censored"] == "true" for r in last)
    for r in last:
        assert r["kmc_mean"] == "" and r["kmc_stderr"] == "" and r["difference"] == ""
        assert int(r["seeds_requested"]) == 4
    reached = [s for s in seeds if s not in capped]
    mu0 = next(r for r in last if r["quantity"] == "mu0_mol")
    assert int(mu0["seeds"]) == len(reached)
    want = [mo.mom_reported_quantities(lines(out / f"seed-{s}" / "observables.jsonl")[-1])["mu0_mol"]
            for s in reached]
    assert rel(float(mu0["censored_kmc_mean"]), math.fsum(want) / len(want)) <= REL


def test_mom_times_match_exactly_and_unambiguously(tmp_path):
    """Item 5."""
    deck, catalog = pe_deck(tmp_path), pe_catalog_file(tmp_path)
    q = {"mu0_mol": 1.0}

    def mom(name, text):
        path = tmp_path / name
        path.write_text(text)
        return path

    near = mom("near.json", json.dumps([{"time_s": 1.0000000001, "quantities": q}]))
    with pytest.raises(kr.RunError, match="time"):
        kr.load_mom(near, [0.0, 1.0])
    for token in ("Infinity", "-Infinity", "NaN"):
        bad = mom(f"bad{token}.json", '[{"time_s": %s, "quantities": {"mu0_mol": 1.0}}]' % token)
        with pytest.raises(kr.RunError, match="finite"):
            kr.load_mom(bad, [0.0, 1.0])
    nan_value = mom("nan_value.json", '[{"time_s": 1.0, "quantities": {"mu0_mol": NaN}}]')
    with pytest.raises(kr.RunError, match="finite"):
        kr.load_mom(nan_value, [0.0, 1.0])
    exact = mom("exact.json", json.dumps([{"time_s": 1.0, "quantities": q}]))
    got, record = kr.load_mom(exact, [0.0, 1.0])
    assert got == [None, q] and record["times"] == 1
    within = mom("within.json", json.dumps([{"time_s": 1.0 + 2.0e-16, "quantities": q}]))
    assert kr.load_mom(within, [0.0, 1.0])[0] == [None, q]
    with pytest.raises(kr.RunError, match="ambiguous"):
        kr.load_mom(exact, [0.0, 1.0, 1.0 + 1.0e-13])
    with pytest.raises(kr.RunError):
        kr.run(deck, [catalog], 700.0, [0.0, 1.0], [1], tmp_path / "near_run", mom=near)
    assert not (tmp_path / "near_run").exists()


def test_a_failed_seed_is_labelled_and_recorded(tmp_path, monkeypatch):
    """Item 6: an invariant failure after the t = 0 row."""
    deck, catalog = pe_deck(tmp_path), pe_catalog_file(tmp_path)
    original = mo.InvariantChecker.check

    def failing(self, state, event=None):
        if self.checks >= 2:
            raise mo.InvariantError("injected invariant failure")
        return original(self, state, event)

    monkeypatch.setattr(mo.InvariantChecker, "check", failing)
    out = tmp_path / "failed"
    m = kr.run(deck, [catalog], 700.0, [0.0, 1.0, 5.0], [1, 2], out)
    on_disk = json.loads((out / "run.json").read_text())
    assert m["status"] == on_disk["status"] == "failed"
    for r in on_disk["runs"]:
        assert r["stop_reason"] == "error"
        assert "InvariantError" in r["error"] and "injected" in r["error"]
        rows = lines(out / f"seed-{r['seed']}" / "observables.jsonl")
        assert 1 <= len(rows) == r["observations"] < 3
        for row in rows:
            assert row["seed"] == r["seed"] and row["catalog_complete"] is True
            assert row["catalog_gate_failures"] == [] and row["catalog_bound"] is True
    assert on_disk["leak"]["excluded_seeds"] == [1, 2] and on_disk["leak"]["per_seed"] == []
    assert {s["seed"] for s in on_disk["censoring"]["stopped_early"]} == {1, 2}
    argv = [str(deck), str(catalog), "--temperature", "700", "--times", "0", "5", "--seeds", "1",
            "--out", str(tmp_path / "cli")]
    assert kr.main(argv) == 1


# --------------------------------------------------------------------------
# rework 2 (spar r110)
# --------------------------------------------------------------------------
def _reordered(executable, order, swap_participants=False):
    """The same executable listed differently: atoms in ``order`` (old
    indices), the bond list reversed, participants optionally exchanged."""
    new_of = {old: new for new, old in enumerate(order)}
    atoms = [dict(executable["atoms"][old]) for old in order]
    if swap_participants:
        for a in atoms:
            a["participant"] = 1 - a["participant"]
    bonds = [[new_of[j], new_of[i], x, y] for i, j, x, y in reversed(executable["bonds"])]
    return {"atoms": atoms, "bonds": bonds}


@pytest.mark.parametrize("variant", ["bond_order_reversed", "atoms_permuted", "participants_swapped"])
def test_an_executable_listed_differently_is_still_a_duplicate(tmp_path, variant):
    """Item 1: identity is the template's structure, not its JSON order."""
    from rmgpy.kmc.template_emission import sha256_json
    payloads = toy_payloads([t for t in S.pe_templates() if t.name in PE_LIVE])
    name = "h_abstraction" if variant == "participants_swapped" else "homolysis"   # asymmetric pair
    (copy,) = [p for p in renamed(payloads, "x") if p["templates"][0]["name"] == name]
    ex = copy["templates"][0]["executable"]
    n = len(ex["atoms"])
    order = list(range(n)) if variant == "bond_order_reversed" else list(reversed(range(n)))
    copy["templates"][0]["executable"] = _reordered(ex, order, variant == "participants_swapped")
    assert sha256_json(copy["templates"][0]["executable"]) != sha256_json(ex)
    first = write_catalog(tmp_path / "first.json", payloads, bind=S.pe_mapping({1: 1}))
    second = write_catalog(tmp_path / "second.json", [copy], bind=S.pe_mapping({1: 1}))
    me.PayloadCatalog([copy]).templates(700.0)            # the reordered copy is a valid template
    with pytest.raises(kr.RunError, match="executable"):
        kr.load_catalogs([first, second], 700.0)
    within = write_catalog(tmp_path / "within.json", payloads + [copy], bind=S.pe_mapping({1: 1}))
    with pytest.raises(kr.RunError, match="executable"):
        kr.load_catalogs([within], 700.0)


def test_distinct_templates_are_not_duplicates(tmp_path):
    """Item 1, the converse: the toy set and the sampled PS catalogs load."""
    kr.load_catalogs([pe_catalog_file(tmp_path)], 700.0)
    kr.load_catalogs(ps_catalog_files(tmp_path), 700.0)


def _full_pe_normalized():
    mapping = S.pe_mapping({10: 1, 20: 1}, small={"methyl": S.METHYL}, monomer=S.ETHYLENE_GRAPH)
    ens = LinearHomopolymerSpec.from_mapping(mapping).ensemble_observables
    mapping["Mn_g_mol"], mapping["Mw_g_mol"] = ens.mn_g_mol, ens.mw_g_mol
    return json.loads(LinearHomopolymerSpec.from_mapping(mapping).normalized_json)


def test_chemistry_digest_covers_every_field_but_the_excluded():
    """Item 2: the digest is the normalized declaration minus an explicit list."""
    normalized = _full_pe_normalized()
    assert kr.CHEMISTRY_EXCLUDED == CHEMISTRY_EXCLUDED
    assert set(CHEMISTRY_EXCLUDED) <= set(normalized)      # every excluded field is perturbed below
    base = kr.chemistry_digest(normalized)
    assert base == chemistry_of_normalized(normalized)
    for key in normalized:
        perturbed = dict(normalized, **{key: {"perturbed": key}})
        assert (kr.chemistry_digest(perturbed) != base) == (key not in CHEMISTRY_EXCLUDED), key
    for key in ("small_species", "monomer_product", "repeat_unit", "end_groups", "locality_policy"):
        assert key in normalized and key not in CHEMISTRY_EXCLUDED
    assert kr.chemistry_digest(dict(normalized, future_field=1)) != base    # new fields count by default
    spec = LinearHomopolymerSpec.from_mapping(S.pe_mapping({10: 1}))
    assert kr.chemistry_sha256(spec) == kr.chemistry_digest(json.loads(spec.normalized_json))


def test_a_declared_small_species_changes_the_binding(tmp_path):
    """Item 2: a zero-population declared radical changes MET bypass, so a
    catalog bound to the plain declaration is refused."""
    plain = S.pe_mapping({10: 1})
    with_radical = S.pe_mapping({10: 1}, small={"ethyl": S.ETHYL})
    assert chemistry_sha256(plain) != chemistry_sha256(with_radical)
    deck = write_deck(tmp_path / "deck.py", with_radical)
    with pytest.raises(kr.RunError, match="chemistry_sha256"):
        kr.run(deck, [pe_catalog_file(tmp_path)], 700.0, [0.0], [1], tmp_path / "out")
    # population, mass, density and transport do not
    other = S.pe_mapping({30: 4}, d_ref=3.0e-17)
    other["initial_mass_kg"], other["density_kg_m3"] = 2.0, 900.0
    assert chemistry_sha256(other) == chemistry_sha256(plain)


def test_a_failing_report_does_not_abort_the_ensemble(tmp_path, monkeypatch):
    """Item 3."""
    deck, catalog = pe_deck(tmp_path), pe_catalog_file(tmp_path)

    def broken(self):
        raise RuntimeError("injected report failure")

    monkeypatch.setattr(me.MeltEngine, "report", broken)
    out = tmp_path / "report"
    m = kr.run(deck, [catalog], 700.0, [0.0, 1.0], [1, 2], out)
    assert m["status"] == "failed"
    assert [r["seed"] for r in m["runs"]] == [1, 2]
    for r in m["runs"]:
        assert r["stop_reason"] == "error" and "injected report failure" in r["error"]
        rows = lines(out / f"seed-{r['seed']}" / "observables.jsonl")
        assert len(rows) == 2 and all(row["seed"] == r["seed"] for row in rows)
    assert m["leak"]["excluded_seeds"] == [1, 2]


@pytest.mark.parametrize("bad", ['"NaN"', '"1.0"', "true", "null", "[1.0]", '{"2": "NaN"}', '{"2": false}',
                                 '{"2": {"x": 1.0}}', "1e999"])
def test_mom_quantities_must_be_finite_json_numbers(tmp_path, bad):
    """Item 4."""
    path = tmp_path / "mom.json"
    path.write_text('[{"time_s": 1.0, "quantities": {"mu0_mol": %s}}]' % bad)
    with pytest.raises(kr.RunError, match="number"):
        kr.load_mom(path, [0.0, 1.0])
    good = tmp_path / "good.json"
    good.write_text('[{"time_s": 1.0, "quantities": {"mu0_mol": 1, "explicit_dp_mol": {"2": 1.5e-6}}}]')
    assert kr.load_mom(good, [0.0, 1.0])[0][1] == {"mu0_mol": 1, "explicit_dp_mol": {"2": 1.5e-6}}


@pytest.mark.parametrize("binding", [
    "def polymer(*args, **kwargs):\n    pass\n",
    "async def polymer():\n    pass\n",
    "class polymer:\n    pass\n",
    "polymer = 1\n",
    "import polymer\n",
    "import os as polymer\n",
    "from os import path as polymer\n",
    "from os import *\n",
    "def f():\n    global polymer\n",
    "def f():\n    x = 1\n    def g():\n        nonlocal x\n        polymer = x\n",
    "f = lambda polymer: 0\n",
    "def f(polymer):\n    pass\n",
    "def f(*, polymer=1):\n    pass\n",
    "def f(*polymer):\n    pass\n",
    "try:\n    pass\nexcept Exception as polymer:\n    pass\n",
    "for polymer in []:\n    pass\n",
    "with open('x') as polymer:\n    pass\n",
    "(polymer := 1)\n",
    "[0 for polymer in []]\n",
])
def test_a_deck_that_binds_the_name_polymer_is_refused(tmp_path, binding):
    """Item 5: any binding of the name, before a literal declaration."""
    mapping = S.pe_mapping({10: 1})
    path = tmp_path / "deck.py"
    path.write_text(DECK_HEAD + binding + f"polymer({mapping['label']!r}, linear_homopolymer_spec={mapping!r})\n"
                    + DECK_TAIL)
    with pytest.raises(kr.RunError, match="polymer"):
        kr.read_deck(path)


# --------------------------------------------------------------------------
# rework 3 (spar r111)
# --------------------------------------------------------------------------
def extended_homolysis(family="toy_homolysis", name="homolysis_ext"):
    """The toy homolysis window grown by one CH2 shell on each side: it
    matches a subset of the toy homolysis's cleavages (nested windows)."""
    from rmgpy.kmc.window_matcher import AtomEdit, BondEdit, Template, WindowAtom
    atoms = [WindowAtom("C"), WindowAtom("C"), WindowAtom("H"), WindowAtom("H"), WindowAtom("H"),
             WindowAtom("H"), WindowAtom("C"), WindowAtom("C"), WindowAtom("H"), WindowAtom("H"),
             WindowAtom("H"), WindowAtom("H"), WindowAtom("C", 0, "port", S.MID),
             WindowAtom("C", 0, "port", S.MID)]
    bonds = [(0, 1, 1), (0, 2, 1), (0, 3, 1), (1, 4, 1), (1, 5, 1), (0, 6, 1), (1, 7, 1), (6, 8, 1),
             (6, 9, 1), (7, 10, 1), (7, 11, 1), (6, 12, 1), (7, 13, 1)]
    return Template(name, family, atoms, bonds,
                    [BondEdit(0, 1, 0), AtomEdit(0, radical=1), AtomEdit(1, radical=1)], weight=0.2)


def test_catalogs_sharing_a_family_are_refused(tmp_path):
    """Item 1 (r111): two chemistry-bound catalogs of one family with nested
    windows would both fire the same cleavages."""
    deck = pe_deck(tmp_path, counts={10: 1, 20: 1})
    toy = pe_catalog_file(tmp_path)
    nested = write_catalog(tmp_path / "nested.json", toy_payloads([extended_homolysis()]),
                           bind=S.pe_mapping({1: 1}))
    with pytest.raises(kr.RunError, match=r"famil.*toy_homolysis"):
        kr.load_catalogs([toy, nested], 700.0)
    with pytest.raises(kr.RunError, match="toy_homolysis"):
        kr.run(deck, [toy, nested], 700.0, [0.0], [1], tmp_path / "out")
    assert not (tmp_path / "out").exists()
    # either catalog alone runs
    kr.load_catalogs([nested], 700.0)
    kr.load_catalogs([toy], 700.0)


def test_a_catalog_without_families_cannot_be_combined(tmp_path):
    """Item 1: provenance.families is required to combine catalogs."""
    first, second = ps_catalog_files(tmp_path)
    doc = json.loads(second.read_text())
    del doc["provenance"]["families"]
    bare = tmp_path / "bare.json"
    bare.write_text(json.dumps(doc))
    with pytest.raises(kr.RunError, match="families"):
        kr.load_catalogs([first, bare], 700.0)
    kr.load_catalogs([bare], 700.0)                     # alone it loads
    kr.load_catalogs([first, second], 700.0)           # disjoint families combine


def test_a_site_matched_from_two_catalogs_is_refused(tmp_path):
    """Item 1: the initial-state check catches nested windows that the family
    policy misses (here declared under different family names)."""
    deck = pe_deck(tmp_path, counts={10: 1, 20: 1})
    toy = pe_catalog_file(tmp_path)
    nested = write_catalog(tmp_path / "nested.json",
                           toy_payloads([extended_homolysis(family="toy_homolysis_shell")]),
                           bind=S.pe_mapping({1: 1}))
    kr.load_catalogs([toy, nested], 700.0)             # disjoint families, distinct executables
    with pytest.raises(kr.RunError, match="site"):
        kr.run(deck, [toy, nested], 700.0, [0.0], [1], tmp_path / "out")
    assert not (tmp_path / "out").exists()


def _beta_variants():
    """Two beta scissions that differ only in which carbon carries which port:
    the same nodes and the same pairs, different pair labels."""
    from rmgpy.kmc.window_matcher import AtomEdit, BondEdit, Template, WindowAtom
    branch = ((1.0, "C", False), (1.0, "C", False), (1.0, "H", False))
    atoms = [WindowAtom("C", 0, radical=1), WindowAtom("H"), WindowAtom("C", 0, "port", S.MID),
             WindowAtom("C"), WindowAtom("H"), WindowAtom("H"), WindowAtom("C"), WindowAtom("H"),
             WindowAtom("H"), WindowAtom("C", 0, "port", branch)]
    common = [(0, 1, 1), (0, 3, 1), (3, 4, 1), (3, 5, 1), (3, 6, 1), (6, 7, 1), (6, 8, 1)]
    edits = [BondEdit(3, 6, 0), BondEdit(0, 3, 2), AtomEdit(0, radical=0), AtomEdit(6, radical=1)]
    x = Template("beta_x", "toy_beta_scission", atoms, common + [(0, 2, 1), (6, 9, 1)], edits, weight=5.0)
    y = Template("beta_y", "toy_beta_scission", atoms, common + [(0, 9, 1), (6, 2, 1)], edits, weight=5.0)
    return x, y


def test_executables_differing_only_in_a_pair_label_are_not_duplicates(tmp_path):
    """Item 2: the identity reads pair labels (bond order before / after)."""
    x, y = _beta_variants()
    gx, gy = kr._template_graph(x), kr._template_graph(y)
    assert sorted(gx.nodes(data=True)) == sorted(gy.nodes(data=True))
    assert set(map(frozenset, gx.edges())) == set(map(frozenset, gy.edges()))
    assert not kr._same_graph(gx, gy)
    assert kr._same_graph(gx, kr._template_graph(x))
    catalog = write_catalog(tmp_path / "beta.json", toy_payloads([x, y]), bind=S.pe_mapping({1: 1}))
    kr.load_catalogs([catalog], 700.0)
    # nested windows are different executables (only the families policy and the site check see them)
    h = next(t for t in S.pe_templates() if t.name == "homolysis")
    assert not kr._same_graph(kr._template_graph(h), kr._template_graph(extended_homolysis()))


def extended_recombination(family="toy_recombination_shell"):
    """The toy end-radical recombination with participant 0's window grown by
    one CH2 shell (nested in the toy recombination)."""
    from rmgpy.kmc.window_matcher import AtomEdit, BondEdit, Template, WindowAtom
    atoms = [WindowAtom("C", 0, radical=1), WindowAtom("H", 0), WindowAtom("H", 0), WindowAtom("C", 0),
             WindowAtom("H", 0), WindowAtom("H", 0), WindowAtom("C", 0, "port", S.MID)] + S._ch2_radical(1)
    bonds = [(0, 1, 1), (0, 2, 1), (0, 3, 1), (3, 4, 1), (3, 5, 1), (3, 6, 1), (7, 8, 1), (7, 9, 1), (7, 10, 1)]
    return Template("recombination_shell", family, atoms, bonds,
                    [BondEdit(0, 7, 1), AtomEdit(0, radical=0), AtomEdit(7, radical=0)], weight=1.0e-2)


def test_a_bimolecular_site_matched_from_two_catalogs_is_found(tmp_path):
    """Item 1, the bimolecular branch of the site check: two end radicals
    (one mid-chain cleavage) are a recombination pair for both the toy
    template and its nested extension in another catalog."""
    import numpy as np
    from rmgpy.kmc.melt_state import materialize
    toy = pe_catalog_file(tmp_path)
    nested = write_catalog(tmp_path / "nested_rec.json", toy_payloads([extended_recombination()]),
                           bind=S.pe_mapping({1: 1}))
    catalog, records = kr.load_catalogs([toy, nested], 700.0)
    owner = {name: i for i, r in enumerate(records) for name in r.template_names}
    spec = LinearHomopolymerSpec.from_mapping(S.pe_mapping({20: 1}))
    state = materialize(spec)
    cs = [c for c in S.carbons(state, 0)
          if sum(1 for b in state.atoms[c].edges if state.atoms[b].element == "C") == 2]
    a = cs[len(cs) // 2]
    b = next(x for x in state.atoms[a].edges if x in cs and x != a)
    from rmgpy.kmc.melt_state import SetBond, SetRadical
    state.apply([SetBond(a, b, 1.0, 0.0), SetRadical(a, 0, 1), SetRadical(b, 0, 1)])
    engine = me.MeltEngine(spec, catalog, 700.0, np.random.default_rng(0), state=state)
    found = kr.cross_catalog_site(engine, owner)
    assert found is not None and "bimolecular" in found and "catalog 0" in found and "catalog 1" in found
    alone = me.MeltEngine(spec, me.PayloadCatalog(json.loads(toy.read_text())["payloads"]), 700.0,
                          np.random.default_rng(0), state=state)
    assert kr.cross_catalog_site(alone, {t.name: 0 for t in alone.templates}) is None


# --------------------------------------------------------------------------
# I-164: the stop reason is the time series' own, passed through to run.json
# --------------------------------------------------------------------------
def test_run_records_the_stop_reason_the_time_series_returns(tmp_path, monkeypatch):
    """run.py passes its event cap to ``write_time_series`` and records the
    reason it returns; it does not infer one (a forged reason shows through)."""
    deck, catalog = pe_deck(tmp_path), pe_catalog_file(tmp_path)
    real = mo.write_time_series
    seen = []

    def forged(engine, times, path, **kw):
        out = real(engine, times, path, **kw)
        seen.append((kw.get("max_events"), out.stop_reason))
        out.stop_reason = "absorbing" if out.stop_reason != "absorbing" else "horizon"
        return out

    monkeypatch.setattr(mo, "write_time_series", forged)
    m = kr.run(deck, [catalog], 700.0, [0.0, 1.0], [1, 2], tmp_path / "forged", max_events=10 ** 6)
    assert [cap for cap, _r in seen] == [10 ** 6, 10 ** 6]
    for r, (_cap, real_reason) in zip(m["runs"], seen):
        assert real_reason == "horizon"
        assert r["stop_reason"] == "absorbing"


def test_each_seed_series_ends_with_its_stop_record(tmp_path):
    deck, catalog = pe_deck(tmp_path), pe_catalog_file(tmp_path)
    times = [0.0, 1.0, 3.0, 50.0]
    m = kr.run(deck, [catalog], 700.0, times, [1, 2], tmp_path / "cap", max_events=5)
    for r in m["runs"]:
        path = tmp_path / "cap" / r["path"]
        stops = records(path)
        assert len(stops) == 1 and json.loads(path.read_text().splitlines()[-1]) == stops[0]
        assert stops[0]["record"] == "stop" and stops[0]["stop_reason"] == r["stop_reason"] == "max_events"
        assert stops[0]["events"] == r["events"] == 5
        assert stops[0]["observations"] == r["observations"] == len(lines(path)) < 4
    for s in m["censoring"]["stopped_early"]:
        assert s["stop_reason"] == "max_events"
        assert s["last_observed_time_s"] == lines(tmp_path / "cap" / f"seed-{s['seed']}" / "observables.jsonl")[-1]["time_s"]
    free = kr.run(deck, [catalog], 700.0, times[:2], [1], tmp_path / "free")
    (r,) = free["runs"]
    (stop,) = records(tmp_path / "free" / r["path"])
    assert stop["stop_reason"] == r["stop_reason"] == "horizon" and r["observations"] == 2
