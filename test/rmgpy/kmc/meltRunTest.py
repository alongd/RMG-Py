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


def write_catalog(path, payloads, *, complete=True, failures=(), tag="test"):
    failures = list(failures)
    doc = {"schema": "kmc-template-catalog/1",
           "provenance": {"artifact_id": f"artifact-{tag}", "families": sorted({p["family"] for p in payloads}),
                          "code_sha": None, "temperatures": [300.0, 1000.0]},
           "payloads": payloads, "refusals": [],
           "gates": {"failures": failures}, "completeness": {}, "complete": complete}
    Path(path).write_text(json.dumps(doc, sort_keys=True, indent=1) + "\n")
    return Path(path)


def pe_catalog_file(tmp_path, name="pe_toy_catalog.json", **kw):
    return write_catalog(tmp_path / name, toy_payloads(
        [t for t in S.pe_templates() if t.name in PE_LIVE]), tag="pe-toy", **kw)


def ps_catalog_files(tmp_path):
    return [write_catalog(tmp_path / f"ps_{i}.json", S.payloads(name), tag=f"ps-{i}")
            for i, name in enumerate(S.PS_CATALOG)]


def pe_deck(tmp_path, counts=None, **kw):
    return write_deck(tmp_path / "pe_deck.py", S.pe_mapping(counts or {10: 1, 20: 1}), **kw)


def ps_deck(tmp_path, counts=None, **kw):
    return write_deck(tmp_path / "ps_deck.py", ps_deck_mapping(counts or {6: 3}), **kw)


def lines(path):
    return [json.loads(x) for x in Path(path).read_text().splitlines()]


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
@pytest.mark.parametrize("via_variable", [False, True])
def test_deck_declaration_equals_from_mapping(tmp_path, polymer, via_variable):
    mapping = ps_deck_mapping({6: 3}) if polymer == "ps" else S.pe_mapping({10: 1, 20: 1})
    path = write_deck(tmp_path / "deck.py", mapping, via_variable=via_variable, cutoff=3)
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


def test_deck_reader_follows_the_input_file_vocabulary():
    """The deck reader stubs exactly the names the RMG input reader exposes and
    binds polymer() arguments in the input reader's order."""
    import rmgpy.rmg.input as rmg_input
    tree = ast.parse(INPUT_PY.read_text())
    reader = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "read_input_file")
    names = None
    for node in ast.walk(reader):
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == "local_context" for t in node.targets):
            names = {k.value for k in node.value.keys}
    assert names is not None
    assert names - {"__builtins__", "True", "False"} == set(kr.DECK_VOCABULARY)
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
