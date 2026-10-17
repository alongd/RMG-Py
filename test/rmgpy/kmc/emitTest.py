"""Tests for rmgpy/kmc/emit.py (I-155): a compiled event-set artifact plus a
roster of declared reference systems becomes one gate-checked template catalog.

Database-free.  The reference systems are declared as a ``LinearHomopolymerSpec``
mapping plus chain variants; ``emit`` resolves the repeat unit to a registered
provider by graph isomorphism.  Assembly runs on stored row fixtures
(``fixtures/i155_rows_dp8.jsonl``: rows of the I-138 acceptance run on the closed
DP8 reference, template maps dropped) with stub services; the database-bound
balance functions (``derive_weights``, ``check_cycles``) are replaced where a
test needs a chosen verdict.
"""

import ast
import copy
import inspect
import json
import re
import sys
import textwrap
from pathlib import Path

import pytest

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import emit
from rmgpy.kmc import template_emission as te
from rmgpy.kmc import windows as win
from rmgpy.polymer_spec import LinearHomopolymerSpec

HERE = Path(__file__).resolve().parent
ROWS = HERE / "fixtures" / "i155_rows_dp8.jsonl"
S = b.PsChainSpec


# --------------------------------------------------------------------------
# declared polymers
# --------------------------------------------------------------------------
def _atom(key, element="C", hydrogens=0):
    return {"key": key, "element": element, "hydrogens": hydrogens,
            "electronic_state": {"radical_electrons": 0}}


def _cap():
    return {"atoms": [_atom("h", "H", 0)], "bonds": [],
            "ports": {"cap": {"atom": "h", "bond_order": 1}},
            "electronic_state": {"multiplicity": 1}}


def polymer_mapping(kind):
    """A ``LinearHomopolymerSpec`` mapping: ``ps`` (styrene repeat), ``pe``
    (ethylene), ``pp`` (propylene: no provider declares it)."""
    ring = [_atom(f"p{i}", hydrogens=0 if i == 1 else 1) for i in range(1, 7)]
    ring_bonds = [{"atoms": [f"p{i}", f"p{i % 6 + 1}"], "order": 1.5} for i in range(1, 7)]
    if kind == "ps":
        atoms = [_atom("c1", hydrogens=2), _atom("c2", hydrogens=1)] + ring
        bonds = [{"atoms": ["c1", "c2"], "order": 1}, {"atoms": ["c2", "p1"], "order": 1}] + ring_bonds
    elif kind == "pe":
        atoms = [_atom("c1", hydrogens=2), _atom("c2", hydrogens=2)]
        bonds = [{"atoms": ["c1", "c2"], "order": 1}]
    elif kind == "pp":
        atoms = [_atom("c1", hydrogens=2), _atom("c2", hydrogens=1), _atom("c3", hydrogens=3)]
        bonds = [{"atoms": ["c1", "c2"], "order": 1}, {"atoms": ["c2", "c3"], "order": 1}]
    else:
        raise ValueError(kind)
    repeat = {"atoms": atoms, "bonds": bonds,
              "ports": {"head": {"atom": "c1", "bond_order": 1},
                        "tail": {"atom": "c2", "bond_order": 1}},
              "electronic_state": {"multiplicity": 1}}
    return {
        "schema_version": 1, "label": kind, "repeat_unit": repeat,
        "end_groups": {"head": _cap(), "tail": _cap()},
        "small_species": {},
        "initial_ensemble": {"dp_counts": {"8": 1}},
        "initial_mobile_counts": {}, "density_basis": "whole_initial_mixture",
        "initial_mass_kg": 1.25, "density_kg_m3": 1050.0,
        "locality_policy": {"radius": 2, "bounds": {
            "elements": ["C", "H"], "isotopes": [-1], "formal_charges": [0],
            "radical_electrons": [0, 1], "lone_pairs": [0],
            "bond_orders": [1, 1.5, 2, 3], "spin_multiplicities": [1, 2, 3],
            "maximum_heavy_atom_degree": 4,
            "maximum_explicit_hydrogens_per_heavy_atom": 4,
            "maximum_total_hydrogens_on_free_small_species": 12,
            "maximum_family_arity": 2, "maximum_root_tuple_size": 4,
            "maximum_free_small_species_heavy_atoms": 8,
            "attachment_classes": [], "topology_metadata_schema": {}}},
        "family_manifest": {}, "library_manifest": {}, "transport_model": {}, "diagnostics": {},
    }


def chain(dp, head=None, tail=None, radical_unit=None):
    return {"chain": {"dp": dp, "head": head, "tail": tail, "radical_unit": radical_unit}}


def write_roster(tmp_path, kind, systems, name="roster.json"):
    path = tmp_path / name
    path.write_text(json.dumps({"schema": emit.ROSTER_SCHEMA, "polymer": polymer_mapping(kind),
                                "systems": systems}))
    return path


def systems_of(tmp_path, kind, systems):
    roster = emit.load_roster(write_roster(tmp_path, kind, systems))
    return emit.build_reference_systems(roster)


def isomorphic(build_a, build_b):
    return build_a.molecule.is_isomorphic(build_b.molecule)


# --------------------------------------------------------------------------
# item 1: reference systems from a declaration
# --------------------------------------------------------------------------
def test_ps_declaration_specs_equal_the_harness_specs(tmp_path):
    systems, refusals = systems_of(tmp_path, "ps", {
        "dp8_closed": [chain(8)],
        "dp4+pe": [chain(4), chain(1, tail="benzylic_radical")],
        "pe2_identical": [chain(2, tail="benzylic_radical"), chain(2, tail="benzylic_radical")]})
    assert refusals == []
    assert [p.spec for p in systems["dp8_closed"]] == [S(8)]
    assert [p.spec for p in systems["dp4+pe"]] == [S(4), S(1, tail_end="benzylic_radical")]
    assert [p.spec for p in systems["pe2_identical"]] == [S(2, tail_end="benzylic_radical")] * 2


def test_ps_declaration_molecules_are_isomorphic_to_the_harness_ones(tmp_path):
    systems, _ = systems_of(tmp_path, "ps", {
        "a": [chain(8)], "b": [chain(4), chain(1, tail="benzylic_radical")],
        "c": [chain(2, tail="benzylic_radical")]})
    pairs = [(systems["a"][0], S(8)), (systems["b"][0], S(4)),
             (systems["b"][1], S(1, tail_end="benzylic_radical")),
             (systems["c"][0], S(2, tail_end="benzylic_radical"))]
    for participant, spec in pairs:
        assert isomorphic(participant.build("p0"), win.build_spec(spec, "p0"))


def test_pe_declaration_is_the_n_alkane(tmp_path):
    systems, refusals = systems_of(tmp_path, "pe", {"c16": [chain(8)]})
    assert refusals == []
    (participant,) = systems["c16"]
    assert participant.spec == ("alkane", 16)
    build = participant.build("p0")
    assert isomorphic(build, win.build_spec(("alkane", 16), "p0"))
    assert build.molecule.get_formula() == "C16H34"


def test_pe_declaration_with_ends_is_the_alkane_end_spec(tmp_path):
    systems, _ = systems_of(tmp_path, "pe", {"r": [chain(3, head="primary_radical", tail="methyl")]})
    assert systems["r"][0].spec == ("alkane", 6, "primary_radical", "methyl")


def test_unregistered_repeat_unit_is_refused_with_a_row(tmp_path):
    systems, refusals = systems_of(tmp_path, "pp", {"pp8": [chain(8)]})
    assert systems == {} or systems.get("pp8") in (None, [])
    assert [r["status"] for r in refusals] == ["refused:no_provider"]
    assert refusals[0]["system"] == "pp8" and refusals[0]["reason"]


def test_unknown_end_variant_is_refused_never_dropped(tmp_path):
    systems, refusals = systems_of(tmp_path, "ps", {"x": [chain(4, tail="bogus")]})
    assert [r["status"] for r in refusals] == ["refused:variant"]
    assert "x" not in systems


def test_interior_radical_on_the_alkane_is_refused_not_built_closed(tmp_path):
    systems, refusals = systems_of(tmp_path, "pe", {"r": [chain(5, radical_unit=2)]})
    assert [r["status"] for r in refusals] == ["refused:variant"] and "r" not in systems


def test_ps_interior_radical_is_the_radical_unit_spec(tmp_path):
    systems, _ = systems_of(tmp_path, "ps", {"r": [chain(5, radical_unit=2)]})
    assert systems["r"][0].spec == S(5, radical_units=(2,))


def test_a_refused_participant_refuses_its_whole_system(tmp_path):
    systems, refusals = systems_of(tmp_path, "ps", {
        "mixed": [chain(4), chain(2, tail="bogus")], "ok": [chain(2)]})
    assert "mixed" not in systems and "ok" in systems
    assert len(refusals) == 1


def test_explicit_small_molecules(tmp_path):
    systems, refusals = systems_of(tmp_path, "ps", {"s": [
        {"smiles": "C"},
        {"adjacency_list": "1 H u1 p0 c0"}]})
    assert refusals == []
    methane, h_atom = systems["s"]
    assert methane.build("p0").molecule.get_formula() == "CH4"
    assert h_atom.build("p1").molecule.get_formula() == "H"
    assert methane.build("p0").metadata.provider == "small"


def test_roster_sha_is_the_file_hash(tmp_path):
    import hashlib
    path = write_roster(tmp_path, "ps", {"a": [chain(2)]})
    assert emit.load_roster(path).sha256 == hashlib.sha256(path.read_bytes()).hexdigest()


def test_no_chemistry_specific_code_in_emit():
    text = Path(emit.__file__).read_text()
    found = re.findall(r"polystyrene|PsChainSpec|alkane|polyethylene|styrene", text, flags=re.IGNORECASE)
    found += re.findall(r"\bPS\b|\bPE\b", text)
    assert not found, found


# --------------------------------------------------------------------------
# harness equivalence of the moved pure functions
# --------------------------------------------------------------------------
def _harness():
    sys.path.insert(0, str(HERE))
    import template_emission_acceptance as acc
    return acc


@pytest.mark.parametrize("name", ["record_key", "bind", "molecularity"])
def test_moved_harness_functions_are_verbatim(name):
    def body(function):
        node = ast.parse(textwrap.dedent(inspect.getsource(function))).body[0]
        statements = node.body
        if isinstance(statements[0], ast.Expr) and isinstance(statements[0].value, ast.Constant):
            statements = statements[1:]       # the docstring
        return ast.dump(node.args), [ast.dump(x) for x in statements]
    assert body(getattr(_harness(), name)) == body(getattr(emit, name))


# --------------------------------------------------------------------------
# record selection
# --------------------------------------------------------------------------
def test_chunks_partition_the_sorted_records():
    records = [{"event_id": f"evt_{k:03d}", "family": "F"} for k in (5, 1, 9, 3, 7, 2, 8)]
    parts = [emit.select_records(records, f"{k}/3") for k in range(3)]
    ids = sorted(r["event_id"] for part in parts for r in part)
    assert ids == sorted(r["event_id"] for r in records)
    assert all(not ({r["event_id"] for r in a} & {r["event_id"] for r in c})
               for i, a in enumerate(parts) for c in parts[i + 1:])
    assert [r["event_id"] for r in parts[0]] == ["evt_001", "evt_005", "evt_009"]


# --------------------------------------------------------------------------
# assembly
# --------------------------------------------------------------------------
class Services:
    counter = None
    kc = None


def ok_midpoint(row, worst=0.0):
    """The row's per-template midpoint entries, as the emit step stores them."""
    return {t["key"]: {"status": "ok", "worst": worst, "worst_T": 650.0, "n": 28}
            for t in (row.get("payload") or {}).get("templates", [])
            if t["status"] == "emitted" and t.get("weight")}


def fixture_rows():
    rows = [json.loads(line) for line in ROWS.read_text().splitlines()]
    rows.sort(key=lambda r: r["event_id"])
    for row in rows:
        if row.get("payload"):
            row["midpoint"] = ok_midpoint(row)
    return rows


def record_id(row):
    return (row.get("payload") or {}).get("record_event_id") or row["event_id"] + "0" * 44


def header(chunk, rows, **over):
    out = {"kind": "header", "schema": emit.ROWS_SCHEMA, "artifact_id": "art",
           "artifact_sha256": "a" * 64, "artifact_records": 6, "families": ["R_Recombination"],
           "database_path": "/db", "database_head": "d" * 40, "code_sha": "c" * 40,
           "roster_sha256": "r" * 64, "temperatures": [600.0, 700.0, 800.0], "chunk": chunk,
           "selection": "all", "assigned": [record_id(r) for r in rows], "roster_refusals": [],
           "systems": {"dp8_closed": ["chain dp=8"]}}
    out.update(over)
    return out


def write_rows(path, head, rows):
    with open(path, "w") as handle:
        handle.write(json.dumps(head) + "\n")
        for row in rows:
            handle.write(json.dumps(dict(row, kind="row")) + "\n")
    return path


def shards(tmp_path, n=2, rows=None, **over):
    rows = rows or fixture_rows()
    paths = []
    for k in range(n):
        mine = rows[k::n]
        paths.append(write_rows(tmp_path / f"rows.{k}of{n}.jsonl", header(f"{k}/{n}", mine, **over), mine))
    return paths


@pytest.fixture(autouse=True)
def balance_stubs(monkeypatch):
    """The database-bound balance functions: identity unless a test sets a verdict."""
    calls = {"derive": 0, "cycles": 0}

    def derive(payloads, counter):
        calls["derive"] += 1
        return payloads

    def cycles(payloads, counter, kc, complete=True, **kw):
        calls["cycles"] += 1
        return payloads, []
    monkeypatch.setattr(te, "derive_weights", derive)
    monkeypatch.setattr(te, "check_cycles", cycles)
    return calls


def run(paths, **kw):
    return emit.assemble([str(p) for p in paths], services=Services(), **kw)


def test_catalog_schema_and_provenance(tmp_path):
    catalog = run(shards(tmp_path))
    assert catalog["schema"] == "kmc-template-catalog/1"
    assert set(catalog) == {"schema", "provenance", "payloads", "refusals", "gates",
                            "completeness", "complete"}
    prov = catalog["provenance"]
    assert prov["artifact_id"] == "art" and prov["artifact_sha256"] == "a" * 64
    assert prov["database_path"] == "/db" and prov["database_head"] == "d" * 40
    assert prov["code_sha"] == "c" * 40 and prov["roster_sha256"] == "r" * 64
    assert prov["families"] == ["R_Recombination"]
    assert prov["temperatures"] == [600.0, 700.0, 800.0]
    assert set(catalog["gates"]) >= {"collisions", "cycles", "balance", "failures"}
    assert set(catalog["gates"]["collisions"]) >= {"keys", "merged", "refused"}
    assert set(catalog["gates"]["cycles"]) >= {"pairs", "refused"}
    assert catalog["complete"] is True and catalog["gates"]["failures"] == []


def test_payloads_hold_exactly_the_fields_loaded_templates_reads(tmp_path):
    catalog = run(shards(tmp_path))
    assert catalog["payloads"]
    for payload in catalog["payloads"]:
        assert set(payload) == {"record_event_id", "family", "reverse_of", "templates"}
        for t in payload["templates"]:
            assert t["status"] == "emitted"
            assert set(t) <= {"key", "name", "family", "status", "direction", "executable",
                              "weight", "derived"}


def test_refusals_and_completeness(tmp_path):
    catalog = run(shards(tmp_path))
    refused = {r["record_event_id"][:24] for r in catalog["refusals"]}
    assert refused == {r["event_id"] for r in fixture_rows() if r["status"] == "refused"}
    assert all(r["status"] and r["reason"] is not None for r in catalog["refusals"])
    totals = {s: sum(c[s] for c in catalog["completeness"].values())
              for s in ("emitted", "partial", "refused", "error")}
    assert totals == {"emitted": 4, "partial": 0, "refused": 2, "error": 0}
    assert all(k.startswith("R_Recombination|") for k in catalog["completeness"])


def test_two_shards_equal_one_unsharded_row_set(tmp_path):
    rows = fixture_rows()
    one = write_rows(tmp_path / "all.jsonl", header(None, rows), rows)
    first = emit.dumps(run([one]))
    two = emit.dumps(run(shards(tmp_path, 2)))
    three = emit.dumps(run(list(reversed(shards(tmp_path, 3)))))
    assert first == two == three


def test_dropped_shard_is_incomplete_and_listed(tmp_path):
    paths = shards(tmp_path, 3)
    catalog = run(paths[:2])
    assert catalog["complete"] is False
    assert any("missing shard" in f and "2/3" in f for f in catalog["gates"]["failures"])


def test_duplicate_shard_is_listed(tmp_path):
    paths = shards(tmp_path, 2)
    catalog = run(paths + [paths[0]])
    assert catalog["complete"] is False
    assert any("duplicate" in f for f in catalog["gates"]["failures"])


def test_mixed_artifacts_are_listed(tmp_path):
    rows = fixture_rows()
    a = write_rows(tmp_path / "a.jsonl", header("0/2", rows[0::2]), rows[0::2])
    c = write_rows(tmp_path / "c.jsonl", header("1/2", rows[1::2], artifact_sha256="b" * 64), rows[1::2])
    catalog = run([a, c])
    assert catalog["complete"] is False
    assert any("artifact_sha256" in f for f in catalog["gates"]["failures"])


def test_error_row_is_incomplete_with_the_failure_listed(tmp_path):
    rows = fixture_rows()
    bad = {"event_id": "evt_ffffffffffffffffffff", "family": "R_Recombination",
           "orientation": "as_generated", "status": "error", "refusal": "error:exception",
           "reason": "Traceback: boom", "seconds": 0.1}
    rows = rows + [bad]
    bad_header = header("0/1", rows, artifact_records=7)
    path = write_rows(tmp_path / "e.jsonl", bad_header, rows)
    catalog = run([path])
    assert catalog["complete"] is False
    assert any("error" in f and "boom" in f for f in catalog["gates"]["failures"])
    assert any(r["status"] == "error" for r in catalog["refusals"])


def test_injected_collision_is_incomplete_with_the_failure_listed(tmp_path):
    rows = fixture_rows()
    emitted = [r for r in rows if r["status"] == "emitted"]
    donor, host = emitted[0], emitted[1]
    clone = copy.deepcopy(donor["payload"]["templates"][0])
    clone["weight"] = dict(clone["weight"], w=[w * 1.0e3 for w in clone["weight"]["w"]])
    host["payload"]["templates"] = host["payload"]["templates"] + [clone]
    host["payload"] = te.finish(host["payload"])
    path = write_rows(tmp_path / "c.jsonl", header(None, rows), rows)
    catalog = run([path])
    assert catalog["gates"]["collisions"]["refused"] >= 1
    assert catalog["complete"] is False
    assert any("collision" in f for f in catalog["gates"]["failures"])


def test_injected_cycle_failure_is_incomplete_with_the_failure_listed(tmp_path, monkeypatch):
    def cycles(payloads, counter, kc, complete=True, **kw):
        return payloads, [{"forward": "evt_" + "a" * 20, "reverse": "evt_" + "b" * 20,
                           "status": "refused:cycle", "reason": "ratio off by 1e3"}]
    monkeypatch.setattr(te, "check_cycles", cycles)
    catalog = run(shards(tmp_path))
    assert catalog["gates"]["cycles"]["pairs"] == 1 and catalog["gates"]["cycles"]["refused"] == 1
    assert catalog["complete"] is False
    assert any("cycle" in f and "ratio off" in f for f in catalog["gates"]["failures"])


def test_assemble_runs_derive_weights_then_check_cycles(tmp_path, balance_stubs):
    run(shards(tmp_path))
    assert balance_stubs == {"derive": 1, "cycles": 1}


def test_check_cycles_receives_completeness_of_coverage(tmp_path, monkeypatch):
    seen = []

    def cycles(payloads, counter, kc, complete=True, **kw):
        seen.append(complete)
        return payloads, []
    monkeypatch.setattr(te, "check_cycles", cycles)
    run(shards(tmp_path, 2))
    run(shards(tmp_path, 3)[:2])
    assert seen == [True, False]


def test_roster_refusal_makes_the_catalog_incomplete(tmp_path):
    refusal = {"kind": "roster_refusal", "system": "pp8", "participant": 0,
               "status": "refused:no_provider", "reason": "no provider declares this repeat unit"}
    catalog = run(shards(tmp_path, roster_refusals=[refusal]))
    assert catalog["complete"] is False
    assert any("refused:no_provider" in f and "pp8" in f for f in catalog["gates"]["failures"])


def test_rows_without_a_shard_header_are_incomplete(tmp_path):
    path = tmp_path / "bare.jsonl"
    path.write_text(ROWS.read_text())
    catalog = run([path])
    assert catalog["complete"] is False
    assert any("header" in f for f in catalog["gates"]["failures"])


def test_round_trip_loads_at_every_provenance_temperature(tmp_path):
    catalog = run(shards(tmp_path))
    expected = sum(1 for r in fixture_rows() if r["payload"]
                   for t in r["payload"]["templates"] if t["status"] == "emitted")
    # the fixture's two rows that emit one executable template are merged into one
    merged = catalog["gates"]["collisions"]["copies"] - catalog["gates"]["collisions"]["keys"]
    assert expected > 0 and merged >= 1
    expected -= merged
    for T in catalog["provenance"]["temperatures"]:
        assert len(te.loaded_templates(catalog["payloads"], T)) == expected
    assert catalog["gates"]["round_trip"]["templates"] == expected
    assert catalog["gates"]["round_trip"]["ok"] is True


def test_round_trip_failure_is_a_gate(tmp_path, monkeypatch):
    def broken(payloads, temperature, **kw):
        raise te.EmissionError("orphaned derived template")
    monkeypatch.setattr(te, "loaded_templates", broken)
    catalog = run(shards(tmp_path))
    assert catalog["complete"] is False
    assert any("round trip" in f for f in catalog["gates"]["failures"])


# --------------------------------------------------------------------------
# CLI end to end
# --------------------------------------------------------------------------
def test_cli_assemble_end_to_end(tmp_path, monkeypatch):
    monkeypatch.setattr(emit, "load_services", lambda database, families: Services())
    paths = [str(p) for p in shards(tmp_path, 2)]
    out_two, out_one = tmp_path / "two.json", tmp_path / "one.json"
    assert emit.main(["assemble", *paths, "--out", str(out_two)]) == 0
    rows = fixture_rows()
    one = write_rows(tmp_path / "all.jsonl", header(None, rows), rows)
    assert emit.main(["assemble", str(one), "--out", str(out_one)]) == 0
    assert out_two.read_bytes() == out_one.read_bytes()
    catalog = json.loads(out_two.read_text())
    assert catalog["schema"] == "kmc-template-catalog/1" and catalog["complete"] is True
    assert te.loaded_templates(catalog["payloads"], 700.0)


def test_cli_incomplete_catalog_is_written_and_exits_nonzero(tmp_path, monkeypatch):
    monkeypatch.setattr(emit, "load_services", lambda database, families: Services())
    paths = [str(p) for p in shards(tmp_path, 2)][:1]
    out = tmp_path / "partial.json"
    assert emit.main(["assemble", *paths, "--out", str(out)]) == 1
    catalog = json.loads(out.read_text())
    assert catalog["complete"] is False and catalog["gates"]["failures"]


def test_module_runs_as_a_script():
    import subprocess
    done = subprocess.run([sys.executable, "-m", "rmgpy.kmc.emit", "assemble", "--help"],
                          capture_output=True, text=True, timeout=120)
    assert done.returncode == 0 and "--out" in done.stdout


# --------------------------------------------------------------------------
# rework 1: the midpoint gate (computed per record at emission)
# --------------------------------------------------------------------------
def _emitted_row(rows, index=0):
    return [r for r in rows if r["status"] == "emitted"][index]


def test_midpoint_gate_aggregates_the_rows(tmp_path):
    catalog = run(shards(tmp_path))
    mid = catalog["gates"]["midpoint"]
    assert mid["evaluated"] is True and mid["threshold"] == 0.01
    assert mid["templates"] == catalog["gates"]["round_trip"]["templates"]
    assert mid["worst"] == 0.0 and mid["failing"] == [] and mid["unverified"] == []
    assert catalog["complete"] is True


def test_injected_midpoint_failure_is_incomplete_with_the_failure_listed(tmp_path):
    rows = fixture_rows()
    row = _emitted_row(rows)
    key = next(iter(row["midpoint"]))
    row["midpoint"][key] = {"status": "ok", "worst": 0.05, "worst_T": 650.0, "n": 28}
    path = write_rows(tmp_path / "m.jsonl", header(None, rows), rows)
    catalog = run([path])
    mid = catalog["gates"]["midpoint"]
    assert mid["worst"] == 0.05 and len(mid["failing"]) == 1
    assert catalog["complete"] is False
    assert any("midpoint" in f and "650" in f for f in catalog["gates"]["failures"])


def test_error_at_the_threshold_passes(tmp_path):
    rows = fixture_rows()
    row = _emitted_row(rows)
    key = next(iter(row["midpoint"]))
    row["midpoint"][key]["worst"] = 0.01
    assert run([write_rows(tmp_path / "m.jsonl", header(None, rows), rows)])["complete"] is True


def test_missing_midpoint_entry_is_unverified_and_incomplete(tmp_path):
    rows = fixture_rows()
    row = _emitted_row(rows)
    row["midpoint"].pop(next(iter(row["midpoint"])))
    catalog = run([write_rows(tmp_path / "m.jsonl", header(None, rows), rows)])
    mid = catalog["gates"]["midpoint"]
    assert mid["evaluated"] is False and len(mid["unverified"]) == 1
    assert catalog["complete"] is False
    assert any("midpoint check not evaluable" in f for f in catalog["gates"]["failures"])


def test_unavailable_midpoint_entry_is_unverified_and_incomplete(tmp_path):
    rows = fixture_rows()
    row = _emitted_row(rows)
    key = next(iter(row["midpoint"]))
    row["midpoint"][key] = {"status": "unavailable", "reason": "no emitted witness"}
    catalog = run([write_rows(tmp_path / "m.jsonl", header(None, rows), rows)])
    assert catalog["gates"]["midpoint"]["evaluated"] is False
    assert catalog["complete"] is False
    assert any("no emitted witness" in f for f in catalog["gates"]["failures"])


def test_a_derived_template_is_not_midpoint_gated(tmp_path, monkeypatch):
    def derive(payloads, counter):
        out = copy.deepcopy(payloads)
        out[0]["templates"][0]["derived"] = {"primary": {"record": "x", "key": "y"}}
        return out
    monkeypatch.setattr(te, "derive_weights", derive)
    rows = fixture_rows()
    for row in rows:
        row["midpoint"] = {}          # no entry at all for any template
    catalog = run([write_rows(tmp_path / "m.jsonl", header(None, rows), rows)])
    mid = catalog["gates"]["midpoint"]
    assert mid["unavailable"] == {"derived": 1}
    assert len(mid["unverified"]) == catalog["gates"]["round_trip"]["templates"] - 1


class FakeRates:
    """Stand-in for ``ProductionRates``: a fixed smooth k(T) on the requested grid."""
    instances = []

    def __init__(self, database, families, *, database_path=None, database_sha=None,
                 temperature_grid=None, **kw):
        self.grid = tuple(temperature_grid)
        FakeRates.instances.append(self)

    def __call__(self, request):
        import math
        from types import SimpleNamespace
        return SimpleNamespace(problems=[], temperatures=self.grid, ssa_multiplier=1.0,
                               k=[1.0e7 * math.exp(-400.0 / t) for t in self.grid])


def _one_row_with_maps():
    return json.loads((HERE / "fixtures" / "i155_row_with_maps.json").read_text())


def test_midpoint_equals_the_harness_midpoint_table(monkeypatch):
    acc = _harness()
    row = _one_row_with_maps()
    payload = row["payload"]
    record = {"event_id": payload["record_event_id"], "family": "R_Recombination",
              "template": "T", "orientation": "as_generated"}
    monkeypatch.setattr(te, "ProductionRates", FakeRates)
    monkeypatch.setitem(acc._S, "by_id", {record["event_id"]: record})
    monkeypatch.setitem(acc._S, "db", None)
    monkeypatch.setitem(acc._S, "records", [record])
    table = acc.midpoint_table([payload])
    entries = emit.midpoint_entries(record, payload, emit.RatesByGrid(None, ["R_Recombination"], "/db"))
    (entry,) = entries.values()
    assert table["templates"] == 1 and entry["status"] == "ok"
    assert entry["worst"] == table["worst"] and entry["worst_T"] == table["worst_T"]
    assert table["worst"] > 0.0


def test_midpoint_without_a_witness_is_unavailable_not_skipped(monkeypatch):
    monkeypatch.setattr(te, "ProductionRates", FakeRates)
    payload = copy.deepcopy(_one_row_with_maps()["payload"])
    for t in payload["templates"]:
        t["maps"] = {}
    record = {"event_id": payload["record_event_id"], "family": "R_Recombination",
              "template": "T", "orientation": "as_generated"}
    entries = emit.midpoint_entries(record, payload, emit.RatesByGrid(None, ["R_Recombination"], "/db"))
    (entry,) = entries.values()
    assert entry["status"] == "unavailable" and "witness" in entry["reason"]


def test_midpoint_rate_failure_is_unavailable_with_the_reason(monkeypatch):
    class Failing(FakeRates):
        def __call__(self, request):
            raise te.EmissionError("no rate")
    monkeypatch.setattr(te, "ProductionRates", Failing)
    payload = _one_row_with_maps()["payload"]
    record = {"event_id": payload["record_event_id"], "family": "R_Recombination",
              "template": "T", "orientation": "as_generated"}
    (entry,) = emit.midpoint_entries(record, payload, emit.RatesByGrid(None, [], "/db")).values()
    assert entry["status"] == "unavailable" and "no rate" in entry["reason"]


# --------------------------------------------------------------------------
# rework 1: coverage is pinned into the cycle check
# --------------------------------------------------------------------------
def test_a_missing_shard_alone_makes_the_cycle_check_partial(tmp_path, monkeypatch):
    """Every record is assigned and has a row, the artifact count matches, no
    header disagrees: the only defect is that shard 1/2 never arrived."""
    seen = []

    def cycles(payloads, counter, kc, complete=True, **kw):
        seen.append(complete)
        status = "ok" if complete else "not_evaluable"
        return payloads, [{"forward": "evt_" + "a" * 20, "reverse": "evt_" + "b" * 20,
                           "status": status, "reason": "reverse record absent"}]
    monkeypatch.setattr(te, "check_cycles", cycles)
    rows = fixture_rows()
    path = write_rows(tmp_path / "only0.jsonl", header("0/2", rows), rows)
    catalog = run([path])
    assert seen == [False]
    assert catalog["gates"]["cycles"]["by_status"] == {"not_evaluable": 1}
    assert catalog["gates"]["coverage"]["complete"] is False
    assert catalog["complete"] is False
    assert any("missing shard 1/2" in f for f in catalog["gates"]["failures"])
    full = run(shards(tmp_path, 2))
    assert seen[-1] is True and full["gates"]["cycles"]["by_status"] == {"ok": 1}
