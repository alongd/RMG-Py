"""Red tests (I-217, S7 of the shared-chemistry design): the comparison's
exact inventory rows, mapped identically from a kMC run, a kernel
reference (PBE, SSA) and a MoM ``--mom`` file, and the join of two row files.

The kMC run is a seeded toy PE melt whose atoms cross classes (C-H
homolysis of the chain ends makes the built-in species H), with short chains
in the explicit ladder from the start.
"""

import copy
import json
import shutil
from fractions import Fraction
from pathlib import Path

import pytest

import melt_engine_support as S
from kernelReferenceTest import RECOMBINATION, HOMOLYSIS, scission_table, toy_table
from meltRunTest import PE_LIVE, lines, toy_payloads, write_catalog, write_deck
from momContractTest import produced, run_cli

from rmgpy.kmc import comparison_rows as cr
from rmgpy.kmc import kernel_rates as kr_
from rmgpy.kmc import kernel_reference as R
from rmgpy.kmc import melt_observables as mo
from rmgpy.kmc import run as kr

T = 700.0
TIMES = [0.0, 1.5, 3.0]
SEEDS = [1, 2, 3]
COUNTS = {3: 4, 10: 3}          # DP 3 chains sit in the ladder (cutoff 4) from the start


def toy_run(directory):
    catalog = write_catalog(directory / "pe_h_loss.json",
                            toy_payloads([t for t in S.pe_templates() if t.name in PE_LIVE]
                                         + S.h_loss_templates()),
                            tag="pe-h-loss", bind=S.pe_mapping({1: 1}))
    deck = write_deck(directory / "deck.py", S.pe_mapping(COUNTS), cutoff=4)
    out = directory / "run"
    kr.run(deck, [catalog], T, TIMES, SEEDS, out)
    return out


@pytest.fixture(scope="module")
def kmc(tmp_path_factory):
    return toy_run(tmp_path_factory.mktemp("kmc"))


def observations(run_dir):
    manifest = json.loads((run_dir / "run.json").read_text())
    return {r["seed"]: lines(run_dir / r["path"]) for r in manifest["runs"]}


def by_name(doc):
    out = {}
    for row in doc["rows"]:
        out[(row["time_s"], row["row"])] = row
    return out


# --------------------------------------------------------------------------
# the kMC observation carries the inputs of the rows
# --------------------------------------------------------------------------
def test_the_observation_partitions_every_atom_by_class(kmc):
    moved = False
    for seed, obs in observations(kmc).items():
        first = obs[0]["tier1"]["species"]["H"]["composition"]
        for o in obs:
            t1 = o["tier1"]
            total = {}
            parts = [t1["polymer_derived"]["composition"], t1["undeclared_light"]["composition"]]
            parts += [s["composition"] for s in t1["species"].values()]
            for part in parts:
                for el, k in part.items():
                    total[el] = total.get(el, 0) + k
            assert {k: v for k, v in total.items() if v} == o["invariants"]["composition"]
            moved = moved or t1["species"]["H"]["composition"] != first
            dist = o["tier2"]["backbone_distribution"]
            assert sum(dist.values()) == t1["polymer_derived"]["molecules"]
            b = o["backbone_atoms_per_unit"]
            assert [sum(k for k in dist.values()), sum((int(x) // b) * k for x, k in dist.items()),
                    sum((int(x) // b) ** 2 * k for x, k in dist.items())] == o["tier2"]["dp_moments_representative"]
    assert moved, "the run must move atoms between classes (C-H homolysis makes H)"


# --------------------------------------------------------------------------
# element inventory: exact zero residual
# --------------------------------------------------------------------------
def test_the_kmc_element_residual_is_exactly_zero(kmc):
    doc = cr.kmc_rows(kmc)
    assert doc["schema"] == cr.SCHEMA and doc["producer"]["kind"] == "kmc"
    rows = by_name(doc)
    for t in TIMES:
        for el in ("C", "H"):
            r = rows[(t, f"element_residual_mol.{el}")]
            assert r["kind"] == "exact" and r["value"] == 0.0 and r["producer"] == "kmc"
            assert rows[(t, f"element_mol.{el}")]["value"] == rows[(0.0, f"element_mol.{el}")]["value"]
    # the species row moved: the inventory really spans classes
    assert rows[(TIMES[-1], "species_mol.H")]["value"] > 0.0


def test_every_row_names_its_producer_and_kind(kmc):
    doc = cr.kmc_rows(kmc)
    for row in doc["rows"]:
        assert row["producer"] == "kmc" and row["kind"] in cr.KINDS
        if row["kind"] == "withheld":
            assert row["reason"]
        if row["kind"] == "interval":
            assert row["lower"] <= row["upper"]
    names = {r["row"] for r in doc["rows"]}
    assert {"total_mass_g", "mass_outside_species_g", "species_mol.H", "species_mol.H2",
            *(f"ladder_mol.{p}" for p in range(1, 9)), *mo.WHOLE_ENSEMBLE_KEYS} <= names


# --------------------------------------------------------------------------
# intervals
# --------------------------------------------------------------------------
def test_a_length_off_the_period_grid_is_an_interval():
    assert cr.length_interval(6, 1) == (6, 6)
    assert cr.length_interval(6, 2) == (3, 3)
    assert cr.length_interval(7, 2) == (3, 4)


def strip_new_fields(run_dir, target):
    """A copy of the run whose observations predate the class compositions
    and the backbone distribution (observables schema 6)."""
    shutil.copytree(run_dir, target)
    manifest = json.loads((target / "run.json").read_text())
    for r in manifest["runs"]:
        path = target / r["path"]
        out = []
        for line in path.read_text().splitlines():
            o = json.loads(line)
            if "record" not in o:
                o["schema"] = "kmc-melt-observables/6"
                del o["tier2"]["backbone_distribution"]
                for part in ("polymer_derived", "undeclared_light"):
                    del o["tier1"][part]["composition"]
                for s in o["tier1"]["species"].values():
                    del s["composition"]
            out.append(json.dumps(o, sort_keys=True))
        path.write_text("\n".join(out) + "\n")
    return target


def test_interval_rows_bracket_the_kmc_value(kmc, tmp_path):
    exact = by_name(cr.kmc_rows(kmc))
    old = cr.kmc_rows(strip_new_fields(kmc, tmp_path / "old"))
    intervals = 0
    for row in old["rows"]:
        if not row["row"].startswith("ladder_mol."):
            continue
        value = exact[(row["time_s"], row["row"])]
        assert value["kind"] == "exact"
        if row["kind"] == "interval":
            intervals += 1
            assert row["lower"] <= value["value"] <= row["upper"], (row, value)
            assert row["value"] is None
        else:
            assert row["kind"] == "exact" and row["value"] == value["value"]
    assert intervals, "the run must give the ladder a molecule whose period length the old record cannot pin"
    # the old record has no length distribution above the ladder and no class split
    assert old["withheld_families"]["mwd_mol"]
    assert any(r["row"] == "element_residual_mol.H" and r["value"] == 0.0 for r in old["rows"])


# --------------------------------------------------------------------------
# kernel references
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def reference(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("ref")
    table = tmp / "toy.kernel_rates.json"
    table.write_text(json.dumps(toy_table(tmp)))
    deck = write_deck(tmp / "deck.py", S.pe_mapping({10: 60}), cutoff=4)
    code = R.main([str(deck), str(table), "--end", "0.2", "--snapshots", "3", "--seeds", "4", "--n-max", "20",
                   "--kernel", HOMOLYSIS, "--kernel", RECOMBINATION])
    assert code == 0
    return tmp / "toy.kernel_rates.reference.json", deck


def test_the_population_balance_rows_carry_the_law_domain_diagnostic(tmp_path):
    """Given the MoM kernels' validity lengths, the balance's rows report,
    per time, the fraction of each family's chains below each (diagnostic,
    never compared); without them, no such row."""
    from rmgpy.kmc import mom_producer as mp
    doc = toy_table(tmp_path)
    table = tmp_path / "toy.kernel_rates.json"
    table.write_text(json.dumps(doc))
    deck = write_deck(tmp_path / "deck.py", S.pe_mapping({10: 60}), cutoff=4)
    terms = mp.kernel_input(doc, kr.read_deck(deck).spec, [HOMOLYSIS, RECOMBINATION], temperature=T, phi=1.0,
                            label="P", radical_label="Re", volume_ratio=1.0)
    producer = tmp_path / "producer.json"
    producer.write_text(json.dumps({"facts": {"kernels": {"record": terms["record"]}}}))
    assert R.main([str(deck), str(table), "--end", "0.2", "--snapshots", "3", "--seeds", "2", "--n-max", "20",
                   "--kernel", HOMOLYSIS, "--kernel", RECOMBINATION, "--validity", str(producer)]) == 0
    path = tmp_path / "toy.kernel_rates.reference.json"
    (h,), (c,) = terms["homolysis"], terms["recombination"]
    rows = [r for r in cr.reference_rows(path, "pbe")["rows"] if r["row"].startswith("law_domain")]
    t0 = [r for r in rows if r["time_s"] == 0.0]
    assert {r["row"] for r in t0} == {f"law_domain_mu0_fraction.{h['template']}@chain.S0"}   # no radicals yet
    assert t0[0]["kind"] == "diagnostic" and t0[0]["value"] == 0.0          # 60-unit chains, all long
    later = [r for r in rows if r["time_s"] > 0.0]
    assert {r["row"] for r in later} == {f"law_domain_mu0_fraction.{h['template']}@chain.S0",
                                         f"law_domain_mu0_fraction.{c['template']}@radical_end.S0",
                                         f"law_domain_mu0_fraction.{c['template']}@radical_end.S1"}
    assert all(0.0 <= r["value"] <= 1.0 for r in later)
    attribution = {r["row"] for r in cr.reference_rows(path, "pbe")["rows"] if r["row"].startswith("law_attr")}
    assert f"law_attribution_mol_s.{h['template']}.S0.mu1" in attribution
    assert not any(r["row"].startswith(("law_domain", "law_attr")) for r in cr.reference_rows(path, "ssa")["rows"])


def test_the_reference_records_its_binding_and_row_observations(reference):
    path, deck = reference
    doc = json.loads(path.read_text())
    assert doc["schema"] == R.SCHEMA
    assert doc["deck"]["sha256"] == kr.read_deck(deck).sha256
    assert len(doc["row_observations"]["pbe"]) == 3 and len(doc["row_observations"]["ssa"]["runs"]) == 4
    # joins outgrow the 20-period grid: the overflow's atoms are part of the inventory
    assert all(run[-1]["overflow_molecules"] > 0 for run in doc["row_observations"]["ssa"]["runs"])


@pytest.mark.parametrize("producer", ["pbe", "ssa"])
def test_the_reference_element_residual_is_exactly_zero(reference, producer):
    path, _deck = reference
    doc = cr.reference_rows(path, producer)
    assert doc["producer"]["kind"] == f"kernel_{producer}"
    rows = by_name(doc)
    times = doc["binding"]["times_s"]
    for t in times:
        for el in ("C", "H"):
            r = rows[(t, f"element_residual_mol.{el}")]
            assert r["kind"] == "exact" and r["value"] == 0.0
    if producer == "pbe":
        assert all(r["kind"] == "diagnostic" for r in doc["rows"] if r["row"].startswith("element_drift_mol."))
    # the run moved: radicals and other lengths appear
    assert rows[(times[-1], "whole_mu0_mol")]["value"] != rows[(times[0], "whole_mu0_mol")]["value"]


def test_an_unbalanced_law_is_named_and_its_residual_is_not_zero():
    """The random-scission toy cuts a chain into two capped chains: it makes
    H atoms.  The model names the term; both references report the residual."""
    spec = S.pe_spec({10: 40})
    model = R.KernelModel(scission_table(spec, 0.7), spec)
    assert [u["template"] for u in model.unbalanced] == ["scission"]
    times = [0.0, 0.5]
    pbe = R.population_balance(model, times)
    ssa, _events = R.stochastic(model, times, 0)
    h0 = model.row_observation(model.initial)["elements"]["H"]
    assert pbe[-1]["row"]["elements"]["H"] > h0 and ssa[-1]["row"]["elements"]["H"] > h0


def test_the_reference_counts_elements_per_state_from_the_table(reference):
    table = toy_table(Path(reference[0]).parent)
    assert table["schema"] == kr_.SCHEMA
    assert table["moment_unit"]["composition"] == {"C": "1", "H": "2"}
    model = R.KernelModel(table, S.pe_spec({10: 5}))
    assert model.composition("chain", 7) == {"C": Fraction(7), "H": Fraction(16)}
    assert model.composition("radical_end", 7) == {"C": Fraction(7), "H": Fraction(15)}
    assert model.unbalanced == []


# --------------------------------------------------------------------------
# MoM
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def mom(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("mom")
    setup, path = produced(tmp)
    assert run_cli(setup, path, tmp / "run") == 0
    return setup, path, tmp / "run"


def test_the_mom_rows_report_moments_and_withhold_what_the_file_lacks(mom):
    setup, path, _run = mom
    doc = cr.mom_rows(path, setup["deck"])
    assert doc["producer"]["kind"] == "mom" and doc["binding"]["pass"] == cr.PASS_WELL_MIXED
    for family in ("element_mol", "element_residual_mol", "mwd_mol"):
        assert doc["withheld_families"][family]
    rows = by_name(doc)
    t0 = doc["binding"]["times_s"][0]
    assert rows[(t0, "whole_mu0_mol")]["kind"] == "exact"
    for p in range(1, 4):
        assert rows[(t0, f"ladder_mol.{p}")]["kind"] == "withheld"


def test_a_withheld_mom_row_carries_the_producers_reason(mom, tmp_path):
    """A whole-ensemble row the producer withheld at a time says why (its
    mom_diagnostics.json), not only that the file lacks it."""
    setup, path, _run = mom
    doc = json.loads(path.read_text())
    t = doc["entries"][1]["time_s"]
    doc["entries"][1]["quantities"].pop("whole_mu2_mol", None)
    shutil.copy(path.parent / "conventions.json", tmp_path / "conventions.json")
    why = "the end-group ledger books 1 g/mol x mol"
    raw = json.dumps([{"time_s": e["time_s"], "whole_withheld": {"whole_mu2_mol": why} if e["time_s"] == t else {}}
                      for e in doc["entries"]]).encode()
    (tmp_path / "mom_diagnostics.json").write_bytes(raw)
    doc["header"]["diagnostics_sha256"] = cr._sha256(raw)
    (tmp_path / "mom.json").write_text(json.dumps(doc))
    rows = by_name(cr.mom_rows(tmp_path / "mom.json", setup["deck"]))
    assert rows[(t, "whole_mu2_mol")]["kind"] == "withheld"
    assert rows[(t, "whole_mu2_mol")]["reason"] == f"withheld by the MoM producer: {why}"


def test_a_mom_file_of_another_deck_is_refused(mom, tmp_path):
    setup, path, _run = mom
    doc = json.loads(path.read_text())
    doc["header"]["producer"]["deck_sha256"] = "0" * 64
    other = tmp_path / "mom.json"
    other.write_text(json.dumps(doc))
    shutil.copy(path.parent / "conventions.json", tmp_path / "conventions.json")
    with pytest.raises(cr.ComparisonRowsError, match="deck"):
        cr.mom_rows(other, setup["deck"])


# --------------------------------------------------------------------------
# the runner and the join
# --------------------------------------------------------------------------
def test_the_runner_writes_the_rows_beside_run_json(kmc):
    manifest = json.loads((kmc / "run.json").read_text())
    assert manifest["rows"] == "rows.json"
    doc = json.loads((kmc / "rows.json").read_text())
    assert doc == json.loads(json.dumps(cr.kmc_rows(kmc)))
    assert doc["binding"]["deck_sha256"] == manifest["deck"]["sha256"]
    assert doc["binding"]["times_s"] == manifest["times_s"] and doc["binding"]["pass"] == cr.PASS_DECLARED


def test_the_join_compares_exact_rows_and_marks_withheld_ones(mom):
    setup, path, run = mom
    a = json.loads((run / "rows.json").read_text())
    b = cr.mom_rows(path, setup["deck"])
    b["binding"]["pass"] = a["binding"]["pass"]           # the comparison under test is kMC vs MoM rows
    table = cr.join(a, b)
    rows = {(r["time_s"], r["row"]): r for r in table}
    t0 = a["binding"]["times_s"][0]
    whole = rows[(t0, "whole_mu0_mol")]
    assert whole["status"] == "compared" and whole["difference"] == whole["a"] - whole["b"]
    assert whole["source_a"] == "kmc" and whole["source_b"] == "mom"
    for (t, name), r in rows.items():
        if name.startswith(("element_mol.", "element_residual_mol.", "ladder_mol.")):
            assert r["status"] == "withheld_b" and r["difference"] is None and r["reason"], r


def test_a_row_withheld_on_either_side_is_never_differenced(kmc):
    a = cr.kmc_rows(kmc)
    b = copy.deepcopy(a)
    target = next(r for r in b["rows"] if r["row"] == "total_mass_g" and r["time_s"] == TIMES[-1])
    target.update(kind="withheld", value=None, stderr=None, reason="test")
    rows = {(r["time_s"], r["row"]): r for r in cr.join(a, b)}
    r = rows[(TIMES[-1], "total_mass_g")]
    assert r["status"] == "withheld_b" and r["difference"] is None and r["reason"] == "test"
    same = rows[(TIMES[0], "total_mass_g")]
    assert same["status"] == "compared" and same["difference"] == 0.0


@pytest.mark.parametrize("field,value,word", [("deck_sha256", "0" * 64, "deck"),
                                              ("temperature_K", 800.0, "temperature"),
                                              ("times_s", [0.0, 1.0, 3.0], "time grid"),
                                              ("pass", cr.PASS_WELL_MIXED, "pass")])
def test_the_join_refuses_a_mismatch_by_name(kmc, field, value, word):
    a = cr.kmc_rows(kmc)
    b = copy.deepcopy(a)
    b["binding"][field] = value
    with pytest.raises(cr.ComparisonRowsError, match=word):
        cr.join(a, b)


def test_the_join_command_line(kmc, tmp_path, capsys):
    a = tmp_path / "a.json"
    a.write_text(json.dumps(cr.kmc_rows(kmc)))
    b = tmp_path / "b.json"
    other = cr.kmc_rows(kmc)
    other["binding"]["temperature_K"] = 800.0
    b.write_text(json.dumps(other))
    assert cr.main(["join", str(a), str(a), "--out", str(tmp_path / "t.csv")]) == 0
    text = (tmp_path / "t.csv").read_text()
    assert "element_residual_mol.H" in text and "compared" in text
    assert cr.main(["join", str(a), str(b)]) == 2
    assert "temperature" in capsys.readouterr().err
    assert cr.main(["kmc", str(kmc), "--out", str(tmp_path / "k.json")]) == 0
    assert json.loads((tmp_path / "k.json").read_text())["schema"] == cr.SCHEMA


def test_the_module_names_no_polymer():
    text = Path(cr.__file__).read_text().lower()
    for word in ("polyethylene", "polystyrene", "polypropylene", "ethylene", "styrene", "ch2"):
        assert word not in text
