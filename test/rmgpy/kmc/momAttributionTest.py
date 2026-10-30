"""Red tests (I-233): what the MoM <-> kMC comparison measures and how it
attributes MoM's difference.

1. kernel mode classes a pool's chains as declared species with the pool's
   family end mass and kMC's identity of the family's state (the table's
   ``declared_species``), not the deck's two caps;
2. the end-group ledger's residual and threshold are exported per time;
3. a term's fit start and each participant's exact-from length are reported
   apart (the law can be exact below where its fit starts, and one
   participant can be exact where the other is not);
4. the law-domain attribution: per MoM term and time, the signed difference
   between the moment law's propensity and the exact one, summed over the
   population balance's chains and weighted by the participant's mu0, mu1,
   mu2 (the number fraction kept as a secondary column);
5. chains past the population balance's grid are bounded or withheld, never
   counted as above a length the grid does not reach;
9. the diagnostics sidecars are bound to the producer's identity and refused
   by name when they belong to another table, temperature or geometry, or
   are partial.
"""

import json
import shutil

import numpy
import pytest

from comparisonRowsTest import by_name, mom  # noqa: F401  (the mom fixture)
from momCrossMethodTest import C, D1, DR, FIXTURES, LABELS, T1 as T1_PE, T, pe
from momKernelProducerTest import H, RADICAL, T1, TABLE3, setup, produce_kernels, terms3  # noqa: F401
from momProducerTest import LABEL, ZERO_GAS, mom_entries, write_run

from rmgpy.kmc import comparison_rows as cr
from rmgpy.kmc import kernel_rates as K
from rmgpy.kmc import kernel_reference as R
from rmgpy.kmc import mom_producer as mp
from rmgpy.molecule import Molecule


# --------------------------------------------------------------------------
# 1. classification
# --------------------------------------------------------------------------
def test_a_pool_chain_is_classed_with_its_familys_end_mass_and_kmcs_identity():
    """A declared n-butyl radical is 4 periods of the end-radical family (one
    H cap): a 2-unit chain of the radical pool, which kMC classes as the
    species.  Read with the deck's two caps (the producer declares every pool
    with them) its mass is off the pool's grid, so it stays in the polymer
    counts and in Mw."""
    butyl = Molecule(smiles="[CH2]CCC")
    m_butyl = butyl.get_molecular_weight() * 1000.0
    a = Molecule(smiles="C=C").get_molecular_weight() * 1000.0
    pools = [{"label": LABEL, "monomer_mw_g_mol": a, "end_groups": ["[H]", "[H]"]},
             {"label": RADICAL, "monomer_mw_g_mol": a, "end_groups": ["[H]", "[H]"]}]
    record = {"moment_unit": {"periods_per_declared_unit": 2, "molar_mass_g_mol": a / 2},
              "pool_family": {LABEL: "chain", RADICAL: "radical_end"},
              "pool_end_mass_g_mol": {LABEL: m_butyl - 2 * a + 1.00794, RADICAL: m_butyl - 2 * a},
              "declared_species_by_pool": {LABEL: [], RADICAL: [["butyl", 4]]}}
    legacy = mp.pool_species_candidates(pools, {"butyl": butyl}, {"butyl": m_butyl}, None)
    assert legacy == {LABEL: [], RADICAL: []}                # the defect: the radical is never a candidate
    got = mp.pool_species_candidates(pools, {"butyl": butyl}, {"butyl": m_butyl}, None, kernel_record=record)
    assert got == {LABEL: [], RADICAL: [("butyl", 2, m_butyl)]}


def test_a_declared_species_off_its_familys_mass_is_refused_by_name():
    a = Molecule(smiles="C=C").get_molecular_weight() * 1000.0
    butyl = Molecule(smiles="[CH2]CCC")
    m_butyl = butyl.get_molecular_weight() * 1000.0
    pools = [{"label": RADICAL, "monomer_mw_g_mol": a, "end_groups": ["[H]", "[H]"]}]
    record = {"moment_unit": {"periods_per_declared_unit": 2, "molar_mass_g_mol": a / 2},
              "pool_family": {RADICAL: "radical_end"}, "pool_end_mass_g_mol": {RADICAL: 2.016},
              "declared_species_by_pool": {RADICAL: [["butyl", 4]]}}
    with pytest.raises(mp.ProducerError, match="butyl"):
        mp.pool_species_candidates(pools, {"butyl": butyl}, {"butyl": m_butyl}, None, kernel_record=record)


def test_kernel_mode_records_each_pools_family_and_its_declared_states(setup):
    out = setup["tmp"] / "out"
    manifest = produce_kernels(setup, out)
    record = manifest["facts"]["kernels"]["record"]
    assert record["pool_family"] == {LABEL: "chain", RADICAL: "radical_end"}
    assert record["declared_species_by_pool"] == {LABEL: [], RADICAL: []}


# --------------------------------------------------------------------------
# 2. the end-group ledger, per time
# --------------------------------------------------------------------------
def test_the_ledger_residual_and_threshold_are_exported_per_time(setup):
    out = setup["tmp"] / "out"
    produce_kernels(setup, out)
    mu = [float(x) for x in setup["spec"].macroscopic_moments_mol]
    d = 1e-3

    def chains(t):
        return mu[0] - d * t, mu[1] - 2 * d * t, mu[2] - 30 * d * t

    def radicals(t):
        return 2 * d * t, 2 * d * t, 30 * d * t
    write_run(out, times=[0.0, 10.0, 20.0, 40.0, 45.0], pools=[(LABEL, chains, 0.0), (RADICAL, radicals, 0.0)],
              gas=ZERO_GAS, rows=[])
    produce_kernels(setup, out, mode="convert")
    diagnostics = json.loads((out / "mom_diagnostics.json").read_text())
    assert len(diagnostics) == len(mom_entries(out))
    for d_ in diagnostics:
        # I-286: kernel-mode pools book e_p mu0_p(t), so none is on the ledger's basis: it is zero, exactly
        ledger, threshold = d_["end_group_ledger_g"], d_["end_group_ledger_threshold_g"]
        assert ledger == 0.0 and threshold == 0.0
    rows = by_name(cr.mom_rows(out / "mom.json", setup["deck"]))
    t = diagnostics[-1]["time_s"]
    assert rows[(t, "mom_end_group_ledger_g")]["kind"] == "diagnostic"
    assert rows[(t, "mom_end_group_ledger_g")]["value"] == diagnostics[-1]["end_group_ledger_g"]
    assert rows[(t, "mom_end_group_ledger_threshold_g")]["value"] == diagnostics[-1]["end_group_ledger_threshold_g"]


# --------------------------------------------------------------------------
# 3. fit start and exact-from, per participant
# --------------------------------------------------------------------------
def test_the_fit_start_and_each_participants_exact_length_are_reported_apart():
    """Without homolysis the reachable floor is the initial chains' 20
    periods, where every fit starts; the toy laws are exact from 4 periods,
    so each participant's exact-from is well below the fit start."""
    out = terms3([C, T1])
    (v,) = [v for v in out["record"]["validity"] if v["kernel"] == T1]
    assert v["fit_from_periods"] == v["valid_from_periods"] == 20
    assert [(p["slot"], p["family"]) for p in v["participants"]] == [(0, "radical_end"), (1, "chain")]
    assert all(p["exact_from_periods"] <= 4 for p in v["participants"])


def test_one_participant_can_be_exact_where_the_other_is_not():
    """The chain donor given 2 sites at 4 periods (its law says 0): the
    joint fit starts at 5; the acceptor's law is exact from 4, the donor's
    from 5."""
    doc = json.loads(TABLE3.read_text())
    (term,) = doc["kernels"][T1]["templates"]
    term["law"]["S1"]["exact"]["4"] = [[2, [["radical_interior:CH", 4]]]]
    far = terms3([H, C, T1], doc=doc, far_law=True)
    (v,) = [v for v in far["record"]["validity"] if v["kernel"] == T1]
    assert v["fit_from_periods"] == 5
    exact = {p["family"]: p["exact_from_periods"] for p in v["participants"]}
    assert exact["chain"] == 5 and exact["radical_end"] <= 4


# --------------------------------------------------------------------------
# 4. the rate-weighted law-domain attribution
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def pe_mom():
    doc, deck = pe()
    table_ = K.kernel_rates(doc, deck.spec, T)
    homolysis = [k for k, v in table_["kernels"].items() if v["type"] == "homolysis"]
    kernels = homolysis + [C, D1, DR, T1_PE]
    terms = mp.kernel_input(table_, deck.spec, kernels, temperature=T, phi=1.0, volume_ratio=1.0, far_law=True,
                            **LABELS)
    model = R.KernelModel(table_, deck.spec, kernels=kernels)
    return {"table": table_, "deck": deck, "record": terms["record"], "terms": terms, "model": model}


def _state(model, chains):
    counts, over = numpy.zeros(model.size), numpy.zeros(3 * len(model.families))
    for (family, n), count in chains.items():
        counts[model.index(family, n)] += count
    return counts, over


def test_on_chains_past_every_exact_length_the_law_is_the_exact_propensity(pe_mom):
    model, record = pe_mom["model"], pe_mom["record"]
    counts, over = _state(model, {("chain", 60): 5, ("radical_end", 60): 3, ("radical_end", 66): 2,
                                  ("alkene_terminal", 64): 4})
    out = R.law_attribution(model, pe_mom["table"], record, counts, over)
    assert {t["template"] for t in out} == {law["template"] for law in record["laws"]}
    seen = 0
    for term in out:
        for slot in term["slots"]:
            for k, m in slot["moments"].items():
                if m["exact_mol_s"]:
                    seen += 1
                    assert abs(m["law_minus_exact_mol_s"]) <= 1e-9 * abs(m["exact_mol_s"]), (term["template"], k)
    assert seen > 0


def test_a_chain_below_the_homolysis_cut_length_shows_the_laws_negative_rate(pe_mom):
    """5 periods, min fragment 4: no cut exists (exact 0), the law a n + b
    gives w (5 a + b) < 0 per chain -- the spurious negative rate, signed."""
    model, terms = pe_mom["model"], pe_mom["terms"]
    (h,) = terms["homolysis"]
    counts, over = _state(model, {("chain", 5): 7})
    (term,) = [t for t in R.law_attribution(model, pe_mom["table"], pe_mom["record"], counts, over)
               if t["template"] == h["template"]]
    (slot,) = term["slots"]
    law = h["w"] * h["phi"] * (h["a"] * 5 + h["b"]) * 7 * model.u
    assert law < 0.0
    for k in range(3):
        m = slot["moments"][str(k)]
        assert m["exact_mol_s"] == 0.0
        assert m["law_minus_exact_mol_s"] == pytest.approx(law * (5 / 2) ** k, rel=1e-12)
    assert slot["mu0_fraction_below"] == 1.0


def test_the_command_line_writes_the_attribution_beside_the_number_fractions(pe_mom, tmp_path):
    manifest = {"facts": {"kernels": {"record": pe_mom["record"]}}}
    (tmp_path / "producer.json").write_text(json.dumps(manifest))
    table_path = tmp_path / "table.json"
    table_path.write_text(json.dumps(pe_mom["table"]))
    deck = FIXTURES / "i213_pe_toy_deck.py"
    kernels = sorted(pe_mom["record"]["kernels"])
    argv = [str(deck), str(table_path), "--end", "1e3", "--snapshots", "2", "--seeds", "1",
            "--validity", str(tmp_path / "producer.json"), "--out", str(tmp_path / "ref.json")]
    for k in kernels:
        argv += ["--kernel", k]
    assert R.main(argv) == 0
    doc = json.loads((tmp_path / "ref.json").read_text())
    assert len(doc["law_domain"]["attribution"]) == len(doc["times_s"])
    assert doc["law_domain"]["identity"]["table_templates_sha256"] == pe_mom["record"]["table_templates_sha256"]


# --------------------------------------------------------------------------
# 5. the grid's overflow
# --------------------------------------------------------------------------
def test_a_length_past_the_grid_withholds_the_number_fraction(pe_mom):
    model = pe_mom["model"]
    record = json.loads(json.dumps(pe_mom["record"]))
    (v,) = [v for v in record["validity"] if v["kernel"] == T1_PE]
    v["participants"][1]["exact_from_periods"] = model.n_max + 5
    counts, over = _state(model, {("chain", 60): 5, ("radical_end", 60): 1})
    over[3 * model.families.index("chain"):3 * model.families.index("chain") + 3] = [2.0, 2.0 * 500, 2.0 * 500 ** 2]
    fractions = R.law_domain_fractions(model, counts, over, record["validity"])
    (f,) = [f for f in fractions if f["template"] == v["template"] and f["family"] == "chain"]
    assert f["mu0_fraction_below"] is None and "grid" in f["withheld"]


def test_the_overflow_bounds_the_moments_its_own_moments_cannot_carry(pe_mom):
    """Chains past the grid: their mu0..mu2 are known, mu3 is not; the
    exact homolysis rate weighted by n^2 needs n^3 there, bounded by
    o2^2 / o1 below and (every period of the system) o2 above -- the value
    is withheld and its bounds bracket the truth (here two chains of 500)."""
    model, terms = pe_mom["model"], pe_mom["terms"]
    (h,) = terms["homolysis"]
    counts, over = _state(model, {("chain", 60): 5})
    j = 3 * model.families.index("chain")
    over[j:j + 3] = [2.0, 2.0 * 500, 2.0 * 500 ** 2]
    (term,) = [t for t in R.law_attribution(model, pe_mom["table"], pe_mom["record"], counts, over)
               if t["template"] == h["template"]]
    (slot,) = term["slots"]
    n_grid = 5 * (h["a"] * 60 + h["b"])
    n_over = h["a"] * 2.0 * 500 + h["b"] * 2.0
    assert slot["moments"]["0"]["exact_mol_s"] == pytest.approx(h["w"] * (n_grid + n_over) * model.u, rel=1e-12)
    assert slot["moments"]["0"]["law_minus_exact_mol_s"] == pytest.approx(0.0, abs=1e-30)
    m2 = slot["moments"]["2"]
    assert m2["exact_mol_s"] is None and "bounded" in slot["withheld"]["2"]
    truth = h["w"] * (5 * (h["a"] * 60 + h["b"]) * 60 ** 2 + 2 * (h["a"] * 500 + h["b"]) * 500 ** 2) / 4 * model.u
    lo, hi = m2["exact_bounds_mol_s"]
    assert lo * (1 - 1e-12) <= truth <= hi * (1 + 1e-12)
    assert m2["law_minus_exact_bounds_mol_s"] == pytest.approx([0.0, 0.0], abs=1e-12 * hi)


# --------------------------------------------------------------------------
# 9. provenance
# --------------------------------------------------------------------------
@pytest.mark.parametrize("edit,word", [
    (lambda r: r.update(temperature_K=r["temperature_K"] + 1.0), "temperature"),
    (lambda r: r["moment_unit"].update(periods_per_declared_unit=3), "geometry"),
    (lambda r: r.update(chemistry_sha256="0" * 64), "chemistry"),
    (lambda r: r.update(table_templates_sha256="0" * 64), "template rows"),
    (lambda r: r.pop("laws"), "no MoM law record"),
    (lambda r: r["laws"].append(dict(r["laws"][0], template="nope:nope", templates=["nope:nope"])), "nope:nope"),
])
def test_a_law_record_of_another_identity_is_refused_by_name(pe_mom, edit, word):
    record = json.loads(json.dumps(pe_mom["record"]))
    edit(record)
    with pytest.raises(R.KernelReferenceError, match=word):
        R.bind_mom_kernels(pe_mom["model"], pe_mom["table"], record)


def _copy_mom(mom_, tmp_path):
    setup_, path, _run = mom_
    for name in ("mom.json", "conventions.json", "mom_diagnostics.json"):
        shutil.copy(path.parent / name, tmp_path / name)
    return setup_, tmp_path / "mom.json"


def test_the_mom_diagnostics_are_bound_to_their_mom_file(mom, tmp_path):  # noqa: F811
    setup_, path = _copy_mom(mom, tmp_path)
    header = json.loads(path.read_text())["header"]
    assert header["diagnostics_sha256"]
    cr.mom_rows(path, setup_["deck"])
    diagnostics = json.loads((tmp_path / "mom_diagnostics.json").read_text())
    diagnostics[0]["whole_withheld"] = {"whole_mu2_mol": "edited"}
    (tmp_path / "mom_diagnostics.json").write_text(json.dumps(diagnostics))
    with pytest.raises(cr.ComparisonRowsError, match="mom_diagnostics.json"):
        cr.mom_rows(path, setup_["deck"])


def test_a_partial_mom_diagnostics_file_is_refused(mom, tmp_path):  # noqa: F811
    setup_, path = _copy_mom(mom, tmp_path)
    doc = json.loads(path.read_text())
    diagnostics = json.loads((tmp_path / "mom_diagnostics.json").read_text())[:-1]
    raw = json.dumps(diagnostics).encode()
    (tmp_path / "mom_diagnostics.json").write_bytes(raw)
    doc["header"]["diagnostics_sha256"] = cr._sha256(raw)
    path.write_text(json.dumps(doc))
    with pytest.raises(cr.ComparisonRowsError, match="no diagnostics at"):
        cr.mom_rows(path, setup_["deck"])
