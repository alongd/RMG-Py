"""The MoM producer's kernel mode (slice S5 a + b of the shared-chemistry
design): a kernel-only MoM input written from a kernel-rates table.

The table fixture is the toy PE catalog's (``meltRunTest.pe_catalog_file``)
at 700 K from ``rmgpy.kmc.kernel_rates`` at kmc-i188 db7660c23
(``kmc-kernel-rates/3``; homolysis a = 1, b = -5 in CH2 periods, fragments
3..n - 3; end recombination per feature 1, zero at 2 periods).  Families are
off; the homolysis and end-recombination kernels carry the table's
coefficients, per template, onto a dead-chain pool and an explicit
end-radical pool.
"""

import ast
import json
from pathlib import Path

import pytest

from momProducerTest import (BASE, CHECKOUT, LABEL, ZERO_GAS, constant_run, make_setup, mom_entries, spec_of,
                             write_run)

from rmgpy.kmc import mom_producer as mp
from rmgpy.kmc import run as kr

TABLE = Path(__file__).parent / "fixtures" / "i219_toy_kernel_rates_700K.json"
H = "homolysis:substrate=chain"
C = "recombination:radical=radical_end:CH2+radical=radical_end:CH2"
RADICAL = LABEL + "_re"


def table():
    return json.loads(TABLE.read_text())


def row(doc, kernel):
    (name,) = [t["template"] for t in doc["kernels"][kernel]["templates"]]
    return next(r for r in doc["templates"] if r["template"] == name)


def terms(doc=None, kernels=(H, C), phi=1.0, counts=BASE, ratio=1.0):
    return mp.kernel_input(table() if doc is None else doc, spec_of(counts), list(kernels), temperature=700.0,
                           phi=phi, label=LABEL, radical_label=RADICAL, volume_ratio=ratio)


def produce_kernels(setup, out, **kw):
    args = dict(checkout=CHECKOUT, temperature=700.0, end=40.0, snapshots=9, families=(), catalogs=(),
                out=out, database=setup["db"], kernel_rates=str(TABLE), kernels=[H, C])
    args.update(kw)
    return mp.produce(setup["deck"], **args)


def polymers(out):
    tree = ast.parse((out / "rmg" / "input.py").read_text())
    calls = [s.value for s in tree.body if isinstance(s, ast.Expr) and isinstance(s.value, ast.Call)
             and getattr(s.value.func, "id", None) == "polymer"]
    out = {}
    for call in calls:
        kwargs = {k.arg: ast.literal_eval(k.value) for k in call.keywords}
        out[kwargs["label"]] = kwargs
    return out


@pytest.fixture
def setup(tmp_path):
    return make_setup(tmp_path)


# --------------------------------------------------------------------------
# the terms
# --------------------------------------------------------------------------
def test_the_terms_are_the_tables_coefficients():
    doc = table()
    out = terms(doc)
    h, c = row(doc, H), row(doc, C)
    assert out["homolysis"] == [{"template": h["template"], "w": h["weight"], "phi": 1.0, "a": 1.0, "b": -5.0,
                                 "min_fragment": 3, "periods_per_unit": 2, "daughter": RADICAL,
                                 "valid_from_periods": 5}]
    assert out["recombination"] == [{"template": c["template"], "k_pair": c["weight"] / c["divisor"],
                                     "sites_per_chain": 1.0, "periods_per_unit": 2, "product": LABEL,
                                     "valid_from_periods": 3}]


def test_phi_thins_the_radical_creating_templates_only():
    out = terms(phi=0.25)
    assert out["homolysis"][0]["phi"] == 0.25
    assert out["recombination"][0]["k_pair"] == terms()["recombination"][0]["k_pair"]


def test_the_pair_coefficient_is_converted_to_moms_volume():
    """MoM's concentration is n / V_poly, kMC's n / V: k_pair = (w / divisor) V_poly / V."""
    c = row(table(), C)
    out = terms(ratio=0.97)
    assert out["recombination"][0]["k_pair"] == pytest.approx(c["weight"] / c["divisor"] * 0.97, rel=1e-15)
    assert out["record"]["volume_ratio_mom_over_kmc"] == 0.97


def test_a_single_kernel_is_written_alone():
    assert terms(kernels=[C])["homolysis"] == [] and len(terms(kernels=[C])["recombination"]) == 1
    assert terms(kernels=[H])["recombination"] == [] and len(terms(kernels=[H])["homolysis"]) == 1


# --------------------------------------------------------------------------
# refusals
# --------------------------------------------------------------------------
def test_a_kernel_moms_kernels_do_not_carry_is_refused():
    with pytest.raises(mp.ProducerError, match="no MoM kernel"):
        terms(kernels=[H, C, "beta_scission:substrate=radical_interior:CH"])


def test_an_unknown_kernel_is_refused():
    with pytest.raises(mp.ProducerError, match="not in the table"):
        terms(kernels=["homolysis:substrate=nothing"])


def test_a_table_of_another_chemistry_or_temperature_or_schema_is_refused():
    doc = table()
    doc["catalog"]["chemistry_sha256"] = "0" * 64
    with pytest.raises(mp.ProducerError, match="chemistry digest"):
        terms(doc)
    with pytest.raises(mp.ProducerError, match="temperature"):
        mp.kernel_input(table(), spec_of(BASE), [H, C], temperature=701.0, phi=1.0, label=LABEL,
                        radical_label=RADICAL, volume_ratio=1.0)
    doc = table()
    doc["schema"] = "kmc-kernel-rates/1"
    with pytest.raises(mp.ProducerError, match="kmc-kernel-rates"):
        terms(doc)


def test_a_homolysis_law_that_is_not_the_uniform_cut_is_refused():
    doc = table()
    law = next(t for t in doc["kernels"][H]["templates"])["law"]
    outcomes = law["exact"]["12"]
    outcomes[0][0] += 1                     # one position counted twice
    with pytest.raises(mp.ProducerError, match="not the uniform cut"):
        terms(doc)


def test_sites_off_the_cut_positions_are_refused():
    doc = table()
    for r in doc["templates"]:
        if r["template"] == row(doc, H)["template"]:
            r["participants"][0]["chain"]["b"] = -4.0
    with pytest.raises(mp.ProducerError, match="b = a"):
        terms(doc)


def test_a_recombination_profile_on_a_fragment_length_is_refused():
    """Zero sites on a radical length homolysis makes: per_feature is then not
    every reachable radical's count."""
    doc = table()
    for r in doc["templates"]:
        if r["template"] == row(doc, C)["template"]:
            for p in r["participants"]:
                p["radical_end:CH2"]["profile"]["fragment_backbone_atoms=4"] = [0]
    with pytest.raises(mp.ProducerError, match="profile"):
        terms(doc)


def test_initial_chains_shorter_than_the_sites_allow_are_refused():
    """a*n + b < 0 at n = 4 periods (2 units): the affine sites are not the count."""
    with pytest.raises(mp.ProducerError, match="shorter than"):
        terms(counts={2: 1, 20: 1})


# --------------------------------------------------------------------------
# the input and its run
# --------------------------------------------------------------------------
def test_kernel_mode_writes_a_kernel_only_input(setup):
    out = setup["tmp"] / "out"
    manifest = produce_kernels(setup, out)
    text = (out / "rmg" / "input.py").read_text()
    assert "kineticsFamilies='none'" in text
    pools = polymers(out)
    assert set(pools) == {LABEL, RADICAL}
    record = manifest["facts"]["kernels"]
    assert pools[LABEL]["kernel_homolysis"] == record["homolysis"]
    assert pools[RADICAL]["kernel_recombination"] == record["recombination"]
    assert pools[RADICAL]["end_radical_site"] == "primary"
    assert pools[RADICAL]["moments"] == [0.0, 0.0, 0.0] and pools[RADICAL]["initial_mass"] == 0.0
    assert pools[RADICAL]["end_groups"] == pools[LABEL]["end_groups"]
    assert record["record"]["volume_ratio_mom_over_kmc"] == manifest["facts"]["volume_ratio_mom_over_kmc"]
    assert record["record"]["table_sha256"] == mp._sha256(TABLE.read_bytes())
    names = {c["name"] for c in json.loads((out / "conventions.json").read_text())["conventions"]}
    assert {"kernels", "kernel_phi", "kernel_volume"} <= names


def test_kernel_mode_writes_no_families(setup):
    with pytest.raises(mp.ProducerError, match="kernel mode writes no families"):
        produce_kernels(setup, setup["tmp"] / "out", families=["toy_homolysis"])


def test_kernel_mode_converts_a_moving_run_without_a_family_census(setup):
    out = setup["tmp"] / "out"
    produce_kernels(setup, out)
    mu = [float(x) for x in setup["spec"].macroscopic_moments_mol]
    d = 1e-3

    def chains(t):
        return mu[0] - d * t, mu[1] - 3 * d * t, mu[2] - 40 * d * t

    def radicals(t):
        return 2 * d * t, 3 * d * t, 40 * d * t
    write_run(out, times=[0.0, 10.0, 20.0, 40.0, 45.0], pools=[(LABEL, chains, 0.0), (RADICAL, radicals, 0.0)],
              gas=ZERO_GAS, rows=[])
    manifest = produce_kernels(setup, out, mode="convert")
    assert manifest["status"] == "converted"
    assert "inactive_chemistry" not in manifest["flags"]
    assert manifest["census"]["mode"] == "kernels" and manifest["census"]["profile_changes"]
    assert len(mom_entries(out)) == len(kr.observation_times(end=40.0, snapshots=9))


def test_kernel_mode_refuses_a_flat_run_as_inactive(setup):
    out = setup["tmp"] / "out"
    produce_kernels(setup, out)
    constant_run(out, setup["spec"], flat=True, rows=[])
    with pytest.raises(mp.ProducerError, match="inactive"):
        produce_kernels(setup, out, mode="convert")


# --------------------------------------------------------------------------
# pair transfer (S5 c + e): end + end disproportionation and its reverse
# --------------------------------------------------------------------------
# The I-221 toy table: the toy catalog plus declared toy templates (a propyl
# donor for D1, end-end reverse disproportionation on a terminal alkene and on
# a propene), from rmgpy.kmc.kernel_rates at kmc-i188 e43ae940f; generator
# /home/alon/runs/i221/toy/make_toy.py.  Each kernel is two templates whose
# sites are uniform only together (as PE's are).
TABLE2 = Path(__file__).parent / "fixtures" / "i221_toy_kernel_rates_700K.json"
D1 = "disproportionation:acceptor=radical_end:CH2+donor=radical_end:CH2"
DR = "reverse_disproportionation:alkene=alkene_terminal:CH2=CH+donor=chain"
ALKENE = LABEL + "_ut"


def table2():
    return json.loads(TABLE2.read_text())


def terms2(kernels, *, doc=None, phi=1.0, counts=BASE, ratio=1.0):
    return mp.kernel_input(table2() if doc is None else doc, spec_of(counts), list(kernels), temperature=700.0,
                           phi=phi, label=LABEL, radical_label=RADICAL, alkene_label=ALKENE, volume_ratio=ratio)


def names(doc, kernel):
    return sorted(t["template"] for t in doc["kernels"][kernel]["templates"])


def pair_coefficient(doc, kernel):
    """w / divisor of the kernel's templates (one value for the toy's pairs)."""
    rows = {r["template"]: r for r in doc["templates"]}
    (k,) = {round(rows[n]["weight"] / rows[n]["divisor"], 18) for n in names(doc, kernel)}
    return k


def only(doc, kernel, keep):
    """The table with ``kernel`` reduced to its templates in ``keep``."""
    doc["kernels"][kernel]["templates"] = [t for t in doc["kernels"][kernel]["templates"] if t["template"] in keep]
    return doc


def test_disproportionation_maps_its_templates_together_onto_the_radical_pool():
    """Acceptor 1 site, donor 2 sites on every length the run makes (the
    propyl donor's own template at 3 periods): one term per chain pair."""
    doc = table2()
    out = terms2([H, C, D1], ratio=0.97)
    (term,) = out["transfer"][RADICAL]
    assert term == {"template": "+".join(names(doc, D1)), "k_pair": pytest.approx(2.0 * pair_coefficient(doc, D1) * 0.97,
                                                                                   rel=1e-12),
                    "phi": 1.0, "partner": RADICAL, "product": LABEL, "partner_product": ALKENE,
                    "valid_from_periods": 3}
    assert set(out["transfer"]) == {RADICAL}
    assert out["record"]["unmapped"] == []
    assert out["record"]["checks"][term["template"]]["reachable_from_periods"] == 3


def test_reverse_disproportionation_maps_onto_the_alkene_pool_times_phi():
    doc = table2()
    out = terms2([H, C, D1, DR], phi=0.25)
    (term,) = out["transfer"][ALKENE]
    assert term == {"template": "+".join(names(doc, DR)), "k_pair": pytest.approx(6.0 * pair_coefficient(doc, DR),
                                                                                   rel=1e-12),
                    "phi": 0.25, "partner": LABEL, "product": RADICAL, "partner_product": RADICAL,
                    "valid_from_periods": 3}
    assert out["transfer"][RADICAL][0]["phi"] == 1.0          # disproportionation creates no radicals


def test_a_template_whose_sites_differ_on_a_reachable_length_is_listed_not_dropped():
    """Without the propyl donor's template the donor has no site at 3
    periods, a length homolysis makes: the kernel has no moment form."""
    doc = table2()
    main = next(n for n in names(doc, D1) if
                next(t for t in doc["kernels"][D1]["templates"] if t["template"] == n)["law"]["S1"]["far"]["outcomes"])
    out = terms2([H, C, D1], doc=only(doc, D1, {main}))
    assert out["transfer"] == {}
    (listed,) = out["record"]["unmapped"]
    assert listed["kernel"] == D1 and listed["templates"] == [main]
    assert "not one constant" in listed["reason"] and "3" in listed["reason"]


def test_the_reachable_floor_is_the_shortest_length_the_kernels_make():
    """Without homolysis no radical is shorter than the initial chains (20
    periods): the one template is uniform there."""
    doc = table2()
    main = next(n for n in names(doc, D1) if
                next(t for t in doc["kernels"][D1]["templates"] if t["template"] == n)["law"]["S1"]["far"]["outcomes"])
    out = terms2([C, D1], doc=only(doc, D1, {main}))
    (term,) = out["transfer"][RADICAL]
    assert term["template"] == main and out["record"]["unmapped"] == []
    assert out["record"]["checks"][main]["reachable_from_periods"] == 20


def test_templates_whose_products_have_no_pool_are_listed():
    """A product class MoM carries no pool for (a mid-chain radical) cannot be
    mapped: listed with its reason, and the sibling it leaves non-uniform too."""
    doc = table2()
    short = next(n for n in names(doc, DR) if
                 next(t for t in doc["kernels"][DR]["templates"] if t["template"] == n)["S0"]["per_feature"] == 0)
    term = next(t for t in doc["kernels"][DR]["templates"] if t["template"] == short)
    law = term["law"]["S1"]
    for n, outs in law["exact"].items():
        law["exact"][n] = [[m, [["radical_interior:CH", int(n)]]] for m, _p in outs]
    for side in ("first", "last"):
        for pos, outs in law["far"][side].items():
            law["far"][side][pos] = [[m, [["radical_interior", 1, 0, 0]]] for m, _p in outs]
    out = terms2([H, C, D1, DR], doc=doc)
    assert ALKENE not in out["transfer"]
    reasons = {tuple(x["templates"]): x["reason"] for x in out["record"]["unmapped"]}
    assert "radical_interior" in reasons[(short,)]
    assert any("not one constant" in r for t, r in reasons.items() if short not in t)


def test_a_kernel_on_a_class_moms_pools_do_not_carry_is_refused():
    doc = table2()
    doc["kernels"][D1]["classes"][1] = "radical_interior:CH"
    with pytest.raises(mp.ProducerError, match="radical_interior"):
        terms2([H, C, D1], doc=doc)


def test_kernel_mode_writes_the_alkene_pool_and_the_transfer_kernels(setup):
    out = setup["tmp"] / "out"
    manifest = produce_kernels(setup, out, kernel_rates=str(TABLE2), kernels=[H, C, D1, DR])
    pools = polymers(out)
    assert set(pools) == {LABEL, RADICAL, ALKENE}
    record = manifest["facts"]["kernels"]
    assert pools[RADICAL]["kernel_transfer"] == record["transfer"][RADICAL]
    assert pools[ALKENE]["kernel_transfer"] == record["transfer"][ALKENE]
    assert pools[ALKENE]["end_alkene_site"] == "primary"
    assert pools[ALKENE]["moments"] == [0.0, 0.0, 0.0] and pools[ALKENE]["initial_mass"] == 0.0
    assert pools[ALKENE]["end_groups"] == pools[LABEL]["end_groups"]
    assert manifest["facts"]["alkene_pool"] == ALKENE
    conventions = {c["name"]: c["value"] for c in json.loads((out / "conventions.json").read_text())["conventions"]}
    assert conventions["kernels"]["transfer"] == record["transfer"]
    assert conventions["kernels"]["record"]["unmapped"] == []


def test_kernel_mode_writes_no_alkene_pool_without_a_transfer_into_it(setup):
    out = setup["tmp"] / "out"
    manifest = produce_kernels(setup, out, kernel_rates=str(TABLE2), kernels=[H, C])
    assert set(polymers(out)) == {LABEL, RADICAL}
    assert manifest["facts"]["alkene_pool"] is None


def test_kernel_mode_converts_a_run_with_an_alkene_pool(setup):
    out = setup["tmp"] / "out"
    produce_kernels(setup, out, kernel_rates=str(TABLE2), kernels=[H, C, D1, DR])
    mu = [float(x) for x in setup["spec"].macroscopic_moments_mol]
    d = 1e-3

    def chains(t):
        return mu[0] - d * t, mu[1] - 3 * d * t, mu[2] - 40 * d * t

    def radicals(t):
        return d * t, 2 * d * t, 20 * d * t

    def alkenes(t):
        return d * t, d * t, 20 * d * t
    write_run(out, times=[0.0, 10.0, 20.0, 40.0, 45.0],
              pools=[(LABEL, chains, 0.0), (RADICAL, radicals, 0.0), (ALKENE, alkenes, 0.0)], gas=ZERO_GAS, rows=[])
    manifest = produce_kernels(setup, out, mode="convert", kernel_rates=str(TABLE2), kernels=[H, C, D1, DR])
    assert manifest["status"] == "converted"
    assert manifest["census"]["mode"] == "kernels" and manifest["census"]["profile_changes"]


# --------------------------------------------------------------------------
# mid-chain radicals (S5 d)
# --------------------------------------------------------------------------
# The I-223 toy table: the toy catalog with a homolysis of min fragment 4
# periods (as PE's) and declared mid-chain templates (end-mid, mid-end and
# mid-mid disproportionation, H transfer to a mid-chain radical's acceptor
# from a CH3 end, a 1,2-shift backbiting and its reverse, internal and end-mid
# reverse disproportionation), from rmgpy.kmc.kernel_rates at kmc-i188
# e43ae940f; generator /home/alon/runs/i223/toy/make_toy.py.  Every law is
# exact from 4 periods; the reverse backbiting is inert (an interior radical
# next to a chain end has no far-field site).
TABLE3 = Path(__file__).parent / "fixtures" / "i223_toy_kernel_rates_700K.json"
D2E = "disproportionation:acceptor=radical_end:CH2+donor=radical_interior:CH"
D2M = "disproportionation:acceptor=radical_interior:CH+donor=radical_end:CH2"
D3 = "disproportionation:acceptor=radical_interior:CH+donor=radical_interior:CH"
T1 = "h_transfer:acceptor=radical_end:CH2+donor=chain"
T2 = "h_transfer:acceptor=radical_interior:CH+donor=chain"
BB = "intra_h_transfer:substrate=radical_end:CH2"
BBR = "intra_h_transfer:substrate=radical_interior:CH"
DRI = "reverse_disproportionation:alkene=alkene_internal:CH=CH+donor=chain"
MIDR = LABEL + "_rm"
INTERNAL = LABEL + "_ui"
ALL3 = [H, C, D1, DR, D2E, D2M, D3, T1, T2, BB, BBR, DRI]


def table3():
    return json.loads(TABLE3.read_text())


def terms3(kernels, *, doc=None, phi=1.0, counts=BASE, ratio=1.0, far_law=False):
    return mp.kernel_input(table3() if doc is None else doc, spec_of(counts), list(kernels), temperature=700.0,
                           phi=phi, label=LABEL, radical_label=RADICAL, alkene_label=ALKENE,
                           interior_radical_label=MIDR, interior_alkene_label=INTERNAL, volume_ratio=ratio,
                           far_law=far_law)


def weight(doc, kernel):
    (name,) = names(doc, kernel)
    r = next(r for r in doc["templates"] if r["template"] == name)
    return r["weight"] / r["divisor"]


def by_template(terms):
    return {t["template"]: t for t in terms}


def test_mid_chain_disproportionation_maps_onto_the_interior_pools():
    """End-mid: the mid-chain donor's 4 sites (two beta CH2) into the
    internal-alkene pool; mid-end: the end donor's 2 into the terminal-alkene
    pool; mid-mid: the pool paired with itself."""
    doc = table3()
    out = terms3([H, C, D2E, D2M, D3], ratio=0.9)
    (d2e,) = out["transfer"][RADICAL]
    assert d2e == {"template": names(doc, D2E)[0], "k_pair": pytest.approx(4.0 * weight(doc, D2E) * 0.9, rel=1e-12),
                   "phi": 1.0, "partner": MIDR, "product": LABEL, "partner_product": INTERNAL,
                   "valid_from_periods": 4}
    mids = by_template(out["transfer"][MIDR])
    assert mids[names(doc, D2M)[0]] == {"template": names(doc, D2M)[0],
                                        "k_pair": pytest.approx(2.0 * weight(doc, D2M) * 0.9, rel=1e-12),
                                        "phi": 1.0, "partner": RADICAL, "product": LABEL, "partner_product": ALKENE,
                                        "valid_from_periods": 4}
    assert mids[names(doc, D3)[0]] == {"template": names(doc, D3)[0],
                                       "k_pair": pytest.approx(4.0 * weight(doc, D3) * 0.9, rel=1e-12),
                                       "phi": 1.0, "partner": MIDR, "product": LABEL, "partner_product": INTERNAL,
                                       "valid_from_periods": 4}
    assert out["record"]["unmapped"] == [] and out["record"]["inert"] == []


def test_h_transfer_maps_with_the_donors_growing_sites():
    """End radical + chain: the chain's interior H, 2m - 8 sites at m periods
    (CH2 between two CH2), the chain picked by them; mid-chain radical + chain:
    the CH3 ends' 6 sites on every chain (constant: the I-221 form)."""
    doc = table3()
    out = terms3([H, C, T1, T2], ratio=0.5)
    (t1,) = out["transfer"][RADICAL]
    w = weight(doc, T1)
    assert {k: v for k, v in t1.items() if k not in ("k_pair", "sites")} == {
        "template": names(doc, T1)[0], "phi": 1.0, "partner": LABEL, "product": LABEL, "partner_product": MIDR,
        "periods_per_unit": 2, "valid_from_periods": 4}
    assert [t1["k_pair"] * c for c in t1["sites"]] == pytest.approx([-8.0 * w * 0.5, 0.0, 2.0 * w * 0.5, 0.0],
                                                                     rel=1e-12, abs=1e-30)
    assert t1["sites"][1] == 0.0 and t1["sites"][3] == 0.0
    (t2,) = out["transfer"][MIDR]
    assert t2 == {"template": names(doc, T2)[0], "k_pair": pytest.approx(6.0 * weight(doc, T2) * 0.5, rel=1e-12),
                  "phi": 1.0, "partner": LABEL, "product": LABEL, "partner_product": RADICAL,
                  "valid_from_periods": 4}
    check = out["record"]["checks"][t1["template"]]
    assert check["law_from_periods"] == check["reachable_from_periods"] == 4


def test_reverse_disproportionation_onto_interior_pools_times_phi():
    doc = table3()
    out = terms3([H, C, D1, DR, DRI], phi=0.25)
    alkene = by_template(out["transfer"][ALKENE])
    (end_mid,) = [t for t in alkene.values() if t["partner_product"] == MIDR]
    w = weight(doc, DR) if len(names(doc, DR)) == 1 else next(
        r["weight"] / r["divisor"] for r in doc["templates"] if r["template"] == end_mid["template"])
    assert end_mid["product"] == RADICAL and end_mid["partner"] == LABEL and end_mid["phi"] == 0.25
    assert [end_mid["k_pair"] * c for c in end_mid["sites"]] == pytest.approx([-8.0 * w, 0.0, 2.0 * w, 0.0],
                                                                               rel=1e-12, abs=1e-30)
    assert any(t["partner_product"] == RADICAL and "sites" not in t for t in alkene.values())
    (internal,) = out["transfer"][INTERNAL]
    wi = weight(doc, DRI)
    assert {k: internal[k] for k in ("phi", "partner", "product", "partner_product")} == {
        "phi": 0.25, "partner": LABEL, "product": MIDR, "partner_product": MIDR}
    assert [internal["k_pair"] * c for c in internal["sites"]] == pytest.approx(
        [-16.0 * wi, 0.0, 4.0 * wi, 0.0], rel=1e-12, abs=1e-30)          # two alkene carbons take the H
    assert out["record"]["unmapped"] == []


def test_backbiting_maps_as_a_move_and_its_inert_reverse_is_listed():
    doc = table3()
    out = terms3([H, C, BB, BBR])
    (move,) = out["move"][RADICAL]
    (name,) = names(doc, BB)
    w = next(r["weight"] for r in doc["templates"] if r["template"] == name)
    assert move == {"template": name, "k": pytest.approx(2.0 * w, rel=1e-12), "a": 0.0, "b": 1.0,
                    "periods_per_unit": 2, "product": MIDR, "valid_from_periods": 4}
    assert set(out["move"]) == {RADICAL}
    (inert,) = out["record"]["inert"]
    assert inert["kernel"] == BBR and inert["templates"] == names(doc, BBR)
    assert "no far-field site" in inert["reason"]
    assert out["record"]["unmapped"] == []


def test_a_law_exact_only_from_a_longer_length_is_listed_unless_the_far_law_is_asked():
    """The chain donor given 2 sites at 4 periods (its law says 0 = 2*4 - 8):
    the bilinear law holds only from 5 periods. By default listed (I-221's
    rule); with far_law mapped from 5 periods, the shorter lengths recorded."""
    doc = table3()
    (term,) = doc["kernels"][T1]["templates"]
    term["law"]["S1"]["exact"]["4"] = [[2, [["radical_interior:CH", 4]]]]
    strict = terms3([H, C, T1], doc=json.loads(json.dumps(doc)))
    assert RADICAL not in strict["transfer"]
    (listed,) = strict["record"]["unmapped"]
    assert listed["kernel"] == T1 and "bilinear" in listed["reason"] and "from 5" in listed["reason"]
    far = terms3([H, C, T1], doc=doc, far_law=True)
    (t1,) = far["transfer"][RADICAL]
    assert far["record"]["unmapped"] == [] and far["record"]["far_law"] is True
    assert far["record"]["checks"][t1["template"]]["law_from_periods"] == 5 == t1["valid_from_periods"]
    (short,) = far["record"]["short_lengths"]
    assert short["kernel"] == T1 and short["from_periods"] == 4 and short["law_from_periods"] == 5
    assert short["worst_relative"] > 0.0


def test_a_move_law_exact_only_from_a_longer_length_follows_the_same_rule():
    doc = table3()
    (term,) = doc["kernels"][BB]["templates"]
    term["law"]["exact"]["4"] = []
    strict = terms3([H, C, BB], doc=json.loads(json.dumps(doc)))
    assert strict["move"] == {}
    assert "affine" in strict["record"]["unmapped"][0]["reason"]
    far = terms3([H, C, BB], doc=doc, far_law=True)
    (move,) = far["move"][RADICAL]
    assert far["record"]["checks"][move["template"]]["law_from_periods"] == 5 == move["valid_from_periods"]


def test_an_intramolecular_move_that_keeps_the_pool_is_recorded_as_identity():
    """A mid-chain radical's own H shift along the chain keeps it a mid-chain
    radical: no MoM term, recorded."""
    doc = table3()
    (term,) = doc["kernels"][BBR]["templates"]
    law = json.loads(json.dumps(doc["kernels"][D3]["templates"][0]["law"]["S0"]))
    law["far"] = {"outcomes": [[3, [["radical_interior", 1, 0, 0]]]]}
    term["law"] = law
    term.pop("law_inert")
    out = terms3([H, C, BBR], doc=doc)
    assert out["move"] == {}
    (same,) = out["record"]["identity"]
    assert same["kernel"] == BBR and same["templates"] == names(doc, BBR)


def test_a_mid_chain_kernel_without_interior_pools_is_refused():
    with pytest.raises(mp.ProducerError, match="radical_interior"):
        mp.kernel_input(table3(), spec_of(BASE), [H, C, D2E], temperature=700.0, phi=1.0, label=LABEL,
                        radical_label=RADICAL, alkene_label=ALKENE, volume_ratio=1.0)


def test_kernel_mode_writes_the_interior_pools_and_their_kernels(setup):
    out = setup["tmp"] / "out"
    manifest = produce_kernels(setup, out, kernel_rates=str(TABLE3), kernels=ALL3, phi=0.5)
    pools = polymers(out)
    assert set(pools) == {LABEL, RADICAL, ALKENE, MIDR, INTERNAL}
    record = manifest["facts"]["kernels"]
    assert pools[MIDR]["interior_feature"] == "radical" and pools[INTERNAL]["interior_feature"] == "alkene"
    for pool in (LABEL, RADICAL, ALKENE, MIDR, INTERNAL):
        assert pools[pool].get("kernel_transfer") == record["transfer"].get(pool)
        assert pools[pool].get("kernel_move") == record["move"].get(pool)
    for pool in (MIDR, INTERNAL):
        assert pools[pool]["moments"] == [0.0, 0.0, 0.0] and pools[pool]["initial_mass"] == 0.0
        assert pools[pool]["end_groups"] == pools[LABEL]["end_groups"]
    assert manifest["facts"]["interior_radical_pool"] == MIDR
    assert manifest["facts"]["interior_alkene_pool"] == INTERNAL
    assert record["record"]["far_law"] is False


def test_kernel_mode_writes_no_interior_pool_without_a_kernel_into_it(setup):
    out = setup["tmp"] / "out"
    manifest = produce_kernels(setup, out, kernel_rates=str(TABLE2), kernels=[H, C, D1, DR])
    assert set(polymers(out)) == {LABEL, RADICAL, ALKENE}
    assert manifest["facts"]["interior_radical_pool"] is None and manifest["facts"]["interior_alkene_pool"] is None


def test_kernel_mode_passes_the_far_law_option_on(setup):
    out = setup["tmp"] / "out"
    manifest = produce_kernels(setup, out, kernel_rates=str(TABLE3), kernels=[H, C, BB], far_law=True)
    assert manifest["facts"]["kernels"]["record"]["far_law"] is True


def test_kernel_mode_converts_a_run_with_interior_pools(setup):
    out = setup["tmp"] / "out"
    produce_kernels(setup, out, kernel_rates=str(TABLE3), kernels=ALL3)
    mu = [float(x) for x in setup["spec"].macroscopic_moments_mol]
    d = 1e-3

    def chains(t):
        return mu[0] - d * t, mu[1] - 3 * d * t, mu[2] - 40 * d * t

    def grown(scale):
        return lambda t: (scale * d * t, scale * d * t, 20 * scale * d * t)
    write_run(out, times=[0.0, 10.0, 20.0, 40.0, 45.0],
              pools=[(LABEL, chains, 0.0), (RADICAL, grown(0.5), 0.0), (ALKENE, grown(0.5), 0.0),
                     (MIDR, grown(0.5), 0.0), (INTERNAL, grown(0.5), 0.0)], gas=ZERO_GAS, rows=[])
    manifest = produce_kernels(setup, out, mode="convert", kernel_rates=str(TABLE3), kernels=ALL3)
    assert manifest["status"] == "converted"
    assert manifest["census"]["mode"] == "kernels" and manifest["census"]["profile_changes"]


def test_a_constant_law_with_rounding_noise_keeps_the_constant_term():
    k = 3.096853570259764e-13
    F = {(n, m): k * (1.0 + 1e-15 * ((7 * n + 13 * m) % 5 - 2)) for n in range(4, 35) for m in range(4, 35)}
    n0, c, worst = mp._bilinear(F, 4, 34)
    assert n0 == 4 and c[1:] == [0.0, 0.0, 0.0] and c[0] == pytest.approx(k, rel=1e-14) and worst == 0.0
    r = {n: k * (1.0 + 1e-15 * (n % 3 - 1)) for n in range(4, 35)}
    n0, c, _w = mp._affine(r, 4, 34)
    assert n0 == 4 and c[1] == 0.0


def test_every_term_declares_its_validity_length_with_its_families():
    """The shortest chain each term's site law counts correctly, for the
    kernel reference to report the fraction of each family's chains below:
    homolysis 2*min_fragment - 1, recombination one past its zero-site
    profile, a pair transfer or move the length its law is exact from."""
    out = terms3([H, C, D2E, T1, BB])
    record = {v["template"]: v for v in out["record"]["validity"]}
    terms = [*out["homolysis"], *out["recombination"], *[t for ts in out["transfer"].values() for t in ts],
             *[t for ts in out["move"].values() for t in ts]]
    assert sorted(record) == sorted(t["template"] for t in terms)
    for t in terms:
        assert record[t["template"]]["valid_from_periods"] == t["valid_from_periods"]
    (h,) = out["homolysis"]
    assert h["valid_from_periods"] == 2 * h["min_fragment"] - 1 and record[h["template"]]["families"] == ["chain"]
    (c,) = out["recombination"]
    assert record[c["template"]]["families"] == ["radical_end"]
    assert c["valid_from_periods"] == max(out["record"]["checks"][c["template"]]["zero_site_periods"]) + 1
    families = {v["kernel"]: v["families"] for v in out["record"]["validity"]}
    assert families[D2E] == ["radical_end", "radical_interior"] and families[T1] == ["chain", "radical_end"]
    assert families[BB] == ["radical_end"]


def test_kernel_mode_exports_mu2_and_mw_with_each_pools_physical_end_mass(setup):
    """Homolysis: a chain (two H caps) becomes two end radicals (one H cap
    each) -- the end-group mass of the chains is conserved, so every chain's
    molar mass is its pool's a n + e_p and the whole mu2 and Mw are exact.
    MoM's polymer() needs mono-radical end groups, so the producer declares
    every pool with the deck's pair; the conversion must take each pool's
    end mass from its kernel-rates family (radical_end: one H), or it books
    an end-group deficit that is not there and withholds mu2 and Mw."""
    out = setup["tmp"] / "out"
    produce_kernels(setup, out)
    mu = [float(x) for x in setup["spec"].macroscopic_moments_mol]
    d = 1e-3

    def chains(t):
        return mu[0] - d * t, mu[1] - 3 * d * t, mu[2] - 40 * d * t

    def radicals(t):
        return 2 * d * t, 3 * d * t, 40 * d * t
    write_run(out, times=[0.0, 10.0, 20.0, 40.0, 45.0], pools=[(LABEL, chains, 0.0), (RADICAL, radicals, 0.0)],
              gas=ZERO_GAS, rows=[])
    manifest = produce_kernels(setup, out, mode="convert")
    offsets = {f: v["molar_mass_offset_g_mol"] for f, v in table()["families"].items()}
    a, e = manifest["facts"]["a_g_mol"], manifest["facts"]["e_g_mol"]
    assert offsets["chain"] == pytest.approx(e, rel=1e-12)
    s = (offsets["radical_end"] - e) / a                         # the radical pool's mass-DP offset
    for entry in mom_entries(out)[1:]:
        q, t = entry["quantities"], entry["time_s"]
        c, r = chains(t), radicals(t)
        mu2 = c[2] + r[2] + 2 * s * r[1] + s * s * r[0]
        assert q["whole_mu2_mol"] == pytest.approx(mu2, rel=1e-9)
        mu1 = c[1] + r[1] + s * r[0]
        assert q["whole_mu1_mol"] == pytest.approx(mu1, rel=1e-9)
        assert q["whole_Mw_g_mol"] == pytest.approx(
            (a * a * mu2 + 2 * a * e * mu1 + e * e * (c[0] + r[0])) / (a * mu1 + e * (c[0] + r[0])), rel=1e-9)
    assert manifest["facts"]["kernels"]["record"]["pool_end_mass_g_mol"] == {
        LABEL: offsets["chain"], RADICAL: offsets["radical_end"]}
