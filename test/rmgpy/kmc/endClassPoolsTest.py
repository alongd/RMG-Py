"""Red tests (I-266): the PE toy's C-C kernels on a period grid, carried by
MoM through per-end-class pools.

The table is ``kernel_rates`` on the PE toy catalog's R_Recombination and
R_Addition_MultipleBond payloads (``fixtures/i261_pe_toy_recombination_payloads.json``,
``fixtures/i262_pe_toy_addition_payloads.json``; ``fixtures/decks/pe_toy_deck.py``)
at 700 K.  Lengths count CH2 periods, two per declared unit.

Before I-266 the table counted pristine chains at whole declared units only,
so the n-propane and n-pentane cuts kMC makes were booked
``small_species_participant`` and the homolysis law was empty at 3 and 5
periods; a class that exists at one length only (the methyl radical at 1
period, ethylene at 2) had no law ("no probe ... at the law's reach"); and the
producer refused, in order: beta scission and addition (no MoM kernel type),
recombination of two end-radical classes, homolysis making a methyl, and
2-period radicals joining below the homolysis law's 6-period validity.

Here: a chain state exists at every period; a short class gets its law at
its lengths; end-radical beta scission releasing the declared monomer maps
onto ``k_depropagation`` of the radical pool, monomer addition onto
``kernel_propagation``, an addition joining two pools onto a
``kernel_recombination`` term with a ``partner``, the methyl class onto a
pool of its own; the templates that act at single sizes only (the methyl cuts
of ethane, propane and n-decane, methyl + ethyl and methyl + nonyl) are listed
as size departures; interior beta scission is refused by name (its fragments'
lengths depend on where the radical sits, which no pool carries).
"""

import json
from collections import Counter
from pathlib import Path

import pytest

from momProducerTest import LABEL

from rmgpy.kmc import kernel_rates as K
from rmgpy.kmc import mom_producer as mp
from rmgpy.kmc import run as kr

FIXTURES = Path(__file__).parent / "fixtures"
DECK = FIXTURES / "decks" / "pe_toy_deck.py"
RADICAL, ALKENE, MID, INTERNAL = (LABEL + s for s in ("_re", "_ut", "_rm", "_ui"))
METHYL_POOL = LABEL + "_re1"

H = "homolysis:substrate=chain"
C = "recombination:radical=radical_end:CH2+radical=radical_end:CH2"
C_METHYL = "recombination:radical=radical_end:CH2+radical=radical_end:CH3"
C_METHYLS = "recombination:radical=radical_end:CH3+radical=radical_end:CH3"
BETA_END = "beta_scission:substrate=radical_end:CH2"
BETA_MID = "beta_scission:substrate=radical_interior:CH"
ADD_MONOMER = "addition:alkene=alkene_terminal:CH2=CH2+radical=radical_end:CH2"
ADD_VINYL = "addition:alkene=alkene_terminal:CH2=CH+radical=radical_end:CH2"

PROPANE, PENTANE = "50e2a74af438:1847838c1a3d", "4824c5f6a98d:9d1848122f70"
METHYL_CUTS = ("72b4c5fe6565:2cec868975dd", PROPANE, "64c3fa6c2701:273afe5ef301")   # ethane, propane, n-decane
METHYL_JOINS = ("e3492f6c27d6:1b7c51ac8b52", "925f4fce5ca2:3c7862ef7633")            # + ethyl, + nonyl


def spec():
    return kr.read_deck(DECK).spec


@pytest.fixture(scope="module")
def table():
    s = spec()
    payloads = []
    for name in ("i261_pe_toy_recombination_payloads.json", "i262_pe_toy_addition_payloads.json"):
        payloads += json.loads((FIXTURES / name).read_text())["payloads"]
    doc = {"schema": "kmc-template-catalog/1", "payloads": payloads,
           "provenance": {"chemistry_sha256": kr.chemistry_sha256(s)}}
    return K.kernel_rates(doc, s, 700.0)


def terms_of(table, kid):
    return {t["template"]: t for t in table["kernels"][kid]["templates"]}


def row(table, name):
    return next(r for r in table["templates"] if r["template"] == name)


def mapped(table, kernels, **kw):
    return mp.kernel_input(table, spec(), list(kernels), temperature=700.0, phi=1.0, label=LABEL,
                           radical_label=RADICAL, alkene_label=ALKENE, interior_radical_label=MID,
                           interior_alkene_label=INTERNAL, volume_ratio=kw.pop("ratio", 1.0), **kw)


def carried(table):
    """Every kernel of the table but interior beta scission."""
    return sorted(set(table["kernels"]) - {BETA_MID})


# --------------------------------------------------------------------------
# kernel_rates: the period grid
# --------------------------------------------------------------------------
def test_a_chain_state_exists_at_every_period(table):
    n_hi = table["laws_exact_to_n"]
    assert table["chain_lengths_probed"] == list(range(1, n_hi + 1))
    assert 5 in table["chain_lengths_probed"]


def test_the_odd_period_cuts_join_the_homolysis_kernel(table):
    """n-propane and n-pentane are chain states now: their cuts are homolysis
    templates, not small-species entries."""
    only = {e["template"] for e in table["kmc_only"]}
    assert PROPANE not in only and PENTANE not in only
    hom = terms_of(table, H)
    assert PROPANE in hom and PENTANE in hom
    union = {n: Counter() for n in (3, 4, 5)}
    for t in hom.values():
        for n in union:
            for m, prods in t["law"]["exact"].get(str(n), []):
                union[n][tuple(sorted(p[1] for p in prods))] += m
    assert union == {3: Counter({(1, 2): 1}), 4: Counter({(2, 2): 1}), 5: Counter({(2, 3): 2})}


def test_a_class_at_its_short_lengths_has_its_law_there(table):
    """The methyl radical (1 period) and ethylene (2 periods) exist at one
    length: the law is read there, a join of the whole reactants."""
    (term,) = table["kernels"][C_METHYLS]["templates"]
    for slot in ("S0", "S1"):
        law = term["law"][slot]
        assert law["class"] == "radical_end:CH3" and law["exact_to_n"] == 1
        assert law["exact"] == {"1": [[1, [["join", "chain", 0]]]]}
    (term,) = table["kernels"][ADD_MONOMER]["templates"]
    law = term["law"]["S0"]
    assert law["class"] == "alkene_terminal:CH2=CH2"
    assert law["exact"] == {"2": [[2, [["join", "radical_end", 0]]]]}
    refused = [t["template"] for k in table["kernels"].values() for t in k["templates"]
               if "law's reach" in (t.get("law_refused") or "")]
    assert refused == []


# --------------------------------------------------------------------------
# the producer: the PE toy's kernels on MoM pools
# --------------------------------------------------------------------------
def test_every_kernel_but_interior_beta_scission_maps(table):
    out = mapped(table, carried(table))
    rec = out["record"]
    assert rec["unmapped"] == []
    assert {k for e in rec["inert"] for k in [e["kernel"]]} <= {C_METHYL}
    assert rec["pools"] == {"chain": LABEL, "radical_end:CH2": RADICAL, "radical_end:CH3": METHYL_POOL,
                            "alkene_terminal:CH2=CH": ALKENE, "radical_interior:CH": MID}


def test_interior_beta_scission_is_refused_by_name(table):
    with pytest.raises(mp.ProducerError, match=r"beta_scission:substrate=radical_interior:CH.*position"):
        mapped(table, [H, C, BETA_MID])


def test_the_methyl_cuts_are_size_departures_and_homolysis_is_the_uniform_cut_from_ethyl(table):
    out = mapped(table, carried(table))
    (h,) = out["homolysis"]
    assert h["min_fragment"] == 2 and h["b"] == pytest.approx(-3.0 * h["a"], rel=1e-12)
    assert h["valid_from_periods"] == 3 and h["daughter"] == RADICAL
    assert set(h["template"].split("+")).isdisjoint(METHYL_CUTS)
    sizes = {(d["kernel"], d["template"]) for d in out["record"]["size_departures"]}
    assert {(H, t) for t in METHYL_CUTS} | {(C_METHYL, t) for t in METHYL_JOINS} <= sizes


def test_two_period_radicals_join_inside_the_homolysis_law(table):
    """(4): radicals from 2 periods join into chains of 4, inside the
    homolysis law's validity (3 periods)."""
    out = mapped(table, [H, C])
    (c,) = out["recombination"]
    assert c["product"] == LABEL and c["valid_from_periods"] <= 2


def test_end_beta_scission_is_the_radical_pools_depropagation(table):
    out = mapped(table, [H, C, BETA_END])
    ((name, t),) = terms_of(table, BETA_END).items()
    assert out["depropagation"] == {RADICAL: {"template": name, "k": t["k_per_feature"], "valid_from_periods": 6}}
    assert out["record"]["monomer"] == {"smiles": "C=C", "class": "alkene_terminal:CH2=CH2", "periods": 2}
    deps = [d for d in out["record"]["depository_departures"] if d["kernel"] == BETA_END]
    assert sorted(d["reactants"][0][1] for d in deps) == [6, 7, 8, 9, 10, 12]


def test_monomer_addition_is_the_radical_pools_propagation(table):
    out = mapped(table, [H, C, BETA_END, ADD_MONOMER], ratio=0.5)
    ((name, t),) = terms_of(table, ADD_MONOMER).items()
    (p,) = out["propagation"][RADICAL]
    # sites: 2 on ethylene (either end), 1 on a radical of 4 periods or more
    assert p == {"template": name, "k_pair": t["w_over_divisor"] * 2 * 1 * 0.5, "valid_from_periods": 4}


def test_vinyl_addition_joins_the_radical_and_alkene_pools_into_the_interior_radical_pool(table):
    out = mapped(table, carried(table), ratio=0.5)
    ((name, t),) = terms_of(table, ADD_VINYL).items()
    (j,) = [x for x in out["joins"][RADICAL] if x["template"] == name]
    assert j == {"template": name, "k_pair": t["w_over_divisor"] * 0.5, "sites_per_chain": 1.0,
                 "partner": ALKENE, "partner_sites_per_chain": 1.0, "periods_per_unit": 2, "product": MID,
                 "valid_from_periods": 5}


def test_methyl_recombination_maps_on_the_methyl_pool_which_nothing_feeds(table):
    out = mapped(table, carried(table))
    ((name, t),) = terms_of(table, C_METHYLS).items()
    assert out["joins"][METHYL_POOL] == [{"template": name, "k_pair": t["w_over_divisor"], "sites_per_chain": 1.0,
                                          "periods_per_unit": 2, "product": LABEL, "valid_from_periods": 1}]
    assert out["record"]["unfed_pools"] == [METHYL_POOL]
