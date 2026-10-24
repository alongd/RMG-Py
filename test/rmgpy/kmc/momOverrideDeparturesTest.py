"""Red tests (I-256): a depository override in a kernel-rates term is a
named, quantified departure of MoM from kMC, never a refused kernel.

The PE toy catalog (``meltRunTest.toy_payloads``) carries one depository
override per kernel type at exact sizes, as RMG's training reactions do on
PE (n-decane, n-dodecane and their radicals): homolysis of n-dodecane into
pentyl + heptyl, their recombination, their disproportionation, and pentyl
abstracting an H from n-dodecane, each at half the template's weight.  kMC
fires those transitions at the override's weight; MoM's kernel carries the
template law at every size, so the producer maps the kernel and lists each
override with its sizes and its weight against the template's.
"""

import copy
import json
from pathlib import Path

import pytest

import melt_engine_support as S
from kernelReferenceTest import _rows_json
from meltRunTest import PE_LIVE, toy_payloads, write_catalog
from momProducerTest import BASE, LABEL, spec_of

from rmgpy.kmc import kernel_rates as K
from rmgpy.kmc import mom_producer as mp

T = 700.0
H = "homolysis:substrate=chain"
C = "recombination:radical=radical_end:CH2+radical=radical_end:CH2"
D = "disproportionation:acceptor=radical_end:CH2+donor=radical_end:CH2"
X = "h_transfer:acceptor=radical_end:CH2+donor=chain"
RATIO = 0.5
OVERRIDES = {
    "homolysis": (["CCCCCCCCCCCC"], ["[CH2]CCCC", "[CH2]CCCCCC"]),
    "recombination": (["[CH2]CCCC", "[CH2]CCCCCC"], ["CCCCCCCCCCCC"]),
    "disproportionation": (["[CH2]CCCC", "[CH2]CCCCCC"], ["CCCCC", "C=CCCCCC"]),
    "h_abstraction": (["[CH2]CCCC", "CCCCCCCCCCCC"], ["CCCCC", "CCC[CH]CCCCCCCC"]),
}


def entry(name):
    return {"depository": f"toy_{name}/training", "index": 1, "label": f"toy {name} at exact sizes"}


@pytest.fixture(scope="module")
def table(tmp_path_factory):
    payloads = toy_payloads([t for t in S.pe_templates() if t.name in PE_LIVE])
    for p in payloads:
        for t in p["templates"]:
            if t["name"] in OVERRIDES:
                r, q = OVERRIDES[t["name"]]
                w = RATIO * S.PE_WEIGHTS[t["name"]]
                t["overrides"] = [{"reaction": {"reactants": _rows_json(r), "products": _rows_json(q)},
                                   "weight": {"T": [300.0, 1000.0], "w": [w, w]}, "entry": entry(t["name"])}]
    path = write_catalog(tmp_path_factory.mktemp("ovr") / "pe_ovr.json", payloads, tag="pe-ovr",
                         bind=S.pe_mapping({1: 1}))
    return K.kernel_rates(json.loads(Path(path).read_text()), S.pe_spec({10: 1}), T)


def without_overrides(table):
    doc = copy.deepcopy(table)
    for kernel in doc["kernels"].values():
        for term in kernel["templates"]:
            term["overrides"] = []
    return doc


def terms(doc, kernels):
    return mp.kernel_input(doc, spec_of(BASE), list(kernels), temperature=T, phi=1.0, label=LABEL,
                           radical_label=LABEL + "_re", volume_ratio=1.0, alkene_label=LABEL + "_at",
                           interior_radical_label=LABEL + "_ri", far_law=True)


def overrides_of(table, kernels):
    return [(kid, term["template"], o) for kid in kernels for term in table["kernels"][kid]["templates"]
            for o in term["overrides"]]


def test_the_toy_carries_one_override_per_kernel(table):
    """The fixture's premise: each kernel's term carries its override, read at
    the override's own molecules, with the catalog's depository entry."""
    found = overrides_of(table, [H, C, D, X])
    assert sorted(kid for kid, _t, _o in found) == sorted([H, C, D, X])
    for _kid, _t, o in found:
        assert o["entry"]["depository"].endswith("/training")
    (_k, _t, h), = overrides_of(table, [H])
    assert h["reactants"] == [["chain", 12, None]]
    assert h["products"] == [["radical_end:CH2", 5], ["radical_end:CH2", 7]]


def test_overridden_homolysis_and_recombination_map_with_the_template_law(table):
    """Refused today ('an inert, unmodelled or overridden law has no moment
    form'): the kernels map, as on the same table without its overrides."""
    got, plain = terms(table, [H, C]), terms(without_overrides(table), [H, C])
    assert got["homolysis"] == plain["homolysis"] and got["recombination"] == plain["recombination"]
    assert got["record"]["laws"] == plain["record"]["laws"]


def test_pair_transfer_terms_with_overrides_map_with_the_template_law(table):
    """An overridden pair-transfer term is mapped, never left in unmapped."""
    got, plain = terms(table, [H, C, D, X]), terms(without_overrides(table), [H, C, D, X])
    assert {law["kernel"] for law in got["record"]["laws"]} == {H, C, D, X}
    assert got["transfer"] == plain["transfer"]
    assert got["record"]["unmapped"] == plain["record"]["unmapped"] == []


def test_every_override_of_a_mapped_term_is_a_named_departure(table):
    """One departure per override: its kernel, template, depository entry,
    exact sizes, multiplicity, the weight kMC fires it at and the template
    weight MoM's kernel carries there, and their ratio."""
    record = terms(table, [H, C, D, X])["record"]
    rows = {r["template"]: r for r in table["templates"]}
    want = overrides_of(table, [H, C, D, X])
    departures = record["depository_departures"]
    assert len(departures) == len(want) == 4
    by_name = {d["override"]: d for d in departures}
    for kid, name, o in want:
        d = by_name[o["name"]]
        assert d["kernel"] == kid and d["template"] == name and d["entry"] == o["entry"]
        assert d["reactants"] == o["reactants"] and d["products"] == o["products"]
        assert d["multiplicity"] == o["multiplicity"]
        assert d["w_kmc"] == o["weight"] and d["w_mom"] == float(rows[name]["weight"])
        assert d["kmc_over_mom"] == pytest.approx(RATIO, rel=1e-12)


def test_only_requested_kernels_list_departures(table):
    record = terms(table, [H, C])["record"]
    assert sorted(d["kernel"] for d in record["depository_departures"]) == [H, C]
