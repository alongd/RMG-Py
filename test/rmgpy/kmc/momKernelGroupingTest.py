"""Red tests (I-258): the MoM producer groups a homolysis or recombination
kernel's templates before it checks MoM's form, as it already does for pair
transfer and moves.

The fixture holds the homolysis and end-recombination kernels of the PE toy
catalog's kernel-rates table at 700 K (kmc-i201 f0a02bdb3,
``fixtures/decks/pe_toy_deck.py``).  RMG splits each over context templates.
Homolysis has an ethyl cut (2 sites), a propyl cut (2), the bulk cut over
[4, n - 4] (n - 7) and the n-hexane cut into two propyls (at n = 6 only).
Their union is the uniform cut over [2, n - 2], n - 3 sites, from n = 6.  No
template cuts n-butane or n-pentane, so those two lengths are outside the
law.  Recombination has propyl + propyl, long + propyl, long + long and long +
ethyl, each a fixed-size profile.  No template joins ethyl + ethyl or ethyl +
propyl.

MoM takes the far-field site-weighted mean weight of a group.  Each template
whose weight departs from that mean is listed in ``position_departures``
with kmc_over_mom.
"""

import copy
import json
from pathlib import Path

import pytest

from momProducerTest import LABEL

from rmgpy.kmc import mom_producer as mp
from rmgpy.kmc import run as kr

FIXTURES = Path(__file__).parent / "fixtures"
TABLE = FIXTURES / "i258_pe_toy_split_kernels_700K.json"
DECK = FIXTURES / "decks" / "pe_toy_deck.py"
H = "homolysis:substrate=chain"
C = "recombination:radical=radical_end:CH2+radical=radical_end:CH2"
RADICAL = LABEL + "_re"
ETHYL, PROPYL, BULK, HEXANE = ("1e40bb1bbf3d:69cd93b3f0a7", "25429a1d318d:2e21ab306d36",
                                "7a78b330d9a0:f5338ad0db2e", "25429a1d318d:e90a01b09296")
LONG_PROPYL, LONG_LONG = "8b190ce1aaeb:1b3b7118bc42", "9a0c689d52bf:e977ec702e62"


def table():
    return json.loads(TABLE.read_text())


def weight(doc, name):
    return next(float(r["weight"]) for r in doc["templates"] if r["template"] == name)


def terms(doc, kernels):
    return mp.kernel_input(doc, kr.read_deck(DECK).spec, list(kernels), temperature=700.0, phi=1.0, label=LABEL,
                           radical_label=RADICAL, volume_ratio=1.0)


def test_split_homolysis_maps_as_one_uniform_cut():
    """The union is the uniform cut from n = 6, the law's validity length
    (n-butane and n-pentane are cut by no template); the deck's chains are 60
    periods."""
    doc = table()
    out = terms(doc, [H])
    joined = "+".join(sorted([ETHYL, PROPYL, BULK, HEXANE]))
    assert out["homolysis"] == [{"template": joined, "w": weight(doc, BULK), "phi": 1.0, "a": 1.0, "b": -3.0,
                                 "min_fragment": 2, "periods_per_unit": 2, "daughter": RADICAL,
                                 "valid_from_periods": 6}]
    assert out["record"]["inert"] == []
    assert out["record"]["checks"][joined]["uniform_from_periods"] == 6


def test_a_template_off_the_mean_weight_is_a_named_position_departure():
    """The ethyl cut's weight is 2.05% above the bulk's; the propyl cut's
    differs in the 14th digit only, which is no departure."""
    doc = table()
    (d,) = terms(doc, [H])["record"]["position_departures"]
    assert d["kernel"] == H and d["template"] == ETHYL
    assert d["w_kmc"] == weight(doc, ETHYL) and d["w_mom"] == weight(doc, BULK)
    assert d["kmc_over_mom"] == pytest.approx(1.0205, abs=1e-4)
    assert d["sites"] == {"a": 0.0, "b": 2.0}


def test_an_override_departs_from_the_grouped_weight():
    doc = table()
    record = terms(doc, [H])["record"]
    ethyl = [d for d in record["depository_departures"] if d["template"] == ETHYL]
    assert ethyl and all(d["w_mom"] == weight(doc, BULK) for d in ethyl)


def test_a_position_cut_twice_is_refused_by_name():
    doc = table()
    term = copy.deepcopy(next(t for t in doc["kernels"][H]["templates"] if t["template"] == ETHYL))
    row = copy.deepcopy(next(r for r in doc["templates"] if r["template"] == ETHYL))
    term["template"] = row["template"] = "copy:" + ETHYL
    doc["kernels"][H]["templates"].append(term)
    doc["templates"].append(row)
    with pytest.raises(mp.ProducerError, match=r"not the uniform cut.*\(2, 6\) covered 4 times, 2 expected"):
        terms(doc, [H])


def test_recombination_alone_maps_from_where_every_pair_is_joined():
    doc = table()
    (out,) = terms(doc, [C])["recombination"]
    assert out["template"] == "+".join(sorted(t["template"] for t in doc["kernels"][C]["templates"]))
    assert out["k_pair"] == weight(doc, LONG_LONG) / 2 and out["sites_per_chain"] == 1.0
    assert out["valid_from_periods"] == 3


def test_recombination_pairs_no_template_joins_are_refused_by_name():
    """Homolysis makes ethyl radicals (min_fragment 2); no template joins
    ethyl + ethyl or ethyl + propyl."""
    with pytest.raises(mp.ProducerError, match=r"profile.*\(2, 2\): 0.*\(2, 3\): 0"):
        terms(table(), [H, C])


def test_a_recombination_template_off_the_mean_weight_is_a_named_departure():
    doc = table()
    for r in doc["templates"]:
        if r["template"] == LONG_PROPYL:
            r["weight"] = 1.1 * float(r["weight"])
    out = terms(doc, [C])
    assert out["recombination"][0]["k_pair"] == weight(doc, LONG_LONG) / 2
    (d,) = out["record"]["position_departures"]
    assert d["template"] == LONG_PROPYL and d["kmc_over_mom"] == pytest.approx(1.1, rel=1e-12)
