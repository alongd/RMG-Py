"""Red tests (I-278): the MoM producer reads the kernel-rates state model
``kmc-kernel-rates/4`` (chain types, oriented positions, codes by class).

A declaration whose backbone is not its own reverse (polypropylene: the
period CH2-CH(CH3)) has more than one chain type (pieces cut inside a period
carry one residue backbone atom), end classes per side of a cut and interior
classes per phase of the period, and its positions are measured from the
head or the tail.  The producer keeps one pool per (family, class, chain
type, oriented position class) and routes every product by the class and
type its code carries:

* homolysis is the uniform cut per cut position within the period: per pair
  of fragment classes and residue (periods the two fragments lose to their
  residues), fragments k and n - residue - k;
* position pools are oriented (``h<d>`` / ``t<d>``) on a type that is not
  its own reverse and folded (``d<d>``) on one that is;
* interior beta scission gives one term per residue;
* a class or type the table declares but no law reaches, a defect and a
  junction are listed by name;
* MoM declares a chain type's pools as end-group variants of the declared
  polymer, each end-featured pool on the site of its side and keyed by its
  class.

The table is built from the I-279 fixture catalog (12 payloads of the final
polypropylene catalog).  Every expectation is the table's own laws read per
length and position (``predict`` / ``predict_near``), never the producer's
moment algebra.
"""

import ast
import gzip
import json
import re
from collections import Counter
from pathlib import Path

import pytest

from rmgpy.kmc import kernel_rates as kr_
from rmgpy.kmc import mom_producer as mp
from rmgpy.kmc import run as kr

T = 700.0
FIXTURES = Path(__file__).parent / "fixtures"
DECK = FIXTURES / "decks" / "pp_deck.py"
LABEL = "PP"
LABELS = dict(label=LABEL, radical_label=LABEL + "_re", alkene_label=LABEL + "_ut",
              interior_radical_label=LABEL + "_rm", interior_alkene_label=LABEL + "_ui")
TOP = 19                                # the fixture table's exact_to_n


@pytest.fixture(scope="module")
def table():
    doc = json.loads(gzip.decompress((FIXTURES / "i279" / "pp_state_model_catalog.json.gz").read_bytes()))
    return kr_.kernel_rates(doc, kr.read_deck(DECK).spec, T)


@pytest.fixture(scope="module")
def out(table):
    return mp.kernel_input(table, kr.read_deck(DECK).spec, sorted(table["kernels"]), temperature=T, phi=1.0,
                           volume_ratio=1.0, far_law=True, positions=True, **LABELS)


def declared_classes(table):
    found = set()
    for t in table["chain_types"].values():
        for f in t["families"].values():
            found.update(f.get("classes") or ())
            for end in (f.get("ends") or {}).values():
                found.add(end["label_far"])
                found.update(end["labels"].values())
    return found


def pool_class(out):
    return {pool: cls for cls, pool in out["record"]["pools_by_class"].items()}


def feature_key(cls):
    return re.sub(r"[^A-Za-z0-9]", "", cls.split(":", 1)[1].split("/")[0])


# --------------------------------------------------------------------------
# the state model is read, or refused by name
# --------------------------------------------------------------------------
def test_a_v4_table_with_chain_types_is_read_in_positions_mode(table, out):
    assert table["chain_types"]
    assert sorted(out["record"]["chain_types"]) == sorted(table["chain_types"])


def test_outside_positions_mode_the_v4_state_model_is_refused_by_name(table):
    with pytest.raises(mp.ProducerError, match="positions mode"):
        mp.kernel_input(table, kr.read_deck(DECK).spec, sorted(table["kernels"]), temperature=T, phi=1.0,
                        volume_ratio=1.0, far_law=True, positions=False, **LABELS)


def test_every_declared_class_has_a_pool_of_its_own(table, out):
    by_class = out["record"]["pools_by_class"]
    assert set(by_class) == declared_classes(table)
    assert len(set(by_class.values())) == len(by_class)
    # the pools a mapped term reads or makes are a class each
    assert set(out["record"]["pools"]) <= set(by_class)
    for cls in ("radical_end:CH", "radical_end:CH2", "radical_end:CH/c2~c2", "radical_end:CH2/c1~c1",
                "chain", "chain/c1~c1", "chain/c2~c2"):
        assert cls in out["record"]["pools"], cls


# --------------------------------------------------------------------------
# homolysis: the uniform cut per cut position within the period
# --------------------------------------------------------------------------
def _cuts_of(term, classes):
    a, e = term["a"], term["min_fragment"]
    first = classes[term["daughter"]]
    second = classes[term.get("partner_daughter", term["daughter"])]
    e2, residue = term.get("partner_min_fragment", e), term.get("residue", 0)

    def at(n):
        got = Counter()
        for k in range(e, n - residue - e2 + 1):
            got[tuple(sorted([(first, k), (second, n - residue - k)]))] += a
        return got
    return at


def test_homolysis_is_the_uniform_cut_per_position_in_the_period(table, out):
    classes = pool_class(out)
    by_pool = out["homolysis_by_pool"]
    for kid, kernel in table["kernels"].items():
        if kernel["type"] != "homolysis":
            continue
        (cls,) = kernel["classes"]
        terms = by_pool[out["record"]["pools_by_class"][cls]]
        assert len(terms) >= 1
        valid = max(t["valid_from_periods"] for t in terms)
        for n in range(valid, TOP + 1):
            want = Counter()
            for term in kernel["templates"]:
                for m, prods in term["law"]["exact"].get(str(n), []):
                    want[tuple(sorted((p[0], int(p[1])) for p in prods))] += m
            got = Counter()
            for t in terms:
                got.update(_cuts_of(t, classes)(n))
            assert got == want, (kid, n)
    # the declared chain cuts at both bonds of the period: one pair of declared ends, one pair of the other types
    pairs = {tuple(sorted([classes[t["daughter"]], classes[t.get("partner_daughter", t["daughter"])]]))
             + (t.get("residue", 0),) for t in by_pool[LABEL]}
    assert pairs == {("radical_end:CH", "radical_end:CH2", 0), ("radical_end:CH/c2~c2", "radical_end:CH2/c1~c1", 1)}


def test_homolysis_weights_carry_every_template_site(table, out):
    rows = {r["template"]: r for r in table["templates"]}
    for kid, kernel in table["kernels"].items():
        if kernel["type"] != "homolysis":
            continue
        (cls,) = kernel["classes"]
        terms = out["homolysis_by_pool"][out["record"]["pools_by_class"][cls]]
        sites = [(float(rows[t["template"]]["weight"]), float(rows[t["template"]]["participants"][0][cls]["a"]))
                 for t in kernel["templates"]]
        assert sum(t["w"] * t["a"] for t in terms) == pytest.approx(sum(w * a for w, a in sites), rel=1e-12)
        assert sum(t["a"] for t in terms) == sum(a for _w, a in sites)


# --------------------------------------------------------------------------
# positions: oriented on an oriented type, folded on a type that is its own reverse
# --------------------------------------------------------------------------
def test_position_pools_are_oriented_only_on_an_oriented_type(table, out):
    pools = out["record"]["position_pools"]
    assert pools
    for pool, v in pools.items():
        symmetric = table["chain_types"][v["type"]]["symmetric"]
        assert re.fullmatch(r"bulk|d\d+" if symmetric else r"bulk|[ht]\d+", v["position"]), (pool, v)
        assert v["class"] in declared_classes(table)
    declared = [v for v in pools.values() if table["chain_types"][v["type"]]["declared"]]
    assert {v["position"][0] for v in declared if v["position"] != "bulk"} == {"h", "t"}


def test_the_end_backbite_lands_at_its_head_distance(table, out):
    kid = "intra_h_transfer:substrate=radical_end:CH2"
    (term,) = table["kernels"][kid]["templates"]
    row = {r["template"]: r for r in table["templates"]}[term["template"]]
    # the table: one radical_interior:C per site, at head distance 3 of the same chain
    got = kr_.predict_placed(term["law"], 12, 2, types=table["chain_types"])
    assert got == Counter({(("radical_interior:C", 12, "head_distance", 3),): 1})
    (move,) = out["move"][out["record"]["pools_by_class"]["radical_end:CH2"]]
    dst = out["record"]["position_pools"][move["product"]]
    assert (dst["class"], dst["position"]) == ("radical_interior:C", "h3")
    per = float(row["participants"][0]["radical_end:CH2"]["per_feature"])
    assert move["k"] * (move["a"] * 10 + move["b"]) == pytest.approx(float(row["weight"]) * per, rel=1e-12)


def test_transfers_route_products_by_the_class_and_type_of_their_codes(table, out):
    by_class = out["record"]["pools_by_class"]
    made = {t["product"] for t in out["transfer"].get(by_class["radical_end:CH/c2~c2"], [])}
    assert made == {by_class["chain/c2~c2"]}
    named = {p for ts in out["transfer"].values() for t in ts
             for p in (t["partner"], t["product"], t["partner_product"])}
    assert named <= set(by_class.values()) | set(out["record"]["position_pools"])


# --------------------------------------------------------------------------
# interior beta scission: one term per residue
# --------------------------------------------------------------------------
def test_interior_scission_gives_one_term_per_residue(table, out):
    classes = pool_class(out)
    positions = out["record"]["position_pools"]
    (bulk,) = [p for p, v in positions.items() if v["class"] == "radical_interior:C" and v["position"] == "bulk"]
    terms = out["interior_scission"][bulk]
    rows = {r["template"]: r for r in table["templates"]}
    (term,) = table["kernels"]["beta_scission:substrate=radical_interior:C"]["templates"]
    w = float(rows[term["template"]]["weight"])
    # the table, far from both ends: two outcomes per feature, one at each beta bond
    by_residue = Counter()
    for outcome, m in kr_.predict_near(term["law"], None, 15, 2, 15, table["chain_types"]).items():
        by_residue[15 - sum(n for _c, n in outcome)] += m
    assert by_residue == Counter({0: 1, 1: 1})
    assert sorted(t.get("residue", 0) for t in terms) == [0, 1]
    for t in terms:
        assert t["k"] == pytest.approx(w * by_residue[t.get("residue", 0)], rel=1e-12)
        assert sum(p["moments"][1][0] for p in t["products"]) == pytest.approx(-t.get("residue", 0), abs=1e-9)
        assert sum(p["moments"][1][1] for p in t["products"]) == pytest.approx(1.0, abs=1e-12)
    made = {classes.get(p["pool"]) or positions[p["pool"]]["class"] for t in terms for p in t["products"]}
    assert made == {"alkene_terminal:CH2=C", "radical_end:CH", "alkene_terminal:CH2=C/c1~c1", "radical_end:CH/c2~c2"}


# --------------------------------------------------------------------------
# listed by name, never dropped
# --------------------------------------------------------------------------
def test_unreached_defects_and_junctions_are_listed_by_name(table, out):
    record = out["record"]
    assert record["unreached"] == table["unreached"]
    defects = {e["kernel"] for e in record["inert"] if "defect" in e["reason"]}
    assert defects >= {"recombination:radical=radical_end:CH+radical=radical_end:CH2/c1~c1",
                       "recombination:radical=radical_end:CH/c2~c2+radical=radical_end:CH2/c1~c1"}
    assert "addition:alkene=alkene_terminal:CH2=C+radical=radical_end:CH" in {
        e["kernel"] for e in record["position_refused"]}


# --------------------------------------------------------------------------
# MoM: chain types as end-group variants, end pools on their side, keyed by class
# --------------------------------------------------------------------------
def _polymers(text):
    found = {}
    for call in ast.parse(text).body:
        node = getattr(call, "value", None)
        if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "polymer":
            kw = {k.arg: ast.literal_eval(k.value) for k in node.keywords}
            found[kw["label"]] = kw
    return found


def test_mom_declares_chain_types_as_end_group_variants(table, tmp_path):
    from rmgpy.molecule import Molecule
    deck = kr.read_deck(DECK)
    source = {"table": table, "ids": sorted(table["kernels"]), "phi": 1.0, "far_law": True, "positions": True,
              "source": {"table": str(tmp_path / "pp.kernel_rates.json"), "table_sha256": "0" * 64}}
    text, facts = mp.build_input(deck, mp.deck_context(DECK), temperature=T, end=10.0, families=[], kernels=source)
    polymers = _polymers(text)
    record = facts["kernels"]["record"]
    by_class = record["pools_by_class"]
    label = facts["label"]

    def smiles(groups):
        return [Molecule(smiles=g).to_smiles() for g in groups]
    head, tail = smiles(polymers[label]["end_groups"])
    assert smiles(polymers[by_class["chain/c1~c1"]]["end_groups"]) == [head, "[CH3]"]
    assert smiles(polymers[by_class["chain/c2~c2"]]["end_groups"]) == ["C[CH2]", tail]
    for cls, site in (("radical_end:CH", "primary"), ("radical_end:CH2", "secondary"),
                      ("radical_end:CH2/c1~c1", "secondary"), ("radical_end:CH/c2~c2", "primary")):
        pool = polymers[by_class[cls]]
        assert pool["end_radical_site"] == site, cls
        assert pool["end_class"] == feature_key(cls), cls
    for pool, v in record["position_pools"].items():
        if pool in polymers:
            assert polymers[pool]["end_class"] == feature_key(v["class"]), pool
    # the conversion reads each typed pool's end mass from its type
    c1 = table["chain_types"]["c1~c1"]["families"]["chain"]["molar_mass_offset_g_mol"]
    assert record["pool_end_mass_g_mol"][by_class["chain/c1~c1"]] == pytest.approx(c1, rel=1e-12)
