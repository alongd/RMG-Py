"""Red tests (I-273): the MoM producer's positions mode
(``mom_producer.kernel_input(..., positions=True)``).

A kernel-rates table states where every mid-chain feature sits (a fixed end
distance d, or the bulk of a chain) and how each interior participant reacts
at every end distance near an end (``law_end_distance`` of the terms with
sites there only, ``law_near_end`` of the law terms).  In positions mode the
producer keeps mid-chain radicals and internal alkenes in position-class
pools: one per end distance d < D and one bulk pool (d >= D, the feature
uniform over those positions), D one past every end distance where a site
count differs from the far field or a product is placed.  Every term then has
one per-chain rate per pool, so:

* the terms the table marks ``law_inert`` (an interior participant with
  sites near an end only: the reverse halves of reciprocal pairs) map, at
  the rate their end-distance law states;
* every transfer and move takes chains from the pool of the class they read
  into the pool of the class they make;
* interior beta scission maps per class pool (``interior_scission``), its
  fragments' length moments fixed coefficients of the pool's moments;
* a product whose position depends on two reactants (a join's junction) has
  no pool: refused by name;
* a bulk-pool product placed relative to the reactant's feature that lands
  below D (the bulk's boundary) is carried in the bulk pool and listed.

The catalog is generic (``melt_engine_support``'s toy linear polymer, the
hand-written templates of the position and near-end tests).  Every
expectation is an explicit sum over chain lengths and every position of the
feature in the class, the table's own laws read per position
(``predict_placed`` / ``predict_near``), never the producer's moment algebra.
"""

import ast
import json
import math
from collections import Counter, defaultdict
from fractions import Fraction

import pytest

import melt_engine_support as S
from meltRunTest import toy_payloads, write_catalog
import kernelRatesNearEndTest as N
import kernelRatesPositionTest as P
from momProducerTest import BASE, CHECKOUT, LABEL, make_setup, spec_of

from rmgpy.kmc import kernel_rates as kr_
from rmgpy.kmc import mom_producer as mp
from rmgpy.polymer_spec import LinearHomopolymerSpec

RADICAL, ALKENE, MID, INTERNAL = (LABEL + s for s in ("_re", "_ut", "_rm", "_ui"))
RATIO = 0.75
INTERIOR = ("radical_interior", "alkene_internal")
BASE_POOL = {"chain": LABEL, "radical_end": RADICAL, "alkene_terminal": ALKENE, "radical_interior": MID,
             "alkene_internal": INTERNAL}
BETA = "beta_scission:substrate=radical_interior:CH"


def templates():
    return [t for t in S.pe_templates() if t.name in ("h_abstraction", "beta_scission")] + \
        [P.end_acceptor(), P.end_donor(), P.end_shift(), P.backbite(), N.mid_donor(), N.alkene_acceptor()]


@pytest.fixture(scope="module")
def table(tmp_path_factory):
    path = write_catalog(tmp_path_factory.mktemp("positions") / "catalog.json", toy_payloads(templates()),
                         tag="positions-toy", bind=S.pe_mapping({1: 1}))
    return kr_.kernel_rates(json.loads(path.read_text()), LinearHomopolymerSpec.from_mapping(S.pe_mapping({10: 1})),
                            700.0)


@pytest.fixture(scope="module")
def out(table):
    return mp.kernel_input(table, spec_of(BASE), sorted(table["kernels"]), temperature=700.0, phi=1.0, label=LABEL,
                           radical_label=RADICAL, alkene_label=ALKENE, interior_radical_label=MID,
                           interior_alkene_label=INTERNAL, volume_ratio=RATIO, far_law=True, positions=True)


# --------------------------------------------------------------------------
# brute force: the table's laws read at every position of every class
# --------------------------------------------------------------------------
def boundary(table):
    """1 + every end distance where an interior participant's sites differ from the far field (the keys of the
    end-distance and near-end laws) or where a law places a product at a fixed end distance."""
    found = {0}

    def scan(x):
        if isinstance(x, dict):
            if x.get("anchor") == "end_distance":
                found.update(int(d) for d in x["at"])
            if "end_distance" in x and isinstance(x["end_distance"], int) and len(x) == 1:
                found.add(x["end_distance"])
            for k, v in x.items():
                if k not in ("exact",):
                    scan(v)
        elif isinstance(x, list):
            for v in x:
                scan(v)
    for k in table["kernels"].values():
        for t in k["templates"]:
            for key in ("law", "law_end_distance", "law_near_end"):
                scan(t.get(key))
    starts = set()

    def bulk(x):
        if isinstance(x, dict):
            if set(x) == {"bulk"} and all(isinstance(v, int) for v in x["bulk"]):
                starts.add(x["bulk"][0])
            for v in x.values():
                bulk(v)
        elif isinstance(x, list):
            for v in x:
                bulk(v)
    for k in table["kernels"].values():
        for t in k["templates"]:
            bulk(t.get("law"))
            bulk(t.get("law_end_distance"))
    return max([1 + d for d in found] + list(starts))


def pool_of(family, d, D):
    base = BASE_POOL[family]
    if family not in INTERIOR:
        return base
    return base if d >= D else f"{base}_d{d}"


def canonical(family, L):
    """{canonical position: how many backbone positions of a chain of L atoms it stands for}."""
    w = kr_.WIDTH[family]
    out = Counter()
    for x in range(1, L - w):
        out[min(x, L - w - x)] += 1
    return out


def slot_law(term, slot, nslots):
    """(law, near-end law) of one slot of a kernel term: the law (or its end-distance law when inert)."""
    if term.get("law") is not None:
        law = term["law"] if nslots == 1 else term["law"][slot]
        near = term.get("law_near_end")
        near = near if nslots == 1 else (near or {}).get(slot)
        return law, near
    led = term["law_end_distance"]
    return (led if nslots == 1 else led[slot]), None


def states(family, n, nbb, D):
    """{pool: [(x or None, weight)]}: the states of a reactant of n periods by pool, each class's positions
    weighted uniformly (the bulk's both halves)."""
    if family not in INTERIOR:
        return {BASE_POOL[family]: [(None, 1.0)]}
    out = defaultdict(list)
    for x, mult in canonical(family, n * nbb).items():
        out[pool_of(family, x, D)].append((x, mult))
    return {p: [(x, m / sum(mm for _x, mm in xs)) for x, m in xs] for p, xs in out.items()}


def outcomes(law, near, n, nbb, x):
    if law["anchor"] == "end_distance":
        return kr_.predict_placed(law, n, nbb, x)
    if law["anchor"] == "interior":
        return kr_.predict_near(law, near, n, nbb, x)
    if law.get("positioned"):
        return kr_.predict_placed(law, n, nbb)
    if n <= int(law["exact_to_n"]):
        return Counter({tuple(sorted((kr_.family_of(p[0]), int(p[1])) for p in prods)): m
                        for m, prods in law["exact"].get(str(n), [])})
    return kr_.predict(law, n, nbb)


def product_pool(product, D, src_pool, src_family, middle):
    """The pool of a product (family, n[, d]); a bulk-pool reactant's product placed relative to its feature
    stays in the bulk pool (``middle``: the products at the bulk's middle position, a fixed one is there too)."""
    family = product[0]
    if family not in INTERIOR:
        return BASE_POOL[family]
    d = product[2]
    if src_family in INTERIOR and src_pool == BASE_POOL[src_family] and d < D and product not in middle:
        return BASE_POOL[family]
    return pool_of(family, d, D)


def flows(law, near, cls, n, nbb, D):
    """{(source pool, product pool): sites per chain of the source pool} of one slot law on n periods."""
    family = kr_.family_of(cls)
    out = Counter()
    for pool, xs in states(family, n, nbb, D).items():
        middle = set()
        if family in INTERIOR and pool == BASE_POOL[family]:
            mid = max(canonical(family, n * nbb))
            middle = {p for o in outcomes(law, near, n, nbb, mid) for p in o}
        for x, w in xs:
            for outcome, m in outcomes(law, near, n, nbb, x).items():
                (product,) = outcome
                out[(pool, product_pool(product, D, pool, family, middle))] += w * m
    return out


def mom_pair(terms, n, m):
    """{(owner, partner, product, partner_product): rate per chain pair} of the transfer terms (m3/(mol s),
    kMC's convention: k_pair over the volume ratio)."""
    out = Counter()
    for owner, ts in terms["transfer"].items():
        for t in ts:
            k = t["k_pair"] / RATIO
            c = t.get("sites", [1.0, 0.0, 0.0, 0.0])
            out[(owner, t["partner"], t["product"], t["partner_product"])] += \
                k * (c[0] + c[1] * n + c[2] * m + c[3] * n * m)
    return out


def valid_from(terms):
    return max([t["valid_from_periods"] for ts in terms["transfer"].values() for t in ts] +
               [t["valid_from_periods"] for ts in terms["move"].values() for t in ts] +
               [t["valid_from_periods"] for ts in terms["interior_scission"].values() for t in ts])


# --------------------------------------------------------------------------
# the pools
# --------------------------------------------------------------------------
def test_mid_chain_features_get_a_pool_per_end_distance_below_the_boundary_and_a_bulk_pool(table, out):
    record = out["record"]
    D = boundary(table)
    assert record["position_boundary"] == D and D > 2
    pools = record["position_pools"]
    for family in INTERIOR:
        mine = {p: v for p, v in pools.items() if v["family"] == family}
        assert mine[BASE_POOL[family]]["position"] == "bulk"
        assert {v["position"] for v in mine.values()} <= {"bulk"} | {f"d{d}" for d in range(1, D)}
    assert {v["position"] for v in pools.values() if v["family"] == "radical_interior"} == \
        {"bulk"} | {f"d{d}" for d in range(1, D)}


def test_every_pool_a_term_reads_or_makes_is_a_position_pool_or_a_family_pool(out):
    known = set(BASE_POOL.values()) | set(out["record"]["position_pools"])
    for owner, ts in out["transfer"].items():
        for t in ts:
            assert {owner, t["partner"], t["product"], t["partner_product"]} <= known
    for owner, ts in out["move"].items():
        for t in ts:
            assert {owner, t["product"]} <= known
    for owner, ts in out["interior_scission"].items():
        for t in ts:
            assert {owner} | {p["pool"] for p in t["products"]} <= known


# --------------------------------------------------------------------------
# transfers and moves against the table on explicit chains
# --------------------------------------------------------------------------
def brute_pairs(table, n, m, D):
    """{(owner, partner, product, partner_product): rate per chain pair} of every two-participant term of the
    table, summed over templates, each reactant's class read at every position."""
    nbb = table["moment_unit"]["backbone_atoms_per_period"]
    rows = {r["template"]: r for r in table["templates"]}
    out = Counter()
    for kid, k in table["kernels"].items():
        if len(k["classes"]) != 2 or k["type"] in ("recombination", "addition"):
            continue
        for t in k["templates"]:
            if t.get("law") is None and not t.get("law_end_distance"):
                continue
            r = rows[t["template"]]
            w = r["weight"] / r["divisor"]
            f0 = flows(*slot_law(t, "S0", 2), k["classes"][0], n, nbb, D)
            f1 = flows(*slot_law(t, "S1", 2), k["classes"][1], m, nbb, D)
            for (a, a2), s0 in f0.items():
                for (b, b2), s1 in f1.items():
                    out[(a, b, a2, b2)] += w * s0 * s1
    return out


def test_every_transfer_equals_the_tables_laws_on_explicit_chains(table, out):
    D, top = boundary(table), valid_from(out) + 6
    checked = 0
    for n in range(valid_from(out), top):
        for m in range(valid_from(out), top):
            want, got = brute_pairs(table, n, m, D), mom_pair(out, n, m)
            for key in set(want) | set(got):
                assert got.get(key, 0.0) == pytest.approx(want.get(key, 0.0), rel=1e-9, abs=1e-30), (key, n, m)
                checked += want.get(key, 0.0) > 0.0
    assert checked


def test_the_reverse_half_of_a_pair_maps_at_its_end_distance_rate(table, out):
    """The inert term (a mid-chain radical next to a chain end takes an H from the bulk of a chain) moves the
    d = 1 pool into the dead-chain pool and the chain into the bulk pool, at the end-distance law's rate."""
    (t,) = [t for k in table["kernels"].values() for t in k["templates"] if t["template"].startswith("end_acceptor")]
    assert t["law"] is None and t["law_inert"]
    names = [x for ts in out["transfer"].values() for x in ts if t["template"] in x["template"].split("@")[0].split("+")]
    assert names and {x["product"] for x in names} == {LABEL}
    owners = {owner for owner, ts in out["transfer"].items() for x in ts if x in names}
    assert owners == {MID + "_d1"}
    assert t["template"] not in [n for i in out["record"]["inert"] for n in i["templates"]]


def brute_moves(table, n, D):
    nbb = table["moment_unit"]["backbone_atoms_per_period"]
    rows = {r["template"]: r for r in table["templates"]}
    out = Counter()
    for k in table["kernels"].values():
        if k["type"] != "intra_h_transfer":
            continue
        for t in k["templates"]:
            if t.get("law") is None and not t.get("law_end_distance"):
                continue
            for (a, b), s in flows(*slot_law(t, None, 1), k["classes"][0], n, nbb, D).items():
                if a != b:
                    out[(a, b)] += rows[t["template"]]["weight"] * s
    return out


def test_every_move_equals_the_tables_laws_on_explicit_chains(table, out):
    D = boundary(table)
    checked = 0
    for n in range(valid_from(out), valid_from(out) + 8):
        want = brute_moves(table, n, D)
        got = Counter()
        for owner, ts in out["move"].items():
            for t in ts:
                got[(owner, t["product"])] += t["k"] * (t["a"] * n + t["b"])
        for key in set(want) | set(got):
            assert got.get(key, 0.0) == pytest.approx(want.get(key, 0.0), rel=1e-9, abs=1e-30), (key, n)
            checked += 1
    assert checked
    assert any(owner.startswith(MID + "_d") and t["product"].startswith(MID + "_d")
               for owner, ts in out["move"].items() for t in ts)        # a shift between fixed classes


# --------------------------------------------------------------------------
# interior beta scission per class
# --------------------------------------------------------------------------
def brute_scission(table, pool, n, D):
    """(rate per chain, {product pool: [E sum L^j per chain-rate unit]}) of interior beta scission on a chain of n
    periods of ``pool``, every position of its class read."""
    nbb = table["moment_unit"]["backbone_atoms_per_period"]
    rows = {r["template"]: r for r in table["templates"]}
    rate, moments_ = Fraction(0), defaultdict(lambda: [Fraction(0)] * 3)
    for t in table["kernels"][BETA]["templates"]:
        law, near = slot_law(t, None, 1)
        w = Fraction(rows[t["template"]]["weight"])
        for x, share in states("radical_interior", n, nbb, D).get(pool, []):
            for outcome, m in outcomes(law, near, n, nbb, x).items():
                rate += w * Fraction(share) * m
                for frag in outcome:
                    q = pool_of(frag[0], frag[2] if len(frag) > 2 else None, D)
                    for j in range(3):
                        moments_[q][j] += w * Fraction(share) * m * Fraction(frag[1]) ** j
    return rate, moments_


def test_interior_beta_scission_maps_per_class_pool_and_equals_the_table(table, out):
    D = boundary(table)
    terms = out["interior_scission"]
    assert set(terms) >= {MID} and any(p.startswith(MID + "_d") for p in terms)
    for pool, ts in terms.items():
        (t,) = ts
        for n in range(t["valid_from_periods"], t["valid_from_periods"] + 10):
            rate, mom = brute_scission(table, pool, n, D)
            assert t["k"] == pytest.approx(float(rate), rel=1e-12), (pool, n)
            got = {p["pool"]: p["moments"] for p in t["products"]}
            assert set(got) == {q for q, v in mom.items() if any(v)}
            for q, q_ji in got.items():
                for j in range(3):
                    want = float(mom[q][j] / rate)
                    assert sum(q_ji[j][i] * n ** i for i in range(j + 1)) == pytest.approx(want, rel=1e-9), \
                        (pool, q, j, n)


def test_a_pool_with_no_scission_feature_gets_no_scission_term(table, out):
    D = boundary(table)
    for d in range(1, D):
        pool = f"{MID}_d{d}"
        rate, _m = brute_scission(table, pool, 30, D)
        assert (pool in out["interior_scission"]) == (rate > 0), pool


# --------------------------------------------------------------------------
# what has no pool, and the bulk's boundary
# --------------------------------------------------------------------------
def test_a_bulk_product_placed_below_the_boundary_is_listed(table, out):
    """The mid-chain radical's own C=C (one bond toward the near end) from a bulk radical at the boundary sits
    at D - 1: carried in the bulk alkene pool, listed."""
    listed = out["record"]["boundary_departures"]
    assert any(x["source"] == MID and x["product"] == INTERNAL for x in listed)
    assert all({"kernel", "templates", "source", "product", "end_distances", "positions_per_chain"} <= set(x)
               for x in listed)


def test_a_join_making_a_mid_chain_feature_is_refused_by_name(table):
    """A join's product feature sits where both reactants meet: no position class carries it."""
    doc = json.loads(json.dumps(table))
    (kid,) = [k for k in doc["kernels"] if k.startswith("disproportionation:acceptor=radical_end")]
    kernel = doc["kernels"][kid]
    kernel["type"], kernel["roles"] = "addition", ["radical", "alkene"]
    for t in kernel["templates"]:
        for slot in ("S0", "S1"):
            law = t["law"][slot]
            law["far"]["outcomes"] = [[1, [["join", "radical_interior", 0]]]]
            law["exact"] = {}
            law["positioned"] = {"far": {"outcomes": [[1, [["join", "radical_interior", 0, {"junction": True}]]]]}}
    terms = mp.kernel_input(doc, spec_of(BASE), [kid], temperature=700.0, phi=1.0, label=LABEL,
                            radical_label=RADICAL, alkene_label=ALKENE, interior_radical_label=MID,
                            interior_alkene_label=INTERNAL, volume_ratio=RATIO, far_law=True, positions=True)
    (refused,) = terms["record"]["position_refused"]
    assert refused["kernel"] == kid and "junction" in refused["reason"]
    assert not terms["joins"]


@pytest.mark.parametrize("drop", ["1", "all"])
def test_a_product_placed_off_the_chain_is_listed_by_name(table, drop):
    """A law read at an end distance its near-end law does not state places its product by the far rule; one
    bond toward the end from end distance 1 is end distance 0, no mid-chain position (I-276: a PP table's
    secondary mid-chain radical, whose near-end law states end distance 2 and skips 1, crashed the producer on
    the missing pool <ui>_d0; with no near-end law at all the boundary is 1 and the product was carried in the
    bulk pool as a boundary departure at end distance 0).  The template is listed in record.unmapped and has
    no term."""
    doc = json.loads(json.dumps(table))
    ((kid, t),) = [(k, t) for k, kernel in doc["kernels"].items() for t in kernel["templates"]
                   if t["template"].startswith("mid_donor:")]
    (slot,) = [s for s, law in t["law"].items() if law["anchor"] == "interior"]
    if drop == "all":
        t.pop("law_near_end")
    else:
        del t["law_near_end"][slot]["at"][drop]
    outcomes = t["law"][slot]["positioned"]["far"]["outcomes"]
    assert any(p[-1].get("offset", 0) < 0 for _m, prods in outcomes for p in prods)
    terms = mp.kernel_input(doc, spec_of(BASE), [kid], temperature=700.0, phi=1.0, label=LABEL,
                            radical_label=RADICAL, alkene_label=ALKENE, interior_radical_label=MID,
                            interior_alkene_label=INTERNAL, volume_ratio=RATIO, far_law=True, positions=True)
    (listed,) = [x for x in terms["record"]["unmapped"] if x["templates"] == [t["template"]]]
    assert listed["kernel"] == kid and "end distance 0" in listed["reason"]
    assert all(t["template"] not in x["template"] for ts in terms["transfer"].values() for x in ts)
    assert all(0 not in x["end_distances"] for x in terms["record"]["boundary_departures"])


@pytest.mark.parametrize("name, slot, kind", [("backbite", None, "move"), ("h_abstractio", "S0", "transfer")])
def test_a_flow_on_short_chains_only_gets_no_term_and_is_listed_with_its_last_length(table, name, slot, kind):
    """An end radical's law with no far-field site (its exact outcomes only, to exact_to_n periods): MoM
    carries no term for it, and the record names the flow and the last length it fires on."""
    doc = json.loads(json.dumps(table))
    ((kid, t),) = [(k, t) for k, kernel in doc["kernels"].items() for t in kernel["templates"]
                   if t["template"].startswith(name + ":")]
    law = t["law"] if slot is None else t["law"][slot]
    law["far"]["outcomes"] = []
    if law.get("positioned"):
        law["positioned"]["far"]["outcomes"] = []
    last = max(int(n) for n, outs in law["exact"].items() if outs)
    terms = mp.kernel_input(doc, spec_of(BASE), [kid], temperature=700.0, phi=1.0, label=LABEL,
                            radical_label=RADICAL, alkene_label=ALKENE, interior_radical_label=MID,
                            interior_alkene_label=INTERNAL, volume_ratio=RATIO, far_law=True, positions=True)
    assert not [x for x in terms[kind].get(RADICAL, []) if t["template"] in x["template"].split("@")[0].split("+")]
    listed = [x for x in terms["record"]["short_lengths"] if x.get("far_law") == 0.0]
    assert listed and all(x["kernel"] == kid and x["templates"] == [t["template"]] for x in listed)
    assert {x["pools"][0] for x in listed} == {RADICAL}
    assert max(x["last_nonzero_periods"] for x in listed) == last


# --------------------------------------------------------------------------
# the MoM input
# --------------------------------------------------------------------------
def test_the_mom_input_declares_every_position_pool_with_its_class_and_its_kernels(table, tmp_path):
    setup = make_setup(tmp_path)
    path = tmp_path / "toy.kernel_rates.json"
    path.write_text(json.dumps(table))
    target = tmp_path / "out"
    manifest = mp.produce(setup["deck"], checkout=CHECKOUT, temperature=700.0, end=40.0, snapshots=9, families=(),
                          catalogs=(), out=target, database=setup["db"], kernel_rates=str(path),
                          kernels=sorted(table["kernels"]), far_law=True, positions=True)
    tree = ast.parse((target / "rmg" / "input.py").read_text())
    pools = {}
    for s in tree.body:
        if isinstance(s, ast.Expr) and isinstance(s.value, ast.Call) and getattr(s.value.func, "id", None) == "polymer":
            kw = {k.arg: ast.literal_eval(k.value) for k in s.value.keywords}
            pools[kw["label"]] = kw
    record = manifest["facts"]["kernels"]
    positions = record["record"]["position_pools"]
    declared = {p for p in positions if p in pools}
    assert declared and declared == {p for p in positions if p in manifest["facts"]["phase_species"]}
    for pool in declared:
        feature = "radical" if positions[pool]["family"] == "radical_interior" else "alkene"
        assert pools[pool]["interior_feature"] == feature
        assert pools[pool]["position_class"] == positions[pool]["position"]
        assert pools[pool].get("kernel_interior_scission") == record["interior_scission"].get(pool)
        assert pools[pool].get("kernel_transfer") == record["transfer"].get(pool)
    assert manifest["facts"]["position_pools"] == {p: [positions[p]["family"], positions[p]["position"]]
                                                   for p in declared}
