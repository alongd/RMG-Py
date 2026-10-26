"""Red tests (I-272): where the interior feature of every law product sits, and
the laws of the terms whose interior participant has sites only near a chain
end (``rmgpy.kmc.kernel_rates``).

The declaration is the toy linear polymer of ``melt_engine_support`` (a
-CH2-CH2- unit, H caps; one CH2 backbone atom per primitive period).  The
catalog is hand-written: the toy H abstraction (a chain donor in the bulk and
at fixed sites near each end) and three context templates, each with sites at
one end distance only: an H transfer onto a mid-chain radical next to a chain
end, a 1,3 H shift of a mid-chain radical three bonds from an end, and a 1,5
H shift of an end radical (its product sits four bonds from the end).

Every expectation is a brute-force enumeration on explicit chains of many
lengths (odd and even, below and past the laws' reach): the melt engine's own
sites of the class on the chain, each fired, and the end distance of the
radical left on the chain read from its backbone.  The table's positioned
laws must give the same sites and positions (``predict_placed``).
"""

import json
from collections import Counter

import numpy
import pytest

import melt_engine_support as S
from meltRunTest import toy_payloads, write_catalog

from rmgpy.kmc import kernel_rates as kr_
from rmgpy.kmc import melt_engine as me
from rmgpy.kmc.window_matcher import AtomEdit, BondEdit, Template, WindowAtom
from rmgpy.polymer_spec import LinearHomopolymerSpec

T = 700.0
MID = S.MID
UNITS = 14                 # brute-force melts: two chains of 28 backbone atoms


def _c(p=0, radical=0):
    return WindowAtom("C", p, radical=radical)


def _h(p=0):
    return WindowAtom("H", p)


def end_acceptor(weight=1.0e-5):
    """H transfer from a CH2 between two CH2 onto a mid-chain radical whose
    other neighbour is the chain-end CH3."""
    atoms = [_c(0), _h(0), _h(0), WindowAtom("C", 0, "port", MID), WindowAtom("C", 0, "port", MID),
             _c(1, radical=1), _h(1), WindowAtom("C", 1, "port", MID), _c(1), _h(1), _h(1), _h(1)]
    bonds = [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (5, 6, 1), (5, 7, 1), (5, 8, 1), (8, 9, 1),
             (8, 10, 1), (8, 11, 1)]
    return Template("end_acceptor", "toy_end_acceptor", atoms, bonds,
                    [BondEdit(0, 1, 0), BondEdit(5, 1, 1), AtomEdit(0, radical=1), AtomEdit(5, radical=0)],
                    weight=weight)


def end_donor(weight=1.0e-5):
    """H transfer onto an end radical from the CH2 next to a chain-end CH3 (a
    chain donor with sites at one edge distance only)."""
    atoms = [_c(0), _h(0), _h(0), WindowAtom("C", 0, "port", MID), _c(0), _h(0), _h(0), _h(0)] + \
        [_c(1, radical=1), _h(1), _h(1), WindowAtom("C", 1, "port", MID)]
    bonds = [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (4, 5, 1), (4, 6, 1), (4, 7, 1), (8, 9, 1),
             (8, 10, 1), (8, 11, 1)]
    return Template("end_donor", "toy_end_donor", atoms, bonds,
                    [BondEdit(0, 1, 0), BondEdit(8, 1, 1), AtomEdit(0, radical=1), AtomEdit(8, radical=0)],
                    weight=weight)


def end_shift(weight=1.0):
    """1,3 H shift: a mid-chain radical three bonds from the chain end takes
    an H of the CH2 next to the end CH3."""
    atoms = [_c(), _h(), _h(), _h(), _c(), _h(), _h(), _c(), _h(), _h(), _c(radical=1), _h(),
             WindowAtom("C", 0, "port", MID)]
    bonds = [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (4, 5, 1), (4, 6, 1), (4, 7, 1), (7, 8, 1),
             (7, 9, 1), (7, 10, 1), (10, 11, 1), (10, 12, 1)]
    return Template("end_shift", "toy_end_shift", atoms, bonds,
                    [BondEdit(4, 5, 0), BondEdit(10, 5, 1), AtomEdit(4, radical=1), AtomEdit(10, radical=0)],
                    weight=weight)


def backbite(weight=1.0):
    """1,5 H shift of an end radical: the radical moves four bonds in."""
    atoms = [_c(radical=1), _h(), _h()]
    bonds = []
    for i in range(4):
        c = len(atoms)
        atoms += [_c(), _h(), _h()]
        bonds += [(c - 3, c, 1), (c, c + 1, 1), (c, c + 2, 1)]
    bonds += [(0, 1, 1), (0, 2, 1)]
    atoms.append(WindowAtom("C", 0, "port", MID))
    bonds.append((12, 15, 1))
    return Template("backbite", "toy_backbite", atoms, bonds,
                    [BondEdit(12, 13, 0), BondEdit(0, 13, 1), AtomEdit(12, radical=1), AtomEdit(0, radical=0)],
                    weight=weight)


def templates():
    return [t for t in S.pe_templates() if t.name == "h_abstraction"] + \
        [end_acceptor(), end_donor(), end_shift(), backbite()]


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    path = write_catalog(tmp_path_factory.mktemp("position") / "catalog.json", toy_payloads(templates()),
                         tag="position-toy", bind=S.pe_mapping({1: 1}))
    doc = json.loads(path.read_text())
    spec = LinearHomopolymerSpec.from_mapping(S.pe_mapping({10: 1}))
    table = kr_.kernel_rates(doc, spec, T)
    catalog = me.PayloadCatalog(doc["payloads"])
    engine = me.MeltEngine(kr_._single_chain_spec(spec, UNITS), catalog, T, numpy.random.default_rng(0))
    probes = kr_._Probes(spec, engine, 2, 1, UNITS)
    return {"table": table, "engine": engine, "probes": probes, "nbb": table["moment_unit"]["backbone_atoms_per_period"]}


def term(table, prefix):
    (t,) = [t for k in table["kernels"].values() for t in k["templates"] if t["template"].startswith(prefix[:12])]
    return t


def class_index(engine, prefix):
    (ci,) = [i for i, t in enumerate(engine.templates) if t.name.startswith(prefix[:12])]
    return ci


# --------------------------------------------------------------------------
# brute force: the engine's own sites on explicit chains, fired one by one
# --------------------------------------------------------------------------
def _radicals_on(state, seg, cids):
    """Counter {(backbone atoms, end distance)} of the radical atoms on the products holding ``seg``'s atoms."""
    out = Counter()
    for cid in cids:
        atoms = set(state.components.atoms(cid))
        mine = [a for a in seg if a in atoms]
        if not mine:
            continue
        path = kr_._walk(state, [a for a in atoms if state.origin[a].kind == "chain" and state.origin[a].is_backbone])
        for i, a in enumerate(path):
            if state.atoms[a].radical:
                out[(len(path), min(i, len(path) - 1 - i))] += 1
    return out


def brute(toy, ci, p, family, n, x=None, partner=None):
    """(sites, Counter {(backbone atoms, end distance): count} of the radicals the sites leave on their own
    product) of participant ``p`` of class ``ci`` on the probe (family, n periods, feature at backbone atom x),
    each site fired once (with the first partner site for a bimolecular class)."""
    engine, probes = toy["engine"], toy["probes"]
    state, paths = probes.melt()
    seg, _feature = probes.put(state, paths[0], family, n, x)
    c0 = state.components.component_of(seg[0])
    cids = [c0]
    if partner is not None:
        pseg, _pf = probes.put(state, paths[1], *partner)
        c1 = state.components.component_of(pseg[0])
        cids.append(c1)
    entries = engine._scan(cids, state=state)
    saved, engine.state = engine.state, state
    out, sites = Counter(), 0
    try:
        if partner is None:
            for tr in kr_._uni_all(entries[c0]).get(ci, ((), 0))[0]:
                sites += 1
                tx, made = engine._begin_apply(tr, me._LocalView(state, [c0]))
                try:
                    out.update(_radicals_on(state, seg, made))
                finally:
                    tx.rollback()
        else:
            mine = [s for (cj, q, _k), ss in sorted(entries[c0].sites.items(), key=repr) if cj == ci and q == p
                    for s in ss]
            theirs = [s for (cj, q, _k), ss in sorted(entries[c1].sites.items(), key=repr)
                      if cj == ci and q == 1 - p for s in ss]
            for s in mine:
                sites += 1
                args = (ci, c0, s, c1, theirs[0]) if p == 0 else (ci, c1, theirs[0], c0, s)
                tr, view = engine._pair_transition(*args)
                tx, made = engine._begin_apply(tr, view)
                try:
                    out.update(_radicals_on(state, seg, made))
                finally:
                    tx.rollback()
    finally:
        engine.state = saved
    return sites, out


def stated(law, n, nbb, x=None):
    """(sites, Counter {(backbone atoms, end distance)} of the interior products) of a positioned law at n
    periods (an end-distance law: its reactant's feature at end distance x)."""
    got = kr_.predict_placed(law, n, nbb, x)
    sites = sum(got.values())
    out = Counter()
    for outcome, m in got.items():
        for product in outcome:
            if product[0] in kr_.INTERIOR_FAMILIES:
                out[(product[1] * nbb, product[2])] += m
    return sites, out


LENGTHS = (4, 5, 7, 8, 11, 13, 16, 19, 22, 25)


# --------------------------------------------------------------------------
# positions on the products of laws that exist today
# --------------------------------------------------------------------------
def _far_positions(far, rows):
    return [p[-1] for _m, prods in rows for p in prods if p[0] == "radical_interior"]


def test_a_chain_donor_in_the_bulk_states_its_radicals_end_distance_range(toy):
    far = term(toy["table"], "h_abstraction")["law"]["S1"]["positioned"]["far"]
    bulk = _far_positions(far, far["bulk"])
    assert bulk and all(set(pos) == {"bulk"} for pos in bulk)


def test_a_chain_donor_at_an_edge_site_states_its_radicals_end_distance(toy):
    far = term(toy["table"], "end_donor")["law"]["S1"]["positioned"]["far"]
    assert far["bulk"] == []
    edges = [pos for side in ("first", "last") for rows in far[side].values() for pos in _far_positions(far, rows)]
    assert edges == [{"end_distance": 1}, {"end_distance": 1}]


@pytest.mark.parametrize("name, donor", [("h_abstraction", 1), ("end_donor", 0)])
def test_brute_force_finds_sites_and_radicals(toy, name, donor):
    """Every length past the shortest has donor sites, each leaving one radical on the chain."""
    ci = class_index(toy["engine"], name)
    found = [brute(toy, ci, donor, "chain", n, partner=("radical_end", 6, None)) for n in LENGTHS]
    assert all(sites > 0 and sum(radicals.values()) == sites for sites, radicals in found[1:])


@pytest.mark.parametrize("n", LENGTHS)
@pytest.mark.parametrize("name, donor", [("h_abstraction", 1), ("end_donor", 0)])
def test_the_chain_donor_positions_equal_a_brute_force_enumeration(toy, name, donor, n):
    t = term(toy["table"], name)
    ci = class_index(toy["engine"], name)
    assert stated(t["law"]["S1"], n, toy["nbb"]) == \
        brute(toy, ci, donor, "chain", n, partner=("radical_end", 6, None))


@pytest.mark.parametrize("n", LENGTHS)
def test_a_backbite_product_sits_at_the_end_distance_brute_force_finds(toy, n):
    t = term(toy["table"], "backbite")
    law = t["law"]
    assert law["anchor"] == "end"
    assert [p[-1] for _m, prods in law["positioned"]["far"]["outcomes"] for p in prods] == [{"end_distance": 4}]
    ci = class_index(toy["engine"], "backbite")
    assert stated(law, n, toy["nbb"]) == brute(toy, ci, 0, "radical_end", n)


# --------------------------------------------------------------------------
# the terms kernel_rates marks law_inert: a law per end distance
# --------------------------------------------------------------------------
def test_an_inert_term_keeps_its_old_keys_and_gains_an_end_distance_law(toy):
    t = term(toy["table"], "end_acceptor")
    assert t["law"] is None and t["law_inert"] == "an interior participant with no far-field site"
    acceptor = t["law_end_distance"]["S0"]
    assert acceptor["anchor"] == "end_distance" and acceptor["class"] == "radical_interior:CH"
    assert sorted(acceptor["at"]) == ["1"]
    assert t["law_end_distance"]["S1"]["anchor"] == "chain"


@pytest.mark.parametrize("n", LENGTHS)
def test_the_end_acceptor_sites_at_every_position_equal_a_brute_force_enumeration(toy, n):
    t = term(toy["table"], "end_acceptor")
    ci = class_index(toy["engine"], "end_acceptor")
    law = t["law_end_distance"]["S0"]
    for x in range(1, n * toy["nbb"] - 1):
        d = min(x, n * toy["nbb"] - 1 - x)
        assert stated(law, n, toy["nbb"], d) == brute(toy, ci, 1, "radical_interior", n, x, ("chain", 12, None)), x


@pytest.mark.parametrize("n", LENGTHS)
def test_the_end_acceptors_donor_positions_equal_a_brute_force_enumeration(toy, n):
    t = term(toy["table"], "end_acceptor")
    ci = class_index(toy["engine"], "end_acceptor")
    assert stated(t["law_end_distance"]["S1"], n, toy["nbb"]) == \
        brute(toy, ci, 0, "chain", n, partner=("radical_interior", 12, 1))


@pytest.mark.parametrize("n", LENGTHS)
def test_an_interior_shift_moves_the_radical_where_brute_force_finds_it(toy, n):
    t = term(toy["table"], "end_shift")
    assert t["law"] is None and t["law_inert"]
    law = t["law_end_distance"]
    assert sorted(law["at"]) == ["3"]
    assert [p[-1] for _m, prods in law["at"]["3"]["outcomes"] for p in prods] == [{"end_distance": 1}]
    ci = class_index(toy["engine"], "end_shift")
    for x in range(1, n * toy["nbb"] - 1):
        d = min(x, n * toy["nbb"] - 1 - x)
        assert stated(law, n, toy["nbb"], d) == brute(toy, ci, 0, "radical_interior", n, x), x


def test_an_end_distance_law_states_its_rate_coefficient_per_radical(toy):
    t = term(toy["table"], "end_shift")
    row = next(r for r in toy["table"]["templates"] if r["template"] == t["template"])
    at = t["law_end_distance"]["at"]["3"]
    assert at["k_per_feature"] == pytest.approx(row["weight"] * at["sites"], rel=1e-12)
