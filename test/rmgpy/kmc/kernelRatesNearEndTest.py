"""Red tests (I-273): the near-end law of every interior participant of a
law term (``rmgpy.kmc.kernel_rates``, ``law_near_end``).

A law term's interior participant (a mid-chain radical, an internal alkene)
has its far-field count and outcomes in the middle of a chain; near an end
some of its sites (or all) are missing, and the ones left need not be the far
law's at that position (a beta scission keeps the bond away from the end, a
donor the H away from it).  The table states, per such participant, every end
distance where the count or the products differ from the far law's (``at``,
zero sites included), and every position of the chains shorter than where
that map holds (``exact``); ``predict_near`` reads it, the far law elsewhere.

The declaration is the toy linear polymer of ``melt_engine_support`` (a
-CH2-CH2- unit, H caps).  The catalog is generic: the toy mid-chain-radical
beta scission, a mid-chain radical giving up the H next to it to an end
radical, and an internal alkene taking an H from a chain.  Every expectation
is a brute-force enumeration on explicit chains of many lengths with the
feature at every position: the melt engine's own sites of the class, each
fired, and the features (radical, C=C) on the products of the probed molecule
read from their backbones.
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
LENGTHS = (4, 5, 7, 8, 11, 13, 16, 19, 22, 25)


def mid_donor(weight=2.0e-5):
    """An end radical takes the H of a CH2 next to a mid-chain radical (the
    radical's other neighbour and the CH2's are chain CH2)."""
    atoms = S._ch2_radical(0) + [WindowAtom("C", 1, radical=1), WindowAtom("H", 1), WindowAtom("C", 1, "port", MID),
                                 WindowAtom("C", 1), WindowAtom("H", 1), WindowAtom("H", 1),
                                 WindowAtom("C", 1, "port", MID)]
    bonds = [(0, 1, 1), (0, 2, 1), (0, 3, 1), (4, 5, 1), (4, 6, 1), (4, 7, 1), (7, 8, 1), (7, 9, 1), (7, 10, 1)]
    return Template("mid_donor", "toy_mid_donor", atoms, bonds,
                    [BondEdit(7, 8, 0), BondEdit(0, 8, 1), BondEdit(4, 7, 2), AtomEdit(0, radical=0),
                     AtomEdit(4, radical=0)], weight=weight)


def alkene_acceptor(weight=3.0e-9):
    """An internal C=C (each carbon one H and a chain CH2 beyond) takes an H
    of a chain CH2 between two CH2: two mid-chain radicals."""
    atoms = [WindowAtom("C"), WindowAtom("C"), WindowAtom("H"), WindowAtom("C", 0, "port", MID), WindowAtom("H"),
             WindowAtom("C", 0, "port", MID),
             WindowAtom("C", 1), WindowAtom("H", 1), WindowAtom("H", 1), WindowAtom("C", 1, "port", MID),
             WindowAtom("C", 1, "port", MID)]
    bonds = [(0, 1, 2), (0, 2, 1), (0, 3, 1), (1, 4, 1), (1, 5, 1), (6, 7, 1), (6, 8, 1), (6, 9, 1), (6, 10, 1)]
    return Template("alkene_acceptor", "toy_alkene_acceptor", atoms, bonds,
                    [BondEdit(6, 7, 0), BondEdit(0, 7, 1), BondEdit(0, 1, 1), AtomEdit(6, radical=1),
                     AtomEdit(1, radical=1)], weight=weight)


def templates():
    return [t for t in S.pe_templates() if t.name == "beta_scission"] + [mid_donor(), alkene_acceptor()]


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    path = write_catalog(tmp_path_factory.mktemp("near_end") / "catalog.json", toy_payloads(templates()),
                         tag="near-end-toy", bind=S.pe_mapping({1: 1}))
    doc = json.loads(path.read_text())
    spec = LinearHomopolymerSpec.from_mapping(S.pe_mapping({10: 1}))
    table = kr_.kernel_rates(doc, spec, T)
    catalog = me.PayloadCatalog(doc["payloads"])
    engine = me.MeltEngine(kr_._single_chain_spec(spec, UNITS), catalog, T, numpy.random.default_rng(0))
    probes = kr_._Probes(spec, engine, 2, 1, UNITS)
    return {"table": table, "engine": engine, "probes": probes, "nbb": table["moment_unit"]["backbone_atoms_per_period"]}


def term_of(table, name):
    (kid, t), = [(kid, t) for kid, k in table["kernels"].items() for t in k["templates"]
                 if t["template"].startswith(name[:12])]
    return kid, t


def class_index(engine, name):
    (ci,) = [i for i, t in enumerate(engine.templates) if t.name.startswith(name[:12])]
    return ci


def interior_slots(table, kid, t):
    """[(slot key or None for one participant, kernel slot index, class)] of the term's interior participants."""
    classes = table["kernels"][kid]["classes"]
    return [(None if len(classes) == 1 else f"S{i}", i, c) for i, c in enumerate(classes)
            if kr_.family_of(c) in kr_.INTERIOR_FAMILIES]


def template_index(table, name, cls):
    """The participant index of class ``cls`` in the template's own order (the engine's)."""
    (row,) = [r for r in table["templates"] if r["template"].startswith(name[:12])]
    (p,) = [p for p, part in enumerate(row["participants"]) if cls in part]
    return p


# --------------------------------------------------------------------------
# brute force: the engine's own sites on explicit chains, fired one by one
# --------------------------------------------------------------------------
def _features_on(state, seg, cids):
    """Counter {(feature, backbone atoms, end distance)} on the products holding ``seg``'s atoms: every radical
    atom (end distance min(i, L - 1 - i)) and backbone C=C (min(i, L - 2 - i) for bond (i, i + 1))."""
    out = Counter()
    for cid in cids:
        atoms = set(state.components.atoms(cid))
        if not any(a in atoms for a in seg):
            continue
        path = kr_._walk(state, [a for a in atoms if state.origin[a].kind == "chain" and state.origin[a].is_backbone])
        L = len(path)
        for i, a in enumerate(path):
            if state.atoms[a].radical:
                out[("radical", L, min(i, L - 1 - i))] += 1
            if i + 1 < L and state.atoms[a].edges.get(path[i + 1]) == 2.0:
                out[("alkene", L, min(i, L - 2 - i))] += 1
    return out


def brute(toy, ci, p, family, n, x, partner=None):
    """(sites, features on the probed molecule's products) of participant ``p`` of class ``ci`` on the probe
    (family, n periods, feature at backbone atom x), each site fired once (with the first partner site)."""
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
                    out.update(_features_on(state, seg, made))
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
                    out.update(_features_on(state, seg, made))
                finally:
                    tx.rollback()
    finally:
        engine.state = saved
    return sites, out


FEATURE = {"radical_end": "radical", "radical_interior": "radical", "alkene_terminal": "alkene",
           "alkene_internal": "alkene"}


def stated(law, near, n, nbb, x):
    """(sites, Counter {(feature, backbone atoms, end distance)}) of a far-field interior law with its near-end
    law at end distance x on n periods (``kr_.predict_near``); an end feature sits at end distance 0."""
    got = kr_.predict_near(law, near, n, nbb, x)
    out = Counter()
    for outcome, m in got.items():
        for product in outcome:
            if product[0] in FEATURE:
                out[(FEATURE[product[0]], product[1] * nbb, product[2] if len(product) > 2 else 0)] += m
    return sum(got.values()), out


PARTNER = {"radical_end": ("radical_end", 6, None), "chain": ("chain", 12, None)}


def cases(table):
    """(name, kernel, term, slot, kernel slot, class) of every interior participant of a law term."""
    out = []
    for name in ("beta_scission", "mid_donor", "alkene_acceptor"):
        kid, t = term_of(table, name)
        for slot, i, cls in interior_slots(table, kid, t):
            out.append((name, kid, t, slot, i, cls))
    return out


# --------------------------------------------------------------------------
# the table
# --------------------------------------------------------------------------
def test_every_interior_participant_of_a_law_term_states_its_near_end_law(toy):
    found = cases(toy["table"])
    assert {c[0] for c in found} == {"beta_scission", "mid_donor", "alkene_acceptor"}
    for name, kid, t, slot, i, cls in found:
        assert t["law"] is not None, name
        near = t["law_near_end"] if slot is None else t["law_near_end"][slot]
        assert near["anchor"] == "end_distance" and near["class"] == cls, name
        assert int(near["from_n"]) >= 1 and near["at"], name


def test_near_end_laws_include_partial_and_zero_sites(toy):
    """Not vacuous: some end distance keeps fewer sites than the far field, and some keeps none."""
    rows = {r["template"]: r for r in toy["table"]["templates"]}
    partial, zero = [], []
    for name, kid, t, slot, i, cls in cases(toy["table"]):
        near = t["law_near_end"] if slot is None else t["law_near_end"][slot]
        per = rows[t["template"]]["participants"][template_index(toy["table"], name, cls)][cls]["per_feature"]
        for d, entry in near["at"].items():
            (partial if 0 < entry["sites"] < per else zero if entry["sites"] == 0 else []).append((name, d))
    assert partial and zero, (partial, zero, rows.keys())


def test_a_unimolecular_near_end_entry_states_its_rate_coefficient(toy):
    _kid, t = term_of(toy["table"], "beta_scission")
    row = next(r for r in toy["table"]["templates"] if r["template"] == t["template"])
    for entry in t["law_near_end"]["at"].values():
        assert entry["k_per_feature"] == pytest.approx(row["weight"] * entry["sites"], rel=1e-12)


@pytest.mark.parametrize("n", LENGTHS)
@pytest.mark.parametrize("name", ["beta_scission", "mid_donor", "alkene_acceptor"])
def test_the_sites_and_products_at_every_position_equal_a_brute_force_enumeration(toy, name, n):
    table, nbb = toy["table"], toy["nbb"]
    ci = class_index(toy["engine"], name)
    checked = 0
    for _name, kid, t, slot, i, cls in [c for c in cases(table) if c[0] == name]:
        family = kr_.family_of(cls)
        law = t["law"] if slot is None else t["law"][slot]
        near = t["law_near_end"] if slot is None else t["law_near_end"][slot]
        p = template_index(table, name, cls)
        partner = None
        if slot is not None:
            other = table["kernels"][kid]["classes"][1 - i]
            partner = PARTNER[kr_.family_of(other)]
        L = n * nbb
        for x in toy["probes"].positions(family, n):
            if x is None or not toy["probes"].carries(cls, n, x):
                continue
            assert stated(law, near, n, nbb, x) == brute(toy, ci, p, family, n, x, partner), (name, slot, n, x)
            checked += 1
    if n >= 7:
        assert checked
