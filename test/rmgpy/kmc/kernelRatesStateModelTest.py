"""Red tests (I-279): the kernel-rates table's state model on a declaration
whose chain is not its own reverse and whose period has two backbone atoms
(``rmgpy.kmc.kernel_rates``, ``kmc-kernel-rates/4``).

The fixture is a reduced PP catalog: the payloads of five templates of the
final PP catalog (pk-i275 restage) with their reciprocal twins and the
primaries their derived templates name, bound to the PP deck.  PP's chain
H-[CH2-CH(CH3)]n-H has an isopropyl head and an n-propyl tail, so:

* a backbone cut leaves a CH end on the head-side piece and a CH2 end on the
  tail-side one (both sides of a cut are states);
* the backbone positions read from the head are not those read from the tail
  (positions are oriented, never folded to the nearer end);
* a cut inside a period (CH2-CH) leaves pieces of an odd number of backbone
  atoms: H-[c1 c2]k-c1-H and H-c2-[c1 c2]k-H, two more chain types;
* along a chain the classes alternate (CH2 and CH positions), so a product
  code names its class and a chain law's bulk is periodic in the anchor.

Every expectation is read independently of the table builder: the melt
engine's own sites on explicit pieces, each fired, and each product read from
its atoms' declaration keys (this module's ``describe``, which knows PP: the
head key is c1, the tail key c2).  PE (the toy of ``melt_engine_support``,
one CH2 per period) is the symmetric control: its table keeps the
``kmc-kernel-rates/3`` state model, with no key of the new one.
"""

import gzip
import json
from collections import Counter
from pathlib import Path

import numpy
import pytest

import melt_engine_support as S
from meltRunTest import toy_payloads, write_catalog

from rmgpy.kmc import kernel_rates as kr_
from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import run as kr
from rmgpy.polymer_spec import LinearHomopolymerSpec

T = 700.0
FIXTURES = Path(__file__).parent / "fixtures"
HOMOLYSIS = "6ff3c29fee36"          # R_Recombination reverse: every backbone bond of a chain
END_BACKBITE = "7966d7d30b15"       # intra_H_migration: the CH2 end (tail-side piece) takes an H
TERTIARY_BACKBITE = "c2d6a938bc13"  # intra_H_migration: a tertiary mid-chain radical moves
CHAIN_DONOR = "503839880460"        # H_Abstraction: a CH end radical takes an H of a chain
BETA = "cb700a132726"               # R_Addition_MultipleBond reverse: a tertiary radical's beta scission
RECOMBINATION = "add96b9d7950"      # R_Recombination forward: CH end + CH2 end


def pp_catalog():
    return json.loads(gzip.decompress((FIXTURES / "i279" / "pp_state_model_catalog.json.gz").read_bytes()))


def pp_spec():
    return kr.read_deck(FIXTURES / "decks" / "pp_deck.py").spec


@pytest.fixture(scope="module")
def pp():
    doc = pp_catalog()
    spec = pp_spec()
    table = kr_.kernel_rates(doc, spec, T)
    n_hi = table["laws_exact_to_n"]
    engine = me.MeltEngine(kr_._single_chain_spec(spec, n_hi), me.PayloadCatalog(doc["payloads"]), T,
                           numpy.random.default_rng(0))
    probes = kr_._Probes(spec, engine, 1, 2, 2 * n_hi + 3)
    return {"table": table, "engine": engine, "probes": probes, "spec": spec, "n_hi": n_hi,
            "types": table.get("chain_types")}


@pytest.fixture(scope="module")
def pe(tmp_path_factory):
    path = write_catalog(tmp_path_factory.mktemp("state_model_pe") / "catalog.json", toy_payloads(S.pe_templates()),
                         tag="state-model-pe", bind=S.pe_mapping({1: 1}))
    return kr_.kernel_rates(json.loads(path.read_text()), LinearHomopolymerSpec.from_mapping(S.pe_mapping({10: 1})), T)


def terms(table, prefix):
    return [(kid, t) for kid, k in table["kernels"].items() for t in k["templates"] if t["template"].startswith(prefix)]


def class_index(engine, prefix, key=None):
    return [i for i, t in enumerate(engine.templates)
            if t.name.startswith(prefix) and (key is None or t.name.endswith(key))]


# --------------------------------------------------------------------------
# the independent reader: a product from its atoms' declaration keys (PP)
# --------------------------------------------------------------------------
def _h(state, a):
    return sum(1 for b in state.atoms[a].edges if state.atoms[b].element == "H")


def _lab(state, atoms):
    return "=".join(state.atoms[a].element + ("H" if _h(state, a) == 1 else f"H{_h(state, a)}" if _h(state, a) else "")
                    for a in atoms)


def _rank(state, a):
    """The atom's place on its own chain, head first: (chain ordinal, unit, key)."""
    o = state.origin[a]
    return o.ordinal, int(o.segment.split()[1]), {"c1": 0, "c2": 1}[o.key]


def describe(state, cid):
    """(label, n[, position kind, d]) of a product component, read from its atoms: its backbone path read head to
    tail (each chain's run in increasing unit order; runs in opposite directions are a defect: None), its chain type
    from the keys of its two end atoms (PP: c1~c2, c2~c1 = c1~c2 read backwards, c1~c1, c2~c2), its feature."""
    atoms = set(state.components.atoms(cid))
    bb = [a for a in atoms if state.origin[a].kind == "chain" and state.origin[a].is_backbone]
    if not bb:
        return None
    path = kr_._walk(state, bb)
    runs = []
    for a in path:
        if runs and state.origin[runs[-1][-1]].ordinal == state.origin[a].ordinal:
            runs[-1].append(a)
        else:
            runs.append([a])
    ups = {_rank(state, r[1]) > _rank(state, r[0]) for r in runs if len(r) > 1}
    if len(ups) > 1:
        return None
    if ups == {False}:
        path = path[::-1]
    keys = [state.origin[a].key for a in path]
    L = len(path)
    name = f"{keys[0]}~{keys[-1]}"
    if name == "c2~c1":
        path, name = path[::-1], "c1~c2"
    symmetric = name != "c1~c2"
    suffix = "" if name == "c1~c2" else f"/{name}"
    rad = [i for i, a in enumerate(path) if state.atoms[a].radical]
    pis = [i for i in range(L - 1) if state.atoms[path[i]].edges.get(path[i + 1]) == 2.0]
    assert len(rad) + len(pis) <= 1
    n = L // 2
    if not rad and not pis:
        return ("chain" + suffix, n)
    i = (rad or pis)[0]
    w = 1 if rad else 2
    if symmetric and i > L - w - i:
        path, i = path[::-1], L - w - i
    if rad:
        if i in (0, L - 1):
            return (f"radical_end:{_lab(state, [path[i]])}" + suffix, n)
        family, feature = "radical_interior", [path[i]]
    else:
        if i == 0 or i + 1 == L - 1:
            end = [path[0], path[1]] if i == 0 else [path[L - 1], path[L - 2]]
            return (f"alkene_terminal:{_lab(state, end)}" + suffix, n)
        family, feature = "alkene_internal", [path[i], path[i + 1]]
    label = f"{family}:{_lab(state, feature)}" + suffix
    if symmetric:
        return (label, n, "end_distance", i)
    return (label, n) + (("head_distance", i) if i <= L - w - i else ("tail_distance", L - w - i))


def brute(pp, ci, p, family, n, x=None, tau=None, partner=None):
    """Counter of the sorted product descriptions of every site of participant p of class ci on the piece
    (family, n, x, tau), each fired once (bimolecular: with the first partner site, the participant's own product
    only)."""
    engine, probes = pp["engine"], pp["probes"]
    state, paths = probes.melt()
    seg, _feature = probes.put(state, paths[0], family, n, x, tau=tau)
    c0 = state.components.component_of(seg[0])
    cids = [c0]
    if partner is not None:
        pseg, _pf = probes.put(state, paths[1], *partner)
        c1 = state.components.component_of(pseg[0])
        cids.append(c1)
    entries = engine._scan(cids, state=state)
    saved, engine.state = engine.state, state
    out = Counter()
    mine = set(seg)
    try:
        if partner is None:
            fires = [(tr, me._LocalView(state, [c0])) for tr in kr_._uni_all(entries[c0]).get(ci, ((), 0))[0]]
        else:
            own = [s for (cj, q, _k), ss in sorted(entries[c0].sites.items(), key=repr) if cj == ci and q == p
                   for s in ss]
            theirs = [s for (cj, q, _k), ss in sorted(entries[c1].sites.items(), key=repr)
                      if cj == ci and q == 1 - p for s in ss]
            fires = [engine._pair_transition(*((ci, c0, s, c1, theirs[0]) if p == 0 else (ci, c1, theirs[0], c0, s)))
                     for s in own]
        for tr, view in fires:
            tx, made = engine._begin_apply(tr, view)
            try:
                got = []
                for cid in made:
                    atoms = set(state.components.atoms(cid))
                    if partner is not None and not atoms & mine:
                        continue
                    d = describe(state, cid)
                    if d is not None:
                        got.append(d)
                out[tuple(sorted(got, key=repr))] += 1
            finally:
                tx.rollback()
    finally:
        engine.state = saved
    return out


def stated(pp, law, n, x=None):
    """Counter of the sorted product tuples a law states at n periods (positioned when it has positions)."""
    types = pp["types"]
    got = kr_.predict_placed(law, n, 2, x, types=types) if law.get("positioned") else kr_.predict(law, n, 2,
                                                                                                    types=types)
    return Counter({tuple(sorted(k, key=repr)): v for k, v in got.items()})


def lengths(pp):
    return [3, 4, 7, pp["n_hi"] - 1, pp["n_hi"] + 1, pp["n_hi"] + 4]


# --------------------------------------------------------------------------
# chain types
# --------------------------------------------------------------------------
def test_the_table_names_its_chain_types_with_their_symmetry(pp):
    assert pp["table"]["schema"] == "kmc-kernel-rates/4"
    types = pp["types"]
    assert sorted(types) == ["c1~c1", "c1~c2", "c2~c2"]
    assert types["c1~c2"]["declared"] and not types["c1~c2"]["symmetric"]
    assert (types["c1~c2"]["start_offset"], types["c1~c2"]["residue"]) == (0, 0)
    assert (types["c1~c1"]["start_offset"], types["c1~c1"]["residue"], types["c1~c1"]["symmetric"]) == (0, 1, True)
    assert (types["c2~c2"]["start_offset"], types["c2~c2"]["residue"], types["c2~c2"]["symmetric"]) == (1, 1, True)


def test_a_sub_period_piece_has_its_own_mass_offset(pp):
    """H-[c1 c2]n-c1-H holds one CH2 more than n periods and the two H caps."""
    types = pp["types"]
    chain = pp["table"]["families"]["chain"]["molar_mass_offset_g_mol"]
    aa = types["c1~c1"]["families"]["chain"]["molar_mass_offset_g_mol"]
    assert aa - chain == pytest.approx(14.02658, abs=1e-3)


# --------------------------------------------------------------------------
# rule 3: homolysis cuts inside a period are states, not refusals
# --------------------------------------------------------------------------
def test_homolysis_has_a_law_on_every_chain_type(pp):
    got = {kid: t for kid, t in terms(pp["table"], HOMOLYSIS)}
    assert sorted(got) == ["homolysis:substrate=chain", "homolysis:substrate=chain/c1~c1",
                           "homolysis:substrate=chain/c2~c2"]
    for kid, t in got.items():
        assert t["law"] is not None, (kid, t.get("law_refused"))


def test_homolysis_bulk_is_periodic_in_the_anchor(pp):
    """A cut of c1-c2 makes c1~c1 and c2~c2 pieces, a cut of c2-c1 two c1~c2 pieces."""
    (t,) = [t for kid, t in terms(pp["table"], HOMOLYSIS) if kid == "homolysis:substrate=chain"]
    far = t["law"]["far"]
    assert "bulk" not in far and sorted(far["bulk_phases"]) == ["0", "1"]
    tokens = {ph: sorted(p[0] for _m, prods in rows for p in prods) for ph, rows in far["bulk_phases"].items()}
    assert tokens == {"0": ["radical_end:CH/c2~c2", "radical_end:CH2/c1~c1"],
                      "1": ["radical_end:CH", "radical_end:CH2"]}


@pytest.mark.parametrize("tau", [None, "c1~c1", "c2~c2"])
def test_homolysis_products_equal_a_brute_force_enumeration(pp, tau):
    kid = "homolysis:substrate=chain" + ("" if tau is None else f"/{tau}")
    (t,) = [t for k, t in terms(pp["table"], HOMOLYSIS) if k == kid]
    (ci,) = class_index(pp["engine"], HOMOLYSIS)
    for n in lengths(pp):
        assert stated(pp, t["law"], n) == brute(pp, ci, 0, "chain", n, tau=tau), (kid, n)


# --------------------------------------------------------------------------
# rule 1: the tail-side end of a cut is probed
# --------------------------------------------------------------------------
def test_the_tail_side_end_radical_has_a_law_read_at_its_own_end(pp):
    ((kid, t),) = terms(pp["table"], END_BACKBITE)
    assert kid == "intra_h_transfer:substrate=radical_end:CH2"
    assert t["law"] is not None, t.get("law_refused")
    assert t["law"]["anchor"] == "end" and t["law"]["side"] == "start"


def test_the_tail_side_end_backbite_equals_a_brute_force_enumeration(pp):
    ((_kid, t),) = terms(pp["table"], END_BACKBITE)
    (ci,) = class_index(pp["engine"], END_BACKBITE)
    for n in lengths(pp):
        assert stated(pp, t["law"], n) == brute(pp, ci, 0, "radical_end", n, "start"), n


def test_the_end_probes_of_both_sides_carry_their_own_class(pp):
    probes = pp["probes"]
    state, paths = probes.melt()
    _seg, head = probes.put(state, paths[0], "radical_end", 6)
    _seg, tail = probes.put(state, paths[1], "radical_end", 6, "start")
    assert (kr_._label(state, head), kr_._label(state, tail)) == ("CH", "CH2")


# --------------------------------------------------------------------------
# rule 2: interior positions are oriented
# --------------------------------------------------------------------------
def test_the_tertiary_backbite_has_a_law(pp):
    got = [(kid, t) for kid, t in terms(pp["table"], TERTIARY_BACKBITE) if t["template"].endswith("837b751137c8")]
    assert got
    for kid, t in got:
        assert t["law"] is not None, (kid, t.get("law_refused"))


def test_ab_positions_are_not_folded(pp):
    probes = pp["probes"]
    assert probes.positions("radical_interior", 10) == list(range(1, 19))
    assert probes.positions("radical_interior", 10, tau="c1~c1") == list(range(1, 10))


def test_the_tertiary_backbite_at_every_position_equals_a_brute_force_enumeration(pp):
    probes = pp["probes"]
    checked = 0
    for kid, t in terms(pp["table"], TERTIARY_BACKBITE):
        if not t["template"].endswith("837b751137c8") or t["law"] is None:
            continue
        law = t["law"]
        tau = law["class"].split("/")[1] if "/" in law["class"] else None
        (ci,) = class_index(pp["engine"], TERTIARY_BACKBITE, "837b751137c8")
        for n in (pp["n_hi"] + 1, pp["n_hi"] + 3):
            for x in probes.positions("radical_interior", n, tau=tau):
                if not probes.carries(law["class"], n, x):
                    continue
                want = brute(pp, ci, 0, "radical_interior", n, x, tau=tau)
                got = Counter({tuple(sorted(k, key=repr)): v
                               for k, v in kr_.predict_near(law, t.get("law_near_end"), n, 2, x,
                                                            types=pp["types"]).items()})
                assert got == want, (kid, n, x)
                checked += 1
    assert checked


# --------------------------------------------------------------------------
# rule 4: codes carry the class; a chain law's bulk is periodic
# --------------------------------------------------------------------------
def test_a_chain_donor_names_the_class_it_makes_per_phase(pp):
    got = [(kid, t) for kid, t in terms(pp["table"], CHAIN_DONOR) if kid.endswith("+donor=chain")]
    assert got
    for kid, t in got:
        assert t["law"] is not None, (kid, t.get("law_refused"))
        far = t["law"]["S1"]["far"]
        tokens = sorted({p[0] for rows in far["bulk_phases"].values() for _m, prods in rows for p in prods})
        assert tokens == ["radical_interior:C", "radical_interior:CH"], kid


def test_a_chain_donor_equals_a_brute_force_enumeration(pp):
    checked = 0
    for kid, t in terms(pp["table"], CHAIN_DONOR):
        if not kid.endswith("+donor=chain") or t["law"] is None:
            continue
        (ci,) = class_index(pp["engine"], CHAIN_DONOR, t["template"].split(":")[1])
        for n in lengths(pp):
            assert stated(pp, t["law"]["S1"], n) == brute(pp, ci, 1, "chain", n, partner=("radical_end", 6, None)), \
                (kid, n)
            checked += 1
    assert checked


def test_a_beta_scission_has_a_law(pp):
    got = terms(pp["table"], BETA)
    assert got and any(t["law"] is not None for _k, t in got), [(k, t.get("law_refused")) for k, t in got]


def test_a_head_to_tail_recombination_joins_into_a_chain(pp):
    got = [(kid, t) for kid, t in terms(pp["table"], RECOMBINATION)
           if kid == "recombination:radical=radical_end:CH+radical=radical_end:CH2"]
    assert got
    for _kid, t in got:
        assert t["law"] is not None, t.get("law_refused")
        joins = {tuple(p[:2]) for s in ("S0", "S1") for _m, prods in t["law"][s]["far"]["outcomes"] for p in prods}
        assert joins == {("join", "chain")}


# --------------------------------------------------------------------------
# nothing silent; readers of the old model refuse by name
# --------------------------------------------------------------------------
def test_every_term_has_a_law_or_says_why(pp):
    for k in pp["table"]["kernels"].values():
        for t in k["templates"]:
            assert t["law"] is not None or t.get("law_refused") or t.get("law_inert"), t


def test_a_class_no_law_reaches_is_listed_by_name(pp):
    unreached = pp["table"]["unreached"]
    participants = {c for k in pp["table"]["kernels"].values() for c in k["classes"]}
    lawful = {c for k in pp["table"]["kernels"].values() for t in k["templates"] if t["law"] is not None
              for c in k["classes"]}
    assert sorted(unreached["classes"]) == sorted(participants - lawful)


def test_the_mom_producer_refuses_a_table_of_the_new_state_model_by_name(pp):
    from rmgpy.kmc import mom_producer as mp
    label = "polypropyl"
    labels = dict(label=label, radical_label=label + mp.RADICAL_POOL_SUFFIX, alkene_label=label + mp.ALKENE_POOL_SUFFIX,
                  interior_radical_label=label + mp.INTERIOR_RADICAL_POOL_SUFFIX,
                  interior_alkene_label=label + mp.INTERIOR_ALKENE_POOL_SUFFIX, volume_ratio=1.0)
    with pytest.raises(mp.ProducerError, match=r"chain types .*c1~c1.*c2~c2"):
        mp.kernel_input(pp["table"], pp["spec"], ["homolysis:substrate=chain"], temperature=T, phi=1.0, **labels)


def test_the_kernel_reference_refuses_a_table_of_the_new_state_model_by_name(pp):
    from rmgpy.kmc import kernel_reference as ref
    with pytest.raises(ref.KernelReferenceError, match=r"chain types .*c1~c1.*c2~c2"):
        ref.KernelModel(pp["table"], pp["spec"])


def test_a_pe_table_keeps_the_old_state_model(pe):
    """One chain type, its own reverse: no key of the new state model, positions folded to the nearer end."""
    assert pe["schema"] == "kmc-kernel-rates/4"
    text = json.dumps(pe)
    for key in ('"chain_types"', '"unreached"', '"bulk_phases"', '"side"', '"head_distance"', '"tail_distance"',
                '"codes"'):
        assert key not in text, key
    assert kr_.state_model_beyond_v3(pe) is None
