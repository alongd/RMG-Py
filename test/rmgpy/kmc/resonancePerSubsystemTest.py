"""I-298: the resonance bound applies to what is enumerated, not to the product
over a component's independent conjugated subsystems.

A component holding several independent delocalized radicals (benzylic,
allylic, phenoxy-type) has a form count equal to the product of its
subsystems' counts.  That product passed ``RESONANCE_FORM_LIMIT`` on cured
networks and stopped the run, although no window evaluates it: a window
touches only the subsystems holding its atoms.  The oracle of every test here
is the matcher as it was before this change (``legacy_component_subs`` and
``legacy_component_forms``, frozen copies of the product-bounded code) run
with the limit raised in the test only: an exact enumeration of the product.
"""

import pickle
import re
import sys
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from rmgpy.kmc import window_matcher as wm  # noqa: E402
from rmgpy.kmc.melt_state import MeltState, Node, Origin  # noqa: E402

EXACT = 10 ** 9   # the oracle's limit: no component here comes near it


# --------------------------------------------------------------------------
# the oracle: the matcher's product-bounded code before I-298, verbatim
# --------------------------------------------------------------------------
def legacy_component_subs(members, value, limit=wm.RESONANCE_FORM_LIMIT):
    subs, product = [], 1
    for group, free in wm._subsystems(members, value):
        bonds, shape = wm._sub_shape(group, free, value)
        budget = limit // product
        try:
            forms = wm._shape_forms(shape, budget)
        except wm.ConjugatedRingError:
            raise
        except wm.ResonanceLimitError as error:
            count = product * (budget + 1) if error.count is not None else None
            raise wm.ResonanceLimitError(
                f"component at atom {members[0]}: more than {limit} resonance forms "
                f"(stopped at form {count})", count) from None
        product *= len(forms)
        subs.append(wm._Sub(group, bonds, {a: k for k, a in enumerate(group)},
                            {b: k for k, b in enumerate(bonds)}, forms))
    return subs, product


def legacy_component_forms(state, cid, limit=wm.RESONANCE_FORM_LIMIT):
    members = state.components.atoms(cid)
    nodes = state.atoms

    def value(a):
        n = nodes[a]
        return n.element, n.charge, n.lone_pairs, n.radical, n.edges

    subs, _total = legacy_component_subs(members, value, limit)
    forms = [((), ())]
    for sub in subs:
        forms = [(bonds + tuple((a, b, o) for (a, b), o in zip(sub.bonds, orders)),
                  rads + tuple(zip(sub.atoms, r)))
                 for bonds, rads in forms for orders, r in sub.forms]
    return tuple(wm.Form(tuple(sorted(b)), tuple(sorted(r))) for b, r in forms)


@contextmanager
def exact_enumeration(limit=EXACT):
    """Run the matcher (and everything calling it) as it was before I-298,
    with ``limit`` in place of ``RESONANCE_FORM_LIMIT`` on the product."""
    saved = wm._component_subs, wm.component_forms
    wm._component_subs = lambda members, value, *_a, **_k: legacy_component_subs(members, value, limit)
    wm.component_forms = lambda state, cid: legacy_component_forms(state, cid, limit)
    wm._VIEWS.clear()
    try:
        yield
    finally:
        wm._component_subs, wm.component_forms = saved
        wm._VIEWS.clear()


def copy(state):
    return pickle.loads(pickle.dumps(state))


# --------------------------------------------------------------------------
# molecules: one saturated backbone carrying independent radical pendants
# --------------------------------------------------------------------------
_RING = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]


def _pendant(kind):
    """(heavy atoms (element, radical, charge, lone pairs), local bonds,
    hydrogens per heavy atom, the radical atom); local atom 0 bonds to the
    backbone."""
    if kind == "benzylic":                    # -CH(*)-C6H5: 7 atoms, 4 forms
        heavy = [("C", 1, 0, 0)] + [("C", 0, 0, 0)] * 6
        bonds = [(0, 1, 1)] + [(a + 1, b + 1, 1.5) for a, b in _RING]
        return heavy, bonds, [1, 0, 1, 1, 1, 1, 1], 0
    if kind == "allylic":                     # -CH(*)-CH=CH2: 3 atoms, 2 forms
        return [("C", 1, 0, 0), ("C", 0, 0, 0), ("C", 0, 0, 0)], [(0, 1, 1), (1, 2, 2)], [1, 1, 2], 0
    if kind == "phenoxy":                     # -C6H4-O(*) (para): 7 atoms, 4 forms
        heavy = [("C", 0, 0, 0)] * 6 + [("O", 1, 0, 2)]
        bonds = [(a, b, 1.5) for a, b in _RING] + [(3, 6, 1)]
        return heavy, bonds, [0, 1, 1, 0, 1, 1, 0], 6
    if kind.startswith("polyenyl"):           # -CH=CH-...-CH2(*): n carbons, (n + 1) / 2 forms
        n = int(kind[len("polyenyl"):])
        heavy = [("C", 1 if i == n - 1 else 0, 0, 0) for i in range(n)]
        bonds = [(i, i + 1, 2 if i % 2 == 0 and i + 1 < n - 1 else 1) for i in range(n - 1)]
        return heavy, bonds, [1] * (n - 1) + [2], n - 1
    raise ValueError(kind)


KINDS = ("benzylic", "allylic", "phenoxy")
FORMS = {"benzylic": 4, "allylic": 2, "phenoxy": 4}


def polyenyl_forms(kind):
    return (int(kind[len("polyenyl"):]) + 1) // 2


def build(kinds, seed=None, extra=True, shared=()):
    """A state: molecule 0 is a saturated backbone with one pendant of each
    of ``kinds`` (each on its own backbone carbon, two backbone carbons
    apart) and the pendants ``shared`` all on one more backbone carbon (the
    hub).  With ``extra``, molecule 1 is a methyl radical and molecule 2
    ethane.  ``seed`` (not None): atom ids are a random permutation.
    Returns (state, named atoms in the final ids): ``radicals`` (pendant
    radical atoms, in pendant order), ``attach`` (each pendant's atom bonded
    to the backbone), per kind its radical atoms, ``backbone`` (two bonded
    backbone carbons at the end, whose window touches no pendant), ``hub``, ``methyl``, ``ethane``."""
    heavy, bonds = [], []
    names = {"radicals": [], "attach": []}
    n_backbone = 2 * (len(kinds) + (1 if shared else 0)) + 3
    for i in range(n_backbone):
        heavy.append(("C", 0, 0, 0))
        if i:
            bonds.append((i - 1, i, 1))
    hydrogens = [3 if i in (0, n_backbone - 1) else 2 for i in range(n_backbone)]
    groups = [(2 * j + 1, [k]) for j, k in enumerate(kinds)]
    if shared:
        groups.append((2 * len(kinds) + 1, list(shared)))
        names["hub"] = 2 * len(kinds) + 1
    for anchor, group in groups:
        for kind in group:
            h, b, hy, rad = _pendant(kind)
            off = len(heavy)
            heavy += h
            bonds += [(i + off, j + off, o) for i, j, o in b] + [(anchor, off, 1)]
            hydrogens += hy
            hydrogens[anchor] -= 1
            names["radicals"].append(rad + off)
            names["attach"].append(off)
            names.setdefault(kind, []).append(rad + off)
    names["backbone"] = (n_backbone - 2, n_backbone - 1)    # no pendant within a bond of either
    mols = [(heavy, bonds, hydrogens)]
    if extra:
        mols += [([("C", 1, 0, 0)], [], [3]), ([("C", 0, 0, 0), ("C", 0, 0, 0)], [(0, 1, 1)], [3, 3])]
    rows, origin, local = [], [], []
    for mol, (h, b, hy) in enumerate(mols):
        base = len(rows)
        local.append(base)
        mine = [list(x) + [{}] for x in h]
        for i, j, o in b:
            mine[i][4][j], mine[j][4][i] = float(o), float(o)
        for i, n in enumerate(hy):
            for _ in range(n):
                mine.append(["H", 0, 0, 0, {i: 1.0}])
                mine[i][4][len(mine) - 1] = 1.0
        for i, row in enumerate(mine):
            rows.append((row[:4], {base + j: o for j, o in row[4].items()}))
            origin.append((mol, i))
    if extra:
        names["methyl"], names["ethane"] = local[1], local[2]
    perm = list(range(len(rows)))
    if seed is not None:
        perm = [int(x) for x in np.random.default_rng(seed).permutation(len(rows))]
    nodes, origins = {}, {}
    for old, ((element, radical, charge, lone_pairs), edges) in enumerate(rows):
        nodes[perm[old]] = Node(element, None, charge, radical, lone_pairs,
                                {perm[b]: o for b, o in edges.items()})
        mol, i = origin[old]
        origins[perm[old]] = Origin("mobile", mol, f"m{mol}", f"a{i}", False)

    def move(v):
        if isinstance(v, list):
            return [perm[a] for a in v]
        if isinstance(v, tuple):
            return tuple(perm[a] for a in v)
        return perm[v]
    return MeltState(nodes, origins), {k: move(v) for k, v in names.items()}


def templates(state, names):
    """Recombination of the first radical with the methyl radical, its
    abstraction of an ethane H (two participants; windows read from the
    state, so they touch the radical's subsystem), and the homolysis of a
    backbone C-C bond (one participant, touching no subsystem)."""
    import networkSameMoleculeContractTest as C
    touched = names["radicals"][0]
    a, b = names["backbone"]
    return [C._recombination(state, touched, names["methyl"], "rec", 1.0e-2),
            C._abstraction(state, touched, [names["ethane"]], names["ethane"], "abs", 1.0e-4),
            C._homolysis(state, a, b, "hom", 1.0e-3)]


def kinds_for(k):
    return tuple(KINDS[i % 3] for i in range(k))


def product_of(kinds):
    total = 1
    for kind in kinds:
        total *= FORMS[kind]
    return total


CASES = [(k, seed) for k in range(2, 7) for seed in (None, 11, 23)]


# --------------------------------------------------------------------------
# k independent radical subsystems: the exact enumeration's transitions
# --------------------------------------------------------------------------
@pytest.mark.parametrize("k,seed", CASES)
def test_independent_radical_subsystems_match_the_exact_enumeration(k, seed):
    kinds = kinds_for(k)
    expected = product_of(kinds)
    state, names = build(kinds, seed)
    cid = state.components.component_of(names["radicals"][0])
    with exact_enumeration():
        tps = templates(copy(state), names)
        oracle = wm.match(tps, copy(state))
        oracle_counted = wm.match(tps, copy(state), joint=False)
        assert len(legacy_component_forms(state, cid, EXACT)) == expected
    got = wm.match(tps, copy(state))
    assert got.transitions == oracle.transitions
    assert got.classes == oracle.classes
    assert wm.match(tps, copy(state), joint=False).classes == oracle_counted.classes
    # the windows evaluate the touched subsystem alone; the others multiply the count
    rec = [t for t in got.transitions if t.template == "rec"]
    hom = [t for t in got.transitions if t.template == "hom"]
    assert rec and hom
    assert all(t.embeddings % (expected // FORMS[kinds[0]]) == 0 for t in rec)
    assert all(t.embeddings % expected == 0 for t in hom)
    with exact_enumeration(wm.RESONANCE_FORM_LIMIT):
        if expected <= wm.RESONANCE_FORM_LIMIT:
            # where the product-bounded code completes, the outputs are its own
            today = wm.match(tps, copy(state))
            assert got.transitions == today.transitions and got.classes == today.classes
        else:
            with pytest.raises(wm.ResonanceLimitError):
                wm.match(tps, copy(state))


@pytest.mark.parametrize("k,seed", [(4, None), (6, 5)])
def test_loaded_transitions_match_the_exact_enumeration(k, seed):
    from rmgpy.kmc import template_emission as te
    state, names = build(kinds_for(k), seed)
    with exact_enumeration():
        tps = templates(copy(state), names)
        oracle = te.loaded_transitions(tps, copy(state))
    assert te.loaded_transitions(tps, copy(state)) == oracle


def _applied(tps, state):
    before = copy(state)
    out = []
    for tr in wm.match(tps, before).transitions:
        after = copy(before)
        wm.apply_transition(after, tr)
        out.append((tr.identity, after.full_hash(), sorted(wm.product_rows(before, tr).items())))
    return out


@pytest.mark.parametrize("k,seed", [(3, None), (5, 7), (6, 13)])
def test_applying_every_transition_gives_the_exact_enumerations_product(k, seed):
    state, names = build(kinds_for(k), seed)
    with exact_enumeration():
        tps = templates(copy(state), names)
        oracle = _applied(tps, state)
    got = _applied(tps, state)
    assert got == oracle and len(got) > 3


# --------------------------------------------------------------------------
# component_forms: factorized, exact count, lazy
# --------------------------------------------------------------------------
@pytest.mark.parametrize("k,seed", [(4, 3), (5, 9), (6, None)])
def test_component_forms_is_the_exact_product_factorized(k, seed):
    state, names = build(kinds_for(k), seed, extra=False)
    cid = state.components.component_of(names["radicals"][0])
    oracle = legacy_component_forms(state, cid, EXACT)
    forms = wm.component_forms(state, cid)
    assert len(forms) == forms.size == len(oracle) > wm.RESONANCE_FORM_LIMIT
    rng = np.random.default_rng(k)
    for i in [0, 1, len(oracle) - 1, -1, -len(oracle)] + [int(x) for x in rng.integers(0, len(oracle), 20)]:
        assert forms[i] == oracle[i]
    for i in (len(oracle), -len(oracle) - 1):
        with pytest.raises(IndexError):
            forms[i]
    # iteration enumerates in the oracle's order, and is bounded like any enumeration
    seen = []
    with pytest.raises(wm.ResonanceLimitError, match=r"more than 64 resonance forms"):
        for f in forms:
            seen.append(f)
    assert tuple(seen) == oracle[:wm.RESONANCE_FORM_LIMIT]
    assert next(f for f in forms if dict(f.radicals)[names["radicals"][0]] == 1) == oracle[0]


def test_component_forms_within_the_bound_is_unchanged():
    for kinds in (("benzylic", "phenoxy", "allylic"), ("benzylic", "phenoxy", "benzylic"), ("allylic",)):
        state, names = build(kinds, seed=4, extra=False)
        cid = state.components.component_of(names["radicals"][0])
        forms = wm.component_forms(state, cid)
        oracle = legacy_component_forms(state, cid)
        assert tuple(forms) == oracle and list(forms) == list(oracle)
        assert len(forms) == len(oracle) == product_of(kinds)
        assert [forms[i] for i in range(len(oracle))] == list(oracle)


def test_a_huge_product_is_counted_exactly_and_never_materialized():
    import networkSameMoleculeContractTest as C
    state, names = build(("benzylic",) * 40)
    cid = state.components.component_of(names["radicals"][0])
    forms = wm.component_forms(state, cid)
    assert forms.size == 4 ** 40 and isinstance(forms.size, int)
    assert wm._view_for(state).comp_forms[cid] == 4 ** 40
    assert all(dict(forms[0].radicals)[r] == 1 for r in names["radicals"])   # every radical on its CH
    a, b = names["backbone"]
    result = wm.match([C._homolysis(state, a, b, "hom", 1.0e-3)], state)
    assert result.transitions
    assert all(t.embeddings % 4 ** 40 == 0 for t in result.transitions)
    assert result.classes["hom"].embeddings == sum(t.embeddings for t in result.transitions)


# --------------------------------------------------------------------------
# what is enumerated stays bounded
# --------------------------------------------------------------------------
_MESSAGE = re.compile(r"^component at atom \d+: more than 64 resonance forms")


def test_a_single_subsystem_over_the_bound_still_raises():
    h_window = wm.Template("h_plus_h", "h", [wm.WindowAtom("H", radical=1), wm.WindowAtom("H", 1, radical=1)],
                           [], [wm.AtomEdit(0, radical=0), wm.AtomEdit(1, radical=0), wm.BondEdit(0, 1, 1.0)],
                           weight=1.0)
    for kinds in (("polyenyl129",), ("benzylic", "allylic", "polyenyl129"), ("polyenyl129", "phenoxy")):
        state, names = build(kinds)
        cid = state.components.component_of(names["radicals"][0])
        with pytest.raises(wm.ResonanceLimitError) as info:
            wm.component_forms(state, cid)
        assert type(info.value) is wm.ResonanceLimitError
        assert _MESSAGE.match(str(info.value)), str(info.value)
        assert info.value.count == wm.RESONANCE_FORM_LIMIT + 1
        with pytest.raises(wm.ResonanceLimitError) as info:
            wm.match([h_window], state)
        assert type(info.value) is wm.ResonanceLimitError and _MESSAGE.match(str(info.value))
    # one subsystem at the bound is fine, whatever the others hold
    state, names = build(("polyenyl127", "benzylic", "phenoxy"), extra=False)
    cid = state.components.component_of(names["radicals"][0])
    assert len(wm.component_forms(state, cid)) == 64 * 4 * 4
    wm.match([h_window], state)


def _window(state, rows, closed_heavy, participant):
    """``closed_heavy`` and their hydrogens closed, every other heavy
    neighbour a port, read from ``rows`` (one form of the component)."""
    from meltEngineTest import _ring_members
    ring = _ring_members(state)
    order = list(closed_heavy) + [y for x in closed_heavy for y, _o in rows[x][5] if rows[y][0] == "H"]
    ports = sorted({y for x in closed_heavy for y, _o in rows[x][5] if y not in order})
    order += ports
    index = {x: i for i, x in enumerate(order)}
    atoms = []
    for x in order:
        element, _iso, charge, radical, lone_pairs, edges = rows[x]
        if x in ports:
            exterior = tuple((o, rows[y][0], y in ring) for y, o in edges if y not in index)
            atoms.append(wm.WindowAtom(element, participant, "port", exterior, charge=charge, radical=radical,
                                       lone_pairs=lone_pairs))
        else:
            atoms.append(wm.WindowAtom(element, participant, charge=charge, radical=radical,
                                       lone_pairs=lone_pairs))
    bonds = sorted({(min(index[x], index[y]), max(index[x], index[y]), o)
                    for x in order for y, o in rows[x][5] if y in index})
    return atoms, bonds, index


@pytest.mark.parametrize("shared,over", [(("polyenyl17", "polyenyl19"), True),
                                         (("polyenyl5", "polyenyl7"), False)])
def test_a_window_over_subsystems_whose_joint_forms_pass_the_bound_raises(shared, over):
    """One participant whose window spans the first atoms of two polyenyl
    pendants on one backbone carbon, written in the form with both radicals
    at the far ends: its atoms allow every form but the one with the radical
    next to the hub (8 x 9 = 72 joint forms for 17 and 19 carbons, 2 x 3 for
    5 and 7)."""
    state, names = build(("allylic",), extra=False, shared=shared)
    hub, firsts, ends = names["hub"], names["attach"][1:], names["radicals"][1:]
    cid = state.components.component_of(hub)
    forms = wm.component_forms(state, cid)
    sizes = [polyenyl_forms(k) for k in shared]
    assert len(forms) == 2 * sizes[0] * sizes[1]
    # the form index with both far ends radical (and the allyl in its first form)
    index = 0
    for sub in forms.subsystems:
        k = 0
        for end in ends:
            if end in sub.apos:
                k = next(i for i, (_o, r) in enumerate(sub.forms) if r[sub.apos[end]] == 1)
        index = index * len(sub.forms) + k
    rows = wm.form_rows(state, cid, forms[index])
    assert all(rows[e][3] == 1 for e in ends)
    atoms, bonds, at = _window(state, rows, [hub] + firsts, 0)
    t = wm.Template("span", "toy", atoms, bonds,
                    [wm.BondEdit(at[hub], at[firsts[0]], 0), wm.AtomEdit(at[hub], radical=1),
                     wm.AtomEdit(at[firsts[0]], radical=1)], weight=1.0)
    if over:
        with pytest.raises(wm.ResonanceLimitError, match=r"^component at atom \d+: more than 64 resonance "
                                                         r"forms over the 2 subsystems") as info:
            wm.match([t], state)
        assert type(info.value) is wm.ResonanceLimitError and info.value.count == (sizes[0] - 1) * (sizes[1] - 1)
    else:
        with exact_enumeration():
            oracle = wm.match([t], copy(state))
        got = wm.match([t], copy(state))
        assert got.transitions == oracle.transitions and got.classes == oracle.classes
        assert got.transitions and got.classes["span"].evaluated % ((sizes[0] - 1) * (sizes[1] - 1)) == 0


def test_the_provenance_report_keeps_the_component_cutoff():
    """``windows.diagnostics`` reports a component's whole form count, or its
    overflow cutoff: unchanged, byte for byte."""
    from rmgpy.kmc import windows

    class Witness:
        def __init__(self, rows):
            self.rows = rows

        def reactant_rows(self):
            return self.rows

        def product_rows(self):
            return self.rows

    for kinds, expected in ((("benzylic", "phenoxy"), {"forms": 16}),
                            (("benzylic", "phenoxy", "allylic"), {"forms": 32}),
                            (("benzylic", "phenoxy", "benzylic"), {"forms": 64}),
                            (("benzylic", "phenoxy", "benzylic", "allylic"), {"cutoff": "> 64 forms"}),
                            (("benzylic",) * 6, {"cutoff": "> 64 forms"})):
        state, names = build(kinds, extra=False)
        cid = state.components.component_of(names["radicals"][0])
        rows = wm.form_rows(state, cid, legacy_component_forms(state, cid, EXACT)[0])
        out = windows.diagnostics(Witness(rows), [names["radicals"][0]])
        got = [{k: v for k, v in r.items() if k in ("forms", "cutoff")} for r in out["resonance"]]
        assert got == [expected, expected]


# --------------------------------------------------------------------------
# the reduced novolac network (law A + phi_seg) that stopped on the product
# --------------------------------------------------------------------------
def _reduced(ncut, homolysis, abstraction, seed, max_events=60):
    """pk-i306 ``gel_run.py reduced law``: the network deck with ``ncut``
    bridges cut, one recombination per cut, the homolysis of a pristine
    bridge and an aryl radical taking a bridge H; returns the event log
    with the exact total propensity after every event, and the report."""
    import melt_engine_support as S
    import networkSameMoleculeContractTest as C
    from meltEngineTest import _recombination
    from network_support import network_spec, with_met
    from rmgpy.kmc import melt_branched as mb
    from rmgpy.kmc import melt_engine as me
    from rmgpy.kmc.melt_state import SetBond, SetRadical, materialize

    def mutate(mapping):
        with_met(mapping)
        C.declare_law(mapping)
    spec = network_spec(mutate)
    state = materialize(spec)
    atom = mb.residue_atoms(state)
    cuts = []
    for (_mol, _unit), keys in sorted(atom.items()):
        if "m" in keys and len(cuts) < ncut:
            m = keys["m"]
            ring = next(b for b in sorted(state.atoms[m].edges) if state.atoms[b].element == "C")
            state.apply([SetBond(m, ring, 1.0, 0.0), SetRadical(m, 0, 1), SetRadical(ring, 0, 1)])
            cuts.append((m, ring))
    tps, seen = [], set()
    for k, (m, r) in enumerate(cuts):
        t = _recombination(state, m, r, f"bridge_ring_{k}")
        key = (t.atoms, tuple(sorted(t.bond.items())), t.edits)
        if key not in seen:
            seen.add(key)
            tps.append(t)
    pristine = next(b for b in sorted(state.atoms) if state.origin[b].key == "m" and not state.atoms[b].radical
                    and sum(state.atoms[x].element == "H" for x in state.atoms[b].edges) == 2)
    ring = next(x for x in sorted(state.atoms[pristine].edges) if state.atoms[x].element == "C")
    if homolysis:
        tps.append(C._homolysis(state, pristine, ring, "bridge_open", 1.0e-3))
    if abstraction:
        tps.append(C._abstraction(state, cuts[0][1], [pristine], pristine, "aryl_takes_bridge_h", 1.0e-4))
    engine = me.MeltEngine(spec, S.ToyCatalog(builder=lambda _t: list(tps)), 700.0,
                           np.random.default_rng(seed), state=state)
    log, stop = [("t0", C.exact_total(engine))], "horizon"
    try:
        while engine.events_fired < max_events and engine.time < 1.0e4:
            event = engine.step(until=1.0e4)
            if event is None:
                break
            log.append((event.time, event.template, event.channel, event.components, event.products,
                        C.exact_total(engine), engine.leak_num, engine.leak_den))
    except wm.ResonanceLimitError as error:
        stop = f"{type(error).__name__}: {error}"
    rep = engine.report()
    return stop, log, {k: rep[k] for k in ("events", "null_proposals", "leak_num", "leak_den", "leak_fraction",
                                           "leak_by_class", "template_counts", "met")}


@pytest.mark.parametrize("ncut,homolysis,abstraction,seed", [
    (2, 0, 1, 1), (2, 1, 1, 1), (3, 1, 0, 2), (4, 1, 0, 1), (6, 0, 0, 1), (8, 0, 0, 3)])
def test_the_reduced_network_runs_as_the_exact_enumeration(ncut, homolysis, abstraction, seed):
    with exact_enumeration():
        oracle = _reduced(ncut, homolysis, abstraction, seed)
    got = _reduced(ncut, homolysis, abstraction, seed)
    assert got == oracle
    assert len(got[1]) > 1
    if (ncut, homolysis, abstraction, seed) == (2, 1, 1, 1):
        # a candidate product holds one 16-atom two-ring subsystem of 113 forms (four radicals, no aryl
        # sigma-radical: one CH2* on ring A, two on ring B, one on the bridge carbon joining them): a
        # single subsystem over the bound, refused by both, after the same six events.  The stop was at
        # five events, on a 14-atom piece of 65 forms that only exceeded the bound while its aryl
        # sigma-radical shifted into the ring; held on its atom, that piece fits and the run goes one
        # event further
        assert got[0] == "ResonanceLimitError: more than 64 resonance forms" and len(got[1]) == 7
    else:
        assert got[0] == "horizon"
