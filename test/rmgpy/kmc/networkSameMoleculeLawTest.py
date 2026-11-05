"""The declared same-molecule law of a network (I-306: design I-297 law A
with the segment escape factor phi_seg), beyond the contract tests of
``networkSameMoleculeContractTest``.

Law A: an allowed pair of class ``c`` whose two sites lie on one molecule
fires at ``f w / (divisor N_A V)``, ``f`` = phi_seg = ``k_D(1) / (k_D(1) +
k_act)`` for a MET class (``k_D(1) = 4 pi N_A s 2 D0 sigma``, ``k_act`` the
molecule's own encounter group) and 1 otherwise.  The oracles are written
out from the declaration (``_k_d_seg``) and from the state (counts, strand
masses by networkx), never read back from the engine's rates.

Covered: the declaration (read, refused); phi_seg with competing products in
both directions (R2), exclusions, overrides and the reciprocal (phi-free)
ratio of the no-split reverse; R3 owner cases (loaded, excluded, deferred,
forbidden, zero weight) and the unimolecular owner of a no-split reverse;
the atom-resolved comparison; zero propensity of every refused or owned pair,
exactly; exact accounting of the leak denominator; sol molecules (R4);
rollback, split, merge, temperature change and full-rebuild equivalence; the
loop-share diagnostic (exact, sampled, trajectory-neutral)."""

import math

import numpy as np
import pytest

import melt_engine_support as S
import networkSameMoleculeContractTest as C

from rmgpy import constants
from rmgpy.kmc import melt_engine as me

D0 = C._NO_SPLIT_D_REF       # declared D0 (m2/s); Ea 0, so every temperature has it
W = C._NO_SPLIT_W            # the ring closure's weight


def _k_d_seg(transport):
    """The segment encounter (m3/mol/s) written out from the declared MET block."""
    return 4.0 * math.pi * constants.Na * transport.spin_factor * (2.0 * D0) * transport.contact_radius_m


def _ci(engine, name):
    return [t.name for t in engine.templates].index(name)


def _self_rate(engine, name):
    return dict(engine.propensities().channels).get(("self", _ci(engine, name)), 0.0)


def _table_of(engine, atom):
    return engine._self_of[engine.state.components.component_of(atom)]


def _reweighted(t, weight):
    from rmgpy.kmc.window_matcher import Template
    bonds = [(i, j, o) for (i, j), o in sorted(t.bond.items())]
    return Template(t.name, t.family, t.atoms, bonds, t.edits, weight=weight, overrides=t.overrides,
                    defers=t.defers)


def _engine(spec, state, templates, seed=3, excluded=(), **kw):
    return me.MeltEngine(spec, S.ToyCatalog(builder=lambda _t: list(templates)), 700.0,
                         np.random.default_rng(seed), state=state, excluded=excluded, **kw)


# --------------------------------------------------------------------------
# the declaration
# --------------------------------------------------------------------------
def test_the_law_is_read_from_the_network_declaration():
    """``transport_model["same_molecule"]`` declares law A with phi_seg; a
    network without it keeps its pairs as same_molecule leak and has no
    same-molecule channel."""
    spec, cut, m, r = C._opened_ring(D0)
    templates = [C._recombination(cut, m, r, "close", W)]
    engine = _engine(spec, cut, templates)
    assert me.SAME_MOLECULE_SCHEMA == C.SAME_MOLECULE_SCHEMA
    assert engine.same_law == me.SameMoleculeLaw("A", "seg")
    assert C.declared_law(engine) == ("A", "seg")
    bare = _engine(C._undeclared(spec), cut, templates)
    assert bare.same_law is None and C.declared_law(bare) == (None, None)
    assert C._classes(bare)["same_molecule"] > 0.0
    assert not any(k[0] == "self" for k, _a in bare.propensities().channels)
    assert "self" not in bare.snapshot()


@pytest.mark.parametrize("raw, match", [
    ({"schema": C.SAME_MOLECULE_SCHEMA, "law": "loop_only", "escape_factor": "segment"}, "not implemented"),
    ({"schema": C.SAME_MOLECULE_SCHEMA, "law": "well_mixed", "escape_factor": "self"}, "not implemented"),
    ({"law": "well_mixed", "escape_factor": "segment"}, "exactly"),
    ({"schema": "other/1", "law": "well_mixed", "escape_factor": "segment"}, "schema"),
    ("A", "exactly"),
])
def test_an_unimplemented_or_malformed_law_is_refused(raw, match):
    def mutate(mapping):
        mapping["transport_model"]["same_molecule"] = raw
    spec = C.law_spec(mutate, law=False)
    with pytest.raises(me.DeclarationError, match=match):
        me.MeltEngine(spec, S.ToyCatalog(builder=lambda _t: []), 700.0, np.random.default_rng(0))


def test_a_same_molecule_law_on_a_linear_declaration_is_refused():
    from rmgpy.polymer_spec import LinearHomopolymerSpec
    mapping = S.pe_mapping({10: 2})
    mapping["transport_model"]["same_molecule"] = {"schema": C.SAME_MOLECULE_SCHEMA, "law": "well_mixed",
                                                   "escape_factor": "segment"}
    with pytest.raises(me.DeclarationError, match="network declaration"):
        me.MeltEngine(LinearHomopolymerSpec.from_mapping(mapping), S.ToyCatalog(), 700.0,
                      np.random.default_rng(0))


# --------------------------------------------------------------------------
# phi_seg, competing products, R2
# --------------------------------------------------------------------------
def _square(key, excluded=()):
    states, templates, _counts, _cycle = C._h_square()
    spec, state = states[key]
    return _engine(spec, state, templates, excluded=excluded)


def test_competing_same_molecule_channels_share_the_molecules_escape_factor():
    """The opened ring of ``_h_square`` with X a radical (state 00): the
    closure (m + r, w) and the competing closure (m + X, 2w) are one encounter
    group, k_act = 3w; each fires at phi_seg x its own w / (N_A V); the
    non-MET H transfer (two H) at 2 w_h / (N_A V)."""
    engine = _square("00")
    k_d = _k_d_seg(engine.transport)
    phi = k_d / (k_d + 3.0 * W)
    assert 0.2 < phi < 0.9, "setup: a finite escape factor"
    assert _self_rate(engine, "close") == pytest.approx(phi * W / engine.nav, rel=1.0e-12)
    assert _self_rate(engine, "close_x") == pytest.approx(phi * 2.0 * W / engine.nav, rel=1.0e-12)
    assert _self_rate(engine, "h_fwd") == pytest.approx(2.0 * 1.0e-2 / engine.nav, rel=1.0e-12)
    assert engine.propensities().leak == 0.0


def test_the_no_split_reverse_reads_the_product_molecules_encounter_group():
    """R2 with competing products: the homolysis of the closed ring (state 10)
    whose product is state 00 is thinned at state 00's phi_seg, k_D / (k_D +
    3w), in the exact open-segment mean and in a fired proposal; the group
    built from the inverse alone (k_D / (k_D + w)) differs, so the assertion
    can fail.  Reciprocal: the forward rate over the physical reverse rate is
    phi-free, w / (N_A V) over the two equivalent homolyses at weight 1."""
    forward, reverse = _square("00"), _square("10")
    k_d = _k_d_seg(forward.transport)
    (table,) = forward._self_of.values()
    assert table.phi == pytest.approx(k_d / (k_d + 3.0 * W), rel=1.0e-12)
    open_ = _ci(reverse, "open")
    mean, sampled = reverse._fallback_phi(open_)
    assert not sampled and mean == pytest.approx(table.phi, rel=1.0e-15)
    assert abs(mean - k_d / (k_d + W)) > 0.05, "an inverse-only group would differ"
    physical = dict(reverse.propensities().channels)[("uni", open_)] * mean
    assert _self_rate(forward, "close") / physical == pytest.approx(W / forward.nav / 2.0, rel=1.0e-12)
    tr, view = reverse._draw_uni(open_)
    reverse._fire(tr, view, "uni", True)
    assert reverse.thinned_proposals == 1 and reverse.phi_sum == pytest.approx(table.phi, rel=1.0e-15)
    assert reverse.met_counts.get("reverse_no_pair", 0) == 0


@pytest.mark.parametrize("excluded, k_act", [(("close_x",), W), (("close",), 2.0 * W)])
def test_an_excluded_class_leaves_the_same_molecule_group_in_both_directions(excluded, k_act):
    """An excluded class is no row of the molecule's group, forward and
    reverse alike (``encounter_group``'s rule): with ``close_x`` excluded
    k_act = w, with ``close`` (the inverse's own class) excluded k_act = 2w;
    the excluded class has no channel and its pair is intrinsic excluded
    propensity."""
    forward, reverse = _square("00", excluded), _square("10", excluded)
    k_d = _k_d_seg(forward.transport)
    phi = k_d / (k_d + k_act)
    (table,) = forward._self_of.values()
    assert table.phi == pytest.approx(phi, rel=1.0e-12)
    (gone,) = excluded
    assert ("self", _ci(forward, gone)) not in dict(forward.propensities().channels)
    weight = {"close": W, "close_x": 2.0 * W}
    assert forward.propensities().excluded == pytest.approx(weight[gone] / forward.nav, rel=1.0e-12)
    mean, sampled = reverse._fallback_phi(_ci(reverse, "open"))
    assert not sampled and mean == pytest.approx(phi, rel=1.0e-12)


def test_with_every_met_class_excluded_a_no_split_reverse_has_no_owner():
    """No MET row left on the product: phi 1, counted ``reverse_no_pair``."""
    reverse = _square("10", ("close", "close_x"))
    open_ = _ci(reverse, "open")
    assert reverse._fallback_phi(open_) == (1.0, False)
    tr, view = reverse._draw_uni(open_)
    assert reverse._fire(tr, view, "uni", True) is not None
    assert reverse.met_counts["reverse_no_pair"] == 1 and reverse.phi_sum == 1.0


def _override_close(cut, m, r, weight):
    """``close`` with a depository override whose reactants are two copies of
    the opened ring and whose product joins m of one to r of the other: the
    opened ring is both participants, so its own (m, r) kind pair takes the
    override's weight (``_pair_override``)."""
    from meltEngineTest import _with_overrides
    from rmgpy.kmc.window_matcher import Override
    cid = cut.components.component_of(m)
    rows = me._first_rows(cut, cid)
    off = max(cut.atoms) + 1
    twin = {a + off: row[:5] + (tuple((b + off, o) for b, o in row[5]),) for a, row in rows.items()}
    reactants = {**rows, **twin}
    products = dict(reactants)
    a, b = products[m], products[r + off]
    products[m] = a[:3] + (0,) + a[4:5] + (tuple(sorted(a[5] + ((r + off, 1.0),))),)
    products[r + off] = b[:3] + (0,) + b[4:5] + (tuple(sorted(b[5] + ((m, 1.0),))),)
    return _with_overrides(C._recombination(cut, m, r, "close", W),
                           [Override("close_o", reactants, products, weight)])


def test_an_override_weight_enters_the_same_molecule_group_in_both_directions():
    """The opened ring's closure at the override's weight W_o: the row carries
    the override, the forward fires at phi_o W_o / (N_A V) with phi_o = k_D /
    (k_D + W_o), and the closed ring's homolysis is thinned at phi_o."""
    from rmgpy.kmc import melt_branched as mb
    w_o = 7.0 * W
    spec, cut, m, r = C._opened_ring(D0)
    close = _override_close(cut, m, r, w_o)
    forward = _engine(spec, cut, [close])
    k_d = _k_d_seg(forward.transport)
    phi = k_d / (k_d + w_o)
    assert _table_of(forward, m).group == ((0, 1.0, 0.0, 0),), "setup: the pair is the override's row"
    assert _self_rate(forward, "close") == pytest.approx(phi * w_o / forward.nav, rel=1.0e-12)
    state = mb.build_molecules(spec, C._RING)
    reverse = _engine(spec, state, [C._homolysis(state, m, r, "open", 1.0), close])
    mean, sampled = reverse._fallback_phi(_ci(reverse, "open"))
    assert not sampled and mean == pytest.approx(phi, rel=1.0e-12)


# --------------------------------------------------------------------------
# R3: a unimolecular closure owns the pair it duplicates, only when it fires
# --------------------------------------------------------------------------
def _ring_close(state, bridge, mine, name="ring_close", weight=1.0e3, defers=()):
    """The contract's R3 closure: a unimolecular template whose window covers
    the strand from the bridge radical to the aryl radical and forms their bond."""
    from collections import deque
    from rmgpy.kmc.window_matcher import AtomEdit, BondEdit, Template
    prev, queue = {bridge: None}, deque([bridge])
    while queue:
        x = queue.popleft()
        for y in sorted(state.atoms[x].edges):
            if y not in prev and state.atoms[y].element != "H":
                prev[y] = x
                queue.append(y)
    path, x = [], mine
    while x is not None:
        path.append(x)
        x = prev[x]
    atoms, bonds, index = C._closed_window(state, path, 0, radical=(bridge, mine))
    return Template(name, "toy_intra_recombination", atoms, bonds,
                    [BondEdit(index[bridge], index[mine], 1), AtomEdit(index[bridge], radical=0),
                     AtomEdit(index[mine], radical=0)], weight=weight, defers=defers)


def _owner_case():
    spec, state, atom = C._melt()
    return spec, state, atom[(0, 3)]["m"], atom[(0, 0)]["r4"]


def _pair_channel(engine):
    """phi_seg w / (N_A V) of the one same-molecule m_p pair (k_act = w)."""
    k_d = _k_d_seg(engine.transport)
    return k_d / (k_d + 1.0e-2) * 1.0e-2 / engine.nav


def test_a_loaded_owner_suppresses_the_pair_exactly():
    spec, state, bridge, mine = _owner_case()
    engine = _engine(spec, state, [_ring_close(state, bridge, mine),
                                   C._recombination(state, bridge, mine, "m_p", 1.0e-2)])
    table = _table_of(engine, bridge)
    assert table.owned == 1 and table.rows == () and table.group == ()
    assert _self_rate(engine, "m_p") == 0.0
    classes = C._classes(engine)
    assert classes["same_molecule"] == classes["short_cycle"] == classes["undeclared_link"] == 0.0


def test_an_excluded_owner_leaves_the_pair_its_chemistry():
    spec, state, bridge, mine = _owner_case()
    engine = _engine(spec, state, [_ring_close(state, bridge, mine),
                                   C._recombination(state, bridge, mine, "m_p", 1.0e-2)],
                     excluded=["ring_close"])
    assert _table_of(engine, bridge).owned == 0
    assert _self_rate(engine, "m_p") == pytest.approx(_pair_channel(engine), rel=1.0e-12)


def test_an_owner_whose_transition_is_deferred_away_owns_only_through_the_class_it_defers_to():
    """``ring_close`` defers to an identical ``ring_close_b``: its transition is
    ``ring_close_b``'s.  Loaded, ``ring_close_b`` owns the pair; excluded, no
    active transition makes the product and the pair is a channel."""
    spec, state, bridge, mine = _owner_case()
    templates = [_ring_close(state, bridge, mine, defers=("ring_close_b",)),
                 _ring_close(state, bridge, mine, "ring_close_b"),
                 C._recombination(state, bridge, mine, "m_p", 1.0e-2)]
    loaded = _engine(spec, state, templates)
    assert loaded.uni_allowed[_ci(loaded, "ring_close")] == 0, "setup: the deferring class leaves it"
    assert _table_of(loaded, bridge).owned == 1 and _self_rate(loaded, "m_p") == 0.0
    out = _engine(spec, state, templates, excluded=["ring_close_b"])
    assert _table_of(out, bridge).owned == 0
    assert _self_rate(out, "m_p") == pytest.approx(_pair_channel(out), rel=1.0e-12)


def _short_pair(state, atom, name):
    p2 = atom[(0, 2)]
    return C._wide_recombination(state, atom[(0, 3)]["m"], [], p2["r4"], [p2["r3"], p2["r5"], p2["r2"], p2["r6"]],
                                 name, 1.0e-2)


def test_a_forbidden_owner_owns_nothing_and_the_pair_stays_refused():
    """A unimolecular closure of the 4-bond cycle (bridge radical onto Ph2's
    para radical) is short_cycle leak; the same-molecule pair it would
    duplicate is not hidden: it is short_cycle leak of its own, at phi_seg w / (N_A V)."""
    spec, state, atom = C._melt(2)
    bridge, short = atom[(0, 3)]["m"], atom[(0, 2)]["r4"]
    engine = _engine(spec, state, [_ring_close(state, bridge, short), _short_pair(state, atom, "m_short")])
    ci = _ci(engine, "ring_close")
    assert engine.uni_allowed[ci] == 0 and engine.uni_leak[ci] == 1, "setup: the owner is forbidden"
    table = _table_of(engine, bridge)
    assert table.owned == 0 and [r.short_cycle for r in table.rows] == [1]
    assert _self_rate(engine, "m_short") == 0.0
    assert C._classes(engine)["short_cycle"] == pytest.approx(1.0e3 + _pair_channel(engine), rel=1.0e-12)


def test_a_zero_weight_owner_owns_nothing():
    """A template refuses a weight <= 0 (so no catalog loads one); the
    engine's own weight table is still checked: an owner at weight 0 does not
    suppress the pair."""
    from rmgpy.kmc.window_matcher import TemplateError
    spec, state, bridge, mine = _owner_case()
    with pytest.raises(TemplateError, match="positive"):
        _ring_close(state, bridge, mine, weight=0.0)
    engine = _engine(spec, state, [_ring_close(state, bridge, mine),
                                   C._recombination(state, bridge, mine, "m_p", 1.0e-2)])
    weights = list(engine.weights)
    weights[_ci(engine, "ring_close")] = 0.0
    engine.weights = tuple(weights)
    table = engine.self_table(engine.entries[engine.state.components.component_of(bridge)])
    assert table.owned == 0 and len(table.rows[0].allowed) == 1


def test_a_no_split_reverse_whose_inverse_a_unimolecular_owner_makes_is_accepted_at_phi_one():
    """R2/R3: the closed ring's homolysis with the closure loaded as a
    unimolecular template (built on the opened ring): the inverse is that
    owner's (atom-resolved), phi 1, counted ``reverse_unimolecular_owner``;
    with the owner excluded the MET pair is the group again (k_D / (k_D + w))."""
    from rmgpy.kmc import melt_branched as mb
    spec, cut, m, r = C._opened_ring(D0)
    owner = _ring_close(cut, m, r, "reclose", 1.0)
    close = C._recombination(cut, m, r, "close", W)
    state = mb.build_molecules(spec, C._RING)
    templates = [C._homolysis(state, m, r, "open", 1.0), owner, close]
    engine = _engine(spec, state, templates)
    open_ = _ci(engine, "open")
    assert engine._fallback_phi(open_) == (1.0, False)
    tr, view = engine._draw_uni(open_)
    assert engine._fire(tr, view, "uni", True) is not None
    assert engine.met_counts["reverse_unimolecular_owner"] == 1
    assert engine.met_counts.get("reverse_no_pair", 0) == 0
    out = _engine(spec, mb.build_molecules(spec, C._RING), templates, excluded=["reclose"])
    k_d = _k_d_seg(out.transport)
    mean, sampled = out._fallback_phi(open_)
    assert not sampled and mean == pytest.approx(k_d / (k_d + W), rel=1.0e-12)


# --------------------------------------------------------------------------
# the atom-resolved comparison
# --------------------------------------------------------------------------
def test_the_atom_resolved_comparison_tells_apart_isomorphic_products_on_other_atoms():
    """The ring's four methylene-ring homolyses make one product species
    (isomorphic graphs: ``_transition_is`` cannot tell them apart) on four
    atom sets: four distinct ``atom_resolved_product`` keys.  In each product
    the inverse key (``inverse_product``) restores the reactant's rows atom
    by atom, and it is the key of the one closure that re-forms that bond."""
    from rmgpy.kmc import melt_branched as mb
    spec, cut, m, r = C._opened_ring(D0)
    state = mb.build_molecules(spec, C._RING)
    engine = _engine(spec, state, [C._homolysis(state, m, r, "open", 1.0),
                                   C._recombination(cut, m, r, "close", W)])
    open_, close = _ci(engine, "open"), _ci(engine, "close")
    (entry,) = engine.entries.values()
    trs = [entry.uni_at(open_, k) for k in range(len(entry.shape.uni[open_][0]))]
    assert len(trs) == 4
    before = {a: me._node_row(n) for a, n in state.atoms.items()}
    keys, graphs, inverses = set(), [], set()
    for tr in trs:
        keys.add(me.atom_resolved_product(tr))
        tx, (product,) = engine._begin_apply(tr, me._LocalView(state, list(tr.components)))
        try:
            graphs.append(C._state_graph(state))
            inverse = me.inverse_product(tr, tx.journal.old_nodes, state.atoms)
            inverses.add(inverse)
            rows = {a: me._node_row(n) for a, n in state.atoms.items()}
            rows.update(dict(inverse))
            assert rows == before, "the inverse key restores every atom's row"
            sites = engine._scan([product])[product].sites
            (s0,) = [s for (ci, p, _k), ss in sites.items() if ci == close and p == 0 for s in ss]
            (s1,) = [s for (ci, p, _k), ss in sites.items() if ci == close and p == 1 for s in ss]
            closing, _view = engine._pair_transition(close, product, s0, product, s1)
            assert me.atom_resolved_product(closing) == inverse
        finally:
            tx.rollback()
    assert len(keys) == 4 and len(inverses) == 4
    assert all(C._isomorphic_states(graphs[0], g) for g in graphs[1:])


# --------------------------------------------------------------------------
# zero propensity, exactly, for every refused pair
# --------------------------------------------------------------------------
def _short_only():
    spec, state, atom = C._melt(2)
    return spec, state, atom[(0, 3)]["m"], [_short_pair(state, atom, "pair")]


def _undeclared_only():
    """The opened ring with an H off M3's methylene (both of its sites
    taken): its recombination with the aryl radical is an undeclared link."""
    from rmgpy.kmc import melt_branched as mb
    spec, cut, m, r = C._opened_ring(D0)
    m3 = mb.residue_atoms(cut)[(0, 3)]["m"]
    S.strip_h(cut, m3)
    return spec, cut, m3, [C._wide_recombination(cut, m3, [], r, [], "pair", 1.0e-2)]


@pytest.mark.parametrize("build, reason", [(_short_only, "short_cycle"), (_undeclared_only, "undeclared_link")])
def test_a_refused_same_molecule_pair_has_zero_propensity_exactly(build, reason):
    """The only pair of the box is refused by declared linearity: its channel
    is 0.0 exactly (no seed is sampled), the box is absorbing, and the pair is
    leak of its own class at phi_seg w / (N_A V)."""
    spec, state, radical, templates = build()
    engine = _engine(spec, state, templates)
    prop = engine.propensities()
    assert prop.enabled == 0.0 and _self_rate(engine, "pair") == 0.0
    assert all(row.allowed == () for row in _table_of(engine, radical).rows)
    classes = C._classes(engine)
    assert classes[reason] == pytest.approx(_pair_channel(engine), rel=1.0e-12)
    assert classes["same_molecule"] == 0.0
    result = engine.run(until=1.0e6)
    assert result.events == 0 and result.stop_reason == "absorbing"
    assert engine.leak_by_class[reason] == pytest.approx(classes[reason] * 1.0e6, rel=1.0e-12)


# --------------------------------------------------------------------------
# exact accounting, sol molecules
# --------------------------------------------------------------------------
def _gel(cuts=3):
    spec = C.law_spec()
    state, out = C._cut_fixture(spec, cuts)
    return spec, state, out


@pytest.mark.parametrize("fixture", ["square", "owner", "short", "gel"])
def test_every_measured_same_molecule_pair_is_a_channel_a_refusal_or_owned(fixture):
    """The table's enumeration of disjoint site pairs (``_self_pairs``) and the
    measured count (``_same_pairs``, ``engine.same``) are two copies of one
    rule: per molecule, allowed + short_cycle + undeclared_link + owned equals
    the measured count (no class excluded here)."""
    if fixture == "square":
        engine = _square("00")
    elif fixture == "owner":
        spec, state, bridge, mine = _owner_case()
        engine = _engine(spec, state, [_ring_close(state, bridge, mine),
                                       C._recombination(state, bridge, mine, "m_p", 1.0e-2)])
    elif fixture == "short":
        spec, state, _radical, templates = _short_only()
        engine = _engine(spec, state, templates)
    else:
        spec, state, templates = _gel_abstraction()
        engine = _engine(spec, state, templates)
    assert engine._self_of, "setup: same-molecule pairs"
    for cid, table in engine._self_of.items():
        measured = sum(engine._same_of.get(cid, {}).values())
        tabled = sum(len(r.allowed) + r.short_cycle + r.undeclared_link for r in table.rows) + table.owned
        assert measured == tabled > 0


def test_the_leak_denominator_integrates_the_exact_same_molecule_rates():
    """The same-molecule channels are exact (no thinning): the run's
    non-thinned integral and leak are the integrals of ``propensities()``
    read between events, and leak_den = leak_num + non_thinned + phi_sum."""
    spec, state, cuts = _gel()
    (m0, r0) = cuts[0]
    engine = _engine(spec, state, [C._recombination(state, m0, r0, "bridge_ring", 1.0e-2)], seed=1)
    first = engine.propensities()
    ci = _ci(engine, "bridge_ring")
    assert dict(first.channels)[("self", ci)] == math.fsum(t.rate(ci) for _c, t in sorted(engine._self_of.items()))
    marks = [(0.0, first)]

    def check(e, _event):
        marks.append((e.time, e.propensities()))
    engine.run(until=1.0e4, check=check)
    assert any(len(set(e.components)) == 1 for e in engine.events), "a same-molecule event"
    ends = [t for t, _p in marks[1:]] + [engine.time]
    non_thinned = math.fsum(p.non_thinned * (t1 - t0) for (t0, p), t1 in zip(marks, ends))
    leak = math.fsum(p.leak * (t1 - t0) for (t0, p), t1 in zip(marks, ends))
    assert non_thinned > 0.0
    assert engine.non_thinned == pytest.approx(non_thinned, rel=1.0e-9)
    assert engine.leak_num == pytest.approx(leak, rel=1.0e-9, abs=1.0e-300)
    assert engine.leak_den == engine.leak_num + engine.non_thinned + engine.phi_sum


def test_a_sol_molecule_fires_by_the_same_law():
    """R4: the far pair on a six-ring chain beside a twelve-ring chain (the
    gel by size) fires at w n / (divisor N_A V)."""
    from rmgpy.kmc import melt_branched as mb
    spec = C.law_spec()

    def names(rings):
        return ["Ph" if u % 2 == 0 else "M" for u in range(2 * rings - 1)]
    state = mb.build_molecules(spec, [(names(6), C._chain_links(6)), (names(12), C._chain_links(12))])
    atom = mb.residue_atoms(state)
    radical = atom[(0, 0)]["r4"]
    S.strip_h(state, radical)
    donor = atom[(0, 9)]["m"]
    far = [atom[(0, 10)][k] for k in ("r1", "r2", "r3", "r4", "r5", "r6", "o")]
    engine = _engine(spec, state, [C._abstraction(state, radical, [donor] + far, donor, "far_h", 1.0e-2)])
    comps = engine.state.components
    assert engine.gel_component() == comps.component_of(atom[(1, 0)]["r1"]) != comps.component_of(radical)
    ci = _ci(engine, "far_h")
    assert _self_rate(engine, "far_h") == pytest.approx(
        engine.weights[ci] * 2 / (engine.classes[ci].divisor * engine.nav), rel=1.0e-12)


# --------------------------------------------------------------------------
# rollback, split, merge, temperature, full rebuild
# --------------------------------------------------------------------------
def test_the_same_molecule_index_after_every_event_equals_a_full_recompute():
    """The ring Ph0-M1-Ph2-M3 beside a chain Ph-M-Ph, the methylene-ring
    homolysis and its recombination: the homolysis opens the ring (no split:
    thinned at phi_seg, often rejected and rolled back) or cuts the chain
    (a split); the recombination merges two molecules or closes one.  After
    every event the snapshot, same-molecule tables included, equals a full
    rebuild's (pk-i306 probe_equiv2: every seed 0-7 covers every kind)."""
    from rmgpy.kmc import melt_branched as mb
    spec, cut, m, r = C._opened_ring(D0)
    state = mb.build_molecules(spec, C._RING + [(["Ph", "M", "Ph"], C._chain_links(2))])
    engine = _engine(spec, state, [C._homolysis(state, m, r, "open", 1.0), C._recombination(cut, m, r, "close", W)],
                     seed=3)
    kinds = {"merge": 0, "split": 0, "same": 0, "no_split": 0}

    def check(e, event):
        assert e.snapshot() == e.fresh_snapshot()
        if event.template == "close":
            kinds["merge" if len(set(event.components)) == 2 else "same"] += 1
        else:
            kinds["split" if len(event.products) == 2 else "no_split"] += 1
    assert engine.snapshot() == engine.fresh_snapshot()
    engine.run(max_events=40, until=1.0e12, check=check)
    assert all(kinds.values()), kinds
    assert engine.met_counts["reverse_rejected"] >= 1, "rolled-back proposals"
    assert engine.met_counts["reverse_self_accepted"] == kinds["no_split"]
    assert engine.met_counts["self_forward"] == kinds["same"]


def _ring_donor(state):
    """An aryl radical on the ring's Ph0 para site and a donor window on M3
    (so the ring holds a same-molecule non-MET pair too)."""
    from rmgpy.kmc import melt_branched as mb
    atom = mb.residue_atoms(state)
    radical = atom[(0, 0)]["r4"]
    S.strip_h(state, radical)
    donor = atom[(0, 3)]["m"]
    return radical, [donor], donor


def test_a_rejected_no_split_proposal_leaves_the_index_and_the_state_unchanged():
    """Rollback (R5): a no-split proposal rejected at phi_seg leaves the
    state, the snapshot and the committed same-molecule tables (the same
    objects) as they were; a table is never read inside the open transaction."""
    from rmgpy.kmc import melt_branched as mb
    spec, cut, m, r = C._opened_ring(D0)
    state = mb.build_molecules(spec, C._RING)
    donor = _ring_donor(state)
    engine = _engine(spec, state, [C._homolysis(state, m, r, "open", 1.0), C._recombination(cut, m, r, "close", W),
                                   C._abstraction(state, *donor, "h", 1.0e-2)])
    assert engine._self_of, "setup: the ring holds a same-molecule pair"
    before, tables, digest = engine.snapshot(), dict(engine._self_of), state.state_hash()
    class Reject:
        """Every acceptance draw rejects (u just below 1)."""
        def __init__(self, rng):
            self.rng = rng

        def random(self):
            return 1.0 - 1.0e-12

        def __getattr__(self, name):
            return getattr(self.rng, name)
    engine.rng = Reject(engine.rng)
    tr, view = engine._draw_uni(_ci(engine, "open"))
    assert engine._fire(tr, view, "uni", True) is None
    assert engine.met_counts["reverse_rejected"] == 1 and 0.0 < engine.phi_sum < 1.0
    assert state.state_hash() == digest and engine.snapshot() == before
    assert set(engine._self_of) == set(tables) and all(engine._self_of[c] is t for c, t in tables.items())
    tx, _products = engine._begin_apply(tr, me._LocalView(state, list(tr.components)))
    try:
        with pytest.raises(me.EngineError, match="committed state"):
            engine.self_table(next(iter(engine.entries.values())))
    finally:
        tx.rollback()


def test_a_temperature_change_moves_phi_seg_with_the_weights():
    """The catalog's weight moves with temperature (w T / 700 K): after
    ``set_temperature`` the channel is phi_seg(w') w' / (N_A V) and the
    snapshot equals a full rebuild's at the new temperature."""
    spec, cut, m, r = C._opened_ring(D0)
    base = C._recombination(cut, m, r, "close", W)
    engine = me.MeltEngine(spec, S.ToyCatalog(builder=lambda t: [_reweighted(base, W * t / 700.0)]), 700.0,
                           np.random.default_rng(0), state=cut, temperature_range=(650.0, 850.0))
    engine.set_temperature(840.0)
    k_d = _k_d_seg(engine.transport)
    w = W * 840.0 / 700.0
    assert _self_rate(engine, "close") == pytest.approx(k_d / (k_d + w) * w / engine.nav, rel=1.0e-12)
    assert engine.snapshot() == engine.fresh_snapshot()


# --------------------------------------------------------------------------
# the loop-share diagnostic
# --------------------------------------------------------------------------
def _oracle_ratio(engine, s, t):
    """c_sp / c_bg of anchors ``s``, ``t``: networkx's weighted shortest path
    over heavy atoms, each bond the mean united-atom mass of its ends."""
    import networkx as nx
    from rmgpy.molecule.element import get_element
    state = engine.state
    cid = state.components.component_of(s)
    heavy = [a for a in state.components.atoms(cid) if state.atoms[a].element != "H"]
    mass = {a: 1000.0 * (get_element(state.atoms[a].element).mass + get_element("H").mass
                         * sum(state.atoms[b].element == "H" for b in state.atoms[a].edges)) for a in heavy}
    g = nx.Graph()
    g.add_weighted_edges_from((a, b, 0.5 * (mass[a] + mass[b])) for a in heavy for b in state.atoms[a].edges
                              if b in mass and a < b)
    r2 = engine.transport.r2_per_molar_mass * nx.dijkstra_path_length(g, s, t)
    return (3.0 / (2.0 * math.pi * r2)) ** 1.5 / constants.Na * engine.nav


def test_the_loop_share_is_exact_over_few_pairs():
    """The far pair (two ordered pairs, one per donor H, both anchored on the
    aryl radical and the donor carbon): exact, share = c_sp / c_bg from the
    oracle; the 24-ring share is smaller."""
    engine, _e = C._far_pair(6)
    s, t = engine._i297_anchors
    out = engine.loop_share(samples=8)
    assert out["exact"] and out["pairs"] == 2 and out["samples"] == 2 and out["unbounded"] == 0
    assert out["standard_error"] == 0.0
    assert out["share"] == pytest.approx(_oracle_ratio(engine, s, t), rel=1.0e-12)
    far, _e = C._far_pair(24)
    assert 0.0 < far.loop_share()["share"] < out["share"]
    assert "unquantified" in out["definition"]


def _gel_abstraction():
    spec, state, cuts = _gel(2)
    (m0, r0) = cuts[0]
    donor = next(b for b in sorted(state.atoms) if state.origin[b].key == "m" and not state.atoms[b].radical
                 and sum(state.atoms[x].element == "H" for x in state.atoms[b].edges) == 2)
    templates = [C._recombination(state, m0, r0, "bridge_ring", 1.0e-2),
                 C._abstraction(state, r0, [donor], donor, "h_gel", 1.0e-4)]
    return spec, state, templates


def test_the_sampled_loop_share_reports_its_size_and_error_and_changes_no_trajectory():
    """On the gel (many abstraction pairs) 16 samples: not exact, a positive
    standard error, reproducible from its seed.  A run observed with the
    diagnostic fires exactly the events of one without it (times, templates,
    components, the final state) and its rows differ only by ``loop_share``."""
    from rmgpy.kmc import melt_network as mn
    runs = []
    for loop_share in (None, 16):
        spec, state, templates = _gel_abstraction()
        engine = _engine(spec, state, templates, seed=2)
        rows = mn.observed_run(engine, [0.0, 1.0, 1.0e4], loop_share=loop_share)
        runs.append((engine, rows))
    (plain, rows0), (seen, rows1) = runs
    assert [(e.time, e.template, e.components, e.products) for e in plain.events] == \
        [(e.time, e.template, e.components, e.products) for e in seen.events]
    assert plain.events_fired >= 1 and plain.state.state_hash() == seen.state.state_hash()
    assert [{k: v for k, v in r.items() if k != "loop_share"} for r in rows1] == rows0
    first = rows1[0]["loop_share"]
    assert first["pairs"] > 16 and not first["exact"] and first["samples"] == 16
    assert first["standard_error"] > 0.0 and first["share"] > 0.0
    assert seen.loop_share(samples=16) == seen.loop_share(samples=16)
    assert seen.loop_share(samples=16, seed=1) != seen.loop_share(samples=16, seed=0)
