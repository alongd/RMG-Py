"""Tests for rmgpy/kmc/windows.py: recertified closed windows (items 1-3),
the span bound (4), chain-end contexts (5), the short-chain ladder (6) and the
spin/resonance diagnostics (7).

Fast part: stub observers on toy witnesses (no RMG database) for the fixed
point; the provider's end reconstruction; context enumeration and merging;
the span bound, L and the ladder plan; support and closure on real
certificates (the RMG database is loaded once).
Slow part (``RMG_KMC_SLOW=1``): the acceptance tables (``windows_acceptance``).
"""

import math
import os
from pathlib import Path

import pytest

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import domain as dom
from rmgpy.kmc import rmg_moves as rm
from rmgpy.kmc import windows as win

SLOW = os.environ.get("RMG_KMC_SLOW") == "1"
HEAD, TAIL = (0, "head"), (0, "tail")


# --------------------------------------------------------------------------
# toy witnesses (no database): an alkane with one edited atom
# --------------------------------------------------------------------------
def toy(carbons=8, centre=3):
    """C<carbons> with both ends artificial; the product side puts a radical
    on ``centre`` (an edit; rows are not checked chemically)."""
    build = b.PROVIDERS["linear_alkane"].build(carbons)
    _mols, parts = b.witness_reactants([build])
    rows = rm.molecule_rows(_mols[0])
    product = dict(rows)
    row = product[centre]
    product[centre] = row[:3] + (1,) + row[4:]
    return b.make_witness(tuple(sorted(rows.items())), tuple(sorted(product.items())), parts,
                          label="toy")


def ball(witness, centre, radius):
    rows = witness.reactant_rows()
    seen = {centre}
    frontier = [centre]
    for _ in range(radius):
        nxt = []
        for a in frontier:
            for n, _o in rows[a][5]:
                if n not in seen:
                    seen.add(n)
                    nxt.append(n)
        frontier = nxt
    return frozenset(seen)


def stub(view_of, radius=lambda w: 1, centre=3, log=None):
    """An observer whose support is the ball of ``radius(w)`` bonds around
    ``centre`` and whose certificate view is ``view_of(w)``."""
    def observe(witness):
        if log is not None:
            log.append(dict(witness.extension_counts))
        support = win.Support(atoms=ball(witness, centre, radius(witness)),
                              reacting=frozenset({centre}))
        value = view_of(witness)
        return win.Observation(support=support, viewer=lambda inside, v=value: v)
    return observe


# --------------------------------------------------------------------------
# item 3: recertified fixed point
# --------------------------------------------------------------------------
def test_fixed_point_pins():
    assert win.LOOP_STATUSES == ("fixed_point", "no_fixed_point", "ceiling_exceeded",
                                 "certificate_unstable")
    assert isinstance(win.PORT_DEPTH_LIMIT, int) and win.PORT_DEPTH_LIMIT >= 2
    assert isinstance(win.FIXED_POINT_RETRIES, int) and win.FIXED_POINT_RETRIES >= 1
    res = win.fixed_point(toy(), stub(lambda w: "same"))
    assert res.status == "fixed_point", res.reason
    # ports in vector order (participant, port id); no extension was needed
    assert [k for k, _d in res.depths] == list(b.vector_ports(toy()))
    assert dict(res.depths) == {HEAD: 0, TAIL: 0}
    # one independent +1 check per artificial port, both passed
    assert sorted(c["port"] for c in res.checks[-1]) == sorted([HEAD, TAIL])
    assert all(c["same_window"] and c["same_certificate"] for c in res.checks[-1])
    assert res.window.atoms >= ball(toy(), 3, 1)


def test_fixed_point_checks_every_port_independently():
    """The certificate changes when the TAIL is extended to depth 1 only: a
    single +1 check (or a joint extension of every port) accepts too early;
    the per-port independent check extends the tail to depth 2."""
    res = win.fixed_point(toy(), stub(lambda w: "s1" if w.depth(TAIL) == 1 else "s0"))
    assert res.status == "fixed_point", res.reason
    assert dict(res.depths) == {HEAD: 0, TAIL: 2}
    failed = [c for round_ in res.checks for c in round_ if not c["same_certificate"]]
    assert {c["port"] for c in failed} == {TAIL}


def test_changing_certificate_is_not_a_fixed_point_until_it_stabilizes():
    seen = []
    res = win.fixed_point(toy(), stub(lambda w: min(w.depth(HEAD) + w.depth(TAIL), 4), log=seen))
    assert res.status == "fixed_point", res.reason
    assert dict(res.depths) == {HEAD: 2, TAIL: 2}
    # every extension was re-observed (recertified), the check witnesses too
    assert {HEAD: 2, TAIL: 2} in seen and {HEAD: 3, TAIL: 2} in seen
    assert res.retries == 2


def test_never_stabilizing_certificate_is_certificate_unstable():
    res = win.fixed_point(toy(), stub(lambda w: w.depth(HEAD) + w.depth(TAIL)))
    assert res.status == "certificate_unstable", (res.status, res.reason)
    assert res.retries == win.FIXED_POINT_RETRIES + 1 or res.retries == win.FIXED_POINT_RETRIES


def test_ceiling_refuses_instead_of_truncating():
    """The support reaches the head port: it must be extended; the ceiling
    forbids it, so the witness is refused, never truncated."""
    w = toy(centre=1)
    res = win.fixed_point(w, stub(lambda w_: "s", radius=lambda w_: 2, centre=1),
                          ceiling=w.heavy_count())
    assert res.status == "ceiling_exceeded", (res.status, res.reason)
    assert res.window is None
    # with room it reaches a fixed point after extending the head
    ok = win.fixed_point(w, stub(lambda w_: "s", radius=lambda w_: 2, centre=1))
    assert ok.status == "fixed_point" and dict(ok.depths)[HEAD] >= 1


def test_growing_window_is_no_fixed_point():
    res = win.fixed_point(toy(centre=1),
                          stub(lambda w: "s", radius=lambda w: 2 + w.depth(HEAD), centre=1))
    assert res.status == "no_fixed_point", (res.status, res.reason)
    assert "limit" in res.reason


# --------------------------------------------------------------------------
# item 5: physical contexts (provider reconstruction), enumeration, merging
# --------------------------------------------------------------------------
def ps_witness(units=3, **kw):
    build = b.PROVIDERS["polystyrene"].build("ps", b.PsChainSpec(units, **kw))
    _mols, parts = b.witness_reactants([build])
    rows = tuple(sorted(rm.molecule_rows(_mols[0]).items()))
    return b.make_witness(rows, rows, parts, label="ps")


def ps_pair():
    builds = [b.PROVIDERS["polystyrene"].build("pe", b.PsChainSpec(1, tail_end="benzylic_radical")),
              b.PROVIDERS["polystyrene"].build("dp3", b.PsChainSpec(3))]
    _mols, parts = b.witness_reactants(builds)
    rows = tuple(sorted(rm.molecule_rows(_mols[0]).items()) + sorted(
        rm.molecule_rows(_mols[1]).items()))
    return b.make_witness(rows, rows, parts, label="pair")


def _partner(witness, key):
    """The backbone neighbour of a port atom in its own unit."""
    meta = witness.participant(key[0]).metadata
    port = witness.port(key)
    unit, position = meta.coordinate(port.atom_id)
    other = 1 - position
    return next(a for a, u, k in meta.coordinates if (u, k) == (unit, other))


def _h_count(rows, atom):
    return sum(1 for n, _o in rows[atom][5] if rows[n][0] == "H")


@pytest.mark.parametrize("key,end,expect", [
    (HEAD, "methyl", "methyl"), (HEAD, "primary_radical", "primary_radical"),
    (HEAD, "vinylidene", "vinylidene"), (TAIL, "benzyl", "benzyl"),
    (TAIL, "benzylic_radical", "benzylic_radical"), (TAIL, "styryl", "styryl")])
def test_physical_end_is_reconstructed_not_relabelled(key, end, expect):
    w = ps_witness()
    port = w.port(key)
    (cap,) = port.cap_atom_ids
    new, record = b.reconstruct_end(w, key, end)
    for rows in (new.reactant_rows(), new.product_rows()):
        assert cap not in rows
    built = new.port(key)
    assert built.physicality == "physical" and built.cap_atom_ids == ()
    assert built.atom_id == port.atom_id
    rows = new.reactant_rows()
    partner = _partner(w, key)
    order = dict(rows[port.atom_id][5]).get(partner)
    if end in ("vinylidene", "styryl"):
        assert order == 2.0
        assert _h_count(rows, partner) == _h_count(w.reactant_rows(), partner) - 1
        assert record.deleted and set(record.deleted) - {cap}
    else:
        assert order == 1.0
    if end in ("methyl", "benzyl"):
        fresh = [a for a in record.fresh_ids]
        assert fresh and all(a not in w.reactant_rows() for a in fresh)
        assert any(n in fresh for n, _o in rows[port.atom_id][5])
    if end.endswith("radical"):
        assert rows[port.atom_id][3] == 1
    verdict = dom.check_rows(dict(rows))
    assert verdict.accepted, verdict.reason
    assert expect in verdict.details["ends"]


def test_end_reconstruction_never_touches_a_reacting_atom():
    w = toy(centre=0)  # the head port atom itself is edited
    with pytest.raises(b.EndConflict):
        b.reconstruct_end(w, HEAD, "primary_radical")
    # a methyl end only replaces the (unreacting) cap by a fresh end H
    built, record = b.reconstruct_end(w, HEAD, "methyl")
    assert b.edits(built) == b.edits(w) and record.fresh_ids
    # a reacting cap is never deleted
    capped = toy(centre=2)
    (cap,) = capped.port(HEAD).cap_atom_ids
    product = dict(capped.product_rows())
    product[cap] = product[cap][:3] + (1,) + product[cap][4:]
    edited = b.make_witness(capped.reactants, tuple(sorted(product.items())),
                            capped.participants)
    with pytest.raises(b.EndConflict):
        b.reconstruct_end(edited, HEAD, "methyl")


def test_context_assignments_policy():
    w = ps_witness()
    alone = win.context_assignments(w, near_ports=())
    assert alone[0] == ((HEAD, None), (TAIL, None))
    singles = [a for a in alone if sum(e is not None for _k, e in a) == 1]
    # far ports: their radical ends only (read whole-species)
    assert sorted(e for a in singles for _k, e in a if e) == ["benzylic_radical",
                                                             "primary_radical"]
    # radical end chemistries interact through the radical inventory: their
    # combinations are always enumerated
    assert ((HEAD, "primary_radical"), (TAIL, "benzylic_radical")) in alone
    assert len(alone) == 4
    one = win.context_assignments(w, near_ports=(HEAD,))
    assert len([a for a in one if sum(e is not None for _k, e in a) == 1]) == 4
    near = win.context_assignments(w, near_ports=(HEAD, TAIL))
    assert len(near) == 16 and len(set(near)) == 16
    assert win.context_assignments(w) == near
    # two participants: per-participant products, cross products only on request
    pair = ps_pair()
    per = win.context_assignments(pair, cross=False)
    crossed = win.context_assignments(pair, cross=True)
    # base, 9 singles, 4 radical combinations, the 3x3 product of the DP3 ends
    # less its radical pair (already a radical combination)
    assert len(per) == 1 + 9 + 4 + (9 - 1) and len(crossed) == 4 * 4 * 4
    mixed = [a for a in per if sum(e is not None for _k, e in a) > 1
             and len({k[0] for k, e in a if e is not None}) > 1]
    assert mixed and all(e in b.RADICAL_ENDS for a in mixed for _k, e in a if e)
    physical_head = ps_witness(radical_end=True)
    only_tail = win.context_assignments(physical_head, near_ports=(TAIL,))
    assert {k for a in only_tail for k, _e in a} == {TAIL}


def _ctx(witness, atoms, view):
    window = win.Window(atoms=frozenset(atoms), support=win.Support(atoms=frozenset(atoms)))
    obs = win.Observation(support=window.support, viewer=lambda inside, v=view: v)
    return win.ContextResult(assignment=(), status="fixed_point", witness=witness, window=window,
                             observation=obs)


def test_contexts_with_equal_windows_but_unequal_certificates_do_not_merge():
    w = toy()
    atoms = ball(w, 3, 1)
    same = win.merge_contexts([_ctx(w, atoms, "a"), _ctx(w, atoms, "a")])
    assert same == ((0, 1),)
    different = win.merge_contexts([_ctx(w, atoms, "a"), _ctx(w, atoms, "b")])
    assert different == ()
    # unequal windows never merge, even with equal certificates
    assert win.merge_contexts([_ctx(w, atoms, "a"), _ctx(w, ball(w, 3, 2), "a")]) == ()
    # refusals never merge
    refused = win.ContextResult(assignment=(), status="refused:domain")
    assert win.merge_contexts([refused, refused]) == ()


# --------------------------------------------------------------------------
# items 4 and 6: span bound, L, ladder
# --------------------------------------------------------------------------
def test_span_bound_formula():
    c = win.FamilyConstants("X", radius=4, reacting_distance=1, free_label_search=True)
    bound = win.span_bound(c)
    assert win.CLOSURE_EXTENSION == 1
    assert bound.reach == 5   # predicate radius + one closure bond
    # a free-label search adds no extent: its candidates are reacting atoms
    # (otherwise the member is refused)
    assert bound.backbone_atoms == 1 + 2 * 5 + 1 == 12
    assert bound.repeats == math.ceil((bound.backbone_atoms + 1) / 2) == 7
    assert bound.L == (bound.backbone_atoms + 1) // 2 == 6
    anchored = win.span_bound(win.FamilyConstants("Y", radius=2, reacting_distance=0,
                                                  free_label_search=False))
    assert anchored.backbone_atoms == 0 + 2 * 3 + 1 and anchored.L == 4
    assert "reach" in bound.formula and "L" in bound.formula


@pytest.mark.parametrize("radius,distance,free", [(2, 1, False), (4, 0, False), (4, 1, True),
                                                  (5, 2, False)])
def test_L_and_L_plus_one(radius, distance, free):
    bound = win.span_bound(win.FamilyConstants("Z", radius, distance, free))
    L, A = bound.L, bound.backbone_atoms
    phases = [dict(head_shift=h, tail_shift=t) for h in (False, True) for t in (False, True)]
    fits = [win.chain_backbone_atoms(b.PsChainSpec(L, **p)) for p in phases]
    assert min(fits) <= A  # some phase of DP L fits in one window component
    over = [win.chain_backbone_atoms(b.PsChainSpec(L + 1, **p)) for p in phases]
    assert min(over) > A   # no phase of DP L+1 does


def test_ladder_covers_cut_phases_and_controls():
    c = win.FamilyConstants("R_Recombination", radius=2, reacting_distance=1,
                            free_label_search=False, bimolecular=True)
    bound = win.span_bound(c)
    plan = win.ladder_plan(c, bound)
    items = plan.items
    assert max(i.dp for i in items) == bound.L + 1
    assert all(i.control for i in items if i.dp == bound.L + 1)
    assert not any(i.control for i in items if i.dp <= bound.L)
    kinds = {k for i in items for k in i.kinds}
    for kind in ("intact", "scission:phase_A", "scission:phase_B", "pendant_loss",
                 "ch_cleavage:benzylic", "ch_cleavage:secondary", "pair:short_short",
                 "pair:short_long"):
        assert kind in kinds, kind
    specs = [s for i in items for s in i.specs]
    assert any(s.head_shift for s in specs) and any(s.tail_shift for s in specs)
    # every declared end chemistry at every DP 1..L+1 of the intact chains
    intact = {(i.dp, i.specs[0].head_end, i.specs[0].tail_end) for i in items
              if "intact" in i.kinds}
    for dp in range(1, bound.L + 2):
        for head in ("methyl", "primary_radical", "vinylidene"):
            for tail in ("benzyl", "benzylic_radical", "styryl"):
                assert (dp, head, tail) in intact or (dp, head, tail) in plan.excluded_keys
    # every planned species is in the declared domain (the rest is excluded
    # with its domain reason, never silently dropped)
    for item in items:
        for spec in item.specs:
            build = b.PROVIDERS["polystyrene"].build("x", spec)
            assert dom.check(build.molecule).accepted, (item.label, dom.check(build.molecule).reason)
    assert all(reason for _key, reason in plan.excluded)


# --------------------------------------------------------------------------
# item 2: closure (no database)
# --------------------------------------------------------------------------
def test_closure_runs_on_both_reaction_sides():
    """The product side makes a 3=4 double bond: the pi system {3, 4} exists
    on the product side only, so atom 4 (and its H) is in the window only if
    the products are closed too."""
    w = toy()
    product = dict(w.product_rows())
    for a, n in ((3, 4), (4, 3)):
        row = product[a]
        edges = tuple(sorted(((x, 2.0 if x == n else o) for x, o in row[5])))
        product[a] = row[:3] + (0,) + row[4:5] + (edges,)
    w2 = b.make_witness(w.reactants, tuple(sorted(product.items())), w.participants)
    window = win.close(w2, win.Support(atoms=frozenset({3}), reacting=frozenset({3, 4})))
    assert window.refusal is None
    assert 4 in window.atoms
    h4 = {n for n, _o in w2.reactant_rows()[4][5] if w2.reactant_rows()[n][0] == "H"}
    assert h4 <= window.atoms
    reactant_only = win.close(w2, win.Support(atoms=frozenset({3}), reacting=frozenset({3})),
                              sides=("reactants",))
    assert 4 not in reactant_only.atoms


def ring_witness(size, centre=0):
    from rmgpy.molecule import Molecule
    m = Molecule(smiles="C1" + "C" * (size - 1) + "1")
    for k, atom in enumerate(m.atoms):
        atom.id = k
    rows = rm.molecule_rows(m)
    product = dict(rows)
    product[centre] = rows[centre][:3] + (1,) + rows[centre][4:]
    return b.make_witness(tuple(sorted(rows.items())), tuple(sorted(product.items())), ())


def test_ring_closure_and_refusals():
    eight = win.close(ring_witness(8), win.Support(atoms=frozenset({0}),
                                                   reacting=frozenset({0})))
    assert eight.refusal is None and set(range(8)) <= eight.atoms
    nine = win.close(ring_witness(9), win.Support(atoms=frozenset({0}),
                                                  reacting=frozenset({0})))
    assert nine.refusal and nine.refusal.startswith("refused:ring")
    assert nine.atoms == frozenset({0})       # refused, never truncated to a partial ring
    touching = win.close(ring_witness(9), win.Support(atoms=frozenset({0, 1}),
                                                      reacting=frozenset({20})))
    assert touching.refusal and touching.refusal.startswith("refused:large_ring")


def test_fixed_point_reports_observation_refusals():
    """A domain (or any certificate) refusal stops the loop by name, without
    a window: never a truncated one."""
    def observe(witness):
        return win.Observation(win.Support(frozenset()), lambda inside: None, None,
                               "refused:domain: an aryl radical product")
    res = win.fixed_point(toy(), observe)
    assert res.status == "refused:domain: an aryl radical product"
    assert res.window is None


# --------------------------------------------------------------------------
# items 1, 3, 5 on real certificates (RMG database)
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def rmg():
    import rmg_moves_harness as harness
    if not harness.DATABASE.is_dir():
        pytest.skip(f"needs the RMG database at {harness.DATABASE}")
    import channelsTest as ct
    from rmgpy.kmc import certificates as cm
    fam = ct.database().kinetics.families
    return fam, cm.Certifier(fam, max_channel_sets=16)


@pytest.fixture(scope="module")
def habs_cases(rmg):
    fam, certifier = rmg
    builds = [b.PROVIDERS["polystyrene"].build("pe", b.PsChainSpec(1, tail_end="benzylic_radical")),
              b.PROVIDERS["polystyrene"].build("dp3", b.PsChainSpec(3))]
    return win.member_cases(certifier, fam["H_Abstraction"], builds, "forward", label="PE+DP3")


def _interior(case, unit=1, position=0):
    coords = {a: (p.index, u, k) for p in case.witness.participants
              for a, u, k in p.metadata.coordinates}
    return coords.get(case.moves[0].atom("*1")) == (1, unit, position) and \
        coords.get(case.moves[0].atom("*3")) == (0, 0, 1)


def test_failed_predicate_outside_the_ball_pulls_its_atoms_into_the_window(rmg, habs_cases):
    fam, certifier = rmg
    case = next(c for c in habs_cases if _interior(c))
    cert = certifier.certify(fam["H_Abstraction"], win.template_molecules(case.witness,
                                                                         "forward"),
                             member_labels=case.wanted)
    full = win.certificate_support(cert, case.wanted, classifier=win.SearchClassifier(
        fam["H_Abstraction"]))
    matched = win.certificate_support(cert, case.wanted, include_failed=False,
                                      classifier=win.SearchClassifier(fam["H_Abstraction"]))
    assert full.refusal is None
    extra = set(full.atoms) - set(matched.atoms)
    assert extra, "no failed predicate reaches beyond the matched balls"
    rows = case.witness.reactant_rows()
    beyond = {a for a in extra if any(n in matched.atoms for n, _o in rows[a][5])}
    assert beyond      # one bond outside the matched predicates' support
    failed_only = {a for node, result, atoms in full.sources if result == "failed"
                   for a in atoms} - set(matched.atoms)
    assert beyond <= failed_only
    window = win.close(case.witness, full)
    assert beyond <= window.atoms
    observed = win.MemberObserver(certifier, fam["H_Abstraction"], "forward", case.wanted)
    assert beyond <= set(observed(case.witness).support.atoms)


def test_real_member_reaches_a_recertified_fixed_point(rmg, habs_cases):
    fam, certifier = rmg
    case = next(c for c in habs_cases if _interior(c))
    bound = win.span_bound(win.family_constants(fam["H_Abstraction"]))
    res = win.run_member(certifier, fam["H_Abstraction"], case, bound=bound, contexts=False)
    assert res["status"] == "fixed_point", res
    assert res["span_bound_held"] and not res["mapping_problems"]
    assert not res["boundary_bond_problems"]
    assert res["observations"] >= 1 + len(b.vector_ports(case.witness))
    for index, span in res["span"].items():
        assert span["backbone_atoms"] <= bound.backbone_atoms


def test_domain_refused_member_is_refused_not_truncated(rmg, habs_cases):
    fam, certifier = rmg
    builds = [b.PROVIDERS["polystyrene"].build("pe", b.PsChainSpec(1, tail_end="benzylic_radical")),
              b.PROVIDERS["polystyrene"].build("dp3", b.PsChainSpec(3))]
    stats = {}
    every = win.member_cases(certifier, fam["H_Abstraction"], builds, "forward",
                             skip_ring_changes=False)
    win.member_cases(certifier, fam["H_Abstraction"], builds, "forward", stats=stats)
    ring = win.ring_atoms(every[0].reactants)
    family = fam["H_Abstraction"]

    def touches_ring(case):
        change = win.connectivity_change(family, next(iter(case.wanted)))
        return any(a in ring for _k, pair in change for a in pair)
    aryls = [c for c in every if touches_ring(c)]
    assert all(win.species_domain_refusal(c.witness) for c in aryls)
    # the ring-change mappings the enumeration skips are exactly these members'
    assert stats["ring_change_mappings"] == sum(len(c.moves) for c in aryls) > 0
    assert len(every) - len(aryls) == len(habs_cases)
    aryl = next(c for c in aryls if c.moves[0].atom("*1") in ring)
    res = win.run_member(certifier, fam["H_Abstraction"], aryl, contexts=False)
    assert res["status"] == "refused:domain", res
    assert "window_atoms" not in res
    # the certificate path refuses it too (eligibility)
    observer = win.MemberObserver(certifier, fam["H_Abstraction"], "forward", aryl.wanted)
    obs = observer(aryl.witness)
    assert obs.refusal and obs.refusal.startswith("refused:domain")
    assert win.fixed_point(aryl.witness, observer).window is None


def test_scoped_certificate_equals_the_full_one_on_the_member(rmg, habs_cases):
    """The scoped certifier keeps every mapping of the member and of any other
    pathway with its outcome; its view of the member equals the full
    certificate's."""
    fam, certifier = rmg
    family = fam["H_Abstraction"]
    for case in [c for c in habs_cases if _interior(c)][:2]:
        mols = win.template_molecules(case.witness, "forward")
        full = certifier.certify(family, mols, member_labels=case.wanted)
        scoped = win.scoped_certifier(certifier, family, case.wanted).certify(
            family, mols, member_labels=case.wanted)
        inside = frozenset(case.witness.reactant_rows())
        assert scoped["complete"] and scoped["eligible"]
        assert win.certificate_view(full, case.wanted, inside) == \
            win.certificate_view(scoped, case.wanted, inside)
        assert win.member_entries(full, case.wanted)[0][1]["mappings"] == \
            win.member_entries(scoped, case.wanted)[0][1]["mappings"]
        def ran(cs):
            return sum(1 for m in cs.moves if m.status != "excluded")
        assert ran(certifier.channel_set(family, mols)) > ran(
            win.scoped_certifier(certifier, family, case.wanted).channel_set(family, mols))


# --------------------------------------------------------------------------
# slow: the required supported cases and the R_Recombination ladder control
# --------------------------------------------------------------------------
@pytest.mark.slow
@pytest.mark.skipif(not SLOW, reason="set RMG_KMC_SLOW=1")
def test_required_supported_cases_reach_recertified_fixed_points(rmg):
    """Every selected member of every required case reaches a fixed point (a
    reacting atom on an artificial port is accounted by its status), every
    mapping's edits lie in the window on both sides, the span bound holds,
    the A+A recombination member covers both template slot orders."""
    import windows_acceptance as wa
    fam, certifier = rmg
    reached = {}
    for label, family_name, direction, specs, selection in wa.required_sets():
        family = fam[family_name]
        bound = win.span_bound(win.family_constants(family))
        cases = [c for c in win.member_cases(certifier, family, wa.builds_of(specs), direction)
                 if wa.select(selection, c)]
        assert cases, label
        rows = [win.run_member(certifier, family, c, bound=bound, contexts=False)
                for c in cases]
        statuses = {r["status"] for r in rows}
        assert statuses <= {"fixed_point", "refused:reacting_port"}, (label, statuses)
        for r in rows:
            if r["status"] == "fixed_point":
                assert not r["mapping_problems"] and r["span_bound_held"], (label, r)
        reached[label] = sum(r["status"] == "fixed_point" for r in rows)
        if "identical" in label:
            assert any(len(r["slot_orders"]) == 2 for r in rows)
        print("REQUIRED", label, len(rows), reached[label])
    assert all(reached.values()), reached
    assert len(reached) == 8


@pytest.mark.slow
@pytest.mark.skipif(not SLOW, reason="set RMG_KMC_SLOW=1")
def test_recombination_ladder_control(rmg):
    """Homolysis of the intact closed chains of DP L and L + 1 (R_Recombination):
    no window of DP L + 1 holds both physical ends (non-overlap control)."""
    fam, certifier = rmg
    family = fam["R_Recombination"]
    bound = win.span_bound(win.family_constants(family))
    both = {}
    for dp in (bound.L, bound.L + 1):
        spec = b.PsChainSpec(dp, head_end="methyl", tail_end="benzyl")
        build = b.PROVIDERS["polystyrene"].build("ladder", spec)
        rows = [win.run_member(certifier, family, c, bound=bound, contexts=False)
                for c in win.member_cases(certifier, family, [build], "reverse")]
        reached = [r for r in rows if r["status"] == "fixed_point"]
        assert reached
        both[dp] = sum(any(s["both_ends"] for s in r["span"].values()) for r in reached)
    print("LADDER CONTROL", bound.L, both)
    assert both[bound.L + 1] == 0
