"""Tests for the provider-generic window plumbing (I-165, slice S2a): the
window machinery runs on declared-chain witnesses through provider objects.

Database-free except the tests marked [DB], which load exactly one kinetics
family (plus thermo and prepared rules) through ``one_family.family_db``."""

import hashlib
import json

import pytest

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import certificates as cm
from rmgpy.kmc import declared_chain as dc
from rmgpy.kmc import rmg_moves as rm
from rmgpy.kmc import windows as win

import declared_parity_support as sup
import interior_radical_specs as irs
import one_family
from declaredChainTest import (_atom, _bond, _declaration, cyclopentane_mapping,
                               pmma_mapping, pp_mapping, spec_of)

E = dc.EditOp
HEAD, TAIL = (0, "head"), (0, "tail")
CEILING = 200
needs_db = pytest.mark.skipif(not one_family.available(), reason="needs the RMG database")
# sha256 of json.dumps(cert, sort_keys=True, default=repr) for hexadecyl + methyl
# under R_Recombination at the base (59246f18b): /home/alon/runs/i165/s2a/scripts
BASE_CERT_DIGEST = "767966c8311b3019d752d7f0791ba65ef9ef314296d220fcbf2c1edca0e1d6fd"


@pytest.fixture(scope="module")
def ps():
    return sup.ps_context()


@pytest.fixture(scope="module")
def pe():
    return dc.DeclaredChainProvider(sup.pe_declaration())


@pytest.fixture(scope="module")
def pp():
    return dc.DeclaredChainProvider(spec_of(pp_mapping()))


def declared(provider, chain):
    return win.DeclaredSpec(provider.name, chain)


def pe_key(pe):
    (key,) = pe.geometry.keys
    return key


# --------------------------------------------------------------------------
# 1. witness_spec refusals
# --------------------------------------------------------------------------
@pytest.mark.parametrize("spec", [
    ("alkane", 16, None, None, 7),
    ("polystyrene", 3, None),
    ("small", "X"),
    ("small",),
    dc.ChainSpec(2),
    "alkane",
    ("alkane", "sixteen"),
])
def test_witness_spec_refuses_legacy_and_bare_shapes(spec):
    with pytest.raises(win.WindowError, match="refused:legacy_spec"):
        win.witness_spec(spec)


def test_witness_spec_refuses_an_unregistered_provider():
    with pytest.raises(win.WindowError, match="refused:legacy_spec"):
        win.witness_spec(win.DeclaredSpec("declared:0000000000000000", dc.ChainSpec(2)))


def test_witness_spec_normalizes_the_accepted_shapes(pe):
    assert win.witness_spec(("alkane", 16)) == ("alkane", 16, None, None)
    assert win.witness_spec(("alkane", 16, "methyl")) == ("alkane", 16, "methyl", None)
    assert win.witness_spec(("small", "H")) == ("small", "H")
    spec = b.PsChainSpec(3)
    assert win.witness_spec(spec) is spec
    d = declared(pe, dc.ChainSpec(3))
    assert win.witness_spec(d) == d


def test_interior_radical_tuple_is_never_built_closed():
    with pytest.raises(win.WindowError, match="refused:legacy_spec"):
        win.build_spec(("alkane", 16, None, None, 7))
    for call in (lambda: win.builds_of([("alkane", 16, None, None, 7)]),
                 lambda: win.longer_specs([("alkane", 16, None, None, 7)], [HEAD]),
                 lambda: win.event_specs([("alkane", 16, None, None, 7)],
                                         ((HEAD, win.INTERIOR),))):
        with pytest.raises(win.WindowError, match="refused:legacy_spec"):
            call()


def test_declared_spec_builds_through_its_provider(pp):
    build = win.build_spec(declared(pp, dc.ChainSpec(3, head="declared")), "x")
    assert build.metadata.provider == pp.name
    assert sup.builds_match(build, pp.build("y", dc.ChainSpec(3, head="declared")))


# --------------------------------------------------------------------------
# 2. declared longer_specs: edits and reversed periods
# --------------------------------------------------------------------------
def _pp_reversed(pp, **kw):
    edits = (E("remove_h", 0, "c2"), E("remove_h", 1, "c2"),
             E("set_bond_order", 0, "c2", 2, other=(1, "c2")),
             E("remove_h", 2, "c3"), E("set_radical", 2, "c3", 1))
    return declared(pp, dc.ChainSpec(4, orientation=(1, -1, 1, -1), edits=edits, **kw))


def _shift(c):
    return lambda coord, c=c: (coord[0] - c, coord[1])


@pytest.mark.parametrize("k", [1, 2])
def test_head_move_shifts_edits_and_prepends_orientation(pp, k):
    spec = _pp_reversed(pp)
    (longer,), shift = win.longer_specs([spec], [HEAD], k)
    assert shift == {0: k}
    chain = longer.chain
    assert chain.units == 4 + k
    assert chain.orientation == (1,) * k + (1, -1, 1, -1)
    assert [e.unit for e in chain.edits] == [0 + k, 1 + k, 0 + k, 2 + k, 2 + k]
    assert chain.edits[2].other == (1 + k, "c2")
    expected = dc.ChainSpec(4 + k, orientation=(1,) * k + (1, -1, 1, -1), edits=tuple(
        E(e.kind, e.unit + k, e.key, e.value,
          None if e.other is None else (e.other[0] + k, e.other[1]))
        for e in spec.chain.edits))
    assert chain == expected
    # the longer spec is the witness extended k periods at its head
    w = sup.identity_witness(win.build_spec(spec))
    extended, _ = b.extend_vector(w, (k, 0), ceiling=CEILING)
    direct = sup.identity_witness(win.build_spec(longer))
    assert sup.witnesses_match(direct, extended, _shift(k))


def test_tail_move_keeps_edits_and_appends_orientation(pp):
    spec = _pp_reversed(pp)
    (longer,), shift = win.longer_specs([spec], [TAIL], 1)
    assert shift == {}
    assert longer.chain.orientation == (1, -1, 1, -1, -1)
    assert longer.chain.edits == spec.chain.edits
    w = sup.identity_witness(win.build_spec(spec))
    extended, _ = b.extend_vector(w, (0, 1), ceiling=CEILING)
    assert sup.witnesses_match(sup.identity_witness(win.build_spec(longer)), extended)


def test_empty_orientation_stays_empty(pp):
    (longer,), _ = win.longer_specs([declared(pp, dc.ChainSpec(3))], [HEAD, TAIL], 2)
    assert longer.chain.orientation == () and longer.chain.units == 7


def test_phase_cut_end_refuses_a_move(ps):
    spec = declared(ps.provider, dc.ChainSpec(3, head_phase=ps.phase, head="declared"))
    with pytest.raises(b.ExtensionRefused):
        win.longer_specs([spec], [HEAD])
    (longer,), _ = win.longer_specs([spec], [TAIL])
    assert longer.chain.units == 4 and longer.chain.head_phase == ps.phase


# --------------------------------------------------------------------------
# 3-4. event_specs parity (PS and the alkane cases)
# --------------------------------------------------------------------------
PS_SPECS = (b.PsChainSpec(4), b.PsChainSpec(4, head_end="primary_radical"),
            b.PsChainSpec(4, tail_end="benzylic_radical"))
PS_ASSIGNMENTS = (((HEAD, win.INTERIOR),), ((TAIL, win.INTERIOR),),
                  ((HEAD, "primary_radical"),), ((TAIL, "benzylic_radical"),))
DECLARED_END = {"primary_radical": "radical", "benzylic_radical": "radical",
                "methyl": "declared", "benzyl": "declared", win.INTERIOR: win.INTERIOR}


@pytest.mark.parametrize("spec", PS_SPECS)
@pytest.mark.parametrize("assignment", PS_ASSIGNMENTS)
def test_event_specs_parity_on_ps(ps, spec, assignment):
    twin = declared(ps.provider, sup.legacy_to_chain(spec, ps))
    sup.identity_witness(win.build_spec(twin))         # the witness builds
    legacy, lshift = win.event_specs([spec], assignment)
    mapped = tuple((port, DECLARED_END[end]) for port, end in assignment)
    ours, dshift = win.event_specs([twin], mapped)
    assert lshift == dshift
    assert sup.builds_match(win.build_spec(legacy[0]), win.build_spec(ours[0]), ps.coord_map)


ALKANE_EVENTS = (
    ([("alkane", 16)], ((HEAD, "primary_radical"),)),
    ([("alkane", 16, None, "methyl")], ((HEAD, win.INTERIOR),)),
    ([("alkane", 16, None, "methyl")], ((TAIL, "methyl"),)),
    ([("alkane", 16, "primary_radical", None)], ((TAIL, win.INTERIOR),)),
)


@pytest.mark.parametrize("specs, assignment", ALKANE_EVENTS)
def test_event_specs_parity_on_the_alkane_cases(pe, specs, assignment):
    twin = [declared(pe, sup.legacy_alkane_to_chain(s, pe_key(pe))) for s in specs]
    legacy, lshift = win.event_specs(specs, assignment)
    ours, dshift = win.event_specs(twin, tuple((p, DECLARED_END[e]) for p, e in assignment))
    assert lshift == dshift
    assert sup.builds_match(win.build_spec(legacy[0]), win.build_spec(ours[0]))


@pytest.mark.parametrize("radical, port, k", [(7, HEAD, 1), (7, TAIL, 1), (3, HEAD, 2)])
def test_interior_radical_survives_a_declared_move(pe, radical, port, k):
    spec = ("alkane", 16, None, None, radical)
    twin = declared(pe, sup.legacy_alkane_to_chain(spec, pe_key(pe)))
    (longer,), shift = win.longer_specs([twin], [port], k)
    moved = radical + shift.get(0, 0)
    expected = irs.build_alkyl_radical(16 + k, moved)
    assert sup.builds_match(expected, win.build_spec(longer))


# --------------------------------------------------------------------------
# 5. geometry
# --------------------------------------------------------------------------
def _span_keys(span):
    return {i: {k: v for k, v in s.items() if k != "provider"} for i, s in span.items()}


def _balls(witness, radius):
    rows = witness.reactant_rows()
    heavy = sorted(a for a, row in rows.items() if row[0] != "H")
    out = []
    for centre in heavy[::3]:
        seen, frontier = {centre}, [centre]
        for _ in range(radius):
            frontier = [n for a in frontier for n, _o in rows[a][5] if n not in seen]
            seen |= set(frontier)
        out.append(seen)
    return out


@pytest.mark.parametrize("spec", [b.PsChainSpec(4), b.PsChainSpec(3, radical_end=True),
                                  b.PsChainSpec(3, head_shift=True, head_end="benzyl")])
def test_window_span_matches_legacy_on_ps(ps, spec):
    legacy = sup.identity_witness(ps.legacy.build("l", spec))
    twin = sup.identity_witness(sup.declared_ps_build(ps, spec))
    for atoms_l in _balls(legacy, 3) + [set(legacy.reactant_rows())]:
        # the same atoms on the twin through one marked isomorphism
        atoms_d = _transport(legacy, twin, atoms_l, ps.coord_map)
        assert _span_keys(win.window_span(legacy, atoms_l)) == \
            _span_keys(win.window_span(twin, atoms_d))


def _transport(legacy, twin, atoms, coord_map):
    from networkx.algorithms import isomorphism
    g1, g2 = sup.witness_graph(legacy, coord_map), sup.witness_graph(twin)
    matcher = isomorphism.GraphMatcher(g1, g2, node_match=lambda x, y: x["label"] == y["label"],
                                       edge_match=lambda x, y: x["label"] == y["label"])
    mapping = next(matcher.isomorphisms_iter())
    return {mapping[a] for a in atoms}


@pytest.mark.parametrize("spec", [("alkane", 12), ("alkane", 9, "methyl", "primary_radical")])
def test_window_span_matches_legacy_on_the_alkane(pe, spec):
    legacy = sup.identity_witness(sup.legacy_alkane_build(spec))
    twin = sup.identity_witness(pe.build("d", sup.legacy_alkane_to_chain(spec, pe_key(pe))))
    for atoms in _balls(legacy, 2) + [set(legacy.reactant_rows())]:
        assert _span_keys(win.window_span(legacy, atoms)) == \
            _span_keys(win.window_span(twin, _transport(legacy, twin, atoms, None)))


def test_window_span_follows_the_built_path_after_a_deletion():
    provider = dc.DeclaredChainProvider(spec_of(cyclopentane_mapping()))
    opened = provider.build("x", dc.ChainSpec(2, edits=(
        E("delete_atom", 1, "r2"), E("set_radical", 1, "r1", 1), E("set_radical", 1, "r3", 1))))
    w = sup.identity_witness(opened)
    span = win.window_span(w, set(w.reactant_rows()))
    assert span[0]["backbone_atoms"] == 7 == provider.backbone_atoms(dc.ChainSpec(2, edits=(
        E("delete_atom", 1, "r2"), E("set_radical", 1, "r1", 1), E("set_radical", 1, "r3", 1))))
    pristine = sup.identity_witness(provider.build("p", dc.ChainSpec(2)))
    assert win.window_span(pristine, set(pristine.reactant_rows()))[0]["backbone_atoms"] == 6


def _constants(a):
    # A = D + 2 * (r + 1) + 1 with r = 2
    return win.FamilyConstants("X", radius=2, reacting_distance=a - 7, free_label_search=False)


TODAY = {8: (5, 4), 11: (6, 6), 12: (7, 6), 15: (8, 8), 20: (11, 10)}


@pytest.mark.parametrize("a", sorted(TODAY))
def test_span_bound_with_and_without_a_provider(ps, a):
    c = _constants(a)
    bound = win.span_bound(c)
    assert bound.backbone_atoms == a and (bound.repeats, bound.L) == TODAY[a]
    legacy_ps = b.PROVIDERS["polystyrene"]
    with_ps = win.span_bound(c, legacy_ps)
    assert (with_ps.repeats, with_ps.L) == TODAY[a]
    with_alkane = win.span_bound(c, b.PROVIDERS["linear_alkane"])
    assert (with_alkane.repeats, with_alkane.L) == (a, a)
    with_declared_ps = win.span_bound(c, ps.provider)
    assert (with_declared_ps.repeats, with_declared_ps.L) == TODAY[a]


def test_repeats_bound_per_provider(pe, pp, ps):
    bound = win.span_bound(_constants(11))
    assert win.repeats_bound(bound, "polystyrene") == 6
    assert win.repeats_bound(bound, "linear_alkane") == 11
    assert win.repeats_bound(bound, b.PROVIDERS["polystyrene"]) == 6
    assert win.repeats_bound(bound, pe) == 11 == win.repeats_bound(bound, pe.name)
    assert pe.coordinate_backbone_atoms == 1 and pe.d_rep == 2
    assert win.repeats_bound(bound, pp) == 6 and pp.coordinate_backbone_atoms == 2
    assert win.repeats_bound(bound, ps.provider) == 6
    with pytest.raises(win.WindowError):
        win.repeats_bound(bound, "small")


def test_check_span_skips_small_species():
    build = win.build_spec(("small", "H2"))
    w = sup.identity_witness(build)
    span = win.check_span(w, set(w.reactant_rows()), win.span_bound(_constants(8)))
    assert span[0]["backbone_atoms"] == 0


# --------------------------------------------------------------------------
# 6. classifier materialization
# --------------------------------------------------------------------------
def _types(atoms):
    return sorted(a.atomtype.label for a in atoms)


def test_declared_interior_atoms_match_the_legacy_ps_unit(ps):
    legacy = win._interior_unit_atoms("polystyrene")
    ours = ps.provider.interior_atoms_rmg()
    assert _types(ours) == _types(legacy)
    assert sum(a.is_hydrogen() for a in ours) == 8


def test_pp_interior_atoms(pp):
    atoms = pp.interior_atoms_rmg()
    heavy = [a for a in atoms if not a.is_hydrogen()]
    assert _types(heavy) == ["Cs", "Cs", "Cs"]
    assert sorted(sum(1 for n in a.bonds if n.is_hydrogen()) for a in heavy) == [1, 2, 3]
    assert _types([a for a in atoms if a.is_hydrogen()]) == ["H0"] * 6    # RMG perceives a bonded neutral H as H0
    alkane = b.PROVIDERS["linear_alkane"].interior_atoms_rmg()
    assert _types(alkane) == _types(win._interior_unit_atoms("linear_alkane"))


@needs_db
def test_classifier_is_equal_for_legacy_and_declared_ps(ps):
    fam = one_family.family_db("R_Recombination").kinetics.families["R_Recombination"]
    from rmgpy.kmc import channels as chm
    legacy = win.SearchClassifier(fam, providers=(b.PROVIDERS["polystyrene"],))
    ours = win.SearchClassifier(fam, providers=(ps.provider,))
    assert _types(legacy.fresh) == _types(ours.fresh) and len(ours.fresh) == 16
    checked = 0
    for node, entry in sorted(fam.groups.entries.items()):
        labels = set()
        for group in chm._group_atoms(fam.groups, entry):
            labels |= set(group.get_all_labeled_atoms())
        for label in sorted(labels):
            assert legacy.extensible(node, label) == ours.extensible(node, label), (node, label)
            checked += 1
    assert checked > 10


# --------------------------------------------------------------------------
# 7. no global-registry dependence; fragment facts per provider
# --------------------------------------------------------------------------
def _ps_witness():
    build = b.PROVIDERS["polystyrene"].build("ps", b.PsChainSpec(3))
    return sup.identity_witness(build)


def test_legacy_facts_and_classifier_ignore_constructed_providers():
    w = _ps_witness()
    win._FACTS.clear()
    before = win._fragment_facts(win.observation_providers(w))
    fresh_before = _types(win.SearchClassifier(None, providers=win.observation_providers(w))
                          .fresh)
    dc.DeclaredChainProvider(spec_of(pmma_mapping()))
    win._FACTS.clear()
    after = win._fragment_facts(win.observation_providers(w))
    fresh_after = _types(win.SearchClassifier(None, providers=win.observation_providers(w))
                         .fresh)
    assert before == after and fresh_before == fresh_after
    assert win._fragment_facts() == after
    assert sorted(after) == sorted(p.name for p in b.LEGACY_PROVIDERS)


def test_declared_facts_come_from_the_declared_fragment(pp):
    w = sup.identity_witness(pp.build("w", dc.ChainSpec(3)))
    providers = win.observation_providers(w)
    assert [p.name for p in providers] == [pp.name]
    facts = win._fragment_facts(providers)
    assert facts == {pp.name: {"backbone_head": {"closed_shell": True, "delta": {"C": 3, "H": 6}},
                               "backbone_tail": {"closed_shell": True,
                                                 "delta": {"C": 3, "H": 6}}}}


class _BrokenProvider(b.BoundaryProvider):
    name = "declared:broken"

    def probe_ports(self):
        raise b.ExtensionRefused("no probe")

    def interior_atoms_rmg(self):
        return []


def _guard_certificate():
    predicate = {"node": "N", "result": "passed",
                 "guards": {"radical_count": {"read": True},
                            "element_count": {"read": True, "group_minimum": {"O": 1},
                                              "structure": {}}}}
    mapping = {"labels": [[0, "*1", 1]], "descent": "d"}
    return {"channels": [{"kind": "kept", "members": [{"mappings": [mapping]}]}],
            "descents": {"d": [{"root": predicate}]}, "global_guards": []}, \
        {frozenset({("*1", 1)})}


def test_undecided_fragment_facts_refuse_by_name():
    facts = win._fragment_facts((_BrokenProvider(),))
    assert "declared:broken" in facts
    cert, wanted = _guard_certificate()
    support = win.certificate_support(cert, wanted, classifier=None, fragment_facts=facts)
    assert support.refusal == "refused:fragment_facts_undecided: declared:broken"
    dispositions = {k: v for k, v, _r in support.dispositions}
    assert dispositions[("guard", "radical_count")] == "refused"
    assert dispositions[("guard", "element_count:O")] == "refused"
    decided = win.certificate_support(cert, wanted, classifier=None,
                                      fragment_facts=win._fragment_facts())
    assert not (decided.refusal or "").startswith("refused:fragment_facts_undecided")


# --------------------------------------------------------------------------
# 8. the radical end
# --------------------------------------------------------------------------
@pytest.mark.parametrize("make", ["ps", "pp", "pmma"])
@pytest.mark.parametrize("port", [HEAD, TAIL])
def test_radical_end_reconstruction(ps, pp, make, port):
    provider = {"ps": ps.provider, "pp": pp,
                "pmma": dc.DeclaredChainProvider(spec_of(pmma_mapping()))}[make]
    assert "radical" in provider.END_CHEMISTRIES["backbone_" + port[1]]
    w = sup.identity_witness(provider.build("o", dc.ChainSpec(3)))
    built, _record = b.reconstruct_end(w, port, "radical")
    direct = provider.build("d", dc.ChainSpec(3, **{port[1]: "radical"}))
    assert sup.witnesses_match(sup.identity_witness(direct), built)
    assert not built.port(port).artificial


def test_context_assignments_read_radical_ends_per_provider(pp):
    w = sup.identity_witness(pp.build("o", dc.ChainSpec(3)))
    alone = win.context_assignments(w, near_ports=())
    assert ((HEAD, "radical"), (TAIL, "radical")) in alone
    singles = sorted(e for a in alone if sum(e is not None for _k, e in a) == 1
                     for _k, e in a if e)
    assert singles == ["declared", "declared", "radical", "radical"]
    # far non-radical joint combination (both declared) is a context
    assert ((HEAD, "declared"), (TAIL, "declared")) in alone
    legacy = _ps_witness()
    far = win.context_assignments(legacy, near_ports=())
    assert ((HEAD, None), (TAIL, "benzylic_radical")) in far
    assert ((HEAD, "primary_radical"), (TAIL, "benzylic_radical")) in far
    assert sum(1 for a in far if a == ((HEAD, None), (TAIL, "benzyl"))) == 1


def test_far_port_verdicts_skip_declared_radical_ends(pp):
    w = sup.identity_witness(pp.build("o", dc.ChainSpec(3)))
    atoms = frozenset(w.reactant_rows())
    window = win.Window(atoms=atoms, support=win.Support(atoms=atoms))

    def at(assignment, view):
        obs = win.Observation(support=window.support, viewer=lambda inside, v=view: v)
        return win.ContextResult(assignment, "fixed_point", "", w, window, obs)
    results = [at(((HEAD, None), (TAIL, None)), "a"), at(((HEAD, "radical"), (TAIL, None)), "b"),
               at(((HEAD, "declared"), (TAIL, None)), "a")]
    assert win.far_port_verdicts(results, far=(HEAD, TAIL)) == {(HEAD, "declared"): True}


# --------------------------------------------------------------------------
# 9. roster refusals
# --------------------------------------------------------------------------
def _rosters(spec):
    c = win.FamilyConstants("X", radius=2, reacting_distance=1, free_label_search=False,
                            classes=())
    pattern = win.Pattern("chain", ("methyl", "benzyl"), ())
    return {
        "kind_specs": lambda: win.kind_specs("closed", 2, provider=spec),
        "LONG_SPECS": lambda: win.LONG_SPECS["closed"](spec),
        "pair_roster": lambda: win.pair_roster(None, None, provider=spec),
        "participant_patterns": lambda: win.participant_patterns(provider=spec),
        "pattern_specs": lambda: win.pattern_specs(pattern, 2, provider=spec),
        "_realizes": lambda: win._realizes(pattern, spec),
        "site_classes": lambda: win.site_classes(provider=spec),
        "class_roster": lambda: win.class_roster(None, None, provider=spec),
        "pattern_site_kinds": lambda: win.pattern_site_kinds(pattern, spec),
        "ladder_plan": lambda: win.ladder_plan(c, provider=spec),
        "close_classes": lambda: win.close_classes(None, None, (), provider=spec),
    }


@pytest.mark.parametrize("name", sorted(_rosters(None)))
def test_rosters_refuse_declared_specs(pp, name):
    spec = declared(pp, dc.ChainSpec(3))
    with pytest.raises(win.WindowError, match="refused:roster: declared providers are not in "
                                              "the PS roster"):
        _rosters(spec)[name]()


def test_rosters_refuse_an_all_declared_witness(pp):
    w = sup.identity_witness(pp.build("o", dc.ChainSpec(3)))
    with pytest.raises(win.WindowError, match="refused:roster"):
        win.kind_specs("closed", 2, provider=w)
    with pytest.raises(win.WindowError, match="refused:roster"):
        win.pair_roster(None, None, provider=pp)


def test_rosters_still_build_ps_without_a_declared_argument():
    assert win.kind_specs("closed", 2) == win.kind_specs("closed", 2, provider=None)
    assert win.LONG_SPECS["closed"]() == (b.PsChainSpec(win.LONG_UNITS),)


# --------------------------------------------------------------------------
# 10. depository growth unit [DB: R_Recombination]
# --------------------------------------------------------------------------
def _c16_methyl():
    c16 = b.PROVIDERS["linear_alkane"].build(16, tail_end="primary_radical")
    mol16 = rm.molecule_from_rows(dict(b.rows_of_molecule(c16.molecule)))
    methyl = rm.molecule_from_rows({
        1000: ("C", None, 0, 1, 0, ((1001, 1.0), (1002, 1.0), (1003, 1.0))),
        1001: ("H", None, 0, 0, 0, ((1000, 1.0),)), 1002: ("H", None, 0, 0, 0, ((1000, 1.0),)),
        1003: ("H", None, 0, 0, 0, ((1000, 1.0),))})
    return [mol16, methyl]


def _digest(cert):
    return hashlib.sha256(json.dumps(cert, sort_keys=True, default=repr).encode()).hexdigest()


@needs_db
def test_growth_unit_is_a_parameter(pe):
    fam = one_family.family_db("R_Recombination").kinetics.families["R_Recombination"]
    certifier = cm.Certifier(fam)
    base = certifier.certify(fam, _c16_methyl())
    assert _digest(base) == BASE_CERT_DIGEST
    assert _digest(cm.Certifier(fam).certify(fam, _c16_methyl(), repeat=None)) == \
        BASE_CERT_DIGEST
    grown = cm.Certifier(fam).certify(fam, _c16_methyl(), repeat=pe.unit_formula())
    exclusions = [ch["depository_exclusion"] for ch in grown["channels"]
                  if "depository_exclusion" in ch]
    assert exclusions and all(x["repeat_unit"] == {"C": 1, "H": 2} for x in exclusions)


class _Stop(Exception):
    pass


class _RecordingCertifier:
    def __init__(self):
        self.calls = []

    def certify(self, family, molecules, **kw):
        self.calls.append(kw)
        raise _Stop


def _observe(witness):
    certifier = _RecordingCertifier()
    observer = win.MemberObserver(certifier, None, "forward", (frozenset(),), scoped=False)
    with pytest.raises(_Stop):
        observer(witness)
    return certifier.calls[0]


def test_member_observer_passes_the_growth_unit(pe, pp):
    """Correction C-1 (I-165 S2b): legacy witnesses take their provider's
    unit too (S2a passed None for both: PS C8H8 unchanged, alkane was C8H8
    by default and is CH2 now)."""
    assert _observe(_ps_witness())["repeat"] == (("C", 8), ("H", 8))
    toy = sup.identity_witness(b.PROVIDERS["linear_alkane"].build(8))
    assert _observe(toy)["repeat"] == (("C", 1), ("H", 2))
    assert _observe(sup.identity_witness(pe.build("w", dc.ChainSpec(5))))["repeat"] == \
        (("C", 1), ("H", 2))
    assert _observe(sup.identity_witness(pp.build("w", dc.ChainSpec(3))))["repeat"] == \
        (("C", 3), ("H", 6))


def test_member_observer_refuses_two_chain_providers(pe, pp):
    molecules = [pe.build("a", dc.ChainSpec(3)), pp.build("b", dc.ChainSpec(2))]
    mols, parts = b.witness_reactants(molecules)
    rows = tuple(sorted(rm.molecule_rows(mols[0]).items()) +
                 sorted(rm.molecule_rows(mols[1]).items()))
    w = b.make_witness(rows, rows, parts)
    observer = win.MemberObserver(_RecordingCertifier(), None, "forward", (frozenset(),),
                                  scoped=False)
    with pytest.raises(win.WindowError, match="two chain providers"):
        observer(w)


# --------------------------------------------------------------------------
# rework 1: a portless declared monomer has a span
# --------------------------------------------------------------------------
def _heavy_of(witness):
    rows = witness.reactant_rows()
    return rows, {a for a, row in rows.items() if row[0] != "H"}


def test_portless_ps_monomer_span_matches_the_legacy_monomer(ps):
    legacy = sup.identity_witness(ps.legacy.build("s", b.PsChainSpec(1, monomer=True)))
    ours = sup.identity_witness(ps.provider.monomer())
    want = _span_keys(win.window_span(legacy, set(legacy.reactant_rows())))
    got = _span_keys(win.window_span(ours, set(ours.reactant_rows())))
    assert got == want
    assert got[0]["backbone_atoms"] == 2 and got[0]["units"] == 1
    bound = win.span_bound(_constants(8))
    win.check_span(ours, set(ours.reactant_rows()), bound)


def _propylene():
    atoms = [{"key": k, "element": "C", "hydrogens": h,
              "electronic_state": {"radical_electrons": 0}}
             for k, h in (("m1", 2), ("m2", 1), ("m3", 3))]
    return {"atoms": atoms, "bonds": [{"atoms": ["m1", "m2"], "order": 2},
                                      {"atoms": ["m2", "m3"], "order": 1}],
            "electronic_state": {"multiplicity": 1}}


def test_portless_pp_monomer_backbone_is_the_vinyl_pair():
    mapping = pp_mapping()
    mapping["monomer_product"] = _propylene()
    provider = dc.DeclaredChainProvider(spec_of(mapping))
    w = sup.identity_witness(provider.monomer())
    rows, heavy = _heavy_of(w)
    span = win.window_span(w, set(rows))
    assert span[0]["backbone_atoms"] == 2 and span[0]["units"] == 1
    assert span[0]["physical_ends"] == ()
    # the methyl carbon (three H) is the pendant, not the backbone
    methyl = next(a for a in heavy if sum(rows[n][0] == "H" for n, _o in rows[a][5]) == 3)
    vinyl = next(a for a in heavy if sum(rows[n][0] == "H" for n, _o in rows[a][5]) == 2)
    assert win.window_span(w, {methyl})[0]["backbone_atoms"] == 0
    assert win.window_span(w, {vinyl})[0]["backbone_atoms"] == 1
    assert win.window_span(w, {vinyl, methyl})[0]["backbone_atoms"] == 1


# --------------------------------------------------------------------------
# rework 2: the monomer backbone is bound once, from the declaration, by the
# vinyl-addition relation, and carried by atom id
# --------------------------------------------------------------------------
def _atom_of(spec):
    """(key, element, hydrogens, lone pairs[, {isotope, charge}])."""
    k, e, h, lp = spec[:4]
    atom = _atom(k, e, hydrogens=h, lone_pairs=lp)
    atom.update(spec[4] if len(spec) > 4 else {})
    return atom


def _graph(atoms, bonds, multiplicity=1):
    """atoms: (key, element, hydrogens, lone pairs[, extra]); bonds: (x, y, order)."""
    return {"atoms": [_atom_of(a) for a in atoms],
            "bonds": [_bond(x, y, o) for x, y, o in bonds],
            "electronic_state": {"multiplicity": multiplicity}}


def _explicit(atoms, bonds):
    """The same graph with every hydrogen an explicit atom."""
    out_atoms, out_bonds = [], list(bonds)
    for k, e, h, lp, *extra in atoms:
        out_atoms.append((k, e, 0, lp, *extra))
        for i in range(h):
            out_atoms.append((f"h_{k}_{i}", "H", 0, 0))
            out_bonds.append((k, f"h_{k}_{i}", 1))
    return out_atoms, out_bonds


def _provider(repeat, head, tail, monomer, elements=("C", "H"), lone_pairs=(0,), **bounds):
    atoms, bonds = repeat
    mapping = _declaration([_atom_of(a) for a in atoms],
                           [_bond(x, y, o) for x, y, o in bonds], head, tail,
                           elements=elements, lone_pairs=lone_pairs)
    mapping["locality_policy"]["bounds"].update(bounds)
    mapping["monomer_product"] = _graph(*monomer)
    return dc.DeclaredChainProvider(spec_of(mapping))


def _signature(rows, a):
    """A labelling-free name of a heavy atom: element, H count and its heavy
    neighbours' (element, H count)."""
    def own(x):
        return rows[x][0], sum(1 for n, _o in rows[x][5] if rows[n][0] == "H")
    return own(a), tuple(sorted(own(n) for n, _o in rows[a][5] if rows[n][0] != "H"))


def _monomer_spans(witness):
    """{signature: backbone atoms} of every one-atom and every bonded
    two-atom heavy window, and the whole span."""
    rows, heavy = _heavy_of(witness)
    out = {"all": win.window_span(witness, set(rows))[0]["backbone_atoms"]}
    for a in heavy:
        out[_signature(rows, a)] = win.window_span(witness, {a})[0]["backbone_atoms"]
        for n, _o in rows[a][5]:
            if n in heavy and n > a:
                key = tuple(sorted((_signature(rows, a), _signature(rows, n))))
                out[key] = win.window_span(witness, {a, n})[0]["backbone_atoms"]
    return out


PP_REPEAT = ([("c1", "C", 2, 0), ("c2", "C", 1, 0), ("c3", "C", 3, 0)],
             [("c1", "c2", 1), ("c2", "c3", 1)])
PROPYLENE = ([("m1", "C", 2, 0), ("m2", "C", 1, 0), ("m3", "C", 3, 0)],
             [("m1", "m2", 2), ("m2", "m3", 1)])
PROPYLENE_RELABELLED = ([("z", "C", 2, 0), ("b", "C", 1, 0), ("a", "C", 3, 0)],
                        [("z", "b", 2), ("b", "a", 1)])


def _pp_spans(repeat, monomer):
    return _monomer_spans(sup.identity_witness(_provider(repeat, "c1", "c2", monomer).monomer()))


def test_pp_monomer_spans_are_independent_of_h_representation_and_labels():
    implicit = _pp_spans(PP_REPEAT, PROPYLENE)
    assert implicit["all"] == 2
    assert sorted(v for k, v in implicit.items() if k != "all") == [0, 1, 1, 1, 2]
    methyl = next(k for k in implicit if k != "all" and k[0] == ("C", 3))
    assert implicit[methyl] == 0
    assert _pp_spans(_explicit(*PP_REPEAT), PROPYLENE) == implicit
    assert _pp_spans(PP_REPEAT, _explicit(*PROPYLENE)) == implicit
    assert _pp_spans(_explicit(*PP_REPEAT), _explicit(*PROPYLENE)) == implicit
    assert _pp_spans(PP_REPEAT, PROPYLENE_RELABELLED) == implicit
    assert _pp_spans(_explicit(*PP_REPEAT), _explicit(*PROPYLENE_RELABELLED)) == implicit
    assert _pp_spans(_explicit(*PP_REPEAT), PROPYLENE_RELABELLED) == implicit


PB12_REPEAT = ([("c1", "C", 2, 0), ("c2", "C", 1, 0), ("c3", "C", 1, 0), ("c4", "C", 2, 0)],
               [("c1", "c2", 1), ("c2", "c3", 1), ("c3", "c4", 2)])


@pytest.mark.parametrize("keys", [("m1", "m2", "m3", "m4"), ("q", "c", "x", "a")])
def test_symmetric_butadiene_monomer_is_refused_under_any_labelling(keys):
    # 1,2-polybutadiene: the repeat unit with c1=c2 is butadiene, whose two
    # vinyl pairs are exchanged by its symmetry, so the backbone pair is not
    # determined; every labelling gives the same named refusal
    w, x, y, z = keys
    monomer = ([(w, "C", 2, 0), (x, "C", 1, 0), (y, "C", 1, 0), (z, "C", 2, 0)],
               [(w, x, 2), (x, y, 1), (y, z, 2)])
    provider = _provider(PB12_REPEAT, "c1", "c2", monomer)
    witness = sup.identity_witness(provider.monomer())
    with pytest.raises(b.MetadataError, match="refused:monomer_backbone_ambiguous"):
        win.window_span(witness, set(witness.reactant_rows()))


def test_pe_ethylene_monomer_backbone_is_both_carbons():
    from melt_engine_support import pe_mapping
    ethylene = _graph([("e1", "C", 2, 0), ("e2", "C", 2, 0)], [("e1", "e2", 2)])
    provider = dc.DeclaredChainProvider(spec_of(pe_mapping({10: 1}, transport=False,
                                                           monomer=ethylene)))
    assert provider.geometry.period == 1        # the primitive period is one CH2
    spans = _monomer_spans(sup.identity_witness(provider.monomer()))
    assert spans["all"] == 2


def _abstract_h(witness, carbon):
    """The witness after one H of ``carbon`` is abstracted (the carbon
    becomes a radical); metadata kept."""
    from dataclasses import replace
    rows = dict(witness.reactant_rows())
    h = min(n for n, _o in rows[carbon][5] if rows[n][0] == "H")
    del rows[h]
    row = rows[carbon]
    rows[carbon] = row[:3] + (row[3] + 1,) + row[4:5] + (
        tuple(e for e in row[5] if e[0] != h),)
    (part,) = witness.participants
    part = replace(part, atom_ids=tuple(a for a in part.atom_ids if a != h))
    rows = tuple(sorted(rows.items()))
    return b.make_witness(rows, rows, (part,))


def _backbone_of(witness):
    (part,) = witness.participants
    provider = b.PROVIDERS[part.metadata.provider]
    rows = {a: r for a, r in witness.reactant_rows().items() if a in set(part.atom_ids)}
    return provider.witness_backbone_indices(rows, part.metadata)


def test_propylene_backbone_survives_methyl_h_abstraction():
    provider = _provider(PP_REPEAT, "c1", "c2", PROPYLENE)
    before = sup.identity_witness(provider.monomer())
    rows, heavy = _heavy_of(before)
    methyl = next(a for a in heavy if _signature(rows, a)[0] == ("C", 3))
    after = _abstract_h(before, methyl)
    assert _backbone_of(after) == _backbone_of(before)
    assert len(_backbone_of(before)) == 2 and methyl not in _backbone_of(before)
    full = set(after.reactant_rows())
    assert win.window_span(after, full)[0]["backbone_atoms"] == 2
    assert win.window_span(after, {methyl})[0]["backbone_atoms"] == 0


MMA_REPEAT = ([("c1", "C", 2, 0), ("c2", "C", 0, 0), ("c3", "C", 3, 0), ("c4", "C", 0, 0),
               ("o1", "O", 0, 2), ("o2", "O", 0, 2), ("c5", "C", 3, 0)],
              [("c1", "c2", 1), ("c2", "c3", 1), ("c2", "c4", 1), ("c4", "o1", 2),
               ("c4", "o2", 1), ("o2", "c5", 1)])
MMA = ([("m1", "C", 2, 0), ("m2", "C", 0, 0), ("m3", "C", 3, 0), ("m4", "C", 0, 0),
        ("q1", "O", 0, 2), ("q2", "O", 0, 2), ("m5", "C", 3, 0)],
       [("m1", "m2", 2), ("m2", "m3", 1), ("m2", "m4", 1), ("m4", "q1", 2),
        ("m4", "q2", 1), ("q2", "m5", 1)])


def test_mma_backbone_survives_o_methyl_h_abstraction():
    provider = _provider(MMA_REPEAT, "c1", "c2", MMA, elements=("C", "H", "O"),
                         lone_pairs=(0, 2))
    before = sup.identity_witness(provider.monomer())
    rows, heavy = _heavy_of(before)
    o_methyl = next(a for a in heavy if _signature(rows, a)[0] == ("C", 3)
                    and _signature(rows, a)[1] == (("O", 0),))
    after = _abstract_h(before, o_methyl)
    assert _backbone_of(after) == _backbone_of(before)
    assert len(_backbone_of(before)) == 2
    assert win.window_span(after, set(after.reactant_rows()))[0]["backbone_atoms"] == 2
    assert win.window_span(after, {o_methyl})[0]["backbone_atoms"] == 0


def test_ring_opening_monomer_has_no_span():
    # poly(ethylene oxide) and ethylene oxide: not related by vinyl addition
    peo = ([("c1", "C", 2, 0), ("c2", "C", 2, 0), ("o1", "O", 0, 2)],
           [("c1", "c2", 1), ("c2", "o1", 1)])
    oxirane = ([("a", "C", 2, 0), ("b", "C", 2, 0), ("o", "O", 0, 2)],
               [("a", "b", 1), ("b", "o", 1), ("o", "a", 1)])
    provider = _provider(peo, "c1", "o1", oxirane, elements=("C", "H", "O"),
                         lone_pairs=(0, 2))
    witness = sup.identity_witness(provider.monomer())
    with pytest.raises(b.MetadataError, match="refused:monomer_backbone_undefined"):
        win.window_span(witness, set(witness.reactant_rows()))


def _with_ports(witness, ports):
    from dataclasses import replace
    (part,) = witness.participants
    meta = b.BoundaryMetadata(part.metadata.provider, ports, part.metadata.coordinates)
    return replace(witness, participants=(replace(part, metadata=meta),))


def test_one_port_and_side_only_participants_are_refused(pp):
    w = sup.identity_witness(pp.build("w", dc.ChainSpec(6)))
    (part,) = w.participants
    head = part.metadata.port("head")
    side = b.Port("side", head.atom_id, (), "pendant", "physical")
    for ports in ((head,), (side,)):
        broken = _with_ports(w, ports)
        with pytest.raises(b.MetadataError, match="refused:backbone_ports"):
            win.window_span(broken, set(broken.reactant_rows()))


# --------------------------------------------------------------------------
# rework 3: the stamp survives extension; the binding matches labelled bonds
# across resonance forms
# --------------------------------------------------------------------------
def _addition_witness():
    """A PP radical chain (head artificial) adding to propylene: the radical
    end bonds the monomer's CH2, whose double bond opens onto its CH."""
    provider = _provider(PP_REPEAT, "c1", "c2", PROPYLENE)
    chain = provider.build("r", dc.ChainSpec(3, tail="radical"))
    molecules, parts = b.witness_reactants([chain, provider.monomer()])
    rows = {}
    for m in molecules:
        rows.update(b.rows_of_molecule(m))
    radical = parts[0].metadata.port("tail").atom_id
    ch2, ch = parts[1].metadata.monomer_backbone
    product = {a: list(r) for a, r in rows.items()}

    def bonds(a, changes):
        edges = dict(product[a][5])
        edges.update(changes)
        product[a][5] = tuple(sorted(edges.items()))
    bonds(radical, {ch2: 1.0})
    bonds(ch2, {radical: 1.0, ch: 1.0})
    bonds(ch, {ch2: 1.0})
    product[radical][3] -= 1
    product[ch][3] += 1
    side = tuple(sorted(rows.items()))
    after = tuple(sorted((a, tuple(r)) for a, r in product.items()))
    return provider, b.make_witness(side, after, parts)


def test_extension_keeps_the_monomer_stamp():
    provider, w = _addition_witness()
    stamp = w.participants[1].metadata.monomer_backbone
    # relabelled with its participant (second participant: non-identity ids)
    assert stamp != provider.monomer().metadata.monomer_backbone
    assert set(stamp) <= set(w.participants[1].atom_ids)
    before = win.window_span(w, set(w.reactant_rows()))[1]
    longer, _record = b.extend(w, (0, "head"))
    assert longer.participants[1].metadata.monomer_backbone == stamp
    after = win.window_span(longer, set(longer.reactant_rows()))
    assert after[1] == before and after[1]["backbone_atoms"] == 2
    assert after[0]["units"] == 4


def test_identity_map_keeps_the_monomer_stamp():
    _provider_, w = _addition_witness()
    same = b.apply_map(w, b.identity_map(w))
    assert [p.metadata for p in same.participants] == [p.metadata for p in w.participants]
    assert same.participants[1].metadata.monomer_backbone == \
        w.participants[1].metadata.monomer_backbone


def test_deleted_backbone_atoms_shrink_the_stamp():
    provider = _provider(PP_REPEAT, "c1", "c2", PROPYLENE)
    build = provider.monomer()
    meta = build.metadata
    first, second = meta.monomer_backbone
    one = meta.carried(meta.ports, meta.coordinates, deleted={first})
    assert one.monomer_backbone == (second,)
    rows = dict(b.rows_of_molecule(build.molecule))
    none = meta.carried(meta.ports, meta.coordinates, deleted={first, second})
    assert none.monomer_backbone == ()
    with pytest.raises(b.MetadataError, match="refused:monomer_backbone_deleted"):
        provider.witness_backbone_indices(rows, none)


def test_monomer_indices_cover_only_the_rows_given():
    provider = _provider(PP_REPEAT, "c1", "c2", PROPYLENE)
    build = provider.monomer()
    rows = dict(b.rows_of_molecule(build.molecule))
    first, second = build.metadata.monomer_backbone
    del rows[first]
    assert provider.witness_backbone_indices(rows, build.metadata) == {second: 1}


def test_monomer_backbone_ids_must_carry_coordinates():
    provider = _provider(PP_REPEAT, "c1", "c2", PROPYLENE)
    build = provider.monomer()
    rows = dict(b.rows_of_molecule(build.molecule))
    hydrogen = min(a for a, r in rows.items() if r[0] == "H")
    meta = build.metadata
    bad = b.BoundaryMetadata(meta.provider, (), meta.coordinates,
                             monomer_backbone=(meta.monomer_backbone[0], hydrogen))
    with pytest.raises(b.MetadataError, match="monomer backbone atoms"):
        bad.check(rows)


# the spar pair: the closed repeat unit is C=CC1=C(C)C=C1, the declared
# monomer C=CC1=CC=C1C (same skeleton and bond-order sums, different bonds)
CBD_REPEAT = ([("c1", "C", 2, 0), ("c2", "C", 1, 0), ("r1", "C", 0, 0), ("r2", "C", 0, 0),
               ("r3", "C", 1, 0), ("r4", "C", 1, 0), ("me", "C", 3, 0)],
              [("c1", "c2", 1), ("c2", "r1", 1), ("r1", "r2", 2), ("r2", "r3", 1),
               ("r3", "r4", 2), ("r4", "r1", 1), ("r2", "me", 1)])
CBD_OTHER = ([("a", "C", 2, 0), ("b", "C", 1, 0), ("s1", "C", 0, 0), ("s2", "C", 1, 0),
              ("s3", "C", 1, 0), ("s4", "C", 0, 0), ("n", "C", 3, 0)],
             [("a", "b", 2), ("b", "s1", 1), ("s1", "s2", 2), ("s2", "s3", 1),
              ("s3", "s4", 2), ("s4", "s1", 1), ("s4", "n", 1)])
CBD_SAME = ([("a", "C", 2, 0), ("b", "C", 1, 0), ("s1", "C", 0, 0), ("s2", "C", 0, 0),
             ("s3", "C", 1, 0), ("s4", "C", 1, 0), ("n", "C", 3, 0)],
            [("a", "b", 2), ("b", "s1", 1), ("s1", "s2", 2), ("s2", "s3", 1),
             ("s3", "s4", 2), ("s4", "s1", 1), ("s2", "n", 1)])


def test_bond_positions_matter_beyond_order_sums():
    other = _provider(CBD_REPEAT, "c1", "c2", CBD_OTHER)
    assert other.monomer_binding.keys is None
    assert other.monomer_binding.refusal.startswith("refused:monomer_backbone_undefined")
    w = sup.identity_witness(other.monomer())
    with pytest.raises(b.MetadataError, match="refused:monomer_backbone_undefined"):
        win.window_span(w, set(w.reactant_rows()))
    assert _provider(CBD_REPEAT, "c1", "c2", CBD_SAME).monomer_binding.keys == ("a", "b")


ACRYL_REPEAT = ([("c1", "C", 2, 0), ("c2", "C", 1, 0), ("c3", "C", 0, 0), ("o", "O", 0, 2),
                 ("n", "N", 2, 1)],
                [("c1", "c2", 1), ("c2", "c3", 1), ("c3", "o", 2), ("c3", "n", 1)])
ACRYLAMIDE = ([("a", "C", 2, 0), ("b", "C", 1, 0), ("c", "C", 0, 0), ("o", "O", 0, 2),
               ("n", "N", 2, 1)],
              [("a", "b", 2), ("b", "c", 1), ("c", "o", 2), ("c", "n", 1)])
ACRYLAMIDE_SEPARATED = ([("a", "C", 2, 0), ("b", "C", 1, 0), ("c", "C", 0, 0),
                         ("o", "O", 0, 3, {"charge": -1}), ("n", "N", 2, 0, {"charge": 1})],
                        [("a", "b", 2), ("b", "c", 1), ("c", "o", 1), ("c", "n", 2)])
CHARGED = dict(elements=("C", "H", "N", "O"), lone_pairs=(0, 1, 2, 3),
               formal_charges=[-1, 0, 1])


def test_charge_separated_acrylamide_binds_like_the_neutral_form():
    neutral = _provider(ACRYL_REPEAT, "c1", "c2", ACRYLAMIDE, **CHARGED)
    separated = _provider(ACRYL_REPEAT, "c1", "c2", ACRYLAMIDE_SEPARATED, **CHARGED)
    assert neutral.monomer_binding.keys == ("a", "b")
    assert separated.monomer_binding.keys == ("a", "b")
    assert _monomer_spans(sup.identity_witness(separated.monomer())) == \
        _monomer_spans(sup.identity_witness(neutral.monomer()))


def _kekule_ps(mapping, graph):
    ring = [k for k in ("p1", "p2", "p3", "p4", "p5", "p6")]
    for bond in mapping[graph]["bonds"]:
        x, y = sorted(bond["atoms"])
        if x in ring and y in ring and bond["order"] == 1.5:
            i, j = ring.index(x), ring.index(y)
            bond["order"] = 2 if {i, j} in ({0, 1}, {2, 3}, {4, 5}) else 1
    return mapping


def test_ps_aromatic_and_kekule_writings_bind_identically():
    providers = []
    for kekule in ((), ("monomer_product",), ("repeat_unit",),
                   ("monomer_product", "repeat_unit")):
        mapping = sup.ps_declaration_mapping()
        for graph in kekule:
            mapping = _kekule_ps(mapping, graph)
        providers.append(dc.DeclaredChainProvider(spec_of(mapping)))
    keys = {p.monomer_binding.keys for p in providers}
    assert keys == {("a", "b")}
    spans = [_monomer_spans(sup.identity_witness(p.monomer())) for p in providers]
    assert all(x == spans[0] for x in spans) and spans[0]["all"] == 2


def _labelled(atoms, key):
    return [a + ({"isotope": 13},) if a[0] == key else a for a in atoms]


def test_isotope_labels_must_agree():
    isotopes = dict(isotopes=[-1, 13])
    labelled_monomer = (_labelled(PROPYLENE[0], "m3"), PROPYLENE[1])
    labelled_repeat = (_labelled(PP_REPEAT[0], "c3"), PP_REPEAT[1])
    mismatch = _provider(PP_REPEAT, "c1", "c2", labelled_monomer, **isotopes)
    assert mismatch.monomer_binding.refusal.startswith("refused:monomer_backbone_undefined")
    both = _provider(labelled_repeat, "c1", "c2", labelled_monomer, **isotopes)
    assert both.monomer_binding.keys == ("m1", "m2")
    wrong_end = (_labelled(PROPYLENE[0], "m1"), PROPYLENE[1])
    assert _provider(labelled_repeat, "c1", "c2", wrong_end, **isotopes
                     ).monomer_binding.keys is None


def test_open_shell_monomer_is_refused_by_name():
    from rmgpy.polymer_spec import GraphDeclaration
    unit = spec_of(pp_mapping()).repeat_unit
    radical = _graph([("m1", "C", 2, 0), ("m2", "C", 1, 0), ("m3", "C", 2, 0)],
                     [("m1", "m2", 2), ("m2", "m3", 1)], multiplicity=2)
    radical["atoms"][2]["electronic_state"] = {"radical_electrons": 1}
    binding = dc.bind_monomer(unit, GraphDeclaration.from_mapping(radical))
    assert binding.keys is None
    assert binding.refusal.startswith("refused:monomer_backbone_open_shell")


# --------------------------------------------------------------------------
# rework 1: the growth unit reaches recertification
# --------------------------------------------------------------------------
class _SpyCertifier:
    """Records the keywords of certify and recertify_certificate; returns a
    minimal certificate over the witness atoms."""

    def __init__(self):
        self.certify_calls, self.recertify_calls = [], []

    @staticmethod
    def _cert(molecules):
        return {"witness": {"atoms": [[a.id] for m in molecules for a in m.atoms]}}

    def certify(self, family, molecules, **kw):
        self.certify_calls.append(kw)
        return self._cert(molecules)

    def recertify_certificate(self, family, old, molecules, atom_map, **kw):
        self.recertify_calls.append(kw)
        return {"certificate": self._cert(molecules), "refusals": []}


def _observe_twice(witness, monkeypatch):
    from types import SimpleNamespace
    monkeypatch.setattr(win, "certificate_refusal", lambda cert: None)
    monkeypatch.setattr(win, "member_refusal", lambda cert, wanted: None)
    monkeypatch.setattr(win, "certificate_support",
                        lambda *a, **k: SimpleNamespace(refusal=None, atoms=frozenset()))
    spy = _SpyCertifier()
    observer = win.MemberObserver(spy, None, "forward", (frozenset(),), scoped=False)
    observer(witness)
    port = next(p for p in b.vector_ports(witness) if witness.port(p).artificial)
    longer, _record = b.extend(witness, port)
    observer(longer)
    assert len(spy.certify_calls) == 1 and len(spy.recertify_calls) == 1
    return spy


def test_recertification_receives_the_declared_growth_unit(pe, monkeypatch):
    spy = _observe_twice(sup.identity_witness(pe.build("w", dc.ChainSpec(5))), monkeypatch)
    assert spy.certify_calls[0]["repeat"] == (("C", 1), ("H", 2))
    assert spy.recertify_calls[0]["repeat"] == (("C", 1), ("H", 2))
    # the declared policy's domain reaches both (I-165 S2b)
    assert spy.certify_calls[0]["domain"] is win.DECLARED_INTERIM.domain
    assert spy.recertify_calls[0]["domain"] is win.DECLARED_INTERIM.domain


def test_legacy_recertification_receives_the_provider_growth_unit(monkeypatch):
    """Correction C-1 (S2a: both None); the legacy policy passes no domain."""
    spy = _observe_twice(_ps_witness(), monkeypatch)
    assert spy.certify_calls[0]["repeat"] == (("C", 8), ("H", 8))
    assert spy.recertify_calls[0]["repeat"] == (("C", 8), ("H", 8))
    assert "domain" not in spy.certify_calls[0] and "domain" not in spy.recertify_calls[0]
    toy = sup.identity_witness(b.PROVIDERS["linear_alkane"].build(8))
    spy = _observe_twice(toy, monkeypatch)
    assert spy.certify_calls[0]["repeat"] == spy.recertify_calls[0]["repeat"] == \
        (("C", 1), ("H", 2))
