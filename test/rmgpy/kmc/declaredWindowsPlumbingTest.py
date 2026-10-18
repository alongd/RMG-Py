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
from declaredChainTest import cyclopentane_mapping, pmma_mapping, pp_mapping, spec_of

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
    assert _observe(_ps_witness()).get("repeat") is None
    toy = sup.identity_witness(b.PROVIDERS["linear_alkane"].build(8))
    assert _observe(toy).get("repeat") is None
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


def test_portless_monomer_without_an_embedding_uses_the_heavy_diameter():
    # an ethylene "monomer" cannot carry the three-carbon PP period
    mapping = pp_mapping()
    graph = _propylene()
    graph["atoms"] = [dict(graph["atoms"][0]), dict(graph["atoms"][1], hydrogens=2)]
    graph["bonds"] = graph["bonds"][:1]
    mapping["monomer_product"] = graph
    provider = dc.DeclaredChainProvider(spec_of(mapping))
    w = sup.identity_witness(provider.monomer())
    rows, heavy = _heavy_of(w)
    assert len(heavy) == 2
    assert win.window_span(w, set(rows))[0]["backbone_atoms"] == 2


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


def test_legacy_recertification_keeps_the_default_growth_unit(monkeypatch):
    spy = _observe_twice(_ps_witness(), monkeypatch)
    assert spy.certify_calls[0].get("repeat") is None
    assert spy.recertify_calls[0].get("repeat") is None
