"""I-165 S2b: an explicit chemistry policy per witness (domain, certificate
eligibility, recertification, the ring-edit generation filter), correction
C-1 (the depository growth unit follows the chain provider), and the
legacy/declared member parity harness (``declared_windows_parity``).

Database-free unless a test says [DB: <family>]; a database test loads that
one kinetics family through ``one_family.family_db`` (``_only_family``
evicts any other first, so one process never holds two)."""

import gc
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import one_family  # noqa: E402

DATA = HERE / "data"
needs_db = pytest.mark.skipif(not one_family.available(),
                              reason=f"needs the RMG database at {one_family.DB}")


def _only_family(name):
    """``family_db(name)`` after evicting every other cached family."""
    for other in [k for k in one_family._DB if k != name]:
        del one_family._DB[other]
    gc.collect()
    return one_family.family_db(name)


def _pp():
    from declaredChainTest import pp_mapping, spec_of
    from rmgpy.kmc import declared_chain as dc
    return dc.DeclaredChainProvider(spec_of(pp_mapping()))


def _pmma():
    from declaredChainTest import pmma_mapping, spec_of
    from rmgpy.kmc import declared_chain as dc
    return dc.DeclaredChainProvider(spec_of(pmma_mapping()))


def _pp_specs(units=4):
    from rmgpy.kmc import declared_chain as dc
    from rmgpy.kmc import windows as win
    return [win.DeclaredSpec(_pp().name, dc.ChainSpec(units, head="artificial",
                                                       tail="artificial")), ("small", "H")]


def _witness_of(specs):
    """A participants-only witness of provider builds (policy resolution
    reads nothing else)."""
    from rmgpy.kmc import boundaries as b
    from rmgpy.kmc import windows as win
    _molecules, participants = b.witness_reactants(win.builds_of(specs))
    return SimpleNamespace(participants=participants)


# --------------------------------------------------------------------------
# 1. the comparator (database-free, fabricated rows)
# --------------------------------------------------------------------------
KEY = "ref|Fam|forward"


def _row(pair, status="fixed_point", **changes):
    row = {"key": KEY, "pair": pair, "run": True, "status": status, "family": "Fam",
           "window": "w-" + pair, "view": "v-" + pair, "span": {"0": [0, 3]},
           "contexts": {"(0,'head'):declared": "fixed_point"}, "deferred": []}
    row.update(changes)
    return row


def _expected(members=2, guarded=(), deferred=0):
    return {KEY: {"members": members, "guarded": sorted(guarded), "deferred": deferred}}


GUARD = "refused:whole_species_search: Node_X label *3 candidates [7, 9]"


@pytest.fixture(scope="module")
def parity():
    import declared_windows_parity as parity
    return parity


def test_comparator_equal_rows_pass(parity):
    rows = [_row("a"), _row("b")]
    report = parity.compare(rows, [dict(r) for r in rows], _expected())
    assert report.ok and report.failures == [] and report.exit_code == 0


@pytest.mark.parametrize("side", ["legacy", "declared"])
def test_comparator_an_unpaired_member_fails(parity, side):
    both = [_row("a"), _row("b")]
    one = [_row("a")]
    legacy, declared = (both, one) if side == "legacy" else (one, both)
    report = parity.compare(legacy, declared, _expected())
    assert not report.ok and report.exit_code == 1
    assert any("without a partner" in f and "b" in f for f in report.failures)


@pytest.mark.parametrize("field,value", [("status", "refused:domain"), ("window", "other"),
                                         ("view", "other"), ("span", {"0": [0, 4]}),
                                         ("contexts", {"(0,'head'):declared": "ceiling"})])
def test_comparator_a_disagreeing_pair_fails(parity, field, value):
    report = parity.compare([_row("a"), _row("b")], [_row("a"), _row("b", **{field: value})],
                            _expected())
    assert not report.ok
    assert any("b" in f and field in f for f in report.failures), report.failures


@pytest.mark.parametrize("side", ["legacy", "declared"])
def test_comparator_a_one_sided_context_fails(parity, side):
    """A context key on one side only is a failure (r124: comparing only the
    shared keys let a side drop every common-end context)."""
    more = _row("a", contexts={"(0,'head'):declared": "fixed_point",
                               "(0,'tail'):radical": "fixed_point"})
    fewer = _row("a", contexts={"(0,'head'):declared": "fixed_point"})
    legacy, declared = (more, fewer) if side == "legacy" else (fewer, more)
    report = parity.compare([legacy, _row("b")], [declared, _row("b")], _expected())
    assert not report.ok
    assert any("(0,'tail'):radical" in f and f"{side} only" in f for f in report.failures), \
        report.failures


def test_comparator_a_one_sided_context_listed_deferred_is_not_failed(parity):
    key = "(0,'head'):vinylidene"
    legacy = _row("a", deferred=[key])
    declared = _row("a", contexts={"(0,'head'):declared": "fixed_point", key: "fixed_point"})
    report = parity.compare([legacy, _row("b")], [declared, _row("b")], _expected(deferred=1))
    assert report.ok, report.failures


def test_comparator_zero_members_compared_fails(parity):
    """Every member unrun on both sides: nothing was compared (r124)."""
    legacy = [_row("a", run=False, status="not_run"), _row("b", run=False, status="not_run")]
    declared = [dict(r) for r in legacy]
    report = parity.compare(legacy, declared, _expected())
    assert not report.ok and report.compared == 0
    assert any("zero members compared" in f for f in report.failures), report.failures


def test_comparator_contexts_requested_and_none_compared_fails(parity):
    rows = [_row("a", contexts={}), _row("b", contexts={})]
    report = parity.compare(rows, [dict(r) for r in rows], _expected())
    assert not report.ok
    assert any("zero contexts compared" in f for f in report.failures), report.failures
    unrequested = [dict(r, contexts_requested=False) for r in rows]
    report = parity.compare(unrequested, [dict(r) for r in unrequested], _expected())
    assert report.ok, report.failures


def test_comparator_a_one_sided_guard_fails(parity):
    report = parity.compare([_row("a", status=GUARD), _row("b")], [_row("a"), _row("b")],
                            _expected(guarded=["a"]))
    assert not report.ok
    assert any("guarded on one side" in f for f in report.failures)


def test_comparator_a_both_sided_guard_is_listed_not_failed(parity):
    report = parity.compare([_row("a", status=GUARD, window=None, view=None), _row("b")],
                            [_row("a", status=GUARD, window=None, view=None), _row("b")],
                            _expected(guarded=["a"]))
    assert report.ok, report.failures
    listed = report.listed["guarded"]
    assert [(g["pair"], g["family"], g["node"], g["label"], g["candidates"]) for g in listed] == \
        [("a", "Fam", "Node_X", "*3", "[7, 9]")]


def test_comparator_a_guarded_member_is_never_compared_as_equal(parity):
    """A guard on both sides is listed; it is not a passing comparison of two
    windows (the guard row's window is never read)."""
    report = parity.compare([_row("a", status=GUARD), _row("b")],
                            [_row("a", status=GUARD, window="different"), _row("b")],
                            _expected(guarded=["a"]))
    assert report.ok
    assert report.compared == 1        # only b


@pytest.mark.parametrize("members", [1, 3])
def test_comparator_a_member_count_off_expected_fails(parity, members):
    rows = [_row("a"), _row("b")]
    report = parity.compare(rows, [dict(r) for r in rows], _expected(members=members))
    assert not report.ok
    assert any("member count" in f for f in report.failures)


def test_comparator_zero_members_fails(parity):
    report = parity.compare([], [], _expected(members=0))
    assert not report.ok
    assert any("zero members" in f for f in report.failures)


def test_comparator_a_reference_missing_from_expected_fails(parity):
    rows = [_row("a"), _row("b")]
    report = parity.compare(rows, [dict(r) for r in rows], {})
    assert not report.ok


def test_comparator_a_guarded_set_off_expected_fails(parity):
    rows = [_row("a", status=GUARD), _row("b")]
    report = parity.compare(rows, [dict(r) for r in rows], _expected(guarded=["b"]))
    assert not report.ok
    assert any("guarded set" in f for f in report.failures)
    report = parity.compare(rows, [dict(r) for r in rows], _expected(guarded=[]))
    assert not report.ok


def test_comparator_a_deferred_context_is_listed_not_failed(parity):
    legacy = _row("a", deferred=["(0,'head'):vinylidene"])
    report = parity.compare([legacy, _row("b")], [_row("a"), _row("b")],
                            _expected(deferred=1))
    assert report.ok, report.failures
    assert report.listed["deferred"] == [{"key": KEY, "pair": "a",
                                          "context": "(0,'head'):vinylidene"}]


def test_comparator_rows_not_run_count_but_are_not_compared(parity):
    """``--limit``: every member is a row (the count is total); members past
    the limit are not run on either side and are not compared, and the
    guarded set is checked on the members run."""
    legacy = [_row("a"), _row("b", run=False, status="not_run")]
    declared = [_row("a"), _row("b", run=False, status="not_run")]
    report = parity.compare(legacy, declared, _expected(guarded=["b"]))
    assert report.ok, report.failures
    assert report.compared == 1


# --------------------------------------------------------------------------
# 2. the interim domain (database-free)
# --------------------------------------------------------------------------
def _mol(smiles):
    from rmgpy.molecule import Molecule
    m = Molecule().from_smiles(smiles)
    for k, atom in enumerate(m.atoms):
        atom.id = k
    return m


@pytest.mark.parametrize("name,smiles", [
    ("PP DP3 closed", "CCC(C)CC(C)CC(C)C"),
    ("PP DP3 secondary radical", "CCC(C)C[C](C)CC(C)C"),
    ("H", "[H]"),
    ("H2", "[H][H]"),
    ("PMMA DP2 closed", "COC(=O)C(C)(C)CC(C)(C)C(=O)OC"),
])
def test_interim_check_accepts(name, smiles):
    from rmgpy.kmc import domain as dom
    verdict = dom.declared_interim_check(_mol(smiles), ("C", "H", "O"))
    assert verdict.accepted, (name, verdict.reason)
    assert verdict.kind == "declared_interim"
    sites = [a.id for a in _mol(smiles).atoms if a.radical_electrons]
    assert sorted(verdict.details["radical_sites"]) == sites


def _h3():
    from rmgpy.molecule import Atom, Bond, Molecule
    atoms = [Atom(element="H", radical_electrons=0, lone_pairs=0) for _ in range(3)]
    m = Molecule(atoms=atoms)
    m.add_bond(Bond(atoms[0], atoms[1], order=1))
    atoms[2].radical_electrons = 1
    m.multiplicity = 2
    return m


def _o_one_lone_pair():
    m = _mol("CO")
    next(a for a in m.atoms if a.symbol == "O").lone_pairs = 1
    return m


def _one_site_triplet():
    m = _mol("CC[CH2]")
    m.multiplicity = 3
    return m


@pytest.mark.parametrize("name,make,reason", [
    ("carbene (u2)", lambda: _mol("[CH2]"), "u2"),
    ("charged", lambda: _mol("C[CH2+]"), "charged"),
    ("hydrogen-only H3", _h3, "hydrogen-only"),
    ("O with one lone pair", _o_one_lone_pair, "lone pair"),
    ("one-site triplet", _one_site_triplet, "multiplicity convention"),
])
def test_interim_check_refuses(name, make, reason):
    from rmgpy.kmc import domain as dom
    verdict = dom.declared_interim_check(make(), ("C", "H", "O"))
    assert not verdict.accepted, name
    assert reason in verdict.reason, (name, verdict.reason)


def test_interim_check_refuses_an_undeclared_element():
    from rmgpy.kmc import domain as dom
    verdict = dom.declared_interim_check(_mol("CCCl"))
    assert not verdict.accepted and "element" in verdict.reason
    assert not dom.declared_interim_check(_mol("CCO")).accepted       # default: C and H


def _oh_pp():
    """PP with a methyl head and a hydroxyl tail end group (O only in the
    tail end group, not in the repeat unit)."""
    from declaredChainTest import _atom, pp_mapping, spec_of
    from rmgpy.kmc import declared_chain as dc
    mapping = pp_mapping()
    methyl = {"atoms": [_atom("m", hydrogens=3)], "bonds": [],
              "ports": {"cap": {"atom": "m", "bond_order": 1}},
              "electronic_state": {"multiplicity": 1}}
    hydroxyl = {"atoms": [_atom("o", "O", hydrogens=1, lone_pairs=2)], "bonds": [],
                "ports": {"cap": {"atom": "o", "bond_order": 1}},
                "electronic_state": {"multiplicity": 1}}
    mapping["end_groups"] = {"head": methyl, "tail": hydroxyl}
    mapping["locality_policy"]["bounds"]["elements"] = ["C", "H", "O"]
    mapping["locality_policy"]["bounds"]["lone_pairs"] = [0, 2]
    return dc.DeclaredChainProvider(spec_of(mapping))


def test_interim_elements_come_from_the_witness_own_declarations(monkeypatch):
    """r124 item 1: the element set is the witness's own declarations (repeat
    unit, end groups, small species, monomer product) plus C and H; a verdict
    never depends on what else is registered."""
    import declared_parity_support as sup
    from rmgpy.kmc import boundaries as b
    from rmgpy.kmc import declared_chain as dc
    from rmgpy.kmc import windows as win
    oh = _oh_pp()
    build = oh.build("oh", dc.ChainSpec(3, head="declared", tail="declared"))
    alone = {name: p for name, p in b.PROVIDERS.items()
             if not isinstance(p, dc.DeclaredChainProvider) or p is oh}
    monkeypatch.setattr(b, "PROVIDERS", alone)
    policy = win.policy_for(sup.identity_witness(build))
    assert policy.name == "declared_interim_v1"
    verdict = policy.domain.check(build.molecule)
    assert verdict.accepted, verdict.reason
    assert policy.domain.check(_mol("CCCC(C)CC(C)CO")).accepted
    monkeypatch.undo()
    # an O-free declaration refuses an O-bearing species, whatever is registered
    _pmma()
    pp_witness = sup.identity_witness(_pp().build("w", dc.ChainSpec(3)))
    refused = win.policy_for(pp_witness).domain.check(_mol("CCCC(C)CC(C)CO"))
    assert not refused.accepted and "element" in refused.reason


def test_the_policy_name_carries_its_element_set():
    from rmgpy.kmc import windows as win
    assert win.declared_interim({"C", "H"}) is win.DECLARED_INTERIM
    assert win.DECLARED_INTERIM.name == "declared_interim_v1"
    cho = win.declared_interim({"O", "C"})
    assert cho.name == cho.domain.name == "declared_interim_v1[C,H,O]"
    assert win.declared_interim(["S", "O"]).name == "declared_interim_v1[C,H,O,S]"
    assert win.declared_interim({"O"}) is cho


def _water_pp(with_water=True):
    """PP whose declaration lists water as a small species (O nowhere else)."""
    from declaredChainTest import _atom, pp_mapping, spec_of
    from rmgpy.kmc import declared_chain as dc
    mapping = pp_mapping()
    if with_water:
        mapping["small_species"] = {"water": {
            "graph": {"atoms": [_atom("o", "O", hydrogens=2, lone_pairs=2)], "bonds": [],
                      "electronic_state": {"multiplicity": 1}},
            "disposition": "mobile"}}
        mapping["locality_policy"]["bounds"]["elements"] = ["C", "H", "O"]
        mapping["locality_policy"]["bounds"]["lone_pairs"] = [0, 2]
    return dc.DeclaredChainProvider(spec_of(mapping))


def test_a_declared_small_species_element_is_admitted():
    """r126 item 3: the small species of the witness's declaration count."""
    import declared_parity_support as sup
    from rmgpy.kmc import declared_chain as dc
    from rmgpy.kmc import windows as win
    for with_water, accepted in ((True, True), (False, False)):
        provider = _water_pp(with_water)
        witness = sup.identity_witness(provider.build("w", dc.ChainSpec(3)))
        verdict = win.policy_for(witness).domain.check(_mol("O"))
        assert verdict.accepted is accepted, (with_water, verdict.reason)


PINNED_STATES = {"H": {(0, 1)}, "C": {(0, 4)}, "Si": {(0, 4)}, "N": {(1, 3)},
                 "P": {(1, 3), (0, 5)}, "O": {(2, 2)}, "S": {(2, 2), (1, 4), (0, 6)},
                 "F": {(3, 1)}, "Cl": {(3, 1)}, "Br": {(3, 1)}, "I": {(3, 1)}}


@pytest.mark.parametrize("element", sorted(PINNED_STATES) + ["B", "Na"])
def test_neutral_state_table_is_pinned(element):
    """M-b: every listed (lone pairs, valence) state is accepted and every
    other state on the grid (0-3 lone pairs, valence 0-7 by halves) refused;
    an unlisted element has no state."""
    from rmgpy.kmc import domain as dom
    allowed = PINNED_STATES.get(element, set())
    for pairs in range(4):
        for twice in range(15):
            valence = twice / 2
            assert dom.in_neutral_state(element, pairs, valence) == \
                ((pairs, valence) in allowed), (element, pairs, valence)


def _vinylnaphthalene_provider():
    """-CH2-CH(2-naphthyl)-, declared in a Kekule form."""
    from declaredChainTest import _atom, _bond, _declaration, spec_of
    from rmgpy.kmc import declared_chain as dc
    ring = [_atom(f"r{i}", hydrogens=0 if i in (2, 5, 10) else 1) for i in range(1, 11)]
    pairs = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 10), (10, 1),
             (5, 10)]
    double = {(1, 2), (3, 4), (5, 10), (6, 7), (8, 9)}
    bonds = [_bond(f"r{x}", f"r{y}", 2 if (x, y) in double else 1) for x, y in pairs]
    mapping = _declaration([_atom("c1", hydrogens=2), _atom("c2", hydrogens=1)] + ring,
                           [_bond("c1", "c2"), _bond("c2", "r2")] + bonds, "c1", "c2")
    return dc.DeclaredChainProvider(spec_of(mapping))


def _aromatic(smiles):
    """``smiles`` with every ring bond written as a benzene bond (B, 1.5)."""
    m = _mol(smiles)
    for bond in m.get_all_edges():
        if m.is_bond_in_cycle(bond):
            bond.order = 1.5
    return m


def test_interim_check_reads_aromatic_bonds_as_rmg_does():
    """r126 item 2: a fused aromatic carbon (three B bonds) is a neutral C,
    as RMG reads it; the declared vinylnaphthalene chain is accepted after
    the production normalization (``aromatic_root``), as are naphthalene,
    pyridine, thiophene and furan written aromatic.  A ring carbon over
    valence in its aromatic form (two B bonds and two H) is refused."""
    import declared_parity_support as sup
    from rmgpy.kmc import certificates as cm
    from rmgpy.kmc import declared_chain as dc
    from rmgpy.kmc import domain as dom
    from rmgpy.kmc import windows as win
    build = _vinylnaphthalene_provider().build("n", dc.ChainSpec(2))
    root = cm.aromatic_root(build.molecule)
    assert any(bond.order == 1.5 for bond in root.get_all_edges())
    verdict = win.policy_for(sup.identity_witness(build)).domain.check(root)
    assert verdict.accepted, verdict.reason
    elements = ("C", "H", "N", "O", "S")
    for smiles in ("c1ccc2ccccc2c1", "c1ccncc1", "c1ccsc1", "c1ccoc1"):
        molecule = _aromatic(smiles)
        assert sum(bond.order == 1.5 for bond in molecule.get_all_edges()) >= 5, smiles
        verdict = dom.declared_interim_check(molecule, elements)
        assert verdict.accepted, (smiles, verdict.reason)
    over = _aromatic("c1ccccc1")
    from rmgpy.molecule import Atom, Bond
    ring_carbon = next(a for a in over.atoms if a.is_carbon())
    extra = Atom(element="H", radical_electrons=0, lone_pairs=0)
    over.add_atom(extra)
    over.add_bond(Bond(ring_carbon, extra, order=1))
    assert not dom.declared_interim_check(over, elements).accepted


def _sulfone_provider():
    """-CH2-CH(S(=O)(=O)CH3)-: a neutral hexavalent S with no lone pair."""
    from declaredChainTest import _atom, _bond, _declaration, spec_of
    from rmgpy.kmc import declared_chain as dc
    mapping = _declaration(
        [_atom("c1", hydrogens=2), _atom("c2", hydrogens=1), _atom("s", "S"),
         _atom("o1", "O", lone_pairs=2), _atom("o2", "O", lone_pairs=2),
         _atom("c3", hydrogens=3)],
        [_bond("c1", "c2"), _bond("c2", "s"), _bond("s", "o1", 2), _bond("s", "o2", 2),
         _bond("s", "c3")], "c1", "c2", elements=("C", "H", "O", "S"), lone_pairs=(0, 2))
    return dc.DeclaredChainProvider(spec_of(mapping))


def test_interim_check_accepts_a_declared_sulfone_chain():
    """r124 item 2: neutral valence states, not one lone-pair count per
    element: S6 (sulfone, 0 lone pairs) is a recognized neutral state."""
    import declared_parity_support as sup
    from rmgpy.kmc import declared_chain as dc
    from rmgpy.kmc import windows as win
    build = _sulfone_provider().build("s", dc.ChainSpec(3))
    sulfur = [a for a in build.molecule.atoms if a.symbol == "S"]
    assert sulfur and all(a.lone_pairs == 0 for a in sulfur)
    verdict = win.policy_for(sup.identity_witness(build)).domain.check(build.molecule)
    assert verdict.accepted, verdict.reason


def _divalent_s_one_lone_pair():
    m = _mol("CSC")
    next(a for a in m.atoms if a.symbol == "S").lone_pairs = 1
    return m


def test_interim_neutral_valence_states():
    """Accepted: thioether S2 (2 lone pairs), sulfoxide S4 (1), sulfone S6
    (0), phosphine P3 (1), phosphine oxide P5 (0).  Refused: a divalent S
    with one lone pair (no neutral state: 2 bonds + 1 pair), a singlet
    carbene (C with 2 bonds and a lone pair), and charges and u2 as before."""
    from rmgpy.kmc import domain as dom
    elements = ("C", "H", "O", "S", "P")
    for smiles in ("CSC", "CS(C)=O", "CS(C)(=O)=O", "CP(C)C", "CP(C)(C)=O"):
        verdict = dom.declared_interim_check(_mol(smiles), elements)
        assert verdict.accepted, (smiles, verdict.reason)
    verdict = dom.declared_interim_check(_divalent_s_one_lone_pair(), elements)
    assert not verdict.accepted and "neutral valence state" in verdict.reason
    carbene = _mol("C[C]C")
    centre = next(a for a in carbene.atoms if a.radical_electrons == 2)
    centre.radical_electrons, centre.lone_pairs = 0, 1
    carbene.multiplicity = 1
    verdict = dom.declared_interim_check(carbene, elements)
    assert not verdict.accepted and "neutral valence state" in verdict.reason
    assert "u2" in dom.declared_interim_check(_mol("[CH2]"), elements).reason
    assert "charged" in dom.declared_interim_check(_mol("C[CH2+]"), elements).reason


def test_interim_check_does_not_apply_the_ps_shape_rules():
    """Branched backbones (PP), closed or with two far-apart radical sites,
    are accepted: the interim widening is the point.  Radical sites closer
    than ``spin_verdict`` allows stay refused (the spin rule is kept)."""
    from rmgpy.kmc import domain as dom
    assert dom.declared_interim_check(_mol("C[CH]C[CH]C")).reason.startswith("spin:")
    for smiles in ("CC(C)(C)CC", "C[C](C)CCCCC[C](C)C"):
        assert dom.declared_interim_check(_mol(smiles)).accepted, smiles
        assert not dom.check(_mol(smiles)).accepted, smiles


def test_check_is_unchanged_on_the_domain_test_molecules():
    """``dom.check`` on the inputs of the existing domain tests equals the
    S2a output frozen in the red commit (``capture_domain_before.py``)."""
    import capture_domain_before as cap
    frozen = json.loads((DATA / "domain_check_before.json").read_text())
    assert cap.verdicts() == frozen


def test_universal_check_is_the_prefix_of_check():
    from rmgpy.kmc import domain as dom
    assert dom.universal_check(_mol("CC(O)C"), ("C", "H")).reason == \
        "an element other than C and H"
    assert dom.universal_check(_mol("CCCC"), ("C", "H")) is None
    assert dom.universal_check(_mol("[H]"), ("C", "H")).kind == "small:H"


# --------------------------------------------------------------------------
# 3. policy resolution (database-free)
# --------------------------------------------------------------------------
def test_policy_for_legacy_witnesses():
    from rmgpy.kmc import boundaries as b
    from rmgpy.kmc import windows as win
    for specs in ([b.PsChainSpec(3), ("small", "H")], [("alkane", 8, None, None)],
                  [("small", "H"), ("small", "H2")]):
        assert win.policy_for(_witness_of(specs)) is win.LEGACY_V1, specs


def test_policy_for_a_declared_witness():
    from rmgpy.kmc import windows as win
    assert win.policy_for(_witness_of(_pp_specs())) is win.DECLARED_INTERIM


def test_policy_for_refuses_a_mixed_witness():
    from rmgpy.kmc import boundaries as b
    from rmgpy.kmc import windows as win
    with pytest.raises(win.WindowError, match="policy"):
        win.policy_for(_witness_of([b.PsChainSpec(2), _pp_specs()[0]]))


def test_policy_objects():
    from rmgpy.kmc import domain as dom
    from rmgpy.kmc import windows as win
    assert (win.LEGACY_V1.name, win.LEGACY_V1.ring_edits) == ("legacy_v1", "exclude")
    from rmgpy.kmc import certificates as cm
    assert win.LEGACY_V1.name == win.LEGACY_V1.domain.name == cm.LEGACY_DOMAIN_POLICY
    assert (win.DECLARED_INTERIM.name, win.DECLARED_INTERIM.ring_edits) == \
        ("declared_interim_v1", "exclude")
    m = _mol("CC(C)(C)CC")
    assert win.LEGACY_V1.domain.check(m) == dom.check(m)
    assert win.DECLARED_INTERIM.domain.check(m) == dom.declared_interim_check(m)
    assert win.DECLARED_INTERIM.domain.spin_verdict(m) == dom.spin_verdict(m)
    with pytest.raises(Exception):
        win.LEGACY_V1.name = "x"


def test_an_explicit_policy_wins(monkeypatch):
    """``species_domain_refusal(witness, policy=)``: the explicit policy is
    used, never the resolved one."""
    from rmgpy.kmc import windows as win
    calls = []

    def spy(molecule):
        calls.append(molecule)
        return win.LEGACY_V1.domain.check(molecule)
    policy = win.ChemistryPolicy("spy", win.PolicyDomain("spy", spy, win.LEGACY_V1.domain
                                                         .spin_verdict), "exclude")
    import declared_parity_support as sup
    from rmgpy.kmc import boundaries as b
    witness = sup.identity_witness(win.builds_of([b.PsChainSpec(2)])[0])
    assert win.species_domain_refusal(witness, policy=policy) is None
    assert calls
    resolved = win.species_domain_refusal(witness)
    assert resolved is None and len(calls) == 2          # one species, both sides; unchanged


def _two_molecule_guard_certificate():
    """A certificate on two molecules (atoms 0-1 and 5) whose one passed
    predicate reads molecule 1 (anchor atom 5) with a carbon minimum it fails."""
    predicate = {"node": "C_rad", "result": "passed", "anchors": {"*3": 5},
                 "guards": {"element_count": {"read": True, "group_minimum": {"C": 1},
                                              "structure": {"H": 1}}}}
    mapping = {"labels": [[1, "*3", 5]], "descent": "d"}
    return {"channels": [{"kind": "kept", "members": [{"mappings": [mapping]}]}],
            "descents": {"d": [{"root": predicate}]}, "global_guards": [],
            "witness": {"molecules": [[0, 1], [5]]}}, {frozenset({("*3", 5)})}


def test_an_element_count_guard_on_a_molecule_no_extension_grows_is_invariant():
    """The guard reads the small participant (an H atom): a chain's step
    delta cannot add carbon to it, so the failed minimum is extension
    invariant, not 'an extension could satisfy it'.  Without ``growable``
    the old reading (every molecule takes every provider's delta) stays."""
    from rmgpy.kmc import windows as win
    facts = {"declared:x": {"backbone_head": {"closed_shell": True, "delta": {"C": 3, "H": 6}}}}
    cert, wanted = _two_molecule_guard_certificate()
    old = win.certificate_support(cert, wanted, fragment_facts=facts)
    assert old.refusal == "refused:guard: element_count C"
    still = win.certificate_support(cert, wanted, fragment_facts=facts,
                                    growable=frozenset({0, 1}))
    assert still.refusal is None
    assert {k: v for k, v, _r in still.dispositions}[("guard", "element_count:C")] == \
        "invariant_failed"
    grows = win.certificate_support(cert, wanted, fragment_facts=facts,
                                    growable=frozenset({0, 1, 5}))
    assert grows.refusal == "refused:guard: element_count C"


# --------------------------------------------------------------------------
# 4-7. [DB: H_Abstraction] eligibility, recertification, PP + H, ring edits
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def habs():
    if not one_family.available():
        pytest.skip(f"needs the RMG database at {one_family.DB}")
    from rmgpy.kmc import certificates as cm
    db = _only_family("H_Abstraction")
    fams = db.kinetics.families
    return fams["H_Abstraction"], cm.Certifier(fams, max_channel_sets=16)


@pytest.fixture(scope="module")
def pp_cases(habs):
    from rmgpy.kmc import windows as win
    family, certifier = habs
    stats = {}
    cases = win.member_cases(certifier, family, win.builds_of(_pp_specs()), "forward",
                             label="PP4+H", stats=stats)
    return cases, stats


def _pp_certificate(habs, case, **kwargs):
    from rmgpy.kmc import windows as win
    family, certifier = habs
    return certifier.certify(family, win.template_molecules(case.witness, "forward"),
                             member_labels=case.wanted, **kwargs)


@needs_db
def test_eligibility_follows_the_declared_policy(habs, pp_cases):
    from rmgpy.kmc import windows as win
    case = pp_cases[0][0]
    cert = _pp_certificate(habs, case, domain=win.DECLARED_INTERIM.domain)
    assert cert["eligible"] is True
    assert cert["domain_policy"] == "declared_interim_v1"
    legacy = _pp_certificate(habs, case)
    assert legacy["eligible"] is False and "domain_policy" not in legacy
    assert {(ch.get("eligibility") or {}).get("disposition") for ch in legacy["channels"]} == \
        {"refused:domain"}


@needs_db
def test_a_legacy_certificate_without_a_domain_is_byte_identical(habs):
    import capture_c1_before as cap
    family, certifier = habs
    frozen = json.loads((DATA / "legacy_ps_dp3_h_certificate.json").read_text())
    live = json.loads(cap.canonical(cap.ps_dp3_h_certificate(certifier, family)))
    assert live["member"] == frozen["member"]
    assert cap.canonical(live["certificate"]) == cap.canonical(frozen["certificate"])


def _longer(case):
    """The member's witness one period longer at its first artificial port,
    and the identity map of the surviving atoms."""
    from rmgpy.kmc import boundaries as b
    from rmgpy.kmc import windows as win
    port = next(p for p in b.vector_ports(case.witness) if case.witness.port(p).artificial)
    longer, _record = b.extend(case.witness, port)
    return win.template_molecules(longer, "forward")


@needs_db
def test_recertification_follows_the_policy(habs, pp_cases):
    from rmgpy.kmc import windows as win
    family, certifier = habs
    for case in pp_cases[0]:
        cert = _pp_certificate(habs, case, domain=win.DECLARED_INTERIM.domain)
        mols = _longer(case)
        present = {a.id for m in mols for a in m.atoms}
        amap = {row[0]: row[0] for row in cert["witness"]["atoms"] if row[0] in present}
        same = certifier.recertify_certificate(family, cert, mols, amap,
                                               domain=win.DECLARED_INTERIM.domain)
        if same.get("same_selection"):
            break
    else:
        pytest.fail("no PP member recertifies with the same selection one period longer")
    assert not same["refusals"]
    changed = certifier.recertify_certificate(family, cert, mols, amap)
    assert not changed.get("same_selection")
    assert any("refused:domain_policy_changed" in r for r in changed["refusals"])


@pytest.fixture(scope="module")
def pp_runs(habs, pp_cases):
    """Every PP + H member through ``run_member(contexts=True)`` with the
    provider's span bound (the declared policy is resolved)."""
    from rmgpy.kmc import windows as win
    family, certifier = habs
    bound = win.span_bound(win.family_constants(family), _pp())
    out = []
    for case in pp_cases[0]:
        out.append((case, win.run_member(certifier, family, case, bound=bound, contexts=True,
                                         specs=_pp_specs())))
    return bound, out


@needs_db
def test_pp_plus_h_is_not_vacuous(habs, pp_cases, pp_runs):
    from rmgpy.kmc import windows as win
    cases, stats = pp_cases
    expected = json.loads((DATA / "declared_parity_expected.json").read_text())
    assert len(cases) > 0
    assert len(cases) == expected["pp4+h|H_Abstraction|forward"]["members"]
    assert stats["ring_change_mappings"] == 0
    bound, runs = pp_runs
    statuses = [r["status"] for _c, r in runs]
    assert not [s for s in statuses if s.startswith("refused:domain")], statuses
    assert not [s for s in statuses if "provider" in s or "roster" in s or "legacy_spec" in s]
    extended = [(c, r) for c, r in runs if r["status"] == "fixed_point"
                and any(d > 0 for _p, d in r["depths"])]
    assert extended, statuses
    case, res = extended[0]
    assert res["span_bound_held"]
    reached = [ctx for _c, r in runs for ctx in r.get("contexts") or ()
               if any(end in ("declared", "radical") for _p, end in ctx["assignment"])
               and not ctx["status"].startswith("refused")]
    assert reached, [r.get("contexts") for _c, r in runs][:3]
    assert all(r.get("seconds") is not None for _c, r in runs)


@needs_db
def test_pp_plus_h_runs_under_the_declared_policy(habs, pp_cases):
    """The certificate path sees the declared policy: a run_member forced to
    LEGACY_V1 is refused by the domain, the resolved one is not."""
    from rmgpy.kmc import windows as win
    family, certifier = habs
    case = pp_cases[0][0]
    legacy = win.run_member(certifier, family, case, contexts=False, policy=win.LEGACY_V1)
    assert legacy["status"] == "refused:domain"


@needs_db
def test_ring_edit_policy(habs):
    from rmgpy.kmc import boundaries as b
    from rmgpy.kmc import windows as win
    family, certifier = habs
    builds = win.builds_of([b.PsChainSpec(3), ("small", "H")])

    def labels(cases):
        return sorted(c.label for c in cases)
    default = win.member_cases(certifier, family, builds, "forward")
    assert labels(win.member_cases(certifier, family, builds, "forward",
                                   policy=win.LEGACY_V1)) == labels(default)
    admit = win.ChemistryPolicy("admit_test", win.LEGACY_V1.domain, "admit")
    stats = {}
    admitted = win.member_cases(certifier, family, builds, "forward", policy=admit)
    every = win.member_cases(certifier, family, builds, "forward", skip_ring_changes=False)
    assert labels(admitted) == labels(every)
    excluded = win.member_cases(certifier, family, builds, "forward", stats=stats)
    assert stats["ring_change_mappings"] > 0
    assert sum(len(c.moves) for c in admitted) > sum(len(c.moves) for c in excluded)


def test_a_changed_element_set_is_a_changed_policy(habs, pp_cases):
    """r126 item 1: certified under the C/H interim domain, recertified under
    C/H/O: the element set is part of the policy identity, so the gate
    refuses (it was accepted with same_selection when both were named
    ``declared_interim_v1``)."""
    from rmgpy.kmc import windows as win
    family, certifier = habs
    cho = win.declared_interim({"O"})
    for case in pp_cases[0]:
        cert = _pp_certificate(habs, case, domain=win.DECLARED_INTERIM.domain)
        mols = _longer(case)
        present = {a.id for m in mols for a in m.atoms}
        amap = {row[0]: row[0] for row in cert["witness"]["atoms"] if row[0] in present}
        same = certifier.recertify_certificate(family, cert, mols, amap,
                                               domain=win.DECLARED_INTERIM.domain)
        if same.get("same_selection"):
            break
    else:
        pytest.fail("no PP member recertifies with the same selection one period longer")
    changed = certifier.recertify_certificate(family, cert, mols, amap, domain=cho.domain)
    assert not changed.get("same_selection")
    assert any("refused:domain_policy_changed" in r and "declared_interim_v1[C,H,O]" in r
               for r in changed["refusals"]), changed["refusals"]
    stamped = _pp_certificate(habs, case, domain=cho.domain)
    assert stamped["domain_policy"] == "declared_interim_v1[C,H,O]"


def test_closing_partners_are_not_training_species(habs):
    """Item 4 (Zeus 4801839): once C-3 let H and H2 close, ``closing_partners``
    picked them; a depository match reads a training species whole, so a
    class closed against one (1-phenylethyl + H2) reaches no fixed point and
    the reduced plan's closings lose its leaf.  A closing partner is never a
    species of the family's training reactions."""
    from rmgpy.kmc import channels as chm
    from rmgpy.kmc import windows as win
    family, certifier = habs
    partners = win.closing_partners(certifier, family, win.family_constants(family))
    training, species = win.training_keys(family), chm.SpeciesKeys()
    assert set(partners["forward"]) == {0, 1}
    for slot, (name, spec) in partners["forward"].items():
        assert species(win.builds_of([spec])[0].molecule) not in training, (slot, name)


def test_cover_cases_are_keyed_by_policy(habs):
    """M-a: the same family, direction and specs under two policies (ring
    edits excluded, then admitted) never share cached cover cases."""
    from rmgpy.kmc import boundaries as b
    from rmgpy.kmc import windows as win
    family, certifier = habs
    specs = (b.PsChainSpec(3), ("small", "H"))
    admit = win.ChemistryPolicy("admit_test", win.LEGACY_V1.domain, "admit")
    win._COVER_CASES.clear()
    excluded = win._cover_cases(certifier, family, specs, "forward", win.LEGACY_V1)
    admitted = win._cover_cases(certifier, family, specs, "forward", admit)
    builds = win.builds_of(specs)
    every = win.member_cases(certifier, family, builds, "forward", label="cover",
                             skip_ring_changes=False)
    assert sorted(c.label for c in admitted) == sorted(c.label for c in every)
    assert sum(len(c.moves) for c in admitted) > sum(len(c.moves) for c in excluded)
    assert win._cover_cases(certifier, family, specs, "forward", win.LEGACY_V1) is excluded


# --------------------------------------------------------------------------
# 8-9. [DB: R_Recombination] C-1 and local parity
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def recombination():
    if not one_family.available():
        pytest.skip(f"needs the RMG database at {one_family.DB}")
    return _only_family("R_Recombination")


@pytest.fixture(scope="module")
def c16_parity(recombination, parity):
    return parity.run_reference("c16", "R_Recombination", "reverse", limit=10)


@needs_db
def test_c1_before_is_frozen_and_after_is_live(c16_parity):
    frozen = json.loads((DATA / "c1_before.json").read_text())["members"]
    rows = c16_parity["corrections"]["C-1"]["rows"]
    assert rows
    for row in rows:
        assert row["before"] == frozen[row["member"]], row["member"]
        assert all(ch["repeat_unit"] == {"C": 1, "H": 2} for ch in row["after"]["channels"])
        assert all(ch["repeat_unit"] == {"C": 8, "H": 8} for ch in row["before"]["channels"])
    changed = [r for r in rows if [c["verdict"] for c in r["before"]["channels"]] !=
               [c["verdict"] for c in r["after"]["channels"]]]
    summary = c16_parity["corrections"]["C-1"]["summary"]
    if not changed:
        assert "no depository_exclusion verdict changes" in summary
        assert "before" in summary and "after" in summary


@needs_db
@pytest.mark.parametrize("reference,limit", [("c16", 10), ("dp8_closed", 20)])
def test_local_parity(recombination, parity, c16_parity, reference, limit):
    result = c16_parity if reference == "c16" else \
        parity.run_reference(reference, "R_Recombination", "reverse", limit=limit)
    report = result["report"]
    assert report["failures"] == [], report["failures"]
    key = f"{reference}|R_Recombination|reverse"
    expected = json.loads((DATA / "declared_parity_expected.json").read_text())[key]
    assert result["members"]["legacy"] == result["members"]["declared"] == \
        expected["members"] > 0
