"""Tests for rmgpy/kmc/rmg_moves.py: RMG's own directed, atom-mapped moves
before any degeneracy collapse, and their agreement with the window matcher.

Fast part: small molecules, homolysis (C16 and closed PS DP8), the DP3
positive control, policies, determinism, the matcher view-cache leak.
Slow part (``RMG_KMC_SLOW=1``): the S2 acceptance states against the S2
hand-written template sets (``rmg_moves_harness``).
"""

import ast
import gc
import hashlib
import json
import os
import subprocess
import sys
import time
import weakref
from collections import Counter
from pathlib import Path

import pytest

from rmgpy.kmc import rmg_moves as rm
from rmgpy.molecule import Molecule

import rmg_moves_harness as harness

TEST_DIR = Path(__file__).resolve().parent
ROOT = TEST_DIR.parents[2]
SLOW = os.environ.get("RMG_KMC_SLOW") == "1"

pytestmark = pytest.mark.skipif(not harness.DATABASE.is_dir(),
                                reason=f"needs the RMG database at {harness.DATABASE}")


@pytest.fixture(scope="module")
def fam():
    return harness.families()


def mol(smiles, start):
    """A Molecule from SMILES with atom ids start, start + 1, ... in RMG atom order."""
    m = Molecule().from_smiles(smiles)
    for k, atom in enumerate(m.atoms):
        atom.id = start + k
    return m


def kept(moves):
    return [m for m in moves if m.status == "kept"]


def h_ids(m, heavy):
    atom = next(a for a in m.atoms if a.id == heavy)
    return sorted(o.id for o in atom.bonds if o.is_hydrogen())


# --------------------------------------------------------------------------
# one move per mapping, atom ids, resonance
# --------------------------------------------------------------------------
def test_one_move_per_mapping_no_product_collapse(fam):
    """CH3* + C2H6: six C-H bonds, six directed moves, six mapped products
    (RMG's degeneracy would fold them into one reaction of degeneracy 6)."""
    methyl, ethane = mol("[CH3]", 0), mol("CC", 10)
    moves = rm.generate_moves(fam["H_Abstraction"], [methyl, ethane])
    assert rm.status_counts(moves)["kept"] == len(moves) == 6
    donors = sorted((m.atom("*1"), m.atom("*2")) for m in moves)
    expected = sorted((c, h) for c in (10, 11) for h in h_ids(ethane, c))
    assert donors == expected
    assert {m.atom("*3") for m in moves} == {0}
    assert {m.order for m in moves} == {(1, 0)}  # ethane in the X_H slot
    assert len({rm.product_key(m) for m in moves}) == 6
    for m in moves:
        rows = m.product_rows()
        assert sorted(rows) == sorted(a.id for x in (methyl, ethane) for a in x.atoms)
        assert (m.atom("*2"), 1.0) in rows[0][5]  # the H now sits on the methyl carbon


def test_atom_ids_survive_relabel(fam):
    """Benzyl + ethane: four reactive benzyl forms (aromatic + three quinoid),
    one acceptor site each; all ids survive the merge, the recipe and RMG's
    relabelling (*1 <-> *3 on the products), which is undone."""
    benzyl, ethane = mol("[CH2]c1ccccc1", 0), mol("CC", 100)
    forms = rm.reactant_forms(benzyl)
    assert len(forms) == 4 and all(f.reactive for f in forms)
    assert all(sorted(a.id for a in f.atoms) == sorted(a.id for a in benzyl.atoms) for f in forms)
    moves = kept(rm.generate_moves(fam["H_Abstraction"], [benzyl, ethane]))
    assert len(moves) == 24
    assert Counter(m.atom("*3") for m in moves) == {0: 6, 2: 6, 4: 6, 6: 6}  # CH2, o, p, o'
    for m in moves:
        plabels = dict(m.product_labels)
        assert plabels["*1"] == m.atom("*3") and plabels["*3"] == m.atom("*1")
        assert plabels["*2"] == m.atom("*2")
        rows = m.product_rows()
        assert sorted(rows) == sorted(a.id for x in (benzyl, ethane) for a in x.atoms)
        assert (m.atom("*2"), 1.0) in rows[m.atom("*3")][5]
        assert all(b != m.atom("*2") for b, _o in rows[m.atom("*1")][5])


def test_same_move_from_two_resonance_forms_has_the_same_atom_ids(fam):
    """CH3* + pent-1-en-3-yl (two forms: radical on C3 or on C1).  Abstraction
    of a terminal methyl H of the donor is found from BOTH donor forms, with
    identical reacting-atom ids and the same product modulo resonance."""
    methyl, donor = mol("[CH3]", 0), mol("CC[CH]C=C", 10)
    assert len(rm.reactant_forms(donor)) == 2
    moves = kept(rm.generate_moves(fam["H_Abstraction"], [methyl, donor]))
    remote = [m for m in moves if m.atom("*1") == 10]
    by_labels = {}
    for m in remote:
        by_labels.setdefault(m.labels, []).append(m)
    assert len(by_labels) == 3
    for labels, same in by_labels.items():
        assert len(same) == 2 and {m.forms for m in same} == {(0, 0), (0, 1)}
        assert same[0].product_rows() != same[1].product_rows()
        assert rm.product_key(same[0]) == rm.product_key(same[1])
        assert rm.same_product(same[0].product_rows(), same[1].product_rows())


def test_reactant_resonance_policy(fam):
    benzyl, ethane = mol("[CH2]c1ccccc1", 0), mol("CC", 100)
    given = rm.generate_moves(fam["H_Abstraction"], [benzyl, ethane], resonance="given")
    assert len(given) == 6 and {m.atom("*3") for m in given} == {0}
    assert {m.forms for m in given} == {(0, 0)}
    every = rm.generate_moves(fam["H_Abstraction"], [benzyl, ethane])
    assert {m.atom("*3") for m in every} == {0, 2, 4, 6}
    with pytest.raises(rm.RmgMovesError):
        rm.generate_moves(fam["H_Abstraction"], [benzyl, ethane], resonance="reactive")
    # the enumerated forms do not depend on the form the reactant is written in
    import windowMatcherTest as wmt
    ring = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]
    kekule = wmt.with_h([("C", 1)] + [("C", 0)] * 6,
                        [(0, 1, 1)] + [(a, b, 2 if k % 2 == 0 else 1)
                                       for k, (a, b) in enumerate(ring)],
                        [2, 0, 1, 1, 1, 1, 1])
    seen = []
    for spec in [wmt.benzyl(), wmt.benzyl("ortho"), wmt.benzyl("para"), kekule]:
        (start_mol,) = rm.molecules_from_state(wmt.molecules(spec))
        seen.append([rm.molecule_rows(f) for f in rm.reactant_forms(start_mol)])
    assert len(seen[0]) == 4 and all(x == seen[0] for x in seen)


# --------------------------------------------------------------------------
# filter accounting
# --------------------------------------------------------------------------
def test_filters_are_statuses_not_drops(fam):
    """CH3* + CH3CH2*: 8 mappings.  Kept: the 3 methyl H of CH3* to the ethyl
    radical centre and the 2 alpha H of Et* to CH3*; forbidden: the 3 beta H of
    Et* (H on a carbon next to the radical: disprop1_base_case, which is
    disproportionation); forbidden moves still carry the product.
    CH3* + CH4: 4 identity moves (unchanged_species).  Cyclohexane bond
    dissociation: 12 C-H cuts x 2 orientations kept, 6 ring C-C x 2 refused by
    the recipe (one product where two are expected)."""
    moves = rm.generate_moves(fam["H_Abstraction"], [mol("[CH3]", 0), mol("C[CH2]", 10)])
    counts = rm.status_counts(moves)
    assert (counts["kept"], counts["forbidden"], len(moves)) == (5, 3, 8)
    for m in moves:
        if m.status == "forbidden":
            assert "disprop1_base_case" in m.reason and m.products is not None
            assert m.template is None
        else:
            assert m.template and m.reason == ""
    identity = rm.generate_moves(fam["H_Abstraction"], [mol("[CH3]", 0), mol("C", 10)])
    assert rm.status_counts(identity)["unchanged_species"] == len(identity) == 4
    ring = rm.generate_moves(fam["R_Recombination"], [mol("C1CCCCC1", 0)], "reverse")
    counts = rm.status_counts(ring)
    assert (counts["kept"], counts["recipe_failed"], len(ring)) == (24, 12, 36)
    for m in ring:
        cut = {a for _s, _l, a in m.labels}
        is_cc = all(a < 6 for a in cut)
        assert (m.status == "recipe_failed") == is_cc
    sel = rm.generate_moves(fam["H_Abstraction"], [mol("[CH3]", 0), mol("C[CH2]", 10)],
                            select=lambda c: c.order == (0, 1))
    assert len(sel) == 8 and all(m.products is None for m in sel if m.status == "excluded")
    assert {m.order for m in sel if m.status == "excluded"} == {(1, 0)}
    assert rm.status_counts(sel)["constraint"] == 0


def test_never_calls_degeneracy_or_constraint_machinery(fam, monkeypatch):
    import rmgpy.data.kinetics.common as common
    import rmgpy.data.kinetics.family as family_module

    def boom(*_a, **_k):
        raise AssertionError("must not be called to produce moves")

    for target, name in [(common, "find_degenerate_reactions"),
                         (family_module, "find_degenerate_reactions"),
                         (family_module, "fails_species_constraints"),
                         (family_module.KineticsFamily, "calculate_degeneracy"),
                         (family_module.KineticsFamily, "_generate_reactions"),
                         (family_module.KineticsFamily, "generate_reactions"),
                         (common, "reduce_same_reactant_degeneracy")]:
        monkeypatch.setattr(target, name, boom)
    assert rm.generate_moves(fam["Disproportionation"], [mol("C[CH2]", 0), mol("C[CH2]", 10)])
    assert rm.generate_moves(fam["R_Recombination"], [mol("CCCC", 0)], "reverse")
    assert rm.generate_moves(fam["intra_H_migration"], [mol("[CH2]CCCC", 0)])
    source = (ROOT / "rmgpy" / "kmc" / "rmg_moves.py").read_text()
    body = source.split("def rmg_degeneracy")[0]
    assert "find_degenerate_reactions(" not in body and "calculate_degeneracy(" not in body


# --------------------------------------------------------------------------
# A + A and A + B
# --------------------------------------------------------------------------
def test_identical_radicals_recombine_in_both_orders(fam):
    """CH3* + CH3* on two objects: RMG tries both template slot orders, so the
    one pair appears as (0, 1) and (1, 0) with the same mapped product."""
    a, b = mol("[CH3]", 0), mol("[CH3]", 10)
    moves = rm.generate_moves(fam["R_Recombination"], [a, b])
    assert [m.order for m in moves] == [(0, 1), (1, 0)]
    assert [m.labels for m in moves] == [((0, "*", 0), (1, "*", 10)), ((0, "*", 10), (1, "*", 0))]
    assert rm.product_key(moves[0]) == rm.product_key(moves[1])
    assert all(m.status == "kept" for m in moves)
    ab = rm.generate_moves(fam["R_Recombination"], [a, mol("C[CH2]", 10)])
    assert [m.order for m in ab] == [(0, 1), (1, 0)]
    with pytest.raises(rm.RmgMovesError):
        rm.generate_moves(fam["R_Recombination"], [a, a])
    with pytest.raises(rm.AtomIdError):
        rm.generate_moves(fam["R_Recombination"], [a, mol("[CH3]", 0)])
    bare = Molecule().from_smiles("[CH3]")
    with pytest.raises(rm.AtomIdError):
        rm.generate_moves(fam["R_Recombination"], [a, bare])


def test_disproportionation_ordered_acceptor_donor(fam):
    """Et* + Et*: 3 moves with A the acceptor (B gives one of its 3 methyl H),
    3 with B the acceptor; the two orders give different mapped products.
    CH3* + Et*: only CH3* can accept (it has no beta H): 3 moves."""
    a, b = mol("C[CH2]", 0), mol("C[CH2]", 10)
    moves = rm.generate_moves(fam["Disproportionation"], [a, b])
    assert len(kept(moves)) == len(moves) == 6
    by_order = Counter(m.order for m in moves)
    assert by_order == {(0, 1): 3, (1, 0): 3}
    acceptors = {m.order: m.atom("*1") for m in moves}
    assert acceptors == {(0, 1): 1, (1, 0): 11}
    keys = {o: {rm.product_key(m) for m in moves if m.order == o} for o in by_order}
    assert len(keys[(0, 1)]) == 3 and len(keys[(1, 0)]) == 3
    assert not keys[(0, 1)] & keys[(1, 0)]
    mixed = rm.generate_moves(fam["Disproportionation"], [mol("[CH3]", 0), b])
    assert len(kept(mixed)) == 3 and {m.atom("*1") for m in kept(mixed)} == {0}


# --------------------------------------------------------------------------
# the bridge to the matcher's identity
# --------------------------------------------------------------------------
def test_product_key_equals_the_matcher_identity_on_benzyl_ethane(fam):
    import windowMatcherTest as wmt
    from rmgpy.kmc import window_matcher as wm

    state = wmt.molecules(wmt.benzyl(), wmt.ethane())
    result = wm.match([wmt.benzyl_acceptor_ethane_donor()], state)
    matcher = {}
    for tr in result.transitions:
        matcher[rm.product_key(wm.product_rows(state, tr))] = wm.product_rows(state, tr)
    assert len(matcher) == 6
    mols = rm.molecules_from_state(state)
    moves = [m for m in kept(rm.generate_moves(fam["H_Abstraction"], list(mols)))
             if m.atom("*3") == 0]
    assert len(moves) == 6  # one reactive form puts the radical on CH2
    assert {rm.product_key(m) for m in moves} == set(matcher)
    for m in moves:
        assert rm.same_product(m.product_rows(), matcher[rm.product_key(m)])
    # the key does not depend on the form the reactant state was written in
    quinoid = wmt.molecules(wmt.benzyl("para"), wmt.ethane())
    qmoves = [m for m in kept(rm.generate_moves(fam["H_Abstraction"],
                                                list(rm.molecules_from_state(quinoid))))
              if m.atom("*3") == 0]
    assert {rm.product_key(m) for m in qmoves} == set(matcher)


def test_same_product_is_strict_and_resonance_aware():
    """Ethylene and the 1,2-ethanediyl biradical share every atom and bond
    (same sigma graph) but not the electrons, and RMG resonance does not join
    them: different products (a non-strict isomorphism would merge them).
    The two allyl forms: one product."""
    ethylene = rm.molecule_rows(mol("C=C", 0))
    biradical = {a: (r[0], r[1], r[2], 1 if r[0] == "C" else 0, r[4],
                     tuple((b, 1.0) for b, _o in r[5])) for a, r in ethylene.items()}
    assert sorted(biradical) == sorted(ethylene)
    assert not rm.same_product(ethylene, biradical)
    assert rm.product_key(ethylene) != rm.product_key(biradical)
    # charge allocation on non-conjugated atoms: two products
    zwitter = rm.molecule_rows(mol("[CH2+]CC[CH2-]", 0))
    ends = [a for a, r in zwitter.items() if r[2] != 0]
    swapped = {a: (r[0], r[1], -r[2], r[3], 1 - r[4], r[5]) if a in ends else r
               for a, r in zwitter.items()}
    assert not rm.same_product(zwitter, swapped)
    assert rm.product_key(zwitter) != rm.product_key(swapped)
    allyl = mol("[CH2]C=C", 0)
    forms = rm.reactant_forms(allyl)
    assert len(forms) == 2
    rows = [rm.molecule_rows(f) for f in forms]
    assert rows[0] != rows[1]
    assert rm.same_product(rows[0], rows[1]) and rm.product_key(rows[0]) == rm.product_key(rows[1])


def _benzyl_start_rows():
    """Atom-id rows of benzyl written in four starting forms (aromatic,
    ortho- and para-quinoid, Kekule benzylic)."""
    import windowMatcherTest as wmt
    ring = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]
    kekule = wmt.with_h([("C", 1)] + [("C", 0)] * 6,
                        [(0, 1, 1)] + [(a, b, 2 if k % 2 == 0 else 1)
                                       for k, (a, b) in enumerate(ring)],
                        [2, 0, 1, 1, 1, 1, 1])
    out = []
    for spec in [wmt.benzyl(), wmt.benzyl("ortho"), wmt.benzyl("para"), kekule]:
        (m,) = rm.molecules_from_state(wmt.molecules(spec))
        out.append(rm.molecule_rows(m))
    return out


def _all_form_rows(start_rows):
    """Every starting rows mapping plus the rows of every reactive RMG form."""
    rows = list(start_rows)
    for start in start_rows:
        rows.extend(rm.molecule_rows(f) for f in rm.reactant_forms(rm.molecule_from_rows(start)))
    return rows


def _check_form_independent(rows):
    distinct = {rm._signature(r) for r in rows}
    assert len(distinct) >= 4  # the forms really are written differently
    fresh = {rm.ProductKeys(cache=False)(r) for r in rows}
    assert len(fresh) == 1
    forward, backward = rm.ProductKeys(cache=True), rm.ProductKeys(cache=True)
    keys = {forward(r) for r in rows} | {backward(r) for r in reversed(rows)}
    assert keys == fresh
    return fresh.pop()


def test_product_key_is_independent_of_the_form_the_rows_arrive_in():
    """No generation.  Benzyl from four starting forms, and a mid-chain
    1-phenylethyl-type radical (CC[C](CC)c1ccccc1) in atom-id rows: every
    reactive RMG form gives one key, with fresh uncached ProductKeys and with
    cached ones fed in two orders.  Moving the radical to the beta carbon of
    the same heavy skeleton (not a resonance move) changes the key."""
    _check_form_independent(_all_form_rows(_benzyl_start_rows()))
    chain = rm.molecule_rows(mol("CC[C](CC)c1ccccc1", 0))
    chain_forms = [rm.molecule_rows(f) for f in rm.reactant_forms(rm.molecule_from_rows(chain))]
    assert len(chain_forms) == 4
    key = _check_form_independent(_all_form_rows([chain] + chain_forms))
    beta = rm.molecule_rows(mol("C[CH]C(CC)c1ccccc1", 0))
    heavy = lambda rows: sorted((a, r[0], tuple(b for b, _o in r[5] if rows[b][0] != "H"))
                                for a, r in rows.items() if r[0] != "H")
    assert heavy(beta) == heavy(chain)
    assert rm.ProductKeys(cache=False)(beta) != key
    assert rm.product_key(beta) != rm.product_key(chain)


def test_converter_uses_rmg_resonance_only():
    """No import of the matcher or its oracle (static and at run time), and a
    nine-ring conjugated radical the matcher refuses is keyed by RMG."""
    tree = ast.parse((ROOT / "rmgpy" / "kmc" / "rmg_moves.py").read_text())
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported.add(node.module or "")
            imported.update(f"{node.module}.{a.name}" for a in node.names)
    assert not {n for n in imported if "window_matcher" in n or "transition_oracle" in n}
    code = (
        "import sys\n"
        "from rmgpy.molecule import Molecule\n"
        "from rmgpy.kmc import rmg_moves as rm\n"
        "m = Molecule().from_smiles('[CH]1C=CC=CC=CC=C1')\n"
        "for k, a in enumerate(m.atoms): a.id = k\n"
        "key = rm.product_key(rm.molecule_rows(m))\n"
        "assert key\n"
        "bad = [n for n in sys.modules if 'window_matcher' in n or 'transition_oracle' in n]\n"
        "assert not bad, bad\n"
        "print('OK')\n")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                         cwd=str(ROOT), timeout=300)
    assert out.stdout.strip().endswith("OK"), out.stderr[-3000:]
    from rmgpy.kmc import window_matcher as wm
    import windowMatcherTest as wmt
    ring = wmt.molecules(wmt.with_h([("C", 1)] + [("C", 0)] * 8,
                                    [(0, 1, 1), (1, 2, 2), (2, 3, 1), (3, 4, 2), (4, 5, 1),
                                     (5, 6, 2), (6, 7, 1), (7, 8, 2), (8, 0, 1)],
                                    [1] * 9))
    with pytest.raises(wm.ConjugatedRingError):
        wm.component_forms(ring, ring.components.ids()[0])
    (ring_mol,) = rm.molecules_from_state(ring)
    assert rm.product_key(rm.molecule_rows(ring_mol))


# --------------------------------------------------------------------------
# homolysis
# --------------------------------------------------------------------------
def test_homolysis_linear_c16_every_bond_and_each_methyl_ch(fam):
    """n-C16H34: 15 C-C + 34 C-H = 49 bonds; R_Recombination reverse (bond
    dissociation) maps each bond in both '*' orientations: 98 moves, 49
    distinct mapped products; the three C-H cuts of each methyl are three
    different moves."""
    c16 = mol("CCCCCCCCCCCCCCCC", 0)
    moves = rm.generate_moves(fam["R_Recombination"], [c16], "reverse")
    assert len(moves) == len(kept(moves)) == 98
    bonds = Counter(tuple(sorted(a for _s, _l, a in m.labels)) for m in moves)
    every = {tuple(sorted((a.id, b.id))) for a in c16.atoms for b in a.bonds}
    assert set(bonds) == every and len(every) == 49 and set(bonds.values()) == {2}
    assert len({rm.product_key(m) for m in moves}) == 49
    for methyl in (0, 15):
        cuts = [m for m in moves if methyl in {a for _s, _l, a in m.labels}
                and any(c16.atoms[a].is_hydrogen() for _s, _l, a in m.labels)]
        hydrogens = h_ids(c16, methyl)
        assert len(hydrogens) == 3
        assert {tuple(sorted(a for _s, _l, a in m.labels)) for m in cuts} == \
            {tuple(sorted((methyl, h))) for h in hydrogens}
        assert len({rm.product_key(m) for m in cuts}) == 3


def test_homolysis_closed_ps_dp8(fam):
    """Closed PS DP8 (H head and tail caps), C64H66, one RMG form.  Single
    bonds by hand: per unit c1-c2, c2-ipso, 3 aliphatic C-H, 5 aryl C-H (10);
    7 backbone links; 2 cap C-H: 8*10 + 7 + 2 = 89.  Aromatic ring bonds do not
    match the single-bond template.  89 x 2 orientations = 178 moves."""
    import windowMatcherTest as wmt
    state = wmt.materialize(wmt.ps_spec({8: 1}))
    (chain,) = rm.molecules_from_state(state)
    assert len(chain.atoms) == 130 and len(rm.reactant_forms(chain)) == 1
    moves = rm.generate_moves(fam["R_Recombination"], [chain], "reverse")
    assert len(moves) == len(kept(moves)) == 178
    bonds = Counter(tuple(sorted(a for _s, _l, a in m.labels)) for m in moves)
    assert len(bonds) == 89 and set(bonds.values()) == {2}
    singles = {tuple(sorted((a, b))) for a, n in state.atoms.items()
               for b, o in n.edges.items() if o == 1.0}
    assert set(bonds) == singles


# --------------------------------------------------------------------------
# DP3 positive control
# --------------------------------------------------------------------------
def dp3_state():
    import windowMatcherTest as wmt
    state = wmt.merge(wmt.materialize(wmt.ps_spec({3: 2})))
    wmt.remove_h(state, wmt.find(state, 0, "unit 3", "c2"))
    return state


def test_positive_control_dp3_h_abstraction_by_hand(fam):
    """DP3 PS: chain 0 lost the H of its tail c2 (benzylic radical, 49 atoms),
    one free H, chain 1 closed (C24H26).  H_Abstraction, chain 0 accepting
    from chain 1.

    Acceptor (*3) sites of chain 0, one RMG form each (4 forms): the benzylic
    carbon, two ortho and the para ring carbon.  Donor (*1-*2) C-H bonds of
    chain 1, by hand:

    =====================  ======  =====
    donor site             C-H     sp3
    =====================  ======  =====
    unit 1 c1 (head CH3)   3       yes
    unit 1 c2 (CH)         1       yes
    unit 2 c1 (CH2)        2       yes
    unit 2 c2 (CH)         1       yes
    unit 3 c1 (CH2)        2       yes
    unit 3 c2 (tail CH2)   2       yes
    3 phenyls, 5 C-H each  15      no
    total                  26      11
    =====================  ======  =====

    Moves of the pair: 4 acceptor sites x 26 = 104, 104 distinct mapped
    products.  No move in the other order (chain 1 is closed).  Two of them
    are RMG ``unchanged_species``: the benzylic acceptor taking a tail-CH2 H
    of chain 1 turns chain 0 into chain 1's species and chain 1 into chain
    0's (an identity reaction for RMG, a real mapped H transfer for the kMC):
    102 kept.  The S2 scope (benzylic acceptor, sp3 donor) has 11 moves, 9
    kept + those 2; the matcher finds all 11 (2 matcher-only, RMG filter).
    """
    import windowMatcherTest as wmt
    state = dp3_state()
    mols = rm.molecules_from_state(state)
    chain0 = next(m for m in mols if len(m.atoms) == 49)
    chain1 = next(m for m in mols if len(m.atoms) == 50)
    moves = rm.generate_moves(fam["H_Abstraction"], [chain0, chain1])
    assert len(moves) == 104 and rm.status_counts(moves)["kept"] == 102
    assert {m.order for m in moves} == {(1, 0)}
    benzylic = wmt.find(state, 0, "unit 3", "c2")
    graph = harness.Graph(mols)
    sites = Counter(m.atom("*3") for m in moves)
    assert len(sites) == 4 and set(sites.values()) == {26} and sites[benzylic] == 26
    donors = Counter()
    for m in moves:
        if m.atom("*3") != benzylic:
            continue
        c = m.atom("*1")
        o = state.origin[c]
        name = "aryl" if c in graph.ring else f"{o.segment} {o.key}"
        donors[name] += 1
    assert donors == {"unit 1 c1": 3, "unit 1 c2": 1, "unit 2 c1": 2, "unit 2 c2": 1,
                      "unit 3 c1": 2, "unit 3 c2": 2, "aryl": 15}
    tail = wmt.find(state, 1, "unit 3", "c2")
    identity = [m for m in moves if m.status == "unchanged_species"]
    assert len(identity) == 2
    assert {(m.atom("*3"), m.atom("*1")) for m in identity} == {(benzylic, tail)}
    assert len({rm.product_key(m) for m in moves}) == 104
    row, mismatches, _r, _m = harness.compare(state, "h_abstraction_benzylic", "DP3")
    assert (row["matched"], row["matcher_only"], row["rmg_only"]) == (9, 2, 0), mismatches
    assert all(x["class"].startswith("RMG filter: unchanged_species") for x in mismatches)


# --------------------------------------------------------------------------
# comparison only: RMG's final degeneracy
# --------------------------------------------------------------------------
DEGENERACY_CASES = [
    ("H_Abstraction", ("[CH3]", "CC"), "forward"),
    ("R_Recombination", ("[CH3]", "[CH3]"), "forward"),
    ("Disproportionation", ("C[CH2]", "C[CH2]"), "forward"),
    ("R_Recombination", ("CCCCCC",), "reverse"),
    ("R_Addition_MultipleBond", ("CC[CH]CC",), "reverse"),
]


def test_degeneracy_comparison_table_is_reported(fam):
    """RMG's final degeneracy beside the directed move count (not asserted)."""
    db = harness.database()
    print("DEGTABLE family | template | products | directed | distinct mapped | RMG degeneracy")
    for family, smiles, direction in DEGENERACY_CASES:
        reactants = [mol(s, 100 * k) for k, s in enumerate(smiles)]
        moves = rm.generate_moves(fam[family], reactants, direction)
        rxns = rm.rmg_degeneracy(db.kinetics, reactants, family, direction)
        table = rm.degeneracy_table(moves, rxns)
        assert sum(r["directed_moves"] for r in table) == len(kept(moves))
        for r in table:
            print("DEGTABLE", family, direction, "|", ",".join(r["template"])[:60], "|",
                  " + ".join(r["products"]), "|", r["directed_moves"], "|",
                  r["distinct_mapped"], "|", r["rmg_degeneracy"])


# --------------------------------------------------------------------------
# determinism
# --------------------------------------------------------------------------
FINGERPRINT = """
import hashlib, sys
sys.path.insert(0, {test_dir!r})
import rmg_moves_harness as h
from rmgpy.molecule import Molecule
from rmgpy.kmc import rmg_moves as rm
def mol(s, start):
    m = Molecule().from_smiles(s)
    for k, a in enumerate(m.atoms): a.id = start + k
    return m
f = h.families()
out = []
out += rm.generate_moves(f["Disproportionation"], [mol("C[CH2]", 0), mol("C[CH2]", 10)])
out += rm.generate_moves(f["H_Abstraction"], [mol("[CH2]c1ccccc1", 0), mol("CC", 100)])
out += rm.generate_moves(f["R_Recombination"], [mol("CCCCCC", 0)], "reverse")
text = repr([(m, rm.product_key(m)) for m in out if m.status == "kept"])
print(hashlib.sha256(text.encode()).hexdigest())
"""


def test_output_is_identical_across_processes_and_hash_seeds():
    code = FINGERPRINT.format(test_dir=str(TEST_DIR))
    digests = set()
    for seed in ("0", "12345"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                             env=env, cwd=str(ROOT), timeout=600)
        assert out.returncode == 0, out.stderr[-3000:]
        digests.add(out.stdout.strip().splitlines()[-1])
    assert len(digests) == 1


def test_moves_are_sorted_by_integer_keys(fam):
    moves = rm.generate_moves(fam["H_Abstraction"], [mol("[CH2]c1ccccc1", 0), mol("CC", 100)])
    assert list(moves) == sorted(moves, key=rm.Move.sort_key)
    for m in moves:
        assert all(isinstance(a, int) for _s, _l, a in m.labels)
        assert all(isinstance(a, int) for a, _r in m.products)


# --------------------------------------------------------------------------
# window matcher: the view cache must not keep states alive
# --------------------------------------------------------------------------
def test_dropped_melt_state_is_collected_after_match():
    import windowMatcherTest as wmt
    from rmgpy.kmc import window_matcher as wm

    state = wmt.molecules(wmt.methyl(), wmt.ethane())
    result = wm.match([wmt.methyl_acceptor_ethane_donor()], state)
    assert len(result.transitions) == 6
    assert state in wm._VIEWS
    ref = weakref.ref(state)
    del state
    gc.collect()
    assert ref() is None
    assert all(k is not None for k in list(wm._VIEWS.keys()))


# --------------------------------------------------------------------------
# slow: S2 acceptance states against the S2 template sets
# --------------------------------------------------------------------------
# Pinned per case: (matched, matcher-only of class "RMG filter").  rmg-only and
# genuine misses must be zero everywhere.  The 43 H-abstraction matcher-only
# transitions are RMG ``forbidden`` (disprop1_base_case and relatives: the
# donor H sits next to a radical centre, which RMG files under
# Disproportionation).
EXPECTED = {
    8: {"h_abstraction_benzylic": (260, 43), "beta_scission_midchain": (6, 0),
        "recombination_chain_end": (5, 0), "disproportionation": (43, 0),
        "h_shift": (5, 0), "h_atom_recombination": (10, 0)},
    12: {"h_abstraction_benzylic": (410, 43), "beta_scission_midchain": (6, 0),
         "recombination_chain_end": (5, 0), "disproportionation": (43, 0),
         "h_shift": (6, 0), "h_shift4": (12, 0), "h_shift6": (12, 0),
         "h_atom_recombination": (10, 0)},
}
S2_CASES = harness.CONTRACT_CASES + ("h_shift4", "h_shift6")
# H abstraction is not run at DP21: its generation took 269 s at DP12 and grows
# with (pairs x mappings), so it would exceed the 5 min per-state budget.
DP21_SKIP = {"h_abstraction_benzylic"}


@pytest.mark.slow
@pytest.mark.skipif(not SLOW, reason="set RMG_KMC_SLOW=1")
@pytest.mark.parametrize("case", S2_CASES)
@pytest.mark.parametrize("dp", [8, 12, 21])
def test_s2_rmg_kept_moves_equal_matcher_transitions(case, dp):
    import windowMatcherTest as wmt
    if dp == 21 and case in DP21_SKIP:
        pytest.skip("DP21 H abstraction exceeds the 5 min generation budget")
    state = wmt.state_at(dp)
    t0 = time.perf_counter()
    row, mismatches, rmg, _mat = harness.compare(state, case, f"DP{dp}")
    row["wall_s"] = round(time.perf_counter() - t0, 1)
    print("S2ROW", json.dumps(row))
    for x in mismatches:
        print("S2MIS", json.dumps(x, default=str))
    assert row["rmg_only"] == 0 and row["genuine_misses"] == 0, mismatches[:5]
    assert row["matcher_only"] == row["matcher_only_rmg_filter"]
    assert row["rmg_generate_s"] < 600
    pinned = EXPECTED.get(dp, {}).get(case)
    if pinned is not None:
        assert (row["matched"], row["matcher_only_rmg_filter"]) == pinned
