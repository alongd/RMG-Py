"""The embedding quotient folds interchangeable branches by the template's own
symmetry, in the matcher and in the engine alike (I-238).

The matcher and the engine enumerate one embedding per class of embeddings
that differ only by a symmetry of the template window that fixes its edited
atoms: the class's first embedding in the full enumeration order, standing
for the class's size.  The hydrogen-leaf quotient was the special case of
hydrogen leaves of one parent; it missed every other symmetric branch (a
central C-C homolysis on linear C22F46 has 2 * (3!)**2 * (2!)**20 ~ 7.5e7
embeddings and one transition), and the engine grouped hydrogen leaves by
parent alone, so a deuterium among a carbon's hydrogens could leave a
template with no embedding at all.
"""

import itertools

import pytest

from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import window_matcher as wm
from rmgpy.kmc.melt_state import MeltState, Node, Origin
from rmgpy.kmc.window_matcher import AtomEdit, BondEdit, Template, WindowAtom


def _state(atom_list, bonds, ids=None):
    """A MeltState of one molecule: atoms (element, isotope), bonds (i, j,
    order) on local indices; ``ids``: the state id of each local index."""
    ids = list(range(len(atom_list))) if ids is None else list(ids)
    edges = {ids[i]: {} for i in range(len(atom_list))}
    for i, j, order in bonds:
        edges[ids[i]][ids[j]] = float(order)
        edges[ids[j]][ids[i]] = float(order)
    atoms, origin = {}, {}
    for i, (element, isotope) in enumerate(atom_list):
        atoms[ids[i]] = Node(element, isotope, 0, 0, LONE_PAIRS.get(element, 0), edges[ids[i]])
        origin[ids[i]] = Origin("mobile", 0, "m0", f"a{i}", False)
    return MeltState(atoms, origin)


LONE_PAIRS = {"F": 3}


def _template(name, atom_list, bonds, edits):
    atoms = [WindowAtom(element, isotope=isotope, lone_pairs=LONE_PAIRS.get(element, 0))
             for element, isotope in atom_list]
    return Template(name, "homolysis", atoms, bonds, edits, weight=1.0)


def _homolysis(atom_list, bonds, a, b):
    return [BondEdit(a, b, 0.0), AtomEdit(a, radical=1), AtomEdit(b, radical=1)]


def _perfluoro_chain(n):
    """Linear CnF(2n+2): carbons 0..n-1, then each carbon's fluorine leaves."""
    atom_list = [("C", None)] * n
    bonds = [(i, i + 1, 1.0) for i in range(n - 1)]
    for c in range(n):
        for _ in range(3 if c in (0, n - 1) else 2):
            atom_list = atom_list + [("F", None)]
            bonds.append((c, len(atom_list) - 1, 1.0))
    return atom_list, bonds


def _classes(template, participant, embeddings):
    """Embeddings up to the template's symmetries that fix its edited atoms:
    the images of the fixed (edited) atoms with each embedding's image set."""
    edited = template.edited[participant]
    return {(e.cid, tuple(e.mapping[t] for t in edited), frozenset(a for a in e.mapping if a is not None),
             e.choices, e.per_choice) for e in embeddings}


def _full_match(template, state):
    """``match`` with the quotient replaced by the full enumeration."""
    original = wm._participant_quotient
    wm._participant_quotient = lambda view, t, p: (wm._participant_embeddings(view, t, p), 1)
    try:
        return wm.match([template], state)
    finally:
        wm._participant_quotient = original


def test_a_perfluoro_chain_homolysis_keeps_one_embedding_per_orientation():
    """Central C-C homolysis on linear C6F14: 2 * (3!)**2 * (2!)**4 = 1152
    embeddings, two orientations of the cut bond, one transition.  The
    quotient keeps the first embedding of each orientation, each standing
    for 576, and ``match`` is the full enumeration's."""
    atom_list, bonds = _perfluoro_chain(6)
    template = _template("cut_c6f14", atom_list, bonds, _homolysis(atom_list, bonds, 2, 3))
    state = _state(atom_list, bonds)
    view = wm._view_for(state)
    full = wm._participant_embeddings(view, template, 0)
    assert len(full) == 1152
    kept, copies = wm._participant_quotient(view, template, 0)
    assert len(kept) == 2 and len(kept) * copies == len(full)
    assert [e for e in full if e in kept] == kept                          # a subsequence ...
    assert all(min(f.mapping for f in full if _classes(template, 0, [f]) == _classes(template, 0, [e]))
               == e.mapping for e in kept)                                 # ... of first embeddings
    fast, slow = wm.match([template], state), _full_match(template, state)
    assert len(fast.transitions) == 1
    assert fast.transitions == slow.transitions
    assert fast.classes == slow.classes


def test_a_phenyl_ring_flip_is_folded():
    """Ethylbenzene, the ethyl C-C bond cut with the whole ring in the
    window: the ring flip (2), the methyl hydrogens (3!) and the CH2
    hydrogens (2!) give 24 embeddings of one class."""
    heavy = [("C", None)] * 8                        # 0-5 ring, 6 CH2, 7 CH3
    bonds = [(i, (i + 1) % 6, 1.5) for i in range(6)] + [(0, 6, 1.0), (6, 7, 1.0)]
    atom_list = list(heavy)
    for c, n in ((1, 1), (2, 1), (3, 1), (4, 1), (5, 1), (6, 2), (7, 3)):
        for _ in range(n):
            atom_list.append(("H", None))
            bonds.append((c, len(atom_list) - 1, 1.0))
    template = _template("cut_ethyl", atom_list, bonds, _homolysis(atom_list, bonds, 6, 7))
    state = _state(atom_list, bonds)
    view = wm._view_for(state)
    full = wm._participant_embeddings(view, template, 0)
    kept, copies = wm._participant_quotient(view, template, 0)
    assert len(full) == 24 and len(kept) == 1 and copies == 24
    assert kept[0] == full[0]
    assert wm.match([template], state).transitions == _full_match(template, state).transitions


def _deuterated_ethane(d_template, d_state):
    """CH2D-CH3 with the D leaf at template position ``d_template`` and state
    position ``d_state`` among carbon 0's three leaves (2, 3, 4)."""
    leaves = [("H", None)] * 3
    leaves[d_template - 2] = ("H", 2)
    atom_list = [("C", None), ("C", None)] + leaves + [("H", None)] * 3
    bonds = [(0, 1, 1.0)] + [(0, k, 1.0) for k in (2, 3, 4)] + [(1, k, 1.0) for k in (5, 6, 7)]
    template = _template("cut_ethane_d", atom_list, bonds, _homolysis(atom_list, bonds, 0, 1))
    state_leaves = [("H", None)] * 3
    state_leaves[d_state - 2] = ("H", 2)
    state_atoms = [("C", None), ("C", None)] + state_leaves + [("H", None)] * 3
    return template, state_atoms, bonds


@pytest.mark.parametrize("d_template", [2, 3, 4])
@pytest.mark.parametrize("d_state", [2, 3, 4])
def test_a_deuterium_leaf_never_empties_the_engine_quotient(d_template, d_state):
    """Whole-ethane C-C homolysis with one D among a carbon's hydrogens:
    under every placement of the D leaf in the template and every
    permutation of the state's atom ids, the engine's quotient keeps exactly
    the full enumeration's classes (not none), and the matcher one
    transition."""
    template, atom_list, bonds = _deuterated_ethane(d_template, d_state)
    for ids in itertools.permutations(range(8)):
        if ids[0] > ids[1] or list(ids[5:]) != sorted(ids[5:]):
            continue                                 # the C-C swap and the CH3 leaves change nothing here
        state = _state(atom_list, bonds, ids)
        view = wm._view_for(state)
        full = me.participant_embeddings(view, template, 0, quotient=False)
        fast = me.participant_embeddings(view, template, 0, quotient=True)
        assert full and [e for e in full if e in fast] == fast
        assert _classes(template, 0, fast) == _classes(template, 0, full), ids
        assert len(wm.match([template], state).transitions) == 1


def test_the_engine_quotient_folds_fluorine_leaves():
    atom_list, bonds = _perfluoro_chain(6)
    template = _template("cut_c6f14", atom_list, bonds, _homolysis(atom_list, bonds, 2, 3))
    view = wm._view_for(_state(atom_list, bonds))
    fast = me.participant_embeddings(view, template, 0, quotient=True)
    full = me.participant_embeddings(view, template, 0, quotient=False)
    assert len(fast) == 2 and _classes(template, 0, fast) == _classes(template, 0, full)
