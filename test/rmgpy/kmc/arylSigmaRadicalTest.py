"""I-313: an aryl sigma-radical stays on its atom in the kMC resonance form model.

A ring carbon whose two subsystem bonds both lie on a conjugated cycle and
whose free valence is 2 puts one free valence into the ring's pi system; its
other free valence is the in-plane sp2 orbital, orthogonal to pi.  A radical
there is a sigma electron: it never shifts through the ring.  Every other
radical delocalizes as before.

The oracle is RMG's ``generate_resonance_structures``, rendered the way this
module counts forms (``_kekulize`` then ``_render``: an alternating 6-cycle
becomes 1.5) and reduced to isomorphism classes of the whole molecule:

- phenyl and benzyl: the kMC classes equal RMG's;
- a sigma-radical beside a pi radical: RMG localizes an aryl radical only when
  every radical of the molecule is aryl (``populate_resonance_algorithms``),
  so here it also lists forms that move the sigma electron into the ring.
  kMC departs from RMG there: every RMG structure that keeps the sigma
  radical on its atom is a kMC class, none that moves it is, and the kMC forms
  are exactly those of the same shape with the sigma electron removed, then
  added back.
"""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from rmgpy.kmc import window_matcher as wm  # noqa: E402

EXACT = 10 ** 6


# --------------------------------------------------------------------------
# the form model before I-313, verbatim (no cache): the frozen reference
# --------------------------------------------------------------------------
def legacy_shape_forms(shape, budget=EXACT):
    free, bonds, orders, radicals = shape
    n = len(free)
    incident = [[] for _ in range(n)]
    for k, (i, j) in enumerate(bonds):
        incident[i].append((k, j))
        incident[j].append((k, i))
    cycles = wm._ring_check_and_cycles(n, bonds, incident)
    six = tuple(c for c in cycles if len(c) == 6)
    even = tuple(c for c in cycles if len(c) % 2 == 0)
    start = (wm._kekulize(free, bonds, orders, radicals, incident), tuple(radicals))
    seen = {start}
    queue = [start]
    rendered = {}
    head = 0
    while head < len(queue):
        pi, rad = queue[head]
        head += 1
        form = wm._render(pi, rad, six)
        if form not in rendered:
            rendered[form] = len(rendered)
            if len(rendered) > budget:
                raise wm.ResonanceLimitError(f"more than {budget} resonance forms", budget + 1)
        for nxt in wm._moves(pi, rad, incident, even):
            if nxt not in seen:
                if len(seen) >= wm.RAW_PLACEMENT_LIMIT:
                    raise wm.ResonanceLimitError(
                        f"more than {wm.RAW_PLACEMENT_LIMIT} pi placements in one subsystem")
                seen.add(nxt)
                queue.append(nxt)
    return tuple(sorted(rendered, key=wm._form_order))


def no_sextet(shape, forms):
    """The I-316 no-sextet exclusion as a filter over a complete form set: in a
    shape with two or more 6-cycles whose best form holds two or more sextets,
    the forms with zero sextets are dropped.  ``_shape_forms`` applies the rule
    during its search; this filter, applied after the frozen legacy search, is
    the independent reference it must equal."""
    free, bonds = shape[0], shape[1]
    incident = [[] for _ in free]
    for k, (i, j) in enumerate(bonds):
        incident[i].append((k, j))
        incident[j].append((k, i))
    six = [c for c in wm._ring_check_and_cycles(len(free), bonds, incident) if len(c) == 6]

    def sextets(form):
        return sum(1 for c in six if all(form[0][k] == 1.5 for k in c))

    if len(six) < 2 or max(sextets(f) for f in forms) < 2:
        return forms
    return tuple(f for f in forms if sextets(f) >= 1)


def forms(shape):
    wm._SHAPE_CACHE.clear()
    return wm._shape_forms(shape, EXACT)


def pinned_reference(shape, sigma):
    """The forms of ``shape`` with one electron held at each atom of
    ``sigma``: the legacy forms of the shape without it, then added back,
    under the no-sextet exclusion (I-316)."""
    free, bonds, orders, radicals = shape
    reduced = (tuple(f - (a in sigma) for a, f in enumerate(free)), bonds, orders,
               tuple(r - (a in sigma) for a, r in enumerate(radicals)))
    back = [(o, tuple(r + (a in sigma) for a, r in enumerate(rad))) for o, rad in legacy_shape_forms(reduced)]
    return no_sextet(shape, tuple(sorted(back, key=wm._form_order)))


# --------------------------------------------------------------------------
# RMG molecules <-> shapes
# --------------------------------------------------------------------------
def _setup(mol):
    free = {a: wm._valence(a.element.symbol) - a.charge - 2 * a.lone_pairs - len(a.bonds) for a in mol.atoms}
    group = [a for a in mol.atoms if free[a] > 0]
    gi = {a: k for k, a in enumerate(group)}
    bonds = tuple(sorted((gi[a], gi[b]) for a in group for b in a.bonds if b in gi and gi[a] < gi[b]))
    n = len(group)
    incident = [[] for _ in range(n)]
    for k, (i, j) in enumerate(bonds):
        incident[i].append((k, j))
        incident[j].append((k, i))
    six = tuple(c for c in wm._ring_check_and_cycles(n, bonds, incident) if len(c) == 6)
    orders = tuple(float(group[i].bonds[group[j]].order) for i, j in bonds)
    radicals = tuple(a.radical_electrons for a in group)
    return group, (tuple(free[a] for a in group), bonds, orders, radicals), incident, six


def _as_molecule(mol, group, bonds, form):
    orders, radicals = form
    m = mol.copy(deep=True)
    at = {a: mol.atoms.index(a) for a in group}
    for (i, j), o in zip(bonds, orders):
        m.atoms[at[group[i]]].bonds[m.atoms[at[group[j]]]].order = o
    for i, r in enumerate(radicals):
        m.atoms[at[group[i]]].radical_electrons = r
    m.update_atomtypes(log_species=False, raise_exception=False)
    m.update_multiplicity()
    return m


def _rendered(mol):
    """``mol`` as this module renders it (one Kekule assignment, alternating
    6-cycles aromatic)."""
    group, (free, bonds, orders, radicals), incident, six = _setup(mol)
    pi = wm._kekulize(free, bonds, orders, radicals, incident)
    return _as_molecule(mol, group, bonds, wm._render(pi, radicals, six))


def _classes(mols):
    out = []
    for m in mols:
        if not any(m.is_isomorphic(x) for x in out):
            out.append(m)
    return out


def _member(m, classes):
    return any(m.is_isomorphic(x) for x in classes)


def _sigma_indices(mol):
    """Molecule atom indices of the aryl sigma-radical atoms (free valence 2,
    radical, two heavy bonds, in a ring).  Here every free-valence-2 atom of
    the test molecules is such an atom in the SMILES form."""
    out = []
    for i, a in enumerate(mol.atoms):
        heavy = [b for b in a.bonds if not b.is_hydrogen()]
        free = wm._valence(a.element.symbol) - a.charge - 2 * a.lone_pairs - len(a.bonds)
        if free == 2 and a.radical_electrons >= 1 and len(heavy) == 2 and mol.is_vertex_in_cycle(a):
            out.append(i)
    return out


def study(smiles):
    """(kMC forms, kMC classes, RMG classes keeping each sigma radical, RMG
    classes moving one, shape, shape-local sigma atoms)."""
    from rmgpy.molecule import Molecule
    mol = Molecule().from_smiles(smiles)
    group, shape, _incident, _six = _setup(mol)
    kforms = forms(shape)
    kmc = _classes([_as_molecule(mol, group, shape[1], f) for f in kforms])
    sigma = _sigma_indices(mol)
    local = {group.index(mol.atoms[i]) for i in sigma}
    keep, move = [], []
    for s in mol.copy(deep=True).generate_resonance_structures(keep_isomorphic=False):
        # RMG reorders atoms, so the sigma atoms are found by structure: a ring
        # atom of free valence 2 on two heavy bonds with no radical (two ring
        # pi bonds, a cumulene) is a sigma atom whose radical moved
        moved = len(_sigma_indices(s)) < len(sigma)
        (move if moved else keep).append(_rendered(s))
    return kforms, kmc, _classes(keep), _classes(move), shape, local


# --------------------------------------------------------------------------
# exact agreement with RMG where RMG localizes the aryl radical
# --------------------------------------------------------------------------
def test_phenyl_radical_has_one_form():
    kforms, kmc, keep, move, _shape, sigma = study("[c]1ccccc1")
    assert len(sigma) == 1 and not move
    assert len(kforms) == 1 and len(kmc) == len(keep) == 1
    assert kmc[0].is_isomorphic(keep[0])


def test_benzyl_radical_still_delocalizes_into_the_ring():
    kforms, kmc, keep, move, shape, sigma = study("[CH2]c1ccccc1")
    assert not sigma and not move
    assert len(kforms) == 4 and len(kmc) == len(keep) == 3
    assert all(_member(m, keep) for m in kmc) and all(_member(m, kmc) for m in keep)
    assert kforms == legacy_shape_forms(shape)
    # the radical reaches the ring: forms with no radical on the CH2 carbon
    assert sum(1 for _o, r in kforms if r[shape[3].index(1)] == 0) == 3


# --------------------------------------------------------------------------
# a sigma-radical beside a pi radical: kMC departs from RMG's whole-molecule switch
# --------------------------------------------------------------------------
def _check_departure(smiles, expected_forms):
    kforms, kmc, keep, move, shape, sigma = study(smiles)
    assert sigma
    assert move, "RMG moves the sigma radical here (its aryl switch is per molecule)"
    assert all(_member(m, kmc) for m in keep), "an RMG form keeping the sigma radical is missing"
    assert not any(_member(m, kmc) for m in move), "a form moving the sigma radical is a kMC form"
    assert kforms == pinned_reference(shape, sigma)
    assert all(rad[a] >= 1 for _o, rad in kforms for a in sigma)
    assert len(kforms) == expected_forms
    return kforms, kmc, keep, move


@pytest.mark.parametrize("smiles", ["[CH2]c1[c]cccc1", "[CH2]c1c[c]ccc1", "[CH2]c1cc[c]cc1"],
                         ids=["ortho", "meta", "para"])
def test_sigma_radical_with_benzylic_radical(smiles):
    _check_departure(smiles, 4)


# the pk-i298 section 6 subsystems, copied from the readiness probe outputs
# (/home/alon/runs/pk-i298/ready/c2_h1_s1.json and c2_h1_s7.json, over[0]:
# free_valence, shape_bonds, shape_radicals); c3_homolysis s1, s3, s8 hold the
# 14-atom shape too
PK_I298_14 = {"free_valence": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1],
              "shape_bonds": [[0, 5, 1.0], [0, 13, 1.0], [1, 9, 1.0], [2, 3, 1.5], [2, 7, 1.5], [3, 4, 1.5],
                              [4, 5, 1.5], [5, 6, 1.5], [6, 7, 1.5], [8, 9, 1.5], [8, 13, 1.5], [9, 10, 1.5],
                              [10, 11, 1.5], [11, 12, 1.5], [12, 13, 1.5]],
              "shape_radicals": [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]}
PK_I298_15 = {"free_valence": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1],
              "shape_bonds": [[0, 6, 1.0], [0, 10, 1.0], [1, 8, 1.0], [2, 14, 1.0], [3, 4, 1.5], [3, 8, 1.5],
                              [4, 5, 1.5], [5, 6, 1.5], [6, 7, 1.5], [7, 8, 1.5], [9, 10, 1.5], [9, 14, 1.5],
                              [10, 11, 1.5], [11, 12, 1.5], [12, 13, 1.5], [13, 14, 1.5]],
              "shape_radicals": [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]}
# the same pieces as molecules: two rings joined by a CH*, a CH2* on ring B and
# the aryl sigma-radical on ring B meta to it (and to the bridge); the
# 15-atom piece carries a second CH2* on ring A
# (21 and 56 forms before the I-316 no-sextet exclusion, 12 and 17 under it)
PIECES = {"14": (PK_I298_14, "[CH](c1ccccc1)c1cc([CH2])c[c]c1", 12),
          "15": (PK_I298_15, "[CH](c1cccc([CH2])c1)c1cc([CH2])c[c]c1", 17)}


def _literal(record):
    bonds = tuple((i, j) for i, j, _o in record["shape_bonds"])
    orders = tuple(float(o) for _i, _j, o in record["shape_bonds"])
    return tuple(record["free_valence"]), bonds, orders, tuple(record["shape_radicals"])


@pytest.mark.parametrize("piece", sorted(PIECES))
def test_pk_i298_subsystem_forms(piece):
    record, smiles, expected = PIECES[piece]
    source = Path("/home/alon/runs/pk-i298/ready") / ("c2_h1_s1.json" if piece == "14" else "c2_h1_s7.json")
    if source.exists():
        over = json.loads(source.read_text())["over"][0]
        assert {k: over[k] for k in record} == record
    shape = _literal(record)
    sigma = {a for a, f in enumerate(shape[0]) if f == 2}
    got = forms(shape)
    assert got == pinned_reference(shape, sigma)
    assert len(got) == expected <= wm.RESONANCE_FORM_LIMIT
    # the molecule built from SMILES has the same subsystem: same counts today and now
    kforms, _kmc, _keep, _move = _check_departure(smiles, expected)
    from rmgpy.molecule import Molecule
    smiles_shape = _setup(Molecule().from_smiles(smiles))[1]
    assert len(legacy_shape_forms(smiles_shape)) == len(legacy_shape_forms(shape)) == (65 if piece == "14" else 131)
    assert sorted(smiles_shape[0]) == sorted(shape[0]) and sorted(smiles_shape[3]) == sorted(shape[3])
    assert len(kforms) == len(got)


def test_every_form_keeps_the_sigma_electron():
    """Each form, stored as the start, gives the same form set: identity does
    not depend on which form a state holds."""
    shapes = [_literal(PK_I298_14), _literal(PK_I298_15)]
    from rmgpy.molecule import Molecule
    for smiles in ("[c]1ccccc1", "[CH2]c1[c]cccc1", "[CH2]c1c[c]ccc1", "[CH2]c1cc[c]cc1"):
        shapes.append(_setup(Molecule().from_smiles(smiles))[1])
    for shape in shapes:
        sigma = {a for a, (f, r) in enumerate(zip(shape[0], shape[3])) if f == 2 and r >= 1}
        base = forms(shape)
        assert all(rad[a] >= 1 for _o, rad in base for a in sigma)
        for orders, rad in base:
            assert forms((shape[0], shape[1], orders, rad)) == base


# --------------------------------------------------------------------------
# every other radical: unchanged
# --------------------------------------------------------------------------
@pytest.mark.parametrize("smiles", ["C=C[CH2]", "C=[CH]", "C=[C]C=C", "[O]c1ccccc1", "[CH](c1ccccc1)c1ccccc1",
                                    "[C]1=CCCCC1", "[CH]1C=CC=CC1", "C=C1C=CC=C1[CH2]",
                                    "[CH](c1ccccc1)c1cccc([CH2])c1", "[CH2]c1ccc2ccccc2c1", "c1ccccc1"])
def test_non_sigma_radicals_are_unchanged(smiles):
    """Unchanged by the sigma rule; the two-sextet shapes here (the 14-atom
    piece, 2-naphthylmethyl) lose their zero-sextet forms to I-316."""
    from rmgpy.molecule import Molecule
    mol = Molecule().from_smiles(smiles)
    _group, shape, _incident, _six = _setup(mol)
    assert forms(shape) == no_sextet(shape, legacy_shape_forms(shape))


def test_a_sigma_atom_without_a_radical_is_not_pinned():
    """Free valence 2 on a ring atom with no radical (a ring allene/aryne
    form) is not a sigma-radical: the shape enumerates as before."""
    shape = ((2, 1, 1, 1, 1, 1), ((0, 1), (0, 5), (1, 2), (2, 3), (3, 4), (4, 5)),
             (2.0, 2.0, 1.0, 2.0, 1.0, 1.0), (0, 0, 0, 0, 1, 0))
    assert forms(shape) == legacy_shape_forms(shape)


# --------------------------------------------------------------------------
# the limits
# --------------------------------------------------------------------------
def test_the_per_subsystem_and_per_window_limits_are_still_64():
    assert wm.RESONANCE_FORM_LIMIT == 64
    assert wm.RAW_PLACEMENT_LIMIT == 4096
    import resonancePerSubsystemTest as R
    R.test_a_single_subsystem_over_the_bound_still_raises()                       # per subsystem
    R.test_a_window_over_subsystems_whose_joint_forms_pass_the_bound_raises(       # per window
        ("polyenyl17", "polyenyl19"), True)
    shape = _literal(PK_I298_15)
    with pytest.raises(wm.ResonanceLimitError):
        wm._SHAPE_CACHE.clear()
        wm._shape_forms(shape, 16)
    wm._SHAPE_CACHE.clear()
    assert len(wm._shape_forms(shape, wm.RESONANCE_FORM_LIMIT)) == 17
