"""I-316: no-sextet exclusion in the kMC resonance form model.

RMG's resonance generation (which the certification keys use) drops every
structure with no aromatic ring once a molecule has two or more aromatic
rings (``filtration.aromaticity_filtration``).  The kMC form model now has the
matching rule, judged per conjugated subsystem: in a subsystem whose best form
holds two or more sextets, a form with zero sextets is not counted, and the
search does not expand it, so the 64-form budget counts only kept forms.

A *sextet* is a 6-cycle drawn at 1.5 (``_render``: its pi values alternate).
A *class* is one isomorphism class of the whole molecule carrying the form's
bond orders and radicals (pk-i315 README section 1): kMC forms are written
onto a copy of the RMG molecule and grouped by ``Molecule.is_isomorphic``.
RMG structures are rendered the module's way (``_kekulize`` then ``_render``)
and grouped the same way.

Pieces are built from SMILES and cut into subsystems by ``_subsystems`` and
``_sub_shape``, the path ``_component_subs`` takes.  The 16-atom piece is the
one that stopped readiness run c2_h1 seed 1 (pk-i313 ``ready/c2_h1_s1.json``
``over[0]``; pk-i315 section 4 shows this SMILES is a form of that shape).
"""

import itertools

import pytest

from rmgpy.kmc import window_matcher as wm
from rmgpy.molecule import Molecule

DPM = "[CH](c1ccccc1)c1ccccc1"
P14 = "[CH](c1ccccc1)c1cccc([CH2])c1"
P15 = "[CH](c1cccc([CH2])c1)c1cccc([CH2])c1"
P16 = "[CH](c1cccc([CH2])c1)c1cc([CH2])cc([CH2])c1"
P16_QUINOID = "C=C1[CH]C(=C)C=C([CH]C2=CC(=C)[CH]C=C2)[CH]1"   # a zero-sextet form of P16
# the 14-atom piece joined by a saturated CH2-CH2 linker to an intact phenyl
P14_PLUS_PHENYL = "[CH](c1ccc(CCc2ccccc2)cc1)c1cccc([CH2])c1"
# two 6-cycles, best form one sextet: a phenyl CH2* on a para-quinodimethane
FALLBACK = "[CH2]c1ccc(C=C2C=CC(=C)C=C2)cc1"
# the 14-atom piece with an aryl sigma radical at the para position of ring A
P14_SIGMA = "[c]1ccc(cc1)[CH]c1cccc([CH2])c1"
BIG = 10 ** 6


def polyene(k):
    """Diphenylmethyl carrying a polyene of k C=C: one radical, so no form can break both rings."""
    return "[C](c1ccccc1)(c1ccccc1)" + "C=C" * k


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def molecule(smiles):
    return Molecule().from_smiles(smiles)


def subsystems(mol):
    """(group, bonds, shape) of every conjugated subsystem of ``mol``, atoms as indices into ``mol.atoms``."""
    idx = {a: k for k, a in enumerate(mol.atoms)}

    def value(k):
        a = mol.atoms[k]
        return (a.element.symbol, a.charge, a.lone_pairs, a.radical_electrons,
                {idx[b]: float(o.order) for b, o in a.bonds.items()})

    out = []
    for group, free in wm._subsystems(list(range(len(mol.atoms))), value):
        bonds, shape = wm._sub_shape(group, free, value)
        out.append((group, bonds, shape))
    return out


def ring_sub(mol):
    """The subsystem with the most atoms."""
    return max(subsystems(mol), key=lambda s: len(s[0]))


def six_cycles(shape):
    free, bonds = shape[0], shape[1]
    n = len(free)
    inc = [[] for _ in range(n)]
    for k, (i, j) in enumerate(bonds):
        inc[i].append((k, j))
        inc[j].append((k, i))
    return inc, tuple(c for c in wm._ring_check_and_cycles(n, bonds, inc) if len(c) == 6)


def sextets(form, six):
    orders = form[0]
    return sum(1 for c in six if all(orders[k] == 1.5 for k in c))


def forms_of(shape, budget=BIG):
    wm._SHAPE_CACHE.clear()
    return wm._shape_forms(shape, budget)


def with_form(mol, group, bonds, form):
    """Copy of ``mol`` with the subsystem's bond orders and radicals set to ``form``."""
    m = mol.copy(deep=True)
    orders, rads = form
    for a, r in zip(group, rads):
        m.atoms[a].radical_electrons = r
    for (a, b), o in zip(bonds, orders):
        m.atoms[a].bonds[m.atoms[b]].order = o
    m.update_atomtypes(log_species=False, raise_exception=False)
    m.update_multiplicity()
    return m


def rendered_rmg(structure):
    """One RMG structure with every subsystem rendered the module's way."""
    m = structure.copy(deep=True)
    for group, bonds, shape in subsystems(m):
        free, sbonds, orders, rads = shape
        inc, six = six_cycles(shape)
        pi = wm._kekulize(free, sbonds, orders, rads, inc)
        m = with_form(m, group, bonds, wm._render(pi, rads, six))
    return m


def classes(mols):
    reps = []
    for m in mols:
        if not any(m.is_isomorphic(x) for x in reps):
            reps.append(m)
    return reps


def same_classes(a, b):
    return (all(any(x.is_isomorphic(y) for y in b) for x in a)
            and all(any(y.is_isomorphic(x) for x in a) for y in b))


def kmc_classes(smiles):
    mol = molecule(smiles)
    group, bonds, shape = ring_sub(mol)
    return classes([with_form(mol, group, bonds, f) for f in forms_of(shape)])


def rmg_structures(smiles, **kw):
    return [rendered_rmg(s) for s in
            molecule(smiles).generate_resonance_structures(keep_isomorphic=False, **kw)]


# --------------------------------------------------------------------------
# the rule on the measured pieces
# --------------------------------------------------------------------------
@pytest.mark.parametrize("smiles,forms,n_classes", [
    (DPM, 7, 3), (P14, 12, 11), (P15, 17, 9), (P16, 25, 18)])
def test_two_ring_pieces_keep_only_forms_with_a_sextet(smiles, forms, n_classes):
    mol = molecule(smiles)
    group, bonds, shape = ring_sub(mol)
    _inc, six = six_cycles(shape)
    got = forms_of(shape)
    assert len(got) == forms
    assert all(sextets(f, six) >= 1 for f in got)
    assert max(sextets(f, six) for f in got) == 2
    assert sextets(got[0], six) == 2    # the first form is still the most aromatic
    assert len(classes([with_form(mol, group, bonds, f) for f in got])) == n_classes


@pytest.mark.parametrize("smiles,n_classes", [(DPM, 3), (P14, 11), (P15, 9)])
def test_kept_classes_equal_rmg_classes(smiles, n_classes):
    kmc = kmc_classes(smiles)
    rmg = classes(rmg_structures(smiles))
    assert len(kmc) == len(rmg) == n_classes
    assert same_classes(kmc, rmg)


def test_zero_sextet_start_reaches_the_same_forms():
    """A subsystem stored in a quinoid placement (the post-edit canonicaliser
    builds its shape from the raw placement) still gives the kept forms."""
    quinoid = molecule(P16_QUINOID)
    q_group, q_bonds, q_shape = ring_sub(quinoid)
    _inc, six = six_cycles(q_shape)
    start = wm._render(wm._kekulize(q_shape[0], q_shape[1], q_shape[2], q_shape[3], six_cycles(q_shape)[0]),
                       q_shape[3], six)
    assert sextets(start, six) == 0
    got = forms_of(q_shape)
    assert len(got) == 25
    assert all(sextets(f, six) >= 1 for f in got)
    mol = molecule(P16)
    group, bonds, shape = ring_sub(mol)
    assert same_classes(classes([with_form(quinoid, q_group, q_bonds, f) for f in got]),
                        classes([with_form(mol, group, bonds, f) for f in forms_of(shape)]))


def test_budget_counts_only_kept_forms():
    _g, _b, shape = ring_sub(molecule(P16))
    assert len(forms_of(shape, 25)) == 25
    with pytest.raises(wm.ResonanceLimitError):
        forms_of(shape, 24)
    assert len(forms_of(shape, wm.RESONANCE_FORM_LIMIT)) == 25
    # the cached answer is bounded the same way
    with pytest.raises(wm.ResonanceLimitError):
        wm._shape_forms(shape, 24)


def test_component_subs_counts_the_kept_forms():
    mol = molecule(P16)
    idx = {a: k for k, a in enumerate(mol.atoms)}

    def value(k):
        a = mol.atoms[k]
        return (a.element.symbol, a.charge, a.lone_pairs, a.radical_electrons,
                {idx[b]: float(o.order) for b, o in a.bonds.items()})

    wm._SHAPE_CACHE.clear()
    subs, total = wm._component_subs(list(range(len(mol.atoms))), value)
    assert total == 25 and [len(s.forms) for s in subs] == [25]


# --------------------------------------------------------------------------
# where the rule does not apply
# --------------------------------------------------------------------------
@pytest.mark.parametrize("smiles,forms,zero", [
    ("[CH2]c1ccccc1", 4, 3), ("[CH2]c1cccc([CH2])c1", 9, 8), ("[CH2]c1ccc([CH2])cc1", 9, 8)])
def test_one_ring_subsystem_is_unchanged(smiles, forms, zero):
    """One 6-cycle: the best form holds one sextet, so zero-sextet forms stay (counts as at efe849fe6)."""
    _g, _b, shape = ring_sub(molecule(smiles))
    _inc, six = six_cycles(shape)
    got = forms_of(shape)
    assert len(got) == forms
    assert sum(1 for f in got if sextets(f, six) == 0) == zero


def test_two_rings_whose_best_form_has_one_sextet_keep_zero_sextet_forms():
    """Two 6-cycles, but the quinodimethane ring can never be a sextet while the
    phenyl is one: the best form holds one sextet, so nothing is dropped."""
    _g, _b, shape = ring_sub(molecule(FALLBACK))
    _inc, six = six_cycles(shape)
    assert len(six) == 2
    got = forms_of(shape)
    assert len(got) == 8
    assert max(sextets(f, six) for f in got) == 1
    assert sum(1 for f in got if sextets(f, six) == 0) == 6


def test_declared_difference_per_subsystem_drops_what_rmg_keeps():
    """The declared difference.  kMC judges each conjugated subsystem on its
    own sextets: the two-ring poly-radical subsystem drops its zero-sextet
    forms although a separate phenyl elsewhere in the molecule is aromatic.
    RMG judges the whole molecule: the phenyl stays aromatic in every
    structure, so RMG keeps the structures in which the poly-radical part has
    no aromatic ring."""
    mol = molecule(P14_PLUS_PHENYL)
    subs = sorted(subsystems(mol), key=lambda s: len(s[0]))
    assert [len(s[0]) for s in subs] == [6, 14]
    (_pg, _pb, phenyl), (group, bonds, shape) = subs
    assert len(forms_of(phenyl)) == 1
    _inc, six = six_cycles(shape)
    kept = forms_of(shape)
    assert len(kept) == 12
    assert all(sextets(f, six) >= 1 for f in kept)
    # RMG: structures whose poly-radical part has no ring drawn at 1.5.  RMG
    # reorders atoms, so each structure is judged on its own subsystem.
    rmg = classes(rmg_structures(P14_PLUS_PHENYL))

    def part_sextets(m):
        (_g2, _b2, part), = [s for s in subsystems(m) if len(s[0]) == 14]
        _i2, six2 = six_cycles(part)
        return sum(1 for c in six2 if all(part[2][k] == 1.5 for k in c))

    zero = [m for m in rmg if part_sextets(m) == 0]
    assert len(rmg) == 17 and len(zero) == 6
    base = rendered_rmg(mol)      # the phenyl subsystem drawn the module's way too
    kmc = classes([with_form(base, group, bonds, f) for f in kept])
    assert len(kmc) == 11
    assert same_classes(kmc, [m for m in rmg if part_sextets(m) >= 1])


def test_aryl_sigma_radical_is_still_held_out():
    """The rule leaves ``_aryl_sigma_atoms`` as it is: the sigma radical stays on its atom in every kept form."""
    mol = molecule(P14_SIGMA)
    group, _bonds, shape = ring_sub(mol)
    free, sbonds, _orders, rads = shape
    inc, six = six_cycles(shape)
    sigma = wm._aryl_sigma_atoms(free, rads, inc, wm._ring_check_and_cycles(len(free), sbonds, inc))
    assert len(sigma) == 1
    (a,) = sigma
    got = forms_of(shape)
    assert len(got) == 12
    assert all(f[1][a] >= 1 for f in got)
    assert all(sextets(f, six) >= 1 for f in got)


# --------------------------------------------------------------------------
# the limits stay
# --------------------------------------------------------------------------
def test_limits_unchanged():
    assert wm.RESONANCE_FORM_LIMIT == 64
    assert wm.RAW_PLACEMENT_LIMIT == 4096


def test_over_limit_subsystem_without_zero_sextet_forms_still_raises():
    """One radical over two rings and a polyene: no form can break both rings,
    so the rule drops nothing.  64 forms pass; 65 raise at the unchanged limit."""
    _g, _b, shape = ring_sub(molecule(polyene(57)))
    _inc, six = six_cycles(shape)
    got = forms_of(shape)
    assert len(got) == 64 and all(sextets(f, six) >= 1 for f in got)
    assert len(forms_of(shape, wm.RESONANCE_FORM_LIMIT)) == 64
    _g, _b, shape = ring_sub(molecule(polyene(58)))
    assert len(forms_of(shape)) == 65
    with pytest.raises(wm.ResonanceLimitError):
        forms_of(shape, wm.RESONANCE_FORM_LIMIT)


# --------------------------------------------------------------------------
# the alternation test is one helper, with the old truth table
# --------------------------------------------------------------------------
def _old_alternates(values):
    return all(v in (0, 1) for v in values) and all(values[i] != values[i - 1] for i in range(len(values)))


def _old_render(pi, radicals, six_cycles):
    """``_render`` as it was at efe849fe6, verbatim."""
    orders = [1.0 + p for p in pi]
    for cycle in six_cycles:
        values = [pi[k] for k in cycle]
        if all(v in (0, 1) for v in values) and all(values[i] != values[i - 1] for i in range(6)):
            for k in cycle:
                orders[k] = 1.5
    return tuple(orders), tuple(radicals)


def test_alternates_truth_table():
    for n in (4, 6, 8):
        for values in itertools.product((0, 1, 2), repeat=n):
            assert wm._alternates(values) == _old_alternates(values), values


@pytest.mark.parametrize("smiles", [DPM, P14, P15, P16])
def test_render_agrees_with_the_old_render_on_every_raw_state(smiles):
    _g, _b, shape = ring_sub(molecule(smiles))
    free, bonds, orders, rads = shape
    inc, six = six_cycles(shape)
    even = tuple(c for c in wm._ring_check_and_cycles(len(free), bonds, inc) if len(c) % 2 == 0)
    start = (wm._kekulize(free, bonds, orders, rads, inc), tuple(rads))
    seen, queue = {start}, [start]
    while queue:
        pi, rad = queue.pop()
        assert wm._render(pi, rad, six) == _old_render(pi, rad, six)
        for nxt in wm._moves(pi, rad, inc, even):
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
    assert len(seen) >= 16
