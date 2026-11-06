"""Product keys skip RMG's Clar optimisation where it cannot change a form (I-307).

RMG's ``generate_resonance_structures`` runs ``generate_clar_structures`` (an
MILP through scipy) instead of ``generate_aromatic_resonance_structure``
whenever a form has more than one aromatic ring, fused or not.  When no atom
lies on two rings, every perceived aromatic ring can be a sextet at once, so
the MILP's optimum is unique and writes exactly the 1.5 bonds the aromatic
branch writes.  ``rmg_moves`` therefore passes ``clar_structures=False`` for
such molecules (and no bond joins two rings, and no atom is charged or has
lone pairs) and keeps RMG's default otherwise.  These tests pin:

- the ring test (``_rings_share_an_atom``) on unfused, fused, spiro and bridged
  graphs;
- that product keys and reactant forms ask RMG for Clar exactly when two rings
  share an atom;
- that the forms, in order, and every key equal RMG's default call on
  multi-phenyl products, radicals, diradicals, cyclohexadienyls, Kekule
  variants and atom-id permutations.
"""

import random

import pytest

from rmgpy.kmc import rmg_moves as rm
from rmgpy.molecule import Molecule
from rmgpy.molecule.resonance import generate_resonance_structures


def rows_of(smiles, ids=None):
    m = Molecule().from_smiles(smiles)
    ids = list(range(len(m.atoms))) if ids is None else list(ids)
    for atom, i in zip(m.atoms, ids):
        atom.id = i
    return rm.molecule_rows(m)


def relabel(rows, mapping):
    return {mapping[i]: r[:5] + (tuple(sorted((mapping[j], o) for j, o in r[5])),)
            for i, r in rows.items()}


def default_forms(rows):
    """RMG's default call, as rmg_moves made it before I-307."""
    mol = rm.molecule_from_rows(rows)
    return tuple(generate_resonance_structures(mol, keep_isomorphic=True, filter_structures=True))


def described(forms):
    return [(rm._signature(rm.molecule_rows(f)), bool(f.reactive)) for f in forms]


def reference_key(rows):
    """The key from RMG's default call (the I-305 key algorithm, uncached)."""
    def component(part):
        sig = rm._signature(part)
        forms = default_forms(part) if len(part) > 1 else (rm.molecule_from_rows(part),)
        sigs = described(forms)
        aromatic = max(rm._aromatic_bonds(f) for f in forms)
        root = min(s for (s, _r), f in zip(sigs, forms) if rm._aromatic_bonds(f) == aromatic)
        if root != sig:
            sigs = described(default_forms(dict((r[0], r[1:]) for r in root)))
        reactive = [s for s, r in sigs if r] or [s for s, _r in sigs]
        return min(reactive)
    return tuple(sorted(component(p) for p in rm.split_rows(rows)))


@pytest.fixture
def clar_requests(monkeypatch):
    """The ``clar_structures`` value of every RMG resonance call rmg_moves makes."""
    seen = []
    real = rm.generate_resonance_structures

    def spy(*args, **kwargs):
        seen.append(kwargs.get("clar_structures", True))
        return real(*args, **kwargs)

    monkeypatch.setattr(rm, "generate_resonance_structures", spy)
    return seen


U = "CC(c1ccccc1)"
UNFUSED = {                                      # no shared ring atom, no lone pairs, neutral
    "ps_benzylic": "C" + U * 3 + "C[CH]c1ccccc1",
    "ps_diradical": "[CH2]C(c1ccccc1)" + U * 2 + "C[CH]c1ccccc1",
    "ps_cyclohexadienyl": "C" + U * 2 + "CC(C1C=C[CH]C=C1)C",
    "ps_alkene": "C" + U * 3 + "C=Cc1ccccc1",
    "ps_aryl": "C" + U * 2 + "CC(c1cc[c]cc1)C",
    "ps_ortho_methylene": "C" + U * 2 + "CC(c1ccccc1[CH2])C",
    "trityl": "[C](c1ccccc1)(c1ccccc1)c1ccccc1",
    "quinodimethane_benzyl": "C=C1C=CC(=C)C=C1C[CH]c1ccccc1",
}
KEEP_CLAR = {                                    # fused or bonded rings, or lone pairs
    "naphthyl_methyl": "[CH2]c1ccc2ccccc2c1",
    "naphthalene_benzyl": "c1ccc2cc([CH]c3ccccc3)ccc2c1",
    "phenoxy_benzyl": "[O]c1ccc(Cc2ccccc2)cc1",
    "biphenyl_benzyl": "[CH2]c1ccc(-c2ccccc2)cc1",
}


# --------------------------------------------------------------------------
# the ring test
# --------------------------------------------------------------------------
@pytest.mark.parametrize("smiles, shared", [
    ("CCCC", False),
    ("c1ccccc1", False),
    ("c1ccc(-c2ccccc2)cc1", False),             # biphenyl: a bridge joins the rings
    ("C" + U * 4 + "C", False),
    ("C1CCC(CC1)Cc1ccccc1", False),
    ("c1ccc2ccccc2c1", True),                    # fused
    ("C1CCC2(CC1)CCCC2", True),                  # spiro: one shared atom
    ("C1CC2CCC1C2", True),                       # bridged
    ("c1ccc2c(c1)CCC2", True),                   # aromatic fused to a saturated ring
])
def test_rings_share_an_atom(smiles, shared):
    rows = rows_of(smiles)
    assert rm._rings_share_an_atom(rm.molecule_from_rows(rows)) is shared


@pytest.mark.parametrize("smiles, skip", [
    ("C" + U * 3 + "C[CH]c1ccccc1", True),
    ("[CH2]c1ccc2ccccc2c1", False),              # fused
    ("[O]c1ccc(Cc2ccccc2)cc1", False),           # lone pairs on O
    ("Nc1ccc(Cc2ccccc2)cc1", False),             # lone pair on N
    ("c1ccc(-c2ccccc2)cc1", False),              # a bond joins the two rings
    ("C=C1C=CC(=C)C=C1Cc1ccccc1", True),
    ("CCCC", True),
])
def test_clar_cannot_matter(smiles, skip):
    assert rm._clar_cannot_matter(rm.molecule_from_rows(rows_of(smiles))) is skip


@pytest.mark.parametrize("smiles", [
    "C=C1C=CC(=C)C=C1",                           # p-quinodimethane: perceived aromatic by RMG
    "C=C1C=CC=CC1=C",                             # o-quinodimethane
    "C=C1C(=C)C(=C)C(=C)C=C1",
    "C=C1C(=C)C(=C)C(=C)C(=C)C1=C",               # hexaradialene
    "C=C1C=C[CH]C=C1",                            # benzyl radical, para form
    "C=C1C=CC=C[CH]1",                            # benzyl radical, ortho form
    "[CH]1C=CCC=C1",                              # cyclohexadienyl
    "C=C1C=C[C](C)C=C1",
    "[CH]1C=C[CH]C=C1",                           # a ring pi diradical
    "C=C1[CH]C=C([CH2])C=C1",
])
def test_every_atom_of_a_perceived_ring_has_a_double_or_benzene_bond(smiles):
    """The premise _clar_cannot_matter relies on: an atom of a ring RMG
    perceives as aromatic carries a double or 1.5 bond, so Clar's MILP can cover
    it (a radical or saturated ring carbon makes the ring non-aromatic)."""
    mol = Molecule().from_smiles(smiles)
    for ring in mol.get_aromatic_rings()[0]:
        for atom in ring:
            assert any(b.is_double() or b.is_benzene() for b in atom.bonds.values()), smiles


def test_a_quinoid_ring_cannot_be_typed_as_benzene():
    """The other premise: RMG perceives a quinodimethane ring as aromatic, but
    a carbon with two benzene bonds and a double bond has no atom type, so the
    aromatic step leaves the ring as it was, as Clar does."""
    from rmgpy.molecule.resonance import generate_aromatic_resonance_structure
    mol = Molecule().from_smiles("C=C1C=CC(=C)C=C1")
    assert len(mol.get_aromatic_rings()[0]) == 1
    assert generate_aromatic_resonance_structure(mol) == []


def test_the_ring_test_does_not_touch_the_molecule():
    mol = rm.molecule_from_rows(rows_of("C" + U * 3 + "C[CH]c1ccccc1"))
    before = rm.molecule_rows(mol), [a.id for a in mol.atoms], mol._sssr
    rm._rings_share_an_atom(mol)
    assert (rm.molecule_rows(mol), [a.id for a in mol.atoms], mol._sssr) == before


# --------------------------------------------------------------------------
# when Clar is asked for
# --------------------------------------------------------------------------
@pytest.mark.parametrize("name", sorted(UNFUSED))
def test_unfused_products_skip_clar(name, clar_requests):
    rm.ProductKeys(cache=False).forms(rows_of(UNFUSED[name]))
    rm.reactant_forms(rm.molecule_from_rows(rows_of(UNFUSED[name])))
    assert clar_requests == [False, False]


@pytest.mark.parametrize("name", sorted(KEEP_CLAR))
def test_fused_or_lone_pair_products_keep_clar(name, clar_requests):
    rm.ProductKeys(cache=False).forms(rows_of(KEEP_CLAR[name]))
    rm.reactant_forms(rm.molecule_from_rows(rows_of(KEEP_CLAR[name])))
    assert clar_requests == [True, True]


# --------------------------------------------------------------------------
# nothing changes
# --------------------------------------------------------------------------
def _inputs(rng):
    """Each product as written, plus every form RMG returns for it, each under
    ids as given, shifted, or randomly permuted."""
    for name, smiles in sorted({**UNFUSED, **KEEP_CLAR}.items()):
        base = rows_of(smiles)
        for rows in [base] + [rm.molecule_rows(f) for f in default_forms(base)]:
            ids = sorted(rows)
            mode = rng.choice(["given", "shifted", "permuted"])
            if mode == "shifted":
                rows = relabel(rows, {i: i + 1000 for i in ids})
            elif mode == "permuted":
                rows = relabel(rows, dict(zip(ids, rng.sample(range(3 * len(ids)), len(ids)))))
            yield name, mode, rows


def test_forms_equal_rmgs_default_call_in_order():
    rng = random.Random(307)
    checked = 0
    for name, mode, rows in _inputs(rng):
        want = default_forms(rows)
        got = rm.ProductKeys(cache=False).forms(rows)
        assert described(got) == described(want), (name, mode)
        assert all(a.is_identical(b) for a, b in zip(got, want)), (name, mode)
        checked += 1
    assert checked > 60


def test_reactant_forms_equal_rmgs_default_call():
    for name, smiles in sorted({**UNFUSED, **KEEP_CLAR}.items()):
        mol = rm.molecule_from_rows(rows_of(smiles))
        want = [f for f in default_forms(rm.molecule_rows(mol)) if f.reactive]
        got = rm.reactant_forms(mol)
        assert sorted(rm._signature(rm.molecule_rows(f)) for f in want) == \
            [rm._signature(rm.molecule_rows(f)) for f in got], name


def test_keys_equal_rmgs_default_call_through_one_cache():
    rng = random.Random(7307)
    inputs = list(_inputs(rng))
    rng.shuffle(inputs)
    keys = rm.ProductKeys(max_entries=64)
    distinct = set()
    for name, mode, rows in inputs:
        want = reference_key(rows)
        assert keys(rows) == want, (name, mode)
        distinct.add(want)
    assert len(distinct) >= len(UNFUSED) + len(KEEP_CLAR)
