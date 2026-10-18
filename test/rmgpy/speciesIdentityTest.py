#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rmgpy.species_identity.same_structure: whether two molecule lists are one
species by RMG's own isomorphism.  This file is byte-identical in the
polymer-mom and kMC checkouts: both compare species with the one module.
"""

import pytest

from rmgpy.molecule import Molecule

try:
    from rmgpy.species_identity import same_structure
except ImportError:                       # before the shared module: each checkout's own copy
    try:
        from rmgpy.polymer import same_structure
    except ImportError:
        from rmgpy.kmc.mom_producer import _same_structure as same_structure


def _kekule_benzene():
    return Molecule(smiles="C1=CC=CC=C1")


def _aromatic_benzene():
    forms = _kekule_benzene().generate_resonance_structures(keep_isomorphic=False)
    aromatic = [f for f in forms if f.is_aromatic()]
    assert aromatic, "RMG gives benzene an aromatic form"
    return aromatic[0]


def _spin(smiles, multiplicity):
    molecule = Molecule(smiles=smiles)
    molecule.multiplicity = multiplicity
    return molecule


# (first, second, one species?)
PAIRS = {
    "benzene_aromatic_kekule": (_aromatic_benzene, _kekule_benzene, True),
    "singlet_diradical_two_drawings": (lambda: _spin("[CH2]C=CC[CH2]", 1), lambda: _spin("C=C[CH]C[CH2]", 1), True),
    "triplet_diradical_two_drawings": (lambda: _spin("[CH2]C=CC[CH2]", 3), lambda: _spin("C=C[CH]C[CH2]", 3), True),
    "singlet_against_triplet": (lambda: _spin("[CH2]C=CC[CH2]", 1), lambda: _spin("C=C[CH]C[CH2]", 3), False),
    "singlet_against_triplet_one_drawing": (lambda: _spin("[CH2]C[CH]C", 1), lambda: _spin("[CH2]C[CH]C", 3), False),
    "tautomers": (lambda: Molecule(smiles="Oc1ccccn1"), lambda: Molecule(smiles="O=c1cccc[nH]1"), False),
    "crotyl_two_drawings": (lambda: Molecule(smiles="[CH2]C=CC"), lambda: Molecule(smiles="C=C[CH]C"), True),
    "one_drawing": (lambda: Molecule(smiles="CCO"), lambda: Molecule(smiles="OCC"), True),
}


@pytest.mark.parametrize("name", sorted(PAIRS))
def test_same_structure_in_both_orders(name):
    """Every pinned pair: the verdict, and the same verdict in both orders."""
    first, second, expected = PAIRS[name]
    assert same_structure([first()], [second()]) is expected
    assert same_structure([second()], [first()]) is expected


def test_the_inputs_are_left_as_they_were():
    """Resonance forms are generated on copies: the inputs keep their bonds
    and multiplicity."""
    a, b = _spin("[CH2]C=CC[CH2]", 1), _spin("C=C[CH]C[CH2]", 1)
    before = (a.to_adjacency_list(), b.to_adjacency_list())
    same_structure([a], [b])
    assert (a.to_adjacency_list(), b.to_adjacency_list()) == before and a.multiplicity == b.multiplicity == 1


def test_nothing_equals_an_empty_list():
    assert same_structure([], [Molecule(smiles="C")]) is False
    assert same_structure([Molecule(smiles="C")], []) is False


def test_this_checkout_compares_species_by_the_one_implementation():
    """MoM's signature and the producer's mask check both call the shared
    function, whichever of them this checkout carries."""
    from rmgpy.species_identity import same_structure as shared
    callers = []
    try:
        from rmgpy.polymer import same_structure as mom
        callers.append(mom)
    except ImportError:
        pass
    try:
        from rmgpy.kmc.mom_producer import _same_structure as producer
        callers.append(producer)
    except ImportError:
        pass
    assert callers and all(c is shared for c in callers)


def _distinct(forms):
    kept = []
    for form in forms:
        if not any(form.is_isomorphic(k, strict=True) for k in kept):
            kept.append(form)
    return kept


def test_the_resonance_forms_compared_are_unfiltered():
    """Default filtering drops octet-violating and charge-separated forms
    (nitromethane's [O-][N+2](C)[O-]); the comparison keeps every form."""
    from rmgpy.species_identity import _forms
    nitromethane = Molecule(smiles="C[N+](=O)[O-]")
    filtered = _distinct(nitromethane.copy(deep=True).generate_resonance_structures(keep_isomorphic=False))
    compared = _distinct(_forms(nitromethane))
    assert len(compared) > len(filtered)
    assert all(any(f.is_isomorphic(c, strict=True) for c in compared) for f in filtered)


def test_the_resonance_forms_keep_the_spin_of_the_input():
    """RMG resets a generated form's multiplicity from its radical count; a
    singlet diradical's forms stay singlets, a triplet's stay triplets."""
    from rmgpy.species_identity import _forms
    for multiplicity in (1, 3):
        forms = _forms(_spin("[CH2]C=CC[CH2]", multiplicity))
        assert len(forms) > 1 and all(f.multiplicity == multiplicity for f in forms)
