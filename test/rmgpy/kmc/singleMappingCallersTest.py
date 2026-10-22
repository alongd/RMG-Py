"""Callers that need one isomorphism, or whether one exists, take the first
one RMG's VF2 reaches and never list the others (I-238).

``find_isomorphism`` lists every mapping before a caller takes one.  On a
polymer proxy with k interchangeable hydrogen leaves (or symmetric heavy
branches: CF3 groups, phenyl rings) that is exponential in k.  Each test
counts the mappings the caller materialises on a small symmetric case,
through a ``Molecule`` subclass whose copies and merges stay counted, and
checks the result is the one the full list gave (its first mapping, or the
first one the caller accepts).
"""

from types import SimpleNamespace

from rmgpy.kmc import compiler as compiler_module
from rmgpy.kmc import declared_proxies
from rmgpy.molecule.molecule import Molecule


class Counting(Molecule):
    """A Molecule whose copies and merges are Counting too, and which
    counts every mapping an isomorphism search hands back."""
    materialised = 0

    def copy(self, deep=False):
        return _counting(Molecule.copy(self, deep))

    def merge(self, other):
        return _counting(Molecule.merge(self, other))

    def find_isomorphism(self, other, initial_map=None, save_order=False, strict=True):
        maps = Molecule.find_isomorphism(self, other, initial_map=initial_map,
                                         save_order=save_order, strict=strict)
        Counting.materialised += len(maps)
        return maps

    def find_first_isomorphism(self, other, initial_map=None, save_order=False, strict=True,
                               accept=None):
        found = Molecule.find_first_isomorphism(self, other, initial_map=initial_map,
                                                save_order=save_order, strict=strict,
                                                accept=accept)
        Counting.materialised += found is not None
        return found


def _counting(molecule):
    counted = Counting(atoms=molecule.atoms)
    counted.multiplicity = molecule.multiplicity
    counted.reactive = molecule.reactive
    counted.is_polymer_proxy = molecule.is_polymer_proxy
    counted.metal, counted.facet = molecule.metal, molecule.facet
    return counted


def _molecule(smiles, first_id):
    molecule = _counting(Molecule().from_smiles(smiles))
    for index, atom in enumerate(molecule.atoms):
        atom.id = first_id + index
    return molecule


def test_the_explicit_h_map_takes_the_first_mapping_without_listing_the_rest():
    """Isobutane pinned on one methyl carbon: 2 * (3!)**3 = 432 full
    isomorphisms, all one heavy-atom map.  The map is the first one listed."""
    source = _molecule("CC(C)C", 0)
    target = _molecule("CC(C)C", 100)
    port_map = {0: 100}
    expected_source, expected_target = source.copy(deep=True), target.copy(deep=True)
    pins = {next(a for a in expected_source.atoms if a.id == 0):
            next(b for b in expected_target.atoms if b.id == 100)}
    listed = Molecule.find_isomorphism(expected_source, expected_target, initial_map=pins,
                                       save_order=True, strict=False)
    assert len(listed) == 432
    expected = {a.id: b.id for a, b in listed[0].items()}
    Counting.materialised = 0
    assert compiler_module._explicit_h_map(source, target, port_map) == expected
    assert Counting.materialised <= 1


def _build(smiles, first_id, ports):
    molecule = _molecule(smiles, first_id)
    return SimpleNamespace(molecule=molecule, metadata=SimpleNamespace(ports=[
        SimpleNamespace(atom_id=first_id + index, artificial=True, cap_atom_ids=())
        for index in ports]))


def test_the_same_structure_check_takes_one_mapping():
    """Two copies of 2,2,3,3-tetramethylbutane with a port on one quaternary
    carbon each: a pinned isomorphism exists (thousands of them do)."""
    first = _build("CC(C)(C)C(C)(C)C", 0, [1])
    second = _build("CC(C)(C)C(C)(C)C", 100, [1])
    Counting.materialised = 0
    assert declared_proxies._same_structure(first, second)
    assert Counting.materialised <= 1
    # a port on a methyl carbon is not the same structure as one on a
    # quaternary carbon
    Counting.materialised = 0
    assert not declared_proxies._same_structure(first, _build("CC(C)(C)C(C)(C)C", 100, [0]))
    assert Counting.materialised == 0


def test_the_transition_automorphism_check_stops_at_the_first_accepted_mapping():
    """Two perfluoroisopropyl radicals recombining, roles swapped: the swap
    is a symmetry of the transition.  The heavy-atom skeletons have
    (2 * (3!)**2)**2 = 5184 automorphisms with the roles pinned; one that
    preserves the species decides."""
    first_radical = _molecule("FC(F)(F)[C](F)C(F)(F)F", 0)
    second_radical = _molecule("FC(F)(F)[C](F)C(F)(F)F", 100)
    centre = next(a.id for a in first_radical.atoms if a.radical_electrons == 1)
    other = next(a.id for a in second_radical.atoms if a.radical_electrons == 1)
    swapped = (((centre, "*1"), (other, "*2")), ((centre, "*2"), (other, "*1")))
    Counting.materialised = 0
    assert compiler_module._transition_automorphism_exists(
        swapped[0], swapped[1], [first_radical, second_radical])
    assert Counting.materialised <= 1
