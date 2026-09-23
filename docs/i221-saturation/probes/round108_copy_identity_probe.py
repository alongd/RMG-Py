#!/usr/bin/env python3
"""
Round 108 probe -- copying the reaction component by component severs its internal
references.

    python docs/i221-saturation/probes/round108_copy_identity_probe.py

**HIGH.** A reaction's `pairs` hold the reaction's *own* `Species` objects --
`Reaction.generate_pairs` appends them straight out of ``self.reactants`` and
``self.products`` -- and its `labeled_atoms` hold `Atom` objects out of those species'
molecules (``family.py``, ``_create_reaction``). Its internal consistency is a property of
the references *between* those lists, not of the lists separately.

`copy()` deep-copied `reactants`, `products`, `pairs` and `labeled_atoms` in four separate
calls, each with its own memo, so on the other side:

* a `pairs` entry is a species that is in no list of the copy. `Species.__eq__` is
  identity, so ``reactants.index(pair[0])`` raises;
* a labelled atom is inside no molecule the copy owns. Relabelling through it --
  ``family.py:2593``, immediately before pairs and templates are regenerated -- changes
  nothing and reports nothing.

Every check here is `is`. Value equality cannot see this: a detached clone has all the
right values, which is why the tests that compare by value were green throughout.

``REPRODUCED`` means the defect is present. Controls must hold in both directions. Exit 1
while any finding is present, 2 if a control breaks.
"""

import pickle
import sys
from copy import deepcopy

from rmgpy.data.base import Entry
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.molecule import Molecule
from rmgpy.species import Species

results = []
controls = []
notes = []


def record(name, reproduced, detail):
    results.append((name, reproduced, detail))
    print("\n{0}  {1}".format("REPRODUCED " if reproduced else "absent     ", name))
    print("    " + detail.replace("\n", "\n    "))


def control(name, ok, detail):
    controls.append((name, ok, detail))
    print("\n{0}  control: {1}".format("holds      " if ok else "BROKEN     ", name))
    print("    " + detail.replace("\n", "\n    "))


def note(name, detail):
    notes.append((name, detail))
    print("\nnote        {0}".format(name))
    print("    " + detail.replace("\n", "\n    "))


def built(shape):
    """A reaction whose `pairs` and `labeled_atoms` point into its own species."""
    reactants = [Species(label="ethane", molecule=[Molecule(smiles="CC")])]
    products = [Species(label="ethyl", molecule=[Molecule(smiles="C[CH2]")])]
    if shape == "template":
        reaction = TemplateReaction(reactants=reactants, products=products,
                                    family="A_Family", entry=Entry(index=4, label="e"))
    else:
        reaction = LibraryReaction(reactants=reactants, products=products,
                                   library="a_library", entry=Entry(index=4, label="e"))
    reaction.pairs = [(reactants[0], products[0])]
    labelled = reactants[0].molecule[0].atoms[0]
    labelled.label = "*1"
    if shape == "template":
        reaction.labeled_atoms = {"reactants": {"*1": labelled}, "products": {}}
    return reaction


def transformed(transform, reaction):
    if transform == "pickle":
        return pickle.loads(pickle.dumps(reaction))
    return reaction.copy()


def atoms_of(reaction):
    return [atom
            for species in list(reaction.reactants) + list(reaction.products)
            for molecule in species.molecule
            for atom in molecule.atoms]


def species_of(reaction):
    return list(reaction.reactants) + list(reaction.products)


def the_pairs_point_outside_the_copy(shape, transform):
    reaction = built(shape)
    after = transformed(transform, reaction)
    owned, original = species_of(after), species_of(reaction)
    stray = [(i, j) for i, pair in enumerate(after.pairs)
             for j, member in enumerate(pair)
             if not any(member is s for s in owned)]
    still_original = any(member is s for pair in after.pairs for member in pair
                         for s in original)
    record("{0} of a {1} reaction leaves pairs pointing outside it".format(
        transform, shape),
        bool(stray),
        "pairs entries not owned by the reaction: {0}; and not the original's either: "
        "{1}. Species.__eq__ is identity, so reactants.index(pair[0]) raises here."
        .format(stray, not still_original) if stray else
        "every pairs member is one of the reaction's own species, and none is the "
        "original's")


def the_labelled_atoms_point_outside_the_copy(transform):
    reaction = built("template")
    after = transformed(transform, reaction)
    owned = atoms_of(after)
    stray = [label for group in after.labeled_atoms.values()
             for label, atom in group.items()
             if not any(atom is a for a in owned)]
    record("{0} leaves labeled_atoms pointing outside the reaction".format(transform),
           bool(stray),
           "labels whose atom is inside no molecule the reaction owns: {0}".format(stray)
           if stray else
           "every labelled atom is an atom of one of the reaction's own molecules")


def relabelling_the_copy_changes_nothing(transform):
    reaction = built("template")
    after = transformed(transform, reaction)
    after.labeled_atoms["reactants"]["*1"].label = "*9"
    visible = any(a.label == "*9" for a in atoms_of(after))
    leaked = any(a.label == "*9" for a in atoms_of(reaction))
    record("relabelling through labeled_atoms after a {0} is a silent no-op".format(
        transform),
        not visible,
        "the label was set and no atom of the reaction carries it; family.py:2593 does "
        "exactly this before regenerating pairs and templates, and would regenerate them "
        "from unlabelled structures" if not visible else
        "the new label is visible on an atom the reaction owns, and did not leak back "
        "into the original ({0})".format("leaked!" if leaked else "no leak"))


def control_the_fixture_holds_the_invariant():
    bad = []
    for shape in ("template", "library"):
        reaction = built(shape)
        if not any(reaction.pairs[0][0] is s for s in reaction.reactants):
            bad.append("{0}: pairs".format(shape))
        if shape == "template" and not any(
                reaction.labeled_atoms["reactants"]["*1"] is a
                for a in atoms_of(reaction)):
            bad.append("{0}: labeled_atoms".format(shape))
    control("the fixture points into its own species before any copy", not bad,
            "; ".join(bad) if bad else
            "pairs and labeled_atoms are the reaction's own objects to begin with, so a "
            "severed reference below is the transform's doing and not the fixture's")


def control_the_copy_is_a_copy():
    bad = []
    for shape in ("template", "library"):
        for transform in ("copy", "pickle"):
            reaction = built(shape)
            after = transformed(transform, reaction)
            if any(a is b for a in species_of(after) for b in species_of(reaction)):
                bad.append("{0}/{1}".format(shape, transform))
    control("each transform really copies the species", not bad,
            "; ".join(bad) if bad else
            "no species of any transformed reaction is an object of the original, so "
            "'points into the copy' is a real property and not aliasing")


def control_the_shallow_half_is_still_shallow():
    bad = []
    for shape in ("template", "library"):
        reaction = built(shape)
        if reaction.copy().entry is not reaction.entry:
            bad.append(shape)
    control("copy() still shares the entry by reference", not bad,
            "; ".join(bad) if bad else
            "entry is the shared database object the gate reads authorship from, and "
            "copy() still hands the copy the same one -- the repair deepened nothing new")


def control_value_equality_cannot_see_any_of_this():
    """The second half of the finding, as a control: the old assertions still pass."""
    reaction = built("template")
    after = reaction.copy()
    by_value = ([tuple(s.label for s in pair) for pair in after.pairs]
                == [tuple(s.label for s in pair) for pair in reaction.pairs]
                and sorted(after.labeled_atoms) == sorted(reaction.labeled_atoms))
    control("value comparison passes whether or not the references survive", by_value,
            "pairs compare equal by label and labeled_atoms by key, at this tip and at "
            "the base alike. A test written that way is structurally unable to fail on "
            "this defect, which is why round 107's tests were green throughout.")


def note_why_pickle_and_not_deepcopy():
    molecule = Molecule(smiles="CC")
    atom = molecule.atoms[0]
    memo = {}
    copied_atom = deepcopy(atom, memo)
    copied_molecule = deepcopy(molecule, memo)
    shared = any(copied_atom is a for a in copied_molecule.atoms)
    note("one shared memo fixes pairs and cannot fix labeled_atoms",
         "deepcopy through one memo shares the atom: {0}\n"
         "`Molecule.__deepcopy__` (rmgpy/molecule/molecule.py:1064) is\n"
         "    return self.copy(deep=True)\n"
         "-- it takes the memo and discards it, so a caller cannot make two references "
         "to one molecule come out as two references to one copy. Species-level sharing "
         "survives because Species has no such override; atom-level sharing cannot. "
         "pickle keeps its own memo, which __deepcopy__ never sees, and is 3x faster "
         "than deepcopy on the same state besides. rmgpy/molecule/ is out of gates."
         .format(shared))


def note_the_same_defect_out_of_gates():
    note("the sites with the same defect that this round cannot touch",
         "  reaction.py:2017   Reaction.copy          `other.pairs = deepcopy(self.pairs)` "
         "beside separately copied reactants -- the same severing, for every subclass "
         "that does not override copy()\n"
         "  family.py:3854     deepcopy(rxn)          deepcopy consults "
         "Molecule.__deepcopy__, so labelled atoms are severed there too\n"
         "  database.py:755    deepcopy(reaction)     the same\n"
         "\n"
         "pickle is unaffected at every one of them, which is why `copy()` now uses it.")


def main():
    for shape in ("template", "library"):
        for transform in ("copy", "pickle"):
            the_pairs_point_outside_the_copy(shape, transform)
    for transform in ("copy", "pickle"):
        the_labelled_atoms_point_outside_the_copy(transform)
        relabelling_the_copy_changes_nothing(transform)

    control_the_fixture_holds_the_invariant()
    control_the_copy_is_a_copy()
    control_the_shallow_half_is_still_shallow()
    control_value_equality_cannot_see_any_of_this()
    note_why_pickle_and_not_deepcopy()
    note_the_same_defect_out_of_gates()

    reproduced = [name for name, flag, _ in results if flag]
    broken = [name for name, ok, _ in controls if not ok]
    print("\n" + "=" * 78)
    print("findings reproduced : {0} of {1}".format(len(reproduced), len(results)))
    print("controls holding    : {0} of {1}".format(
        len(controls) - len(broken), len(controls)))
    print("noted, not findings : {0}".format(len(notes)))
    for name in reproduced:
        print("  REPRODUCED  {0}".format(name))
    for name in broken:
        print("  CONTROL BROKEN  {0}".format(name))
    if broken:
        return 2
    return 1 if reproduced else 0


if __name__ == "__main__":
    sys.exit(main())
