#!/usr/bin/env python3
"""
Round 105 probe -- a carry that rescales, an enumeration taken from the wrong class, and a
containment flag that can evaporate.

    python docs/i221-saturation/probes/round105_degeneracy_and_enumeration_probe.py

**HIGH.** `_carry_entry_fields` assigns `degeneracy` through the property, and
`Reaction.degeneracy`'s setter is not an assignment: when kinetics are already attached it
multiplies the rate by a ratio and edits the kinetics comment. The kinetics object at that
moment is ``entry.data`` itself -- the shared database object -- so the damage is not
confined to the reaction being built.

The ratio has two branches (`reaction.py:359`), and which one fires depends on the value
the constructor was given:

* the template shape passes ``degeneracy=`` to the constructor, so the old value equals
  the new one and the ratio is 1 for every degeneracy >= 2 -- which is why a probe using
  ``degeneracy=3`` sees nothing;
* the two `LibraryReaction` shapes pass no ``degeneracy`` at all, so the old value is 1,
  the ``< 2`` branch fires, and the ratio is the **whole new degeneracy**.

**MEDIUM.** `as_library_reaction` enumerates `LibraryReaction.__init__`'s parameters, so
state the class holds but the constructor does not take is dropped by construction.

**MEDIUM.** `_read_manifest` takes `O_NOFOLLOW` with `getattr(os, 'O_NOFOLLOW', 0)`. Where
the constant is absent the flag becomes a no-op and every component below the anchor is
followed -- the containment is gone and nothing says so.

``REPRODUCED`` means the defect is present. Controls must hold in both directions. Exit 1
while any finding is present, 2 if a control breaks.
"""

import os
import shutil
import sys
import tempfile

from rmgpy import settings
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import KineticsLibrary, LibraryReaction
from rmgpy.kinetics.arrhenius import Arrhenius, Marcus
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import as_library_reaction
from rmgpy.species import Species
from rmgpy.data.kinetics.quarantine import QUARANTINE_FILENAME, resolve_quarantine

MANIFEST = """
name = "{0}/quarantine"
state = "QUARANTINED FOR TESTING"
appliesToKineticsClass = "Marcus"
reason = "{1}"
"""

results = []
controls = []
notes = []


def record(name, reproduced, detail):
    results.append((name, reproduced, detail))
    print("\n[{0}] {1}\n    {2}".format(
        "REPRODUCED" if reproduced else "NOT REPRODUCED", name,
        str(detail).replace("\n", "\n    ")))


def control(name, ok, detail):
    controls.append((name, ok, detail))
    print("\n[{0}] control: {1}\n    {2}".format(
        "HOLDS" if ok else "BROKEN", name, str(detail).replace("\n", "\n    ")))


def note(name, detail):
    notes.append((name, detail))
    print("\n[NOTED, not a finding] {0}\n    {1}".format(
        name, str(detail).replace("\n", "\n    ")))


def clear_caches():
    import rmgpy.data.kinetics.quarantine as q
    for name in ("_DISK_QUARANTINE_CACHE", "_DISK_ANY_QUARANTINE_CACHE",
                 "_UNATTRIBUTED_WARNED", "_UNANSWERED_WARNED", "_UNSAFE_LABELS_WARNED",
                 "_UNSAFE_MANIFESTS_WARNED", "_LEGACY_CALL_SITES_WARNED"):
        marked = getattr(q, name, None)
        if marked is not None:
            marked.clear()


class pinned_database(object):
    def __init__(self, root):
        self.root = root

    def __enter__(self):
        self.previous = settings.get("database.directory")
        settings["database.directory"] = self.root
        clear_caches()
        return self.root

    def __exit__(self, *exc):
        settings["database.directory"] = self.previous
        clear_caches()
        return False


#: The rate every sweep below starts from. A round number so a ratio is unmistakable.
A_BEFORE = 10.0


def arrhenius(comment="Estimated from node Root"):
    return Arrhenius(A=(A_BEFORE, "m^3/(mol*s)"), n=0, Ea=(0, "kJ/mol"), T0=(1, "K"),
                     comment=comment)


def library_declaring(shape, degeneracy, data=None):
    """
    A library whose single entry declares `degeneracy` ON ``entry.item``.

    `shape` selects which of the three branches of `get_library_reactions` runs:
    ``"template"`` (auto-generated, rate rule), ``"originally_from"`` (auto-generated,
    re-exported library) or ``"ordinary"``.
    """
    if shape == "template":
        long_desc = "\n".join([
            "Matched reaction 3 Lip + CH3 <=> CH3Li in A_Family/rate rule [Root]",
            "Euclidian distance = 0",
            "family: A_Family",
        ])
    elif shape == "originally_from":
        long_desc = "Originally from reaction library: some_other_library"
    else:
        long_desc = ""
    library = KineticsLibrary(label="a_seed", name="a_seed")
    library.auto_generated = shape != "ordinary"
    library.entries = {
        1: Entry(
            index=1, label="Lip + CH3 <=> CH3Li",
            item=Reaction(
                reactants=[Species(label="Lip", molecule=[Molecule(smiles="[Li+]")],
                                   reactive=False),
                           Species(label="CH3", molecule=[Molecule(smiles="[CH3]")],
                                   reactive=False)],
                products=[Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")],
                                  reactive=False)],
                reversible=False, electrons=1, degeneracy=degeneracy),
            data=arrhenius() if data is None else data,
            long_desc=long_desc,
        )
    }
    return library


SHAPES = ("template", "originally_from", "ordinary")
SWEEP = (0.5, 1.0, 1.5, 2.0, 3.0)


# ---------------------------------------------------------------------------
# HIGH -- carrying `degeneracy` rescales the rate and mutates the database object
# ---------------------------------------------------------------------------


def carrying_degeneracy_rescales_the_rate():
    damage = []
    for shape in SHAPES:
        for degeneracy in SWEEP:
            library = library_declaring(shape, degeneracy)
            entry = library.entries[1]
            reaction = library.get_library_reactions()[0]
            after = reaction.kinetics.A.value_si
            shared = entry.data.A.value_si
            if abs(after - A_BEFORE) > 1e-9 or abs(shared - A_BEFORE) > 1e-9:
                damage.append("{0:>16s} degeneracy={1:<4g} A {2:g} -> {3:g}"
                              "  (x{4:g}), entry.data.A -> {5:g}".format(
                                  shape, degeneracy, A_BEFORE, after,
                                  after / A_BEFORE, shared))
    record("carrying degeneracy rescales the rate and the shared entry",
           bool(damage),
           "\n".join(damage) if damage else
           "every shape x degeneracy left A at {0:g} on both the reaction and the "
           "entry".format(A_BEFORE))


def carrying_degeneracy_edits_the_shared_comment():
    """
    The second half of the setter, which fires even where the ratio is exactly 1.

    `entry.data` is the library's own kinetics object. This campaign reads authorship out
    of a kinetics comment (`authoring_families`), so appending to it is not cosmetic.
    """
    damage = []
    for shape in SHAPES:
        for degeneracy in SWEEP:
            library = library_declaring(shape, degeneracy)
            entry = library.entries[1]
            before = entry.data.comment
            library.get_library_reactions()
            if entry.data.comment != before:
                damage.append("{0:>16s} degeneracy={1:<4g} entry.data.comment {2!r} -> "
                              "{3!r}".format(shape, degeneracy, before,
                                             entry.data.comment))
    record("carrying degeneracy appends to the shared kinetics comment",
           bool(damage),
           "\n".join(damage) if damage else
           "no shape x degeneracy edited the comment on the database object")


# ---------------------------------------------------------------------------
# MEDIUM -- the conversion enumerates the destination constructor
# ---------------------------------------------------------------------------


def the_conversion_enumerates_the_destination_constructor():
    source = TemplateReaction(
        reactants=[Species(label="Lip", molecule=[Molecule(smiles="[Li+]")])],
        products=[Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")])],
        kinetics=arrhenius(comment="a rate"), family="A_Family", template="Root",
        degeneracy=1, electrons=1)
    source.is_forward = True
    source.rank = 7
    source.comment = "kept?"
    source.label = "Lip + CH3 <=> CH3Li"

    converted = as_library_reaction(source, "a_library")
    lost = []
    for name in ("is_forward", "rank", "comment", "label"):
        was = getattr(source, name)
        now = getattr(converted, name, "<absent>")
        if was != now:
            lost.append("{0:<12s} {1!r} -> {2!r}".format(name, was, now))
    record("the conversion drops state the destination constructor does not take",
           bool(lost),
           "\n".join(lost) if lost else
           "every field checked survived the conversion")


def every_carried_field_actually_arrives():
    """
    A carry that cannot carry a field is a drop, whatever the partition says.

    `_carry_entry_fields` swallows `AttributeError` from *both* sides, so a field the
    target refuses to accept is classified "carried" and silently is not.
    """
    from rmgpy.data.kinetics import library as library_module

    fields = getattr(library_module, "REACTION_STATE_FIELDS", None)
    if fields is None:
        fields = library_module._REACTION_FIELDS
    excluded = library_module._NOT_CARRIED_FROM_ENTRY
    refused = []
    probe = Reaction()
    for name in sorted(set(fields) - set(excluded)):
        try:
            setattr(probe, name, getattr(probe, name))
        except AttributeError as error:
            refused.append("{0:<12s} {1}".format(name, error))
    record("a field the partition calls carried cannot be assigned at all",
           bool(refused),
           "\n".join(refused) if refused else
           "every field the partition calls carried accepts an assignment")


# ---------------------------------------------------------------------------
# MEDIUM -- containment degrades silently when `O_NOFOLLOW` is absent
# ---------------------------------------------------------------------------


def escape_at(tmpdir, name, component):
    """A family reached through a symlinked `component`, exactly as the tests build it."""
    label = "A_Family_Reached_Through_A_Linked_Parent"
    base = os.path.join(tmpdir, name)
    outside = os.path.join(base, "outside_the_database")
    os.makedirs(os.path.join(outside, label))
    with open(os.path.join(outside, label, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST.format(label, "reached from outside the database"))
    root = os.path.join(base, "database")
    os.makedirs(os.path.join(root, "kinetics"))
    os.symlink(outside, os.path.join(root, "kinetics", component))
    return root, label


class without_o_nofollow(object):
    """Delete the constant, as a platform that never defined it would have it."""

    def __enter__(self):
        self.value = getattr(os, "O_NOFOLLOW", None)
        if self.value is not None:
            del os.O_NOFOLLOW
        return self

    def __exit__(self, *exc):
        if self.value is not None:
            os.O_NOFOLLOW = self.value
        return False


def a_missing_o_nofollow_follows_every_link(tmpdir):
    root, label = escape_at(tmpdir, "nofollow_absent", "families")
    with pinned_database(root):
        with without_o_nofollow():
            answer = resolve_quarantine(label)
    record("without O_NOFOLLOW the descent follows links out of the database",
           answer != (None, False),
           "resolve_quarantine({0!r}) = {1!r}; the family directory was reached through "
           "a symlinked `kinetics/families` and the manifest behind it was executed, "
           "because the flag degraded to 0 and nothing refused".format(label, answer))


# ---------------------------------------------------------------------------
# controls
# ---------------------------------------------------------------------------


def control_the_degeneracy_value_is_forwarded():
    forwarded = []
    for shape in SHAPES:
        for degeneracy in SWEEP:
            library = library_declaring(shape, degeneracy)
            reaction = library.get_library_reactions()[0]
            if abs(reaction.degeneracy - degeneracy) > 1e-9:
                forwarded.append("{0} degeneracy={1:g} -> {2:g}".format(
                    shape, degeneracy, reaction.degeneracy))
    control("the degeneracy value itself reaches the reaction", not forwarded,
            "\n".join(forwarded) if forwarded else
            "all {0} shape x degeneracy combinations forwarded the scalar".format(
                len(SHAPES) * len(SWEEP)))


def control_a_marcus_rate_is_untouched_at_degeneracy_one():
    library = library_declaring("template", 1.0, data=Marcus(
        A=(1.73e06, "m^3/(mol*s)"), n=2,
        lmbd_i_coefs=[21824.5, -0.0341626, -0.0013254, 4.92966e-07],
        beta=(1.2e10, "1/m"), wr=(0, "kJ/mol"), wp=(0, "kJ/mol"),
        lmbd_o=(0, "J/mol"), comment="Estimated from node Root"))
    entry = library.entries[1]
    before = list(entry.data.lmbd_i_coefs.value_si) + [entry.data.A.value_si]
    library.get_library_reactions()
    after = list(entry.data.lmbd_i_coefs.value_si) + [entry.data.A.value_si]
    control("a degeneracy of 1 leaves the Marcus fixture alone", before == after,
            "{0} -> {1}".format(before, after))


def control_the_conversion_still_replaces_the_library():
    source = TemplateReaction(
        reactants=[Species(label="Lip", molecule=[Molecule(smiles="[Li+]")])],
        products=[Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")])],
        kinetics=arrhenius(), family="A_Family", template="Root", electrons=1)
    converted = as_library_reaction(source, "a_library")
    control("the conversion still replaces library and family",
            isinstance(converted, LibraryReaction)
            and converted.library == "a_library" and converted.family == "a_library"
            and converted.electrons == 1,
            "library={0!r} family={1!r} electrons={2!r}".format(
                converted.library, converted.family, converted.electrons))


def control_the_conversion_does_not_share_the_reactant_list():
    source = TemplateReaction(
        reactants=[Species(label="Lip", molecule=[Molecule(smiles="[Li+]")])],
        products=[Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")])],
        kinetics=arrhenius(), family="A_Family", template="Root")
    converted = as_library_reaction(source, "a_library")
    control("the conversion still copies the reactant and product lists",
            converted.reactants is not source.reactants
            and converted.products is not source.products,
            "reactants copied={0} products copied={1}".format(
                converted.reactants is not source.reactants,
                converted.products is not source.products))


def control_the_escape_is_refused_with_o_nofollow(tmpdir):
    root, label = escape_at(tmpdir, "nofollow_present", "families")
    with pinned_database(root):
        answer = resolve_quarantine(label)
    control("with O_NOFOLLOW the same escape is refused", answer == (None, False),
            "resolve_quarantine({0!r}) = {1!r}".format(label, answer))


def control_an_in_tree_family_still_answers(tmpdir):
    label = "A_Quarantined_Family"
    root = os.path.join(tmpdir, "in_tree")
    family = os.path.join(root, "kinetics", "families", label)
    os.makedirs(family)
    with open(os.path.join(family, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST.format(label, "an ordinary in-tree manifest"))
    with pinned_database(root):
        quarantine, answered = resolve_quarantine(label)
    control("an ordinary in-tree family is still answered",
            answered and quarantine is not None,
            "answered={0} quarantine={1}".format(answered, bool(quarantine)))


def note_the_ratio_branch():
    note("the `< 2` branch is why a degeneracy of 3 looks harmless",
         "reaction.py:359 uses `new` as the ratio when the OLD degeneracy is below 2, "
         "and `new / old` otherwise. Re-assigning the value a reaction already holds is "
         "therefore idempotent only for old >= 2: at old == 1 (the two LibraryReaction "
         "shapes, which pass no degeneracy to the constructor) the rate is multiplied by "
         "the whole new degeneracy, and at old == 0.5 by the new value again.")


def note_the_ordering_rule_and_the_census():
    """
    The property that decides every other site, measured rather than read.

    Assigning `degeneracy` is safe exactly while `kinetics` is None. Every carry in the
    engine is safe or unsafe by that one positional fact, and none of them says so.
    """
    before = Reaction(degeneracy=1)
    before.degeneracy = 4.0
    before.kinetics = arrhenius(comment="")
    after = Reaction(degeneracy=1, kinetics=arrhenius(comment=""))
    after.degeneracy = 4.0
    note("the ordering rule: assign degeneracy before the kinetics, or it is a rescale",
         "degeneracy assigned BEFORE attaching kinetics -> A = {0:g}\n"
         "degeneracy assigned AFTER  attaching kinetics -> A = {1:g}\n"
         "\n"
         "Census of `.degeneracy = ` in rmgpy/ (9 sites):\n"
         "  species.py:1018            not a Reaction -- TransitionState.__init__\n"
         "  common.py:397, :411        an application at generation time, before any\n"
         "                             kinetics exist; this is what the setter is for\n"
         "  chemkin.pyx:895            an application that then undoes the rescale by\n"
         "                             hand with change_rate(1/degen) -- the same repair\n"
         "                             as this round's, done lossily, and its comment is\n"
         "                             still appended\n"
         "  reaction.py:2033           a carry in Reaction.copy(), SAFE: assigned two\n"
         "                             lines before `other.kinetics` is attached\n"
         "  family.py:210              a carry in TemplateReaction.copy(), SAFE, same\n"
         "                             ordering\n"
         "  model.py:1152              a carry from the reverse reaction, SAFE: four\n"
         "                             lines before `reaction.kinetics = kinetics`\n"
         "  isotopes.py:452            a carry, SAFE: before `rxn.kinetics = ...`\n"
         "  isotopes.py:458            a carry, **UNSAFE**: three lines later, after the\n"
         "                             kinetics were attached, so the isotopologue's rate\n"
         "                             is multiplied by the reverse degeneracy. Out of\n"
         "                             this round's gates; named, not fixed\n"
         "  library.py                 the carry this round repairs\n"
         "\n"
         "Four of the five carries are safe only because of where the line sits, and "
         "nothing marks it.".format(before.kinetics.A.value_si, after.kinetics.A.value_si))


def main():
    tmpdir = tempfile.mkdtemp(prefix="round105-probe-")
    try:
        carrying_degeneracy_rescales_the_rate()
        carrying_degeneracy_edits_the_shared_comment()
        the_conversion_enumerates_the_destination_constructor()
        every_carried_field_actually_arrives()
        a_missing_o_nofollow_follows_every_link(tmpdir)

        control_the_degeneracy_value_is_forwarded()
        control_a_marcus_rate_is_untouched_at_degeneracy_one()
        control_the_conversion_still_replaces_the_library()
        control_the_conversion_does_not_share_the_reactant_list()
        control_the_escape_is_refused_with_o_nofollow(tmpdir)
        control_an_in_tree_family_still_answers(tmpdir)
        note_the_ratio_branch()
        note_the_ordering_rule_and_the_census()
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

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
