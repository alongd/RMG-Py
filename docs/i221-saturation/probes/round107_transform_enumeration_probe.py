#!/usr/bin/env python3
"""
Round 107 probe -- the fields survive the loader and die in the transforms.

    python docs/i221-saturation/probes/round107_transform_enumeration_probe.py

**HIGH.** Rounds 102 and 105 made the *loader* preserve every field a `Reaction` carries,
by discovering the field set from the class instead of writing it down. The transforms
still enumerated by hand, so the guarantee ended at the first ``pickle`` or ``copy()``:

* `TemplateReaction.__reduce__` rebuilt through a constructor that has no parameter for
  `elementary_high_p`, `allow_pdep_route` or `allow_max_rate_violation`, so all three
  came back ``False`` -- and `rank`, `comment` and `label` came back empty too;
* `TemplateReaction.copy()` dropped the same three plus `rank`, `network_kinetics` and
  `labeled_atoms`, and left `comment` as ``None`` where the class declares a ``str``;
* `LibraryReaction.__reduce__` dropped `rank`, `comment` and `label` and turned
  `is_forward` from ``True`` into ``False`` -- round 105's three fields, one transform
  over;
* `LibraryReaction` had no `copy()` at all, so it inherited `Reaction.copy`, which builds
  a **base `Reaction`**: `library`, `family` and `entry` gone, carrier and all.

`elementary_high_p` is the one with teeth. A library reaction that loses it is not
explored for pressure-dependent routing, and nothing reports that it was dropped.

**MEDIUM.** `carry_reaction_state` caught `AttributeError` from both the source read and
the target write and continued, so a field the partition called carried could be dropped
in silence -- the opposite of what its docstring promised.

``REPRODUCED`` means the defect is present. Controls must hold in both directions. Exit 1
while any finding is present, 2 if a control breaks.
"""

import pickle
import sys

from rmgpy.data.base import Entry
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import LibraryReaction, carry_reaction_state
from rmgpy.data.kinetics.library import _NOT_CARRIED_FROM_ENTRY
from rmgpy.kinetics.arrhenius import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species, TransitionState

results = []
controls = []
notes = []

A_BEFORE = 10.0


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


#: Every field this probe sets, with the value it sets, so a drop is visible by name.
#: `kinetics` is attached last on purpose: assigning `degeneracy` while kinetics are
#: attached is round 105's transformation, and a fixture that triggered it would be
#: measuring itself.
def dressed(shape):
    reactants = [Species(label="Lip", molecule=[Molecule(smiles="[Li+]")])]
    products = [Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")])]
    entry = Entry(index=4, label="an entry", long_desc="family: A_Family")
    kinetics = Arrhenius(A=(A_BEFORE, "m^3/(mol*s)"), n=0, Ea=(0, "kJ/mol"), T0=(1, "K"))

    if shape == "template":
        reaction = TemplateReaction(reactants=reactants, products=products,
                                    family="A_Family", entry=entry)
        reaction.template = ["Root"]
        reaction.estimator = "rate rules"
        reaction.labeled_atoms = {"reactants": {"*1": "Lip"}, "products": {}}
    else:
        reaction = LibraryReaction(reactants=reactants, products=products,
                                   library="a_library", entry=entry)

    reaction.elementary_high_p = True
    reaction.allow_pdep_route = True
    reaction.allow_max_rate_violation = True
    reaction.rank = 7
    reaction.comment = "a comment worth keeping"
    reaction.label = "Lip + CH3 <=> CH3Li"
    reaction.is_forward = True
    reaction.network_kinetics = Arrhenius(A=(3.0, "m^3/(mol*s)"), n=0, Ea=(0, "kJ/mol"),
                                          T0=(1, "K"))
    reaction.transition_state = TransitionState()
    reaction._degeneracy = 3.0
    reaction.kinetics = kinetics
    return reaction


#: name -> (expected value, how to read it off a transformed reaction)
EXPECTED = (
    ("elementary_high_p", True, lambda r: getattr(r, "elementary_high_p", "<missing>")),
    ("allow_pdep_route", True, lambda r: getattr(r, "allow_pdep_route", "<missing>")),
    ("allow_max_rate_violation", True,
     lambda r: getattr(r, "allow_max_rate_violation", "<missing>")),
    ("rank", 7, lambda r: getattr(r, "rank", "<missing>")),
    ("comment", "a comment worth keeping", lambda r: getattr(r, "comment", "<missing>")),
    ("label", "Lip + CH3 <=> CH3Li", lambda r: getattr(r, "label", "<missing>")),
    ("is_forward", True, lambda r: getattr(r, "is_forward", "<missing>")),
    ("degeneracy", 3.0, lambda r: getattr(r, "degeneracy", "<missing>")),
    ("network_kinetics", 3.0,
     lambda r: None if getattr(r, "network_kinetics", None) is None
     else r.network_kinetics.A.value_si),
    ("entry", "an entry",
     lambda r: None if getattr(r, "entry", None) is None else r.entry.label),
)

TEMPLATE_EXTRA = (
    ("template", ["Root"], lambda r: getattr(r, "template", "<missing>")),
    ("estimator", "rate rules", lambda r: getattr(r, "estimator", "<missing>")),
    ("labeled_atoms", ["products", "reactants"],
     lambda r: sorted(getattr(r, "labeled_atoms", {}))),
)

LIBRARY_EXTRA = (
    ("library", "a_library", lambda r: getattr(r, "library", "<missing>")),
    ("family", "a_library", lambda r: getattr(r, "family", "<missing>")),
)


def transformed(transform, reaction):
    if transform == "pickle":
        return pickle.loads(pickle.dumps(reaction))
    return reaction.copy()


def losses(shape, transform):
    reaction = dressed(shape)
    after = transformed(transform, reaction)
    expected = EXPECTED + (TEMPLATE_EXTRA if shape == "template" else LIBRARY_EXTRA)
    lost = []
    for name, want, read in expected:
        try:
            got = read(after)
        except AttributeError:
            got = "<unset>"
        if got != want:
            lost.append("{0:>26s}  {1!r} -> {2!r}".format(name, want, got))
    return after, lost


def a_transform_drops_what_the_loader_preserved(shape, transform):
    after, lost = losses(shape, transform)
    kind = "" if type(after).__name__ == type(dressed(shape)).__name__ else \
        "\n{0:>26s}  {1} -> {2}".format("(class)", shape, type(after).__name__)
    record("{0} of a {1} reaction drops fields the loader preserved".format(
        transform, shape),
        bool(lost) or bool(kind),
        ("\n".join(lost) + kind) if (lost or kind) else
        "every field set on the reaction came back with the value it was given")


def the_library_copy_returns_a_different_class():
    reaction = dressed("library")
    other = reaction.copy()
    wrong = type(other) is not LibraryReaction
    record("copy() of a library reaction returns a different class",
           wrong,
           "LibraryReaction.copy() -> {0}; `Reaction.copy` builds "
           "`Reaction.__new__(Reaction)`, so library, family and entry have nowhere to "
           "land. rmgpy/tools/isotopes.py:394 copies whatever Reaction it is handed."
           .format(type(other).__name__)
           if wrong else "LibraryReaction.copy() -> LibraryReaction")


def a_field_the_target_refuses_is_dropped_in_silence():
    """`protons` is read-only. Claim it as carried and watch the carry not complain."""
    claimed = {name: reason for name, reason in _NOT_CARRIED_FROM_ENTRY.items()
               if name != "protons"}
    try:
        carry_reaction_state(LibraryReaction(), Reaction(), claimed)
    except Exception as error:  # noqa: BLE001 -- the refusal is what we are measuring
        record("a field the target cannot take is dropped in silence", False,
               "refused with {0}: {1}".format(type(error).__name__, error))
        return
    record("a field the target cannot take is dropped in silence", True,
           "'protons' was classified as carried and the carry returned normally; the "
           "AttributeError from setattr was caught and the loop moved on")


def a_source_holding_no_state_carries_nothing_in_silence():
    class HoldsNoReactionState:
        pass

    target = LibraryReaction()
    try:
        carry_reaction_state(target, HoldsNoReactionState(), _NOT_CARRIED_FROM_ENTRY)
    except Exception as error:  # noqa: BLE001
        record("a source holding none of the state carries nothing in silence", False,
               "refused with {0}: {1}".format(type(error).__name__, error))
        return
    record("a source holding none of the state carries nothing in silence", True,
           "every getattr raised AttributeError, every one was swallowed, and the carry "
           "returned a reaction it had written nothing to")


def control_the_transform_really_happened():
    for shape in ("template", "library"):
        for transform in ("pickle", "copy"):
            reaction = dressed(shape)
            after = transformed(transform, reaction)
            if after is reaction or after.reactants is reaction.reactants:
                control("{0} of a {1} reaction produces a new object".format(
                    transform, shape), False,
                    "the transform returned the original, so every comparison above is "
                    "vacuous")
                return
    control("each transform produces a new object", True,
            "all four shape x transform combinations returned a distinct reaction with a "
            "distinct reactant list, so a surviving field really survived a round trip")


def control_the_fixture_sets_what_it_checks():
    missing = []
    for shape in ("template", "library"):
        reaction = dressed(shape)
        expected = EXPECTED + (TEMPLATE_EXTRA if shape == "template" else LIBRARY_EXTRA)
        for name, want, read in expected:
            if read(reaction) != want:
                missing.append("{0} on the {1} shape".format(name, shape))
    control("the fixture holds every field before the transform", not missing,
            "; ".join(missing) if missing else
            "both shapes carry all {0} marked fields before anything is transformed, so "
            "a loss below is the transform's".format(len(EXPECTED) + len(TEMPLATE_EXTRA)))


def control_the_entry_survives_the_template_transforms():
    """Round 95 added `entry` to both TemplateReaction transforms. It must still hold."""
    bad = []
    for transform in ("pickle", "copy"):
        after = transformed(transform, dressed("template"))
        if getattr(after, "entry", None) is None:
            bad.append(transform)
    control("round 95's carrier fix is still in place on the template shape", not bad,
            "entry is gone after: {0}".format(", ".join(bad)) if bad else
            "entry survives both template transforms, at this tip and at the base -- so "
            "this probe is measuring the fields round 102 and 105 added, not a fixture "
            "that never set them")


def control_no_transform_restates_the_rate():
    bad = []
    for shape in ("template", "library"):
        for transform in ("pickle", "copy"):
            after = transformed(transform, dressed(shape))
            if abs(after.kinetics.A.value_si - dressed(shape).kinetics.A.value_si) > 1e-12:
                bad.append("{0}/{1} A -> {2:g}".format(
                    shape, transform, after.kinetics.A.value_si))
    control("no transform restates the rate while carrying the degeneracy", not bad,
            "; ".join(bad) if bad else
            "A is unchanged in all four combinations at degeneracy 3.0 -- round 105's "
            "storage rule survives being routed through the transforms")


def control_protons_really_is_read_only():
    try:
        Reaction().protons = 1
    except AttributeError:
        control("'protons' really is read-only", True,
                "assigning it raises AttributeError, so the finding above is about a "
                "real field and not an invented one")
        return
    control("'protons' really is read-only", False,
            "it accepted an assignment, so the silent-drop finding needs another field")


def control_the_discovery_still_finds_fields():
    from rmgpy.data.kinetics.library import REACTION_STATE_FIELDS

    control("the field discovery has not gone empty", len(REACTION_STATE_FIELDS) >= 20,
            "{0} fields discovered from Reaction; an empty discovery would make every "
            "check above pass vacuously".format(len(REACTION_STATE_FIELDS)))


def note_the_census():
    note("every site in rmgpy/ that enumerates reaction state: 9, of which 4 are in "
         "this round's gates",
         "  reaction.py:322      Reaction.__reduce__        out of gates; drops "
         "allow_max_rate_violation and is_forward\n"
         "  reaction.py:2017     Reaction.copy              out of gates; drops the "
         "same two, and returns a base Reaction for every subclass that does not "
         "override it\n"
         "  family.py:__reduce__ TemplateReaction           repaired this round\n"
         "  family.py:copy       TemplateReaction           repaired this round\n"
         "  library.py:__reduce__ LibraryReaction           repaired this round\n"
         "  library.py:copy      LibraryReaction            added this round\n"
         "  library.py:_carry_entry_fields                  round 102/105's loader\n"
         "  depository.py:83     DepositoryReaction.__reduce__  out of gates; drops the "
         "three flags, rank, comment, label and network_kinetics\n"
         "  rmg/pdep.py:90       PDepReaction.__reduce__    out of gates; drops the "
         "three flags, rank and comment\n"
         "  rmg/model.py:79      as_library_reaction        round 105's conversion, "
         "already derived\n"
         "\n"
         "The review's census named three. The two it missed inside the gates are the "
         "`LibraryReaction` pair, and one of those is where round 105's own three "
         "fields were being dropped again.")


def note_the_pickle_format():
    note("the reduce tuple changed shape, deliberately",
         "`__reduce__` returns (class, (), state) now rather than (class, (17 args)). A "
         "pickle written before this tip still loads here -- no third element, no "
         "__setstate__ call, the old constructor path. A pickle written here needs the "
         "new __setstate__, because the default unpickling path would put every Reaction "
         "field in the instance dictionary where the class's descriptor shadows it. RMG "
         "pickles reactions within a run (multiprocessing, deepcopy), not across "
         "versions, so this is stated rather than worked around.")


def main():
    for shape in ("template", "library"):
        for transform in ("pickle", "copy"):
            a_transform_drops_what_the_loader_preserved(shape, transform)
    the_library_copy_returns_a_different_class()
    a_field_the_target_refuses_is_dropped_in_silence()
    a_source_holding_no_state_carries_nothing_in_silence()

    control_the_transform_really_happened()
    control_the_fixture_sets_what_it_checks()
    control_the_entry_survives_the_template_transforms()
    control_no_transform_restates_the_rate()
    control_protons_really_is_read_only()
    control_the_discovery_still_finds_fields()
    note_the_census()
    note_the_pickle_format()

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
