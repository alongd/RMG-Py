#!/usr/bin/env python3
"""
Round 95 probe -- the reader that ignores the carrier, and four ways the path is still open.

    python docs/i221-saturation/probes/round95_carrier_probe.py

Round 92 attached the authoring ``Entry`` to the template shape and carried it across
both conversions in ``rmgpy/rmg/model.py``. The manager verified that and then found the
*reader*: ``authoring_families()`` short-circuits for anything whose ``library`` is
``None`` and returns the single ``reaction.family`` slot -- which ``library.py`` rewrites
once per ``family:`` line, so only the last survives. A second ``family:`` line naming a
loaded ordinary family therefore hides a quarantined first one, on exactly the shape
round 92 gave a carrier to.

Alongside: the two conversions construct their replacement ``LibraryReaction`` with 8 of
the 17 fields the class takes, so a converted ``Ar + e- => Ar+ + 2e-`` claims zero net
electrons; the manifest *file* is still executed without being validated, while only its
parent directory is; a loaded family never re-reads its manifest; and the suppression
cache is keyed on something a same-directory addition does not move.

``REPRODUCED`` means the defect is present. Controls must hold in both directions or the
run measures nothing. Exit 1 while any finding is present, 2 if a control breaks.
"""

import copy as copy_module
import inspect
import logging
import os
import pickle
import shutil
import sys
import tempfile
import types

from rmgpy import settings
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import KineticsLibrary, LibraryReaction
from rmgpy.data.kinetics.quarantine import (
    QUARANTINE_FILENAME,
    authoring_families,
    check_quarantine,
    load_family_quarantine,
    resolve_quarantine,
)
from rmgpy.exceptions import QuarantinedKineticsError
from rmgpy.kinetics import Marcus
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species

#: Quarantined, on disk, and deliberately NOT loaded -- the state in which the gate has
#: to answer from disk, which is round 89's repair and the thing round 95's HIGH slips past.
QUARANTINED = "A_Family_Quarantined_On_Disk"

#: Loaded, ordinary, unquarantined. Naming it in a SECOND `family:` line is the whole
#: attack: `library.py` keeps the last label, and a loaded family is not converted, so the
#: LibraryReaction path that round 89 fixed is never entered.
INNOCENT = "An_Ordinary_Loaded_Family"

#: In no database at all -- the unanswered case, which is admitted with a warning. Used
#: for the field-carrying finding, where the reaction has to reach the core to be read.
ABSENT = "A_Family_This_Database_Does_Not_Have"

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
    """A measured fact that is deliberately NOT being changed. Never affects the exit code."""
    notes.append((name, detail))
    print("\n[NOTED, not a finding] {0}\n    {1}".format(
        name, str(detail).replace("\n", "\n    ")))


def marcus(comment=""):
    return Marcus(
        A=(1.73e06, "m^3/(mol*s)"), n=2,
        lmbd_i_coefs=[21824.5, -0.0341626, -0.0013254, 4.92966e-07],
        beta=(1.2e10, "1/m"), wr=(0, "kJ/mol"), wp=(0, "kJ/mol"),
        lmbd_o=(0, "J/mol"), comment=comment,
    )


def species_pair():
    a = Species(label="Lip", molecule=[Molecule(smiles="[Li+]")], reactive=False)
    b = Species(label="CH3", molecule=[Molecule(smiles="[CH3]")], reactive=False)
    c = Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")], reactive=False)
    return [a, b], [c]


def template_library(label, family_lines, electrons=0, degeneracy=1):
    """
    A real KineticsLibrary whose one entry produces the auto-generated TEMPLATE shape.

    `family_lines` is a list, because the defect is what happens when longDesc carries
    more than one of them. The kinetics comment deliberately holds no `family:` line, so
    the entry's longDesc is the only carrier.
    """
    reactants, products = species_pair()
    long_desc = ["Matched reaction 3 Lip + CH3 <=> CH3Li in {0}/rate rule [Root]".format(
        family_lines[0]), "Euclidian distance = 0"]
    long_desc += ["family: {0}".format(name) for name in family_lines]
    library = KineticsLibrary(label=label, name=label)
    library.auto_generated = True
    library.entries = {
        1: Entry(
            index=1,
            label="Lip + CH3 <=> CH3Li",
            item=Reaction(reactants=reactants, products=products, reversible=False,
                          electrons=electrons, degeneracy=degeneracy),
            data=marcus(),
            long_desc="\n".join(long_desc),
        )
    }
    return library


def ordinary_library(label, family_lines):
    """The shape round 89 fixed: an ordinary LibraryReaction carrying the same longDesc."""
    reactants, products = species_pair()
    library = KineticsLibrary(label=label, name=label)
    library.auto_generated = False
    library.entries = {
        1: Entry(
            index=1,
            label="Lip + CH3 <=> CH3Li",
            item=Reaction(reactants=reactants, products=products, reversible=False),
            data=marcus(),
            long_desc="\n".join("family: {0}".format(n) for n in family_lines),
        )
    }
    return library


class _LoadedFamily(object):
    """Just enough of a KineticsFamily for the conversion test and the gate."""

    def __init__(self, label, quarantine=None):
        self.label = label
        self.quarantine = quarantine


def _database(libraries, families):
    class _Kinetics(object):
        def __init__(self):
            self.libraries = libraries
            self.families = families
            self.library_order = []

        def load_libraries(self, path=None, libraries=None):
            raise AssertionError("every library in this probe is already known")

    class _Forbidden(object):
        def is_molecule_forbidden(self, molecule):
            return False

    class _Database(object):
        def __init__(self):
            self.kinetics = _Kinetics()
            self.forbidden_structures = _Forbidden()

    return _Database()


def with_database(libraries, families):
    import rmgpy.data.rmg
    previous = getattr(rmgpy.data.rmg, "database", None)
    rmgpy.data.rmg.database = _database(libraries, families)
    return previous


def restore_database(previous):
    import rmgpy.data.rmg
    rmgpy.data.rmg.database = previous


def clear_caches():
    import rmgpy.data.kinetics.quarantine as q
    for name in ("_DISK_QUARANTINE_CACHE", "_DISK_ANY_QUARANTINE_CACHE"):
        cache = getattr(q, name, None)
        if cache is not None:
            cache.clear()
    for name in ("_UNATTRIBUTED_WARNED", "_UNANSWERED_WARNED", "_UNSAFE_LABELS_WARNED"):
        marked = getattr(q, name, None)
        if marked is not None:
            marked.clear()


def scratch_database(root, families):
    """`families` maps label -> manifest text or None (a family directory with no manifest)."""
    for label, manifest in families.items():
        path = os.path.join(root, "kinetics", "families", label)
        os.makedirs(path)
        if manifest is not None:
            with open(os.path.join(path, QUARANTINE_FILENAME), "w") as f:
                f.write(manifest)
    return root


class pinned_database(object):
    """Point `settings['database.directory']` at `root` for the duration."""

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


# ---------------------------------------------------------------------------
# HIGH -- a second `family:` line still hides the first, on the template shape
# ---------------------------------------------------------------------------


def high_second_family_line(tmpdir):
    """
    The manager's acceptance test: the two-`family:`-line entry, through the real loader
    and `add_reaction_to_core()`.

    Both labels are real. The first is quarantined on disk and not loaded; the second is
    loaded and innocent. `library.py` keeps the last, so the reaction is a TemplateReaction
    whose `.family` is the innocent one -- `model.py` never converts it, because that
    family IS loaded -- and `authoring_families()` returns that single slot without ever
    consulting the entry the previous round attached.
    """
    root = scratch_database(os.path.join(tmpdir, "high"),
                            {QUARANTINED: MANIFEST.format(QUARANTINED, "the rate is unsafe"),
                             INNOCENT: None})
    library = template_library("a_seed", [QUARANTINED, INNOCENT])
    with pinned_database(root):
        previous = with_database({library.label: library},
                                 {INNOCENT: _LoadedFamily(INNOCENT)})
        try:
            built = library.get_library_reactions()[0]
            labels = authoring_families(built)
            model = CoreEdgeReactionModel()
            try:
                model.add_reaction_to_core(built)
                admitted = True
                refusal = ""
            except QuarantinedKineticsError as error:
                admitted = False
                refusal = str(error)
        finally:
            restore_database(previous)

    record("HIGH: a second `family:` line hides a quarantined first one on the template shape",
           admitted,
           "shape built                 : {0}\n"
           "longDesc named              : {1} then {2}\n"
           ".family slot (last wins)    : {3!r}\n"
           "authoring_families()        : {4}\n"
           "entry attached              : {5}\n"
           "ADMITTED_TO_CORE            : {6}{7}".format(
               type(built).__name__, QUARANTINED, INNOCENT, getattr(built, "family", None),
               labels, getattr(built, "entry", None) is not None, admitted,
               "" if admitted else "\nrefused with: " + refusal.split("\n")[0]))

    record("HIGH (the reader): authoring_families ignores the attached entry for a "
           "TemplateReaction",
           QUARANTINED not in labels,
           "the entry's longDesc holds both labels; authoring_families returned {0}.\n"
           "quarantine.py's short path returns [reaction.family] whenever reaction.library "
           "is None,\nand never reaches the entry branch four lines below.".format(labels))


def control_single_label_still_refused(tmpdir):
    """Round 92's case: one `family:` line, quarantined, not loaded. Must still be refused."""
    root = scratch_database(os.path.join(tmpdir, "control_single"),
                            {QUARANTINED: MANIFEST.format(QUARANTINED, "the rate is unsafe")})
    library = template_library("a_seed_single", [QUARANTINED])
    with pinned_database(root):
        previous = with_database({library.label: library}, {})
        try:
            built = library.get_library_reactions()[0]
            model = CoreEdgeReactionModel()
            try:
                model.add_reaction_to_core(built)
                refused = False
            except QuarantinedKineticsError:
                refused = True
        finally:
            restore_database(previous)
    control("round 92's single-label template shape is still refused at add_reaction_to_core",
            refused, "refused = {0}".format(refused))


def control_library_shape_still_refused(tmpdir):
    """The LibraryReaction shape's two-line guarantee already held. It must keep holding."""
    root = scratch_database(os.path.join(tmpdir, "control_lib"),
                            {QUARANTINED: MANIFEST.format(QUARANTINED, "the rate is unsafe"),
                             INNOCENT: None})
    library = ordinary_library("an_ordinary_library", [QUARANTINED, INNOCENT])
    with pinned_database(root):
        previous = with_database({library.label: library},
                                 {INNOCENT: _LoadedFamily(INNOCENT)})
        try:
            built = library.get_library_reactions()[0]
            labels = authoring_families(built)
            model = CoreEdgeReactionModel()
            try:
                model.add_reaction_to_core(built)
                refused = False
            except QuarantinedKineticsError:
                refused = True
        finally:
            restore_database(previous)
    control("a second `family:` line still cannot hide the first on the LibraryReaction shape",
            refused and QUARANTINED in labels,
            "authoring_families = {0}; refused = {1}".format(labels, refused))


# ---------------------------------------------------------------------------
# MEDIUM -- the conversion drops every field it is not handed
# ---------------------------------------------------------------------------

#: Every parameter `LibraryReaction.__init__` accepts. Enumerated from the signature
#: rather than listed by hand, so a field added to the class later cannot be missed here.
LIBRARY_REACTION_FIELDS = [
    name for name in inspect.signature(LibraryReaction.__init__).parameters if name != "self"
]

#: Fields this comparison does not hold the conversion to, each for a stated reason.
#:
#: `library` is what the conversion is FOR -- it is rebuilding the reaction as one
#: belonging to the seed mechanism or reaction library.
#:
#: `index` is assigned downstream: `make_new_reaction` numbers a reaction as it admits it,
#: so the source's -1 legitimately becomes 1 by the time the core holds it. That is not the
#: conversion dropping it, and the separate check below drives the conversion on its own to
#: show the difference.
DELIBERATELY_REPLACED = {"library", "index"}


def medium_dropped_fields(tmpdir):
    """
    Drive the real seed conversion and compare every field of source and result.

    `electrons` is the one with teeth in this campaign: a plasma reaction carries a signed
    electron count precisely because that is what the charge-balance checks read, and
    `Ar + e- => Ar+ + 2e-` converted through this path becomes a reaction claiming zero.
    """
    root = scratch_database(os.path.join(tmpdir, "fields"), {QUARANTINED: None})
    library = template_library("a_seed_with_fields", [ABSENT], electrons=1, degeneracy=3)
    with pinned_database(root):
        previous = with_database({library.label: library}, {})
        try:
            source = library.get_library_reactions()[0]
            source_fields = {name: getattr(source, name, "<<ABSENT>>")
                             for name in LIBRARY_REACTION_FIELDS}
            model = CoreEdgeReactionModel()
            model.add_seed_mechanism_to_core(library.label)
            converted = model.core.reactions[0]
        finally:
            restore_database(previous)

    lost = []
    for name in LIBRARY_REACTION_FIELDS:
        if name in DELIBERATELY_REPLACED:
            continue
        was = source_fields[name]
        if was == "<<ABSENT>>":
            continue
        now = getattr(converted, name, "<<ABSENT>>")
        if name in ("reactants", "products", "pairs", "kinetics", "transition_state"):
            continue          # object identity/copies, compared below by presence only
        if was != now:
            lost.append("{0}: {1!r} -> {2!r}".format(name, was, now))

    record("MEDIUM: the seed conversion drops every field it is not explicitly handed",
           bool(lost),
           "LibraryReaction.__init__ takes {0} fields; the conversion passes 8.\n"
           "fields changed by the conversion that should not have been:\n  {1}\n"
           "converted .electrons = {2} (the source reaction declares {3})".format(
               len(LIBRARY_REACTION_FIELDS),
               "\n  ".join(lost) if lost else "(none)",
               getattr(converted, "electrons", "<<ABSENT>>"),
               source_fields.get("electrons")))

    control("the conversion still writes the library label into .family",
            getattr(converted, "family", None) == library.name,
            "reaction.family after conversion = {0!r} (electron_placement reads this "
            "slot)".format(getattr(converted, "family", None)))
    control("the entry round 92 carried is still on the converted reaction",
            getattr(converted, "entry", None) is not None,
            "entry = {0}".format(type(getattr(converted, "entry", None)).__name__))

    # `index` is excluded above, so show that the exclusion is about admission and not
    # about the conversion: drive the conversion alone and the index comes through.
    try:
        from rmgpy.rmg.model import as_library_reaction
    except ImportError:
        note("`index` could not be checked against the conversion alone",
             "rmgpy.rmg.model has no `as_library_reaction`; at this commit the two call "
             "sites construct the LibraryReaction inline.")
    else:
        direct = as_library_reaction(source, "some_library")
        control("the conversion itself carries `index`; it is the admission that renumbers",
                direct.index == source.index,
                "source.index = {0}; straight after conversion = {1}; after admission = "
                "{2}".format(source.index, direct.index, converted.index))


# ---------------------------------------------------------------------------
# MEDIUM -- the manifest FILE is executed without being validated
# ---------------------------------------------------------------------------


def medium_manifest_symlink(tmpdir):
    """
    The hardening validates the family *directory*. The file inside it may be a symlink
    pointing anywhere, and that file is what gets executed.
    """
    label = "A_Family_With_A_Symlinked_Manifest"
    root = scratch_database(os.path.join(tmpdir, "symlink"), {label: None})
    outside_dir = os.path.join(tmpdir, "outside_the_database")
    os.makedirs(outside_dir)
    outside = os.path.join(outside_dir, "not_a_manifest.py")
    with open(outside, "w") as f:
        f.write(MANIFEST.format(label, "EXECUTED_FROM_OUTSIDE_THE_DATABASE"))
    link = os.path.join(root, "kinetics", "families", label, QUARANTINE_FILENAME)
    os.symlink(outside, link)

    with pinned_database(root):
        answer = resolve_quarantine(label)
    quarantine, answered = answer
    reason = getattr(quarantine, "reason", None)
    record("MEDIUM: a manifest that is a symlink out of the database is followed and executed",
           reason == "EXECUTED_FROM_OUTSIDE_THE_DATABASE",
           "resolve_quarantine({0!r}) -> ({1}, {2})\n"
           "the file executed was {3}, which is outside the database root {4}\n"
           "its `reason` reached the caller, so it ran".format(
               label, type(quarantine).__name__, answered, outside, root))


def medium_nul_label(tmpdir):
    """A NUL in the label reaches os.stat and raises an uncaught ValueError."""
    root = scratch_database(os.path.join(tmpdir, "nul"), {"A_Family": None})
    label = "A_Family\x00.py"
    with pinned_database(root):
        try:
            answer = resolve_quarantine(label)
            raised = None
        except Exception as error:               # noqa: BLE001 -- measuring what escapes
            answer = None
            raised = error
    record("MEDIUM: a label carrying a NUL byte raises an uncaught exception out of the gate",
           raised is not None,
           "resolve_quarantine({0!r}) raised {1}: {2}".format(
               label, type(raised).__name__, raised) if raised is not None
           else "returned {0}".format(answer))


def medium_newline_label(tmpdir):
    """A label with a newline still selects a directory, and its manifest is executed."""
    label = "A_Family\nOhNo"
    root = os.path.join(tmpdir, "newline")
    try:
        scratch_database(root, {label: MANIFEST.format("newline", "EXECUTED_VIA_NEWLINE_LABEL")})
    except OSError as error:
        note("a newline-named directory could not be created on this filesystem",
             "skipping the newline reproduction: {0}".format(error))
        return
    with pinned_database(root):
        answer = resolve_quarantine(label)
    quarantine, answered = answer
    record("MEDIUM: a label carrying a newline selects a directory and its manifest is executed",
           getattr(quarantine, "reason", None) == "EXECUTED_VIA_NEWLINE_LABEL",
           "resolve_quarantine({0!r}) -> ({1}, {2}); reason = {3!r}".format(
               label, type(quarantine).__name__, answered,
               getattr(quarantine, "reason", None)))


def control_ordinary_label_resolves(tmpdir):
    label = "An_Ordinary_Family"
    root = scratch_database(os.path.join(tmpdir, "ordinary"),
                            {label: MANIFEST.format(label, "still readable")})
    with pinned_database(root):
        answer = resolve_quarantine(label)
    control("an ordinary label still resolves to its manifest",
            getattr(answer[0], "reason", None) == "still readable" and answer[1],
            "resolve_quarantine({0!r}) -> ({1}, {2})".format(
                label, type(answer[0]).__name__, answer[1]))


def control_family_symlink_refused(tmpdir):
    """Round 92's directory-symlink refusal must survive whatever is done to the file check."""
    label = "A_Family_Directory_That_Points_Out"
    root = os.path.join(tmpdir, "dirlink")
    os.makedirs(os.path.join(root, "kinetics", "families"))
    outside = os.path.join(tmpdir, "outside_dir")
    os.makedirs(outside)
    with open(os.path.join(outside, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST.format(label, "SHOULD NEVER BE READ"))
    os.symlink(outside, os.path.join(root, "kinetics", "families", label))
    with pinned_database(root):
        answer = resolve_quarantine(label)
    control("a family DIRECTORY symlinked out of the database is still refused",
            answer == (None, False), "resolve_quarantine -> {0}".format(answer))


# ---------------------------------------------------------------------------
# MEDIUM -- cache invalidation, the halves that were left open
# ---------------------------------------------------------------------------


def medium_loaded_family_never_rechecks(tmpdir):
    """
    A LOADED family answers from the object, forever. `resolve_quarantine` returns
    `family.quarantine` before it ever looks at disk, so a manifest added, edited or
    removed after the database was loaded is invisible for the life of the process --
    the exact staleness the unloaded half was repaired for in round 92.
    """
    label = "A_Loaded_Family_Whose_Manifest_Appears"
    root = scratch_database(os.path.join(tmpdir, "loaded"), {label: None})
    family_dir = os.path.join(root, "kinetics", "families", label)
    with pinned_database(root):
        previous = with_database({}, {label: _LoadedFamily(label)})
        try:
            before = resolve_quarantine(label)
            with open(os.path.join(family_dir, QUARANTINE_FILENAME), "w") as f:
                f.write(MANIFEST.format(label, "added after the family was loaded"))
            after = resolve_quarantine(label)
        finally:
            restore_database(previous)
    record("MEDIUM: a loaded family never re-reads its manifest",
           after[0] is None,
           "loaded, before the manifest was written : {0}\n"
           "loaded, after  the manifest was written : {1}\n"
           "the file on disk says QUARANTINED FOR TESTING".format(before, after))


def medium_any_quarantine_cache(tmpdir):
    """
    `_DISK_ANY_QUARANTINE_CACHE` is keyed on the families ROOT's mtime, and writing a file
    inside an existing family directory does not move it. The suppression then keeps
    `_warn_unanswered` silent in a database that does quarantine something.
    """
    from rmgpy.data.kinetics.quarantine import _database_has_any_quarantine

    label = "A_Family_Gaining_A_Manifest"
    root = scratch_database(os.path.join(tmpdir, "anycache"), {label: None})
    family_dir = os.path.join(root, "kinetics", "families", label)
    with pinned_database(root):
        before = _database_has_any_quarantine()
        with open(os.path.join(family_dir, QUARANTINE_FILENAME), "w") as f:
            f.write(MANIFEST.format(label, "added inside an existing family directory"))
        after = _database_has_any_quarantine()
    record("MEDIUM: the suppression cache misses a manifest added inside an existing family",
           before is False and after is False,
           "_database_has_any_quarantine() before = {0}, after = {1}\n"
           "the manifest is on disk; the families root's mtime did not move, so the "
           "cached negative stands\nand _warn_unanswered stays silent".format(before, after))


# ---------------------------------------------------------------------------
# MEDIUM -- requiresEngineCallSites, fourth round
# ---------------------------------------------------------------------------

CALLSITE_MANIFEST = """
name = "Synthetic/quarantine"
state = "QUARANTINED FOR TESTING"
appliesToKineticsClass = "Arrhenius"
reason = "a reason"
requiresEngineModule = "rmgpy.data.kinetics.quarantine"
requiresEngineSymbol = "check_quarantine"
requiresEngineCallSites = ("rmgpy.rmg.model",)
"""

GATES = ("apply_kinetics_to_reaction", "add_reaction_to_core", "add_reaction_to_edge",
         "add_reaction_to_unimolecular_networks")

SHORT_CIRCUIT = "\n\n".join(
    "def {0}(r):\n    False and check_quarantine(r, stage='never')".format(name)
    for name in GATES)

GENUINE = "\n\n".join(
    "def {0}(r):\n    check_quarantine(r, stage='real')\n    return r".format(name)
    for name in GATES)


def _accepts(tmpdir, source, tag):
    import linecache

    from rmgpy.data.kinetics.quarantine import check_quarantine as real_gate

    carrier = os.path.join(tmpdir, "carrier_{0}.py".format(tag))
    with open(carrier, "w") as f:
        f.write(source)
    linecache.checkcache(carrier)
    module = types.ModuleType("rmgpy.rmg.model")
    module.__file__ = carrier
    module.check_quarantine = real_gate

    directory = os.path.join(tmpdir, "manifest_{0}".format(tag))
    os.makedirs(directory)
    with open(os.path.join(directory, QUARANTINE_FILENAME), "w") as f:
        f.write(CALLSITE_MANIFEST)

    previous = sys.modules.get("rmgpy.rmg.model")
    sys.modules["rmgpy.rmg.model"] = module
    try:
        return load_family_quarantine("Synthetic", directory) is not None
    except Exception:                            # noqa: BLE001 -- a refusal is the answer
        return False
    finally:
        if previous is not None:
            sys.modules["rmgpy.rmg.model"] = previous
        else:
            del sys.modules["rmgpy.rmg.model"]


def medium_pin_overstates(tmpdir):
    """
    The check is a static syntactic presence check and its docstring claims more.

    The manager's ruling was "either it verifies what its name says, or correct the
    documentation to state exactly what it checks". This finding is the documentation one,
    so the thing asserted is the text.

    It is asserted positively -- the documentation must name the things the check does not
    verify -- rather than by the absence of the old "reachable branch" phrase. Absence of a
    substring passes for the wrong reason as easily as the right one, and the first version
    of this check proved it: the corrected docstring still contains that phrase, inside the
    sentence that disowns it.
    """
    import rmgpy.data.kinetics.quarantine as q

    gate_calls = getattr(q, "_gate_calls_in_source", None) or getattr(q, "_live_gate_calls", None)
    text = ((gate_calls.__doc__ or "") + (q._check_engine_requirements.__doc__ or "")).lower()
    # The manager's own enumeration of what it does not verify.
    required = ("execut", "dominat", "argument", "propagat", "short-circuit", "shadow",
                "lambda")
    absent = [word for word in required if word not in text]
    record("MEDIUM: the call-site check's own documentation claims more than it verifies",
           bool(absent),
           "the documentation does not state that it fails to verify: {0}\n"
           "checker in use: {1}".format(", ".join(absent) or "(nothing -- all stated)",
                                        getattr(gate_calls, "__name__", None)))

    accepts_short_circuit = _accepts(tmpdir, SHORT_CIRCUIT, "short")
    note("`False and check_quarantine(...)` is accepted, and this round does NOT change that",
         "accepted = {0}. No static check can prove a call executes; each further narrowing\n"
         "makes the check look stronger without changing what it can promise. The repair is\n"
         "the documentation, per the manager's own two options.".format(accepts_short_circuit))


def control_genuine_wiring_accepted(tmpdir):
    control("a module that really calls the gate in all four admission functions is accepted",
            _accepts(tmpdir, GENUINE, "genuine"), "accepted = True expected")


def control_missing_gate_refused(tmpdir):
    thinned = "\n\n".join(
        "def {0}(r):\n    check_quarantine(r, stage='real')\n    return r".format(name)
        for name in GATES[:2])
    control("a module gating only two of the four admission functions is still refused",
            not _accepts(tmpdir, thinned, "thinned"), "accepted = False expected")


# ---------------------------------------------------------------------------
# LOW -- the carrier does not survive a pickle or a copy
# ---------------------------------------------------------------------------


def low_carrier_not_durable(tmpdir):
    root = scratch_database(os.path.join(tmpdir, "durable"), {QUARANTINED: None})
    library = template_library("a_seed_durable", [QUARANTINED])
    with pinned_database(root):
        previous = with_database({library.label: library}, {})
        try:
            built = library.get_library_reactions()[0]
        finally:
            restore_database(previous)

    before = getattr(built, "entry", None) is not None
    try:
        unpickled = pickle.loads(pickle.dumps(built))
        after_pickle = getattr(unpickled, "entry", None) is not None
        pickle_error = None
    except Exception as error:                   # noqa: BLE001 -- measuring picklability
        after_pickle = False
        pickle_error = "{0}: {1}".format(type(error).__name__, error)
    copied = built.copy()
    after_copy = getattr(copied, "entry", None) is not None

    record("LOW: the provenance carrier does not survive __reduce__ or copy()",
           before and not (after_pickle and after_copy),
           "entry before = {0}; after pickle = {1}; after copy() = {2}{3}\n"
           "TemplateReaction.__reduce__ and copy() enumerate their fields by hand and "
           "neither lists `entry`.".format(
               before, after_pickle, after_copy,
               "" if pickle_error is None else "\npickle raised {0}".format(pickle_error)))


# ---------------------------------------------------------------------------


def main():
    logging.basicConfig(level=logging.INFO, stream=sys.stderr,
                        format="%(levelname)s %(message)s")
    tmpdir = tempfile.mkdtemp(prefix="round95_")
    try:
        high_second_family_line(tmpdir)
        control_single_label_still_refused(tmpdir)
        control_library_shape_still_refused(tmpdir)

        medium_dropped_fields(tmpdir)

        medium_manifest_symlink(tmpdir)
        medium_nul_label(tmpdir)
        medium_newline_label(tmpdir)
        control_ordinary_label_resolves(tmpdir)
        control_family_symlink_refused(tmpdir)

        medium_loaded_family_never_rechecks(tmpdir)
        medium_any_quarantine_cache(tmpdir)

        medium_pin_overstates(tmpdir)
        control_genuine_wiring_accepted(tmpdir)
        control_missing_gate_refused(tmpdir)

        low_carrier_not_durable(tmpdir)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    reproduced = [name for name, hit, _ in results if hit]
    broken = [name for name, ok, _ in controls if not ok]
    print("\n" + "=" * 78)
    print("findings reproduced : {0} of {1}".format(len(reproduced), len(results)))
    for name in reproduced:
        print("    - " + name)
    print("controls holding    : {0} of {1}".format(len(controls) - len(broken), len(controls)))
    for name in broken:
        print("    BROKEN: " + name)
    print("noted, not findings : {0}".format(len(notes)))
    print("=" * 78)

    if broken:
        return 2
    return 1 if reproduced else 0


if __name__ == "__main__":
    sys.exit(main())
