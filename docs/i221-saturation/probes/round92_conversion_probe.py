#!/usr/bin/env python3
"""
Round 92 probe -- the third reaction shape, and four ways the gate fails open.

    python docs/i221-saturation/probes/round92_conversion_probe.py

Round 89 delivered the authoring entry on the two ``LibraryReaction`` shapes that
``KineticsLibrary.get_library_reactions`` builds. The third shape is a
``TemplateReaction``, and the two places in ``rmgpy/rmg/model.py`` that convert one
into a ``LibraryReaction`` carry neither the entry nor the parsed family across --
one statement after logging that family to the user.

Everything about the HIGH runs through production objects: a real ``KineticsLibrary``,
a real ``CoreEdgeReactionModel``, ``add_seed_mechanism_to_core()`` and
``add_reaction_library_to_edge()``, and the real ``Cation_R_Recombination`` manifest
in the database on disk.

``REPRODUCED`` means the defect is present. Controls must hold in both directions or
the run measures nothing. Exit 1 while any finding is present, 2 if a control breaks.
"""

import logging
import os
import shutil
import sys
import tempfile
import types

from rmgpy import settings
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.library import KineticsLibrary, LibraryReaction
from rmgpy.data.kinetics.quarantine import (
    QUARANTINE_FILENAME,
    authoring_families,
    check_quarantine,
    load_family_quarantine,
    resolve_quarantine,
)
from rmgpy.exceptions import DatabaseError, QuarantinedKineticsError
from rmgpy.kinetics import Marcus
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species

#: A family that really is quarantined in the database on disk, and that a test run
#: does not load. Both halves matter: the gate must be able to answer from disk, and
#: the conversion must be the thing that loses the answer.
QUARANTINED_FAMILY = "Cation_R_Recombination"

results = []
controls = []


def record(name, reproduced, detail):
    results.append((name, reproduced, detail))
    print("\n[{0}] {1}\n    {2}".format(
        "REPRODUCED" if reproduced else "NOT REPRODUCED", name,
        str(detail).replace("\n", "\n    ")))


def control(name, ok, detail):
    controls.append((name, ok, detail))
    print("\n[{0}] control: {1}\n    {2}".format(
        "HOLDS" if ok else "BROKEN", name, str(detail).replace("\n", "\n    ")))


def marcus(comment=""):
    return Marcus(
        A=(1.73e06, "m^3/(mol*s)"), n=2,
        lmbd_i_coefs=[21824.5, -0.0341626, -0.0013254, 4.92966e-07],
        beta=(1.2e10, "1/m"), wr=(0, "kJ/mol"), wp=(0, "kJ/mol"),
        lmbd_o=(0, "J/mol"), comment=comment,
    )


# ---------------------------------------------------------------------------
# A real library, of the shape that produces a TemplateReaction
# ---------------------------------------------------------------------------

#: What RMG writes into an auto-generated library entry's longDesc for an estimated
#: rate. `rate rule` is what selects shape 2 in get_library_reactions; `family:` is
#: the authorship the gate needs. The kinetics comment deliberately carries no family
#: line, so the longDesc is the only carrier -- the case round 89 closed for the other
#: two shapes.
LONG_DESC = ("Matched reaction 3 Lip + CH3 <=> CH3Li in {0}/rate rule [Root]\n"
             "Euclidian distance = 0\n"
             "family: {0}".format(QUARANTINED_FAMILY))


def species_pair():
    a = Species(label="Lip", molecule=[Molecule(smiles="[Li+]")], reactive=False)
    b = Species(label="CH3", molecule=[Molecule(smiles="[CH3]")], reactive=False)
    c = Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")], reactive=False)
    return [a, b], [c]


def template_library(label="a_seed_written_elsewhere"):
    """A real KineticsLibrary whose single entry produces shape 2."""
    reactants, products = species_pair()
    library = KineticsLibrary(label=label, name=label)
    library.auto_generated = True
    library.entries = {
        1: Entry(
            index=1,
            label="Lip + CH3 <=> CH3Li",
            item=Reaction(reactants=reactants, products=products, reversible=False),
            data=marcus(comment="Matched reaction 3 in {0}".format(QUARANTINED_FAMILY)),
            long_desc=LONG_DESC,
        )
    }
    return library


def library_library(label="an_ordinary_library"):
    """A real KineticsLibrary whose single entry produces shape 3 (LibraryReaction)."""
    reactants, products = species_pair()
    library = KineticsLibrary(label=label, name=label)
    library.auto_generated = False
    library.entries = {
        1: Entry(
            index=1,
            label="Lip + CH3 <=> CH3Li",
            item=Reaction(reactants=reactants, products=products, reversible=False),
            data=marcus(comment=""),
            long_desc="family: {0}".format(QUARANTINED_FAMILY),
        )
    }
    return library


class _Kinetics(object):
    def __init__(self, libraries):
        self.libraries = libraries
        self.families = {}          # the quarantined family is NOT loaded
        self.library_order = []

    def load_libraries(self, path=None, libraries=None):
        raise AssertionError("the probe's libraries are all known; nothing should load")


class _Forbidden(object):
    def is_molecule_forbidden(self, molecule):
        return False


class _Database(object):
    def __init__(self, libraries):
        self.kinetics = _Kinetics(libraries)
        self.forbidden_structures = _Forbidden()


def with_database(libraries):
    import rmgpy.data.rmg
    previous = getattr(rmgpy.data.rmg, "database", None)
    rmgpy.data.rmg.database = _Database(libraries)
    return previous


def restore_database(previous):
    import rmgpy.data.rmg
    rmgpy.data.rmg.database = previous


def clear_caches():
    import rmgpy.data.kinetics.quarantine as q
    for name in ("_DISK_QUARANTINE_CACHE",):
        cache = getattr(q, name, None)
        if cache is not None:
            cache.clear()
    for name in ("_UNATTRIBUTED_WARNED", "_UNANSWERED_WARNED"):
        marked = getattr(q, name, None)
        if marked is not None:
            marked.clear()


# ---------------------------------------------------------------------------
# HIGH -- the conversion sites
# ---------------------------------------------------------------------------

def drive(method_name, library):
    """
    Run one real admission path end to end and report what the gate did.

    Returns (authorship_before, refused_before, authorship_after, entry_after,
             admitted, family_slot_after).
    """
    clear_caches()
    previous = with_database({library.label: library})
    try:
        built = library.get_library_reactions()[0]
        before = authoring_families(built)
        try:
            check_quarantine(built, stage="before conversion", kinetics=built.kinetics)
            refused_before = False
        except QuarantinedKineticsError:
            refused_before = True

        model = CoreEdgeReactionModel()
        admitted, after, entry_after, slot = True, None, None, None
        try:
            getattr(model, method_name)(library.label)
        except QuarantinedKineticsError:
            admitted = False
        finally:
            pool = list(model.core.reactions) + list(model.edge.reactions)
            if pool:
                after = authoring_families(pool[0])
                entry_after = getattr(pool[0], "entry", None)
                slot = getattr(pool[0], "family", None)
        return before, refused_before, after, entry_after, admitted, slot
    finally:
        restore_database(previous)


def high_seed():
    library = template_library()
    before, refused_before, after, entry_after, admitted, slot = drive(
        "add_seed_mechanism_to_core", library)
    detail = (
        "before_conversion_authorship = {0}\n"
        "BEFORE_CONVERSION            = {1}\n"
        "after_conversion_authorship  = {2}\n"
        "after_conversion_entry       = {3}\n"
        "ADMITTED_TO_CORE             = {4}".format(
            before, "REFUSED_FROM_DISK" if refused_before else "ADMITTED",
            after, entry_after, admitted))
    record("HIGH: add_seed_mechanism_to_core converts away the authorship it just logged",
           refused_before and admitted, detail)
    control("the conversion leaves the library label in .family for electron_placement",
            slot in (library.label, None),
            "reaction.family after conversion = {0!r}".format(slot))


def high_library():
    library = template_library(label="a_library_written_elsewhere")
    before, refused_before, after, entry_after, admitted, slot = drive(
        "add_reaction_library_to_edge", library)
    detail = (
        "before_conversion_authorship = {0}\n"
        "BEFORE_CONVERSION            = {1}\n"
        "after_conversion_authorship  = {2}\n"
        "ADMITTED_TO_EDGE             = {3}".format(
            before, "REFUSED_FROM_DISK" if refused_before else "ADMITTED", after, admitted))
    record("HIGH (second site): add_reaction_library_to_edge drops it the same way",
           refused_before and admitted, detail)


def control_untouched_shape():
    library = library_library()
    before, refused_before, after, entry_after, admitted, slot = drive(
        "add_reaction_library_to_edge", library)
    control("the shape round 89 fixed is still refused on the same path",
            refused_before and not admitted,
            "authorship {0}; refused before conversion = {1}; refused on the path = {2}".format(
                before, refused_before, not admitted))


# ---------------------------------------------------------------------------
# MEDIUM 1 -- the disk cache
# ---------------------------------------------------------------------------

MANIFEST = """
name = "{0}/quarantine"
state = "QUARANTINED FOR TESTING"
appliesToKineticsClass = "Marcus"
reason = "a reason that must reach the error message"
"""


def scratch_database(root, families):
    """Build a database directory holding `families` -> manifest text or None."""
    for label, manifest in families.items():
        path = os.path.join(root, "kinetics", "families", label)
        os.makedirs(path)
        if manifest is not None:
            with open(os.path.join(path, QUARANTINE_FILENAME), "w") as f:
                f.write(manifest)
    return root


def medium_stale_cache(tmpdir):
    label = "A_Family_That_Changes"
    root = scratch_database(os.path.join(tmpdir, "stale"), {label: None})
    previous = settings.get("database.directory")
    settings["database.directory"] = root
    clear_caches()
    try:
        first = resolve_quarantine(label)
        manifest_path = os.path.join(root, "kinetics", "families", label, QUARANTINE_FILENAME)
        with open(manifest_path, "w") as f:
            f.write(MANIFEST.format(label))
        second = resolve_quarantine(label)
        record("MEDIUM: the disk cache answers from a database state that no longer exists",
               first[0] is None and second[0] is None,
               "before the manifest was written: {0}\n"
               "after  the manifest was written: {1}\n"
               "the file on disk now says QUARANTINED FOR TESTING".format(first, second))
    finally:
        settings["database.directory"] = previous
        clear_caches()


def medium_toctou(tmpdir):
    """
    The window between resolve_quarantine's existence check and the loader's own.

    Reproduced by deleting the manifest from inside that existence check, which is
    what a concurrent database edit does at an unpredictable moment. The shipped code
    then reads a missing file, gets ``None`` back, and records it as a clean bill of
    health -- and caches that.
    """
    label = "A_Family_Whose_Manifest_Vanishes"
    root = scratch_database(os.path.join(tmpdir, "toctou"), {label: MANIFEST.format(label)})
    manifest_path = os.path.join(root, "kinetics", "families", label, QUARANTINE_FILENAME)
    previous = settings.get("database.directory")
    settings["database.directory"] = root
    clear_caches()

    import rmgpy.data.kinetics.quarantine as q
    real_exists = q.os.path.exists
    state = {"fired": False}

    def racing_exists(path):
        answer = real_exists(path)
        if answer and path == manifest_path and not state["fired"]:
            state["fired"] = True
            os.remove(path)          # the race, made deterministic
        return answer

    try:
        q.os.path.exists = racing_exists
        answer = resolve_quarantine(label)
        cached = dict(getattr(q, "_DISK_QUARANTINE_CACHE", {}))
        record("MEDIUM: a manifest that vanishes mid-read is recorded as a clean bill of health",
               answer == (None, True),
               "resolve_quarantine returned {0} -- answered=True means 'this family has no "
               "manifest'\ncache now holds {1}\nthe manifest was there when the check ran and "
               "gone when the loader looked".format(answer, cached))
    finally:
        q.os.path.exists = real_exists
        settings["database.directory"] = previous
        clear_caches()


# ---------------------------------------------------------------------------
# MEDIUM 2 -- the label is an unsanitised path
# ---------------------------------------------------------------------------

def medium_path_escape(tmpdir):
    outside = os.path.join(tmpdir, "outside")
    os.makedirs(outside)
    with open(os.path.join(outside, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST.format("Planted"))

    root = scratch_database(os.path.join(tmpdir, "escape"), {"An_Ordinary_Family": None})
    previous = settings.get("database.directory")
    settings["database.directory"] = root
    clear_caches()
    try:
        absolute = resolve_quarantine(outside)
        record("MEDIUM: an absolute family label discards the database prefix",
               absolute[0] is not None,
               "resolve_quarantine({0!r}) -> {1}\nthe manifest executed is the one planted "
               "outside the database".format(outside, absolute))

        clear_caches()
        depth = os.path.join("..", "..", "..", "outside")
        relative = resolve_quarantine(depth)
        record("MEDIUM: a '..' family label escapes the families directory",
               relative[0] is not None,
               "resolve_quarantine({0!r}) -> {1}".format(depth, relative))

        clear_caches()
        ordinary = resolve_quarantine("An_Ordinary_Family")
        control("an ordinary label still resolves inside the database",
                ordinary == (None, True),
                "resolve_quarantine('An_Ordinary_Family') -> {0}".format(ordinary))
    finally:
        settings["database.directory"] = previous
        clear_caches()


# ---------------------------------------------------------------------------
# MEDIUM 3 -- requiresEngineCallSites still does not prove wiring
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


def _fake_model(tmpdir, source, tag):
    """
    A module object carrying `source`, under the name the manifest declares.

    Nothing is executed: the check reads the module's source with
    ``inspect.getsource`` and compares its binding of the symbol, so a module with
    ``__file__`` set and the real gate bound is exactly the shape it inspects. That
    lets the probe ask what the check accepts without editing rmgpy/rmg/model.py in
    the worktree.
    """
    import linecache
    from rmgpy.data.kinetics.quarantine import check_quarantine as real_gate
    path = os.path.join(tmpdir, "fake_model_{0}.py".format(tag))
    with open(path, "w") as f:
        f.write(source)
    linecache.checkcache(path)
    module = types.ModuleType("rmgpy.rmg.model")
    module.__file__ = path
    module.check_quarantine = real_gate
    return module


def _try_manifest(tmpdir, tag, module):
    directory = os.path.join(tmpdir, "manifest_" + tag)
    os.makedirs(directory)
    with open(os.path.join(directory, QUARANTINE_FILENAME), "w") as f:
        f.write(CALLSITE_MANIFEST)
    previous = sys.modules.get("rmgpy.rmg.model")
    sys.modules["rmgpy.rmg.model"] = module
    try:
        load_family_quarantine("Synthetic", directory)
        return "loaded"
    except DatabaseError as exc:
        return "refused: {0}".format(str(exc)[:120])
    finally:
        if previous is not None:
            sys.modules["rmgpy.rmg.model"] = previous
        else:
            del sys.modules["rmgpy.rmg.model"]


DEAD_CALL = '''\
from rmgpy.data.kinetics.quarantine import check_quarantine


class CoreEdgeReactionModel(object):
    def add_reaction_to_core(self, rxn):
        if False:
            check_quarantine(rxn, stage="never runs")
        self.core.reactions.append(rxn)
'''

UNRELATED_CALL = '''\
from rmgpy.data.kinetics.quarantine import check_quarantine


class CoreEdgeReactionModel(object):
    def add_reaction_to_core(self, rxn):
        somebody_else.check_quarantine(rxn)
        self.core.reactions.append(rxn)
'''


def medium_callsites(tmpdir):
    dead = _try_manifest(tmpdir, "dead", _fake_model(tmpdir, DEAD_CALL, "dead"))
    record("MEDIUM: a call the interpreter can never reach satisfies requiresEngineCallSites",
           dead == "loaded",
           "a module whose only call is inside `if False:` -> {0}".format(dead))

    unrelated = _try_manifest(tmpdir, "unrelated",
                              _fake_model(tmpdir, UNRELATED_CALL, "unrelated"))
    record("MEDIUM: an unrelated attribute call of the same name satisfies it too",
           unrelated == "loaded",
           "a module calling somebody_else.check_quarantine(...) -> {0}".format(unrelated))

    import inspect as _inspect
    import rmgpy.rmg.model as real_model
    real_source = _inspect.getsource(real_model)
    thinned = real_source
    removed = 0
    for stage in ("admission to the model core",
                  "admission to the model edge",
                  "admission to a pressure-dependent network"):
        for line in real_source.split("\n"):
            if "check_quarantine(" in line and stage in line:
                thinned = thinned.replace(line + "\n", "")
                removed += 1
    one_left = _try_manifest(tmpdir, "thinned", _fake_model(tmpdir, thinned, "thinned"))
    record("MEDIUM: rmgpy.rmg.model with {0} of its 4 gates deleted still satisfies it".format(
               removed),
           one_left == "loaded",
           "gate calls removed: {0}; what remains -> {1}".format(removed, one_left))

    real = _try_manifest(tmpdir, "real", real_model)
    control("the real rmgpy.rmg.model satisfies the check",
            real == "loaded", "rmgpy.rmg.model -> {0}".format(real))


# ---------------------------------------------------------------------------
# MEDIUM 4 -- the suppression looks only at loaded families
# ---------------------------------------------------------------------------

def medium_silent_admission(tmpdir):
    """
    A database that DOES quarantine something, with nothing loaded yet, admits a rate
    whose family it cannot find -- and says nothing, because the suppression test
    iterates loaded families only.
    """
    root = scratch_database(os.path.join(tmpdir, "silent"),
                            {"A_Quarantined_Family": MANIFEST.format("A_Quarantined_Family")})
    previous = settings.get("database.directory")
    settings["database.directory"] = root
    clear_caches()
    previous_db = with_database({})

    records = []

    class _Catch(logging.Handler):
        def emit(self, rec):
            records.append(rec.getMessage())

    handler = _Catch()
    logging.getLogger().addHandler(handler)
    try:
        reactants, products = species_pair()
        rxn = LibraryReaction(
            reactants=reactants, products=products, library="a_foreign_seed",
            kinetics=marcus(comment="family: A_Family_This_Database_Does_Not_Have"))
        admitted = True
        try:
            check_quarantine(rxn, stage="the round-92 probe", kinetics=rxn.kinetics)
        except QuarantinedKineticsError:
            admitted = False
        warned = [m for m in records if "Cannot tell whether" in m]
        record("MEDIUM: a database with a manifest on disk admits an unanswerable rate in silence",
               admitted and not warned,
               "admitted = {0}; warnings emitted = {1}\n"
               "the database DOES carry a quarantine manifest, at\n{2}\n"
               "but nothing is loaded, and the suppression test iterates loaded families "
               "only".format(admitted, warned,
                             os.path.join(root, "kinetics", "families",
                                          "A_Quarantined_Family", QUARANTINE_FILENAME)))
    finally:
        logging.getLogger().removeHandler(handler)
        restore_database(previous_db)
        settings["database.directory"] = previous
        clear_caches()


# ---------------------------------------------------------------------------

def main():
    tmpdir = tempfile.mkdtemp(prefix="round92-probe-")
    print("probe scratch: {0}".format(tmpdir))
    try:
        high_seed()
        high_library()
        control_untouched_shape()
        medium_stale_cache(tmpdir)
        medium_toctou(tmpdir)
        medium_path_escape(tmpdir)
        medium_callsites(tmpdir)
        medium_silent_admission(tmpdir)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    print("\n" + "=" * 72)
    print("RESULT")
    print("=" * 72)
    for name, ok, _ in controls:
        print("  {0:<14} control: {1}".format("HOLDS" if ok else "BROKEN", name))
    for name, reproduced, _ in results:
        print("  {0:<14} {1}".format("REPRODUCED" if reproduced else "not reproduced", name))
    reproduced = sum(1 for _, r, _ in results if r)
    holding = sum(1 for _, ok, _ in controls if ok)
    print("\n{0} of {1} findings reproduced; {2} of {3} controls hold.".format(
        reproduced, len(results), holding, len(controls)))
    if holding != len(controls):
        return 2
    return 1 if reproduced else 0


if __name__ == "__main__":
    sys.exit(main())
