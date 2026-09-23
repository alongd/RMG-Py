#!/usr/bin/env python3
"""
Round 89 probe -- the authorship key is right; three paths never deliver it.

    python docs/i221-saturation/probes/round89_delivery_probe.py

Round 87 fixed what the gate reads. Round 89 is about whether anything puts the
value there. Every check below goes through the PRODUCTION object -- a real
``KineticsLibrary``, a real ``CoreEdgeReactionModel`` -- because the defect this
round is built on is a fixture that hand-built the field production omits.

`REPRODUCED` means the defect is present; controls must hold in both directions
or the run measures nothing.
"""

import ast
import inspect
import os
import shutil
import sys
import tempfile
import traceback

from rmgpy.data.base import Entry
from rmgpy.data.kinetics.library import KineticsLibrary, LibraryReaction
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.quarantine import (
    QUARANTINE_FILENAME,
    authoring_family,
    check_quarantine,
    load_family_quarantine,
)
from rmgpy.exceptions import QuarantinedKineticsError
from rmgpy.kinetics import Arrhenius, Marcus
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species

FAMILY = "Fake_Quarantined_Family"

MANIFEST = """
name = "Fake_Quarantined_Family/quarantine"
state = "QUARANTINED FOR TESTING"
appliesToKineticsClass = "Marcus"
reason = "a reason that must reach the error message"
"""

#: The provenance RMG writes into an estimated rate's comment, and saves into an
#: entry's longDesc when the rate is written to a library.
PROVENANCE = ("Estimated using template [Root_2R->C] for rate rule [Root_2R->C]\n"
              "Euclidian distance = 0\n"
              "family: {0}".format(FAMILY))

results = []
controls = []


def record(name, reproduced, detail):
    results.append((name, reproduced, detail))
    print("\n[{0}] {1}\n    {2}".format(
        "REPRODUCED" if reproduced else "NOT REPRODUCED", name,
        detail.replace("\n", "\n    ")))


def control(name, ok, detail):
    controls.append((name, ok, detail))
    print("\n[{0}] control: {1}\n    {2}".format(
        "HOLDS" if ok else "BROKEN", name, detail.replace("\n", "\n    ")))


def marcus(comment=""):
    return Marcus(
        A=(1.73e06, "m^3/(mol*s)"), n=2,
        lmbd_i_coefs=[21824.5, -0.0341626, -0.0013254, 4.92966e-07],
        beta=(1.2e10, "1/m"), wr=(0, "kJ/mol"), wp=(0, "kJ/mol"),
        lmbd_o=(0, "J/mol"), comment=comment,
    )


def reactants_products():
    return ([Species(label="Lip", molecule=[Molecule(smiles="[Li+]")]),
             Species(label="CH3", molecule=[Molecule(smiles="[CH3]")])],
            [Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")])])


def register(quarantine, extra_families=()):
    class _Family(object):
        def __init__(self, q, label):
            self.quarantine = q
            self.label = label

    class _Kinetics(object):
        def __init__(self, families):
            self.families = families

    class _Database(object):
        def __init__(self, families):
            self.kinetics = _Kinetics(families)

    import rmgpy.data.rmg
    families = {FAMILY: _Family(quarantine, FAMILY)}
    for label in extra_families:
        families[label] = _Family(None, label)
    rmgpy.data.rmg.database = _Database(families)


def load_manifest(tmpdir):
    directory = os.path.join(tmpdir, "family")
    if not os.path.isdir(directory):
        os.makedirs(directory)
    with open(os.path.join(directory, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST)
    return load_family_quarantine(FAMILY, directory)


def refused(reaction, kinetics=None):
    try:
        check_quarantine(reaction, stage="the round-89 probe",
                         kinetics=kinetics if kinetics is not None else reaction.kinetics)
    except QuarantinedKineticsError:
        return True
    return False


# --------------------------------------------------------------------------
# HIGH 3 -- the production loader drops the entry
# --------------------------------------------------------------------------
def build_library(tmpdir, label="copied_seed"):
    """A real KineticsLibrary whose one entry carries the provenance in longDesc."""
    reactants, products = reactants_products()
    library = KineticsLibrary(label=label)
    library.entries = {}
    entry = Entry(
        index=1,
        label="Lip + CH3 <=> CH3Li",
        item=Reaction(reactants=reactants, products=products, reversible=True),
        data=marcus(),          # comment deliberately EMPTY: longDesc is the only carrier
        long_desc=PROVENANCE,
    )
    library.entries[1] = entry
    return library, entry


def probe_loader_drops_the_entry(tmpdir):
    register(load_manifest(tmpdir))
    library, entry = build_library(tmpdir)

    reactions = library.get_library_reactions()
    assert len(reactions) == 1
    reaction = reactions[0]

    attached = getattr(reaction, "entry", None) is not None
    recovered = authoring_family(reaction)
    admitted = not refused(reaction)

    control("the entry really did carry the provenance",
            "family:" in (entry.long_desc or ""),
            "entry.long_desc = {0!r}".format((entry.long_desc or "").splitlines()[-1]))

    record("HIGH 3: get_library_reactions() drops the entry, so longDesc is unreachable",
           admitted,
           "entry attached to the LibraryReaction? {0}\n"
           "reaction.family                      = {1!r}\n"
           "authoring_family(reaction)           = {2!r}\n"
           "admitted by the gate                 = {3}".format(
               attached, reaction.family, recovered, admitted))

    # And the half that DOES work today, so the finding is located precisely: when the
    # provenance is in the kinetics comment rather than the longDesc, the gate sees it.
    library2, _ = build_library(tmpdir, label="copied_seed_2")
    library2.entries[1].data.comment = PROVENANCE
    reaction2 = library2.get_library_reactions()[0]
    control("the comment carrier is unaffected -- only the longDesc half is lost",
            refused(reaction2),
            "same library with the provenance in kinetics.comment instead: refused")


# --------------------------------------------------------------------------
# HIGH 2 -- authorship known, family not loaded
# --------------------------------------------------------------------------
def probe_unloaded_family(tmpdir):
    """
    The situation `add_seed_mechanism_to_core` creates on purpose: a reaction whose
    authoring family is NOT among the loaded ones. Authorship is perfectly recoverable;
    the family object is not there to ask.

    Two shapes, and they are not the same question. (a) The family exists in this
    database and merely was not loaded -- answerable, and silence is a defect. (b) The
    family is not in this database at all -- unanswerable, and silence is still a defect
    because the run cannot know it is missing something.
    """
    import logging as _logging

    quarantine = load_manifest(tmpdir)
    register(quarantine)
    reactants, products = reactants_products()
    reaction = LibraryReaction(reactants=reactants, products=products,
                               library="a_seed_from_another_chemistry",
                               kinetics=marcus(PROVENANCE), reversible=True)
    control("with the family loaded, this very reaction is refused",
            refused(reaction),
            "authorship {0!r} resolves while the family is registered".format(
                authoring_family(reaction)))

    # (a) A real database directory holding the family and its manifest -- but a loaded
    # database that does NOT include the family. Only some OTHER family is loaded, and it
    # carries a quarantine, which is the production shape: the run knows about
    # quarantines, just not about this family.
    from rmgpy import settings

    db_root = os.path.join(tmpdir, "db")
    family_dir = os.path.join(db_root, "kinetics", "families", FAMILY)
    os.makedirs(family_dir)
    with open(os.path.join(family_dir, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST)

    other = load_manifest(tmpdir)
    import rmgpy.data.rmg

    class _F(object):
        def __init__(self, q):
            self.quarantine = q
            self.label = "Some_Other_Loaded_Family"

    rmgpy.data.rmg.database.kinetics.families = {"Some_Other_Loaded_Family": _F(other)}

    previous = settings.get("database.directory")
    settings["database.directory"] = db_root
    try:
        on_disk_refused = refused(reaction)
    finally:
        settings["database.directory"] = previous

    record("HIGH 2a: the family is in the database, merely not loaded -- and is not consulted",
           not on_disk_refused,
           "authoring_family {0!r}; loaded families {1}; the manifest is on disk at\n"
           "{2}\nand the rate was {3}".format(
               authoring_family(reaction),
               sorted(rmgpy.data.rmg.database.kinetics.families),
               os.path.join(family_dir, QUARANTINE_FILENAME),
               "ADMITTED without the manifest being read" if not on_disk_refused
               else "refused, the manifest having been read from disk"))

    # (b) The family is nowhere this run can reach. The answer is genuinely unavailable;
    # what is measured is whether the run is TOLD.
    settings["database.directory"] = os.path.join(tmpdir, "an-empty-database")
    records = []

    class _Capture(_logging.Handler):
        def emit(self, record_):
            records.append(record_.getMessage())

    handler = _Capture()
    root = _logging.getLogger()
    root.addHandler(handler)
    previous_level = root.level
    root.setLevel(_logging.WARNING)
    try:
        admitted = not refused(reaction)
    finally:
        root.removeHandler(handler)
        root.setLevel(previous_level)
        settings["database.directory"] = previous

    told = any("cannot know" in message.lower() or "cannot tell" in message.lower()
               for message in records)
    record("HIGH 2b: an unanswerable question is admitted in silence",
           admitted and not told,
           "admitted {0}; the run was told {1}\n"
           "warnings emitted: {2}\n"
           "None from 'no authorship' and None from 'family not found' are the same "
           "value arriving from opposite situations".format(
               admitted, told, records or "none"))


# --------------------------------------------------------------------------
# HIGH 1 -- the pdep network path is gated by nothing
# --------------------------------------------------------------------------
def probe_pdep_bypass(tmpdir):
    register(load_manifest(tmpdir))

    model = CoreEdgeReactionModel()
    reactants, products = reactants_products()
    reaction = TemplateReaction(reactants=reactants, products=products,
                                family=FAMILY, kinetics=marcus(), reversible=True)
    reaction.template = ["Root_2R->C"]

    gated = False
    note = "no exception"
    try:
        model.add_reaction_to_unimolecular_networks(reaction, new_species=reactants[0])
    except QuarantinedKineticsError:
        gated = True
        note = "refused by the gate"
    except Exception as exc:                      # the network machinery, not the gate
        note = "{0}: {1}".format(type(exc).__name__, str(exc).split("\n")[0])

    record("HIGH 1: add_reaction_to_unimolecular_networks consults no gate",
           not gated,
           "core {0}   edge {1}   refused by a gate: {2}\n"
           "outcome: {3}\n"
           "a quarantined rate reaches a pdep network with nothing having looked at it".format(
               len(model.core.reactions), len(model.edge.reactions), gated, note))

    source = inspect.getsource(CoreEdgeReactionModel.add_reaction_to_unimolecular_networks)
    record("HIGH 1 (static): the network method names no gate at all",
           "check_quarantine" not in source,
           "check_quarantine appears in the method's source: {0}\n"
           "the behavioural check above and this one must agree; a behavioural pass with "
           "no call in the source would mean the probe measured the wrong object".format(
               "check_quarantine" in source))


# --------------------------------------------------------------------------
# MEDIUM -- the call-site check pins an import, and authorship is forgeable
# --------------------------------------------------------------------------
def probe_call_site_check(tmpdir):
    """A module that imports the gate and never calls it must not satisfy the check."""
    directory = os.path.join(tmpdir, "callsite")
    os.makedirs(directory)
    sys.path.insert(0, directory)
    with open(os.path.join(directory, "imports_but_never_calls.py"), "w") as f:
        f.write("from rmgpy.data.kinetics.quarantine import check_quarantine\n"
                "# every call site deleted; the import is all that is left\n")

    manifest_dir = os.path.join(tmpdir, "callsite-manifest")
    os.makedirs(manifest_dir)
    with open(os.path.join(manifest_dir, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST
                + 'requiresEngineModule = "rmgpy.data.kinetics.quarantine"\n'
                + 'requiresEngineSymbol = "check_quarantine"\n'
                + 'requiresEngineCallSites = ("imports_but_never_calls",)\n')
    try:
        load_family_quarantine(FAMILY, manifest_dir)
        loaded = True
    except Exception as exc:
        loaded = False
        detail = str(exc).split("\n")[0][:120]
    finally:
        sys.path.remove(directory)

    record("MEDIUM: requiresEngineCallSites accepts a module that never calls the gate",
           loaded,
           "a module whose entire body is the import {0}".format(
               "loaded clean" if loaded else "was refused: " + detail))

    calls = sum(1 for node in ast.walk(ast.parse(inspect.getsource(
        sys.modules["rmgpy.rmg.model"])))
        if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "check_quarantine")
    control("the real call site does contain real calls",
            calls > 0,
            "rmgpy.rmg.model contains {0} call(s) to check_quarantine".format(calls))


def probe_forgeable_authorship(tmpdir):
    register(load_manifest(tmpdir))
    reactants, products = reactants_products()

    shadowed = LibraryReaction(
        reactants=reactants, products=products, library="lib",
        kinetics=marcus("family: An_Innocent_Family\n" + PROVENANCE), reversible=True)
    record("MEDIUM: one prepended line shadows genuine provenance",
           not refused(shadowed),
           "the comment contains BOTH 'family: An_Innocent_Family' and 'family: {0}';\n"
           "authoring_family returns {1!r} because it takes the first line and stops".format(
               FAMILY, authoring_family(shadowed)))

    forged = LibraryReaction(
        reactants=reactants, products=products, library="lib",
        kinetics=Arrhenius(A=(1e13, "cm^3/(mol*s)"), n=0, Ea=(0, "kJ/mol"),
                           comment="family: {0}".format(FAMILY)),
        reversible=True)
    control("a forged line cannot quarantine a rate of an unaffected class",
            not refused(forged),
            "Arrhenius does not match the manifest's criterion, so the criterion still "
            "bounds what a forged label can do")


def main():
    tmpdir = tempfile.mkdtemp(prefix="round89-probe-")
    print("probe scratch: {0}".format(tmpdir))
    for fn in (probe_loader_drops_the_entry, probe_unloaded_family, probe_pdep_bypass,
               probe_call_site_check, probe_forgeable_authorship):
        try:
            fn(tmpdir)
        except Exception:
            traceback.print_exc()
            record(fn.__name__, False, "the probe itself raised -- see traceback above")
    shutil.rmtree(tmpdir, ignore_errors=True)

    print("\n" + "=" * 72 + "\nRESULT\n" + "=" * 72)
    for name, ok, _ in controls:
        print("  {0:<14} control: {1}".format("HOLDS" if ok else "BROKEN", name))
    for name, reproduced, _ in results:
        print("  {0:<14} {1}".format("REPRODUCED" if reproduced else "not reproduced", name))

    broken = [n for n, ok, _ in controls if not ok]
    reproduced = [n for n, r, _ in results if r]
    print("\n{0} of {1} findings reproduced; {2} of {3} controls hold.".format(
        len(reproduced), len(results), len(controls) - len(broken), len(controls)))
    if broken:
        print("CONTROLS BROKEN -- this run measures nothing until they are restored.")
        return 2
    return 1 if reproduced else 0


if __name__ == "__main__":
    sys.exit(main())
