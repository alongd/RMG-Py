#!/usr/bin/env python3
"""
Round 101 probe -- the descent protects the last component, not the path.

    python docs/i221-saturation/probes/round101_descent_probe.py

Round 99 replaced a validate-the-name/open-the-name pair with a descent: open the parent,
open the family label through it with ``O_NOFOLLOW``, open the manifest through that. Two
of the three opens are pinned. The first is not, and the first is
``<database>/kinetics/families``.

What makes this worth a round of its own rather than a footnote is that
``_family_directory``'s containment check *cannot* catch it. It compares
``realpath(family_path)`` against ``realpath(families_root)`` -- and when ``families`` is
itself a link, both resolve into the same foreign tree, so the prefix test passes and
reports containment in a directory that is not the database. A check and the thing it
protects, disagreeing about which tree they are talking about.

Five findings, plus a LOW:

1. ``kinetics/families`` replaced by a link reaches a manifest outside the database. No
   race: the link is simply there, and the containment check approves.
2. The same one level up, at ``kinetics`` -- so a repair naming one level is not a repair.
3. The ``dir_fd``-unsupported fallback opens the joined pathname and follows every link in
   it. A documented hole is still a hole.
4. An empty ``family_path`` makes ``os.path.join('', 'quarantine.py')`` a *relative* name,
   so the fallback executes whatever ``quarantine.py`` sits in the working directory.
5. A permission error on ``kinetics/families`` makes ``os.path.isdir`` return False, and
   ``resolve_quarantine`` tested existence before ``kind`` -- so a **loaded** family
   answered from its own attribute with ``answered=True``. Round 99's defect on the branch
   round 99 did not test.
6. (LOW) the both-spellings refusal tests truthiness, not presence.

``REPRODUCED`` means the defect is present. Controls must hold in both directions. Exit 1
while any finding is present, 2 if a control breaks.
"""

import os
import shutil
import sys
import tempfile

from rmgpy import settings
from rmgpy.data.kinetics.quarantine import (
    QUARANTINE_FILENAME,
    load_family_quarantine,
    resolve_quarantine,
)
from rmgpy.exceptions import DatabaseError

MANIFEST = """
name = "{0}/quarantine"
state = "QUARANTINED FOR TESTING"
appliesToKineticsClass = "Marcus"
reason = "{1}"
"""

GATE = ('requiresEngineModule = "rmgpy.data.kinetics.quarantine"\n'
        'requiresEngineSymbol = "check_quarantine"\n')

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
    """A measured fact that is deliberately NOT a finding. Never affects the exit code."""
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


class loaded_families(object):
    """Register `families` (label -> quarantine or None) as the loaded kinetics database."""

    def __init__(self, families):
        self.families = families

    def __enter__(self):
        from rmgpy.data import rmg as rmg_database

        class _Family(object):
            def __init__(self, label, quarantine):
                self.label = label
                self.quarantine = quarantine
                self.quarantine_path = None

        class _Kinetics(object):
            def __init__(self, families):
                self.families = families

        class _Database(object):
            def __init__(self, families):
                self.kinetics = _Kinetics(families)

        self.previous = rmg_database.database
        rmg_database.database = _Database(
            {label: _Family(label, q) for label, q in self.families.items()})
        return self

    def __exit__(self, *exc):
        from rmgpy.data import rmg as rmg_database
        rmg_database.database = self.previous
        return False


def escape_at(tmpdir, tag, component):
    """
    Build a database whose `component` of the descent is a link to a family outside it.

    Returns (database root, label, the manifest that must never be executed).
    """
    label = "A_Family_Reached_Through_A_Linked_Parent"
    outside = os.path.join(tmpdir, "outside_" + tag)
    root = os.path.join(tmpdir, "database_" + tag)

    if component == "families":
        os.makedirs(os.path.join(outside, label))
        target = os.path.join(outside, label, QUARANTINE_FILENAME)
        os.makedirs(os.path.join(root, "kinetics"))
        os.symlink(outside, os.path.join(root, "kinetics", "families"))
    else:
        os.makedirs(os.path.join(outside, "families", label))
        target = os.path.join(outside, "families", label, QUARANTINE_FILENAME)
        os.makedirs(root)
        os.symlink(outside, os.path.join(root, "kinetics"))

    with open(target, "w") as f:
        f.write(MANIFEST.format(label, "EXECUTED_FROM_OUTSIDE_VIA_" + component.upper()))
    return root, label, target


def linked_parent_reaches_outside(tmpdir, component):
    root, label, target = escape_at(tmpdir, component, component)
    with pinned_database(root):
        quarantine, answered = resolve_quarantine(label)
    reason = getattr(quarantine, "reason", None)
    record("HIGH: `{0}` replaced by a link reaches a manifest outside the database"
           .format("kinetics/families" if component == "families" else "kinetics"),
           reason == "EXECUTED_FROM_OUTSIDE_VIA_" + component.upper(),
           "resolve_quarantine({0!r}) -> ({1}, {2}); reason = {3!r}\n"
           "the file executed is {4}\n"
           "the database it must have come from is {5}\n"
           "the containment check approves: with the link in place the root and the "
           "family resolve into the same foreign tree".format(
               label, type(quarantine).__name__, answered, reason, target, root))


def fallback_follows_everything(tmpdir):
    root, label, target = escape_at(tmpdir, "nodirfd", "families")
    original = os.supports_dir_fd
    os.supports_dir_fd = set()
    try:
        with pinned_database(root):
            quarantine, answered = resolve_quarantine(label)
    finally:
        os.supports_dir_fd = original
    record("HIGH: with no directory descriptors the descent falls back to opening the "
           "joined pathname, following every link in it",
           getattr(quarantine, "reason", None) == "EXECUTED_FROM_OUTSIDE_VIA_FAMILIES",
           "os.supports_dir_fd emptied; resolve_quarantine({0!r}) -> ({1}, {2})\n"
           "the file executed is {3}".format(
               label, type(quarantine).__name__, answered, target))


def empty_path_reads_the_working_directory(tmpdir):
    here = os.path.join(tmpdir, "working_directory")
    os.makedirs(here)
    with open(os.path.join(here, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST.format("cwd", "EXECUTED_FROM_THE_WORKING_DIRECTORY"))
    was = os.getcwd()
    os.chdir(here)
    try:
        loaded, raised = load_family_quarantine("A_Family", ""), None
    except Exception as error:                   # noqa: BLE001 -- measuring what escapes
        loaded, raised = None, error
    finally:
        os.chdir(was)
    record("HIGH: an empty family path reads `quarantine.py` out of the process's "
           "working directory and executes it",
           getattr(loaded, "reason", None) == "EXECUTED_FROM_THE_WORKING_DIRECTORY",
           "os.path.split('') is ('', '') and os.path.join('', 'quarantine.py') is a "
           "RELATIVE name\nload_family_quarantine('A_Family', '') -> {0}{1}".format(
               loaded, "" if raised is None else "; raised {0}: {1}".format(
                   type(raised).__name__, raised)))
    if raised is not None:
        note("the empty-path refusal raised instead of refusing",
             "{0}: {1} -- a gate that raises out of the path it is refusing has failed "
             "in the direction this module exists to avoid".format(
                 type(raised).__name__, raised))


def loaded_family_behind_an_unreadable_parent(tmpdir):
    label = "A_Loaded_Family_Behind_An_Unreadable_Parent"
    root = os.path.join(tmpdir, "unreadable_parent")
    families = os.path.join(root, "kinetics", "families")
    os.makedirs(os.path.join(families, label))
    with open(os.path.join(families, label, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST.format(label, "must not read as clean"))
    os.chmod(families, 0o000)
    try:
        if os.path.isdir(os.path.join(families, label)):
            note("this process reads through a 0o000 directory (running as root?)",
                 "the loaded-family reproduction could not be built here")
            return
        with pinned_database(root), loaded_families({label: None}):
            answer = resolve_quarantine(label)
    finally:
        os.chmod(families, 0o755)

    record("HIGH: a permission error on `kinetics/families` gives a LOADED family a clean "
           "bill of health, because existence is consulted before readability",
           answer == (None, True),
           "resolve_quarantine({0!r}) -> {1}\n"
           "os.path.isdir cannot look, so it returns False, so the family reads as 'not "
           "where this database would put it' and the loaded object answers instead\n"
           "the family's own attribute is None, so the answer is 'carries no quarantine'"
           .format(label, answer))


def both_spellings_with_one_empty(tmpdir):
    directory = os.path.join(tmpdir, "both_spellings")
    os.makedirs(directory)
    with open(os.path.join(directory, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST.format("both", "a reason") + GATE
                + 'requiresEngineCallsInSource = ()\n'
                + 'requiresEngineCallSites = ("rmgpy.rmg.model",)\n')
    try:
        loaded, raised = load_family_quarantine("A_Family", directory), None
    except DatabaseError as error:
        loaded, raised = None, error
    record("LOW: the both-spellings refusal tests truthiness, not presence, so declaring "
           "the new field EMPTY beside a populated old one is not seen as declaring both",
           raised is None,
           "requiresEngineCallsInSource = () with requiresEngineCallSites populated -> "
           "{0}".format("loaded, unrefused" if raised is None else "refused"))


# ---------------------------------------------------------------------------
# Controls -- each must hold before AND after
# ---------------------------------------------------------------------------


def control_ordinary_label_resolves(tmpdir):
    label = "An_Ordinary_Family"
    root = os.path.join(tmpdir, "ordinary")
    family = os.path.join(root, "kinetics", "families", label)
    os.makedirs(family)
    with open(os.path.join(family, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST.format(label, "still readable"))
    with pinned_database(root):
        answer = resolve_quarantine(label)
    control("an ordinary label still resolves to its manifest and the manifest is read",
            getattr(answer[0], "reason", None) == "still readable" and answer[1],
            "resolve_quarantine({0!r}) -> ({1}, {2})".format(
                label, type(answer[0]).__name__, answer[1]))


def control_the_anchor_may_be_a_link(tmpdir):
    """
    The anchor has to stay followable, and this is what stops the repair overshooting.

    A path cannot be descended without something to descend from, and that first open
    always resolves what it is given. Putting the anchor on `database.directory` puts the
    one followed open on local configuration rather than on anything shipped data steers.
    """
    label = "A_Family_Under_A_Linked_Root"
    real_root = os.path.join(tmpdir, "real_root")
    family = os.path.join(real_root, "kinetics", "families", label)
    os.makedirs(family)
    with open(os.path.join(family, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST.format(label, "reached under a linked root"))
    linked = os.path.join(tmpdir, "linked_root")
    os.symlink(real_root, linked)
    with pinned_database(linked):
        answer = resolve_quarantine(label)
    control("a database root that is itself a link still resolves its manifests",
            answer[1] and answer[0] is not None,
            "resolve_quarantine({0!r}) -> ({1}, {2})".format(
                label, type(answer[0]).__name__, answer[1]))


def control_a_family_without_a_manifest_is_answered(tmpdir):
    label = "A_Family_With_No_Manifest"
    root = os.path.join(tmpdir, "clean")
    os.makedirs(os.path.join(root, "kinetics", "families", label))
    with pinned_database(root):
        answer = resolve_quarantine(label)
    control("a family that carries no manifest is still ANSWERED, not unanswered",
            answer == (None, True),
            "resolve_quarantine({0!r}) -> {1}".format(label, answer))


def control_the_real_pin_still_loads(tmpdir):
    directory = os.path.join(tmpdir, "real_pin")
    os.makedirs(directory)
    with open(os.path.join(directory, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST.format("real", "a reason") + GATE
                + 'requiresEngineCallSites = ("rmgpy.rmg.model",)\n')
    try:
        loaded, raised = load_family_quarantine("A_Family", directory), None
    except DatabaseError as error:
        loaded, raised = None, error
    control("a manifest pinning the real call site under the legacy spelling still loads",
            raised is None and loaded is not None,
            "load_family_quarantine -> {0}".format(
                "loaded" if raised is None else "refused: {0}".format(raised)))


def control_a_loaded_family_still_rereads_a_readable_disk(tmpdir):
    """Round 95's repair must survive the reordering: the disk still wins over the object."""
    label = "A_Loaded_Family_With_A_Readable_Manifest"
    root = os.path.join(tmpdir, "loaded_readable")
    family = os.path.join(root, "kinetics", "families", label)
    os.makedirs(family)
    with open(os.path.join(family, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST.format(label, "read from disk, not from the object"))
    with pinned_database(root), loaded_families({label: None}):
        answer = resolve_quarantine(label)
    control("a loaded family whose attribute is None still picks up a manifest on disk",
            getattr(answer[0], "reason", None) == "read from disk, not from the object",
            "resolve_quarantine({0!r}) -> ({1}, {2})".format(
                label, type(answer[0]).__name__, answer[1]))


def main():
    tmpdir = tempfile.mkdtemp(prefix="round101-descent-")
    try:
        linked_parent_reaches_outside(tmpdir, "families")
        linked_parent_reaches_outside(tmpdir, "kinetics")
        fallback_follows_everything(tmpdir)
        empty_path_reads_the_working_directory(tmpdir)
        loaded_family_behind_an_unreadable_parent(tmpdir)
        both_spellings_with_one_empty(tmpdir)

        control_ordinary_label_resolves(tmpdir)
        control_the_anchor_may_be_a_link(tmpdir)
        control_a_family_without_a_manifest_is_answered(tmpdir)
        control_the_real_pin_still_loads(tmpdir)
        control_a_loaded_family_still_rereads_a_readable_disk(tmpdir)
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
