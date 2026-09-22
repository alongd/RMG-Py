#!/usr/bin/env python3
"""
Round 99 probe -- four ways the manifest is addressed by a name instead of by a fact.

    python docs/i221-saturation/probes/round99_manifest_toctou_probe.py

Four findings, one shape between them: something the gate needs to KNOW is taken from a
name, a timestamp, or a field spelling, each of which can be true of two different states
of the world.

1. The directory is validated by name and opened by name, so a symlink swapped into the
   window between the two is followed all the way to ``exec()`` (below).
2. ``_manifest_signature`` maps every ``OSError`` from its ``lstat`` onto ``'absent'``, so
   a manifest that exists and cannot be examined is reported as one that was never there
   -- a clean bill of health issued by the failure of the check itself, then cached.
3. The cache identity ``(mtime_ns, size, inode)`` is entirely restorable: an equal-length
   in-place edit with ``os.utime`` putting the mtime back leaves an edited manifest
   invisible for the rest of the run.
4. ``requiresEngineCallSites`` claims to verify call sites and performs a static syntactic
   presence check; the corrected spelling ``requiresEngineCallsInSource`` is not read at
   all at the base, so a manifest pinning the gate under it pins nothing.

Round 95 closed the symlinked-manifest hole by opening the family *directory* and reading
``quarantine.py`` through that descriptor with ``O_NOFOLLOW``. It left the directory itself
addressed by name twice, independently:

* ``_family_directory()`` resolves the path with ``os.path.realpath`` and checks containment
  with ``startswith`` -- then returns the **unresolved** ``os.path.join(root, label)``;
* ``_read_manifest()`` afterwards calls ``os.open(family_path, O_RDONLY|O_DIRECTORY)`` on
  that string, with no ``O_NOFOLLOW``, so every component including the last is followed.

Check and use are two separate resolutions of the same name. A directory symlink swapped
into the window between them is followed by the open and reaches ``exec()``. That is the
finding: round 95 converted a straightforward escape into a race rather than closing it.

The window is made deterministic here by swapping from inside the check itself, which is
the same technique round 92's vanishing-manifest test uses. The harness asserts that its
own hook fired -- a probe whose instrument falls off the code path reports NOT REPRODUCED
and proves nothing.

``REPRODUCED`` means the defect is present. Controls must hold in both directions or the
run measures nothing. Exit 1 while any finding is present, 2 if a control breaks.
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
    for name in ("_DISK_QUARANTINE_CACHE", "_DISK_ANY_QUARANTINE_CACHE"):
        cache = getattr(q, name, None)
        if cache is not None:
            cache.clear()
    for name in ("_UNATTRIBUTED_WARNED", "_UNANSWERED_WARNED", "_UNSAFE_LABELS_WARNED",
                 "_UNSAFE_MANIFESTS_WARNED", "_LEGACY_CALL_SITES_WARNED"):
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


class swap_during_the_check(object):
    """
    Make the check-then-use window deterministic.

    `_family_directory` calls `os.path.realpath` on the family path; that call IS the
    containment check. We let it compute the true, in-tree answer, then replace the
    directory with a symlink pointing outside the database, then hand the true answer
    back. The check therefore passes on the directory that WAS there, and every later
    user of the returned *string* reaches the one that is there now.

    Nothing about this is exotic to the probe: `os.rename` and `os.symlink` are two
    ordinary syscalls, and a real attacker races them rather than being invited in. The
    hook only decides *when*, not *whether*, the swap is possible.
    """

    def __init__(self, family_path, replacement):
        self.family_path = os.path.abspath(family_path)
        self.replacement = replacement
        self.state = {"fired": False}

    def __enter__(self):
        self.original = os.path.realpath
        state = self.state
        family_path = self.family_path
        replacement = self.replacement
        original = self.original

        def realpath(path, *args, **kwargs):
            answer = original(path, *args, **kwargs)
            if not state["fired"] and os.path.abspath(path) == family_path:
                state["fired"] = True
                os.rename(family_path, family_path + ".was_here")
                os.symlink(replacement, family_path)
            return answer

        os.path.realpath = realpath
        return self.state

    def __exit__(self, *exc):
        os.path.realpath = self.original
        return False


# ---------------------------------------------------------------------------
# The finding -- the family directory is checked by name and opened by name
# ---------------------------------------------------------------------------


def toctou_directory_swap(tmpdir):
    label = "A_Family_Swapped_Under_The_Check"
    root = scratch_database(os.path.join(tmpdir, "race"), {label: None})
    family_path = os.path.join(root, "kinetics", "families", label)

    outside = os.path.join(tmpdir, "outside_the_database")
    os.makedirs(outside)
    with open(os.path.join(outside, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST.format(label, "EXECUTED_FROM_OUTSIDE_THE_DATABASE_VIA_RACE"))

    with pinned_database(root):
        with swap_during_the_check(family_path, outside) as state:
            answer = resolve_quarantine(label)

    assert state["fired"] is True, (
        "this probe did not exercise the race it is named for: the containment check "
        "never called os.path.realpath on the family path, so nothing was swapped")

    quarantine, answered = answer
    reason = getattr(quarantine, "reason", None)
    record("HIGH: a family directory swapped between the containment check and the open is "
           "followed, and the manifest outside the database is executed",
           reason == "EXECUTED_FROM_OUTSIDE_THE_DATABASE_VIA_RACE",
           "resolve_quarantine({0!r}) -> ({1}, {2}); reason = {3!r}\n"
           "the file that would have been executed is {4}\n"
           "the database root it must have come from is {5}\n"
           "the swap fired inside the containment check: {6}".format(
               label, type(quarantine).__name__, answered, reason,
               os.path.join(outside, QUARANTINE_FILENAME), root, state["fired"]))


# ---------------------------------------------------------------------------
# Controls -- each must hold before AND after, or the run measures nothing
# ---------------------------------------------------------------------------


def control_ordinary_label_resolves(tmpdir):
    label = "An_Ordinary_Family"
    root = scratch_database(os.path.join(tmpdir, "ordinary"),
                            {label: MANIFEST.format(label, "still readable")})
    with pinned_database(root):
        answer = resolve_quarantine(label)
    control("an ordinary label still resolves to its manifest and the manifest is read",
            getattr(answer[0], "reason", None) == "still readable" and answer[1],
            "resolve_quarantine({0!r}) -> ({1}, {2})".format(
                label, type(answer[0]).__name__, answer[1]))


def control_family_without_manifest_is_answered_clean(tmpdir):
    label = "A_Family_With_No_Manifest"
    root = scratch_database(os.path.join(tmpdir, "clean"), {label: None})
    with pinned_database(root):
        answer = resolve_quarantine(label)
    control("a family that carries no manifest is still ANSWERED, not unanswered",
            answer == (None, True), "resolve_quarantine({0!r}) -> {1}".format(label, answer))


def control_family_symlink_out_refused(tmpdir):
    """Round 92's static directory-symlink refusal, with no race involved."""
    label = "A_Family_Directory_That_Points_Out"
    root = os.path.join(tmpdir, "dirlink")
    os.makedirs(os.path.join(root, "kinetics", "families"))
    outside = os.path.join(tmpdir, "dirlink_outside")
    os.makedirs(outside)
    with open(os.path.join(outside, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST.format(label, "SHOULD NEVER BE READ"))
    os.symlink(outside, os.path.join(root, "kinetics", "families", label))
    with pinned_database(root):
        answer = resolve_quarantine(label)
    control("a family DIRECTORY symlinked out of the database is refused, unanswered",
            answer == (None, False), "resolve_quarantine({0!r}) -> {1}".format(label, answer))


def control_manifest_symlink_out_refused(tmpdir):
    """Round 95's manifest-file refusal must survive whatever is done to the directory open."""
    label = "A_Family_With_A_Symlinked_Manifest"
    root = scratch_database(os.path.join(tmpdir, "filelink"), {label: None})
    outside_dir = os.path.join(tmpdir, "filelink_outside")
    os.makedirs(outside_dir)
    outside = os.path.join(outside_dir, "not_a_manifest.py")
    with open(outside, "w") as f:
        f.write(MANIFEST.format(label, "SHOULD NEVER BE READ"))
    os.symlink(outside, os.path.join(root, "kinetics", "families", label,
                                     QUARANTINE_FILENAME))
    with pinned_database(root):
        answer = resolve_quarantine(label)
    control("a manifest FILE symlinked out of the database is refused, unanswered",
            answer == (None, False), "resolve_quarantine({0!r}) -> {1}".format(label, answer))


def control_nul_label_refused(tmpdir):
    """Verifier point 2: refused by NAME, and without raising out of the gate."""
    root = scratch_database(os.path.join(tmpdir, "nul"), {"A_Family": None})
    label = "A_Family\x00.py"
    with pinned_database(root):
        try:
            answer, raised = resolve_quarantine(label), None
        except Exception as error:               # noqa: BLE001 -- measuring what escapes
            answer, raised = None, error
    control("a label carrying a NUL byte is refused by name, with nothing raised",
            raised is None and answer == (None, False),
            "resolve_quarantine({0!r}) -> {1}{2}".format(
                label, answer,
                "" if raised is None else "; raised {0}: {1}".format(
                    type(raised).__name__, raised)))


def control_newline_label_refused(tmpdir):
    """Verifier point 2: a newline-named directory really can exist, so the name must be refused."""
    label = "A_Family\nOhNo"
    root = os.path.join(tmpdir, "newline")
    try:
        scratch_database(root, {label: MANIFEST.format("newline", "SHOULD NEVER BE READ")})
    except OSError as error:
        note("a newline-named directory could not be created on this filesystem",
             "the newline control could not be run: {0}".format(error))
        return
    with pinned_database(root):
        answer = resolve_quarantine(label)
    control("a label carrying a newline is refused by name, and its manifest is not read",
            answer == (None, False),
            "resolve_quarantine({0!r}) -> {1}; the directory exists and holds a manifest".format(
                label, answer))


# ---------------------------------------------------------------------------
# Measured, deliberately not a finding
# ---------------------------------------------------------------------------


def note_in_tree_family_symlink(tmpdir):
    """
    A family directory that is a symlink to somewhere else INSIDE the database.

    `_family_directory` admits it -- it resolves under the families root, which is all
    the containment check asks. Refusing the symlink outright is the only form that makes
    containment structural rather than a second resolution of a string, so this case
    changes answer with the repair. Recorded on both sides rather than asserted, because
    it is a deliberate narrowing and not a defect either way.
    """
    label = "A_Family_Linked_To_Its_Sibling"
    real = "A_Real_Family"
    root = scratch_database(os.path.join(tmpdir, "intree"),
                            {real: MANIFEST.format(real, "read through an in-tree symlink")})
    families = os.path.join(root, "kinetics", "families")
    os.symlink(os.path.join(families, real), os.path.join(families, label))
    with pinned_database(root):
        answer = resolve_quarantine(label)
    note("a family directory symlinked to a sibling INSIDE the database",
         "resolve_quarantine({0!r}) -> ({1}, {2}); reason = {3!r}".format(
             label, type(answer[0]).__name__, answer[1],
             getattr(answer[0], "reason", None)))


def note_families_root_is_configuration(tmpdir):
    """
    The families root itself is reached by name, following symlinks, and that is on purpose.

    It is built from ``settings['database.directory']`` -- local configuration -- not from
    a `family:` line in shipped data. Symlinking a whole database checkout is ordinary and
    supported; the label is the component an entry controls, and the repair pins that one.
    """
    label = "A_Family_Under_A_Symlinked_Root"
    real_root = scratch_database(os.path.join(tmpdir, "realroot"),
                                 {label: MANIFEST.format(label, "reached under a linked root")})
    linked_root = os.path.join(tmpdir, "linkedroot")
    os.symlink(real_root, linked_root)
    with pinned_database(linked_root):
        answer = resolve_quarantine(label)
    note("a database root that is itself a symlink still resolves",
         "resolve_quarantine({0!r}) -> ({1}, {2}); reason = {3!r}".format(
             label, type(answer[0]).__name__, answer[1],
             getattr(answer[0], "reason", None)))


# ---------------------------------------------------------------------------
# HIGH -- an unreadable manifest is reported as an absent one
# ---------------------------------------------------------------------------


def unreadable_manifest_reads_as_clean(tmpdir):
    """
    `_manifest_signature` maps every `OSError` from its `lstat` onto `'absent'`. With the
    family directory still visible, `resolve_quarantine` turns that into `(None, True)` --
    "this family carries no manifest" -- and caches it for the life of the run.

    Demonstrated with a permission error rather than a deleted file: a deleted file really
    is absent, and would show nothing about the collapse of "no" into "I could not look".
    """
    label = "A_Family_Whose_Manifest_Cannot_Be_Read"
    root = scratch_database(os.path.join(tmpdir, "unreadable"),
                            {label: MANIFEST.format(label, "must not be read as clean")})
    family = os.path.join(root, "kinetics", "families", label)
    os.chmod(family, 0o000)
    try:
        if os.access(os.path.join(family, QUARANTINE_FILENAME), os.F_OK):
            note("this process reads through a 0o000 directory (running as root?)",
                 "the permission reproduction could not be built here")
            return
        import rmgpy.data.kinetics.quarantine as q
        with pinned_database(root):
            first = resolve_quarantine(label)
            # Read INSIDE the block: `pinned_database.__exit__` clears the caches, and a
            # harness that looks afterwards measures the teardown rather than the code.
            cached = dict(q._DISK_QUARANTINE_CACHE).get((root, label))
            os.chmod(family, 0o755)
            second = resolve_quarantine(label)
    finally:
        os.chmod(family, 0o755)

    record("HIGH: a manifest that exists and cannot be examined is reported as absent, "
           "which is a clean bill of health, and is then cached",
           first == (None, True),
           "with the family directory unreadable : resolve_quarantine -> {0}\n"
           "                      readable again : resolve_quarantine -> ({1}, {2})\n"
           "the first answer says the family carries no manifest; the manifest is there".format(
               first, type(second[0]).__name__, second[1]))

    record("HIGH (second half): the answer taken from a failure is entered into the cache, "
           "under a key that says nothing about the failure that produced it",
           cached is not None and cached[1] == (None, True),
           "_DISK_QUARANTINE_CACHE[({0!r}, {1!r})] = {2}".format(root, label, cached))

    note("the cached clean answer does NOT survive the condition clearing",
         "after the permission was restored: resolve_quarantine -> ({0}, {1}). The review "
         "said the clean answer is 'cached for the rest of the run'; measured, it is not. "
         "Restoring the permission moves the signature from (True, None, 'absent') to "
         "(True, <identity>, 'regular'), which misses the cache and forces a re-read. The "
         "entry is still wrong while the error lasts, and is still removed here on the "
         "principle that an answer produced by a failure is not a property of the key -- "
         "but the permanence in the review does not reproduce.".format(
             type(second[0]).__name__, second[1]))


# ---------------------------------------------------------------------------
# MEDIUM -- the cache identity is restorable
# ---------------------------------------------------------------------------


def equal_length_edit_is_invisible(tmpdir):
    """`(mtime_ns, size, inode)` all survive an in-place edit with the mtime put back."""
    label = "A_Family_Edited_In_Place"
    before = MANIFEST.format(label, "reason A -- xxxxxxxxxxxxxxxxxxxxxxxxxxxx")
    after = MANIFEST.format(label, "reason B -- xxxxxxxxxxxxxxxxxxxxxxxxxxxx")
    assert len(before) == len(after), "the probe's own edit must not change the length"

    root = scratch_database(os.path.join(tmpdir, "inplace"), {label: before})
    manifest = os.path.join(root, "kinetics", "families", label, QUARANTINE_FILENAME)

    with pinned_database(root):
        first = resolve_quarantine(label)
        stamps = os.stat(manifest)
        with open(manifest, "r+") as handle:
            handle.seek(0)
            handle.write(after)
        os.utime(manifest, ns=(stamps.st_atime_ns, stamps.st_mtime_ns))
        repeat = os.stat(manifest)
        assert (repeat.st_mtime_ns, repeat.st_size, repeat.st_ino) == (
            stamps.st_mtime_ns, stamps.st_size, stamps.st_ino), (
            "this probe did not exercise the case it is named for: the edit moved one of "
            "the three fields the old identity was built from")
        second = resolve_quarantine(label)

    record("MEDIUM: an equal-length in-place edit with the mtime restored leaves the "
           "cached answer in place, so an edited manifest never takes effect",
           getattr(second[0], "reason", "").startswith("reason A"),
           "before the edit : reason = {0!r}\n"
           "after the edit  : reason = {1!r}\n"
           "(mtime_ns, size, inode) identical across the edit; st_ctime_ns moved".format(
               getattr(first[0], "reason", None), getattr(second[0], "reason", None)))


# ---------------------------------------------------------------------------
# MEDIUM -- the manifest field claims to verify call sites
# ---------------------------------------------------------------------------


def _loads(tmpdir, tag, extra):
    """Load a manifest declaring `extra`; return (loaded, raised)."""
    directory = os.path.join(tmpdir, "field_" + tag)
    os.makedirs(directory)
    with open(os.path.join(directory, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST.format(tag, "a reason") + extra)
    try:
        return load_family_quarantine("A_Family", directory), None
    except DatabaseError as error:
        return None, error


GATE = ('requiresEngineModule = "rmgpy.data.kinetics.quarantine"\n'
        'requiresEngineSymbol = "check_quarantine"\n')


def the_renamed_field_is_not_read(tmpdir):
    """
    The rename, measured as behaviour rather than as prose.

    `os` binds no gate, so a manifest pinning the gate to it must be refused. Declared
    under the new spelling the field is not read at all at the base, so the manifest
    loads: the pin is inert, which is this campaign's oldest defect -- a pin nothing reads
    is a comment.
    """
    loaded, raised = _loads(tmpdir, "new",
                            GATE + 'requiresEngineCallsInSource = ("os",)\n')
    record("MEDIUM: a manifest pinning the gate under the field's corrected name is not "
           "read, so the pin is inert",
           raised is None and loaded is not None,
           "load_family_quarantine with requiresEngineCallsInSource = (\"os\",) -> "
           "{0}\n`os` does not bind the gate, so this must be refused".format(
               "loaded, unrefused" if raised is None else "refused: {0}".format(raised)))


def control_the_old_spelling_still_pins(tmpdir):
    loaded, raised = _loads(tmpdir, "old", GATE + 'requiresEngineCallSites = ("os",)\n')
    control("the old spelling of the field still refuses a module that binds no gate",
            raised is not None and "wired into" in str(raised),
            "load_family_quarantine with requiresEngineCallSites = (\"os\",) -> {0}".format(
                "loaded, unrefused" if raised is None else "refused"))


def control_a_real_pin_still_loads(tmpdir):
    loaded, raised = _loads(tmpdir, "real",
                            GATE + 'requiresEngineCallSites = ("rmgpy.rmg.model",)\n')
    control("a manifest pinning the real call site still loads",
            raised is None and loaded is not None,
            "load_family_quarantine -> {0}".format(
                "loaded" if raised is None else "refused: {0}".format(raised)))


def note_the_fifth_counterexample(tmpdir):
    """
    Local shadowing satisfies the check. Recorded, not repaired: the ruling was to rename
    the field rather than strengthen the check a sixth time, and each narrowing removes one
    example while leaving the class untouched.
    """
    module_dir = os.path.join(tmpdir, "shadow_site")
    os.makedirs(module_dir)
    with open(os.path.join(module_dir, "shadows_the_gate.py"), "w") as f:
        f.write("from rmgpy.data.kinetics.quarantine import check_quarantine\n"
                "\n"
                "def add_reaction_to_core(rxn):\n"
                "    check_quarantine = lambda *a, **k: None\n"
                "    check_quarantine(rxn)\n")
    sys.path.insert(0, module_dir)
    try:
        loaded, raised = _loads(tmpdir, "shadow",
                                GATE + 'requiresEngineCallSites = ("shadows_the_gate",)\n')
    finally:
        sys.path.remove(module_dir)
    note("the fifth construct that satisfies the call check while gating nothing",
         "a local `check_quarantine = lambda: None` shadowing the import, then called: "
         "{0}. Not repaired -- the ruling was to rename the field, not to narrow the "
         "check again.".format("accepted" if raised is None else "refused"))


def main():
    tmpdir = tempfile.mkdtemp(prefix="round99-toctou-")
    try:
        toctou_directory_swap(tmpdir)
        unreadable_manifest_reads_as_clean(tmpdir)
        equal_length_edit_is_invisible(tmpdir)
        the_renamed_field_is_not_read(tmpdir)
        control_the_old_spelling_still_pins(tmpdir)
        control_a_real_pin_still_loads(tmpdir)
        note_the_fifth_counterexample(tmpdir)

        control_ordinary_label_resolves(tmpdir)
        control_family_without_manifest_is_answered_clean(tmpdir)
        control_family_symlink_out_refused(tmpdir)
        control_manifest_symlink_out_refused(tmpdir)
        control_nul_label_refused(tmpdir)
        control_newline_label_refused(tmpdir)

        note_in_tree_family_symlink(tmpdir)
        note_families_root_is_configuration(tmpdir)
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
