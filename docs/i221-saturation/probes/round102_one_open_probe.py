#!/usr/bin/env python3
"""
Round 102 probe -- one open for the signature and the content, one enumeration for the fields.

    python docs/i221-saturation/probes/round102_one_open_probe.py

Two findings, and both are patterns this campaign has already written down.

**HIGH 1 is round 99's lesson one level up.** Round 99 fixed validate-the-name /
open-the-name for the *path*. The **signature** and the **content** are still two
independent resolutions of one name: `_manifest_signature` stats it, then
`load_family_quarantine` opens it again. Rename between them and B's quarantine is stored
**under A's signature** -- after which restoring A is a cache HIT returning B's answer,
and A's criterion is bypassed without A ever being read. `st_ctime_ns` cannot help; the
cached one is A's.

**HIGH 2 is the fixture-assigns-what-production-never-sets shape.** The template shape in
`get_library_reactions` forwards `degeneracy` and `electrons` from `entry.item` and not
`elementary_high_p`, `allow_pdep_route` or `allow_max_rate_violation`; all three shapes
drop the last. A reaction that loses `elementary_high_p` silently misses pressure-
dependent routing. The round-95 test could not catch it because it assigned those three
onto the loader's *output*, so the converter was only asked to carry what the fixture had
put there.

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
from rmgpy.kinetics.arrhenius import Marcus
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
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


def marcus():
    return Marcus(A=(1.73e06, "m^3/(mol*s)"), n=2,
                  lmbd_i_coefs=[21824.5, -0.0341626, -0.0013254, 4.92966e-07],
                  beta=(1.2e10, "1/m"), wr=(0, "kJ/mol"), wp=(0, "kJ/mol"),
                  lmbd_o=(0, "J/mol"), comment="Estimated from node Root")


def library_declaring(auto=True, **flags):
    """
    A library whose single entry declares `flags` ON `entry.item`.

    That is where `KineticsLibrary.load_entry` puts them when a library file declares
    them, and the only place they can come from. Nothing here touches the loader's
    output -- a probe that set these on the built reaction would be measuring itself.
    """
    long_desc = "\n".join([
        "Matched reaction 3 Lip + CH3 <=> CH3Li in A_Family/rate rule [Root]",
        "Euclidian distance = 0",
        "family: A_Family",
    ])
    library = KineticsLibrary(label="a_seed", name="a_seed")
    library.auto_generated = auto
    library.entries = {
        1: Entry(
            index=1, label="Lip + CH3 <=> CH3Li",
            item=Reaction(
                reactants=[Species(label="Lip", molecule=[Molecule(smiles="[Li+]")],
                                   reactive=False)],
                products=[Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")],
                                  reactive=False)],
                reversible=False, electrons=1, degeneracy=3, **flags),
            data=marcus(),
            long_desc=long_desc if auto else "Originally from reaction library: other",
        )
    }
    return library


FLAGS = ("elementary_high_p", "allow_pdep_route", "allow_max_rate_violation")


# ---------------------------------------------------------------------------
# HIGH 1 -- the signature and the content come from two different opens
# ---------------------------------------------------------------------------


def rename_between_the_stat_and_the_read(tmpdir):
    label = "A_Family_Whose_Manifest_Is_Renamed_Mid_Lookup"
    root = os.path.join(tmpdir, "rename")
    family = os.path.join(root, "kinetics", "families", label)
    os.makedirs(family)
    manifest = os.path.join(family, QUARANTINE_FILENAME)
    spare = os.path.join(family, "manifest_b.py")
    stashed = os.path.join(family, "manifest_a.py")
    with open(manifest, "w") as f:
        f.write(MANIFEST.format(label, "manifest A"))
    with open(spare, "w") as f:
        f.write(MANIFEST.format(label, "manifest B"))

    state = {"fired": False}
    original = os.lstat
    target = os.path.abspath(manifest)

    def lstat(path, *args, **kwargs):
        answer = original(path, *args, **kwargs)
        if not state["fired"] and os.path.abspath(str(path)) == target:
            state["fired"] = True
            os.rename(target, stashed)
            os.rename(spare, target)
        return answer

    with pinned_database(root):
        os.lstat = lstat
        try:
            first = resolve_quarantine(label)
        finally:
            os.lstat = original

        assert state["fired"] is True, (
            "this probe did not exercise the rename it is named for: the signature was "
            "never taken on the manifest path")

        import rmgpy.data.kinetics.quarantine as q
        cached = q._DISK_QUARANTINE_CACHE.get((root, label))
        parsed = os.stat(target)          # `quarantine.py` is B now; B is what was read
        parsed_identity = (parsed.st_mtime_ns, parsed.st_ctime_ns, parsed.st_size,
                           parsed.st_ino)

        # Put A back and ask again, which is the exploitation the brief describes.
        os.rename(target, spare)
        os.rename(stashed, target)
        second = resolve_quarantine(label)

    cached_identity = cached[0][1] if cached else None
    record("HIGH: the cached answer is stored under the identity of a file that was "
           "never parsed -- the signature comes from one open and the content from "
           "another",
           cached_identity != parsed_identity,
           "the read picked up B (that is the race, not the defect):\n"
           "    first -> reason = {0!r}\n"
           "the cache now holds B's quarantine under:\n"
           "    cached identity = {1}\n"
           "    B's identity    = {2}\n"
           "a key and a value that describe two different files".format(
               getattr(first[0], "reason", None), cached_identity, parsed_identity))

    note("the brief's exploitation does NOT reproduce at this base, and the reason is "
         "round 99's own repair",
         "with the ORIGINAL manifest restored, resolve_quarantine -> reason = {0!r}, not "
         "'manifest B'. A rename moves the inode's st_ctime_ns, and round 99 put "
         "st_ctime_ns into the identity, so the restored file no longer matches the "
         "stale key and the cache MISSES. The brief says ctime 'does not help'; measured, "
         "it is what makes this unexploitable by rename on an ordinary POSIX filesystem. "
         "The key still describes a file that was never parsed, which is the finding "
         "above, and is worth closing on its own terms.".format(
             getattr(second[0], "reason", None)))


def control_the_cached_key_describes_the_parsed_file(tmpdir):
    import rmgpy.data.kinetics.quarantine as q

    label = "An_Ordinary_Family"
    root = os.path.join(tmpdir, "identity")
    family = os.path.join(root, "kinetics", "families", label)
    os.makedirs(family)
    with open(os.path.join(family, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST.format(label, "still readable"))
    with pinned_database(root):
        answer = resolve_quarantine(label)
        cached = q._DISK_QUARANTINE_CACHE.get((root, label))
    control("an ordinary manifest still resolves and is still cached",
            getattr(answer[0], "reason", None) == "still readable" and cached is not None,
            "resolve_quarantine({0!r}) -> ({1}, {2}); cached = {3}".format(
                label, type(answer[0]).__name__, answer[1], cached is not None))


# ---------------------------------------------------------------------------
# HIGH 2 -- fields discarded before the converter ever sees them
# ---------------------------------------------------------------------------


def template_shape_drops_three_flags(tmpdir):
    reaction = library_declaring(**{flag: True for flag in FLAGS}).get_library_reactions()[0]
    dropped = [flag for flag in FLAGS if getattr(reaction, flag) is not True]
    record("HIGH: the template shape drops every flag it does not name by hand, so both "
           "conversions in model.py carry false defaults onward",
           bool(dropped),
           "entry.item declares {0}\n"
           "the built {1} carries {2}\n"
           "dropped: {3}\n"
           "`elementary_high_p` lost here means the reaction silently misses "
           "pressure-dependent routing".format(
               {flag: True for flag in FLAGS}, type(reaction).__name__,
               {flag: getattr(reaction, flag) for flag in FLAGS}, dropped or "none"))


def every_shape_drops_allow_max_rate_violation(tmpdir):
    reaction = library_declaring(auto=False, allow_max_rate_violation=True
                                 ).get_library_reactions()[0]
    record("HIGH (second half): the ordinary library shape names two of the three flags "
           "and drops `allow_max_rate_violation`, so no shape carries it",
           reaction.allow_max_rate_violation is not True,
           "entry.item declares allow_max_rate_violation = True\n"
           "the built {0} carries {1}".format(
               type(reaction).__name__, reaction.allow_max_rate_violation))


def control_electrons_and_degeneracy_still_carried(tmpdir):
    reaction = library_declaring().get_library_reactions()[0]
    control("round 95's repair is intact: electrons and degeneracy still arrive",
            reaction.electrons == 1 and reaction.degeneracy == 3,
            "electrons = {0}, degeneracy = {1}".format(
                reaction.electrons, reaction.degeneracy))


def control_a_false_flag_stays_false(tmpdir):
    reaction = library_declaring().get_library_reactions()[0]
    control("carrying is not setting: a flag the entry leaves false stays false",
            not any(getattr(reaction, flag) for flag in FLAGS),
            {flag: getattr(reaction, flag) for flag in FLAGS})


def control_the_entry_still_travels(tmpdir):
    reaction = library_declaring().get_library_reactions()[0]
    control("round 92's carrier is intact: the entry still travels with the reaction",
            getattr(reaction, "entry", None) is not None,
            "entry = {0}".format(type(getattr(reaction, "entry", None)).__name__))


def note_the_out_of_gates_sibling(tmpdir):
    """
    The census's fourth site, in a file this round may not edit.

    `KineticsDatabase.generate_reactions_from_library` builds a `LibraryReaction` from
    `entry.item` naming eight fields and dropping the rest -- including `electrons`,
    which is the charge-balance field round 95 was about.
    """
    import inspect

    from rmgpy.data.kinetics.database import KineticsDatabase

    source = inspect.getsource(KineticsDatabase.generate_reactions_from_library)
    named = [field for field in ("electrons", "allow_pdep_route", "elementary_high_p",
                                 "allow_max_rate_violation", "network_kinetics", "pairs",
                                 "transition_state")
             if field + "=" in source]
    note("rmgpy/data/kinetics/database.py: generate_reactions_from_library builds from "
         "entry.item and drops the same class of field",
         "of (electrons, allow_pdep_route, elementary_high_p, allow_max_rate_violation, "
         "network_kinetics, pairs, transition_state) it forwards: {0}\n"
         "out of this round's gates -- named, not fixed".format(named or "none"))


def main():
    tmpdir = tempfile.mkdtemp(prefix="round102-one-open-")
    try:
        rename_between_the_stat_and_the_read(tmpdir)
        template_shape_drops_three_flags(tmpdir)
        every_shape_drops_allow_max_rate_violation(tmpdir)

        control_the_cached_key_describes_the_parsed_file(tmpdir)
        control_electrons_and_degeneracy_still_carried(tmpdir)
        control_a_false_flag_stays_false(tmpdir)
        control_the_entry_still_travels(tmpdir)
        note_the_out_of_gates_sibling(tmpdir)
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
