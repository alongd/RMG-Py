"""Round-75 RED arms: break one repair at a time and confirm a guard fails.

Same two traps as round 74, both of which have cost this campaign a round:

* **Never restore with ``git checkout -- <file>``.** It restores from the INDEX, and
  these repairs may be unstaged, so it silently reverts them to HEAD. Every arm then
  measures the same unrepaired engine while appearing to measure different breakages.
  Restoration here is from a byte backup taken before the first mutation.

* **Never restore with ``shutil.copy2``.** It preserves the backup's mtime, so Cython
  decides the source is older than the ``.so`` and skips the rebuild, leaving the BROKEN
  extension loaded after the restore. ``shutil.copy`` plus an explicit ``os.utime`` is
  what makes the rebuild happen.

A third trap is specific to this round and is why ``must_fail`` is asserted per arm rather
than "some test went red": **an inert mutation looks exactly like a caught one.**
``rmgpy/reaction.py``'s marking of expansion leaves is now overwritten by
``yaml_cantera1.reaction_to_dicts``, so breaking it is invisible from every writer. Its
arm is guarded by a test that calls ``Reaction.to_cantera`` directly. An arm whose guard
does not go red is a repair with no test pointed at it, which is the finding, not an error
in the harness.
"""

import os
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PYX = os.path.join(ROOT, "rmgpy", "chemkin.pyx")
YAML2 = os.path.join(ROOT, "rmgpy", "yaml_cantera2.py")
YAML1 = os.path.join(ROOT, "rmgpy", "yaml_cantera1.py")
REACTION = os.path.join(ROOT, "rmgpy", "reaction.py")
LOGS = os.path.join(ROOT, "docs", "i244-chemkin-duplicate-electron-aware", "logs")
PY = "/home/alon/anaconda3/envs/rmg_env/bin/python"

I244 = os.path.join(ROOT, "test", "rmgpy", "i244ChemkinDuplicateElectronTest.py")
PLASMA = os.path.join(ROOT, "test", "rmgpy", "plasmaExportTest.py")

WRITER1_CORE = "test_writer1_core_only_export_is_loadable_after_production_marking"
PAIRWISE_AGREES = "test_the_pairwise_marking_production_uses_agrees_with_the_group_key"
CHEMKIN_ONE_LEAF = "test_the_chemkin_deck_cantera_converts_also_loads"
WRITERS_ONE_LEAF = "test_both_cantera_writers_load"
TO_CANTERA_ONE_LEAF = "test_to_cantera_alone_does_not_mark_a_one_leaf_wrapper"
WRITER1_COLLIDER = "test_writer1_refuses_a_collider_it_cannot_put_in_the_equation"
WRITER2_COLLIDER = "test_a_grouped_reaction_with_an_unrenderable_collider_is_refused"

ARMS = [
    {
        # HIGH 1: put Writer1 back to serialising the object's own flag.
        "name": "D5-writer1-reads-the-flag",
        "file": YAML1,
        "old": "    duplicate_flags = chemkin_duplicate_flags(rxn_list)",
        "new": "    duplicate_flags = [rxn.duplicate for rxn in rxn_list]",
        "rebuild": False,
        "tests": [I244],
        "must_fail": [WRITER1_CORE],
    },
    {
        # HIGH 2: put the class comparison back into the PAIRWISE marker. The writers
        # still recompute, so only the test that reaches the pairwise path goes red --
        # which is the point: that is the test round 75 said did not exist.
        "name": "D6-class-back-in-pairwise",
        "file": PYX,
        "old": "    for reaction2 in reaction_list:\n        same_dir_match =",
        "new": ("    for reaction2 in reaction_list:\n"
                "        if reaction1.__class__ != reaction2.__class__:\n"
                "            continue\n"
                "        same_dir_match ="),
        "rebuild": True,
        "tests": [I244],
        "must_fail": [PAIRWISE_AGREES],
    },
    {
        # HIGH 3a: Chemkin marks every leaf unconditionally again.
        "name": "D7-chemkin-marks-every-leaf",
        "file": PYX,
        "old": "        leaf_duplicate = len(reaction.kinetics.arrhenius) > 1 or duplicate",
        "new": "        leaf_duplicate = True",
        "rebuild": True,
        "tests": [I244],
        "must_fail": [CHEMKIN_ONE_LEAF],
    },
    {
        # HIGH 3b: Writer2 marks every leaf unconditionally again.
        "name": "D8-writer2-marks-every-leaf",
        "file": YAML2,
        "old": ("        sub_duplicate = len(sub_kinetics_list) > 1 or bool(\n"
                "            reaction.duplicate if duplicate is None else duplicate)"),
        "new": "        sub_duplicate = True",
        "rebuild": False,
        "tests": [I244],
        "must_fail": [WRITERS_ONE_LEAF],
    },
    {
        # HIGH 3c: to_cantera marks every leaf unconditionally again. Invisible through
        # every writer, because Writer1 overwrites the flag -- so this arm is the one
        # that proves its guard reaches the site rather than a writer downstream of it.
        "name": "D9-to-cantera-marks-every-leaf",
        "file": REACTION,
        "old": "            duplicate = self.duplicate or len(ct_reaction) > 1",
        "new": "            duplicate = True",
        "rebuild": True,
        "tests": [I244],
        "must_fail": [TO_CANTERA_ONE_LEAF],
    },
    {
        # MEDIUM a: Writer2 drops an unrenderable collider silently again.
        "name": "D10-writer2-drops-the-collider",
        "file": YAML2,
        "old": "    _collider = getattr(reaction, 'specific_collider', None)\n    if _collider is not None",
        "new": "    _collider = None\n    if _collider is not None",
        "rebuild": False,
        "tests": [PLASMA],
        "must_fail": [WRITER2_COLLIDER],
    },
    {
        # MEDIUM b: Writer1 drops an unrenderable collider silently again.
        "name": "D11-writer1-drops-the-collider",
        "file": YAML1,
        "old": "    if collider is not None and not isinstance(kin, (ThirdBody, Lindemann, Troe)):",
        "new": "    if False:",
        "rebuild": False,
        "tests": [I244],
        "must_fail": [WRITER1_COLLIDER],
    },
]


def run(cmd, log_name):
    with open(os.path.join(LOGS, log_name + "-stdout.log"), "w") as out, \
            open(os.path.join(LOGS, log_name + "-stderr.log"), "w") as err:
        proc = subprocess.run(cmd, cwd=ROOT, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, text=True)
        out.write(proc.stdout)
        err.write(proc.stderr)
    return proc


def build(log_name):
    env = dict(os.environ, PATH="/home/alon/anaconda3/envs/rmg_env/bin:" + os.environ["PATH"])
    with open(os.path.join(LOGS, log_name + "-stdout.log"), "w") as out, \
            open(os.path.join(LOGS, log_name + "-stderr.log"), "w") as err:
        subprocess.run([PY, "setup.py", "build_ext", "--inplace"], cwd=ROOT,
                       stdout=out, stderr=err, env=env)


DISCRIMINATOR_CODE = r"""
import inspect, json
from rmgpy.chemkin import mark_duplicate_reaction, write_kinetics_entry
from rmgpy.kinetics import Arrhenius, MultiArrhenius, PDepArrhenius
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.species import Species
from rmgpy.molecule import Molecule
from rmgpy.thermo import NASA, NASAPolynomial
import rmgpy.yaml_cantera1 as y1
import rmgpy.yaml_cantera2 as y2
from rmgpy.exceptions import MechanismWriterError

def spc(label, index, smiles):
    s = Species(label=label, molecule=[Molecule(smiles=smiles)]); s.index = index
    c = [2.5, 0, 0, 0, 0, -745.375, -11.7246]
    s.thermo = NASA(polynomials=[NASAPolynomial(coeffs=c, Tmin=(200,'K'), Tmax=(1000,'K')),
                                 NASAPolynomial(coeffs=c, Tmin=(1000,'K'), Tmax=(6000,'K'))],
                    Tmin=(200,'K'), Tmax=(6000,'K'))
    return s

# Two reactants, because write_kinetics_entry asserts that the rate's units match the
# molecularity: a unimolecular reaction carrying a bimolecular A aborts inside the writer
# with a units assertion, long before it renders a DUPLICATE line. That is how the first
# version of this probe failed -- silently, printing "<probe failed>" for every arm while
# the arms themselves ran with no build check at all.
A, B, C, AR = (spc('A', 1, '[H]'), spc('B', 2, '[OH]'),
               spc('C', 3, '[H][H]'), spc('Ar', 4, '[Ar]'))
SPCS = [A, B, C, AR]
rate = Arrhenius(A=(1e12, 'cm^3/(mol*s)'), n=0.0, Ea=(0.0, 'kcal/mol'))

lib = LibraryReaction(reactants=[A, B], products=[C], library='L', reversible=False, kinetics=rate)
tpl = TemplateReaction(reactants=[A, B], products=[C], family='F', reversible=False, kinetics=rate)
mark_duplicate_reaction(tpl, [lib])

one = LibraryReaction(reactants=[A, B], products=[C], library='L', reversible=False,
                      kinetics=MultiArrhenius(arrhenius=[rate]))
two = LibraryReaction(reactants=[A, B], products=[C], library='L', reversible=False,
                      kinetics=MultiArrhenius(arrhenius=[rate, rate]))

plog = LibraryReaction(reactants=[A, B], products=[C], library='L', reversible=False,
                       kinetics=PDepArrhenius(pressures=([0.1, 10.0], 'bar'),
                                              arrhenius=[rate, rate]))
plog.specific_collider = AR

multi_with_collider = LibraryReaction(
    reactants=[A, B], products=[C], library='L', reversible=False,
    kinetics=MultiArrhenius(arrhenius=[rate, rate]))
multi_with_collider.specific_collider = AR

def refuses(fn):
    try:
        fn()
    except MechanismWriterError:
        return True
    except Exception:
        return 'wrong-exception'
    return False

# Writer2 refuses an unrenderable collider in TWO places, so removing either one leaves
# the other still refusing -- a plain 'does it raise' discriminator cannot see the
# difference. What changes is WHICH type the message names: refusing only downstream
# reports the Arrhenius of a leaf the caller never constructed.
def refusal_names(fn, wanted):
    try:
        fn()
    except MechanismWriterError as exc:
        return wanted in str(exc)
    except Exception:
        return 'wrong-exception'
    return False

print(json.dumps({
 'writer1_keys_its_own_list':
     'chemkin_duplicate_flags(rxn_list)' in inspect.getsource(y1._collect_reactions),
 'pairwise_marker_ignores_class': [lib.duplicate, tpl.duplicate] == [True, True],
 'chemkin_one_leaf_writes_no_duplicate':
     'DUPLICATE' not in write_kinetics_entry(one, SPCS, verbose=False, duplicate=False),
 'chemkin_two_leaf_still_writes_duplicate':
     write_kinetics_entry(two, SPCS, verbose=False, duplicate=False).count('DUPLICATE') == 2,
 'writer2_one_leaf_not_marked':
     y2.reaction_to_dict_list(one, SPCS, duplicate=False)[0].get('duplicate') is None,
 'to_cantera_one_leaf_not_marked':
     [r.duplicate for r in one.to_cantera(SPCS, use_chemkin_identifier=True)] == [False],
 'writer2_refuses_unrenderable_collider':
     refuses(lambda: y2.reaction_to_dict_list(plog, SPCS)),
 'writer2_refusal_names_the_wrapper':
     refusal_names(lambda: y2.reaction_to_dict_list(multi_with_collider, SPCS),
                   'MultiArrhenius'),
 'writer1_refuses_unrenderable_collider':
     refuses(lambda: y1.reaction_to_dicts(plog, SPCS)),
}))
"""


def discriminators():
    """Eight values read off the BUILT artifacts, so a build that did not take shows up
    as the wrong line moving rather than as a mysteriously identical failure set."""
    proc = subprocess.run([PY, "-c", DISCRIMINATOR_CODE], cwd=ROOT,
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    return proc.stdout.strip().splitlines()[-1] if proc.stdout.strip() else "<probe failed>"


def main():
    files = sorted({arm["file"] for arm in ARMS})
    backups = {path: path + ".round75-backup" for path in files}
    for path, backup in backups.items():
        shutil.copy(path, backup)
    print("byte backups taken:", ", ".join(os.path.basename(v) for v in backups.values()))

    uncaught = []
    try:
        for arm in ARMS:
            path = arm["file"]
            source = open(path).read()
            assert arm["old"] in source, "arm {0}: anchor not found".format(arm["name"])
            open(path, "w").write(source.replace(arm["old"], arm["new"], 1))
            os.utime(path, None)
            if arm["rebuild"]:
                build("R75-" + arm["name"] + "-build")

            print("\n=== {0} ===".format(arm["name"]))
            print("discriminators:", discriminators())
            proc = run([PY, "-m", "pytest"] + arm["tests"] + ["--no-cov", "-q", "-p", "no:randomly"],
                       "R75-RED-" + arm["name"])
            failed = [l for l in proc.stdout.splitlines() if l.startswith("FAILED")]
            print("summary:", proc.stdout.strip().splitlines()[-1])
            for line in failed:
                print("  ", line.split("::")[-1])
            for guard in arm["must_fail"]:
                if not any(guard in l for l in failed):
                    uncaught.append((arm["name"], guard))
                    print("  !! UNCAUGHT: {0} did not go red".format(guard))

            shutil.copy(backups[path], path)
            os.utime(path, None)
            if arm["rebuild"]:
                build("R75-" + arm["name"] + "-restore-build")
    finally:
        for path, backup in backups.items():
            shutil.copy(backup, path)
            os.utime(path, None)
        build("R75-final-restore-build")
        for backup in backups.values():
            os.remove(backup)

    print("\n=== restored ===")
    print("discriminators:", discriminators())
    proc = run([PY, "-m", "pytest", I244, PLASMA, "--no-cov", "-q", "-p", "no:randomly"],
               "R75-GREEN-restored")
    print("summary:", proc.stdout.strip().splitlines()[-1])

    if uncaught:
        print("\nARMS WITH NO GUARD ({0}):".format(len(uncaught)))
        for name, guard in uncaught:
            print("  {0} -> {1}".format(name, guard))
        sys.exit(1)
    print("\nevery arm took its guard with it")


if __name__ == "__main__":
    main()
