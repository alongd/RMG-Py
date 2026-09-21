"""Round-74 RED arms: break one repair at a time and confirm the guard fails.

Two traps are designed around here, both of which have already cost this campaign a
round:

* **Never restore with ``git checkout -- <file>``.** It restores from the INDEX, and
  these repairs are unstaged, so it silently reverts them to HEAD. Three arms then
  measure the same unrepaired engine while appearing to measure three different
  breakages. Restoration here is from a byte backup taken before the first mutation.

* **Never restore with ``shutil.copy2``.** It preserves the backup's mtime, so Cython
  decides the source is older than the ``.so`` and skips the rebuild, leaving the
  BROKEN extension loaded after the restore. ``shutil.copy`` plus an explicit
  ``os.utime`` is what makes the rebuild happen.

Each arm also prints a discriminator read off the BUILT artifact, one value per arm, so
a build that did not take shows up as the wrong line moving rather than as a
mysteriously identical failure set.
"""

import os
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PYX = os.path.join(ROOT, "rmgpy", "chemkin.pyx")
YAML = os.path.join(ROOT, "rmgpy", "yaml_cantera2.py")
LOGS = os.path.join(ROOT, "docs", "i244-chemkin-duplicate-electron-aware", "logs")
PY = "/home/alon/anaconda3/envs/rmg_env/bin/python"

TESTS = os.path.join(ROOT, "test", "rmgpy", "i244ChemkinDuplicateElectronTest.py")

CROSS_CLASS = "test_a_cross_class_pair_is_not_cleared_into_a_deck_cantera_rejects"
PRODUCTION = "test_a_core_only_cantera_export_is_loadable_after_production_marking"
NONMUTATION = "test_a_render_leaves_every_reactions_flag_exactly_as_it_found_it"
THIRD_BODY = "test_a_three_body_is_not_grouped_with_a_falloff_that_writes_another_equation"

ARMS = [
    {
        "name": "D1-class-back-in-key",
        "file": PYX,
        "old": "    return (id(reaction.specific_collider),",
        "new": "    return (reaction.__class__, id(reaction.specific_collider),",
        "rebuild": True,
        "guards": [CROSS_CLASS],
    },
    {
        "name": "D2-cantera-reads-the-flag",
        "file": YAML,
        "old": "    duplicate_flags = chemkin_duplicate_flags(rxns)",
        "new": "    duplicate_flags = [rxn.duplicate for rxn in rxns]",
        "rebuild": False,
        "guards": [PRODUCTION],
    },
    {
        "name": "D3-third-body-bool-only",
        "file": PYX,
        "old": "            pressure_dependent, chemkin_third_body_token_shape(kinetics),",
        "new": "            pressure_dependent,",
        "rebuild": True,
        "guards": [THIRD_BODY],
    },
    {
        "name": "D4-noop-authority",
        "file": PYX,
        "old": "    flags = [False] * len(reactions)",
        "new": "    return [reaction.duplicate for reaction in reactions]\n    flags = [False] * len(reactions)",
        "rebuild": True,
        "guards": [CROSS_CLASS, NONMUTATION],
    },
]


def run(cmd, log_name):
    """Run `cmd`, tee both streams into the ticket's logs directory, return the output."""
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


def discriminators():
    """Four values read off the BUILT extension and the writer module, one per repair."""
    code = (
        "import json;"
        "from rmgpy.chemkin import chemkin_duplicate_group_key, chemkin_duplicate_flags;"
        "from rmgpy.kinetics import Arrhenius, ThirdBody, Troe;"
        "from rmgpy.data.kinetics.library import LibraryReaction;"
        "from rmgpy.data.kinetics.family import TemplateReaction;"
        "from rmgpy.species import Species;"
        "from rmgpy.molecule import Molecule;"
        "import inspect, rmgpy.yaml_cantera2 as y;"
        "a=Species(index=1,label='A',molecule=[Molecule(smiles='[H]')]);"
        "b=Species(index=2,label='B',molecule=[Molecule(smiles='[OH]')]);"
        "lib=LibraryReaction(reactants=[a],products=[b],library='L',reversible=False);"
        "tpl=TemplateReaction(reactants=[a],products=[b],family='F',reversible=False);"
        "marked=LibraryReaction(reactants=[a],products=[b],library='M',reversible=False,duplicate=True);"
        "_lo=Arrhenius(A=(1e6,'m^6/(mol^2*s)'),n=0.0,Ea=(10.0,'kJ/mol'));"
        "_hi=Arrhenius(A=(1e6,'m^3/(mol*s)'),n=0.0,Ea=(10.0,'kJ/mol'));"
        "tb=LibraryReaction(reactants=[a],products=[b],library='T',reversible=False,kinetics=ThirdBody(arrheniusLow=_lo));"
        "fo=LibraryReaction(reactants=[a],products=[b],library='F',reversible=False,kinetics=Troe(arrheniusHigh=_hi,arrheniusLow=_lo,alpha=0.5,T3=(100.0,'K'),T1=(200.0,'K'),T2=(300.0,'K')));"
        "print(json.dumps({"
        "'key_ignores_class': chemkin_duplicate_group_key(lib)==chemkin_duplicate_group_key(tpl),"
        "'flags_recompute_clears_a_lone_mark': chemkin_duplicate_flags([marked])==[False],"
        "'cantera_writer_keys_its_list': 'chemkin_duplicate_flags(rxns)' in inspect.getsource(y._collect_reaction_entries),"
        "'cantera_entry_takes_a_duplicate_arg': 'duplicate' in inspect.signature(y.reaction_to_dict_list).parameters,"
        "'key_separates_third_body_from_falloff': chemkin_duplicate_group_key(tb)!=chemkin_duplicate_group_key(fo),"
        "}))"
    )
    proc = subprocess.run([PY, "-c", code], cwd=ROOT, stdout=subprocess.PIPE,
                          stderr=subprocess.DEVNULL, text=True)
    return proc.stdout.strip().splitlines()[-1] if proc.stdout.strip() else "<probe failed>"


def main():
    backups = {}
    for path in (PYX, YAML):
        backups[path] = path + ".round74-backup"
        shutil.copy(path, backups[path])
    print("byte backups taken:", ", ".join(os.path.basename(v) for v in backups.values()))

    try:
        for arm in ARMS:
            path = arm["file"]
            source = open(path).read()
            assert arm["old"] in source, "arm {0}: anchor not found".format(arm["name"])
            open(path, "w").write(source.replace(arm["old"], arm["new"], 1))
            os.utime(path, None)
            if arm["rebuild"]:
                build("R74-" + arm["name"] + "-build")

            print("\n=== {0} ===".format(arm["name"]))
            print("discriminators:", discriminators())
            proc = run([PY, "-m", "pytest", TESTS, "--no-cov", "-q"],
                       "R74-RED-" + arm["name"])
            tail = [l for l in proc.stdout.splitlines() if l.startswith("FAILED")]
            summary = proc.stdout.strip().splitlines()[-1]
            print("summary:", summary)
            for line in tail:
                print("  ", line.split("::")[-1])
            for guard in arm["guards"]:
                assert any(guard in l for l in tail), \
                    "arm {0} did NOT fail its guard {1}".format(arm["name"], guard)

            # restore THIS arm before the next one
            shutil.copy(backups[path], path)
            os.utime(path, None)
            if arm["rebuild"]:
                build("R74-" + arm["name"] + "-restore-build")
    finally:
        for path, backup in backups.items():
            shutil.copy(backup, path)
            os.utime(path, None)
        build("R74-final-restore-build")
        for backup in backups.values():
            os.remove(backup)

    print("\n=== restored ===")
    print("discriminators:", discriminators())
    proc = run([PY, "-m", "pytest", TESTS, "--no-cov", "-q"], "R74-GREEN-restored")
    print("summary:", proc.stdout.strip().splitlines()[-1])


if __name__ == "__main__":
    main()
