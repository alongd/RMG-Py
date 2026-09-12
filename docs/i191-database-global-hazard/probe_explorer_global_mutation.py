#!/usr/bin/env python3
"""
Probe: does ``arkane/explorer.py`` mutate the process-global database in a way that
changes what a later caller resolves, and does that mutation accumulate across jobs?

Mechanism under test (verbatim from ``arkane/explorer.py`` lines 111-118)::

    kinetics_database = get_db('kinetics')
    thermo_database = get_db('thermo')

    thermo_database.libraries['thermojobs'] = thermo_library
    thermo_database.library_order.insert(0, 'thermojobs')

    kinetics_database.libraries['kineticsjobs'] = kinetics_library
    kinetics_database.library_order.insert(0, ('kineticsjobs', 'Reaction Library'))

Nothing undoes this.  ``arkane/main.py`` lines 293-296 run that block once per
``ExplorerJob`` in the input file, in a single process.

This probe does NOT construct a second RMGDatabase.  It is deliberately a test of the
*other* half of the process-global hazard: in-place mutation of the one global.

To keep it honest with real data instead of a fabricated library, the stand-in for the
explorer's generated ``thermojobs`` library is a real, already-loaded RMG thermo library
(GRI-Mech3.0) re-registered under the key ``thermojobs`` -- structurally identical to what
the explorer does (a ThermoLibrary object bound to that key and pushed to the front of
``library_order``).

Run from the repository root of the worktree, so that ./rmgrc is found:

    python docs/i191-database-global-hazard/probe_explorer_global_mutation.py
"""

import os
import sys

from rmgpy import settings

print("database.directory = {}".format(settings["database.directory"]))
print("resolved = {}".format(os.path.abspath(settings["database.directory"])))
print("cwd = {}".format(os.getcwd()))
sys.stdout.flush()

from rmgpy.data.rmg import RMGDatabase, get_db  # noqa: E402
from rmgpy.species import Species  # noqa: E402

BASE_LIB = "primaryThermoLibrary"
DONOR_LIB = "GRI-Mech3.0"
SMILES = "C"  # methane


def resolve(label):
    """Resolve a *freshly created* species through the current global thermo database."""
    spc = Species(smiles=SMILES)
    spc.generate_resonance_structures()
    thermo = get_db("thermo").get_thermo_data(spc)
    print(
        "  [{}] H298 = {:.3f} kJ/mol | comment = {}".format(
            label, thermo.get_enthalpy(298.0) / 1000.0, thermo.comment
        )
    )
    return thermo.get_enthalpy(298.0), thermo.comment


def main():
    db = RMGDatabase()
    db.load_thermo(
        os.path.join(settings["database.directory"], "thermo"),
        thermo_libraries=[BASE_LIB, DONOR_LIB],
        depository=False,
        surface=False,
    )
    tdb = get_db("thermo")
    print("\nloaded libraries      : {}".format(sorted(tdb.libraries.keys())))
    print("library_order (t=0)   : {}".format(tdb.library_order))

    print("\n-- baseline, before any explorer job --")
    h_base, c_base = resolve("baseline")

    donor = tdb.libraries[DONOR_LIB]

    # ---- explorer job #1: the two lines from arkane/explorer.py, verbatim in shape ----
    tdb.libraries["thermojobs"] = donor
    tdb.library_order.insert(0, "thermojobs")
    print("\n-- after explorer job #1 --")
    print("library_order         : {}".format(tdb.library_order))
    h_1, c_1 = resolve("after explorer #1")

    # ---- explorer job #2: same input file, second explorer() block ----
    tdb.libraries["thermojobs"] = donor
    tdb.library_order.insert(0, "thermojobs")
    print("\n-- after explorer job #2 --")
    print("library_order         : {}".format(tdb.library_order))
    print(
        "'thermojobs' appears {} time(s) in library_order".format(
            tdb.library_order.count("thermojobs")
        )
    )
    h_2, c_2 = resolve("after explorer #2")

    print("\n=== findings ===")
    changed = (h_base, c_base) != (h_1, c_1)
    print("resolution changed by explorer job #1 : {}".format(changed))
    print("  baseline           : {:.3f} kJ/mol | {}".format(h_base / 1000.0, c_base))
    print("  after explorer #1  : {:.3f} kJ/mol | {}".format(h_1 / 1000.0, c_1))
    print(
        "library_order duplicate accumulation  : {} entries named 'thermojobs' "
        "after 2 explorer jobs".format(tdb.library_order.count("thermojobs"))
    )
    print("mutation ever undone                  : no (no reset exists; grep 'thermojobs')")

    # Non-zero exit if the probe failed to demonstrate anything, so the exit code is
    # informative rather than always-0.
    if not changed:
        print("\nPROBE INCONCLUSIVE: the two libraries agree on this species; "
              "pick a species where they differ.")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
