#!/usr/bin/env python
"""I-262 -- reproduce the anonymous-atomtype refusal and show the named replacement.

Run from the worktree root (rmgrc in cwd):
    /home/alon/anaconda3/envs/rmg_env/bin/python docs/i262-anonymous-atomtype/repro_i262.py

Drives three of the seven I-247 survey crash cells, each in a DIFFERENT family building a
DIFFERENT impossible species, and prints the full exception (class name + message) that reaction
generation raises. Against the base plasma module this is a bare AtomTypeError naming only the
impossible atom's shape; against the i262 module it is an UntypeableStructureError naming the
family, the matched template, and the reacting species as adjacency lists.

Then it feeds each reactant adjacency list the message printed back through
Molecule().from_adjacency_list(...) to prove it round-trips.

Exit code is non-zero if any cell FAILED TO CRASH (a silent skip would be the defect this hunts);
the crashes themselves are the expected, correct behaviour and are caught here only so the driver
can print all cells and check every round-trip in one run.
"""
import sys, os, logging, traceback
logging.disable(logging.WARNING)
from rmgpy.molecule import Molecule
from rmgpy import settings

ADJ = {
    "Ar_ground":   "1 Ar u0 p4 c0",
    "Ar_meta_3P2": "1 Ar u2 p3 c0",
    "Ar+":         "1 Ar u1 p3 c+1",
    "Ar2+":        "1 Ar u0 p3 c+1 {2,S}\n2 Ar u1 p3 c0 {1,S}",
    "H":           "1 H u1 p0 c0",
}
def mol(adj):
    return Molecule().from_adjacency_list(adj, saturate_h=False)

# (label, family, [reactant adj keys]) -- three distinct families, distinct species.
CELLS = [
    ("Ar_ground unimolecular",  "Plasma_Radiative_Recombination", ["Ar_ground"]),
    ("Ar2+ + H",                "R_Recombination",                ["Ar2+", "H"]),
    ("Ar2+ + Ar_meta_3P2",      "Birad_R_Recombination",          ["Ar2+", "Ar_meta_3P2"]),
]
FAMILIES = sorted({fam for _, fam, _ in CELLS})


def main():
    from rmgpy.data.kinetics import KineticsDatabase
    kdb = KineticsDatabase()
    kdb.load_families(settings['database.directory'] + '/kinetics/families', families=FAMILIES)
    print("Loaded families:", FAMILIES)
    print("family.py in use:", sys.modules['rmgpy.data.kinetics.family'].__file__)

    any_did_not_crash = False
    round_trip_failures = []
    for label, fam, keys in CELLS:
        reactants = [mol(ADJ[k]) for k in keys]
        print("\n" + "=" * 78)
        print(f"CELL: {label}   (expected family: {fam})")
        print("=" * 78)
        try:
            rxns = kdb.generate_reactions_from_families(reactants, only_families=[fam])
            print(f"  DID NOT CRASH -- generated {len(rxns)} reactions (unexpected; a silent skip "
                  f"would be the defect this ticket hunts)")
            any_did_not_crash = True
        except Exception as e:
            print(f"  EXCEPTION CLASS: {type(e).__module__}.{type(e).__name__}")
            print(f"  __cause__ CLASS: {type(e.__cause__).__name__ if e.__cause__ else None}")
            print("  MESSAGE:")
            for line in str(e).splitlines():
                print("    | " + line)
            # Round-trip: re-parse every reactant adjacency list the exception carries (i262),
            # else the reactants we fed in (base module has no such attribute).
            adjlists = getattr(e, 'reactants', None)
            source = "exception.reactants (i262)"
            if not adjlists:
                adjlists = [r.to_adjacency_list() for r in reactants]
                source = "driver-side reactants (base module has no attribute)"
            print(f"  ROUND-TRIP of adjacency lists from {source}:")
            for i, adj in enumerate(adjlists):
                try:
                    reparsed = Molecule().from_adjacency_list(adj, saturate_h=False)
                    print(f"    reactant[{i}] re-parsed OK -> {reparsed.get_formula()} "
                          f"charge {reparsed.get_net_charge()} mult {reparsed.multiplicity}")
                except Exception as re:
                    print(f"    reactant[{i}] RE-PARSE FAILED: {type(re).__name__}: {re}")
                    round_trip_failures.append((label, i, str(re)))

    print("\n" + "=" * 78)
    print("SUMMARY")
    print(f"  cells that did NOT crash (bad): {any_did_not_crash}")
    print(f"  round-trip failures: {round_trip_failures}")
    # Non-zero exit only if the refusal stopped being a refusal, or a printed adj list is unparseable.
    if any_did_not_crash or round_trip_failures:
        sys.exit(2)
    print("  OK: every cell refused loudly and every printed adjacency list round-tripped.")


if __name__ == "__main__":
    main()
