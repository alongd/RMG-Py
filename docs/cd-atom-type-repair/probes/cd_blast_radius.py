#!/usr/bin/env python3
"""Size the Cd/CO/CS/Cdd repair across everything the two checks are given.

``databaseTest.py`` carries the same dead predicate twice::

    num_of_d_bonds = sum([1 if x.order[0] == "D" and len(x.order) == 1 else 0 for x in atom.bonds.values()])

* ``kinetics_check_cd_atom_type`` (:1164 pre-repair) -- the 140 kinetics families;
* ``general_check_cd_atom_type`` (:2024 pre-repair) -- called from ``test_thermo``,
  ``test_solvation``, ``test_statmech`` and ``test_transport``, so it covers those four group
  databases.

A ``Group`` bond stores ``order`` as a list of *numbers*, so ``order[0] == "D"`` is never true,
``num_of_d_bonds`` is always ``0``, and neither check could fail for anything.

This script does not import the repaired file.  It re-implements the shared predicate twice --
verbatim and corrected -- and runs both over every group either check is given, so the size of the
repair is a number rather than a guess.  It reports *reach* alongside the complaint counts, because
a zero from a check that has never run is otherwise indistinguishable from a check that still does
not run, and it refuses to report anything until a negative control shows the corrected predicate
does fire.

Usage:  python cd_blast_radius.py
"""

import sys
from collections import defaultdict

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.molecule import Group
from rmgpy.molecule.atomtype import ATOMTYPES

TARGET_LABELS = ["Cd", "CO", "CS", "Cdd"]
TARGET_ATOM_TYPES = [ATOMTYPES[x] for x in TARGET_LABELS]

REACH = defaultdict(int)


def count_double_bonds(atom, corrected):
    """The checks' ``num_of_d_bonds``, verbatim or with the string comparison repaired."""
    if corrected:
        return sum(
            1 for bond in atom.bonds.values()
            if len(bond.order) == 1 and abs(2 - bond.order[0]) < 1e-7
        )
    return sum(
        1 if bond.order[0] == "D" and len(bond.order) == 1 else 0
        for bond in atom.bonds.values()
    )


def scan_entries(entries, ignore, corrected, count_reach):
    """The body both checks share: return the list of (entry_name, missing_atomtype) complaints."""
    complaints = []
    for entry_name, entry in entries.items():
        if entry in ignore:
            continue
        if not isinstance(entry.item, Group):
            continue
        if count_reach:
            REACH["entries_scanned"] += 1
        for atom in entry.item.atoms:
            if not any(at in TARGET_ATOM_TYPES for at in atom.atomtype):
                continue
            if count_reach:
                REACH["target_typed_atoms"] += 1
            num_of_d_bonds = count_double_bonds(atom, corrected)
            if count_reach and num_of_d_bonds:
                REACH["atoms_with_a_double_bond"] += 1
            correct_atom_list = []
            if num_of_d_bonds == 2:
                correct_atom_list.append("Cdd")
            elif num_of_d_bonds == 1:
                for ligand, bond in atom.bonds.items():
                    if any(abs(2 - order) < 1e-7 for order in bond.order):
                        if ligand.atomtype[0].is_specific_case_of(ATOMTYPES["O"]):
                            correct_atom_list.append("CO")
                        elif ligand.atomtype[0].is_specific_case_of(ATOMTYPES["S"]):
                            correct_atom_list.append("CS")
            for correct_atom in set(correct_atom_list):
                if ATOMTYPES[correct_atom] not in atom.atomtype:
                    complaints.append((entry_name, correct_atom))
    return complaints


def kinetics_ignore(family):
    """`kinetics_check_cd_atom_type`'s ignore list. `general_check_cd_atom_type` has none."""
    if family.own_reverse:
        return []
    ignore = []
    for product in family.forward_template.products:
        ignore.append(product)
        ignore.extend(product.children)
    return ignore


def negative_control():
    """Prove the corrected predicate can fire at all before believing a zero from it."""

    class _Entry:
        def __init__(self, item):
            self.item = item
            self.children = []

    entries = {"BAD": _Entry(Group().from_adjacency_list("1 Cd u0 {2,D}\n2 O2d u0 {1,D}\n"))}
    shipped = scan_entries(entries, [], corrected=False, count_reach=False)
    repaired = scan_entries(entries, [], corrected=True, count_reach=False)
    print("--- negative control: Cd double-bonded to O, typed Cd but not CO ---")
    print(f"  shipped predicate complains:   {shipped}")
    print(f"  repaired predicate complains:  {repaired}")
    ok = not shipped and repaired
    print(f"  control {'PASSES' if ok else 'FAILS'} "
          "(shipped must be silent, repaired must complain)\n")
    return ok


def main():
    if not negative_control():
        print("negative control failed -- the numbers below mean nothing", file=sys.stderr)
        return 1

    database = RMGDatabase()
    database.load(settings["database.directory"], kinetics_families="all")

    totals = defaultdict(int)
    failing = defaultdict(list)

    # kinetics_check_cd_atom_type: the 140 families.
    n_families = len(database.kinetics.families)
    for name in sorted(database.kinetics.families):
        family = database.kinetics.families[name]
        ignore = kinetics_ignore(family)
        for corrected in (False, True):
            got = scan_entries(family.groups.entries, ignore, corrected, count_reach=corrected)
            totals["corrected" if corrected else "verbatim"] += len(got)
            if got and corrected:
                failing[f"kinetics/{name}"] = got

    # general_check_cd_atom_type: the four group databases its callers walk.
    general_sources = [
        ("thermo", getattr(database.thermo, "groups", {}) or {}),
        ("solvation", getattr(database.solvation, "groups", {}) or {}),
        ("statmech", getattr(database.statmech, "groups", {}) or {}),
        ("transport", getattr(database.transport, "groups", {}) or {}),
    ]
    n_groups = 0
    for source_name, groups in general_sources:
        for group_name, group in groups.items():
            n_groups += 1
            for corrected in (False, True):
                got = scan_entries(group.entries, [], corrected, count_reach=corrected)
                totals["corrected" if corrected else "verbatim"] += len(got)
                if got and corrected:
                    failing[f"{source_name}/{group_name}"] = got

    print(f"kinetics families scanned:           {n_families}")
    print(f"thermo/solvation/statmech/transport")
    print(f"  group databases scanned:           {n_groups}")
    print(f"complaints, checks AS SHIPPED:       {totals['verbatim']}")
    print(f"complaints, checks WITH REPAIR:      {totals['corrected']}")
    print(f"group databases that would FAIL:     {len(failing)}")
    for name, complaints in sorted(failing.items()):
        nodes = sorted({c[0] for c in complaints})
        print(f"  {name}: {len(complaints)} complaints across {len(nodes)} nodes")
        for node, wanted in sorted(complaints)[:5]:
            print(f"      node {node!r} is missing atomtype {wanted}")
        if len(complaints) > 5:
            print(f"      ... and {len(complaints) - 5} more")
    print()
    print("reach (a zero above is only meaningful if these are not zero):")
    print(f"  Group entries scanned:             {REACH['entries_scanned']}")
    print(f"  atoms typed Cd/CO/CS/Cdd:          {REACH['target_typed_atoms']}")
    print(f"  ...of those, with a double bond:   {REACH['atoms_with_a_double_bond']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
