#!/usr/bin/env python3
"""Size the Cd/CO/CS/Cdd repair across everything the two checks are given.

``databaseTest.py`` carried the same dead predicate twice::

    num_of_d_bonds = sum([1 if x.order[0] == "D" and len(x.order) == 1 else 0 for x in atom.bonds.values()])

* ``kinetics_check_cd_atom_type`` -- the 140 kinetics families;
* ``general_check_cd_atom_type`` -- called from ``test_thermo``, ``test_solvation``,
  ``test_statmech`` and ``test_transport``, so it covers those four group databases.

A ``Group`` bond stores ``order`` as a list of *numbers*, so ``order[0] == "D"`` was never true,
``num_of_d_bonds`` was always ``0``, and neither check could fail for anything.

The two arms here are deliberately asymmetric, so this script cannot become a mirror that drifts
away from the thing it measures:

* the **as-shipped** arm is re-implemented below, because that code no longer exists in the file;
* the **repaired** arm calls the *real* checks by path import, so it cannot disagree with what the
  suite runs.

It reports *reach* alongside the complaint counts, because a zero from a check that has never run
is otherwise indistinguishable from a check that still does not run, and it refuses to report
anything until a negative control confirms the as-shipped predicate is silent exactly where the
real repaired check complains.

Usage:  python cd_blast_radius.py
"""

import importlib.util
import os
import sys
from collections import defaultdict

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.molecule import Group
from rmgpy.molecule.atomtype import ATOMTYPES

TARGET_ATOM_TYPES = [ATOMTYPES[x] for x in ["Cd", "CO", "CS", "Cdd"]]

REACH = defaultdict(int)


def real_checks():
    """Import ``test/database/databaseTest.py`` by path and return its ``TestDatabase``."""
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.normpath(os.path.join(here, "..", "..", "..", "test", "database", "databaseTest.py"))
    spec = importlib.util.spec_from_file_location("_databaseTest_for_blast_radius", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    print(f"repaired arm drives the real checks from {path}\n")
    return module.TestDatabase


CHECKS = real_checks()


def scan_as_shipped(entries, ignore, count_reach):
    """The check body as it stood before the repair, re-implemented because it is gone from HEAD."""
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
            # The dead line: a Group bond order is a number, so this is never true.
            num_of_d_bonds = sum(
                1 if bond.order[0] == "D" and len(bond.order) == 1 else 0
                for bond in atom.bonds.values()
            )
            if count_reach and any(
                len(b.order) == 1 and abs(2 - b.order[0]) < 1e-7 for b in atom.bonds.values()
            ):
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


def repaired_general_complains(group_name, group):
    """Drive the real `general_check_cd_atom_type`; True means it raised."""
    try:
        CHECKS.general_check_cd_atom_type(None, group_name, group)
        return False
    except ValueError:
        return True


def repaired_kinetics_complains(family_name, family):
    """Drive the real `kinetics_check_cd_atom_type`; True means it raised."""
    import types

    harness = types.SimpleNamespace(
        database=types.SimpleNamespace(kinetics=types.SimpleNamespace(families={family_name: family}))
    )
    try:
        CHECKS.kinetics_check_cd_atom_type(harness, family_name)
        return False
    except ValueError:
        return True


def negative_control():
    """The as-shipped predicate must be silent exactly where the real repaired check complains."""
    import types

    def entries_for(adjlist):
        return {"BAD": types.SimpleNamespace(
            label="BAD", item=Group().from_adjacency_list(adjlist), children=[])}

    misused_cd = "1 Cd u0 {2,D}\n2 O2d u0 {1,D}\n"
    cdd_subgraph = "1 * O2d u0 {2,D}\n2   Cdd u0 {1,D}\n"

    shipped_on_defect = scan_as_shipped(entries_for(misused_cd), [], count_reach=False)
    repaired_on_defect = repaired_general_complains(
        "MisusedCd", types.SimpleNamespace(entries=entries_for(misused_cd)))
    repaired_on_subgraph = repaired_general_complains(
        "O2d-Cdd", types.SimpleNamespace(entries=entries_for(cdd_subgraph)))

    print("--- negative control ---")
    print(f"  Cd double-bonded to O, typed Cd but not CO:")
    print(f"      as shipped complains:  {bool(shipped_on_defect)}   (must be False -- the bug)")
    print(f"      repaired complains:    {repaired_on_defect}   (must be True)")
    print(f"  Cdd carbon with one of its two double bonds drawn (the solvation O2d-Cdd shape):")
    print(f"      repaired complains:    {repaired_on_subgraph}  (must be False -- CO is impossible there)")
    ok = (not shipped_on_defect) and repaired_on_defect and (not repaired_on_subgraph)
    print(f"  control {'PASSES' if ok else 'FAILS'}\n")
    return ok


def kinetics_ignore(family):
    """`kinetics_check_cd_atom_type`'s ignore list. `general_check_cd_atom_type` has none."""
    if family.own_reverse:
        return []
    ignore = []
    for product in family.forward_template.products:
        ignore.append(product)
        ignore.extend(product.children)
    return ignore


def main():
    if not negative_control():
        print("negative control failed -- the numbers below mean nothing", file=sys.stderr)
        return 1

    database = RMGDatabase()
    database.load(settings["database.directory"], kinetics_families="all")

    shipped_total = 0
    repaired_failing = []

    n_families = len(database.kinetics.families)
    for name in sorted(database.kinetics.families):
        family = database.kinetics.families[name]
        shipped_total += len(scan_as_shipped(family.groups.entries, kinetics_ignore(family), True))
        if repaired_kinetics_complains(name, family):
            repaired_failing.append(f"kinetics/{name}")

    n_groups = 0
    for source_name, groups in [
        ("thermo", getattr(database.thermo, "groups", {}) or {}),
        ("solvation", getattr(database.solvation, "groups", {}) or {}),
        ("statmech", getattr(database.statmech, "groups", {}) or {}),
        ("transport", getattr(database.transport, "groups", {}) or {}),
    ]:
        for group_name, group in groups.items():
            n_groups += 1
            shipped_total += len(scan_as_shipped(group.entries, [], True))
            if repaired_general_complains(group_name, group):
                repaired_failing.append(f"{source_name}/{group_name}")

    print(f"kinetics families scanned:           {n_families}")
    print(f"group databases scanned:             {n_groups}")
    print(f"complaints, checks AS SHIPPED:       {shipped_total}")
    print(f"targets the REAL repaired checks fail: {len(repaired_failing)}")
    for name in repaired_failing:
        print(f"  {name}")
    print()
    print("reach (a zero above is only meaningful if these are not zero):")
    print(f"  Group entries scanned:             {REACH['entries_scanned']}")
    print(f"  atoms typed Cd/CO/CS/Cdd:          {REACH['target_typed_atoms']}")
    print(f"  ...of those, with a double bond:   {REACH['atoms_with_a_double_bond']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
