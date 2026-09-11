#!/usr/bin/env python3
"""
Q2 kinetics teeth test - separate process, same discipline as the thermo P3
experiment. Load DB_A with kinetics_families=['H_Abstraction']. Obtain a
reaction's family label string under DB_A ('H_Abstraction'). Construct DB_B
(the second construction) with a kinetics_families set that EXCLUDES
H_Abstraction. Then evaluate the production expression used at
rmgpy/rmg/model.py:1073 -- get_db('kinetics').families[reaction.family] --
against the label obtained under DB_A.

MEASUREMENT: separate-process.
"""
import sys

import rmgpy.data.rmg as rmg_data_module
from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase, get_db

print("database.directory =", settings['database.directory'])
assert settings['database.directory'] == '/home/alon/Code/RMG-database-plasma/input', \
    "wrong database.directory - run from the i191 worktree root"

# --- DB_A: only H_Abstraction ---
db_a = RMGDatabase()
db_a.load_kinetics(
    settings['database.directory'] + '/kinetics',
    reaction_libraries=[],
    seed_mechanisms=[],
    kinetics_families=['H_Abstraction'],
    kinetics_depositories=['!training'],
)
print("id(db_a) =", id(db_a))
print("db_a.kinetics.families.keys() =", sorted(db_a.kinetics.families.keys()))

# A reaction's family label string, obtained while DB_A is live -- this is
# exactly the shape of `reaction.family`, a plain str attribute that survives
# independently of which RMGDatabase object produced it.
reaction_family_label = 'H_Abstraction'
assert reaction_family_label in db_a.kinetics.families

# --- DB_B: the second construction, a DIFFERENT family set that excludes H_Abstraction ---
db_b = RMGDatabase()
db_b.load_kinetics(
    settings['database.directory'] + '/kinetics',
    reaction_libraries=[],
    seed_mechanisms=[],
    kinetics_families=['R_Recombination'],
    kinetics_depositories=['!training'],
)
print("id(db_b) =", id(db_b))
print("db_b.kinetics.families.keys() =", sorted(db_b.kinetics.families.keys()))

print()
print("id(rmgpy.data.rmg.database) =", id(rmg_data_module.database))
print("id(get_db('kinetics'))       =", id(get_db('kinetics')))
print("get_db('kinetics') is db_a.kinetics:", get_db('kinetics') is db_a.kinetics)
print("get_db('kinetics') is db_b.kinetics:", get_db('kinetics') is db_b.kinetics)

print()
print("--- reproducing rmgpy/rmg/model.py:1073  get_db('kinetics').families[reaction.family] ---")
try:
    family = get_db('kinetics').families[reaction_family_label]
    print("Lookup SUCCEEDED (unexpected):", family)
except KeyError as e:
    print("Lookup RAISED KeyError:", e)
    print("This is the production code path -- a reaction whose .family was set")
    print("while DB_A was live now fails a dict lookup against DB_B's families,")
    print("because DB_B never loaded H_Abstraction at all.")

sys.exit(0)
