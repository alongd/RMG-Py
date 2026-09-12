#!/usr/bin/env python3
"""
Q2 P3 - separate process, the experiment. Load DB_A (with primaryThermoLibrary).
Create methane. Do NOT resolve yet. Construct DB_B (without the library) - the
second construction, exactly the shape RMGDatabase() calls take in
rmgpy/tools/simulate.py's per-reaction-system loop. THEN resolve the species'
thermo via generate_thermo_data(spc), which internally calls get_db('thermo').

Which answer does it get: DB_A's (library, P1) or DB_B's (group additivity, P2)?

MEASUREMENT: separate-process (single process, but containing the two
RMGDatabase constructions under test - the species creation happens before
either evaluation, matching the brief's requested harness).
"""
import sys

import rmgpy.data.rmg as rmg_data_module
from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase, get_db
from rmgpy.species import Species
from rmgpy.thermo.thermoengine import generate_thermo_data

print("database.directory =", settings['database.directory'])
assert settings['database.directory'] == '/home/alon/Code/RMG-database-plasma/input', \
    "wrong database.directory - run from the i191 worktree root"

# --- construct DB_A, with the library ---
db_a = RMGDatabase()
db_a.load_thermo(settings['database.directory'] + '/thermo', thermo_libraries=['primaryThermoLibrary'])
print("id(db_a) =", id(db_a))

# --- create the species while DB_A is the global, but do NOT resolve thermo yet ---
spc = Species(smiles='C')
spc.generate_resonance_structures()
print("species created; spc.thermo =", spc.thermo)

# --- construct DB_B, the second construction, WITHOUT the library ---
db_b = RMGDatabase()
db_b.load_thermo(settings['database.directory'] + '/thermo', thermo_libraries=[])
print("id(db_b) =", id(db_b))

print()
print("id(rmgpy.data.rmg.database) =", id(rmg_data_module.database))
print("id(get_db('thermo'))        =", id(get_db('thermo')))
print("id(db_a.thermo)              =", id(db_a.thermo))
print("id(db_b.thermo)              =", id(db_b.thermo))
print("get_db('thermo') is db_a.thermo:", get_db('thermo') is db_a.thermo)
print("get_db('thermo') is db_b.thermo:", get_db('thermo') is db_b.thermo)

# --- now resolve, AFTER the second construction ---
thermo = generate_thermo_data(spc)
print()
print("H298 (kJ/mol) =", thermo.get_enthalpy(298) / 1000.0)
print("comment =", repr(thermo.comment))
print()
print("Matches P1 (library, -74.534 kJ/mol, 'Thermo library: primaryThermoLibrary')"
      " or P2 (group additivity, -74.894 kJ/mol, 'Thermo group additivity"
      " estimation: group(Cs-HHHH)')? See comment string above.")
sys.exit(0)
