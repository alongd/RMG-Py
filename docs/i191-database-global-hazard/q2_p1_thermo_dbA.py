#!/usr/bin/env python3
"""
Q2 P1 - separate process. Load DB_A (WITH primaryThermoLibrary) only.
Create methane, resolve thermo via rmgpy.thermo.thermoengine.generate_thermo_data.
Record H298 and comment.

MEASUREMENT: separate-process.
"""
import sys

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.species import Species
from rmgpy.thermo.thermoengine import generate_thermo_data

print("database.directory =", settings['database.directory'])
assert settings['database.directory'] == '/home/alon/Code/RMG-database-plasma/input', \
    "wrong database.directory - run from the i191 worktree root"

db = RMGDatabase()
db.load_thermo(settings['database.directory'] + '/thermo', thermo_libraries=['primaryThermoLibrary'])
print("id(db) =", id(db))

spc = Species(smiles='C')
spc.generate_resonance_structures()
thermo = generate_thermo_data(spc)
print("H298 (kJ/mol) =", thermo.get_enthalpy(298) / 1000.0)
print("comment =", repr(thermo.comment))
sys.exit(0)
