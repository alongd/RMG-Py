#!/usr/bin/env python3
"""
Q2 P4 - separate process, control. Load DB_A (with primaryThermoLibrary),
create methane, resolve immediately, NO second construction. Must equal P1.

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

db_a = RMGDatabase()
db_a.load_thermo(settings['database.directory'] + '/thermo', thermo_libraries=['primaryThermoLibrary'])
print("id(db_a) =", id(db_a))

spc = Species(smiles='C')
spc.generate_resonance_structures()
thermo = generate_thermo_data(spc)
print("H298 (kJ/mol) =", thermo.get_enthalpy(298) / 1000.0)
print("comment =", repr(thermo.comment))
print()
print("Expect: H298 = -74.53377599999996, comment = 'Thermo library: primaryThermoLibrary' (must equal P1)")
sys.exit(0)
