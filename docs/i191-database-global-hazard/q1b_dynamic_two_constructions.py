#!/usr/bin/env python3
"""
Q1b - dynamic, in-process demonstration that a single process can construct
RMGDatabase twice via unmodified production code, the way
rmgpy/tools/simulate.py:83 and :88 do inside its per-reaction-system loop.

We drive rmgpy.rmg.main.RMG.load_database() twice on the same RMG object,
narrowing the load to a single kinetics family and empty library lists so the
load is cheap. This narrowing only affects how much data is loaded, not the
identity behaviour under test: RMGDatabase.__init__ unconditionally rebinds
the rmgpy.data.rmg.database module global regardless of what is subsequently
loaded into the new instance.

MEASUREMENT: in-process (single Python process, two constructions).
"""
import io
import logging
import sys

import rmgpy.data.rmg as rmg_data_module
from rmgpy import settings
from rmgpy.rmg.main import RMG

print("database.directory =", settings['database.directory'])
assert settings['database.directory'] == '/home/alon/Code/RMG-database-plasma/input', \
    "wrong database.directory - run from the i191 worktree root"

# Capture logging.warning output so we can show the 're-initializing' warning
log_stream = io.StringIO()
handler = logging.StreamHandler(log_stream)
handler.setLevel(logging.WARNING)
logging.getLogger().addHandler(handler)
logging.getLogger().setLevel(logging.INFO)

rmg = RMG()
rmg.database_directory = settings['database.directory']
# Minimal, cheap load: one family, no extra libraries/depositories.
rmg.thermo_libraries = []
rmg.transport_libraries = []
rmg.reaction_libraries = []
rmg.seed_mechanisms = []
rmg.kinetics_families = ['H_Abstraction']
rmg.kinetics_depositories = ['!training']
rmg.statmech_libraries = []
rmg.trimolecular_product_reversible = True
rmg.kinetics_estimator = 'rate rules'
rmg.kinetics_datastore = False
rmg.species_constraints = {}
rmg.forbidden_structures = []
rmg.binding_energies = None
rmg.solvent = None
rmg.adsorption_groups = 'adsorptionPt111'
rmg.verbose_comments = False

print("--- first load_database() call (mimics simulate.py:83, diffusion-limited liquid branch) ---")
rmg.load_database()
id_after_first = id(rmg_data_module.database)
print("id(rmgpy.data.rmg.database) after call #1 =", id_after_first)
print("id(rmg.database) after call #1            =", id(rmg.database))
assert rmg_data_module.database is rmg.database

print("--- second load_database() call (mimics simulate.py:87-88, uncertainty branch) ---")
rmg.verbose_comments = True  # exactly what simulate.py line 87 does before the reload
rmg.load_database()
id_after_second = id(rmg_data_module.database)
print("id(rmgpy.data.rmg.database) after call #2 =", id_after_second)
print("id(rmg.database) after call #2            =", id(rmg.database))
assert rmg_data_module.database is rmg.database

print()
print("RESULT: identity changed across the two calls in ONE process:",
      id_after_first != id_after_second)

logging.getLogger().removeHandler(handler)
captured = log_stream.getvalue()
print()
print("--- captured logging.warning output (filtered to the re-initialization line) ---")
for line in captured.splitlines():
    if 'Re-initializing' in line or 'already exists' in line:
        print(line)

sys.exit(0)
