#!/usr/bin/env python3
"""
Q1c - does get_db('transport')/get_db('statmech') raise DatabaseError, or
return None, when the global database exists but those attributes were never
populated? And what does that mean for arkane/input.py:571-588
load_necessary_databases(), whose guard is
    try:
        get_db('transport')
        get_db('statmech')
    except DatabaseError:
        ... construct and load them ...

MEASUREMENT: in-process (single Python process; no second process needed -
this is about get_db's return-vs-raise behaviour on one already-constructed,
partially-populated database).
"""
import sys

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase, get_db
from rmgpy.exceptions import DatabaseError

print("database.directory =", settings['database.directory'])
assert settings['database.directory'] == '/home/alon/Code/RMG-database-plasma/input', \
    "wrong database.directory - run from the i191 worktree root"

# Construct a database carrying ONLY thermo (mimics e.g. arkane/input.py:111
# `get_db() or RMGDatabase()` having been reached by some earlier path that
# loaded thermo but never called load_necessary_databases' transport/statmech
# loaders).
db = RMGDatabase()
db.load_thermo(settings['database.directory'] + '/thermo', thermo_libraries=['primaryThermoLibrary'])
print("db.transport is", db.transport)
print("db.statmech is", db.statmech)

print()
print("--- calling get_db('transport') on a database whose .transport is None ---")
try:
    result = get_db('transport')
    print("get_db('transport') returned (no exception):", result)
    raised = False
except DatabaseError as e:
    print("get_db('transport') RAISED DatabaseError:", e)
    raised = True

print()
print("--- calling get_db('statmech') on a database whose .statmech is None ---")
try:
    result2 = get_db('statmech')
    print("get_db('statmech') returned (no exception):", result2)
    raised2 = False
except DatabaseError as e:
    print("get_db('statmech') RAISED DatabaseError:", e)
    raised2 = True

print()
print("RESULT: get_db returns None silently (does not raise) when the attribute")
print("is None but the global `database` object itself exists:", (not raised) and (not raised2))

print()
print("--- consequence for arkane/input.py:load_necessary_databases() ---")
print("Reproducing its exact guard logic against this partially-populated database:")
try:
    get_db('transport')
    get_db('statmech')
    print("No exception raised -> the `except DatabaseError` branch that would")
    print("construct+load transport/statmech databases is NEVER entered.")
    print("transport/statmech remain None: db.transport =", db.transport,
          " db.statmech =", db.statmech)
except DatabaseError:
    print("DatabaseError raised -> transport/statmech WOULD be loaded.")

sys.exit(0)
