#!/usr/bin/env python3
"""
Q1a follow-up. arkane/input.py:111 is cited as a Q1a static lead:

    rmg_database = get_db() or RMGDatabase()

Read as a fallback idiom: "use the existing global database, or make one if
there isn't one." Source-reading get_db() (rmgpy/data/rmg.py:247-282) shows
that when the global `database` is falsy (None), get_db() does not return a
falsy value for the `or` to fall through on -- it raises DatabaseError
unconditionally, for every value of `name` including the default ''. So the
`RMGDatabase()` right-hand side of the `or` can never execute via normal
short-circuit-on-falsy evaluation; it is unreachable except if get_db()
itself is later changed to return falsy instead of raising.

PART 1 of this script calls the production `database(...)` directive
DIRECTLY (bypassing arkane/input.py's load_input_file wrapper), in a fresh
process where rmgpy.data.rmg.database starts at None. This demonstrates
that database(...) alone, unprotected, crashes on a fresh process: the
RMGDatabase() half of the `or` is unreachable via short-circuit-on-falsy
evaluation, because get_db() raises DatabaseError rather than returning
anything falsy.

PART 2 resolves why real Arkane runs do not hit this crash: arkane/input.py
:653 (`load_input_file`) calls `load_necessary_databases()` BEFORE `exec()`
ever runs the input file's statements. `load_necessary_databases()`
(arkane/input.py:571-589) constructs an RMGDatabase() itself, in the
`except DatabaseError:` branch, populating db.statmech and db.transport --
and this construction rebinds the module global. By the time the input
file's own `database(...)` directive later executes, `get_db()` no longer
raises, so `or RMGDatabase()` is reached only vacuously (get_db() already
returns a truthy object; the right-hand side is never evaluated). Part 2
reproduces this exact ordering and shows database(...) succeeding the
second time, immediately after load_necessary_databases() has run.

MEASUREMENT: separate-process (single process containing two sequential
calls, matching load_input_file's real ordering: load_necessary_databases()
then the exec'd `database(...)` directive).
"""
import sys

import rmgpy.data.rmg as rmg_data_module
from rmgpy import settings

print("database.directory =", settings['database.directory'])
assert settings['database.directory'] == '/home/alon/Code/RMG-database-plasma/input', \
    "wrong database.directory - run from the i191 worktree root"

print("rmgpy.data.rmg.database at process start:", rmg_data_module.database)
assert rmg_data_module.database is None, \
    "harness invalid: some earlier import already constructed a database in this process"

from arkane.input import database, load_necessary_databases  # noqa: E402  (exact production functions)

print()
print("=== PART 1: database(...) called directly, unprotected, fresh process ===")
try:
    database(thermoLibraries=['primaryThermoLibrary'])
    print("database(...) call SUCCEEDED (unexpected)")
    print("rmgpy.data.rmg.database now:", rmg_data_module.database)
except Exception as e:
    print("database(...) call RAISED:", type(e).__name__, "-", e)
    print("Confirms: unprotected, on a fresh process, arkane/input.py:111's")
    print("`get_db() or RMGDatabase()` crashes with an uncaught DatabaseError")
    print("on the FIRST call -- get_db() raises rather than returning a")
    print("falsy value for `or` to fall through on, so the RMGDatabase()")
    print("fallback branch is unreachable via short-circuit evaluation.")

print()
print("=== PART 2: real ordering -- load_necessary_databases() first, matching")
print("    arkane/input.py:653's load_input_file(), THEN database(...) ===")
load_necessary_databases()
print("rmgpy.data.rmg.database after load_necessary_databases():", rmg_data_module.database)
try:
    database(thermoLibraries=['primaryThermoLibrary'])
    print("database(...) call SUCCEEDED")
    print("rmgpy.data.rmg.database now:", rmg_data_module.database)
    print()
    print("Confirms: in the real call order, load_necessary_databases() (called")
    print("by load_input_file BEFORE the input file's statements are exec'd)")
    print("constructs the global database first (in its except DatabaseError:")
    print("branch), so get_db() no longer raises by the time the input file's")
    print("own database(...) directive runs. The `or RMGDatabase()` at")
    print("arkane/input.py:111 is reached in this ordering, but vacuously --")
    print("get_db() already returns a truthy object, so RMGDatabase() (the")
    print("right-hand side) is never evaluated. Net: line 111's OWN")
    print("RMGDatabase() branch is unreachable in normal Arkane usage; not")
    print("because it doesn't need to be (get_db() truly can raise), but")
    print("because load_necessary_databases() -- a different construction")
    print("site -- always runs first and prevents get_db() from raising.")
except Exception as e:
    print("database(...) call RAISED (unexpected):", type(e).__name__, "-", e)
sys.exit(0)
