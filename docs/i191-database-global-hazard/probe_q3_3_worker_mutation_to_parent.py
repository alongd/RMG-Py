#!/usr/bin/env python3
"""
Q3-3: does a worker's mutation of the inherited database reach the parent?

Two distinct mutations, tested in separate pools so they don't interfere:
  (a) in-place mutation: insert a sentinel key into get_db('thermo').libraries
  (b) rebind mutation: rmgpy.data.rmg.database = RMGDatabase() (a brand new object)

Parent-side checks after pool close/join are in-process measurements of
the PARENT's own memory; the mutation itself happens inside a forked
worker process (separate-process).
"""
import multiprocessing
import os

import rmgpy
import rmgpy.data.rmg as rmg_mod
from rmgpy.data.rmg import RMGDatabase, get_db

print('database.directory =', rmgpy.settings['database.directory'])
print('start_method =', multiprocessing.get_start_method())

THERMO_PATH = '/home/alon/Code/RMG-database-plasma/input/thermo'
assert os.path.isdir(THERMO_PATH), THERMO_PATH

SENTINEL_KEY = '__i191_sentinel_from_worker__'


def mutate_inplace(i):
    """Insert a sentinel key into the inherited thermo libraries dict."""
    thermo_db = get_db('thermo')
    thermo_db.libraries[SENTINEL_KEY] = thermo_db.libraries['primaryThermoLibrary']
    return {
        'pid': os.getpid(),
        'task_index': i,
        'sentinel_present_in_worker_immediately_after_mutation': SENTINEL_KEY in get_db('thermo').libraries,
        'id(database)': id(rmg_mod.database),
    }


def mutate_rebind(i):
    """Rebind the module-level global to a brand new RMGDatabase instance."""
    old_id = id(rmg_mod.database)
    rmg_mod.database = RMGDatabase()  # RMGDatabase.__init__ also does `global database; database = self`
    new_id = id(rmg_mod.database)
    return {
        'pid': os.getpid(),
        'task_index': i,
        'old_id(database)': old_id,
        'new_id(database)_in_worker': new_id,
        'rebind_took_effect_in_worker': old_id != new_id,
    }


if __name__ == '__main__':
    # ---- Experiment (a): in-place mutation of the inherited dict ----
    db = RMGDatabase()
    db.load_thermo(THERMO_PATH, thermo_libraries=['primaryThermoLibrary'])

    parent_id_before_a = id(rmg_mod.database)
    parent_has_sentinel_before_a = SENTINEL_KEY in get_db('thermo').libraries
    print('(a) PARENT id(database) before pool =', parent_id_before_a)
    print('(a) PARENT has sentinel BEFORE pool =', parent_has_sentinel_before_a)

    with multiprocessing.Pool(processes=3) as pool:
        results_a = pool.map(mutate_inplace, range(6))
    for r in results_a:
        print('(a) WORKER result:', r)

    parent_id_after_a = id(rmg_mod.database)
    parent_has_sentinel_after_a = SENTINEL_KEY in get_db('thermo').libraries
    print('(a) PARENT id(database) after pool =', parent_id_after_a)
    print('(a) PARENT has sentinel AFTER pool =', parent_has_sentinel_after_a)
    print('(a) RESULT in_place_worker_mutation_visible_to_parent =', parent_has_sentinel_after_a)

    # ---- Experiment (b): rebind the module-level global itself ----
    # Fresh parent-side database object again, so (a)'s sentinel doesn't confuse (b).
    db2 = RMGDatabase()
    db2.load_thermo(THERMO_PATH, thermo_libraries=['primaryThermoLibrary'])
    parent_id_before_b = id(rmg_mod.database)
    print('(b) PARENT id(database) before pool =', parent_id_before_b)

    with multiprocessing.Pool(processes=3) as pool:
        results_b = pool.map(mutate_rebind, range(6))
    for r in results_b:
        print('(b) WORKER result:', r)

    parent_id_after_b = id(rmg_mod.database)
    print('(b) PARENT id(database) after pool =', parent_id_after_b)
    print('(b) RESULT rebind_in_worker_visible_to_parent =', parent_id_after_b != parent_id_before_b)
