#!/usr/bin/env python3
"""
Q3-2: is a database loaded in the parent process visible to workers of a
multiprocessing.Pool started AFTER the load, under the 'fork' start method?

Parent load is in-process. Each pool worker's fingerprint is a
separate-process measurement (each worker is a distinct forked PID).
"""
import multiprocessing
import os

import rmgpy
from rmgpy.data.rmg import RMGDatabase, database as module_database, get_db

print('database.directory =', rmgpy.settings['database.directory'])
print('start_method =', multiprocessing.get_start_method())

THERMO_PATH = '/home/alon/Code/RMG-database-plasma/input/thermo'
assert os.path.isdir(THERMO_PATH), THERMO_PATH


def fingerprint(tag):
    import rmgpy.data.rmg as rmg_mod
    db = rmg_mod.database
    out = {
        'tag': tag,
        'pid': os.getpid(),
        'database_is_none': db is None,
    }
    if db is not None:
        out['id(database)'] = id(db)
        thermo_db = get_db('thermo')
        out['thermo_is_none'] = thermo_db is None
        if thermo_db is not None:
            libs = sorted(thermo_db.libraries.keys())
            out['library_keys'] = libs
            if 'primaryThermoLibrary' in thermo_db.libraries:
                out['n_entries_primaryThermoLibrary'] = len(
                    thermo_db.libraries['primaryThermoLibrary'].entries)
    return out


def worker_task(i):
    fp = fingerprint('worker')
    fp['task_index'] = i
    return fp


if __name__ == '__main__':
    # Load a real, small database in the PARENT, before creating the Pool.
    db = RMGDatabase()
    db.load_thermo(THERMO_PATH, thermo_libraries=['primaryThermoLibrary'])

    parent_fp_before_pool = fingerprint('parent-before-pool')
    print('PARENT fingerprint (before creating Pool):', parent_fp_before_pool)

    with multiprocessing.Pool(processes=3) as pool:
        results = pool.map(worker_task, range(6))

    for r in results:
        print('WORKER fingerprint:', r)

    parent_fp_after_pool = fingerprint('parent-after-pool')
    print('PARENT fingerprint (after pool closed/joined):', parent_fp_after_pool)

    # Decisive comparison
    same_id_as_any_worker = any(
        (not w['database_is_none']) and w.get('id(database)') == parent_fp_before_pool.get('id(database)')
        for w in results
    )
    same_libs = all(
        (not w['database_is_none']) and w.get('library_keys') == parent_fp_before_pool.get('library_keys')
        for w in results
    )
    print('RESULT any_worker_id_matches_parent_id =', same_id_as_any_worker)
    print('RESULT all_workers_library_keys_match_parent =', same_libs)
