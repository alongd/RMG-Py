#!/usr/bin/env python3
"""
Q3-4: does a worker's mutation reach (a) a DIFFERENT worker process in the
same pool, or (b) a LATER task dispatched to the SAME worker process?

We submit far more tasks (40) than pool workers (4), with
maxtasksperchild=None (the default, same as production react.py), and have
every task report os.getpid() plus whether it sees a sentinel key inserted
by an earlier task. Each task, before checking, inserts its OWN
task-indexed sentinel key, so later tasks on the same pid accumulate
sentinels from all prior tasks handled by that pid, while tasks on a
different pid never see them (unless COW/fork gave them a common start,
which experiment 2/3 already showed is not retained after the pool ends
-- here we check DURING the pool's lifetime, across tasks).

This whole script is separate-process measurement (the interesting values
come from inside forked workers); the final parent-side summary is
in-process bookkeeping over what the workers reported.
"""
import multiprocessing
import os
import re

from rmgpy.data.rmg import RMGDatabase, get_db
import rmgpy

print('database.directory =', rmgpy.settings['database.directory'])
print('start_method =', multiprocessing.get_start_method())

THERMO_PATH = '/home/alon/Code/RMG-database-plasma/input/thermo'
assert os.path.isdir(THERMO_PATH)

N_TASKS = 40
N_WORKERS = 4


def task(i):
    thermo_db = get_db('thermo')
    libs = thermo_db.libraries

    # sentinel keys seen that were inserted by EARLIER tasks (possibly on
    # this same pid, possibly inherited via fork before any mutation)
    seen_before = sorted(k for k in libs.keys() if k.startswith('__sentinel_'))

    # this task's own sentinel, so later tasks (same or different pid) can detect it
    my_key = '__sentinel_task_{0:02d}__'.format(i)
    libs[my_key] = libs['primaryThermoLibrary']

    return {
        'task_index': i,
        'pid': os.getpid(),
        'sentinels_seen_before_my_mutation': seen_before,
    }


if __name__ == '__main__':
    db = RMGDatabase()
    db.load_thermo(THERMO_PATH, thermo_libraries=['primaryThermoLibrary'])

    with multiprocessing.Pool(processes=N_WORKERS, maxtasksperchild=None) as pool:
        results = list(pool.map(task, range(N_TASKS), chunksize=1))

    for r in results:
        print('TASK result:', r)

    # ---- Analysis ----
    by_pid = {}
    for r in results:
        by_pid.setdefault(r['pid'], []).append(r['task_index'])
    print('RESULT distinct_worker_pids =', sorted(by_pid.keys()))
    print('RESULT n_distinct_worker_pids =', len(by_pid))
    for pid, idxs in sorted(by_pid.items()):
        print('RESULT pid {0} handled task_indices (in submission order) = {1}'.format(pid, idxs))

    # (a) cross-pid leakage: did any task ever see a sentinel whose task_index
    # was handled by a DIFFERENT pid, at the time this task ran?
    idx_to_pid = {r['task_index']: r['pid'] for r in results}
    cross_pid_leak_examples = []
    for r in results:
        for skey in r['sentinels_seen_before_my_mutation']:
            src_idx = int(re.match(r'__sentinel_task_(\d+)__', skey).group(1))
            src_pid = idx_to_pid.get(src_idx)
            if src_pid is not None and src_pid != r['pid']:
                cross_pid_leak_examples.append(
                    (r['task_index'], r['pid'], src_idx, src_pid))
    print('RESULT cross_pid_leak_examples (task_idx, task_pid, source_task_idx, source_pid) =',
          cross_pid_leak_examples)
    print('RESULT any_cross_pid_leak =', bool(cross_pid_leak_examples))

    # (b) same-pid persistence: for each pid, did a later task on that pid see
    # sentinels from an earlier task on that SAME pid?
    same_pid_persistence_examples = []
    for r in results:
        for skey in r['sentinels_seen_before_my_mutation']:
            src_idx = int(re.match(r'__sentinel_task_(\d+)__', skey).group(1))
            src_pid = idx_to_pid.get(src_idx)
            if src_pid == r['pid'] and src_idx != r['task_index']:
                same_pid_persistence_examples.append(
                    (r['task_index'], r['pid'], src_idx))
    print('RESULT same_pid_persistence_examples (task_idx, pid, earlier_task_idx_same_pid) =',
          same_pid_persistence_examples)
    print('RESULT any_same_pid_persistence =', bool(same_pid_persistence_examples))
