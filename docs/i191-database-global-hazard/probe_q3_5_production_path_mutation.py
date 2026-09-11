#!/usr/bin/env python3
"""
Q3-5: does the REAL production reaction-generation code path -- the one
p.map(_react_species_star, ...) in rmgpy/rmg/react.py actually calls --
mutate any shared, module-global-reachable database state in a worker, in a
way that would (a) be silently lost, or (b) silently persist and leak to
later tasks on the same worker pid (per Q3-4's finding that same-pid
persistence is real for synthetic mutations)?

Source reading done before writing this probe (see report for line numbers):
  - family.py:1845 generate_reactions()'s own docstring states it "does not
    estimate the kinetics of these reactions at this time" -- and grep
    confirms estimate_kinetics/estimate_kinetics_using_rate_rules are never
    called from generate_reactions, _generate_reactions, or
    generate_reactions_from_families. They ARE called from
    family.get_kinetics, which is only invoked from
    rmgpy/rmg/model.py:1095 CoreEdgeReactionModel.generate_kinetics -- i.e.
    AFTER react()/react_all() returns to the PARENT, not inside a worker.
    So KineticsRules.estimate_kinetics (rules.py:361), and any mutation it
    might perform on self.rules.entries, is not reachable from the
    multiprocessing worker path at all.
  - database.py:509-553 generate_reactions_from_families() reads
    self.families[...], calls react_molecules -> family.generate_reactions
    (read-only against family.forward_template/reverse_template/groups in
    the code inspected), then for own_reverse families calls
    family.add_reverse_attribute(rxn) UNCONDITIONALLY for every reaction
    from such a family (not just on error).
  - family.py:1894-1981 add_reverse_attribute() contains exactly one write
    to a shared, persistent-looking attribute: `self.forbidden =
    ForbiddenStructures()` (family.py:1940), but ONLY on the rare
    zero-reverse-reactions-found error/debug branch, and it is unconditionally
    restored via try/finally (family.py:1941-1947) before the function
    returns, even if an exception occurs.

This probe empirically checks the *normal* (non-error) path: real families,
real reactants, run generate_reactions_from_families(...) (the exact call
react_species() makes) many times across a Pool with maxtasksperchild=None
(production's default) and more tasks than workers, and checks whether
family.forbidden identity/content, family.rules.entries count, and
family.groups.entries count ever change, either within a worker or across
tasks on the same worker pid.

Parent-side database load is in-process. Each task's generate_reactions_from_families
call and its fingerprint are separate-process measurements (inside a forked
worker). The final comparison/analysis is in-process (parent), over values
reported back by workers.
"""
import multiprocessing
import os

import rmgpy
from rmgpy.data.rmg import RMGDatabase, get_db
from rmgpy.molecule import Molecule

print('database.directory =', rmgpy.settings['database.directory'])
print('start_method =', multiprocessing.get_start_method())

KINETICS_FAMILIES = ['R_Recombination', 'H_Abstraction']
N_TASKS = 24
N_WORKERS = 4


def fingerprint(tag, i=None):
    kdb = get_db('kinetics')
    fp = {'tag': tag, 'pid': os.getpid()}
    if i is not None:
        fp['task_index'] = i
    for fam_label in KINETICS_FAMILIES:
        fam = kdb.families[fam_label]
        fp[fam_label] = {
            'id(forbidden)': id(fam.forbidden),
            'forbidden_n_entries': len(fam.forbidden.entries) if fam.forbidden is not None else None,
            'rules_n_entries': len(fam.rules.entries) if fam.rules is not None else None,
            'groups_n_entries': len(fam.groups.entries) if fam.groups is not None else None,
        }
    return fp


def task(i):
    before = fingerprint('worker-before-call', i)

    # This is exactly the call react_species() makes (rmgpy/rmg/react.py:94),
    # against real molecules and real families, using the production API.
    ch3 = Molecule(smiles='[CH3]')
    ch4 = Molecule(smiles='C')
    reactions_a = get_db('kinetics').generate_reactions_from_families(
        [ch3, ch3], only_families=KINETICS_FAMILIES)
    reactions_b = get_db('kinetics').generate_reactions_from_families(
        [ch4], only_families=KINETICS_FAMILIES)

    after = fingerprint('worker-after-call', i)

    return {
        'pid': os.getpid(),
        'task_index': i,
        'n_reactions_generated': len(reactions_a) + len(reactions_b),
        'before': before,
        'after': after,
    }


if __name__ == '__main__':
    db = RMGDatabase()
    db.load(
        rmgpy.settings['database.directory'],
        thermo_libraries=['primaryThermoLibrary'],
        kinetics_families=KINETICS_FAMILIES,
        reaction_libraries=[],
        seed_mechanisms=[],
        surface=False,
    )

    parent_fp_before = fingerprint('parent-before-pool')
    print('PARENT fingerprint before pool:', parent_fp_before)

    with multiprocessing.Pool(processes=N_WORKERS, maxtasksperchild=None) as pool:
        results = list(pool.map(task, range(N_TASKS), chunksize=1))

    parent_fp_after = fingerprint('parent-after-pool')
    print('PARENT fingerprint after pool:', parent_fp_after)

    for r in results:
        print('TASK result: pid={0} task_index={1} n_reactions_generated={2}'.format(
            r['pid'], r['task_index'], r['n_reactions_generated']))

    # ---- Analysis ----
    # NOTE: fingerprint() includes 'tag'/'task_index'/'pid' bookkeeping keys
    # that trivially differ between calls; only compare the per-family
    # sub-dicts (the actually content-bearing measurement) below.

    # (1) Within-task: did the call change this worker's own fingerprint?
    within_task_changes = []
    for r in results:
        for fam_label in KINETICS_FAMILIES:
            b = r['before'][fam_label]
            a = r['after'][fam_label]
            if b != a:
                within_task_changes.append((r['task_index'], r['pid'], fam_label, b, a))
    print('RESULT within_task_fingerprint_changes =', within_task_changes)
    print('RESULT any_within_task_change =', bool(within_task_changes))

    # (2) Across tasks on the SAME pid: did a later task's "before" family
    # fingerprint differ from an earlier task's "before" family fingerprint
    # on the same pid?
    by_pid_in_order = {}
    for r in results:
        by_pid_in_order.setdefault(r['pid'], []).append(r)
    same_pid_drift = []
    for pid, task_results in by_pid_in_order.items():
        first_before = {fam: task_results[0]['before'][fam] for fam in KINETICS_FAMILIES}
        for r in task_results[1:]:
            this_before = {fam: r['before'][fam] for fam in KINETICS_FAMILIES}
            if this_before != first_before:
                same_pid_drift.append(
                    (pid, task_results[0]['task_index'], r['task_index'], first_before, this_before))
    print('RESULT same_pid_drift_examples (pid, first_task_idx, drifted_task_idx, first_fp, drifted_fp) =',
          same_pid_drift)
    print('RESULT any_same_pid_drift =', bool(same_pid_drift))

    # (3) Parent before vs after pool: any change visible in the parent's own
    # family fingerprints (excluding the 'tag'/'pid' bookkeeping keys) from
    # all this worker activity?
    parent_before_fams = {fam: parent_fp_before[fam] for fam in KINETICS_FAMILIES}
    parent_after_fams = {fam: parent_fp_after[fam] for fam in KINETICS_FAMILIES}
    parent_changed = parent_before_fams != parent_after_fams
    print('RESULT parent_fingerprint_changed_after_pool =', parent_changed)

    # (4) Sanity: did we actually generate any reactions (i.e. is this a live
    # exercise of the code path, not a no-op)?
    total_reactions = sum(r['n_reactions_generated'] for r in results)
    print('RESULT total_reactions_generated_across_all_tasks =', total_reactions)
    print('RESULT probe_is_live (total_reactions_generated > 0) =', total_reactions > 0)
