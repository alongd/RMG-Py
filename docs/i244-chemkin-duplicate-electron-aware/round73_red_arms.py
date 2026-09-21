"""
Show each new or repaired check RED by breaking exactly the thing it guards.

For every arm: patch ``rmgpy/chemkin.pyx``, rebuild the extension, run the named
tests, restore the file, rebuild. A check that cannot be made to fail is not a
check, and this is what makes that visible rather than asserted.

    python docs/i244-chemkin-duplicate-electron-aware/round73_red_arms.py <arm>

with ``<arm>`` one of the keys of ``ARMS``, or ``restore`` to put the file back.
Each arm prints the tests it expects to fail; the run's own pytest output is the
evidence, captured by the caller.
"""
import os
import shutil
import sys

SOURCE = 'rmgpy/chemkin.pyx'
#: A byte copy of the file, taken before the first patch and restored after each arm.
#: NOT ``git checkout -- <file>``: that restores from the INDEX, so on a worktree whose
#: repairs are not yet staged it silently reverts them to HEAD. It did, here, and cost a
#: full reapply -- and the four arms it ran were all measuring the unrepaired engine
#: while appearing to measure four different breakages.
BACKUP = SOURCE + '.round73-backup'

ARMS = {
    # Defect 1: the group key keys each side in list order again, so a pair
    # differing only in reactant order splits into two singletons and is cleared.
    'D1-ordered-key': (
        [("""    reactant_side = (tuple(sorted([id(spc) for spc in reaction.reactants])), reactant_electrons)
    product_side = (tuple(sorted([id(spc) for spc in reaction.products])), product_electrons)""",
          """    reactant_side = (tuple([id(spc) for spc in reaction.reactants]), reactant_electrons)
    product_side = (tuple([id(spc) for spc in reaction.products]), product_electrons)""")],
        ['TestTheAnswerBelongsToTheDeck::test_a_permuted_pair_is_not_cleared_into_a_deck_cantera_rejects'],
    ),
    # Defect 2: the render stores its answer on the reactions again, so it
    # outlives the deck it was computed for.
    'D2-leaking-recompute': (
        [("""    if check_for_duplicates:
        duplicate_flags = chemkin_duplicate_flags(reactions)
    else:
        duplicate_flags = [reaction.duplicate for reaction in reactions]

    if elements_in_use is None:""",
          """    if check_for_duplicates:
        mark_duplicate_reactions(reactions)
    duplicate_flags = [reaction.duplicate for reaction in reactions]

    if elements_in_use is None:""")],
        ['TestTheAnswerBelongsToTheDeck::test_a_render_leaves_every_reactions_flag_exactly_as_it_found_it',
         'TestTheAnswerBelongsToTheDeck::test_a_core_plus_edge_save_does_not_mark_a_core_only_cantera_entry'],
    ),
    # Defect 3: the render stops materializing its argument, so a generator is
    # exhausted by the keying pass and the writer emits nothing.
    'D3-consumed-generator': (
        [("""    reactions = list(reactions)
    # Check for duplicate
    if check_for_duplicates:
        duplicate_flags = chemkin_duplicate_flags(reactions)
    else:
        duplicate_flags = [reaction.duplicate for reaction in reactions]

    if elements_in_use is None:""",
          """    # Check for duplicate
    if check_for_duplicates:
        duplicate_flags = chemkin_duplicate_flags(reactions)
    else:
        duplicate_flags = [reaction.duplicate for reaction in reactions]

    if elements_in_use is None:""")],
        ['TestTheAnswerBelongsToTheDeck::test_a_generator_of_reactions_is_written_not_consumed'],
    ),
    # Defect 4: the authority is made a no-op -- it hands back the flags it was
    # given. Any check whose inputs already hold the values it asserts still
    # passes; only a check that requires a flag to CHANGE can see this.
    'D4-noop-authority': (
        [("""    reactions = list(reactions)
    groups = {}
    for index, reaction in enumerate(reactions):""",
          """    reactions = list(reactions)
    return [reaction.duplicate for reaction in reactions]
    groups = {}
    for index, reaction in enumerate(reactions):""")],
        ['TestDuplicateIsAGroupPredicate::test_a_cross_group_comparison_does_not_clear_another_groups_flags'],
    ),
}


def patch(arm):
    if not os.path.exists(BACKUP):
        shutil.copy2(SOURCE, BACKUP)
        print('backed up %s -> %s' % (SOURCE, BACKUP))
    text = open(SOURCE).read()
    for old, new in ARMS[arm][0]:
        if old not in text:
            raise SystemExit('arm %s: anchor not found in %s' % (arm, SOURCE))
        text = text.replace(old, new, 1)
    open(SOURCE, 'w').write(text)
    print('patched %s for arm %s' % (SOURCE, arm))


def restore():
    if not os.path.exists(BACKUP):
        raise SystemExit('no backup at %s; refusing to guess what to restore' % BACKUP)
    shutil.copy(BACKUP, SOURCE)          # copy, not copy2: copy2 preserves the backup's
    os.utime(SOURCE, None)               # mtime, and Cython then skips the rebuild, so the
    os.remove(BACKUP)                    # .so keeps running the broken arm after restore.
    print('restored %s from its byte backup' % SOURCE)


if __name__ == '__main__':
    what = sys.argv[1]
    if what == 'restore':
        restore()
    else:
        patch(what)
        print('expected RED:')
        for node in ARMS[what][1]:
            print('   ', node)
