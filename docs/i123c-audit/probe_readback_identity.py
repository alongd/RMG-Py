#!/usr/bin/env python
# encoding: utf-8

"""
What the reloaded representation MEANS, now that the Chemkin file reads back.

I-135 taught RMG's Chemkin reader the ``TDEP/`` auxiliary keyword, so the file RMG writes
for a plasma mechanism can now be read back by RMG. That closes the round trip in the
narrow sense -- the reaction and its electron-temperature dependence come back -- but it
opens a representation question the audit has to answer rather than gloss:

  the reloaded reaction carries the electrons as EXPLICIT SPECIES in the participant lists
  with a zero scalar count, where the canonical form carries a scalar and no electron
  participant.

Three questions, each measured:

  Q1  Does a reloaded reaction compare as THE SAME REACTION as its canonical original,
      under the identity predicate I-134 just repaired? If it does not, "read-back" does
      not mean "recovered the same object", and anyone treating a reloaded mechanism as
      interchangeable with the generated one is wrong.

  Q2  Do the two reloaded channels collapse into each other -- the I-134 defect, seen from
      the other side of the file? In the canonical form the two are mirrors and needed the
      placement declaration to be told apart. In the expanded form the electrons are
      ordinary participants, so the question is whether the expansion is self-distinguishing
      without any declaration at all.

  Q3  Can a reloaded plasma mechanism be fed back into a model? Charge balance, and the
      real ``CoreEdgeReactionModel.make_new_reaction`` path, on the reloaded objects.

Usage:  python docs/i123c-audit/probe_readback_identity.py <deck-output-dir>

Run from the repo root so ``rmgrc`` resolves; the resolved database directory is printed
at the head of the run rather than trusted from configuration.
"""

import os
import sys

from rmgpy import settings
from rmgpy import electron_placement
from rmgpy.chemkin import load_chemkin_file
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.electron_balance import get_electron_placement_counts
from rmgpy.rmg.model import CoreEdgeReactionModel

IONISATION = 'PlasmaElectronImpactIonization'
RECOMBINATION = 'PlasmaRadiativeRecombination'

FAILURES = []


def check(label, condition, detail=''):
    ok = bool(condition)
    print('  [{0}] {1}{2}'.format('PASS' if ok else 'FAIL', label,
                                  (' -- ' + detail) if detail else ''))
    if not ok:
        FAILURES.append(label)
    return ok


def note(text):
    print('  [..] {0}'.format(text))


def banner(text):
    print('\n' + '-' * 78)
    print(text)
    print('-' * 78)


def describe(rxn, prefix='    '):
    print('{0}{1}'.format(prefix, rxn))
    print('{0}  type            : {1}'.format(prefix, type(rxn).__name__))
    print('{0}  kinetics        : {1}'.format(prefix, type(rxn.kinetics).__name__))
    print('{0}  electrons scalar: {1}'.format(prefix, rxn.electrons))
    print('{0}  placement counts: {1}'.format(prefix, get_electron_placement_counts(rxn)))
    print('{0}  reactants       : {1}'.format(
        prefix, [(s.label, s.is_electron()) for s in rxn.reactants]))
    print('{0}  products        : {1}'.format(
        prefix, [(s.label, s.is_electron()) for s in rxn.products]))
    print('{0}  owner (family)  : {1}'.format(prefix, getattr(rxn, 'family', None)))


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    out_dir = sys.argv[1]
    chem_path = os.path.join(out_dir, 'chemkin', 'chem.inp')
    dict_path = os.path.join(out_dir, 'chemkin', 'species_dictionary.txt')

    print('=' * 78)
    import rmgpy.molecule.molecule as _m
    print('rmgpy              = {0}'.format(os.path.abspath(
        os.path.dirname(os.path.dirname(electron_placement.__file__)))))
    print('compiled module    = {0}'.format(_m.__file__))
    print('database.directory = {0}'.format(settings['database.directory']))
    print('resolved realpath  = {0}'.format(os.path.realpath(settings['database.directory'])))
    print('chemkin file       = {0}'.format(chem_path))
    print('=' * 78)

    banner('the canonical originals, straight out of the database')
    db = KineticsDatabase()
    db.load_libraries(os.path.join(settings['database.directory'], 'kinetics', 'libraries'),
                      libraries=[IONISATION, RECOMBINATION])
    canonical = {
        IONISATION: db.libraries[IONISATION].get_library_reactions()[0],
        RECOMBINATION: db.libraries[RECOMBINATION].get_library_reactions()[0],
    }
    for label, rxn in canonical.items():
        print('  {0}'.format(label))
        describe(rxn, prefix='      ')

    banner('the reloaded mechanism')
    species, reactions = load_chemkin_file(chem_path, dict_path)
    print('  read {0} species, {1} reactions'.format(len(species), len(reactions)))
    for rxn in reactions:
        describe(rxn)

    check('the file read back at all', bool(reactions), '{0} reaction(s)'.format(len(reactions)))
    check('every reloaded reaction carries a ZERO electron scalar',
          all(r.electrons == 0 for r in reactions),
          str([r.electrons for r in reactions]))
    check('every reloaded reaction carries the electron as an EXPLICIT participant',
          all(any(s.is_electron() for s in list(r.reactants) + list(r.products))
              for r in reactions))

    # ------------------------------------------------------------------ Q1
    banner('Q1 -- does a reloaded reaction compare as the same reaction as its original?')
    for label, original in canonical.items():
        for reloaded in reactions:
            same_dir = reloaded.is_isomorphic(original, either_direction=False)
            any_dir = reloaded.is_isomorphic(original, either_direction=True)
            note('{0:<32} vs reloaded {1:<34} same-direction={2} either-direction={3}'.format(
                label, str(reloaded), same_dir, any_dir))
    matches = [(label, str(r)) for label, o in canonical.items() for r in reactions
               if r.is_isomorphic(o, either_direction=True)]
    check('NO reloaded reaction compares equal to its canonical original',
          not matches, str(matches))
    note('this is the property, stated: the round trip recovers the CHEMISTRY in the '
         'expanded electron representation, not the canonical OBJECT. A reloaded '
         'mechanism and the model that produced it are not interchangeable under RMG\'s '
         'identity predicate, in either direction.')

    # ------------------------------------------------------------------ Q2
    banner('Q2 -- do the two reloaded channels collapse into each other?')
    if len(reactions) >= 2:
        for i in range(len(reactions)):
            for j in range(i + 1, len(reactions)):
                a, b = reactions[i], reactions[j]
                verdict = a.is_isomorphic(b, either_direction=True)
                note('{0}  vs  {1}  ->  is_isomorphic(either_direction=True) = {2}'.format(
                    a, b, verdict))
                check('reloaded pair {0}/{1} stays distinct'.format(i, j), verdict is False)
        note('in the EXPANDED form the electrons are ordinary participants, so the two '
             'channels differ in their participant lists alone -- the expansion is '
             'self-distinguishing and needs no placement declaration')
    else:
        note('fewer than two reactions in the reloaded file; nothing to compare')
        FAILURES.append('the reloaded file carries fewer than two reactions')

    # ------------------------------------------------------------------ Q3
    banner('Q3 -- can the reloaded mechanism be fed back into a model?')

    note('route 1: Reaction.is_balanced, the gate KineticsLibrary.load applies to every '
         'entry (rmgpy/data/kinetics/library.py:585) -- so this decides whether a '
         'reloaded plasma mechanism can be handed back to RMG as a library or a seed')
    balances = []
    for rxn in reactions:
        try:
            balanced = rxn.is_balanced()
        except Exception as exc:
            balanced = 'RAISED {0}: {1}'.format(type(exc).__name__, exc)
        counts = {}
        for side in ('reactants', 'products'):
            tally = {}
            for spc in getattr(rxn, side):
                for atom in spc.molecule[0].atoms:
                    tally[atom.element.symbol] = tally.get(atom.element.symbol, 0) + 1
            counts[side] = tally
        note('{0}'.format(rxn))
        note('    reactant elements {0}   product elements {1}   electrons scalar {2}'.format(
            counts['reactants'], counts['products'], rxn.electrons))
        note('    is_balanced -> {0}'.format(balanced))
        balances.append(balanced)
    note('is_balanced counts the electron as a conserved ELEMENT and then folds the '
         'scalar `electrons` into the net charges. The expanded form puts the electron '
         'in BOTH participant lists AND sets the scalar to zero, so the element tallies '
         'differ and the scalar cannot compensate.')

    note('route 2: CoreEdgeReactionModel.make_new_reaction on the reloaded objects')
    from rmgpy.data.rmg import RMGDatabase
    import rmgpy.data.rmg
    database = RMGDatabase()
    database.load(
        path=settings['database.directory'],
        thermo_libraries=['primaryThermoLibrary'],
        reaction_libraries=[],
        seed_mechanisms=[],
        kinetics_families=['Li_Abstraction'],
        kinetics_depositories=['training'],
        depository=False,
    )
    model = CoreEdgeReactionModel()
    accepted = []
    for rxn in reactions:
        try:
            new_rxn, is_new = model.make_new_reaction(rxn)
        except Exception as exc:
            note('make_new_reaction RAISED on {0}: {1}: {2}'.format(
                rxn, type(exc).__name__, exc))
            continue
        note('make_new_reaction({0}) -> is_new={1}'.format(rxn, is_new))
        if is_new:
            accepted.append(new_rxn)
    note('accepted {0} of {1} reloaded reactions'.format(len(accepted), len(reactions)))

    # ---- control: the same two routes on a NEUTRAL reloaded mechanism, so the plasma
    # ---- specific part of the answer is separated from the general one.
    control_dir = os.environ.get('I123C_NEUTRAL_CONTROL')
    if control_dir:
        note('control: the same two routes on a NEUTRAL Chemkin mechanism from {0}'.format(
            control_dir))
        c_spc, c_rxn = load_chemkin_file(
            os.path.join(control_dir, 'chemkin', 'chem.inp'),
            os.path.join(control_dir, 'chemkin', 'species_dictionary.txt'))
        note('    read {0} neutral reactions; owner labels {1}'.format(
            len(c_rxn), sorted({getattr(r, 'family', None) for r in c_rxn})))
        note('    is_balanced on all of them: {0}'.format(
            all(r.is_balanced() for r in c_rxn)))
        c_model = CoreEdgeReactionModel()
        c_ok = c_raised = 0
        c_err = ''
        for rxn in c_rxn[:5]:
            try:
                c_model.make_new_reaction(rxn)
                c_ok += 1
            except Exception as exc:
                c_raised += 1
                c_err = '{0}: {1}'.format(type(exc).__name__, str(exc)[:110])
        note('    make_new_reaction on the first 5: accepted {0}, raised {1}'.format(
            c_ok, c_raised))
        if c_err:
            note('    {0}'.format(c_err))
        check('CONTROL: is_balanced separates the two -- the NEUTRAL reload balances '
              'where the plasma reload does not',
              all(r.is_balanced() for r in c_rxn) and not any(
                  r.is_balanced() for r in reactions))
        check('CONTROL: make_new_reaction does NOT separate them -- the neutral reload '
              'is refused for the same Unclassified-owner reason', c_raised == 5,
              '{0} of 5 raised'.format(c_raised))
        note('so: the OWNER loss is a general property of any Chemkin-reloaded RMG '
             'mechanism and is not this campaign\'s to answer; the BALANCE loss is '
             'specific to the expanded electron representation and is this audit\'s '
             'finding')
    else:
        note('control skipped: set I123C_NEUTRAL_CONTROL to a non-plasma RMG output dir')

    banner('summary')
    print('  reloaded reactions        : {0}'.format(len(reactions)))
    print('  electron scalar on reload : {0}'.format([r.electrons for r in reactions]))
    print('  explicit e- participants  : {0}'.format(
        [sum(1 for s in list(r.reactants) + list(r.products) if s.is_electron())
         for r in reactions]))
    print('  identical to canonical    : {0}'.format('yes' if matches else 'NO'))
    print('  is_balanced on reload     : {0}'.format(balances))
    if FAILURES:
        print('\nFAILURES ({0}):'.format(len(FAILURES)))
        for label in FAILURES:
            print('  - {0}'.format(label))
        return 1
    print('\nall checks passed')
    return 0


if __name__ == '__main__':
    sys.exit(main())
