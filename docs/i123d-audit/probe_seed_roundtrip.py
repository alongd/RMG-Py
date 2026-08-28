#!/usr/bin/env python
# encoding: utf-8

"""
The seed mechanism RMG writes for a plasma model, and what happens when it is read back.

    python docs/i123d-audit/probe_seed_roundtrip.py <deck-output-dir>

Every RMG run writes ``seed/`` and a ready-to-run ``restart_from_seed.py`` beside its
Chemkin output, unconditionally and without being asked. That is a third serialisation of
the same model, alongside Chemkin and Cantera, and no pass of this audit has exercised it.
Unlike Chemkin it keeps the CANONICAL representation -- the original ``VoronovEIArrhenius``
and ``BadnellRRArrhenius`` objects, with the scalar ``electrons`` metadata -- so on the face
of it, it is the lossless one.

What it does NOT keep is the reaction's OWNER. The seed library is written under the label
``seed`` (``restart`` once loaded), not under the label of the library the reaction came
from, and the electron-placement declaration is keyed on exactly that label. So every
plasma reaction that passes through a seed loses its declaration and silently falls back
to the net-derived rule.

This probe measures that, in four steps, on the reaction objects themselves rather than
through a whole RMG run -- so the arithmetic is visible and not just the crash:

  S1  what the seed file carries: the rate law classes and the scalar electron counts.
  S2  the owner each reloaded seed reaction reports, and whether that owner is a key in
      the placement registry.
  S3  the placement counts the canonical library entry gets and the placement counts the
      same chemistry gets after the seed round trip, side by side. Where those differ,
      state the consequence in the same line.
  S4  the identity verdict between the seed copy and the canonical library copy of the
      SAME reaction, per channel. A pair that fails to compare equal is a reaction RMG
      will add twice.

The two whole-run measurements this probe explains -- restart with the libraries declared
(the file RMG itself writes) and restart with them removed -- are run separately and their
logs are beside this one in ``evidence/``.

The resolved database directory is printed at the head of the run rather than trusted from
configuration.
"""

import os
import sys

from rmgpy import settings
from rmgpy import electron_placement
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.library import KineticsLibrary
from rmgpy.electron_balance import get_electron_placement_counts

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


def owner_of(rxn):
    """The label the placement registry is keyed on, for whatever kind of reaction."""
    for attribute in ('family', 'library'):
        value = getattr(rxn, attribute, None)
        if value is not None:
            return getattr(value, 'label', value)
    return None


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    deck = os.path.abspath(sys.argv[1])
    seed_dir = os.path.join(deck, 'seed', 'seed')

    print('=' * 78)
    print('database.directory = {0}'.format(settings['database.directory']))
    print('resolved           = {0}'.format(
        os.path.realpath(settings['database.directory'])))
    print('seed directory     = {0}'.format(seed_dir))
    print('=' * 78)

    if not os.path.isfile(os.path.join(seed_dir, 'reactions.py')):
        print('no seed mechanism at {0}'.format(seed_dir))
        return 2

    banner('S0 -- the canonical originals, from the database')
    db = KineticsDatabase()
    db.load_libraries(os.path.join(settings['database.directory'], 'kinetics', 'libraries'),
                      libraries=[IONISATION, RECOMBINATION])
    canonical = {}
    for label in (IONISATION, RECOMBINATION):
        rxn = db.libraries[label].get_library_reactions()[0]
        canonical[label] = rxn
        note('{0:<32} {1:<18} electrons={2:+d}  owner={3}  placement={4}'.format(
            label, str(rxn), rxn.electrons, owner_of(rxn),
            get_electron_placement_counts(rxn)))

    banner('S1 -- what the seed file carries')
    seed_db = KineticsDatabase()
    seed_db.load_libraries(os.path.join(deck, 'seed'), libraries=['seed'])
    seed_lib = seed_db.libraries.get('seed')
    check('the seed mechanism loads at all', seed_lib is not None)
    if seed_lib is None:
        return 1
    seed_rxns = seed_lib.get_library_reactions()
    check('the seed carries both charge channels', len(seed_rxns) == 2,
          '{0} reaction(s)'.format(len(seed_rxns)))
    for rxn in seed_rxns:
        note('{0:<28} kinetics={1:<22} electrons={2:+d}'.format(
            str(rxn), type(rxn.kinetics).__name__, rxn.electrons))
    check('the seed kept the ORIGINAL rate laws, not a Chemkin reduction',
          sorted(type(r.kinetics).__name__ for r in seed_rxns)
          == ['BadnellRRArrhenius', 'VoronovEIArrhenius'],
          str(sorted(type(r.kinetics).__name__ for r in seed_rxns)))
    check('the seed kept the scalar electron counts',
          sorted(r.electrons for r in seed_rxns) == [-1, 1],
          str(sorted(r.electrons for r in seed_rxns)))

    banner('S2 -- the owner each reloaded seed reaction reports')
    registry = electron_placement.FAMILY_ELECTRON_PLACEMENT
    note('FAMILY_ELECTRON_PLACEMENT keys = {0}'.format(sorted(registry)))
    owners = []
    for rxn in seed_rxns:
        owner = owner_of(rxn)
        owners.append(owner)
        note('{0:<28} owner={1!r}   declared in the registry: {2}'.format(
            str(rxn), owner, owner in registry))
    check('no reloaded seed reaction names an owner the registry knows',
          all(o not in registry for o in owners), str(owners))
    note('the originating library label is NOT recoverable from the seed entry: '
         'it survives only as free text in longDesc, which nothing parses')

    banner('S3 -- placement counts, canonical versus after the seed round trip')
    note('{0:<24} {1:<12} {2:<12} {3}'.format(
        'channel', 'canonical', 'via seed', 'consequence'))
    pairs = []
    for label, want in ((IONISATION, (1, 2)), (RECOMBINATION, (1, 0))):
        canon = canonical[label]
        # pair each seed reaction with the canonical one it is the same chemistry as
        match = None
        for rxn in seed_rxns:
            if (rxn.reactants[0].is_isomorphic(canon.reactants[0])
                    and rxn.products[0].is_isomorphic(canon.products[0])):
                match = rxn
                break
        if match is None:
            check('the seed carries the {0} channel'.format(label), False)
            continue
        canon_counts = get_electron_placement_counts(canon)
        seed_counts = get_electron_placement_counts(match)
        same = canon_counts == seed_counts
        note('{0:<24} {1:<12} {2:<12} {3}'.format(
            label.replace('Plasma', ''), str(canon_counts), str(seed_counts),
            'declaration preserved by luck' if same
            else 'DECLARATION LOST -- falls back to the net-derived rule'))
        check('the canonical {0} entry still declares {1}'.format(label, want),
              canon_counts == want, str(canon_counts))
        pairs.append((label, canon, match, canon_counts, seed_counts))

    note('')
    note('the fallback rule for an undeclared owner is: electrons<0 -> (-electrons, 0), '
         'electrons>0 -> (0, electrons). For the recombination (electrons=-1) that gives '
         '(1, 0), which is what the declaration says anyway, so nothing changes. For the '
         'ionisation (electrons=+1) it gives (0, 1) where the declaration says (1, 2): '
         'the electron that is CONSUMED disappears from the reactant side, which is '
         'exactly the incident-order information a net scalar cannot carry.')

    banner('S4 -- identity: is the seed copy the same reaction as the canonical copy?')
    for label, canon, match, canon_counts, seed_counts in pairs:
        same_direction = canon.is_isomorphic(match)
        either = canon.is_isomorphic(match, either_direction=True)
        note('{0:<32} same-direction={1}  either-direction={2}'.format(
            label, same_direction, either))
        if canon_counts == seed_counts:
            check('{0}: the seed copy IS the canonical reaction'.format(label), either,
                  'placement agrees, so identity agrees')
        else:
            check('{0}: the seed copy is NOT the canonical reaction'.format(label),
                  not either,
                  'placement disagrees {0} vs {1}, so RMG adds this channel TWICE when '
                  'the seed and the library are both declared'.format(
                      canon_counts, seed_counts))

    print('\n' + '=' * 78)
    if FAILURES:
        print('SEED ROUND TRIP: {0} CHECK(S) FAILED'.format(len(FAILURES)))
        for f in FAILURES:
            print('  - {0}'.format(f))
        return 1
    print('SEED ROUND TRIP: every step measured as described above')
    print('=' * 78)
    return 0


if __name__ == '__main__':
    sys.exit(main())
