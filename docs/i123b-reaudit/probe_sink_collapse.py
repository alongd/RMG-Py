#!/usr/bin/env python
# encoding: utf-8

"""
Why does the generated lithium plasma mechanism have a source and no sink?

The I-123 deck runs to completion and writes a Chemkin file containing exactly one
reaction -- the electron-impact ionisation. The radiative recombination that the deck
also asks for never reaches the model. RMG.log records the reason in one informational
line:

    Adding reaction library PlasmaRadiativeRecombination to model edge...
    This library reaction was not new: [Lip](3) => Li(2)

This probe establishes, against the running code rather than by reading it, that:

  1. the two library reactions are physically distinct -- one creates a free electron,
     the other consumes one, and `Reaction.electrons` says so;
  2. `CoreEdgeReactionModel.make_new_reaction` nevertheless reports the second one as
     not new, so it is dropped;
  3. the drop is order-dependent -- whichever library is added second is the one that
     disappears -- which rules out "the sink was rejected on its merits";
  4. the duplicate check never consults `Reaction.electrons`, which is the only field
     that distinguishes them in the canonical representation.

Run from the RMG-Py worktree root so `rmgrc` resolves:

    python docs/i123b-reaudit/probe_sink_collapse.py
"""

import logging
import os

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.rmg.model import CoreEdgeReactionModel

logging.basicConfig(level=logging.CRITICAL)

IONISATION = 'PlasmaElectronImpactIonization'
RECOMBINATION = 'PlasmaRadiativeRecombination'

FAILURES = []


def check(label, condition, detail=''):
    mark = 'PASS' if condition else 'FAIL'
    if not condition:
        FAILURES.append(label)
    print('  [{0}] {1}{2}'.format(mark, label, ' -- {0}'.format(detail) if detail else ''))


def load_libraries():
    db = RMGDatabase()
    db.kinetics = db.kinetics or None
    from rmgpy.data.kinetics import KineticsDatabase

    db.kinetics = KineticsDatabase()
    db.kinetics.load_libraries(
        os.path.join(settings['database.directory'], 'kinetics', 'libraries'),
        libraries=[IONISATION, RECOMBINATION],
    )
    return db


def sole_reaction(db, name):
    library = db.kinetics.libraries[name]
    reactions = library.get_library_reactions()
    assert len(reactions) == 1, '{0} has {1} entries, expected 1'.format(name, len(reactions))
    return reactions[0]


def fresh_model(db):
    model = CoreEdgeReactionModel()
    model.kinetics_database = db.kinetics
    return model


def add(model, rxn):
    """Return (reaction, is_new) exactly as add_reaction_library_to_edge does."""
    return model.make_new_reaction(rxn, generate_thermo=False, generate_kinetics=False)


def main():
    print('=' * 78)
    print('database.directory = {0}'.format(settings['database.directory']))
    print('rmgpy              = {0}'.format(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))))))
    print('=' * 78)

    db = load_libraries()

    print('\nSTAGE 1 -- the two library reactions are physically distinct')
    ion = sole_reaction(db, IONISATION)
    rec = sole_reaction(db, RECOMBINATION)
    print('    ionisation    : {0}   electrons={1}'.format(ion, ion.electrons))
    print('    recombination : {0}   electrons={1}'.format(rec, rec.electrons))
    check('ionisation nets an electron GAIN', ion.electrons == 1, ion.electrons)
    check('recombination nets an electron LOSS', rec.electrons == -1, rec.electrons)
    check('the two differ in electron count', ion.electrons != rec.electrons,
          '{0} vs {1}'.format(ion.electrons, rec.electrons))
    check('their heavy-species halves are exact mirrors',
          [s.label for s in ion.reactants] == [s.label for s in rec.products]
          and [s.label for s in ion.products] == [s.label for s in rec.reactants],
          '{0} => {1}  vs  {2} => {3}'.format(
              [s.label for s in ion.reactants], [s.label for s in ion.products],
              [s.label for s in rec.reactants], [s.label for s in rec.products]))

    print('\nSTAGE 2 -- source first, then sink: the SINK is dropped')
    model = fresh_model(db)
    _, ion_new = add(model, sole_reaction(db, IONISATION))
    _, rec_new = add(model, sole_reaction(db, RECOMBINATION))
    check('ionisation entered the model', ion_new is True, 'is_new={0}'.format(ion_new))
    check('recombination was judged NOT NEW and dropped', rec_new is False,
          'is_new={0}'.format(rec_new))
    order_a = (ion_new, rec_new)

    print('\nSTAGE 3 -- sink first, then source: now the SOURCE is dropped')
    model = fresh_model(db)
    _, rec_new2 = add(model, sole_reaction(db, RECOMBINATION))
    _, ion_new2 = add(model, sole_reaction(db, IONISATION))
    check('recombination entered the model', rec_new2 is True, 'is_new={0}'.format(rec_new2))
    check('ionisation was judged NOT NEW and dropped', ion_new2 is False,
          'is_new={0}'.format(ion_new2))
    order_b = (rec_new2, ion_new2)
    check('the outcome is order-dependent, not a judgement on the sink',
          order_a == order_b == (True, False),
          'first-added survives in both orders')

    print('\nSTAGE 4 -- the verdict does not depend on Reaction.electrons at all')
    # `check_for_existing_reaction` must not be called directly with library-owned
    # Species objects: `are_identical_species_references` compares object identity, and
    # only `make_new_reaction` performs the conversion to the model's own Species. An
    # earlier revision of this probe called it directly, got False, and would have
    # reported a contradiction that was purely its own calling error. So drive the same
    # path STAGE 2 drives, and vary only the field under test.
    for forced in (0, -1, 12345):
        model = fresh_model(db)
        add(model, sole_reaction(db, IONISATION))
        probe = sole_reaction(db, RECOMBINATION)
        probe.electrons = forced
        _, is_new = add(model, probe)
        check('electrons={0}: sink still judged not new'.format(forced), is_new is False,
              'is_new={0}'.format(is_new))

    print('\n' + '=' * 78)
    if FAILURES:
        print('SINK-COLLAPSE PROBE: {0} CHECK(S) FAILED'.format(len(FAILURES)))
        for f in FAILURES:
            print('   - {0}'.format(f))
    else:
        print('SINK-COLLAPSE PROBE: the sink is discarded as a duplicate of the source.')
        print('The generated mechanism has an ionisation source and no recombination sink.')
    print('=' * 78)
    return 1 if FAILURES else 0


if __name__ == '__main__':
    raise SystemExit(main())
