#!/usr/bin/env python
# encoding: utf-8

"""
I-148: the EDGE seed's electron-placement round trip, driven through the real
writer and the real reader.

    python docs/i148-seed-placement/probe_edge_seed_roundtrip.py <scratch-dir>

The audit named the edge seed as an artifact it did not exercise. No shipped
plasma deck can produce a CHARGED edge seed to exercise it with: the plasma
solver promotes the cation to the core at t=0 "to avoid singularity"
(measured on the toleranceMoveToCore=1e6 variant of the i123-integration deck,
whose final edge is 0 species / 0 reactions), so the charged channels never
remain in the edge. This probe therefore drives the same code path the deck
cannot reach: it places the two charged library reactions in the EDGE of a
hand-built model, writes the seed with the real writer (``RMG.make_seed_mech``,
which serialises the edge through its own loop with its own provenance
``try/except``), reloads ``seed_edge`` with the real reader
(``KineticsDatabase.load_libraries`` + ``get_library_reactions``), and measures
the placement counts off the reloaded objects.

  E1  the written seed_edge carries both channels with their rate-law classes
      and scalar electron counts intact (asserted on the CLASS, not the number).
  E2  each reloaded edge reaction reports family = the container label and
      library = the original owner, and the placement counts equal the
      canonical declarations -- (1, 2) for the ionisation, (1, 0) for the
      recombination. Before I-148 the ionisation read (0, 1).
  E3  no SEED LOSES ELECTRON PLACEMENT warning fired: the edge writer's
      provenance line is sufficient for both.

Exit 0 only if all checks hold.
"""

import logging
import os
import shutil
import sys

import rmgpy
from rmgpy import settings
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.electron_balance import get_electron_placement_counts, get_placement_owner
from rmgpy.kinetics import BadnellRRArrhenius, VoronovEIArrhenius
from rmgpy.molecule import Molecule
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species

IONISATION = 'PlasmaElectronImpactIonization'
RECOMBINATION = 'PlasmaRadiativeRecombination'
WARNING_TOKEN = 'SEED LOSES ELECTRON PLACEMENT'

FAILURES = []


def check(label, condition, detail=''):
    ok = bool(condition)
    print('  [{0}] {1}{2}'.format('PASS' if ok else 'FAIL', label,
                                  (' -- ' + detail) if detail else ''))
    if not ok:
        FAILURES.append(label)
    return ok


class _Capture(logging.Handler):
    def __init__(self):
        super(_Capture, self).__init__(level=logging.WARNING)
        self.records = []

    def emit(self, record):
        self.records.append(record.getMessage())


def species():
    lithium = Species(label='Li', molecule=[
        Molecule().from_adjacency_list('multiplicity 2\n1 Li u1 p0 c0')])
    cation = Species(label='Lip', molecule=[
        Molecule().from_adjacency_list('1 Li u0 p0 c+1')])
    return lithium, cation


def charged_channels():
    lithium, cation = species()
    # The exact shipped parameters, copied from the seed the i123-integration
    # deck writes (Voronov 1997 Z=3 N=3; Badnell 2006 Z=3 N=2) -- copied, not
    # authored: this probe measures serialisation, not kinetics.
    ionisation = LibraryReaction(
        reactants=[lithium], products=[cation], electrons=1, reversible=False,
        library=IONISATION,
        kinetics=VoronovEIArrhenius(A=(1.39e-07, 'cm^3/(molecule*s)'), P=0,
                                    X=0.438, K=0.41, dE=(5.4, 'eV'), electrons=1,
                                    Tmin=(11604.5, 'K'), Tmax=(2.3209e+08, 'K'),
                                    comment='Voronov (1997) e-impact ionization fit, Z=3, N=3'))
    recombination = LibraryReaction(
        reactants=[cation], products=[lithium], electrons=-1, reversible=False,
        library=RECOMBINATION,
        kinetics=BadnellRRArrhenius(A=(8.7e-12, 'cm^3/(molecule*s)'), B=0.364,
                                    T0=(147, 'K'), T1=(7.153e+06, 'K'), C=0.1508,
                                    T2=(715400, 'K'), electrons=-1,
                                    Tmin=(10, 'K'), Tmax=(1e+07, 'K'),
                                    comment='Badnell (2006) RR fit, Z=3, N=2'))
    return lithium, cation, ionisation, recombination


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    scratch = os.path.abspath(sys.argv[1])
    workdir = os.path.join(scratch, 'edge-seed')
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)

    print('=' * 78)
    print('rmgpy              = {0}'.format(rmgpy.__file__))
    print('database.directory = {0}'.format(
        os.path.realpath(settings['database.directory'])))
    print('(the database itself is not loaded; the reactions are built in place)')
    print('=' * 78)

    lithium, cation, ionisation, recombination = charged_channels()

    input_file = os.path.join(workdir, 'input.py')
    with open(input_file, 'w') as f:
        f.write('# placeholder input for the I-148 edge-seed probe\n')
    rmg = RMG(input_file=input_file, output_directory=workdir)
    rmg.name = 'Seed'
    rmg.save_seed_to_database = False
    rmg.save_seed_modulus = -1
    rmg.restart = False
    rmg.reaction_model = CoreEdgeReactionModel()
    rmg.reaction_model.core.species = [lithium]
    rmg.reaction_model.core.reactions = []
    rmg.reaction_model.edge.species = [cation]
    rmg.reaction_model.edge.reactions = [ionisation, recombination]

    capture = _Capture()
    logging.getLogger().addHandler(capture)
    try:
        rmg.initialize_seed_mech()
        rmg.make_seed_mech()
    finally:
        logging.getLogger().removeHandler(capture)

    print('\nE1 -- what the written seed_edge carries')
    edge_file = os.path.join(workdir, 'seed', 'seed_edge', 'reactions.py')
    check('seed_edge/reactions.py was written', os.path.isfile(edge_file))
    if FAILURES:
        return 1

    seed_db = KineticsDatabase()
    seed_db.load_libraries(os.path.join(workdir, 'seed'), libraries=['seed_edge'])
    edge_lib = seed_db.libraries.get('seed_edge')
    check('the edge seed loads', edge_lib is not None)
    if edge_lib is None:
        return 1
    rxns = edge_lib.get_library_reactions()
    check('the edge seed carries both charge channels', len(rxns) == 2,
          '{0} reaction(s)'.format(len(rxns)))
    kin_classes = sorted(type(r.kinetics).__name__ for r in rxns)
    check('the ORIGINAL rate-law classes survived (asserted on the class)',
          kin_classes == ['BadnellRRArrhenius', 'VoronovEIArrhenius'], str(kin_classes))
    check('the scalar electron counts survived',
          sorted(r.electrons for r in rxns) == [-1, 1],
          str(sorted(r.electrons for r in rxns)))

    print('\nE2 -- owners and placements off the reloaded edge objects')
    want = {1: (IONISATION, (1, 2)), -1: (RECOMBINATION, (1, 0))}
    for rxn in rxns:
        owner_want, counts_want = want[rxn.electrons]
        counts = get_electron_placement_counts(rxn)
        print('  {0!s:<14} family={1!r} library={2!r} owner={3!r} placement={4}'.format(
            rxn, getattr(rxn, 'family', None), getattr(rxn, 'library', None),
            get_placement_owner(rxn), counts))
        check('{0}: family is the container label'.format(owner_want),
              rxn.family == 'seed_edge', repr(rxn.family))
        check('{0}: library is the original owner'.format(owner_want),
              rxn.library == owner_want, repr(rxn.library))
        check('{0}: placement is the canonical declaration'.format(owner_want),
              counts == counts_want, '{0} (canonical {1})'.format(counts, counts_want))

    print('\nE3 -- the seed-write warning stayed silent (provenance sufficed)')
    warnings = [msg for msg in capture.records if WARNING_TOKEN in msg]
    check('no placement warning fired for the edge entries', not warnings,
          '{0} warning(s)'.format(len(warnings)))

    print('\n' + '=' * 78)
    if FAILURES:
        print('EDGE SEED ROUND TRIP: {0} CHECK(S) FAILED'.format(len(FAILURES)))
        for f in FAILURES:
            print('  - {0}'.format(f))
        return 1
    print('EDGE SEED ROUND TRIP: written by the real writer, read by the real '
          'reader, placement preserved')
    print('=' * 78)
    return 0


if __name__ == '__main__':
    sys.exit(main())
