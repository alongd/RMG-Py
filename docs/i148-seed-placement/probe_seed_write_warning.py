#!/usr/bin/env python
# encoding: utf-8

"""
I-148: the loud warning at seed-WRITE time, demonstrated on the real writer.

    python docs/i148-seed-placement/probe_seed_write_warning.py <scratch-dir>

``RMG.make_seed_mech`` is the function every RMG run calls, unconditionally, to
write the ``seed/`` directory. Since I-148 it warns, per entry, when a reaction
whose owner carries an electron-placement declaration is about to be written
WITHOUT the parsed provenance the seed reader recovers the owner from -- so a
seed that cannot restart announces itself in the run that writes it, not in the
run that fails to restart from it days later.

This probe drives that writer twice, on a minimal hand-built model (no database
load), through the real ``initialize_seed_mech`` + ``make_seed_mech`` path:

  P1  a reaction that carries its owner in ``Reaction.library`` -- the shipped
      shape, whose provenance the writer emits as an ``Originally from reaction
      library:`` line -- must produce NO warning;
  P2  the same reaction with its ``library`` provenance stripped (owner known
      only through ``family``, nothing for the writer to emit) must produce the
      SEED LOSES ELECTRON PLACEMENT warning, loudly, at write time.

The probe exits 0 only if P1 is silent, P2 fires, and both seeds were actually
written.

The resolved database directory is printed at the head of the run rather than
trusted from configuration (no database is loaded; the print is the campaign's
standing discipline).
"""

import logging
import os
import shutil
import sys

import rmgpy
from rmgpy import settings
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species

IONISATION = 'PlasmaElectronImpactIonization'
WARNING_TOKEN = 'SEED LOSES ELECTRON PLACEMENT'


class _Capture(logging.Handler):
    def __init__(self):
        super(_Capture, self).__init__(level=logging.WARNING)
        self.records = []

    def emit(self, record):
        self.records.append(record.getMessage())


def build_rmg(workdir, reaction):
    """A minimal RMG object carrying one core reaction, enough for the seed writer."""
    input_file = os.path.join(workdir, 'input.py')
    with open(input_file, 'w') as f:
        f.write('# placeholder input for the I-148 seed-write warning probe\n')
    rmg = RMG(input_file=input_file, output_directory=workdir)
    rmg.name = 'Seed'
    rmg.save_seed_to_database = False
    rmg.save_seed_modulus = -1
    rmg.restart = False
    rmg.reaction_model = CoreEdgeReactionModel()
    rmg.reaction_model.core.species = list(reaction.reactants) + list(reaction.products)
    rmg.reaction_model.core.reactions = [reaction]
    rmg.reaction_model.edge.species = []
    rmg.reaction_model.edge.reactions = []
    return rmg


def ionisation_reaction():
    lithium = Species(label='Li', molecule=[
        Molecule().from_adjacency_list('multiplicity 2\n1 Li u1 p0 c0')])
    cation = Species(label='Lip', molecule=[
        Molecule().from_adjacency_list('1 Li u0 p0 c+1')])
    kinetics = Arrhenius(A=(1.0, 'm^3/(mol*s)'), n=0, Ea=(0, 'kJ/mol'), T0=(1, 'K'),
                         comment='')
    return LibraryReaction(reactants=[lithium], products=[cation], electrons=1,
                           reversible=False, library=IONISATION, kinetics=kinetics)


def write_seed(workdir, reaction):
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    rmg = build_rmg(workdir, reaction)
    capture = _Capture()
    logging.getLogger().addHandler(capture)
    try:
        rmg.initialize_seed_mech()
        rmg.make_seed_mech()
    finally:
        logging.getLogger().removeHandler(capture)
    wrote = os.path.isfile(os.path.join(workdir, 'seed', 'seed', 'reactions.py'))
    warnings = [msg for msg in capture.records if WARNING_TOKEN in msg]
    return wrote, warnings


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    scratch = os.path.abspath(sys.argv[1])

    print('=' * 78)
    print('rmgpy              = {0}'.format(rmgpy.__file__))
    print('database.directory = {0}'.format(
        os.path.realpath(settings['database.directory'])))
    print('(no database is loaded by this probe)')
    print('=' * 78)

    failures = []

    print('\nP1 -- provenance present (Reaction.library = {0!r}): expect silence'.format(
        IONISATION))
    wrote, warnings = write_seed(os.path.join(scratch, 'p1'), ionisation_reaction())
    print('  seed written: {0};  {1} placement warning(s)'.format(wrote, len(warnings)))
    if not wrote:
        failures.append('P1: seed was not written')
    if warnings:
        failures.append('P1: warning fired for a preservable reaction')

    print('\nP2 -- provenance stripped (Reaction.library = None, owner only in '
          'Reaction.family): expect the warning')
    rxn = ionisation_reaction()
    rxn.library = None  # family still names the owner; the writer has nothing to emit
    wrote, warnings = write_seed(os.path.join(scratch, 'p2'), rxn)
    print('  seed written: {0};  {1} placement warning(s)'.format(wrote, len(warnings)))
    for msg in warnings:
        print('  WARNING: {0}'.format(msg.split('\n')[0]))
    if not wrote:
        failures.append('P2: seed was not written')
    if not warnings:
        failures.append('P2: the warning did NOT fire for an unpreservable reaction')

    print('\n' + '=' * 78)
    if failures:
        print('SEED-WRITE WARNING PROBE: {0} CHECK(S) FAILED'.format(len(failures)))
        for f in failures:
            print('  - {0}'.format(f))
        return 1
    print('SEED-WRITE WARNING PROBE: P1 silent, P2 fired, both seeds written')
    print('=' * 78)
    return 0


if __name__ == '__main__':
    sys.exit(main())
