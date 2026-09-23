#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
Round 112 end-to-end close gate (I-221, brief verifier item 10).

A really generated electron-impact reaction, taken the whole way:
family generation -> ``copy()`` -> pair reconstruction -> ``saturate_radicals`` (through
``saturate_for_estimation``, the production caller) -> atom typing -> model construction
(``CoreEdgeReactionModel.make_new_reaction``) -> group-additivity thermo of every model
species. Run for three reactants: Li (the campaign's own ionisation), CH3 (an ordinary
organic radical) and metastable argon ``Ar u2 p3`` (the species I-221 is about, whose
saturated form no atom type owns).

The gate: no stage raises a bare ``AtomTypeError``. A ``SaturatedStructureError`` naming the
species is the intended answer for metastable argon and is recorded, not failed; so is a
``DatabaseError`` from the thermo group tree (Ar+ has no group; no libraries are loaded). Pairs are
checked by POSITION on each side, before and after the copy and after the model rebuilt the
lists, on the reactor-facing view where the electron occurs on both sides and twice on one.

Exit 0 when the gate holds, 1 otherwise.
"""

import logging
import os
import sys

logging.disable(logging.INFO)

import rmgpy.data.rmg
from rmgpy import settings
from rmgpy.data.base import saturate_for_estimation
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.thermo import ThermoDatabase
from rmgpy.electron_placement import resolve_electron_placement
from rmgpy.exceptions import AtomTypeError, DatabaseError, SaturatedStructureError
from rmgpy.molecule import Molecule
from rmgpy.reaction import pair_occurrences
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species

IONISATION = 'Plasma_Electron_Impact_Ionization'
REACTANTS = {
    'Li': Molecule(smiles='[Li]'),
    'CH3': Molecule(smiles='[CH3]'),
    'Ar metastable': Molecule().from_adjacency_list('multiplicity 3\n1 Ar u2 p3 c0\n'),
}

failures = []


def check(name, ok, detail=''):
    print('{0}  {1}{2}'.format('ok  ' if ok else 'FAIL', name, (' -- ' + detail) if detail else ''))
    if not ok:
        failures.append(name)


def stage(label, name, function):
    """Run one stage; a bare AtomTypeError is the gate failure."""
    try:
        value = function()
        check('{0}: {1}'.format(label, name), True)
        return value
    except SaturatedStructureError as exc:
        check('{0}: {1} refused by name (SaturatedStructureError)'.format(label, name), True,
              str(exc).splitlines()[0][:100])
        return 'refused'
    except DatabaseError as exc:
        # The thermo group tree has no node for the structure (Ar+ has none; production
        # supplies it from a library, and this probe loads no libraries on purpose, so
        # the estimator itself is reached). Loud and named -- not the gate's failure.
        check('{0}: {1} refused by the database (DatabaseError)'.format(label, name), True,
              str(exc).splitlines()[0][:100])
        return 'refused'
    except AtomTypeError as exc:
        check('{0}: {1} raised a BARE AtomTypeError'.format(label, name), False, str(exc)[:160])
        return 'bare'


def molecules(species_or_molecule):
    if isinstance(species_or_molecule, Species):
        return species_or_molecule.molecule
    return [species_or_molecule]


def main():
    database_directory = settings['database.directory']
    kinetics = KineticsDatabase()
    kinetics.load_families(path=os.path.join(database_directory, 'kinetics', 'families'),
                           families=[IONISATION])
    family = kinetics.families[IONISATION]
    thermo = ThermoDatabase()
    thermo.load(os.path.join(database_directory, 'thermo'), libraries=[])

    saved = rmgpy.data.rmg.database
    database = RMGDatabase()
    database.kinetics = kinetics
    database.thermo = thermo
    rmgpy.data.rmg.database = database
    electron = Species(label='e', molecule=[Molecule(smiles='e')])
    try:
        for label, molecule in REACTANTS.items():
            print('\n== {0}'.format(label))
            species = Species(label=label, molecule=[molecule.copy(deep=True)])
            species.generate_resonance_structures()
            # 1. family generation, default arguments
            generated = family.generate_reactions([species.molecule[0]])
            check('{0}: the family generates a reaction'.format(label), len(generated) >= 1,
                  '{0} reaction(s)'.format(len(generated)))
            if not generated:
                continue
            reaction = generated[0]
            reaction.kinetics = family.get_kinetics(
                reaction, template_labels=reaction.template, degeneracy=reaction.degeneracy)[0][0]
            reaction.ensure_species()
            view = resolve_electron_placement(
                reaction, [electron] + list(reaction.reactants) + list(reaction.products))
            view.generate_pairs()
            print('    view: {0}'.format(view))
            electrons_left = sum(1 for s in view.reactants if s is electron)
            electrons_right = sum(1 for s in view.products if s is electron)
            check('{0}: the view carries the electron on both sides, twice on the right'.format(label),
                  (electrons_left, electrons_right) == (1, 2),
                  'placement ({0}, {1})'.format(electrons_left, electrons_right))

            # 2. copy() and 3. pair reconstruction
            before = pair_occurrences(view.pairs, view.reactants, view.products)
            copied = view.copy()
            after = pair_occurrences(copied.pairs, copied.reactants, copied.products)
            check('{0}: copy() carries every pair to the same positions'.format(label),
                  before == after, '{0} -> {1}'.format(before, after))
            originals = {id(s) for s in list(view.reactants) + list(view.products)}
            aliased = [m for pair in copied.pairs for m in pair if id(m) in originals]
            check('{0}: no copied pair member is an object of the original'.format(label),
                  not aliased, '{0} aliased'.format(len(aliased)))

            # 4. saturate_radicals, through its production caller, on every structure
            for structure in list(copied.reactants) + list(copied.products):
                for m in molecules(structure):
                    if m.is_radical():
                        stage(label, 'saturating {0}'.format(m.to_smiles()),
                              lambda m=m: saturate_for_estimation(m, 'thermodynamic data'))

            # 5. atom typing of every structure the copy holds
            for structure in list(copied.reactants) + list(copied.products):
                for m in molecules(structure):
                    stage(label, 'typing {0}'.format(m.to_smiles()),
                          lambda m=m: m.update_atomtypes())

            # 6. model construction. The model takes the canonical (family) reaction -- the
            # electron lives in `electrons`, not in the lists; the view above is what the
            # reactor sees. So the model is fed a copy of the canonical reaction.
            reaction.generate_pairs()
            canonical = pair_occurrences(reaction.pairs, reaction.reactants, reaction.products)
            forward = reaction.copy()
            check('{0}: copy() keeps the canonical reaction\'s class and electrons'.format(label),
                  type(forward) is type(reaction) and forward.electrons == reaction.electrons,
                  '{0}, electrons {1}'.format(type(forward).__name__, forward.electrons))
            model = CoreEdgeReactionModel()
            result = stage(label, 'make_new_reaction',
                           lambda: model.make_new_reaction(forward, generate_thermo=False,
                                                           generate_kinetics=False))
            if result in ('refused', 'bare'):
                continue
            built = result[0]
            built_occurrences = pair_occurrences(built.pairs, built.reactants, built.products)
            check('{0}: the model reaction keeps the pairs at the canonical positions'.format(label),
                  built_occurrences == canonical, '{0} -> {1}'.format(canonical, built_occurrences))
            check('{0}: the model\'s pair members are the model\'s own species'.format(label),
                  all(pair[0] in built.reactants and pair[1] in built.products
                      for pair in built.pairs))

            # 7. group-additivity thermo of every heavy model species
            for s in list(built.reactants) + list(built.products):
                if s.molecule[0].is_electron():
                    continue
                stage(label, 'thermo for {0}'.format(s.molecule[0].to_smiles()),
                      lambda s=s: thermo.get_thermo_data(s))
    finally:
        rmgpy.data.rmg.database = saved

    print('\n{0} failure(s)'.format(len(failures)))
    for name in failures:
        print('  ' + name)
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
