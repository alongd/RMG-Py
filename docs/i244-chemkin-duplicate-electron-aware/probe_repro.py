#!/usr/bin/env python3
"""
I-244 reproduction probe: does the Chemkin writer's duplicate check collapse two
reactions whose electron placements differ per side?

Drives the REAL export path -- the shipped libraries through the real loader, a
real CoreEdgeReactionModel so that make_new_species unifies isomorphic species
into single objects, then render_chemkin_file with check_for_duplicates=True.
The unification is the part that matters: at library-load time the two libraries
hold distinct Species objects and the identity comparison in
mark_duplicate_reaction is False for unrelated reasons.

Run from the worktree root with rmgrc pointing at RMG-database-plasma.
"""
import logging
import sys

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
import rmgpy.data.rmg
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.chemkin import render_chemkin_file
from rmgpy.electron_balance import get_electron_placement_counts

logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format='%(levelname)s %(message)s')

LIBRARIES = ['PlasmaElectronImpactIonization', 'PlasmaRadiativeRecombination']
THERMO = ['LithiumPrimaryThermo', 'LithiumAdditionalThermo',
          'primaryThermoLibrary', 'electrocatThermo']


def main():
    db = RMGDatabase()
    db.load(
        path=settings['database.directory'],
        thermo_libraries=THERMO,
        transport_libraries=[],
        reaction_libraries=LIBRARIES,
        seed_mechanisms=[],
        kinetics_families=['Li_Abstraction', 'Li_Addition_MultipleBond'],
        kinetics_depositories=['training'],
        depository=False,
        statmech_libraries=[],
        solvation=False,
        surface=False,
    )
    rmgpy.data.rmg.database = db

    model = CoreEdgeReactionModel()
    # The electron pseudo-species, exactly as the deck declares it, through the
    # same unification the library species go through.
    from rmgpy.species import Species
    from rmgpy.molecule import Molecule
    electron = Species(label='e-', molecule=[Molecule().from_adjacency_list('1 e u1 p0 c-1')])
    electron, _ = model.make_new_species(electron, label='e-')
    electron.thermo = db.thermo.get_thermo_data(electron)
    model.add_species_to_edge(electron)

    for lib in LIBRARIES:
        model.add_reaction_library_to_edge(lib)

    # A SECOND owner for the same ionisation channel, carrying the same net
    # electron count but NO placement declaration -- so the net-derived rule
    # gives it (0, 1) where PlasmaElectronImpactIonization's declaration gives
    # (1, 2). This is the shape I-148 recorded from a seed round trip ("the
    # ionisation channel's (1, 2) collapsed to (0, 1), the restarted core
    # carried the channel twice") and the shape any user library that writes the
    # channel without declaring its incident order produces. The object is built
    # exactly the way KineticsLibrary.load builds one: LibraryReaction over the
    # model's already-unified Species objects, with `electrons` assigned after
    # construction, which is the line library.py:628 runs.
    from rmgpy.data.kinetics.library import LibraryReaction
    from rmgpy.kinetics import Arrhenius

    li = [s for s in model.edge.species if str(s).startswith('[Li](')][0]
    lip = [s for s in model.edge.species if str(s).startswith('[Lip](')][0]
    undeclared = LibraryReaction(
        reactants=[li],
        products=[lip],
        library='UndeclaredIonizationLibrary',
        # First order as the net rule writes it -- one reactant, no incident
        # electron. The units have to match the equation the writer emits.
        kinetics=Arrhenius(A=(1e10, 's^-1'), n=0, Ea=(0, 'kcal/mol')),
        reversible=False,
        duplicate=False,
    )
    undeclared.electrons = 1
    model.edge.reactions.append(undeclared)

    reactions = model.edge.reactions
    species = model.edge.species

    print('\n==== loaded reactions, BEFORE any duplicate marking ====')
    for rxn in reactions:
        print('  {0!s:40s} class={1:16s} family={2!r} library={3!r} '
              'reversible={4} electrons={5:+d} placement={6} duplicate={7}'.format(
                  rxn, type(rxn).__name__, getattr(rxn, 'family', None),
                  getattr(rxn, 'library', None), rxn.reversible,
                  getattr(rxn, 'electrons', 0) or 0,
                  get_electron_placement_counts(rxn), rxn.duplicate))

    print('\n==== species object identity after unification ====')
    for spc in species:
        print('  {0!s:12s} id={1}'.format(spc, id(spc)))
    for rxn in reactions:
        print('  {0!s:30s} reactant ids={1} product ids={2}'.format(
            rxn, [id(s) for s in rxn.reactants], [id(s) for s in rxn.products]))

    print('\n==== render_chemkin_file(check_for_duplicates=True) ====')
    for i, spc in enumerate(species):
        spc.index = i + 1
    text = render_chemkin_file(species, reactions, verbose=False,
                               check_for_duplicates=True)
    print('\n---- REACTIONS block of the emitted deck ----')
    in_rxns = False
    for line in text.splitlines():
        if line.strip().startswith('REACTIONS'):
            in_rxns = True
        if in_rxns:
            print(line)

    print('\n==== duplicate flags AFTER marking ====')
    for rxn in reactions:
        print('  {0!s:40s} duplicate={1} placement={2}'.format(
            rxn, rxn.duplicate, get_electron_placement_counts(rxn)))

    n_dup = sum(1 for r in reactions if r.duplicate)
    print('\nRESULT: {0} of {1} reactions marked DUPLICATE'.format(n_dup, len(reactions)))

    # ---- the Cantera leg ----------------------------------------------------
    # yaml_cantera2 never calls mark_duplicate_reactions itself -- it reads
    # reaction.duplicate, which the Chemkin writer above has just set. Measure
    # that rather than assume it: render the Cantera YAML from the SAME reaction
    # objects and show which entries carry `duplicate: true`.
    from rmgpy.yaml_cantera2 import generate_cantera_data
    from rmgpy.rmg.model import ReactionModel

    elements = ReactionModel(species=species).get_elements()
    data = generate_cantera_data(species, reactions, elements_in_use=elements,
                                 is_plasma=True)
    print('\n==== Cantera YAML reaction entries ====')
    for entry in data['reactions']:
        print('  equation={0!r:55s} duplicate={1}'.format(
            entry.get('equation'), entry.get('duplicate', False)))


if __name__ == '__main__':
    main()
