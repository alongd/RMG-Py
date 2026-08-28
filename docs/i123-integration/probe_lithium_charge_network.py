#!/usr/bin/env python
# encoding: utf-8

"""
I-123 -- the lithium charge network, end to end on the integrated tree.

The backlog exists to deliver one thing: a neutral lithium feed that reaches its first
cation through electron-impact ionisation, and a cation that has a radiative
recombination sink. Every piece of that has been tested on its own branch. This probe is
the first time all of them are present at once, so it walks the network stage by stage
and prints what each stage produced, rather than asserting a single pass/fail at the end.

Stages, for BOTH channels:

  1. LOAD       -- the two kinetics libraries load out of the integration database.
  2. BALANCE    -- each canonical reaction closes its charge gap via Reaction.is_balanced,
                   which only became a real check on i093-isbalanced.
  3. PLACEMENT  -- resolve_electron_placement turns the canonical `A => A+` form into the
                   reactor-facing `A + e- => A+ + 2 e-`, from the library's declaration in
                   FAMILY_ELECTRON_PLACEMENT. Both library labels must be registry keys,
                   which is what i116 (ionisation) and i119 (recombination) added.
  4. REACTOR    -- one PlasmaReactor holds BOTH reactions at once, resolves placement
                   itself, and evaluates each at its own rate law's electron-temperature
                   rate. Source and sink in the same model is the integration claim.

Run from the RMG-Py integration worktree root so rmgrc resolves:

    conda activate rmg_env
    python docs/i123-integration/probe_lithium_charge_network.py
"""

import os
import sys

from rmgpy import settings
from rmgpy import electron_placement
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial

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


def flat_thermo():
    """Constant-Cp NASA polynomial. The reactor needs thermo to pack a model; both
    reactions here are irreversible, so no Keq is ever formed and the numbers below do
    not depend on this being the real thermochemistry."""
    coeffs = [2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    return NASA(
        polynomials=[NASAPolynomial(coeffs=coeffs, Tmin=(100, 'K'), Tmax=(5000, 'K')),
                     NASAPolynomial(coeffs=coeffs, Tmin=(5000, 'K'), Tmax=(20000, 'K'))],
        Tmin=(100, 'K'), Tmax=(20000, 'K'))


def main():
    print('=' * 78)
    print('database.directory = {0}'.format(settings['database.directory']))
    print('rmgpy              = {0}'.format(
        os.path.dirname(os.path.dirname(os.path.abspath(electron_placement.__file__)))))
    print('=' * 78)

    # ---------------------------------------------------------------- stage 1: load
    print('\nSTAGE 1 -- LOAD both kinetics libraries out of the integration database')
    db = KineticsDatabase()
    db.load_libraries(os.path.join(settings['database.directory'], 'kinetics', 'libraries'),
                      libraries=[IONISATION, RECOMBINATION])
    check('ionisation library loaded', IONISATION in db.libraries)
    check('recombination library loaded', RECOMBINATION in db.libraries)

    ion_rxns = db.libraries[IONISATION].get_library_reactions()
    rec_rxns = db.libraries[RECOMBINATION].get_library_reactions()
    check('ionisation has exactly one entry', len(ion_rxns) == 1, str(len(ion_rxns)))
    check('recombination has exactly one entry', len(rec_rxns) == 1, str(len(rec_rxns)))

    ionisation, recombination = ion_rxns[0], rec_rxns[0]
    print('    ionisation    : {0}   electrons={1}  kinetics={2}'.format(
        ionisation, ionisation.electrons, type(ionisation.kinetics).__name__))
    print('    recombination : {0}   electrons={1}  kinetics={2}'.format(
        recombination, recombination.electrons, type(recombination.kinetics).__name__))

    check('ionisation declares a net electron GAIN', ionisation.electrons == +1,
          str(ionisation.electrons))
    check('recombination declares a net electron LOSS', recombination.electrons == -1,
          str(recombination.electrons))
    check('the neutral feed and the cation are the same two species in both channels',
          ionisation.reactants[0].is_isomorphic(recombination.products[0])
          and ionisation.products[0].is_isomorphic(recombination.reactants[0]),
          'source produces exactly what the sink consumes')

    # ------------------------------------------------------------ stage 2: balance
    print('\nSTAGE 2 -- BALANCE: each canonical reaction closes its own charge gap')
    check('ionisation is_balanced', ionisation.is_balanced())
    check('recombination is_balanced', recombination.is_balanced())

    # ---------------------------------------------------------- stage 3: placement
    print('\nSTAGE 3 -- PLACEMENT: the resolver restores the electron participants')
    registry = electron_placement.FAMILY_ELECTRON_PLACEMENT
    print('    FAMILY_ELECTRON_PLACEMENT = {0}'.format(registry))
    check('the ionisation library is a declared owner', IONISATION in registry)
    check('the recombination library is a declared owner', RECOMBINATION in registry)
    check('ionisation declares one electron in, two out',
          registry.get(IONISATION) == (1, 2), str(registry.get(IONISATION)))
    check('recombination declares one electron in, none out',
          registry.get(RECOMBINATION) == (1, 0), str(registry.get(RECOMBINATION)))

    electron = Species(label='e', molecule=[Molecule().from_adjacency_list('1 e u0 p0 c-1')])

    def resolve(rxn):
        return electron_placement.resolve_electron_placement(
            rxn, list(rxn.reactants) + list(rxn.products) + [electron])

    def electron_sides(view):
        return (sum(1 for s in view.reactants if s.is_electron()),
                sum(1 for s in view.products if s.is_electron()))

    ion_view, rec_view = resolve(ionisation), resolve(recombination)
    print('    ionisation view    : {0}  ->  reactant/product electrons {1}, order {2}'
          .format(ion_view, electron_sides(ion_view), len(ion_view.reactants)))
    print('    recombination view : {0}  ->  reactant/product electrons {1}, order {2}'
          .format(rec_view, electron_sides(rec_view), len(rec_view.reactants)))
    check('ionisation resolves to Li + e- => Li+ + 2 e-',
          electron_sides(ion_view) == (1, 2), str(electron_sides(ion_view)))
    check('recombination resolves to Li+ + e- => Li + hv',
          electron_sides(rec_view) == (1, 0), str(electron_sides(rec_view)))
    check('both views carry no residual electron metadata',
          ion_view.electrons == 0 and rec_view.electrons == 0)
    check('both resolved views are second order in the reactor',
          len(ion_view.reactants) == 2 and len(rec_view.reactants) == 2)
    check('resolution did not mutate the canonical reactions',
          len(ionisation.reactants) == 1 and len(recombination.reactants) == 1)

    # ------------------------------------------------------------ stage 4: reactor
    print('\nSTAGE 4 -- REACTOR: one PlasmaReactor holds the source AND the sink')
    from rmgpy.solver.plasma import PlasmaReactor

    Te = 10000.0
    lithium = ionisation.reactants[0]
    cation = ionisation.products[0]
    # The sink's species objects must BE the source's, or the reactor sees four species
    # where the chemistry has three. This is the join the whole backlog is about.
    recombination.reactants[0] = cation
    recombination.products[0] = lithium
    for species in (lithium, cation, electron):
        species.thermo = flat_thermo()

    reactor = PlasmaReactor(
        T=(1000, 'K'), P=(1, 'bar'), Te=(Te, 'K'),
        initial_mole_fractions={lithium: 1.0, cation: 1e-12, electron: 1e-12},
        termination=[])
    reactor.initialize_model(core_species=[lithium, cation, electron],
                            core_reactions=[ionisation, recombination],
                            edge_species=[], edge_reactions=[])

    check('the reactor accepted both reactions', len(reactor.kf) == 2, str(len(reactor.kf)))
    check('the reactor located the electron', reactor.electron_index == 2,
          str(reactor.electron_index))

    kf_ion = ionisation.kinetics.get_rate_coefficient_electron_temp(Te)
    kf_rec = recombination.kinetics.get_rate_coefficient_electron_temp(Te)
    print('    kf[source] = {0:.6e}   expected (Voronov) = {1:.6e}'.format(reactor.kf[0], kf_ion))
    print('    kf[sink]   = {0:.6e}   expected (Badnell) = {1:.6e}'.format(reactor.kf[1], kf_rec))
    check('source evaluated at the Voronov rate',
          abs(reactor.kf[0] - kf_ion) <= 1e-9 * abs(kf_ion))
    check('sink evaluated at the Badnell rate',
          abs(reactor.kf[1] - kf_rec) <= 1e-9 * abs(kf_rec))
    check('the cation has a non-zero loss channel', reactor.kf[1] > 0.0)

    print('\n' + '=' * 78)
    if FAILURES:
        print('LITHIUM CHARGE NETWORK: {0} STAGE CHECK(S) FAILED'.format(len(FAILURES)))
        for f in FAILURES:
            print('  - {0}'.format(f))
        return 1
    print('LITHIUM CHARGE NETWORK: ALL STAGES PASSED -- neutral feed -> cation -> sink')
    print('=' * 78)
    return 0


if __name__ == '__main__':
    sys.exit(main())
