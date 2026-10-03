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

"""Excited-state reverse-rate enforcement through library and reactor APIs."""

from pathlib import Path

import numpy as np
import pytest

from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.exceptions import NonEquilibriumReverseRateError
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.species import Species
from rmgpy.thermo import ThermoData

FIXTURE_ROOT = Path(__file__).resolve().parents[2] / 'test_data' / 'superelastic'


def load_reactions(label):
    database = KineticsDatabase()
    database.load_libraries(str(FIXTURE_ROOT), libraries=[label])
    return database.libraries[label].get_library_reactions()


def reactor_for(reactions):
    species = {sp.label: sp for rxn in reactions for sp in rxn.reactants + rxn.products}
    ion = Species(label='Ar+').from_adjacency_list('1 Ar u1 p3 c+1')
    species[ion.label] = ion
    for sp in species.values():
        if not sp.is_electron():
            sp.thermo = ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
                Cpdata=([29.1] * 7, 'J/(mol*K)'),
                H298=(27.87 if sp.label == 'N2v1' else 0.0, 'kJ/mol'),
                S298=(191.6, 'J/(mol*K)'),
            )
    if 'e-' not in species:
        species['e-'] = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    fractions = {sp: (1.0e-4 if sp.label in ('e-', 'Ar+') else
                     0.1 if sp.label == 'N2v1' else 0.8998) for sp in species.values()}
    reactor = PlasmaReactor(300.0, 1.0e5, fractions, (11604.5, 'K'),
                            termination=[], thermo_source_assertions=['Ar+'])
    return reactor, list(species.values())


def test_probe_reversible_excited_collision():
    """Record the earliest refusal while attempting a real library-to-reactor run."""
    try:
        reactions = load_reactions('reversible')
    except NonEquilibriumReverseRateError as error:
        print('PROBE: refused at library load:', error)
        assert 'N2v1' in str(error)
        return
    print('PROBE: library load accepted')
    reactor, species = reactor_for(reactions)
    stage = 'reactor initialization'
    with pytest.raises(NonEquilibriumReverseRateError) as caught:
        reactor.initialize_model(species, reactions, [], [])
        stage = 'reactor advance'
        reactor.advance(1.0e-9)
    print('PROBE: refused at', stage + ':', caught.value)
    assert 'N2v1' in str(caught.value)


def test_reversible_excited_collision_refused_at_library_load():
    with pytest.raises(NonEquilibriumReverseRateError) as caught:
        load_reactions('reversible')
    message = str(caught.value)
    assert 'e- + N2v1 <=> e- + N2' in message
    assert 'N2v1' in message
    assert 'vibrationallevel 1' in message
    assert 'irreversible' in message


def test_explicit_irreversible_pair_loads_and_integrates():
    reactions = load_reactions('explicit-pair')
    assert len(reactions) == 2
    assert all(not rxn.reversible for rxn in reactions)
    assert reactions[0].is_isomorphic(reactions[1], either_direction=True)
    assert not reactions[0].is_isomorphic(reactions[1], either_direction=False)
    reactor, species = reactor_for(reactions)
    reactor.initialize_model(species, reactions, [], [])
    assert list(reactor.kb) == [0.0, 0.0]
    assert list(reactor.kf) == pytest.approx([1.0e3, 1.0e2])
    excited_index = next(i for i, sp in enumerate(species) if sp.label == 'N2v1')
    before = reactor.y[excited_index]
    reactor.advance(1.0e-6)
    assert np.all(np.isfinite(reactor.y))
    assert 0.0 < reactor.y[excited_index] < before


def test_thermal_resolved_reaction_keeps_equilibrium_reverse():
    reactions = load_reactions('thermal')
    rxn = reactions[0]
    assert rxn.reversible
    assert rxn.get_reverse_from_equilibrium_refusal() is None
    reactor, species = reactor_for(reactions)
    reactor.initialize_model(species, reactions, [], [])
    # Independently known level spacing: equal Cp/S, DeltaH=-27.87 kJ/mol,
    # and no change in particle count. Kc=exp(27870/(R*300)).
    # Allow the small difference between modern R and RMG's older CODATA value.
    expected_keq = np.exp(27870.0 / (8.31446261815324 * 300.0))
    assert float(reactor.Keq[0]) == pytest.approx(expected_keq, rel=2.0e-5)
    assert float(reactor.kb[0]) == pytest.approx(1.0e3 / expected_keq, rel=2.0e-5)
    reverse = rxn.generate_reverse_rate_coefficient(Tmin=rxn.kinetics.Tmin, Tmax=rxn.kinetics.Tmax)
    assert reverse.get_rate_coefficient(300.0) == pytest.approx(1.0e3 / expected_keq, rel=2.0e-5)


def test_unresolved_te_dependent_entry_keeps_baseline_behavior():
    reactions = load_reactions('unresolved')
    rxn = reactions[0]
    assert rxn.reversible
    assert not any(mol.has_resolved_state() for sp in rxn.reactants + rxn.products for mol in sp.molecule)
    with pytest.raises(NonEquilibriumReverseRateError):
        rxn.check_reverse_from_equilibrium_supported()
    reactor, species = reactor_for(reactions)
    with pytest.raises(NonEquilibriumReverseRateError, match=r'Te-dependent reaction e- \+ N2 <=> e- \+ N \+ N'):
        reactor.initialize_model(species, reactions, [], [])


TE_KINETICS = [
    "TwoTemperaturePlasma(A=(1.0e3, 'm^3/(mol*s)'))",
    "ElectronCollisionPlasma(energies=([0.0, 1.0, 2.0], 'eV/molecule'), sigma=([0.0, 1.0e-20, 1.0e-20], 'm^2'))",
    "BadnellRRArrhenius(A=(1.0e3, 'm^3/(mol*s)'), T0=(1.0, 'K'), T1=(1.0, 'K'), electrons=0)",
    "VoronovEIArrhenius(A=(1.0e3, 'm^3/(mol*s)'), dE=(1.0, 'eV'), electrons=0)",
]


@pytest.mark.parametrize('kinetics', TE_KINETICS, ids=['two-temperature', 'collision', 'badnell', 'voronov'])
@pytest.mark.parametrize('state', ['vibrationallevel 0', 'electronicstate A3Su+',
                                   'electronicstate A3Su+\nvibrationallevel 1'])
@pytest.mark.parametrize('equation', ['e- + N2v1 <=> e- + N2', 'e- + N2 <=> e- + N2v1'])
def test_resolved_state_and_te_kinetics_variants_refused_at_load(tmp_path, kinetics, state, equation):
    """Both sides, both state headers (including level zero), and every current Te rate law."""
    from rmgpy.data.kinetics.library import KineticsLibrary

    template = FIXTURE_ROOT / 'reversible'
    dictionary = (template / 'dictionary.txt').read_text().replace('vibrationallevel 1', state)
    (tmp_path / 'dictionary.txt').write_text(dictionary)
    (tmp_path / 'reactions.py').write_text(
        "name = 'state variants'\nshortDesc = ''\nlongDesc = ''\n"
        "entry(index=1, label={!r}, kinetics={})\n".format(equation, kinetics))
    database = KineticsDatabase()
    with pytest.raises(NonEquilibriumReverseRateError) as caught:
        KineticsLibrary(label='variants').load(str(tmp_path / 'reactions.py'),
                                             database.local_context, database.global_context)
    message = str(caught.value)
    assert equation in message
    assert 'Resolved species: N2v1' in message
    for header in state.splitlines():
        assert header in message


def test_resolved_collision_reactor_refusal_remains():
    """Bypass library loading to retain the reactor's independent safety net."""
    rxn = load_reactions('explicit-pair')[0]
    rxn.reversible = True
    reactor, species = reactor_for([rxn])
    with pytest.raises(NonEquilibriumReverseRateError, match='N2v1'):
        reactor.initialize_model(species, [rxn], [], [])


def test_resolved_collision_reverse_guard_refusal_remains():
    rxn = load_reactions('explicit-pair')[0]
    rxn.reversible = True
    with pytest.raises(NonEquilibriumReverseRateError, match='N2v1'):
        rxn.check_reverse_from_equilibrium_supported()
