#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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

from pathlib import Path
from types import SimpleNamespace

import pytest

import rmgpy.data.rmg as rmg_data_module
from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary
from rmgpy.exceptions import PlasmaStateError
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.species import Species

FIXTURES = Path(__file__).resolve().parents[1] / 'test_data' / 'excited_states'
ELECTRON = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
N2 = '1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}'


@pytest.fixture
def database(monkeypatch):
    database = ThermoDatabase()
    library = ThermoLibrary(label='ExcitedStateFixture')
    library.load(str(FIXTURES / 'thermo.py'), database.local_context, {})
    database.libraries = {library.label: library}
    database.library_order = [library.label]
    monkeypatch.setattr(rmg_data_module, 'database', SimpleNamespace(thermo=database, solvation=None))
    return database


def reactor_for(species, **kwargs):
    return PlasmaReactor(T=300, P=1e5, initial_mole_fractions={species: 1, ELECTRON: 1e-12}, Te=(10000,'K'), termination=[], **kwargs)


@pytest.mark.parametrize('state,label', [('vibrationallevel 1\n','N2v1'), ('electronicstate A3Su+\n','N2A')])
def test_plasma_accepts_exact_resolved_library_thermo(database, state, label):
    species = Species(label=label).from_adjacency_list(state + N2)
    species.thermo = database.get_thermo_data(species)
    reactor = reactor_for(species)
    reactor.initialize_model([species, ELECTRON], [], [], [])
    assert 'ExcitedStateFixture' in reactor.thermo_provenance_diagnostics[label]
    assert label in reactor.thermo_provenance_diagnostics[label]


@pytest.mark.parametrize('state,label', [('vibrationallevel 1\n','N2v1'), ('electronicstate A3Su+\n','N2A')])
def test_plasma_refuses_altered_resolved_thermo(database, state, label):
    species = Species(label=label).from_adjacency_list(state + N2)
    species.thermo = database.get_thermo_data(species)
    species.thermo.change_base_enthalpy(1000)
    with pytest.raises(PlasmaStateError, match=label):
        reactor_for(species).initialize_model([species, ELECTRON], [], [], [])


def test_plasma_ground_library_entry_cannot_source_excited_state(database):
    species = Species(label='N2v1').from_adjacency_list('vibrationallevel 1\n' + N2)
    species.thermo = database.get_thermo_data(Species(label='N2').from_adjacency_list(N2))
    with pytest.raises(PlasmaStateError, match='N2v1'):
        reactor_for(species).initialize_model([species, ELECTRON], [], [], [])


def test_plasma_refuses_missing_resolved_thermo(database):
    species = Species(label='N2A').from_adjacency_list('electronicstate A3Su+\n' + N2)
    with pytest.raises(PlasmaStateError, match='N2A'):
        reactor_for(species).initialize_model([species, ELECTRON], [], [], [])


def test_plasma_refuses_resolved_state_without_loaded_library_even_with_assertion(monkeypatch):
    monkeypatch.setattr(rmg_data_module, 'database', None)
    species = Species(label='N2A').from_adjacency_list('electronicstate A3Su+\n' + N2)
    from rmgpy.thermo import NASA, NASAPolynomial
    species.thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5,0,0,0,0,70000,3], Tmin=(200,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'))
    with pytest.raises(PlasmaStateError, match='N2A'):
        reactor_for(species, thermo_source_assertions=['N2A']).initialize_model([species, ELECTRON], [], [], [])
