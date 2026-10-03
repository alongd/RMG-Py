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

"""State-matched electron records must be independent of entry/library order."""

import pytest

from rmgpy.data.base import Entry
from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary
from rmgpy.exceptions import StateProvenanceError
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.thermo import ThermoData

GROUND = ('', -1)
STATES = [('A', -1), ('', 1), ('A', 1)]


def electron(state, radicals=0):
    molecule = Molecule().from_adjacency_list('1 e u{} p0 c-1'.format(radicals))
    molecule.electronic_state, molecule.vibrational_level = state
    return Species(molecule=[molecule])


def record(label, state, enthalpy):
    data = ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
                      Cpdata=([0] * 7, 'J/(mol*K)'), H298=(enthalpy, 'J/mol'),
                      S298=(0, 'J/(mol*K)'))
    return Entry(label=label, item=electron(state).molecule[0], data=data)


def lookup(database, species, api, library=None):
    if api == 'direct':
        return database.get_thermo_data_from_library(species, library)[0]
    if api == 'libraries':
        return database.get_thermo_data_from_libraries(species)[0]
    return database.get_thermo_data(species)


@pytest.mark.parametrize('record_state', STATES)
@pytest.mark.parametrize('query_state', [GROUND] + STATES)
@pytest.mark.parametrize('order', [('resolved', 'ground'), ('ground', 'resolved')])
@pytest.mark.parametrize('radicals', [0, 1])
@pytest.mark.parametrize('api', ['direct', 'libraries', 'full'])
def test_mixed_electron_entries_match_the_query_in_either_order(record_state, query_state, order, radicals, api):
    entries = {'ground': record('ground', GROUND, 0),
               'resolved': record('resolved', record_state, 12345)}
    library = ThermoLibrary(label='MixedElectron')
    library.entries = {label: entries[label] for label in order}
    database = ThermoDatabase()
    database.library_order = [library.label]
    database.libraries = {library.label: library}
    species = electron(query_state, radicals)
    if query_state not in (GROUND, record_state):
        with pytest.raises(StateProvenanceError, match='electron.*state'):
            lookup(database, species, api, library)
    else:
        result = lookup(database, species, api, library)
        expected = entries['ground' if query_state == GROUND else 'resolved']
        assert result.H298.value_si == expected.data.H298.value_si
        assert result is not expected.data
        assert result.label == expected.label
    assert (species.molecule[0].electronic_state, species.molecule[0].vibrational_level) == query_state


@pytest.mark.parametrize('record_state', STATES)
@pytest.mark.parametrize('query_state', [GROUND] + STATES)
@pytest.mark.parametrize('order', [('ResolvedElectron', 'GroundElectron'), ('GroundElectron', 'ResolvedElectron')])
@pytest.mark.parametrize('api', ['libraries', 'full'])
@pytest.mark.parametrize('solvent', [None, 'water'])
def test_electron_search_continues_to_a_matching_later_library(record_state, query_state, order, api, solvent, monkeypatch):
    import rmgpy.rmg.main
    monkeypatch.setattr(rmgpy.rmg.main, 'solvent', solvent)
    libraries = {}
    for name, state, enthalpy in [('GroundElectron', GROUND, 0), ('ResolvedElectron', record_state, 12345)]:
        library = ThermoLibrary(label=name, solvent=solvent)
        library.entries = {name: record(name, state, enthalpy)}
        libraries[name] = library
    database = ThermoDatabase()
    database.libraries = libraries
    database.library_order = list(order)
    species = electron(query_state, 1)
    if query_state not in (GROUND, record_state):
        with pytest.raises(StateProvenanceError, match='electron.*state'):
            lookup(database, species, api)
    else:
        result = lookup(database, species, api)
        assert result.H298.value_si == (0 if query_state == GROUND else 12345)
        assert result.label == ('GroundElectron' if query_state == GROUND else 'ResolvedElectron')


@pytest.mark.parametrize('state', STATES)
@pytest.mark.parametrize('api', ['direct', 'libraries', 'full'])
def test_missing_resolved_electron_record_refuses_even_an_empty_library(state, api):
    library = ThermoLibrary(label='Empty')
    database = ThermoDatabase()
    database.libraries = {library.label: library}
    database.library_order = [library.label]
    with pytest.raises(StateProvenanceError, match='electron.*state'):
        lookup(database, electron(state), api, library)


@pytest.mark.parametrize('state', STATES)
def test_matched_resolved_thermo_reaches_cp_and_symmetry_arithmetic(state):
    molecule = Molecule(smiles='C', electronic_state=state[0], vibrational_level=state[1])
    species = Species(molecule=[molecule])
    library = ThermoLibrary(label='ResolvedMethane')
    entry = record('matched', GROUND, 12345)
    entry.item = molecule.copy(deep=True)
    library.entries = {entry.label: entry}
    database = ThermoDatabase()
    database.libraries = {library.label: library}
    database.library_order = [library.label]
    result = database.get_thermo_data(species)
    assert result.H298.value_si == 12345
    assert result.Cp0.value_si == molecule.calculate_cp0() == species.calculate_cp0()
    assert result.CpInf.value_si == molecule.calculate_cpinf() == species.calculate_cpinf()
    assert species.get_symmetry_number() == molecule.calculate_symmetry_number() == 12
    assert (molecule.electronic_state, molecule.vibrational_level) == state
