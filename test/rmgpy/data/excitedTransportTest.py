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

import yaml

import pytest

from rmgpy.chemkin import save_chemkin_file, save_transport_file
from rmgpy.data.transport import TransportDatabase, TransportLibrary
from rmgpy.species import Species
from rmgpy.rmg.model import ReactionModel
from rmgpy.thermo import NASA, NASAPolynomial
from rmgpy.yaml_cantera2 import save_cantera_model, species_to_dict

FIXTURES = Path(__file__).resolve().parents[1] / 'test_data' / 'excited_states'
N2 = '1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}'
BORROWED = 'Transport borrowed from the ground state for vibrationally resolved species'
FALLBACK = 'Transport fallback to the ground state for electronically resolved species'


def nitrogen(state='', label='N2'):
    return Species(label=label).from_adjacency_list(state + N2)


def transport_database(electronic=False):
    database = TransportDatabase()
    for filename in (['transport_ground.py', 'transport_electronic.py'] if electronic else ['transport_ground.py']):
        library = TransportLibrary(label=filename[:-3])
        library.load(str(FIXTURES / filename), database.local_context, {})
        database.libraries[library.label] = library
        database.library_order.append(library.label)
    return database


def assert_ground_numbers(data):
    assert data.epsilon.value_si == pytest.approx(831.4462618, rel=0.00001)
    assert data.sigma.value_si == pytest.approx(3.7e-10)
    assert data.shapeIndex == 1
    assert data.rotrelaxcollnum == 4


def test_vibrational_transport_borrows_ground_without_mutation():
    database = transport_database()
    species = nitrogen('vibrationallevel 1\n', 'N2v1')
    data = database.get_transport_properties(species)[0]
    assert_ground_numbers(data)
    assert BORROWED in data.comment
    assert 'vibrationallevel 1' in data.comment
    assert species.molecule[0].vibrational_level == 1
    ground = database.get_transport_properties(nitrogen())[0]
    assert ground.comment == 'transport_ground'
    assert BORROWED not in database.libraries['transport_ground'].entries['N2'].data.comment


def test_electronic_transport_uses_exact_library_entry(caplog):
    data = transport_database(electronic=True).get_transport_properties(nitrogen('electronicstate A3Su+\n', 'N2A'))[0]
    assert data.epsilon.value_si == pytest.approx(2078.6156545, rel=0.00001)
    assert data.sigma.value_si == pytest.approx(4.1e-10)
    assert data.comment == 'transport_electronic'
    assert not caplog.records


def test_electronic_transport_fallback_warns_and_reaches_exports(caplog, tmp_path, monkeypatch):
    from types import SimpleNamespace
    from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary
    import rmgpy.data.rmg as data_module
    database = ThermoDatabase()
    library = ThermoLibrary(label='ExcitedStateFixture')
    library.load(str(FIXTURES / 'thermo.py'), database.local_context, {})
    database.libraries = {library.label: library}
    database.library_order = [library.label]
    monkeypatch.setattr(data_module, 'database', SimpleNamespace(thermo=database, solvation=None))
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    species.transport_data = transport_database().get_transport_properties(species)[0]
    data = species.transport_data
    assert_ground_numbers(data)
    assert FALLBACK in data.comment
    assert 'electronicstate A3Su+' in data.comment
    assert 'N2A' in caplog.text
    assert FALLBACK in caplog.text
    species.thermo = species.get_thermo_data()
    assert FALLBACK in species_to_dict(species, [species])['transport']['note']
    cantera = tmp_path / 'chem_annotated.yaml'
    save_cantera_model(ReactionModel(species=[species], reactions=[]), str(cantera))
    exported = yaml.safe_load(cantera.read_text())
    assert FALLBACK in exported['species'][0]['transport']['note']
    annotated = tmp_path / 'chem_annotated.inp'
    save_chemkin_file(str(annotated), [species], [], verbose=True)
    assert FALLBACK in annotated.read_text()
    transport = tmp_path / 'tran.dat'
    save_transport_file(str(transport), [species])
    assert FALLBACK in transport.read_text()


@pytest.mark.parametrize('state,comment', [('vibrationallevel 1\n', BORROWED), ('electronicstate A3Su+\n', FALLBACK)])
def test_all_transport_query_obeys_state_rule(state, comment):
    data = transport_database().get_all_transport_properties(nitrogen(state))[0][0]
    assert_ground_numbers(data)
    assert comment in data.comment


def test_vibrational_transport_ignores_a_state_specific_library_entry():
    database = transport_database()
    from rmgpy.data.base import Entry
    from rmgpy.transport import TransportData
    level = nitrogen('vibrationallevel 1\n', 'N2v1')
    library = database.libraries['transport_ground']
    library.entries['N2v1'] = Entry(index=1, label='N2v1', item=level.molecule[0],
                                   data=TransportData(shapeIndex=1, epsilon=(900,'K'), sigma=(8,'angstroms')))
    data = database.get_transport_properties(level)[0]
    assert_ground_numbers(data)
    assert BORROWED in data.comment


@pytest.mark.parametrize('radicals', [0, 1])
@pytest.mark.parametrize('state', ['electronic', 'vibrational'])
@pytest.mark.parametrize('route', ['attached', 'generate', 'properties', 'all', 'library', 'groups', 'lennard_jones', 'critical'])
def test_mutated_canonical_electron_cannot_acquire_transport(state, route, radicals):
    from rmgpy.data.base import Entry
    from rmgpy.exceptions import ExcitedSpeciesThermoError, StateProvenanceError
    from rmgpy.transport import TransportData

    species = Species(label='e-').from_adjacency_list('1 e u{} p0 c-1'.format(radicals))
    if state == 'electronic':
        species.molecule[0].electronic_state = 'TEST'
    else:
        species.molecule[0].vibrational_level = 1
    species.transport_data = TransportData(shapeIndex=0, epsilon=(100, 'K'), sigma=(3, 'angstroms'))
    database = TransportDatabase()
    library = TransportLibrary(label='ElectronMutationProbe')
    library.entries['e-'] = Entry(index=1, label='e-', item=species.molecule[0].copy(deep=True),
                                  data=species.transport_data)
    database.libraries[library.label] = library
    database.library_order = [library.label]
    calls = {
        'attached': species.get_transport_data,
        'generate': species.generate_transport_data,
        'properties': lambda: database.get_transport_properties(species),
        'all': lambda: database.get_all_transport_properties(species),
        'library': lambda: database.get_transport_properties_from_library(species, library),
        'groups': lambda: database.get_transport_properties_via_group_estimates(species),
        'critical': lambda: database.estimate_critical_properties_via_group_additivity(species.molecule[0]),
        'lennard_jones': lambda: database.get_transport_properties_via_lennard_jones_parameters(species),
    }
    expected = StateProvenanceError if route == 'critical' and radicals else ExcitedSpeciesThermoError
    message = 'resolved-state' if expected is StateProvenanceError else 'electron'
    with pytest.raises(expected, match=message):
        calls[route]()
