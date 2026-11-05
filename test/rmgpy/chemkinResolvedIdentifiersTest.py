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

"""Compact resolved identifiers retain exact molecular identity."""

import csv
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from rmgpy.chemkin import (get_species_identifier, load_chemkin_file,
                          save_chemkin_file, save_species_dictionary,
                          validate_species_identifiers)
from rmgpy.exceptions import ChemkinError, ChemkinIdentifierCollisionError
from rmgpy.export import SpeciesReferences
from rmgpy.kinetics import Arrhenius
from rmgpy.reaction import Reaction
from rmgpy.rmg.listener import SimulationProfileWriter
from rmgpy.species import Species

ROOT = Path(__file__).resolve().parents[2]
DATA = json.loads((ROOT / 'test/rmgpy/test_data/chemkin/resolved_identifiers.json').read_text())
spec = importlib.util.spec_from_file_location(
    'compact_identifier_helpers', ROOT / 'test/rmgpy/yaml_writer/excitedExportHelpers.py')
helpers = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helpers)


def census_species(row, index=9999):
    species = Species(label=row['labels'][0]).from_adjacency_list(row['adjacency'])
    species.index = index
    return species


@pytest.mark.parametrize('index', [-1, 0, 99, 9999])
def test_all_pinned_database_states_fit_without_collisions(index):
    assert DATA['database_pin'] == '5f566a3a93fe06a79a89a0be69525e25ec26331d'
    assert len(DATA['rows']) == 53
    species = [census_species(row, index) for row in DATA['rows']]
    identifiers = [get_species_identifier(spc) for spc in species]
    assert max(map(len, identifiers)) <= 16
    assert len(set(identifiers)) == 53
    validate_species_identifiers(species)


@pytest.mark.parametrize('state,tag', [
    ('3Sigma_u+', 't3Sup'), ('1Sigma_u-', 't1Sum'),
    ('2Sigma_u+', 't2Sup'), ('3Pi_g', 't3Pg'),
    ('3Pi_u', 't3Pu'), ('1Pi_g', 't1Pg'), ('1Delta_u', 't1Du'),
    ('3Sup', 'e3Sup'), ('3S_u+', 'e3S__u_p'),
    ('Sigma', 'eSigma'), ('3Sigma_u', 't3Su'),
])
def test_term_encoding_and_literal_namespace(state, tag):
    spc = helpers.nitrogen()
    spc.molecule[0].electronic_state = state
    spc.index = 7
    assert get_species_identifier(spc) == 'N2(' + tag + ')(7)'


def test_declaration_boundary_refuses_constructed_collision():
    first, second = helpers.nitrogen(level=1), helpers.nitrogen(level=2)
    with pytest.raises(ChemkinIdentifierCollisionError, match='collid'):
        SpeciesReferences([first, second], identifiers=['N2(v1)', 'N2(v1)'],
                          collision_type=ChemkinIdentifierCollisionError, context='Chemkin')


def test_term_and_literal_states_are_distinct_in_declaration_set():
    first, second = helpers.nitrogen(), helpers.nitrogen()
    first.molecule[0].electronic_state = '3Sigma_u+'
    second.molecule[0].electronic_state = '3Sup'
    first.index = second.index = 99
    validate_species_identifiers([first, second])
    assert get_species_identifier(first) != get_species_identifier(second)


@pytest.mark.parametrize('electronic,level', [('3Sigma_u+', 12), ('VeryLongStateToken', -1)])
def test_overlong_states_still_refuse(electronic, level):
    spc = helpers.nitrogen(level=level)
    spc.molecule[0].electronic_state = electronic
    spc.index = 9999
    with pytest.raises(ChemkinError, match='state cannot be dropped or truncated'):
        get_species_identifier(spc)


def test_chemkin_mechanism_and_dictionary_round_trip(tmp_path):
    selected = [('N2', '3Sigma_u+', -1), ('Ar', '3P2', -1), ('N2', '', 1)]
    species = [census_species(next(row for row in DATA['rows']
               if (row['formula'], row['electronic'], row['vibrational']) == key), i)
               for i, key in enumerate(selected, 9997)]
    for spc in species:
        spc.thermo = helpers.nitrogen().thermo
    reaction = Reaction(reactants=[species[0], species[1]], products=[species[2], species[1]],
                        reversible=False, kinetics=Arrhenius(A=(1e6, 'm^3/(mol*s)'),
                                                            n=0, Ea=(0, 'J/mol')))
    chemkin = str(tmp_path / 'chem.inp')
    dictionary = str(tmp_path / 'species_dictionary.txt')
    with helpers.unit_thermo_library(species):
        save_chemkin_file(chemkin, species, [reaction])
        save_species_dictionary(dictionary, species)
    loaded, reactions = load_chemkin_file(chemkin, dictionary, use_chemkin_names=True)
    assert len(loaded) == 3
    assert len(reactions) == 1
    assert [spc.index for spc in loaded] == [9997, 9998, 9999]
    assert [(spc.molecule[0].get_formula(), spc.molecule[0].electronic_state,
             spc.molecule[0].vibrational_level) for spc in loaded] == selected
    assert all(first.is_isomorphic(second) for first, second in zip(species, loaded))
    assert reactions[0].is_isomorphic(reaction)
    assert reactions[0].reactants[0] is loaded[0]
    assert reactions[0].products[0] is loaded[2]
    assert 'electronicstate 3Sigma_u+' in Path(dictionary).read_text()
    assert 'N2(t3Sup)(9997)' in Path(chemkin).read_text()


def test_native_profile_header_contains_nitrogen_electronic_state(tmp_path):
    row = next(row for row in DATA['rows'] if row['electronic'] == '3Sigma_u+')
    spc = census_species(row)
    (tmp_path / 'solver').mkdir()
    SimulationProfileWriter(str(tmp_path), 0, [spc]).update(
        SimpleNamespace(snapshots=[[0.0, 1.0, 1.0]]))
    with next((tmp_path / 'solver').glob('*.csv')).open() as stream:
        header = next(csv.reader(stream))
    assert header == ['Time (s)', 'Volume (m^3)', 'N2(t3Sup)(9999)']


def test_chemkin_declarations_refuse_actual_formula_state_index_collision():
    # Two nonisomorphic graphs can share formula, state token, and index.
    first = helpers.nitrogen()
    first.molecule[0].electronic_state = '3Pi_g'
    first.index = 0
    row = next(row for row in DATA['rows'] if row['electronic'] == '3Pi_g')
    second = census_species(row, 0)
    assert not first.is_isomorphic(second)
    assert get_species_identifier(first) == get_species_identifier(second)
    with pytest.raises(ChemkinIdentifierCollisionError, match='collid'):
        validate_species_identifiers([first, second])
