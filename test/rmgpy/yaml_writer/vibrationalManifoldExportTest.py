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


import importlib.util
import json
from pathlib import Path

import pytest
from rmgpy import chemkin, yaml_rms, yaml_cantera1, yaml_cantera2
from rmgpy.exceptions import SpeciesIdentityError
from rmgpy.export import SpeciesReferences, export_molecule, resolve_species_reference
from rmgpy.molecule import Molecule
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species

_spec = importlib.util.spec_from_file_location('manifold_source_fixture',
    Path(__file__).resolve().parents[1] / 'data/excitedThermoSourceTest.py')
_source = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_source)
database = _source.database

_spec = importlib.util.spec_from_file_location('manifold_export_census',
    Path(__file__).with_name('excitedExportCensusTest.py'))
_census = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_census)
# These two formats explicitly carry the header-free declaration while still
# refusing explicitly resolved input. Their existing round-trip tests cover it.
REFUSED = [key for key in _census.N_SITES if key not in {
    'arkane/output.py:save_thermo_lib', 'rmgpy/rmg/input.py:save_input_file'}]


@pytest.fixture
def manifold(database):
    species = _source.nitrogen()
    CoreEdgeReactionModel().declare_vibrational_manifold(species)
    species.get_thermo_data()
    return species


@pytest.mark.parametrize('route', ['cantera1', 'cantera2', 'rms', 'native', 'adjacency', 'repr'])
def test_manifold_export_carries_v0_without_changing_source(manifold, tmp_path, route):
    original = dict(manifold.props), dict(manifold.molecule[0].props)
    if route == 'cantera1':
        record = yaml_cantera1.species_to_dict(manifold, [manifold])
        adjacency = record['note']
    elif route == 'cantera2':
        record = yaml_cantera2.species_to_dict(manifold, [manifold])
        adjacency = record['note']
    elif route == 'rms':
        path = tmp_path / 'manifold.rms'
        yaml_rms.write_rms([manifold], [], path=str(path))
        loaded = yaml_rms.load_rms_species(str(path))[0]
        adjacency = loaded.to_adjacency_list()
    elif route == 'native':
        native = manifold.to_cantera(use_chemkin_identifier=True, all_species=[manifold])
        assert '(v0)' in native.name
        adjacency = native.input_data['note']
    elif route == 'adjacency':
        adjacency = manifold.to_adjacency_list()
    else:
        molecule = eval(repr(manifold.molecule[0]), {'Molecule': Molecule})
        adjacency = molecule.to_adjacency_list()
    loaded = Species().from_adjacency_list(adjacency)
    assert loaded.molecule[0].vibrational_level == 0
    assert not loaded.is_isomorphic(_source.nitrogen())
    assert manifold.molecule[0].vibrational_level == -1
    assert (manifold.props, manifold.molecule[0].props) == original


def test_manifold_reference_resolves_v0_and_refuses_ground_collisions(manifold):
    fixed = _source.nitrogen('vibrationallevel 0\n')
    declarations = SpeciesReferences([manifold], identifiers=chemkin.get_species_identifier)
    assert resolve_species_reference(fixed, declarations) == chemkin.get_species_identifier(manifold)
    with pytest.raises(SpeciesIdentityError):
        resolve_species_reference(_source.nitrogen(), declarations)
    with pytest.raises(SpeciesIdentityError):
        SpeciesReferences([manifold, _source.nitrogen()])
    with pytest.raises(SpeciesIdentityError):
        resolve_species_reference(manifold, [], allow_missing_efficiency=True)


def test_manifold_legacy_dictionary_and_unqualified_native_refuse(manifold):
    with pytest.raises(SpeciesIdentityError):
        chemkin.render_species_dictionary([manifold], old_style=True)
    with pytest.raises(SpeciesIdentityError):
        manifold.to_cantera()


def test_species_owner_declaration_is_carried_without_molecule_marker(manifold):
    manifold.molecule[0].props.pop('vibrational_manifold')
    assert export_molecule(manifold).vibrational_level == 0
    assert 'vibrationallevel 0' in manifold.to_adjacency_list()


@pytest.mark.parametrize('key', REFUSED, ids=REFUSED)
def test_manifold_format_refuses_before_output(manifold, key, tmp_path):
    _census.refusal_witness(key, tmp_path, spc=manifold)


def test_manifold_qm_geometry_refuses_before_embedding_or_writing(manifold, tmp_path):
    from types import SimpleNamespace
    from rmgpy.qm.molecule import Geometry
    with pytest.raises(SpeciesIdentityError, match='Geometry.rd_embed'):
        Geometry.rd_embed(SimpleNamespace(molecule=manifold.molecule[0]), None, 1)
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize('kind', ['solute', 'solvent', 'mixture'])
def test_manifold_solvation_records_carry_v0(manifold, kind, tmp_path):
    import io
    from rmgpy.data.base import Entry
    from rmgpy.data.solvation import save_entry, SoluteLibrary, SolventLibrary, SoluteData, SolventData
    library = SoluteLibrary() if kind == 'solute' else SolventLibrary()
    item = manifold if kind == 'solute' else [manifold] if kind == 'solvent' else [manifold, _source.nitrogen(label='other')]
    data = SoluteData() if kind == 'solute' else SolventData()
    stream = io.StringIO()
    save_entry(stream, Entry(index=1, label='manifold', item=item, data=data))
    assert 'vibrationallevel 0' in stream.getvalue()
    path = tmp_path / 'solvation.py'
    path.write_text(stream.getvalue())
    library.load(str(path))
    item = library.entries['manifold'].item
    restored = item if isinstance(item, Species) else item[0]
    assert restored.molecule[0].vibrational_level == 0
    assert manifold.molecule[0].vibrational_level == -1


def test_manifold_smiles_keyed_efficiency_refuses_before_serialization(manifold):
    from rmgpy.data.kinetics.common import library_serializable_kinetics
    from rmgpy.kinetics import ThirdBody, Arrhenius
    rate = ThirdBody(arrheniusLow=Arrhenius(A=(1., 'm^3/(mol*s)')),
                     efficiencies={manifold.molecule[0]: 2.})
    with pytest.raises(SpeciesIdentityError, match='SMILES-keyed'):
        library_serializable_kinetics(rate, [manifold])


@pytest.mark.parametrize('record_type', ['kinetics', 'forbidden', 'thermo', 'statmech', 'transport'])
def test_manifold_library_dictionary_families_carry_v0(manifold, record_type, tmp_path):
    import io
    from rmgpy.data.base import Entry, ForbiddenStructures
    from rmgpy.data.kinetics.library import KineticsLibrary
    from rmgpy.kinetics import Arrhenius
    from rmgpy.reaction import Reaction
    if record_type == 'kinetics':
        library = KineticsLibrary()
        reaction = Reaction(reactants=[manifold], products=[_source.nitrogen(label='ensemble')])
        library.entries[1] = Entry(index=1, label='N2 <=> ensemble', item=reaction,
                                  data=Arrhenius(A=(1., 's^-1')))
        path = tmp_path / 'reactions.py'
        library.save(str(path))
        restored = KineticsLibrary()
        restored.load(str(path), local_context={'Arrhenius': Arrhenius})
        owner = next(iter(restored.entries.values())).item.reactants[0]
        assert owner.molecule[0].vibrational_level == 0
    elif record_type == 'forbidden':
        library = ForbiddenStructures()
        library.entries['N2'] = Entry(index=1, label='N2', item=manifold)
        path = tmp_path / 'forbidden.py'
        library.save(str(path))
        restored = ForbiddenStructures().load(str(path))
        assert restored.entries['N2'].item.molecule[0].vibrational_level == 0
    else:
        module = __import__('rmgpy.data.' + record_type, fromlist=['save_entry'])
        stream = io.StringIO()
        module.save_entry(stream, Entry(index=1, label='N2', item=manifold.molecule[0]))
        assert 'vibrationallevel 0' in stream.getvalue()
    assert manifold.molecule[0].vibrational_level == -1
