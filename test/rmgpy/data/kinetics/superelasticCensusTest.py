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

"""Keep source discovery and the reviewed reversal census in lockstep."""
import ast
import importlib.util
import json
from pathlib import Path
import re
import shutil

import pytest

ROOT = Path(__file__).resolve().parents[4]
spec = importlib.util.spec_from_file_location(
    'reversal_census_generator', ROOT / 'scripts/generate_reversal_census.py')
generator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generator)
VERDICTS = json.loads(
    (ROOT / 'documentation/source/users/rmg/reversal_census_verdicts.json').read_text())


def test_generated_census_matches_reviewed_sites():
    sites = generator.discover(ROOT)
    generator.validate(sites, VERDICTS)
    checked_in = json.loads(
        (ROOT / 'documentation/source/users/rmg/reversal_census.json').read_text())
    assert checked_in == [dict(row, **VERDICTS[row['site']]) for row in sites]


def test_arkane_reaction_emitters_use_serialization_chokepoint():
    inventory = json.loads(
        (ROOT / 'documentation/source/users/rmg/arkane_writer_inventory.json').read_text())
    emitters = generator.discover_arkane_emitters(ROOT)
    generator.validate_arkane_emitters(emitters, inventory)


def test_arkane_emitter_gate_rejects_chokepoint_bypass(tmp_path):
    inventory = json.loads(
        (ROOT / 'documentation/source/users/rmg/arkane_writer_inventory.json').read_text())
    shutil.copytree(ROOT / 'arkane', tmp_path / 'arkane')
    path = tmp_path / 'arkane/pdep.py'
    source = path.read_text()
    tree = ast.parse(source)
    target = next(
        function
        for cls in tree.body if isinstance(cls, ast.ClassDef) and cls.name == 'PressureDependenceJob'
        for function in cls.body
        if isinstance(function, ast.FunctionDef) and function.name == 'save_input_file'
    )
    lines = source.splitlines(keepends=True)
    body = ''.join(lines[target.lineno - 1:target.end_lineno])
    checked = 'check_reaction_rate_serialization(reaction, serialized_kinetics, True)'
    mutated = body.replace(checked, 'bypassed_serialization_check(reaction, serialized_kinetics, True)')
    assert mutated != body
    lines[target.lineno - 1:target.end_lineno] = [mutated]
    path.write_text(''.join(lines))

    emitters = generator.discover_arkane_emitters(tmp_path)
    with pytest.raises(ValueError, match=(
            'Emitter checks the wrong rate or direction at chokepoint: '
            'arkane/pdep.py::PressureDependenceJob.save_input_file')):
        generator.validate_arkane_emitters(emitters, inventory)


def test_static_gate_rejects_removed_verdict():
    verdicts = dict(VERDICTS)
    verdicts.pop(next(iter(verdicts)))
    with pytest.raises(ValueError, match='Unclassified sites:'):
        generator.validate(generator.discover(ROOT), verdicts)


def test_static_gate_rejects_new_unclassified_site(tmp_path):
    (tmp_path / 'rmgpy').mkdir()
    (tmp_path / 'arkane').mkdir()
    (tmp_path / 'rmgpy/new_exporter.py').write_text(
        'def export(reaction):\n    return reaction.generate_reverse_rate_coefficient()\n')
    discovered = generator.discover(tmp_path)
    assert any(row['kind'] == 'reverse' for row in discovered)
    with pytest.raises(ValueError, match='Unclassified sites:.*new_exporter'):
        generator.validate(generator.discover(ROOT) + discovered, VERDICTS)


def test_static_gate_rejects_disappeared_site():
    sites = generator.discover(ROOT)
    with pytest.raises(ValueError, match='Disappeared sites:'):
        generator.validate(sites[1:], VERDICTS)


@pytest.mark.parametrize('site', sorted(VERDICTS))
def test_census_verdict_has_executable_witness(site):
    verdict = VERDICTS[site]
    for selector in verdict['tests']:
        filename, function = selector.split('::')
        source = (ROOT / filename).read_text()
        assert re.search(r'^\s*def ' + re.escape(function) + r'\s*\(', source, re.M), selector
    assert verdict['reason'].strip()


def test_scanner_covers_cython_templates_and_io(tmp_path):
    (tmp_path / 'rmgpy/tools').mkdir(parents=True)
    (tmp_path / 'arkane').mkdir()
    (tmp_path / 'rmgpy/tools/example.pyx').write_text(
        'cpdef import_reaction(reaction):\n'
        '    # reaction.generate_reverse_rate_coefficient() is prose\n'
        '    return reaction.is_same_reaction(other)\n')
    (tmp_path / 'arkane/render.py').write_text(
        'def write_output(reaction):\n'
        '    template = "{{ reaction.generate_reverse_rate_coefficient() }}"\n'
        '    return template\n')
    rows = generator.discover(tmp_path)
    assert any('example.pyx::import_reaction::is_same_reaction' in row['site'] for row in rows)
    assert any(row['kind'] == 'template-expression' for row in rows)
    assert any(row['site'].endswith('write_output::reaction-io') for row in rows)
    assert sum('generate_reverse_rate_coefficient' in row['site'] for row in rows) == 1


def test_scanner_covers_pressure_dependence_rate_primitives(tmp_path):
    (tmp_path / 'rmgpy/pdep').mkdir(parents=True)
    (tmp_path / 'arkane').mkdir()
    (tmp_path / 'rmgpy/pdep/routes.py').write_text(
        'def routes(reaction, network):\n'
        '    reaction.calculate_microcanonical_rate_coefficient()\n'
        '    calculate_microcanonical_rate_coefficient(reaction)\n'
        '    network.calculate_rate_coefficients()\n'
        '    get_rate_coefficients_CSE_Advanced(network)\n'
        '    get_rate_coefficients_SLS(network)\n')
    sites = {row['site'] for row in generator.discover(tmp_path)}
    for primitive in (
            'calculate_microcanonical_rate_coefficient',
            'calculate_rate_coefficients',
            'get_rate_coefficients_CSE_Advanced',
            'get_rate_coefficients_SLS'):
        assert any(f'::{primitive}#' in site for site in sites), primitive


def test_structural_identity_has_no_reaction_direction():
    from rmgpy.species import Species
    from rmgpy.molecule import Molecule, Group
    resolved = Species(label='N2v1').from_adjacency_list(
        'vibrationallevel 1\n1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}')
    ground = Species(label='N2').from_smiles('N#N')
    assert not resolved.is_isomorphic(ground)
    assert resolved.is_isomorphic(resolved.copy(deep=True))
    for structural_object in (resolved, resolved.molecule[0], Group().from_adjacency_list('1 N u0')):
        assert not hasattr(structural_object, 'generate_reverse_rate_coefficient')
        assert not hasattr(structural_object, 'reversible')


def test_resolved_species_cannot_enter_ordinary_family_recipe():
    from rmgpy.species import Species
    from rmgpy.molecule import Group
    from rmgpy.data.base import Entry
    from rmgpy.data.kinetics.family import KineticsFamily
    species = Species(label='N2v1').from_adjacency_list(
        'vibrationallevel 1\n1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}')
    root = Group().from_adjacency_list('1 *1 N u0')
    assert not KineticsFamily(label='test')._match_reactant_to_template(
        species.molecule[0], Entry(item=root))


def test_scanner_discovers_an_exporter_without_known_name(tmp_path):
    (tmp_path / 'rmgpy').mkdir()
    (tmp_path / 'arkane').mkdir()
    (tmp_path / 'rmgpy/new_writer.py').write_text(
        'def custom(r):\n'
        '    return {"reversible": r.reversible, "kinetics": r.kinetics}\n')
    rows = generator.discover(tmp_path)
    assert len(rows) == 1
    assert rows[0]['site'] == 'rmgpy/new_writer.py::custom::reaction-io'
    with pytest.raises(ValueError, match='Unclassified sites:'):
        generator.validate(rows, {})


def test_nonreaction_serialized_data_has_no_direction():
    from rmgpy.quantity import ScalarQuantity, ArrayQuantity
    from rmgpy.thermo import NASAPolynomial
    for obj in (ScalarQuantity(300, 'K'), ArrayQuantity([300, 1000], 'K'),
                NASAPolynomial(coeffs=[3.5, 0, 0, 0, 0, 0, 1], Tmin=(200, 'K'), Tmax=(3000, 'K'))):
        assert not hasattr(obj, 'generate_reverse_rate_coefficient')
        assert 'reversible' not in obj.as_dict()


def test_generic_dictionary_expansion_retains_reaction_objects():
    from rmgpy.reaction import Reaction
    from rmgpy.rmgobject import RMGObject, expand_to_dict
    reaction = Reaction(reversible=False)
    assert not isinstance(reaction, RMGObject)
    assert not hasattr(reaction, 'as_dict')
    assert not hasattr(reaction, 'make_object')
    assert expand_to_dict(reaction) is reaction
    assert expand_to_dict({'reaction': reaction})['reaction'] is reaction


def test_input_deck_io_does_not_serialize_reaction_rates(tmp_path, monkeypatch):
    from rmgpy.rmg.main import RMG
    from rmgpy.rmg.input import read_input_file, read_thermo_input_file, save_input_file
    from rmgpy.pdep.collision import SingleExponentialDown
    from rmgpy.thermo import ThermoData

    def supply_state_matched_data(model, species, rename=False):
        species.thermo = ThermoData(
            Tdata=([300, 1000], 'K'),
            Cpdata=([30, 30], 'J/(mol*K)'),
            H298=(0, 'kJ/mol'),
            S298=(200, 'J/(mol*K)'))
        species.energy_transfer_model = SingleExponentialDown(alpha0=(1, 'kJ/mol'))

    monkeypatch.setattr('rmgpy.rmg.model.CoreEdgeReactionModel.generate_thermo',
                        supply_state_matched_data)
    species = (
        "species(label='N2v1', reactive=True, structure=adjacencyList("
        "'vibrationallevel 1\\n1 N u0 p1 c0 {2,T}\\n2 N u0 p1 c0 {1,T}'))\n"
        "species(label='N2', reactive=True, structure=adjacencyList("
        "'1 N u0 p1 c0 {2,T}\\n2 N u0 p1 c0 {1,T}'))\n"
        "vibrationalManifold(species='N2')\n"
    )
    database = "database(thermoLibraries=[], reactionLibraries=[], kineticsFamilies=[])\n"
    path = tmp_path / 'input.py'
    path.write_text(database + species +
                    "simpleReactor(temperature=(300, 'K'), pressure=(1, 'bar'), "
                    "initialMoleFractions={'N2v1': 1}, terminationTime=(1e-6, 's'))\n"
                    "simulator(atol=1e-16, rtol=1e-8)\n"
                    "model(toleranceMoveToCore=0.1, toleranceKeepInEdge=0)\n")
    job = RMG()
    read_input_file(str(path), job)
    assert job.initial_species[0].molecule[0].has_resolved_state()
    assert job.reaction_model.core.reactions == []
    saved = tmp_path / 'saved.py'
    from rmgpy.exceptions import SpeciesIdentityError
    with pytest.raises(SpeciesIdentityError, match='save_input_file.*N2v1'):
        save_input_file(str(saved), job)
    assert not saved.exists()
    thermo_path = tmp_path / 'thermo.py'
    thermo_path.write_text(database + species)
    thermo_job = RMG()
    from rmgpy.exceptions import VibrationalManifoldError
    with pytest.raises(VibrationalManifoldError, match='N2v1'):
        read_thermo_input_file(str(thermo_path), thermo_job)
