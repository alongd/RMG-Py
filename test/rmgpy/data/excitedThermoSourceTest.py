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

import rmgpy.data.rmg as data_module
import rmgpy.rmg.input as inp
from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary
from rmgpy.exceptions import NetworkError, SpeciesIdentityError, ExcitedSpeciesThermoError, StateProvenanceError, VibrationalManifoldError
from rmgpy.species import Species

FIXTURES = Path(__file__).resolve().parents[1] / 'test_data' / 'excited_states'
N2 = '1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}'


def nitrogen(state='', label='N2'):
    return Species(label=label).from_adjacency_list(state + N2)


@pytest.fixture
def database(monkeypatch):
    db = ThermoDatabase()
    lib = ThermoLibrary(label='ExactStateFixture')
    lib.load(str(FIXTURES / 'thermo.py'), db.local_context, {})
    db.libraries = {lib.label: lib}
    db.library_order = [lib.label]
    monkeypatch.setattr(data_module, 'database', SimpleNamespace(thermo=db, solvation=None))
    monkeypatch.setattr(inp, 'rmg', None)
    return db


@pytest.mark.parametrize('field, value', [
    ('E0', (1234, 'kJ/mol')),
    ('Cp0', (20000, 'J/(mol*K)')),
    ('CpInf', (10000, 'J/(mol*K)')),
])
def test_attached_energy_and_limits_must_match_whole_library_source(database, field, value):
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    thermo = species.get_thermo_data()
    setattr(thermo, field, value)
    with pytest.raises(ExcitedSpeciesThermoError, match='N2A'):
        species.get_thermo_data()


@pytest.mark.parametrize('state,entry,field', [
    ('electronicstate A3Su+\n', 'N2A', 'CpInf'),
    ('vibrationallevel 1\n', 'N2v1', 'Cp0'),
])
def test_missing_library_limits_refuse_without_formula(database, state, entry, field):
    setattr(database.libraries['ExactStateFixture'].entries[entry].data, field, None)
    with pytest.raises(ExcitedSpeciesThermoError, match=field):
        nitrogen(state, entry).get_thermo_data()


def test_configuration_refreshes_energy_after_state_mutation(database):
    from rmgpy.pdep.configuration import Configuration
    species = nitrogen()
    species.get_thermo_data()
    initial = Configuration(species).E0
    species.molecule[0].electronic_state = 'A3Su+'
    actual = Configuration(species).E0
    assert actual != initial
    assert actual == pytest.approx(582013.0788860451)


@pytest.mark.parametrize('state', ['electronic', 'vibrational', 'declared'])
@pytest.mark.parametrize('consumer', ['cantera', 'configuration'])
def test_consumers_refresh_ground_thermo_after_state_is_cleared(database, state, consumer):
    from rmgpy.pdep.configuration import Configuration
    from rmgpy.rmg.model import CoreEdgeReactionModel
    ground = nitrogen()
    ground.get_thermo_data()
    expected_energy = Configuration(ground).E0
    species = nitrogen()
    if state == 'electronic':
        species.molecule[0].electronic_state = 'A3Su+'
    elif state == 'vibrational':
        species.molecule[0].vibrational_level = 1
    else:
        CoreEdgeReactionModel().declare_vibrational_manifold(species)
    species.get_thermo_data()
    if state == 'electronic':
        species.molecule[0].electronic_state = ''
    elif state == 'vibrational':
        species.molecule[0].vibrational_level = -1
    else:
        species.props.pop('vibrational_manifold')
        for molecule in species.molecule:
            molecule.props.pop('vibrational_manifold', None)
    if consumer == 'cantera':
        exported = species.to_cantera()
        assert exported.thermo.cp(298.15) / 1000 == pytest.approx(37.415081781)
    else:
        assert Configuration(species).E0 == pytest.approx(expected_energy)


def test_cantera_refuses_unmatched_attached_state(database):
    species = nitrogen('electronicstate absent\n', 'absent')
    species.thermo = database.libraries['ExactStateFixture'].entries['N2'].data
    with pytest.raises(ExcitedSpeciesThermoError, match='absent'):
        species.to_cantera()


def test_arkane_thermo_refuses_before_changing_attached_data(database):
    from arkane.thermo import ThermoJob
    from rmgpy.statmech import Conformer
    species = Species(label='O1D').from_adjacency_list('electronicstate 1D\n1 O u0 p3 c0')
    species.conformer = Conformer(E0=(200, 'kJ/mol'), spin_multiplicity=1)
    with pytest.raises(ExcitedSpeciesThermoError, match='O1D'):
        ThermoJob.generate_thermo(SimpleNamespace(species=species, thermo_class='NASA'))
    assert species.thermo is None


def test_isotope_generation_and_direct_entropy_correction_refuse(database):
    from rmgpy.tools.isotopes import generate_isotopomers, correct_entropy
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    thermo = species.get_thermo_data()
    before = thermo.get_entropy(298.15)
    with pytest.raises(ExcitedSpeciesThermoError, match='N2A'):
        generate_isotopomers(species)
    with pytest.raises(ExcitedSpeciesThermoError, match='N2A'):
        correct_entropy(species, nitrogen())
    assert species.thermo.get_entropy(298.15) == before


def test_exact_solute_does_not_authorize_derived_solvation(database, monkeypatch):
    from rmgpy.thermo.thermoengine import process_thermo_data
    # A correction-capable exact-state solute provider still cannot supply liquid thermo.
    solvation = SimpleNamespace(
        get_solute_data=lambda species: SimpleNamespace(comment='exact-state solute'),
        get_solvent_data=lambda name: object(),
        get_solvation_correction=lambda *args: SimpleNamespace(enthalpy=-29754.926, entropy=-59.53))
    data_module.database.solvation = solvation
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    thermo = species.get_thermo_data()
    before = thermo.get_enthalpy(298.15)
    with pytest.raises(ExcitedSpeciesThermoError, match='N2A'):
        process_thermo_data(species, thermo, solvent_name='water')
    assert species.thermo.get_enthalpy(298.15) == before


def test_molecule_aware_ring_helper_refuses_named_error(database):
    species = Species(label='CycloExcited').from_smiles('C1CCCCC1')
    species.molecule[0].electronic_state = 'TEST'
    ring = species.molecule[0].get_smallest_set_of_smallest_rings()[0]
    with pytest.raises(ExcitedSpeciesThermoError):
        database._add_ring_correction_thermo_data_from_tree(None, None, species.molecule[0], ring)


def test_atom_only_ring_helper_retains_resolved_owner(database):
    species = Species(label='BicycloExcited').from_smiles('C1=CC2CCCCC2C=C1')
    species.molecule[0].electronic_state = 'TEST'
    _, poly = species.molecule[0].get_disparate_cycles()
    with pytest.raises(ExcitedSpeciesThermoError):
        database.get_bicyclic_correction_thermo_data_from_heuristic(poly[0])


def test_manifold_library_export_preserves_v0_identity(database, tmp_path):
    from arkane.output import save_thermo_lib
    from rmgpy.rmg.model import CoreEdgeReactionModel
    species = nitrogen()
    model = CoreEdgeReactionModel()
    model.declare_vibrational_manifold(species)
    expected = species.get_thermo_data().get_heat_capacity(298.15)
    ensemble = nitrogen().get_thermo_data().get_heat_capacity(298.15)
    save_thermo_lib([species], str(tmp_path), 'saved', 'test')
    loaded = ThermoLibrary(label='Saved')
    loaded.load(str(tmp_path / 'saved.py'), database.local_context, {})
    assert loaded.entries['N2'].item.vibrational_level == 0
    database.libraries['Saved'] = loaded
    database.library_order.insert(0, 'Saved')
    assert nitrogen().get_thermo_data().get_heat_capacity(298.15) == ensemble
    fixed = nitrogen('vibrationallevel 0\n')
    assert fixed.get_thermo_data().get_heat_capacity(298.15) == expected


def test_thermodata_library_reload_preserves_energy_and_limits(database, tmp_path):
    from rmgpy.thermo import ThermoData
    library = ThermoLibrary(label='Table')
    library.load_entry(index=1, label='fixed', molecule='vibrationallevel 1\n' + N2,
                       thermo=ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
                                        Cpdata=([29.1] * 7, 'J/(mol*K)'), H298=(20, 'kJ/mol'),
                                        S298=(200, 'J/(mol*K)'), Cp0=(29.1, 'J/(mol*K)'),
                                        CpInf=(29.1, 'J/(mol*K)'), E0=(10, 'kJ/mol')))
    library.save(str(tmp_path / 'table.py'))
    loaded = ThermoLibrary(label='Table')
    loaded.load(str(tmp_path / 'table.py'), database.local_context, {})
    original = library.entries['fixed'].data
    actual = loaded.entries['fixed'].data
    for field in ('E0', 'Cp0', 'CpInf'):
        assert getattr(actual, field) is not None
        assert getattr(actual, field).value_si == getattr(original, field).value_si
    database.libraries = {'Table': loaded}
    database.library_order = ['Table']
    assert nitrogen('vibrationallevel 1\n').get_thermo_data() is not None


def admit(model, species):
    phase = model.edge.phase_system.phases['Default']
    phase.names.append(species.label)
    phase.species.append(SimpleNamespace(name=species.label))
    phase.rmg_species.append(species)
    model.edge.phase_system.species_dict[species.label] = species
    model.edge.species.append(species)


def test_ground_rename_keeps_own_phase_name_and_promotion(database):
    from rmgpy.rmg.model import CoreEdgeReactionModel
    model = CoreEdgeReactionModel()
    species, _ = model.make_new_species(nitrogen(), label='N2')
    admit(model, species)
    species.thermo = None
    model.generate_thermo(species, rename=True)
    assert species.label == 'N2'
    assert model.edge.phase_system.phases['Default'].names == ['N2']
    model.edge.phase_system.pass_species(species.label, model.core.phase_system)
    assert model.core.phase_system.species_dict['N2'] is species


def test_admitted_automatic_label_preemption_is_atomic(database):
    from rmgpy.rmg.model import CoreEdgeReactionModel
    from rmgpy.exceptions import DuplicateSpeciesLabelError
    model = CoreEdgeReactionModel()
    excited, _ = model.make_new_species(nitrogen('electronicstate A3Su+\n', '').molecule[0], generate_thermo=False)
    admit(model, excited)
    name = excited.label
    with pytest.raises(DuplicateSpeciesLabelError, match=name.replace('+', '\\+')):
        model.make_new_species(Species().from_smiles('C'), label=name, generate_thermo=False)
    assert excited.label == name
    assert model.edge.phase_system.phases['Default'].names == [name]
    model.edge.phase_system.pass_species(name, model.core.phase_system)
    assert model.core.phase_system.species_dict[name] is excited


def api_call(module, scope, **values):
    """Call a public boundary with only the inputs needed to observe its refusal."""
    import importlib
    import inspect
    function = importlib.import_module(module)
    for name in scope.split('.'):
        function = getattr(function, name)
    signature = inspect.signature(function)
    args = []
    kwargs = {}
    for name, parameter in signature.parameters.items():
        if parameter.kind in (parameter.VAR_POSITIONAL, parameter.VAR_KEYWORD):
            continue
        if name in values:
            value = values[name]
        elif parameter.default is not parameter.empty:
            continue
        else:
            value = None
        if parameter.kind == parameter.POSITIONAL_ONLY:
            args.append(value)
        else:
            kwargs[name] = value
    return function(*args, **kwargs)


@pytest.mark.parametrize('module,scope,owner', [
    ('rmgpy.data.thermo', 'ThermoDatabase._add_adsorption_correction', 'molecule'),
    ('rmgpy.data.thermo', 'ThermoDatabase._add_polycyclic_correction_thermo_data', 'molecule'),
    ('rmgpy.data.thermo', 'ThermoDatabase._add_group_thermo_data', 'molecule'),
    ('rmgpy.data.thermo', 'ThermoDatabase._remove_group_thermo_data', 'molecule'),
    ('rmgpy.data.solvation', 'SolvationDatabase._add_polycyclic_correction_solute_data', 'molecule'),
    ('rmgpy.data.solvation', 'SolvationDatabase._add_ring_correction_solute_data_from_tree', 'molecule'),
    ('rmgpy.data.solvation', 'SolvationDatabase._add_group_solute_data', 'molecule'),
    ('rmgpy.data.solvation', 'SolvationDatabase._remove_group_solute_data', 'molecule'),
    ('rmgpy.data.statmech', 'StatmechGroups.get_statmech_data', 'molecule'),
    ('rmgpy.data.statmech', 'StatmechDatabase.get_statmech_data', 'molecule'),
    ('rmgpy.data.statmech', 'StatmechDatabase.get_statmech_data_from_depository', 'molecule'),
    ('rmgpy.data.statmech', 'StatmechDatabase.get_statmech_data_from_library', 'molecule'),
    ('rmgpy.data.statmech', 'StatmechDatabase.get_statmech_data_from_groups', 'molecule'),
])
def test_lower_acquisition_census_refuses(module, scope, owner, database):
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    upstream_refusals = {
        'ThermoDatabase._add_adsorption_correction',
        'ThermoDatabase._add_polycyclic_correction_thermo_data',
        'SolvationDatabase._add_polycyclic_correction_solute_data',
        'StatmechGroups.get_statmech_data',
    }
    expected = StateProvenanceError if scope in upstream_refusals else ExcitedSpeciesThermoError
    with pytest.raises(expected):
        api_call(module, scope, self=database, **{owner: species.molecule[0]})


@pytest.mark.parametrize('module,scope', [
    ('arkane.thermo', 'ThermoJob.generate_thermo'),
    ('arkane.thermo', 'ThermoJob.plot'),
    ('arkane.statmech', 'StatMechJob.load'),
    ('arkane.statmech', 'StatMechJob.write_output'),
    ('arkane.common', 'ArkaneSpecies.update_species_attributes'),
    ('rmgpy.qm.molecule', 'QMMolecule.load_thermo_data'),
    ('rmgpy.qm.molecule', 'QMMolecule.save_thermo_data'),
    ('rmgpy.qm.molecule', 'QMMolecule.calculate_thermo_data'),
    ('arkane.sensitivity', 'KineticsSensitivity.perturb'),
    ('arkane.sensitivity', 'KineticsSensitivity.unperturb'),
])
def test_jobs_and_mutators_census_refuses_before_use(module, scope, database):
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    species.get_thermo_data()
    job = SimpleNamespace(species=species, molecule=species.molecule[0])
    upstream_refusals = {
        'StatMechJob.load', 'QMMolecule.load_thermo_data', 'QMMolecule.calculate_thermo_data',
    }
    expected = (SpeciesIdentityError if scope == "ThermoJob.plot" else
                StateProvenanceError if scope in upstream_refusals else ExcitedSpeciesThermoError)
    with pytest.raises(expected):
        api_call(module, scope, self=job, species=species)


@pytest.mark.parametrize('module,scope', [
    ('arkane.pdep', 'PressureDependenceJob.execute'),
    ('arkane.pdep', 'PressureDependenceJob.initialize'),
    ('arkane.pdep', 'PressureDependenceJob.save'),
    ('arkane.pdep', 'PressureDependenceJob.save_input_file'),
    ('rmgpy.pdep.network', 'Network.initialize'),
    ('rmgpy.pdep.network', 'Network.calculate_rate_coefficients'),
    ('rmgpy.pdep.network', 'Network.calculate_densities_of_states'),
    ('rmgpy.pdep.network', 'Network.calculate_microcanonical_rates'),
    ('rmgpy.pdep.network', 'Network.select_energy_grains'),
    ('rmgpy.pdep.network', 'Network.set_conditions'),
    ('rmgpy.pdep.network', 'Network.log_summary'),
    ('rmgpy.rmg.pdep', 'PDepNetwork.update'),
])
def test_network_census_refuses_even_with_cached_energy(module, scope, database):
    from rmgpy.pdep.configuration import Configuration
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    species.get_thermo_data()
    network = SimpleNamespace(isomers=[Configuration(species)], reactants=[], products=[],
                              path_reactions=[], net_reactions=[], get_all_species=lambda: [species])
    if module == 'rmgpy.pdep.network':
        from rmgpy.pdep.network import Network
        network = Network(isomers=[Configuration(species)])
    expected = (SpeciesIdentityError if module == "arkane.pdep" and scope != "PressureDependenceJob.initialize"
                else NetworkError if scope == 'PDepNetwork.update' else ExcitedSpeciesThermoError)
    match = 'unsupported network type' if scope == 'PDepNetwork.update' else None
    with pytest.raises(expected, match=match):
        api_call(module, scope, self=SimpleNamespace(network=network) if module == 'arkane.pdep' else network)


@pytest.mark.parametrize('module,scope', [
    ('rmgpy.yaml_cantera1', 'species_to_dict'),
    ('rmgpy.yaml_cantera2', 'species_to_dict'),
    ('rmgpy.yaml_cantera1', 'write_cantera'),
    ('rmgpy.yaml_cantera2', 'generate_cantera_data'),
    ('rmgpy.yaml_rms', 'obj_to_dict'),
    ('rmgpy.rmg.reactionmechanismsimulator_reactors', 'to_rms'),
    ('rmgpy.chemkin', 'write_thermo_entry'),
    ('arkane.output', 'save_thermo_lib'),
    ('rmgpy.rmg.output', 'save_output_html'),
    ('rmgpy.rmg.output', 'save_diff_html'),
    ('rmgpy.tools.diffmodels', 'enthalpy_diff'),
    ('rmgpy.tools.diffmodels', 'identical_thermo'),
    ('rmgpy.rmg.model', 'CoreEdgeReactionModel.thermo_filter_species'),
    ('rmgpy.rmg.model', 'CoreEdgeReactionModel.set_thermodynamic_filtering_parameters'),
    ('rmgpy.rmg.model', 'ReactionModel.merge'),
])
def test_export_and_filter_census_rejects_altered_source(module, scope, database, tmp_path):
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    species.get_thermo_data().CpInf = (10000, 'J/(mol*K)')
    model = SimpleNamespace(species=[species], core=SimpleNamespace(species=[species], reactions=[]),
                            edge=SimpleNamespace(species=[species], reactions=[]),
                            output_species_list=[], Tmax=1000.)
    if scope == 'ReactionModel.merge':
        from rmgpy.rmg.model import ReactionModel
        model = ReactionModel(species=[species])
    expected = (SpeciesIdentityError if scope in {"to_rms", "save_thermo_lib", "save_diff_html", "save_output_html"}
                else ExcitedSpeciesThermoError)
    with pytest.raises(expected, match='N2A'):
        api_call(module, scope, self=model, species=[species, species] if scope == 'enthalpy_diff' else species, obj=species, spcs=[species],
                 species_list=[species], species_list1=[species], species_list2=[],
                 common_species_list=[], species_pair=[species, species],
                 reaction_model=model, other=model, path=str(tmp_path / 'refused'),
                 part_core_edge='core', name='refused', common_reactions=[],
                 unique_reactions1=[], unique_reactions2=[])


def test_atom_ownership_does_not_leak_into_ground_graphs(database):
    import gc
    import weakref
    from rmgpy.thermo.state import require_atom_thermo_allowed
    resolved = Species().from_smiles('C1CCCCC1').molecule[0]
    resolved.electronic_state = 'TEST'
    atoms = resolved.atoms
    alias = resolved.copy(deep=False)
    alias.electronic_state = ''
    with pytest.raises(ExcitedSpeciesThermoError):
        require_atom_thermo_allowed(alias.atoms)
    reference = weakref.ref(resolved)
    del resolved
    gc.collect()
    assert reference() is None
    require_atom_thermo_allowed(atoms)
    # State changes are observed even when the atom list was extracted earlier.
    alias.electronic_state = 'TEST'
    with pytest.raises(ExcitedSpeciesThermoError):
        require_atom_thermo_allowed(atoms)
    alias.electronic_state = ''
    require_atom_thermo_allowed(atoms)


def test_owner_free_thermo_models_preserve_fields_in_round_trip():
    from copy import deepcopy
    from rmgpy.thermo import NASA, NASAPolynomial
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5, 0, 0, 0, 0, 100, 3],
                                             Tmin=(200, 'K'), Tmax=(6000, 'K'))],
                  Tmin=(200, 'K'), Tmax=(6000, 'K'), E0=(10, 'kJ/mol'),
                  Cp0=(29.1, 'J/(mol*K)'), CpInf=(29.1, 'J/(mol*K)'))
    restored = deepcopy(thermo)
    for field in ('E0', 'Cp0', 'CpInf'):
        assert getattr(restored, field).value_si == getattr(thermo, field).value_si
    assert restored.get_enthalpy(298.15) == thermo.get_enthalpy(298.15)


@pytest.mark.parametrize('path', ['tst', 'arkane', 'well', 'reconstruct', 'solute_ts'])
def test_kinetics_derivation_census_refuses(database, path):
    from rmgpy.reaction import Reaction
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    species.get_thermo_data()
    reaction = Reaction(reactants=[species], products=[nitrogen()])
    with pytest.raises(ExcitedSpeciesThermoError):
        if path == 'tst':
            reaction.calculate_tst_rate_coefficient(1000)
        elif path == 'arkane':
            api_call('arkane.kinetics', 'KineticsJob.generate_kinetics', self=SimpleNamespace(reaction=reaction))
        elif path == 'well':
            api_call('arkane.kinetics', 'Well.__init__', self=SimpleNamespace(), species_list=[species])
        elif path == 'reconstruct':
            api_call('rmgpy.data.kinetics.database', 'KineticsDatabase.reconstruct_kinetics_from_source',
                     self=None, reaction=reaction)
        else:
            api_call('rmgpy.data.kinetics.family', 'get_site_solute_data', rxn=reaction)


@pytest.mark.parametrize('path', ['partition', 'density', 'sum', 'generate', 'configuration'])
def test_statmech_consumers_census_refuses(database, path):
    import numpy as np
    from rmgpy.pdep.configuration import Configuration
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    species.get_thermo_data()
    with pytest.raises(ExcitedSpeciesThermoError):
        if path == 'partition':
            species.get_partition_function(1000)
        elif path == 'density':
            species.get_density_of_states(np.array([0., 1000.]))
        elif path == 'sum':
            species.get_sum_of_states(np.array([0., 1000.]))
        elif path == 'generate':
            species.generate_statmech()
        else:
            Configuration(species).calculate_density_of_states(np.array([0., 1000.]))


def test_surface_census_refuses_derived_coverage(database):
    from rmgpy.solver.surface import SurfaceReactor
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    species.get_thermo_data()
    reactor = SurfaceReactor(T=(600, 'K'), P_initial=(1, 'bar'),
                             initial_gas_mole_fractions={species: 1.0},
                             initial_surface_coverages={}, surface_volume_ratio=(1, 'm^-1'),
                             surface_site_density=(2.5e-5, 'mol/m^2'), thermo_coverage_dependence=True, n_sims=1)
    with pytest.raises(ExcitedSpeciesThermoError, match='N2A'):
        reactor.initialize_model([species], [], [], [])


def test_processing_preserves_explicit_library_energy(database):
    library = database.libraries['ExactStateFixture']
    library.entries['N2A'].data.E0 = (600, 'kJ/mol')
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    assert species.get_thermo_data().E0.value_si == 600000
    from rmgpy.pdep.configuration import Configuration
    assert Configuration(species).E0 == 600000


def test_processing_rejects_altered_attached_energy_before_mutation(database):
    from rmgpy.thermo.thermoengine import process_thermo_data
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    thermo = species.get_thermo_data()
    initial = species.conformer.E0.value_si
    thermo.E0 = (1234, 'kJ/mol')
    with pytest.raises(ExcitedSpeciesThermoError, match='N2A'):
        process_thermo_data(species, thermo)
    assert species.conformer.E0.value_si == initial


def test_ground_library_export_reload_preserves_thermo(database, tmp_path):
    from arkane.output import save_thermo_lib
    species = nitrogen()
    original = species.get_thermo_data()
    values = [(original.get_heat_capacity(t), original.get_enthalpy(t), original.get_entropy(t))
              for t in (298.15, 1000., 3000.)]
    save_thermo_lib([species], str(tmp_path), 'ground', 'test')
    library = ThermoLibrary(label='Ground')
    library.load(str(tmp_path / 'ground.py'), database.local_context, {})
    assert library.entries['N2'].item.vibrational_level == -1
    database.libraries = {'Ground': library}
    database.library_order = ['Ground']
    restored = nitrogen().get_thermo_data()
    assert values == [(restored.get_heat_capacity(t), restored.get_enthalpy(t), restored.get_entropy(t))
                      for t in (298.15, 1000., 3000.)]


def test_declared_atom_list_keeps_manifold_refusal(database):
    from rmgpy.rmg.model import CoreEdgeReactionModel
    from rmgpy.thermo.state import require_atom_thermo_allowed
    from rmgpy.exceptions import VibrationalManifoldError
    species = Species(label='cycle').from_smiles('C1CCCCC1')
    CoreEdgeReactionModel().declare_vibrational_manifold(species)
    atoms = species.molecule[0].get_smallest_set_of_smallest_rings()[0]
    with pytest.raises(VibrationalManifoldError, match='cycle'):
        require_atom_thermo_allowed(atoms)


def test_arkane_yaml_does_not_prefer_unkeyed_smiles_over_resolved_header(database, tmp_path):
    import yaml
    from arkane.common import ArkaneSpecies
    species = nitrogen()
    species.get_thermo_data()
    path = tmp_path / 'species' / 'N2.yml'
    record = ArkaneSpecies(species=species)
    record.save_yaml(str(tmp_path))
    # save_yaml accepts an output directory and uses species/<label>.yml.
    if not path.is_file():
        path = next(tmp_path.rglob('*.yml'))
    payload = yaml.safe_load(path.read_text())
    payload['adjacency_list'] = 'electronicstate A3Su+\n' + N2
    payload['smiles'] = 'N#N'
    path.write_text(yaml.safe_dump(payload))
    with pytest.raises(ExcitedSpeciesThermoError):
        from rmgpy.statmech import Conformer
        ArkaneSpecies(conformer=Conformer()).load_yaml(str(path))


def test_cached_arkane_well_refuses_after_state_mutation(database):
    from arkane.kinetics import Well
    species = nitrogen()
    species.get_thermo_data()
    well = Well([species])
    before = well.E0
    well.E0 = before + 1
    assert well.E0 == before + 1
    species.molecule[0].electronic_state = 'A3Su+'
    with pytest.raises(ExcitedSpeciesThermoError):
        _ = well.E0


def test_isodesmic_acquisition_refuses_resolved_graph(database):
    from arkane.encorr.isodesmic import ErrorCancelingSpecies
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    with pytest.raises(ExcitedSpeciesThermoError):
        ErrorCancelingSpecies(species.molecule[0], (10, 'kJ/mol'), None)


def test_configuration_accepts_exact_attached_thermo_without_conformer(database):
    from copy import deepcopy
    from rmgpy.pdep.configuration import Configuration
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    species.thermo = deepcopy(database.libraries['ExactStateFixture'].entries['N2A'].data)
    assert species.conformer is None
    assert Configuration(species).E0 == pytest.approx(582013.0788860451)


def test_refused_solvent_request_does_not_poison_gas_source(database):
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    thermo = species.get_thermo_data()
    with pytest.raises(ExcitedSpeciesThermoError):
        species.get_thermo_data('water')
    assert species.get_thermo_data() is thermo


def test_inconsistent_resonance_states_refuse_before_library_matching(database):
    species = nitrogen()
    species.molecule.append(species.molecule[0].copy(deep=True))
    species.molecule[1].electronic_state = 'A3Su+'
    for call in (species.get_thermo_data, lambda: database.get_thermo_data(species)):
        with pytest.raises(ExcitedSpeciesThermoError, match='inconsistent'):
            call()


@pytest.mark.parametrize('deep', [False, True])
def test_declared_graph_copy_keeps_estimator_refusal(database, deep):
    from rmgpy.thermo.state import require_atom_thermo_allowed, require_thermo_estimation_allowed
    molecule = nitrogen().molecule[0]
    molecule.props['vibrational_manifold'] = 'N2'
    copied = molecule.copy(deep=deep)
    assert copied.props['vibrational_manifold'] == 'N2'
    with pytest.raises(VibrationalManifoldError):
        require_thermo_estimation_allowed(copied)
    with pytest.raises(VibrationalManifoldError):
        require_atom_thermo_allowed(copied.vertices)
    wrapped = Species(molecule=[copied])
    with pytest.raises(VibrationalManifoldError):
        wrapped.get_thermo_data()


def test_atom_registration_survives_graph_population_after_state_assignment(database):
    from rmgpy.molecule import Atom, Molecule
    from rmgpy.thermo.state import require_atom_thermo_allowed
    molecule = Molecule(electronic_state='excited')
    atom = Atom('C')
    molecule.add_atom(atom)
    with pytest.raises(ExcitedSpeciesThermoError):
        require_atom_thermo_allowed([atom])
    ring = Molecule(smiles='C1CCCCC1')
    ring.electronic_state = 'excited'
    cycles = ring.get_smallest_set_of_smallest_rings()
    with pytest.raises(ExcitedSpeciesThermoError):
        require_atom_thermo_allowed(cycles[0])


@pytest.mark.parametrize('module, name', [
    ('rmgpy.pdep.msc', 'apply_modified_strong_collision_method'),
    ('rmgpy.pdep.rs', 'apply_reservoir_state_method'),
    ('rmgpy.pdep.cse', 'apply_chemically_significant_eigenvalues_method'),
    ('rmgpy.pdep.cse', 'get_rate_coefficients_CSE_Advanced'),
    ('rmgpy.pdep.cse', 'apply_chemically_significant_eigenvalues_method_georgievskii'),
    ('rmgpy.pdep.me', 'generate_full_me_matrix'),
    ('rmgpy.pdep.me', 'states_to_configurations'),
    ('rmgpy.pdep.sls', 'apply_simulation_least_squares_method'),
])
def test_compiled_network_consumers_refuse_cached_ground_arrays(database, module, name):
    import importlib
    from rmgpy.pdep.configuration import Configuration
    species = nitrogen()
    species.molecule[0].electronic_state = 'A3Su+'
    network = SimpleNamespace(isomers=[Configuration(species)], reactants=[], products=[], path_reactions=[])
    function = getattr(importlib.import_module(module), name)
    args = [network]
    if name == 'get_rate_coefficients_CSE_Advanced':
        args += [1000., 1e5]
    elif name == 'states_to_configurations':
        args += [None, None]
    with pytest.raises(NetworkError, match='unsupported network type'):
        function(*args)


@pytest.mark.parametrize('method', ['get_reference_enthalpy', 'to_error_canceling_spcs'])
def test_reference_enthalpy_refuses_mutated_stored_identity(database, method):
    from arkane.encorr.reference import ReferenceSpecies
    reference = ReferenceSpecies(species=nitrogen())
    reference.adjacency_list = 'electronicstate A3Su+\n' + N2
    with pytest.raises(ExcitedSpeciesThermoError):
        getattr(reference, method)(*([None] if method == 'to_error_canceling_spcs' else []))


@pytest.mark.parametrize('method', ['checkSpecies', 'printThermo', 'printSpeciesComments'])
def test_checkmodels_uses_checked_source(database, method):
    from scripts import checkModels
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    species.get_thermo_data().CpInf = (10000, 'J/(mol*K)')
    with pytest.raises(ExcitedSpeciesThermoError):
        if method == 'checkSpecies':
            checkModels.checkSpecies([(species, species)], [], [])
        else:
            getattr(checkModels, method)(species)


def test_transport_refreshes_after_public_state_mutation(database, monkeypatch):
    species = nitrogen()
    ground = object()
    excited = object()
    species.transport_data = ground
    assert species.get_transport_data() is ground
    transport = SimpleNamespace(get_transport_properties=lambda owner: (excited, None, None))
    monkeypatch.setattr(data_module.database, 'transport', transport, raising=False)
    species.molecule[0].electronic_state = 'A3Su+'
    assert species.get_transport_data() is excited


@pytest.mark.parametrize('attribute', ['ref_data', 'calc_data', 'bac_data'])
def test_bac_energy_caches_refuse_after_reference_state_mutation(database, attribute):
    from arkane.encorr.data import BACDatapoint, BACDataset
    from arkane.encorr.reference import ReferenceSpecies
    reference = ReferenceSpecies(species=nitrogen())
    datapoint = BACDatapoint(reference, level_of_theory=object())
    setattr(datapoint, '_' + attribute, 10.)
    dataset = BACDataset([datapoint])
    assert getattr(dataset, attribute)[0] == 10.
    reference.adjacency_list = 'electronicstate A3Su+\n' + N2
    with pytest.raises(ExcitedSpeciesThermoError):
        getattr(datapoint, attribute)
    with pytest.raises(ExcitedSpeciesThermoError):
        getattr(dataset, attribute)


@pytest.mark.parametrize('method', ['get_correction', '_get_petersson_correction', '_get_melius_correction'])
def test_bac_corrections_refuse_resolved_reference(database, method):
    from arkane.encorr.bac import BAC
    reference = SimpleNamespace(adjacency_list='electronicstate A3Su+\n' + N2)
    with pytest.raises(ExcitedSpeciesThermoError):
        getattr(BAC, method)(SimpleNamespace(), datapoint=SimpleNamespace(spc=reference))


def test_wilhoit_processing_preserves_energy_or_refuses_named(database):
    from copy import deepcopy
    from rmgpy.thermo import Wilhoit
    from rmgpy.thermo.thermoengine import process_thermo_data
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    entry = database.libraries['ExactStateFixture'].entries['N2A']
    entry.data = entry.data.to_wilhoit()
    result = process_thermo_data(species, deepcopy(entry.data), thermo_class=Wilhoit)
    assert isinstance(result, Wilhoit)
    assert result.E0.value_si == pytest.approx(entry.data.E0.value_si)


def test_cantera_refreshes_transport_after_state_mutation(database, monkeypatch):
    from rmgpy.transport import TransportData
    ground = TransportData(shapeIndex=1, epsilon=(100., 'K'), sigma=(3., 'angstrom'),
                           dipoleMoment=(0., 'De'), polarizability=(0., 'angstrom^3'), rotrelaxcollnum=1.)
    excited = TransportData(shapeIndex=1, epsilon=(200., 'K'), sigma=(4., 'angstrom'),
                            dipoleMoment=(0., 'De'), polarizability=(0., 'angstrom^3'), rotrelaxcollnum=1.)
    species = nitrogen()
    species.transport_data = ground
    species.get_transport_data()
    monkeypatch.setattr(data_module.database, 'transport',
                        SimpleNamespace(get_transport_properties=lambda owner: (excited, None, None)), raising=False)
    species.molecule[0].electronic_state = 'A3Su+'
    result = species.to_cantera(use_chemkin_identifier=True)
    assert species.transport_data is excited
    assert result.transport.diameter == pytest.approx(4e-10)


@pytest.mark.parametrize('route', ['database', 'library', 'numeric', 'wilhoit'])
def test_legacy_thermo_writers_refuse_lossy_resolved_records(database, tmp_path, route):
    from copy import deepcopy
    from io import StringIO
    from rmgpy.data.thermo import save_entry
    library = database.libraries['ExactStateFixture']
    library.entries = {'N2A': deepcopy(library.entries['N2A'])}
    with pytest.raises(ExcitedSpeciesThermoError):
        if route == 'database':
            database.save_old(str(tmp_path / 'legacy'))
        elif route == 'library':
            library.save_old(str(tmp_path / 'dictionary'), '', str(tmp_path / 'library'))
        elif route == 'numeric':
            library.save_old_library(str(tmp_path / 'library'))
        else:
            entry = library.entries['N2A']
            entry.data = entry.data.to_wilhoit()
            save_entry(StringIO(), entry)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('record_type', ['arkane', 'reference'])
def test_arkane_record_writers_refuse_mutated_stored_identity(database, tmp_path, record_type):
    from arkane.common import ArkaneSpecies
    from arkane.encorr.reference import ReferenceSpecies
    cls = ArkaneSpecies if record_type == 'arkane' else ReferenceSpecies
    record = cls(species=nitrogen())
    record.adjacency_list = 'electronicstate A3Su+\n' + N2
    expected = SpeciesIdentityError if record_type == 'reference' else ExcitedSpeciesThermoError
    with pytest.raises(expected):
        record.save_yaml(str(tmp_path))
    assert list(tmp_path.iterdir()) == []


def test_native_weakrefs_stay_out_of_complete_transport(database):
    import weakref
    from rmgpy.data.kinetics.family import object_state, complete_round_trip
    molecule = nitrogen('vibrationallevel 1\n').molecule[0]
    reference = weakref.ref(molecule)
    assert '__weakref__' not in object_state(molecule)
    copied = complete_round_trip(molecule)
    assert copied is not molecule
    assert copied.vibrational_level == 1
    assert reference() is molecule
    assert all(value is not reference for value in weakref.getweakrefs(copied))


@pytest.mark.parametrize('electronic', ['', 'A3Su+'])
def test_unnamed_manifold_refuses_before_mutating_declaration(electronic):
    from rmgpy.molecule import Molecule
    from rmgpy.rmg.model import CoreEdgeReactionModel
    from rmgpy.thermo import ThermoData
    mol = Molecule(smiles='N#N')
    mol.electronic_state = electronic
    thermo = ThermoData()
    species = Species(molecule=[mol], thermo=thermo)
    before_species = dict(species.props)
    before_molecule = dict(mol.props)
    model = CoreEdgeReactionModel()
    with pytest.raises(VibrationalManifoldError, match='non-empty input label'):
        model.declare_vibrational_manifold(species)
    assert model.vibrational_manifolds == []
    assert species.props == before_species and mol.props == before_molecule
    assert species.thermo is thermo


def test_nasa9_attachment_cannot_hide_between_five_cp_samples(database):
    from copy import deepcopy
    import numpy as np
    from rmgpy.thermo import NASAPolynomial
    from rmgpy.solver.plasma import _thermo_comparison_plan
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    reference = species.get_thermo_data()
    actual = deepcopy(reference)
    _, _, points = _thermo_comparison_plan(actual, reference)['segments'][0]
    perturbation = np.polynomial.polynomial.polyfromroots(np.array(points) / 1000.)
    coefficients = np.r_[0., 0., reference.polynomials[0].coeffs]
    for i, coefficient in enumerate(perturbation):
        coefficients[i] += -0.02 * coefficient * 1000.**(2 - i)
    polynomial = NASAPolynomial(coeffs=coefficients, Tmin=reference.polynomials[0].Tmin,
                               Tmax=reference.polynomials[0].Tmax)
    midpoint = points[len(points) // 2]
    polynomial.change_base_enthalpy(reference.get_enthalpy(midpoint) - polynomial.get_enthalpy(midpoint))
    polynomial.change_base_entropy(reference.get_entropy(midpoint) - polynomial.get_entropy(midpoint))
    actual.polynomials = [polynomial]
    assert actual.get_heat_capacity(298.15) == pytest.approx(178.1756377716229)
    species.thermo = actual
    with pytest.raises(ExcitedSpeciesThermoError):
        species.get_heat_capacity(298.15)


def test_chemkin_manifold_roundtrip_retains_library_only_v0(database, tmp_path):
    from rmgpy import chemkin
    from rmgpy.rmg.model import CoreEdgeReactionModel
    species = nitrogen()
    species.index = 1
    CoreEdgeReactionModel().declare_vibrational_manifold(species)
    fixed = species.get_heat_capacity(298.15)
    chemkin.save_chemkin_file(str(tmp_path / 'chem.inp'), [species], [])
    chemkin.save_species_dictionary(str(tmp_path / 'dictionary.txt'), [species])
    loaded = chemkin.load_chemkin_file(str(tmp_path / 'chem.inp'), str(tmp_path / 'dictionary.txt'))[0][0]
    assert loaded.molecule[0].vibrational_level == 0
    assert not loaded.is_isomorphic(nitrogen())
    assert loaded.get_heat_capacity(298.15) == pytest.approx(fixed, rel=1e-6)
    database.libraries['ExactStateFixture'].entries.pop('N2v0')
    with pytest.raises(ExcitedSpeciesThermoError):
        loaded.get_heat_capacity(298.15)


@pytest.mark.parametrize('change', ['class', 'count', 'range', 'coefficient', 'zero', 'coverage'])
def test_exact_thermo_structure_rejects_changed_library_fields(database, change):
    from copy import deepcopy
    from rmgpy.thermo.state import thermo_fields_match
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    reference = species.get_thermo_data()
    actual = deepcopy(reference)
    assert thermo_fields_match(actual, reference)
    if change == 'class':
        actual = reference.to_wilhoit()
    elif change == 'count':
        actual.polynomials = actual.polynomials * 2
    elif change == 'range':
        actual.polynomials[0].Tmin = (210, 'K')
    elif change in ('coefficient', 'zero'):
        coefficients = actual.polynomials[0].coeffs.copy()
        coefficients[0 if change == 'coefficient' else 1] += 1e-7
        actual.polynomials[0].coeffs = coefficients
    else:
        actual.thermo_coverage_dependence = {N2: {'model': 'polynomial', 'enthalpy-coefficients': [(1., 'J/mol')],
                                                 'entropy-coefficients': [(0., 'J/(mol*K)')]}}
    assert not thermo_fields_match(actual, reference)


def test_exact_thermo_structure_has_coefficient_tolerance_and_named_unknown_refusal(database):
    from copy import deepcopy
    from rmgpy.thermo.state import thermo_fields_match
    reference = nitrogen('electronicstate A3Su+\n', 'N2A').get_thermo_data()
    actual = deepcopy(reference)
    coefficients = actual.polynomials[0].coeffs.copy()
    coefficients[0] *= 1 + 5e-10
    actual.polynomials[0].coeffs = coefficients
    assert thermo_fields_match(actual, reference)
    coefficients[0] *= 1 + 2e-9
    actual.polynomials[0].coeffs = coefficients
    assert not thermo_fields_match(actual, reference)
    with pytest.raises(ExcitedSpeciesThermoError, match='Cannot inspect thermo model classes'):
        thermo_fields_match(object(), reference)


@pytest.mark.parametrize('attachment', ['unknown', 'polynomial', 'future_unknown', 'future_library'])
def test_attached_model_inspection_refuses_unknown_classes_by_name(database, attachment):
    from concurrent.futures import Future
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    reference = species.get_thermo_data()
    model = reference.polynomials[0] if attachment == 'polynomial' else object()
    if attachment.startswith('future'):
        future = Future()
        future.set_result(reference if attachment == 'future_library' else model)
        model = future
    species.thermo = model
    if attachment == 'future_library':
        assert species.get_thermo_data() is reference
    else:
        with pytest.raises(ExcitedSpeciesThermoError):
            species.get_thermo_data()
