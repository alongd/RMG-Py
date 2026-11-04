"""Output-only bounded Arrhenius fits for tabulated kinetics."""

import pytest

from rmgpy.kinetics import KineticsData
from rmgpy.kinetics.export import fit_kinetics_data_for_export


def test_fit_uses_only_source_grid_sets_bounds_and_reports_error_without_mutating_table():
    table = KineticsData(
        Tdata=([300, 400, 500, 700, 1000], "K"),
        kdata=([1e-18, 1.5e-18, 2.0e-18, 3.2e-18, 5e-18], "m^3/(mol*s)"),
        Tmin=(250, "K"), Tmax=(1500, "K"),
        comment="Source: fixture table; SHA-256 abc",
    )

    fit = fit_kinetics_data_for_export(table, n_gas=2)

    assert (fit.Tmin.value_si, fit.Tmax.value_si) == (300, 1000)
    assert "Source: fixture table" in fit.comment
    assert "maximum relative fit error" in fit.comment
    assert "The source table is the record" in fit.comment
    assert isinstance(table, KineticsData)
    assert (table.Tmin.value_si, table.Tmax.value_si) == (250, 1500)


def test_non_arrhenius_table_is_refused_when_fit_error_exceeds_25_percent():
    table = KineticsData(
        Tdata=([300, 400, 500, 600, 700, 800], "K"),
        kdata=([1e-18, 1e-14, 1e-18, 1e-14, 1e-18, 1e-14], "m^3/(mol*s)"),
        Tmin=(300, "K"), Tmax=(800, "K"),
    )

    import pickle
    before = pickle.dumps(table)
    with pytest.raises(ValueError, match="exceeds 25%"):
        fit_kinetics_data_for_export(table, n_gas=2)
    assert pickle.dumps(table) == before


def test_positive_sum_fit_is_deterministic_bounded_and_preserves_source_bytes():
    import pickle
    import numpy as np
    from rmgpy.kinetics import Arrhenius, MultiArrhenius
    temperatures = np.linspace(300., 3500., 65)
    terms = [Arrhenius(A=(2e-12, 'm^3/(mol*s)'), n=1.2, Ea=(1000., 'J/mol')),
             Arrhenius(A=(1e-5, 'm^3/(mol*s)'), n=-0.5, Ea=(40000., 'J/mol'))]
    values = [sum(term.get_rate_coefficient(float(t)) for term in terms) for t in temperatures]
    table = KineticsData(Tdata=(temperatures, 'K'), kdata=(values, 'm^3/(mol*s)'))
    before = pickle.dumps(table)
    fit = fit_kinetics_data_for_export(table, 2)
    assert isinstance(fit, MultiArrhenius)
    assert len(fit.arrhenius) == 2
    assert repr(fit_kinetics_data_for_export(table, 2)) == repr(fit)
    assert pickle.dumps(table) == before
    assert all(term.A.value_si > 0 for term in fit.arrhenius)
    assert all((term.Tmin.value_si, term.Tmax.value_si) == (300., 3500.) for term in fit.arrhenius)
    assert max(abs(fit.get_rate_coefficient(float(t)) / k - 1) for t, k in zip(temperatures, values)) <= .1


@pytest.fixture(params=[30, 0])
def admitted_export_reaction(monkeypatch, request):
    """Native writer/reader coverage: V-V single term and V-T DUPLICATE sum."""
    index = request.param
    from pathlib import Path
    import pickle
    from rmgpy import settings
    import rmgpy.data.rmg as data_module
    from rmgpy.data.rmg import RMGDatabase
    from rmgpy.data.kinetics.database import KineticsDatabase
    from rmgpy.data.thermo import ThermoDatabase
    from rmgpy.data.solvation import SolvationDatabase
    from rmgpy.rmg.model import CoreEdgeReactionModel
    from rmgpy.species import Species
    from rmgpy.kinetics import MultiArrhenius
    from rmgpy.thermo.thermoengine import process_thermo_data

    database_path = Path(settings['database.directory'])
    hub = RMGDatabase()
    hub.kinetics = KineticsDatabase()
    hub.kinetics.load_libraries(str(database_path / 'kinetics/libraries'), ['N2HeavyParticleRelaxation'])
    hub.thermo = ThermoDatabase()
    hub.thermo.load_libraries(str(database_path / 'thermo/libraries'), ['N2Vibrational'])
    hub.solvation = SolvationDatabase()
    monkeypatch.setattr(data_module, 'database', hub)
    model = CoreEdgeReactionModel()
    model.declare_vibrational_manifold(Species(label='N2').from_smiles('N#N'))
    reaction = hub.kinetics.libraries['N2HeavyParticleRelaxation'].get_library_reactions()[index]
    for species in reaction.reactants + reaction.products:
        model._apply_vibrational_manifold(species)
        species.thermo = hub.thermo.get_thermo_data(species)
    reaction, _ = model.make_new_reaction(reaction, generate_thermo=False)
    species_list = list(dict.fromkeys(reaction.reactants + reaction.products))
    for species in species_list:
        species.thermo = process_thermo_data(species, species.thermo)
    reaction.generate_pairs()
    assert reaction.pairs
    return index, reaction, species_list


def test_real_admitted_annotated_chemkin_round_trip_with_flux_pairs(tmp_path, admitted_export_reaction):
    import pickle
    from rmgpy.kinetics import MultiArrhenius
    from rmgpy.chemkin import save_chemkin_file, save_species_dictionary, load_chemkin_file
    index, reaction, species_list = admitted_export_reaction
    snapshot = pickle.dumps(reaction.kinetics)
    fit = fit_kinetics_data_for_export(reaction.kinetics, 2)
    if index == 0:
        assert isinstance(fit, MultiArrhenius)
    mechanism, dictionary = tmp_path / 'chem.inp', tmp_path / 'dictionary.txt'
    save_chemkin_file(str(mechanism), species_list, [reaction])
    save_species_dictionary(str(dictionary), species_list)
    text = mechanism.read_text()
    assert 'Flux pairs:' in text and 'Output-only' in text
    if index == 0:
        assert 'DUPLICATE' in text
    loaded_species, loaded_reactions = load_chemkin_file(str(mechanism), str(dictionary))
    assert loaded_species and len(loaded_reactions) == 1
    assert loaded_reactions[0].pairs
    if index == 0:
        assert isinstance(loaded_reactions[0].kinetics, MultiArrhenius)
    for temperature in (300., 1000., 3500.):
        assert loaded_reactions[0].get_rate_coefficient(temperature) == pytest.approx(
            fit.get_rate_coefficient(temperature), rel=.01, abs=0)
    assert pickle.dumps(reaction.kinetics) == snapshot


@pytest.mark.parametrize('writer', ['cantera1', 'cantera2'])
def test_real_admitted_cantera_round_trip(tmp_path, admitted_export_reaction, writer):
    import pickle
    import cantera as ct
    import math
    from rmgpy.rmg.model import ReactionModel
    from rmgpy.data.transport import TransportDatabase
    from rmgpy.yaml_cantera1 import write_cantera
    from rmgpy.yaml_cantera2 import save_cantera_model
    from rmgpy import settings
    from pathlib import Path
    index, reaction, species_list = admitted_export_reaction
    transport = TransportDatabase()
    transport.load(str(Path(settings['database.directory']) / 'transport'), ['PrimaryTransportLibrary'])
    # Writer1 requires gas transport entries; request them from the loaded DB.
    for species in species_list:
        species.transport_data = transport.get_transport_properties(species)[0]
    model = ReactionModel(species=species_list, reactions=[reaction])
    path = tmp_path / 'chem.yaml'
    snapshot = pickle.dumps(reaction.kinetics)
    if writer == 'cantera1':
        write_cantera(species_list, [reaction], model.get_elements(), path=str(path))
    else:
        save_cantera_model(model, str(path))
    gas = ct.Solution(str(path))
    fit = fit_kinetics_data_for_export(reaction.kinetics, 2)
    assert gas.n_reactions == len(getattr(fit, 'arrhenius', [fit]))
    for temperature in (300., 1000., 3500.):
        gas.TP = temperature, 1e5
        # RMG and Cantera use slightly different numerical gas constants.
        # Compare the emitted Arrhenius parameters using Cantera's own R.
        expected = sum(term.A.value_si * (temperature / term.T0.value_si) ** term.n.value_si
                       * math.exp(-term.Ea.value_si * 1000. / (ct.gas_constant * temperature))
                       for term in getattr(fit, 'arrhenius', [fit])) * 1000.
        assert sum(gas.forward_rate_constants) == pytest.approx(expected, rel=1e-10, abs=0)
    assert 'The source table is the record' in path.read_text()
    assert pickle.dumps(reaction.kinetics) == snapshot


@pytest.mark.parametrize('conflict_at', [None, 1, 2])
def test_duplicate_merge_preserves_agreed_flux_pairs_and_refuses_conflicts(conflict_at):
    from rmgpy.chemkin import _process_duplicate_reactions
    from rmgpy.data.kinetics.library import LibraryReaction
    from rmgpy.exceptions import ChemkinError
    from rmgpy.kinetics import Arrhenius, MultiArrhenius
    from rmgpy.species import Species

    reactant = Species(label='H').from_smiles('[H]')
    product = Species(label='H2').from_smiles('[H][H]')
    pairs = [(reactant, product)]
    leaves = [LibraryReaction(reactants=[reactant], products=[product],
                              kinetics=Arrhenius(A=(1., 's^-1')), library='fixture',
                              pairs=list(pairs), duplicate=True) for _ in range(3)]
    if conflict_at is not None:
        leaves[conflict_at].pairs = []
        with pytest.raises(ChemkinError, match='Conflicting flux pairs for duplicate reaction'):
            _process_duplicate_reactions(leaves)
    else:
        _process_duplicate_reactions(leaves)
        assert len(leaves) == 1
        assert leaves[0].pairs == pairs
        assert isinstance(leaves[0].kinetics, MultiArrhenius)
        assert len(leaves[0].kinetics.arrhenius) == 3


def test_default_translation_retains_resolved_state_refusal(tmp_path, admitted_export_reaction):
    from rmgpy.exceptions import MechanismWriterError, SpeciesIdentityError
    from rmgpy.rmg.main import RMG
    from rmgpy.rmg.model import CoreEdgeReactionModel, ReactionModel
    from rmgpy.rmg.settings import WriterConfig

    _, reaction, species_list = admitted_export_reaction
    job = RMG(output_directory=str(tmp_path))
    assert job.generate_cantera_from_chemkin is True
    job.chemkin_writer_config = WriterConfig(save_interval=1)
    job.reaction_model = CoreEdgeReactionModel()
    job.reaction_model.core = ReactionModel(species=species_list, reactions=[reaction])
    job.generate_end_of_run_cantera_files()
    assert len(job.export_failures) == 1
    step, error = job.export_failures[0]
    assert step == 'Chemkin-to-Cantera translation (cantera_from_ck/)'
    assert type(error) is SpeciesIdentityError
    assert str(error).startswith(
        'Resolved electronic or vibrational states are not supported by '
        'rmgpy/rmg/main.py:RMG.generate_cantera_files_from_chemkin:')
    with pytest.raises(MechanismWriterError, match='Chemkin-to-Cantera translation'):
        job.report_export_failures()
