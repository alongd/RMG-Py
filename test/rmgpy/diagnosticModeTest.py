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

"""
The prescribed-Te transport diagnostic must declare itself in its artifacts, and
must refuse the one claim it cannot make.

Every artifact assertion here reads the file back **off disk**. Asserting that a
writer was handed the right string tests the writer's argument list, not whether
the string reached the file -- and the artifact is the whole point: the two
external reports this work exists because of were drawn from files, not from a
Python session.
"""

import csv
import os
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from rmgpy.chemkin import save_chemkin_file, save_species_dictionary, save_transport_file
from rmgpy.diagnostic_mode import (PRESCRIBED_TE_DISCLAIMER, PRESCRIBED_TE_DISCLAIMER_ASCII,
                                   PRESCRIBED_TE_MODE, SOLVED_TE_MODE, SOLVED_TE_DISCLAIMER,
                                   EEDF_MODE, EEDF_DISCLAIMER, MAXWELLIAN_MODE, MAXWELLIAN_DISCLAIMER,
                                   comment_block, run_provenance_lines)
from rmgpy.exceptions import PlasmaStateError, PrescribedTeDiagnosticError
from rmgpy.kinetics import Arrhenius, TwoTemperaturePlasma
from rmgpy.reaction import Reaction
from rmgpy.rmg.listener import SimulationProfileWriter
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.solver.simple import SimpleReactor
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial
from rmgpy.tools.plot import parse_csv_data, parse_csv_provenance

T_GAS = 300.0
T_E = 10000.0
P0 = (1.0, 'bar')


def _thermo(a0, a6):
    """A NASA7 polynomial; Chemkin and Cantera export both require NASA, not ThermoData."""
    coeffs = [a0, 0.0, 0.0, 0.0, 0.0, 0.0, a6]
    return NASA(polynomials=[
        NASAPolynomial(coeffs=coeffs, Tmin=(298.0, 'K'), Tmax=(1000.0, 'K')),
        NASAPolynomial(coeffs=coeffs, Tmin=(1000.0, 'K'), Tmax=(3000.0, 'K')),
    ], Tmin=(298.0, 'K'), Tmax=(3000.0, 'K'))


def _species():
    electron = Species(label='e-', index=1).from_adjacency_list('1 e u1 p0 c-1')
    electron.thermo = _thermo(2.5, 12.0)
    ar = Species(label='Ar', index=2).from_adjacency_list('1 Ar u0 p4 c0')
    ar.thermo = _thermo(2.5, 4.4)
    ar_ion = Species(label='Ar+', index=3).from_adjacency_list('1 Ar u1 p3 c+1')
    ar_ion.thermo = _thermo(2.5, 4.9)
    return electron, ar, ar_ion


def _plasma_reactor(termination=None):
    electron, ar, ar_ion = _species()
    imf = {electron: 1.0e-8, ar: 1.0, ar_ion: 1.0e-8}
    return PlasmaReactor(T_GAS, P0, imf, (T_E, 'K'), n_sims=1,
                         termination=termination if termination is not None else [])


def _simple_reactor():
    _, ar, _ = _species()
    return SimpleReactor(T=T_GAS, P=P0, initial_mole_fractions={ar: 1.0}, n_sims=1, termination=[])


class _FakeReactionSystem(object):
    """A reaction system carrying nothing but snapshots, for the csv writer."""

    def __init__(self, snapshots):
        self.snapshots = snapshots


class _PrescribedTeSystem(_FakeReactionSystem):
    diagnostic_mode = PRESCRIBED_TE_MODE


def _write_profile(tmp_path, reaction_system, core_species):
    os.makedirs(os.path.join(str(tmp_path), 'solver'), exist_ok=True)
    SimulationProfileWriter(str(tmp_path), 0, core_species).update(reaction_system)
    return os.path.join(str(tmp_path), 'solver',
                        'simulation_1_{0:d}.csv'.format(len(core_species)))


class DiagnosticModeStringsTest:
    """The canonical strings, and the query every artifact writer uses."""

    def test_mode_name_is_the_exact_ruling_string(self):
        assert PRESCRIBED_TE_MODE == 'LEGACY: prescribed-Te transport diagnostic'

    def test_disclaimer_is_verbatim_with_the_u2091_subscript(self):
        assert PRESCRIBED_TE_DISCLAIMER == (
            'Species evolution conditional on prescribed Te and nₑ; '
            'absolute electron density and ionisation degree are not predicted.')

    def test_ascii_variant_differs_only_in_the_subscript(self):
        assert PRESCRIBED_TE_DISCLAIMER_ASCII == PRESCRIBED_TE_DISCLAIMER.replace('nₑ', 'n_e')
        PRESCRIBED_TE_DISCLAIMER_ASCII.encode('ascii')  # raises if it is not ASCII

    def test_plasma_reactor_declares_the_mode_itself(self):
        assert _plasma_reactor().diagnostic_mode == PRESCRIBED_TE_MODE

    def test_plasma_reactor_identity_names_the_mode(self):
        assert PRESCRIBED_TE_MODE in _plasma_reactor()._identity()

    def test_non_plasma_reactor_declares_nothing(self):
        assert run_provenance_lines(_simple_reactor()) == []
        assert comment_block(_simple_reactor(), '#') == ''

    def test_provenance_lines_carry_mode_then_disclaimer(self):
        lines = run_provenance_lines(_plasma_reactor())
        assert lines == ['RMG diagnostic mode: ' + PRESCRIBED_TE_MODE, PRESCRIBED_TE_DISCLAIMER]


class SimulationCsvLabelTest:
    """The artifact the 0.9988 figure was drawn from."""

    def test_plasma_csv_on_disk_carries_mode_and_disclaimer(self, tmp_path):
        electron, ar, ar_ion = _species()
        system = _PrescribedTeSystem([[0.0, 1.0, 1e-8, 1.0, 1e-8], [1.0, 1.0, 2e-8, 1.0, 2e-8]])
        path = _write_profile(tmp_path, system, [electron, ar, ar_ion])

        with open(path, encoding='utf-8') as f:
            text = f.read()

        assert PRESCRIBED_TE_MODE in text
        assert PRESCRIBED_TE_DISCLAIMER in text
        assert text.startswith('#')

    def test_non_plasma_csv_on_disk_is_exactly_the_old_bytes(self, tmp_path):
        electron, ar, ar_ion = _species()
        system = _FakeReactionSystem([[0.0, 1.0, 1e-8, 1.0, 1e-8]])
        path = _write_profile(tmp_path, system, [electron, ar, ar_ion])

        with open(path, encoding='utf-8') as f:
            text = f.read()

        assert '#' not in text
        assert text.splitlines()[0] == 'Time (s),Volume (m^3),e-(1),Ar(2),Ar(3)'

    def test_existing_reader_still_returns_the_numbers(self, tmp_path):
        electron, ar, ar_ion = _species()
        system = _PrescribedTeSystem([[0.0, 1.0, 1e-8, 1.0, 1e-8], [1.0, 2.0, 2e-8, 3.0, 4e-8]])
        path = _write_profile(tmp_path, system, [electron, ar, ar_ion])

        time, data_list = parse_csv_data(path)

        # The whole table, not one column collapsed by zip() truncation.
        assert len(data_list) == 4
        assert np.allclose(time.data, [0.0, 1.0])
        assert np.allclose(data_list[0].data, [1.0, 2.0])
        assert np.allclose(data_list[-1].data, [1e-8, 4e-8])

    def test_provenance_is_recoverable_from_the_file_alone(self, tmp_path):
        electron, ar, ar_ion = _species()
        system = _PrescribedTeSystem([[0.0, 1.0, 1e-8, 1.0, 1e-8]])
        path = _write_profile(tmp_path, system, [electron, ar, ar_ion])

        assert parse_csv_provenance(path) == run_provenance_lines(system)

    def test_a_csv_without_provenance_reports_none(self, tmp_path):
        electron, ar, ar_ion = _species()
        system = _FakeReactionSystem([[0.0, 1.0, 1e-8, 1.0, 1e-8]])
        path = _write_profile(tmp_path, system, [electron, ar, ar_ion])

        assert parse_csv_provenance(path) == []


class ChemkinAndDictionaryLabelTest:
    """Chemkin decks are ASCII, and are re-read by RMG's own parsers."""

    def _model(self):
        electron, ar, ar_ion = _species()
        rxn = Reaction(reactants=[electron, ar], products=[ar_ion, electron, electron],
                       reversible=False,
                       kinetics=TwoTemperaturePlasma(A=(1.0e-3, 'm^3/(mol*s)'), n=0.5,
                                                     Ea_g=(4000.0, 'J/mol'),
                                                     Ea_e=(60000.0, 'J/mol')))
        return [electron, ar, ar_ion], [rxn]

    def test_chemkin_deck_on_disk_carries_the_ascii_disclaimer(self, tmp_path):
        species, reactions = self._model()
        path = os.path.join(str(tmp_path), 'chem.inp')
        save_chemkin_file(path, species, reactions,
                          provenance=comment_block(_plasma_reactor(), '!', ascii_only=True))

        with open(path, encoding='utf-8') as f:
            text = f.read()

        assert PRESCRIBED_TE_MODE in text
        assert PRESCRIBED_TE_DISCLAIMER_ASCII in text
        # A Chemkin deck must stay ASCII: this is why the subscript is spelled n_e.
        text.encode('ascii')
        assert text.startswith('! RMG diagnostic mode:')

    def test_chemkin_deck_without_provenance_is_unchanged(self, tmp_path):
        species, reactions = self._model()
        path = os.path.join(str(tmp_path), 'chem.inp')
        save_chemkin_file(path, species, reactions)

        with open(path, encoding='utf-8') as f:
            assert f.read().startswith('ELEMENTS')

    def test_species_dictionary_label_survives_being_read_back(self, tmp_path):
        from rmgpy.chemkin import load_species_dictionary

        species, _ = self._model()
        path = os.path.join(str(tmp_path), 'species_dictionary.txt')
        save_species_dictionary(path, species,
                                provenance=comment_block(_plasma_reactor(), '//', ascii_only=True))

        with open(path, encoding='utf-8') as f:
            assert PRESCRIBED_TE_DISCLAIMER_ASCII in f.read()

        # The label must not turn the dictionary into a file RMG can no longer load.
        loaded = load_species_dictionary(path, generate_resonance_structures=False)
        assert set(loaded) == {'e-(1)', 'Ar(2)', 'Ar(3)'}

    def test_transport_file_carries_the_label(self, tmp_path):
        from rmgpy.transport import TransportData

        species, _ = self._model()
        for spc in species:
            spc.transport_data = TransportData(shapeIndex=0, sigma=(3.758, 'angstrom'),
                                               epsilon=(148.6, 'K'), dipoleMoment=(0.0, 'C*m'),
                                               polarizability=(0.0, 'angstrom^3'),
                                               rotrelaxcollnum=0.0)
        path = os.path.join(str(tmp_path), 'tran.dat')
        save_transport_file(path, species,
                            provenance=comment_block(_plasma_reactor(), '!', ascii_only=True))

        with open(path, encoding='utf-8') as f:
            text = f.read()
        assert PRESCRIBED_TE_DISCLAIMER_ASCII in text
        text.encode('ascii')


class CanteraLabelTest:
    def test_cantera_yaml_on_disk_carries_the_label_in_and_above_the_mapping(self, tmp_path):
        import yaml

        from rmgpy.rmg.model import ReactionModel
        from rmgpy.yaml_cantera2 import save_cantera_model

        electron, ar, ar_ion = _species()
        rxn = Reaction(reactants=[electron, ar], products=[ar_ion, electron, electron],
                       reversible=False,
                       kinetics=TwoTemperaturePlasma(A=(1.0e-3, 'm^3/(mol*s)'), n=0.5,
                                                     Ea_g=(4000.0, 'J/mol'),
                                                     Ea_e=(60000.0, 'J/mol')))
        model = ReactionModel(species=[electron, ar, ar_ion], reactions=[rxn])
        path = os.path.join(str(tmp_path), 'chem.yaml')
        save_cantera_model(model, path, provenance_source=_plasma_reactor())

        with open(path, encoding='utf-8') as f:
            text = f.read()

        # Above the mapping, for the human who opens the file ...
        assert text.startswith('# RMG diagnostic mode:')
        assert PRESCRIBED_TE_DISCLAIMER in text
        # ... and inside it, for the tool that reparses and drops comments.
        assert PRESCRIBED_TE_DISCLAIMER in yaml.safe_load(text)['description']


class PredictiveElectronDensityRefusalTest:
    """The claim this mode may not make."""

    def test_absolute_electron_density_raises(self):
        with pytest.raises(PrescribedTeDiagnosticError):
            _plasma_reactor().get_electron_density()

    def test_ionisation_degree_raises(self):
        with pytest.raises(PrescribedTeDiagnosticError):
            _plasma_reactor().get_ionisation_degree()

    def test_us_spelling_raises_too(self):
        with pytest.raises(PrescribedTeDiagnosticError):
            _plasma_reactor().get_ionization_degree()

    def test_refusal_teaches_the_distinction(self):
        with pytest.raises(PrescribedTeDiagnosticError) as excinfo:
            _plasma_reactor().get_ionisation_degree()
        message = str(excinfo.value)

        # names the mode ...
        assert PRESCRIBED_TE_MODE in message
        # ... carries the disclaimer verbatim ...
        assert PRESCRIBED_TE_DISCLAIMER in message
        # ... and says what the run MAY claim, so the reader learns rather than loses.
        assert 'may legitimately claim' in message
        for permitted in ('extinction versus sustainment', 'wall flux', 'geometry sensitivity',
                          'particle-balance residual'):
            assert permitted in message

    def test_refusal_is_an_error_not_a_warning(self, recwarn):
        with pytest.raises(PrescribedTeDiagnosticError):
            _plasma_reactor().get_electron_density()
        # A warning that can be ignored is insufficient; nothing here is ignorable.
        assert len(recwarn) == 0


class PermittedPredictionsStillWorkTest:
    """The cases the discriminator deliberately leaves alone."""

    def test_running_to_a_stationary_state_is_still_configurable(self):
        from rmgpy.solver.termination import TerminationSteadyState, TerminationTime

        reactor = _plasma_reactor(termination=[TerminationSteadyState(tolerance=1e-6),
                                               TerminationTime((1.0, 's'))])
        assert any(isinstance(t, TerminationSteadyState) for t in reactor.termination)

    def test_the_steady_state_verdict_is_readable_without_raising(self):
        reactor = _plasma_reactor()
        assert reactor.steady_state_reached is False
        assert np.isnan(reactor.steady_state_residual)

    def test_prescribed_inputs_are_readable_without_raising(self):
        reactor = _plasma_reactor()
        assert reactor.Te.value_si == pytest.approx(T_E)
        assert reactor.T.value_si == pytest.approx(T_GAS)

    def test_the_electron_trajectory_is_readable_from_the_labelled_csv(self, tmp_path):
        electron, ar, ar_ion = _species()
        system = _PrescribedTeSystem([[0.0, 1.0, 1e-8, 1.0, 1e-8], [1.0, 1.0, 5e-8, 1.0, 5e-8]])
        path = _write_profile(tmp_path, system, [electron, ar, ar_ion])

        _, data_list = parse_csv_data(path)
        electron_column = [d for d in data_list if d.label == 'e-(1)'][0]

        # Growth versus decay -- sustainment versus extinction -- is a permitted read.
        assert electron_column.data[-1] > electron_column.data[0]
        # And it arrives with the qualification attached, not stripped.
        assert PRESCRIBED_TE_DISCLAIMER in '\n'.join(parse_csv_provenance(path))

    def test_a_non_plasma_reactor_is_untouched_by_all_of_this(self):
        reactor = _simple_reactor()
        assert getattr(reactor, 'diagnostic_mode', None) is None
        assert not hasattr(reactor, 'get_ionisation_degree')


def _mode_reactor(mode):
    """Real constructors, without an EEDF integration or table load."""
    if mode == PRESCRIBED_TE_MODE:
        return _plasma_reactor()
    electron, ar, ion = _species()
    energy = {'absorbed_power': (0.5, 'W'), 'chamber_volume': (1e-4, 'm^3'),
              'sheath': 'floating_wall', 'electron_energies': {},
              'elastic_collisions': {'Ar': {'ignore': 'constructor-only fixture'}}}
    kwargs = {}
    if mode in (EEDF_MODE, MAXWELLIAN_MODE):
        kwargs['electron_kinetics'] = {
            'provider': 'maxwellian-table' if mode == MAXWELLIAN_MODE else 'loki-table', 'table': ('unused.json', 'a' * 64),
            'branch': 'branch_0', 'initial_reduced_field': (10., 'Td'), 'empirical_laws': {},
        }
        energy.pop('elastic_collisions')
    return PlasmaReactor((300., 'K'), P0, {electron: 1e-6, ion: 1e-6, ar: .999998},
                         (T_E, 'K'), termination=[], diffusion_length=(.01, 'm'),
                         ion_reduced_mobility=(1e-4, 'm^2/(V*s)'),
                         electron_energy_balance=energy, **kwargs)


MODES = [(PRESCRIBED_TE_MODE, PRESCRIBED_TE_DISCLAIMER),
         (SOLVED_TE_MODE, SOLVED_TE_DISCLAIMER), (EEDF_MODE, EEDF_DISCLAIMER),
         (MAXWELLIAN_MODE, MAXWELLIAN_DISCLAIMER)]


@pytest.mark.parametrize('energy_balance,eedf,expected_mode,expected_qualification', [
    pytest.param(False, False, 'LEGACY: prescribed-Te transport diagnostic',
                 'Species evolution conditional on prescribed Te and nₑ; '
                 'absolute electron density and ionisation degree are not predicted.', id='P'),
    pytest.param(True, False, 'LEGACY: solved-Te energy-balance engineering intermediate',
                 'Electron temperature and species evolution are solved conditional on specified '
                 'absorbed power and the declared Maxwellian, wall and sheath closures; coupling '
                 'of an external power source to the plasma is not predicted.', id='S'),
    pytest.param(True, 'loki-table', 'EEDF energy-balance engineering intermediate',
                 'Electron kinetics use the declared EEDF table; reduced field and species '
                 'evolution are solved conditional on specified absorbed power, empirical laws '
                 'and wall/sheath closures. This label does not establish terminal qualification '
                 'or predict external power-source coupling.', id='E'),
    pytest.param(True, 'maxwellian-table',
                 'MODEL-FORM COMPARATOR: Maxwellian-table energy-balance diagnostic',
                 'Electron kinetics use a prescribed Maxwellian on the source-table grid. '
                 'This is a model-form comparator; it is not a production EEDF result and '
                 'does not establish EEDF self-consistency or a physical reduced electric field.', id='M'),
])
def test_configured_reactor_has_literal_mode_and_qualification(
        energy_balance, eedf, expected_mode, expected_qualification):
    # Neither the physical configuration nor this oracle uses production mode
    # constants: swapping their values must not swap the test's expectations.
    electron, ar, ion = _species()
    kwargs = {}
    if energy_balance:
        energy = {'absorbed_power': (0.5, 'W'), 'chamber_volume': (1e-4, 'm^3'),
                  'sheath': 'floating_wall', 'electron_energies': {}}
        if not eedf:
            energy['elastic_collisions'] = {'Ar': {'ignore': 'constructor-only fixture'}}
        kwargs.update(electron_energy_balance=energy, diffusion_length=(.01, 'm'),
                      ion_reduced_mobility=(1e-4, 'm^2/(V*s)'))
    if eedf:
        kwargs['electron_kinetics'] = {
            'provider': eedf, 'table': ('unused.json', 'a' * 64),
            'branch': 'branch_0', 'initial_reduced_field': (10., 'Td'), 'empirical_laws': {},
        }
    reactor = PlasmaReactor((300., 'K'), (1., 'bar'),
                            {electron: 1e-6, ion: 1e-6, ar: .999998}, (10000., 'K'),
                            termination=[], **kwargs)
    assert reactor.energy_balance == energy_balance
    assert reactor.eedf_mode == bool(eedf)
    assert reactor.diagnostic_mode == expected_mode
    expected = ['RMG diagnostic mode: ' + expected_mode, expected_qualification]
    if energy_balance:
        expected.append('Mechanism files exclude reactor-level operators: charged-particle wall loss.')
    assert reactor.get_diagnostic_provenance() == expected
    assert run_provenance_lines(reactor) == expected


@pytest.mark.parametrize('mode,disclaimer', MODES)
def test_real_reactor_classifies_mode_independently_of_te(mode, disclaimer):
    reactor = _mode_reactor(mode)
    expected = ['RMG diagnostic mode: ' + mode, disclaimer]
    if reactor.has_wall:
        expected.append('Mechanism files exclude reactor-level operators: charged-particle wall loss.')
    assert reactor.diagnostic_mode == mode
    assert reactor.get_diagnostic_provenance() == expected
    reactor.Te.value_si *= 2
    assert reactor.get_diagnostic_provenance() == expected
    assert mode in reactor._identity()
    for other_mode, other_disclaimer in MODES:
        if other_mode != mode:
            assert other_disclaimer not in expected


@pytest.mark.parametrize('mode', [SOLVED_TE_MODE, EEDF_MODE, MAXWELLIAN_MODE])
@pytest.mark.parametrize('accessor', ['get_electron_density', 'get_ionisation_degree', 'get_ionization_degree'])
def test_solved_accessors_do_not_claim_prescribed_te(mode, accessor):
    with pytest.raises(NotImplementedError, match=mode):
        getattr(_mode_reactor(mode), accessor)()


@pytest.mark.parametrize('bad', ['', 'unknown', 1, [], False])
def test_unknown_declared_mode_refuses(bad):
    with pytest.raises(ValueError, match='Unknown diagnostic mode'):
        run_provenance_lines(SimpleNamespace(diagnostic_mode=bad))


@pytest.mark.parametrize('mode', [EEDF_MODE, MAXWELLIAN_MODE])
def test_inconsistent_eedf_state_refuses_provenance(mode):
    reactor = _mode_reactor(mode)
    reactor.energy_balance = False
    assert "mode='unconfigured electron energy balance'" in reactor._identity()
    with pytest.raises(PlasmaStateError, match='EEDF provenance requires'):
        run_provenance_lines(reactor)


def test_mixed_modes_are_ordered_and_do_not_leak_to_next_run():
    reactors = [_mode_reactor(mode) for mode, _ in MODES]
    expected = [line for mode, disclaimer in MODES for line in ['RMG diagnostic mode: ' + mode, disclaimer]]
    excluded = 'Mechanism files exclude reactor-level operators: charged-particle wall loss.'
    expected.append(excluded)
    assert run_provenance_lines(reactors + reactors) == expected
    assert run_provenance_lines(_simple_reactor()) == []
    assert run_provenance_lines(reactors[1]) == expected[2:4] + [excluded]


def _seed_job(tmp_path, mode):
    from rmgpy.rmg.main import RMG
    from rmgpy.rmg.model import CoreEdgeReactionModel
    job = RMG(output_directory=str(tmp_path))
    job.reaction_model = CoreEdgeReactionModel()
    job.reaction_systems = [_mode_reactor(mode)]
    job.name = 'provenance'
    species, reactions = ChemkinAndDictionaryLabelTest()._model()
    # Both seed parts must be nonempty: opposite channels have different structures.
    reverse = Reaction(reactants=[species[2], species[0]], products=[species[1]],
                       reversible=False, kinetics=Arrhenius(A=(1., 'm^3/(mol*s)'), n=0, Ea=(0., 'J/mol')))
    job.reaction_model.core.species = list(species)
    job.reaction_model.core.reactions = reactions
    job.reaction_model.edge.species = []
    job.reaction_model.edge.reactions = [reverse]
    job.save_seed_to_database = False
    job.save_seed_modulus = -1
    return job, species


@pytest.mark.parametrize('mode,disclaimer', MODES)
def test_both_labelled_seed_dictionaries_round_trip_by_structure(tmp_path, mode, disclaimer):
    from rmgpy.data.base import Database
    from rmgpy.data.kinetics.library import KineticsLibrary
    from rmgpy.util import generation_token, read_generation_pair
    job, expected = _seed_job(tmp_path, mode)
    input_path = tmp_path / 'input.py'
    input_path.write_text('# original input\n')
    job.input_file = str(input_path)
    job.initialize_seed_mech()
    job.make_seed_mech()
    assert input_path.read_text() == '# original input\n'
    assert disclaimer in (tmp_path / 'restart_from_seed.py').read_text()
    for name in ('seed', 'seed_edge'):
        directory = tmp_path / 'seed' / name
        dictionary = directory / 'dictionary.txt'
        reactions = directory / 'reactions.py'
        assert disclaimer in dictionary.read_text()
        assert mode in reactions.read_text()
        assert generation_token(dictionary.read_text()) == generation_token(reactions.read_text())
        read_generation_pair(str(reactions), str(dictionary))
        loaded = Database().get_species(str(dictionary))
        assert len(loaded) == len(expected) > 0
        for species in expected:
            assert any(species.is_isomorphic(candidate) for candidate in loaded.values()), species
        assert not expected[0].is_isomorphic(expected[1])
        library = KineticsLibrary()
        from rmgpy.data.kinetics.database import KineticsDatabase
        library.load(str(reactions), local_context=dict(KineticsDatabase().local_context), global_context={})
        assert len(library.entries) == 1
    import h5py
    import yaml
    with h5py.File(tmp_path / 'seed/filters/filters.h5') as stream:
        assert disclaimer in stream.attrs['run_provenance']
    mapping = tmp_path / 'seed/filters/species_map.yml'
    assert disclaimer in mapping.read_text()
    assert len(yaml.safe_load(mapping.read_text())) == len(expected)


def test_failed_labelled_seed_pair_save_restores_both_files(tmp_path, monkeypatch):
    from rmgpy.data.base import Entry
    from rmgpy.data.kinetics.library import KineticsLibrary
    from rmgpy import util
    _, reactions = ChemkinAndDictionaryLabelTest()._model()
    library = KineticsLibrary(name='seed')
    library.entries = {1: Entry(index=1, item=reactions[0], data=reactions[0].kinetics)}
    path = tmp_path / 'reactions.py'
    library.save(str(path), provenance_source=_plasma_reactor())
    before = {p: p.read_bytes() for p in (path, tmp_path / 'dictionary.txt')}
    rewrite = util._rewrite_save
    calls = []
    def fail_second(record, rendered):
        calls.append(record)
        if len(calls) == 2:
            raise OSError('injected dictionary write failure')
        return rewrite(record, rendered)
    monkeypatch.setattr(util, '_rewrite_save', fail_second)
    with pytest.raises(OSError, match='injected dictionary'):
        library.save(str(path), provenance_source=_mode_reactor(SOLVED_TE_MODE))
    assert len(calls) == 2
    assert all(path.read_bytes() == content for path, content in before.items())


@pytest.mark.parametrize('mode,disclaimer', MODES)
def test_profile_plot_stats_html_carry_the_matching_mode(tmp_path, mode, disclaimer):
    from PIL import Image
    import xlrd
    from rmgpy.rmg.output import save_output
    from rmgpy.stats import ExecutionStatsWriter
    from rmgpy.tools.plot import SimulationPlot
    job, species = _seed_job(tmp_path, mode)
    reactor = job.reaction_systems[0]
    reactor.snapshots = [[0., 1., .01, .98, .01], [1., 1., .02, .96, .02]]
    # EEDF manifests are tested with a fully initialized provider below.
    if mode in (EEDF_MODE, MAXWELLIAN_MODE):
        proxy = SimpleNamespace(diagnostic_mode=mode, snapshots=reactor.snapshots)
    else:
        proxy = reactor
    profile = _write_profile(tmp_path, proxy, species)
    time, columns = parse_csv_data(profile)
    np.testing.assert_array_equal(time.data, [0., 1.])
    assert len(columns) == 4
    image_path = tmp_path / 'simulation.png'
    SimulationPlot(csv_file=profile, title='existing title').plot(str(image_path))
    with Image.open(image_path) as image:
        assert image.info['Title'] == 'existing title'
        assert disclaimer in image.info['Description']
    job.exec_time = [1.]
    writer = ExecutionStatsWriter(str(tmp_path))
    writer.coreSpeciesCount = [3]
    writer.coreReactionCount = [1]
    writer.edgeSpeciesCount = [0]
    writer.edgeReactionCount = [1]
    writer.memoryUse = [1.]
    writer.save_execution_statistics(job)
    writer.generate_execution_plots(job)
    book = xlrd.open_workbook(tmp_path / 'statistics.xls')
    assert book.sheet_by_name('Diagnostic mode').cell_value(1, 0) == disclaimer
    assert book.sheet_by_name('Statistics').cell_value(1, 1) == 3
    for name in ('coreSize', 'edgeSize', 'memoryUse'):
        assert disclaimer in (tmp_path / 'plot' / (name + '.svg')).read_text()
    # The HTML notice applies even to a species-only mechanism.
    job.reaction_model.core.reactions = []
    job.reaction_model.edge.reactions = []
    save_output(job, save_edge=True)
    for name in ('output.html', 'output_edge.html'):
        assert disclaimer in (tmp_path / name).read_text()
    # Scope exemption: structure thumbnails are produced normally, without run metadata.
    for path in (tmp_path / 'species').glob('*.png'):
        with Image.open(path) as image:
            assert 'RMG diagnostic mode:' not in repr(image.info)


def test_real_eedf_json_preserves_provider_manifest(tmp_path):
    import runpy
    build_reactor = runpy.run_path(str(Path(__file__).parent / 'solver/plasmaEEDFCoreTest.py'))['build_reactor']
    reactor, core, _ = build_reactor(tmp_path)
    reactor.snapshots = [[0., 1.] + [0.] * len(core)]
    before = reactor.eedf_run_manifest()
    path = Path(_write_profile(tmp_path, reactor, core)).with_suffix('.eedf.json')
    record = json.loads(path.read_text())
    assert record['run_provenance'] == run_provenance_lines(reactor)
    assert record == before
    assert 'artifact_sha256' in record


@pytest.mark.parametrize('mode,disclaimer', MODES)
def test_optional_yaml_and_comparison_writers(tmp_path, mode, disclaimer):
    import yaml
    from rmgpy.yaml_cantera1 import write_cantera
    from rmgpy.yaml_rms import write_rms
    from rmgpy.tools.compare_cantera_yaml import compare_yaml_files_and_report
    reactor = _mode_reactor(mode)
    # Direct renderers of an admissible neutral mechanism receive the enclosing run's provenance.
    ar = _species()[1]
    from rmgpy.transport import TransportData
    ar.transport_data = TransportData(shapeIndex=0, epsilon=(100., 'K'), sigma=(3.4, 'angstroms'),
                                      dipoleMoment=(0., 'C*m'), polarizability=(0., 'angstroms^3'),
                                      rotrelaxcollnum=0.)
    from rmgpy.rmg.model import ReactionModel
    model = ReactionModel(species=[ar], reactions=[])
    cantera = tmp_path / 'cantera.yaml'
    write_cantera([ar], [], model.get_elements(), path=str(cantera), provenance_source=reactor)
    rms = tmp_path / 'chem.rms'
    write_rms([ar], [], path=str(rms), provenance_source=reactor)
    for path in (cantera, rms):
        assert disclaimer in path.read_text()
        assert yaml.safe_load(path.read_text())
    report = tmp_path / 'comparison.txt'
    assert compare_yaml_files_and_report(str(cantera), str(cantera), output=str(report),
                                         provenance=comment_block(reactor, '#'))
    assert disclaimer in report.read_text()


def test_wall_json_adds_provenance_without_mutating_closure_record(tmp_path):
    mode, disclaimer = MODES[1]
    manifest = {'closure': 'fixture', 'geometry_arm': 'cylinder', 'scientific_status': 'conditional',
                'gates': None, 'floating_potential_e_over_kTe': None, 'time': 0.}
    wall_record = {'h': 1., 'floating_potential_e_over_kTe': None, 'gates': None}
    reactor = SimpleNamespace(diagnostic_mode=mode, snapshots=[[0., 1., .1, .8, .1]],
                              electronegative_wall_model='fixture',
                              electronegative_wall_manifest=lambda: manifest,
                              electronegative_wall_history={0.: wall_record})
    profile = Path(_write_profile(tmp_path, reactor, _species()))
    record = json.loads(profile.with_suffix('.electronegative-wall.json').read_text())
    assert record.pop('run_provenance') == ['RMG diagnostic mode: ' + mode, disclaimer]
    assert record == manifest
    assert 'run_provenance' not in manifest
    from rmgpy.rmg.listener import SimulationProfilePlotter
    SimulationProfilePlotter(str(tmp_path), 0, _species()).update(reactor)
    from PIL import Image
    with Image.open(profile.with_suffix('.png')) as image:
        assert mode in image.info['Description']


def test_empty_edge_dictionary_keeps_label_and_loads_empty(tmp_path):
    from rmgpy.data.base import Database
    from rmgpy.data.kinetics.library import KineticsLibrary
    library = KineticsLibrary(name='seed_edge')
    library.save(str(tmp_path / 'reactions.py'), provenance_source=_plasma_reactor())
    dictionary = tmp_path / 'dictionary.txt'
    assert PRESCRIBED_TE_DISCLAIMER in dictionary.read_text()
    assert Database().get_species(str(dictionary)) == {}


@pytest.mark.parametrize('mode,disclaimer', MODES[:2])
def test_converted_yaml_preserves_chemkin_provenance(tmp_path, mode, disclaimer):
    import yaml
    from rmgpy.chemkin import save_chemkin
    from rmgpy.rmg.main import RMG
    from rmgpy.rmg.model import CoreEdgeReactionModel
    job = RMG(output_directory=str(tmp_path))
    job.reaction_systems = [_mode_reactor(mode)]
    job.reaction_model = CoreEdgeReactionModel()
    job.reaction_model.core.species = [_species()[1]]
    from rmgpy.transport import TransportData
    job.reaction_model.core.species[0].transport_data = TransportData(
        shapeIndex=0, epsilon=(100., 'K'), sigma=(3.4, 'angstrom'),
        dipoleMoment=(0., 'C*m'), polarizability=(0., 'angstrom^3'), rotrelaxcollnum=0.)
    directory = tmp_path / 'chemkin'
    directory.mkdir()
    save_chemkin(job.reaction_model, str(directory / 'chem.inp'), str(directory / 'chem_annotated.inp'),
                 dictionary_path=str(directory / 'species_dictionary.txt'),
                 transport_path=str(directory / 'tran.dat'), provenance_source=job)
    path = job.generate_cantera_files_from_chemkin(str(directory / 'chem.inp'))
    data = yaml.safe_load(Path(path).read_text())
    assert mode in data['description']
    assert disclaimer.replace('nₑ', 'n_e') in data['description']


@pytest.fixture(autouse=True)
def comparator_outputs_never_claim_production(request, tmp_path):
    """Inspect the artifacts written by each existing writer's M case."""
    yield
    params = getattr(getattr(request.node, 'callspec', None), 'params', {})
    if params.get('mode') != 'MODEL-FORM COMPARATOR: Maxwellian-table energy-balance diagnostic':
        return
    for path in tmp_path.rglob('*'):
        if not path.is_file():
            continue
        if path.suffix in ('.csv', '.json', '.py', '.txt', '.inp', '.yaml', '.yml', '.rms', '.svg', '.html'):
            assert 'EEDF energy-balance engineering intermediate' not in path.read_text()
        elif path.suffix == '.png' and path.parent.name != 'species':
            from PIL import Image
            with Image.open(path) as picture:
                assert 'EEDF energy-balance engineering intermediate' not in repr(picture.info)


def test_mechanism_writers_name_excluded_reactor_operators(tmp_path):
    from rmgpy.yaml_cantera1 import write_cantera
    from rmgpy.yaml_cantera2 import save_cantera_model
    from rmgpy.yaml_rms import write_rms
    from rmgpy.rmg.model import ReactionModel
    from rmgpy.transport import TransportData
    system = _plasma_reactor()
    # Public declarations also represent jobs with wall loss alone or trapping alone.
    source = SimpleNamespace(diagnostic_mode=system.diagnostic_mode, has_wall=True,
                             wall_neutral_diffusion={'Ars': {}},
                             resonance_radiation={'model':'complete','geometryModel':'slab'})
    ar = _species()[1]
    ar.transport_data = TransportData(shapeIndex=0, epsilon=(100.,'K'), sigma=(3.4,'angstroms'),
        dipoleMoment=(0.,'C*m'), polarizability=(0.,'angstroms^3'), rotrelaxcollnum=0.)
    model = ReactionModel(species=[ar],reactions=[])
    ck, ct1, ct2, rms = [tmp_path/name for name in ('chem.inp','ct1.yaml','ct2.yaml','chem.rms')]
    save_chemkin_file(str(ck), [ar], [], provenance=comment_block(source,'!',ascii_only=True))
    write_cantera([ar], [], model.get_elements(), path=str(ct1), provenance_source=source)
    save_cantera_model(model, str(ct2), provenance_source=source)
    write_rms([ar], [], path=str(rms), provenance_source=source)
    for path in (ck,ct1,ct2,rms):
        text = path.read_text()
        assert 'Mechanism files exclude reactor-level operators:' in text
        assert 'charged-particle wall loss' in text
        assert 'neutral wall diffusion' in text
        assert 'Ar resonance radiation (complete)' in text
    assert len(model.reactions) == 0


def test_radiation_profile_writes_accepted_manifest(tmp_path):
    manifest = {'accepted': {'time_s':1.,'photon_power_W':2.}, 'declaration':{'model':'holsteinWalsh'}}
    system = SimpleNamespace(diagnostic_mode=PRESCRIBED_TE_MODE,
        resonance_radiation=manifest['declaration'], radiation_manifest=lambda: manifest,
        snapshots=[[1.,1.,1e-8,1.,1e-8]], radiation_history={1.:manifest['accepted']})
    path = Path(_write_profile(tmp_path, system, _species())).with_suffix('.radiation.json')
    saved = json.loads(path.read_text())
    assert saved['accepted'] == manifest['accepted']
    assert saved['trajectory'] == [manifest['accepted']]
    assert 'run_provenance' in saved
    assert 'trajectory' not in manifest
