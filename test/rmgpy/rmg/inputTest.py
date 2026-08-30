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

from unittest.mock import patch

import numpy as np

import rmgpy.constants as constants
import rmgpy.rmg.input as inp
from rmgpy.exceptions import InputError
from rmgpy.rmg.input import _parse_writer_config, _writer_config_to_input
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.rmg.settings import WriterConfig
from rmgpy.ml.estimator import ADMONITION

import pytest


def setup_module(module):
    """
    A method that is run before the class.
    """
    # set-up RMG object and get global rmg object in input.py file
    # so methods can be tested
    global rmg
    rmg = RMG()
    inp.set_global_rmg(rmg)


def teardown_module(module):
    # remove RMG object
    global rmg
    rmg = None


class TestInputDatabase:
    """
    Contains unit tests rmgpy.rmg.input.database
    """

    def teardown_class(self):
        # remove the reactionLibraries value
        global rmg
        rmg.reaction_libraries = None

    def test_importing_database_reaction_libraries_from_string(self):
        """
        Test that reaction libraries given as plain strings are stored as strings.
        """
        global rmg
        inp.database(reactionLibraries=["test"])
        assert rmg.reaction_libraries == ["test"]
        assert "test" not in rmg.reaction_libraries_output_edge

    def test_importing_database_reaction_libraries_from_false_tuple(self):
        """
        Test that (name, False) tuples are stored as plain strings without output_edge.
        """
        global rmg
        inp.database(reactionLibraries=[("test", False)])
        assert rmg.reaction_libraries == ["test"]
        assert "test" not in rmg.reaction_libraries_output_edge

    def test_importing_database_reaction_libraries_from_true_tuple(self):
        """
        Test that (name, True) tuples are stored as plain strings with the name in output_edge.
        """
        global rmg
        inp.database(reactionLibraries=[("test", True)])
        assert rmg.reaction_libraries == ["test"]
        assert "test" in rmg.reaction_libraries_output_edge

class TestInputDatabaseAutoSelection:
    """Tests for 'auto' and '<PAH_libs>' token handling in database()."""

    def test_thermo_auto_string(self):
        global rmg
        inp.database(thermoLibraries='auto')
        assert rmg.thermo_libraries == 'auto'

    def test_thermo_auto_in_list(self):
        global rmg
        inp.database(thermoLibraries=['myLib', 'auto'])
        assert rmg.thermo_libraries == ['myLib', 'auto']

    def test_thermo_pah_libs_in_list(self):
        global rmg
        inp.database(thermoLibraries=['auto', '<PAH_libs>'])
        assert rmg.thermo_libraries == ['auto', '<PAH_libs>']

    def test_transport_auto_string(self):
        global rmg
        inp.database(transportLibraries='auto')
        assert rmg.transport_libraries == 'auto'

    def test_reaction_libraries_auto_string(self):
        global rmg
        inp.database(reactionLibraries='auto')
        assert rmg.reaction_libraries == 'auto'

    def test_reaction_libraries_auto_in_list(self):
        global rmg
        inp.database(reactionLibraries=['auto', '<PAH_libs>'])
        assert 'auto' in rmg.reaction_libraries
        assert '<PAH_libs>' in rmg.reaction_libraries

    def test_reaction_libraries_mixed_auto_and_tuple(self):
        global rmg
        inp.database(reactionLibraries=[('myLib', True), 'auto'])
        assert rmg.reaction_libraries == ['myLib', 'auto']
        assert 'myLib' in rmg.reaction_libraries_output_edge

    def test_seed_mechanisms_auto(self):
        global rmg
        inp.database(seedMechanisms='auto')
        assert rmg.seed_mechanisms == 'auto'

    def test_kinetics_families_auto(self):
        global rmg
        inp.database(kineticsFamilies='auto')
        assert rmg.kinetics_families == 'auto'

    def test_kinetics_families_auto_with_exclusion(self):
        global rmg
        inp.database(kineticsFamilies=['!H_Abstraction', 'auto'])
        assert rmg.kinetics_families == ['!H_Abstraction', 'auto']

    def test_default_no_auto(self):
        """Without 'auto', fields should behave as before."""
        global rmg
        inp.database(
            thermoLibraries=['primaryThermoLibrary'],
            reactionLibraries=[('lib1', False)],
            seedMechanisms=[],
            transportLibraries=None,
            kineticsFamilies='default',
        )
        assert rmg.thermo_libraries == ['primaryThermoLibrary']
        assert rmg.reaction_libraries == ['lib1']
        assert 'lib1' not in rmg.reaction_libraries_output_edge
        assert rmg.seed_mechanisms == []
        assert rmg.transport_libraries is None
        assert rmg.kinetics_families == 'default'

    def test_pah_libs_standalone_thermo_raises(self):
        with pytest.raises(InputError):
            inp.database(thermoLibraries='<PAH_libs>')

    def test_pah_libs_standalone_reaction_raises(self):
        with pytest.raises(InputError):
            inp.database(reactionLibraries='<PAH_libs>')

    def test_pah_libs_standalone_transport_raises(self):
        with pytest.raises(InputError):
            inp.database(transportLibraries='<PAH_libs>')

    def test_pah_libs_standalone_seeds_raises(self):
        with pytest.raises(InputError):
            inp.database(seedMechanisms='<PAH_libs>')

    def test_pah_libs_tuple_in_reaction_libs_raises(self):
        with pytest.raises(InputError):
            inp.database(reactionLibraries=[('<PAH_libs>', True)])


@pytest.mark.skip(reason=ADMONITION)
class TestInputMLEstimator:
    """
    Contains unit tests rmgpy.rmg.input.mlEstimator
    """

    def teardown_class(self):
        # remove the reactionLibraries value
        global rmg
        rmg.ml_estimator = None

    def test_ml_estimator(self):
        """
        Test that we can input.
        """
        from rmgpy.ml.estimator import MLEstimator

        global rmg
        # add database properties to RMG
        inp.ml_estimator(thermo=True)
        assert isinstance(rmg.ml_estimator, MLEstimator)
        assert isinstance(rmg.ml_settings, dict)


class TestInputThermoCentralDatabase:
    """
    Contains unit tests rmgpy.rmg.input.thermo_central_database
    """

    def teardown_class(self):
        # remove the reactionLibraries value
        global rmg
        rmg.thermo_central_database = None

    def test_thermo_central_database(self):
        """
        Test that we can input.
        """
        global rmg
        # add database properties to RMG
        inp.thermo_central_database(
            host="some_host",
            port=0,
            username="some_usr",
            password="some_pw",
            application="some_app",
        )
        assert rmg.thermo_central_database.host == "some_host"
        assert rmg.thermo_central_database.port == 0
        assert rmg.thermo_central_database.username == "some_usr"
        assert rmg.thermo_central_database.password == "some_pw"
        assert rmg.thermo_central_database.application == "some_app"
        assert rmg.thermo_central_database.client == None


class TestInputReactors:
    """
    Contains unit tests for reactor input classes
    """

    @pytest.fixture(autouse=True)
    def setup_rmg(self):
        """This method is run before every test in this class"""
        # Create a mock species dictionary
        # In reality, the values would be Species objects, but it doesn't matter for testing
        species_dict = {
            "A": "A",
            "B": "B",
            "C": "C",
            "X": "X",
        }

        # Assign to global variable in the input module
        inp.species_dict = species_dict

        # Initialize the rmg.reaction_systems attribute
        global rmg
        rmg.reaction_systems = []
        yield

    def teardown_class(self):
        """This method is run after every test in this class"""
        # Reset the global species_dict variable in the input module
        inp.species_dict = {}

        # Reset the rmg.reaction_systems attribute
        global rmg
        rmg.reaction_systems = []

    def test_simple_reactor_mole_fractions(self):
        """Test that SimpleReactor mole fractions are set properly"""
        inp.simple_reactor(
            temperature=(1000, "K"),
            pressure=(1, "atm"),
            initialMoleFractions={
                "A": 0.5,
                "B": 0.3,
                "C": 0.2,
            },
            terminationTime=(1, "s"),
        )

        global rmg
        reactor = rmg.reaction_systems[0]
        assert reactor.initial_mole_fractions["A"] == 0.5
        assert reactor.initial_mole_fractions["B"] == 0.3
        assert reactor.initial_mole_fractions["C"] == 0.2

    @patch("rmgpy.rmg.input.logging")
    def test_simple_reactor_mole_fractions_normalize_1(self, mock_logging):
        """Test that SimpleReactor mole fractions are normalized properly"""
        inp.simple_reactor(
            temperature=(1000, "K"),
            pressure=(1, "atm"),
            initialMoleFractions={
                "A": 5,
                "B": 3,
                "C": 2,
            },
            terminationTime=(1, "s"),
        )

        global rmg
        reactor = rmg.reaction_systems[0]
        assert reactor.initial_mole_fractions["A"] == 0.5
        assert reactor.initial_mole_fractions["B"] == 0.3
        assert reactor.initial_mole_fractions["C"] == 0.2

        mock_logging.warning.assert_called_with("Initial mole fractions do not sum to one; normalizing.")

    @patch("rmgpy.rmg.input.logging")
    def test_simple_reactor_mole_fractions_normalize_2(self, mock_logging):
        """Test that SimpleReactor mole fractions are normalized properly"""
        inp.simple_reactor(
            temperature=[(1000, "K"), (2000, "K")],
            pressure=[(1, "atm"), (10, "atm")],
            initialMoleFractions={
                "A": 5,
                "B": 3,
                "C": 2,
            },
            terminationTime=(1, "s"),
        )

        global rmg
        reactor = rmg.reaction_systems[0]
        assert reactor.initial_mole_fractions["A"] == 0.5
        assert reactor.initial_mole_fractions["B"] == 0.3
        assert reactor.initial_mole_fractions["C"] == 0.2

        mock_logging.warning.assert_called_with("Initial mole fractions do not sum to one; normalizing.")

    def test_simple_reactor_mole_fractions_ranged(self):
        """Test that SimpleReactor ranged mole fractions are not normalized"""
        inp.simple_reactor(
            temperature=[(1000, "K"), (2000, "K")],
            pressure=[(1, "atm"), (10, "atm")],
            initialMoleFractions={
                "A": [5, 8],
                "B": 3,
                "C": 2,
            },
            terminationTime=(1, "s"),
        )

        global rmg
        reactor = rmg.reaction_systems[0]
        assert reactor.initial_mole_fractions["A"] == [5, 8]
        assert reactor.initial_mole_fractions["B"] == 3
        assert reactor.initial_mole_fractions["C"] == 2

    def test_liquid_reactor_concentrations(self):
        """Test that LiquidReactor concentrations are set properly"""
        inp.liquid_reactor(
            temperature=(1000, "K"),
            initialConcentrations={
                "A": (0.3, "mol/L"),
                "B": (0.2, "mol/L"),
                "C": (0.1, "mol/L"),
            },
            terminationTime=(1, "s"),
        )

        global rmg
        reactor = rmg.reaction_systems[0]

        # Values get converted to default SI units, mol/m^3
        assert reactor.initial_concentrations["A"] == 300
        assert reactor.initial_concentrations["B"] == 200
        assert reactor.initial_concentrations["C"] == 100

    def test_surface_reactor_mole_fractions(self):
        """Test that SurfaceReactor mole fractions are set properly"""
        inp.surface_reactor(
            temperature=(1000, "K"),
            initialPressure=(1, "atm"),
            initialGasMoleFractions={
                "A": 0.5,
                "B": 0.3,
                "C": 0.2,
            },
            initialSurfaceCoverages={"X": 1.0},
            surfaceVolumeRatio=(1e1, "m^-1"),
            terminationTime=(1, "s"),
        )

        global rmg
        reactor = rmg.reaction_systems[0]
        assert reactor.initial_gas_mole_fractions["A"] == 0.5
        assert reactor.initial_gas_mole_fractions["B"] == 0.3
        assert reactor.initial_gas_mole_fractions["C"] == 0.2

    @patch("rmgpy.rmg.input.logging")
    def test_surface_reactor_mole_fractions_normalize_1(self, mock_logging):
        """Test that SurfaceReactor mole fractions are normalized properly"""
        inp.surface_reactor(
            temperature=(1000, "K"),
            initialPressure=(1, "atm"),
            initialGasMoleFractions={
                "A": 5,
                "B": 3,
                "C": 2,
            },
            initialSurfaceCoverages={"X": 1.0},
            surfaceVolumeRatio=(1e1, "m^-1"),
            terminationTime=(1, "s"),
        )

        global rmg
        reactor = rmg.reaction_systems[0]
        assert reactor.initial_gas_mole_fractions["A"] == 0.5
        assert reactor.initial_gas_mole_fractions["B"] == 0.3
        assert reactor.initial_gas_mole_fractions["C"] == 0.2

        mock_logging.warning.assert_called_with("Initial gas mole fractions do not sum to one; renormalizing.")

    def test_mb_sampled_reactor_mole_fractions(self):
        """Test that MBSampledReactor mole fractions are set properly"""
        inp.mb_sampled_reactor(
            temperature=(1000, "K"),
            pressure=(1, "atm"),
            initialMoleFractions={
                "A": 0.5,
                "B": 0.3,
                "C": 0.2,
            },
            mbsamplingRate=3500,
            terminationTime=(1, "s"),
            constantSpecies=["B", "C"],
        )

        global rmg
        reactor = rmg.reaction_systems[0]
        assert reactor.initial_mole_fractions["A"] == 0.5
        assert reactor.initial_mole_fractions["B"] == 0.3
        assert reactor.initial_mole_fractions["C"] == 0.2

    @patch("rmgpy.rmg.input.logging")
    def test_mb_sampled_reactor_mole_fractions_normalize_1(self, mock_logging):
        """Test that MBSampledReactor mole fractions are normalized properly"""
        inp.mb_sampled_reactor(
            temperature=(1000, "K"),
            pressure=(1, "atm"),
            initialMoleFractions={
                "A": 5,
                "B": 3,
                "C": 2,
            },
            mbsamplingRate=3500,
            terminationTime=(1, "s"),
            constantSpecies=["B", "C"],
        )

        global rmg
        reactor = rmg.reaction_systems[0]
        assert reactor.initial_mole_fractions["A"] == 0.5
        assert reactor.initial_mole_fractions["B"] == 0.3
        assert reactor.initial_mole_fractions["C"] == 0.2

        mock_logging.warning.assert_called_with("Initial mole fractions do not sum to one; normalizing.")


class TestInputPressureDependence:
    """
    Contains unit tests for pressure dependence input, including completedNetworks
    """

    def setup_method(self):
        """This method is run before every test in this class"""
        global rmg
        # Reset the completed networks set before each test
        rmg.reaction_model = CoreEdgeReactionModel()
        rmg.reaction_model.completed_pdep_networks = set()

    def test_completed_networks_single(self):
        """Test that a single completedNetwork can be added via pressure_dependence"""
        global rmg
        
        inp.pressure_dependence(
            method='modified strong collision',
            temperatures=(300, 2000, 'K', 8),
            pressures=(0.01, 100, 'bar', 5),
            maximumGrainSize=(0.5, 'kcal/mol'),
            minimumNumberOfGrains=250,
            interpolation=('Chebyshev', 6, 4),
            maximumAtoms=16,
            completedNetworks=['CH2O2'],
        )
        
        # Check that the network was added
        assert len(rmg.reaction_model.completed_pdep_networks) == 1
        # The formula CH2O2 should be converted to a sorted tuple of elements
        expected_key = (('C', 1), ('H', 2), ('O', 2))
        assert expected_key in rmg.reaction_model.completed_pdep_networks

    def test_completed_networks_multiple(self):
        """Test that multiple completedNetworks can be added via pressure_dependence"""
        global rmg
        
        inp.pressure_dependence(
            method='modified strong collision',
            temperatures=(300, 2000, 'K', 8),
            pressures=(0.01, 100, 'bar', 5),
            maximumGrainSize=(0.5, 'kcal/mol'),
            minimumNumberOfGrains=250,
            interpolation=('Chebyshev', 6, 4),
            maximumAtoms=16,
            completedNetworks=['CH2O2', 'C2H6'],
        )
        
        # Check that both networks were added
        assert len(rmg.reaction_model.completed_pdep_networks) == 2
        expected_key1 = (('C', 1), ('H', 2), ('O', 2))
        expected_key2 = (('C', 2), ('H', 6))
        assert expected_key1 in rmg.reaction_model.completed_pdep_networks
        assert expected_key2 in rmg.reaction_model.completed_pdep_networks

    def test_completed_networks_none(self):
        """Test that pressure_dependence works without completedNetworks"""
        global rmg
        
        inp.pressure_dependence(
            method='modified strong collision',
            temperatures=(300, 2000, 'K', 8),
            pressures=(0.01, 100, 'bar', 5),
            maximumGrainSize=(0.5, 'kcal/mol'),
            minimumNumberOfGrains=250,
            interpolation=('Chebyshev', 6, 4),
            maximumAtoms=16,
        )
        
        # Check that no networks were added
        assert len(rmg.reaction_model.completed_pdep_networks) == 0


class TestWriterConfig:
    """Unit tests for WriterConfig and the _parse_writer_config / _writer_config_to_input helpers."""

    def test_parse_false_disables(self):
        cfg = _parse_writer_config(False)
        assert not cfg.enabled
        assert cfg.save_interval == 0

    def test_parse_true_enables_default_interval(self):
        cfg = _parse_writer_config(True)
        assert cfg.enabled
        assert cfg.save_interval == 1
        assert cfg.save_edge is None

    def test_parse_true_custom_default_interval(self):
        cfg = _parse_writer_config(True, default_save_interval=5)
        assert cfg.save_interval == 5

    def test_parse_dict_full(self):
        cfg = _parse_writer_config({'saveInterval': -1, 'saveEdge': False})
        assert cfg.enabled
        assert cfg.save_interval == -1
        assert cfg.save_edge is False

    def test_parse_dict_partial(self):
        cfg = _parse_writer_config({'saveInterval': 3})
        assert cfg.save_interval == 3
        assert cfg.save_edge is None

    def test_parse_invalid_raises(self):
        with pytest.raises(InputError):
            _parse_writer_config(42)

    def test_should_write_every_iteration(self):
        cfg = WriterConfig(save_interval=1)
        assert cfg.should_write(0, False)
        assert cfg.should_write(1, False)
        assert cfg.should_write(5, False)

    def test_should_write_every_n_iterations(self):
        cfg = WriterConfig(save_interval=3)
        assert cfg.should_write(0, False)
        assert not cfg.should_write(1, False)
        assert not cfg.should_write(2, False)
        assert cfg.should_write(3, False)
        assert cfg.should_write(6, False)

    def test_should_write_end_only(self):
        cfg = WriterConfig(save_interval=-1)
        assert not cfg.should_write(0, False)
        assert not cfg.should_write(5, False)
        assert cfg.should_write(5, True)

    def test_should_write_disabled(self):
        cfg = WriterConfig(save_interval=0)
        assert not cfg.should_write(0, False)
        assert not cfg.should_write(0, True)

    def test_should_write_final_always_writes(self):
        cfg = WriterConfig(save_interval=5)
        # Iteration 7 was not written (not multiple of 5)
        assert cfg.should_write(7, True)

    def test_should_write_final_no_double_write(self):
        cfg = WriterConfig(save_interval=5)
        # Iteration 5 is a multiple of 5 — written during loop
        cfg.should_write(5, False)
        # Final call at same iteration should be skipped
        assert not cfg.should_write(5, True)

    def test_writer_config_to_input_false(self):
        cfg = WriterConfig(save_interval=0)
        assert _writer_config_to_input(cfg) is False

    def test_writer_config_to_input_true(self):
        cfg = WriterConfig(save_interval=1)
        assert _writer_config_to_input(cfg) is True

    def test_writer_config_to_input_dict(self):
        cfg = WriterConfig(save_interval=-1, save_edge=False)
        result = _writer_config_to_input(cfg)
        assert "'saveInterval': -1" in result
        assert "'saveEdge': False" in result

    def test_writer_config_to_input_none(self):
        assert _writer_config_to_input(None) is False


class TestInputPlasmaReactor:
    """
    Functional tests for the ``plasmaReactor`` input directive and the electron
    number-density conversion (I-067 increment 1).

    Every test drives the *public* :func:`rmgpy.rmg.input.read_input_file` so that
    the ``convert_initial_keys_to_species_objects`` hook at input.py:1737-1739 is
    genuinely exercised. None of these tests constructs ``plasma_reactor`` or
    ``PlasmaReactor`` directly: doing so would bypass that hook and hide whether a
    user can actually run a plasma simulation from an input file.
    """

    T = 1000.0        # K, gas temperature
    TE = 11604.5      # K, electron temperature (~1 eV)
    P_ATM_SI = 101325.0
    N_E = 1.0e23      # m^-3, requested electron number density

    # ---- helpers -----------------------------------------------------------

    def _read(self, tmp_path, body):
        """Write ``body`` as an input file and parse it through read_input_file."""
        path = tmp_path / "input.py"
        path.write_text(body)
        rmg = RMG()
        inp.read_input_file(str(path), rmg)
        return rmg

    @staticmethod
    def _e_species():
        return "species(label='e-', structure=adjacencyList(\"1 e u1 p0 c-1\"))\n"

    @staticmethod
    def _heavy_species():
        return (
            "species(label='Ar', structure=adjacencyList(\"1 Ar u0 p4 c0\"))\n"
            "species(label='He', structure=adjacencyList(\"1 He u0 p1 c0\"))\n"
        )

    def _preamble(self):
        return self._e_species() + self._heavy_species()

    def _expected_conversion(self, n_e, T, Te, P_si):
        # The builder uses the reactor's own kB = R/Na (not constants.kB) so the EOS it
        # inverts and the conversion agree by construction; the oracle must match.
        boltzmann = constants.R / constants.Na
        q = n_e * boltzmann * Te / P_si
        a = q * T / Te
        x_e = a / (1.0 - q + a)
        heavy_scale = (1.0 - q) / (1.0 - q + a)
        return q, x_e, heavy_scale

    def _plasma_block(self, extra, mole_fractions="{'Ar': 0.7, 'He': 0.3}",
                      temperature="(1000,'K')", pressure="(1,'atm')",
                      electron_temperature="(11604.5,'K')"):
        return (
            "plasmaReactor(\n"
            f"    temperature={temperature},\n"
            f"    pressure={pressure},\n"
            f"    electronTemperature={electron_temperature},\n"
            f"{extra}"
            f"    initialMoleFractions={mole_fractions},\n"
            "    terminationTime=(1,'s'),\n"
            ")\n"
        )

    # ---- positive path -----------------------------------------------------

    def test_electron_density_conversion(self, tmp_path):
        """Assertion 1 & 2: electron at x_e, heavy at input*heavy_scale, sum == 1."""
        body = self._preamble() + self._plasma_block("    electronDensity=(1e23,'m^-3'),\n")
        rmg = self._read(tmp_path, body)
        reactor = rmg.reaction_systems[0]
        imf = {k.label: v for k, v in reactor.initial_mole_fractions.items()}
        _, x_e, heavy_scale = self._expected_conversion(self.N_E, self.T, self.TE, self.P_ATM_SI)
        assert imf['e-'] == pytest.approx(x_e, rel=1e-12)
        assert imf['Ar'] == pytest.approx(0.7 * heavy_scale, rel=1e-12)
        assert imf['He'] == pytest.approx(0.3 * heavy_scale, rel=1e-12)
        assert sum(reactor.initial_mole_fractions.values()) == pytest.approx(1.0, abs=1e-12)

    def test_electron_density_round_trip(self, tmp_path):
        """Assertion 3: the reactor's own compute_volume reproduces the requested n_e."""
        body = self._preamble() + self._plasma_block("    electronDensity=(1e23,'m^-3'),\n")
        rmg = self._read(tmp_path, body)
        reactor = rmg.reaction_systems[0]
        core = list(reactor.initial_mole_fractions.keys())
        reactor.initialize_model(core, [], [], [])
        volume = reactor.compute_volume(reactor.y0)
        n_e_recovered = reactor.y0[reactor.electron_index] * constants.Na / volume
        # The conversion uses the reactor's own kB = R/Na (not constants.kB), so the
        # two sides of the equation of state agree by construction and the round-trip
        # is exact to floating point.
        assert n_e_recovered == pytest.approx(self.N_E, rel=1e-12)

    def test_directive_optional_electron_named_directly(self, tmp_path):
        """Assertion 4: naming the electron directly, with no density directive, still works."""
        body = self._preamble() + self._plasma_block(
            "", mole_fractions="{'Ar': 0.7, 'He': 0.2999, 'e-': 1e-4}")
        rmg = self._read(tmp_path, body)
        reactor = rmg.reaction_systems[0]
        imf = {k.label: v for k, v in reactor.initial_mole_fractions.items()}
        assert 'e-' in imf and imf['e-'] > 0.0
        assert sum(reactor.initial_mole_fractions.values()) == pytest.approx(1.0, abs=1e-12)
        core = list(reactor.initial_mole_fractions.keys())
        reactor.initialize_model(core, [], [], [])
        assert reactor.y0[reactor.electron_index] > 0.0

    # ---- guard refusals (assertion 5: one per guard) -----------------------

    def test_guard1_ranged_conditions_rejected(self, tmp_path):
        body = self._preamble() + self._plasma_block(
            "    electronDensity=(1e23,'m^-3'),\n", mole_fractions="{'Ar': 1.0}",
            temperature="[(1000,'K'),(2000,'K')]")
        with pytest.raises(InputError, match='electronDensity'):
            self._read(tmp_path, body)

    def test_guard2_te_units_rejected(self, tmp_path):
        body = self._preamble() + self._plasma_block(
            "    electronDensity=(1e23,'m^-3'),\n", mole_fractions="{'Ar': 1.0}",
            electron_temperature="(1,'eV')")
        with pytest.raises(InputError, match='electronTemperature'):
            self._read(tmp_path, body)

    def test_guard3_q_at_or_above_one_rejected(self, tmp_path):
        # n_e large enough that P_e >= P (q >= 1): physically impossible.
        body = self._preamble() + self._plasma_block(
            "    electronDensity=(1e25,'m^-3'),\n", mole_fractions="{'Ar': 1.0}")
        with pytest.raises(InputError, match='electronDensity'):
            self._read(tmp_path, body)

    def test_guard4_near_degenerate_q_rejected(self, tmp_path):
        from rmgpy.rmg.input import PLASMA_ELECTRON_PRESSURE_MARGIN as margin
        q_target = 1.0 - 0.5 * margin  # inside the degeneracy band, still q < 1
        # Use the builder's kB = R/Na so the recomputed q lands exactly at q_target.
        n_e = q_target * self.P_ATM_SI / ((constants.R / constants.Na) * self.TE)
        body = self._preamble() + self._plasma_block(
            f"    electronDensity=({n_e!r},'m^-3'),\n", mole_fractions="{'Ar': 1.0}")
        with pytest.raises(InputError, match='electronDensity'):
            self._read(tmp_path, body)

    def test_guard5_duplicate_specification_rejected(self, tmp_path):
        body = self._preamble() + self._plasma_block(
            "    electronDensity=(1e23,'m^-3'),\n",
            mole_fractions="{'Ar': 1.0, 'e-': 1e-4}")
        with pytest.raises(InputError, match='electronDensity'):
            self._read(tmp_path, body)

    def test_guard6_zero_density_rejected(self, tmp_path):
        body = self._preamble() + self._plasma_block(
            "    electronDensity=(0,'m^-3'),\n", mole_fractions="{'Ar': 1.0}")
        with pytest.raises(InputError, match='electronDensity'):
            self._read(tmp_path, body)

    def test_guard6_negative_density_rejected(self, tmp_path):
        body = self._preamble() + self._plasma_block(
            "    electronDensity=(-1.0,'m^-3'),\n", mole_fractions="{'Ar': 1.0}")
        with pytest.raises(InputError, match='electronDensity'):
            self._read(tmp_path, body)

    def test_guard6_infinite_density_rejected(self, tmp_path):
        # 1e400 evaluates to inf in the input file (no builtins needed).
        body = self._preamble() + self._plasma_block(
            "    electronDensity=(1e400,'m^-3'),\n", mole_fractions="{'Ar': 1.0}")
        with pytest.raises(InputError, match='electronDensity'):
            self._read(tmp_path, body)

    def test_guard6_nan_density_rejected(self, tmp_path):
        # 1e400 - 1e400 evaluates to nan in the input file (no builtins needed).
        body = self._preamble() + self._plasma_block(
            "    electronDensity=(1e400 - 1e400,'m^-3'),\n", mole_fractions="{'Ar': 1.0}")
        with pytest.raises(InputError, match='electronDensity'):
            self._read(tmp_path, body)

    def test_guard7_missing_electron_species_rejected(self, tmp_path):
        # Electron species is NOT declared; the density directive must refuse.
        body = self._heavy_species() + self._plasma_block(
            "    electronDensity=(1e23,'m^-3'),\n", mole_fractions="{'Ar': 1.0}")
        with pytest.raises(InputError, match='electron'):
            self._read(tmp_path, body)

    # ---- review-fix regression guards (I-067 round 2) ----------------------

    def test_guard_density_bare_number_rejected(self, tmp_path):
        # B1: a bare number carries no dimension; it must not be read as m^-3.
        body = self._preamble() + self._plasma_block(
            "    electronDensity=1e23,\n", mole_fractions="{'Ar': 1.0}")
        with pytest.raises(InputError, match='electronDensity'):
            self._read(tmp_path, body)

    def test_density_molecule_count_units_accepted_with_avogadro(self, tmp_path):
        # 'molecule/cm^3' is a count density used interchangeably with cm^-3 in the
        # plasma literature. RMG defines 'molecule' = mol/Na, so it must be accepted and
        # converted with the Avogadro factor: 1e11 molecule/cm^3 == 1e17 m^-3. Proving
        # x_e equality (not merely that it parsed) is what would catch a missing or
        # doubled Na.
        body_molecule = self._preamble() + self._plasma_block(
            "    electronDensity=(1e11,'molecule/cm^3'),\n", mole_fractions="{'Ar': 1.0}")
        body_si = self._preamble() + self._plasma_block(
            "    electronDensity=(1e17,'m^-3'),\n", mole_fractions="{'Ar': 1.0}")
        rmg_molecule = self._read(tmp_path, body_molecule)
        rmg_si = self._read(tmp_path, body_si)
        xe_molecule = {k.label: v for k, v in
                       rmg_molecule.reaction_systems[0].initial_mole_fractions.items()}['e-']
        xe_si = {k.label: v for k, v in
                 rmg_si.reaction_systems[0].initial_mole_fractions.items()}['e-']
        assert xe_molecule == pytest.approx(xe_si, rel=1e-12)

    def test_density_plural_molecules_units_accepted(self, tmp_path):
        # The plural 'molecules/cm^3' must be treated identically to 'molecule/cm^3'.
        body_plural = self._preamble() + self._plasma_block(
            "    electronDensity=(1e11,'molecules/cm^3'),\n", mole_fractions="{'Ar': 1.0}")
        body_si = self._preamble() + self._plasma_block(
            "    electronDensity=(1e17,'m^-3'),\n", mole_fractions="{'Ar': 1.0}")
        rmg_plural = self._read(tmp_path, body_plural)
        rmg_si = self._read(tmp_path, body_si)
        xe_plural = {k.label: v for k, v in
                     rmg_plural.reaction_systems[0].initial_mole_fractions.items()}['e-']
        xe_si = {k.label: v for k, v in
                 rmg_si.reaction_systems[0].initial_mole_fractions.items()}['e-']
        assert xe_plural == pytest.approx(xe_si, rel=1e-12)

    def test_density_reciprocal_volume_units_accepted(self, tmp_path):
        # '1/cm^3' is dimensionally inverse-volume (dim 1/m^3) and must be accepted as a
        # raw count density, equal to the cm^-3 spelling.
        body_recip = self._preamble() + self._plasma_block(
            "    electronDensity=(1e17,'1/cm^3'),\n", mole_fractions="{'Ar': 1.0}")
        body_cm = self._preamble() + self._plasma_block(
            "    electronDensity=(1e17,'cm^-3'),\n", mole_fractions="{'Ar': 1.0}")
        rmg_recip = self._read(tmp_path, body_recip)
        rmg_cm = self._read(tmp_path, body_cm)
        xe_recip = {k.label: v for k, v in
                    rmg_recip.reaction_systems[0].initial_mole_fractions.items()}['e-']
        xe_cm = {k.label: v for k, v in
                 rmg_cm.reaction_systems[0].initial_mole_fractions.items()}['e-']
        assert xe_recip == pytest.approx(xe_cm, rel=1e-12)

    def test_guard_density_molar_concentration_still_rejected(self, tmp_path):
        # A molar concentration is NOT a count density; it must stay rejected even though
        # it shares the molecule dimensionality.
        body = self._preamble() + self._plasma_block(
            "    electronDensity=(1e17,'mol/m^3'),\n", mole_fractions="{'Ar': 1.0}")
        with pytest.raises(InputError, match='electronDensity'):
            self._read(tmp_path, body)

    def test_density_cm3_units_accepted(self, tmp_path):
        # B1: cm^-3 is a valid number density; the guard must not reject it, and the
        # conversion must round-trip (1e17 cm^-3 == 1e23 m^-3).
        body = self._preamble() + self._plasma_block(
            "    electronDensity=(1e17,'cm^-3'),\n", mole_fractions="{'Ar': 1.0}")
        rmg = self._read(tmp_path, body)
        reactor = rmg.reaction_systems[0]
        core = list(reactor.initial_mole_fractions.keys())
        reactor.initialize_model(core, [], [], [])
        volume = reactor.compute_volume(reactor.y0)
        n_e_recovered = reactor.y0[reactor.electron_index] * constants.Na / volume
        assert n_e_recovered == pytest.approx(1e23, rel=1e-9)

    def test_shared_mole_fraction_dict_not_mutated(self, tmp_path):
        # B2: one dict bound once and passed to two reactor blocks must not have the
        # first reactor rewritten by the second (which would trip the duplicate guard).
        body = (
            self._preamble() +
            "imf = {'Ar': 0.7, 'He': 0.3}\n"
            "plasmaReactor(temperature=(1000,'K'), pressure=(1,'atm'), "
            "electronTemperature=(11604.5,'K'), electronDensity=(1e23,'m^-3'), "
            "initialMoleFractions=imf, terminationTime=(1,'s'))\n"
            "plasmaReactor(temperature=(1500,'K'), pressure=(1,'atm'), "
            "electronTemperature=(11604.5,'K'), electronDensity=(1e23,'m^-3'), "
            "initialMoleFractions=imf, terminationTime=(1,'s'))\n"
        )
        rmg = self._read(tmp_path, body)
        assert len(rmg.reaction_systems) == 2
        for reactor in rmg.reaction_systems:
            labels = {k.label for k in reactor.initial_mole_fractions}
            assert 'e-' in labels
            assert sum(reactor.initial_mole_fractions.values()) == pytest.approx(1.0, abs=1e-12)

    def test_guard_multiple_electron_species_rejected(self, tmp_path):
        # B3: more than one electron species is ambiguous; refuse before mutating,
        # rather than let the reactor reject it later at the wrong layer.
        body = (
            "species(label='e-', structure=adjacencyList(\"1 e u1 p0 c-1\"))\n"
            "species(label='eX', structure=adjacencyList(\"1 e u0 p0 c-1\"))\n"
            "species(label='Ar', structure=adjacencyList(\"1 Ar u0 p4 c0\"))\n"
            + self._plasma_block("    electronDensity=(1e23,'m^-3'),\n",
                                 mole_fractions="{'Ar': 1.0}")
        )
        with pytest.raises(InputError, match='electron'):
            self._read(tmp_path, body)

    def test_zero_pressure_reported_as_pressure(self, tmp_path):
        # B5: P=0 must raise a clear 'pressure' InputError, not a raw ZeroDivisionError
        # and not be misattributed to electronDensity.
        body = self._preamble() + self._plasma_block(
            "    electronDensity=(1e23,'m^-3'),\n", mole_fractions="{'Ar': 1.0}",
            pressure="(0,'Pa')")
        with pytest.raises(InputError, match='pressure'):
            self._read(tmp_path, body)

    def test_negative_temperature_reported_as_temperature(self, tmp_path):
        # B5: T is validated before the conversion arithmetic, with its own message.
        body = self._preamble() + self._plasma_block(
            "    electronDensity=(1e23,'m^-3'),\n", mole_fractions="{'Ar': 1.0}",
            temperature="(-5,'K')")
        with pytest.raises(InputError, match='temperature'):
            self._read(tmp_path, body)

    def test_bare_number_electron_temperature_rejected(self, tmp_path):
        # Minor: Te must be given as (value, 'K'); a bare number is refused.
        body = self._preamble() + self._plasma_block(
            "    electronDensity=(1e23,'m^-3'),\n", mole_fractions="{'Ar': 1.0}",
            electron_temperature="11604.5")
        with pytest.raises(InputError, match='electronTemperature'):
            self._read(tmp_path, body)

    # ---- non-finite mole fractions (I-067 round 4) -------------------------

    def test_nan_mole_fraction_rejected_plain_path(self, tmp_path):
        # nan < 0 is False, so without an explicit finiteness check a nan slips past
        # every guard. Plain path (no electronDensity), electron named directly.
        body = self._preamble() + self._plasma_block(
            "", mole_fractions="{'Ar': (1e400 - 1e400), 'He': 0.3, 'e-': 1e-4}")
        with pytest.raises(InputError, match='finite'):
            self._read(tmp_path, body)

    def test_inf_mole_fraction_rejected_plain_path(self, tmp_path):
        body = self._preamble() + self._plasma_block(
            "", mole_fractions="{'Ar': 1e400, 'He': 0.3, 'e-': 1e-4}")
        with pytest.raises(InputError, match='finite'):
            self._read(tmp_path, body)

    def test_nan_mole_fraction_rejected_density_path(self, tmp_path):
        # electronDensity path normalizes differently (heavy_total then scale); a nan
        # heavy fraction must be refused, not divided through.
        body = self._preamble() + self._plasma_block(
            "    electronDensity=(1e23,'m^-3'),\n",
            mole_fractions="{'Ar': (1e400 - 1e400), 'He': 0.3}")
        with pytest.raises(InputError, match='finite'):
            self._read(tmp_path, body)

    def test_inf_mole_fraction_rejected_density_path(self, tmp_path):
        body = self._preamble() + self._plasma_block(
            "    electronDensity=(1e23,'m^-3'),\n",
            mole_fractions="{'Ar': 1e400, 'He': 0.3}")
        with pytest.raises(InputError, match='finite'):
            self._read(tmp_path, body)

    # ---- save_input_file round-trip (I-067 round 4) ------------------------

    def test_save_input_file_round_trips_plasma_reactor(self, tmp_path):
        # A plasma reactor written by save_input_file must come back as a PlasmaReactor,
        # not a simpleReactor with electronTemperature dropped. Real round-trip: read ->
        # save -> RE-READ the whole saved file -> assert PlasmaReactor with T, P, Te and
        # every mole fraction intact (to the writer's %g precision).
        from rmgpy.solver.plasma import PlasmaReactor
        body = (
            "database(thermoLibraries=['primaryThermoLibrary'], reactionLibraries=[], "
            "seedMechanisms=[], kineticsFamilies='default')\n"
            + self._preamble()
            + self._plasma_block("    electronDensity=(1e23,'m^-3'),\n",
                                 mole_fractions="{'Ar': 0.7, 'He': 0.3}")
            + "simulator(atol=1e-16, rtol=1e-8)\n"
            + "model(toleranceMoveToCore=0.1, toleranceInterruptSimulation=0.1)\n"
        )
        rmg1 = self._read(tmp_path, body)
        reactor1 = rmg1.reaction_systems[0]
        assert isinstance(reactor1, PlasmaReactor)

        saved = tmp_path / "saved.py"
        inp.save_input_file(str(saved), rmg1)

        rmg2 = RMG()
        inp.read_input_file(str(saved), rmg2)
        reactor2 = rmg2.reaction_systems[0]
        assert isinstance(reactor2, PlasmaReactor)
        assert reactor2.T.value_si == pytest.approx(reactor1.T.value_si, rel=1e-5)
        assert reactor2.P.value_si == pytest.approx(reactor1.P.value_si, rel=1e-5)
        assert reactor2.Te.value_si == pytest.approx(reactor1.Te.value_si, rel=1e-5)
        mf1 = {k.label: v for k, v in reactor1.initial_mole_fractions.items()}
        mf2 = {k.label: v for k, v in reactor2.initial_mole_fractions.items()}
        assert set(mf1) == set(mf2)
        for label in mf1:
            assert mf2[label] == pytest.approx(mf1[label], rel=1e-5)

    def test_save_input_file_round_trips_kinetics_depositories_all(self, tmp_path):
        # The reader maps kineticsDepositories='all' to rmg.kinetics_depositories = None;
        # the writer must serialize None back to 'all' or the saved file cannot be re-read.
        # The list form already round-trips, so 'all' is the case that proves the sentinel.
        from rmgpy.solver.plasma import PlasmaReactor
        body = (
            "database(thermoLibraries=['primaryThermoLibrary'], reactionLibraries=[], "
            "seedMechanisms=[], kineticsFamilies='default', kineticsDepositories='all')\n"
            + self._preamble()
            + self._plasma_block("    electronDensity=(1e23,'m^-3'),\n",
                                 mole_fractions="{'Ar': 0.7, 'He': 0.3}")
            + "simulator(atol=1e-16, rtol=1e-8)\n"
            + "model(toleranceMoveToCore=0.1, toleranceInterruptSimulation=0.1)\n"
        )
        rmg1 = self._read(tmp_path, body)
        assert rmg1.kinetics_depositories is None  # 'all' -> None in the reader
        saved = tmp_path / "saved.py"
        inp.save_input_file(str(saved), rmg1)

        rmg2 = RMG()
        inp.read_input_file(str(saved), rmg2)  # must not raise on the 'all' sentinel
        assert isinstance(rmg2.reaction_systems[0], PlasmaReactor)
        assert rmg2.kinetics_depositories is None  # the sentinel survived the round-trip


class TestInputPlasmaChargeBalance:
    """
    Tests for ``chargeBalanceSpecies`` on the ``plasmaReactor`` directive (I-169).

    Like :class:`TestInputPlasmaReactor` above, every test drives the public
    :func:`rmgpy.rmg.input.read_input_file`: the whole complaint this keyword answers
    is about what a modeller can express in a deck, so a test that called
    ``plasma_reactor`` directly would not be testing it.

    The oracle these tests hold the implementation to is stated independently of it:
    a composition is charge balanced when ``sum(x_i * z_i) == 0``, computed from the
    resolved reactor composition and each species' own ``get_net_charge()``.
    """

    # ---- helpers -----------------------------------------------------------

    def _read(self, tmp_path, body):
        path = tmp_path / "input.py"
        path.write_text(body)
        rmg = RMG()
        inp.read_input_file(str(path), rmg)
        return rmg

    @staticmethod
    def _preamble():
        return (
            "species(label='e-', structure=adjacencyList(\"1 e u1 p0 c-1\"))\n"
            "species(label='Ar', structure=adjacencyList(\"1 Ar u0 p4 c0\"))\n"
            "species(label='Arp', structure=adjacencyList(\"\"\"\n"
            "multiplicity 2\n"
            "1 Ar u1 p3 c+1\n"
            "\"\"\"))\n"
        )

    @staticmethod
    def _block(extra, mole_fractions="{'Ar': 1.0}"):
        return (
            "plasmaReactor(\n"
            "    temperature=(1000,'K'),\n"
            "    pressure=(1,'atm'),\n"
            "    electronTemperature=(11604.5,'K'),\n"
            f"{extra}"
            f"    initialMoleFractions={mole_fractions},\n"
            "    terminationTime=(1,'s'),\n"
            ")\n"
        )

    @staticmethod
    def _net_charge(reactor):
        """sum(x_i * z_i) over the resolved composition -- the oracle, stated here."""
        return sum(value * spc.get_net_charge()
                   for spc, value in reactor.initial_mole_fractions.items())

    # ---- the neutral-seeding path ------------------------------------------

    def test_charge_balance_with_electron_density_is_neutral(self, tmp_path):
        # The whole point: no hand-derived constant in the deck, and the resulting
        # composition is neutral to roundoff rather than to the modeller's arithmetic.
        body = self._preamble() + self._block(
            "    electronDensity=(1e23,'m^-3'),\n"
            "    chargeBalanceSpecies='Arp',\n")
        reactor = self._read(tmp_path, body).reaction_systems[0]
        imf = {k.label: v for k, v in reactor.initial_mole_fractions.items()}
        assert 'Arp' in imf                                   # RMG created the entry
        assert imf['Arp'] == pytest.approx(imf['e-'], rel=1e-12)  # singly charged: x+ == x_e
        assert sum(reactor.initial_mole_fractions.values()) == pytest.approx(1.0, abs=1e-12)
        assert abs(self._net_charge(reactor)) < 1e-15

    def test_charge_balance_preserves_electron_density_round_trip(self, tmp_path):
        # Balancing must not disturb the electron density the deck asked for: the
        # cation is inserted on the heavy side, and the EOS must still recover n_e.
        body = self._preamble() + self._block(
            "    electronDensity=(1e23,'m^-3'),\n"
            "    chargeBalanceSpecies='Arp',\n")
        reactor = self._read(tmp_path, body).reaction_systems[0]
        core = list(reactor.initial_mole_fractions.keys())
        reactor.initialize_model(core, [], [], [])
        volume = reactor.compute_volume(reactor.y0)
        n_e_recovered = reactor.y0[reactor.electron_index] * constants.Na / volume
        assert n_e_recovered == pytest.approx(1e23, rel=1e-9)

    def test_charge_balance_with_explicit_electron_fraction(self, tmp_path):
        # The keyword must work in the OTHER branch too -- an explicitly typed electron
        # fraction, no electronDensity. A keyword that no-ops in one branch is a defect.
        body = self._preamble() + self._block(
            "    chargeBalanceSpecies='Arp',\n",
            mole_fractions="{'Ar': 0.9, 'e-': 1e-6}")
        reactor = self._read(tmp_path, body).reaction_systems[0]
        imf = {k.label: v for k, v in reactor.initial_mole_fractions.items()}
        assert imf['Arp'] == pytest.approx(imf['e-'], rel=1e-12)
        assert sum(reactor.initial_mole_fractions.values()) == pytest.approx(1.0, abs=1e-12)
        assert abs(self._net_charge(reactor)) < 1e-15

    def test_charge_balance_accounts_for_charge_already_present(self, tmp_path):
        # A deck that already types SOME cation must have only the remainder supplied,
        # not x_e again -- otherwise the "balance" overshoots into net positive.
        body = self._preamble() + self._block(
            "    chargeBalanceSpecies='Arp',\n",
            mole_fractions="{'Ar': 0.9, 'e-': 1e-6, 'Nep2': 2e-7}")
        body = body.replace(
            "plasmaReactor(",
            "species(label='Nep2', structure=adjacencyList(\"\"\"\n"
            "1 Ne u0 p3 c+2\n"
            "\"\"\"))\n"
            "plasmaReactor(", 1)
        reactor = self._read(tmp_path, body).reaction_systems[0]
        imf = {k.label: v for k, v in reactor.initial_mole_fractions.items()}
        # Arp must carry x_e - 2*x_Nep2 worth of charge, not x_e.
        assert imf['Arp'] == pytest.approx(imf['e-'] - 2.0 * imf['Nep2'], rel=1e-9)
        assert abs(self._net_charge(reactor)) < 1e-15

    def test_multiply_charged_balance_species(self, tmp_path):
        # z_b = +2: the balancing fraction is x_e/2, not x_e. Guards against an
        # implementation that assumes singly-charged ions.
        body = self._preamble() + self._block(
            "    electronDensity=(1e23,'m^-3'),\n"
            "    chargeBalanceSpecies='Nep2',\n")
        body = body.replace(
            "plasmaReactor(",
            "species(label='Nep2', structure=adjacencyList(\"\"\"\n"
            "1 Ne u0 p3 c+2\n"
            "\"\"\"))\n"
            "plasmaReactor(", 1)
        reactor = self._read(tmp_path, body).reaction_systems[0]
        imf = {k.label: v for k, v in reactor.initial_mole_fractions.items()}
        assert imf['Nep2'] == pytest.approx(imf['e-'] / 2.0, rel=1e-12)
        assert abs(self._net_charge(reactor)) < 1e-15

    def test_omitting_the_keyword_leaves_the_composition_untouched(self, tmp_path):
        # Existing decks must keep working unchanged -- bit for bit, not approximately.
        (tmp_path / 'a').mkdir()
        (tmp_path / 'b').mkdir()
        body = self._preamble() + self._block("    electronDensity=(1e23,'m^-3'),\n")
        r1 = self._read(tmp_path / 'a', body).reaction_systems[0]
        r2 = self._read(tmp_path / 'b', body).reaction_systems[0]
        imf1 = {k.label: v for k, v in r1.initial_mole_fractions.items()}
        imf2 = {k.label: v for k, v in r2.initial_mole_fractions.items()}
        assert imf1 == imf2
        assert 'Arp' not in imf1          # no cation invented behind the modeller's back
        assert imf1['Ar'] + imf1['e-'] == pytest.approx(1.0, abs=1e-15)

    # ---- parse-time guards -------------------------------------------------

    def test_undeclared_balance_species_rejected_at_parse_time(self, tmp_path):
        body = self._preamble() + self._block(
            "    electronDensity=(1e23,'m^-3'),\n"
            "    chargeBalanceSpecies='Krp',\n")
        with pytest.raises(InputError, match="Krp"):
            self._read(tmp_path, body)

    def test_uncharged_balance_species_rejected(self, tmp_path):
        body = self._preamble() + self._block(
            "    electronDensity=(1e23,'m^-3'),\n"
            "    chargeBalanceSpecies='Ar',\n")
        with pytest.raises(InputError, match="no net charge"):
            self._read(tmp_path, body)

    def test_electron_as_balance_species_rejected(self, tmp_path):
        body = self._preamble() + self._block(
            "    electronDensity=(1e23,'m^-3'),\n"
            "    chargeBalanceSpecies='e-',\n")
        with pytest.raises(InputError, match="electron pseudo-species"):
            self._read(tmp_path, body)

    def test_balance_species_also_in_mole_fractions_rejected(self, tmp_path):
        # Two sources of truth for one number -- the same rule electronDensity applies
        # to the electron.
        body = self._preamble() + self._block(
            "    electronDensity=(1e23,'m^-3'),\n"
            "    chargeBalanceSpecies='Arp',\n",
            mole_fractions="{'Ar': 1.0, 'Arp': 1e-7}")
        with pytest.raises(InputError, match="two sources of truth"):
            self._read(tmp_path, body)

    def test_non_string_balance_species_rejected(self, tmp_path):
        body = self._preamble() + self._block(
            "    electronDensity=(1e23,'m^-3'),\n"
            "    chargeBalanceSpecies=['Arp'],\n")
        with pytest.raises(InputError, match="string label"):
            self._read(tmp_path, body)

    def test_wrong_sign_balance_species_rejected(self, tmp_path):
        # An anion cannot cancel an excess of negative charge; the required fraction is
        # negative, which must be a named error rather than a negative mole fraction.
        body = self._preamble() + self._block(
            "    electronDensity=(1e23,'m^-3'),\n"
            "    chargeBalanceSpecies='Om',\n")
        body = body.replace(
            "plasmaReactor(",
            "species(label='Om', structure=adjacencyList(\"\"\"\n"
            "multiplicity 2\n"
            "1 O u1 p3 c-1\n"
            "\"\"\"))\n"
            "plasmaReactor(", 1)
        with pytest.raises(InputError, match="not a usable mole fraction"):
            self._read(tmp_path, body)
