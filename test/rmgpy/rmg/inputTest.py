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

import rmgpy.polymer
import rmgpy.rmg.input
import rmgpy.rmg.input as inp
from rmgpy.exceptions import InputError
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import CoreEdgeReactionModel
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
        Test that we can import Reaction Libraries using the non-tuple form.
        """
        global rmg
        # add database properties to RMG
        inp.database(reactionLibraries=["test"])
        assert isinstance(rmg.reaction_libraries[0], tuple)
        assert not rmg.reaction_libraries[0][1]

    def test_importing_database_reaction_libraries_from_false_tuple(self):
        """
        Test that we can import Reaction Libraries using the Tuple False form.
        """
        global rmg
        # add database properties to RMG
        inp.database(reactionLibraries=[("test", False)])
        assert isinstance(rmg.reaction_libraries[0], tuple)
        assert not rmg.reaction_libraries[0][1]

    def test_importing_database_reaction_libraries_from_true_tuple(self):
        """
        Test that we can import Reaction Libraries using the Tuple True form.
        """
        global rmg
        # add database properties to RMG
        inp.database(reactionLibraries=[("test", True)])
        assert isinstance(rmg.reaction_libraries[0], tuple)
        assert rmg.reaction_libraries[0][1]

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


class TestGeneratePolymerConstraints:
    """Unit tests for the generatePolymerConstraints input block parser."""

    def setup_method(self):
        self.rmg = RMG()
        inp.set_global_rmg(self.rmg)

    def teardown_method(self):
        inp.set_global_rmg(None)

    def test_valid_block_stores(self):
        inp.generate_polymer_constraints(maximumCarbonAtoms=30)
        assert self.rmg.polymer_constraints == {"maximumCarbonAtoms": 30}

    def test_mixed_block_with_nonbounding_key_accepted(self):
        inp.generate_polymer_constraints(maximumCarbonAtoms=30, allowSingletO2=True)
        assert self.rmg.polymer_constraints == {"maximumCarbonAtoms": 30, "allowSingletO2": True}

    def test_unknown_key_raises(self):
        with pytest.raises(InputError):
            inp.generate_polymer_constraints(bogusKey=3)

    def test_zero_keys_raises(self):
        with pytest.raises(InputError, match="at least one finite"):
            inp.generate_polymer_constraints()

    def test_only_nonbounding_keys_raises(self):
        with pytest.raises(InputError, match="at least one finite"):
            inp.generate_polymer_constraints(allowSingletO2=True, speciesCuttingThreshold=20)

    def test_all_maxima_unlimited_raises(self):
        with pytest.raises(InputError, match="at least one finite"):
            inp.generate_polymer_constraints(maximumCarbonAtoms=-1, maximumHeavyAtoms=-1)


class TestPolymerSizeThresholdParsing:
    def setup_method(self):
        self.rmg = RMG()
        rmgpy.rmg.input.rmg = self.rmg

    def teardown_method(self):
        rmgpy.rmg.input.rmg = None

    def test_polymer_size_threshold_parsed(self):
        rmgpy.rmg.input.generate_polymer_constraints(
            maximumHeavyAtoms=30, polymerSizeThreshold=15,
        )
        assert self.rmg.polymer_constraints["polymerSizeThreshold"] == 15
        assert self.rmg.polymer_constraints["maximumHeavyAtoms"] == 30

    def test_polymer_size_threshold_alone_still_requires_a_bound(self):
        with pytest.raises(InputError):
            rmgpy.rmg.input.generate_polymer_constraints(polymerSizeThreshold=15)


class TestPolymerConduitAdmissionParsing:
    """M18.4 opt-in: options(polymerConduitAdmission=...) parsing + the
    fail-closed type guard. Input files are exec'd, so a truthy non-bool
    (bool("off") is True) must RAISE, never silently enable this
    prediction-changing feature."""

    def setup_method(self):
        self.rmg = RMG()
        rmgpy.rmg.input.rmg = self.rmg

    def teardown_method(self):
        rmgpy.rmg.input.rmg = None

    def test_default_is_none_inherits_fallback(self):
        rmgpy.rmg.input.options()
        assert self.rmg.polymer_conduit_admission is None

    def test_true_and_false_parsed(self):
        rmgpy.rmg.input.options(polymerConduitAdmission=True)
        assert self.rmg.polymer_conduit_admission is True
        rmgpy.rmg.input.options(polymerConduitAdmission=False)
        assert self.rmg.polymer_conduit_admission is False

    def test_string_value_rejected_fails_closed(self):
        # bool("off") is True -> would silently ENABLE admission; must raise.
        with pytest.raises(InputError, match="polymerConduitAdmission"):
            rmgpy.rmg.input.options(polymerConduitAdmission="off")

    def test_int_value_rejected(self):
        with pytest.raises(InputError, match="polymerConduitAdmission"):
            rmgpy.rmg.input.options(polymerConduitAdmission=1)


class TestPolymerCopolymerDeck:
    """
    Deck-level tests for the copolymer composition (``polymer(monomers=[...])``).

    The unit tests in polymerTest.py pin the composition maths and the dyad
    proxy GRAPHS; these pin the thing only the deck path can establish -- that
    the dyad proxies actually reach the core as reactive species. That is the
    load-bearing step: reaction families only ever see species in the core, so
    a dyad proxy that is built but never registered would leave the pool
    carrying a copolymer's mass on a homopolymer's chemistry, silently.
    """

    ETHYLENE = '[CH2][CH2]'
    PROPYLENE = '[CH2][CH](C)'
    DIENE = "CC=C1CC2[CH][CH]C1C2"

    def setup_method(self):
        self.rmg = RMG()
        self.rmg.reaction_model = CoreEdgeReactionModel()
        self.rmg.initial_species = []
        inp.set_global_rmg(self.rmg)
        inp.species_dict = {}

    def teardown_method(self):
        inp.set_global_rmg(None)
        inp.species_dict = {}

    def _declare(self, units, **kwargs):
        params = dict(end_groups=['[CH3]', '[H]'], cutoff=3,
                      Mn=5000.0, Mw=10000.0, initial_mass=1.0)
        params.update(kwargs)
        return inp.polymer(
            label='EPDM',
            monomers=[dict(monomer=m, fraction=f) for m, f in units],
            **params)

    def test_dyad_proxies_are_registered_as_core_species(self):
        pool = self._declare([(self.ETHYLENE, 0.60), (self.PROPYLENE, 0.35),
                              (self.DIENE, 0.05)])
        # 3 units -> 6 unordered pairs, minus the dominant homo-dyad (which is
        # the baseline proxy, already registered as the pool itself).
        assert len(pool.dyad_proxy_species) == 5
        for spc in pool.dyad_proxy_species:
            assert spc in self.rmg.initial_species, \
                'dyad proxy never reached initial_species; families cannot react it'
            assert spc.reactive
            assert spc.props['copolymer_dyad_origin'][0] == pool.label

    def test_homopolymer_deck_registers_no_dyads(self):
        """The legacy declaration must be untouched: no extra core species."""
        pool = inp.polymer(label='PE', monomer=self.ETHYLENE,
                           end_groups=['[CH3]', '[H]'], cutoff=3,
                           Mn=5000.0, Mw=10000.0, initial_mass=1.0)
        assert pool.dyad_proxy_species == []
        assert pool.comonomers is None

    def test_deck_rejects_monomer_and_monomers_together(self):
        with pytest.raises(InputError, match='mutually exclusive'):
            inp.polymer(label='both', monomer=self.ETHYLENE,
                        monomers=[dict(monomer=self.PROPYLENE, fraction=1.0)],
                        end_groups=['[CH3]', '[H]'], Mn=5000.0, Mw=10000.0)

    def test_per_comonomer_monomer_product_is_registered(self):
        """A unit that unzips needs a live gas destination for its volatile."""
        pool = inp.polymer(
            label='EPDM_products',
            monomers=[dict(monomer=self.ETHYLENE, fraction=0.7,
                           monomer_product='C=C'),
                      dict(monomer=self.PROPYLENE, fraction=0.3,
                           monomer_product='C=CC')],
            end_groups=['[CH3]', '[H]'], cutoff=3, Mn=5000.0, Mw=10000.0)
        assert len(pool.comonomer_product_species) == 2
        for spc in pool.comonomer_product_species:
            assert spc in self.rmg.initial_species
            assert spc.reactive

    def test_composition_reaches_the_sidecar_block(self):
        """The artifact must carry the whole composition, not just the
        dominant unit -- a consumer reconstructing condensed mass from
        `monomer_smiles` alone would use the wrong repeat mass."""
        pool = self._declare([(self.ETHYLENE, 0.60), (self.PROPYLENE, 0.35),
                              (self.DIENE, 0.05)])
        block = rmgpy.polymer._serialize_pool_for_sidecar(pool)
        comp = block['composition']
        assert comp['statistics'] == 'bernoullian_random'
        assert len(comp['units']) == 3
        assert sum(u['fraction'] for u in comp['units']) == pytest.approx(1.0)
        assert comp['monomer_mw_g_mol'] == pytest.approx(pool.monomer_mw_g_mol)

    def test_dyad_proxies_route_to_the_polymer_constraint_tier(self):
        """
        Dyad proxies are much larger than any volatile the deck should generate
        (an ENB--ENB dyad is C28), so they must be bounded by the POLYMER tier,
        not the gas tier. Routing is by heavy-atom count, so the smallest dyad
        proxy sets the usable polymerSizeThreshold -- a threshold above it would
        drop that proxy onto the gas tier, where it would be refused (or force a
        gas bound loose enough for RMG to generate C28 volatiles).
        """
        from rmgpy.constraints import is_polymer_constraint_member
        pool = self._declare([(self.ETHYLENE, 0.7101), (self.PROPYLENE, 0.2761),
                              (self.DIENE, 0.0138)])
        heavy = [sum(1 for a in spc.molecule[0].atoms if a.symbol != 'H')
                 for spc in pool.dyad_proxy_species]
        threshold = min(heavy)
        constraints = {'polymerSizeThreshold': threshold}
        for spc in pool.dyad_proxy_species:
            assert is_polymer_constraint_member(spc, polymer_constraints=constraints), \
                'a dyad proxy fell through to the gas constraint tier'
        # The example deck pins this threshold; if the proxy construction ever
        # shrinks the smallest dyad, examples/rmg/epdm/input.py must move too.
        assert threshold == 8

    def test_composition_is_declarable_in_a_real_deck_context(self):
        """
        Input files are exec'd with ``__builtins__`` set to None, so a deck
        cannot call ``dict(...)`` -- only dict LITERALS survive. Calling
        inp.polymer() from Python (as the tests above do) hides that, which is
        exactly how a deck that cannot be read got written; this test executes
        the declaration the way read_input_file does.
        """
        deck = (
            "polymer(label='EPDM',\n"
            "        monomers=[{'monomer': '[CH2][CH2]', 'fraction': 0.7101},\n"
            "                  {'monomer': '[CH2][CH](C)', 'fraction': 0.2761},\n"
            "                  {'monomer': 'CC=C1CC2[CH][CH]C1C2', 'fraction': 0.0138}],\n"
            "        end_groups=['[CH3]', '[H]'], cutoff=3,\n"
            "        Mn=5000.0, Mw=10000.0, initial_mass=1.0)\n")
        global_context = {'__builtins__': None}
        local_context = {'__builtins__': None, 'True': True, 'False': False,
                         'polymer': inp.polymer}
        exec(deck, global_context, local_context)   # must not raise
        pool = self.rmg.reaction_model.core.species  # deck path registered species
        assert pool is not None
        assert len(self.rmg.initial_species) >= 5   # 5 dyad proxies at minimum

    def test_dict_builtin_is_not_available_to_decks(self):
        """Negative control: pins WHY the literals above are required."""
        with pytest.raises(TypeError):
            exec("x = dict(a=1)", {'__builtins__': None}, {'__builtins__': None})
