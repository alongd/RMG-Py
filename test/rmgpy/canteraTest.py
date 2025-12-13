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


import cantera as ct
import os
import shutil
import numpy as np
import pytest
import yaml

import rmgpy.constants as constants
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.kinetics import (
    Arrhenius,
    PDepArrhenius,
    MultiArrhenius,
    Chebyshev,
    Troe,
    Lindemann,
    ThirdBody,
    TwoTemperaturePlasma,
    ElectronCollisionPlasma
)
from rmgpy.thermo import NASA, NASAPolynomial
from rmgpy.transport import TransportData
from rmgpy.molecule import Molecule
from rmgpy.cantera import (
    CanteraWriter,
    save_cantera_files,
    species_to_dict,
    reaction_to_dict_list,
    generate_cantera_data
)


class TestCanteraWriter:

    def setup_method(self):
        """
        Create a temporary directory for file I/O tests.
        """
        base_dir = os.path.dirname(os.path.abspath(__file__))
        self.tmp_dir = os.path.join(base_dir, 'tmp')

        # Ensure a clean start: delete if exists, then create
        if os.path.exists(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.makedirs(self.tmp_dir)

    def teardown_method(self):
        """
        Clean up the temporary directory after tests.
        """
        shutil.rmtree(self.tmp_dir)

    # --------------------------------------------------------------------------
    #    Helper: Create Dummy RMG Objects
    # --------------------------------------------------------------------------

    def _create_dummy_species(self, label, formula, is_electron=False):
        """Helper to create a functional RMG Species object with thermo/transport"""
        if is_electron:
            # RMG represents electron usually as 'e' with Element 'E' or just charge
            # For the writer test, the label is the critical part
            mol = Molecule().from_adjacency_list("1 e u1 p0 c-1")
            sp = Species(label=label, molecule=[mol])
        else:
            sp = Species(label=label).from_smiles(formula)# Create dummy NASA7 Thermo
        coeffs = [1.0, 0.0, 0.0, 0.0, 0.0, -100.0, 1.0]
        poly_low = NASAPolynomial(coeffs=coeffs, Tmin=(200, 'K'), Tmax=(1000, 'K'))
        poly_high = NASAPolynomial(coeffs=coeffs, Tmin=(1000, 'K'), Tmax=(6000, 'K'))
        sp.thermo = NASA(polynomials=[poly_low, poly_high], Tmin=(200, 'K'), Tmax=(6000, 'K'))

        # Create dummy Transport
        if not is_electron:
            # Determine correct geometry index to satisfy Cantera validation
            # 0: atom, 1: linear, 2: nonlinear
            num_atoms = len(sp.molecule[0].atoms)
            if num_atoms == 1:
                shape_idx = 0
            elif num_atoms == 2:
                shape_idx = 1
            else:
                shape_idx = 2

            sp.transport_data = TransportData(
                shapeIndex=shape_idx,
                sigma=(3.0, 'angstrom'),
                epsilon=(100.0, 'K'),
                dipoleMoment=(0.0, 'De'),
                polarizability=(0.0, 'angstrom^3'),
                rotrelaxcollnum=1.0
            )
        return sp

    # --------------------------------------------------------------------------
    #    Unit Tests: Species Translation
    # --------------------------------------------------------------------------

    def test_species_to_dict_standard(self):
        """Test conversion of a standard gas species."""
        sp = self._create_dummy_species("H2", "[H][H]")
        d = species_to_dict(sp)

        assert d['name'] == "H2"
        assert 'composition' in d
        assert d['thermo']['model'] == 'NASA7'
        assert len(d['thermo']['data']) == 2

        # Verify Transport
        assert 'transport' in d
        assert d['transport']['model'] == 'gas'
        assert d['transport']['geometry'] == 'linear'
        # Diameter should be in meters (SI)
        assert np.isclose(d['transport']['diameter'], 3.0e-10)

    def test_species_to_dict_electron(self):
        """Test conversion of electron species (special thermo/transport)."""
        # Test variations of labels
        for label in ['e', 'E', 'e-']:
            sp = self._create_dummy_species(label, "", is_electron=True)
            d = species_to_dict(sp)

            assert d['name'] == label
            assert d['composition'] == {'E': 1}
            assert d['thermo']['model'] == 'constant-cp'
            assert d['transport']['model'] == 'gas'

    # --------------------------------------------------------------------------
    #    Unit Tests: Reaction Translation
    # --------------------------------------------------------------------------

    def test_reaction_to_dict_arrhenius(self):
        """Test standard Arrhenius kinetics."""
        r = self._create_dummy_species("R", "[H]")
        p = self._create_dummy_species("P", "[H]")
        rxn = Reaction(
            reactants=[r], products=[p],
            kinetics=Arrhenius(A=(1e10, "s^-1"), n=0.5, Ea=(10, "kJ/mol"), T0=(1, "K"))
        )

        entries = reaction_to_dict_list(rxn)
        assert len(entries) == 1
        data = entries[0]

        assert data['equation'] == "R <=> P"  # Reversible default
        assert 'rate-constant' in data
        assert np.isclose(data['rate-constant']['A'], 1e10)
        assert np.isclose(data['rate-constant']['b'], 0.5)
        assert np.isclose(data['rate-constant']['Ea'], 10000.0)

    def test_reaction_to_dict_two_temp_plasma(self):
        """Test TwoTemperaturePlasma kinetics structure."""
        e = self._create_dummy_species("e", "", is_electron=True)
        r = self._create_dummy_species("R", "[H]")

        rxn = Reaction(
            reactants=[r, e], products=[r, e],
            kinetics=TwoTemperaturePlasma(
                A=(1e-10, "cm^3/(molecule*s)"), n=1.0,
                Ea_g=(1000, "J/mol"), Ea_e=(2.0, "eV/molecule")
            )
        )

        entries = reaction_to_dict_list(rxn)
        data = entries[0]

        assert data['type'] == 'two-temperature-plasma'
        rc = data['rate-constant']

        # Check units conversion
        # A: 1e-10 cm3/s * 1e-6 (m3/cm3) * Na ~= 6.022e7
        assert np.isclose(rc['A'], 1e-10 * 1e-6 * 6.02214e23, rtol=1e-3)
        assert rc['Ea-gas'] == 1000.0
        # Ea-e: stored as J/mol in RMG
        assert np.isclose(rc['Ea-electron'], 2.0 * 96485.33, rtol=1e-3)

    def test_reaction_to_dict_electron_collision(self):
        """Test ElectronCollisionPlasma kinetics and unit conversion."""
        e = self._create_dummy_species("e", "", is_electron=True)
        r = self._create_dummy_species("R", "[H]")

        energies_eV = np.array([0.0, 10.0, 50.0])
        sigma_m2 = np.array([0.0, 1e-20, 0.5e-20])

        rxn = Reaction(
            reactants=[r, e], products=[r, e],
            kinetics=ElectronCollisionPlasma(
                energies=(energies_eV, "eV/molecule"),
                sigma=(sigma_m2, "m^2")
            )
        )

        entries = reaction_to_dict_list(rxn)
        data = entries[0]

        assert data['type'] == 'electron-collision-plasma'
        assert 'energy-levels' in data
        assert 'cross-sections' in data

        # Verify J/mol -> eV conversion happened correctly
        assert np.allclose(data['energy-levels'], energies_eV, atol=1e-6)
        assert np.allclose(data['cross-sections'], sigma_m2)

    def test_reaction_to_dict_duplicates(self):
        """Test that MultiKinetics objects result in multiple YAML entries."""
        r = self._create_dummy_species("R", "[H]")
        k1 = Arrhenius(A=(1e10, "s^-1"), n=0, Ea=(0, "J/mol"), T0=(1, "K"))
        k2 = Arrhenius(A=(2e10, "s^-1"), n=0, Ea=(0, "J/mol"), T0=(1, "K"))

        rxn = Reaction(
            reactants=[r], products=[r],
            kinetics=MultiArrhenius(arrhenius=[k1, k2]),
            duplicate=True
        )

        entries = reaction_to_dict_list(rxn)
        assert len(entries) == 2
        assert entries[0]['rate-constant']['A'] == 1e10
        assert entries[1]['rate-constant']['A'] == 2e10
        assert entries[0].get('duplicate') is True

    def test_reaction_to_dict_troe(self):
        """Test Falloff/Troe serialization."""
        r = self._create_dummy_species("R", "[H]")
        M = self._create_dummy_species("M", "[Ar]")

        # Troe
        k_high = Arrhenius(A=(1e14, "s^-1"), n=0, Ea=(10, "kJ/mol"), T0=(1, "K"))
        k_low = Arrhenius(A=(1e20, "cm^3/(mol*s)"), n=0, Ea=(10, "kJ/mol"), T0=(1, "K"))

        troe = Troe(
            arrheniusHigh=k_high, arrheniusLow=k_low,
            alpha=0.5, T3=(100, "K"), T1=(200, "K"), T2=(300, "K"),
            efficiencies={M.molecule[0]: 2.0}
        )

        rxn = Reaction(reactants=[r], products=[r], kinetics=troe)
        entries = reaction_to_dict_list(rxn, species_list=[r, M])
        data = entries[0]

        assert data['type'] == 'falloff'
        assert 'Troe' in data
        assert data['Troe']['A'] == 0.5
        assert data['Troe']['T2'] == 300.0
        # Efficiencies should map label -> val
        assert data['efficiencies'] == {"M": 2.0}

    # --------------------------------------------------------------------------
    #    Unit Tests: Phase Generation
    # --------------------------------------------------------------------------

    def test_generate_cantera_data_detects_plasma(self):
        """Test that the writer detects 'e' and sets thermo: plasma."""
        h2 = self._create_dummy_species("H2", "[H][H]")
        e = self._create_dummy_species("e", "", is_electron=True)

        # Case 1: No Electron
        data = generate_cantera_data([h2], [], is_plasma=False)
        phase = data['phases'][0]
        assert phase['thermo'] == 'ideal-gas'
        assert phase['transport'] == 'mixture-averaged'

        # Case 2: With Electron (passed explicitly via is_plasma flag logic in save_cantera_files)
        # Note: generate_cantera_data takes the boolean, the logic to determine it is in save_cantera_model
        data = generate_cantera_data([h2, e], [], is_plasma=True)
        phase = data['phases'][0]
        assert phase['thermo'] == 'plasma'
        assert phase['transport'] == 'ionized-gas'
        assert 'electron-energy-distribution' in phase
        assert 'E' in phase['elements']

    # --------------------------------------------------------------------------
    #    Integration Test: Full Cycle (Write -> Read)
    # --------------------------------------------------------------------------

    def test_full_integration_plasma_model(self):
        """
        Create a comprehensive RMG model, write it to disk, and load it in Cantera
        to ensure all fields are valid and parsed correctly.
        """

        # 1. Create Model Components
        H2 = self._create_dummy_species("H2", "[H][H]")
        H = self._create_dummy_species("H", "[H]")
        e = self._create_dummy_species("e", "", is_electron=True)
        Ar = self._create_dummy_species("Ar", "[Ar]")

        species = [H2, H, e, Ar]

        # Reaction 1: Arrhenius
        r1 = Reaction(
            reactants=[H2], products=[H, H],
            kinetics=Arrhenius(A=(1e13, "s^-1"), n=0, Ea=(400, "kJ/mol"), T0=(1, "K"))
        )

        # Reaction 2: TwoTemp Plasma
        r2 = Reaction(
            reactants=[H, e], products=[H, e],
            kinetics=TwoTemperaturePlasma(
                A=(2e-10, "cm^3/(molecule*s)"), n=0.5, Ea_g=(0, "J/mol"), Ea_e=(1.0, "eV/molecule")
            )
        )

        # Reaction 3: Electron Collision
        energies = np.array([0.0, 15.0, 100.0])
        sigma = np.array([0.0, 2e-20, 0.5e-20])
        r3 = Reaction(
            reactants=[H2, e], products=[H, H, e],
            kinetics=ElectronCollisionPlasma(energies=(energies, "eV/molecule"), sigma=(sigma, "m^2"))
        )

        # Reaction 4: Third Body
        r4 = Reaction(
            reactants=[H, H], products=[H2],
            kinetics=ThirdBody(
                arrheniusLow=Arrhenius(A=(1e18, "cm^6/(mol^2*s)"), n=-1, Ea=(0, "J/mol"), T0=(1, "K")),
                efficiencies={Ar.molecule[0]: 0.7}
            )
        )

        reactions = [r1, r2, r3, r4]

        # 2. Mock RMG Object Structure
        # The writer expects: rmg.output_directory and rmg.reaction_model.core
        class MockCore:
            def __init__(self):
                self.species = species
                self.reactions = reactions

        class MockModel:
            def __init__(self):
                self.core = MockCore()
                self.edge = MockCore()  # Empty for now

        class MockRMG:
            def __init__(self, out_dir):
                self.output_directory = out_dir
                self.reaction_model = MockModel()
                self.save_edge_species = False

        mock_rmg = MockRMG(self.tmp_dir)

        # 3. Run Writer
        save_cantera_files(mock_rmg)

        # Verify files exist
        yaml_file = os.path.join(self.tmp_dir, "cantera", "chem.yaml")
        versioned_file = os.path.join(self.tmp_dir, "cantera", "chem0004.yaml")
        assert os.path.exists(yaml_file)
        assert os.path.exists(versioned_file)

        # 4. Validate with Cantera
        # Note: Cantera 3.0+ is required for the plasma features.
        # If running on 2.6, this might fail on the 'thermo: plasma' field.
        try:
            # Suppress warnings about thermodynamic inconsistencies in dummy data
            sol = ct.Solution(yaml_file)
        except Exception as e:
            pytest.fail(f"Cantera failed to load the generated YAML: {e}")

        # 5. Assertions on Loaded Solution

        # Check Species
        assert sol.n_species == 4
        assert sol.species_index("e") >= 0

        # Check Phase
        # In Cantera 3.0+, we can verify it's a plasma phase
        # Note: 'name' usually defaults to 'gas' unless specified
        # Check transport model if accessible
        assert sol.transport_model == "ionized-gas"

        # Check Reactions
        assert sol.n_reactions == 4

        # Check Reaction 2 (TwoTemp)
        ct_r2 = sol.reaction(1)
        # Verify it loaded as TwoTempPlasma type (name varies by Cantera version)
        # Cantera 3.0: reaction_type is string "TwoTempPlasma" or similar
        assert "TwoTemp" in ct_r2.reaction_type or isinstance(ct_r2.rate, ct.TwoTempPlasmaRate)

        # Check Reaction 3 (ElectronCollision)
        ct_r3 = sol.reaction(2)
        assert "ElectronCollision" in ct_r3.reaction_type or isinstance(ct_r3.rate, ct.ElectronCollisionPlasmaRate)

        # Check Reaction 4 (Efficiencies)
        ct_r4 = sol.reaction(3)
        assert "three-body" in ct_r4.reaction_type or "ThreeBody" in ct_r4.reaction_type
        # Verify Ar efficiency
        assert np.isclose(ct_r4.third_body.efficiencies["Ar"], 0.7)

    def test_reaction_to_dict_pdep_arrhenius(self):
        """Test Pressure-Dependent Arrhenius (PLOG) structure."""
        r = self._create_dummy_species("R", "[H]")
        p = self._create_dummy_species("P", "[H]")

        # Create two Arrhenius expressions for different pressures
        k_low = Arrhenius(A=(1e10, "s^-1"), n=0, Ea=(10, "kJ/mol"), T0=(1, "K"))
        k_high = Arrhenius(A=(1e12, "s^-1"), n=0, Ea=(15, "kJ/mol"), T0=(1, "K"))

        # Pressures: 0.1 atm and 1.0 atm
        # RMG stores pressures in SI (Pa) internally
        pdep = PDepArrhenius(
            pressures=([0.1, 1.0], "atm"),
            arrhenius=[k_low, k_high]
        )

        rxn = Reaction(reactants=[r], products=[p], kinetics=pdep)

        entries = reaction_to_dict_list(rxn)
        data = entries[0]

        assert data['type'] == 'pressure-dependent-Arrhenius'
        rates = data['rate-constants']
        assert len(rates) == 2

        # Verify first pressure point (0.1 atm -> 10132.5 Pa)
        assert np.isclose(rates[0]['P'], 0.1 * 101325.0)
        assert np.isclose(rates[0]['A'], 1e10)
        assert np.isclose(rates[0]['Ea'], 10000.0)

        # Verify second pressure point (1.0 atm -> 101325.0 Pa)
        assert np.isclose(rates[1]['P'], 1.0 * 101325.0)
        assert np.isclose(rates[1]['A'], 1e12)
        assert np.isclose(rates[1]['Ea'], 15000.0)

    def test_reaction_to_dict_chebyshev(self):
        """Test Chebyshev kinetics structure."""
        r = self._create_dummy_species("R", "[H]")

        # 2x2 Coefficients matrix
        coeffs = np.array([[1.0, 2.0], [3.0, 4.0]])
        cheb = Chebyshev(
            Tmin=(300, "K"), Tmax=(2000, "K"),
            Pmin=(0.01, "atm"), Pmax=(100, "atm"),
            coeffs=coeffs,
            kunits="s^-1"
        )

        rxn = Reaction(reactants=[r], products=[r], kinetics=cheb)

        entries = reaction_to_dict_list(rxn)
        data = entries[0]

        assert data['type'] == 'Chebyshev'

        # Check Ranges
        assert np.allclose(data['temperature-range'], [300.0, 2000.0])
        assert np.allclose(data['pressure-range'], [0.01 * 101325.0, 100 * 101325.0])

        # Check Coefficients
        # RMG stores coeffs in .value_si as numpy array, Writer converts to list
        assert np.allclose(data['data'], coeffs)

    def test_reaction_to_dict_lindemann(self):
        """Test Lindemann (Falloff without Troe parameters)."""
        r = self._create_dummy_species("R", "[H]")
        M = self._create_dummy_species("M", "[Ar]")

        # High Pressure Limit
        k_high = Arrhenius(A=(1e14, "s^-1"), n=0, Ea=(10, "kJ/mol"), T0=(1, "K"))
        # Low Pressure Limit (Third body type units)
        k_low = Arrhenius(A=(1e20, "cm^3/(mol*s)"), n=0, Ea=(10, "kJ/mol"), T0=(1, "K"))

        lind = Lindemann(
            arrheniusHigh=k_high,
            arrheniusLow=k_low,
            efficiencies={M.molecule[0]: 5.0}
        )

        rxn = Reaction(reactants=[r], products=[r], kinetics=lind)

        entries = reaction_to_dict_list(rxn, species_list=[r, M])
        data = entries[0]

        assert data['type'] == 'falloff'
        assert 'high-P-rate-constant' in data
        assert 'low-P-rate-constant' in data

        # Verify High P
        assert np.isclose(data['high-P-rate-constant']['A'], 1e14)

        # Verify Low P
        # 1e20 cm^3/(mol*s) = 1e20 * 1e-6 m^3/(mol*s) = 1e14
        assert np.isclose(data['low-P-rate-constant']['A'], 1e14)

        # Verify Efficiencies
        assert data['efficiencies'] == {"M": 5.0}

        # Ensure no Troe parameters leaked in
        assert 'Troe' not in data

    def test_cantera_writer_class_listener(self):
        """
        Test the CanteraWriter class directly to ensure it correctly initializes
        subdirectories and triggers the save on update().
        """
        # 1. Initialize Writer
        # This should create the 'cantera' subdirectory immediately
        writer = CanteraWriter(self.tmp_dir)
        cantera_dir = os.path.join(self.tmp_dir, 'cantera')
        assert os.path.exists(cantera_dir)
        assert os.path.isdir(cantera_dir)

        # 2. Setup Mock RMG Subject
        mock_rmg = self._create_dummy_model()

        # 3. Trigger Update
        # This simulates the RMG solver calling the listener after an iteration
        writer.update(mock_rmg)

        # 4. Verify Files
        # Based on _create_dummy_model, we have 3 species
        versioned_file = os.path.join(cantera_dir, 'chem0003.yaml')
        latest_file = os.path.join(cantera_dir, 'chem.yaml')

        assert os.path.exists(versioned_file)
        assert os.path.exists(latest_file)

        # Verify content briefly (ensure it's not empty)
        with open(latest_file, 'r') as f:
            content = f.read()
            assert "generator: RMG-Py CanteraWriter" in content
            assert "phases:" in content
            assert "species:" in content

    def _create_dummy_model(self):
        """Creates a mock object structure resembling RMG.reaction_model"""

        # 1. Species
        sp_H2 = self._create_dummy_species("H2", "[H][H]")
        sp_H = self._create_dummy_species("H", "[H]")
        sp_e = self._create_dummy_species("e", "", is_electron=True)

        species_list = [sp_H2, sp_H, sp_e]

        # 2. Reactions
        rxn_arr = Reaction(
            reactants=[sp_H2], products=[sp_H, sp_H],
            kinetics=Arrhenius(A=(1e13, "s^-1"), n=0.0, Ea=(200, "kJ/mol"), T0=(1, "K"))
        )

        reaction_list = [rxn_arr]

        # Mock Object Structure
        class MockCore:
            def __init__(self, s, r):
                self.species = s
                self.reactions = r

        class MockModel:
            def __init__(self, core):
                self.core = core
                self.edge = MockCore([], [])
                self.output_species_list = []
                self.output_reaction_list = []

        class MockRMG:
            def __init__(self, out_dir, model):
                self.output_directory = out_dir
                self.reaction_model = model
                self.save_edge_species = False

        return MockRMG(self.tmp_dir, MockModel(MockCore(species_list, reaction_list)))


class TestOxygenPlasmaReplication:

    def setup_method(self):
        """
        Setup paths.
        The reference file is expected to be in ../test_data/cantera/ relative to this file.
        """
        # Define paths
        self.base_dir = os.path.dirname(os.path.abspath(__file__))
        self.ref_yaml_path = os.path.join(self.base_dir, 'test_data', 'cantera', 'oxygen-plasma-itikawa.yaml')

        # Temp dir for generated file
        self.tmp_dir = os.path.join(self.base_dir, 'tmp_oxygen_test')
        if os.path.exists(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.makedirs(self.tmp_dir)

    def teardown_method(self):
        shutil.rmtree(self.tmp_dir)

    def _create_species_objects(self):
        """Create RMG Species objects matching the provided data."""

        def make_nasa(coeffs_low, coeffs_high, Tmin, Tmid=1000.0, Tmax=5000.0):
            return NASA(
                polynomials=[
                    NASAPolynomial(coeffs=coeffs_low, Tmin=(Tmin, 'K'), Tmax=(Tmid, 'K')),
                    NASAPolynomial(coeffs=coeffs_high, Tmin=(Tmid, 'K'), Tmax=(Tmax, 'K'))
                ],
                Tmin=(Tmin, 'K'), Tmax=(Tmax, 'K')
            )

        # 1. Electron (e)
        e = Species(label="e", molecule=[Molecule().from_adjacency_list("1 e u1 p0 c-1")])

        # 2. O
        coeffs_low_O = [2.542060E+00, -2.755062E-05, -3.102803E-09, 4.551067E-12, -4.368052E-16, 2.923080E+04,
                        4.920308E+00]
        coeffs_high_O = [2.946429E+00, -1.638166E-03, 2.421032E-06, -1.602843E-09, 3.890696E-13, 2.914764E+04,
                         2.963995E+00]
        O = Species(label="O", molecule=[Molecule().from_adjacency_list("multiplicity 3\n1 O u2 p2 c0")])
        O.thermo = make_nasa(coeffs_low_O, coeffs_high_O, Tmin=298.0)
        O.transport_data = TransportData(shapeIndex=0, sigma=(2.75, 'angstrom'), epsilon=(665.16 / 8.314, 'K'))

        # 3. O2
        coeffs_low_O2 = [3.697578E+00, 6.135197E-04, -1.258842E-07, 1.775281E-11, -1.136435E-15, -1.233930E+03,
                         3.189166E+00]
        coeffs_high_O2 = [3.212936E+00, 1.127486E-03, -5.756150E-07, 1.313877E-09, -8.768554E-13, -1.005249E+03,
                          6.034738E+00]
        O2 = Species(label="O2", molecule=[
            Molecule().from_adjacency_list("multiplicity 3\n1 O u1 p2 c0 {2,S}\n2 O u1 p2 c0 {1,S}")])
        O2.thermo = make_nasa(coeffs_low_O2, coeffs_high_O2, Tmin=298.0)
        O2.transport_data = TransportData(shapeIndex=1, sigma=(3.46, 'angstrom'), epsilon=(892.97 / 8.314, 'K'),
                                          polarizability=(1.60, 'angstrom^3'), rotrelaxcollnum=3.80)

        # 4. O2-
        coeffs_low_O2m = [3.95666294E+00, 5.98141823E-04, -2.12133905E-07, 3.63267581E-11, -2.24989228E-15,
                          -7.06287229E+03, 2.27871017E+00]
        coeffs_high_O2m = [3.66442522E+00, -9.28741138E-04, 6.45477082E-06, -7.74703380E-09, 2.93332662E-12,
                           -6.87076983E+03, 4.35140681E+00]
        O2m = Species(label="O2-", molecule=[
            Molecule().from_adjacency_list("multiplicity 2\n1 O u1 p2 c0 {2,S}\n2 O u0 p3 c-1 {1,S}")])
        O2m.thermo = make_nasa(coeffs_low_O2m, coeffs_high_O2m, Tmin=298.15, Tmax=6000.0)
        O2m.transport_data = TransportData(shapeIndex=1, sigma=(3.46, 'angstrom'), epsilon=(892.97 / 8.314, 'K'),
                                           polarizability=(1.60, 'angstrom^3'), rotrelaxcollnum=3.80)

        # 5. O-
        coeffs_low_Om = [2.54474869E+00, -4.66695513E-05, 1.84912357E-08, -3.18159223E-12, 1.98962956E-16,
                         1.15042089E+04, 4.52131015E+00]
        coeffs_high_Om = [2.90805921E+00, -1.69804907E-03, 2.98069955E-06, -2.43835127E-09, 7.61229311E-13,
                          1.14357717E+04, 2.80339097E+00]
        Om = Species(label="O-", molecule=[Molecule().from_adjacency_list("1 O u1 p3 c-1")])
        Om.thermo = make_nasa(coeffs_low_Om, coeffs_high_Om, Tmin=298.15, Tmax=6000.0)
        Om.transport_data = TransportData(shapeIndex=0, sigma=(2.75, 'angstrom'), epsilon=(665.16 / 8.314, 'K'))

        # 6. O2+
        coeffs_low_O2p = [4.61017, -0.00635952, 1.42426e-05, -1.20998e-08, 3.70957e-12, 139742, -0.201327]
        coeffs_high_O2p = [3.31676, 0.00111522, -3.83493e-07, 5.72785e-11, -2.77648e-15, 139877, 5.44726]
        O2p = Species(label="O2+", molecule=[Molecule().from_adjacency_list("1 O u0 p2 c+1 {2,S}\n2 O u1 p2 c0 {1,S}")])
        O2p.thermo = make_nasa(coeffs_low_O2p, coeffs_high_O2p, Tmin=298.15, Tmax=6000.0)
        O2p.transport_data = TransportData(shapeIndex=1, sigma=(3.46, 'angstrom'), epsilon=(892.97 / 8.314, 'K'),
                                           polarizability=(1.60, 'angstrom^3'), rotrelaxcollnum=3.80)

        return [e, O, O2, O2m, Om, O2p]

    def _create_reaction_objects(self, sp_map):
        e, O, O2, O2m, Om, O2p = sp_map["e"], sp_map["O"], sp_map["O2"], sp_map["O2-"], sp_map["O-"], sp_map["O2+"]
        rxns = []

        # 1. O2+ + e => O + O (TwoTemp)
        r1 = Reaction(reactants=[O2p, e], products=[O, O], reversible=False,
                      kinetics=TwoTemperaturePlasma(
                          A=(6.0e-5, "cm^3/(molecule*s)"), n=-1.0, Ea_g=(0, 'K'), Ea_e=(0, 'K')))
        rxns.append(r1)

        # 2. e + O2 + O2 => O2- + O2 (TwoTemp)
        r2 = Reaction(reactants=[e, O2, O2], products=[O2m, O2], reversible=False,
                      kinetics=TwoTemperaturePlasma(
                          A=(4.2e-27, "cm^6/(molecule^2*s)"), n=-1.0, Ea_g=(600, 'K'), Ea_e=(700, 'K')))
        rxns.append(r2)

        # 3. O- + O => O2 + e (Arrhenius)
        r3 = Reaction(reactants=[Om, O], products=[O2, e], reversible=False,
                      kinetics=Arrhenius(A=(5.0e-10, "cm^3/(molecule*s)"), n=0.0, Ea=(0, 'K'), T0=(1, 'K')))
        rxns.append(r3)

        # 4. O2- + O2 => O2 + e + O2
        r4 = Reaction(reactants=[O2m, O2], products=[O2, e, O2], reversible=False,
                      kinetics=Arrhenius(A=(1.559e-11, "cm^3/(molecule*s)"), n=0.5, Ea=(5590, 'K'), T0=(1, 'K')))
        rxns.append(r4)

        # 5. O2- + O => O2 + O-
        r5 = Reaction(reactants=[O2m, O], products=[O2, Om], reversible=False,
                      kinetics=Arrhenius(A=(3.3e-10, "cm^3/(molecule*s)"), n=0.0, Ea=(0, 'K'), T0=(1, 'K')))
        rxns.append(r5)

        # 6. O2- + O2+ => O2 + O2
        r6 = Reaction(reactants=[O2m, O2p], products=[O2, O2], reversible=False,
                      kinetics=Arrhenius(A=(3.464e-06, "cm^3/(molecule*s)"), n=-0.5, Ea=(0, 'K'), T0=(1, 'K')))
        rxns.append(r6)

        # 7. O- + O2+ => O + O2
        r7 = Reaction(reactants=[Om, O2p], products=[O, O2], reversible=False,
                      kinetics=Arrhenius(A=(3.464e-06, "cm^3/(molecule*s)"), n=-0.5, Ea=(0, 'K'), T0=(1, 'K')))
        rxns.append(r7)

        # 8. O2- + O2+ => O2 + O + O
        r8 = Reaction(reactants=[O2m, O2p], products=[O2, O, O], reversible=False,
                      kinetics=Arrhenius(A=(1.0e-07, "cm^3/(molecule*s)"), n=0.0, Ea=(0, 'K'), T0=(1, 'K')))
        rxns.append(r8)

        # 9. O2- + O2+ + M => O2 + O2 + M (ThreeBody)
        r9 = Reaction(reactants=[O2m, O2p], products=[O2, O2], reversible=False,
                      kinetics=ThirdBody(
                          arrheniusLow=Arrhenius(A=(3.118e-19, "cm^6/(molecule^2*s)"), n=-2.5, Ea=(0, 'K'),
                                                 T0=(1, 'K')),
                          efficiencies={O2.molecule[0]: 1.0, O.molecule[0]: 1.0}
                      ))
        rxns.append(r9)

        # 10. O2 + e => O- + O (Collision)
        energies_c1 = np.array([4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2,
                                5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.5, 6.6, 6.7, 6.8, 6.9,
                                7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5,
                                8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.8, 9.9])
        sigma_c1 = np.array([0.0, 8.8e-25, 2.64e-24, 4.4e-24, 7.04e-24, 9.68e-24, 1.32e-23,
                             1.76e-23, 2.2e-23, 2.9e-23, 3.61e-23, 4.49e-23, 5.37e-23, 6.33e-23, 7.48e-23,
                             8.53e-23, 9.59e-23, 1.05e-22, 1.14e-22, 1.23e-22, 1.31e-22, 1.36e-22, 1.41e-22,
                             1.4e-22, 1.37e-22, 1.34e-22, 1.28e-22, 1.22e-22, 1.14e-22, 1.06e-22, 9.85e-23,
                             8.97e-23, 8.18e-23, 7.39e-23, 6.42e-23, 5.72e-23, 5.01e-23, 4.49e-23, 3.87e-23,
                             3.34e-23, 2.82e-23, 2.38e-23, 2.02e-23, 1.67e-23, 1.41e-23, 1.23e-23, 1.06e-23,
                             8.8e-24, 7.04e-24, 7.04e-24, 6.16e-24, 5.28e-24, 4.4e-24, 4.4e-24, 3.52e-24, 3.52e-24])
        c1 = Reaction(reactants=[O2, e], products=[Om, O], reversible=False,
                      kinetics=ElectronCollisionPlasma(energies=(energies_c1, 'eV/molecule'), sigma=(sigma_c1, 'cm^2')))
        rxns.append(c1)

        # 11. O2 + e => e + e + O2+ (Collision)
        energies_c2 = np.array([13.0, 15.5, 18, 23])
        sigma_c2 = np.array([1.17e-22, 7.3e-22, 1.64e-21, 3.66e-21])
        c2 = Reaction(reactants=[O2, e], products=[e, e, O2p], reversible=False,
                      kinetics=ElectronCollisionPlasma(energies=(energies_c2, 'eV/molecule'), sigma=(sigma_c2, 'cm^2')))
        rxns.append(c2)

        # 12. O2 + e => e + O + O (Collision)
        energies_c3 = np.array([13.5, 18.5, 21.0, 23.5])
        sigma_c3 = np.array([2.2e-21, 5.29e-21, 5.65e-21, 5.25e-21])
        c3 = Reaction(reactants=[O2, e], products=[e, O, O], reversible=False,
                      kinetics=ElectronCollisionPlasma(energies=(energies_c3, 'eV/molecule'), sigma=(sigma_c3, 'cm^2')))
        rxns.append(c3)

        return rxns

    def test_replicate_oxygen_plasma(self):
        """
        1. Build RMG object model for Oxygen plasma.
        2. Write to Cantera YAML (chem.yaml).
        3. Load existing 'reference.yaml' from disk.
        4. Compare species, reactions, and rates.
        """
        # --- A. Build Model ---
        species_list = self._create_species_objects()
        sp_map = {s.label: s for s in species_list}
        reaction_list = self._create_reaction_objects(sp_map)

        # Mock RMG Structure
        class MockCore:
            def __init__(self, s, r): self.species, self.reactions = s, r

        class MockModel:
            def __init__(self, core): self.core, self.edge = core, MockCore([], [])

        class MockRMG:
            def __init__(self, out_dir, model):
                self.output_directory = out_dir
                self.reaction_model = model
                self.save_edge_species = False

        mock_rmg = MockRMG(self.tmp_dir, MockModel(MockCore(species_list, reaction_list)))

        # --- B. Write Generated File ---
        save_cantera_files(mock_rmg)
        gen_yaml_path = os.path.join(self.tmp_dir, 'cantera', 'chem.yaml')
        assert os.path.exists(gen_yaml_path)

        # --- C. Load Generated File in Cantera (Validation) ---
        try:
            sol_gen = ct.Solution(gen_yaml_path)
        except Exception as e:
            pytest.fail(f"Generated Cantera file is invalid: {e}")

        # Check basic properties of loaded solution
        assert sol_gen.n_species == 6
        assert sol_gen.n_reactions == 12
        assert sol_gen.transport_model == 'ionized-gas'

        # --- D. Compare Contents (YAML vs YAML) ---
        # We read the reference file manually since it can't be loaded by Cantera 3.0
        with open(self.ref_yaml_path, 'r') as f:
            ref_data = yaml.safe_load(f)

        with open(gen_yaml_path, 'r') as f:
            gen_data = yaml.safe_load(f)

        # 1. Compare Species counts
        assert len(gen_data['species']) == 6
        assert len(ref_data['species']) == 1  # only defines an electron, refers to an external file for the other species

        # 2. Compare Specific Reactions
        # We assume the order is preserved because we created them in that order.

        # Reaction 0: O2+ + e => O + O (TwoTemp)
        # Ref: rate-constant: {A: 6.0e-5, b: -1.0, Ea-gas: 0.0, Ea-electron: 0.0}
        r_gen = gen_data['reactions'][0]
        r_ref = ref_data['reactions'][0]

        assert r_gen['type'] == 'two-temperature-plasma'
        print(r_gen['rate-constant'], r_ref['rate-constant'])
        assert np.isclose(r_gen['rate-constant']['A'], r_ref['rate-constant']['A'] * constants.Na * 1e-6)
        assert r_gen['rate-constant']['b'] == r_ref['rate-constant']['b']

        # Reaction 8: O2- + O2+ + M => O2 + O2 + M (ThreeBody)
        # Ref efficiencies: {O2: 1.0, O: 1.0}
        r_gen = gen_data['reactions'][8]
        r_ref = ref_data['reactions'][8]

        assert r_gen['type'] == 'three-body'
        assert r_gen['efficiencies']['O2'] == r_ref['efficiencies']['O2']
        assert r_gen['efficiencies']['O'] == r_ref['efficiencies']['O']

        # Reaction 10: O2 + e => O- + O (ElectronCollision)
        # Compare a few cross-sections to ensure array integrity
        r_gen = gen_data['reactions'][9]  # Index 9 is the 10th reaction
        r_ref = ref_data['collisions'][0]  # Note: Ref file separates 'collisions', our writer puts them in 'reactions'

        assert r_gen['type'] == 'electron-collision-plasma'

        # Compare Energy Levels
        # Note: RMG Writer converts to float list, Ref file has them as list.
        gen_energies = np.array(r_gen['energy-levels'])
        ref_energies = np.array(r_ref['energy-levels'])
        np.testing.assert_allclose(gen_energies, ref_energies, rtol=1e-4, err_msg="Energy levels mismatch")

        # Compare Cross-Sections
        gen_sigma = np.array(r_gen['cross-sections'])
        ref_sigma = np.array(r_ref['cross-sections'])
        # Convert Reference (cm^2) to SI (m^2) for comparison
        # 1 cm^2 = 1e-4 m^2
        ref_sigma_si = ref_sigma * 1e-4
        np.testing.assert_allclose(gen_sigma, ref_sigma_si, rtol=1e-4, err_msg="Cross sections mismatch")
