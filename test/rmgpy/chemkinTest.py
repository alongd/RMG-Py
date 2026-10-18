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

import os

from unittest import mock

import rmgpy
from rmgpy.chemkin import (
    get_species_identifier,
    load_chemkin_file,
    load_transport_file,
    mark_duplicate_reactions,
    read_kinetics_entry,
    read_reaction_comments,
    read_thermo_entry,
    save_chemkin_file,
    save_chemkin_surface_file,
    save_species_dictionary,
    save_transport_file,
    read_reactions_block,
    write_kinetics_entry,
    write_thermo_entry,
)
from rmgpy.chemkin import _remove_line_breaks, _process_duplicate_reactions
from rmgpy.data.kinetics import LibraryReaction
from rmgpy.exceptions import ChemkinError
from rmgpy.kinetics.arrhenius import Arrhenius, MultiArrhenius
from rmgpy.kinetics.chebyshev import Chebyshev
from rmgpy.kinetics.falloff import Lindemann, ThirdBody, Troe
from rmgpy.kinetics.arrhenius import PDepArrhenius
from rmgpy import constants
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial
from rmgpy.transport import TransportData
from rmgpy.quantity import Quantity
from rmgpy.kinetics.surface import SurfaceArrhenius, StickingCoefficient
import pytest


class ChemkinTest:
    @mock.patch("rmgpy.chemkin.logging")
    def test_read_thermo_entry_bad_element_count(self, mock_logging):
        """
        Test that invalid element count logs the appropriate warning.

        This test uses the `mock` module in order to test calls to logging.
        The `mock.patch` decorator replaces the logging module instance in
        rmgpy.chemkin with a mock instance that can be accessed by this
        unit test. By using the mock instance, it is possible to assert that
        the expected logging statements are being created.
        """
        entry = """C2H6                    H   XC   X          L   100.000  5000.000  827.28      1
                2.44813916E+00 1.83377834E-02-7.25714119E-06 1.35300042E-09-9.60327447E-14    2
                -1.19655244E+04 8.07917520E+00 3.50507145E+00-3.65219841E-03 6.32200490E-05    3
                -8.01049582E-08 3.19734088E-11-1.15627878E+04 6.67152939E+00                   4
                """
        with pytest.raises(ValueError):
            read_thermo_entry(entry)

        mock_logging.info.assert_called_with(
            "Trouble reading line 'C2H6                    H   XC   X          L   100.000  5000.000  827.28      1' element segment 'H   X'"
        )

    @mock.patch("rmgpy.chemkin.logging")
    def test_read_thermo_entry_not_gas_phase(self, mock_logging):
        """
        Test that non gas phase data logs the appropriate warning.

        This test uses the `mock` module in order to test calls to logging.
        The `mock.patch` decorator replaces the logging module instance in
        rmgpy.chemkin with a mock instance that can be accessed by this
        unit test. By using the mock instance, it is possible to assert that
        the expected logging statements are being created.
        """
        entry = """C2H6                    H   6C   2          L   100.000  5000.000  827.28      1
                2.44813916E+00 1.83377834E-02-7.25714119E-06 1.35300042E-09-9.60327447E-14    2
                -1.19655244E+04 8.07917520E+00 3.50507145E+00-3.65219841E-03 6.32200490E-05    3
                -8.01049582E-08 3.19734088E-11-1.15627878E+04 6.67152939E+00                   4
                """
        species, thermo, formula = read_thermo_entry(entry)

        mock_logging.warning.assert_called_with("Was expecting gas phase thermo data for C2H6. Skipping thermo data.")
        assert species == "C2H6"
        assert formula is None
        assert thermo is None

    @mock.patch("rmgpy.chemkin.logging")
    def test_read_thermo_entry_not_float(self, mock_logging):
        """
        Test that non-float parameters log the appropriate warning.

        This test uses the `mock` module in order to test calls to logging.
        The `mock.patch` decorator replaces the logging module instance in
        rmgpy.chemkin with a mock instance that can be accessed by this
        unit test. By using the mock instance, it is possible to assert that
        the expected logging statements are being created.
        """
        entry = """C2H6                    H   6C   2          G   100.000  5000.000  827.28      1
 X.44813916E+00 1.83377834E-02-7.25714119E-06 1.35300042E-09-9.60327447E-14    2
-1.19655244E+04 8.07917520E+00 3.50507145E+00-3.65219841E-03 6.32200490E-05    3
-8.01049582E-08 3.19734088E-11-1.15627878E+04 6.67152939E+00                   4
"""
        species, thermo, formula = read_thermo_entry(entry)

        mock_logging.warning.assert_called_with("could not convert string to float: 'X.44813916E+00'")
        assert species == "C2H6"
        assert formula is None
        assert thermo is None

    def test_read_thermo_entry_no_temperature_range(self):
        """Test that missing temperature range can be handled for thermo entry."""
        entry = """C2H6                    H   6C   2          G                                  1
 2.44813916E+00 1.83377834E-02-7.25714119E-06 1.35300042E-09-9.60327447E-14    2
-1.19655244E+04 8.07917520E+00 3.50507145E+00-3.65219841E-03 6.32200490E-05    3
-8.01049582E-08 3.19734088E-11-1.15627878E+04 6.67152939E+00                   4
"""
        species, thermo, formula = read_thermo_entry(entry, Tmin=100.0, Tint=827.28, Tmax=5000.0)

        assert species == "C2H6"
        assert formula == {"H": 6, "C": 2}
        assert isinstance(thermo, NASA)

    def test_read_and_write_and_read_template_reaction_family_for_minimal_example(self):
        """
        This example tests if family and templates info can be correctly
        parsed from comments like '!Template reaction: R_Recombination'.
        """
        folder = os.path.join(os.path.dirname(__file__), "test_data", "chemkin", "chemkin_py")

        chemkin_path = os.path.join(folder, "minimal", "chem.inp")
        dictionary_path = os.path.join(folder, "minimal", "species_dictionary.txt")

        # read original chemkin file
        species, reactions = load_chemkin_file(chemkin_path, dictionary_path)

        # ensure correct reading
        reaction1 = reactions[0]
        assert reaction1.family == "R_Recombination"
        assert frozenset("C_methyl;C_methyl".split(";")) == frozenset(reaction1.template)
        reaction2 = reactions[1]
        assert reaction2.family == "H_Abstraction"
        assert frozenset(r"C/H3/Cs\H3;C_methyl".split(";")) == frozenset(reaction2.template)
        # save_chemkin_file
        chemkin_save_path = os.path.join(folder, "minimal", "chem_new.inp")
        dictionary_save_path = os.path.join(folder, "minimal", "species_dictionary_new.txt")

        save_chemkin_file(
            chemkin_save_path,
            species,
            reactions,
            verbose=True,
            check_for_duplicates=True,
        )
        save_species_dictionary(dictionary_save_path, species, old_style=False)

        assert os.path.isfile(chemkin_save_path)
        assert os.path.isfile(dictionary_save_path)

        # read newly written chemkin file to make sure the entire cycle works
        _, reactions2 = load_chemkin_file(chemkin_save_path, dictionary_save_path)

        reaction1_new = reactions2[0]
        assert reaction1_new.family == reaction1_new.family
        assert reaction1_new.template == reaction1_new.template
        assert reaction1_new.degeneracy == reaction1_new.degeneracy

        reaction2_new = reactions2[1]
        assert reaction2_new.family == reaction2_new.family
        assert reaction2_new.template == reaction2_new.template
        assert reaction2_new.degeneracy == reaction2_new.degeneracy

        # clean up
        os.remove(chemkin_save_path)
        os.remove(dictionary_save_path)

    def test_read_and_write_template_reaction_family_for_pdd_example(self):
        """
        This example is mainly to ensure comments like
        '! Kinetics were estimated in this direction instead
        of the reverse because:' or '! This direction matched
        an entry in H_Abstraction, the other was just an estimate.'
        won't interfere reaction family info retrival.
        """
        folder = os.path.join(os.path.dirname(__file__), "test_data", "chemkin", "chemkin_py")

        chemkin_path = os.path.join(folder, "pdd", "chem.inp")
        dictionary_path = os.path.join(folder, "pdd", "species_dictionary.txt")

        # load_chemkin_file
        species, reactions = load_chemkin_file(chemkin_path, dictionary_path)

        reaction1 = reactions[0]
        assert reaction1.family == "H_Abstraction"

        reaction2 = reactions[1]
        assert reaction2.family == "H_Abstraction"

        # save_chemkin_file
        chemkin_save_path = os.path.join(folder, "minimal", "chem_new.inp")
        dictionary_save_path = os.path.join(folder, "minimal", "species_dictionary_new.txt")

        save_chemkin_file(
            chemkin_save_path,
            species,
            reactions,
            verbose=False,
            check_for_duplicates=False,
        )
        save_species_dictionary(dictionary_save_path, species, old_style=False)

        assert os.path.isfile(chemkin_save_path)
        assert os.path.isfile(dictionary_save_path)

        # clean up
        os.remove(chemkin_save_path)
        os.remove(dictionary_save_path)

    def test_transport_data_read_and_write(self):
        """
        Test that we can write to chemkin and recreate the same transport object
        """
        Ar = Species(
            label="Ar",
            transport_data=TransportData(
                shapeIndex=0,
                epsilon=(1134.93, "J/mol"),
                sigma=(3.33, "angstrom"),
                dipoleMoment=(0, "De"),
                polarizability=(0, "angstrom^3"),
                rotrelaxcollnum=0.0,
                comment="""GRI-Mech""",
            ),
        )

        Ar_write = Species(label="Ar")
        folder = os.path.join(os.path.dirname(__file__), "test_data", "chemkin", "chemkin_py")

        temp_transport_path = os.path.join(folder, "tran_temp.dat")

        save_transport_file(temp_transport_path, [Ar])
        species_dict = {"Ar": Ar_write}
        load_transport_file(temp_transport_path, species_dict)
        assert repr(Ar) == repr(Ar_write)

        os.remove(temp_transport_path)

    def test_use_chemkin_names(self):
        """
        Test that the official chemkin names are used as labels for the created Species objects.
        """

        folder = os.path.join(os.path.dirname(__file__), "test_data", "chemkin", "chemkin_py")

        chemkin_path = os.path.join(folder, "minimal", "chem.inp")
        dictionary_path = os.path.join(folder, "minimal", "species_dictionary.txt")

        # load_chemkin_file
        species, reactions = load_chemkin_file(chemkin_path, dictionary_path, use_chemkin_names=True)

        expected = ["Ar", "He", "Ne", "N2", "ethane", "CH3", "C2H5", "C"]

        for spc, label in zip(species, expected):
            assert spc.label == label

    def test_reactant_n2_is_reactive_and_gets_right_species_identifier(self):
        """
        Test that after loading chemkin files, species such as N2, which is in the default
        inert list of RMG, should be treated as reactive species and given right species
        Identifier when it's reacting in reactions.
        """
        folder = os.path.join(os.path.dirname(__file__), "test_data", "chemkin", "chemkin_py")

        chemkin_path = os.path.join(folder, "NC", "chem.inp")
        dictionary_path = os.path.join(folder, "NC", "species_dictionary.txt")

        # load_chemkin_file
        species, reactions = load_chemkin_file(chemkin_path, dictionary_path, use_chemkin_names=True)

        for n2 in species:
            if n2.label == "N2":
                break
        assert n2.reactive

        assert get_species_identifier(n2) == "N2(35)"

    def test_read_specific_collider(self):
        """
        Test that a Chemkin reaction with a specific species as a third body collider can be properly read
        even if the species name contains parenthesis
        """
        entry = """O2(4)+H(5)(+N2(5))<=>HO2(10)(+N2(5))                          4.651e+12 0.440     0.000"""
        species_dict = {}
        s1 = Species().from_adjacency_list(
            """O2(4)
                                         multiplicity 3
                                         1 O u1 p2 c0 {2,S}
                                         2 O u1 p2 c0 {1,S}"""
        )
        s2 = Species().from_adjacency_list(
            """H(5)
                                         multiplicity 2
                                         1 H u1 p0 c0"""
        )
        s3 = Species().from_adjacency_list(
            """N2(5)
                                         1 N u0 p1 c0 {2,T}
                                         2 N u0 p1 c0 {1,T}"""
        )
        s4 = Species().from_adjacency_list(
            """HO2(10)
                                         multiplicity 2
                                         1 O u0 p2 c0 {2,S} {3,S}
                                         2 O u1 p2 c0 {1,S}
                                         3 H u0 p0 c0 {1,S}"""
        )
        species_dict["O2(4)"] = s1
        species_dict["H(5)"] = s2
        species_dict["N2(5)"] = s3
        species_dict["HO2(10)"] = s4
        A_units = ["", "s^-1", "cm^3/(mol*s)", "cm^6/(mol^2*s)", "cm^9/(mol^3*s)"]
        A_units_surf = []
        E_units = "kcal/mol"
        reaction = read_kinetics_entry(entry, species_dict, A_units, A_units_surf, E_units)

        assert reaction.specific_collider.label == "N2(5)"

    def test_process_duplicate_reactions(self):
        """
        Test that duplicate reactions are handled correctly when
        loading a Chemkin file.
        """
        s1 = Species().from_smiles("CC")
        s2 = Species().from_smiles("[CH3]")
        s3 = Species().from_smiles("[OH]")
        s4 = Species().from_smiles("C[CH2]")
        s5 = Species().from_smiles("O")
        r1 = Reaction(reactants=[s1], products=[s2, s2], duplicate=False, kinetics=Arrhenius())
        r2 = Reaction(reactants=[s1, s3], products=[s4, s5], duplicate=True, kinetics=Arrhenius())
        r3 = Reaction(reactants=[s1, s3], products=[s4, s5], duplicate=True, kinetics=Arrhenius())
        r4 = Reaction(reactants=[s1, s3], products=[s4, s5], duplicate=False, kinetics=Arrhenius())
        r5 = LibraryReaction(
            reactants=[s1, s3],
            products=[s4, s5],
            duplicate=True,
            kinetics=Arrhenius(),
            library="lib1",
        )
        r6 = LibraryReaction(
            reactants=[s1, s3],
            products=[s4, s5],
            duplicate=True,
            kinetics=Arrhenius(),
            library="lib2",
        )
        r7 = LibraryReaction(
            reactants=[s1, s3],
            products=[s4, s5],
            duplicate=True,
            kinetics=Chebyshev(),
            library="lib1",
        )
        r8 = LibraryReaction(
            reactants=[s1, s3],
            products=[s4, s5],
            duplicate=True,
            kinetics=Arrhenius(),
            library="lib1",
        )
        r9 = LibraryReaction(
            reactants=[s1, s3],
            products=[s4, s5],
            duplicate=False,
            kinetics=MultiArrhenius(arrhenius=[Arrhenius(), Arrhenius()]),
            library="lib1",
        )
        reaction_list_with_duplicate = [r1, r2, r3]
        reaction_list_with_duplicate2 = [r1, r2, r3]
        reaction_list_unmarked_duplicate = [r1, r2, r4]
        reaction_list_unequal_libraries = [r1, r5, r6]
        reaction_list_mixed_kinetics = [r1, r5, r7]
        reaction_list_mergeable = [r1, r5, r8]
        reaction_list_merged = [r1, r9]

        # Test that duplicates are not removed for non-library reactions
        _process_duplicate_reactions(reaction_list_with_duplicate)
        assert reaction_list_with_duplicate == reaction_list_with_duplicate2

        # Test that unmarked duplicate reactions are detected if both
        # reactions are p-dep or p-indep
        with pytest.raises(ChemkinError, match="Encountered unmarked duplicate reaction"):
            _process_duplicate_reactions(
                reaction_list_unmarked_duplicate,
            )

        # Test that unequal libraries are recognized
        with pytest.raises(ChemkinError, match="from different libraries"):
            _process_duplicate_reactions(
                reaction_list_unequal_libraries,
            )

        # Test that an error is raised for reactions with kinetics
        # that cannot be merged
        with pytest.raises(ChemkinError, match="Mixed kinetics for duplicate reaction"):
            _process_duplicate_reactions(
                reaction_list_mixed_kinetics,
            )

        # Test that duplicate library reactions are merged successfully
        _process_duplicate_reactions(reaction_list_mergeable)
        assert len(reaction_list_mergeable) == len(reaction_list_merged)
        assert reaction_list_mergeable[0] == reaction_list_merged[0]
        rtest = reaction_list_mergeable[1]
        rtrue = reaction_list_merged[1]
        assert rtest.reactants == rtrue.reactants
        assert rtest.products == rtrue.products
        assert rtest.duplicate == rtrue.duplicate
        assert rtest.library == rtrue.library
        assert isinstance(rtest.kinetics, MultiArrhenius)
        assert all(isinstance(k, Arrhenius) for k in rtest.kinetics.arrhenius)

    def test_mark_duplicate_reactions(self):
        """Test that we can properly mark duplicate reactions for Chemkin."""
        s1 = Species().from_smiles("CC")
        s2 = Species().from_smiles("[CH3]")
        s3 = Species().from_smiles("[OH]")
        s4 = Species().from_smiles("C[CH2]")
        s5 = Species().from_smiles("O")
        s6 = Species().from_smiles("[H]")

        # Try initializing with duplicate=False
        reaction_list = [
            Reaction(reactants=[s1], products=[s2, s2], duplicate=False, kinetics=Arrhenius()),
            Reaction(reactants=[s1], products=[s2, s2], duplicate=False, kinetics=Arrhenius()),
            Reaction(
                reactants=[s1, s3],
                products=[s4, s5],
                duplicate=False,
                kinetics=Arrhenius(),
            ),
            Reaction(
                reactants=[s1, s3],
                products=[s4, s5],
                duplicate=False,
                kinetics=Chebyshev(),
            ),
            Reaction(
                reactants=[s1],
                products=[s4, s6],
                duplicate=False,
                kinetics=Arrhenius(),
                reversible=False,
            ),
            Reaction(
                reactants=[s1],
                products=[s4, s6],
                duplicate=False,
                kinetics=Arrhenius(),
                reversible=False,
            ),
            Reaction(
                reactants=[s5],
                products=[s3, s6],
                duplicate=False,
                kinetics=Arrhenius(),
                reversible=False,
            ),
            Reaction(
                reactants=[s3, s6],
                products=[s5],
                duplicate=False,
                kinetics=Arrhenius(),
                reversible=False,
            ),
        ]

        expected_flags = [True, True, False, False, True, True, False, False]

        mark_duplicate_reactions(reaction_list)
        duplicate_flags = [rxn.duplicate for rxn in reaction_list]

        assert duplicate_flags == expected_flags

        # Try initializing with duplicate=True
        reaction_list = [
            Reaction(reactants=[s1], products=[s2, s2], duplicate=True, kinetics=Arrhenius()),
            Reaction(reactants=[s1], products=[s2, s2], duplicate=True, kinetics=Arrhenius()),
            Reaction(
                reactants=[s1, s3],
                products=[s4, s5],
                duplicate=True,
                kinetics=Arrhenius(),
            ),
            Reaction(
                reactants=[s1, s3],
                products=[s4, s5],
                duplicate=True,
                kinetics=Chebyshev(),
            ),
            Reaction(
                reactants=[s1],
                products=[s4, s6],
                duplicate=True,
                kinetics=Arrhenius(),
                reversible=False,
            ),
            Reaction(
                reactants=[s1],
                products=[s4, s6],
                duplicate=True,
                kinetics=Arrhenius(),
                reversible=False,
            ),
            Reaction(
                reactants=[s5],
                products=[s3, s6],
                duplicate=True,
                kinetics=Arrhenius(),
                reversible=False,
            ),
            Reaction(
                reactants=[s3, s6],
                products=[s5],
                duplicate=True,
                kinetics=Arrhenius(),
                reversible=False,
            ),
        ]

        mark_duplicate_reactions(reaction_list)
        duplicate_flags = [rxn.duplicate for rxn in reaction_list]

        assert duplicate_flags == expected_flags

    def test_read_write_coverage_dependence(self):
        """Test that we can properly read and write coverage dependent parameters"""

        folder = os.path.join(os.path.dirname(__file__), "test_data", "chemkin", "chemkin_py")

        s_x_entry = """X                       X   1               G   100.000  5000.000 1554.80      1
 1.60299900E-01-2.52235409E-04 1.14181275E-07-1.21471653E-11 3.85790025E-16    2
-7.08100885E+01-9.09527530E-01 7.10139498E-03-4.25619522E-05 8.98533016E-08    3
-7.80193649E-11 2.32465471E-14-8.76101712E-01-3.11211229E-02                   4"""
        s_hx_entry = """H*                      H   1X   1          G   100.000  5000.000  952.91      1
 2.80339655E+00-5.41047017E-04 4.99507978E-07-7.54963647E-11 3.06772366E-15    2
-2.34636021E+03-1.59436787E+01-3.80965452E-01 5.47228709E-03 2.60912778E-06    3
-9.64961980E-09 4.63946753E-12-1.40561079E+03 1.01725550E+00                   4"""
        s_h2_entry = """H2                      H   2               G   100.000  5000.000 1959.08      1
 2.78816619E+00 5.87640475E-04 1.59010635E-07-5.52739465E-11 4.34311304E-15    2
-5.96144481E+02 1.12730527E-01 3.43536411E+00 2.12710383E-04-2.78625110E-07    3
 3.40267219E-10-7.76032129E-14-1.03135984E+03-3.90841731E+00                   4"""
        s_ohx_entry = """OH*                     H   1O   1X   1     G   100.000  5000.000  914.54      1
 2.43541504E+00 4.64599933E-03-2.39987608E-06 4.26351871E-10-2.60607840E-14    2
-2.17972457E+04-1.03879495E+01-1.29522792E+00 3.36487597E-02-7.07760603E-05    3
 6.54375105E-08-2.19437885E-11-2.16453879E+04 4.37657050E+00                   4"""
        s_ox_entry = """O*                      O   1X   1          G   100.000  5000.000  888.26      1
 1.89893631E+00 2.03295404E-03-1.19976562E-06 2.32680628E-10-1.53508256E-14    2
-2.21111565E+04-9.64104147E+00-7.59013115E-01 1.89868504E-02-3.82473769E-05    3
 3.43558429E-08-1.13974388E-11-2.18356104E+04 1.76017413E+00                   4"""

        s_h2 = Species().from_smiles("[H][H]")
        s_x = Species().from_adjacency_list("1 X u0 p0")
        s_x.label = "X"
        s_hx = Species().from_adjacency_list("1 H u0 p0 {2,S} \n 2 X u0 p0 {1,S}")
        s_hx.label = "HX"
        s_ohx = Species().from_smiles("O[*]")
        s_ohx.label = "OHX"
        s_ox = Species().from_smiles("O=[*]")
        s_ox.label = "OX"

        species, thermo, formula = read_thermo_entry(s_x_entry)
        s_x.thermo = thermo
        species, thermo, formula = read_thermo_entry(s_hx_entry)
        s_hx.thermo = thermo
        species, thermo, formula = read_thermo_entry(s_h2_entry)
        s_h2.thermo = thermo
        species, thermo, formula = read_thermo_entry(s_ohx_entry)
        s_ohx.thermo = thermo
        species, thermo, formula = read_thermo_entry(s_ox_entry)
        s_ox.thermo = thermo

        species = [s_h2, s_x, s_hx, s_ohx, s_ox]

        # test a coverage dependent arrhenius and a coverage dependent sticking coefficient
        rxn1_entry = """H2+X+X=H*+H*                                        4.000e-02 0.000     7.098
    STICK
    COV / H*                                        0         0         1.099     /"""

        rxn2_entry = """X+OH*=O*+H*                                         7.390000e+19 0.000     18.475   
    COV / O*                                        0         0         -17.500   /"""

        species_dict = {}
        species_dict["H2"] = s_h2
        species_dict["X"] = s_x
        species_dict["H*"] = s_hx
        species_dict["OH*"] = s_ohx
        species_dict["O*"] = s_ox

        A_units = ["", "s^-1", "cm^3/(mol*s)", "cm^6/(mol^2*s)", "cm^9/(mol^3*s)"]
        A_units_surf = ["", "s^-1", "cm^2/(mol*s)", "cm^4/(mol^2*s)", "cm^6/(mol^3*s)"]

        E_units = "kcal/mol"
        self.rxn1 = read_kinetics_entry(rxn1_entry, species_dict, A_units, A_units_surf, E_units)

        self.rxn2 = read_kinetics_entry(rxn2_entry, species_dict, A_units, A_units_surf, E_units)

        reactions = [self.rxn1, self.rxn2]

        # save_chemkin_file
        chemkin_save_path = os.path.join(folder, "surface", "chem-surface_new.inp")
        dictionary_save_path = os.path.join(folder, "surface", "species_dictionary_new.txt")
        save_chemkin_file(
            chemkin_save_path,
            species,
            reactions,
            verbose=True,
            check_for_duplicates=True,
        )
        save_species_dictionary(dictionary_save_path, species, old_style=False)

        # load chemkin file
        species_load, reactions_load = load_chemkin_file(chemkin_save_path, dictionary_save_path)

        # compare only labels because the objects are different
        species_label = []
        species_load_label = []
        for s in range(len(species)):
            species_label.append(species[s].label)
            species_load_label.append(species_load[s].label)

        # check the written chemkin file matches what is expected
        assert all(i in species_load_label for i in species_label)

        for r in range(len(reactions)):
            for s in range(len(reactions[r].products)):
                assert reactions[r].products[s].label == reactions_load[r].products[s].label
            for s in range(len(reactions[r].reactants)):
                assert reactions[r].reactants[s].label == reactions_load[r].reactants[s].label
            assert reactions[r].kinetics.coverage_dependence is not None
            assert reactions_load[r].kinetics.coverage_dependence is not None

        # clean up
        os.remove(chemkin_save_path)
        os.remove(dictionary_save_path)


class TestThermoReadWrite:
    def setup_class(self):
        """This method is run once before each test."""
        coeffs_low = [
            4.03055,
            -0.00214171,
            4.90611e-05,
            -5.99027e-08,
            2.38945e-11,
            -11257.6,
            3.5613,
        ]
        coeffs_high = [
            -0.307954,
            0.0245269,
            -1.2413e-05,
            3.07724e-09,
            -3.01467e-13,
            -10693,
            22.628,
        ]
        Tmin = 300.0
        Tmax = 3000.0
        Tint = 650.73
        E0 = -782292.0  # J/mol.
        comment = "C2H6"
        self.nasa = NASA(
            polynomials=[
                NASAPolynomial(coeffs=coeffs_low, Tmin=(Tmin, "K"), Tmax=(Tint, "K")),
                NASAPolynomial(coeffs=coeffs_high, Tmin=(Tint, "K"), Tmax=(Tmax, "K")),
            ],
            Tmin=(Tmin, "K"),
            Tmax=(Tmax, "K"),
            E0=(E0, "J/mol"),
            comment=comment,
        )

        # Chemkin entries for testing - note that the values are all the same
        self.entry1 = """C2H6                    C   2H   6          G   300.000  3000.000  650.73      1
-3.07954000E-01 2.45269000E-02-1.24130000E-05 3.07724000E-09-3.01467000E-13    2
-1.06930000E+04 2.26280000E+01 4.03055000E+00-2.14171000E-03 4.90611000E-05    3
-5.99027000E-08 2.38945000E-11-1.12576000E+04 3.56130000E+00                   4
"""

        self.entry2 = """CH3NO2X                                     G   300.000  3000.000  650.73      1&
C 1 H 3 N 1 O 2 X 1
-3.07954000E-01 2.45269000E-02-1.24130000E-05 3.07724000E-09-3.01467000E-13    2
-1.06930000E+04 2.26280000E+01 4.03055000E+00-2.14171000E-03 4.90611000E-05    3
-5.99027000E-08 2.38945000E-11-1.12576000E+04 3.56130000E+00                   4
"""

        self.entry3 = """CH3NO2SX                                    G   300.000  3000.000  650.73      1&
C 1 H 3 N 1 O 2 S 1 X 1
-3.07954000E-01 2.45269000E-02-1.24130000E-05 3.07724000E-09-3.01467000E-13    2
-1.06930000E+04 2.26280000E+01 4.03055000E+00-2.14171000E-03 4.90611000E-05    3
-5.99027000E-08 2.38945000E-11-1.12576000E+04 3.56130000E+00                   4
"""

    def test_write_thermo_block(self):
        """Test that we can write a normal thermo block"""
        species = Species(smiles="CC")
        species.thermo = self.nasa

        result = write_thermo_entry(species, verbose=False)

        assert result == self.entry1

    def test_read_thermo_block(self):
        """Test that we can read a normal thermo block"""
        species, thermo, formula = read_thermo_entry(self.entry1)

        assert species == "C2H6"
        assert formula == {"H": 6, "C": 2}
        assert self.nasa.is_identical_to(thermo)

    def test_write_thermo_block_5_elem(self):
        """Test that we can write a thermo block for a species with 5 elements"""
        species = Species().from_adjacency_list(
            """
1 O u0 p3 c-1 {3,S}
2 O u0 p2 c0 {3,D}
3 N u0 p0 c+1 {1,S} {2,D} {4,S}
4 C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 X u0 p0 c0
"""
        )
        species.thermo = self.nasa

        result = write_thermo_entry(species, verbose=False)

        assert result == self.entry2

    def test_read_thermo_block_5_elem(self):
        """Test that we can read a thermo block with 5 elements"""
        species, thermo, formula = read_thermo_entry(self.entry2)

        assert species == "CH3NO2X"
        assert formula == {"X": 1, "C": 1, "O": 2, "H": 3, "N": 1}
        assert self.nasa.is_identical_to(thermo)

    def test_write_thermo_block_6_elem(self):
        """Test that we can write a thermo block for a species with 6 elements"""
        species = Species().from_adjacency_list(
            """
1 O u0 p3 c-1 {2,S}
2 N u0 p0 c+1 {1,S} {3,D} {4,S}
3 O u0 p2 c0 {2,D}
4 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
5 S u0 p2 c0 {4,S} {8,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
9 X u0 p0 c0
"""
        )
        species.thermo = self.nasa

        result = write_thermo_entry(species, verbose=False)

        assert result == self.entry3

    def test_read_thermo_block_6_elem(self):
        """Test that we can read a thermo block with 6 elements"""
        species, thermo, formula = read_thermo_entry(self.entry3)

        assert species == "CH3NO2SX"
        assert formula == {"X": 1, "C": 1, "O": 2, "H": 3, "N": 1, "S": 1}
        assert self.nasa.is_identical_to(thermo)

    def test_write_bidentate_species(self):
        """Test that species with 2 or more surface sites get proper formatting"""

        folder = os.path.join(os.path.dirname(__file__), "test_data", "chemkin", "chemkin_py")
        chemkin_path = os.path.join(folder, "surface", "chem-surface.inp")
        dictionary_path = os.path.join(folder, "surface", "species_dictionary.txt")
        chemkin_save_path = os.path.join(folder, "surface", "chem-surface-test.inp")
        species, reactions = load_chemkin_file(chemkin_path, dictionary_path, use_chemkin_names=True)
        chox3 = next(iter(s for s in species if s.label=="CHOX3"))
        surface_atom_count = chox3.molecule[0].get_num_atoms("X")
        assert surface_atom_count == 3
        save_chemkin_surface_file(
            chemkin_save_path,
            species,
            reactions,
            verbose=False,
            check_for_duplicates=False,
        )

        bidentate_test = "CH2OX2/2/"
        tridentate_test = "CHOX3/3/"
        with open(chemkin_save_path, "r") as f:
            lines = [line.strip() for line in f]
        assert bidentate_test in lines
        assert tridentate_test in lines

        os.remove(chemkin_save_path)


class TestReadReactionComments:
    @classmethod
    def setup_class(cls):
        r = Species().from_smiles("[CH3]")
        r.label = "[CH3]"
        p = Species().from_smiles("CC")
        p.label = "CC"

        cls.reaction = Reaction(
            reactants=[r, r],
            products=[p],
            kinetics=Arrhenius(A=(8.26e17, "cm^3/(mol*s)"), n=-1.4, Ea=(1, "kcal/mol"), T0=(1, "K")),
        )
        cls.comments_list = [
            """
                              Reaction index: Chemkin #1; RMG #1
                              Template reaction: R_Recombination
                              Exact match found for rate rule (C_methyl;C_methyl)
                              Multiplied by reaction path degeneracy 0.5
                              """,
            """
                             Reaction index: Chemkin #2; RMG #4
                             Template reaction: H_Abstraction
                             Estimated using template (C/H3/Cs;C_methyl) for rate rule (C/H3/Cs\\H3;C_methyl)
                             Multiplied by reaction path degeneracy 6
                             """,
            """
                              Reaction index: Chemkin #13; RMG #8
                              Template reaction: H_Abstraction
                              Flux pairs: [CH3], CC; [CH3], CC; 
                              Estimated using an average for rate rule [C/H3/Cs\\H3;C_rad/H2/Cs]
                              Multiplied by reaction path degeneracy 6.0
                              """,
            """
                              Reaction index: Chemkin #17; RMG #31
                              Template reaction: H_Abstraction
                              Flux pairs: [CH3], CC; [CH3], CC; 
                              Estimated using average of templates [C/H3/Cs;H_rad] + [C/H3/Cs\\H3;Y_rad] for rate rule [C/H3/Cs\\H3;H_rad]
                              Multiplied by reaction path degeneracy 6.0
                              """,
            """
                              Reaction index: Chemkin #69; RMG #171
                              Template reaction: intra_H_migration
                              Flux pairs: [CH3], CC; [CH3], CC; 
                              Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] for rate rule
                              [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
                              Multiplied by reaction path degeneracy 3.0
                              """,
            """
                              Reaction index: Chemkin #3; RMG #243
                              Template reaction: Disproportionation
                              Flux pairs: [CH3], CC; [CH3], CC; 
                              Average of [Average of [O2b;O_Csrad] + Average of [O_atom_triplet;O_Csrad + CH2_triplet;O_Csrad] + Average of [Average of [Ct_rad/Ct;O_Csrad from
                              training reaction 0] + Average of [O_pri_rad;O_Csrad + Average of [O_rad/NonDeC;O_Csrad + O_rad/NonDeO;O_Csrad]] + Average of [Cd_pri_rad;O_Csrad] +
                              Average of [CO_pri_rad;O_Csrad] + Average of [C_methyl;O_Csrad + Average of [C_rad/H2/Cs;O_Csrad + C_rad/H2/Cd;O_Csrad + C_rad/H2/O;O_Csrad] + Average
                              of [C_rad/H/NonDeC;O_Csrad] + Average of [Average of [C_rad/Cs3;O_Csrad]]] + H_rad;O_Csrad]]
                              Estimated using template [Y_rad_birad_trirad_quadrad;O_Csrad] for rate rule [CH_quartet;O_Csrad]
                              """,
            """
                              Reaction index: Chemkin #4; RMG #303
                              Template reaction: Disproportionation
                              Flux pairs: [CH3], CC; [CH3], CC; 
                              Matched reaction 0 C2H + CH3O <=> C2H2 + CH2O in Disproportionation/training
                              """,
            """
                              Reaction index: Chemkin #51; RMG #136
                              Template reaction: H_Abstraction
                              Flux pairs: [CH3], CC; [CH3], CC; 
                              Estimated using an average for rate rule [C/H3/Cd\\H_Cd\\H2;C_rad/H2/Cs]
                              Euclidian distance = 0
                              Multiplied by reaction path degeneracy 3.0
                              """,
            """
                             Reaction index: Chemkin #32; RMG #27
                             Template reaction: R_Recombination
                             Flux pairs: [CH3], CC; [CH3], CC; 
                             Matched reaction 20 CH3 + CH3 <=> C2H6 in R_Recombination/training
                             This reaction matched rate rule [C_methyl;C_methyl]
                             """,
            """
                             Reaction index: Chemkin #2; RMG #4
                             Template reaction: R_Recombination
                             Flux pairs: [CH3], CC; [CH3], CC; 
                             From training reaction 21 used for C_rad/H2/Cs;C_methyl
                             Exact match found for rate rule [C_rad/H2/Cs;C_methyl]
                             Euclidian distance = 0
                             """,
        ]
        cls.template_list = [
            ["C_methyl", "C_methyl"],
            [r"C/H3/Cs\H3", "C_methyl"],
            [r"C/H3/Cs\H3", "C_rad/H2/Cs"],
            [r"C/H3/Cs\H3", "H_rad"],
            ["R3H_SS_Cs", "O_rad_out", "Cs_H_out_2H"],
            ["CH_quartet", "O_Csrad"],
            None,
            [r"C/H3/Cd\H_Cd\H2", "C_rad/H2/Cs"],
            ["C_methyl", "C_methyl"],
            ["C_rad/H2/Cs", "C_methyl"],
        ]
        cls.family_list = [
            "R_Recombination",
            "H_Abstraction",
            "H_Abstraction",
            "H_Abstraction",
            "intra_H_migration",
            "Disproportionation",
            "Disproportionation",
            "H_Abstraction",
            "R_Recombination",
            "R_Recombination",
        ]
        cls.degeneracy_list = [0.5, 6, 6, 6, 3, 1, 1, 3, 1, 1]
        cls.expected_lines = [4, 4, 5, 5, 5, 5, 4, 6, 5, 6]

    def test_read_reaction_comments_template(self):
        """
        Test that the template is picked up from reading reaction comments.
        """
        for index, comment in enumerate(self.comments_list):
            new_rxn = read_reaction_comments(self.reaction, comment)

            # only check template if meant to find one
            if self.template_list[index]:
                assert new_rxn.template, "The template was not saved from the reaction comment {}".format(comment)
                assert frozenset(new_rxn.template) == frozenset(self.template_list[index]), "The reaction template does not match"
            else:
                assert not new_rxn.template

    def test_read_reaction_comments_family(self):
        """
        Test that the family is picked up from reading reaction comments.
        """
        for index, comment in enumerate(self.comments_list):
            new_rxn = read_reaction_comments(self.reaction, comment)

            assert new_rxn.family == self.family_list[index], "wrong reaction family stored"

    def test_read_reaction_comments_degeneracy(self):
        """
        Test that the degeneracy is picked up from reading reaction comments.

        Also checks that reaction rate was not modified in the process.
        """
        for index, comment in enumerate(self.comments_list):
            # Clear any leftover kinetics comments
            self.reaction.kinetics.comment = ""
            previous_rate = self.reaction.kinetics.A.value_si
            new_rxn = read_reaction_comments(self.reaction, comment)
            new_rate = new_rxn.kinetics.A.value_si

            assert new_rxn.degeneracy == self.degeneracy_list[index], "wrong degeneracy was stored"
            assert previous_rate == new_rate

            # Check that the comment only appears once in the kinetics comment
            if new_rxn.degeneracy != 1:
                assert (
                    new_rxn.kinetics.comment.count("Multiplied by reaction path degeneracy {}".format(new_rxn.degeneracy)) == 1
                ), "Reaction degeneracy comment duplicated while reading Chemkin comments"
            else:
                assert "Multiplied by reaction path degeneracy" not in new_rxn.kinetics.comment

    def test_remove_line_breaks(self):
        """
        tests that _remove_line_breaks functions properly
        """
        for index, comment in enumerate(self.comments_list):
            new_comment = _remove_line_breaks(comment)
            new_comment_lines = len(new_comment.strip().splitlines())
            assert (
                new_comment_lines == self.expected_lines[index]
            ), 'Found {} more lines than expected for comment \n\n""{}""\n\n which converted to \n\n""{}""'.format(
                new_comment_lines - self.expected_lines[index],
                comment.strip(),
                new_comment.strip(),
            )


class TestLongLabelSpeciesIdentifier:
    """Regression for the >16-char moment-dummy Chemkin identifier overflow.

    Moment dummies are created with index = -1 (rmgpy/rmg/model.py), so
    get_species_identifier used to return the raw label even when longer than
    16 characters, overflowing the fixed-column THERMO format and breaking
    load_chemkin_file on the EPDM artifact (ValueError: could not convert
    string to float: 'E+0').
    """

    @staticmethod
    def _mu_dummy(label):
        from rmgpy.molecule import Molecule
        spc = Species(label=label, reactive=False)
        spc.molecule = [Molecule().from_smiles("[Ne]")]
        spc.index = -1
        spc.is_moment_dummy = True
        return spc

    @staticmethod
    def _ne_thermo():
        poly_low = NASAPolynomial(
            coeffs=[2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 3.35532],
            Tmin=(200.0, "K"), Tmax=(1000.0, "K"))
        poly_high = NASAPolynomial(
            coeffs=[2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 3.35532],
            Tmin=(1000.0, "K"), Tmax=(6000.0, "K"))
        return NASA(polynomials=[poly_low, poly_high],
                    Tmin=(200.0, "K"), Tmax=(6000.0, "K"))

    def test_long_label_is_compressed_to_16_chars(self):
        spc = self._mu_dummy("epdm_scission_tail_mu0")
        ident = get_species_identifier(spc)
        assert len(ident) <= 16
        assert ident.endswith("_mu0")  # readable moment tail preserved

    def test_short_labels_unchanged(self):
        spc = self._mu_dummy("epdm_mu0")
        assert get_species_identifier(spc) == "epdm_mu0"

    def test_long_non_polymer_label_roundtrips_through_reaction_block(self, tmpdir):
        """Long ordinary labels must agree between the reaction writer and reader."""
        from rmgpy.molecule import Molecule

        reactant = Species(label="hydroperoxylvinoxy", reactive=True)
        reactant.molecule = [Molecule().from_smiles("[CH2]C(=O)OO")]
        product = Species(label="acetylperoxy", reactive=True)
        product.molecule = [Molecule().from_smiles("CC(=O)O[O]")]
        reaction = Reaction(
            reactants=[reactant], products=[product],
            kinetics=Arrhenius(A=(1.0, "s^-1"), n=0, Ea=(0, "kcal/mol")),
        )
        chem_path = os.path.join(str(tmpdir), "chem.inp")
        with open(chem_path, "w") as chemkin_file:
            chemkin_file.write("REACTIONS\n")
            chemkin_file.write(write_kinetics_entry(reaction, [reactant, product], verbose=False))
            chemkin_file.write("END\n")

        with open(chem_path) as chemkin_file:
            loaded_reactions = read_reactions_block(
                chemkin_file, {reactant.label: reactant, product.label: product}
            )
        assert loaded_reactions[0].reactants == [reactant]

    def test_distinct_long_labels_get_distinct_identifiers(self):
        labels = ["epdm_scission_tail_mu0", "epdm_scission_head_mu0",
                  "epdm_scission_tail_mu1", "epdm_scission_tail_mu2"]
        idents = [get_species_identifier(self._mu_dummy(l)) for l in labels]
        assert len(set(idents)) == len(labels)
        assert all(len(i) <= 16 for i in idents)

    def test_identifier_is_deterministic(self):
        a = get_species_identifier(self._mu_dummy("epdm_scission_tail_mu0"))
        b = get_species_identifier(self._mu_dummy("epdm_scission_tail_mu0"))
        assert a == b

    def test_multi_digit_moment_tail_preserved(self):
        """_mu<digits> regex must match multi-digit moments (e.g. _mu10)."""
        spc = self._mu_dummy("epdm_scission_tail_mu10")
        ident = get_species_identifier(spc)
        assert len(ident) <= 16
        assert ident.endswith("_mu10"), (
            "Multi-digit moment tail '_mu10' should be preserved in compressed "
            "identifier; got {!r}".format(ident)
        )

    def test_near_collision_pair_distinct_with_6hex(self):
        """Widening the digest from 4 to 6 hex chars must resolve the
        reviewer-reported collision between epdm_scission_70_mu0 and
        epdm_scission_306_mu0 (both mapped to 'epdm_sc-ae04_mu0' at 4 hex).
        With 6 hex the digests are ae04c7 vs ae0420 — distinct."""
        spc_70 = self._mu_dummy("epdm_scission_70_mu0")
        spc_306 = self._mu_dummy("epdm_scission_306_mu0")
        ident_70 = get_species_identifier(spc_70)
        ident_306 = get_species_identifier(spc_306)
        assert ident_70 != ident_306, (
            "Near-collision pair should get distinct identifiers with 6-hex "
            "digest; both got {!r}".format(ident_70)
        )
        assert len(ident_70) <= 16
        assert len(ident_306) <= 16

    def test_chemkin_roundtrip_with_long_moment_dummy_labels(self, tmpdir):
        """save_chemkin_file + save_species_dictionary then load_chemkin_file
        must succeed with >16-char moment-dummy labels (EPDM failure mode)."""
        long_dummy = self._mu_dummy("epdm_scission_tail_mu0")
        long_dummy.thermo = self._ne_thermo()
        short_dummy = self._mu_dummy("epdm_mu0")
        short_dummy.thermo = self._ne_thermo()
        species = [long_dummy, short_dummy]

        chem_path = os.path.join(str(tmpdir), "chem.inp")
        dict_path = os.path.join(str(tmpdir), "species_dictionary.txt")
        save_chemkin_file(chem_path, species, [], verbose=False, check_for_duplicates=False)
        save_species_dictionary(dict_path, species)

        loaded_species, loaded_reactions = load_chemkin_file(chem_path, dict_path)
        assert len(loaded_species) == 2
        labels = {s.label for s in loaded_species}
        # the compressed identifier round-trips as the species label
        assert any(l.endswith("_mu0") and len(l) <= 16 for l in labels)


class TestEndKeyword:
    """A Chemkin block ends only on a line whose first token is END (any case,
    comments removed); a species identifier that merely holds 'end'
    (``blend(3)``, ``end_cap(4)``) keeps its reactions and its thermo."""

    THERM_ENTRY = [
        "{0:<24}C   1H   4          G   200.000  6000.000 1000.00      1",
        " 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2",
        "-7.45375000E+02 3.35532000E+00 2.50000000E+00 0.00000000E+00 0.00000000E+00    3",
        " 0.00000000E+00 0.00000000E+00-7.45375000E+02 3.35532000E+00                   4",
    ]
    NAMES = {"blend(3)": "CC", "end_cap(4)": "C", "H(1)": "[H]"}
    REACTIONS = ["blend(3)=>end_cap(4)+end_cap(4)    1.0 0.0 0.0",
                 "end_cap(4)+H(1)=>blend(3)    2.0 0.0 0.0",
                 "H(1)+H(1)=>blend(3)    3.0 0.0 0.0"]
    # each closer, then a line that the block would read had it not closed
    CLOSERS = {"lowercase": ["end"], "comment_after": ["END ! closes the block"],
               "comment_before": ["! a comment line", "  END"]}

    def load(self, tmp_path, closer):
        therm = ["THERM ALL", "   300.000  1000.000  5000.000"]
        for name in self.NAMES:
            therm += [line.format(name) for line in self.THERM_ENTRY]
        text = "\n".join(["ELEMENTS H C END", "SPECIES", *self.NAMES, "END", *therm, *closer,
                          "REACTIONS    KCAL/MOLE   MOLES", *self.REACTIONS, *closer,
                          "H(1)+blend(3)=>end_cap(4)+end_cap(4)+H(1)    9.0 0.0 0.0", ""])
        (tmp_path / "chem.inp").write_text(text)
        (tmp_path / "dict.txt").write_text("\n".join(f"{n}\n{Molecule(smiles=s).to_adjacency_list()}"
                                                     for n, s in self.NAMES.items()))
        return load_chemkin_file(str(tmp_path / "chem.inp"), str(tmp_path / "dict.txt"),
                                 check_duplicates=False, generate_resonance_structures=False)

    @pytest.mark.parametrize("closer", sorted(CLOSERS))
    def test_every_reaction_and_thermo_entry_survives_and_end_still_closes(self, tmp_path, closer):
        species, reactions = self.load(tmp_path, self.CLOSERS[closer])
        assert [str(r) for r in reactions] == ["blend(3) => end_cap(4) + end_cap(4)",
                                               "end_cap(4) + H(1) => blend(3)", "H(1) + H(1) => blend(3)"]
        assert [r.kinetics.A.value_si for r in reactions] == pytest.approx([1.0, 2.0e-6, 3.0e-6], rel=1e-12)
        assert sorted(s.label for s in species if s.thermo is not None) == ["H", "blend", "end_cap"]


class TestKeywordsAreTokens:
    """An auxiliary line's keyword (DUPLICATE, LOW, HIGH, TROE, SRI, CHEB,
    PLOG, REV, COV) is its first token, and a block header is a line's first
    token: a collider label that merely holds one (``blow(4)``, ``dup_x(5)``)
    on an efficiency line, or a thermo entry named ``species(9)``, is read as
    what it is."""

    LABELS = ["blow(4)", "dup_x(5)", "rev(6)", "troex(7)", "highT(8)", "cov1(9)", "cheb1(10)", "plog_x(11)",
              "SRIx(12)"]
    SMILES = ["C", "CC", "CCC", "CCCC", "CCCCC", "CCCCCC", "CCCCCCC", "CCCCCCCC", "CCCCCCCCC"]

    def write(self, tmp_path, names, body, thermo=()):
        text = "\n".join(["ELEMENTS H C N END", "SPECIES", *names, "END", *thermo,
                          "REACTIONS    KCAL/MOLE   MOLES", *body, "END", ""])
        (tmp_path / "chem.inp").write_text(text)
        return str(tmp_path / "chem.inp")

    def dictionary(self, tmp_path, structures):
        (tmp_path / "dict.txt").write_text("\n".join(f"{n}\n{Molecule(smiles=s).to_adjacency_list()}"
                                                     for n, s in structures.items()))
        return str(tmp_path / "dict.txt")

    @pytest.mark.parametrize("label", LABELS)
    def test_a_collider_label_holding_a_keyword_is_an_efficiency(self, tmp_path, label):
        structures = {"H(1)": "[H]", "H2(2)": "[H][H]", "N2": "N#N", label: self.SMILES[self.LABELS.index(label)]}
        path = self.write(tmp_path, structures, ["H(1)+H(1)+M<=>H2(2)+M    1.0 0.0 0.0", f"{label}/2.0/ N2/0.5/"])
        _species, reactions = load_chemkin_file(path, self.dictionary(tmp_path, structures),
                                                check_duplicates=False, generate_resonance_structures=False)
        (reaction,) = reactions
        assert reaction.reversible and not reaction.duplicate
        assert type(reaction.kinetics).__name__ == "ThirdBody"
        efficiencies = {m.to_smiles(): e for m, e in reaction.kinetics.efficiencies.items()}
        assert efficiencies == {Molecule(smiles=structures[label]).to_smiles(): 2.0, "N#N": 0.5}

    def test_a_thermo_entry_named_like_a_header_is_no_block(self, tmp_path):
        """With a separate thermo file the mechanism's own THERM block is not
        read; its entry lines must not open a SPECIES or SITE block."""
        entry = ["{0:<24}C   1H   4          G   200.000  6000.000 1000.00      1",
                 " 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2",
                 "-7.45375000E+02 3.35532000E+00 2.50000000E+00 0.00000000E+00 0.00000000E+00    3",
                 " 0.00000000E+00 0.00000000E+00-7.45375000E+02 3.35532000E+00                   4"]
        structures = {"species(9)": "C", "SITE_a(2)": "CC", "H(1)": "[H]"}
        therm = ["THERM ALL", "   300.000  1000.000  5000.000"]
        for name in structures:
            therm += [line.format(name) for line in entry]
        therm += ["END"]
        path = self.write(tmp_path, structures, ["species(9)+H(1)=>SITE_a(2)    1.0 0.0 0.0"], thermo=therm)
        (tmp_path / "therm.dat").write_text("\n".join(therm) + "\n")
        species, reactions = load_chemkin_file(path, self.dictionary(tmp_path, structures),
                                               thermo_path=str(tmp_path / "therm.dat"), check_duplicates=False,
                                               generate_resonance_structures=False)
        assert sorted(s.label for s in species) == ["H", "SITE_a", "species"]
        assert [str(r) for r in reactions] == ["species(9) + H(1) => SITE_a(2)"]


class TestAuxiliaryRecordsAndHeaders:
    """Chemkin allows several auxiliary records on one line, each with the
    meaning it has on a line of its own, and the abbreviated block headers
    ELEM, SPEC, THER and REAC."""

    THERM_ENTRY = TestEndKeyword.THERM_ENTRY
    NAMES = {"N2": "N#N", "H(1)": "[H]", "H2(2)": "[H][H]", "CH4(3)": "C"}

    def load(self, tmp_path, reactions, headers=("ELEMENTS", "SPECIES", "THERM ALL", "REACTIONS    KCAL/MOLE   MOLES"),
             name="chem"):
        therm = [headers[2], "   300.000  1000.000  5000.000"]
        for species in self.NAMES:
            therm += [line.format(species) for line in self.THERM_ENTRY]
        text = "\n".join([f"{headers[0]} H C N END", headers[1], *self.NAMES, "END", *therm, "END",
                          headers[3], *reactions, "END", ""])
        path = tmp_path / f"{name}.inp"
        path.write_text(text)
        (tmp_path / "dict.txt").write_text("\n".join(f"{n}\n{Molecule(smiles=s).to_adjacency_list()}"
                                                     for n, s in self.NAMES.items()))
        return load_chemkin_file(str(path), str(tmp_path / "dict.txt"), check_duplicates=False,
                                 generate_resonance_structures=False)

    @staticmethod
    def described(kinetics):
        efficiencies = sorted((m.to_smiles(), e) for m, e in getattr(kinetics, "efficiencies", {}).items())
        pressures = list(kinetics.pressures.value_si) if hasattr(kinetics, "pressures") else None
        return type(kinetics).__name__, repr(kinetics), efficiencies, pressures

    CASES = {
        "low_troe": (["H(1)+H(1)(+M)=>H2(2)(+M)    1.0e13 0.0 0.0", "LOW/1e24 -2 0/ TROE/0.5 100 1000/"],
                     ["H(1)+H(1)(+M)=>H2(2)(+M)    1.0e13 0.0 0.0", "LOW/1e24 -2 0/", "TROE/0.5 100 1000/"],
                     "Troe"),
        "plog_plog": (["CH4(3)=>H(1)+H(1)+H2(2)    1.0 0.0 0.0", "PLOG/0.1 1e10 0 0/ PLOG/1.0 1e11 0 0/"],
                      ["CH4(3)=>H(1)+H(1)+H2(2)    1.0 0.0 0.0", "PLOG/0.1 1e10 0 0/", "PLOG/1.0 1e11 0 0/"],
                      "PDepArrhenius"),
        "low_efficiency": (["H(1)+H(1)(+M)=>H2(2)(+M)    1.0e13 0.0 0.0", "LOW/1e24 -2 0/ N2/0.5/"],
                           ["H(1)+H(1)(+M)=>H2(2)(+M)    1.0e13 0.0 0.0", "LOW/1e24 -2 0/", "N2/0.5/"],
                           "Lindemann"),
        "efficiency_low": (["H(1)+H(1)(+M)=>H2(2)(+M)    1.0e13 0.0 0.0", "N2/0.5/ LOW/1e24 -2 0/"],
                           ["H(1)+H(1)(+M)=>H2(2)(+M)    1.0e13 0.0 0.0", "N2/0.5/", "LOW/1e24 -2 0/"],
                           "Lindemann"),
    }

    @pytest.mark.parametrize("case", sorted(CASES))
    def test_records_on_one_line_read_as_on_separate_lines(self, tmp_path, case):
        together, apart, kind = self.CASES[case]
        _s, (one,) = self.load(tmp_path, together, name="together")
        _s, (two,) = self.load(tmp_path, apart, name="apart")
        assert self.described(one.kinetics) == self.described(two.kinetics)
        assert type(one.kinetics).__name__ == kind
        if case == "plog_plog":
            assert len(one.kinetics.pressures.value_si) == 2
        if case.startswith(("low_e", "efficiency")):
            assert self.described(one.kinetics)[2] == [("N#N", 0.5)]

    def test_a_line_with_an_unreadable_record_is_an_error(self, tmp_path):
        with pytest.raises(Exception):
            self.load(tmp_path, ["H(1)+H(1)(+M)=>H2(2)(+M)    1.0e13 0.0 0.0", "LOW/1e24 -2 0/ TROE/0.5/"])

    @pytest.mark.parametrize("which", [0, 1, 2, 3])
    def test_abbreviated_headers_load_the_same_model(self, tmp_path, which):
        full = ["ELEMENTS", "SPECIES", "THERM ALL", "REACTIONS    KCAL/MOLE   MOLES"]
        short = list(full)
        short[which] = ["ELEM", "SPEC", "THER", "REAC    KCAL/MOLE   MOLES"][which]
        body = ["CH4(3)=>H(1)+H(1)+H2(2)    1.0 0.0 3.0", "H(1)+H(1)(+M)=>H2(2)(+M)    1.0e13 0.0 0.0",
                "LOW/1e24 -2 0/", "N2/0.5/"]
        sp1, rx1 = self.load(tmp_path, body, headers=full, name="full")
        sp2, rx2 = self.load(tmp_path, body, headers=short, name="short")
        assert [s.label for s in sp2] == [s.label for s in sp1] and all(s.thermo is not None for s in sp2)
        assert [str(r) for r in rx2] == [str(r) for r in rx1] and len(rx2) == 2
        assert [self.described(r.kinetics) for r in rx2] == [self.described(r.kinetics) for r in rx1]

    def test_reac_keeps_its_units(self, tmp_path):
        body = ["CH4(3)=>H(1)+H(1)+H2(2)    1.0 0.0 1000.0"]
        _s, (full,) = self.load(tmp_path, body, headers=("ELEMENTS", "SPECIES", "THERM ALL",
                                                         "REACTIONS MOLES KELVINS"), name="full")
        _s, (short,) = self.load(tmp_path, body, headers=("ELEMENTS", "SPECIES", "THERM ALL", "REAC MOLES KELVINS"),
                                 name="short")
        assert short.kinetics.Ea.value_si == pytest.approx(1000.0 * constants.R, rel=1e-12)
        assert repr(short.kinetics) == repr(full.kinetics)


# a SITE header whose site density puts a 1 in column 80 (spar r158 item 4)
SITE_DENSITY_IN_COLUMN_80 = "SITE/SURF/   SDEN/2.7" + "0" * (79 - len("SITE/SURF/   SDEN/2.7")) + "1E-09/"
assert SITE_DENSITY_IN_COLUMN_80[79] == "1"


class TestRecordShapes:
    """A record whose name is both a keyword and a declared species is read by
    its shape; an unknown bare word, values with no name, and an unreadable
    record fail the line; comments are not records; an inline THERM block
    skipped for a separate thermo file is skipped whole."""

    THERM_ENTRY = TestEndKeyword.THERM_ENTRY
    BASE = {"N2": "N#N", "H(1)": "[H]", "H2(2)": "[H][H]"}
    EXTRA = {"REV": "CC", "DUP": "CCC", "LOW": "CCCC", "CHEB": "CCCCC"}
    EQUATION = "H(1)+H(1)(+M)<=>H2(2)(+M)    1.0e13 0.0 0.0"

    def load(self, tmp_path, aux, names=None, thermo_header="THERM ALL", thermo_path=False, name="chem",
             equation=None):
        names = dict(self.BASE, **self.EXTRA) if names is None else names
        therm = [thermo_header, "   300.000  1000.000  5000.000"]
        for species in names:
            therm += [line.format(species) for line in self.THERM_ENTRY]
        therm += ["END"]
        text = "\n".join(["ELEMENTS H C N END", "SPECIES", *names, "END", *therm,
                          "REACTIONS    KCAL/MOLE   MOLES", equation or self.EQUATION, *aux, "END", ""])
        path = tmp_path / f"{name}.inp"
        path.write_text(text)
        (tmp_path / "dict.txt").write_text("\n".join(f"{n}\n{Molecule(smiles=s).to_adjacency_list()}"
                                                     for n, s in names.items()))
        kwargs = {}
        if thermo_path:
            (tmp_path / "therm.dat").write_text("\n".join(therm) + "\n")
            kwargs["thermo_path"] = str(tmp_path / "therm.dat")
        return load_chemkin_file(str(path), str(tmp_path / "dict.txt"), check_duplicates=False,
                                 generate_resonance_structures=False, **kwargs)

    @staticmethod
    def efficiencies(kinetics):
        return sorted((m.to_smiles(), e) for m, e in kinetics.efficiencies.items())

    @pytest.mark.parametrize("label, value", [("REV", "0.00"), ("DUP", "2"), ("LOW", "2")])
    def test_a_species_named_like_a_keyword_keeps_its_efficiency(self, tmp_path, label, value):
        _s, (reaction,) = self.load(tmp_path, [f"N2/0.50/ {label}/{value}/", "LOW/1e24 -2 0/"])
        assert reaction.reversible and not reaction.duplicate
        assert type(reaction.kinetics).__name__ == "Lindemann"
        assert self.efficiencies(reaction.kinetics) == sorted(
            [("N#N", 0.5), (Molecule(smiles=self.EXTRA[label]).to_smiles(), float(value))])
        assert reaction.kinetics.arrheniusLow.A.value_si == pytest.approx(1e24 * 1e-12, rel=1e-12)

    def test_a_keyword_record_of_keyword_shape_is_the_keyword(self, tmp_path):
        _s, (reaction,) = self.load(tmp_path, ["N2/0.50/ REV/1e13 0 0/", "LOW/1e24 -2 0/"])
        assert self.efficiencies(reaction.kinetics) == [("N#N", 0.5)]
        assert "explicit reverse" in reaction.kinetics.comment

    def test_a_single_cheb_value_on_a_falloff_reaction_is_an_efficiency(self, tmp_path):
        """Spar r153 P2: one value is a Chebyshev coefficient only on a
        reaction already read as Chebyshev; on a falloff reaction it is the
        declared species' efficiency."""
        _s, (reaction,) = self.load(tmp_path, ["N2/0.50/ CHEB/2/", "LOW/1e24 -2 0/"])
        assert type(reaction.kinetics).__name__ == "Lindemann"
        assert self.efficiencies(reaction.kinetics) == sorted([("N#N", 0.5), ("CCCCC", 2.0)])

    def test_a_record_that_fits_no_reading_names_both(self, tmp_path):
        """LOW on a third-body reaction: no falloff collider for the keyword,
        three values for an efficiency of the declared species LOW."""
        with pytest.raises(ChemkinError, match="keyword LOW.*efficiency"):
            self.load(tmp_path, ["N2/0.50/ LOW/1e24 -2 0/"], equation="H(1)+H(1)+M<=>H2(2)+M    1.0e13 0.0 0.0")

    def test_a_lowercase_record_names_the_capital_species(self, tmp_path):
        """Species names are matched as written, else upper-cased (the old
        reader's fallback): ``cheb/2/`` is the efficiency of the declared
        CHEB; the keyword test uses the same lookup."""
        _s, (reaction,) = self.load(tmp_path, ["N2/0.50/ cheb/2/", "LOW/1e24 -2 0/"])
        assert self.efficiencies(reaction.kinetics) == sorted([("N#N", 0.5), ("CCCCC", 2.0)])

    def test_a_capital_record_is_no_species_declared_lowercase(self, tmp_path):
        """Only ``cheb`` is declared: ``CHEB/2/`` names no species (never
        lower-cased), so it is the keyword, which a falloff reaction cannot
        take."""
        names = dict(self.BASE, cheb="CCCCC")
        with pytest.raises(ChemkinError, match="CHEB"):
            self.load(tmp_path, ["N2/0.50/ CHEB/2/", "LOW/1e24 -2 0/"], names=names)

    def test_rmgs_writer_output_round_trips(self, tmp_path):
        polys = [NASAPolynomial(coeffs=[2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 3.35], Tmin=lo, Tmax=hi)
                 for lo, hi in (((200, "K"), (1000, "K")), ((1000, "K"), (6000, "K")))]
        species = []
        for label, smiles in [("N2", "N#N"), ("H", "[H]"), ("H2", "[H][H]"), ("REV", "CC"), ("DUP", "CCC"),
                              ("LOW", "CCCC")]:
            spc = Species(label=label).from_smiles(smiles)
            spc.label = label
            spc.thermo = NASA(polynomials=polys, Tmin=(200, "K"), Tmax=(6000, "K"))
            species.append(spc)
        n2, h, h2, rev, dup, low = species
        kinetics = Lindemann(arrheniusHigh=Arrhenius(A=(1e7, "m^3/(mol*s)"), n=0, Ea=(0, "J/mol"), T0=(1, "K")),
                             arrheniusLow=Arrhenius(A=(1e12, "m^6/(mol^2*s)"), n=-2, Ea=(0, "J/mol"), T0=(1, "K")),
                             efficiencies={n2.molecule[0]: 0.5, rev.molecule[0]: 0.0, dup.molecule[0]: 2.0,
                                           low.molecule[0]: 2.0})
        reaction = Reaction(reactants=[h, h], products=[h2], kinetics=kinetics, reversible=True)
        save_chemkin_file(str(tmp_path / "chem.inp"), species, [reaction], verbose=False, check_for_duplicates=False)
        save_species_dictionary(str(tmp_path / "dict.txt"), species)
        _s, (back,) = load_chemkin_file(str(tmp_path / "chem.inp"), str(tmp_path / "dict.txt"),
                                        check_duplicates=False, generate_resonance_structures=False)
        assert back.reversible and not back.duplicate
        assert self.efficiencies(back.kinetics) == self.efficiencies(kinetics)

    @pytest.mark.parametrize("thermo_path", [False, True])
    @pytest.mark.parametrize("block", [["REACTIONS    KCAL/MOLE   MOLES"], ["REAC"], ["SPECIES", "N2", "END"],
                                       ["ELEMENTS H C N", "END"], ["THERMO", "END"], ["SITE/SURF/", "END"],
                                       [SITE_DENSITY_IN_COLUMN_80, "END"]],
                             ids=["REACTIONS", "REAC", "SPECIES", "ELEMENTS", "THERMO", "SITE", "SITE_column_80"])
    def test_an_unterminated_thermo_block_fails_by_name(self, tmp_path, thermo_path, block):
        """Spar r153 P1: a THERM block without its END, read or skipped for a
        separate thermo file, stops at the next block's header with an
        error naming both; it never consumes another block."""
        therm = ["THERM ALL", "   300.000  1000.000  5000.000"]
        for species in self.BASE:
            therm += [line.format(species) for line in self.THERM_ENTRY]
        reactions = ["REACTIONS    KCAL/MOLE   MOLES", self.EQUATION, "LOW/1e24 -2 0/", "END"]
        text = "\n".join(["ELEMENTS H C N END", "SPECIES", *self.BASE, "END", *therm,
                          *(block if block[0] not in ("REACTIONS    KCAL/MOLE   MOLES", "REAC") else []),
                          *(["REAC"] + reactions[1:] if block == ["REAC"] else reactions), ""])
        (tmp_path / "chem.inp").write_text(text)
        (tmp_path / "dict.txt").write_text("\n".join(f"{n}\n{Molecule(smiles=s).to_adjacency_list()}"
                                                     for n, s in self.BASE.items()))
        kwargs = {}
        if thermo_path:
            (tmp_path / "therm.dat").write_text("\n".join(therm + ["END"]) + "\n")
            kwargs["thermo_path"] = str(tmp_path / "therm.dat")
        with pytest.raises(ChemkinError, match=r"THERM.*END.*" + block[0].split()[0].split("/")[0]):
            load_chemkin_file(str(tmp_path / "chem.inp"), str(tmp_path / "dict.txt"), check_duplicates=False,
                              generate_resonance_structures=False, **kwargs)

    @pytest.mark.parametrize("thermo_path", [False, True])
    def test_an_inline_thermo_entry_named_reac_loads(self, tmp_path, thermo_path):
        names = dict(self.BASE, REAC="CC")
        species, (reaction,) = self.load(tmp_path, ["LOW/1e24 -2 0/"], names=names, thermo_path=thermo_path)
        assert sorted(s.label for s in species) == ["H", "H2", "N2", "REAC"]
        assert all(s.thermo is not None for s in species)
        assert type(reaction.kinetics).__name__ == "Lindemann"

    def test_a_separate_thermo_file_headed_ther_is_read(self, tmp_path):
        species, _r = self.load(tmp_path, ["LOW/1e24 -2 0/"], names=dict(self.BASE), thermo_header="THER",
                                thermo_path=True)
        assert all(s.thermo is not None for s in species)

    @pytest.mark.parametrize("aux", [["TYPO N2/0.5/", "LOW/1e24 -2 0/"], ["N2/0.5/ GARBAGE", "LOW/1e24 -2 0/"],
                                     ["LOW/1e24 -2 0/ /0.5/"]], ids=["typo", "garbage", "values_with_no_name"])
    def test_an_unreadable_record_fails_the_line(self, tmp_path, aux):
        with pytest.raises(ChemkinError):
            self.load(tmp_path, aux, names=dict(self.BASE))

    def test_a_bare_keyword_before_a_record_is_read(self, tmp_path):
        _s, (together,) = self.load(tmp_path, ["DUP LOW/1e24 -2 0/"], names=dict(self.BASE), name="together")
        _s, (apart,) = self.load(tmp_path, ["DUP", "LOW/1e24 -2 0/"], names=dict(self.BASE), name="apart")
        assert together.duplicate and apart.duplicate
        assert repr(together.kinetics) == repr(apart.kinetics)

    def test_a_comment_holding_a_slash_is_no_record(self):
        species = {label: Species(label=label).from_smiles(smiles) for label, smiles in self.BASE.items()}
        aunits = ['', 's^-1', 'cm^3/(mol*s)', 'cm^6/(mol^2*s)', 'cm^9/(mol^3*s)']
        aunits_surf = ['', 's^-1', 'cm^2/(mol*s)', 'cm^4/(mol^2*s)']
        entry = "H(1)+H(1)(+M)=>H2(2)(+M)    1.0e13 0.0 0.0\nLOW/1e24 -2 0/ ! comment / slash"
        reaction = read_kinetics_entry(entry, species, aunits, aunits_surf, "kcal/mol")
        assert type(reaction.kinetics).__name__ == "Lindemann"
        assert reaction.kinetics.arrheniusLow.A.value_si == pytest.approx(1e24 * 1e-12, rel=1e-12)

    # every keyword the reader knows, as a species name, and a lowercase cheb beside CHEB
    KEYWORD_SPECIES = ["CHEB", "cheb", "REV", "DUP", "DUPL", "DUPLICATE", "LOW", "HIGH", "TROE", "SRI", "PLOG",
                       "TCHEB", "PCHEB", "COV", "COVERAGE", "STICK", "MOME", "XSMI"]

    @staticmethod
    def _close(a, b):
        return a == pytest.approx(b, rel=1e-9, abs=1e-12)

    def _same_arrhenius(self, a, b):
        return all(self._close(x.value_si, y.value_si) for x, y in ((a.A, b.A), (a.Ea, b.Ea), (a.T0, b.T0))) \
            and self._close(a.n.value_si, b.n.value_si)

    def _same_kinetics(self, a, b):
        assert type(a) is type(b), (type(a), type(b))
        if isinstance(a, Chebyshev):
            return (a.coeffs.value_si.shape == b.coeffs.value_si.shape
                    and all(self._close(x, y) for x, y in zip(a.coeffs.value_si.flat, b.coeffs.value_si.flat))
                    and all(self._close(getattr(a, k).value_si, getattr(b, k).value_si)
                            for k in ("Tmin", "Tmax", "Pmin", "Pmax")))
        if isinstance(a, PDepArrhenius):
            return (all(self._close(x, y) for x, y in zip(a.pressures.value_si, b.pressures.value_si))
                    and len(a.arrhenius) == len(b.arrhenius)
                    and all(self._same_kinetics(x, y) for x, y in zip(a.arrhenius, b.arrhenius)))
        if isinstance(a, MultiArrhenius):
            return len(a.arrhenius) == len(b.arrhenius) and all(
                self._same_arrhenius(x, y) for x, y in zip(a.arrhenius, b.arrhenius))
        if isinstance(a, Arrhenius):
            return self._same_arrhenius(a, b)
        same = self.efficiencies(a) == self.efficiencies(b) and self._same_arrhenius(a.arrheniusLow, b.arrheniusLow)
        if isinstance(a, (Lindemann, Troe)):
            same = same and self._same_arrhenius(a.arrheniusHigh, b.arrheniusHigh)
        if isinstance(a, Troe):
            same = same and self._close(a.alpha, b.alpha) and all(
                (getattr(a, k) is None and getattr(b, k) is None) or self._close(getattr(a, k).value_si,
                                                                                  getattr(b, k).value_si)
                for k in ("T1", "T2", "T3"))
        return same

    def test_every_entry_rmgs_writer_emits_reads_back(self, tmp_path):
        """Spar r153 oracle: a model whose species are named like every
        keyword (and cheb beside CHEB), with Chebyshev (one, three and six
        pressure terms), Lindemann, Troe, third-body, PLOG and duplicate
        reactions, written by RMG's own writer, reads back exactly."""
        polys = [NASAPolynomial(coeffs=[2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 3.35], Tmin=lo, Tmax=hi)
                 for lo, hi in (((200, "K"), (1000, "K")), ((1000, "K"), (6000, "K")))]
        species = {}
        for label, smiles in [("N2", "N#N"), ("H", "[H]"), ("H2", "[H][H]")] + [
                (x, "C" * (k + 1)) for k, x in enumerate(self.KEYWORD_SPECIES)]:
            spc = Species(label=label).from_smiles(smiles)
            spc.label = label
            spc.thermo = NASA(polynomials=polys, Tmin=(200, "K"), Tmax=(6000, "K"))
            species[label] = spc
        h, h2 = species["H"], species["H2"]
        eff = {species[x].molecule[0]: 0.25 * (k + 1) for k, x in enumerate(["N2"] + self.KEYWORD_SPECIES)}

        def arr(a, units, n=0.5, ea=2.0):
            return Arrhenius(A=(a, units), n=n, Ea=(ea, "kcal/mol"), T0=(1, "K"))

        def cheb(rows, cols, units):
            coeffs = [[1.25 - 0.5 * i - 0.125 * j for j in range(cols)] for i in range(rows)]
            return Chebyshev(coeffs=coeffs, kunits=units, Tmin=(300, "K"), Tmax=(2000, "K"),
                             Pmin=(0.01, "atm"), Pmax=(100, "atm"))
        bi, ter = "cm^3/(mol*s)", "cm^6/(mol^2*s)"
        kinetics = [
            Lindemann(arrheniusHigh=arr(1.5e13, bi), arrheniusLow=arr(2.5e18, ter, n=-1.5), efficiencies=eff),
            Troe(arrheniusHigh=arr(1.5e13, bi), arrheniusLow=arr(2.5e18, ter, n=-1.5), alpha=0.5, T3=(100, "K"),
                 T1=(1000, "K"), efficiencies=eff),
            Troe(arrheniusHigh=arr(1.5e13, bi), arrheniusLow=arr(2.5e18, ter, n=-1.5), alpha=0.5, T3=(100, "K"),
                 T1=(1000, "K"), T2=(5000, "K"), efficiencies=eff),
            ThirdBody(arrheniusLow=arr(2.5e18, ter, n=-1.5), efficiencies=eff),
            cheb(3, 1, bi), cheb(2, 3, bi), cheb(5, 6, bi),
            PDepArrhenius(pressures=([0.1, 1.0, 10.0], "atm"), arrhenius=[arr(1.5e13, bi), arr(2.5e13, bi),
                                                                          arr(3.5e13, bi)]),
            arr(1.5e13, bi), arr(2.5e13, bi, n=1.5),
        ]
        reactions = [Reaction(reactants=[h, h], products=[h2], kinetics=k, reversible=True,
                              duplicate=j >= len(kinetics) - 2) for j, k in enumerate(kinetics)]
        save_chemkin_file(str(tmp_path / "chem.inp"), list(species.values()), reactions, verbose=False,
                          check_for_duplicates=False)
        save_species_dictionary(str(tmp_path / "dict.txt"), list(species.values()))
        _s, back = load_chemkin_file(str(tmp_path / "chem.inp"), str(tmp_path / "dict.txt"),
                                     check_duplicates=False, generate_resonance_structures=False)
        assert len(back) == len(reactions)
        for original, read in zip(reactions, back):
            assert [s.label for s in read.reactants] == ["H", "H"] and [s.label for s in read.products] == ["H2"]
            assert read.reversible and read.duplicate == original.duplicate
            assert self._same_kinetics(original.kinetics, read.kinetics), (original.kinetics, read.kinetics)

    # ---- I-177 rework 3 (spar r158) ----
    SURFACE_UNITS = (['', 's^-1', 'cm^3/(mol*s)', 'cm^6/(mol^2*s)', 'cm^9/(mol^3*s)'],
                     ['', 's^-1', 'cm^2/(mol*s)', 'cm^4/(mol^2*s)', 'cm^6/(mol^3*s)'])

    @pytest.mark.parametrize("cov_first", [False, True])
    def test_stick_keeps_the_coverage_of_an_earlier_cov(self, cov_first):
        """Spar r158 item 1: STICK turns the surface Arrhenius into a sticking
        coefficient carrying every field earlier records set: COV before
        STICK keeps its coverage dependence, as after it."""
        x = Species().from_adjacency_list("1 X u0 p0")
        hx = Species().from_adjacency_list("1 H u0 p0 {2,S} \n 2 X u0 p0 {1,S}")
        h2 = Species().from_smiles("[H][H]")
        species = {"X": x, "H*": hx, "H2": h2}
        lines = ["    STICK", "    COV / H*    0.5    1.0    1.099 /"]
        entry = "\n".join(["H2+X+X=H*+H*    4.000e-02 0.000 7.098"] + (lines[::-1] if cov_first else lines))
        reaction = read_kinetics_entry(entry, species, *self.SURFACE_UNITS, "kcal/mol")
        assert isinstance(reaction.kinetics, StickingCoefficient)
        (spc, cov), = reaction.kinetics.coverage_dependence.items()
        # I-179 item 3: a, m and E are quantities, as the coverage_dependence setter makes them, in both orders
        assert spc is hx and cov["a"].value_si == 0.5 and cov["m"].value_si == 1.0
        assert cov["E"].value_si == pytest.approx(1.099 * 4184.0, rel=1e-12)
        assert reaction.kinetics.A.value_si == pytest.approx(4.0e-2, rel=1e-12)
        assert reaction.kinetics.Ea.value_si == pytest.approx(7.098 * 4184.0, rel=1e-12)
        # the read SurfaceReactor makes of a coverage dependence (rmgpy/solver/surface.pyx initialize_model)
        assert [(p['a'].value_si, p['m'].value_si, p['E'].value_si)
                for p in reaction.kinetics.coverage_dependence.values()] == [(0.5, 1.0, cov["E"].value_si)]

    @pytest.mark.parametrize("cov_first", [False, True])
    def test_a_read_coverage_dependence_writes_and_reads_back(self, tmp_path, cov_first):
        """I-179 item 3: the quantities RMG's writer gets from the reader
        write and read back unchanged, in both orders."""
        x = Species().from_adjacency_list("1 X u0 p0")
        x.label = "X"
        hx = Species().from_adjacency_list("1 H u0 p0 {2,S} \n 2 X u0 p0 {1,S}")
        hx.label = "H*"
        h2 = Species().from_smiles("[H][H]")
        h2.label = "H2"
        species = {"X": x, "H*": hx, "H2": h2}
        lines = ["    STICK", "    COV / H*    0.5    1.0    1.099 /"]
        entry = "\n".join(["H2+X+X=H*+H*    4.000e-02 0.000 7.098"] + (lines[::-1] if cov_first else lines))
        reaction = read_kinetics_entry(entry, species, *self.SURFACE_UNITS, "kcal/mol")
        text = write_kinetics_entry(reaction, list(species.values()), verbose=False)
        back = read_kinetics_entry(text, species, *self.SURFACE_UNITS, "kcal/mol")
        (_s, cov), = back.kinetics.coverage_dependence.items()
        assert (cov["a"].value_si, cov["m"].value_si) == (0.5, 1.0)
        assert cov["E"].value_si == pytest.approx(1.099 * 4184.0, rel=1e-6)

    def test_duplicate_chebyshev_fits_write_and_read_back(self, tmp_path):
        """I-179 item 1: duplicate Chebyshev pairs of 2x6 (twelve coefficients,
        the last CHEB line unfinished), 3x1 and 5x6 write DUPLICATE on its own
        line and read back equal."""
        polys = [NASAPolynomial(coeffs=[2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 3.35], Tmin=lo, Tmax=hi)
                 for lo, hi in (((200, "K"), (1000, "K")), ((1000, "K"), (6000, "K")))]
        species = []
        for label, smiles in [("N2", "N#N"), ("H", "[H]"), ("H2", "[H][H]")]:
            spc = Species(label=label).from_smiles(smiles)
            spc.label = label
            spc.thermo = NASA(polynomials=polys, Tmin=(200, "K"), Tmax=(6000, "K"))
            species.append(spc)
        _n2, h, h2 = species

        def cheb(rows, cols, shift):
            coeffs = [[1.25 - 0.5 * i - 0.125 * j + shift for j in range(cols)] for i in range(rows)]
            return Chebyshev(coeffs=coeffs, kunits="cm^3/(mol*s)", Tmin=(300, "K"), Tmax=(2000, "K"),
                             Pmin=(0.01, "atm"), Pmax=(100, "atm"))
        reactions = [Reaction(reactants=[h, h], products=[h2], kinetics=cheb(rows, cols, shift), reversible=True,
                              duplicate=True)
                     for rows, cols in ((2, 6), (3, 1), (5, 6)) for shift in (0.0, 0.25)]
        save_chemkin_file(str(tmp_path / "chem.inp"), species, reactions, verbose=False, check_for_duplicates=False)
        save_species_dictionary(str(tmp_path / "dict.txt"), species)
        text = (tmp_path / "chem.inp").read_text()
        assert all(line.strip() == "DUPLICATE" for line in text.splitlines() if "DUPLICATE" in line)
        _s, back = load_chemkin_file(str(tmp_path / "chem.inp"), str(tmp_path / "dict.txt"),
                                     check_duplicates=False, generate_resonance_structures=False)
        assert len(back) == len(reactions)
        for original, read in zip(reactions, back):
            assert read.duplicate and self._same_kinetics(original.kinetics, read.kinetics)

    def test_a_therm_block_ending_in_a_comment_fails_by_name(self, tmp_path):
        """I-179 item 2: THERM ALL followed only by a comment line: end of
        file stops the reader with an error naming the block, never a hang
        (guarded: the read runs in a daemon thread given 60 s)."""
        import threading
        (tmp_path / "chem.inp").write_text("ELEMENTS H END\nSPECIES\nH2\nEND\nTHERM ALL\n! only a comment\n")
        (tmp_path / "dict.txt").write_text("H2\n" + Molecule(smiles="[H][H]").to_adjacency_list())
        outcome = []

        def read():
            try:
                load_chemkin_file(str(tmp_path / "chem.inp"), str(tmp_path / "dict.txt"),
                                  generate_resonance_structures=False)
                outcome.append(None)
            except Exception as error:
                outcome.append(error)
        worker = threading.Thread(target=read, daemon=True)
        worker.start()
        worker.join(60.0)
        assert not worker.is_alive(), "the reader hangs on a THERM block ending at end of file"
        assert isinstance(outcome[0], ChemkinError) and "THERM" in str(outcome[0])

    @pytest.mark.parametrize("order", [0, 1])
    def test_plog_with_an_optional_collider_is_plog_beside_a_species_plog(self, tmp_path, order):
        """Spar r158 item 2: PLOG on a (+M) reaction (Cantera strips the
        optional collider) stays PLOG with a species named PLOG declared."""
        names = dict(self.BASE, PLOG="CC")
        aux = ["PLOG/1.0 1e13 0 0/", "PLOG/10.0 2e13 0 0/"]
        _s, (reaction,) = self.load(tmp_path, aux if order == 0 else aux[::-1], names=names)
        assert type(reaction.kinetics).__name__ == "PDepArrhenius"
        assert [p for p in reaction.kinetics.pressures.value_si] == pytest.approx([101325.0, 1013250.0])

    CHEB_LINES = ["TCHEB/ 300.0 2000.0 /", "PCHEB/ 0.01 100.0 /", "CHEB/ 2 1 /", "CHEB/ 1.25 /", "CHEB/ 0.5 /"]

    @pytest.mark.parametrize("position", range(6))
    def test_low_beside_chebyshev_reads_in_every_order(self, tmp_path, position):
        """Spar r158 item 3: LOW alongside Chebyshev records (the manual
        allows it) reads in every line order with a species named LOW
        declared: the reaction is Chebyshev, the LOW record the keyword."""
        names = dict(self.BASE, LOW="CC")
        aux = self.CHEB_LINES[:position] + ["LOW/1e24 -2 0/"] + self.CHEB_LINES[position:]
        _s, (reaction,) = self.load(tmp_path, aux, names=names)
        assert type(reaction.kinetics).__name__ == "Chebyshev"
        assert reaction.kinetics.coeffs.value_si.shape == (2, 1)

