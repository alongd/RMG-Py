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

"""Regression contract for the second independent state-identity review."""
import os
import pickle
import runpy
from pathlib import Path
from multiprocessing.reduction import ForkingPickler

import pytest
import numpy as np

from rmgpy.data.kinetics.family import complete_round_trip
from rmgpy.data.thermo import ThermoDatabase, bicyclic_decomposition_for_polyring
from rmgpy.thermo import ThermoData
from rmgpy.exceptions import SpeciesError, StateProvenanceError, SpeciesIdentityError
from rmgpy.kinetics.arrhenius import get_w0
from rmgpy.molecule import Molecule
from rmgpy.molecule.draw import MoleculeDrawer
from rmgpy.reaction import Reaction
from rmgpy.species import Species
assert_isolated_reader = runpy.run_path(str(Path(__file__).with_name('stateReworkTest.py')))['assert_isolated_reader']


class TestR89Finding1:
    @pytest.mark.parametrize('transport', ['copy', 'forking_pickle'])
    @pytest.mark.parametrize('mutation', ['electronic', 'vibrational', 'unresolved'])
    @pytest.mark.parametrize('first_read', ['augmented', 'fingerprint'])
    def test_transported_cache_tracks_mutation_before_any_first_read(self, transport, mutation, first_read):
        source = Species(molecule=[Molecule(smiles='N#N', electronic_state='A', vibrational_level=1)])
        source.get_augmented_inchi()
        source.fingerprint
        copied = complete_round_trip(source) if transport == 'copy' else pickle.loads(ForkingPickler.dumps(source))
        molecule = copied.molecule[0]
        if mutation == 'electronic':
            molecule.electronic_state = 'B'
        elif mutation == 'vibrational':
            molecule.vibrational_level = 2
        else:
            molecule.electronic_state, molecule.vibrational_level = '', -1
        # Do not prime any copied cache before changing its state.
        if first_read == 'augmented':
            assert copied.aug_inchi is None
        else:
            assert copied.fingerprint == molecule.fingerprint
        assert copied.get_augmented_inchi() == molecule.to_augmented_inchi()
        assert copied.fingerprint == molecule.fingerprint
        assert hash(copied) == hash(Species(molecule=[molecule.copy(deep=True)]))


class TestR89Finding2:
    @pytest.mark.parametrize('smiles', ['N#N', 'C', 'CC', 'O', 'O=C=O', 'C=C', 'c1ccccc1'])
    @pytest.mark.parametrize('state', [('', 1), ('A', -1), ('A', 2)])
    @pytest.mark.parametrize('species', [False, True])
    def test_resolved_symmetry_matches_unresolved_structure(self, smiles, state, species):
        ground = Molecule(smiles=smiles)
        excited = Molecule(smiles=smiles, electronic_state=state[0], vibrational_level=state[1])
        before = excited.to_adjacency_list()
        if species:
            expected = Species(molecule=[ground]).get_symmetry_number()
            got = Species(molecule=[excited]).get_symmetry_number()
        else:
            expected = ground.get_symmetry_number()
            got = excited.get_symmetry_number()
        assert got == expected
        # Resonance generation and isomorphism may sort atoms. Preserve identity,
        # rather than requiring a bookkeeping order that symmetry does not promise.
        assert excited.is_isomorphic(Molecule().from_adjacency_list(before))
        assert (excited.electronic_state, excited.vibrational_level) == state

    def test_ring_perception_can_split_a_resolved_temporary_graph(self):
        ground = Molecule(smiles='c1ccccc1Cc1ccccc1')
        excited = Molecule(smiles='c1ccccc1Cc1ccccc1', vibrational_level=1)
        before = excited.to_adjacency_list()
        assert len(ground.get_all_cycles_of_size(6)) == 2
        assert len(excited.get_all_cycles_of_size(6)) == 2
        assert excited.to_adjacency_list() == before

    def test_resonance_generation_preserves_resolved_state_with_separate_ring_cores(self):
        excited = Molecule(smiles='c1ccccc1Cc1ccccc1', electronic_state='A', vibrational_level=1)
        forms = excited.generate_resonance_structures()
        assert forms
        assert all((mol.electronic_state, mol.vibrational_level) == ('A', 1) for mol in forms)
        assert all(len(mol.get_all_cycles_of_size(6)) == 2 for mol in forms)

    @pytest.mark.parametrize('smiles', ['N#N', '[O][O]', 'c1ccccc1Cc1ccccc1'])
    def test_drawing_refuses_resolved_identity(self, tmp_path, smiles):
        excited = Molecule(smiles=smiles, vibrational_level=1)
        before = excited.to_adjacency_list()
        ground = MoleculeDrawer().draw(Molecule(smiles=smiles), 'svg', str(tmp_path / 'ground.svg'))
        assert ground[0] is not None
        with pytest.raises(SpeciesIdentityError, match='MoleculeDrawer.draw'):
            MoleculeDrawer().draw(excited, 'svg', str(tmp_path / 'excited.svg'))
        assert not (tmp_path / 'excited.svg').exists()
        assert excited.to_adjacency_list() == before

    def test_thermo_ring_decomposition_uses_atoms_without_changing_resolved_identity(self):
        ground = Molecule(smiles='C1CC2CCC1C2')
        excited = Molecule(smiles='C1CC2CCC1C2', electronic_state='A', vibrational_level=1)
        before = excited.to_adjacency_list()
        expected, expected_occurrences = bicyclic_decomposition_for_polyring(ground.atoms)
        result, occurrences = bicyclic_decomposition_for_polyring(excited.atoms)
        assert len(result) == len(expected) == 1
        assert sorted(occurrences.values()) == sorted(expected_occurrences.values())
        assert excited.to_adjacency_list() == before

    @pytest.mark.parametrize('smiles,expected_energy', [('[CH2]=*', 1000.0), ('CC.*', 2000.0)])
    def test_thermo_surface_connectivity_keeps_identity_and_refuses_resolved_values(self, smiles, expected_energy):
        # Replace the external trained estimator with known data; exercise the real
        # decomposition/routing API without loading a machine-specific SIDT model.
        class Estimator:
            nodes = {'Root': None}

            def __init__(self, energy):
                self.energy = energy

            def evaluate(self, molecule, **kwargs):
                return np.array([self.energy, 0.0] + [0.0] * 7), np.zeros(9), 'test estimator'

        database = ThermoDatabase()
        database.adsorption_groups = 'SIDT'
        names = ['Pt111_monodentate_adsorption_corrections', 'Pt111_vdw_adsorption_corrections']
        database.sidts = {name: Estimator(index + 1) for index, name in enumerate(names)}
        database.sidt_taggings_and_decompositions = {name: lambda molecule: None for name in names}
        molecule = Molecule(smiles=smiles, vibrational_level=1)
        before = molecule.to_adjacency_list()
        thermo = ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
                            Cpdata=([0.0] * 7, 'J/(mol*K)'), H298=(0.0, 'J/mol'), S298=(0.0, 'J/(mol*K)'))
        with pytest.raises(StateProvenanceError, match='Adsorption'):
            database._add_adsorption_correction(thermo, None, molecule, molecule.get_surface_sites())
        assert thermo.H298.value_si == 0.0
        assert molecule.to_adjacency_list() == before
        # The unchanged unresolved routing still distinguishes connected and vdw adsorbates.
        ground = Molecule(smiles=smiles)
        assert database._add_adsorption_correction(thermo, None, ground, ground.get_surface_sites())
        assert thermo.H298.value_si == expected_energy

    @pytest.mark.parametrize('actions', [[], [['CHANGE_BOND', '*1', -1, '*2']]])
    def test_bond_energy_structural_calculation_can_merge_resolved_reactants(self, actions):
        def reaction(resolved):
            mol = Molecule().from_adjacency_list('1 *1 N u0 p1 c0 {2,T}\n2 *2 N u0 p1 c0 {1,T}')
            if resolved:
                mol.electronic_state, mol.vibrational_level = 'A', 1
            return Reaction(reactants=[Species(molecule=[mol]), Species(smiles='[He]')])
        ground, excited = reaction(False), reaction(True)
        before = excited.reactants[0].to_adjacency_list()
        assert get_w0(actions, excited) == get_w0(actions, ground)
        assert excited.reactants[0].to_adjacency_list() == before


    def test_bond_energy_structural_calculation_can_merge_resolved_products(self):
        def reaction(resolved):
            reactant = next(mol for mol in Molecule(smiles='c1ccccc1').generate_resonance_structures()
                            if any(bond.is_benzene() for atom in mol.atoms for bond in atom.edges.values()))
            atom1, atom2 = next((atom, neighbor) for atom in reactant.atoms
                                for neighbor, bond in atom.edges.items() if bond.is_benzene())
            atom1.label, atom2.label = '*1', '*2'
            product = reactant.copy(deep=True)
            if resolved:
                product.electronic_state, product.vibrational_level = 'A', 1
            return Reaction(reactants=[Species(molecule=[reactant])],
                            products=[Species(molecule=[product]), Species(smiles='[He]')])
        ground, excited = reaction(False), reaction(True)
        before = excited.products[0].to_adjacency_list()
        # The aromatic 1.5 -> 0.5 recipe uses actual product bonds, hence its
        # separate product-merge path rather than the ordinary reactant merge.
        actions = [['CHANGE_BOND', '*1', -1, '*2']]
        assert get_w0(actions, excited) == get_w0(actions, ground) == 0.0
        assert excited.products[0].to_adjacency_list() == before


class TestR89Finding3:
    @staticmethod
    def mixed_species():
        species = Species(molecule=[Molecule(smiles='[CH2]C=C', electronic_state='A', vibrational_level=1)])
        species.generate_resonance_structures()
        assert len(species.molecule) > 1
        species.molecule[0].vibrational_level = 2
        return species

    @pytest.mark.parametrize('access', ['fingerprint', 'aug_inchi', 'inchi', 'smiles', 'multiplicity',
                                      'hash', 'get_augmented_inchi', 'to_adjacency_list'])
    def test_identity_access_refuses_a_mixed_resonance_state(self, access):
        species = self.mixed_species()
        with pytest.raises(SpeciesError, match='state.*match'):
            if access == 'hash':
                hash(species)
            elif access in ('get_augmented_inchi', 'to_adjacency_list'):
                getattr(species, access)()
            else:
                getattr(species, access)

    @pytest.mark.parametrize('level', [1, 2])
    @pytest.mark.parametrize('reverse', [False, True])
    @pytest.mark.parametrize('strict', [False, True])
    @pytest.mark.parametrize('comparison', ['species_isomorphic', 'species_identical', 'reaction_isomorphic'])
    def test_species_and_reaction_comparisons_refuse_both_mixed_states(self, level, reverse, strict, comparison):
        mixed = self.mixed_species()
        candidate = Species(molecule=[Molecule(smiles='[CH2]C=C', electronic_state='A', vibrational_level=level)])
        candidate.generate_resonance_structures()
        left, right = (candidate, mixed) if reverse else (mixed, candidate)
        with pytest.raises(SpeciesError, match='state.*match'):
            if comparison == 'species_isomorphic':
                left.is_isomorphic(right, strict=strict)
            elif comparison == 'species_identical':
                left.is_identical(right, strict=strict)
            else:
                product = Species(smiles='[He]')
                Reaction(reactants=[left], products=[product]).is_isomorphic(Reaction(reactants=[right], products=[product]), strict=strict)


class TestR89Finding5:
    @pytest.mark.parametrize('mode', ['unset', 'multiple', 'alternate_database'])
    def test_isolated_reader_checks_the_checkout_without_environment_pins(self, monkeypatch, tmp_path, mode):
        tree = Path(__file__).resolve().parents[3]
        if mode == 'unset':
            monkeypatch.delenv('PYTHONPATH', raising=False)
        elif mode == 'multiple':
            monkeypatch.setenv('PYTHONPATH', str(tree) + os.pathsep + str(tmp_path))
        else:
            database = tmp_path / 'other_input'
            database.mkdir()
            (tmp_path / 'rmgrc').write_text('database.directory = ' + str(database) + '\n')
            monkeypatch.chdir(tmp_path)
            monkeypatch.setenv('PYTHONPATH', str(tree))
        assert_isolated_reader('assert Molecule(smiles="N#N").vibrational_level == -1', {})
