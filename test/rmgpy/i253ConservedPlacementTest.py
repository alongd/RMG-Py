#!/usr/bin/env python3

"""Conserved owner electron placement survives mechanism export."""

import pytest

from rmgpy.chemkin import save_chemkin_file, save_species_dictionary, write_reaction_string
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.electron_balance import expand_electrons
from rmgpy.electron_placement import FAMILY_ELECTRON_PLACEMENT
from rmgpy.exceptions import MechanismWriterError
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial


OWNER = 'Test_Conserved_Electron_Placement'


def _species():
    electron = Species(label='e-', molecule=[
        Molecule().from_adjacency_list('1 e u1 p0 c-1')])
    argon = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    # Keep the export fixture independent of the excited-state thermo gate;
    # the reaction's Ar* label is enough to exercise the equation round trip.
    excited_argon = Species(label='Ar*').from_adjacency_list('1 Ar u0 p4 c0')
    electron.index, argon.index, excited_argon.index = 1, 2, 3
    return electron, argon, excited_argon


def _declared_reaction(argon, excited_argon):
    reaction = LibraryReaction(
        reactants=[argon], products=[excited_argon], library=OWNER,
        reversible=False, electrons=0,
        kinetics=Arrhenius(A=(1.0e-12, 'cm^3/(mol*s)'), n=0,
                           Ea=(0.0, 'kcal/mol'), T0=(1.0, 'K')))
    return reaction


def _flat_nasa():
    coeffs = [2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    return NASA(polynomials=[
        NASAPolynomial(coeffs=coeffs, Tmin=(100, 'K'), Tmax=(1000, 'K')),
        NASAPolynomial(coeffs=coeffs, Tmin=(1000, 'K'), Tmax=(5000, 'K')),
    ], Tmin=(100, 'K'), Tmax=(5000, 'K'))


def _electron_counts(participants):
    return sum(1 for species in participants if species.is_electron())


class TestConservedElectronPlacementExport:
    @pytest.fixture(autouse=True)
    def declare_test_owner(self, monkeypatch):
        # pytest restores the registry entry during fixture teardown.
        monkeypatch.setitem(FAMILY_ELECTRON_PLACEMENT, OWNER, (1, 1))

    def test_declared_net_zero_places_one_electron_on_each_side(self):
        electron, argon, excited_argon = _species()
        reaction = _declared_reaction(argon, excited_argon)

        reactants, products = expand_electrons(reaction, [electron, argon, excited_argon])

        assert _electron_counts(reactants) == 1
        assert _electron_counts(products) == 1
        assert len(reactants) == 2

    def test_already_explicit_declared_pair_is_not_duplicated(self):
        electron, argon, excited_argon = _species()
        reaction = _declared_reaction(argon, excited_argon)
        reaction.reactants.append(electron)
        reaction.products.append(electron)

        reactants, products = expand_electrons(reaction, [electron, argon, excited_argon])

        assert len(reactants) == len(reaction.reactants) == 2
        assert len(products) == len(reaction.products) == 2
        assert all(actual is expected for actual, expected in
                   zip(reactants, reaction.reactants))
        assert all(actual is expected for actual, expected in
                   zip(products, reaction.products))

    def test_partial_explicit_pair_is_refused(self):
        electron, argon, excited_argon = _species()
        reaction = _declared_reaction(argon, excited_argon)
        reaction.reactants.append(electron)

        with pytest.raises(MechanismWriterError, match='explicit electron'):
            expand_electrons(reaction, [electron, argon, excited_argon])

    def test_undeclared_net_zero_reaction_is_unchanged(self):
        electron, argon, excited_argon = _species()
        reaction = Reaction(reactants=[argon], products=[excited_argon], electrons=0)

        reactants, products = expand_electrons(reaction, [electron, argon, excited_argon])

        assert reactants == [argon]
        assert products == [excited_argon]

    def test_declared_net_zero_chemkin_without_species_list_is_refused(self):
        _electron, argon, excited_argon = _species()
        reaction = _declared_reaction(argon, excited_argon)

        with pytest.raises(MechanismWriterError, match='conserved electron pair.*species_list'):
            write_reaction_string(reaction)

        with pytest.raises(MechanismWriterError, match='conserved electron pair.*species_list'):
            reaction.to_chemkin(kinetics=False)

    def test_undeclared_and_already_explicit_reactions_keep_no_list_output(self):
        electron, argon, excited_argon = _species()
        undeclared = Reaction(reactants=[argon], products=[excited_argon], electrons=0)
        assert write_reaction_string(undeclared) == 'Ar(2) <=> Ar*(3)'

        explicit = _declared_reaction(argon, excited_argon)
        explicit.reactants.append(electron)
        explicit.products.append(electron)
        assert write_reaction_string(explicit) == 'Ar(2)+e-(1)=>Ar*(3)+e-(1)'

    def test_chemkin_round_trip_keeps_conserved_pair_and_is_stable(self, tmp_path):
        from rmgpy.chemkin import load_chemkin_file

        electron, argon, excited_argon = _species()
        species = [electron, argon, excited_argon]
        for item in species:
            item.thermo = _flat_nasa()
        reaction = _declared_reaction(argon, excited_argon)
        chemkin_path = str(tmp_path / 'chem.inp')
        dictionary_path = str(tmp_path / 'species_dictionary.txt')

        save_chemkin_file(chemkin_path, species, [reaction], verbose=False,
                          check_for_duplicates=False)
        save_species_dictionary(dictionary_path, species)

        _read_species, read_reactions = load_chemkin_file(chemkin_path, dictionary_path)
        assert len(read_reactions) == 1
        read_reaction = read_reactions[0]
        assert _electron_counts(read_reaction.reactants) == 1
        assert _electron_counts(read_reaction.products) == 1

        once = expand_electrons(read_reaction, _read_species)
        twice = expand_electrons(read_reaction, _read_species)
        assert once == twice
        assert _electron_counts(once[0]) == 1
        assert _electron_counts(once[1]) == 1
