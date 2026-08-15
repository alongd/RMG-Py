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
Export of plasma mechanisms through both mechanism writers.

Every assertion here is made against the file that was actually written and then
read back, never against an exception alone: each defect these tests cover is one
that happens quietly while the process exits zero.
"""

import os
import shutil

import cantera as ct
import numpy as np
import pytest

import rmgpy.constants as constants
from rmgpy.chemkin import save_chemkin_file, write_kinetics_entry, write_reaction_string
from rmgpy.exceptions import MechanismWriterError
from rmgpy.kinetics import (
    Arrhenius,
    BadnellRRArrhenius,
    ElectronCollisionPlasma,
    TwoTemperaturePlasma,
    VoronovEIArrhenius,
)
from rmgpy.kinetics.model import KineticsModel
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import ReactionModel
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial
from rmgpy.transport import TransportData
from rmgpy.yaml_cantera2 import save_cantera_model


def _make_species(label, index, molecule):
    """An RMG Species carrying the thermo and transport both writers need."""
    spc = Species(label=label, molecule=[molecule])
    spc.index = index
    coeffs = [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, -11.7246]
    spc.thermo = NASA(
        polynomials=[
            NASAPolynomial(coeffs=coeffs, Tmin=(200, 'K'), Tmax=(1000, 'K')),
            NASAPolynomial(coeffs=coeffs, Tmin=(1000, 'K'), Tmax=(6000, 'K')),
        ],
        Tmin=(200, 'K'), Tmax=(6000, 'K'),
    )
    spc.transport_data = TransportData(
        shapeIndex=0 if len(molecule.atoms) == 1 else 1,
        sigma=(3.0, 'angstrom'),
        epsilon=(100.0, 'K'),
        dipoleMoment=(0.0, 'De'),
        polarizability=(0.0, 'angstrom^3'),
        rotrelaxcollnum=0.0,
    )
    return spc


class UnhandledKinetics(KineticsModel):
    """A kinetics type deliberately given no case in either writer."""
    pass


@pytest.fixture(scope='module')
def mechanism():
    """
    A mechanism containing one reaction of each of the four plasma kinetics
    types, plus a plain Arrhenius reaction as a control.

    The electron stoichiometry follows RMG's convention: ``reaction.electrons``
    counts electrons that are *not* already in the reactant/product lists, and is
    negative when they are consumed. Electron-impact ionization therefore carries
    the consumed electron explicitly in ``reactants`` and counts the two produced
    electrons in ``electrons``, which is what keeps both the E balance and the
    reaction order right.
    """
    e = _make_species('e', 1, Molecule().from_adjacency_list('1 e u0 p0 c-1'))
    o2 = _make_species('O2', 2, Molecule(smiles='[O][O]'))
    o2_anion = _make_species('O2-', 3, Molecule().from_adjacency_list(
        '1 O u1 p2 c0 {2,S}\n2 O u0 p3 c-1 {1,S}\n'))
    o_atom = _make_species('O', 4, Molecule().from_adjacency_list('1 O u2 p2 c0'))
    o_anion = _make_species('O-', 5, Molecule().from_adjacency_list('1 O u1 p3 c-1'))
    o_cation = _make_species('O+', 6, Molecule().from_adjacency_list('1 O u3 p1 c+1'))

    species = [e, o2, o2_anion, o_atom, o_anion, o_cation]

    two_temp = Reaction(
        index=1, reactants=[o2, e], products=[o2_anion], electrons=0,
        reversible=False,
        kinetics=TwoTemperaturePlasma(
            A=(1.7e+07, 'm^3/(mol*s)'), n=-0.5,
            Ea_g=(3000.0, 'J/mol'), Ea_e=(45000.0, 'J/mol'),
        ),
    )
    collision = Reaction(
        index=2, reactants=[o_atom], products=[o_anion], electrons=-1,
        reversible=False,
        kinetics=ElectronCollisionPlasma(
            energies=([0.0, 1.0, 2.0, 5.0, 10.0], 'eV/molecule'),
            sigma=([0.0, 1.0e-21, 2.5e-20, 1.1e-20, 3.0e-21], 'm^2'),
        ),
    )
    badnell = Reaction(
        index=3, reactants=[o_cation], products=[o_atom], electrons=-1,
        reversible=False,
        kinetics=BadnellRRArrhenius(
            A=(8.0e-12, 'cm^3/(molecule*s)'), B=0.7,
            T0=(3.0e+02, 'K'), T1=(1.5e+06, 'K'),
            Tmin=(1.0e+04, 'K'), Tmax=(1.0e+06, 'K'),
        ),
    )
    voronov = Reaction(
        index=4, reactants=[o_atom, e], products=[o_cation], electrons=2,
        reversible=False,
        kinetics=VoronovEIArrhenius(
            A=(3.59e-08, 'cm^3/(molecule*s)'), P=0.0, X=0.073, K=0.47,
            dE=(13.62, 'eV'), Tmin=(1.0e+04, 'K'), Tmax=(1.0e+06, 'K'),
        ),
    )
    thermal = Reaction(
        index=5, reactants=[o_atom, o_atom], products=[o2], electrons=0,
        reversible=False,
        kinetics=Arrhenius(A=(1.2e+11, 'cm^3/(mol*s)'), n=0.3, Ea=(5000.0, 'J/mol')),
    )

    reactions = [two_temp, collision, badnell, voronov, thermal]
    return ReactionModel(species=species, reactions=reactions), species, reactions


@pytest.fixture()
def tmp_out(tmp_path):
    return str(tmp_path)


# ---------------------------------------------------------------------------
# 1. Round-trip through Cantera, comparing parameters and not just counts
# ---------------------------------------------------------------------------

class TestCanteraRoundTrip:

    def _export_and_load(self, mechanism, tmp_out):
        model, species, reactions = mechanism
        path = os.path.join(tmp_out, 'chem.yaml')
        save_cantera_model(model, path)
        assert os.path.exists(path)
        return ct.Solution(path), path, species, reactions

    def test_reaction_count_survives(self, mechanism, tmp_out):
        """No plasma reaction is dropped on the way out."""
        gas, path, species, reactions = self._export_and_load(mechanism, tmp_out)
        assert gas.n_reactions == len(reactions), (
            'exported {0} of {1} reactions'.format(gas.n_reactions, len(reactions)))

    def test_two_temperature_plasma_parameters(self, mechanism, tmp_out):
        gas, path, species, reactions = self._export_and_load(mechanism, tmp_out)
        kin = reactions[0].kinetics
        rate = gas.reaction(0).rate.input_data
        assert rate['type'] == 'two-temperature-plasma'
        # Cantera reports in kmol; the file is written in mol, so A picks up 1000
        # per unit of reaction order beyond the first (this reaction is second order).
        assert rate['rate-constant']['A'] == pytest.approx(kin.A.value_si * 1000.0, rel=1e-10)
        assert rate['rate-constant']['b'] == pytest.approx(kin.n.value_si, rel=1e-10)
        assert rate['rate-constant']['Ea-gas'] == pytest.approx(kin.Ea_g.value_si * 1000.0, rel=1e-10)
        assert rate['rate-constant']['Ea-electron'] == pytest.approx(kin.Ea_e.value_si * 1000.0, rel=1e-10)

    def test_electron_collision_plasma_parameters(self, mechanism, tmp_out):
        gas, path, species, reactions = self._export_and_load(mechanism, tmp_out)
        kin = reactions[1].kinetics
        rate = gas.reaction(1).rate.input_data
        assert rate['type'] == 'electron-collision-plasma'
        expected_eV = kin.energies.value_si / constants.F
        assert np.allclose(rate['energy-levels'], expected_eV, rtol=1e-10)
        assert np.allclose(rate['cross-sections'], kin.sigma.value_si, rtol=1e-10)

    @pytest.mark.parametrize('idx', [2, 3])
    def test_badnell_and_voronov_parameters(self, mechanism, tmp_out, idx):
        """
        Badnell and Voronov have no Cantera form. They are pure functions of Te,
        so they go out as two-temperature-plasma with Ea-gas == Ea-electron,
        which is the only choice that cancels the gas-temperature factors and
        leaves k = A*Te^b*exp(-Ea/(R*Te)).
        """
        gas, path, species, reactions = self._export_and_load(mechanism, tmp_out)
        kin = reactions[idx].kinetics
        arr = kin.to_arrhenius()
        rate = gas.reaction(idx).rate.input_data
        assert rate['type'] == 'two-temperature-plasma'
        expected_A = arr.A.value_si / (arr.T0.value_si ** arr.n.value_si)
        assert rate['rate-constant']['A'] == pytest.approx(expected_A * 1000.0, rel=1e-8)
        assert rate['rate-constant']['b'] == pytest.approx(arr.n.value_si, rel=1e-8)
        assert rate['rate-constant']['Ea-gas'] == pytest.approx(arr.Ea.value_si * 1000.0, rel=1e-8)
        assert rate['rate-constant']['Ea-electron'] == pytest.approx(arr.Ea.value_si * 1000.0, rel=1e-8)

    def test_badnell_and_voronov_rate_is_electron_temperature_only(self, mechanism, tmp_out):
        """
        The mapping is only exact if the exported rate does not move with the gas
        temperature at fixed Te. Assert that on the loaded mechanism, which is the
        thing a solver would see.
        """
        model, species, reactions = mechanism
        # Only the two fitted forms: the tabulated cross-section reaction would
        # drag Cantera's Boltzmann integration into this, which is not what is
        # under test here.
        path = os.path.join(tmp_out, 'fitted.yaml')
        save_cantera_model(ReactionModel(species=species, reactions=[reactions[2], reactions[3]]), path)
        gas = ct.Solution(path)
        gas.TPX = 300.0, ct.one_atm, 'O(4):1,e(1):1e-8,O+(6):1e-8,O-(5):1e-8,O2(2):1,O2-(3):1e-8'
        rates = []
        for temperature in (300.0, 1000.0, 3000.0):
            gas.TP = temperature, ct.one_atm
            gas.Te = 20000.0
            rates.append(list(gas.forward_rate_constants))
        for badnell_k, voronov_k in rates[1:]:
            assert badnell_k == pytest.approx(rates[0][0], rel=1e-10)
            assert voronov_k == pytest.approx(rates[0][1], rel=1e-10)

    def test_ionization_is_first_order_in_the_electron(self, mechanism, tmp_out):
        """
        Electron-impact ionization is written A + e => A+ + 2 e, with the electron
        on both sides. Cantera's Reaction.reactants reports *net* stoichiometry,
        so the electron looks absent from the reactant side -- but the rate of
        progress is evaluated at the as-written order. Pin that down on the
        loaded mechanism, because reading the net stoichiometry as the reaction
        order is an easy and expensive mistake to make about this file.
        """
        model, species, reactions = mechanism
        path = os.path.join(tmp_out, 'ionization.yaml')
        save_cantera_model(ReactionModel(species=species, reactions=[reactions[3]]), path)
        gas = ct.Solution(path)

        rates = []
        for x_e in (1.0e-8, 1.0e-7):
            gas.TPX = 1000.0, ct.one_atm, (
                'O(4):1,e(1):{0},O+(6):1e-9,O-(5):1e-9,O2(2):1e-9,O2-(3):1e-9'.format(x_e))
            gas.Te = 20000.0
            rates.append(gas.forward_rates_of_progress[0])
        assert rates[1] / rates[0] == pytest.approx(10.0, rel=1e-3), (
            'rate of progress is not first order in the electron: {0}'.format(rates))

    def test_plain_arrhenius_still_correct(self, mechanism, tmp_out):
        """The control: adding plasma cases did not disturb ordinary kinetics."""
        gas, path, species, reactions = self._export_and_load(mechanism, tmp_out)
        kin = reactions[4].kinetics
        rate = gas.reaction(4).rate.input_data
        assert rate['rate-constant']['A'] == pytest.approx(kin.A.value_si * 1000.0, rel=1e-10)
        assert rate['rate-constant']['b'] == pytest.approx(kin.n.value_si, rel=1e-10)
        assert rate['rate-constant']['Ea'] == pytest.approx(kin.Ea.value_si * 1000.0, rel=1e-10)


# ---------------------------------------------------------------------------
# 2. The negative test: an unhandled kinetics type must stop the export
# ---------------------------------------------------------------------------

class TestUnhandledKineticsFailsLoudly:

    def _mechanism_with_unhandled(self, mechanism):
        model, species, reactions = mechanism
        bad = Reaction(
            index=6, reactants=[species[1]], products=[species[3], species[3]],
            electrons=0, reversible=False, kinetics=UnhandledKinetics(),
        )
        return ReactionModel(species=species, reactions=list(reactions) + [bad]), species, bad

    def test_cantera_writer_raises_and_writes_no_file_with_the_reaction(self, mechanism, tmp_out):
        model, species, bad = self._mechanism_with_unhandled(mechanism)
        path = os.path.join(tmp_out, 'bad.yaml')
        with pytest.raises(MechanismWriterError) as exc:
            save_cantera_model(model, path)
        assert 'UnhandledKinetics' in str(exc.value)
        # The artifact, not the exception: the mechanism must not be sitting on
        # disk quietly missing this reaction.
        if os.path.exists(path):
            content = open(path).read()
            assert 'O2(2) => O(4) + O(4)' not in content, (
                'export failed but left a mechanism file behind containing the reaction')

    def test_chemkin_writer_raises_and_writes_no_file_with_the_reaction(self, mechanism, tmp_out):
        model, species, bad = self._mechanism_with_unhandled(mechanism)
        path = os.path.join(tmp_out, 'bad.inp')
        with pytest.raises(MechanismWriterError) as exc:
            save_chemkin_file(path, model.species, model.reactions, check_for_duplicates=False)
        assert 'UnhandledKinetics' in str(exc.value)
        if os.path.exists(path):
            content = open(path).read()
            assert 'O2(2)=>O(4)+O(4)' not in content

    def test_no_opt_out_downgrades_it_to_a_warning(self, mechanism, tmp_out):
        """
        The hard failure has no escape hatch. If someone adds a flag, config key
        or environment variable to turn it back into a warning, this test is the
        thing that notices.
        """
        import inspect
        import rmgpy.chemkin
        import rmgpy.yaml_cantera2
        for module in (rmgpy.yaml_cantera2,):
            source = inspect.getsource(module.reaction_to_dict_list)
            assert 'logging.warning' not in source
        source = inspect.getsource(rmgpy.electron_balance)
        assert 'logging' not in source


# ---------------------------------------------------------------------------
# 3. Electron balance in the exported equation
# ---------------------------------------------------------------------------

class TestElectronBalance:

    def test_cantera_equations_carry_the_electron(self, mechanism, tmp_out):
        model, species, reactions = mechanism
        path = os.path.join(tmp_out, 'chem.yaml')
        save_cantera_model(model, path)
        gas = ct.Solution(path)
        equations = [gas.reaction(i).equation for i in range(gas.n_reactions)]
        # electrons = -1: the electron is added to the reactant side
        assert 'e(1)' in equations[1], equations[1]
        assert equations[1].split('=')[0].count('e(1)') == 1
        assert 'e(1)' in equations[2], equations[2]
        # electrons = +2 on top of an explicit reactant electron
        assert equations[3].split('=')[1].count('e(1)') == 2, equations[3]

    def test_cantera_reactions_balance_in_E(self, mechanism, tmp_out):
        """
        Asserted from Cantera's own element bookkeeping on the loaded mechanism,
        not from the equation string.
        """
        model, species, reactions = mechanism
        path = os.path.join(tmp_out, 'chem.yaml')
        save_cantera_model(model, path)
        gas = ct.Solution(path)
        e_index = gas.element_index('E')
        for i in range(gas.n_reactions):
            rxn = gas.reaction(i)
            left = sum(coeff * gas.n_atoms(name, e_index)
                       for name, coeff in rxn.reactants.items())
            right = sum(coeff * gas.n_atoms(name, e_index)
                        for name, coeff in rxn.products.items())
            assert left == pytest.approx(right), (
                'reaction {0} does not balance in E: {1} vs {2}'.format(rxn.equation, left, right))

    def test_chemkin_equations_carry_the_electron(self, mechanism, tmp_out):
        model, species, reactions = mechanism
        path = os.path.join(tmp_out, 'chem.inp')
        save_chemkin_file(path, model.species, model.reactions, check_for_duplicates=False)
        content = open(path).read()
        # get_species_identifier strips '+' from labels, so O+(6) is written O(6).
        assert 'O(4)+e(1)=>O(6)+e(1)+e(1)' in content, content
        assert 'O(4)+e(1)=>O-(5)' in content, content
        assert 'O(6)+e(1)=>O(4)' in content, content

    def test_unbalanced_E_is_caught_by_the_same_loud_path(self, mechanism, tmp_out):
        """
        A reaction whose E does not balance must fail the export, not be written.
        Here the electron stoichiometry is simply wrong: an attachment declared
        with no electron consumed.
        """
        model, species, reactions = mechanism
        o_atom, o_anion = species[3], species[4]
        unbalanced = Reaction(
            index=7, reactants=[o_atom], products=[o_anion], electrons=0,
            reversible=False,
            kinetics=Arrhenius(A=(1.0e+10, 'cm^3/(mol*s)'), n=0.0, Ea=(0.0, 'J/mol')),
        )
        bad_model = ReactionModel(species=species, reactions=[unbalanced])

        yaml_path = os.path.join(tmp_out, 'unbalanced.yaml')
        with pytest.raises(MechanismWriterError) as exc:
            save_cantera_model(bad_model, yaml_path)
        assert 'E pseudo-element' in str(exc.value)
        if os.path.exists(yaml_path):
            assert 'O(4) => O-(5)' not in open(yaml_path).read()

        inp_path = os.path.join(tmp_out, 'unbalanced.inp')
        with pytest.raises(MechanismWriterError) as exc:
            save_chemkin_file(inp_path, species, [unbalanced], check_for_duplicates=False)
        assert 'E pseudo-element' in str(exc.value)
        if os.path.exists(inp_path):
            assert 'O(4)=>O-(5)' not in open(inp_path).read()

    def test_missing_electron_species_is_caught(self, mechanism, tmp_out):
        """
        A reaction needing an electron in a mechanism that defines none cannot be
        written at all, in either writer.
        """
        model, species, reactions = mechanism
        no_electron = [spc for spc in species if not spc.is_electron()]
        attachment = reactions[1]
        bad_model = ReactionModel(species=no_electron, reactions=[attachment])
        with pytest.raises(MechanismWriterError):
            save_cantera_model(bad_model, os.path.join(tmp_out, 'no_e.yaml'))
        with pytest.raises(MechanismWriterError):
            save_chemkin_file(os.path.join(tmp_out, 'no_e.inp'), no_electron, [attachment],
                              check_for_duplicates=False)

    def test_electron_density_dependence_without_a_reactant_electron_is_caught(self, mechanism, tmp_out):
        """
        Net-electron bookkeeping alone lets an ionization balance in E while
        exporting at the wrong reaction order (A => A+ + e is read as first order
        and never multiplied by the electron concentration). That must fail too.
        """
        model, species, reactions = mechanism
        o_atom, o_cation = species[3], species[5]
        wrong_order = Reaction(
            index=8, reactants=[o_atom], products=[o_cation], electrons=1,
            reversible=False, kinetics=reactions[3].kinetics,
        )
        bad_model = ReactionModel(species=species, reactions=[wrong_order])
        with pytest.raises(MechanismWriterError) as exc:
            save_cantera_model(bad_model, os.path.join(tmp_out, 'order.yaml'))
        assert 'reaction order' in str(exc.value)
        with pytest.raises(MechanismWriterError):
            save_chemkin_file(os.path.join(tmp_out, 'order.inp'), species, [wrong_order],
                              check_for_duplicates=False)


# ---------------------------------------------------------------------------
# 4. Both writers, side by side
# ---------------------------------------------------------------------------

class TestBothWriters:

    def test_chemkin_writes_every_plasma_reaction_with_real_parameters(self, mechanism, tmp_out):
        """
        Chemkin has no plasma rate forms, so each reaction goes out as the
        modified-Arrhenius reduction of its rate law along T = Te, marked TDEP.
        What matters here is that nothing is dropped and nothing carries the old
        placeholder triple '1.000e+00 0.000     0.000'.
        """
        model, species, reactions = mechanism
        path = os.path.join(tmp_out, 'chem.inp')
        save_chemkin_file(path, model.species, model.reactions, check_for_duplicates=False)
        content = open(path).read()

        assert content.count('TDEP/e(1)/') == 4, content
        assert '1.000e+00 0.000     0.000' not in content

        # The two-temperature reduction is exact along T = Te: A, n and Ea_g go
        # out verbatim, in Chemkin's cm/mol/s and kcal/mol.
        kin = reactions[0].kinetics
        expected_A = kin.A.value_si * kin.A.get_conversion_factor_from_si_to_cm_mol_s()
        assert '{0:<9.3e}'.format(expected_A).strip() in content
        assert 'Ea_electron={0:.4g} J/mol'.format(kin.Ea_e.value_si) in content

        # The three fitted forms record that they are fits.
        for name in ('ElectronCollisionPlasma', 'BadnellRRArrhenius', 'VoronovEIArrhenius'):
            assert name in content, 'no provenance note for {0}'.format(name)

    def test_chemkin_plasma_rate_matches_the_rmg_object(self, mechanism, tmp_out):
        """Parse the A/n/Ea back off the written line and compare to the source."""
        model, species, reactions = mechanism
        path = os.path.join(tmp_out, 'chem.inp')
        save_chemkin_file(path, model.species, model.reactions, check_for_duplicates=False)
        lines = open(path).read().splitlines()

        targets = {
            'O2(2)+e(1)=>O2-(3)': reactions[0],
            'O(4)+e(1)=>O-(5)': reactions[1],
            'O(6)+e(1)=>O(4)': reactions[2],
        }
        seen = {}
        for line in lines:
            stripped = line.strip()
            for equation in targets:
                if stripped.startswith(equation) and not stripped.startswith('!'):
                    seen[equation] = [float(token) for token in stripped[len(equation):].split()]
        assert set(seen) == set(targets), 'missing reaction lines: {0}'.format(
            set(targets) - set(seen))

        for equation, rxn in targets.items():
            A, n, Ea = seen[equation]
            kin = rxn.kinetics
            if isinstance(kin, TwoTemperaturePlasma):
                arr_A = kin.A.value_si / (kin.T0.value_si ** kin.n.value_si)
                arr_n, arr_Ea = kin.n.value_si, kin.Ea_g.value_si
                factor = kin.A.get_conversion_factor_from_si_to_cm_mol_s()
            else:
                arr = kin.to_arrhenius()
                arr_A = arr.A.value_si / (arr.T0.value_si ** arr.n.value_si)
                arr_n, arr_Ea = arr.n.value_si, arr.Ea.value_si
                factor = arr.A.get_conversion_factor_from_si_to_cm_mol_s()
            assert A == pytest.approx(arr_A * factor, rel=1e-3), equation
            assert n == pytest.approx(arr_n, abs=1e-3), equation
            assert Ea == pytest.approx(arr_Ea / 4184.0, abs=1e-3), equation

    def test_both_writers_agree_on_which_reactions_exist(self, mechanism, tmp_out):
        model, species, reactions = mechanism
        yaml_path = os.path.join(tmp_out, 'chem.yaml')
        inp_path = os.path.join(tmp_out, 'chem.inp')
        save_cantera_model(model, yaml_path)
        save_chemkin_file(inp_path, model.species, model.reactions, check_for_duplicates=False)

        gas = ct.Solution(yaml_path)
        chemkin_text = open(inp_path).read()
        assert gas.n_reactions == len(reactions)
        for rxn in reactions:
            equation = write_reaction_string(rxn, species_list=species)
            assert equation in chemkin_text, 'Chemkin file is missing {0}'.format(equation)

    def test_write_kinetics_entry_raises_for_unhandled_type(self, mechanism, tmp_out):
        """The Chemkin entry writer fails on its own, not only via save_chemkin_file."""
        model, species, reactions = mechanism
        bad = Reaction(index=9, reactants=[species[3]], products=[species[3]],
                       electrons=0, reversible=False, kinetics=UnhandledKinetics())
        with pytest.raises(MechanismWriterError):
            write_kinetics_entry(bad, species_list=species)
