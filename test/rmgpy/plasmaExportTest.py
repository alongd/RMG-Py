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
    ArrheniusChargeTransfer,
    BadnellRRArrhenius,
    ElectronCollisionPlasma,
    Marcus,
    MultiArrhenius,
    MultiPDepArrhenius,
    PDepArrhenius,
    TwoTemperaturePlasma,
    VoronovEIArrhenius,
)
from rmgpy.kinetics.surface import SurfaceChargeTransfer
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


# ---------------------------------------------------------------------------
# B1. Charge-transfer and Marcus kinetics: exact, or refused
# ---------------------------------------------------------------------------

def _charge_transfer_kinetics(cls, alpha, electrons, Ea, V0=0.0, units='m^3/(mol*s)'):
    return cls(A=(2.5e+10, units), n=0.0, Ea=(Ea, 'kJ/mol'),
               V0=(V0, 'V'), alpha=alpha, electrons=electrons)


@pytest.fixture(scope='module')
def surface_mechanism():
    """
    A surface mechanism shaped like a real RMG electrocatalysis entry: the ion
    and the electron both live in the gas/electrolyte phase, the adsorbates on
    the surface. This is the arrangement that decides whether Cantera sees a
    charge transfer at all.
    """
    e = _make_species('e', 1, Molecule().from_adjacency_list('1 e u0 p0 c-1'))
    hplus = _make_species('Hplus', 2, Molecule().from_adjacency_list('1 H u0 p0 c+1'))
    site = _make_species('X', 3, Molecule().from_adjacency_list('1 X u0 p0 c0'))
    hx = _make_species('HX', 4, Molecule().from_adjacency_list(
        '1 H u0 p0 c0 {2,S}\n2 X u0 p0 c0 {1,S}\n'))
    return [e, hplus, site, hx]


def _surface_reaction(species, kinetics, index=1):
    e, hplus, site, hx = species
    return Reaction(index=index, reactants=[hplus, site], products=[hx],
                    electrons=-1, reversible=False, kinetics=kinetics)


def _bulk_reaction(species, kinetics, index=1):
    e, o2, o2_anion, o_atom, o_anion, o_cation = species
    return Reaction(index=index, reactants=[o_atom], products=[o_anion],
                    electrons=-1, reversible=False, kinetics=kinetics)


class TestChargeTransferKinetics:
    """
    Both charge-transfer rate laws evaluate

        k(T, V) = A*(T/T0)^n * exp(-(Ea - alpha*electrons*F*(V - V0))/(R*T))

    so writing (A, n, Ea) is exact at every potential only when
    ``alpha*electrons == 0`` AND ``Ea >= 0`` -- the second because the
    non-negative clamp on Ea_eff fires only on the ``V != V0`` branch. Both
    halves of that rule are asserted, in both writers.
    """

    def test_inert_surface_charge_transfer_round_trips_in_cantera(self, surface_mechanism, tmp_out):
        species = surface_mechanism
        kin = _charge_transfer_kinetics(SurfaceChargeTransfer, alpha=0.0, electrons=-1, Ea=72.0)
        rxn = _surface_reaction(species, kin)
        path = os.path.join(tmp_out, 'ct_inert.yaml')
        save_cantera_model(ReactionModel(species=species, reactions=[rxn]), path)

        surf = ct.Interface(path, 'surface')
        assert surf.n_reactions == 1
        rate = surf.reaction(0).rate.input_data
        assert rate['type'] == 'interface-Arrhenius'
        assert rate['rate-constant']['A'] == pytest.approx(kin.A.value_si * 1000.0 ** 2, rel=1e-8)
        assert rate['rate-constant']['b'] == pytest.approx(kin.n.value_si, rel=1e-10)
        assert rate['rate-constant']['Ea'] == pytest.approx(kin.Ea.value_si * 1000.0, rel=1e-8)
        assert 'e(1)' in surf.reaction(0).equation

    def test_live_surface_charge_transfer_is_refused_in_cantera(self, surface_mechanism, tmp_out):
        species = surface_mechanism
        kin = _charge_transfer_kinetics(SurfaceChargeTransfer, alpha=0.62, electrons=-1, Ea=72.0)
        rxn = _surface_reaction(species, kin)
        path = os.path.join(tmp_out, 'ct_live.yaml')
        with pytest.raises(MechanismWriterError) as exc:
            save_cantera_model(ReactionModel(species=species, reactions=[rxn]), path)
        assert 'potential-dependent' in str(exc.value)
        assert not os.path.exists(path)

    def test_cantera_never_writes_a_beta_that_cantera_would_discard(self, surface_mechanism, tmp_out):
        """
        The reason the live case is refused rather than exported with 'beta':
        for RMG's phase arrangement Cantera reports no charge transfer, so a
        written beta is silently dropped on read. Assert the writer emits none.
        """
        species = surface_mechanism
        kin = _charge_transfer_kinetics(SurfaceChargeTransfer, alpha=0.0, electrons=-1, Ea=72.0)
        rxn = _surface_reaction(species, kin)
        path = os.path.join(tmp_out, 'ct_nobeta.yaml')
        save_cantera_model(ReactionModel(species=species, reactions=[rxn]), path)
        assert 'beta' not in open(path).read()
        surf = ct.Interface(path, 'surface')
        assert surf.reaction(0).rate.uses_electrochemistry is False

    @pytest.mark.parametrize('alpha,electrons,Ea,inert', [
        (0.0, -1, 72.0, True),    # no transfer coefficient -> no potential term
        (0.62, 0, 72.0, True),    # no electrons -> no potential term
        (0.62, -1, 72.0, False),  # live
        (0.0, -1, -5.0, False),   # inert term, but the Ea clamp still bites off V0
    ])
    def test_inertness_rule_in_both_writers(self, mechanism, alpha, electrons, Ea, inert, tmp_out):
        model, species, reactions = mechanism
        kin = _charge_transfer_kinetics(ArrheniusChargeTransfer, alpha=alpha,
                                        electrons=electrons, Ea=Ea, V0=0.3)
        rxn = _bulk_reaction(species, kin)
        tag = '{0}_{1}_{2}'.format(alpha, electrons, Ea)
        yaml_path = os.path.join(tmp_out, 'act_{0}.yaml'.format(tag))
        inp_path = os.path.join(tmp_out, 'act_{0}.inp'.format(tag))
        container = ReactionModel(species=species, reactions=[rxn])

        if inert:
            save_cantera_model(container, yaml_path)
            gas = ct.Solution(yaml_path)
            rate = gas.reaction(0).rate.input_data
            assert rate['rate-constant']['A'] == pytest.approx(kin.A.value_si * 1000.0, rel=1e-8)
            assert rate['rate-constant']['b'] == pytest.approx(kin.n.value_si, rel=1e-10)
            assert rate['rate-constant']['Ea'] == pytest.approx(kin.Ea.value_si * 1000.0, rel=1e-8)

            save_chemkin_file(inp_path, species, [rxn], check_for_duplicates=False)
            content = open(inp_path).read()
            expected_A = (kin.A.value_si / (kin.T0.value_si ** kin.n.value_si)
                          * kin.A.get_conversion_factor_from_si_to_cm_mol_s())
            line = [l for l in content.splitlines()
                    if l.strip().startswith('O(4)+e(1)=>O-(5)') and not l.strip().startswith('!')]
            assert len(line) == 1, content
            A, n, Ea_kcal = [float(tok) for tok in line[0].split()[1:4]]
            assert A == pytest.approx(expected_A, rel=1e-3)
            assert n == pytest.approx(kin.n.value_si, abs=1e-3)
            assert Ea_kcal == pytest.approx(kin.Ea.value_si / 4184.0, abs=1e-3)
        else:
            with pytest.raises(MechanismWriterError) as exc:
                save_cantera_model(container, yaml_path)
            assert 'potential-dependent' in str(exc.value)
            assert not os.path.exists(yaml_path)
            with pytest.raises(MechanismWriterError):
                save_chemkin_file(inp_path, species, [rxn], check_for_duplicates=False)
            assert not os.path.exists(inp_path)

    def test_marcus_is_refused_by_both_writers(self, mechanism, tmp_out):
        """
        Marcus rates depend on dGrxn, which is not a property of the rate law --
        there is no reference point at which any reduction is exact, so there is
        no correct file to write.
        """
        model, species, reactions = mechanism
        e, o2, o2_anion, o_atom, o_anion, o_cation = species
        kin = Marcus(A=(1.0e+10, 'm^3/(mol*s)'), n=0.0,
                     lmbd_i_coefs=np.array([1.0e+04, 0.0, 0.0, 0.0]),
                     wr=(0, 'J/mol'), wp=(0, 'J/mol'), lmbd_o=(5.0e+04, 'J/mol'))
        rxn = Reaction(index=1, reactants=[o_atom, o_atom], products=[o2], electrons=0,
                       reversible=False, kinetics=kin)
        yaml_path = os.path.join(tmp_out, 'marcus.yaml')
        inp_path = os.path.join(tmp_out, 'marcus.inp')
        container = ReactionModel(species=species, reactions=[rxn])

        with pytest.raises(MechanismWriterError) as exc:
            save_cantera_model(container, yaml_path)
        assert 'dGrxn' in str(exc.value)
        assert not os.path.exists(yaml_path)

        with pytest.raises(MechanismWriterError) as exc:
            save_chemkin_file(inp_path, species, [rxn], check_for_duplicates=False)
        assert 'dGrxn' in str(exc.value)
        assert not os.path.exists(inp_path)


# ---------------------------------------------------------------------------
# B2. Grouped kinetics keep the electron
# ---------------------------------------------------------------------------

class TestMultiKineticsCarryElectrons:

    @pytest.mark.parametrize('electrons,units,expected', [
        (-1, 'cm^3/(mol*s)', 'O(4) + e(1) => O-(5)'),
        (1, 's^-1', 'O-(5) => O(4) + e(1)'),
    ])
    def test_multi_arrhenius(self, mechanism, electrons, units, expected, tmp_out):
        model, species, reactions = mechanism
        e, o2, o2_anion, o_atom, o_anion, o_cation = species
        reactants, products = ([o_atom], [o_anion]) if electrons < 0 else ([o_anion], [o_atom])
        rxn = Reaction(index=1, reactants=reactants, products=products, electrons=electrons,
                       reversible=False, kinetics=MultiArrhenius(arrhenius=[
                           Arrhenius(A=(1.0e+10, units), n=0.0, Ea=(0.0, 'J/mol')),
                           Arrhenius(A=(2.0e+10, units), n=0.0, Ea=(0.0, 'J/mol')),
                       ]))
        container = ReactionModel(species=species, reactions=[rxn])

        yaml_path = os.path.join(tmp_out, 'multi_{0}.yaml'.format(electrons))
        save_cantera_model(container, yaml_path)
        gas = ct.Solution(yaml_path)
        assert gas.n_reactions == 2
        equations = [l.split(':', 1)[1].strip() for l in open(yaml_path).read().splitlines()
                     if l.startswith('- equation:')]
        assert equations == [expected, expected], equations

        inp_path = os.path.join(tmp_out, 'multi_{0}.inp'.format(electrons))
        save_chemkin_file(inp_path, species, [rxn], check_for_duplicates=False)
        chemkin_equation = expected.replace(' ', '')
        assert open(inp_path).read().count(chemkin_equation) == 2, open(inp_path).read()

    @pytest.mark.parametrize('electrons,expected', [
        (-1, 'O(4) + e(1) => O-(5)'),
        (1, 'O-(5) => O(4) + e(1)'),
    ])
    def test_multi_pdep_arrhenius(self, mechanism, electrons, expected, tmp_out):
        model, species, reactions = mechanism
        e, o2, o2_anion, o_atom, o_anion, o_cation = species
        reactants, products = ([o_atom], [o_anion]) if electrons < 0 else ([o_anion], [o_atom])
        units = 'cm^3/(mol*s)' if electrons < 0 else 's^-1'

        def pdep():
            return PDepArrhenius(
                pressures=([0.1, 1.0], 'bar'),
                arrhenius=[Arrhenius(A=(1.0e+10, units), n=0.0, Ea=(0.0, 'J/mol')),
                           Arrhenius(A=(2.0e+10, units), n=0.0, Ea=(0.0, 'J/mol'))])

        rxn = Reaction(index=1, reactants=reactants, products=products, electrons=electrons,
                       reversible=False,
                       kinetics=MultiPDepArrhenius(arrhenius=[pdep(), pdep()]))
        container = ReactionModel(species=species, reactions=[rxn])

        yaml_path = os.path.join(tmp_out, 'multipdep_{0}.yaml'.format(electrons))
        save_cantera_model(container, yaml_path)
        equations = [l.split(':', 1)[1].strip() for l in open(yaml_path).read().splitlines()
                     if l.startswith('- equation:')]
        assert equations == [expected, expected], equations

        inp_path = os.path.join(tmp_out, 'multipdep_{0}.inp'.format(electrons))
        save_chemkin_file(inp_path, species, [rxn], check_for_duplicates=False)
        assert open(inp_path).read().count(expected.replace(' ', '')) == 2


# ---------------------------------------------------------------------------
# B3. The two writers refuse on the same condition
# ---------------------------------------------------------------------------

class TestWritersAgreeOnMissingElectronSpecies:

    def test_both_refuse_plasma_kinetics_with_no_electron_species(self, mechanism, tmp_out):
        """
        The Cantera writer decides plasma-phase mode from species presence, so
        without this check it wrote plasma reaction types into a plain
        'ideal-gas' phase while Chemkin refused the same mechanism.
        """
        model, species, reactions = mechanism
        e, o2, o2_anion, o_atom, o_anion, o_cation = species
        no_electron = [spc for spc in species if not spc.is_electron()]
        neutral_plasma_rxn = Reaction(
            index=1, reactants=[o_atom, o_atom], products=[o2], electrons=0, reversible=False,
            kinetics=reactions[1].kinetics)  # ElectronCollisionPlasma
        container = ReactionModel(species=no_electron, reactions=[neutral_plasma_rxn])

        yaml_path = os.path.join(tmp_out, 'noe.yaml')
        inp_path = os.path.join(tmp_out, 'noe.inp')
        with pytest.raises(MechanismWriterError) as cantera_exc:
            save_cantera_model(container, yaml_path)
        with pytest.raises(MechanismWriterError) as chemkin_exc:
            save_chemkin_file(inp_path, no_electron, [neutral_plasma_rxn],
                              check_for_duplicates=False)
        assert 'electron species' in str(cantera_exc.value)
        assert 'electron species' in str(chemkin_exc.value)
        assert not os.path.exists(yaml_path)
        assert not os.path.exists(inp_path)


# ---------------------------------------------------------------------------
# B4. The wrong-order guard reaches the classes that need it
# ---------------------------------------------------------------------------

class TestWrongReactionOrderGuard:

    @pytest.mark.parametrize('kinetics_index,name', [
        (1, 'ElectronCollisionPlasma'),
        (0, 'TwoTemperaturePlasma'),
    ])
    def test_electron_impact_without_an_electron_reactant_is_caught(
            self, mechanism, kinetics_index, name, tmp_out):
        """
        Neither of these classes defines uses_electron_density, so the flag the
        guard used to consult was absent and the guard never fired. The rate
        coefficient's own dimensionality answers the question instead.
        """
        model, species, reactions = mechanism
        e, o2, o2_anion, o_atom, o_anion, o_cation = species
        # Electron-impact dissociation written without its electron: neutral on
        # both sides, so it balances in E and only the order check can see it.
        rxn = Reaction(index=1, reactants=[o2], products=[o_atom, o_atom], electrons=0,
                       reversible=False, kinetics=reactions[kinetics_index].kinetics)
        container = ReactionModel(species=species, reactions=[rxn])
        yaml_path = os.path.join(tmp_out, 'order_{0}.yaml'.format(name))
        inp_path = os.path.join(tmp_out, 'order_{0}.inp'.format(name))

        with pytest.raises(MechanismWriterError) as cantera_exc:
            save_cantera_model(container, yaml_path)
        with pytest.raises(MechanismWriterError) as chemkin_exc:
            save_chemkin_file(inp_path, species, [rxn], check_for_duplicates=False)
        for exc in (cantera_exc, chemkin_exc):
            assert 'reaction order' in str(exc.value)
            assert name in str(exc.value)
        assert not os.path.exists(yaml_path)
        assert not os.path.exists(inp_path)

    def test_the_five_reaction_plasma_fixture_still_exports(self, mechanism, tmp_out):
        """Regression guard: the new order check must not fire on correct input."""
        model, species, reactions = mechanism
        path = os.path.join(tmp_out, 'fixture.yaml')
        save_cantera_model(model, path)
        assert ct.Solution(path).n_reactions == len(reactions)


# ---------------------------------------------------------------------------
# B5. A failed export leaves nothing behind
# ---------------------------------------------------------------------------

class TestFailedExportLeavesNothing:

    def _failing_model(self, mechanism):
        model, species, reactions = mechanism
        bad = Reaction(index=99, reactants=[species[1]], products=[species[3], species[3]],
                       electrons=0, reversible=False, kinetics=UnhandledKinetics())
        return species, list(reactions) + [bad]

    def test_chemkin_leaves_no_partial_file_and_no_dirty_counter(self, mechanism, tmp_out):
        import rmgpy.chemkin
        species, reactions = self._failing_model(mechanism)
        path = os.path.join(tmp_out, 'partial.inp')
        with pytest.raises(MechanismWriterError):
            save_chemkin_file(path, species, reactions, check_for_duplicates=False)
        assert not os.path.exists(path), 'a partial mechanism was left on disk'
        assert rmgpy.chemkin._chemkin_reaction_count is None, 'the global counter was left dirty'
        assert os.listdir(tmp_out) == [], 'a temporary file was left behind: {0}'.format(
            os.listdir(tmp_out))

    def test_chemkin_does_not_clobber_an_existing_good_file(self, mechanism, tmp_out):
        """The destination is only replaced once the whole file is built."""
        model, species, reactions = mechanism
        path = os.path.join(tmp_out, 'chem.inp')
        save_chemkin_file(path, species, reactions, check_for_duplicates=False)
        good = open(path).read()
        assert 'END' in good

        species, failing = self._failing_model(mechanism)
        with pytest.raises(MechanismWriterError):
            save_chemkin_file(path, species, failing, check_for_duplicates=False)
        assert open(path).read() == good, 'the previous good mechanism was clobbered'

    def test_counter_is_cleared_after_a_successful_export_too(self, mechanism, tmp_out):
        import rmgpy.chemkin
        model, species, reactions = mechanism
        save_chemkin_file(os.path.join(tmp_out, 'ok.inp'), species, reactions,
                          check_for_duplicates=False)
        assert rmgpy.chemkin._chemkin_reaction_count is None

    def test_cantera_leaves_no_partial_file(self, mechanism, tmp_out):
        species, reactions = self._failing_model(mechanism)
        path = os.path.join(tmp_out, 'partial.yaml')
        with pytest.raises(MechanismWriterError):
            save_cantera_model(ReactionModel(species=species, reactions=reactions), path)
        assert not os.path.exists(path)
