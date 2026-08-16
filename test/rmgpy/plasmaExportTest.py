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
from rmgpy.chemkin import (save_chemkin, save_chemkin_file, save_chemkin_surface_file,
                           write_kinetics_entry, write_reaction_string)
from rmgpy.exceptions import MechanismWriterError
from rmgpy.kinetics import (
    Arrhenius,
    ArrheniusChargeTransfer,
    BadnellRRArrhenius,
    ElectronCollisionPlasma,
    Lindemann,
    Marcus,
    MultiArrhenius,
    MultiPDepArrhenius,
    PDepArrhenius,
    ThirdBody,
    Troe,
    TwoTemperaturePlasma,
    VoronovEIArrhenius,
)
from rmgpy.kinetics.surface import StickingCoefficient, SurfaceArrhenius, SurfaceChargeTransfer
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


# ---------------------------------------------------------------------------
# 9. Reference temperature: the two writers must agree on the *rate*
# ---------------------------------------------------------------------------

#: A reference temperature that is not RMG's conventional 1 K, together with a
#: temperature exponent that is not zero, so that ``T0**n`` is a factor of 9e4 --
#: large enough that no rounding can hide it.
_T0 = (300.0, 'K')
_N = 2.0
_EA = (10.0, 'kJ/mol')
#: Site density shared by both writers. Their defaults differ (2.72e-9 mol/cm^2
#: in the Chemkin writer, 2.5e-5 mol/m^2 here), which would otherwise show up as
#: a rate difference that has nothing to do with the reference temperature.
_SITE_DENSITY = 2.5e-5  # mol/m^2


@pytest.fixture(scope='module')
def gas_forms_species():
    """Neutral and charged gas species enough to build one of every gas rate form."""
    e = _make_species('e', 1, Molecule().from_adjacency_list('1 e u0 p0 c-1'))
    o_atom = _make_species('O', 2, Molecule().from_adjacency_list('1 O u2 p2 c0'))
    o2 = _make_species('O2', 3, Molecule(smiles='[O][O]'))
    o3 = _make_species('O3', 4, Molecule(smiles='[O-][O+]=O'))
    o_anion = _make_species('O-', 5, Molecule().from_adjacency_list('1 O u1 p3 c-1'))
    return [e, o_atom, o2, o3, o_anion]


def _gas_form_reaction(name, species, T0, n):
    """
    One reaction per Cantera writer case that emits an ``{A, b, Ea}`` mapping
    from a gas-phase RMG rate, built at the given reference temperature.
    """
    e, o_atom, o2, o3, o_anion = species

    def arr(A, units):
        return Arrhenius(A=(A, units), n=n, T0=T0, Ea=_EA)

    if name == 'Arrhenius':
        return Reaction(index=1, reactants=[o_atom, o2], products=[o3], reversible=False,
                        kinetics=arr(1.0e6, 'm^3/(mol*s)'))
    if name == 'ThirdBody':
        return Reaction(index=1, reactants=[o_atom, o_atom], products=[o2], reversible=False,
                        kinetics=ThirdBody(arrheniusLow=arr(1.0e6, 'm^6/(mol^2*s)')))
    if name == 'Troe':
        return Reaction(index=1, reactants=[o_atom, o_atom], products=[o2], reversible=False,
                        kinetics=Troe(arrheniusHigh=arr(1.0e6, 'm^3/(mol*s)'),
                                      arrheniusLow=arr(2.0e6, 'm^6/(mol^2*s)'),
                                      alpha=0.5, T3=(100.0, 'K'), T1=(200.0, 'K'), T2=(300.0, 'K')))
    if name == 'Lindemann':
        return Reaction(index=1, reactants=[o_atom, o_atom], products=[o2], reversible=False,
                        kinetics=Lindemann(arrheniusHigh=arr(1.0e6, 'm^3/(mol*s)'),
                                           arrheniusLow=arr(2.0e6, 'm^6/(mol^2*s)')))
    if name == 'PDepArrhenius':
        return Reaction(index=1, reactants=[o_atom, o2], products=[o3], reversible=False,
                        kinetics=PDepArrhenius(pressures=([0.1, 10.0], 'bar'),
                                               arrhenius=[arr(1.0e6, 'm^3/(mol*s)'),
                                                          arr(4.0e6, 'm^3/(mol*s)')]))
    if name == 'MultiArrhenius':
        return Reaction(index=1, reactants=[o_atom, o2], products=[o3], reversible=False,
                        kinetics=MultiArrhenius(arrhenius=[arr(1.0e6, 'm^3/(mol*s)'),
                                                           arr(3.0e6, 'm^3/(mol*s)')]))
    if name == 'ArrheniusChargeTransfer':
        return Reaction(index=1, reactants=[o_atom], products=[o_anion], electrons=-1,
                        reversible=False,
                        kinetics=ArrheniusChargeTransfer(
                            A=(1.0e6, 'm^3/(mol*s)'), n=n, T0=T0, Ea=_EA,
                            alpha=0.0, electrons=-1, V0=(0.0, 'V')))
    raise AssertionError('no such gas rate form: {0}'.format(name))


GAS_FORMS = ['Arrhenius', 'ThirdBody', 'Troe', 'Lindemann', 'PDepArrhenius', 'MultiArrhenius',
             'ArrheniusChargeTransfer']


@pytest.fixture(scope='module')
def surface_forms_species():
    """Gas and surface species enough to build one of every surface rate form."""
    e = _make_species('e', 1, Molecule().from_adjacency_list('1 e u0 p0 c-1'))
    hplus = _make_species('Hplus', 2, Molecule().from_adjacency_list('1 H u0 p0 c+1'))
    h2 = _make_species('H2', 3, Molecule(smiles='[H][H]'))
    site = _make_species('X', 4, Molecule().from_adjacency_list('1 X u0 p0 c0'))
    hx = _make_species('HX', 5, Molecule().from_adjacency_list(
        '1 H u0 p0 c0 {2,S}\n2 X u0 p0 c0 {1,S}\n'))
    h2x = _make_species('H2X', 6, Molecule().from_adjacency_list(
        '1 H u0 p0 c0 {2,S}\n2 H u0 p0 c0 {1,S}\n3 X u0 p0 c0\n'))
    return [e, hplus, h2, site, hx, h2x]


def _surface_form_reaction(name, species, T0, n):
    """One reaction per surface-phase Cantera writer case that emits an A-factor."""
    e, hplus, h2, site, hx, h2x = species
    if name == 'StickingCoefficient':
        return Reaction(index=1, reactants=[h2, site], products=[h2x], reversible=False,
                        kinetics=StickingCoefficient(A=0.1, n=n, T0=T0, Ea=_EA))
    if name == 'SurfaceArrhenius':
        return Reaction(index=1, reactants=[h2, site], products=[h2x], reversible=False,
                        kinetics=SurfaceArrhenius(A=(1.0e3, 'm^3/(mol*s)'), n=n, T0=T0, Ea=_EA))
    if name == 'SurfaceChargeTransfer':
        return Reaction(index=1, reactants=[hplus, site], products=[hx], electrons=-1,
                        reversible=False,
                        kinetics=SurfaceChargeTransfer(
                            A=(1.0e3, 'm^3/(mol*s)'), n=n, T0=T0, Ea=_EA,
                            alpha=0.0, electrons=-1, V0=(0.0, 'V')))
    raise AssertionError('no such surface rate form: {0}'.format(name))


SURFACE_FORMS = ['StickingCoefficient', 'SurfaceArrhenius', 'SurfaceChargeTransfer']

_RATE_TEMPERATURES = (300.0, 500.0, 800.0, 1500.0, 2500.0)
_RATE_PRESSURES = (0.2 * ct.one_atm, ct.one_atm, 8.0 * ct.one_atm)


def _cantera_gas_from_chemkin(species, reactions, directory):
    """
    Write the Chemkin artifact, convert it with Cantera's own ``ck2yaml`` and
    return the resulting :class:`cantera.Solution`.

    The conversion is what makes this a *rate* comparison rather than a field
    comparison: ck2yaml, not the test, does the ``cm^3/(mol*s)`` to
    ``m^3/(kmol*s)`` conversion, so the test cannot pass by replaying the
    Chemkin writer's own arithmetic.
    """
    from cantera import ck2yaml
    inp_path = os.path.join(directory, 'chem.inp')
    converted = os.path.join(directory, 'from_chemkin.yaml')
    save_chemkin_file(inp_path, species, reactions, check_for_duplicates=False)
    ck2yaml.convert(input_file=inp_path, out_name=converted, permissive=True, quiet=True)
    return ct.Solution(converted)


def _mole_fractions(species):
    return ', '.join('{0}({1}):1'.format(spc.label, spc.index) for spc in species)


class TestWritersAgreeOnTheRate:
    """
    RMG evaluates ``k = A*(T/T0)^n*exp(-Ea/RT)``; Cantera's ``{A, b, Ea}`` mapping
    evaluates ``k = A*T^b*exp(-Ea/RT)`` and has nowhere to put ``T0``. The Chemkin
    writer has always folded ``T0**n`` into the A-factor it writes; the Cantera
    writer did not, so the same reaction exported to the two formats evaluated to
    two different rates -- by a factor of ``T0**n``, silently, in the deliverable
    format.

    These tests compare the *rates* the two artifacts evaluate to. Comparing
    fields or strings would pass while the physics diverges, which is how this
    survived a previous review.
    """

    @pytest.mark.parametrize('form', GAS_FORMS)
    def test_both_writers_evaluate_to_the_same_rate(self, gas_forms_species, form, tmp_out):
        species = gas_forms_species
        reactions = [_gas_form_reaction(form, species, _T0, _N)]
        directory = os.path.join(tmp_out, form)
        os.makedirs(directory)

        yaml_path = os.path.join(directory, 'chem.yaml')
        save_cantera_model(ReactionModel(species=species, reactions=reactions), yaml_path)
        from_cantera = ct.Solution(yaml_path)
        from_chemkin = _cantera_gas_from_chemkin(species, reactions, directory)

        composition = _mole_fractions(species)
        assert from_cantera.n_reactions == from_chemkin.n_reactions
        for temperature in _RATE_TEMPERATURES:
            for pressure in _RATE_PRESSURES:
                from_cantera.TPX = temperature, pressure, composition
                from_chemkin.TPX = temperature, pressure, composition
                for k_cantera, k_chemkin in zip(from_cantera.forward_rate_constants,
                                                from_chemkin.forward_rate_constants):
                    # 1e-3 is the Chemkin *format's* precision -- it carries four
                    # significant figures of A and three decimals of Ea in
                    # kcal/mol. The defect this guards against is a factor of
                    # 9e4, five orders of magnitude clear of that.
                    assert k_cantera == pytest.approx(k_chemkin, rel=1e-3), (
                        '{0} at T={1} K, P={2} Pa: Cantera {3!r} vs Chemkin {4!r}'.format(
                            form, temperature, pressure, k_cantera, k_chemkin))

    @pytest.mark.parametrize('form', GAS_FORMS + SURFACE_FORMS)
    def test_exported_rate_does_not_depend_on_the_reference_temperature(
            self, gas_forms_species, surface_forms_species, form, tmp_out):
        """
        ``change_t0`` rewrites A so the rate law is unchanged. Two kinetics
        objects that differ only in ``T0`` are therefore the same physical rate,
        and must export to the same evaluated rate. This is what covers the
        surface cases, whose Chemkin artifacts ck2yaml cannot read back (it emits
        no definition for the ``X`` site element).
        """
        surface = form in SURFACE_FORMS
        species = surface_forms_species if surface else gas_forms_species
        build = _surface_form_reaction if surface else _gas_form_reaction

        solutions = []
        for label, rebase in (('shifted', False), ('conventional', True)):
            reaction = build(form, species, _T0, _N)
            if rebase:
                _rebase_to_t0_one(reaction.kinetics)
            directory = os.path.join(tmp_out, form + '_' + label)
            os.makedirs(directory)
            path = os.path.join(directory, 'chem.yaml')
            save_cantera_model(ReactionModel(species=species, reactions=[reaction]), path,
                               site_density=_SITE_DENSITY if surface else None)
            solutions.append(_load_for_rates(path, species, surface))

        shifted, conventional = solutions
        for temperature in _RATE_TEMPERATURES:
            for pressure in _RATE_PRESSURES:
                _set_state(shifted, temperature, pressure, species, surface)
                _set_state(conventional, temperature, pressure, species, surface)
                for k_shifted, k_conventional in zip(shifted.forward_rate_constants,
                                                     conventional.forward_rate_constants):
                    assert k_shifted == pytest.approx(k_conventional, rel=1e-12), (
                        '{0} at T={1} K, P={2} Pa: T0=300 K export {3!r} but the same rate '
                        'law written against T0=1 K exports {4!r}'.format(
                            form, temperature, pressure, k_shifted, k_conventional))


def _rebase_to_t0_one(kinetics):
    """
    Rewrite every A-factor inside ``kinetics`` in place so the rate law is
    expressed against ``T0 = 1 K`` without changing the rate it evaluates:
    ``change_t0`` folds the old reference temperature into A.
    """
    inner = getattr(kinetics, 'arrhenius', None)
    if inner is not None:  # MultiArrhenius, PDepArrhenius
        for component in inner:
            _rebase_to_t0_one(component)
        return
    rebased = False
    for attribute in ('arrheniusHigh', 'arrheniusLow'):  # ThirdBody, Lindemann, Troe
        component = getattr(kinetics, attribute, None)
        if component is not None:
            _rebase_to_t0_one(component)
            rebased = True
    if not rebased:
        kinetics.change_t0(1.0)


def _load_for_rates(path, species, surface):
    if not surface:
        return ct.Solution(path)
    import yaml as pyyaml
    with open(path) as handle:
        phase_names = [phase['name'] for phase in pyyaml.safe_load(handle)['phases']]
    return ct.Interface(path, phase_names[-1])


def _set_state(solution, temperature, pressure, species, surface):
    gas_species = [spc for spc in species if not spc.contains_surface_site()]
    if not surface:
        solution.TPX = temperature, pressure, _mole_fractions(species)
        return
    surface_species = [spc for spc in species if spc.contains_surface_site()]
    solution.TP = temperature, pressure
    solution.adjacent['gas'].TPX = temperature, pressure, _mole_fractions(gas_species)
    solution.coverages = ', '.join(
        '{0}({1}):{2}'.format(spc.label, spc.index, 1.0 / len(surface_species))
        for spc in surface_species)


# ---------------------------------------------------------------------------
# 10. The conventional T0 = 1 K export is provably unchanged
# ---------------------------------------------------------------------------

#: Golden Cantera YAML produced by the Cantera writer as it stood at commit
#: 331614fe1, i.e. before ``T0`` normalization was added. RMG's own convention is
#: ``T0 = 1 K``, so ``T0**n == 1`` and the normalization must change nothing;
#: these files are what makes that claim checkable rather than asserted. See
#: test/rmgpy/test_data/plasma_export/README.md for how they were generated.
_GOLDEN_DIR = os.path.join(os.path.dirname(__file__), 'test_data', 'plasma_export')


def _normalize_generator(text):
    """
    Blank out the one entry of a Cantera YAML file that cannot be stable across
    commits: ``generator`` embeds the writer module's path and the current git
    commit, and is long enough that the dumper wraps it over several lines.
    Everything else -- every parameter, every note, every key order -- is
    compared byte for byte.
    """
    lines, out, skipping = text.split('\n'), [], False
    for line in lines:
        if line.startswith('generator:'):
            out.append('generator: <normalized>')
            skipping = True
        elif skipping and (line.startswith(' ') or line.startswith('\t')):
            continue
        else:
            skipping = False
            out.append(line)
    return '\n'.join(out)


def _t0_one_gas_model(species):
    return ReactionModel(
        species=species,
        reactions=[_gas_form_reaction(form, species, (1.0, 'K'), _N) for form in GAS_FORMS],
    )


def _t0_one_surface_model(species):
    return ReactionModel(
        species=species,
        reactions=[_surface_form_reaction(form, species, (1.0, 'K'), _N)
                   for form in SURFACE_FORMS],
    )


class TestConventionalReferenceTemperatureIsUnchanged:
    """
    ``T0`` normalization was applied writer-wide, which is only safe because
    RMG writes ``T0 = 1 K``. These tests are the proof that nothing moved for
    that case: the exported files still match, byte for byte, what the writer
    produced before the change.
    """

    def test_gas_export_is_byte_identical_to_the_base_commit(self, gas_forms_species, tmp_out):
        path = os.path.join(tmp_out, 'chem.yaml')
        save_cantera_model(_t0_one_gas_model(gas_forms_species), path)
        with open(path) as handle:
            written = handle.read()
        with open(os.path.join(_GOLDEN_DIR, 't0_one_gas.yaml')) as handle:
            golden = handle.read()
        assert _normalize_generator(written) == _normalize_generator(golden)

    def test_surface_export_is_byte_identical_to_the_base_commit(self, surface_forms_species,
                                                                 tmp_out):
        path = os.path.join(tmp_out, 'chem.yaml')
        save_cantera_model(_t0_one_surface_model(surface_forms_species), path,
                           site_density=_SITE_DENSITY)
        with open(path) as handle:
            written = handle.read()
        with open(os.path.join(_GOLDEN_DIR, 't0_one_surface.yaml')) as handle:
            golden = handle.read()
        assert _normalize_generator(written) == _normalize_generator(golden)

    def test_the_goldens_really_do_use_the_conventional_reference_temperature(
            self, gas_forms_species, surface_forms_species):
        """
        A golden generated from ``T0 != 1`` kinetics would make the two tests
        above vacuous -- they would compare an un-normalized export to another
        un-normalized export.
        """
        models = [_t0_one_gas_model(gas_forms_species),
                  _t0_one_surface_model(surface_forms_species)]
        checked = 0
        for model in models:
            for reaction in model.reactions:
                for kinetics in _every_rate_object(reaction.kinetics):
                    assert kinetics.T0.value_si == 1.0
                    assert kinetics.n.value_si != 0.0
                    checked += 1
        # 11 gas rate objects (Troe, Lindemann, PDepArrhenius and MultiArrhenius
        # each carry two) and 3 surface ones.
        assert checked == 14


def _every_rate_object(kinetics):
    """Yield every rate object inside ``kinetics`` that carries an A-factor."""
    inner = getattr(kinetics, 'arrhenius', None)
    if inner is not None:
        for component in inner:
            for rate in _every_rate_object(component):
                yield rate
        return
    nested = False
    for attribute in ('arrheniusHigh', 'arrheniusLow'):
        component = getattr(kinetics, attribute, None)
        if component is not None:
            nested = True
            for rate in _every_rate_object(component):
                yield rate
    if not nested:
        yield kinetics


# ---------------------------------------------------------------------------
# 11. Failure behaviour: the artifacts already on disk survive it
# ---------------------------------------------------------------------------

_SENTINEL = 'PREVIOUS GOOD EXPORT - MUST SURVIVE A FAILED ONE\n'


def _leftovers(directory, keep):
    """Files in ``directory`` other than the expected outputs -- i.e. temp files."""
    return sorted(set(os.listdir(directory)) - set(keep))


class TestCanteraWriteIsAtomic:
    """
    The Cantera writer used to hand the destination file straight to
    ``yaml.dump``. A failure inside the dump therefore truncated whatever was
    already there. The previous test only covered a failure raised *before* the
    file was opened, which is a case the old code also survived -- so it read as
    closed while the real one was open.
    """

    def test_a_failure_inside_the_dump_leaves_the_previous_file_intact(
            self, mechanism, tmp_out, monkeypatch):
        import rmgpy.yaml_cantera2 as yaml_cantera2

        model, species, reactions = mechanism
        path = os.path.join(tmp_out, 'chem_annotated.yaml')
        with open(path, 'w') as handle:
            handle.write(_SENTINEL)

        def exploding_dump(data, stream, **kwargs):
            # Write real content first, then fail: this is what a mid-write
            # failure looks like, and what truncates a streamed destination.
            stream.write('phases:\n- name: gas\n  thermo: ideal-gas\n')
            raise RuntimeError('dump failed half way through')

        monkeypatch.setattr(yaml_cantera2.yaml, 'dump', exploding_dump)
        with pytest.raises(RuntimeError, match='half way through'):
            save_cantera_model(model, path)

        with open(path) as handle:
            assert handle.read() == _SENTINEL
        assert _leftovers(tmp_out, ['chem_annotated.yaml']) == []

    def test_a_failure_inside_the_dump_creates_no_file_where_there_was_none(
            self, mechanism, tmp_out, monkeypatch):
        import rmgpy.yaml_cantera2 as yaml_cantera2

        model, species, reactions = mechanism
        path = os.path.join(tmp_out, 'chem_annotated.yaml')

        def exploding_dump(data, stream, **kwargs):
            stream.write('phases:\n')
            raise RuntimeError('dump failed half way through')

        monkeypatch.setattr(yaml_cantera2.yaml, 'dump', exploding_dump)
        with pytest.raises(RuntimeError):
            save_cantera_model(model, path)

        assert os.listdir(tmp_out) == []


class TestChemkinOutputSetIsCoherent:
    """
    Each Chemkin file is written atomically, but ``save_chemkin`` writes several.
    Landing them one at a time let a failure on a later file leave a new gas file
    beside a stale surface file: a mechanism split across two generations of the
    model, which reads as a valid export and is not one.
    """

    @pytest.fixture()
    def surface_model(self, surface_forms_species):
        species = surface_forms_species
        e, hplus, h2, site, hx, h2x = species
        good = Reaction(index=1, reactants=[h2, site], products=[h2x], reversible=False,
                        kinetics=SurfaceArrhenius(A=(1.0e3, 'm^3/(mol*s)'), n=0.0,
                                                  Ea=(10.0, 'kJ/mol')))
        gas_only = Reaction(index=2, reactants=[h2], products=[h2], reversible=False,
                            kinetics=Arrhenius(A=(1.0e3, 's^-1'), n=0.0, Ea=(10.0, 'kJ/mol')))

        class _Core(object):
            pass

        core = _Core()
        core.species = species
        core.reactions = [gas_only, good]
        model = _Core()
        model.core = core
        model.output_species_list = []
        model.output_reaction_list = []
        model.surface_site_density = None
        return model, species

    def _paths(self, tmp_out):
        return dict(
            path=os.path.join(tmp_out, 'chem.inp'),
            verbose_path=os.path.join(tmp_out, 'chem_annotated.inp'),
            dictionary_path=os.path.join(tmp_out, 'species_dictionary.txt'),
            transport_path=os.path.join(tmp_out, 'tran.dat'),
        )

    def _expected_files(self):
        return ['chem-gas.inp', 'chem-surface.inp', 'chem_annotated-gas.inp',
                'chem_annotated-surface.inp', 'species_dictionary.txt', 'tran.dat']

    def test_a_full_export_writes_the_whole_set(self, surface_model, tmp_out):
        model, species = surface_model
        save_chemkin(model, **self._paths(tmp_out))
        assert sorted(os.listdir(tmp_out)) == sorted(self._expected_files())

    def test_a_failure_on_a_later_file_leaves_the_earlier_ones_untouched(
            self, surface_model, tmp_out):
        model, species = surface_model
        paths = self._paths(tmp_out)

        # A first, good export, so there is a real previous generation on disk.
        save_chemkin(model, **paths)
        before = {}
        for name in os.listdir(tmp_out):
            with open(os.path.join(tmp_out, name)) as handle:
                before[name] = handle.read()

        # Move the gas rate *and* break the surface reaction. The gas file is
        # rendered first, so a writer that landed as it went would already have
        # replaced it with the new rate by the time the surface file failed --
        # leaving a new gas file beside a stale surface file. Moving the gas rate
        # is what makes that visible: without it the rewritten gas file would
        # have identical bytes and the check would pass either way.
        model.core.reactions[0].kinetics = Arrhenius(A=(5.0e3, 's^-1'), n=0.0,
                                                     Ea=(10.0, 'kJ/mol'))
        model.core.reactions[1].kinetics = UnhandledKinetics()
        with pytest.raises(MechanismWriterError):
            save_chemkin(model, **paths)

        after = {}
        for name in os.listdir(tmp_out):
            with open(os.path.join(tmp_out, name)) as handle:
                after[name] = handle.read()
        assert after == before
        assert _leftovers(tmp_out, self._expected_files()) == []

    def test_a_failure_on_a_later_file_writes_nothing_at_all_the_first_time(
            self, surface_model, tmp_out):
        model, species = surface_model
        model.core.reactions[1].kinetics = UnhandledKinetics()
        with pytest.raises(MechanismWriterError):
            save_chemkin(model, **self._paths(tmp_out))
        assert os.listdir(tmp_out) == []


class TestMultiKineticsReconstructionKeepsTheReaction:
    """
    Each component of a grouped ``MultiArrhenius``/``MultiPDepArrhenius`` is
    re-entered into the writer as a freshly built ``Reaction``. Anything not
    named in that constructor call is gone by the time the entry is built --
    which is how the third-body collider and the flux pairs were being dropped
    from grouped reactions only.
    """

    @pytest.fixture()
    def grouped_reaction(self, gas_forms_species):
        e, o_atom, o2, o3, o_anion = gas_forms_species
        kinetics = MultiArrhenius(arrhenius=[
            Arrhenius(A=(1.0e6, 'm^3/(mol*s)'), n=0.0, Ea=(10.0, 'kJ/mol')),
            Arrhenius(A=(3.0e6, 'm^3/(mol*s)'), n=0.0, Ea=(20.0, 'kJ/mol')),
        ])
        return Reaction(index=17, reactants=[o_atom], products=[o_anion], electrons=-1,
                        reversible=False, specific_collider=o2,
                        pairs=[(o_atom, o_anion)], kinetics=kinetics)

    def test_every_component_keeps_index_collider_pairs_and_electrons(
            self, gas_forms_species, grouped_reaction, monkeypatch):
        import rmgpy.yaml_cantera2 as yaml_cantera2

        built = []
        real_reaction = yaml_cantera2.Reaction

        class RecordingReaction(real_reaction):
            def __init__(self, *args, **kwargs):
                super(RecordingReaction, self).__init__(*args, **kwargs)
                built.append(self)

        monkeypatch.setattr(yaml_cantera2, 'Reaction', RecordingReaction)
        entries = yaml_cantera2.reaction_to_dict_list(grouped_reaction, gas_forms_species)

        assert len(built) == 2, 'expected one sub-reaction per grouped component'
        for sub in built:
            assert sub.index == 17
            assert sub.specific_collider is grouped_reaction.specific_collider
            assert sub.pairs == grouped_reaction.pairs
            assert sub.electrons == -1

        # ...and the artifact carries what the artifact can carry.
        assert len(entries) == 2
        for entry in entries:
            assert 'Specific third body collider: O2' in entry['note']
            assert 'Flux pairs: O(2), O-(5)' in entry['note']
            assert 'e(1)' in entry['equation']

    def test_a_grouped_third_body_reaction_keeps_its_collider_in_the_equation(
            self, gas_forms_species):
        """
        ``MultiPDepArrhenius`` groups ``PDepArrhenius`` components, but the same
        reconstruction feeds ``ThirdBody``-shaped equations, where the collider
        is part of the equation string rather than the note. Written as a direct
        check on ``get_reaction_equation`` so the coupling is explicit.
        """
        import rmgpy.yaml_cantera2 as yaml_cantera2

        e, o_atom, o2, o3, o_anion = gas_forms_species
        reaction = Reaction(index=1, reactants=[o_atom, o_atom], products=[o2], reversible=False,
                            specific_collider=o3,
                            kinetics=ThirdBody(arrheniusLow=Arrhenius(
                                A=(1.0e6, 'm^6/(mol^2*s)'), n=0.0, Ea=(10.0, 'kJ/mol'))))
        equation = yaml_cantera2.get_reaction_equation(reaction, gas_forms_species)
        assert equation.count('O3(4)') == 2
