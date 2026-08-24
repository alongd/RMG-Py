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

"""
Routing of ``Marcus`` kinetics through :meth:`Reaction.get_rate_coefficient`.

``Reaction.get_rate_coefficient(T, P, ...)`` used to end in a bare
``self.kinetics.get_rate_coefficient(T, P)``, which assumes every kinetics class reads its second
positional argument as a pressure in Pa. ``Marcus.get_rate_coefficient(double T, double dGrxn=0.0)``
does not: it reads a free energy of reaction in J/mol. A 1 bar reactor therefore evaluated

    dG_act = (lmbd_i(T) + lmbd_o)/4 * (1 + 1e5/(lmbd_i(T) + lmbd_o))^2

and with no pressure supplied at all -- the far commoner case -- ``P`` defaulted to 0, so every
Marcus reaction was evaluated as *exactly thermoneutral*, its barrier collapsed to the intrinsic
``(lmbd_i + lmbd_o)/4`` no matter what its real thermochemistry said. That failure is silent: no
exception, no warning, a physically plausible number. ``Cation_R_Recombination`` is 12/12 Marcus,
so it is the whole family that was affected.

These tests pin both halves of the fix, because either half alone is satisfiable by a wrong fix:

* no pressure ever reaches the rate law -- but a fix hardcoding ``dGrxn = 0.0`` also achieves that;
* the rate answers to the reaction's own thermochemistry, along the Marcus parabola specifically --
  which a linear Evans-Polanyi-shaped substitution would not reproduce.
"""

import os

import numpy as np
import pytest

from rmgpy import constants, settings
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.exceptions import KineticsError
from rmgpy.kinetics.arrhenius import Marcus
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.thermo import ThermoData

#: Pressures spanning four decades. ``Marcus`` carries no pressure dependence whatsoever, so a test
#: at a single pressure would also pass on a fix that merely substituted a different wrong number.
PRESSURES = [1.0e3, 1.0e4, 1.0e5, 1.0e6]

#: Total reorganization energy of the fixtures, in J/mol. Held in one place because the fixtures
#: below are positioned relative to it: the Marcus parabola is symmetric about ``dGrxn = -LAMBDA``.
LAMBDA = 1.0e5

T = 500.0


def make_marcus_kinetics():
    """
    A temperature-independent ``Marcus`` model with total reorganization energy ``LAMBDA``.

    ``lmbd_i_coefs`` is a cubic in T; only the constant term is populated so that ``lmbd_i`` does
    not move with temperature and the fixtures' arithmetic stays checkable by hand.
    """
    return Marcus(
        A=(1.0e10, "m^3/(mol*s)"), n=0.0,
        lmbd_i_coefs=np.array([LAMBDA, 0.0, 0.0, 0.0]),
        wr=(0, "J/mol"), wp=(0, "J/mol"), lmbd_o=(0.0, "J/mol"),
    )


def _species(smiles, h298_kj):
    """A species whose only distinguishing thermochemistry is its enthalpy of formation."""
    return Species(
        molecule=[Molecule().from_smiles(smiles)],
        thermo=ThermoData(
            Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
            Cpdata=([8.0] * 7, "cal/(mol*K)"),
            H298=(h298_kj, "kJ/mol"), S298=(50.0, "cal/(mol*K)"),
        ),
    )


def make_reaction(dg_kj, electrons=0):
    """
    A one-to-one gas-phase reaction with free energy of reaction ``dg_kj`` kJ/mol at every T.

    Reactant and product carry identical Cp and S298, so the entropy and heat-capacity terms cancel
    exactly and ``get_free_energy_of_reaction(T)`` returns the enthalpy difference at any T. That
    makes the fixture's driving force an input to the test rather than something it has to measure.
    """
    return Reaction(
        reactants=[_species("[CH3]", 0.0)],
        products=[_species("[NH2]", dg_kj)],
        electrons=electrons,
        kinetics=make_marcus_kinetics(),
    )


def test_fixture_free_energy_is_the_number_it_claims_to_be():
    """Guard on the fixture itself: every assertion below is stated in terms of these values."""
    for dg_kj in (-200.0, -40.0, 0.0, 40.0):
        assert make_reaction(dg_kj).get_free_energy_of_reaction(T) == pytest.approx(
            dg_kj * 1000.0, rel=1e-9, abs=1e-6
        )


class TestPressureDoesNotReachTheRateLaw:
    """Half one: the pressure argument is not a free energy and must never be read as one."""

    def test_pressure_does_not_reach_the_rate_law(self):
        """
        The regression itself. At 1 bar the old path evaluated ``dGrxn = 1e5 J/mol``; at 10 bar it
        underflowed to a denormal. Neither number has anything to do with the reaction.
        """
        rxn = make_reaction(0.0)
        fed_the_pressure = rxn.kinetics.get_rate_coefficient(T, 1.0e5)

        assert rxn.get_rate_coefficient(T, 1.0e5) != fed_the_pressure

    def test_rate_is_invariant_to_pressure(self):
        rxn = make_reaction(40.0)
        rates = [rxn.get_rate_coefficient(T, P) for P in PRESSURES]

        assert rates[0] > 0.0, "invariant but annihilated is not a fix"
        assert len(set(rates)) == 1, f"rate varied with pressure: {list(zip(PRESSURES, rates))}"

    def test_omitted_pressure_and_supplied_pressure_agree(self):
        """
        ``P`` defaults to 0, which used to mean "thermoneutral" rather than "no pressure given".
        Leaving it out and passing 1 bar must now be the same call.
        """
        rxn = make_reaction(40.0)
        assert rxn.get_rate_coefficient(T) == rxn.get_rate_coefficient(T, 1.0e5)


class TestRateRespondsToTheReactionsOwnThermo:
    """
    Half two: the tripwire. Every test in the class above passes on a fix that hardcodes
    ``dGrxn = 0.0`` -- which would be the same class of error as passing the pressure, just quieter.
    """

    def test_free_energy_of_reaction_reaches_the_rate_law(self):
        rxn = make_reaction(40.0)
        dGrxn = rxn.get_free_energy_of_reaction(T)

        assert dGrxn != 0.0, "fixture is thermoneutral, so it cannot detect a laundered dGrxn"
        assert rxn.get_rate_coefficient(T, 1.0e5) != rxn.kinetics.get_rate_coefficient(T, 0.0)
        assert rxn.get_rate_coefficient(T, 1.0e5) == rxn.kinetics.get_rate_coefficient(T, dGrxn)

    def test_endothermic_is_slower_and_exothermic_is_faster(self):
        """
        In the normal Marcus region (|dGrxn| < lmbd) the barrier rises monotonically with dGrxn,
        so k must fall monotonically. Two reactions differing only in the product's free energy.
        """
        exothermic = make_reaction(-40.0).get_rate_coefficient(T, 1.0e5)
        thermoneutral = make_reaction(0.0).get_rate_coefficient(T, 1.0e5)
        endothermic = make_reaction(40.0).get_rate_coefficient(T, 1.0e5)

        assert exothermic > thermoneutral > endothermic

    def test_barriers_match_the_marcus_parabola_by_hand(self):
        """
        dG_act = lmbd/4 * (1 + dGrxn/lmbd)^2, computed here from the fixture's own numbers rather
        than from the object under test, so an error inside ``get_gibbs_activation_energy`` cannot
        cancel out of both sides.
        """
        for dg_kj in (-40.0, 0.0, 40.0):
            dGrxn = dg_kj * 1000.0
            expected_barrier = LAMBDA / 4.0 * (1.0 + dGrxn / LAMBDA) ** 2
            expected_k = 1.0e10 * np.exp(-expected_barrier / (constants.R * T))

            assert make_reaction(dg_kj).get_rate_coefficient(T, 1.0e5) == pytest.approx(
                expected_k, rel=1e-12
            )

    def test_inverted_region_is_reproduced(self):
        """
        The signature that separates Marcus from any linear barrier-vs-driving-force model: the
        parabola is symmetric about ``dGrxn = -lmbd``, so a reaction at ``dGrxn = -2*lmbd`` has
        *exactly* the barrier of a thermoneutral one, despite being 200 kJ/mol downhill.

        A fix that hardcoded ``dGrxn = 0`` would pass this test and fail the two above; a fix that
        substituted the driving force into a linear (Evans-Polanyi) form would pass those two and
        fail this one. Only the real dGrxn on the real rate law passes all three.
        """
        inverted = make_reaction(-2.0 * LAMBDA / 1000.0)
        thermoneutral = make_reaction(0.0)

        assert inverted.get_free_energy_of_reaction(T) == pytest.approx(-2.0 * LAMBDA, rel=1e-9)
        assert inverted.get_rate_coefficient(T, 1.0e5) == pytest.approx(
            thermoneutral.get_rate_coefficient(T, 1.0e5), rel=1e-12
        )
        # ...and it is genuinely far downhill, not a degenerate fixture.
        assert inverted.get_rate_coefficient(T, 1.0e5) < make_reaction(-40.0).get_rate_coefficient(T, 1.0e5)


class TestPotentialHandling:
    """
    ``Marcus`` has no ``V0``: its constructor takes none, no property exposes one, and its rate law
    has no Butler-Volmer term. What an applied field does reach is the *thermodynamics* -- dGrxn of
    a charge transfer reaction shifts by ``-n*F*V`` -- and that is what this branch passes through.
    """

    def test_no_electrode_means_no_field(self):
        """The gas-phase default: ``potential=None`` is "no electrode", i.e. exactly 0 V."""
        rxn = make_reaction(40.0, electrons=-1)

        assert rxn.get_rate_coefficient(T, 1.0e5) == \
            rxn.kinetics.get_rate_coefficient(T, rxn.get_free_energy_of_reaction(T, potential=0.0))

    def test_applied_potential_shifts_the_driving_force(self):
        """An explicitly requested field is honoured, through dGrxn and only through dGrxn."""
        rxn = make_reaction(40.0, electrons=-1)
        shifted = rxn.get_free_energy_of_reaction(T, potential=0.1)

        assert shifted == pytest.approx(40.0e3 + constants.F * 0.1, rel=1e-12)
        assert rxn.get_rate_coefficient(T, 1.0e5, potential=0.1) == \
            rxn.kinetics.get_rate_coefficient(T, shifted)
        # electrons = -1 raises dGrxn with V, and dGrxn = +40 kJ/mol is in the normal region.
        assert rxn.get_rate_coefficient(T, 1.0e5, potential=0.1) < rxn.get_rate_coefficient(T, 1.0e5)

    def test_potential_is_ignored_when_no_electrons_are_transferred(self):
        """A neutral Marcus reaction has no ``-n*F*V`` term, so a field cannot move its rate."""
        rxn = make_reaction(40.0)

        assert rxn.get_rate_coefficient(T, 1.0e5, potential=0.5) == rxn.get_rate_coefficient(T, 1.0e5)


class TestMissingThermoRaises:
    """Falling back on a thermoneutral 0.0 is the bug, not the fallback."""

    def test_missing_thermo_raises_instead_of_inventing_a_free_energy(self):
        rxn = Reaction(
            reactants=[Species(molecule=[Molecule().from_smiles("[CH3]")])],
            products=[Species(molecule=[Molecule().from_smiles("[NH2]")])],
            kinetics=make_marcus_kinetics(),
        )

        with pytest.raises(KineticsError) as exc:
            rxn.get_rate_coefficient(T, 1.0e5)
        assert "Marcus" in str(exc.value)


class TestMarcusCarriesNoReferencePotential:
    """The ``V0`` row in ``Marcus``'s class docstring described a parameter that does not exist."""

    def test_marcus_has_no_v0(self):
        kinetics = make_marcus_kinetics()

        assert not hasattr(kinetics, "V0")
        assert "V0" not in repr(kinetics)
        assert "`V0`" not in Marcus.__doc__


@pytest.mark.database
class TestRealDatabaseMarcusFamily:
    """
    The same guarantee, on the family that actually ships Marcus kinetics.

    ``Cation_R_Recombination`` has no non-Marcus path through it, so a defect on this branch is a
    defect on every reaction the family generates.
    """

    FAMILY = "Cation_R_Recombination"

    @classmethod
    def setup_class(cls):
        database = KineticsDatabase()
        database.load_families(
            path=os.path.join(settings["database.directory"], "kinetics", "families"),
            families=[cls.FAMILY],
        )
        family = database.families[cls.FAMILY]

        cls.rules = [entry.data
                     for entries in family.rules.entries.values()
                     for entry in entries]
        cls.training = [entry.data for entry in family.get_training_depository().entries.values()]

    def test_the_family_really_is_all_marcus(self):
        """Guard on the fixture: a non-Marcus entry would make the tests below vacuously pass."""
        assert self.rules, "no rate rules loaded"
        assert self.training, "no training reactions loaded"
        assert all(isinstance(k, Marcus) for k in self.rules)
        assert all(isinstance(k, Marcus) for k in self.training)

    def test_real_reorganization_energies_are_small(self):
        """
        Recorded because it sets the driving forces the test below may use: this family's total
        reorganization energy is only ~16-54 kJ/mol, so a routine +-40 kJ/mol reaction can sit past
        ``dGrxn = -lmbd`` and into the Marcus inverted region, where k falls again as the reaction
        gets *more* exothermic. A monotone "more exothermic is faster" assertion is therefore not a
        safe way to test these entries -- it fails on correct physics.
        """
        for kinetics in self.rules + self.training:
            lmbd = kinetics.get_lmbd_i(T) + kinetics.lmbd_o.value_si
            assert 1.0e4 < lmbd < 6.0e4, f"{kinetics!r} has lmbd = {lmbd} J/mol at {T} K"

    def test_real_rate_rules_ignore_pressure_and_answer_to_thermo(self):
        for kinetics in self.rules + self.training:
            # Driving forces scaled to each entry's own lmbd, so every fixture stays in the normal
            # region (|dGrxn| < lmbd) and the ordering below is the physics, not an accident.
            step_kj = (kinetics.get_lmbd_i(T) + kinetics.lmbd_o.value_si) / 4.0 / 1000.0
            rates = {}
            for dg_kj in (-step_kj, 0.0, step_kj):
                rxn = make_reaction(dg_kj)
                rxn.kinetics = kinetics
                at_each_pressure = [rxn.get_rate_coefficient(T, P) for P in PRESSURES]

                assert at_each_pressure[0] > 0.0, f"{kinetics!r} annihilated at dGrxn={dg_kj}"
                assert len(set(at_each_pressure)) == 1, (
                    f"{kinetics!r} varied with pressure: "
                    f"{list(zip(PRESSURES, at_each_pressure))}"
                )
                assert at_each_pressure[0] == \
                    kinetics.get_rate_coefficient(T, rxn.get_free_energy_of_reaction(T))
                rates[dg_kj] = at_each_pressure[0]

            assert rates[-step_kj] > rates[0.0] > rates[step_kj], f"{kinetics!r} ignored its dGrxn"
