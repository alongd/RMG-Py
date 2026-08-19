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
Routing of gas-phase charge transfer kinetics through :meth:`Reaction.get_rate_coefficient`.

``Reaction.get_rate_coefficient(T, P, ...)`` used to end in a bare
``self.kinetics.get_rate_coefficient(T, P)``, which assumes every kinetics class reads its second
positional argument as a pressure in Pa. The two gas-phase charge transfer classes do not:
``ArrheniusChargeTransfer`` reads it as an electrode potential in V, and
``ArrheniusChargeTransferBM`` reads it as an enthalpy of reaction in J/mol. Handing them a
pressure of 1e5 shifted the barrier by ~4.8e9 J/mol and returned exactly 0.0, silently deleting
these reactions from every gas-phase simulation.

These tests pin the fix: no pressure ever reaches those two classes, the default evaluation is at
the reference potential ``V0``, an explicitly requested potential is honoured, and the surface path
is left exactly where it was.
"""

import os

import pytest

from rmgpy import settings
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.exceptions import KineticsError
from rmgpy.kinetics.arrhenius import ArrheniusChargeTransfer, ArrheniusChargeTransferBM
from rmgpy.kinetics.diffusionLimited import diffusion_limiter
from rmgpy.kinetics.surface import StickingCoefficient, SurfaceChargeTransfer
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.thermo import ThermoData

#: Pressures spanning four decades. These classes carry no pressure dependence at all, so a test
#: at a single pressure would also pass on a fix that merely substituted a different wrong number.
PRESSURES = [1.0e3, 1.0e4, 1.0e5, 1.0e6]


def make_charge_transfer_kinetics(V0=0.0):
    """A gas-phase ArrheniusChargeTransfer with volumetric units, as the cation families use."""
    return ArrheniusChargeTransfer(
        A=(1.0e10, "m^3/(mol*s)"), n=0, Ea=(20.0, "kJ/mol"), V0=(V0, "V"),
        alpha=0.5, electrons=-1, T0=(1.0, "K"),
    )


def make_charge_transfer_bm_kinetics(V0=0.0):
    """The Blowers-Masel form of the same, as the cation families' rate rules use."""
    return ArrheniusChargeTransferBM(
        A=(1.0e10, "m^3/(mol*s)"), n=0, w0=(100.0, "kJ/mol"), E0=(20.0, "kJ/mol"),
        V0=(V0, "V"), alpha=0.5, electrons=-1,
    )


class TestArrheniusChargeTransferRouting:
    """``ArrheniusChargeTransfer`` is evaluated at V0, never at the caller's pressure."""

    def test_pressure_does_not_reach_the_rate_law(self):
        """The regression itself: a 1 bar reactor used to get exactly 0.0 out of this."""
        kinetics = make_charge_transfer_kinetics()
        rxn = Reaction(kinetics=kinetics)
        expected = kinetics.get_rate_coefficient(1000.0, kinetics.V0.value_si)

        assert expected > 0.0
        assert rxn.get_rate_coefficient(1000.0, 1.0e5) == expected

    def test_evaluated_at_v0_not_at_zero_volts(self):
        """A non-zero V0 shows the default is really V0 and not 0 V by coincidence."""
        kinetics = make_charge_transfer_kinetics(V0=0.5)
        rxn = Reaction(kinetics=kinetics)

        at_v0 = kinetics.get_rate_coefficient(1000.0, 0.5)
        at_zero_volts = kinetics.get_rate_coefficient(1000.0, 0.0)
        assert at_v0 != at_zero_volts, "fixture is degenerate: V0 and 0 V give the same rate"

        assert rxn.get_rate_coefficient(1000.0, 1.0e5) == at_v0

    def test_rate_is_invariant_to_pressure(self):
        """No pressure dependence exists in this model, so none may appear in the answer."""
        rxn = Reaction(kinetics=make_charge_transfer_kinetics(V0=0.5))
        rates = [rxn.get_rate_coefficient(1000.0, P) for P in PRESSURES]

        # Without this, the old behaviour -- 0.0 at every pressure -- would satisfy the test.
        assert rates[0] > 0.0, "invariant but annihilated is not a fix"
        assert len(set(rates)) == 1, f"rate varied with pressure: {list(zip(PRESSURES, rates))}"

    def test_explicit_potential_is_honoured(self):
        """A deliberately applied field is legitimate and must still be obeyed."""
        kinetics = make_charge_transfer_kinetics()
        rxn = Reaction(kinetics=kinetics)

        assert rxn.get_rate_coefficient(1000.0, 1.0e5, potential=0.1) == \
            kinetics.get_rate_coefficient(1000.0, 0.1)

    def test_explicit_potential_moves_the_rate_the_right_way(self):
        """Ea -= alpha * electrons * F * (V - V0); with electrons = -1, raising V raises Ea."""
        kinetics = make_charge_transfer_kinetics()
        rxn = Reaction(kinetics=kinetics)

        at_v0 = rxn.get_rate_coefficient(1000.0, 1.0e5)
        above = rxn.get_rate_coefficient(1000.0, 1.0e5, potential=0.1)
        below = rxn.get_rate_coefficient(1000.0, 1.0e5, potential=-0.1)

        assert below > at_v0 > above

    def test_explicit_zero_volts_is_distinguishable_from_no_potential(self):
        """``potential=0.0`` means an applied 0 V, which is not the same as leaving it unset."""
        kinetics = make_charge_transfer_kinetics(V0=0.5)
        rxn = Reaction(kinetics=kinetics)

        assert rxn.get_rate_coefficient(1000.0, 1.0e5, potential=0.0) == \
            kinetics.get_rate_coefficient(1000.0, 0.0)
        assert rxn.get_rate_coefficient(1000.0, 1.0e5) == \
            kinetics.get_rate_coefficient(1000.0, 0.5)


def make_endothermic_reaction(kinetics):
    """
    An endothermic gas-phase reaction carrying `kinetics`.

    The Blowers-Masel barrier depends on the enthalpy of reaction, so a thermoneutral fixture would
    not distinguish "resolved against the reaction's own enthalpy" from "substituted 0.0".
    """
    def species(smiles, h298):
        return Species(
            molecule=[Molecule().from_smiles(smiles)],
            thermo=ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                Cpdata=([8.0] * 7, "cal/(mol*K)"),
                H298=(h298, "kJ/mol"), S298=(50.0, "cal/(mol*K)"),
            ),
        )

    return Reaction(
        reactants=[species("[CH3]", 0.0)], products=[species("[NH2]", 40.0)], kinetics=kinetics,
    )


class TestArrheniusChargeTransferBMRouting:
    """``ArrheniusChargeTransferBM`` reads its second argument as dHrxn, so it too must be spared."""

    def test_pressure_does_not_reach_the_rate_law(self):
        """1e5 Pa read as 1e5 J/mol is a ~1.5e4x error -- wrong, but not zero, so easy to miss."""
        rxn = make_endothermic_reaction(make_charge_transfer_bm_kinetics())

        fed_the_pressure = rxn.kinetics.get_rate_coefficient(1000.0, 1.0e5)
        assert rxn.get_rate_coefficient(1000.0, 1.0e5) != fed_the_pressure

    def test_barrier_is_resolved_against_the_reactions_own_enthalpy(self):
        """
        The tripwire against trading one substituted quantity for another: a fix that passed
        ``dHrxn = 0.0`` instead of the pressure would be just as wrong, and would pass every other
        test in this class.
        """
        rxn = make_endothermic_reaction(make_charge_transfer_bm_kinetics())
        h298 = rxn.get_enthalpy_of_reaction(298.0)

        assert h298 > 0.0, "fixture is thermoneutral, so it cannot detect a laundered dHrxn"
        assert rxn.get_rate_coefficient(1000.0, 1.0e5) != rxn.kinetics.get_rate_coefficient(1000.0, 0.0)
        # Not bit-identical: to_arrhenius_charge_transfer round-trips Ea through kJ/mol.
        assert rxn.get_rate_coefficient(1000.0, 1.0e5) == \
            pytest.approx(rxn.kinetics.get_rate_coefficient(1000.0, h298), rel=1e-12)

    def test_agrees_with_the_conversion_fix_barrier_height_would_have_done(self):
        """
        ``fix_barrier_height`` converts this model to an ``ArrheniusChargeTransfer`` using H298.
        Routing must not disagree with the pipeline's own conversion, or the rate would change
        depending on whether that conversion had happened yet.
        """
        rxn = make_endothermic_reaction(make_charge_transfer_bm_kinetics())
        routed = rxn.get_rate_coefficient(1000.0, 1.0e5)

        converted = rxn.copy()
        converted.kinetics = rxn.kinetics.to_arrhenius_charge_transfer(rxn.get_enthalpy_of_reaction(298.0))
        assert isinstance(converted.kinetics, ArrheniusChargeTransfer)
        assert converted.get_rate_coefficient(1000.0, 1.0e5) == routed

    def test_rate_is_invariant_to_pressure(self):
        rxn = make_endothermic_reaction(make_charge_transfer_bm_kinetics())
        rates = [rxn.get_rate_coefficient(1000.0, P) for P in PRESSURES]

        assert rates[0] > 0.0, "invariant but annihilated is not a fix"
        assert len(set(rates)) == 1, f"rate varied with pressure: {list(zip(PRESSURES, rates))}"

    def test_explicit_potential_is_honoured(self):
        rxn = make_endothermic_reaction(make_charge_transfer_bm_kinetics())

        at_v0 = rxn.get_rate_coefficient(1000.0, 1.0e5)
        above = rxn.get_rate_coefficient(1000.0, 1.0e5, potential=0.1)
        below = rxn.get_rate_coefficient(1000.0, 1.0e5, potential=-0.1)

        assert below > at_v0 > above

    def test_applied_potential_uses_the_records_own_alpha(self):
        """
        ``to_arrhenius_charge_transfer`` silently resets alpha to its 0.5 default, so routing
        through it has to put the record's alpha back. alpha cancels at V0 but not under a field.
        """
        rxn = make_endothermic_reaction(make_charge_transfer_bm_kinetics())
        rxn.kinetics.alpha = 0.9
        default_alpha = make_endothermic_reaction(make_charge_transfer_bm_kinetics())

        assert default_alpha.kinetics.alpha.value_si == 0.5
        assert rxn.get_rate_coefficient(1000.0, 1.0e5) == \
            default_alpha.get_rate_coefficient(1000.0, 1.0e5), "alpha must cancel at V0"
        assert rxn.get_rate_coefficient(1000.0, 1.0e5, potential=0.1) != \
            default_alpha.get_rate_coefficient(1000.0, 1.0e5, potential=0.1)

        # Any non-default alpha would satisfy the inequality above, so pin the actual value: the
        # rate must equal that of the equivalent ArrheniusChargeTransfer carrying alpha = 0.9.
        equivalent = rxn.kinetics.to_arrhenius_charge_transfer(rxn.get_enthalpy_of_reaction(298.0))
        equivalent.alpha = 0.9
        assert rxn.get_rate_coefficient(1000.0, 1.0e5, potential=0.1) == \
            equivalent.get_rate_coefficient(1000.0, 0.1)
        # ...and that value is not what the lossy converter's own 0.5 default would give.
        assert equivalent.get_rate_coefficient(1000.0, 0.1) != \
            rxn.kinetics.to_arrhenius_charge_transfer(
                rxn.get_enthalpy_of_reaction(298.0)).get_rate_coefficient(1000.0, 0.1)

    def test_missing_thermo_raises_instead_of_inventing_an_enthalpy(self):
        """
        The Blowers-Masel barrier is undefined without an enthalpy of reaction. Returning a number
        anyway -- from a thermoneutral 0.0, say -- is the failure mode this whole module is about,
        so the routing raises with the cause named.

        This is stricter than the base commit, which happily returned a rate computed from the
        pressure. A species-less Reaction is covered too: get_enthalpy_of_reaction raises TypeError
        on it rather than summing to zero.
        """
        thermoless = Reaction(
            reactants=[Species(molecule=[Molecule().from_smiles("[CH3]")])],
            products=[Species(molecule=[Molecule().from_smiles("[NH2]")])],
            kinetics=make_charge_transfer_bm_kinetics(),
        )
        with pytest.raises(KineticsError, match="enthalpy of reaction"):
            thermoless.get_rate_coefficient(1000.0, 1.0e5)

        speciesless = Reaction(kinetics=make_charge_transfer_bm_kinetics())
        with pytest.raises(KineticsError, match="enthalpy of reaction"):
            speciesless.get_rate_coefficient(1000.0, 1.0e5)

    def test_potential_equal_to_v0_matches_the_default(self):
        rxn = make_endothermic_reaction(make_charge_transfer_bm_kinetics(V0=0.3))

        assert rxn.get_rate_coefficient(1000.0, 1.0e5, potential=0.3) == \
            rxn.get_rate_coefficient(1000.0, 1.0e5)


class TestSurfacePathIsUnchanged:
    """The surface charge transfer path, sticking coefficients and the diffusion limiter."""

    @classmethod
    def setup_class(cls):
        m_proton = Molecule().from_smiles("[H+]")
        m_ch2x = Molecule().from_adjacency_list(
            """
            1 C u0 p0 c0 {2,S} {3,S} {4,D}
            2 H u0 p0 c0 {1,S}
            3 H u0 p0 c0 {1,S}
            4 X u0 p0 c0 {1,D}
            """
        )
        m_ch3x = Molecule().from_adjacency_list(
            """
            1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
            2 H u0 p0 c0 {1,S}
            3 H u0 p0 c0 {1,S}
            4 H u0 p0 c0 {1,S}
            5 X u0 p0 c0 {1,S}
            """
        )

        s_proton = Species(
            molecule=[m_proton],
            thermo=ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                Cpdata=([3.4475, 3.4875, 3.497, 3.5045, 3.5405, 3.6095, 3.86], "cal/(mol*K)"),
                H298=(0, "kcal/mol"), S298=(15.6165, "cal/(mol*K)"),
            ),
        )
        s_ch2x = Species(
            molecule=[m_ch2x],
            thermo=ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                Cpdata=([28.4959, 36.3588, 42.0219, 46.2428, 52.3978, 56.921, 64.1119], "J/(mol*K)"),
                H298=(0.654731, "kJ/mol"), S298=(19.8341, "J/(mol*K)"),
            ),
        )
        s_ch3x = Species(
            molecule=[m_ch3x],
            thermo=ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                Cpdata=([37.3325, 44.9406, 51.3613, 56.8115, 65.537, 72.3287, 83.3007], "J/(mol*K)"),
                H298=(-45.8036, "kJ/mol"), S298=(57.7449, "J/(mol*K)"),
            ),
        )

        # X=CH2 + H+ + e- <--> X-CH3
        cls.surface_rxn = Reaction(
            reactants=[s_proton, s_ch2x], products=[s_ch3x], electrons=-1,
            kinetics=SurfaceChargeTransfer(
                A=(2.483e21, "cm^3/(mol*s)"), V0=(0, "V"), Ea=(10, "kJ/mol"),
                Tmin=(200, "K"), Tmax=(3000, "K"), electrons=-1,
            ),
        )
        cls.sticking_rxn = Reaction(
            reactants=[s_proton, s_ch2x], products=[s_ch3x],
            kinetics=StickingCoefficient(A=0.1, n=0, Ea=(20, "kJ/mol"), T0=(1, "K")),
        )

    #: Frozen rates for the fixtures below, in the [m, mol, s] combination
    #: get_surface_rate_coefficient returns. A self-comparison between the dispatcher and the
    #: surface method would move together if both moved, so pin the absolute numbers too.
    SURFACE_RATE_AT_0V = 745832162915775.5
    SURFACE_RATE_AT_0_3V = 130820327780032.66
    STICKING_RATE = 3800697019.8819513

    def test_surface_charge_transfer_still_routes_to_the_surface_method(self):
        """Same value by the two routes, and the same absolute number as before the change."""
        rxn = self.surface_rxn

        assert rxn.get_rate_coefficient(1000.0, 1.0e5) == \
            rxn.get_surface_rate_coefficient(1000.0, 0.0, 0.0)
        assert rxn.get_rate_coefficient(1000.0, 1.0e5) == self.SURFACE_RATE_AT_0V

    def test_surface_charge_transfer_default_potential_is_still_zero_volts(self):
        """Surface kinetics keep their historical 0 V default; only the gas-phase default moved."""
        rxn = self.surface_rxn

        assert rxn.get_rate_coefficient(1000.0, 1.0e5) == \
            rxn.get_rate_coefficient(1000.0, 1.0e5, potential=0.0)

    def test_surface_charge_transfer_still_carries_its_potential(self):
        rxn = self.surface_rxn

        assert rxn.get_rate_coefficient(1000.0, 1.0e5, potential=0.3) == \
            rxn.get_surface_rate_coefficient(1000.0, 0.0, 0.3)
        assert rxn.get_rate_coefficient(1000.0, 1.0e5, potential=0.3) != \
            rxn.get_rate_coefficient(1000.0, 1.0e5)
        assert rxn.get_rate_coefficient(1000.0, 1.0e5, potential=0.3) == self.SURFACE_RATE_AT_0_3V

    def test_sticking_coefficient_still_demands_a_site_density(self):
        with pytest.raises(ValueError):
            self.sticking_rxn.get_rate_coefficient(1000.0, 1.0e5)

    def test_sticking_coefficient_still_routes_to_the_surface_method(self):
        rxn = self.sticking_rxn

        assert rxn.get_rate_coefficient(1000.0, 1.0e5, surface_site_density=2.72e-9) == \
            rxn.get_surface_rate_coefficient(1000.0, 2.72e-9)
        assert rxn.get_rate_coefficient(1000.0, 1.0e5, surface_site_density=2.72e-9) == \
            self.STICKING_RATE

    def test_diffusion_limiter_still_outranks_charge_transfer(self):
        """
        The new charge transfer branch sits *after* the diffusion limiter, so an enabled limiter
        still claims the reaction. If that order ever inverts, liquid-phase behaviour changes
        silently.

        The limiter is stubbed because its real internals cannot run on these kinetics at all:
        diffusionLimited.py calls ``kinetics.get_rate_coefficient(T, P=...)`` by keyword, and the
        charge transfer classes name that parameter V or dHrxn, so the real call raises TypeError.
        That is a separate, pre-existing defect; this test only pins the dispatch order.
        """
        rxn = Reaction(kinetics=make_charge_transfer_kinetics())
        sentinel = 1.2345e6
        calls = []

        original_enabled = diffusion_limiter.enabled
        original_get_effective_rate = diffusion_limiter.get_effective_rate
        try:
            diffusion_limiter.enabled = True
            diffusion_limiter.get_effective_rate = lambda reaction, T: calls.append((reaction, T)) or sentinel

            assert rxn.get_rate_coefficient(1000.0, 1.0e5) == sentinel
            assert len(calls) == 1
        finally:
            diffusion_limiter.enabled = original_enabled
            diffusion_limiter.get_effective_rate = original_get_effective_rate


@pytest.mark.database
class TestRealDatabaseChargeTransferReaction:
    """The same guarantee, on a reaction built from the real charge transfer families."""

    FAMILIES = [
        "Cation_Addition_MultipleBond",
        "Cation_Addition_MultipleBond_Disprop",
        "Cation_Li_Abstraction",
        "Cation_NO_Ring_Opening",
        "Cation_NO_Substitution",
    ]

    @classmethod
    def setup_class(cls):
        database = KineticsDatabase()
        database.load_families(
            path=os.path.join(settings["database.directory"], "kinetics", "families"),
            families=cls.FAMILIES,
        )

        cls.reactions = []
        for label in cls.FAMILIES:
            depository = database.families[label].get_training_depository()
            for entry in depository.entries.values():
                if isinstance(entry.data, ArrheniusChargeTransfer):
                    cls.reactions.append(
                        Reaction(
                            reactants=entry.item.reactants,
                            products=entry.item.products,
                            kinetics=entry.data,
                        )
                    )

        # The rate rules, as distinct from the training reactions, carry the Blowers-Masel form.
        cls.rules = []
        for label in cls.FAMILIES:
            for entries in database.families[label].rules.entries.values():
                for entry in entries:
                    if isinstance(entry.data, ArrheniusChargeTransferBM):
                        cls.rules.append(entry.data)

    def test_the_families_really_supply_charge_transfer_kinetics(self):
        """Guard on the fixture: an empty list would make every test below vacuously pass."""
        assert len(self.reactions) >= 24, (
            f"expected the cation families to supply ArrheniusChargeTransfer training reactions, "
            f"found {len(self.reactions)}"
        )
        assert len(self.rules) >= 28, (
            f"expected the cation families to supply ArrheniusChargeTransferBM rate rules, "
            f"found {len(self.rules)}"
        )

    def test_real_rate_rules_evaluate_and_ignore_pressure(self):
        """The Blowers-Masel rate rules, driven through a reaction that carries real thermo."""
        for rule in self.rules:
            rxn = make_endothermic_reaction(rule)
            rates = [rxn.get_rate_coefficient(1000.0, P) for P in PRESSURES]

            assert len(set(rates)) == 1, f"{rule} moved with pressure: {list(zip(PRESSURES, rates))}"
            # Not bit-identical: to_arrhenius_charge_transfer round-trips Ea through kJ/mol.
            assert rates[0] == pytest.approx(
                rule.get_rate_coefficient(1000.0, rxn.get_enthalpy_of_reaction(298.0)), rel=1e-12
            )

    def test_evaluation_does_not_mutate_the_stored_records(self):
        """
        The fix is call routing, not a data change: alpha, electrons, V0, Ea and E0 on the loaded
        records must read back identically after the rate has been evaluated.
        """
        before = [
            (k.alpha.value_si, k.electrons.value_si, k.V0.value_si,
             getattr(k, "Ea", None).value_si if hasattr(k, "Ea") else k.E0.value_si)
            for k in [r.kinetics for r in self.reactions] + self.rules
        ]

        for rxn in self.reactions:
            rxn.get_rate_coefficient(1000.0, 1.0e5)
            rxn.get_rate_coefficient(1000.0, 1.0e5, potential=0.1)
        for rule in self.rules:
            bm_rxn = make_endothermic_reaction(rule)
            bm_rxn.get_rate_coefficient(1000.0, 1.0e5)
            bm_rxn.get_rate_coefficient(1000.0, 1.0e5, potential=0.1)

        after = [
            (k.alpha.value_si, k.electrons.value_si, k.V0.value_si,
             getattr(k, "Ea", None).value_si if hasattr(k, "Ea") else k.E0.value_si)
            for k in [r.kinetics for r in self.reactions] + self.rules
        ]
        assert before == after

    def test_stored_electron_sign_and_alpha_are_what_the_records_say(self):
        """Guard that this ticket did not re-sign the fields a separate, merged patch corrected."""
        for kinetics in [r.kinetics for r in self.reactions] + self.rules:
            assert kinetics.alpha.value_si == 0.5
            assert kinetics.electrons.value_si == -1.0

    def test_real_reactions_evaluate_at_v0_and_ignore_pressure(self):
        for rxn in self.reactions:
            expected = rxn.kinetics.get_rate_coefficient(1000.0, rxn.kinetics.V0.value_si)
            assert expected > 0.0, f"degenerate fixture, zero rate at V0: {rxn}"
            for P in PRESSURES:
                assert rxn.get_rate_coefficient(1000.0, P) == expected, (
                    f"{rxn} moved with pressure {P} Pa"
                )

    def test_real_reactions_were_broken_before_the_fix(self):
        """
        Confirms the RED condition is really what this fix removes: feeding the pressure straight
        to the rate law -- which is what the old fall-through did -- annihilates the rate.
        """
        annihilated = [
            rxn for rxn in self.reactions
            if rxn.kinetics.get_rate_coefficient(1000.0, 1.0e5) == 0.0
        ]
        assert len(annihilated) == len(self.reactions), (
            f"only {len(annihilated)} of {len(self.reactions)} real reactions returned 0.0 when "
            "handed the pressure as a potential"
        )
