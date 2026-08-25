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
The reverse-reconstruction guard.

Every reactor reconstructs the reverse rate of a nominally reversible reaction as
``kb = kf(Tgas) / Keq(Tgas)``. That reconstruction is only defined when ``Keq(Tgas)``
prices every participant of the reaction at the state the forward rate was evaluated
at. These tests pin the two halves of that contract:

* a reaction for which it is *not* defined must hard-fail before integration, naming
  itself, rather than silently producing a ``kb``;
* ordinary thermal reversible chemistry must be numerically unchanged, to the bit.

The literals in :func:`test_thermal_control_is_bit_identical` were captured from the
unguarded base branch and must never move.
"""

import numpy as np
import pytest

import rmgpy.constants as constants
from rmgpy.exceptions import NonEquilibriumReverseRateError
from rmgpy.kinetics import Arrhenius, Marcus, SurfaceArrhenius, SurfaceChargeTransfer
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.solver.liquid import LiquidReactor
from rmgpy.solver.mbSampled import MBSampledReactor
from rmgpy.solver.simple import SimpleReactor
from rmgpy.solver.surface import SurfaceReactor
from rmgpy.species import Species
from rmgpy.thermo import ThermoData


def _thermo(h298_kJ, s298, cp=30.0):
    """A constant-Cp surrogate, so these tests do not depend on database thermo."""
    return ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
        Cpdata=([cp] * 7, "J/(mol*K)"),
        H298=(h298_kJ, "kJ/mol"),
        S298=(s298, "J/(mol*K)"),
    )


def _species(label, adjlist=None, smiles=None, h298_kJ=0.0, s298=100.0):
    if adjlist is not None:
        mol = Molecule().from_adjacency_list(adjlist)
    else:
        mol = Molecule(smiles=smiles)
    spc = Species(label=label, molecule=[mol])
    spc.thermo = _thermo(h298_kJ, s298)
    return spc


def _marcus():
    """The rate rule carried by the Cation_R_Recombination root node."""
    return Marcus(
        A=(1.73e06, "m^3/(mol*s)"),
        n=2,
        lmbd_i_coefs=[51902.4, -1.10249, -0.00813508, 3.26585e-06],
        beta=(1.2e10, "m^-1"),
        wr=(0, "J/mol"),
        wp=(0, "J/mol"),
        lmbd_o=(0, "J/mol"),
    )


def _cation_recombination():
    """
    Li+ + CH3 <=> CH3Li -- the shape of the twelve Cation_R_Recombination entries:
    Marcus kinetics, a metadata-only electron (``electrons = -1``), and the family's
    ``reversible = True``.
    """
    li_cation = _species("Li+", adjlist="1 Li u0 p0 c+1", h298_kJ=680.0, s298=133.0)
    ch3 = _species("CH3", smiles="[CH3]", h298_kJ=146.0, s298=194.0)
    ch3li = _species("CH3Li", smiles="C[Li]", h298_kJ=120.0, s298=230.0)
    rxn = Reaction(
        reactants=[li_cation, ch3],
        products=[ch3li],
        kinetics=_marcus(),
        reversible=True,
        electrons=-1,
    )
    return rxn, [li_cation, ch3, ch3li]


def _thermal_control():
    """CH4 + H <=> CH3 + H2 -- ordinary, uncharged, Tgas-only reversible chemistry."""
    ch3 = _species("CH3", smiles="[CH3]", h298_kJ=146.0, s298=194.0)
    h = _species("H", smiles="[H]", h298_kJ=218.0, s298=114.7)
    h2 = _species("H2", smiles="[H][H]", h298_kJ=0.0, s298=130.7)
    ch4 = _species("CH4", smiles="C", h298_kJ=-74.6, s298=186.3)
    rxn = Reaction(
        reactants=[ch4, h],
        products=[ch3, h2],
        kinetics=Arrhenius(A=(4.1e03, "m^3/(mol*s)"), n=1.5, Ea=(35.0, "kJ/mol"), T0=(1, "K")),
        reversible=True,
    )
    return rxn, [ch3, h, h2, ch4]


class ReverseReconstructionPredicateTest:
    """The predicate itself, independent of any reactor."""

    def test_metadata_electron_is_refused(self):
        rxn, _ = _cation_recombination()
        reason = rxn.get_reverse_from_equilibrium_refusal()
        assert reason is not None
        assert "electron" in reason.lower()

    def test_explicit_electron_is_refused(self):
        electron = Species(label="e", molecule=[Molecule().from_adjacency_list("1 e u0 p0 c-1")])
        electron.thermo = _thermo(0.0, 20.0)
        ar = _species("Ar", smiles="[Ar]", h298_kJ=0.0, s298=154.8)
        ar_cation = _species("Ar+", adjlist="1 Ar u1 p3 c+1", h298_kJ=1520.0, s298=160.0)
        rxn = Reaction(
            reactants=[ar_cation, electron],
            products=[ar],
            kinetics=Arrhenius(A=(1e06, "m^3/(mol*s)"), n=0, Ea=(0, "kJ/mol"), T0=(1, "K")),
            reversible=True,
        )
        reason = rxn.get_reverse_from_equilibrium_refusal()
        assert reason is not None
        assert "electron" in reason.lower()

    def test_thermal_reaction_is_permitted(self):
        rxn, _ = _thermal_control()
        assert rxn.get_reverse_from_equilibrium_refusal() is None

    def test_irreversible_declaration_is_not_the_guard_s_business(self):
        """The refusal is about the reversal, so an irreversible entry has nothing to refuse."""
        rxn, _ = _cation_recombination()
        rxn.reversible = False
        assert rxn.get_reverse_from_equilibrium_refusal() is None

    def test_charge_transfer_kinetics_with_a_potential_reference_is_permitted(self):
        """
        Electrochemistry is the one supported charge-transfer reverse: the kinetics
        declares the electrode potential at which the electron is priced, so the
        reversal is referenced rather than unreferenced.
        """
        site = _species("X", adjlist="1 X u0 p0 c0", h298_kJ=0.0, s298=0.0)
        adsorbate = _species(
            "HX",
            adjlist="1 H u0 p0 c0 {2,S}\n2 X u0 p0 c0 {1,S}",
            h298_kJ=-40.0,
            s298=20.0,
        )
        proton = _species("H+", adjlist="1 H u0 p0 c+1", h298_kJ=0.0, s298=20.0)
        rxn = Reaction(
            reactants=[proton, site],
            products=[adsorbate],
            kinetics=SurfaceChargeTransfer(
                A=(1e06, "m^3/(mol*s)"), n=0, Ea=(0, "kJ/mol"), electrons=1, alpha=0.5
            ),
            reversible=True,
            electrons=-1,
        )
        assert rxn.get_reverse_from_equilibrium_refusal() is None

    def test_electron_temperature_kinetics_is_refused(self):
        """
        The predicate is on the property, not the class: a Tgas-only-invalid rate law
        is refused even with no electron anywhere in the reaction.
        """
        rxn, _ = _thermal_control()

        class _TwoTemperatureArrhenius(Arrhenius):
            pass

        kin = _TwoTemperatureArrhenius(A=(1.0, "m^3/(mol*s)"), n=0, Ea=(0, "kJ/mol"), T0=(1, "K"))
        kin.uses_electron_temperature = True
        rxn.kinetics = kin
        reason = rxn.get_reverse_from_equilibrium_refusal()
        assert reason is not None
        assert "Tgas" in reason or "temperature" in reason.lower()


class ReverseReconstructionGuardInReactorsTest:
    """The refusal must fire in every reactor that performs the reconstruction."""

    T = (1000, "K")
    P = (1e5, "Pa")

    def _forbidden_and_control(self):
        forbidden, forbidden_species = _cation_recombination()
        control, control_species = _thermal_control()
        # _cation_recombination and _thermal_control each build their own CH3; merge
        # them so the reactor sees a single consistent species list.
        species = forbidden_species + [s for s in control_species if s.label != "CH3"]
        control.reactants = [s for s in species if s.label == "CH4"] + [
            s for s in species if s.label == "H"
        ]
        control.products = [s for s in species if s.label == "CH3"] + [
            s for s in species if s.label == "H2"
        ]
        return forbidden, control, species

    def test_simple_reactor_refuses(self):
        forbidden, _, species = self._forbidden_and_control()
        reactor = SimpleReactor(
            T=self.T,
            P=self.P,
            initial_mole_fractions={s: 1.0 / len(species) for s in species},
            termination=[],
        )
        with pytest.raises(NonEquilibriumReverseRateError) as exc:
            reactor.initialize_model(species, [forbidden], [], [])
        assert "Li+" in str(exc.value) and "CH3Li" in str(exc.value)

    def test_liquid_reactor_refuses(self):
        forbidden, _, species = self._forbidden_and_control()
        reactor = LiquidReactor(
            T=self.T,
            initial_concentrations={s: 1.0 for s in species},  # LiquidReactor takes SI
            termination=[],
        )
        with pytest.raises(NonEquilibriumReverseRateError) as exc:
            reactor.initialize_model(species, [forbidden], [], [])
        assert "Li+" in str(exc.value) and "CH3Li" in str(exc.value)

    def test_mb_sampled_reactor_refuses(self):
        forbidden, _, species = self._forbidden_and_control()
        reactor = MBSampledReactor(
            T=self.T,
            P=self.P,
            initial_mole_fractions={s: 1.0 / len(species) for s in species},
            k_sampling=(1e-3, "s^-1"),
            constantSpeciesList=[],
            termination=[],
        )
        with pytest.raises(NonEquilibriumReverseRateError) as exc:
            reactor.initialize_model(species, [forbidden], [], [])
        assert "Li+" in str(exc.value) and "CH3Li" in str(exc.value)

    def test_surface_reactor_refuses(self):
        forbidden, _, species = self._forbidden_and_control()
        site = _species("X", adjlist="1 X u0 p0 c0", h298_kJ=0.0, s298=0.0)
        reactor = SurfaceReactor(
            T=self.T,
            P_initial=self.P,
            n_sims=1,
            initial_gas_mole_fractions={s: 1.0 / len(species) for s in species},
            initial_surface_coverages={site: 1.0},
            surface_volume_ratio=(1e1, "m^-1"),
            surface_site_density=(2.72e-9, "mol/cm^2"),
            termination=[],
        )
        with pytest.raises(NonEquilibriumReverseRateError) as exc:
            reactor.initialize_model(species + [site], [forbidden], [], [])
        assert "Li+" in str(exc.value) and "CH3Li" in str(exc.value)

    def test_thermal_control_is_bit_identical(self):
        """
        The negative control, and the reason this guard is narrow: an ordinary thermal
        reversible reaction must integrate exactly as it did before the guard existed.
        These literals were captured on the unguarded base branch.
        """
        _, control, species = self._forbidden_and_control()
        reactor = SimpleReactor(
            T=self.T,
            P=self.P,
            initial_mole_fractions={s: 1.0 / len(species) for s in species},
            termination=[],
        )
        reactor.initialize_model(species, [control], [], [])
        j = reactor.reaction_index[control]
        assert reactor.kf[j] == 1925791.2777157801
        assert reactor.Keq[j] == 12.651095698063447
        assert reactor.kb[j] == 152223.2796018268

    def test_irreversible_cation_reaction_still_integrates(self):
        """
        The refusal is not a ban on the chemistry. Declared irreversible, the same
        reaction goes through untouched with kb = 0 -- which is what the ruling
        requires the *mechanism* to say, rather than the integrator to assume.
        """
        forbidden, _, species = self._forbidden_and_control()
        forbidden.reversible = False
        reactor = SimpleReactor(
            T=self.T,
            P=self.P,
            initial_mole_fractions={s: 1.0 / len(species) for s in species},
            termination=[],
        )
        reactor.initialize_model(species, [forbidden], [], [])
        j = reactor.reaction_index[forbidden]
        assert reactor.kb[j] == 0.0
        assert np.isinf(reactor.Keq[j])
        assert reactor.kf[j] == 8.2990941663987537e-74
