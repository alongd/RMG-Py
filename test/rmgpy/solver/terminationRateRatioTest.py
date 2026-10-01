#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)    #
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
Regression tests for :class:`TerminationRateRatio`.

Defect under test
-----------------
``char_rate`` (``base.pyx``) is the L2 norm of the *net* core species production
rates, and ``max_char_rate`` is its running maximum. The main ``simulate()``
loop sets ``first_time = True`` and skips ``self.step()`` on its first pass, so
``max_char_rate`` is seeded from ``char_rate`` evaluated at ``t = 0`` exactly --
on the user's initial condition, before any integration has happened.

When a simulation starts from a composition that is not at partial equilibrium
and has no radical pool to drive an induction-then-ignition rise, the largest
characteristic flux of the entire trajectory is that ``t = 0`` relaxation of the
reactant into its own fast reversible channels. That relaxation is closed and
mass-conserving, so it converts essentially nothing. ``max_char_rate`` is then
pinned to a boundary value that is not a chemical event, the ratio decreases
monotonically from 1.0, and the criterion is guaranteed to fire during
induction -- reporting a converged simulation at ~0 conversion.

The documented meaning (``documentation/source/users/rmg/input.rst``) is the
flux "relative to the main chemical process", so the reference maximum must be
an *interior* maximum: evidence that a chemical event actually occurred.
"""

import numpy as np

from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
from rmgpy.solver.base import TerminationRateRatio, TerminationTime
from rmgpy.solver.simple import SimpleReactor
from rmgpy.species import Species
from rmgpy.thermo import ThermoData


TDATA = ([300, 400, 500, 600, 800, 1000, 1500], "K")
CPDATA = ([20.0, 24.0, 28.0, 31.0, 36.0, 39.0, 44.0], "cal/(mol*K)")


def _species(smiles, h298, s298=70.0):
    return Species(
        molecule=[Molecule().from_smiles(smiles)],
        thermo=ThermoData(Tdata=TDATA, Cpdata=CPDATA,
                          H298=(h298, "kcal/mol"), S298=(s298, "cal/(mol*K)")),
    )


def _build_induction_system(termination):
    """A closed-shell feed whose characteristic rate is maximal at t = 0.

    ``parent`` isomerises reversibly and very fast to a much less stable
    tautomer (``iso``), and separately decomposes slowly and irreversibly to
    ``prod`` -- the "real" chemistry. At t = 0 the tautomer is absent, so the
    isomerisation runs at full forward rate with zero reverse flux and
    dominates ``char_rate``. It equilibrates within picoseconds at negligible
    conversion (the tautomer is +40 kcal/mol, so K_eq is ~1e-7), after which
    the net flux collapses and only the slow decomposition remains. This is the
    structure of a neat-polymer pyrolysis deck, and of no ignition deck.
    """
    parent = _species("c1ccccc1C", 12.0)
    iso = _species("C=C1C=CC=CC1", 52.0)          # +40 kcal/mol tautomer
    prod = _species("c1ccccc1", 19.8)
    n2 = _species("N#N", 0.0)

    fast_isomerisation = Reaction(
        reactants=[parent], products=[iso], reversible=True,
        kinetics=Arrhenius(A=(1.0e13, "s^-1"), n=0, Ea=(40.0, "kcal/mol"), T0=(1, "K")),
    )
    slow_decomposition = Reaction(
        reactants=[parent], products=[prod], reversible=False,
        kinetics=Arrhenius(A=(1.0e13, "s^-1"), n=0, Ea=(77.0, "kcal/mol"), T0=(1, "K")),
    )

    core_species = [parent, iso, prod, n2]
    core_reactions = [fast_isomerisation, slow_decomposition]
    rxn_system = SimpleReactor(
        T=1300.0, P=1.0e5,
        initial_mole_fractions={parent: 0.1, iso: 0.0, prod: 0.0, n2: 0.9},
        n_sims=1,
        termination=termination,
    )
    rxn_system.initialize_model(core_species, core_reactions, [], [])
    return rxn_system, core_species, core_reactions, parent


def _run(termination):
    rxn_system, core_species, core_reactions, parent = _build_induction_system(termination)
    i_parent = core_species.index(parent)
    y0 = rxn_system.y.copy()
    rxn_system.simulate(
        core_species, core_reactions, [], [], [], [],
        model_settings=ModelSettings(tol_keep_in_edge=0, tol_move_to_core=1, tol_interrupt_simulation=1e100),
        simulator_settings=SimulatorSettings(),
    )
    conversion = 1.0 - rxn_system.y[i_parent] / y0[i_parent]
    return rxn_system.t, conversion


class TestTerminationRateRatioInductionRegime:
    """The characteristic rate is maximal at t = 0 here, so its decay is not
    evidence that the chemistry finished."""

    def test_char_rate_maximum_is_at_t_zero(self):
        """Premise check: this system's char_rate really does peak at t = 0.

        If this ever fails the fixture has stopped exercising the defect and
        the assertions below are meaningless.
        """
        rxn_system, core_species, core_reactions, _ = _build_induction_system([])
        rates = []
        for t in [0.0] + list(np.logspace(-14, -3, 60)):
            if t > 0.0:
                rxn_system.advance(t)
            else:
                rxn_system.y = rxn_system.y  # evaluate the initial condition as-is
            rxn_system.residual(rxn_system.t, rxn_system.y, np.zeros(rxn_system.y.shape))
            rates.append(float(np.sqrt(np.sum(rxn_system.core_species_rates ** 2))))
        assert rates[0] == max(rates), (
            "fixture no longer peaks at t=0: max is at index "
            f"{int(np.argmax(rates))} of {len(rates)}"
        )

    def test_rate_ratio_does_not_terminate_during_induction(self):
        """The defect: with a time backstop of 1 s and a rate ratio of 0.1, the
        rate ratio must not pre-empt the time criterion at ~zero conversion.

        Before the fix this terminates in the picosecond range having converted
        essentially none of the parent -- the tautomer equilibration alone.
        """
        t_end, conversion = _run([TerminationTime((1.0, "s")), TerminationRateRatio(0.1)])
        assert t_end >= 1.0, (
            f"simulation stopped at t={t_end:.4e} s, before the 1 s termination time: "
            "TerminationRateRatio fired during the induction period"
        )
        assert conversion > 0.5, (
            f"parent conversion was only {conversion:.6e}; the simulation was cut off "
            "before the chemistry it was meant to explore"
        )

    def test_rate_ratio_still_fires_after_a_real_peak(self):
        """The fix must not disable the criterion where it is meaningful.

        An ordinary radical-chain pyrolysis: slow homolysis of the fuel feeds a
        radical that is consumed much faster than it is made. The radical is
        absent at t = 0, so the characteristic rate *rises* as its pool fills
        and the fast step switches on, giving a genuine interior maximum. It
        then decays as the fuel depletes, and the criterion must still fire on
        that decay -- this is the case terminationRateRatio exists for.
        """
        fuel = _species("CCCC", -30.0)          # n-butane
        rad = _species("C[CH2]", 28.4)          # ethyl
        olefin = _species("C=C", 12.5)          # ethylene
        h_atom = _species("[H]", 52.1, s298=27.4)
        n2 = _species("N#N", 0.0)
        homolysis = Reaction(                    # C4H10 -> 2 C2H5, slow initiation
            reactants=[fuel], products=[rad, rad], reversible=False,
            kinetics=Arrhenius(A=(1.0e13, "s^-1"), n=0, Ea=(82.0, "kcal/mol"), T0=(1, "K")),
        )
        beta_scission = Reaction(                # C2H5 -> C2H4 + H, fast
            reactants=[rad], products=[olefin, h_atom], reversible=False,
            kinetics=Arrhenius(A=(1.0e13, "s^-1"), n=0, Ea=(40.0, "kcal/mol"), T0=(1, "K")),
        )
        core_species = [fuel, rad, olefin, h_atom, n2]
        core_reactions = [homolysis, beta_scission]
        rxn_system = SimpleReactor(
            T=1300.0, P=1.0e5,
            initial_mole_fractions={fuel: 0.1, rad: 0.0, olefin: 0.0, h_atom: 0.0, n2: 0.9},
            n_sims=1,
            termination=[TerminationTime((1.0e4, "s")), TerminationRateRatio(0.1)],
        )
        rxn_system.initialize_model(core_species, core_reactions, [], [])
        rxn_system.simulate(
            core_species, core_reactions, [], [], [], [],
            model_settings=ModelSettings(tol_keep_in_edge=0, tol_move_to_core=1, tol_interrupt_simulation=1e100),
            simulator_settings=SimulatorSettings(),
        )
        assert rxn_system.t < 1.0e4, (
            "TerminationRateRatio failed to fire on a genuine post-peak decay; "
            f"ran to the {rxn_system.t:.4e} s time backstop instead"
        )

    def test_rate_ratio_only_deck_is_unchanged(self):
        """The suppression must NOT apply when terminationRateRatio is the only
        criterion.

        `test/regression/oxidation/input.py` and `test/regression/nitrogen/input.py`
        are exactly this shape. For them the rate ratio is the criterion of last
        resort: suppressing it would leave `simulate()` with no reachable bound and
        would also keep `terminated` False, which switches off edge pruning
        (`main.py`, gated on `all_terminated`). So with no TerminationTime backstop
        the old behaviour must survive verbatim, boundary maximum or not.
        """
        t_end, conversion = _run([TerminationRateRatio(0.1)])
        assert t_end < 1.0e-6, (
            f"a rate-ratio-only system ran to t={t_end:.4e} s; the criterion of last "
            "resort was suppressed and nothing else can stop this simulation"
        )
