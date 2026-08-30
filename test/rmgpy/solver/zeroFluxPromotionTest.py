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
The zero-flux escape hatch in :meth:`ReactionSystem.simulate`.

When the characteristic rate is exactly zero every edge and network rate ratio is a
0/0 division, so no object can ever exceed a tolerance and the model could never
grow.  ``simulate`` breaks that singularity by force-promoting the edge species with
the largest rate.  That is a defensible choice while the edge carries flux, and no
choice at all when it does not: ``np.argmax`` over an all-zero array returns 0, so
the species promoted is whichever one happens to sit first in the edge list.

These tests pin the halves of the boundary.  The edge list is deliberately ordered so
that the largest-rate species is *not* the first one, which is what makes "chosen by
rate" distinguishable from "chosen by position".  A third test guards the F5 defect:
the inert test above uses an edge rate of *exactly* zero, but the promotion guard used
to test ``max_species_rate == 0`` exactly, so an edge species carrying subnormal
numerical dust (a rate in ``(0, np.finfo(float).tiny)``, physically inert but not
``== 0``) walked past it and was promoted.  The guard now rejects the whole subnormal
band, and that test constructs a rate inside it.

A fourth test guards the char_rate sibling of the F5 defect, one layer up: ``char_rate``
(the L2 norm of the core species rates) used to gate promotion with ``char_rate == 0``
exactly.  A core carrying a single reaction with a tiny enough A-factor produces a
``char_rate`` that is nonzero yet still numerical dust -- normal-magnitude dust, not
F5's subnormal band, because ``sqrt`` of any tiny positive sum-of-squares returns a
normal double.  The old exact-zero test let that case fall through to ratio-based
promotion, dividing by the same near-nothing ``char_rate`` and promoting an arbitrary
species on the resulting garbage ratio.  The guard now rejects the whole dust band
below ``_CHAR_RATE_FLOOR``, and that test constructs a ``char_rate`` inside it.
"""

import logging

import numpy as np

from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
from rmgpy.solver.base import TerminationTime, _CHAR_RATE_FLOOR
from rmgpy.solver.simple import SimpleReactor
from rmgpy.species import Species
from rmgpy.thermo import ThermoData

TDATA = ([300, 400, 500, 600, 800, 1000, 1500], "K")


def _species(smiles, cpdata, h298, s298):
    return Species(
        molecule=[Molecule().from_smiles(smiles)],
        thermo=ThermoData(
            Tdata=TDATA,
            Cpdata=(cpdata, "cal/(mol*K)"),
            H298=(h298, "kcal/mol"),
            S298=(s298, "cal/(mol*K)"),
        ),
    )


class ZeroFluxPromotionTest:
    def _build(self, a_unimolecular, a_bimolecular):
        """
        A reactor whose core carries no reactions at all, so the characteristic rate
        is zero by construction, and whose two edge reactions have their A factors
        supplied by the caller.  Setting both to zero makes the whole system inert.
        """
        ch4 = _species("C", [8.615, 9.687, 10.963, 12.301, 14.841, 16.976, 20.528], -17.714, 44.472)
        c2h6 = _species("CC", [12.684, 15.506, 18.326, 20.971, 25.500, 29.016, 34.595], -19.521, 54.799)
        ch3 = _species("[CH3]", [9.397, 10.123, 10.856, 11.571, 12.899, 14.055, 16.195], 9.357, 45.174)
        c2h5 = _species("C[CH2]", [11.635, 13.744, 16.085, 18.246, 21.885, 24.676, 29.107], 29.496, 56.687)
        h2 = _species("[H][H]", [6.895, 6.975, 6.994, 7.009, 7.081, 7.219, 7.720], 0.0, 31.233)

        # C2H6 -> 2 CH3 makes CH3 the fastest-produced edge species when it runs.
        edge_rxn1 = Reaction(
            reactants=[c2h6],
            products=[ch3, ch3],
            kinetics=Arrhenius(A=(a_unimolecular, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
        )
        # CH4 + C2H6 -> CH3 + C2H5 + H2 brings the other two species onto the edge.
        edge_rxn2 = Reaction(
            reactants=[ch4, c2h6],
            products=[ch3, c2h5, h2],
            kinetics=Arrhenius(A=(a_bimolecular, "m^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
        )

        core_species = [ch4, c2h6]
        core_reactions = []
        # CH3 is placed last on purpose: an argmax over zeros would return C2H5.
        edge_species = [c2h5, h2, ch3]
        edge_reactions = [edge_rxn1, edge_rxn2]

        rxn_system = SimpleReactor(
            T=1000.0,
            P=1.0e5,
            initial_mole_fractions={ch4: 0.5, c2h6: 0.5},
            n_sims=1,
            termination=[TerminationTime((1e-6, "s"))],
        )
        rxn_system.initialize_model(core_species, core_reactions, edge_species, edge_reactions)
        return rxn_system, core_species, core_reactions, edge_species, edge_reactions

    @staticmethod
    def _simulate(rxn_system, core_species, core_reactions, edge_species, edge_reactions):
        return rxn_system.simulate(
            core_species,
            core_reactions,
            edge_species,
            edge_reactions,
            [],
            [],
            model_settings=ModelSettings(tol_keep_in_edge=0, tol_move_to_core=1e5, tol_interrupt_simulation=1e5),
            simulator_settings=SimulatorSettings(),
        )

    def test_no_promotion_when_the_whole_system_is_inert(self, caplog):
        """
        With every rate coefficient zero there is no flux anywhere, so there is no
        rate comparison to make and nothing may be promoted into the core.  The
        condition is reported rather than papered over.
        """
        system = self._build(a_unimolecular=0.0, a_bimolecular=0.0)
        rxn_system, core_species, core_reactions, edge_species, edge_reactions = system

        with caplog.at_level(logging.INFO):
            terminated, resurrected, invalid_objects, _surf_spc, _surf_rxn, _t, _x = self._simulate(*system)

        assert np.all(rxn_system.core_species_rates == 0.0), "the core is meant to be inert here"
        assert np.all(rxn_system.edge_species_rates == 0.0), "the edge is meant to be inert here"
        assert invalid_objects == [], (
            "a system with no flux anywhere promoted {0!r} into the core; with every edge rate "
            "at zero there is no largest-rate species, so this is the first edge species chosen "
            "on no physical basis".format(invalid_objects)
        )
        assert not any(
            "added to model core to avoid singularity" in record.getMessage() for record in caplog.records
        ), "the singularity escape hatch reported a promotion it should not have made"
        zero_flux_records = [
            record
            for record in caplog.records
            if record.levelno >= logging.ERROR and record.getMessage().startswith("ZERO FLUX:")
        ]
        assert zero_flux_records, "the zero-flux condition was never reported: {0!r}".format(
            [r.getMessage() for r in caplog.records]
        )
        message = zero_flux_records[0].getMessage()
        assert "NO species was added to the model core" in message, message
        assert "SimpleReactor" in message, message

    def test_promotion_still_picks_the_fastest_species_when_the_edge_has_flux(self):
        """
        The escape hatch keeps its job in the ordinary case: the characteristic rate is
        still zero, but the edge carries flux, so the largest-rate species is promoted.
        Only C2H6 -> 2 CH3 runs, so that species is CH3 -- last in the edge list, which
        is what distinguishes a rate comparison from an argmax over zeros.
        """
        system = self._build(a_unimolecular=1.0e6, a_bimolecular=0.0)
        rxn_system, core_species, core_reactions, edge_species, edge_reactions = system

        terminated, resurrected, invalid_objects, _surf_spc, _surf_rxn, _t, _x = self._simulate(*system)

        assert np.all(rxn_system.core_species_rates == 0.0), "the core is meant to be inert here"
        assert rxn_system.edge_species_rates.max() > 0.0, "the edge is meant to carry flux here"
        assert invalid_objects == [edge_species[2]], (
            "expected the fastest edge species (CH3) to be promoted, got {0!r}".format(invalid_objects)
        )

    def test_no_promotion_when_the_only_edge_flux_is_subnormal_dust(self, caplog):
        """
        The F5 defect: the characteristic rate is exactly zero and the edge carries a
        rate that is *nonzero* but subnormal -- underflow/cancellation dust, physically
        indistinguishable from inert.  The old guard tested ``max_species_rate == 0``
        exactly, so this rate slipped past and the species was force-promoted on no
        physical basis.  With a tiny unimolecular A factor, C2H6 -> 2 CH3 produces a CH3
        rate down in the subnormal band ``(0, np.finfo(float).tiny)``; the guard must now
        recognise it as dust, promote nothing, and report the zero-flux condition.
        """
        system = self._build(a_unimolecular=1.0e-315, a_bimolecular=0.0)
        rxn_system, core_species, core_reactions, edge_species, edge_reactions = system

        with caplog.at_level(logging.INFO):
            terminated, resurrected, invalid_objects, _surf_spc, _surf_rxn, _t, _x = self._simulate(*system)

        max_edge_rate = float(np.abs(rxn_system.edge_species_rates).max())
        assert 0.0 < max_edge_rate < np.finfo(np.float64).tiny, (
            "this test needs the largest edge rate to be nonzero subnormal dust so that the "
            "old exact-zero guard would have promoted it; got {0!r} (tiny={1!r})".format(
                max_edge_rate, np.finfo(np.float64).tiny
            )
        )
        assert np.all(rxn_system.core_species_rates == 0.0), "the core is meant to be inert here"
        assert invalid_objects == [], (
            "a system whose only edge flux is subnormal dust promoted {0!r} into the core; a "
            "rate below the smallest normal float is numerical dust, not a physical basis for "
            "promotion".format(invalid_objects)
        )
        assert not any(
            "added to model core to avoid singularity" in record.getMessage() for record in caplog.records
        ), "the singularity escape hatch promoted a species whose rate is numerical dust"
        zero_flux_records = [
            record
            for record in caplog.records
            if record.levelno >= logging.ERROR and record.getMessage().startswith("ZERO FLUX:")
        ]
        assert zero_flux_records, "the zero-flux condition was never reported: {0!r}".format(
            [r.getMessage() for r in caplog.records]
        )
        assert "numerical dust" in zero_flux_records[0].getMessage(), zero_flux_records[0].getMessage()

    def test_no_promotion_when_core_char_rate_is_normal_magnitude_dust(self, caplog):
        """
        The char_rate sibling of the F5 defect: ``char_rate`` is nonzero but is itself
        numerical dust, not a physical rate.  A single core reaction (C2H6 -> 2 CH3) with
        an A-factor of 1e-105 produces a core rate small enough that the L2 norm
        ``char_rate = sqrt(sum(rate_i**2))`` lands strictly inside ``(0, _CHAR_RATE_FLOOR)``
        -- nonzero, so the old ``char_rate == 0`` gate would have skipped it and fallen
        through to ratio-based promotion (dividing by that same near-nothing char_rate).
        The edge is fully inert (no edge reactions), so with the widened gate this must
        take the same zero-flux path as a truly all-zero core: nothing promoted, and the
        condition reported rather than papered over.
        """
        ch4 = _species("C", [8.615, 9.687, 10.963, 12.301, 14.841, 16.976, 20.528], -17.714, 44.472)
        c2h6 = _species("CC", [12.684, 15.506, 18.326, 20.971, 25.500, 29.016, 34.595], -19.521, 54.799)
        ch3 = _species("[CH3]", [9.397, 10.123, 10.856, 11.571, 12.899, 14.055, 16.195], 9.357, 45.174)
        c2h5 = _species("C[CH2]", [11.635, 13.744, 16.085, 18.246, 21.885, 24.676, 29.107], 29.496, 56.687)
        h2 = _species("[H][H]", [6.895, 6.975, 6.994, 7.009, 7.081, 7.219, 7.720], 0.0, 31.233)

        core_rxn = Reaction(
            reactants=[c2h6],
            products=[ch3, ch3],
            kinetics=Arrhenius(A=(1e-105, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
        )

        core_species = [ch4, c2h6, ch3]
        core_reactions = [core_rxn]
        edge_species = [c2h5, h2]
        edge_reactions = []

        rxn_system = SimpleReactor(
            T=1000.0,
            P=1.0e5,
            initial_mole_fractions={ch4: 0.5, c2h6: 0.5},
            n_sims=1,
            termination=[TerminationTime((1e-6, "s"))],
        )
        rxn_system.initialize_model(core_species, core_reactions, edge_species, edge_reactions)

        with caplog.at_level(logging.INFO):
            terminated, resurrected, invalid_objects, _surf_spc, _surf_rxn, _t, _x = rxn_system.simulate(
                core_species,
                core_reactions,
                edge_species,
                edge_reactions,
                [],
                [],
                model_settings=ModelSettings(tol_keep_in_edge=0, tol_move_to_core=0.1, tol_interrupt_simulation=1e5),
                simulator_settings=SimulatorSettings(),
            )

        char_rate = float(np.sqrt(np.sum(rxn_system.core_species_rates ** 2)))
        assert 0.0 < char_rate < _CHAR_RATE_FLOOR, (
            "this test needs char_rate to be nonzero normal-magnitude dust so that the old "
            "exact-zero gate would have let it fall through to ratio-based promotion; got "
            "{0!r} (floor={1!r})".format(char_rate, _CHAR_RATE_FLOOR)
        )
        assert invalid_objects == [], (
            "a core whose only nonzero rate is dust-magnitude char_rate promoted {0!r} into the "
            "core; a char_rate this small is numerical dust, not a physical basis for ratio-based "
            "promotion".format(invalid_objects)
        )
        assert not any(
            "added to model core to avoid singularity" in record.getMessage() for record in caplog.records
        ), "the singularity escape hatch promoted a species on a dust-magnitude char_rate"
        zero_flux_records = [
            record
            for record in caplog.records
            if record.levelno >= logging.ERROR and record.getMessage().startswith("ZERO FLUX:")
        ]
        assert zero_flux_records, "the zero-flux condition was never reported: {0!r}".format(
            [r.getMessage() for r in caplog.records]
        )
