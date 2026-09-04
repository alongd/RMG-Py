#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@northeastern.edu) and the RMG Team            #
# (rmg_dev@mit.edu)                                                           #
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
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL     #
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER  #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                    #
#                                                                             #
###############################################################################

"""
Characterization tests for issue i209 -- "source priority is an accident of
ordering, not a guarantee".

RMG decides two generated reactions are *the same reaction* by comparing only
the reactant/product species references and a per-side electron-placement count
(:func:`rmgpy.rmg.model.are_identical_species_references`). Kinetics, provenance,
and measured-vs-estimated status are never consulted. On a collision,
:meth:`CoreEdgeReactionModel.check_for_existing_reaction` returns the reaction
that was registered FIRST (the incumbent) and
:meth:`CoreEdgeReactionModel.make_new_reaction` hands the caller that incumbent,
discarding the newcomer. ``are_identical_species_references``' own docstring
states it: "the model then silently keeps whichever was offered first."

These tests RECORD that behaviour as it stands today, so that a future change to
the duplicate path (e.g. one that started preferring a sourced value over an
estimate) is caught by a failing characterization test rather than shipping
silently. They assert nothing about what the code *should* do.

The load-bearing finding is the pair
``test_sourced_value_wins_when_registered_first`` /
``test_estimate_wins_when_registered_first``: the SAME predicate keeps the
sourced value in one ordering and the estimate in the other. Which value
survives is decided purely by registration order -- and RMG's initialize()
sequence registers seed mechanisms (core) and reaction libraries (edge) before
any family reaction is generated, so a family estimate normally loses to a
library. A seed mechanism carrying a previously-estimated rate is registered
*ahead* of reaction libraries, which is the ordering the second test pins.
"""

import os

import rmgpy.data.rmg
from rmgpy import settings
from rmgpy.data.base import ForbiddenStructures
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.electron_balance import get_electron_placement_counts
from rmgpy.kinetics import Arrhenius
from rmgpy.rmg.model import CoreEdgeReactionModel, are_identical_species_references
from rmgpy.species import Species

#: A deliberately "measured/sourced" rate: fast, near-zero barrier.
SOURCED = Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(5.0, "kcal/mol"))
#: A deliberately "family-estimated" rate: different everywhere from SOURCED.
ESTIMATE = Arrhenius(A=(9.9e9, "cm^3/(mol*s)"), n=0.5, Ea=(30.0, "kcal/mol"))

# Two library labels that resolve through get_family_library_object, and one
# family that generates the same heavy-atom transformation. These live in the
# committed testing_database, so the test needs no RMG-database checkout.
SOURCED_LIBRARY = "GRI-Mech3.0"
SEED_LIBRARY = "ethane-oxidation"
FAMILY = "H_Abstraction"


class _SingletonKineticsDatabase:
    """Register a kinetics database as the process-global singleton and restore
    whatever was there before.

    ``check_for_existing_reaction`` resolves a reaction's owner via
    ``get_family_library_object``, which reads ``rmgpy.data.rmg.database`` -- not
    the model's own attribute -- so a test exercising the real model path must
    register one. Mirrors the helper in ``i134DuplicateElectronsTest``.
    """

    def __init__(self, kinetics_database):
        self._kinetics_database = kinetics_database
        self._saved = None

    def __enter__(self):
        self._saved = rmgpy.data.rmg.database
        database = RMGDatabase()
        database.kinetics = self._kinetics_database
        database.forbidden_structures = ForbiddenStructures()
        rmgpy.data.rmg.database = database
        return database

    def __exit__(self, *exc_info):
        rmgpy.data.rmg.database = self._saved
        return False


class TestSourcePriorityIsOrdering:
    """The collapse recorded end to end through the model's real duplicate path."""

    @classmethod
    def setup_class(cls):
        path = os.path.join(settings["test_data.directory"], "testing_database")
        cls.database = KineticsDatabase()
        cls.database.load(
            os.path.join(path, "kinetics"),
            families=[FAMILY],
            libraries=[SOURCED_LIBRARY, SEED_LIBRARY],
        )

    def _species(self):
        """Fresh Species objects per test. ``make_new_reaction`` /
        ``check_for_existing_reaction`` compare species by reference, so sharing
        objects across tests would leak one test's identities into the next."""
        a = Species().from_smiles("[H]")
        b = Species().from_smiles("C=C[CH2]C")
        c = Species().from_smiles("C=C=CC")
        d = Species().from_smiles("[H][H]")
        a.label, b.label, c.label, d.label = "[H]", "C=C[CH2]C", "C=C=CC", "[H][H]"
        b.generate_resonance_structures()
        return a, b, c, d

    def _model(self, *species):
        cerm = CoreEdgeReactionModel()
        for spc in species:
            cerm.add_species_to_core(spc)
        return cerm

    @staticmethod
    def _register(cerm, rxn):
        rxn.reactants.sort()
        rxn.products.sort()
        cerm.add_reaction_to_core(rxn)
        cerm.register_reaction(rxn)

    @staticmethod
    def _library(reactants, products, library, kinetics):
        return LibraryReaction(reactants=list(reactants), products=list(products),
                               library=library, kinetics=kinetics, reversible=True)

    @staticmethod
    def _template(reactants, products, kinetics):
        return TemplateReaction(reactants=list(reactants), products=list(products),
                                family=FAMILY, template=["Csd", "H"], kinetics=kinetics)

    def _singleton(self):
        return _SingletonKineticsDatabase(self.database)

    def test_identity_ignores_kinetics_and_provenance(self):
        """A library reaction (sourced) and a family reaction (estimated) over
        the same heavy species, with wildly different kinetics and different
        owners, are declared THE SAME reaction. The neutral electron count is
        ``(0, 0)`` on both sides, so identity reduces to the heavy-species
        reference comparison."""
        a, b, c, d = self._species()
        sourced = self._library([a, b], [c, d], SOURCED_LIBRARY, SOURCED)
        estimate = self._template([a, b], [c, d], ESTIMATE)
        for rxn in (sourced, estimate):
            rxn.reactants.sort()
            rxn.products.sort()
        assert get_electron_placement_counts(sourced) == (0, 0)
        assert get_electron_placement_counts(estimate) == (0, 0)
        assert are_identical_species_references(sourced, estimate)
        # The predicate is blind to which rate is which.
        assert sourced.kinetics.A.value_si != estimate.kinetics.A.value_si

    def test_sourced_value_wins_when_registered_first(self):
        """Q1 normal ordering: a reaction library reaches the model before family
        enlargement, so the sourced library reaction is the incumbent and the
        family estimate is discarded. The model keeps the SOURCED rate."""
        with self._singleton():
            a, b, c, d = self._species()
            cerm = self._model(a, b, c, d)
            sourced = self._library([a, b], [c, d], SOURCED_LIBRARY, SOURCED)
            self._register(cerm, sourced)

            estimate = self._template([a, b], [c, d], ESTIMATE)
            found, kept = cerm.check_for_existing_reaction(estimate)

            assert found
            assert kept is sourced
            assert kept.kinetics.A.value_si == SOURCED.A.value_si

    def test_estimate_wins_when_registered_first(self):
        """Q3 inversion: a seed mechanism carrying a previously family-estimated
        rate is registered ahead of reaction libraries. Now the estimate is the
        incumbent, and the later sourced library reaction is discarded. The model
        keeps the ESTIMATE -- the sourced value never becomes operative. Nothing
        about the kinetics or provenance changes this verdict; only the order
        did."""
        with self._singleton():
            a, b, c, d = self._species()
            cerm = self._model(a, b, c, d)
            # A seed reaction is a LibraryReaction whose numbers were estimated in
            # a prior run and written to the seed mechanism.
            seed_estimate = self._library([a, b], [c, d], SEED_LIBRARY, ESTIMATE)
            self._register(cerm, seed_estimate)

            sourced = self._library([a, b], [c, d], SOURCED_LIBRARY, SOURCED)
            found, kept = cerm.check_for_existing_reaction(sourced)

            assert found
            assert kept is seed_estimate
            assert kept.kinetics.A.value_si == ESTIMATE.A.value_si
            # The sourced value was not kept.
            assert kept.kinetics.A.value_si != SOURCED.A.value_si

    def test_make_new_reaction_returns_the_incumbent_estimate(self):
        """The same inversion end to end through the true production entry point
        ``make_new_reaction`` (what both ``add_seed_mechanism_to_core`` and
        ``add_reaction_library_to_edge`` call). Both reactions enter the model the
        way a real run introduces them, and the second reaction is built from
        FRESH species objects, so the collapse is resolved through isomorphism
        de-duplication exactly as in a run -- not by sharing Python references.
        The estimate is registered first; the later sourced reaction comes back
        with ``is_new=False`` and the estimate's kinetics, so the caller silently
        drops the sourced newcomer."""
        with self._singleton():
            cerm = CoreEdgeReactionModel()

            # 1) A seed reaction carrying a previously-estimated rate enters first.
            a, b, c, d = self._species()
            seed_estimate = self._library([a, b], [c, d], SEED_LIBRARY, ESTIMATE)
            registered_seed, seed_is_new = cerm.make_new_reaction(
                seed_estimate, generate_thermo=False, generate_kinetics=False)
            assert seed_is_new is True

            # 2) The sourced library reaction enters later, over independently
            #    constructed species with the same structures.
            a2, b2, c2, d2 = self._species()
            sourced = self._library([a2, b2], [c2, d2], SOURCED_LIBRARY, SOURCED)
            returned, is_new = cerm.make_new_reaction(
                sourced, generate_thermo=False, generate_kinetics=False)

            assert is_new is False
            assert returned is registered_seed
            assert returned.kinetics.A.value_si == ESTIMATE.A.value_si
            assert returned.kinetics.A.value_si != SOURCED.A.value_si
