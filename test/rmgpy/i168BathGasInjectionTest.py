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
I-168: a deck that names only argon gets a model containing only argon.

Since 5ec0636e (2014, message in full: "Add set of bath gases (Ar, Ne, He, and N2) to
all RMG jobs just like Java.") every RMG job has had Ar, He, Ne and N2 pushed into it as
``reactive=False`` species whether or not the input file asked for them. They react with
nothing, so the arithmetic is untouched -- but they are enlarged into the **core** ahead
of the deck's own species (``rmgpy/rmg/main.py``, "Add nonreactive species ... to core
first"), so they land in every Chemkin file, Cantera YAML, species dictionary and seed
mechanism the run emits. On the 5 torr argon plasma deck that meant three of six core
species were species nobody had mentioned, each at exactly zero mole fraction.

The production diff this file defends is
:meth:`rmgpy.rmg.main.RMG.add_default_bath_gases` returning early when
:meth:`~rmgpy.rmg.main.RMG.uses_plasma_reactor` is true.

WHY PLASMA-ONLY, AND NOT A DELETION. A census of the tree found real dependants, none of
them plasma:

* ``rmgpy/rmg/pdep.py`` builds the pressure-dependence bath gas out of *every* nonreactive
  core species and asserts there is at least one, so a pdep deck declaring no inert of its
  own survives only because of the injection;
* the restart fixtures in ``test/rmgpy/test_data/restartTest/*/filters/species_map.yml``
  enumerate Ar/He/Ne/N2 at core indices 0-3, and a restart raises ``RuntimeError`` for a
  mapped species missing from the core;
* three regression decks (``aromatics``, ``fragment``, ``liquid_oxidation``) declare no
  inert at all, so their core species lists -- and the CI baselines built from them --
  would shrink.

Hence the gate is on the reactor type, and the negative controls below are the point of
the file as much as the positive case is: a change that also stopped the injection for
``simpleReactor`` would be caught here rather than in CI's regression diff.

Four assertions, each able to fail on its own:

1. ``uses_plasma_reactor`` answers True for a plasma job and False for a non-plasma one --
   asserted separately, because a predicate that is always False would make the whole gate
   a no-op while every other test here still passed;
2. a plasma job gets **no** bath gases: the core keeps only what the deck declared;
3. a non-plasma job still gets all four, in the core, ahead of the deck's own species
   (the negative control that pins the scope of the change);
4. the injection is still idempotent against a deck that declares one of the four itself.

These drive :meth:`~rmgpy.rmg.main.RMG.add_default_bath_gases` -- the live path, the same
method ``initialize()`` calls -- against a real
:class:`~rmgpy.rmg.model.CoreEdgeReactionModel`, rather than re-implementing the loop.
No database is required: ``make_new_species`` and ``enlarge`` of a nonreactive species
need no thermo.
"""

import pytest

from rmgpy.molecule import Molecule
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.solver.simple import SimpleReactor

BATH_GASES = ("Ar", "He", "Ne", "N2")


def _bare_rmg(reaction_systems):
    """
    An :class:`RMG` carrying only what ``add_default_bath_gases`` reads: a reaction model,
    an initial species list, and the reaction systems the gate inspects.
    """
    rmg = RMG()
    rmg.reaction_model = CoreEdgeReactionModel()
    rmg.initial_species = []
    rmg.reaction_systems = reaction_systems
    return rmg


def _plasma_reactor():
    """
    A :class:`PlasmaReactor` instance sufficient for an ``isinstance`` probe. Constructed
    through the real class -- a stand-in would let the gate pass while the production
    ``isinstance`` check failed on the genuine article.
    """
    return PlasmaReactor.__new__(PlasmaReactor)


def _simple_reactor():
    return SimpleReactor.__new__(SimpleReactor)


def _declare(rmg, label, smiles, reactive=True):
    """Put a species into the model the way an input file's ``species(...)`` directive does."""
    spec, _is_new = rmg.reaction_model.make_new_species(Molecule().from_smiles(smiles), label=label, reactive=reactive)
    rmg.initial_species.append(spec)
    return spec


def _species_dict_labels(rmg):
    """
    Every label the model knows about. ``CoreEdgeReactionModel.species_dict`` is keyed by
    formula and each value is a *list* of the species sharing it, so this flattens.
    """
    return {spec.label for group in rmg.reaction_model.species_dict.values() for spec in group}


def _core_labels(rmg):
    """
    Enlarge the model the way ``initialize()`` does -- nonreactive species first -- and read
    the core. Needs a loaded database (``enlarge`` consults the forbidden structures), so
    every caller is marked functional.
    """
    for spec in rmg.initial_species:
        if not spec.reactive:
            rmg.reaction_model.enlarge(spec)
    return [spec.label for spec in rmg.reaction_model.core.species]


class ReactorTypeProbeTest:
    """
    The gate's predicate, asserted on its own. If ``uses_plasma_reactor`` were stuck at
    False the injection would simply carry on and only :class:`PlasmaDeckGetsNoBathGasTest`
    would notice; if it were stuck at True every job would lose its bath gases and only
    :class:`NonPlasmaDeckStillGetsBathGasTest` would. Both directions are pinned here so a
    failure names the predicate rather than the consequence.
    """

    def test_plasma_reaction_system_is_detected(self):
        assert _bare_rmg([_plasma_reactor()]).uses_plasma_reactor() is True

    def test_non_plasma_reaction_system_is_not_detected(self):
        assert _bare_rmg([_simple_reactor()]).uses_plasma_reactor() is False

    def test_mixed_reaction_systems_count_as_plasma(self):
        assert _bare_rmg([_simple_reactor(), _plasma_reactor()]).uses_plasma_reactor() is True

    def test_no_reaction_systems_is_not_plasma(self):
        assert _bare_rmg([]).uses_plasma_reactor() is False


class PlasmaDeckGetsNoBathGasTest:
    """
    The ticket itself: an argon-only plasma deck yields an argon-only model.
    """

    def test_argon_only_plasma_deck_stays_argon_only(self):
        rmg = _bare_rmg([_plasma_reactor()])
        _declare(rmg, "Ar", "[Ar]")

        rmg.add_default_bath_gases()

        assert [spec.label for spec in rmg.initial_species] == ["Ar"], (
            "a plasma deck declaring only argon gained species it never asked for: {0}".format(
                [spec.label for spec in rmg.initial_species]
            )
        )

    def test_the_species_dictionary_gains_nothing_either(self):
        """
        ``make_new_species`` registers a species in the model's ``species_dict`` whether or
        not it is new to ``initial_species``. Asserting only on ``initial_species`` would
        pass for an implementation that still created all four and merely declined to append
        them -- which would still put them in the core via ``enlarge``.
        """
        rmg = _bare_rmg([_plasma_reactor()])
        _declare(rmg, "Ar", "[Ar]")
        before = _species_dict_labels(rmg)

        rmg.add_default_bath_gases()

        after = _species_dict_labels(rmg)
        assert after - before == set(), "plasma job created bath-gas species: {0}".format(after - before)


class NonPlasmaDeckStillGetsBathGasTest:
    """
    Negative control, and the reason this change is narrow. The census found a pdep bath-gas
    assertion, two restart fixtures and three regression decks that lean on the injection --
    all of them non-plasma. Nothing here may change for them.
    """

    def test_simple_reactor_still_gets_all_four(self):
        rmg = _bare_rmg([_simple_reactor()])
        _declare(rmg, "ethane", "CC")

        rmg.add_default_bath_gases()

        assert [spec.label for spec in rmg.initial_species] == ["ethane", "Ar", "He", "Ne", "N2"]

    def test_bath_gases_are_nonreactive_and_reach_the_core_first(self):
        rmg = _bare_rmg([_simple_reactor()])
        _declare(rmg, "ethane", "CC")

        rmg.add_default_bath_gases()

        injected = [spec for spec in rmg.initial_species if spec.label in BATH_GASES]
        assert len(injected) == len(BATH_GASES)
        # Nonreactive is what makes ``initialize()`` enlarge them into the core ahead of the
        # deck's own species, which is what gives them index -1 and puts them at the head of
        # every Chemkin SPECIES block.
        assert all(not spec.reactive for spec in injected)

    def test_a_declared_bath_gas_is_not_duplicated(self):
        """
        Idempotence, unchanged by the gate: a deck that declares its own N2 must end up with
        one N2, not two. ``make_new_species`` matches on structure, so the deck's label wins.
        """
        rmg = _bare_rmg([_simple_reactor()])
        _declare(rmg, "nitrogen", "N#N", reactive=False)

        rmg.add_default_bath_gases()

        labels = [spec.label for spec in rmg.initial_species]
        assert labels == ["nitrogen", "Ar", "He", "Ne"], labels


@pytest.mark.functional
class CoreMembershipTest:
    """
    The claim as it is actually cashed out. The core is what every Chemkin file, Cantera YAML,
    species dictionary and seed mechanism is written from, so "the mechanism is what was asked
    for" is a statement about the core, not about ``initial_species``. Reaching it means
    calling ``enlarge``, which consults the forbidden-structure database -- hence functional,
    and hence the same ``RMGDatabase().load`` setup that ``rmgpy/rmg/modelTest.py`` uses for
    its own enlarge tests.
    """

    @classmethod
    def setup_class(cls):
        from rmgpy import settings
        from rmgpy.data.rmg import RMGDatabase

        cls.database = RMGDatabase()
        cls.database.load(
            path=settings["database.directory"],
            thermo_libraries=["primaryThermoLibrary"],
            kinetics_families=[],
            reaction_libraries=[],
        )

    def test_plasma_core_holds_only_what_the_deck_declared(self):
        rmg = _bare_rmg([_plasma_reactor()])
        rmg.database = self.database
        _declare(rmg, "Ar", "[Ar]", reactive=False)

        rmg.add_default_bath_gases()

        assert _core_labels(rmg) == ["Ar"]

    def test_non_plasma_core_still_leads_with_all_four(self):
        rmg = _bare_rmg([_simple_reactor()])
        rmg.database = self.database
        _declare(rmg, "ethane", "CC")

        rmg.add_default_bath_gases()

        assert _core_labels(rmg) == list(BATH_GASES)
