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
I-148: the electron-placement declaration must survive the seed round trip.

The placement declaration is keyed on a reaction's owner label
(``FAMILY_ELECTRON_PLACEMENT``). The seed mechanism renames its container to the
fixed label ``seed`` on write and ``restart`` on read, so before I-148 every
plasma reaction reloaded from a seed lost its declaration and fell back to the
net-derived rule. The original owner is NOT lost, though: the seed writer emits
``Originally from reaction library: <label>`` into each entry's ``longDesc``,
and ``KineticsLibrary.get_library_reactions`` parses it back into
``LibraryReaction.library`` for auto-generated libraries -- while overwriting
``family`` with the container label. I-148 makes the placement-owner lookup read
that preserved per-reaction provenance (``family`` first, then ``library``)
instead of the renamed container alone, and adds a loud seed-write-time warning
for any reaction whose declaration would NOT survive the round trip.

These tests pin both halves on constructed objects, with no database loaded.
The measured asymmetry they encode, from ``docs/i123d-audit.md`` section 9:

    channel                      canonical    via seed (pre-I-148)
    electron-impact ionisation   (1, 2)       (0, 1)   declaration lost
    radiative recombination      (1, 0)       (1, 0)   preserved by coincidence
"""

import logging

import pytest

from rmgpy.data.base import Entry
from rmgpy.data.kinetics.library import KineticsLibrary, LibraryReaction, \
    seed_placement_survives
from rmgpy.electron_balance import expand_electrons, get_electron_placement_counts, \
    get_placement_owner
from rmgpy.molecule import Molecule
from rmgpy.rmg.main import warn_if_seed_loses_placement
from rmgpy.species import Species

IONISATION = 'PlasmaElectronImpactIonization'
RECOMBINATION = 'PlasmaRadiativeRecombination'
ATTACHMENT = 'Plasma_Electron_Attachment'


def electron():
    return Species(label='e-', molecule=[Molecule().from_adjacency_list('1 e u1 p0 c-1')])


def lithium():
    return Species(label='Li', molecule=[
        Molecule().from_adjacency_list('multiplicity 2\n1 Li u1 p0 c0')])


def lithium_cation():
    return Species(label='Lip', molecule=[
        Molecule().from_adjacency_list('1 Li u0 p0 c+1')])


def seed_shaped(reactants, products, electrons, original_library, container='seed'):
    """
    A LibraryReaction exactly as ``KineticsLibrary.get_library_reactions``
    rebuilds it from an auto-generated (seed) library: ``library`` carries the
    parsed original owner, ``family`` is then overwritten with the container
    label (``seed`` on a plain reload, ``restart`` on a restart).
    """
    rxn = LibraryReaction(reactants=reactants, products=products,
                          electrons=electrons, library=original_library)
    rxn.family = container
    return rxn


class OwnerRecoveryTest:
    """The placement owner is read from the reaction, not from its container."""

    def test_seed_shaped_ionisation_recovers_its_declaration(self):
        rxn = seed_shaped([lithium()], [lithium_cation()], +1, IONISATION)
        assert get_placement_owner(rxn) == IONISATION
        assert get_electron_placement_counts(rxn) == (1, 2)

    def test_restart_shaped_ionisation_recovers_its_declaration(self):
        rxn = seed_shaped([lithium()], [lithium_cation()], +1, IONISATION,
                          container='restart')
        assert get_placement_owner(rxn) == IONISATION
        assert get_electron_placement_counts(rxn) == (1, 2)

    def test_seed_shaped_recombination_recovers_its_declaration(self):
        rxn = seed_shaped([lithium_cation()], [lithium()], -1, RECOMBINATION)
        assert get_placement_owner(rxn) == RECOMBINATION
        assert get_electron_placement_counts(rxn) == (1, 0)

    def test_family_wins_over_library_provenance(self):
        # The current attribution governs; the parsed provenance is consulted
        # only when the current attribution carries no declaration.
        rxn = seed_shaped([lithium_cation()], [lithium()], -1, IONISATION,
                          container=RECOMBINATION)
        assert get_placement_owner(rxn) == RECOMBINATION
        assert get_electron_placement_counts(rxn) == (1, 0)

    def test_undeclared_owner_still_uses_the_net_rule(self):
        rxn = seed_shaped([lithium_cation()], [lithium()], -1, 'SomeNeutralLibrary')
        assert get_placement_owner(rxn) is None
        assert get_electron_placement_counts(rxn) == (1, 0)
        rxn2 = seed_shaped([lithium()], [lithium_cation()], +1, 'SomeNeutralLibrary')
        assert get_placement_owner(rxn2) is None
        assert get_electron_placement_counts(rxn2) == (0, 1)

    def test_no_owner_at_all_is_unchanged(self):
        rxn = LibraryReaction(reactants=[lithium_cation()], products=[lithium()],
                              electrons=-1, library=None)
        assert get_placement_owner(rxn) is None
        assert get_electron_placement_counts(rxn) == (1, 0)


class SeedIdentityTest:
    """
    The two consequences measured in the audit: the seed copy of the ionisation
    channel compared as a DIFFERENT reaction from its canonical original (so
    restart added the channel twice), while the recombination pair compared
    equal by arithmetic coincidence.
    """

    def test_seed_copy_of_ionisation_is_the_canonical_reaction(self):
        canonical = LibraryReaction(reactants=[lithium()], products=[lithium_cation()],
                                    electrons=+1, library=IONISATION)
        via_seed = seed_shaped([lithium()], [lithium_cation()], +1, IONISATION)
        assert canonical.is_isomorphic(via_seed)
        assert canonical.is_isomorphic(via_seed, either_direction=True)

    def test_seed_copy_of_recombination_is_the_canonical_reaction(self):
        canonical = LibraryReaction(reactants=[lithium_cation()], products=[lithium()],
                                    electrons=-1, library=RECOMBINATION)
        via_seed = seed_shaped([lithium_cation()], [lithium()], -1, RECOMBINATION)
        assert canonical.is_isomorphic(via_seed, either_direction=True)

    def test_the_two_channels_do_not_collapse_into_one_reaction(self):
        # The I-134 hazard must not resurface through the recovered owners:
        # ionisation reversed is (2, 1), which is not recombination's (1, 0).
        ionisation = seed_shaped([lithium()], [lithium_cation()], +1, IONISATION)
        recombination = seed_shaped([lithium_cation()], [lithium()], -1, RECOMBINATION)
        assert not ionisation.is_isomorphic(recombination, either_direction=True)


class SeedExportTest:
    """The writers place the recovered owner's declaration, not the net rule."""

    def test_expand_electrons_places_both_sides_for_seed_shaped_ionisation(self):
        e = electron()
        rxn = seed_shaped([lithium()], [lithium_cation()], +1, IONISATION)
        reactants, products = expand_electrons(rxn, [rxn.reactants[0], rxn.products[0], e])
        assert sum(1 for spc in reactants if spc.is_electron()) == 1
        assert sum(1 for spc in products if spc.is_electron()) == 2

    def test_expand_electrons_net_rule_untouched_for_undeclared_owner(self):
        e = electron()
        rxn = seed_shaped([lithium()], [lithium_cation()], +1, 'SomeNeutralLibrary')
        reactants, products = expand_electrons(rxn, [rxn.reactants[0], rxn.products[0], e])
        assert sum(1 for spc in reactants if spc.is_electron()) == 0
        assert sum(1 for spc in products if spc.is_electron()) == 1


def entry_for(rxn, long_desc):
    entry = Entry(index=1, label=str(rxn), item=rxn, data=rxn.kinetics)
    entry.long_desc = long_desc
    return entry


class SeedWriteWarningTest:
    """
    The loud warning at seed-write time: a reaction whose placement declaration
    would not survive the seed round trip must announce itself in the run that
    writes the seed, not in the run that reloads it.
    """

    def test_survives_when_the_library_provenance_line_is_present(self):
        rxn = LibraryReaction(reactants=[lithium()], products=[lithium_cation()],
                              electrons=+1, library=IONISATION)
        entry = entry_for(rxn, 'Originally from reaction library: {0}\n'.format(IONISATION))
        assert seed_placement_survives(entry, IONISATION)

    def test_survives_when_the_rate_rule_family_lines_are_present(self):
        rxn = LibraryReaction(reactants=[lithium()], products=[lithium_cation()],
                              electrons=-1, library=ATTACHMENT)
        entry = entry_for(rxn, 'Estimated using template [X] for rate rule [Y]\n'
                               'family: {0}\n'.format(ATTACHMENT))
        assert seed_placement_survives(entry, ATTACHMENT)

    def test_does_not_survive_without_parsed_provenance(self):
        rxn = LibraryReaction(reactants=[lithium()], products=[lithium_cation()],
                              electrons=+1, library=IONISATION)
        assert not seed_placement_survives(entry_for(rxn, ''), IONISATION)
        assert not seed_placement_survives(
            entry_for(rxn, 'Originally from reaction library: SomeOtherLibrary\n'),
            IONISATION)

    def test_warning_fires_for_an_unpreservable_declared_reaction(self, caplog):
        rxn = LibraryReaction(reactants=[lithium()], products=[lithium_cation()],
                              electrons=+1, library=IONISATION)
        entry = entry_for(rxn, 'a comment carrying no parsed provenance at all')
        with caplog.at_level(logging.WARNING):
            warn_if_seed_loses_placement(rxn, entry, 'core')
        assert any(IONISATION in rec.message and 'placement' in rec.message
                   for rec in caplog.records if rec.levelno >= logging.WARNING)

    def test_warning_silent_when_provenance_survives(self, caplog):
        rxn = LibraryReaction(reactants=[lithium()], products=[lithium_cation()],
                              electrons=+1, library=IONISATION)
        entry = entry_for(rxn, 'Originally from reaction library: {0}\n'.format(IONISATION))
        with caplog.at_level(logging.WARNING):
            warn_if_seed_loses_placement(rxn, entry, 'core')
        assert not caplog.records

    def test_warning_silent_for_an_undeclared_neutral_reaction(self, caplog):
        rxn = LibraryReaction(reactants=[lithium()], products=[lithium()],
                              electrons=0, library='SomeNeutralLibrary')
        entry = entry_for(rxn, '')
        with caplog.at_level(logging.WARNING):
            warn_if_seed_loses_placement(rxn, entry, 'core')
        assert not caplog.records
