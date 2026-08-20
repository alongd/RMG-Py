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
RMG can already *generate* electron attachment from the real
``Plasma_Electron_Attachment`` family, but until now nothing in the suite proved it. The
family exporting the right classes, and the guard refusing charged reactants, were each pinned
elsewhere -- but the one claim the milestone is actually about, "start from neutral O2 and the
engine produces ``O2 + e- -> O2-`` with the right bookkeeping", was only ever measured by hand.
This file commits that measurement.

Everything here runs through the ordinary family-generation path -- ``KineticsFamily``'s own
``generate_reactions`` and ``KineticsDatabase.generate_reactions_from_families``, the same calls
the model builder makes -- against the real database on ``settings['database.directory']``, not a
fixture. The point of the milestone is that the *shipped* family does this, so a fixture would
prove the wrong thing.

Provenance, not a coincidence of names. The reaction that comes back is a ``TemplateReaction``
carrying ``family = 'Plasma_Electron_Attachment'``, produced by family generation, which never
consults a kinetics library. That it is *not* a ``LibraryReaction`` is asserted directly: a
library lookup and a generated reaction would both read as "an O2- reaction exists", and only the
type and the generation path tell them apart.

The degeneracy question is pinned at BOTH levels on purpose. The family level returns two
reactions that are *identical* to each other; the database level returns one. A test that checked
only the database level would not notice the family-level pair reappearing as a genuine second
channel (the rate silently doubling); a test that checked only the family level would not notice
the collapse breaking. The collapse is a *duplicate drop*, not a sum that happens to reach 1.0:
in ``rmgpy.data.kinetics.common.find_degenerate_reactions`` the second reaction is found identical
to the first (``check_identical=True``) and is discarded without ever joining a sublist, so the
surviving one-member sublist sums to degeneracy 1.0. Whether 1.0 is the physically correct value
for two symmetry-equivalent oxygens is a scientific question this file does not settle: it asserts
what the code produces and leaves the question open.

One assertion is characterization, and says so at the assertion: the generated anion ``[O][O-]``
re-attaching through the same family returns zero reactions. Today that is true only because the
family's narrow root union happens to exclude ``O2-`` -- it documents the current group shape and
would keep passing if the engine's charged-reactant control were deleted. The real test of that
control lives in ``neutralReactantTest.py`` against purpose-built fixture families, and must not be
confused with this one.
"""

import os.path

import pytest

from rmgpy import settings
from rmgpy.data.kinetics.common import get_molecularity
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.molecule import Molecule
from rmgpy.species import Species

# [O][O] -- neutral molecular oxygen, both oxygens radicals, net charge 0.
NEUTRAL_O2 = "[O][O]"
# [O][O-] -- the superoxide anion the family produces, net charge -1.
ANION_O2 = "[O][O-]"

FAMILY = "Plasma_Electron_Attachment"


@pytest.mark.database
class TestPlasmaAttachmentGeneration:
    """
    Generation of ``O2 + e- -> O2-`` from the real ``Plasma_Electron_Attachment`` family.
    """

    @classmethod
    def setup_class(cls):
        """
        Load the real family (and the recommended sets, so a caller's view of it is complete) once,
        then generate the attachment reactions both ways.
        """
        cls.families_path = os.path.join(settings["database.directory"], "kinetics", "families")
        database = KineticsDatabase()
        database.load_recommended_families(os.path.join(cls.families_path, "recommended.py"))
        database.load_families(cls.families_path, families=[FAMILY])
        cls.database = database
        cls.family = database.families[FAMILY]

        # Family-level generation: the direct KineticsFamily path.
        cls.family_reactions = cls.family.generate_reactions(
            [Molecule().from_smiles(NEUTRAL_O2)]
        )

        # Database-level generation: the same call the model builder makes, on a Species whose
        # resonance structures have been generated. only_families keeps the list to this family so
        # the count is unambiguous.
        species = Species().from_smiles(NEUTRAL_O2)
        species.generate_resonance_structures()
        cls.database_reactions = cls.database.generate_reactions_from_families(
            [species], only_families=[FAMILY]
        )

    # -- item 1: the reaction is generated, from this family, and is not a library lookup --------

    def test_attachment_reaction_is_generated(self):
        """A reaction comes back at all; generation did not silently return nothing."""
        assert len(self.family_reactions) > 0
        assert len(self.database_reactions) == 1

    def test_reaction_is_o2_plus_electron_to_superoxide(self):
        """
        The generated reaction really is O2 -> O2- (electron implicit), not some other channel.
        Asserted by structural isomorphism against reference molecules, not by SMILES-string
        equality, so it does not turn on the SMILES backend's canonical ordering.
        """
        reactant_ref = Molecule().from_smiles(NEUTRAL_O2)
        product_ref = Molecule().from_smiles(ANION_O2)
        reaction = self.database_reactions[0]
        assert len(reaction.reactants) == 1 and len(reaction.products) == 1
        assert reaction.reactants[0].is_isomorphic(reactant_ref)
        assert reaction.products[0].is_isomorphic(product_ref)

    def test_provenance_is_the_family_not_a_library(self):
        """
        The reaction is a family-generated ``TemplateReaction`` from this family, not a
        ``LibraryReaction``. The proof is the *generation path*, not any single attribute: both
        ``KineticsFamily.generate_reactions`` and ``generate_reactions_from_families`` go through
        families only and never search a kinetics library (that happens only in
        ``KineticsDatabase.generate_reactions`` when no families are named). The type and family
        label are then what a later reader checks to confirm it. Both family-level reactions are
        checked, not only the first, so a mislabelled duplicate cannot slip through.
        """
        for reaction in (self.family_reactions + self.database_reactions):
            assert isinstance(reaction, TemplateReaction)
            assert not isinstance(reaction, LibraryReaction)
            assert reaction.family == FAMILY
            assert reaction.template == ["O_in_O2"]
            # get_source() on a TemplateReaction returns its family label (it is not independent of
            # .family, but it is the accessor the model builder itself reads for provenance).
            assert reaction.get_source() == FAMILY

    # -- item 2: electron bookkeeping ------------------------------------------------------------

    def test_reaction_consumes_one_electron(self):
        """The family carries ``electrons = -1``; every generated reaction inherits it."""
        for reaction in (self.family_reactions + self.database_reactions):
            assert reaction.electrons == -1

    # -- item 3: molecularity --------------------------------------------------------------------

    def test_molecularity_is_correct(self):
        """
        The reaction is written with one heavy reactant and one heavy product, the electron left
        implicit (``electrons = -1``). Counting the implicit electron the reaction is bimolecular
        -- A + e- -- which is what fixes the rate-coefficient units below. Molecularity is taken
        from the production helper ``get_molecularity`` (``len(reactants) + max(-electrons, 0)``),
        not a hand-rolled count: a released electron (``electrons > 0``) must NOT add to it, and
        ``max(-electrons, 0)`` is the convention that gets that right where ``abs`` would not.
        """
        reaction = self.database_reactions[0]
        assert len(reaction.reactants) == 1
        assert len(reaction.products) == 1
        assert get_molecularity(reaction) == 2

    # -- item 4: rate-coefficient units ----------------------------------------------------------

    def test_rate_coefficient_units_match_the_molecularity(self):
        """
        The rate the family would give the generated reaction carries second-order
        ``cm^3/(mol*s)`` units, matching the bimolecular A + e- molecularity: a volume per mole per
        time, the A-factor dimensions of an order-2 rate law. First-order ``s^-1`` units would
        signal the electron had been dropped from the rate law.

        The template retrieved is the generated reaction's *own* ``['O_in_O2']`` template, so the
        units are tied to that reaction, not to an arbitrary rule. This checks the units the
        molecularity implies; it is not itself a second proof of the electron bookkeeping, which
        the molecularity test above pins directly.

        A fresh family is loaded here on purpose: ``add_rules_from_training`` and
        ``fill_rules_by_averaging_up`` mutate the family's rule table, and doing that to the
        class-shared ``self.family`` would leave an order-dependent side effect for any later test.
        """
        database = KineticsDatabase()
        database.load_families(self.families_path, families=[FAMILY])
        family = database.families[FAMILY]
        family.add_rules_from_training(thermo_database=None)
        family.fill_rules_by_averaging_up()

        template = family.retrieve_template(self.database_reactions[0].template)
        kinetics = family.get_kinetics_for_template(template, degeneracy=1)[0]

        assert hasattr(kinetics, "A"), f"kinetics {type(kinetics).__name__} has no rate coefficient A"
        assert kinetics.A.units == "cm^3/(mol*s)"
        # The units above are the second-order form the bimolecular molecularity demands.
        assert get_molecularity(self.database_reactions[0]) == 2

    # -- item 5: charge behaviour ----------------------------------------------------------------

    def test_reactant_is_neutral_and_product_is_anionic(self):
        """
        Neutral O2 in, superoxide anion out: net charge 0 -> -1. Every generated reaction is
        checked -- both family-level duplicates and the database-level one -- so a duplicate that
        collapsed the charge wrong could not hide behind its twin.
        """
        for reaction in (self.family_reactions + self.database_reactions):
            assert [self._net_charge(spc) for spc in reaction.reactants] == [0]
            assert [self._net_charge(spc) for spc in reaction.products] == [-1]

    # -- item 6: charged-reactant rejection is armed on this family ------------------------------

    def test_charged_reactant_rejection_is_active(self):
        """
        The guard the milestone relies on is declared on the shipped family: it accepts charged
        *species* in general but refuses charged *reactants*. Silence would have made these equal;
        the asymmetry is what makes it a declaration.
        """
        assert self.family.allow_charged_species is True
        assert self.family.allow_charged_reactants is False

    # -- item 7: the anion does not re-attach (CHARACTERIZATION, not a guard test) ----------------

    def test_generated_anion_does_not_reattach(self):
        """
        CHARACTERIZATION, not a test of the charged-reactant control. The generated anion
        ``[O][O-]`` returns zero reactions through this family, but *today* that holds only because
        the family's narrow root union excludes ``O2-`` -- the anion does not even match the
        ``O_in_O2`` group. This assertion documents the current root-group shape and would keep
        passing if ``allow_charged_reactants`` were deleted. The real test of that engine control
        lives in ``neutralReactantTest.py`` (``TestNeutralReactantGeneration`` and
        ``TestChargedReactantFamiliesUnaffected``) against purpose-built fixture families with a
        deliberately wide root; that coverage is STRONGER than anything the real family can offer
        and must not be deleted as redundant to this. Do not mistake this for it, and do not widen
        the root here to make it mean more.
        """
        anion = Molecule().from_smiles(ANION_O2)
        assert anion.get_net_charge() == -1
        assert self.family.generate_reactions([anion]) == []
        # The reason it is empty: the anion does not match the root group at all. If this stops
        # being true (the root is widened), the assertion above starts depending on the engine
        # control, and this comment stops applying.
        assert not anion.is_subgraph_isomorphic(self.family.groups.entries["O_in_O2"].item)

    # -- the degeneracy question: pin BOTH levels ------------------------------------------------

    def test_family_level_returns_two_identical_reactions(self):
        """
        The direct family path returns two reactions that are *identical* to one another -- not
        merely isomorphic. Both carry degeneracy 1.0. Pinning the pair here is what would catch the
        duplicate becoming a genuine second channel that doubles the attachment rate.
        """
        assert len(self.family_reactions) == 2
        first, second = self.family_reactions
        assert first.is_isomorphic(second)
        assert first.is_isomorphic(
            second, check_identical=True, strict=False, check_template_rxn_products=True
        )
        assert first.degeneracy == 1.0
        assert second.degeneracy == 1.0

    def test_database_level_collapses_to_one(self):
        """
        The database path collapses the identical pair to a single reaction of degeneracy 1.0: the
        attachment rate is not doubled. Pinning this level is what would catch the collapse
        breaking. The value 1.0 -- not 2.0 -- is itself what distinguishes a *duplicate drop* from
        a degeneracy *sum*: had ``find_degenerate_reactions`` summed the two degeneracy-1.0
        reactions instead of discarding the identical copy, this would read 2.0.
        """
        assert len(self.database_reactions) == 1
        assert self.database_reactions[0].degeneracy == 1.0
        # A sum of the two family-level reactions (1.0 + 1.0) would be 2.0; the drop keeps it at 1.0.
        assert self.database_reactions[0].degeneracy != 2.0

    # -- the family is not reachable without being named -----------------------------------------

    def test_family_is_a_member_of_no_recommended_set(self):
        """
        ``Plasma_Electron_Attachment`` is a member of no recommended set -- not ``default``, and
        there is no ``plasma`` set at all -- so a caller reaches it only by naming it. This asserts
        set *membership* (the shipped ``recommended.py`` is a flat module of label sets); it is not
        a claim about every reachability path. Reported, not fixed: adding a recommended set is out
        of scope for this ticket.
        """
        recommended = self.database.recommended_families
        assert "plasma" not in recommended
        for name, members in recommended.items():
            assert FAMILY not in members, f"unexpectedly found {FAMILY} in recommended set {name!r}"

    @staticmethod
    def _net_charge(species_or_molecule):
        """Net charge whether generation handed back a Molecule (family level) or Species (db)."""
        if isinstance(species_or_molecule, Molecule):
            return species_or_molecule.get_net_charge()
        return species_or_molecule.molecule[0].get_net_charge()
