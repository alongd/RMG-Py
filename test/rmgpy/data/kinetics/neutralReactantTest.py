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
A reaction family could not say "my reactants are neutral molecules".

``allowChargedSpecies`` is one flag filtering ``reaction.reactants + reaction.products``
together, so a family that must make a charged product -- electron attachment, ``A + e- => A-``,
is the plain case -- was forced to also accept a charged reactant, and would happily go on to
generate ``A- + e- => A(2-)``. Nothing in the group language can close that: RMG groups match
subgraphs, so a ``c0`` constrains the atom it is written on and never the molecule, and
``generatedSpeciesConstraints`` has no charge key.

``allowChargedReactants`` is the one-sided companion. Undeclared it inherits
``allowChargedSpecies``, so it is invisible to every family that does not ask for it.

These tests go through ``KineticsFamily.generate_reactions`` -- family generation, the same call
the model builder makes -- rather than constructing a ``Reaction`` by hand, because the claim
under test is about what the engine *generates*.

Two fixture families under ``neutral_reactant_data/`` carry the argument. They are identical
except for the one declaration, and both have a deliberately wide root, ``1 *1 O u1 p2 c0``, that
matches the radical oxygen of ``[O-][O]`` as well as either oxygen of ``[O][O]``. The width is
not incidental: the committed ``Plasma_Electron_Attachment`` family generates no
anion re-attachment today, but only because its narrow root union happens to exclude ``O2-``.
A test written against that family would pass on the group shape and would keep passing if the
engine control were deleted. These fixtures make the declaration the only thing standing between
the anion and the reaction.
"""

import logging
import os.path

import pytest

from rmgpy import settings
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.molecule import Molecule

FIXTURE_FAMILIES = os.path.join(
    os.path.dirname(__file__), "neutral_reactant_data", "kinetics", "families"
)

# [O][O] -- neutral, both oxygens are radicals, net charge 0.
NEUTRAL_O2 = """
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""

# [O-][O] -- already an anion, net charge -1, and its second oxygen is an unremarkable
# `O u1 p2 c0` that the wide root matches.
ANION_O2 = """
1 O u0 p3 c-1 {2,S}
2 O u1 p2 c0  {1,S}
"""

# The same anion with the two atoms written in the other order. Nothing about the refusal may
# depend on which atom the adjacency list happens to list first.
ANION_O2_REORDERED = """
1 O u1 p2 c0  {2,S}
2 O u0 p3 c-1 {1,S}
"""


def _load(*labels):
    """
    Load the named fixture families through the production loader and return them in order.
    """
    database = KineticsDatabase()
    database.load_families(path=FIXTURE_FAMILIES, families=list(labels))
    return [database.families[label] for label in labels]


class TestNeutralReactantDeclaration:
    """
    What the declaration is, before what it does: it must be read off the family, and its default
    must be the old behaviour.
    """

    @classmethod
    def setup_class(cls):
        cls.neutral_only, cls.agnostic = _load(
            "Neutral_Only_Attachment", "Charge_Agnostic_Attachment"
        )

    def test_declaration_is_read_from_the_family(self):
        """A family that declares the control gets it, independently of allowChargedSpecies."""
        assert self.neutral_only.allow_charged_species is True
        assert self.neutral_only.allow_charged_reactants is False

    def test_undeclared_control_inherits_allow_charged_species(self):
        """
        Silence must mean "as before". Every family in the database is silent about this, so a
        default of anything other than allow_charged_species would change all of them.
        """
        assert self.agnostic.allow_charged_species is True
        assert self.agnostic.allow_charged_reactants is True

    def test_predicate_answers_for_the_declared_reason(self):
        """
        The refusal is available as a cause, not only as an absence. This is the difference
        between "the family declines this reactant" and "no reaction came back".
        """
        anion = Molecule().from_adjacency_list(ANION_O2)
        neutral = Molecule().from_adjacency_list(NEUTRAL_O2)

        assert self.neutral_only.is_charged_reactant_forbidden(anion) is True
        assert self.neutral_only.is_charged_reactant_forbidden(neutral) is False
        # The control is off here, so the same anion is fine.
        assert self.agnostic.is_charged_reactant_forbidden(anion) is False

    def test_declaration_survives_a_save_load_round_trip(self, tmp_path):
        """
        A declaration the writer drops is a declaration that vanishes the first time a family is
        regenerated. Equally, a family that never declared it must not acquire one.
        """
        armed = tmp_path / "armed_groups.py"
        self.neutral_only.save_groups(str(armed))
        assert "allowChargedReactants = False" in armed.read_text()

        silent = tmp_path / "silent_groups.py"
        self.agnostic.save_groups(str(silent))
        # The assignment, not the word: this fixture's own longDesc discusses the declaration.
        assert "allowChargedReactants = " not in silent.read_text()


class TestNeutralReactantGeneration:
    """
    The behaviour, asserted on the generated reaction list.
    """

    @classmethod
    def setup_class(cls):
        cls.neutral_only, cls.agnostic = _load(
            "Neutral_Only_Attachment", "Charge_Agnostic_Attachment"
        )
        cls.root = cls.neutral_only.groups.entries["Attacher"].item

    def test_the_anion_really_does_match_the_group(self):
        """
        The premise of every assertion below. If the anion did not match the structural root
        there would be nothing for the engine control to refuse, and an empty reaction list
        would prove only that the group shape excluded it -- the exact fragile arrangement this
        work exists to replace.
        """
        for adjacency_list in (ANION_O2, ANION_O2_REORDERED):
            anion = Molecule().from_adjacency_list(adjacency_list)
            assert anion.get_net_charge() == -1
            assert anion.is_subgraph_isomorphic(self.root), (
                "fixture root stopped matching the anion; the test below would then pass for "
                "the wrong reason"
            )

    def test_neutral_attachment_still_generates(self):
        """The intended chemistry is untouched: A + e- => A-."""
        reactions = self.neutral_only.generate_reactions([Molecule().from_adjacency_list(NEUTRAL_O2)])

        assert len(reactions) > 0
        for reaction in reactions:
            assert [spc.get_net_charge() for spc in reaction.reactants] == [0]
            assert [spc.get_net_charge() for spc in reaction.products] == [-1]

    def test_anionic_reactant_generates_nothing(self):
        """
        An empty list, returned normally. Not an exception escaping generate_reactions, which
        the model builder would have to catch and could not distinguish from a real failure.
        """
        reactions = self.neutral_only.generate_reactions([Molecule().from_adjacency_list(ANION_O2)])

        assert isinstance(reactions, list)
        assert reactions == []

    def test_refusal_does_not_depend_on_atom_ordering(self):
        """Same molecule, other adjacency-list order."""
        reactions = self.neutral_only.generate_reactions(
            [Molecule().from_adjacency_list(ANION_O2_REORDERED)]
        )

        assert reactions == []

    def test_refusal_does_not_depend_on_where_the_charge_sits(self):
        """
        A second anion, matching the root through a structurally different neighbourhood. In
        `[O]C=C[O-]` the matched atom is the same `O u1 p2 c0` as in the neutral case and the
        negative charge is three bonds away from it, so nothing local to the matched group
        distinguishes this reactant -- only the charge of the molecule as a whole, which is the
        property no group can see and the declaration exists to reach.
        """
        anion = Molecule().from_smiles("[O]C=C[O-]")
        assert anion.get_net_charge() == -1
        assert anion.is_subgraph_isomorphic(self.root)

        assert self.neutral_only.generate_reactions([anion]) == []
        # And the control family, differing by the one declaration, still generates it.
        assert len(self.agnostic.generate_reactions([Molecule().from_smiles("[O]C=C[O-]")])) > 0

    def test_refusal_holds_for_every_resonance_form(self):
        """
        The refusal reads net charge, which every resonance form of a molecule shares, so no
        form can slip through.

        Stated rather than hidden: RMG returns a single structure for both anions used here --
        `generate_resonance_structures` does not explore charge-separated forms for these -- so
        this loop currently has one member each. It is the assertion that would catch a
        regression if that changes, not evidence that it has been exercised broadly.
        """
        for adjacency_list in (ANION_O2, ANION_O2_REORDERED):
            for structure in Molecule().from_adjacency_list(adjacency_list).generate_resonance_structures():
                assert structure.get_net_charge() == -1
                assert self.neutral_only.generate_reactions([structure]) == []

    def test_refusal_is_attributed_to_the_charge_rule(self, caplog):
        """
        A dropped reaction and a declined reactant look identical from the outside -- both are an
        absent reaction. The family says which one happened, naming itself, the reactant and the
        charge, so a silent drop elsewhere in generation is not mistaken for this.
        """
        with caplog.at_level(logging.DEBUG, logger="root"):
            self.neutral_only.generate_reactions([Molecule().from_adjacency_list(ANION_O2)])

        attributed = [
            record.getMessage()
            for record in caplog.records
            if "allowChargedReactants" in record.getMessage()
        ]
        assert attributed, "the family refused the reactant without saying why"
        assert "Neutral_Only_Attachment" in attributed[0]
        assert "-1" in attributed[0]

    def test_control_is_what_refuses_it(self):
        """
        The control family and this one differ in one declaration and nothing else. This one
        still generates the dianion, so the declaration -- not the group, the recipe or the
        anion -- is what refused it next door.
        """
        reactions = self.agnostic.generate_reactions([Molecule().from_adjacency_list(ANION_O2)])

        assert len(reactions) > 0
        for reaction in reactions:
            assert [spc.get_net_charge() for spc in reaction.reactants] == [-1]
            assert [spc.get_net_charge() for spc in reaction.products] == [-2]


@pytest.mark.database
class TestChargedReactantFamiliesUnaffected:
    """
    Families that legitimately consume an ion must be untouched. They are silent about the new
    declaration, so this is the default under load from the real database rather than a fixture.
    """

    CHARGE_CONSUMING = [
        "Cation_R_Recombination",
        "Cation_Addition_MultipleBond",
        "Surface_Proton_Electron_Reduction_Alpha",
        "Surface_Proton_Electron_Reduction_Beta",
    ]

    @classmethod
    def setup_class(cls):
        cls.families_path = os.path.join(settings["database.directory"], "kinetics", "families")
        database = KineticsDatabase()
        database.load_families(
            path=cls.families_path,
            families=cls.CHARGE_CONSUMING + ["Plasma_Electron_Attachment"],
        )
        cls.families = database.families

    @pytest.mark.parametrize("label", CHARGE_CONSUMING)
    def test_charge_consuming_family_still_accepts_charged_reactants(self, label):
        family = self.families[label]

        assert family.allow_charged_species is True
        assert family.allow_charged_reactants is True
        assert family.is_charged_reactant_forbidden(Molecule().from_adjacency_list(ANION_O2)) is False

    def test_shipped_attachment_family_is_not_armed_by_this_change(self):
        """
        The engine control ships unarmed. `Plasma_Electron_Attachment` declares nothing, so it
        inherits allowChargedSpecies exactly as before; arming it is a one-line declaration in
        the database and is deliberately not part of this change.
        """
        family = self.families["Plasma_Electron_Attachment"]

        assert family.allow_charged_reactants is True

    def test_shipped_attachment_family_still_generates_neutral_attachment(self):
        """The chemistry the family exists for, from the real database, unchanged."""
        family = self.families["Plasma_Electron_Attachment"]
        reactions = family.generate_reactions([Molecule().from_adjacency_list(NEUTRAL_O2)])

        assert len(reactions) > 0
        for reaction in reactions:
            assert [spc.get_net_charge() for spc in reaction.reactants] == [0]
            assert [spc.get_net_charge() for spc in reaction.products] == [-1]

    def test_shipped_attachment_family_generates_no_dianion(self):
        """
        Zero anion re-attachment from the real family -- but read the reason, not the number.
        Today the anion does not even match the narrow root union, so this passes on group shape
        and would keep passing with the engine control deleted. That is the whole defect the
        control exists to make impossible to depend on, and it is why the fixture families above
        carry the real argument. When the root is widened, this assertion starts depending on the
        declaration, and the declaration has to be there.
        """
        family = self.families["Plasma_Electron_Attachment"]
        anion = Molecule().from_adjacency_list(ANION_O2)

        assert family.generate_reactions([anion]) == []
        assert not anion.is_subgraph_isomorphic(family.groups.entries["O_in_O2"].item)
