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
An electron transferred by a reaction is a reactant particle on one side of the load boundary and a
piece of metadata on the other. A reaction family declares its electron count once, in ``groups.py``;
generated reactions receive it (``KineticsFamily._create_reaction`` passes ``electrons=self.electrons``
into every ``TemplateReaction``) but reactions loaded from the family's own training depository used
to receive it only when the *kinetics data object* carried one, which is true for surface charge
transfer and for nothing else. A family's declaration therefore died silently at the depository load.

Two consequences, both tested here:

1. The count must reach the loaded training reactions, so that anything reading ``rxn.electrons``
   after a load sees what the family declared.
2. Molecularity must count electrons that are consumed. ``len(reaction.reactants)`` is not the
   molecularity of a reaction whose electron is carried by a flag, which is why the units checker in
   ``test/database/databaseTest.py`` demanded first-order units from a bimolecular electron
   attachment and rejected a correct rate.

Everything here is deliberately non-plasma: methyl radical recombination and hand-built reactions
with no family at all. The defect is in how an electron count crosses a boundary, not in plasma
chemistry, and a proof that only holds for the plasma family would not be a proof of the fix.
"""

import os.path

from rmgpy.data.kinetics.common import get_molecularity
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.depository import KineticsDepository
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.thermo import ThermoData


class TestMolecularityCountsElectrons:
    """
    ``get_molecularity`` on hand-built reactions: no family, no database, no plasma.
    """

    @classmethod
    def setup_class(cls):
        cls.methyl = Species().from_smiles("[CH3]")
        cls.ethane = Species().from_smiles("CC")
        cls.oxygen = Species().from_smiles("[O][O]")

    def test_electrons_zero_is_unchanged(self):
        """The ordinary case must be untouched: molecularity is the number of reactants."""
        unimolecular = Reaction(reactants=[self.ethane], products=[self.methyl, self.methyl], electrons=0)
        bimolecular = Reaction(reactants=[self.methyl, self.methyl], products=[self.ethane], electrons=0)
        assert get_molecularity(unimolecular) == 1
        assert get_molecularity(bimolecular) == 2

    def test_consumed_electron_counts_as_a_reactant(self):
        """
        One reactant plus one consumed electron is a bimolecular reaction. This is the shape of
        non-dissociative electron attachment, and the reason its rate is in cm^3/(mol*s).
        """
        attachment = Reaction(reactants=[self.oxygen], products=[self.oxygen], electrons=-1)
        assert get_molecularity(attachment) == 2

    def test_consumed_electrons_count_one_each(self):
        """A two-electron reduction of a single species is termolecular."""
        two_electron = Reaction(reactants=[self.oxygen], products=[self.oxygen], electrons=-2)
        assert get_molecularity(two_electron) == 3

    def test_consumed_electron_adds_to_a_bimolecular_reaction(self):
        """The correction is additive, not a special case for one reactant."""
        rxn = Reaction(reactants=[self.methyl, self.methyl], products=[self.ethane], electrons=-1)
        assert get_molecularity(rxn) == 3

    def test_released_electron_does_not_count(self):
        """
        A positive count means the electron is released, i.e. it is a product. Products do not
        enter the rate expression, so an ionization is unimolecular in its single reactant.
        """
        ionization = Reaction(reactants=[self.oxygen], products=[self.oxygen], electrons=1)
        assert get_molecularity(ionization) == 1

    def test_missing_attribute_falls_back_to_the_reactant_count(self):
        """
        Anything reaction-shaped but lacking an electron count keeps the old behaviour rather than
        raising. Entries whose item is a template reaction over groups take this path.
        """

        class ReactionWithoutElectrons:
            reactants = [object(), object()]

        assert get_molecularity(ReactionWithoutElectrons()) == 2


class TestFamilyPropagatesElectronsToItsTrainingDepository:
    """
    A family's declared electron count must reach the reactions loaded from its training
    depository. Both fixture families are ordinary methyl recombination; they differ only in
    whether ``groups.py`` declares ``electrons``.
    """

    @classmethod
    def setup_class(cls):
        cls.families_root = os.path.join(os.path.dirname(__file__), "electron_propagation_data")
        database = KineticsDatabase()
        database.load_families(
            path=cls.families_root,
            families=["Electron_Carrying_Recombination", "Plain_Recombination"],
        )
        cls.families = database.families

    def _training_entry(self, family_label, index):
        """Return the loaded training entry `index` of `family_label`, checking that it loaded."""
        family = self.families[family_label]
        depository = family.get_training_depository()
        assert index in depository.entries, (
            f"{family_label} training entry {index} did not load; got {sorted(depository.entries)}"
        )
        return depository.entries[index]

    def test_family_declaration_is_read(self):
        """Guard on the fixtures themselves: the two families really do differ in what they declare."""
        assert self.families["Electron_Carrying_Recombination"].electrons == -1
        assert self.families["Plain_Recombination"].electrons == 0

    def test_declared_electrons_reach_the_loaded_training_reactions(self):
        """
        The assertion is on the loaded object. Reading the fixture's groups.py would pass while the
        propagation stayed broken, which is exactly how this defect survived.
        """
        entry = self._training_entry("Electron_Carrying_Recombination", 0)
        assert entry.item.electrons == -1, f"entry {entry.label}: electrons={entry.item.electrons}"

    def test_training_reactions_of_an_ordinary_family_are_unchanged(self):
        """The control: a family that declares nothing leaves its training reactions at zero."""
        entry = self._training_entry("Plain_Recombination", 0)
        assert entry.item.electrons == 0, f"entry {entry.label}: electrons={entry.item.electrons}"

    def test_loaded_training_reactions_have_the_physical_molecularity(self):
        """
        End to end: the value the units checker consumes. The same reaction, ``CH3 + CH3 <=> C2H6``,
        is termolecular in the family that consumes an electron and bimolecular in the one that does
        not, while ``len(reactants)`` is 2 in both.
        """
        with_electron = self._training_entry("Electron_Carrying_Recombination", 0).item
        without_electron = self._training_entry("Plain_Recombination", 0).item

        assert len(with_electron.reactants) == 2
        assert len(without_electron.reactants) == 2

        assert get_molecularity(with_electron) == 3
        assert get_molecularity(without_electron) == 2

    def test_data_borne_count_still_beats_the_family_declaration(self):
        """
        The family's count is only a default. Kinetics that carry their own electron count -- charge
        transfer, i.e. every electrochemistry entry in the real database -- must keep it, so this fix
        cannot silently rewrite them. The fixture's entry 1 declares +1 against the family's -1.
        """
        entry = self._training_entry("Electron_Carrying_Recombination", 1)
        assert entry.data.electrons.value == 1, "fixture no longer sets a conflicting data count"
        assert entry.item.electrons == 1, f"family default overwrote the data count: {entry.item.electrons}"

    def test_depository_carries_the_count_handed_down_by_its_family(self):
        """The depository is where the family's declaration lands; it has no other route to it."""
        assert KineticsDepository(label="fixture", electrons=-1).electrons == -1
        assert KineticsDepository(label="fixture").electrons == 0


class TestReverseTrainingReactionNegatesElectronCount:
    """
    A training reaction stored in the reverse (dissociation) direction is flipped to the family's
    forward direction by get_training_set. That flip rebuilt the reaction with a bare Reaction()
    that never passed ``electrons``, so ``electrons`` defaulted to 0 -- the count was *erased*, not
    merely left un-negated. The family declares ``electrons = -1``, so the flipped forward reaction
    must instead carry the negated count, +1.
    """

    @classmethod
    def setup_class(cls):
        root = os.path.join(os.path.dirname(__file__), "electron_reverse_data")
        database = KineticsDatabase()
        database.load_families(path=root, families=["Electron_Reverse_Recombination"])
        cls.family = database.families["Electron_Reverse_Recombination"]

    @staticmethod
    def _placeholder_thermo():
        # The flip computes a reverse rate coefficient, which needs species free energies; supply a
        # placeholder so the flip runs at all. Its numeric value is irrelevant to the electron
        # count, which is all this test asserts on.
        return ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
                          Cpdata=([10, 10, 10, 10, 10, 10, 10], 'cal/(mol*K)'),
                          H298=(0, 'kcal/mol'), S298=(50, 'cal/(mol*K)'),
                          Cp0=(4.0, 'cal/(mol*K)'), CpInf=(20.0, 'cal/(mol*K)'))

    def test_flipped_reverse_training_reaction_negates_electrons(self):
        depository = self.family.get_training_depository()
        for entry in depository.entries.values():
            for spc in entry.item.reactants + entry.item.products:
                spc.thermo = self._placeholder_thermo()

        rxns = self.family.get_training_set(get_reverse=True, estimate_thermo=False)

        assert len(rxns) == 1
        flipped = rxns[0]
        # The one entry was stored reversed relative to the recombination template, so what comes
        # back is the forward-oriented recombination.
        assert len(flipped.reactants) == 2 and len(flipped.products) == 1, str(flipped)
        assert flipped.electrons == 1, f"expected the negated count +1, got {flipped.electrons}"


class TestCreateReactionNegatesElectronsWhenReversed:
    """
    KineticsFamily._create_reaction swaps the participant lists when is_forward is False but copied
    the family-forward electron declaration onto the reaction unchanged. Reaction.electrons is
    signed relative to the reaction's current orientation, so a reversed reaction must carry the
    negated count.
    """

    @classmethod
    def setup_class(cls):
        from rmgpy.molecule import Molecule

        root = os.path.join(os.path.dirname(__file__), "electron_propagation_data")
        database = KineticsDatabase()
        database.load_families(path=root, families=["Electron_Carrying_Recombination"])
        cls.family = database.families["Electron_Carrying_Recombination"]
        assert cls.family.electrons == -1
        cls.reactants = [Molecule().from_smiles("[CH3]"), Molecule().from_smiles("[CH3]")]
        cls.products = [Molecule().from_smiles("CC")]

    def test_forward_reaction_keeps_the_family_declaration(self):
        forward = self.family._create_reaction(self.reactants, self.products, is_forward=True)
        assert forward is not None
        assert forward.electrons == -1

    def test_reverse_reaction_negates_the_family_declaration(self):
        reverse = self.family._create_reaction(self.reactants, self.products, is_forward=False)
        assert reverse is not None
        assert reverse.electrons == 1


class TestDepositoryDataBorneElectronPrecedence:
    """
    KineticsDepository.load kept the data-borne electron count only for an isinstance allowlist
    (SurfaceChargeTransfer, SurfaceArrheniusBEP) that missed SurfaceChargeTransferBEP,
    ArrheniusChargeTransfer and ArrheniusChargeTransferBM -- all of which carry their own electrons
    field -- so a training entry using one of those had its count overwritten by the family default.
    Keying off the presence of the attribute fixes all of them, and any future class, at once.
    """

    def test_data_borne_count_wins_for_every_charge_transfer_class(self):
        from rmgpy.data.kinetics.depository import KineticsDepository
        from rmgpy.data.kinetics.database import KineticsDatabase

        root = os.path.join(os.path.dirname(__file__), "electron_precedence_data")
        depository = KineticsDepository(electrons=-1)  # family default the entries must override
        depository.load(os.path.join(root, "reactions.py"), KineticsDatabase().local_context, {})

        by_index = {entry.index: entry for entry in depository.entries.values()}
        assert by_index[0].item.electrons == 2, "SurfaceChargeTransferBEP count overwritten"
        assert by_index[1].item.electrons == -3, "ArrheniusChargeTransfer count overwritten"
        assert by_index[2].item.electrons == 2, "ArrheniusChargeTransferBM count overwritten"
