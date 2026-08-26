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
Regression lock for I-086: the electron sign of a *reverse*-generated cation reaction.

The setting. ``Cation_R_Recombination`` is reversible and NOT its own reverse; its forward
template is ``Li+ + R.  (+ e-)  ->  R-Li`` with a family-forward declaration ``electrons = -1``
(one electron consumed on the reactant side). When RMG has only the neutral product ``CH3Li`` in
the model, it reaches this family in the *reverse* generation direction: it matches ``CH3Li``
against the product template and reconstructs the reactant side, ``Li+ + CH3``. Critically, the
engine stores that reaction back in family-forward molecular orientation --
``reactants=[Li+, CH3], products=[CH3Li]`` -- while flagging it ``is_forward=False`` to record how
it was found. Because the stored orientation is family-forward, the electron is still a *reactant*
(the equation only balances as ``Li+ + CH3 + e- -> CH3Li``), so ``Reaction.electrons`` must stay at
the family-forward value ``-1``.

The defect this locks out: ``_create_reaction`` negated ``electrons`` whenever ``is_forward`` was
false, treating the ``is_forward`` flag as if it also reversed the molecular orientation. It does
not -- the reactant/product lists were swapped back to family-forward in the same constructor. The
negation therefore produced ``electrons = +1`` for a reaction whose electron is physically a
reactant, and the Chemkin writer correctly refused it as unbalanced in the ``E`` pseudo-element.

Two quantities, asserted separately (the point of the ticket). They are computed by two
independent routes so a single wrong scalar cannot satisfy both:

  * NET electron change  -- read straight off ``Reaction.electrons``; must be ``-1``.
  * REACTANT-side electron participation -- obtained by actually expanding the electron into the
    equation the writer builds (``expand_electrons``) and counting the electrons that land on the
    reactant side; must be exactly ``1`` on the reactant side and ``0`` on the product side.

A reaction whose net change said "+1 produced" but whose reactant participation was "1 consumed"
would be the exact conflation the electron-representation work exists to prevent; here both must
agree on consumption. The final assertion runs the writer's own ``check_electron_balance`` on the
expanded equation: with the correct sign it does not raise, so the guard passes *because the
reaction balances*, not because the guard was touched.

Runs against the real plasma database on ``settings['database.directory']`` (point it at
RMG-database-plasma), through the same family-generation calls the model builder makes.
"""

import os.path

import pytest

from rmgpy import settings
from rmgpy.data.kinetics.common import get_molecularity
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.electron_balance import check_electron_balance, expand_electrons
from rmgpy.molecule import Molecule
from rmgpy.species import Species

# C[Li] -- neutral methyllithium, net charge 0; the only lithium-bearing species the deck feeds.
NEUTRAL_CH3LI = "C[Li]"

FAMILY = "Cation_R_Recombination"

# The canonical electron pseudo-species (u0, charge -1), as the plasma database and the SMILES
# translator define it. Built once for the whole module.
ELECTRON = Species(label="e-", molecule=[Molecule().from_adjacency_list("1 e u0 p0 c-1")])


@pytest.mark.database
class TestCationReverseElectronSign:
    """The reverse-generated ``Li+ + CH3 <=> CH3Li`` carries the correct reactant-side electron."""

    @classmethod
    def setup_class(cls):
        families_path = os.path.join(settings["database.directory"], "kinetics", "families")
        database = KineticsDatabase()
        database.load_recommended_families(os.path.join(families_path, "recommended.py"))
        database.load_families(families_path, families=[FAMILY])
        cls.database = database
        cls.family = database.families[FAMILY]

        # Database-level generation, the same call the model builder makes: feed the neutral
        # product and let the family reach the cation via reverse generation.
        species = Species().from_smiles(NEUTRAL_CH3LI)
        species.generate_resonance_structures()
        cls.reactions = cls.database.generate_reactions_from_families(
            [species], only_families=[FAMILY]
        )

        # Select the reverse recombination Li+ + CH3 <=> CH3Li: two reactants (one of them the
        # lithium cation) and a single neutral product. There must be exactly one.
        matches = [
            rxn for rxn in cls.reactions
            if len(rxn.reactants) == 2 and len(rxn.products) == 1
            and any(cls._net_charge(spc) == 1 for spc in rxn.reactants)
        ]
        assert len(matches) == 1, (
            "expected exactly one reverse Li+ + CH3 <=> CH3Li reaction, got "
            "{0}: {1}".format(len(matches), [str(r) for r in cls.reactions])
        )
        cls.reaction = matches[0]

    # -- the reaction really is the reverse-generated cation recombination -----------------------

    def test_reaction_shape_and_direction(self):
        """It is a reverse-generated (is_forward=False) recombination of this family, charge +1 -> 0
        on the heavy species, so the electron must supply the missing charge on the reactant side."""
        rxn = self.reaction
        assert isinstance(rxn, TemplateReaction)
        assert rxn.family == FAMILY
        assert rxn.is_forward is False
        assert sorted(self._net_charge(spc) for spc in rxn.reactants) == [0, 1]
        assert [self._net_charge(spc) for spc in rxn.products] == [0]

    # -- quantity 1: NET electron change (read off the scalar) -----------------------------------

    def test_net_electron_change_is_minus_one(self):
        """The family declares ``electrons = -1`` for its forward orientation; the reverse-generated
        reaction is stored in that same orientation, so its net electron change is ``-1`` (one
        electron net consumed), NOT ``+1``. This is the assertion that fails RED before the fix."""
        assert self.reaction.electrons == -1

    # -- quantity 2: REACTANT-side electron participation (expanded, counted, not derived) --------

    def test_reactant_side_electron_participation_is_one(self):
        """Expand the electron into the equation the writer builds and count where it lands. It must
        appear once among the reactants and never among the products -- the electron is consumed."""
        reactants, products = expand_electrons(self.reaction, [ELECTRON])
        reactant_electrons = sum(1 for spc in reactants if spc.is_electron())
        product_electrons = sum(1 for spc in products if spc.is_electron())
        assert reactant_electrons == 1
        assert product_electrons == 0

    def test_molecularity_counts_the_consumed_electron(self):
        """Independently of the expansion, the production helper ``get_molecularity`` counts the
        consumed electron: two heavy reactants plus one incident electron -> molecularity 3. A
        released-electron (+1) sign would drop it to 2."""
        assert get_molecularity(self.reaction) == 3

    # -- the two quantities together let the writer's guard pass because it BALANCES --------------

    def test_expanded_equation_balances_in_E_and_charge(self):
        """With the correct sign, the writer's own ``check_electron_balance`` does not raise on the
        expanded equation, and the expanded reactant/product charges are equal (0 = 0). The guard
        passes because the reaction balances -- the guard itself is unchanged."""
        reactants, products = expand_electrons(self.reaction, [ELECTRON])
        equation = "{0} <=> {1}".format(
            " + ".join(str(s) for s in reactants), " + ".join(str(s) for s in products)
        )
        # Does not raise -> the E pseudo-element balances on what the writer would emit.
        check_electron_balance(self.reaction, reactants, products, equation)

        reactant_charge = sum(spc.molecule[0].get_net_charge() for spc in reactants)
        product_charge = sum(spc.molecule[0].get_net_charge() for spc in products)
        assert reactant_charge == 0
        assert product_charge == 0

    @staticmethod
    def _net_charge(species_or_molecule):
        """Net charge whether generation handed back a Molecule or a Species."""
        if isinstance(species_or_molecule, Molecule):
            return species_or_molecule.get_net_charge()
        return species_or_molecule.molecule[0].get_net_charge()
