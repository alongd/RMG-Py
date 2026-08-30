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
I-190 -- species resolution must be charge-aware.

An anion and a cation with the same heavy-atom skeleton (O- and O+, a single O
atom differing only in charge and lone pairs) are isomorphic under
``is_isomorphic(..., strict=False)``, which is charge-blind on purpose so that
resonance structures match. Several lookups across RMG use ``strict=False`` as a
species-identity test and so silently confuse opposite-charge species. The shared
guard ``is_isomorphic_same_charge`` requires equal NET charge around the
strict-False match; because net charge is conserved across resonance structures,
it never rejects a genuine resonance match.

These tests pin the helper directly. The call sites that use it are covered by
their own tests (``modelTest`` for check_for_existing_species) and by the fact
that the i186 argon deck no longer flips O- into O+.
"""

import pytest

from rmgpy.electron_balance import is_isomorphic_same_charge
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction, same_species_lists
from rmgpy.species import Species


O_CATION_ADJ = "multiplicity 2\n1 O u1 p2 c+1"
O_ANION_ADJ = "multiplicity 2\n1 O u1 p3 c-1"
O_NEUTRAL_ADJ = "multiplicity 3\n1 O u2 p2 c0"
N_CATION_ADJ = "multiplicity 3\n1 N u2 p1 c+1"
N_ANION_ADJ = "multiplicity 3\n1 N u2 p2 c-1"
N_NEUTRAL_ADJ = "multiplicity 4\n1 N u3 p1 c0"


class TestIsIsomorphicSameCharge:
    """Unit tests for rmgpy.electron_balance.is_isomorphic_same_charge."""

    def test_opposite_charges_do_not_match(self):
        """O+ and O- are isomorphic under strict=False but carry opposite net
        charge, so the guard must reject them."""
        cation = Species().from_adjacency_list(O_CATION_ADJ)
        anion = Species().from_adjacency_list(O_ANION_ADJ)
        # Precondition: charge-blind isomorphism would have matched them.
        assert cation.is_isomorphic(anion, strict=False)
        assert is_isomorphic_same_charge(cation, anion, strict=False) is False

    def test_cation_and_neutral_do_not_match(self):
        """A cation and its neutral atom differ in net charge and must not
        match."""
        cation = Species().from_adjacency_list(O_CATION_ADJ)
        neutral = Species().from_adjacency_list(O_NEUTRAL_ADJ)
        assert is_isomorphic_same_charge(cation, neutral, strict=False) is False

    def test_same_charge_isomorph_matches(self):
        """Two representations of the same charged species still match."""
        cation_a = Species().from_adjacency_list(O_CATION_ADJ)
        cation_b = Species().from_adjacency_list(O_CATION_ADJ)
        assert is_isomorphic_same_charge(cation_a, cation_b, strict=False) is True

    def test_accepts_molecule_on_right(self):
        """``other`` may be a Molecule (the coverage-dependence call sites pass a
        Molecule built from an adjacency-list key), not only a Species."""
        cation = Species().from_adjacency_list(O_CATION_ADJ)
        anion_mol = Molecule().from_adjacency_list(O_ANION_ADJ)
        cation_mol = Molecule().from_adjacency_list(O_CATION_ADJ)
        assert is_isomorphic_same_charge(cation, anion_mol, strict=False) is False
        assert is_isomorphic_same_charge(cation, cation_mol, strict=False) is True

    def test_different_skeleton_does_not_match_even_at_equal_charge(self):
        """Equal net charge is necessary but not sufficient: the isomorphism
        check still runs, so two different neutral species do not match."""
        neutral_o = Species().from_adjacency_list(O_NEUTRAL_ADJ)
        methane = Species().from_smiles("C")
        assert neutral_o.get_net_charge() == methane.get_net_charge() == 0
        assert is_isomorphic_same_charge(neutral_o, methane, strict=False) is False


class TestSameSpeciesListsChargeAware:
    """
    Reaction-level matching (same_species_lists / _same_object) must not treat
    two lists that differ only by which species carries which charge as the same.
    [O+, N-] and [O-, N+] have the same net-charge MULTISET, so a side-total check
    is insufficient: the guard has to bind charge to the individual species
    pairing the isomorphism chooses. Without it, reaction degeneracy silently sums
    two distinct charged channels and drops one (rmgpy/data/kinetics/common.py).
    """

    @staticmethod
    def _spc(adj):
        return Species().from_adjacency_list(adj)

    def test_charge_swapped_lists_are_not_the_same(self):
        """[O+, N-] vs [O-, N+]: charge-blind pairing would match them; the
        per-pair net-charge guard must reject every pairing."""
        op, nm = self._spc(O_CATION_ADJ), self._spc(N_ANION_ADJ)
        om, np_ = self._spc(O_ANION_ADJ), self._spc(N_CATION_ADJ)
        # Precondition: the O skeletons are isomorphic ignoring charge, so the
        # (O+~O-, N-~N+) pairing is exactly what a charge-blind match accepts.
        assert op.is_isomorphic(om, strict=False)
        assert nm.is_isomorphic(np_, strict=False)
        assert same_species_lists([op, nm], [om, np_], strict=False) is False

    def test_reordered_same_charge_lists_still_match(self):
        """A genuine reordering of the same charged species must still match:
        the guard rejects only charge-mismatched pairings, not valid ones."""
        op, nm = self._spc(O_CATION_ADJ), self._spc(N_ANION_ADJ)
        op2, nm2 = self._spc(O_CATION_ADJ), self._spc(N_ANION_ADJ)
        assert same_species_lists([op, nm], [nm2, op2], strict=False) is True

    def test_neutral_lists_are_unaffected(self):
        """Ordinary neutral chemistry is untouched: both sides are net-charge 0,
        so the guard is a pass-through and reordered neutral lists still match."""
        o, n = self._spc(O_NEUTRAL_ADJ), self._spc(N_NEUTRAL_ADJ)
        o2, n2 = self._spc(O_NEUTRAL_ADJ), self._spc(N_NEUTRAL_ADJ)
        assert o.get_net_charge() == n.get_net_charge() == 0
        assert same_species_lists([o, n], [n2, o2], strict=False) is True

    def test_reaction_isomorphic_distinguishes_charge_swapped_channels(self):
        """The real defect path: Reaction.is_isomorphic -> same_species_lists.
        Two reactions whose reactant charges are swapped are distinct channels,
        so degeneracy must not fold them together."""
        product = self._spc(O_NEUTRAL_ADJ)
        rxn1 = Reaction(reactants=[self._spc(O_CATION_ADJ), self._spc(N_ANION_ADJ)],
                        products=[product])
        rxn2 = Reaction(reactants=[self._spc(O_ANION_ADJ), self._spc(N_CATION_ADJ)],
                        products=[product])
        assert rxn1.is_isomorphic(rxn2, either_direction=False, strict=False) is False
