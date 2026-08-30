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
from rmgpy.species import Species


O_CATION_ADJ = "multiplicity 2\n1 O u1 p2 c+1"
O_ANION_ADJ = "multiplicity 2\n1 O u1 p3 c-1"
O_NEUTRAL_ADJ = "multiplicity 3\n1 O u2 p2 c0"


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
