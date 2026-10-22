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


import pytest
from numpy import testing

from rmgpy.molecule.graph import Edge, Graph, Vertex, get_vertex_connectivity_value
from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.vf2 import VF2


class TestVF2:
    """
    Contains unit tests of the methods for computing symmetry numbers for a
    given Molecule object.
    """

    def setup_class(self):
        self.vf2 = VF2()
        self.mol = Molecule().from_smiles("CC(=O)C[CH2]")
        self.mol2 = self.mol.copy(deep=True)

    def test_import_graph(self):
        """Test that we can add graphs to the object and that they are sorted"""

        self.mol.sort_vertices()

        ordered_original_connectivity_order = [get_vertex_connectivity_value(atom) for atom in self.mol.atoms]

        self.vf2.graphA = self.mol
        self.vf2.graphB = self.mol2

        final_connectivity_order = [get_vertex_connectivity_value(atom) for atom in self.vf2.graphA.atoms]
        final_connectivity_order2 = [get_vertex_connectivity_value(atom) for atom in self.vf2.graphB.atoms]

        testing.assert_array_equal(final_connectivity_order, final_connectivity_order2)
        testing.assert_array_equal(final_connectivity_order, ordered_original_connectivity_order)

    def test_feasible(self):
        """
        Test that feasibility returns correct values on highly functional molecule

        `feasible` method isn't perfect in assigning values but it should do a good
        job on highly functional values
        """

        self.vf2.graphA = self.mol
        self.vf2.graphB = self.mol2

        for atom1 in self.vf2.graphA.atoms:
            for atom2 in self.vf2.graphB.atoms:
                # same connectivity values should result in `feasible` being true
                if get_vertex_connectivity_value(atom1) == get_vertex_connectivity_value(atom2):
                    assert self.vf2.feasible(atom1, atom2)
                else:  # different connectivity values should return false
                    assert not self.vf2.feasible(atom1, atom2)

    def test_clear_mapping(self):
        """Test that vertex mapping is cleared after isomorphism."""
        self.vf2.is_isomorphic(self.mol, self.mol2, None)

        for atom in self.mol.atoms:
            assert atom.mapping is None
            assert not atom.terminal

    def test_save_order_on_graphs_compared_before(self):
        """A save_order comparison of graphs the matcher compared just before
        (it keeps them and does not re-sort) restores their order, never
        raises on an order it did not save."""
        a = Molecule().from_smiles("[CH2]CCC")
        b = Molecule().from_smiles("[CH2]CCC")
        assert self.vf2.is_isomorphic(a, b, None)
        order_a, order_b = list(a.atoms), list(b.atoms)
        assert self.vf2.is_isomorphic(a, b, None, save_order=True)
        assert a.atoms == order_a and b.atoms == order_b

    def test_a_kept_graph_whose_connectivity_was_reset_is_re_sorted(self):
        """A graph the matcher kept from its previous call whose connectivity
        values were reset since is re-sorted (its connectivity recomputed),
        never matched with the reset values."""
        a = Molecule().from_smiles("CCC")
        b = Molecule().from_smiles("CCC")
        assert self.vf2.is_isomorphic(a, b, None)
        a.reset_connectivity_values()
        assert self.vf2.is_isomorphic(a, b, None)

    def test_a_kept_graph_whose_vertex_count_changed_is_re_sorted(self):
        """A graph the matcher kept from its previous call that lost a vertex
        since (its sorting labels still valid) is re-sorted, never matched
        with its connectivity from before the change."""
        def path(n):
            graph = Graph()
            vertices = [Vertex() for _ in range(n)]
            for vertex in vertices:
                graph.add_vertex(vertex)
            for v1, v2 in zip(vertices, vertices[1:]):
                graph.add_edge(Edge(v1, v2))
            return graph, vertices

        four, vertices = path(4)
        assert self.vf2.is_isomorphic(four, path(4)[0], None)
        four.remove_vertex(vertices[-1])
        assert self.vf2.is_isomorphic(four, path(3)[0], None)

    @pytest.mark.parametrize("subgraph", [False, True])
    @pytest.mark.parametrize("empty", ["both", "first", "second"])
    def test_save_order_on_empty_graphs(self, subgraph, empty):
        """``save_order`` on an empty graph (both, or either one) gives the
        answer of the plain call and keeps both orders, never refuses to
        restore the (empty) order it saved."""
        def pair():
            return (Graph() if empty in ("both", "first") else Molecule().from_smiles("CC"),
                    Graph() if empty in ("both", "second") else Molecule().from_smiles("CC"))

        def compare(graph1, graph2, save_order):
            if subgraph:
                return VF2().is_subgraph_isomorphic(graph1, graph2, None, save_order=save_order)
            return VF2().is_isomorphic(graph1, graph2, None, save_order=save_order)

        plain = compare(*pair(), False)
        graph1, graph2 = pair()
        order = [list(graph1.vertices), list(graph2.vertices)]
        assert compare(graph1, graph2, True) == plain
        assert [list(graph1.vertices), list(graph2.vertices)] == order

    @pytest.mark.parametrize("subgraph", [False, True])
    def test_save_order_restores_the_order_when_the_sizes_reject(self, subgraph):
        """``save_order`` restores both graphs' order also when the sizes
        alone reject the pair (an early return after sorting)."""
        def pair():
            big, small = Molecule().from_smiles("CCO"), Molecule().from_smiles("CC")
            return (small, big) if subgraph else (big, small)

        def compare(graph1, graph2, save_order):
            if subgraph:
                return VF2().is_subgraph_isomorphic(graph1, graph2, None, save_order=save_order)
            return VF2().is_isomorphic(graph1, graph2, None, save_order=save_order)

        graph1, graph2 = pair()
        order = [list(graph1.vertices), list(graph2.vertices)]
        assert not compare(graph1, graph2, False)
        assert [list(graph1.vertices), list(graph2.vertices)] != order   # the call sorts them
        graph1, graph2 = pair()
        order = [list(graph1.vertices), list(graph2.vertices)]
        assert not compare(graph1, graph2, True)
        assert [list(graph1.vertices), list(graph2.vertices)] == order


def _numbered(smiles):
    molecule = Molecule().from_smiles(smiles)
    molecule.assign_atom_ids()
    return molecule


def _ids(mapping):
    return None if mapping is None else {a.id: b.id for a, b in mapping.items()}


FIRST_CASES = [
    # (smiles a, smiles b, pinned atom indices, strict)
    ("CC(C)(C)C", "CC(C)(C)C", (), False),
    ("CC(C)(C)C", "CC(C)(C)C", ((1, 1),), True),
    ("CC(C)C", "CC(C)C", ((0, 2),), False),
    ("CC(C)C", "CC(C)C", ((0, 1),), False),          # a methyl pinned on the centre: none
    ("[CH2]C=CC", "C=C[CH]C", (), False),
    ("[CH2]C=CC", "C=C[CH]C", (), True),             # other resonance drawings, strict: none
    ("c1ccccc1CC", "c1ccccc1CC", ((6, 6),), False),
    ("FC(F)(F)C(F)(F)[C](F)C(F)(F)F", "FC(F)(F)C(F)(F)[C](F)C(F)(F)F", (), False),
]


class TestFindFirstIsomorphism:
    """``find_first_isomorphism`` is ``find_isomorphism``'s first mapping (the
    same search order: listing every mapping only continues past the first),
    or with ``accept`` the first listed mapping it accepts, without listing
    the others."""

    @pytest.mark.parametrize("save_order", [False, True])
    @pytest.mark.parametrize("smiles_a,smiles_b,pins,strict", FIRST_CASES)
    def test_the_first_mapping_is_the_first_one_find_isomorphism_lists(
            self, smiles_a, smiles_b, pins, strict, save_order):
        a, b = _numbered(smiles_a), _numbered(smiles_b)

        def initial(x, y):
            return {x.atoms[i]: y.atoms[j] for i, j in pins} or None

        x, y = a.copy(deep=True), b.copy(deep=True)
        listed = VF2().find_isomorphism(x, y, initial(x, y), save_order=save_order, strict=strict)
        x, y = a.copy(deep=True), b.copy(deep=True)
        first = VF2().find_first_isomorphism(x, y, initial(x, y), save_order=save_order, strict=strict)
        expected = _ids(listed[0]) if listed else None
        assert _ids(first) == expected
        x, y = a.copy(deep=True), b.copy(deep=True)
        assert _ids(x.find_first_isomorphism(y, initial_map=initial(x, y), save_order=save_order,
                                             strict=strict)) == expected

    @pytest.mark.parametrize("smiles_a,smiles_b,pins,strict", FIRST_CASES)
    def test_with_accept_it_is_the_first_listed_mapping_accepted(self, smiles_a, smiles_b, pins, strict):
        a, b = _numbered(smiles_a), _numbered(smiles_b)
        listed = a.copy(deep=True).find_isomorphism(b.copy(deep=True), strict=strict)
        if not listed:
            assert a.copy(deep=True).find_first_isomorphism(
                b.copy(deep=True), strict=strict, accept=lambda m: True) is None
            return
        heavy = [atom.id for atom in a.atoms if not atom.is_hydrogen()]
        # accept only mappings sending the last heavy atom onto the largest
        # image any listed mapping gives it
        target = max(_ids(m)[heavy[-1]] for m in listed)
        seen = []

        def accept(mapping):
            seen.append(1)
            return _ids(mapping)[heavy[-1]] == target

        expected = next(_ids(m) for m in listed if _ids(m)[heavy[-1]] == target)
        found = a.copy(deep=True).find_first_isomorphism(b.copy(deep=True), strict=strict, accept=accept)
        assert _ids(found) == expected
        # it stops at the accepted mapping: the ones before it, and that one
        assert len(seen) == [_ids(m) for m in listed].index(expected) + 1
        assert a.copy(deep=True).find_first_isomorphism(
            b.copy(deep=True), strict=strict, accept=lambda m: False) is None

    def test_none_when_not_isomorphic(self):
        a, b = _numbered("CCO"), _numbered("COC")
        assert VF2().find_first_isomorphism(a, b, None) is None
        assert a.find_first_isomorphism(b) is None
        assert _numbered("[CH2]C").find_first_isomorphism(_numbered("C=C")) is None   # the fingerprint differs

    def test_the_vertex_order_is_kept_with_save_order(self):
        a, b = _numbered("CC(C)(C)C"), _numbered("CC(C)(C)C")
        order = [list(a.atoms), list(b.atoms)]
        assert a.find_first_isomorphism(b, save_order=True) is not None
        assert [list(a.atoms), list(b.atoms)] == order
