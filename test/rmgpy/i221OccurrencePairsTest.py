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
I-221 round 112 -- a reaction is a list of OCCURRENCES, and `pairs` must be carried by them.

`Species.__eq__` is identity, and one species object can occur twice on a side (``A + A``)
and on both sides -- every electron-impact ionisation holds the one electron object three
times, ``Li + e- => Li+ + e- + e-``. Every remapper of `pairs` used to key on the species
(``id(species)`` in `Reaction.copy`, a global ``list.index`` in
`CoreEdgeReactionModel.make_new_reaction`, an isomorphism search in `ensure_species`), so
all occurrences of one species collapsed onto whichever one the map met first or last.

Each assertion here is about POSITIONS: which reactant slot and which product slot a pair
names, by identity. A test that compared labels or species sets could not see the defect,
because every collapsed answer still names a species with the right label.

The five controls the round's verifier names: the same species twice on one side; a species
on both sides; duplicated electron products; a pair member that is not on its side (which
must raise); and an ordinary thermal reaction (unchanged).
"""

import pytest

from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.molecule.molecule import Molecule
from rmgpy.reaction import Reaction, pair_occurrences, pairs_from_occurrences
from rmgpy.rmg.model import CoreEdgeReactionModel, ReactionModel
from rmgpy.species import Species
from rmgpy.thermo import ThermoData


def _spc(label, adjlist=None, smiles=None):
    spc = Species(label=label)
    return spc.from_adjacency_list(adjlist) if adjlist else spc.from_smiles(smiles)


def _ionisation_species():
    return (_spc('Li', adjlist='1 Li u1 p0 c0'),
            _spc('Li+', adjlist='1 Li u0 p0 c+1'),
            _spc('e', adjlist='1 e u1 p0 c-1'))


def _slots(side, member):
    """Every position on `side` that holds `member` itself."""
    return [i for i, occupant in enumerate(side) if occupant is member]


def _positions(reaction):
    """`pairs` as (reactant position, product position), refusing an ambiguous answer."""
    positions = []
    for reactant, product in reaction.pairs:
        r, p = _slots(reaction.reactants, reactant), _slots(reaction.products, product)
        assert len(r) == 1, 'pair reactant {0!r} occupies reactant slots {1}'.format(reactant, r)
        assert len(p) == 1, 'pair product {0!r} occupies product slots {1}'.format(product, p)
        positions.append((r[0], p[0]))
    return positions


class TestCopyCarriesOccurrences:
    """`Reaction.copy` deepens each occurrence separately, so every pair must name one slot."""

    def test_same_species_twice_on_one_side(self):
        s, p = _spc('S', smiles='[CH3]'), _spc('P', smiles='CC')
        rxn = Reaction(reactants=[s, s], products=[p], pairs=[(s, p), (s, p)])
        other = rxn.copy()
        assert other.reactants[0] is not other.reactants[1]
        # At 71ae97bd5 both pairs named reactant slot 1: `id(s)` mapped to the last copy.
        assert _positions(other) == [(0, 0), (1, 0)]

    def test_species_on_both_sides(self):
        x, y, z = _spc('X', smiles='[H]'), _spc('Y', smiles='[CH3]'), _spc('Z', smiles='C')
        rxn = Reaction(reactants=[x, y], products=[x, z], pairs=[(x, x), (y, z)])
        other = rxn.copy()
        assert other.reactants[0] is not other.products[0]
        assert _positions(other) == [(0, 0), (1, 1)]

    def test_duplicated_electron_products(self):
        li, li_plus, e = _ionisation_species()
        rxn = Reaction(reactants=[li, e], products=[li_plus, e, e])
        rxn.generate_pairs()
        # generate_pairs' own shape, by slot on the original: e-e, Li-Li+, Li-e.
        mentioned = [(rxn.reactants.index(a), [i for i, x in enumerate(rxn.products) if x is b])
                     for a, b in rxn.pairs]
        assert mentioned == [(1, [1, 2]), (0, [0]), (0, [1, 2])]
        other = rxn.copy()
        assert len({id(spc) for spc in other.products}) == 3
        # The first mention of the product electron is slot 1 and the second slot 2. At
        # 71ae97bd5 both named slot 2, and product slot 1 was named by no pair at all.
        assert _positions(other) == [(1, 1), (0, 0), (0, 2)]
        assert sorted(p for _, p in _positions(other)) == [0, 1, 2]

    def test_missing_pair_member_raises(self):
        a, b, alien = _spc('A', smiles='C'), _spc('B', smiles='[CH3]'), _spc('Q', smiles='O')
        rxn = Reaction(reactants=[a], products=[b], pairs=[(a, alien)])
        # At 71ae97bd5 the alien passed through aliased to the original object.
        with pytest.raises(ValueError, match='not a product of this reaction'):
            rxn.copy()

    def test_member_on_the_wrong_side_raises(self):
        a, b = _spc('A', smiles='C'), _spc('B', smiles='[CH3]')
        rxn = Reaction(reactants=[a], products=[b], pairs=[(b, a)])
        with pytest.raises(ValueError, match='not a reactant of this reaction'):
            rxn.copy()

    def test_ordinary_thermal_reaction_unchanged(self):
        a, b, c = _spc('A', smiles='CC'), _spc('B', smiles='[CH3]'), _spc('C', smiles='[CH3]')
        rxn = Reaction(reactants=[a], products=[b, c])
        rxn.generate_pairs()
        before = _positions(rxn)
        other = rxn.copy()
        assert _positions(other) == before == [(0, 0), (0, 1)]
        assert all(o is not s for o, s in zip(other.reactants + other.products,
                                               rxn.reactants + rxn.products))

    def test_no_pairs_stays_none(self):
        a, b = _spc('A', smiles='C'), _spc('B', smiles='[CH3]')
        assert Reaction(reactants=[a], products=[b]).copy().pairs is None


class TestTemplateCopyCarriesOccurrences:
    """`TemplateReaction.copy` is one pickle round trip; the memo keeps the identity
    structure, so a shared electron stays shared and every pair names the copy's own objects."""

    def test_electron_identity_structure_survives(self):
        li, li_plus, e = _ionisation_species()
        rxn = TemplateReaction(reactants=[li, e], products=[li_plus, e, e], family='F')
        rxn.generate_pairs()
        other = rxn.copy()
        assert other.reactants[1] is other.products[1] is other.products[2]
        assert other.reactants[1] is not e
        own = {id(spc) for spc in other.reactants + other.products}
        assert all(id(a) in own and id(b) in own for a, b in other.pairs)
        assert pair_occurrences(other.pairs, other.reactants, other.products) == \
            pair_occurrences(rxn.pairs, rxn.reactants, rxn.products)

    def test_alien_member_raises(self):
        a, b, alien = _spc('A', smiles='C'), _spc('B', smiles='[CH3]'), _spc('Q', smiles='O')
        rxn = TemplateReaction(reactants=[a], products=[b], pairs=[(a, alien)], family='F')
        # At 71ae97bd5 the pickle memo handed the copy a fresh detached alien, no error.
        with pytest.raises(ValueError, match='not a product of this reaction'):
            rxn.copy()


class TestPairOccurrences:
    """The resolver every remapper goes through."""

    def test_kth_mention_is_kth_occurrence(self):
        li, li_plus, e = _ionisation_species()
        rxn = Reaction(reactants=[li, e], products=[li_plus, e, e])
        rxn.generate_pairs()
        assert pair_occurrences(rxn.pairs, rxn.reactants, rxn.products) == \
            [(1, 1), (0, 0), (0, 2)]

    def test_rebuild_on_new_lists(self):
        new_r = [object(), object()]
        new_p = [object(), object(), object()]
        pairs = pairs_from_occurrences([(1, 1), (0, 0), (0, 2)], new_r, new_p)
        assert all(a is new_r[i] and b is new_p[j]
                   for (a, b), (i, j) in zip(pairs, [(1, 1), (0, 0), (0, 2)]))

    def test_malformed_pair_raises(self):
        a = object()
        with pytest.raises(ValueError):
            pair_occurrences([(a,)], [a], [a])

    def test_none_is_none(self):
        assert pair_occurrences(None, [], []) is None
        assert pairs_from_occurrences(None, [], []) is None


class TestModelCarriesPairsByPosition:
    """HIGH 2: `make_new_reaction` swaps each participant for its model species. The pairs
    follow by position; a global ``list.index`` answered the first slot for every mention."""

    def _model(self, monkeypatch):
        model = CoreEdgeReactionModel()
        # One distinct model object per OCCURRENCE, so a collapsed mapping is visible. The
        # real make_new_species maps both electrons to one model species, which hides it.
        monkeypatch.setattr(model, 'make_new_species',
                            lambda spc, **kwargs: (_spc(spc.label, adjlist=spc.molecule[0].to_adjacency_list()), True))
        monkeypatch.setattr(model, 'register_reaction', lambda rxn: None)
        return model

    def test_duplicated_electron_products(self, monkeypatch):
        li, li_plus, e = _ionisation_species()
        forward = TemplateReaction(reactants=[li, e], products=[li_plus, e, e],
                                   family='F', is_forward=False)
        forward.generate_pairs()
        rxn, new = self._model(monkeypatch).make_new_reaction(
            forward, check_existing=False, generate_thermo=False, generate_kinetics=False)
        assert new
        assert all(spc is not orig for spc, orig in zip(rxn.products, [li_plus, e, e]))
        # At 71ae97bd5: `forward.products.index(e)` is 1 for both electron mentions.
        assert _positions(rxn) == [(1, 1), (0, 0), (0, 2)]

    def test_reverse_pairs_follow(self, monkeypatch):
        li, li_plus, e = _ionisation_species()
        forward = TemplateReaction(reactants=[li, e], products=[li_plus, e, e],
                                   family='F', is_forward=False)
        forward.generate_pairs()
        forward.reverse = TemplateReaction(reactants=[li_plus, e, e], products=[li, e],
                                           family='F', is_forward=True)
        forward.reverse.pairs = [(b, a) for a, b in forward.pairs]
        rxn, _ = self._model(monkeypatch).make_new_reaction(
            forward, check_existing=False, generate_thermo=False, generate_kinetics=False)
        assert rxn.reverse.pairs == [(b, a) for a, b in rxn.pairs]
        # Round 113: the reverse pairs must name occurrences of the reverse's OWN sides. At
        # fc60e5ba4 the sides were never rebuilt, so every member was detached and this raised.
        assert pair_occurrences(rxn.reverse.pairs, rxn.reverse.reactants, rxn.reverse.products) == \
            [(1, 1), (0, 0), (2, 0)]
        assert _positions(rxn.reverse) == [(1, 1), (0, 0), (2, 0)]

    def test_reverse_sides_are_the_model_species(self, monkeypatch):
        """The reverse `add_reverse_attribute` attaches is built from Molecules, in the
        template's order. Its sides become the forward's model species, swapped, one slot per
        occurrence, in lists of their own."""
        li, li_plus, e = _ionisation_species()
        forward = TemplateReaction(reactants=[li, e], products=[li_plus, e, e],
                                   family='F', is_forward=False)
        forward.generate_pairs()
        e1, e2, e3 = (e.molecule[0].copy(deep=True) for _ in range(3))
        m_li_plus, m_li = li_plus.molecule[0].copy(deep=True), li.molecule[0].copy(deep=True)
        forward.reverse = TemplateReaction(reactants=[e1, m_li_plus, e2], products=[e3, m_li],
                                           family='F', is_forward=True)
        forward.reverse.pairs = [(e1, e3), (m_li_plus, m_li), (e2, m_li)]
        rxn, _ = self._model(monkeypatch).make_new_reaction(
            forward, check_existing=False, generate_thermo=False, generate_kinetics=False)
        assert all(a is b for a, b in zip(rxn.reverse.reactants, rxn.products))
        assert all(a is b for a, b in zip(rxn.reverse.products, rxn.reactants))
        assert len(rxn.reverse.reactants) == 3 and len(rxn.reverse.products) == 2
        assert rxn.reverse.reactants is not rxn.products and rxn.reverse.products is not rxn.reactants
        assert _positions(rxn.reverse) == [(b, a) for a, b in _positions(rxn)]

    def test_ordinary_thermal_reaction_unchanged(self, monkeypatch):
        a, b, c = _spc('A', smiles='CC'), _spc('B', smiles='[CH3]'), _spc('C', smiles='[H]')
        forward = TemplateReaction(reactants=[a], products=[b, c], family='F', is_forward=False)
        forward.generate_pairs()
        rxn, _ = self._model(monkeypatch).make_new_reaction(
            forward, check_existing=False, generate_thermo=False, generate_kinetics=False)
        assert _positions(rxn) == [(0, 0), (0, 1)]


class TestEnsureSpeciesCarriesOccurrences:
    """`ensure_species` wraps Molecule participants in Species; pairs follow by position."""

    def test_duplicated_electron_products(self):
        li = Molecule().from_adjacency_list('1 Li u1 p0 c0')
        li_plus = Molecule().from_adjacency_list('1 Li u0 p0 c+1')
        e = Molecule().from_adjacency_list('1 e u1 p0 c-1')
        rxn = Reaction(reactants=[li, e], products=[li_plus, e, e],
                       pairs=[(e, e), (li, li_plus), (li, e)])
        rxn.ensure_species()
        assert all(isinstance(spc, Species) for spc in rxn.reactants + rxn.products)
        assert len({id(spc) for spc in rxn.products}) == 3
        # At 71ae97bd5 an isomorphism search mapped both electron mentions to slot 1.
        assert _positions(rxn) == [(1, 1), (0, 0), (0, 2)]


class TestMergeModelsCarriesPairsByPosition:
    """`merge_models` swaps species for the final model's; pairs follow by position."""

    def test_species_on_both_sides_and_twice(self):
        li, li_plus, e = _ionisation_species()
        common_e = _spc('e', adjlist='1 e u1 p0 c-1')
        # merge() compares the thermo of every species the two models share.
        e.thermo = common_e.thermo = ThermoData(
            Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'), Cpdata=([20.8] * 7, 'J/(mol*K)'),
            H298=(0, 'kJ/mol'), S298=(20, 'J/(mol*K)'), Cp0=(20.8, 'J/(mol*K)'),
            CpInf=(20.8, 'J/(mol*K)'))
        first = ReactionModel(species=[common_e], reactions=[])
        rxn = Reaction(reactants=[li, e], products=[li_plus, e, e])
        rxn.generate_pairs()
        second = ReactionModel(species=[li, li_plus, e], reactions=[rxn])
        merged = first.merge(second)
        out = merged.reactions[-1]
        assert out.reactants[1] is common_e and out.products[1] is common_e \
            and out.products[2] is common_e
        assert pair_occurrences(out.pairs, out.reactants, out.products) == \
            [(1, 1), (0, 0), (0, 2)]
        assert all(b is common_e for a, b in (out.pairs[0], out.pairs[2]))
