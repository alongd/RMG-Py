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
A scission daughter's pool proxy keeps the radical of the product it is built
from. A backbone homolysis piece carries one radical at the atom the bond
broke; the wing removal that turns the piece into an end group cuts a bond at
that same atom, and the cut must open its own stitch site. The daughter's
proxy is then neutral, with the piece's radical count at the piece's radical
site, and the producer never hands back a charged proxy for a neutral piece.

The parents are styrene-, methyl methacrylate- and ethylene-like test proxies
(test data only). A setup step that fails raises RuntimeError, never an
AssertionError.
"""

import pytest

from rmgpy.molecule import Bond, Molecule
from rmgpy.polymer import Polymer
import rmgpy.polymer as polymer

PARENTS = {
    'PS': '[CH2][CH]c1ccccc1',
    'PMMA': '[CH2][C](C)C(=O)OC',
    'PE': '[CH2][CH2]',
}


def _setup(ok, what):
    """A setup step that fails raises RuntimeError, never an AssertionError."""
    if not ok:
        raise RuntimeError('setup failed: ' + what)


def _parent(name):
    return Polymer(label=name, monomer=PARENTS[name], end_groups=['[CH3]', '[H]'], cutoff=3,
                   Mn=5000.0, Mw=6000.0, initial_mass=1.0)


def _backbone_cuts(base):
    for i, a in enumerate(base.atoms):
        for b in list(a.bonds):
            j = base.atoms.index(b)
            bond = base.get_bond(a, b)
            if (j > i and a.is_carbon() and b.is_carbon() and bond.is_single()
                    and not base.is_bond_in_cycle(bond)):
                yield i, j


def _split(base, i, j, concerted):
    """Cut bond i-j: a homolysis (two radicals) or a concerted split (an H
    from a carbon neighbour of i moves onto j and closes a C=C at i)."""
    m = base.copy(deep=True)
    ai, aj = m.atoms[i], m.atoms[j]
    m.remove_bond(m.get_bond(ai, aj))
    if concerted:
        n = next((x for x in ai.bonds if x.is_carbon() and any(h.is_hydrogen() for h in x.bonds)), None)
        if n is None:
            return []
        h = next(x for x in n.bonds if x.is_hydrogen())
        m.remove_bond(m.get_bond(n, h))
        m.add_bond(Bond(h, aj, order=1))
        m.get_bond(ai, n).increment_order()
    else:
        ai.increment_radical()
        aj.increment_radical()
    pieces = m.split()
    for p in pieces:
        p.update()
    return pieces if len(pieces) == 2 else []


def _pieces(parent, concerted):
    """Every distinct piece of every backbone cut of the parent's proxy."""
    base = parent.molecule[0]
    seen, out = set(), []
    for i, j in _backbone_cuts(base):
        for piece in _split(base, i, j, concerted):
            key = piece.to_smiles()
            if key not in seen:
                seen.add(key)
                out.append(piece)
    return out


def _radical_sites(mol):
    """Where the radicals sit: (element, radical electrons, H neighbours,
    heavy neighbours) of every radical atom, sorted."""
    return sorted((a.symbol, a.radical_electrons,
                   sum(1 for b in a.bonds if b.is_hydrogen()),
                   sum(1 for b in a.bonds if not b.is_hydrogen()))
                  for a in mol.atoms if a.radical_electrons)


def _scission_daughters(name):
    parent = _parent(name)
    out = []
    for piece in _pieces(parent, concerted=False):
        d = parent.create_reacted_copy(piece.copy(deep=True))
        if isinstance(d, Polymer) and (getattr(d, 'spawn_metadata', None) or {}).get('source') == 'scission':
            out.append((piece, d))
    _setup(out, f'no homolysis piece of {name} makes a scission daughter')
    return out


@pytest.mark.parametrize('name', sorted(PARENTS))
def test_a_homolysis_pieces_scission_daughter_keeps_its_radical(name):
    """Neutral, with the piece's radical count, at the piece's radical site."""
    for piece, daughter in _scission_daughters(name):
        proxy = daughter.molecule[0]
        assert (proxy.get_net_charge(), proxy.get_radical_count(), _radical_sites(proxy)) == \
            (0, piece.get_radical_count(), _radical_sites(piece)), (piece.to_smiles(), proxy.to_smiles())


@pytest.mark.parametrize('name', sorted(PARENTS))
def test_the_producer_never_returns_a_charged_proxy_for_a_neutral_piece(name):
    parent = _parent(name)
    charged = []
    for concerted in (False, True):
        for piece in _pieces(parent, concerted):
            _setup(piece.get_net_charge() == 0, 'neutral piece')
            d = parent.create_reacted_copy(piece.copy(deep=True))
            if isinstance(d, Polymer) and any(m.get_net_charge() != 0 for m in d.molecule):
                charged.append((piece.to_smiles(), d.label, d.molecule[0].to_smiles()))
    assert charged == []


def test_a_cut_at_a_radical_atom_opens_one_more_site():
    """The ethyl radical with its methyl cut away as a head wing: the kept
    CH2 already carries the piece's radical and gains the cut's stitch site."""
    ethyl = Molecule().from_smiles('[CH2]C')
    ch2 = next(a for a in ethyl.atoms if a.radical_electrons == 1)
    ch3 = next(a for a in ethyl.atoms if a.is_carbon() and a is not ch2)
    removed = {ch3} | {h for h in ch3.bonds if h.is_hydrogen()}
    kept = polymer.Polymer._extract_remainder(ethyl, removed)
    polymer.Polymer._restore_labels(kept, original_mol=ethyl, removed_atoms=removed,
                                    head_match_atoms=removed, tail_match_atoms=None)
    site = next(a for a in kept.atoms if a.label == '*2')
    assert site.radical_electrons == 2


def test_the_generated_homolysis_rows_tail_is_a_neutral_radical():
    """Through the production generation seam (make_new_reaction): the
    scission tail of a backbone homolysis of the styrene-like parent."""
    from rmgpy.data.kinetics.family import TemplateReaction
    from rmgpy.kinetics import Arrhenius
    from rmgpy.rmg.model import CoreEdgeReactionModel
    from rmgpy.species import Species
    parent = _parent('PS')
    model = CoreEdgeReactionModel()
    model._register_polymer(parent, generate_thermo=False)
    gas = Species(molecule=[Molecule().from_smiles('[CH2]C(C)c1ccccc1')])
    gas.generate_resonance_structures()
    base = parent.molecule[0]
    pieces = next((p for i, j in _backbone_cuts(base) for p in [_split(base, i, j, False)]
                   if p and any(gas.is_isomorphic(Molecule().from_smiles(x.to_smiles())) for x in p)), None)
    _setup(pieces is not None, 'the cut releasing the end unit')
    raw = TemplateReaction(reactants=[base.copy(deep=True)], products=pieces, family='test_family',
                           template=['a'], is_forward=True, degeneracy=1,
                           kinetics=Arrhenius(A=(1.0, 's^-1'), n=0.0, Ea=(0.0, 'kcal/mol'), T0=(298.15, 'K')))
    out, _ = model.make_new_reaction(raw, check_existing=False, generate_thermo=False, generate_kinetics=False)
    _setup(out is not None, 'row generated')
    tails = [p for p in out.products if isinstance(p, Polymer)]
    _setup(len(tails) == 1 and tails[0].label == 'PS_scission_tail', 'one scission tail product')
    proxy = tails[0].molecule[0]
    assert (proxy.get_net_charge(), proxy.get_radical_count()) == (0, 1)


@pytest.mark.parametrize('name,piece', [
    ('PS', '[CH2]C(CC(C)c1ccccc1)c1ccccc1'),            # the radical at the end the wing does not cut
])
def test_a_piece_whose_radical_is_away_from_the_cut_is_not_booked(name, piece):
    """The end group keeps the radicals the product left on the cut atom
    only; a piece with its radical elsewhere is refused as before."""
    assert _parent(name).create_reacted_copy(Molecule().from_smiles(piece)) is None


@pytest.mark.parametrize('name,piece', [
    ('PS', 'C[CH]CC(CC(C)c1ccccc1)c1ccccc1'),          # a phenyl loss: the backbone stays whole
    ('PMMA', 'COC(=O)C(C)(C)CC(C)(C[C](C)C)C(=O)OC'),   # an ester loss
])
def test_a_side_group_loss_is_never_a_scission_piece(name, piece):
    """A piece that keeps the whole backbone lost a side group: never a
    half-length scission population, but its X-loss feature (I-257)."""
    d = _parent(name).create_reacted_copy(Molecule().from_smiles(piece))
    assert (getattr(d, 'spawn_metadata', None) or {}).get('source') == 'side_group_loss'
