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
The identity of an H-loss feature daughter (``{label}_mod``). Pools register
by fingerprint, and a radical feature unit's fingerprint read each radical
atom only as (element, radicals, H count, heavy-neighbour count). Two
H-loss sites with that same reading but different structures (a methyl on an
ester oxygen and a methyl on a backbone carbon; the ortho, meta and para ring
positions) were registered as ONE pool: the second site's rows were booked to
the first site's structure, and RMG kept one rate for both.

- a pool is one per isomorphism class of the feature proxies: distinct sites
  register distinct pools;
- a row books the structure that reacted: its daughter's feature proxy is the
  row's raw product, rendered at the center;
- positional twins (the center and an end unit rendering of one monomer site)
  still book one pool -- the same structure, read at another position.

Monomers are a methyl methacrylate-, a styrene- and an ethylene-like repeat
unit (test data only). A setup step that fails raises RuntimeError, never an
AssertionError.
"""

import itertools

import pytest

from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.polymer import Polymer, stitch_molecules_by_labeled_atoms
from rmgpy.rmg.model import CoreEdgeReactionModel

MONOMERS = {
    'methacrylate': ('[CH2][C](C)C(=O)OC', ['[CH3]', '[H]']),
    'styrene': ('[CH2][CH]c1ccccc1', ['[CH3]', '[H]']),
    'ethylene': ('[CH2][CH2]', ['[H]', '[H]']),
}
SOURCE = 'radical_feature_h_loss'
KIN = Arrhenius(A=(1.0, 'cm^3/(mol*s)'), n=0.0, Ea=(0.0, 'kcal/mol'), T0=(298.15, 'K'))
RED = pytest.mark.xfail(strict=True, raises=AssertionError,
                        reason='I-259 red: same-reading H-loss sites share one pool')


def _setup(ok, what):
    """A setup step that fails raises RuntimeError, never an AssertionError."""
    if not ok:
        raise RuntimeError('setup failed: ' + what)


def _parent(name):
    monomer, end_groups = MONOMERS[name]
    return Polymer(label='X', monomer=monomer, end_groups=end_groups, cutoff=3, Mn=5000.0, Mw=6000.0,
                   initial_mass=1.0)


def _unit_daughters(parent):
    """One born-at-zero daughter per H-loss feature unit of the monomer."""
    units = parent._h_loss_feature_units()
    _setup(units, 'no H-loss feature units')
    return [parent._born_at_zero_mod_daughter(u, source=SOURCE) for u in units]


def _classes(daughters):
    """The isomorphism classes of the daughters' feature proxies."""
    classes = []
    for d in daughters:
        for c in classes:
            if c[0].feature_proxy.is_isomorphic(d.feature_proxy):
                c.append(d)
                break
        else:
            classes.append([d])
    return classes


def _pool_count(name):
    """(registered pools, isomorphism classes) over every H-loss unit."""
    parent = _parent(name)
    model = CoreEdgeReactionModel()
    daughters = _unit_daughters(parent)
    registered = [model._register_polymer(d, generate_thermo=False)[0] for d in daughters]
    return len({id(p) for p in registered}), len(_classes(daughters))


@RED
def test_a_monomer_with_one_site_class_registers_one_pool():
    """Every H of an ethylene-like unit is one site: one pool, never a split
    (the H lost at either stitch atom read as two radical-site labels)."""
    assert _pool_count('ethylene') == (1, 1)


@RED
def test_methacrylate_like_sites_register_one_pool_each():
    """A methyl on the ester oxygen and the methyl on the backbone carbon
    read alike (C, one radical, 2 H, one heavy neighbour): two structures."""
    n_pools, n_classes = _pool_count('methacrylate')
    assert n_pools == n_classes, (n_pools, n_classes)


@RED
def test_styrene_like_ring_sites_register_one_pool_each():
    """The ortho, meta and para ring H losses read alike: three structures."""
    n_pools, n_classes = _pool_count('styrene')
    assert n_pools == n_classes, (n_pools, n_classes)


@RED
@pytest.mark.parametrize('name', ['methacrylate', 'styrene'])
def test_each_registered_pool_carries_the_structure_it_was_registered_for(name):
    """The pool a daughter registers to has that daughter's feature proxy."""
    parent = _parent(name)
    model = CoreEdgeReactionModel()
    for d in _unit_daughters(parent):
        pool, _ = model._register_polymer(d, generate_thermo=False)
        assert pool.feature_proxy.is_isomorphic(d.feature_proxy), (
            d.feature_proxy.molecule[0].to_smiles(), pool.feature_proxy.molecule[0].to_smiles())


@pytest.mark.parametrize('name', [
    'methacrylate', 'styrene',
    pytest.param('ethylene', marks=RED)])
def test_twin_units_share_one_fingerprint(name):
    """Units whose feature proxies are isomorphic share one fingerprint."""
    daughters = _unit_daughters(_parent(name))
    for a, b in itertools.combinations(daughters, 2):
        same = a.feature_proxy.is_isomorphic(b.feature_proxy)
        if same:
            assert a.fingerprint == b.fingerprint


@RED
@pytest.mark.parametrize('name', ['methacrylate', 'styrene'])
def test_non_isomorphic_units_have_distinct_fingerprints(name):
    daughters = _unit_daughters(_parent(name))
    for a, b in itertools.combinations(daughters, 2):
        if not a.feature_proxy.is_isomorphic(b.feature_proxy):
            assert a.fingerprint != b.fingerprint, (a.feature_proxy.molecule[0].to_smiles(),
                                                    b.feature_proxy.molecule[0].to_smiles())


# --------------------------------------------------------------------------
# at the generation seam: one H abstraction per H site of the proxy
# --------------------------------------------------------------------------
def _tagged_proxy(parent):
    """The capped trimer with each atom's unit position (0..2, None for an
    end group) in props."""
    units = []
    for k in range(3):
        u = parent.monomer.copy(deep=True)
        for a in u.atoms:
            a.props = dict(a.props or {}, unit=k)
        units.append(u)
    e0, e1 = [eg.copy(deep=True) for eg in parent.end_groups]
    for a in e0.atoms + e1.atoms:
        a.props = dict(a.props or {}, unit=None)
    mol = e0
    for u in units:
        mol = stitch_molecules_by_labeled_atoms(mol, u)
    mol = stitch_molecules_by_labeled_atoms(mol, e1)
    _setup(mol is not None, 'trimer stitch')
    mol.clear_labeled_atoms()
    mol.update()
    return mol


def _h_abstraction_rows(name):
    """[(unit position, raw product, the row's polymer product)] for one H
    per heavy atom of the proxy, abstracted by CH3 on ONE model, so daughters
    register as in a run."""
    parent = _parent(name)
    model = CoreEdgeReactionModel()
    pool, _ = model._register_polymer(parent.copy(deep=True), generate_thermo=False)
    proxy = _tagged_proxy(parent)
    rows = []
    for i, atom in enumerate(proxy.atoms):
        if atom.is_hydrogen() or not any(b.is_hydrogen() for b in atom.bonds):
            continue
        product = proxy.copy(deep=True)
        x = product.atoms[i]
        product.remove_atom(next(b for b in x.bonds if b.is_hydrogen()))
        x.increment_radical()
        for a in product.atoms:
            a.props = {}
        product.update()
        raw = TemplateReaction(reactants=[pool.molecule[0].copy(deep=True), Molecule().from_smiles('[CH3]')],
                               products=[product.copy(deep=True), Molecule().from_smiles('C')],
                               family='H_Abstraction', template=['a'], is_forward=True, degeneracy=1,
                               kinetics=KIN)
        rxn, _ = model.make_new_reaction(raw, check_existing=False, generate_thermo=False,
                                         generate_kinetics=False)
        _setup(rxn is not None, 'H abstraction row')
        booked = [p for p in rxn.products if isinstance(p, Polymer)]
        _setup(len(booked) == 1, 'one polymer product')
        rows.append((atom.props['unit'], product, booked[0]))
    return parent, rows


def _renders(daughter, product):
    """True if ``product`` is one of the daughter's positional renderings
    (its feature unit at the center or either flank)."""
    spc = Polymer._h_loss_normalized_product_species(product)
    unit = daughter.feature_monomer
    return any(v is not None and v.is_isomorphic(spc)
               for v in (daughter._h_loss_positional_species(unit, p) for p in (1, 0, 2)))


@RED
@pytest.mark.parametrize('name', ['methacrylate', 'styrene'])
def test_a_center_h_loss_row_books_the_structure_that_reacted(name):
    """Every H lost from the center unit: the row's daughter renders, at the
    center, the row's own raw product."""
    _, rows = _h_abstraction_rows(name)
    centre = [(product, d) for unit, product, d in rows if unit == 1]
    _setup(centre, 'center-unit rows')
    for product, d in centre:
        assert d.feature_proxy.is_isomorphic(Polymer._h_loss_normalized_product_species(product)), (
            product.to_smiles(), d.label, d.feature_proxy.molecule[0].to_smiles())


@pytest.mark.parametrize('name', sorted(MONOMERS))
def test_positional_twins_book_one_pool(name):
    """A unit-0 H loss that renders a center site's daughter at a flank is
    that site read at another position: it books the center row's pool."""
    _, rows = _h_abstraction_rows(name)
    centre = {id(d): d for unit, _, d in rows if unit == 1}
    flank = [(product, d) for unit, product, d in rows if unit == 0]
    for product, d in flank:
        twin = [c for c in centre.values() if _renders(c, product)]
        if twin:
            assert any(c is d for c in twin), (product.to_smiles(), d.label)
