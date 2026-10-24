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
Side-group losses of a polymer pool, read from structure.

A side group is the far side of a single, acyclic bond of the repeat unit that
holds no stitch atom; the leaving group is that side with one radical at its
root. The side-group kernel takes any such leaving group whose element and
mass ledger closes (the unit minus the group, plus the group, is the unit),
and a carrier whose chains another kernel moves carries its losses as an
X-loss feature daughter with the exact per-chain defect M(group).

At the generation seam, a product that is the pool's proxy with one side
group cut off (the backbone whole, the radical at the cut) is the X-loss
feature with the exact defect, never a half-length scission tail; an H lost
from an end unit books its defect like any other unit's.

The parents are styrene- and methyl methacrylate-like test proxies (test data
only). A setup step that fails raises RuntimeError, never an AssertionError.
"""

import pytest

from rmgpy.molecule import Molecule
from rmgpy.polymer import Polymer
from rmgpy.species import Species

MONOMERS = {
    'PS': '[CH2][CH]c1ccccc1',
    'PMMA': '[CH2][C](C)C(=O)OC',
}

# (parent, channel label, leaving-group radical, site selector)
LEAVING_GROUPS = [
    ('PMMA', 'ester', 'CO[C]=O', 'aliphatic'),
    ('PMMA', 'methoxy', 'C[O]', 'aliphatic'),
    ('PMMA', 'alpha_methyl', '[CH3]', 'aliphatic'),
    ('PMMA', 'ester_methyl', '[CH3]', 'hetero'),
    ('PS', 'phenyl', '[c]1ccccc1', 'aliphatic'),
]

# Homolysis products of the trimer proxies (end groups CH3 / H) that lose one
# side group and keep the whole backbone, with the radical at the cut:
# (parent, product, leaving group).
SIDE_LOSS_PRODUCTS = [
    ('PMMA', 'COC(=O)C(C)(C)C[C](C)CC(C)(C)C(=O)OC', 'CO[C]=O'),         # center ester
    ('PMMA', 'COC(=O)[C](CC(C)(C)C(=O)OC)CC(C)(C)C(=O)OC', '[CH3]'),     # center alpha-methyl
    ('PMMA', 'COC(=O)C(C)(C)CC(C)(CC(C)(C)C(=O)OC)C([O])=O', '[CH3]'),   # center ester methyl
    ('PMMA', 'COC(=O)C(C)(C)CC(C)([C]=O)CC(C)(C)C(=O)OC', 'C[O]'),       # center methoxy
    ('PMMA', 'COC(=O)C(C)(C)CC(C)(C[C](C)C)C(=O)OC', 'CO[C]=O'),         # flank ester
    ('PS', 'CC(C[CH]CC(C)c1ccccc1)c1ccccc1', '[c]1ccccc1'),              # center phenyl
    ('PS', 'C[CH]CC(CC(C)c1ccccc1)c1ccccc1', '[c]1ccccc1'),              # flank phenyl
]

# An H lost from an end unit's side group (the H-abstraction shape): the
# product the producer folded into the parent.
END_UNIT_H_LOSS = [
    ('PMMA', '[CH2]OC(=O)C(C)(C)CC(C)(CC(C)(C)C(=O)OC)C(=O)OC'),         # unit 0 ester methyl
    ('PS', 'CC(CC(CC(C)c1ccccc1)c1ccccc1)c1cc[c]cc1'),                   # unit 0 para H
]


def _setup(ok, what):
    """A setup step that fails raises RuntimeError, never an AssertionError."""
    if not ok:
        raise RuntimeError('setup failed: ' + what)


def _parent(name, **kw):
    return Polymer(label=name, monomer=MONOMERS[name], end_groups=['[CH3]', '[H]'], cutoff=3,
                   Mn=5000.0, Mw=6000.0, initial_mass=1.0, **kw)


def _channel(label, gas, selector, sites=1.0):
    return dict(label=label, A=1.0e16, n=0.0, Ea=3.0e5, site_selector=selector,
                sites_per_unit=sites, gas_product=gas)


def _qssa():
    return dict(initiation=dict(A=1.0e15, n=0.0, Ea=3.0e5),
                depropagation=dict(A=1.0e13, n=0.0, Ea=8.0e4),
                termination=dict(A=1.0e8, n=0.0, Ea=1.0e4))


def _census(mol):
    out = {}
    for a in mol.atoms:
        out[a.symbol] = out.get(a.symbol, 0) + 1
    return out


def _mw(smiles_or_mol):
    mol = Molecule().from_smiles(smiles_or_mol) if isinstance(smiles_or_mol, str) else smiles_or_mol
    return mol.get_molecular_weight() * 1000.0


def _attempt(fn):
    """fn() or the ValueError it raised (a refusal is an outcome to assert
    on, not a setup failure)."""
    try:
        return fn()
    except ValueError as e:
        return e


# --- the kernel reads a polyatomic leaving group from structure ------------

@pytest.mark.parametrize('name,label,gas,selector', LEAVING_GROUPS)
def test_the_side_group_kernel_accepts_a_polyatomic_leaving_group(name, label, gas, selector):
    """validate_side_group_homolysis with the repeat unit resolves the group's
    one site per unit (sites_per_unit is checked against it)."""
    from rmgpy.solver.polymer import validate_side_group_homolysis
    poly = _parent(name)
    out = _attempt(lambda: validate_side_group_homolysis(name, [_channel(label, gas, selector)],
                                                         monomer=poly.monomer))
    assert not isinstance(out, Exception), f'{name} {label} {gas}: refused: {out}'
    assert [ch['gas_product'] for ch in out] == [gas]


@pytest.mark.parametrize('name,label,gas,selector', LEAVING_GROUPS)
def test_a_polyatomic_side_loss_unit_closes_its_ledger(name, label, gas, selector):
    """The X-loss feature unit is the repeat unit minus the group: its
    elements plus the group's are the unit's, its mass is the unit's less
    M(group), and it carries one defect radical beside the two stitch
    sites."""
    poly = _parent(name)
    unit = _attempt(lambda: poly._side_loss_feature_unit(label, gas, selector))
    assert not isinstance(unit, Exception), f'{name} {label} {gas}: refused: {unit}'
    total = _census(unit)
    for el, n in _census(Molecule().from_smiles(gas)).items():
        total[el] = total.get(el, 0) + n
    assert total == _census(poly.monomer)
    assert _mw(poly.monomer) - _mw(unit) == pytest.approx(_mw(gas), rel=1e-12)
    assert unit.get_radical_count() == 3


def test_a_carrier_whose_chains_another_kernel_moves_gets_an_x_loss_daughter():
    """Beside an unzip kernel the inventory accounting (a stationary carrier)
    cannot hold, so each channel spawns its X-loss feature daughter carrying
    the exact per-chain defect M(group)."""
    poly = _parent('PMMA', radical_qssa_unzip=_qssa(),
                   side_group_homolysis=[_channel('ester', 'CO[C]=O', 'aliphatic')])
    daughters = _attempt(lambda: poly.generate_side_loss_daughters())
    assert not isinstance(daughters, Exception), f'refused: {daughters}'
    assert [d.label for d in daughters] == ['PMMA_sidegrp_ester']
    d = daughters[0]
    assert d.chain_mass_defect_g_mol == pytest.approx(_mw('CO[C]=O'), rel=1e-12)
    assert d.parent_pool_label == 'PMMA'
    assert (d.spawn_metadata or {}).get('source') == 'side_group_homolysis'
    want = poly._side_loss_feature_unit('ester', 'CO[C]=O', 'aliphatic')
    assert d.feature_monomer.is_isomorphic(want)


def test_a_stationary_carrier_spawns_no_daughter():
    """Guard (passes today): without another kernel moving its chains the
    carrier keeps the inventory accounting and spawns nothing."""
    poly = _parent('PMMA', side_group_homolysis=[_channel('ester', 'CO[C]=O', 'aliphatic')])
    assert poly.generate_side_loss_daughters() == ()


def test_the_deck_takes_a_polyatomic_side_group_beside_an_unzip_kernel():
    """The deck helper accepts the ester channel beside radical_qssa_unzip,
    registers the gas radical and the X-loss feature daughter."""
    from unittest.mock import MagicMock
    from rmgpy.rmg import input as rmg_input

    def _make_new_species(obj, **kwargs):
        if isinstance(obj, Species):
            return obj, True
        return Species(label=obj.to_smiles(), molecule=[obj]), True

    mock_rmg = MagicMock()
    mock_rmg.initial_species = []
    mock_rmg.reaction_model.iteration_num = 0
    mock_rmg.reaction_model.new_species_list = []
    mock_rmg.reaction_model.make_new_species.side_effect = _make_new_species
    old_rmg, old_sd = rmg_input.rmg, rmg_input.species_dict
    rmg_input.rmg, rmg_input.species_dict = mock_rmg, {}
    try:
        poly = _attempt(lambda: rmg_input.polymer(
            label='PMMA', monomer=MONOMERS['PMMA'], end_groups=['[CH3]', '[H]'], cutoff=3,
            Mn=5000.0, Mw=6000.0, initial_mass=0.001, monomer_product='C=C(C)C(=O)OC',
            radical_qssa_unzip=_qssa(),
            side_group_homolysis=[_channel('ester', 'CO[C]=O', 'aliphatic')]))
        sd = dict(rmg_input.species_dict)
    except Exception as e:  # InputError: the deck refused the pair
        poly, sd = e, {}
    finally:
        rmg_input.rmg, rmg_input.species_dict = old_rmg, old_sd
    assert not isinstance(poly, Exception), f'refused: {poly}'
    daughter = sd.get('PMMA_sidegrp_ester')
    assert isinstance(daughter, Polymer)
    assert daughter in mock_rmg.initial_species
    assert daughter.chain_mass_defect_g_mol == pytest.approx(_mw('CO[C]=O'), rel=1e-12)


# --- the producer routes a backbone-preserving loss to its X-loss feature ---

@pytest.mark.parametrize('name,product,gas', SIDE_LOSS_PRODUCTS)
def test_a_side_group_loss_is_the_x_loss_feature_never_a_scission_tail(name, product, gas):
    """create_reacted_copy of a product that keeps the whole backbone and
    lost one side group: a feature daughter with the exact per-chain defect,
    never a scission daughter."""
    parent = _parent(name)
    d = parent.create_reacted_copy(Molecule().from_smiles(product))
    source = (getattr(d, 'spawn_metadata', None) or {}).get('source')
    assert source != 'scission', f'{product}: booked as {d.label}'
    assert isinstance(d, Polymer) and source == 'side_group_loss', \
        f'{product}: {getattr(d, "label", d)} ({source})'
    assert d.chain_mass_defect_g_mol == pytest.approx(_mw(gas), rel=1e-12)


@pytest.mark.parametrize('name,product', END_UNIT_H_LOSS)
def test_an_end_unit_h_loss_books_its_defect(name, product):
    """An H lost from an end unit (handshake evidence: exactly one H) is the
    H-loss feature with the 1.008 g/mol defect, not the parent pool."""
    parent = _parent(name)
    d = parent.create_reacted_copy(Molecule().from_smiles(product), h_loss_feature=True)
    assert isinstance(d, Polymer), f'{product}: refused'
    assert d.label != name, f'{product}: folded into the parent'
    assert d.chain_mass_defect_g_mol == pytest.approx(_mw('[H]'), rel=1e-12)


def test_the_generated_side_loss_row_moves_no_repeat_units():
    """Through make_new_reaction: the ester homolysis of the center unit is a
    volatile ejection into the X-loss feature whose moment shift is zero (the
    chain keeps its length; its mass leaves as the defect), so it reads as a
    side-group loss, the side-group kernel's channel."""
    from rmgpy.data.kinetics.family import TemplateReaction
    from rmgpy.kinetics import Arrhenius
    from rmgpy.polymer import (CHANNEL_SIDE_GROUP_LOSS, kernel_carried_channel, reaction_channel,
                               volatile_ejection_moment_units)
    from rmgpy.rmg.model import CoreEdgeReactionModel
    from rmgpy.solver.polymer import FLUX_VOLATILE_EJECTION
    parent = _parent('PMMA')
    model = CoreEdgeReactionModel()
    model._register_polymer(parent, generate_thermo=False)
    pieces = [Molecule().from_smiles('COC(=O)C(C)(C)C[C](C)CC(C)(C)C(=O)OC'), Molecule().from_smiles('CO[C]=O')]
    raw = TemplateReaction(reactants=[parent.molecule[0].copy(deep=True)], products=pieces, family='test_family',
                           template=['a'], is_forward=True, degeneracy=1,
                           kinetics=Arrhenius(A=(1.0, 's^-1'), n=0.0, Ea=(0.0, 'kcal/mol'), T0=(298.15, 'K')))
    out, _ = model.make_new_reaction(raw, check_existing=False, generate_thermo=False, generate_kinetics=False)
    _setup(out is not None, 'row generated')
    pools = [p for p in out.products if isinstance(p, Polymer)]
    assert [(getattr(p, 'spawn_metadata', None) or {}).get('source') for p in pools] == ['side_group_loss']
    assert int(out.polymer_flux_archetype) == FLUX_VOLATILE_EJECTION and not out.polymer_refused
    a_mom = volatile_ejection_moment_units(float(out.polymer_eject_units), 0.0, pools[0].chain_mass_defect_g_mol,
                                           parent.monomer_mw_g_mol)
    assert a_mom == 0.0
    assert kernel_carried_channel(reaction_channel(out), a_mom != 0.0) == CHANNEL_SIDE_GROUP_LOSS
