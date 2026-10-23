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
Side-group triage on a methacrylate-like pool (methyl ester side group, test
data only): three of the symptoms a methacrylate deck prints, each reproduced
through the production seams (the Polymer producer, model registration, the
enlarge boundary, the solver rebuild).

Each symptom has a pin that passes today (the path is reached and the symptom
is what it is) and an expected-to-fail test naming it. The xfails are strict
and accept only an AssertionError; a setup that fails raises RuntimeError, so
an xfail can only pass by the symptom's own assertion failing.

- symptom 2: a spawned scission tail is never promoted by the item-16
  enlarge-boundary pass (it carries no spawn marker), so its row is refused
  "item-16 pending" on every rebuild;
- symptom 3: a born-empty spawned pool is warned "no initial conditions";
- symptom 4: a hydrogen-abstraction volatile-ejection row, which cuts no
  backbone bond, trips the QSSA double-count tripwire.
"""

import logging

import numpy as np
import pytest

from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species

MONOMER = '[CH2][C](C)C(=O)OC'
MONOMER_CLOSED = 'C=C(C)C(=O)OC'
RADICAL_CHANNEL = {
    'initiation': {'A': 1.0e16, 'n': 0.0, 'Ea': 300000.0},
    'depropagation': {'A': 1.0e13, 'n': 0.0, 'Ea': 100000.0},
    'termination': {'A': 1.0e7, 'n': 0.0, 'Ea': 0.0},
    'transfer': {'A': 1.0e-10, 'n': 0.0, 'Ea': 0.0},
    'efficiency': 1.0,
    'monomer_yield': 1.0,
    'basis': 'backbone_bonds_mu1_minus_mu0',
}
KIN_UNI = dict(kinetics=Arrhenius(A=(1.0, 's^-1'), n=0.0, Ea=(0.0, 'kcal/mol'), T0=(298.15, 'K')),
               reversible=False)
KIN_BI = dict(kinetics=Arrhenius(A=(2.0, 'm^3/(mol*s)'), n=0.0, Ea=(0.0, 'kcal/mol'), T0=(298.15, 'K')),
              reversible=False)


def _spc(smiles, label):
    s = Species(molecule=[Molecule().from_smiles(smiles)])
    s.label = label
    return s


def _setup(ok, what):
    """A setup step that fails must not satisfy an xfail(raises=AssertionError)."""
    if not ok:
        raise RuntimeError('setup failed: ' + what)


def _methyl_radical(label):
    return Molecule().from_adjacency_list(f"""multiplicity 2
        1 {label} C u1 p0 c0 {{2,S}} {{3,S}} {{4,S}}
        2 H u0 p0 c0 {{1,S}}
        3 H u0 p0 c0 {{1,S}}
        4 H u0 p0 c0 {{1,S}}""")


def _parent(channel=None):
    from rmgpy.polymer import Polymer
    return Polymer(label='PMMA', monomer=MONOMER, end_groups=['[CH3]', '[H]'], cutoff=3,
                   Mn=5000.0, Mw=6000.0, initial_mass=0.05,
                   radical_qssa_unzip=channel)


def _ester_h_loss_daughter(parent):
    """The radical-feature daughter of an H-loss on an ester methyl. Only the
    middle repeat unit's routes to ``PMMA_mod``: the producer books the same
    H-loss on either end unit as the parent pool itself."""
    from rmgpy.polymer import Polymer
    base = parent.baseline_proxy.molecule[0]
    for i, atom in enumerate(base.atoms):
        if not (atom.is_carbon() and any(n.is_oxygen() for n in atom.bonds)
                and sum(1 for n in atom.bonds if n.is_hydrogen()) == 3):
            continue
        mol = base.copy(deep=True)
        target = mol.atoms[i]
        mol.remove_atom(next(n for n in target.bonds if n.is_hydrogen()))
        target.increment_radical()
        mol.update()
        daughter = parent.create_reacted_copy(mol, h_loss_feature=True)
        if isinstance(daughter, Polymer) and daughter.label == 'PMMA_mod':
            return daughter
    _setup(False, 'no ester-methyl H-loss routes to PMMA_mod')


def _model_with_parent(parent, *gas):
    """The parent pool core-resident (proxy and moment dummies), plus core gas."""
    from rmgpy.rmg.model import CoreEdgeReactionModel
    model = CoreEdgeReactionModel()
    model._register_polymer(parent, generate_thermo=False)
    model.core.species.append(parent)
    for m in list(model.new_species_list):
        if getattr(m, 'is_moment_dummy', False):
            model.core.species.append(m)
    for g in gas:
        model.core.species.append(g)
    model.new_species_list = []
    return model


def _place_edge_row(model, rxn):
    for spec in rxn.products:
        if spec not in model.core.species and spec not in model.edge.species:
            model.add_species_to_edge(spec)
    model.add_reaction_to_edge(rxn)


def _row(model, rxn):
    """The solver row of ``rxn``: core rows first, then edge rows. A row moves
    to the core when its last edge participant is promoted."""
    rows = list(model.core.reactions) + list(model.edge.reactions)
    _setup(sum(1 for x in rows if x is rxn) == 1, 'row present once')
    return next(i for i, x in enumerate(rows) if x is rxn)


def _enlarge_boundary(model, monkeypatch):
    """Only the end-of-enlarge passes run (all react flags False)."""
    import rmgpy.rmg.model as model_mod

    class _NoForbidden:
        @staticmethod
        def is_molecule_forbidden(mol):
            return False

    monkeypatch.setattr(model_mod, 'get_db', lambda name: _NoForbidden())
    n = len(model.core.species)
    model.enlarge(react_edge=True, unimolecular_react=np.zeros(n, dtype=bool),
                  bimolecular_react=np.zeros((n, n), dtype=bool))


def _rebuild(model, parent, gas, channel=None):
    """The production rebuild over the live model lists; only the parent is
    deck-configured, every daughter config comes from the engine."""
    from rmgpy.quantity import Quantity
    from rmgpy.rmg.polymer_input import HybridPolymerReactor, PolymerPhase, PolymerPool
    mu = [next(s for s in model.core.species if s.label == f'PMMA_mu{k}') for k in (0, 1, 2)]
    pool = PolymerPool(label='PMMA', xs=3, monomer=Molecule().from_smiles(MONOMER_CLOSED),
                       explicit_map={}, mu_species=mu, proxy_species=parent,
                       radical_qssa_unzip=channel,
                       monomer_product=(None if channel is None else
                                        next(s for s in model.core.species if s.label == 'MMA')))
    phase = PolymerPhase(density=Quantity(1180.0, 'kg/m^3'),
                         initial_moments={'PMMA': (0.01, 0.5, 30.0)},
                         initial_explicit={}, pools=[pool])
    reactor = HybridPolymerReactor(temperature=(700.0, 'K'), pressure=(1.0e5, 'Pa'),
                                   initialMoles={gas: 1.0}, polymerPhase=phase)
    solver = reactor.to_solver_object(model.core.species, model.core.reactions,
                                      model.edge.species, model.edge.reactions)
    solver.initialize_model(model.core.species, model.core.reactions,
                            model.edge.species, model.edge.reactions)
    return solver


class TestScissionTailPromotion:
    """Symptom 2. The scission producer builds ``<label>_scission_tail`` with no
    ``parent_pool_label`` and no ``spawn_metadata``, so
    ``is_engine_spawned_pool_daughter`` is False and the item-16 pass never
    promotes it; the r91 rebuild then refuses its stamped row "item-16 pending",
    and the refused row's zero flux can never promote it either."""

    def _spawn(self, monkeypatch, caplog):
        import rmgpy.solver.polymer as sp_mod
        from rmgpy.polymer import Polymer, stitch_molecules_by_labeled_atoms
        parent = _parent()
        ch3 = _spc('[CH3]', 'CH3')
        model = _model_with_parent(parent, ch3)
        fragment = stitch_molecules_by_labeled_atoms(parent._stitch_wing('head'), _methyl_radical('*2'))
        tail = parent.create_reacted_copy(fragment)
        _setup(isinstance(tail, Polymer) and tail.label == 'PMMA_scission_tail', 'scission tail producer')
        tail, is_new = model.make_new_species(tail, generate_thermo=False)
        _setup(is_new, 'tail registration')
        rxn = Reaction(reactants=[parent], products=[tail, ch3], **KIN_UNI)
        rxn.polymer_flux_archetype = sp_mod.FLUX_SCISSION_FRAGMENT
        _place_edge_row(model, rxn)
        with caplog.at_level(logging.WARNING):
            _enlarge_boundary(model, monkeypatch)
            solver = _rebuild(model, parent, ch3)
        return model, solver, tail, _row(model, rxn)

    def test_tail_is_unmarked_unpromoted_and_refused_today(self, monkeypatch, caplog):
        from rmgpy.polymer import is_engine_spawned_pool_daughter
        model, solver, tail, r = self._spawn(monkeypatch, caplog)
        assert getattr(tail, 'parent_pool_label', None) is None
        assert not getattr(tail, 'spawn_metadata', None)
        assert not is_engine_spawned_pool_daughter(tail)
        assert tail not in model.core.species
        assert 'PMMA_scission_tail' not in {p.label for p in solver.polymer_pools}
        assert solver.reaction_refused[r] == 1
        assert [c['pools'] for c in solver.spawned_pool_refusal_census] == [['PMMA_scission_tail']]
        assert any('item-16 pending' in rec.getMessage() for rec in caplog.records)

    @pytest.mark.xfail(strict=True, raises=AssertionError,
                       reason='symptom 2: a spawned scission tail carries no spawn marker, is never '
                              'promoted by item 16, and its row stays refused "item-16 pending"')
    def test_tail_is_promoted_configured_and_live(self, monkeypatch, caplog):
        model, solver, tail, r = self._spawn(monkeypatch, caplog)
        assert tail in model.core.species
        assert 'PMMA_scission_tail' in {p.label for p in solver.polymer_pools}
        assert solver.reaction_refused[r] == 0
        assert not solver.spawned_pool_refusal_census


class TestBornEmptyPoolInitialConditions:
    """Symptom 3. A spawned daughter pool is born empty by design (moments
    0, 0, 0; the parent keeps its mass), so it has no entry in the deck's
    initial moments; set_initial_conditions warns "no initial conditions" for
    it on every rebuild while its moment slots start at exactly zero."""

    def _spawn(self, monkeypatch, caplog):
        import rmgpy.solver.polymer as sp_mod
        from rmgpy.polymer import Polymer
        parent = _parent()
        h = _spc('[H]', 'H')
        model = _model_with_parent(parent, h)
        daughter, is_new = model.make_new_species(_ester_h_loss_daughter(parent), generate_thermo=False)
        _setup(is_new, 'daughter registration')
        h2 = _spc('[H][H]', 'H2')
        rxn = Reaction(reactants=[h, parent], products=[h2, daughter], **KIN_BI)
        rxn.polymer_flux_archetype = sp_mod.FLUX_VOLATILE_EJECTION
        rxn.polymer_eject_units = ((Molecule(smiles='[H][H]').get_molecular_weight()
                                    - Molecule(smiles='[H]').get_molecular_weight()) * 1000.0
                                   / parent.monomer_mw_g_mol)
        _place_edge_row(model, rxn)
        with caplog.at_level(logging.WARNING):
            _enlarge_boundary(model, monkeypatch)
            solver = _rebuild(model, parent, h)
        warned = [rec.getMessage() for rec in caplog.records
                  if 'has no initial conditions provided' in rec.getMessage()]
        return solver, warned

    def test_born_empty_pool_is_warned_and_starts_at_zero(self, monkeypatch, caplog):
        solver, warned = self._spawn(monkeypatch, caplog)
        pools = {p.label: p for p in solver.polymer_pools}
        assert set(pools) == {'PMMA', 'PMMA_mod'}
        assert [float(solver.y0[i]) for i in pools['PMMA_mod'].mu_indices] == [0.0, 0.0, 0.0]
        assert [float(solver.y0[i]) for i in pools['PMMA'].mu_indices] == pytest.approx([0.01, 0.5, 30.0])
        assert warned == ["Polymer pool 'PMMA_mod' has no initial conditions provided (Explicit or Moments)."]

    @pytest.mark.xfail(strict=True, raises=AssertionError,
                       reason='symptom 3: a born-empty spawned pool is warned "no initial conditions" '
                              'although zero moments are its declared initial condition')
    def test_born_empty_spawned_pool_is_not_warned(self, monkeypatch, caplog):
        solver, warned = self._spawn(monkeypatch, caplog)
        _setup('PMMA_mod' in {p.label for p in solver.polymer_pools}, 'daughter configured')
        assert warned == []


class TestQssaDoubleCountOnHydrogenAbstraction:
    """Symptom 4. The QSSA tripwire counts any live SCISSION_FRAGMENT or
    VOLATILE_EJECTION row sourced from a QSSA pool as overlapping its random
    backbone homolysis. A hydrogen abstraction into the pool's radical-feature
    daughter (volatile ejection of one H) cuts no backbone bond, yet trips it:
    on the triage deck this row shape is the tripwire's only trigger."""

    def _build(self, monkeypatch, caplog):
        import rmgpy.solver.polymer as sp_mod
        from rmgpy.polymer import Polymer
        parent = _parent(RADICAL_CHANNEL)
        mma = _spc(MONOMER_CLOSED, 'MMA')
        ch3 = _spc('[CH3]', 'CH3')
        model = _model_with_parent(parent, mma, ch3)
        parent.monomer_product_species = mma
        daughter, _ = model.make_new_species(_ester_h_loss_daughter(parent), generate_thermo=False)
        ch4 = _spc('C', 'CH4')
        rxn = Reaction(reactants=[ch3, parent], products=[ch4, daughter], **KIN_BI)
        rxn.polymer_flux_archetype = sp_mod.FLUX_VOLATILE_EJECTION
        rxn.polymer_eject_units = ((Molecule(smiles='C').get_molecular_weight()
                                    - Molecule(smiles='[CH3]').get_molecular_weight()) * 1000.0
                                   / parent.monomer_mw_g_mol)
        _place_edge_row(model, rxn)
        with caplog.at_level(logging.WARNING):
            _enlarge_boundary(model, monkeypatch)
            solver = _rebuild(model, parent, ch3, channel=RADICAL_CHANNEL)
        r = _row(model, rxn)
        _setup(solver.reaction_refused[r] == 0
               and solver.reaction_flux_archetype[r] == sp_mod.FLUX_VOLATILE_EJECTION,
               'live hydrogen-abstraction row')
        return solver

    def test_hydrogen_abstraction_row_trips_the_tripwire_today(self, monkeypatch, caplog):
        solver = self._build(monkeypatch, caplog)
        assert solver.qssa_double_count_census == [{'pool': 'PMMA', 'overlap': 'generated_scission_ve'}]
        assert any('QSSA DOUBLE-COUNT TRIPWIRE' in rec.getMessage() for rec in caplog.records)

    @pytest.mark.xfail(strict=True, raises=AssertionError,
                       reason='symptom 4: a hydrogen-abstraction volatile-ejection row cuts no '
                              'backbone bond yet trips the QSSA double-count tripwire')
    def test_hydrogen_abstraction_row_does_not_trip(self, monkeypatch, caplog):
        solver = self._build(monkeypatch, caplog)
        assert solver.qssa_double_count_census == []
