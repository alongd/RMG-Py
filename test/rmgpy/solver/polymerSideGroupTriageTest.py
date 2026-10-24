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
data only): three symptoms a methacrylate deck printed, and the double-count
ruling that closing the first one needs, each exercised through the
production seams (the Polymer producer, model registration, the enlarge
boundary, the solver rebuild). A setup step that fails raises RuntimeError,
never an AssertionError.

- symptom 2: a spawned scission tail carries its spawn markers, so the
  item-16 enlarge-boundary pass promotes it and its row runs live; a spawned
  pool the pass does not reach is refused naming that cause;
- the ruling, by channel identity (polymerChannelIdentityTest has the
  table): a generated backbone homolysis beside the pool's QSSA kernel is
  refused superseded and counted; a row the kernel does not carry stays live;
- symptom 3: a born-empty spawned pool starts at zero without a "no initial
  conditions" warning; a deck pool without one is still warned;
- symptom 4: neither a hydrogen-abstraction volatile-ejection row nor a
  concerted chain-end ejection is a channel the QSSA kernel carries: both
  stay live and no double-count census counts them.
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
    """A setup step that fails raises RuntimeError, never an AssertionError."""
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


def _spawn_tail(monkeypatch, caplog, channel=None, marked=True, reaction_channel=None):
    """A scission tail of the parent, registered, its row placed on the edge,
    then the enlarge boundary and the production rebuild. ``marked=False``
    strips the spawn markers: a daughter the item-16 promotion does not reach.
    ``reaction_channel`` is the stamp generation writes from the raw
    molecules (the tail proxy here is a capped stand-in, not the radical
    piece a homolysis makes); for a concerted stamp the small product is
    closed-shell, so the row's species read its stamp (I-255)."""
    import rmgpy.solver.polymer as sp_mod
    from rmgpy.polymer import Polymer, stitch_molecules_by_labeled_atoms
    parent = _parent(channel)
    ch3 = _spc('[CH3]', 'CH3')
    gas = [ch3]
    small = ch3
    if reaction_channel == 'concerted_molecular':
        small = _spc('C', 'CH4')
        gas.append(small)
    if channel is not None:
        mma = _spc(MONOMER_CLOSED, 'MMA')
        gas.append(mma)
        parent.monomer_product_species = mma
    model = _model_with_parent(parent, *gas)
    fragment = stitch_molecules_by_labeled_atoms(parent._stitch_wing('head'), _methyl_radical('*2'))
    tail = parent.create_reacted_copy(fragment)
    _setup(isinstance(tail, Polymer) and tail.label == 'PMMA_scission_tail', 'scission tail producer')
    if not marked:
        tail.parent_pool_label = None
        tail.spawn_metadata = None
    tail, is_new = model.make_new_species(tail, generate_thermo=False)
    _setup(is_new, 'tail registration')
    rxn = Reaction(reactants=[parent], products=[tail, small], **KIN_UNI)
    rxn.polymer_flux_archetype = sp_mod.FLUX_SCISSION_FRAGMENT
    rxn.polymer_reaction_channel = reaction_channel
    _place_edge_row(model, rxn)
    with caplog.at_level(logging.WARNING):
        _enlarge_boundary(model, monkeypatch)
        solver = _rebuild(model, parent, ch3, channel=channel)
    return model, solver, tail, rxn, _row(model, rxn)


class TestScissionTailPromotion:
    """Symptom 2. The scission producer stamps ``<label>_scission_tail`` with
    ``parent_pool_label`` and ``spawn_metadata``, so the item-16 pass promotes
    it, the rebuild derives its pool config, and its stamped row runs live. A
    spawned pool the pass does not reach (no markers) still has its row
    refused, and the refusal names that cause."""

    def test_tail_is_promoted_configured_and_live(self, monkeypatch, caplog):
        from rmgpy.polymer import is_engine_spawned_pool_daughter
        model, solver, tail, _, r = _spawn_tail(monkeypatch, caplog)
        assert tail.parent_pool_label == 'PMMA'
        assert tail.spawn_metadata == {'source': 'scission'}
        assert is_engine_spawned_pool_daughter(tail)
        assert tail in model.core.species
        assert 'PMMA_scission_tail' in {p.label for p in solver.polymer_pools}
        assert solver.reaction_refused[r] == 0
        assert not solver.spawned_pool_refusal_census
        assert solver.kernel_supersession_census == []

    def test_unconfigured_spawned_tail_is_refused_naming_the_cause(self, monkeypatch, caplog):
        model, solver, tail, _, r = _spawn_tail(monkeypatch, caplog, marked=False)
        _setup(tail not in model.core.species, 'an unmarked tail is not promoted')
        assert solver.reaction_refused[r] == 1
        assert [c['pools'] for c in solver.spawned_pool_refusal_census] == [['PMMA_scission_tail']]
        refusals = [rec.getMessage() for rec in caplog.records
                    if 'SPAWNED-POOL DEMOTION REFUSAL' in rec.getMessage()]
        assert len(refusals) == 1
        assert 'no solver pool config' in refusals[0]
        assert 'item-16 pending' not in refusals[0]


class TestKernelSupersessionByChannel:
    """The double-count ruling, by channel identity. A generated backbone
    homolysis cutting a pool whose QSSA kernel carries that initiation is
    refused by name (zero flux) and counted; its tail, whose only row that
    is, stays on the edge. A row whose channel the kernel does not carry
    stays live, whatever its archetype."""

    def test_a_backbone_homolysis_beside_qssa_is_refused_by_name(self, monkeypatch, caplog):
        import rmgpy.solver.polymer as sp_mod
        from rmgpy.polymer import _kernel_supersession_warned
        _kernel_supersession_warned.clear()
        model, solver, tail, rxn, r = _spawn_tail(monkeypatch, caplog, channel=RADICAL_CHANNEL,
                                                  reaction_channel='bond_fission')
        assert tail not in model.core.species
        assert 'PMMA_scission_tail' not in {p.label for p in solver.polymer_pools}
        assert not solver.spawned_pool_refusal_census
        assert solver.reaction_refused[r] == 1
        assert solver.reaction_flux_archetype[r] == sp_mod.FLUX_NONE
        assert solver.kernel_supersession_census == [
            {'reaction': str(rxn), 'pool': 'PMMA', 'channel': 'backbone_bond_fission',
             'kernels': ['radical_qssa_unzip'], 'reason': 'superseded:backbone_bond_fission'}]
        assert solver.kernel_double_count_census == []
        assert solver.qssa_double_count_census == []
        lines = [rec.getMessage() for rec in caplog.records
                 if 'KERNEL SUPERSESSION CENSUS' in rec.getMessage()]
        assert len(lines) == 1 and "'PMMA'" in lines[0] and str(rxn) in lines[0]

    def test_a_scission_row_qssa_does_not_carry_stays_live(self, monkeypatch, caplog):
        model, solver, tail, rxn, r = _spawn_tail(monkeypatch, caplog, channel=RADICAL_CHANNEL,
                                                  reaction_channel='concerted_molecular')
        assert tail in model.core.species
        assert solver.reaction_refused[r] == 0
        assert solver.kernel_supersession_census == []
        assert solver.kernel_double_count_census == []

    def test_the_kernel_channel_table(self):
        import rmgpy.solver.polymer as sp_mod
        from rmgpy.polymer import kernels_carrying_channel

        def pool(**kw):
            return sp_mod.PolymerPoolConfig(label='P', xs=3, explicit_dp_to_species_index={},
                                            mu_indices=(0, 1, 2), **kw)

        trip = {'A': 1.0e15, 'n': 0.0, 'Ea': 3.0e5}
        backbone, beta, side = 'backbone_bond_fission', 'radical_beta_scission', 'side_group_loss'
        sgh = [{'label': 'aliphatic_C-Br', 'gas_product': '[Br]'}]
        br, cl = _spc('[Br]', 'Br'), _spc('[Cl]', 'Cl')
        for channel in (backbone, beta, side, None):
            assert kernels_carrying_channel(pool(), channel, [br]) == ()
        assert kernels_carrying_channel(pool(radical_qssa_unzip=RADICAL_CHANNEL), backbone) == ('radical_qssa_unzip',)
        assert kernels_carrying_channel(pool(radical_qssa_unzip=RADICAL_CHANNEL), beta) == ('radical_qssa_unzip',)
        assert kernels_carrying_channel(pool(k_homolysis=trip), backbone) == ('k_homolysis',)
        assert kernels_carrying_channel(pool(k_homolysis=trip), beta) == ()
        assert kernels_carrying_channel(pool(kernel_homolysis=[{'template': 't'}]), backbone) == ('kernel_homolysis',)
        assert kernels_carrying_channel(pool(kernel_homolysis=[]), backbone) == ()
        assert kernels_carrying_channel(pool(k_scission=1.0), backbone) == ('k_scission',)
        assert kernels_carrying_channel(pool(k_scission=0.0), backbone) == ()
        assert kernels_carrying_channel(pool(k_depropagation=trip), beta) == ('k_depropagation',)
        assert kernels_carrying_channel(pool(k_depropagation=trip), backbone) == ()
        assert kernels_carrying_channel(pool(k_unzip=0.5), beta) == ('k_unzip',)
        assert kernels_carrying_channel(pool(side_group_homolysis=sgh), side, [br]) == ('side_group_homolysis',)
        assert kernels_carrying_channel(pool(side_group_homolysis=sgh), side, [cl]) == ()
        assert kernels_carrying_channel(pool(side_group_homolysis=sgh), backbone, [br]) == ()


class TestBornEmptyPoolInitialConditions:
    """Symptom 3. A spawned daughter pool is born empty by design (moments
    0, 0, 0; the parent keeps its mass), so it has no entry in the deck's
    initial moments; set_initial_conditions does not warn for it. A deck pool
    without initial conditions is still warned."""

    def _spawn(self, monkeypatch, caplog):
        import rmgpy.solver.polymer as sp_mod
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

    def test_born_empty_spawned_pool_starts_at_zero_unwarned(self, monkeypatch, caplog):
        solver, warned = self._spawn(monkeypatch, caplog)
        pools = {p.label: p for p in solver.polymer_pools}
        _setup(set(pools) == {'PMMA', 'PMMA_mod'}, 'daughter configured')
        assert pools['PMMA_mod'].engine_spawned and not pools['PMMA'].engine_spawned
        assert [float(solver.y0[i]) for i in pools['PMMA_mod'].mu_indices] == [0.0, 0.0, 0.0]
        assert [float(solver.y0[i]) for i in pools['PMMA'].mu_indices] == pytest.approx([0.01, 0.5, 30.0])
        assert warned == []

    def test_configured_pool_without_initial_conditions_still_warns(self, caplog):
        """Two deck-configured pools, one without initial conditions: that one
        is warned (a deck omission, not a spawned daughter)."""
        from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig
        core = [_spc('N#N', 'N2'),
                _spc('CO', 'A_mu0'), _spc('C=O', 'A_mu1'), _spc('C#N', 'A_mu2'),
                _spc('CCO', 'B_mu0'), _spc('CC=O', 'B_mu1'), _spc('CC#N', 'B_mu2')]
        pools = [PolymerPoolConfig(label=lbl, xs=3, explicit_dp_to_species_index={},
                                   mu_indices=mu, monomer_mw_g_mol=100.0)
                 for lbl, mu in (('A', (1, 2, 3)), ('B', (4, 5, 6)))]
        solver = HybridPolymerSystem(
            T=700.0, P=1.0e5, initial_mole_fractions={core[0]: 1.0}, V_poly=1.0,
            polymer_pools=pools, mass_transfer=[],
            gas_species_mask=np.array([True] + [False] * 6), constant_gas_volume=False,
            initial_polymer_moments={'A': (0.01, 0.5, 30.0)}, termination=[])
        with caplog.at_level(logging.WARNING):
            solver.initialize_model(core, [], [], [])
        assert [rec.getMessage() for rec in caplog.records
                if 'has no initial conditions provided' in rec.getMessage()] == [
            "Polymer pool 'B' has no initial conditions provided (Explicit or Moments)."]


class TestQssaDoubleCountOnChainMovingRows:
    """Symptom 4. Neither a hydrogen abstraction into the pool's
    radical-feature daughter (a bimolecular row: one H leaves, chains move
    intact) nor a same-pool concerted chain-end ejection of one repeat unit
    (closed-shell into closed-shell) is a channel the QSSA kernel carries
    (random backbone homolysis, radical depropagation): both stay live, and
    no double-count census or tripwire counts them."""

    def _build(self, monkeypatch, caplog, row):
        """``row``: the hydrogen abstraction into the radical-feature daughter
        (one H leaves, chains move intact), or a same-pool chain-end ejection
        of one repeat unit (closed-shell pieces: concerted)."""
        import rmgpy.solver.polymer as sp_mod
        parent = _parent(RADICAL_CHANNEL)
        mma = _spc(MONOMER_CLOSED, 'MMA')
        ch3 = _spc('[CH3]', 'CH3')
        model = _model_with_parent(parent, mma, ch3)
        parent.monomer_product_species = mma
        if row == 'h_abstraction':
            daughter, _ = model.make_new_species(_ester_h_loss_daughter(parent), generate_thermo=False)
            ch4 = _spc('C', 'CH4')
            rxn = Reaction(reactants=[ch3, parent], products=[ch4, daughter], **KIN_BI)
            rxn.polymer_eject_units = ((Molecule(smiles='C').get_molecular_weight()
                                        - Molecule(smiles='[CH3]').get_molecular_weight()) * 1000.0
                                       / parent.monomer_mw_g_mol)
        else:
            rxn = Reaction(reactants=[parent], products=[parent, mma], **KIN_UNI)
            rxn.polymer_eject_units = 1.0
        rxn.polymer_flux_archetype = sp_mod.FLUX_VOLATILE_EJECTION
        _place_edge_row(model, rxn)
        with caplog.at_level(logging.WARNING):
            _enlarge_boundary(model, monkeypatch)
            solver = _rebuild(model, parent, ch3, channel=RADICAL_CHANNEL)
        r = _row(model, rxn)
        _setup(solver.reaction_refused[r] == 0
               and solver.reaction_flux_archetype[r] == sp_mod.FLUX_VOLATILE_EJECTION,
               'live volatile-ejection row')
        return solver, r

    def test_a_concerted_chain_end_ejection_is_not_counted(self, monkeypatch, caplog):
        solver, r = self._build(monkeypatch, caplog, row='chain_end')
        _setup(solver.reaction_eject_units_moment[r] == 1.0, 'one repeat unit leaves')
        assert solver.qssa_double_count_census == []
        assert solver.kernel_double_count_census == []
        assert solver.kernel_supersession_census == []
        assert not any('DOUBLE-COUNT' in rec.getMessage() for rec in caplog.records)

    def test_hydrogen_abstraction_row_is_not_counted(self, monkeypatch, caplog):
        solver, r = self._build(monkeypatch, caplog, row='h_abstraction')
        _setup(solver.reaction_eject_units_moment[r] == 0.0, 'chains move intact')
        assert solver.qssa_double_count_census == []
        assert solver.kernel_double_count_census == []
        assert not any('DOUBLE-COUNT' in rec.getMessage() for rec in caplog.records)
