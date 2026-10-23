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
Supersession by channel identity. A generated row is refused as superseded
only when a kernel on the pool it decomposes carries the same reaction
channel; every other row stays live and is not counted as a double count.

- a row's channel is read from its real molecules where it is generated
  (bond fission, radical beta-scission, concerted molecular, ...), in its
  decomposition direction, so it does not depend on how the row is written;
- the kernels carry: homolysis (k_homolysis, kernel_homolysis, k_scission)
  backbone bond fission; radical_qssa_unzip backbone bond fission (its
  initiation) and radical beta-scission; k_depropagation and k_unzip radical
  beta-scission; side_group_homolysis the fission that releases its declared
  gas product;
- one census flags a live row whose channel a kernel on its pool carries;
- a pool that never held more than its floor (born empty) is not announced
  as exhausted; a spawned daughter is promoted only once it has a live row.

The generated rows are built from a polystyrene-like test proxy (test data
only). A setup step that fails raises RuntimeError, never an AssertionError.
"""

import logging

import numpy as np
import pytest

import rmgpy.solver.polymer as sp_mod
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Bond, Molecule
from rmgpy.reaction import Reaction
from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig
from rmgpy.species import Species

from polymerKernelsTest import h_term

MONOMER = '[CH2][CH]c1ccccc1'
MONOMER_CLOSED = 'C=Cc1ccccc1'
QSSA = {
    'initiation': {'A': 1.0e15, 'n': 0.0, 'Ea': 3.0e5},
    'depropagation': {'A': 1.0e13, 'n': 0.0, 'Ea': 8.0e4},
    'termination': {'A': 1.0e8, 'n': 0.0, 'Ea': 1.0e4},
}
TRIPLET = {'A': 1.0e15, 'n': 0.0, 'Ea': 3.0e5}
KIN = Arrhenius(A=(1.0, 's^-1'), n=0.0, Ea=(0.0, 'kcal/mol'), T0=(298.15, 'K'))


def _spc(smiles, label):
    s = Species(molecule=[Molecule().from_smiles(smiles)])
    s.label = label
    return s


def _setup(ok, what):
    """A setup step that fails raises RuntimeError, never an AssertionError."""
    if not ok:
        raise RuntimeError('setup failed: ' + what)


def _stamp(rxn, channel):
    """The channel stamp generation writes (a hand-built row stands in for a
    generated one). Before the ruling the attribute does not exist and the
    row stays unstamped."""
    try:
        rxn.polymer_reaction_channel = channel
    except AttributeError:
        pass


def _census(rs, name):
    return getattr(rs, name, None)


# --------------------------------------------------------------------------
# the channel of a row, from its real molecules
# --------------------------------------------------------------------------
def _mols(*smiles):
    return [Molecule().from_smiles(s) for s in smiles]


CHANNEL_CASES = [
    # (reactants, products, channel)
    (('CCCC',), ('[CH2]C', '[CH2]C'), 'bond_fission'),
    (('[CH2]C', '[CH2]C'), ('CCCC',), 'bond_fission'),
    (('CC[CH]CC',), ('C=CCC', '[CH3]'), 'radical_beta_scission'),
    (('C=CCC', '[CH3]'), ('CC[CH]CC',), 'radical_beta_scission'),
    (('C=CCCCC',), ('C=CC', 'C=CC'), 'concerted_molecular'),
    (('[CH2]CCCC',), ('CC[CH]CC',), 'isomerization'),
    (('C', '[CH3]'), ('[CH3]', 'C'), 'bimolecular'),
    (('CCCC',), ('[CH2]CC', 'C'), 'unidentified'),
]


@pytest.mark.parametrize('reactants,products,channel', CHANNEL_CASES)
def test_the_channel_is_read_from_the_molecules_in_the_decomposition_direction(reactants, products, channel):
    from rmgpy.polymer import classify_reaction_channel
    assert classify_reaction_channel(_mols(*reactants), _mols(*products)) == channel


# --------------------------------------------------------------------------
# generated rows: the production generation seam on a styrene-like proxy
# --------------------------------------------------------------------------
def _parent(channel=None):
    from rmgpy.polymer import Polymer
    return Polymer(label='PS', monomer=MONOMER, end_groups=['[CH3]', '[H]'], cutoff=3,
                   Mn=5000.0, Mw=6000.0, initial_mass=0.05, radical_qssa_unzip=channel)


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


def _backbone_cuts(base):
    for i, a in enumerate(base.atoms):
        for b in list(a.bonds):
            j = base.atoms.index(b)
            bond = base.get_bond(a, b)
            if (j > i and a.is_carbon() and b.is_carbon() and bond.is_single()
                    and not base.is_bond_in_cycle(bond)):
                yield i, j


def _split(base, i, j, concerted):
    """Cut bond i-j of the proxy: as a homolysis (two radicals), or as a
    concerted split (an H from a carbon neighbour of i moves onto j, which
    closes a C=C at i; both pieces closed-shell)."""
    m = base.copy(deep=True)
    ai, aj = m.atoms[i], m.atoms[j]
    m.remove_bond(m.get_bond(ai, aj))
    if concerted:
        n = next((x for x in ai.bonds if x.is_carbon() and any(h.is_hydrogen() for h in x.bonds)), None)
        if n is None:
            return None
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
    return pieces


def _raw_row(parent, gas_smiles, concerted):
    """The raw generated row (molecules, before the handshake) whose small
    product is ``gas_smiles``: the proxy cut at the bond that releases it."""
    ref = Species(molecule=[Molecule().from_smiles(gas_smiles)])
    ref.generate_resonance_structures()
    base = parent.molecule[0]
    for i, j in _backbone_cuts(base):
        pieces = _split(base, i, j, concerted)
        if pieces and len(pieces) == 2 and any(ref.is_isomorphic(Molecule().from_smiles(p.to_smiles()))
                                               for p in pieces):
            return TemplateReaction(reactants=[base.copy(deep=True)], products=list(pieces),
                                    family='test_family', template=['a'], is_forward=True,
                                    degeneracy=1, kinetics=KIN)
    _setup(False, 'no cut releases ' + gas_smiles)


HOMOLYSIS_GAS = '[CH2]C(C)c1ccccc1'      # one end unit leaves as a radical
CONCERTED_GAS = MONOMER_CLOSED           # one unit leaves as the closed monomer


def _generate(model, parent, gas_smiles, concerted):
    """The row through CoreEdgeReactionModel.make_new_reaction, then placed on
    the edge: parent -> scission tail + a small gas."""
    from rmgpy.polymer import Polymer
    out, _ = model.make_new_reaction(_raw_row(parent, gas_smiles, concerted), check_existing=False,
                                     generate_thermo=False, generate_kinetics=False)
    _setup(out is not None and out.reactants == [parent], 'row generated from the parent')
    tail = [p for p in out.products if isinstance(p, Polymer)]
    _setup(len(tail) == 1 and tail[0].label == 'PS_scission_tail', 'scission tail product')
    out.reversible = False   # the solver needs no thermo for an irreversible row
    _setup(out.polymer_flux_archetype == sp_mod.FLUX_VOLATILE_EJECTION and not out.polymer_refused,
           'live volatile-ejection row')
    for spec in out.products:
        if spec not in model.core.species and spec not in model.edge.species:
            model.add_species_to_edge(spec)
    model.add_reaction_to_edge(out)
    return out, tail[0]


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


def _rebuild(model, parent, monomer, channel=None):
    """The production rebuild over the live model lists; only the parent is
    deck-configured."""
    from rmgpy.quantity import Quantity
    from rmgpy.rmg.polymer_input import HybridPolymerReactor, PolymerPhase, PolymerPool
    mu = [next(s for s in model.core.species if s.label == f'PS_mu{k}') for k in (0, 1, 2)]
    pool = PolymerPool(label='PS', xs=3, monomer=Molecule().from_smiles(MONOMER_CLOSED),
                       explicit_map={}, mu_species=mu, proxy_species=parent,
                       radical_qssa_unzip=channel, monomer_product=monomer)
    phase = PolymerPhase(density=Quantity(1050.0, 'kg/m^3'),
                         initial_moments={'PS': (0.01, 0.5, 30.0)},
                         initial_explicit={}, pools=[pool])
    gas = next(s for s in model.core.species if s.label == 'N2')
    reactor = HybridPolymerReactor(temperature=(700.0, 'K'), pressure=(1.0e5, 'Pa'),
                                   initialMoles={gas: 1.0}, polymerPhase=phase)
    solver = reactor.to_solver_object(model.core.species, model.core.reactions,
                                      model.edge.species, model.edge.reactions)
    solver.initialize_model(model.core.species, model.core.reactions,
                            model.edge.species, model.edge.reactions)
    return solver


def _row(model, rxn):
    rows = list(model.core.reactions) + list(model.edge.reactions)
    _setup(sum(1 for x in rows if x is rxn) == 1, 'row present once')
    return next(i for i, x in enumerate(rows) if x is rxn)


def _qssa_deck(monkeypatch, caplog, gas_smiles, concerted, warm_qssa_key=True):
    from rmgpy.polymer import _qssa_double_count_warned
    if warm_qssa_key:
        _qssa_double_count_warned.discard(('PS', 'generated_scission_ve'))
    parent = _parent(QSSA)
    styrene = _spc(MONOMER_CLOSED, 'styrene')
    parent.monomer_product_species = styrene
    model = _model_with_parent(parent, _spc('N#N', 'N2'), styrene)
    rxn, tail = _generate(model, parent, gas_smiles, concerted)
    with caplog.at_level(logging.WARNING):
        _enlarge_boundary(model, monkeypatch)
        solver = _rebuild(model, parent, styrene, channel=QSSA)
    return model, solver, rxn, tail


class TestGeneratedRowsBesideTheQssaKernel:
    """The PS deck's case. Chain-end concerted rows (Retroene there; a
    concerted split here) make the scission tails beside a radical_qssa_unzip
    kernel. No radical kernel carries a concerted molecular channel, so the
    row stays live and no census counts it. A backbone homolysis is the
    kernel's own initiation: refused by name."""

    def test_a_generated_row_is_stamped_with_its_channel(self):
        stamped = {}
        for gas, concerted in ((HOMOLYSIS_GAS, False), (CONCERTED_GAS, True)):
            parent = _parent()
            model = _model_with_parent(parent, _spc('N#N', 'N2'))
            rxn, _ = _generate(model, parent, gas, concerted=concerted)
            stamped[concerted] = getattr(rxn, 'polymer_reaction_channel', None)
        assert stamped == {False: 'bond_fission', True: 'concerted_molecular'}

    def test_a_concerted_row_beside_qssa_stays_live_and_is_not_flagged(self, monkeypatch, caplog):
        model, solver, rxn, tail = _qssa_deck(monkeypatch, caplog, CONCERTED_GAS, concerted=True)
        r = _row(model, rxn)
        assert tail in model.core.species
        assert solver.reaction_refused[r] == 0
        assert solver.reaction_flux_archetype[r] == sp_mod.FLUX_VOLATILE_EJECTION
        assert solver.reaction_eject_units_moment[r] == pytest.approx(1.0)
        assert solver.qssa_double_count_census == []
        assert not any('DOUBLE-COUNT' in rec.getMessage() for rec in caplog.records)
        assert _census(solver, 'kernel_supersession_census') == []
        assert _census(solver, 'kernel_double_count_census') == []

    def test_a_backbone_homolysis_beside_qssa_is_refused_superseded(self, monkeypatch, caplog):
        model, solver, rxn, tail = _qssa_deck(monkeypatch, caplog, HOMOLYSIS_GAS, concerted=False)
        r = _row(model, rxn)
        assert solver.reaction_refused[r] == 1
        assert solver.reaction_flux_archetype[r] == sp_mod.FLUX_NONE
        assert _census(solver, 'kernel_supersession_census') == [
            {'reaction': str(rxn), 'pool': 'PS', 'channel': 'backbone_bond_fission',
             'kernels': ['radical_qssa_unzip'], 'reason': 'superseded:backbone_bond_fission'}]
        assert _census(solver, 'kernel_double_count_census') == []
        lines = [rec.getMessage() for rec in caplog.records
                 if 'KERNEL SUPERSESSION CENSUS' in rec.getMessage()]
        assert len(lines) == 1 and "'PS'" in lines[0] and str(rxn) in lines[0]


# --------------------------------------------------------------------------
# the rule, kernel by kernel (hand-built rows stamped as generation stamps)
# --------------------------------------------------------------------------
def _pool(label, mu0, **kw):
    return PolymerPoolConfig(label=label, xs=2, explicit_dp_to_species_index={},
                             mu_indices=(mu0, mu0 + 1, mu0 + 2), monomer_mw_g_mol=104.15, **kw)


def _kernel_system(kernel, channel, archetype=sp_mod.FLUX_VOLATILE_EJECTION, cls=HybridPolymerSystem,
                   t_proxy='CCCCCC'):
    """Pool A carries ``kernel``; one row decomposes A: A -> T + X (a
    volatile ejection moving one repeat unit) or, as a chip, A -> A + X.
    T, Re (a kernel_homolysis daughter) and the k_homolysis end-radical pools
    are born empty."""
    core = [_spc('N#N', 'N2'), _spc('[CH2]CC', 'X'), _spc('C=C', 'M')]
    labels = ['A', 'T', 'Re', 'A_rad_primary_end', 'A_rad_secondary_end']
    pools = []
    for n, label in enumerate(labels):
        stem = 'C' * (n + 1)
        base = len(core)
        proxy = t_proxy if label == 'T' else 'C' * (n + 5)
        core += [_spc(proxy, label), _spc(stem + 'O', f'{label}_mu0'),
                 _spc(stem + 'N', f'{label}_mu1'), _spc(stem + 'S', f'{label}_mu2')]
        kw = {}
        if label == 'A':
            kw = {'k_homolysis': dict(k_homolysis=TRIPLET),
                  'kernel_homolysis': dict(kernel_homolysis=[h_term(daughter='Re')]),
                  'k_scission': dict(k_scission=0.5),
                  'radical_qssa_unzip': dict(radical_qssa_unzip=QSSA, monomer_poly_index=2),
                  'k_depropagation': dict(k_depropagation=TRIPLET, monomer_poly_index=2),
                  'k_unzip': dict(k_unzip=0.5, monomer_poly_index=2),
                  None: {}}[kernel]
        pools.append(_pool(label, base + 1, **kw))
    sp = {s.label: s for s in core}
    if archetype == sp_mod.FLUX_DISCRETE_CHIP:
        rxn = Reaction(reactants=[sp['A']], products=[sp['A'], sp['X']], kinetics=KIN, reversible=False)
        rxn.polymer_chip_units = 1
    else:
        rxn = Reaction(reactants=[sp['A']], products=[sp['T'], sp['X']], kinetics=KIN, reversible=False)
        rxn.polymer_eject_units = 1.0
    rxn.polymer_flux_archetype = archetype
    if channel is not None:
        _stamp(rxn, channel)
    mask = np.array([s.label in ('N2', 'X', 'M') for s in core], dtype=bool)
    rs = cls(T=700.0, P=1.0e5, initial_mole_fractions={core[0]: 1.0}, V_poly=1.0,
             polymer_pools=pools, mass_transfer=[], gas_species_mask=mask, constant_gas_volume=False,
             initial_polymer_moments={'A': (1.0, 5.0, 30.0)}, termination=[])
    rs.initialize_model(core, [rxn], [], [])
    return rs, rxn


SUPERSEDED = [
    ('k_homolysis', 'bond_fission', 'backbone_bond_fission'),
    ('kernel_homolysis', 'bond_fission', 'backbone_bond_fission'),
    ('k_scission', 'bond_fission', 'backbone_bond_fission'),
    ('radical_qssa_unzip', 'bond_fission', 'backbone_bond_fission'),
    ('radical_qssa_unzip', 'radical_beta_scission', 'radical_beta_scission'),
    ('k_depropagation', 'radical_beta_scission', 'radical_beta_scission'),
    ('k_unzip', 'radical_beta_scission', 'radical_beta_scission'),
]

LIVE = [
    ('radical_qssa_unzip', 'concerted_molecular'),
    ('k_homolysis', 'radical_beta_scission'),
    ('k_homolysis', 'concerted_molecular'),
    ('k_depropagation', 'bond_fission'),
    ('k_scission', 'concerted_molecular'),
    ('k_scission', 'bimolecular'),
    (None, 'bond_fission'),
    (None, 'radical_beta_scission'),
]


@pytest.mark.parametrize('kernel,channel,carried', SUPERSEDED)
def test_a_row_whose_channel_a_kernel_on_its_pool_carries_is_superseded(kernel, channel, carried):
    rs, rxn = _kernel_system(kernel, channel)
    assert rs.reaction_refused[0] == 1
    assert rs.reaction_flux_archetype[0] == sp_mod.FLUX_NONE
    assert _census(rs, 'kernel_supersession_census') == [
        {'reaction': str(rxn), 'pool': 'A', 'channel': carried, 'kernels': [kernel],
         'reason': 'superseded:' + carried}]
    assert _census(rs, 'kernel_double_count_census') == []


@pytest.mark.parametrize('kernel,channel', LIVE)
def test_a_row_whose_channel_no_kernel_on_its_pool_carries_stays_live(kernel, channel):
    rs, rxn = _kernel_system(kernel, channel)
    assert rs.reaction_refused[0] == 0
    assert rs.reaction_flux_archetype[0] == sp_mod.FLUX_VOLATILE_EJECTION
    assert not _census(rs, 'kernel_supersession_census')
    assert not _census(rs, 'kernel_double_count_census')


def test_a_chip_row_follows_the_same_rule():
    """A discrete chip is a row like any other: a fission chip beside
    k_scission is superseded; a concerted chip stays live and the item-18
    tripwire no longer counts it."""
    rs_fission, rxn = _kernel_system('k_scission', 'bond_fission', archetype=sp_mod.FLUX_DISCRETE_CHIP)
    assert rs_fission.reaction_refused[0] == 1
    assert _census(rs_fission, 'kernel_supersession_census') == [
        {'reaction': str(rxn), 'pool': 'A', 'channel': 'backbone_bond_fission', 'kernels': ['k_scission'],
         'reason': 'superseded:backbone_bond_fission'}]
    rs_split, _ = _kernel_system('k_scission', 'concerted_molecular', archetype=sp_mod.FLUX_DISCRETE_CHIP)
    assert rs_split.reaction_refused[0] == 0
    assert rs_split.reaction_flux_archetype[0] == sp_mod.FLUX_DISCRETE_CHIP
    assert _census(rs_split, 'kernel_double_count_census') == []
    assert not getattr(rs_split, 'double_count_census', None)


def test_an_unstamped_row_is_read_from_its_own_species():
    """A row built without the generation stamp (a hand-built or unpickled
    row) is classified from its own species: a closed-shell proxy into two
    radicals is a bond fission, superseded beside k_scission."""
    rs, rxn = _kernel_system('k_scission', None, t_proxy='[CH2]CCCCC')
    assert rs.reaction_refused[0] == 1
    assert [c['channel'] for c in _census(rs, 'kernel_supersession_census') or []] == ['backbone_bond_fission']


def test_an_unstamped_row_that_is_no_fission_stays_live():
    """A closed T and one radical X: no channel a kernel carries."""
    rs, _ = _kernel_system('k_scission', None)
    assert rs.reaction_refused[0] == 0


class _NoSupersession(HybridPolymerSystem):
    """Injects a live duplicate: the supersession pass does nothing."""

    def _refuse_kernel_superseded_rows(self, *args, **kwargs):
        return None


def test_an_injected_live_duplicate_of_a_kernels_channel_is_flagged(caplog):
    from rmgpy import polymer
    getattr(polymer, '_kernel_double_count_warned', set()).clear()
    with caplog.at_level(logging.WARNING):
        rs, rxn = _kernel_system('radical_qssa_unzip', 'bond_fission', cls=_NoSupersession)
    _setup(rs.reaction_refused[0] == 0, 'the duplicate is live')
    assert _census(rs, 'kernel_double_count_census') == [
        {'reaction': str(rxn), 'pool': 'A', 'channel': 'backbone_bond_fission',
         'kernels': ['radical_qssa_unzip']}]
    lines = [rec.getMessage() for rec in caplog.records
             if 'KERNEL CHANNEL DOUBLE-COUNT CENSUS' in rec.getMessage()]
    assert len(lines) == 1 and "'A'" in lines[0] and str(rxn) in lines[0]


# --------------------------------------------------------------------------
# side-group loss: the fission that releases the declared gas product
# --------------------------------------------------------------------------
def _sgh_system(gas_product_row):
    """An SGH carrier 'poly' declaring [Br]; one fission row poly -> F + gas,
    F an X-loss feature pool (its per-chain defect is the gas's mass, so the
    chains move with no moment shift)."""
    channel = dict(label='aliphatic_C-Br', A=1.0e13, n=0.5, Ea=1.2e5, site_selector='aliphatic',
                   sites_per_unit=2.0, gas_product='[Br]')
    core = [_spc('N#N', 'N2'), _spc('[Br]', 'Br_rad'), _spc('[Cl]', 'Cl_rad'),
            _spc('CCCC', 'poly'), _spc('CO', 'poly_mu0'), _spc('C=O', 'poly_mu1'), _spc('C#N', 'poly_mu2'),
            _spc('CCCCC', 'F'), _spc('CCO', 'F_mu0'), _spc('CC=O', 'F_mu1'), _spc('CC#N', 'F_mu2')]
    gas = {'[Br]': core[1], '[Cl]': core[2]}[gas_product_row]
    m_x = gas.molecule[0].get_molecular_weight() * 1000.0
    pools = [PolymerPoolConfig(label='poly', xs=2, explicit_dp_to_species_index={}, mu_indices=(4, 5, 6),
                               monomer_mw_g_mol=106.949, side_group_homolysis=[channel],
                               side_group_gas_indices=(1,)),
             PolymerPoolConfig(label='F', xs=2, explicit_dp_to_species_index={}, mu_indices=(8, 9, 10),
                               monomer_mw_g_mol=106.949, chain_mass_defect_g_mol=m_x)]
    rxn = Reaction(reactants=[core[3]], products=[core[7], gas], kinetics=KIN, reversible=False)
    rxn.polymer_flux_archetype = sp_mod.FLUX_VOLATILE_EJECTION
    rxn.polymer_eject_units = m_x / 106.949
    _stamp(rxn, 'bond_fission')
    mask = np.array([s.label in ('N2', 'Br_rad', 'Cl_rad') for s in core], dtype=bool)
    rs = HybridPolymerSystem(T=800.0, P=1.0e5, initial_mole_fractions={core[0]: 1.0}, V_poly=1.0,
                             polymer_pools=pools, mass_transfer=[], gas_species_mask=mask,
                             constant_gas_volume=False,
                             initial_polymer_moments={'poly': (1.0, 5.0, 30.0)}, termination=[])
    rs.initialize_model(core, [rxn], [], [])
    return rs, rxn


def test_a_side_group_loss_of_the_declared_gas_is_superseded():
    """Before the ruling the row is live, so the stationary-carrier guard
    refuses the build (ValueError naming the carrier)."""
    rs, rxn = _sgh_system('[Br]')
    assert rs.reaction_eject_units_moment[0] == 0.0
    assert rs.reaction_refused[0] == 1
    assert _census(rs, 'kernel_supersession_census') == [
        {'reaction': str(rxn), 'pool': 'poly', 'channel': 'side_group_loss',
         'kernels': ['side_group_homolysis'], 'reason': 'superseded:side_group_loss'}]


def test_a_side_group_loss_of_another_gas_is_not_the_kernels_channel():
    """Not superseded, so the row still moves the carrier's chains, which the
    stationary-carrier guard refuses."""
    with pytest.raises(ValueError, match=r"poly"):
        _sgh_system('[Cl]')


# --------------------------------------------------------------------------
# follow-ups: the born-empty exhaustion line; promotion on a live row
# --------------------------------------------------------------------------
def _two_pool_exhaustion(caplog):
    core = [_spc('N#N', 'N2'),
            _spc('CO', 'A_mu0'), _spc('C=O', 'A_mu1'), _spc('C#N', 'A_mu2'),
            _spc('CCO', 'B_mu0'), _spc('CC=O', 'B_mu1'), _spc('CC#N', 'B_mu2')]
    pools = [PolymerPoolConfig(label=lbl, xs=3, explicit_dp_to_species_index={},
                               mu_indices=mu, monomer_mw_g_mol=100.0, engine_spawned=(lbl == 'B'))
             for lbl, mu in (('A', (1, 2, 3)), ('B', (4, 5, 6)))]
    rs = HybridPolymerSystem(T=700.0, P=1.0e5, initial_mole_fractions={core[0]: 1.0}, V_poly=1.0,
                             polymer_pools=pools, mass_transfer=[],
                             gas_species_mask=np.array([True] + [False] * 6), constant_gas_volume=False,
                             initial_polymer_moments={'A': (0.01, 0.5, 30.0)}, termination=[])
    with caplog.at_level(logging.WARNING):
        rs.initialize_model(core, [], [], [])
        rs.residual(0.0, rs.y, np.zeros_like(rs.y))
    return rs


def _exhausted_lines(caplog, label):
    return [rec.getMessage() for rec in caplog.records
            if rec.getMessage().startswith(f'POOL EXHAUSTION CENSUS: pool {label} exhausted')]


def test_a_pool_born_empty_is_not_announced_exhausted_at_t0(caplog):
    _two_pool_exhaustion(caplog)
    assert _exhausted_lines(caplog, 'B') == []


def test_a_pool_that_held_chains_and_drained_is_still_announced(caplog):
    rs = _two_pool_exhaustion(caplog)
    y = np.array(rs.y, dtype=float)
    y[1:4] = 0.0
    with caplog.at_level(logging.WARNING):
        rs.residual(1.0, y, np.zeros_like(y))
    assert len(_exhausted_lines(caplog, 'A')) == 1


class TestPromotionOnALiveRow:
    """A spawned daughter whose only row is superseded or refused has nothing
    to integrate: it is not promoted. It is promoted at the enlarge boundary
    after a live row reaches it."""

    def _deck(self, monkeypatch):
        parent = _parent(QSSA)
        styrene = _spc(MONOMER_CLOSED, 'styrene')
        parent.monomer_product_species = styrene
        model = _model_with_parent(parent, _spc('N#N', 'N2'), styrene)
        return model, parent

    def test_a_daughter_whose_only_row_is_superseded_is_not_promoted(self, monkeypatch):
        model, parent = self._deck(monkeypatch)
        rxn, tail = _generate(model, parent, HOMOLYSIS_GAS, concerted=False)
        _enlarge_boundary(model, monkeypatch)
        assert tail not in model.core.species

    def test_a_daughter_whose_only_row_is_refused_is_not_promoted(self, monkeypatch):
        model, parent = self._deck(monkeypatch)
        rxn, tail = _generate(model, parent, CONCERTED_GAS, concerted=True)
        rxn.polymer_refused = True
        _enlarge_boundary(model, monkeypatch)
        assert tail not in model.core.species

    def test_it_is_promoted_when_its_first_live_row_appears(self, monkeypatch):
        model, parent = self._deck(monkeypatch)
        rxn, tail = _generate(model, parent, HOMOLYSIS_GAS, concerted=False)
        _enlarge_boundary(model, monkeypatch)
        assert tail not in model.core.species
        styrene = next(s for s in model.core.species if s.label == 'styrene')
        live = Reaction(reactants=[parent], products=[tail, styrene], kinetics=KIN, reversible=False)
        live.polymer_flux_archetype = sp_mod.FLUX_VOLATILE_EJECTION
        live.polymer_eject_units = 1.0
        _stamp(live, 'concerted_molecular')
        model.add_reaction_to_edge(live)
        _enlarge_boundary(model, monkeypatch)
        assert tail in model.core.species

    def test_a_daughter_with_a_live_row_is_promoted(self, monkeypatch):
        model, parent = self._deck(monkeypatch)
        rxn, tail = _generate(model, parent, CONCERTED_GAS, concerted=True)
        _enlarge_boundary(model, monkeypatch)
        assert tail in model.core.species
