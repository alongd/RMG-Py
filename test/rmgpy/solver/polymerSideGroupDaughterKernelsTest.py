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
A side-group X-loss daughter carries its parent's backbone kernels, and every
side group a chain loses is counted.

A chain that lost a side group keeps the parent's repeat units, ends and
backbone (the group holds no stitch atom), so the parent's backbone kernels
apply to it: its unzip releases a monomer, or the modified monomer (the unit
less its group) when the unzipping unit lost its group, and its remaining
units still lose their groups.

Producer: the daughter inherits the parent's unzip and scission kernels and its
own channel, and names its modified monomer, for an ester, a phenyl, a chlorine
and a methyl leaving group (test-data parents).

Solver: a methyl methacrylate-like carrier with the radical QSSA unzip and the
ester channel (the constants of a kernels-only test deck; test data only), its
daughter configured through derive_daughter_pool_configs as a run configures
it. Once each chain has lost one group, the monomer yield (monomer plus
modified monomer) is the parent-only carrier's; the gas holds every lost group,
s*mu1(0) at full loss and s*mu1(0)*(1 - exp(-k t)) on the way; elements and
mass close.

A setup step that fails raises RuntimeError, never an AssertionError.
"""
import contextlib
import copy
import io
import math

import numpy as np
import pytest

from rmgpy.molecule import Molecule
from rmgpy.polymer import Polymer
from rmgpy.rmg.polymer_input import derive_daughter_pool_configs
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
from rmgpy.solver.base import TerminationTime
from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig
from rmgpy.species import Species

# (parent, repeat unit, channel, leaving group, selector, modified monomer)
LEAVING_GROUPS = [
    ("PMMA", "[CH2][C](C)C(=O)OC", "ester", "CO[C]=O", "aliphatic", "C=[C]C"),
    ("PS", "[CH2][CH]c1ccccc1", "phenyl", "[c]1ccccc1", "aliphatic", "C=[CH]"),
    ("PVC", "[CH2][CH]Cl", "chlorine", "[Cl]", "aliphatic", "C=[CH]"),
    ("PP", "[CH2][CH]C", "methyl", "[CH3]", "aliphatic", "C=[CH]"),
]

# The kernels-only test deck's carrier (test data): 700 K, 0.01 mol of chains
# of Mn 5000 g/mol in 50 g of melt at 1180 kg/m^3.
T_K = 700.0
V_POLY = 0.05 / 1180.0
MOMENTS = (0.01, 0.49942182628652415, 29.93065926856405)
QSSA = dict(initiation=dict(A=6.082608239264e23, n=-2.45832700751, Ea=320216.873028),
            depropagation=dict(A=3.173648881804e18, n=-1.69158135443, Ea=114750.084903),
            termination=dict(A=2.655663993173e8, n=-0.553042197084, Ea=728.80426424),
            transfer=dict(A=1.057374272215e-10, n=5.72127747789, Ea=26918.0295221),
            efficiency=1.0, monomer_yield=1.0, basis="backbone_bonds_mu1_minus_mu0")
MONOMER = "[CH2][C](C)C(=O)OC"
MONOMER_PRODUCT = "C=C(C)C(=O)OC"
GAS = "CO[C]=O"
GAS_ATOMS = {"C": 2, "H": 3, "O": 2}
MODIFIED = "C=[C]C"
REPEAT = (("C", 5), ("H", 8), ("O", 2))
ENDS = (("C", 1), ("H", 3), ("H", 1))
DAUGHTER = "PMMA_sidegrp_ester"


def _setup(ok, what):
    if not ok:
        raise RuntimeError("setup failed: " + what)


def _spc(smiles, label):
    s = Species(molecule=[Molecule(smiles=smiles)])
    s.label = label
    return s


def _channel(label, gas, selector, k):
    return dict(label=label, A=float(k), n=0.0, Ea=0.0, site_selector=selector, sites_per_unit=1.0,
                gas_product=gas)


def _parent(label, unit, channels, **kw):
    return Polymer(label=label, monomer=unit, end_groups=["[CH3]", "[H]"], cutoff=3, Mn=5000.0, Mw=6000.0,
                   initial_mass=1.0, side_group_homolysis=channels, **kw)


def _unit_mw():
    return Molecule().from_smiles(MONOMER_PRODUCT).get_molecular_weight() * 1000.0


class TestProducer:

    @pytest.mark.parametrize("parent,unit,ch,gas,sel,modified", LEAVING_GROUPS)
    def test_the_daughter_inherits_the_unzip_and_its_own_channel(self, parent, unit, ch, gas, sel, modified):
        """The daughter of a radical-QSSA carrier carries a copy of the
        parent's unzip channel and the channel it was spawned by, and names
        the modified monomer its defect units release."""
        channel = _channel(ch, gas, sel, 1.0e3)
        p = _parent(parent, unit, [channel], radical_qssa_unzip=copy.deepcopy(QSSA))
        daughters = p.generate_side_loss_daughters()
        _setup(len(daughters) == 1, f"one daughter for {parent}")
        d = daughters[0]
        assert d.radical_qssa_unzip == p.radical_qssa_unzip
        assert d.radical_qssa_unzip is not p.radical_qssa_unzip
        assert d.side_group_homolysis == [channel]
        mod = getattr(d, "side_loss_monomer", None)
        assert mod is not None and mod.is_isomorphic(Molecule().from_smiles(modified)), mod

    def test_the_daughter_inherits_the_legacy_unzip_and_scission(self):
        """k_unzip and k_scission are backbone kernels too."""
        p = _parent("PMMA", MONOMER, [_channel("ester", GAS, "aliphatic", 1.0e3)], k_unzip=2.0e-3,
                    k_scission=5.0e-4)
        d = p.generate_side_loss_daughters()[0]
        assert (d.k_unzip, d.k_scission) == (2.0e-3, 5.0e-4)


def _systems(k_sgh, qssa=QSSA, sgh=True):
    """The carrier (and, with the channel, its X-loss daughter configured
    from the daughter Polymer through derive_daughter_pool_configs, with the
    gas and modified-monomer species registered as the deck setup does)."""
    channel = _channel("ester", GAS, "aliphatic", k_sgh)
    n2, monomer, gas, modified = (_spc("N#N", "N2"), _spc(MONOMER_PRODUCT, "monomer"), _spc(GAS, "ester_rad"),
                                  _spc(MODIFIED, "modified_monomer"))
    core = [n2, monomer, gas, modified, _spc("CCCC", "PMMA"), _spc("CO", "PMMA_mu0"), _spc("C=O", "PMMA_mu1"),
            _spc("C#N", "PMMA_mu2")]
    kw = dict(side_group_homolysis=[channel], side_group_gas_indices=(2,)) if sgh else {}
    parent_cfg = PolymerPoolConfig(label="PMMA", mu_indices=(5, 6, 7), xs=2, explicit_dp_to_species_index={},
                                   k_scission=0.0, k_unzip=0.0, tail_kinetics=None, monomer_poly_index=1,
                                   monomer_mw_g_mol=_unit_mw(), repeat_unit_atoms=REPEAT, end_group_atoms=ENDS,
                                   radical_qssa_unzip=copy.deepcopy(qssa), **kw)
    pools, moments = [parent_cfg], {"PMMA": MOMENTS}
    if sgh:
        parent = _parent("PMMA", MONOMER, [channel], radical_qssa_unzip=copy.deepcopy(qssa))
        parent.monomer_product_species = monomer
        parent.side_group_gas_species = [gas]
        daughters = parent.generate_side_loss_daughters()
        _setup(len(daughters) == 1, "one daughter")
        d = daughters[0]
        if getattr(d, "side_loss_monomer", None) is not None:
            d.side_loss_monomer_species = modified
        core += [d, _spc("CCO", DAUGHTER + "_mu0"), _spc("CC=O", DAUGHTER + "_mu1"), _spc("CC#N", DAUGHTER + "_mu2")]
        spc_map = {s: i for i, s in enumerate(core)}
        derived = derive_daughter_pool_configs(core, spc_map, ["PMMA"])
        _setup([c.label for c in derived] == [DAUGHTER], f"derived {[c.label for c in derived]}")
        pools.append(derived[0])
        moments[DAUGHTER] = (0.0, 0.0, 0.0)
    mask = np.array([not any(s.label == p.label or s.label.startswith(p.label + "_mu") for p in pools)
                     for s in core], dtype=bool)
    rs = HybridPolymerSystem(T=T_K, P=1.0e5, initial_mole_fractions={n2: 1.0}, V_poly=V_POLY,
                             polymer_pools=pools, mass_transfer=[], gas_species_mask=mask,
                             constant_gas_volume=False, initial_polymer_moments=moments, termination=[])
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            rs.initialize_model(core, [], [], [])
    except Exception as e:
        raise RuntimeError(f"setup failed: {e!r}")
    return rs, core


def _simulate(rs, core, t_end):
    rs.termination.append(TerminationTime((t_end, "s")))
    ms = ModelSettings(tol_keep_in_edge=0.0, tol_move_to_core=1.0e8, tol_interrupt_simulation=1.0e8)
    with contextlib.redirect_stdout(io.StringIO()):
        rs.simulate(list(core), [], [], [], [], [], model_settings=ms,
                    simulator_settings=SimulatorSettings(atol=1e-16, rtol=1e-8))


def _index(core, label):
    return next(i for i, s in enumerate(core) if s.label == label)


def _released(rs, core):
    """Monomer plus modified monomer released [mol]."""
    return float(rs.y[_index(core, "monomer")]) + float(rs.y[_index(core, "modified_monomer")])


def _mass(rs, core, y):
    gas = sum(float(y[i]) * core[i].molecule[0].get_molecular_weight() * 1000.0
              for i in range(len(core)) if rs.gas_species_mask[i])
    residue = sum(rs.exhausted_residue_mass_g().values()) if y is rs.y else 0.0
    return rs.get_total_polymer_condensed_mass_g(y) + gas + residue


def _inventory(rs, y):
    defects = {DAUGHTER: GAS_ATOMS}
    out = dict(rs.code_convention_inventory(y, defect_atoms=defects))
    if y is rs.y:
        for el, n in rs.exhausted_residue_atoms(defect_atoms=defects).items():
            out[el] = out.get(el, 0.0) + n
    return out


T_END = 2.0


def _release_rate(rs, core):
    """Released units per second to the end of the integration [mol/s]: the
    integrator stops past T_END at a time of its own, and with B nearly
    constant the release is steady."""
    return _released(rs, core) / float(rs.t)


def _parent_only_rate(qssa):
    rs, core = _systems(0.0, qssa=qssa, sgh=False)
    _simulate(rs, core, T_END)
    return float(rs.y[_index(core, "monomer")]) / float(rs.t)


# Termination-limited radicals (no transfer): R ~ sqrt(B), where one radical
# population and one per pool differ (by up to sqrt(2) at an even split).
NO_TRANSFER = dict(QSSA, transfer=None)


@pytest.fixture(scope="module")
def parent_only():
    return _parent_only_rate(QSSA)


class TestDaughterKernels:

    def test_the_monomer_yield_survives_the_loss_of_one_group_per_chain(self, parent_only):
        """Every chain loses its first group within milliseconds; the units
        it releases per second are the parent-only carrier's."""
        _setup(parent_only > 0.0, "parent-only release")
        rs, core = _systems(360.0)
        _simulate(rs, core, T_END)
        assert _release_rate(rs, core) == pytest.approx(parent_only, rel=1.0e-2), (_release_rate(rs, core), parent_only)

    def test_the_gas_holds_every_lost_group(self):
        """At k*t >> 1 every unit has lost its ester: s*mu1(0) of gas, not
        one per chain (mu0)."""
        rs, core = _systems(360.0)
        _simulate(rs, core, T_END)
        lost = float(rs.y[_index(core, "ester_rad")])
        # the units released before their group left never lose it: at most
        # the released amount short of s*mu1(0)
        assert MOMENTS[1] - _released(rs, core) - 1.0e-9 <= lost <= MOMENTS[1] * (1.0 + 1.0e-9), lost

    def test_every_site_is_lost_at_its_own_rate(self):
        """Without unzip, each site is lost independently at k: the gas is
        s*mu1(0)*(1 - exp(-k t)) exactly, whichever pool holds the chain."""
        frozen = copy.deepcopy(QSSA)
        frozen["initiation"] = dict(A=1.0e-30, n=0.0, Ea=0.0)
        k = 0.4
        rs, core = _systems(k, qssa=frozen)
        _simulate(rs, core, T_END)
        lost = float(rs.y[_index(core, "ester_rad")])
        assert lost == pytest.approx(MOMENTS[1] * -math.expm1(-k * float(rs.t)), rel=1.0e-6), (lost, rs.t)

    def test_a_partial_loss_releases_what_the_parent_alone_would(self):
        """With the chains split between parent and daughter through the
        run (k*s*DP*t ~ 1) and termination-limited radicals, the two are one
        melt with one radical population: together they release what the
        parent alone would, not up to sqrt(2) more."""
        rs, core = _systems(0.01, qssa=NO_TRANSFER)
        _simulate(rs, core, T_END)
        moved = float(rs.y[rs.polymer_pools[1].mu_indices[0]]) / MOMENTS[0]
        _setup(0.2 < moved < 0.8, f"a chain-level split ({moved:.2f} of the chains moved)")
        alone = _parent_only_rate(NO_TRANSFER)
        assert _release_rate(rs, core) == pytest.approx(alone, rel=1.0e-3), (_release_rate(rs, core), alone)

    def test_elements_and_mass_close(self):
        """Through the loss and the unzip, the lost groups (the first per
        chain on the defect, the rest on the inventory) and the released
        units close the element and mass ledger."""
        rs, core = _systems(0.4)
        inv0, m0 = _inventory(rs, rs.y0), _mass(rs, core, rs.y0)
        _simulate(rs, core, T_END)
        now = _inventory(rs, rs.y)
        for el in ("C", "H", "O"):
            assert now.get(el, 0.0) == pytest.approx(inv0[el], rel=1.0e-9), (el, now, inv0)
        assert _mass(rs, core, rs.y) == pytest.approx(m0, rel=1.0e-9)
