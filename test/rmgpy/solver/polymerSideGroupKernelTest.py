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
The side-group kernel with a polyatomic leaving group, in the solver.

A methyl methacrylate-like carrier (test data only; repeat unit C5H8O2)
loses its ester group as the C2H3O2 radical:

- a stationary carrier keeps the removable-group inventory Z; integrated
  through full side-group loss (Z -> 0), the element inventory (code
  convention, the lost groups booked off the carrier) and the mass (condensed
  plus gas) stay what they were at t = 0;
- a carrier whose chains another kernel moves carries the loss as its X-loss
  feature daughter: parent dmu_j -= k*s*mu_{j+1}, the daughter gains exactly
  that, the gas gains k*s*mu1, and there is no inventory slot. Integrated
  until the carrier is gone, elements (the daughter's per-chain defect atoms
  off its chains) and mass close.

A setup step that fails raises RuntimeError, never an AssertionError.
"""
import contextlib
import io

import numpy as np
import pytest

from rmgpy.molecule import Molecule
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
from rmgpy.solver.base import TerminationTime
from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig
from rmgpy.species import Species

T_K = 800.0
DP = 50.0
MOMENTS = (1.0e-3, DP * 1.0e-3, 2.0 * DP * DP * 1.0e-3)
REPEAT = (("C", 5), ("H", 8), ("O", 2))
ENDS = (("C", 1), ("H", 3), ("H", 1))
UNIT_MW = Molecule().from_smiles("C=C(C)C(=O)OC").get_molecular_weight() * 1000.0
GAS = "CO[C]=O"
GAS_ATOMS = {"C": 2, "H": 3, "O": 2}
DAUGHTER = "poly_sidegrp_ester"


def _setup(ok, what):
    if not ok:
        raise RuntimeError("setup failed: " + what)


def _spc(smiles, label):
    s = Species(molecule=[Molecule(smiles=smiles)])
    s.label = label
    return s


def _rate(k):
    """An Arrhenius triplet with k(T) == k exactly."""
    return dict(A=float(k), n=0.0, Ea=0.0)


def _gas_mw():
    return Molecule().from_smiles(GAS).get_molecular_weight() * 1000.0


def _channel(k):
    return dict(label="ester", site_selector="aliphatic", sites_per_unit=1.0, gas_product=GAS, **_rate(k))


def _pool(label, mu_indices, **kw):
    cfg = dict(xs=2, explicit_dp_to_species_index={}, k_scission=0.0, k_unzip=0.0, tail_kinetics=None,
               monomer_mw_g_mol=UNIT_MW, repeat_unit_atoms=REPEAT, end_group_atoms=ENDS)
    cfg.update(kw)
    return PolymerPoolConfig(label=label, mu_indices=mu_indices, **cfg)


def stationary(k_sgh):
    """The carrier alone: inventory accounting."""
    core = [_spc("N#N", "N2"), _spc(GAS, "ester_rad"), _spc("CCCC", "poly"),
            _spc("CO", "poly_mu0"), _spc("C=O", "poly_mu1"), _spc("C#N", "poly_mu2")]
    pools = [_pool("poly", (3, 4, 5), monomer_poly_index=None,
                   side_group_homolysis=[_channel(k_sgh)], side_group_gas_indices=(1,))]
    return core, pools, {"poly": MOMENTS}


def moved(k_sgh, k_dep=None, sgh=True):
    """The carrier unzips from its chain ends (end-radical depropagation)
    beside the side-group kernel: the X-loss feature daughter."""
    core = [_spc("N#N", "N2"), _spc("C=C(C)C(=O)OC", "monomer"), _spc(GAS, "ester_rad"),
            _spc("CCCC", "poly"), _spc("CO", "poly_mu0"), _spc("C=O", "poly_mu1"), _spc("C#N", "poly_mu2"),
            _spc("CCC", DAUGHTER), _spc("CCO", DAUGHTER + "_mu0"), _spc("CC=O", DAUGHTER + "_mu1"),
            _spc("CC#N", DAUGHTER + "_mu2")]
    kw = dict(side_group_homolysis=[_channel(k_sgh)], side_group_gas_indices=(2,)) if sgh else {}
    if k_dep is not None:
        kw["k_depropagation"] = _rate(k_dep)
    else:
        kw["k_unzip"] = 1.0e-3
    pools = [_pool("poly", (4, 5, 6), monomer_poly_index=1, **kw),
             _pool(DAUGHTER, (8, 9, 10), monomer_poly_index=None, chain_mass_defect_g_mol=_gas_mw())]
    return core, pools, {"poly": MOMENTS, DAUGHTER: (0.0, 0.0, 0.0)}


def _system(core, pools, moments):
    mask = np.array([not any(s.label == p.label or s.label.startswith(p.label + "_mu") for p in pools)
                     for s in core], dtype=bool)
    rs = HybridPolymerSystem(T=T_K, P=1.0e5, initial_mole_fractions={core[0]: 1.0}, V_poly=1.0,
                             polymer_pools=pools, mass_transfer=[], gas_species_mask=mask,
                             constant_gas_volume=False, initial_polymer_moments=moments, termination=[])
    with contextlib.redirect_stdout(io.StringIO()):
        rs.initialize_model(core, [], [], [])
    return rs


def _attempt_system(core, pools, moments):
    """The built system, or the ValueError the build raised (a refusal is an
    outcome to assert on)."""
    try:
        return _system(core, pools, moments)
    except ValueError as e:
        return e


def _simulate(rs, core, t_end):
    rs.termination.append(TerminationTime((t_end, "s")))
    ms = ModelSettings(tol_keep_in_edge=0.0, tol_move_to_core=1.0e8, tol_interrupt_simulation=1.0e8)
    with contextlib.redirect_stdout(io.StringIO()):
        rs.simulate(list(core), [], [], [], [], [], model_settings=ms, simulator_settings=SimulatorSettings())


def _mass(rs, core, y):
    """Condensed mass plus every gas species' mass plus the exhausted
    residue [g]."""
    gas = sum(float(y[i]) * core[i].molecule[0].get_molecular_weight() * 1000.0
              for i in range(len(core)) if rs.gas_species_mask[i])
    residue = sum(rs.exhausted_residue_mass_g().values()) if y is rs.y else 0.0
    return rs.get_total_polymer_condensed_mass_g(y) + gas + residue


def _inventory(rs, y, defect_atoms):
    out = dict(rs.code_convention_inventory(y, defect_atoms=defect_atoms))
    if y is rs.y:
        for el, n in rs.exhausted_residue_atoms(defect_atoms=defect_atoms).items():
            out[el] = out.get(el, 0.0) + n
    return out


class TestInventoryAccounting:

    def test_a_stationary_carrier_loses_its_polyatomic_groups_and_the_ledger_closes(self):
        """Through full side-group loss: Z goes to zero, the gas holds one
        ester radical per repeat unit, and elements and mass close."""
        core, pools, moments = stationary(10.0)
        rs = _attempt_system(core, pools, moments)
        assert not isinstance(rs, Exception), f"refused: {rs}"
        inv0, m0 = _inventory(rs, rs.y0, None), _mass(rs, core, rs.y0)
        _simulate(rs, core, 3.0)
        z = float(rs.y[int(rs.sgh_z_slot[0])])
        assert z <= 1.0e-9 * MOMENTS[1], z
        assert float(rs.y[1]) == pytest.approx(MOMENTS[1], rel=1e-6)
        now = _inventory(rs, rs.y, None)
        for el in ("C", "H", "O"):
            assert now.get(el, 0.0) == pytest.approx(inv0[el], rel=1e-7), (el, now, inv0)
        assert _mass(rs, core, rs.y) == pytest.approx(m0, rel=1e-7)


class TestDaughterAccounting:

    def test_the_x_loss_daughter_takes_exactly_the_carriers_debit(self):
        """The kernel's share of the residual (with the channel less without
        it): parent dmu0 = -k*s*mu1, dmu1 = -k*s*mu2, dmu2 < 0; the daughter
        gains exactly that; the gas gains k*s*mu1; no inventory slot."""
        k = 10.0
        rs = _attempt_system(*moved(k))
        assert not isinstance(rs, Exception), f"refused: {rs}"
        base = _system(*moved(k, sgh=False))
        assert rs.num_sgh_z == 0
        d = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0] - base.residual(0.0, base.y, np.zeros_like(base.y))[0]
        mu0, mu1, mu2 = MOMENTS
        assert d[4] == pytest.approx(-k * mu1, rel=1e-12)
        assert d[5] == pytest.approx(-k * mu2, rel=1e-12)
        assert d[6] < 0.0
        assert list(d[8:11]) == pytest.approx([-x for x in d[4:7]], rel=1e-9)
        assert d[2] == pytest.approx(k * mu1, rel=1e-12)

    def test_a_moved_carrier_loses_its_groups_into_the_daughter_and_the_ledger_closes(self):
        """Until the carrier is gone (side-group loss and its own unzip):
        elements, with the daughter's defect atoms off each of its chains,
        and mass close; the daughter holds one chain per lost group."""
        core, pools, moments = moved(10.0, k_dep=1.0)
        rs = _attempt_system(core, pools, moments)
        assert not isinstance(rs, Exception), f"refused: {rs}"
        defects = {DAUGHTER: GAS_ATOMS}
        inv0, m0 = _inventory(rs, rs.y0, defects), _mass(rs, core, rs.y0)
        _simulate(rs, core, 5.0)
        daughter_mu0 = float(rs.y[8])
        assert daughter_mu0 == pytest.approx(float(rs.y[2]), rel=1e-6)
        assert daughter_mu0 > 0.5 * MOMENTS[0]
        now = _inventory(rs, rs.y, defects)
        for el in ("C", "H", "O"):
            assert now.get(el, 0.0) == pytest.approx(inv0[el], rel=1e-7), (el, now, inv0)
        assert _mass(rs, core, rs.y) == pytest.approx(m0, rel=1e-7)
