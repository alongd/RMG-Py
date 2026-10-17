"""
Declaration -> composition path of the chain end-group bookkeeping: the
producers (PolymerPool.to_config, derive_daughter_pool_configs, the artifact
reader) must turn declared end groups and the monomer into the element tuples
the ledger uses. The closure tests build configs by hand; these go through
the producers.
"""
import dataclasses
import os
import sys

import numpy as np
import pytest

from rmgpy.molecule import Molecule
from rmgpy.polymer import Polymer
from rmgpy.rmg.polymer_input import (PolymerPool, _end_group_atoms,
                                     derive_daughter_pool_configs)
from rmgpy.solver.polymer import HybridPolymerSystem
from rmgpy.species import Species

sys.path.insert(0, os.path.dirname(__file__))
from solverPolymerEndGroupLedgerTest import _hook, _spc  # noqa: E402

PS_MONOMER = "[CH2][CH]c1ccccc1"


def _polymer(end_groups, label="X", monomer="[CH2][CH2]"):
    return Polymer(label=label, monomer=monomer, end_groups=end_groups, cutoff=3,
                   Mn=5000.0, Mw=6000.0, initial_mass=0.001)


@pytest.mark.parametrize("groups, expected", [
    (["[H]", "[H]"], (("H", 2),)),
    (["[CH3]", "[H]"], (("C", 1), ("H", 4))),
    (["[CH2]c1ccccc1", "[H]"], (("C", 7), ("H", 8))),
])
def test_end_group_atoms_from_declared_polymer(groups, expected):
    assert _end_group_atoms(_polymer(groups).end_groups) == expected


def test_end_group_atoms_none_when_nothing_declared():
    assert _end_group_atoms(None) is None
    assert _end_group_atoms([]) is None


def _pool_config(polymer, monomer_smiles="[CH2][CH2]"):
    proxy = polymer   # a Polymer is the pool's proxy Species and carries end_groups
    mus = [_spc("CO", "X_mu0"), _spc("C=O", "X_mu1"), _spc("C#N", "X_mu2")]
    pool = PolymerPool(label="X", xs=3, monomer=Molecule().from_smiles(monomer_smiles),
                       explicit_map={}, mu_species=mus, proxy_species=proxy)
    core = [proxy] + mus
    return pool.to_config({s: i for i, s in enumerate(core)})


def test_to_config_carries_declared_end_groups_and_repeat_unit():
    cfg = _pool_config(_polymer(["[CH3]", "[H]"]))
    assert cfg.end_group_atoms == (("C", 1), ("H", 4))
    assert cfg.repeat_unit_atoms == (("C", 2), ("H", 4))
    cfg = _pool_config(_polymer(["[H]", "[H]"], monomer=PS_MONOMER), PS_MONOMER)
    assert cfg.end_group_atoms == (("H", 2),)
    assert cfg.repeat_unit_atoms == (("C", 8), ("H", 8))


def test_derived_daughter_inherits_end_groups_and_repeat_unit():
    parent = _spc("CCCC", "PS")
    pm = [_spc("CO", "PS_mu0"), _spc("C=O", "PS_mu1"), _spc("C#N", "PS_mu2")]
    daughter = _polymer(["[CH3]", "[H]"], label="PS_d1", monomer=PS_MONOMER)
    dm = []
    for lab in ("PS_d1_mu0", "PS_d1_mu1", "PS_d1_mu2"):
        s = _spc("[Ne]", lab)
        s.reactive = False
        s.is_moment_dummy = True
        dm.append(s)
    core = [parent] + pm + [daughter] + dm
    derived = derive_daughter_pool_configs(
        core, {s: i for i, s in enumerate(core)}, existing_pool_labels={"PS"})
    assert [c.label for c in derived] == ["PS_d1"]
    assert derived[0].end_group_atoms == (("C", 1), ("H", 4))
    assert derived[0].repeat_unit_atoms == (("C", 8), ("H", 8))


def test_artifact_reader_parses_end_groups_and_monomer():
    from rmgpy.tools.polymer_moments_runner import (
        _end_group_atoms_from_pool_entry, _monomer_atoms_from_pool_entry)
    entry = {
        "end_groups": ["[CH3]", "[H]"],
        "monomer_adj_list": ("multiplicity 3\n1 *1 C u1 p0 c0 {2,S} {3,S} {4,S}\n"
                             "2 *2 C u1 p0 c0 {1,S} {5,S} {6,S}\n3    H u0 p0 c0 {1,S}\n"
                             "4    H u0 p0 c0 {1,S}\n5    H u0 p0 c0 {2,S}\n"
                             "6    H u0 p0 c0 {2,S}\n"),
    }
    assert _end_group_atoms_from_pool_entry(entry) == (("C", 1), ("H", 4))
    assert _monomer_atoms_from_pool_entry(entry) == (("C", 2), ("H", 4))
    assert _end_group_atoms_from_pool_entry({"end_groups": ["[H]", "[H]"]}) == (("H", 2),)
    assert _end_group_atoms_from_pool_entry({}) is None
    assert _monomer_atoms_from_pool_entry({"monomer_smiles": "[CH2][CH2]"}) == (("C", 2), ("H", 4))


def test_declared_pool_drives_physical_closure_end_to_end():
    """A pool built by the producer from a declaration (CH3/H ends, so the
    ledger must count C as well as H) closes scission + recombination + gas."""
    cfg = _pool_config(_polymer(["[CH3]", "[H]"]))
    mon = _spc("C=C", "mon")
    pool_sp = _spc("CCCC", "X")  # core placeholder at index 0 (config indices fixed above)
    mus = [_spc("CO", "X_mu0"), _spc("C=O", "X_mu1"), _spc("C#N", "X_mu2")]
    core = [pool_sp] + mus + [mon]
    cfg = dataclasses.replace(cfg, monomer_poly_index=4,
                              tail_kinetics=_hook(True, True, True, 4))
    rs = HybridPolymerSystem(
        T=800.0, P=1.0e5, initial_mole_fractions={mon: 0.5}, V_poly=1.0,
        polymer_pools=[cfg], mass_transfer=[],
        gas_species_mask=np.array([False] * 4 + [True]), constant_gas_volume=False,
        initial_polymer_moments={"X": (1.0, 5.0, 30.0)}, termination=[])
    rs.initialize_model(core, [], [], [])
    inv0 = rs.physical_elemental_inventory(rs.y0)
    rs.advance(3.0)
    inv = rs.physical_elemental_inventory(rs.y)
    for el in ("C", "H"):
        assert abs(inv[el] - inv0[el]) <= 1e-10 * inv0[el], (el, inv, inv0)
    orph = rs.orphaned_end_atoms(rs.y)
    dm0 = float(rs.y0[1] - rs.y[1])
    assert orph["C"] == pytest.approx(1.0 * dm0, rel=1e-6, abs=1e-12)
    assert orph["H"] == pytest.approx(4.0 * dm0, rel=1e-6, abs=1e-12)
    assert abs(dm0) > 1e-3
