"""
Chain end-group bookkeeping of the method-of-moments polymer solver.

The moment closure books condensed atoms as mu1 * (repeat unit); the two end
groups of every chain are not in that count. A chain that vanishes (recombined
away, or fully converted to gas) takes its end atoms with it, and a scission
creates two ends out of nothing, so the PHYSICAL elemental balance (units plus
mu0 * end groups plus gas) drifts although the code convention closes. The
solver therefore carries one non-feeding "orphaned end groups" accumulator per
pool whose declared ``end_group_atoms`` make the physical balance close.
"""
import numpy as np
import pytest

from rmgpy.molecule import Molecule
from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig
from rmgpy.species import Species

S_SCISSION, S_RECOMB, S_GAS = 0.3, 0.2, 0.1


def _spc(smiles, label):
    s = Species(molecule=[Molecule(smiles=smiles)])
    s.label = label
    return s


def _hook(scission, recomb, gas, mon_idx):
    def f(T, P, mu0, mu1, mu2, mu3):
        d0 = 0.0
        d1 = 0.0
        d2 = 0.0
        src = {}
        if scission:
            d0 += S_SCISSION * (mu1 - mu0)
        if recomb:
            d0 -= S_RECOMB * mu0
        if gas:
            # whole chains leave as gas monomer units: mu0 and mu1 shrink together
            d0 -= S_GAS * mu0
            d1 -= S_GAS * mu1
            d2 -= S_GAS * mu2
            src[mon_idx] = S_GAS * mu1
        return d0, d1, d2, src
    return f


def _run(scission, recomb, gas, tend=3.0):
    pool_sp = _spc("CCCCCC", "A")
    mu = [_spc("CO", "A_mu0"), _spc("C=O", "A_mu1"), _spc("C#N", "A_mu2")]
    mon = _spc("C=C", "mon")
    core = [pool_sp] + mu + [mon]
    mask = np.array([False] * 4 + [True], dtype=bool)
    pool = PolymerPoolConfig(
        label="A", xs=2, explicit_dp_to_species_index={}, mu_indices=(1, 2, 3),
        monomer_poly_index=4, k_scission=0.0, k_unzip=0.0,
        tail_kinetics=_hook(scission, recomb, gas, 4),
        repeat_unit_atoms=(("C", 2), ("H", 4)),
        end_group_atoms=(("H", 1), ("H", 1)),
    )
    rs = HybridPolymerSystem(
        T=800.0, P=1.0e5, initial_mole_fractions={mon: 0.5}, V_poly=1.0,
        polymer_pools=[pool], mass_transfer=[], gas_species_mask=mask.copy(),
        constant_gas_volume=False,
        initial_polymer_moments={"A": (1.0, 5.0, 30.0)}, termination=[])
    rs.initialize_model(core, [], [], [])
    inv0 = rs.physical_elemental_inventory(rs.y0)
    mu0_0 = float(rs.y0[1])
    rs.advance(tend)
    return rs, inv0, mu0_0


CASES = {"scission": (True, False, False), "recombination": (False, True, False),
         "gas": (False, False, True), "all": (True, True, True)}


@pytest.mark.parametrize("name", list(CASES))
def test_physical_elements_close_with_end_group_ledger(name):
    rs, inv0, mu0_0 = _run(*CASES[name])
    inv = rs.physical_elemental_inventory(rs.y)
    for el in ("C", "H"):
        assert abs(inv[el] - inv0[el]) <= 1e-10 * inv0[el], (name, el, inv, inv0)
    # the ledger is the named output; H it holds is 2 per orphaned chain
    mu0 = float(rs.y[1])
    orphans = rs.orphaned_end_atoms(rs.y)
    assert orphans.get("H", 0.0) == pytest.approx(2.0 * (mu0_0 - mu0), rel=1e-6, abs=1e-12)
    assert abs(orphans.get("H", 0.0)) > 1e-3, "synthetic run must move mu0"


def test_code_convention_still_closes_and_excludes_ledger():
    rs, inv0, _ = _run(*CASES["all"])
    code = rs.code_convention_inventory(rs.y)
    code0 = rs.code_convention_inventory(rs.y0)
    assert abs(code["H"] - code0["H"]) <= 1e-10 * code0["H"]
    assert abs(code["C"] - code0["C"]) <= 1e-10 * code0["C"]


def test_pool_without_declared_end_groups_is_refused_by_the_reporters():
    pool_sp = _spc("CCCCCC", "A")
    mu = [_spc("CO", "A_mu0"), _spc("C=O", "A_mu1"), _spc("C#N", "A_mu2")]
    pool_sp2 = _spc("CCCCC", "B")
    mu2 = [_spc("CCO", "B_mu0"), _spc("CC=O", "B_mu1"), _spc("CC#N", "B_mu2")]
    core = [pool_sp] + mu + [pool_sp2] + mu2
    mask = np.array([False] * 8, dtype=bool)
    pools = [PolymerPoolConfig(label="A", xs=2, explicit_dp_to_species_index={},
                               mu_indices=(1, 2, 3), end_group_atoms=(("H", 1), ("H", 1)),
                               repeat_unit_atoms=(("C", 2), ("H", 4))),
             PolymerPoolConfig(label="B", xs=2, explicit_dp_to_species_index={},
                               mu_indices=(5, 6, 7))]
    rs = HybridPolymerSystem(
        T=800.0, P=1.0e5, initial_mole_fractions={}, V_poly=1.0,
        polymer_pools=pools, mass_transfer=[], gas_species_mask=mask.copy(),
        constant_gas_volume=False,
        initial_polymer_moments={"A": (1.0, 5.0, 30.0), "B": (1.0, 5.0, 30.0)},
        termination=[])
    rs.initialize_model(core, [], [], [])
    with pytest.raises(ValueError, match="end_group_atoms"):
        rs.orphaned_end_atoms(rs.y)
    with pytest.raises(ValueError, match="_atoms"):
        rs.physical_elemental_inventory(rs.y)


def test_explicit_oligomer_handoff_orphans_nothing():
    """A chain handed to the explicit DP=xs species keeps its end groups (the
    species carries them): the physical balance closes and the ledger books
    only the other mu0 changes, not the handed-off chains."""
    pool_sp = _spc("CCCCCC", "A")
    mu = [_spc("CO", "A_mu0"), _spc("C=O", "A_mu1"), _spc("C#N", "A_mu2")]
    mon = _spc("C=C", "mon")
    oligo = _spc("CCCCCC", "A_x3")    # H-(C2H4)3-H = C6H14
    core = [pool_sp] + mu + [mon, oligo]
    mask = np.array([False] * 4 + [True, False], dtype=bool)
    pool = PolymerPoolConfig(
        label="A", xs=3, explicit_dp_to_species_index={3: 5}, mu_indices=(1, 2, 3),
        monomer_poly_index=4, k_scission=0.0, k_unzip=0.5, tail_kinetics=None,
        repeat_unit_atoms=(("C", 2), ("H", 4)), end_group_atoms=(("H", 1), ("H", 1)))
    rs = HybridPolymerSystem(
        T=800.0, P=1.0e5, initial_mole_fractions={mon: 0.5}, V_poly=1.0,
        polymer_pools=[pool], mass_transfer=[], gas_species_mask=mask.copy(),
        constant_gas_volume=False,
        initial_polymer_moments={"A": (1.0, 12.0, 200.0)}, termination=[])
    rs.initialize_model(core, [], [], [])
    inv0 = rs.physical_elemental_inventory(rs.y0)
    rs.advance(2.0)
    assert rs.y[5] > 1e-6, "handshake must have produced explicit oligomer"
    inv = rs.physical_elemental_inventory(rs.y)
    for el in ("C", "H"):
        assert abs(inv[el] - inv0[el]) <= 1e-10 * inv0[el], (el, inv, inv0)
    # ends of the handed-off chains are NOT in the ledger
    led = rs.end_group_ledger(rs.y)["A"]
    handed = float(rs.y[5])
    assert led == pytest.approx(float(rs.y0[1] - rs.y[1]) - handed, abs=1e-9)
