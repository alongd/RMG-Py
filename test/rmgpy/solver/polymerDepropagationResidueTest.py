"""Red tests (I-286): the end groups a depropagation leaves behind.

``k_depropagation`` releases one monomer per unzip event (dmu1 = -R, the
monomer +R) and drains the chains it leaves with no unit, dmu0 = -k_dep N1
(N1 = mu0 P(DP = 1), the existing closure).  What is left of such a chain
is its end groups: in a kernel-mode deck a pool's chain weighs M_u n + e_p,
so the drain removes e_p per chain from the system unless something
receives it.  ``depropagation_residue`` names that small species: the drain
credits it one molecule per drained chain, so the pool's chain count and
the residue's amount move by the same float (mass conserved when the
residue weighs e_p).  Without it the drain is what it was.
"""

import numpy as np
import pytest

from rmgpy.exceptions import InputError
from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig

from polymerEndClassKernelsTest import KDEP, _deck
from polymerKernelsTest import _spc

V = 2.0
# short chains (mol of chains by DP in units): a sizeable DP-1 fraction
CHAINS = {1: 0.4, 2: 0.3, 5: 0.2}
MU = tuple(sum(N * n ** k for n, N in CHAINS.items()) for k in range(3))


def _core(residue=True):
    """N2 (gas), the pool's proxy and mu slots 2-4, the gas monomer M (5), the residue R (6)."""
    core = [_spc("N#N", "N2"), _spc("CCCC", "Re"), _spc("CO", "Re_mu0"), _spc("CN", "Re_mu1"),
            _spc("CS", "Re_mu2"), _spc("C=C", "M")]
    return core + [_spc("[H]", "R")] if residue else core


def _system(core=None, **kw):
    core = core or _core()
    kw.setdefault("monomer_poly_index", 5)
    kw.setdefault("k_depropagation", KDEP)
    pool = PolymerPoolConfig(label="Re", xs=2, explicit_dp_to_species_index={}, mu_indices=(2, 3, 4), **kw)
    mask = np.array([s.label in ("N2", "M") for s in core], dtype=bool)
    rs = HybridPolymerSystem(T=700.0, P=1.0e5, initial_mole_fractions={core[0]: 1.0}, V_poly=V,
                             polymer_pools=[pool], mass_transfer=[], gas_species_mask=mask,
                             constant_gas_volume=False, initial_polymer_moments={"Re": MU}, termination=[])
    rs.initialize_model(core, [], [], [])
    return rs


def _rates(rs):
    return rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]


def test_the_drained_chains_end_groups_go_to_the_residue():
    dn = _rates(_system(depropagation_residue_index=6))
    assert dn[2] < 0.0                                  # chains leave
    assert dn[6] == -dn[2]                              # one residue per drained chain, the same float
    assert dn[5] == -dn[3]                              # one monomer per unit, as before


def test_without_a_residue_the_drain_is_unchanged():
    a = _rates(_system(depropagation_residue_index=6))
    b = _rates(_system())
    assert b[6] == 0.0
    assert list(a[:6]) == list(b[:6])


def test_a_residue_without_depropagation_is_refused():
    with pytest.raises(ValueError, match="depropagation_residue.*k_depropagation"):
        _system(k_depropagation=None, monomer_poly_index=None, depropagation_residue_index=6)


@pytest.mark.parametrize("index", [7, -1, 2, 5])
def test_a_residue_off_the_core_on_a_moment_slot_or_on_the_monomer_is_refused(index):
    with pytest.raises(ValueError, match="depropagation_residue"):
        _system(depropagation_residue_index=index)


def test_the_deck_registers_the_residue_species():
    r = _deck(k_depropagation=KDEP, monomer_product="C=C", depropagation_residue="[H]")
    assert r.depropagation_residue_species is not None
    assert r.depropagation_residue_species.molecule[0].to_smiles() == "[H]"
    assert _deck(k_depropagation=KDEP, monomer_product="C=C").depropagation_residue_species is None


def test_a_deck_residue_without_depropagation_is_refused():
    with pytest.raises(InputError, match="depropagation_residue.*k_depropagation"):
        _deck(depropagation_residue="[H]")


def test_a_copy_of_the_pool_keeps_its_residue_by_reference():
    """The core holds a copy of the deck's Polymer, and the residue resolves
    by object (spc_map): a copy that dropped it would wire no credit."""
    r = _deck(k_depropagation=KDEP, monomer_product="C=C", depropagation_residue="[H]")
    assert r.copy(deep=True).depropagation_residue_species is r.depropagation_residue_species


def test_a_born_empty_pool_keeps_its_residue():
    """A born-empty end-radical pool (no initialMoles entry) is configured
    by derive_daughter_pool_configs, not compile_polymer_phase: its residue
    must resolve on that path too."""
    from rmgpy.rmg.polymer_input import derive_daughter_pool_configs
    r = _deck(k_depropagation=KDEP, monomer_product="C=C", depropagation_residue="[H]")
    dummies = []
    for k in range(3):
        s = _spc("[Ne]", f"Re_mu{k}")
        s.reactive = False
        s.is_moment_dummy = True
        dummies.append(s)
    core = [r, r.monomer_product_species, r.depropagation_residue_species] + dummies
    (derived,) = derive_daughter_pool_configs(core, {s: i for i, s in enumerate(core)}, existing_pool_labels=set())
    assert derived.label == "Re" and derived.k_depropagation is not None
    assert derived.depropagation_residue_index == 2
    off = core[:2] + core[3:]                         # the residue is not a core species
    with pytest.raises(ValueError, match="depropagation_residue"):
        derive_daughter_pool_configs(off, {s: i for i, s in enumerate(off)}, existing_pool_labels=set())


def test_the_pool_config_resolves_the_residue_index():
    from rmgpy.rmg.polymer_input import PolymerPool
    mus = [_spc("CO", "Re_mu0"), _spc("C=O", "Re_mu1"), _spc("C#N", "Re_mu2")]
    m, r = _spc("C=C", "M"), _spc("[H]", "R")
    pool = PolymerPool(label="Re", xs=2, monomer=_spc("C=C", "Re_repeat"), explicit_map={}, mu_species=mus,
                       monomer_product=m, k_depropagation=KDEP, depropagation_residue=r)
    assert pool.to_config({mus[0]: 1, mus[1]: 2, mus[2]: 3, m: 4, r: 5}).depropagation_residue_index == 5
    with pytest.raises(ValueError, match="depropagation_residue"):
        pool.to_config({mus[0]: 1, mus[1]: 2, mus[2]: 3, m: 4})
