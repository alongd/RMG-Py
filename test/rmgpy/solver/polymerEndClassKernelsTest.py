"""Red tests (I-266): what MoM needs to carry a kinetic Monte Carlo kernel
set whose chemistry has more than one end type, on the existing moments.

* A join of two pools: a ``kernel_recombination`` term with a ``partner``
  pool takes one chain of this pool (``sites_per_chain`` sites) and one of
  the partner (``partner_sites_per_chain``) into one chain of n + m periods
  in ``product`` at ``k_pair*s*s'*c*c'`` per volume (two end-radical classes
  recombining, an addition joining a radical and an alkene pool).  The
  convolution moments are exact (no closure); mass and the radical count are
  conserved: each event takes one chain of each pool and makes one product.
* ``end_class`` keys an explicit end-radical pool by its class, so two
  end-radical pools of one site are two pools.
* ``k_depropagation`` on an explicit end-radical pool (the deck declares it
  with its ``monomer_product``): the solver's existing depropagation law.
* ``kernel_propagation``: a chain of the pool adds one monomer unit at
  ``k_pair*c_M`` per chain, ``c_M`` the pool's monomer_product amount over
  the polymer volume (kMC's convention: every small species at its melt
  concentration), the monomer consumed.

Every expected value is an explicit sum over the chains of a discrete
distribution, not the moment algebra under test.
"""

import math

import numpy as np
import pytest

from rmgpy.exceptions import InputError
from rmgpy.molecule import Molecule
from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig
from rmgpy.species import Species
import rmgpy.solver.polymer as SP

from polymerKernelsTest import PPU, _spc, moments, r_term

V = 2.0          # m3 of polymer


def _system(pools, chains, *, monomer=None):
    """Core: N2 (gas), then per pool its proxy and three moment slots
    (pool i at 1 + 4i .. 4 + 4i), then the gas monomer M when given."""
    core = [_spc("N#N", "N2")]
    for i, pool in enumerate(pools):
        stem = "C" * (i + 1)
        core += [_spc("C" * (i + 4), pool.label), _spc(stem + "O", f"{pool.label}_mu0"),
                 _spc(stem + "N", f"{pool.label}_mu1"), _spc(stem + "S", f"{pool.label}_mu2")]
    if monomer is not None:
        core.append(_spc("C=C", "M"))
    mask = np.array([s.label in ("N2", "M") for s in core], dtype=bool)
    rs = HybridPolymerSystem(T=700.0, P=1.0e5, initial_mole_fractions={core[0]: 1.0}, V_poly=V,
                             polymer_pools=pools, mass_transfer=[], gas_species_mask=mask,
                             constant_gas_volume=False,
                             initial_polymer_moments={p.label: moments(chains.get(p.label, {})) for p in pools},
                             termination=[])
    rs.initialize_model(core, [], [], [])
    if monomer is not None:
        rs.y[len(core) - 1] = monomer
    return rs


def _pool(label, i, **kw):
    return PolymerPoolConfig(label=label, xs=2, explicit_dp_to_species_index={},
                             mu_indices=(2 + 4 * i, 3 + 4 * i, 4 + 4 * i), **kw)


def _slots(dn, i):
    return np.array(dn[2 + 4 * i:5 + 4 * i])


def join_term(**kw):
    term = dict(template="join_t", k_pair=3.0, sites_per_chain=1.0, partner="B", partner_sites_per_chain=2.0,
                periods_per_unit=PPU, product="P")
    term.update(kw)
    return term


def join_sums(a, b, k, s, s2):
    """d(moments)/dt (pool units, mol/s) of A, B and the product by explicit
    sums over every pair: rate k s s' N_i N_j / V."""
    da, db, dp = np.zeros(3), np.zeros(3), np.zeros(3)
    for i, Ni in a.items():
        for j, Nj in b.items():
            rate = k * s * s2 * Ni * Nj / V
            for m in range(3):
                da[m] -= rate * (i / PPU) ** m
                db[m] -= rate * (j / PPU) ** m
                dp[m] += rate * ((i + j) / PPU) ** m
    return da, db, dp


# --------------------------------------------------------------------------
# joins of two pools
# --------------------------------------------------------------------------
@pytest.mark.parametrize("b_chains", [{1: 0.7}, {3: 0.2, 8: 0.5, 21: 0.1}], ids=["fixed-length", "spread"])
def test_a_join_of_two_pools_is_the_pair_sum(b_chains):
    a_chains = {4: 0.3, 9: 0.6, 30: 0.2}
    pools = [_pool("A", 0, kernel_recombination=[join_term()]), _pool("B", 1), _pool("P", 2)]
    rs = _system(pools, {"A": a_chains, "B": b_chains})
    dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
    da, db, dp = join_sums(a_chains, b_chains, 3.0, 1.0, 2.0)
    assert _slots(dn, 0) == pytest.approx(da, rel=1e-12)
    assert _slots(dn, 1) == pytest.approx(db, rel=1e-12)
    assert _slots(dn, 2) == pytest.approx(dp, rel=1e-12)


def test_a_two_end_class_join_conserves_mass_and_the_radical_count():
    """Each event takes one chain of each end-radical pool and makes one
    product chain; the periods (mass) are conserved."""
    pools = [_pool("A", 0, kernel_recombination=[join_term()]), _pool("B", 1), _pool("P", 2)]
    rs = _system(pools, {"A": {6: 0.4, 40: 0.3}, "B": {1: 0.9}})
    dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
    a, b, p = _slots(dn, 0), _slots(dn, 1), _slots(dn, 2)
    assert a[0] == pytest.approx(b[0], rel=1e-14) and p[0] == pytest.approx(-a[0], rel=1e-14) and p[0] > 0.0
    assert a[0] + b[0] + 2.0 * p[0] == pytest.approx(0.0, abs=1e-14 * p[0])   # radicals: two ends per event
    assert a[1] + b[1] + p[1] == pytest.approx(0.0, abs=1e-13 * p[1])          # units (mass)


def test_a_term_without_a_partner_is_the_same_pool_recombination():
    """The existing law: two chains of this pool."""
    pools = [_pool("A", 0, kernel_recombination=[r_term(product="P")]), _pool("P", 1)]
    chains = {6: 0.4, 40: 0.3}
    rs = _system(pools, {"A": chains})
    dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
    rad, prod = np.zeros(3), np.zeros(3)
    for i, Ni in chains.items():
        for j, Nj in chains.items():
            rate = 5.0 * Ni * Nj / V
            for m in range(3):
                rad[m] -= rate * ((i / PPU) ** m + (j / PPU) ** m)
                prod[m] += rate * ((i + j) / PPU) ** m
    assert _slots(dn, 0) == pytest.approx(rad, rel=1e-12)
    assert _slots(dn, 1) == pytest.approx(prod, rel=1e-12)


@pytest.mark.parametrize("term, match", [
    (join_term(partner_sites_per_chain=None), r"or with all of \['partner', 'partner_sites_per_chain'\]"),
    (join_term(partner=None), r"or with all of \['partner', 'partner_sites_per_chain'\]"),
    (join_term(partner="A"), "partner='A' must name another pool"),
    (join_term(partner_sites_per_chain=0.0), r"partner_sites_per_chain=0.0 must be > 0"),
])
def test_a_malformed_join_is_refused(term, match):
    term = {k: v for k, v in term.items() if v is not None}
    with pytest.raises(ValueError, match=match):
        SP.validate_kernel_recombination("A", [term])


def test_a_join_naming_an_unconfigured_partner_is_refused():
    pools = [_pool("A", 0, kernel_recombination=[join_term(partner="Nowhere")]), _pool("P", 1)]
    with pytest.raises(ValueError, match="partner pool 'Nowhere', which is not among the configured"):
        _system(pools, {"A": {6: 1.0}})


def test_the_join_connects_the_partner_to_the_product():
    pool = _pool("A", 0, kernel_recombination=[join_term()])
    assert ("kernel_recombination", "join_t", "B", "P") in SP.catalog_kernel_connections(pool)
    assert ("kernel_recombination", "join_t", "A", "P") in SP.catalog_kernel_connections(pool)


# --------------------------------------------------------------------------
# propagation
# --------------------------------------------------------------------------
def p_term(**kw):
    term = dict(template="propagation_t", k_pair=0.8)
    term.update(kw)
    return term


def test_propagation_adds_one_unit_per_event_and_consumes_the_monomer():
    chains = {4: 0.3, 9: 0.6, 30: 0.2}
    n_m = 0.05
    pools = [_pool("Re", 0, kernel_propagation=[p_term()], monomer_poly_index=5)]
    rs = _system(pools, {"Re": chains}, monomer=n_m)
    dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
    want, used = np.zeros(3), 0.0
    for n, N in chains.items():
        rate = 0.8 * (n_m / V) * N                  # per chain k_pair c_M, c_M on the polymer volume
        used += rate
        for m in range(3):
            want[m] += rate * ((n / PPU + 1) ** m - (n / PPU) ** m)
    assert _slots(dn, 0) == pytest.approx(want, rel=1e-12)
    assert dn[5] == pytest.approx(-used, rel=1e-12)
    assert _slots(dn, 0)[1] == pytest.approx(-dn[5], rel=1e-12)       # one unit per monomer


def test_propagation_without_a_monomer_target_is_refused():
    pools = [_pool("Re", 0, kernel_propagation=[p_term()])]
    with pytest.raises(ValueError, match="kernel_propagation.*monomer"):
        _system(pools, {"Re": {6: 1.0}})


@pytest.mark.parametrize("term, match", [(p_term(k_pair=0.0), "k_pair"), (dict(p_term(), extra=1), "keys")])
def test_a_malformed_propagation_is_refused(term, match):
    with pytest.raises(ValueError, match=match):
        SP.validate_kernel_propagation("Re", [term])


# --------------------------------------------------------------------------
# the deck surface
# --------------------------------------------------------------------------
def _deck(**kw):
    from unittest.mock import MagicMock
    from rmgpy.rmg import input as rmg_input

    def make(obj, **_k):
        if isinstance(obj, Molecule):
            obj = Species(label="M", molecule=[obj])
        return obj, True
    mock_rmg = MagicMock()
    mock_rmg.initial_species = []
    mock_rmg.reaction_model.iteration_num = 0
    mock_rmg.reaction_model.new_species_list = []
    mock_rmg.reaction_model.make_new_species.side_effect = make
    old_rmg, old_sd = rmg_input.rmg, rmg_input.species_dict
    rmg_input.rmg, rmg_input.species_dict = mock_rmg, {}
    args = dict(label="Re", monomer="[CH2][CH2]", end_groups=["[H]", "[H]"], cutoff=4, moments=[0.0, 0.0, 0.0],
                initial_mass=0.0, end_radical_site="primary")
    args.update(kw)
    try:
        return rmg_input.polymer(**args)
    finally:
        rmg_input.rmg, rmg_input.species_dict = old_rmg, old_sd


KDEP = {"A": 2930.0, "n": 0.0, "Ea": 0.0}


def test_an_explicit_end_radical_pool_carries_depropagation_and_propagation():
    r = _deck(k_depropagation=KDEP, monomer_product="C=C", kernel_propagation=[p_term()])
    assert r.k_depropagation == SP.validate_k_depropagation("Re", KDEP)
    assert r.kernel_propagation == SP.validate_kernel_propagation("Re", [p_term()])
    assert r.monomer_product_species is not None


def test_depropagation_on_a_pool_that_is_no_explicit_end_radical_pool_still_needs_homolysis():
    with pytest.raises(InputError, match="k_homolysis"):
        _deck(end_radical_site=None, k_depropagation=KDEP, monomer_product="C=C")


def test_propagation_without_a_monomer_product_is_refused_at_the_deck():
    with pytest.raises(InputError, match="kernel_propagation.*monomer_product"):
        _deck(kernel_propagation=[p_term()])


def test_an_end_class_keys_a_second_end_radical_pool():
    plain, methyl = _deck(), _deck(label="Re1", end_class="CH3")
    assert methyl.end_class == "CH3" and plain.end_class is None
    assert methyl.fingerprint != plain.fingerprint
    assert methyl.copy(deep=True).fingerprint == methyl.fingerprint


@pytest.mark.parametrize("kw, match", [(dict(end_radical_site=None, end_class="CH3"), "end_class"),
                                       (dict(end_class="C H3"), "end_class"), (dict(end_class=""), "end_class")])
def test_a_malformed_end_class_is_refused(kw, match):
    with pytest.raises(InputError, match=match):
        _deck(**kw)


def test_the_pool_config_carries_an_explicit_pools_depropagation_and_propagation():
    from rmgpy.rmg.polymer_input import PolymerPool
    mus = [_spc("CO", "Re_mu0"), _spc("C=O", "Re_mu1"), _spc("C#N", "Re_mu2")]
    m = _spc("C=C", "M")
    pool = PolymerPool(label="Re", xs=2, monomer=_spc("C=C", "Re_repeat"), explicit_map={}, mu_species=mus,
                       monomer_product=m, k_depropagation=KDEP, kernel_propagation=[p_term()])
    config = pool.to_config({mus[0]: 1, mus[1]: 2, mus[2]: 3, m: 4})
    assert config.k_depropagation == SP.validate_k_depropagation("Re", KDEP)
    assert config.kernel_propagation == SP.validate_kernel_propagation("Re", [p_term()])
    assert config.monomer_poly_index == 4
