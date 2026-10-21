"""Catalog-driven pool kernels (slice S5 a + b of the shared-chemistry design):
backbone homolysis of a dead-chain pool into an explicit end-radical pool,
``kernel_homolysis``, and end-radical recombination back into a dead-chain
pool, ``kernel_recombination``.

Lengths in a kernel term count the term's period, ``periods_per_unit`` of
them per pool monomer unit (the pools keep their monomer unit; the closure
is homogeneous in the length scale).  Homolysis: a chain of n periods has
``a*n + b`` sites, each cut leaving fragments k and n - k with k uniform on
[min_fragment, n - min_fragment] (so b = a*(1 - 2*min_fragment)); event rate
``w*phi*(a*mu1 + b*mu0)``.  Recombination: two chains of the radical pool,
each with ``sites_per_chain`` sites, join into one chain of n0 + n1 periods
in the product pool at ``k_pair * S**2`` per volume (S the radical pool's
site concentration).  Every expected value below is an explicit sum over the
chains of a discrete distribution, not the moment algebra under test.
"""

import math

import numpy as np
import pytest

from rmgpy.exceptions import InputError
import rmgpy.solver.polymer as SP
from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig
from rmgpy.species import Species
from rmgpy.molecule import Molecule

PPU = 2          # periods per monomer unit (CH2 periods of a C2 repeat unit)
E = 4            # min_fragment, periods
A, B = 1.0, 1.0 - 2 * E


def _spc(smiles, label):
    return Species(label=label, molecule=[Molecule(smiles=smiles)])


def h_term(**kw):
    term = dict(template="homolysis_t", w=0.3, phi=0.5, a=A, b=B, min_fragment=E,
                periods_per_unit=PPU, daughter="Re")
    term.update(kw)
    return term


def r_term(**kw):
    term = dict(template="recombination_t", k_pair=5.0, sites_per_chain=1.0,
                periods_per_unit=PPU, product="P")
    term.update(kw)
    return term


def moments(chains):
    """Pool-unit moments (mol) of {length in periods: mol}."""
    return tuple(math.fsum(N * (n / PPU) ** k for n, N in chains.items()) for k in range(3))


def system(p_chains, r_chains, *, homolysis=None, recombination=None, V_poly=1.0,
           p_extra=None):
    """Core: N2 (gas), P proxy + mu slots 2-4, Re proxy + mu slots 6-8."""
    core = [_spc("N#N", "N2"),
            _spc("CCCC", "P"), _spc("CO", "P_mu0"), _spc("C=O", "P_mu1"), _spc("C#N", "P_mu2"),
            _spc("[CH2]CCC", "Re"), _spc("CCO", "Re_mu0"), _spc("CC=O", "Re_mu1"), _spc("CC#N", "Re_mu2")]
    pools = [PolymerPoolConfig(label="P", xs=2, explicit_dp_to_species_index={}, mu_indices=(2, 3, 4),
                               kernel_homolysis=homolysis, **(p_extra or {})),
             PolymerPoolConfig(label="Re", xs=2, explicit_dp_to_species_index={}, mu_indices=(6, 7, 8),
                               kernel_recombination=recombination)]
    mask = np.array([s.label == "N2" for s in core], dtype=bool)
    rs = HybridPolymerSystem(T=700.0, P=1.0e5, initial_mole_fractions={core[0]: 1.0}, V_poly=V_poly,
                             polymer_pools=pools, mass_transfer=[], gas_species_mask=mask,
                             constant_gas_volume=False,
                             initial_polymer_moments={"P": moments(p_chains), "Re": moments(r_chains)},
                             termination=[])
    rs.initialize_model(core, [], [], [])
    return rs


def rates(rs):
    dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
    return np.array(dn[2:5]), np.array(dn[6:9]), dn


def homolysis_sums(chains, w, phi):
    """d(moments)/dt (pool units, mol/s) of parent and daughter by explicit
    sums over every chain and cut."""
    parent, daughter = [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]
    for n, N in chains.items():
        for k in range(E, n - E + 1):
            rate = w * phi * A * N                  # one cut position, a sites
            for j in range(3):
                parent[j] -= rate * (n / PPU) ** j
                daughter[j] += rate * ((k / PPU) ** j + ((n - k) / PPU) ** j)
    return np.array(parent), np.array(daughter)


def recombination_sums(chains, k_pair, s, V):
    """d(moments)/dt of the radical and product pool by explicit sums over
    every ordered site pair: rate k_pair*s^2*N_i*N_j/V (mol/s)."""
    rad, prod = [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]
    for i, Ni in chains.items():
        for j, Nj in chains.items():
            rate = k_pair * s * s * Ni * Nj / V
            for m in range(3):
                rad[m] -= rate * ((i / PPU) ** m + (j / PPU) ** m)
                prod[m] += rate * ((i + j) / PPU) ** m
    return np.array(rad), np.array(prod)


# --------------------------------------------------------------------------
# the moment laws
# --------------------------------------------------------------------------
def test_homolysis_moment_law_is_the_cut_sum_on_a_monodisperse_pool():
    chains = {60: 1.5}                              # the log-Lagrange mu3 is exact here
    rs = system(chains, {}, homolysis=[h_term()])
    p, r, dn = rates(rs)
    p_exp, r_exp = homolysis_sums(chains, 0.3, 0.5)
    assert p == pytest.approx(p_exp, rel=1e-12)
    assert r == pytest.approx(r_exp, rel=1e-12)
    assert p[1] + r[1] == pytest.approx(0.0, abs=1e-12 * abs(p[1]))   # periods are conserved
    assert r[0] == pytest.approx(-2.0 * p[0], rel=1e-14)               # two fragments per cut
    assert dn[0] == 0.0 and dn[1] == 0.0 and dn[5] == 0.0


def test_homolysis_rate_counts_the_end_exclusion():
    """At n = 9 periods only k = 4, 5 cut: 2 sites, not the 9 - 1 bonds."""
    chains = {9: 1.0}
    rs = system(chains, {}, homolysis=[h_term(phi=1.0)])
    p, _r, _dn = rates(rs)
    assert p[0] == pytest.approx(-0.3 * 2.0, rel=1e-14)


def test_homolysis_sums_its_templates_each_with_its_own_phi():
    chains = {60: 1.0}
    rs = system(chains, {}, homolysis=[h_term(phi=0.25), h_term(template="other", w=0.1, phi=1.0)])
    p, r, _dn = rates(rs)
    p1, r1 = homolysis_sums(chains, 0.3, 0.25)
    p2, r2 = homolysis_sums(chains, 0.1, 1.0)
    assert p == pytest.approx(p1 + p2, rel=1e-12)
    assert r == pytest.approx(r1 + r2, rel=1e-12)


def test_recombination_moment_law_is_the_pair_sum():
    chains = {8: 0.05, 14: 0.12, 31: 0.03}
    V = 2.0
    rs = system({60: 1.0}, chains, recombination=[r_term(sites_per_chain=1.5)], V_poly=V)
    p, r, _dn = rates(rs)
    r_exp, p_exp = recombination_sums(chains, 5.0, 1.5, V)
    assert r == pytest.approx(r_exp, rel=1e-12)
    assert p == pytest.approx(p_exp, rel=1e-12)
    assert p[1] + r[1] == pytest.approx(0.0, abs=1e-12 * abs(p[1]))
    assert r[0] == pytest.approx(-2.0 * p[0], rel=1e-14)


def test_recombination_gains_mu2_beyond_the_transfer():
    """The product's mu2 gain exceeds what the radical pool loses by
    2*k_pair*s^2*mu1^2/V: the cross term of (n0 + n1)^2."""
    chains = {10: 0.1, 20: 0.1}
    rs = system({60: 1.0}, chains, recombination=[r_term()], V_poly=1.0)
    p, r, _dn = rates(rs)
    mu1 = moments(chains)[1]
    assert p[2] + r[2] == pytest.approx(2.0 * 5.0 * mu1 * mu1, rel=1e-12)


# --------------------------------------------------------------------------
# exact trajectories
# --------------------------------------------------------------------------
def test_homolysis_from_a_monodisperse_pool_follows_its_exact_solution():
    """Parent chains of n0 periods stay monodisperse (a cut removes a whole
    chain), so N(t) = N0 exp(-w phi s(n0) t) and every daughter moment is
    (N0 - N) times the per-cut fragment moment."""
    n0, N0, w, phi = 40, 1.0, 0.3, 0.5
    s = A * n0 + B
    frag = [math.fsum((k / PPU) ** j + ((n0 - k) / PPU) ** j for k in range(E, n0 - E + 1)) / (n0 - 2 * E + 1)
            for j in range(3)]
    rs = system({n0: N0}, {}, homolysis=[h_term(w=w, phi=phi)])
    for t in (0.05, 0.2, 0.6):
        rs.advance(t)
        y = np.asarray(rs.y)
        N = N0 * math.exp(-w * phi * s * t)
        assert y[2:5] == pytest.approx([N * (n0 / PPU) ** j for j in range(3)], rel=1e-5)
        assert y[6:9] == pytest.approx([(N0 - N) * f for f in frag], rel=1e-5)


def test_recombination_alone_follows_its_exact_solution():
    """Random pairing leaves the radical pool's shape: mu_j(t) = mu_j(0) f,
    f = 1/(1 + alpha t), alpha = 2 k s^2 mu0(0)/V; the product's mu2 gains
    (k s^2/V)(2 mu0 mu2 + 2 mu1^2)(0) t/(1 + alpha t)."""
    chains = {6: 0.2, 11: 0.3, 25: 0.1}
    k, s, V = 0.4, 1.0, 1.0
    m0 = moments(chains)
    rs = system({}, chains, recombination=[r_term(k_pair=k, sites_per_chain=s)], V_poly=V)
    alpha = 2.0 * k * s * s * m0[0] / V
    for t in (0.5, 2.0, 6.0):
        rs.advance(t)
        y = np.asarray(rs.y)
        f = 1.0 / (1.0 + alpha * t)
        assert y[6:9] == pytest.approx([m * f for m in m0], rel=1e-5)
        mu2_p = (k * s * s / V) * (2.0 * m0[0] * m0[2] + 2.0 * m0[1] ** 2) * t / (1.0 + alpha * t)
        assert y[2:5] == pytest.approx([m0[0] * (1 - f) / 2.0, m0[1] * (1 - f), mu2_p], rel=1e-5)


# --------------------------------------------------------------------------
# nulls
# --------------------------------------------------------------------------
def test_recombination_only_without_radicals_is_flat():
    rs = system({60: 1.0, 70: 0.5}, {}, recombination=[r_term()])
    y0 = np.array(rs.y)
    assert np.all(rates(rs)[2] == 0.0)
    rs.advance(10.0)
    assert np.array_equal(np.asarray(rs.y), y0)


def test_phi_zero_is_flat():
    rs = system({60: 1.0, 70: 0.5}, {}, homolysis=[h_term(phi=0.0)], recombination=[r_term()])
    y0 = np.array(rs.y)
    assert np.all(rates(rs)[2] == 0.0)
    rs.advance(10.0)
    assert np.array_equal(np.asarray(rs.y), y0)


# --------------------------------------------------------------------------
# refusals
# --------------------------------------------------------------------------
@pytest.mark.parametrize("change, match", [
    (dict(b=0.0), r"b = a\*\(1 - 2\*min_fragment\)"),
    (dict(phi=1.5), "phi=1.5"),
    (dict(phi=-0.1), "phi=-0.1"),
    (dict(w=0.0), "w=0.0"),
    (dict(a=-1.0, b=7.0), "a=-1.0"),
    (dict(min_fragment=0, b=1.0), "min_fragment=0"),
    (dict(periods_per_unit=0), "periods_per_unit=0"),
    (dict(daughter="P"), "daughter"),
])
def test_a_malformed_homolysis_term_is_refused(change, match):
    with pytest.raises(ValueError, match=match):
        SP.validate_kernel_homolysis("P", [h_term(**change)])


def test_homolysis_terms_need_exactly_the_keys_and_distinct_templates():
    term = h_term()
    del term["phi"]
    with pytest.raises(ValueError, match="keys"):
        SP.validate_kernel_homolysis("P", [term])
    with pytest.raises(ValueError, match="duplicate"):
        SP.validate_kernel_homolysis("P", [h_term(), h_term()])
    with pytest.raises(ValueError, match="list"):
        SP.validate_kernel_homolysis("P", h_term())


@pytest.mark.parametrize("change, match", [
    (dict(k_pair=0.0), "k_pair=0.0"),
    (dict(sites_per_chain=-1.0), "sites_per_chain=-1.0"),
    (dict(product="Re"), "product"),
])
def test_a_malformed_recombination_term_is_refused(change, match):
    with pytest.raises(ValueError, match=match):
        SP.validate_kernel_recombination("Re", [r_term(**change)])


def test_a_kernel_naming_an_unconfigured_pool_is_refused():
    with pytest.raises(ValueError, match="'Nowhere'"):
        system({60: 1.0}, {}, homolysis=[h_term(daughter="Nowhere")])
    with pytest.raises(ValueError, match="'Nowhere'"):
        system({60: 1.0}, {}, recombination=[r_term(product="Nowhere")])


def test_kernel_homolysis_and_k_homolysis_on_one_pool_are_refused():
    with pytest.raises(ValueError, match="kernel_homolysis.*k_homolysis"):
        system({60: 1.0}, {}, homolysis=[h_term()],
               p_extra=dict(k_homolysis=dict(A=1.0e13, n=0.0, Ea=1.0e5)))


# --------------------------------------------------------------------------
# the deck surface
# --------------------------------------------------------------------------
def _deck_polymers(**kw):
    """polymer() twice under a mocked RMG: a dead-chain pool and its explicit
    end-radical pool."""
    from unittest.mock import MagicMock
    from rmgpy.rmg import input as rmg_input

    mock_rmg = MagicMock()
    mock_rmg.initial_species = []
    mock_rmg.reaction_model.iteration_num = 0
    mock_rmg.reaction_model.new_species_list = []
    mock_rmg.reaction_model.make_new_species.side_effect = lambda obj, **_k: (obj, True)
    old_rmg, old_sd = rmg_input.rmg, rmg_input.species_dict
    rmg_input.rmg, rmg_input.species_dict = mock_rmg, {}
    try:
        p = rmg_input.polymer(label="P", monomer="[CH2][CH2]", end_groups=["[H]", "[H]"], cutoff=4,
                              moments=[1.0, 30.0, 900.0], initial_mass=0.001,
                              kernel_homolysis=kw.get("homolysis", [h_term()]))
        r = rmg_input.polymer(label="Re", monomer="[CH2][CH2]", end_groups=["[H]", "[H]"], cutoff=4,
                              moments=[0.0, 0.0, 0.0], initial_mass=0.0, end_radical_site="primary",
                              kernel_recombination=kw.get("recombination", [r_term()]))
    finally:
        rmg_input.rmg, rmg_input.species_dict = old_rmg, old_sd
    return p, r


def test_the_deck_declares_an_explicit_end_radical_pool_with_its_kernels():
    p, r = _deck_polymers()
    assert p.kernel_homolysis == SP.validate_kernel_homolysis("P", [h_term()])
    assert r.kernel_recombination == SP.validate_kernel_recombination("Re", [r_term()])
    assert r.end_radical_site == "primary"
    assert r.molecule[0].multiplicity == 2 and p.molecule[0].multiplicity == 1
    copied = p.copy(deep=True)
    assert copied.kernel_homolysis == p.kernel_homolysis
    assert r.copy(deep=True).kernel_recombination == r.kernel_recombination


def test_the_deck_refuses_a_malformed_kernel():
    with pytest.raises(InputError, match=r"b = a\*\(1 - 2\*min_fragment\)"):
        _deck_polymers(homolysis=[h_term(b=0.0)])


def test_the_artifact_records_the_kernels_and_a_consumer_without_them_refuses_it():
    from rmgpy.polymer import build_polymer_moments_artifact
    from rmgpy.tools.polymer_moments_runner import _check_schema_version_known
    p, r = _deck_polymers()
    artifact = build_polymer_moments_artifact([p, r], core_species=None, core_reactions=[],
                                              configured_pool_labels=["P", "Re"], monomer_routing_by_pool={})
    entries = {e["label"]: e for e in artifact["pools"]}
    assert entries["P"]["catalog_kernels"] == {"homolysis": p.kernel_homolysis}
    assert entries["Re"]["catalog_kernels"] == {"recombination": r.kernel_recombination}
    assert artifact["schema_version"] == "3.2"
    with pytest.raises(Exception, match="3.2"):
        _check_schema_version_known(artifact)
    plain = build_polymer_moments_artifact([_deck_polymers(homolysis=None, recombination=None)[0]],
                                           core_species=None, core_reactions=[], configured_pool_labels=["P"],
                                           monomer_routing_by_pool={})
    assert "catalog_kernels" not in plain["pools"][0] and plain["schema_version"] != "3.2"


def test_a_pool_configured_without_initial_moles_keeps_its_kernels():
    """A born-empty end-radical pool has no initialMoles entry, so it is
    configured by derive_daughter_pool_configs, not compile_polymer_phase:
    its kernels must survive that path (a dropped term is a silent no-op)."""
    from rmgpy.rmg.polymer_input import derive_daughter_pool_configs
    p, r = _deck_polymers()
    dummies = []
    for k in range(3):
        s = _spc("[Ne]", f"Re_mu{k}")
        s.reactive = False
        s.is_moment_dummy = True
        dummies.append(s)
    core = [_spc("CCCC", "P"), r] + dummies
    (derived,) = derive_daughter_pool_configs(core, {s: i for i, s in enumerate(core)}, existing_pool_labels={"P"})
    assert derived.label == "Re"
    assert derived.kernel_recombination == r.kernel_recombination
    assert derived.kernel_homolysis is None


def test_the_pool_config_carries_the_kernels():
    from rmgpy.rmg.polymer_input import PolymerPool
    mu = [_spc("CO", "P_mu0"), _spc("C=O", "P_mu1"), _spc("C#N", "P_mu2")]
    proxy = _spc("CCCC", "P")
    pool = PolymerPool(label="P", xs=4, monomer=_spc("C=C", "P_repeat"), explicit_map={}, mu_species=mu,
                       proxy_species=proxy, kernel_homolysis=[h_term()], kernel_recombination=None)
    spc_map = {s: i for i, s in enumerate([proxy] + mu)}
    config = pool.to_config(spc_map)
    assert config.kernel_homolysis == SP.validate_kernel_homolysis("P", [h_term()])
    assert config.kernel_recombination is None
