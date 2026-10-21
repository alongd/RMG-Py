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
site concentration).

Pair transfer (slice S5 c + e), ``kernel_transfer``: one chain of this pool
and one of a ``partner`` pool (this pool too) change pool, lengths unchanged,
into ``product`` and ``partner_product`` at ``k_pair*phi*N*N_partner/V``
(sites folded into k_pair): end + end disproportionation R_e + R_e -> P + U_t
on the end-radical pool, and its reverse U_t + P -> R_e + R_e (x phi) on an
explicit terminal-alkene pool (``end_alkene_site``).

Mid-chain radicals (slice S5 d): explicit interior-radical and internal-alkene
pools (``interior_feature``); pair transfer whose sites grow with a chain
(``sites`` [c00, c10, c01, c11]: rate ``k_pair*phi*(c00 + c10 n + c01 m + c11
n m)`` per ordered pair, n, m in periods), e.g. H transfer from a dead chain
picked by its sites; and ``kernel_move``, a chain changing pool with its length
at ``k*(a n + b)`` per chain (backbiting and its reverse).

Every expected value below is an explicit sum over the chains of a discrete
distribution, not the moment algebra under test.
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


# --------------------------------------------------------------------------
# pair transfer (S5 c + e): end + end disproportionation and its reverse
# --------------------------------------------------------------------------
def t_term(**kw):
    """End + end disproportionation on the end-radical pool: the acceptor
    becomes a dead chain, the donor a terminal alkene."""
    term = dict(template="disproportionation_t", k_pair=3.0, phi=1.0, partner="Re", product="P",
                partner_product="Ut")
    term.update(kw)
    return term


def rev_term(**kw):
    """Its reverse on the terminal-alkene pool: an alkene and a dead chain
    both become end radicals, x phi."""
    term = dict(template="reverse_t", k_pair=0.7, phi=0.4, partner="P", product="Re", partner_product="Re")
    term.update(kw)
    return term


def system3(p_chains, r_chains, u_chains, *, r_transfer=None, u_transfer=None, recombination=None, V_poly=1.0):
    """Core: N2, P (slots 2-4), Re (6-8), Ut (10-12)."""
    core = [_spc("N#N", "N2"),
            _spc("CCCC", "P"), _spc("CO", "P_mu0"), _spc("C=O", "P_mu1"), _spc("C#N", "P_mu2"),
            _spc("[CH2]CCC", "Re"), _spc("CCO", "Re_mu0"), _spc("CC=O", "Re_mu1"), _spc("CC#N", "Re_mu2"),
            _spc("C=CCC", "Ut"), _spc("CCCO", "Ut_mu0"), _spc("CCC=O", "Ut_mu1"), _spc("CCC#N", "Ut_mu2")]
    pools = [PolymerPoolConfig(label="P", xs=2, explicit_dp_to_species_index={}, mu_indices=(2, 3, 4)),
             PolymerPoolConfig(label="Re", xs=2, explicit_dp_to_species_index={}, mu_indices=(6, 7, 8),
                               kernel_recombination=recombination, kernel_transfer=r_transfer),
             PolymerPoolConfig(label="Ut", xs=2, explicit_dp_to_species_index={}, mu_indices=(10, 11, 12),
                               kernel_transfer=u_transfer)]
    mask = np.array([s.label == "N2" for s in core], dtype=bool)
    rs = HybridPolymerSystem(T=700.0, P=1.0e5, initial_mole_fractions={core[0]: 1.0}, V_poly=V_poly,
                             polymer_pools=pools, mass_transfer=[], gas_species_mask=mask,
                             constant_gas_volume=False,
                             initial_polymer_moments={"P": moments(p_chains), "Re": moments(r_chains),
                                                      "Ut": moments(u_chains)},
                             termination=[])
    rs.initialize_model(core, [], [], [])
    return rs


def rates3(rs):
    dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
    return np.array(dn[2:5]), np.array(dn[6:9]), np.array(dn[10:13]), dn


def pair_sums(owner, partner, k, V, same):
    """d(moments)/dt of (owner, partner, product, partner_product) by explicit
    sums over every ordered (owner chain, partner chain) pair, rate
    k*N_i*N_j/V, lengths unchanged."""
    out = [np.zeros(3) for _ in range(4)]
    for i, Ni in owner.items():
        for j, Nj in partner.items():
            rate = k * Ni * Nj / V
            for m in range(3):
                out[0][m] -= rate * (i / PPU) ** m
                out[1][m] -= rate * (j / PPU) ** m
                out[2][m] += rate * (i / PPU) ** m
                out[3][m] += rate * (j / PPU) ** m
    if same:
        out[0], out[1] = out[0] + out[1], out[0] + out[1]
    return out


def test_disproportionation_moment_law_is_the_pair_sum():
    """Every ordered pair of end radicals: the acceptor's length into P, the
    donor's into Ut (lengths unchanged, so no closure)."""
    chains = {8: 0.05, 14: 0.12, 31: 0.03}
    V = 2.0
    rs = system3({}, chains, {}, r_transfer=[t_term()], V_poly=V)
    p, r, u, dn = rates3(rs)
    loss, _same, to_p, to_u = pair_sums(chains, chains, 3.0, V, same=True)
    assert r == pytest.approx(loss, rel=1e-12)
    assert p == pytest.approx(to_p, rel=1e-12)
    assert u == pytest.approx(to_u, rel=1e-12)
    assert u[1] > 0.0 and u == pytest.approx(p, rel=1e-12)       # one alkene per dead chain
    assert p[1] + r[1] + u[1] == pytest.approx(0.0, abs=1e-12 * abs(r[1]))
    assert dn[0] == 0.0


def test_reverse_disproportionation_moment_law_is_the_pair_sum_times_phi():
    """Every (alkene, dead chain) pair at k*phi: both become end radicals."""
    alkenes, dead = {10: 0.02, 26: 0.05}, {40: 0.3, 60: 0.1}
    V = 1.5
    rs = system3(dead, {}, alkenes, u_transfer=[rev_term()], V_poly=V)
    p, r, u, _dn = rates3(rs)
    u_loss, p_loss, from_u, from_p = pair_sums(alkenes, dead, 0.7 * 0.4, V, same=False)
    assert u == pytest.approx(u_loss, rel=1e-12)
    assert p == pytest.approx(p_loss, rel=1e-12)
    assert r == pytest.approx(from_u + from_p, rel=1e-12)
    assert p[1] + r[1] + u[1] == pytest.approx(0.0, abs=1e-12 * abs(r[1]))


def test_transfer_sums_its_terms():
    chains = {9: 0.1, 20: 0.2}
    rs = system3({}, chains, {}, r_transfer=[t_term(), t_term(template="other", k_pair=1.0)])
    _p, r, _u, _dn = rates3(rs)
    assert r == pytest.approx(pair_sums(chains, chains, 4.0, 1.0, same=True)[0], rel=1e-12)


def test_disproportionation_alone_follows_its_exact_solution():
    """Random pairing leaves the radical pool's shape: mu_j(t) = mu_j(0) f,
    f = 1/(1 + alpha t), alpha = 2 k mu0(0)/V; P and Ut each hold half of
    what left, with the radicals' shape."""
    chains = {6: 0.2, 11: 0.3, 25: 0.1}
    k, V = 0.4, 1.0
    m0 = moments(chains)
    rs = system3({}, chains, {}, r_transfer=[t_term(k_pair=k)], V_poly=V)
    alpha = 2.0 * k * m0[0] / V
    for t in (0.5, 2.0, 6.0):
        rs.advance(t)
        y = np.asarray(rs.y)
        f = 1.0 / (1.0 + alpha * t)
        assert y[6:9] == pytest.approx([m * f for m in m0], rel=1e-5)
        assert y[2:5] == pytest.approx([m * (1.0 - f) / 2.0 for m in m0], rel=1e-5)
        assert y[10:13] == pytest.approx([m * (1.0 - f) / 2.0 for m in m0], rel=1e-5)


def test_reverse_disproportionation_alone_follows_its_exact_solution():
    """dN_u/dt = dN_p/dt = -k phi N_u N_p / V, every chain alike: Δ = N_p - N_u
    is constant, N_u(t) = Δ N_u0 / (N_p0 exp(k phi Δ t / V) - N_u0), and each
    pool keeps its shape."""
    alkenes, dead = {10: 0.02, 26: 0.05}, {40: 0.3, 60: 0.1}
    k, phi, V = 2.0, 0.5, 1.0
    mu, mp = moments(alkenes), moments(dead)
    rs = system3(dead, {}, alkenes, u_transfer=[rev_term(k_pair=k, phi=phi)], V_poly=V)
    d = mp[0] - mu[0]
    for t in (1.0, 5.0, 20.0):
        rs.advance(t)
        y = np.asarray(rs.y)
        nu = d * mu[0] / (mp[0] * math.exp(k * phi * d * t / V) - mu[0])
        npp = nu + d
        assert y[10:13] == pytest.approx([m * nu / mu[0] for m in mu], rel=1e-5)
        assert y[2:5] == pytest.approx([m * npp / mp[0] for m in mp], rel=1e-5)
        gained = [mu[j] * (1 - nu / mu[0]) + mp[j] * (1 - npp / mp[0]) for j in range(3)]
        assert y[6:9] == pytest.approx(gained, rel=1e-5)


def test_no_radicals_and_no_alkenes_is_flat():
    rs = system3({60: 1.0, 70: 0.5}, {}, {}, r_transfer=[t_term()], u_transfer=[rev_term()],
                 recombination=[r_term()])
    y0 = np.array(rs.y)
    assert np.all(rates3(rs)[3] == 0.0)
    rs.advance(10.0)
    assert np.array_equal(np.asarray(rs.y), y0)


def test_phi_zero_with_only_the_reverse_is_flat():
    rs = system3({60: 1.0, 70: 0.5}, {}, {20: 0.3, 30: 0.1}, u_transfer=[rev_term(phi=0.0)])
    y0 = np.array(rs.y)
    assert np.all(rates3(rs)[3] == 0.0)
    rs.advance(10.0)
    assert np.array_equal(np.asarray(rs.y), y0)


@pytest.mark.parametrize("change, match", [
    (dict(k_pair=0.0), "k_pair=0.0"),
    (dict(phi=1.5), "phi=1.5"),
    (dict(product="Re"), "product"),
    (dict(partner=""), "partner"),
    (dict(partner_product="Re"), "partner_product"),
])
def test_a_malformed_transfer_term_is_refused(change, match):
    with pytest.raises(ValueError, match=match):
        SP.validate_kernel_transfer("Re", [t_term(**change)])


def test_transfer_terms_need_exactly_the_keys():
    term = t_term()
    del term["phi"]
    with pytest.raises(ValueError, match="keys"):
        SP.validate_kernel_transfer("Re", [term])


@pytest.mark.parametrize("key", ["partner", "product", "partner_product"])
def test_a_transfer_naming_an_unconfigured_pool_is_refused(key):
    with pytest.raises(ValueError, match="'Nowhere'"):
        system3({}, {10: 1.0}, {}, r_transfer=[t_term(**{key: "Nowhere"})])


def _deck_alkene_pool(**kw):
    from unittest.mock import MagicMock
    from rmgpy.rmg import input as rmg_input

    mock_rmg = MagicMock()
    mock_rmg.initial_species = []
    mock_rmg.reaction_model.iteration_num = 0
    mock_rmg.reaction_model.new_species_list = []
    mock_rmg.reaction_model.make_new_species.side_effect = lambda obj, **_k: (obj, True)
    old_rmg, old_sd = rmg_input.rmg, rmg_input.species_dict
    rmg_input.rmg, rmg_input.species_dict = mock_rmg, {}
    args = dict(label="Ut", monomer="[CH2][CH2]", end_groups=["[H]", "[H]"], cutoff=4, moments=[0.0, 0.0, 0.0],
                initial_mass=0.0, end_alkene_site="primary", kernel_transfer=[rev_term()])
    args.update(kw)
    try:
        return rmg_input.polymer(**args)
    finally:
        rmg_input.rmg, rmg_input.species_dict = old_rmg, old_sd


def test_the_deck_declares_a_terminal_alkene_pool_with_its_kernel():
    """Its proxy is the end-radical oligomer with the radical end turned into
    a terminal C=C (the disproportionation donor's product), so the pool is
    distinct from the dead-chain and end-radical pools of the same monomer."""
    u = _deck_alkene_pool()
    p, r = _deck_polymers()
    assert u.kernel_transfer == SP.validate_kernel_transfer("Ut", [rev_term()])
    assert u.end_alkene_site == "primary"
    assert u.molecule[0].is_isomorphic(Molecule(smiles="C=CCCCC"))
    assert len({u.fingerprint, p.fingerprint, r.fingerprint}) == 3
    copied = u.copy(deep=True)
    assert copied.kernel_transfer == u.kernel_transfer and copied.end_alkene_site == "primary"
    assert copied.fingerprint == u.fingerprint


@pytest.mark.parametrize("change, match", [
    (dict(end_alkene_site="middle"), "end_alkene_site"),
    (dict(end_radical_site="primary"), "end_alkene_site.*end_radical_site"),
    (dict(kernel_transfer=[rev_term(phi=2.0)]), "phi=2.0"),
])
def test_the_deck_refuses_a_malformed_alkene_pool(change, match):
    with pytest.raises(InputError, match=match):
        _deck_alkene_pool(**change)


def test_the_artifact_records_the_transfer_kernels_at_3_3():
    from rmgpy.polymer import build_polymer_moments_artifact
    from rmgpy.tools.polymer_moments_runner import _check_schema_version_known
    p, r = _deck_polymers()
    u = _deck_alkene_pool()
    artifact = build_polymer_moments_artifact([p, r, u], core_species=None, core_reactions=[],
                                              configured_pool_labels=["P", "Re", "Ut"], monomer_routing_by_pool={})
    entries = {e["label"]: e for e in artifact["pools"]}
    assert entries["Ut"]["catalog_kernels"] == {"transfer": u.kernel_transfer}
    assert artifact["schema_version"] == "3.3"
    with pytest.raises(Exception, match="3.3"):
        _check_schema_version_known(artifact)


def test_a_born_empty_alkene_pool_keeps_its_transfer_kernel():
    from rmgpy.rmg.polymer_input import derive_daughter_pool_configs
    u = _deck_alkene_pool()
    dummies = []
    for k in range(3):
        s = _spc("[Ne]", f"Ut_mu{k}")
        s.reactive = False
        s.is_moment_dummy = True
        dummies.append(s)
    core = [_spc("CCCC", "P"), u] + dummies
    (derived,) = derive_daughter_pool_configs(core, {s: i for i, s in enumerate(core)}, existing_pool_labels={"P"})
    assert derived.label == "Ut"
    assert derived.kernel_transfer == u.kernel_transfer


def test_the_pool_config_carries_the_transfer_kernel():
    from rmgpy.rmg.polymer_input import PolymerPool
    mu = [_spc("CO", "Re_mu0"), _spc("C=O", "Re_mu1"), _spc("C#N", "Re_mu2")]
    proxy = _spc("[CH2]CCC", "Re")
    pool = PolymerPool(label="Re", xs=4, monomer=_spc("C=C", "Re_repeat"), explicit_map={}, mu_species=mu,
                       proxy_species=proxy, kernel_transfer=[t_term()])
    config = pool.to_config({s: i for i, s in enumerate([proxy] + mu)})
    assert config.kernel_transfer == SP.validate_kernel_transfer("Re", [t_term()])


# --------------------------------------------------------------------------
# mid-chain radicals (S5 d): sites that grow with a chain, pool moves
# --------------------------------------------------------------------------
LABELS5 = ("P", "Re", "Ut", "Rm", "Ui")
PROXIES5 = {"P": "CCCC", "Re": "[CH2]CCC", "Ut": "C=CCC", "Rm": "C[CH]CC", "Ui": "CC=CC"}


def system5(chains, *, transfer=None, move=None, V_poly=1.0):
    """Core: N2, then per pool of LABELS5 its proxy and three mu slots (P 2-4,
    Re 6-8, Ut 10-12, Rm 14-16, Ui 18-20); ``transfer`` and ``move`` map a
    pool label to its term list."""
    core, pools = [_spc("N#N", "N2")], []
    for i, label in enumerate(LABELS5):
        stem = "C" * (i + 1)
        core += [_spc(PROXIES5[label], label), _spc(stem + "O", f"{label}_mu0"), _spc(stem + "N", f"{label}_mu1"),
                 _spc(stem + "S", f"{label}_mu2")]
        pools.append(PolymerPoolConfig(label=label, xs=2, explicit_dp_to_species_index={},
                                       mu_indices=(2 + 4 * i, 3 + 4 * i, 4 + 4 * i),
                                       kernel_transfer=(transfer or {}).get(label),
                                       kernel_move=(move or {}).get(label)))
    mask = np.array([s.label == "N2" for s in core], dtype=bool)
    rs = HybridPolymerSystem(T=700.0, P=1.0e5, initial_mole_fractions={core[0]: 1.0}, V_poly=V_poly,
                             polymer_pools=pools, mass_transfer=[], gas_species_mask=mask,
                             constant_gas_volume=False,
                             initial_polymer_moments={label: moments(chains.get(label, {})) for label in LABELS5},
                             termination=[])
    rs.initialize_model(core, [], [], [])
    return rs


def rates5(rs):
    dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
    return {label: np.array(dn[2 + 4 * i:5 + 4 * i]) for i, label in enumerate(LABELS5)}, dn


def state5(rs):
    y = np.asarray(rs.y)
    return {label: np.array(y[2 + 4 * i:5 + 4 * i]) for i, label in enumerate(LABELS5)}


def ht_term(**kw):
    """H transfer on the end-radical pool: an end radical takes an H from a
    dead chain, which becomes a mid-chain radical; the dead chain of m
    periods has 2m - 4 sites."""
    term = dict(template="h_transfer_t", k_pair=0.6, phi=1.0, partner="P", product="P", partner_product="Rm",
                sites=[-4.0, 0.0, 2.0, 0.0], periods_per_unit=PPU)
    term.update(kw)
    return term


def bb_term(**kw):
    """Backbiting on the end-radical pool: an end radical becomes a mid-chain
    radical of the same length at k*(a n + b) per chain."""
    term = dict(template="backbiting_t", k=0.5, a=0.0, b=2.0, periods_per_unit=PPU, product="Rm")
    term.update(kw)
    return term


def bilinear_sums(owner, partner, k, sites, V, same):
    """pair_sums with the rate of an ordered pair (n, m periods)
    k*(c00 + c10 n + c01 m + c11 n m)*N_n*N_m/V."""
    c00, c10, c01, c11 = sites
    out = [np.zeros(3) for _ in range(4)]
    for i, Ni in owner.items():
        for j, Nj in partner.items():
            rate = k * (c00 + c10 * i + c01 * j + c11 * i * j) * Ni * Nj / V
            for m in range(3):
                out[0][m] -= rate * (i / PPU) ** m
                out[1][m] -= rate * (j / PPU) ** m
                out[2][m] += rate * (i / PPU) ** m
                out[3][m] += rate * (j / PPU) ** m
    if same:
        out[0], out[1] = out[0] + out[1], out[0] + out[1]
    return out


def move_sums(chains, k, a, b):
    """d(moments)/dt of the source pool by explicit sums: each chain of n
    periods leaves at k*(a n + b)."""
    return -np.array([math.fsum(k * (a * n + b) * N * (n / PPU) ** j for n, N in chains.items())
                      for j in range(3)])


def test_h_transfer_moment_law_is_the_site_weighted_pair_sum():
    """A monodisperse donor pool (its mu3 closure exact): every pool's flow is
    the explicit sum over ordered pairs at the donor's sites."""
    radicals, dead = {8: 0.05, 14: 0.1, 23: 0.02}, {40: 0.5}
    V = 1.5
    rs = system5({"Re": radicals, "P": dead}, transfer={"Re": [ht_term()]}, V_poly=V)
    got, dn = rates5(rs)
    r_loss, p_loss, to_p, to_rm = bilinear_sums(radicals, dead, 0.6, (-4.0, 0.0, 2.0, 0.0), V, same=False)
    assert got["Re"] == pytest.approx(r_loss, rel=1e-12)
    assert got["P"] == pytest.approx(p_loss + to_p, rel=1e-12)
    assert got["Rm"] == pytest.approx(to_rm, rel=1e-12)
    assert np.all(got["Ut"] == 0.0) and np.all(got["Ui"] == 0.0)
    assert sum(g[1] for g in got.values()) == pytest.approx(0.0, abs=1e-12 * abs(got["Re"][1]))
    assert dn[0] == 0.0


def test_h_transfer_picks_each_donor_chain_by_its_sites():
    """A polydisperse donor pool: the chains that become mid-chain radicals
    are picked in proportion to their sites 2m - 4, so the radical pool gains
    sum_m (2m - 4) m^j N_m, not the pool's mean sites times its moments (exact
    for j = 0, 1; j = 2 needs mu3, the log-Lagrange closure)."""
    radicals, dead = {9: 0.1, 30: 0.04}, {12: 0.3, 60: 0.1, 33: 0.2}
    rs = system5({"Re": radicals, "P": dead}, transfer={"Re": [ht_term()]})
    got, _dn = rates5(rs)
    _r, _p, _to_p, to_rm = bilinear_sums(radicals, dead, 0.6, (-4.0, 0.0, 2.0, 0.0), 1.0, same=False)
    assert got["Rm"][:2] == pytest.approx(to_rm[:2], rel=1e-12)
    m = [math.fsum(N * n ** j for n, N in dead.items()) for j in range(3)]       # periods
    m3 = m[0] * (m[2] / m[1]) ** 3
    k0 = 0.6 * moments(radicals)[0]
    assert got["Rm"][2] == pytest.approx(k0 * (-4.0 * m[2] + 2.0 * m3) / PPU ** 2, rel=1e-12)
    uniform = k0 * (-4.0 * m[0] + 2.0 * m[1]) / m[0] * m[1] / PPU
    assert abs(got["Rm"][1] - uniform) > 1e-3 * abs(uniform)


def test_sites_on_both_participants_sum_the_bilinear_law():
    radicals, dead = {12: 0.2}, {30: 0.4}
    sites = [1.5, 0.25, 0.5, 0.125]
    rs = system5({"Re": radicals, "P": dead}, transfer={"Re": [ht_term(sites=sites, k_pair=0.3)]}, V_poly=0.7)
    got, _dn = rates5(rs)
    r_loss, p_loss, to_p, to_rm = bilinear_sums(radicals, dead, 0.3, sites, 0.7, same=False)
    assert got["Re"] == pytest.approx(r_loss, rel=1e-12)
    assert got["P"] == pytest.approx(p_loss + to_p, rel=1e-12)
    assert got["Rm"] == pytest.approx(to_rm, rel=1e-12)


def test_constant_sites_are_the_plain_pair_transfer():
    chains = {"Re": {9: 0.1, 20: 0.2}, "P": {15: 0.3, 44: 0.2}}
    plain = {k: v for k, v in ht_term(k_pair=1.2).items() if k not in ("sites", "periods_per_unit")}
    a = rates5(system5(chains, transfer={"Re": [plain]}))[0]
    b = rates5(system5(chains, transfer={"Re": [ht_term(k_pair=0.4, sites=[3.0, 0.0, 0.0, 0.0])]}))[0]
    for label in LABELS5:
        assert b[label] == pytest.approx(a[label], rel=1e-14, abs=0.0)


def test_mid_mid_disproportionation_pairs_the_pool_with_itself():
    """R_m + R_m -> P + U_i: the donor's length into the internal-alkene pool,
    nothing into the terminal-alkene pool."""
    mids = {7: 0.05, 18: 0.1, 40: 0.02}
    term = dict(template="d3_t", k_pair=2.0, phi=1.0, partner="Rm", product="P", partner_product="Ui")
    rs = system5({"Rm": mids}, transfer={"Rm": [term]}, V_poly=1.3)
    got, _dn = rates5(rs)
    loss, _same, to_p, to_ui = pair_sums(mids, mids, 2.0, 1.3, same=True)
    assert got["Rm"] == pytest.approx(loss, rel=1e-12)
    assert got["P"] == pytest.approx(to_p, rel=1e-12)
    assert got["Ui"] == pytest.approx(to_ui, rel=1e-12)
    assert np.all(got["Ut"] == 0.0)


def test_end_mid_reverse_disproportionation_grows_with_the_donor():
    """U_t + P -> R_e + R_m (x phi), the dead chain's sites 2m - 4."""
    alkenes, dead = {10: 0.02, 26: 0.05}, {48: 0.3}
    term = ht_term(template="reverse_em_t", k_pair=0.9, phi=0.25, product="Re", partner_product="Rm")
    rs = system5({"Ut": alkenes, "P": dead}, transfer={"Ut": [term]})
    got, _dn = rates5(rs)
    u_loss, p_loss, to_re, to_rm = bilinear_sums(alkenes, dead, 0.9 * 0.25, (-4.0, 0.0, 2.0, 0.0), 1.0, same=False)
    assert got["Ut"] == pytest.approx(u_loss, rel=1e-12)
    assert got["P"] == pytest.approx(p_loss, rel=1e-12)
    assert got["Re"] == pytest.approx(to_re, rel=1e-12)
    assert got["Rm"] == pytest.approx(to_rm, rel=1e-12)


def test_mid_mid_disproportionation_alone_follows_its_exact_solution():
    mids = {6: 0.2, 11: 0.3, 25: 0.1}
    k = 0.4
    m0 = moments(mids)
    term = dict(template="d3_t", k_pair=k, phi=1.0, partner="Rm", product="P", partner_product="Ui")
    rs = system5({"Rm": mids}, transfer={"Rm": [term]})
    alpha = 2.0 * k * m0[0]
    for t in (0.5, 2.0, 6.0):
        rs.advance(t)
        y, f = state5(rs), 1.0 / (1.0 + alpha * t)
        assert y["Rm"] == pytest.approx([m * f for m in m0], rel=1e-5)
        assert y["P"] == pytest.approx([m * (1.0 - f) / 2.0 for m in m0], rel=1e-5)
        assert y["Ui"] == pytest.approx([m * (1.0 - f) / 2.0 for m in m0], rel=1e-5)


def test_backbiting_moment_law_is_the_per_chain_sum():
    radicals = {8: 0.05, 14: 0.1, 23: 0.02}
    rs = system5({"Re": radicals}, move={"Re": [bb_term()]})
    got, dn = rates5(rs)
    loss = move_sums(radicals, 0.5, 0.0, 2.0)
    assert got["Re"] == pytest.approx(loss, rel=1e-12)
    assert got["Rm"] == pytest.approx(-loss, rel=1e-12)
    assert sum(np.abs(got[k]).sum() for k in ("P", "Ut", "Ui")) == 0.0 and dn[0] == 0.0


def test_a_length_dependent_move_weights_each_chain_by_its_rate():
    radicals = {8: 0.05, 14: 0.1, 23: 0.02}
    rs = system5({"Re": radicals}, move={"Re": [bb_term(a=0.5, b=-1.0)]})
    got, _dn = rates5(rs)
    loss = move_sums(radicals, 0.5, 0.5, -1.0)
    assert got["Re"][:2] == pytest.approx(loss[:2], rel=1e-12)
    assert got["Rm"][:2] == pytest.approx(-loss[:2], rel=1e-12)
    mono = {30: 0.2}
    got = rates5(system5({"Re": mono}, move={"Re": [bb_term(a=0.5, b=-1.0)]}))[0]
    assert got["Re"] == pytest.approx(move_sums(mono, 0.5, 0.5, -1.0), rel=1e-12)


def test_backbiting_alone_follows_its_exact_solution():
    radicals = {8: 0.05, 14: 0.1, 23: 0.02}
    m0 = moments(radicals)
    rs = system5({"Re": radicals}, move={"Re": [bb_term(k=0.3, b=2.0)]})
    for t in (0.5, 2.0, 6.0):
        rs.advance(t)
        y, f = state5(rs), math.exp(-0.6 * t)
        assert y["Re"] == pytest.approx([m * f for m in m0], rel=1e-5)
        assert y["Rm"] == pytest.approx([m * (1.0 - f) for m in m0], rel=1e-5)


def test_backbiting_and_its_reverse_reach_the_two_state_balance():
    """R_e -> R_m at kf, R_m -> R_e at kr per chain: every chain flips
    independently, R_e(t) = R_e(0) (kr + kf exp(-(kf + kr) t)) / (kf + kr)."""
    radicals = {8: 0.05, 14: 0.1, 23: 0.02}
    m0 = moments(radicals)
    kf, kr = 0.6, 0.2
    rs = system5({"Re": radicals}, move={"Re": [bb_term(k=0.3, b=2.0)],
                                         "Rm": [bb_term(template="reverse_backbiting_t", k=0.1, b=2.0,
                                                        product="Re")]})
    for t in (0.5, 2.0, 6.0, 30.0):
        rs.advance(t)
        y = state5(rs)
        f = (kr + kf * math.exp(-(kf + kr) * t)) / (kf + kr)
        assert y["Re"] == pytest.approx([m * f for m in m0], rel=1e-5)
        assert y["Rm"] == pytest.approx([m * (1.0 - f) for m in m0], rel=1e-5)


def test_backbiting_only_leaves_the_dead_chains_alone():
    chains = {"P": {60: 1.0, 70: 0.5}, "Re": {8: 0.05, 14: 0.1}}
    rs = system5(chains, move={"Re": [bb_term()], "Rm": [bb_term(template="rev", k=0.1, product="Re")]})
    p0, total0 = state5(rs)["P"], state5(rs)["Re"] + state5(rs)["Rm"]
    rs.advance(10.0)
    y = state5(rs)
    assert np.array_equal(y["P"], p0)
    assert y["Re"] + y["Rm"] == pytest.approx(total0, rel=1e-9)


def test_no_radicals_and_no_unsaturation_is_flat_with_every_kernel():
    d3 = dict(template="d3_t", k_pair=2.0, phi=1.0, partner="Rm", product="P", partner_product="Ui")
    rev_i = ht_term(template="reverse_int_t", phi=0.5, product="Rm", partner_product="Rm")
    rs = system5({"P": {60: 1.0, 70: 0.5}},
                 transfer={"Re": [ht_term()], "Rm": [d3], "Ui": [rev_i],
                           "Ut": [ht_term(template="reverse_em_t", product="Re", partner_product="Rm")]},
                 move={"Re": [bb_term()], "Rm": [bb_term(template="rev", product="Re")]})
    y0 = np.array(rs.y)
    assert np.all(rates5(rs)[1] == 0.0)
    rs.advance(10.0)
    assert np.array_equal(np.asarray(rs.y), y0)


def test_phi_zero_with_only_the_internal_reverse_is_flat():
    rev_i = ht_term(template="reverse_int_t", phi=0.0, product="Rm", partner_product="Rm")
    rs = system5({"P": {60: 1.0, 70: 0.5}, "Ui": {20: 0.3, 30: 0.1}}, transfer={"Ui": [rev_i]})
    y0 = np.array(rs.y)
    assert np.all(rates5(rs)[1] == 0.0)
    rs.advance(10.0)
    assert np.array_equal(np.asarray(rs.y), y0)


@pytest.mark.parametrize("change, match", [
    (dict(sites=[1.0, 2.0, 3.0]), "sites="),
    (dict(sites=[1.0, float("nan"), 0.0, 0.0]), "sites="),
    (dict(sites=[0.0, 0.0, 0.0, 0.0]), "sites="),
    (dict(sites="1"), "sites="),
    (dict(periods_per_unit=0), "periods_per_unit=0"),
])
def test_a_malformed_site_law_is_refused(change, match):
    with pytest.raises(ValueError, match=match):
        SP.validate_kernel_transfer("Re", [ht_term(**change)])


def test_sites_need_their_period():
    term = ht_term()
    del term["periods_per_unit"]
    with pytest.raises(ValueError, match="keys"):
        SP.validate_kernel_transfer("Re", [term])


@pytest.mark.parametrize("change, match", [
    (dict(k=0.0), "k=0.0"),
    (dict(a=0.0, b=0.0), "a.*b"),
    (dict(a=float("inf")), "a=inf"),
    (dict(product="Re"), "product"),
    (dict(periods_per_unit=0), "periods_per_unit=0"),
])
def test_a_malformed_move_term_is_refused(change, match):
    with pytest.raises(ValueError, match=match):
        SP.validate_kernel_move("Re", [bb_term(**change)])


def test_move_terms_need_exactly_the_keys():
    term = bb_term()
    del term["a"]
    with pytest.raises(ValueError, match="keys"):
        SP.validate_kernel_move("Re", [term])


def test_a_move_naming_an_unconfigured_pool_is_refused():
    with pytest.raises(ValueError, match="'Nowhere'"):
        system5({"Re": {10: 1.0}}, move={"Re": [bb_term(product="Nowhere")]})


def _deck_pool(**kw):
    from unittest.mock import MagicMock
    from rmgpy.rmg import input as rmg_input

    mock_rmg = MagicMock()
    mock_rmg.initial_species = []
    mock_rmg.reaction_model.iteration_num = 0
    mock_rmg.reaction_model.new_species_list = []
    mock_rmg.reaction_model.make_new_species.side_effect = lambda obj, **_k: (obj, True)
    old_rmg, old_sd = rmg_input.rmg, rmg_input.species_dict
    rmg_input.rmg, rmg_input.species_dict = mock_rmg, {}
    args = dict(label="Rm", monomer="[CH2][CH2]", end_groups=["[H]", "[H]"], cutoff=4, moments=[0.0, 0.0, 0.0],
                initial_mass=0.0, interior_feature="radical",
                kernel_move=[bb_term(template="rev", product="Re")])
    args.update(kw)
    try:
        return rmg_input.polymer(**args)
    finally:
        rmg_input.rmg, rmg_input.species_dict = old_rmg, old_sd


def test_the_deck_declares_interior_radical_and_alkene_pools():
    """Their proxies are the dead chain's oligomer with the radical (or a C=C)
    at its most central backbone position, so each pool is distinct from the
    dead-chain, end-radical and terminal-alkene pools of the same monomer."""
    rm = _deck_pool()
    ui = _deck_pool(label="Ui", interior_feature="alkene", kernel_move=None)
    p, r = _deck_polymers()
    u = _deck_alkene_pool()
    assert rm.interior_feature == "radical" and ui.interior_feature == "alkene"
    assert rm.molecule[0].is_isomorphic(Molecule(smiles="CC[CH]CCC"))
    assert ui.molecule[0].is_isomorphic(Molecule(smiles="CCC=CCC"))
    assert len({x.fingerprint for x in (p, r, u, rm, ui)}) == 5
    assert rm.kernel_move == SP.validate_kernel_move("Rm", [bb_term(template="rev", product="Re")])
    copied = rm.copy(deep=True)
    assert copied.kernel_move == rm.kernel_move and copied.interior_feature == "radical"
    assert copied.fingerprint == rm.fingerprint


@pytest.mark.parametrize("change, match", [
    (dict(interior_feature="middle"), "interior_feature"),
    (dict(end_radical_site="primary"), "interior_feature.*end_radical_site"),
    (dict(end_alkene_site="primary"), "interior_feature.*end_alkene_site"),
    (dict(kernel_move=[bb_term(k=-1.0)]), "k=-1.0"),
])
def test_the_deck_refuses_a_malformed_interior_pool(change, match):
    with pytest.raises(InputError, match=match):
        _deck_pool(**change)


def test_the_artifact_records_moves_and_growing_sites_at_3_4():
    from rmgpy.polymer import build_polymer_moments_artifact
    from rmgpy.tools.polymer_moments_runner import _check_schema_version_known
    p, r = _deck_polymers()
    rm = _deck_pool()
    artifact = build_polymer_moments_artifact([p, r, rm], core_species=None, core_reactions=[],
                                              configured_pool_labels=["P", "Re", "Rm"], monomer_routing_by_pool={})
    entries = {e["label"]: e for e in artifact["pools"]}
    assert entries["Rm"]["catalog_kernels"] == {"move": rm.kernel_move}
    assert artifact["schema_version"] == "3.4"
    with pytest.raises(Exception, match="3.4"):
        _check_schema_version_known(artifact)
    u = _deck_alkene_pool(kernel_transfer=[ht_term(product="Re", partner_product="Rm")])
    grown = build_polymer_moments_artifact([p, r, u], core_species=None, core_reactions=[],
                                           configured_pool_labels=["P", "Re", "Ut"], monomer_routing_by_pool={})
    assert grown["schema_version"] == "3.4"


def test_a_born_empty_interior_pool_keeps_its_move_kernel():
    from rmgpy.rmg.polymer_input import derive_daughter_pool_configs
    rm = _deck_pool()
    dummies = []
    for k in range(3):
        s = _spc("[Ne]", f"Rm_mu{k}")
        s.reactive = False
        s.is_moment_dummy = True
        dummies.append(s)
    core = [_spc("CCCC", "P"), rm] + dummies
    (derived,) = derive_daughter_pool_configs(core, {s: i for i, s in enumerate(core)}, existing_pool_labels={"P"})
    assert derived.label == "Rm"
    assert derived.kernel_move == rm.kernel_move


def test_the_pool_config_carries_the_move_kernel():
    from rmgpy.rmg.polymer_input import PolymerPool
    mu = [_spc("CO", "Re_mu0"), _spc("C=O", "Re_mu1"), _spc("C#N", "Re_mu2")]
    proxy = _spc("[CH2]CCC", "Re")
    pool = PolymerPool(label="Re", xs=4, monomer=_spc("C=C", "Re_repeat"), explicit_map={}, mu_species=mu,
                       proxy_species=proxy, kernel_move=[bb_term()])
    config = pool.to_config({s: i for i, s in enumerate([proxy] + mu)})
    assert config.kernel_move == SP.validate_kernel_move("Re", [bb_term()])
