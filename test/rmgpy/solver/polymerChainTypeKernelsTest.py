"""Red tests (I-278): what MoM needs to carry a kinetic Monte Carlo kernel
set whose chemistry has chain types with a residue, on the existing moments.

* ``kernel_homolysis`` with two daughters: optional keys ``partner_daughter``,
  ``partner_min_fragment`` (e') and ``residue`` (r), all three or none.  A
  chain of n periods has a*n + b sites, b = a*(1 - r - e - e'); a cut leaves k
  periods in ``daughter`` and n - r - k in ``partner_daughter``, k uniform over
  [e, n - r - e'], a sites per position.  The credits are exact polynomials of
  the moments (m3 from the existing closure); without the keys the legacy law
  runs untouched.
* ``kernel_interior_scission`` with ``residue``: the products' lengths sum to
  n - residue periods per event.
* ``end_class`` also keys a terminal-alkene pool and an interior pool.
* ``end_alkene_proxy`` takes the BACKBONE neighbour of the end atom (the
  branch holding the most heavy atoms), and ``end_alkene_proxy`` and
  ``interior_proxy`` honour ``end_class``.
* Chain types are end-group variants: pools of one repeat unit, different end
  groups, may exchange chains.

Every expected value is an explicit sum over the chains of a discrete
distribution (monodisperse where the log-Lagrange mu3 enters).
"""

import logging

import numpy as np
import pytest

from rmgpy.exceptions import InputError
import rmgpy.solver.polymer as SP

from polymerKernelsTest import PPU, _spc, moments, t_term
from polymerEndClassKernelsTest import _artifact, _deck, _pool, _slots, _system
from polymerInteriorScissionTest import s_term

W, PHI, A_SITES = 0.3, 0.5, 1.5


def h2_term(e=2, e2=3, residue=1, a=A_SITES, **kw):
    """A two-daughter homolysis term: sites a*n + a*(1 - r - e - e')."""
    term = dict(template="homolysis_two", w=W, phi=PHI, a=a, b=a * (1 - residue - e - e2), min_fragment=e,
                periods_per_unit=PPU, daughter="D", partner_daughter="Q", partner_min_fragment=e2, residue=residue)
    term.update(kw)
    return term


def brute_two(chains, e, e2, residue, a=A_SITES):
    """d(moments)/dt (pool units, mol/s) of the parent, daughter and partner daughter by explicit sums over every
    chain and cut position k in [e, n - r - e'], a sites per position."""
    par, dau, prt = np.zeros(3), np.zeros(3), np.zeros(3)
    for n, N in chains.items():
        for k in range(e, n - residue - e2 + 1):
            rate = W * PHI * a * N
            for j in range(3):
                par[j] -= rate * (n / PPU) ** j
                dau[j] += rate * (k / PPU) ** j
                prt[j] += rate * ((n - residue - k) / PPU) ** j
    return par, dau, prt


def run_h(term, chains):
    pools = [_pool("P", 0, kernel_homolysis=[term]), _pool("D", 1), _pool("Q", 2)]
    rs = _system(pools, {"P": chains})
    dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
    return _slots(dn, 0), _slots(dn, 1), _slots(dn, 2)


# --------------------------------------------------------------------------
# kernel_homolysis: two daughters
# --------------------------------------------------------------------------
@pytest.mark.parametrize("e, e2, residue", [(2, 3, 1), (3, 2, 1), (2, 3, 0), (4, 1, 2), (2, 2, 1)])
@pytest.mark.parametrize("n", [9, 14, 60])
def test_two_daughters_credit_the_cut_sums_on_a_monodisperse_pool(e, e2, residue, n):
    chains = {n: 0.7}                    # the log-Lagrange mu3 is exact here
    got = run_h(h2_term(e, e2, residue), chains)
    for g, want in zip(got, brute_two(chains, e, e2, residue)):
        assert g == pytest.approx(want, rel=1e-12)


@pytest.mark.parametrize("e, e2, residue", [(2, 3, 1), (3, 2, 0)])
def test_two_daughters_credit_events_and_units_on_a_spread_pool(e, e2, residue):
    """mu3 enters only the third-moment rows: the event count, the parent's lengths and the daughters' mu1 are
    exact sums on any distribution."""
    chains = {n: 0.1 * (1 + n % 3) for n in range(6, 13)}
    got = run_h(h2_term(e, e2, residue), chains)
    want = brute_two(chains, e, e2, residue)
    for g, w in zip(got, want):
        assert g[0] == pytest.approx(w[0], rel=1e-12)
        assert g[1] == pytest.approx(w[1], rel=1e-12)
    assert got[0][2] < 0.0 and got[1][2] > 0.0 and got[2][2] > 0.0


def test_two_daughters_conserve_events_and_all_but_the_residue_in_periods():
    chains = {20: 0.7}
    par, dau, prt = run_h(h2_term(2, 3, 1), chains)
    assert dau[0] == pytest.approx(-par[0], rel=1e-13) and prt[0] == pytest.approx(-par[0], rel=1e-13)
    # per event the fragments hold n - residue periods: the lost period is 1/PPU of a unit per event
    assert par[1] + dau[1] + prt[1] == pytest.approx(par[0] * 1.0 / PPU, rel=1e-12)


def test_two_daughters_may_share_one_pool():
    term = h2_term(2, 3, 1, partner_daughter="D")
    pools = [_pool("P", 0, kernel_homolysis=[term]), _pool("D", 1)]
    rs = _system(pools, {"P": {20: 0.7}})
    dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
    par, dau, prt = brute_two({20: 0.7}, 2, 3, 1)
    assert _slots(dn, 1) == pytest.approx(dau + prt, rel=1e-12)


def test_a_term_without_the_optional_keys_credits_as_it_did():
    from polymerKernelsTest import homolysis_sums
    legacy = dict(template="homolysis_t", w=W, phi=PHI, a=1.0, b=-7.0, min_fragment=4, periods_per_unit=PPU,
                  daughter="D")
    pools = [_pool("P", 0, kernel_homolysis=[legacy]), _pool("D", 1)]
    rs = _system(pools, {"P": {60: 1.5}})
    dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
    par, dau = homolysis_sums({60: 1.5}, W, PHI)
    assert _slots(dn, 0) == pytest.approx(par, rel=1e-12) and _slots(dn, 1) == pytest.approx(dau, rel=1e-12)


def test_a_two_daughter_term_outside_its_law_domain_is_skipped_loudly(caplog):
    """Chains of 5 periods hold fewer than e + e' + r: negative sites."""
    with caplog.at_level(logging.WARNING):
        got = run_h(h2_term(2, 3, 1), {5: 0.7})
    assert all(np.all(g == 0.0) for g in got)
    assert "kernel_homolysis" in caplog.text


# --------------------------------------------------------------------------
# validation
# --------------------------------------------------------------------------
def test_a_valid_two_daughter_term_is_accepted_and_normalized():
    out = SP.validate_kernel_homolysis("P", [h2_term()])
    assert out[0]["partner_daughter"] == "Q" and out[0]["partner_min_fragment"] == 3 and out[0]["residue"] == 1
    assert SP.validate_kernel_homolysis("P", [h2_term(residue=0)])[0]["residue"] == 0


@pytest.mark.parametrize("term, match", [
    (h2_term(b=-5.0), r"b = a\*\(1 - residue - min_fragment - partner_min_fragment\)"),
    (h2_term(partner_min_fragment=0), "partner_min_fragment=0 must be an integer >= 1"),
    (h2_term(residue=-1), "residue=-1 must be an integer >= 0"),
    (h2_term(residue=1.0), "residue=1.0 must be an integer >= 0"),
    (h2_term(partner_daughter="P"), "partner_daughter='P' must name another pool"),
    (h2_term(valid_from_periods=3), "valid_from_periods=3 is below"),
])
def test_a_malformed_two_daughter_term_is_refused(term, match):
    with pytest.raises(ValueError, match=match):
        SP.validate_kernel_homolysis("P", [term])


@pytest.mark.parametrize("drop", ["partner_daughter", "partner_min_fragment", "residue"])
def test_the_optional_keys_come_all_together_or_not_at_all(drop):
    term = h2_term()
    del term[drop]
    with pytest.raises(ValueError, match=r"or with all of \['partner_daughter', 'partner_min_fragment', 'residue'\]"):
        SP.validate_kernel_homolysis("P", [term])


def test_a_terms_validity_floor_is_e_plus_e_prime_plus_residue_minus_one():
    pools = [_pool("P", 0, kernel_homolysis=[h2_term(2, 3, 1)]), _pool("D", 1), _pool("Q", 2)]
    rs = _system(pools, {"P": {20: 0.7}})
    entry = [v for v in rs.kernel_validity if v["kernel"] == "kernel_homolysis"][0]
    assert entry["valid_from_periods"] == 5 and entry["source"] == "derived"
    assert SP.validate_kernel_homolysis("P", [h2_term(2, 3, 1, valid_from_periods=5)])[0]["valid_from_periods"] == 5


def test_the_partner_daughter_is_connected_and_must_be_configured():
    pool = _pool("P", 0, kernel_homolysis=[h2_term()])
    links = SP.catalog_kernel_connections(pool)
    assert ("kernel_homolysis", "homolysis_two", "P", "D") in links
    assert ("kernel_homolysis", "homolysis_two", "P", "Q") in links
    pools = [_pool("P", 0, kernel_homolysis=[h2_term(partner_daughter="Nowhere")]), _pool("D", 1)]
    with pytest.raises(ValueError, match="partner_daughter pool 'Nowhere', which is not among the configured"):
        _system(pools, {"P": {20: 0.7}})


def test_a_legacy_term_has_only_its_one_connection():
    legacy = dict(template="homolysis_t", w=W, phi=PHI, a=1.0, b=-7.0, min_fragment=4, periods_per_unit=PPU,
                  daughter="D")
    links = SP.catalog_kernel_connections(_pool("P", 0, kernel_homolysis=[legacy]))
    assert [l for l in links if l[0] == "kernel_homolysis"] == [("kernel_homolysis", "homolysis_t", "P", "D")]


def test_a_two_daughter_term_carries_through_the_deck_and_the_pool_config():
    from rmgpy.rmg.polymer_input import PolymerPool
    r = _deck(label="P", end_radical_site=None, kernel_homolysis=[h2_term()])
    assert r.kernel_homolysis == SP.validate_kernel_homolysis("P", [h2_term()])
    assert r.copy(deep=True).kernel_homolysis == r.kernel_homolysis
    mus = [_spc("CO", "P_mu0"), _spc("C=O", "P_mu1"), _spc("C#N", "P_mu2")]
    pool = PolymerPool(label="P", xs=2, monomer=_spc("C=C", "P_repeat"), explicit_map={}, mu_species=mus,
                       kernel_homolysis=[h2_term()])
    config = pool.to_config({mus[0]: 1, mus[1]: 2, mus[2]: 3})
    assert config.kernel_homolysis == SP.validate_kernel_homolysis("P", [h2_term()])
    with pytest.raises(InputError, match="partner_min_fragment=0"):
        _deck(label="P", end_radical_site=None, kernel_homolysis=[h2_term(partner_min_fragment=0)])


# --------------------------------------------------------------------------
# kernel_interior_scission: residue
# --------------------------------------------------------------------------
def _residue_term(residue):
    """Two products whose lengths sum to n - residue: A gets n - residue - 5 periods, R gets 5."""
    return dict(s_term(), residue=residue,
                products=[{"pool": "A", "moments": [[1.0], [-5.0 - residue, 1.0], [(5.0 + residue) ** 2, -2 * (5.0 + residue), 1.0]]},
                          {"pool": "R", "moments": [[1.0], [5.0, 0.0], [25.0, 0.0, 0.0]]}])


def test_an_interior_scission_with_a_residue_conserves_all_but_the_residue():
    out = SP.validate_kernel_interior_scission("P", [_residue_term(1)])
    assert out[0]["residue"] == 1
    assert SP.validate_kernel_interior_scission("P", [s_term()])[0].get("residue") is None


def test_an_interior_scission_with_a_residue_credits_n_minus_residue_periods():
    pools = [_pool("P", 0, kernel_interior_scission=[_residue_term(1)]), _pool("A", 1), _pool("R", 2)]
    rs = _system(pools, {"P": {24: 0.3}})
    dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
    parent, a, r = _slots(dn, 0), _slots(dn, 1), _slots(dn, 2)
    assert parent[1] + a[1] + r[1] == pytest.approx(parent[0] / PPU, rel=1e-12)   # one period lost per event


@pytest.mark.parametrize("term, match", [
    (dict(_residue_term(1), residue=0), r"residue 0"),
    (dict(_residue_term(1), residue=2), r"residue 2"),
    (dict(_residue_term(1), residue=-1), "residue=-1 must be an integer >= 0"),
    (dict(_residue_term(0)), r"residue 0"),
    (dict(s_term(), residue=1), r"residue 1"),
])
def test_an_interior_scission_whose_lengths_do_not_match_its_residue_is_refused(term, match):
    with pytest.raises(ValueError, match=match):
        SP.validate_kernel_interior_scission("P", [term])


# --------------------------------------------------------------------------
# end_class on terminal-alkene and interior pools; the backbone neighbour
# --------------------------------------------------------------------------
PP = "[CH](C)[CH2]"           # head port *2 on the CH2, tail port *1 on the CH


def _alkene(label="Ut", **kw):
    args = dict(label=label, monomer=PP, end_radical_site=None, end_alkene_site="primary")
    args.update(kw)
    return _deck(**args)


def _interior(label="Rm", **kw):
    args = dict(label=label, monomer=PP, end_radical_site=None, interior_feature="radical")
    args.update(kw)
    return _deck(**args)


def _label(atom):
    h = sum(1 for b in atom.edges if b.is_hydrogen())
    return atom.symbol + ("" if h == 0 else "H" if h == 1 else f"H{h}")


def _alkene_labels(spc):
    mol = spc.molecule[0]
    bonds = [(a, b) for a in mol.atoms for b in a.edges if a is not b and not a.is_hydrogen() and not b.is_hydrogen()
             and mol.get_bond(a, b).is_double() and mol.atoms.index(a) < mol.atoms.index(b)]
    assert len(bonds) == 1
    return {_label(bonds[0][0]), _label(bonds[0][1])}, bonds[0]


def test_an_end_class_keys_alkene_and_interior_pools():
    plain, keyed = _alkene(), _alkene(label="Ut1", end_class="CHCH")
    assert keyed.end_class == "CHCH" and plain.end_class is None
    assert keyed.fingerprint != plain.fingerprint and "-class" not in plain.fingerprint
    assert "_EndAlkene-primary-classCHCH" in keyed.fingerprint
    assert keyed.copy(deep=True).fingerprint == keyed.fingerprint and keyed.copy(deep=True).end_class == "CHCH"
    ip, ik = _interior(position_class="d3"), _interior(label="Rm1", position_class="d3", end_class="C")
    assert ik.fingerprint != ip.fingerprint and "-class" not in ip.fingerprint
    assert ik.fingerprint.count("-classC") == 1 and ik.fingerprint.index("-pos") < ik.fingerprint.index("-classC")
    assert ik.copy(deep=True).fingerprint == ik.fingerprint


@pytest.mark.parametrize("kw", [dict(end_class="C H"), dict(end_class="")])
def test_a_malformed_end_class_on_an_alkene_or_interior_pool_is_refused(kw):
    with pytest.raises(InputError, match="end_class"):
        _alkene(**kw)
    with pytest.raises(InputError, match="end_class"):
        _interior(**kw)


def test_an_end_class_still_needs_a_site_or_a_feature():
    with pytest.raises(InputError, match="end_class"):
        _deck(end_radical_site=None, end_class="CH3")


@pytest.mark.parametrize("site, labels", [("primary", {"CH"}), ("secondary", {"CH2", "C"})])
def test_a_polypropylene_end_alkene_takes_its_backbone_neighbour(site, labels):
    pool = _alkene(end_alkene_site=site)
    got, (a, b) = _alkene_labels(pool.end_alkene_proxy)
    assert got == labels
    mol = pool.end_alkene_proxy.molecule[0]
    # the pi bond is to the backbone neighbour: the methyl side group keeps its three hydrogens and no double bond
    methyls = [x for x in mol.atoms if x.symbol == "C" and _label(x) == "CH3"]
    assert all(not any(mol.get_bond(x, y).is_double() for y in x.edges) for x in methyls)


def test_the_polyethylene_end_alkene_is_what_it_was():
    from rmgpy.molecule import Molecule
    pe = _alkene(monomer="[CH2][CH2]")
    assert pe.molecule[0].is_isomorphic(Molecule(smiles="C=CCCCC"))
    assert pe.end_alkene_proxy.molecule[0].is_isomorphic(Molecule(smiles="C=CCCCC"))


@pytest.mark.parametrize("site, end_class", [("primary", "CHCH"), ("secondary", "CH2C")])
def test_an_alkene_end_class_that_matches_the_formed_feature_is_accepted(site, end_class):
    pool = _alkene(end_alkene_site=site, end_class=end_class)
    assert pool.end_alkene_proxy is not None


@pytest.mark.parametrize("site, end_class", [("primary", "CH2C"), ("secondary", "CHCH")])
def test_an_alkene_end_class_that_does_not_match_is_refused_naming_both(site, end_class):
    pool = _alkene(end_alkene_site=site, end_class=end_class)
    formed = "CHCH" if site == "primary" else "CH2C"
    with pytest.raises(ValueError, match=f"{end_class}.*{formed}|{formed}.*{end_class}"):
        pool.end_alkene_proxy


def test_an_interior_radical_end_class_picks_the_matching_backbone_atom():
    tertiary = _interior(end_class="C").interior_proxy
    secondary = _interior(label="Rm1", end_class="CH").interior_proxy
    for spc, label in ((tertiary, "C"), (secondary, "CH")):
        radicals = [a for a in spc.molecule[0].atoms if a.radical_electrons > 0]
        assert len(radicals) == 1 and _label(radicals[0]) == label
    assert not tertiary.molecule[0].is_isomorphic(secondary.molecule[0])
    # no end_class: the pool is what it was (the central tertiary atom)
    assert _interior(label="Rm2").interior_proxy.molecule[0].is_isomorphic(tertiary.molecule[0])


def test_an_interior_alkene_end_class_picks_the_matching_pair_in_either_order():
    plain = _interior(label="Ia", interior_feature="alkene").interior_proxy
    assert _alkene_labels(plain)[0] == {"C", "CH"}
    for end_class in ("CCH", "CHC"):
        got = _interior(label="Ib", interior_feature="alkene", end_class=end_class).interior_proxy
        assert got.molecule[0].is_isomorphic(plain.molecule[0])
    other = _interior(label="Ic", interior_feature="alkene", end_class="CHCH").interior_proxy
    assert _alkene_labels(other)[0] == {"CH"} and not other.molecule[0].is_isomorphic(plain.molecule[0])


def test_an_interior_end_class_no_atom_matches_is_refused():
    with pytest.raises(ValueError, match="CH3"):
        _interior(end_class="CH3").interior_proxy
    with pytest.raises(ValueError, match="CH3"):
        _interior(interior_feature="alkene", end_class="CH3").interior_proxy


# --------------------------------------------------------------------------
# chain types are end-group variants
# --------------------------------------------------------------------------
@pytest.mark.parametrize("end_groups, site", [(["[H]", "[CH3]"], "secondary"), (["C[CH2]", "[H]"], "primary")])
def test_a_chain_type_pool_is_an_end_group_variant_of_the_declared_pool(end_groups, site):
    plain = _deck(label="Pd", monomer=PP, end_radical_site=None)
    typed = _deck(label="Pt", monomer=PP, end_radical_site=None, end_groups=end_groups)
    assert typed.fingerprint != plain.fingerprint
    rad = _deck(label="Rt", monomer=PP, end_groups=end_groups, end_radical_site=site)
    assert rad.end_radical_proxy is not None and rad.end_radical_proxy.molecule[0].multiplicity == 2
    assert rad.fingerprint != _deck(label="Rd", monomer=PP, end_radical_site=site).fingerprint


def test_chains_move_between_a_chain_type_pool_and_the_plain_pools():
    """Same repeat unit, different end groups: the repeat-unit check passes for a transfer and a homolysis."""
    from rmgpy.rmg.polymer_input import _element_count_tuple
    plain = _deck(label="Pd", monomer=PP, end_radical_site=None)
    typed = _deck(label="Pt", monomer=PP, end_radical_site=None, end_groups=["[H]", "[CH3]"])
    unit = dict(repeat_unit_atoms=_element_count_tuple(plain.monomer),
                monomer_mw_g_mol=plain.monomer.get_molecular_weight() * 1000.0)
    assert unit == dict(repeat_unit_atoms=_element_count_tuple(typed.monomer),
                        monomer_mw_g_mol=typed.monomer.get_molecular_weight() * 1000.0)
    transfer = t_term(partner="Pt", product="D", partner_product="Q")
    pools = [_pool("P", 0, kernel_transfer=[transfer], **unit), _pool("Pt", 1, **unit),
             _pool("D", 2, **unit), _pool("Q", 3, **unit)]
    _system(pools, {"P": {24: 0.3}, "Pt": {20: 0.2}})
    pools = [_pool("P", 0, kernel_homolysis=[h2_term()], **unit), _pool("D", 1, **unit), _pool("Q", 2, **unit)]
    _system(pools, {"P": {24: 0.3}})


def test_the_artifact_records_a_two_daughter_homolysis_in_its_homolysis_terms():
    p = _deck(label="P", end_radical_site=None, kernel_homolysis=[h2_term()])
    artifact = _artifact([p])
    assert artifact["pools"][0]["catalog_kernels"] == {"homolysis": p.kernel_homolysis}
