"""Red tests (I-273): what MoM needs to carry mid-chain radicals by where
they sit, on the existing moments.

* ``position_class`` keys an interior pool (``interior_feature``) by the
  position class of its feature (a fixed end distance, or the bulk), so two
  mid-chain-radical pools of one feature are two pools.
* ``kernel_interior_scission``: every chain of the pool breaks at k per
  chain [1/s] into fragments; the pool loses dmu_j = -k mu_j, and each
  product pool gains dm_j = k sum_{i<=j} q_ji m_i in the term's period
  (m_j = c^j mu_j, c = periods_per_unit), the q fixed coefficients of
  E[sum of the fragment lengths^j it receives | parent length n] = sum_i
  q_ji n^i.  That is exact whenever the radical's position is fixed or
  uniform over a range [p, n - q] (no closure); the products must conserve
  the parent's periods.
* The pool table records both at schema 3.6; without them nothing changes.

Every expected value is an explicit sum over chains of explicit (length,
radical position), the fragments of each event enumerated.  The q of each
term are obtained by exact interpolation from the same enumeration at three
lengths, never from a closed form.
"""

import logging
from fractions import Fraction

import numpy as np
import pytest

from rmgpy.exceptions import InputError
import rmgpy.solver.polymer as SP

from polymerKernelsTest import PPU, _spc, moments
from polymerEndClassKernelsTest import V, _artifact, _deck, _pool, _slots, _system


def outcomes(n, j):
    """The events of a mid-chain radical at backbone position j on n periods: [(weight, [(pool, length)])]: the
    beta bond away from the near end (terminal alkene j + 2, end radical n - j - 2) and the other one (alkene
    n - j + 1, radical j - 1), as the far law of a beta scission states them."""
    return [(1, [("A", j + 2), ("R", n - j - 2)]), (1, [("A", n - j + 1), ("R", j - 1)])]


def near_outcomes(n, j):
    """Only the bond away from the near end (a radical near an end keeps one feature)."""
    return outcomes(n, j)[:1]


def per_event(law, n, positions, pool, power):
    """E[sum of the lengths^power of the fragments ``pool`` receives] per event on a parent of n periods with
    its radical uniform over ``positions``."""
    total, events = Fraction(0), Fraction(0)
    for j in positions:
        for w, frags in law(n, j):
            events += w
            total += w * sum(Fraction(L) ** power for p, L in frags if p == pool)
    return total / events


def q_matrix(law, positions_of, pool, lengths=(30, 31, 32)):
    """[[q00], [q10, q11], [q20, q21, q22]] by exact interpolation of per_event in n on three lengths."""
    out = []
    for power in range(3):
        rows = [[Fraction(n) ** i for i in range(power + 1)] for n in lengths[:power + 1]]
        rhs = [per_event(law, n, positions_of(n), pool, power) for n in lengths[:power + 1]]
        # solve the (power + 1) x (power + 1) Vandermonde system exactly
        m = [r + [b] for r, b in zip(rows, rhs)]
        size = power + 1
        for c in range(size):
            piv = next(r for r in range(c, size) if m[r][c] != 0)
            m[c], m[piv] = m[piv], m[c]
            m[c] = [x / m[c][c] for x in m[c]]
            for r in range(size):
                if r != c:
                    m[r] = [a - m[r][c] * b for a, b in zip(m[r], m[c])]
        out.append([float(m[r][-1]) for r in range(size)])
    return out


def s_term(law=outcomes, positions_of=lambda n: [5], k=7.0, **kw):
    term = dict(template="scission_t", k=k, periods_per_unit=PPU,
                products=[{"pool": p, "moments": q_matrix(law, positions_of, p)} for p in ("A", "R")])
    term.update(kw)
    return term


def brute(chains, law, k):
    """d(moments)/dt (pool units, mol/s) of the parent, A and R by explicit sums over {(n, j): mol}."""
    dp, da, dr = np.zeros(3), np.zeros(3), np.zeros(3)
    for (n, j), N in chains.items():
        events = law(n, j)
        W = sum(w for w, _f in events)
        for m in range(3):
            dp[m] -= k * N * (n / PPU) ** m
        for w, frags in events:
            for pool, L in frags:
                for m in range(3):
                    (da if pool == "A" else dr)[m] += k * N * (w / W) * (L / PPU) ** m
    return dp, da, dr


def by_length(chains):
    out = {}
    for (n, _j), N in chains.items():
        out[n] = out.get(n, 0.0) + N
    return out


def run(term, chains):
    pools = [_pool("P", 0, kernel_interior_scission=[term]), _pool("A", 1), _pool("R", 2)]
    rs = _system(pools, {"P": by_length(chains)})
    dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
    return _slots(dn, 0), _slots(dn, 1), _slots(dn, 2)


# --------------------------------------------------------------------------
# the moment law against explicit chains
# --------------------------------------------------------------------------
@pytest.mark.parametrize("d", [3, 5, 9])
def test_a_fixed_position_pool_breaks_as_its_chains_do(d):
    chains = {(n, d): N for n, N in ((24, 0.3), (31, 0.5), (60, 0.2), (75, 0.05))}
    got = run(s_term(positions_of=lambda n, d=d: [d]), chains)
    for g, want in zip(got, brute(chains, outcomes, 7.0)):
        assert g == pytest.approx(want, rel=1e-12)


@pytest.mark.parametrize("p, q", [(6, 6), (4, 7)])
def test_a_bulk_pool_breaks_as_its_chains_do_with_the_radical_uniform_over_its_range(p, q):
    chains = {}
    for n, N in ((24, 0.3), (31, 0.5), (60, 0.2), (75, 0.05)):
        js = range(p, n - q + 1)
        for j in js:
            chains[(n, j)] = N / len(js)
    got = run(s_term(positions_of=lambda n, p=p, q=q: range(p, n - q + 1)), chains)
    for g, want in zip(got, brute(chains, outcomes, 7.0)):
        assert g == pytest.approx(want, rel=1e-12)


def test_a_pool_with_one_feature_left_near_the_end_breaks_by_that_one():
    chains = {(n, 3): N for n, N in ((24, 0.3), (60, 0.2))}
    got = run(s_term(law=near_outcomes, positions_of=lambda n: [3], k=3.0), chains)
    for g, want in zip(got, brute(chains, near_outcomes, 3.0)):
        assert g == pytest.approx(want, rel=1e-12)


def test_each_event_conserves_periods_and_makes_one_more_chain():
    chains = {(n, 5): N for n, N in ((24, 0.3), (31, 0.5), (60, 0.2))}
    parent, a, r = run(s_term(), chains)
    assert parent[0] < 0.0
    assert a[0] + r[0] == pytest.approx(-2.0 * parent[0], rel=1e-13)          # two fragments per event
    assert parent[1] + a[1] + r[1] == pytest.approx(0.0, abs=1e-13 * a[1])    # units (mass)


def test_a_term_whose_fragments_leave_its_law_domain_is_skipped_loudly(caplog):
    """A fragment law n - 40 on chains of 24 periods: negative lengths; no flux, and a warning."""
    term = s_term(products=[{"pool": "A", "moments": [[1.0], [-40.0, 1.0], [1600.0, -80.0, 1.0]]},
                            {"pool": "R", "moments": [[1.0], [40.0, 0.0], [1600.0, 0.0, 0.0]]}])
    with caplog.at_level(logging.WARNING):
        got = run(term, {(24, 5): 0.3})
    assert all(np.all(g == 0.0) for g in got)
    assert "kernel_interior_scission" in caplog.text


# --------------------------------------------------------------------------
# validation
# --------------------------------------------------------------------------
def _bad(**kw):
    term = s_term()
    term.update(kw)
    return term


@pytest.mark.parametrize("term, match", [
    (_bad(k=0.0), "k=0.0 must be > 0"),
    (_bad(products=[]), "products"),
    (_bad(products=[{"pool": "P", "moments": [[1.0], [0.0, 1.0], [0.0, 0.0, 1.0]]}]), "another pool"),
    (_bad(products=[{"pool": "A", "moments": [[1.0], [0.0, 1.0]]}]), "moments"),
    (_bad(products=[{"pool": "A", "moments": [[2.0], [0.0, 2.0], [0.0, 0.0, 2.0]]}]), "periods"),
    (dict(s_term(), extra=1), "keys"),
])
def test_a_malformed_interior_scission_is_refused(term, match):
    with pytest.raises(ValueError, match=match):
        SP.validate_kernel_interior_scission("P", [term])


def test_a_term_naming_an_unconfigured_product_is_refused():
    term = s_term()
    term["products"][1]["pool"] = "Nowhere"
    pools = [_pool("P", 0, kernel_interior_scission=[term]), _pool("A", 1)]
    with pytest.raises(ValueError, match="'Nowhere', which is not among the configured"):
        _system(pools, {"P": {24: 0.3}})


def test_the_term_connects_the_pool_to_each_product():
    pool = _pool("P", 0, kernel_interior_scission=[s_term()])
    links = SP.catalog_kernel_connections(pool)
    assert ("kernel_interior_scission", "scission_t", "P", "A") in links
    assert ("kernel_interior_scission", "scission_t", "P", "R") in links


# --------------------------------------------------------------------------
# the deck surface and the artifact
# --------------------------------------------------------------------------
def _interior(**kw):
    args = dict(label="Rm", end_radical_site=None, interior_feature="radical")
    args.update(kw)
    return _deck(**args)


def test_a_position_class_keys_a_second_interior_pool():
    bulk, near = _interior(label="Rm", position_class="bulk"), _interior(label="Rmd3", position_class="d3")
    plain = _interior(label="Rmx")
    assert near.position_class == "d3" and plain.position_class is None
    assert len({bulk.fingerprint, near.fingerprint, plain.fingerprint}) == 3
    assert near.copy(deep=True).fingerprint == near.fingerprint


@pytest.mark.parametrize("kw, match", [(dict(interior_feature=None, position_class="d3"), "position_class"),
                                       (dict(position_class="d 3"), "position_class"),
                                       (dict(position_class=""), "position_class")])
def test_a_malformed_position_class_is_refused(kw, match):
    with pytest.raises(InputError, match=match):
        _interior(**kw)


def test_an_interior_pool_carries_its_scission_through_the_deck_and_the_pool_config():
    from rmgpy.rmg.polymer_input import PolymerPool
    r = _interior(position_class="d3", kernel_interior_scission=[s_term()])
    assert r.kernel_interior_scission == SP.validate_kernel_interior_scission("Rm", [s_term()])
    assert r.copy(deep=True).kernel_interior_scission == r.kernel_interior_scission
    mus = [_spc("CO", "Rm_mu0"), _spc("C=O", "Rm_mu1"), _spc("C#N", "Rm_mu2")]
    pool = PolymerPool(label="Rm", xs=2, monomer=_spc("C=C", "Rm_repeat"), explicit_map={}, mu_species=mus,
                       kernel_interior_scission=[s_term()])
    config = pool.to_config({mus[0]: 1, mus[1]: 2, mus[2]: 3})
    assert config.kernel_interior_scission == SP.validate_kernel_interior_scission("Rm", [s_term()])


def test_the_artifact_records_an_interior_scission_at_3_6():
    r = _interior(position_class="d3", kernel_interior_scission=[s_term()])
    artifact = _artifact([r])
    assert artifact["pools"][0]["catalog_kernels"] == {"interior_scission": r.kernel_interior_scission}
    assert artifact["schema_version"] == "3.6"


def test_a_position_class_alone_leaves_the_artifact_as_it_was():
    plain, keyed = _artifact([_interior()]), _artifact([_interior(position_class="d3")])
    assert keyed["schema_version"] == plain["schema_version"]
    assert keyed["pools"][0].keys() == plain["pools"][0].keys()
