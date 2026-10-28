"""
An empty pool's moments are integrator noise, and the accepted-state
negative-moment check must judge them at the scale that noise has at each
moment order.

On PP, a born-empty pair of mid-chain-radical pools exchanges chains by fast
moves (k ~ 1e7 s^-1) and is fed far below atol.  The accepted floor of a
born-empty pool is the absolute 100*atol alone, the same for every moment
order -- ONE error weight of its moment slots -- so the check refused an
accepted state of mu0 = -9.3e-17, mu2 = -1.27e-14 (CH pool) and -5.64e-14 mol
(C pool) against a 1e-14 floor.  The noise is chains passed through the pool
at mu0's error scale, so mu_k noise ~ mean^k * mu0 noise (mu2/mu0 = 11.7^2
here, 28.4^2 on a position pool later), with the deck's chains 30 long.  Once
negative, the pair's draining terms are skipped (a flow outside its law's
domain), so the excursion is held, not pulled back.

The refusal floor is ``max(F_k, L^k * A_0)``, F the accepted floors, A_0 the
absolute mu0 floor (an empty pool's noise scale) and L the largest
weight-average length mu2/mu1 of a trusted pool.  The floors that end
a pool (pool_moment_accepted_floors) are unchanged, so an empty pool between
the two is carried; a negative excursion beyond the refusal floor is still a
hard stop.
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

DP = 30.0
GAS = ["C#C", "O", "N", "S", "C=O", "CO", "CC", "OO", "NN", "C#N", "[He]", "[Ne]", "[Ar]",
       "CF", "CCl", "FF", "ClCl", "O=O", "C=C", "CN"]
LABELS = ("P", "Re", "Rm", "Ra", "Ut")
PROXIES = {"P": "CCCC", "Re": "[CH2]CCC", "Rm": "C[CH]CC", "Ra": "CC[CH]CC", "Ut": "C=CCC"}
# beta scission of an interior radical into an end radical and a terminal
# alkene, PP's fixed coefficients
CUT = [{"pool": "Re", "moments": [[1.0], [-1.0, 0.5], [23.0 / 3.0, -2.5, 1.0 / 3.0]]},
       {"pool": "Ut", "moments": [[1.0], [1.0, 0.5], [23.0 / 3.0, -0.5, 1.0 / 3.0]]}]
# the refused accepted state on PP (fixture table, atol 1e-16, t = 2.69e-5 s)
PP_RM = (-9.270516014958691e-17, -1.513098406075246e-15, -1.2712029181068988e-14)
PP_RA = (-4.042e-16, -6.573e-15, -5.640e-14)


def _spc(smiles, label):
    s = Species(molecule=[Molecule(smiles=smiles)])
    s.label = label
    return s


def pair_deck(n_gas=20, cls=HybridPolymerSystem):
    """A dead-chain pool P cut by homolysis into a born-empty end-radical pool
    Re; Re takes an H from P, making a born-empty mid-chain radical pool Rm;
    Rm and Ra exchange chains by fast moves and Ra breaks by beta scission
    (the PP c2~c2 pair, its rates).  The pair holds ~1e-21 mol at 1e-4 s,
    five decades below atol.  ``n_gas`` inert gas species set the norm's N."""
    core = [_spc("N#N", "N2")]
    core += [_spc(s, f"gas{i}") for i, s in enumerate(GAS[:n_gas])]
    mu = {}
    for j, label in enumerate(LABELS):
        base = len(core)
        core += [_spc(PROXIES[label], label)]
        core += [_spc("C" * (j + 1) + s, f"{label}_mu{k}") for k, s in enumerate(("CO", "C=O", "C#N"))]
        mu[label] = (base + 1, base + 2, base + 3)
    move = lambda k, owner, product: [dict(template=f"move_{owner}", k=k, a=0.0, b=1.0, periods_per_unit=1,
                                           product=product)]
    kw = {"P": dict(kernel_homolysis=[dict(template="homolysis", w=1.65e-8, phi=1.0, a=1.0, b=-3.0,
                                           min_fragment=2, periods_per_unit=1, daughter="Re")]),
          "Re": dict(kernel_transfer=[dict(template="h_transfer", k_pair=1.0e-4, phi=1.0, partner="P",
                                           product="P", partner_product="Rm")]),
          "Rm": dict(kernel_move=move(3.7e7, "Rm", "Ra")),
          "Ra": dict(kernel_move=move(1.2e7, "Ra", "Rm"),
                     kernel_interior_scission=[dict(template="beta", k=3.3e6, periods_per_unit=1, products=CUT)]),
          "Ut": {}}
    pools = [PolymerPoolConfig(label=label, xs=2, explicit_dp_to_species_index={}, mu_indices=mu[label],
                               **kw[label]) for label in LABELS]
    moments = {label: (0.0, 0.0, 0.0) for label in LABELS}
    moments["P"] = (1.0, DP, DP * DP + DP)
    mask = np.array([s.label == "N2" or s.label.startswith("gas") for s in core], dtype=bool)
    rs = cls(T=700.0, P=1.0e5, initial_mole_fractions={core[0]: 1.0}, V_poly=1.0, polymer_pools=pools,
             mass_transfer=[], gas_species_mask=mask, constant_gas_volume=False,
             initial_polymer_moments=moments, termination=[])
    with contextlib.redirect_stdout(io.StringIO()):
        rs.initialize_model(core, [], [], [])
    rs.test_core = core
    return rs


def _simulate(rs, t_end):
    rs.termination.append(TerminationTime((t_end, "s")))
    ms = ModelSettings(tol_keep_in_edge=0.0, tol_move_to_core=1.0e8, tol_interrupt_simulation=1.0e8)
    with contextlib.redirect_stdout(io.StringIO()):
        rs.simulate(list(rs.test_core), [], [], [], [], [], model_settings=ms,
                    simulator_settings=SimulatorSettings())


def _index(rs, label):
    return [p.label for p in rs.polymer_pools].index(label)


def _set(rs, label, moments):
    for i, m in zip(rs.polymer_pools[_index(rs, label)].mu_indices, moments):
        rs.y[i] = m


def _moments(rs, label):
    return [float(rs.y[i]) for i in rs.polymer_pools[_index(rs, label)].mu_indices]


class _KickedPairSystem(HybridPolymerSystem):
    """The deck with PP's noise put into the empty pair by hand: over the
    first T_KICK s a source drives Rm and Ra to the refused PP state, then
    stops.  (On PP, DASPK puts it there itself; on five pools it does not.)"""

    T_KICK = 1.0e-6

    def residual(self, t, y, dydt, senpar=None):
        if senpar is None:
            delta, ires = HybridPolymerSystem.residual(self, t, y, dydt)
        else:
            delta, ires = HybridPolymerSystem.residual(self, t, y, dydt, senpar)
        if t < self.T_KICK:
            for label, target in (("Rm", PP_RM), ("Ra", PP_RA)):
                for i, m in zip(self.polymer_pools[_index(self, label)].mu_indices, target):
                    delta[i] += m / self.T_KICK
        return delta, ires


class _SinkOnEmptyPoolSystem(HybridPolymerSystem):
    """A sign error on the empty pool Ra: a constant sink that ignores what is
    left, past every error scale and every feed (-1e-6 mol/s in each
    moment).  Ra's linear drains hold a sink at about -RATE/k (k ~ 2.5e6
    s^-1 for the pair), so it must outrun k times the 1e-14 mol floor."""

    RATE = 1.0e-6

    def residual(self, t, y, dydt, senpar=None):
        if senpar is None:
            delta, ires = HybridPolymerSystem.residual(self, t, y, dydt)
        else:
            delta, ires = HybridPolymerSystem.residual(self, t, y, dydt, senpar)
        for i in self.polymer_pools[_index(self, "Ra")].mu_indices:
            delta[i] -= self.RATE
        return delta, ires


class TestEmptyPairNoise:

    def test_the_refused_pp_state_on_the_empty_pair_is_accepted(self):
        rs = pair_deck()
        _set(rs, "Rm", PP_RM)
        _set(rs, "Ra", PP_RA)
        rs._assert_pool_moments_accepted()

    def test_the_kicked_empty_pair_does_not_stop_the_run(self):
        """The pair holds PP's refused state from T_KICK on, its drains
        skipped while negative; the run raised there.  It reaches its end
        time, and the feed refills the pair."""
        rs = pair_deck(cls=_KickedPairSystem)
        _simulate(rs, 1.0)
        assert rs.t >= 1.0
        assert all(m > 0.0 for m in _moments(rs, "Rm") + _moments(rs, "Ra"))

    def test_the_empty_pair_is_carried_not_ended(self):
        """Between the floor that ends a pool and the error scale, an empty
        pool is carried: not ended, not refused."""
        rs = pair_deck()
        _set(rs, "Rm", PP_RM)
        _set(rs, "Ra", PP_RA)
        rs._update_pool_mu_peak()
        assert rs.end_exhausted_pools() == []
        assert _moments(rs, "Ra") == list(PP_RA)
        rs._assert_pool_moments_accepted()


class TestRefusalFloor:

    def test_the_refusal_floor_grows_with_the_moment_order_by_the_decks_chain_length(self):
        rs = pair_deck()
        # P is the trusted pool: mu2/mu1 = (DP^2 + DP)/DP
        assert rs._pool_length_scale == pytest.approx(DP + 1.0, rel=1e-15)
        floors = rs.pool_moment_refusal_floors()
        accepted = rs.pool_moment_accepted_floors()
        length = rs._pool_length_scale
        for p in range(len(rs.polymer_pools)):
            for k in range(3):
                assert floors[p, k] == max(accepted[p, k], rs._pool_mu_floors[p, 0] * length ** k)
        # a held pool keeps its relative floor at every order: P's mu2 floor is rtol*peak|mu2|
        assert floors[_index(rs, "P"), 2] == accepted[_index(rs, "P"), 2]
        assert list(floors[_index(rs, "Ra")]) == [1.0e-14, 1.0e-14 * length, 1.0e-14 * length * length]

    def test_an_untrusted_pool_does_not_set_the_length_scale(self):
        """Noise with a large mu2/mu1 in an empty pool, or a non-realizable
        triple above the floor, leaves L alone."""
        rs = pair_deck()
        length = rs._pool_length_scale
        _set(rs, "Rm", (5.0e-15, 1.0e-14, 1.0e-10))         # mu0 inside its floor
        _set(rs, "Ra", (1.0e-12, 2.0e-12, 3.0e-12))         # mu0*mu2 < mu1^2
        rs._update_pool_mu_peak()
        assert rs._pool_length_scale == length
        _set(rs, "Ra", (1.0e-12, 5.0e-11, 1.0e-8))          # a trusted pool of 200-mers
        rs._update_pool_mu_peak()
        assert rs._pool_length_scale == pytest.approx(200.0, rel=1e-15)

    def test_inside_the_refusal_floor_is_accepted_beyond_it_refused(self):
        rs = pair_deck()
        p = _index(rs, "Ra")
        floors = rs.pool_moment_refusal_floors()[p]
        for k in range(3):
            i = rs.polymer_pools[p].mu_indices[k]
            rs.y[i] = -0.99 * floors[k]
            rs._assert_pool_moments_accepted()
            rs.y[i] = -1.01 * floors[k]
            with pytest.raises(ValueError, match="beyond the exhaustion floor"):
                rs._assert_pool_moments_accepted()
            rs.y[i] = 0.0

    def test_the_floors_that_end_a_pool_are_unchanged(self):
        rs = pair_deck()
        _set(rs, "Ra", PP_RA)
        rs._update_pool_mu_peak()
        accepted = rs.pool_moment_accepted_floors()
        expected = np.maximum(rs._pool_mu_floors, rs._pool_mu_rtol * rs._pool_mu_peak)
        assert np.array_equal(accepted, expected)

    def test_a_sign_error_on_the_empty_pool_is_still_refused(self):
        rs = pair_deck(cls=_SinkOnEmptyPoolSystem)
        with pytest.raises(ValueError, match="beyond the exhaustion floor"):
            _simulate(rs, 1.0)
        # the sink takes mu0 past 1e-14 mol within ~1e-7 s
        assert rs.t < 1.0e-3
