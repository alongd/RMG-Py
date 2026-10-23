"""
A polymer pool that depolymerizes completely ends at zero moments, and the
integrator's noise there is relative to the pool's own scale, not to the
absolute tolerance.  The accepted-state negative-moment check refused any
moment below -100*atol (1e-14 mol at the default atol) whatever the pool's
initial moments (mu2(0) ~ 2e3 mol): DASPK's RMS error norm lets one moment
slot carry several error weights (the moment slots' DASPK atol is that same
1e-14), so a pool that empties completely while the reactor still holds other
states raised a ValueError on round-off.

The accepted floor is relative (per moment order k, -max(100*atol,
rtol*peak|mu_k|)), and a pool whose every moment falls inside it is ended by a
named ``pool_exhausted`` event: moments clamped to zero, the residue booked on
a ledger so the element and mass balances stay exact, the reactor continuing
with the rest.  A truly negative state is still refused.
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
DP = 1000.0
# a deck-scale pool: mu2(0) = 2e3 mol
SMALL = (1.0e-3, DP * 1.0e-3, 2.0 * DP * DP * 1.0e-3)
LARGE = (1.0, DP, 2.0 * DP * DP)
GAS_SPECIES = ["C#C", "O", "N", "S", "C=O", "CO", "CC", "OO", "NN", "C#N", "[He]", "[Ne]", "[Ar]",
               "CF", "CCl", "FF", "ClCl", "O=O", "C=C", "CN"]
REPEAT = (("C", 3), ("H", 6))
ENDS = (("H", 1), ("H", 1))


def _spc(smiles, label):
    s = Species(molecule=[Molecule(smiles=smiles)])
    s.label = label
    return s


def _rate(k):
    """An Arrhenius triplet with k(T) == k exactly."""
    return dict(A=float(k), n=0.0, Ea=0.0)


def _pool(label, mu_indices, **kw):
    cfg = dict(xs=2, explicit_dp_to_species_index={}, k_scission=0.0, k_unzip=0.0, tail_kinetics=None,
               monomer_mw_g_mol=42.08, repeat_unit_atoms=REPEAT, end_group_atoms=ENDS)
    cfg.update(kw)
    return PolymerPoolConfig(label=label, mu_indices=mu_indices, **cfg)


def end_unzip(k_dep, moments, n_gas=0):
    """One pool unzipping from its chain ends (end-radical depropagation;
    the last unit of a chain leaves with its chain), releasing gas monomer,
    beside ``n_gas`` other gas species."""
    core = [_spc("N#N", "N2"), _spc("C=CC", "monomer"), _spc("[CH2]CC", "poly"),
            _spc("CO", "poly_mu0"), _spc("C=O", "poly_mu1"), _spc("C#N", "poly_mu2")]
    core += [_spc(s, f"gas{i}") for i, s in enumerate(GAS_SPECIES[:n_gas])]
    pools = [_pool("poly", (3, 4, 5), monomer_poly_index=1, k_depropagation=_rate(k_dep))]
    return core, pools, {"poly": tuple(moments)}


def scission_unzip(k_hom, k_dep, moments):
    """A pool cut at random backbone bonds into two end-radical daughter
    pools that unzip (the parent empties first, its daughters still full)."""
    core = [_spc("N#N", "N2"), _spc("C=CC", "monomer"), _spc("CCCC", "poly"),
            _spc("CO", "poly_mu0"), _spc("C=O", "poly_mu1"), _spc("C#N", "poly_mu2"),
            _spc("[CH2]CC", "poly_rad_primary_end"), _spc("CCO", "poly_rad_primary_end_mu0"),
            _spc("CC=O", "poly_rad_primary_end_mu1"), _spc("CC#N", "poly_rad_primary_end_mu2"),
            _spc("C[CH]C", "poly_rad_secondary_end"), _spc("CCCO", "poly_rad_secondary_end_mu0"),
            _spc("CCC=O", "poly_rad_secondary_end_mu1"), _spc("CCC#N", "poly_rad_secondary_end_mu2")]
    pools = [_pool("poly", (3, 4, 5), monomer_poly_index=None, k_homolysis=_rate(k_hom))]
    for label, mu in (("poly_rad_primary_end", (7, 8, 9)), ("poly_rad_secondary_end", (11, 12, 13))):
        pools.append(_pool(label, mu, monomer_poly_index=1, k_depropagation=_rate(k_dep)))
    moments = {"poly": tuple(moments), "poly_rad_primary_end": (0.0, 0.0, 0.0),
               "poly_rad_secondary_end": (0.0, 0.0, 0.0)}
    return core, pools, moments


def _system(core, pools, moments, cls=HybridPolymerSystem):
    mask = np.array([not any(s.label == p.label or s.label.startswith(p.label + "_mu") for p in pools)
                     for s in core], dtype=bool)
    rs = cls(T=T_K, P=1.0e5, initial_mole_fractions={core[0]: 1.0}, V_poly=1.0, polymer_pools=pools,
             mass_transfer=[], gas_species_mask=mask, constant_gas_volume=False,
             initial_polymer_moments=moments, termination=[])
    with contextlib.redirect_stdout(io.StringIO()):
        rs.initialize_model(core, [], [], [])
    return rs


def _simulate(rs, core, t_end):
    rs.termination.append(TerminationTime((t_end, "s")))
    ms = ModelSettings(tol_keep_in_edge=0.0, tol_move_to_core=1.0e8, tol_interrupt_simulation=1.0e8)
    with contextlib.redirect_stdout(io.StringIO()):
        rs.simulate(list(core), [], [], [], [], [], model_settings=ms, simulator_settings=SimulatorSettings())


def _moments(rs, label):
    pool = next(p for p in rs.polymer_pools if p.label == label)
    return [float(rs.y[i]) for i in pool.mu_indices]


def _closes(rs, y0_inventory):
    """The code-convention elements now, plus the exhausted residue, equal
    those at t = 0 (gas monomer released keeps the repeat units)."""
    now = rs.code_convention_inventory()
    for el, n in rs.exhausted_residue_atoms().items():
        now[el] = now.get(el, 0.0) + n
    for el, n0 in y0_inventory.items():
        assert now.get(el, 0.0) == pytest.approx(n0, rel=1e-7), (el, now, y0_inventory)


class TestCompleteDepolymerization:

    def test_a_parent_emptied_by_scission_ends_by_event_and_its_daughters_go_on(self):
        """The parent is cut into two unzipping daughters and empties at
        t ~ 0.045 s, its round-off (mu2 ~ -2e-14 mol) below the absolute
        floor: the run raised there.  Its moments are inside the relative
        floor, so it ends by pool_exhausted; the daughters keep unzipping
        until they are exhausted in turn."""
        core, pools, moments = scission_unzip(10.0, 1.0e3, SMALL)
        rs = _system(core, pools, moments)
        inv0 = rs.code_convention_inventory(rs.y0)
        _simulate(rs, core, 20.0)
        assert rs.t >= 20.0
        events = {e["pool"]: e for e in rs.pool_exhaustion_events}
        assert set(events) == {"poly", "poly_rad_primary_end", "poly_rad_secondary_end"}, events
        parent = events["poly"]
        assert parent["event"] == "pool_exhausted"
        assert 0.0 < parent["time_s"] < 0.1
        assert all(abs(m) <= f for m, f in zip(parent["residual_moments_mol"], parent["floors_mol"]))
        assert all(e["time_s"] > parent["time_s"] for label, e in events.items() if label != "poly")
        assert _moments(rs, "poly") == [0.0, 0.0, 0.0]
        _closes(rs, inv0)

    def test_an_end_unzipping_pool_beside_many_gas_species_ends_by_event(self):
        """Twenty more gas species: the RMS error norm lets the emptied pool's
        mu2 slot reach -1.3e-14 mol on an accepted state (raised at t ~ 51 s).
        One event, at the pool's end, and the run reaches its end time."""
        core, pools, moments = end_unzip(1.0e3, LARGE, n_gas=20)
        rs = _system(core, pools, moments)
        inv0 = rs.code_convention_inventory(rs.y0)
        _simulate(rs, core, 100.0)
        assert rs.t >= 100.0
        assert [e["pool"] for e in rs.pool_exhaustion_events] == ["poly"]
        assert _moments(rs, "poly") == [0.0, 0.0, 0.0]
        _closes(rs, inv0)


class _SignErrorSystem(HybridPolymerSystem):
    """A pool drained by a constant sink that ignores what is left: once its
    moments reach zero it keeps going, to -mu_k(0) at t = 2 TAU."""

    TAU = 0.5

    def residual(self, t, y, dydt, senpar=None):
        if senpar is None:
            delta, ires = HybridPolymerSystem.residual(self, t, y, dydt)
        else:
            delta, ires = HybridPolymerSystem.residual(self, t, y, dydt, senpar)
        for i in self.polymer_pools[0].mu_indices:
            delta[i] -= self.y0[i] / self.TAU
        return delta, ires


class TestTrulyNegativeStateRefused:

    def test_a_sign_error_is_refused_on_the_accepted_state(self):
        core, pools, moments = end_unzip(1.0e-6, SMALL)
        rs = _system(core, pools, moments, cls=_SignErrorSystem)
        with pytest.raises(ValueError, match="beyond the exhaustion floor"):
            _simulate(rs, core, 5.0)
        assert rs.t < 2.0 * _SignErrorSystem.TAU


class TestRelativeFloor:

    def test_the_accepted_floor_is_the_larger_of_the_absolute_and_relative_terms(self):
        core, pools, moments = end_unzip(1.0e3, LARGE)
        rs = _system(core, pools, moments)
        floors = rs.pool_moment_accepted_floors()
        absolute = np.asarray(rs._pool_mu_floors)
        rtol = 1.0e-8   # the system's rtol (initialize_model's default)
        assert floors.shape == (1, 3)
        for k, mu in enumerate(LARGE):
            assert floors[0, k] == max(absolute[0, k], rtol * mu)
        # a born-empty pool keeps the absolute floor
        core, pools, moments = scission_unzip(10.0, 1.0e3, SMALL)
        rs = _system(core, pools, moments)
        floors = rs.pool_moment_accepted_floors()
        assert list(floors[1]) == list(np.asarray(rs._pool_mu_floors)[1])

    def test_inside_the_relative_floor_is_accepted_beyond_it_refused(self):
        core, pools, moments = end_unzip(1.0e3, SMALL)
        rs = _system(core, pools, moments)
        mu2 = pools[0].mu_indices[2]
        f2 = rs.pool_moment_accepted_floors()[0, 2]
        rs.y[mu2] = -0.5 * f2
        rs._assert_pool_moments_accepted()
        rs.y[mu2] = -2.0 * f2
        with pytest.raises(ValueError, match="beyond the exhaustion floor"):
            rs._assert_pool_moments_accepted()

    def test_the_absolute_floors_the_other_consumers_read_are_unchanged(self):
        core, pools, moments = end_unzip(1.0e3, LARGE)
        rs = _system(core, pools, moments)
        assert list(np.asarray(rs._pool_mu_floors)[0]) == [1.0e-14] * 3


class TestExhaustionEvent:

    def _in_band(self, rs, label="poly"):
        pool = next(p for p in rs.polymer_pools if p.label == label)
        floors = rs.pool_moment_accepted_floors()[rs.polymer_pools.index(pool)]
        residue = (0.25 * floors[0], -0.5 * floors[1], 0.75 * floors[2])
        for i, m in zip(pool.mu_indices, residue):
            rs.y[i] = m
        return pool, residue

    def test_the_event_clamps_the_pool_and_books_its_residue(self):
        core, pools, moments = end_unzip(1.0e3, LARGE)
        rs = _system(core, pools, moments)
        pool, residue = self._in_band(rs)
        before = rs.physical_elemental_inventory()
        rs.end_exhausted_pools()
        assert [float(rs.y[i]) for i in pool.mu_indices] == [0.0, 0.0, 0.0]
        (event,) = rs.pool_exhaustion_events
        assert event["event"] == "pool_exhausted" and event["pool"] == "poly"
        assert event["time_s"] == rs.t
        assert event["residual_moments_mol"] == list(residue)
        assert event["residual_mass_g"] == pool.condensed_mass_g(residue[0], residue[1])
        assert rs.exhausted_residue_mass_g() == {"poly": event["residual_mass_g"]}
        after = rs.physical_elemental_inventory()
        for el, n in rs.exhausted_residue_atoms().items():
            after[el] = after.get(el, 0.0) + n
        for el in before:
            assert after[el] == pytest.approx(before[el], rel=1e-15, abs=1e-30), el

    def test_an_ended_pool_is_ended_once_and_a_live_pool_never(self):
        core, pools, moments = end_unzip(1.0e3, LARGE)
        rs = _system(core, pools, moments)
        rs.end_exhausted_pools()
        assert rs.pool_exhaustion_events == []
        self._in_band(rs)
        rs.end_exhausted_pools()
        self._in_band(rs)
        rs.end_exhausted_pools()
        assert len(rs.pool_exhaustion_events) == 1

    def test_a_born_empty_pool_is_not_exhausted(self):
        core, pools, moments = scission_unzip(10.0, 1.0e3, SMALL)
        rs = _system(core, pools, moments)
        rs.end_exhausted_pools()
        assert rs.pool_exhaustion_events == []
