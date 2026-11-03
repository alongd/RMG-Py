"""The kMC initial state, leak classes and observables of a network
declaration (S7).

``materialize`` runs ``random_bond_v1`` above the gel point: one giant
component plus sol, the realized links equal to the declared ones, the gel
fraction and branch points reported.  At runtime a network reports leak per
class: an undeclared link, a bond closing a cycle shorter than
``minimum_cycle_bonds``, a bimolecular pair on one molecule (measured, never
fired) and the gel component's whole-molecule mobility (it never moves as a
whole: its diffusivity enters no escape factor, and what it would add is
measured).  The observables read the gel as the largest component, the rest
of the condensed melt as sol, and an inert filler as residue mass."""

import math
from collections import Counter

import numpy as np
import pytest

import melt_engine_support as S
from branched_support import branched_spec
from meltEngineTest import _pair_counts, _recombination
from branchedMeltTest import _closure, _cycle_through, _inter_residue_bonds, _uni
from network_support import counts, network_spec, no_filler, with_met

from rmgpy.kmc import melt_engine as me


# --------------------------------------------------------------------------
# random_bond_v1 above the gel point
# --------------------------------------------------------------------------
def test_the_generator_gives_one_giant_component_plus_sol():
    from rmgpy.kmc import melt_branched as mb
    spec = network_spec()
    report = mb.realized(spec)
    sizes = mb.random_bond(spec).molecules
    largest = max(len(m) for m in sizes)
    assert report["gel_residues"] == largest
    assert report["gel_fraction"] == pytest.approx(largest / 76)
    assert report["gel_fraction"] > 0.5
    assert report["sol_molecules"] == len(sizes) - 1
    assert report["link_counts"] == report["declared_link_counts"] == {"m_o": 52, "m_p": 20}
    assert report["gel_branching_ratio"] > 1.0
    assert isinstance(report["branch_points"], int) and report["branch_points"] > 0
    assert report["expected_branch_points"] > 0


def test_materialize_realizes_the_declared_links_exactly():
    from rmgpy.kmc.melt_state import materialize
    state = materialize(network_spec())
    assert len(_inter_residue_bonds(state)) == 72


def test_no_cycle_in_the_network_is_shorter_than_declared():
    from rmgpy.kmc.melt_state import materialize
    spec = network_spec()
    state = materialize(spec)
    cycles = [_cycle_through(state, a, b) for a, b in _inter_residue_bonds(state)]
    closed = [c for c in cycles if c is not None]
    assert closed, "a network above the gel point closes cycles"
    assert min(closed) >= spec.minimum_cycle_bonds


def test_the_branched_report_is_unchanged():
    from rmgpy.kmc import melt_branched as mb
    assert "gel_fraction" not in mb.realized(branched_spec())


# --------------------------------------------------------------------------
# leak per class
# --------------------------------------------------------------------------
def _engine(spec, state, templates):
    return me.MeltEngine(spec, S.ToyCatalog(builder=lambda _t: list(templates)), 700.0,
                         np.random.default_rng(3), state=state)


def _melt(spec_of=network_spec):
    """mol0 Ph-M-Ph-M-Ph (o2-o2, o6-o2) with M3-Ph4 cut: a bridge radical on
    {Ph0 M1 Ph2 M3} and an aryl radical alone; mol1 Ph-M-Ph (o2-o2).  Then
    an H off the vacant para site of Ph0 on both molecules: two aryl radicals
    of one window, one on the bridge radical's molecule."""
    from rmgpy.kmc import melt_branched as mb
    from rmgpy.kmc.melt_state import SetBond, SetRadical

    def mutate(mapping):
        with_met(mapping)
        no_filler(mapping)
    spec = spec_of(mutate)
    state = mb.build_molecules(spec, [
        (["Ph", "M", "Ph", "M", "Ph"], [(1, "a", 0, "o2"), (1, "b", 2, "o2"),
                                        (3, "a", 2, "o6"), (3, "b", 4, "o2")]),
        (["Ph", "M", "Ph"], [(1, "a", 0, "o2"), (1, "b", 2, "o2")])])
    atom = mb.residue_atoms(state)
    m, r = atom[(0, 3)]["m"], atom[(0, 4)]["r2"]
    state.apply([SetBond(m, r, 1.0, 0.0), SetRadical(m, 0, 1), SetRadical(r, 0, 1)])
    for molecule in (0, 1):
        S.strip_h(state, atom[(molecule, 0)]["r4"])
    return spec, state, atom


def _classes(engine):
    return dict(engine.propensities().leak_classes)


def test_a_bimolecular_pair_on_one_molecule_is_leak_and_never_fires():
    spec, state, atom = _melt()
    bridge, mine, other = atom[(0, 3)]["m"], atom[(0, 0)]["r4"], atom[(1, 0)]["r4"]
    engine = _engine(spec, state, [_recombination(state, bridge, other, "m_p")])
    allowed, leak = _pair_counts(engine, bridge, other, "m_p")
    assert allowed > 0 and leak == 0
    index = [t.name for t in engine.templates].index("m_p")
    classes = _classes(engine)
    assert classes["same_molecule"] == pytest.approx(engine.weights[index] / engine.nav)
    # the same-molecule pair is no channel: the MET channel holds the other molecule's pair only
    comps = engine.state.components
    assert set(engine.pairs) == {tuple(sorted((comps.component_of(bridge), comps.component_of(other))))}
    assert comps.component_of(mine) == comps.component_of(bridge)
    engine.run(max_events=1)
    assert engine.template_counts.get("m_p", 0) == 1
    assert engine.state.components.component_of(bridge) == engine.state.components.component_of(other)


def test_a_bond_closing_a_short_cycle_is_leak_of_its_own_class():
    """The bridge radical closing on Ph2's vacant para site: Ph2-M3 ring of
    fewer than eight bonds, at two vacant sites of a declared link."""
    spec, state, atom = _melt()
    ring = atom[(0, 2)]["r4"]
    S.strip_h(state, ring)
    engine = _engine(spec, state, [_closure(state, ring, atom[(0, 3)]["m"], "short")])
    allowed, leak = _uni(engine, "short")
    assert allowed == 0 and leak > 0
    classes = _classes(engine)
    assert classes["short_cycle"] > 0 and classes["undeclared_link"] == 0


def test_a_bond_at_an_undeclared_site_is_leak_of_its_own_class():
    spec, state, atom = _melt()
    ring = atom[(1, 2)]["r3"]
    S.strip_h(state, ring)
    engine = _engine(spec, state, [_recombination(state, atom[(0, 3)]["m"], ring, "off_site")])
    assert _pair_counts(engine, atom[(0, 3)]["m"], ring, "off_site") == (0, 1)
    classes = _classes(engine)
    assert classes["undeclared_link"] > 0 and classes["short_cycle"] == 0


def test_the_gel_component_never_moves_as_a_whole():
    """The largest component (Ph0 M1 Ph2 M3) has no whole-molecule
    diffusivity in the escape factor; what it would add is mobility leak."""
    spec, state, atom = _melt()
    bridge, other = atom[(0, 3)]["m"], atom[(1, 0)]["r4"]
    engine = _engine(spec, state, [_recombination(state, bridge, other, "m_p")])
    comps = engine.state.components
    key = tuple(sorted((comps.component_of(bridge), comps.component_of(other))))
    rate = engine.pairs[key]
    gel, sol = engine.entries[comps.component_of(bridge)], engine.entries[comps.component_of(other)]
    assert engine.gel_component() == gel.cid
    still = engine.transport.k_diffusion_still(engine.temperature, sol.met_length, gel.met_length)
    full = engine.transport.k_diffusion(engine.temperature, sol.met_length, gel.met_length)
    assert 0.0 < still < full
    assert rate.k_d == pytest.approx(still)
    phi, phi_full = me.escape_factor(rate.k_act, still), me.escape_factor(rate.k_act, full)
    assert rate.allowed == pytest.approx(phi * rate.k_act_allowed / engine.nav)
    assert rate.mobility == pytest.approx((phi_full - phi) * rate.k_act_allowed / engine.nav)
    assert _classes(engine)["gel_mobility"] == pytest.approx(rate.mobility)


def test_leak_classes_integrate_to_the_leak_number():
    spec, state, atom = _melt()
    bridge, other = atom[(0, 3)]["m"], atom[(1, 0)]["r4"]
    engine = _engine(spec, state, [_recombination(state, bridge, other, "m_p")])
    engine.run(until=1.0e-9)
    report = engine.report()
    by_class = report["leak_by_class"]
    assert set(by_class) == {"undeclared_link", "short_cycle", "same_molecule", "gel_mobility"}
    assert math.fsum(by_class.values()) == pytest.approx(report["leak_num"])
    assert by_class["same_molecule"] > 0


def test_a_branched_engine_reports_no_leak_classes():
    spec, state, atom = _melt(branched_spec)
    bridge, other = atom[(0, 3)]["m"], atom[(1, 0)]["r4"]
    engine = _engine(spec, state, [_recombination(state, bridge, other, "m_p")])
    assert "leak_by_class" not in engine.report()
    assert engine.propensities().leak_classes == ()
    comps = engine.state.components
    key = tuple(sorted((comps.component_of(bridge), comps.component_of(other))))
    rate = engine.pairs[key]
    assert rate.mobility == 0.0
    assert rate.k_d == pytest.approx(engine.transport.k_diffusion(
        engine.temperature, engine.entries[key[0]].met_length, engine.entries[key[1]].met_length))


# --------------------------------------------------------------------------
# observables
# --------------------------------------------------------------------------
def test_sol_and_gel_mass_fractions_at_time_zero():
    from rmgpy.kmc import melt_branched as mb
    from rmgpy.kmc import melt_network as mn
    from rmgpy.kmc.melt_state import materialize
    spec = network_spec()
    obs = mn.observe(spec, materialize(spec))
    assert obs["gel_mass_fraction"] + obs["sol_mass_fraction"] == pytest.approx(1.0)
    assert 0.5 < obs["gel_mass_fraction"] < 1.0
    assert obs["gel_residues"] == mb.realized(spec)["gel_residues"]
    assert obs["mass_loss_fraction"] == pytest.approx(0.0, abs=1e-12)
    assert obs["residue_mass_fraction"] == pytest.approx(1.0)


def test_an_inert_filler_stays_as_residue_mass():
    """One H lost from the resin: the filler (half the composite's mass)
    stays as residue, so the composite loses half the resin's fraction."""
    from rmgpy.kmc import melt_network as mn
    from rmgpy.kmc.melt_state import materialize
    losses = {}
    for name, mutate in (("filled", None), ("bare", no_filler)):
        spec = network_spec(mutate)
        state = materialize(spec)
        obs0 = mn.observe(spec, state)
        ring = next(a for a, o in state.origin.items() if o.key == "r3" and state.atoms[a].element == "C")
        S.strip_h(state, ring)
        obs = mn.observe(spec, state)
        losses[name] = obs["mass_loss_fraction"]
        assert obs["inert_mass_g_mol"] == obs0["inert_mass_g_mol"]
    assert losses["bare"] > 0
    assert losses["filled"] == pytest.approx(losses["bare"] / 2)
    spec = network_spec()
    obs = mn.observe(spec, materialize(spec))
    assert obs["inert_mass_g_mol"] == pytest.approx(obs["initial_resin_mass_g_mol"])
    assert obs["inert_composition"] == {"filler": {"C": 1}}


def test_the_observation_names_its_masses():
    from rmgpy.kmc import melt_network as mn
    spec, state, _atom = _melt()
    obs = mn.observe(spec, state)
    # the two H atoms stripped to make the aryl radicals are the built-in species H
    assert set(obs["volatile_by_species_g_mol"]) == {"H"}
    assert obs["volatile_mass_g_mol"] == pytest.approx(2 * 1.00794, rel=1e-3)
    assert set(obs) >= {"gel_mass_fraction", "sol_mass_fraction", "mass_loss_fraction",
                        "residue_mass_fraction", "volatile_mass_g_mol", "sol_molecule_sizes"}


def test_an_observed_run_records_the_gel_fraction_and_leak_over_time():
    from rmgpy.kmc import melt_network as mn
    spec, state, atom = _melt()
    bridge, other = atom[(0, 3)]["m"], atom[(1, 0)]["r4"]
    engine = _engine(spec, state, [_recombination(state, bridge, other, "m_p")])
    rows = mn.observed_run(engine, [0.0, 1.0e-12, 1.0e3])
    assert [r["time"] for r in rows] == [0.0, 1.0e-12, 1.0e3]
    assert all(set(r["leak_by_class"]) == {"undeclared_link", "short_cycle", "same_molecule",
                                           "gel_mobility"} for r in rows)
    assert rows[-1]["events"] >= 1
    assert all(r["gel_mass_fraction"] + r["sol_mass_fraction"] == pytest.approx(1.0) for r in rows)


def test_the_network_index_after_events_equals_a_full_recompute():
    """Same-molecule pairs, short-cycle counts, the gel and its mobility are
    kept incrementally: after every event they equal a fresh rebuild's."""
    spec, state, atom = _melt()
    bridge, other = atom[(0, 3)]["m"], atom[(1, 0)]["r4"]
    ring = atom[(0, 2)]["r4"]
    S.strip_h(state, ring)
    templates = [_recombination(state, bridge, other, "m_p"), _closure(state, ring, bridge, "short")]
    engine = _engine(spec, state, templates)
    checked = []

    def check(e, _event):
        assert e.snapshot() == e.fresh_snapshot()
        checked.append(1)
    assert engine.snapshot() == engine.fresh_snapshot()
    engine.run(until=1.0e3, check=check)
    assert checked
