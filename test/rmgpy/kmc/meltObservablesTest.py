"""S5 red tests: melt observables comparable with the method-of-moments solver (I-154).

Design: I-127 v3 §5 (release and species identity) and §7 (observables, the
MoM comparison, exact invariants), §9 (S5 red tests).  Expected numbers are
derived independently of the observables module: the declaration's own
initialization (``macroscopic_moments_mol``, ``scale_factor``,
``ensemble_observables``), integer DP counts of the declared ensemble, exact
masses read from the state's component index, and hand-applied edits.
"""

import json
import math
import re
from collections import Counter
from dataclasses import replace
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

from rmgpy import constants
from rmgpy.polymer_mass import Composition

import melt_engine_support as S

from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import melt_observables as mo
from rmgpy.kmc.melt_state import AncestryError, MeltStateError, SetBond, SetRadical, materialize

PE_LIVE = ("homolysis", "recombination", "h_abstraction", "disproportionation", "beta_scission",
           "mid_end_recombination", "intra_mid_end")
REL = 1.0e-12


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def chain_carbons_in_order(state, ordinal):
    """Backbone carbons of chain ``ordinal`` from the head end to the tail end."""
    cs = [a for a in S.chain_atoms(state, ordinal) if state.origin[a].is_backbone]
    ends = [a for a in cs if sum(1 for b in state.atoms[a].edges if b in cs) == 1]
    path, prev = [min(ends)], None
    while len(path) < len(cs):
        nxt = [b for b in state.atoms[path[-1]].edges if b in cs and b != prev]
        prev = path[-1]
        path.append(nxt[0])
    return path


def homolyse(state, a, b):
    """Break the a-b single bond; both atoms gain a radical electron."""
    state.apply([SetBond(a, b, 1.0, 0.0), SetRadical(a, state.atoms[a].radical, state.atoms[a].radical + 1),
                 SetRadical(b, state.atoms[b].radical, state.atoms[b].radical + 1)])


def recombine(state, a, b):
    state.apply([SetBond(a, b, 0.0, 1.0), SetRadical(a, state.atoms[a].radical, state.atoms[a].radical - 1),
                 SetRadical(b, state.atoms[b].radical, state.atoms[b].radical - 1)])


# --------------------------------------------------------------------------
# C5': t = 0 export, physical mass moments, integer n_c
# --------------------------------------------------------------------------
def test_c5_t0_export_equals_the_declared_macroscopic_moments():
    counts = {10: 1, 20: 1, 30: 1}
    spec = S.ps_spec(counts)
    engine = me.MeltEngine(spec, S.ps_catalog(), 700.0, np.random.default_rng(0))
    o = mo.MeltObserver(spec, xs=0).observe(engine)
    init = spec.initialization()
    assert o["schema"] == mo.OBSERVABLES_SCHEMA and o["time_s"] == 0.0
    assert o["scale_factor"] == init.scale_factor and o["backbone_atoms_per_unit"] == 2
    t2 = o["tier2"]
    rep = t2["dp_moments_representative"]
    assert all(type(x) is int for x in rep)
    assert rep == [3, 60, sum(dp * dp for dp in counts)]
    mu = spec.macroscopic_moments_mol
    for got, want in zip(t2["dp_moments_mol"], mu):
        assert rel(got, want) <= REL, (got, want)
    # tier 1: polymer-derived mass and molecule count, kMC and MoM definition, equal the declaration
    p = o["tier1"]["polymer_derived"]
    assert p["molecules"] == 3
    assert rel(p["amount_mol"], mu[0]) <= REL
    assert rel(p["mass_g"], spec.initial_mass_kg * 1000.0) <= REL
    assert rel(p["mom_definition"]["amount_mol"], p["amount_mol"]) <= REL
    assert rel(p["mom_definition"]["mass_g"], p["mass_g"]) <= REL
    # mass moments are physical and equal the declared ensemble at t = 0
    assert t2["M1_g"] > 0 and t2["M0_mol"] * t2["M2_g2_mol"] >= t2["M1_g"] ** 2
    assert t2["Mw_g_mol"] >= t2["Mn_g_mol"] > 0
    ens = spec.ensemble_observables
    assert rel(t2["Mn_g_mol"], ens.mn_g_mol) <= REL and rel(t2["Mw_g_mol"], ens.mw_g_mol) <= REL
    assert rel(t2["mom_convention"]["Mn_g_mol"], ens.mn_g_mol) <= REL
    assert t2["end_mass_convention"]["M1_difference_g"] == 0.0
    assert t2["end_mass_convention"]["M1_bound_g"] == 0.0


def test_c5_t0_declared_mobile_species_amounts_are_scale_times_count_over_avogadro():
    spec = S.ps_spec({10: 1, 20: 1, 30: 1}, mobile={"styrene": 2}, monomer=True)
    state = materialize(spec)
    o = mo.MeltObserver(spec, xs=0).observe(state)
    scale = spec.initialization().scale_factor
    sty = o["tier1"]["species"]["styrene"]
    styrene_graph = dict((n, g) for n, g, _d in spec.small_species)["styrene"]
    assert sty["molecules"] == 2
    assert rel(sty["amount_mol"], scale * 2 / constants.Na) <= REL
    assert rel(sty["mass_g"], scale * 2 / constants.Na * styrene_graph.composition.molar_mass_g_mol()) <= REL
    # the mobile monomer carries backbone ancestry but is a declared species, not polymer-derived
    p = o["tier1"]["polymer_derived"]
    assert p["molecules"] == 3 and rel(p["mass_g"], spec.initial_mass_kg * 1000.0) <= REL
    assert rel(o["tier1"]["total_mass_g"], p["mass_g"] + sty["mass_g"]) <= REL


def test_n_c_is_the_integer_floor_of_backbone_atoms_over_b():
    spec = S.ps_spec({10: 1})
    state = materialize(spec)
    path = chain_carbons_in_order(state, 0)
    # break the CH2-CH(Ph) bond inside unit 4: both fragments carry an odd backbone count
    homolyse(state, path[6], path[7])
    obs = mo.MeltObserver(spec, xs=3)
    records = obs.components(state)
    assert sorted(r.backbone for r in records) == [7, 13]
    for r in records:
        assert type(r.n) is int and r.n == r.backbone // 2
    o = obs.observe(state)
    assert o["tier2"]["dp_moments_representative"] == [2, 3 + 6, 9 + 36]
    assert all(type(x) is int for x in o["tier2"]["dp_moments_representative"])
    t3 = o["tier3"]
    assert t3["diagnostic_only"] is True and t3["xs"] == 3
    assert t3["tail"]["dp_moments_representative"] == [1, 6, 36]
    assert t3["explicit"]["dp_molecules"] == {"3": 1}
    mu0, mu1, _ = t3["tail"]["dp_moments_mol"]
    assert mu1 >= (3 + 1) * mu0


def test_tier3_partition_uses_integer_n_and_xs():
    spec = S.ps_spec({10: 1, 20: 1, 30: 1})
    o = mo.MeltObserver(spec, xs=15).observe(materialize(spec))
    scale = spec.initialization().scale_factor
    t3 = o["tier3"]
    assert t3["tail"]["molecules"] == 2 and t3["tail"]["dp_moments_representative"] == [2, 50, 1300]
    assert t3["explicit"]["dp_molecules"] == {"10": 1}
    assert rel(t3["explicit"]["dp_amount_mol"]["10"], scale / constants.Na) <= REL
    m = mo.mom_reported_quantities(o)
    assert rel(m["polymer_derived_molecules_mol"], o["tier1"]["polymer_derived"]["amount_mol"]) <= REL
    assert rel(m["mu0_mol"] + m["explicit_dp_mol"]["10"], spec.macroscopic_moments_mol[0]) <= REL


def test_end_mass_convention_difference_is_reported_with_a_bound():
    spec = S.pe_spec({6: 1})
    state = materialize(spec)
    end = S.end_carbons(state, 0)[0]
    nxt = next(b for b in state.atoms[end].edges if state.atoms[b].element == "C")
    homolyse(state, end, nxt)       # CH3. (B = 1, n = 0) + C11H23. (B = 11, n = 5)
    o = mo.MeltObserver(spec, xs=0).observe(state)
    conv = o["tier2"]["end_mass_convention"]
    scale = spec.initialization().scale_factor
    ch = Composition.from_mapping({"C": 1, "H": 1}).molar_mass_g_mol()
    assert rel(conv["delta_max_g_mol"], ch) <= 1e-12 and rel(conv["delta_min_g_mol"], ch) <= 1e-12
    assert rel(conv["M1_difference_g"], scale / constants.Na * ch) <= 1e-12
    assert abs(conv["M1_difference_g"]) <= conv["M1_bound_g"]
    assert abs(conv["M2_difference_g2_mol"]) <= conv["M2_bound_g2_mol"]
    t2 = o["tier2"]
    assert rel(t2["M1_g"] - t2["mom_convention"]["mass_g"], conv["M1_difference_g"]) <= 1e-9


# --------------------------------------------------------------------------
# exact invariants after every event
# --------------------------------------------------------------------------
def test_invariants_hold_after_every_event_of_a_200_event_ps_trajectory():
    spec = S.ps_spec({10: 2, 6: 2})
    engine = me.MeltEngine(spec, S.ps_catalog(), 700.0, np.random.default_rng(3))
    checker = mo.InvariantChecker(engine.state)
    result = engine.run(max_events=200, check=checker)
    assert result.events == 200 and checker.checks == 200
    assert checker.counts["split"] >= 1 and checker.counts["merge"] >= 1, checker.counts
    assert sum(checker.counts.values()) == 200


def test_bridge_scission_and_recombination_moment_increments_are_exact():
    spec = S.pe_spec({10: 1, 6: 1})
    state = materialize(spec)
    path = chain_carbons_in_order(state, 0)
    before = mo.exact_totals(state)
    checker = mo.InvariantChecker(state)
    (cid,) = {state.components.component_of(a) for a in path}
    homolyse(state, path[7], path[8])
    after = mo.exact_totals(state)
    ca, cb = state.components.component_of(path[7]), state.components.component_of(path[8])
    ma, mb = state.components.molar_mass_g_mol(ca), state.components.molar_mass_g_mol(cb)
    assert isinstance(ma, Fraction)
    assert after.count == before.count + 1
    assert after.mass2 - before.mass2 == -2 * ma * mb
    assert after.mass == before.mass and after.composition == before.composition
    assert after.backbone == before.backbone
    checker.check(state, me.Event(0.0, "hand", "hand", "uni", (cid,), tuple(sorted((ca, cb))), state.revision))
    recombine(state, path[7], path[8])
    back = mo.exact_totals(state)
    assert back.count == after.count - 1 and back.mass2 - after.mass2 == 2 * ma * mb
    assert back == before
    checker.check(state, me.Event(0.0, "hand", "hand", "met", tuple(sorted((ca, cb))), (cid,), state.revision))
    assert checker.counts == Counter({"split": 1, "merge": 1})


def test_invariant_checker_rejects_a_wrong_scission_increment():
    spec = S.pe_spec({10: 1})
    state = materialize(spec)
    path = chain_carbons_in_order(state, 0)
    checker = mo.InvariantChecker(state)
    cid = state.components.component_of(path[0])
    homolyse(state, path[3], path[4])
    products = tuple(sorted({state.components.component_of(a) for a in path}))
    # an event that claims an unchanged molecule count is caught
    with pytest.raises(mo.InvariantError):
        checker.check(state, me.Event(0.0, "x", "x", "uni", (cid,), products[:1], state.revision))


def test_invariant_checker_catches_a_wrong_sum_m2_increment():
    """With exact mass conservation and the right component sets the bridge
    increment -2 M_a M_b is implied; the checker's M^2 line catches index mass
    that moved between the event's products (total mass unchanged)."""
    spec = S.pe_spec({10: 1})
    state = materialize(spec)
    path = chain_carbons_in_order(state, 0)
    checker = mo.InvariantChecker(state)
    cid = state.components.component_of(path[0])
    homolyse(state, path[3], path[4])
    a, b = sorted({state.components.component_of(x) for x in path})
    rec = state.components._rec
    shift = Fraction(1, 1000)
    rec[a] = replace(rec[a], mass=rec[a].mass + shift)
    rec[b] = replace(rec[b], mass=rec[b].mass - shift)
    with pytest.raises(mo.InvariantError, match=r"M\^2"):
        checker.check(state, me.Event(0.0, "x", "x", "uni", (cid,), (a, b), state.revision))


@pytest.mark.parametrize("field", ["backbone", "composition"])
def test_invariant_checker_catches_index_drift(field):
    spec = S.pe_spec({10: 1, 6: 1})
    state = materialize(spec)
    checker = mo.InvariantChecker(state)
    checker.check(state)
    cid = state.components.ids()[0]
    rec = state.components._rec[cid]
    if field == "backbone":
        state.components._rec[cid] = replace(rec, backbone=rec.backbone + 1)
    else:
        comp = dict(rec.composition)
        comp[("H", None)] += 1
        state.components._rec[cid] = replace(rec, composition=tuple(sorted(comp.items(), key=str)))
    with pytest.raises(mo.InvariantError, match=field):
        checker.check(state)


# --------------------------------------------------------------------------
# mobile-monomer ancestry
# --------------------------------------------------------------------------
def test_initially_mobile_monomer_gets_backbone_origin_through_the_unique_map():
    spec = S.ps_spec({4: 1}, mobile={"styrene": 2}, monomer=True)
    state = materialize(spec)
    mobiles = [c for c in state.components.ids() if state.origin_kinds(c) == {"mobile"}]
    assert len(mobiles) == 2
    for c in mobiles:
        flagged = sorted(state.origin[a].key for a in state.components.atoms(c) if state.origin[a].is_backbone)
        assert flagged == ["a", "b"]                      # the vinyl carbons, not the ring
        assert state.components.backbone_count(c) == 2
    records = {r.cid: r for r in mo.MeltObserver(spec, xs=0).components(state)}
    assert all(records[c].n == 1 and records[c].species == "styrene" for c in mobiles)
    # without a declared monomer_product there is no ancestry (and no refusal)
    plain = materialize(S.ps_spec({4: 1}, mobile={"styrene": 2}))
    assert not any(o.is_backbone for o in plain.origin.values() if o.is_mobile)
    # a non-PS monomer (ethylene: two isomorphisms, one ancestry) is mapped the same way
    pe = S.pe_spec({4: 1}, small={"ethylene": S.ETHYLENE}, mobile={"ethylene": 1}, monomer=S.ETHYLENE_GRAPH)
    pstate = materialize(pe)
    (eth,) = [c for c in pstate.components.ids() if pstate.origin_kinds(c) == {"mobile"}]
    assert pstate.components.backbone_count(eth) == 2


def test_ambiguous_monomer_ancestry_refuses_initial_mobile_monomer():
    with pytest.raises(AncestryError, match="ambiguous"):
        materialize(S.vinyl_butadiene_spec({4: 1}, mobile={"butadiene": 1}))
    spec = S.vinyl_butadiene_spec({4: 1}, mobile={"butadiene": 1})
    with pytest.raises(me.DeclarationError, match="ambiguous"):
        me.MeltEngine(spec, S.ToyCatalog(), 700.0, np.random.default_rng(0))
    # declared but not initially present: nothing to map, nothing refused
    state = materialize(S.vinyl_butadiene_spec({4: 1}))
    assert all(not o.is_mobile for o in state.origin.values())


# --------------------------------------------------------------------------
# release reporting
# --------------------------------------------------------------------------
def test_declared_mobile_species_in_no_loaded_template_is_refused():
    spec = S.ps_spec({6: 1}, declare_styrene=True)
    with pytest.raises(me.DeclarationError, match="styrene"):
        me.MeltEngine(spec, S.ps_catalog(), 700.0, np.random.default_rng(0))
    methyl = S.pe_spec({6: 1}, small={"methyl": S.METHYL})
    with pytest.raises(me.DeclarationError, match="methyl"):
        me.MeltEngine(methyl, S.ToyCatalog(names=PE_LIVE), 700.0, np.random.default_rng(0))
    engine = me.MeltEngine(methyl, S.ToyCatalog(names=PE_LIVE + ("methyl_recombination",)), 700.0,
                           np.random.default_rng(0))
    assert engine.report()["declared_mobile_participation"] == {"methyl": ["methyl_recombination"]}


def test_undeclared_light_products_are_listed_by_canonical_key():
    spec = S.pe_spec({6: 1})
    state = materialize(spec)
    end = S.end_carbons(state, 0)[0]
    nxt = next(b for b in state.atoms[end].edges if state.atoms[b].element == "C")
    homolyse(state, end, nxt)
    o = mo.MeltObserver(spec, xs=0).observe(state)
    (light,) = o["release"]["undeclared"]
    assert light["key"] == "[CH3]" and light["molecules"] == 1 and light["n"] == 0
    assert light["category"] == "light"
    assert o["tier1"]["undeclared_light"]["molecules"] == 1
    assert o["tier1"]["polymer_derived"]["molecules"] == 1
    declared = S.pe_spec({6: 1}, small={"methyl": S.METHYL})
    o2 = mo.MeltObserver(declared, xs=0).observe(state)
    assert o2["release"]["undeclared"] == [] and o2["tier1"]["species"]["methyl"]["molecules"] == 1
    total = sum(o["tier1"][k]["mass_g"] for k in ("polymer_derived", "undeclared_light"))
    assert rel(o["tier1"]["total_mass_g"], total) <= REL


# --------------------------------------------------------------------------
# non-PS polymer, export
# --------------------------------------------------------------------------
def test_non_ps_linear_polymer_through_the_same_observables():
    spec = S.pe_spec({10: 1, 20: 1, 30: 1})
    engine = me.MeltEngine(spec, S.ToyCatalog(names=PE_LIVE), 700.0, np.random.default_rng(5))
    observer = mo.MeltObserver(spec, xs=5)
    o0 = observer.observe(engine)
    assert o0["backbone_atoms_per_unit"] == 2
    assert o0["tier2"]["dp_moments_representative"] == [3, 60, 1400]
    for got, want in zip(o0["tier2"]["dp_moments_mol"], spec.macroscopic_moments_mol):
        assert rel(got, want) <= REL
    checker = mo.InvariantChecker(engine.state)
    engine.run(max_events=200, check=checker)
    assert checker.checks == 200 and checker.counts["split"] and checker.counts["merge"]
    o = observer.observe(engine)
    p = o["tier1"]["polymer_derived"]
    assert rel(p["mom_definition"]["mass_g"], p["mass_g"]) <= REL
    assert rel(p["mom_definition"]["amount_mol"], p["amount_mol"]) <= REL
    assert rel(o["tier1"]["total_mass_g"], o0["tier1"]["total_mass_g"]) <= REL
    for r in observer.components(engine.state):
        assert r.n == engine.state.components.backbone_count(r.cid) // 2


def test_observables_name_no_polymer():
    source = Path(mo.__file__).read_text()
    code = re.sub(r'"""[\s\S]*?"""', "", source)
    code = re.sub(r"#.*", "", code)
    assert not re.search(r"\b(PS|PE|PP|styrene|polystyrene|polyethylene|phenyl)\b", code, re.I)
    for name in ("PS_M0", "PS_C_R2", "rmgpy.kmc.met", "rmgpy.solver"):
        assert name not in source


def test_time_series_is_written_at_requested_times_without_perturbing_the_run(tmp_path):
    spec = S.pe_spec({10: 1, 20: 1})
    times = [0.0, 0.5, 2.0, 5.0]
    engine = me.MeltEngine(spec, S.ToyCatalog(names=PE_LIVE), 700.0, np.random.default_rng(42))
    path = tmp_path / "series.jsonl"
    out = mo.write_time_series(engine, times, path, xs=3)
    lines = [json.loads(x) for x in path.read_text().splitlines()]
    assert len(lines) == len(out) == 4
    assert [x["time_s"] for x in lines] == times
    reference = me.MeltEngine(spec, S.ToyCatalog(names=PE_LIVE), 700.0, np.random.default_rng(42))
    reference.run(until=5.0)
    assert reference.state.state_hash() == engine.state.state_hash()
    assert [(e.time, e.template) for e in reference.events] == [(e.time, e.template) for e in engine.events]
    last = lines[-1]
    assert last["counts"]["events"] == engine.events_fired > 0
    assert last["counts"]["family"] == dict(engine.family_counts)
    for part in ("leak_num", "leak_den", "non_thinned_integral", "phi_sum", "thinned_proposals"):
        assert part in last["leak_estimate"]
    assert last["invariant_checks"] == engine.events_fired
    for key in ("tier1", "tier2", "tier3", "release"):
        assert key in last
    mom = mo.mom_reported_quantities(last)
    for key in ("mu0_mol", "mu1_mol", "mu2_mol", "explicit_dp_mol", "end_inventory_mol",
                "polymer_condensed_mass_g", "Mn_g_mol", "Mw_g_mol", "polymer_conversion", "species_mol",
                "polymer_derived_molecules_mol"):
        assert key in mom
    # a MoM run that reports exactly the kMC state agrees on every compared row
    rows = mo.side_by_side(last, mom)
    assert {r["tier"] for r in rows} == {1, 3}
    assert all(r["agreement"] in ("equal", "not_a_comparison") for r in rows), rows


# --------------------------------------------------------------------------
# rework 1
# --------------------------------------------------------------------------
def _mass_law_replay(spec, o):
    """MoM's shared mass law on the exported quantities, written out here:
    the tail pool (repeat x mu1 + end-group composition x E) plus every
    explicit DP-n amount at the capped DP-n molar mass."""
    from rmgpy.polymer_mass import EvolvedPoolMass, InventoryTerm
    q = mo.mom_reported_quantities(o)
    pool = EvolvedPoolMass("tail", spec.repeat_unit.composition, (q["mu0_mol"], q["mu1_mol"], q["mu2_mol"]),
                           end_inventory=(InventoryTerm("ends", spec.end_composition, q["end_inventory_mol"],
                                                        True),))
    a = spec.repeat_unit.composition.molar_mass_g_mol()
    e = spec.end_composition.molar_mass_g_mol()
    return pool.mass_g + sum(amount * (a * int(n) + e) for n, amount in q["explicit_dp_mol"].items()), q


def test_end_inventory_is_physical_and_bridge_scission_keeps_the_mass_law_exact():
    spec = S.pe_spec({6: 1})
    state = materialize(spec)
    path = chain_carbons_in_order(state, 0)
    u = spec.initialization().scale_factor / constants.Na
    homolyse(state, path[5], path[6])                       # two DP3 radicals
    o = mo.MeltObserver(spec, xs=0).observe(state)
    replay, q = _mass_law_replay(spec, o)
    kmc = o["tier1"]["polymer_derived"]["mass_g"]
    assert q["mu0_mol"] == pytest.approx(2 * u, rel=REL)
    assert q["end_inventory_mol"] == pytest.approx(u, rel=REL)  # one end pair: no cap atom appeared
    assert rel(replay, kmc) <= REL, (replay, kmc)
    assert rel(q["polymer_condensed_mass_g"], kmc) <= REL
    law = o["mom_mass_law"]
    assert law["difference_g"] == 0.0 and law["bound_g"] == 0.0
    recombine(state, path[5], path[6])
    o = mo.MeltObserver(spec, xs=0).observe(state)
    replay, q = _mass_law_replay(spec, o)
    assert q["end_inventory_mol"] == pytest.approx(u, rel=REL)
    assert rel(replay, o["tier1"]["polymer_derived"]["mass_g"]) <= REL


def test_mass_law_gap_of_radical_explicit_chains_is_covered_by_the_bound():
    """The two DP3 radicals in the explicit ladder (xs = 3): MoM values each at
    the capped DP3 mass, so the gap is one cap H per radical end."""
    spec = S.pe_spec({6: 1})
    state = materialize(spec)
    path = chain_carbons_in_order(state, 0)
    u = spec.initialization().scale_factor / constants.Na
    homolyse(state, path[5], path[6])
    o = mo.MeltObserver(spec, xs=3).observe(state)
    replay, q = _mass_law_replay(spec, o)
    kmc = o["tier1"]["polymer_derived"]["mass_g"]
    h = Composition.from_mapping({"H": 1}).molar_mass_g_mol()
    law = o["mom_mass_law"]
    assert law["difference_g"] == pytest.approx(kmc - replay, rel=1e-9)
    assert law["difference_g"] == pytest.approx(-2 * u * h, rel=1e-9)
    assert abs(kmc - replay) <= law["bound_g"]


def test_ps_chain_end_homolysis_mass_law_gap_is_covered_by_the_bound():
    spec = S.ps_spec({10: 1})
    state = materialize(spec)
    path = chain_carbons_in_order(state, 0)
    u = spec.initialization().scale_factor / constants.Na
    homolyse(state, path[2], path[3])          # 2-phenylpropyl (B = 3, n = 1) + the rest (B = 17, n = 8)
    o = mo.MeltObserver(spec, xs=0).observe(state)
    replay, q = _mass_law_replay(spec, o)
    kmc = o["tier1"]["polymer_derived"]["mass_g"]
    a = spec.repeat_unit.composition.molar_mass_g_mol()
    assert q["end_inventory_mol"] == pytest.approx(u, rel=REL)
    assert q["mu1_mol"] == pytest.approx(9 * u, rel=REL)
    assert kmc - replay == pytest.approx(u * a, rel=1e-9)         # one repeat unit split over two ends
    assert o["mom_mass_law"]["difference_g"] == pytest.approx(kmc - replay, rel=1e-9)
    assert abs(kmc - replay) <= o["mom_mass_law"]["bound_g"]


def _add_monomer_to_chain_end(state, mobile_cid):
    """Strip an H from the chain's backbone end and add the monomer's
    backbone C=C to it (the monomer becomes the chain's next unit)."""
    path = chain_carbons_in_order(state, 0)
    S.strip_h(state, path[0])
    flagged = [x for x in state.components.atoms(mobile_cid) if state.origin[x].is_backbone]
    assert len(flagged) == 2
    a, b = flagged
    if sum(1 for y in state.atoms[a].edges if state.atoms[y].element == "H") < 2:
        a, b = b, a
    state.apply([SetBond(path[0], a, 0.0, 1.0), SetBond(a, b, 2.0, 1.0),
                 SetRadical(path[0], 1, 0), SetRadical(b, 0, 1)])


@pytest.mark.parametrize("case", ["explicit_h_species", "explicit_h_both", "kekule_species", "kekule_both"])
def test_monomer_ancestry_is_matched_after_normalization(case):
    if case.startswith("explicit_h"):
        monomer = S.ETHYLENE_EXPLICIT_H if case == "explicit_h_both" else S.ETHYLENE_GRAPH
        spec = S.pe_spec({4: 1}, small={"ethylene": {"graph": S.ETHYLENE_EXPLICIT_H, "disposition": "mobile"}},
                         mobile={"ethylene": 1}, monomer=monomer)
    else:
        monomer = S.styrene_graph(kekule=(case == "kekule_both"))
        spec = S.ps_spec_declared({4: 1}, S.styrene_graph(kekule=True), monomer, mobile={"styrene": 1})
    state = materialize(spec)
    (mob,) = [c for c in state.components.ids() if state.origin_kinds(c) == {"mobile"}]
    assert state.components.backbone_count(mob) == 2
    _add_monomer_to_chain_end(state, mob)
    o = mo.MeltObserver(spec, xs=0).observe(state)
    assert o["tier2"]["dp_moments_representative"] == [1, 5, 25]


def test_declared_small_species_isomorphic_to_an_initial_chain_is_refused():
    spec = S.pe_spec({2: 1}, small={"butane": S.BUTANE})
    with pytest.raises(me.DeclarationError, match="butane"):
        mo.MeltObserver(spec, xs=0)
    with pytest.raises(me.DeclarationError, match="butane"):
        me.MeltEngine(spec, S.ToyCatalog(), 700.0, np.random.default_rng(0))
    ok = S.pe_spec({3: 1}, small={"butane": S.BUTANE})
    o = mo.MeltObserver(ok, xs=0).observe(materialize(ok))
    for got, want in zip(o["tier2"]["dp_moments_mol"], ok.macroscopic_moments_mol):
        assert rel(got, want) <= REL
    assert o["tier1"]["species"]["butane"]["molecules"] == 0


def test_a_fragment_matching_a_declared_species_is_that_species():
    """Run-time rule: a component isomorphic (under resonance) to a declared
    species is that species whatever its ancestry, so a polymer-derived
    fragment (here n_c = 1) that matches one leaves the polymer class."""
    spec = S.pe_spec({6: 1}, small={"ethyl": S.ETHYL})
    state = materialize(spec)
    path = chain_carbons_in_order(state, 0)
    homolyse(state, path[1], path[2])          # ethyl (B = 2, n = 1) + C10H21 radical
    o = mo.MeltObserver(spec, xs=0).observe(state)
    u = spec.initialization().scale_factor / constants.Na
    assert o["tier1"]["species"]["ethyl"]["molecules"] == 1
    assert o["tier1"]["polymer_derived"]["molecules"] == 1
    assert mo.mom_reported_quantities(o)["species_mol"]["ethyl"] == pytest.approx(u, rel=REL)
    assert "declared species" in mo.SPECIES_RULE


def test_polymer_light_boundary_is_n_c_equal_one():
    ps = S.ps_spec({10: 1})
    state = materialize(ps)
    path = chain_carbons_in_order(state, 0)
    homolyse(state, path[2], path[3])          # 2-phenylpropyl: n_c = 1 -> polymer-derived
    o = mo.MeltObserver(ps, xs=0).observe(state)
    assert o["tier1"]["polymer_derived"]["molecules"] == 2 and o["tier1"]["undeclared_light"]["molecules"] == 0
    pe = S.pe_spec({6: 1})
    state = materialize(pe)
    path = chain_carbons_in_order(state, 0)
    homolyse(state, path[0], path[1])          # methyl: n_c = 0 -> light
    o = mo.MeltObserver(pe, xs=0).observe(state)
    assert o["tier1"]["polymer_derived"]["molecules"] == 1 and o["tier1"]["undeclared_light"]["molecules"] == 1


def test_side_by_side_differences_come_from_the_supplied_mom_values():
    spec = S.pe_spec({6: 1}, small={"methyl": S.METHYL})
    o = mo.MeltObserver(spec, xs=0).observe(materialize(spec))
    assert "accounting identity" in o["tier1"]["mom_definition_note"]
    assert "accounting identity" in mo.__doc__
    rows = {r["quantity"]: r for r in mo.side_by_side(o, {"polymer_condensed_mass_g": 1.0e5,
                                                          "species_mol": {"methyl": 1.0}})}
    mass = rows["polymer_condensed_mass_g"]
    assert mass["kmc"] == pytest.approx(1000.0, rel=REL) and mass["mom"] == 1.0e5
    assert mass["difference"] == pytest.approx(1000.0 - 1.0e5, rel=REL) and mass["agreement"] == "differs"
    methyl = rows["species_mol.methyl"]
    assert methyl["kmc"] == 0.0 and methyl["mom"] == 1.0 and methyl["agreement"] == "differs"
    assert rows["mu0_mol"]["agreement"] == "missing_in_mom"
    for key in ("mu0_mol", "mu1_mol", "mu2_mol", "end_inventory_mol", "polymer_condensed_mass_g",
                "convention_gap_g", "Mn_g_mol", "Mw_g_mol", "polymer_conversion",
                "polymer_derived_molecules_mol", "species_mol.methyl", "species_mol.H", "species_mol.H2"):
        assert key in rows, key
    # a MoM quantity the kMC run never produced still gets a row
    rows = {r["quantity"]: r for r in mo.side_by_side(o, {"explicit_dp_mol": {"7": 0.5}})}
    assert rows["explicit_dp_mol.7"]["kmc"] == 0.0 and rows["explicit_dp_mol.7"]["agreement"] == "differs"
    # no MoM values: no difference is claimed
    plain = mo.side_by_side(o)
    assert all(r["difference"] is None and r["agreement"] == "no_mom_values"
               for r in plain if r["quantity"] != "convention_gap_g")
    assert [r["agreement"] for r in plain if r["quantity"] == "convention_gap_g"] == ["not_a_comparison"]


@pytest.mark.parametrize("label", ["H", "H2", "monomer_product"])
def test_a_declared_label_colliding_with_a_builtin_is_refused(label):
    spec = S.pe_spec({6: 1}, small={label: S.METHYL})
    with pytest.raises(me.DeclarationError, match="built-in"):
        mo.MeltObserver(spec, xs=0)
    with pytest.raises(me.DeclarationError, match="built-in"):
        me.MeltEngine(spec, S.ToyCatalog(builder=lambda T: S.pe_templates(T)), 700.0, np.random.default_rng(0))


# --------------------------------------------------------------------------
# rework 2
# --------------------------------------------------------------------------
def _remove_br(state, carbon):
    br = next(x for x in state.atoms[carbon].edges if state.atoms[x].element == "Br")
    homolyse(state, carbon, br)
    return br


M_BR = Composition.from_mapping({"Br": 1}).molar_mass_g_mol()
PVB_CHANNEL = [S.sgh_channel("aliphatic_C-Br", "aliphatic", 1.0)]


def _channel(o, label):
    (row,) = [r for r in o["side_group_channels"] if r["label"] == label]
    return row


def _observe(spec, xs, channels, edit):
    """Observer + release tracker on the declared melt; ``edit(state)`` is
    applied, the tracker updated from the state (rescan), and one
    observation returned with the state and the per-molecule unit."""
    state = materialize(spec)
    observer = mo.MeltObserver(spec, xs=xs, side_group_homolysis=channels)
    tracker = observer.release_tracker(state)
    edit(state, tracker)
    tracker.update(state)
    u = spec.initialization().scale_factor / constants.Na
    return observer.observe(state, releases=tracker), state, u


def test_a_true_release_on_a_tail_chain_is_counted_once_from_the_event():
    spec = S.pvb_spec({4: 1})
    state = materialize(spec)
    observer = mo.MeltObserver(spec, xs=0, side_group_homolysis=PVB_CHANNEL)
    tracker = observer.release_tracker(state)
    path = chain_carbons_in_order(state, 0)
    cid = state.components.component_of(path[0])
    br = _remove_br(state, path[1])
    products = tuple(sorted({cid, state.components.component_of(br)}))
    tracker.update(state, me.Event(0.0, "hand", "hand", "uni", (cid,), products, state.revision))
    tracker.update(state, me.Event(0.0, "hand", "hand", "uni", (cid,), products, state.revision))  # replay
    o = observer.observe(state, releases=tracker)
    u = spec.initialization().scale_factor / constants.Na
    row = _channel(o, "aliphatic_C-Br")
    assert row["released_count"] == 1 and rel(row["released_mol"], u) <= REL
    assert rel(row["released_mass_g"], u * M_BR) <= REL and rel(row["free_now_mol"], u) <= REL
    assert rel(row["inventory_tail_mol"], 3 * u) <= REL and row["inventory_explicit_mol"] == 0.0
    law = o["mom_mass_law"]
    assert "side_group_depletion" not in law
    assert rel(law["value_g"], 1000.0) <= REL                                # the no-depletion law
    kmc = o["tier1"]["polymer_derived"]["mass_g"]
    assert rel(kmc, 1000.0 - u * M_BR) <= REL
    assert rel(law["difference_g"], -u * M_BR) <= 1e-9 and abs(law["difference_g"]) <= law["bound_g"]
    assert law["undeclared_loss"]["mass_g"] == 0.0                           # Br left through the channel
    q = mo.mom_reported_quantities(o)
    assert rel(q["polymer_condensed_mass_g"], kmc) <= REL                    # the physical condensed mass
    assert rel(q["released_mol"]["aliphatic_C-Br"], u) <= REL
    assert "side_group_Z_mol" not in q


def test_a_true_release_on_an_explicit_chain_is_counted_once():
    def edit(state, _t):
        _remove_br(state, chain_carbons_in_order(state, 0)[3])
    o, _s, u = _observe(S.pvb_spec({4: 1}), 5, PVB_CHANNEL, edit)
    row = _channel(o, "aliphatic_C-Br")
    assert row["released_count"] == 1
    assert row["inventory_tail_mol"] == 0.0 and rel(row["inventory_explicit_mol"], 3 * u) <= REL


def test_a_partition_crossing_books_no_release():
    """PVB DP6 cut into two DP3 at xs = 3: the chains leave the tail, no Br leaves the melt."""
    def edit(state, _t):
        path = chain_carbons_in_order(state, 0)
        homolyse(state, path[5], path[6])
    o, _s, u = _observe(S.pvb_spec({6: 1}), 3, PVB_CHANNEL, edit)
    row = _channel(o, "aliphatic_C-Br")
    assert row["released_count"] == 0 and row["released_mol"] == 0.0
    assert row["inventory_tail_mol"] == 0.0 and rel(row["inventory_explicit_mol"], 6 * u) <= REL
    q = mo.mom_reported_quantities(o)
    assert rel(q["polymer_condensed_mass_g"], 1000.0) <= REL


def test_growth_then_a_real_loss_books_exactly_one_release():
    spec = S.pvb_spec_with_monomer({4: 1}, {"vinyl_bromide": 1})

    def edit(state, tracker):
        (mob,) = [c for c in state.components.ids() if state.origin_kinds(c) == {"mobile"}]
        _add_monomer_to_chain_end(state, mob)
        tracker.update(state)
        _remove_br(state, chain_carbons_in_order(state, 0)[3])
    o, _s, u = _observe(spec, 0, PVB_CHANNEL, edit)
    row = _channel(o, "aliphatic_C-Br")
    assert row["released_count"] == 1
    assert rel(row["inventory_tail_mol"], 4 * u) <= REL                       # 4 + 1 grown - 1 released


def test_equivalent_br_swap_books_nothing():
    def edit(state, _t):
        path = chain_carbons_in_order(state, 0)
        ca, cb = path[1], path[3]
        ba = next(x for x in state.atoms[ca].edges if state.atoms[x].element == "Br")
        bb = next(x for x in state.atoms[cb].edges if state.atoms[x].element == "Br")
        state.apply([SetBond(ca, ba, 1.0, 0.0), SetBond(cb, bb, 1.0, 0.0),
                     SetBond(ca, bb, 0.0, 1.0), SetBond(cb, ba, 0.0, 1.0)])
    o, _s, u = _observe(S.pvb_spec({4: 1}), 0, PVB_CHANNEL, edit)
    row = _channel(o, "aliphatic_C-Br")
    assert row["released_count"] == 0 and rel(row["inventory_tail_mol"], 4 * u) <= REL


def test_identity_h_transfer_books_nothing():
    """A C6H13 primary radical takes an H from its own CH3 end: the graph is
    unchanged (C6H13 primary radical), no H leaves the melt."""
    channel = [S.sgh_channel("aliphatic_C-H", "aliphatic", 4.0, gas="[H]")]
    counts = {}

    def edit(state, tracker):
        path = chain_carbons_in_order(state, 0)
        homolyse(state, path[5], path[6])
        tracker.update(state)
        counts["before"] = _inventory(state, tracker)
        head, rad = path[0], path[5]
        h = min(x for x in state.atoms[head].edges if state.atoms[x].element == "H")
        state.apply([SetBond(head, h, 1.0, 0.0), SetBond(rad, h, 0.0, 1.0),
                     SetRadical(rad, 1, 0), SetRadical(head, 0, 1)])
    o, _s, u = _observe(S.pe_spec({6: 1}), 0, channel, edit)
    row = _channel(o, "aliphatic_C-H")
    assert row["released_count"] == 0
    assert rel(row["inventory_tail_mol"], 26 * u) <= REL and counts["before"] == 26


def _inventory(state, tracker):
    return sum(tracker.observer.site_inventory(state).get(0, {}).values())


def test_c_h_channel_counts_hydrogen_release_and_inventory_by_current_graph():
    channel = [S.sgh_channel("aliphatic_C-H", "aliphatic", 4.0, gas="[H]")]

    def edit(state, _t):
        S.strip_h(state, chain_carbons_in_order(state, 0)[5])
    o, _s, u = _observe(S.pe_spec({6: 1}), 0, channel, edit)
    row = _channel(o, "aliphatic_C-H")
    assert row["released_count"] == 1 and rel(row["inventory_tail_mol"], 25 * u) <= REL


def test_aryl_channel_inventory_is_the_br_on_aromatic_carbons():
    def edit(state, _t):
        br = next(x for x in state.atoms if state.atoms[x].element == "Br")
        homolyse(state, next(iter(state.atoms[br].edges)), br)
    o, _s, u = _observe(S.para_bromo_ps_spec({4: 1}), 0, [S.sgh_channel("aryl_C-Br", "aryl", 1.0)], edit)
    row = _channel(o, "aryl_C-Br")
    assert row["released_count"] == 1 and rel(row["inventory_tail_mol"], 3 * u) <= REL
    assert rel(o["mom_mass_law"]["value_g"], 1250.0) <= REL
    assert mo.mom_reported_quantities(o)["polymer_condensed_mass_g"] == pytest.approx(1113.960, abs=5e-4)


def test_no_declared_channel_reports_the_loss_as_undeclared():
    def edit(state, _t):
        _remove_br(state, chain_carbons_in_order(state, 0)[1])
    o, _s, u = _observe(S.pvb_spec({4: 1}), 0, None, edit)
    law = o["mom_mass_law"]
    assert o["side_group_channels"] == [] and rel(law["value_g"], 1000.0) <= REL
    assert rel(law["undeclared_loss"]["mass_g"], u * M_BR) <= REL and law["undeclared_loss"]["atoms"] == {"Br": 1}


def test_pendant_phenyl_loss_has_no_mom_channel_and_is_an_undeclared_loss():
    def edit(state, _t):
        c2 = chain_carbons_in_order(state, 0)[9]
        ring = next(x for x in state.atoms[c2].edges
                    if state.atoms[x].element == "C" and not state.origin[x].is_backbone)
        homolyse(state, c2, ring)
    o, _s, u = _observe(S.ps_spec({10: 1}), 0, None, edit)
    m_ph = Composition.from_mapping({"C": 6, "H": 5}).molar_mass_g_mol()
    law = o["mom_mass_law"]
    assert o["side_group_channels"] == [] and rel(law["value_g"], 1250.0) <= REL
    assert rel(law["undeclared_loss"]["mass_g"], u * m_ph) <= REL
    assert law["undeclared_loss"]["atoms"] == {"C": 6, "H": 5}


def test_declared_channels_are_separate_and_compared_per_label():
    channels = [S.sgh_channel("aliphatic_C-Br", "aliphatic", 1.0), S.sgh_channel("benzylic_C-Br", "benzylic", 1.0)]

    def edit(state, _t):
        _remove_br(state, chain_carbons_in_order(state, 0)[2])     # the aliphatic Br of unit 2
    o, _s, u = _observe(S.two_bromine_spec({4: 1}), 0, channels, edit)
    ali, ben = _channel(o, "aliphatic_C-Br"), _channel(o, "benzylic_C-Br")
    assert ali["released_count"] == 1 and ben["released_count"] == 0
    assert rel(ali["inventory_tail_mol"], 3 * u) <= REL and rel(ben["inventory_tail_mol"], 4 * u) <= REL
    mom = mo.mom_reported_quantities(o)
    mom["released_mol"] = {"aliphatic_C-Br": u, "benzylic_C-Br": u}
    rows = {r["quantity"]: r for r in mo.side_by_side(o, mom)}
    assert rows["released_mol.aliphatic_C-Br"]["agreement"] == "equal"
    assert rows["released_mol.benzylic_C-Br"]["agreement"] == "differs"


def test_untracked_observation_reports_no_release():
    spec = S.pvb_spec({4: 1})
    o = mo.MeltObserver(spec, xs=0, side_group_homolysis=PVB_CHANNEL).observe(materialize(spec))
    row = _channel(o, "aliphatic_C-Br")
    assert row["released_tracked"] is False and row["released_count"] is None
    rows = {r["quantity"]: r for r in mo.side_by_side(o, {"released_mol": {"aliphatic_C-Br": 0.0}})}
    assert rows["released_mol.aliphatic_C-Br"]["agreement"] == "missing_in_kmc"


def test_declared_key_colliding_with_a_generated_hydrogen_key_is_refused():
    """An explicit H declared with key ``c1.H0`` (bonded to c2) while c1
    carries implicit hydrogens: the generated key would collide."""
    mapping = S.pe_mapping({4: 1})
    repeat = mapping["repeat_unit"]
    repeat["atoms"] = [S._atom("c1", hydrogens=2), S._atom("c2", hydrogens=1), {"key": "c1.H0", "element": "H"}]
    repeat["bonds"].append({"atoms": ["c2", "c1.H0"], "order": 1})
    from rmgpy.polymer_spec import LinearHomopolymerSpec
    spec = LinearHomopolymerSpec.from_mapping(mapping)
    with pytest.raises(MeltStateError, match="collide"):
        materialize(spec)
    with pytest.raises((MeltStateError, me.DeclarationError), match="collide"):
        mo.MeltObserver(spec, xs=0)


@pytest.mark.parametrize("bad, match", [
    ([S.sgh_channel("x", "aliphatic", 1.0, gas="[CH3]")], "monoatomic"),
    ([S.sgh_channel("x", "aliphatic", 2.0)], "sites_per_unit"),
    ([S.sgh_channel("x", "aryl", 1.0)], "matches NO"),
    ([S.sgh_channel("x", "aliphatic", 1.0), S.sgh_channel("y", "aliphatic", 1.0)], "SAME"),
    ([S.sgh_channel("x", "aliphatic", 1.0), S.sgh_channel("x", "aliphatic", 1.0)], "duplicate"),
    ([dict(S.sgh_channel("x", "aliphatic", 1.0), extra=1)], "exactly"),
    ([], "empty"),
])
def test_only_channels_mom_can_represent_are_accepted(bad, match):
    with pytest.raises(mo.ObservablesError, match=match):
        mo.MeltObserver(S.pvb_spec({4: 1}), xs=0, side_group_homolysis=bad)


def test_explicit_ladder_term_of_the_mass_law_is_the_capped_dp_mass():
    spec = S.pe_spec({6: 1, 10: 1})
    state = materialize(spec)
    path = chain_carbons_in_order(state, 0)       # the DP6 chain
    assert len(path) == 12
    homolyse(state, path[5], path[6])             # two DP3 radicals: explicit at xs = 3
    o = mo.MeltObserver(spec, xs=3).observe(state)
    u = spec.initialization().scale_factor / constants.Na
    a = spec.repeat_unit.composition.molar_mass_g_mol()
    e = spec.end_composition.molar_mass_g_mol()
    law = o["mom_mass_law"]
    assert rel(law["explicit_capped_mass_g"], 2 * u * (3 * a + e)) <= REL
    hand = a * 10 * u + e * u + 2 * u * (3 * a + e)
    assert rel(law["value_g"], hand) <= REL


def test_mom_numbers_are_compared_only_in_the_same_convention():
    """A MoM run that released 49.9 g of styrene the kMC melt did not: every
    mass row differs, and the convention gap is never an agreement."""
    spec = S.ps_spec({10: 1}, declare_styrene=True)
    state = materialize(spec)
    path = chain_carbons_in_order(state, 0)
    homolyse(state, path[2], path[3])             # convention gap of one repeat unit (~125 g)
    o = mo.MeltObserver(spec, xs=0).observe(state)
    assert o["mom_mass_law"]["bound_g"] > 49.9
    mom = mo.mom_reported_quantities(o)
    m_sty = Composition.from_mapping({"C": 8, "H": 8}).molar_mass_g_mol()
    mom["polymer_condensed_mass_g"] -= 49.9
    mom["polymer_conversion"] += 49.9 / o["initial_polymer_mass_g"]
    mom["species_mol"]["styrene"] += 49.9 / m_sty
    rows = {r["quantity"]: r for r in mo.side_by_side(o, mom)}
    for key in ("polymer_condensed_mass_g", "polymer_conversion", "species_mol.styrene"):
        assert rows[key]["agreement"] == "differs", rows[key]
    gap = rows["convention_gap_g"]
    assert gap["agreement"] == "not_a_comparison" and gap["mom"] is None and gap["difference"] is None
    assert abs(gap["kmc"]) <= gap["bound"]
    assert not any(r["agreement"] == "within_bound" for r in rows.values())
    assert "polymer_derived_mass_exact_g" not in rows


@pytest.mark.parametrize("relative, expected", [(1.0e-6, "differs"), (1.0e-14, "equal")])
def test_equal_means_a_relative_difference_of_at_most_1e_12(relative, expected):
    spec = S.pe_spec({6: 1})
    o = mo.MeltObserver(spec, xs=0).observe(materialize(spec))
    mom = mo.mom_reported_quantities(o)
    mom["polymer_condensed_mass_g"] *= 1.0 + relative
    rows = {r["quantity"]: r for r in mo.side_by_side(o, mom)}
    assert rows["polymer_condensed_mass_g"]["agreement"] == expected


def test_channel_rules_mirror_the_solver():
    """Parity with the MoM solver's own site selection and channel validation
    (read here, never imported by the observables module)."""
    solver = pytest.importorskip("rmgpy.solver.polymer")
    specs = [S.pvb_spec({4: 1}), S.para_bromo_ps_spec({4: 1}), S.two_bromine_spec({4: 1}), S.pe_spec({4: 1})]
    for spec in specs:
        repeat, _ends = spec.legacy_molecules()
        for element in ("Br", "H"):
            for selector in ("aryl", "benzylic", "aliphatic"):
                assert (mo.side_group_site_indices(repeat, element, selector)
                        == solver.side_group_site_atom_indices(repeat, element, selector))
    cases = [[S.sgh_channel("x", "aliphatic", 1.0, gas="[CH3]")], [S.sgh_channel("x", "aliphatic", 2.0)],
             [S.sgh_channel("x", "aryl", 1.0)], [S.sgh_channel("x", "aliphatic", 1.0)],
             [S.sgh_channel("x", "aliphatic", 1.0), S.sgh_channel("y", "aliphatic", 1.0)]]
    spec = S.pvb_spec({4: 1})
    for channels in cases:
        repeat, _ends = spec.legacy_molecules()
        try:
            solver.validate_side_group_homolysis("pool", channels, monomer=repeat)
            ok_solver = True
        except ValueError:
            ok_solver = False
        try:
            mo.validate_side_group_channels(spec, channels)
            ok_mine = True
        except mo.ObservablesError:
            ok_mine = False
        assert ok_solver == ok_mine, channels


# --------------------------------------------------------------------------
# I-158 rework 1
# --------------------------------------------------------------------------
H_CHANNEL = [S.sgh_channel("aliphatic_C-H", "aliphatic", 4.0, gas="[H]")]


def _h_loss_engine(seed=3, **kw):
    spec = S.pe_spec({6: 10})                         # 20 CH3 ends: up to 20 H-loss events
    return me.MeltEngine(spec, S.ToyCatalog(builder=lambda T: S.h_loss_templates(T, **kw)), 700.0,
                         np.random.default_rng(seed))


def _balance_holds(observer, state, tracker):
    for element, b in observer.x_balance(state, tracker).items():
        parts = (b["condensed_at_declared_sites_count"] + b["condensed_elsewhere_count"] + b["free_count"]
                 + b["other_count"])
        assert b["initial_count"] == parts, (element, b)


def test_appended_export_keeps_cumulative_releases(tmp_path):
    engine = _h_loss_engine()
    path = tmp_path / "series.jsonl"
    first = mo.write_time_series(engine, [10.0], path, xs=0, side_group_homolysis=H_CHANNEL)[-1]
    second = mo.write_time_series(engine, [10.0], path, xs=0, side_group_homolysis=H_CHANNEL, append=True)[-1]
    a, b = _channel(first, "aliphatic_C-H"), _channel(second, "aliphatic_C-H")
    assert engine.events_fired >= 1 and a["released_count"] == engine.events_fired
    assert b["released_count"] == a["released_count"] and a["released_tracked"] and b["released_tracked"]
    third = mo.write_time_series(engine, [10.0, 20.0], path, xs=0, side_group_homolysis=H_CHANNEL, append=True)[-1]
    assert _channel(third, "aliphatic_C-H")["released_count"] == engine.events_fired


def test_an_export_that_cannot_prove_continuity_reports_untracked(tmp_path):
    engine = _h_loss_engine()
    engine.run(until=10.0)
    assert engine.events_fired >= 1
    obs = mo.write_time_series(engine, [10.0], tmp_path / "s.jsonl", xs=0, side_group_homolysis=H_CHANNEL)[-1]
    row = _channel(obs, "aliphatic_C-H")
    assert row["released_tracked"] is False


def test_tracker_recovers_a_release_committed_before_a_post_commit_failure():
    engine = _h_loss_engine()
    observer = mo.MeltObserver(engine.spec, xs=0, side_group_homolysis=H_CHANNEL)
    tracker = observer.release_tracker(engine)
    original = engine._add
    calls = {"n": 0}

    def once(*args, **kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            raise RuntimeError("injected in _add")
        return original(*args, **kwargs)

    engine._add = once
    with pytest.raises(RuntimeError, match="injected"):
        engine.run(max_events=1, check=tracker)
    del engine._add
    engine.check_consistency()
    assert engine.events_fired == 1 and engine.events[0].template == "h_loss"
    row = _channel(observer.observe(engine, releases=tracker), "aliphatic_C-H")
    assert row["released_tracked"] is True and row["released_count"] == 1
    # two events the tracker never saw: continuity cannot be proven
    engine.run(max_events=2)
    row = _channel(observer.observe(engine, releases=tracker), "aliphatic_C-H")
    assert row["released_tracked"] is False


def test_an_unattributed_polymer_release_is_an_undeclared_loss():
    def edit(state, _t):
        _remove_br(state, chain_carbons_in_order(state, 0)[3])     # the benzylic Br of unit 2
    o, _s, u = _observe(S.two_bromine_spec({4: 1}), 0, [S.sgh_channel("aliphatic_C-Br", "aliphatic", 1.0)], edit)
    assert _channel(o, "aliphatic_C-Br")["released_count"] == 0
    assert o["side_group_unattributed_releases"] == {"polymer": 1, "other": 0}
    loss = o["mom_mass_law"]["undeclared_loss"]
    assert loss["atoms"] == {"Br": 1} and rel(loss["mass_g"], u * M_BR) <= REL


def test_mobile_species_cleavage_is_not_a_polymer_release():
    spec = S.pvb_with_methyl_bromide({4: 1})

    def edit(state, _t):
        (mob,) = [c for c in state.components.ids() if state.origin_kinds(c) == {"mobile"}]
        c = next(x for x in state.components.atoms(mob) if state.atoms[x].element == "C")
        _remove_br(state, c)
    o, _s, u = _observe(spec, 0, PVB_CHANNEL, edit)
    row = _channel(o, "aliphatic_C-Br")
    assert row["released_count"] == 0 and row["released_from_other_count"] == 1
    assert rel(row["inventory_tail_mol"], 4 * u) <= REL
    assert rel(o["tier1"]["polymer_derived"]["mass_g"], 1000.0) <= REL


@pytest.mark.parametrize("kekule", ["ring_b_alternating", "ring_b_not_alternating"])
def test_site_classification_is_resonance_invariant(kekule):
    spec = S.bromonaphthyl_spec({3: 1}, kekule)
    channel = [S.sgh_channel("aryl_C-Br", "aryl", 1.0)]
    state = materialize(spec)
    observer = mo.MeltObserver(spec, xs=0, side_group_homolysis=channel)
    assert observer.site_inventory(state)[0] == {"tail": 3, "explicit": 0}

    def edit(state, _t):
        br = next(x for x in state.atoms if state.atoms[x].element == "Br")
        homolyse(state, next(iter(state.atoms[br].edges)), br)
    o, _s, u = _observe(spec, 0, channel, edit)
    row = _channel(o, "aryl_C-Br")
    assert row["released_count"] == 1 and rel(row["inventory_tail_mol"], 2 * u) <= REL
    assert o["side_group_unattributed_releases"] == {"polymer": 0, "other": 0}


def test_a_charged_gas_product_is_matched_exactly():
    spec = S.oxy_polymer_spec({4: 1}, -1)
    channel = [S.sgh_channel("aliphatic_C-O", "aliphatic", 1.0, gas="[O-]")]

    def edit(state, _t):
        c2 = chain_carbons_in_order(state, 0)[1]
        o_atom = next(x for x in state.atoms[c2].edges if state.atoms[x].element == "O")
        homolyse(state, c2, o_atom)
    o, _s, u = _observe(spec, 0, channel, edit)
    row = _channel(o, "aliphatic_C-O")
    assert row["released_count"] == 1 and rel(row["free_now_mol"], u) <= REL
    assert o["mom_mass_law"]["undeclared_loss"]["mass_g"] == 0.0


def test_a_gas_product_that_homolysis_cannot_produce_is_refused():
    spec = S.oxy_polymer_spec({4: 1}, +1)              # O+ site: homolysis gives O+., not [O-]
    with pytest.raises(mo.ObservablesError, match="cannot be matched"):
        mo.MeltObserver(spec, xs=0, side_group_homolysis=[S.sgh_channel("x", "aliphatic", 1.0, gas="[O-]")])


def test_release_count_is_gross_and_recaptures_are_reported():
    spec = S.pvb_spec({4: 1})
    state = materialize(spec)
    observer = mo.MeltObserver(spec, xs=0, side_group_homolysis=PVB_CHANNEL)
    tracker = observer.release_tracker(state)
    path = chain_carbons_in_order(state, 0)
    _balance_holds(observer, state, tracker)
    br = _remove_br(state, path[1])
    tracker.update(state)
    _balance_holds(observer, state, tracker)
    recombine(state, path[1], br)                       # the free Br re-adds at its site
    tracker.update(state)
    _balance_holds(observer, state, tracker)
    _remove_br(state, path[1])
    tracker.update(state)
    _balance_holds(observer, state, tracker)
    row = _channel(observer.observe(state, releases=tracker), "aliphatic_C-Br")
    assert row["released_count"] == 2 and row["recaptured_count"] == 1 and row["released_net_count"] == 1


def test_bromine_balance_holds_through_release_transfer_and_capture():
    spec = S.pvb_with_methyl_bromide({4: 2})
    state = materialize(spec)
    observer = mo.MeltObserver(spec, xs=0, side_group_homolysis=PVB_CHANNEL)
    tracker = observer.release_tracker(state)
    a, b = chain_carbons_in_order(state, 0), chain_carbons_in_order(state, 1)
    steps = []

    def step(edits):
        state.apply(edits)
        tracker.update(state)
        _balance_holds(observer, state, tracker)
        steps.append(observer.x_balance(state, tracker)["Br"])

    br = next(x for x in state.atoms[a[1]].edges if state.atoms[x].element == "Br")
    step([SetBond(a[1], br, 1.0, 0.0), SetRadical(a[1], 0, 1), SetRadical(br, 0, 1)])        # release
    h = min(x for x in state.atoms[b[2]].edges if state.atoms[x].element == "H")
    step([SetBond(b[2], h, 1.0, 0.0), SetBond(br, h, 0.0, 1.0), SetRadical(br, 1, 0),
          SetRadical(b[2], 0, 1)])                                                              # Br + H -> HBr
    assert steps[-1]["other_count"] >= 2 and steps[-1]["free_count"] == 0                      # HBr and CH3Br
    assert steps[0]["condensed_at_declared_sites_count"] == 7
    row = _channel(observer.observe(state, releases=tracker), "aliphatic_C-Br")
    assert row["released_count"] == 1
    assert observer.observe(state, releases=tracker)["mom_mass_law"]["undeclared_loss"]["mass_g"] > 0.0  # the H


def test_balance_holds_after_every_engine_event():
    engine = _h_loss_engine(seed=5)
    observer = mo.MeltObserver(engine.spec, xs=0, side_group_homolysis=H_CHANNEL)
    tracker = observer.release_tracker(engine)

    def check(eng, event):
        tracker(eng, event)
        _balance_holds(observer, eng.state, tracker)
    engine.run(max_events=20, check=check)
    assert engine.events_fired == 20
    assert _channel(observer.observe(engine, releases=tracker), "aliphatic_C-H")["released_count"] == 20



# --------------------------------------------------------------------------
# I-158 rework 2
# --------------------------------------------------------------------------
TWO_BR_CHANNELS = [S.sgh_channel("aliphatic_C-Br", "aliphatic", 1.0), S.sgh_channel("benzylic_C-Br", "benzylic", 1.0)]


def test_reversed_channel_order_gives_the_same_per_label_counts():
    spec = S.two_bromine_spec({4: 1})
    state = materialize(spec)
    forward = mo.MeltObserver(spec, xs=0, side_group_homolysis=TWO_BR_CHANNELS)
    tracker = forward.release_tracker(state)
    _remove_br(state, chain_carbons_in_order(state, 0)[2])      # the aliphatic Br of unit 2
    tracker.update(state)
    backward = mo.MeltObserver(spec, xs=0, side_group_homolysis=list(reversed(TWO_BR_CHANNELS)))
    for observer in (forward, backward):
        o = observer.observe(state, releases=tracker)
        ali, ben = _channel(o, "aliphatic_C-Br"), _channel(o, "benzylic_C-Br")
        assert ali["released_tracked"] and ben["released_tracked"]
        assert ali["released_count"] == 1 and ben["released_count"] == 0
    other = mo.MeltObserver(spec, xs=0, side_group_homolysis=TWO_BR_CHANNELS[:1])
    assert _channel(other.observe(state, releases=tracker), "aliphatic_C-Br")["released_tracked"] is False


def _failing_once(engine, name="_add"):
    original = getattr(engine, name)
    calls = {"n": 0}

    def once(*args, **kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            raise RuntimeError(f"injected in {name}")
        return original(*args, **kwargs)
    setattr(engine, name, once)


def test_export_resumes_after_a_recovered_post_commit_failure(tmp_path):
    engine = _h_loss_engine()
    path = tmp_path / "series.jsonl"
    mo.write_time_series(engine, [0.0], path, xs=0, side_group_homolysis=H_CHANNEL)
    _failing_once(engine)
    with pytest.raises(RuntimeError, match="injected"):
        mo.write_time_series(engine, [1.0e6], path, xs=0, side_group_homolysis=H_CHANNEL, append=True)
    del engine._add
    assert engine.events_fired == 1
    last = mo.write_time_series(engine, [1.0e6], path, xs=0, side_group_homolysis=H_CHANNEL, append=True)[-1]
    row = _channel(last, "aliphatic_C-H")
    assert engine.events_fired == 20
    assert row["released_tracked"] is True and row["released_count"] == 20


def test_observing_before_the_event_callback_keeps_continuity():
    engine = _h_loss_engine()
    observer = mo.MeltObserver(engine.spec, xs=0, side_group_homolysis=H_CHANNEL)
    tracker = observer.release_tracker(engine)

    def hook(eng, event):
        observer.observe(eng, releases=tracker)          # syncs the event first
        tracker(eng, event)
    engine.run(max_events=5, check=hook)
    row = _channel(observer.observe(engine, releases=tracker), "aliphatic_C-H")
    assert row["released_tracked"] is True and row["released_count"] == 5


def test_a_callback_after_missed_events_reports_untracked():
    engine = _h_loss_engine()
    observer = mo.MeltObserver(engine.spec, xs=0, side_group_homolysis=H_CHANNEL)
    tracker = observer.release_tracker(engine)
    engine.run(max_events=2)                             # two events the tracker never saw
    engine.run(max_events=1, check=tracker)
    row = _channel(observer.observe(engine, releases=tracker), "aliphatic_C-H")
    assert row["released_tracked"] is False and row["released_count"] is None


def test_an_attached_tracker_is_the_exports_tracker(tmp_path):
    engine = _h_loss_engine()
    observer = mo.MeltObserver(engine.spec, xs=0, side_group_homolysis=H_CHANNEL)
    tracker = observer.release_tracker(engine)
    engine.run(max_events=1, check=tracker)
    obs = mo.write_time_series(engine, [engine.time], tmp_path / "s.jsonl", xs=0, observer=observer)[-1]
    row = _channel(obs, "aliphatic_C-H")
    assert row["released_tracked"] is True and row["released_count"] == 1
    assert observer.release_tracker(engine) is tracker


def test_polymer_inventory_excludes_declared_sites_in_small_species():
    spec = S.pvb_with_methyl_bromide({4: 1})
    observer = mo.MeltObserver(spec, xs=0, side_group_homolysis=PVB_CHANNEL)
    state = materialize(spec)
    assert observer.site_inventory(state)[0] == {"tail": 4, "explicit": 0}      # not the CH3Br bromine
    balance = observer.x_balance(state)["Br"]
    assert balance["condensed_at_declared_sites_count"] == 4 and balance["other_count"] == 1


@pytest.mark.parametrize("case", ["naphthyl_b_alternating", "naphthyl_b_not_alternating", "para_bromo_ps",
                                  "pyridyl", "biphenylyl"])
def test_ring_system_aromaticity_equals_whole_molecule_perception(case):
    specs = {"naphthyl_b_alternating": lambda: S.bromonaphthyl_spec({3: 1}, "ring_b_alternating"),
             "naphthyl_b_not_alternating": lambda: S.bromonaphthyl_spec({3: 1}, "ring_b_not_alternating"),
             "para_bromo_ps": lambda: S.para_bromo_ps_spec({4: 1}),
             "pyridyl": lambda: S.pyridyl_spec({3: 1}),
             "biphenylyl": lambda: S.biphenylyl_bromide_spec({3: 1})}
    state = materialize(specs[case]())
    for cid in state.components.ids():
        assert mo.aromatic_atoms(state, cid) == mo.aromatic_atoms_whole(state, cid)
    if case != "pyridyl":
        assert mo.aromatic_atoms(state, state.components.ids()[0])
    br = [x for x in state.atoms if state.atoms[x].element == "Br"]
    if br:
        homolyse(state, next(iter(state.atoms[br[0]].edges)), br[0])        # an aryl sigma radical
        for cid in state.components.ids():
            assert mo.aromatic_atoms(state, cid) == mo.aromatic_atoms_whole(state, cid)


def test_aromatic_perception_cost_does_not_grow_with_chain_length():
    """Perception runs per ring system and is cached by its graph: the atoms
    handed to RMG for a whole chain's inventory do not depend on DP."""
    perceived = []
    for dp in (8, 16):
        spec = S.para_bromo_ps_spec({dp: 1})
        observer = mo.MeltObserver(spec, xs=0, side_group_homolysis=[S.sgh_channel("aryl_C-Br", "aryl", 1.0)])
        state = materialize(spec)
        mo.clear_aromaticity_cache()
        before = mo.PERCEPTION_STATS["atoms"]
        assert observer.site_inventory(state)[0]["tail"] == dp
        perceived.append(mo.PERCEPTION_STATS["atoms"] - before)
    assert perceived[0] == perceived[1] and perceived[0] < 40
