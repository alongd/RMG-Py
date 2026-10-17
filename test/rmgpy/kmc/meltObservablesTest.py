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
from rmgpy.kmc.melt_state import AncestryError, SetBond, SetRadical, materialize

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
    # a MoM run that reports exactly the kMC state agrees on every row
    rows = mo.side_by_side(last, mom)
    assert {r["tier"] for r in rows} == {1, 3}
    assert all(r["agreement"] in ("equal", "within_bound") for r in rows), rows


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
                "polymer_derived_mass_exact_g", "Mn_g_mol", "Mw_g_mol", "polymer_conversion",
                "polymer_derived_molecules_mol", "species_mol.methyl", "species_mol.H", "species_mol.H2"):
        assert key in rows, key
    # a MoM quantity the kMC run never produced still gets a row
    rows = {r["quantity"]: r for r in mo.side_by_side(o, {"explicit_dp_mol": {"7": 0.5}})}
    assert rows["explicit_dp_mol.7"]["kmc"] == 0.0 and rows["explicit_dp_mol.7"]["agreement"] == "differs"
    # no MoM values: no difference is claimed
    assert all(r["difference"] is None and r["agreement"] == "no_mom_values" for r in mo.side_by_side(o))


@pytest.mark.parametrize("label", ["H", "H2", "monomer_product"])
def test_a_declared_label_colliding_with_a_builtin_is_refused(label):
    spec = S.pe_spec({6: 1}, small={label: S.METHYL})
    with pytest.raises(me.DeclarationError, match="built-in"):
        mo.MeltObserver(spec, xs=0)
    with pytest.raises(me.DeclarationError, match="built-in"):
        me.MeltEngine(spec, S.ToyCatalog(builder=lambda T: S.pe_templates(T)), 700.0, np.random.default_rng(0))
