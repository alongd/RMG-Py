"""S4 red tests: the melt engine over a loaded template catalog (I-151).

Design: I-127 v3 §3 (propensity, exact pair sampler, role_swap), §4 (engine),
§7 (MET length and invariants), §8 (scope), §9 (S4 red tests).  Every
expected number below is derived independently of the engine: exact pair
enumeration, the declared transport formula written out here, whole-species
Kc from the catalog's stored thermo (``template_emission.kc_at``), and atom
counts read from the state.
"""

import math
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

from rmgpy import constants
from rmgpy.kmc import template_emission as te

import melt_engine_support as S

from rmgpy.kmc import melt_engine as me

ENGINE_SOURCE = Path(me.__file__).read_text()
PE_LIVE = ("homolysis", "recombination", "h_abstraction", "disproportionation", "beta_scission",
           "mid_end_recombination", "intra_mid_end")


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def totals(state):
    """(sum of compositions, sum of backbone counts), read from the atoms."""
    comp = Counter((n.element, n.isotope) for n in state.atoms.values())
    backbone = sum(1 for o in state.origin.values() if o.is_backbone)
    per_component = Counter()
    bb = 0
    for cid in state.components.ids():
        per_component.update(state.components.composition(cid))
        bb += state.components.backbone_count(cid)
    assert per_component == comp and bb == backbone
    return comp, backbone


def max_backbone_degree(state):
    origin = state.origin
    return max(sum(1 for b in state.atoms[a].edges if origin[b].is_backbone)
               for a in state.atoms if origin[a].is_backbone)


def k_diffusion(met_block, temperature, m0, n1, n2):
    """The declared encounter law, written out independently of the engine."""
    d = met_block["segment_diffusivity"]
    if d["law"] == "arrhenius":
        d0 = d["D_ref_m2_s"] * math.exp(-d["Ea_J_mol"] / constants.R * (1 / temperature - 1 / d["T_ref_K"]))
    else:
        def shift(x):
            return -d["C1"] * (x - d["T_wlf_K"]) / (d["C2_K"] + x - d["T_wlf_K"])
        d0 = d["D_ref_m2_s"] * 10 ** (shift(d["T_ref_K"]) - shift(temperature))
    scaling = met_block["chain_scaling"]

    def diff(n):
        if scaling["law"] == "rouse" or n <= scaling["entanglement_units"]:
            return d0 / n
        return d0 * scaling["entanglement_units"] / n ** 2
    rg = math.sqrt(met_block["mean_square_end_to_end_per_molar_mass_m2_mol_per_g"] * m0 * min(n1, n2) / 6)
    radius = max(met_block["contact_radius_m"], 2 * rg)
    return 4 * math.pi * constants.Na * met_block["spin_factor"] * (diff(n1) + diff(n2)) * radius


def two_state_ratio(engine, events):
    """Run ``events`` events of the dimer <-> radical-pair system; return
    (time in dimer / time in pair, completed sojourns)."""
    comps = len(engine.state.components)
    time_in = Counter()
    last = engine.time
    sojourns = 0
    for _ in range(events):
        before = len(engine.state.components)
        event = engine.step()
        assert event is not None
        time_in[before] += event.time - last
        last = event.time
        sojourns += 1
        assert len(engine.state.components) != before
    assert comps in time_in
    return time_in[1] / time_in[2], sojourns


# --------------------------------------------------------------------------
# C3': declared contents only, positive initiation
# --------------------------------------------------------------------------
@pytest.mark.parametrize("mobile", [{}, {"styrene": 2}])
def test_c3_closed_melt_from_declared_contents_has_positive_initiation(mobile):
    spec = S.ps_spec({10: 1, 20: 1, 30: 1}, mobile=mobile)
    catalog = S.ps_catalog()
    if mobile:
        # S5 (design §5): styrene is declared mobile but no template of the
        # sampled PS catalogs involves it, so the run is refused at start; the
        # declared contents are still materialized and identified as styrene.
        with pytest.raises(me.DeclarationError, match="styrene"):
            me.MeltEngine(spec, catalog, 700.0, np.random.default_rng(0))
        catalog = me.PayloadCatalog(catalog.payloads)
        state = me.materialize(spec)
        identity = me.SpeciesIdentity(spec)
        others = [c for c in state.components.ids() if state.origin_kinds(c) == {"mobile"}]
        assert len(others) == 2 and all(identity.identify(state, c) == "styrene" for c in others)
        spec = S.ps_spec({10: 1, 20: 1, 30: 1})
        mobile = {}
    engine = me.MeltEngine(spec, catalog, 700.0, np.random.default_rng(0))
    state = engine.state
    chains = [c for c in state.components.ids() if state.origin_kinds(c) == {"chain"}]
    others = [c for c in state.components.ids() if c not in chains]
    formulas = sorted(tuple(sorted(state.components.composition(c).items(), key=str)) for c in chains)
    assert formulas == sorted(((("C", None), 8 * n), (("H", None), 8 * n + 2)) for n in (10, 20, 30))
    assert len(others) == sum(mobile.values())
    assert all(engine.entries[c].species == "styrene" for c in others)
    assert all(state.components.radical_count(c) == 0 for c in state.components.ids())
    prop = engine.propensities()
    assert prop.enabled > 0.0
    live = [k for k, a in prop.channels if a > 0]
    assert live and all(k[0] == "uni" and engine.classes[k[1]].radical_delta > 0 for k in live)
    # the MET length unit: n = max(1, floor(B_c / b)), b = 2 backbone atoms per PS unit
    assert engine.transport.backbone_atoms_per_unit == 2
    assert sorted(engine.entries[c].met_length for c in chains) == [10, 20, 30]


# --------------------------------------------------------------------------
# index: incremental == full recompute
# --------------------------------------------------------------------------
def test_incremental_equals_full_recompute_after_200_events():
    spec = S.pe_spec({10: 1, 20: 1, 30: 1})
    engine = me.MeltEngine(spec, S.ToyCatalog(names=PE_LIVE), 700.0, np.random.default_rng(11))
    hab = next(c.index for c in engine.classes if c.name == "h_abstraction")
    split_x = merge_x = 0
    previous = engine.snapshot()
    for _ in range(200):
        event = engine.step()
        assert event is not None, "the toy melt went absorbing before 200 events"
        engine.check_consistency()
        mine, fresh = engine.snapshot(), engine.fresh_snapshot()
        assert mine == fresh
        x_before = {k: v for k, v in previous["X"] if k[0] == hab}
        x_after = {k: v for k, v in mine["X"] if k[0] == hab}
        if x_before != x_after:
            if len(event.products) > len(event.components):
                split_x += 1
            elif len(event.products) < len(event.components):
                merge_x += 1
        previous = mine
    assert split_x >= 1 and merge_x >= 1, (split_x, merge_x)


def test_incremental_equals_full_recompute_on_the_loaded_ps_catalog():
    """The same check on emitted templates (sampled closed-PS-DP8 R_Recombination
    and H_Abstraction catalogs) over a closed PS melt."""
    engine = me.MeltEngine(S.ps_spec({6: 3}), S.ps_catalog(), 700.0, np.random.default_rng(8))
    fired = Counter()
    for _ in range(40):
        event = engine.step()
        assert event is not None
        fired[event.channel] += 1
        engine.check_consistency()
        assert engine.snapshot() == engine.fresh_snapshot()
    assert fired["uni"] and fired["met"], fired


@pytest.mark.parametrize("quotient_catalog", ["ps", "pe"])
def test_hydrogen_quotient_equals_full_enumeration(quotient_catalog):
    if quotient_catalog == "ps":
        spec, catalog = S.ps_spec({6: 2}), S.ps_catalog()
    else:
        spec, catalog = S.pe_spec({8: 2}), S.ToyCatalog()
    engine = me.MeltEngine(spec, catalog, 700.0, np.random.default_rng(4))
    for _ in range(6):
        if engine.step() is None:
            break
    ids = engine.state.components.ids()
    fast = engine._scan(ids)
    engine.quotient = False
    full = engine._scan(ids)
    assert {c: e.canonical() for c, e in fast.items()} == {c: e.canonical() for c, e in full.items()}


# --------------------------------------------------------------------------
# two-state equilibrium and the MET escape factor
# --------------------------------------------------------------------------
def _pe2_engine(met_block, seed, temperature=700.0):
    catalog, state, raw = S.pe2_two_state()
    derived = [t for p in raw for t in p["templates"] if t.get("derived")]
    forward = next(t for t in derived if t["key"].startswith("165c69761932"))
    kc = te.kc_at(forward["derived"]["thermo"], temperature)   # dimer -> 2 radicals, mol/m^3
    assert len(forward["derived"]["thermo"]["reactants"]) == 1
    assert len(forward["derived"]["thermo"]["products"]) == 2
    volume = 2.0 / (kc * constants.Na)          # expected pi(dimer)/pi(pair) = 2/(Kc N_A V) = 1
    spec = S.ps_spec({2: 1}, met_block=met_block)
    engine = me.MeltEngine(spec, catalog, temperature, np.random.default_rng(seed),
                           state=state, volume_m3=volume)
    return engine, 2.0 / (kc * constants.Na * volume), kc


def test_two_state_equilibrium_reproduces_kc():
    """2 R. <-> R-R (identical radicals, role_swap).  The stationary ratio of
    time in R-R over time in 2 R. is 2 / (Kc N_A V) with Kc the whole-species
    dissociation constant from the catalog's stored thermo.  Tolerance: four
    standard errors of a ratio of two means of n exponential sojourns each,
    4 * sqrt(2 / n)."""
    fast = dict(S.ps_met(), segment_diffusivity={"law": "arrhenius", "D_ref_m2_s": 1.0,
                                                 "T_ref_K": 700.0, "Ea_J_mol": 0.0})
    engine, expected, _kc = _pe2_engine(fast, seed=20261017)
    (pair,) = engine.pairs.values()
    assert pair.phi > 1 - 1e-6
    ratio, n = two_state_ratio(engine, 3000)
    tolerance = 4 * math.sqrt(2.0 / (n / 2))
    assert abs(ratio / expected - 1) < tolerance, (ratio, expected, tolerance)


def test_met_escape_factor_preserves_k():
    """Same system with transport declared so that phi = 0.3: association is
    slowed by phi (Collins-Kimball on the summed activation), homolysis is
    thinned by the phi of its product pair, and K is unchanged."""
    catalog, state, _raw = S.pe2_two_state()
    probe = me.MeltEngine(S.ps_spec({2: 1}), catalog, 700.0, np.random.default_rng(0), state=state)
    (k_act,) = [p.k_act for p in probe.pairs.values()]
    m0 = probe.transport.repeat_unit_molar_mass_g_mol
    base = dict(S.ps_met(), segment_diffusivity={"law": "arrhenius", "D_ref_m2_s": 1.0,
                                                 "T_ref_K": 700.0, "Ea_J_mol": 0.0})
    unit = k_diffusion(base, 700.0, m0, 1, 1)
    target_kd = k_act * 0.3 / 0.7
    slow = dict(base, segment_diffusivity=dict(base["segment_diffusivity"], D_ref_m2_s=target_kd / unit))
    engine, expected, _kc = _pe2_engine(slow, seed=7)
    (pair,) = engine.pairs.values()
    kd = k_diffusion(slow, 700.0, m0, 1, 1)
    assert pair.k_d == pytest.approx(kd, rel=1e-12)
    assert pair.phi == pytest.approx(kd / (kd + k_act), rel=1e-12)
    assert pair.phi == pytest.approx(0.3, rel=1e-9)
    assert engine.propensities().enabled == pytest.approx(pair.phi * k_act / engine.nav, rel=1e-12)
    ratio, n = two_state_ratio(engine, 3000)
    tolerance = 4 * math.sqrt(2.0 / (n / 2))
    assert abs(ratio / expected - 1) < tolerance, (ratio, expected, tolerance)
    rejected = engine.met_counts["reverse_rejected"]
    accepted = engine.met_counts["reverse_accepted"]
    assert accepted / (accepted + rejected) == pytest.approx(0.3, abs=4 * math.sqrt(0.21 / (accepted + rejected)))
    assert engine.null_proposals == rejected


def test_ck_applied_once_per_component_pair():
    """Two radical chains, two MET classes between them (recombination,
    role_swap: T = 1; disproportionation, ordered, either chain the donor,
    either of the donor's two beta H: T = 4): the pair rate is
    CK of the SUMMED activation, and a class is drawn proportional to w T."""
    spec = S.pe_spec({12: 2})
    catalog = S.ToyCatalog(builder=lambda t: S.pe_templates(t, weights={"disproportionation": 5.0e-3}),
                           names=("recombination", "disproportionation"))
    engine = me.MeltEngine(spec, catalog, 700.0, np.random.default_rng(3))
    state = engine.state
    for chain in (0, 1):
        S.strip_h(state, S.end_carbons(state, chain)[0])
    engine.rebuild()
    w = {t.name: t.weight for t in engine.templates}
    k_act = w["recombination"] * 1 + w["disproportionation"] * 4
    chains = sorted(c for c in state.components.ids() if state.components.radical_count(c))
    chains = [c for c in chains if state.components.backbone_count(c)]
    assert len(chains) == 2
    (key, pair), = [(k, p) for k, p in engine.pairs.items() if set(k) == set(chains)]
    m0 = engine.transport.repeat_unit_molar_mass_g_mol
    kd = k_diffusion(S.pe_met(), 700.0, m0, 12, 12)
    assert pair.k_act == pytest.approx(k_act, rel=1e-12)
    assert pair.k_d == pytest.approx(kd, rel=1e-12)
    phi = kd / (kd + k_act)
    assert 0.2 < phi < 0.8   # the per-pair and per-site-type sums differ materially here
    assert pair.allowed == pytest.approx(phi * k_act / engine.nav, rel=1e-12)
    met = dict(engine.propensities().channels)[("met",)]
    assert met == pytest.approx(phi * k_act / engine.nav, rel=1e-12)
    # within the pair: class ~ w T
    rng_draws = Counter()
    for _ in range(4000):
        tr, _view, _rate = engine._draw_met(met)
        rng_draws[tr.template] += 1
    p_rec = w["recombination"] / k_act
    n = sum(rng_draws.values())
    assert abs(rng_draws["recombination"] / n - p_rec) < 4 * math.sqrt(p_rec * (1 - p_rec) / n)


def test_declared_small_species_bypass_met_and_is_counted():
    spec = S.pe_spec({12: 1}, small={"methyl": S.METHYL}, mobile={"methyl": 1})
    engine = me.MeltEngine(spec, S.ToyCatalog(names=("methyl_recombination",)), 700.0,
                           np.random.default_rng(1))
    state = engine.state
    S.strip_h(state, S.end_carbons(state, 0)[0])
    engine.rebuild()
    methyl = next(c for c in state.components.ids() if engine.entries[c].species == "methyl")
    (pair,) = [p for k, p in engine.pairs.items() if methyl in k]
    assert pair.bypass and pair.phi == 1.0
    assert engine.step().template == "methyl_recombination"
    assert engine.met_counts["forward_bypass"] == 1


# --------------------------------------------------------------------------
# branching leak
# --------------------------------------------------------------------------
def test_branch_creating_transitions_accrue_as_leak():
    spec = S.pe_spec({14: 2})
    names = ("recombination", "mid_end_recombination", "intra_mid_end")
    engine = me.MeltEngine(spec, S.ToyCatalog(names=names), 700.0, np.random.default_rng(9))
    state = engine.state
    c0 = S.carbons(state, 0)
    mid = sorted(c0, key=lambda a: a)[len(c0) // 2]
    S.strip_h(state, mid)
    S.strip_h(state, S.end_carbons(state, 0)[0])
    S.strip_h(state, S.end_carbons(state, 1)[0])
    engine.rebuild()
    prop = engine.propensities()
    assert prop.leak > 0.0
    w = {t.name: t.weight for t in engine.templates}
    # intramolecular: exactly one mid/end transition on chain 0 (uni, leak)
    intra = engine.classes[[t.name for t in engine.templates].index("intra_mid_end")].index
    assert engine.uni_allowed[intra] == 0 and engine.uni_leak[intra] == 1
    assert dict(prop.channels).get(("uni", intra)) == 0.0
    met_leak = sum(p.leak for p in engine.pairs.values())
    assert prop.leak == pytest.approx(w["intra_mid_end"] * 1 + met_leak, rel=1e-12)
    assert met_leak > 0.0
    fired = Counter()

    def check(eng, event):
        fired[event.template] += 1
        assert max_backbone_degree(eng.state) <= 2

    result = engine.run(max_events=1, check=check)
    assert fired == {"recombination": 1}
    assert result.leak_num > 0.0 and 0.0 < result.leak < 1.0
    assert engine.report()["leak_fraction"] == result.leak
    # a live toy melt: branch-creating classes never fire, their propensity accrues
    live = me.MeltEngine(S.pe_spec({10: 1, 20: 1, 30: 1}), S.ToyCatalog(names=PE_LIVE), 700.0,
                         np.random.default_rng(5))
    fired.clear()
    result = live.run(max_events=200, check=check)
    assert result.events == 200
    assert fired["mid_end_recombination"] == 0 and fired["intra_mid_end"] == 0
    assert fired["h_abstraction"] > 0 and fired["beta_scission"] > 0   # mid radicals existed
    assert result.leak_num > 0.0


# --------------------------------------------------------------------------
# exact pair sampler
# --------------------------------------------------------------------------
def _sites(prefix, cid, n):
    return tuple((prefix, cid, k) for k in range(n))


@pytest.mark.parametrize("layout", [
    {0: (2, 1), 5: (0, 2), 9: (3, 0), 14: (1, 2), 20: (2, 3), 27: (1, 1)},
    {1: (4, 4), 3: (1, 0), 4: (0, 1), 8: (1, 1), 11: (2, 0), 30: (0, 3)},
])
def test_pair_sampler_chi2_against_exact_enumeration(layout):
    sampler = me.ComponentPairSampler(32)
    for cid in sorted(layout):
        sampler.set(cid, (), _sites("x", cid, 5))      # dummy content first: updates must replace it
    for cid, (n0, n1) in sorted(layout.items()):
        sampler.set(cid, _sites("a", cid, n0), _sites("b", cid, n1))
    pairs = [(s0, s1) for c0, (n0, _x) in layout.items() for c1, (_y, n1) in layout.items() if c0 != c1
             for s0 in _sites("a", c0, n0) for s1 in _sites("b", c1, n1)]
    assert sampler.ordered_total == len(pairs)
    n0 = sum(v[0] for v in layout.values())
    n1 = sum(v[1] for v in layout.values())
    assert len(pairs) == n0 * n1 - sum(a * b for a, b in layout.values())
    rng = np.random.default_rng(12345)
    draws = 200 * len(pairs)
    counts = Counter()
    for _ in range(draws):
        c0, s0, c1, s1 = sampler.draw(rng)
        assert c0 != c1 and s0[1] == c0 and s1[1] == c1
        counts[(s0, s1)] += 1
    assert set(counts) <= set(pairs)
    observed = np.array([counts[p] for p in pairs], dtype=float)
    chi2 = float(((observed - 200.0) ** 2 / 200.0).sum())
    assert chi2 < stats.chi2.ppf(0.999, len(pairs) - 1), chi2


# --------------------------------------------------------------------------
# declaration: generic transport, refused when absent
# --------------------------------------------------------------------------
def test_run_start_refuses_a_declaration_without_transport():
    spec = S.ps_spec({10: 1}, transport=False)
    with pytest.raises(me.DeclarationError, match="transport_model"):
        me.MeltEngine(spec, S.ps_catalog(), 700.0, np.random.default_rng(0))
    broken = dict(S.ps_met())
    del broken["contact_radius_m"]
    with pytest.raises(me.DeclarationError, match="contact_radius_m"):
        me.MeltEngine(S.ps_spec({10: 1}, met_block=broken), S.ps_catalog(), 700.0,
                      np.random.default_rng(0))


def test_engine_names_no_polymer_and_reads_no_ps_constants():
    code = re.sub(r'"""[\s\S]*?"""', "", ENGINE_SOURCE)
    code = re.sub(r"#.*", "", code)
    assert not re.search(r"\b(PS|PE|PP|styrene|polystyrene|polyethylene|phenyl)\b", code, re.I)
    for name in ("PS_M0", "PS_C_R2", "SIGMA_CONTACT", "TRANSPORT_ARMS", "_D_FERRY", "rmgpy.kmc.met"):
        assert name not in ENGINE_SOURCE


def test_non_ps_linear_polymer_runs_200_events():
    spec = S.pe_spec({10: 1, 20: 1, 30: 1})
    engine = me.MeltEngine(spec, S.ToyCatalog(names=PE_LIVE), 700.0, np.random.default_rng(5))
    assert engine.transport.backbone_atoms_per_unit == 2
    assert engine.transport.repeat_unit_molar_mass_g_mol == pytest.approx(28.05, abs=0.01)
    start = totals(engine.state)

    def check(eng, event):
        assert totals(eng.state) == start
        for cid, entry in eng.entries.items():
            assert entry.met_length == max(1, eng.state.components.backbone_count(cid) // 2)

    result = engine.run(max_events=200, check=check)
    assert result.events == 200 and result.stop_reason == "max_events"
    assert engine.conservation_checks == 200
    assert result.family_counts["toy_homolysis"] > 0 and result.family_counts["toy_recombination"] > 0
    assert sum(result.family_counts.values()) == 200


# --------------------------------------------------------------------------
# clock and temperature
# --------------------------------------------------------------------------
def _trace(engine):
    return [(e.time, e.template, e.components, e.products) for e in engine.events]


def test_horizon_schedule_does_not_change_the_trajectory():
    spec = S.pe_spec({10: 1, 20: 1})
    whole = me.MeltEngine(spec, S.ToyCatalog(names=PE_LIVE), 700.0, np.random.default_rng(42))
    whole.run(until=20.0)
    pieces = me.MeltEngine(spec, S.ToyCatalog(names=PE_LIVE), 700.0, np.random.default_rng(42))
    for horizon in np.linspace(0.37, 20.0, 41):
        pieces.run(until=float(horizon))
    assert _trace(whole) == _trace(pieces) and len(whole.events) > 20
    assert whole.state.state_hash() == pieces.state.state_hash()
    assert whole.null_proposals == pieces.null_proposals > 0
    assert whole.leak_num == pytest.approx(pieces.leak_num, rel=1e-12)


def test_temperature_refreshes_weights_through_the_catalog():
    catalog = S.ToyCatalog(builder=lambda t: S.pe_templates(t, scale=t / 700.0), names=PE_LIVE)
    engine = me.MeltEngine(S.pe_spec({10: 2}), catalog, 700.0, np.random.default_rng(2),
                           temperature_range=(650.0, 850.0))
    before = engine.propensities().enabled
    engine.set_temperature(840.0)
    assert catalog.requested == [700.0, 840.0]
    assert engine.temperature == 840.0
    assert engine.propensities().enabled == pytest.approx(before * 840.0 / 700.0, rel=1e-12)
    assert engine.snapshot() == engine.fresh_snapshot()
    with pytest.raises(me.DeclarationError, match="outside the run range"):
        engine.set_temperature(900.0)


def test_non_natural_isotopes_are_refused_at_run_start():
    from rmgpy.kmc.melt_state import MeltState, Node
    spec = S.pe_spec({4: 1})
    state = S.witness_state(["CCCC"])
    atoms = dict(state.atoms)
    first = min(a for a in atoms if atoms[a].element == "C")
    n = atoms[first]
    atoms[first] = Node(n.element, 13, n.charge, n.radical, n.lone_pairs, n.edges)
    labelled = MeltState(atoms, state.origin)
    with pytest.raises(me.DeclarationError, match="isotope"):
        me.MeltEngine(spec, S.ToyCatalog(), 700.0, np.random.default_rng(0), state=labelled)


# --------------------------------------------------------------------------
# rework 1
# --------------------------------------------------------------------------
def test_declared_species_identity_is_graph_isomorphism_not_an_attribute_multiset():
    """Crotyl [CH2]C=CC is declared; its isomer [CH2]CC=C (but-3-en-1-yl) is
    not.  Identity must be exact isomorphism under resonance, on the cached
    path (crotyl first, then the isomer) and on a fresh identity."""
    spec = S.pe_spec({6: 1}, small={"crotyl": S.CROTYL}, mobile={"crotyl": 1})
    crotyl_state = S.witness_state(["[CH2]C=CC", "C=CC[CH2]"])
    first = crotyl_state.components.ids()
    assert len(first) == 2
    cached = me.SpeciesIdentity(spec)
    assert cached.identify(crotyl_state, first[0]) == "crotyl"
    assert cached.identify(crotyl_state, first[1]) is None           # the isomer, after a cache fill
    fresh = me.SpeciesIdentity(spec)
    assert fresh.identify(crotyl_state, first[1]) is None            # the isomer alone
    assert fresh.identify(crotyl_state, first[0]) == "crotyl"
    # the other resonance form of crotyl is crotyl
    other = S.witness_state(["C=C[CH]C"])
    assert me.SpeciesIdentity(spec).identify(other, other.components.ids()[0]) == "crotyl"
    # the engine: only the declared species bypasses MET
    # (S5: a declared mobile species must take part in a loaded template)
    engine = me.MeltEngine(spec, S.ToyCatalog(builder=S.allyl_recombination), 700.0,
                           np.random.default_rng(0), state=crotyl_state)
    assert [engine.entries[c].species for c in first] == ["crotyl", None]


def _leak_fraction_state():
    """PE DP12 with a mid radical and an end radical on one chain (the
    reviewer's case): proposals of homolysis are thinned hard, the intra
    mid/end branching class carries the leak; slow declared transport makes
    phi of the homolysis products small)."""
    spec = S.pe_spec({12: 1}, d_ref=3.0e-20)
    engine = me.MeltEngine(spec, S.ToyCatalog(names=("homolysis", "recombination", "intra_mid_end")),
                           700.0, np.random.default_rng(31))
    state = engine.state
    c = S.carbons(state, 0)
    S.strip_h(state, sorted(c)[len(c) // 2])
    S.strip_h(state, S.end_carbons(state, 0)[0])
    engine.rebuild()
    return engine


def _physical_reverse_hazard(engine):
    """Exact sum over the phi-thinned (radical-increasing) transitions of w * phi
    of their product pair: each transition applied tentatively and rolled back."""
    total = 0.0
    state = engine.state
    for cls in engine.classes:
        if not cls.met_reverse or cls.bimolecular:
            continue
        for entry in list(engine.entries.values()):
            for tr in entry.uni.get(cls.index, ((), 0))[0]:
                view = me._LocalView(state, list(tr.components))
                tx, products = engine._begin_apply(tr, view)
                try:
                    phi = 1.0
                    if len(products) == 2:
                        new = engine._scan(products)
                        rate = engine.pair_rate(new[products[0]], new[products[1]])
                        phi = 1.0 if rate is None else rate.phi
                finally:
                    tx.rollback()
                total += engine.weights[cls.index] * phi
    return total


def test_leak_fraction_is_measured_against_physical_hazards():
    engine = _leak_fraction_state()
    before = engine.state.state_hash()
    prop = engine.propensities()
    proposal = prop.enabled
    physical = _physical_reverse_hazard(engine)
    assert engine.state.state_hash() == before
    assert physical < 0.05 * proposal          # the thinning is strong here
    exact = prop.leak / (prop.leak + physical)
    naive = prop.leak / (prop.leak + proposal)
    result = engine.run(until=40.0, max_events=1)
    assert engine.null_proposals >= 40
    # the reported fraction estimates the physical one, not the proposal-based one
    assert abs(result.leak - exact) < 0.03, (result.leak, exact, naive)
    assert abs(result.leak - naive) > 0.5
    assert engine.report()["leak_fraction"] == result.leak


def test_time_is_monotone_across_a_resync_after_a_horizon():
    """Observation horizon, then an external edit: the resync redraws from
    the current time (the reviewer's case: PE DP6 homolysis only, seed 1)."""
    engine = me.MeltEngine(S.pe_spec({6: 1}), S.ToyCatalog(names=("homolysis",)), 700.0,
                           np.random.default_rng(1))
    horizon = 0.5748369784
    result = engine.run(until=horizon)
    assert result.events == 0 and engine.time == horizon
    leak_den = engine.leak_den
    S.strip_h(engine.state, S.end_carbons(engine.state, 0)[0])
    event = engine.step()
    assert event is not None and event.time > horizon
    assert engine.leak_den >= leak_den
    times = [e.time for e in engine.events]
    engine.run(max_events=20)
    times = [e.time for e in engine.events]
    assert times == sorted(times) and times[0] > horizon


def test_reverse_thinning_is_transactional():
    """An exception inside the tentative apply / product scan / transport
    evaluation of a thinned reverse leaves the state byte-identical and no
    transaction open, and the exception propagates."""
    for target in ("pair_rate", "_scan"):
        engine = me.MeltEngine(S.pe_spec({10: 1}), S.ToyCatalog(names=("homolysis", "recombination")),
                               700.0, np.random.default_rng(2))
        before = engine.state.state_hash()
        revision = engine.state.revision

        def boom(*args, **kwargs):
            raise RuntimeError(f"injected in {target}")

        setattr(engine, target, boom)
        with pytest.raises(RuntimeError, match="injected"):
            engine.step()
        assert engine.state.state_hash() == before
        assert engine.state.revision == revision
        assert engine.state._tx is None


@pytest.mark.parametrize("arm", ["wlf", "arrhenius_nonfinite"])
def test_singular_transport_is_refused_at_run_start(arm):
    if arm == "wlf":
        # T_wlf - C2 = 714.7 K: the WLF shift is singular between it and the run temperature
        block = dict(S.ps_met(), segment_diffusivity={"law": "wlf", "D_ref_m2_s": 1e-12, "T_ref_K": 800.0,
                                                      "C1": 13.39, "C2_K": 45.28, "T_wlf_K": 760.0})
    else:
        # exp(Ea/R (1/T_ref - 1/T)) overflows at 700 K
        block = dict(S.ps_met(), segment_diffusivity={"law": "arrhenius", "D_ref_m2_s": 1e-12,
                                                      "T_ref_K": 600.0, "Ea_J_mol": 5.0e7})
    spec = S.ps_spec({10: 1}, met_block=block)
    with pytest.raises(me.DeclarationError, match="singular|finite"):
        me.MeltEngine(spec, S.ps_catalog(), 700.0, np.random.default_rng(0))


@pytest.mark.parametrize("scaling", [{"law": "rouse"}, {"law": "crossover", "entanglement_units": 10.0}])
@pytest.mark.parametrize("diffusivity", [
    {"law": "arrhenius", "D_ref_m2_s": 3.0e-11, "T_ref_K": 650.0, "Ea_J_mol": 4.0e4},
    {"law": "wlf", "D_ref_m2_s": 2.0e-12, "T_ref_K": 473.0, "C1": 13.39, "C2_K": 45.28, "T_wlf_K": 378.15},
])
def test_k_diffusion_is_the_declared_law_for_unequal_lengths(scaling, diffusivity):
    block = dict(S.pe_met(), segment_diffusivity=diffusivity, chain_scaling=scaling)
    law = me.TransportLaw.from_spec(S.pe_spec({6: 1}, met_block=block))
    m0 = law.repeat_unit_molar_mass_g_mol
    for n1, n2 in ((3, 40), (40, 3), (7, 25), (12, 12)):
        expected = k_diffusion(block, 700.0, m0, n1, n2)
        assert law.k_diffusion(700.0, n1, n2) == pytest.approx(expected, rel=1e-12)
    # the capture radius is set by the SHORTER chain: max-length Rg would give a different value here
    rg_short = math.sqrt(block["mean_square_end_to_end_per_molar_mass_m2_mol_per_g"] * m0 * 3 / 6)
    assert 2 * rg_short > block["contact_radius_m"]


def test_phi_uses_the_whole_activation_including_leak_site_pairs():
    """A component pair with an allowed MET site pair (end + end) and a
    branch-creating one (mid + end): phi = k_D / (k_D + k_act) with k_act over
    BOTH, and the allowed/leak split is phi * w T / (N_A V) each."""
    spec = S.pe_spec({14: 2})
    engine = me.MeltEngine(spec, S.ToyCatalog(names=("recombination", "mid_end_recombination")), 700.0,
                           np.random.default_rng(9))
    state = engine.state
    c0 = S.carbons(state, 0)
    S.strip_h(state, sorted(c0)[len(c0) // 2])
    S.strip_h(state, S.end_carbons(state, 0)[0])
    S.strip_h(state, S.end_carbons(state, 1)[0])
    engine.rebuild()
    w = {t.name: t.weight for t in engine.templates}
    chains = [c for c in state.components.ids() if state.components.backbone_count(c)]
    (pair,) = [p for k, p in engine.pairs.items() if set(k) == set(chains)]
    m0 = engine.transport.repeat_unit_molar_mass_g_mol
    kd = k_diffusion(S.pe_met(), 700.0, m0, 14, 14)
    k_allowed = w["recombination"] * 1
    k_leak = w["mid_end_recombination"] * 1
    k_act = k_allowed + k_leak
    phi = kd / (kd + k_act)
    assert pair.k_act == pytest.approx(k_act, rel=1e-12)
    assert pair.phi == pytest.approx(phi, rel=1e-12)
    assert abs(pair.phi - kd / (kd + k_allowed)) > 0.05 * pair.phi
    assert pair.allowed == pytest.approx(phi * k_allowed / engine.nav, rel=1e-12)
    assert pair.leak == pytest.approx(phi * k_leak / engine.nav, rel=1e-12)


def test_hydrogen_order_cache_holds_templates_weakly():
    import gc
    templates = S.pe_templates()
    for t in templates:
        me._h_order(t)
    alive = len(me._H_ORDER)
    del templates, t
    gc.collect()
    assert len(me._H_ORDER) <= alive - len(S.PE_WEIGHTS)


# --------------------------------------------------------------------------
# rework 2
# --------------------------------------------------------------------------
@pytest.mark.parametrize("declared_kekule", [True, False])
def test_declared_benzyl_is_identified_in_any_representation_and_bypasses_met(declared_kekule):
    """The query is normalized like the declared forms: a Kekule-materialized,
    an aromatic-materialized and a ring-radical resonance form of benzyl are
    all 'benzyl', on a cached identity and on fresh ones; the engine pair of
    two declared benzyls bypasses MET (phi = 1)."""
    declared = {"benzyl": S.benzyl(declared_kekule)}
    queries = []
    for kekule in (True, False):
        state = me.materialize(S.pe_spec({6: 1}, small={"benzyl": S.benzyl(kekule)}, mobile={"benzyl": 1}))
        queries.append((state, [c for c in state.components.ids() if state.origin_kinds(c) == {"mobile"}][0]))
    for smiles in ("[CH2]C1=CC=CC=C1", "[CH2]c1ccccc1", "C=C1C=C[CH]C=C1", "C=C1[CH]C=CC=C1"):
        state = S.witness_state([smiles])
        queries.append((state, state.components.ids()[0]))
    spec = S.pe_spec({6: 1}, small=declared, mobile={"benzyl": 2})
    cached = me.SpeciesIdentity(spec)
    for state, cid in queries:
        assert cached.identify(state, cid) == "benzyl"
        assert me.SpeciesIdentity(spec).identify(state, cid) == "benzyl"
    tolyl = S.witness_state(["Cc1ccc[c]c1"])          # an isomer (m-tolyl, sigma radical): not benzyl
    assert cached.identify(tolyl, tolyl.components.ids()[0]) is None
    engine = me.MeltEngine(spec, S.ToyCatalog(builder=S.benzyl_recombination), 700.0,
                           np.random.default_rng(0))
    mobiles = [c for c in engine.state.components.ids() if engine.state.origin_kinds(c) == {"mobile"}]
    assert [engine.entries[c].species for c in mobiles] == ["benzyl", "benzyl"]
    (pair,) = engine.pairs.values()
    assert pair.bypass and pair.phi == 1.0


def _two_end_radicals(seed=4):
    engine = me.MeltEngine(S.pe_spec({10: 2}), S.ToyCatalog(names=("recombination",)), 700.0,
                           np.random.default_rng(seed))
    for chain in (0, 1):
        S.strip_h(engine.state, S.end_carbons(engine.state, chain)[0])
    engine.rebuild()
    return engine


@pytest.mark.parametrize("target", ["_scan", "_add"])
def test_post_commit_failure_logs_the_event_and_leaves_a_consistent_index(target):
    engine = _two_end_radicals()
    original = getattr(engine, target)
    calls = {"n": 0}

    def once(*args, **kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            raise RuntimeError(f"injected in {target}")
        return original(*args, **kwargs)

    before = engine.state.revision
    setattr(engine, target, once)
    with pytest.raises(RuntimeError, match="injected"):
        engine.step()
    delattr(engine, target)
    assert engine.state.revision == before + 1                  # the chemistry is committed
    assert len(engine.events) == 1 and engine.events_fired == 1
    assert engine.events[0].template == "recombination"
    assert set(engine.entries) == set(engine.state.components.ids())
    engine.check_consistency()
    assert engine.snapshot() == engine.fresh_snapshot()


def _fast_leak_state():
    spec = S.pe_spec({12: 1})       # fast declared transport: phi of the homolysis products ~ 1
    engine = me.MeltEngine(spec, S.ToyCatalog(names=("homolysis", "recombination", "intra_mid_end")),
                           700.0, np.random.default_rng(31))
    state = engine.state
    c = S.carbons(state, 0)
    S.strip_h(state, sorted(c)[len(c) // 2])
    S.strip_h(state, S.end_carbons(state, 0)[0])
    engine.rebuild()
    return engine


def test_short_horizon_leak_report_states_its_bias_and_keeps_the_parts():
    """No proposal yet: the thinned part of the denominator is 0, so the fraction is
    high here (a single-run ratio estimate, bias sign unspecified); the open-segment hint is reported but not in leak_den."""
    engine = _fast_leak_state()
    prop = engine.propensities()
    exact = prop.leak / (prop.leak + _physical_reverse_hazard(engine))
    result = engine.run(until=1.0e-6)
    assert result.events == 0 and engine.null_proposals == 0
    assert result.leak > exact
    assert result.leak == pytest.approx(prop.leak / (prop.leak + (prop.enabled - prop.thinned)),
                                        rel=1e-9)
    parts = engine.report()["leak_estimate"]
    assert parts["estimate"] is True and parts["estimator"] == me.LEAK_ESTIMATOR
    assert parts["leak_num"] == engine.leak_num and parts["leak_den"] == engine.leak_den
    assert set(parts) >= {"non_thinned_integral", "phi_sum", "thinned_proposals", "open_thinned"}
    assert parts["open_thinned"] == pytest.approx(
        _physical_reverse_hazard(engine) * 1.0e-6, rel=1e-9)     # the hint is the exact value here
    assert "estimate" in me.LEAK_DEFINITION and "unspecified bias sign" in me.LEAK_DEFINITION


def test_long_run_leak_estimate_matches_the_physical_fraction():
    engine = _leak_fraction_state()
    prop = engine.propensities()
    exact = prop.leak / (prop.leak + _physical_reverse_hazard(engine))
    result = engine.run(until=40.0, max_events=1)
    assert engine.null_proposals >= 40
    assert abs(result.leak - exact) < 0.03
    parts = engine.report()["leak_estimate"]
    assert parts["thinned_proposals"] == engine.null_proposals + result.events


def test_unchanged_rebuild_keeps_the_pending_draw():
    def engine():
        return me.MeltEngine(S.pe_spec({6: 1}), S.ToyCatalog(names=("homolysis",)), 700.0,
                             np.random.default_rng(1))
    horizon = 0.5748369784
    plain, rebuilt = engine(), engine()
    plain.run(until=horizon)
    rebuilt.run(until=horizon)
    rebuilt.rebuild()
    a, b = plain.step(), rebuilt.step()
    assert a.time == b.time and a.time > horizon
    assert rebuilt.snapshot() == plain.snapshot()


# --------------------------------------------------------------------------
# I-152: one consistent leak estimator, labelled samples, counters vs the log
# --------------------------------------------------------------------------
def _exact_current_phi(engine, ci):
    """Exhaustive mean phi over every current transition of unimolecular class ci."""
    phis = []
    for c in sorted(engine.entries):
        for tr in engine.entries[c].uni.get(ci, ((), 0))[0]:
            tx, products = engine._begin_apply(tr, me._LocalView(engine.state, list(tr.components)))
            try:
                phi = 1.0
                if len(products) == 2:
                    new = engine._scan(products)
                    rate = engine.pair_rate(new[products[0]], new[products[1]])
                    phi = 1.0 if rate is None else rate.phi
            finally:
                tx.rollback()
            phis.append(phi)
    return math.fsum(phis) / len(phis)


def test_horizon_view_adds_no_biased_term_to_the_unbiased_phi_sum():
    """leak_den = leak_num + non-thinned integral + sum of phi over proposals, always:
    the open-segment estimate is a reported hint, not a term of the denominator."""
    engine = _leak_fraction_state()
    engine.run(until=0.7)
    parts = engine.leak_estimate()
    assert parts["open_thinned"] > 0.0
    assert engine.leak_den == pytest.approx(
        engine.leak_num + parts["non_thinned_integral"] + parts["phi_sum"], rel=1e-12, abs=0.0)


def test_leak_estimator_is_unbiased_in_the_denominator_over_seeds():
    """Mean over seeds of the physical part of the denominator at a horizon equals the
    exact integral (apply/rollback enumeration) of the physical firing hazard.  A
    one-shot regime (about one proposal per run), where an open-tail term shows."""
    probe = _leak_fraction_state()
    prop = probe.propensities()
    physical = _physical_reverse_hazard(probe)
    rate = (prop.enabled - prop.thinned) + physical
    horizon = 1.0 / prop.thinned
    dens = []
    for seed in range(240):
        engine = _leak_fraction_state()
        engine.rng = np.random.default_rng(1000 + seed)
        engine.run(until=horizon)
        dens.append(engine.leak_den - engine.leak_num)
    # state changes by a rare accepted event bias this by well under the tolerance
    assert np.mean(dens) == pytest.approx(rate * horizon, rel=0.12)


def test_phi_on_a_stride_subsample_is_labelled_a_sample():
    engine = me.MeltEngine(S.pe_spec({5: 128}), S.ToyCatalog(names=("homolysis",)), 700.0,
                           np.random.default_rng(1))
    ci = next(c.index for c in engine.classes if c.met_reverse)
    ntr = sum(len(en.uni.get(ci, ((), 0))[0]) for en in engine.entries.values())
    assert ntr > me.PHI_FALLBACK_LIMIT
    engine.run(until=1.0e-12)
    parts = engine.leak_estimate()
    assert parts["thinned_proposals"] == 0 and parts["open_thinned"] > 0.0
    exhaustive = _exact_current_phi(engine, ci) * engine.propensities().thinned * 1.0e-12
    if "exact" in parts["open_phi_source"]:
        assert parts["open_thinned"] == pytest.approx(exhaustive, rel=1e-12)
    else:
        assert "sample" in parts["open_phi_source"]
    assert "exact" not in parts["open_phi_source"]          # 512+ transitions are sampled here


def test_phi_over_few_transitions_is_exhaustive_and_labelled_exact():
    engine = me.MeltEngine(S.pe_spec({5: 4}), S.ToyCatalog(names=("homolysis",)), 700.0,
                           np.random.default_rng(1))
    ci = next(c.index for c in engine.classes if c.met_reverse)
    engine.run(until=1.0e-12)
    parts = engine.leak_estimate()
    assert parts["open_phi_source"] == "exact_current_transitions"
    assert parts["open_thinned"] == pytest.approx(
        _exact_current_phi(engine, ci) * engine.propensities().thinned * 1.0e-12, rel=1e-12)


def test_counters_agree_with_the_event_log_under_recovery():
    """An event committed and then lost to a post-commit failure is in met_counts
    and the family/template counters exactly when it is in the log."""
    engine = me.MeltEngine(S.pe_spec({10: 2}), S.ToyCatalog(names=("recombination",)),
                           700.0, np.random.default_rng(4))
    for chain in (0, 1):
        S.strip_h(engine.state, S.end_carbons(engine.state, chain)[0])
    engine.rebuild()
    def failing(*args, **kwargs):
        raise RuntimeError("injected")

    engine._update = failing
    with pytest.raises(RuntimeError, match="injected"):
        engine.step()
    del engine._update
    assert len(engine.events) == 1
    ev = engine.events[0]
    assert sum(engine.template_counts.values()) == sum(engine.family_counts.values()) == 1
    logged = Counter(e.channel for e in engine.events)
    forward = engine.met_counts["forward"] + engine.met_counts["forward_bypass"]
    assert forward == logged["met"]
    assert ev.channel == "met"


def test_running_mean_open_term_is_phi_sum_over_n_not_n_plus_one():
    engine = _leak_fraction_state()
    engine.run(until=3.0, max_events=10 ** 6)
    n = engine.thinned_proposals
    if engine.events_fired:
        pytest.skip("state moved")
    assert n >= 2
    assert len(engine._proposals_by_class) == 1
    engine.run(until=3.0 + 1.0e-3)
    prop = engine.propensities()
    elapsed = engine.time - engine._segment_start
    expected = prop.thinned * (engine.phi_sum / n) * elapsed
    assert engine.leak_estimate()["open_thinned"] == pytest.approx(expected, rel=1e-9)


GOLDEN_TRAJECTORY = "fa6f80838f4de91a3ec81910b60e9d932906e8cf0a6be91c6472a524c5a0f4b6"


def test_fixed_seed_trajectory_is_byte_identical_to_the_pre_i152_engine():
    import hashlib
    spec = S.pe_spec({10: 1, 6: 1})
    e = me.MeltEngine(spec, S.ToyCatalog(names=PE_LIVE), 700.0, np.random.default_rng(42))
    for h in np.linspace(0.37, 20.0, 41):
        e.run(until=float(h))
    tr = [(repr(x.time), x.template, x.components, x.products) for x in e.events]
    digest = repr((tr, e.state.state_hash(), e.null_proposals, repr(e.leak_num), e.events_fired,
                   dict(e.template_counts)))
    assert len(tr) == 136 and e.null_proposals == 18
    assert hashlib.sha256(digest.encode()).hexdigest() == GOLDEN_TRAJECTORY


# --------------------------------------------------------------------------
# I-152 rework 1: no cancellation in the non-thinned integral; pooled ensemble estimate
# --------------------------------------------------------------------------
def test_non_thinned_integral_does_not_cancel_against_a_huge_thinned_hazard():
    weights = {"homolysis": 1.0e17, "recombination": 1.0e-2}
    engine = me.MeltEngine(
        S.pe_spec({10: 2}),
        S.ToyCatalog(builder=lambda T: S.pe_templates(T, weights=weights),
                     names=("homolysis", "recombination")),
        700.0, np.random.default_rng(4))
    for chain in (0, 1):
        S.strip_h(engine.state, S.end_carbons(engine.state, chain)[0])
    engine.rebuild()
    prop = engine.propensities()
    direct = math.fsum(a for k, a in prop.channels
                       if k[0] == "met" or not engine.classes[k[1]].met_reverse)
    assert direct > 0.0 and prop.thinned > 1.0e15 * direct
    assert (prop.enabled - prop.thinned) < 0.5 * direct     # the subtraction has lost it
    horizon = 1.0e-18
    engine.run(until=horizon)
    assert engine.events_fired == 0
    assert engine.leak_estimate()["non_thinned_integral"] == pytest.approx(direct * horizon, rel=1e-9, abs=0.0)


def _balanced_leak_state(seed):
    """_leak_fraction_state with the leak weight set so that leak is about as large as the
    physical hazard (exact fraction near 0.5); the leak never fires, so dynamics are unchanged."""
    spec = S.pe_spec({12: 1}, d_ref=3.0e-20)
    engine = me.MeltEngine(
        spec, S.ToyCatalog(builder=lambda T: S.pe_templates(T, weights={"intra_mid_end": 0.0125}),
                           names=("homolysis", "recombination", "intra_mid_end")),
        700.0, np.random.default_rng(seed))
    state = engine.state
    c = S.carbons(state, 0)
    S.strip_h(state, sorted(c)[len(c) // 2])
    S.strip_h(state, S.end_carbons(state, 0)[0])
    engine.rebuild()
    return engine


def test_pooled_estimate_is_consistent_in_runs_and_a_single_run_is_a_biased_ratio():
    probe = _balanced_leak_state(0)
    prop = probe.propensities()
    exact = prop.leak / (prop.leak + prop.non_thinned + _physical_reverse_hazard(probe))
    assert 0.4 < exact < 0.6
    horizon = 1.0 / prop.thinned          # about one proposal per run
    reports, singles = [], []
    for seed in range(240):
        engine = _balanced_leak_state(5000 + seed)
        engine.run(until=horizon)
        reports.append(engine.report())
        singles.append(reports[-1]["leak_fraction"])
    pooled = me.pooled_leak_fraction(reports)
    assert pooled["runs"] == 240
    assert pooled["fraction"] == pytest.approx(exact, rel=0.05)
    # the per-run fraction is a ratio estimate: its expectation is not the exact value
    assert abs(np.mean(singles) - exact) > 0.05 * exact
    assert "high" not in me.LEAK_DEFINITION
    assert "pooled_leak_fraction" in me.LEAK_DEFINITION


def test_pooled_fraction_is_the_ratio_of_ensemble_sums():
    runs = [{"leak_estimate": {"leak_num": 1.0, "leak_den": 4.0}},
            {"leak_estimate": {"leak_num": 0.0, "leak_den": 2.0}},
            {"leak_estimate": {"leak_num": 3.0, "leak_den": 6.0}}]
    pooled = me.pooled_leak_fraction(runs)
    assert pooled["fraction"] == pytest.approx(4.0 / 12.0) and pooled["runs"] == 3
    assert pooled["leak_num"] == 4.0 and pooled["leak_den"] == 12.0


# --------------------------------------------------------------------------
# I-164: the fallback-phi cache, the phi-only product scan, the thinning report
# --------------------------------------------------------------------------
def _twin_homolysis(temperature):
    """The toy PE catalog plus a second, identically shaped homolysis class:
    two thinned unimolecular classes evaluated at the same state."""
    out = S.pe_templates(temperature)
    twin = next(t for t in S.pe_templates(temperature) if t.name == "homolysis")
    twin.name = "homolysis_twin"
    return out + [twin]


def _twin_engine(seed=1):
    catalog = S.ToyCatalog(builder=_twin_homolysis, names=("homolysis", "homolysis_twin", "recombination"))
    return me.MeltEngine(S.pe_spec({5: 4}), catalog, 700.0, np.random.default_rng(seed))


def _count_applies(engine):
    calls = []
    original = engine._begin_apply

    def counted(tr, view):
        calls.append(tr.class_index)
        return original(tr, view)

    engine._begin_apply = counted
    return calls


def _round_trip(engine):
    """Two commits that restore the graph: the revision moves by 2, the state
    hash and the index do not."""
    from rmgpy.kmc.melt_state import SetBond, SetRadical
    state = engine.state
    before = (state.revision, state.state_hash())
    carbon = S.end_carbons(state, 0)[0]
    h = S.strip_h(state, carbon)
    state.apply([SetBond(carbon, h, 0.0, 1.0), SetRadical(carbon, 1, 0), SetRadical(h, 1, 0)])
    assert state.revision == before[0] + 2 and state.state_hash() == before[1]


def test_fallback_phi_cache_holds_every_class_at_one_revision():
    """The cache is keyed per class: evaluating one thinned class does not
    evict another, so a repeated horizon view at an unchanged state re-applies
    nothing (the single-entry cache recomputed each class at every view)."""
    engine = _twin_engine()
    a, b = sorted(c.index for c in engine.classes if c.met_reverse)
    calls = _count_applies(engine)
    va = engine._fallback_phi(a)
    vb = engine._fallback_phi(b)
    assert calls and set(calls) == {a, b}
    n = len(calls)
    assert engine._fallback_phi(a) == va and engine._fallback_phi(b) == vb
    assert len(calls) == n
    assert va[0] == _exact_current_phi(engine, a) and vb[0] == _exact_current_phi(engine, b)


def test_fallback_phi_cache_is_keyed_on_the_state_revision():
    """A commit (even one that restores the graph) moves the revision and the
    cached value is recomputed; after a real change the served value is the
    exhaustive mean over the new state's transitions, never the old one."""
    engine = _twin_engine()
    a = min(c.index for c in engine.classes if c.met_reverse)
    calls = _count_applies(engine)
    v0 = engine._fallback_phi(a)
    n = len(calls)
    _round_trip(engine)
    assert engine._fallback_phi(a) == v0
    assert len(calls) == 2 * n                    # recomputed at the new revision
    S.strip_h(engine.state, S.end_carbons(engine.state, 1)[0])
    engine.rebuild()
    v1 = engine._fallback_phi(a)
    assert v1[0] == _exact_current_phi(engine, a)
    assert v1 != v0


def test_fallback_phi_cache_is_bounded():
    assert isinstance(me.PHI_FALLBACK_CACHE_LIMIT, int) and me.PHI_FALLBACK_CACHE_LIMIT >= 1
    engine = _twin_engine()
    a, b = sorted(c.index for c in engine.classes if c.met_reverse)
    limit = 3
    original = me.PHI_FALLBACK_CACHE_LIMIT
    me.PHI_FALLBACK_CACHE_LIMIT = limit
    try:
        sizes = []
        for _ in range(4):
            engine._fallback_phi(a)
            engine._fallback_phi(b)
            sizes.append(len(engine._phi_fallback))
            _round_trip(engine)
    finally:
        me.PHI_FALLBACK_CACHE_LIMIT = original
    assert max(sizes) == limit                       # the bound is reached and held
    calls = _count_applies(engine)
    engine._fallback_phi(a)
    engine._fallback_phi(a)
    assert calls and len(set(calls)) == 1


def test_phi_scan_gives_the_full_scan_pair_rate():
    """phi needs only the MET-forward sites, the MET lengths and the species of
    the two products: the restricted scan gives the same PairRate, bit for bit,
    for every current transition of every thinned class, and scans no
    unimolecular class."""
    engine = _twin_engine()
    S.strip_h(engine.state, S.end_carbons(engine.state, 2)[0])
    engine.rebuild()
    checked = 0
    for ci in (c.index for c in engine.classes if c.met_reverse):
        for c in sorted(engine.entries):
            for tr in engine.entries[c].uni.get(ci, ((), 0))[0]:
                tx, products = engine._begin_apply(tr, me._LocalView(engine.state, list(tr.components)))
                try:
                    full = engine._scan(products)
                    part = engine._phi_scan(products)
                finally:
                    tx.rollback()
                assert set(part) == set(full) == set(products)
                for p in products:
                    assert part[p].uni == {}
                    assert (part[p].met_length, part[p].species) == (full[p].met_length, full[p].species)
                assert engine.pair_rate(part[products[0]], part[products[1]]) == \
                    engine.pair_rate(full[products[0]], full[products[1]])
                checked += 1
    assert checked > 20


def test_report_decomposes_the_thinned_acceptance():
    """The run report states, per thinned class and overall, the proposals,
    the accepted ones, the summed escape factor and its mean (the escape
    probability per proposal), so zero events read as a probability."""
    engine = me.MeltEngine(S.pe_spec({10: 2, 20: 2}), S.ToyCatalog(names=("homolysis", "recombination")),
                           700.0, np.random.default_rng(3))
    engine.run(until=5.0, max_events=400)
    assert engine.thinned_proposals > 0
    th = engine.report()["thinning"]
    assert th["proposals"] == engine.thinned_proposals
    assert th["phi_sum"] == engine.phi_sum
    assert th["escape_probability"] == engine.phi_sum / engine.thinned_proposals
    # phi_sum is the compensator of the accepted count, not a conditional expectation
    assert "expected_accepted" not in th
    assert "compensator" in th["phi_sum_definition"] and "event cap" in th["phi_sum_definition"]
    fired = Counter(e.template for e in engine.events)
    names = {c.name for c in engine.classes if c.met_reverse}
    assert {row["template"] for row in th["classes"]} == names
    for row in th["classes"]:
        ci = next(c.index for c in engine.classes if c.name == row["template"])
        assert row["proposals"] == engine._proposals_by_class[ci]
        assert row["phi_sum"] == engine._phi_by_class[ci]
        assert row["accepted"] == fired[row["template"]]
        assert row["escape_probability"] == (row["phi_sum"] / row["proposals"] if row["proposals"] else None)
    assert th["accepted"] == sum(fired[n] for n in names)


# --------------------------------------------------------------------------
# I-164 rework 1: every key component and invalidation of the fallback-phi
# cache is pinned by a case in which dropping it serves a stale value
# --------------------------------------------------------------------------
def _mid_strip_engine():
    """An engine indexed with one mid-chain carbon radical (its two C-C bonds
    are no homolysis site), whose state then re-attaches the H outside the
    engine: one commit, after which the index is stale."""
    from rmgpy.kmc.melt_state import SetBond, SetRadical
    engine = me.MeltEngine(S.pe_spec({5: 4}), S.ToyCatalog(names=("homolysis", "recombination")), 700.0,
                           np.random.default_rng(1))
    state = engine.state
    mid = [c for c in S.carbons(state, 0)
           if sum(1 for b in state.atoms[c].edges if state.atoms[b].element == "C") == 2][1]
    h = S.strip_h(state, mid)
    engine.rebuild()
    old = sorted({state.components.component_of(mid), state.components.component_of(h)})
    state.apply([SetBond(mid, h, 0.0, 1.0), SetRadical(mid, 1, 0), SetRadical(h, 1, 0)])
    added = [state.components.component_of(mid)]
    return engine, old, added


@pytest.mark.parametrize("refresh", ["rebuild", "update"])
def test_fallback_phi_cache_is_keyed_on_the_index_revision(refresh):
    """At one state revision the index changes (rebuild, or _update of the
    changed components): the cached value of the old index is not served."""
    engine, old, added = _mid_strip_engine()
    ci = min(c.index for c in engine.classes if c.met_reverse)
    revision = engine.state.revision
    stale = engine._fallback_phi(ci)
    assert stale[0] == _exact_current_phi(engine, ci)       # the old index's own value
    if refresh == "rebuild":
        engine.rebuild()
    else:
        engine._update(old, added)
    assert engine.state.revision == revision
    fresh = _exact_current_phi(engine, ci)
    assert fresh != stale[0]                                # the index change moves the value
    assert engine._fallback_phi(ci)[0] == fresh


def test_update_bumps_the_index_revision():
    engine, old, added = _mid_strip_engine()
    before = engine._index_revision
    engine._update(old, added)
    assert engine._index_revision == before + 1


def _transport_engine(temperature=700.0):
    """Constant template weights (the toy catalog ignores T) and a transport
    with an activation energy: only k_D, hence phi, moves with temperature."""
    met = dict(S.pe_met(), segment_diffusivity={"law": "arrhenius", "D_ref_m2_s": 1.0e-17, "T_ref_K": 700.0,
                                                "Ea_J_mol": 6.0e4})
    engine = me.MeltEngine(S.pe_spec({5: 4}, met_block=met), S.ToyCatalog(names=("homolysis", "recombination")),
                           temperature, np.random.default_rng(1), temperature_range=(600.0, 900.0))
    S.strip_h(engine.state, S.end_carbons(engine.state, 0)[0])
    engine.rebuild()
    return engine


def test_fallback_phi_cache_is_keyed_on_the_temperature():
    engine = _transport_engine()
    ci = min(c.index for c in engine.classes if c.met_reverse)
    weights, index, revision = engine.weights, engine._index_revision, engine.state.revision
    at_700 = engine._fallback_phi(ci)
    engine.set_temperature(840.0)
    assert (engine.weights, engine._index_revision, engine.state.revision) == (weights, index, revision)
    fresh = _exact_current_phi(engine, ci)
    assert fresh > at_700[0] + 0.1                         # faster diffusion, more escape
    assert engine._fallback_phi(ci)[0] == fresh
    assert fresh == pytest.approx(_exact_current_phi(_transport_engine(840.0), ci), rel=1e-12)


def _strip_last_h(state, carbon):
    from rmgpy.kmc.melt_state import SetBond, SetRadical
    h = max(x for x in state.atoms[carbon].edges if state.atoms[x].element == "H")
    state.apply([SetBond(carbon, h, 1.0, 0.0), SetRadical(carbon, 0, 1), SetRadical(h, 0, 1)])


def test_fallback_phi_cache_belongs_to_one_state_object():
    """The engine given another state object at the same revision number (and
    an unchanged index revision) does not serve the first state's value."""
    def make():
        return me.MeltEngine(S.pe_spec({5: 2, 8: 2}), S.ToyCatalog(names=("homolysis", "recombination")),
                             700.0, np.random.default_rng(1))

    engine, other = make(), make()
    _strip_last_h(engine.state, S.end_carbons(engine.state, 0)[0])
    engine.rebuild()
    _strip_last_h(other.state, S.end_carbons(other.state, 3)[0])
    ci = min(c.index for c in engine.classes if c.met_reverse)
    first = engine._fallback_phi(ci)
    index = engine._index_revision
    engine.state = other.state
    assert engine.state.revision == 1 and engine._index_revision == index
    fresh = _exact_current_phi(engine, ci)
    assert fresh != first[0]
    assert engine._fallback_phi(ci)[0] == fresh
