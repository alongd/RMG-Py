"""The contract of a same-molecule law on a network (I-297, design only).

A network declaration (``NetworkPolymerSpec``) measures a bimolecular pair whose
two sites lie on one molecule as leak (``same_molecule``) and never fires it.
On a gel holding radicals that pair is the only chemistry left, so the run is
absorbed by leak.  The design (``/home/alon/runs/pk-i297/design.md``, revision
2) lists the candidate laws and leaves the choice to the owner.  The tests
assert only observables the engine already exposes (propensities, leak
classes, the event log, the MET counters) and fall in two groups.

Implementation invariants (any law):

1. a gel holding radicals evolves: the same-molecule pair fires;
2. classification: no same-molecule pair is left in ``same_molecule``; a pair
   closing a cycle shorter than ``minimum_cycle_bonds`` is ``short_cycle``
   leak and never fires; a pair whose product a loaded unimolecular closure
   already makes is that closure's (counted once);
3. the no-split reverse: a homolysis that leaves both radicals on one molecule
   is accepted at exactly the escape factor of the same-molecule pair it
   makes (R2), measured on that transition alone;
4. the fired rate is the reported rate (first-event waiting time);
5. the rate is a graph invariant (atom renumbering);
6. a far pair's excess over the two-molecule rate does not grow with distance.

Physical assumptions (marked ``law(...)`` with the options of design §3 that
make them, or ``phi_closure(...)`` with the escape-factor closures of design
§2; with ``I297_LAW`` / ``I297_PHI`` set to the chosen option the others
skip):

* the far pair fires at the two-molecule rate (a background term: A, C, D);
* replicating the box leaves the per-molecule rate unchanged (A, C, D);
* Kolmogorov's criterion on a closed reaction cycle (A, B_R);
* the pair rate is the option's own formula (one parameter per option);
* the no-split reverse's escape factor is the chosen closure's (self, seg).

The contract tests (``contract``) are skipped until the law is implemented
(set ``I297_RUN_RED=1`` to run them red against the engine as it is: each
must then fail on its contract assertion, not on its setup).  An unknown or
empty ``I297_LAW`` or ``I297_PHI`` is a collection error, never a silent
skip of every option.  ``declare_law`` (a declaration field, if the chosen
law needs one), ``declared_law`` (reads it back), ``exact_total`` (if the
law fires by thinning) and ``class_rate`` (a class's physical rate) are the
hooks the implementation edits.

The groundwork tests (I-301: the encounter group, the one escape-factor
evaluator, the phi_self mapping) assert today's engine and always run."""

import heapq
import math
import os

import numpy as np
import pytest

import melt_engine_support as S
from meltEngineTest import _radical_window, _ring_members
from network_support import network_spec, no_filler

from rmgpy import constants
from rmgpy.kmc import melt_engine as me

contract = pytest.mark.skipif(
    os.environ.get("I297_RUN_RED") != "1",
    reason="I-297 same-molecule law: design only (/home/alon/runs/pk-i297/design.md); these contract tests "
           "turn on with the owner's chosen law")

# the options of design §3 (revision 2) and the phi closures of §2; the implementation sets the chosen ones
LAWS = ("A", "B_sp", "B_R", "C_sp", "C_R", "C_img", "D")
PHI_CLOSURES = ("self", "seg")


def _selection(variable, options):
    """The option ``variable`` selects (None: unset, every option's tests
    run).  Set to anything not in ``options``, the empty string included, it
    is refused here, at collection: a typo must not skip every option."""
    value = os.environ.get(variable)
    if value is not None and value not in options:
        raise pytest.UsageError(f"{variable}={value!r} is not one of {options} (unset it to run every option)")
    return value


CHOSEN_LAW = _selection("I297_LAW", LAWS)
CHOSEN_PHI = _selection("I297_PHI", PHI_CLOSURES)


def law(*options):
    """A physical-assumption test: it binds only the options that make the
    assumption (all run while no option is chosen)."""
    assert set(options) <= set(LAWS), options
    return pytest.mark.skipif(CHOSEN_LAW is not None and CHOSEN_LAW not in options,
                              reason=f"a physical assumption of options {options}, not of {CHOSEN_LAW}")


def phi_closure(*options):
    """As ``law``, for the escape-factor closure of a same-molecule MET pair."""
    assert set(options) <= set(PHI_CLOSURES), options
    return pytest.mark.skipif(CHOSEN_PHI is not None and CHOSEN_PHI not in options,
                              reason=f"a phi closure of {options}, not of {CHOSEN_PHI}")


def declare_law(mapping):
    """The owner's chosen same-molecule law, as the implementation declares it
    on the network (no-op while the law is a design)."""
    return mapping


def declared_law(engine):
    """(law, phi closure) the engine fires, read from its declaration; (None,
    None) while no law is declared.  The implementation points this at the
    declaration field ``declare_law`` writes."""
    return None, None


def exact_total(engine):
    """The engine's exact total propensity (fired channels).  An implementation
    that fires same-molecule pairs by bounded thinning (design §5) reports an
    envelope in ``propensities()``; it points this hook at its exact
    evaluation."""
    return engine.propensities().total


def class_rate(engine, name):
    """The physical rate (1/s) of class ``name``'s fired transitions in the
    current state: its unimolecular and two-molecule channels, a thinned
    reverse class's times its exact mean escape factor (``_fallback_phi``,
    exhaustive below ``PHI_FALLBACK_LIMIT``), a MET class's share of every
    component pair's allowed rate.  Same-molecule pairs are leak today (0
    here); the implementation adds their fired rate."""
    ci = [t.name for t in engine.templates].index(name)
    cls = engine.classes[ci]
    rate = math.fsum(a for key, a in engine.propensities().channels if key[0] != "met" and key[1] == ci)
    if cls.met_reverse and rate:
        phi, sampled = engine._fallback_phi(ci)
        assert not sampled, "the class's mean phi is exact"
        rate *= phi
    rate += math.fsum(r.phi * engine._row_weight(row) * row[1] / engine.nav
                      for r in engine.pairs.values() for row in r.classes if row[0] == ci)
    return rate


def law_spec(*mutations, d_ref=1.0e-17):
    def mutate(mapping):
        mapping["transport_model"] = {"met": S.pe_met(d_ref)}
        no_filler(mapping)
        for m in mutations:
            m(mapping)
        declare_law(mapping)
    return network_spec(mutate)


def _engine(spec, state, templates, seed=3):
    return me.MeltEngine(spec, S.ToyCatalog(builder=lambda _t: list(templates)), 700.0,
                         np.random.default_rng(seed), state=state)


def _classes(engine):
    return dict(engine.propensities().leak_classes)


def _carbon_radicals(state):
    return sorted(a for a, n in state.atoms.items() if n.radical and n.element != "H")


# --------------------------------------------------------------------------
# hand-built windows (read from the state, so a non-match cannot pass silently)
# --------------------------------------------------------------------------
def _closed_window(state, closed_heavy, participant, radical=None):
    """``closed_heavy`` and their hydrogens closed, every other heavy
    neighbour a port whose exterior is read from the state, in the first
    resonance form (that keeps a radical on each atom of ``radical``, if given); returns
    the atoms, the induced bonds and the window index of each state atom."""
    from rmgpy.kmc import window_matcher as wm
    from rmgpy.kmc.window_matcher import WindowAtom
    cid = state.components.component_of(closed_heavy[0])
    rows = next(r for r in (wm.form_rows(state, cid, f) for f in wm.component_forms(state, cid))
                if radical is None or all(r[x][3] == 1 for x in (radical if isinstance(radical, (list, tuple))
                                                                   else (radical,))))
    ring = _ring_members(state)
    order = list(closed_heavy)
    order += [y for x in closed_heavy for y, _o in rows[x][5] if rows[y][0] == "H"]
    ports = sorted({y for x in closed_heavy for y, _o in rows[x][5] if y not in order})
    order += ports
    index = {x: i for i, x in enumerate(order)}
    atoms = []
    for x in order:
        element, _iso, charge, radical, lone_pairs, edges = rows[x]
        if x in ports:
            exterior = tuple((o, rows[y][0], y in ring) for y, o in edges if y not in index)
            atoms.append(WindowAtom(element, participant, "port", exterior, charge=charge, radical=radical,
                                    lone_pairs=lone_pairs))
        else:
            atoms.append(WindowAtom(element, participant, charge=charge, radical=radical, lone_pairs=lone_pairs))
    bonds = sorted({(min(index[x], index[y]), max(index[x], index[y]), o)
                    for x in order for y, o in rows[x][5] if y in index})
    return atoms, bonds, index


def _recombination(state, left, right, name, weight):
    """Radical ``left`` (participant 0) + radical ``right`` (participant 1) -> bond."""
    from rmgpy.kmc.window_matcher import AtomEdit, BondEdit, Template
    a0, b0 = _radical_window(state, left, 0)
    a1, b1 = _radical_window(state, right, 1)
    n = len(a0)
    return Template(name, "toy_recombination", a0 + a1, b0 + [(i + n, j + n, o) for i, j, o in b1],
                    [BondEdit(0, n, 1), AtomEdit(0, radical=0), AtomEdit(n, radical=0)], weight=weight)


def _wide_recombination(state, left, left_closed, right, right_closed, name, weight):
    """As ``_recombination``, with each radical's window closing
    ``*_closed`` (the radical first): a pair-specific template where the
    one-atom windows of two radicals coincide."""
    from rmgpy.kmc.window_matcher import AtomEdit, BondEdit, Template
    a0, b0, i0 = _closed_window(state, [left] + list(left_closed), 0, radical=left)
    a1, b1, i1 = _closed_window(state, [right] + list(right_closed), 1, radical=right)
    n = len(a0)
    return Template(name, "toy_recombination", a0 + a1, b0 + [(i + n, j + n, o) for i, j, o in b1],
                    [BondEdit(i0[left], i1[right] + n, 1), AtomEdit(i0[left], radical=0),
                     AtomEdit(i1[right] + n, radical=0)], weight=weight)


def _abstraction(state, radical, donor_closed, donor_carbon, name, weight):
    """Radical ``radical`` (participant 0) takes an H of ``donor_carbon``
    (participant 1, whose window closes ``donor_closed``): the radical count
    is conserved, so an activation-controlled (non-MET) class."""
    from rmgpy.kmc.window_matcher import AtomEdit, BondEdit, Template
    a0, b0 = _radical_window(state, radical, 0)
    a1, b1, index = _closed_window(state, donor_closed, 1)
    n = len(a0)
    c = index[donor_carbon]
    h = next(i for i, atom in enumerate(a1)
             if atom.element == "H" and any({i, c} == {p, q} for p, q, _o in b1))
    return Template(name, "toy_abstraction", a0 + a1, b0 + [(i + n, j + n, o) for i, j, o in b1],
                    [BondEdit(c + n, h + n, 0), BondEdit(0, h + n, 1), AtomEdit(0, radical=0),
                     AtomEdit(c + n, radical=1)], weight=weight)


def _homolysis(state, a, b, name, weight):
    """Bond ``a``-``b`` -> two radicals (unimolecular, radical-increasing:
    the reverse of a MET class, fired by thinning)."""
    from rmgpy.kmc.window_matcher import AtomEdit, BondEdit, Template
    atoms, bonds, index = _closed_window(state, [a, b], 0)
    return Template(name, "toy_homolysis", atoms, bonds,
                    [BondEdit(index[a], index[b], 0), AtomEdit(index[a], radical=1),
                     AtomEdit(index[b], radical=1)], weight=weight)


# --------------------------------------------------------------------------
# strand metric and the options' effective concentrations (design §2)
# --------------------------------------------------------------------------
_MASS = {"C": 12.011, "H": 1.008, "O": 15.999}


def _strand_graph(state, cid):
    """Heavy atoms of component ``cid`` and their bonds, each weighted by the
    mean united-atom mass (atom + its bonded H) of its two ends."""
    heavy = [a for a, n in state.atoms.items()
             if n.element != "H" and state.components.component_of(a) == cid]
    m = {a: _MASS[state.atoms[a].element]
         + sum(_MASS["H"] for b in state.atoms[a].edges if state.atoms[b].element == "H") for a in heavy}
    w = {(a, b): 0.5 * (m[a] + m[b]) for a in heavy for b in state.atoms[a].edges if b in m and a < b}
    return heavy, m, w


def _strand_mass(heavy, w, s, t):
    """The least strand mass between anchors ``s`` and ``t`` (a value: ties in the path are immaterial)."""
    adj = {a: [] for a in heavy}
    for (a, b), x in w.items():
        adj[a].append((b, x))
        adj[b].append((a, x))
    dist, queue = {s: 0.0}, [(0.0, s)]
    while queue:
        d, x = heapq.heappop(queue)
        if x == t:
            return d
        if d > dist[x]:
            continue
        for y, c in adj[x]:
            if d + c < dist.get(y, math.inf):
                dist[y] = d + c
                heapq.heappush(queue, (d + c, y))
    raise AssertionError("the anchors are not connected")


def _resistance(heavy, w, s, t):
    """Resistance distance between ``s`` and ``t`` with each bond a resistor of its strand mass."""
    index = {a: i for i, a in enumerate(heavy)}
    lap = np.zeros((len(heavy), len(heavy)))
    for (a, b), x in w.items():
        i, j = index[a], index[b]
        lap[i, i] += 1 / x
        lap[j, j] += 1 / x
        lap[i, j] -= 1 / x
        lap[j, i] -= 1 / x
    p = np.linalg.pinv(lap)
    i, j = index[s], index[t]
    return p[i, i] + p[j, j] - 2 * p[i, j]


def _gauss(r2):
    """Contact concentration (mol m^-3) of a Gaussian end-to-end vector of mean square ``r2``."""
    return (3.0 / (2.0 * math.pi * r2)) ** 1.5 / constants.Na


def _image_sum(r2, box):
    """The same, normalised in a periodic cube of side ``box`` (sum over images)."""
    s2 = r2 / 3.0
    if s2 < (box / 3.0) ** 2:
        n = range(-3, 4)
        tot = sum(math.exp(-(i * i + j * j + k * k) * box * box / (2 * s2)) for i in n for j in n for k in n)
        return tot * (2 * math.pi * s2) ** -1.5 / constants.Na
    q = range(-6, 7)
    tot = sum(math.exp(-(2 * math.pi / box) ** 2 * (i * i + j * j + k * k) * s2 / 2)
              for i in q for j in q for k in q)
    return tot / box ** 3 / constants.Na


def _c_eff(option, engine, s, t):
    """The effective concentration of anchors ``s``, ``t`` (one component) under ``option``."""
    state = engine.state
    heavy, m, w = _strand_graph(state, state.components.component_of(s))
    r2pm = engine.transport.r2_per_molar_mass
    c_bg = 1.0 / engine.nav
    if option == "A":
        return c_bg
    sp = _gauss(r2pm * _strand_mass(heavy, w, s, t))
    res = _gauss(r2pm * (_resistance(heavy, w, s, t) + 0.5 * (m[s] + m[t])))
    return {"B_sp": sp, "B_R": res, "C_sp": c_bg + sp, "C_R": c_bg + res, "D": c_bg + sp,
            "C_img": _image_sum(r2pm * _strand_mass(heavy, w, s, t), engine.volume ** (1.0 / 3.0))}[option]


# --------------------------------------------------------------------------
# a far pair on one chain (activation controlled: no escape factor)
# --------------------------------------------------------------------------
def _far_chain(rings, copies=1, reverse=False):
    """Ph0-M-Ph1-M-...-Ph(rings-1), each M joining o6 of the ring before to
    o2 of the ring after, ``copies`` times; Ph0's para H stripped on every
    copy (an aryl radical at one end).  ``reverse`` numbers the units from the
    other end (the same graph, other atom ids).  Returns the spec, the state
    and (radical, donor carbon, far-ring atoms) of copy 0."""
    from rmgpy.kmc import melt_branched as mb
    spec = law_spec()
    n = 2 * rings - 1
    unit = (lambda u: n - 1 - u) if reverse else (lambda u: u)
    names = ["Ph" if unit(u) % 2 == 0 else "M" for u in range(n)]
    links = []
    for k in range(rings - 1):
        links += [(unit(2 * k + 1), "a", unit(2 * k), "o6"), (unit(2 * k + 1), "b", unit(2 * k + 2), "o2")]
    state = mb.build_molecules(spec, [(names, links)] * copies)
    atom = mb.residue_atoms(state)
    for c in range(copies):
        S.strip_h(state, atom[(c, unit(0))]["r4"])
    last = 2 * (rings - 1)
    donor = atom[(0, unit(last - 1))]["m"]
    far = [atom[(0, unit(last))][k] for k in ("r1", "r2", "r3", "r4", "r5", "r6", "o")]
    return spec, state, (atom[(0, unit(0))]["r4"], donor, far)


def _far_pair(rings, copies=1, reverse=False, seed=3):
    """The far chain with an abstraction of the far-end methylene's H by the
    aryl radical, the donor window closing the far end ring (so no other
    methylene matches).  Returns the engine and the two-molecule law's
    propensity of one copy's own pair, ``w n / (divisor N_A V)`` (n = 2: the
    donor's two H), in the volume of ``copies`` copies."""
    spec, state, (radical, donor, far) = _far_chain(rings, copies, reverse)
    template = _abstraction(state, radical, [donor] + far, donor, "far_h", 1.0e-2)
    volume = me.declared_initialization(spec).volume_m3 * copies
    engine = me.MeltEngine(spec, S.ToyCatalog(builder=lambda _t: [template]), 700.0,
                           np.random.default_rng(seed), state=state, volume_m3=volume)
    ci = [t.name for t in engine.templates].index("far_h")
    engine._i297_anchors = (radical, donor)
    return engine, engine.weights[ci] * 2 / (engine.classes[ci].divisor * engine.nav)


@contract
@law("A", "C_sp", "C_R", "C_img", "D")
@pytest.mark.parametrize("rings", [6, 24])
def test_a_far_pair_on_one_molecule_fires_at_the_two_molecule_rate(rings):
    """A background term (options A, C, D; B has none): the pair's propensity
    is never below the two-molecule law's and is within 10% of it at 24 rings
    (strand 1.14e3 g/mol: the loop term is 5.0% (C_sp, D), 6.3% (C_R) and
    0% (C_img) of the background there; probe_oracle)."""
    engine, expected = _far_pair(rings)
    prop = engine.propensities()
    assert _classes(engine)["same_molecule"] == 0.0
    assert prop.leak == 0.0
    ratio = exact_total(engine) / expected
    assert ratio >= 1.0 - 1.0e-9
    if rings == 24:
        assert ratio < 1.10


@contract
def test_the_far_pair_excess_over_the_two_molecule_rate_does_not_grow_with_distance():
    """Invariant: every option's effective concentration is non-increasing in strand length."""
    (near, e_near), (far, e_far) = _far_pair(6), _far_pair(24)
    r_near, r_far = exact_total(near) / e_near, exact_total(far) / e_far
    assert r_near > 0.0 and r_far > 0.0
    assert r_far - 1.0 <= r_near - 1.0 + 1.0e-12


@contract
def test_the_pair_rate_is_a_graph_invariant():
    """Invariant: numbering the units from the other end (other atom ids, the
    same graph) leaves the pair's propensity unchanged."""
    forward, _e = _far_pair(6)
    backward, _e = _far_pair(6, reverse=True)
    a, b = exact_total(forward), exact_total(backward)
    assert a > 0.0
    assert b == pytest.approx(a, rel=1.0e-12, abs=0.0)


@contract
def test_the_first_event_waits_one_over_the_reported_propensity():
    """Invariant (a waiting-time oracle): with the far pair the only channel,
    the mean first-event time over 100 seeds is 1/P within 4 standard
    errors (an exponential's: 1/sqrt(100)), P the reported propensity."""
    times, rates = [], []
    for seed in range(100):
        engine, _e = _far_pair(6, seed=seed)
        prop = engine.propensities()
        assert prop.leak == 0.0
        rates.append(exact_total(engine))
        assert exact_total(engine) > 0.0, "the far pair is a channel"
        engine.run(max_events=1)
        assert engine.events_fired == 1 and engine.template_counts["far_h"] == 1
        times.append(engine.events[0].time)
    assert max(rates) == pytest.approx(min(rates), rel=1.0e-12)
    assert abs(np.mean(times) * rates[0] - 1.0) < 4.0 / math.sqrt(len(times))


@contract
@pytest.mark.parametrize("option", [pytest.param(o, marks=law(o)) for o in LAWS])
def test_the_pair_rate_is_the_options_formula(option):
    """A quantitative oracle per option (a constant background passes the
    far-rate tests; it fails here unless chosen): the 6-ring pair's
    propensity is ``w n c_eff / divisor`` with c_eff from the declared
    <R^2>/M, the least united-atom strand mass (or the resistance distance
    plus the formed contact's weight) between the radical and the donor
    carbon, and the declared volume (design §2-§3)."""
    engine, expected = _far_pair(6)
    s, t = engine._i297_anchors
    want = expected * engine.nav * _c_eff(option, engine, s, t)
    assert want > 0.0
    assert exact_total(engine) == pytest.approx(want, rel=1.0e-9)


@contract
@law("A", "C_sp", "C_R", "D")
def test_replicating_the_box_leaves_the_per_molecule_rate_unchanged():
    """k copies of the far chain in k times the volume: one copy's own pair
    plus its pairs with the other copies fire at the same total per copy.
    A, C_sp, C_R and D hold this for an activation-controlled class (their
    background is the replicas' share, design §3); B does not (at k = 1 the
    box has no replica: +c_bg/2 per copy at k = 2), nor C_img (its own
    term depends on the box side: +10% at k = 2, 6 rings)."""
    one, _e = _far_pair(6, copies=1)
    two, _e = _far_pair(6, copies=2)
    p1, p2 = one.propensities(), two.propensities()
    assert p1.leak == 0.0 and p2.leak == 0.0
    assert exact_total(one) > 0.0
    assert exact_total(two) / 2 == pytest.approx(exact_total(one), rel=1.0e-9)


# --------------------------------------------------------------------------
# 1. a gel holding radicals evolves
# --------------------------------------------------------------------------
def _cut_fixture(spec, cuts=2):
    """The network fixture with ``cuts`` methylene-ring bridges cut (as
    pk-i243 §9): a benzylic and an aryl radical per cut."""
    from rmgpy.kmc import melt_branched as mb
    from rmgpy.kmc.melt_state import SetBond, SetRadical, materialize
    state = materialize(spec)
    atom = mb.residue_atoms(state)
    out = []
    for (_mol, _unit), keys in sorted(atom.items()):
        if "m" in keys and len(out) < cuts:
            m = keys["m"]
            ring = next(b for b in sorted(state.atoms[m].edges) if state.atoms[b].element == "C")
            state.apply([SetBond(m, ring, 1.0, 0.0), SetRadical(m, 0, 1), SetRadical(ring, 0, 1)])
            out.append((m, ring))
    return state, out


@contract
def test_a_gel_holding_radicals_evolves():
    """Without a law the run stalls after one recombination with both remaining
    radicals on the gel (pk-i243 §9); with it a same-molecule pair fires (on
    the gel it closes a 20-bond cycle, above the declared 8) and no radical
    is left."""
    spec = law_spec()
    state, cuts = _cut_fixture(spec)
    engine = _engine(spec, state, [_recombination(state, cuts[0][0], cuts[0][1], "bridge_ring", 1.0e-2)],
                     seed=1)
    engine.run(until=1.0e4, max_events=10)
    assert engine.template_counts.get("bridge_ring", 0) == 2
    assert _carbon_radicals(engine.state) == []
    # radicals m0, m3 (benzylic) and r98, r170 (aryl): 0, 3 and 170 start on one sol fragment, 98 on
    # the gel, so every complete pairing has one recombination within one molecule (first or second)
    assert any(len(set(e.components)) == 1 for e in engine.events), "a recombination within one molecule"
    assert engine.leak_by_class["same_molecule"] == 0.0


# --------------------------------------------------------------------------
# 2. classification
# --------------------------------------------------------------------------
def _melt(*strip):
    """networkMeltTest's melt (mol0 Ph0-M1-Ph2-M3 with a bridge radical on M3,
    mol1 Ph-M-Ph, an aryl radical at Ph0's para site on both) under the law,
    plus a para radical on mol0 at each unit in ``strip``."""
    from networkMeltTest import _melt as network_melt
    spec, state, atom = network_melt(lambda mutate: law_spec(mutate))
    for unit in strip:
        S.strip_h(state, atom[(0, unit)]["r4"])
    return spec, state, atom


@contract
def test_no_same_molecule_pair_is_left_as_same_molecule_leak():
    """The bridge radical and Ph0's aryl radical of one molecule (an 8-bond
    cycle: the declared minimum) are a channel, not leak."""
    spec, state, atom = _melt()
    bridge, mine = atom[(0, 3)]["m"], atom[(0, 0)]["r4"]
    engine = _engine(spec, state, [_recombination(state, bridge, mine, "m_p", 1.0e-2)])
    prop = engine.propensities()
    classes = _classes(engine)
    assert classes["same_molecule"] == 0.0
    # a law may still measure what a declared approximation omits (gel_mobility), never same_molecule
    assert classes["undeclared_link"] == classes["short_cycle"] == 0.0
    inter = math.fsum(r.allowed for (a, b), r in engine.pairs.items() if a != b)
    assert exact_total(engine) > inter, "the same-molecule pair is a channel of its own"
    engine.run(until=1.0e3)
    assert engine.leak_by_class["same_molecule"] == 0.0
    assert math.fsum(engine.leak_by_class.values()) == pytest.approx(engine.leak_num)


@contract
def test_a_same_molecule_pair_closing_a_short_cycle_is_short_cycle_leak_and_never_fires():
    """Ph2's para radical closing on the bridge radical M3 of the same
    molecule: a 4-bond cycle, below the declared 8.  It is short_cycle leak,
    and over 12 seeds the first event never makes that bond."""
    fired_short = 0
    for seed in range(12):
        spec, state, atom = _melt(2)
        bridge, mine, short = atom[(0, 3)]["m"], atom[(0, 0)]["r4"], atom[(0, 2)]["r4"]
        engine = _engine(spec, state, [_recombination(state, bridge, mine, "m_p", 1.0e-2)], seed=seed)
        if seed == 0:
            classes = _classes(engine)
            assert classes["same_molecule"] == 0.0
            assert classes["short_cycle"] > 0.0
            assert classes["undeclared_link"] == 0.0
        engine.run(max_events=1)
        assert engine.events_fired == 1
        fired_short += short in engine.state.atoms[bridge].edges
    assert fired_short == 0


@contract
def test_a_unimolecular_closure_owns_the_pair_it_duplicates():
    """R3: a loaded unimolecular template that closes the same 8-bond ring
    (its window covers the strand) makes exactly the product of the
    same-molecule pair; the pair is that closure's and is counted once
    (never as a second channel, never as same_molecule leak)."""
    from collections import deque
    from rmgpy.kmc.window_matcher import AtomEdit, BondEdit, Template
    spec, state, atom = _melt()
    bridge, mine = atom[(0, 3)]["m"], atom[(0, 0)]["r4"]
    prev, queue = {bridge: None}, deque([bridge])
    while queue:
        x = queue.popleft()
        for y in sorted(state.atoms[x].edges):
            if y not in prev and state.atoms[y].element != "H":
                prev[y] = x
                queue.append(y)
    path, x = [], mine
    while x is not None:
        path.append(x)
        x = prev[x]
    atoms, bonds, index = _closed_window(state, path, 0, radical=(bridge, mine))
    uni = Template("ring_close", "toy_intra_recombination", atoms, bonds,
                   [BondEdit(index[bridge], index[mine], 1), AtomEdit(index[bridge], radical=0),
                    AtomEdit(index[mine], radical=0)], weight=1.0e3)
    alone = _engine(spec, state, [uni])
    p_uni = alone.propensities()
    assert p_uni.leak == 0.0 and exact_total(alone) > 0.0, "setup: the closure is a unimolecular channel"
    both = _engine(spec, state, [uni, _recombination(state, bridge, mine, "m_p", 1.0e-2)])
    prop = both.propensities()
    assert _classes(both)["same_molecule"] == 0.0
    inter = math.fsum(r.allowed for (a, b), r in both.pairs.items() if a != b)
    assert exact_total(both) - inter == pytest.approx(exact_total(alone), rel=1.0e-12)


# --------------------------------------------------------------------------
# 3. the no-split reverse (R2)
# --------------------------------------------------------------------------
_RING = [(["Ph", "M", "Ph", "M"], [(1, "a", 0, "o2"), (1, "b", 2, "o2"), (3, "a", 2, "o6"), (3, "b", 0, "o6")])]


def _opened_ring(d_ref):
    """The ring Ph0-M1-Ph2-M3 (one 8-bond cycle) opened at M1-Ph0: one
    molecule with a methylene and an aryl radical; the engine loads only
    their recombination."""
    from rmgpy.kmc import melt_branched as mb
    from rmgpy.kmc.melt_state import SetBond, SetRadical
    spec = law_spec(d_ref=d_ref)
    cut = mb.build_molecules(spec, _RING)
    atom = mb.residue_atoms(cut)
    m, r = atom[(0, 1)]["m"], atom[(0, 0)]["r2"]
    cut.apply([SetBond(m, r, 1.0, 0.0), SetRadical(m, 0, 1), SetRadical(r, 0, 1)])
    return spec, cut, m, r


_NO_SPLIT_D_REF = 1.0e-17   # declared D0 (m2/s) at 700 K (S.pe_met: Arrhenius, Ea 0, T_ref 700 K)
_NO_SPLIT_W = 5.0e-2        # the recombination's weight: k_act = w (one ordered site pair, divisor 1)


def _analytic_k_d(closure, transport, met_length):
    """k_D (m3/mol/s) of a same-molecule MET pair under ``closure``, from the
    declared MET block alone (design §2; no engine rate is read):
    ``self``, ``k_diffusion(T, L, L)`` written out (Rouse D(L) = D0/L, radius
    max(sigma, 2 Rg(L)), Rg^2 = <R^2>/M * M0 * L / 6); ``seg``, the segment
    closure with design §2's recommended radius: D(1) = D0, radius sigma."""
    d0 = _NO_SPLIT_D_REF
    sigma, spin = transport.contact_radius_m, transport.spin_factor
    if closure == "self":
        rg = math.sqrt(transport.r2_per_molar_mass * transport.repeat_unit_molar_mass_g_mol * met_length / 6.0)
        return 4.0 * math.pi * constants.Na * spin * (2.0 * d0 / met_length) * max(sigma, 2.0 * rg)
    return 4.0 * math.pi * constants.Na * spin * (2.0 * d0) * sigma


@contract
@pytest.mark.parametrize("closure", [pytest.param(c, marks=phi_closure(c)) for c in PHI_CLOSURES])
def test_the_no_split_reverse_is_accepted_at_the_same_molecule_pairs_escape_factor(closure):
    """R2: every methylene-ring homolysis of the ring leaves both radicals on
    one molecule (no split), and the four are equivalent (one product species),
    so every proposal before the first event has one acceptance factor and
    ``phi_sum / thinned_proposals`` is it.  It must be the escape factor of
    the same-molecule pair the homolysis makes, k_D / (k_D + k_act), from an
    analytic oracle: k_act = w (the methylene radical as participant 0, the
    aryl radical as participant 1: one ordered site pair, divisor 1) and k_D
    the chosen closure's (``_analytic_k_d``).  The declared D0 and w put phi
    at 0.5-0.6 under either closure, so the run needs a few proposals."""
    from rmgpy.kmc import melt_branched as mb
    spec, cut, m, r = _opened_ring(_NO_SPLIT_D_REF)
    state = mb.build_molecules(spec, _RING)
    engine = _engine(spec, state, [_homolysis(state, m, r, "open", 1.0),
                                   _recombination(cut, m, r, "close", _NO_SPLIT_W)])
    (entry,) = engine.entries.values()
    close = [t.name for t in engine.templates].index("close")
    assert engine.classes[close].divisor == 1, "setup: one ordered site pair is one pair"
    k_d = _analytic_k_d(closure, engine.transport, entry.met_length)
    phi = k_d / (k_d + _NO_SPLIT_W)
    assert 0.2 < phi < 0.8, f"setup: moderate acceptance ({phi})"
    engine.run(max_events=1)
    assert engine.template_counts.get("open", 0) == 1 and engine.events_fired == 1
    assert engine.thinned_proposals >= 1
    assert engine.phi_sum / engine.thinned_proposals == pytest.approx(phi, rel=1.0e-12)
    assert engine.met_counts.get("reverse_no_pair", 0) == 0


# --------------------------------------------------------------------------
# 4. Kolmogorov's criterion on closed reaction cycles (one catalog, finite phi)
# --------------------------------------------------------------------------
def _wide_homolysis(state, a, b, closed, name, weight):
    """As ``_homolysis``, the window closing ``closed`` as well (a bond told
    from the equivalent bonds of the molecule)."""
    from rmgpy.kmc.window_matcher import AtomEdit, BondEdit, Template
    atoms, bonds, index = _closed_window(state, [a, b] + list(closed), 0)
    return Template(name, "toy_homolysis", atoms, bonds,
                    [BondEdit(index[a], index[b], 0), AtomEdit(index[a], radical=1),
                     AtomEdit(index[b], radical=1)], weight=weight)


def _move_h(state, donor, acceptor):
    """An H of ``donor`` to the radical ``acceptor``: the radical moves."""
    from rmgpy.kmc.melt_state import SetBond, SetRadical
    h = min(b for b in state.atoms[donor].edges if state.atoms[b].element == "H")
    state.apply([SetBond(donor, h, 1.0, 0.0), SetBond(acceptor, h, 0.0, 1.0),
                 SetRadical(donor, 0, 1), SetRadical(acceptor, 1, 0)])


_H_SQUARE = [(["Ph", "M", "Ph", "M", "M", "Ph", "M", "Ph"],
              [(1, "a", 0, "o2"), (1, "b", 2, "o2"), (3, "a", 2, "o6"), (3, "b", 0, "o6"),
               (4, "a", 2, "p4"), (4, "b", 5, "o2"), (6, "a", 5, "o6"), (6, "b", 7, "o2")])]


def _h_square():
    """A ring Ph0-M1-Ph2-M3 (8 bonds) with a tail M4-Ph5-M6-Ph7 on Ph2's para
    site and an aryl radical X at Ph7's para site.  Four states: the ring
    closed or opened at M1-Ph0 (``1``/``0`` first), X's radical in place or
    moved by an H transfer from M4 to X (``0``/``1`` second).  One catalog:
    the ring closure (methylene + aryl radical, MET) and its homolysis; a
    competing MET closure of the methylene radical on X (its row is in the
    opened ring's encounter group only while X is a radical, so phi differs
    between the two opened states); X abstracting an H of M4 and the benzylic
    M4 radical abstracting it back from Ph7 (non-MET, both one molecule).
    The ring is mirror-symmetric (M1 <-> M3), so the homolysis has two
    equivalent bonds in each closed state, and M4 two equivalent H: the cycle
    is on the lumped (species) states, each multiplicity on both sides."""
    from rmgpy.kmc import melt_branched as mb
    from rmgpy.kmc.melt_state import SetBond, SetRadical
    spec = law_spec(d_ref=_NO_SPLIT_D_REF)
    states, at = {}, None
    for key in ("00", "10", "01", "11"):
        state = mb.build_molecules(spec, _H_SQUARE)
        atom = mb.residue_atoms(state)
        at = {"m": atom[(0, 1)]["m"], "r": atom[(0, 0)]["r2"], "x": atom[(0, 7)]["r4"], "d": atom[(0, 4)]["m"],
              "r1": atom[(0, 0)]["r1"], "r3": atom[(0, 0)]["r3"], "p2": atom[(0, 2)]["r2"],
              "p4": atom[(0, 2)]["r4"], "x3": atom[(0, 7)]["r3"], "x5": atom[(0, 7)]["r5"]}
        S.strip_h(state, at["x"])
        if key[0] == "0":
            state.apply([SetBond(at["m"], at["r"], 1.0, 0.0), SetRadical(at["m"], 0, 1),
                         SetRadical(at["r"], 0, 1)])
        if key[1] == "1":
            _move_h(state, at["d"], at["x"])
        states[key] = (spec, state)
    s00, s10, s01 = states["00"][1], states["10"][1], states["01"][1]
    templates = [
        _recombination(s00, at["m"], at["r"], "close", _NO_SPLIT_W),
        _recombination(s00, at["m"], at["x"], "close_x", 2.0 * _NO_SPLIT_W),
        _wide_homolysis(s10, at["m"], at["r"], [at["r1"], at["r3"], at["p2"]], "open", 1.0),
        _abstraction(s00, at["x"], [at["d"], at["p4"]], at["d"], "h_fwd", 1.0e-2),
        _abstraction(s01, at["d"], [at["x"], at["x3"], at["x5"]], at["x"], "h_rev", 1.0e-2),
    ]
    counts = {"00": {"close": 1, "close_x": 1, "h_fwd": 2}, "10": {"open": 2, "h_fwd": 2},
              "01": {"close": 1, "h_rev": 1}, "11": {"open": 2, "h_rev": 1}}
    cycle = [("00", "close", "10", "open"), ("10", "h_fwd", "11", "h_rev"),
             ("11", "open", "01", "close"), ("01", "h_rev", "00", "h_fwd")]
    return states, templates, counts, cycle


def _transition_counts(engine):
    """Class name -> its transitions in the engine's state: a unimolecular
    class's allowed transitions, a bimolecular class's ordered site pairs on
    one component (the measured same-molecule count, over the divisor) and
    across components."""
    names = [t.name for t in engine.templates]
    out = {}
    for ci, cls in enumerate(engine.classes):
        if cls.bimolecular:
            n = sum(v for (cj, _k0, _k1), v in engine.same.items() if cj == ci) // cls.divisor
            n += sum(round(row[1] + row[2]) for r in engine.pairs.values() for row in r.classes if row[0] == ci)
        else:
            n = sum(len(e.shape.uni.get(ci, ((), 0))[0]) for e in engine.entries.values())
        if n:
            out[names[ci]] = n
    return out


def _state_graph(state):
    import networkx as nx
    g = nx.Graph()
    for a, n in state.atoms.items():
        g.add_node(a, element=n.element, radical=n.radical)
    g.add_edges_from((a, b, {"order": o}) for a, n in state.atoms.items() for b, o in n.edges.items() if a < b)
    return g


def _isomorphic_states(g, h):
    import networkx as nx
    return nx.is_isomorphic(g, h, node_match=lambda x, y: x == y, edge_match=lambda x, y: x == y)


def _uni_products(engine, name):
    """The state graph each transition of unimolecular class ``name`` makes
    (applied tentatively, rolled back)."""
    ci = [t.name for t in engine.templates].index(name)
    out = []
    for cid in sorted(engine.entries):
        entry = engine.entries[cid]
        for k in range(len(entry.shape.uni.get(ci, ((), 0))[0])):
            tr = entry.uni_at(ci, k)
            tx, _products = engine._begin_apply(tr, me._LocalView(engine.state, list(tr.components)))
            try:
                out.append(_state_graph(engine.state))
            finally:
                tx.rollback()
    return out


@contract
@law("A", "B_R")
def test_the_closed_reaction_cycle_satisfies_kolmogorovs_criterion():
    """``_h_square``: one catalog holds every forward and reverse class of the
    cycle (closure and its homolysis, H transfer and its reverse), at a finite
    escape factor (declared D0 1e-17 m2/s).  Setup checks that each state has
    exactly the transitions ``counts`` lists and that every unimolecular
    transition of an edge makes the edge's other state (graph isomorphism),
    so a window that matches elsewhere fails setup, not the criterion.  The
    product of the physical rates (``class_rate``: a thinned reverse at its
    escape factor) along the cycle equals the product against it.  The
    competing closure's row is in one opened state's encounter group and not
    the other's, so this holds only when a reverse is thinned by the group the
    forward reads (thinned by an inverse-only group instead, the ratio is
    phi_00/phi_01 = 0.459 under phi_self: probe_square_oracle); c_eff here is a state function under every
    option (the H transfer is off the ring's strand), so the c_eff half is
    ``test_the_nested_closure_rates_are_path_independent``."""
    states, templates, counts, cycle = _h_square()
    engines = {key: _engine(spec, state, templates) for key, (spec, state) in states.items()}
    for key, engine in engines.items():
        assert _transition_counts(engine) == counts[key], f"setup: the transitions of state {key}"
    graphs = {key: _state_graph(state) for key, (_spec, state) in states.items()}
    for src, name, dst, rev in cycle:
        for at, cls, target in ((src, name, dst), (dst, rev, src)):
            engine = engines[at]
            if not engine.classes[[t.name for t in engine.templates].index(cls)].bimolecular:
                assert all(_isomorphic_states(g, graphs[target]) for g in _uni_products(engine, cls)), \
                    f"setup: every {cls} transition of state {at} makes state {target}"
    # each edge (src, class, dst, reverse class): along the cycle src -> dst, against it dst -> src
    forward = [class_rate(engines[src], name) for src, name, _dst, _rev in cycle]
    backward = [class_rate(engines[dst], rev) for _src, _name, dst, rev in cycle]
    assert all(r > 0.0 for r in forward + backward), "every edge of the cycle is a fired channel"
    assert math.prod(forward) == pytest.approx(math.prod(backward), rel=1.0e-9)


def _chords(closed, radicals):
    """Ph0-M1-Ph2-M3-Ph4-M5-Ph6 (o6 -> o2) with chord A (MA: Ph0 o2 -> Ph6 o6)
    and chord B (MB: Ph2 p4 -> Ph6 p4); the chords not in ``closed`` are
    opened at Ph6 (a methylene and an aryl radical each).  B shortens A's
    strand (187 -> 139 g/mol); A leaves B's unchanged (139 g/mol).  Returns
    the engine loading one recombination specific to ``radicals``' two ends
    (the methylene's window closes its ring anchor, so it tells MA from MB)."""
    from rmgpy.kmc import melt_branched as mb
    from rmgpy.kmc.melt_state import SetBond, SetRadical
    spec = law_spec(d_ref=1.0e-6)
    names, links = [], []
    for k in range(4):
        names.append("Ph")
        if k < 3:
            names.append("M")
            links += [(2 * k + 1, "a", 2 * k, "o6"), (2 * k + 1, "b", 2 * k + 2, "o2")]
    ma, mb_ = len(names), len(names) + 1
    names += ["M", "M"]
    links += [(ma, "a", 0, "o2"), (ma, "b", 6, "o6"), (mb_, "a", 2, "p4"), (mb_, "b", 6, "p4")]
    end = {"A": (ma, "r6", (0, "r2")), "B": (mb_, "r4", (2, "r4"))}
    state = mb.build_molecules(spec, [(names, links)])
    atom = mb.residue_atoms(state)
    for c in "AB":
        if c not in closed:
            m, r = atom[(0, end[c][0])]["m"], atom[(0, 6)][end[c][1]]
            state.apply([SetBond(m, r, 1.0, 0.0), SetRadical(m, 0, 1), SetRadical(r, 0, 1)])
    mu, key, (anchor_unit, anchor_key) = end[radicals]
    m, r = atom[(0, mu)]["m"], atom[(0, 6)][key]
    template = _wide_recombination(state, m, [atom[(0, anchor_unit)][anchor_key]], r, [], "close", 1.0e-2)
    return _engine(spec, state, [template])


@contract
@law("A", "B_R")
def test_the_nested_closure_rates_are_path_independent(monkeypatch):
    """The c_eff half of Kolmogorov's criterion, forward rates only (the
    reversible one-catalog check at finite phi is
    ``test_the_closed_reaction_cycle_satisfies_kolmogorovs_criterion``; this
    fixture's homolyses cannot be told apart by a window in every state, so
    each forward rate is read with its own one-template catalog at phi -> 1).
    G0 -A-> GA -B-> G2 and G0 -B-> GB -A-> G2: were the reverse homolyses
    state-independent unimolecular rates, detailed balance would need
    P_A(G0) P_B(GA) = P_B(G0) P_A(GB), i.e. c_eff a state function.  Options
    A and B_R hold it (design §3: tau(G+e) = tau(G)(1 + R/w_e)); B_sp, C_sp,
    C_R and C_img miss it by 0.641, 0.807, 1.063 and 0.816 on this cycle
    (probe_cycle_nested)."""
    from rmgpy.kmc import window_matcher as wm
    # benzylic radicals on a 4-ring chain exceed the matcher's 64 resonance forms
    monkeypatch.setattr(wm._component_subs, "__defaults__", (4096,))
    rate = {}
    for key, closed, radicals in (("A|G0", (), "A"), ("B|GA", ("A",), "B"),
                                  ("B|G0", (), "B"), ("A|GB", ("B",), "A")):
        engine = _chords(closed, radicals)
        assert engine.propensities().leak == 0.0
        assert exact_total(engine) > 0.0, f"{key}: the pair is a channel"
        rate[key] = exact_total(engine)
    assert rate["A|G0"] * rate["B|GA"] == pytest.approx(rate["B|G0"] * rate["A|GB"], rel=1.0e-9)


@contract
def test_the_selected_options_are_the_engines_declared_law():
    """The selection is validated against the declaration, not only against
    the option names: with ``I297_LAW`` / ``I297_PHI`` set, the engine fires
    that law and that closure (``declared_law``)."""
    if CHOSEN_LAW is None and CHOSEN_PHI is None:
        pytest.skip("no option selected (I297_LAW / I297_PHI unset)")
    engine, _e = _far_pair(6)
    law_, phi_ = declared_law(engine)
    if CHOSEN_LAW is not None:
        assert law_ == CHOSEN_LAW
    if CHOSEN_PHI is not None:
        assert phi_ == CHOSEN_PHI


# --------------------------------------------------------------------------
# groundwork (I-301): today's engine, always run
# --------------------------------------------------------------------------
def _chain_links(rings):
    return [x for k in range(rings - 1) for x in ((2 * k + 1, "a", 2 * k, "o6"), (2 * k + 1, "b", 2 * k + 2, "o2"))]


def _competing(bonded, excluded=()):
    """X = Ph-M-Ph with an aryl radical A (Ph0's para site); Y =
    Ph-M-Ph-M-Ph with an aryl radical C (Ph0's para site) and a benzylic
    radical B (M3).  Two MET classes pair X with Y and make different
    products: ``ab`` (B + A, weight 1) and ``ac`` (A + C, weight 9; a role
    swap: its two ordered site pairs over divisor 2).  ``bonded``: A-B already
    bonded (one molecule XY with C left) and the homolysis ``open`` of that
    bond loaded, whose products are X and Y again."""
    from rmgpy.kmc import melt_branched as mb
    from rmgpy.kmc.melt_state import SetBond, SetRadical
    spec = law_spec()
    state = mb.build_molecules(spec, [(["Ph", "M", "Ph"], _chain_links(2)),
                                      (["Ph", "M", "Ph", "M", "Ph"], _chain_links(3))])
    atom = mb.residue_atoms(state)
    a, c, b = atom[(0, 0)]["r4"], atom[(1, 0)]["r4"], atom[(1, 3)]["m"]
    for x in (a, c, b):
        S.strip_h(state, x)
    templates = [_recombination(state, b, a, "ab", 1.0), _recombination(state, a, c, "ac", 9.0)]
    if bonded:
        state.apply([SetBond(a, b, 0.0, 1.0), SetRadical(a, 1, 0), SetRadical(b, 1, 0)])
        templates.append(_homolysis(state, b, a, "open", 1.0))
    return me.MeltEngine(spec, S.ToyCatalog(builder=lambda _t: list(templates)), 700.0,
                         np.random.default_rng(0), state=state, excluded=excluded)


def _unit_k_d(monkeypatch):
    monkeypatch.setattr(me.TransportLaw, "k_diffusion", lambda self, t, n1, n2: 1.0)
    monkeypatch.setattr(me.TransportLaw, "k_diffusion_still", lambda self, t, n1, n2: 1.0)


@pytest.mark.parametrize("excluded, phi", [((), 1.0 / 11.0), (("ac",), 1.0 / 2.0)])
def test_competing_channels_share_one_escape_factor_forward_and_in_the_reverse_thinning(monkeypatch, excluded, phi):
    """The encounter group (``encounter_group``) is the component pair's, not
    the inverse's: with k_D = 1 and activations 1 (``ab``) and 9 (``ac``, a
    channel to other products), the pair's forward phi is 1/11, and the
    homolysis of the A-B bond (the inverse of ``ab`` alone) is thinned at
    1/11 too (an inverse-only group would give 1/2), in the fired proposal and
    in the open-segment mean (``_fallback_phi``).  An excluded class leaves
    the group in both directions (1/2).  Today's engine already does this for
    two-molecule pairs (I-301: green before the refactor, kept green)."""
    _unit_k_d(monkeypatch)
    forward = _competing(False, excluded)
    (rate,) = forward.pairs.values()
    names = [t.name for t in forward.templates]
    activation = {names[row[0]]: forward._row_weight(row) * (row[1] + row[2]) for row in rate.classes}
    assert activation == {k: v for k, v in {"ab": 1.0, "ac": 9.0}.items() if k not in excluded}
    assert rate.k_d == 1.0
    assert rate.phi == pytest.approx(phi, rel=1.0e-15)
    reverse = _competing(True, excluded)
    ci = [t.name for t in reverse.templates].index("open")
    assert sum(len(e.shape.uni.get(ci, ((), 0))[0]) for e in reverse.entries.values()) == 1
    mean, sampled = reverse._fallback_phi(ci)
    assert not sampled and mean == pytest.approx(phi, rel=1.0e-15)
    reverse.run(max_events=1)
    assert reverse.template_counts["open"] == 1 and reverse.met_counts["reverse_accepted"] == 1
    assert reverse.phi_sum / reverse.thinned_proposals == pytest.approx(phi, rel=1.0e-15)


def test_both_thinning_paths_take_the_escape_factor_from_one_evaluator(monkeypatch):
    """``_fire``'s reverse thinning and ``_fallback_phi`` read phi from
    ``_reverse_phi`` and nowhere else: replaced by a constant 0.25, both
    report 0.25 (a path computing its own phi would report 1/11 here)."""
    _unit_k_d(monkeypatch)
    engine = _competing(True)
    calls = []

    def constant(self, tr, products):
        calls.append(len(products))
        return 0.25, False
    monkeypatch.setattr(me.MeltEngine, "_reverse_phi", constant)
    ci = [t.name for t in engine.templates].index("open")
    assert engine._fallback_phi(ci) == (0.25, False)
    engine.run(max_events=1)
    assert engine.template_counts["open"] == 1
    assert engine.phi_sum == pytest.approx(0.25 * engine.thinned_proposals, rel=1.0e-15)
    assert set(calls) == {2} and len(calls) >= 1 + engine.thinned_proposals


def _two_aryl_chain():
    """Ph0-M1-Ph2-M3-Ph4 with an aryl radical at both end rings' para sites
    (one window up to a mirror: an aryl-aryl recombination is a role swap)."""
    from rmgpy.kmc import melt_branched as mb
    spec = law_spec(d_ref=_NO_SPLIT_D_REF)
    state = mb.build_molecules(spec, [(["Ph", "M", "Ph", "M", "Ph"], _chain_links(3))])
    atom = mb.residue_atoms(state)
    a, c = atom[(0, 0)]["r4"], atom[(0, 4)]["r4"]
    for x in (a, c):
        S.strip_h(state, x)
    return spec, state, a, c


def test_the_self_pair_rows_count_each_ordered_site_pair_once():
    """phi_self's mapping (``self_pair_rows``): each ordered site pair of one
    molecule, sharing no atom, is one pair of its class's row over the
    class's divisor; ``pair_rate(entry, entry)`` counts every ordered pair
    twice (both role assignments are the one molecule) and keeps the pairs
    that share an atom.  Opened ring (methylene + aryl radical, divisor 1):
    k_act 0.01 against pair_rate's 0.02 (round 58 item 3).  Two aryl radicals
    on one chain (role swap, divisor 2): ordered pairs (A, C) and (C, A) over
    2 is one pair, k_act w; pair_rate(e, e) gives 4 w ((A, C), (C, A), (A, A),
    (C, C), twice, over 2).  The rows also price the measured same_molecule
    leak (same count, same weights).  phi_self = k_D(L, L) / (k_D + k_act)
    with k_D written out (``_analytic_k_d``)."""
    spec, cut, m, r = _opened_ring(_NO_SPLIT_D_REF)
    engine = _engine(spec, cut, [_recombination(cut, m, r, "close", 1.0e-2)])
    (entry,) = engine.entries.values()
    rows = engine.self_pair_rows(entry)
    assert rows.classes == ((0, 1.0, 0.0),)
    assert rows.k_act == rows.k_act_allowed == pytest.approx(1.0e-2, rel=1.0e-15)
    assert engine.pair_rate(entry, entry).k_act == pytest.approx(2.0e-2, rel=1.0e-15)
    k_d = _analytic_k_d("self", engine.transport, entry.met_length)
    assert rows.k_d_self == pytest.approx(k_d, rel=1.0e-12)
    assert rows.phi_self == pytest.approx(k_d / (k_d + 1.0e-2), rel=1.0e-12)
    assert _classes(engine)["same_molecule"] == pytest.approx(rows.k_act / engine.nav, rel=1.0e-12)

    spec, state, a, c = _two_aryl_chain()
    engine = _engine(spec, state, [_recombination(state, a, c, "aa", 3.0e-2)])
    assert engine.classes[0].divisor == 2, "setup: the aryl-aryl recombination is a role swap"
    (entry,) = [e for e in engine.entries.values() if e.species is None]
    rows = engine.self_pair_rows(entry)
    assert sum(r[1] + r[2] for r in rows.classes) == 1.0
    assert rows.k_act == pytest.approx(3.0e-2, rel=1.0e-15)
    assert engine.pair_rate(entry, entry).k_act == pytest.approx(4 * 3.0e-2, rel=1.0e-15)
    assert _classes(engine)["same_molecule"] == pytest.approx(rows.k_act / engine.nav, rel=1.0e-12)


def test_the_self_pair_rows_leave_excluded_classes_out():
    """An excluded class is no row of the same-molecule group either (as in
    ``encounter_group``); with only it loaded there are no rows."""
    spec, cut, m, r = _opened_ring(_NO_SPLIT_D_REF)
    engine = me.MeltEngine(spec, S.ToyCatalog(builder=lambda _t: [_recombination(cut, m, r, "close", 1.0e-2)]),
                           700.0, np.random.default_rng(0), state=cut, excluded=["close"])
    (entry,) = engine.entries.values()
    assert engine.self_pair_rows(entry) is None
