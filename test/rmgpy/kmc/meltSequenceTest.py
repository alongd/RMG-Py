"""The melt kMC of a linear copolymer declaration (schema v2): the initial
state sampled from the declared sequence statistics (``melt_sequence``,
``melt_state.materialize``), the realized composition and its declared
tolerance, the mass-based MET length, and the engine, observer and runner
on a v2 deck.

Expected values come from the declaration alone: the sampler is re-derived
here from its definition (``markov_chain_v1``), masses from the residue
graphs' atoms, and the toy PE catalog (``melt_engine_support``) stands in
for an emitted E/P catalog."""

import copy
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

import melt_engine_support as S
from meltRunTest import DECK_HEAD, DECK_TAIL, toy_payloads, write_catalog

from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import melt_observables as mo
from rmgpy.kmc import melt_sequence as msq
from rmgpy.kmc import melt_state as ms
from rmgpy.kmc import run as kr
from rmgpy.kmc.melt_state import _element_mass
from rmgpy.polymer_spec import PolymerSpec

from polymer_spec_v2_support import DECKS, polymer_keyword, q

EP_DECK = DECKS / "ep_v2_deck.py"
PE_LIVE = ("homolysis", "recombination", "h_abstraction", "disproportionation", "beta_scission",
           "mid_end_recombination", "intra_mid_end")
TOLERANCE = q(1, 10)
# a Markov chain whose stationary distribution is the fixture's (1/2, 1/2)
STICKY = {"E": {"E": q(3, 4), "P": q(1, 4)}, "P": {"E": q(1, 4), "P": q(3, 4)}}


def ep_mapping(*, tolerance=TOLERANCE, met=True, statistics=None, counts=None, seed=None):
    """The E/P fixture's declaration with the runtime blocks a melt run reads
    (labelled test values: tolerance, MET parameters)."""
    mapping = copy.deepcopy(polymer_keyword(EP_DECK, "polymer_spec"))
    ensemble = mapping["initial_ensemble"]
    if tolerance is not None:
        ensemble["composition_relative_tolerance"] = tolerance
    if statistics is not None:
        ensemble["sequence_statistics"] = statistics
    if counts is not None:
        ensemble["dp_counts"] = dict(counts)
    if seed is not None:
        ensemble["seed"] = seed
    if met:
        mapping["transport_model"] = {"met": S.pe_met()}
    return mapping


def ep_spec(**kw):
    return PolymerSpec.from_mapping(ep_mapping(**kw))


# --------------------------------------------------------------------------
# the sampler, re-derived from its definition
# --------------------------------------------------------------------------
def reference_sequences(spec):
    """``markov_chain_v1``: the uniform u of draw (chain, step) is the first 8
    bytes (big-endian) of sha256("markov_chain_v1:<seed>:<chain>:<step>")
    over 2^64; the state drawn is the first, in state order, whose cumulative
    weight exceeds u.  Step 0 draws from pi, step k from the row of step k-1."""
    ensemble = spec.initial_ensemble
    composition = dict(ensemble.composition)
    present = sorted(r for r, x in composition.items() if x > 0)
    statistics = ensemble.statistics
    if statistics.model == "alternating":
        cycle = statistics.cycle
        states = list(range(len(cycle)))
        pi = {k: Fraction(1, len(cycle)) for k in states}
        rows = {k: {(k + 1) % len(cycle): Fraction(1)} for k in states}
        name = dict(enumerate(cycle))
    else:
        states = present
        pi = {r: composition[r] for r in present}
        rows = ({i: dict(pi) for i in present} if statistics.model == "ideal_random"
                else {i: dict(row) for i, row in statistics.transition})
        name = {r: r for r in present}

    def draw(chain, step, weights):
        digest = hashlib.sha256(f"markov_chain_v1:{ensemble.seed}:{chain}:{step}".encode()).digest()
        u = Fraction(int.from_bytes(digest[:8], "big"), 2 ** 64)
        total = Fraction(0)
        for state in states:
            total += weights.get(state, Fraction(0))
            if u < total:
                return state
        raise AssertionError("weights do not sum to one")

    out = []
    for dp, count in sorted(ensemble.chains.counts.items()):
        for _ in range(count):
            chain = len(out)
            state = draw(chain, 0, pi)
            sequence = [state]
            for step in range(1, dp):
                state = draw(chain, step, rows[state])
                sequence.append(state)
            out.append(tuple(name[s] for s in sequence))
    return tuple(out)


def first_difference(a, b):
    """The first chain where two samples differ (None: equal); keeps a failing
    comparison of whole melts cheap to report."""
    if len(a) != len(b):
        return f"{len(a)} chains != {len(b)} chains"
    return next(((i, x, y) for i, (x, y) in enumerate(zip(a, b)) if x != y), None)


def residue_mass(residue) -> Fraction:
    """g/mol, from the residue graph's atoms (every hydrogen counted)."""
    total = Fraction(0)
    for atom in residue.graph.atoms:
        total += _element_mass(atom.element, None if atom.isotope == -1 else atom.isotope)
        total += atom.hydrogens * _element_mass("H", None)
    return total * 1000


def mean_residue_mass(spec) -> Fraction:
    residues = dict(spec.residues)
    return sum((x * residue_mass(residues[r]) for r, x in spec.initial_ensemble.composition), Fraction(0))


def chain_sequences(state):
    """Per chain ordinal, the residue names of its units in unit order (from
    the atoms' origin)."""
    units = {}
    for origin in state.origin.values():
        if origin.kind != "chain":
            continue
        residue = ms.segment_residue(origin.segment)
        if residue is not None:
            unit = int(origin.segment.split()[1])
            units.setdefault(origin.ordinal, {})[unit] = residue
    return tuple(tuple(units[c][k] for k in sorted(units[c])) for c in sorted(units))


# --------------------------------------------------------------------------
# the tolerance field
# --------------------------------------------------------------------------
def test_the_composition_tolerance_is_an_optional_exact_rational_outside_the_chemistry_digest():
    with_tolerance, without = ep_spec(), ep_spec(tolerance=None)
    assert with_tolerance.initial_ensemble.composition_tolerance == Fraction(1, 10)
    assert without.initial_ensemble.composition_tolerance is None
    data = json.loads(with_tolerance.normalized_json)["initial_ensemble"]
    assert data["composition_relative_tolerance"] == {"numerator": 1, "denominator": 10}
    assert "composition_relative_tolerance" not in json.loads(without.normalized_json)["initial_ensemble"]
    # a deck that declares none keeps its bytes
    assert ep_spec(tolerance=None, met=False).normalized_json == kr.read_deck(EP_DECK).spec.normalized_json
    assert kr.chemistry_sha256(with_tolerance) == kr.chemistry_sha256(without)
    assert kr.declaration_sha256(with_tolerance) != kr.declaration_sha256(without)
    assert PolymerSpec.from_mapping(json.loads(with_tolerance.normalized_json)).normalized_json \
        == with_tolerance.normalized_json
    for bad in (0.1, q(-1, 10), "1/10"):
        with pytest.raises(ValueError, match="composition_relative_tolerance"):
            ep_spec(tolerance=bad)


# --------------------------------------------------------------------------
# sampling
# --------------------------------------------------------------------------
@pytest.mark.parametrize("statistics", [None, {"model": "markov1", "transition": STICKY},
                                        {"model": "alternating", "cycle": ["E", "P"]}],
                         ids=["ideal_random", "markov1", "alternating"])
def test_the_sampler_is_markov_chain_v1_as_defined(statistics):
    spec = ep_spec(statistics=statistics)
    sequences = msq.sample_sequences(spec)
    assert first_difference(sequences, reference_sequences(spec)) is None
    assert [len(s) for s in sequences] == [30] * 100


def test_markov1_transitions_follow_the_rows_not_the_stationary_distribution():
    """E->E is 3/4 under the declared rows; i.i.d. draws from pi give 1/2."""
    sequences = msq.sample_sequences(ep_spec(statistics={"model": "markov1", "transition": STICKY}))
    pairs = [(a, b) for s in sequences for a, b in zip(s, s[1:])]
    stay = sum(1 for a, b in pairs if a == b) / len(pairs)
    assert abs(stay - 0.75) < 0.03


def test_alternating_chains_alternate():
    sequences = msq.sample_sequences(ep_spec(statistics={"model": "alternating", "cycle": ["E", "P"]}))
    assert {s[0] for s in sequences} == {"E", "P"}
    assert all(a != b for s in sequences for a, b in zip(s, s[1:]))


def test_the_declared_seed_decides_the_sample():
    assert first_difference(msq.sample_sequences(ep_spec()), msq.sample_sequences(ep_spec())) is None
    assert first_difference(msq.sample_sequences(ep_spec()), msq.sample_sequences(ep_spec(seed=2))) is not None


# --------------------------------------------------------------------------
# the materialized melt
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def ep_state():
    return ms.materialize(ep_spec())


def test_materialize_builds_the_sampled_chains_head_to_tail(ep_state):
    spec = ep_spec()
    assert first_difference(chain_sequences(ep_state), msq.sample_sequences(spec)) is None
    assert len(ep_state.components) == 100
    assert sorted(ep_state.atoms) == list(range(len(ep_state.atoms)))
    assert ms.materialize(ep_spec()).full_hash() == ep_state.full_hash()
    residues = dict(spec.residues)
    junctions = spec.junctions
    sequence = msq.sample_sequences(spec)[0]
    atom = {(o.segment, o.key): i for i, o in ep_state.origin.items() if o.ordinal == 0 and o.kind == "chain"}
    for k, (a, b) in enumerate(zip(sequence, sequence[1:]), start=1):
        exit_atom, entry_atom, order = junctions[(a, b)]
        x, y = atom[(f"unit {k} {a}", exit_atom)], atom[(f"unit {k + 1} {b}", entry_atom)]
        assert ep_state.atoms[x].edges[y] == order
    # the backbone is each residue's entry-site -> exit-site path
    for (segment, key), i in atom.items():
        name = ms.segment_residue(segment)
        if name is not None:
            assert ep_state.origin[i].is_backbone == (key in msq.residue_backbone_keys(spec, name))
    assert msq.residue_backbone_keys(spec, "P") == frozenset({"p1", "p2"})
    # every chain's composition is its residues plus the two caps
    for cid in ep_state.components.ids():
        ordinal = ep_state.origin[cid].ordinal
        mass = sum((residue_mass(residues[r]) for r in msq.sample_sequences(spec)[ordinal]), Fraction(0))
        mass += 2 * _element_mass("H", None) * 1000
        assert ep_state.components.molar_mass_g_mol(cid) == mass


def test_segment_names_the_residue():
    assert ms.segment_residue("unit 3 E") == "E"
    assert ms.segment_residue("unit 3") is None
    assert ms.segment_residue("head_end") is None


def test_the_realized_composition_is_refused_above_the_declared_tolerance():
    spec = ep_spec()
    sequences = msq.sample_sequences(spec)
    residues = [r for s in sequences for r in s]
    realized = {r: Fraction(residues.count(r), len(residues)) for r in ("E", "P")}
    assert msq.realized_composition(spec) == realized
    deviation = max(abs(realized[r] - Fraction(1, 2)) / Fraction(1, 2) for r in realized)
    assert deviation > 0
    tolerance = {"numerator": deviation.numerator, "denominator": deviation.denominator}
    ms.materialize(ep_spec(tolerance=tolerance))         # at the tolerance: accepted
    half = deviation / 2
    with pytest.raises(ms.MeltStateError, match="^refused:composition_tolerance: .*E.*P"):
        ms.materialize(ep_spec(tolerance={"numerator": half.numerator, "denominator": half.denominator}))
    with pytest.raises(ms.MeltStateError, match="^refused:composition_tolerance: .*declares no"):
        ms.materialize(ep_spec(tolerance=None))


ETHYLENE = {"atoms": [S._atom("a", hydrogens=2), S._atom("b", hydrogens=2)],
            "bonds": [{"atoms": ["a", "b"], "order": 2}], "electronic_state": {"multiplicity": 1}}


def test_an_initially_mobile_monomer_product_gets_backbone_origin():
    mapping = ep_mapping(counts={5: 4}, tolerance=q(1, 1))
    mapping["monomer_products"] = {"E": ETHYLENE}
    mapping["small_species"] = {"ethylene": {"graph": ETHYLENE, "disposition": "mobile"}}
    mapping["initial_mobile_counts"] = {"ethylene": 2}
    mapping["density_basis"] = "whole_initial_mixture"
    state = ms.materialize(PolymerSpec.from_mapping(mapping))
    mobile = [i for i, o in state.origin.items() if o.kind == "mobile"]
    assert len([i for i in mobile if state.atoms[i].element == "C"]) == 4
    assert all(state.origin[i].is_backbone == (state.atoms[i].element == "C") for i in mobile)
    assert len(state.components) == 6


# --------------------------------------------------------------------------
# engine, transport and observer
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def ep_engine():
    return me.MeltEngine(ep_spec(), S.ToyCatalog(names=PE_LIVE), 700.0, np.random.default_rng(3))


def test_the_met_length_is_the_chain_mass_over_the_mean_residue_mass(ep_engine):
    spec = ep_spec()
    m0 = mean_residue_mass(spec)
    assert m0 == msq.mean_residue_molar_mass(spec)
    assert ep_engine.transport.repeat_unit_molar_mass_g_mol == float(m0)
    lengths = set()
    for cid, entry in ep_engine.entries.items():
        mass = ep_engine.state.components.molar_mass_g_mol(cid)
        assert entry.met_length == max(1, math.floor(mass / m0))
        lengths.add(entry.met_length)
    assert lengths != {30}          # not the residue count: the composition varies per chain


def test_the_volume_is_the_sampled_melt_mass_over_the_density(ep_engine):
    spec = ep_spec()
    total = sum((ep_engine.state.components.mass_kg_per_mol(c) for c in ep_engine.state.components.ids()),
                Fraction(0))
    from rmgpy import constants
    expected = float(total) / constants.Na / spec.density_kg_m3
    assert ep_engine.volume == pytest.approx(expected, rel=1e-12)
    init = msq.initialization(spec)
    assert init.scale_factor == pytest.approx(spec.initial_mass_kg * constants.Na / float(total), rel=1e-12)


def test_the_observer_reports_composition_beside_the_declared_one(ep_engine):
    spec = ep_spec()
    observer = mo.MeltObserver(spec, xs=4, identity=ep_engine.identity)
    obs = observer.observe(ep_engine)
    realized = msq.realized_composition(spec)
    block = obs["composition"]
    assert block["declared"] == {"E": 0.5, "P": 0.5}
    assert block["realized_initial"] == {r: float(x) for r, x in realized.items()}
    assert block["polymer_derived"] == block["realized_initial"]
    assert block["relative_tolerance"] == 0.1
    assert obs["length_rule"] == "mean_residue_molar_mass"
    assert obs["mom_mass_law"] is None and obs["tier3"] is None
    masses = [float(ep_engine.state.components.molar_mass_g_mol(c)) for c in ep_engine.state.components.ids()]
    assert obs["tier2"]["Mn_g_mol"] == pytest.approx(sum(masses) / len(masses), rel=1e-12)
    with pytest.raises(mo.ObservablesError, match="refused:schema_v2_mom"):
        mo.mom_reported_quantities(obs)


def test_side_group_channels_are_refused_for_a_v2_declaration():
    with pytest.raises(mo.ObservablesError, match="refused:schema_v2_mom"):
        mo.MeltObserver(ep_spec(), xs=4, side_group_homolysis=[{"label": "x"}])


def test_the_engine_fires_on_the_copolymer_and_conserves_atoms():
    engine = me.MeltEngine(ep_spec(counts={20: 6}, tolerance=q(1, 2)), S.ToyCatalog(names=PE_LIVE), 700.0,
                           np.random.default_rng(1))
    before = mo.exact_totals(engine.state)
    checker = mo.InvariantChecker(engine.state)
    result = engine.run(max_events=60, check=checker)
    assert result.events == 60
    assert mo.exact_totals(engine.state).composition == before.composition


# --------------------------------------------------------------------------
# the runner
# --------------------------------------------------------------------------
def write_v2_deck(path, mapping):
    Path(path).write_text(DECK_HEAD + f"polymer({mapping['label']!r}, polymer_spec={mapping!r}, cutoff=4)\n"
                          + DECK_TAIL)
    return Path(path)


def ep_catalog(path, spec):
    write_catalog(path, toy_payloads([t for t in S.pe_templates() if t.name in PE_LIVE]), tag="ep-toy")
    doc = json.loads(Path(path).read_text())
    doc["provenance"]["chemistry_sha256"] = kr.chemistry_sha256(spec)
    Path(path).write_text(json.dumps(doc, sort_keys=True, indent=1) + "\n")
    return Path(path)


def test_the_runner_runs_a_v2_deck(tmp_path):
    mapping = ep_mapping(counts={20: 6}, tolerance=q(1, 2))
    deck = write_v2_deck(tmp_path / "ep.py", mapping)
    spec = PolymerSpec.from_mapping(mapping)
    catalog = ep_catalog(tmp_path / "ep_catalog.json", spec)
    manifest = kr.run(deck, [catalog], 700.0, [0.0, 1.0e-3], [1, 2], tmp_path / "out", max_events=50)
    assert manifest["status"] == "complete"
    realized = msq.realized_composition(spec)
    assert manifest["composition"]["declared"] == {"E": "1/2", "P": "1/2"}
    assert manifest["composition"]["realized"] == {r: str(x) for r, x in realized.items()}
    assert manifest["composition"]["relative_tolerance"] == "1/2"
    assert manifest["comparison"] is None and manifest["comparison_refused"].startswith("refused:schema_v2_mom")
    for run in manifest["runs"]:
        rows = [json.loads(x) for x in (tmp_path / "out" / run["path"]).read_text().splitlines()]
        observations = [r for r in rows if "record" not in r]
        assert observations[0]["composition"]["realized_initial"] == {r: float(x) for r, x in realized.items()}
        first, last = observations[0]["invariants"], observations[-1]["invariants"]
        assert first["composition"] == last["composition"]     # C and H close


def test_the_runner_refuses_a_mom_file_for_a_v2_deck(tmp_path):
    mapping = ep_mapping(counts={20: 2}, tolerance=q(1, 1))
    deck = write_v2_deck(tmp_path / "ep.py", mapping)
    catalog = ep_catalog(tmp_path / "ep_catalog.json", PolymerSpec.from_mapping(mapping))
    mom = tmp_path / "mom.json"
    mom.write_text("{}")
    with pytest.raises(kr.RunError, match="refused:schema_v2_mom"):
        kr.run(deck, [catalog], 700.0, [0.0], [1], tmp_path / "out", mom=mom)
    assert not (tmp_path / "out").exists()
