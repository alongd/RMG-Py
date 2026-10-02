"""Contract tests for the isothermal type-aggregated SSA engine."""

import math
import os
import hashlib
import importlib.util
import json
import multiprocessing
import subprocess
import sys
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import chisquare

from rmgpy.kmc.ssa import (
    E_PROPENSITY_NEGATIVE,
    E_PROPENSITY_NONFINITE,
    ConstantPairKernel,
    IsothermalSSA,
    LengthBinnedThinningSampler,
    METPairKernel,
    PairItem,
    PropensityError,
    SiteIndex,
    SumPairKernel,
    ThinningBoundError,
    _forward_rate_bound,
    build_met_population,
    canonical_propensities,
    epsilon_propensity,
)
from rmgpy.kmc.event_record import EventRecord
from rmgpy.kmc.met import (
    N_A,
    PAIR_CLASSES,
    PS_C_R2,
    PS_M0,
    SIGMA_CONTACT,
    TRANSPORT_ARMS,
    collins_kimball,
    compile_bulk_table,
    diffusion_rate,
)
from rmgpy.kmc.compiler import compiler_source_hash, validate_artifact
from rmgpy.kmc.state import AtomRef, KMCState, Site, Strand


PROPENSITY_INPUT = (
    ("event-00", "1/8"),
    ("event-01", "3/8"),
    ("event-02", "1/2"),
    ("zero-event-00", "0"),
    ("zero-event-01", "0"),
    ("invalid-nan", "NaN"),
    ("invalid-inf", "+Inf"),
    ("invalid-negative", "-1/8"),
)
PROPENSITY_EXPECTED = (
    ("event_ids", ("event-00", "event-01", "event-02")),
    ("exact_total", "1"),
    ("pairwise_binary64", repr(math.fsum((1 / 8, 3 / 8, 1 / 2)))),
    ("zero_fixture", "absorbing:no-draw"),
    (
        "invalid",
        "E_PROPENSITY_NONFINITE,E_PROPENSITY_NONFINITE,E_PROPENSITY_NEGATIVE",
    ),
)
REPO_ROOT = Path(__file__).resolve().parents[3]
REAL_DATABASE_PATH = Path(
    os.environ.get("RMG_DATABASE_PATH", "/home/alon/Code/RMG-database")
)
REAL_CACHE_ROOT = Path(
    os.environ.get(
        "RMG_KMC_CACHE_ROOT", str(Path(__file__).with_name(".real-event-cache"))
    )
)


def slow(test_function):
    test_function = pytest.mark.slow(test_function)
    return pytest.mark.skipif(
        os.environ.get("RMG_KMC_SLOW") != "1",
        reason="set RMG_KMC_SLOW=1",
    )(test_function)


def _git_head(path):
    return subprocess.check_output(
        ["git", "-C", str(path), "rev-parse", "HEAD"], text=True
    ).strip()


@pytest.fixture(scope="session")
def real_ps_artifact():
    cache_key = f"{_git_head(REPO_ROOT)}-{_git_head(REAL_DATABASE_PATH)}-{compiler_source_hash()}"
    cache_dir = REAL_CACHE_ROOT / cache_key
    cache_dir.mkdir(parents=True, exist_ok=True)
    artifacts = sorted(cache_dir.glob("*.json"))
    if not artifacts:
        fixture_script = Path(__file__).with_name("compile_event_set_fixture.py")
        environment = os.environ.copy()
        environment.update(
            {
                "PYTHONPATH": str(REPO_ROOT),
                "PYTHONHASHSEED": "0",
                "MPLCONFIGDIR": str(cache_dir / "matplotlib"),
            }
        )
        with (cache_dir / "compile.stdout.log").open("w") as stdout, (
            cache_dir / "compile.stderr.log"
        ).open("w") as stderr:
            subprocess.run(
                [
                    sys.executable,
                    str(fixture_script),
                    str(REAL_DATABASE_PATH),
                    str(cache_dir),
                ],
                cwd=REPO_ROOT,
                env=environment,
                stdout=stdout,
                stderr=stderr,
                check=True,
                timeout=4 * 60 * 60,
            )
        artifacts = sorted(cache_dir.glob("*.json"))
    assert len(artifacts) == 1
    payload = artifacts[0].read_bytes()
    assert hashlib.sha256(payload).hexdigest() == artifacts[0].stem
    artifact = json.loads(payload)
    validate_artifact(artifact)
    return artifact


@pytest.fixture(scope="session")
def independent_state_oracle():
    path = Path(__file__).with_name("stateTest.py")
    spec = importlib.util.spec_from_file_location("kmc_state_independent_oracle", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_propensity_numerics_fixture_and_absorbing_state_do_not_advance_rng():
    expected = dict(PROPENSITY_EXPECTED)
    vector = canonical_propensities(
        [(event_id, float(Fraction(value))) for event_id, value in PROPENSITY_INPUT[:3]]
    )
    assert vector.event_ids == expected["event_ids"]
    assert repr(vector.total) == expected["pairwise_binary64"]
    assert abs(vector.total - math.fsum(vector.values)) <= epsilon_propensity(3)
    exact_probabilities = tuple(
        float(Fraction(value)) for _, value in PROPENSITY_INPUT[:3]
    )
    assert all(
        abs(value / vector.total - exact_probability)
        <= epsilon_propensity(3) * exact_probability
        for value, exact_probability in zip(vector.values, exact_probabilities)
    )

    zero = canonical_propensities(
        [(event_id, float(value)) for event_id, value in PROPENSITY_INPUT[3:5]]
    )
    rng = np.random.default_rng(921_013)
    before = rng.bit_generator.state
    assert zero.draw(rng) is None
    assert rng.bit_generator.state == before
    assert zero.absorbing
    assert expected["zero_fixture"] == "absorbing:no-draw"


@pytest.mark.parametrize(
    ("value", "code"),
    [
        (math.nan, E_PROPENSITY_NONFINITE),
        (math.inf, E_PROPENSITY_NONFINITE),
        (-1 / 8, E_PROPENSITY_NEGATIVE),
    ],
)
def test_propensity_numerics_reject_invalid_values_before_summation(value, code):
    with pytest.raises(PropensityError) as caught:
        canonical_propensities((("event", value), ("other", 1.0)))
    assert caught.value.code == code


def _one_atom_radical_record(status="enabled"):
    reason = "test frontier" if status != "enabled" else ""
    return EventRecord(
        family="toy",
        arity=1,
        participant_site_types=["radical"],
        reactant_multiplicities=[1],
        rate_order=1,
        rate_units="s^-1",
        k_table={
            "T": [600.0, 800.0],
            "k": [2.0, 2.0],
            "interpolation": "linear-ln-k",
            "extrapolation": "refuse",
        },
        bond_ops=[{"action": "set_radical", "atom": 0, "value": 0}],
        radical_delta=-1,
        status=status,
        status_reason=reason,
        reactant_graphs=["1 *1 C u1 p0 c0"],
        product_graphs=["1 *1 C u0 p0 c0"],
    )


def _leak_records(k1, k2):
    enabled = EventRecord(
        family="toy-enabled",
        arity=1,
        participant_site_types=["A"],
        reactant_multiplicities=[1],
        rate_order=1,
        rate_units="s^-1",
        k_table={
            "T": [600.0, 800.0],
            "k": [k1, k1],
            "interpolation": "linear-ln-k",
            "extrapolation": "refuse",
        },
        bond_ops=[{"action": "set_radical", "atom": 0, "value": 0}],
        radical_delta=-1,
        reactant_graphs=["1 *1 C u1 p0 c0"],
        product_graphs=["1 *1 C u0 p0 c0"],
    )
    refused = EventRecord(
        family="toy-refused",
        arity=1,
        participant_site_types=["B"],
        reactant_multiplicities=[1],
        rate_order=1,
        rate_units="s^-1",
        k_table={
            "T": [600.0, 800.0],
            "k": [k2, k2],
            "interpolation": "linear-ln-k",
            "extrapolation": "refuse",
        },
        status="refused",
        status_reason="toy leak frontier",
        reactant_graphs=["1 *1 C u0 p0 c0"],
        product_graphs=["1 *1 C u0 p0 c0"],
    )
    return enabled, refused


def _run_leak(k1, k2, horizon, seed):
    state = KMCState((_radical_strand("strand", "strand-C", 5),))
    engine = IsothermalSSA(
        state,
        _leak_records(k1, k2),
        temperature=700.0,
        volume=1.0,
        rng=np.random.default_rng(seed),
    )
    result = engine.run(until=horizon)
    fired_at = result.events[0].time if result.events else None
    expected_refused = k2 * (horizon - fired_at) if fired_at is not None else 0.0
    expected_enabled = k1 * (fired_at if fired_at is not None else horizon)
    assert result.leak_num == pytest.approx(expected_refused, rel=1e-12, abs=1e-14)
    assert result.leak_den == pytest.approx(
        expected_refused + expected_enabled, rel=1e-12, abs=1e-14
    )
    return result


def test_leak_toy_integrates_refused_hazard_over_realized_event_intervals():
    _run_leak(1.0, 0.2, 2.0, 6_013_700)


def _radical_strand(strand_id, atom_uuid, length):
    ref = AtomRef(atom_uuid, strand_id, 0, "chain_end")
    return Strand(
        strand_id,
        length,
        formula={"C": 1},
        radical_count=1,
        atom_refs={0: ref},
        atom_graph={
            atom_uuid: {
                "element": "C",
                "radical": 1,
                "charge": 0,
                "lone_pairs": 0,
                "implicit_hydrogens": 0,
                "edges": {},
            }
        },
    )


def test_site_index_updates_only_touched_graph_components_after_apply():
    record = _one_atom_radical_record()
    state = KMCState(
        (
            _radical_strand("left", "left-C", 4),
            _radical_strand("right", "right-C", 9),
        )
    )
    index = SiteIndex(state, (record,))
    assert len(index.candidates(record.event_id)) == 2
    untouched = deepcopy(state.strands["right"].atom_graph)

    selected = next(
        candidate
        for candidate in index.candidates(record.event_id)
        if candidate[0].strand_id == "left"
    )
    index.apply(record, selected)

    assert state.strands["right"].atom_graph == untouched
    assert index.last_rescanned_atom_uuids == frozenset({"left-C"})
    assert [
        candidate[0].strand_id for candidate in index.candidates(record.event_id)
    ] == ["right"]
    index.assert_index_consistency()


def _indexed_candidate_keys(index, records):
    return {
        record.event_id: {
            tuple(site.graph[site.atom]["atom_ref"].uuid for site in candidate)
            for candidate in index.candidates(record.event_id)
        }
        for record in records
    }


@slow
def test_real_site_index_matches_independent_brute_force_rescan(
    real_ps_artifact, independent_state_oracle
):
    seed = 7_013_700
    state = independent_state_oracle._seed_evolving_state(
        real_ps_artifact, "R_Addition_MultipleBond", seed
    )
    records = tuple(EventRecord.from_dict(data) for data in real_ps_artifact["records"])
    index = SiteIndex(state, records)
    representatives = {}
    for data in real_ps_artifact["records"]:
        key = (
            tuple(data["reactant_graphs"]),
            data["orientation"],
            data["inventory_class"],
            bool(data["reverse_of"]),
        )
        representatives.setdefault(key, data)
    representative_records = tuple(
        EventRecord.from_dict(data) for data in representatives.values()
    )
    prepared = tuple(
        independent_state_oracle._prepare_stress_record(data)
        for data in representatives.values()
    )
    brute = independent_state_oracle._current_candidates(state, prepared, {})
    brute_keys = {record.event_id: set() for record in representative_records}
    for record, candidate in brute:
        brute_keys[record.event_id].add(
            tuple(site.graph[site.atom]["atom_ref"].uuid for site in candidate)
        )
    indexed = _indexed_candidate_keys(index, representative_records)
    assert indexed == brute_keys
    for record in records:
        key = (
            tuple(record.reactant_graphs),
            record.orientation,
            record.inventory_class,
            bool(record.reverse_of),
        )
        representative = representatives[key]["event_id"]
        if record.orientation != "reversed":
            assert (
                indexed[representative]
                == _indexed_candidate_keys(index, (record,))[record.event_id]
            )


@slow
def test_real_k_act_is_pair_specific_at_700_k(real_ps_artifact):
    from rmgpy.data.rmg import RMGDatabase
    from met_rate_oracle import independent_termination_rates

    records = real_ps_artifact["records"]
    by_event = {record["event_id"]: record for record in records}
    database = RMGDatabase()
    database.load_kinetics(
        str(REAL_DATABASE_PATH / "input/kinetics"),
        reaction_libraries=[],
        seed_mechanisms=None,
        kinetics_families=["R_Recombination", "Disproportionation"],
        kinetics_depositories=["training"],
    )
    database.load_thermo(
        str(REAL_DATABASE_PATH / "input/thermo"),
        thermo_libraries=["primaryThermoLibrary"],
        depository=True,
    )
    independent_rates = independent_termination_rates(real_ps_artifact, database, 700.0)
    expected = {}
    for inventory in ("R0", "R1"):
        groups = Counter()
        for event_id, rate in independent_rates.items():
            record = by_event[event_id]
            if inventory == "R0" and record["inventory_class"] == "R1:J_ring":
                continue
            site_types = tuple(
                sorted(
                    site_type
                    for site_type, count in zip(
                        record["participant_site_types"],
                        record["reactant_multiplicities"],
                    )
                    for _ in range(count)
                )
            )
            groups[site_types] += rate
        expected[inventory] = dict(groups)
    pack_temperature = 700.0
    pack_exponent = -1.00291
    j_para_forward = 1.0 * 1.76793e10 * pack_temperature**pack_exponent
    j_ortho_forward = 3.53586e10 * pack_temperature**pack_exponent
    ring_total = sum(
        rate
        for event_id, rate in independent_rates.items()
        if by_event[event_id]["inventory_class"] == "R1:J_ring"
    )
    assert ring_total == pytest.approx(j_ortho_forward + j_para_forward, rel=1e-9)
    rows = []
    for inventory in ("R0", "R1"):
        table = compile_bulk_table(
            records, TRANSPORT_ARMS["A0_REF_H_CROSS_NEc"], inventory
        )
        for channel in table.channels:
            bound = _forward_rate_bound(channel, 700.0, 725.0)
            dense_maximum = max(
                channel.forward_rate(temperature)
                for temperature in np.linspace(700.0, 725.0, 1001)
            )
            assert bound >= dense_maximum * (1.0 - 1e-14)
        per_pair = Counter()
        for channel in table.channels:
            record = by_event[channel.event_id]
            site_types = tuple(
                sorted(
                    site_type
                    for site_type, count in zip(
                        record["participant_site_types"],
                        record["reactant_multiplicities"],
                    )
                    for _ in range(count)
                )
            )
            per_pair[site_types] += channel.forward_rate(700.0)
        lumped = sum(channel.forward_rate(700.0) for channel in table.channels)
        assert dict(per_pair) == pytest.approx(expected[inventory], rel=1e-12)
        assert lumped == pytest.approx(sum(expected[inventory].values()), rel=1e-12)
        rows.append((inventory, lumped, dict(per_pair)))
    print(f"k_act 700 K: {rows}")


def _real_cycle_state(artifact, oracle, seed):
    state = oracle._seed_evolving_state(artifact, "R_Addition_MultipleBond", seed)
    owners = state._uuid_owners()
    active_strands = {
        owners[atom_uuid] for edge in state.junctions for atom_uuid in edge
    }
    for strand_id, strand in state.strands.items():
        if strand_id not in active_strands:
            strand.atom_refs = {}
            strand.atom_graph = {}
            strand.radical_count = 0
    state._refresh_components()
    return state


def _active_r1_record_data(artifact):
    return tuple(
        data for data in artifact["records"] if data["inventory_class"] == "R1:J_ring"
    )


def _candidate_anchor_keys(candidates):
    return {record.event_id: set() for record, _ in candidates} | {
        record.event_id: {
            tuple(site.graph[site.atom]["atom_ref"].uuid for site in sites)
            for candidate_record, sites in candidates
            if candidate_record.event_id == record.event_id
        }
        for record, _ in candidates
    }


def _brute_r1_met_total(state, oracle, prepared, table, volume):
    brute = oracle._current_candidates(state, prepared, {})
    by_event = {}
    for record, sites in brute:
        if record.arity != 2 or record.radical_delta >= 0:
            continue
        components = tuple(sorted(state.components[site.strand_id] for site in sites))
        if len(set(components)) != 2:
            continue
        by_event.setdefault(record.event_id, []).append((components, sites))
    ortho_ids = {
        prepared_record["record"].event_id
        for prepared_record in prepared
        if prepared_record["record"].junction_ops
        and str(
            prepared_record["record"].junction_ops[0].get("junction_kind", "")
        ).startswith("J_ortho")
        and prepared_record["record"].arity == 2
    }
    pairs = {}
    for channel in table.channels:
        event_ids = ortho_ids if channel.channel_id == "J_ortho" else {channel.event_id}
        matches = [
            match for event_id in event_ids for match in by_event.get(event_id, ())
        ]
        for components, sites in matches:
            pairs.setdefault(components, {"sites": sites, "channels": {}})["channels"][
                channel.channel_id
            ] = channel
    total = 0.0
    for components, data in pairs.items():
        activation = sum(
            channel.forward_rate(700.0) for channel in data["channels"].values()
        )
        site_classes = []
        for site in data["sites"]:
            radical_refs = [
                node["atom_ref"]
                for node in site.graph.values()
                if int(node.get("radical", 0)) > 0
            ]
            site_classes.append(
                "end"
                if any(
                    ref.position in {0, state.strands[ref.strand_id].length - 1}
                    for ref in radical_refs
                )
                else "mid"
            )
        pair_class = "/".join(sorted(site_classes))
        left_length = state.component_lengths[components[0]]
        right_length = state.component_lengths[components[1]]
        total += collins_kimball(
            activation,
            diffusion_rate(
                table.arm,
                700.0,
                left_length,
                right_length,
                pair_class,
                spin_factor=1.0,
            ),
        ) / (N_A * volume)
    return total


def test_met_ssa_exact_real_r1_evolving_trajectory(
    real_ps_artifact, independent_state_oracle
):
    """Compare independent MET totals on the archived ring capture/inverse cycle."""
    seed = 8_013_700
    event_count = 20_000 if os.environ.get("RMG_KMC_SLOW") == "1" else 2_000
    relative_tolerance = 1e-12
    volume = 1e-50
    state = _real_cycle_state(real_ps_artifact, independent_state_oracle, seed)
    initial_ledger = state.total_ledger
    active_data = _active_r1_record_data(real_ps_artifact)
    active_records = tuple(EventRecord.from_dict(data) for data in active_data)
    records = active_records
    table = compile_bulk_table(records, TRANSPORT_ARMS["A0_REF_H_CROSS_NEc"], "R1")
    prepared = tuple(
        independent_state_oracle._prepare_stress_record(data) for data in active_data
    )
    engine = IsothermalSSA(
        state,
        records,
        temperature=700.0,
        volume=volume,
        rng=np.random.default_rng(seed),
        termination_table=table,
        spin_factor=1.0,
    )
    bin_migrations = 0
    previous_bins = sorted(
        length.bit_length() - 1 for length in state.component_lengths.values()
    )
    oracle_cache = {}
    trajectory = []
    for step in range(event_count):
        population = engine._met_population()
        assert population is not None
        if step == 0:
            report = engine.channel_propensities()
            brute_total = math.fsum(report.propensities.values())
            assert abs(report.total_enabled - brute_total) <= (
                epsilon_propensity(len(report.propensities)) * brute_total
            )
        state_hash = state.state_hash()
        if state_hash not in oracle_cache:
            brute_total = _brute_r1_met_total(
                state, independent_state_oracle, prepared, table, volume
            )
            brute_candidates = independent_state_oracle._current_candidates(
                state, prepared, {}
            )
            oracle_cache[state_hash] = (
                brute_total,
                _candidate_anchor_keys(brute_candidates),
            )
        brute_total, brute_keys = oracle_cache[state_hash]
        assert population.sampler.exact_total_propensity == pytest.approx(
            brute_total, rel=relative_tolerance, abs=1e-300
        )
        indexed_keys = _indexed_candidate_keys(engine.index, active_records)
        for record in active_records:
            assert indexed_keys[record.event_id] == brute_keys.get(
                record.event_id, set()
            )
        event = engine.step()
        assert event is not None, f"absorbed at step {step}"
        trajectory.append((event.event_id, event.state_hash))
        state.assert_component_consistency()
        assert state.moments == state.brute_force_moments()
        assert state.total_ledger == initial_ledger
        current_bins = sorted(
            length.bit_length() - 1 for length in state.component_lengths.values()
        )
        bin_migrations += current_bins != previous_bins
        previous_bins = current_bins
    result = engine.run(max_events=0)
    assert len(result.events) == event_count
    assert bin_migrations > 0
    assert result.bound_violations == 0
    if os.environ.get("RMG_KMC_SLOW") == "1":
        replay_state = _real_cycle_state(
            real_ps_artifact, independent_state_oracle, seed
        )
        replay = IsothermalSSA(
            replay_state,
            records,
            temperature=700.0,
            volume=volume,
            rng=np.random.default_rng(seed),
            termination_table=table,
            spin_factor=1.0,
        )
        for expected_event_id, expected_hash in trajectory[:10_000]:
            replay_event = replay.step()
            assert replay_event is not None
            assert (replay_event.event_id, replay_event.state_hash) == (
                expected_event_id,
                expected_hash,
            )
        replay_result = replay.run(max_events=0)
        assert replay_state.total_ledger == initial_ledger
        assert replay_state.moments == replay_state.brute_force_moments()
        print(
            f"ARCHIVED-RING-SMOKE R1: events={len(replay_result.events)} seed={seed} "
            f"L={replay_result.leak:.12g} "
            f"irreversible={replay_result.irreversible_fraction:.12g}"
        )
    print(
        f"MET-SSA-EXACT R1: events={event_count} seed={seed} "
        f"bin_migrations={bin_migrations} bounds={result.bound_checks} "
        f"L={result.leak:.12g} irreversible={result.irreversible_fraction:.12g}"
    )


@slow
def test_met_ssa_exact_real_r0_frozen_state(real_ps_artifact, independent_state_oracle):
    state = independent_state_oracle._seed_evolving_state(
        real_ps_artifact, "R_Addition_MultipleBond", 8_013_701
    )
    records = tuple(EventRecord.from_dict(data) for data in real_ps_artifact["records"])
    index = SiteIndex(state, records)
    table = compile_bulk_table(records, TRANSPORT_ARMS["A0_REF_H_CROSS_NEc"], "R0")
    assert table.channels
    assert all(channel.reversible for channel in table.channels)
    population = build_met_population(
        state,
        index,
        records,
        table,
        temperature=700.0,
        volume=1e-21,
        spin_factor=1.0,
    )
    assert population.sampler.bound_violations == 0
    assert population.sampler.exact_total_propensity >= 0.0


@slow
def test_r0_pristine_chain_smoke_fires_10000_events_and_reports_irreversibility(
    real_ps_artifact, independent_state_oracle
):
    records = tuple(
        EventRecord.from_dict(data)
        for data in real_ps_artifact["records"]
        if data["inventory_class"] != "R1:J_ring"
    )
    homolysis = [
        record
        for record in records
        if record.family == "R_Recombination"
        and record.arity == 1
        and record.radical_delta == 2
        and record.site_type == "pristine"
    ]
    assert homolysis
    selected = None
    for candidate in homolysis:
        graph = independent_state_oracle._parse_adjacency(candidate.reactant_graphs[0])
        for operation in candidate.bond_ops:
            if operation["action"] != "break":
                continue
            left, right = operation["atoms"]
            remaining = set(graph)
            groups = []
            while remaining:
                members = set()
                pending = [min(remaining)]
                while pending:
                    atom = pending.pop()
                    if atom in members:
                        continue
                    members.add(atom)
                    pending.extend(
                        neighbor
                        for neighbor in graph[atom].get("edges", {})
                        if {atom, neighbor} != {left, right} and neighbor not in members
                    )
                remaining -= members
                groups.append(members)
            if (
                len(groups) != 2
                or min(
                    sum(graph[atom]["element"] == "C" for atom in members)
                    for members in groups
                )
                < 8
            ):
                continue
            selected = candidate, graph
            break
        if selected is not None:
            break
    assert selected is not None
    graph = selected[1]
    backbone = {
        atom
        for atom, node in graph.items()
        if node["element"] == "C"
        and all(float(order) == 1.0 for order in node.get("edges", {}).values())
    }
    endpoints = [
        atom
        for atom in backbone
        if len(set(graph[atom].get("edges", {})) & backbone) == 1
    ]
    assert len(endpoints) == 2
    ordered = [min(endpoints)]
    while len(ordered) < len(backbone):
        following = set(graph[ordered[-1]].get("edges", {})) & backbone - set(ordered)
        assert len(following) == 1
        ordered.append(following.pop())
    roles = {atom: "backbone" for atom in backbone}
    positions = {atom: position for position, atom in enumerate(ordered)}
    for atom in graph:
        visited = {atom}
        pending = {atom}
        while not pending & backbone:
            pending = {
                neighbor
                for current in pending
                for neighbor in graph[current].get("edges", {})
                if neighbor not in visited
            }
            assert pending
            visited |= pending
        anchors = pending & backbone
        assert len(anchors) == 1
        roles[(atom, "position")] = positions[anchors.pop()]
    state, _, _ = independent_state_oracle._harness(selected[0], roles=roles)
    assert state.total_radicals == 0
    initial_ledger = state.total_ledger
    carbon_by_graph = {}
    carbon_by_event = {}
    for record in records:
        carbon = 0
        for adjacency in record.reactant_graphs:
            if adjacency not in carbon_by_graph:
                carbon_by_graph[adjacency] = sum(
                    node["element"] == "C"
                    for node in independent_state_oracle._parse_adjacency(
                        adjacency
                    ).values()
                )
            carbon += carbon_by_graph[adjacency]
        carbon_by_event[record.event_id] = carbon
    full_index = SiteIndex(state, records)
    impossible = tuple(
        record
        for record in records
        if carbon_by_event[record.event_id] > initial_ledger["C"]
    )
    assert all(not full_index.candidates(record.event_id) for record in impossible)
    possible = tuple(
        record
        for record in records
        if carbon_by_event[record.event_id] <= initial_ledger["C"]
    )
    print(
        f"R0 full inventory admitted={len(records)}; stoichiometrically possible={len(possible)}; excluded only by conserved carbon budget={len(impossible)}",
        flush=True,
    )
    engine = IsothermalSSA(
        state,
        possible,
        temperature=800.0,
        volume=1e-50,
        rng=np.random.default_rng(35_000_001),
    )
    result = engine.run(max_events=10_000)
    assert len(result.events) == 10_000
    fired = {event.event_id for event in result.events}
    assert fired & {record.event_id for record in homolysis}
    assert result.events[0].event_id in {record.event_id for record in homolysis}
    assert result.irreversible_fraction == sum(
        event.irreversible for event in result.events
    ) / len(result.events)
    assert state.total_ledger == initial_ledger
    assert not state.sink
    print(
        f"REAL-SET-SMOKE R0 pristine: events={len(result.events)} homolysis={sum(event.event_id in {record.event_id for record in homolysis} for event in result.events)} irreversible={result.irreversible_fraction:.12g} seed=35000001"
    )


@slow
def test_r0_melt_density_mixed_pristine_chains_fire_non_inverse_chemistry(
    real_ps_artifact, independent_state_oracle
):
    records = tuple(
        EventRecord.from_dict(data)
        for data in real_ps_artifact["records"]
        if data["inventory_class"] != "R1:J_ring"
    )
    pristine = {}
    for record in records:
        if (
            record.family == "R_Recombination"
            and record.arity == 1
            and record.radical_delta == 2
            and record.participant_site_types == ["pristine"]
        ):
            graph = independent_state_oracle._parse_adjacency(record.reactant_graphs[0])
            carbon = sum(node["element"] == "C" for node in graph.values())
            pristine.setdefault(carbon // 8, (record, graph))
    lengths = (3, 5, 3, 5, 3)
    assert set(lengths) <= pristine.keys()
    strands = []
    for ordinal, length in enumerate(lengths):
        record, graph = pristine[length]
        assert sum(node["element"] == "C" for node in graph.values()) == 8 * length
        backbone = {
            atom
            for atom, node in graph.items()
            if node["element"] == "C"
            and all(float(order) == 1.0 for order in node.get("edges", {}).values())
        }
        endpoints = [
            atom
            for atom in backbone
            if len(set(graph[atom].get("edges", {})) & backbone) == 1
        ]
        assert len(endpoints) == 2
        ordered = [min(endpoints)]
        while len(ordered) < len(backbone):
            following = set(graph[ordered[-1]].get("edges", {})) & backbone - set(
                ordered
            )
            assert len(following) == 1
            ordered.append(following.pop())
        assert len(ordered) == 2 * length
        positions = {atom: position for position, atom in enumerate(ordered)}
        roles = {atom: "backbone" for atom in backbone}
        for atom in graph:
            visited = {atom}
            pending = {atom}
            while not pending & backbone:
                pending = {
                    neighbor
                    for current in pending
                    for neighbor in graph[current].get("edges", {})
                    if neighbor not in visited
                }
                assert pending
                visited |= pending
            anchors = pending & backbone
            assert len(anchors) == 1
            roles[(atom, "position")] = positions[anchors.pop()]
        chain, _, _ = independent_state_oracle._harness(
            record, roles=roles, namespace=f"melt-{ordinal}-"
        )
        assert len(chain.strands) == 1
        strand = next(iter(chain.strands.values()))
        strand.length = len(ordered)
        strands.append(strand)
    state = KMCState(strands)
    assert len(state.strands) == 5
    assert len({strand.length for strand in state.strands.values()}) == 2
    assert state.total_radicals == 0
    initial_ledger = state.total_ledger
    carbon = initial_ledger["C"]
    repeat_equivalents = carbon / 8.0
    density_g_cm3 = 1.05
    repeat_molar_mass_g_mol = 104.15
    mass_g = repeat_equivalents * repeat_molar_mass_g_mol / N_A
    volume_m3 = mass_g / density_g_cm3 * 1e-6
    print(
        f"R0 MELT SETUP chain_units={lengths} carbon={carbon} "
        f"repeat_equivalents=C/8={repeat_equivalents:g} "
        f"mass_g=(C/8)*104.15/N_A={mass_g:.12g} "
        f"volume_m3=mass_g/1.05*1e-6={volume_m3:.12g} "
        "density_g_cm3=1.05 temperature_K=800 seed=35000002; "
        "PS repeat-equivalent density neglects terminal end-cap mass",
        flush=True,
    )
    engine = IsothermalSSA(
        state,
        records,
        temperature=800.0,
        volume=volume_m3,
        rng=np.random.default_rng(35_000_002),
    )
    by_id = {record.event_id: record for record in records}
    histogram = Counter()
    qualifying = Counter()
    previous_id = None
    homolyses = 0
    for ordinal in range(10_000):
        before = state.total_ledger, state.total_radicals, state.total_mass
        event = engine.step()
        if event is None:
            break
        record = by_id[event.event_id]
        independent_state_oracle._assert_record_balance(state, record, *before)
        independent_state_oracle._assert_strand_radicals_match_graph(state)
        assert state.total_ledger == initial_ledger
        histogram[record.family] += 1
        homolyses += record.family == "R_Recombination" and record.radical_delta > 0
        beta_scission = (
            record.family == "R_Addition_MultipleBond"
            and record.arity == 1
            and any(operation["action"] == "break" for operation in record.bond_ops)
        )
        if record.reverse_of != previous_id and (
            record.family in {"H_Abstraction", "Disproportionation"} or beta_scission
        ):
            qualifying[record.family] += 1
        previous_id = record.event_id
        if (ordinal + 1) % 1_000 == 0:
            print(
                f"R0 MELT PROGRESS events={ordinal + 1} "
                f"families={dict(sorted(histogram.items()))} "
                f"non_inverse={dict(sorted(qualifying.items()))}",
                flush=True,
            )
    fired = len(engine.events)
    irreversible_fraction = engine.irreversible_fired / fired if fired else 0.0
    print(
        f"R0 MELT RESULT events={fired} homolyses={homolyses} "
        f"families={dict(sorted(histogram.items()))} "
        f"non_inverse={dict(sorted(qualifying.items()))} "
        f"irreversible_fraction={irreversible_fraction:.12g} "
        f"ledger={state.total_ledger} sink={state.sink} "
        f"time_s={engine.time:.12g}",
        flush=True,
    )
    assert fired == 10_000
    assert homolyses > 0
    assert (
        qualifying
    ), "no non-inverse chemistry fired within the fixed 10,000-event budget"
    assert not state.sink
    assert engine.bound_violations == 0


def test_production_thinning_sampler_counts_unordered_pairs_and_checks_bound():
    items = [
        PairItem(f"cluster-{index}", f"cluster-{index}", 1, "cluster")
        for index in range(8)
    ]
    sampler = LengthBinnedThinningSampler(
        items,
        ConstantPairKernel(3.0),
        (("cluster", "cluster"),),
        volume=2.0,
    )
    assert len(sampler.cells) == 1
    assert sampler.cells[0].pair_count == 8 * 7 // 2
    assert sampler.total_bound_propensity == 3.0 * 28 / 2.0
    proposal = sampler.propose(np.random.default_rng(13_004))
    assert proposal is not None and proposal.accepted
    assert proposal.left.item_id != proposal.right.item_id
    assert sampler.bound_checks == 1
    assert sampler.bound_violations == 0

    broken = LengthBinnedThinningSampler(
        items,
        SumPairKernel(2.0, bound_scale=0.999),
        (("cluster", "cluster"),),
        volume=2.0,
    )
    with pytest.raises(ThinningBoundError):
        broken.propose(np.random.default_rng(13_005))


def _run_coagulation(kernel, target_times, seed, initial_count=50_000):
    items = (
        PairItem(f"initial-{index}", f"initial-{index}", 1, "cluster")
        for index in range(initial_count)
    )
    sampler = LengthBinnedThinningSampler(
        items,
        kernel,
        (("cluster", "cluster"),),
        volume=float(initial_count),
    )
    rng = np.random.default_rng(seed)
    now = 0.0
    serial = 0
    pending = sampler.propose(rng)
    event_time = math.inf if pending is None else pending.elapsed
    observations = []
    for target in target_times:
        while event_time <= target:
            now = event_time
            if pending.accepted:
                serial += 1
                sampler.replace_pair(
                    pending.left.item_id,
                    pending.right.item_id,
                    PairItem(
                        f"merged-{serial}",
                        f"merged-{serial}",
                        pending.left.length + pending.right.length,
                        "cluster",
                    ),
                )
            pending = sampler.propose(rng)
            event_time = math.inf if pending is None else now + pending.elapsed
        counts = Counter(item.length for item in sampler.items.values())
        observations.append(
            (len(sampler.items), tuple(counts.get(size, 0) for size in range(1, 6)))
        )
    assert sampler.bound_violations == 0
    return observations, sampler.bound_checks


def _run_coagulation_case(arguments):
    return _run_coagulation(*arguments)


def _constant_solution(initial_count, time, rate=2.0):
    # Our convention uses unordered component pairs, giving
    # tau=K*N0*t/(2V); V=N0 and K=2 in this verifier.
    tau = rate * initial_count * time / (2.0 * initial_count)
    total = initial_count / (1.0 + tau)
    sizes = tuple(
        initial_count * tau ** (size - 1) / (1.0 + tau) ** (size + 1)
        for size in range(1, 6)
    )
    return total, sizes


def _sum_solution(initial_count, time, coefficient=1.0):
    # With K(i,j)=b(i+j), the exponential generating function obeys the
    # rooted-tree equation.  Lagrange inversion gives the Borel cluster-number
    # law p_k=exp(-q*k)*(q*k)^(k-1)/k!, q=1-exp(-g).  Multiplying p_k by
    # N=N0*exp(-g) yields n_k in this unordered-pair convention.
    g = coefficient * initial_count * time / initial_count
    fraction = math.exp(-g)
    q = 1.0 - fraction
    total = initial_count * fraction
    sizes = tuple(
        total * math.exp(-q * size) * (q * size) ** (size - 1) / math.factorial(size)
        for size in range(1, 6)
    )
    return total, sizes


@slow
@pytest.mark.parametrize(
    ("name", "kernel", "solution", "target_times"),
    [
        (
            "constant",
            ConstantPairKernel(2.0),
            _constant_solution,
            (1.0, 1.5, 7.0 / 3.0, 4.0),
        ),
        (
            "sum",
            SumPairKernel(1.0),
            _sum_solution,
            tuple(-math.log(fraction) for fraction in (0.5, 0.4, 0.3, 0.2)),
        ),
    ],
)
def test_c4_production_thinning_matches_smoluchowski(
    name, kernel, solution, target_times
):
    initial_count = 50_000
    replicas = 64 if name == "sum" else 32
    relative_tolerance = 0.01
    sigma_tolerance = 4.0
    seed = 4_013_700
    arguments = [
        (kernel, target_times, seed + replica, initial_count)
        for replica in range(replicas)
    ]
    with multiprocessing.get_context("fork").Pool(processes=8) as pool:
        runs = pool.map(_run_coagulation_case, arguments)
    checks = sum(run[1] for run in runs)
    values = np.asarray(
        [[total, *sizes] for run, _ in runs for total, sizes in run], dtype=float
    ).reshape(replicas, len(target_times), 6)
    rows = []
    for time_index, target in enumerate(target_times):
        expected_total, expected_sizes = solution(initial_count, target)
        expected = np.asarray((expected_total, *expected_sizes))
        mean = values[:, time_index, :].mean(axis=0)
        standard_error = values[:, time_index, :].std(axis=0, ddof=1) / math.sqrt(
            replicas
        )
        relative = np.abs(mean / expected - 1.0)
        assert np.all(relative <= relative_tolerance)
        assert np.all(np.abs(mean - expected) <= sigma_tolerance * standard_error)
        rows.append((target, float(relative.max())))
    print(f"C4 {name}: replicas={replicas} checks={checks} max-errors={rows}")


def _normalized_expected(weights, draws):
    values = np.asarray(weights, dtype=float)
    return draws * values / values.sum()


@slow
def test_c5_heterogeneous_component_pair_sampling():
    draws = 100_000
    seed = 5_013_700
    kinds = {
        (name, position): (name, position)
        for name in "ABC"
        for position in ("end", "mid")
    }
    raw_items = [
        ("c0-A", "c0", 2, kinds["A", "end"], 2),
        ("c1-A", "c1", 3, kinds["A", "end"], 1),
        ("c2-A", "c2", 5, kinds["A", "end"], 1),
        ("c3-A", "c3", 8, kinds["A", "end"], 1),
        ("c0-B", "c0", 2, kinds["B", "mid"], 1),  # biradical c0
        ("c4-B", "c4", 4, kinds["B", "mid"], 1),
        ("c5-B", "c5", 6, kinds["B", "mid"], 1),
        ("c6-B", "c6", 12, kinds["B", "mid"], 1),
        ("c7-C", "c7", 2, kinds["C", "end"], 1),
        ("c8-C", "c8", 5, kinds["C", "end"], 1),
        ("c9-C", "c9", 9, kinds["C", "end"], 1),
    ]
    items = [PairItem(*item) for item in raw_items]
    activation = {
        ("A", "A"): 1.1e7,
        ("A", "B"): 1.3e7,
        ("A", "C"): 1.7e7,
        ("B", "C"): 2.1e7,
        ("C", "C"): 2.3e7,
    }
    kernel = METPairKernel(
        activation,
        TRANSPORT_ARMS["A0_REF_H_CROSS_NEc"],
        700.0,
        spin_factor=1.0,
    )
    pair_kinds = (
        (kinds["A", "end"], kinds["A", "end"]),
        (kinds["A", "end"], kinds["B", "mid"]),
        (kinds["A", "end"], kinds["C", "end"]),
        (kinds["B", "mid"], kinds["C", "end"]),
        (kinds["C", "end"], kinds["C", "end"]),
    )
    sampler = LengthBinnedThinningSampler(
        items, kernel, pair_kinds, volume=1.0, normalization=1.0
    )
    assert {cell.left_bin for cell in sampler.cells} | {
        cell.right_bin for cell in sampler.cells
    } >= {1, 2, 3}
    assert all(cell.left_bin != 4 and cell.right_bin != 4 for cell in sampler.cells)
    assert any(
        cell.left_kind == kinds["A", "end"]
        and cell.left_bin == 3
        and len(
            [item for item in items if item.kind == cell.left_kind and item.length == 8]
        )
        == 1
        for cell in sampler.cells
    )
    assert all(
        left.component_id != right.component_id
        for cell in sampler.cells
        for left, right in sampler.cell_pairs(cell)
    )

    observed = Counter()
    rng = np.random.default_rng(seed)
    for _ in range(draws):
        proposal = sampler.draw_accepted(rng)
        assert proposal is not None
        observed[proposal.cell.cell_id] += 1
    exact = sampler.exact_cell_propensities()
    cell_ids = [cell.cell_id for cell in sampler.cells]
    actual = np.asarray([observed[cell_id] for cell_id in cell_ids])
    expected = _normalized_expected([exact[cell_id] for cell_id in cell_ids], draws)
    correct_p = float(chisquare(actual, expected).pvalue)
    assert correct_p > 1e-3

    site_pair_weights = []
    internal_pair_weights = []
    identical_twice_weights = []
    for cell in sampler.cells:
        correct = exact[cell.cell_id]
        site_pair_weights.append(
            sum(
                kernel.exact_rate(left, right) * int(left.payload) * int(right.payload)
                for left, right in sampler.cell_pairs(cell)
            )
        )
        identical_twice_weights.append(
            correct / 2.0 if cell.left_kind == cell.right_kind else correct
        )
        internal = correct
        if (
            cell.left_kind == kinds["A", "end"]
            and cell.right_kind == kinds["B", "mid"]
            and cell.left_bin == 1
            and cell.right_bin == 1
        ):
            internal += kernel.exact_rate(items[0], items[4])
        internal_pair_weights.append(internal)
    control_ps = {
        "site-pair": float(
            chisquare(actual, _normalized_expected(site_pair_weights, draws)).pvalue
        ),
        "identical-factor-twice": float(
            chisquare(
                actual, _normalized_expected(identical_twice_weights, draws)
            ).pvalue
        ),
        "biradical-internal": float(
            chisquare(actual, _normalized_expected(internal_pair_weights, draws)).pvalue
        ),
    }
    assert all(value < 1e-3 for value in control_ps.values())
    print(
        f"C5: draws={draws} seed={seed} cells={len(cell_ids)} "
        f"p={correct_p:.6g} controls={control_ps}"
    )


@slow
def test_met_bound_ramp_all_arms_classes_spins_and_length_pairs():
    temperature = 700.0
    refresh_delta = 25.0
    maximum_length = 320
    activation = {("radical", "radical"): 1.0e8}
    kind_by_class = {"end": ("radical", "end"), "mid": ("radical", "mid")}
    bins = []
    lower = 1
    while lower <= maximum_length:
        upper = min(2 * lower - 1, maximum_length)
        bins.append((lower, upper))
        lower *= 2

    capture_transition = 6.0 * (SIGMA_CONTACT / 2.0) ** 2 / (PS_C_R2 * PS_M0)
    assert 1.0 < capture_transition < maximum_length
    cross_transitions = [
        arm.n_e for arm in TRANSPORT_ARMS.values() if arm.n_e is not None
    ]
    assert min(cross_transitions) > 1.0
    assert max(cross_transitions) < maximum_length

    # Every D0 package is monotone on the complete refresh interval.  For
    # fixed i,j the remaining factors are temperature-independent, and
    # Collins--Kimball is monotone in k_D and k_act.  Thus endpoint maxima
    # prove the bound for every intervening temperature (the activation bound
    # is supplied separately, as production does at a refresh).
    temperature_grid = np.linspace(temperature, temperature + refresh_delta, 101)
    for arm in TRANSPORT_ARMS.values():
        d0_values = np.asarray([arm.d0(value) for value in temperature_grid])
        assert np.all(np.diff(d0_values) >= 0.0)

    checks = 0
    for arm in TRANSPORT_ARMS.values():
        for spin_factor in (0.25, 1.0):
            for pair_class in sorted(PAIR_CLASSES):
                left_class, right_class = pair_class.split("/")
                left_kind = kind_by_class[left_class]
                right_kind = kind_by_class[right_class]
                kernel = METPairKernel(
                    activation,
                    arm,
                    temperature + refresh_delta,
                    spin_factor=spin_factor,
                    bound_temperature=temperature + refresh_delta,
                    bound_activation_rates=activation,
                )
                for left_bin_index, left_range in enumerate(bins):
                    for right_range in bins[left_bin_index:]:
                        bound = kernel.bound_rate(
                            left_kind, right_kind, left_range, right_range
                        )
                        for left_length in range(left_range[0], left_range[1] + 1):
                            for right_length in range(
                                right_range[0], right_range[1] + 1
                            ):
                                exact = kernel.exact_rate(
                                    PairItem("left", "left", left_length, left_kind),
                                    PairItem(
                                        "right", "right", right_length, right_kind
                                    ),
                                )
                                assert exact <= bound + 8.0 * math.ulp(
                                    max(exact, bound, 1.0)
                                )
                                checks += 1

                # The corner maximum is valid here: for i<=j, D(j)<=D(i);
                # multiplying (D(i)+D(j)) by max(sigma,c*sqrt(i)) has a
                # non-positive derivative in i for both N^-1 and N^-2 laws,
                # and increasing j only lowers D(j).  A dense discrete scan
                # across both transitions independently checks that proof.
                diagonal = [
                    kernel.exact_rate(
                        PairItem("left", "left", length, left_kind),
                        PairItem("right", "right", length, right_kind),
                    )
                    for length in range(1, maximum_length + 1)
                ]
                assert np.all(np.diff(diagonal) <= 0.0)
    print(
        f"MET-BOUND-RAMP: T=[{temperature},{temperature + refresh_delta}] "
        f"bins={bins} checks={checks} arms={len(TRANSPORT_ARMS)}"
    )


@slow
@pytest.mark.parametrize(("k1", "k2", "horizon"), [(1.0, 0.2, 2.0), (0.5, 2.0, 1.0)])
def test_leak_toy_replica_ratio_matches_analytic_value(k1, k2, horizon):
    replicas = 20_000
    seed = 6_013_700
    sigma_tolerance = 4.0
    results = [
        _run_leak(k1, k2, horizon, seed + replica) for replica in range(replicas)
    ]
    numerators = np.asarray([result.leak_num for result in results])
    denominators = np.asarray([result.leak_den for result in results])
    actual = float(numerators.sum() / denominators.sum())
    enabled_integral = 1.0 - math.exp(-k1 * horizon)
    refused_integral = k2 * (horizon - (1.0 - math.exp(-k1 * horizon)) / k1)
    expected = refused_integral / (enabled_integral + refused_integral)
    centered = numerators - expected * denominators
    standard_error = float(
        centered.std(ddof=1) / math.sqrt(replicas) / denominators.mean()
    )
    assert abs(actual - expected) <= sigma_tolerance * standard_error
    if (k1, k2, horizon) == (1.0, 0.2, 2.0):
        assert 0.05 < expected < 0.5
    print(
        f"Leak toy k1={k1} k2={k2} T={horizon}: replicas={replicas} "
        f"L={actual:.12g} analytic={expected:.12g} se={standard_error:.3g}"
    )


class _Fenwick:
    def __init__(self, values, capacity):
        self.tree = np.zeros(capacity + 1, dtype=np.int64)
        self.size = len(values)
        for index, value in enumerate(values, start=1):
            self.add(index, int(value))

    def add(self, index, delta):
        while index < len(self.tree):
            self.tree[index] += int(delta)
            index += index & -index

    @property
    def total(self):
        value = 0
        index = self.size
        while index:
            value += int(self.tree[index])
            index -= index & -index
        return value

    def append(self, value):
        self.size += 1
        self.add(self.size, int(value))
        return self.size

    def find(self, ordinal):
        index = 0
        bit = 1 << (self.size.bit_length() - 1)
        remaining = int(ordinal)
        while bit:
            candidate = index + bit
            if candidate <= self.size and int(self.tree[candidate]) <= remaining:
                remaining -= int(self.tree[candidate])
                index = candidate
            bit >>= 1
        return index + 1


def _random_scission_record():
    reactant = "\n".join(
        (
            "1 *1 C u0 p0 c0 {2,S}",
            "2 *2 C u0 p0 c0 {1,S}",
        )
    )
    return EventRecord(
        family="random_scission",
        arity=1,
        participant_site_types=["backbone_bond"],
        reactant_multiplicities=[1],
        rate_order=1,
        rate_units="s^-1",
        k_table={
            "T": [600.0, 800.0],
            "k": [1.0, 1.0],
            "interpolation": "linear-ln-k",
            "extrapolation": "refuse",
        },
        atom_map={0: 0, 1: 1},
        bond_ops=[{"action": "break", "atoms": [0, 1], "order": "1.0"}],
        radical_delta=0,
        status="irreversible",
        status_reason="synthetic isolated random-scission verifier",
        reactant_graphs=[reactant],
        product_graphs=["1 *1 C u0 p0 c0", "1 *2 C u0 p0 c0"],
    )


def _apply_executor_scission(state, record, strand_id, cut, serial):
    strand = state.strands[strand_id]
    left_ref = AtomRef(f"cut-{serial}-L", strand_id, cut, "backbone")
    right_ref = AtomRef(f"cut-{serial}-R", strand_id, cut + 1, "backbone")
    strand.atom_refs = {0: left_ref, 1: right_ref}
    strand.atom_graph = {
        left_ref.uuid: {
            "element": "C",
            "radical": 0,
            "charge": 0,
            "lone_pairs": 0,
            "implicit_hydrogens": 0,
            "edges": {right_ref.uuid: "1.0"},
        },
        right_ref.uuid: {
            "element": "C",
            "radical": 0,
            "charge": 0,
            "lone_pairs": 0,
            "implicit_hydrogens": 0,
            "edges": {left_ref.uuid: "1.0"},
        },
    }
    state._issued_uuids.update((left_ref.uuid, right_ref.uuid))
    graph = {
        0: {
            "element": "C",
            "radical": 0,
            "charge": 0,
            "lone_pairs": 0,
            "implicit_hydrogens": 0,
            "edges": {1: "1.0"},
            "atom_ref": left_ref,
            "position": cut,
            "role": "backbone",
        },
        1: {
            "element": "C",
            "radical": 0,
            "charge": 0,
            "lone_pairs": 0,
            "implicit_hydrogens": 0,
            "edges": {0: "1.0"},
            "atom_ref": right_ref,
            "position": cut + 1,
            "role": "backbone",
        },
    }
    entry = state.apply(record, (Site("backbone_bond", strand_id, cut, 0, graph),))
    assert entry.derived_cut_offsets == (cut,)
    affected = set(entry.inheritance)
    for affected_id in affected:
        state.strands[affected_id].atom_refs = {}
        state.strands[affected_id].atom_graph = {}
    new_id = next(item for item in affected if item != strand_id)
    return new_id


def _schulz_zimm_lengths(rng, chain_count, dispersity):
    dpn = 100_000.0 / 104.15
    support = np.arange(1, 40_001, dtype=float)
    shape = 1.0 / (dispersity - 1.0)
    logp = (shape - 1.0) * np.log(support) - support * shape / dpn - math.lgamma(shape)
    probabilities = np.exp(logp - logp.max())
    probabilities /= probabilities.sum()
    # Randomized stratified sampling is still a sample from the report's
    # discrete Schulz--Zimm law, but resolves its long tail much better than
    # iid draws at a tractable executor population size.
    quantiles = (np.arange(chain_count, dtype=float) + rng.random()) / chain_count
    return support[np.searchsorted(np.cumsum(probabilities), quantiles)].astype(int)


def _exact_scission_moments(initial_lengths, conversion):
    lengths = np.asarray(initial_lengths, dtype=float)
    q = 1.0 - conversion
    if conversion == 0.0:
        pair_sum = lengths * (lengths - 1.0) / 2.0
    else:
        pair_sum = (lengths * q * (1.0 - q) - q * (1.0 - q**lengths)) / (1.0 - q) ** 2
    mu0 = np.sum(1.0 + (lengths - 1.0) * conversion)
    mu1 = np.sum(lengths)
    mu2 = np.sum(lengths + 2.0 * pair_sum)
    return mu0, mu1, mu2


def _mw_dispersity(lengths):
    lengths = np.asarray(lengths, dtype=float)
    mu0 = float(len(lengths))
    mu1 = float(lengths.sum())
    mu2 = float(np.dot(lengths, lengths))
    return 104.15 * mu2 / mu1 / 1e3, mu2 * mu0 / mu1**2


def _run_random_scission_case(arguments):
    dispersity, replica, chain_count, conversions, seed = arguments
    rng = np.random.default_rng(seed + replica)
    initial = _schulz_zimm_lengths(rng, chain_count, dispersity)
    states = [KMCState((Strand("chain", int(length)),)) for length in initial]
    record = _random_scission_record()
    ids = [(index, "chain") for index in range(chain_count)]
    weights = [int(length) - 1 for length in initial]
    fenwick = _Fenwick(weights, int(initial.sum()) + chain_count + 1)
    targets = [-math.log1p(-value) for value in conversions]
    now = 0.0
    serial = 0
    pending = float(rng.exponential(1.0 / fenwick.total))
    observations = []
    for target in targets:
        while now + pending <= target:
            now += pending
            remaining = fenwick.total
            selected_index = fenwick.find(int(rng.integers(remaining)))
            state_index, strand_id = ids[selected_index - 1]
            state = states[state_index]
            old_weight = state.strands[strand_id].length - 1
            cut = int(rng.integers(old_weight))
            serial += 1
            new_id = _apply_executor_scission(state, record, strand_id, cut, serial)
            left_weight = state.strands[strand_id].length - 1
            right_weight = state.strands[new_id].length - 1
            fenwick.add(selected_index, left_weight - old_weight)
            ids.append((state_index, new_id))
            fenwick.append(right_weight)
            pending = (
                math.inf
                if fenwick.total == 0
                else float(rng.exponential(1.0 / fenwick.total))
            )
        actual = _mw_dispersity(
            [
                strand.length
                for molecule_state in states
                for strand in molecule_state.strands.values()
            ]
        )
        exact_mu = _exact_scission_moments(initial, conversions[len(observations)])
        expected = (
            104.15 * exact_mu[2] / exact_mu[1] / 1e3,
            exact_mu[2] * exact_mu[0] / exact_mu[1] ** 2,
        )
        observations.append((actual, expected))
    return observations, serial


@slow
@pytest.mark.parametrize("dispersity", [2.0, 1.05])
def test_random_scission_executor_matches_exact_population_balance(dispersity):
    dpn = 100_000.0 / 104.15
    conversions = (1e-3, 1.0 / (dpn - 1.0), 3.0 / (dpn - 1.0), 1e-2)
    chain_count = 5_000
    replicas = 16
    seed = 9_013_700 + int(100 * dispersity)
    arguments = [
        (dispersity, replica, chain_count, conversions, seed)
        for replica in range(replicas)
    ]
    with multiprocessing.get_context("fork").Pool(processes=8) as pool:
        runs = pool.map(_run_random_scission_case, arguments)
    actual = np.asarray(
        [[*observation[0]] for run, _ in runs for observation in run]
    ).reshape(replicas, len(conversions), 2)
    expected = np.asarray(
        [[*observation[1]] for run, _ in runs for observation in run]
    ).reshape(replicas, len(conversions), 2)
    differences = actual - expected
    standard_error = differences.std(axis=0, ddof=1) / math.sqrt(replicas)
    assert np.all(np.abs(differences.mean(axis=0)) <= 3.0 * standard_error)
    printed = {
        2.0: {
            conversions[1]: (99.92, 1.998),
            conversions[2]: (49.90, 1.996),
            1e-2: (18.78, 1.989),
        },
        1.05: {
            conversions[1]: (75.34, None),
            conversions[2]: (45.73, None),
        },
    }[dispersity]
    rows = []
    for index, conversion in enumerate(conversions):
        mean = actual[:, index, :].mean(axis=0)
        spread_error = actual[:, index, :].std(axis=0, ddof=1) / math.sqrt(replicas)
        if conversion in printed:
            mw_printed, dispersity_printed = printed[conversion]
            assert abs(mean[0] - mw_printed) <= 3.0 * spread_error[0]
            if dispersity_printed is not None:
                assert abs(mean[1] - dispersity_printed) <= 3.0 * spread_error[1]
        rows.append((conversion, *mean, *expected[:, index, :].mean(axis=0)))
    print(
        f"Random scission D0={dispersity}: chains={chain_count} "
        f"replicas={replicas} events={sum(run[1] for run in runs)} rows={rows}"
    )
