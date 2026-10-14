"""Observation independence of the SSA clock (I-126).

A schedule of ``run(until=...)`` calls must reproduce the trajectory of one run to the
last horizon: same event times, selections, final state, RNG state and integrated
diagnostics.  Fixtures are hand-authored (fixed atom ids, no database, no uuid use).
"""

import hashlib
from copy import deepcopy
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

from rmgpy.kmc.event_record import EventRecord
from rmgpy.kmc.ssa import IsothermalSSA
from rmgpy.kmc.state import AtomRef, KMCState, Strand

SEEDS = (11, 12, 13)
HORIZON = 6.0


def _strand(strand_id, atom_uuid, length):
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


def _toggle(k):
    return EventRecord(
        family="toy",
        arity=1,
        participant_site_types=["radical"],
        reactant_multiplicities=[1],
        rate_order=1,
        rate_units="s^-1",
        k_table={
            "T": [600.0, 800.0],
            "k": [k, k],
            "interpolation": "linear-ln-k",
            "extrapolation": "refuse",
        },
        bond_ops=[{"action": "set_radical", "atom": 0, "value": 0}],
        radical_delta=-1,
        reactant_graphs=["1 *1 C u1 p0 c0"],
        product_graphs=["1 *1 C u0 p0 c0"],
    )


def _refused(k):
    return EventRecord(
        family="toy-refused",
        arity=1,
        participant_site_types=["B"],
        reactant_multiplicities=[1],
        rate_order=1,
        rate_units="s^-1",
        k_table={
            "T": [600.0, 800.0],
            "k": [k, k],
            "interpolation": "linear-ln-k",
            "extrapolation": "refuse",
        },
        status="refused",
        status_reason="toy leak frontier",
        reactant_graphs=["1 *1 C u0 p0 c0"],
        product_graphs=["1 *1 C u0 p0 c0"],
    )


def _records(reversible=True):
    forward = _toggle(2.0)
    if not reversible:
        return (forward, _refused(0.3))
    reverse = replace(
        forward,
        event_id="",
        participant_site_types=["pristine"],
        radical_delta=1,
        bond_ops=[{"action": "set_radical", "atom": 0, "value": 1}],
        reactant_graphs=forward.product_graphs,
        product_graphs=forward.reactant_graphs,
        k_table={**forward.k_table, "k": [0.5, 8.0]},
    )
    return (forward, reverse, _refused(0.3))


def _state():
    return KMCState(tuple(_strand(f"chain-{i}", f"carbon-{i}", 5 + i) for i in range(3)))


def _engine(seed, *, reversible=True, rng=None, **kwargs):
    return IsothermalSSA(
        _state(),
        _records(reversible),
        temperature=700.0,
        volume=1.0,
        rng=rng if rng is not None else np.random.default_rng(seed),
        **kwargs,
    )


def _snapshot(engine):
    """Everything an observer could compare: trace, state, RNG and diagnostics."""
    return {
        "trace": [
            (event.time.hex(), event.event_id, event.state_hash, event.irreversible)
            for event in engine.events
        ],
        "applied": [
            (entry.event_id, entry.participants, entry.derived_cut_offsets)
            for entry in engine.state.event_log
        ],
        "signature": engine.state.state_hash(),
        "rng": engine.rng.bit_generator.state,
        "time": engine.time.hex(),
        "leak_num": engine.leak_num.hex(),
        "leak_den": engine.leak_den.hex(),
        "irreversible_fired": engine.irreversible_fired,
        "bound_checks": engine.bound_checks,
        "bound_violations": engine.bound_violations,
    }


def _run_schedule(engine, horizons, check_each=None):
    for horizon in horizons:
        engine.run(until=horizon)
        if check_each is not None:
            check_each(engine, horizon)


def _reference_events(seed, horizon=HORIZON):
    engine = _engine(seed)
    engine.run(until=horizon)
    return [event.time for event in engine.events]


def _schedules(seed):
    times = _reference_events(seed)
    assert len(times) > 8, "fixture must fire enough events to be a test"
    rng = np.random.default_rng(1000 + seed)
    uneven = sorted(float(x) for x in rng.uniform(0.0, HORIZON, 96)) + [HORIZON]
    return {
        "one": [HORIZON],
        "ten": [HORIZON * i / 10 for i in range(1, 11)],
        "ninety_seven_uneven": uneven,
        "several_before_first_event": [times[0] * f for f in (0.2, 0.4, 0.6, 0.8)]
        + [HORIZON],
        "exactly_at_event_times": [*times, HORIZON],
        "just_after_each_event": [float(np.nextafter(t, np.inf)) for t in times]
        + [HORIZON],
        "just_before_each_event": [float(np.nextafter(t, 0.0)) for t in times]
        + [HORIZON],
    }


@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize(
    "name",
    (
        "one",
        "ten",
        "ninety_seven_uneven",
        "several_before_first_event",
        "exactly_at_event_times",
        "just_after_each_event",
        "just_before_each_event",
    ),
)
def test_any_horizon_schedule_reproduces_the_single_run_trace(seed, name):
    horizons = _schedules(seed)[name]
    single = _engine(seed)
    single.run(until=HORIZON)
    piecewise = _engine(seed)
    _run_schedule(piecewise, horizons)
    assert _snapshot(piecewise) == _snapshot(single)


@pytest.mark.parametrize("seed", SEEDS)
def test_diagnostics_agree_with_a_fresh_run_at_every_common_horizon(seed):
    horizons = _schedules(seed)["ninety_seven_uneven"]
    piecewise = _engine(seed)
    for horizon in horizons:
        piecewise.run(until=horizon)
        fresh = _engine(seed)
        fresh.run(until=horizon)
        assert _snapshot(piecewise) == _snapshot(fresh), horizon


def test_single_uninterrupted_run_is_byte_identical_to_the_base_engine():
    """Digest recorded from the pre-fix engine: consumption order must not change."""
    engine = _engine(11)
    engine.run(until=HORIZON)
    digest = hashlib.sha256(repr(_snapshot(engine)).encode()).hexdigest()
    assert digest == "cf51e5ad2098f7db6cc8ef90e6aa2e5e4fdf028d606549603ab920419ca2c660"


def test_run_to_the_same_horizon_twice_is_a_no_op():
    engine = _engine(11)
    engine.run(until=2.0)
    before = _snapshot(engine)
    engine.run(until=2.0)
    assert _snapshot(engine) == before


class _CountingRng:
    """Delegating generator that counts waiting-time draws."""

    def __init__(self, seed):
        self.inner = np.random.default_rng(seed)
        self.exponentials = 0

    def exponential(self, *args, **kwargs):
        self.exponentials += 1
        return self.inner.exponential(*args, **kwargs)

    def __getattr__(self, name):
        return getattr(self.inner, name)


def _engine_with_pending_draw(seed=11):
    rng = _CountingRng(seed)
    engine = _engine(seed, rng=rng)
    first_event = _reference_events(seed)[0]
    engine.run(until=first_event * 0.5)
    assert not engine.events and rng.exponentials == 1
    return engine, rng, first_event


def test_valid_pending_draw_is_reused_without_redrawing():
    engine, rng, first_event = _engine_with_pending_draw()
    engine.run(until=first_event)
    assert engine.events[0].time == first_event
    # one draw for the first proposal, then one for the proposal after the event
    assert rng.exponentials == 2


def test_state_mutation_between_runs_discards_the_pending_draw():
    engine, rng, first_event = _engine_with_pending_draw()
    record = next(r for r in engine.records if r.participant_site_types == ["radical"])
    participants = engine.index.candidates(record.event_id)[0]
    engine.index.apply(record, participants)  # outside step
    engine.run(until=10.0)
    assert rng.exponentials >= 2  # the stale draw was replaced
    assert engine.events[0].time != first_event


def test_direct_state_edit_between_runs_discards_the_pending_draw():
    engine, rng, first_event = _engine_with_pending_draw()
    engine.state.sink["X"] = 1  # not an index event; still a state change
    engine.run(until=10.0)
    assert rng.exponentials >= 2
    assert engine.events[0].time != first_event


@pytest.mark.parametrize("change", ("temperature", "volume"))
def test_environment_change_between_runs_redraws_from_the_horizon(change):
    engine, rng, first_event = _engine_with_pending_draw()
    horizon = engine.time
    setattr(engine, change, 760.0 if change == "temperature" else 2.0)
    expected_rng = deepcopy(engine.rng.inner)
    total = engine.channel_propensities().total_enabled
    expected = horizon + expected_rng.exponential(1.0 / total)
    engine.run(until=10.0)
    assert rng.exponentials >= 2
    assert engine.events[0].time == expected


def test_invalidation_integrates_the_stale_interval_with_the_old_rates():
    engine, rng, first_event = _engine_with_pending_draw()
    horizon = engine.time
    leak_den = engine.leak_den
    assert leak_den > 0.0
    engine.temperature = 760.0
    engine.run(until=horizon)  # nothing to do, the draw is discarded lazily
    engine.run(until=first_event * 0.75)
    new_den = engine.leak_den
    assert new_den > leak_den
    assert engine.time == first_event * 0.75


class _RejectingSampler:
    """Smallest thinning construct: every proposal is rejected and consumes one draw."""

    cells = ()
    total_bound_propensity = 3.0
    bound_checks = 0
    bound_violations = 0

    def select_pair(self, rng):
        rng.random()
        self.bound_checks += 1
        return SimpleNamespace(accepted=False)


class _ThinningEngine(IsothermalSSA):
    """Real non-MET chemistry plus a MET leaf whose proposals are always rejected."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, incremental=False, **kwargs)
        self._fake = SimpleNamespace(
            channels={},
            component_channels={},
            sampler=_RejectingSampler(),
        )

    def _met_population(self):
        return self._fake


def _thinning_engine(seed):
    return _ThinningEngine(
        _state(),
        _records(),
        temperature=700.0,
        volume=1.0,
        rng=np.random.default_rng(seed),
    )


@pytest.mark.parametrize("seed", SEEDS)
def test_rejected_thinning_proposals_are_horizon_independent(seed):
    single = _thinning_engine(seed)
    single.run(until=HORIZON)
    assert single.bound_checks > 5, "fixture must exercise rejected proposals"
    assert single.events
    for horizons in (
        [HORIZON * i / 10 for i in range(1, 11)],
        sorted(float(x) for x in np.random.default_rng(seed).uniform(0, HORIZON, 96))
        + [HORIZON],
    ):
        piecewise = _thinning_engine(seed)
        _run_schedule(piecewise, horizons)
        assert _snapshot(piecewise) == _snapshot(single)


@pytest.mark.parametrize("seed", SEEDS)
def test_horizon_between_a_rejected_proposal_and_the_next_one(seed):
    # Probe where proposals fall by stepping a reference engine proposal by proposal.
    reference = _thinning_engine(seed)
    times = []
    original = reference._fake.sampler.select_pair

    def spy(rng):
        times.append(reference.time)
        return original(rng)

    reference._fake.sampler.select_pair = spy
    reference.run(until=HORIZON)
    assert len(times) > 5
    horizons = [float(np.nextafter(t, np.inf)) for t in times] + [HORIZON]
    piecewise = _thinning_engine(seed)
    _run_schedule(piecewise, horizons)
    assert _snapshot(piecewise) == _snapshot(reference)


def test_stop_reasons_and_global_event_counter():
    engine = _engine(11)
    result = engine.run(until=2.0)
    assert result.stop_reason == "horizon"
    assert result.total_events == len(engine.events)

    capped = _engine(11)
    first = capped.run(max_events=3)
    assert first.stop_reason == "max_events"
    assert len(first.events) == 3 and first.total_events == 3
    second = capped.run(max_events=3)  # per-call cap is kept
    assert second.stop_reason == "max_events"
    assert len(second.events) == 6 and second.total_events == 6
    assert capped.run(until=HORIZON).total_events == len(capped.events)


def test_absorbing_state_is_distinguishable_from_the_horizon():
    engine = _engine(11, reversible=False)
    result = engine.run(until=1000.0)
    assert result.stop_reason == "absorbing"
    assert len(result.events) == 3 and result.total_events == 3
    assert result.time == 1000.0  # the clock still advances to the horizon
    again = engine.run(until=2000.0)
    assert again.stop_reason == "absorbing" and again.time == 2000.0

    open_ended = _engine(11, reversible=False).run(max_events=10)
    assert open_ended.stop_reason == "absorbing"
    assert len(open_ended.events) == 3


def test_absorbing_runs_are_horizon_independent():
    single = _engine(11, reversible=False)
    single.run(until=1000.0)
    piecewise = _engine(11, reversible=False)
    _run_schedule(piecewise, [0.1, 0.5, 3.0, 40.0, 1000.0])
    assert _snapshot(piecewise) == _snapshot(single)
