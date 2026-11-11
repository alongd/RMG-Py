"""Red tests for the long-chain PS melt (I-125).

Each test asserts the behaviour a runnable DP 10-30 melt needs and is marked
``xfail(strict=True)`` with the claim it demonstrates; it turns green (and, being
strict, fails the run until the marker is removed) when the defect is fixed.
Evidence and numbers: ``long_chain_runtime_probe.py report`` and
/home/alon/runs/i125/REPORT.md.  The reduced artifact is built once per
module (about 3 minutes cold) into I125_ARTIFACT_DIR.
"""


import os

import numpy as np
import pytest

import long_chain_runtime_probe as probe
from rmgpy.kmc.ssa import IsothermalSSA, SiteIndex
from rmgpy.kmc.state import LedgerError


@pytest.fixture(scope="module")
def records():
    return probe.load_records()


def test_probe_harness_matches_on_a_dp3_chain(records):
    """Control (must pass): the same artifact and builder do match a DP3 chain."""
    state = probe.build_state([probe.ChainGraph(3, "t", ("end",)), probe.ChainGraph(3, "p")])
    index = SiteIndex(state, records)
    assert sum(bool(index.candidates(r.event_id)) for r in records) > 0


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="C1: DP3 capped reactant graphs embed in no position of a DP20 chain "
    "(cap hydrogens are not chain continuation), even where the chain-end "
    "environment is chemically valid",
)
def test_c1_records_match_chemically_valid_sites_of_a_dp20_chain(records):
    subset = [r for r in records
              if r.arity == 1 and r.family in ("R_Addition_MultipleBond", "intra_H_migration")]
    state = probe.build_state(probe.scenario_chains("interior+end"))
    index = SiteIndex(state, subset)
    graph = probe._nx_chain(state)
    expecting = [r for r in subset if probe.expected_sites(r, graph)[0]]
    assert expecting, "probe premise: the DP20 radical chain has chemically valid sites"
    missing = [r.event_id for r in expecting if not index.candidates(r.event_id)]
    assert not missing, f"{len(missing)} of {len(expecting)} records with valid sites match nowhere"


@pytest.mark.xfail(
    strict=True,
    raises=LedgerError,
    reason="C2: _split_strand's element ledger needs every strand atom mapped or an "
    "unmapped repeat-unit formula, so a cut inside a reaction window of a DP20 "
    "chain raises LedgerError",
)
def test_c2_backbone_cut_inside_a_reaction_window_splits_a_dp20_chain():
    chain = probe.ChainGraph(probe.DP, "t", ("interior", "end"))
    state = probe.build_state([chain], with_styrene=False, repeat_unit_formula={"C": 4, "H": 4})
    strand = state.strands["t"]
    by_position = {r.position: r.uuid for r in strand.atom_refs.values() if r.role == "backbone"}
    cut = len(chain.backbone) // 2 - 1
    first, second = by_position[cut], by_position[cut + 1]
    window = probe._window(state, "t", first, second)
    elements = {u: strand.atom_graph[u]["element"] for u in window}
    right_id, _ = state._split_strand(first, second, elements)
    for piece in (state.strands["t"], state.strands[right_id]):
        assert dict(piece.formula) == probe._graph_formula(piece)


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="C3: a closed DP10-30 melt has zero total propensity (absorbing at t=0) "
    "although the artifact holds radical-free-substrate homolysis records that "
    "fire on a closed DP3 melt",
)
def test_c3_closed_long_chain_melt_has_initiation_propensity(records):
    state = probe.build_state(probe.closed_melt_chains({10: 1, 20: 1, 30: 1}))
    engine = IsothermalSSA(state, records, temperature=700.0, volume=1e-24,
                           rng=np.random.default_rng(1),
                           inventory="R1")
    assert engine.channel_propensities().total_enabled > 0.0


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="C4: every compiled record has reverse_of set, so the coproduct debit is "
    "skipped; a released H atom stays in the polymer graph and is matched by later events",
)
def test_c4_released_coproduct_atoms_leave_the_active_graph(records):
    result = probe.c4_coproduct(records, {"H": 1}, clear_reverse_of=False)
    assert result["coproduct_atoms"] == 1
    assert result["sites_touching_coproduct_atoms"] == 0
    assert result["sink"] == {"H": 1}


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="C5: KMCState.moments is (count, sum nominal mass, sum mass^2), not the "
    "(count, sum DP, sum DP^2) of the declared realized ensemble",
)
def test_c5_state_moments_are_the_declared_dp_moments():
    result = probe.c5_report({10: 3, 20: 4, 30: 3})
    assert result["state_moments"] == result["declared_dp_moments"]


def test_c6_trajectory_is_independent_of_the_snapshot_schedule(records):
    single = probe.run_trajectory(records, "named", 1, "single")
    piecewise = probe.run_trajectory(records, "named", 1, "ten")
    assert [e[1] for e in single["events"]] == [e[1] for e in piecewise["events"]]


C7_BUILDS = 5
C7_MIN_PERMUTED_WIDTH = 4


def _require_permuted_width(orders, label):
    """Precondition of the C7 bound; fails as pytest.fail (not AssertionError), so never xfails."""
    widths = [probe.permuted_width(order) for order in orders]
    if len(widths) != C7_BUILDS or min(widths) < C7_MIN_PERMUTED_WIDTH:
        pytest.fail(f"C7 premise broken in {label}: permuted widths {widths}, need "
                    f"{C7_BUILDS} builds of width >= {C7_MIN_PERMUTED_WIDTH}; the 24^-5 bound no longer holds")


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="C7: atom ids come from uuid4 and candidates are ordered by uuid, so the "
    "same seed in two fresh processes indexes a different candidate order",
)
def test_c7_same_seed_in_two_fresh_processes_sees_the_same_candidate_order():
    """The list a seeded draw indexes into must not depend on the process.

    Each of two fresh processes builds the DP3 control melt 5 times through the uuid4
    path and records the t=0 candidate order of every record, ids left out.  The four
    closed chains c0..c3 each give one candidate (a site carrying the strand's 50 ids)
    to the same records; SiteIndex orders them by each strand's smallest id, the minimum
    of 50 independent uuid4s per strand, so their order is a uniform permutation of 4:
    two builds agree with probability 1/4! = 1/24 (id collisions, 2^-122, ignored).  Builds draw fresh ids, so
    all 5 pairs agree by chance with probability 24^-5 = 1.3e-7: the verdict is fixed.
    The width (>= 4) is checked first and fails loudly, outside the xfail, if the melt
    stops providing it.  The single trajectory pair this test used to compare agrees
    about half the time (test_c7_single_pair_verdict_is_a_coin_flip).
    """
    first = probe.candidate_orders_in_subprocess("uuid4", C7_BUILDS)
    second = probe.candidate_orders_in_subprocess("uuid4", C7_BUILDS)
    _require_permuted_width(first, "process 1")
    _require_permuted_width(second, "process 2")
    assert first == second


def test_c7_named_candidate_order_is_reproducible_across_processes():
    """Control (must pass): with deterministic ids the same order comes back in a fresh process."""
    first = probe.candidate_orders_in_subprocess("named", 2)
    second = probe.candidate_orders_in_subprocess("named", 2)
    assert first == second
    assert min(probe.permuted_width(order) for order in first) >= C7_MIN_PERMUTED_WIDTH


SLOW_RATE_CHECKS = "RMG_KMC_SLOW_RATE_CHECKS"
SINGLE_PAIR_RATE_PAIRS = 20


@pytest.mark.skipif(not os.environ.get(SLOW_RATE_CHECKS),
                    reason=f"slow (about 3 min): set {SLOW_RATE_CHECKS}=1 to run")
def test_c7_single_pair_verdict_is_a_coin_flip():
    """One pair of uuid4 trajectories decides the old C7 verdict, and that verdict is chance.

    Two fresh processes each run 20 seed-1 uuid4 trajectories (every one builds its own
    melt, so its ids are fresh draws); trajectory i of one is compared with trajectory i
    of the other, exactly as the single-pair C7 assertion compares its two.  Measured
    same-pick rate p = 40/81 (about 0.49), so both outcomes occur; this check misses one
    with probability p^20 + (1 - p)^20, about 2e-6 at p = 0.5 and 3.7e-5 at p = 0.6.
    """
    first = probe.trace_many_in_subprocess("uuid4", 1, SINGLE_PAIR_RATE_PAIRS)
    second = probe.trace_many_in_subprocess("uuid4", 1, SINGLE_PAIR_RATE_PAIRS)
    same = [a["applied"] == b["applied"] for a, b in zip(first, second)]
    assert len(same) == SINGLE_PAIR_RATE_PAIRS
    assert True in same and False in same, f"single-pair verdicts {same}"
