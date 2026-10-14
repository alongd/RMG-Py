"""Red tests for the long-chain PS melt (I-125).

Each test asserts the behaviour a runnable DP 10-30 melt needs and is marked
``xfail(strict=True)`` with the claim it demonstrates; it turns green (and, being
strict, fails the run until the marker is removed) when the defect is fixed.
Evidence and numbers: ``long_chain_runtime_probe.py report`` and
/home/alon/runs/i125/REPORT.md.  The reduced artifact is built once per
module (about 3 minutes cold) into I125_ARTIFACT_DIR.
"""


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


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="C6: run(until=T) draws a waiting time and discards it at the horizon, so "
    "ten snapshot calls consume a different random stream than one call",
)
def test_c6_trajectory_is_independent_of_the_snapshot_schedule(records):
    single = probe.run_trajectory(records, "named", 1, "single")
    piecewise = probe.run_trajectory(records, "named", 1, "ten")
    assert [e[1] for e in single["events"]] == [e[1] for e in piecewise["events"]]


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="C7: atom ids come from uuid4 and candidates are ordered by uuid, so the "
    "same seed in two fresh processes picks different sites",
)
def test_c7_same_seed_in_two_fresh_processes_picks_the_same_sites():
    first = probe.trace_in_subprocess("uuid4", 1)
    second = probe.trace_in_subprocess("uuid4", 1)
    assert first["applied"] == second["applied"]
