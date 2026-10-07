"""Generate and verify the immutable Step-1 legacy PS kMC baseline.

This module is test-only.  It observes commit-238930123 behavior without
changing production code or compiling chemistry.  The large event artifact is
read from its content-addressed phase-2b location.
"""

from __future__ import annotations

import argparse
import base64
import copy
import hashlib
import importlib.util
import json
import platform
import subprocess
import sys
from collections import Counter
from dataclasses import asdict
from pathlib import Path
from typing import Any, Mapping, Sequence

import networkx
import numpy as np
import scipy

import rmgpy
from rmgpy.kmc.event_record import EventRecord
from rmgpy.kmc.kinetics_library import matches_head_to_tail
from rmgpy.kmc.met import TRANSPORT_ARMS, compile_bulk_table, diffusion_rate
from rmgpy.kmc.ssa import (
    IsothermalSSA,
    SiteIndex,
    build_met_population,
    canonical_propensities,
    channel_propensities,
    eligible_candidates,
)
from rmgpy.kmc.state import AtomRef, KMCState, Site, Strand


ROOT = Path(__file__).resolve().parents[3]
DATABASE = Path("/home/alon/runs/phase2b/database")
ARTIFACT = Path(
    "/home/alon/runs/phase2b/cache/artifact/"
    "69c8023b80e9c0d58475650b8a1e773513bb958f6d5b69afca6aefee133d500c.json"
)
DESIGN = Path("/home/alon/runs/i058-generic-design/generic_linear_kmc_design.md")
FIXTURE = Path(__file__).with_name("fixtures") / "i052_legacy_baseline.json"
MANIFEST = Path(__file__).with_name("fixtures") / "i052_legacy_baseline_manifest.json"
RMG_PY_PIN = "23893012393d946efc8d70529bfc58780b8cdf86"
STEP1_FREEZE_COMMIT = "d1459e1b9cd1f9d877ac820b0152ee931673b7c7"
RMG_DATABASE_PIN = "cd86d4e1c187a132109e16cd86f624ed9fb217df"
MOM_PIN = "c220aa83c3e6f3016c9998e5efa69a10c4de6bb5"
STEP1_FREEZE_ONLY_FILES = frozenset({
    "test/rmgpy/kmc/fixtures/i052_legacy_baseline.json",
    "test/rmgpy/kmc/i052BaselineFreezeTest.py",
    "test/rmgpy/kmc/i052_baseline_freeze.py",
})
TEMPERATURE = 700.0
VOLUME = 1.0e-24


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _git(root: Path, *arguments: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(root), *arguments], text=True
    ).strip()


def _require_git_commit(root: Path, commit: str) -> None:
    """Fail loudly when git or an immutable baseline commit is unavailable."""
    subprocess.check_call(
        ["git", "-C", str(root), "cat-file", "-e", f"{commit}^{{commit}}"]
    )


def _git_blob(root: Path, commit: str, relative: str) -> bytes:
    """Read an exact committed blob; missing git, commit, or path is fatal."""
    return subprocess.check_output(
        ["git", "-C", str(root), "cat-file", "blob", f"{commit}:{relative}"]
    )


def _verify_frozen_json_bytes() -> None:
    """The record files themselves remain exactly as committed by Step 1."""
    for path in (MANIFEST, FIXTURE):
        relative = str(path.relative_to(ROOT))
        expected = _git_blob(ROOT, STEP1_FREEZE_COMMIT, relative)
        assert path.read_bytes() == expected, relative


def _state_oracle():
    path = Path(__file__).with_name("stateTest.py")
    spec = importlib.util.spec_from_file_location("i052_state_oracle", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def canonical_state_bytes(state: KMCState) -> bytes:
    """Return the exact bytes fed to legacy ``KMCState.state_hash``."""
    payload = {
        "strands": {
            key: {
                "length": strand.length,
                "features": strand.features,
                "formula": strand.formula,
                "radicals": strand.radical_count,
                "atoms": {
                    index: ref.__dict__
                    for index, ref in sorted(strand.atom_refs.items())
                },
                "atom_graph": strand.atom_graph,
                "repeat_unit_formula": strand.repeat_unit_formula,
            }
            for key, strand in sorted(state.strands.items())
        },
        "junctions": sorted(state.junctions),
        "sink": state.sink,
        "sink_radicals": state.sink_radicals,
    }
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), default=str
    ).encode()
    assert hashlib.sha256(encoded).hexdigest() == state.state_hash()
    return encoded


def state_snapshot(state: KMCState) -> dict[str, Any]:
    encoded = canonical_state_bytes(state)
    return {
        "canonical_bytes_b64": base64.b64encode(encoded).decode("ascii"),
        "state_hash": state.state_hash(),
        "components": dict(sorted(state.components.items())),
        "component_lengths": dict(sorted(state.component_lengths.items())),
        "component_ledgers": dict(sorted(state.component_ledgers.items())),
        "component_radicals": dict(sorted(state.component_radicals.items())),
        "moments": list(state.moments),
        "total_ledger": state.total_ledger,
        "total_mass": state.total_mass,
        "total_radicals": state.total_radicals,
        "sink": state.sink,
        "sink_radicals": state.sink_radicals,
        "junctions": [list(pair) for pair in sorted(state.junctions)],
        "event_ids": [entry.event_id for entry in state.event_log],
    }


def _site_snapshot(site: Site) -> dict[str, Any]:
    refs = sorted(
        node["atom_ref"].uuid
        for node in site.graph.values()
        if isinstance(node.get("atom_ref"), AtomRef)
    )
    root = site.graph[site.atom]["atom_ref"]
    return {
        "site_type": site.site_type,
        "strand_id": site.strand_id,
        "position": site.position,
        "atom": site.atom,
        "root_uuid": root.uuid,
        "binding_uuids": refs,
    }


def candidate_snapshot(state: KMCState, records: Sequence[Any]) -> dict[str, Any]:
    index = SiteIndex(state, records)
    rows = {}
    for record in records:
        event_id = str(record.event_id)
        candidates = eligible_candidates(state, index, record)
        rows[event_id] = {
            "record": {
                "event_id": event_id,
                "status": record.status,
                "rate_order": record.rate_order,
                "rate_units": record.rate_units,
                "degeneracy": record.degeneracy,
                "raw_path_degeneracy": record.raw_path_degeneracy,
                "ssa_multiplier": record.ssa_multiplier,
                "reactant_pair_convention": record.reactant_pair_convention,
            },
            "candidates": [
                [_site_snapshot(site) for site in candidate] for candidate in candidates
            ],
        }
    return rows


def propensity_snapshot(state: KMCState, records: Sequence[Any]) -> dict[str, Any]:
    index = SiteIndex(state, records)
    report = channel_propensities(state, index, records, TEMPERATURE, VOLUME)
    return {
        "enabled_channel_ids": sorted(report.enabled_channel_ids),
        "propensities": dict(report.propensities),
        "refused_propensities": dict(report.refused_propensities),
        "total_enabled": report.total_enabled,
        "total_refused": report.total_refused,
    }


def _entry_snapshot(entry: Any) -> dict[str, Any]:
    return {
        "event_id": entry.event_id,
        "participants": list(entry.participants),
        "bindings": dict(entry.bindings),
        "derived_cut_offsets": list(entry.derived_cut_offsets),
        "inheritance": dict(entry.inheritance),
        "formed_bond": list(entry.formed_bond) if entry.formed_bond else None,
        "product_graphs": entry.product_graphs,
    }


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return [_jsonable(item) for item in sorted(value, key=repr)]
    if isinstance(value, np.generic):
        return value.item()
    return value


class TracingRNG:
    """Record exactly the random values requested by the legacy runtime."""

    def __init__(self, seed: int):
        self.generator = np.random.default_rng(seed)
        self.trace: list[dict[str, Any]] = []

    @staticmethod
    def _plain(value: Any) -> Any:
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        return value

    def _call(self, method: str, *args: Any, **kwargs: Any) -> Any:
        value = getattr(self.generator, method)(*args, **kwargs)
        self.trace.append(
            {
                "method": method,
                "args": [self._plain(item) for item in args],
                "kwargs": {key: self._plain(item) for key, item in kwargs.items()},
                "value": self._plain(value),
            }
        )
        return value

    def random(self, *args: Any, **kwargs: Any) -> Any:
        return self._call("random", *args, **kwargs)

    def exponential(self, *args: Any, **kwargs: Any) -> Any:
        return self._call("exponential", *args, **kwargs)

    def integers(self, *args: Any, **kwargs: Any) -> Any:
        return self._call("integers", *args, **kwargs)

    def choice(self, *args: Any, **kwargs: Any) -> Any:
        return self._call("choice", *args, **kwargs)


def _mapped_fixtures(artifact: Mapping[str, Any], records: dict[str, Any]):
    by_id = {item["event_id"]: item for item in artifact["records"]}
    unary_data = by_id[
        "evt_0040e2cc0127d4797d6b0c6d37f033428667c836482ebc1f9fc85b4c1392d65c"
    ]
    record = EventRecord.from_dict(unary_data)
    records[record.event_id] = record.to_dict()
    oracle = _state_oracle()
    fixtures = []

    def mapped_state(namespaces, lengths, *, doubly_featured=False):
        strands = []
        for namespace, length in zip(namespaces, lengths):
            harness, _, _ = oracle._harness(unary_data, namespace=namespace)
            strand = copy.deepcopy(next(iter(harness.strands.values())))
            strand.length = length
            if namespace == "interior-":
                strand.atom_refs = {
                    key: AtomRef(ref.uuid, ref.strand_id, ref.position + 7, "interior")
                    for key, ref in strand.atom_refs.items()
                }
            if doubly_featured and namespace == namespaces[-1]:
                strand.features = {
                    1: {"legacy_observation_tag": "alpha"},
                    length - 2: {"legacy_observation_tag": "beta"},
                }
            strands.append(strand)
        return KMCState(strands)

    states = {
        "ends_interiors_doubly_featured": mapped_state(
            ("end-", "interior-", "double-"),
            (52, 67, 83),
            doubly_featured=True,
        ),
        "same_local_class_different_dp": mapped_state(
            ("dp-short-", "dp-long-"), (52, 91)
        ),
    }
    for name, state in states.items():
        fixtures.append(
            {
                "category": "mapped_state",
                "name": name,
                "record_ids": [record.event_id],
                "outcome": {"status": "observed", "exception": None},
                "pre_state": state_snapshot(state),
                "post_state": state_snapshot(state),
                "candidates": candidate_snapshot(state, (record,)),
                "propensity_vector": propensity_snapshot(state, (record,)),
            }
        )

    equivalent_a, _, _ = oracle._harness(unary_data, namespace="equivalent-a-")
    equivalent_b, _, _ = oracle._harness(unary_data, namespace="equivalent-b-")
    other, _, _ = oracle._harness(unary_data, namespace="inequivalent-")
    root_sources = (
        next(iter(equivalent_a.strands.values())),
        next(iter(equivalent_b.strands.values())),
        next(iter(other.strands.values())),
    )
    root_refs = {}
    root_graph = {}
    root_formula = Counter()
    root_radicals = 0
    for ordinal, source in enumerate(root_sources):
        position_shift = 70 if ordinal == 2 else 0
        for ref in source.atom_refs.values():
            root_refs[len(root_refs)] = AtomRef(
                ref.uuid,
                "root-chain",
                ref.position + position_shift,
                "interior" if ordinal == 2 else ref.role,
            )
        root_graph.update(copy.deepcopy(source.atom_graph))
        root_formula.update(source.formula)
        root_radicals += source.radical_count
    root_state = KMCState(
        (
            Strand(
                "root-chain",
                140,
                formula=dict(root_formula),
                radical_count=root_radicals,
                atom_refs=root_refs,
                atom_graph=root_graph,
            ),
        )
    )
    fixtures.append(
        {
            "category": "mapped_state",
            "name": "two_equivalent_plus_one_inequivalent_root",
            "record_ids": [record.event_id],
            "outcome": {"status": "observed", "exception": None},
            "pre_state": state_snapshot(root_state),
            "post_state": state_snapshot(root_state),
            "candidates": candidate_snapshot(root_state, (record,)),
            "propensity_vector": propensity_snapshot(root_state, (record,)),
            "root_classes": {
                "equivalent": [
                    "equivalent-a-" + record.event_id + "-0-0",
                    "equivalent-b-" + record.event_id + "-0-0",
                ],
                "inequivalent": ["inequivalent-" + record.event_id + "-0-0"],
            },
        }
    )

    identical_data = next(
        item
        for item in artifact["records"]
        if item.get("reactant_pair_convention") == "unordered-identical-pair N(N-1)/2"
    )
    pair_state, pair_record, _ = oracle._harness(identical_data)
    records[pair_record.event_id] = pair_record.to_dict()
    fixtures.append(
        {
            "category": "mapped_state",
            "name": "full_species_identical_pair",
            "record_ids": [pair_record.event_id],
            "outcome": {"status": "observed", "exception": None},
            "pre_state": state_snapshot(pair_state),
            "post_state": state_snapshot(pair_state),
            "candidates": candidate_snapshot(pair_state, (pair_record,)),
            "propensity_vector": propensity_snapshot(pair_state, (pair_record,)),
        }
    )
    return fixtures


def _attachment_records(artifact: Mapping[str, Any], kind: str):
    by_id = {item["event_id"]: item for item in artifact["records"]}
    capture_data = next(
        item
        for item in artifact["records"]
        if item.get("junction_ops")
        and item["junction_ops"][0].get("junction_kind") == kind
        and item["junction_ops"][0].get("action") == "create"
    )
    return capture_data, by_id[capture_data["reverse_of"]]


def _capture_harness(oracle: Any, capture_data: Mapping[str, Any], namespace=""):
    form_atoms = next(
        operation["atoms"]
        for operation in capture_data["bond_ops"]
        if operation["action"] == "form"
    )
    junction = capture_data["junction_ops"][0]
    labels = junction["formed_bond"]["label_pair"]
    return oracle._harness(
        capture_data,
        roles={form_atoms[0]: "attacked_ring", form_atoms[1]: "attacker"},
        labels={form_atoms[0]: labels[0], form_atoms[1]: labels[1]},
        namespace=namespace,
    )


def _reverse_site_for_entry(state: KMCState, reverse: EventRecord, entry: Any) -> Site:
    index = SiteIndex(state, (reverse,))
    expected = set(entry.bindings.values())
    for candidate in eligible_candidates(state, index, reverse):
        supplied = {
            node["atom_ref"].uuid
            for site in candidate
            for node in site.graph.values()
            if isinstance(node.get("atom_ref"), AtomRef)
        }
        if expected <= supplied:
            assert len(candidate) == 1
            return candidate[0]
    raise AssertionError("legacy index did not expose the bound reverse candidate")


def _attachment_fixtures(artifact: Mapping[str, Any], records: dict[str, Any]):
    oracle = _state_oracle()
    by_id = {item["event_id"]: item for item in artifact["records"]}
    downstream_data = by_id[
        "evt_0040e2cc0127d4797d6b0c6d37f033428667c836482ebc1f9fc85b4c1392d65c"
    ]
    downstream_record = EventRecord.from_dict(downstream_data)
    records[downstream_record.event_id] = downstream_record.to_dict()
    fixtures = []
    for kind in ("J_para", "J_ortho_S6", "J_ortho_S7"):
        capture_data, reverse_data = _attachment_records(artifact, kind)
        capture = EventRecord.from_dict(capture_data)
        reverse = EventRecord.from_dict(reverse_data)
        records[capture.event_id] = capture.to_dict()
        records[reverse.event_id] = reverse.to_dict()

        state, _, sites = _capture_harness(oracle, capture_data)
        before = state_snapshot(state)
        fixtures.append(
            {
                "category": "attachment",
                "name": f"{kind}_before_association",
                "variant": kind,
                "record_ids": [capture.event_id, reverse.event_id],
                "outcome": {"status": "not_fired"},
                "pre_state": before,
                "post_state": before,
                "candidates": candidate_snapshot(state, (capture, reverse)),
                "propensity_vector": propensity_snapshot(state, (capture, reverse)),
            }
        )

        state, _, sites = _capture_harness(oracle, capture_data)
        pre = state_snapshot(state)
        entry = state.apply(capture, sites)
        fixtures.append(
            {
                "category": "attachment",
                "name": f"{kind}_after_association",
                "variant": kind,
                "record_ids": [capture.event_id, reverse.event_id],
                "outcome": {"status": "success", "event": _entry_snapshot(entry)},
                "pre_state": pre,
                "post_state": state_snapshot(state),
                "candidates": candidate_snapshot(state, (capture, reverse)),
                "propensity_vector": propensity_snapshot(state, (capture, reverse)),
            }
        )

        state, _, sites = _capture_harness(oracle, capture_data)
        pre = state_snapshot(state)
        capture_entry = state.apply(capture, sites)
        reverse_site = _reverse_site_for_entry(state, reverse, capture_entry)
        reverse_entry = state.apply(reverse, [reverse_site])
        fixtures.append(
            {
                "category": "attachment",
                "name": f"{kind}_immediate_dissociation",
                "variant": kind,
                "record_ids": [capture.event_id, reverse.event_id],
                "outcome": {
                    "status": "success",
                    "events": [
                        _entry_snapshot(capture_entry),
                        _entry_snapshot(reverse_entry),
                    ],
                },
                "pre_state": pre,
                "post_state": state_snapshot(state),
                "candidates": candidate_snapshot(state, (capture, reverse)),
                "propensity_vector": propensity_snapshot(state, (capture, reverse)),
            }
        )

        state, _, sites = _capture_harness(oracle, capture_data)
        capture_entry = state.apply(capture, sites)
        reverse_site = _reverse_site_for_entry(state, reverse, capture_entry)
        bound_positions = {
            (state._ref(atom_uuid).strand_id, state._ref(atom_uuid).position)
            for atom_uuid in capture_entry.bindings.values()
        }
        distant_strand = next(
            strand
            for strand in state.strands.values()
            if any(
                (strand.strand_id, position) not in bound_positions
                for position in range(strand.length)
            )
        )
        distant_position = next(
            position
            for position in range(distant_strand.length)
            if (distant_strand.strand_id, position) not in bound_positions
        )
        distant_strand.features[distant_position] = {
            "legacy_distant_component_edit": kind
        }
        distant_strand.features = dict(sorted(distant_strand.features.items()))
        pre_failure = state_snapshot(state)
        try:
            state.apply(reverse, [reverse_site])
        except Exception as error:  # legacy exception type/message are frozen
            failure = {"type": type(error).__name__, "message": str(error)}
        else:  # pragma: no cover - an unexpected legacy change
            failure = {"type": None, "message": None}
        fixtures.append(
            {
                "category": "attachment",
                "name": f"{kind}_distant_edit_then_dissociation",
                "variant": kind,
                "record_ids": [capture.event_id, reverse.event_id],
                "outcome": {
                    "status": "failure",
                    "exception": failure,
                    "edit": {
                        "kind": "direct_legacy_component_feature_edit",
                        "strand_id": distant_strand.strand_id,
                        "position": distant_position,
                    },
                },
                "pre_state": pre_failure,
                "post_state": state_snapshot(state),
                "candidates": candidate_snapshot(state, (capture, reverse)),
                "propensity_vector": propensity_snapshot(state, (capture, reverse)),
            }
        )

        state, _, sites = _capture_harness(oracle, capture_data)
        capture_entry = state.apply(capture, sites)
        reverse_site = _reverse_site_for_entry(state, reverse, capture_entry)
        endpoint = capture_entry.formed_bond[0]
        endpoint_ref = state._ref(endpoint)
        endpoint_node = state.strands[endpoint_ref.strand_id].atom_graph[endpoint]
        old_radical = int(endpoint_node.get("radical", 0))
        new_radical = 1 if old_radical == 0 else 0
        state._feature_atom(endpoint, radical=new_radical)
        state.strands[endpoint_ref.strand_id].radical_count += new_radical - old_radical
        pre_failure = state_snapshot(state)
        try:
            state.apply(reverse, [reverse_site])
        except Exception as error:
            failure = {"type": type(error).__name__, "message": str(error)}
        else:  # pragma: no cover
            failure = {"type": None, "message": None}
        fixtures.append(
            {
                "category": "attachment",
                "name": f"{kind}_conflicting_edit",
                "variant": kind,
                "record_ids": [capture.event_id, reverse.event_id],
                "outcome": {
                    "status": "failure",
                    "exception": failure,
                    "edit": {
                        "kind": "direct_legacy_binding_edit",
                        "atom_uuid": endpoint,
                        "field": "radical",
                        "before": old_radical,
                        "after": new_radical,
                    },
                },
                "pre_state": pre_failure,
                "post_state": state_snapshot(state),
                "candidates": candidate_snapshot(state, (capture, reverse)),
                "propensity_vector": propensity_snapshot(state, (capture, reverse)),
            }
        )

        left, _, left_sites = _capture_harness(oracle, capture_data, "left-")
        right, _, right_sites = _capture_harness(oracle, capture_data, "right-")
        downstream, _, downstream_sites = oracle._harness(
            downstream_data, namespace="downstream-"
        )
        combined = KMCState(
            tuple(left.strands.values())
            + tuple(right.strands.values())
            + tuple(downstream.strands.values())
        )
        left_owner = next(iter(left.strands.values()))
        downstream_owner = next(iter(downstream.strands.values()))
        ownership_edge = tuple(
            sorted(
                (
                    sorted(left_owner.atom_graph)[0],
                    sorted(downstream_owner.atom_graph)[0],
                )
            )
        )
        combined.junctions.add(ownership_edge)
        combined._refresh_components()
        pre = state_snapshot(combined)
        first = combined.apply(capture, left_sites)
        assert (
            combined.components[downstream_sites[0].strand_id]
            == combined.components[left_sites[0].strand_id]
        )
        chemistry = combined.apply(downstream_record, downstream_sites)
        second = combined.apply(capture, right_sites)
        fixtures.append(
            {
                "category": "attachment",
                "name": f"{kind}_downstream_further_attachment",
                "variant": kind,
                "record_ids": [
                    capture.event_id,
                    downstream_record.event_id,
                    reverse.event_id,
                ],
                "outcome": {
                    "status": "success",
                    "events": [
                        _entry_snapshot(first),
                        _entry_snapshot(chemistry),
                        _entry_snapshot(second),
                    ],
                    "downstream_component_edge": list(ownership_edge),
                    "open_capture_count": len(
                        combined._open_captures[capture.event_id]
                    ),
                },
                "pre_state": pre,
                "post_state": state_snapshot(combined),
                "candidates": candidate_snapshot(combined, (capture, reverse)),
                "propensity_vector": propensity_snapshot(combined, (capture, reverse)),
            }
        )
    return fixtures


def _plp_fixtures(artifact: Mapping[str, Any], records: dict[str, Any]):
    library_records = [
        item
        for item in artifact["records"]
        if item.get("rate_source", {}).get("kind") == "kMC kinetics library"
    ]
    intact = library_records[0]
    records[intact["event_id"]] = intact
    copied = copy.deepcopy(intact)
    oracle = _state_oracle()
    first_state, intact_record, _ = oracle._harness(intact, namespace="ancestry-a-")
    second_state, _, _ = oracle._harness(intact, namespace="ancestry-b-")
    defective_strands = [
        copy.deepcopy(strand) for strand in first_state.strands.values()
    ]
    defective_strand = defective_strands[0]
    defect_position = defective_strand.length - 1
    central_uuid = "defective-branch-central-C"
    branch_graph = {
        central_uuid: {
            "element": "C",
            "radical": 0,
            "charge": 0,
            "lone_pairs": 0,
            "implicit_hydrogens": 0,
            "edges": {},
        }
    }
    central_h = "defective-branch-central-H"
    branch_graph[central_h] = {
        "element": "H",
        "radical": 0,
        "charge": 0,
        "lone_pairs": 0,
        "implicit_hydrogens": 0,
        "edges": {central_uuid: "1.0"},
    }
    branch_graph[central_uuid]["edges"][central_h] = "1.0"
    for arm in range(3):
        carbon_uuid = f"defective-branch-arm-{arm}-C"
        branch_graph[carbon_uuid] = {
            "element": "C",
            "radical": 0,
            "charge": 0,
            "lone_pairs": 0,
            "implicit_hydrogens": 0,
            "edges": {central_uuid: "1.0"},
        }
        branch_graph[central_uuid]["edges"][carbon_uuid] = "1.0"
        for hydrogen in range(3):
            hydrogen_uuid = f"defective-branch-arm-{arm}-H-{hydrogen}"
            branch_graph[hydrogen_uuid] = {
                "element": "H",
                "radical": 0,
                "charge": 0,
                "lone_pairs": 0,
                "implicit_hydrogens": 0,
                "edges": {carbon_uuid: "1.0"},
            }
            branch_graph[carbon_uuid]["edges"][hydrogen_uuid] = "1.0"
    assert all(
        sum(float(order) for order in node["edges"].values())
        == (4.0 if node["element"] == "C" else 1.0)
        for node in branch_graph.values()
    )
    assert (
        sum(
            branch_graph[neighbor]["element"] == "C"
            for neighbor in branch_graph[central_uuid]["edges"]
        )
        == 3
    )
    for atom_uuid in sorted(branch_graph):
        defective_strand.atom_refs[max(defective_strand.atom_refs) + 1] = AtomRef(
            atom_uuid,
            defective_strand.strand_id,
            defect_position,
            "branch_defect",
        )
    defective_strand.atom_graph.update(branch_graph)
    defective_strand.formula["C"] = defective_strand.formula.get("C", 0) + 4
    defective_strand.formula["H"] = defective_strand.formula.get("H", 0) + 10
    defective_strand.features[defect_position] = {
        "legacy_global_defect": "branched_chain",
        "branch_points": 1,
        "topology": "branched",
    }
    defective_state = KMCState(defective_strands)
    first_candidates = candidate_snapshot(first_state, (intact_record,))
    second_candidates = candidate_snapshot(second_state, (intact_record,))
    defective_candidates = candidate_snapshot(defective_state, (intact_record,))
    first_propensity = propensity_snapshot(first_state, (intact_record,))
    second_propensity = propensity_snapshot(second_state, (intact_record,))
    defective_propensity = propensity_snapshot(defective_state, (intact_record,))
    ancestry_fields = sorted(
        key for key in intact if "ancestry" in key.lower() or "origin" in key.lower()
    )
    return [
        {
            "category": "plp_scope",
            "name": "plp_intact_chain",
            "record_ids": [intact["event_id"]],
            "outcome": {"status": "matched", "exception": None},
            "pre_state": state_snapshot(first_state),
            "post_state": state_snapshot(first_state),
            "candidates": first_candidates,
            "propensity_vector": first_propensity,
            "matcher": matches_head_to_tail(intact),
            "rate_source": intact["rate_source"],
            "k_table": intact["k_table"],
        },
        {
            "category": "plp_scope",
            "name": "plp_locally_matching_globally_defective_chain",
            "record_ids": [intact["event_id"]],
            "outcome": {"status": "matched", "exception": None},
            "pre_state": state_snapshot(defective_state),
            "post_state": state_snapshot(defective_state),
            "candidates": defective_candidates,
            "propensity_vector": defective_propensity,
            "matcher": matches_head_to_tail(intact),
            "global_defect": {
                "kind": "branch",
                "feature_position": defect_position,
                "central_atom_uuid": central_uuid,
                "formula": {"C": 4, "H": 10},
                "topology": "branched",
            },
        },
        {
            "category": "plp_scope",
            "name": "plp_isomorphic_different_ancestry",
            "record_ids": [intact["event_id"]],
            "outcome": {"status": "matched", "exception": None},
            "pre_state": state_snapshot(first_state),
            "post_state": state_snapshot(second_state),
            "candidates": {
                "first": first_candidates,
                "second": second_candidates,
            },
            "propensity_vector": {
                "first": first_propensity,
                "second": second_propensity,
            },
            "first_matcher": matches_head_to_tail(intact),
            "second_matcher": matches_head_to_tail(copied),
            "record_bytes_equal": json.dumps(intact, sort_keys=True)
            == json.dumps(copied, sort_keys=True),
            "legacy_ancestry_fields": ancestry_fields,
            "candidate_counts": [
                len(first_candidates[intact_record.event_id]["candidates"]),
                len(second_candidates[intact_record.event_id]["candidates"]),
            ],
            "propensity_totals": [
                first_propensity["total_enabled"],
                second_propensity["total_enabled"],
            ],
        },
    ]


def _met_fixture(artifact: Mapping[str, Any], records: dict[str, Any]):
    by_id = {item["event_id"]: item for item in artifact["records"]}
    met_data = (
        by_id["evt_01509d2c793aee254f85a53b02979318b715df2f84124107442a4aa0b04f69b3"],
        by_id["evt_0770369a122c8f0e948f8a2bf3f0098cecdae5c3bf6c234905c116ba88bec71b"],
    )
    met_records = tuple(EventRecord.from_dict(item) for item in met_data)
    records.update({record.event_id: record.to_dict() for record in met_records})
    oracle = _state_oracle()

    def met_state():
        strands = []
        for ordinal, data in enumerate(met_data):
            harness, _, _ = oracle._harness(data, namespace=f"met-{ordinal}-")
            for strand_index, strand in enumerate(harness.strands.values()):
                strand = copy.deepcopy(strand)
                strand.length += 11 * (ordinal + strand_index + 1)
                if ordinal == 1 and strand_index == 0:
                    root_key = min(strand.atom_refs)
                    root = strand.atom_refs[root_key]
                    strand.atom_refs[root_key] = AtomRef(
                        root.uuid,
                        root.strand_id,
                        strand.length // 2,
                        "interior",
                    )
                strands.append(strand)
        result = KMCState(strands)
        first_group = (strands[0], strands[2])
        second_group = (strands[1], strands[3])
        for group in (first_group, second_group):
            endpoints = tuple(sorted(strand.atom_graph)[0] for strand in group)
            result.junctions.add(tuple(sorted(endpoints)))
        result._refresh_components()
        return result

    state = met_state()
    index = SiteIndex(state, met_records)
    table = compile_bulk_table(met_records, TRANSPORT_ARMS["A0_REF_H_CROSS_NEc"], "R0")
    population = build_met_population(
        state,
        index,
        met_records,
        table,
        temperature=TEMPERATURE,
        volume=VOLUME,
        spin_factor=1.0,
    )
    proposal_rng = TracingRNG(52_001)
    proposals = []
    for _ in range(4):
        proposal = population.sampler.propose(proposal_rng)
        proposals.append(None if proposal is None else asdict(proposal))
    report = channel_propensities(
        state,
        index,
        met_records,
        TEMPERATURE,
        VOLUME,
        met_population=population,
    )

    clock_state = met_state()
    clock_rng = TracingRNG(52_002)
    engine = IsothermalSSA(
        clock_state,
        met_records,
        temperature=TEMPERATURE,
        volume=VOLUME,
        rng=clock_rng,
        termination_table=table,
        spin_factor=1.0,
        incremental=False,
    )
    event = engine.step()
    assert event is not None
    result = engine.run(max_events=0)

    def options(values):
        return [
            {
                "channel": _jsonable(asdict(option.channel)),
                "record_ids": [record.event_id for record in option.records],
            }
            for option in values
        ]

    return {
        "category": "met_rng_clock",
        "name": "mixed_site_types_classes_and_thinning",
        "record_ids": [record.event_id for record in met_records],
        "outcome": {
            "status": "success",
            "exception": None,
            "record": asdict(event),
        },
        "pre_state": state_snapshot(state),
        "candidates": candidate_snapshot(state, met_records),
        "met": {
            "channels": {
                repr(key): options(value)
                for key, value in sorted(
                    population.channels.items(), key=lambda row: repr(row[0])
                )
            },
            "component_channels": {
                repr(key): options(value)
                for key, value in sorted(
                    population.component_channels.items(), key=lambda row: repr(row[0])
                )
            },
            "item_sites": {
                key: [_site_snapshot(site) for site in value]
                for key, value in sorted(population.item_sites.items())
            },
            "cells": [asdict(cell) for cell in population.sampler.cells],
            "exact_cell_propensities": population.sampler.exact_cell_propensities(),
            "exact_total_propensity": population.sampler.exact_total_propensity,
            "total_bound_propensity": population.sampler.total_bound_propensity,
        },
        "thinning": {"proposals": proposals, "rng_variates": proposal_rng.trace},
        "propensity_vector": {
            "enabled_channel_ids": sorted(report.enabled_channel_ids),
            "propensities": dict(report.propensities),
            "refused_propensities": dict(report.refused_propensities),
            "total_enabled": report.total_enabled,
            "total_refused": report.total_refused,
        },
        "clock": {
            "rng_variates": clock_rng.trace,
            "event_time": event.time,
            "event": asdict(event),
        },
        "observables": asdict(result),
        "post_state": state_snapshot(clock_state),
    }


def _spin_fixture():
    arm = TRANSPORT_ARMS["A0_REF_H_CROSS_NEc"]
    rows = []
    for spin_factor in (1.0, 0.25):
        for pair_class in ("end/end", "end/mid", "mid/mid"):
            rows.append(
                {
                    "spin_factor": spin_factor,
                    "pair_class": pair_class,
                    "diffusion_rate": diffusion_rate(
                        arm,
                        TEMPERATURE,
                        7,
                        19,
                        pair_class,
                        spin_factor=spin_factor,
                    ),
                }
            )
    return {"category": "spin", "name": "met_spin_factor", "rates": rows}


def _ghost_fixture(records: dict[str, Any]):
    oracle = _state_oracle()
    state, record, sites, captured = oracle._captured_coproduct_radical_case()
    records[record.event_id] = record.to_dict()
    pre = state_snapshot(state)
    candidates = candidate_snapshot(state, (record,))
    propensity = propensity_snapshot(state, (record,))
    entry = state.apply(record, sites)
    active_atoms = sorted(
        atom_uuid
        for strand in state.strands.values()
        for atom_uuid in strand.atom_graph
    )
    return {
        "category": "ghost_coproduct",
        "name": "archived_coproduct_atoms_remain_active",
        "record_ids": [record.event_id],
        "captured_record_artifact_sha256": captured["record_artifact_sha256"],
        "outcome": {"status": "success", "event": _entry_snapshot(entry)},
        "pre_state": pre,
        "post_state": state_snapshot(state),
        "candidates": candidates,
        "propensity_vector": propensity,
        "active_atom_uuids": active_atoms,
    }


def _mass_fixture():
    formulas = {
        "ps_repeat": {"C": 8, "H": 8},
        "unknown_element": {"Xe": 1},
        "isotope_key": {"C13": 1},
        "combined": {"C": 8, "H": 8, "Xe": 1, "C13": 1},
    }
    rows = {}
    for name, formula in formulas.items():
        state = KMCState((Strand(name, 1, formula=formula),))
        rows[name] = {
            "formula": formula,
            "strand_mass": state.strands[name].mass,
            "total_mass": state.total_mass,
            "moments": list(state.moments),
        }
    return {
        "category": "nominal_mass",
        "name": "nominal_and_unknown_mass",
        "rows": rows,
    }


def _rollback_fixture(artifact: Mapping[str, Any], records: dict[str, Any]):
    data = next(item for item in artifact["records"] if len(item["bond_ops"]) > 2)
    oracle = _state_oracle()
    state, record, sites = oracle._harness(data)
    records[record.event_id] = record.to_dict()
    pre = state_snapshot(state)
    candidates = candidate_snapshot(state, (record,))
    propensity = propensity_snapshot(state, (record,))
    calls = 0

    def fail_after_first(operation):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise RuntimeError("injected mid-rewrite failure")

    state._after_operation = fail_after_first
    try:
        state.apply(record, sites)
    except Exception as error:
        failure = {"type": type(error).__name__, "message": str(error)}
    else:  # pragma: no cover
        failure = {"type": None, "message": None}
    return {
        "category": "atomic_rollback",
        "name": "mid_rewrite_failure_is_atomic",
        "record_ids": [record.event_id],
        "outcome": {
            "status": "failure",
            "exception": failure,
            "operation_calls": calls,
        },
        "pre_state": pre,
        "post_state": state_snapshot(state),
        "candidates": candidates,
        "propensity_vector": propensity,
    }


def _legacy_clock(vector: Mapping[str, Any], seed: int) -> dict[str, Any]:
    if vector.get("availability") == "not_applicable":
        return {
            "availability": "not_applicable",
            "rng_variates": [],
            "event_time": None,
            "event": None,
        }
    propensities = vector["propensities"]
    canonical = canonical_propensities(propensities.items())
    if canonical.absorbing:
        return {
            "availability": "absorbing",
            "propensity_vector": copy.deepcopy(vector),
            "rng_variates": [],
            "event_time": None,
            "event": None,
        }
    rng = TracingRNG(seed)
    event_time = float(rng.exponential(1.0 / canonical.total))
    event_id = canonical.draw(rng)
    return {
        "availability": "legacy_clock_draw",
        "propensity_vector": copy.deepcopy(vector),
        "rng_variates": rng.trace,
        "event_time": event_time,
        "event": {"event_id": event_id, "status": "selected"},
    }


def _clock_surface(name: str, vector: Mapping[str, Any]) -> dict[str, Any]:
    seed = int.from_bytes(hashlib.sha256(name.encode()).digest()[:8], "big")
    if "propensities" in vector:
        return _legacy_clock(vector, seed)
    states = {
        key: _legacy_clock(value, seed + ordinal)
        for ordinal, (key, value) in enumerate(sorted(vector.items()))
    }
    return {"availability": "legacy_multi_state_clock", "states": states}


def _complete_legacy_surface(fixture: dict[str, Any]) -> dict[str, Any]:
    """Make unavailable legacy observations explicit rather than silently absent."""
    outcome = fixture.setdefault("outcome", {"status": "observed"})
    outcome.setdefault("exception", None)
    fixture.setdefault("record_ids", [])
    fixture.setdefault("pre_state", None)
    fixture.setdefault("post_state", fixture["pre_state"])
    fixture.setdefault("candidates", {})
    fixture.setdefault(
        "propensity_vector",
        {
            "enabled_channel_ids": [],
            "propensities": {},
            "refused_propensities": {},
            "total_enabled": 0.0,
            "total_refused": 0.0,
            "availability": "not_applicable",
        },
    )
    fixture.setdefault(
        "met",
        {
            "availability": "not_applicable",
            "channels": {},
            "component_channels": {},
            "item_sites": {},
        },
    )
    fixture.setdefault(
        "thinning",
        {
            "availability": "not_applicable",
            "proposals": [],
            "rng_variates": [],
        },
    )
    fixture.setdefault(
        "clock", _clock_surface(fixture["name"], fixture["propensity_vector"])
    )
    if "observables" not in fixture:
        post = fixture["post_state"]
        fixture["observables"] = (
            {
                key: post[key]
                for key in (
                    "moments",
                    "total_ledger",
                    "total_mass",
                    "total_radicals",
                    "sink",
                    "sink_radicals",
                )
            }
            if isinstance(post, Mapping)
            else {"availability": "not_applicable"}
        )
    return fixture


def _validate_observation_surface(
    fixtures: Sequence[Mapping[str, Any]], records: Mapping[str, Any]
) -> None:
    required = {
        "outcome",
        "pre_state",
        "post_state",
        "candidates",
        "propensity_vector",
        "met",
        "thinning",
        "clock",
        "observables",
    }
    for fixture in fixtures:
        assert required <= set(fixture), fixture["name"]
        assert "status" in fixture["outcome"], fixture["name"]
        assert "exception" in fixture["outcome"], fixture["name"]
        assert set(fixture["record_ids"]) <= set(records), fixture["name"]
    assert not any(event_id.startswith("i052-") for event_id in records)


def generate_observations(*, remove_attachment: bool = False) -> dict[str, Any]:
    with ARTIFACT.open() as stream:
        artifact = json.load(stream)
    assert artifact["provenance"]["rmgpy_sha"] == RMG_PY_PIN
    assert artifact["provenance"]["rmg_database_sha"] == RMG_DATABASE_PIN
    records: dict[str, Any] = {}
    fixtures = []
    fixtures.extend(_mapped_fixtures(artifact, records))
    fixtures.extend(_attachment_fixtures(artifact, records))
    fixtures.extend(_plp_fixtures(artifact, records))
    fixtures.append(_met_fixture(artifact, records))
    fixtures.append(_spin_fixture())
    fixtures.append(_ghost_fixture(records))
    fixtures.append(_mass_fixture())
    fixtures.append(_rollback_fixture(artifact, records))
    fixtures = [_complete_legacy_surface(fixture) for fixture in fixtures]
    _validate_observation_surface(fixtures, records)
    if remove_attachment:
        target = next(
            fixture
            for fixture in fixtures
            if fixture["name"] == "J_para_after_association"
        )
        target["post_state"]["junctions"] = []
    counts = Counter(fixture["category"] for fixture in fixtures)
    return _jsonable(
        {
            "schema_version": "kmc_legacy_baseline/1",
            "pins": {
                "rmgpy": RMG_PY_PIN,
                "rmg_database": RMG_DATABASE_PIN,
                "mom": MOM_PIN,
                "artifact_sha256": ARTIFACT.stem,
            },
            "fixture_counts": dict(sorted(counts.items())),
            "records": {key: records[key] for key in sorted(records)},
            "fixtures": fixtures,
        }
    )


def _tracked_files(root: Path, *, include_untracked: bool = False):
    arguments = ["ls-files", "-z"]
    if include_untracked:
        arguments = ["ls-files", "-z", "--cached", "--others", "--exclude-standard"]
    output = subprocess.check_output(["git", "-C", str(root), *arguments])
    return [Path(value.decode()) for value in output.split(b"\0") if value]


def _materialized_files(root: Path):
    return [
        path.relative_to(root)
        for path in root.rglob("*")
        if path.is_file() and ".git" not in path.relative_to(root).parts
    ]


def build_manifest() -> dict[str, Any]:
    repository_hashes = {}
    for relative in sorted(_tracked_files(ROOT, include_untracked=True)):
        path = ROOT / relative
        if not path.is_file() or path == MANIFEST:
            continue
        repository_hashes[str(relative)] = _sha256(path)
    database_hashes = {
        str(relative): _sha256(DATABASE / relative)
        for relative in sorted(_materialized_files(DATABASE))
        if (DATABASE / relative).is_file()
    }
    return {
        "schema_version": "kmc_legacy_baseline_manifest/1",
        "pins": {
            "rmgpy_commit": RMG_PY_PIN,
            "rmg_database_commit": RMG_DATABASE_PIN,
            "mom_commit": MOM_PIN,
        },
        "artifact": {"path": str(ARTIFACT), "sha256": _sha256(ARTIFACT)},
        "design": {"path": str(DESIGN), "sha256": _sha256(DESIGN)},
        "repository_root": str(ROOT),
        "repository_files": repository_hashes,
        "database_root": str(DATABASE),
        "database_files": database_hashes,
        "environment": environment_pins(),
        "required_runtime_sources": {
            name: repository_hashes[f"rmgpy/kmc/{name}"]
            for name in (
                "state.py",
                "ssa.py",
                "met.py",
                "compiler.py",
                "event_record.py",
                "kinetics_library.py",
            )
        },
        "fixture": {
            "path": str(FIXTURE.relative_to(ROOT)),
            "sha256": repository_hashes[str(FIXTURE.relative_to(ROOT))],
        },
    }


def sysconfig_soabi():
    import sysconfig

    return sysconfig.get_config_var("SOABI")


def environment_pins():
    return {
        "python_executable": sys.executable,
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "python_soabi": sysconfig_soabi(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "networkx": networkx.__version__,
        "rmgpy_version": getattr(rmgpy, "__version__", None),
        "platform": platform.platform(),
        "byteorder": sys.byteorder,
    }


def verify_manifest(manifest: Mapping[str, Any]) -> None:
    _require_git_commit(ROOT, RMG_PY_PIN)
    _require_git_commit(ROOT, STEP1_FREEZE_COMMIT)
    _verify_frozen_json_bytes()
    assert manifest["pins"] == {
        "rmgpy_commit": RMG_PY_PIN,
        "rmg_database_commit": RMG_DATABASE_PIN,
        "mom_commit": MOM_PIN,
    }
    assert _sha256(Path(manifest["artifact"]["path"])) == manifest["artifact"]["sha256"]
    assert _sha256(Path(manifest["design"]["path"])) == manifest["design"]["sha256"]
    assert _git(ROOT, "merge-base", RMG_PY_PIN, "HEAD") == RMG_PY_PIN
    with ARTIFACT.open() as stream:
        artifact_provenance = json.load(stream)["provenance"]
    assert artifact_provenance["rmg_database_sha"] == RMG_DATABASE_PIN
    assert artifact_provenance["rmgpy_sha"] == RMG_PY_PIN
    assert manifest["environment"] == environment_pins()
    for relative, expected in manifest["repository_files"].items():
        # The freeze commit adds exactly three Step-1 record/harness files to
        # the pinned source tree. Preserve their recorded hashes from that
        # immutable overlay; every pre-existing input comes directly from the
        # RMG-Py pin rather than the evolving checkout.
        source_commit = (
            STEP1_FREEZE_COMMIT
            if relative in STEP1_FREEZE_ONLY_FILES
            else RMG_PY_PIN
        )
        assert _sha256_bytes(_git_blob(ROOT, source_commit, relative)) == expected, relative
    for relative, expected in manifest["database_files"].items():
        assert _sha256(DATABASE / relative) == expected, relative


def compare_frozen(expected: Any, actual: Any, path: str = "root") -> None:
    """Exact fields, exact zero, and relative 1e-12 for nonzero numerics."""
    if isinstance(expected, bool) or expected is None or isinstance(expected, str):
        assert actual == expected, f"{path}: {actual!r} != {expected!r}"
        return
    if isinstance(expected, (int, float)) and not isinstance(expected, bool):
        assert isinstance(actual, (int, float)) and not isinstance(actual, bool), path
        if expected == 0:
            assert actual == 0, f"{path}: baseline zero changed to {actual!r}"
        else:
            relative = abs((actual - expected) / expected)
            assert relative <= 1.0e-12, f"{path}: relative delta {relative:.17g}"
        return
    if isinstance(expected, list):
        assert isinstance(actual, list), path
        assert len(actual) == len(expected), f"{path}: list length changed"
        for index, (left, right) in enumerate(zip(expected, actual)):
            compare_frozen(left, right, f"{path}[{index}]")
        return
    if isinstance(expected, dict):
        assert isinstance(actual, dict), path
        assert set(actual) == set(expected), f"{path}: keys changed"
        for key in expected:
            compare_frozen(expected[key], actual[key], f"{path}.{key}")
        return
    assert actual == expected, f"{path}: {actual!r} != {expected!r}"


def write_freeze() -> None:
    FIXTURE.parent.mkdir(parents=True, exist_ok=True)
    observations = generate_observations()
    FIXTURE.write_text(json.dumps(observations, indent=2, sort_keys=True) + "\n")
    manifest = build_manifest()
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--negative-control", action="store_true")
    arguments = parser.parse_args(argv)
    if arguments.write:
        write_freeze()
        return 0
    if arguments.check or arguments.negative_control:
        expected = json.loads(FIXTURE.read_text())
        actual = generate_observations(remove_attachment=arguments.negative_control)
        compare_frozen(expected, actual)
        verify_manifest(json.loads(MANIFEST.read_text()))
        return 0
    parser.error("choose --write, --check, or --negative-control")


if __name__ == "__main__":
    raise SystemExit(main())
