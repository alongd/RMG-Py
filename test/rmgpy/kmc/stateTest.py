"""K2 acceptance checks for sparse strand/component record execution."""

import copy
import hashlib
import json
import multiprocessing
import os
import random
import subprocess
import sys
from collections import Counter
from itertools import product
from pathlib import Path

import networkx as nx
import pytest

from rmgpy.kmc.compiler import apply_record, validate_artifact
from rmgpy.kmc.event_record import EventRecord
from rmgpy.kmc.state import (
    ArityError,
    AtomRef,
    GraphMismatchError,
    KMCState,
    LedgerError,
    ReverseBindingError,
    Site,
    Strand,
    UUIDReuseError,
    _parse_adjacency,
)
from rmgpy.molecule.molecule import Molecule


ROOT = Path(__file__).resolve().parents[3]
DATABASE = Path(os.environ.get("RMG_DATABASE_PATH", ROOT.parent / "RMG-database"))
CACHE = (
    Path(os.environ.get("RMG_KMC_CACHE_ROOT", str(ROOT / ".kmc-cache"))) / "event-set"
)
_STRESS_ARTIFACT = None


@pytest.fixture(scope="session")
def ps_artifact():
    """Compile once per session, cached by both repository inputs."""
    from rmgpy.kmc.compiler import compiler_source_hash
    from cache_provenance import supplied_artifact

    supplied = supplied_artifact()
    if supplied:
        return supplied[1]

    source_hash = compiler_source_hash()
    head = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    database_head = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=DATABASE, text=True
    ).strip()
    CACHE.mkdir(parents=True, exist_ok=True)
    artifact_path = None
    for candidate in CACHE.glob("*.json"):
        artifact = json.loads(candidate.read_bytes())
        provenance = artifact.get("provenance", {})
        if (
            provenance.get("rmgpy_sha") == head
            and provenance.get("rmg_database_sha") == database_head
            and provenance.get("compiler_sources_sha256") == source_hash
        ):
            artifact_path = candidate
            break
    if artifact_path is None:
        fixture = Path(__file__).with_name("compile_event_set_fixture.py")
        environment = os.environ.copy()
        environment.update(
            {
                "PYTHONPATH": str(ROOT),
                "MPLCONFIGDIR": str(CACHE / "matplotlib"),
            }
        )
        subprocess.run(
            [sys.executable, str(fixture), str(DATABASE), str(CACHE)],
            cwd=ROOT,
            env=environment,
            check=True,
        )
        candidates = [
            path
            for path in CACHE.glob("*.json")
            if (
                json.loads(path.read_bytes()).get("provenance", {}).get("rmgpy_sha")
                == head
                and json.loads(path.read_bytes())
                .get("provenance", {})
                .get("rmg_database_sha")
                == database_head
                and json.loads(path.read_bytes())
                .get("provenance", {})
                .get("compiler_sources_sha256")
                == source_hash
            )
        ]
        assert len(candidates) == 1
        artifact_path = candidates[0]
    artifact = json.loads(artifact_path.read_bytes())
    validate_artifact(artifact)
    assert hashlib.sha256(artifact_path.read_bytes()).hexdigest() == artifact_path.stem
    return artifact


def _expanded_types(record):
    return [
        site_type
        for site_type, count in zip(
            record.participant_site_types, record.reactant_multiplicities
        )
        for _ in range(count)
    ]


def _harness(data, *, roles=None, labels=None, namespace=""):
    """Build current-state participant sites from an actual record reactant."""
    record = EventRecord.from_dict(data) if isinstance(data, dict) else data
    roles, labels = roles or {}, labels or {}
    graphs = [_parse_adjacency(text) for text in record.reactant_graphs]
    strands, sites = [], []
    offset = 0
    types = _expanded_types(record)
    for ordinal, (graph, site_type) in enumerate(zip(graphs, types)):
        strand_id = f"{namespace}strand-{ordinal}"
        length = max(8, len(graph) + 1)
        refs = {}
        formula = Counter()
        radicals = 0
        for local, node in sorted(graph.items()):
            role = roles.get(offset + local, "feature")
            position = roles.get((offset + local, "position"), local)
            ref = AtomRef(
                f"{namespace}{record.event_id}-{ordinal}-{local}",
                strand_id,
                position,
                role,
            )
            refs[local] = ref
            formula[node["element"]] += 1
            radicals += int(node.get("radical", 0))
        atom_graph = {}
        for local, node in graph.items():
            atom_graph[refs[local].uuid] = {
                "element": node["element"],
                "radical": int(node.get("radical", 0)),
                "charge": int(node.get("charge", 0)),
                "lone_pairs": int(node.get("lone_pairs", 0)),
                "implicit_hydrogens": int(node.get("implicit_hydrogens", 0)),
                "edges": {
                    refs[neighbor].uuid: order
                    for neighbor, order in node.get("edges", {}).items()
                },
            }
        strands.append(
            Strand(strand_id, length, {}, dict(formula), radicals, refs, atom_graph)
        )
        actual = copy.deepcopy(graph)
        for local, node in actual.items():
            ref = refs[local]
            node.update(
                {
                    "atom_ref": ref,
                    "position": ref.position,
                    "role": ref.role,
                }
            )
            if offset + local in labels:
                node["label"] = labels[offset + local]
        sites.append(Site(site_type, strand_id, 0, 0, actual))
        offset += len(graph)
    state = KMCState(strands)
    current_sites = []
    for site in sites:
        graph = copy.deepcopy(site.graph)
        for node in graph.values():
            node["atom_ref"] = state._ref(node["atom_ref"].uuid)
        current_sites.append(
            Site(site.site_type, site.strand_id, site.position, site.atom, graph)
        )
    return state, record, current_sites


def _nx_executor(component):
    graph = nx.Graph()
    for (
        atom_uuid,
        element,
        radical,
        charge,
        lone_pairs,
        implicit_hydrogens,
        edges,
    ) in component:
        graph.add_node(
            atom_uuid,
            element=element,
            radical=radical,
            charge=charge,
            lone_pairs=lone_pairs,
            implicit_hydrogens=implicit_hydrogens,
        )
        for neighbor, order in edges:
            graph.add_edge(atom_uuid, neighbor, order=float(order))
    return graph


def _nx_molecule(molecule):
    graph = nx.Graph()
    for atom in molecule.atoms:
        graph.add_node(
            id(atom),
            element=atom.element.symbol,
            radical=atom.radical_electrons,
            charge=atom.charge,
            lone_pairs=atom.lone_pairs,
            implicit_hydrogens=getattr(atom, "implicit_hydrogens", 0),
        )
    for bond in molecule.get_all_edges():
        graph.add_edge(id(bond.vertex1), id(bond.vertex2), order=float(bond.order))
    return graph


def _same_products(entry, record):
    actual = [_nx_executor(component) for component in entry.product_graphs]
    participants = [
        Molecule().from_adjacency_list(text) for text in record.reactant_graphs
    ]
    expected = [
        _nx_molecule(molecule) for molecule in apply_record(record, participants)
    ]
    node_match = nx.algorithms.isomorphism.categorical_node_match(
        ["element", "radical", "charge", "lone_pairs", "implicit_hydrogens"],
        [None, 0, 0, 0, 0],
    )
    edge_match = nx.algorithms.isomorphism.numerical_edge_match("order", 1.0)
    unmatched = list(expected)
    for graph in actual:
        for index, candidate in enumerate(unmatched):
            if nx.is_isomorphic(
                graph, candidate, node_match=node_match, edge_match=edge_match
            ):
                del unmatched[index]
                break
        else:
            return False
    return not unmatched


def _assert_record_balance(state, record, before_formula, before_radicals, before_mass):
    expected_formula = Counter(before_formula)
    expected_formula.update(record.formula_delta)
    assert state.total_ledger == {
        key: value for key, value in sorted(expected_formula.items()) if value
    }
    assert state.total_radicals == before_radicals + record.radical_delta
    expected_mass = before_mass + record.mass_delta
    assert abs(state.total_mass - expected_mass) <= 1e-12 * max(1.0, abs(expected_mass))
    state.assert_component_consistency()
    assert state.moments == state.brute_force_moments()


def _assert_strand_radicals_match_graph(state):
    ledger = {
        strand_id: strand.radical_count
        for strand_id, strand in sorted(state.strands.items())
    }
    graph = {
        strand_id: sum(
            int(node.get("radical", 0)) for node in strand.atom_graph.values()
        )
        for strand_id, strand in sorted(state.strands.items())
    }
    assert (
        ledger == graph
    ), f"strand radical recount mismatch: ledger={ledger}, graph={graph}"


def _formula_with(formula, element, amount):
    updated = Counter(formula)
    updated[element] += amount
    return {key: value for key, value in sorted(updated.items()) if value}


_EDGE_MATCH = nx.algorithms.isomorphism.numerical_edge_match("order", 1.0)
_ANCHORED_NODE_MATCH = nx.algorithms.isomorphism.categorical_node_match(
    [
        "element",
        "radical",
        "charge",
        "lone_pairs",
        "implicit_hydrogens",
        "anchor",
    ],
    [None, 0, 0, 0, 0, False],
)


def _state_graph_components(state):
    nodes = {}
    for strand in state.strands.values():
        for atom_uuid, node in strand.atom_graph.items():
            nodes[atom_uuid] = node
    unseen = set(nodes)
    components = []
    while unseen:
        root = min(unseen)
        members = set()
        stack = [root]
        while stack:
            current = stack.pop()
            if current in members:
                continue
            members.add(current)
            stack.extend(set(nodes[current].get("edges", {})) - members)
        unseen -= members
        components.append({atom_uuid: nodes[atom_uuid] for atom_uuid in members})
    return components


def _nx_current(component):
    graph = nx.Graph()
    for atom_uuid, node in component.items():
        graph.add_node(
            atom_uuid,
            element=node["element"],
            radical=int(node.get("radical", 0)),
            charge=int(node.get("charge", 0)),
            lone_pairs=int(node.get("lone_pairs", 0)),
            implicit_hydrogens=int(node.get("implicit_hydrogens", 0)),
        )
        for neighbor, order in node.get("edges", {}).items():
            if neighbor in component:
                graph.add_edge(atom_uuid, neighbor, order=float(order))
    return graph


def _nx_expected(graph):
    expected = nx.Graph()
    for index, node in graph.items():
        expected.add_node(
            index,
            element=node["element"],
            radical=int(node.get("radical", 0)),
            charge=int(node.get("charge", 0)),
            lone_pairs=int(node.get("lone_pairs", 0)),
            implicit_hydrogens=int(node.get("implicit_hydrogens", 0)),
        )
        for neighbor, order in node.get("edges", {}).items():
            expected.add_edge(index, neighbor, order=float(order))
    return expected


def _raw_graph_signature(graph):
    nodes = Counter(
        (
            data["element"],
            int(data.get("radical", 0)),
            int(data.get("charge", 0)),
            int(data.get("lone_pairs", 0)),
            int(data.get("implicit_hydrogens", 0)),
        )
        for data in graph.values()
    )
    edges = Counter()
    seen = set()
    for node, data in graph.items():
        for neighbor, order in data.get("edges", {}).items():
            if neighbor in graph and frozenset((node, neighbor)) not in seen:
                seen.add(frozenset((node, neighbor)))
                edges[float(order)] += 1
    return (
        len(graph),
        sum(edges.values()),
        tuple(sorted(nodes.items())),
        tuple(sorted(edges.items())),
    )


def _signature_can_contain(container, expected):
    container_nodes, container_edges = dict(container[2]), dict(container[3])
    return (
        container[0] >= expected[0]
        and container[1] >= expected[1]
        and all(container_nodes.get(key, 0) >= count for key, count in expected[2])
        and all(container_edges.get(key, 0) >= count for key, count in expected[3])
    )


def _family_seed_record(artifact, family, seed):
    if family == "R_Addition_MultipleBond":
        return _scission_record(artifact)
    candidates = [
        data
        for data in artifact["records"]
        if data["family"] == family and data["inventory_class"] != "R1:J_ring"
    ]
    if not candidates:
        return None
    return random.Random(f"{seed}:{family}").choice(
        sorted(candidates, key=lambda data: data["event_id"])
    )


def _stress_record_data(artifact, family, seed):
    capture_data, reverse_data = _r1_pair(artifact)
    selected = [capture_data, reverse_data]
    family_record = _family_seed_record(artifact, family, seed)
    if family_record is not None:
        selected.append(family_record)
    return tuple({data["event_id"]: data for data in selected}.values())


def _prepare_stress_record(data):
    record = EventRecord.from_dict(data)
    graphs = [_parse_adjacency(text) for text in record.reactant_graphs]
    expected = [_nx_expected(graph) for graph in graphs]
    offsets = []
    offset = 0
    for graph in graphs:
        offsets.append(offset)
        offset += len(graph)
    return {
        "record": record,
        "graphs": graphs,
        "expected": expected,
        "signatures": [_raw_graph_signature(graph) for graph in graphs],
        "labels": _record_labels(record),
        "types": _expanded_types(record),
        "offsets": offsets,
    }


def _stress_records(artifact, family, seed, sampled, cache):
    selected = list(_stress_record_data(artifact, family, seed)) + list(sampled)
    prepared = []
    seen = set()
    for data in selected:
        if data["event_id"] in seen:
            continue
        seen.add(data["event_id"])
        if data["event_id"] not in cache:
            cache[data["event_id"]] = _prepare_stress_record(data)
        prepared.append(cache[data["event_id"]])
    return tuple(prepared)


def _seed_evolving_state(artifact, family, seed):
    rng = random.Random(seed)
    capture_data, _ = _r1_pair(artifact)
    form_atoms = next(
        operation["atoms"]
        for operation in capture_data["bond_ops"]
        if operation["action"] == "form"
    )
    capture_roles = {form_atoms[0]: "attacked_ring", form_atoms[1]: "attacker"}
    capture_labels = {form_atoms[0]: "S9", form_atoms[1]: "P9"}
    capture_state, capture, capture_sites = _harness(
        capture_data,
        roles=capture_roles,
        labels=capture_labels,
        namespace=f"chain-{seed}-r1-",
    )
    strands = [copy.deepcopy(strand) for strand in capture_state.strands.values()]
    family_record = _family_seed_record(artifact, family, seed)
    if family_record is not None:
        roles = {}
        if family == "R_Addition_MultipleBond":
            cut = next(
                operation["atoms"]
                for operation in family_record["bond_ops"]
                if operation.get("action") == "break"
                and operation.get("atoms") == [46, 47]
            )
            roles = {
                cut[0]: "backbone",
                cut[1]: "backbone",
                (cut[0], "position"): 3,
                (cut[1], "position"): 4,
            }
        family_state, _, _ = _harness(
            family_record,
            roles=roles,
            namespace=f"chain-{seed}-{family}-",
        )
        strands.extend(
            copy.deepcopy(strand) for strand in family_state.strands.values()
        )
    previous_length = 0
    for strand in strands:
        seeded_length = strand.length + rng.randrange(1, 12)
        strand.length = max(seeded_length, previous_length + 1)
        previous_length = strand.length

    # Embed each R1 site in seeded local context so candidate discovery must
    # find anchored subgraphs rather than exact record-shaped graph islands.
    context_strand = strands[rng.randrange(2)]
    context_position = context_strand.length - 1
    attach_uuid = rng.choice(sorted(context_strand.atom_graph))
    for context_index in range(2 + rng.randrange(3)):
        context_uuid = f"chain-{seed}-context-{context_index}"
        context_ref = AtomRef(
            context_uuid, context_strand.strand_id, context_position, "feature"
        )
        context_strand.atom_refs[max(context_strand.atom_refs) + 1] = context_ref
        context_strand.atom_graph[context_uuid] = {
            "element": "H",
            "radical": 0,
            "charge": 0,
            "lone_pairs": 0,
            "implicit_hydrogens": 0,
            "edges": {attach_uuid: "1.0"},
        }
        context_strand.atom_graph[attach_uuid]["edges"][context_uuid] = "1.0"
        context_strand.formula = _formula_with(context_strand.formula, "H", 1)
        attach_uuid = context_uuid

    state = KMCState(strands)
    state.apply(capture, capture_sites)
    assert len(state.strands) >= 3
    assert len({strand.length for strand in state.strands.values()}) >= 2
    assert sum(strand.radical_count for strand in state.strands.values()) > 0 or (
        family_record is not None and family_record["radical_delta"] > 0
    )
    assert state.junctions
    return state


def _record_labels(record):
    labels = {}
    for junction in record.junction_ops or []:
        action = "form" if junction.get("action") == "create" else "break"
        operation = next(
            (
                item
                for item in record.bond_ops
                if item.get("action") == action and len(item.get("atoms", [])) == 2
            ),
            None,
        )
        pair = junction.get("formed_bond", {}).get("label_pair", [])
        if operation is not None and len(pair) == 2:
            labels.update(zip(operation["atoms"], pair))
    return labels


def _mapping_still_matches(expected, component, mapping):
    if set(mapping) != set(expected) or not set(mapping.values()) <= set(component):
        return False
    mapped_uuids = set(mapping.values())
    for index, expected_node in expected.items():
        current_node = component[mapping[index]]
        if any(
            int(current_node.get(field, 0)) != int(expected_node.get(field, 0))
            for field in (
                "radical",
                "charge",
                "lone_pairs",
                "implicit_hydrogens",
            )
        ) or current_node.get("element") != expected_node.get("element"):
            return False
        expected_edges = {
            mapping[neighbor]: float(order)
            for neighbor, order in expected_node.get("edges", {}).items()
        }
        current_edges = {
            neighbor: float(order)
            for neighbor, order in current_node.get("edges", {}).items()
            if neighbor in mapped_uuids
        }
        if expected_edges != current_edges:
            return False
    return True


def _anchored_mappings(expected_graph, expected, component, reverse_captures, labels):
    anchor = min(expected_graph)
    anchor_node = expected_graph[anchor]
    current = _nx_current(component)
    nx.set_node_attributes(current, False, "anchor")
    nx.set_node_attributes(expected, False, "anchor")
    expected.nodes[anchor]["anchor"] = True
    mappings = []
    for atom_uuid, node in sorted(component.items()):
        if node.get("element") != anchor_node.get("element") or any(
            int(node.get(field, 0)) != int(anchor_node.get(field, 0))
            for field in (
                "radical",
                "charge",
                "lone_pairs",
                "implicit_hydrogens",
            )
        ):
            continue
        current.nodes[atom_uuid]["anchor"] = True
        matcher = nx.algorithms.isomorphism.GraphMatcher(
            current,
            expected,
            node_match=_ANCHORED_NODE_MATCH,
            edge_match=_EDGE_MATCH,
        )
        for current_mapping in matcher.subgraph_isomorphisms_iter():
            candidate = {
                expected_atom: current_atom
                for current_atom, expected_atom in current_mapping.items()
            }
            if not _mapping_still_matches(expected_graph, component, candidate):
                continue
            if reverse_captures and not any(
                all(
                    candidate[index] == capture["bindings"].get(label)
                    for index, label in labels.items()
                )
                for capture in reverse_captures
            ):
                continue
            mappings.append(candidate)
            break
        current.nodes[atom_uuid]["anchor"] = False
    expected.nodes[anchor]["anchor"] = False
    return mappings


def _component_match_key(component):
    return tuple(
        (
            atom_uuid,
            node["element"],
            int(node.get("radical", 0)),
            int(node.get("charge", 0)),
            int(node.get("lone_pairs", 0)),
            int(node.get("implicit_hydrogens", 0)),
            tuple(
                sorted(
                    (neighbor, float(order))
                    for neighbor, order in node.get("edges", {}).items()
                    if neighbor in component
                )
            ),
        )
        for atom_uuid, node in sorted(component.items())
    )


def _current_candidates(state, prepared_records, match_cache):
    """Brute-force anchored record/site matches over state-owned current graphs."""
    components = _state_graph_components(state)
    signatures = [_raw_graph_signature(component) for component in components]
    component_keys = []
    for component in components:
        component_key = _component_match_key(component)
        component_keys.append(
            match_cache.setdefault(("component", component_key), component_key)
        )
    candidates = []
    for prepared in prepared_records:
        record = prepared["record"]
        expected_graphs = prepared["graphs"]
        labels = prepared["labels"]
        reverse_captures = (
            state._open_captures.get(record.reverse_of, [])
            if record.orientation == "reversed"
            and record.inventory_class == "R1:J_ring"
            else ()
        )
        match_groups = []
        for ordinal, expected in enumerate(prepared["expected"]):
            signature = prepared["signatures"][ordinal]
            matches = []
            for component_index, component in enumerate(components):
                if not _signature_can_contain(signatures[component_index], signature):
                    continue
                cache_key = (record.event_id, ordinal, tuple(sorted(component)))
                mappings = [
                    mapping
                    for mapping in match_cache.get(cache_key, [])
                    if _mapping_still_matches(
                        expected_graphs[ordinal], component, mapping
                    )
                ]
                missing_key = (
                    "no-match",
                    record.event_id,
                    ordinal,
                    component_keys[component_index],
                    tuple(
                        tuple(sorted(capture["bindings"].items()))
                        for capture in reverse_captures
                    ),
                )
                if not mappings and missing_key not in match_cache:
                    mappings = _anchored_mappings(
                        expected_graphs[ordinal],
                        expected,
                        component,
                        reverse_captures,
                        labels,
                    )
                    if mappings:
                        match_cache[cache_key] = mappings
                    else:
                        match_cache[missing_key] = None
                matches.extend((component_index, mapping) for mapping in mappings)
            match_groups.append(matches)
        if any(not matches for matches in match_groups):
            continue
        types = prepared["types"]
        offsets = prepared["offsets"]
        for combination in product(*match_groups):
            mapped_groups = [set(item[1].values()) for item in combination]
            if any(
                left & right
                for index, left in enumerate(mapped_groups)
                for right in mapped_groups[index + 1 :]
            ):
                continue
            global_mapping = {
                offsets[ordinal] + local: atom_uuid
                for ordinal, (_, mapping) in enumerate(combination)
                for local, atom_uuid in mapping.items()
            }
            if (
                record.orientation == "reversed"
                and record.inventory_class == "R1:J_ring"
            ):
                captures = state._open_captures.get(record.reverse_of, [])
                if not any(
                    all(
                        global_mapping[index] == capture["bindings"].get(label)
                        for index, label in labels.items()
                    )
                    for capture in captures
                ):
                    continue
            sites = []
            for ordinal, expected_graph in enumerate(expected_graphs):
                mapping = combination[ordinal][1]
                site_graph = copy.deepcopy(expected_graph)
                for local, node in site_graph.items():
                    ref = state._ref(mapping[local])
                    node.update(
                        {
                            "atom_ref": ref,
                            "position": ref.position,
                            "role": ref.role,
                        }
                    )
                    label = labels.get(offsets[ordinal] + local)
                    if label:
                        node["label"] = label
                anchor = min(site_graph)
                anchor_ref = site_graph[anchor]["atom_ref"]
                sites.append(
                    Site(
                        types[ordinal],
                        anchor_ref.strand_id,
                        anchor_ref.position,
                        anchor,
                        site_graph,
                    )
                )
            candidates.append((record, sites))
    return candidates


def _exercise_evolving_chain(artifact, family, count, seed, before_apply=None):
    rng = random.Random(seed)
    record_cache = {}
    match_cache = {}
    state = _seed_evolving_state(artifact, family, seed)
    stats = Counter(chains=1)
    for _ in range(count):
        sampled = rng.sample(artifact["records"], 1)
        records = _stress_records(
            artifact, family, seed + stats["reseeds"], sampled, record_cache
        )
        stats["record_samples"] += len(records)
        candidates = _current_candidates(state, records, match_cache)
        if not candidates:
            stats["reseeds"] += 1
            stats["chains"] += 1
            state = _seed_evolving_state(artifact, family, seed + stats["reseeds"])
            match_cache = {}
            records = _stress_records(
                artifact, family, seed + stats["reseeds"], sampled, record_cache
            )
            candidates = _current_candidates(state, records, match_cache)
        assert candidates
        family_record = _family_seed_record(artifact, family, seed + stats["reseeds"])
        target_candidates = [
            candidate
            for candidate in candidates
            if family_record is not None
            and candidate[0].event_id == family_record["event_id"]
        ]
        record, sites = rng.choice(
            target_candidates
            if not stats["target_events"] and target_candidates
            else candidates
        )
        before_formula = state.total_ledger
        before_radicals = state.total_radicals
        before_mass = state.total_mass
        before_components = len(set(state.components.values()))
        if before_apply is not None:
            before_apply(state, record, sites)
        entry = state.apply(record, sites)
        _assert_record_balance(
            state, record, before_formula, before_radicals, before_mass
        )
        _assert_strand_radicals_match_graph(state)
        after_components = len(set(state.components.values()))
        stats["steps"] += 1
        stats["splits"] += len(entry.derived_cut_offsets)
        stats["component_splits"] += max(0, after_components - before_components)
        stats["merges"] += max(0, before_components - after_components)
        stats[f"family:{record.family}"] += 1
        stats[f"event:{record.event_id}"] += 1
        if family_record is not None and record.event_id == family_record["event_id"]:
            stats["target_events"] += 1
        if record.inventory_class == "R1:J_ring":
            action = record.junction_ops[0]["action"]
            stats["r1_captures" if action == "create" else "r1_releases"] += 1
    return dict(stats)


def _slow_stress_worker(arguments):
    worker, count, family = arguments
    return _exercise_evolving_chain(_STRESS_ARTIFACT, family, count, 8675309 + worker)


def _slow_stress_jobs(total_steps, families, shards=24):
    """Return the stable million-step stress-test shard plan."""
    base, remainder = divmod(total_steps, shards)
    return [
        (worker, base + (worker < remainder), families[worker % len(families)])
        for worker in range(shards)
    ]


def _slow_stress_processes(shards):
    """Return the bounded process-pool size for the slow stress test."""
    configured = os.environ.get("RMG_KMC_STRESS_PROCESSES")
    processes = (
        min(4, multiprocessing.cpu_count()) if configured is None else int(configured)
    )
    return min(processes, shards)


def test_slow_stress_shard_plan_and_pool_size(monkeypatch):
    families = ["alpha", "beta", "gamma"]
    assert _slow_stress_jobs(1_000_000, families) == [
        (worker, 41_667 if worker < 16 else 41_666, families[worker % 3])
        for worker in range(24)
    ]

    monkeypatch.setattr(multiprocessing, "cpu_count", lambda: 8)
    monkeypatch.delenv("RMG_KMC_STRESS_PROCESSES", raising=False)
    assert _slow_stress_processes(24) == 4

    monkeypatch.setenv("RMG_KMC_STRESS_PROCESSES", "3")
    assert _slow_stress_processes(24) == 3

    monkeypatch.setenv("RMG_KMC_STRESS_PROCESSES", "30")
    assert _slow_stress_processes(24) == 24


def test_zero_radical_seed_fires_real_selected_initiator(ps_artifact, monkeypatch):
    record_data = next(
        record
        for record in ps_artifact["records"]
        if record["family"] == "R_Recombination"
        and record["arity"] == 1
        and record["radical_delta"] == 2
        and record["participant_site_types"] == ["pristine"]
    )
    monkeypatch.setattr(
        sys.modules[__name__], "_family_seed_record", lambda *arguments: record_data
    )
    state = _seed_evolving_state(ps_artifact, "R_Recombination", 8_675_317)
    assert state.total_radicals == 0
    candidates = _current_candidates(state, [_prepare_stress_record(record_data)], {})
    assert candidates
    record, sites = candidates[0]
    before_formula = state.total_ledger
    before_radicals = state.total_radicals
    before_mass = state.total_mass
    state.apply(record, sites)
    _assert_record_balance(state, record, before_formula, before_radicals, before_mass)
    _assert_strand_radicals_match_graph(state)
    assert state.total_radicals == 2


@pytest.mark.parametrize(
    "field",
    ["element", "radical", "charge", "lone_pairs", "implicit_hydrogens", "edges"],
)
def test_failed_match_key_tracks_every_matching_field(ps_artifact, field):
    state = _seed_evolving_state(ps_artifact, "Disproportionation", 8_675_309)
    component = _state_graph_components(state)[0]
    changed = copy.deepcopy(component)
    atom_uuid = next(atom_uuid for atom_uuid, node in changed.items() if node["edges"])
    if field == "element":
        changed[atom_uuid][field] = "Si"
    elif field == "edges":
        neighbor = next(iter(changed[atom_uuid]["edges"]))
        changed[atom_uuid]["edges"][neighbor] = "3.0"
        changed[neighbor]["edges"][atom_uuid] = "3.0"
    else:
        changed[atom_uuid][field] = int(changed[atom_uuid].get(field, 0)) + 1
    assert set(component) == set(changed)
    assert _component_match_key(component) != _component_match_key(changed)


def test_failed_match_cache_tracks_real_capture_bindings(ps_artifact, monkeypatch):
    state = _seed_evolving_state(ps_artifact, "Disproportionation", 8_675_309)
    capture_data, reverse_data = _r1_pair(ps_artifact)
    prepared = _prepare_stress_record(reverse_data)
    bindings = state._open_captures[capture_data["event_id"]][0]["bindings"]
    saved = dict(bindings)
    for label in bindings:
        bindings[label] = "not-an-owned-atom"
    original = _anchored_mappings
    attempts = []

    def counted(*arguments):
        attempts.append(1)
        return original(*arguments)

    monkeypatch.setattr(sys.modules[__name__], "_anchored_mappings", counted)
    cache = {}
    assert _current_candidates(state, [prepared], cache) == []
    assert attempts
    assert any(key[0] == "no-match" for key in cache)
    for key in cache:
        if key[0] == "no-match":
            assert cache[("component", key[3])] is key[3]
    attempts.clear()
    assert _current_candidates(state, [prepared], cache) == []
    assert not attempts
    bindings.update(saved)
    assert _current_candidates(state, [prepared], cache)
    assert attempts


def test_real_corpus_is_valid_and_self_authenticating(ps_artifact):
    assert ps_artifact["records"]
    assert {record["family"] for record in ps_artifact["records"]} == set(
        ps_artifact["families"]
    )


def test_executor_agrees_with_compiler_oracle_across_every_family(ps_artifact):
    samples = {}
    for data in ps_artifact["records"]:
        if data["inventory_class"] != "R1:J_ring":
            samples.setdefault(data["family"], data)
    assert set(samples) == set(ps_artifact["families"])
    for data in samples.values():
        state, record, sites = _harness(data)
        entry = state.apply(record, sites)
        assert _same_products(entry, record), record.event_id
        state.assert_component_consistency()
        assert state.moments == state.brute_force_moments()


def test_evolving_trajectory_discovers_sites_from_current_state(ps_artifact):
    trajectory = _exercise_evolving_chain(
        ps_artifact, "R_Addition_MultipleBond", 10, 8675309
    )
    assert trajectory["steps"] == 10
    assert trajectory["splits"] > 0
    assert trajectory["merges"] > 0
    assert trajectory["r1_captures"] > 0
    assert trajectory["r1_releases"] > 0


def _captured_coproduct_radical_case():
    path = Path(__file__).with_name("fixtures") / "i033_coproduct_radical_state.json"
    captured = json.loads(path.read_bytes())
    strands = []
    for data in captured["strands"]:
        strands.append(
            Strand(
                strand_id=data["strand_id"],
                length=data["length"],
                features={int(key): value for key, value in data["features"].items()},
                formula=data["formula"],
                radical_count=data["radical_count"],
                atom_refs={
                    int(key): AtomRef(**value)
                    for key, value in data["atom_refs"].items()
                },
                atom_graph=data["atom_graph"],
                repeat_unit_formula=data["repeat_unit_formula"],
            )
        )
    state = KMCState(strands)
    state.junctions = {tuple(pair) for pair in captured["junctions"]}
    state.sink = captured["sink"]
    state.sink_radicals = captured["sink_radicals"]
    state._strand_serial = captured["strand_serial"]
    state._refresh_components()

    provenance = captured["provenance"]
    assert captured["record_artifact_sha256"] == (
        "517e3decc530ae435b7a2e2dab647cb0a7047cd4b7b2269ae61b08ef0e2075a6"
    )
    record = EventRecord.from_dict(captured["record"])
    assert record.event_id == provenance["event_id"]
    graphs = [_parse_adjacency(text) for text in record.reactant_graphs]
    labels = _record_labels(record)
    sites = []
    offset = 0
    for saved, graph in zip(captured["sites"], graphs):
        for local, node in graph.items():
            ref = state._ref(saved["atom_uuids"][str(local)])
            node.update({"atom_ref": ref, "position": ref.position, "role": ref.role})
            label = labels.get(offset + local)
            if label:
                node["label"] = label
        sites.append(
            Site(
                saved["site_type"],
                saved["strand_id"],
                saved["position"],
                saved["atom"],
                graph,
            )
        )
        offset += len(graph)
    return state, record, sites, captured


def test_evolving_coproduct_radical_debit_matches_atom_graph():
    """Replay the debit that corrupted the graph, not a split/merge or bad record."""
    state, record, sites, captured = _captured_coproduct_radical_case()
    provenance = captured["provenance"]
    assert provenance == {
        "worker": 2,
        "seed": 8675311,
        "step_one_based": 414,
        "reseed_count": 1,
        "event_id": "evt_31c1877561e0abce21b079cb6e3c41349daef0fb4fb8410ae4c8ed866bdd85ff",
        "family": "R_Addition_MultipleBond",
        "inventory_class": None,
        "orientation": "reversed",
        "later_guard": {
            "step_one_based": 13472,
            "event_id": "evt_cdcfeb42b7a40698042ac7dbe31d5ae27eb91f4addb578cccb4c86853669fc06",
            "family": "R_Addition_MultipleBond",
            "inventory_class": None,
            "orientation": "reversed",
            "site": {
                "site_type": "end_radical",
                "strand_id": "chain-8675418-r1-strand-0",
                "position": 0,
                "atom": 0,
                "all_mapped_atoms_owner": "chain-8675418-r1-strand-0",
                "radical_atom_local": 25,
                "radical_atom_uuid": "chain-8675418-r1-evt_4131af12d1b5000b117b391c60415b9c82e1bd789875d727b36e16ec9d56c400-0-25",
            },
        },
    }
    assert record.family == provenance["family"]
    assert record.inventory_class is provenance["inventory_class"]
    assert [
        (site.site_type, site.strand_id, site.position, site.atom) for site in sites
    ] == [("end_radical", "chain-8675312-r1-strand-0", 1, 0)]
    assignment = json.dumps(
        captured["sites"], sort_keys=True, separators=(",", ":")
    ).encode()
    assert hashlib.sha256(assignment).hexdigest() == (
        "393d83ba6b965894041e95a29b69f978eec24cd415f4454f4898039f4e7f642b"
    )
    before = {
        strand_id: (
            strand.radical_count,
            sum(int(node.get("radical", 0)) for node in strand.atom_graph.values()),
        )
        for strand_id, strand in sorted(state.strands.items())
    }
    assert before == {
        "chain-8675312-R_Recombination-strand-0": (0, 0),
        "chain-8675312-R_Recombination-strand-1": (0, 0),
        "chain-8675312-r1-strand-0": (1, 1),
        "chain-8675312-r1-strand-1": (1, 1),
    }

    state.apply(record, sites)

    _assert_strand_radicals_match_graph(state)
    assert state.strands["chain-8675312-r1-strand-0"].radical_count == 0
    assert state.sink_radicals == 1


def test_evolving_recount_catches_wrong_coproduct_owner(monkeypatch):
    state, record, sites, _ = _captured_coproduct_radical_case()
    actual_owner = "chain-8675312-r1-strand-0"
    wrong_owner = "chain-8675312-r1-strand-1"
    assert state.strands[actual_owner].radical_count == 1
    assert state.strands[wrong_owner].radical_count == 1
    monkeypatch.setattr(state, "_coproduct_atom_owner", lambda atom_uuid: wrong_owner)

    state.apply(record, sites)

    with pytest.raises(AssertionError, match="strand radical recount mismatch"):
        _assert_strand_radicals_match_graph(state)


def test_slow_chain_seed_starts_with_distinct_strand_lengths(ps_artifact):
    state = _seed_evolving_state(ps_artifact, "intra_H_migration", 8675470)
    assert len({strand.length for strand in state.strands.values()}) >= 2


def test_ledger_components_and_moments_over_seeded_real_draws(ps_artifact):
    by_family = {}
    for data in ps_artifact["records"]:
        is_unbound_reverse = data["inventory_class"] == "R1:J_ring" and any(
            operation.get("action") == "dissociate"
            for operation in data.get("junction_ops", [])
        )
        if not is_unbound_reverse:
            by_family.setdefault(data["family"], []).append(data)
    assert set(by_family) == set(ps_artifact["families"])
    # Force the topology-changing corpus cases into this same stress path.
    scission_data = _scission_record(ps_artifact)
    scission_roles = {
        46: "backbone",
        47: "backbone",
        (46, "position"): 3,
        (47, "position"): 4,
    }
    scission_state, scission_record, scission_sites = _harness(
        scission_data, roles=scission_roles
    )
    scission_before = (
        scission_state.total_ledger,
        scission_state.total_radicals,
        scission_state.total_mass,
    )
    scission_state.apply(scission_record, scission_sites)
    _assert_record_balance(scission_state, scission_record, *scission_before)
    assert len(scission_state.strands) == 2

    capture_data, reverse_data = _r1_pair(ps_artifact)
    form_atoms = next(
        operation["atoms"]
        for operation in capture_data["bond_ops"]
        if operation["action"] == "form"
    )
    capture_roles = {form_atoms[0]: "attacked_ring", form_atoms[1]: "attacker"}
    capture_labels = {form_atoms[0]: "S9", form_atoms[1]: "P9"}
    r1_state, capture, capture_sites = _harness(
        capture_data, roles=capture_roles, labels=capture_labels
    )
    r1_before = r1_state.state_hash()
    capture_before = (
        r1_state.total_ledger,
        r1_state.total_radicals,
        r1_state.total_mass,
    )
    capture_entry = r1_state.apply(capture, capture_sites)
    _assert_record_balance(r1_state, capture, *capture_before)
    reverse = EventRecord.from_dict(reverse_data)
    reverse_site = _reverse_site(r1_state, reverse_data, capture_entry)
    reverse_before = (
        r1_state.total_ledger,
        r1_state.total_radicals,
        r1_state.total_mass,
    )
    r1_state.apply(reverse, [reverse_site])
    _assert_record_balance(r1_state, reverse, *reverse_before)
    assert r1_state.state_hash() == r1_before

    total_steps = 1_000_000 if os.environ.get("RMG_KMC_SLOW") == "1" else 10_000
    families = sorted(by_family)
    if os.environ.get("RMG_KMC_SLOW") == "1":
        global _STRESS_ARTIFACT
        _STRESS_ARTIFACT = ps_artifact
        jobs = _slow_stress_jobs(total_steps, families)
        processes = _slow_stress_processes(len(jobs))
        with multiprocessing.get_context("fork").Pool(
            processes, maxtasksperchild=1
        ) as pool:
            results = pool.map(_slow_stress_worker, jobs)
    else:
        base, remainder = divmod(total_steps, len(families))
        results = [
            _exercise_evolving_chain(
                ps_artifact,
                family,
                base + (ordinal < remainder),
                8675309 + ordinal,
            )
            for ordinal, family in enumerate(families)
        ]
    trajectory = Counter()
    for result in results:
        trajectory.update(result)
    assert trajectory["steps"] == total_steps
    assert trajectory["splits"] > 0
    assert trajectory["component_splits"] > 0
    assert trajectory["merges"] > 0
    assert trajectory["r1_captures"] > 0
    assert trajectory["r1_releases"] > 0
    assert {family for family in families if trajectory[f"family:{family}"]} == set(
        families
    )
    print(
        "K2 evolving trajectory: "
        + " ".join(
            f"{key}={trajectory[key]}"
            for key in (
                "chains",
                "steps",
                "reseeds",
                "splits",
                "component_splits",
                "merges",
                "r1_captures",
                "r1_releases",
                "record_samples",
            )
        )
        + f" families={len(families)}"
    )


def test_atomicity_when_failure_is_injected_mid_rewrite(ps_artifact, monkeypatch):
    data = next(
        record for record in ps_artifact["records"] if len(record["bond_ops"]) > 2
    )
    state, record, sites = _harness(data)
    before = state.state_hash()
    calls = 0

    def fail_after_first(operation):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise RuntimeError("injected mid-rewrite failure")

    monkeypatch.setattr(state, "_after_operation", fail_after_first)
    with pytest.raises(RuntimeError, match="injected"):
        state.apply(record, sites)
    assert calls == 2
    assert state.state_hash() == before


def _r1_pair(artifact):
    records = {record["event_id"]: record for record in artifact["records"]}
    capture = next(
        record
        for record in artifact["records"]
        if record["inventory_class"] == "R1:J_ring"
        and record["junction_ops"][0]["action"] == "create"
    )
    return capture, records[capture["reverse_of"]]


def _scission_record(artifact):
    return next(
        record
        for record in artifact["records"]
        if record["family"] == "R_Addition_MultipleBond"
        and record["arity"] == 1
        and any(
            operation.get("action") == "break" and operation.get("atoms") == [46, 47]
            for operation in record["bond_ops"]
        )
    )


def _reverse_site(state, reverse, entry, *, swap_labels=False):
    graph = _parse_adjacency(reverse["reactant_graphs"][0])
    expected = nx.Graph()
    for index, node in graph.items():
        expected.add_node(
            index,
            element=node["element"],
            radical=node.get("radical", 0),
            charge=node.get("charge", 0),
            lone_pairs=node.get("lone_pairs", 0),
            implicit_hydrogens=node.get("implicit_hydrogens", 0),
        )
        for neighbor, order in node.get("edges", {}).items():
            expected.add_edge(index, neighbor, order=float(order))
    actual = nx.compose_all(
        [_nx_executor(component) for component in entry.product_graphs]
    )
    node_match = nx.algorithms.isomorphism.categorical_node_match(
        ["element", "radical", "charge", "lone_pairs", "implicit_hydrogens"],
        [None, 0, 0, 0, 0],
    )
    edge_match = nx.algorithms.isomorphism.numerical_edge_match("order", 1.0)
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        expected, actual, node_match=node_match, edge_match=edge_match
    )
    break_atoms = next(
        op["atoms"] for op in reverse["bond_ops"] if op["action"] == "break"
    )
    bindings = dict(entry.bindings)
    assigned = next(
        mapping
        for mapping in matcher.isomorphisms_iter()
        if mapping[break_atoms[0]] == bindings["S9"]
        and mapping[break_atoms[1]] == bindings["P9"]
    )
    for index, node in sorted(graph.items()):
        atom_uuid = assigned[index]
        ref = state._ref(atom_uuid)
        node.update({"atom_ref": ref, "position": ref.position, "role": ref.role})
    first_label, second_label = ("P9", "S9") if swap_labels else ("S9", "P9")
    graph[break_atoms[0]]["label"] = first_label
    graph[break_atoms[1]]["label"] = second_label
    strand_id = state._ref(assigned[break_atoms[0]]).strand_id
    return Site("J_ring", strand_id, 0, 0, graph)


def test_real_j_para_round_trip_and_wrong_binding_refusal(ps_artifact):
    capture_data, reverse_data = _r1_pair(ps_artifact)
    form_atoms = next(
        op["atoms"] for op in capture_data["bond_ops"] if op["action"] == "form"
    )
    roles = {form_atoms[0]: "attacked_ring", form_atoms[1]: "attacker"}
    labels = {form_atoms[0]: "S9", form_atoms[1]: "P9"}
    state, capture, sites = _harness(capture_data, roles=roles, labels=labels)
    before = state.state_hash()
    capture_entry = state.apply(capture, sites)
    assert capture_entry.formed_bond in state.junctions
    reverse = EventRecord.from_dict(reverse_data)
    wrong = _reverse_site(state, reverse_data, capture_entry, swap_labels=True)
    captured_hash = state.state_hash()
    with pytest.raises(ReverseBindingError):
        state.apply(reverse, [wrong])
    assert state.state_hash() == captured_hash
    correct = _reverse_site(state, reverse_data, capture_entry)
    state.apply(reverse, [correct])
    assert state.state_hash() == before
    assert [event.event_id for event in state.event_log] == [
        capture.event_id,
        reverse.event_id,
    ]


def test_r1_scoped_restore_preserves_unrelated_later_sink_edit(ps_artifact):
    capture_data, reverse_data = _r1_pair(ps_artifact)
    form_atoms = next(
        op["atoms"] for op in capture_data["bond_ops"] if op["action"] == "form"
    )
    state, capture, sites = _harness(
        capture_data,
        roles={form_atoms[0]: "attacked_ring", form_atoms[1]: "attacker"},
        labels={form_atoms[0]: "S9", form_atoms[1]: "P9"},
    )
    spectator_ref = AtomRef("spectator-h", "spectator", 0, "feature")
    state.strands["spectator"] = Strand(
        "spectator",
        1,
        {},
        {"H": 1},
        0,
        {0: spectator_ref},
        {
            spectator_ref.uuid: {
                "element": "H",
                "radical": 0,
                "charge": 0,
                "lone_pairs": 0,
                "implicit_hydrogens": 0,
                "edges": {},
            }
        },
    )
    state._issued_uuids.add(spectator_ref.uuid)
    state._refresh_components({"spectator"})
    capture_entry = state.apply(capture, sites)

    sink_record = EventRecord(
        event_id="evt_sink_spectator",
        arity=1,
        participant_site_types=["spectator"],
        reactant_multiplicities=[1],
        reactant_graphs=["1 H u0 p0 c0"],
        product_graphs=["", "1 H u0 p0 c0"],
        coproducts=[{"formula": {"H": 1}}],
    )
    sink_graph = {
        0: {
            "element": "H",
            "radical": 0,
            "charge": 0,
            "lone_pairs": 0,
            "implicit_hydrogens": 0,
            "edges": {},
            "atom_ref": spectator_ref,
            "position": 0,
            "role": "feature",
        }
    }
    state.apply(sink_record, [Site("spectator", "spectator", 0, 0, sink_graph)])
    assert state.sink == {"H": 1}

    reverse = EventRecord.from_dict(reverse_data)
    reverse_site = _reverse_site(state, reverse_data, capture_entry)
    state.apply(reverse, [reverse_site])
    assert state.sink == {"H": 1}


def test_r1_reverse_refuses_an_intervening_edit_to_captured_component(ps_artifact):
    capture_data, reverse_data = _r1_pair(ps_artifact)
    form_atoms = next(
        op["atoms"] for op in capture_data["bond_ops"] if op["action"] == "form"
    )
    state, capture, sites = _harness(
        capture_data,
        roles={form_atoms[0]: "attacked_ring", form_atoms[1]: "attacker"},
        labels={form_atoms[0]: "S9", form_atoms[1]: "P9"},
    )
    capture_entry = state.apply(capture, sites)
    attacked = state._ref(capture_entry.bindings["S9"])
    edit = EventRecord(
        event_id="evt_intervening_formula",
        arity=1,
        participant_site_types=["captured"],
        reactant_multiplicities=[1],
        reactant_graphs=["1 C u0 p0 c0"],
        product_graphs=["1 C u0 p0 c0"],
        formula_delta={"H": 1},
    )
    graph = {
        0: {
            "element": "C",
            "radical": 0,
            "charge": 0,
            "lone_pairs": 0,
            "implicit_hydrogens": 0,
            "edges": {},
            "atom_ref": attacked,
            "position": attacked.position,
            "role": attacked.role,
        }
    }
    state.apply(
        edit, [Site("captured", attacked.strand_id, attacked.position, 0, graph)]
    )
    reverse = EventRecord.from_dict(reverse_data)
    reverse_site = _reverse_site(state, reverse_data, capture_entry)
    before = state.state_hash()
    with pytest.raises(ReverseBindingError, match="changed after"):
        state.apply(reverse, [reverse_site])
    assert state.state_hash() == before


def test_real_backbone_break_splits_features_and_preserves_uuids(ps_artifact):
    data = _scission_record(ps_artifact)
    break_op = next(
        operation
        for operation in data["bond_ops"]
        if operation["action"] == "break" and operation["atoms"] == [46, 47]
    )
    roles = {
        break_op["atoms"][0]: "backbone",
        break_op["atoms"][1]: "backbone",
        (break_op["atoms"][0], "position"): 3,
        (break_op["atoms"][1], "position"): 4,
    }
    state, record, sites = _harness(data, roles=roles)
    before_uuids = set(state._uuid_owners())
    entry = state.apply(record, sites)
    assert entry.derived_cut_offsets == (3,)
    assert len(state.strands) == 2
    assert set(state._uuid_owners()) == before_uuids
    state.assert_uuid_uniqueness()


def test_terminal_backbone_form_joins_two_sparse_strands():
    record = EventRecord(
        event_id="evt_terminal_join",
        arity=2,
        participant_site_types=["left", "right"],
        reactant_multiplicities=[1, 1],
        atom_map={0: 0, 1: 1},
        reactant_graphs=["1 C u0 p0 c0", "1 C u0 p0 c0"],
        product_graphs=["1 C u0 p0 c0 {2,S}\n2 C u0 p0 c0 {1,S}"],
        bond_ops=[{"action": "form", "atoms": [0, 1], "order": "1.0"}],
    )
    state, record, sites = _harness(
        record,
        roles={
            0: "backbone",
            1: "backbone",
            (0, "position"): 7,
            (1, "position"): 0,
        },
    )
    state.apply(record, sites)
    assert len(state.strands) == 1
    joined = next(iter(state.strands.values()))
    assert joined.length == 16
    assert not state.junctions


def test_negative_controls_inject_executor_defects(ps_artifact, monkeypatch):
    data = _scission_record(ps_artifact)
    roles = {
        46: "backbone",
        47: "backbone",
        (46, "position"): 3,
        (47, "position"): 4,
    }

    dropped, record, sites = _harness(data, roles=roles)
    original_operations = dropped._record_operations
    monkeypatch.setattr(
        dropped,
        "_record_operations",
        lambda item: [
            op
            for op in original_operations(item)
            if not (op["action"] == "break" and op["atoms"] == [46, 47])
        ],
    )
    with pytest.raises(LedgerError, match="share a mapped component"):
        dropped.apply(record, sites)

    capture_data, reverse_data = _r1_pair(ps_artifact)
    form_atoms = next(
        operation["atoms"]
        for operation in capture_data["bond_ops"]
        if operation["action"] == "form"
    )
    capture_roles = {form_atoms[0]: "attacked_ring", form_atoms[1]: "attacker"}
    capture_labels = {form_atoms[0]: "S9", form_atoms[1]: "P9"}
    skipped, capture, capture_sites = _harness(
        capture_data, roles=capture_roles, labels=capture_labels
    )
    capture_entry = skipped.apply(capture, capture_sites)
    reverse_site = _reverse_site(skipped, reverse_data, capture_entry)
    reverse = EventRecord.from_dict(reverse_data)
    monkeypatch.setattr(skipped, "_refresh_components", lambda affected=None: None)
    skipped.apply(reverse, [reverse_site])
    with pytest.raises(AssertionError, match="component ids"):
        skipped.assert_component_consistency()

    reused, record, sites = _harness(data, roles=roles)
    original_move = reused._moved_ref
    duplicate_uuid = next(iter(reused._uuid_owners()))

    def reuse_uuid(ref, strand_id, position):
        moved = original_move(ref, strand_id, position)
        if ref.uuid != duplicate_uuid:
            return AtomRef(duplicate_uuid, moved.strand_id, moved.position, moved.role)
        return moved

    monkeypatch.setattr(reused, "_moved_ref", reuse_uuid)
    with pytest.raises(UUIDReuseError):
        reused.apply(record, sites)


def test_precondition_failure_is_typed_and_atomic(ps_artifact):
    data = ps_artifact["records"][0]
    state, record, sites = _harness(data)
    before = state.state_hash()
    with pytest.raises(ArityError):
        state.apply(record, [])
    assert state.state_hash() == before
    bad = copy.deepcopy(sites[0].graph)
    source = next(index for index, node in bad.items() if node["edges"])
    neighbor = next(iter(bad[source]["edges"]))
    del bad[source]["edges"][neighbor]
    del bad[neighbor]["edges"][source]
    bad_site = Site(sites[0].site_type, sites[0].strand_id, 0, 0, bad)
    with pytest.raises(GraphMismatchError):
        state.apply(record, [bad_site] + sites[1:])
    assert state.state_hash() == before

    outsider = AtomRef("outsider-atom", "outsider", 0, "feature")
    element = sites[0].graph[0]["element"]
    state.strands["outsider"] = Strand(
        "outsider", 1, {}, {element: 1}, 0, {0: outsider}
    )
    state._issued_uuids.add(outsider.uuid)
    state._refresh_components({"outsider"})
    forged = copy.deepcopy(sites[0].graph)
    forged[0]["atom_ref"] = outsider
    forged_site = Site(sites[0].site_type, sites[0].strand_id, 0, 0, forged)
    forged_before = state.state_hash()
    with pytest.raises(GraphMismatchError, match="outside the participant component"):
        state.apply(record, [forged_site] + sites[1:])
    assert state.state_hash() == forged_before
