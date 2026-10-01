"""Sparse strand/component state and atomic compiled-rewrite executor.

F2 deliberately derives its execution plan rather than extending F1. At an
application site the record reactant graph is matched to the transient local
graph with the site atom anchored. Ties use ``atom_map`` and then mapped
``(position, role)``. Backbone cuts and fragment inheritance are derived from
mapped bond operations, and R1 labels bind through that same mapping to
persistent UUIDs. F1 therefore needs no strand positions, cut offsets, or
inheritance fields.

Strands store only length, sorted sparse features, and provenance atoms. A
transaction snapshots participant components and their incident junctions,
never the whole state.
"""

from __future__ import annotations

import copy
import hashlib
import json
import uuid
from collections import Counter
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Mapping, Sequence


ATOMIC_MASSES = {"C": 12, "H": 1, "N": 14, "O": 16, "S": 32}
ORDERS = {"S": "1.0", "D": "2.0", "T": "3.0", "B": "1.5"}


class StateError(RuntimeError):
    """Base class for state-executor failures."""


class ArityError(StateError):
    pass


class SiteTypeError(StateError):
    pass


class GraphMismatchError(StateError):
    pass


class LedgerError(StateError):
    pass


class ReverseBindingError(StateError):
    pass


class UUIDReuseError(StateError):
    pass


@dataclass(frozen=True)
class AtomRef:
    uuid: str
    strand_id: str
    position: int
    role: str = "feature"


@dataclass
class Strand:
    """A backbone length plus sorted sparse position-to-feature data."""

    strand_id: str
    length: int
    features: dict[int, dict[str, Any]] = field(default_factory=dict)
    formula: dict[str, int] = field(default_factory=dict)
    radical_count: int = 0
    atom_refs: dict[int, AtomRef] = field(default_factory=dict)
    atom_graph: dict[str, dict[str, Any]] = field(default_factory=dict)
    repeat_unit_formula: dict[str, int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.length < 1:
            raise ValueError("a strand must have a positive backbone length")
        if any(position < 0 or position >= self.length for position in self.features):
            raise ValueError("feature position lies outside its strand")
        self.features = dict(sorted(self.features.items()))
        self.formula = _clean_formula(self.formula)
        self.atom_graph = copy.deepcopy(self.atom_graph)
        self.repeat_unit_formula = _clean_formula(self.repeat_unit_formula)
        if not set(self.atom_graph) <= {ref.uuid for ref in self.atom_refs.values()}:
            raise ValueError("atom graph references an atom outside the strand")

    @property
    def mass(self) -> int:
        return _formula_mass(self.formula)


@dataclass(frozen=True)
class Site:
    site_type: str
    strand_id: str
    position: int
    atom: int
    graph: Mapping[int, Mapping[str, Any]]


@dataclass(frozen=True)
class AppliedEvent:
    event_id: str
    participants: tuple[str, ...]
    bindings: Mapping[str, str]
    derived_cut_offsets: tuple[int, ...]
    inheritance: Mapping[str, Any]
    formed_bond: tuple[str, str] | None = None
    product_graphs: tuple[tuple[tuple[Any, ...], ...], ...] = ()


@dataclass
class _Transaction:
    strands: dict[str, Strand]
    incident_junctions: set[tuple[str, str]]
    affected_uuids: set[str]
    sink: dict[str, int]
    sink_radicals: int
    event_log_length: int
    open_captures: dict[str, list[dict[str, Any]] | None]
    strand_serial: int


def _clean_formula(formula: Mapping[str, int]) -> dict[str, int]:
    return dict(
        sorted((key, int(value)) for key, value in formula.items() if int(value))
    )


def _formula_add(left: Mapping[str, int], right: Mapping[str, int]) -> dict[str, int]:
    result = Counter({key: int(value) for key, value in left.items()})
    result.update({key: int(value) for key, value in right.items()})
    return _clean_formula(result)


def _formula_subtract(
    left: Mapping[str, int], right: Mapping[str, int]
) -> dict[str, int]:
    result = Counter({key: int(value) for key, value in left.items()})
    result.subtract({key: int(value) for key, value in right.items()})
    if any(value < 0 for value in result.values()):
        raise LedgerError("coproduct debit exceeds the owning strand ledger")
    return _clean_formula(result)


def _formula_mass(formula: Mapping[str, int]) -> int:
    return sum(ATOMIC_MASSES.get(key, 0) * int(value) for key, value in formula.items())


def _order(value: Any) -> str:
    text = str(value)
    if text in ORDERS:
        return ORDERS[text]
    try:
        return f"{float(text):.1f}"
    except ValueError:
        return text


def _parse_adjacency(text: str) -> dict[int, dict[str, Any]]:
    """Parse the explicit-atom adjacency subset used by compiled records."""
    graph: dict[int, dict[str, Any]] = {}
    pending: list[tuple[int, int, str]] = []
    for line in text.splitlines():
        fields = line.split()
        if not fields or not fields[0].isdigit():
            continue
        index = int(fields[0]) - 1
        cursor = 1
        label = None
        if fields[cursor].startswith("*"):
            label = fields[cursor].lstrip("*")
            cursor += 1
        node = {
            "element": fields[cursor],
            "radical": 0,
            "charge": 0,
            "lone_pairs": 0,
            "implicit_hydrogens": 0,
            "edges": {},
        }
        if label:
            node["label"] = label
        for token in fields[cursor + 1 :]:
            if token.startswith("u") and token[1:].lstrip("-").isdigit():
                node["radical"] = int(token[1:])
            elif token.startswith("p") and token[1:].lstrip("-").isdigit():
                node["lone_pairs"] = int(token[1:])
            elif token.startswith("c") and token[1:].lstrip("+-").isdigit():
                node["charge"] = int(token[1:])
            elif token.startswith("{") and token.endswith("}"):
                other, bond_order = token[1:-1].split(",", 1)
                pending.append((index, int(other) - 1, _order(bond_order)))
        graph[index] = node
    for index, other, bond_order in pending:
        if other in graph:
            graph[index]["edges"][other] = bond_order
            graph[other]["edges"][index] = bond_order
    return graph


def _record_graphs(
    record: Any, field: str = "reactant_graphs"
) -> list[dict[int, dict[str, Any]]]:
    return [_parse_adjacency(item) for item in _value(record, field, ())]


def _value(record: Any, name: str, default: Any = None) -> Any:
    return (
        record.get(name, default)
        if isinstance(record, Mapping)
        else getattr(record, name, default)
    )


class KMCState:
    """Sparse connected-component state with all-or-nothing exact rewrites."""

    def __init__(self, strands: Sequence[Strand] = ()):
        self.strands = {strand.strand_id: copy.deepcopy(strand) for strand in strands}
        self.junctions: set[tuple[str, str]] = set()  # UUID endpoint pairs
        self.components: dict[str, str] = {}
        self.component_ledgers: dict[str, dict[str, int]] = {}
        self.component_radicals: dict[str, int] = {}
        self.component_lengths: dict[str, int] = {}
        self.component_moments: dict[str, tuple[int, int, int]] = {}
        self.moments: tuple[int, int, int] = (0, 0, 0)
        self.sink: dict[str, int] = {}
        self.sink_radicals = 0
        self.event_log: list[AppliedEvent] = []
        self._open_captures: dict[str, list[dict[str, Any]]] = {}
        self._issued_uuids = {
            ref.uuid
            for strand in self.strands.values()
            for ref in strand.atom_refs.values()
        }
        self._strand_serial = 0
        self.assert_uuid_uniqueness()
        self._refresh_components()

    def _allocate_uuid(self) -> str:
        value = str(uuid.uuid4())
        while value in self._issued_uuids:
            value = str(uuid.uuid4())
        self._issued_uuids.add(value)
        return value

    def new_atom(self, strand_id: str, position: int, role: str = "feature") -> AtomRef:
        ref = AtomRef(self._allocate_uuid(), strand_id, position, role)
        strand = self.strands[strand_id]
        strand.atom_refs[max(strand.atom_refs, default=-1) + 1] = ref
        return ref

    @property
    def total_ledger(self) -> dict[str, int]:
        result: dict[str, int] = {}
        for ledger in self.component_ledgers.values():
            result = _formula_add(result, ledger)
        return _formula_add(result, self.sink)

    @property
    def total_radicals(self) -> int:
        return sum(self.component_radicals.values()) + self.sink_radicals

    @property
    def total_mass(self) -> float:
        return float(_formula_mass(self.total_ledger))

    def state_hash(self) -> str:
        payload = {
            "strands": {
                key: {
                    "length": strand.length,
                    "features": strand.features,
                    "formula": strand.formula,
                    "radicals": strand.radical_count,
                    "atoms": {
                        i: ref.__dict__ for i, ref in sorted(strand.atom_refs.items())
                    },
                    "atom_graph": strand.atom_graph,
                    "repeat_unit_formula": strand.repeat_unit_formula,
                }
                for key, strand in sorted(self.strands.items())
            },
            "junctions": sorted(self.junctions),
            "sink": self.sink,
            "sink_radicals": self.sink_radicals,
        }
        encoded = json.dumps(
            payload, sort_keys=True, separators=(",", ":"), default=str
        )
        return hashlib.sha256(encoded.encode()).hexdigest()

    def _uuid_owners(self) -> dict[str, str]:
        return {
            ref.uuid: strand_id
            for strand_id, strand in self.strands.items()
            for ref in strand.atom_refs.values()
        }

    def _ref(self, atom_uuid: str) -> AtomRef:
        for strand in self.strands.values():
            for ref in strand.atom_refs.values():
                if ref.uuid == atom_uuid:
                    return ref
        raise GraphMismatchError(f"persistent atom {atom_uuid} is absent")

    def _component_partition(self) -> dict[str, str]:
        owners = self._uuid_owners()
        adjacency = {strand_id: set() for strand_id in self.strands}
        for left_uuid, right_uuid in self.junctions:
            if left_uuid not in owners or right_uuid not in owners:
                raise GraphMismatchError(
                    "junction references an absent persistent atom"
                )
            left, right = owners[left_uuid], owners[right_uuid]
            if left != right:
                adjacency[left].add(right)
                adjacency[right].add(left)
        result: dict[str, str] = {}
        for root in sorted(adjacency):
            if root in result:
                continue
            members: set[str] = set()
            stack = [root]
            while stack:
                current = stack.pop()
                if current in members:
                    continue
                members.add(current)
                stack.extend(adjacency[current] - members)
            component_id = min(members)
            result.update({member: component_id for member in members})
        return result

    def _refresh_components(self, affected: set[str] | None = None) -> None:
        if affected is None or not self.components:
            partition = self._component_partition()
            old_ids = set(self.component_ledgers)
            impacted = set(self.strands)
        else:
            old_ids = {self.components[s] for s in affected if s in self.components}
            impacted = {s for s, cid in self.components.items() if cid in old_ids}
            impacted.update(s for s in affected if s in self.strands)
            owners = self._uuid_owners()
            adjacency: dict[str, set[str]] = {}
            for left_uuid, right_uuid in self.junctions:
                left, right = owners[left_uuid], owners[right_uuid]
                if left != right:
                    adjacency.setdefault(left, set()).add(right)
                    adjacency.setdefault(right, set()).add(left)
            stack = list(impacted)
            while stack:
                current = stack.pop()
                for neighbor in adjacency.get(current, set()):
                    if neighbor not in impacted:
                        impacted.add(neighbor)
                        stack.append(neighbor)
            partition = {}
            unseen = impacted & set(self.strands)
            while unseen:
                root = min(unseen)
                members, stack = set(), [root]
                while stack:
                    current = stack.pop()
                    if current in members or current not in unseen:
                        continue
                    members.add(current)
                    stack.extend(adjacency.get(current, set()))
                unseen -= members
                component_id = min(members)
                partition.update({member: component_id for member in members})
            old_ids.update(self.components[s] for s in impacted if s in self.components)
        for component_id in old_ids:
            row = self.component_moments.pop(component_id, None)
            self.component_ledgers.pop(component_id, None)
            self.component_radicals.pop(component_id, None)
            self.component_lengths.pop(component_id, None)
            if row:
                self.moments = tuple(
                    left - right for left, right in zip(self.moments, row)
                )
        for strand_id in list(self.components):
            if strand_id not in self.strands or strand_id in impacted:
                self.components.pop(strand_id, None)
        for component_id in sorted({partition[s] for s in impacted if s in partition}):
            members = [s for s, cid in partition.items() if cid == component_id]
            ledger: dict[str, int] = {}
            radicals = length = 0
            for strand_id in members:
                strand = self.strands[strand_id]
                self.components[strand_id] = component_id
                ledger = _formula_add(ledger, strand.formula)
                radicals += strand.radical_count
                length += strand.length
            mass = _formula_mass(ledger)
            row = (1, mass, mass * mass)
            self.component_ledgers[component_id] = ledger
            self.component_radicals[component_id] = radicals
            self.component_lengths[component_id] = length
            self.component_moments[component_id] = row
            self.moments = tuple(left + right for left, right in zip(self.moments, row))

    def recompute_components(self) -> dict[str, str]:
        return self._component_partition()

    def brute_force_moments(self) -> tuple[int, int, int]:
        partition = self._component_partition()
        masses = [
            sum(self.strands[s].mass for s, value in partition.items() if value == cid)
            for cid in sorted(set(partition.values()))
        ]
        return len(masses), sum(masses), sum(mass * mass for mass in masses)

    def assert_component_consistency(self) -> None:
        if self.components != self.recompute_components():
            raise AssertionError("incremental component ids differ from full recompute")

    def assert_uuid_uniqueness(self) -> None:
        values = [
            ref.uuid
            for strand in self.strands.values()
            for ref in strand.atom_refs.values()
        ]
        if len(values) != len(set(values)):
            raise UUIDReuseError("persistent UUID reuse")
        if hasattr(self, "_issued_uuids") and not set(values) <= self._issued_uuids:
            raise UUIDReuseError("persistent UUID bypassed the allocator")

    def _match(
        self,
        expected: Mapping[int, Mapping[str, Any]],
        actual: Mapping[int, Mapping[str, Any]],
        anchor: int,
        ranks: Mapping[int, int],
    ) -> dict[int, int]:
        if anchor not in expected or anchor not in actual:
            raise GraphMismatchError("missing anchored atom")
        if len(expected) != len(actual):
            raise GraphMismatchError("local graph has the wrong atom count")

        def compatible(left: Mapping[str, Any], right: Mapping[str, Any]) -> bool:
            return (
                left.get("element") == right.get("element")
                and int(left.get("radical", 0)) == int(right.get("radical", 0))
                and int(left.get("charge", 0)) == int(right.get("charge", 0))
                and int(left.get("lone_pairs", 0)) == int(right.get("lone_pairs", 0))
                and int(left.get("implicit_hydrogens", 0))
                == int(right.get("implicit_hydrogens", 0))
            )

        # Normal site indexing preserves canonical record order. Validate the
        # fast path fully before the general isomorphism fallback.
        if set(expected) == set(actual) and all(
            compatible(expected[i], actual[i]) for i in expected
        ):
            if all(
                {n: _order(o) for n, o in expected[i].get("edges", {}).items()}
                == {n: _order(o) for n, o in actual[i].get("edges", {}).items()}
                for i in expected
            ):
                return {i: i for i in expected}

        candidates = {
            source: [target for target, node in actual.items() if compatible(exp, node)]
            for source, exp in expected.items()
        }
        if anchor not in candidates[anchor]:
            raise GraphMismatchError(
                "anchored atom does not match record reactant graph"
            )
        candidates[anchor] = [anchor]
        visit_order = sorted(
            expected, key=lambda i: (i != anchor, len(candidates[i]), ranks.get(i, i))
        )
        found: list[dict[int, int]] = []

        def visit(offset: int, mapping: dict[int, int]) -> None:
            if offset == len(visit_order):
                if all(
                    {
                        mapping[neighbor]: _order(edge)
                        for neighbor, edge in expected[source].get("edges", {}).items()
                    }
                    == {
                        neighbor: _order(edge)
                        for neighbor, edge in actual[mapping[source]]
                        .get("edges", {})
                        .items()
                    }
                    for source in expected
                ):
                    found.append(dict(mapping))
                return
            source = visit_order[offset]
            for target in candidates[source]:
                if target in mapping.values():
                    continue
                if all(
                    _order(actual[target].get("edges", {}).get(mapping[neighbor]))
                    == _order(edge)
                    for neighbor, edge in expected[source].get("edges", {}).items()
                    if neighbor in mapping
                ):
                    mapping[source] = target
                    visit(offset + 1, mapping)
                    del mapping[source]

        visit(0, {})
        if not found:
            raise GraphMismatchError("local graph does not match record reactant graph")
        ranked = sorted(expected, key=lambda i: (ranks.get(i, i), i))
        return min(
            found,
            key=lambda mapping: tuple(
                (
                    actual[mapping[i]].get("position", 0),
                    actual[mapping[i]].get("role", ""),
                    mapping[i],
                )
                for i in ranked
            ),
        )

    def _expected_site_types(self, record: Any) -> list[str]:
        types = list(_value(record, "participant_site_types", []))
        multiplicities = list(_value(record, "reactant_multiplicities", []))
        if multiplicities and len(types) == len(multiplicities):
            return [
                kind
                for kind, count in zip(types, multiplicities)
                for _ in range(int(count))
            ]
        return types

    def _preflight(self, record: Any, participants: Sequence[Site]) -> tuple[
        dict[int, str],
        dict[str, str],
        dict[str, dict[str, Any]],
        list[dict[int, dict[str, Any]]],
    ]:
        if len(participants) != int(_value(record, "arity", 0)):
            raise ArityError("participant arity does not match record")
        if [
            participant.site_type for participant in participants
        ] != self._expected_site_types(record):
            raise SiteTypeError("participant site types do not match record")
        graphs = _record_graphs(record)
        if len(graphs) != len(participants):
            raise GraphMismatchError("record graph count does not match participants")

        raw_atom_map = _value(record, "atom_map", {})
        atom_map = {int(key): int(value) for key, value in raw_atom_map.items()}
        heavy_rank: dict[int, int] = {}
        heavy = global_offset = 0
        for graph in graphs:
            for local, node in sorted(graph.items()):
                if node.get("element") != "H":
                    heavy_rank[global_offset + local] = atom_map.get(heavy, heavy)
                    heavy += 1
            global_offset += len(graph)

        mapping: dict[int, str] = {}
        bindings: dict[str, str] = {}
        local_graph: dict[str, dict[str, Any]] = {}
        offset = 0
        for expected, participant in zip(graphs, participants):
            if participant.strand_id not in self.strands:
                raise GraphMismatchError("participant strand is absent")
            ranks = {i: heavy_rank.get(offset + i, offset + i) for i in expected}
            current_refs: dict[int, AtomRef] = {}
            for local, supplied in participant.graph.items():
                supplied_ref = supplied.get("atom_ref")
                if not isinstance(supplied_ref, AtomRef):
                    raise GraphMismatchError("local graph atom lacks persistent UUID")
                ref = self._ref(supplied_ref.uuid)
                if (
                    self.components[ref.strand_id]
                    != self.components[participant.strand_id]
                ):
                    raise GraphMismatchError(
                        "local graph atom lies outside the participant component"
                    )
                current_refs[local] = ref
            uuid_to_local = {ref.uuid: local for local, ref in current_refs.items()}
            if len(uuid_to_local) != len(current_refs):
                raise GraphMismatchError("local graph repeats a persistent UUID")
            current_graph: dict[int, dict[str, Any]] = {}
            for local, ref in current_refs.items():
                state_node = self.strands[ref.strand_id].atom_graph.get(ref.uuid)
                if state_node is None:
                    raise GraphMismatchError(
                        "persistent atom has no state-owned local graph data"
                    )
                current_edges = {
                    uuid_to_local[neighbor]: _order(bond_order)
                    for neighbor, bond_order in state_node.get("edges", {}).items()
                    if neighbor in uuid_to_local
                }
                supplied_edges = {
                    neighbor: _order(bond_order)
                    for neighbor, bond_order in participant.graph[local]
                    .get("edges", {})
                    .items()
                    if neighbor in current_refs
                }
                if supplied_edges != current_edges:
                    raise GraphMismatchError(
                        "site-index graph disagrees with state-owned adjacency"
                    )
                current_graph[local] = {
                    "element": state_node.get("element"),
                    "radical": int(state_node.get("radical", 0)),
                    "charge": int(state_node.get("charge", 0)),
                    "lone_pairs": int(state_node.get("lone_pairs", 0)),
                    "implicit_hydrogens": int(state_node.get("implicit_hydrogens", 0)),
                    "edges": current_edges,
                    "position": ref.position,
                    "role": ref.role,
                }
                label = participant.graph[local].get("label")
                if isinstance(label, str):
                    current_graph[local]["label"] = label
            matched = self._match(expected, current_graph, participant.atom, ranks)
            inverse_match = {actual: source for source, actual in matched.items()}
            for source, actual in matched.items():
                node = current_graph[actual]
                ref = current_refs[actual]
                mapping[offset + source] = ref.uuid
                local_graph[ref.uuid] = {
                    "element": node.get("element"),
                    "radical": int(node.get("radical", 0)),
                    "charge": int(node.get("charge", 0)),
                    "lone_pairs": int(node.get("lone_pairs", 0)),
                    "implicit_hydrogens": int(node.get("implicit_hydrogens", 0)),
                    "edges": {},
                }
                if isinstance(node.get("label"), str):
                    bindings[node["label"].lstrip("*")] = ref.uuid
            for actual, node in current_graph.items():
                if actual not in inverse_match:
                    continue
                source_uuid = mapping[offset + inverse_match[actual]]
                for neighbor, bond_order in node.get("edges", {}).items():
                    if neighbor in inverse_match:
                        local_graph[source_uuid]["edges"][
                            mapping[offset + inverse_match[neighbor]]
                        ] = _order(bond_order)
            offset += len(expected)
        self._derive_r1_bindings(record, mapping, bindings)
        return mapping, bindings, local_graph, graphs

    def _derive_r1_bindings(
        self, record: Any, mapping: Mapping[int, str], bindings: dict[str, str]
    ) -> None:
        for junction in _value(record, "junction_ops", None) or []:
            attacked = junction.get("attacked_atom_label")
            attacker = junction.get("attacker_site_label")
            if attacked in bindings and attacker in bindings:
                continue
            action = "form" if junction.get("action") == "create" else "break"
            candidates = [
                op
                for op in _value(record, "bond_ops", [])
                if op.get("action") == action and len(op.get("atoms", [])) == 2
            ]
            labels = junction.get("formed_bond", {}).get("label_pair", [])
            if len(candidates) != 1 or len(labels) != 2:
                continue
            first, second = (
                self._ref(mapping[index]) for index in candidates[0]["atoms"]
            )
            if "ring" in first.role or "attacked" in first.role:
                bindings[attacked], bindings[attacker] = first.uuid, second.uuid
            elif "ring" in second.role or "attacked" in second.role:
                bindings[attacked], bindings[attacker] = second.uuid, first.uuid
            else:
                bindings[labels[0]], bindings[labels[1]] = first.uuid, second.uuid

    def _affected_strands(self, participants: Sequence[Site]) -> set[str]:
        component_ids = {
            self.components[p.strand_id]
            for p in participants
            if p.strand_id in self.components
        }
        return {
            strand_id
            for strand_id, cid in self.components.items()
            if cid in component_ids
        }

    def _begin(self, participants: Sequence[Site], record: Any) -> _Transaction:
        affected = self._affected_strands(participants)
        atom_uuids = {
            ref.uuid
            for strand_id in affected
            for ref in self.strands[strand_id].atom_refs.values()
        }
        capture_keys = {
            str(value)
            for value in (_value(record, "event_id"), _value(record, "reverse_of"))
            if value
        }
        return _Transaction(
            {
                strand_id: copy.deepcopy(self.strands[strand_id])
                for strand_id in affected
            },
            {
                edge
                for edge in self.junctions
                if edge[0] in atom_uuids or edge[1] in atom_uuids
            },
            atom_uuids,
            dict(self.sink),
            self.sink_radicals,
            len(self.event_log),
            {key: copy.deepcopy(self._open_captures.get(key)) for key in capture_keys},
            self._strand_serial,
        )

    def _rollback(self, snapshot: _Transaction) -> None:
        current_keys = {
            strand_id
            for strand_id, strand in self.strands.items()
            if any(
                ref.uuid in snapshot.affected_uuids for ref in strand.atom_refs.values()
            )
        }
        for key in current_keys:
            del self.strands[key]
        self.strands.update(snapshot.strands)
        self.junctions = {
            edge
            for edge in self.junctions
            if edge[0] not in snapshot.affected_uuids
            and edge[1] not in snapshot.affected_uuids
        } | snapshot.incident_junctions
        self.sink = snapshot.sink
        self.sink_radicals = snapshot.sink_radicals
        del self.event_log[snapshot.event_log_length :]
        for key, captures in snapshot.open_captures.items():
            if captures is None:
                self._open_captures.pop(key, None)
            else:
                self._open_captures[key] = captures
        self._strand_serial = snapshot.strand_serial
        self._refresh_components()

    def _record_operations(self, record: Any) -> list[dict[str, Any]]:
        operations = [dict(operation) for operation in _value(record, "bond_ops", [])]
        seen = {json.dumps(operation, sort_keys=True) for operation in operations}
        for operation in _value(record, "feature_ops", None) or []:
            encoded = json.dumps(operation, sort_keys=True)
            if encoded not in seen:
                operations.append(dict(operation))
                seen.add(encoded)
        return operations

    def _feature_atom(self, atom_uuid: str, **changes: Any) -> None:
        ref = self._ref(atom_uuid)
        strand = self.strands[ref.strand_id]
        atoms = strand.features.setdefault(ref.position, {}).setdefault("atoms", {})
        atoms.setdefault(atom_uuid, {}).update(changes)
        state_node = strand.atom_graph.get(atom_uuid)
        if state_node is None:
            raise GraphMismatchError("persistent atom has no state-owned graph node")
        state_node.update(changes)
        strand.features = dict(sorted(strand.features.items()))

    def _feature_bond(self, first: str, second: str, bond_order: str | None) -> None:
        key = "|".join(sorted((first, second)))
        for atom_uuid in (first, second):
            ref = self._ref(atom_uuid)
            feature = self.strands[ref.strand_id].features.setdefault(ref.position, {})
            bonds = feature.setdefault("bonds", {})
            if bond_order is None:
                bonds.pop(key, None)
            else:
                bonds[key] = bond_order
            if not bonds:
                feature.pop("bonds", None)
        first_ref, second_ref = self._ref(first), self._ref(second)
        first_node = self.strands[first_ref.strand_id].atom_graph[first]
        second_node = self.strands[second_ref.strand_id].atom_graph[second]
        if bond_order is None:
            first_node["edges"].pop(second, None)
            second_node["edges"].pop(first, None)
        else:
            first_node["edges"][second] = bond_order
            second_node["edges"][first] = bond_order

    def _rewrite_graph(
        self,
        operations: Sequence[Mapping[str, Any]],
        mapping: Mapping[int, str],
        graph: dict[str, dict[str, Any]],
    ) -> tuple[list[tuple[str, str]], int]:
        cuts: list[tuple[str, str]] = []
        radical_change = 0
        replacement_pairs = {
            tuple(sorted(operation["atoms"]))
            for operation in operations
            if operation.get("action") in {"form", "change"}
            and len(operation.get("atoms", [])) == 2
        }
        for operation in operations:
            action = operation.get("action")
            if action in {"break", "form", "change"}:
                try:
                    first, second = (
                        mapping[int(index)] for index in operation["atoms"]
                    )
                except (KeyError, TypeError) as error:
                    raise GraphMismatchError(
                        "bond operation references an unmapped atom"
                    ) from error
                existing = graph[first]["edges"].get(second)
                if action == "break":
                    if existing is None:
                        raise GraphMismatchError("cannot break a missing mapped bond")
                    graph[first]["edges"].pop(second, None)
                    graph[second]["edges"].pop(first, None)
                    self._feature_bond(first, second, None)
                    first_ref, second_ref = self._ref(first), self._ref(second)
                    if (
                        tuple(sorted(operation["atoms"])) not in replacement_pairs
                        and first_ref.strand_id == second_ref.strand_id
                        and first_ref.role == second_ref.role == "backbone"
                        and abs(first_ref.position - second_ref.position) == 1
                    ):
                        cuts.append((first, second))
                    self.junctions.discard(tuple(sorted((first, second))))
                else:
                    if action == "form" and existing is not None:
                        raise GraphMismatchError("cannot form an existing mapped bond")
                    if action == "change" and existing is None:
                        raise GraphMismatchError("cannot change a missing mapped bond")
                    bond_order = _order(operation["order"])
                    graph[first]["edges"][second] = bond_order
                    graph[second]["edges"][first] = bond_order
                    self._feature_bond(first, second, bond_order)
                    first_ref, second_ref = self._ref(first), self._ref(second)
                    if (
                        action == "form"
                        and first_ref.strand_id != second_ref.strand_id
                        and graph[first].get("element") != "H"
                        and graph[second].get("element") != "H"
                    ):
                        if not self._join_terminal_strands(first, second):
                            self.junctions.add(tuple(sorted((first, second))))
            elif action.startswith("set_"):
                try:
                    atom_uuid = mapping[int(operation["atom"])]
                except (KeyError, TypeError) as error:
                    raise GraphMismatchError(
                        "feature operation references an unmapped atom"
                    ) from error
                field_name = {
                    "set_radical": "radical",
                    "set_charge": "charge",
                    "set_lone_pairs": "lone_pairs",
                    "set_implicit_hydrogens": "implicit_hydrogens",
                }.get(action)
                if field_name is None:
                    raise GraphMismatchError(f"unknown feature rewrite {action!r}")
                new_value = int(operation["value"])
                old_value = int(graph[atom_uuid].get(field_name, 0))
                graph[atom_uuid][field_name] = new_value
                self._feature_atom(atom_uuid, **{field_name: new_value})
                if field_name == "radical":
                    change = new_value - old_value
                    self.strands[self._ref(atom_uuid).strand_id].radical_count += change
                    radical_change += change
            else:
                raise GraphMismatchError(f"unknown rewrite action {action!r}")
            self._after_operation(operation)
        return cuts, radical_change

    def _after_operation(self, operation: Mapping[str, Any]) -> None:
        """Fault-injection seam used to prove mid-rewrite atomicity."""

    def _allocate_strand_id(self, base: str) -> str:
        while True:
            self._strand_serial += 1
            candidate = f"{base}~{self._strand_serial}"
            if candidate not in self.strands:
                return candidate

    def _orient_strand_endpoint(
        self, strand: Strand, endpoint_uuid: str, want_end: bool
    ) -> Strand:
        endpoint = next(
            ref for ref in strand.atom_refs.values() if ref.uuid == endpoint_uuid
        )
        wanted = strand.length - 1 if want_end else 0
        if endpoint.position == wanted:
            return copy.deepcopy(strand)
        opposite = 0 if want_end else strand.length - 1
        if endpoint.position != opposite:
            raise GraphMismatchError(
                "terminal join references an interior backbone atom"
            )
        features = {
            strand.length - 1 - position: copy.deepcopy(value)
            for position, value in strand.features.items()
        }
        refs = {
            key: self._moved_ref(
                ref, strand.strand_id, strand.length - 1 - ref.position
            )
            for key, ref in strand.atom_refs.items()
        }
        return Strand(
            strand.strand_id,
            strand.length,
            features,
            strand.formula,
            strand.radical_count,
            refs,
            strand.atom_graph,
            strand.repeat_unit_formula,
        )

    def _join_terminal_strands(self, first_uuid: str, second_uuid: str) -> bool:
        first_ref, second_ref = self._ref(first_uuid), self._ref(second_uuid)
        if first_ref.role != "backbone" or second_ref.role != "backbone":
            return False
        first_strand, second_strand = (
            self.strands[first_ref.strand_id],
            self.strands[second_ref.strand_id],
        )
        if first_ref.position not in {0, first_strand.length - 1}:
            return False
        if second_ref.position not in {0, second_strand.length - 1}:
            return False
        if (
            first_strand.repeat_unit_formula
            and second_strand.repeat_unit_formula
            and first_strand.repeat_unit_formula != second_strand.repeat_unit_formula
        ):
            return False
        left = self._orient_strand_endpoint(first_strand, first_uuid, True)
        right = self._orient_strand_endpoint(second_strand, second_uuid, False)
        shift = left.length
        features = copy.deepcopy(left.features)
        features.update(
            {
                position + shift: copy.deepcopy(value)
                for position, value in right.features.items()
            }
        )
        refs = dict(left.atom_refs)
        next_key = max(refs, default=-1) + 1
        for ref in right.atom_refs.values():
            refs[next_key] = self._moved_ref(ref, left.strand_id, ref.position + shift)
            next_key += 1
        atom_graph = copy.deepcopy(left.atom_graph)
        atom_graph.update(copy.deepcopy(right.atom_graph))
        self.strands[left.strand_id] = Strand(
            left.strand_id,
            left.length + right.length,
            features,
            _formula_add(left.formula, right.formula),
            left.radical_count + right.radical_count,
            refs,
            atom_graph,
            left.repeat_unit_formula or right.repeat_unit_formula,
        )
        del self.strands[right.strand_id]
        return True

    def _moved_ref(self, ref: AtomRef, strand_id: str, position: int) -> AtomRef:
        """Preserve identity while changing location (UUID mutation-test seam)."""
        return AtomRef(ref.uuid, strand_id, position, ref.role)

    def _split_strand(
        self,
        first_uuid: str,
        second_uuid: str,
        atom_elements: Mapping[str, str] | None = None,
    ) -> tuple[str, int]:
        first, second = self._ref(first_uuid), self._ref(second_uuid)
        if first.strand_id != second.strand_id:
            raise GraphMismatchError("backbone cut endpoints do not share a strand")
        strand = self.strands[first.strand_id]
        cut = min(first.position, second.position)
        if cut < 0 or cut >= strand.length - 1:
            raise GraphMismatchError("mapped backbone cut lies outside the strand")
        right_id = self._allocate_strand_id(strand.strand_id)
        left_length, right_length = cut + 1, strand.length - cut - 1
        left_features = {
            position: copy.deepcopy(value)
            for position, value in strand.features.items()
            if position <= cut
        }
        right_features = {
            position - left_length: copy.deepcopy(value)
            for position, value in strand.features.items()
            if position > cut
        }
        left_refs: dict[int, AtomRef] = {}
        right_refs: dict[int, AtomRef] = {}
        for key, ref in strand.atom_refs.items():
            if ref.position <= cut:
                left_refs[key] = self._moved_ref(ref, strand.strand_id, ref.position)
            else:
                right_refs[key] = self._moved_ref(
                    ref, right_id, ref.position - left_length
                )
        moved_uuids = [ref.uuid for ref in left_refs.values()] + [
            ref.uuid for ref in right_refs.values()
        ]
        if len(moved_uuids) != len(set(moved_uuids)):
            raise UUIDReuseError("persistent UUID reuse during strand split")
        atom_elements = atom_elements or {}
        left_known = Counter(
            atom_elements[ref.uuid]
            for ref in left_refs.values()
            if ref.uuid in atom_elements
        )
        right_known = Counter(
            atom_elements[ref.uuid]
            for ref in right_refs.values()
            if ref.uuid in atom_elements
        )
        left_formula: dict[str, int] = {}
        right_formula: dict[str, int] = {}
        for element, count in strand.formula.items():
            known = left_known[element] + right_known[element]
            if known > int(count):
                raise LedgerError("mapped atoms exceed the strand element ledger")
            residual = int(count) - known
            unit_count = strand.repeat_unit_formula.get(element, 0)
            if residual != unit_count * strand.length:
                raise LedgerError(
                    "exact split needs mapped atoms or an unmapped repeat-unit formula"
                )
            left_residual = unit_count * left_length
            left_formula[element] = left_known[element] + left_residual
            right_formula[element] = right_known[element] + residual - left_residual
        left_graph = {
            ref.uuid: copy.deepcopy(strand.atom_graph[ref.uuid])
            for ref in left_refs.values()
        }
        right_graph = {
            ref.uuid: copy.deepcopy(strand.atom_graph[ref.uuid])
            for ref in right_refs.values()
        }
        right_radicals = min(
            strand.radical_count,
            sum(
                int(feature.get("radical", False))
                + sum(
                    int(atom.get("radical", 0))
                    for atom in feature.get("atoms", {}).values()
                )
                for feature in right_features.values()
            ),
        )
        self.strands[strand.strand_id] = Strand(
            strand.strand_id,
            left_length,
            left_features,
            left_formula,
            strand.radical_count - right_radicals,
            left_refs,
            left_graph,
            strand.repeat_unit_formula,
        )
        self.strands[right_id] = Strand(
            right_id,
            right_length,
            right_features,
            right_formula,
            right_radicals,
            right_refs,
            right_graph,
            strand.repeat_unit_formula,
        )
        self.assert_uuid_uniqueness()
        return right_id, cut

    def _transfer_mobile_atoms(
        self, graph: Mapping[str, Mapping[str, Any]]
    ) -> set[str]:
        """Rehome transferred explicit H atoms on the strand owning the new bond."""
        affected: set[str] = set()
        for atom_uuid, node in graph.items():
            if node.get("element") != "H":
                continue
            heavy_neighbors = [
                neighbor
                for neighbor in node.get("edges", {})
                if graph[neighbor].get("element") != "H"
            ]
            if len(heavy_neighbors) != 1:
                continue
            old_ref = self._ref(atom_uuid)
            new_ref = self._ref(heavy_neighbors[0])
            if old_ref.strand_id == new_ref.strand_id:
                continue
            old_strand, new_strand = (
                self.strands[old_ref.strand_id],
                self.strands[new_ref.strand_id],
            )
            old_key = next(
                key
                for key, ref in old_strand.atom_refs.items()
                if ref.uuid == atom_uuid
            )
            del old_strand.atom_refs[old_key]
            new_key = max(new_strand.atom_refs, default=-1) + 1
            new_strand.atom_refs[new_key] = self._moved_ref(
                old_ref, new_ref.strand_id, new_ref.position
            )
            new_strand.atom_graph[atom_uuid] = old_strand.atom_graph.pop(atom_uuid)
            old_strand.formula = _formula_subtract(old_strand.formula, {"H": 1})
            new_strand.formula = _formula_add(new_strand.formula, {"H": 1})
            radical = int(node.get("radical", 0))
            old_strand.radical_count -= radical
            new_strand.radical_count += radical
            affected.update((old_ref.strand_id, new_ref.strand_id))
        return affected

    def _apply_formula_delta(
        self, delta: Mapping[str, int], owners: Sequence[str]
    ) -> None:
        distinct = sorted(set(owners))
        if not distinct:
            if any(int(value) for value in delta.values()):
                raise LedgerError("record delta has no mapped owning strand")
            return
        if len(distinct) > 1 and any(int(value) for value in delta.values()):
            raise LedgerError(
                "multi-participant formula delta lacks per-atom ownership"
            )
        for element, value in sorted(delta.items()):
            amount, sign = abs(int(value)), 1 if int(value) >= 0 else -1
            for index in range(amount):
                owner = distinct[0]
                before = self.strands[owner].formula
                if sign < 0:
                    self.strands[owner].formula = _formula_subtract(
                        before, {element: 1}
                    )
                else:
                    self.strands[owner].formula = _formula_add(before, {element: 1})

    def _debit_coproducts(
        self,
        record: Any,
        mapping: Mapping[int, str],
        reactant_graphs: Sequence[Mapping[int, Mapping[str, Any]]],
        local_graph: Mapping[str, Mapping[str, Any]],
    ) -> None:
        coproducts = list(_value(record, "coproducts", []))
        if not coproducts:
            return
        product_graphs = _record_graphs(record, "product_graphs")
        if len(product_graphs) != len(coproducts) + 1:
            raise LedgerError("coproducts do not align with product graphs")
        ranges: list[range] = []
        cursor = 0
        for graph in product_graphs:
            count = sum(node.get("element") != "H" for node in graph.values())
            ranges.append(range(cursor, cursor + count))
            cursor += count
        atom_map = {
            int(key): int(value)
            for key, value in _value(record, "atom_map", {}).items()
        }
        inverse_map = {product: reactant for reactant, product in atom_map.items()}
        heavy_globals: list[int] = []
        offset = 0
        for graph in reactant_graphs:
            heavy_globals.extend(
                offset + index
                for index, node in sorted(graph.items())
                if node.get("element") != "H"
            )
            offset += len(graph)

        components = self._graph_component_members(local_graph)
        component_for = {
            atom_uuid: index
            for index, members in enumerate(components)
            for atom_uuid in members
        }

        assignments: dict[int, int | None] = {}
        claimed: set[int] = set()
        deferred: list[int] = []
        for product_index, (product_graph, product_range) in enumerate(
            zip(product_graphs, ranges)
        ):
            anchors = []
            for product_atom in product_range:
                reactant_atom = inverse_map.get(product_atom)
                if reactant_atom is None or reactant_atom >= len(heavy_globals):
                    raise LedgerError("product atom map is incomplete")
                try:
                    anchors.append(mapping[heavy_globals[reactant_atom]])
                except KeyError as error:
                    raise LedgerError("product atom map is incomplete") from error
            if anchors:
                component_indexes = {
                    component_for.get(atom_uuid) for atom_uuid in anchors
                }
                if None in component_indexes or len(component_indexes) != 1:
                    raise LedgerError("mapped product atoms span product components")
                component_index = component_indexes.pop()
                if component_index in claimed:
                    raise LedgerError("product graphs share a mapped component")
                assignments[product_index] = component_index
                claimed.add(component_index)
            elif not product_graph:
                assignments[product_index] = None
            else:
                deferred.append(product_index)

        candidates = {
            product_index: [
                component_index
                for component_index, members in enumerate(components)
                if component_index not in claimed
                and self._component_matches_product(
                    product_graphs[product_index], local_graph, members
                )
            ]
            for product_index in deferred
        }
        solutions: list[dict[int, int]] = []

        def assign_deferred(
            position: int, used: set[int], current: dict[int, int]
        ) -> None:
            if len(solutions) > 1:
                return
            if position == len(deferred):
                solutions.append(dict(current))
                return
            product_index = deferred[position]
            for component_index in candidates[product_index]:
                if component_index in used:
                    continue
                current[product_index] = component_index
                assign_deferred(position + 1, used | {component_index}, current)
                del current[product_index]

        assign_deferred(0, set(claimed), {})
        if len(solutions) != 1:
            raise LedgerError("unmapped product components are ambiguous")
        assignments.update(solutions[0])
        assigned_components = {
            component_index
            for component_index in assignments.values()
            if component_index is not None
        }
        if assigned_components != set(range(len(components))):
            raise LedgerError("product component attribution is incomplete")

        for product_index, coproduct in enumerate(coproducts, start=1):
            component_index = assignments[product_index]
            if component_index is None:
                raise LedgerError("coproduct graph has no attributable atoms")
            members = components[component_index]
            formula = coproduct.get("formula", coproduct.get("element_delta", {}))
            formula_by_owner: dict[str, Counter] = {}
            radicals_by_owner: Counter = Counter()
            actual_formula: Counter = Counter()
            for atom_uuid in members:
                node = local_graph[atom_uuid]
                owner = self._coproduct_atom_owner(atom_uuid)
                formula_by_owner.setdefault(owner, Counter())[node.get("element")] += 1
                radical = int(node.get("radical", 0))
                radicals_by_owner[owner] += radical
                actual_formula[node.get("element")] += 1
            if _clean_formula(actual_formula) != _clean_formula(formula):
                raise LedgerError("coproduct graph formula disagrees with its ledger")
            for owner, owner_formula in formula_by_owner.items():
                self.strands[owner].formula = _formula_subtract(
                    self.strands[owner].formula, owner_formula
                )
            self.sink = _formula_add(self.sink, formula)
            radicals = sum(radicals_by_owner.values())
            for owner, owner_radicals in radicals_by_owner.items():
                if self.strands[owner].radical_count < owner_radicals:
                    raise LedgerError(
                        "coproduct radical debit exceeds its owning strand"
                    )
                self.strands[owner].radical_count -= owner_radicals
            self.sink_radicals += radicals
            for atom_uuid in members:
                if int(local_graph[atom_uuid].get("radical", 0)):
                    self._feature_atom(atom_uuid, radical=0)

    def _component_matches_product(
        self,
        expected: Mapping[int, Mapping[str, Any]],
        actual: Mapping[str, Mapping[str, Any]],
        members: set[str],
    ) -> bool:
        if len(expected) != len(members):
            return False

        def signature(node: Mapping[str, Any]) -> tuple[Any, ...]:
            return (
                node.get("element"),
                int(node.get("radical", 0)),
                int(node.get("charge", 0)),
                int(node.get("lone_pairs", 0)),
                int(node.get("implicit_hydrogens", 0)),
            )

        candidates = {
            source: [
                target
                for target in members
                if signature(expected_node) == signature(actual[target])
                and len(expected_node.get("edges", {}))
                == len(set(actual[target].get("edges", {})) & members)
            ]
            for source, expected_node in expected.items()
        }
        if any(not options for options in candidates.values()):
            return False
        order = sorted(expected, key=lambda source: len(candidates[source]))

        def visit(position: int, matched: dict[int, str]) -> bool:
            if position == len(order):
                return True
            source = order[position]
            for target in candidates[source]:
                if target in matched.values():
                    continue
                if all(
                    _order(actual[target].get("edges", {}).get(matched[neighbor]))
                    == _order(bond_order)
                    for neighbor, bond_order in expected[source]
                    .get("edges", {})
                    .items()
                    if neighbor in matched
                ):
                    matched[source] = target
                    if visit(position + 1, matched):
                        return True
                    del matched[source]
            return False

        return visit(0, {})

    def _coproduct_atom_owner(self, atom_uuid: str) -> str:
        return self._ref(atom_uuid).strand_id

    def _graph_component_members(
        self, graph: Mapping[str, Mapping[str, Any]]
    ) -> list[set[str]]:
        unseen = set(graph)
        components = []
        while unseen:
            root = min(unseen)
            stack, members = [root], set()
            while stack:
                current = stack.pop()
                if current in members:
                    continue
                members.add(current)
                stack.extend(set(graph[current]["edges"]) - members)
            unseen -= members
            components.append(members)
        return components

    def _graph_components(
        self, graph: Mapping[str, Mapping[str, Any]]
    ) -> tuple[tuple[tuple[Any, ...], ...], ...]:
        components = []
        for members in self._graph_component_members(graph):
            rows = []
            for atom_uuid in sorted(members):
                node = graph[atom_uuid]
                rows.append(
                    (
                        atom_uuid,
                        node.get("element"),
                        int(node.get("radical", 0)),
                        int(node.get("charge", 0)),
                        int(node.get("lone_pairs", 0)),
                        int(node.get("implicit_hydrogens", 0)),
                        tuple(
                            sorted(
                                (neighbor, _order(bond_order))
                                for neighbor, bond_order in node["edges"].items()
                                if neighbor in members
                            )
                        ),
                    )
                )
            components.append(tuple(rows))
        return tuple(sorted(components, key=lambda component: component[0][0]))

    def _capture_payload(
        self, snapshot: _Transaction, bindings: Mapping[str, str]
    ) -> dict[str, Any]:
        payload = {
            "before_strands": copy.deepcopy(snapshot.strands),
            "before_junctions": set(snapshot.incident_junctions),
            "before_ids": set(snapshot.strands),
            "bindings": dict(bindings),
        }
        payload["post_scope_hash"] = self._binding_scope_hash(bindings)
        return payload

    def _binding_scope_hash(self, bindings: Mapping[str, str]) -> str:
        owners = self._uuid_owners()
        starts = {owners[value] for value in bindings.values() if value in owners}
        component_ids = {self.components[strand_id] for strand_id in starts}
        strand_ids = {
            strand_id
            for strand_id, component_id in self.components.items()
            if component_id in component_ids
        }
        atom_uuids = {
            ref.uuid
            for strand_id in strand_ids
            for ref in self.strands[strand_id].atom_refs.values()
        }
        payload = {
            "strands": {
                strand_id: {
                    "length": self.strands[strand_id].length,
                    "features": self.strands[strand_id].features,
                    "formula": self.strands[strand_id].formula,
                    "radicals": self.strands[strand_id].radical_count,
                    "atoms": {
                        key: ref.__dict__
                        for key, ref in sorted(
                            self.strands[strand_id].atom_refs.items()
                        )
                    },
                    "atom_graph": self.strands[strand_id].atom_graph,
                    "repeat_unit_formula": self.strands[strand_id].repeat_unit_formula,
                }
                for strand_id in sorted(strand_ids)
            },
            "junctions": sorted(
                edge
                for edge in self.junctions
                if edge[0] in atom_uuids or edge[1] in atom_uuids
            ),
        }
        encoded = json.dumps(
            payload, sort_keys=True, separators=(",", ":"), default=str
        )
        return hashlib.sha256(encoded.encode()).hexdigest()

    def _restore_capture(self, capture: Mapping[str, Any]) -> set[str]:
        owners = self._uuid_owners()
        starts = {
            owners[atom_uuid]
            for atom_uuid in capture["bindings"].values()
            if atom_uuid in owners
        }
        component_ids = {self.components[strand_id] for strand_id in starts}
        current_ids = {
            strand_id
            for strand_id, component_id in self.components.items()
            if component_id in component_ids
        }
        current_uuids = {
            ref.uuid
            for strand_id in current_ids
            for ref in self.strands[strand_id].atom_refs.values()
        }
        self.junctions = {
            edge
            for edge in self.junctions
            if edge[0] not in current_uuids and edge[1] not in current_uuids
        }
        for strand_id in current_ids:
            self.strands.pop(strand_id, None)
        for strand_id, strand in capture["before_strands"].items():
            self.strands[strand_id] = copy.deepcopy(strand)
        self.junctions |= set(capture["before_junctions"])
        return current_ids | set(capture["before_ids"])

    def _find_reverse_capture(
        self,
        record: Any,
        participants: Sequence[Site],
        reverse_bindings: Mapping[str, str],
    ) -> tuple[int, dict[str, Any]]:
        reverse_of = _value(record, "reverse_of")
        candidates = self._open_captures.get(reverse_of, [])
        supplied = {
            node["atom_ref"].uuid
            for participant in participants
            for node in participant.graph.values()
            if isinstance(node.get("atom_ref"), AtomRef)
        }
        for index in range(len(candidates) - 1, -1, -1):
            capture = candidates[index]
            exact_labels = all(
                reverse_bindings.get(label) == atom_uuid
                for label, atom_uuid in capture["bindings"].items()
            )
            if exact_labels and set(capture["bindings"].values()) <= supplied:
                if (
                    self._binding_scope_hash(capture["bindings"])
                    != capture["post_scope_hash"]
                ):
                    raise ReverseBindingError(
                        "captured component changed after the R1 binding was created"
                    )
                return index, capture
        raise ReverseBindingError("reverse does not match an open R1 UUID binding")

    def apply(self, record: Any, participants: Sequence[Site]) -> AppliedEvent:
        """Apply one compiled record exactly, or leave state byte-identical."""
        snapshot = self._begin(participants, record)
        try:
            mapping, bindings, local_graph, reactant_graphs = self._preflight(
                record, participants
            )
            is_reverse = bool(
                _value(record, "reverse_of")
                and _value(record, "orientation", "") == "reversed"
                and any(
                    operation.get("action") == "dissociate"
                    for operation in (_value(record, "junction_ops", None) or [])
                )
            )
            reverse_capture = None
            if is_reverse:
                reverse_capture = self._find_reverse_capture(
                    record, participants, bindings
                )
            operations = self._record_operations(record)
            backbone_cuts, radical_change = self._rewrite_graph(
                operations, mapping, local_graph
            )
            declared_radical_delta = int(_value(record, "radical_delta", 0))
            if radical_change != declared_radical_delta:
                raise LedgerError(
                    f"mapped radical edits sum to {radical_change}, "
                    f"record declares {declared_radical_delta}"
                )

            cuts: list[int] = []
            inheritance: dict[str, Any] = {}
            affected = self._affected_strands(participants)
            affected |= self._transfer_mobile_atoms(local_graph)
            if is_reverse:
                assert reverse_capture is not None
                capture_index, capture = reverse_capture
                bindings = dict(capture["bindings"])
                affected |= self._restore_capture(capture)
                self._open_captures[_value(record, "reverse_of")].pop(capture_index)
            else:
                for first, second in sorted(
                    backbone_cuts,
                    key=lambda pair: min(
                        self._ref(pair[0]).position, self._ref(pair[1]).position
                    ),
                    reverse=True,
                ):
                    old_id = self._ref(first).strand_id
                    new_id, cut = self._split_strand(
                        first,
                        second,
                        {
                            atom_uuid: node.get("element")
                            for atom_uuid, node in local_graph.items()
                        },
                    )
                    affected.update((old_id, new_id))
                    cuts.append(cut)
                for first, second in backbone_cuts:
                    if self._ref(first).strand_id == self._ref(second).strand_id:
                        raise AssertionError(
                            "mapped backbone break did not split its strand"
                        )
                # A different form operation may reconnect the new fragments;
                # represent that mapped heavy-atom edge as a UUID junction.
                for first, node in local_graph.items():
                    for second in node["edges"]:
                        if first >= second:
                            continue
                        first_ref, second_ref = self._ref(first), self._ref(second)
                        if (
                            first_ref.strand_id != second_ref.strand_id
                            and node.get("element") != "H"
                            and local_graph[second].get("element") != "H"
                        ):
                            self.junctions.add(tuple(sorted((first, second))))
                formula_delta = _value(record, "formula_delta", None) or _value(
                    record, "element_delta", {}
                )
                owners = [
                    self._ref(atom_uuid).strand_id for atom_uuid in mapping.values()
                ]
                self._apply_formula_delta(formula_delta, owners)
                self._debit_coproducts(record, mapping, reactant_graphs, local_graph)

            self.assert_uuid_uniqueness()
            self._refresh_components(affected)

            for strand_id in sorted(affected & set(self.strands)):
                strand = self.strands[strand_id]
                inheritance[strand_id] = {
                    "length": strand.length,
                    "features": tuple(sorted(strand.features)),
                }
            formed_bond = None
            for operation in _value(record, "junction_ops", None) or []:
                labels = operation.get("formed_bond", {}).get("label_pair", [])
                if len(labels) == 2 and all(label in bindings for label in labels):
                    formed_bond = tuple(
                        sorted((bindings[labels[0]], bindings[labels[1]]))
                    )
                    break
            entry = AppliedEvent(
                str(_value(record, "event_id", "")),
                tuple(participant.strand_id for participant in participants),
                MappingProxyType(dict(bindings)),
                tuple(sorted(cuts)),
                MappingProxyType(inheritance),
                formed_bond,
                self._graph_components(local_graph),
            )
            self.event_log.append(entry)
            if any(
                operation.get("action") == "create"
                for operation in (_value(record, "junction_ops", None) or [])
            ):
                payload = self._capture_payload(snapshot, bindings)
                self._open_captures.setdefault(entry.event_id, []).append(payload)
            return entry
        except Exception:
            self._rollback(snapshot)
            raise
