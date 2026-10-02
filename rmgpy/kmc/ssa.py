"""Isothermal type-aggregated stochastic simulation for sparse KMC states.

The module owns propensity reduction, production site discovery, length-binned
pair thinning, and refused-channel accounting.  A caller supplies the sole
``numpy.random.Generator`` used by a trajectory.
"""

from __future__ import annotations

import math
import copy
from collections import Counter
from collections import OrderedDict
from collections.abc import Mapping
from dataclasses import dataclass
from itertools import product
from typing import Any, Hashable, Iterable, Sequence

from rmgpy.kmc.state import AtomRef, KMCState, Site, _parse_adjacency
from rmgpy.kmc.met import (
    N_A,
    PS_C_R2,
    PS_M0,
    R_RMG,
    SIGMA_CONTACT,
    CompiledChannel,
    CompiledTerminationTable,
    RateTable,
    TransportArm,
    collins_kimball,
    diffusion_rate,
    draw_ortho_attack_site,
)


E_PROPENSITY_NONFINITE = "E_PROPENSITY_NONFINITE"
E_PROPENSITY_NEGATIVE = "E_PROPENSITY_NEGATIVE"


class PropensityError(ValueError):
    """A canonical propensity is outside the registered numeric domain."""

    def __init__(self, code: str, event_id: str, value: float):
        self.code = code
        self.event_id = event_id
        self.value = value
        super().__init__(f"{code}: {event_id}={value!r}")


def epsilon_propensity(count: int) -> float:
    """Return the registered binary64 reduction tolerance ε_PROP(n)."""
    if count < 0:
        raise ValueError("propensity count must be non-negative")
    unit_roundoff = 2.0**-53
    depth = math.ceil(math.log2(count)) if count > 1 else 0
    gamma = depth * unit_roundoff / (1.0 - depth * unit_roundoff)
    return 8.0 * unit_roundoff + 2.0 * gamma


def pairwise_sum(values: Iterable[float]) -> float:
    """Sum using one fixed left-to-right binary reduction tree."""
    level = [float(value) for value in values]
    while len(level) > 1:
        level = [
            level[index] + level[index + 1] if index + 1 < len(level) else level[index]
            for index in range(0, len(level), 2)
        ]
    return level[0] if level else 0.0


@dataclass(frozen=True)
class CanonicalPropensities:
    """Validated propensities in canonical event-ID order."""

    event_ids: tuple[str, ...]
    values: tuple[float, ...]
    total: float

    @property
    def absorbing(self) -> bool:
        return self.total == 0.0

    def draw(self, rng: Any) -> str | None:
        """Draw an event ID, leaving ``rng`` untouched when absorbing."""
        if self.absorbing:
            return None
        threshold = float(rng.random()) * self.total
        # Traverse the same fixed binary tree used for the total.  A linear
        # cumulative walk would reintroduce O(n*u) probability-bound error.
        nodes = [
            (value, event_id, None, None)
            for event_id, value in zip(self.event_ids, self.values)
        ]
        while len(nodes) > 1:
            next_level = []
            for index in range(0, len(nodes), 2):
                left = nodes[index]
                if index + 1 == len(nodes):
                    next_level.append(left)
                    continue
                right = nodes[index + 1]
                next_level.append((left[0] + right[0], None, left, right))
            nodes = next_level
        node = nodes[0]
        if threshold >= node[0]:
            threshold = math.nextafter(node[0], 0.0)
        while node[1] is None:
            left = node[2]
            right = node[3]
            assert left is not None and right is not None
            if threshold < left[0]:
                node = left
            else:
                threshold -= left[0]
                node = right
        return node[1]


class _PropensityTree:
    """Update the original binary reduction tree without changing its leaves."""

    def __init__(self, values: Mapping[str, float]):
        vector = canonical_propensities(values.items())
        self.event_ids = vector.event_ids
        self.positions = {key: index for index, key in enumerate(self.event_ids)}
        self.levels = [list(vector.values)]
        while len(self.levels[-1]) > 1:
            level = self.levels[-1]
            self.levels.append(
                [
                    (
                        level[index] + level[index + 1]
                        if index + 1 < len(level)
                        else level[index]
                    )
                    for index in range(0, len(level), 2)
                ]
            )

    @property
    def total(self) -> float:
        return self.levels[-1][0] if self.event_ids else 0.0

    @property
    def absorbing(self) -> bool:
        return self.total == 0.0

    def update(self, event_id: str, value: float) -> None:
        if not math.isfinite(value):
            raise PropensityError(E_PROPENSITY_NONFINITE, event_id, value)
        if value < 0.0:
            raise PropensityError(E_PROPENSITY_NEGATIVE, event_id, value)
        position = self.positions[event_id]
        self.levels[0][position] = value
        for depth in range(1, len(self.levels)):
            position //= 2
            lower = self.levels[depth - 1]
            left = 2 * position
            self.levels[depth][position] = (
                lower[left] + lower[left + 1] if left + 1 < len(lower) else lower[left]
            )

    def draw(self, rng: Any) -> str | None:
        if self.absorbing:
            return None
        threshold = float(rng.random()) * self.total
        if threshold >= self.total:
            threshold = math.nextafter(self.total, 0.0)
        position = 0
        for level in reversed(self.levels[:-1]):
            position *= 2
            # An unpaired node is carried unchanged to its parent.
            if position + 1 < len(level) and threshold >= level[position]:
                threshold -= level[position]
                position += 1
        return self.event_ids[position]


def canonical_propensities(
    propensities: Iterable[tuple[str, float]],
) -> CanonicalPropensities:
    """Validate, order, and reduce canonical propensities.

    Every value is checked before the first addition.  Duplicate identifiers
    are rejected because an event must occupy exactly one leaf of the fixed
    reduction tree.
    """
    items = sorted((str(event_id), float(value)) for event_id, value in propensities)
    identifiers = [event_id for event_id, _ in items]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("duplicate canonical propensity event ID")
    for event_id, value in items:
        if not math.isfinite(value):
            raise PropensityError(E_PROPENSITY_NONFINITE, event_id, value)
        if value < 0.0:
            raise PropensityError(E_PROPENSITY_NEGATIVE, event_id, value)
    values = tuple(value for _, value in items)
    total = pairwise_sum(values)
    return CanonicalPropensities(tuple(identifiers), values, total)


def _field(record: Any, name: str, default: Any = None) -> Any:
    return (
        record.get(name, default)
        if isinstance(record, Mapping)
        else getattr(record, name, default)
    )


def _expanded_site_types(record: Any) -> tuple[str, ...]:
    types = tuple(_field(record, "participant_site_types", ()))
    multiplicities = tuple(_field(record, "reactant_multiplicities", ()))
    if multiplicities and len(types) == len(multiplicities):
        return tuple(
            site_type
            for site_type, count in zip(types, multiplicities)
            for _ in range(int(count))
        )
    return types


def _record_labels(record: Any) -> dict[int, str]:
    labels: dict[int, str] = {}
    for junction in _field(record, "junction_ops", ()) or ():
        action = "form" if junction.get("action") == "create" else "break"
        operation = next(
            (
                item
                for item in _field(record, "bond_ops", ())
                if item.get("action") == action and len(item.get("atoms", ())) == 2
            ),
            None,
        )
        pair = junction.get("formed_bond", {}).get("label_pair", ())
        if operation is not None and len(pair) == 2:
            labels.update(zip(operation["atoms"], pair))
    return labels


def _node_matches(expected: Mapping[str, Any], actual: Mapping[str, Any]) -> bool:
    return expected.get("element") == actual.get("element") and all(
        int(expected.get(field, 0)) == int(actual.get(field, 0))
        for field in ("radical", "charge", "lone_pairs", "implicit_hydrogens")
    )


def _node_signature(node: Mapping[str, Any]) -> tuple[Any, ...]:
    return (
        node.get("element"),
        int(node.get("radical", 0)),
        int(node.get("charge", 0)),
        int(node.get("lone_pairs", 0)),
        int(node.get("implicit_hydrogens", 0)),
    )


def _neighborhood_can_contain(
    source: Any,
    target: Any,
    expected: Mapping[Any, Mapping[str, Any]],
    actual: Mapping[Any, Mapping[str, Any]],
) -> bool:
    required = Counter(
        (float(order), _node_signature(expected[neighbor]))
        for neighbor, order in expected[source].get("edges", {}).items()
    )
    available = Counter(
        (float(order), _node_signature(actual[neighbor]))
        for neighbor, order in actual[target].get("edges", {}).items()
        if neighbor in actual
    )
    return all(available[key] >= count for key, count in required.items())


def _edge_order(node: Mapping[str, Any], neighbor: Any) -> float | None:
    value = node.get("edges", {}).get(neighbor)
    return None if value is None else float(value)


def _graph_components(state: KMCState) -> dict[tuple[str, ...], dict[str, Any]]:
    nodes = {
        atom_uuid: node
        for strand in state.strands.values()
        for atom_uuid, node in strand.atom_graph.items()
    }
    unseen = set(nodes)
    components: dict[tuple[str, ...], dict[str, Any]] = {}
    while unseen:
        root = min(unseen)
        members: set[str] = set()
        stack = [root]
        while stack:
            current = stack.pop()
            if current in members:
                continue
            members.add(current)
            stack.extend(
                neighbor
                for neighbor in nodes[current].get("edges", {})
                if neighbor in nodes and neighbor not in members
            )
        unseen -= members
        key = tuple(sorted(members))
        components[key] = {atom_uuid: nodes[atom_uuid] for atom_uuid in key}
    return components


def _subgraph_mappings(
    expected: Mapping[int, Mapping[str, Any]],
    component: Mapping[str, Mapping[str, Any]],
) -> tuple[dict[int, str], ...]:
    """Enumerate one deterministic induced match for every possible anchor."""
    if not expected:
        return ()
    if len(expected) > len(component):
        return ()
    expected_nodes = Counter(
        (
            node.get("element"),
            int(node.get("radical", 0)),
            int(node.get("charge", 0)),
            int(node.get("lone_pairs", 0)),
            int(node.get("implicit_hydrogens", 0)),
        )
        for node in expected.values()
    )
    actual_nodes = Counter(
        (
            node.get("element"),
            int(node.get("radical", 0)),
            int(node.get("charge", 0)),
            int(node.get("lone_pairs", 0)),
            int(node.get("implicit_hydrogens", 0)),
        )
        for node in component.values()
    )
    if any(actual_nodes[key] < count for key, count in expected_nodes.items()):
        return ()
    expected_edges = Counter(
        float(order)
        for source, node in expected.items()
        for neighbor, order in node.get("edges", {}).items()
        if source < neighbor
    )
    actual_edges = Counter(
        float(order)
        for source, node in component.items()
        for neighbor, order in node.get("edges", {}).items()
        if source < neighbor
    )
    if any(actual_edges[key] < count for key, count in expected_edges.items()):
        return ()
    anchor = min(expected)
    candidates = {
        source: tuple(
            target
            for target in sorted(component)
            if _node_matches(expected_node, component[target])
            and _neighborhood_can_contain(source, target, expected, component)
        )
        for source, expected_node in expected.items()
    }
    if any(not values for values in candidates.values()):
        return ()
    results = []
    for anchor_target in candidates[anchor]:
        found: list[dict[int, str]] = []

        def visit(mapping: dict[int, str]) -> None:
            if found:
                return
            if len(mapping) == len(expected):
                found.append(dict(mapping))
                return
            remaining = set(expected) - set(mapping)
            source = (
                anchor
                if not mapping
                else min(
                    remaining,
                    key=lambda item: (
                        -sum(
                            neighbor in mapping
                            for neighbor in expected[item].get("edges", {})
                        ),
                        len(candidates[item]),
                        item,
                    ),
                )
            )
            targets = (anchor_target,) if source == anchor else candidates[source]
            for target in targets:
                if target in mapping.values():
                    continue
                compatible = True
                for other_source, other_target in mapping.items():
                    expected_order = _edge_order(expected[source], other_source)
                    actual_order = _edge_order(component[target], other_target)
                    if expected_order != actual_order:
                        compatible = False
                        break
                if compatible:
                    mapping[source] = target
                    visit(mapping)
                    del mapping[source]
                    if found:
                        return

        visit({})
        if found:
            results.append(
                min(
                    found,
                    key=lambda mapping: tuple(
                        mapping[source] for source in sorted(mapping)
                    ),
                )
            )
    return tuple(results)


@dataclass(frozen=True)
class _PreparedRecord:
    record: Any
    event_id: str
    graphs: tuple[dict[int, dict[str, Any]], ...]
    graph_keys: tuple[Hashable, ...]
    offsets: tuple[int, ...]
    labels: Mapping[int, str]
    site_types: tuple[str, ...]


def _prepare_record(record: Any, graph_cache: dict, graph_ids: dict) -> _PreparedRecord:
    graphs, graph_keys = [], []
    for text in _field(record, "reactant_graphs", ()):
        if text not in graph_cache:
            graph = _parse_adjacency(text)
            key = _graph_key(graph)
            graph_cache[text] = graph, graph_ids.setdefault(key, len(graph_ids))
        graph, key = graph_cache[text]
        graphs.append(graph)
        graph_keys.append(key)
    offsets = []
    offset = 0
    for graph in graphs:
        offsets.append(offset)
        offset += len(graph)
    site_types = _expanded_site_types(record)
    if len(graphs) != len(site_types):
        raise ValueError(
            f"record {_field(record, 'event_id', '')} graph/site arity disagrees"
        )
    return _PreparedRecord(
        record,
        str(_field(record, "event_id", "")),
        tuple(graphs),
        tuple(graph_keys),
        tuple(offsets),
        _record_labels(record),
        site_types,
    )


def _graph_key(graph: Mapping[int, Mapping[str, Any]]) -> Hashable:
    return tuple(
        (
            index,
            node.get("element"),
            int(node.get("radical", 0)),
            int(node.get("charge", 0)),
            int(node.get("lone_pairs", 0)),
            int(node.get("implicit_hydrogens", 0)),
            tuple(
                sorted(
                    (int(neighbor), float(order))
                    for neighbor, order in node.get("edges", {}).items()
                )
            ),
        )
        for index, node in sorted(graph.items())
    )


def _site_key(site: Site) -> tuple[str, ...]:
    return tuple(
        sorted(
            node["atom_ref"].uuid
            for node in site.graph.values()
            if isinstance(node.get("atom_ref"), AtomRef)
        )
    )


class SiteIndex:
    """Incremental production index of compiled-record participant tuples.

    "Touched" means every atom-graph component owned by a polymer connected
    component containing a participant immediately before the rewrite, plus
    every resulting strand reported in ``AppliedEvent.inheritance``.  Only
    atom-graph components intersecting that set are structurally re-matched;
    candidate products are then rebuilt from the retained and refreshed local
    matches.
    """

    def __init__(self, state: KMCState, records: Sequence[Any]):
        self.state = state
        self.records = tuple(
            sorted(records, key=lambda record: str(_field(record, "event_id", "")))
        )
        graph_cache, graph_ids = {}, {}
        self._prepared = tuple(
            _prepare_record(record, graph_cache, graph_ids) for record in self.records
        )
        identifiers = [prepared.event_id for prepared in self._prepared]
        if not all(identifiers) or len(identifiers) != len(set(identifiers)):
            raise ValueError("site-index records need unique event IDs")
        self._graphs = {}
        self._records_for_graph: dict[Hashable, set[int]] = {}
        for ordinal, prepared in enumerate(self._prepared):
            for key, graph in zip(prepared.graph_keys, prepared.graphs):
                self._graphs[key] = graph
                self._records_for_graph.setdefault(key, set()).add(ordinal)
        self._components: dict[tuple[str, ...], dict[str, Any]] = {}
        self._matches: dict[
            tuple[Hashable, tuple[str, ...]], tuple[dict[int, str], ...]
        ] = {}
        self._candidates: dict[str, tuple[tuple[Site, ...], ...]] = dict.fromkeys(
            identifiers, ()
        )
        self._state_cache: OrderedDict[str, tuple[Any, Any, Any]] = OrderedDict()
        self.last_rescanned_atom_uuids: frozenset[str] = frozenset()
        self.revision = 0
        self.changed_event_ids: frozenset[str] = frozenset()
        self.active_event_ids: frozenset[str] = frozenset()
        self._population_signatures: dict[str, Any] = {}
        self._eligible: dict[str, tuple[tuple[Site, ...], ...]] = {}
        self._full_rescan()
        self._update_populations()
        self._remember_state()

    def _update_populations(self) -> None:
        signatures = {}
        for event_id, candidates in self._candidates.items():
            if candidates:
                signatures[event_id] = tuple(
                    tuple(
                        (
                            _site_key(site),
                            site.strand_id,
                            site.position,
                            tuple(
                                (
                                    local,
                                    node["atom_ref"].uuid,
                                    node["atom_ref"].position,
                                    node["atom_ref"].role,
                                )
                                for local, node in sorted(site.graph.items())
                            ),
                            self.state.components[site.strand_id],
                            self.state.component_lengths[
                                self.state.components[site.strand_id]
                            ],
                        )
                        for site in candidate
                    )
                    for candidate in candidates
                )
        self.changed_event_ids = frozenset(
            key
            for key in self._population_signatures.keys() | signatures.keys()
            if self._population_signatures.get(key) != signatures.get(key)
        )
        self.active_event_ids = frozenset(signatures)
        self._population_signatures = signatures
        for key in self.changed_event_ids:
            self._eligible.pop(key, None)
        self.revision += 1

    def _remember_state(self) -> None:
        state_hash = self.state.state_hash()
        self._state_cache[state_hash] = (
            dict(self._components),
            dict(self._matches),
            self._candidates,
        )
        self._state_cache.move_to_end(state_hash)
        while len(self._state_cache) > 32:
            self._state_cache.popitem(last=False)

    def _scan_components(self, component_keys: Iterable[tuple[str, ...]]) -> None:
        keys = tuple(component_keys)
        for graph_key, graph in self._graphs.items():
            for key in keys:
                self._matches[(graph_key, key)] = _subgraph_mappings(
                    graph, self._components[key]
                )

    def _full_rescan(self) -> None:
        self._components = _graph_components(self.state)
        self._matches.clear()
        self._scan_components(self._components)
        self.last_rescanned_atom_uuids = frozenset(
            atom_uuid for key in self._components for atom_uuid in key
        )
        self._rebuild_candidates()

    def _reverse_mapping_allowed(
        self, prepared: _PreparedRecord, mapping: Mapping[int, str]
    ) -> bool:
        record = prepared.record
        if not (
            _field(record, "orientation", "") == "reversed"
            and _field(record, "inventory_class") == "R1:J_ring"
        ):
            return True
        captures = self.state._open_captures.get(_field(record, "reverse_of"), ())
        return any(
            all(
                mapping[index] == capture["bindings"].get(label)
                for index, label in prepared.labels.items()
            )
            for capture in captures
        )

    def _make_site(
        self,
        prepared: _PreparedRecord,
        ordinal: int,
        mapping: Mapping[int, str],
    ) -> Site:
        graph = copy.deepcopy(prepared.graphs[ordinal])
        offset = prepared.offsets[ordinal]
        for local, node in graph.items():
            ref = self.state._ref(mapping[local])
            node.update({"atom_ref": ref, "position": ref.position, "role": ref.role})
            label = prepared.labels.get(offset + local)
            if label:
                node["label"] = label
        anchor = min(graph)
        anchor_ref = graph[anchor]["atom_ref"]
        return Site(
            prepared.site_types[ordinal],
            anchor_ref.strand_id,
            anchor_ref.position,
            anchor,
            graph,
        )

    def _rebuild_candidates(self) -> None:
        candidates = self._candidates.copy()
        for event_id in self.active_event_ids:
            candidates[event_id] = ()
        component_keys = tuple(sorted(self._components))
        groups_by_graph = {}
        possible = set()
        for graph_key in self._graphs:
            matches = tuple(
                (key, mapping)
                for key in component_keys
                for mapping in self._matches.get((graph_key, key), ())
            )
            if matches:
                groups_by_graph[graph_key] = matches
                possible.update(self._records_for_graph[graph_key])
        for record_ordinal in sorted(possible):
            prepared = self._prepared[record_ordinal]
            groups = [groups_by_graph.get(key, ()) for key in prepared.graph_keys]
            found = []
            if groups and all(groups):
                for combination in product(*groups):
                    mapped = [set(mapping.values()) for _, mapping in combination]
                    if any(
                        left & right
                        for index, left in enumerate(mapped)
                        for right in mapped[index + 1 :]
                    ):
                        continue
                    global_mapping = {
                        prepared.offsets[ordinal] + local: atom_uuid
                        for ordinal, (_, mapping) in enumerate(combination)
                        for local, atom_uuid in mapping.items()
                    }
                    if not self._reverse_mapping_allowed(prepared, global_mapping):
                        continue
                    found.append(
                        tuple(
                            self._make_site(prepared, ordinal, mapping)
                            for ordinal, (_, mapping) in enumerate(combination)
                        )
                    )
            candidates[prepared.event_id] = tuple(
                sorted(found, key=lambda item: tuple(_site_key(site) for site in item))
            )
        self._candidates = candidates

    def candidates(self, event_id: str) -> tuple[tuple[Site, ...], ...]:
        return self._candidates[str(event_id)]

    def apply(self, record: Any, participants: Sequence[Site]) -> Any:
        """Apply through the product executor and refresh only touched matches."""
        before_component_ids = {
            self.state.components[site.strand_id]
            for site in participants
            if site.strand_id in self.state.components
        }
        before_strands = {
            strand_id
            for strand_id, component_id in self.state.components.items()
            if component_id in before_component_ids
        }
        before_uuids = {
            ref.uuid
            for strand_id in before_strands
            for ref in self.state.strands[strand_id].atom_refs.values()
        }
        entry = self.state.apply(record, participants)
        after_strands = set(entry.inheritance) & set(self.state.strands)
        after_uuids = {
            ref.uuid
            for strand_id in after_strands
            for ref in self.state.strands[strand_id].atom_refs.values()
        }
        touched = before_uuids | after_uuids
        state_hash = self.state.state_hash()
        cached = self._state_cache.get(state_hash)
        if cached is not None:
            _, matches, candidates = cached
            self._components = _graph_components(self.state)
            self._matches = dict(matches)
            self._candidates = candidates
            self._state_cache.move_to_end(state_hash)
            self.last_rescanned_atom_uuids = frozenset()
            self._update_populations()
            return entry
        current = _graph_components(self.state)
        stale_keys = {
            key for key in self._components if set(key) & touched or key not in current
        }
        for cache_key in tuple(self._matches):
            if cache_key[1] in stale_keys:
                del self._matches[cache_key]
        refreshed = {
            key for key in current if set(key) & touched or key not in self._components
        }
        self._components = current
        self._scan_components(sorted(refreshed))
        self.last_rescanned_atom_uuids = frozenset(
            atom_uuid for key in refreshed for atom_uuid in key
        )
        self._rebuild_candidates()
        self._update_populations()
        self._remember_state()
        return entry

    def _candidate_keys(self) -> dict[str, tuple[tuple[tuple[str, ...], ...], ...]]:
        return {
            event_id: tuple(
                tuple(_site_key(site) for site in candidate) for candidate in candidates
            )
            for event_id, candidates in self._candidates.items()
        }

    def assert_index_consistency(self) -> None:
        """Assert equality with a new full structural rescan."""
        full = SiteIndex(self.state, self.records)
        if self._candidate_keys() != full._candidate_keys():
            raise AssertionError("incremental site index differs from full rescan")


class ThinningBoundError(RuntimeError):
    """An exact pair rate exceeded its declared thinning bound."""


@dataclass(frozen=True)
class PairItem:
    """One component-level participant in a binned pair population."""

    item_id: str
    component_id: str
    length: int
    kind: Hashable
    payload: Any = None

    def __post_init__(self) -> None:
        if self.length < 1:
            raise ValueError("pair-item length must be positive")


@dataclass(frozen=True)
class PairCell:
    cell_id: str
    left_kind: Hashable
    right_kind: Hashable
    left_bin: int
    right_bin: int
    left_range: tuple[int, int]
    right_range: tuple[int, int]
    pair_count: int
    rate_bound: float
    propensity_bound: float


@dataclass(frozen=True)
class PairProposal:
    elapsed: float
    cell: PairCell
    left: PairItem
    right: PairItem
    exact_rate: float
    accepted: bool


class ConstantPairKernel:
    """Constant coagulation kernel used by the registered C4 verifier."""

    def __init__(self, rate: float, *, bound_scale: float = 1.0):
        self.rate = float(rate)
        self.bound_scale = float(bound_scale)

    def exact_rate(self, left: PairItem, right: PairItem) -> float:
        return self.rate

    def bound_rate(
        self,
        left_kind: Hashable,
        right_kind: Hashable,
        left_range: tuple[int, int],
        right_range: tuple[int, int],
    ) -> float:
        return self.rate * self.bound_scale


class SumPairKernel:
    """Additive ``K(i,j)=b(i+j)`` coagulation kernel for C4."""

    def __init__(self, coefficient: float, *, bound_scale: float = 1.0):
        self.coefficient = float(coefficient)
        self.bound_scale = float(bound_scale)

    def exact_rate(self, left: PairItem, right: PairItem) -> float:
        return self.coefficient * (left.length + right.length)

    def bound_rate(
        self,
        left_kind: Hashable,
        right_kind: Hashable,
        left_range: tuple[int, int],
        right_range: tuple[int, int],
    ) -> float:
        return self.coefficient * (left_range[1] + right_range[1]) * self.bound_scale


class _Bucket:
    def __init__(self) -> None:
        self.items: list[PairItem] = []
        self.positions: dict[str, int] = {}
        self.components: set[str] = set()

    def add(self, item: PairItem) -> None:
        if item.item_id in self.positions:
            raise ValueError(f"duplicate pair item {item.item_id}")
        if item.component_id in self.components:
            raise ValueError("one pair kind may occur only once per component")
        self.positions[item.item_id] = len(self.items)
        self.items.append(item)
        self.components.add(item.component_id)

    def remove(self, item: PairItem) -> None:
        position = self.positions.pop(item.item_id)
        last = self.items.pop()
        if position < len(self.items):
            self.items[position] = last
            self.positions[last.item_id] = position
        self.components.remove(item.component_id)


def length_bin(length: int) -> tuple[int, int, int]:
    """Return the base-two bin index and inclusive range for a length."""
    if length < 1:
        raise ValueError("length must be positive")
    index = int(length).bit_length() - 1
    return index, 1 << index, (1 << (index + 1)) - 1


class LengthBinnedThinningSampler:
    """Production component-pair sampler using exact accept/reject thinning.

    Pair populations are component-level: an item kind can occur at most once
    on a component, and cross-kind cells subtract component IDs present on both
    sides.  Thus an intracomponent biradical pair is never proposed.
    """

    def __init__(
        self,
        items: Iterable[PairItem],
        kernel: Any,
        pair_kinds: Iterable[tuple[Hashable, Hashable]],
        *,
        volume: float,
        normalization: float = 1.0,
        binning: Any = length_bin,
    ):
        self.kernel = kernel
        self.volume = float(volume)
        self.normalization = float(normalization)
        self.binning = binning
        if not math.isfinite(self.volume) or self.volume <= 0.0:
            raise ValueError("volume must be finite and positive")
        if not math.isfinite(self.normalization) or self.normalization <= 0.0:
            raise ValueError("normalization must be finite and positive")
        self.pair_kinds = tuple(pair_kinds)
        if len(self.pair_kinds) != len(set(self.pair_kinds)):
            raise ValueError("duplicate pair-kind convention")
        self.items: dict[str, PairItem] = {}
        self._buckets: dict[tuple[Hashable, int], _Bucket] = {}
        self._bin_ranges: dict[tuple[Hashable, int], tuple[int, int]] = {}
        self.bound_checks = 0
        self.bound_violations = 0
        for item in items:
            self.add(item, rebuild=False)
        self._cells = self._build_cells()

    def _bucket(self, kind: Hashable, bin_index: int) -> _Bucket:
        return self._buckets.setdefault((kind, bin_index), _Bucket())

    def add(self, item: PairItem, *, rebuild: bool = True) -> None:
        if item.item_id in self.items:
            raise ValueError(f"duplicate pair item {item.item_id}")
        bin_index, lower, upper = self.binning(item.length)
        self._bucket(item.kind, bin_index).add(item)
        self._bin_ranges[(item.kind, bin_index)] = (lower, upper)
        self.items[item.item_id] = item
        if rebuild:
            self._cells = self._build_cells()

    def remove(self, item_id: str, *, rebuild: bool = True) -> PairItem:
        item = self.items.pop(item_id)
        bin_index, _, _ = self.binning(item.length)
        key = (item.kind, bin_index)
        bucket = self._buckets[key]
        bucket.remove(item)
        if not bucket.items:
            del self._buckets[key]
            del self._bin_ranges[key]
        if rebuild:
            self._cells = self._build_cells()
        return item

    def replace_pair(self, left_id: str, right_id: str, product_item: PairItem) -> None:
        if left_id == right_id:
            raise ValueError("a pair replacement needs two items")
        self.remove(left_id, rebuild=False)
        self.remove(right_id, rebuild=False)
        self.add(product_item, rebuild=False)
        self._cells = self._build_cells()

    def _build_cells(self) -> tuple[PairCell, ...]:
        cells = []
        serial = 0
        for left_kind, right_kind in self.pair_kinds:
            left_bins = sorted(
                bin_index for kind, bin_index in self._buckets if kind == left_kind
            )
            right_bins = sorted(
                bin_index for kind, bin_index in self._buckets if kind == right_kind
            )
            for left_bin in left_bins:
                for right_bin in right_bins:
                    if left_kind == right_kind and right_bin < left_bin:
                        continue
                    left = self._buckets[(left_kind, left_bin)]
                    right = self._buckets[(right_kind, right_bin)]
                    if left is right:
                        pair_count = len(left.items) * (len(left.items) - 1) // 2
                    else:
                        pair_count = len(left.items) * len(right.items) - len(
                            left.components & right.components
                        )
                    if pair_count <= 0:
                        continue
                    left_range = self._bin_ranges[(left_kind, left_bin)]
                    right_range = self._bin_ranges[(right_kind, right_bin)]
                    bound = float(
                        self.kernel.bound_rate(
                            left_kind,
                            right_kind,
                            left_range,
                            right_range,
                        )
                    )
                    if not math.isfinite(bound):
                        raise PropensityError(
                            E_PROPENSITY_NONFINITE, f"pair-cell-{serial}", bound
                        )
                    if bound < 0.0:
                        raise PropensityError(
                            E_PROPENSITY_NEGATIVE, f"pair-cell-{serial}", bound
                        )
                    propensity = bound * pair_count / (self.normalization * self.volume)
                    cells.append(
                        PairCell(
                            f"pair-cell-{serial:08d}",
                            left_kind,
                            right_kind,
                            left_bin,
                            right_bin,
                            left_range,
                            right_range,
                            pair_count,
                            bound,
                            propensity,
                        )
                    )
                    serial += 1
        return tuple(cells)

    @property
    def cells(self) -> tuple[PairCell, ...]:
        return self._cells

    @property
    def total_bound_propensity(self) -> float:
        return canonical_propensities(
            (cell.cell_id, cell.propensity_bound) for cell in self.cells
        ).total

    def _draw_item_pair(self, cell: PairCell, rng: Any) -> tuple[PairItem, PairItem]:
        left_bucket = self._buckets[(cell.left_kind, cell.left_bin)]
        right_bucket = self._buckets[(cell.right_kind, cell.right_bin)]
        if left_bucket is right_bucket:
            size = len(left_bucket.items)
            first = int(rng.integers(size))
            second = int(rng.integers(size - 1))
            if second >= first:
                second += 1
            if second < first:
                first, second = second, first
            return left_bucket.items[first], left_bucket.items[second]
        while True:
            left = left_bucket.items[int(rng.integers(len(left_bucket.items)))]
            right = right_bucket.items[int(rng.integers(len(right_bucket.items)))]
            if left.component_id != right.component_id:
                return left, right

    def select_pair(
        self,
        rng: Any,
        vector: CanonicalPropensities | None = None,
    ) -> PairProposal | None:
        """Select and test one bound proposal without drawing elapsed time."""
        if vector is None:
            vector = canonical_propensities(
                (cell.cell_id, cell.propensity_bound) for cell in self.cells
            )
        selected_id = vector.draw(rng)
        if selected_id is None:
            return None
        cell = next(cell for cell in self.cells if cell.cell_id == selected_id)
        left, right = self._draw_item_pair(cell, rng)
        exact = float(self.kernel.exact_rate(left, right))
        if not math.isfinite(exact):
            raise PropensityError(E_PROPENSITY_NONFINITE, cell.cell_id, exact)
        if exact < 0.0:
            raise PropensityError(E_PROPENSITY_NEGATIVE, cell.cell_id, exact)
        self.bound_checks += 1
        tolerance = 8.0 * math.ulp(max(exact, cell.rate_bound, 1.0))
        if exact > cell.rate_bound + tolerance:
            self.bound_violations += 1
            raise ThinningBoundError(
                f"exact pair rate {exact:.17g} exceeds bound "
                f"{cell.rate_bound:.17g} in {cell.cell_id}"
            )
        accepted = bool(rng.random() * cell.rate_bound < exact)
        return PairProposal(0.0, cell, left, right, exact, accepted)

    def propose(self, rng: Any) -> PairProposal | None:
        vector = canonical_propensities(
            (cell.cell_id, cell.propensity_bound) for cell in self.cells
        )
        if vector.absorbing:
            return None
        elapsed = float(rng.exponential(1.0 / vector.total))
        proposal = self.select_pair(rng, vector)
        assert proposal is not None
        return PairProposal(
            elapsed,
            proposal.cell,
            proposal.left,
            proposal.right,
            proposal.exact_rate,
            proposal.accepted,
        )

    def draw_accepted(self, rng: Any) -> PairProposal | None:
        elapsed = 0.0
        while True:
            proposal = self.propose(rng)
            if proposal is None:
                return None
            elapsed += proposal.elapsed
            if proposal.accepted:
                return PairProposal(
                    elapsed,
                    proposal.cell,
                    proposal.left,
                    proposal.right,
                    proposal.exact_rate,
                    True,
                )

    def cell_pairs(self, cell: PairCell) -> Iterable[tuple[PairItem, PairItem]]:
        """Iterate valid component pairs in one cell for diagnostics/oracles."""
        left = self._buckets[(cell.left_kind, cell.left_bin)].items
        right = self._buckets[(cell.right_kind, cell.right_bin)].items
        if left is right:
            for left_index, left_item in enumerate(left):
                for right_item in left[left_index + 1 :]:
                    yield left_item, right_item
        else:
            for left_item in left:
                for right_item in right:
                    if left_item.component_id != right_item.component_id:
                        yield left_item, right_item

    def exact_cell_propensities(self) -> dict[str, float]:
        """Return diagnostic exact cell propensities by explicit pair sum."""
        result = {}
        for cell in self.cells:
            rates = [
                float(self.kernel.exact_rate(left, right))
                for left, right in self.cell_pairs(cell)
            ]
            result[cell.cell_id] = pairwise_sum(rates) / (
                self.normalization * self.volume
            )
        return result

    @property
    def exact_total_propensity(self) -> float:
        return pairwise_sum(self.exact_cell_propensities().values())


def radical_site_class(site: Site, state: KMCState) -> str:
    """Classify a radical site as ``end`` or ``mid`` from persistent positions."""
    radical_refs = [
        node.get("atom_ref")
        for node in site.graph.values()
        if int(node.get("radical", 0)) > 0
    ]
    radical_refs = [ref for ref in radical_refs if isinstance(ref, AtomRef)]
    if not radical_refs:
        raise ValueError("radical site has no radical atom")
    for ref in radical_refs:
        strand = state.strands[ref.strand_id]
        if ref.position in {0, strand.length - 1}:
            return "end"
    return "mid"


def radical_pair_class(left: Site, right: Site, state: KMCState) -> str:
    classes = sorted(
        (radical_site_class(left, state), radical_site_class(right, state))
    )
    return f"{classes[0]}/{classes[1]}".replace("end/mid", "end/mid")


def _pair_type(left: str, right: str) -> tuple[str, str]:
    return tuple(sorted((str(left), str(right))))


class METPairKernel:
    """Exact Collins–Kimball rate and a proven decoupled length-bin bound."""

    def __init__(
        self,
        activation_rates: Mapping[tuple[str, str], float],
        arm: TransportArm,
        temperature: float,
        *,
        spin_factor: float,
        bound_temperature: float | None = None,
        bound_activation_rates: Mapping[tuple[str, str], float] | None = None,
        component_activation_rates: (
            Mapping[tuple[str, str, tuple[str, str]], float] | None
        ) = None,
        component_bound_activation_rates: (
            Mapping[tuple[str, str, tuple[str, str]], float] | None
        ) = None,
    ):
        self.activation_rates = {
            _pair_type(*pair): float(rate) for pair, rate in activation_rates.items()
        }
        self.bound_activation_rates = {
            _pair_type(*pair): float(rate)
            for pair, rate in (
                activation_rates
                if bound_activation_rates is None
                else bound_activation_rates
            ).items()
        }
        self.arm = arm
        self.temperature = float(temperature)
        self.bound_temperature = float(
            temperature if bound_temperature is None else bound_temperature
        )
        self.spin_factor = float(spin_factor)
        self.component_activation_rates = (
            dict(component_activation_rates)
            if component_activation_rates is not None
            else None
        )
        self.component_bound_activation_rates = (
            dict(component_activation_rates)
            if component_bound_activation_rates is None
            and component_activation_rates is not None
            else (
                dict(component_bound_activation_rates)
                if component_bound_activation_rates is not None
                else None
            )
        )
        if not 0.25 <= self.spin_factor <= 1.0:
            raise ValueError("spin factor must lie in [1/4, 1]")
        if self.bound_temperature < self.temperature:
            raise ValueError("bound temperature must not be below run temperature")

    @staticmethod
    def _site_type(kind: Hashable) -> str:
        return str(kind[0] if isinstance(kind, tuple) else kind)

    @staticmethod
    def _site_class(kind: Hashable) -> str:
        if not isinstance(kind, tuple) or len(kind) < 2:
            raise ValueError("MET item kind must carry (site_type, end_or_mid)")
        value = str(kind[1])
        if value not in {"end", "mid"}:
            raise ValueError(f"unknown radical site class {value!r}")
        return value

    def _activation(self, left_kind: Hashable, right_kind: Hashable) -> float:
        pair = _pair_type(self._site_type(left_kind), self._site_type(right_kind))
        try:
            return self.activation_rates[pair]
        except KeyError as error:
            raise ValueError(f"no MET activation rate for pair {pair}") from error

    def _class(self, left_kind: Hashable, right_kind: Hashable) -> str:
        values = sorted((self._site_class(left_kind), self._site_class(right_kind)))
        return f"{values[0]}/{values[1]}"

    def exact_rate(self, left: PairItem, right: PairItem) -> float:
        activation = self._activation(left.kind, right.kind)
        if self.component_activation_rates is not None:
            components = tuple(sorted((left.component_id, right.component_id)))
            pair = _pair_type(self._site_type(left.kind), self._site_type(right.kind))
            activation = self.component_activation_rates.get(
                (components[0], components[1], pair), 0.0
            )
            if activation == 0.0:
                return 0.0
        return collins_kimball(
            activation,
            diffusion_rate(
                self.arm,
                self.temperature,
                left.length,
                right.length,
                self._class(left.kind, right.kind),
                spin_factor=self.spin_factor,
            ),
        )

    def bound_rate(
        self,
        left_kind: Hashable,
        right_kind: Hashable,
        left_range: tuple[int, int],
        right_range: tuple[int, int],
    ) -> float:
        # D(N) is non-increasing in each product scaling law, while the capture
        # radius is non-decreasing in min(i,j).  Taking their separate extrema
        # and multiplying is therefore an upper bound even where their product
        # is non-monotone across N_e or the 2Rg/sigma transition.
        temperature = self.bound_temperature
        diffusivity_bound = self.arm.chain_diffusivity(
            temperature, left_range[0]
        ) + self.arm.chain_diffusivity(temperature, right_range[0])
        capture_bound = max(
            SIGMA_CONTACT,
            2.0 * math.sqrt(PS_C_R2 * PS_M0 * min(left_range[1], right_range[1]) / 6.0),
        )
        diffusion_bound = (
            4.0 * math.pi * N_A * self.spin_factor * diffusivity_bound * capture_bound
        )
        pair = _pair_type(self._site_type(left_kind), self._site_type(right_kind))
        try:
            activation = self.bound_activation_rates[pair]
        except KeyError as error:
            raise ValueError(f"no MET activation bound for pair {pair}") from error
        if (
            self.component_bound_activation_rates is not None
            and isinstance(left_kind, tuple)
            and isinstance(right_kind, tuple)
            and len(left_kind) >= 3
            and len(right_kind) >= 3
        ):
            components = tuple(sorted((str(left_kind[2]), str(right_kind[2]))))
            activation = self.component_bound_activation_rates.get(
                (components[0], components[1], pair), 0.0
            )
            if activation == 0.0:
                return 0.0
        if activation == 0.0:
            return 0.0
        return collins_kimball(activation, diffusion_bound)


@dataclass(frozen=True)
class METChannelOption:
    channel: CompiledChannel
    records: tuple[Any, ...]


@dataclass(frozen=True)
class METPopulation:
    sampler: LengthBinnedThinningSampler
    channels: Mapping[tuple[str, str], tuple[METChannelOption, ...]]
    component_channels: Mapping[
        tuple[str, str, tuple[str, str]], tuple[METChannelOption, ...]
    ]
    item_sites: Mapping[str, tuple[Site, ...]]


def _forward_rate_bound(
    channel: CompiledChannel, lower_temperature: float, upper_temperature: float
) -> float:
    """Bound one channel over a refresh interval at all analytic extrema."""
    candidates = {float(lower_temperature), float(upper_temperature)}
    if channel.continuous_forward is None:
        candidates.update(
            value
            for value in channel.rate.temperatures
            if lower_temperature < value < upper_temperature
        )
    else:
        exponent = channel.continuous_forward.temperature_exponent
        if exponent != 0.0:
            stationary = -channel.continuous_forward.activation_energy / (
                exponent * R_RMG
            )
            if lower_temperature < stationary < upper_temperature:
                candidates.add(stationary)
    return max(channel.forward_rate(value) for value in candidates)


def _record_is_met(record: Any) -> bool:
    return (
        int(_field(record, "arity", 0)) == 2
        and int(_field(record, "radical_delta", 0)) < 0
        and str(_field(record, "family", ""))
        in {"R_Recombination", "Disproportionation"}
    )


def _ortho_label(record: Any) -> str | None:
    operations = _field(record, "junction_ops", ()) or ()
    if len(operations) != 1:
        return None
    operation = operations[0]
    if operation.get("action") != "create" or not str(
        operation.get("junction_kind", "")
    ).startswith("J_ortho"):
        return None
    label = str(operation.get("attacked_atom_label", ""))
    return label if label in {"S6", "S7"} else None


def build_met_population(
    state: KMCState,
    index: SiteIndex,
    records: Sequence[Any],
    table: CompiledTerminationTable,
    *,
    temperature: float,
    volume: float,
    spin_factor: float,
    bound_temperature: float | None = None,
    _layout: Any = None,
    _rates: Any = None,
) -> METPopulation:
    """Build component-level MET bins from current indexed radical sites."""
    if table.kernel != "bulk":
        raise ValueError("SSA termination needs a bulk MET table")
    if _layout is not None:
        channels, grouped_sites, component_channels = _layout.populations(state, index)
        activation = {pair: 0.0 for pair in channels}
    else:
        channels, grouped_sites, component_channels, activation = _full_met_populations(
            state, index, records, table, temperature
        )

    forward = (
        (lambda channel: channel.forward_rate(temperature))
        if _rates is None
        else _rates.forward
    )
    maximum_temperature = (
        temperature if bound_temperature is None else bound_temperature
    )
    bound = (
        (lambda channel: _forward_rate_bound(channel, temperature, maximum_temperature))
        if _rates is None
        else (lambda channel: _rates.bound(channel, maximum_temperature))
    )
    items = []
    item_sites = {}
    for (component, site_type, site_class), sites_by_key in sorted(
        grouped_sites.items()
    ):
        item_id = f"{component}:{site_type}:{site_class}"
        sites = tuple(sites_by_key[key] for key in sorted(sites_by_key))
        item_sites[item_id] = sites
        items.append(
            PairItem(
                item_id,
                component,
                int(state.component_lengths[component]),
                (site_type, site_class),
            )
        )
    available_kinds = {item.kind for item in items}
    pair_kinds = []
    for left in sorted(available_kinds, key=repr):
        for right in sorted(available_kinds, key=repr):
            if repr(right) < repr(left):
                continue
            if _pair_type(str(left[0]), str(right[0])) in channels:
                pair_kinds.append((left, right))
    component_activation = {
        key: pairwise_sum(forward(option.channel) for option in options)
        for key, options in component_channels.items()
    }
    component_bound_activation = {
        key: pairwise_sum(bound(option.channel) for option in options)
        for key, options in component_channels.items()
    }
    bound_activation = {
        pair: max(
            (
                rate
                for (*_, component_pair), rate in component_bound_activation.items()
                if component_pair == pair
            ),
            default=0.0,
        )
        for pair in channels
    }
    kernel = METPairKernel(
        activation,
        table.arm,
        temperature,
        spin_factor=spin_factor,
        bound_temperature=bound_temperature,
        bound_activation_rates=bound_activation,
        component_activation_rates=component_activation,
        component_bound_activation_rates=component_bound_activation,
    )
    sampler = LengthBinnedThinningSampler(
        items,
        kernel,
        pair_kinds,
        volume=volume,
        normalization=N_A,
    )
    return METPopulation(
        sampler,
        {pair: tuple(options) for pair, options in channels.items()},
        {key: tuple(options) for key, options in component_channels.items()},
        item_sites,
    )


def _full_met_populations(state, index, records, table, temperature):
    """The uncached full-catalogue oracle's population discovery."""
    by_event = {str(_field(record, "event_id", "")): record for record in records}
    ortho_records = tuple(
        record for record in records if _ortho_label(record) is not None
    )
    channels: dict[tuple[str, str], list[METChannelOption]] = {}
    activation: dict[tuple[str, str], float] = {}
    for channel in table.channels:
        option_records = (
            ortho_records
            if channel.channel_id == "J_ortho"
            else (by_event[channel.event_id],)
        )
        representative = option_records[0]
        site_types = _expanded_site_types(representative)
        if len(site_types) != 2:
            raise ValueError(f"MET channel {channel.channel_id} is not bimolecular")
        pair = _pair_type(*site_types)
        channels.setdefault(pair, []).append(METChannelOption(channel, option_records))
        activation[pair] = activation.get(pair, 0.0) + channel.forward_rate(temperature)

    grouped_sites: dict[tuple[str, str, str], dict[tuple[str, ...], Site]] = {}
    for record in records:
        if (
            not _record_is_met(record)
            or _field(record, "status", "enabled") == "refused"
        ):
            continue
        for candidate in index.candidates(str(_field(record, "event_id", ""))):
            for site in candidate:
                component = state.components[site.strand_id]
                site_class = radical_site_class(site, state)
                grouped_sites.setdefault((component, site.site_type, site_class), {})[
                    _site_key(site)
                ] = site

    component_channels: dict[
        tuple[str, str, tuple[str, str]], list[METChannelOption]
    ] = {}
    for pair, options in channels.items():
        for option in options:
            component_pairs = set()
            for record in option.records:
                for candidate in eligible_candidates(state, index, record):
                    components = tuple(
                        sorted(state.components[site.strand_id] for site in candidate)
                    )
                    if len(set(components)) == 2:
                        component_pairs.add(components)
            for components in component_pairs:
                component_channels.setdefault(
                    (components[0], components[1], pair), []
                ).append(option)

    return channels, grouped_sites, component_channels, activation


def record_rate(record: Any, temperature: float) -> float:
    """Evaluate one compiled rate table with its one recorded SSA multiplier."""
    table = _field(record, "k_table")
    if not table:
        return 0.0
    return RateTable.from_mapping(table)(temperature) * float(
        _field(record, "ssa_multiplier", 1.0)
    )


class _RateCache:
    """Lazy fixed-temperature rates; a new temperature gets a fresh cache."""

    def __init__(self, temperature: float):
        self.temperature = temperature
        self.records: dict[str, float] = {}
        self.channels: dict[str, float] = {}
        self.bounds: dict[tuple[str, float], float] = {}

    def record(self, record: Any) -> float:
        key = str(_field(record, "event_id", ""))
        if key not in self.records:
            self.records[key] = record_rate(record, self.temperature)
        return self.records[key]

    def forward(self, channel: CompiledChannel) -> float:
        if channel.channel_id not in self.channels:
            self.channels[channel.channel_id] = channel.forward_rate(self.temperature)
        return self.channels[channel.channel_id]

    def bound(self, channel: CompiledChannel, maximum_temperature: float) -> float:
        key = channel.channel_id, maximum_temperature
        if key not in self.bounds:
            self.bounds[key] = (
                self.forward(channel)
                if maximum_temperature == self.temperature
                else _forward_rate_bound(channel, self.temperature, maximum_temperature)
            )
        return self.bounds[key]


def _indexed_eligible(state: KMCState, index: SiteIndex, record: Any):
    key = str(_field(record, "event_id", ""))
    if key not in index._eligible:
        index._eligible[key] = eligible_candidates(state, index, record)
    return index._eligible[key]


class _METLayout:
    """Static channel layout with changed-record site/pair contributions."""

    def __init__(self, records: Sequence[Any], table: CompiledTerminationTable):
        if table.kernel != "bulk":
            raise ValueError("SSA termination needs a bulk MET table")
        by_event = {str(_field(record, "event_id", "")): record for record in records}
        self.record_order = {key: ordinal for ordinal, key in enumerate(by_event)}
        ortho = tuple(record for record in records if _ortho_label(record) is not None)
        self.records = {
            key: record
            for key, record in by_event.items()
            if _record_is_met(record)
            and _field(record, "status", "enabled") != "refused"
        }
        self.channels: dict[tuple[str, str], list[METChannelOption]] = {}
        self.options_by_event: dict[str, list[tuple[Any, int, METChannelOption]]] = {}
        self.met_event_ids: set[str] = set()
        self.channel_ids: set[str] = set()
        for ordinal, channel in enumerate(table.channels):
            option_records = (
                ortho
                if channel.channel_id == "J_ortho"
                else (by_event[channel.event_id],)
            )
            types = _expanded_site_types(option_records[0])
            if len(types) != 2:
                raise ValueError(f"MET channel {channel.channel_id} is not bimolecular")
            pair = _pair_type(*types)
            option = METChannelOption(channel, option_records)
            self.channels.setdefault(pair, []).append(option)
            self.channel_ids.add(channel.channel_id)
            for record in option_records:
                key = str(_field(record, "event_id", ""))
                self.met_event_ids.add(key)
                self.options_by_event.setdefault(key, []).append(
                    (pair, ordinal, option)
                )
                self.records[key] = record
        self.contributions: dict[str, tuple[Any, Any]] = {}
        self.revision = -1

    def populations(self, state: KMCState, index: SiteIndex):
        changed = (
            index.changed_event_ids
            if self.revision == index.revision - 1
            else self.records.keys()
        )
        for key in changed:
            record = self.records.get(key)
            if record is None:
                continue
            sites = {}
            if (
                _record_is_met(record)
                and _field(record, "status", "enabled") != "refused"
            ):
                for candidate in index.candidates(key):
                    for site in candidate:
                        component = state.components[site.strand_id]
                        group = (
                            component,
                            site.site_type,
                            radical_site_class(site, state),
                        )
                        sites.setdefault(group, {})[_site_key(site)] = site
            pairs = (
                {
                    tuple(
                        sorted(state.components[site.strand_id] for site in candidate)
                    )
                    for candidate in _indexed_eligible(state, index, record)
                    if len({state.components[site.strand_id] for site in candidate})
                    == 2
                }
                if key in self.options_by_event
                else set()
            )
            if sites or pairs:
                self.contributions[key] = sites, pairs
            else:
                self.contributions.pop(key, None)
        self.revision = index.revision
        grouped_sites = {}
        # Preserve the full builder's overwrite order for equivalent sites.
        for key in sorted(self.contributions, key=self.record_order.__getitem__):
            sites, _ = self.contributions[key]
            for group, entries in sites.items():
                grouped_sites.setdefault(group, {}).update(entries)
        indexed_options = {}
        for key, (_, pairs) in self.contributions.items():
            for pair, ordinal, option in self.options_by_event.get(key, ()):
                for left, right in pairs:
                    indexed_options.setdefault((left, right, pair), {})[
                        ordinal
                    ] = option
        component_channels = {
            key: [options[ordinal] for ordinal in sorted(options)]
            for key, options in indexed_options.items()
        }
        return self.channels, grouped_sites, component_channels


def _candidate_key(candidate: Sequence[Site]) -> tuple[tuple[str, ...], ...]:
    return tuple(_site_key(site) for site in candidate)


def eligible_candidates(
    state: KMCState, index: SiteIndex, record: Any
) -> tuple[tuple[Site, ...], ...]:
    """Return executor candidates after §M1 and same-component filtering."""
    candidates = index.candidates(str(_field(record, "event_id", "")))
    if int(_field(record, "arity", 0)) != 2:
        return candidates
    identical = (
        _field(record, "reactant_pair_convention", "")
        == "unordered-identical-pair N(N-1)/2"
    )
    found = {}
    for candidate in candidates:
        if (
            state.components[candidate[0].strand_id]
            == state.components[candidate[1].strand_id]
        ):
            continue
        key = _candidate_key(candidate)
        if identical and key != min(key, tuple(reversed(key))):
            continue
        found[key] = candidate
    return tuple(found[key] for key in sorted(found))


def record_propensity(
    state: KMCState,
    index: SiteIndex,
    record: Any,
    temperature: float,
    volume: float,
) -> float:
    """Return a non-MET record propensity under the compiled convention."""
    candidates = eligible_candidates(state, index, record)
    rate = record_rate(record, temperature)
    arity = int(_field(record, "arity", 0))
    if arity == 1:
        return rate * len(candidates)
    if arity == 2:
        if not math.isfinite(volume) or volume <= 0.0:
            raise ValueError("volume must be finite and positive")
        return rate * len(candidates) / (N_A * volume)
    raise ValueError(f"SSA supports only arity one or two, got {arity}")


@dataclass(frozen=True)
class ChannelPropensityReport:
    enabled_channel_ids: frozenset[str]
    propensities: Mapping[str, float]
    refused_propensities: Mapping[str, float]
    total_enabled: float
    total_refused: float


def _met_channel_propensities(
    population: METPopulation, temperature: float, _rates: _RateCache | None = None
) -> dict[str, float]:
    result: dict[str, float] = {}
    for cell in population.sampler.cells:
        pair = _pair_type(str(cell.left_kind[0]), str(cell.right_kind[0]))
        for left, right in population.sampler.cell_pairs(cell):
            components = tuple(sorted((left.component_id, right.component_id)))
            options = population.component_channels.get(
                (components[0], components[1], pair), ()
            )
            if not options:
                continue
            option_rates = [
                (
                    option.channel.forward_rate(temperature)
                    if _rates is None
                    else _rates.forward(option.channel)
                )
                for option in options
            ]
            activation = pairwise_sum(option_rates)
            exact = population.sampler.kernel.exact_rate(left, right) / (
                population.sampler.normalization * population.sampler.volume
            )
            for option, rate in zip(options, option_rates):
                result[option.channel.channel_id] = (
                    result.get(option.channel.channel_id, 0.0)
                    + exact * rate / activation
                )
    return result


def channel_propensities(
    state: KMCState,
    index: SiteIndex,
    records: Sequence[Any],
    temperature: float,
    volume: float,
    *,
    met_population: METPopulation | None = None,
) -> ChannelPropensityReport:
    """Expose canonical enabled channel IDs and per-channel propensities."""
    enabled: dict[str, float] = {}
    refused: dict[str, float] = {}
    met_event_ids = set()
    if met_population is not None:
        enabled.update(_met_channel_propensities(met_population, temperature))
        met_event_ids = {
            str(_field(record, "event_id", ""))
            for options in met_population.channels.values()
            for option in options
            for record in option.records
        }
    for record in records:
        event_id = str(_field(record, "event_id", ""))
        if event_id in met_event_ids:
            continue
        value = record_propensity(state, index, record, temperature, volume)
        target = (
            refused if _field(record, "status", "enabled") == "refused" else enabled
        )
        target[event_id] = value
    enabled_vector = canonical_propensities(enabled.items())
    refused_vector = canonical_propensities(refused.items())
    return ChannelPropensityReport(
        frozenset(enabled),
        dict(zip(enabled_vector.event_ids, enabled_vector.values)),
        dict(zip(refused_vector.event_ids, refused_vector.values)),
        enabled_vector.total,
        refused_vector.total,
    )


@dataclass(frozen=True)
class SSAEvent:
    time: float
    event_id: str
    state_hash: str
    irreversible: bool


@dataclass(frozen=True)
class SSARunResult:
    time: float
    events: tuple[SSAEvent, ...]
    leak_num: float
    leak_den: float
    leak: float
    irreversible_fraction: float
    bound_checks: int
    bound_violations: int


class _IncrementalPropensities:
    """Per-engine maintenance, independent even when engines share an index."""

    def __init__(self, engine: Any):
        self.engine = engine
        self.rates = _RateCache(engine.temperature)
        self.layout = (
            _METLayout(engine.records, engine.termination_table)
            if engine.termination_table is not None
            else None
        )
        met_ids = self.layout.met_event_ids if self.layout else set()
        self.records = {
            str(_field(record, "event_id", "")): record
            for record in engine.records
            if str(_field(record, "event_id", "")) not in met_ids
        }
        self.enabled = {
            key: 0.0
            for key, record in self.records.items()
            if _field(record, "status", "enabled") != "refused"
        }
        self.refused = {key: 0.0 for key in self.records if key not in self.enabled}
        self.refused_tree = _PropensityTree(self.refused)
        sampling = {
            key: value
            for key, value in self.enabled.items()
            if self.layout is None or key not in self.layout.channel_ids
        }
        if self.layout is not None:
            sampling[engine._MET_LEAF] = 0.0
        self.sampling_tree = _PropensityTree(sampling)
        self.enabled_trees: OrderedDict[tuple[str, ...], _PropensityTree] = (
            OrderedDict()
        )
        self.enabled_tree = _PropensityTree(self.enabled)
        self.population: METPopulation | None = None
        self.met_values: dict[str, float] = {}
        self.revision = -1

    def sync(self) -> None:
        engine = self.engine
        index = engine.index
        if self.revision == index.revision:
            return
        consecutive = self.revision == index.revision - 1
        changed = index.changed_event_ids if consecutive else self.records.keys()
        updates = {}
        for key in changed:
            record = self.records.get(key)
            if record is None:
                continue
            arity = int(_field(record, "arity", 0))
            if arity not in {1, 2}:
                raise ValueError(f"SSA supports only arity one or two, got {arity}")
            candidates = _indexed_eligible(engine.state, index, record)
            # Ineligible records do not construct or interpolate rate tables.
            value = self.rates.record(record) * len(candidates) if candidates else 0.0
            if arity == 2:
                value /= N_A * engine.volume
            if key in self.refused:
                self.refused[key] = value
                self.refused_tree.update(key, value)
            else:
                self.enabled[key] = value
                updates[key] = value
                if key in self.sampling_tree.positions:
                    self.sampling_tree.update(key, value)
        met_changed = self.layout is not None and (
            not consecutive
            or bool(index.changed_event_ids & self.layout.records.keys())
        )
        if met_changed:
            self.population = build_met_population(
                engine.state,
                index,
                engine.records,
                engine.termination_table,
                temperature=engine.temperature,
                volume=engine.volume,
                spin_factor=engine.spin_factor,
                bound_temperature=engine.bound_temperature,
                _layout=self.layout,
                _rates=self.rates,
            )
            self.met_values = _met_channel_propensities(
                self.population, engine.temperature, self.rates
            )
            self.sampling_tree.update(
                engine._MET_LEAF, self.population.sampler.total_bound_propensity
            )
        # MET report leaves exist only when eligible. Cache each exact tree
        # shape; zero padding or dropping other zero leaves would change rounding.
        shape = tuple(sorted(self.met_values))
        tree = self.enabled_trees.get(shape)
        for cached in self.enabled_trees.values():
            for key, value in updates.items():
                cached.update(key, value)
        if tree is None:
            tree = _PropensityTree({**self.enabled, **self.met_values})
            self.enabled_trees[shape] = tree
        else:
            for key, value in self.met_values.items():
                tree.update(key, value)
        self.enabled_trees.move_to_end(shape)
        while len(self.enabled_trees) > 8:
            self.enabled_trees.popitem(last=False)
        self.enabled_tree = tree
        self.revision = index.revision

    def report(self) -> ChannelPropensityReport:
        self.sync()
        values = {**self.enabled, **self.met_values}
        return ChannelPropensityReport(
            frozenset(values),
            {key: values[key] for key in self.enabled_tree.event_ids},
            {key: self.refused[key] for key in self.refused_tree.event_ids},
            self.enabled_tree.total,
            self.refused_tree.total,
        )


class IsothermalSSA:
    """Direct-method isothermal SSA with MET thinning and exact leak clocks."""

    _MET_LEAF = "~MET-THINNING-BOUND"

    def __init__(
        self,
        state: KMCState,
        records: Sequence[Any],
        *,
        temperature: float,
        volume: float,
        rng: Any,
        termination_table: CompiledTerminationTable | None = None,
        spin_factor: float = 1.0,
        bound_temperature: float | None = None,
        incremental: bool = True,
    ):
        self.state = state
        self.records = tuple(records)
        self.temperature = float(temperature)
        self.volume = float(volume)
        if not math.isfinite(self.volume) or self.volume <= 0.0:
            raise ValueError("volume must be finite and positive")
        self.rng = rng
        self.termination_table = termination_table
        self.spin_factor = float(spin_factor)
        self.bound_temperature = bound_temperature
        self.incremental = bool(incremental)
        self._maintenance: _IncrementalPropensities | None = None
        self._maintenance_configuration: tuple[Any, ...] | None = None
        self.index = SiteIndex(state, self.records)
        self.time = 0.0
        self.leak_num = 0.0
        self.leak_den = 0.0
        self.events: list[SSAEvent] = []
        self.irreversible_fired = 0
        self.bound_checks = 0
        self.bound_violations = 0
        self._records_by_id = {
            str(_field(record, "event_id", "")): record for record in self.records
        }

    def _met_population(self) -> METPopulation | None:
        if self.incremental:
            maintenance = self._incremental_propensities()
            maintenance.sync()
            return maintenance.population
        if self.termination_table is None:
            return None
        return build_met_population(
            self.state,
            self.index,
            self.records,
            self.termination_table,
            temperature=self.temperature,
            volume=self.volume,
            spin_factor=self.spin_factor,
            bound_temperature=self.bound_temperature,
        )

    def channel_propensities(self) -> ChannelPropensityReport:
        if self.incremental:
            return self._incremental_propensities().report()
        return channel_propensities(
            self.state,
            self.index,
            self.records,
            self.temperature,
            self.volume,
            met_population=self._met_population(),
        )

    def _incremental_propensities(self) -> _IncrementalPropensities:
        configuration = (
            id(self.index),
            id(self.records),
            id(self.termination_table),
            self.temperature,
            self.volume,
            self.spin_factor,
            self.bound_temperature,
        )
        if configuration != self._maintenance_configuration or (
            self._maintenance is not None and self._maintenance.engine is not self
        ):
            self._maintenance = _IncrementalPropensities(self)
            self._maintenance_configuration = configuration
        assert self._maintenance is not None
        return self._maintenance

    def _integrate(self, elapsed: float, enabled: float, refused: float) -> None:
        self.leak_num += refused * elapsed
        self.leak_den += (enabled + refused) * elapsed
        self.time += elapsed

    def _record_event(self, record: Any) -> SSAEvent:
        irreversible = _field(record, "status", "enabled") == "irreversible"
        if irreversible:
            self.irreversible_fired += 1
        event = SSAEvent(
            self.time,
            str(_field(record, "event_id", "")),
            self.state.state_hash(),
            irreversible,
        )
        self.events.append(event)
        return event

    def _fire_non_met(self, event_id: str) -> SSAEvent:
        record = self._records_by_id[event_id]
        candidates = (
            _indexed_eligible(self.state, self.index, record)
            if self.incremental
            else eligible_candidates(self.state, self.index, record)
        )
        if not candidates:
            raise AssertionError("positive record propensity has no candidate")
        participants = candidates[int(self.rng.integers(len(candidates)))]
        self.index.apply(record, participants)
        return self._record_event(record)

    def _choose_met_record(
        self, population: METPopulation, proposal: PairProposal
    ) -> tuple[Any, tuple[Site, ...]]:
        pair = _pair_type(str(proposal.left.kind[0]), str(proposal.right.kind[0]))
        components = tuple(
            sorted((proposal.left.component_id, proposal.right.component_id))
        )
        options = population.component_channels[(components[0], components[1], pair)]
        vector = canonical_propensities(
            (
                option.channel.channel_id,
                (
                    self._incremental_propensities().rates.forward(option.channel)
                    if self.incremental
                    else option.channel.forward_rate(self.temperature)
                ),
            )
            for option in options
        )
        channel_id = vector.draw(self.rng)
        option = next(
            option for option in options if option.channel.channel_id == channel_id
        )
        if option.channel.channel_id == "J_ortho":
            label = draw_ortho_attack_site(self.rng)
            record = next(
                record for record in option.records if _ortho_label(record) == label
            )
        else:
            record = option.records[0]
        component_set = {proposal.left.component_id, proposal.right.component_id}
        site_keys = {
            *(_site_key(site) for site in population.item_sites[proposal.left.item_id]),
            *(
                _site_key(site)
                for site in population.item_sites[proposal.right.item_id]
            ),
        }
        candidates = []
        for candidate in eligible_candidates(self.state, self.index, record):
            candidate_components = {
                self.state.components[site.strand_id] for site in candidate
            }
            if candidate_components == component_set and all(
                _site_key(site) in site_keys for site in candidate
            ):
                candidates.append(candidate)
        if not candidates:
            raise AssertionError("accepted MET component pair has no record mapping")
        return record, candidates[int(self.rng.integers(len(candidates)))]

    def step(self, *, until: float | None = None) -> SSAEvent | None:
        """Advance through null proposals until one event fires or a horizon ends."""
        while until is None or self.time < until:
            if self.incremental:
                maintenance = self._incremental_propensities()
                maintenance.sync()
                population = maintenance.population
                vector = maintenance.sampling_tree
                enabled = maintenance.enabled_tree.total
                refused = maintenance.refused_tree.total
            else:
                population = self._met_population()
                report = channel_propensities(
                    self.state,
                    self.index,
                    self.records,
                    self.temperature,
                    self.volume,
                    met_population=population,
                )
                met_ids = (
                    {
                        option.channel.channel_id
                        for options in population.channels.values()
                        for option in options
                    }
                    if population is not None
                    else set()
                )
                leaves = [
                    (event_id, value)
                    for event_id, value in report.propensities.items()
                    if event_id not in met_ids
                ]
                if population is not None:
                    leaves.append(
                        (self._MET_LEAF, population.sampler.total_bound_propensity)
                    )
                vector = canonical_propensities(leaves)
                enabled, refused = report.total_enabled, report.total_refused
            if vector.absorbing:
                if until is not None:
                    self._integrate(
                        until - self.time,
                        enabled,
                        refused,
                    )
                return None
            elapsed = float(self.rng.exponential(1.0 / vector.total))
            if until is not None and self.time + elapsed > until:
                self._integrate(until - self.time, enabled, refused)
                return None
            self._integrate(elapsed, enabled, refused)
            selected = vector.draw(self.rng)
            if selected != self._MET_LEAF:
                assert selected is not None
                return self._fire_non_met(selected)
            assert population is not None
            checks, violations = (
                population.sampler.bound_checks,
                population.sampler.bound_violations,
            )
            proposal = population.sampler.select_pair(self.rng)
            self.bound_checks += population.sampler.bound_checks - checks
            self.bound_violations += population.sampler.bound_violations - violations
            if proposal is None or not proposal.accepted:
                continue
            record, participants = self._choose_met_record(population, proposal)
            self.index.apply(record, participants)
            return self._record_event(record)
        return None

    def run(
        self, *, until: float | None = None, max_events: int | None = None
    ) -> SSARunResult:
        if until is None and max_events is None:
            raise ValueError("run needs a time or event limit")
        start_events = len(self.events)
        while (until is None or self.time < until) and (
            max_events is None or len(self.events) - start_events < max_events
        ):
            event = self.step(until=until)
            if event is None:
                break
        leak = self.leak_num / self.leak_den if self.leak_den else 0.0
        irreversible_fraction = (
            self.irreversible_fired / len(self.events) if self.events else 0.0
        )
        return SSARunResult(
            self.time,
            tuple(self.events),
            self.leak_num,
            self.leak_den,
            leak,
            irreversible_fraction,
            self.bound_checks,
            self.bound_violations,
        )


__all__ = [
    "CanonicalPropensities",
    "ChannelPropensityReport",
    "ConstantPairKernel",
    "E_PROPENSITY_NEGATIVE",
    "E_PROPENSITY_NONFINITE",
    "PropensityError",
    "LengthBinnedThinningSampler",
    "IsothermalSSA",
    "METChannelOption",
    "METPairKernel",
    "METPopulation",
    "PairCell",
    "PairItem",
    "PairProposal",
    "SSAEvent",
    "SSARunResult",
    "SiteIndex",
    "SumPairKernel",
    "ThinningBoundError",
    "build_met_population",
    "canonical_propensities",
    "channel_propensities",
    "eligible_candidates",
    "epsilon_propensity",
    "length_bin",
    "pairwise_sum",
    "radical_pair_class",
    "radical_site_class",
    "record_propensity",
    "record_rate",
]
