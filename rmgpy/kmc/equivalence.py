"""Semantic comparison of compiled kMC event-set artifacts.

Event identifiers deliberately include compiler provenance.  This module instead
compares the rooted, atom-mapped reaction graph represented by each record and
then compares the physical-channel data attached to that graph.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from functools import lru_cache
import json
import math
from pathlib import Path
import sys
import time
from typing import Any, Iterable

import networkx as nx


class EquivalenceError(ValueError):
    """An artifact lacks information required for semantic comparison."""


@dataclass(frozen=True)
class _ParsedMolecule:
    atom_labels: tuple[str, ...]
    bonds: tuple[tuple[int, int, str], ...]
    heavy_atoms: tuple[int, ...]
    formula: tuple[tuple[str, int], ...]
    multiplicity: int


@dataclass(frozen=True)
class _CompactGraph:
    nodes: tuple[str, ...]
    edges: tuple[tuple[int, int, str], ...]
    fingerprint: str

    def materialize(self) -> nx.Graph:
        graph = nx.Graph()
        graph.add_nodes_from(
            (index, {"label": label}) for index, label in enumerate(self.nodes)
        )
        graph.add_edges_from(
            (first, second, {"label": label}) for first, second, label in self.edges
        )
        return graph


def _canonical_label(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


@lru_cache(maxsize=512)
def _parse_molecule(adjacency: str) -> _ParsedMolecule:
    from rmgpy.molecule.molecule import Molecule

    molecule = Molecule().from_adjacency_list(adjacency)
    atom_index = {atom: index for index, atom in enumerate(molecule.atoms)}
    atom_labels = tuple(
        _canonical_label(
            {
                "element": atom.element.symbol,
                "isotope": atom.element.isotope,
                "radicals": int(atom.radical_electrons),
                "lone_pairs": int(atom.lone_pairs),
                "charge": int(atom.charge),
                "implicit_hydrogens": int(getattr(atom, "implicit_hydrogens", 0)),
            }
        )
        for atom in molecule.atoms
    )
    bonds = tuple(
        sorted(
            (
                min(atom_index[bond.vertex1], atom_index[bond.vertex2]),
                max(atom_index[bond.vertex1], atom_index[bond.vertex2]),
                _canonical_label({"kind": "bond", "order": str(bond.order)}),
            )
            for bond in molecule.get_all_edges()
        )
    )
    return _ParsedMolecule(
        atom_labels=atom_labels,
        bonds=bonds,
        heavy_atoms=tuple(
            index
            for index, atom in enumerate(molecule.atoms)
            if atom.element.number != 1
        ),
        formula=tuple(
            sorted(
                (key, int(value)) for key, value in molecule.get_element_count().items()
            )
        ),
        multiplicity=int(molecule.multiplicity),
    )


def _compact_graph(
    nodes: list[str], edges: list[tuple[int, int, str]]
) -> _CompactGraph:
    graph = nx.Graph()
    graph.add_nodes_from((index, {"label": label}) for index, label in enumerate(nodes))
    graph.add_edges_from(
        (first, second, {"label": label}) for first, second, label in edges
    )
    fingerprint = nx.weisfeiler_lehman_graph_hash(
        graph,
        node_attr="label",
        edge_attr="label",
        iterations=4,
        digest_size=32,
    )
    return _CompactGraph(tuple(nodes), tuple(edges), fingerprint)


def _isomorphic(first: _CompactGraph, second: _CompactGraph) -> bool:
    if first == second:
        return True
    if first.fingerprint != second.fingerprint:
        return False
    if len(first.nodes) != len(second.nodes) or len(first.edges) != len(second.edges):
        return False
    return nx.is_isomorphic(
        first.materialize(),
        second.materialize(),
        node_match=lambda left, right: left["label"] == right["label"],
        edge_match=lambda left, right: left["label"] == right["label"],
    )


def _coproduct_formulas(record: dict[str, Any]) -> Counter:
    formulas = Counter()
    for coproduct in record.get("coproducts", []):
        formula = coproduct.get("formula") if isinstance(coproduct, dict) else None
        if not isinstance(formula, dict):
            raise EquivalenceError("coproduct entries must contain a formula")
        formulas[
            tuple(sorted((key, int(value)) for key, value in formula.items()))
        ] += 1
    return formulas


def _product_fates(
    record: dict[str, Any], products: list[_ParsedMolecule]
) -> list[str]:
    """Infer permutation-stable melt/leave labels from coproduct formula counts.

    The current schema stores coproduct formulae rather than graph handles.  A
    partial formula match is therefore resolved by the writer's convention that
    the first product stays in the melt and later products are coproducts.
    """
    remaining = _coproduct_formulas(record)
    fates = ["melt"] * len(products)
    by_formula: dict[tuple[tuple[str, int], ...], list[int]] = defaultdict(list)
    for index, molecule in enumerate(products):
        by_formula[molecule.formula].append(index)
    for formula, indices in by_formula.items():
        count = remaining.pop(formula, 0)
        if count > len(indices):
            raise EquivalenceError("coproduct formula count exceeds matching products")
        for index in indices[-count:] if count else ():
            fates[index] = "leave"
    if remaining:
        raise EquivalenceError("coproduct formula does not match a product graph")
    return fates


def _reaction_graphs(record: dict[str, Any]) -> tuple[_CompactGraph, _CompactGraph]:
    reactants = [_parse_molecule(item) for item in record.get("reactant_graphs", [])]
    products = [_parse_molecule(item) for item in record.get("product_graphs", [])]
    if not reactants or not products:
        raise EquivalenceError("records require reactant_graphs and product_graphs")

    nodes: list[str] = []
    edges: list[tuple[int, int, str]] = []
    reactant_atoms: list[int] = []
    reactant_heavy: list[int] = []
    product_heavy: list[int] = []

    def add_side(side: str, molecules: list[_ParsedMolecule]) -> None:
        for molecule in molecules:
            participant = len(nodes)
            nodes.append(
                _canonical_label(
                    {
                        "kind": "participant",
                        "side": side,
                        "multiplicity": molecule.multiplicity,
                    }
                )
            )
            offset = len(nodes)
            nodes.extend(
                _canonical_label({"kind": "atom", "side": side, "state": atom_label})
                for atom_label in molecule.atom_labels
            )
            for local_index in range(len(molecule.atom_labels)):
                atom_node = offset + local_index
                edges.append(
                    (participant, atom_node, _canonical_label({"kind": "membership"}))
                )
                if side == "reactant":
                    reactant_atoms.append(atom_node)
            for first, second, bond_label in molecule.bonds:
                edges.append((offset + first, offset + second, bond_label))
            heavy = [offset + local_index for local_index in molecule.heavy_atoms]
            (reactant_heavy if side == "reactant" else product_heavy).extend(heavy)

    add_side("reactant", reactants)
    add_side("product", products)

    atom_map = {
        int(key): int(value) for key, value in record.get("atom_map", {}).items()
    }
    if set(atom_map) != set(range(len(reactant_heavy))):
        raise EquivalenceError("atom_map reactant indices must cover all heavy atoms")
    if set(atom_map.values()) != set(range(len(product_heavy))):
        raise EquivalenceError("atom_map product indices must cover all heavy atoms")
    for reactant_index, product_index in atom_map.items():
        edges.append(
            (
                reactant_heavy[reactant_index],
                product_heavy[product_index],
                _canonical_label({"kind": "atom_map"}),
            )
        )

    for operation in record.get("bond_ops", []):
        if not isinstance(operation, dict):
            raise EquivalenceError("bond_ops entries must be objects")
        operation_data = {
            key: value
            for key, value in operation.items()
            if key not in {"atom", "atoms"}
        }
        operation_node = len(nodes)
        nodes.append(
            _canonical_label({"kind": "rewrite_operation", "data": operation_data})
        )
        references: Iterable[Any]
        references = operation.get("atoms", [operation.get("atom")])
        for reference in references:
            if not isinstance(reference, int) or not 0 <= reference < len(
                reactant_atoms
            ):
                raise EquivalenceError(
                    "bond operation references an unknown reactant atom"
                )
            edges.append(
                (
                    operation_node,
                    reactant_atoms[reference],
                    _canonical_label({"kind": "rewrite_root"}),
                )
            )

    chemistry = _compact_graph(nodes, edges)
    fate_nodes: list[str] = []
    fate_edges: list[tuple[int, int, str]] = []
    for molecule, product_fate in zip(products, _product_fates(record, products)):
        participant = len(fate_nodes)
        fate_nodes.append(
            _canonical_label(
                {
                    "kind": "participant",
                    "multiplicity": molecule.multiplicity,
                    "fate": product_fate,
                }
            )
        )
        offset = len(fate_nodes)
        fate_nodes.extend(molecule.atom_labels)
        fate_edges.extend(
            (
                participant,
                offset + local_index,
                _canonical_label({"kind": "membership"}),
            )
            for local_index in range(len(molecule.atom_labels))
        )
        fate_edges.extend(
            (offset + first, offset + second, bond_label)
            for first, second, bond_label in molecule.bonds
        )
    fate = _compact_graph(fate_nodes, fate_edges)
    return chemistry, fate


@dataclass
class _Channel:
    index: int
    chemistry: _CompactGraph
    record_count: int = 0
    event_ids: list[str] = field(default_factory=list)
    reverse_ids: list[str] = field(default_factory=list)
    reverse_channels: set[int] = field(default_factory=set)
    unresolved_reverse_ids: set[str] = field(default_factory=set)
    nonreciprocal_reverse: bool = False
    statuses: set[str] = field(default_factory=set)
    rate_orders: set[int] = field(default_factory=set)
    rate_units: set[str] = field(default_factory=set)
    pair_conventions: set[str] = field(default_factory=set)
    ssa_multipliers: set[float] = field(default_factory=set)
    degeneracy: float = 0.0
    raw_path_degeneracy: float = 0.0
    symmetry_degeneracy: float = 0.0
    rates: dict[tuple[float, ...], list[float]] = field(default_factory=dict)
    propensities: dict[tuple[float, ...], list[float]] = field(default_factory=dict)
    has_missing_rate: bool = False
    fate_variants: dict[str, list[_CompactGraph]] = field(default_factory=dict)

    def add(self, record: dict[str, Any], fate: _CompactGraph) -> None:
        self.record_count += 1
        event_id = record.get("event_id")
        if isinstance(event_id, str):
            self.event_ids.append(event_id)
        reverse_id = record.get("reverse_of")
        if isinstance(reverse_id, str):
            self.reverse_ids.append(reverse_id)
        self.statuses.add(str(record.get("status", "")))
        self.rate_orders.add(int(record.get("rate_order", 0)))
        self.rate_units.add(str(record.get("rate_units", "")))
        self.pair_conventions.add(str(record.get("reactant_pair_convention", "")))
        ssa = float(record.get("ssa_multiplier", 1.0))
        degeneracy = float(record.get("degeneracy", 1.0))
        raw = float(record.get("raw_path_degeneracy", degeneracy))
        self.ssa_multipliers.add(ssa)
        self.degeneracy += degeneracy
        self.raw_path_degeneracy += raw
        self.symmetry_degeneracy += degeneracy * ssa

        variants = self.fate_variants.setdefault(fate.fingerprint, [])
        if not any(_isomorphic(fate, candidate) for candidate in variants):
            variants.append(fate)

        table = record.get("k_table")
        if table is None:
            self.has_missing_rate = True
            return
        temperatures = tuple(float(value) for value in table.get("T", []))
        rates = [float(value) for value in table.get("k", [])]
        if not temperatures or len(temperatures) != len(rates):
            raise EquivalenceError("k_table must contain aligned T and k arrays")
        rate_sum = self.rates.setdefault(temperatures, [0.0] * len(rates))
        propensity_sum = self.propensities.setdefault(temperatures, [0.0] * len(rates))
        for index, value in enumerate(rates):
            rate_sum[index] += value
            propensity_sum[index] += value * ssa


@dataclass
class _ArtifactSummary:
    temperature_grid: tuple[float, ...]
    record_count: int
    channels: list[_Channel]
    by_fingerprint: dict[str, list[int]]


def _graph_cache_key(record: dict[str, Any]) -> tuple[Any, ...]:
    """Return an exact-representation key for graph and fate construction."""
    return (
        tuple(record.get("reactant_graphs", [])),
        tuple(record.get("product_graphs", [])),
        tuple(
            sorted(
                (int(key), int(value))
                for key, value in record.get("atom_map", {}).items()
            )
        ),
        _canonical_label(record.get("bond_ops", [])),
        _canonical_label(record.get("coproducts", [])),
    )


def _summarize(
    artifact: dict[str, Any],
    graph_cache: (
        dict[tuple[Any, ...], tuple[_CompactGraph, _CompactGraph]] | None
    ) = None,
) -> _ArtifactSummary:
    if artifact.get("schema_version") != "kmc_event_set/0.1":
        raise EquivalenceError("unknown event-set schema")
    grid = tuple(
        float(value) for value in artifact.get("inputs", {}).get("temperature_grid", [])
    )
    channels: list[_Channel] = []
    by_fingerprint: dict[str, list[int]] = defaultdict(list)
    event_to_channel: dict[str, int] = {}
    records = artifact.get("records")
    if not isinstance(records, list):
        raise EquivalenceError("artifact records must be a list")

    for record in records:
        cache_key = _graph_cache_key(record)
        cached = graph_cache.get(cache_key) if graph_cache is not None else None
        if cached is None:
            chemistry, fate = _reaction_graphs(record)
            if graph_cache is not None:
                graph_cache[cache_key] = chemistry, fate
        else:
            chemistry, fate = cached
        channel = None
        for candidate_index in by_fingerprint[chemistry.fingerprint]:
            candidate = channels[candidate_index]
            if _isomorphic(chemistry, candidate.chemistry):
                channel = candidate
                break
        if channel is None:
            channel = _Channel(len(channels), chemistry)
            channels.append(channel)
            by_fingerprint[chemistry.fingerprint].append(channel.index)
        channel.add(record, fate)
        event_id = record.get("event_id")
        if isinstance(event_id, str):
            event_to_channel[event_id] = channel.index

    for channel in channels:
        for reverse_id in channel.reverse_ids:
            target = event_to_channel.get(reverse_id)
            if target is None:
                channel.unresolved_reverse_ids.add(reverse_id)
            else:
                channel.reverse_channels.add(target)
    for channel in channels:
        for target in channel.reverse_channels:
            if channel.index not in channels[target].reverse_channels:
                channel.nonreciprocal_reverse = True

    return _ArtifactSummary(grid, len(records), channels, dict(by_fingerprint))


def _load_artifact(value: str | Path | dict[str, Any]) -> dict[str, Any]:
    if isinstance(value, dict):
        return value
    with Path(value).open(encoding="utf-8") as handle:
        return json.load(handle)


def _close(first: float, second: float, rtol: float) -> bool:
    return math.isclose(first, second, rel_tol=rtol, abs_tol=0.0)


def _sequences_close(first: list[float], second: list[float], rtol: float) -> bool:
    return len(first) == len(second) and all(
        _close(left, right, rtol) for left, right in zip(first, second)
    )


def _variant_sets_equal(
    first: dict[str, list[_CompactGraph]], second: dict[str, list[_CompactGraph]]
) -> bool:
    if set(first) != set(second):
        return False
    for fingerprint in first:
        unmatched = list(second[fingerprint])
        for graph in first[fingerprint]:
            for index, candidate in enumerate(unmatched):
                if _isomorphic(graph, candidate):
                    del unmatched[index]
                    break
            else:
                return False
        if unmatched:
            return False
    return True


def _channel_ref(channel: _Channel) -> dict[str, Any]:
    return {
        "channel": channel.chemistry.fingerprint,
        "records": channel.record_count,
    }


def _compare_channel_data(
    left: _Channel,
    right: _Channel,
    left_summary: _ArtifactSummary,
    right_summary: _ArtifactSummary,
    channel_mapping: dict[int, int],
    rtol: float,
) -> list[dict[str, Any]]:
    if left is right and left_summary is right_summary:
        return []
    reasons: list[dict[str, Any]] = []

    def mismatch(reason: str, left_value: Any, right_value: Any) -> None:
        reasons.append({"reason": reason, "left": left_value, "right": right_value})

    if left_summary.temperature_grid != right_summary.temperature_grid:
        mismatch(
            "artifact temperature grid",
            left_summary.temperature_grid,
            right_summary.temperature_grid,
        )
    if left.statuses != right.statuses:
        mismatch("reversibility/status", sorted(left.statuses), sorted(right.statuses))
    if left.rate_orders != right.rate_orders:
        mismatch("rate order", sorted(left.rate_orders), sorted(right.rate_orders))
    if left.rate_units != right.rate_units:
        mismatch("rate units", sorted(left.rate_units), sorted(right.rate_units))
    if left.pair_conventions != right.pair_conventions:
        mismatch(
            "symmetry pair convention",
            sorted(left.pair_conventions),
            sorted(right.pair_conventions),
        )
    if left.ssa_multipliers != right.ssa_multipliers:
        mismatch(
            "SSA symmetry multiplier",
            sorted(left.ssa_multipliers),
            sorted(right.ssa_multipliers),
        )
    if not _close(left.degeneracy, right.degeneracy, rtol):
        mismatch("degeneracy", left.degeneracy, right.degeneracy)
    if not _close(left.raw_path_degeneracy, right.raw_path_degeneracy, rtol):
        mismatch(
            "raw path degeneracy", left.raw_path_degeneracy, right.raw_path_degeneracy
        )
    if not _close(left.symmetry_degeneracy, right.symmetry_degeneracy, rtol):
        mismatch(
            "symmetry-adjusted degeneracy",
            left.symmetry_degeneracy,
            right.symmetry_degeneracy,
        )
    if left.has_missing_rate != right.has_missing_rate:
        mismatch("rate availability", left.has_missing_rate, right.has_missing_rate)
    if set(left.rates) != set(right.rates):
        mismatch("rate temperature grid", sorted(left.rates), sorted(right.rates))
    else:
        for grid in sorted(left.rates):
            if not _sequences_close(left.rates[grid], right.rates[grid], rtol):
                mismatch("summed channel rate", left.rates[grid], right.rates[grid])
            if not _sequences_close(
                left.propensities[grid], right.propensities[grid], rtol
            ):
                mismatch(
                    "summed physical propensity",
                    left.propensities[grid],
                    right.propensities[grid],
                )
    if not _variant_sets_equal(left.fate_variants, right.fate_variants):
        mismatch(
            "product fate", sorted(left.fate_variants), sorted(right.fate_variants)
        )

    mapped_reverse = {
        channel_mapping[target]
        for target in left.reverse_channels
        if target in channel_mapping
    }
    if mapped_reverse != right.reverse_channels:
        mismatch(
            "reverse linkage", sorted(mapped_reverse), sorted(right.reverse_channels)
        )
    left_missing_targets = {
        target for target in left.reverse_channels if target not in channel_mapping
    }
    mapped_right_channels = set(channel_mapping.values())
    right_added_targets = {
        target
        for target in right.reverse_channels
        if target not in mapped_right_channels
    }
    left_link_errors = {
        "unresolved": max(
            0, len(left.unresolved_reverse_ids) - len(right_added_targets)
        ),
        "nonreciprocal": left.nonreciprocal_reverse,
    }
    right_link_errors = {
        "unresolved": max(
            0, len(right.unresolved_reverse_ids) - len(left_missing_targets)
        ),
        "nonreciprocal": right.nonreciprocal_reverse,
    }
    if left_link_errors != right_link_errors:
        mismatch("reverse linkage validity", left_link_errors, right_link_errors)
    return reasons


def compare_artifacts(
    left: str | Path | dict[str, Any],
    right: str | Path | dict[str, Any],
    *,
    rtol: float = 1e-9,
) -> dict[str, Any]:
    """Compare two artifacts and return a JSON-serializable report."""
    if not math.isfinite(rtol) or rtol < 0:
        raise ValueError("rtol must be a finite non-negative number")
    started = time.perf_counter()
    same_input = left is right
    if not same_input and not isinstance(left, dict) and not isinstance(right, dict):
        same_input = Path(left).resolve() == Path(right).resolve()
    graph_cache: dict[tuple[Any, ...], tuple[_CompactGraph, _CompactGraph]] = {}
    left_summary = _summarize(_load_artifact(left), graph_cache)
    right_summary = (
        left_summary if same_input else _summarize(_load_artifact(right), graph_cache)
    )

    mapping: dict[int, int] = {}
    used_right: set[int] = set()
    missing: list[dict[str, Any]] = []
    for left_channel in left_summary.channels:
        match = None
        for candidate_index in right_summary.by_fingerprint.get(
            left_channel.chemistry.fingerprint, []
        ):
            if candidate_index in used_right:
                continue
            candidate = right_summary.channels[candidate_index]
            if same_input or _isomorphic(left_channel.chemistry, candidate.chemistry):
                match = candidate_index
                break
        if match is None:
            missing.append(
                {
                    **_channel_ref(left_channel),
                    "reason": "rooted atom-mapped rewrite absent from right artifact",
                }
            )
        else:
            mapping[left_channel.index] = match
            used_right.add(match)

    added = [
        {
            **_channel_ref(channel),
            "reason": "rooted atom-mapped rewrite absent from left artifact",
        }
        for channel in right_summary.channels
        if channel.index not in used_right
    ]
    changed = []
    for left_index, right_index in mapping.items():
        left_channel = left_summary.channels[left_index]
        right_channel = right_summary.channels[right_index]
        reasons = _compare_channel_data(
            left_channel,
            right_channel,
            left_summary,
            right_summary,
            mapping,
            rtol,
        )
        if reasons:
            changed.append({**_channel_ref(left_channel), "reasons": reasons})

    report = {
        "equal": not (missing or added or changed),
        "relative_tolerance": rtol,
        "left": {
            "records": left_summary.record_count,
            "channels": len(left_summary.channels),
        },
        "right": {
            "records": right_summary.record_count,
            "channels": len(right_summary.channels),
        },
        "missing": missing,
        "added": added,
        "changed": changed,
        "elapsed_seconds": time.perf_counter() - started,
    }
    return report


def human_summary(report: dict[str, Any]) -> str:
    """Render a concise human-readable result."""
    left = report["left"]
    right = report["right"]
    if report["equal"]:
        return (
            f"EQUAL: {left['channels']} physical channels "
            f"({left['records']} vs {right['records']} records)"
        )
    lines = [
        "NOT EQUAL: "
        f"{len(report['missing'])} missing, {len(report['added'])} added, "
        f"{len(report['changed'])} changed channels"
    ]
    for category in ("missing", "added"):
        for item in report[category]:
            lines.append(f"{category.upper()} {item['channel']}: {item['reason']}")
    for item in report["changed"]:
        reasons = ", ".join(reason["reason"] for reason in item["reasons"])
        lines.append(f"CHANGED {item['channel']}: {reasons}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("left", help="first compiled event-set JSON artifact")
    parser.add_argument("right", help="second compiled event-set JSON artifact")
    parser.add_argument(
        "--rtol",
        type=float,
        default=1e-9,
        help="relative numeric tolerance (default: 1e-9)",
    )
    parser.add_argument(
        "--json",
        metavar="PATH",
        default="-",
        help="write the JSON report to PATH (default: stdout)",
    )
    arguments = parser.parse_args(argv)
    try:
        report = compare_artifacts(arguments.left, arguments.right, rtol=arguments.rtol)
    except (OSError, json.JSONDecodeError, EquivalenceError, ValueError) as error:
        parser.error(str(error))
    summary = human_summary(report)
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.json == "-":
        sys.stdout.write(payload)
        print(summary, file=sys.stderr)
    else:
        Path(arguments.json).write_text(payload, encoding="utf-8")
        print(summary)
    return 0 if report["equal"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
