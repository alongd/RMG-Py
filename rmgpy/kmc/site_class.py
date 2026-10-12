"""Canonical rooted local graph identities for kMC reaction sites."""

from __future__ import annotations

from collections import deque
from hashlib import sha256
import json
from typing import Any, Iterable, Mapping, Sequence


SITE_CLASS_SCHEMA = "SiteClassKeyV1"


class SiteClassError(ValueError):
    """A fail-closed site-class normalization error with a stable code."""

    def __init__(self, code: str, message: str):
        super().__init__(f"{code}: {message}")
        self.code = code


_REQUIRED_POLICY_FIELDS = {
    "allowed_atom_stereo",
    "allowed_bond_orders",
    "allowed_bond_stereo",
    "allowed_charges",
    "allowed_elements",
    "allowed_isotopes",
    "allowed_lone_pairs",
    "allowed_multiplicities",
    "allowed_radicals",
    "canonicalizer_version",
    "max_degree",
    "max_explicit_h",
    "max_heavy_atoms",
    "max_radius",
    "max_root_size",
    "max_valence_by_element",
    "orientation",
    "policy_id",
    "radius",
}

_IGNORED_ATOM_PROPS = {"display_label", "inRing", "persistent_uuid", "uuid"}
_CONTROL_ATOM_PROPS = {"artificial_cap", "boundary_port", "reacting", "stereo"}
_ALLOWED_GRAPH_PROPS = {
    "atom_stereo",
    "bond_stereo",
    "display_label",
    "persistent_uuid",
    "uuid",
}


def _canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _fail(code: str, message: str) -> None:
    raise SiteClassError(code, message)


def _normalize_policy(locality_policy: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(locality_policy, Mapping):
        _fail("invalid_policy", "locality_policy must be a mapping")
    missing = sorted(_REQUIRED_POLICY_FIELDS - set(locality_policy))
    if missing:
        _fail("invalid_policy", f"missing fields: {', '.join(missing)}")
    unknown = sorted(set(locality_policy) - _REQUIRED_POLICY_FIELDS)
    if unknown:
        _fail("invalid_policy", f"unsupported fields: {', '.join(unknown)}")

    policy = dict(locality_policy)
    for name in (
        "allowed_atom_stereo",
        "allowed_bond_orders",
        "allowed_bond_stereo",
        "allowed_charges",
        "allowed_elements",
        "allowed_isotopes",
        "allowed_lone_pairs",
        "allowed_multiplicities",
        "allowed_radicals",
    ):
        value = policy[name]
        if not isinstance(value, (list, tuple)) or not value:
            _fail("invalid_policy", f"{name} must be a nonempty finite sequence")
        policy[name] = sorted(set(value))

    valence_table = policy["max_valence_by_element"]
    if not isinstance(valence_table, Mapping) or not valence_table:
        _fail("invalid_policy", "max_valence_by_element must be a nonempty mapping")
    if not set(policy["allowed_elements"]).issubset(valence_table):
        _fail(
            "invalid_policy",
            "max_valence_by_element must define every allowed element",
        )
    normalized_valences = {}
    for element, maximum in valence_table.items():
        if (
            not isinstance(maximum, (int, float))
            or isinstance(maximum, bool)
            or maximum <= 0
        ):
            _fail("invalid_policy", f"invalid maximum valence for {element!r}")
        normalized_valences[element] = float(maximum)
    policy["max_valence_by_element"] = normalized_valences

    for name in (
        "max_degree",
        "max_explicit_h",
        "max_heavy_atoms",
        "max_radius",
        "max_root_size",
        "radius",
    ):
        value = policy[name]
        if not isinstance(value, int) or isinstance(value, bool) or value < 0:
            _fail("invalid_policy", f"{name} must be a nonnegative integer")
    if policy["radius"] > policy["max_radius"]:
        _fail("domain_exceeded", "radius exceeds max_radius")
    if not isinstance(policy["policy_id"], str) or not policy["policy_id"]:
        _fail("invalid_policy", "policy_id must be a nonempty string")
    if (
        not isinstance(policy["canonicalizer_version"], str)
        or not policy["canonicalizer_version"]
    ):
        _fail("invalid_policy", "canonicalizer_version must be a nonempty string")
    if not isinstance(policy["orientation"], str) or not policy["orientation"]:
        _fail("invalid_policy", "orientation must be a nonempty string")
    try:
        _canonical_json(policy)
    except (TypeError, ValueError) as exc:
        _fail("invalid_policy", f"policy is not canonical-JSON serializable: {exc}")
    return policy


def _atoms(graph: Any) -> list[Any]:
    atoms = getattr(graph, "atoms", None)
    if atoms is None:
        atoms = getattr(graph, "vertices", None)
    if atoms is None:
        _fail("unsupported_field", "graph exposes neither atoms nor vertices")
    return list(atoms)


def _is_hydrogen(atom: Any) -> bool:
    return getattr(atom, "symbol", None) == "H"


def _props(atom: Any) -> dict[str, Any]:
    value = getattr(atom, "props", None)
    return dict(value or {})


def _is_cap(atom: Any) -> bool:
    return _props(atom).get("artificial_cap") is True


def _bond_between(atom1: Any, atom2: Any) -> Any:
    try:
        return atom1.bonds[atom2]
    except (AttributeError, KeyError):
        _fail("malformed_root", "root bond endpoints are not bonded in graph")


def _normalize_roots(
    graph_atoms: Sequence[Any], ordered_root: Iterable[Any], max_root_size: int
) -> tuple[list[dict[str, Any]], list[Any]]:
    if ordered_root is None:
        _fail("missing_root", "ordered_root is required")
    try:
        items = list(ordered_root)
    except TypeError:
        _fail("malformed_root", "ordered_root must be iterable")
    if not items:
        _fail("missing_root", "ordered_root is empty")
    if len(items) > max_root_size:
        _fail("domain_exceeded", "root tuple exceeds max_root_size")

    graph_atom_ids = {id(atom) for atom in graph_atoms}
    normalized: list[dict[str, Any]] = []
    root_atoms: list[Any] = []
    roles: set[str] = set()
    for index, item in enumerate(items):
        if not isinstance(item, (tuple, list)) or len(item) != 2:
            _fail("malformed_root", f"root item {index} must be (role, atom_or_bond)")
        role, entity = item
        if not isinstance(role, str) or not role:
            _fail("malformed_root", f"root item {index} has no role label")
        if role in roles:
            _fail("malformed_root", f"duplicate root role {role!r}")
        roles.add(role)

        if id(entity) in graph_atom_ids:
            endpoints = (entity,)
            kind = "atom"
        elif hasattr(entity, "atom1") and hasattr(entity, "atom2"):
            endpoints = (entity.atom1, entity.atom2)
            if any(id(atom) not in graph_atom_ids for atom in endpoints):
                _fail("malformed_root", f"root bond {role!r} is not in graph")
            if _bond_between(*endpoints) is not entity:
                _fail("malformed_root", f"root bond {role!r} is not graph-owned")
            kind = "bond"
        else:
            _fail(
                "malformed_root", f"root item {role!r} is neither a graph atom nor bond"
            )

        for atom in endpoints:
            if _is_hydrogen(atom):
                _fail("invalid_root", f"root {role!r} contains hydrogen")
            props = _props(atom)
            if props.get("artificial_cap") is True:
                _fail("invalid_root", f"root {role!r} contains an artificial cap")
            if props.get("reacting") is True:
                _fail(
                    "invalid_root", f"root {role!r} contains a pre-marked reacting atom"
                )
            if all(existing is not atom for existing in root_atoms):
                root_atoms.append(atom)
        normalized.append({"role": role, "kind": kind, "endpoints": endpoints})
    return normalized, root_atoms


def _allowed(value: Any, policy: Mapping[str, Any], field: str, name: str) -> None:
    if value not in policy[field]:
        _fail("domain_exceeded", f"{name}={value!r} is outside {field}")


def _atom_stereo(graph: Any, atom: Any) -> Any:
    graph_props = getattr(graph, "props", None) or {}
    annotations = graph_props.get("atom_stereo", {})
    if annotations and not isinstance(annotations, Mapping):
        _fail("unsupported_field", "graph.props['atom_stereo'] must be a mapping")
    explicit = annotations.get(atom)
    prop_value = _props(atom).get("stereo")
    if explicit is not None and prop_value is not None and explicit != prop_value:
        _fail("unsupported_field", "conflicting atom stereo annotations")
    return explicit if explicit is not None else prop_value


def _validate_atom(atom: Any, graph: Any, policy: Mapping[str, Any]) -> None:
    props = _props(atom)
    unknown_props = set(props) - _IGNORED_ATOM_PROPS - _CONTROL_ATOM_PROPS
    if unknown_props:
        _fail("unsupported_field", f"unsupported atom props: {sorted(unknown_props)}")
    if getattr(atom, "site", "") or getattr(atom, "morphology", ""):
        _fail("unsupported_field", "surface site and morphology are unsupported")

    element = getattr(atom, "element", None)
    if element is None:
        _fail("unsupported_field", "atom has no element")
    _allowed(atom.symbol, policy, "allowed_elements", "element")
    _allowed(element.isotope, policy, "allowed_isotopes", "isotope")
    _allowed(atom.charge, policy, "allowed_charges", "charge")
    _allowed(atom.radical_electrons, policy, "allowed_radicals", "radical_electrons")
    _allowed(atom.lone_pairs, policy, "allowed_lone_pairs", "lone_pairs")
    valence = 0.0
    for bond in atom.bonds.values():
        try:
            order = float(bond.order)
        except (TypeError, ValueError):
            _fail("unsupported_field", f"nonnumeric bond order {bond.order!r}")
        _allowed(order, policy, "allowed_bond_orders", "bond_order")
        valence += order
    if valence > policy["max_valence_by_element"][atom.symbol]:
        _fail("domain_exceeded", f"valence {valence!r} exceeds finite policy bound")
    stereo = _atom_stereo(graph, atom)
    if stereo is not None:
        if stereo not in policy["allowed_atom_stereo"]:
            _fail("unsupported_stereo", f"unsupported atom stereo {stereo!r}")


def _bond_stereo(graph: Any, bond: Any) -> Any:
    graph_props = getattr(graph, "props", None) or {}
    annotations = graph_props.get("bond_stereo", {})
    if annotations and not isinstance(annotations, Mapping):
        _fail("unsupported_field", "graph.props['bond_stereo'] must be a mapping")
    return annotations.get(bond)


def _bond_descriptor(
    graph: Any, bond: Any, policy: Mapping[str, Any]
) -> dict[str, Any]:
    order = float(bond.order)
    _allowed(order, policy, "allowed_bond_orders", "bond_order")
    stereo = _bond_stereo(graph, bond)
    if stereo is not None and stereo not in policy["allowed_bond_stereo"]:
        _fail("unsupported_stereo", f"unsupported bond stereo {stereo!r}")
    return {"order": order, "stereo": stereo}


def _hydrogen_descriptor(
    graph: Any, atom: Any, bond: Any, policy: Mapping[str, Any]
) -> dict[str, Any]:
    _validate_atom(atom, graph, policy)
    return {
        "bond": _bond_descriptor(graph, bond, policy),
        **_chemical_state_descriptor(graph, atom),
    }


def _chemical_state_descriptor(graph: Any, atom: Any) -> dict[str, Any]:
    return {
        "charge": atom.charge,
        "isotope": atom.element.isotope,
        "lone_pairs": atom.lone_pairs,
        "radical_electrons": atom.radical_electrons,
        "stereo": _atom_stereo(graph, atom),
    }


def _atom_descriptor(
    graph: Any, atom: Any, policy: Mapping[str, Any]
) -> dict[str, Any]:
    _validate_atom(atom, graph, policy)
    hydrogens = []
    for neighbor, bond in atom.bonds.items():
        if _is_hydrogen(neighbor) and not _is_cap(neighbor):
            hydrogens.append(_hydrogen_descriptor(graph, neighbor, bond, policy))
    if len(hydrogens) > policy["max_explicit_h"]:
        _fail("domain_exceeded", "explicit hydrogen count exceeds max_explicit_h")
    hydrogens.sort(key=_canonical_json)
    return {
        **_chemical_state_descriptor(graph, atom),
        "element": atom.symbol,
        "explicit_hydrogens": hydrogens,
    }


def _heavy_neighbors(atom: Any) -> list[Any]:
    return [
        neighbor
        for neighbor in atom.bonds
        if not _is_hydrogen(neighbor) and not _is_cap(neighbor)
    ]


def _context_atoms(root_atoms: Sequence[Any], radius: int) -> list[Any]:
    distances = {id(atom): 0 for atom in root_atoms}
    by_id = {id(atom): atom for atom in root_atoms}
    queue = deque(root_atoms)
    while queue:
        atom = queue.popleft()
        distance = distances[id(atom)]
        if distance == radius:
            continue
        for neighbor in _heavy_neighbors(atom):
            neighbor_id = id(neighbor)
            if neighbor_id not in distances:
                distances[neighbor_id] = distance + 1
                by_id[neighbor_id] = neighbor
                queue.append(neighbor)
    return list(by_id.values())


def _refine(
    seed: Sequence[Any], adjacency: Sequence[list[tuple[str, int]]]
) -> list[int]:
    tokens = [_canonical_json(value) for value in seed]
    while True:
        palette = {value: index for index, value in enumerate(sorted(set(tokens)))}
        colors = [palette[value] for value in tokens]
        signatures = []
        for index, color in enumerate(colors):
            neighbors = sorted(
                (edge, colors[other]) for edge, other in adjacency[index]
            )
            signatures.append(_canonical_json([color, neighbors]))
        next_palette = {
            value: index for index, value in enumerate(sorted(set(signatures)))
        }
        next_colors = [next_palette[value] for value in signatures]
        if next_colors == colors:
            return colors
        tokens = signatures


def _canonical_order(
    node_seeds: Sequence[Any],
    adjacency: Sequence[list[tuple[str, int]]],
    serialize,
) -> tuple[list[int], dict[str, Any]]:
    best_text: str | None = None
    best_order: list[int] | None = None
    best_pattern: dict[str, Any] | None = None

    def search(seed: Sequence[Any]) -> None:
        nonlocal best_text, best_order, best_pattern
        colors = _refine(seed, adjacency)
        cells: dict[int, list[int]] = {}
        for index, color in enumerate(colors):
            cells.setdefault(color, []).append(index)
        unresolved = [cell for cell in cells.values() if len(cell) > 1]
        if not unresolved:
            order = sorted(range(len(colors)), key=colors.__getitem__)
            pattern = serialize(order)
            text = _canonical_json(pattern)
            if (
                best_text is None
                or text < best_text
                or (text == best_text and tuple(order) < tuple(best_order or ()))
            ):
                best_text = text
                best_order = order
                best_pattern = pattern
            return
        cell = min(unresolved, key=lambda value: (len(value), colors[value[0]]))
        for chosen in cell:
            individualized = [
                [colors[index], 0 if index == chosen else 1]
                for index in range(len(colors))
            ]
            search(individualized)

    search(node_seeds)
    assert best_order is not None and best_pattern is not None
    return best_order, best_pattern


def extract_site_class(
    graph: Any,
    ordered_root: Iterable[Any],
    locality_policy: Mapping[str, Any],
) -> tuple[str, dict[int, Any], dict[str, Any]]:
    """Return a canonical finite site-class key, input mapping, and port pattern."""

    policy = _normalize_policy(locality_policy)
    graph_atoms = _atoms(graph)
    graph_props = getattr(graph, "props", None) or {}
    if not isinstance(graph_props, Mapping):
        _fail("unsupported_field", "graph props must be a mapping")
    unknown_graph_props = set(graph_props) - _ALLOWED_GRAPH_PROPS
    if unknown_graph_props:
        _fail(
            "unsupported_field",
            f"unsupported graph props: {sorted(unknown_graph_props)}",
        )
    roots, root_atoms = _normalize_roots(
        graph_atoms, ordered_root, policy["max_root_size"]
    )
    for atom in graph_atoms:
        props = _props(atom)
        is_cap = props.get("artificial_cap") is True
        has_port = "boundary_port" in props
        if has_port and not is_cap:
            _fail("malformed_port", "boundary_port annotation requires artificial_cap")
        if is_cap:
            if not _is_hydrogen(atom):
                _fail("malformed_port", "an artificial cap must be hydrogen")
            _validate_atom(atom, graph, policy)
    multiplicity = getattr(graph, "multiplicity", None)
    if multiplicity in (None, -187):
        _fail("missing_multiplicity", "graph multiplicity is not authoritative")
    _allowed(multiplicity, policy, "allowed_multiplicities", "multiplicity")
    if getattr(graph, "metal", "") or getattr(graph, "facet", ""):
        _fail("unsupported_field", "surface metal and facet are unsupported")

    context = _context_atoms(root_atoms, policy["radius"])
    if len(context) > policy["max_heavy_atoms"]:
        _fail("domain_exceeded", "context exceeds max_heavy_atoms")
    context_ids = {id(atom) for atom in context}
    local_index = {id(atom): index for index, atom in enumerate(context)}
    atom_descriptors = [_atom_descriptor(graph, atom, policy) for atom in context]

    for atom in context:
        degree = len(_heavy_neighbors(atom))
        if degree > policy["max_degree"]:
            _fail("domain_exceeded", "heavy-atom degree exceeds max_degree")

    root_entries = []
    atom_root_marks: dict[int, list[list[Any]]] = {}
    root_bonds: dict[tuple[int, int], list[list[Any]]] = {}
    for ordinal, root in enumerate(roots):
        endpoints = [local_index[id(atom)] for atom in root["endpoints"]]
        if root["kind"] == "atom":
            atom_root_marks.setdefault(endpoints[0], []).append([ordinal, root["role"]])
        else:
            edge = tuple(sorted(endpoints))
            root_bonds.setdefault(edge, []).append([ordinal, root["role"]])
        root_entries.append(
            {"kind": root["kind"], "role": root["role"], "endpoints": endpoints}
        )

    edges: dict[tuple[int, int], dict[str, Any]] = {}
    adjacency: list[list[tuple[str, int]]] = [[] for _ in context]
    for left, atom in enumerate(context):
        for neighbor, bond in atom.bonds.items():
            if id(neighbor) not in context_ids:
                continue
            right = local_index[id(neighbor)]
            if left >= right:
                continue
            descriptor = _bond_descriptor(graph, bond, policy)
            marks = root_bonds.get((left, right), [])
            edge_value = {"bond": descriptor, "root_roles": marks}
            edges[(left, right)] = edge_value
            edge_token = _canonical_json(edge_value)
            adjacency[left].append((edge_token, right))
            adjacency[right].append((edge_token, left))

    ports: list[dict[str, Any]] = []
    for source, atom in enumerate(context):
        for neighbor, bond in atom.bonds.items():
            if id(neighbor) in context_ids:
                continue
            if _is_cap(neighbor):
                annotation = _props(neighbor).get("boundary_port")
                if annotation is None:
                    _fail(
                        "malformed_port",
                        "artificial cap lacks boundary_port annotation",
                    )
                signature = annotation
            elif _is_hydrogen(neighbor):
                continue
            else:
                _validate_atom(neighbor, graph, policy)
                signature = {
                    "bond_orders": sorted(
                        float(other_bond.order)
                        for other, other_bond in neighbor.bonds.items()
                        if not _is_hydrogen(other) and not _is_cap(other)
                    ),
                    "element": neighbor.symbol,
                }
            if not isinstance(signature, Mapping):
                _fail("malformed_port", "boundary_port must be a mapping")
            if set(signature) != {"bond_orders", "element"}:
                _fail(
                    "malformed_port", "boundary_port requires element and bond_orders"
                )
            if signature["element"] not in policy["allowed_elements"]:
                _fail("domain_exceeded", "boundary element is outside allowed_elements")
            if not isinstance(signature["bond_orders"], (list, tuple)):
                _fail("malformed_port", "boundary bond_orders must be a sequence")
            try:
                orders = sorted(float(order) for order in signature["bond_orders"])
            except (TypeError, ValueError):
                _fail("malformed_port", "boundary bond_orders must be numeric")
            for order in orders:
                _allowed(order, policy, "allowed_bond_orders", "boundary_bond_order")
            ports.append(
                {
                    "bond": _bond_descriptor(graph, bond, policy),
                    "excluded": {
                        "bond_orders": orders,
                        "element": signature["element"],
                    },
                    "source": source,
                }
            )

    ports_by_source: dict[int, list[dict[str, Any]]] = {}
    for port in ports:
        ports_by_source.setdefault(port["source"], []).append(
            {"bond": port["bond"], "excluded": port["excluded"]}
        )
    node_seeds = []
    for index, descriptor in enumerate(atom_descriptors):
        node_seeds.append(
            {
                "atom": descriptor,
                "ports": sorted(ports_by_source.get(index, []), key=_canonical_json),
                "root_roles": atom_root_marks.get(index, []),
            }
        )

    policy_text = _canonical_json(policy)
    policy_digest = sha256(policy_text.encode("ascii")).hexdigest()

    def serialize(order: Sequence[int]) -> dict[str, Any]:
        position = {local: canonical for canonical, local in enumerate(order)}
        canonical_bonds = []
        for (left, right), value in edges.items():
            endpoints = sorted([position[left], position[right]])
            canonical_bonds.append(
                {
                    "atoms": endpoints,
                    "bond": value["bond"],
                    "root_roles": value["root_roles"],
                }
            )
        canonical_bonds.sort(key=_canonical_json)
        canonical_ports = []
        for port in ports:
            canonical_ports.append(
                {
                    "bond": port["bond"],
                    "excluded": port["excluded"],
                    "source": position[port["source"]],
                }
            )
        canonical_ports.sort(key=_canonical_json)
        canonical_roots = []
        for root in root_entries:
            endpoints = [position[endpoint] for endpoint in root["endpoints"]]
            entry = {"kind": root["kind"], "role": root["role"]}
            if root["kind"] == "atom":
                entry["position"] = endpoints[0]
            else:
                entry["positions"] = sorted(endpoints)
            canonical_roots.append(entry)
        return {
            "atoms": [atom_descriptors[local] for local in order],
            "bonds": canonical_bonds,
            "boundary_ports": canonical_ports,
            "multiplicity": multiplicity,
            "ordered_root": canonical_roots,
            "orientation": policy["orientation"],
            "radius": policy["radius"],
            "schema": "RootedPortPatternV1",
        }

    order, pattern = _canonical_order(node_seeds, adjacency, serialize)
    key_document = {
        "canonicalizer_version": policy["canonicalizer_version"],
        "pattern": pattern,
        "policy_digest": policy_digest,
        "policy_id": policy["policy_id"],
        "schema": SITE_CLASS_SCHEMA,
    }
    key = sha256(_canonical_json(key_document).encode("ascii")).hexdigest()
    mapping = {position: context[local] for position, local in enumerate(order)}
    return key, mapping, pattern
