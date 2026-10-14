"""Atom-native melt state for the kMC runtime.

Atoms are the only authority.  A molecule is a connected component of the
atom graph, so chains of any length, cuts anywhere, transfers and free small
fragments (H, H2, ...) need no side ledgers.  Every hydrogen is an explicit
atom node.

The module covers: the state itself, deterministic materialization from a
``LinearHomopolymerSpec`` declaration, the component index, journaled atomic
edits with exact rollback, a canonical state hash, and a few observables.
It does not match templates, resolve resonance, identify species, release or
sink fragments, or touch the SSA; those are later layers.

Edits never create or delete atoms.  An edit set is applied in a transaction:
preconditions (each operation states the expected old value) are checked, the
mutation support (every changed node, both endpoints of every changed edge) is
journaled, the edits are written, the components of the union of touched
components are recomputed, and postconditions (valence per element, symmetric
adjacency, unchanged atom-id set, canonical component ids) are verified.  Any
failure replays the journal and leaves the state byte-identical.
"""

import hashlib
from collections import Counter
from dataclasses import dataclass, replace
from fractions import Fraction
from types import MappingProxyType
from typing import Callable, Dict, Iterable, Mapping, Optional, Tuple

from rmgpy.molecule.element import PeriodicSystem, get_element

# Stages at which a fault can be injected (test seam of ``MeltState.apply``).
STAGES = ("after_apply", "after_components", "after_postcheck")

# Elements whose valence is a closed shell: bond orders + radical electrons +
# 2 * lone pairs == valence electrons - charge.  Hypervalent elements are not
# supported and are refused with a clear error.
SUPPORTED_ELEMENTS = ("H", "C", "N", "O", "F", "Cl", "Br", "I")
BOND_ORDERS = (1.0, 1.5, 2.0, 3.0)


class MeltStateError(ValueError):
    """Base class for every refusal of this module."""


class ValenceError(MeltStateError):
    """A node violates the valence rule of its element."""


class EditError(MeltStateError):
    """An edit set is malformed or its stated old values do not hold."""


class PostconditionError(MeltStateError):
    """An applied edit set broke an invariant (and was rolled back)."""


# --------------------------------------------------------------------------
# nodes and origin
# --------------------------------------------------------------------------
class _Edges(Mapping):
    """Immutable, picklable id -> bond-order mapping (sorted by id)."""

    __slots__ = ("_items", "_dict")

    def __init__(self, items):
        pairs = tuple(sorted((int(j), float(o)) for j, o in dict(items).items()))
        object.__setattr__(self, "_items", pairs)
        object.__setattr__(self, "_dict", dict(pairs))

    def __setattr__(self, name, value):
        raise AttributeError("edges are immutable")

    def __getitem__(self, key):
        return self._dict[key]

    def __iter__(self):
        return iter(self._dict)

    def __len__(self):
        return len(self._dict)

    def __hash__(self):
        return hash(self._items)

    def __eq__(self, other):
        if isinstance(other, _Edges):
            return self._items == other._items
        return Mapping.__eq__(self, other)

    def __reduce__(self):
        return (_Edges, (self._items,))

    def __repr__(self):
        return f"_Edges({dict(self._items)!r})"


@dataclass(frozen=True)
class Node:
    """One atom: element, isotope (mass number or None), charge, radical
    electrons, lone pairs and its bonds (neighbour id -> bond order)."""

    element: str
    isotope: Optional[int]
    charge: int
    radical: int
    lone_pairs: int
    edges: Mapping

    def __post_init__(self):
        if not isinstance(self.element, str) or not self.element:
            raise MeltStateError("a node needs an element symbol")
        if self.isotope is not None and (
                isinstance(self.isotope, bool) or not isinstance(self.isotope, int)
                or self.isotope < 1):
            raise MeltStateError(f"isotope must be a mass number or None, got {self.isotope!r}")
        for name in ("charge", "radical", "lone_pairs"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int):
                raise MeltStateError(f"{name} must be an integer, got {value!r}")
        if self.radical < 0 or self.lone_pairs < 0:
            raise MeltStateError("radical electrons and lone pairs cannot be negative")
        edges = self.edges if isinstance(self.edges, _Edges) else _Edges(self.edges)
        for neighbour, order in edges.items():
            if order not in BOND_ORDERS:
                raise MeltStateError(f"bond order {order!r} is not one of {BOND_ORDERS}")
        object.__setattr__(self, "edges", edges)


@dataclass(frozen=True)
class Origin:
    """Immutable provenance of an atom.

    ``kind`` is ``chain`` or ``mobile``; ``ordinal`` counts molecules of that
    kind in materialization order; ``segment`` is ``head_end`` / ``unit k`` /
    ``tail_end`` for chains and the species key for mobile molecules; ``key``
    is the declaration atom key (hydrogens: ``<parent key>.H<i>``);
    ``is_backbone`` flags atoms on the repeat-unit path between its head and
    tail ports.
    """

    kind: str
    ordinal: int
    segment: str
    key: str
    is_backbone: bool

    def __post_init__(self):
        if self.kind not in ("chain", "mobile"):
            raise MeltStateError(f"origin kind must be chain or mobile, got {self.kind!r}")


def _check_valence(atom_id: int, node: Node):
    if node.element not in SUPPORTED_ELEMENTS:
        raise ValenceError(
            f"atom {atom_id}: element {node.element!r} is not supported "
            f"(supported: {', '.join(SUPPORTED_ELEMENTS)})")
    expected = PeriodicSystem.valence_electrons[node.element] - node.charge
    used = sum(node.edges.values()) + node.radical + 2 * node.lone_pairs
    if used != expected:
        raise ValenceError(
            f"atom {atom_id} ({node.element}): bond orders {sum(node.edges.values())} + "
            f"radicals {node.radical} + 2*lone pairs {node.lone_pairs} = {used}, "
            f"valence data require {expected} at charge {node.charge}")


# --------------------------------------------------------------------------
# components
# --------------------------------------------------------------------------
_MASS_CACHE: Dict[Tuple[str, Optional[int]], Fraction] = {}


def _element_mass(element: str, isotope: Optional[int]) -> Fraction:
    key = (element, isotope)
    if key not in _MASS_CACHE:
        _MASS_CACHE[key] = Fraction(get_element(element, -1 if isotope is None else isotope).mass)
    return _MASS_CACHE[key]


@dataclass(frozen=True)
class _Record:
    atoms: Tuple[int, ...]
    composition: Tuple[Tuple[Tuple[str, Optional[int]], int], ...]
    radicals: int
    mass: Fraction
    backbone: int


def _make_record(members: Tuple[int, ...], atoms: Mapping, origin: Mapping) -> _Record:
    counts = Counter()
    radicals = backbone = 0
    for a in members:
        node = atoms[a]
        counts[(node.element, node.isotope)] += 1
        radicals += node.radical
        backbone += 1 if origin[a].is_backbone else 0
    composition = tuple(sorted(counts.items(), key=lambda kv: (kv[0][0], kv[0][1] or 0)))
    mass = sum((_element_mass(el, iso) * n for (el, iso), n in composition), Fraction(0))
    return _Record(members, composition, radicals, mass, backbone)


class ComponentIndex:
    """atom -> component id (= minimum atom id of the component), with sorted
    atom tuples, composition, radical count, exact mass and backbone count."""

    def __init__(self):
        self._of: Dict[int, int] = {}
        self._rec: Dict[int, _Record] = {}

    # queries
    def ids(self) -> Tuple[int, ...]:
        return tuple(sorted(self._rec))

    def __len__(self):
        return len(self._rec)

    def component_of(self, atom: int) -> int:
        return self._of[atom]

    def atoms(self, cid: int) -> Tuple[int, ...]:
        return self._rec[cid].atoms

    def composition(self, cid: int) -> Dict[Tuple[str, Optional[int]], int]:
        return dict(self._rec[cid].composition)

    def radical_count(self, cid: int) -> int:
        return self._rec[cid].radicals

    def mass(self, cid: int) -> Fraction:
        """Exact molecular mass (kg/mol) from the isotope-aware composition."""
        return self._rec[cid].mass

    def backbone_count(self, cid: int) -> int:
        return self._rec[cid].backbone

    # mutation (internal)
    def _install(self, remove: Iterable[int], records: Iterable[_Record]):
        for cid in remove:
            del self._rec[cid]
        for record in records:
            cid = record.atoms[0]
            self._rec[cid] = record
            for a in record.atoms:
                self._of[a] = cid

    def __getstate__(self):
        return (tuple(sorted(self._of.items())),
                tuple(sorted(self._rec.items())))

    def __setstate__(self, state):
        of, rec = state
        self._of = dict(of)
        self._rec = dict(rec)


def _partition(universe: Iterable[int], atoms: Mapping) -> Tuple[Tuple[int, ...], ...]:
    """Connected components of ``universe`` (sorted atom tuples, sorted by id)."""
    universe = sorted(universe)
    allowed = set(universe)
    seen = set()
    parts = []
    for start in universe:
        if start in seen:
            continue
        seen.add(start)
        stack, members = [start], [start]
        while stack:
            for j in atoms[stack.pop()].edges:
                if j not in allowed:
                    raise PostconditionError(f"bond leaves the recomputed region at atom {j}")
                if j not in seen:
                    seen.add(j)
                    stack.append(j)
                    members.append(j)
        parts.append(tuple(sorted(members)))
    return tuple(parts)


# --------------------------------------------------------------------------
# edits
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class SetBond:
    """Change the bond order between two existing atoms (0 = none)."""

    a: int
    b: int
    old: float
    new: float


@dataclass(frozen=True)
class SetRadical:
    atom: int
    old: int
    new: int


@dataclass(frozen=True)
class SetCharge:
    atom: int
    old: int
    new: int


@dataclass(frozen=True)
class SetLonePairs:
    atom: int
    old: int
    new: int


_FIELD_OPS = {SetRadical: "radical", SetCharge: "charge", SetLonePairs: "lone_pairs"}


class Journal:
    """Record of one transaction: the full mutation support and old membership."""

    def __init__(self, edits):
        self.edits = edits
        self.old_nodes: Dict[int, Node] = {}
        self.old_components: Tuple[Tuple[int, Tuple[int, ...]], ...] = ()
        self.new_components: Tuple[Tuple[int, Tuple[int, ...]], ...] = ()
        self._old_records: Dict[int, _Record] = {}
        self._installed = False


def _is_int(value) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


# --------------------------------------------------------------------------
# state
# --------------------------------------------------------------------------
class MeltState:
    """Atoms (authoritative), their immutable origin, and the component index."""

    def __init__(self, atoms: Mapping[int, Node], origin: Mapping[int, Origin]):
        if not atoms:
            raise MeltStateError("a melt state cannot be empty")
        if set(atoms) != set(origin):
            raise MeltStateError("atoms and origin must cover the same atom ids")
        for i, node in atoms.items():
            if not _is_int(i) or i < 0:
                raise MeltStateError(f"atom ids must be non-negative integers, got {i!r}")
            if not isinstance(node, Node):
                raise MeltStateError(f"atom {i} is not a Node")
            if not isinstance(origin[i], Origin):
                raise MeltStateError(f"atom {i} has no Origin")
        self._atoms: Dict[int, Node] = {i: atoms[i] for i in sorted(atoms)}
        self._origin: Dict[int, Origin] = {i: origin[i] for i in sorted(origin)}
        self._ids = frozenset(self._atoms)
        for i, node in self._atoms.items():
            self._check_node(i, node)
        self.components = ComponentIndex()
        parts = _partition(self._atoms, self._atoms)
        self.components._install((), (_make_record(p, self._atoms, self._origin) for p in parts))

    # -- read access -------------------------------------------------------
    @property
    def atoms(self) -> Mapping[int, Node]:
        return MappingProxyType(self._atoms)

    @property
    def origin(self) -> Mapping[int, Origin]:
        return MappingProxyType(self._origin)

    # -- validation --------------------------------------------------------
    def _check_node(self, i: int, node: Node):
        _check_valence(i, node)
        for j, order in node.edges.items():
            if j == i:
                raise MeltStateError(f"atom {i} is bonded to itself")
            other = self._atoms.get(j)
            if other is None:
                raise MeltStateError(f"atom {i} is bonded to nonexistent atom {j}")
            if other.edges.get(i) != order:
                raise MeltStateError(
                    f"asymmetric adjacency: {i}->{j} has order {order} but {j}->{i} "
                    f"has {other.edges.get(i)}")

    # -- hash and observables ---------------------------------------------
    def state_hash(self) -> str:
        """SHA-256 over atoms and edges in canonical (id-sorted) form."""
        digest = hashlib.sha256()
        for i in sorted(self._atoms):
            n = self._atoms[i]
            edges = ",".join(f"{j}:{o!r}" for j, o in sorted(n.edges.items()))
            digest.update(
                f"{i}|{n.element}|{n.isotope}|{n.charge}|{n.radical}|{n.lone_pairs}|{edges}\n"
                .encode("ascii"))
        return digest.hexdigest()

    def backbone_count(self, cid: int) -> int:
        """Number of atoms of the component whose origin is_backbone (an integer)."""
        return self.components.backbone_count(cid)

    def mass_moments(self, population: Callable[[int], bool]):
        """(count, sum m, sum m^2, sum m^3) over the components for which
        ``population(component id)`` is true; masses are exact Fractions."""
        count, m1, m2, m3 = 0, Fraction(0), Fraction(0), Fraction(0)
        for cid in self.components.ids():
            if population(cid):
                m = self.components.mass(cid)
                count += 1
                m1 += m
                m2 += m * m
                m3 += m * m * m
        return count, m1, m2, m3

    # -- pickling: canonical, order-independent ----------------------------
    def __getstate__(self):
        return (tuple(sorted(self._atoms.items())),
                tuple(sorted(self._origin.items())),
                self.components)

    def __setstate__(self, state):
        atoms, origin, components = state
        self._atoms = dict(atoms)
        self._origin = dict(origin)
        self._ids = frozenset(self._atoms)
        self.components = components

    # -- edits -------------------------------------------------------------
    def apply(self, edits, *, fault: Optional[Callable[[str], None]] = None) -> Journal:
        """Apply an edit set atomically.  ``fault(stage)`` is a test seam called
        at each stage of ``STAGES``; anything it raises triggers the rollback."""
        edits = tuple(edits)
        plan = self._validate(edits)
        journal = Journal(edits)
        try:
            support = sorted(plan["support"])
            self._journal_nodes(journal, support)
            touched = sorted({self.components.component_of(a) for a in support})
            journal._old_records = {c: self.components._rec[c] for c in touched}
            journal.old_components = tuple((c, self.components.atoms(c)) for c in touched)
            self._write(plan)
            if fault:
                fault("after_apply")
            self._recompute(journal, touched)
            if fault:
                fault("after_components")
            self._postcheck(journal, support, touched)
            if fault:
                fault("after_postcheck")
        except BaseException:
            self._rollback(journal)
            raise
        return journal

    def _validate(self, edits):
        bonds, fields, support = {}, {}, set()
        seen = set()
        for edit in edits:
            kind = type(edit)
            if kind is SetBond:
                a, b = edit.a, edit.b
                self._require_atoms(edit, a, b)
                if a == b:
                    raise EditError(f"{edit}: a bond needs two distinct atoms")
                for value in (edit.old, edit.new):
                    if value != 0 and float(value) not in BOND_ORDERS:
                        raise EditError(f"{edit}: bond order {value!r} is not 0 or one of {BOND_ORDERS}")
                key = ("B",) + tuple(sorted((a, b)))
                if key in seen:
                    raise EditError(f"{edit}: this bond is edited twice in one set")
                seen.add(key)
                current = self._atoms[a].edges.get(b, 0.0)
                if current != float(edit.old):
                    raise EditError(f"{edit}: expected old order {edit.old}, state has {current}")
                bonds[tuple(sorted((a, b)))] = float(edit.new)
                support.update((a, b))
            elif kind in _FIELD_OPS:
                name = _FIELD_OPS[kind]
                self._require_atoms(edit, edit.atom)
                for value in (edit.old, edit.new):
                    if not _is_int(value):
                        raise EditError(f"{edit}: {name} values must be integers")
                if name != "charge" and edit.new < 0:
                    raise EditError(f"{edit}: {name} cannot be negative")
                key = (name, edit.atom)
                if key in seen:
                    raise EditError(f"{edit}: {name} of atom {edit.atom} is edited twice in one set")
                seen.add(key)
                current = getattr(self._atoms[edit.atom], name)
                if current != edit.old:
                    raise EditError(f"{edit}: expected old {name} {edit.old}, state has {current}")
                fields[key] = edit.new
                support.add(edit.atom)
            else:
                raise EditError(f"unsupported operation {edit!r}: edits only change existing atoms")
        return {"bonds": bonds, "fields": fields, "support": support}

    def _require_atoms(self, edit, *ids):
        for i in ids:
            if not _is_int(i) or i not in self._ids:
                raise EditError(f"{edit}: atom {i!r} does not exist (atoms are never created)")

    def _journal_nodes(self, journal: Journal, ids):
        """Record the old node of every atom in the mutation support."""
        for i in ids:
            journal.old_nodes.setdefault(i, self._atoms[i])

    def _write(self, plan):
        edges = {a: dict(self._atoms[a].edges) for a in plan["support"]}
        for (a, b), order in plan["bonds"].items():
            for x, y in ((a, b), (b, a)):
                if order == 0.0:
                    edges[x].pop(y, None)
                else:
                    edges[x][y] = order
        for a in sorted(plan["support"]):
            node = self._atoms[a]
            changes = {name: value for (name, atom), value in plan["fields"].items() if atom == a}
            self._atoms[a] = replace(node, edges=_Edges(edges[a]), **changes)

    def _recompute(self, journal: Journal, touched):
        universe = [a for c in touched for a in self.components.atoms(c)]
        parts = _partition(universe, self._atoms)
        records = [_make_record(p, self._atoms, self._origin) for p in parts]
        self.components._install(touched, records)
        journal._installed = True
        journal.new_components = tuple((r.atoms[0], r.atoms) for r in records)

    def _postcheck(self, journal: Journal, support, touched):
        if self._atoms.keys() != self._ids:
            raise PostconditionError("the atom-id set changed")
        for i in support:
            self._check_node(i, self._atoms[i])
        before = sorted(a for _, atoms in journal.old_components for a in atoms)
        after = sorted(a for _, atoms in journal.new_components for a in atoms)
        if before != after:
            raise PostconditionError("recomputed components do not cover the touched atoms")
        for cid, atoms in journal.new_components:
            if cid != min(atoms) or any(self.components.component_of(a) != cid for a in atoms):
                raise PostconditionError(f"component {cid} is not indexed by its minimum atom id")

    def _rollback(self, journal: Journal):
        for i, node in journal.old_nodes.items():
            self._atoms[i] = node
        if journal._installed:
            self.components._install(
                [c for c, _ in journal.new_components], journal._old_records.values())
            journal._installed = False


# --------------------------------------------------------------------------
# deterministic materialization from a LinearHomopolymerSpec
# --------------------------------------------------------------------------
class _Builder:
    def __init__(self):
        self.atoms: Dict[int, dict] = {}
        self.origin: Dict[int, Origin] = {}

    def add_segment(self, graph, kind, ordinal, segment, backbone_keys=()):
        """Add one declared graph (heavy atoms by key, each followed by its H)."""
        ids = {}
        for decl in graph.atoms:  # GraphDeclaration sorts atoms by key
            heavy = len(self.atoms)
            ids[decl.key] = heavy
            self._new(heavy, decl.element, None if decl.isotope == -1 else decl.isotope,
                      decl.charge, decl.radical_electrons, decl.lone_pairs)
            self.origin[heavy] = Origin(kind, ordinal, segment, decl.key,
                                        decl.key in backbone_keys)
            for h in range(decl.hydrogens):
                hid = len(self.atoms)
                self._new(hid, "H", None, 0, 0, 0)
                self.origin[hid] = Origin(kind, ordinal, segment, f"{decl.key}.H{h}", False)
                self.bond(heavy, hid, 1.0)
        for bond in graph.bonds:
            self.bond(ids[bond.atoms[0]], ids[bond.atoms[1]], bond.order)
        return ids

    def _new(self, i, element, isotope, charge, radical, lone_pairs):
        self.atoms[i] = {"element": element, "isotope": isotope, "charge": charge,
                         "radical": radical, "lone_pairs": lone_pairs, "edges": {}}

    def bond(self, a, b, order):
        if b in self.atoms[a]["edges"]:
            raise MeltStateError(f"duplicate bond between atoms {a} and {b}")
        self.atoms[a]["edges"][b] = float(order)
        self.atoms[b]["edges"][a] = float(order)

    def state(self) -> MeltState:
        nodes = {i: Node(a["element"], a["isotope"], a["charge"], a["radical"],
                         a["lone_pairs"], a["edges"]) for i, a in self.atoms.items()}
        return MeltState(nodes, self.origin)


def _backbone_keys(repeat) -> frozenset:
    """Atom keys on the shortest head-port -> tail-port path of the repeat unit."""
    ports = dict(repeat.ports)
    start, goal = ports["head"].atom, ports["tail"].atom
    neighbours = {atom.key: [] for atom in repeat.atoms}
    for bond in repeat.bonds:
        a, b = bond.atoms
        neighbours[a].append(b)
        neighbours[b].append(a)
    previous = {start: None}
    frontier = [start]
    while frontier and goal not in previous:
        following = []
        for key in frontier:
            for nxt in sorted(neighbours[key]):
                if nxt not in previous:
                    previous[nxt] = key
                    following.append(nxt)
        frontier = following
    if goal not in previous:
        raise MeltStateError("repeat unit has no path between its head and tail ports")
    path, key = [], goal
    while key is not None:
        path.append(key)
        key = previous[key]
    return frozenset(path)


def materialize(spec) -> MeltState:
    """Materialize the declared melt: chains by (DP, declaration ordinal), then
    the initial mobile molecules in declaration order.  Atom ids are 0..N-1 in
    canonical order, so the same declaration gives identical ids and bytes in
    any process."""
    from rmgpy.polymer_spec import LinearHomopolymerSpec

    if not isinstance(spec, LinearHomopolymerSpec):
        raise MeltStateError(
            f"only LinearHomopolymerSpec declarations are supported, got {type(spec).__name__}")
    repeat, head, tail = spec.repeat_unit, spec.head_end_group, spec.tail_end_group
    repeat_ports = dict(repeat.ports)
    backbone = _backbone_keys(repeat)
    builder = _Builder()
    molecules = 0
    chain = 0
    for dp, count in sorted(spec.initial_ensemble.counts.items()):
        if dp < 1:
            raise MeltStateError(f"chains need DP >= 1, got {dp}")
        for _ in range(count):
            previous = builder.add_segment(head, "chain", chain, "head_end")
            joint = (previous[head.ports[0][1].atom], head.ports[0][1].bond_order)
            for unit in range(1, dp + 1):
                ids = builder.add_segment(repeat, "chain", chain, f"unit {unit}", backbone)
                builder.bond(joint[0], ids[repeat_ports["head"].atom], joint[1])
                joint = (ids[repeat_ports["tail"].atom], repeat_ports["tail"].bond_order)
            ids = builder.add_segment(tail, "chain", chain, "tail_end")
            builder.bond(joint[0], ids[tail.ports[0][1].atom], joint[1])
            chain += 1
            molecules += 1
    small = {name: graph for name, graph, _disposition in spec.small_species}
    mobile = 0
    for name, count in spec.initial_mobile_counts:
        for _ in range(count):
            builder.add_segment(small[name], "mobile", mobile, name)
            mobile += 1
            molecules += 1
    state = builder.state()
    if len(state.components) != molecules:
        raise MeltStateError(
            f"declaration yields {len(state.components)} connected molecules, expected {molecules}")
    return state
