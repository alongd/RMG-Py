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

Several applies can share one public transaction (``state.transaction()``);
the journal then accumulates, an exception or ``rollback()`` restores the state
byte-identically, and ``commit()`` finalizes and bumps ``state.revision`` once.

Identity for caches: ``state_hash()`` is the graph identity (atoms and edges);
``full_hash()`` also covers the immutable ancestry (origin), which determines
backbone counts and chain-length observables.  Key caches on ``full_hash()``;
``revision`` is a monotone commit counter (never decremented, never bumped by
a rollback) for ordering and staleness checks.
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
# Fault stage inside the component install, after the old records are removed
# and before the new ones are added (so the indexes are half-written).
INSTALL_STAGE = "during_install"

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


class AncestryError(MeltStateError):
    """An initially mobile monomer cannot be given a unique backbone ancestry."""


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

    @property
    def is_mobile(self) -> bool:
        """True for atoms of an initially mobile molecule (``kind == "mobile"``)."""
        return self.kind == "mobile"


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

    def mass_kg_per_mol(self, cid: int) -> Fraction:
        """Exact molecular mass in kg/mol (RMG's element-mass unit), isotope-aware."""
        return self._rec[cid].mass

    def molar_mass_g_mol(self, cid: int) -> Fraction:
        """Exact molar mass in g/mol (the unit downstream polymer code uses)."""
        return self._rec[cid].mass * 1000

    def mass(self, cid: int) -> Fraction:
        """Alias of ``mass_kg_per_mol`` (kg/mol); prefer the explicit names."""
        return self.mass_kg_per_mol(cid)

    def backbone_count(self, cid: int) -> int:
        return self._rec[cid].backbone

    # mutation (internal)
    def _install(self, remove: Iterable[int], records: Iterable[_Record], fault=None):
        for cid in remove:
            del self._rec[cid]
        if fault:
            fault(INSTALL_STAGE)
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
    """Record of a transaction: the full mutation support and old membership.

    ``status`` is ``open``, ``committed`` or ``rolled_back``.  A committed or
    rolled-back journal can never be rolled back again.
    """

    def __init__(self, edits=(), state=None):
        self.edits = tuple(edits)
        self.old_nodes: Dict[int, Node] = {}
        self.old_components: Tuple[Tuple[int, Tuple[int, ...]], ...] = ()
        self.new_components: Tuple[Tuple[int, Tuple[int, ...]], ...] = ()
        self.status = "open"
        self.revision: Optional[int] = None  # state revision after the commit
        self.radical_delta = 0
        self._state = state
        self._old_records: Dict[int, _Record] = {}
        self._affected: set = set()
        self._installing = False  # set BEFORE the component indexes are mutated


class Transaction:
    """Several applies sharing one journal; commit finalizes, rollback undoes all.

    Use ``with state.transaction() as tx:`` (commit on success, rollback on any
    exception) or ``tx = state.begin()`` with explicit ``commit()`` /
    ``rollback()``.  Nested transactions are refused.
    """

    def __init__(self, state, expected_radical_delta=None):
        self.state = state
        self.expected_radical_delta = expected_radical_delta
        self.journal = Journal((), state)
        self.journal.base_revision = state._revision
        self.open = True

    def _require_open(self):
        if not self.open:
            raise MeltStateError(f"transaction is {self.journal.status}, not open")

    def apply(self, edits, *, fault=None, expected_radical_delta=None) -> Journal:
        """Apply an edit set inside the transaction (atomic on its own: a failed
        apply undoes only itself and raises; the transaction stays open)."""
        self._require_open()
        return self.state._run(edits, fault, expected_radical_delta, self)

    def _absorb(self, call: Journal):
        j = self.journal
        for i, node in call.old_nodes.items():
            j.old_nodes.setdefault(i, node)
        for cid, record in call._old_records.items():
            if record.atoms[0] not in j._affected:
                j._old_records[cid] = record
        j._affected |= call._affected
        j._installing = True
        j.edits += call.edits
        j.radical_delta += call.radical_delta
        j.old_components = tuple((c, r.atoms) for c, r in sorted(j._old_records.items()))
        rec = self.state.components._rec
        j.new_components = tuple((a, rec[a].atoms) for a in sorted(j._affected) if a in rec)

    def commit(self) -> Journal:
        self._require_open()
        j = self.journal
        if self.expected_radical_delta is not None and j.radical_delta != self.expected_radical_delta:
            self.rollback()
            raise PostconditionError(
                f"transaction radical delta {j.radical_delta} != expected {self.expected_radical_delta}")
        self.state._revision += 1
        j.status, j.revision = "committed", self.state._revision
        self.state._tx = None
        self.open = False
        return j

    def rollback(self):
        self._require_open()
        self.state._undo(self.journal)
        self.journal.status = "rolled_back"
        self.state._tx = None
        self.open = False

    def __enter__(self):
        self._require_open()
        return self

    def __exit__(self, exc_type, exc, tb):
        if self.open:
            if exc_type is None:
                self.commit()
            else:
                self.rollback()
        return False


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
        self._revision = 0
        self._tx = None
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

    @property
    def revision(self) -> int:
        """Monotone commit counter: +1 per committed transaction (or one-shot
        apply), unchanged by rollback and by failed applies."""
        return self._revision

    def origin_kinds(self, cid: int) -> frozenset:
        """Origin kinds ("chain" / "mobile") present in a component."""
        return frozenset(self._origin[a].kind for a in self.components.atoms(cid))

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

    def full_hash(self) -> str:
        """SHA-256 over atoms, edges AND origin (ancestry).  Use this, not
        ``state_hash()``, as the cache key: backbone counts and chain-length
        observables depend on origin."""
        digest = hashlib.sha256(self.state_hash().encode("ascii"))
        for i in sorted(self._origin):
            o = self._origin[i]
            digest.update(
                f"\n{i}|{o.kind}|{o.ordinal}|{o.segment}|{o.key}|{int(o.is_backbone)}"
                .encode("utf-8"))
        return digest.hexdigest()

    def backbone_count(self, cid: int) -> int:
        """Number of atoms of the component whose origin is_backbone (an integer).

        An initially mobile molecule that is the declared ``monomer_product``
        carries backbone ancestry on the atoms its unique map onto the repeat
        unit sends to the head-tail path (``monomer_backbone_keys``); any other
        initially mobile molecule has none (0): consumers must not read zero as
        "not polymer-derived"; test ``origin[a].is_mobile`` /
        ``origin_kinds(cid)`` for that.  No DP or comparison basis is defined
        here.
        """
        return self.components.backbone_count(cid)

    def mass_moments(self, population: Callable[[int], bool], unit: str = "kg_per_mol"):
        """(count, sum m, sum m^2, sum m^3) over the components for which
        ``population(component id)`` is true; masses are exact Fractions in
        ``unit`` (``kg_per_mol`` or ``g_per_mol``)."""
        if unit not in ("kg_per_mol", "g_per_mol"):
            raise MeltStateError(f"unknown mass unit {unit!r}")
        count, m1, m2, m3 = 0, Fraction(0), Fraction(0), Fraction(0)
        for cid in self.components.ids():
            if population(cid):
                m = (self.components.mass_kg_per_mol(cid) if unit == "kg_per_mol"
                     else self.components.molar_mass_g_mol(cid))
                count += 1
                m1 += m
                m2 += m * m
                m3 += m * m * m
        return count, m1, m2, m3

    # -- pickling: canonical, order-independent ----------------------------
    def __getstate__(self):
        return (tuple(sorted(self._atoms.items())),
                tuple(sorted(self._origin.items())),
                self.components, self._revision)

    def __setstate__(self, state):
        atoms, origin, components, revision = state
        self._revision = revision
        self._tx = None
        self._atoms = dict(atoms)
        self._origin = dict(origin)
        self._ids = frozenset(self._atoms)
        self.components = components

    # -- transactions and edits ---------------------------------------------
    def begin(self, expected_radical_delta: Optional[int] = None) -> Transaction:
        """Open a transaction (refused while another is open).  If given, the
        net radical-count change of the whole transaction must equal
        ``expected_radical_delta`` at commit."""
        if self._tx is not None:
            raise MeltStateError("a transaction is already open (nested transactions are refused)")
        self._tx = Transaction(self, expected_radical_delta)
        return self._tx

    transaction = begin

    def rollback(self, journal: Journal):
        """Roll back an OPEN transaction through its journal.  Refused for a
        journal of another state, of a committed or already rolled-back
        transaction, or of anything but the currently open transaction."""
        if journal._state is not self:
            raise MeltStateError("this journal belongs to a different state")
        if journal.status != "open":
            raise MeltStateError(f"journal is {journal.status}; only an open transaction can be rolled back")
        tx = self._tx
        if tx is None or tx.journal is not journal or journal.base_revision != self._revision:
            raise MeltStateError("journal is not the open transaction of this state")
        tx.rollback()

    def apply(self, edits, *, fault: Optional[Callable[[str], None]] = None,
              expected_radical_delta: Optional[int] = None) -> Journal:
        """Apply an edit set atomically.  Outside a transaction this is a one-shot
        transaction (commits on return, +1 revision).  Inside one it accumulates
        into the transaction.  ``fault(stage)`` is a test seam called at each
        stage of ``STAGES`` and at ``INSTALL_STAGE``; anything it raises triggers
        the rollback.  ``expected_radical_delta`` (optional) must equal the net
        change in radical electrons of the touched atoms."""
        if self._tx is not None:
            return self._tx.apply(edits, fault=fault, expected_radical_delta=expected_radical_delta)
        return self._run(edits, fault, expected_radical_delta, None)

    def _run(self, edits, fault, expected, tx) -> Journal:
        edits = tuple(edits)
        plan = self._validate(edits)
        call = Journal(edits, self)
        support = sorted(plan["support"])
        charge_before = sum(self._atoms[a].charge for a in support)
        radical_before = sum(self._atoms[a].radical for a in support)
        try:
            self._journal_nodes(call, support)
            touched = sorted({self.components.component_of(a) for a in support})
            call._old_records = {c: self.components._rec[c] for c in touched}
            call._affected = {a for c in touched for a in self.components.atoms(c)}
            call.old_components = tuple((c, self.components.atoms(c)) for c in touched)
            self._write(plan)
            if fault:
                fault("after_apply")
            self._recompute(call, touched, fault)
            if fault:
                fault("after_components")
            self._postcheck(call, support, touched)
            if sum(self._atoms[a].charge for a in support) != charge_before:
                raise PostconditionError(
                    "total charge of the touched atoms is not conserved "
                    f"({charge_before} -> {sum(self._atoms[a].charge for a in support)})")
            call.radical_delta = sum(self._atoms[a].radical for a in support) - radical_before
            if expected is not None and call.radical_delta != expected:
                raise PostconditionError(
                    f"radical delta {call.radical_delta} != expected {expected}")
            if fault:
                fault("after_postcheck")
        except BaseException:
            self._undo(call)
            raise
        if tx is not None:
            tx._absorb(call)
        else:
            self._revision += 1
            call.status, call.revision = "committed", self._revision
        return call

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

    def _recompute(self, journal: Journal, touched, fault=None):
        universe = [a for c in touched for a in self.components.atoms(c)]
        parts = _partition(universe, self._atoms)
        records = [_make_record(p, self._atoms, self._origin) for p in parts]
        journal._installing = True  # before the first mutation of the indexes
        self.components._install(touched, records, fault)
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

    def _undo(self, journal: Journal):
        """Replay a journal: old nodes back, and (if the indexes may have been
        touched) every component inside the affected region replaced by its old
        record.  Idempotent, and valid for a half-written install."""
        for i, node in journal.old_nodes.items():
            self._atoms[i] = node
        if journal._installing:
            rec = self.components._rec
            for a in journal._affected:
                rec.pop(a, None)  # component ids are atom ids of their own region
            self.components._install((), journal._old_records.values())
            journal._installing = False


# --------------------------------------------------------------------------
# deterministic materialization from a LinearHomopolymerSpec
# --------------------------------------------------------------------------
class _Builder:
    def __init__(self):
        self.atoms: Dict[int, dict] = {}
        self.origin: Dict[int, Origin] = {}

    def add_segment(self, graph, kind, ordinal, segment, backbone_keys=()):
        """Add one declared graph (heavy atoms by key, each followed by its H).
        Implicit hydrogens get origin keys ``<parent key>.H<i>``; a declaration
        whose own atom keys collide with a generated key is refused (two atoms
        would share one origin key)."""
        declared = {decl.key for decl in graph.atoms}
        collide = sorted(f"{decl.key}.H{h}" for decl in graph.atoms for h in range(decl.hydrogens)
                         if f"{decl.key}.H{h}" in declared)
        if collide:
            raise MeltStateError(
                f"declared atom keys {collide} collide with generated implicit-hydrogen keys "
                "(<parent key>.H<i>); rename the declared atoms")
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


def _normalized_forms(graph, closing: bool = False):
    """(origin keys by atom id, [labelled networkx graph per resonance form],
    backbone atom ids) of one declaration materialized like the melt: every
    hydrogen explicit (declared as atoms or as counts), aromatic rings in
    every supported form (Kekule and 1.5 declarations give the same forms).
    ``closing``: ``graph`` is a repeat unit whose head-tail bond is raised by
    the port bond order (the polymerized C=C restored)."""
    import networkx as nx
    from rmgpy.kmc import window_matcher as wm
    builder = _Builder()
    backbone = _backbone_keys(graph) if closing else frozenset()
    ids = builder.add_segment(graph, "mobile", 0, "declaration", backbone)
    if closing:
        ports = dict(graph.ports)
        head, tail = ports["head"], ports["tail"]
        a, b = ids[head.atom], ids[tail.atom]
        if b not in builder.atoms[a]["edges"] or head.bond_order != tail.bond_order:
            raise AncestryError(
                "initial mobile monomer refused: the repeat unit's head and tail atoms are not bonded "
                "by one bond with equal port orders, so the monomer's C=C has no head-tail bond to map onto")
        for x, y in ((a, b), (b, a)):
            builder.atoms[x]["edges"][y] += float(head.bond_order)
    state = builder.state()
    if len(state.components) != 1:
        raise AncestryError("a monomer declaration must be one connected molecule")
    cid = state.components.ids()[0]
    graphs = []
    for form in wm.component_forms(state, cid):
        g = nx.Graph()
        for atom, row in wm.form_rows(state, cid, form).items():
            g.add_node(atom, label=row[:5])
            for other, order in row[5]:
                if atom < other:
                    g.add_edge(atom, other, order=float(order))
        graphs.append(g)
    keys = {i: o.key for i, o in state.origin.items()}
    flagged = frozenset(i for i, o in state.origin.items() if o.is_backbone)
    return keys, graphs, flagged


def _isomorphisms(g, h):
    from networkx.algorithms.isomorphism import GraphMatcher
    matcher = GraphMatcher(g, h, node_match=lambda x, y: x["label"] == y["label"],
                           edge_match=lambda x, y: x["order"] == y["order"])
    return matcher.isomorphisms_iter()


def monomer_backbone_keys(spec, graph) -> Optional[frozenset]:
    """Atom keys of ``graph`` (a mobile species declaration) that receive
    backbone origin, or None when ``graph`` is not the declared
    ``monomer_product`` (or none is declared).

    Both sides are normalized the way species identity is: every hydrogen
    explicit (declared as atoms or as counts) and resonance forms (a Kekule
    and an aromatic ring are the same molecule); the species' first form is
    matched against every form of the other side.  The repeat unit is closed
    by raising its head-tail bond by the port bond order (the vinyl C=C
    polymerized into the head-tail bond); every isomorphism of ``graph`` onto
    that closed unit sends some atoms onto the head-port -> tail-port path.
    The ancestry must be unique: all isomorphisms must flag the same atoms
    (automorphisms that leave the flagged set unchanged, such as a ring flip,
    hydrogen permutations or the two ends of ethylene, are not an ambiguity).
    Refused (``AncestryError``) when the head and tail atoms are not bonded,
    the ports differ in order, no isomorphism exists, or two isomorphisms
    flag different atoms (ambiguous)."""
    monomer = getattr(spec, "monomer_product", None)
    if monomer is None:
        return None
    keys, species_forms, _none = _normalized_forms(graph)
    query = species_forms[0]
    _k, monomer_forms, _none = _normalized_forms(monomer)
    if not any(next(_isomorphisms(query, form), None) is not None for form in monomer_forms):
        return None
    _k, closed_forms, backbone = _normalized_forms(spec.repeat_unit, closing=True)
    flagged = {frozenset(keys[k] for k, v in iso.items() if v in backbone)
               for form in closed_forms for iso in _isomorphisms(query, form)}
    if not flagged:
        raise AncestryError(
            "initial mobile monomer refused: the declared monomer_product is not isomorphic to the "
            "repeat unit with its head-tail bond closed")
    if len(flagged) != 1:
        raise AncestryError(
            "initial mobile monomer refused: its map onto the repeat unit is ambiguous "
            f"(backbone atom sets {sorted(sorted(f) for f in flagged)})")
    return next(iter(flagged))


def materialize(spec) -> MeltState:
    """Materialize the declared melt: chains by (DP, declaration ordinal), then
    the initial mobile molecules in declaration order.  Atom ids are 0..N-1 in
    canonical order, so the same declaration gives identical ids and bytes in
    any process.  An initially mobile molecule that is the declared
    ``monomer_product`` gets backbone origin through its unique map onto the
    repeat unit (``monomer_backbone_keys``; ambiguous: ``AncestryError``)."""
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
        if not count:
            continue
        ancestry = monomer_backbone_keys(spec, small[name]) or frozenset()
        for _ in range(count):
            builder.add_segment(small[name], "mobile", mobile, name, ancestry)
            mobile += 1
            molecules += 1
    state = builder.state()
    if len(state.components) != molecules:
        raise MeltStateError(
            f"declaration yields {len(state.components)} connected molecules, expected {molecules}")
    return state
