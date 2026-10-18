"""Chain witnesses derived from a polymer declaration alone (I-159).

``DeclaredChainProvider`` builds every chain proxy of a
``LinearHomopolymerSpec`` (``rmgpy/polymer_spec.py``) from the declared
repeat unit, end groups and monomer product; nothing here names a polymer.

Derived once per declaration (``ChainGeometry``):

* rule 2: every inter-residue bond of the joined chain ``J_3`` is a bridge
  (ladder and spiro joins are refused);
* rule 2b: no pi system of the pristine ``J_3`` contains an inter-residue
  bond (``windows._pi_atoms``/``_pi_component`` semantics);
* ``n_bb = |melt_state._backbone_keys(unit)|``; ``d_rep``: the bond distance
  between corresponding backbone atoms (the head-port atoms) of consecutive
  units in ``J_3``; ``period``: the smallest ``p`` dividing ``d_rep`` such
  that cutting the unit's backbone path every ``p`` bonds gives isomorphic
  pieces joined like the units (the shift by ``p`` is an automorphism of the
  infinite chain);
* the primitive period (one piece) is the provider's unit: its atom keys,
  sorted, give the ``position`` of a ``(unit, position)`` coordinate;
* ``phases``: one per bridge bond on the backbone path inside the primitive
  period (the places a scission can cut a unit: head and tail pieces);
* ``atom_classes``: orbits of the middle-period atoms of ``J_{2q+1}`` under
  ``site_class.extract_site_class`` keys (declared bond orders), through a
  policy adapter from the declaration's ``locality_policy``.

Builds (``ChainSpec``): primitive periods in a row, each oriented +1 (head
port toward the chain head) or -1; the first and last period optionally
trimmed to one side of a phase bond; ends ``artificial`` (H caps on the
port, port artificial), ``declared`` (the declared end group attached, port
physical) or ``radical`` (port physical, the port atom keeps the cap's
electron); then explicit, coordinate-addressed edits (``EditOp``), after
which the valence of every atom is checked.  Aromatic (1.5) declarations
are built in one deterministic Kekule form.  Nothing is silently dropped:
an unknown field, an out-of-range unit, an absent key or an edit that
breaks valence raises a ``BoundaryError`` subclass.

Coordinates are ``(unit, position)``: position = index of the key in the
period's sorted keys, plus ``width`` for a period oriented -1 (so the
orientation survives into extension and backbone indexing).  Period
detection compares pieces as species (``piece_signature``: RDKit canonical
SMILES through RMG's conversion), so implicit and explicit hydrogens and
every writing of an aromatic ring read alike while every other bond order
compares exactly.

Scope limit: closed-shell declarations only.  A repeat unit, end group or
monomer product that declares radical electrons or a multiplicity other than
1 is refused (``DeclarationRefused``, rule ``closed_shell``): witness rows
carry no multiplicity, so a declared spin could not survive extension or
template matching.  Radicals in chain builds come only from ends and edits
and take RMG's default multiplicity (radical electrons + 1).  Open-shell
declarations are a later slice that needs a multiplicity slot in the rows.

Ids follow construction order: heavy atoms (head end group, then period by
period in sorted-key order, then the tail end group), then hydrogens grouped
by their heavy atom (declared, added, end-group, then caps: the cap is the
last H created for its port)."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, fields
from math import ceil
from types import MappingProxyType
from typing import Any, Mapping, Optional

from rmgpy.kmc.boundaries import (
    PROVIDERS,
    BoundaryError,
    BoundaryMetadata,
    BoundaryProvider,
    EndEdit,
    ExtensionRefused,
    Fragment,
    MetadataError,
    Port,
    ProxyBuild,
)

END_KINDS = ("artificial", "declared", "radical")
EDIT_KINDS = ("set_radical", "set_bond_order", "delete_atom", "remove_h", "add_h")
BOND_ORDERS = (1.0, 1.5, 2.0, 3.0)
# per-element maximum valence for the site-class policy adapter (the same
# table the declaration parser uses for its overfill check)
MAX_VALENCE = {"H": 1.0, "B": 3.0, "C": 4.0, "N": 4.0, "O": 2.0, "F": 1.0, "Si": 4.0,
               "P": 5.0, "S": 6.0, "Cl": 1.0, "Br": 1.0, "I": 1.0}
SEARCH_NODES = 200_000


class DeclarationRefused(BoundaryError):
    """The declaration is outside the scope of a linear chain (rule 2: a
    non-bridge inter-residue bond; rule 2b: a pi system across residues)."""

    def __init__(self, rule: str, message: str, component=()):
        super().__init__(f"rule {rule}: {message}")
        self.rule = rule
        self.component = tuple(component)


class ChainSpecError(BoundaryError):
    """A ChainSpec or EditOp is malformed or does not fit the declaration."""


class EditError(ChainSpecError):
    """An edit names an absent site, or the edited graph breaks valence."""


# ---------------------------------------------------------------------------
# chain specifications
# ---------------------------------------------------------------------------
def _is_int(value) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


@dataclass(frozen=True)
class EditOp:
    """One explicit, coordinate-addressed edit.  ``unit`` and ``key``
    address an atom of the chain (``key``: a declaration atom key of the
    primitive period).  ``set_radical``: the atom's radical electrons become
    ``value``; ``remove_h`` / ``add_h``: ``value`` hydrogens (0 means one;
    caps are never removed); ``delete_atom``: the atom and its hydrogens;
    ``set_bond_order``: the bond to ``other`` = (unit, key) gets order
    ``value``."""

    kind: str
    unit: int
    key: str
    value: Any = 0
    other: Optional[tuple] = None

    def __post_init__(self):
        if self.kind not in EDIT_KINDS:
            raise ChainSpecError(f"unknown edit kind {self.kind!r} (known: {EDIT_KINDS})")
        if not _is_int(self.unit):
            raise ChainSpecError(f"edit {self.kind}: unit must be an integer")
        if not isinstance(self.key, str) or not self.key:
            raise ChainSpecError(f"edit {self.kind}: key must be a non-empty string")
        if self.kind == "set_bond_order":
            other = self.other
            if not isinstance(other, (tuple, list)) or len(other) != 2 or \
                    not _is_int(other[0]) or not isinstance(other[1], str):
                raise ChainSpecError("set_bond_order needs other=(unit, key)")
            object.__setattr__(self, "other", (int(other[0]), str(other[1])))
            if float(self.value) not in BOND_ORDERS:
                raise ChainSpecError(f"set_bond_order: order {self.value!r} not in {BOND_ORDERS}")
        elif self.other is not None:
            raise ChainSpecError(f"edit {self.kind} takes no other atom")
        if self.kind == "set_radical" and (not _is_int(self.value) or self.value < 0):
            raise ChainSpecError("set_radical needs a non-negative integer value")
        if self.kind in ("remove_h", "add_h") and (not _is_int(self.value) or self.value < 0):
            raise ChainSpecError(f"{self.kind} needs a non-negative integer count")
        if self.kind == "delete_atom" and self.value not in (0, None):
            raise ChainSpecError("delete_atom takes no value")

    @property
    def count(self) -> int:
        return int(self.value) or 1

    @classmethod
    def from_mapping(cls, value: Mapping) -> "EditOp":
        return _from_mapping(cls, value)


def _from_mapping(cls, value):
    if isinstance(value, cls):
        return value
    if not isinstance(value, Mapping):
        raise ChainSpecError(f"{cls.__name__} needs a mapping, got {type(value).__name__}")
    known = {f.name for f in fields(cls)}
    unknown = sorted(set(value) - known)
    if unknown:
        raise ChainSpecError(f"unknown {cls.__name__} fields: {unknown}")
    try:
        return cls(**dict(value))
    except TypeError as error:
        raise ChainSpecError(f"{cls.__name__}: {error}") from error


@dataclass(frozen=True)
class ChainSpec:
    """A chain of ``units`` primitive periods.  ``orientation``: per unit +1
    (head port toward the chain head) or -1, () = all +1.  ``head`` /
    ``tail``: an END_KINDS value.  ``head_phase`` / ``tail_phase``: 0 the
    whole period, k the k-th phase (``ChainGeometry.phases[k - 1]``): the
    first (last) period keeps only the side of the phase bond toward the
    rest of the chain.  ``edits``: EditOp, applied in order."""

    units: int
    orientation: tuple = ()
    head: str = "artificial"
    tail: str = "artificial"
    head_phase: int = 0
    tail_phase: int = 0
    edits: tuple = ()

    def __post_init__(self):
        if not _is_int(self.units) or self.units < 1:
            raise ChainSpecError(f"units must be an integer >= 1, got {self.units!r}")
        orientation = tuple(self.orientation)
        if orientation and len(orientation) != self.units:
            raise ChainSpecError(f"orientation has {len(orientation)} entries for "
                                 f"{self.units} units")
        if any(o not in (1, -1) or isinstance(o, bool) for o in orientation):
            raise ChainSpecError("orientation entries must be +1 or -1")
        object.__setattr__(self, "orientation", orientation)
        for end in (self.head, self.tail):
            if end not in END_KINDS:
                raise ChainSpecError(f"end {end!r} is not one of {END_KINDS}")
        for phase in (self.head_phase, self.tail_phase):
            if not _is_int(phase) or phase < 0:
                raise ChainSpecError(f"phase must be a non-negative integer, got {phase!r}")
        edits = tuple(EditOp.from_mapping(e) if isinstance(e, Mapping) else e
                      for e in self.edits)
        if any(not isinstance(e, EditOp) for e in edits):
            raise ChainSpecError("edits must be EditOp entries")
        object.__setattr__(self, "edits", edits)

    def orientations(self) -> tuple:
        return self.orientation or (1,) * self.units

    @classmethod
    def from_mapping(cls, value: Mapping) -> "ChainSpec":
        return _from_mapping(cls, value)


# ---------------------------------------------------------------------------
# declaration graph helpers
# ---------------------------------------------------------------------------
def _valence_electrons(element: str) -> int:
    from rmgpy.molecule.element import PeriodicSystem

    return PeriodicSystem.valence_electrons[element]


def _target_order(element, charge, radical, lone_pairs) -> float:
    """Total bond order RMG's consistency rule requires of an atom."""
    return _valence_electrons(element) - charge - radical - 2 * lone_pairs


def _isotope(declaration) -> Optional[int]:
    return None if declaration.isotope in (-1, None) else int(declaration.isotope)


def _attrs(declaration) -> tuple:
    return (declaration.element, _isotope(declaration), int(declaration.charge),
            int(declaration.radical_electrons), int(declaration.lone_pairs))


def _adjacency(atoms, bonds) -> dict:
    out = {k: {} for k in atoms}
    for (x, y), order in bonds.items():
        out[x][y] = order
        out[y][x] = order
    return out


def _bonds_of(graph) -> dict:
    return {tuple(sorted(bond.atoms)): float(bond.order) for bond in graph.bonds}


def _kekule(atoms: Mapping, bonds: Mapping, extra: Mapping) -> dict:
    """One deterministic Kekule assignment of the 1.5 bonds: every atom of
    the aromatic subgraph whose valence (aromatic bonds counted 1, plus the
    ``extra`` bond order it carries outside the graph) lacks one unit gets
    exactly one double bond."""
    aromatic = sorted(pair for pair, order in bonds.items() if order == 1.5)
    out = dict(bonds)
    if not aromatic:
        return out
    involved = sorted({k for pair in aromatic for k in pair})
    need = {}
    for key in involved:
        a = atoms[key]
        used = a.hydrogens + extra.get(key, 0.0) + sum(
            1.0 if order == 1.5 else order for pair, order in bonds.items() if key in pair)
        slack = _target_order(a.element, a.charge, a.radical_electrons, a.lone_pairs) - used
        if slack not in (0, 1):
            raise MetadataError(f"aromatic atom {key!r} cannot take a Kekule form "
                                f"(valence slack {slack})")
        need[key] = int(slack)
    partners = {k: sorted(o for pair in aromatic if k in pair for o in pair if o != k)
                for k in involved}

    def solve(matched):
        free = [k for k in involved if need[k] and k not in matched]
        if not free:
            return {}
        key = free[0]
        for other in partners[key]:
            if need[other] and other not in matched:
                rest = solve(matched | {key, other})
                if rest is not None:
                    rest[tuple(sorted((key, other)))] = 2.0
                    return rest
        return None

    doubles = solve(frozenset())
    if doubles is None:
        raise MetadataError(f"the aromatic system {involved} has no Kekule form")
    for pair in aromatic:
        out[pair] = doubles.get(pair, 1.0)
    return out


def _joins(unit) -> list:
    """(tail key, head key, order) per inter-residue join: the declared
    head/tail ports first, then any ``head:<s>``/``tail:<s>`` port pair
    (a second join between consecutive residues)."""
    ports = dict(unit.ports)
    if "head" not in ports or "tail" not in ports:
        raise MetadataError("a repeat unit needs head and tail ports")
    joins = [(ports["tail"].atom, ports["head"].atom, float(ports["tail"].bond_order))]
    for name in sorted(ports):
        if name in ("head", "tail"):
            continue
        side, _sep, suffix = name.partition(":")
        if side == "head" and suffix:
            partner = ports.get(f"tail:{suffix}")
            if partner is None:
                raise MetadataError(f"port {name!r} has no tail:{suffix} partner")
            joins.append((partner.atom, ports[name].atom, float(partner.bond_order)))
        elif side == "tail" and suffix:
            if f"head:{suffix}" not in ports:
                raise MetadataError(f"port {name!r} has no head:{suffix} partner")
        else:
            raise MetadataError(f"unknown repeat-unit port {name!r}")
    return joins


def _joined_rows(atoms, bonds, joins, copies: int):
    """J_copies of the declared unit as windows-style rows (implicit H
    omitted): (id -> row, id -> (unit, key), inter-residue id pairs)."""
    keys = sorted(atoms)
    index = {k: i for i, k in enumerate(keys)}

    def ident(u, k):
        return u * len(keys) + index[k]
    edges = {ident(u, k): [] for u in range(copies) for k in keys}
    inter = []
    for u in range(copies):
        for (x, y), order in bonds.items():
            edges[ident(u, x)].append((ident(u, y), order))
            edges[ident(u, y)].append((ident(u, x), order))
        if u + 1 < copies:
            for tail, head, order in joins:
                a, c = ident(u, tail), ident(u + 1, head)
                edges[a].append((c, order))
                edges[c].append((a, order))
                inter.append((a, c))
    rows, where = {}, {}
    for u in range(copies):
        for k in keys:
            i = ident(u, k)
            a = atoms[k]
            rows[i] = (a.element, None, a.charge, a.radical_electrons, a.lone_pairs,
                       tuple(sorted(edges[i])))
            where[i] = (u, k)
    return rows, where, inter


def _name(where, i) -> str:
    u, k = where[i]
    return f"{u}:{k}"


def _bfs(adjacency, start, allowed=None) -> dict:
    seen = {start: 0}
    frontier = [start]
    while frontier:
        following = []
        for key in frontier:
            for nxt in sorted(adjacency[key]):
                if nxt not in seen and (allowed is None or nxt in allowed):
                    seen[nxt] = seen[key] + 1
                    following.append(nxt)
        frontier = following
    return seen


def _backbone_path(adjacency, start, goal) -> list:
    """The head -> tail path ``melt_state._backbone_keys`` selects (same
    breadth-first order and tie-break), in path order."""
    previous = {start: None}
    frontier = [start]
    while frontier and goal not in previous:
        following = []
        for key in frontier:
            for nxt in sorted(adjacency[key]):
                if nxt not in previous:
                    previous[nxt] = key
                    following.append(nxt)
        frontier = following
    if goal not in previous:
        raise MetadataError("the repeat unit has no head-to-tail path")
    path, key = [], goal
    while key is not None:
        path.append(key)
        key = previous[key]
    return path[::-1]


def _components_without(adjacency, keys, removed) -> list:
    removed = {tuple(sorted(pair)) for pair in removed}
    seen, parts = set(), []
    for start in sorted(keys):
        if start in seen:
            continue
        part, stack = set(), [start]
        seen.add(start)
        while stack:
            key = stack.pop()
            part.add(key)
            for nxt in adjacency[key]:
                if nxt in keys and nxt not in seen and tuple(sorted((key, nxt))) not in removed:
                    seen.add(nxt)
                    stack.append(nxt)
        parts.append(frozenset(part))
    return parts


HEAD_CAP_ISOTOPE, TAIL_CAP_ISOTOPE = 901, 902   # signature-only port marks


def piece_signature(atoms, bonds, keys, head, tail, join_order) -> str:
    """The identity of a piece as a species: RDKit's canonical isomeric SMILES
    (RMG's own ``to_rdkit_mol`` conversion, which sanitizes and perceives
    aromaticity) of the piece in a Kekule form, every hydrogen explicit and
    each port capped with ``join_order`` hydrogens whose isotope labels mark
    head and tail.  Two writings of one aromatic ring (1.5 or any Kekule
    form, six- or five-membered, hetero or not) give one signature; every
    other bond order, charge, radical and isotope is compared exactly."""
    from rdkit import Chem
    from rmgpy.kmc.rmg_moves import molecule_from_rows

    keys = sorted(keys)
    inside = {pair: o for pair, o in bonds.items() if set(pair) <= set(keys)}
    extra = {head: join_order}
    extra[tail] = extra.get(tail, 0.0) + join_order
    asm = _Assembly()
    handles = asm.add_period(atoms, keys, {k: i for i, k in enumerate(keys)},
                             _kekule(atoms, inside, extra), 0, set(keys))
    marks = {}
    for _ in range(int(join_order)):
        marks[asm.add_h(handles[head], 4)] = HEAD_CAP_ISOTOPE
        marks[asm.add_h(handles[tail], 4)] = TAIL_CAP_ISOTOPE
    rows, ids = asm.rows()
    molecule = molecule_from_rows(rows)
    try:
        rdmol, index = molecule.to_rdkit_mol(remove_h=False, return_mapping=True)
    except (ValueError, RuntimeError) as error:
        raise MetadataError(f"piece {keys} has no RDKit form: {error}") from error
    by_id = {atom.id: index[atom] for atom in molecule.atoms}
    for handle, isotope in marks.items():
        rdmol.GetAtomWithIdx(by_id[ids[handle]]).SetIsotope(isotope)
    return Chem.MolToSmiles(rdmol, isomericSmiles=True, canonical=True)


def _is_bridge(adjacency, keys, pair) -> bool:
    return len(_components_without(adjacency, keys, [pair])) > \
        len(_components_without(adjacency, keys, []))


# ---------------------------------------------------------------------------
# site-class policy adapter
# ---------------------------------------------------------------------------
def site_class_policy(locality_policy: Mapping, max_heavy_atoms: int) -> dict:
    """``site_class`` policy fields from a declaration ``locality_policy``
    (the schemas differ): bounds map one to one, stereo is unsupported
    (``"none"`` only), the radius is the declared locality radius."""
    bounds = locality_policy["bounds"]
    radius = int(locality_policy["radius"])
    elements = sorted(bounds["elements"])
    missing = [e for e in elements if e not in MAX_VALENCE]
    if missing:
        raise MetadataError(f"no maximum valence known for elements {missing}")
    return {
        "allowed_atom_stereo": ["none"],
        "allowed_bond_orders": sorted(float(o) for o in bounds["bond_orders"]),
        "allowed_bond_stereo": ["none"],
        "allowed_charges": sorted(int(c) for c in bounds["formal_charges"]),
        "allowed_elements": elements,
        "allowed_isotopes": sorted(int(i) for i in bounds["isotopes"]),
        "allowed_lone_pairs": sorted(int(i) for i in bounds["lone_pairs"]),
        "allowed_multiplicities": sorted(int(i) for i in bounds["spin_multiplicities"]),
        "allowed_radicals": sorted(int(i) for i in bounds["radical_electrons"]),
        "canonicalizer_version": "declared_chain.v1",
        "max_degree": int(bounds["maximum_heavy_atom_degree"]),
        "max_explicit_h": int(bounds["maximum_explicit_hydrogens_per_heavy_atom"]),
        "max_heavy_atoms": int(max_heavy_atoms),
        "max_radius": radius,
        "max_root_size": int(bounds["maximum_root_tuple_size"]),
        "max_search_nodes": SEARCH_NODES,
        "max_valence_by_element": {e: MAX_VALENCE[e] for e in elements},
        "orientation": "unoriented",
        "policy_id": "declared_chain.atom_classes",
        "radius": radius,
    }


def _default_locality(unit) -> dict:
    elements = sorted({a.element for a in unit.atoms} | {"H"})
    return {"radius": 2, "bounds": {
        "elements": elements, "isotopes": sorted({a.isotope for a in unit.atoms} | {-1}),
        "formal_charges": sorted({a.charge for a in unit.atoms} | {0}),
        "radical_electrons": [0, 1], "lone_pairs": sorted({a.lone_pairs for a in unit.atoms}),
        "bond_orders": [1.0, 1.5, 2.0, 3.0], "spin_multiplicities": [1, 2, 3],
        "maximum_heavy_atom_degree": 4, "maximum_explicit_hydrogens_per_heavy_atom": 4,
        "maximum_root_tuple_size": 4}}


# ---------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Phase:
    """A bridge bond on the backbone path inside the primitive period:
    ``bond`` (head-side key, tail-side key) and the keys of each side."""

    index: int
    bond: tuple
    head_piece: tuple
    tail_piece: tuple


_GEOMETRY_CACHE: dict = {}


@dataclass(frozen=True)
class ChainGeometry:
    """Everything derived once from a declared repeat unit (module doc)."""

    digest: str
    keys: tuple                 # primitive-period keys, sorted (position = index)
    positions: Mapping          # key -> position
    atoms: Mapping              # key -> AtomDeclaration
    bonds: Mapping              # (key, key) -> build (Kekule) order, inside one period
    declared_bonds: Mapping     # (key, key) -> declared order, inside one period
    head_key: str
    tail_key: str
    join_order: float
    bridge_ok: bool
    pi_ok: bool
    largest_pi: int
    backbone_keys: frozenset    # melt_state._backbone_keys of the declared unit
    backbone_path: tuple        # primitive-period backbone path, head first
    n_bb: int
    d_rep: int
    period: int
    phases: tuple
    backbone_offsets: Mapping   # key -> distance from the period's head (shortest-path atoms)
    atom_classes: tuple         # ((site-class key, positions), ...)
    # the period is palindromic (head and tail ports swapped give the same
    # rooted graph): a reversed period equals a forward one and orientation is
    # not a discriminator.  A chain may still have a global mirror symmetry
    # that maps period boundaries elsewhere (CH2-C6H4 para): reversing one
    # period there makes a new junction, so it is not counted here.
    reflection_symmetric: bool
    unit_formula: Mapping

    @property
    def width(self) -> int:
        """Positions per orientation: a period oriented -1 (``ChainSpec``)
        takes positions ``width + index`` so the coordinate carries it."""
        return len(self.keys)

    @property
    def oriented(self) -> bool:
        """Whether a reversed period differs in its keys (head atom is not the
        tail atom)."""
        return self.head_key != self.tail_key

    @property
    def heavy_classes(self) -> int:
        return sum(1 for _key, positions in self.atom_classes
                   if self.atoms[self.keys[positions[0]]].element != "H")

    @staticmethod
    def of(spec) -> "ChainGeometry":
        _require_closed_shell(spec)
        digest = hashlib.sha256(spec.normalized_json).hexdigest()
        cached = _GEOMETRY_CACHE.get(digest)
        if cached is None:
            cached = _derive(spec.repeat_unit, spec.locality_policy, digest)
            _GEOMETRY_CACHE[digest] = cached
        return cached

    @staticmethod
    def of_unit(unit, locality_policy: Optional[Mapping] = None) -> "ChainGeometry":
        """The geometry of a bare repeat-unit graph (a ``GraphDeclaration``
        whose extra ``head:<s>``/``tail:<s>`` ports declare further joins)."""
        _closed_shell(unit, "repeat unit")
        digest = hashlib.sha256(json.dumps(unit.as_mapping(), sort_keys=True).encode()
                                ).hexdigest()
        return _derive(unit, locality_policy or _default_locality(unit), digest)


def _derive(unit, locality_policy, digest) -> ChainGeometry:
    atoms = {a.key: a for a in unit.atoms}
    unit_keys = frozenset(atoms)
    declared = _bonds_of(unit)
    adjacency = _adjacency(atoms, declared)
    joins = _joins(unit)

    # rules 2 and 2b on the pristine J_3
    import networkx as nx
    from rmgpy.kmc.windows import _pi_atoms, _pi_component

    rows, where, inter = _joined_rows(atoms, declared, joins, 3)
    graph = nx.Graph()
    graph.add_nodes_from(rows)
    for i, row in rows.items():
        for j, _o in row[5]:
            graph.add_edge(i, j)
    bridges = {tuple(sorted(e)) for e in nx.bridges(graph)}
    for pair in inter:
        if tuple(sorted(pair)) not in bridges:
            rest = graph.copy()
            rest.remove_edge(*pair)
            cycle = sorted(_name(where, i) for i in nx.shortest_path(rest, *pair))
            raise DeclarationRefused(
                "2", f"the inter-residue bond {_name(where, pair[0])}-{_name(where, pair[1])} "
                     "lies on a cycle of the joined chain (ladder or spiro)", cycle)
    pi = _pi_atoms(rows)
    largest_pi, seen = 0, set()
    for start in sorted(pi):
        if start not in seen:
            component = _pi_component(rows, pi, start)
            seen |= component
            largest_pi = max(largest_pi, len(component))
    for a, c in inter:
        if a in pi and c in pi and c in _pi_component(rows, pi, a):
            component = sorted(_name(where, i) for i in _pi_component(rows, pi, a))
            raise DeclarationRefused(
                "2b", f"the pristine pi system {component} crosses the residue bond "
                      f"{_name(where, a)}-{_name(where, c)}", component)
    if len(joins) != 1:  # pragma: no cover - two joins always close a cycle
        raise DeclarationRefused("2", "more than one join between residues")
    tail, head, join_order = joins[0]
    ports = dict(unit.ports)
    if float(ports["head"].bond_order) != join_order:
        raise MetadataError("head and tail port bond orders differ")
    if join_order not in (1.0, 2.0, 3.0):
        raise MetadataError(f"a fractional port bond order ({join_order}) cannot be capped")

    # backbone quantities
    from rmgpy.kmc.melt_state import _backbone_keys

    backbone_keys = frozenset(_backbone_keys(unit))
    n_bb = len(backbone_keys)
    distances = nx.single_source_shortest_path_length(
        graph, next(i for i, w in where.items() if w == (1, head)))
    d_rep = distances[next(i for i, w in where.items() if w == (2, head))]
    path = _backbone_path(adjacency, head, tail)
    if frozenset(path) != backbone_keys:  # pragma: no cover - same search
        raise MetadataError("backbone path differs from melt_state._backbone_keys")

    # primitive period
    period, piece, piece_tail = d_rep, unit_keys, tail
    if len(path) == d_rep:
        for p in range(1, d_rep):
            if d_rep % p:
                continue
            cuts = [(path[j * p - 1], path[j * p]) for j in range(1, d_rep // p)]
            if any(declared[tuple(sorted(c))] != join_order or
                   not _is_bridge(adjacency, unit_keys, c) for c in cuts):
                continue
            parts = _components_without(adjacency, unit_keys, cuts)
            pieces = []
            for j in range(d_rep // p):
                part = next(x for x in parts if path[j * p] in x)
                pieces.append(piece_signature(atoms, declared, part, path[j * p],
                                              path[(j + 1) * p - 1], join_order))
            if len(set(pieces)) == 1:
                period = p
                piece = next(x for x in parts if head in x)
                piece_tail = path[p - 1]
                break
    keys = tuple(sorted(piece))
    positions = {k: i for i, k in enumerate(keys)}
    piece_declared = {pair: o for pair, o in declared.items() if set(pair) <= piece}
    extra = {head: join_order}
    extra[piece_tail] = extra.get(piece_tail, 0.0) + join_order
    piece_bonds = _kekule(atoms, piece_declared, extra)
    piece_adjacency = _adjacency({k: atoms[k] for k in keys}, piece_declared)
    piece_path = tuple(path[:period])

    # phases: bridge bonds on the backbone path inside the period
    phases = []
    for x, y in zip(piece_path, piece_path[1:]):
        if _is_bridge(piece_adjacency, piece, (x, y)):
            parts = _components_without(piece_adjacency, piece, [(x, y)])
            head_side = next(p for p in parts if x in p)
            tail_side = next(p for p in parts if y in p)
            phases.append(Phase(len(phases) + 1, (x, y), tuple(sorted(head_side)),
                                tuple(sorted(tail_side))))

    # backbone offsets: atoms on any shortest head -> tail path of the period
    heavy = {k for k in keys if atoms[k].element != "H"}
    from_head = _bfs(piece_adjacency, head, heavy)
    from_tail = _bfs(piece_adjacency, piece_tail, heavy)
    span = from_head[piece_tail]
    offsets = {k: from_head[k] for k in heavy
               if k in from_head and k in from_tail and from_head[k] + from_tail[k] == span}

    reflection = piece_signature(atoms, piece_declared, piece, head, piece_tail, join_order) \
        == piece_signature(atoms, piece_declared, piece, piece_tail, head, join_order)
    formula = {}
    for a in unit.atoms:
        formula[a.element] = formula.get(a.element, 0) + 1
        if a.hydrogens:
            formula["H"] = formula.get("H", 0) + a.hydrogens
    classes = _atom_classes(atoms, keys, positions, piece_declared, head, piece_tail,
                            join_order, period, locality_policy)
    return ChainGeometry(
        digest=digest, keys=keys, positions=MappingProxyType(positions),
        atoms=MappingProxyType({k: atoms[k] for k in keys}),
        bonds=MappingProxyType(piece_bonds), declared_bonds=MappingProxyType(piece_declared),
        head_key=head, tail_key=piece_tail, join_order=join_order, bridge_ok=True,
        pi_ok=True, largest_pi=largest_pi, backbone_keys=backbone_keys,
        backbone_path=piece_path, n_bb=n_bb, d_rep=d_rep, period=period,
        phases=tuple(phases), backbone_offsets=MappingProxyType(offsets),
        atom_classes=classes, reflection_symmetric=reflection,
        unit_formula=MappingProxyType(formula))


def _closed_shell(graph, where) -> None:
    radicals = sorted(a.key for a in graph.atoms if a.radical_electrons)
    if radicals or int(graph.multiplicity) != 1:
        raise DeclarationRefused(
            "closed_shell", f"the {where} is open-shell (multiplicity {graph.multiplicity}, "
                            f"radical electrons on {radicals}); S1 builds closed-shell "
                            "declarations only", radicals)


def _require_closed_shell(spec) -> None:
    _closed_shell(spec.repeat_unit, "repeat unit")
    _closed_shell(spec.head_end_group, "head end group")
    _closed_shell(spec.tail_end_group, "tail end group")
    if spec.monomer_product is not None:
        _closed_shell(spec.monomer_product, "monomer_product")


def _atom_classes(atoms, keys, positions, bonds, head, tail, join_order, period,
                  locality_policy) -> tuple:
    """Orbits of the middle-period atoms of J_{2q+1} (declared bond orders,
    H-capped ends) under site-class keys.  q = ceil((radius + 2) / period)
    keeps every context, its hydrogens and its boundary signatures clear of
    the chain ends."""
    from rmgpy.kmc.rmg_moves import molecule_from_rows
    from rmgpy.kmc.site_class import extract_site_class

    radius = int(locality_policy["radius"])
    q = ceil((radius + 2) / period)
    copies = 2 * q + 1
    asm = _Assembly()
    handles = []
    for u in range(copies):
        handles.append(asm.add_period(atoms, keys, positions, bonds, u, set(keys)))
        if u:
            asm.bond(handles[u - 1][tail], handles[u][head], join_order)
    for _ in range(int(join_order)):
        asm.add_h(handles[0][head], slot=4)
        asm.add_h(handles[-1][tail], slot=4)
    rows, ids = asm.rows()
    molecule = molecule_from_rows(rows)
    by_id = {atom.id: atom for atom in molecule.atoms}
    heavy_count = sum(1 for atom in molecule.atoms if atom.symbol != "H")
    policy = site_class_policy(locality_policy, heavy_count)
    middle = handles[q]
    found = {}
    for key in keys:
        if atoms[key].element == "H":
            continue
        site, _mapping, _pattern = extract_site_class(
            molecule, [("site", by_id[ids[middle[key]]])], policy)
        found.setdefault(site, []).append(positions[key])
    return tuple(sorted((site, tuple(sorted(p))) for site, p in found.items()))


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------
class _Assembly:
    """A mutable explicit-H graph whose atoms carry construction sort keys:
    heavy atoms (0, seq); hydrogens (1, owner seq, slot, seq) with slots
    0 declared, 1 explicit declared H, 2 added, 3 end group, 4 cap."""

    def __init__(self):
        self.attrs = {}
        self.sort = {}
        self.coord = {}
        self.adj = {}
        self.label = {}
        self.heavy_seq = 0
        self.h_seq = 0
        self.next = 0

    def _new(self, attrs, sort, label):
        handle = self.next
        self.next += 1
        self.attrs[handle] = list(attrs)
        self.sort[handle] = sort
        self.adj[handle] = {}
        self.label[handle] = label
        return handle

    def heavy(self, attrs, label, coord=None, seq=None):
        """A heavy atom; ``seq`` (negative) places it before every period."""
        if seq is None:
            seq = self.heavy_seq
            self.heavy_seq += 1
        handle = self._new(attrs, (0, seq), label)
        if coord is not None:
            self.coord[handle] = coord
        return handle

    def add_h(self, owner, slot, attrs=("H", None, 0, 0, 0), label=None):
        self.h_seq += 1
        handle = self._new(attrs, (1, self.sort[owner][1], slot, self.h_seq),
                           label or ("H", self.label[owner]))
        self.bond(owner, handle, 1.0)
        return handle

    def bond(self, a, c, order):
        if c in self.adj[a]:
            raise MetadataError(f"duplicate bond {self.label[a]}-{self.label[c]}")
        self.adj[a][c] = float(order)
        self.adj[c][a] = float(order)

    def remove(self, handle):
        for other in list(self.adj[handle]):
            del self.adj[other][handle]
        for table in (self.attrs, self.sort, self.coord, self.adj, self.label):
            table.pop(handle, None)

    def hydrogens(self, handle):
        return sorted((h for h in self.adj[handle] if self.attrs[h][0] == "H"),
                      key=lambda h: self.sort[h])

    def add_period(self, atoms, keys, positions, bonds, unit, present, offset=0) -> dict:
        """One primitive period: heavy atoms (sorted keys), internal bonds,
        declared and explicit hydrogens; coordinates ``(unit, position +
        offset)``.  key -> handle."""
        handles = {}
        for key in keys:
            if key in present and atoms[key].element != "H":
                handles[key] = self.heavy(_attrs(atoms[key]), (unit, key),
                                          (unit, positions[key] + offset))
        for (x, y), order in sorted(bonds.items()):
            if x in handles and y in handles:
                self.bond(handles[x], handles[y], order)
        for key in keys:
            if key in handles:
                for _ in range(atoms[key].hydrogens):
                    self.add_h(handles[key], 0)
        for key in keys:
            if key in present and atoms[key].element == "H":
                owners = [o for pair in bonds for o in pair if key in pair and o != key]
                if len(owners) != 1 or owners[0] not in handles:
                    raise MetadataError(f"explicit hydrogen {key!r} needs exactly one "
                                        "heavy neighbour in the period")
                handles[key] = self.add_h(handles[owners[0]], 1, _attrs(atoms[key]),
                                          (unit, key))
        return handles

    def ordered(self) -> list:
        return sorted(self.attrs, key=lambda h: self.sort[h])

    def rows(self):
        order = self.ordered()
        ids = {h: i for i, h in enumerate(order)}
        rows = {}
        for h in order:
            edges = tuple(sorted((ids[o], order_) for o, order_ in self.adj[h].items()))
            rows[ids[h]] = tuple(self.attrs[h]) + (edges,)
        return rows, ids

    def check_valence(self, where) -> None:
        for h, attrs in self.attrs.items():
            element, _iso, charge, radical, lone_pairs = attrs
            total = sum(self.adj[h].values())
            target = _target_order(element, charge, radical, lone_pairs)
            if abs(total - target) > 1e-9:
                raise EditError(f"{where}: atom {self.label[h]} ({element}) has bond order "
                                f"{total:g} with {radical} radical electrons; its valence "
                                f"needs {target:g}")

    def check_connected(self, where) -> None:
        handles = list(self.attrs)
        seen, stack = {handles[0]}, [handles[0]]
        while stack:
            for nxt in self.adj[stack.pop()]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        if len(seen) != len(handles):
            raise EditError(f"{where}: the edited graph is disconnected")


def _canonical_digest(spec) -> str:
    from rmgpy.polymer_spec import normalized_linear_homopolymer_spec

    canonical = normalized_linear_homopolymer_spec(spec.normalized_json)
    return hashlib.sha256(canonical.normalized_json).hexdigest()


# ---------------------------------------------------------------------------
# provider
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class _Plan:
    units: int
    orientation: tuple
    entry: tuple
    exit: tuple
    present: tuple
    head_atom: str
    head_order: float
    tail_atom: str
    tail_order: float


class DeclaredChainProvider(BoundaryProvider):
    """Chain proxies of one declaration (module doc).  Registered in
    ``boundaries.PROVIDERS`` under ``declared:<sha256[:16]>`` of the
    normalized declaration."""

    def __init__(self, spec):
        self.spec = spec
        self.geometry = ChainGeometry.of(spec)
        self.name = "declared:" + _canonical_digest(spec)[:16]
        self.END_CHEMISTRIES = {"backbone_head": ("declared", "radical"),
                                "backbone_tail": ("declared", "radical")}
        PROVIDERS[self.name] = self

    RADICAL_END_CHEMISTRIES = ("radical",)

    @property
    def coordinate_backbone_atoms(self) -> int:
        """Backbone atoms per coordinate unit: one primitive period (never
        ``d_rep``)."""
        return self.geometry.period

    # -- derived quantities --------------------------------------------------
    @property
    def n_bb(self) -> int:
        return self.geometry.n_bb

    @property
    def d_rep(self) -> int:
        return self.geometry.d_rep

    @property
    def period(self) -> int:
        return self.geometry.period

    def key_of(self, position: int) -> tuple:
        """(declaration key, orientation) of a coordinate position."""
        g = self.geometry
        if not _is_int(position) or not 0 <= position < 2 * g.width:
            raise KeyError(position)
        if position >= g.width:
            return g.keys[position - g.width], -1
        return g.keys[position], 1

    def backbone_index(self, unit: int, position: int, build: Optional[ProxyBuild] = None
                       ) -> Optional[int]:
        """The backbone index of a coordinate.  With ``build``: the index on
        that built chain (``backbone_indices``: edits included).  Without:
        the pristine index, unit * period + distance from the period's
        chain-head-side end (the head atom of a +1 period, the tail atom of a
        -1 period), for atoms on a shortest head-to-tail path of the period
        (equidistant ring atoms share an index: an extent measure only), None
        otherwise; consecutive along any orientation sequence."""
        if build is not None:
            for a, u, k in build.metadata.coordinates:
                if (u, k) == (unit, position):
                    return self.backbone_indices(build).get(a)
            return None
        g = self.geometry
        try:
            key, orientation = self.key_of(position)
        except KeyError:
            return None
        offset = g.backbone_offsets.get(key)
        if offset is None:
            return None
        if orientation == -1:
            offset = g.backbone_offsets[g.tail_key] - offset
        return unit * g.period + offset

    @staticmethod
    def _path_over_rows(rows, metadata):
        """(distances from the head port, from the tail port, span) over the
        heavy atoms of ``rows`` (one participant)."""
        heavy = {a: [c for c, _o in row[5] if c in rows and rows[c][0] != "H"]
                 for a, row in rows.items() if row[0] != "H"}
        head = metadata.port("head").atom_id
        tail = metadata.port("tail").atom_id
        from_head, from_tail = _bfs(heavy, head), _bfs(heavy, tail)
        if tail not in from_head:
            raise MetadataError("the head and tail ports are not connected")
        return from_head, from_tail, from_head[tail]

    @classmethod
    def _built_path(cls, build):
        from rmgpy.kmc.boundaries import rows_of_molecule

        return cls._path_over_rows(dict(rows_of_molecule(build.molecule)), build.metadata)

    def backbone_indices(self, build: ProxyBuild) -> dict:
        """Atom id -> backbone index on the BUILT chain
        (``witness_backbone_indices`` over its rows)."""
        from rmgpy.kmc.boundaries import rows_of_molecule

        return self.witness_backbone_indices(dict(rows_of_molecule(build.molecule)),
                                             build.metadata)

    def witness_backbone_indices(self, rows, metadata: BoundaryMetadata) -> dict:
        """Atom id -> backbone index for the heavy atoms on a shortest
        head-port -> tail-port path of the participant's actual graph
        (``rows``: edits and extensions included; a deletion reroutes it): the
        head port's pristine index plus the distance from the head port.
        Consecutive along the path; atoms off every shortest path have no
        entry.  A participant without both ports (the portless monomer) uses
        ``portless_backbone_indices``."""
        if {port.port_id for port in metadata.ports} < {"head", "tail"}:
            return self.portless_backbone_indices(rows)
        from_head, from_tail, span = self._path_over_rows(rows, metadata)
        head = metadata.port("head").atom_id
        base = self.backbone_index(*metadata.coordinate(head)) or 0
        return {a: base + d for a, d in from_head.items()
                if a in from_tail and d + from_tail[a] == span}

    def portless_backbone_indices(self, rows) -> dict:
        """Backbone indices of a participant without ports (the declared
        monomer): the primitive period's heavy atoms embedded in its heavy
        atoms (a monomorphism matching elements, each image carrying at least
        the declared hydrogens), the period's backbone path at its offsets
        from the head key (a vinyl monomer: its vinyl pair, as the legacy
        monomer).
        Without an embedding: every heavy atom on a shortest path between the
        two ends of a heavy-atom diameter, indexed from one end."""
        import networkx as nx
        from networkx.algorithms import isomorphism

        g = self.geometry
        heavy = sorted(a for a, row in rows.items() if row[0] != "H")
        if not heavy:
            return {}
        adjacency = {a: sorted(c for c, _o in rows[a][5] if c in rows and rows[c][0] != "H")
                     for a in heavy}
        mol = nx.Graph()
        for a in heavy:
            mol.add_node(a, element=rows[a][0],
                         hydrogens=sum(1 for c, _o in rows[a][5] if rows[c][0] == "H"))
        mol.add_edges_from((a, c) for a in heavy for c in adjacency[a] if a < c)
        unit = nx.Graph()
        for key in g.keys:
            if g.atoms[key].element != "H":
                unit.add_node(key, element=g.atoms[key].element,
                              hydrogens=g.atoms[key].hydrogens)
        unit.add_edges_from(pair for pair in sorted(g.declared_bonds)
                            if all(k in unit for k in pair))
        matcher = isomorphism.GraphMatcher(
            mol, unit, node_match=lambda m, u: m["element"] == u["element"]
            and m["hydrogens"] >= u["hydrogens"])
        embedding = next(matcher.subgraph_monomorphisms_iter(), None)
        if embedding is not None:
            image = {key: atom for atom, key in embedding.items()}
            return {image[key]: g.backbone_offsets[key] for key in g.backbone_path}
        start = heavy[0]
        seen = _bfs(adjacency, start)
        first = min(seen, key=lambda a: (-seen[a], a))
        from_first = _bfs(adjacency, first)
        last = min(from_first, key=lambda a: (-from_first[a], a))
        from_last = _bfs(adjacency, last)
        span = from_first[last]
        return {a: d for a, d in from_first.items()
                if a in from_last and d + from_last[a] == span}

    def backbone_atoms(self, chain: ChainSpec) -> int:
        """Heavy atoms on the shortest head-port -> tail-port path of the built
        chain (the path ``backbone_indices`` indexes)."""
        return self._built_path(self.build("backbone", chain))[2] + 1

    def interior_atoms_rmg(self) -> list:
        """Perceived RMG atoms of the middle period of a 3-period artificial
        chain: its heavy atoms and their hydrogens."""
        build = self.build("interior", ChainSpec(3))
        molecule = build.molecule
        molecule.update()
        middle = {a for a, u, _k in build.metadata.coordinates if u == 1}
        atoms = [a for a in molecule.atoms if a.id in middle]
        return atoms + [h for a in atoms for h in a.bonds if h.is_hydrogen()]

    def probe_ports(self) -> tuple:
        """(port, coordinate) of the head and tail ports of a 2-period
        artificial chain."""
        build = self.build("probe", ChainSpec(2))
        return tuple((port, build.metadata.coordinate(port.atom_id))
                     for port in build.metadata.ports)

    def unit_formula(self) -> tuple:
        """One primitive period's atoms and their hydrogens, from the
        declaration: ((element, count), ...) sorted."""
        counts = {}
        for atom in self.geometry.atoms.values():
            counts[atom.element] = counts.get(atom.element, 0) + 1
            if atom.hydrogens:
                counts["H"] = counts.get("H", 0) + atom.hydrogens
        return tuple(sorted(counts.items()))

    def interior_unit_atoms(self) -> tuple:
        """(position, key) of the heavy atoms of one primitive period."""
        g = self.geometry
        return tuple((g.positions[k], k) for k in g.keys if g.atoms[k].element != "H")

    # -- builds ----------------------------------------------------------------
    def _plan(self, chain: ChainSpec) -> _Plan:
        if not isinstance(chain, ChainSpec):
            raise ChainSpecError(f"a declared chain needs a ChainSpec, got {chain!r}")
        g = self.geometry
        units = chain.units
        orientation = chain.orientations()
        entry = tuple(g.head_key if o == 1 else g.tail_key for o in orientation)
        exit_ = tuple(g.tail_key if o == 1 else g.head_key for o in orientation)
        present = [frozenset(g.keys) for _ in range(units)]
        for phase in (chain.head_phase, chain.tail_phase):
            if phase > len(g.phases):
                raise ChainSpecError(f"phase {phase} is not declared (the period has "
                                     f"{len(g.phases)} phases)")
        head_atom, head_order = entry[0], g.join_order
        tail_atom, tail_order = exit_[-1], g.join_order
        if chain.head_phase:
            ph = g.phases[chain.head_phase - 1]
            keep = ph.tail_piece if exit_[0] in ph.tail_piece else ph.head_piece
            present[0] = present[0] & frozenset(keep)
            head_atom = ph.bond[1] if ph.bond[1] in keep else ph.bond[0]
            head_order = g.bonds[tuple(sorted(ph.bond))]
        if chain.tail_phase:
            ph = g.phases[chain.tail_phase - 1]
            keep = ph.head_piece if entry[-1] in ph.head_piece else ph.tail_piece
            present[-1] = present[-1] & frozenset(keep)
            tail_atom = ph.bond[0] if ph.bond[0] in keep else ph.bond[1]
            tail_order = g.bonds[tuple(sorted(ph.bond))]
        if not present[0] or head_atom not in present[0] or tail_atom not in present[-1]:
            raise ChainSpecError("the head and tail phases leave no period between them")
        for order in (head_order, tail_order):
            if order not in (1.0, 2.0, 3.0):
                raise ChainSpecError(f"a port of bond order {order} cannot be capped")
        return _Plan(units, orientation, entry, exit_, tuple(present), head_atom, head_order,
                     tail_atom, tail_order)

    def _end_group(self, asm, side, port_handle, order, unit):
        """Attach the declared end group of ``side`` at ``port_handle``.  Its
        heavy atoms take coordinates (unit, position in its sorted keys) with
        unit -1 (head) or ``units`` (tail); a head group precedes the first
        period in construction order."""
        group = self.spec.head_end_group if side == "head" else self.spec.tail_end_group
        (_name_, port), = group.ports
        if float(port.bond_order) != order:
            raise ChainSpecError(f"the declared {side} end group binds with order "
                                 f"{port.bond_order:g}, the {side} port needs {order:g}")
        atoms = {a.key: a for a in group.atoms}
        keys = sorted(atoms)
        bonds = _kekule(atoms, _bonds_of(group), {port.atom: order})
        heavy = [k for k in keys if atoms[k].element != "H"]
        handles = {}
        for rank, key in enumerate(heavy):
            handles[key] = asm.heavy(_attrs(atoms[key]), (side, key), (unit, keys.index(key)),
                                     seq=rank - len(heavy) if side == "head" else None)
        for (x, y), o in sorted(bonds.items()):
            if x in handles and y in handles:
                asm.bond(handles[x], handles[y], o)
        if port.atom in handles:
            asm.bond(port_handle, handles[port.atom], order)
        for key in heavy:
            for _ in range(atoms[key].hydrogens):
                asm.add_h(handles[key], 0)
        for key in keys:
            if atoms[key].element != "H":
                continue
            if key == port.atom:
                owner = port_handle
            else:
                neighbours = [o for pair in bonds for o in pair if key in pair and o != key]
                if len(neighbours) != 1 or neighbours[0] not in handles:
                    raise MetadataError(f"end-group hydrogen {key!r} needs one heavy neighbour")
                owner = handles[neighbours[0]]
            handles[key] = asm.add_h(owner, 3, _attrs(atoms[key]), (side, key))

    def build(self, name: str, chain: ChainSpec) -> ProxyBuild:
        """The proxy of ``chain`` (see ChainSpec); RMG's default multiplicity
        (the declaration is closed-shell: radicals come from ends and edits)."""
        from rmgpy.kmc.rmg_moves import molecule_from_rows

        g = self.geometry
        plan = self._plan(chain)
        asm = _Assembly()
        units = plan.units
        sites, periods = {}, []
        for u in range(units):
            offset = g.width if plan.orientation[u] == -1 and g.oriented else 0
            periods.append(asm.add_period(g.atoms, g.keys, g.positions, g.bonds, u,
                                          plan.present[u], offset))
            for key, handle in periods[u].items():
                sites[(u, key)] = handle
        for u in range(units - 1):
            asm.bond(periods[u][plan.exit[u]], periods[u + 1][plan.entry[u + 1]], g.join_order)
        ports = []
        for side, end, handle, order, tag in (
                ("head", chain.head, periods[0][plan.head_atom], plan.head_order,
                 "backbone_head"),
                ("tail", chain.tail, periods[-1][plan.tail_atom], plan.tail_order,
                 "backbone_tail")):
            caps = []
            if end == "artificial":
                caps = [asm.add_h(handle, 4) for _ in range(int(order))]
            elif end == "radical":
                asm.attrs[handle][3] += int(order)
            else:
                self._end_group(asm, side, handle, order, -1 if side == "head" else units)
            ports.append((side, handle, caps, tag,
                          "artificial" if end == "artificial" else "physical"))
        port_handles = {h for _s, h, _c, _t, _p in ports}
        cap_handles = {c for _s, _h, caps, _t, _p in ports for c in caps}
        for edit in chain.edits:
            self._apply(asm, edit, units, sites, port_handles, cap_handles)
        asm.check_valence(name)
        asm.check_connected(name)
        rows, ids = asm.rows()
        built_ports = tuple(Port(side, ids[h], tuple(ids[c] for c in caps), tag, physicality)
                            for side, h, caps, tag, physicality in ports)
        coordinates = tuple((ids[h], u, k) for h, (u, k) in asm.coord.items())
        metadata = BoundaryMetadata(self.name, built_ports, coordinates)
        metadata.check(rows, name)
        return ProxyBuild(name, molecule_from_rows(rows), metadata, tuple(range(len(rows))))

    def _apply(self, asm, edit, units, sites, port_handles, cap_handles):
        g = self.geometry

        def locate(unit, key):
            if not 0 <= unit < units:
                raise EditError(f"{edit}: unit {unit} is outside the chain (0..{units - 1})")
            if key not in g.positions:
                raise EditError(f"{edit}: key {key!r} is not an atom of the period "
                                f"{list(g.keys)}")
            handle = sites.get((unit, key))
            if handle is None or handle not in asm.attrs:
                raise EditError(f"{edit}: atom ({unit}, {key!r}) is not in this chain "
                                "(trimmed by a phase or deleted)")
            return handle

        handle = locate(edit.unit, edit.key)
        if edit.kind == "set_radical":
            asm.attrs[handle][3] = int(edit.value)
        elif edit.kind == "add_h":
            for _ in range(edit.count):
                asm.add_h(handle, 2)
        elif edit.kind == "remove_h":
            free = [h for h in asm.hydrogens(handle) if h not in cap_handles]
            if len(free) < edit.count:
                raise EditError(f"{edit}: atom ({edit.unit}, {edit.key!r}) has {len(free)} "
                                f"removable hydrogens")
            for h in free[len(free) - edit.count:]:
                asm.remove(h)
        elif edit.kind == "set_bond_order":
            other = locate(*edit.other)
            if other not in asm.adj[handle]:
                raise EditError(f"{edit}: ({edit.unit}, {edit.key!r}) and {edit.other} "
                                "are not bonded")
            asm.adj[handle][other] = asm.adj[other][handle] = float(edit.value)
        elif edit.kind == "delete_atom":
            if handle in port_handles:
                raise EditError(f"{edit}: ({edit.unit}, {edit.key!r}) carries a port")
            for h in [h for h in asm.adj[handle] if asm.attrs[h][0] == "H"
                      and len(asm.adj[h]) == 1]:
                asm.remove(h)
            asm.remove(handle)

    def monomer(self, name: str = "monomer") -> ProxyBuild:
        """The declared ``monomer_product`` as a portless proxy (unit 0,
        positions from its own sorted keys)."""
        from rmgpy.kmc.rmg_moves import molecule_from_rows

        graph = self.spec.monomer_product
        if graph is None:
            raise ChainSpecError("the declaration has no monomer_product")
        atoms = {a.key: a for a in graph.atoms}
        keys = tuple(sorted(atoms))
        positions = {k: i for i, k in enumerate(keys)}
        asm = _Assembly()
        asm.add_period(atoms, keys, positions, _kekule(atoms, _bonds_of(graph), {}), 0,
                       set(keys))
        asm.check_valence(name)
        rows, ids = asm.rows()
        metadata = BoundaryMetadata(self.name, (),
                                    tuple((ids[h], u, k) for h, (u, k) in asm.coord.items()))
        metadata.check(rows, name)
        return ProxyBuild(name, molecule_from_rows(rows), metadata, tuple(range(len(rows))))

    # -- extension -------------------------------------------------------------
    def fragment(self, port: Port, coordinate: tuple) -> Fragment:
        """One primitive period at ``port``, in declared bond orders, oriented
        like the port's period (the orientation is in the coordinate): at a
        tail port (on the period's chain-tail-side atom) the new period
        follows (unit + 1), at a head port it precedes (unit - 1); it bonds
        through its port atom facing the old one and its other port atom
        carries the new caps.  A port on any other atom is a phase cut (part
        of its period is absent) and is refused."""
        g = self.geometry
        unit, position = coordinate
        try:
            key, orientation = self.key_of(position)
        except KeyError:
            raise ExtensionRefused(f"{self.name}: position {position} is not a period position")
        entry = g.head_key if orientation == 1 else g.tail_key
        exit_ = g.tail_key if orientation == 1 else g.head_key
        if port.tag == "backbone_tail":
            expected, attach, new_port, new_unit = exit_, entry, exit_, unit + 1
        elif port.tag == "backbone_head":
            expected, attach, new_port, new_unit = entry, exit_, entry, unit - 1
        else:
            raise ExtensionRefused(f"{self.name}: no extension path for a {port.tag} port")
        if key != expected:
            raise ExtensionRefused(f"{self.name}: the {port.tag} port on {key!r} is on a phase "
                                   f"cut (its period's {expected!r} end is absent); a "
                                   "phase-cut end is not extended")
        asm = _Assembly()
        offset = g.width if orientation == -1 else 0
        handles = asm.add_period(g.atoms, g.keys, g.positions, g.declared_bonds, new_unit,
                                 set(g.keys), offset)
        caps = [asm.add_h(handles[new_port], 4) for _ in range(int(g.join_order))]
        rows, ids = asm.rows()
        bonds = tuple(sorted((i, j, o) for i, row in rows.items() for j, o in row[5] if i < j))
        return Fragment(tuple(row[:5] for _i, row in sorted(rows.items())), bonds,
                        ((ids[handles[attach]], g.join_order),), ids[handles[new_port]],
                        tuple(ids[c] for c in caps),
                        tuple(sorted((ids[h], u, k) for h, (u, k) in asm.coord.items())))

    def end_edit(self, metadata: BoundaryMetadata, port: Port, end: str, rows,
                 protected) -> EndEdit:
        """An end at an artificial port, exactly as a build with that end:
        ``"declared"`` deletes the caps and attaches the whole end-group graph
        of the port's side (heavy atoms with coordinates (port unit -/+ 1,
        position in the group's sorted keys)); ``"radical"`` deletes the caps
        and puts their electrons on the port atom."""
        if end not in self.END_CHEMISTRIES.get(port.tag, ()):
            raise ExtensionRefused(f"{self.name}: {end!r} is not a declared end of a "
                                   f"{port.tag} port")
        if end in self.RADICAL_END_CHEMISTRIES:
            # the radical end: the caps go and the port atom keeps their
            # electrons, as ``build`` with ``head``/``tail`` "radical"
            return EndEdit(end, tuple(port.cap_atom_ids),
                           radicals=((port.atom_id, len(port.cap_atom_ids)),))
        side = "head" if port.tag == "backbone_head" else "tail"
        group = self.spec.head_end_group if side == "head" else self.spec.tail_end_group
        (_name_, group_port), = group.ports
        order = float(group_port.bond_order)
        if len(port.cap_atom_ids) != int(order):
            raise ExtensionRefused(f"{self.name}: the {side} end group binds with order "
                                   f"{order:g}; the port has {len(port.cap_atom_ids)} caps")
        unit, _position = metadata.coordinate(port.atom_id)
        asm = _Assembly()
        anchor = asm.heavy(("C", None, 0, 0, 0), ("anchor",))
        self._end_group(asm, side, anchor, order, unit - 1 if side == "head" else unit + 1)
        fresh = [h for h in asm.ordered() if h != anchor]
        local = {h: i for i, h in enumerate(fresh)}
        internal = tuple(sorted((local[a], local[c], o) for a in fresh
                                for c, o in asm.adj[a].items() if c != anchor and
                                local[a] < local[c]))
        attach = tuple((local[h], port.atom_id, o) for h, o in sorted(asm.adj[anchor].items()))
        return EndEdit(end, tuple(port.cap_atom_ids), tuple(tuple(asm.attrs[h]) for h in fresh),
                       attach, fresh_internal_bonds=internal,
                       fresh_coordinates=tuple((local[h], u, k) for h, (u, k)
                                               in sorted(asm.coord.items()) if h != anchor))
