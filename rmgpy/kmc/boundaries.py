# I-136 red commit: API stub; every algorithm raises NotImplementedError.
"""Polymer boundary providers: construction-carried port identity, exact
explicit-H cap/continuation maps, and per-port extension vectors.

A finite proxy molecule stands in for a piece of an infinite polymer.  Where
the proxy was cut out of the chain it carries an *artificial* port: a heavy
atom whose declared cap atom(s) (usually one H) stand in for the rest of the
chain.  Where the proxy ends because the polymer really ends (a chain-end
radical, a monomer) the port is *physical* and is never extended.

Pinned policies
---------------
* **Identity by construction.**  A provider builds a proxy together with its
  ``BoundaryMetadata``: per port a stable ``port_id``, the port atom id, the
  exact cap atom id(s), a tag (``backbone_head | backbone_tail | pendant |
  other``) and ``physical | artificial``; per heavy atom its repeat
  coordinate (unit, position in the unit).  Nothing is inferred later from
  graph symmetry, phenyl attachment or proxy names.
* **Witnesses are explicit-H row graphs** (``rmg_moves`` rows: element,
  isotope, charge, radical electrons, lone pairs, sorted edges) with one atom
  id namespace for both reaction sides.  ``witness_from_move`` takes a move
  from ``rmg_moves.generate_moves`` and asserts that every metadata id
  survived RMG generation (``MetadataLostError`` names the family).
* **One core, many chemistries.**  The core never assumes two ports, one
  repeat chemistry, a unique backbone path, carbon-only continuation or an H
  cap.  A provider supplies one ``Fragment`` per extension step (fresh atoms,
  their bonds, the attachment bonds to the old port atom, the new port and
  its cap); the core deletes the old port's DECLARED caps and adds the
  fragment identically on both reaction sides.
* **Exact maps.**  Every extension returns an ``ExtensionMap``: surviving
  atoms (identity map, element and isotope recorded), deleted declared caps
  (one set for both sides), fresh atoms with new ids and their bonds, the new
  port with its new cap, and the fresh repeat coordinates.  Deleting any atom
  of the edit support (atoms whose attributes or incident bonds, H included,
  differ across the reaction) is forbidden.  ``verify_extension`` re-checks a
  record against the original and extended witnesses from scratch: full
  bijection between the two sides, survivors unchanged, fresh atoms exactly
  as recorded, caps exactly the declared ones on both sides, edits
  transported unchanged, every step equal to its provider fragment.
* **Order independence.**  Fresh ids are a function of (port rank, depth,
  fragment position) only, so applying one extension vector in any port
  order yields the same witness and the same map; composite maps keep their
  steps sorted by (port, depth).
* **Vectors.**  One integer per artificial port, ports ordered by
  (participant index, port id).  ``enumerate_vectors`` walks a stated box
  under a hard heavy-atom ceiling (``PROXY_PADDING_MAX_HEAVY_ATOMS`` unless
  given) and reports ``ok | refused | ceiling_exceeded`` per vector.
* **Marked distinctness.**  Reacting atoms (the edit support, pinned by atom
  id) and ports (pinned by port key) are marks; two extended witnesses are
  duplicates only if isomorphic by a map that fixes every mark.
* **Determinism.**  No object identity, set order or hash seed reaches an
  output; inputs are keyed by atom id, never by atom list order.
"""
from __future__ import annotations
import itertools
from collections import deque
from dataclasses import dataclass, field, replace
from types import MappingProxyType
from typing import Any, Iterable, Mapping, Optional, Sequence
TAGS = ('backbone_head', 'backbone_tail', 'pendant', 'other')
PHYSICALITIES = ('physical', 'artificial')
VECTOR_STATUSES = ('ok', 'refused', 'ceiling_exceeded')
FRAGMENT_STRIDE = 1000
PORT_STRIDE = 1000000
MAX_PORTS = 1000

class BoundaryError(ValueError):
    """Base class for every refusal of this module."""

class MetadataError(BoundaryError):
    """Boundary metadata is malformed or does not fit its graph."""

class MetadataLostError(MetadataError):
    """RMG generation did not preserve a metadata atom id."""

class MapCheckError(BoundaryError):
    """An extension map failed an exactness check."""

class ExtensionRefused(BoundaryError):
    """A port may not be extended in this witness (structured refusal)."""

class CeilingExceeded(ExtensionRefused):
    """An extension would exceed the hard heavy-atom ceiling."""

@dataclass(frozen=True)
class Port:
    """One boundary of a proxy, declared at construction."""
    port_id: str
    atom_id: int
    cap_atom_ids: tuple[int, ...]
    tag: str
    physicality: str

    def __post_init__(self):
        if not isinstance(self.port_id, str) or not self.port_id:
            raise MetadataError('a port needs a non-empty string port_id')
        if self.tag not in TAGS:
            raise MetadataError(f'port tag must be one of {TAGS}, got {self.tag!r}')
        if self.physicality not in PHYSICALITIES:
            raise MetadataError(f'port physicality must be one of {PHYSICALITIES}, got {self.physicality!r}')
        caps = tuple((int(a) for a in self.cap_atom_ids))
        if len(set(caps)) != len(caps) or self.atom_id in caps:
            raise MetadataError(f'port {self.port_id}: cap atom ids must be distinct and differ from the port atom')
        if self.physicality == 'artificial' and (not caps):
            raise MetadataError(f'artificial port {self.port_id} has no declared cap')
        if self.physicality == 'physical' and caps:
            raise MetadataError(f'physical port {self.port_id} cannot carry a cap')
        object.__setattr__(self, 'cap_atom_ids', tuple(sorted(caps)))
        object.__setattr__(self, 'atom_id', int(self.atom_id))

    @property
    def artificial(self) -> bool:
        return self.physicality == 'artificial'

    def relabel(self, mapping: Mapping[int, int]) -> 'Port':
        return replace(self, atom_id=mapping[self.atom_id], cap_atom_ids=tuple((mapping[a] for a in self.cap_atom_ids)))

@dataclass(frozen=True)
class BoundaryMetadata:
    """Ports and repeat coordinates of one participant graph.

    ``coordinates`` holds ``(atom id, unit, position)`` for every heavy atom
    (sorted by atom id); ``ports`` is sorted by ``port_id``."""
    provider: str
    ports: tuple[Port, ...]
    coordinates: tuple[tuple[int, int, int], ...]

    def __post_init__(self):
        ports = tuple(sorted(self.ports, key=lambda p: p.port_id))
        if len({p.port_id for p in ports}) != len(ports):
            raise MetadataError('port ids must be unique within a participant')
        coords = tuple(sorted(((int(a), int(u), int(k)) for a, u, k in self.coordinates)))
        if len({a for a, _u, _k in coords}) != len(coords):
            raise MetadataError('an atom has two repeat coordinates')
        object.__setattr__(self, 'ports', ports)
        object.__setattr__(self, 'coordinates', coords)

    def port(self, port_id: str) -> Port:
        for port in self.ports:
            if port.port_id == port_id:
                return port
        raise KeyError(port_id)

    def coordinate(self, atom_id: int) -> tuple[int, int]:
        for a, unit, position in self.coordinates:
            if a == atom_id:
                return (unit, position)
        raise KeyError(atom_id)

    def atom_ids(self) -> frozenset:
        ids = {a for a, _u, _k in self.coordinates}
        for port in self.ports:
            ids.add(port.atom_id)
            ids.update(port.cap_atom_ids)
        return frozenset(ids)

    def relabel(self, mapping: Mapping[int, int]) -> 'BoundaryMetadata':
        return BoundaryMetadata(self.provider, tuple((p.relabel(mapping) for p in self.ports)), tuple(((mapping[a], u, k) for a, u, k in self.coordinates)))

    def check(self, rows: Mapping[int, tuple], where: str='graph', atom_ids: Optional[Iterable[int]]=None) -> None:
        """The metadata fits ``rows`` (atom id -> rmg_moves row); ``atom_ids``
        restricts the heavy-atom coverage check to one participant."""
        owned = set(rows) if atom_ids is None else set(atom_ids)
        heavy = {a for a, row in rows.items() if row[0] != 'H' and a in owned}
        coords = {a for a, _u, _k in self.coordinates}
        missing = sorted(self.atom_ids() - set(rows))
        if missing:
            raise MetadataLostError(f'{where}: metadata atom ids {missing} are absent')
        if not self.atom_ids() <= owned:
            raise MetadataError(f'{where}: metadata names atoms outside its participant')
        if coords != heavy:
            raise MetadataError(f'{where}: repeat coordinates do not cover exactly the heavy atoms (extra {sorted(coords - heavy)}, missing {sorted(heavy - coords)})')
        seen = set()
        for a, unit, position in self.coordinates:
            if (unit, position) in seen:
                raise MetadataError(f'{where}: repeat coordinate {(unit, position)} is reused')
            seen.add((unit, position))
        for port in self.ports:
            if rows[port.atom_id][0] == 'H':
                raise MetadataError(f'{where}: port {port.port_id} is a hydrogen')
            _check_cap_group(rows, port, where)

def _check_cap_group(rows, port: Port, where: str) -> None:
    """The cap atoms are bonded to the rest of the graph only through the port."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

@dataclass(frozen=True)
class ProxyBuild:
    """A provider-built proxy: the RMG molecule, its metadata, and every atom
    id in construction order (the stable relabelling key)."""
    name: str
    molecule: Any
    metadata: BoundaryMetadata
    construction_order: tuple[int, ...]

    def canonical_metadata(self) -> BoundaryMetadata:
        """The metadata with atom ids replaced by construction indices (free
        of RMG's process-global ``assign_atom_ids`` counter)."""
        return self.metadata.relabel({a: k for k, a in enumerate(self.construction_order)})

@dataclass(frozen=True)
class Fragment:
    """One repeat added at a port, in local indices.

    ``atoms``: (element, isotope, charge, radical, lone pairs) per local atom;
    ``bonds``: (i, j, order) inside the fragment; ``attachments``: (local
    index, order) bonded to the OLD port atom; ``new_port`` and ``new_caps``:
    local indices; ``coordinates``: (local index, unit, position) per heavy
    fragment atom."""
    atoms: tuple[tuple, ...]
    bonds: tuple[tuple[int, int, float], ...]
    attachments: tuple[tuple[int, float], ...]
    new_port: int
    new_caps: tuple[int, ...]
    coordinates: tuple[tuple[int, int, int], ...]

    def heavy_count(self) -> int:
        return sum((1 for atom in self.atoms if atom[0] != 'H'))

class BoundaryProvider:
    """One polymer chemistry: build proxies with metadata; extend one
    artificial port by one repeat (``fragment``)."""
    name = 'abstract'

    def build(self, *args, **kwargs) -> ProxyBuild:
        raise NotImplementedError

    def fragment(self, port: Port, coordinate: tuple[int, int]) -> Fragment:
        raise NotImplementedError

def _attrs(row) -> tuple:
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def _row_map(side) -> dict:
    """Atom id -> row from a sorted ((id, row), ...) side; duplicates refused."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def _bonds(rows) -> dict:
    """(i, j) -> order with i < j; asymmetric or dangling edges refused."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def _side_from(atoms: Mapping[int, tuple], bonds: Mapping[tuple, float]) -> tuple:
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def rows_of_molecule(molecule) -> tuple:
    """Sorted ((atom id, row), ...) of an RMG molecule (keyed by atom id)."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

@dataclass(frozen=True)
class Participant:
    """Boundary metadata of one reactant participant, in witness ids, and the
    participant's atom ids (sorted)."""
    index: int
    metadata: BoundaryMetadata
    atom_ids: tuple[int, ...]

    def __post_init__(self):
        object.__setattr__(self, 'atom_ids', tuple(sorted((int(a) for a in self.atom_ids))))

@dataclass(frozen=True)
class Witness:
    """A mapped reaction witness on explicit-H rows plus its boundaries."""
    reactants: tuple
    products: tuple
    participants: tuple[Participant, ...]
    id_base: int
    extension_counts: tuple[tuple[tuple[int, str], int], ...] = ()
    label: str = ''

    def _memo(self, name, compute):
        cached = self.__dict__.get(name)
        if cached is None:
            cached = compute()
            object.__setattr__(self, name, cached)
        return cached

    def reactant_rows(self) -> Mapping:
        """Atom id -> row of the reactant side (read-only view)."""
        return self._memo('_reactant_rows', lambda: MappingProxyType(_row_map(self.reactants)))

    def product_rows(self) -> Mapping:
        return self._memo('_product_rows', lambda: MappingProxyType(_row_map(self.products)))

    def reactant_bonds(self) -> Mapping:
        """(i, j) -> bond order (i < j) of the reactant side (read-only view)."""
        return self._memo('_reactant_bonds', lambda: MappingProxyType(_bonds(self.reactant_rows())))

    def product_bonds(self) -> Mapping:
        return self._memo('_product_bonds', lambda: MappingProxyType(_bonds(self.product_rows())))

    def port_keys(self) -> tuple[tuple[int, str], ...]:
        return tuple(sorted(((p.index, port.port_id) for p in self.participants for port in p.metadata.ports)))

    def artificial_port_keys(self) -> tuple[tuple[int, str], ...]:
        return tuple((key for key in self.port_keys() if self.port(key).artificial))

    def participant(self, index: int) -> Participant:
        for participant in self.participants:
            if participant.index == index:
                return participant
        raise KeyError(index)

    def port(self, key: tuple[int, str]) -> Port:
        return self.participant(key[0]).metadata.port(key[1])

    def depth(self, key) -> int:
        return dict(self.extension_counts).get(tuple(key), 0)

    def heavy_count(self) -> int:
        return sum((1 for _a, row in self.reactants if row[0] != 'H'))

def edit_support(witness: Witness) -> frozenset:
    """Atoms whose attributes or incident bonds (H included) differ across sides."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def _edit_support(witness: Witness) -> frozenset:
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def edits(witness: Witness) -> tuple:
    """The reaction's edits: changed atom attributes and changed bonds."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def _edits(witness: Witness) -> tuple:
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def make_witness(reactant_side, product_side, participants: Sequence[Participant], label: str='', where: str='witness') -> Witness:
    """A witness from sorted row sides; checks the bijection and metadata."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def _port_rank_keys(participants) -> tuple:
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def relabel_build(build: ProxyBuild, offset: int):
    """A copy of the proxy molecule with ids ``offset + construction index``
    and the metadata relabelled by the same map."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def witness_reactants(builds: Sequence[ProxyBuild]):
    """Relabelled reactant molecules (disjoint ids) and their participants."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def witness_from_move(move, reactant_molecules, participants, *, resonance: str='all', forms=None) -> Witness:
    """The witness of one ``rmg_moves`` move.  The reactant side is the
    resonance form the move used; metadata survival is asserted."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

@dataclass(frozen=True)
class ExtensionStep:
    """One repeat added at one port."""
    port_key: tuple[int, str]
    depth: int
    provider: str
    old_port: Port
    new_port: Port
    deleted: tuple[int, ...]
    fresh_ids: tuple[int, ...]

    def sort_key(self):
        return (self.port_key, self.depth)

@dataclass(frozen=True)
class ExtensionMap:
    """Exact record of an extension of a witness (one or several steps)."""
    steps: tuple[ExtensionStep, ...]
    surviving: tuple[tuple[int, str, Optional[int]], ...]
    deleted: tuple[int, ...]
    fresh_atoms: tuple[tuple[int, tuple], ...]
    fresh_bonds: tuple[tuple[int, int, float], ...]
    fresh_coordinates: tuple[tuple[int, int, int, int], ...]
    ports: tuple[tuple[tuple[int, str], Port], ...]

    @property
    def fresh_ids(self) -> tuple[int, ...]:
        return tuple((a for a, _attrs in self.fresh_atoms))

def _fresh_ids(witness: Witness, port_key, depth: int, size: int) -> tuple[int, ...]:
    """Order-independent fresh ids: a function of (port rank, depth, local index)."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def _forbid_support_deletion(support, deleted, where: str) -> None:
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def plan_extension(witness: Witness, port_key, provider: BoundaryProvider, support=None) -> ExtensionMap:
    """The map that extends ``port_key`` by one repeat (not yet applied)."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def _apply_side(side, deleted, fresh_atoms, fresh_bonds, where, rows=None, side_bonds=None):
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def apply_map(witness: Witness, record: ExtensionMap, *, support=None) -> Witness:
    """Apply an extension map to both sides of a witness identically."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def compose(first: ExtensionMap, second: ExtensionMap) -> ExtensionMap:
    """The single map equal to applying ``first`` then ``second``."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')
IDENTITY_MAP = ExtensionMap((), (), (), (), (), (), ())

def identity_map(witness: Witness) -> ExtensionMap:
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def _sides(witness: Witness):
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def verify_extension(original: Witness, record: ExtensionMap, extended: Witness, providers: Optional[Mapping[str, BoundaryProvider]]=None) -> None:
    """Re-check an extension map from scratch; raise ``MapCheckError``."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def _verify_steps(original, record, extended, providers) -> None:
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def _ceiling(ceiling):
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def _provider_for(witness, port_key, providers):
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def extend(witness: Witness, port_key, providers=None, *, ceiling=None, verify: bool=True, support=None):
    """Extend one artificial port by one repeat; returns (witness, map)."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def vector_ports(witness: Witness) -> tuple[tuple[int, str], ...]:
    """The vector coordinates: artificial ports ordered by (participant, port id)."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def extend_vector(witness: Witness, vector: Sequence[int], providers=None, *, order: Optional[Sequence[int]]=None, ceiling=None, verify: bool=True, verify_steps: bool=False):
    """Apply an extension vector port by port (in ``order``, default
    ascending); returns (witness, composite map).  With ``verify`` the
    composite map is verified against the original witness (which re-checks
    every step against its provider fragment) and must equal the sequential
    application; ``verify_steps`` also verifies each single step."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

@dataclass(frozen=True)
class VectorResult:
    vector: tuple[int, ...]
    status: str
    reason: str
    witness: Optional[Witness] = field(default=None, compare=False)
    record: Optional[ExtensionMap] = field(default=None, compare=False)

def _vector_plan(witness: Witness, providers=None) -> tuple[dict, dict]:
    """Per vector position: the refusal reason, or the heavy atoms one step adds."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def evaluate_vector(witness: Witness, vector: Sequence[int], providers=None, *, ceiling=None, check_order: bool=True, verify: bool=True, plan=None) -> VectorResult:
    """Status of one vector: ``refused`` (a non-zero entry on a port that may
    not be extended), ``ceiling_exceeded`` (the hard heavy-atom ceiling,
    projected before building and enforced per step) or ``ok``.  For an
    ``ok`` vector with several extended ports the reversed port order is
    applied too and must give the same witness and map."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def enumerate_vectors(witness: Witness, box, providers=None, *, ceiling=None, check_order: bool=True, verify: bool=True) -> tuple:
    """``evaluate_vector`` on every vector of the box ``0..box[i]`` (an int
    applies to every port), in lexicographic order."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def marks(witness: Witness) -> dict:
    """Pinned marks: reacting atoms by atom id, ports by port key."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def marked_graph(witness: Witness):
    """networkx graph of the whole witness: node label (element, isotope,
    attributes on both sides, mark), edge label (order before, order after)."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def _invariant(graph) -> tuple:
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def marked_isomorphic(first: Witness, second: Witness) -> bool:
    """Isomorphic by a map that fixes every mark (reacting atoms, ports)."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def marked_duplicates(items: Sequence[tuple[Any, Witness]]) -> tuple:
    """Pairs ``(key_a, key_b)`` (first occurrence first) of marked-isomorphic
    witnesses; every item is compared with the first member of each class."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def vector_duplicates(results: Sequence[VectorResult]) -> tuple:
    """Marked duplicates among the ``ok`` vectors of one enumeration."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

def _construction_metadata(provider, objects, roles, ports, name):
    """Metadata from atom objects recorded at construction.

    ``objects``: atoms in construction order (heavy atoms in declaration
    order, then H); ``roles``: heavy declaration index -> (unit, position);
    ``ports``: (port_id, heavy declaration index, tag, physicality)."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')

@dataclass(frozen=True)
class PsChainSpec:
    """A linear head-to-tail PS proxy, declared unit by unit."""
    units: int
    radical_units: tuple[int, ...] = ()
    radical_end: bool = False
    benzylic_end: bool = False
    ring_radical_unit: Optional[int] = None
    monomer: bool = False
PS_UNIT_HEAVY = 8

class PolystyreneProvider(BoundaryProvider):
    """Head-to-tail PS.  Unit u: position 0 the backbone CH2 (head side),
    position 1 the phenyl-bearing CH (tail side), positions 2..7 the ring
    (ipso first, in ring order).  Head port: unit 0 position 0; tail port:
    last unit position 1."""
    name = 'polystyrene'

    @staticmethod
    def smiles(spec: PsChainSpec) -> str:
        raise NotImplementedError('I-136 boundaries: not implemented yet')

    def build(self, name: str, spec: PsChainSpec) -> ProxyBuild:
        """Build exactly as ``ps_proxy_set`` always has (SMILES, ``update``,
        ``assign_atom_ids``) while recording each atom's role from the
        declaration: RMG keeps SMILES heavy-atom order, then the hydrogens."""
        raise NotImplementedError('I-136 boundaries: not implemented yet')

    def fragment(self, port: Port, coordinate: tuple[int, int]) -> Fragment:
        raise NotImplementedError('I-136 boundaries: not implemented yet')

def ps_structure_specs(units: int) -> dict:
    raise NotImplementedError('I-136 boundaries: not implemented yet')

class LinearAlkaneProvider(BoundaryProvider):
    """Unbranched CnH2n+2 as an explicit graph, one-carbon CH2 repeat, methyl
    ends.  Ids: carbons 0..n-1 along the chain, then each carbon's H in
    carbon order (the last H of each end carbon is its declared cap).  Head
    port carbon 0, tail port carbon n-1; carbon i has coordinate (i, 0)."""
    name = 'linear_alkane'

    def build(self, carbons: int, name: Optional[str]=None) -> ProxyBuild:
        raise NotImplementedError('I-136 boundaries: not implemented yet')

    def fragment(self, port: Port, coordinate: tuple[int, int]) -> Fragment:
        """One-atom repeat: the attached carbon IS the new port (head and tail
        coincide), so it carries the new cap itself."""
        raise NotImplementedError('I-136 boundaries: not implemented yet')
PROVIDERS = {PolystyreneProvider.name: PolystyreneProvider(), LinearAlkaneProvider.name: LinearAlkaneProvider()}

def ps_build(name: str, units: int) -> ProxyBuild:
    """Build one named ``ps_proxy_set`` structure with its metadata."""
    raise NotImplementedError('I-136 boundaries: not implemented yet')
