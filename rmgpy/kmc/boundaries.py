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

from collections import deque
from dataclasses import dataclass, field, replace
from types import MappingProxyType
from typing import Any, Iterable, Mapping, Optional, Sequence

TAGS = ("backbone_head", "backbone_tail", "pendant", "other")
PHYSICALITIES = ("physical", "artificial")
VECTOR_STATUSES = ("ok", "refused", "ceiling_exceeded")
# Fresh ids: id_base + port_rank * PORT_STRIDE + depth * FRAGMENT_STRIDE + local.
FRAGMENT_STRIDE = 1_000
PORT_STRIDE = 1_000_000
MAX_PORTS = 1_000


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


class BudgetExceeded(BoundaryError):
    """An enumeration would build more vectors than its budget allows."""


class CeilingExceeded(ExtensionRefused):
    """An extension would exceed the hard heavy-atom ceiling."""


# ---------------------------------------------------------------------------
# metadata
# ---------------------------------------------------------------------------
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
            raise MetadataError("a port needs a non-empty string port_id")
        if self.tag not in TAGS:
            raise MetadataError(f"port tag must be one of {TAGS}, got {self.tag!r}")
        if self.physicality not in PHYSICALITIES:
            raise MetadataError(
                f"port physicality must be one of {PHYSICALITIES}, got {self.physicality!r}")
        caps = tuple(int(a) for a in self.cap_atom_ids)
        if len(set(caps)) != len(caps) or self.atom_id in caps:
            raise MetadataError(f"port {self.port_id}: cap atom ids must be distinct "
                                "and differ from the port atom")
        if self.physicality == "artificial" and not caps:
            raise MetadataError(f"artificial port {self.port_id} has no declared cap")
        if self.physicality == "physical" and caps:
            raise MetadataError(f"physical port {self.port_id} cannot carry a cap")
        object.__setattr__(self, "cap_atom_ids", tuple(sorted(caps)))
        object.__setattr__(self, "atom_id", int(self.atom_id))

    @property
    def artificial(self) -> bool:
        return self.physicality == "artificial"

    def relabel(self, mapping: Mapping[int, int]) -> "Port":
        return replace(self, atom_id=mapping[self.atom_id],
                       cap_atom_ids=tuple(mapping[a] for a in self.cap_atom_ids))


@dataclass(frozen=True)
class BoundaryMetadata:
    """Ports and repeat coordinates of one participant graph.

    ``coordinates`` holds ``(atom id, unit, position)`` for every heavy atom
    (sorted by atom id); ``ports`` is sorted by ``port_id``.
    ``monomer_backbone``: None, or (on a portless declared monomer) the atom
    ids of its backbone in path order, bound once from the declaration
    (empty: the declaration defines none); carried by atom id.
    ``residues``: None, or (on a residue-sequence chain) the residue letter
    of every coordinate unit, in unit order."""

    provider: str
    ports: tuple[Port, ...]
    coordinates: tuple[tuple[int, int, int], ...]
    monomer_backbone: Optional[tuple] = field(default=None, repr=False)
    residues: Optional[tuple] = field(default=None, repr=False)

    def __post_init__(self):
        if self.residues is not None:
            residues = tuple(self.residues)
            if not all(isinstance(r, str) and r for r in residues):
                raise MetadataError("residue letters must be non-empty strings")
            object.__setattr__(self, "residues", residues)
        if self.monomer_backbone is not None:
            if self.ports:
                raise MetadataError("a monomer backbone annotation needs a portless participant")
            object.__setattr__(self, "monomer_backbone",
                               tuple(int(a) for a in self.monomer_backbone))
        ports = tuple(sorted(self.ports, key=lambda p: p.port_id))
        if len({p.port_id for p in ports}) != len(ports):
            raise MetadataError("port ids must be unique within a participant")
        coords = tuple(sorted((int(a), int(u), int(k)) for a, u, k in self.coordinates))
        if len({a for a, _u, _k in coords}) != len(coords):
            raise MetadataError("an atom has two repeat coordinates")
        object.__setattr__(self, "ports", ports)
        object.__setattr__(self, "coordinates", coords)

    def port(self, port_id: str) -> Port:
        for port in self.ports:
            if port.port_id == port_id:
                return port
        raise KeyError(port_id)

    def coordinate(self, atom_id: int) -> tuple[int, int]:
        for a, unit, position in self.coordinates:
            if a == atom_id:
                return unit, position
        raise KeyError(atom_id)

    def atom_ids(self) -> frozenset:
        ids = {a for a, _u, _k in self.coordinates}
        for port in self.ports:
            ids.add(port.atom_id)
            ids.update(port.cap_atom_ids)
        return frozenset(ids)

    def carried(self, ports, coordinates, deleted=frozenset()) -> "BoundaryMetadata":
        """The same participant's metadata with new ports and coordinates:
        every other field kept (``dataclasses.replace``), the monomer
        backbone annotation minus the ``deleted`` atom ids."""
        stamp = self.monomer_backbone
        if stamp is not None:
            stamp = tuple(a for a in stamp if a not in deleted)
        return replace(self, ports=tuple(ports), coordinates=tuple(coordinates),
                       monomer_backbone=stamp)

    def relabel(self, mapping: Mapping[int, int]) -> "BoundaryMetadata":
        return replace(
            self, ports=tuple(p.relabel(mapping) for p in self.ports),
            coordinates=tuple((mapping[a], u, k) for a, u, k in self.coordinates),
            monomer_backbone=None if self.monomer_backbone is None else
            tuple(mapping[a] for a in self.monomer_backbone))

    def check(self, rows: Mapping[int, tuple], where: str = "graph",
              atom_ids: Optional[Iterable[int]] = None) -> None:
        """The metadata fits ``rows`` (atom id -> rmg_moves row); ``atom_ids``
        restricts the heavy-atom coverage check to one participant."""
        owned = set(rows) if atom_ids is None else set(atom_ids)
        heavy = {a for a, row in rows.items() if row[0] != "H" and a in owned}
        coords = {a for a, _u, _k in self.coordinates}
        missing = sorted(self.atom_ids() - set(rows))
        if missing:
            raise MetadataLostError(f"{where}: metadata atom ids {missing} are absent")
        if not self.atom_ids() <= owned:
            raise MetadataError(f"{where}: metadata names atoms outside its participant")
        if coords != heavy:
            raise MetadataError(f"{where}: repeat coordinates do not cover exactly the "
                                f"heavy atoms (extra {sorted(coords - heavy)}, "
                                f"missing {sorted(heavy - coords)})")
        if self.monomer_backbone and not set(self.monomer_backbone) <= coords:
            raise MetadataError(f"{where}: monomer backbone atoms "
                                f"{sorted(set(self.monomer_backbone) - coords)} carry no "
                                "repeat coordinate")
        seen = set()
        for a, unit, position in self.coordinates:
            if (unit, position) in seen:
                raise MetadataError(f"{where}: repeat coordinate {(unit, position)} is reused")
            seen.add((unit, position))
        for port in self.ports:
            if rows[port.atom_id][0] == "H":
                raise MetadataError(f"{where}: port {port.port_id} is a hydrogen")
            _check_cap_group(rows, port, where)


def _check_cap_group(rows, port: Port, where: str) -> None:
    """The cap atoms are bonded to the rest of the graph only through the port."""
    caps = set(port.cap_atom_ids)
    attachment = []
    for cap in sorted(caps):
        for neighbour, order in rows[cap][5]:
            if neighbour in caps:
                continue
            if neighbour != port.atom_id:
                raise MetadataError(f"{where}: cap {cap} of port {port.port_id} is bonded "
                                    f"to {neighbour}, not only to the port")
            attachment.append(order)
    if port.cap_atom_ids and not attachment:
        raise MetadataError(f"{where}: cap of port {port.port_id} is not bonded to its port")


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
        return sum(1 for atom in self.atoms if atom[0] != "H")


class BoundaryProvider:
    """One polymer chemistry: build proxies with metadata; extend one
    artificial port by one repeat (``fragment``).

    The window machinery (``windows``) reads a provider only through this
    interface (I-165): ``coordinate_backbone_atoms`` (backbone atoms per
    coordinate unit; None: no backbone, no span check), ``backbone_index`` /
    ``witness_backbone_indices`` (backbone positions of a participant),
    ``RADICAL_END_CHEMISTRIES`` (the end chemistries that put an unpaired
    electron on the port atom), ``interior_atoms_rmg`` (perceived RMG atoms
    of an interior repeat, hydrogens included), ``probe_ports`` (port and
    coordinate pairs whose fragments give the extension facts),
    ``unit_formula`` (the depository growth unit; None: no chain) and
    ``growth_sites`` (where one unit inserts in a bare molecule; None: the
    provider cannot say)."""

    name = "abstract"
    coordinate_backbone_atoms: Optional[int] = None
    RADICAL_END_CHEMISTRIES: tuple = ()
    ALKENE_END_CHEMISTRIES: tuple = ()

    def build(self, *args, **kwargs) -> ProxyBuild:  # pragma: no cover - interface
        raise NotImplementedError

    def fragment(self, port: Port, coordinate: tuple[int, int]) -> Fragment:  # pragma: no cover
        raise NotImplementedError

    def backbone_index(self, unit: int, position: int) -> Optional[int]:
        return None

    def witness_backbone_indices(self, rows, metadata: BoundaryMetadata) -> dict:
        """Atom id -> backbone index of the participant's atoms in ``rows``
        (by repeat coordinate: ``backbone_index``)."""
        out = {}
        for a, unit, position in metadata.coordinates:
            if a in rows:
                index = self.backbone_index(unit, position)
                if index is not None:
                    out[a] = index
        return out

    def interior_atoms_rmg(self) -> list:
        return []

    def probe_ports(self) -> tuple:
        return ()

    def unit_formula(self) -> Optional[tuple]:
        return None

    def growth_sites(self, molecule) -> Optional[tuple]:
        """((Port, coordinate), ...): the artificial ends of ``molecule`` (its
        atoms carry ids) at which ``fragment`` inserts one unit; None when the
        provider cannot say where a unit inserts in a molecule without
        metadata."""
        return None

    def reference_special_declarations(self) -> tuple:
        """Extra ``(label, dp, head, tail, interior_radical)`` chain phases
        needed to represent a provider's shorter-than-one-repeat species."""
        return ()

    def spec_for_declared_provider(self, declared, dp, head=None, tail=None,
                                   interior_radical=None):
        """Translate this provider's reference vocabulary to ``declared``.

        Closed and radical end names have declaration-owned equivalents.
        Providers with phase-shifted or unsaturated ends override this adapter.
        """
        radicals = set(self.RADICAL_END_CHEMISTRIES)
        alkenes = set(self.ALKENE_END_CHEMISTRIES)

        def end_kind(end, port):
            if end is None or end in ("declared", "radical"):
                return end
            allowed = set(self.END_CHEMISTRIES.get(port, ()))
            if end not in allowed:
                raise ValueError(f"{end!r} is not a declared {port} end "
                                 f"(declared: {sorted(allowed)})")
            if end in radicals:
                return "radical"
            if end in alkenes:
                raise ExtensionRefused(f"declared provider has no adapter for alkene end {end!r}")
            return "declared"

        return declared.spec_from_declaration(
            dp, end_kind(head, "backbone_head"), end_kind(tail, "backbone_tail"),
            interior_radical)

    def reference_monomer(self, spec):
        """Build the explicitly declared addition monomer, or derive it by
        raising the repeat-unit port bond when no monomer is declared."""
        from rmgpy.kmc.declared_chain import bind_monomer, _declaration_molecule
        if spec.monomer_product is not None:
            binding = bind_monomer(spec.repeat_unit, spec.monomer_product)
            if binding.keys is None:
                raise ExtensionRefused(binding.refusal)
            return _declaration_molecule(spec.monomer_product)[0]
        ports = dict(spec.repeat_unit.ports)
        built = _declaration_molecule(
            spec.repeat_unit, (ports["head"].atom, ports["tail"].atom))
        return None if built is None else built[0]


def explicit_species_build(name, *, smiles=None, adjacency_list=None):
    """Build an explicit, port-free species while preserving declared spin."""
    from rmgpy.molecule.molecule import Molecule
    molecule = Molecule(smiles=smiles) if smiles is not None else \
        Molecule().from_adjacency_list(adjacency_list)
    multiplicity = molecule.multiplicity
    objects = list(molecule.atoms)
    molecule.update()
    molecule.multiplicity = multiplicity
    molecule.assign_atom_ids()
    rows = dict(rows_of_molecule(molecule))
    heavy = sorted(a for a, row in rows.items() if row[0] != "H")
    metadata = BoundaryMetadata("small", (),
                                tuple((a, k, 0) for k, a in enumerate(heavy)))
    metadata.check(rows, name)
    return ProxyBuild(name, molecule, metadata, tuple(atom.id for atom in objects))


# the two extension probes of the named providers: a tail port on unit 0
# position 1 and a head port on unit 0 position 0 (``windows`` facts)
LEGACY_PROBE_PORTS = ((Port("p", 0, (1,), "backbone_tail", "artificial"), (0, 1)),
                      (Port("p", 0, (1,), "backbone_head", "artificial"), (0, 0)))


# ---------------------------------------------------------------------------
# rows
# ---------------------------------------------------------------------------
def _attrs(row) -> tuple:
    return tuple(row[:5])


def _row_map(side) -> dict:
    """Atom id -> row from a sorted ((id, row), ...) side; duplicates refused."""
    rows = {}
    for atom_id, row in side:
        if atom_id in rows:
            raise MapCheckError(f"atom id {atom_id} occurs twice on one reaction side")
        rows[atom_id] = row
    return rows


def _bonds(rows) -> dict:
    """(i, j) -> order with i < j; asymmetric or dangling edges refused."""
    directed = {}
    for a, row in rows.items():
        for b, order in row[5]:
            if (a, b) in directed:
                raise MapCheckError(f"atom {a} lists neighbour {b} twice")
            directed[(a, b)] = order
    bonds = {}
    for (a, b), order in directed.items():
        if b not in rows:
            raise MapCheckError(f"atom {a} is bonded to absent atom {b}")
        if directed.get((b, a)) != order:
            raise MapCheckError(f"bond {a}-{b} is not symmetric")
        if a < b:
            bonds[(a, b)] = float(order)
    return bonds


def _side_from(atoms: Mapping[int, tuple], bonds: Mapping[tuple, float]) -> tuple:
    edges = {a: [] for a in atoms}
    for (a, b), order in bonds.items():
        edges[a].append((b, float(order)))
        edges[b].append((a, float(order)))
    return tuple((a, tuple(atoms[a]) + (tuple(sorted(edges[a])),)) for a in sorted(atoms))


def rows_of_molecule(molecule) -> tuple:
    """Sorted ((atom id, row), ...) of an RMG molecule (keyed by atom id)."""
    from rmgpy.kmc.rmg_moves import molecule_rows

    return tuple(sorted(molecule_rows(molecule).items()))


# ---------------------------------------------------------------------------
# witnesses
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Participant:
    """Boundary metadata of one reactant participant, in witness ids, and the
    participant's atom ids (sorted)."""

    index: int
    metadata: BoundaryMetadata
    atom_ids: tuple[int, ...]

    def __post_init__(self):
        object.__setattr__(self, "atom_ids", tuple(sorted(int(a) for a in self.atom_ids)))


@dataclass(frozen=True)
class Witness:
    """A mapped reaction witness on explicit-H rows plus its boundaries."""

    reactants: tuple
    products: tuple
    participants: tuple[Participant, ...]
    id_base: int
    extension_counts: tuple[tuple[tuple[int, str], int], ...] = ()
    label: str = ""

    def _memo(self, name, compute):
        # derived views are cached on the (immutable) instance; they are not
        # fields, so equality, hashing and repr ignore them
        cached = self.__dict__.get(name)
        if cached is None:
            cached = compute()
            object.__setattr__(self, name, cached)
        return cached

    def reactant_rows(self) -> Mapping:
        """Atom id -> row of the reactant side (read-only view)."""
        return self._memo("_reactant_rows",
                          lambda: MappingProxyType(_row_map(self.reactants)))

    def product_rows(self) -> Mapping:
        return self._memo("_product_rows", lambda: MappingProxyType(_row_map(self.products)))

    def reactant_bonds(self) -> Mapping:
        """(i, j) -> bond order (i < j) of the reactant side (read-only view)."""
        return self._memo("_reactant_bonds",
                          lambda: MappingProxyType(_bonds(self.reactant_rows())))

    def product_bonds(self) -> Mapping:
        return self._memo("_product_bonds",
                          lambda: MappingProxyType(_bonds(self.product_rows())))

    def port_keys(self) -> tuple[tuple[int, str], ...]:
        return tuple(sorted((p.index, port.port_id) for p in self.participants
                            for port in p.metadata.ports))

    def artificial_port_keys(self) -> tuple[tuple[int, str], ...]:
        return tuple(key for key in self.port_keys() if self.port(key).artificial)

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
        return sum(1 for _a, row in self.reactants if row[0] != "H")


def edit_support(witness: Witness) -> frozenset:
    """Atoms whose attributes or incident bonds (H included) differ across sides."""
    return witness._memo("_edit_support", lambda: _edit_support(witness))


def _edit_support(witness: Witness) -> frozenset:
    before, after = witness.reactant_rows(), witness.product_rows()
    support = set()
    for atom_id in set(before) | set(after):
        if atom_id not in before or atom_id not in after:
            raise MapCheckError(f"atom {atom_id} is on one reaction side only")
        if before[atom_id] != after[atom_id]:
            support.add(atom_id)
    return frozenset(support)


def edits(witness: Witness) -> tuple:
    """The reaction's edits: changed atom attributes and changed bonds."""
    return witness._memo("_edits", lambda: _edits(witness))


def _edits(witness: Witness) -> tuple:
    before, after = witness.reactant_rows(), witness.product_rows()
    bonds_before, bonds_after = witness.reactant_bonds(), witness.product_bonds()
    atoms = tuple(sorted((a, _attrs(before[a]), _attrs(after[a])) for a in before
                         if _attrs(before[a]) != _attrs(after[a])))
    bonds = tuple(sorted((pair, bonds_before.get(pair, 0.0), bonds_after.get(pair, 0.0))
                         for pair in set(bonds_before) | set(bonds_after)
                         if bonds_before.get(pair, 0.0) != bonds_after.get(pair, 0.0)))
    return atoms, bonds


def make_witness(reactant_side, product_side, participants: Sequence[Participant],
                 label: str = "", where: str = "witness") -> Witness:
    """A witness from sorted row sides; checks the bijection and metadata."""
    reactants, products = tuple(reactant_side), tuple(product_side)
    before, after = _row_map(reactants), _row_map(products)
    _bonds(before), _bonds(after)
    if set(before) != set(after):
        raise MapCheckError(f"{where}: reactant and product atom ids differ")
    for atom_id in before:
        if before[atom_id][:2] != after[atom_id][:2]:
            raise MapCheckError(f"{where}: atom {atom_id} changes element or isotope")
    participants = tuple(sorted(participants, key=lambda p: p.index))
    if len({p.index for p in participants}) != len(participants):
        raise MetadataError(f"{where}: two participants share an index")
    owned = set()
    for participant in participants:
        ids = set(participant.atom_ids)
        if owned & ids:
            raise MetadataError(f"{where}: two participants claim one atom id")
        owned |= ids
        missing = sorted(participant.metadata.atom_ids() - set(after))
        if missing:
            raise MetadataLostError(f"{where}: metadata atoms {missing} are absent from "
                                    "the products")
        participant.metadata.check(before, f"{where} reactants (participant "
                                   f"{participant.index})", participant.atom_ids)
    if participants and owned != set(before):
        raise MetadataError(f"{where}: participants do not cover the witness atoms")
    if len(_port_rank_keys(participants)) > MAX_PORTS:
        raise MetadataError("too many ports")
    id_base = max(before, default=-1) + 1
    return Witness(reactants, products, participants, id_base, (), label)


def _port_rank_keys(participants) -> tuple:
    return tuple(sorted((p.index, port.port_id) for p in participants
                        for port in p.metadata.ports))


def relabel_build(build: ProxyBuild, offset: int):
    """A copy of the proxy molecule with ids ``offset + construction index``
    and the metadata relabelled by the same map."""
    mapping = {old: offset + k for k, old in enumerate(build.construction_order)}
    if len(mapping) != len(build.construction_order):
        raise MetadataError(f"{build.name}: construction order repeats an atom id")
    molecule = build.molecule.copy(deep=True)
    ids = [atom.id for atom in build.molecule.atoms]
    if sorted(ids) != sorted(mapping):
        raise MetadataError(f"{build.name}: construction order does not cover the molecule")
    for atom, old in zip(molecule.atoms, ids):
        atom.id = mapping[old]
    return molecule, build.metadata.relabel(mapping)


def witness_reactants(builds: Sequence[ProxyBuild]):
    """Relabelled reactant molecules (disjoint ids) and their participants."""
    molecules, participants, offset = [], [], 0
    for index, build in enumerate(builds):
        molecule, metadata = relabel_build(build, offset)
        offset += len(build.construction_order)
        molecules.append(molecule)
        participants.append(Participant(index, metadata,
                                        tuple(atom.id for atom in molecule.atoms)))
    return tuple(molecules), tuple(participants)


def witness_from_move(move, reactant_molecules, participants, *, resonance: str = "all",
                      forms=None) -> Witness:
    """The witness of one ``rmg_moves`` move.  The reactant side is the
    resonance form the move used; metadata survival is asserted."""
    from rmgpy.kmc import rmg_moves as rm

    if move.products is None:
        raise BoundaryError(f"a {move.status} move has no product")
    if forms is None:
        forms = [rm.reactant_forms(m, resonance) for m in reactant_molecules]
    reactant_rows = {}
    for k, molecule_forms in enumerate(forms):
        reactant_rows.update(rm.molecule_rows(molecule_forms[move.forms[k]]))
    label = f"{move.family}/{move.direction}/{move.order}/{move.forms}/{move.labels}"
    where = f"{move.family} ({move.direction})"
    product_rows = dict(move.products)
    if set(product_rows) != set(reactant_rows):
        raise MetadataLostError(f"{where}: RMG generation did not preserve the atom ids")
    try:
        return make_witness(tuple(sorted(reactant_rows.items())),
                            tuple(sorted(product_rows.items())), participants, label, where)
    except MetadataLostError as error:
        raise MetadataLostError(f"{where}: {error}") from error


# ---------------------------------------------------------------------------
# extension maps
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ExtensionStep:
    """One repeat added at one port."""

    port_key: tuple[int, str]
    depth: int
    provider: str
    old_port: Port
    new_port: Port
    deleted: tuple[int, ...]
    fresh_ids: tuple[int, ...]  # in fragment-local order

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
    fresh_coordinates: tuple[tuple[int, int, int, int], ...]  # (atom, participant, unit, pos)
    ports: tuple[tuple[tuple[int, str], Port], ...]

    @property
    def fresh_ids(self) -> tuple[int, ...]:
        return tuple(a for a, _attrs in self.fresh_atoms)


def _fresh_ids(witness: Witness, port_key, depth: int, size: int) -> tuple[int, ...]:
    """Order-independent fresh ids: a function of (port rank, depth, local index)."""
    if size > FRAGMENT_STRIDE:
        raise BoundaryError("fragment is larger than the fresh-id stride")
    if depth * FRAGMENT_STRIDE + size > PORT_STRIDE:
        raise BoundaryError("extension depth exceeds the fresh-id stride")
    rank = witness.port_keys().index(tuple(port_key))
    start = witness.id_base + rank * PORT_STRIDE + depth * FRAGMENT_STRIDE
    return tuple(range(start, start + size))


def _forbid_support_deletion(support, deleted, where: str) -> None:
    hit = sorted(set(deleted) & set(support))
    if hit:
        raise MapCheckError(f"{where}: deleting edit-support atoms {hit} is forbidden")


def plan_extension(witness: Witness, port_key, provider: BoundaryProvider,
                   support=None) -> ExtensionMap:
    """The map that extends ``port_key`` by one repeat (not yet applied)."""
    port_key = tuple(port_key)
    if port_key not in witness.port_keys():
        raise BoundaryError(f"unknown port {port_key}")
    port = witness.port(port_key)
    if not port.artificial:
        raise ExtensionRefused(f"port {port_key} is physical")
    metadata = witness.participant(port_key[0]).metadata
    if provider.name != metadata.provider:
        raise BoundaryError(f"port {port_key} belongs to provider {metadata.provider!r}, "
                            f"not {provider.name!r}")
    if support is None:
        support = edit_support(witness)
    if port.atom_id in support:
        raise ExtensionRefused(f"port {port_key} atom {port.atom_id} is a reacting atom")
    caps = tuple(port.cap_atom_ids)
    reactant_rows = witness.reactant_rows()
    _forbid_support_deletion(support, caps, f"port {port_key}")
    fragment = provider.fragment(port, metadata.coordinate(port.atom_id))
    depth = witness.depth(port_key)
    ids = _fresh_ids(witness, port_key, depth, len(fragment.atoms))
    bonds = [(ids[i], ids[j], float(o)) for i, j, o in fragment.bonds]
    bonds += [(port.atom_id, ids[i], float(o)) for i, o in fragment.attachments]
    bonds = tuple(sorted((min(a, b), max(a, b), o) for a, b, o in bonds))
    new_port = Port(port.port_id, ids[fragment.new_port],
                    tuple(ids[i] for i in fragment.new_caps), port.tag, port.physicality)
    step = ExtensionStep(port_key, depth, provider.name, port, new_port, caps, ids)
    deleted = set(caps)
    return ExtensionMap(
        steps=(step,),
        surviving=tuple((a, row[0], row[1]) for a, row in sorted(reactant_rows.items())
                        if a not in deleted),
        deleted=tuple(sorted(caps)),
        fresh_atoms=tuple((ids[i], tuple(atom)) for i, atom in enumerate(fragment.atoms)),
        fresh_bonds=bonds,
        fresh_coordinates=tuple(sorted((ids[i], port_key[0], u, k)
                                       for i, u, k in fragment.coordinates)),
        ports=((port_key, new_port),),
    )


def _apply_side(side, deleted, fresh_atoms, fresh_bonds, where, rows=None, side_bonds=None):
    rows = _row_map(side) if rows is None else rows
    side_bonds = _bonds(rows) if side_bonds is None else side_bonds
    missing = sorted(set(deleted) - set(rows))
    if missing:
        raise MapCheckError(f"{where}: atoms to delete {missing} are absent")
    atoms = {a: _attrs(row) for a, row in rows.items() if a not in deleted}
    bonds = {pair: o for pair, o in side_bonds.items()
             if pair[0] not in deleted and pair[1] not in deleted}
    for atom_id, attrs in fresh_atoms:
        if atom_id in atoms or atom_id in rows:
            raise MapCheckError(f"{where}: fresh atom id {atom_id} is already in use")
        atoms[atom_id] = tuple(attrs)
    for a, b, order in fresh_bonds:
        if a not in atoms or b not in atoms:
            raise MapCheckError(f"{where}: fresh bond {a}-{b} has an absent endpoint")
        pair = (min(a, b), max(a, b))
        if pair in bonds:
            raise MapCheckError(f"{where}: fresh bond {a}-{b} already exists")
        bonds[pair] = float(order)
    return _side_from(atoms, bonds)


def apply_map(witness: Witness, record: ExtensionMap, *, support=None) -> Witness:
    """Apply an extension map to both sides of a witness identically."""
    if support is None:
        support = edit_support(witness)
    _forbid_support_deletion(support, record.deleted, "apply")
    fresh = [a for a, _attrs in record.fresh_atoms]
    if len(set(fresh)) != len(fresh):
        raise MapCheckError("apply: a fresh atom id is duplicated")
    reactants = _apply_side(witness.reactants, set(record.deleted), record.fresh_atoms,
                            record.fresh_bonds, "apply reactants", witness.reactant_rows(),
                            witness.reactant_bonds())
    products = _apply_side(witness.products, set(record.deleted), record.fresh_atoms,
                           record.fresh_bonds, "apply products", witness.product_rows(),
                           witness.product_bonds())
    ports = dict(record.ports)
    coordinates = {}
    for atom_id, participant, unit, position in record.fresh_coordinates:
        coordinates.setdefault(participant, []).append((atom_id, unit, position))
    deleted = set(record.deleted)
    members = {}
    for step in record.steps:
        members.setdefault(step.port_key[0], []).extend(step.fresh_ids)
    fresh_set = set(fresh)
    members = {k: [a for a in v if a in fresh_set] for k, v in members.items()}
    participants = []
    for participant in witness.participants:
        metadata = participant.metadata
        new_ports = tuple(ports.get((participant.index, p.port_id), p) for p in metadata.ports)
        new_coords = tuple(c for c in metadata.coordinates if c[0] not in deleted) + tuple(
            coordinates.get(participant.index, ()))
        new_atoms = tuple(a for a in participant.atom_ids if a not in deleted) + tuple(
            members.get(participant.index, ()))
        participants.append(Participant(participant.index, metadata.carried(
            new_ports, new_coords, deleted), new_atoms))
    counts = dict(witness.extension_counts)
    for step in record.steps:
        counts[step.port_key] = counts.get(step.port_key, 0) + 1
    return Witness(reactants, products, tuple(participants), witness.id_base,
                   tuple(sorted(counts.items())), witness.label)


def compose(first: ExtensionMap, second: ExtensionMap) -> ExtensionMap:
    """The single map equal to applying ``first`` then ``second``."""
    fresh_first = {a for a, _attrs in first.fresh_atoms}
    removed = set(second.deleted)
    deleted = set(first.deleted) | (removed - fresh_first)
    fresh_atoms = [item for item in first.fresh_atoms if item[0] not in removed]
    fresh_atoms += list(second.fresh_atoms)
    bonds = [b for b in first.fresh_bonds if b[0] not in removed and b[1] not in removed]
    bonds += list(second.fresh_bonds)
    coordinates = [c for c in first.fresh_coordinates if c[0] not in removed]
    coordinates += list(second.fresh_coordinates)
    ports = dict(first.ports)
    ports.update(dict(second.ports))
    return ExtensionMap(
        steps=tuple(sorted(first.steps + second.steps, key=ExtensionStep.sort_key)),
        surviving=tuple(s for s in first.surviving if s[0] not in removed),
        deleted=tuple(sorted(deleted)),
        fresh_atoms=tuple(sorted(fresh_atoms)),
        fresh_bonds=tuple(sorted(bonds)),
        fresh_coordinates=tuple(sorted(coordinates)),
        ports=tuple(sorted(ports.items())),
    )


def identity_map(witness: Witness) -> ExtensionMap:
    """The empty extension of ``witness`` (every atom survives)."""
    rows = witness.reactant_rows()
    return ExtensionMap((), tuple((a, row[0], row[1]) for a, row in sorted(rows.items())),
                        (), (), (), (), ())


# ---------------------------------------------------------------------------
# verification
# ---------------------------------------------------------------------------
def verify_extension(original: Witness, record: ExtensionMap, extended: Witness,
                     providers: Optional[Mapping[str, BoundaryProvider]] = None) -> None:
    """Re-check an extension map from scratch; raise ``MapCheckError``."""
    providers = PROVIDERS if providers is None else providers
    support = edit_support(original)
    _forbid_support_deletion(support, record.deleted, "verify")
    _verify_history(original, record)
    _verify_bookkeeping(original, record, extended)
    deleted = set(record.deleted)
    fresh = [a for a, _attrs in record.fresh_atoms]
    if len(set(fresh)) != len(fresh):
        raise MapCheckError("a fresh atom id is duplicated in the map")
    fresh_set = set(fresh)
    fresh_attrs = dict(record.fresh_atoms)
    surviving = {a: (element, isotope) for a, element, isotope in record.surviving}
    if len(surviving) != len(record.surviving):
        raise MapCheckError("a surviving atom is listed twice")
    attach_ports = {step.old_port.atom_id for step in record.steps}
    record_bonds = {}
    for a, b, order in record.fresh_bonds:
        pair = (min(a, b), max(a, b))
        if pair in record_bonds:
            raise MapCheckError(f"fresh bond {pair} is listed twice")
        if a not in fresh_set and b not in fresh_set:
            raise MapCheckError(f"fresh bond {pair} touches no fresh atom")
        for end in pair:
            if end not in fresh_set and end not in attach_ports:
                raise MapCheckError(f"fresh bond {pair} attaches to {end}, which is not "
                                    "an extended port atom")
        record_bonds[pair] = float(order)
    # full bijection between the extended sides
    ext_before, ext_after = extended.reactant_rows(), extended.product_rows()
    if set(ext_before) != set(ext_after):
        raise MapCheckError("extended sides do not share one atom-id set")
    for atom_id in ext_before:
        if ext_before[atom_id][:2] != ext_after[atom_id][:2]:
            raise MapCheckError(f"atom {atom_id} changes element or isotope across sides")
    views = {"reactants": (original.reactant_rows, original.reactant_bonds,
                           extended.reactant_rows, extended.reactant_bonds),
             "products": (original.product_rows, original.product_bonds,
                          extended.product_rows, extended.product_bonds)}
    for name, (orig_rows, orig_bonds_of, ext_rows, ext_bonds_of) in views.items():
        orig, ext = orig_rows(), ext_rows()
        if not deleted <= set(orig):
            raise MapCheckError(f"{name}: deleted atoms {sorted(deleted - set(orig))} "
                                "were never there")
        if deleted & set(ext):
            raise MapCheckError(f"{name}: deleted atoms {sorted(deleted & set(ext))} "
                                "are still present")
        if fresh_set & set(orig):
            raise MapCheckError(f"{name}: fresh ids {sorted(fresh_set & set(orig))} "
                                "collide with original atoms")
        if set(surviving) != set(orig) - deleted:
            raise MapCheckError(f"{name}: the surviving set is not original minus deleted")
        if set(ext) != set(surviving) | fresh_set:
            raise MapCheckError(f"{name}: extended atoms are not survivors plus fresh atoms "
                                f"(missing {sorted((set(surviving) | fresh_set) - set(ext))}, "
                                f"extra {sorted(set(ext) - set(surviving) - fresh_set)})")
        for atom_id, (element, isotope) in surviving.items():
            if orig[atom_id][:2] != (element, isotope) or ext[atom_id][:2] != (element, isotope):
                raise MapCheckError(f"{name}: surviving atom {atom_id} changed element "
                                    "or isotope")
            if _attrs(orig[atom_id]) != _attrs(ext[atom_id]):
                raise MapCheckError(f"{name}: surviving atom {atom_id} changed attributes")
        for atom_id in fresh_set:
            if _attrs(ext[atom_id]) != tuple(fresh_attrs[atom_id]):
                raise MapCheckError(f"{name}: fresh atom {atom_id} differs from the map")
        orig_bonds, ext_bonds = orig_bonds_of(), ext_bonds_of()
        kept = {p: o for p, o in orig_bonds.items() if p[0] not in deleted and p[1] not in deleted}
        ext_kept = {p: o for p, o in ext_bonds.items()
                    if p[0] not in fresh_set and p[1] not in fresh_set}
        if kept != ext_kept:
            diff = sorted(set(kept.items()) ^ set(ext_kept.items()))
            raise MapCheckError(f"{name}: bonds among surviving atoms changed {diff[:4]}")
        ext_fresh = {p: o for p, o in ext_bonds.items()
                     if p[0] in fresh_set or p[1] in fresh_set}
        if ext_fresh != record_bonds:
            diff = sorted(set(ext_fresh.items()) ^ set(record_bonds.items()))
            raise MapCheckError(f"{name}: bonds of fresh atoms differ from the map {diff[:4]}")
    _verify_steps(original, record, extended, providers)
    _verify_replay(original, record, extended, providers)
    # Implied by the checks above (kept as defence in depth): survivors keep
    # their per-side attributes and bonds, fresh atoms and their bonds equal
    # the map on BOTH sides (so carry no edit), and deleted caps lie outside
    # the edit support (so carried none); hence the edits are transported.
    # The cases it would catch are tested through those checks.
    if edits(extended) != edits(original):
        raise MapCheckError("the edit difference is not transported unchanged")


def _verify_history(original: Witness, record: ExtensionMap) -> None:
    """Every id any step ever introduced (including ids a later step deletes)
    is new: not an original atom, not another step's, never a recycled one."""
    used = set(original.reactant_rows()) | set(original.product_rows())
    for step in sorted(record.steps, key=ExtensionStep.sort_key):
        for atom_id in step.fresh_ids:
            if atom_id in used:
                raise MapCheckError(f"step {step.port_key} depth {step.depth} reuses "
                                    f"historical atom id {atom_id}")
            used.add(atom_id)


def _verify_bookkeeping(original: Witness, record: ExtensionMap, extended: Witness) -> None:
    """Provider, id base, extension counts and label are carried exactly."""
    if [p.index for p in extended.participants] != [p.index for p in original.participants]:
        raise MapCheckError("the participants of the extended witness differ")
    for before, after in zip(original.participants, extended.participants):
        if after.metadata.provider != before.metadata.provider:
            raise MapCheckError(f"participant {before.index}: provider "
                                f"{after.metadata.provider!r} is not the original "
                                f"{before.metadata.provider!r}")
    if extended.id_base != original.id_base:
        raise MapCheckError("id_base is not preserved by the extension")
    expected = dict(original.extension_counts)
    for step in record.steps:
        expected[step.port_key] = expected.get(step.port_key, 0) + 1
    expected = {k: v for k, v in expected.items() if v}
    if {k: v for k, v in extended.extension_counts if v} != expected:
        raise MapCheckError("extension count is not original plus the steps' ports "
                            f"(expected {sorted(expected.items())}, "
                            f"got {list(extended.extension_counts)})")
    if extended.label != original.label:
        raise MapCheckError("the witness label changed")


def _verify_replay(original: Witness, record: ExtensionMap, extended: Witness,
                   providers) -> None:
    """Replay every step from the original witness: each must be exactly what
    ``plan_extension`` produces on the intermediate witness."""
    current, maps = original, []
    for step in sorted(record.steps, key=ExtensionStep.sort_key):
        try:
            provider = providers[step.provider]
            planned = plan_extension(current, step.port_key, provider)
        except (KeyError, BoundaryError) as error:
            raise MapCheckError(f"replay of step {step.port_key} depth {step.depth} "
                                f"failed: {error!r}") from error
        if planned.steps != (step,):
            raise MapCheckError(f"replay: step {step.port_key} depth {step.depth} differs "
                                "from the planned step on the intermediate witness")
        current = apply_map(current, planned)
        maps.append(planned)
    if current != extended:
        raise MapCheckError("replay: the replayed witness differs from the extended witness")
    replayed = identity_map(original)
    for planned in maps:
        replayed = compose(replayed, planned)
    if replayed != record:
        raise MapCheckError("replay: the composed single-step maps differ from the record")


def _verify_steps(original, record, extended, providers) -> None:
    deleted = set(record.deleted)
    fresh_set = set(record.fresh_ids)
    fresh_attrs = dict(record.fresh_atoms)
    record_bonds = {(min(a, b), max(a, b)): o for a, b, o in record.fresh_bonds}
    by_port = {}
    for step in record.steps:
        by_port.setdefault(step.port_key, []).append(step)
    ext_coords = {}
    for participant in extended.participants:
        for a, unit, position in participant.metadata.coordinates:
            ext_coords[a] = (participant.index, unit, position)
    orig_coords = {}
    for participant in original.participants:
        for a, unit, position in participant.metadata.coordinates:
            orig_coords[a] = (participant.index, unit, position)
    recorded_coords = {a: (p, u, k) for a, p, u, k in record.fresh_coordinates}
    final_ports = dict(record.ports)
    claimed, stepped, step_deleted = [], [], set()
    for step in record.steps:
        step_deleted.update(step.deleted)
    for port_key, steps in sorted(by_port.items()):
        steps = sorted(steps, key=ExtensionStep.sort_key)
        current = original.port(port_key)
        depth = original.depth(port_key)
        for step in steps:
            if step.old_port != current:
                raise MapCheckError(f"port {port_key}: step starts from {step.old_port}, "
                                    f"not the declared {current}")
            if step.depth != depth:
                raise MapCheckError(f"port {port_key}: step depth {step.depth} != {depth}")
            if tuple(step.deleted) != tuple(current.cap_atom_ids):
                raise MapCheckError(f"port {port_key}: deleted {step.deleted} are not the "
                                    f"declared caps {current.cap_atom_ids}")
            provider = providers[step.provider]
            participant = original.participant(port_key[0]).metadata.provider
            if step.provider != participant:
                raise MapCheckError(f"port {port_key}: step provider {step.provider!r} "
                                    f"is not the participant's {participant!r}")
            coord = orig_coords.get(current.atom_id) or recorded_coords.get(current.atom_id)
            if coord is None:
                raise MapCheckError(f"port {port_key}: port atom has no repeat coordinate")
            fragment = provider.fragment(current, coord[1:])
            ids = tuple(step.fresh_ids)
            if len(ids) != len(fragment.atoms):
                raise MapCheckError(f"port {port_key}: step size differs from its fragment")
            for i, j, o in fragment.bonds:
                claimed.append(((min(ids[i], ids[j]), max(ids[i], ids[j])), float(o)))
            for i, o in fragment.attachments:
                claimed.append(((min(current.atom_id, ids[i]), max(current.atom_id, ids[i])),
                                float(o)))
            for k, atom_id in enumerate(ids):
                stepped.append(atom_id)
                if atom_id in fresh_set and tuple(fresh_attrs[atom_id]) != tuple(
                        fragment.atoms[k]):
                    raise MapCheckError(f"port {port_key}: fresh atom {atom_id} is not "
                                        "the fragment atom")
            expected_port = Port(current.port_id, ids[fragment.new_port],
                                 tuple(ids[i] for i in fragment.new_caps), current.tag,
                                 current.physicality)
            if step.new_port != expected_port:
                raise MapCheckError(f"port {port_key}: new port {step.new_port} is not the "
                                    f"fragment's {expected_port}")
            for local, unit, position in fragment.coordinates:
                atom_id = ids[local]
                if atom_id in fresh_set and recorded_coords.get(atom_id) != (
                        port_key[0], unit, position):
                    raise MapCheckError(f"port {port_key}: fresh atom {atom_id} has the "
                                        "wrong repeat coordinate")
            current = step.new_port
            depth += 1
        if final_ports.get(port_key) != current:
            raise MapCheckError(f"port {port_key}: the map's final port is not the last step's")
        if extended.port(port_key) != current:
            raise MapCheckError(f"port {port_key}: extended metadata does not carry the new port")
        if extended.depth(port_key) != depth:
            raise MapCheckError(f"port {port_key}: extension count is not carried")
    for port_key in original.port_keys():
        if port_key not in by_port and extended.port(port_key) != original.port(port_key):
            raise MapCheckError(f"port {port_key}: an unextended port changed")
    if set(final_ports) != set(by_port):
        raise MapCheckError("the map lists a port it never extended")
    if len(set(stepped)) != len(stepped):
        raise MapCheckError("two steps share a fresh atom id")
    if fresh_set != set(stepped) - step_deleted:
        raise MapCheckError("fresh atoms are not the steps' atoms minus later-deleted caps")
    if deleted != step_deleted - set(stepped):
        raise MapCheckError("deleted atoms are not the steps' declared original caps")
    expected_bonds = {pair: order for pair, order in claimed
                      if pair[0] not in step_deleted and pair[1] not in step_deleted}
    if len(expected_bonds) != len([c for c in claimed if c[0][0] not in step_deleted
                                   and c[0][1] not in step_deleted]):
        raise MapCheckError("two steps claim one bond")
    if expected_bonds != record_bonds:
        diff = sorted(set(expected_bonds.items()) ^ set(record_bonds.items()))
        raise MapCheckError(f"map bonds are not the steps' fragment bonds {diff[:4]}")
    expected_coords = {a: c for a, c in orig_coords.items() if a not in deleted}
    expected_coords.update(recorded_coords)
    if ext_coords != expected_coords:
        raise MapCheckError("extended repeat coordinates differ from original plus fresh")
    fresh_owner = {}
    for step in record.steps:
        for atom_id in step.fresh_ids:
            fresh_owner[atom_id] = step.port_key[0]
    for participant in extended.participants:
        original_atoms = set(original.participant(participant.index).atom_ids)
        expected = (original_atoms - deleted) | {a for a in fresh_set
                                                 if fresh_owner.get(a) == participant.index}
        if set(participant.atom_ids) != expected:
            raise MapCheckError(f"participant {participant.index}: atom set is not original "
                                "minus deleted plus its fresh atoms")
        participant.metadata.check(extended.reactant_rows(),
                                   f"extended reactants (participant {participant.index})",
                                   participant.atom_ids)


# ---------------------------------------------------------------------------
# extension, vectors, distinctness
# ---------------------------------------------------------------------------
def _ceiling(ceiling):
    if ceiling is None:
        from rmgpy.kmc.compiler import PROXY_PADDING_MAX_HEAVY_ATOMS

        return PROXY_PADDING_MAX_HEAVY_ATOMS
    return int(ceiling)


def _provider_for(witness, port_key, providers):
    providers = PROVIDERS if providers is None else providers
    name = witness.participant(port_key[0]).metadata.provider
    try:
        return providers[name]
    except KeyError as error:
        raise BoundaryError(f"no provider registered as {name!r}") from error


def extend(witness: Witness, port_key, providers=None, *, ceiling=None, verify: bool = True,
           support=None):
    """Extend one artificial port by one repeat; returns (witness, map)."""
    provider = _provider_for(witness, tuple(port_key), providers)
    if support is None:
        support = edit_support(witness)
    record = plan_extension(witness, port_key, provider, support)
    limit = _ceiling(ceiling)
    heavy = witness.heavy_count() + sum(
        1 for _a, attrs in record.fresh_atoms if attrs[0] != "H") - sum(
        1 for a in record.deleted if witness.reactant_rows()[a][0] != "H")
    if heavy > limit:
        raise CeilingExceeded(f"extension reaches {heavy} heavy atoms > ceiling {limit}")
    extended = apply_map(witness, record, support=support)
    if verify:
        verify_extension(witness, record, extended, providers)
    return extended, record


def vector_ports(witness: Witness) -> tuple[tuple[int, str], ...]:
    """The vector coordinates: artificial ports ordered by (participant, port id)."""
    return witness.artificial_port_keys()


def extend_vector(witness: Witness, vector: Sequence[int], providers=None, *,
                  order: Optional[Sequence[int]] = None, ceiling=None, verify: bool = True,
                  verify_steps: bool = False):
    """Apply an extension vector port by port (in ``order``, default
    ascending); returns (witness, composite map).  With ``verify`` the
    composite map is verified against the original witness (which re-checks
    every step against its provider fragment) and must equal the sequential
    application; ``verify_steps`` also verifies each single step."""
    ports = vector_ports(witness)
    if len(vector) != len(ports):
        raise BoundaryError(f"vector has {len(vector)} entries for {len(ports)} ports")
    if any(isinstance(v, bool) or not isinstance(v, int) or v < 0 for v in vector):
        raise BoundaryError("vector entries must be non-negative integers")
    order = tuple(range(len(ports))) if order is None else tuple(order)
    if sorted(order) != list(range(len(ports))):
        raise BoundaryError("order must be a permutation of the port positions")
    # extension never touches the edit support (fresh atoms are identical on
    # both sides), so the original support holds for every intermediate
    support = edit_support(witness)
    current, record = witness, identity_map(witness)
    for position in order:
        for _ in range(vector[position]):
            current, step = extend(current, ports[position], providers, ceiling=ceiling,
                                   verify=verify_steps, support=support)
            record = compose(record, step)
    if verify and any(vector):
        once = apply_map(witness, record, support=support)
        if once != current:
            raise MapCheckError("composite map differs from sequential application")
        verify_extension(witness, record, current, providers)
    return current, record


@dataclass(frozen=True)
class VectorResult:
    vector: tuple[int, ...]
    status: str
    reason: str
    witness: Optional[Witness] = field(default=None, compare=False)
    record: Optional[ExtensionMap] = field(default=None, compare=False)


def _vector_plan(witness: Witness, providers=None) -> tuple[dict, dict]:
    """Per vector position: the refusal reason, or the heavy atoms one step adds."""
    support = edit_support(witness)
    rows = witness.reactant_rows()
    refusals, per_step = {}, {}
    for position, key in enumerate(vector_ports(witness)):
        port = witness.port(key)
        if port.atom_id in support:
            refusals[position] = f"port {key} atom {port.atom_id} is a reacting atom"
            continue
        provider = _provider_for(witness, key, providers)
        coordinate = witness.participant(key[0]).metadata.coordinate(port.atom_id)
        try:
            fragment = provider.fragment(port, coordinate)
        except ExtensionRefused as error:
            refusals[position] = str(error)
            continue
        per_step[position] = fragment.heavy_count() - sum(
            1 for a in port.cap_atom_ids if rows[a][0] != "H")
    return refusals, per_step


def evaluate_vector(witness: Witness, vector: Sequence[int], providers=None, *,
                    ceiling=None, check_order: bool = True, verify: bool = True,
                    plan=None) -> VectorResult:
    """Status of one vector: ``refused`` (a non-zero entry on a port that may
    not be extended), ``ceiling_exceeded`` (the hard heavy-atom ceiling,
    projected before building and enforced per step) or ``ok``.  For an
    ``ok`` vector with several extended ports the reversed port order is
    applied too and must give the same witness and map."""
    vector = tuple(vector)
    ports = vector_ports(witness)
    if len(vector) != len(ports):
        raise BoundaryError(f"vector has {len(vector)} entries for {len(ports)} ports")
    limit = _ceiling(ceiling)
    refusals, per_step = _vector_plan(witness, providers) if plan is None else plan
    refused = [refusals[p] for p, v in enumerate(vector) if v and p in refusals]
    if refused:
        return VectorResult(vector, "refused", "; ".join(refused))
    projected = witness.heavy_count() + sum(v * per_step[p] for p, v in enumerate(vector) if v)
    if projected > limit:
        return VectorResult(vector, "ceiling_exceeded",
                            f"{projected} heavy atoms > ceiling {limit}")
    try:
        extended, record = extend_vector(witness, vector, providers, ceiling=limit,
                                         verify=verify)
    except CeilingExceeded as error:
        return VectorResult(vector, "ceiling_exceeded", str(error))
    if check_order and sum(1 for v in vector if v) > 1:
        other, other_record = extend_vector(
            witness, vector, providers, order=tuple(reversed(range(len(ports)))),
            ceiling=limit, verify=False)
        if other != extended or other_record != record:
            raise MapCheckError(f"vector {vector}: port order changes the extension")
    return VectorResult(vector, "ok", "", extended, record)


def iter_vectors(witness: Witness, box, providers=None, *, ceiling=None,
                 check_order: bool = True, verify: bool = True, budget=None,
                 prune_ceiling: bool = False, stats=None):
    """Generator of ``evaluate_vector`` over the box ``0..box[i]`` (an int
    applies to every port), in lexicographic order, by a lazy odometer
    (nothing is pooled: memory is independent of the box).  ``budget``
    caps every VISITED vector (yielded, whatever its status: refusals and
    ceiling verdicts count as well as built vectors); the vector after the
    budget raises ``BudgetExceeded``.  With ``prune_ceiling`` sub-boxes
    whose projected heavy count already exceeds the ceiling are skipped
    without visiting them (per-step increments are non-negative, so every
    extension of such a prefix exceeds it too; a coordinate's scan stops at
    the first value over the ceiling).  ``stats`` receives the counts
    ``visited``, ``yielded``, ``built`` and ``pruned``."""
    ports = vector_ports(witness)
    bounds = (box,) * len(ports) if isinstance(box, int) else tuple(box)
    if len(bounds) != len(ports) or any(
            isinstance(bound, bool) or not isinstance(bound, int) or bound < 0
            for bound in bounds):
        raise BoundaryError("box must give a non-negative integer bound per port")
    if budget is not None and (isinstance(budget, bool) or not isinstance(budget, int)
                               or budget < 0):
        raise BoundaryError("budget must be a non-negative integer")
    stats = {} if stats is None else stats
    for name in ("visited", "yielded", "built", "pruned"):
        stats.setdefault(name, 0)
    limit = _ceiling(ceiling)
    plan = _vector_plan(witness, providers)
    refusals, per_step = plan
    base = witness.heavy_count()

    def leaf(vector):
        if budget is not None and stats["visited"] >= budget:
            raise BudgetExceeded(f"vector {vector}: more than {budget} vectors to visit")
        stats["visited"] += 1
        refused = any(v and p in refusals for p, v in enumerate(vector))
        projected = base + sum(v * per_step[p] for p, v in enumerate(vector)
                               if v and p in per_step)
        if not refused and projected <= limit:
            stats["built"] += 1
        stats["yielded"] += 1
        return evaluate_vector(witness, vector, providers, ceiling=limit,
                               check_order=check_order, verify=verify, plan=plan)

    n = len(bounds)
    if not prune_ceiling:
        vector = [0] * n
        while True:
            yield leaf(tuple(vector))
            position = n - 1
            while position >= 0 and vector[position] == bounds[position]:
                vector[position] = 0
                position -= 1
            if position < 0:
                return
            vector[position] += 1

    # vectors under a prefix of each length (Python ints: exact for any box)
    below = [1] * (n + 1)
    for p in range(n - 1, -1, -1):
        below[p] = below[p + 1] * (bounds[p] + 1)

    def walk(position, prefix, projected):
        if position == n:
            yield leaf(tuple(prefix))
            return
        step = per_step.get(position, 0)
        for v in range(bounds[position] + 1):
            added = projected + v * step
            if added > limit:  # increments are non-negative: so is every larger v
                stats["pruned"] += (bounds[position] + 1 - v) * below[position + 1]
                break
            yield from walk(position + 1, prefix + [v], added)

    yield from walk(0, [], base)


def enumerate_vectors(witness: Witness, box, providers=None, *, ceiling=None,
                      check_order: bool = True, verify: bool = True, budget=None,
                      prune_ceiling: bool = False) -> tuple:
    """``evaluate_vector`` on every vector of the box (see ``iter_vectors``;
    ``budget`` bounds the visited vectors)."""
    return tuple(iter_vectors(witness, box, providers, ceiling=ceiling,
                              check_order=check_order, verify=verify, budget=budget,
                              prune_ceiling=prune_ceiling))


def _roles(witness: Witness) -> dict:
    return {a: p.index for p in witness.participants for a in p.atom_ids}


def marks(witness: Witness) -> dict:
    """Pinned marks: reacting atoms by atom id, port atoms by (participant,
    port id, tag, physicality, provider), declared cap atoms by (participant,
    port id)."""
    labels = {a: (f"r{a}",) for a in edit_support(witness)}
    for participant in witness.participants:
        metadata = participant.metadata
        for port in metadata.ports:
            labels[port.atom_id] = labels.get(port.atom_id, ()) + (
                f"p{participant.index}:{port.port_id}:{port.tag}:{port.physicality}:"
                f"{metadata.provider}",)
            for cap in port.cap_atom_ids:
                labels[cap] = labels.get(cap, ()) + (f"c{participant.index}:{port.port_id}",)
    return {a: "|".join(sorted(v)) for a, v in labels.items()}


def marked_graph(witness: Witness):
    """networkx graph of the whole witness: node label (element, isotope,
    attributes on both sides, mark), edge label (order before, order after).
    The mark of a node is its participant role plus its pinned marks; the
    node's ``mark`` attribute holds the pinned marks only."""
    import networkx as nx

    before, after = witness.reactant_rows(), witness.product_rows()
    pinned = marks(witness)
    roles = _roles(witness)
    graph = nx.Graph()
    for atom_id in sorted(before):
        pin = pinned.get(atom_id, "")
        mark = "|".join(sorted(filter(None, (f"q{roles.get(atom_id, -1)}", pin))))
        graph.add_node(atom_id, mark=pin, label=repr((_attrs(before[atom_id]),
                                                       _attrs(after[atom_id]), mark)))
    bonds_before, bonds_after = witness.reactant_bonds(), witness.product_bonds()
    for pair in sorted(set(bonds_before) | set(bonds_after)):
        graph.add_edge(*pair, label=repr((bonds_before.get(pair, 0.0),
                                          bonds_after.get(pair, 0.0))))
    return graph


def _invariant(graph) -> tuple:
    labels = tuple(sorted(data["label"] for _n, data in graph.nodes(data=True)))
    marked = sorted((data["label"], n) for n, data in graph.nodes(data=True) if data["mark"])
    distances = []
    for label, node in marked:
        seen = {node: 0}
        frontier = deque([node])
        while frontier:
            current = frontier.popleft()
            for neighbour in graph[current]:
                if neighbour not in seen:
                    seen[neighbour] = seen[current] + 1
                    frontier.append(neighbour)
        distances.append((label, tuple(sorted(
            (graph.nodes[m]["label"], seen.get(m, -1)) for _l, m in marked))))
    return labels, tuple(distances)


def marked_isomorphic(first: Witness, second: Witness) -> bool:
    """Isomorphic by a map that fixes every mark (reacting atoms, ports)."""
    import networkx as nx
    from networkx.algorithms import isomorphism

    g1, g2 = marked_graph(first), marked_graph(second)
    if _invariant(g1) != _invariant(g2):
        return False
    if nx.weisfeiler_lehman_graph_hash(g1, node_attr="label", edge_attr="label") != \
            nx.weisfeiler_lehman_graph_hash(g2, node_attr="label", edge_attr="label"):
        return False
    matcher = isomorphism.GraphMatcher(
        g1, g2, node_match=lambda a, b: a["label"] == b["label"],
        edge_match=lambda a, b: a["label"] == b["label"])
    return matcher.is_isomorphic()


def marked_duplicates(items: Sequence[tuple[Any, Witness]]) -> tuple:
    """Pairs ``(key_a, key_b)`` (first occurrence first) of marked-isomorphic
    witnesses; every item is compared with the first member of each class."""
    classes = []
    duplicates = []
    # only witnesses with equal atom-label multisets can be isomorphic
    composition = {}
    for key, witness in items:
        signature = tuple(sorted(
            (row[:5], row[5] and len(row[5])) for row in witness.reactant_rows().values()))
        composition.setdefault(signature, []).append(key)
    shared = {key for keys in composition.values() if len(keys) > 1 for key in keys}
    for key, witness in items:
        if key not in shared:
            continue
        graph = marked_graph(witness)
        invariant = _invariant(graph)
        for representative_key, representative, rep_invariant in classes:
            if rep_invariant == invariant and marked_isomorphic(representative, witness):
                duplicates.append((representative_key, key))
                break
        else:
            classes.append((key, witness, invariant))
    return tuple(duplicates)


def vector_duplicates(results: Sequence[VectorResult]) -> tuple:
    """Marked duplicates among the ``ok`` vectors of one enumeration."""
    return marked_duplicates([(r.vector, r.witness) for r in results if r.status == "ok"])


# ---------------------------------------------------------------------------
# providers
# ---------------------------------------------------------------------------
def _construction_metadata(provider, objects, roles, ports, name):
    """Metadata from atom objects recorded at construction.

    ``objects``: atoms in construction order (heavy atoms in declaration
    order, then H); ``roles``: heavy declaration index -> (unit, position);
    ``ports``: (port_id, heavy declaration index, tag, physicality)."""
    heavy = [atom for atom in objects if atom.element.number != 1]
    if len(heavy) != len(roles):
        raise MetadataError(f"{name}: {len(heavy)} heavy atoms for {len(roles)} roles")
    coordinates = tuple((atom.id, *roles[k]) for k, atom in enumerate(heavy))
    built = []
    for port_id, index, tag, physicality in ports:
        atom = heavy[index]
        caps = ()
        if physicality == "artificial":
            hydrogens = [h for h in objects if h.element.number == 1 and atom in h.edges]
            if not hydrogens:
                raise MetadataError(f"{name}: artificial port {port_id} has no hydrogen to "
                                    "declare as its cap")
            caps = (hydrogens[-1].id,)  # the last H added at construction
        built.append(Port(port_id, atom.id, caps, tag, physicality))
    return BoundaryMetadata(provider, tuple(built), coordinates)


@dataclass(frozen=True)
class PsChainSpec:
    """A linear head-to-tail PS proxy, declared unit by unit.

    Declared ends (additive, I-137; ``declared``): ``head_end`` / ``tail_end``
    build a PHYSICAL end with that end chemistry (``None``: the artificial,
    H-capped cut point).  In phase the head is a unit position 0 terminal
    (``HEAD_ENDS``) and the tail a position 1 terminal (``TAIL_ENDS``);
    ``head_shift`` starts the chain at position 1 of unit 0 (the head then
    carries a position-1 chemistry, ``TAIL_ENDS``) and ``tail_shift`` ends it
    at position 0 of an extra half unit ``units`` (a ``HEAD_ENDS``
    chemistry): the two phases a backbone cut leaves.  A shifted end is
    always physical.  ``secondary_radical_units`` put a radical on position 0
    of those interior units; ``pendant_loss_units`` drop the phenyl of
    position 1 and leave a radical there (the chain fragment of a pendant
    loss).  ``radical_units`` keeps its meaning (a benzylic radical on
    position 1).  A requested site must exist on the declared backbone
    (refused otherwise, never silently dropped).

    One-carbon chains (both ends on one backbone carbon): ``units=1,
    head_shift=True`` is a lone position-1 carbon (toluene, benzyl radical:
    ``TAIL_ENDS`` at both ends); ``units=0, tail_shift=True`` a lone
    position-0 carbon (methane, methyl radical: ``HEAD_ENDS``).  Two radical
    ends on one carbon build its u2 form (the domain refuses it); an alkene
    end needs two backbone carbons (refused)."""

    units: int
    radical_units: tuple[int, ...] = ()
    radical_end: bool = False
    benzylic_end: bool = False
    ring_radical_unit: Optional[int] = None
    monomer: bool = False
    head_end: Optional[str] = None
    tail_end: Optional[str] = None
    head_shift: bool = False
    tail_shift: bool = False
    secondary_radical_units: tuple[int, ...] = ()
    pendant_loss_units: tuple[int, ...] = ()

    @property
    def declared(self) -> bool:
        """Whether the spec uses the declared-end fields (I-137)."""
        return bool(self.head_end or self.tail_end or self.head_shift or self.tail_shift
                    or self.secondary_radical_units or self.pendant_loss_units)


PS_UNIT_HEAVY = 8  # backbone CH2 (0), backbone CH(Ph) (1), ring carbons (2..7)
# declared PS end chemistries (``domain.END_GROUPS``) by the unit position of
# the terminal backbone carbon: position 0 (CH2-type) and position 1 (the
# phenyl-bearing carbon)
HEAD_ENDS = ("methyl", "primary_radical", "vinylidene")
TAIL_ENDS = ("benzyl", "benzylic_radical", "styryl")
RADICAL_ENDS = ("primary_radical", "benzylic_radical")
ALKENE_ENDS = ("vinylidene", "styryl")
# the end-reconstruction fresh-id slot of a port (the last fragment slot of
# its stride: never reached by an extension under any heavy-atom ceiling the
# stride admits, and checked)
END_SLOT = PORT_STRIDE // FRAGMENT_STRIDE - 1


class EndConflict(ExtensionRefused):
    """A declared end cannot be built without touching the reaction."""


@dataclass(frozen=True)
class EndEdit:
    """One declared end chemistry at a port, in witness ids: ``deleted`` (the
    port's declared caps and, for an alkene end, one H of the partner carbon),
    ``fresh_atoms`` (attributes, in local order) bonded by ``fresh_bonds``
    ((local index, witness atom, order)), ``radicals`` ((atom, radical
    electrons)) and ``bonds`` ((atom, atom, new order)).  Additive (I-159):
    ``fresh_internal_bonds`` ((local, local, order)) bond fresh atoms to each
    other and ``fresh_coordinates`` ((local, unit, position)) give fresh heavy
    atoms their repeat coordinates (a multi-atom end group)."""

    end: str
    deleted: tuple[int, ...]
    fresh_atoms: tuple[tuple, ...] = ()
    fresh_bonds: tuple[tuple[int, int, float], ...] = ()
    radicals: tuple[tuple[int, int], ...] = ()
    bonds: tuple[tuple[int, int, float], ...] = ()
    fresh_internal_bonds: tuple[tuple[int, int, float], ...] = ()
    fresh_coordinates: tuple[tuple[int, int, int], ...] = ()


def _free_hydrogen(rows, atom, protected) -> int:
    hydrogens = sorted(n for n, _o in rows[atom][5] if rows[n][0] == "H")
    free = [h for h in hydrogens if h not in protected]
    if not free:
        raise EndConflict(f"atom {atom} has no hydrogen outside the reaction's edit support")
    return free[0]


def standard_end_edit(port: Port, end: str, partner: Optional[int], rows,
                      protected) -> EndEdit:
    """The edit of one declared end on a CH-capped carbon chain: saturated ends
    (methyl, benzyl) replace the cap by a FRESH end hydrogen; radical ends
    remove the cap and put the unpaired electron on the port atom; alkene ends
    remove the cap and one partner H and make the port-partner bond double."""
    caps = tuple(port.cap_atom_ids)
    if end in ("methyl", "benzyl"):
        return EndEdit(end, caps, (("H", None, 0, 0, 0),), ((0, port.atom_id, 1.0),))
    if end in RADICAL_ENDS:
        return EndEdit(end, caps, radicals=((port.atom_id, 1),))
    if end in ALKENE_ENDS:
        if partner is None:
            raise ExtensionRefused(f"{end}: the port atom has no backbone partner")
        if partner in protected:
            raise EndConflict(f"{end}: the alkene partner {partner} is a reacting atom")
        h = _free_hydrogen(rows, partner, protected)
        return EndEdit(end, caps + (h,), bonds=((port.atom_id, partner, 2.0),))
    raise ExtensionRefused(f"undeclared end chemistry {end!r}")


class PolystyreneProvider(BoundaryProvider):
    """Head-to-tail PS.  Unit u: position 0 the backbone CH2 (head side),
    position 1 the phenyl-bearing CH (tail side), positions 2..7 the ring
    (ipso first, in ring order).  Head port: unit 0 position 0; tail port:
    last unit position 1."""

    name = "polystyrene"
    # declared end chemistries per port tag (additive, I-137)
    END_CHEMISTRIES = {"backbone_head": HEAD_ENDS, "backbone_tail": TAIL_ENDS}
    RADICAL_END_CHEMISTRIES = RADICAL_ENDS
    ALKENE_END_CHEMISTRIES = ALKENE_ENDS
    coordinate_backbone_atoms = 2      # backbone CH2 (position 0) and CH(Ph) (1)

    def backbone_index(self, unit: int, position: int) -> Optional[int]:
        return 2 * unit + position if position in (0, 1) else None

    def interior_atoms_rmg(self) -> list:
        """RMG atoms (atom types perceived) of unit 1 of a pristine DP3 and
        their H."""
        build = self.build("interior", PsChainSpec(3))
        coords = {a: (u, k) for a, u, k in build.metadata.coordinates}
        heavy = {a for a, (u, _k) in coords.items() if u == 1}
        atoms = [a for a in build.molecule.atoms if a.id in heavy]
        return atoms + [h for a in atoms for h in a.bonds if h.is_hydrogen()]

    def probe_ports(self) -> tuple:
        return LEGACY_PROBE_PORTS

    def unit_formula(self) -> tuple:
        return (("C", 8), ("H", 8))

    @staticmethod
    def normalize_spec(spec):
        """The spec if it is this provider's (a ``PsChainSpec``), else None."""
        return spec if isinstance(spec, PsChainSpec) else None

    # The repeat unit this provider builds, in the ``polymer_spec`` graph
    # declaration schema (heavy atoms with their hydrogen counts, aromatic
    # bonds 1.5, head/tail ports): a roster matches a declared repeat unit to
    # a provider by graph isomorphism against this (additive, I-155).
    REPEAT_UNIT = {
        "atoms": [{"key": "c0", "element": "C", "hydrogens": 2},
                  {"key": "c1", "element": "C", "hydrogens": 1},
                  {"key": "r0", "element": "C", "hydrogens": 0}] +
                 [{"key": f"r{k}", "element": "C", "hydrogens": 1} for k in range(1, 6)],
        "bonds": [{"atoms": ["c0", "c1"], "order": 1}, {"atoms": ["c1", "r0"], "order": 1}] +
                 [{"atoms": [f"r{k}", f"r{(k + 1) % 6}"], "order": 1.5} for k in range(6)],
        "ports": {"head": {"atom": "c0", "bond_order": 1}, "tail": {"atom": "c1", "bond_order": 1}},
        "electronic_state": {"multiplicity": 1},
    }

    def spec_from_declaration(self, dp, head=None, tail=None, interior_radical=None):
        """The chain spec of ``dp`` declared repeat units with the head/tail
        end variants named from ``END_CHEMISTRIES`` (``None``: the artificial
        H-capped end) and, optionally, a radical on the phenyl-bearing carbon
        of unit ``interior_radical``.  ValueError on a name or count this
        provider does not declare."""
        if not isinstance(dp, int) or dp < 0:
            raise ValueError(f"degree of polymerization must be a non-negative integer, got {dp!r}")
        for end, allowed in ((head, HEAD_ENDS + TAIL_ENDS), (tail, HEAD_ENDS + TAIL_ENDS)):
            if end is not None and end not in allowed:
                raise ValueError(f"{end!r} is not a declared end (declared: {allowed})")
        extra = {} if interior_radical is None else {"radical_units": (int(interior_radical),)}
        return PsChainSpec(dp, head_end=head, tail_end=tail, head_shift=head in TAIL_ENDS,
                           tail_shift=tail in HEAD_ENDS, **extra)

    def spec_for_declared_provider(self, declared, dp, head=None, tail=None,
                                   interior_radical=None):
        """Translate a generated PS reference shape to the declaration's
        primitive-period ``ChainSpec``; the resulting build is owned entirely
        by ``declared`` even when this provider supplied the phase vocabulary."""
        from rmgpy.kmc.declared_chain import ChainSpec, EditOp
        from rmgpy.kmc.windows import DeclaredSpec

        generic_head = {"declared": "methyl", "radical": "primary_radical"}
        generic_tail = {"declared": "benzyl", "radical": "benzylic_radical"}
        head = generic_head.get(head, head)
        tail = generic_tail.get(tail, tail)
        spec = self.spec_from_declaration(dp, head, tail, interior_radical)
        if spec.ring_radical_unit is not None or spec.secondary_radical_units or \
                spec.pendant_loss_units or spec.monomer:
            raise ExtensionRefused(f"unsupported generated PS reference shape {spec!r}")

        ports = dict(declared.spec.repeat_unit.ports)
        keys = (ports["head"].atom, ports["tail"].atom)
        phase = next((item.index for item in declared.geometry.phases
                      if set(item.bond) == set(keys)), None)
        if (spec.head_shift or spec.tail_shift) and phase is None:
            raise ExtensionRefused("the declaration has no phase between its port atoms")

        units = spec.units + int(spec.tail_shift)
        head_kind = "radical" if spec.radical_end or spec.head_end in RADICAL_ENDS else \
            "declared" if spec.head_end else "artificial"
        tail_kind = "radical" if spec.benzylic_end or spec.tail_end in RADICAL_ENDS else \
            "declared" if spec.tail_end else "artificial"
        edits = []
        for unit in sorted(set(spec.radical_units)):
            edits += [EditOp("remove_h", unit, keys[1]),
                      EditOp("set_radical", unit, keys[1], 1)]

        backbone = []
        for unit in range(spec.units):
            if not (unit == 0 and spec.head_shift):
                backbone.append((unit, 0))
            backbone.append((unit, 1))
        if spec.tail_shift:
            backbone.append((spec.units, 0))
        doubles = []
        for end, pair in ((spec.head_end, backbone[:2]),
                          (spec.tail_end, backbone[-2:])):
            if end in ALKENE_ENDS:
                if len(pair) < 2:
                    raise ExtensionRefused(f"{end}: one backbone carbon")
                pair = tuple(sorted(pair))
                if pair not in doubles:
                    doubles.append(pair)
        if len(doubles) == 2 and set(doubles[0]) & set(doubles[1]):
            raise ExtensionRefused("two alkene ends share a backbone carbon")
        for (u1, k1), (u2, k2) in doubles:
            edits += [EditOp("remove_h", u1, keys[k1]),
                      EditOp("remove_h", u2, keys[k2]),
                      EditOp("set_bond_order", u1, keys[k1], 2,
                             other=(u2, keys[k2]))]
        return DeclaredSpec(declared.name, ChainSpec(
            units, head=head_kind, tail=tail_kind,
            head_phase=phase if spec.head_shift else 0,
            tail_phase=phase if spec.tail_shift else 0, edits=tuple(edits)))

    def reference_special_declarations(self) -> tuple:
        return (("methane", 0, "methyl", "methyl", None),
                ("methyl-radical", 0, "primary_radical", "methyl", None),
                ("toluene", 1, "benzyl", "benzyl", None),
                ("benzyl-radical", 1, "benzyl", "benzylic_radical", None))

    @staticmethod
    def declared_layout(spec: PsChainSpec):
        """(SMILES, heavy roles in SMILES order, ports) of a declared-end spec.
        Each backbone carbon is followed by its ring (ipso first)."""
        if spec.monomer or spec.radical_end or spec.benzylic_end or \
                spec.ring_radical_unit is not None:
            raise MetadataError("declared ends do not combine with the legacy end flags")
        if spec.units < 1 and not (spec.units == 0 and spec.tail_shift and
                                   not spec.head_shift):
            raise MetadataError("a PS chain needs at least one unit (units=0 only as the "
                                "lone position-0 carbon of a tail-shifted spec)")
        if spec.head_shift and spec.head_end not in TAIL_ENDS:
            raise MetadataError(f"a shifted head needs one of {TAIL_ENDS}, got {spec.head_end!r}")
        if not spec.head_shift and spec.head_end not in (None,) + HEAD_ENDS:
            raise MetadataError(f"an in-phase head needs one of {HEAD_ENDS}, got "
                                f"{spec.head_end!r}")
        if spec.tail_shift and spec.tail_end not in HEAD_ENDS:
            raise MetadataError(f"a shifted tail needs one of {HEAD_ENDS}, got {spec.tail_end!r}")
        if not spec.tail_shift and spec.tail_end not in (None,) + TAIL_ENDS:
            raise MetadataError(f"an in-phase tail needs one of {TAIL_ENDS}, got "
                                f"{spec.tail_end!r}")
        backbone = []
        for u in range(spec.units):
            if not (u == 0 and spec.head_shift):
                backbone.append((u, 0))
            backbone.append((u, 1))
        if spec.tail_shift:
            backbone.append((spec.units, 0))
        last = len(backbone) - 1
        if last < 0:
            raise MetadataError("a declared chain needs a backbone carbon")
        if last == 0 and (spec.head_end in ALKENE_ENDS or spec.tail_end in ALKENE_ENDS):
            raise MetadataError("an alkene end needs two backbone carbons")
        sites = set(backbone)
        for flag, position in ((spec.radical_units, 1), (spec.pendant_loss_units, 1),
                               (spec.secondary_radical_units, 0)):
            missing = sorted((u, position) for u in flag if (u, position) not in sites)
            if missing:
                raise MetadataError(f"requested sites {missing} are not on the declared "
                                    "backbone")
        double = set()   # i: the bond i-(i+1) is double
        if spec.head_end in ALKENE_ENDS:
            double.add(0)
        if spec.tail_end in ALKENE_ENDS:
            double.add(last - 1)
        # two backbone carbons with both ends unsaturated: one double bond read
        # from either end (styrene is vinylidene at one end, styryl at the other)
        same_bond = last == 1 and len(double) == 1
        if len(double) == 2 and max(double) - min(double) < 2 or \
                spec.head_end in ALKENE_ENDS and spec.tail_end in ALKENE_ENDS and last < 3 and \
                not same_bond:
            raise MetadataError("two alkene ends share a backbone carbon")
        interior = set(range(1, last))
        ring, radical = {}, {}
        for i, (u, k) in enumerate(backbone):
            ring[i] = k == 1 and u not in spec.pendant_loss_units
            flagged = (k == 1 and (u in spec.radical_units or u in spec.pendant_loss_units)) or \
                (k == 0 and u in spec.secondary_radical_units)
            if flagged and i not in interior:
                raise MetadataError(f"a mid-chain radical flag names the terminal {(u, k)}")
            radical[i] = int(flagged)
        for flag in (spec.radical_units, spec.pendant_loss_units, spec.secondary_radical_units):
            if any(not 0 <= u < spec.units + 1 for u in flag):
                raise MetadataError("a unit flag is outside the chain")
        if spec.head_end in RADICAL_ENDS:
            radical[0] += 1
        if spec.tail_end in RADICAL_ENDS:
            radical[last] += 1
        tokens, roles = [], []
        for i, (u, k) in enumerate(backbone):
            heavy = (i > 0) + (i < last) + ring[i]
            pi = (i in double) + (i - 1 in double)
            h = 4 - heavy - pi - radical[i]
            if h < 0:
                raise MetadataError(f"backbone carbon {(u, k)} has no valence left")
            token = ("[C]" if h == 0 else "[CH]" if h == 1 else f"[CH{h}]") if radical[i] else "C"
            if ring[i]:
                token += "(c1ccccc1)"
            if i in double:
                token += "="
            tokens.append(token)
            roles.append((u, k))
            if ring[i]:
                roles.extend((u, 2 + r) for r in range(6))
        tail_index = len(roles) - (7 if ring[last] else 1)
        ports = [("head", 0, "backbone_head", "physical" if spec.head_end else "artificial"),
                 ("tail", tail_index, "backbone_tail",
                  "physical" if spec.tail_end else "artificial")]
        return "".join(tokens), tuple(roles), tuple(ports)

    @staticmethod
    def backbone_partner(metadata: BoundaryMetadata, atom_id: int) -> Optional[int]:
        """The backbone carbon bonded to ``atom_id`` in its own repeat unit
        (positions 0 and 1 of one unit), or None."""
        unit, position = metadata.coordinate(atom_id)
        if position not in (0, 1):
            return None
        for a, u, k in metadata.coordinates:
            if u == unit and k == 1 - position:
                return a
        return None

    def end_edit(self, metadata: BoundaryMetadata, port: Port, end: str, rows,
                 protected) -> EndEdit:
        """The declared end chemistry ``end`` built at ``port`` (see
        ``standard_end_edit``); only the chemistries of the port's tag in its
        unit phase (``END_CHEMISTRIES``)."""
        _unit, position = metadata.coordinate(port.atom_id)
        allowed = HEAD_ENDS if position == 0 else TAIL_ENDS if position == 1 else ()
        if end not in allowed:
            raise ExtensionRefused(f"polystyrene: {end!r} is not a declared end of a unit "
                                   f"position {position} terminal (declared: {allowed})")
        partner = self.backbone_partner(metadata, port.atom_id)
        if end in ALKENE_ENDS and partner is not None:
            bonded = {n for n, _o in rows[port.atom_id][5]}
            if partner not in bonded:
                partner = None
        return standard_end_edit(port, end, partner, rows, protected)

    @staticmethod
    def smiles(spec: PsChainSpec) -> str:
        if getattr(spec, "declared", False):
            return PolystyreneProvider.declared_layout(spec)[0]
        if spec.monomer:
            return "C=Cc1ccccc1"
        if spec.units < 1:
            raise MetadataError("a PS chain needs at least one unit")
        pieces = []
        for u in range(spec.units):
            pieces.append("[CH2]" if spec.radical_end and u == 0 else "C")
            ring = "c1cc[c]cc1" if spec.ring_radical_unit == u else "c1ccccc1"
            if spec.benzylic_end and u == spec.units - 1:
                pieces.append(f"[CH]({ring})")
            elif u in spec.radical_units:
                pieces.append(f"[C]({ring})")
            else:
                pieces.append(f"C({ring})")
        return "".join(pieces)

    def build(self, name: str, spec: PsChainSpec) -> ProxyBuild:
        """Build exactly as ``ps_proxy_set`` always has (SMILES, ``update``,
        ``assign_atom_ids``) while recording each atom's role from the
        declaration: RMG keeps SMILES heavy-atom order, then the hydrogens."""
        from rmgpy.molecule.molecule import Molecule

        molecule = Molecule(smiles=self.smiles(spec))
        objects = list(molecule.atoms)
        molecule.update()
        molecule.assign_atom_ids()
        units = 1 if spec.monomer else spec.units
        roles = [(u, k) for u in range(units) for k in range(PS_UNIT_HEAVY)]
        ports = []
        if spec.declared:
            _smiles, roles, ports = self.declared_layout(spec)
        elif not spec.monomer:
            ports.append(("head", 0, "backbone_head",
                          "physical" if spec.radical_end else "artificial"))
            ports.append(("tail", (units - 1) * PS_UNIT_HEAVY + 1, "backbone_tail",
                          "physical" if spec.benzylic_end else "artificial"))
        metadata = _construction_metadata(self.name, objects, roles, ports, name)
        build = ProxyBuild(name, molecule, metadata, tuple(atom.id for atom in objects))
        metadata.check(dict(rows_of_molecule(molecule)), name)
        return build

    def fragment(self, port: Port, coordinate: tuple[int, int]) -> Fragment:
        unit, position = coordinate
        carbon = ("C", None, 0, 0, 0)
        hydrogen = ("H", None, 0, 0, 0)
        ring = [(2 + k, 2 + (k + 1) % 6, 1.5) for k in range(6)]
        if port.tag == "backbone_tail" and position == 1:
            # 0 CH2, 1 CH(Ph)+cap, 2..7 ring; H: 8, 9 on 0; 10 on 1; 11 cap on 1; 12..16 ring
            atoms = (carbon,) * 8 + (hydrogen,) * 9
            bonds = [(0, 1, 1.0), (1, 2, 1.0)] + ring + [
                (0, 8, 1.0), (0, 9, 1.0), (1, 10, 1.0), (1, 11, 1.0)] + [
                (3 + k, 12 + k, 1.0) for k in range(5)]
            return Fragment(atoms, tuple(bonds), ((0, 1.0),), 1, (11,),
                            tuple((k, unit + 1, k) for k in range(8)))
        if port.tag == "backbone_head" and position == 0:
            # 0 CH2+cap, 1 CH(Ph), 2..7 ring; H: 8, 9 on 0; 10 cap on 0; 11 on 1; 12..16 ring
            atoms = (carbon,) * 8 + (hydrogen,) * 9
            bonds = [(0, 1, 1.0), (1, 2, 1.0)] + ring + [
                (0, 8, 1.0), (0, 9, 1.0), (0, 10, 1.0), (1, 11, 1.0)] + [
                (3 + k, 12 + k, 1.0) for k in range(5)]
            return Fragment(atoms, tuple(bonds), ((1, 1.0),), 0, (10,),
                            tuple((k, unit - 1, k) for k in range(8)))
        raise ExtensionRefused(f"polystyrene has no extension path for a {port.tag} port "
                               f"at unit position {position}")


# The PS proxy structures of ``compiler.ps_proxy_set`` (name -> spec builder).
def ps_structure_specs(units: int) -> dict:
    return {
        "pristine": PsChainSpec(units),
        "interior_radical": PsChainSpec(units, radical_units=(units // 2,)),
        "doubly_featured": PsChainSpec(units, radical_units=(0, 1)),
        "end_radical": PsChainSpec(units, radical_end=True),
        "end_radical_short": PsChainSpec(units - 1, radical_end=True),
        "benzylic_end_radical": PsChainSpec(units, benzylic_end=True),
        "benzylic_end_radical_short": PsChainSpec(units - 1, benzylic_end=True),
        "junction_radical": PsChainSpec(units, ring_radical_unit=0),
        "styrene": PsChainSpec(1, monomer=True),
    }


class LinearAlkaneProvider(BoundaryProvider):
    """Unbranched CnH2n+2 as an explicit graph, one-carbon CH2 repeat, methyl
    ends.  Ids: carbons 0..n-1 along the chain, then each carbon's H in
    carbon order (the last H of each end carbon is its declared cap).  Head
    port carbon 0, tail port carbon n-1; carbon i has coordinate (i, 0)."""

    name = "linear_alkane"
    # declared end chemistries (additive, I-137): the saturated methyl end and
    # the primary end radical (both in the declared domain for an unbranched
    # chain)
    END_CHEMISTRIES = {"backbone_head": ("methyl", "primary_radical"),
                       "backbone_tail": ("methyl", "primary_radical")}
    RADICAL_END_CHEMISTRIES = ("primary_radical",)
    coordinate_backbone_atoms = 1      # the one-carbon repeat (``fragment``)

    def backbone_index(self, unit: int, position: int) -> Optional[int]:
        return unit

    def interior_atoms_rmg(self) -> list:
        """RMG atoms (atom types perceived) of carbon 2 of a C5 and its H."""
        build = self.build(5)
        atoms = [a for a in build.molecule.atoms if a.id in {2}]
        return atoms + [h for a in atoms for h in a.bonds if h.is_hydrogen()]

    def probe_ports(self) -> tuple:
        return LEGACY_PROBE_PORTS

    def unit_formula(self) -> tuple:
        return (("C", 1), ("H", 2))

    @staticmethod
    def normalize_spec(spec):
        """("alkane", n, head end, tail end, interior radical) of an
        ("alkane", n[, head[, tail[, interior radical]]]) spec; None for
        another shape; ValueError for a malformed one."""
        if not isinstance(spec, tuple) or not spec or spec[0] != "alkane":
            return None
        if not 2 <= len(spec) <= 5:
            raise ValueError(f"an alkane spec has 2 to 5 entries, got {len(spec)}: {spec!r}")
        if isinstance(spec[1], bool) or not isinstance(spec[1], int):
            raise ValueError(f"an alkane spec needs an integer carbon count: {spec!r}")
        return ("alkane", spec[1]) + tuple(spec[2:]) + (None,) * (5 - len(spec))

    # the two-carbon repeat of the declared chain (additive, I-155; see
    # ``PolystyreneProvider.REPEAT_UNIT``)
    REPEAT_UNIT = {
        "atoms": [{"key": "c0", "element": "C", "hydrogens": 2},
                  {"key": "c1", "element": "C", "hydrogens": 2}],
        "bonds": [{"atoms": ["c0", "c1"], "order": 1}],
        "ports": {"head": {"atom": "c0", "bond_order": 1}, "tail": {"atom": "c1", "bond_order": 1}},
        "electronic_state": {"multiplicity": 1},
    }
    REPEAT_CARBONS = 2

    def spec_from_declaration(self, dp, head=None, tail=None, interior_radical=None):
        """The chain spec of ``dp`` declared repeat units (``2 * dp`` carbons)
        with end variants from ``END_CHEMISTRIES``.  ValueError on a name, count
        or interior radical this provider does not declare."""
        if not isinstance(dp, int) or dp < 1:
            raise ValueError(f"degree of polymerization must be a positive integer, got {dp!r}")
        carbons = self.REPEAT_CARBONS * dp
        for end in (head, tail):
            if end is not None and end not in self.END_CHEMISTRIES["backbone_head"]:
                raise ValueError(f"{end!r} is not a declared end "
                                 f"(declared: {self.END_CHEMISTRIES['backbone_head']})")
        radical = None
        if interior_radical is not None:
            if isinstance(interior_radical, bool) or not isinstance(interior_radical, int) or \
                    not 0 <= interior_radical < dp:
                raise ValueError(f"interior radical unit must be in 0..{dp - 1}, got "
                                 f"{interior_radical!r}")
            radical = self.REPEAT_CARBONS * interior_radical
        if head is None and tail is None and radical is None:
            return ("alkane", carbons)
        return ("alkane", carbons, head, tail, radical)

    def end_edit(self, metadata: BoundaryMetadata, port: Port, end: str, rows,
                 protected) -> EndEdit:
        if end not in self.END_CHEMISTRIES.get(port.tag, ()):
            raise ExtensionRefused(f"linear alkane: {end!r} is not a declared end")
        return standard_end_edit(port, end, None, rows, protected)

    def build(self, carbons: int, name: Optional[str] = None, head_end: Optional[str] = None,
              tail_end: Optional[str] = None,
              interior_radical: Optional[int] = None) -> ProxyBuild:
        """``head_end`` / ``tail_end`` (additive, I-137 rework 2): build that end
        PHYSICAL with a declared end chemistry (``methyl``: the CH3 terminus,
        no cap; ``primary_radical``: a .CH2 terminus) instead of the
        artificial H-capped cut point."""
        from rmgpy.kmc.rmg_moves import molecule_from_rows

        if carbons < 2:
            raise MetadataError("a linear alkane proxy needs at least two carbons")
        for end in (head_end, tail_end):
            if end is not None and end not in self.END_CHEMISTRIES["backbone_head"]:
                raise MetadataError(f"linear alkane: {end!r} is not a declared end")
        if interior_radical is not None and (isinstance(interior_radical, bool) or
                                             not isinstance(interior_radical, int) or
                                             not 0 <= interior_radical < carbons):
            raise MetadataError(f"linear alkane: interior radical carbon {interior_radical!r} "
                                f"is outside 0..{carbons - 1}")
        name = name or f"C{carbons}"
        neighbours = {i: [] for i in range(carbons)}
        for i in range(carbons - 1):
            neighbours[i].append(i + 1)
            neighbours[i + 1].append(i)
        radical = {0: int(head_end == "primary_radical"),
                   carbons - 1: int(tail_end == "primary_radical")}
        if interior_radical is not None:
            radical[interior_radical] = radical.get(interior_radical, 0) + 1
        h_of, next_id = {}, carbons
        for i in range(carbons):
            count = 4 - len(neighbours[i]) - radical.get(i, 0)
            h_of[i] = tuple(range(next_id, next_id + count))
            next_id += count
        rows = {}
        for i in range(carbons):
            edges = [(j, 1.0) for j in neighbours[i]] + [(h, 1.0) for h in h_of[i]]
            rows[i] = ("C", None, 0, radical.get(i, 0), 0, tuple(sorted(edges)))
            for h in h_of[i]:
                rows[h] = ("H", None, 0, 0, 0, ((i, 1.0),))
        molecule = molecule_from_rows(rows)

        def port(port_id, atom, tag, end):
            if end is None:
                return Port(port_id, atom, (h_of[atom][-1],), tag, "artificial")
            return Port(port_id, atom, (), tag, "physical")
        ports = (port("head", 0, "backbone_head", head_end),
                 port("tail", carbons - 1, "backbone_tail", tail_end))
        metadata = BoundaryMetadata(self.name, ports, tuple((i, i, 0) for i in range(carbons)))
        metadata.check(rows, name)
        return ProxyBuild(name, molecule, metadata, tuple(range(next_id)))

    def fragment(self, port: Port, coordinate: tuple[int, int]) -> Fragment:
        """One-atom repeat: the attached carbon IS the new port (head and tail
        coincide), so it carries the new cap itself."""
        unit, _position = coordinate
        step = {"backbone_head": -1, "backbone_tail": 1}.get(port.tag)
        if step is None:
            raise ExtensionRefused(f"linear alkane has no extension path for a {port.tag} port")
        atoms = (("C", None, 0, 0, 0),) + (("H", None, 0, 0, 0),) * 3
        bonds = ((0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0))
        return Fragment(atoms, bonds, ((0, 1.0),), 0, (3,), ((0, unit + step, 0),))


PROVIDERS = {
    PolystyreneProvider.name: PolystyreneProvider(),
    LinearAlkaneProvider.name: LinearAlkaneProvider(),
}
# the named providers, in registry order: the classifier and fragment-fact
# providers of a witness without a declared participant (I-165), and the
# owners of the named spec shapes (``normalize_spec``)
LEGACY_PROVIDERS = (PROVIDERS[PolystyreneProvider.name], PROVIDERS[LinearAlkaneProvider.name])


# ---------------------------------------------------------------------------
# physical ends (additive, I-137)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class EndMap:
    """Exact record of one end reconstruction: the port made physical, the
    deleted atoms (declared caps plus any H the end chemistry removes), the
    fresh atoms with their bonds, the radical increments and the bond orders
    set (identically on both reaction sides)."""

    port_key: tuple[int, str]
    end: str
    old_port: Port
    new_port: Port
    deleted: tuple[int, ...]
    fresh_atoms: tuple[tuple[int, tuple], ...]
    fresh_bonds: tuple[tuple[int, int, float], ...]
    radicals: tuple[tuple[int, int], ...]
    bonds: tuple[tuple[int, int, float], ...]
    fresh_coordinates: tuple[tuple[int, int, int], ...] = ()

    @property
    def fresh_ids(self) -> tuple[int, ...]:
        return tuple(a for a, _attrs in self.fresh_atoms)


def end_chemistries(witness: Witness, port_key, providers=None) -> tuple:
    """The declared end chemistries a port can be reconstructed with."""
    provider = _provider_for(witness, tuple(port_key), providers)
    return tuple(getattr(provider, "END_CHEMISTRIES", {}).get(witness.port(port_key).tag, ()))


def reconstruct_end(witness: Witness, port_key, end: str, providers=None):
    """Make an artificial port PHYSICAL by building the declared end
    chemistry ``end`` with the participant's provider (``end_edit``): its
    declared caps are DELETED, fresh atoms (e.g. the end H of a methyl) get
    new ids from the port's end slot, radicals and bond orders are set, all
    identically on both reaction sides.  Never a relabelled cap.  Refused
    (``EndConflict``) when an atom to delete or a bond to set belongs to the
    reaction's edit support.  Returns (witness, ``EndMap``)."""
    port_key = tuple(port_key)
    port = witness.port(port_key)
    if not port.artificial:
        raise ExtensionRefused(f"port {port_key} is already physical")
    provider = _provider_for(witness, port_key, providers)
    if not hasattr(provider, "end_edit"):
        raise ExtensionRefused(f"provider {provider.name!r} declares no end chemistry")
    metadata = witness.participant(port_key[0]).metadata
    support = edit_support(witness)
    rows = witness.reactant_rows()
    edit = provider.end_edit(metadata, port, end, rows, support)
    hit = sorted(set(edit.deleted) & support)
    if hit:
        raise EndConflict(f"{end} at {port_key}: deleting reacting atoms {hit}")
    hit = sorted({a for a, _r in edit.radicals} & support)
    if hit:
        raise EndConflict(f"{end} at {port_key}: the end radical would sit on reacting atoms "
                          f"{hit} (a different event: build that end on the reactants)")
    before, after = witness.reactant_bonds(), witness.product_bonds()
    for a, b, _order in edit.bonds:
        pair = (min(a, b), max(a, b))
        if before.get(pair) != after.get(pair):
            raise EndConflict(f"{end} at {port_key}: the bond {pair} is edited by the reaction")
    rank = witness.port_keys().index(port_key)
    start = witness.id_base + rank * PORT_STRIDE + END_SLOT * FRAGMENT_STRIDE
    ids = tuple(range(start, start + len(edit.fresh_atoms)))
    if len(ids) > FRAGMENT_STRIDE:
        raise BoundaryError("end fragment is larger than the fresh-id stride")
    fresh_atoms = tuple((ids[k], tuple(attrs)) for k, attrs in enumerate(edit.fresh_atoms))
    fresh_bonds = tuple(sorted([(min(ids[k], a), max(ids[k], a), float(o))
                                for k, a, o in edit.fresh_bonds] +
                               [(min(ids[i], ids[j]), max(ids[i], ids[j]), float(o))
                                for i, j, o in edit.fresh_internal_bonds]))
    fresh_coordinates = tuple(sorted((ids[k], u, q) for k, u, q in edit.fresh_coordinates))
    deleted = set(edit.deleted)
    radicals = dict(edit.radicals)
    orders = {(min(a, b), max(a, b)): float(o) for a, b, o in edit.bonds}

    def side(rows_, bonds_):
        atoms = {}
        for a, row in rows_.items():
            if a in deleted:
                continue
            attrs = list(row[:5])
            attrs[3] = attrs[3] + radicals.get(a, 0)
            atoms[a] = tuple(attrs)
        bonds = {p: o for p, o in bonds_.items() if p[0] not in deleted and p[1] not in deleted}
        for pair, order in orders.items():
            if pair not in bonds:
                raise MapCheckError(f"end {end}: bond {pair} to set is absent")
            bonds[pair] = order
        for a, attrs in fresh_atoms:
            atoms[a] = attrs
        for a, b, order in fresh_bonds:
            bonds[(a, b)] = order
        return _side_from(atoms, bonds)

    reactants = side(witness.reactant_rows(), before)
    products = side(witness.product_rows(), after)
    new_port = replace(port, cap_atom_ids=(), physicality="physical")
    participants = []
    for participant in witness.participants:
        if participant.index != port_key[0]:
            participants.append(participant)
            continue
        meta = participant.metadata
        ports = tuple(new_port if p.port_id == port.port_id else p for p in meta.ports)
        atoms = tuple(a for a in participant.atom_ids if a not in deleted) + ids
        participants.append(Participant(participant.index,
                                        meta.carried(ports, meta.coordinates + fresh_coordinates,
                                                     set(deleted)),
                                        atoms))
    extended = Witness(reactants, products, tuple(participants), witness.id_base,
                       witness.extension_counts, witness.label)
    record = EndMap(port_key, end, port, new_port, tuple(sorted(deleted)), fresh_atoms,
                    fresh_bonds, tuple(sorted(radicals.items())), tuple(sorted(
                        (a, b, o) for (a, b), o in orders.items())), fresh_coordinates)
    verify_end(witness, record, extended)
    return extended, record


def verify_end(original: Witness, record: EndMap, built: Witness) -> None:
    """Re-check an end reconstruction from scratch: survivors keep their ids,
    elements, electrons (charge, radical, lone pairs: only the recorded
    radical increments) and bonds (only the recorded orders set; no other
    bond appears or vanishes), the deleted atoms are the port's declared caps
    plus atoms outside the edit support, fresh atoms are new ids with the
    recorded attributes, bonded as recorded, the edits are transported
    unchanged and the metadata fits."""
    deleted = set(record.deleted)
    if not set(record.old_port.cap_atom_ids) <= deleted:
        raise MapCheckError("an end reconstruction must delete the port's declared caps")
    if deleted & edit_support(original):
        raise MapCheckError("an end reconstruction deleted a reacting atom")
    fresh = set(record.fresh_ids)
    for name, orig, new in (("reactants", original.reactant_rows(), built.reactant_rows()),
                            ("products", original.product_rows(), built.product_rows())):
        if fresh & set(orig):
            raise MapCheckError(f"{name}: fresh end ids collide with original atoms")
        if set(new) != (set(orig) - deleted) | fresh:
            raise MapCheckError(f"{name}: atoms are not survivors plus fresh end atoms")
        radicals = dict(record.radicals)
        for a in set(orig) - deleted:
            if orig[a][:2] != new[a][:2]:
                raise MapCheckError(f"{name}: atom {a} changed element or isotope")
            want = list(_attrs(orig[a]))
            want[3] += radicals.get(a, 0)
            if tuple(want) != _attrs(new[a]):
                raise MapCheckError(f"{name}: surviving atom {a} changed its electrons "
                                    f"(charge, radical, lone pairs) beyond the recorded end")
        for a, attrs in record.fresh_atoms:
            if _attrs(new[a]) != tuple(attrs):
                raise MapCheckError(f"{name}: fresh end atom {a} is not as recorded")
        set_orders = {(min(a, c), max(a, c)): float(o) for a, c, o in record.bonds}
        old_bonds, new_bonds = _bonds(orig), _bonds(new)
        want_bonds = {pair: set_orders.get(pair, o) for pair, o in old_bonds.items()
                      if pair[0] not in deleted and pair[1] not in deleted}
        if set(set_orders) - set(want_bonds):
            raise MapCheckError(f"{name}: a bond the end sets is absent")
        want_bonds.update({(min(a, c), max(a, c)): float(o) for a, c, o in record.fresh_bonds})
        if new_bonds != want_bonds:
            changed = sorted(set(new_bonds.items()) ^ set(want_bonds.items()))[:4]
            raise MapCheckError(f"{name}: bonds differ from the survivors' plus the recorded "
                                f"end bonds: {changed}")
    if edits(built) != edits(original):
        raise MapCheckError("the end reconstruction changed the reaction's edits")
    for participant in built.participants:
        participant.metadata.check(built.reactant_rows(),
                                   f"end reconstruction (participant {participant.index})",
                                   participant.atom_ids)
    if built.port(record.port_key) != record.new_port or record.new_port.artificial:
        raise MapCheckError("the reconstructed port is not the recorded physical port")
    before = original.participant(record.port_key[0]).metadata.coordinates
    after = built.participant(record.port_key[0]).metadata.coordinates
    if set(after) != set(before) | set(record.fresh_coordinates) or \
            {a for a, _u, _k in record.fresh_coordinates} - fresh:
        raise MapCheckError("coordinates are not the survivors' plus the recorded fresh ones")


def ps_build(name: str, units: int) -> ProxyBuild:
    """Build one named ``ps_proxy_set`` structure with its metadata."""
    return PROVIDERS["polystyrene"].build(name, ps_structure_specs(units)[name])
