"""Graph-native completion of explicitly annotated polymer proxy boundaries."""

from __future__ import annotations

import copy
from collections import deque
from dataclasses import dataclass
from typing import Any, Iterable, Sequence

from rmgpy.kmc.resource_errors import reraise_resource


@dataclass(frozen=True)
class RepeatUnitGraph:
    """A repeat-unit graph with oriented head and tail continuation ports."""

    name: str
    molecule: Any
    head_atom_index: int
    tail_atom_index: int

    @classmethod
    def from_smiles(
        cls, name: str, smiles: str, *, head_atom_index: int, tail_atom_index: int
    ) -> "RepeatUnitGraph":
        from rmgpy.molecule.molecule import Molecule

        return cls(
            name=name,
            molecule=Molecule(smiles=smiles),
            head_atom_index=head_atom_index,
            tail_atom_index=tail_atom_index,
        )

    def provenance(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "adjacency_list": self.molecule.to_adjacency_list(remove_h=False),
            "head_atom_index": self.head_atom_index,
            "tail_atom_index": self.tail_atom_index,
        }


@dataclass(frozen=True)
class BoundaryPort:
    """One annotated artificial continuation point on a participant graph.

    ``cap_atom_ids`` are the cap atoms declared at construction (see
    ``rmgpy.kmc.boundaries``); empty means the legacy first-neighbouring-H
    choice.  ``source`` flags where the port came from: ``construction``
    (provider metadata) or ``inferred`` (legacy graph/name inference)."""

    atom_id: int
    orientation: str
    participant_index: int = 0
    kind: str = "artificial"
    cap_atom_ids: tuple = ()
    source: str = "inferred"


@dataclass(frozen=True)
class PaddedMolecule:
    molecule: Any
    boundaries: tuple[BoundaryPort, ...]
    extensions: int


@dataclass(frozen=True)
class PaddedReactionWitness:
    """A padded rate/thermo reaction plus auditable completion metadata."""

    reaction: Any
    boundaries: tuple[BoundaryPort, ...]
    reacting_atom_ids: tuple[int, ...]
    extensions: int
    original_heavy_atoms: int
    padded_heavy_atoms: int
    executable_to_witness_projection: dict[str, list[dict[str, int]]]

    def provenance(
        self, repeat_unit: RepeatUnitGraph, min_distance: int, tolerance_log10: float
    ) -> dict[str, Any]:
        return {
            "status": "padded",
            "minimum_heavy_bond_distance": min_distance,
            "tolerance_max_abs_log10": tolerance_log10,
            "repeat_unit": repeat_unit.provenance(),
            "reacting_atoms": {
                side: _canonical_atom_references(
                    getattr(self.reaction, side), self.reacting_atom_ids
                )
                for side in ("reactants", "products")
            },
            "artificial_boundaries": [
                {
                    side: _canonical_atom_reference(
                        getattr(self.reaction, side), item.atom_id
                    )
                    for side in ("reactants", "products")
                } | {
                    "orientation": item.orientation,
                }
                for item in self.boundaries
            ],
            "executable_to_witness_projection": copy.deepcopy(
                self.executable_to_witness_projection
            ),
            "repeat_extensions": self.extensions,
            "original_heavy_atoms": self.original_heavy_atoms,
            "padded_heavy_atoms": self.padded_heavy_atoms,
        }


class PaddingLimitExceeded(ValueError):
    """The requested completion exceeded its pinned graph-size ceiling."""


class ReactionNotPaddable(ValueError):
    """An annotated artificial boundary cannot be replaced for this event."""


def _canonical_atoms(molecule) -> list[Any]:
    """Return an atom order determined by graph structure, never RMG atom IDs."""
    if molecule.__class__.__module__.startswith("rmgpy."):
        from rdkit import Chem
        from rmgpy.molecule.converter import to_rdkit_mol

        rdkit_molecule, mapping = to_rdkit_mol(
            molecule, remove_h=False, return_mapping=True, save_order=True
        )
        ranks = Chem.CanonicalRankAtoms(
            rdkit_molecule, breakTies=True, includeChirality=True,
            includeIsotopes=True,
        )
        return sorted(molecule.atoms, key=lambda atom: ranks[mapping[atom]])
    return list(molecule.atoms)


def _canonical_heavy_index(molecule, atom_id: int) -> int:
    heavy = [atom for atom in _canonical_atoms(molecule) if atom.element.number != 1]
    matches = [index for index, atom in enumerate(heavy) if atom.id == atom_id]
    if len(matches) != 1:
        raise ValueError("canonical atom reference is absent or ambiguous")
    return matches[0]


def _canonical_atom_references(participants, atom_ids) -> list[dict[str, int]]:
    references = []
    for atom_id in atom_ids:
        participant_index = _participant_containing(participants, atom_id)
        references.append({
            "participant_index": participant_index,
            "atom_index": _canonical_heavy_index(
                _molecule(participants[participant_index]), atom_id
            ),
        })
    return sorted(
        references,
        key=lambda item: (item["participant_index"], item["atom_index"]),
    )


def _canonical_atom_reference(participants, atom_id) -> dict[str, int]:
    participant_index = _participant_containing(participants, atom_id)
    return {
        "participant_index": participant_index,
        "atom_index": _canonical_heavy_index(
            _molecule(participants[participant_index]), atom_id
        ),
    }


def _projection(executable, witness) -> dict[str, list[dict[str, int]]]:
    result = {}
    for side_name in ("reactants", "products"):
        source_side = getattr(executable, side_name)
        witness_side = getattr(witness, side_name)
        entries = []
        for participant_index, source_participant in enumerate(source_side):
            source_molecule = _molecule(source_participant)
            witness_molecule = _molecule(witness_side[participant_index])
            witness_ids = {
                atom.id for atom in witness_molecule.atoms
                if atom.element.number != 1
            }
            for source_atom in _canonical_atoms(source_molecule):
                if source_atom.element.number == 1:
                    continue
                if source_atom.id not in witness_ids:
                    raise ValueError("executable atom is absent from padded witness")
                entries.append({
                    "participant_index": participant_index,
                    "witness_participant_index": participant_index,
                    "executable_atom_index": _canonical_heavy_index(
                        source_molecule, source_atom.id
                    ),
                    "witness_atom_index": _canonical_heavy_index(
                        witness_molecule, source_atom.id
                    ),
                })
        result[side_name] = entries
    return result


def heavy_atom_distance(molecule, first_id: int, second_id: int) -> int:
    """Return the heavy-bond distance between two atom IDs."""
    by_id = {atom.id: atom for atom in molecule.atoms if atom.element.number != 1}
    if first_id not in by_id or second_id not in by_id:
        raise ValueError("distance endpoint is absent from the heavy-atom graph")
    target = by_id[second_id]
    seen = {by_id[first_id]}
    frontier = deque([(by_id[first_id], 0)])
    while frontier:
        atom, distance = frontier.popleft()
        if atom is target:
            return distance
        for neighbor in atom.edges:
            if neighbor.element.number != 1 and neighbor not in seen:
                seen.add(neighbor)
                frontier.append((neighbor, distance + 1))
    raise ValueError("distance endpoints are disconnected")


def _remove_cap_hydrogen(molecule, atom, cap_atom_ids=(), support=None) -> None:
    """Remove the cap of continuation port ``atom``: the declared cap atoms
    (each must be a hydrogen bonded to the port by a single bond and only to
    it) or, without a declaration, the first hydrogen of the port.  A cap in
    the reaction's edit ``support`` (atom ids whose attributes or bonds
    change, hydrogens included) is never removed."""
    support = set(support or ())
    if cap_atom_ids:
        wanted = set(cap_atom_ids)
        caps = [neighbor for neighbor in atom.edges if neighbor.id in wanted]
        if len(caps) != len(wanted):
            raise ReactionNotPaddable(
                "declared cap atom is not bonded to its continuation port"
            )
        for cap in caps:
            if cap.element.number != 1:
                raise ReactionNotPaddable(
                    f"declared cap atom {cap.id} is not a hydrogen"
                )
            if len(cap.edges) != 1 or float(cap.edges[atom].order) != 1.0:
                raise ReactionNotPaddable(
                    f"declared cap hydrogen {cap.id} is not singly bonded to its port only"
                )
            if cap.id in support:
                raise ReactionNotPaddable(
                    f"declared cap hydrogen {cap.id} is in the reaction's edit support"
                )
        for cap in caps:
            molecule.remove_atom(cap)
        return
    hydrogens = [
        neighbor
        for neighbor in atom.edges
        if neighbor.element.number == 1 and neighbor.id not in support
    ]
    if not hydrogens:
        raise ValueError("continuation port has no hydrogen cap to replace")
    molecule.remove_atom(hydrogens[0])


def _append_repeat(
    molecule,
    boundary: BoundaryPort,
    repeat_unit: RepeatUnitGraph,
    *,
    atom_ids: Sequence[int] | None = None,
    support=None,
) -> tuple[Any, BoundaryPort]:
    from rmgpy.molecule.molecule import Bond

    by_id = {atom.id: atom for atom in molecule.atoms}
    try:
        old_port = by_id[boundary.atom_id]
    except KeyError as error:
        raise ValueError("annotated continuation port is absent from graph") from error
    if old_port.element.number == 1:
        raise ValueError("continuation port must be a heavy atom")

    unit = repeat_unit.molecule.copy(deep=True)
    heavy = [atom for atom in unit.atoms if atom.element.number != 1]
    try:
        head = heavy[repeat_unit.head_atom_index]
        tail = heavy[repeat_unit.tail_atom_index]
    except IndexError as error:
        raise ValueError("repeat-unit continuation index is outside its heavy graph") from error
    if head is tail:
        raise ValueError("repeat-unit head and tail ports must be distinct")
    if boundary.orientation == "tail":
        attach, new_port = head, tail
    elif boundary.orientation == "head":
        attach, new_port = tail, head
    else:
        raise ValueError("boundary orientation must be 'head' or 'tail'")

    used_ids = {atom.id for atom in molecule.atoms}
    if atom_ids is None:
        next_id = max(used_ids, default=-1) + 1
        assigned_ids = []
        for _ in unit.atoms:
            while next_id in used_ids:
                next_id += 1
            assigned_ids.append(next_id)
            used_ids.add(next_id)
            next_id += 1
    else:
        assigned_ids = list(atom_ids)
        if len(assigned_ids) != len(unit.atoms):
            raise ValueError("repeat atom-ID assignment has the wrong length")
        if len(set(assigned_ids)) != len(assigned_ids) or used_ids.intersection(
            assigned_ids
        ):
            raise ValueError("repeat atom-ID assignment is not fresh")
    for atom, atom_id in zip(unit.atoms, assigned_ids):
        atom.id = atom_id

    _remove_cap_hydrogen(molecule, old_port, boundary.cap_atom_ids, support)
    _remove_cap_hydrogen(unit, attach)
    new_caps = ()
    if boundary.cap_atom_ids:
        # declared caps continue as declared: the last H of the new port in
        # repeat-unit atom order (identical on both reaction sides)
        unit_hydrogens = [
            neighbor.id for neighbor in new_port.edges if neighbor.element.number == 1
        ]
        if not unit_hydrogens:
            raise ReactionNotPaddable("new continuation port has no hydrogen cap")
        new_caps = (unit_hydrogens[-1],)
    multiplicity = molecule.multiplicity
    molecule = molecule.merge(unit)
    molecule.multiplicity = multiplicity
    molecule.add_bond(Bond(old_port, attach, order="S"))
    molecule.update(sort_atoms=False)
    molecule.multiplicity = multiplicity
    return molecule, BoundaryPort(
        atom_id=new_port.id,
        orientation=boundary.orientation,
        participant_index=boundary.participant_index,
        kind=boundary.kind,
        cap_atom_ids=new_caps,
        source=boundary.source,
    )


def pad_molecule(
    molecule,
    boundaries: Sequence[BoundaryPort],
    repeat_unit: RepeatUnitGraph,
    reacting_atom_ids: Iterable[int],
    min_distance: int,
) -> PaddedMolecule:
    """Extend artificial ports until each is sufficiently far from all reacting atoms."""
    if min_distance < 0:
        raise ValueError("minimum boundary distance must be non-negative")
    padded = copy.deepcopy(molecule)
    reacting_ids = tuple(reacting_atom_ids)
    available = {
        atom.id for atom in padded.atoms if atom.element.number != 1
    }
    if not set(reacting_ids).issubset(available):
        raise ValueError("reacting atom is absent from the padded molecule")
    updated = []
    extensions = 0
    for boundary in boundaries:
        current = boundary
        if current.kind == "artificial" and reacting_ids:
            while min(
                heavy_atom_distance(padded, atom_id, current.atom_id)
                for atom_id in reacting_ids
            ) < min_distance:
                padded, current = _append_repeat(padded, current, repeat_unit)
                extensions += 1
        updated.append(current)
    return PaddedMolecule(padded, tuple(updated), extensions)


def _molecule(participant):
    molecules = getattr(participant, "molecule", None)
    return molecules[0] if molecules else participant


def _set_molecule(participant, molecule):
    if getattr(participant, "molecule", None) is not None:
        participant.molecule[0] = molecule
    else:
        return molecule
    return participant


def _side_atoms(participants: Sequence[Any]) -> dict[int, Any]:
    atoms = {}
    for participant in participants:
        for atom in _molecule(participant).atoms:
            if atom.element.number == 1:
                continue
            if atom.id in atoms:
                raise ValueError("reaction side contains duplicate heavy-atom IDs")
            atoms[atom.id] = atom
    return atoms


def reaction_center_atom_ids(reaction) -> tuple[int, ...]:
    """Identify atoms whose bonds or electronic attributes change across a reaction."""
    reactants = _side_atoms(reaction.reactants)
    products = _side_atoms(reaction.products)
    if set(reactants) != set(products):
        raise ValueError("reaction witness does not have a bijective heavy-atom map")

    def bonds(participants):
        result = {}
        for participant in participants:
            for bond in _molecule(participant).get_all_edges():
                if bond.vertex1.element.number == 1 or bond.vertex2.element.number == 1:
                    continue
                pair = tuple(sorted((bond.vertex1.id, bond.vertex2.id)))
                result[pair] = float(bond.order)
        return result

    before = bonds(reaction.reactants)
    after = bonds(reaction.products)
    reacting = {
        atom_id
        for pair in set(before) | set(after)
        if before.get(pair) != after.get(pair)
        for atom_id in pair
    }
    attributes = ("radical_electrons", "charge", "lone_pairs")
    for atom_id in reactants:
        if any(
            getattr(reactants[atom_id], name, 0)
            != getattr(products[atom_id], name, 0)
            for name in attributes
        ):
            reacting.add(atom_id)
    return tuple(sorted(reacting))


def _participant_containing(participants: Sequence[Any], atom_id: int) -> int:
    matches = [
        index
        for index, participant in enumerate(participants)
        if any(
            atom.id == atom_id and atom.element.number != 1
            for atom in _molecule(participant).atoms
        )
    ]
    if len(matches) != 1:
        raise ValueError(
            "artificial boundary must identify exactly one participant graph"
        )
    return matches[0]


def _heavy_count(participants: Sequence[Any]) -> int:
    return sum(
        atom.element.number != 1
        for participant in participants
        for atom in _molecule(participant).atoms
    )


def pad_reaction_witness(
    reaction,
    boundaries: Sequence[BoundaryPort],
    repeat_unit: RepeatUnitGraph,
    min_distance: int,
    *,
    reacting_atom_ids: Iterable[int] | None = None,
    max_heavy_atoms: int = 512,
    allow_unverified: bool = False,
) -> PaddedReactionWitness:
    """Pad a reaction bijectively while keeping its executable counterpart separate.

    Legacy (unverified) path for inferred ports.  Construction-annotated
    ports carry declared caps and metadata and must be padded through the
    verified explicit-H extension (``pad_reaction_witness_verified``);
    they are refused here unless ``allow_unverified`` (unit tests of this
    helper only).  Caps in the reaction's edit support are never removed."""
    if not allow_unverified and any(
        getattr(boundary, "source", "inferred") == "construction" for boundary in boundaries
    ):
        raise ReactionNotPaddable(
            "construction-annotated ports require the verified extension path "
            "(pad_reaction_witness_verified)"
        )
    support = _edit_support_ids(reaction)
    witness = copy.deepcopy(reaction)
    centers = tuple(
        sorted(
            reaction_center_atom_ids(witness)
            if reacting_atom_ids is None
            else set(reacting_atom_ids)
        )
    )
    reactant_atoms = _side_atoms(witness.reactants)
    product_atoms = _side_atoms(witness.products)
    if set(reactant_atoms) != set(product_atoms):
        raise ValueError("reaction witness does not have a bijective heavy-atom map")
    if not set(centers).issubset(reactant_atoms):
        raise ValueError("reacting atom is absent from the reaction witness")
    original_heavy_atoms = len(reactant_atoms)
    if original_heavy_atoms > max_heavy_atoms:
        raise PaddingLimitExceeded("unpadded witness exceeds maximum graph size")
    reacting_boundaries = sorted(
        boundary.atom_id for boundary in boundaries
        if boundary.kind == "artificial" and boundary.atom_id in centers
    )
    if reacting_boundaries:
        raise ReactionNotPaddable(
            "artificial continuation port is a reacting atom"
        )
    updated = []
    extensions = 0
    all_ids = set(reactant_atoms) | {
        atom.id
        for participant in tuple(witness.reactants) + tuple(witness.products)
        for atom in _molecule(participant).atoms
    }
    next_id = max(all_ids, default=-1) + 1
    for boundary in boundaries:
        reactant_index = _participant_containing(
            witness.reactants, boundary.atom_id
        )
        product_index = _participant_containing(witness.products, boundary.atom_id)
        reactant_molecule = _molecule(witness.reactants[reactant_index])
        product_molecule = _molecule(witness.products[product_index])
        reactant_centers = tuple(
            atom_id
            for atom_id in centers
            if any(atom.id == atom_id for atom in reactant_molecule.atoms)
        )
        product_centers = tuple(
            atom_id
            for atom_id in centers
            if any(atom.id == atom_id for atom in product_molecule.atoms)
        )
        current = BoundaryPort(
            boundary.atom_id,
            boundary.orientation,
            reactant_index,
            boundary.kind,
            boundary.cap_atom_ids,
            boundary.source,
        )
        if current.kind == "artificial" and (
            reactant_centers or product_centers
        ):
            def side_too_close():
                return (
                    reactant_centers
                    and min(
                        heavy_atom_distance(
                            reactant_molecule, atom_id, current.atom_id
                        )
                        for atom_id in reactant_centers
                    )
                    < min_distance
                ) or (
                    product_centers
                    and min(
                        heavy_atom_distance(
                            product_molecule, atom_id, current.atom_id
                        )
                        for atom_id in product_centers
                    )
                    < min_distance
                )

            while side_too_close():
                repeat_size = len(repeat_unit.molecule.atoms)
                assigned_ids = tuple(range(next_id, next_id + repeat_size))
                next_id += repeat_size
                if (
                    _heavy_count(witness.reactants)
                    + sum(
                        atom.element.number != 1
                        for atom in repeat_unit.molecule.atoms
                    )
                    > max_heavy_atoms
                ):
                    raise PaddingLimitExceeded(
                        "padding did not converge within maximum graph size"
                    )
                reactant_molecule, reactant_port = _append_repeat(
                    reactant_molecule,
                    current,
                    repeat_unit,
                    atom_ids=assigned_ids,
                    support=support,
                )
                product_molecule, product_port = _append_repeat(
                    product_molecule,
                    BoundaryPort(
                        current.atom_id,
                        current.orientation,
                        product_index,
                        current.kind,
                        current.cap_atom_ids,
                        current.source,
                    ),
                    repeat_unit,
                    atom_ids=assigned_ids,
                    support=support,
                )
                if (reactant_port.atom_id, reactant_port.cap_atom_ids) != (
                    product_port.atom_id, product_port.cap_atom_ids
                ):
                    raise ValueError("padded atom correspondence diverged")
                witness.reactants[reactant_index] = _set_molecule(
                    witness.reactants[reactant_index], reactant_molecule
                )
                witness.products[product_index] = _set_molecule(
                    witness.products[product_index], product_molecule
                )
                current = BoundaryPort(
                    reactant_port.atom_id,
                    current.orientation,
                    reactant_index,
                    current.kind,
                    reactant_port.cap_atom_ids,
                    current.source,
                )
                extensions += 1
        updated.append(current)

    padded_heavy_atoms = _heavy_count(witness.reactants)
    if set(_side_atoms(witness.reactants)) != set(_side_atoms(witness.products)):
        raise ValueError("padding did not preserve the heavy-atom bijection")
    return PaddedReactionWitness(
        witness,
        tuple(updated),
        centers,
        extensions,
        original_heavy_atoms,
        padded_heavy_atoms,
        _projection(reaction, witness),
    )


def _edit_support_ids(reaction) -> frozenset:
    """Atom ids (hydrogens included) whose attributes or bonds differ across
    the reaction (one id namespace on both sides); empty when the sides do
    not share ids."""
    from rmgpy.kmc.rmg_moves import molecule_rows

    def rows(side):
        out = {}
        for participant in side:
            out.update(molecule_rows(_molecule(participant)))
        return out

    try:
        before, after = rows(reaction.reactants), rows(reaction.products)
    except Exception as error:
        reraise_resource(error)
        return frozenset()
    if set(before) != set(after):
        return frozenset()
    return frozenset(a for a in before if before[a] != after[a])


@dataclass(frozen=True)
class VerifiedPaddedReactionWitness(PaddedReactionWitness):
    """A padded witness built by the verified explicit-H extension path:
    ``extension_map`` is the composite ``boundaries.ExtensionMap`` and
    ``verified`` records that ``boundaries.verify_extension`` (history,
    bookkeeping, replay of every step) accepted it."""

    extension_map: Any = None
    verified: bool = False

    def provenance(
        self, repeat_unit: RepeatUnitGraph, min_distance: int, tolerance_log10: float
    ) -> dict[str, Any]:
        out = super().provenance(repeat_unit, min_distance, tolerance_log10)
        record = self.extension_map
        out["extension"] = {
            "method": "rmgpy.kmc.boundaries ExtensionMap (explicit H), verify_extension",
            "verified": bool(self.verified),
            "steps": [
                [step.port_key[0], step.port_key[1], step.depth, step.provider]
                for step in record.steps
            ],
            "deleted_caps": len(record.deleted),
            "fresh_atoms": len(record.fresh_atoms),
        }
        return out


_ORIENTATION_OF_TAG = {"backbone_head": "head", "backbone_tail": "tail"}


def _fragment_molecule(fragment):
    """The closed repeat unit one provider step adds: the fragment's atoms
    and bonds (local ids) with one hydrogen standing in for each attachment
    bond to the old port."""
    from rmgpy.kmc.rmg_moves import molecule_from_rows

    edges = {k: [] for k in range(len(fragment.atoms))}
    for i, j, order in fragment.bonds:
        edges[i].append((j, float(order)))
        edges[j].append((i, float(order)))
    next_id = len(fragment.atoms)
    extra = {}
    for local, _order in fragment.attachments:
        edges[local].append((next_id, 1.0))
        extra[next_id] = ("H", None, 0, 0, 0, ((local, 1.0),))
        next_id += 1
    rows = {k: tuple(atom) + (tuple(sorted(edges[k])),)
            for k, atom in enumerate(fragment.atoms)}
    rows.update(extra)
    return molecule_from_rows(rows)


def check_padding_declarations(original, boundaries, repeat_unit, providers) -> None:
    """The declared padding (``boundaries``: orientation, caps, source;
    ``repeat_unit``: the repeat chemistry and its head/tail atoms named in
    the provenance) must be what the transported construction metadata of
    ``original`` (a ``boundaries.Witness``) and its providers will execute;
    raise ``ReactionNotPaddable`` on any disagreement.  Per artificial port:
    a construction-sourced boundary on the port atom, the orientation of the
    port's tag, exactly the declared caps, and the provider's one-step
    fragment, closed with a hydrogen on each attachment, isomorphic to the
    declared repeat unit (constitution, hydrogens explicit) by a map taking
    the attachment atom to the unit's continuing end and the new port to its
    other end (tail extension: attachment = head atom, new port = tail atom;
    head extension: the reverse)."""
    if repeat_unit is None:
        raise ReactionNotPaddable("no declared repeat unit to check the providers against")
    by_atom = {}
    for key in original.artificial_port_keys():
        by_atom[original.port(key).atom_id] = key
    unit = repeat_unit.molecule
    unit_heavy = [atom for atom in unit.atoms if atom.element.number != 1]
    try:
        unit_head = unit_heavy[repeat_unit.head_atom_index]
        unit_tail = unit_heavy[repeat_unit.tail_atom_index]
    except IndexError:
        raise ReactionNotPaddable(
            f"declared repeat unit {repeat_unit.name!r}: head/tail atom index out of range"
        ) from None
    for boundary in boundaries:
        if boundary.kind != "artificial":
            continue
        key = by_atom.get(boundary.atom_id)
        if key is None:
            raise ReactionNotPaddable(
                f"declared port atom {boundary.atom_id} is not a transported artificial port"
            )
        port = original.port(key)
        if getattr(boundary, "source", "inferred") != "construction":
            raise ReactionNotPaddable(f"port {key}: the declared boundary is not construction-sourced")
        expected = _ORIENTATION_OF_TAG.get(port.tag)
        if boundary.orientation != expected:
            raise ReactionNotPaddable(
                f"port {key}: declared orientation {boundary.orientation!r} is not the "
                f"{expected!r} of its transported tag {port.tag!r}"
            )
        if tuple(sorted(boundary.cap_atom_ids)) != tuple(port.cap_atom_ids):
            raise ReactionNotPaddable(
                f"port {key}: declared cap {tuple(boundary.cap_atom_ids)} is not the transported "
                f"cap {port.cap_atom_ids}"
            )
        metadata = original.participant(key[0]).metadata
        try:
            provider = providers[metadata.provider]
            fragment = provider.fragment(port, metadata.coordinate(port.atom_id))
        except Exception as error:
            reraise_resource(error)
            raise ReactionNotPaddable(f"port {key}: provider fragment unavailable: {error}") \
                from error
        if len(fragment.attachments) != 1:
            raise ReactionNotPaddable(
                f"port {key}: a repeat unit graph declares one attachment, the provider has "
                f"{len(fragment.attachments)}"
            )
        closed = _fragment_molecule(fragment)
        attach = closed.atoms[fragment.attachments[0][0]]
        new_port = closed.atoms[fragment.new_port]
        if expected == "tail":
            initial = {attach: unit_head, new_port: unit_tail}
        else:
            initial = {attach: unit_tail, new_port: unit_head}
        if attach is new_port:
            initial = {attach: unit_head} if unit_head is unit_tail else None
        # save_order: ``repeat_unit`` is shared across proxies and its head and
        # tail are atom INDICES; RMG's default sorts both graphs' atoms
        if initial is None or not closed.is_isomorphic(unit, initial_map=initial, strict=False,
                                                        save_order=True):
            raise ReactionNotPaddable(
                f"port {key}: declared repeat unit {repeat_unit.name!r} is not the chemistry "
                f"(or orientation) its provider {metadata.provider!r} adds"
            )


def pad_reaction_witness_verified(
    reaction,
    boundaries: Sequence[BoundaryPort],
    participants: Sequence[Any],
    min_distance: int,
    *,
    max_heavy_atoms: int = 512,
    providers=None,
    repeat_unit: RepeatUnitGraph | None = None,
) -> VerifiedPaddedReactionWitness:
    """Pad a reaction through ``rmgpy.kmc.boundaries``: the reaction is read
    as an explicit-H witness (both sides one id namespace) carrying the
    transported construction metadata (``participants``); each artificial
    port is extended one verified step at a time (declared caps deleted on
    both sides, edit support protected) until it is ``min_distance`` heavy
    bonds from every reacting heavy atom on both sides; the composite map is
    re-verified from scratch (``verify_extension``) and the witness is
    converted back onto the reaction's participants.  Before anything is
    extended the declarations are checked against the transported metadata
    (``check_padding_declarations``): ``repeat_unit`` (the chemistry the
    provenance names), every port's orientation and caps."""
    from rmgpy.kmc import boundaries as bnd
    from rmgpy.kmc.rmg_moves import molecule_from_rows, molecule_rows, split_rows

    if min_distance < 0:
        raise ValueError("minimum boundary distance must be non-negative")
    if not participants:
        raise ReactionNotPaddable(
            "construction metadata was not transported: the verified extension path "
            "needs boundary participants"
        )
    witness_reaction = copy.deepcopy(reaction)
    sides = {
        "reactants": list(witness_reaction.reactants),
        "products": list(witness_reaction.products),
    }
    rows = {
        name: [molecule_rows(_molecule(p)) for p in side] for name, side in sides.items()
    }
    owned = {p.index: frozenset(p.atom_ids) for p in participants}
    side_a = None
    for name in ("reactants", "products"):
        sets = [frozenset(r) for r in rows[name]]
        if sorted(sets, key=sorted) == sorted(owned.values(), key=sorted):
            side_a = name
            break
    if side_a is None:
        raise ReactionNotPaddable(
            "the transported construction participants do not match a reaction side"
        )
    side_b = "products" if side_a == "reactants" else "reactants"
    flat_a = {a: r for part in rows[side_a] for a, r in part.items()}
    flat_b = {a: r for part in rows[side_b] for a, r in part.items()}
    if set(flat_a) != set(flat_b):
        raise ReactionNotPaddable(
            "explicit-H witness: the reaction sides do not share one atom-id namespace"
        )
    try:
        original = bnd.make_witness(
            tuple(sorted(flat_a.items())), tuple(sorted(flat_b.items())),
            participants, label="production padding", where="production padding",
        )
    except bnd.BoundaryError as error:
        raise ReactionNotPaddable(f"explicit-H witness refused: {error}") from error
    centers = tuple(sorted(reaction_center_atom_ids(witness_reaction)))
    original_heavy = original.heavy_count()
    if original_heavy > max_heavy_atoms:
        raise PaddingLimitExceeded("unpadded witness exceeds maximum graph size")
    by_atom = {}
    for key in original.artificial_port_keys():
        by_atom[original.port(key).atom_id] = key
    if providers is None:
        providers = bnd.PROVIDERS
    check_padding_declarations(original, boundaries, repeat_unit, providers)
    if set(by_atom) & set(centers):
        raise ReactionNotPaddable("artificial continuation port is a reacting atom")
    ports = {
        (boundary.atom_id): boundary
        for boundary in boundaries
        if boundary.kind == "artificial"
    }
    if set(ports) != set(by_atom):
        raise ReactionNotPaddable(
            "the artificial ports do not match the transported construction metadata"
        )

    def heavy_distance(side_rows, source):
        distance = {source: 0}
        frontier = deque([source])
        while frontier:
            atom = frontier.popleft()
            for neighbour, _order in side_rows[atom][5]:
                if neighbour not in distance and side_rows[neighbour][0] != "H":
                    distance[neighbour] = distance[atom] + 1
                    frontier.append(neighbour)
        return distance

    def too_close(witness, key):
        port_atom = witness.port(key).atom_id
        for side_rows in (witness.reactant_rows(), witness.product_rows()):
            distance = heavy_distance(side_rows, port_atom)
            reach = [distance[c] for c in centers if c in distance]
            if reach and min(reach) < min_distance:
                return True
        return False

    current = original
    steps = []
    final_ports = {}
    for key in original.artificial_port_keys():
        while too_close(current, key):
            try:
                current, record = bnd.extend(
                    current, key, providers, ceiling=max_heavy_atoms, verify=True
                )
            except bnd.CeilingExceeded as error:
                raise PaddingLimitExceeded(
                    f"padding did not converge within maximum graph size: {error}"
                ) from error
            except bnd.BoundaryError as error:
                raise ReactionNotPaddable(f"verified extension refused: {error}") from error
            steps.append(record)
        final_ports[key] = current.port(key)
    composite = bnd.identity_map(original)
    for record in steps:
        composite = record if not composite.steps else bnd.compose(composite, record)
    try:
        bnd.verify_extension(original, composite, current, providers)
    except bnd.BoundaryError as error:
        raise ReactionNotPaddable(f"extension map verification failed: {error}") from error

    def rebuild(name, side_rows):
        components = split_rows(dict(side_rows))
        used = set()
        for position, participant in enumerate(sides[name]):
            molecule = _molecule(participant)
            ids = {atom.id for atom in molecule.atoms}
            hits = [k for k, component in enumerate(components) if ids & set(component)]
            if len(hits) != 1 or hits[0] in used:
                raise ReactionNotPaddable("padded witness does not split onto the participants")
            used.add(hits[0])
            rebuilt = molecule_from_rows(components[hits[0]])
            rebuilt.multiplicity = molecule.multiplicity
            sides[name][position] = _set_molecule(participant, rebuilt)
        if len(used) != len(components):
            raise ReactionNotPaddable("padded witness has a component without a participant")

    rebuild(side_a, current.reactant_rows())
    rebuild(side_b, current.product_rows())
    witness_reaction.reactants = sides["reactants"]
    witness_reaction.products = sides["products"]
    updated = []
    for key in original.artificial_port_keys():
        boundary = ports[original.port(key).atom_id]
        port = final_ports[key]
        position = _participant_containing(witness_reaction.reactants, port.atom_id)
        updated.append(
            BoundaryPort(
                port.atom_id, boundary.orientation, position, boundary.kind,
                tuple(port.cap_atom_ids), boundary.source,
            )
        )
    if set(_side_atoms(witness_reaction.reactants)) != set(
        _side_atoms(witness_reaction.products)
    ):
        raise ValueError("padding did not preserve the heavy-atom bijection")
    return VerifiedPaddedReactionWitness(
        witness_reaction,
        tuple(updated),
        centers,
        len(steps),
        original_heavy,
        current.heavy_count(),
        _projection(reaction, witness_reaction),
        extension_map=composite,
        verified=True,
    )
