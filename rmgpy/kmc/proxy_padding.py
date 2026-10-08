"""Graph-native completion of explicitly annotated polymer proxy boundaries."""

from __future__ import annotations

import copy
import os
from collections import deque
from dataclasses import dataclass
from typing import Any, Iterable, Sequence


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
    """One annotated artificial continuation point on a participant graph."""

    atom_id: int
    orientation: str
    participant_index: int = 0
    kind: str = "artificial"


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

    def provenance(
        self, repeat_unit: RepeatUnitGraph, min_distance: int, tolerance_log10: float
    ) -> dict[str, Any]:
        return {
            "status": "padded",
            "minimum_heavy_bond_distance": min_distance,
            "tolerance_max_abs_log10": tolerance_log10,
            "repeat_unit": repeat_unit.provenance(),
            "reacting_atom_ids": list(self.reacting_atom_ids),
            "artificial_boundaries": [
                {
                    "participant_index": item.participant_index,
                    "atom_id": item.atom_id,
                    "orientation": item.orientation,
                }
                for item in self.boundaries
            ],
            "repeat_extensions": self.extensions,
            "original_heavy_atoms": self.original_heavy_atoms,
            "padded_heavy_atoms": self.padded_heavy_atoms,
        }


class PaddingLimitExceeded(ValueError):
    """The requested completion exceeded its pinned graph-size ceiling."""


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


def _remove_cap_hydrogen(molecule, atom) -> None:
    hydrogens = [neighbor for neighbor in atom.edges if neighbor.element.number == 1]
    if not hydrogens:
        raise ValueError("continuation port has no hydrogen cap to replace")
    molecule.remove_atom(hydrogens[0])


def _append_repeat(
    molecule,
    boundary: BoundaryPort,
    repeat_unit: RepeatUnitGraph,
    *,
    atom_ids: Sequence[int] | None = None,
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

    _remove_cap_hydrogen(molecule, old_port)
    _remove_cap_hydrogen(unit, attach)
    multiplicity = molecule.multiplicity
    molecule = molecule.merge(unit)
    molecule.multiplicity = multiplicity
    molecule.add_bond(Bond(old_port, attach, order="S"))
    molecule.update(sort_atoms=False)
    return molecule, BoundaryPort(
        atom_id=new_port.id,
        orientation=boundary.orientation,
        participant_index=boundary.participant_index,
        kind=boundary.kind,
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
    if os.environ.get("RMG_KMC_DISABLE_PROXY_PADDING") == "1":
        return PaddedMolecule(padded, tuple(boundaries), 0)
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
) -> PaddedReactionWitness:
    """Pad a reaction bijectively while keeping its executable counterpart separate."""
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
    if os.environ.get("RMG_KMC_DISABLE_PROXY_PADDING") == "1":
        return PaddedReactionWitness(
            witness,
            tuple(boundaries),
            centers,
            0,
            original_heavy_atoms,
            original_heavy_atoms,
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
                )
                product_molecule, product_port = _append_repeat(
                    product_molecule,
                    BoundaryPort(
                        current.atom_id,
                        current.orientation,
                        product_index,
                        current.kind,
                    ),
                    repeat_unit,
                    atom_ids=assigned_ids,
                )
                if reactant_port.atom_id != product_port.atom_id:
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
    )
