"""Observational accounting for RMG reaction-path merging.

This module deliberately owns no chemistry.  It drives the same generation,
merge, and correction operations as
``KineticsDatabase.generate_reactions_from_families`` while retaining the
raw-path evidence that RMG normally discards.
"""

from collections import OrderedDict
import copy
from dataclasses import dataclass
from typing import Any, Iterable, Sequence, Tuple

from rmgpy.data.kinetics.common import (
    check_for_same_reactants,
    ensure_independent_atom_ids,
    find_degenerate_reactions,
    generate_molecule_combos,
)


@dataclass(frozen=True)
class MappedGraphV1:
    """One participant graph plus the atom IDs aligned with its atom order."""

    participant_index: int
    adjacency_list: str
    atom_ids: Tuple[int, ...]


@dataclass(frozen=True, order=True)
class RootAtomV1:
    """A reaction-centre atom rooted in one reactant participant."""

    participant_index: int
    atom_id: int


@dataclass(frozen=True)
class RootedPathV1:
    """One raw RMG path before merge accounting mutates its degeneracy."""

    path_id: str
    degeneracy: float
    root: Tuple[RootAtomV1, ...]
    contributes_to_pre_correction_sum: bool = True


@dataclass(frozen=True)
class RootWeightV1:
    """The final degeneracy allocated to one distinct rooted event."""

    root: Tuple[RootAtomV1, ...]
    weight: float
    contributing_path_ids: Tuple[str, ...]


@dataclass(frozen=True)
class RMGMergeAccountingV1:
    """Auditable evidence for one reaction returned by RMG generation."""

    schema: str
    family: str
    direction: str
    family_own_reverse: bool
    reactant_graphs: Tuple[MappedGraphV1, ...]
    product_graphs: Tuple[MappedGraphV1, ...]
    raw_paths: Tuple[RootedPathV1, ...]
    merge_sublist_path_ids: Tuple[str, ...]
    contributing_path_ids: Tuple[str, ...]
    pre_correction_degeneracy: float
    generation_same_reactants: int
    same_reactants: int
    full_species_reactants_identical: bool
    identical_reactant_divisor: float
    reverse_degeneracy_recalculated: bool
    recalculated_rooted_paths: Tuple[RootedPathV1, ...]
    G_final: float
    root_weights: Tuple[RootWeightV1, ...]


@dataclass(frozen=True)
class MergeAccountingResult:
    """Generated reactions and their optional observational records."""

    reactions: list
    records: Tuple[RMGMergeAccountingV1, ...]


class MergeAccountingError(ValueError):
    """Raised when RMG's final degeneracy cannot be rooted without guessing."""


@dataclass
class _MergeGroup:
    contributors: list
    members: list


def _molecule(participant):
    molecules = getattr(participant, "molecule", None)
    return molecules[0] if molecules is not None else participant


def _mapped_graphs(participants: Sequence[Any]) -> Tuple[MappedGraphV1, ...]:
    return tuple(
        MappedGraphV1(
            participant_index=index,
            adjacency_list=_molecule(participant).to_adjacency_list(remove_h=False),
            atom_ids=tuple(atom.id for atom in _molecule(participant).atoms),
        )
        for index, participant in enumerate(participants)
    )


def _side_atoms(participants: Sequence[Any]) -> dict:
    result = {}
    for participant_index, participant in enumerate(participants):
        for atom in _molecule(participant).atoms:
            if atom.id in result:
                raise MergeAccountingError(
                    "reaction side contains duplicate atom IDs; rooted mapping is ambiguous"
                )
            result[atom.id] = (participant_index, atom)
    return result


def _reaction_center(reaction) -> Tuple[RootAtomV1, ...]:
    reactants = _side_atoms(reaction.reactants)
    products = _side_atoms(reaction.products)
    if set(reactants) != set(products):
        raise MergeAccountingError(
            "reaction does not preserve a bijective atom-ID mapping"
        )

    def bonds(participants):
        result = {}
        for participant in participants:
            molecule = _molecule(participant)
            for bond in molecule.get_all_edges():
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
    attributes = (
        "radical_electrons",
        "charge",
        "lone_pairs",
    )
    for atom_id, (_, reactant_atom) in reactants.items():
        product_atom = products[atom_id][1]
        if any(
            getattr(reactant_atom, attribute, 0)
            != getattr(product_atom, attribute, 0)
            for attribute in attributes
        ):
            reacting.add(atom_id)

    if not reacting:
        raise MergeAccountingError("reaction path has no mapped reaction centre")
    return tuple(
        RootAtomV1(reactants[atom_id][0], atom_id)
        for atom_id in sorted(reacting)
    )


def _group_raw_reactions(reactions: Sequence[Any]) -> list:
    """Mirror the grouping half of ``find_degenerate_reactions``.

    Exact duplicate paths remain visible as group members but are marked as
    non-contributors, matching RMG's line 348-350 accounting.
    """

    groups = []
    for reaction0 in reactions:
        reaction0.ensure_species(save_order=False)
        if not groups:
            groups.append(_MergeGroup([reaction0], [reaction0]))
            continue
        for group in groups:
            isomorphic = False
            identical = False
            same_template = True
            for reaction in group.contributors:
                isomorphic = reaction0.is_isomorphic(
                    reaction,
                    check_identical=False,
                    strict=False,
                    check_template_rxn_products=True,
                    save_order=False,
                )
                if isomorphic:
                    identical = reaction0.is_isomorphic(
                        reaction,
                        check_identical=True,
                        strict=False,
                        check_template_rxn_products=True,
                        save_order=False,
                    )
                    if identical:
                        break
                    same_template = frozenset(reaction.template) == frozenset(
                        reaction0.template
                    )
                else:
                    break
            if identical:
                group.members.append(reaction0)
                break
            if isomorphic and same_template:
                group.contributors.append(reaction0)
                group.members.append(reaction0)
                break
            if isomorphic:
                reaction0.duplicate = True
                group.contributors[0].duplicate = True
                continue
        else:
            groups.append(_MergeGroup([reaction0], [reaction0]))
    return groups


def _path_snapshots(reactions: Sequence[Any], prefix="raw") -> dict:
    return {
        id(reaction): RootedPathV1(
            path_id="{0}-{1:06d}".format(prefix, index),
            degeneracy=float(reaction.degeneracy),
            root=_reaction_center(reaction),
        )
        for index, reaction in enumerate(reactions)
    }


def _corrected_divisor(
    reaction, pre_correction_degeneracy: float, final_degeneracy: float
) -> float:
    if not reaction.is_forward:
        return 1.0
    if final_degeneracy == 0:
        raise MergeAccountingError("final degeneracy is zero")
    divisor = pre_correction_degeneracy / final_degeneracy
    if divisor not in (1.0, 2.0, 6.0):
        raise MergeAccountingError(
            "RMG applied an unsupported identical-reactant divisor {0!r}".format(
                divisor
            )
        )
    return divisor


def _root_weights(
    paths: Iterable[RootedPathV1], divisor: float, final_degeneracy: float
) -> Tuple[RootWeightV1, ...]:
    grouped = OrderedDict()
    for path in paths:
        if not path.contributes_to_pre_correction_sum:
            continue
        item = grouped.setdefault(path.root, [0.0, []])
        item[0] += path.degeneracy / divisor
        item[1].append(path.path_id)
    if not grouped:
        raise MergeAccountingError("no rooted paths contribute to final degeneracy")

    entries = list(grouped.items())
    weights = []
    allocated = 0.0
    for index, (root, (weight, path_ids)) in enumerate(entries):
        if index == len(entries) - 1:
            weight = final_degeneracy - allocated
        allocated += weight
        weights.append(RootWeightV1(root, weight, tuple(path_ids)))
    if sum(weight.weight for weight in weights) != final_degeneracy:
        raise MergeAccountingError(
            "root weights do not sum exactly to final degeneracy"
        )
    return tuple(weights)


def _trace_recalculated_degeneracy(family, reaction, resonance):
    """Re-run ``calculate_degeneracy``'s isolated generation for evidence."""

    traced_reaction = copy.deepcopy(reaction)
    reactants, same_reactants = check_for_same_reactants(
        traced_reaction.reactants
    )
    ensure_independent_atom_ids(reactants, resonance=resonance)
    if any(
        participant.is_polymer_proxy
        for participant in traced_reaction.reactants + traced_reaction.products
    ):
        for participant in reactants:
            participant.is_polymer_proxy = True

    raw_reactions = []
    for combo in generate_molecule_combos(reactants):
        raw_reactions.extend(
            family._generate_reactions(
                combo,
                products=traced_reaction.products,
                forward=True,
                prod_resonance=resonance,
                react_non_reactive=True,
                apply_species_constraints=False,
            )
        )
    template = frozenset(traced_reaction.template)
    selected = [
        candidate
        for candidate in raw_reactions
        if frozenset(candidate.template) == template
    ]
    if not selected:
        raise MergeAccountingError(
            "reverse degeneracy recalculation produced no paths for the returned template"
        )

    groups = _group_raw_reactions(selected)
    snapshots = _path_snapshots(selected, prefix="recalculated")
    group_by_representative = {
        id(group.contributors[0]): group for group in groups
    }
    merged = find_degenerate_reactions(
        raw_reactions,
        same_reactants,
        template=traced_reaction.template,
        kinetics_family=family,
        resonance=resonance,
    )
    if len(merged) != 1:
        raise MergeAccountingError(
            "reverse degeneracy recalculation returned {0} merged channels, expected one".format(
                len(merged)
            )
        )
    representative = merged[0]
    if float(representative.degeneracy) != float(reaction.degeneracy):
        raise MergeAccountingError(
            "recalculated trace degeneracy {0!r} does not match G_final {1!r}".format(
                representative.degeneracy, reaction.degeneracy
            )
        )
    group = group_by_representative[id(representative)]
    contributors = set(map(id, group.contributors))
    paths = tuple(
        RootedPathV1(
            path_id=snapshots[id(path)].path_id,
            degeneracy=snapshots[id(path)].degeneracy,
            root=snapshots[id(path)].root,
            contributes_to_pre_correction_sum=id(path) in contributors,
        )
        for path in group.members
    )
    contributing = tuple(
        path for path in paths if path.contributes_to_pre_correction_sum
    )
    pre_correction = sum(path.degeneracy for path in contributing)
    divisor = _corrected_divisor(
        representative, pre_correction, float(representative.degeneracy)
    )
    return paths, contributing, int(same_reactants), divisor


def _generate_traced(
    kinetics_database,
    reactants,
    products=None,
    only_families=None,
    resonance=True,
) -> MergeAccountingResult:
    if isinstance(reactants, tuple):
        reactants = list(reactants)
    reactants, same_reactants = check_for_same_reactants(reactants)
    ensure_independent_atom_ids(reactants, resonance=resonance)

    raw_reactions = []
    for combo in generate_molecule_combos(reactants):
        raw_reactions.extend(
            kinetics_database.react_molecules(
                combo,
                products=products,
                only_families=only_families,
                prod_resonance=resonance,
            )
        )

    groups = _group_raw_reactions(raw_reactions)
    snapshots = _path_snapshots(raw_reactions)
    group_by_representative = {id(group.contributors[0]): group for group in groups}

    reactions = find_degenerate_reactions(
        raw_reactions,
        same_reactants,
        kinetics_database=kinetics_database,
        resonance=resonance,
    )

    to_delete = []
    for index, reaction in enumerate(reactions):
        family = kinetics_database.families[reaction.family]
        if family.own_reverse and not family.add_reverse_attribute(reaction):
            to_delete.append(index)
    for index in reversed(to_delete):
        del reactions[index]

    records = []
    for reaction in reactions:
        family = kinetics_database.families[reaction.family]
        group = group_by_representative[id(reaction)]
        contributors = set(map(id, group.contributors))
        raw_paths = tuple(
            RootedPathV1(
                path_id=snapshots[id(path)].path_id,
                degeneracy=snapshots[id(path)].degeneracy,
                root=snapshots[id(path)].root,
                contributes_to_pre_correction_sum=id(path) in contributors,
            )
            for path in group.members
        )
        contributing_paths = tuple(
            path for path in raw_paths if path.contributes_to_pre_correction_sum
        )
        pre_correction = sum(path.degeneracy for path in contributing_paths)
        final_degeneracy = float(reaction.degeneracy)
        reverse_recalculated = not reaction.is_forward and not family.own_reverse
        if reverse_recalculated:
            (
                recalculated_paths,
                allocation_paths,
                accounting_same_reactants,
                divisor,
            ) = _trace_recalculated_degeneracy(
                family, reaction, resonance
            )
        else:
            recalculated_paths = ()
            allocation_paths = contributing_paths
            accounting_same_reactants = int(same_reactants)
            divisor = _corrected_divisor(
                reaction, pre_correction, final_degeneracy
            )
        records.append(
            RMGMergeAccountingV1(
                schema="RMGMergeAccountingV1",
                family=reaction.family,
                direction="forward" if reaction.is_forward else "reverse",
                family_own_reverse=bool(family.own_reverse),
                reactant_graphs=_mapped_graphs(reaction.reactants),
                product_graphs=_mapped_graphs(reaction.products),
                raw_paths=raw_paths,
                merge_sublist_path_ids=tuple(path.path_id for path in raw_paths),
                contributing_path_ids=tuple(
                    path.path_id for path in contributing_paths
                ),
                pre_correction_degeneracy=pre_correction,
                generation_same_reactants=int(same_reactants),
                same_reactants=accounting_same_reactants,
                full_species_reactants_identical=accounting_same_reactants
                not in (0, 1),
                identical_reactant_divisor=divisor,
                reverse_degeneracy_recalculated=reverse_recalculated,
                recalculated_rooted_paths=recalculated_paths,
                G_final=final_degeneracy,
                root_weights=_root_weights(
                    allocation_paths, divisor, final_degeneracy
                ),
            )
        )
    return MergeAccountingResult(reactions, tuple(records))


def generate_reactions_with_merge_accounting(
    kinetics_database,
    reactants,
    products=None,
    only_families=None,
    resonance=True,
    trace=True,
) -> MergeAccountingResult:
    """Generate family reactions with an optional immutable accounting trace.

    ``trace=False`` delegates directly to RMG.  ``trace=True`` drives the same
    isolated sequence while retaining evidence before RMG mutates degeneracy.
    """

    if not trace:
        reactions = kinetics_database.generate_reactions_from_families(
            reactants,
            products=products,
            only_families=only_families,
            resonance=resonance,
        )
        return MergeAccountingResult(reactions, ())
    return _generate_traced(
        kinetics_database,
        reactants,
        products=products,
        only_families=only_families,
        resonance=resonance,
    )


__all__ = [
    "MappedGraphV1",
    "MergeAccountingError",
    "MergeAccountingResult",
    "RMGMergeAccountingV1",
    "RootAtomV1",
    "RootedPathV1",
    "RootWeightV1",
    "generate_reactions_with_merge_accounting",
]
