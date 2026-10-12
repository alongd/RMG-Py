"""Observational accounting for RMG reaction-path merging.

The traced backend drives RMG's family generation, grouping, correction, and
base-family reverse-degeneracy algorithm while retaining discarded evidence.

Tracing performs extra generation. Family objects, optional profiling counters,
and the process-wide atom-ID allocator are shared mutable RMG state. Traced
generation is therefore serialized here and must not share a backend with any
other thread. Cold ON/OFF parity tests use separately loaded databases.
"""

from collections import OrderedDict
from contextlib import contextmanager
from dataclasses import dataclass
from fractions import Fraction
import itertools
import math
import threading
from typing import Any, Iterable, Sequence, Tuple

from rmgpy.data.kinetics.common import (
    check_for_same_reactants,
    ensure_independent_atom_ids,
    generate_molecule_combos,
    reduce_same_reactant_degeneracy,
)
from rmgpy.data.kinetics.family import KineticsFamily
from rmgpy.exceptions import KineticsError


@dataclass(frozen=True)
class MappedGraphV2:
    participant_index: int
    adjacency_list: str
    atom_ids: Tuple[int, ...]


@dataclass(frozen=True, order=True)
class RootBindingV2:
    role_label: str
    participant_index: int
    atom_id: int


@dataclass(frozen=True)
class RecipeEditV2:
    action: str
    arguments: Tuple[Any, ...]


@dataclass(frozen=True, order=True)
class PathAtomMapV2:
    raw_participant_index: int
    raw_atom_index: int
    raw_atom_id: int
    serialized_participant_index: int
    serialized_atom_index: int
    serialized_atom_id: int
    raw_element: str
    raw_isotope: int
    serialized_element: str
    serialized_isotope: int


@dataclass(frozen=True)
class RootedPathV2:
    path_id: str
    degeneracy: Any
    root: Tuple[RootBindingV2, ...]
    recipe: Tuple[RecipeEditV2, ...]
    atom_map: Tuple[PathAtomMapV2, ...] = ()
    contributes_to_pre_correction_sum: bool = True


@dataclass(frozen=True)
class RootWeightV2:
    root: Tuple[RootBindingV2, ...]
    recipe: Tuple[RecipeEditV2, ...]
    weight: Fraction
    contributing_path_ids: Tuple[str, ...]


@dataclass(frozen=True)
class RMGMergeAccountingV2:
    schema: str
    family: str
    direction: str
    family_own_reverse: bool
    reactant_graphs: Tuple[MappedGraphV2, ...]
    product_graphs: Tuple[MappedGraphV2, ...]
    raw_paths: Tuple[RootedPathV2, ...]
    merge_sublist_path_ids: Tuple[str, ...]
    contributing_path_ids: Tuple[str, ...]
    pre_correction_degeneracy: Fraction
    allocation_pre_correction_degeneracy: Fraction
    generation_same_reactants: int
    same_reactants: int
    full_species_reactants_identical: bool
    identical_reactant_divisor: int
    reverse_degeneracy_recalculated: bool
    recalculated_rooted_paths: Tuple[RootedPathV2, ...]
    G_final: Fraction
    root_weights: Tuple[RootWeightV2, ...]


@dataclass(frozen=True)
class MergeAccountingResult:
    reactions: list
    records: Tuple[RMGMergeAccountingV2, ...]


class MergeAccountingError(ValueError):
    """RMG's final degeneracy cannot be rooted without guessing."""


@dataclass
class _MergeGroup:
    contributor_indices: list
    member_indices: list


@dataclass(frozen=True)
class _PathSnapshot:
    path_id: str
    degeneracy: Fraction
    root: Tuple[RootBindingV2, ...]
    recipe: Tuple[RecipeEditV2, ...]


@dataclass
class _MergedEvidence:
    reaction: Any
    group: _MergeGroup
    raw_pre_correction: Fraction
    recalculated_paths: Tuple[RootedPathV2, ...]
    allocation_pre_correction: Fraction
    accounting_same_reactants: int
    divisor: int


_trace_lock = threading.Lock()


@contextmanager
def _exclusive_trace():
    if not _trace_lock.acquire(blocking=False):
        raise MergeAccountingError(
            "traced generation shares process-wide RMG state and may not "
            "run concurrently"
        )
    try:
        yield
    finally:
        _trace_lock.release()


def _molecule(participant):
    molecules = getattr(participant, "molecule", None)
    return molecules[0] if molecules is not None else participant


def _mapped_graphs(participants: Sequence[Any]) -> Tuple[MappedGraphV2, ...]:
    return tuple(
        MappedGraphV2(
            index,
            _molecule(participant).to_adjacency_list(remove_h=False),
            tuple(atom.id for atom in _molecule(participant).atoms),
        )
        for index, participant in enumerate(participants)
    )


def _positive_fraction(value, field):
    try:
        scalar = float(value)
    except (TypeError, ValueError, OverflowError):
        raise MergeAccountingError(
            "{0} is not a finite positive number".format(field)
        )
    if not math.isfinite(scalar) or scalar <= 0:
        raise MergeAccountingError(
            "{0} is not a finite positive number".format(field)
        )
    return Fraction(str(value))


def _recipe_for(family, is_forward):
    recipe = family.forward_recipe if is_forward else family.reverse_recipe
    if recipe is None:
        raise MergeAccountingError(
            "family {0!r} has no recipe for the generated direction".format(
                family.label
            )
        )
    return tuple(
        RecipeEditV2(str(action[0]), tuple(action[1:]))
        for action in recipe.actions
    )


def _root_bindings(reaction, recipe):
    labeled = getattr(reaction, "labeled_atoms", None)
    if not labeled or "reactants" not in labeled:
        raise MergeAccountingError(
            "raw reaction did not retain recipe-role bindings"
        )
    required = {
        argument
        for edit in recipe
        for argument in edit.arguments
        if isinstance(argument, str) and argument.startswith("*")
    }
    bindings = []
    present = set()
    for label, value in labeled["reactants"].items():
        for atom in value if isinstance(value, list) else [value]:
            locations = [
                i
                for i, participant in enumerate(reaction.reactants)
                if any(
                    candidate is atom
                    for candidate in _molecule(participant).atoms
                )
            ]
            if len(locations) != 1:
                raise MergeAccountingError(
                    "recipe role {0!r} does not identify one raw reactant "
                    "atom".format(label)
                )
            bindings.append(RootBindingV2(label, locations[0], atom.id))
            present.add(label)
    missing = sorted(required - present)
    if missing:
        raise MergeAccountingError(
            "raw reaction is missing recipe roles {0!r}".format(missing)
        )
    if not bindings:
        raise MergeAccountingError("reaction path has no recipe-role bindings")
    return tuple(sorted(bindings))


def _snapshot_reaction(reaction, family, path_id):
    recipe = _recipe_for(family, reaction.is_forward)
    return _PathSnapshot(
        path_id,
        _positive_fraction(reaction.degeneracy, "raw path degeneracy"),
        _root_bindings(reaction, recipe),
        recipe,
    )


def _clear_generation_labels(reactions, input_molecules):
    for molecule in input_molecules:
        molecule.clear_labeled_atoms()
    for reaction in reactions:
        for participant in itertools.chain(
            reaction.reactants, reaction.products
        ):
            participant.clear_labeled_atoms()
        if hasattr(reaction, "labeled_atoms"):
            del reaction.labeled_atoms


def _assert_no_reaction_aliases(reactions):
    seen = set()
    for index, reaction in enumerate(reactions):
        marker = id(reaction)
        if marker in seen:
            raise MergeAccountingError(
                "raw reaction occurrence {0} aliases an earlier occurrence; "
                "refusing ambiguous accounting".format(index)
            )
        seen.add(marker)


def _generate_database_paths(
    database, combos, products, only_families, resonance
):
    reactions, snapshots = [], []
    for combo in combos:
        combo_reactions = []
        for label, family in database.families.items():
            if only_families is not None and label not in only_families:
                continue
            generated = family.generate_reactions(
                combo,
                products=products,
                prod_resonance=resonance,
                delete_labels=False,
            )
            for reaction in generated:
                snapshots.append(
                    _snapshot_reaction(
                        reaction, family, "raw-{0:06d}".format(len(reactions))
                    )
                )
                reactions.append(reaction)
                combo_reactions.append(reaction)
        _clear_generation_labels(combo_reactions, combo)
    _assert_no_reaction_aliases(reactions)
    return reactions, snapshots


def _generate_family_paths(family, combos, products, resonance):
    reactions, snapshots = [], []
    for combo in combos:
        generated = family._generate_reactions(
            combo,
            products=products,
            forward=True,
            prod_resonance=resonance,
            react_non_reactive=True,
            delete_labels=False,
            apply_species_constraints=False,
        )
        for reaction in generated:
            snapshots.append(
                _snapshot_reaction(
                    reaction,
                    family,
                    "recalculated-{0:06d}".format(len(reactions)),
                )
            )
            reactions.append(reaction)
        _clear_generation_labels(generated, combo)
    _assert_no_reaction_aliases(reactions)
    return reactions, snapshots


def _group_raw_reactions(reactions: Sequence[Any]) -> list:
    """Mirror RMG grouping while tracking raw occurrences by position."""
    groups = []
    for index, reaction0 in enumerate(reactions):
        reaction0.ensure_species(save_order=False)
        if not groups:
            groups.append(_MergeGroup([index], [index]))
            continue
        for group in groups:
            isomorphic = identical = False
            same_template = True
            for contributor_index in group.contributor_indices:
                reaction = reactions[contributor_index]
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
                group.member_indices.append(index)
                break
            if isomorphic and same_template:
                group.contributor_indices.append(index)
                group.member_indices.append(index)
                break
            if isomorphic:
                reaction0.duplicate = True
                reactions[group.contributor_indices[0]].duplicate = True
                continue
        else:
            groups.append(_MergeGroup([index], [index]))
    return groups


def _derived_divisor(reaction, same_reactants):
    if not reaction.is_forward or same_reactants in (0, 1):
        return 1
    if len(reaction.reactants) == 2:
        return 2 if same_reactants == 2 else 1
    if len(reaction.reactants) == 3:
        return 6 if same_reactants == 3 else 2 if same_reactants == 2 else 1
    return 1


def _corrected_divisor(
    reaction, same_reactants, pre_correction_degeneracy, final_degeneracy
):
    """Derive RMG's correction branch, then verify its observed result."""
    pre = _positive_fraction(
        pre_correction_degeneracy, "pre-correction degeneracy"
    )
    final = float(_positive_fraction(final_degeneracy, "final degeneracy"))
    divisor = _derived_divisor(reaction, same_reactants)
    expected = pre / divisor
    if float(expected) != final:
        raise MergeAccountingError(
            "derived divisor {0} predicts {1!r}, but RMG produced "
            "{2!r}".format(divisor, float(expected), final)
        )
    return divisor


def _path_atom_mapping(path_reaction, representative):
    targets = {}
    for participant_index, participant in enumerate(representative.reactants):
        for atom_index, atom in enumerate(_molecule(participant).atoms):
            if atom.id in targets:
                raise MergeAccountingError(
                    "serialized reactant graphs contain duplicate atom IDs"
                )
            targets[atom.id] = (participant_index, atom_index, atom)
    result = []
    for raw_participant_index, participant in enumerate(
        path_reaction.reactants
    ):
        for raw_atom_index, atom in enumerate(_molecule(participant).atoms):
            if atom.id not in targets:
                raise MergeAccountingError(
                    "raw path atom {0} is absent from serialized reactant "
                    "graphs".format(atom.id)
                )
            target_participant_index, target_atom_index, target = targets[
                atom.id
            ]
            raw_identity = (atom.element.symbol, int(atom.element.isotope))
            target_identity = (
                target.element.symbol,
                int(target.element.isotope),
            )
            if raw_identity != target_identity:
                raise MergeAccountingError(
                    "raw-to-serialized mapping changes element or isotope"
                )
            result.append(
                PathAtomMapV2(
                    raw_participant_index,
                    raw_atom_index,
                    atom.id,
                    target_participant_index,
                    target_atom_index,
                    target.id,
                    raw_identity[0],
                    raw_identity[1],
                    target_identity[0],
                    target_identity[1],
                )
            )
    if len(result) != len(targets):
        raise MergeAccountingError(
            "raw path does not map bijectively into serialized reactants"
        )
    return tuple(result)


def _serialize_group_paths(group, reactions, snapshots, representative):
    contributors = set(group.contributor_indices)
    paths = []
    for index in group.member_indices:
        snapshot = snapshots[index]
        atom_map = _path_atom_mapping(reactions[index], representative)
        targets = {
            (item.raw_participant_index, item.raw_atom_id): (
                item.serialized_participant_index,
                item.serialized_atom_id,
            )
            for item in atom_map
        }
        root = tuple(
            sorted(
                RootBindingV2(
                    binding.role_label,
                    *targets[(binding.participant_index, binding.atom_id)]
                )
                for binding in snapshot.root
            )
        )
        paths.append(
            RootedPathV2(
                snapshot.path_id,
                snapshot.degeneracy,
                root,
                snapshot.recipe,
                atom_map,
                index in contributors,
            )
        )
    return tuple(paths)


def _root_weights(
    paths: Iterable[RootedPathV2], divisor: int, final_degeneracy
):
    divisor_fraction = _positive_fraction(
        divisor, "identical-reactant divisor"
    )
    final = _positive_fraction(final_degeneracy, "G_final")
    grouped = OrderedDict()
    for path in paths:
        if not path.contributes_to_pre_correction_sum:
            continue
        item = grouped.setdefault((path.root, path.recipe), [Fraction(0), []])
        item[0] += _positive_fraction(path.degeneracy, "raw path degeneracy")
        item[1].append(path.path_id)
    if not grouped:
        raise MergeAccountingError(
            "no rooted paths contribute to final degeneracy"
        )
    natural_total = (
        sum((item[0] for item in grouped.values()), Fraction(0))
        / divisor_fraction
    )
    if natural_total != final:
        raise MergeAccountingError(
            "rooted-path allocation {0!r} does not match G_final {1!r}".format(
                natural_total, final
            )
        )
    weights = tuple(
        RootWeightV2(root, recipe, total / divisor_fraction, tuple(path_ids))
        for (root, recipe), (total, path_ids) in grouped.items()
    )
    if sum((weight.weight for weight in weights), Fraction(0)) != final:
        raise MergeAccountingError(
            "root weights do not sum exactly to final degeneracy"
        )
    return weights


def _merge_forward_groups(reactions, snapshots, groups, same_reactants):
    merged = []
    for group in groups:
        representative = reactions[group.contributor_indices[0]]
        representative.degeneracy = sum(
            reactions[i].degeneracy for i in group.contributor_indices
        )
        pre = sum(
            (snapshots[i].degeneracy for i in group.contributor_indices),
            Fraction(0),
        )
        reduce_same_reactant_degeneracy(representative, same_reactants)
        divisor = _corrected_divisor(
            representative, same_reactants, pre, representative.degeneracy
        )
        merged.append((representative, group, pre, divisor))
    return merged


def _trace_base_degeneracy(family, reaction, resonance):
    """Capture roots in the base calculation that returns the scalar."""
    if reaction.is_charge_transfer_reaction():
        raise MergeAccountingError(
            "charge-transfer calculate_degeneracy returns 1 without a rooted "
            "decomposition"
        )
    if (
        type(family).calculate_degeneracy
        is not KineticsFamily.calculate_degeneracy
    ):
        raise MergeAccountingError(
            "family calculate_degeneracy override cannot be rooted without "
            "an evidence hook"
        )
    reactants, same_reactants = check_for_same_reactants(reaction.reactants)
    reaction.reactants = reactants
    ensure_independent_atom_ids(reactants, resonance=resonance)
    if any(p.is_polymer_proxy for p in reaction.reactants + reaction.products):
        for participant in reactants:
            participant.is_polymer_proxy = True
    raw, snapshots = _generate_family_paths(
        family,
        generate_molecule_combos(reactants),
        reaction.products,
        resonance,
    )
    template = frozenset(reaction.template)
    selected = [
        i
        for i, candidate in enumerate(raw)
        if frozenset(candidate.template) == template
    ]
    reactions = [raw[i] for i in selected]
    selected_snapshots = [snapshots[i] for i in selected]
    merged = (
        _merge_forward_groups(
            reactions,
            selected_snapshots,
            _group_raw_reactions(reactions),
            same_reactants,
        )
        if reactions
        else []
    )
    if len(merged) != 1:
        raise KineticsError(
            "Unable to calculate degeneracy for reaction {0} in reaction "
            "family {1}. "
            "Expected 1 reaction but generated {2}".format(
                reaction, family.label, len(merged)
            )
        )
    recalculated, group, pre, divisor = merged[0]
    paths = _serialize_group_paths(
        group, reactions, selected_snapshots, reaction
    )
    return recalculated.degeneracy, paths, int(same_reactants), pre, divisor


def _generate_traced(
    database, reactants, products=None, only_families=None, resonance=True
):
    if isinstance(reactants, tuple):
        reactants = list(reactants)
    reactants, same_reactants = check_for_same_reactants(reactants)
    ensure_independent_atom_ids(reactants, resonance=resonance)
    raw, snapshots = _generate_database_paths(
        database,
        generate_molecule_combos(reactants),
        products,
        only_families,
        resonance,
    )
    groups = _group_raw_reactions(raw)
    merged = []
    for group in groups:
        reaction = raw[group.contributor_indices[0]]
        reaction.degeneracy = sum(
            raw[i].degeneracy for i in group.contributor_indices
        )
        raw_pre = sum(
            (snapshots[i].degeneracy for i in group.contributor_indices),
            Fraction(0),
        )
        family = database.families[reaction.family]
        if reaction.is_forward:
            reduce_same_reactant_degeneracy(reaction, same_reactants)
            divisor = _corrected_divisor(
                reaction, same_reactants, raw_pre, reaction.degeneracy
            )
            recalculated, allocation_pre, accounting_same = (
                (),
                raw_pre,
                int(same_reactants),
            )
        elif family.own_reverse:
            divisor = _corrected_divisor(
                reaction, same_reactants, raw_pre, reaction.degeneracy
            )
            recalculated, allocation_pre, accounting_same = (
                (),
                raw_pre,
                int(same_reactants),
            )
        else:
            (
                reaction.degeneracy,
                recalculated,
                accounting_same,
                allocation_pre,
                divisor,
            ) = _trace_base_degeneracy(family, reaction, resonance)
        merged.append(
            _MergedEvidence(
                reaction,
                group,
                raw_pre,
                recalculated,
                allocation_pre,
                accounting_same,
                divisor,
            )
        )

    reactions = [item.reaction for item in merged]
    to_delete = []
    for index, reaction in enumerate(reactions):
        family = database.families[reaction.family]
        if family.own_reverse and not family.add_reverse_attribute(reaction):
            to_delete.append(index)
    for index in reversed(to_delete):
        del reactions[index]
        del merged[index]

    records = []
    for evidence in merged:
        reaction = evidence.reaction
        group = evidence.group
        raw_pre = evidence.raw_pre_correction
        recalculated = evidence.recalculated_paths
        allocation_pre = evidence.allocation_pre_correction
        accounting_same = evidence.accounting_same_reactants
        divisor = evidence.divisor
        family = database.families[reaction.family]
        raw_paths = _serialize_group_paths(group, raw, snapshots, reaction)
        contributing = tuple(
            p for p in raw_paths if p.contributes_to_pre_correction_sum
        )
        reverse_recalculated = (
            not reaction.is_forward and not family.own_reverse
        )
        allocation_paths = (
            recalculated if reverse_recalculated else contributing
        )
        final = allocation_pre / divisor
        if float(final) != float(reaction.degeneracy):
            raise MergeAccountingError(
                "captured degeneracy does not match returned reaction "
                "degeneracy"
            )
        records.append(
            RMGMergeAccountingV2(
                "RMGMergeAccountingV2",
                reaction.family,
                "forward" if reaction.is_forward else "reverse",
                bool(family.own_reverse),
                _mapped_graphs(reaction.reactants),
                _mapped_graphs(reaction.products),
                raw_paths,
                tuple(p.path_id for p in raw_paths),
                tuple(p.path_id for p in contributing),
                raw_pre,
                allocation_pre,
                int(same_reactants),
                accounting_same,
                accounting_same not in (0, 1),
                divisor,
                reverse_recalculated,
                recalculated,
                final,
                _root_weights(allocation_paths, divisor, final),
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
):
    """Generate reactions and, when requested, immutable merge evidence."""
    if not trace:
        reactions = kinetics_database.generate_reactions_from_families(
            reactants,
            products=products,
            only_families=only_families,
            resonance=resonance,
        )
        return MergeAccountingResult(reactions, ())
    with _exclusive_trace():
        return _generate_traced(
            kinetics_database,
            reactants,
            products=products,
            only_families=only_families,
            resonance=resonance,
        )


__all__ = [
    "MappedGraphV2",
    "MergeAccountingError",
    "MergeAccountingResult",
    "PathAtomMapV2",
    "RecipeEditV2",
    "RMGMergeAccountingV2",
    "RootBindingV2",
    "RootedPathV2",
    "RootWeightV2",
    "generate_reactions_with_merge_accounting",
]
