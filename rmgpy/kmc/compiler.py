"""Offline compilation of RMG family reactions into kMC rewrite records.

The compiler asks RMG's public ``generate_reactions_from_families`` API for
merged reaction paths and serializes their graphs and kinetics as a stable,
content-addressed JSON artifact.  Owner-approved J-ring capture channels use
frozen reactions and kinetics reconstructed from their archived RMG runs.
"""

from __future__ import annotations

import copy
from contextvars import ContextVar
import hashlib
import json
import logging
import math
import os
import subprocess
from dataclasses import dataclass, field, replace
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable, Sequence
from types import SimpleNamespace

from rmgpy.kmc.atom_map import extract_atom_map
from rmgpy.kmc.event_record import EventRecord, ssa_multiplier_for
from rmgpy.kmc.kinetics_library import (
    load_plpsec_entry, matches_head_to_tail, plpsec_rate_table,
)
from rmgpy.kmc.reference_thermo import (
    FROZEN_THERMO_PROPERTY,
    GasPhaseRMGReferenceThermo,
    ReferenceThermoProvider,
    SharedThermoAssignment,
    ThermoUnavailable,
)
from rmgpy.kinetics.arrhenius import ArrheniusBM, ArrheniusEP
from rmgpy.kinetics.model import get_rate_coefficient_units_from_reaction_order
from rmgpy.exceptions import ActionError


_PAIR_MOLECULE_CACHE = ContextVar("pair_molecule_cache", default=None)
SCHEMA_VERSION = "kmc_event_set/0.1"
DEFAULT_T_GRID = tuple(float(temperature) for temperature in range(300, 1001, 25))
ATOMIC_MASS_NUMBERS = {"H": 1, "C": 12, "N": 14, "O": 16, "S": 32}
PS_PROXY_UNITS = 3
PS_FAMILY_CANDIDATES = (
    "Disproportionation",
    "H_Abstraction",
    "R_Addition_MultipleBond",
    "R_Recombination",
    "intra_H_migration",
)
PS_FAMILY_FILTER_REASON = (
    "outside the documented PS radical-growth family filter "
    "(H transfer, radical addition/recombination/disproportionation, and migration)"
)
LINK_PLACEHOLDER = "evt_" + "0" * 64
PERSISTENT_CARBENE_POLICY_VERSION = "persistent-carbene-applicability/1"
R1_JUNCTION_LABELS = {
    "J_para": {
        "attacked_atom_label": "S9",
        "chosen_resonance_localisation": "para:S9",
        "reflection": None,
    },
    "J_ortho_S7": {
        "attacked_atom_label": "S7",
        "chosen_resonance_localisation": "ortho:S7",
        "reflection": None,
    },
    "J_ortho_S6": {
        "attacked_atom_label": "S6",
        "chosen_resonance_localisation": "ortho:S6-reflection",
        "reflection": {"S6": "S7", "S7": "S6", "S8": "S10", "S10": "S8"},
    },
}
ARCHIVED_J_PARA_RULE = (
    "Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_1C-inRing_"
    "Ext-1C-R_Sp-3R!H-1C_N-2R-inRing_Ext-2R-R"
)
ARCHIVED_J_PARA_RATE_PROVENANCE = {
    "channel_id": "J_para",
    "reaction_label": "para J_ring",
    "family": "R_Recombination",
    "rule_entry_index": 176,
    "rule_entry_label": ARCHIVED_J_PARA_RULE,
    "rule_rank": 11,
    "direction": "forward",
    "degeneracy": 1.0,
    "arrhenius": {
        "A_m3_mol_s": 1.76793e10,
        "n": -1.00291,
        "Ea_J_mol": 0.0,
        "T0_K": 1.0,
    },
    "database_sha": "cd86d4e1c187a132109e16cd86f624ed9fb217df",
    "extraction_script": "polymer-pm/reports/R-009_sources/R-009_v8_rmg_archive.py",
    "archive_record": "polymer-pm/reports/R-009_sources/R-009_v8_rmg_archive.txt",
}


def canonical_json_bytes(value: Any) -> bytes:
    """Return the one canonical representation used for artifact identity."""
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def compiler_source_hash() -> str:
    """Fingerprint the code loaded by this interpreter, not later file edits."""
    return _LOADED_SOURCE_HASH


_LOADED_SOURCE_HASH = hashlib.sha256(
    b"".join(
        (Path(__file__).parent / filename).read_bytes()
        for filename in (
            "compiler.py",
            "reference_thermo.py",
            "event_record.py",
            "atom_map.py",
            "kinetics_library.py",
        )
    )
).hexdigest()
_LOADED_COMPILER_HASH = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def _git_sha(path: str | Path | None) -> str | None:
    if not path:
        return None
    try:
        return subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def prepare_rate_rules(kinetics_database, thermo_database, *,
                       kinetics_depositories=("training",), verbose=True):
    """Apply RMG.load_database's non-ATG rate-rule preparation once per family.

    Training reactions must not inherit model-growth species/polymer limits.
    Keep the caller's global input and constraints intact, including on failure.
    Existing family methods perform all rule construction and averaging.
    """
    from rmgpy.rmg import input as rmg_input

    families = getattr(kinetics_database, "families", {})
    add_training = "!training" not in kinetics_depositories
    policy = {"training_requested": add_training, "verbose": verbose}
    prepared = {}
    original_input = rmg_input.rmg
    owner = original_input or SimpleNamespace(
        species_constraints={}, polymer_constraints=None, quantum_mechanics=None,
        ml_estimator=None, ml_settings=None, thermo_central_database=None,
    )
    species_constraints = owner.species_constraints
    polymer_constraints = owner.polymer_constraints
    try:
        rmg_input.rmg = owner
        owner.species_constraints = {}
        owner.polymer_constraints = None
        for label, family in sorted(families.items()):
            # Lightweight compiler fixtures may have no RMG family objects.
            if not hasattr(family, "auto_generated") or family.auto_generated:
                continue
            previous = getattr(family, "_kmc_rate_rule_preparation", None)
            if previous is not None:
                if previous["policy"] != policy:
                    raise ValueError("rate-rule preparation policy changed on a loaded family")
                prepared[label] = copy.deepcopy(previous)
                continue
            if add_training and thermo_database is None:
                raise ValueError("training rate rules require a thermo database")
            count = lambda: sum(len(entries) for entries in family.rules.entries.values())
            before = count()
            if add_training:
                logging.getLogger(__name__).info("adding training rate rules: %s", label)
                family.add_rules_from_training(thermo_database=thermo_database)
            after_training = count()
            # RMG restores constraints before averaging the resulting rules.
            owner.species_constraints = species_constraints
            owner.polymer_constraints = polymer_constraints
            family.fill_rules_by_averaging_up(verbose=verbose)
            owner.species_constraints = {}
            owner.polymer_constraints = None
            result = {
                "policy": policy, "training_rules_added": add_training,
                "rules_before": before, "rules_after_training": after_training,
                "rules_after_averaging": count(),
            }
            family._kmc_rate_rule_preparation = copy.deepcopy(result)
            prepared[label] = result
    finally:
        owner.species_constraints = species_constraints
        owner.polymer_constraints = polymer_constraints
        rmg_input.rmg = original_input
    return {
        "procedure": "RMG add_rules_from_training then fill_rules_by_averaging_up",
        "families": prepared,
        "auto_generated_families_untouched": sorted(
            label for label, family in families.items()
            if getattr(family, "auto_generated", False)
        ),
    }


def _rate_rule_source(family, reaction):
    """Serialize RMG's own source extraction without replacing its comment."""
    training, source = family.extract_source_from_comments(reaction)

    def entry_data(entry):
        return {"index": entry.index, "label": entry.label, "rank": entry.rank,
                "short_desc": entry.short_desc}

    result = {"template": _template(reaction),
              "comment": reaction.kinetics.comment}
    if training:
        result["training"] = [{"entry": entry_data(source[1]),
                                "reverse": source[2], "weight": 1.0}]
        result["rules"] = []
    else:
        details = source[1]
        result["exact"] = details["exact"]
        result["rules"] = [
            {"entry": entry_data(entry), "weight": weight}
            for entry, weight in details["rules"]
        ]
        result["training"] = [
            {"rule": entry_data(rule), "entry": entry_data(entry), "weight": weight}
            for rule, entry, weight in details["training"]
        ]
    return result


def _persistent_labeled_u2(molecule, atom) -> tuple[bool, int]:
    """Return whether one labeled triplet u2 carbon survives every resonance form."""
    from rmgpy.species import Species

    if not (
        atom.element.number == 6
        and atom.charge == 0
        and atom.radical_electrons == 2
        and math.isclose(sum(bond.order for bond in atom.edges.values()), 2.0)
    ):
        return False, 0
    copied = molecule.copy(deep=True)
    marker = "*kmc_applicability_root"
    copied.atoms[molecule.atoms.index(atom)].label = marker
    species = Species(molecule=[copied])
    species.generate_resonance_structures()
    roots = []
    for form in species.molecule:
        matches = [candidate for candidate in form.atoms if candidate.label == marker]
        if len(matches) != 1:
            raise ValueError("resonance generation did not preserve mapped root identity")
        roots.append(matches[0])
    return (
        all(
            root.element.number == 6
            and root.charge == 0
            and root.radical_electrons == 2
            and math.isclose(sum(bond.order for bond in root.edges.values()), 2.0)
            for root in roots
        ),
        len(species.molecule),
    )


def _touched_atom_indices(operations) -> set[int]:
    """Return every stored atom index named by rewrite operations."""
    return {
        index
        for operation in operations
        for index in (
            operation.get("atoms", [])
            + ([operation["atom"]] if "atom" in operation else [])
        )
    }


def _select_family_root_candidate(
    family_candidates: set[int], preferred: set[int], label: str
) -> tuple[int, str | None]:
    """Select one jointly bound root, or retain a deterministic unresolved witness."""
    candidates = family_candidates & preferred if preferred else family_candidates
    if len(candidates) == 1:
        return next(iter(candidates)), None
    witness_pool = candidates or family_candidates or preferred
    witness = min(witness_pool) if witness_pool else -1
    if family_candidates and preferred and not candidates:
        reason = "stored rewrite does not touch the family-derived recipe root"
    else:
        reason = f"mapped recipe root {label} is ambiguous in the stored reaction"
    return witness, reason


def _mapped_reaction_u2_roots(
    family, reaction, record: EventRecord | dict[str, Any] | None = None
) -> list[dict[str, Any]]:
    """Map recipe-labelled u2 centres without replacing the fired resonance form."""
    mapping_reaction = _reaction_from_record(record) if record is not None else reaction
    if not any(
        atom.element.number == 6 and atom.radical_electrons == 2
        for participant in mapping_reaction.reactants + mapping_reaction.products
        for atom in _molecule(participant).atoms
    ):
        return []
    original_reactants = [
        (
            participant.molecule[0]
            if hasattr(participant, "molecule")
            else participant
        ).copy(deep=True)
        for participant in mapping_reaction.reactants
    ]
    labeled = copy.deepcopy(mapping_reaction)
    family.add_atom_labels_for_reaction(
        labeled, output_with_resonance=False, save_order=True
    )
    recipe_labels = {
        token
        for action in family.forward_recipe.actions
        for token in action
        if isinstance(token, str) and token.startswith("*")
    }
    side_molecules = {
        side: [
            participant.molecule[0]
            if hasattr(participant, "molecule")
            else participant
            for participant in getattr(labeled, f"{side}s")
        ]
        for side in ("reactant", "product")
    }
    reactant_atoms = [atom for molecule in original_reactants for atom in molecule.atoms]
    reactant_root_candidates: dict[str, set[int]] = {}
    for original, labeled_molecule in zip(
        original_reactants, side_molecules["reactant"]
    ):
        mappings = original.find_isomorphism(labeled_molecule, save_order=True)
        if not mappings:
            raise ValueError("family labeling changed a stored reactant graph")
        for mapping in mappings:
            for original_atom, labeled_atom in mapping.items():
                if labeled_atom.label in recipe_labels:
                    reactant_root_candidates.setdefault(
                        labeled_atom.label, set()
                    ).add(reactant_atoms.index(original_atom))
    preferred_by_side = {"reactant": set(), "product": set()}
    if record is not None:
        data = record.to_dict() if isinstance(record, EventRecord) else record
        touched_indices = _touched_atom_indices(data.get("bond_ops", []))
        preferred_by_side["product"] = {
            operation["atom"]
            for operation in data.get("bond_ops", [])
            if operation.get("action") == "set_radical"
            and operation.get("value") == 2
        }
        preferred_by_side["reactant"] = {
            index
            for index in touched_indices
            if reactant_atoms[index].element.number == 6
            and reactant_atoms[index].charge == 0
            and reactant_atoms[index].radical_electrons == 2
        }
    by_side = {}
    for side in ("reactant", "product"):
        by_side[side] = {
            atom.label: (molecule, atom)
            for molecule in side_molecules[side]
            for atom in molecule.atoms
            if atom.label in recipe_labels
        }
    roots = []
    for label in sorted(recipe_labels):
        reactant_item = by_side["reactant"].get(label)
        product_item = by_side["product"].get(label)
        if reactant_item is None or product_item is None:
            continue
        selected = []
        for side, item in (("reactant", reactant_item), ("product", product_item)):
            persistent, resonance_count = _persistent_labeled_u2(*item)
            if item[1].element.number == 6 and item[1].radical_electrons == 2:
                selected.append((side, persistent, resonance_count))
        for side, persistent, resonance_count in selected:
            family_candidates = reactant_root_candidates.get(label, set())
            preferred = preferred_by_side[side]
            root_index, mapping_error = _select_family_root_candidate(
                family_candidates, preferred, label
            )
            roots.append(
                {
                    "reactant_atom_index": root_index,
                    "recipe_label": label,
                    "record_role": side,
                    "family_forward_role": side,
                    "persistent_neutral_divalent_carbon": persistent,
                    "resonance_form_count": resonance_count,
                    **(
                        {
                            "mapping_verified": False,
                            "mapping_error": mapping_error,
                        }
                        if mapping_error is not None
                        else {}
                    ),
                }
            )
    return roots


def _training_source_domain(family, reaction, cache) -> list[dict[str, Any]] | None:
    """Serialize the actual calibrated contributors selected by an RMG estimate."""
    if getattr(family, "auto_generated", True):
        return None
    training, source = family.extract_source_from_comments(reaction)
    rule_entries = []
    if training:
        entries = [(source[1], 1.0, source[2], None)]
    else:
        details = source[1]
        entries = [
            (entry, float(weight), False, rule)
            for rule, entry, weight in details.get("training", [])
        ]
        rule_entries = details.get("rules", [])
    domain = []
    for entry, weight, reverse, rule in entries:
        key = (family.label, entry.index, bool(reverse))
        if key not in cache:
            source_reaction = copy.deepcopy(entry.item)
            if reverse:
                source_reaction.reactants, source_reaction.products = (
                    source_reaction.products,
                    source_reaction.reactants,
                )
            cache[key] = _mapped_reaction_u2_roots(family, source_reaction)
        mapped_roots = copy.deepcopy(cache[key])
        domain.append(
            {
                "source": {
                    "kind": "training",
                    "entry_index": entry.index,
                    "entry_label": entry.label,
                    "rule_index": getattr(rule, "index", None),
                    "weight": weight,
                    "stored_reverse": bool(reverse),
                },
                "mapped_roots": mapped_roots,
            }
        )
    domain.extend(
        {
            "source": {
                "kind": "rule",
                "entry_index": entry.index,
                "entry_label": entry.label,
                "weight": float(weight),
            },
            "mapped_roots": None,
        }
        for entry, weight in rule_entries
    )
    return domain or None


def _reaction_from_record(record: EventRecord | dict[str, Any]):
    """Reconstruct the exact stored direction for family root labeling."""
    from rmgpy.molecule.molecule import Molecule
    from rmgpy.reaction import Reaction
    from rmgpy.species import Species

    def make_side(graphs):
        return [
            Species(molecule=[Molecule().from_adjacency_list(graph)])
            for graph in graphs
        ]
    data = record.to_dict() if isinstance(record, EventRecord) else record
    return Reaction(
        reactants=make_side(data["reactant_graphs"]),
        products=make_side(data["product_graphs"]),
    )


def _stored_u2_root_fallback(record: EventRecord | dict[str, Any]) -> list[dict[str, Any]]:
    """Bind touched u2 atoms when family relabeling cannot survive normalization."""
    from rmgpy.molecule.molecule import Molecule

    data = record.to_dict() if isinstance(record, EventRecord) else record
    reactant_atoms = [
        atom
        for graph in data.get("reactant_graphs", [])
        for atom in Molecule().from_adjacency_list(graph).atoms
    ]
    touched = _touched_atom_indices(data.get("bond_ops", []))
    roots = []
    reactant_u2 = {
        index
        for index in touched
        if reactant_atoms[index].element.number == 6
        and reactant_atoms[index].charge == 0
        and reactant_atoms[index].radical_electrons == 2
    }
    product_u2 = {
        operation["atom"]
        for operation in data.get("bond_ops", [])
        if operation.get("action") == "set_radical"
        and operation.get("value") == 2
    }
    for side, indices in (("reactant", reactant_u2), ("product", product_u2)):
        for index in sorted(indices):
            roots.append(
                {
                    "reactant_atom_index": index,
                    "recipe_label": "__runtime_recipe_root_unresolved__",
                    "record_role": side,
                    "family_forward_role": side,
                    "mapping_caveat": (
                        "family relabeling failed; exact touched atom is bound to "
                        "the stored rewrite"
                    ),
                }
            )
    return roots


def _atom_total_order(atom) -> tuple[Any, ...]:
    """Explicit order for RMG atoms; never inherit container iteration order."""
    identifier = getattr(atom, "id", None)
    if not isinstance(identifier, int):
        identifier = -1
    neighbors = tuple(
        sorted(
            (
                getattr(neighbor.element, "number", 0),
                getattr(neighbor, "id", -1),
                str(bond.order),
            )
            for neighbor, bond in getattr(atom, "edges", {}).items()
        )
    )
    return (
        getattr(atom.element, "number", 0),
        getattr(atom.element, "isotope", -1),
        getattr(atom, "radical_electrons", 0),
        getattr(atom, "lone_pairs", 0),
        getattr(atom, "charge", 0),
        getattr(atom, "label", ""),
        neighbors,
        identifier,
    )


def _ordered_atoms(molecule) -> list[Any]:
    if molecule.__class__.__module__.startswith("rmgpy."):
        from rdkit import Chem
        from rmgpy.molecule.converter import to_rdkit_mol

        rdkit_molecule, mapping = to_rdkit_mol(
            molecule,
            remove_h=False,
            return_mapping=True,
            save_order=True,
        )
        ranks = list(
            Chem.CanonicalRankAtoms(
                rdkit_molecule,
                breakTies=True,
                includeChirality=True,
                includeIsotopes=True,
            )
        )
        if len(ranks) != len(set(ranks)):
            raise ValueError("RDKit did not produce a canonical total atom order")
        return sorted(molecule.atoms, key=lambda atom: ranks[mapping[atom]])
    return sorted(molecule.atoms, key=_atom_total_order)


def _canonical_adjacency(molecule) -> str:
    clone = molecule.copy(deep=True)
    clone.atoms = _ordered_atoms(clone)
    return clone.to_adjacency_list(remove_h=False)


def _molecule_total_order(molecule) -> tuple[str, str]:
    clone = molecule.copy(deep=True)
    smiles = clone.to_smiles() if hasattr(clone, "to_smiles") else ""
    return smiles, _canonical_adjacency(molecule)


def _canonical_resonance_representative(molecule):
    """Collapse arbitrary Kekule localization to RMG's aromatic graph form."""
    if not molecule.__class__.__module__.startswith("rmgpy."):
        return molecule
    from rmgpy.molecule.resonance import generate_aromatic_resonance_structure

    aromatic = generate_aromatic_resonance_structure(
        molecule, copy=True, save_order=True
    )
    return aromatic[0] if aromatic else molecule


def _molecule(species_or_molecule):
    cache = _PAIR_MOLECULE_CACHE.get()
    identity = id(species_or_molecule)
    if cache is not None and identity in cache:
        return cache[identity][1]
    if not hasattr(species_or_molecule, "molecule"):
        selected = _canonical_resonance_representative(species_or_molecule)
    else:
        molecules = species_or_molecule.molecule
        if not molecules:
            raise ValueError("compiled species has no molecular graph")
        candidates = [_canonical_resonance_representative(item) for item in molecules]
        selected = min(candidates, key=_molecule_total_order)
    if cache is not None:
        cache[identity] = species_or_molecule, selected
    return selected


def _formula(molecules: Iterable) -> dict[str, int]:
    result: dict[str, int] = {}
    for item in molecules:
        for element, count in _molecule(item).get_element_count().items():
            result[element] = result.get(element, 0) + int(count)
    return dict(sorted(result.items()))


def _radicals(molecules: Iterable) -> int:
    return sum(
        atom.radical_electrons for item in molecules for atom in _molecule(item).atoms
    )


def _implicit_h(molecules: Iterable) -> int:
    return sum(
        getattr(atom, "implicit_hydrogens", 0)
        for item in molecules
        for atom in _molecule(item).atoms
    )


def _bonds(molecules: Iterable) -> dict[tuple[int, int], str]:
    result = {}
    for item in molecules:
        for bond in _molecule(item).get_all_edges():
            a, b = sorted((bond.vertex1.id, bond.vertex2.id))
            result[(a, b)] = str(bond.order)
    return result


def _canonical_atom_map(reaction) -> dict[int, int]:
    reactant_atoms = [
        atom
        for item in reaction.reactants
        for atom in _ordered_atoms(_molecule(item))
        if atom.element.number != 1
    ]
    product_atoms = [
        atom
        for item in reaction.products
        for atom in _ordered_atoms(_molecule(item))
        if atom.element.number != 1
    ]
    reactant_index = {atom.id: index for index, atom in enumerate(reactant_atoms)}
    product_index = {atom.id: index for index, atom in enumerate(product_atoms)}
    if set(reactant_index) != set(product_index):
        raise ValueError("heavy-atom ids are not conserved by the generated reaction")
    return {
        reactant_index[identifier]: product_index[identifier]
        for identifier in sorted(reactant_index)
    }


def _bond_ops(reaction) -> list[dict[str, Any]]:
    reactant_atoms = [
        atom for item in reaction.reactants for atom in _ordered_atoms(_molecule(item))
    ]
    reactant_index = {atom.id: index for index, atom in enumerate(reactant_atoms)}
    product_atoms = [
        atom for item in reaction.products for atom in _ordered_atoms(_molecule(item))
    ]
    product_by_id = {atom.id: atom for atom in product_atoms}
    before, after = _bonds(reaction.reactants), _bonds(reaction.products)
    ops = []
    for pair in sorted(set(before) | set(after)):
        old, new = before.get(pair), after.get(pair)
        if old == new:
            continue
        if not all(identifier in reactant_index for identifier in pair):
            raise ValueError(
                "bond rewrite references an atom absent from the reactants"
            )
        atoms = [reactant_index[identifier] for identifier in pair]
        if old is not None:
            ops.append({"action": "break", "atoms": atoms, "order": old})
        if new is not None:
            ops.append({"action": "form", "atoms": atoms, "order": new})
    attributes = (
        ("radical_electrons", "set_radical"),
        ("charge", "set_charge"),
        ("lone_pairs", "set_lone_pairs"),
        ("implicit_hydrogens", "set_implicit_hydrogens"),
    )
    for identifier, index in sorted(reactant_index.items(), key=lambda item: item[1]):
        reactant_atom = reactant_atoms[index]
        product_atom = product_by_id[identifier]
        for attribute, action in attributes:
            old = getattr(reactant_atom, attribute, 0)
            new = getattr(product_atom, attribute, 0)
            if old != new:
                ops.append({"action": action, "atom": index, "value": int(new)})
    return ops


def _template(reaction) -> str:
    template = getattr(reaction, "template", None) or []
    return ";".join(getattr(item, "label", str(item)) for item in template)


def _proxy_fingerprint(proxy) -> dict[str, Any]:
    """Structural, rather than object-identity, provenance for a proxy input."""
    return {
        "site_type": proxy.site_type,
        "frontier": proxy.frontier,
        "metadata": proxy.metadata,
        "participant_site_types": list(proxy.participant_site_types),
        "reactants": [
            _canonical_adjacency(_molecule(item)) for item in proxy.reactants
        ],
    }


def _graph_adjacencies(participants: Iterable[Any]) -> list[str]:
    """Serialize participant graphs as stable, explicit-H adjacency lists."""
    return [_canonical_adjacency(_molecule(item)) for item in participants]


def _copy_participant(item):
    """Give the public generator a pristine participant on every invocation."""
    return item.copy(deep=True) if hasattr(item, "copy") else item


def _is_resonance_derived(reaction) -> bool:
    """RMG exposes resonance products as multi-structure Species objects."""
    return any(
        len(getattr(product, "molecule", ())) > 1 for product in reaction.products
    )


def _cyclic_carbon_ids(participants: Iterable[Any]) -> set[int]:
    result = set()
    for participant in participants:
        molecule = _molecule(participant)
        if not hasattr(molecule, "is_atom_in_cycle"):
            continue
        result.update(
            atom.id
            for atom in molecule.atoms
            if atom.element.number == 6 and molecule.is_atom_in_cycle(atom)
        )
    return result


def _inventory_class(reaction) -> str | None:
    family = getattr(reaction, "family", "")
    before = _bonds(reaction.reactants)
    aromatic_atoms = _cyclic_carbon_ids(reaction.reactants)
    changed_aromatic = False
    formed_at_aromatic = False
    after = _bonds(reaction.products)
    for pair in set(before) | set(after):
        if before.get(pair) != after.get(pair) and any(
            atom in aromatic_atoms for atom in pair
        ):
            changed_aromatic = True
        if (
            pair not in before
            and pair in after
            and any(atom in aromatic_atoms for atom in pair)
        ):
            formed_at_aromatic = True
    if family == "R_Recombination" and formed_at_aromatic:
        return "R1:J_ring"
    if (
        family == "Disproportionation"
        and changed_aromatic
        and _is_resonance_derived(reaction)
    ):
        return "R1:quinoid_disproportionation"
    return None


def _junction_rewrite(
    action: str, attacked: str, attacker: str
) -> list[dict[str, Any]]:
    if action == "create":
        return [
            {"action": "form", "atoms": [attacked, attacker], "order": "1.0"},
            {"action": "set_radical", "atom": attacked, "value": 0},
            {"action": "set_radical", "atom": attacker, "value": 0},
        ]
    if action == "dissociate":
        return [
            {"action": "break", "atoms": [attacked, attacker], "order": "1.0"},
            {"action": "set_radical", "atom": attacked, "value": 1},
            {"action": "set_radical", "atom": attacker, "value": 1},
        ]
    raise ValueError(f"unknown junction action: {action}")


def _junction_ops(
    inventory_class: str | None,
    metadata: dict[str, Any],
    *,
    action: str = "create",
    reverse_event_handle: str = LINK_PLACEHOLDER,
) -> list[dict[str, Any]] | None:
    if inventory_class != "R1:J_ring":
        return None
    junction_kind = metadata.get("junction_kind", "J_para")
    try:
        localisation = R1_JUNCTION_LABELS[junction_kind]
    except KeyError as error:
        raise ValueError(f"unknown R1 junction kind: {junction_kind}") from error
    attacked = localisation["attacked_atom_label"]
    attacker = metadata.get("attacker_site_label", "P9")
    inverse_action = "dissociate" if action == "create" else "create"
    return [
        {
            "action": action,
            "junction_kind": junction_kind,
            "attacker_site_label": attacker,
            "attacked_atom_label": attacked,
            "formed_bond": {"label_pair": [attacked, attacker], "order": "1.0"},
            "chosen_resonance_localisation": localisation[
                "chosen_resonance_localisation"
            ],
            "reflection": localisation["reflection"],
            "reverse_event_handle": reverse_event_handle,
            "exact_inverse_rewrite": _junction_rewrite(
                inverse_action, attacked, attacker
            ),
            "aromaticity": "dearomatised" if action == "create" else "restore",
            "eligibility_refresh": True,
            "subsequent_chemistry_eligible": True,
        }
    ]


def _partition_degeneracy(record: EventRecord, divisor: int) -> EventRecord:
    """Partition one degeneracy-inclusive channel among equivalent records."""
    if divisor < 1:
        raise ValueError("degeneracy divisor must be positive")
    k_table = copy.deepcopy(record.k_table)
    if k_table is not None:
        k_table["k"] = [value / divisor for value in k_table["k"]]
    return replace(
        record,
        raw_path_degeneracy=record.raw_path_degeneracy / divisor,
        degeneracy=record.degeneracy / divisor,
        k_table=k_table,
        event_id="",
    )


def _reverse_view(reaction):
    source = getattr(reaction, "source_reaction", reaction)
    return SimpleNamespace(
        reactants=reaction.products,
        products=reaction.reactants,
        family=getattr(reaction, "family", ""),
        template=getattr(reaction, "template", []),
        degeneracy=getattr(reaction, "degeneracy", 1.0),
        reversible=True,
        kinetics=None,
        source_reaction=source,
    )


def _pair_convention(reactants: Sequence[Any]) -> str:
    if len(reactants) != 2:
        return "ordered-unimolecular"
    first, second = reactants
    same = first is second or (
        hasattr(first, "is_isomorphic") and first.is_isomorphic(second)
    )
    return "unordered-identical-pair N(N-1)/2" if same else "distinct-pair N_A*N_B"


def _participant_inventory(proxy, reaction) -> tuple[list[str], list[int]]:
    declared = list(
        proxy.participant_site_types or [proxy.site_type] * len(reaction.reactants)
    )
    if len(declared) != len(reaction.reactants):
        raise ValueError("proxy participant site types do not match reaction arity")
    site_types = []
    multiplicities = []
    for site_type in declared:
        if site_type in site_types:
            multiplicities[site_types.index(site_type)] += 1
        else:
            site_types.append(site_type)
            multiplicities.append(1)
    return site_types, multiplicities


def _same_participants(first: Sequence[Any], second: Sequence[Any]) -> bool:
    if len(first) != len(second):
        return False
    remaining = list(second)
    for item in first:
        for index, candidate in enumerate(remaining):
            left, right = _molecule(item), _molecule(candidate)
            same = item is candidate or left is right
            if not same and hasattr(item, "is_isomorphic"):
                same = (
                    item.is_isomorphic(candidate, save_order=True)
                    if item.__class__.__module__.startswith("rmgpy.")
                    else item.is_isomorphic(candidate)
                )
            if not same and hasattr(left, "is_isomorphic"):
                same = (
                    left.is_isomorphic(right, save_order=True)
                    if left.__class__.__module__.startswith("rmgpy.")
                    else left.is_isomorphic(right)
                )
            if (
                not same
                and hasattr(left, "to_adjacency_list")
                and hasattr(right, "to_adjacency_list")
            ):
                same = left.to_adjacency_list(
                    remove_h=False
                ) == right.to_adjacency_list(remove_h=False)
            if same:
                del remaining[index]
                break
        else:
            return False
    return True


def _graph_lists_isomorphic(first: Sequence[str], second: Sequence[str]) -> bool:
    from rmgpy.molecule.molecule import Molecule

    if len(first) != len(second):
        return False
    remaining = [Molecule().from_adjacency_list(adjacency) for adjacency in second]
    for adjacency in first:
        molecule = Molecule().from_adjacency_list(adjacency)
        for index, candidate in enumerate(remaining):
            if molecule.is_isomorphic(candidate):
                del remaining[index]
                break
        else:
            return False
    return True


def _graph_list_key(graphs: Sequence[str]) -> tuple[str, ...]:
    """Canonical lookup key; candidate matches are still checked isomorphically."""
    from rmgpy.molecule.molecule import Molecule

    return tuple(
        sorted(Molecule().from_adjacency_list(graph).to_smiles() for graph in graphs)
    )


def _orient_to_proxy(proxy, reaction):
    """Return the generated reaction in the direction initiated by the proxy."""
    if _same_participants(proxy.reactants, reaction.reactants):
        return reaction, "as_generated"
    reverse = _same_participants(proxy.reactants, reaction.products)
    if not reverse and getattr(reaction, "is_forward", None) is False:
        # RMG canonicalizes own_reverse=False families.  A resonance structure
        # of the initiating proxy may then be non-isomorphic to either stored
        # side, but is_forward still records that the input fired in reverse.
        reverse = True
    if reverse:
        attributes = {
            name: getattr(reaction, name)
            for name in ("family", "template", "degeneracy", "reversible", "is_forward")
            if hasattr(reaction, name)
        }
        return (
            SimpleNamespace(
                reactants=reaction.products,
                products=reaction.reactants,
                kinetics=None,
                source_reaction=reaction,
                **attributes,
            ),
            "reversed",
        )
    if getattr(reaction, "is_forward", None) is True:
        return reaction, "as_generated"
    raise ValueError(
        f"generated {reaction.family} reaction does not contain its centred proxy"
    )


def ps_proxy_set(units: int = PS_PROXY_UNITS) -> tuple[SiteProxy, ...]:
    """Centred PS structures; discover firing families through RMG generation.

    Retain the original three-unit structures and add full five-unit contexts
    for the radius-one completeness bound.  All five candidate families are
    tried on every structure, including radical/closed-chain donor pairs.
    """
    from rmgpy.molecule.molecule import Molecule
    from rmgpy.species import Species

    if units < 3:
        raise ValueError("PS proxy contexts require at least three repeat units")
    structures = {
        "pristine": linear_ps_smiles(units),
        "interior_radical": linear_ps_smiles(units, radical_unit=units // 2),
        "doubly_featured": linear_ps_smiles(units, radical_units=(0, 1)),
        "end_radical": linear_ps_smiles(units, radical_end=True),
        "end_radical_short": linear_ps_smiles(units - 1, radical_end=True),
        "benzylic_end_radical": benzylic_ps_end_smiles(units),
        "benzylic_end_radical_short": benzylic_ps_end_smiles(units - 1),
        "junction_radical": linear_ps_smiles(units).replace(
            "c1ccccc1", "c1cc[c]cc1", 1
        ),
        "styrene": "C=Cc1ccccc1",
    }
    species = {}
    for name, smiles in structures.items():
        molecule = Molecule(smiles=smiles)
        molecule.update()
        species[name] = Species(molecule=[molecule])
    unimolecular_candidates = PS_FAMILY_CANDIDATES
    radical_pair_candidates = PS_FAMILY_CANDIDATES
    declarations = (
        (
            "pristine",
            ("pristine",),
            ("pristine",),
            "pristine",
            PS_FAMILY_CANDIDATES,
        ),
        (
            "interior_radical",
            ("interior_radical",),
            ("interior_radical",),
            "featured",
            unimolecular_candidates,
        ),
        (
            "doubly_featured",
            ("doubly_featured",),
            ("doubly_featured",),
            "featured",
            unimolecular_candidates,
        ),
        (
            "end_radical",
            ("end_radical",),
            ("end_radical",),
            "end-proximal",
            unimolecular_candidates,
        ),
        (
            "junction_radical",
            ("junction_radical",),
            ("junction_radical",),
            "junction",
            PS_FAMILY_CANDIDATES,
        ),
        (
            "end_radical+end_radical",
            ("end_radical", "end_radical"),
            ("end_radical", "end_radical"),
            "end-proximal",
            radical_pair_candidates,
        ),
        (
            "junction_radical+end_radical",
            ("junction_radical", "end_radical"),
            ("junction_radical", "end_radical"),
            "junction",
            PS_FAMILY_CANDIDATES,
        ),
        (
            "end_radical+styrene",
            ("end_radical_short" if units == 3 else "end_radical", "styrene"),
            ("end_radical", "styrene"),
            "end-proximal",
            PS_FAMILY_CANDIDATES,
        ),
    ) + (
        (
            "benzylic_end_radical",
            ("benzylic_end_radical",),
            ("benzylic_end_radical",),
            "end-proximal",
            PS_FAMILY_CANDIDATES,
        ),
        (
            "benzylic_end_radical+styrene",
            ("benzylic_end_radical_short", "styrene"),
            ("benzylic_end_radical", "styrene"),
            "end-proximal",
            PS_FAMILY_CANDIDATES,
        ),
        (
            "benzylic_end_radical+pristine",
            ("benzylic_end_radical", "pristine"),
            ("benzylic_end_radical", "pristine"),
            "end-proximal",
            PS_FAMILY_CANDIDATES,
        ),
        (
            "benzylic_end_radical+end_radical",
            ("benzylic_end_radical", "end_radical"),
            ("benzylic_end_radical", "end_radical"),
            "end-proximal",
            PS_FAMILY_CANDIDATES,
        ),
        (
            "benzylic_end_radical+benzylic_end_radical",
            ("benzylic_end_radical", "benzylic_end_radical"),
            ("benzylic_end_radical", "benzylic_end_radical"),
            "end-proximal",
            PS_FAMILY_CANDIDATES,
        ),
        (
            "junction_radical+benzylic_end_radical",
            ("junction_radical", "benzylic_end_radical"),
            ("junction_radical", "benzylic_end_radical"),
            "junction",
            PS_FAMILY_CANDIDATES,
        ),
    ) + tuple(
        (
            f"{radical}+pristine",
            (radical, "pristine"),
            (radical, "pristine"),
            "junction" if radical == "junction_radical" else "featured",
            PS_FAMILY_CANDIDATES,
        )
        for radical in (
            "interior_radical",
            "end_radical",
            "junction_radical",
        )
    )
    proxies = tuple(
        SiteProxy(
            site_type,
            tuple(species[name] for name in participants),
            participant_site_types=participant_site_types,
            metadata={
                "centred": True,
                "context_class": context_class,
                "proxy_units": units,
                "family_candidates": list(proxy_candidates),
                **(
                    {
                        "junction_kind": "J_para",
                        "attacker_site_label": "P9",
                    }
                    if site_type == "junction_radical+end_radical"
                    else {}
                ),
            },
        )
        for (
            site_type,
            participants,
            participant_site_types,
            context_class,
            proxy_candidates,
        ) in declarations
    )
    if units == 3:
        extended = tuple(
            replace(
                proxy,
                site_type=f"{proxy.site_type}@5",
                metadata={**proxy.metadata, "coverage_site_type": proxy.site_type},
            )
            for proxy in ps_proxy_set(5)
        )
        return proxies + extended
    return proxies


def _archived_junction_proxy(junction_kind: str) -> SiteProxy:
    """Return the mapped cage pair used by an archived R-009 junction oracle."""
    from rmgpy.molecule.molecule import Molecule
    from rmgpy.species import Species

    secondary_smiles = (
        "CCCC=C1C=C[CH]C=C1" if junction_kind == "J_para" else "CCCC=C1[CH]C=CC=C1"
    )
    reactants = []
    for smiles in (secondary_smiles, "[CH2]C(C)c1ccccc1"):
        molecule = Molecule(smiles=smiles)
        molecule.update()
        molecule.assign_atom_ids()
        reactants.append(Species(molecule=[molecule]))
    return SiteProxy(
        "junction_radical+end_radical",
        tuple(reactants),
        participant_site_types=("junction_radical", "end_radical"),
        metadata={
            "centred": True,
            "context_class": "junction",
            "proxy_units": PS_PROXY_UNITS,
            "family_candidates": ["R_Recombination"],
            "junction_kind": junction_kind,
            "attacker_site_label": "P9",
        },
    )


def _ortho_junction_proxy(junction_kind: str = "J_ortho_S7") -> SiteProxy:
    """Return the mapped cage pair used by the archived J_ortho oracle."""
    return _archived_junction_proxy(junction_kind)


def _para_junction_proxy() -> SiteProxy:
    """Return the mapped cage pair used by the archived J_para oracle."""
    return _archived_junction_proxy("J_para")


def _archived_junction_reaction(
    proxy: SiteProxy,
    *,
    preexponential: float,
    degeneracy: float,
    template: Sequence[str],
    product_join: float,
    product_low: Sequence[float],
    product_high: Sequence[float],
):
    """Reconstruct one frozen RMG-mapped R-009 junction reaction."""
    from rmgpy.data.kinetics.family import TemplateReaction
    from rmgpy.kinetics.arrhenius import Arrhenius
    from rmgpy.molecule.molecule import Bond
    from rmgpy.species import Species
    from rmgpy.thermo.nasa import NASA, NASAPolynomial

    def nasa(join, low, high):
        return NASA(
            polynomials=[
                NASAPolynomial(coeffs=low, Tmin=(300.0, "K"), Tmax=(join, "K")),
                NASAPolynomial(coeffs=high, Tmin=(join, "K"), Tmax=(3000.0, "K")),
            ],
            Tmin=(300.0, "K"),
            Tmax=(3000.0, "K"),
            comment="R-009 archived RMG group-additivity thermo",
        )

    secondary = _molecule(proxy.reactants[0]).copy(deep=True)
    primary = _molecule(proxy.reactants[1]).copy(deep=True)
    attacked = next(
        atom
        for atom in secondary.atoms
        if atom.radical_electrons and secondary.is_atom_in_cycle(atom)
    )
    attacker = next(atom for atom in primary.atoms if atom.radical_electrons)
    product = secondary.merge(primary)
    product.add_bond(Bond(attacked, attacker, order=1.0))
    attacked.radical_electrons = 0
    attacker.radical_electrons = 0
    product.update(sort_atoms=False)
    proxy.reactants[0].thermo = nasa(
        1254.2536410728535,
        (
            -0.9486839628215864,
            0.08486284068971323,
            -3.938637982344858e-05,
            2.2192613672944045e-09,
            1.9058827184182995e-12,
            13488.681882038003,
            32.230444564008096,
        ),
        (
            2.8546259360984383,
            0.08172943981695815,
            -4.6397539366542894e-05,
            1.1664266596876454e-08,
            -1.1165136655519307e-12,
            11827.021589717222,
            10.199278347265654,
        ),
    )
    proxy.reactants[1].thermo = nasa(
        1177.6065961907934,
        (
            -5.607484013799355,
            0.09820532625745784,
            -7.478953987950637e-05,
            2.784763043057718e-08,
            -3.6642297619144356e-12,
            22890.317829865562,
            55.06864760491035,
        ),
        (
            12.837928688294944,
            0.04618929187592041,
            -2.208327352961384e-05,
            5.680597353073413e-09,
            -5.868026777213613e-13,
            17808.42272068456,
            -40.067143782909014,
        ),
    )
    product_species = Species(molecule=[product])
    product_species.thermo = nasa(product_join, product_low, product_high)
    for species in [*proxy.reactants, product_species]:
        species.props[FROZEN_THERMO_PROPERTY] = True
    return TemplateReaction(
        reactants=list(proxy.reactants),
        products=[product_species],
        kinetics=Arrhenius(
            A=(preexponential, "m^3/(mol*s)"),
            n=-1.00291,
            Ea=(0.0, "J/mol"),
            T0=(1.0, "K"),
        ),
        reversible=True,
        degeneracy=degeneracy,
        family="R_Recombination",
        template=list(template),
        is_forward=True,
    )


def _archived_para_reaction(proxy: SiteProxy):
    """Reconstruct the frozen RMG-mapped g=1 para capture reaction."""
    return _archived_junction_reaction(
        proxy,
        preexponential=1.76793e10,
        degeneracy=1.0,
        template=(ARCHIVED_J_PARA_RULE,),
        product_join=1403.6826169506783,
        product_low=(
            -10.949627764777754,
            0.20969652550530965,
            -0.00015924632577030666,
            6.232354330618766e-08,
            -9.94128604215258e-12,
            11636.982987189213,
            85.42480168403608,
        ),
        product_high=(
            22.36559012247251,
            0.11475767839669494,
            -5.779061059436066e-05,
            1.4136912533315018e-08,
            -1.3589033989204612e-12,
            2284.404766054439,
            -86.59824668710979,
        ),
    )


def _archived_ortho_reaction(proxy: SiteProxy):
    """Reconstruct the frozen RMG-mapped g=2 ortho capture reaction."""
    return _archived_junction_reaction(
        proxy,
        preexponential=3.53586e10,
        degeneracy=2.0,
        template=(),
        product_join=1171.6898321871104,
        product_low=(
            -8.455190642047764,
            0.19501882679191176,
            -0.00013017705248095492,
            3.885484757972701e-08,
            -3.261608985996416e-12,
            13316.25635917201,
            75.03280013132829,
        ),
        product_high=(
            18.46952197543009,
            0.12343415887301927,
            -6.45643901319227e-05,
            1.633326197511283e-08,
            -1.6163446936598018e-12,
            5611.067403002024,
            -65.08547018371374,
        ),
    )


def _reflect_ortho_reaction(reaction):
    """Apply the archived S6/S7 and S8/S10 graph automorphism."""
    from rmgpy.molecule.molecule import Bond

    secondary = _molecule(reaction.reactants[0])
    attacked = next(
        atom
        for atom in secondary.atoms
        if atom.radical_electrons and secondary.is_atom_in_cycle(atom)
    )
    ring = set(
        next(
            cycle
            for cycle in secondary.get_smallest_set_of_smallest_rings()
            if attacked in cycle
        )
    )
    attacked_ring_neighbors = [atom for atom in attacked.edges if atom in ring]
    ipso = next(
        atom
        for atom in attacked_ring_neighbors
        if any(
            neighbor not in ring and neighbor.element.number != 1
            for neighbor in atom.edges
        )
    )
    attacked_meta = next(atom for atom in attacked_ring_neighbors if atom is not ipso)
    reflected_ortho = next(
        atom for atom in ipso.edges if atom in ring and atom is not attacked
    )
    reflected_meta = next(
        atom for atom in reflected_ortho.edges if atom in ring and atom is not ipso
    )
    reflection = {
        attacked.id: reflected_ortho.id,
        reflected_ortho.id: attacked.id,
        attacked_meta.id: reflected_meta.id,
        reflected_meta.id: attacked_meta.id,
    }
    reflected = reaction.copy()
    for participant in [*reflected.reactants, *reflected.products]:
        reflected_molecules = []
        for molecule in participant.molecule:
            source_by_id = {atom.id: atom for atom in molecule.atoms}
            edge_specs = []
            for bond in molecule.get_all_edges():
                edge_specs.append(
                    (
                        reflection.get(bond.vertex1.id, bond.vertex1.id),
                        reflection.get(bond.vertex2.id, bond.vertex2.id),
                        float(bond.order),
                    )
                )
            clone = molecule.copy(deep=True)
            clone_by_id = {atom.id: atom for atom in clone.atoms}
            for bond in list(clone.get_all_edges()):
                clone.remove_bond(bond)
            for first, second, order in edge_specs:
                clone.add_bond(Bond(clone_by_id[first], clone_by_id[second], order))
            for source_id, source_atom in source_by_id.items():
                target = clone_by_id[reflection.get(source_id, source_id)]
                for attribute in (
                    "radical_electrons",
                    "charge",
                    "lone_pairs",
                ):
                    setattr(target, attribute, getattr(source_atom, attribute, 0))
            clone.update(sort_atoms=False)
            reflected_molecules.append(clone)
        participant.molecule = reflected_molecules
    return reflected


def linear_ps_smiles(
    units: int,
    *,
    radical_unit: int | None = None,
    radical_units: Sequence[int] = (),
    radical_end: bool = False,
) -> str:
    """Build an H-capped constitution-only linear PS oligomer."""
    if units < 1:
        raise ValueError("a PS oligomer must contain at least one repeat unit")
    if radical_unit is not None and radical_units:
        raise ValueError("use one multiple-feature specification at a time")
    if radical_end and (radical_unit is not None or radical_units):
        raise ValueError("primary end and unit radical choices are incompatible")
    radical_sites = set(radical_units)
    if radical_unit is not None:
        radical_sites.add(radical_unit)
    if any(index < 0 or index >= units for index in radical_sites):
        raise ValueError("radical unit lies outside the oligomer")
    pieces = []
    for index in range(units):
        pieces.append("[CH2]" if radical_end and index == 0 else "C")
        pieces.append("[C](c1ccccc1)" if index in radical_sites else "C(c1ccccc1)")
    return "".join(pieces)


def benzylic_ps_end_smiles(units: int) -> str:
    """Build an H-capped head-to-tail PS chain with terminal CH•(Ph)."""
    if units < 1:
        raise ValueError("a PS oligomer must contain at least one repeat unit")
    return "CC(c1ccccc1)" * (units - 1) + "C[CH](c1ccccc1)"


def short_ps_molecule_catalogue(radius: int, feature_cap: int = 2) -> dict[str, Any]:
    """Finite whole-molecule catalogue for the PS homopolymer under span ``r``."""
    from rmgpy.molecule.molecule import Molecule

    if radius < 0:
        raise ValueError("radius must be non-negative")
    if feature_cap < 0:
        raise ValueError("feature cap must be non-negative")
    maximum = 2 * radius + 1
    molecules = []
    termination_bound = 0
    for units in range(1, maximum + 1):
        base = Molecule(smiles=linear_ps_smiles(units))
        backbone = [
            atom
            for atom in base.atoms
            if atom.element.symbol == "C" and not base.is_atom_in_cycle(atom)
        ]
        backbone_sites = len(backbone)
        if backbone_sites != 2 * units:
            raise ValueError("linear PS builder produced the wrong backbone size")
        for feature_count in range(min(feature_cap, backbone_sites) + 1):
            candidates = list(combinations(range(backbone_sites), feature_count))
            termination_bound += len(candidates)
            seen: dict[str, list[Any]] = {}
            for sites in candidates:
                molecule = base.copy(deep=True)
                molecule_backbone = [
                    atom
                    for atom in molecule.atoms
                    if atom.element.symbol == "C"
                    and not molecule.is_atom_in_cycle(atom)
                ]
                for site in sites:
                    atom = molecule_backbone[site]
                    hydrogens = [
                        neighbor
                        for neighbor in atom.edges
                        if neighbor.element.number == 1
                    ]
                    if not hydrogens:
                        raise ValueError(
                            "radical feature exceeds the site's hydrogen inventory"
                        )
                    molecule.remove_atom(hydrogens[0])
                    atom.radical_electrons += 1
                molecule.update(sort_atoms=False)
                # Index reflection swaps phenyl-bearing and plain backbone C.
                # Quotient only actual featured graphs, never index tuples.
                key = molecule.to_smiles()
                equivalent = seen.setdefault(key, [])
                if any(molecule.is_isomorphic(other) for other in equivalent):
                    continue
                equivalent.append(molecule)
                molecules.append(
                    {
                        "repeat_units": units,
                        "formula": dict(sorted(molecule.get_element_count().items())),
                        "radical_sites": list(sites),
                        "state": (
                            "closed_shell_linear"
                            if feature_count == 0
                            else "radical_linear"
                        ),
                        "adjacency_list": molecule.to_adjacency_list(remove_h=False),
                    }
                )
    return {
        "radius": radius,
        "feature_cap": feature_cap,
        "max_units": maximum,
        "termination_bound": termination_bound,
        "terminates": len(molecules) <= termination_bound,
        "molecules": molecules,
        "size": len(molecules),
    }


def ceiling_temperature(
    propagation: dict[str, Any],
    depropagation: dict[str, Any],
    monomer_concentration: float,
) -> float | None:
    """Return the tabulated crossing of k_prop[M] and k_deprop, if bracketed."""
    prop, dep = propagation.get("k_table"), depropagation.get("k_table")
    if not prop or not dep or prop["T"] != dep["T"]:
        return None
    if monomer_concentration <= 0:
        raise ValueError("monomer concentration must be positive")
    residual = [
        math.log(kp * monomer_concentration / kd) for kp, kd in zip(prop["k"], dep["k"])
    ]
    for index in range(len(residual) - 1):
        if residual[index] == 0:
            return prop["T"][index]
        if residual[index] * residual[index + 1] < 0:
            fraction = -residual[index] / (residual[index + 1] - residual[index])
            return prop["T"][index] + fraction * (
                prop["T"][index + 1] - prop["T"][index]
            )
    return None


def validate_artifact(artifact: dict[str, Any]) -> None:
    """Reject structural corruption in a compiled event-set artifact."""
    if artifact.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("unknown event-set schema")
    records = artifact.get("records", [])
    by_id = {record.get("event_id"): record for record in records}
    ortho_records = []
    for record in records:
        operations = record.get("junction_ops")
        if not isinstance(operations, list) or not operations:
            continue
        operation = operations[0]
        if not isinstance(operation, dict):
            continue
        if operation.get("junction_kind") in {"J_ortho_S6", "J_ortho_S7"}:
            ortho_records.append((record, operation))
    if "R_Recombination" in artifact.get("families", []) or ortho_records:
        expected_labels = {"S6", "S7"}
        forwards_by_label = {
            label: [
                (record, operation)
                for record, operation in ortho_records
                if operation.get("action") == "create"
                and operation.get("attacked_atom_label") == label
            ]
            for label in expected_labels
        }
        incomplete = [
            label
            for label, forwards in sorted(forwards_by_label.items())
            if len(forwards) != 1
        ]
        if incomplete:
            raise ValueError(
                "J_ortho completeness requires exactly one forward record for "
                f"S6 and S7; invalid or missing {', '.join(incomplete)}"
            )
        for label, forwards in sorted(forwards_by_label.items()):
            forward, _ = forwards[0]
            partner = by_id.get(forward.get("reverse_of"))
            partner_operations = partner.get("junction_ops") if partner else None
            partner_operation = (
                partner_operations[0]
                if isinstance(partner_operations, list) and partner_operations
                else {}
            )
            if (
                partner is None
                or partner.get("reverse_of") != forward.get("event_id")
                or partner_operation.get("action") != "dissociate"
                or partner_operation.get("attacked_atom_label") != label
            ):
                raise ValueError(
                    f"J_ortho completeness requires S6 and S7 reverse_of partners; "
                    f"invalid or missing {label} reverse"
                )
    event_ids = [record.get("event_id") for record in records]
    legacy_ids = [
        record.get("event_id")
        for record in records
        if not (
            record.get("junction_ops")
            and record["junction_ops"][0].get("junction_kind", "").startswith("J_ortho")
        )
    ]
    ortho_ids = [event_id for event_id in event_ids if event_id not in legacy_ids]
    expected_order = sorted(legacy_ids) + sorted(ortho_ids)
    if event_ids != expected_order or len(event_ids) != len(set(event_ids)):
        raise ValueError("records must have unique canonical event-id ordering")
    provenance = artifact.get("provenance", {})
    for index, data in enumerate(records):
        record = EventRecord.from_dict(data)
        record.validate()
        if record.canonical_index != index:
            raise ValueError("record canonical_index does not match canonical ordering")
        if record.provenance != provenance:
            raise ValueError("record provenance does not match artifact provenance")
    for record in records:
        partner = record.get("reverse_of")
        if partner is not None and (
            partner not in by_id
            or by_id[partner].get("reverse_of") != record["event_id"]
        ):
            raise ValueError("reverse_of links must be reciprocal")
        if record.get("inventory_class") == "R1:J_ring":
            if partner is None or partner not in by_id:
                raise ValueError("R1 J_ring reverse-event handle does not resolve")
            operations = record.get("junction_ops")
            if not isinstance(operations, list) or len(operations) != 1:
                raise ValueError("R1 J_ring record requires one provenance operation")
            operation = operations[0]
            required = {
                "junction_kind",
                "attacker_site_label",
                "attacked_atom_label",
                "formed_bond",
                "chosen_resonance_localisation",
                "reverse_event_handle",
                "exact_inverse_rewrite",
            }
            if not required <= set(operation):
                raise ValueError("R1 J_ring provenance is incomplete")
            kind = operation["junction_kind"]
            if kind not in R1_JUNCTION_LABELS:
                raise ValueError("R1 J_ring provenance has an unknown junction kind")
            expected_attacked = R1_JUNCTION_LABELS[kind]["attacked_atom_label"]
            if operation["attacked_atom_label"] != expected_attacked:
                raise ValueError("R1 J_ring attacked atom label is inconsistent")
            attacker = operation["attacker_site_label"]
            if attacker != "P9" or attacker == expected_attacked:
                raise ValueError("R1 J_ring attacker must be external label P9")
            if operation["formed_bond"] != {
                "label_pair": [expected_attacked, attacker],
                "order": "1.0",
            }:
                raise ValueError("R1 J_ring formed bond labels are inconsistent")
            if operation["reverse_event_handle"] != partner:
                raise ValueError("R1 J_ring reverse-event handle does not resolve")
            action = operation.get("action")
            inverse_action = "dissociate" if action == "create" else "create"
            if action not in {"create", "dissociate"} or operation[
                "exact_inverse_rewrite"
            ] != _junction_rewrite(inverse_action, expected_attacked, attacker):
                raise ValueError("R1 J_ring inverse rewrite is inconsistent")
            partner_operations = by_id[partner].get("junction_ops")
            if not isinstance(partner_operations, list) or len(partner_operations) != 1:
                raise ValueError("R1 J_ring partner provenance is incomplete")
            partner_operation = partner_operations[0]
            if partner_operation.get("action") != inverse_action:
                raise ValueError("R1 J_ring linked actions are not inverse")
    linked_records = [
        record for record in records if record.get("reverse_of") in by_id
    ]
    duplicates = duplicate_transition_groups(linked_records)
    if duplicates:
        raise ValueError(f"duplicate reciprocal transitions: {duplicates}")
    catalogue = artifact.get("short_molecule_catalogue", {})
    if catalogue.get("size") != len(catalogue.get("molecules", [])):
        raise ValueError("short-molecule catalogue size is inconsistent")
    if catalogue.get("termination_bound", -1) < catalogue.get(
        "size", 0
    ) or not catalogue.get("terminates"):
        raise ValueError("short-molecule catalogue has no valid termination bound")
    if provenance.get("family_list_sha256") != sha256_json(
        tuple(artifact.get("families", []))
    ):
        raise ValueError("family-list provenance hash mismatch")


def apply_record(
    record: EventRecord | dict[str, Any], participants: Sequence[Any]
) -> list[Any]:
    """Apply a record's atom/bond rewrite to participant molecules.

    Atom references in rewrite operations are zero-based positions in the
    concatenated participant atom order.  The returned list contains the
    connected product components, independently reconstructed from the record.
    """
    from rmgpy.molecule.molecule import Bond

    data = record.to_dict() if isinstance(record, EventRecord) else record
    atom_map = {int(key): int(value) for key, value in data.get("atom_map", {}).items()}
    if len(atom_map) != len(set(atom_map.values())):
        raise ValueError("atom_map must be a bijection")

    molecules = [_molecule(item).copy(deep=True) for item in participants]
    atoms = [atom for molecule in molecules for atom in molecule.atoms]
    heavy_atoms = [atom for atom in atoms if atom.element.number != 1]
    if atom_map and set(atom_map) != set(range(len(heavy_atoms))):
        raise ValueError("atom_map reactant indices must be canonical")
    if atom_map and set(atom_map.values()) != set(range(len(heavy_atoms))):
        raise ValueError("atom_map product indices must be canonical")
    combined = molecules[0]
    for molecule in molecules[1:]:
        combined = combined.merge(molecule)

    for operation in data.get("bond_ops", []):
        action = operation.get("action")
        if action in {"break", "form", "change"}:
            first, second = (atoms[index] for index in operation["atoms"])
            existing = first.edges.get(second)
            if action == "break":
                if existing is None:
                    raise ValueError("cannot break a missing bond")
                combined.remove_bond(existing)
            elif action == "form":
                if existing is not None:
                    raise ValueError("cannot form an existing bond")
                combined.add_bond(Bond(first, second, order=float(operation["order"])))
            else:
                if existing is None:
                    raise ValueError("cannot change a missing bond")
                existing.order = float(operation["order"])
        elif action == "set_radical":
            atoms[operation["atom"]].radical_electrons = int(operation["value"])
        elif action == "set_charge":
            atoms[operation["atom"]].charge = int(operation["value"])
        elif action == "set_lone_pairs":
            atoms[operation["atom"]].lone_pairs = int(operation["value"])
        elif action == "set_implicit_hydrogens":
            atoms[operation["atom"]].implicit_hydrogens = int(operation["value"])
        else:
            raise ValueError(f"unknown rewrite action: {action!r}")

    products = combined.split()
    try:
        for product in products:
            product.update(sort_atoms=False)
    except Exception as error:
        raise ValueError("rewrite produced an invalid molecular graph") from error
    return sorted(
        products, key=lambda molecule: molecule.to_adjacency_list(remove_h=False)
    )


def _record_mapping_root(record, mapped_root: dict[str, Any]) -> dict[str, Any]:
    """Bind an applicability root to the stored rewrite and its resonance set."""
    from rmgpy.molecule.molecule import Molecule

    attribution_unresolved = mapped_root.get("mapping_verified") is False
    data = record.to_dict() if isinstance(record, EventRecord) else record
    reactants = [
        Molecule().from_adjacency_list(graph)
        for graph in data.get("reactant_graphs", [])
    ]
    root_index = int(mapped_root["reactant_atom_index"])
    reactant_atoms = [atom for molecule in reactants for atom in molecule.atoms]
    if root_index < 0 or root_index >= len(reactant_atoms):
        if attribution_unresolved:
            return {
                "persistent_neutral_divalent_carbon": False,
                "resonance_form_count": 0,
                **mapped_root,
            }
        return {
            **mapped_root,
            "mapping_verified": False,
            "mapping_error": "mapped applicability root is outside the reactant atom map",
            "structural_inconsistency": True,
            "persistent_neutral_divalent_carbon": False,
            "resonance_form_count": 0,
        }
    root_before = reactant_atoms[root_index]
    root_side = mapped_root.get(
        "record_role", mapped_root["family_forward_role"]
    )
    if root_side not in {"reactant", "product"}:
        return {
            **mapped_root,
            "mapping_verified": False,
            "mapping_error": "mapped applicability root role must be reactant or product",
            "persistent_neutral_divalent_carbon": False,
            "resonance_form_count": 0,
        }

    touched_indices = _touched_atom_indices(data.get("bond_ops", []))
    if root_index not in touched_indices:
        return {
            **mapped_root,
            "mapping_verified": False,
            "mapping_error": "mapped recipe root is not touched by the stored rewrite",
            "persistent_neutral_divalent_carbon": False,
            "resonance_form_count": 0,
        }

    marker = "*kmc_stored_rewrite_root"
    original_label = root_before.label
    root_before.label = marker
    expected_products = data.get("product_graphs", [])
    try:
        rewritten = apply_record(data, reactants)
    except ValueError as error:
        return {
            **mapped_root,
            "mapping_verified": False,
            "mapping_error": f"stored rewrite is invalid: {error}",
            "structural_inconsistency": True,
            "product_rewrite_verified": False,
            "persistent_neutral_divalent_carbon": False,
            "resonance_form_count": 0,
        }
    rewritten_roots = [
        atom
        for molecule in rewritten
        for atom in molecule.atoms
        if atom.label == marker
    ]
    if len(rewritten_roots) != 1:
        return {
            **mapped_root,
            "mapping_verified": False,
            "mapping_error": "stored rewrite did not preserve the mapped root identity",
            "structural_inconsistency": True,
            "persistent_neutral_divalent_carbon": False,
            "resonance_form_count": 0,
        }
    rewritten_root = rewritten_roots[0]
    rewritten_root.label = original_label
    rewritten_graphs = [
        molecule.to_adjacency_list(remove_h=False) for molecule in rewritten
    ]
    rewrite_verified = _graph_lists_isomorphic(rewritten_graphs, expected_products)
    if not rewrite_verified:
        return {
            **mapped_root,
            "mapping_verified": False,
            "mapping_error": "mapped applicability root does not reproduce the stored products",
            "structural_inconsistency": True,
            "product_rewrite_verified": False,
            "persistent_neutral_divalent_carbon": False,
            "resonance_form_count": 0,
        }

    selected_atom = root_before
    selected_molecule = next(
        molecule for molecule in reactants if root_before in molecule.atoms
    )
    if root_before.element.number == 1:
        return {
            **mapped_root,
            "mapping_verified": False,
            "mapping_error": "persistent-carbene roots must be heavy atoms",
            "product_rewrite_verified": True,
            "persistent_neutral_divalent_carbon": False,
            "resonance_form_count": 0,
        }
    heavy_before = [
        atom
        for molecule in reactants
        for atom in molecule.atoms
        if atom.element.number != 1
    ]
    reactant_heavy_index = heavy_before.index(root_before)
    atom_map = {
        int(key): int(value) for key, value in data.get("atom_map", {}).items()
    }
    product_heavy_index = atom_map.get(reactant_heavy_index)
    product_molecules = [
        Molecule().from_adjacency_list(graph) for graph in expected_products
    ]
    heavy_after = [
        atom
        for molecule in product_molecules
        for atom in molecule.atoms
        if atom.element.number != 1
    ]
    if product_heavy_index is None or not 0 <= product_heavy_index < len(heavy_after):
        return {
            **mapped_root,
            "mapping_verified": False,
            "mapping_error": "stored atom map does not contain the mapped root",
            "structural_inconsistency": True,
            "product_rewrite_verified": True,
            "persistent_neutral_divalent_carbon": False,
            "resonance_form_count": 0,
        }
    mapped_product_atom = heavy_after[product_heavy_index]
    mapped_product_molecule = next(
        molecule
        for molecule in product_molecules
        if mapped_product_atom in molecule.atoms
    )
    if mapped_product_atom.element.number != root_before.element.number:
        return {
            **mapped_root,
            "mapping_verified": False,
            "mapping_error": "stored atom map changes the mapped root element",
            "structural_inconsistency": True,
            "product_rewrite_verified": True,
            "persistent_neutral_divalent_carbon": False,
            "resonance_form_count": 0,
        }
    rewritten_molecule = next(
        molecule for molecule in rewritten if rewritten_root in molecule.atoms
    )
    correspondence_verified = any(
        mapping.get(rewritten_root) is mapped_product_atom
        for mapping in rewritten_molecule.find_isomorphism(
            mapped_product_molecule, save_order=True
        )
    )
    if not correspondence_verified:
        return {
            **mapped_root,
            "mapping_verified": False,
            "mapping_error": (
                "stored atom map disagrees with the recipe-labelled rewrite root"
            ),
            "structural_inconsistency": True,
            "product_rewrite_verified": True,
            "persistent_neutral_divalent_carbon": False,
            "resonance_form_count": 0,
        }
    if root_side == "product":
        selected_atom = mapped_product_atom
        selected_molecule = mapped_product_molecule

    persistent, resonance_form_count = False, 0
    if root_index in touched_indices:
        persistent, resonance_form_count = _persistent_labeled_u2(
            selected_molecule, selected_atom
        )

    return {
        **mapped_root,
        "mapping_verified": not attribution_unresolved,
        "product_rewrite_verified": rewrite_verified,
        "resonance_form_count": resonance_form_count,
        "persistent_neutral_divalent_carbon": persistent,
    }


def _reciprocal_mapping_root(
    record, mapped_root: dict[str, Any]
) -> dict[str, Any] | None:
    """Carry a mapped reactant root into the reciprocal reactant ordering."""
    from rmgpy.molecule.molecule import Molecule

    data = record.to_dict() if isinstance(record, EventRecord) else record
    reactants = [
        Molecule().from_adjacency_list(graph)
        for graph in data.get("reactant_graphs", [])
    ]
    reactant_atoms = [atom for molecule in reactants for atom in molecule.atoms]
    root_index = int(mapped_root["reactant_atom_index"])
    if root_index < 0 or root_index >= len(reactant_atoms):
        return None
    root_before = reactant_atoms[root_index]
    if root_before.element.number == 1:
        return None
    heavy_before = [
        atom for atom in reactant_atoms if atom.element.number != 1
    ]
    atom_map = {
        int(key): int(value) for key, value in data.get("atom_map", {}).items()
    }
    product_heavy_index = atom_map.get(heavy_before.index(root_before))
    product_molecules = [
        Molecule().from_adjacency_list(graph)
        for graph in data.get("product_graphs", [])
    ]
    product_atoms = [
        atom for molecule in product_molecules for atom in molecule.atoms
    ]
    heavy_after = [
        atom for atom in product_atoms if atom.element.number != 1
    ]
    if product_heavy_index is None or not 0 <= product_heavy_index < len(heavy_after):
        return None
    record_role = mapped_root.get(
        "record_role", mapped_root["family_forward_role"]
    )
    reciprocal_role = {
        "reactant": "product",
        "product": "reactant",
    }.get(record_role, record_role)
    return {
        **mapped_root,
        "reactant_atom_index": product_atoms.index(heavy_after[product_heavy_index]),
        "record_role": reciprocal_role,
    }


def classify_persistent_carbene(
    record,
    mapped_root: dict[str, Any],
    rate_source_domain: list[dict[str, Any]] | None,
) -> dict[str, Any]:
    """Classify one mapped reacting centre against selected-rate contributors."""
    root = _record_mapping_root(record, mapped_root)
    result = {
        "policy_version": PERSISTENT_CARBENE_POLICY_VERSION,
        "mapped_root": root,
        "rate_source": copy.deepcopy(rate_source_domain),
    }
    if not root.get("mapping_verified", True):
        if root.get("structural_inconsistency"):
            return {
                **result,
                "disposition": "refused-structural-inconsistency",
                "reason": f"structural-inconsistency: {root['mapping_error']}",
            }
        return {
            **result,
            "disposition": "retained-unresolved-applicability",
            "reason": root["mapping_error"],
        }
    if not root["persistent_neutral_divalent_carbon"]:
        return {**result, "disposition": "not-applicable"}
    if root.get("mapping_caveat"):
        return {
            **result,
            "disposition": "retained-unresolved-applicability",
            "reason": root["mapping_caveat"],
        }
    if rate_source_domain is None:
        return {
            **result,
            "disposition": "retained-unresolved-applicability",
            "reason": "selected rate has no serialized calibrated contributor domain",
        }
    contributors = [
        contributor
        for contributor in rate_source_domain
        if float(
            contributor.get("source", {}).get("weight", 1.0)
            if isinstance(contributor.get("source"), dict)
            else 1.0
        ) > 0.0
    ]
    matches = [
        contributor
        for contributor in contributors
        if any(
            source_root.get("recipe_label") == root["recipe_label"]
            and source_root.get("family_forward_role")
            == root["family_forward_role"]
            and source_root.get("persistent_neutral_divalent_carbon") is True
            for source_root in contributor.get("mapped_roots") or []
        )
    ]
    if contributors and len(matches) == len(contributors):
        return {
            **result,
            "disposition": "retained-supported-transfer",
            "supporting_contributors": copy.deepcopy(matches),
            "caveat": "same reacting role; donor/acceptor environment may differ",
        }
    if matches or any(
        contributor.get("mapped_roots") is None for contributor in contributors
    ):
        return {
            **result,
            "disposition": "retained-unresolved-applicability",
            "reason": "the whole selected contributor domain is not approved",
        }
    return {
        **result,
        "disposition": "refused-unsupported-transfer",
        "reason": "no selected rate contributor supports persistent u2 in the mapped role",
    }


def apply_persistent_carbene_policy(
    pair: Sequence[Any],
    mapped_root: dict[str, Any],
    rate_source_domain: list[dict[str, Any]] | None,
) -> tuple[list[Any], dict[str, Any] | None]:
    """Apply one decision atomically before or after reciprocal pair linking."""
    if len(pair) != 2:
        raise ValueError("applicability policy requires exactly one reciprocal pair")
    records = [
        record.to_dict() if isinstance(record, EventRecord) else copy.deepcopy(record)
        for record in pair
    ]
    by_id = {record["event_id"]: record for record in records}
    prepublication = all(record.get("reverse_of") is None for record in records)
    linked = all(
        record.get("reverse_of") in by_id
        and by_id[record["reverse_of"]].get("reverse_of") == record["event_id"]
        for record in records
    )
    if len(by_id) != 2 or not (prepublication or linked):
        raise ValueError("applicability policy requires reciprocal reverse links")
    decision = classify_persistent_carbene(
        records[0], mapped_root, rate_source_domain
    )
    if decision["disposition"] != "refused-structural-inconsistency":
        reciprocal_root = _reciprocal_mapping_root(records[0], mapped_root)
        if reciprocal_root is not None:
            reciprocal_decision = classify_persistent_carbene(
                records[1], reciprocal_root, rate_source_domain
            )
            if (
                reciprocal_decision["disposition"]
                == "refused-structural-inconsistency"
            ):
                decision = reciprocal_decision
    if decision["disposition"] == "not-applicable":
        return records, None
    if decision["disposition"] in {
        "refused-unsupported-transfer",
        "refused-structural-inconsistency",
    }:
        return [], {
            **decision,
            "record_ids": [record["event_id"] for record in records],
        }
    for record in records:
        record["rate_source"] = {
            **record.get("rate_source", {}),
            "applicability": copy.deepcopy(decision),
        }
    return records, None


def _pair_relative_record(record: dict[str, Any]) -> dict[str, Any]:
    """Remove concrete link identities while retaining reciprocal semantics."""
    ignored = {"event_id", "reverse_of", "canonical_index"}
    result = {
        key: copy.deepcopy(value)
        for key, value in record.items()
        if key not in ignored
    }
    for operation in result.get("junction_ops") or []:
        if operation.get("reverse_event_handle") is not None:
            operation["reverse_event_handle"] = "__pair_partner__"
    return result


_TRANSITION_SHAPE_FIELDS = {
    "reactant_graphs",
    "bond_ops",
    "inventory_class",
    "junction_ops",
}


def _transition_shape_record(record: dict[str, Any]) -> dict[str, Any]:
    """Return a compact prefilter for potentially duplicate transitions."""
    return _pair_relative_record({
        key: value
        for key, value in record.items()
        if key in _TRANSITION_SHAPE_FIELDS
    })


def duplicate_transition_groups(records: Sequence[dict[str, Any]]) -> list[list[tuple[str, str]]]:
    """Return duplicate reciprocal transitions using pair-relative identities."""
    by_id = {record["event_id"]: record for record in records}
    seen = set()
    groups: dict[str, list[tuple[str, str]]] = {}
    for event_id in sorted(by_id):
        if event_id in seen:
            continue
        record = by_id[event_id]
        partner_id = record.get("reverse_of")
        partner = by_id.get(partner_id)
        if partner is None or partner.get("reverse_of") != event_id:
            raise ValueError("duplicate audit requires complete reciprocal pairs")
        pair_ids = tuple(sorted((event_id, partner_id)))
        seen.update(pair_ids)
        shape_views = sorted(
            canonical_json_bytes(_transition_shape_record(item)).decode("ascii")
            for item in (record, partner)
        )
        groups.setdefault(sha256_json(shape_views), []).append(pair_ids)
    return sorted(
        [sorted(group) for group in groups.values() if len(group) > 1]
    )


@dataclass(frozen=True)
class SiteProxy:
    """A centred finite proxy and its declared kMC site type."""

    site_type: str
    reactants: Sequence[Any]
    frontier: bool = False
    metadata: dict[str, Any] = field(default_factory=dict)
    participant_site_types: Sequence[str] = field(default_factory=tuple)


class EventSetCompiler:
    """Compile declared family reactions on finite, centred polymer proxies."""

    def __init__(
        self,
        kinetics_database,
        proxies: Iterable[SiteProxy],
        families: Iterable[str],
        *,
        rmgpy_path: str | Path | None = None,
        database_path: str | Path | None = None,
        temperature_grid: Iterable[float] = DEFAULT_T_GRID,
        excluded_families: dict[str, str] | None = None,
        span_radius: int = 1,
        rmgpy_sha: str | None = None,
        rmg_database_sha: str | None = None,
        thermo_database=None,
        reference_thermo_provider: ReferenceThermoProvider | None = None,
        ceiling_monomer_concentration_mol_m3: float = 1000.0,
        reaction_cache: dict[str, Sequence[Any]] | None = None,
        family_candidates: Iterable[str] = PS_FAMILY_CANDIDATES,
        kinetics_depositories: Iterable[str] = ("training",),
        use_plpsec_library: bool | None = None,
    ):
        # Constructor selection takes precedence over the environment. Invalid
        # values must fail rather than silently selecting a sensitivity arm.
        if use_plpsec_library is None:
            selection = os.environ.get("RMG_KMC_PLPSEC_LIBRARY", "1")
            if selection not in {"0", "1"}:
                raise ValueError("RMG_KMC_PLPSEC_LIBRARY must be 0 or 1")
            use_plpsec_library = selection == "1"
        self.use_plpsec_library = bool(use_plpsec_library)
        self.kinetics_database = kinetics_database
        self.proxies = tuple(copy.deepcopy(tuple(proxies)))
        self.families = tuple(sorted(set(families)))
        self.rmgpy_path = (
            Path(rmgpy_path) if rmgpy_path else Path(__file__).resolve().parents[2]
        )
        self.database_path = Path(database_path) if database_path else None
        self.temperature_grid = tuple(float(t) for t in temperature_grid)
        self.excluded_families = dict(sorted((excluded_families or {}).items()))
        self.span_radius = span_radius
        self.rmgpy_sha = rmgpy_sha
        self.rmg_database_sha = rmg_database_sha
        self.thermo_database = thermo_database
        self.rate_rule_preparation = prepare_rate_rules(
            kinetics_database, thermo_database,
            kinetics_depositories=tuple(kinetics_depositories), verbose=True,
        )
        database_commit = rmg_database_sha or _git_sha(self.database_path)
        self.thermo_assignment = getattr(
            reference_thermo_provider, "assignment", None
        ) or SharedThermoAssignment(thermo_database, database_commit)
        self.reference_thermo_provider = reference_thermo_provider or (
            GasPhaseRMGReferenceThermo(
                thermo_database,
                database_commit,
                assignment=self.thermo_assignment,
            )
        )
        self.ceiling_monomer_concentration_mol_m3 = float(
            ceiling_monomer_concentration_mol_m3
        )
        self.reaction_cache = dict(reaction_cache or {})
        self.family_candidates = tuple(sorted(set(family_candidates)))
        self._applicability_source_cache: dict[
            tuple[str, int, bool], list[dict[str, Any]]
        ] = {}
        self._applicability_refusals: list[dict[str, Any]] = []
        self._compiled_artifact: dict[str, Any] | None = None

    def _generate(self, proxy: SiteProxy):
        """Use RMG's public merged-path pipeline (the C2 degeneracy oracle)."""
        if proxy.site_type in self.reaction_cache:
            return (
                reaction
                for reaction in self.reaction_cache[proxy.site_type]
                if reaction.family in self.families
            )
        reactions = []
        for family in self.families:
            reactions.extend(
                self.kinetics_database.generate_reactions_from_families(
                    [_copy_participant(item) for item in proxy.reactants],
                    only_families=[family],
                    resonance=True,
                )
            )
        return reactions

    def _rate_table(self, reaction) -> tuple[dict[str, Any] | None, dict[str, Any]]:
        source_reaction = getattr(reaction, "source_reaction", reaction)
        kinetics = getattr(source_reaction, "kinetics", None)
        source = None
        entry = None
        if kinetics is None:
            family = self.kinetics_database.families[
                getattr(source_reaction, "family", "")
            ]
            try:
                kinetics, source, entry, _ = family.get_kinetics(
                    source_reaction,
                    template_labels=source_reaction.template,
                    degeneracy=source_reaction.degeneracy,
                    return_all_kinetics=False,
                )
            except Exception as error:
                return None, {
                    "kind": "RMG family estimate",
                    "available": False,
                    "error": str(error),
                }
        if kinetics is None:
            return None, {"kind": "RMG family estimate", "available": False}
        evaluated_kinetics = kinetics
        kinetics_conversion = None
        if isinstance(kinetics, (ArrheniusBM, ArrheniusEP)):
            try:
                species_thermo_assignments = (
                    self.thermo_assignment.assign_reaction(source_reaction)
                )
                reaction_enthalpy = float(
                    source_reaction.get_enthalpy_of_reaction(298)
                )
                model_generation_reaction = copy.deepcopy(source_reaction)
                model_generation_reaction.kinetics = copy.deepcopy(kinetics)
                model_generation_reaction.fix_barrier_height()
                evaluated_kinetics = model_generation_reaction.kinetics
                kinetics_conversion = {
                    "input_model": type(kinetics).__name__,
                    "method": "to_arrhenius(reaction.get_enthalpy_of_reaction(298))",
                    "post_conversion": "reaction.fix_barrier_height()",
                    "output_model": type(evaluated_kinetics).__name__,
                    "reaction_enthalpy_J_per_mol": reaction_enthalpy,
                    "species_thermo_assignments": species_thermo_assignments,
                    "activation_energy_J_per_mol": float(
                        evaluated_kinetics.Ea.value_si
                    ),
                }
            except Exception as error:
                return None, {
                    "kind": "RMG family estimate",
                    "available": False,
                    "error": f"enthalpy-dependent kinetics conversion failed: {error}",
                }
        values = [
            float(evaluated_kinetics.get_rate_coefficient(t))
            for t in self.temperature_grid
        ]
        reference = "forward RMG family estimate"
        equilibrium_constants = None
        if source_reaction is not reaction:
            if self.thermo_database is None:
                return None, {
                    "kind": "RMG family estimate",
                    "available": False,
                    "error": "reverse rate requires gas-phase thermochemistry",
                }
            try:
                self.thermo_assignment.assign_reaction(source_reaction)
                equilibrium_constants = [
                    float(
                        source_reaction.get_equilibrium_constant(temperature, type="Kc")
                    )
                    for temperature in self.temperature_grid
                ]
                values = [
                    value / equilibrium_constant
                    for value, equilibrium_constant in zip(
                        values, equilibrium_constants
                    )
                ]
                reference = "reverse from RMG gas-phase Kc"
            except Exception as error:
                return None, {
                    "kind": "RMG gas-phase reverse",
                    "available": False,
                    "error": str(error),
                }
        units = get_rate_coefficient_units_from_reaction_order(len(reaction.reactants))
        if any(not math.isfinite(value) or value <= 0 for value in values):
            return None, {
                "kind": "RMG family estimate",
                "available": False,
                "error": "non-finite or non-positive rate coefficient",
            }
        return (
            {
                "T": list(self.temperature_grid),
                "k": values,
                "interpolation": "linear-ln-k",
                "extrapolation": "refuse",
            },
            {
                "kind": "RMG family estimate",
                "available": True,
                "units": units,
                "reference": reference,
                "source": str(source) if source is not None else "attached",
                "entry": str(entry) if entry else None,
                "rank": getattr(entry, "rank", None),
                "uncertainty": str(getattr(kinetics, "uncertainty", None) or "UNKNOWN"),
                **(
                    {"kinetics_conversion": kinetics_conversion}
                    if kinetics_conversion is not None
                    else {}
                ),
                "source_temperature_range_K": [
                    getattr(getattr(kinetics, bound, None), "value_si", None)
                    for bound in ("Tmin", "Tmax")
                ],
                "equilibrium_constant_table": (
                    {"T": list(self.temperature_grid), "Kc": equilibrium_constants}
                    if equilibrium_constants is not None
                    else None
                ),
                "lifetime_s": (
                    [1.0 / value if value > 0 else None for value in values]
                    if source_reaction is not reaction and len(reaction.reactants) == 1
                    else None
                ),
            },
        )

    def _record(
        self,
        proxy: SiteProxy,
        reaction,
        provenance: dict[str, Any],
        orientation: str = "as_generated",
        rate_override=None,
    ) -> EventRecord:
        # The public pipeline returns Species.  atom_map accepts molecule lists.
        graph_reaction = SimpleNamespace(
            reactants=[_molecule(item) for item in reaction.reactants],
            products=[_molecule(item) for item in reaction.products],
        )
        mapping = extract_atom_map(graph_reaction)
        reactant_formula, product_formula = (
            mapping["reactant_element_counts"],
            mapping["product_element_counts"],
        )
        formula_delta = {
            key: product_formula.get(key, 0) - reactant_formula.get(key, 0)
            for key in sorted(set(reactant_formula) | set(product_formula))
        }
        family = getattr(reaction, "family", "")
        participant_site_types, reactant_multiplicities = _participant_inventory(
            proxy, reaction
        )
        inventory_class = _inventory_class(reaction)
        if (
            proxy.metadata.get("generic_reference_pair")
            and inventory_class == "R1:J_ring"
        ):
            inventory_class = None
        is_ring = inventory_class == "R1:J_ring"
        template = _template(reaction)
        arity = len(reaction.reactants)
        degeneracy = float(getattr(reaction, "degeneracy", 1.0))
        atom_map = _canonical_atom_map(reaction)
        bond_ops = _bond_ops(reaction)
        radical_delta = _radicals(reaction.products) - _radicals(reaction.reactants)
        implicit_h_delta = _implicit_h(reaction.products) - _implicit_h(
            reaction.reactants
        )
        reactant_graphs = _graph_adjacencies(reaction.reactants)
        product_graphs = _graph_adjacencies(reaction.products)
        coproducts = [{"formula": _formula([p])} for p in reaction.products[1:]]
        junction_ops = _junction_ops(inventory_class, proxy.metadata)
        pair_convention = _pair_convention(reaction.reactants)
        ssa_multiplier = ssa_multiplier_for(list(reaction.reactants))
        resonance_derived = _is_resonance_derived(reaction)
        # Kinetics estimation is allowed to normalize resonance forms in-place.
        # Capture the generated structural oracle first so the record always
        # represents the public pipeline reaction that actually fired.
        k_table, rate_source = (
            rate_override if rate_override is not None else self._rate_table(reaction)
        )
        status, reason = (
            ("refused", "frontier channel") if proxy.frontier else ("enabled", "")
        )
        if not proxy.frontier and not getattr(reaction, "reversible", True):
            status, reason = "irreversible", "RMG reaction is declared irreversible"
        if not rate_source["available"]:
            status, reason = "refused", "RMG family kinetics unavailable"
        return EventRecord(
            family=family,
            template=template,
            arity=arity,
            participant_site_types=participant_site_types,
            reactant_multiplicities=reactant_multiplicities,
            raw_path_degeneracy=degeneracy,
            degeneracy=degeneracy,
            ssa_multiplier=ssa_multiplier,
            rate_order=arity,
            rate_units=rate_source.get("units", ""),
            k_table=k_table,
            atom_map=atom_map,
            bond_ops=bond_ops,
            element_delta=formula_delta,
            formula_delta=formula_delta,
            radical_delta=radical_delta,
            implicit_h_delta=implicit_h_delta,
            mass_delta=sum(
                ATOMIC_MASS_NUMBERS.get(k, 0) * v for k, v in formula_delta.items()
            ),
            reactant_pair_convention=pair_convention,
            resonance_derived=resonance_derived,
            status=status,
            status_reason=reason,
            orientation=orientation,
            rate_source=rate_source,
            site_type=proxy.site_type,
            coproducts=coproducts,
            inventory_class=inventory_class,
            thermo_provenance=(
                {
                    "gas_phase": "RMG Kc reference",
                    "condensed_phase_constraint": "UNKNOWN",
                }
                if is_ring
                else self.reference_thermo_provider.provenance
            ),
            provenance=provenance,
            reactant_graphs=reactant_graphs,
            product_graphs=product_graphs,
            feature_ops=[
                operation
                for operation in bond_ops
                if operation["action"].startswith("set_")
            ],
            junction_ops=junction_ops,
        )

    def _direction_proxy(self, proxy, reaction):
        if _same_participants(proxy.reactants, reaction.reactants):
            return proxy
        labels = []
        for participant in reaction.reactants:
            radical_count = _radicals([participant])
            molecule = _molecule(participant)
            radicals = [atom for atom in molecule.atoms if atom.radical_electrons]
            if radical_count > 1:
                label = "doubly_featured"
            elif (
                radical_count
                and hasattr(molecule, "is_atom_in_cycle")
                and molecule.is_atom_in_cycle(radicals[0])
            ):
                label = "junction_radical"
            elif radical_count:
                carbons = [
                    neighbor for neighbor in getattr(radicals[0], "edges", {})
                    if neighbor.element.number == 6
                ]
                ring_neighbors = [
                    neighbor for neighbor in carbons
                    if hasattr(molecule, "is_atom_in_cycle")
                    and molecule.is_atom_in_cycle(neighbor)
                ]
                backbone_degree = len(carbons) - len(ring_neighbors)
                label = (
                    "interior_radical" if backbone_degree >= 2
                    else "benzylic_end_radical" if ring_neighbors
                    else "end_radical"
                )
            else:
                label = (
                    "styrene"
                    if hasattr(molecule, "to_smiles")
                    and molecule.to_smiles() == "C=Cc1ccccc1"
                    else "pristine"
                )
            labels.append(label)
        return replace(
            proxy, reactants=reaction.reactants, participant_site_types=tuple(labels)
        )

    def _linked_family_pair(self, proxy, reaction, provenance):
        """Reuse immutable structural representatives only within this pair."""
        token = _PAIR_MOLECULE_CACHE.set({})
        try:
            return self._build_linked_family_pair(
                proxy, copy.deepcopy(reaction), provenance
            )
        finally:
            _PAIR_MOLECULE_CACHE.reset(token)

    def _build_linked_family_pair(self, proxy, reaction, provenance):
        """Estimate exactly one direction, then invert its reference-state Kc."""
        source = getattr(reaction, "source_reaction", reaction)
        proxy = replace(
            proxy, metadata={**proxy.metadata, "generic_reference_pair": True}
        )
        structural_forward = self._record(
            self._direction_proxy(proxy, source),
            source,
            provenance,
            rate_override=(None, {"available": True}),
        )
        estimate = copy.deepcopy(source)
        kinetics = getattr(estimate, "kinetics", None)
        estimated_forward = True
        source_name, entry = "attached", None
        if kinetics is None:
            family = self.kinetics_database.families[estimate.family]
            try:
                kinetics, source_name, entry, estimated_forward = family.get_kinetics(
                    estimate,
                    template_labels=estimate.template,
                    degeneracy=estimate.degeneracy,
                    return_all_kinetics=False,
                )
            except Exception as error:
                refused = replace(
                    structural_forward,
                    status="refused",
                    status_reason=f"RMG family kinetics unavailable: {error}",
                    event_id="",
                )
                return [refused], refused
        estimate.kinetics = kinetics
        if not estimated_forward:
            estimate.reactants, estimate.products = (
                estimate.products,
                estimate.reactants,
            )
        family = getattr(self.kinetics_database, "families", {}).get(estimate.family)
        reverse_view = _reverse_view(source)
        structural_reverse = self._record(
            self._direction_proxy(proxy, reverse_view),
            reverse_view,
            provenance,
            "reversed",
            rate_override=(None, {"available": True}),
        )
        forward, reverse = (
            (structural_forward, structural_reverse)
            if estimated_forward
            else (structural_reverse, structural_forward)
        )
        applicability = None
        has_selected_u2 = any(
            " u2 " in graph
            for graph in forward.reactant_graphs + forward.product_graphs
        )
        if family is not None and has_selected_u2:
            try:
                mapped_roots = _mapped_reaction_u2_roots(
                    family, _reaction_from_record(forward), forward
                )
            except ActionError:
                if not getattr(family, "auto_generated", False):
                    raise
                mapped_roots = _stored_u2_root_fallback(forward)
            if mapped_roots:
                rate_source_domain = _training_source_domain(
                    family, estimate, self._applicability_source_cache
                )
                policy_results = [
                    apply_persistent_carbene_policy(
                        (forward, reverse), mapped_root, rate_source_domain
                    )
                    for mapped_root in mapped_roots
                ]
                refused = next(
                    (
                        refusal for _, refusal in policy_results
                        if refusal is not None
                    ),
                    None,
                )
                if refused is not None:
                    refusal = {
                        **refused,
                        "family": estimate.family,
                        "template": _template(estimate),
                    }
                    self._applicability_refusals.append(refusal)
                    return [], None
                applicability = next(
                    (
                        published[0]["rate_source"].get("applicability")
                        for published, _ in policy_results
                        if published[0]["rate_source"].get("applicability")
                        is not None
                    ),
                    None,
                )

        # Applicability is decided before evaluating kinetics or reference thermo.
        forward_table, rate_source = self._rate_table(estimate)
        rate_source = {
            **rate_source,
            "source": str(source_name),
            "entry": str(entry) if entry else None,
            "rank": getattr(entry, "rank", None),
            "generated_is_forward": getattr(source, "is_forward", None),
            "family_template_direction": "forward" if estimated_forward else "reverse",
            "reference_thermo": self.reference_thermo_provider.provenance,
            **(
                {"applicability": applicability}
                if applicability is not None
                else {}
            ),
        }
        if family is not None and getattr(family, "auto_generated", True) is False:
            rate_source.update(_rate_rule_source(family, estimate))
        forward = replace(
            forward,
            k_table=forward_table,
            rate_source=rate_source,
            rate_units=rate_source.get("units", ""),
            event_id="",
        )
        initiated_forward = _same_participants(proxy.reactants, estimate.reactants)
        if forward_table is None:
            refused = replace(
                forward,
                status="refused",
                status_reason="RMG family kinetics unavailable",
                event_id="",
            )
            return [refused], refused
        if not getattr(source, "reversible", True):
            irreversible = replace(
                forward,
                status="refused" if proxy.frontier else "irreversible",
                status_reason=(
                    "frontier channel"
                    if proxy.frontier
                    else "RMG reaction is declared irreversible"
                ),
                event_id="",
            )
            return [irreversible], irreversible
        reference_result = None
        try:
            if hasattr(self.reference_thermo_provider, "evaluate"):
                reference_result = self.reference_thermo_provider.evaluate(
                    estimate, self.temperature_grid
                )
                constants = list(reference_result.equilibrium_constants)
            else:
                constants = self.reference_thermo_provider.equilibrium_constants(
                    estimate, self.temperature_grid
                )
        except (ThermoUnavailable, AttributeError) as error:
            irreversible = replace(
                forward,
                status="refused" if proxy.frontier else "irreversible",
                status_reason="frontier channel" if proxy.frontier else str(error),
                event_id="",
            )
            return [irreversible], irreversible
        if self.use_plpsec_library and (
            matches_head_to_tail(forward) or matches_head_to_tail(reverse)
        ):
            replaced_table = copy.deepcopy(forward_table)
            replaced_source = copy.deepcopy(rate_source)
            # RMG may estimate the unimolecular family direction. Normalize to
            # physical propagation before installing the library coefficient.
            if not matches_head_to_tail(forward):
                forward, reverse = reverse, forward
                if reference_result is not None:
                    reference_result = reference_result.reversed_direction()
                    constants = list(reference_result.equilibrium_constants)
                else:
                    constants = [1.0 / constant for constant in constants]
                replaced_table["k"] = [
                    rate * constant for rate, constant in zip(forward_table["k"], constants)
                ]
                initiated_forward = not initiated_forward
            library_entry = load_plpsec_entry()
            rate_source = {
                "kind": "kMC kinetics library",
                "available": True,
                "units": library_entry["output_units"],
                "library_entry": library_entry["entry_id"],
                "citation": library_entry["citation"],
                "entry": library_entry,
                "replaced_rmg_estimate": {
                    "template": forward.template,
                    "source": replaced_source,
                    "k_table": copy.deepcopy(forward_table),
                    "units": replaced_source["units"],
                    "propagation_k_table": replaced_table,
                    "propagation_units": library_entry["output_units"],
                },
                "reference_thermo": self.reference_thermo_provider.provenance,
            }
            forward_table = plpsec_rate_table(self.temperature_grid)
            forward = replace(
                forward, k_table=forward_table, rate_source=rate_source,
                rate_units=library_entry["output_units"], event_id="",
            )
        reverse_table = {
            **forward_table,
            "k": [
                rate / constant for rate, constant in zip(forward_table["k"], constants)
            ],
        }
        if any(not math.isfinite(rate) or rate <= 0 for rate in reverse_table["k"]):
            raise ValueError(
                "reference-thermo reverse rates must be positive and finite"
            )
        thermo = {
            **self.reference_thermo_provider.provenance,
            **(
                {
                    "reaction_enthalpy_298_J_per_mol": (
                        reference_result.reaction_enthalpy_298_J_per_mol
                    ),
                    "species_thermo_assignments": (
                        list(reference_result.species_thermo_assignments)
                    ),
                }
                if reference_result is not None
                else {}
            ),
            "equilibrium_constant_table": {
                "T": list(self.temperature_grid),
                "Kc": constants,
            },
        }
        forward = replace(
            forward, thermo_provenance=thermo, reverse_of=LINK_PLACEHOLDER, event_id=""
        )
        reverse = replace(
            reverse,
            k_table=reverse_table,
            rate_source={
                "kind": "reference-thermo reverse",
                "available": True,
                "units": reverse.rate_units
                or get_rate_coefficient_units_from_reaction_order(reverse.arity),
                "reference": ("k_library / Kc" if rate_source["kind"] == "kMC kinetics library"
                              else "k_family / Kc"),
                **self.reference_thermo_provider.provenance,
                **(
                    {"forward_kinetics_conversion": rate_source["kinetics_conversion"]}
                    if "kinetics_conversion" in rate_source
                    else {}
                ),
                **({"forward_rate_source": rate_source}
                   if "comment" in rate_source or rate_source["kind"] == "kMC kinetics library" else {}),
                **(
                    {"applicability": applicability}
                    if applicability is not None
                    else {}
                ),
            },
            rate_units=get_rate_coefficient_units_from_reaction_order(reverse.arity),
            thermo_provenance=thermo,
            reverse_of=forward.event_id,
            event_id="",
        )
        forward = replace(forward, reverse_of=reverse.event_id, event_id="")
        reverse = replace(reverse, reverse_of=forward.event_id, event_id="")
        return [forward, reverse], forward if initiated_forward else reverse

    def _linked_ortho_pair(
        self,
        proxy: SiteProxy,
        reaction,
        provenance: dict[str, Any],
    ) -> tuple[EventRecord, EventRecord]:
        """Compile one half-rate ortho capture and its exact dissociation."""
        forward = _partition_degeneracy(self._record(proxy, reaction, provenance), 2)
        if forward.inventory_class != "R1:J_ring":
            raise ValueError("archived J_ortho reaction is not a ring capture")
        forward = replace(forward, reverse_of=LINK_PLACEHOLDER, event_id="")
        reverse_proxy = SiteProxy(
            "J_ring",
            reaction.products,
            participant_site_types=("J_ring",),
            metadata={
                "inventory_class": "R1",
                "reference_thermo": "RMG gas-phase Kc",
                "junction_kind": proxy.metadata["junction_kind"],
                "attacker_site_label": "P9",
            },
        )
        reverse = self._record(
            reverse_proxy,
            _reverse_view(reaction),
            provenance,
            "reversed",
        )
        reverse = replace(
            reverse,
            inventory_class="R1:J_ring",
            junction_ops=_junction_ops(
                "R1:J_ring",
                reverse_proxy.metadata,
                action="dissociate",
                reverse_event_handle=forward.event_id,
            ),
            status="enabled" if reverse.k_table is not None else reverse.status,
            status_reason="" if reverse.k_table is not None else reverse.status_reason,
            thermo_provenance={
                "gas_phase": "RMG Kc reference",
                "condensed_phase_constraint": "UNKNOWN",
            },
            reverse_of=forward.event_id,
            event_id="",
        )
        forward = replace(
            forward,
            junction_ops=_junction_ops(
                "R1:J_ring",
                proxy.metadata,
                action="create",
                reverse_event_handle=reverse.event_id,
            ),
            reverse_of=reverse.event_id,
            event_id="",
        )
        reverse = replace(
            reverse,
            junction_ops=_junction_ops(
                "R1:J_ring",
                reverse_proxy.metadata,
                action="dissociate",
                reverse_event_handle=forward.event_id,
            ),
            reverse_of=forward.event_id,
            event_id="",
        )
        return forward, reverse

    def _linked_para_pair(
        self,
        proxy: SiteProxy,
        reaction,
        provenance: dict[str, Any],
    ) -> tuple[EventRecord, EventRecord]:
        """Compile the archived para capture and its exact dissociation."""
        forward = self._record(proxy, reaction, provenance)
        if forward.inventory_class != "R1:J_ring":
            raise ValueError("archived J_para reaction is not a ring capture")
        forward_source = {
            **forward.rate_source,
            **ARCHIVED_J_PARA_RATE_PROVENANCE,
            "kind": "R-009 archived RMG rate rule",
            "source": ARCHIVED_J_PARA_RATE_PROVENANCE["archive_record"],
            "entry": ARCHIVED_J_PARA_RULE,
            "rank": ARCHIVED_J_PARA_RATE_PROVENANCE["rule_rank"],
        }
        forward = replace(
            forward,
            rate_source=forward_source,
            reverse_of=LINK_PLACEHOLDER,
            event_id="",
        )
        reverse_proxy = SiteProxy(
            "J_ring",
            reaction.products,
            participant_site_types=("J_ring",),
            metadata={
                "inventory_class": "R1",
                "reference_thermo": "RMG gas-phase Kc",
                "junction_kind": "J_para",
                "attacker_site_label": "P9",
            },
        )
        reverse = self._record(
            reverse_proxy,
            _reverse_view(reaction),
            provenance,
            "reversed",
        )
        reverse_source = {
            **reverse.rate_source,
            **ARCHIVED_J_PARA_RATE_PROVENANCE,
            "kind": "R-009 archived RMG gas-phase reverse",
            "reference": "reverse from archived RMG gas-phase Kc",
            "source": ARCHIVED_J_PARA_RATE_PROVENANCE["archive_record"],
            "entry": ARCHIVED_J_PARA_RULE,
            "rank": ARCHIVED_J_PARA_RATE_PROVENANCE["rule_rank"],
            "direction": "reverse through Kc",
        }
        reverse = replace(
            reverse,
            rate_source=reverse_source,
            inventory_class="R1:J_ring",
            junction_ops=_junction_ops(
                "R1:J_ring",
                reverse_proxy.metadata,
                action="dissociate",
                reverse_event_handle=forward.event_id,
            ),
            status="enabled" if reverse.k_table is not None else reverse.status,
            status_reason="" if reverse.k_table is not None else reverse.status_reason,
            thermo_provenance={
                "gas_phase": "R-009 archived RMG NASA Kc reference",
                "condensed_phase_constraint": "UNKNOWN",
            },
            reverse_of=forward.event_id,
            event_id="",
        )
        forward = replace(
            forward,
            junction_ops=_junction_ops(
                "R1:J_ring",
                proxy.metadata,
                action="create",
                reverse_event_handle=reverse.event_id,
            ),
            thermo_provenance={
                "gas_phase": "R-009 archived RMG NASA Kc reference",
                "condensed_phase_constraint": "UNKNOWN",
            },
            reverse_of=reverse.event_id,
            event_id="",
        )
        reverse = replace(
            reverse,
            junction_ops=_junction_ops(
                "R1:J_ring",
                reverse_proxy.metadata,
                action="dissociate",
                reverse_event_handle=forward.event_id,
            ),
            reverse_of=forward.event_id,
            event_id="",
        )
        return forward, reverse

    def _para_junction_records(self, provenance: dict[str, Any]) -> list[EventRecord]:
        """Compile the archived g=1 para capture and exact reverse."""
        if "R_Recombination" not in self.families:
            return []
        proxy = _para_junction_proxy()
        return list(
            self._linked_para_pair(proxy, _archived_para_reaction(proxy), provenance)
        )

    def _ortho_junction_records(self, provenance: dict[str, Any]) -> list[EventRecord]:
        """Compile the archived g=2 ortho channel as S6/S7 half-rate records."""
        if "R_Recombination" not in self.families:
            return []
        proxy = _ortho_junction_proxy()
        reaction = _archived_ortho_reaction(proxy)

        records = []
        for junction_kind in ("J_ortho_S7", "J_ortho_S6"):
            localised_proxy = replace(
                proxy,
                metadata={**proxy.metadata, "junction_kind": junction_kind},
            )
            records.extend(
                self._linked_ortho_pair(
                    localised_proxy,
                    (
                        reaction
                        if junction_kind == "J_ortho_S7"
                        else _reflect_ortho_reaction(reaction)
                    ),
                    provenance,
                )
            )
        return records

    @classmethod
    def discover_family_reactions(
        cls,
        kinetics_database,
        proxies: Iterable[SiteProxy],
        candidate_families: Iterable[str] = PS_FAMILY_CANDIDATES,
        family_universe: Iterable[str] | None = None,
    ) -> tuple[list[str], dict[str, str], dict[str, Sequence[Any]]]:
        """Account for every database family, generating only the loaded PS filter.

        A family is active only if the public merged-path pipeline produces a
        reaction on the bounded proxy set.  Loaded families outside the explicit
        PS radical-growth filter remain visible with a stable exclusion reason;
        callers may cheaply enumerate that full universe from the database tree.
        """
        loaded = set(kinetics_database.families)
        universe = set(family_universe) if family_universe is not None else loaded
        declared_candidates = set(candidate_families)
        missing_candidates = (universe & declared_candidates) - loaded
        if missing_candidates:
            raise ValueError(
                "PS family candidates were enumerated but not loaded: "
                + ", ".join(sorted(missing_candidates))
            )
        candidates = sorted(loaded & declared_candidates)
        active_set: set[str] = set()
        scheduled_set: set[str] = set()
        reaction_cache = {}
        for proxy in proxies:
            logging.getLogger(__name__).info(
                "generating all candidate families on %s", proxy.site_type
            )
            proxy_candidates = sorted(
                set(proxy.metadata.get("family_candidates", candidates))
                & set(candidates)
            )
            scheduled_set.update(proxy_candidates)
            reactions = []
            for family in proxy_candidates:
                logging.getLogger(__name__).info(
                    "generating %s on %s", family, proxy.site_type
                )
                generated = kinetics_database.generate_reactions_from_families(
                    [_copy_participant(item) for item in proxy.reactants],
                    only_families=[family],
                    resonance=True,
                )
                reactions.extend(generated)
                logging.getLogger(__name__).info(
                    "generated %d %s reactions on %s",
                    len(generated),
                    family,
                    proxy.site_type,
                )
            proxy.metadata["generated_families"] = sorted(
                {reaction.family for reaction in reactions}
            )
            reaction_cache[proxy.site_type] = reactions
            logging.getLogger(__name__).info(
                "generated %d reactions on %s", len(reactions), proxy.site_type
            )
            active_set.update(reaction.family for reaction in reactions)
        active = sorted(active_set)
        excluded = {}
        for family in sorted(universe - active_set):
            excluded[family] = (
                PS_FAMILY_FILTER_REASON
                if family not in candidates
                else (
                    "no compatible bounded L=3 proxy declaration in the "
                    "PS family filter"
                    if family not in scheduled_set
                    else "no reaction on bounded L=3 PS proxy set with "
                    "resonance enabled"
                )
            )
        return active, excluded, reaction_cache

    @classmethod
    def discover_families(
        cls,
        kinetics_database,
        proxies: Iterable[SiteProxy],
        candidate_families: Iterable[str] = PS_FAMILY_CANDIDATES,
        family_universe: Iterable[str] | None = None,
    ) -> tuple[list[str], dict[str, str]]:
        """Return firing filtered families plus reasons for every loaded exclusion."""
        active, excluded, _ = cls.discover_family_reactions(
            kinetics_database, proxies, candidate_families, family_universe
        )
        return active, excluded

    def compile(self) -> dict[str, Any]:
        if self._compiled_artifact is not None:
            return copy.deepcopy(self._compiled_artifact)
        self._applicability_refusals = []

        proxy_inputs = [
            _proxy_fingerprint(p)
            for p in sorted(self.proxies, key=lambda p: p.site_type)
        ]
        provenance = {
            "compiler_sha256": _LOADED_COMPILER_HASH,
            "compiler_sources_sha256": compiler_source_hash(),
            "family_list_sha256": sha256_json(self.families),
            "family_filter_sha256": sha256_json(self.family_candidates),
            "rmgpy_sha": self.rmgpy_sha or _git_sha(self.rmgpy_path),
            "rmg_database_sha": self.rmg_database_sha or _git_sha(self.database_path),
            "proxy_set_sha256": sha256_json(proxy_inputs),
            "kinetics_libraries": {
                "styrene_plpsec": {
                    "enabled": self.use_plpsec_library,
                    "entry": load_plpsec_entry(),
                    "entry_sha256": sha256_json(load_plpsec_entry()),
                    "precedence": "exact forward chemical rewrite replaces RMG family estimate; inverse remains k_forward/Kc",
                },
            },
            "rate_rule_preparation": {
                **self.rate_rule_preparation,
                "database_sha": self.rmg_database_sha or _git_sha(self.database_path),
            },
            "applicability_policy_version": PERSISTENT_CARBENE_POLICY_VERSION,
        }
        if "R_Recombination" in self.families:
            provenance["archived_j_para_rate"] = copy.deepcopy(
                ARCHIVED_J_PARA_RATE_PROVENANCE
            )
        records = []
        discovery = []
        excluded_channels = []
        pairs = {}
        for proxy in sorted(self.proxies, key=lambda p: p.site_type):
            logging.getLogger(__name__).info("compiling pairs on %s", proxy.site_type)
            candidates = []
            for reaction in self._generate(proxy):
                reactant_side = tuple(sorted(_graph_adjacencies(reaction.reactants)))
                product_side = tuple(sorted(_graph_adjacencies(reaction.products)))
                candidates.append(
                    (
                        (
                            reaction.family,
                            _template(reaction),
                            reactant_side,
                            product_side,
                            float(reaction.degeneracy),
                        ),
                        reaction,
                        tuple(sorted((reactant_side, product_side))),
                    )
                )
            for _, reaction, sides in sorted(candidates, key=lambda item: item[0]):
                oriented, _ = _orient_to_proxy(proxy, reaction)
                coverage_site = proxy.metadata.get(
                    "coverage_site_type", proxy.site_type
                )
                if (
                    _inventory_class(oriented) == "R1:J_ring"
                    and coverage_site == "junction_radical+end_radical"
                ):
                    excluded_channels.append(
                        {
                            "site_type": coverage_site,
                            "family": reaction.family,
                            "template": _template(reaction),
                            "reason": "bounded J_ring capture replaced by immutable R-009 archived para/ortho pairs",
                        }
                    )
                    continue
                pair_key = (reaction.family, proxy.frontier, sides)
                if pair_key not in pairs:
                    pair_records, initiating = self._linked_family_pair(
                        proxy, reaction, provenance
                    )
                    pairs[pair_key] = pair_records
                    records.extend(pair_records)
                else:
                    pair_records = pairs[pair_key]
                    initiating = next(
                        (
                            record
                            for record in pair_records
                            if _graph_lists_isomorphic(
                                record.reactant_graphs,
                                _graph_adjacencies(oriented.reactants),
                            )
                        ),
                        pair_records[0] if pair_records else None,
                    )
                discovery.append(
                    {
                        "site_type": coverage_site,
                        "proxy_site_type": proxy.site_type,
                        "proxy_units": proxy.metadata.get("proxy_units"),
                        "family": reaction.family,
                        "template": _template(reaction),
                        "raw_path_degeneracy": float(reaction.degeneracy),
                        "event_id": initiating.event_id if initiating else None,
                    }
                )
        records.extend(self._para_junction_records(provenance))
        records.sort(key=lambda record: record.event_id)
        records = [
            replace(record, canonical_index=index)
            for index, record in enumerate(records)
        ]
        ortho_records = sorted(
            self._ortho_junction_records(provenance), key=lambda record: record.event_id
        )
        first_ortho_index = len(records)
        records.extend(
            replace(record, canonical_index=first_ortho_index + index)
            for index, record in enumerate(ortho_records)
        )
        artifact = {
            "schema_version": SCHEMA_VERSION,
            "families": list(self.families),
            "excluded_families": self.excluded_families,
            "provenance": provenance,
            "inventory_arms": {
                "selected": None,
                "R0": "resonance-derived ring-site capture omitted",
                "R1": "full resonance with reversible J_ring",
            },
            "inputs": {
                "families": list(self.families),
                "family_filter": {
                    "candidates": list(self.family_candidates),
                    "excluded_loaded_family_reason": PS_FAMILY_FILTER_REASON,
                },
                "proxies": proxy_inputs,
                "temperature_grid": list(self.temperature_grid),
                "span_radius": self.span_radius,
            },
            "records": [record.to_dict() for record in records],
            "discovery": discovery,
            "excluded_channels": excluded_channels,
            "applicability_refusals": copy.deepcopy(
                self._applicability_refusals
            ),
            "irreversible_pairs": [
                {
                    "event_id": r.event_id,
                    "family": r.family,
                    "reactant_graphs_sha256": sha256_json(r.reactant_graphs),
                    "product_graphs_sha256": sha256_json(r.product_graphs),
                    "reason": r.status_reason,
                }
                for r in records
                if r.status == "irreversible"
            ],
            "short_molecule_catalogue": short_ps_molecule_catalogue(self.span_radius),
            "ps_ceiling_temperature_K": None,
        }
        by_id = {record.event_id: record for record in records}
        pair_fields = {
            "benzylic_end_radical": "ps_ceiling_pairs",
            "end_radical": "ps_primary_end_ceiling_pairs",
        }
        from rmgpy.molecule.molecule import Molecule

        end_contexts = {label: [] for label in pair_fields}
        for proxy in self.proxies:
            label = proxy.metadata.get("coverage_site_type", proxy.site_type)
            units = proxy.metadata.get("proxy_units")
            if label not in end_contexts or units is None:
                continue
            products = _graph_adjacencies(proxy.reactants)
            reactants = (
                _graph_adjacencies([
                    Molecule(smiles=benzylic_ps_end_smiles(units - 1)),
                    Molecule(smiles="C=Cc1ccccc1"),
                ]) if label == "benzylic_end_radical" else None
            )
            end_contexts[label].append((units, products, reactants, _graph_list_key(products)))
        for end_label, field in pair_fields.items():
            pairs = []
            for prop in records:
                if (prop.family != "R_Addition_MultipleBond" or prop.arity != 2
                        or not prop.k_table or prop.status != "enabled"):
                    continue
                if "C=Cc1ccccc1" not in _graph_list_key(prop.reactant_graphs):
                    continue
                product_key = _graph_list_key(prop.product_graphs)
                product_units = next((
                    units for units, products, reactants, key in end_contexts[end_label]
                    if product_key == key
                    and _graph_lists_isomorphic(prop.product_graphs, products)
                    and (reactants is None or _graph_lists_isomorphic(prop.reactant_graphs, reactants))
                ), None)
                if product_units is None:
                    continue
                dep = by_id.get(prop.reverse_of)
                if (dep is None or dep.reverse_of != prop.event_id or dep.arity != 1
                        or not dep.k_table or dep.status != "enabled"
                        or not _graph_lists_isomorphic(prop.reactant_graphs, dep.product_graphs)
                        or not _graph_lists_isomorphic(prop.product_graphs, dep.reactant_graphs)):
                    continue
                temperature = ceiling_temperature(
                    prop.to_dict(), dep.to_dict(),
                    self.ceiling_monomer_concentration_mol_m3,
                )
                if temperature is None and end_label == "end_radical":
                    continue
                pairs.append({
                    "propagation_event_id": prop.event_id,
                    "depropagation_event_id": dep.event_id,
                    "monomer_concentration_mol_m3": self.ceiling_monomer_concentration_mol_m3,
                    "temperature_K": temperature,
                    **({"reactant_repeat_units": product_units - 1,
                        "product_repeat_units": product_units}
                       if end_label == "benzylic_end_radical" else {}),
                })
            pairs.sort(key=(
                (lambda pair: (pair["product_repeat_units"], pair["propagation_event_id"]))
                if end_label == "benzylic_end_radical" else
                (lambda pair: (pair["temperature_K"], pair["propagation_event_id"]))
            ))
            artifact[field] = pairs
        # Explicit chemical anchor: n=2→3 benzylic growth, never the minimum
        # ceiling across chain lengths or unrelated primary-end channels.
        anchor = next((pair for pair in artifact["ps_ceiling_pairs"]
                       if pair["product_repeat_units"] == PS_PROXY_UNITS), None)
        artifact["ps_ceiling_anchor_event_id"] = (
            anchor["propagation_event_id"] if anchor else None
        )
        artifact["ps_ceiling_temperature_K"] = anchor["temperature_K"] if anchor else None
        validate_artifact(artifact)
        self._compiled_artifact = copy.deepcopy(artifact)
        return copy.deepcopy(self._compiled_artifact)

    def write_artifact(self, directory: str | Path) -> tuple[Path, dict[str, Any]]:
        artifact = self.compile()
        payload = canonical_json_bytes(artifact)
        artifact_id = hashlib.sha256(payload).hexdigest()
        target = Path(directory) / f"{artifact_id}.json"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(payload)
        return target, artifact
