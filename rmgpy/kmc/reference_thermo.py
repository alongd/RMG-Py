"""Reference-state equilibrium constants for compiled reversible event pairs."""

from __future__ import annotations

import copy
import hashlib
import json
import math
from dataclasses import dataclass
from typing import Protocol, Sequence


ASSIGNMENT_VERSION = "kmc_shared_thermo/1"
FROZEN_THERMO_PROPERTY = "kmc_frozen_thermo"


@dataclass(frozen=True)
class ReferenceThermoResult:
    """One atomic equilibrium calculation and the assignment that produced it."""

    equilibrium_constants: tuple[float, ...]
    reaction_enthalpy_298_J_per_mol: float
    species_thermo_assignments: tuple[dict, ...]


class ReferenceThermoProvider(Protocol):
    """A provider returns Kc in the displayed reaction's concentration units."""

    @property
    def provenance(self) -> dict: ...

    def equilibrium_constants(
        self, reaction, temperatures: Sequence[float]
    ) -> list[float]: ...


class ThermoUnavailable(ValueError):
    """A named species or reference-state calculation is unavailable."""


class SharedThermoAssignment:
    """Assign one pinned, resonance-aware thermo view without touching graphs."""

    def __init__(self, thermo_database, database_commit: str | None):
        self.thermo_database = thermo_database
        self.database_commit = database_commit
        self._cache = {}

    @property
    def library_order(self) -> list[str]:
        return list(getattr(self.thermo_database, "library_order", ()) or ())

    @property
    def provenance(self) -> dict:
        provider = self.thermo_database
        return {
            "assignment_version": ASSIGNMENT_VERSION,
            "thermo_provider": (
                f"{type(provider).__module__}.{type(provider).__qualname__}"
                if provider is not None
                else None
            ),
            "rmg_database_sha": self.database_commit,
            "thermo_library_order": self.library_order,
        }

    @staticmethod
    def _name(molecule) -> str:
        isolated = (
            molecule.copy(deep=True) if hasattr(molecule, "copy") else molecule
        )
        return (
            isolated.to_smiles()
            if hasattr(isolated, "to_smiles")
            else isolated.to_adjacency_list(remove_h=False)
        )

    @staticmethod
    def _thermo_source(thermo) -> dict:
        comment = str(getattr(thermo, "comment", "") or "")
        marker = "Thermo library:"
        if marker in comment:
            library = comment.split(marker, 1)[1].split("\n", 1)[0].strip()
            return {"kind": "library", "library": library}
        if "R-009 archived" in comment:
            return {"kind": "archived", "archive": "R-009"}
        if "group additivity" in comment.lower():
            return {"kind": "group-additivity"}
        return {"kind": "unknown", "comment": comment or None}

    @staticmethod
    def _value_si(thermo, name):
        value = getattr(thermo, name, None)
        return float(value.value_si) if value is not None else None

    def _isolated_reference(self, species):
        from rmgpy.molecule.molecule import Molecule
        from rmgpy.species import Species

        molecule = species.molecule[0]
        if not molecule.__class__.__module__.startswith("rmgpy."):
            return copy.deepcopy(species)
        reference = Species(
            label=getattr(species, "label", ""),
            molecule=[
                Molecule().from_adjacency_list(
                    molecule.to_adjacency_list(remove_h=False)
                )
            ],
        )
        reference.generate_resonance_structures()
        return reference

    def _identity(self, reference) -> tuple[dict, str]:
        resonance_smiles = sorted(
            {
                (
                    molecule.to_smiles()
                    if hasattr(molecule, "to_smiles")
                    else molecule.to_adjacency_list(remove_h=False)
                )
                for molecule in reference.molecule
            }
        )
        multiplicity = int(getattr(reference.molecule[0], "multiplicity", 1))
        chemical_identity = {
            "resonance_smiles": resonance_smiles,
            "multiplicity": multiplicity,
        }
        key_payload = {
            **self.provenance,
            "chemical_identity": chemical_identity,
        }
        key = hashlib.sha256(
            json.dumps(key_payload, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        return chemical_identity, key

    def assign_reaction(self, reaction) -> list[dict]:
        """Assign cached thermo to every occurrence and return its provenance."""
        assignments = []
        participants = (
            [
                ("reactant", index, species)
                for index, species in enumerate(reaction.reactants)
            ]
            + [
                ("product", index, species)
                for index, species in enumerate(reaction.products)
            ]
        )
        for role, index, species in participants:
            molecule = species.molecule[0]
            name = self._name(molecule)
            if self.thermo_database is None:
                raise ThermoUnavailable(
                    f"thermo unavailable for species {name}: no thermo database"
                )
            source_graph = molecule.to_adjacency_list(remove_h=False)
            try:
                reference = self._isolated_reference(species)
                chemical_identity, key = self._identity(reference)
                frozen = bool(
                    getattr(species, "props", {}).get(FROZEN_THERMO_PROPERTY, False)
                )
                if frozen:
                    thermo = getattr(species, "thermo", None)
                else:
                    if key not in self._cache:
                        self._cache[key] = self.thermo_database.get_thermo_data(
                            reference
                        )
                    thermo = self._cache[key]
                    species.thermo = thermo
                if thermo is None:
                    raise ValueError("RMG returned no thermochemistry")
            except Exception as error:
                raise ThermoUnavailable(
                    f"thermo unavailable for species {name}: {error}"
                ) from error
            assignments.append(
                {
                    "role": role,
                    "index": index,
                    "label": getattr(species, "label", ""),
                    "source_graph_sha256": hashlib.sha256(
                        source_graph.encode()
                    ).hexdigest(),
                    "chemical_identity_sha256": key,
                    "multiplicity": chemical_identity["multiplicity"],
                    "resonance_structure_count": len(reference.molecule),
                    "thermo_source": self._thermo_source(thermo),
                    "H298_J_per_mol": self._value_si(thermo, "H298"),
                    "E0_J_per_mol": self._value_si(thermo, "E0"),
                    "frozen_archive": frozen,
                }
            )
        return assignments


class GasPhaseRMGReferenceThermo:
    """Gas-phase RMG thermochemistry, without a melt/solvation correction."""

    def __init__(
        self,
        thermo_database,
        database_commit: str | None,
        assignment: SharedThermoAssignment | None = None,
    ):
        self.thermo_database = thermo_database
        self.database_commit = database_commit
        self.assignment = assignment or SharedThermoAssignment(
            thermo_database, database_commit
        )

    @property
    def provenance(self) -> dict:
        return {
            "reference_thermo": "RMG gas-phase Kc",
            "rmg_database_sha": self.database_commit,
            "condensed_phase_constraint": "UNKNOWN",
            **self.assignment.provenance,
        }

    def evaluate(
        self, reaction, temperatures: Sequence[float]
    ) -> ReferenceThermoResult:
        assignments = self.assignment.assign_reaction(reaction)
        reaction_enthalpy = float(reaction.get_enthalpy_of_reaction(298))
        try:
            values = [
                float(reaction.get_equilibrium_constant(temperature, type="Kc"))
                for temperature in temperatures
            ]
        except Exception as error:
            names = ", ".join(
                species.molecule[0].to_smiles()
                for species in list(reaction.reactants) + list(reaction.products)
            )
            raise ThermoUnavailable(
                f"Kc unavailable for species [{names}]: {error}"
            ) from error
        if any(not math.isfinite(value) or value <= 0 for value in values):
            names = ", ".join(
                species.molecule[0].to_smiles()
                for species in list(reaction.reactants) + list(reaction.products)
            )
            raise ThermoUnavailable(
                f"Kc unavailable for species [{names}]: non-positive or non-finite values at {list(temperatures)} K"
            )
        return ReferenceThermoResult(
            equilibrium_constants=tuple(values),
            reaction_enthalpy_298_J_per_mol=reaction_enthalpy,
            species_thermo_assignments=tuple(assignments),
        )

    def equilibrium_constants(
        self, reaction, temperatures: Sequence[float]
    ) -> list[float]:
        return list(self.evaluate(reaction, temperatures).equilibrium_constants)
