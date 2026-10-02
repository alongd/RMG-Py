"""Reference-state equilibrium constants for compiled reversible event pairs."""

from __future__ import annotations

import math
from typing import Protocol, Sequence


class ReferenceThermoProvider(Protocol):
    """A provider returns Kc in the displayed reaction's concentration units."""

    @property
    def provenance(self) -> dict: ...

    def equilibrium_constants(
        self, reaction, temperatures: Sequence[float]
    ) -> list[float]: ...


class ThermoUnavailable(ValueError):
    """A named species or reference-state calculation is unavailable."""


class GasPhaseRMGReferenceThermo:
    """Gas-phase RMG thermochemistry, without a melt/solvation correction."""

    def __init__(self, thermo_database, database_commit: str | None):
        self.thermo_database = thermo_database
        self.database_commit = database_commit
        self._cache = {}

    @property
    def provenance(self) -> dict:
        return {
            "reference_thermo": "RMG gas-phase Kc",
            "rmg_database_sha": self.database_commit,
            "condensed_phase_constraint": "UNKNOWN",
        }

    def equilibrium_constants(
        self, reaction, temperatures: Sequence[float]
    ) -> list[float]:
        from rmgpy.molecule.molecule import Molecule
        from rmgpy.species import Species

        for species in list(reaction.reactants) + list(reaction.products):
            molecule = species.molecule[0]
            name = (
                molecule.to_smiles()
                if hasattr(molecule, "to_smiles")
                else molecule.to_adjacency_list(remove_h=False)
            )
            if self.thermo_database is None:
                raise ThermoUnavailable(
                    f"thermo unavailable for species {name}: no thermo database"
                )
            key = name
            try:
                if key not in self._cache:
                    reference = Species(
                        molecule=[
                            Molecule().from_adjacency_list(
                                molecule.to_adjacency_list(remove_h=False)
                            )
                        ]
                    )
                    reference.generate_resonance_structures()
                    self._cache[key] = self.thermo_database.get_thermo_data(reference)
                species.thermo = self._cache[key]
                if species.thermo is None:
                    raise ValueError("RMG returned no thermochemistry")
            except Exception as error:
                raise ThermoUnavailable(
                    f"thermo unavailable for species {name}: {error}"
                ) from error
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
        return values
