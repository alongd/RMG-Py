"""Compiler-local zero-K enthalpy providers for the barrier-height path."""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass

import numpy as np


PROVIDER_NAME = "fixed-b-wilhoit"
PROVIDER_VERSION = "1"


class BarrierE0Error(ValueError):
    """A fixed-B provider input or fit cannot produce trustworthy E0."""


@dataclass(frozen=True)
class FixedBBarrierE0Provider:
    """Populate missing E0 on isolated thermo copies with one fixed B."""

    B: float

    def __post_init__(self):
        if isinstance(self.B, bool) or not isinstance(self.B, (int, float)):
            raise TypeError("B must be a finite positive number")
        value = float(self.B)
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError("B must be a finite positive number")
        object.__setattr__(self, "B", value)

    @property
    def provenance(self) -> dict:
        return {
            "enabled": True,
            "name": PROVIDER_NAME,
            "version": PROVIDER_VERSION,
            "B_K": self.B,
            "fit_temperature_grid": "participant ThermoData.Tdata",
            "fit_weights": "uniform least squares",
        }

    @staticmethod
    def _fit_inputs(thermo, B):
        try:
            temperatures = [float(value) for value in thermo.Tdata.value_si]
            capacities = [float(value) for value in thermo.Cpdata.value_si]
            cp0 = float(thermo.Cp0.value_si)
            cp_inf = float(thermo.CpInf.value_si)
            h298 = float(thermo.H298.value_si)
        except (AttributeError, TypeError) as error:
            raise BarrierE0Error(
                "fixed-B E0 provider requires Tdata, Cpdata, H298, "
                "Cp0, and CpInf"
            ) from error
        if len(temperatures) != len(capacities) or len(temperatures) < 4:
            raise BarrierE0Error(
                "fixed-B E0 provider requires at least four paired "
                "Tdata/Cpdata values"
            )
        values = temperatures + capacities + [cp0, cp_inf, h298]
        if not all(math.isfinite(value) for value in values):
            raise BarrierE0Error("fixed-B E0 provider inputs must be finite")
        if any(value <= 0.0 for value in temperatures):
            raise BarrierE0Error(
                "fixed-B E0 provider temperatures must be positive"
            )
        if any(
            second <= first
            for first, second in zip(temperatures, temperatures[1:])
        ):
            raise BarrierE0Error(
                "fixed-B E0 provider temperatures must be strictly increasing"
            )
        if any(value < 0.0 for value in capacities) or min(cp0, cp_inf) < 0.0:
            raise BarrierE0Error(
                "fixed-B E0 provider heat capacities must be non-negative"
            )
        if not math.isclose(cp0, cp_inf, rel_tol=1e-12, abs_tol=1e-12):
            y = np.asarray(temperatures) / (
                np.asarray(temperatures) + B
            )
            design = np.column_stack(
                [(y**3 - y**2) * y**power for power in range(4)]
            )
            if np.linalg.matrix_rank(design) < 4:
                raise BarrierE0Error(
                    "fixed-B E0 provider fit design is rank deficient"
                )
        return temperatures

    def prepare_reaction(self, reaction):
        """Return a barrier-only reaction and per-participant origins."""
        prepared = copy.deepcopy(reaction)
        prepared.reactants = [
            copy.deepcopy(species) for species in reaction.reactants
        ]
        prepared.products = [
            copy.deepcopy(species) for species in reaction.products
        ]
        assignments = []
        participants = [
            ("reactant", index, species)
            for index, species in enumerate(prepared.reactants)
        ] + [
            ("product", index, species)
            for index, species in enumerate(prepared.products)
        ]
        for role, index, species in participants:
            thermo = getattr(species, "thermo", None)
            if thermo is None:
                raise BarrierE0Error(
                    f"fixed-B E0 provider requires thermo for {role} {index}"
                )
            isolated = copy.deepcopy(thermo)
            supplied = getattr(isolated, "E0", None)
            if supplied is None:
                temperatures = self._fit_inputs(isolated, self.B)
                try:
                    fitted = isolated.to_wilhoit(B=self.B)
                except Exception as error:
                    raise BarrierE0Error(
                        f"fixed-B E0 fit failed for {role} {index}: {error}"
                    ) from error
                fitted_e0 = getattr(fitted, "E0", None)
                if fitted_e0 is None or not math.isfinite(
                    float(fitted_e0.value_si)
                ):
                    raise BarrierE0Error(
                        "fixed-B E0 fit returned no finite E0 for "
                        f"{role} {index}"
                    )
                isolated.E0 = (float(fitted_e0.value_si), "J/mol")
                origin = "provider"
                fit_temperatures = temperatures
                weights = [1.0] * len(temperatures)
            else:
                if not math.isfinite(float(supplied.value_si)):
                    raise BarrierE0Error(
                        f"supplied E0 is not finite for {role} {index}"
                    )
                origin = "supplied"
                fit_temperatures = None
                weights = None
            species.thermo = isolated
            assignments.append(
                {
                    "role": role,
                    "index": index,
                    "origin": origin,
                    "E0_J_per_mol": float(isolated.E0.value_si),
                    "temperature_grid_K": fit_temperatures,
                    "weights": weights,
                }
            )
        return prepared, assignments
