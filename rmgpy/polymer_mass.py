"""Shared physical mass and moment accounting for linear homopolymers.

This module is intentionally independent of both the method-of-moments solver
and the kMC runtime.  It is the common accounting surface consumed by their
thin adapters.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import Dict, Iterable, Mapping, Optional, Sequence, Tuple

from rmgpy.exceptions import ElementError
from rmgpy.molecule.element import get_element


_LEGACY_NOMINAL_MASS = {"C": 12, "H": 1, "N": 14, "O": 16, "S": 32}


def _isotope_number(value) -> int:
    if value in (None, -1, 0, ""):
        return -1
    if isinstance(value, bool):
        raise ValueError("an isotope number cannot be boolean")
    try:
        isotope = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("an isotope number must be an integer") from exc
    if isinstance(value, str):
        if value.strip() != str(isotope):
            raise ValueError("an isotope number must be an integer")
    elif isotope != value:
        raise ValueError("an isotope number must be an integer")
    if isotope <= 0:
        raise ValueError("an isotope number must be positive")
    return isotope


@dataclass(frozen=True)
class Composition:
    """An exact elemental/isotopic inventory.

    Counts may be signed so the same type can represent a defect relative to
    the repeat-unit inventory.  Physical graphs use only nonnegative counts.
    """

    entries: Tuple[Tuple[str, int, int], ...]

    @classmethod
    def from_mapping(cls, values: Mapping) -> "Composition":
        counts: Dict[Tuple[str, int], int] = {}
        for raw_key, raw_count in values.items():
            if isinstance(raw_key, tuple):
                if len(raw_key) != 2:
                    raise ValueError("composition tuple keys must be (element, isotope)")
                symbol, isotope = raw_key
            else:
                symbol, isotope = raw_key, -1
            if not isinstance(symbol, str) or not symbol:
                raise ValueError("composition element symbols must be non-empty strings")
            if isinstance(raw_count, bool):
                raise ValueError("composition counts cannot be boolean")
            count = int(raw_count)
            if count != raw_count:
                raise ValueError("composition counts must be integers")
            isotope = _isotope_number(isotope)
            if count:
                key = (symbol, isotope)
                counts[key] = counts.get(key, 0) + count
        return cls(tuple(
            (symbol, isotope, count)
            for (symbol, isotope), count in sorted(counts.items())
            if count
        ))

    @classmethod
    def from_atoms(cls, atoms: Sequence[Mapping]) -> "Composition":
        counts: Dict[Tuple[str, int], int] = {}
        for atom in atoms:
            symbol = atom.get("element")
            isotope = _isotope_number(atom.get("isotope", -1))
            if not isinstance(symbol, str) or not symbol:
                raise ValueError("every atom needs a non-empty element symbol")
            counts[(symbol, isotope)] = counts.get((symbol, isotope), 0) + 1
            hydrogens = atom.get("hydrogens", 0)
            if isinstance(hydrogens, bool) or int(hydrogens) != hydrogens or hydrogens < 0:
                raise ValueError("atom hydrogens must be a nonnegative integer")
            if hydrogens:
                counts[("H", -1)] = counts.get(("H", -1), 0) + int(hydrogens)
        return cls.from_mapping(counts)

    def as_mapping(self) -> Dict[Tuple[str, int], int]:
        return {(symbol, isotope): count for symbol, isotope, count in self.entries}

    def plus(self, other: "Composition") -> "Composition":
        values = self.as_mapping()
        for key, count in other.as_mapping().items():
            values[key] = values.get(key, 0) + count
        return Composition.from_mapping(values)

    def minus(self, other: "Composition") -> "Composition":
        return self.plus(other.scaled(-1))

    def scaled(self, multiplier: int) -> "Composition":
        if isinstance(multiplier, bool) or int(multiplier) != multiplier:
            raise ValueError("composition multipliers must be integers")
        return Composition.from_mapping({
            (symbol, isotope): count * int(multiplier)
            for symbol, isotope, count in self.entries
        })

    def molar_mass_kg_mol(self) -> float:
        total = 0.0
        for symbol, isotope, count in self.entries:
            try:
                element = get_element(symbol, isotope)
            except ElementError as exc:
                raise ValueError(
                    "unknown element/isotope in physical mass inventory: "
                    f"{symbol!r}, isotope={isotope}"
                ) from exc
            if element.number <= 0 or element.mass <= 0.0:
                raise ValueError(
                    f"element {symbol!r} is not a physical positive-mass atom"
                )
            total += count * element.mass
        return total

    def molar_mass_g_mol(self) -> float:
        return 1000.0 * self.molar_mass_kg_mol()

    def legacy_nominal_mass_g_mol(self) -> float:
        """Reproduce the pre-Correction-5 integer/unknown-zero diagnostic."""
        return float(sum(
            _LEGACY_NOMINAL_MASS.get(symbol, 0) * count
            for symbol, _isotope, count in self.entries
        ))


@dataclass(frozen=True)
class MomentObservables:
    mass_g: float
    mn_g_mol: float
    mw_g_mol: float


def observables_from_moments(
    moments: Sequence[float],
    repeat_mass_g_mol: float,
    end_mass_g_mol: float,
) -> MomentObservables:
    """Return end-corrected mass, Mn, and Mw from DP moments."""
    mu0, mu1, mu2 = _validate_moments(moments)
    a = _positive(repeat_mass_g_mol, "repeat molar mass")
    b = _nonnegative(end_mass_g_mol, "end-group molar mass")
    if mu0 == 0.0:
        return MomentObservables(0.0, 0.0, 0.0)
    mass = a * mu1 + b * mu0
    mn = mass / mu0
    mw = (a * a * mu2 + 2.0 * a * b * mu1 + b * b * mu0) / mass
    return MomentObservables(mass, mn, mw)


def moments_from_observables(
    mass_g: float,
    mn_g_mol: float,
    mw_g_mol: float,
    repeat_mass_g_mol: float,
    end_mass_g_mol: float,
) -> Tuple[float, float, float]:
    """Reconstruct end-corrected DP moments from mass, Mn, and Mw."""
    mass = _nonnegative(mass_g, "mass")
    if mass == 0.0:
        return 0.0, 0.0, 0.0
    mn = _positive(mn_g_mol, "Mn")
    mw = _positive(mw_g_mol, "Mw")
    a = _positive(repeat_mass_g_mol, "repeat molar mass")
    b = _nonnegative(end_mass_g_mol, "end-group molar mass")
    if mn < b:
        raise ValueError("Mn is smaller than the two-end-group molar mass")
    mu0 = mass / mn
    mu1 = mu0 * (mn - b) / a
    mu2 = (mw * (a * mu1 + b * mu0) - 2.0 * a * b * mu1 - b * b * mu0) / (a * a)
    return _validate_moments((mu0, mu1, mu2))


def explicit_moment_bundle(dp_moles: Mapping[int, float]) -> Tuple[float, float, float]:
    mu0 = mu1 = mu2 = 0.0
    for raw_dp, raw_moles in dp_moles.items():
        if isinstance(raw_dp, bool) or int(raw_dp) != raw_dp or raw_dp <= 0:
            raise ValueError("explicit-chain DP values must be positive integers")
        moles = _nonnegative(raw_moles, "explicit-chain amount")
        dp = int(raw_dp)
        mu0 += moles
        mu1 += moles * dp
        mu2 += moles * dp * dp
    return mu0, mu1, mu2


def subtract_explicit_moments(
    moments: Sequence[float], dp_moles: Mapping[int, float]
) -> Tuple[float, float, float]:
    values = _validate_moments(moments)
    debit = explicit_moment_bundle(dp_moles)
    result = tuple(values[index] - debit[index] for index in range(3))
    if any(value < -1.0e-15 for value in result):
        raise ValueError("explicit-chain moment debit exceeds the pool moments")
    return _validate_moments(result)


@dataclass(frozen=True)
class InventoryTerm:
    provenance_id: str
    composition: Composition
    amount_mol: float
    chain_coupled: bool = False

    def __post_init__(self):
        if not isinstance(self.provenance_id, str) or not self.provenance_id:
            raise ValueError("inventory provenance IDs must be non-empty strings")
        _nonnegative(self.amount_mol, "inventory amount")

    @property
    def mass_g(self) -> float:
        return self.composition.molar_mass_g_mol() * self.amount_mol


@dataclass(frozen=True)
class ExplicitInventory:
    provenance_id: str
    composition: Composition
    amount_mol: float

    def __post_init__(self):
        if not isinstance(self.provenance_id, str) or not self.provenance_id:
            raise ValueError("explicit-inventory provenance IDs must be non-empty strings")
        _nonnegative(self.amount_mol, "explicit-inventory amount")

    @property
    def mass_g(self) -> float:
        return self.composition.molar_mass_g_mol() * self.amount_mol


@dataclass(frozen=True)
class EvolvedPoolMass:
    """Exclusive evolved-pool mass law from design section 2.3."""

    label: str
    repeat: Composition
    moments: Tuple[float, float, float]
    end_inventory: Tuple[InventoryTerm, ...] = ()
    defects: Tuple[InventoryTerm, ...] = ()
    losses: Tuple[InventoryTerm, ...] = ()

    def __post_init__(self):
        _validate_moments(self.moments)
        ids = [term.provenance_id for term in self.end_inventory + self.defects + self.losses]
        if len(ids) != len(set(ids)):
            raise ValueError("end, defect, and loss provenance IDs must be disjoint")

    @property
    def mass_g(self) -> float:
        _mu0, mu1, _mu2 = self.moments
        return (
            self.repeat.molar_mass_g_mol() * mu1
            + sum(term.mass_g for term in self.end_inventory)
            + sum(term.mass_g for term in self.defects)
            - sum(term.mass_g for term in self.losses)
        )

    def depropagate(
        self,
        units_mol: float,
        destination_id: str,
        mu0_debit: float = 0.0,
        mu2_debit: Optional[float] = None,
    ) -> Tuple["EvolvedPoolMass", ExplicitInventory]:
        """Debit repeat-unit moments and credit monomer, never the loss ledger."""
        units = _nonnegative(units_mol, "depropagated units")
        mu0_delta = _nonnegative(mu0_debit, "depropagated chain debit")
        if mu2_debit is None:
            raise ValueError(
                "mu2_debit must be explicit or derived from the actual depropagation operator"
            )
        mu2_delta = _nonnegative(mu2_debit, "depropagated second-moment debit")
        mu0, mu1, mu2 = self.moments
        new_moments = (mu0 - mu0_delta, mu1 - units, mu2 - mu2_delta)
        if any(value < -1.0e-15 for value in new_moments):
            raise ValueError("depropagation debit exceeds the pool moments")
        updated = replace(self, moments=tuple(max(0.0, value) for value in new_moments))
        return updated, ExplicitInventory(destination_id, self.repeat, units)

    def transfer_complete_chain_bundle(
        self, dp: int, chains_mol: float, destination_id: str
    ) -> Tuple["EvolvedPoolMass", ExplicitInventory]:
        """Move moments plus chain-coupled end/defect inventory exactly once."""
        if isinstance(dp, bool) or int(dp) != dp or dp <= 0:
            raise ValueError("transferred-chain DP must be a positive integer")
        amount = _nonnegative(chains_mol, "transferred chains")
        mu0, mu1, mu2 = self.moments
        new_moments = (mu0 - amount, mu1 - dp * amount, mu2 - dp * dp * amount)
        if any(value < -1.0e-15 for value in new_moments):
            raise ValueError("complete chain bundle exceeds the pool moments")

        explicit_composition = self.repeat.scaled(dp)

        def debit_coupled(terms: Tuple[InventoryTerm, ...]):
            nonlocal explicit_composition
            updated = []
            for term in terms:
                if term.chain_coupled:
                    if term.amount_mol + 1.0e-15 < amount:
                        raise ValueError("chain-coupled inventory is smaller than the transfer")
                    explicit_composition = explicit_composition.plus(term.composition)
                    updated.append(replace(term, amount_mol=max(0.0, term.amount_mol - amount)))
                else:
                    updated.append(term)
            return tuple(updated)

        ends = debit_coupled(self.end_inventory)
        defects = debit_coupled(self.defects)
        if any(term.chain_coupled for term in self.losses):
            raise ValueError("stationary loss inventory cannot be chain-coupled")
        updated = replace(
            self,
            moments=tuple(max(0.0, value) for value in new_moments),
            end_inventory=ends,
            defects=defects,
        )
        return updated, ExplicitInventory(destination_id, explicit_composition, amount)

    def deplete_stationary_inventory(
        self,
        provenance_id: str,
        composition: Composition,
        amount_mol: float,
        destination_id: str,
    ) -> Tuple["EvolvedPoolMass", ExplicitInventory]:
        """Record a stationary-inventory loss and its destination credit once."""
        amount = _nonnegative(amount_mol, "stationary depletion")
        if provenance_id in {
            term.provenance_id
            for term in self.end_inventory + self.defects + self.losses
        }:
            raise ValueError("stationary-depletion provenance must be unique")
        requested = composition.as_mapping()
        if any(count < 0 for count in requested.values()):
            raise ValueError("stationary depletion must name a nonnegative physical inventory")
        available = self._elemental_inventory_mol()
        for key, count in requested.items():
            needed = count * amount
            present = available.get(key, 0.0)
            tolerance = 1.0e-12 * max(abs(needed), abs(present), 1.0e-30)
            if present + tolerance < needed:
                symbol, isotope = key
                raise ValueError(
                    "stationary depletion exceeds available elemental inventory for "
                    f"{symbol} isotope={isotope}: requested {needed:g} mol-atoms, "
                    f"available {present:g} mol-atoms"
                )
        term = InventoryTerm(provenance_id, composition, amount, False)
        updated = replace(self, losses=self.losses + (term,))
        return updated, ExplicitInventory(destination_id, composition, amount)

    def _elemental_inventory_mol(self) -> Dict[Tuple[str, int], float]:
        """Return the live atom inventory represented by moments/E/D/L."""
        inventory: Dict[Tuple[str, int], float] = {}

        def add(composition: Composition, amount: float):
            for key, count in composition.as_mapping().items():
                inventory[key] = inventory.get(key, 0.0) + count * amount

        add(self.repeat, self.moments[1])
        for term in self.end_inventory + self.defects:
            add(term.composition, term.amount_mol)
        for term in self.losses:
            add(term.composition, -term.amount_mol)
        return inventory


@dataclass(frozen=True)
class MassVolumeLedger:
    pools: Tuple[EvolvedPoolMass, ...]
    explicit: Tuple[ExplicitInventory, ...]
    density_kg_m3: float

    def __post_init__(self):
        _positive(self.density_kg_m3, "mixture density")

    @property
    def mass_kg(self) -> float:
        return (
            sum(pool.mass_g for pool in self.pools)
            + sum(item.mass_g for item in self.explicit)
        ) / 1000.0

    @property
    def volume_m3(self) -> float:
        return self.mass_kg / self.density_kg_m3


def _validate_moments(values: Sequence[float]) -> Tuple[float, float, float]:
    if len(values) < 3:
        raise ValueError("moments must contain mu0, mu1, and mu2")
    mu0, mu1, mu2 = (float(values[index]) for index in range(3))
    if not all(math.isfinite(value) for value in (mu0, mu1, mu2)):
        raise ValueError("moments must be finite")
    if min(mu0, mu1, mu2) < -1.0e-15:
        raise ValueError("moments cannot be negative")
    mu0, mu1, mu2 = max(0.0, mu0), max(0.0, mu1), max(0.0, mu2)
    if mu0 == 0.0:
        if mu1 or mu2:
            raise ValueError("nonzero higher moments require nonzero mu0")
        return mu0, mu1, mu2
    if mu1 + 1.0e-15 < mu0:
        raise ValueError("a positive-DP distribution requires mu1 >= mu0")
    if mu2 + 1.0e-15 < mu1:
        raise ValueError("a positive-DP distribution requires mu2 >= mu1")
    lhs = mu0 * mu2
    rhs = mu1 * mu1
    if lhs + 1.0e-12 * max(abs(lhs), abs(rhs)) < rhs:
        raise ValueError("moments violate mu0*mu2 >= mu1^2")
    return mu0, mu1, mu2


def _nonnegative(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return value


def _positive(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value
