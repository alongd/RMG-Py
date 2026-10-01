"""Melt encounter-limited radical termination and reversible geminate cages.

Transport is always injected through one of the pre-registered arm objects.
This module consumes compiled rewrite records; it does not generate chemistry,
choose an arm or inventory, or run a stochastic simulation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from types import MappingProxyType
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np


R_RMG = 8.314472
K_B = 1.380649e-23
N_A = 6.02214076e23
R_SI = K_B * N_A
P_STANDARD = 1.0e5
SIGMA_CONTACT = 6.766081442101e-10
PS_M0 = 104.15
PS_C_R2 = 0.434e-20
RATE_TOLERANCE = 1.0e-11
RATE_GRID = (600.0, 650.0, 700.0, 750.0, 800.0)
PAIR_CLASSES = frozenset({"end/end", "end/mid", "mid/mid"})
FALLBACK_STRESS_GRID = (0.0, 0.25, 0.5, 0.75, 1.0)
ESCAPED_RADICAL_PAIR = (
    "CCC[CH]c1ccccc1",
    "CC(c1ccccc1)[CH2]",
)
R0_DISPROPORTIONATION_PRODUCTS = (
    "C=C(C)C1=CC=CC=C1 + CCCCC1=CC=CC=C1",
    "CC(C)C1=CC=CC=C1 + CCC=CC1=CC=CC=C1",
)


class UnsupportedTopologyError(ValueError):
    """The MET transport law is undefined for this component topology."""


class CoverageError(ValueError):
    """The compiled radical-decreasing inventory is not fully classified."""


class CageBalanceError(ValueError):
    """A cage cannot absorb into the required permanent fate set."""

    def __init__(self, required_terminal_fate_mass: float):
        self.required_terminal_fate_mass = float(required_terminal_fate_mass)
        super().__init__(
            "required terminal fate mass is "
            f"{self.required_terminal_fate_mass:.13g}, not 1"
        )


class MissingMETChannelError(ValueError):
    """A requested inventory is absent from the compiled event set."""

    def __init__(self, channel_id: str, kernel: str):
        self.channel_id = channel_id
        self.kernel = kernel
        super().__init__(
            f"{kernel} R1 compilation is missing required MET channel {channel_id}"
        )


@dataclass(frozen=True)
class RateTable:
    """A positive rate table with the compiler's linear-in-ln(k) contract."""

    temperatures: tuple[float, ...]
    rates: tuple[float, ...]

    def __post_init__(self) -> None:
        temperatures = tuple(float(value) for value in self.temperatures)
        rates = tuple(float(value) for value in self.rates)
        if len(temperatures) != len(rates) or not temperatures:
            raise ValueError("rate table temperatures and rates must align")
        if temperatures != tuple(sorted(set(temperatures))):
            raise ValueError("rate table temperatures must be strictly increasing")
        if any(not math.isfinite(value) or value <= 0.0 for value in rates):
            raise ValueError("rate table rates must be finite and positive")
        object.__setattr__(self, "temperatures", temperatures)
        object.__setattr__(self, "rates", rates)

    @classmethod
    def from_mapping(cls, table: Mapping[str, Any]) -> "RateTable":
        if table.get("interpolation") != "linear-ln-k":
            raise ValueError("MET requires linear-ln-k interpolation")
        if table.get("extrapolation") != "refuse":
            raise ValueError("MET refuses rate-table extrapolation")
        return cls(tuple(table["T"]), tuple(table["k"]))

    def to_mapping(self) -> dict[str, Any]:
        return {
            "T": list(self.temperatures),
            "k": list(self.rates),
            "interpolation": "linear-ln-k",
            "extrapolation": "refuse",
        }

    def __call__(self, temperature: float) -> float:
        temperature = float(temperature)
        if not math.isfinite(temperature):
            raise ValueError("temperature must be finite")
        if temperature < self.temperatures[0] or temperature > self.temperatures[-1]:
            raise ValueError("temperature is outside the compiled rate table")
        index = int(np.searchsorted(self.temperatures, temperature))
        if index < len(self.temperatures) and self.temperatures[index] == temperature:
            return self.rates[index]
        lower = index - 1
        fraction = (temperature - self.temperatures[lower]) / (
            self.temperatures[index] - self.temperatures[lower]
        )
        log_rate = math.log(self.rates[lower]) + fraction * (
            math.log(self.rates[index]) - math.log(self.rates[lower])
        )
        return math.exp(log_rate)


@dataclass(frozen=True)
class ArrheniusRate:
    """One archived RMG Arrhenius object evaluated with R_RMG."""

    preexponential: float
    temperature_exponent: float
    activation_energy: float
    reference_temperature: float = 1.0

    def __call__(self, temperature: float) -> float:
        temperature = _positive_finite(temperature, "temperature")
        return (
            self.preexponential
            * (temperature / self.reference_temperature) ** self.temperature_exponent
            * math.exp(-self.activation_energy / (R_RMG * temperature))
        )


ARCHIVED_CHANNEL_RATES: Mapping[str, ArrheniusRate] = MappingProxyType(
    {
        "D1": ArrheniusRate(1.051628e8, -0.55, 11935.1976504764516),
        "D2": ArrheniusRate(5.25814e7, -0.55, 8686.37878983243354),
        "D3": ArrheniusRate(5.25814e7, -0.55, 0.0),
        "D4": ArrheniusRate(2.9e6, 0.0, 0.0),
        "D5": ArrheniusRate(168600.0, 0.0, 113223.849350046265),
        "null": ArrheniusRate(2.4834e8, -0.557189, 0.0),
        "J_para": ArrheniusRate(1.76793e10, -1.00291, 0.0),
        "J_ortho": ArrheniusRate(3.53586e10, -1.00291, 0.0),
    }
)
ARCHIVED_JUNCTION_RATES = MappingProxyType(
    {
        channel_id: ARCHIVED_CHANNEL_RATES[channel_id]
        for channel_id in ("J_para", "J_ortho")
    }
)

# Two-range NASA records archived by the binding transport pack.  These are
# the labelled gas-phase RMG reference model, not condensed-phase thermochemistry.
_NASA: Mapping[str, tuple[tuple[float, float, tuple[float, ...]], ...]] = (
    MappingProxyType(
        {
            "prim": (
                (
                    300.0,
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
                ),
                (
                    1177.6065961907934,
                    3000.0,
                    (
                        12.837928688294944,
                        0.04618929187592041,
                        -2.208327352961384e-05,
                        5.680597353073413e-09,
                        -5.868026777213613e-13,
                        17808.42272068456,
                        -40.067143782909014,
                    ),
                ),
            ),
            "sec": (
                (
                    300.0,
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
                ),
                (
                    1254.2536410728535,
                    3000.0,
                    (
                        2.8546259360984383,
                        0.08172943981695815,
                        -4.6397539366542894e-05,
                        1.1664266596876454e-08,
                        -1.1165136655519307e-12,
                        11827.021589717222,
                        10.199278347265654,
                    ),
                ),
            ),
            "para": (
                (
                    300.0,
                    1403.6826169506783,
                    (
                        -10.949627764777754,
                        0.20969652550530965,
                        -0.00015924632577030666,
                        6.232354330618766e-08,
                        -9.94128604215258e-12,
                        11636.982987189213,
                        85.42480168403608,
                    ),
                ),
                (
                    1403.6826169506783,
                    3000.0,
                    (
                        22.36559012247251,
                        0.11475767839669494,
                        -5.779061059436066e-05,
                        1.4136912533315018e-08,
                        -1.3589033989204612e-12,
                        2284.404766054439,
                        -86.59824668710979,
                    ),
                ),
            ),
            "ortho": (
                (
                    300.0,
                    1171.6898321871104,
                    (
                        -8.455190642047764,
                        0.19501882679191176,
                        -0.00013017705248095492,
                        3.885484757972701e-08,
                        -3.261608985996416e-12,
                        13316.25635917201,
                        75.03280013132829,
                    ),
                ),
                (
                    1171.6898321871104,
                    3000.0,
                    (
                        18.46952197543009,
                        0.12343415887301927,
                        -6.45643901319227e-05,
                        1.633326197511283e-08,
                        -1.6163446936598018e-12,
                        5611.067403002024,
                        -65.08547018371374,
                    ),
                ),
            ),
        }
    )
)


def _nasa_gibbs(temperature: float, species: str) -> float:
    temperature = _positive_finite(temperature, "temperature")
    try:
        record = next(
            item for item in _NASA[species] if item[0] <= temperature <= item[1]
        )
    except (KeyError, StopIteration) as error:
        raise ValueError(
            f"temperature is outside the archived NASA model for {species}"
        ) from error
    c0, c1, c2, c3, c4, c5, c6 = record[2]
    h_rt = (
        c0
        + c1 * temperature / 2.0
        + c2 * temperature**2 / 3.0
        + c3 * temperature**3 / 4.0
        + c4 * temperature**4 / 5.0
        + c5 / temperature
    )
    s_r = (
        c0 * math.log(temperature)
        + c1 * temperature
        + c2 * temperature**2 / 2.0
        + c3 * temperature**3 / 3.0
        + c4 * temperature**4 / 4.0
        + c6
    )
    return R_RMG * temperature * (h_rt - s_r)


def rmg_reference_equilibrium_constant(temperature: float, junction: str) -> float:
    """Continuous Kc for the archived gas-phase RMG junction reference."""
    product = {"J_para": "para", "J_ortho": "ortho"}.get(junction)
    if product is None:
        raise ValueError(f"no archived RMG reference model for {junction}")
    delta_g = (
        _nasa_gibbs(temperature, product)
        - _nasa_gibbs(temperature, "sec")
        - _nasa_gibbs(temperature, "prim")
    )
    return math.exp(-delta_g / (R_RMG * temperature)) / (
        P_STANDARD / (R_RMG * temperature)
    )


@dataclass(frozen=True)
class TransportArm:
    """One frozen D0 package and its one permitted chain-scaling law."""

    arm_id: str
    package: str
    scaling: str
    n_e: float | None

    def d0(self, temperature: float) -> float:
        temperature = _positive_finite(temperature, "temperature")
        if self.package == "H":
            return _d0_h(temperature)
        if self.package == "U":
            return _d0_u(temperature)
        if self.package == "B":
            return _d0_b(temperature)
        raise ValueError(f"unknown transport package: {self.package}")

    def chain_diffusivity(self, temperature: float, units: float) -> float:
        units = _positive_finite(units, "chain length")
        d0 = self.d0(temperature)
        if self.scaling == "rouse":
            return d0 / units
        if self.scaling == "cross" and self.n_e is not None:
            return d0 / units if units <= self.n_e else d0 * self.n_e / units**2
        raise ValueError(f"invalid transport scaling: {self.scaling}")

    def signature(self) -> tuple[Any, ...]:
        return (
            self.arm_id,
            self.package,
            self.scaling,
            self.n_e,
            tuple(self.d0(temperature) for temperature in RATE_GRID),
        )


def _positive_finite(value: float, label: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{label} must be finite and positive")
    return value


def _log10_shift(temperature: float, c1: float, c2: float, reference: float) -> float:
    return -c1 * (temperature - reference) / (c2 + temperature - reference)


def _wlf(
    temperature: float,
    anchor_d: float,
    anchor_t: float,
    c1: float,
    c2: float,
    reference: float,
) -> float:
    return anchor_d * 10.0 ** (
        _log10_shift(anchor_t, c1, c2, reference)
        - _log10_shift(temperature, c1, c2, reference)
    )


_T_FERRY = 473.0
_D_FERRY = K_B * _T_FERRY / (10.0**-6.95 * 1.0e-3)
_D0_U = 17.5 * 10.0**-5.900749063670412 * 1.0e-4
_D0_B = (3600.0 / PS_M0) * 7.0e-11


def _d0_h(temperature: float) -> float:
    return _wlf(temperature, _D_FERRY, _T_FERRY, 13.39, 45.28, 378.15)


def _d0_u(temperature: float) -> float:
    return _wlf(temperature, _D0_U, 581.15, 9.76, 60.45, 59.0 + 273.15 + 12.2)


def _d0_b(temperature: float) -> float:
    return _D0_B * math.exp(-62.2e3 / R_SI * (1.0 / temperature - 1.0 / 503.15))


_NE_C = 14800.0 / PS_M0
_NE_LO = 14800.0 * 0.195 / 0.210 / PS_M0
_NE_HI = 14800.0 * 0.195 / 0.180 / PS_M0
TRANSPORT_ARMS: Mapping[str, TransportArm] = MappingProxyType(
    {
        "A0_REF_H_CROSS_NEc": TransportArm("A0_REF_H_CROSS_NEc", "H", "cross", _NE_C),
        "A1_H_ROUSE": TransportArm("A1_H_ROUSE", "H", "rouse", None),
        "A2_H_CROSS_NElo": TransportArm("A2_H_CROSS_NElo", "H", "cross", _NE_LO),
        "A3_H_CROSS_NEhi": TransportArm("A3_H_CROSS_NEhi", "H", "cross", _NE_HI),
        "A4_LU_ROUSE": TransportArm("A4_LU_ROUSE", "U", "rouse", None),
        "A5_LU_CROSS_NElo": TransportArm("A5_LU_CROSS_NElo", "U", "cross", _NE_LO),
        "A6_LU_CROSS_NEhi": TransportArm("A6_LU_CROSS_NEhi", "U", "cross", _NE_HI),
        "A7_LB_ROUSE": TransportArm("A7_LB_ROUSE", "B", "rouse", None),
        "A8_LB_CROSS_NElo": TransportArm("A8_LB_CROSS_NElo", "B", "cross", _NE_LO),
        "A9_LB_CROSS_NEhi": TransportArm("A9_LB_CROSS_NEhi", "B", "cross", _NE_HI),
    }
)


def diffusion_rate(
    arm: TransportArm,
    temperature: float,
    units_i: float,
    units_j: float,
    pair_class: str,
    *,
    spin_factor: float,
) -> float:
    """Return the report's long-time k_D for one of its three pair classes.

    The report defines no class-specific multiplier; the class is explicit so a
    future reference test cannot silently conflate end/end, end/mid and mid/mid.
    """
    if pair_class not in PAIR_CLASSES:
        raise ValueError(f"unknown MET pair class: {pair_class!r}")
    spin_factor = _positive_finite(spin_factor, "spin factor")
    if spin_factor > 1.0:
        raise ValueError("spin factor must not exceed one")
    units_i = _positive_finite(units_i, "chain length i")
    units_j = _positive_finite(units_j, "chain length j")
    radius_gyration = math.sqrt(PS_C_R2 * PS_M0 * min(units_i, units_j) / 6.0)
    capture_radius = max(SIGMA_CONTACT, 2.0 * radius_gyration)
    diffusivity = arm.chain_diffusivity(temperature, units_i) + arm.chain_diffusivity(
        temperature, units_j
    )
    return 4.0 * math.pi * N_A * spin_factor * diffusivity * capture_radius


def collins_kimball(k_act: float, k_diffusion: float) -> float:
    """Combine finite positive activation and diffusion rates in series."""
    k_act = float(k_act)
    k_diffusion = float(k_diffusion)
    if math.isinf(k_act) and k_act > 0.0:
        return _positive_finite(k_diffusion, "diffusion rate")
    if math.isinf(k_diffusion) and k_diffusion > 0.0:
        return _positive_finite(k_act, "activation rate")
    k_act = _positive_finite(k_act, "activation rate")
    k_diffusion = _positive_finite(k_diffusion, "diffusion rate")
    return 1.0 / (1.0 / k_act + 1.0 / k_diffusion)


def _field(record: Any, name: str, default: Any = None) -> Any:
    if isinstance(record, Mapping):
        return record.get(name, default)
    return getattr(record, name, default)


def record_activation_rate(record: Any, temperature: float) -> float:
    """Evaluate a compiled RMG rate with its recorded SSA pair correction."""
    table = _field(record, "k_table")
    if not table:
        raise ValueError("compiled record has no RMG rate table")
    multiplier = _positive_finite(
        _field(record, "ssa_multiplier", 1.0), "SSA multiplier"
    )
    return RateTable.from_mapping(table)(temperature) * multiplier


class CoverageClass(str, Enum):
    BIMolecular_MET = "bimolecular MET"
    UNIMOLECULAR = "unimolecular event"
    HELD_OUT = "held out"
    OUT_OF_SCOPE = "out of scope"


@dataclass(frozen=True)
class FamilyClassification:
    category: CoverageClass
    reason: str


FAMILY_CLASSIFICATION: Mapping[str, FamilyClassification] = MappingProxyType(
    {
        "R_Recombination": FamilyClassification(
            CoverageClass.BIMolecular_MET,
            "C/H radical recombination is a candidate RMG contact-reactivity source",
        ),
        "Disproportionation": FamilyClassification(
            CoverageClass.BIMolecular_MET,
            "C/H radical disproportionation is a candidate RMG contact-reactivity source",
        ),
        "R_Addition_MultipleBond_Disprop": FamilyClassification(
            CoverageClass.HELD_OUT,
            "not established for C/H polymers; two lithium training reactions",
        ),
        "Intra_Disproportionation": FamilyClassification(
            CoverageClass.UNIMOLECULAR,
            "same-component biradical disproportionation is not a pair encounter",
        ),
        "Birad_recombination": FamilyClassification(
            CoverageClass.UNIMOLECULAR,
            "same-component biradical ring closure is not a pair encounter",
        ),
        "1,2-Birad_to_alkene": FamilyClassification(
            CoverageClass.UNIMOLECULAR,
            "adjacent biradical collapse is a separate unimolecular event",
        ),
        "Disproportionation-Y": FamilyClassification(
            CoverageClass.OUT_OF_SCOPE,
            "the root requires a halogen donor absent from scoped PS chemistry",
        ),
        "CO_Disproportionation": FamilyClassification(
            CoverageClass.OUT_OF_SCOPE,
            "formyl/CO chemistry is absent from the scoped radical inventory",
        ),
        "Peroxyl_Termination": FamilyClassification(
            CoverageClass.OUT_OF_SCOPE,
            "oxygen/peroxyl chemistry is outside inert-atmosphere pyrolysis scope",
        ),
        "Birad_R_Recombination": FamilyClassification(
            CoverageClass.OUT_OF_SCOPE,
            "reserved for O/S/N triplet biradicals and has no fitted rules",
        ),
        "1,4_Cyclic_birad_scission": FamilyClassification(
            CoverageClass.UNIMOLECULAR,
            "cyclic biradical scission is a separate unimolecular event",
        ),
        "1,4_Linear_birad_scission": FamilyClassification(
            CoverageClass.UNIMOLECULAR,
            "linear biradical scission is a separate unimolecular event",
        ),
        "Cation_Addition_MultipleBond_Disprop": FamilyClassification(
            CoverageClass.OUT_OF_SCOPE,
            "cation chemistry is absent from the neutral PS radical inventory",
        ),
        "Cation_R_Recombination": FamilyClassification(
            CoverageClass.OUT_OF_SCOPE,
            "cation chemistry is absent from the neutral PS radical inventory",
        ),
        "Surface_Adsorption_Double": FamilyClassification(
            CoverageClass.OUT_OF_SCOPE,
            "surface adsorption is outside the homogeneous polymer melt model",
        ),
        "Surface_Adsorption_Single": FamilyClassification(
            CoverageClass.OUT_OF_SCOPE,
            "surface adsorption is outside the homogeneous polymer melt model",
        ),
    }
)
RADICAL_DECREASING_FAMILY_INVENTORY = tuple(FAMILY_CLASSIFICATION)


@dataclass(frozen=True)
class CoverageReport:
    classified_families: Mapping[str, FamilyClassification]
    library_reactions: tuple[str, ...]
    admitted_event_ids: tuple[str, ...]


def enumerate_radical_decreasing_families(kinetics_database: Any) -> tuple[str, ...]:
    """Inspect every loaded family recipe for a net radical decrease."""
    result = []
    for name, family in _field(kinetics_database, "families", {}).items():
        recipe = _field(family, "forward_recipe") or _field(family, "recipe")
        actions = _field(recipe, "actions", ()) if recipe is not None else ()
        balance = 0
        for action in actions:
            kind = str(action[0]).upper() if action else ""
            amount = int(action[2]) if len(action) > 2 else 1
            if kind == "LOSE_RADICAL":
                balance -= amount
            elif kind == "GAIN_RADICAL":
                balance += amount
        if balance < 0:
            result.append(str(name))
    return tuple(sorted(result))


def _participant_radicals(participants: Iterable[Any]) -> int:
    total = 0
    for participant in participants:
        molecules = _field(participant, "molecule")
        molecule = molecules[0] if molecules else participant
        total += sum(
            int(_field(atom, "radical_electrons", 0))
            for atom in _field(molecule, "atoms", ())
        )
    return total


def enumerate_radical_decreasing_library_reactions(
    kinetics_database: Any,
) -> tuple[str, ...]:
    """Inspect every entry in every loaded reaction library."""
    result = []
    for library_name, library in _field(kinetics_database, "libraries", {}).items():
        for entry_key, entry in _field(library, "entries", {}).items():
            reaction = _field(entry, "item", entry)
            reactants = _field(reaction, "reactants", ())
            products = _field(reaction, "products", ())
            if _participant_radicals(products) < _participant_radicals(reactants):
                result.append(f"{library_name}:{entry_key}")
    return tuple(sorted(result))


def audit_met_coverage(
    records: Iterable[Any],
    loaded_families: Iterable[str] | Any,
    library_classification: Mapping[str, FamilyClassification] | None = None,
    *,
    classification: Mapping[str, FamilyClassification] = FAMILY_CLASSIFICATION,
) -> CoverageReport:
    """Cross-check the static inventory against loaded and compiled chemistry."""
    database_families = _field(loaded_families, "families")
    if isinstance(database_families, Mapping):
        expected_loaded = set(enumerate_radical_decreasing_families(loaded_families))
        loaded_library_inventory = set(
            enumerate_radical_decreasing_library_reactions(loaded_families)
        )
    else:
        loaded_universe = set(loaded_families)
        expected_loaded = set(RADICAL_DECREASING_FAMILY_INVENTORY) & loaded_universe
        loaded_library_inventory = set()
    missing = sorted(expected_loaded - set(classification))
    if missing:
        raise CoverageError(
            "unclassified radical-decreasing families: " + ", ".join(missing)
        )
    library_entries = dict(library_classification or {})
    missing_libraries = sorted(loaded_library_inventory - set(library_entries))
    if missing_libraries:
        raise CoverageError(
            "unclassified radical-decreasing library reactions: "
            + ", ".join(missing_libraries)
        )
    admitted = []
    observed_libraries = set()
    for record in records:
        if int(_field(record, "radical_delta", 0)) >= 0:
            continue
        family = str(_field(record, "family", ""))
        event_id = str(_field(record, "event_id", ""))
        rate_source = _field(record, "rate_source", {}) or {}
        library_id = (
            rate_source.get("library_reaction")
            if isinstance(rate_source, Mapping)
            else None
        )
        if library_id is not None:
            observed_libraries.add(library_id)
            if library_id not in library_entries:
                raise CoverageError(
                    f"unclassified radical-decreasing library reaction: {library_id}"
                )
            library_entry = library_entries[library_id]
            if (
                library_entry.category
                in {CoverageClass.HELD_OUT, CoverageClass.OUT_OF_SCOPE}
                and _field(record, "status", "enabled") == "enabled"
            ):
                raise CoverageError(
                    f"{library_entry.category.value} kinetics-library reaction admitted: {library_id}"
                )
            if library_entry.category is CoverageClass.BIMolecular_MET:
                if int(_field(record, "arity", 0)) != 2:
                    raise CoverageError(
                        f"MET library reaction {library_id} is not bimolecular"
                    )
                admitted.append(event_id)
            elif (
                library_entry.category is CoverageClass.UNIMOLECULAR
                and int(_field(record, "arity", 0)) != 1
            ):
                raise CoverageError(
                    f"unimolecular library reaction {library_id} has wrong arity"
                )
            continue
        if family not in classification:
            raise CoverageError(
                f"compiled radical-decreasing record has unlisted family: {family}"
            )
        entry = classification[family]
        arity = int(_field(record, "arity", 0))
        if entry.category is CoverageClass.BIMolecular_MET:
            if arity != 2:
                raise CoverageError(f"MET family {family} compiled with arity {arity}")
            admitted.append(event_id)
        elif entry.category is CoverageClass.UNIMOLECULAR:
            if arity != 1:
                raise CoverageError(
                    f"unimolecular family {family} compiled with arity {arity}"
                )
        elif _field(record, "status", "enabled") == "enabled":
            raise CoverageError(f"{entry.category.value} family admitted: {family}")
    return CoverageReport(
        MappingProxyType(
            {name: classification[name] for name in sorted(expected_loaded)}
        ),
        tuple(sorted(loaded_library_inventory | observed_libraries)),
        tuple(sorted(admitted)),
    )


def validate_met_topology(component: Any) -> None:
    """Reject topologies outside the pre-registered linear-melt domain."""
    is_gel = bool(_field(component, "is_gel", False))
    branch_points = int(_field(component, "branch_points", 0))
    topology = str(_field(component, "topology", "linear"))
    if is_gel:
        raise UnsupportedTopologyError("MET does not support gel components")
    if branch_points > 0 or topology not in {"linear", "unbranched"}:
        raise UnsupportedTopologyError("MET does not support branched components")


@dataclass(frozen=True)
class MappedGraph:
    """Minimal atom-mapped graph used to verify compiled rewrite behavior."""

    nodes: Mapping[Any, tuple[str, int, int]]
    edges: Mapping[tuple[Any, Any], float]

    def __post_init__(self) -> None:
        nodes = dict(self.nodes)
        edges = {
            tuple(sorted(pair)): float(order) for pair, order in self.edges.items()
        }
        if any(a not in nodes or b not in nodes for a, b in edges):
            raise ValueError("graph edge references an unknown mapped atom")
        object.__setattr__(self, "nodes", MappingProxyType(nodes))
        object.__setattr__(self, "edges", MappingProxyType(edges))

    def canonical(self) -> tuple[Any, ...]:
        adjacency = {label: set() for label in self.nodes}
        for a, b in self.edges:
            adjacency[a].add(b)
            adjacency[b].add(a)
        unseen = set(adjacency)
        components = []
        while unseen:
            stack = [min(unseen)]
            component = set()
            while stack:
                atom = stack.pop()
                if atom in component:
                    continue
                component.add(atom)
                unseen.discard(atom)
                stack.extend(adjacency[atom] - component)
            components.append(tuple(sorted(component)))
        return (
            tuple(sorted((label,) + state for label, state in self.nodes.items())),
            tuple(sorted((a, b, order) for (a, b), order in self.edges.items())),
            tuple(sorted(components)),
        )

    def rewrite(
        self,
        bond_changes: Iterable[tuple[Any, Any, float]],
        node_changes: Iterable[tuple[Any, tuple[str, int, int]]],
    ) -> "MappedGraph":
        nodes = dict(self.nodes)
        edges = dict(self.edges)
        for a, b, order in bond_changes:
            pair = tuple(sorted((a, b)))
            if float(order) == 0.0:
                edges.pop(pair, None)
            else:
                edges[pair] = float(order)
        for label, state in node_changes:
            nodes[label] = state
        return MappedGraph(nodes, edges)

    def relabel(self, permutation: Mapping[Any, Any]) -> "MappedGraph":
        """Return an independently owned graph with mapped atom labels changed."""

        def mapped(label):
            return permutation.get(label, label)

        return MappedGraph(
            {mapped(label): state for label, state in self.nodes.items()},
            {
                tuple(sorted((mapped(a), mapped(b)))): order
                for (a, b), order in self.edges.items()
            },
        )


def apply_junction_operation(
    graph: MappedGraph, operation: Mapping[str, Any]
) -> MappedGraph:
    """Apply the compiler's mapped create/dissociate junction operation."""
    attacked = operation["attacked_atom_label"]
    attacker = operation["attacker_site_label"]
    action = operation["action"]
    order = float(operation["formed_bond"]["order"])
    if action == "create":
        bond_order = order
        radical = 0
    elif action == "dissociate":
        bond_order = 0.0
        radical = 1
    else:
        raise ValueError(f"unsupported junction action: {action}")
    node_changes = []
    for label in (attacked, attacker):
        element, _, hydrogens = graph.nodes[label]
        node_changes.append((label, (element, radical, hydrogens)))
    return graph.rewrite(((attacked, attacker, bond_order),), node_changes)


def graph_postcondition_failures(
    observable: str, actual: MappedGraph, expected: MappedGraph
) -> tuple[str, ...]:
    """Name the one graph postcondition whose canonical behavior differs."""
    return () if actual.canonical() == expected.canonical() else (observable,)


def apply_compiled_graph_rewrite(
    record: Any,
    graph: MappedGraph,
    *,
    atom_labels: Mapping[int, Any] | Sequence[Any] | None = None,
    attacked_atom_id: Any | None = None,
    inverse: bool = False,
) -> MappedGraph:
    """Apply an EventRecord rewrite, resolving compiler indices to persistent IDs."""

    def resolve(atom: Any) -> Any:
        if not isinstance(atom, int):
            return atom
        if atom_labels is None:
            raise ValueError(
                "integer compiler atom indices require persistent atom_labels"
            )
        try:
            return atom_labels[atom]
        except (KeyError, IndexError) as error:
            raise ValueError(
                f"no persistent label for compiler atom index {atom}"
            ) from error

    def apply_operations(operations: Iterable[Mapping[str, Any]]) -> MappedGraph:
        nodes = dict(graph.nodes)
        edges = dict(graph.edges)
        for operation in operations:
            action = operation.get("action")
            if action in {"break", "form"}:
                a, b = (resolve(atom) for atom in operation["atoms"])
                if a not in nodes or b not in nodes:
                    raise ValueError(
                        "compiled rewrite references an unknown persistent atom"
                    )
                pair = tuple(sorted((a, b)))
                if action == "break":
                    edges.pop(pair, None)
                else:
                    edges[pair] = float(operation["order"])
                continue
            atom = resolve(operation.get("atom"))
            if atom not in nodes:
                raise ValueError(
                    f"compiled rewrite references unknown mapped atom {atom}"
                )
            element, radical, hydrogens = nodes[atom]
            if action == "set_radical":
                radical = int(operation["value"])
            elif action == "set_implicit_hydrogens":
                hydrogens = int(operation["value"])
            else:
                raise ValueError(f"unsupported compiled graph operation: {action}")
            nodes[atom] = (element, radical, hydrogens)
        return MappedGraph(nodes, edges)

    junction_ops = _field(record, "junction_ops", ()) or ()
    if junction_ops:
        if len(junction_ops) != 1:
            raise ValueError("a compiled junction record must contain one operation")
        operation = dict(junction_ops[0])
        source_attacked = operation.get("attacked_atom_label")
        declared = str(operation.get("junction_kind", ""))
        channel_id = "J_ortho" if declared.startswith("J_ortho") else declared
        if channel_id not in {"J_para", "J_ortho"}:
            raise ValueError(f"unknown compiled junction channel: {declared}")
        if attacked_atom_id is not None:
            if channel_id != "J_ortho" or attacked_atom_id not in {
                "S6",
                "S7",
            }:
                raise ValueError(
                    "attacked atom override is valid only for J_ortho S6/S7"
                )
            operation["attacked_atom_label"] = attacked_atom_id
            formed_bond = dict(operation["formed_bond"])
            formed_bond["label_pair"] = [
                attacked_atom_id,
                operation["attacker_site_label"],
            ]
            operation["formed_bond"] = formed_bond
        if inverse:
            inverse_operations = [
                dict(item) for item in operation.get("exact_inverse_rewrite", ())
            ]
            if not inverse_operations:
                raise ValueError("junction provenance has no exact inverse rewrite")
            if attacked_atom_id is not None:
                for item in inverse_operations:
                    if "atoms" in item:
                        item["atoms"] = [
                            attacked_atom_id if atom == source_attacked else atom
                            for atom in item["atoms"]
                        ]
                    if item.get("atom") == source_attacked:
                        item["atom"] = attacked_atom_id
            return apply_operations(inverse_operations)
        return apply_junction_operation(graph, operation)
    if inverse:
        raise ValueError("non-junction EventRecord has no compiled inverse rewrite")
    return apply_operations(_field(record, "bond_ops", ()) or ())


def mapped_product_graph(
    record: Any, reactant_heavy_labels: Mapping[int, Any] | Sequence[Any]
) -> MappedGraph:
    """Decode EventRecord.product_graphs using its heavy-atom map."""
    from rmgpy.molecule.molecule import Molecule

    atom_map = {
        int(key): int(value) for key, value in _field(record, "atom_map", {}).items()
    }
    labels = (
        dict(reactant_heavy_labels)
        if isinstance(reactant_heavy_labels, Mapping)
        else dict(enumerate(reactant_heavy_labels))
    )
    if set(labels) != set(atom_map):
        raise ValueError("persistent heavy-atom labels must cover EventRecord.atom_map")
    product_labels = {
        product: labels[reactant] for reactant, product in atom_map.items()
    }
    graphs = _field(record, "product_graphs", ()) or ()
    if not graphs:
        raise ValueError("EventRecord has no product_graphs postcondition")
    nodes: dict[Any, tuple[str, int, int]] = {}
    edges: dict[tuple[Any, Any], float] = {}
    product_index = 0
    for adjacency in graphs:
        molecule = Molecule().from_adjacency_list(adjacency, saturate_h=True)
        heavy = [atom for atom in molecule.atoms if atom.element.number != 1]
        atom_labels_by_object = {}
        for atom in heavy:
            if product_index not in product_labels:
                raise ValueError("product graph has an unmapped heavy atom")
            label = product_labels[product_index]
            atom_labels_by_object[atom] = label
            explicit_h = sum(neighbor.element.number == 1 for neighbor in atom.edges)
            nodes[label] = (
                atom.element.symbol,
                int(atom.radical_electrons),
                int(getattr(atom, "implicit_hydrogens", 0)) + explicit_h,
            )
            product_index += 1
        for atom in heavy:
            for neighbor, bond in atom.edges.items():
                if (
                    neighbor.element.number == 1
                    or atom not in atom_labels_by_object
                    or neighbor not in atom_labels_by_object
                ):
                    continue
                pair = tuple(
                    sorted(
                        (atom_labels_by_object[atom], atom_labels_by_object[neighbor])
                    )
                )
                edges[pair] = float(bond.order)
    if product_index != len(product_labels):
        raise ValueError("EventRecord.atom_map and product_graphs disagree")
    return MappedGraph(nodes, edges)


def compiled_graph_postcondition_failures(
    observable: str,
    record: Any,
    actual: MappedGraph,
    reactant_heavy_labels: Mapping[int, Any] | Sequence[Any],
) -> tuple[str, ...]:
    """Compare an executed rewrite with the EventRecord product graph itself."""
    expected = mapped_product_graph(record, reactant_heavy_labels)
    return graph_postcondition_failures(observable, actual, expected)


def draw_ortho_attack_site(rng: Any) -> str:
    """Draw one of the two symmetry-equivalent ortho sites at no extra rate factor."""
    return str(rng.choice(("S6", "S7")))


@dataclass(frozen=True)
class CompiledChannel:
    channel_id: str
    event_id: str
    family: str
    inventory_class: str
    rate: RateTable
    ssa_multiplier: float
    reversible: bool
    reverse_link: str | None
    reverse_rate: RateTable | None
    continuous_forward: ArrheniusRate | None
    reference_junction: str | None
    rewrite_graph: tuple[Any, ...]
    restrictions: frozenset[str] = field(default_factory=frozenset)

    def forward_rate(self, temperature: float) -> float:
        if self.continuous_forward is not None:
            return self.continuous_forward(temperature)
        return self.rate(temperature) * self.ssa_multiplier

    def dissociation_rate(self, temperature: float) -> float:
        if self.reference_junction is not None:
            return self.forward_rate(temperature) / rmg_reference_equilibrium_constant(
                temperature, self.reference_junction
            )
        if self.reverse_rate is None:
            raise ValueError(f"channel {self.channel_id} has no compiled reverse rate")
        return self.reverse_rate(temperature)

    def equilibrium_constant(self, temperature: float) -> float:
        if self.reference_junction is not None:
            return rmg_reference_equilibrium_constant(
                temperature, self.reference_junction
            )
        if self.reverse_rate is None:
            raise ValueError(
                f"channel {self.channel_id} has no compiled equilibrium object"
            )
        return self.forward_rate(temperature) / self.reverse_rate(temperature)


@dataclass(frozen=True)
class CompiledTerminationTable:
    kernel: str
    arm: TransportArm
    inventory: str
    channels: tuple[CompiledChannel, ...]

    def indexed_channels(self) -> dict[str, CompiledChannel]:
        identifiers = [channel.channel_id for channel in self.channels]
        if len(identifiers) != len(set(identifiers)):
            raise ValueError(f"duplicate channel id in {self.kernel} table")
        return dict(zip(identifiers, self.channels))

    def bulk_rate(
        self,
        temperature: float,
        units_i: float,
        units_j: float,
        pair_class: str,
        *,
        spin_factor: float,
    ) -> float:
        if self.kernel != "bulk":
            raise ValueError("bulk rate is available only from a bulk table")
        activation = sum(
            channel.forward_rate(temperature)
            for channel in self.channels
            if channel.family in {"R_Recombination", "Disproportionation"}
        )
        return collins_kimball(
            activation,
            diffusion_rate(
                self.arm,
                temperature,
                units_i,
                units_j,
                pair_class,
                spin_factor=spin_factor,
            ),
        )


def _validate_inventory(inventory: str) -> None:
    if inventory not in {"R0", "R1"}:
        raise ValueError("inventory must be 'R0' or 'R1'")


def _rewrite_signature(record: Any, channel_id: str) -> tuple[Any, ...]:
    return (
        _freeze(_field(record, "atom_map", {})),
        _freeze(_field(record, "bond_ops", ())),
        _freeze(_field(record, "junction_ops", ())),
        _freeze(_field(record, "product_graphs", ())),
    )


def _freeze(value: Any) -> Any:
    if isinstance(value, Mapping):
        return tuple(sorted((str(key), _freeze(item)) for key, item in value.items()))
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, set):
        return tuple(sorted(_freeze(item) for item in value))
    return value


def _channel_id(record: Any) -> str:
    source = _field(record, "rate_source", {}) or {}
    junction_ops = _field(record, "junction_ops", ()) or ()
    declared_junction = junction_ops[0].get("junction_kind") if junction_ops else None
    source_channel = (
        str(source["channel_id"])
        if isinstance(source, Mapping) and source.get("channel_id")
        else None
    )
    if source_channel is not None and declared_junction is None:
        return source_channel
    table = _field(record, "k_table")
    if not table:
        raise ValueError("compiled MET record has no rate object for channel identity")
    compiled_rate = RateTable.from_mapping(table)
    multiplier = _positive_finite(
        _field(record, "ssa_multiplier", 1.0), "SSA multiplier"
    )
    family = str(_field(record, "family", ""))
    candidates = (
        ("J_para", "J_ortho")
        if _field(record, "inventory_class") == "R1:J_ring"
        and int(_field(record, "arity", 0)) == 2
        else (
            ("null",)
            if family == "R_Recombination"
            else (
                ("D1", "D2", "D3", "D4", "D5") if family == "Disproportionation" else ()
            )
        )
    )
    matches = []
    for channel_id in candidates:
        expected = ARCHIVED_CHANNEL_RATES[channel_id]
        actual_values = [
            compiled_rate(temperature) * multiplier for temperature in RATE_GRID
        ]
        expected_values = [expected(temperature) for temperature in RATE_GRID]
        if _relative_error(actual_values, expected_values) <= RATE_TOLERANCE:
            matches.append(channel_id)
    if len(matches) > 1:
        raise ValueError(
            "compiled MET record does not match exactly one archived channel: "
            f"event={_field(record, 'event_id', '')}, matches={matches}"
        )
    if matches:
        matched = matches[0]
    elif declared_junction is not None:
        declared_channel = (
            "J_ortho"
            if str(declared_junction).startswith("J_ortho")
            else str(declared_junction)
        )
        if declared_channel not in {"J_para", "J_ortho"}:
            raise ValueError(f"unknown compiled junction channel: {declared_junction}")
        raise ValueError(
            f"junction provenance {declared_junction!r} does not match an "
            "archived junction rate"
        )
    elif family in {"R_Recombination", "Disproportionation"}:
        matched = f"{family}:{_field(record, 'event_id', '')}"
    else:
        raise ValueError(
            "compiled MET record has no stable channel identity: "
            f"event={_field(record, 'event_id', '')}"
        )
    if declared_junction is not None:
        declared_channel = (
            "J_ortho"
            if str(declared_junction).startswith("J_ortho")
            else declared_junction
        )
        if declared_channel != matched:
            raise ValueError(
                f"junction provenance {declared_junction!r} disagrees with "
                f"rate-inferred channel {matched}"
            )
        if source_channel is not None and source_channel != matched:
            raise ValueError(
                f"rate-source channel {source_channel!r} disagrees with "
                f"rate-inferred channel {matched}"
            )
    return matched


def _compiled_channel(record: Any, by_id: Mapping[str, Any]) -> CompiledChannel:
    """Construct one owned behavioral row from a compiled forward record."""
    family = str(_field(record, "family", ""))
    inventory_class = str(_field(record, "inventory_class", None) or "base")
    table = _field(record, "k_table")
    reverse_link = _field(record, "reverse_of")
    reverse = by_id.get(str(reverse_link)) if reverse_link else None
    reverse_table = _field(reverse, "k_table") if reverse is not None else None
    forward_rate = RateTable.from_mapping(table)
    reverse_rate = RateTable.from_mapping(reverse_table) if reverse_table else None
    multiplier = _positive_finite(
        _field(record, "ssa_multiplier", 1.0), "SSA multiplier"
    )
    channel_id = _channel_id(record)
    continuous_forward = ARCHIVED_JUNCTION_RATES.get(channel_id)
    if continuous_forward is not None:
        archived_values = [continuous_forward(temperature) for temperature in RATE_GRID]
        compiled_values = [
            forward_rate(temperature) * multiplier for temperature in RATE_GRID
        ]
        if _relative_error(archived_values, compiled_values) > RATE_TOLERANCE:
            continuous_forward = None
    return CompiledChannel(
        channel_id,
        str(_field(record, "event_id", "")),
        family,
        inventory_class,
        forward_rate,
        multiplier,
        reverse is not None,
        str(reverse_link) if reverse_link else None,
        reverse_rate,
        continuous_forward,
        channel_id if continuous_forward is not None else None,
        _rewrite_signature(record, channel_id),
        frozenset(_field(record, "met_restrictions", ()) or ()),
    )


def _is_forward_met_record(record: Any, inventory: str) -> bool:
    if _field(record, "status", "enabled") == "refused":
        return False
    if (
        int(_field(record, "radical_delta", 0)) >= 0
        or int(_field(record, "arity", 0)) != 2
    ):
        return False
    if str(_field(record, "family", "")) not in {
        "R_Recombination",
        "Disproportionation",
    }:
        return False
    if inventory == "R0" and _field(record, "inventory_class") == "R1:J_ring":
        return False
    return bool(_field(record, "k_table"))


def _ortho_forward_label(record: Any) -> str | None:
    junction_ops = _field(record, "junction_ops", ()) or ()
    if len(junction_ops) != 1:
        return None
    operation = junction_ops[0]
    if operation.get("action") != "create" or not str(
        operation.get("junction_kind", "")
    ).startswith("J_ortho"):
        return None
    label = str(operation.get("attacked_atom_label", ""))
    return label if label in {"S6", "S7"} else None


def _compiled_ortho_channel(
    records: Sequence[Any], by_id: Mapping[str, Any], kernel: str
) -> CompiledChannel:
    localisations = {label: [] for label in ("S6", "S7")}
    for record in records:
        label = _ortho_forward_label(record)
        if label is not None and _is_forward_met_record(record, "R1"):
            localisations[label].append(record)
    if any(len(matches) != 1 for matches in localisations.values()):
        raise MissingMETChannelError("J_ortho", kernel)

    representative = localisations["S7"][0]
    rate_tables = [
        RateTable.from_mapping(_field(record, "k_table"))
        for matches in localisations.values()
        for record in matches
    ]
    temperatures = rate_tables[0].temperatures
    if any(table.temperatures != temperatures for table in rate_tables[1:]):
        raise ValueError("J_ortho localisations use different temperature grids")
    rates = tuple(
        sum(
            record_activation_rate(matches[0], temperature)
            for matches in localisations.values()
        )
        for temperature in temperatures
    )
    combined = (
        dict(representative)
        if isinstance(representative, Mapping)
        else dict(vars(representative))
    )
    combined["k_table"] = RateTable(temperatures, rates).to_mapping()
    combined["ssa_multiplier"] = 1.0
    return _compiled_channel(combined, by_id)


def _compile_cage_channels(
    records: Sequence[Any], inventory: str
) -> tuple[CompiledChannel, ...]:
    """Cage-owned compilation path."""
    by_id = {str(_field(record, "event_id", "")): record for record in records}
    channels = []
    for record in records:
        if (
            _is_forward_met_record(record, inventory)
            and _ortho_forward_label(record) is None
        ):
            channels.append(_compiled_channel(record, by_id))
    if inventory == "R1":
        channels.append(_compiled_ortho_channel(records, by_id, "cage"))
    channels.sort(key=lambda channel: channel.channel_id)
    return tuple(channels)


def _compile_bulk_channels(
    records: Sequence[Any], inventory: str
) -> tuple[CompiledChannel, ...]:
    """Bulk-owned compilation path; it never consumes cage rows."""
    by_id = {str(_field(record, "event_id", "")): record for record in records}
    channels = [
        _compiled_channel(record, by_id)
        for record in records
        if _is_forward_met_record(record, inventory)
        and _ortho_forward_label(record) is None
    ]
    if inventory == "R1":
        channels.append(_compiled_ortho_channel(records, by_id, "bulk"))
    return tuple(sorted(channels, key=lambda channel: channel.channel_id))


def _require_r1_channels(channels: Sequence[CompiledChannel], kernel: str) -> None:
    available = {channel.channel_id for channel in channels}
    for required in ("J_para", "J_ortho"):
        if required not in available:
            raise MissingMETChannelError(required, kernel)


def compile_cage_table(
    records: Iterable[Any], arm: TransportArm, inventory: str
) -> CompiledTerminationTable:
    """Compile the cage path directly from records."""
    _validate_inventory(inventory)
    owned_records = tuple(records)
    channels = _compile_cage_channels(owned_records, inventory)
    if inventory == "R1":
        _require_r1_channels(channels, "cage")
    return CompiledTerminationTable("cage", arm, inventory, channels)


def compile_bulk_table(
    records: Iterable[Any], arm: TransportArm, inventory: str
) -> CompiledTerminationTable:
    """Compile the bulk path independently from records."""
    _validate_inventory(inventory)
    owned_records = tuple(records)
    channels = _compile_bulk_channels(owned_records, inventory)
    if inventory == "R1":
        _require_r1_channels(channels, "bulk")
    return CompiledTerminationTable("bulk", arm, inventory, tuple(channels))


def _relative_error(left: Sequence[float], right: Sequence[float]) -> float:
    left_array = np.asarray(left, dtype=float)
    right_array = np.asarray(right, dtype=float)
    if not np.all(np.isfinite(left_array)) or not np.all(np.isfinite(right_array)):
        return math.inf
    return float(np.max(np.abs(left_array / right_array - 1.0)))


def compare_compiled_tables(
    cage: CompiledTerminationTable,
    bulk: CompiledTerminationTable,
    *,
    allowed_restrictions: Mapping[tuple[str, str, str], str] | None = None,
) -> tuple[str, ...]:
    """Return named behavioral mismatches between independently compiled paths."""
    allowed = allowed_restrictions or {}
    failures = []
    arm_id = cage.arm.arm_id
    left_signature = cage.arm.signature()
    right_signature = bulk.arm.signature()
    for index, label in enumerate(("arm_id", "package", "scaling", "Ne", "D0")):
        left = left_signature[index]
        right = right_signature[index]
        if label == "D0":
            mismatch = _relative_error(left, right) > RATE_TOLERANCE
        elif label == "Ne" and left is not None and right is not None:
            mismatch = abs(left / right - 1.0) > RATE_TOLERANCE
        else:
            mismatch = left != right
        if mismatch:
            failures.append(f"transport_{label}:{arm_id}")
    left_channels = cage.indexed_channels()
    right_channels = bulk.indexed_channels()
    for channel_id in sorted(set(left_channels) - set(right_channels)):
        if ("bulk", channel_id, "unavailable") not in allowed:
            failures.append(f"missing_bulk_channel:{channel_id}")
    for channel_id in sorted(set(right_channels) - set(left_channels)):
        if ("cage", channel_id, "unavailable") not in allowed:
            failures.append(f"missing_cage_channel:{channel_id}")
    for channel_id in sorted(set(left_channels) & set(right_channels)):
        left = left_channels[channel_id]
        right = right_channels[channel_id]
        if left.inventory_class != right.inventory_class:
            failures.append(f"inventory_class:{channel_id}")
        left_rates = [
            left.rate(temperature) * left.ssa_multiplier for temperature in RATE_GRID
        ]
        right_rates = [
            right.rate(temperature) * right.ssa_multiplier for temperature in RATE_GRID
        ]
        if _relative_error(left_rates, right_rates) > RATE_TOLERANCE:
            failures.append(f"rate_object:{channel_id}")
        if left.reversible != right.reversible:
            failures.append(f"reversibility:{channel_id}")
        if left.reverse_link != right.reverse_link:
            failures.append(f"reverse_link:{channel_id}")
        for label, evaluator in (
            ("Kc_object", CompiledChannel.equilibrium_constant),
            ("kdiss_object", CompiledChannel.dissociation_rate),
        ):
            if (left.reverse_rate is None) != (right.reverse_rate is None):
                failures.append(f"{label}:{channel_id}")
            elif left.reverse_rate is not None:
                left_values = [
                    evaluator(left, temperature) for temperature in RATE_GRID
                ]
                right_values = [
                    evaluator(right, temperature) for temperature in RATE_GRID
                ]
                if _relative_error(left_values, right_values) > RATE_TOLERANCE:
                    failures.append(f"{label}:{channel_id}")
        if left.rewrite_graph != right.rewrite_graph:
            failures.append(f"rewrite_graph:{channel_id}")
        for restriction in sorted(left.restrictions - right.restrictions):
            if ("cage", channel_id, restriction) not in allowed:
                failures.append(f"unlisted_restriction:{channel_id}:{restriction}")
        for restriction in sorted(right.restrictions - left.restrictions):
            if ("bulk", channel_id, restriction) not in allowed:
                failures.append(f"unlisted_restriction:{channel_id}:{restriction}")
    return tuple(failures)


@dataclass(frozen=True)
class CageObservables:
    terminal_sequence: tuple[str, ...]
    terminal_probabilities: tuple[float, ...]
    occupancies: Mapping[str, float]
    mean_absorption_time: float
    expected_cycles: float
    f_capture: float
    f_redissociate: float
    f_escape: float
    f_permanent: float


@dataclass(frozen=True)
class CageReconciliation:
    """Binding response to a cage-balance result for one inventory arm."""

    status: str
    inventory: str
    fallback_grid: tuple[float, ...] = ()


@dataclass(frozen=True)
class FrozenMetricVerdict:
    """One arm's pre-frozen scientific interpretation."""

    passed: bool
    sign: str
    direction: str
    ranking: tuple[str, ...]
    interpretation: str

    def signature(self) -> tuple[Any, ...]:
        return (self.sign, self.direction, self.ranking, self.interpretation)


@dataclass(frozen=True)
class FallbackStressPoint:
    """Executable R0 full-k_fiss fate split at one pre-registered f_const."""

    f_const: float
    k_fiss: float
    escaped_pair_topology: tuple[str, str]
    escaped_pair_probability: float
    disproportionation_probabilities: Mapping[str, float]


def compile_r0_fallback_stress(
    k_fiss: float, disproportionation_rates: Mapping[str, float]
) -> tuple[FallbackStressPoint, ...]:
    """Build the literal full-k_fiss, no-null, topology-preserving stress grid."""
    k_fiss = _positive_finite(k_fiss, "fission rate")
    rates = {
        product: _positive_finite(rate, f"disproportionation rate {product}")
        for product, rate in disproportionation_rates.items()
    }
    if set(rates) != set(R0_DISPROPORTIONATION_PRODUCTS):
        raise ValueError(
            "R0 fallback requires exactly the two archived disproportionation "
            "product topologies"
        )
    total = sum(rates.values())
    return tuple(
        FallbackStressPoint(
            f_const,
            k_fiss,
            ESCAPED_RADICAL_PAIR,
            f_const,
            MappingProxyType(
                {
                    product: (1.0 - f_const) * rate / total
                    for product, rate in rates.items()
                }
            ),
        )
        for f_const in FALLBACK_STRESS_GRID
    )


def reconcile_cage_balance(
    inventory: str,
    balance_passed: bool,
    *,
    fallback_k_fiss: float | None = None,
    disproportionation_rates: Mapping[str, float] | None = None,
    run_fallback_case: (
        Callable[[FallbackStressPoint], FrozenMetricVerdict] | None
    ) = None,
) -> CageReconciliation:
    """Apply the manager reconciliation without choosing a fallback value."""
    _validate_inventory(inventory)
    if balance_passed:
        return CageReconciliation("valid", inventory)
    if inventory == "R1":
        return CageReconciliation("B4-unresolved", inventory)
    supplied = (
        fallback_k_fiss is not None,
        disproportionation_rates is not None,
        run_fallback_case is not None,
    )
    if not any(supplied):
        return CageReconciliation(
            "fallback-stress-required", inventory, FALLBACK_STRESS_GRID
        )
    if not all(supplied):
        raise ValueError(
            "R0 fallback execution requires k_fiss, archived disproportionation "
            "rates, and a case runner"
        )
    points = compile_r0_fallback_stress(fallback_k_fiss, disproportionation_rates)
    verdicts = [run_fallback_case(point) for point in points]
    if any(not isinstance(verdict, FrozenMetricVerdict) for verdict in verdicts):
        raise TypeError("R0 fallback runner must return FrozenMetricVerdict")
    reference = verdicts[0].signature()
    status = "fallback-stress"
    if any(not verdict.passed for verdict in verdicts) or any(
        verdict.signature() != reference for verdict in verdicts[1:]
    ):
        status = "transport-nonidentified"
    return CageReconciliation(status, inventory, FALLBACK_STRESS_GRID)


def all_arm_conjunction(results: Mapping[str, FrozenMetricVerdict]) -> bool:
    """Enforce pass plus identical sign/direction/ranking/interpretation on A0-A9."""
    if set(results) != set(TRANSPORT_ARMS):
        missing = sorted(set(TRANSPORT_ARMS) - set(results))
        extra = sorted(set(results) - set(TRANSPORT_ARMS))
        raise ValueError(
            f"all-arm conjunction requires A0-A9; missing={missing}, extra={extra}"
        )
    verdicts = [results[arm_id] for arm_id in TRANSPORT_ARMS]
    reference = verdicts[0].signature()
    return all(verdict.passed for verdict in verdicts) and all(
        verdict.signature() == reference for verdict in verdicts[1:]
    )


def solve_cage(
    table: CompiledTerminationTable,
    temperature: float,
    *,
    concentration: float = 1.0,
    escape_rate: float | None = None,
    dissociation_scale: float | Mapping[str, float] = 1.0,
) -> CageObservables:
    """Solve the cage's absorbing CTMC exactly; no SSA is performed here."""
    if table.kernel != "cage":
        raise ValueError("cage observables require a cage table")
    concentration = _positive_finite(concentration, "contact concentration")
    channels = table.indexed_channels()
    terminal_sequence = ("null", "D1", "D2", "D3", "D4", "D5", "escape")
    missing = [channel for channel in terminal_sequence[:-1] if channel not in channels]
    if missing:
        raise ValueError(
            "cage table is missing terminal channels: " + ", ".join(missing)
        )
    terminal_rates = [
        channels[name].forward_rate(temperature) * concentration
        for name in terminal_sequence[:-1]
    ]
    if escape_rate is None:
        escape_rate = (
            4.0 * math.pi * N_A * SIGMA_CONTACT * 2.0 * table.arm.d0(temperature)
        )
    terminal_rates.append(_positive_finite(escape_rate, "escape rate"))
    capture_names = tuple(name for name in ("J_para", "J_ortho") if name in channels)
    capture_rates = [
        channels[name].forward_rate(temperature) * concentration
        for name in capture_names
    ]
    if isinstance(dissociation_scale, Mapping):
        scales = [float(dissociation_scale.get(name, 1.0)) for name in capture_names]
    else:
        scales = [float(dissociation_scale)] * len(capture_names)
    if any(not math.isfinite(scale) or scale < 0.0 for scale in scales):
        raise ValueError("dissociation scale must be finite and nonnegative")
    reverse_rates = [
        channels[name].dissociation_rate(temperature) * scale
        for name, scale in zip(capture_names, scales)
    ]
    if capture_rates and any(rate == 0.0 for rate in reverse_rates):
        required_mass = sum(terminal_rates) / (sum(terminal_rates) + sum(capture_rates))
        raise CageBalanceError(required_mass)
    state_names = ("G",) + capture_names
    size = len(state_names)
    q_matrix = np.zeros((size, size), dtype=float)
    total_g = sum(terminal_rates) + sum(capture_rates)
    q_matrix[0, 0] = -total_g
    for index, (capture, reverse) in enumerate(
        zip(capture_rates, reverse_rates), start=1
    ):
        q_matrix[0, index] = capture
        q_matrix[index, 0] = reverse
        q_matrix[index, index] = -reverse
    absorption = np.zeros((size, len(terminal_sequence)), dtype=float)
    absorption[0, :] = terminal_rates
    fundamental = np.linalg.inv(-q_matrix)
    fates = tuple(float(value) for value in (fundamental @ absorption)[0])
    occupancies = MappingProxyType(
        {name: float(fundamental[0, index]) for index, name in enumerate(state_names)}
    )
    mean_time = float(sum(occupancies.values()))
    capture_probability = sum(capture_rates) / total_g
    cycles = capture_probability / (1.0 - capture_probability)
    permanent = sum(fates[1:6])
    return CageObservables(
        terminal_sequence,
        fates,
        occupancies,
        mean_time,
        cycles,
        capture_probability,
        1.0 if capture_names else 0.0,
        fates[-1],
        permanent,
    )


__all__ = [
    "ARCHIVED_JUNCTION_RATES",
    "ARCHIVED_CHANNEL_RATES",
    "ArrheniusRate",
    "CageReconciliation",
    "CageObservables",
    "CageBalanceError",
    "CompiledChannel",
    "CompiledTerminationTable",
    "CoverageClass",
    "CoverageError",
    "CoverageReport",
    "ESCAPED_RADICAL_PAIR",
    "FAMILY_CLASSIFICATION",
    "FALLBACK_STRESS_GRID",
    "FamilyClassification",
    "FallbackStressPoint",
    "FrozenMetricVerdict",
    "MappedGraph",
    "MissingMETChannelError",
    "PAIR_CLASSES",
    "RATE_GRID",
    "RATE_TOLERANCE",
    "RADICAL_DECREASING_FAMILY_INVENTORY",
    "R0_DISPROPORTIONATION_PRODUCTS",
    "RateTable",
    "SIGMA_CONTACT",
    "TRANSPORT_ARMS",
    "TransportArm",
    "UnsupportedTopologyError",
    "apply_junction_operation",
    "apply_compiled_graph_rewrite",
    "all_arm_conjunction",
    "audit_met_coverage",
    "collins_kimball",
    "compare_compiled_tables",
    "compiled_graph_postcondition_failures",
    "compile_bulk_table",
    "compile_cage_table",
    "compile_r0_fallback_stress",
    "diffusion_rate",
    "draw_ortho_attack_site",
    "enumerate_radical_decreasing_families",
    "enumerate_radical_decreasing_library_reactions",
    "graph_postcondition_failures",
    "mapped_product_graph",
    "record_activation_rate",
    "reconcile_cage_balance",
    "rmg_reference_equilibrium_constant",
    "solve_cage",
    "validate_met_topology",
]
