"""Versioned shared declaration for generic linear homopolymers.

The normalized JSON emitted here is the intermethod contract.  kMC and the
method-of-moments adapter retain these exact bytes instead of independently
interpreting polymer-specific input.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple, Union

from rmgpy import constants
from rmgpy.molecule import Atom, Bond, Molecule
from rmgpy.molecule.adjlist import ConsistencyChecker
from rmgpy.molecule.element import get_element
from rmgpy.polymer_mass import (
    Composition,
    MomentObservables,
    moments_from_observables,
    observables_from_moments,
)


SCHEMA_VERSION = 1
ROUNDING_V1 = "largest_remainder_dp_ascending_v1"
# Where a species lives (small_species.<name>.residence).  "melt": retained in
# the closed melt, as the kinetic Monte Carlo keeps every species; the default
# for every small species, H, H2 and every undeclared light.  "volatile"
# (removed on formation) is a later extension.
SMALL_SPECIES_RESIDENCES = ("melt",)
DEFAULT_RESIDENCE = "melt"

_LOCALITY_BOUND_FIELDS = {
    "elements", "isotopes", "formal_charges", "radical_electrons",
    "lone_pairs", "bond_orders", "spin_multiplicities",
    "maximum_heavy_atom_degree", "maximum_explicit_hydrogens_per_heavy_atom",
    "maximum_total_hydrogens_on_free_small_species", "maximum_family_arity",
    "maximum_root_tuple_size", "maximum_free_small_species_heavy_atoms",
    "attachment_classes", "topology_metadata_schema",
}


def _integer(value: Any, name: str, minimum: Optional[int] = None) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{name} cannot be boolean")
    result = int(value)
    if isinstance(value, str):
        if value.strip() != str(result):
            raise ValueError(f"{name} must be an integer")
    elif result != value:
        raise ValueError(f"{name} must be an integer")
    if minimum is not None and result < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return result


def _finite(value: Any, name: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _canonical_data(value: Any) -> Any:
    """Return JSON-safe data with deterministic mapping order and no NaN."""
    if isinstance(value, Mapping):
        return {
            str(key): _canonical_data(value[key])
            for key in sorted(value, key=lambda item: str(item))
        }
    if isinstance(value, (list, tuple)):
        return [_canonical_data(item) for item in value]
    if isinstance(value, Decimal):
        if not value.is_finite():
            raise ValueError("normalized declaration values must be finite")
        return int(value) if value == value.to_integral_value() else float(value)
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("normalized declaration values must be finite")
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise ValueError(f"value {value!r} is not JSON-normalizable")


def _reject_duplicate_json_pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _normalized_dp_mapping(value: Mapping, name: str) -> Dict[int, Any]:
    normalized = {}
    for raw_dp, item in value.items():
        dp = _integer(raw_dp, "DP", 1)
        if dp in normalized:
            raise ValueError(f"duplicate DP {dp} in {name}")
        normalized[dp] = item
    return normalized


def _validated_locality_policy(value: Any) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or set(value) != {"radius", "bounds"}:
        raise ValueError("locality_policy must declare exactly radius and finite bounds")
    radius = _integer(value["radius"], "locality radius", 0)
    bounds = value["bounds"]
    if not isinstance(bounds, Mapping):
        raise ValueError("locality_policy bounds must be a mapping")
    missing = _LOCALITY_BOUND_FIELDS - set(bounds)
    unknown = set(bounds) - _LOCALITY_BOUND_FIELDS
    if missing:
        raise ValueError(f"missing locality bounds: {sorted(missing)}")
    if unknown:
        raise ValueError(f"unknown locality bounds: {sorted(unknown)}")

    normalized = {}
    for field in ("elements", "attachment_classes"):
        items = bounds[field]
        if not isinstance(items, Sequence) or isinstance(items, str):
            raise ValueError(f"locality bound {field} must be a finite list")
        if any(not isinstance(item, str) or not item for item in items):
            raise ValueError(f"locality bound {field} entries must be non-empty strings")
        if len(items) != len(set(items)):
            raise ValueError(f"locality bound {field} cannot contain duplicates")
        normalized[field] = sorted(items)
    if not normalized["elements"]:
        raise ValueError("locality bound elements cannot be empty")

    integer_lists = ("isotopes", "formal_charges", "radical_electrons",
                     "lone_pairs", "spin_multiplicities")
    for field in integer_lists:
        items = bounds[field]
        if not isinstance(items, Sequence) or isinstance(items, str) or not items:
            raise ValueError(f"locality bound {field} must be a non-empty finite list")
        converted = [_integer(item, f"locality bound {field}") for item in items]
        if len(converted) != len(set(converted)):
            raise ValueError(f"locality bound {field} cannot contain duplicates")
        normalized[field] = sorted(converted)

    if any(item != -1 and item <= 0 for item in normalized["isotopes"]):
        raise ValueError("locality bound isotopes must be -1 or positive integers")
    for field in ("radical_electrons", "lone_pairs"):
        if any(item < 0 for item in normalized[field]):
            raise ValueError(f"locality bound {field} cannot contain negative values")
    if any(item < 1 for item in normalized["spin_multiplicities"]):
        raise ValueError("locality bound spin_multiplicities must be positive")

    bond_orders = bounds["bond_orders"]
    if not isinstance(bond_orders, Sequence) or isinstance(bond_orders, str) or not bond_orders:
        raise ValueError("locality bound bond_orders must be a non-empty finite list")
    normalized["bond_orders"] = sorted({_finite(item, "locality bond order") for item in bond_orders})
    if any(order not in (1.0, 1.5, 2.0, 3.0)
           for order in normalized["bond_orders"]):
        raise ValueError(
            "locality bound bond_orders may contain only 1, 1.5, 2, or 3"
        )

    maxima = (
        "maximum_heavy_atom_degree", "maximum_explicit_hydrogens_per_heavy_atom",
        "maximum_total_hydrogens_on_free_small_species", "maximum_family_arity",
        "maximum_root_tuple_size", "maximum_free_small_species_heavy_atoms",
    )
    for field in maxima:
        normalized[field] = _integer(bounds[field], f"locality bound {field}", 0)
    if normalized["maximum_family_arity"] != 2:
        raise ValueError("locality bound maximum_family_arity must be 2 in schema v1")
    if not isinstance(bounds["topology_metadata_schema"], Mapping):
        raise ValueError("locality bound topology_metadata_schema must be a mapping")
    normalized["topology_metadata_schema"] = _canonical_data(bounds["topology_metadata_schema"])
    return {"radius": radius, "bounds": _canonical_data(normalized)}


@dataclass(frozen=True)
class AtomDeclaration:
    key: str
    element: str
    isotope: int = -1
    charge: int = 0
    radical_electrons: int = 0
    lone_pairs: int = 0
    hydrogens: int = 0

    @classmethod
    def from_mapping(cls, value: Mapping) -> "AtomDeclaration":
        allowed = {"key", "element", "isotope", "charge", "electronic_state", "hydrogens", "lone_pairs"}
        unknown = set(value) - allowed
        if unknown:
            raise ValueError(f"unknown atom fields: {sorted(unknown)}")
        key = value.get("key")
        element = value.get("element")
        if not isinstance(key, str) or not key:
            raise ValueError("atom keys must be non-empty strings")
        if not isinstance(element, str) or not element:
            raise ValueError("atom elements must be non-empty strings")
        isotope = value.get("isotope", -1)
        isotope = -1 if isotope in (None, -1, 0, "") else _integer(isotope, "isotope", 1)
        try:
            physical_element = get_element(element, isotope)
        except Exception as exc:
            raise ValueError(f"unknown element/isotope: {element!r}, isotope={isotope}") from exc
        if physical_element.number <= 0 or physical_element.mass <= 0.0:
            raise ValueError(f"{element!r} is not a physical positive-mass atom")
        state = value.get("electronic_state", {})
        if not isinstance(state, Mapping):
            raise ValueError("electronic_state must be a mapping")
        if set(state) - {"radical_electrons"}:
            raise ValueError("electronic_state only supports radical_electrons in schema v1")
        return cls(
            key=key,
            element=element,
            isotope=isotope,
            charge=_integer(value.get("charge", 0), "formal charge"),
            radical_electrons=_integer(state.get("radical_electrons", 0), "radical electrons", 0),
            lone_pairs=_integer(value.get("lone_pairs", 0), "lone pairs", 0),
            hydrogens=_integer(value.get("hydrogens", 0), "hydrogens", 0),
        )

    def as_mapping(self) -> Dict[str, Any]:
        result = {
            "key": self.key,
            "element": self.element,
            "charge": self.charge,
            "electronic_state": {"radical_electrons": self.radical_electrons},
            "lone_pairs": self.lone_pairs,
            "hydrogens": self.hydrogens,
        }
        if self.isotope != -1:
            result["isotope"] = self.isotope
        return result


@dataclass(frozen=True)
class BondDeclaration:
    atoms: Tuple[str, str]
    order: float

    @classmethod
    def from_mapping(cls, value: Mapping) -> "BondDeclaration":
        if set(value) != {"atoms", "order"}:
            raise ValueError("bonds require exactly atoms and order")
        atoms = value["atoms"]
        if not isinstance(atoms, Sequence) or isinstance(atoms, str) or len(atoms) != 2:
            raise ValueError("bond atoms must contain exactly two atom keys")
        left, right = str(atoms[0]), str(atoms[1])
        if not left or not right or left == right:
            raise ValueError("bond atom keys must be distinct and non-empty")
        order = _finite(value["order"], "bond order")
        if order not in (1.0, 1.5, 2.0, 3.0):
            raise ValueError("schema v1 bond order must be 1, 1.5, 2, or 3")
        return cls(tuple(sorted((left, right))), order)

    def as_mapping(self) -> Dict[str, Any]:
        order = int(self.order) if self.order.is_integer() else self.order
        return {"atoms": list(self.atoms), "order": order}


@dataclass(frozen=True)
class PortDeclaration:
    atom: str
    bond_order: float

    @classmethod
    def from_mapping(cls, value: Mapping) -> "PortDeclaration":
        if set(value) != {"atom", "bond_order"}:
            raise ValueError("ports require exactly atom and bond_order")
        atom = value["atom"]
        if not isinstance(atom, str) or not atom:
            raise ValueError("port atom must be a non-empty atom key")
        order = _finite(value["bond_order"], "port bond order")
        if order not in (1.0, 1.5, 2.0, 3.0):
            raise ValueError("schema v1 port bond order must be 1, 1.5, 2, or 3")
        return cls(atom, order)

    def as_mapping(self) -> Dict[str, Any]:
        order = int(self.bond_order) if self.bond_order.is_integer() else self.bond_order
        return {"atom": self.atom, "bond_order": order}


@dataclass(frozen=True)
class GraphDeclaration:
    atoms: Tuple[AtomDeclaration, ...]
    bonds: Tuple[BondDeclaration, ...]
    ports: Tuple[Tuple[str, PortDeclaration], ...] = ()
    multiplicity: int = 1

    @classmethod
    def from_mapping(
        cls,
        value: Mapping,
        required_ports: Sequence[str] = (),
        single_connection_port: bool = False,
    ) -> "GraphDeclaration":
        if not isinstance(value, Mapping):
            raise ValueError("graph declarations must be mappings")
        allowed = {"atoms", "bonds", "ports", "electronic_state"}
        if set(value) - allowed:
            raise ValueError(f"unknown graph fields: {sorted(set(value) - allowed)}")
        atoms_value = value.get("atoms")
        bonds_value = value.get("bonds", [])
        ports_value = value.get("ports", {})
        electronic_state = value.get("electronic_state")
        if not isinstance(electronic_state, Mapping) \
                or set(electronic_state) != {"multiplicity"}:
            raise ValueError(
                "graph electronic_state must declare exactly multiplicity"
            )
        multiplicity = _integer(
            electronic_state["multiplicity"], "graph multiplicity", 1
        )
        if not isinstance(atoms_value, Sequence) or isinstance(atoms_value, str) or not atoms_value:
            raise ValueError("a physical graph needs a non-empty atom list")
        if not isinstance(bonds_value, Sequence) or isinstance(bonds_value, str):
            raise ValueError("graph bonds must be a list")
        if not isinstance(ports_value, Mapping):
            raise ValueError("graph ports must be a mapping")
        atoms = tuple(sorted(
            (AtomDeclaration.from_mapping(item) for item in atoms_value),
            key=lambda item: item.key,
        ))
        keys = [atom.key for atom in atoms]
        if len(keys) != len(set(keys)):
            raise ValueError("graph atom keys must be unique")
        bonds = tuple(sorted(
            (BondDeclaration.from_mapping(item) for item in bonds_value),
            key=lambda item: (item.atoms, item.order),
        ))
        bond_pairs = [bond.atoms for bond in bonds]
        if len(bond_pairs) != len(set(bond_pairs)):
            raise ValueError("a graph cannot declare a bond pair more than once")
        for bond in bonds:
            if not set(bond.atoms) <= set(keys):
                raise ValueError("bond refers to an unknown atom key")
        ports = tuple(sorted(
            ((str(name), PortDeclaration.from_mapping(item)) for name, item in ports_value.items()),
            key=lambda item: item[0],
        ))
        if any(port.atom not in set(keys) for _name, port in ports):
            raise ValueError("port refers to an unknown atom key")
        if required_ports and {name for name, _port in ports} != set(required_ports):
            raise ValueError(f"graph requires exactly ports {sorted(required_ports)}")
        if single_connection_port and len(ports) != 1:
            raise ValueError("an end-group graph requires exactly one connection port")
        graph = cls(atoms, bonds, ports, multiplicity)
        graph._validate_connected()
        graph._validate_valence()
        return graph

    @property
    def composition(self) -> Composition:
        return Composition.from_atoms([atom.as_mapping() for atom in self.atoms])

    def as_mapping(self, canonical_end_port: bool = False) -> Dict[str, Any]:
        ports = self.ports
        if canonical_end_port and len(ports) == 1:
            ports = (("connection", ports[0][1]),)
        return {
            "atoms": [atom.as_mapping() for atom in self.atoms],
            "bonds": [bond.as_mapping() for bond in self.bonds],
            "ports": {name: port.as_mapping() for name, port in ports},
            "electronic_state": {"multiplicity": self.multiplicity},
        }

    def _validate_connected(self):
        if len(self.atoms) == 1:
            return
        adjacency = {atom.key: set() for atom in self.atoms}
        for bond in self.bonds:
            left, right = bond.atoms
            adjacency[left].add(right)
            adjacency[right].add(left)
        seen = set()
        stack = [self.atoms[0].key]
        while stack:
            key = stack.pop()
            if key not in seen:
                seen.add(key)
                stack.extend(adjacency[key] - seen)
        if len(seen) != len(self.atoms):
            raise ValueError("physical graph atoms must form one connected component")

    def _validate_valence(self):
        maximum = {"H": 1.0, "B": 3.0, "C": 4.0, "N": 4.0, "O": 2.0,
                   "F": 1.0, "Si": 4.0, "P": 5.0, "S": 6.0, "Cl": 1.0,
                   "Br": 1.0, "I": 1.0}
        used = {atom.key: float(atom.hydrogens) for atom in self.atoms}
        for bond in self.bonds:
            for key in bond.atoms:
                used[key] += bond.order
        for _name, port in self.ports:
            used[port.atom] += port.bond_order
        for atom in self.atoms:
            limit = maximum.get(atom.element)
            if limit is not None and used[atom.key] > limit + 1.0e-12:
                raise ValueError(
                    f"joining ports overfills atom {atom.key!r}: valence {used[atom.key]} > {limit}"
                )

    def to_legacy_molecule(self, port_labels: Mapping[str, str]) -> Molecule:
        """Create a disposable radical/label stitching view for legacy MoM."""
        by_key: Dict[str, Atom] = {}
        atoms = []
        port_by_atom: Dict[str, str] = {}
        for name, port in self.ports:
            if name in port_labels:
                port_by_atom[port.atom] = port_labels[name]
        for declaration in self.atoms:
            element = get_element(declaration.element, declaration.isotope)
            atom = Atom(
                element=element,
                radical_electrons=declaration.radical_electrons + (1 if declaration.key in port_by_atom else 0),
                charge=declaration.charge,
                label=port_by_atom.get(declaration.key, ""),
                lone_pairs=declaration.lone_pairs,
            )
            by_key[declaration.key] = atom
            atoms.append(atom)
        molecule = Molecule(atoms=atoms)
        for declaration in self.bonds:
            left, right = declaration.atoms
            molecule.add_bond(Bond(by_key[left], by_key[right], declaration.order))
        for declaration in self.atoms:
            for _index in range(declaration.hydrogens):
                hydrogen = Atom(element=get_element("H"), lone_pairs=0)
                molecule.add_atom(hydrogen)
                molecule.add_bond(Bond(by_key[declaration.key], hydrogen, 1.0))
        molecule.update_multiplicity()
        # Ports add temporary stitching radicals, so their legacy view needs
        # the inferred multiplicity.  A closed physical declaration has no
        # such representation radicals: retain its authoritative electronic
        # state exactly, including singlet diradicals used for mobile binding.
        if not port_by_atom:
            molecule.multiplicity = self.multiplicity
        molecule.update_atomtypes()
        return molecule


@dataclass(frozen=True)
class InitialEnsemble:
    fractions: Tuple[Tuple[int, Fraction], ...]
    representative_chain_count: int
    rounding: str = ROUNDING_V1

    @classmethod
    def from_mapping(cls, value: Mapping) -> "InitialEnsemble":
        if not isinstance(value, Mapping):
            raise ValueError("initial_ensemble must be a mapping")
        allowed = {"dp_number_fractions", "dp_counts", "representative_chain_count", "rounding"}
        if set(value) - allowed:
            raise ValueError(f"unknown initial_ensemble fields: {sorted(set(value) - allowed)}")
        has_fractions = "dp_number_fractions" in value
        has_counts = "dp_counts" in value
        if has_fractions == has_counts:
            raise ValueError("provide exactly one of dp_number_fractions or dp_counts")
        rounding = value.get("rounding", ROUNDING_V1)
        if rounding != ROUNDING_V1:
            raise ValueError(f"rounding must be the literal {ROUNDING_V1!r}")
        if has_counts:
            raw_counts = value["dp_counts"]
            if not isinstance(raw_counts, Mapping) or not raw_counts:
                raise ValueError("dp_counts must be a non-empty mapping")
            raw_counts = _normalized_dp_mapping(raw_counts, "dp_counts")
            counts = {dp: _integer(count, "DP count", 0)
                      for dp, count in raw_counts.items()}
            total = sum(counts.values())
            if total <= 0:
                raise ValueError("dp_counts must contain at least one chain")
            if "representative_chain_count" in value:
                declared = _integer(value["representative_chain_count"], "representative chain count", 1)
                if declared != total:
                    raise ValueError("representative_chain_count must equal sum(dp_counts)")
            fractions = tuple(sorted((dp, Fraction(count, total)) for dp, count in counts.items()))
            ensemble = cls(fractions, total, rounding)
            if ensemble.counts != {dp: count for dp, count in counts.items() if count}:
                raise ValueError("dp_counts do not reconstruct under largest-remainder v1")
            return ensemble

        raw_fractions = value["dp_number_fractions"]
        if not isinstance(raw_fractions, Mapping) or not raw_fractions:
            raise ValueError("dp_number_fractions must be a non-empty mapping")
        raw_fractions = _normalized_dp_mapping(raw_fractions, "dp_number_fractions")
        fractions = tuple(sorted(
            (dp, _fraction(item))
            for dp, item in raw_fractions.items()
        ))
        if any(fraction < 0 for _dp, fraction in fractions):
            raise ValueError("DP number fractions cannot be negative")
        if sum((fraction for _dp, fraction in fractions), Fraction(0, 1)) != 1:
            raise ValueError("DP number fractions must sum exactly to one")
        count = _integer(value.get("representative_chain_count"), "representative chain count", 1)
        return cls(fractions, count, rounding)

    @property
    def counts(self) -> Dict[int, int]:
        floors = {}
        remainders = []
        assigned = 0
        for dp, fraction in self.fractions:
            quota = fraction * self.representative_chain_count
            floor = quota.numerator // quota.denominator
            floors[dp] = floor
            assigned += floor
            remainders.append((quota - floor, dp))
        for _remainder, dp in sorted(remainders, key=lambda item: (-item[0], item[1]))[
            :self.representative_chain_count - assigned
        ]:
            floors[dp] += 1
        return {dp: count for dp, count in sorted(floors.items()) if count}

    def as_mapping(self) -> Dict[str, Any]:
        return {
            "dp_number_fractions": {
                str(dp): {"numerator": fraction.numerator, "denominator": fraction.denominator}
                for dp, fraction in self.fractions
            },
            "representative_chain_count": self.representative_chain_count,
            "rounding": self.rounding,
        }


def _fraction(value: Any) -> Fraction:
    if isinstance(value, Mapping):
        if set(value) != {"numerator", "denominator"}:
            raise ValueError("rational fractions require numerator and denominator")
        numerator = _integer(value["numerator"], "fraction numerator")
        denominator = _integer(value["denominator"], "fraction denominator", 1)
        return Fraction(numerator, denominator)
    if isinstance(value, bool):
        raise ValueError("a fraction cannot be boolean")
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("a decimal fraction must be finite")
        return Fraction(Decimal(str(value)))
    return Fraction(value)


@dataclass(frozen=True)
class SharedInitialization:
    normalized_json: bytes
    dp_counts: Tuple[Tuple[int, int], ...]
    moments_mol: Tuple[float, float, float]
    representative_polymer_mass_kg: float
    representative_mobile_mass_kg: float
    volume_m3: float
    macroscopic_polymer_mass_kg: float
    scale_factor: float


@dataclass(frozen=True)
class LinearHomopolymerSpec:
    schema_version: int
    label: str
    repeat_unit: GraphDeclaration
    head_end_group: GraphDeclaration
    tail_end_group: GraphDeclaration
    monomer_product: Optional[GraphDeclaration]
    small_species: Tuple[Tuple[str, GraphDeclaration, str], ...]
    initial_ensemble: InitialEnsemble
    mn_g_mol: Optional[float]
    mw_g_mol: Optional[float]
    observable_relative_tolerance: float
    initial_mass_kg: float
    density_kg_m3: float
    density_basis: str
    initial_mobile_counts: Tuple[Tuple[str, int], ...]
    locality_policy: Mapping[str, Any]
    family_manifest: Mapping[str, Any]
    library_manifest: Mapping[str, Any]
    transport_model: Mapping[str, Any]
    diagnostics: Mapping[str, Any]

    @classmethod
    def from_mapping(cls, value: Mapping) -> "LinearHomopolymerSpec":
        if not isinstance(value, Mapping):
            raise ValueError("LinearHomopolymerSpec input must be a mapping")
        required = {
            "schema_version", "label", "repeat_unit", "end_groups", "small_species",
            "initial_ensemble", "initial_mass_kg", "density_kg_m3", "locality_policy",
            "family_manifest", "library_manifest", "transport_model", "diagnostics",
        }
        optional = {
            "monomer_product", "Mn_g_mol", "Mw_g_mol", "observable_relative_tolerance",
            "density_basis", "initial_mobile_counts",
        }
        missing = required - set(value)
        unknown = set(value) - required - optional
        if missing:
            raise ValueError(f"missing LinearHomopolymerSpec fields: {sorted(missing)}")
        if unknown:
            raise ValueError(f"unknown LinearHomopolymerSpec fields: {sorted(unknown)}")
        version = _integer(value["schema_version"], "schema_version", 1)
        if version != SCHEMA_VERSION:
            raise ValueError(f"unknown LinearHomopolymerSpec schema_version {version}")
        label = value["label"]
        if not isinstance(label, str) or not label:
            raise ValueError("spec label must be a non-empty display string")
        ends = value["end_groups"]
        if not isinstance(ends, Mapping) or set(ends) != {"head", "tail"}:
            raise ValueError("end_groups must contain exactly head and tail")
        repeat = GraphDeclaration.from_mapping(value["repeat_unit"], required_ports=("head", "tail"))
        if dict(repeat.ports)["head"].atom == dict(repeat.ports)["tail"].atom:
            raise ValueError("repeat-unit head and tail ports must attach to distinct atoms")
        head = GraphDeclaration.from_mapping(ends["head"], single_connection_port=True)
        tail = GraphDeclaration.from_mapping(ends["tail"], single_connection_port=True)
        repeat_ports = dict(repeat.ports)
        if head.ports[0][1].bond_order != repeat_ports["head"].bond_order:
            raise ValueError("head end-group port bond order does not complement repeat head")
        if tail.ports[0][1].bond_order != repeat_ports["tail"].bond_order:
            raise ValueError("tail end-group port bond order does not complement repeat tail")
        monomer_product = value.get("monomer_product")
        if monomer_product is not None:
            monomer_product = GraphDeclaration.from_mapping(monomer_product)
            if monomer_product.ports:
                raise ValueError("monomer_product is a physical species and cannot contain ports")
        small = value["small_species"]
        if not isinstance(small, Mapping):
            raise ValueError("small_species must be a mapping")
        small_species = [_small_species_entry(name, declaration) for name, declaration in sorted(small.items())]
        mobile_counts_value = value.get("initial_mobile_counts", {})
        if not isinstance(mobile_counts_value, Mapping):
            raise ValueError("initial_mobile_counts must be a mapping")
        known_small = {name for name, _graph, _disposition in small_species}
        mobile_counts = tuple(sorted(
            (str(name), _integer(count, "representative mobile count", 0))
            for name, count in mobile_counts_value.items()
        ))
        if not {name for name, _count in mobile_counts} <= known_small:
            raise ValueError("initial_mobile_counts names a species absent from small_species")
        density_basis = value.get("density_basis", "dry_polymer")
        if density_basis not in ("dry_polymer", "whole_initial_mixture"):
            raise ValueError("density_basis must be dry_polymer or whole_initial_mixture")
        if any(count for _name, count in mobile_counts) and density_basis != "whole_initial_mixture":
            raise ValueError("mobile representative counts require whole_initial_mixture density")
        locality = _validated_locality_policy(value["locality_policy"])
        policies = []
        for name in ("family_manifest", "library_manifest", "transport_model", "diagnostics"):
            policy = value[name]
            if not isinstance(policy, Mapping):
                raise ValueError(f"{name} must be a mapping")
            policies.append(_canonical_data(policy))
        mn = value.get("Mn_g_mol")
        mw = value.get("Mw_g_mol")
        mn = None if mn is None else _finite(mn, "Mn_g_mol")
        mw = None if mw is None else _finite(mw, "Mw_g_mol")
        if (mn is None) != (mw is None):
            raise ValueError("Mn_g_mol and Mw_g_mol must be declared together")
        tolerance = _finite(value.get("observable_relative_tolerance", 1.0e-6), "observable tolerance")
        if tolerance < 0.0:
            raise ValueError("observable tolerance cannot be negative")
        spec = cls(
            version, label, repeat, head, tail, monomer_product, tuple(small_species),
            InitialEnsemble.from_mapping(value["initial_ensemble"]), mn, mw, tolerance,
            _finite(value["initial_mass_kg"], "initial_mass_kg"),
            _finite(value["density_kg_m3"], "density_kg_m3"), density_basis,
            mobile_counts, locality, *policies,
        )
        if spec.initial_mass_kg < 0.0:
            raise ValueError("initial_mass_kg cannot be negative")
        if spec.density_kg_m3 <= 0.0:
            raise ValueError("density_kg_m3 must be positive")
        spec._validate_physical_electronic_states()
        spec._validate_locality_coverage()
        spec._validate_observables()
        return spec

    def _validate_locality_coverage(self):
        """Require the finite locality domain to contain every declaration."""
        bounds = self.locality_policy["bounds"]
        graphs = [
            ("repeat_unit", self.repeat_unit, False),
            ("head_end_group", self.head_end_group, False),
            ("tail_end_group", self.tail_end_group, False),
        ]
        if self.monomer_product is not None:
            graphs.append(("monomer_product", self.monomer_product, True))
        graphs.extend(
            (f"small_species[{name}]", graph, True)
            for name, graph, _disposition in self.small_species
        )

        for name, graph, free_small in graphs:
            atom_by_key = {atom.key: atom for atom in graph.atoms}
            heavy_degree = {
                atom.key: 0 for atom in graph.atoms if atom.element != "H"
            }
            explicit_hydrogens = {key: 0 for key in heavy_degree}
            for bond in graph.bonds:
                left, right = bond.atoms
                if left in heavy_degree and atom_by_key[right].element != "H":
                    heavy_degree[left] += 1
                if right in heavy_degree and atom_by_key[left].element != "H":
                    heavy_degree[right] += 1
                if left in heavy_degree and atom_by_key[right].element == "H":
                    explicit_hydrogens[left] += 1
                if right in heavy_degree and atom_by_key[left].element == "H":
                    explicit_hydrogens[right] += 1
                if bond.order not in bounds["bond_orders"]:
                    raise ValueError(
                        f"locality bond_orders exclude {bond.order:g} used by {name}"
                    )
            for _port_name, port in graph.ports:
                if port.atom in heavy_degree:
                    heavy_degree[port.atom] += 1
                if port.bond_order not in bounds["bond_orders"]:
                    raise ValueError(
                        f"locality bond_orders exclude {port.bond_order:g} used by {name}"
                    )

            for atom in graph.atoms:
                if atom.element not in bounds["elements"]:
                    raise ValueError(
                        f"locality elements exclude {atom.element!r} used by {name}"
                    )
                for field, value in (
                    ("isotopes", atom.isotope),
                    ("formal_charges", atom.charge),
                    ("radical_electrons", atom.radical_electrons),
                    ("lone_pairs", atom.lone_pairs),
                ):
                    if value not in bounds[field]:
                        raise ValueError(
                            f"locality {field} exclude {value} used by {name}"
                        )
                if (atom.element != "H" and atom.hydrogens
                        + explicit_hydrogens[atom.key]
                        > bounds["maximum_explicit_hydrogens_per_heavy_atom"]):
                    raise ValueError(
                        "locality maximum_explicit_hydrogens_per_heavy_atom "
                        f"is smaller than {name}'s declaration"
                    )
            if heavy_degree and max(heavy_degree.values()) > bounds["maximum_heavy_atom_degree"]:
                raise ValueError(
                    f"locality maximum_heavy_atom_degree is smaller than {name}'s declaration"
                )

            # Ports are non-atomic topology metadata. Electronic-state
            # coverage therefore follows only the declared physical radical
            # electrons; the capped-dimer consistency check owns the joined
            # molecule separately.
            if graph.multiplicity not in bounds["spin_multiplicities"]:
                raise ValueError(
                    "locality spin_multiplicities exclude "
                    f"{graph.multiplicity} used by {name}"
                )

            if free_small:
                heavy_atoms = sum(atom.element != "H" for atom in graph.atoms)
                total_hydrogens = (
                    sum(atom.hydrogens for atom in graph.atoms)
                    + sum(atom.element == "H" for atom in graph.atoms)
                )
                if heavy_atoms > bounds["maximum_free_small_species_heavy_atoms"]:
                    raise ValueError(
                        "locality maximum_free_small_species_heavy_atoms "
                        f"is smaller than {name}'s declaration"
                    )
                if total_hydrogens > bounds["maximum_total_hydrogens_on_free_small_species"]:
                    raise ValueError(
                        "locality maximum_total_hydrogens_on_free_small_species "
                        f"is smaller than {name}'s declaration"
                    )

    def _validate_physical_electronic_states(self):
        """Apply RMG's valence/electron rules to a physical capped dimer."""
        portless = []
        if self.monomer_product is not None:
            portless.append(("monomer_product", self.monomer_product))
        portless.extend(
            (f"small_species[{name}]", graph)
            for name, graph, _disposition in self.small_species
        )
        for name, graph in portless:
            molecule = graph.to_legacy_molecule({})
            molecule.multiplicity = graph.multiplicity
            try:
                for atom in molecule.atoms:
                    ConsistencyChecker.check_partial_charge(atom)
                    ConsistencyChecker.check_hund_rule(atom, molecule.multiplicity)
                ConsistencyChecker.check_multiplicity(
                    molecule.get_radical_count(), molecule.multiplicity
                )
            except Exception as exc:
                raise ValueError(
                    f"declared electronic state is inconsistent in {name}"
                ) from exc

        fragments = (self.head_end_group, self.repeat_unit,
                     self.repeat_unit, self.tail_end_group)
        atom_maps = []
        molecule = Molecule()
        for fragment in fragments:
            by_key = {}
            for declaration in fragment.atoms:
                atom = Atom(
                    element=get_element(declaration.element, declaration.isotope),
                    radical_electrons=declaration.radical_electrons,
                    charge=declaration.charge,
                    lone_pairs=declaration.lone_pairs,
                )
                molecule.add_atom(atom)
                by_key[declaration.key] = atom
                for _index in range(declaration.hydrogens):
                    hydrogen = Atom(element=get_element("H"), lone_pairs=0)
                    molecule.add_atom(hydrogen)
                    molecule.add_bond(Bond(atom, hydrogen, 1.0))
            for bond in fragment.bonds:
                left, right = bond.atoms
                molecule.add_bond(Bond(by_key[left], by_key[right], bond.order))
            atom_maps.append(by_key)

        head_port = self.head_end_group.ports[0][1]
        repeat_ports = dict(self.repeat_unit.ports)
        tail_port = self.tail_end_group.ports[0][1]
        joins = (
            (0, head_port, 1, repeat_ports["head"]),
            (1, repeat_ports["tail"], 2, repeat_ports["head"]),
            (2, repeat_ports["tail"], 3, tail_port),
        )
        for left_index, left_port, right_index, right_port in joins:
            molecule.add_bond(Bond(
                atom_maps[left_index][left_port.atom],
                atom_maps[right_index][right_port.atom],
                left_port.bond_order,
            ))
        molecule.update_multiplicity()
        try:
            for atom in molecule.atoms:
                ConsistencyChecker.check_partial_charge(atom)
                ConsistencyChecker.check_hund_rule(atom, molecule.multiplicity)
            ConsistencyChecker.check_multiplicity(
                molecule.get_radical_count(), molecule.multiplicity
            )
        except Exception as exc:
            raise ValueError(
                "declared electronic state is inconsistent in the physical capped dimer"
            ) from exc

    @property
    def end_composition(self) -> Composition:
        return self.head_end_group.composition.plus(self.tail_end_group.composition)

    @property
    def normalized_json(self) -> bytes:
        return json.dumps(
            self.as_mapping(), sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")

    @property
    def ensemble_observables(self) -> MomentObservables:
        counts = self.initial_ensemble.counts
        n = sum(counts.values())
        masses = {
            dp: self.repeat_unit.composition.molar_mass_g_mol() * dp
                + self.end_composition.molar_mass_g_mol()
            for dp in counts
        }
        total = sum(counts[dp] * masses[dp] for dp in counts)
        mn = total / n
        mw = sum(counts[dp] * masses[dp] * masses[dp] for dp in counts) / total
        return MomentObservables(total / constants.Na, mn, mw)

    @property
    def macroscopic_moments_mol(self) -> Tuple[float, float, float]:
        observed = self.ensemble_observables
        return moments_from_observables(
            self.initial_mass_kg * 1000.0,
            observed.mn_g_mol,
            observed.mw_g_mol,
            self.repeat_unit.composition.molar_mass_g_mol(),
            self.end_composition.molar_mass_g_mol(),
        )

    def initialization(self) -> SharedInitialization:
        counts = self.initial_ensemble.counts
        polymer_mass_kg = sum(
            count * (
                dp * self.repeat_unit.composition.molar_mass_kg_mol()
                + self.end_composition.molar_mass_kg_mol()
            ) / constants.Na
            for dp, count in counts.items()
        )
        small_by_name = {name: graph for name, graph, _disposition in self.small_species}
        mobile_mass_kg = sum(
            count * small_by_name[name].composition.molar_mass_kg_mol() / constants.Na
            for name, count in self.initial_mobile_counts
        )
        total = polymer_mass_kg + mobile_mass_kg
        scale = self.initial_mass_kg / polymer_mass_kg if polymer_mass_kg else 0.0
        return SharedInitialization(
            normalized_json=self.normalized_json,
            dp_counts=tuple(sorted(counts.items())),
            moments_mol=self.macroscopic_moments_mol,
            representative_polymer_mass_kg=polymer_mass_kg,
            representative_mobile_mass_kg=mobile_mass_kg,
            volume_m3=total / self.density_kg_m3,
            macroscopic_polymer_mass_kg=self.initial_mass_kg,
            scale_factor=scale,
        )

    def for_kmc(self) -> SharedInitialization:
        return self.initialization()

    def for_method_of_moments(self) -> SharedInitialization:
        return self.initialization()

    def legacy_molecules(self) -> Tuple[Molecule, Tuple[Molecule, Molecule]]:
        repeat_ports = dict(self.repeat_unit.ports)
        if any(port.bond_order != 1.0 for port in repeat_ports.values()):
            raise ValueError("the legacy MoM stitch adapter only supports single-bond ports")
        if self.head_end_group.ports[0][1].bond_order != 1.0 or self.tail_end_group.ports[0][1].bond_order != 1.0:
            raise ValueError("the legacy MoM stitch adapter only supports single-bond end ports")
        # The legacy stitcher joins the open *1 of the chain built so far to
        # the *2 of the next piece.  The physical joins are head cap -> repeat
        # head, repeat tail -> next repeat head, repeat tail -> tail cap, so a
        # repeat unit enters through its HEAD (*2) and leaves through its
        # TAIL (*1); the head cap offers *1 and the tail cap receives *2.
        repeat = self.repeat_unit.to_legacy_molecule({"head": "*2", "tail": "*1"})
        head_port_name = self.head_end_group.ports[0][0]
        tail_port_name = self.tail_end_group.ports[0][0]
        head = self.head_end_group.to_legacy_molecule({head_port_name: "*1"})
        tail = self.tail_end_group.to_legacy_molecule({tail_port_name: "*2"})
        return repeat, (head, tail)

    def residence(self, name: str) -> str:
        """Where species ``name`` lives: every small species, H, H2 and every
        undeclared light is melt-resident (``DEFAULT_RESIDENCE``, the only
        residence declared so far; ``small_species.<name>.residence``)."""
        return DEFAULT_RESIDENCE

    def as_mapping(self) -> Dict[str, Any]:
        result = {
            "schema_version": self.schema_version,
            "label": self.label,
            "repeat_unit": self.repeat_unit.as_mapping(),
            "end_groups": {
                "head": self.head_end_group.as_mapping(canonical_end_port=True),
                "tail": self.tail_end_group.as_mapping(canonical_end_port=True),
            },
            "small_species": {
                name: {"graph": graph.as_mapping(), "disposition": disposition}
                for name, graph, disposition in self.small_species
            },
            "initial_ensemble": self.initial_ensemble.as_mapping(),
            "initial_mass_kg": self.initial_mass_kg,
            "density_kg_m3": self.density_kg_m3,
            "density_basis": self.density_basis,
            "initial_mobile_counts": {name: count for name, count in self.initial_mobile_counts},
            "observable_relative_tolerance": self.observable_relative_tolerance,
            "locality_policy": self.locality_policy,
            "family_manifest": self.family_manifest,
            "library_manifest": self.library_manifest,
            "transport_model": self.transport_model,
            "diagnostics": self.diagnostics,
        }
        if self.monomer_product is not None:
            result["monomer_product"] = self.monomer_product.as_mapping()
        if self.mn_g_mol is not None:
            result["Mn_g_mol"] = self.mn_g_mol
            result["Mw_g_mol"] = self.mw_g_mol
        return result

    def _validate_observables(self):
        if self.mn_g_mol is None:
            return
        observed = self.ensemble_observables
        for name, declared, actual in (
            ("Mn_g_mol", self.mn_g_mol, observed.mn_g_mol),
            ("Mw_g_mol", self.mw_g_mol, observed.mw_g_mol),
        ):
            if declared <= 0.0:
                raise ValueError(f"{name} must be positive")
            relative = abs(declared - actual) / max(abs(declared), abs(actual))
            if relative > self.observable_relative_tolerance:
                raise ValueError(
                    f"declared {name}={declared} disagrees with realized ensemble {actual}"
                )


def normalized_linear_homopolymer_spec(value: Any) -> LinearHomopolymerSpec:
    if isinstance(value, LinearHomopolymerSpec):
        return value
    if isinstance(value, (str, bytes, bytearray)):
        value = json.loads(
            value,
            parse_float=Decimal,
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    return LinearHomopolymerSpec.from_mapping(value)


# ---------------------------------------------------------------------------
# schema v2: linear statistical copolymers (``PolymerSpec``)
# ---------------------------------------------------------------------------
# A v2 declaration names residues with typed sites, the links that join a
# site of one class to a site of another, and the sequence support (which
# residue may follow which along a chain).  The support is chemistry: it
# decides which local environments exist.  Composition and sequence
# statistics are the initial state and live in ``initial_ensemble``, outside
# the chemistry digest.  ``LinearHomopolymerSpec`` (v1) is untouched; every
# v1 consumer refuses a v2 declaration by name (``refuse_schema_v2``).
SCHEMA_VERSION_V2 = 2
TOPOLOGY_VOCABULARY = {
    "architecture": ("linear", "branched", "network"),
    "sequence_models": ("ideal_random", "alternating", "markov1"),
    "network_generators": ("random_bond_v1",),
}
SAMPLER_V1 = "markov_chain_v1"
VACANCIES = ("forbidden", "H")
_V1_ENSEMBLE_FIELDS = ("dp_number_fractions", "dp_counts", "representative_chain_count", "rounding")


class SchemaV2Refusal(ValueError):
    """A consumer that reads v1 declarations only was given a v2 one."""

    status = "refused:schema_v2"

    def __init__(self, consumer: str):
        super().__init__(f"{self.status}: {consumer} reads v1 declarations only")
        self.consumer = consumer


def is_schema_v2(value: Any) -> bool:
    """Whether ``value`` is a v2 declaration: a ``PolymerSpec``, or a mapping
    (or its JSON text) declaring schema_version 2."""
    if isinstance(value, (PolymerSpec, BranchedPolymerSpec)):
        return True
    if isinstance(value, (str, bytes, bytearray)):
        try:
            value = json.loads(value)
        except ValueError:
            return False
    if not isinstance(value, Mapping):
        return False
    try:
        return _integer(value.get("schema_version"), "schema_version") == SCHEMA_VERSION_V2
    except (TypeError, ValueError):
        return False


def refuse_schema_v2(value: Any, consumer: str) -> Any:
    """``value`` unchanged, or ``SchemaV2Refusal`` naming ``consumer`` if it
    is a v2 declaration."""
    if is_schema_v2(value):
        raise SchemaV2Refusal(consumer)
    return value


def _exact(value: Any, name: str) -> Fraction:
    """An exact rational: ``{"numerator", "denominator"}`` integers or an
    integer.  Floats are refused, including inside the mapping."""
    if isinstance(value, Mapping) and set(value) == {"numerator", "denominator"}:
        numerator, denominator = value["numerator"], value["denominator"]
        if all(isinstance(x, int) and not isinstance(x, bool) for x in (numerator, denominator)) \
                and denominator > 0:
            return Fraction(numerator, denominator)
    elif isinstance(value, int) and not isinstance(value, bool):
        return Fraction(value)
    raise ValueError(f"{name} must be an exact rational ({{numerator, denominator}} integers, "
                     f"denominator > 0) or an integer, got {value!r}")


def _rational(value: Fraction) -> Dict[str, int]:
    return {"numerator": value.numerator, "denominator": value.denominator}


def _names(value: Any, name: str) -> Tuple[str, ...]:
    if not isinstance(value, Sequence) or isinstance(value, str):
        raise ValueError(f"{name} must be a list of names")
    if any(not isinstance(item, str) or not item for item in value):
        raise ValueError(f"{name} entries must be non-empty strings")
    return tuple(value)


def _link_order(value: Any, name: str) -> float:
    order = _finite(value, name)
    if order not in (1.0, 2.0, 3.0):
        raise ValueError(f"{name} must be 1, 2 or 3 (a vacant end is capped by hydrogens)")
    return order


def _order_data(order: float):
    return int(order) if order.is_integer() else order


@dataclass(frozen=True)
class SiteDeclaration:
    atom: str
    bond_order: float
    site_class: str
    vacancy: str

    @classmethod
    def from_mapping(cls, value: Mapping, name: str, vacancies=("forbidden",)) -> "SiteDeclaration":
        """``vacancies``: the vacancies the architecture admits (a linear chain
        occupies every site; a branched one may leave a site H-capped)."""
        if not isinstance(value, Mapping) or set(value) != {"atom", "bond_order", "class", "vacancy"}:
            raise ValueError(f"site {name} requires exactly atom, bond_order, class and vacancy")
        atom, site_class, vacancy = value["atom"], value["class"], value["vacancy"]
        if not isinstance(atom, str) or not atom:
            raise ValueError(f"site {name}: atom must be a non-empty atom key")
        if not isinstance(site_class, str) or not site_class:
            raise ValueError(f"site {name}: class must be a non-empty string")
        if vacancy not in VACANCIES:
            raise ValueError(f"site {name}: vacancy must be one of {list(VACANCIES)}")
        if vacancy not in vacancies:
            raise ValueError(f"refused:vacancy: site {name} declares vacancy {vacancy!r}; a linear "
                             "chain occupies every site, capped sites arrive with branched chains")
        return cls(atom, _link_order(value["bond_order"], f"site {name} bond_order"), site_class, vacancy)

    def as_mapping(self) -> Dict[str, Any]:
        return {"atom": self.atom, "bond_order": _order_data(self.bond_order),
                "class": self.site_class, "vacancy": self.vacancy}


@dataclass(frozen=True)
class ResidueDeclaration:
    """A residue: its graph with every site occupied (hydrogens counted as
    a v1 repeat unit counts them with its ports joined) and its sites.
    ``graph.ports`` are the sites (atom, bond order), so the v1 graph rules
    (connectivity, valence with every site occupied) apply unchanged."""

    graph: GraphDeclaration
    sites: Tuple[Tuple[str, SiteDeclaration], ...]

    @classmethod
    def from_mapping(cls, value: Mapping, name: str = "residue",
                     vacancies=("forbidden",)) -> "ResidueDeclaration":
        if not isinstance(value, Mapping) or set(value) != {"graph", "sites"}:
            raise ValueError(f"residue {name} requires exactly graph and sites")
        graph, sites = value["graph"], value["sites"]
        if not isinstance(graph, Mapping) or "ports" in graph:
            raise ValueError(f"residue {name}: the graph is a mapping without ports (sites declare them)")
        if not isinstance(sites, Mapping) or not sites:
            raise ValueError(f"residue {name} declares no sites")
        parsed = tuple(sorted((str(site), SiteDeclaration.from_mapping(item, f"{name}.{site}", vacancies))
                              for site, item in sites.items()))
        ports = {site: {"atom": item.atom, "bond_order": item.bond_order} for site, item in parsed}
        return cls(GraphDeclaration.from_mapping(dict(graph, ports=ports)), parsed)

    def as_mapping(self) -> Dict[str, Any]:
        graph = self.graph.as_mapping()
        del graph["ports"]
        return {"graph": graph, "sites": {name: site.as_mapping() for name, site in self.sites}}


@dataclass(frozen=True)
class LinkDeclaration:
    """A link joins a site of class ``source`` on one residue to a site of
    class ``target`` on the next one along a chain."""

    source: str
    target: str
    bond_order: float

    @classmethod
    def from_mapping(cls, value: Mapping, name: str = "link") -> "LinkDeclaration":
        if not isinstance(value, Mapping) or not {"from", "to", "bond_order"} <= set(value) \
                or set(value) - {"from", "to", "bond_order", "bridge"}:
            raise ValueError(f"link {name} requires from, to and bond_order (and optionally bridge)")
        if value.get("bridge") is not None:
            raise ValueError(f"refused:bridge: link {name} declares a bridge graph; bridged links "
                             "arrive with networks")
        source, target = value["from"], value["to"]
        if not all(isinstance(item, str) and item for item in (source, target)):
            raise ValueError(f"link {name}: from and to must be site classes")
        return cls(source, target, _link_order(value["bond_order"], f"link {name} bond_order"))

    def as_mapping(self) -> Dict[str, Any]:
        return {"from": self.source, "to": self.target, "bond_order": _order_data(self.bond_order)}


def _entry_and_exit(residue: ResidueDeclaration, links: Mapping[str, LinkDeclaration]):
    targets = {link.target for link in links.values()}
    sources = {link.source for link in links.values()}
    entries = [(name, site) for name, site in residue.sites if site.site_class in targets]
    exits = [(name, site) for name, site in residue.sites if site.site_class in sources]
    return entries, exits


def _junction_rules(residues, a, b, bonds) -> None:
    """Rules 2 and 2b of the declared chain on the joined pair ``a`` -> ``b``:
    no inter-residue bond on a cycle, no pi system across a link."""
    first, second = residues[a].graph, residues[b].graph
    adjacency, orders = {}, {}

    def connect(x, y, order):
        adjacency.setdefault(x, set()).add(y)
        adjacency.setdefault(y, set()).add(x)
        orders[x] = max(orders.get(x, 0.0), order)
        orders[y] = max(orders.get(y, 0.0), order)

    for side, graph in ((0, first), (1, second)):
        for atom in graph.atoms:
            adjacency.setdefault((side, atom.key), set())
        for bond in graph.bonds:
            connect((side, bond.atoms[0]), (side, bond.atoms[1]), bond.order)
    inter = [((0, x), (1, y), order) for x, y, order in bonds]
    for x, y, order in inter:
        connect(x, y, order)
    for x, y, _order in inter:
        seen, stack = {x}, [x]
        while stack:
            node = stack.pop()
            for nxt in adjacency[node]:
                if {node, nxt} != {x, y} and nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        if y in seen:
            raise ValueError(f"rule 2: the inter-residue bond {a}.{x[1]}-{b}.{y[1]} lies on a cycle "
                             "of the joined pair (ladder or spiro)")
    radical = {(side, atom.key) for side, graph in ((0, first), (1, second))
               for atom in graph.atoms if atom.radical_electrons}
    for x, y, _order in inter:
        if all(node in radical or orders.get(node, 0.0) > 1.0 for node in (x, y)):
            raise ValueError(f"rule 2b: the link {a}.{x[1]}-{b}.{y[1]} joins two atoms of one pi "
                             "system across residues")


def linear_junctions(residues: Mapping[str, ResidueDeclaration], links: Mapping[str, LinkDeclaration],
                     support: Mapping[str, Sequence[str]]) -> Dict[Tuple[str, str], Tuple[str, str, float]]:
    """``(a, b) -> (exit atom of a, entry atom of b, bond order)`` for every
    pair of the support: the one bond a declared link forms when residue
    ``b`` follows residue ``a`` along a chain."""
    by_classes = {}
    for name, link in sorted(links.items()):
        key = (link.source, link.target)
        if key in by_classes:
            raise ValueError(f"links {by_classes[key][0]} and {name} both join {key[0]} to {key[1]}")
        by_classes[key] = (name, link)
    both = sorted({link.source for link in links.values()} & {link.target for link in links.values()})
    if both:
        raise ValueError(f"site classes {both} both leave and enter a link (a network link); a linear "
                         "chain enters each residue through one class and leaves through another")
    classes = {c for key in by_classes for c in key}
    for rname, residue in sorted(residues.items()):
        for sname, site in residue.sites:
            if site.site_class not in classes:
                raise ValueError(f"site {rname}.{sname} of class {site.site_class!r} takes no declared link")
    out = {}
    for a in sorted(support):
        for b in sorted(support[a]):
            _entries, exits = _entry_and_exit(residues[a], links)
            entries, _exits = _entry_and_exit(residues[b], links)
            bonds, used = [], set()
            for sname, site in exits:
                partners = [(tname, target) for tname, target in entries
                            if (site.site_class, target.site_class) in by_classes]
                if len(partners) != 1:
                    raise ValueError(f"junction {a}->{b}: site {a}.{sname} matches {len(partners)} "
                                     f"entry sites of {b}; a junction needs exactly one")
                tname, target = partners[0]
                lname, link = by_classes[(site.site_class, target.site_class)]
                if not site.bond_order == target.bond_order == link.bond_order:
                    raise ValueError(f"junction {a}->{b}: link {lname} has bond order "
                                     f"{link.bond_order:g}, sites {site.bond_order:g} and "
                                     f"{target.bond_order:g}")
                used.add(tname)
                bonds.append((site.atom, target.atom, link.bond_order))
            if not bonds or used != {tname for tname, _site in entries}:
                raise ValueError(f"junction {a}->{b}: the exit sites of {a} and the entry sites of "
                                 f"{b} do not pair one to one through the declared links")
            _junction_rules(residues, a, b, bonds)
            out[(a, b)] = bonds[0]
    return out


def _strongly_connected(nodes, successors) -> bool:
    nodes = sorted(nodes)
    if not nodes:
        return False

    def reach(start, step):
        seen, stack = {start}, [start]
        while stack:
            for nxt in step(stack.pop()):
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        return seen
    predecessors = {n: [m for m in nodes if n in successors(m)] for n in nodes}
    return reach(nodes[0], successors) == set(nodes) == reach(nodes[0], predecessors.__getitem__)


@dataclass(frozen=True)
class SequenceStatistics:
    """``ideal_random`` (Bernoulli: P(next = j) = x_j), ``alternating`` (a
    fixed ``cycle``) or ``markov1`` (first-order ``transition`` rows over the
    residues present, zero entries dropped)."""

    model: str
    cycle: Tuple[str, ...] = ()
    transition: Tuple[Tuple[str, Tuple[Tuple[str, Fraction], ...]], ...] = ()

    def as_mapping(self) -> Dict[str, Any]:
        result: Dict[str, Any] = {"model": self.model}
        if self.model == "alternating":
            result["cycle"] = list(self.cycle)
        if self.model == "markov1":
            result["transition"] = {row: {col: _rational(p) for col, p in entries}
                                    for row, entries in self.transition}
        return result


def _statistics(value, composition: Mapping[str, Fraction],
                support: Mapping[str, Sequence[str]]) -> SequenceStatistics:
    if not isinstance(value, Mapping) or "model" not in value:
        raise ValueError("sequence_statistics must be a mapping with a model")
    model = value["model"]
    present = sorted(r for r, x in composition.items() if x > 0)
    fields = {"ideal_random": {"model"}, "alternating": {"model", "cycle"},
              "markov1": {"model", "transition"}}
    if model not in fields:
        raise ValueError(f"unknown sequence model {model!r}; known: {sorted(fields)}")
    if set(value) != fields[model]:
        raise ValueError(f"sequence model {model} takes exactly {sorted(fields[model])}")
    if model == "ideal_random":
        for i in present:
            for j in present:
                if j not in support[i]:
                    raise ValueError(f"ideal_random needs every pair of present residues in "
                                     f"sequence_support; {i}->{j} is missing")
        return SequenceStatistics(model)
    if model == "alternating":
        cycle = _names(value["cycle"], "alternating cycle")
        if len(cycle) < 2 or set(cycle) - set(composition):
            raise ValueError("an alternating cycle lists at least two declared residues")
        realized = {r: Fraction(cycle.count(r), len(cycle)) for r in composition}
        if realized != dict(composition):
            raise ValueError(f"the alternating cycle {list(cycle)} has composition "
                             f"{ {r: str(x) for r, x in sorted(realized.items())} }, not the declared one")
        for i, j in zip(cycle, cycle[1:] + cycle[:1]):
            if j not in support[i]:
                raise ValueError(f"the alternating cycle step {i}->{j} is outside sequence_support")
        return SequenceStatistics(model, cycle=cycle)
    rows = value["transition"]
    if not isinstance(rows, Mapping) or set(rows) != set(present):
        raise ValueError(f"markov1 transition rows must be exactly the present residues {present}")
    matrix = {}
    for i in present:
        row = rows[i]
        if not isinstance(row, Mapping) or set(row) - set(composition):
            raise ValueError(f"markov1 row {i} must map declared residues to rationals")
        entries = {j: _exact(p, f"transition {i}->{j}") for j, p in row.items()}
        if any(p < 0 for p in entries.values()):
            raise ValueError(f"markov1 row {i} has a negative entry")
        total = sum(entries.values(), Fraction(0))
        if total != 1:
            raise ValueError(f"markov1 row {i} sums to {total}, not 1")
        for j, p in entries.items():
            if p and j not in support[i]:
                raise ValueError(f"markov1 transition {i}->{j} > 0 is outside sequence_support")
            if p and j not in present:
                raise ValueError(f"markov1 transition {i}->{j} > 0 leads to an absent residue")
        matrix[i] = {j: p for j, p in sorted(entries.items()) if p}
    if not _strongly_connected(present, lambda i: matrix[i]):
        raise ValueError("the markov1 chain is not irreducible over the present residues")
    for j in present:
        flow = sum((composition[i] * matrix[i].get(j, Fraction(0)) for i in present), Fraction(0))
        if flow != composition[j]:
            raise ValueError(f"the markov1 stationary distribution differs from the composition at "
                             f"{j} ({flow} != {composition[j]})")
    return SequenceStatistics(model, transition=tuple(
        (i, tuple(matrix[i].items())) for i in present))


@dataclass(frozen=True)
class CopolymerEnsemble:
    """The v1 chain-length ensemble plus composition, sequence statistics,
    sampler and seed (every number an exact rational or integer), and the
    optional ``composition_relative_tolerance``: the largest relative
    deviation of the sampled melt's composition from the declared one that
    the melt kMC accepts (``rmgpy.kmc.melt_sequence``; absent, the kMC
    refuses to sample)."""

    chains: InitialEnsemble
    composition: Tuple[Tuple[str, Fraction], ...]
    statistics: SequenceStatistics
    sampler: str
    seed: int
    composition_tolerance: Optional[Fraction] = None

    @classmethod
    def from_mapping(cls, value: Mapping, residues: Sequence[str],
                     support: Mapping[str, Sequence[str]]) -> "CopolymerEnsemble":
        if not isinstance(value, Mapping):
            raise ValueError("initial_ensemble must be a mapping")
        new = {"composition", "sequence_statistics", "sampler", "seed"}
        unknown = set(value) - set(_V1_ENSEMBLE_FIELDS) - new - {"composition_relative_tolerance"}
        if unknown:
            raise ValueError(f"unknown initial_ensemble fields: {sorted(unknown)}")
        if new - set(value):
            raise ValueError(f"a v2 initial_ensemble requires {sorted(new - set(value))}")
        chains = InitialEnsemble.from_mapping({k: value[k] for k in _V1_ENSEMBLE_FIELDS if k in value})
        raw = value["composition"]
        if not isinstance(raw, Mapping) or set(raw) != set(residues):
            raise ValueError(f"composition must give every residue {sorted(residues)} a mole fraction")
        composition = {r: _exact(x, f"composition {r}") for r, x in raw.items()}
        if any(x < 0 for x in composition.values()):
            raise ValueError("composition mole fractions cannot be negative")
        if sum(composition.values(), Fraction(0)) != 1:
            raise ValueError("composition mole fractions must sum exactly to one")
        statistics = _statistics(value["sequence_statistics"], composition, support)
        if value["sampler"] != SAMPLER_V1:
            raise ValueError(f"sampler must be the literal {SAMPLER_V1!r}")
        seed = value["seed"]
        if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
            raise ValueError(f"seed must be a non-negative integer, got {seed!r}")
        tolerance = value.get("composition_relative_tolerance")
        if tolerance is not None:
            tolerance = _exact(tolerance, "composition_relative_tolerance")
            if tolerance < 0:
                raise ValueError("composition_relative_tolerance cannot be negative")
        return cls(chains, tuple(sorted(composition.items())), statistics, SAMPLER_V1, seed, tolerance)

    def as_mapping(self) -> Dict[str, Any]:
        result = dict(self.chains.as_mapping(),
                      composition={r: _rational(x) for r, x in self.composition},
                      sequence_statistics=self.statistics.as_mapping(),
                      sampler=self.sampler, seed=self.seed)
        if self.composition_tolerance is not None:
            result["composition_relative_tolerance"] = _rational(self.composition_tolerance)
        return result


def _capped_chain_molecule(pieces, joins) -> Molecule:
    """``pieces`` (graphs) joined by ``joins`` ((piece, atom key, piece,
    atom key, order)), every declared hydrogen explicit."""
    molecule, maps = Molecule(), []
    for graph in pieces:
        by_key = {}
        for declaration in graph.atoms:
            atom = Atom(element=get_element(declaration.element, declaration.isotope),
                        radical_electrons=declaration.radical_electrons,
                        charge=declaration.charge, lone_pairs=declaration.lone_pairs)
            molecule.add_atom(atom)
            by_key[declaration.key] = atom
            for _index in range(declaration.hydrogens):
                hydrogen = Atom(element=get_element("H"), lone_pairs=0)
                molecule.add_atom(hydrogen)
                molecule.add_bond(Bond(atom, hydrogen, 1.0))
        for bond in graph.bonds:
            molecule.add_bond(Bond(by_key[bond.atoms[0]], by_key[bond.atoms[1]], bond.order))
        maps.append(by_key)
    for left, left_key, right, right_key, order in joins:
        molecule.add_bond(Bond(maps[left][left_key], maps[right][right_key], order))
    molecule.update_multiplicity()
    return molecule


def _check_electronic_state(molecule: Molecule, where: str) -> None:
    try:
        for atom in molecule.atoms:
            ConsistencyChecker.check_partial_charge(atom)
            ConsistencyChecker.check_hund_rule(atom, molecule.multiplicity)
        ConsistencyChecker.check_multiplicity(molecule.get_radical_count(), molecule.multiplicity)
    except Exception as exc:
        raise ValueError(f"declared electronic state is inconsistent in {where}") from exc


def _small_species_entry(name, declaration) -> Tuple[str, GraphDeclaration, str]:
    """One ``small_species`` entry: graph, disposition and the optional
    ``residence`` (``SMALL_SPECIES_RESIDENCES``; the default is canonical and
    not stored, so stating it changes no digest)."""
    if not isinstance(name, str) or not name:
        raise ValueError("small-species labels must be non-empty strings")
    if not isinstance(declaration, Mapping) or not (
            {"graph", "disposition"} <= set(declaration) <= {"graph", "disposition", "residence"}):
        raise ValueError("each small species requires graph and disposition (and may declare residence)")
    residence = declaration.get("residence", DEFAULT_RESIDENCE)
    if residence == "volatile":
        raise ValueError(f"small species {name!r}: residence 'volatile' (removed from the melt on formation) "
                         "is a later extension")
    if residence not in SMALL_SPECIES_RESIDENCES:
        raise ValueError(f"small species {name!r}: residence {residence!r} is not one of "
                         f"{SMALL_SPECIES_RESIDENCES}")
    graph = GraphDeclaration.from_mapping(declaration["graph"])
    if graph.ports:
        raise ValueError("small physical species cannot contain connection ports")
    disposition = declaration["disposition"]
    if not isinstance(disposition, str) or not disposition:
        raise ValueError("small-species disposition must be a non-empty string")
    return name, graph, disposition


def _small_species(value) -> Tuple[Tuple[str, GraphDeclaration, str], ...]:
    if not isinstance(value, Mapping):
        raise ValueError("small_species must be a mapping")
    return tuple(_small_species_entry(name, declaration) for name, declaration in sorted(value.items()))


def _v1_locality(locality, repeat, head, tail, monomer=None, small=()) -> None:
    """The v1 locality-coverage rules on one residue (as the repeat unit),
    the end groups and the free species."""
    from types import SimpleNamespace
    LinearHomopolymerSpec._validate_locality_coverage(SimpleNamespace(
        locality_policy=locality, repeat_unit=repeat, head_end_group=head, tail_end_group=tail,
        monomer_product=monomer, small_species=small))


@dataclass(frozen=True)
class PolymerSpec:
    schema_version: int
    label: str
    architecture: str
    residues: Tuple[Tuple[str, ResidueDeclaration], ...]
    links: Tuple[Tuple[str, LinkDeclaration], ...]
    sequence_support: Tuple[Tuple[str, Tuple[str, ...]], ...]
    head_end_group: GraphDeclaration
    tail_end_group: GraphDeclaration
    monomer_products: Tuple[Tuple[str, GraphDeclaration], ...]
    small_species: Tuple[Tuple[str, GraphDeclaration, str], ...]
    initial_ensemble: CopolymerEnsemble
    initial_mass_kg: float
    density_kg_m3: float
    density_basis: str
    initial_mobile_counts: Tuple[Tuple[str, int], ...]
    locality_policy: Mapping[str, Any]
    family_manifest: Mapping[str, Any]
    library_manifest: Mapping[str, Any]
    transport_model: Mapping[str, Any]
    diagnostics: Mapping[str, Any]
    junction_table: Tuple[Tuple[Tuple[str, str], Tuple[str, str, float]], ...]

    @property
    def junctions(self) -> Dict[Tuple[str, str], Tuple[str, str, float]]:
        """``(a, b) -> (exit atom of a, entry atom of b, bond order)``."""
        return dict(self.junction_table)

    @classmethod
    def from_mapping(cls, value: Mapping) -> "PolymerSpec":
        if not isinstance(value, Mapping):
            raise ValueError("PolymerSpec input must be a mapping")
        version = _integer(value.get("schema_version"), "schema_version", 1)
        if version != SCHEMA_VERSION_V2:
            raise ValueError(f"PolymerSpec reads schema_version {SCHEMA_VERSION_V2}, got {version} "
                             "(schema_version 1 is a LinearHomopolymerSpec)")
        required = {
            "schema_version", "label", "architecture", "residues", "links", "sequence_support",
            "end_groups", "small_species", "initial_ensemble", "initial_mass_kg", "density_kg_m3",
            "locality_policy", "family_manifest", "library_manifest", "transport_model", "diagnostics",
        }
        optional = {"monomer_products", "density_basis", "initial_mobile_counts"}
        missing = required - set(value)
        unknown = set(value) - required - optional
        if missing:
            raise ValueError(f"missing PolymerSpec fields: {sorted(missing)}")
        if unknown:
            raise ValueError(f"unknown PolymerSpec fields: {sorted(unknown)}")
        label = value["label"]
        if not isinstance(label, str) or not label:
            raise ValueError("spec label must be a non-empty display string")
        architecture = value["architecture"]
        if architecture not in TOPOLOGY_VOCABULARY["architecture"]:
            raise ValueError(f"unknown architecture {architecture!r}")
        if architecture != "linear":
            raise ValueError(f"refused:architecture: {architecture} chains are not declarable yet; "
                             "schema v2 declares linear chains")

        raw_residues = value["residues"]
        if not isinstance(raw_residues, Mapping) or not raw_residues:
            raise ValueError("residues must be a non-empty mapping")
        residues = {str(name): ResidueDeclaration.from_mapping(item, str(name))
                    for name, item in raw_residues.items()}
        if len(residues) == 1:
            raise ValueError("a linear chain of one residue is a homopolymer: declare it as v1 "
                             "(schema_version 1, linear_homopolymer_spec)")
        raw_links = value["links"]
        if not isinstance(raw_links, Mapping) or not raw_links:
            raise ValueError("links must be a non-empty mapping")
        links = {str(name): LinkDeclaration.from_mapping(item, str(name)) for name, item in raw_links.items()}

        locality = _validated_locality_policy(value["locality_policy"])
        bounds = locality["bounds"]
        classes = set(bounds["attachment_classes"])
        for rname, residue in sorted(residues.items()):
            for sname, site in residue.sites:
                if site.site_class not in classes:
                    raise ValueError(f"site class {site.site_class!r} of {rname}.{sname} is not in "
                                     "locality_policy.bounds.attachment_classes")
        for lname, link in sorted(links.items()):
            for site_class in (link.source, link.target):
                if site_class not in classes:
                    raise ValueError(f"link {lname} class {site_class!r} is not in "
                                     "locality_policy.bounds.attachment_classes")
        vocabulary = bounds["topology_metadata_schema"]
        if set(vocabulary) - set(TOPOLOGY_VOCABULARY):
            raise ValueError(f"unknown topology_metadata_schema keys: "
                             f"{sorted(set(vocabulary) - set(TOPOLOGY_VOCABULARY))}")
        for key, names in vocabulary.items():
            if set(_names(names, f"topology_metadata_schema {key}")) - set(TOPOLOGY_VOCABULARY[key]):
                raise ValueError(f"topology_metadata_schema {key} names outside "
                                 f"{list(TOPOLOGY_VOCABULARY[key])}")
        if architecture not in vocabulary.get("architecture", ()):
            raise ValueError(f"architecture {architecture!r} is not in "
                             "locality_policy.bounds.topology_metadata_schema")

        raw_support = value["sequence_support"]
        if not isinstance(raw_support, Mapping) or set(raw_support) != set(residues):
            raise ValueError(f"sequence_support must list the successors of every residue {sorted(residues)}")
        support = {}
        for name, successors in raw_support.items():
            successors = _names(successors, f"sequence_support {name}")
            if len(successors) != len(set(successors)) or set(successors) - set(residues):
                raise ValueError(f"sequence_support {name} must name distinct declared residues")
            support[name] = tuple(sorted(successors))
        if not _strongly_connected(residues, support.__getitem__):
            raise ValueError("sequence_support must be strongly connected: every residue reachable "
                             "from every other along a chain")
        junctions = linear_junctions(residues, links, support)
        for name, residue in sorted(residues.items()):
            entries, exits = _entry_and_exit(residue, links)
            if len(entries) != 1 or len(exits) != 1:
                raise ValueError(f"residue {name}: a linear residue has exactly one entry site and "
                                 "one exit site")

        ends = value["end_groups"]
        if not isinstance(ends, Mapping) or set(ends) != {"head", "tail"}:
            raise ValueError("end_groups must contain exactly head and tail")
        head = GraphDeclaration.from_mapping(ends["head"], single_connection_port=True)
        tail = GraphDeclaration.from_mapping(ends["tail"], single_connection_port=True)
        for name, residue in sorted(residues.items()):
            entries, exits = _entry_and_exit(residue, links)
            if head.ports[0][1].bond_order != entries[0][1].bond_order:
                raise ValueError(f"head end-group port bond order does not complement {name}'s entry site")
            if tail.ports[0][1].bond_order != exits[0][1].bond_order:
                raise ValueError(f"tail end-group port bond order does not complement {name}'s exit site")

        raw_products = value.get("monomer_products", {})
        if not isinstance(raw_products, Mapping) or set(raw_products) - set(residues):
            raise ValueError("monomer_products maps declared residues to physical species")
        products = []
        for name, graph in sorted(raw_products.items()):
            graph = GraphDeclaration.from_mapping(graph)
            if graph.ports:
                raise ValueError("a monomer product is a physical species and cannot contain ports")
            products.append((name, graph))
        small_species = _small_species(value["small_species"])
        mobile = value.get("initial_mobile_counts", {})
        if not isinstance(mobile, Mapping):
            raise ValueError("initial_mobile_counts must be a mapping")
        mobile_counts = tuple(sorted((str(name), _integer(count, "representative mobile count", 0))
                                     for name, count in mobile.items()))
        if not {name for name, _count in mobile_counts} <= {name for name, _g, _d in small_species}:
            raise ValueError("initial_mobile_counts names a species absent from small_species")
        density_basis = value.get("density_basis", "dry_polymer")
        if density_basis not in ("dry_polymer", "whole_initial_mixture"):
            raise ValueError("density_basis must be dry_polymer or whole_initial_mixture")
        if any(count for _name, count in mobile_counts) and density_basis != "whole_initial_mixture":
            raise ValueError("mobile representative counts require whole_initial_mixture density")
        ensemble = CopolymerEnsemble.from_mapping(value["initial_ensemble"], tuple(residues), support)
        if ensemble.statistics.model not in vocabulary.get("sequence_models", ()):
            raise ValueError(f"sequence model {ensemble.statistics.model!r} is not in "
                             "locality_policy.bounds.topology_metadata_schema")
        policies = []
        for name in ("family_manifest", "library_manifest", "transport_model", "diagnostics"):
            if not isinstance(value[name], Mapping):
                raise ValueError(f"{name} must be a mapping")
            policies.append(_canonical_data(value[name]))
        spec = cls(
            SCHEMA_VERSION_V2, label, architecture, tuple(sorted(residues.items())),
            tuple(sorted(links.items())), tuple(sorted(support.items())), head, tail, tuple(products),
            small_species, ensemble, _finite(value["initial_mass_kg"], "initial_mass_kg"),
            _finite(value["density_kg_m3"], "density_kg_m3"), density_basis, mobile_counts, locality,
            *policies, tuple(sorted(junctions.items())),
        )
        if spec.initial_mass_kg < 0.0:
            raise ValueError("initial_mass_kg cannot be negative")
        if spec.density_kg_m3 <= 0.0:
            raise ValueError("density_kg_m3 must be positive")
        spec._validate_coverage_and_states()
        return spec

    def _validate_coverage_and_states(self) -> None:
        residues = dict(self.residues)
        for _name, residue in self.residues:
            _v1_locality(self.locality_policy, residue.graph, self.head_end_group, self.tail_end_group)
        first = self.residues[0][1].graph
        for _name, product in self.monomer_products:
            _v1_locality(self.locality_policy, first, self.head_end_group, self.tail_end_group,
                         monomer=product)
        _v1_locality(self.locality_policy, first, self.head_end_group, self.tail_end_group,
                     small=self.small_species)
        for name, graph in list(self.monomer_products) + [(n, g) for n, g, _d in self.small_species]:
            molecule = graph.to_legacy_molecule({})
            molecule.multiplicity = graph.multiplicity
            _check_electronic_state(molecule, f"species {name}")
        links = dict(self.links)
        sites = {name: _entry_and_exit(residue, links) for name, residue in self.residues}
        _head_name, head_port = self.head_end_group.ports[0]
        _tail_name, tail_port = self.tail_end_group.ports[0]
        for (a, b), (exit_atom, entry_atom, order) in self.junction_table:
            molecule = _capped_chain_molecule(
                (self.head_end_group, residues[a].graph, residues[b].graph, self.tail_end_group),
                ((0, head_port.atom, 1, sites[a][0][0][1].atom, head_port.bond_order),
                 (1, exit_atom, 2, entry_atom, order),
                 (2, sites[b][1][0][1].atom, 3, tail_port.atom, tail_port.bond_order)))
            _check_electronic_state(molecule, f"the capped pair {a}-{b}")

    @property
    def normalized_json(self) -> bytes:
        return json.dumps(
            self.as_mapping(), sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")

    def residence(self, name: str) -> str:
        """Where species ``name`` lives: every small species, H, H2 and every
        undeclared light is melt-resident (``DEFAULT_RESIDENCE``, the only
        residence declared so far; ``small_species.<name>.residence``)."""
        return DEFAULT_RESIDENCE

    def as_mapping(self) -> Dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "label": self.label,
            "architecture": self.architecture,
            "residues": {name: residue.as_mapping() for name, residue in self.residues},
            "links": {name: link.as_mapping() for name, link in self.links},
            "sequence_support": {name: list(successors) for name, successors in self.sequence_support},
            "end_groups": {
                "head": self.head_end_group.as_mapping(canonical_end_port=True),
                "tail": self.tail_end_group.as_mapping(canonical_end_port=True),
            },
            "monomer_products": {name: graph.as_mapping() for name, graph in self.monomer_products},
            "small_species": {
                name: {"graph": graph.as_mapping(), "disposition": disposition}
                for name, graph, disposition in self.small_species
            },
            "initial_ensemble": self.initial_ensemble.as_mapping(),
            "initial_mass_kg": self.initial_mass_kg,
            "density_kg_m3": self.density_kg_m3,
            "density_basis": self.density_basis,
            "initial_mobile_counts": {name: count for name, count in self.initial_mobile_counts},
            "locality_policy": self.locality_policy,
            "family_manifest": self.family_manifest,
            "library_manifest": self.library_manifest,
            "transport_model": self.transport_model,
            "diagnostics": self.diagnostics,
        }


def normalized_polymer_spec(value: Any) -> Union["PolymerSpec", "BranchedPolymerSpec"]:
    """A v2 declaration from a ``PolymerSpec`` / ``BranchedPolymerSpec``, a
    mapping or its JSON text: ``architecture`` ``branched`` reads as a
    ``BranchedPolymerSpec``, ``linear`` as a ``PolymerSpec``; ``network`` is
    refused by name.  JSON floats parse as floats, so ``normalized_json`` is
    a fixed point."""
    if isinstance(value, (PolymerSpec, BranchedPolymerSpec)):
        return value
    if isinstance(value, (str, bytes, bytearray)):
        value = json.loads(value, object_pairs_hook=_reject_duplicate_json_pairs)
    architecture = value.get("architecture") if isinstance(value, Mapping) else None
    if architecture == "network":
        raise ValueError("refused:architecture: network chains are not declarable yet (S7); schema v2 "
                         "declares linear and branched chains")
    if architecture == "branched":
        return BranchedPolymerSpec.from_mapping(value)
    return PolymerSpec.from_mapping(value)


# ---------------------------------------------------------------------------
# schema v2, architecture "branched": finite branched molecules below the gel
# point (``BranchedPolymerSpec``)
# ---------------------------------------------------------------------------
# Residues carry sites that a link may occupy or leave vacant (``vacancy: H``:
# an H cap) and sites that are always linked (``forbidden``: a linking
# residue).  A link joins a site of class ``from`` to a site of class ``to``
# on two different residues; which residue follows which is pure site
# occupancy, so there is no sequence support and no end group (a chain end is
# a vacant site).  ``minimum_cycle_bonds`` (chemistry) is the shortest cycle
# through inter-residue bonds any molecule may hold; at 2(r+1)+2 or more every
# environment is a tree (an environment reads r + 1 bonds on each side of its
# atom), so shorter values are refused.  The initial state is the seeded
# ``random_bond_v1`` matching of the declared link counts over the declared
# residue counts (``initial_ensemble``, outside the chemistry digest);
# ``rmgpy.kmc.melt_branched`` builds it.  Counts at or above the gel point are
# refused (networks are S7): the gel point is where the half-edge branching
# process of the random matching becomes critical (``gel_criterion``).
GENERATOR_V1 = "random_bond_v1"
_BRANCHED_ENSEMBLE_FIELDS = {"generator", "representative_residue_counts", "link_counts", "seed"}


def default_minimum_cycle_bonds(radius: int) -> int:
    """2(r + 1) + 2: the shortest cycle no environment of radius ``radius``
    (with its boundary shell) can hold whole."""
    return 2 * (int(radius) + 1) + 2


@dataclass(frozen=True)
class BranchedEnsemble:
    """``random_bond_v1``: residue counts, link counts (links of each type)
    and the seed, exact integers."""

    generator: str
    residue_counts: Tuple[Tuple[str, int], ...]
    link_counts: Tuple[Tuple[str, int], ...]
    seed: int

    @classmethod
    def from_mapping(cls, value: Mapping, residues: Sequence[str], links: Sequence[str]) -> "BranchedEnsemble":
        if not isinstance(value, Mapping) or set(value) != _BRANCHED_ENSEMBLE_FIELDS:
            raise ValueError(f"a branched initial_ensemble has exactly {sorted(_BRANCHED_ENSEMBLE_FIELDS)}")
        if value["generator"] != GENERATOR_V1:
            raise ValueError(f"generator must be the literal {GENERATOR_V1!r}")
        out = []
        for field, names in (("representative_residue_counts", residues), ("link_counts", links)):
            raw = value[field]
            if not isinstance(raw, Mapping) or set(raw) != set(names):
                raise ValueError(f"{field} must give every declared name {sorted(names)} a count")
            out.append(tuple(sorted((str(name), _integer(count, f"{field} {name}", 0))
                                    for name, count in raw.items())))
        seed = value["seed"]
        if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
            raise ValueError(f"seed must be a non-negative integer, got {seed!r}")
        return cls(GENERATOR_V1, out[0], out[1], seed)

    def as_mapping(self) -> Dict[str, Any]:
        return {"generator": self.generator,
                "representative_residue_counts": dict(self.residue_counts),
                "link_counts": dict(self.link_counts), "seed": self.seed}


def _class_tables(residues, links, residue_counts, link_counts):
    """Per site class: (sites in the melt, link endpoints), and per pair of
    classes the endpoints that a link of that pair joins (both directions)."""
    sites, ends, pairs = {}, {}, {}
    for name, residue in residues:
        for _s, site in residue.sites:
            sites[site.site_class] = sites.get(site.site_class, 0) + residue_counts[name]
    for name, link in links:
        n = link_counts[name]
        for x, y in ((link.source, link.target), (link.target, link.source)):
            ends[x] = ends.get(x, 0) + n
            pairs[(x, y)] = pairs.get((x, y), 0) + n
    return sites, ends, pairs


def gel_criterion(residues, links, residue_counts: Mapping[str, int],
                  link_counts: Mapping[str, int]) -> Tuple[bool, list]:
    """(below the gel point, the half-edge branching matrix A) of a random
    matching of the declared links over the declared residues.

    A half edge arrives at a site of class c; it lands on residue R with
    probability N_R k_Rc / S_c (k_Rc: R's sites of class c, S_c: all sites of
    class c), every other site s' of R is linked with probability p = E/S of
    its class (E: link endpoints), and a link from class c' arrives at class
    c'' with probability q(c' -> c'') (its share of c' endpoints).  So
    A[c][c''] = sum_R P(R | c) sum_{s' != the arrival site} p(c') q(c' -> c'').
    Finite molecules (below the gel point) iff the spectral radius of A is
    below 1, iff I - A is a nonsingular M-matrix, iff every leading principal
    minor of I - A is positive (exact in rationals).  Equal reactivity on three
    sites per residue gives Flory's alpha_c = 1/2."""
    sites, ends, pairs = _class_tables(residues, links, residue_counts, link_counts)
    classes = sorted(c for c in ends if ends[c] and sites.get(c))
    index = {c: i for i, c in enumerate(classes)}
    p = {c: Fraction(ends[c], sites[c]) for c in classes}
    matrix = [[Fraction(0)] * len(classes) for _ in classes]
    for c in classes:
        for name, residue in residues:
            k = {}
            for _s, site in residue.sites:
                k[site.site_class] = k.get(site.site_class, 0) + 1
            if not k.get(c) or not residue_counts[name]:
                continue
            land = Fraction(residue_counts[name] * k[c], sites[c])
            for other, count in k.items():
                count -= other == c
                if other not in index or not count:
                    continue
                for (x, y), n in pairs.items():
                    if x == other and y in index:
                        matrix[index[c]][index[y]] += land * count * p[other] * Fraction(n, ends[other])
    m = [[(Fraction(1) if i == j else Fraction(0)) - matrix[i][j] for j in range(len(classes))]
         for i in range(len(classes))]
    for size in range(1, len(classes) + 1):
        if _determinant([row[:size] for row in m[:size]]) <= 0:
            return False, matrix
    return True, matrix


def _determinant(rows) -> Fraction:
    a = [list(row) for row in rows]
    det = Fraction(1)
    for col in range(len(a)):
        pivot = next((r for r in range(col, len(a)) if a[r][col] != 0), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det = -det
        det *= a[col][col]
        for r in range(col + 1, len(a)):
            factor = a[r][col] / a[col][col]
            for c in range(col, len(a)):
                a[r][c] -= factor * a[col][c]
    return det


def _vacancy_cap() -> "GraphDeclaration":
    return GraphDeclaration.from_mapping(
        {"atoms": [{"key": "h", "element": "H", "hydrogens": 0, "electronic_state": {"radical_electrons": 0}}],
         "bonds": [], "ports": {"cap": {"atom": "h", "bond_order": 1}},
         "electronic_state": {"multiplicity": 1}}, single_connection_port=True)


@dataclass(frozen=True)
class BranchedPolymerSpec:
    """A branched v2 declaration (comment above)."""

    schema_version: int
    label: str
    architecture: str
    residues: Tuple[Tuple[str, ResidueDeclaration], ...]
    links: Tuple[Tuple[str, LinkDeclaration], ...]
    minimum_cycle_bonds: int
    monomer_products: Tuple[Tuple[str, GraphDeclaration], ...]
    small_species: Tuple[Tuple[str, GraphDeclaration, str], ...]
    initial_ensemble: BranchedEnsemble
    initial_mass_kg: float
    density_kg_m3: float
    density_basis: str
    initial_mobile_counts: Tuple[Tuple[str, int], ...]
    locality_policy: Mapping[str, Any]
    family_manifest: Mapping[str, Any]
    library_manifest: Mapping[str, Any]
    transport_model: Mapping[str, Any]
    diagnostics: Mapping[str, Any]

    @property
    def vacancy_cap(self) -> GraphDeclaration:
        """The H that caps a vacant site (as a one-port end group)."""
        return _vacancy_cap()

    @classmethod
    def from_mapping(cls, value: Mapping) -> "BranchedPolymerSpec":
        if not isinstance(value, Mapping):
            raise ValueError("BranchedPolymerSpec input must be a mapping")
        version = _integer(value.get("schema_version"), "schema_version", 1)
        if version != SCHEMA_VERSION_V2:
            raise ValueError(f"BranchedPolymerSpec reads schema_version {SCHEMA_VERSION_V2}, got {version}")
        if value.get("architecture") != "branched":
            raise ValueError(f"a BranchedPolymerSpec declares architecture 'branched', got "
                             f"{value.get('architecture')!r}")
        required = {
            "schema_version", "label", "architecture", "residues", "links", "small_species",
            "initial_ensemble", "initial_mass_kg", "density_kg_m3", "locality_policy",
            "family_manifest", "library_manifest", "transport_model", "diagnostics",
        }
        optional = {"minimum_cycle_bonds", "monomer_products", "density_basis", "initial_mobile_counts"}
        missing, unknown = required - set(value), set(value) - required - optional
        if missing:
            raise ValueError(f"missing BranchedPolymerSpec fields: {sorted(missing)}")
        if unknown:
            raise ValueError(f"unknown BranchedPolymerSpec fields: {sorted(unknown)} (a branched "
                             "declaration has no sequence_support or end_groups: a chain end is a vacant site)")
        label = value["label"]
        if not isinstance(label, str) or not label:
            raise ValueError("spec label must be a non-empty display string")

        raw_residues = value["residues"]
        if not isinstance(raw_residues, Mapping) or not raw_residues:
            raise ValueError("residues must be a non-empty mapping")
        residues = {str(name): ResidueDeclaration.from_mapping(item, str(name), VACANCIES)
                    for name, item in raw_residues.items()}
        raw_links = value["links"]
        if not isinstance(raw_links, Mapping) or not raw_links:
            raise ValueError("links must be a non-empty mapping")
        links = {str(name): LinkDeclaration.from_mapping(item, str(name)) for name, item in raw_links.items()}

        locality = _validated_locality_policy(value["locality_policy"])
        bounds = locality["bounds"]
        classes = set(bounds["attachment_classes"])
        vacancy = {}
        for rname, residue in sorted(residues.items()):
            for sname, site in residue.sites:
                if site.site_class not in classes:
                    raise ValueError(f"site class {site.site_class!r} of {rname}.{sname} is not in "
                                     "locality_policy.bounds.attachment_classes")
                if site.bond_order != 1.0:
                    raise ValueError(f"refused:site_order: site {rname}.{sname} has bond order "
                                     f"{site.bond_order:g}; a branched site is joined (or capped) by one "
                                     "single bond")
                vacancy.setdefault(site.site_class, set()).add(site.vacancy)
        mixed = sorted(c for c, kinds in vacancy.items() if len(kinds) > 1)
        if mixed:
            raise ValueError(f"refused:mixed_vacancy: site classes {mixed} are forbidden-vacant on one "
                             "residue and vacant-capable on another; the random matching cannot keep "
                             "every forbidden site linked")
        if not any("H" in kinds for kinds in vacancy.values()):
            raise ValueError("a branched declaration without a vacant-capable site branches nowhere: "
                             "declare it as linear (architecture 'linear')")
        linked = set()
        for lname, link in sorted(links.items()):
            for site_class in (link.source, link.target):
                if site_class not in classes:
                    raise ValueError(f"link {lname} class {site_class!r} is not in "
                                     "locality_policy.bounds.attachment_classes")
                if site_class not in vacancy:
                    raise ValueError(f"link {lname} class {site_class!r} is no residue's site")
                linked.add(site_class)
            if link.bond_order != 1.0:
                raise ValueError(f"refused:site_order: link {lname} has bond order {link.bond_order:g}")
        unlinked = sorted(set(vacancy) - linked)
        if unlinked:
            raise ValueError(f"site classes {unlinked} take no declared link")
        for lname, link in sorted(links.items()):
            for a, ra in sorted(residues.items()):
                for _x, xs in ra.sites:
                    if xs.site_class != link.source:
                        continue
                    for b, rb in sorted(residues.items()):
                        for _y, ys in rb.sites:
                            if ys.site_class == link.target:
                                _junction_rules(residues, a, b, [(xs.atom, ys.atom, link.bond_order)])
        vocabulary = bounds["topology_metadata_schema"]
        if set(vocabulary) - set(TOPOLOGY_VOCABULARY):
            raise ValueError(f"unknown topology_metadata_schema keys: "
                             f"{sorted(set(vocabulary) - set(TOPOLOGY_VOCABULARY))}")
        for key, names in vocabulary.items():
            if set(_names(names, f"topology_metadata_schema {key}")) - set(TOPOLOGY_VOCABULARY[key]):
                raise ValueError(f"topology_metadata_schema {key} names outside "
                                 f"{list(TOPOLOGY_VOCABULARY[key])}")
        if "branched" not in vocabulary.get("architecture", ()):
            raise ValueError("architecture 'branched' is not in locality_policy.bounds.topology_metadata_schema")
        radius = int(locality["radius"])
        floor = default_minimum_cycle_bonds(radius)
        minimum = value.get("minimum_cycle_bonds", floor)
        if isinstance(minimum, bool) or not isinstance(minimum, int) or minimum < floor:
            raise ValueError(f"refused:minimum_cycle_bonds: {minimum!r} is below 2(r+1)+2 = {floor} at "
                             f"radius {radius}: a shorter cycle fits inside one environment, which the "
                             "environment set (a tree set) does not hold")

        raw_products = value.get("monomer_products", {})
        if not isinstance(raw_products, Mapping) or set(raw_products) - set(residues):
            raise ValueError("monomer_products maps declared residues to physical species")
        products = []
        for name, graph in sorted(raw_products.items()):
            graph = GraphDeclaration.from_mapping(graph)
            if graph.ports:
                raise ValueError("a monomer product is a physical species and cannot contain ports")
            products.append((name, graph))
        small_species = _small_species(value["small_species"])
        mobile = value.get("initial_mobile_counts", {})
        if not isinstance(mobile, Mapping):
            raise ValueError("initial_mobile_counts must be a mapping")
        mobile_counts = tuple(sorted((str(name), _integer(count, "representative mobile count", 0))
                                     for name, count in mobile.items()))
        if not {name for name, _count in mobile_counts} <= {name for name, _g, _d in small_species}:
            raise ValueError("initial_mobile_counts names a species absent from small_species")
        density_basis = value.get("density_basis", "dry_polymer")
        if density_basis not in ("dry_polymer", "whole_initial_mixture"):
            raise ValueError("density_basis must be dry_polymer or whole_initial_mixture")
        if any(count for _name, count in mobile_counts) and density_basis != "whole_initial_mixture":
            raise ValueError("mobile representative counts require whole_initial_mixture density")
        ensemble = BranchedEnsemble.from_mapping(value["initial_ensemble"], tuple(residues), tuple(links))
        if ensemble.generator not in vocabulary.get("network_generators", ()):
            raise ValueError(f"generator {ensemble.generator!r} is not in "
                             "locality_policy.bounds.topology_metadata_schema network_generators")
        policies = []
        for name in ("family_manifest", "library_manifest", "transport_model", "diagnostics"):
            if not isinstance(value[name], Mapping):
                raise ValueError(f"{name} must be a mapping")
            policies.append(_canonical_data(value[name]))
        spec = cls(SCHEMA_VERSION_V2, label, "branched", tuple(sorted(residues.items())),
                   tuple(sorted(links.items())), minimum, tuple(products), small_species, ensemble,
                   _finite(value["initial_mass_kg"], "initial_mass_kg"),
                   _finite(value["density_kg_m3"], "density_kg_m3"), density_basis, mobile_counts,
                   locality, *policies)
        if spec.initial_mass_kg < 0.0:
            raise ValueError("initial_mass_kg cannot be negative")
        if spec.density_kg_m3 <= 0.0:
            raise ValueError("density_kg_m3 must be positive")
        spec._validate_counts()
        spec._validate_coverage_and_states()
        return spec

    def _validate_counts(self) -> None:
        residue_counts, link_counts = dict(self.initial_ensemble.residue_counts), \
            dict(self.initial_ensemble.link_counts)
        sites, ends, _pairs = _class_tables(self.residues, self.links, residue_counts, link_counts)
        forbidden = {site.site_class for _n, residue in self.residues for _s, site in residue.sites
                     if site.vacancy == "forbidden"}
        for c in sorted(sites):
            if ends.get(c, 0) > sites[c]:
                raise ValueError(f"refused:link_counts: {ends[c]} link ends of class {c!r} but only "
                                 f"{sites[c]} such sites in the declared residues")
            if c in forbidden and ends.get(c, 0) != sites[c]:
                raise ValueError(f"refused:link_counts: the {sites[c]} forbidden-vacant sites of class "
                                 f"{c!r} need exactly that many link ends, the counts give {ends.get(c, 0)}")
        below, _matrix = gel_criterion(self.residues, self.links, residue_counts, link_counts)
        if not below:
            raise ValueError("refused:gel_point: the declared counts are at or above the gel point (the "
                             "half-edge branching process of the random matching is critical or "
                             "supercritical): a network is S7")

    def _validate_coverage_and_states(self) -> None:
        cap = self.vacancy_cap
        for _name, residue in self.residues:
            _v1_locality(self.locality_policy, residue.graph, cap, cap)
        first = self.residues[0][1].graph
        for _name, product in self.monomer_products:
            _v1_locality(self.locality_policy, first, cap, cap, monomer=product)
        _v1_locality(self.locality_policy, first, cap, cap, small=self.small_species)
        for name, graph in list(self.monomer_products) + [(n, g) for n, g, _d in self.small_species]:
            molecule = graph.to_legacy_molecule({})
            molecule.multiplicity = graph.multiplicity
            _check_electronic_state(molecule, f"species {name}")

    @property
    def normalized_json(self) -> bytes:
        return json.dumps(
            self.as_mapping(), sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")

    def residence(self, name: str) -> str:
        """Every species is melt-resident (``DEFAULT_RESIDENCE``)."""
        return DEFAULT_RESIDENCE

    def as_mapping(self) -> Dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "label": self.label,
            "architecture": self.architecture,
            "residues": {name: residue.as_mapping() for name, residue in self.residues},
            "links": {name: link.as_mapping() for name, link in self.links},
            "minimum_cycle_bonds": self.minimum_cycle_bonds,
            "monomer_products": {name: graph.as_mapping() for name, graph in self.monomer_products},
            "small_species": {
                name: {"graph": graph.as_mapping(), "disposition": disposition}
                for name, graph, disposition in self.small_species
            },
            "initial_ensemble": self.initial_ensemble.as_mapping(),
            "initial_mass_kg": self.initial_mass_kg,
            "density_kg_m3": self.density_kg_m3,
            "density_basis": self.density_basis,
            "initial_mobile_counts": {name: count for name, count in self.initial_mobile_counts},
            "locality_policy": self.locality_policy,
            "family_manifest": self.family_manifest,
            "library_manifest": self.library_manifest,
            "transport_model": self.transport_model,
            "diagnostics": self.diagnostics,
        }


class BranchedMomRefusal(ValueError):
    """The method of moments was given a branched declaration (deferred: a
    finite branched moment seed needs branch-point moments)."""

    status = "refused:branched_mom"

    def __init__(self, consumer: str):
        super().__init__(f"{self.status}: {consumer}: the method of moments reads linear declarations; "
                         "branched chains need branch-point moments (deferred)")
        self.consumer = consumer


def is_branched(value: Any) -> bool:
    """Whether ``value`` is a branched v2 declaration (an object, a mapping or its JSON)."""
    if isinstance(value, BranchedPolymerSpec):
        return True
    if isinstance(value, (str, bytes, bytearray)):
        try:
            value = json.loads(value)
        except ValueError:
            return False
    return isinstance(value, Mapping) and is_schema_v2(value) and value.get("architecture") == "branched"
