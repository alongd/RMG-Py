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
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

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
        small_species = []
        for name, declaration in sorted(small.items()):
            if not isinstance(name, str) or not name:
                raise ValueError("small-species labels must be non-empty strings")
            if not isinstance(declaration, Mapping) or set(declaration) != {"graph", "disposition"}:
                raise ValueError("each small species requires graph and disposition")
            graph = GraphDeclaration.from_mapping(declaration["graph"])
            if graph.ports:
                raise ValueError("small physical species cannot contain connection ports")
            disposition = declaration["disposition"]
            if not isinstance(disposition, str) or not disposition:
                raise ValueError("small-species disposition must be a non-empty string")
            small_species.append((name, graph, disposition))
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
        repeat = self.repeat_unit.to_legacy_molecule({"head": "*1", "tail": "*2"})
        head_port_name = self.head_end_group.ports[0][0]
        tail_port_name = self.tail_end_group.ports[0][0]
        head = self.head_end_group.to_legacy_molecule({head_port_name: "*1"})
        tail = self.tail_end_group.to_legacy_molecule({tail_port_name: "*2"})
        return repeat, (head, tail)

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
