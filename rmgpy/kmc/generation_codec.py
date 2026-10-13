"""Lossless, bounded JSON codec for the output of one public family generation.

``encode_generated_reactions(reactions) -> bytes`` serializes the object graph
that ``generate_reactions_from_families`` returns; ``decode_generated_reactions``
rebuilds it.  The format is explicit and closed:

* Reaction, Species, Molecule, Atom, Bond, kinetics and thermo objects live in
  reference tables.  A shared instance is one table row referenced many times;
  two isomorphic but distinct instances are two rows.  ``reverse`` links
  (including cycles) and ``pairs`` endpoints are table references, never rebuilt
  by graph matching.
* Encoding raises ``EncodeError`` for any populated field or type the schema
  does not carry (the caller treats that as a cache miss).  Nothing is encoded
  with ``repr``, ``eval``, pickle or import-by-name.
* Decoding never estimates anything: atom types, connectivity values, element
  singletons and quantities are restored, not recomputed.  Every rejection is a
  ``DecodeError``, raised before any RMG object is constructed whenever the
  rejection is structural (bytes, nesting, tokens, duplicate keys, non-finite
  numbers, types, counts, references, bonds).  The only checks that need an RMG
  object, unit-string resolution and element lookup, are done in the same
  validation pass on read-only tables.

Not carried (typed rejection on encode): Species conformer, transport data,
energy transfer model, solvation data, props other than ``kmc_frozen_thermo``;
Molecule props, metal, facet, Fragment/Polymer subclasses; Atom coords; any
Reaction specific collider, transition state, network kinetics, polymer
flux-handshake stamps, non-empty rate cache or ``labeled_atoms``; kinetics
other than Arrhenius, ArrheniusEP, ArrheniusBM (and ``solute``); thermo other
than ThermoData.  Python list aliasing (a reaction's ``products`` list being the
same list object as its reverse's ``reactants``) is not preserved; element
sharing is.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from typing import Any, Callable, Optional, Sequence

SCHEMA = "kmc.generated-reactions/1"


class CodecError(Exception):
    """Base class of every codec failure."""


class EncodeError(CodecError):
    """The object graph holds a populated field or type the schema cannot carry."""


class DecodeError(CodecError):
    """The payload is not a well-formed, in-bounds encoding."""


@dataclass(frozen=True)
class DecodeLimits:
    """Hard bounds enforced before and during decoding."""

    max_bytes: int = 256 * 1024 * 1024
    max_depth: int = 24
    max_tokens: int = 100_000_000
    max_string_length: int = 8192
    max_reactions: int = 1_000_000
    max_species: int = 2_000_000
    max_molecules: int = 2_000_000
    max_atoms: int = 20_000_000
    max_bonds: int = 20_000_000
    max_kinetics: int = 2_000_000
    max_thermo: int = 2_000_000
    max_references: int = 200_000_000


DEFAULT_LIMITS = DecodeLimits()

_INT_LIMIT_DIGITS = 18
_FLOAT_LIMIT_CHARS = 64

# ---------------------------------------------------------------------------
# Class and default tables
# ---------------------------------------------------------------------------

_REACTION_MUST_BE_DEFAULT = (
    ("specific_collider", None),
    ("transition_state", None),
    ("network_kinetics", None),
    ("polymer_refused", False),
    ("polymer_refused_accumulating", False),
    ("polymer_refused_reason", None),
    ("polymer_flux_archetype", 0),
    ("polymer_chip_units", 0),
    ("polymer_eject_units", 0.0),
    ("polymer_conduit_params", None),
    ("polymer_conduit_dst_pool", None),
    ("polymer_conduit_admission_pending", False),
    ("polymer_conduit_admission_readjudicated", False),
)
_SPECIES_MUST_BE_NONE = (
    "conformer",
    "transport_data",
    "energy_transfer_model",
    "liquid_volumetric_mass_transfer_coefficient_data",
    "henry_law_constant_data",
)
FROZEN_THERMO_PROPERTY = "kmc_frozen_thermo"
_TEMPLATE_ATTRIBUTES = {"family", "template", "estimator", "reverse", "labeled_atoms"}

_SCALAR_QUANTITY_KEYS = {"units", "uncertainty_si", "uncertainty_type", "value_si"}
_ARRAY_QUANTITY_KEYS = _SCALAR_QUANTITY_KEYS
_UNCERTAINTY_TYPES = {"+|-", "*|/"}
_KINETICS_FIELDS = {
    "Arrhenius": ("A", "n", "Ea", "T0"),
    "ArrheniusEP": ("A", "n", "alpha", "E0"),
    "ArrheniusBM": ("A", "n", "w0", "E0"),
}
_KINETICS_OPTIONAL = ("Tmin", "Tmax", "Pmin", "Pmax")
_THERMO_REQUIRED = ("H298", "S298")
_THERMO_ARRAYS = ("Tdata", "Cpdata")
_THERMO_OPTIONAL = ("Cp0", "CpInf", "Tmin", "Tmax", "E0")
_PROP_SCALARS = (bool, int, float, str, type(None))

_C_INT = (-(2**31), 2**31 - 1)
_C_SHORT = (-(2**15), 2**15 - 1)


def _modules():
    """Import RMG classes lazily so importing this module stays cheap."""
    from rmgpy.data.kinetics.family import TemplateReaction
    from rmgpy.kinetics.arrhenius import Arrhenius, ArrheniusBM, ArrheniusEP
    from rmgpy.kinetics.uncertainties import RateUncertainty
    from rmgpy.molecule.atomtype import ATOMTYPES
    from rmgpy.molecule.element import get_element
    from rmgpy.molecule.molecule import Atom, Bond, Molecule
    from rmgpy.quantity import ArrayQuantity, ScalarQuantity
    from rmgpy.reaction import Reaction
    from rmgpy.species import Species
    from rmgpy.thermo.thermodata import ThermoData

    return {
        "Arrhenius": Arrhenius,
        "ArrheniusBM": ArrheniusBM,
        "ArrheniusEP": ArrheniusEP,
        "ArrayQuantity": ArrayQuantity,
        "ATOMTYPES": ATOMTYPES,
        "Atom": Atom,
        "Bond": Bond,
        "Molecule": Molecule,
        "RateUncertainty": RateUncertainty,
        "Reaction": Reaction,
        "ScalarQuantity": ScalarQuantity,
        "Species": Species,
        "TemplateReaction": TemplateReaction,
        "ThermoData": ThermoData,
        "get_element": get_element,
    }


# ---------------------------------------------------------------------------
# Encoding
# ---------------------------------------------------------------------------


class _Table:
    """Identity-keyed reference table with deterministic first-visit order."""

    def __init__(self) -> None:
        self.index: dict[int, int] = {}
        self.keep: list[Any] = []  # pin objects so ids stay unique while encoding
        self.rows: list[Any] = []

    def reserve(self, obj: Any) -> tuple[int, bool]:
        key = id(obj)
        if key in self.index:
            return self.index[key], False
        position = len(self.rows)
        self.index[key] = position
        self.keep.append(obj)
        self.rows.append(None)
        return position, True


def _is_int(value: Any) -> bool:
    return type(value) is int


def _enc_int(value: Any, name: str, bounds: tuple[int, int] = _C_INT) -> int:
    if type(value) is not int or not bounds[0] <= value <= bounds[1]:
        raise EncodeError(f"{name} must be an int in range, got {value!r}")
    return value


def _enc_float(value: Any, name: str) -> float:
    if type(value) is bool or not isinstance(value, (int, float)) and not _is_numpy_real(value):
        raise EncodeError(f"{name} must be a real number, got {type(value).__name__}")
    result = float(value)
    if not math.isfinite(result):
        raise EncodeError(f"{name} must be finite")
    return result


def _is_numpy_real(value: Any) -> bool:
    import numpy as np

    return isinstance(value, np.floating)


def _enc_bool(value: Any, name: str) -> bool:
    if type(value) is not bool:
        raise EncodeError(f"{name} must be a bool, got {value!r}")
    return value


def _enc_str(value: Any, name: str) -> str:
    if type(value) is not str:
        raise EncodeError(f"{name} must be a str, got {type(value).__name__}")
    return value


def _enc_optional_str(value: Any, name: str) -> Optional[str]:
    return None if value is None else _enc_str(value, name)


class _Encoder:
    def __init__(self) -> None:
        self.m = _modules()
        self.reactions = _Table()
        self.species = _Table()
        self.molecules = _Table()
        self.atoms = _Table()
        self.bonds = _Table()
        self.kinetics = _Table()
        self.thermo = _Table()
        self.pending_atoms: list[tuple[int, Any]] = []

    # -- top level ---------------------------------------------------------

    def encode(self, reactions: Any) -> dict[str, Any]:
        if type(reactions) is not list:
            raise EncodeError("generated reactions must be a list")
        result = [self.reaction(item) for item in reactions]
        while self.pending_atoms:
            position, atom = self.pending_atoms.pop()
            self.atoms.rows[position]["edges"] = [
                [self.atom_ref(neighbour), self.bond_ref(bond)]
                for neighbour, bond in atom.edges.items()
            ]
        return {
            "atoms": self.atoms.rows,
            "bonds": self.bonds.rows,
            "kinetics": self.kinetics.rows,
            "molecules": self.molecules.rows,
            "reactions": self.reactions.rows,
            "result": result,
            "schema": SCHEMA,
            "species": self.species.rows,
            "thermo": self.thermo.rows,
        }

    # -- reactions ---------------------------------------------------------

    def reaction(self, reaction: Any) -> int:
        m = self.m
        cls = type(reaction)
        if cls is m["TemplateReaction"]:
            template = True
        elif cls is m["Reaction"]:
            template = False
        else:
            raise EncodeError(f"unsupported reaction class {cls.__name__}")
        position, new = self.reactions.reserve(reaction)
        if not new:
            return position
        for name, default in _REACTION_MUST_BE_DEFAULT:
            value = getattr(reaction, name, default)
            if value != default or type(value) is not type(default):
                raise EncodeError(f"reaction.{name} is populated: {value!r}")
        if getattr(reaction, "k_effective_cache", {}) != {}:
            raise EncodeError("reaction.k_effective_cache is populated")
        extra = set(getattr(reaction, "__dict__", {}))
        if template:
            extra -= _TEMPLATE_ATTRIBUTES
        if extra:
            raise EncodeError(f"unsupported reaction attributes {sorted(extra)}")
        record: dict[str, Any] = {
            "allow_max_rate_violation": _enc_bool(
                reaction.allow_max_rate_violation, "allow_max_rate_violation"
            ),
            "allow_pdep_route": _enc_bool(reaction.allow_pdep_route, "allow_pdep_route"),
            "class": cls.__name__,
            "comment": _enc_str(reaction.comment, "comment"),
            "degeneracy": _enc_float(reaction.degeneracy, "degeneracy"),
            "duplicate": _enc_bool(reaction.duplicate, "duplicate"),
            "electrons": _enc_int(reaction.electrons, "electrons"),
            "elementary_high_p": _enc_bool(reaction.elementary_high_p, "elementary_high_p"),
            "index": _enc_int(reaction.index, "index"),
            "is_end_group_reaction": _enc_bool(
                reaction.is_end_group_reaction, "is_end_group_reaction"
            ),
            "is_forward": _enc_bool(reaction.is_forward, "is_forward"),
            "kinetics": self.kinetics_ref(reaction.kinetics),
            "label": _enc_str(reaction.label, "label"),
            "pairs": self.pairs(reaction.pairs),
            "products": self.species_list(reaction.products, "products"),
            "protons": _enc_int(reaction._protons, "protons"),
            "rank": self.rank(reaction.rank),
            "reactants": self.species_list(reaction.reactants, "reactants"),
            "reversible": _enc_bool(reaction.reversible, "reversible"),
        }
        if template:
            labeled = reaction.__dict__.get("labeled_atoms", None)
            has_labeled = "labeled_atoms" in reaction.__dict__
            if has_labeled and labeled != {"reactants": {}, "products": {}}:
                raise EncodeError("reaction.labeled_atoms is populated")
            template_value = reaction.template
            if template_value is not None:
                if type(template_value) is not list:
                    raise EncodeError("reaction.template must be a list of str")
                template_value = [_enc_str(item, "template item") for item in template_value]
            reverse = reaction.reverse
            record.update(
                {
                    "estimator": _enc_optional_str(reaction.estimator, "estimator"),
                    "family": _enc_optional_str(reaction.family, "family"),
                    "has_labeled_atoms": has_labeled,
                    "reverse": None if reverse is None else self.reaction(reverse),
                    "template": template_value,
                }
            )
        self.reactions.rows[position] = record
        return position

    @staticmethod
    def rank(value: Any) -> Optional[int]:
        if value is None:
            return None
        return _enc_int(value, "rank")

    def species_list(self, species: Any, name: str) -> list[int]:
        if type(species) is not list:
            raise EncodeError(f"reaction.{name} must be a list")
        return [self.species_ref(item) for item in species]

    def pairs(self, pairs: Any) -> Optional[list[dict[str, Any]]]:
        if pairs is None:
            return None
        if type(pairs) is not list:
            raise EncodeError("reaction.pairs must be a list or None")
        encoded = []
        for pair in pairs:
            kind = {list: "list", tuple: "tuple"}.get(type(pair))
            if kind is None or len(pair) != 2:
                raise EncodeError("each pair must be a 2-element list or tuple")
            encoded.append({"items": [self.species_ref(item) for item in pair], "kind": kind})
        return encoded

    # -- species / molecules / atoms --------------------------------------

    def species_ref(self, species: Any) -> int:
        m = self.m
        if type(species) is not m["Species"]:
            raise EncodeError(f"unsupported species class {type(species).__name__}")
        position, new = self.species.reserve(species)
        if not new:
            return position
        for name in _SPECIES_MUST_BE_NONE:
            if getattr(species, name, None) is not None:
                raise EncodeError(f"species.{name} is populated")
        props = species.props
        if type(props) is not dict or set(props) - {FROZEN_THERMO_PROPERTY}:
            raise EncodeError(f"unsupported species props {sorted(map(str, props))}")
        if FROZEN_THERMO_PROPERTY in props:
            _enc_bool(props[FROZEN_THERMO_PROPERTY], FROZEN_THERMO_PROPERTY)
        if type(species.molecule) is not list:
            raise EncodeError("species.molecule must be a list")
        weight = species._molecular_weight
        record = {
            "aug_inchi": _enc_optional_str(species.aug_inchi, "aug_inchi"),
            "creation_iteration": _enc_int(species.creation_iteration, "creation_iteration"),
            "explicitly_allowed": _enc_bool(species.explicitly_allowed, "explicitly_allowed"),
            "index": _enc_int(species.index, "index"),
            "is_moment_dummy": _enc_bool(species.is_moment_dummy, "is_moment_dummy"),
            "is_polymer_proxy": _enc_bool(species._is_polymer_proxy, "is_polymer_proxy"),
            "is_solvent": _enc_bool(species.is_solvent, "is_solvent"),
            "label": _enc_str(species.label, "label"),
            "molecular_weight": None if weight is None else self.scalar_quantity(weight, "molecular_weight"),
            "molecules": [self.molecule_ref(item) for item in species.molecule],
            "props": dict(props),
            "reactive": _enc_bool(species.reactive, "reactive"),
            "symmetry_number": _enc_float(species.symmetry_number, "symmetry_number"),
            "thermo": self.thermo_ref(species.thermo),
        }
        self.species.rows[position] = record
        return position

    def molecule_ref(self, molecule: Any) -> int:
        m = self.m
        if type(molecule) is not m["Molecule"]:
            raise EncodeError(f"unsupported molecule class {type(molecule).__name__}")
        position, new = self.molecules.reserve(molecule)
        if not new:
            return position
        if molecule.props or molecule.metal or molecule.facet:
            raise EncodeError("molecule props, metal or facet are populated")
        record = {
            "atoms": [self.atom_ref(atom) for atom in molecule.atoms],
            "is_polymer_proxy": _enc_bool(molecule.is_polymer_proxy, "is_polymer_proxy"),
            "multiplicity": _enc_int(molecule.multiplicity, "multiplicity"),
            "reactive": _enc_bool(molecule.reactive, "reactive"),
            "symmetry_number": _enc_float(molecule.symmetry_number, "symmetry_number"),
        }
        self.molecules.rows[position] = record
        return position

    def atom_ref(self, atom: Any) -> int:
        m = self.m
        if type(atom) is not m["Atom"]:
            raise EncodeError(f"unsupported atom class {type(atom).__name__}")
        position, new = self.atoms.reserve(atom)
        if not new:
            return position
        if atom.coords is not None and getattr(atom.coords, "size", 0):
            raise EncodeError("atom.coords is populated")
        if atom.mapping is not None:
            raise EncodeError("atom.mapping is populated")
        props = atom.props
        if type(props) is not dict:
            raise EncodeError("atom.props must be a dict")
        for key, value in props.items():
            if type(key) is not str or type(value) not in _PROP_SCALARS:
                raise EncodeError(f"unsupported atom prop {key!r}={value!r}")
            if type(value) is float and not math.isfinite(value):
                raise EncodeError(f"non-finite atom prop {key!r}")
        record: dict[str, Any] = {
            "atomtype": None if atom.atomtype is None else _enc_str(atom.atomtype.label, "atomtype"),
            "charge": _enc_int(atom.charge, "charge", _C_SHORT),
            "connectivity1": _enc_int(atom.connectivity1, "connectivity1", _C_SHORT),
            "connectivity2": _enc_int(atom.connectivity2, "connectivity2", _C_SHORT),
            "connectivity3": _enc_int(atom.connectivity3, "connectivity3", _C_SHORT),
            "edges": None,
            "element": _enc_str(atom.element.symbol, "element"),
            "id": _enc_int(atom.id, "id"),
            "ignore": _enc_bool(atom.ignore, "ignore"),
            "isotope": _enc_int(atom.element.isotope, "isotope"),
            "label": _enc_str(atom.label, "label"),
            "lone_pairs": _enc_int(atom.lone_pairs, "lone_pairs", _C_SHORT),
            "morphology": _enc_str(atom.morphology, "morphology"),
            "props": dict(props),
            "radical_electrons": _enc_int(atom.radical_electrons, "radical_electrons", _C_SHORT),
            "site": _enc_str(atom.site, "site"),
            "sorting_label": _enc_int(atom.sorting_label, "sorting_label", _C_SHORT),
            "terminal": _enc_bool(atom.terminal, "terminal"),
        }
        self.atoms.rows[position] = record
        self.pending_atoms.append((position, atom))
        return position

    def bond_ref(self, bond: Any) -> int:
        m = self.m
        if type(bond) is not m["Bond"]:
            raise EncodeError(f"unsupported bond class {type(bond).__name__}")
        position, new = self.bonds.reserve(bond)
        if not new:
            return position
        order = bond.order
        if type(order) not in (int, float) or not math.isfinite(order):
            raise EncodeError(f"unsupported bond order {order!r}")
        self.bonds.rows[position] = {
            "order": order,
            "v1": self.atom_ref(bond.vertex1),
            "v2": self.atom_ref(bond.vertex2),
        }
        return position

    # -- quantities, kinetics, thermo --------------------------------------

    def scalar_quantity(self, quantity: Any, name: str) -> dict[str, Any]:
        if type(quantity) is not self.m["ScalarQuantity"]:
            raise EncodeError(f"{name} must be a ScalarQuantity")
        uncertainty_type = _enc_str(quantity.uncertainty_type, name)
        if uncertainty_type not in _UNCERTAINTY_TYPES:
            raise EncodeError(f"{name} has unsupported uncertainty type")
        return {
            "uncertainty_si": _enc_float(quantity.uncertainty_si, name),
            "uncertainty_type": uncertainty_type,
            "units": _enc_str(quantity.units, name),
            "value_si": _enc_float(quantity.value_si, name),
        }

    def array_quantity(self, quantity: Any, name: str) -> dict[str, Any]:
        if type(quantity) is not self.m["ArrayQuantity"]:
            raise EncodeError(f"{name} must be an ArrayQuantity")
        value, uncertainty = quantity.value_si, quantity.uncertainty_si
        if getattr(value, "ndim", None) != 1 or getattr(uncertainty, "shape", None) != value.shape:
            raise EncodeError(f"{name} must be a one-dimensional array")
        uncertainty_type = _enc_str(quantity.uncertainty_type, name)
        if uncertainty_type not in _UNCERTAINTY_TYPES:
            raise EncodeError(f"{name} has unsupported uncertainty type")
        return {
            "uncertainty_si": [_enc_float(item, name) for item in uncertainty.tolist()],
            "uncertainty_type": uncertainty_type,
            "units": _enc_str(quantity.units, name),
            "value_si": [_enc_float(item, name) for item in value.tolist()],
        }

    def optional_scalar(self, quantity: Any, name: str) -> Optional[dict[str, Any]]:
        return None if quantity is None else self.scalar_quantity(quantity, name)

    def kinetics_ref(self, kinetics: Any) -> Optional[int]:
        if kinetics is None:
            return None
        name = type(kinetics).__name__
        if type(kinetics) is not self.m.get(name) or name not in _KINETICS_FIELDS:
            raise EncodeError(f"unsupported kinetics class {name}")
        position, new = self.kinetics.reserve(kinetics)
        if not new:
            return position
        if kinetics.solute is not None:
            raise EncodeError("kinetics.solute is populated")
        uncertainty = kinetics.uncertainty
        encoded_uncertainty = None
        if uncertainty is not None:
            if type(uncertainty) is not self.m["RateUncertainty"]:
                raise EncodeError("unsupported kinetics uncertainty class")
            correlation = uncertainty.correlation
            encoded_uncertainty = {
                "correlation": _enc_optional_str(correlation, "correlation"),
                "data_mean": _enc_float(uncertainty.data_mean, "data_mean"),
                "N": _enc_int(uncertainty.N, "N"),
                "Tref": _enc_float(uncertainty.Tref, "Tref"),
                "mu": _enc_float(uncertainty.mu, "mu"),
                "var": _enc_float(uncertainty.var, "var"),
            }
        record: dict[str, Any] = {
            "class": name,
            "comment": _enc_str(kinetics.comment, "comment"),
            "uncertainty": encoded_uncertainty,
        }
        for field in _KINETICS_FIELDS[name]:
            value = getattr(kinetics, "_" + field)
            if value is None:
                raise EncodeError(f"kinetics.{field} is missing")
            record[field] = self.scalar_quantity(value, field)
        for field in _KINETICS_OPTIONAL:
            record[field] = self.optional_scalar(getattr(kinetics, "_" + field), field)
        self.kinetics.rows[position] = record
        return position

    def thermo_ref(self, thermo: Any) -> Optional[int]:
        if thermo is None:
            return None
        if type(thermo) is not self.m["ThermoData"]:
            raise EncodeError(f"unsupported thermo class {type(thermo).__name__}")
        position, new = self.thermo.reserve(thermo)
        if not new:
            return position
        record: dict[str, Any] = {
            "comment": _enc_str(thermo.comment, "comment"),
            "label": _enc_str(thermo.label, "label"),
        }
        for field in _THERMO_ARRAYS:
            record[field] = self.array_quantity(getattr(thermo, "_" + field), field)
        for field in _THERMO_REQUIRED:
            record[field] = self.scalar_quantity(getattr(thermo, "_" + field), field)
        for field in _THERMO_OPTIONAL:
            record[field] = self.optional_scalar(getattr(thermo, "_" + field), field)
        self.thermo.rows[position] = record
        return position


def encode_generated_reactions(reactions: Sequence[Any]) -> bytes:
    """Serialize a generator result; raise ``EncodeError`` on anything unsupported."""
    try:
        document = _Encoder().encode(reactions)
        return json.dumps(
            document, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False
        ).encode("ascii")
    except CodecError:
        raise
    except (ValueError, TypeError, AttributeError, OverflowError, RecursionError) as error:
        raise EncodeError(f"cannot encode generated reactions: {error}") from error


# ---------------------------------------------------------------------------
# Decoding: byte preflight, strict parse, whole-structure validation, build
# ---------------------------------------------------------------------------

_TOKEN = re.compile(rb'(?s)"[^"\\]*(?:\\.[^"\\]*)*"|[\[\]{}:,]|"')


def _preflight(data: bytes, limits: DecodeLimits) -> None:
    """Single pass over the raw bytes: nesting depth, token count, string length."""
    depth = 0
    tokens = 0
    for match in _TOKEN.finditer(data):
        tokens += 1
        if tokens > limits.max_tokens:
            raise DecodeError("token count exceeds the limit")
        token = match.group()
        first = token[:1]
        if first == b'"':
            if len(token) < 2:
                raise DecodeError("unterminated string")
            if len(token) > 6 * limits.max_string_length + 2:
                raise DecodeError("string exceeds the length limit")
        elif first == b"[" or first == b"{":
            depth += 1
            if depth > limits.max_depth:
                raise DecodeError("nesting depth exceeds the limit")
        elif first == b"]" or first == b"}":
            depth -= 1
            if depth < 0:
                raise DecodeError("unbalanced brackets")
    if depth != 0:
        raise DecodeError("unbalanced brackets")


def _pairs_hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise DecodeError(f"duplicate key {key!r}")
        result[key] = value
    return result


def _reject_constant(token: str) -> Any:
    raise DecodeError(f"non-finite constant {token}")


def _parse_float(token: str) -> float:
    if len(token) > _FLOAT_LIMIT_CHARS:
        raise DecodeError("float literal too long")
    value = float(token)
    if not math.isfinite(value):
        raise DecodeError(f"out-of-range float {token}")
    return value


def _parse_int(token: str) -> int:
    if len(token.lstrip("-")) > _INT_LIMIT_DIGITS:
        raise DecodeError("integer literal out of range")
    return int(token)


class _Validator:
    """Check the whole decoded JSON document without constructing RMG objects."""

    def __init__(self, document: Any, limits: DecodeLimits) -> None:
        self.doc = document
        self.limits = limits
        self.references = 0
        self.m = None  # filled lazily for element/atomtype lookups only

    # -- primitives --------------------------------------------------------

    def fail(self, where: str, message: str) -> Any:
        raise DecodeError(f"{where}: {message}")

    def obj(self, value: Any, keys: set[str], where: str) -> dict[str, Any]:
        if type(value) is not dict:
            self.fail(where, "expected an object")
        if set(value) != keys:
            self.fail(where, f"keys {sorted(value)} differ from {sorted(keys)}")
        return value

    def seq(self, value: Any, where: str) -> list:
        if type(value) is not list:
            self.fail(where, "expected an array")
        self.references += len(value)
        if self.references > self.limits.max_references:
            self.fail(where, "reference count exceeds the limit")
        return value

    def integer(self, value: Any, where: str, bounds: tuple[int, int] = _C_INT) -> int:
        if type(value) is not int or not bounds[0] <= value <= bounds[1]:
            self.fail(where, f"expected an in-range integer, got {value!r}")
        return value

    def real(self, value: Any, where: str) -> float:
        if type(value) is not float or not math.isfinite(value):
            self.fail(where, f"expected a finite float, got {value!r}")
        return value

    def boolean(self, value: Any, where: str) -> bool:
        if type(value) is not bool:
            self.fail(where, f"expected a bool, got {value!r}")
        return value

    def string(self, value: Any, where: str, optional: bool = False) -> Optional[str]:
        if value is None and optional:
            return None
        if type(value) is not str:
            self.fail(where, "expected a string")
        if len(value) > self.limits.max_string_length:
            self.fail(where, "string exceeds the length limit")
        return value

    def ref(self, value: Any, size: int, where: str, optional: bool = False) -> Optional[int]:
        if value is None and optional:
            return None
        if type(value) is not int or not 0 <= value < size:
            self.fail(where, f"invalid reference {value!r}")
        return value

    # -- quantities --------------------------------------------------------

    def units(self, value: Any, where: str) -> str:
        text = self.string(value, where)
        if not re.fullmatch(r"[A-Za-z0-9_/*^()\-\. %]{0,128}", text):
            self.fail(where, f"unsupported unit string {text!r}")
        try:
            import quantities as pq

            pq.Quantity(1.0, text)
        except Exception:
            self.fail(where, f"unresolvable unit string {text!r}")
        return text

    def scalar_quantity(self, value: Any, where: str) -> None:
        item = self.obj(value, _SCALAR_QUANTITY_KEYS, where)
        self.units(item["units"], where)
        self.real(item["value_si"], where)
        self.real(item["uncertainty_si"], where)
        if item["uncertainty_type"] not in _UNCERTAINTY_TYPES:
            self.fail(where, "unsupported uncertainty type")

    def array_quantity(self, value: Any, where: str) -> None:
        item = self.obj(value, _ARRAY_QUANTITY_KEYS, where)
        self.units(item["units"], where)
        values = self.seq(item["value_si"], where)
        uncertainties = self.seq(item["uncertainty_si"], where)
        if not values or len(values) != len(uncertainties):
            self.fail(where, "array and uncertainty lengths differ or are empty")
        for entry in (*values, *uncertainties):
            self.real(entry, where)
        if item["uncertainty_type"] not in _UNCERTAINTY_TYPES:
            self.fail(where, "unsupported uncertainty type")

    def optional_scalar(self, value: Any, where: str) -> None:
        if value is not None:
            self.scalar_quantity(value, where)

    # -- document ----------------------------------------------------------

    def validate(self) -> None:
        limits = self.limits
        doc = self.obj(
            self.doc,
            {"atoms", "bonds", "kinetics", "molecules", "reactions", "result", "schema", "species", "thermo"},
            "document",
        )
        if doc["schema"] != SCHEMA:
            self.fail("document", f"unsupported schema {doc['schema']!r}")
        tables = {}
        for name, bound in (
            ("reactions", limits.max_reactions),
            ("species", limits.max_species),
            ("molecules", limits.max_molecules),
            ("atoms", limits.max_atoms),
            ("bonds", limits.max_bonds),
            ("kinetics", limits.max_kinetics),
            ("thermo", limits.max_thermo),
        ):
            rows = doc[name]
            if type(rows) is not list:
                self.fail(name, "expected an array")
            if len(rows) > bound:
                self.fail(name, f"{len(rows)} rows exceed the limit {bound}")
            tables[name] = rows
        self.sizes = {name: len(rows) for name, rows in tables.items()}
        for index in self.seq(doc["result"], "result"):
            self.ref(index, self.sizes["reactions"], "result")
        for row in tables["thermo"]:
            self.thermo(row)
        for row in tables["kinetics"]:
            self.kinetics(row)
        self.graph(tables["atoms"], tables["bonds"], tables["molecules"])
        for row in tables["species"]:
            self.species(row)
        for row in tables["reactions"]:
            self.reaction(row)

    def thermo(self, row: Any) -> None:
        keys = {"comment", "label", *_THERMO_ARRAYS, *_THERMO_REQUIRED, *_THERMO_OPTIONAL}
        item = self.obj(row, keys, "thermo")
        self.string(item["comment"], "thermo.comment")
        self.string(item["label"], "thermo.label")
        for field in _THERMO_ARRAYS:
            self.array_quantity(item[field], f"thermo.{field}")
        for field in _THERMO_REQUIRED:
            self.scalar_quantity(item[field], f"thermo.{field}")
        for field in _THERMO_OPTIONAL:
            self.optional_scalar(item[field], f"thermo.{field}")

    def kinetics(self, row: Any) -> None:
        if type(row) is not dict or row.get("class") not in _KINETICS_FIELDS:
            self.fail("kinetics", "unsupported kinetics class")
        name = row["class"]
        keys = {"class", "comment", "uncertainty", *_KINETICS_FIELDS[name], *_KINETICS_OPTIONAL}
        item = self.obj(row, keys, "kinetics")
        self.string(item["comment"], "kinetics.comment")
        for field in _KINETICS_FIELDS[name]:
            self.scalar_quantity(item[field], f"kinetics.{field}")
        for field in _KINETICS_OPTIONAL:
            self.optional_scalar(item[field], f"kinetics.{field}")
        if item["uncertainty"] is not None:
            u = self.obj(
                item["uncertainty"],
                {"correlation", "data_mean", "N", "Tref", "mu", "var"},
                "kinetics.uncertainty",
            )
            self.string(u["correlation"], "kinetics.uncertainty.correlation", optional=True)
            for field in ("data_mean", "Tref", "mu", "var"):
                self.real(u[field], f"kinetics.uncertainty.{field}")
            self.integer(u["N"], "kinetics.uncertainty.N")

    def graph(self, atoms: list, bonds: list, molecules: list) -> None:
        atom_total, bond_total = len(atoms), len(bonds)
        endpoints: dict[int, tuple[int, int]] = {}
        seen_pairs: set[frozenset] = set()
        for index, row in enumerate(bonds):
            item = self.obj(row, {"order", "v1", "v2"}, "bond")
            order = item["order"]
            if type(order) not in (int, float) or not 0 < order <= 4:
                self.fail("bond.order", f"invalid bond order {order!r}")
            v1 = self.ref(item["v1"], atom_total, "bond.v1")
            v2 = self.ref(item["v2"], atom_total, "bond.v2")
            if v1 == v2:
                self.fail("bond", f"self bond on atom {v1}")
            pair = frozenset((v1, v2))
            if pair in seen_pairs:
                self.fail("bond", f"duplicate bond between atoms {v1} and {v2}")
            seen_pairs.add(pair)
            endpoints[index] = (v1, v2)
        referenced_by: dict[int, set[int]] = {}
        neighbours: list[dict[int, int]] = []
        atom_keys = {
            "atomtype", "charge", "connectivity1", "connectivity2", "connectivity3", "edges",
            "element", "id", "ignore", "isotope", "label", "lone_pairs", "morphology", "props",
            "radical_electrons", "site", "sorting_label", "terminal",
        }
        from rmgpy.molecule.atomtype import ATOMTYPES
        from rmgpy.molecule.element import get_element

        for index, row in enumerate(atoms):
            item = self.obj(row, atom_keys, "atom")
            symbol = self.string(item["element"], "atom.element")
            isotope = self.integer(item["isotope"], "atom.isotope")
            try:
                get_element(symbol, isotope)
            except Exception:
                self.fail("atom", f"unknown element {symbol!r} isotope {isotope}")
            label = self.string(item["atomtype"], "atom.atomtype", optional=True)
            if label is not None and label not in ATOMTYPES:
                self.fail("atom", f"unknown atom type {label!r}")
            for field in ("charge", "connectivity1", "connectivity2", "connectivity3",
                          "lone_pairs", "radical_electrons", "sorting_label"):
                self.integer(item[field], f"atom.{field}", _C_SHORT)
            self.integer(item["id"], "atom.id")
            for field in ("ignore", "terminal"):
                self.boolean(item[field], f"atom.{field}")
            for field in ("label", "morphology", "site"):
                self.string(item[field], f"atom.{field}")
            props = item["props"]
            if type(props) is not dict or len(props) > 64:
                self.fail("atom.props", "expected a small object")
            for key, value in props.items():
                self.string(key, "atom.props key")
                if type(value) not in _PROP_SCALARS:
                    self.fail("atom.props", f"unsupported value for {key!r}")
                if type(value) is str:
                    self.string(value, "atom.props value")
            mine: dict[int, int] = {}
            for edge in self.seq(item["edges"], "atom.edges"):
                if type(edge) is not list or len(edge) != 2:
                    self.fail("atom.edges", "expected [neighbour, bond] pairs")
                neighbour = self.ref(edge[0], atom_total, "atom.edges neighbour")
                bond = self.ref(edge[1], bond_total, "atom.edges bond")
                if neighbour == index:
                    self.fail("atom.edges", f"atom {index} bonds to itself")
                if neighbour in mine:
                    self.fail("atom.edges", f"duplicate neighbour {neighbour} on atom {index}")
                if set(endpoints[bond]) != {index, neighbour}:
                    self.fail("atom.edges", f"bond {bond} does not join atoms {index} and {neighbour}")
                mine[neighbour] = bond
                referenced_by.setdefault(bond, set()).add(index)
            neighbours.append(mine)
        for bond, (v1, v2) in endpoints.items():
            if referenced_by.get(bond, set()) != {v1, v2}:
                self.fail("bond", f"bond {bond} is not listed symmetrically by both atoms")
        for index, mine in enumerate(neighbours):
            for neighbour, bond in mine.items():
                if neighbours[neighbour].get(index) != bond:
                    self.fail("atom.edges", f"asymmetric bond between atoms {index} and {neighbour}")
        for row in molecules:
            item = self.obj(
                row,
                {"atoms", "is_polymer_proxy", "multiplicity", "reactive", "symmetry_number"},
                "molecule",
            )
            members = [self.ref(a, atom_total, "molecule.atoms") for a in self.seq(item["atoms"], "molecule.atoms")]
            if not members or len(set(members)) != len(members):
                self.fail("molecule", "atoms must be non-empty and unique")
            member_set = set(members)
            for atom in members:
                if not set(neighbours[atom]) <= member_set:
                    self.fail("molecule", f"atom {atom} bonds outside its molecule")
            self.boolean(item["is_polymer_proxy"], "molecule.is_polymer_proxy")
            self.integer(item["multiplicity"], "molecule.multiplicity")
            self.boolean(item["reactive"], "molecule.reactive")
            self.real(item["symmetry_number"], "molecule.symmetry_number")

    def species(self, row: Any) -> None:
        item = self.obj(
            row,
            {
                "aug_inchi", "creation_iteration", "explicitly_allowed", "index", "is_moment_dummy",
                "is_polymer_proxy", "is_solvent", "label", "molecular_weight", "molecules", "props",
                "reactive", "symmetry_number", "thermo",
            },
            "species",
        )
        self.string(item["aug_inchi"], "species.aug_inchi", optional=True)
        self.integer(item["creation_iteration"], "species.creation_iteration")
        self.integer(item["index"], "species.index")
        for field in ("explicitly_allowed", "is_moment_dummy", "is_polymer_proxy", "is_solvent", "reactive"):
            self.boolean(item[field], f"species.{field}")
        self.string(item["label"], "species.label")
        self.real(item["symmetry_number"], "species.symmetry_number")
        self.optional_scalar(item["molecular_weight"], "species.molecular_weight")
        for index in self.seq(item["molecules"], "species.molecules"):
            self.ref(index, self.sizes["molecules"], "species.molecules")
        props = item["props"]
        if type(props) is not dict or set(props) - {FROZEN_THERMO_PROPERTY}:
            self.fail("species.props", "unsupported properties")
        for value in props.values():
            self.boolean(value, "species.props")
        self.ref(item["thermo"], self.sizes["thermo"], "species.thermo", optional=True)

    def reaction(self, row: Any) -> None:
        base = {
            "allow_max_rate_violation", "allow_pdep_route", "class", "comment", "degeneracy", "duplicate",
            "electrons", "elementary_high_p", "index", "is_end_group_reaction", "is_forward", "kinetics",
            "label", "pairs", "products", "protons", "rank", "reactants", "reversible",
        }
        if type(row) is not dict or row.get("class") not in ("Reaction", "TemplateReaction"):
            self.fail("reaction", "unsupported reaction class")
        template = row["class"] == "TemplateReaction"
        keys = base | ({"estimator", "family", "has_labeled_atoms", "reverse", "template"} if template else set())
        item = self.obj(row, keys, "reaction")
        for field in ("allow_max_rate_violation", "allow_pdep_route", "duplicate", "elementary_high_p",
                      "is_end_group_reaction", "is_forward", "reversible"):
            self.boolean(item[field], f"reaction.{field}")
        self.string(item["comment"], "reaction.comment")
        self.string(item["label"], "reaction.label")
        self.real(item["degeneracy"], "reaction.degeneracy")
        for field in ("electrons", "index", "protons"):
            self.integer(item[field], f"reaction.{field}")
        if item["rank"] is not None:
            self.integer(item["rank"], "reaction.rank")
        self.ref(item["kinetics"], self.sizes["kinetics"], "reaction.kinetics", optional=True)
        for field in ("reactants", "products"):
            for index in self.seq(item[field], f"reaction.{field}"):
                self.ref(index, self.sizes["species"], f"reaction.{field}")
        if item["pairs"] is not None:
            for pair in self.seq(item["pairs"], "reaction.pairs"):
                pair = self.obj(pair, {"items", "kind"}, "reaction.pair")
                if pair["kind"] not in ("list", "tuple"):
                    self.fail("reaction.pair", "unsupported pair container")
                items = self.seq(pair["items"], "reaction.pair")
                if len(items) != 2:
                    self.fail("reaction.pair", "expected two endpoints")
                for index in items:
                    self.ref(index, self.sizes["species"], "reaction.pair")
        if template:
            self.string(item["estimator"], "reaction.estimator", optional=True)
            self.string(item["family"], "reaction.family", optional=True)
            self.boolean(item["has_labeled_atoms"], "reaction.has_labeled_atoms")
            self.ref(item["reverse"], self.sizes["reactions"], "reaction.reverse", optional=True)
            if item["template"] is not None:
                for entry in self.seq(item["template"], "reaction.template"):
                    self.string(entry, "reaction.template")


class _Builder:
    """Construct the RMG object graph from a validated document."""

    def __init__(self, document: dict[str, Any]) -> None:
        self.doc = document
        self.m = _modules()

    def scalar(self, item: Optional[dict[str, Any]]) -> Any:
        if item is None:
            return None
        quantity = self.m["ScalarQuantity"](0.0, item["units"], 0.0, item["uncertainty_type"])
        quantity.value_si = item["value_si"]
        quantity.uncertainty_si = item["uncertainty_si"]
        return quantity

    def array(self, item: dict[str, Any]) -> Any:
        import numpy as np

        value = np.array(item["value_si"], dtype=float)
        quantity = self.m["ArrayQuantity"](np.zeros_like(value), item["units"], None, item["uncertainty_type"])
        quantity.value_si = value
        quantity.uncertainty_si = np.array(item["uncertainty_si"], dtype=float)
        return quantity

    def build(self) -> list:
        m, doc = self.m, self.doc
        elements = m["get_element"]
        atoms = []
        for row in doc["atoms"]:
            atom = m["Atom"](
                element=elements(row["element"], row["isotope"]),
                radical_electrons=row["radical_electrons"],
                charge=row["charge"],
                label=row["label"],
                lone_pairs=row["lone_pairs"],
                site=row["site"],
                morphology=row["morphology"],
                id=row["id"],
                props=dict(row["props"]),
            )
            atom.atomtype = None if row["atomtype"] is None else m["ATOMTYPES"][row["atomtype"]]
            atom.connectivity1 = row["connectivity1"]
            atom.connectivity2 = row["connectivity2"]
            atom.connectivity3 = row["connectivity3"]
            atom.sorting_label = row["sorting_label"]
            atom.terminal = row["terminal"]
            atom.ignore = row["ignore"]
            atoms.append(atom)
        bonds = [
            m["Bond"](atoms[row["v1"]], atoms[row["v2"]], row["order"]) for row in doc["bonds"]
        ]
        for atom, row in zip(atoms, doc["atoms"]):
            atom.edges = {}
            for neighbour, bond in row["edges"]:
                atom.edges[atoms[neighbour]] = bonds[bond]
        molecules = []
        for row in doc["molecules"]:
            molecule = m["Molecule"](
                [atoms[i] for i in row["atoms"]],
                symmetry=row["symmetry_number"],
                multiplicity=row["multiplicity"],
                reactive=row["reactive"],
                is_polymer_proxy=row["is_polymer_proxy"],
            )
            molecules.append(molecule)
        thermo = []
        for row in doc["thermo"]:
            optional = {name: self.scalar(row[name]) for name in _THERMO_OPTIONAL}
            item = m["ThermoData"](
                Tdata=self.array(row["Tdata"]),
                Cpdata=self.array(row["Cpdata"]),
                H298=self.scalar(row["H298"]),
                S298=self.scalar(row["S298"]),
                label=row["label"],
                comment=row["comment"],
                **optional,
            )
            # The setters convert through the display value; restore the SI bits.
            item._Tdata = self.array(row["Tdata"])
            item._Cpdata = self.array(row["Cpdata"])
            item._H298 = self.scalar(row["H298"])
            item._S298 = self.scalar(row["S298"])
            for name in _THERMO_OPTIONAL:
                setattr(item, "_" + name, self.scalar(row[name]))
            thermo.append(item)
        kinetics = [self.kinetic(row) for row in doc["kinetics"]]
        species = []
        for row in doc["species"]:
            item = m["Species"](index=row["index"], label=row["label"], reactive=row["reactive"])
            item.molecule = [molecules[i] for i in row["molecules"]]
            item.thermo = None if row["thermo"] is None else thermo[row["thermo"]]
            item._molecular_weight = self.scalar(row["molecular_weight"])
            item.props = dict(row["props"])
            item.aug_inchi = row["aug_inchi"]
            item.symmetry_number = row["symmetry_number"]
            item.is_solvent = row["is_solvent"]
            item.is_moment_dummy = row["is_moment_dummy"]
            item._is_polymer_proxy = row["is_polymer_proxy"]
            item.creation_iteration = row["creation_iteration"]
            item.explicitly_allowed = row["explicitly_allowed"]
            species.append(item)
        reactions = []
        for row in doc["reactions"]:
            if row["class"] == "TemplateReaction":
                reaction = m["TemplateReaction"]()
                if not row["has_labeled_atoms"]:
                    del reaction.labeled_atoms
                reaction.family = row["family"]
                reaction.estimator = row["estimator"]
                reaction.template = None if row["template"] is None else list(row["template"])
            else:
                reaction = m["Reaction"]()
            reactions.append(reaction)
        for reaction, row in zip(reactions, doc["reactions"]):
            reaction.index = row["index"]
            reaction.label = row["label"]
            reaction.reactants = [species[i] for i in row["reactants"]]
            reaction.products = [species[i] for i in row["products"]]
            reaction.kinetics = None if row["kinetics"] is None else kinetics[row["kinetics"]]
            reaction.reversible = row["reversible"]
            reaction.duplicate = row["duplicate"]
            reaction._degeneracy = row["degeneracy"]
            reaction.electrons = row["electrons"]
            reaction._protons = row["protons"]
            reaction.allow_pdep_route = row["allow_pdep_route"]
            reaction.elementary_high_p = row["elementary_high_p"]
            reaction.comment = row["comment"]
            reaction.is_forward = row["is_forward"]
            reaction.allow_max_rate_violation = row["allow_max_rate_violation"]
            reaction.is_end_group_reaction = row["is_end_group_reaction"]
            reaction.rank = row["rank"]
            if row["pairs"] is None:
                reaction.pairs = None
            else:
                reaction.pairs = [
                    (tuple if pair["kind"] == "tuple" else list)(species[i] for i in pair["items"])
                    for pair in row["pairs"]
                ]
            if row["class"] == "TemplateReaction":
                reaction.reverse = None if row["reverse"] is None else reactions[row["reverse"]]
        return [reactions[i] for i in doc["result"]]

    def kinetic(self, row: dict[str, Any]) -> Any:
        m = self.m
        name = row["class"]
        values = {field: self.scalar(row[field]) for field in _KINETICS_FIELDS[name]}
        values.update({field: self.scalar(row[field]) for field in _KINETICS_OPTIONAL})
        uncertainty = None
        if row["uncertainty"] is not None:
            u = row["uncertainty"]
            uncertainty = m["RateUncertainty"](
                mu=u["mu"], var=u["var"], Tref=u["Tref"], N=u["N"], correlation=u["correlation"]
            )
            uncertainty.data_mean = u["data_mean"]
        kinetics = m[name](uncertainty=uncertainty, comment=row["comment"], **values)
        # Quantities were restored from SI values; assign them without re-conversion.
        for field in _KINETICS_FIELDS[name]:
            setattr(kinetics, "_" + field, values[field])
        for field in _KINETICS_OPTIONAL:
            setattr(kinetics, "_" + field, values[field])
        return kinetics


def decode_generated_reactions(payload: bytes, *, limits: DecodeLimits) -> list:
    """Rebuild generator output; every rejection is a ``DecodeError``."""
    if not isinstance(limits, DecodeLimits):
        raise TypeError("limits must be a DecodeLimits")
    try:
        if not isinstance(payload, (bytes, bytearray, memoryview)):
            raise DecodeError("payload must be bytes")
        data = bytes(payload)
        if len(data) > limits.max_bytes:
            raise DecodeError("payload exceeds the byte cap")
        _preflight(data, limits)
        try:
            text = data.decode("utf-8")
        except UnicodeDecodeError as error:
            raise DecodeError(f"payload is not UTF-8: {error}") from error
        document = json.loads(
            text,
            object_pairs_hook=_pairs_hook,
            parse_constant=_reject_constant,
            parse_float=_parse_float,
            parse_int=_parse_int,
        )
        _Validator(document, limits).validate()
        return _Builder(document).build()
    except DecodeError:
        raise
    except (ValueError, TypeError, KeyError, IndexError, AttributeError, RecursionError,
            OverflowError, MemoryError) as error:
        raise DecodeError(f"cannot decode generated reactions: {error}") from error


# ---------------------------------------------------------------------------
# Parity hook
# ---------------------------------------------------------------------------


def roundtrip_generation_provider_factory(
    limits: DecodeLimits = DEFAULT_LIMITS,
) -> Callable[[Callable], Callable]:
    """Build a ``generation_provider_factory`` that round-trips every result.

    The returned provider calls the database's own generator and then returns
    ``decode(encode(result))``, so a compile run through it proves the codec is
    transparent to every consumer of generation output.
    """

    def factory(default_generate: Callable) -> Callable:
        def provider(reactants, products=None, only_families=None, resonance=True):
            generated = default_generate(
                reactants, products=products, only_families=only_families, resonance=resonance
            )
            return decode_generated_reactions(encode_generated_reactions(generated), limits=limits)

        return provider

    return factory
