from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field, fields
from typing import Any, Literal, Optional
from hashlib import sha256


@dataclass(frozen=True)
class EventRecord:
    """
    Draft schema rewrite_record/0.1 for a compiled kMC event.

    Strand-specific fields (cut_offset, feature_ops, inheritance, junction_ops)
    are documented null placeholders; the compiler will define them later.
    """

    schema_version: str = "rewrite_record/0.1"
    event_id: str = ""
    canonical_index: int = -1
    family: str = ""
    template: str = ""
    arity: int = 0
    participant_site_types: list[str] = field(default_factory=list)
    reactant_multiplicities: list[int] = field(default_factory=list)
    raw_path_degeneracy: float = 1.0
    degeneracy: float = 1.0
    ssa_multiplier: float = 1.0
    rate_order: int = 0
    rate_units: str = ""
    k_table: Optional[dict[str, list[float]]] = None
    atom_map: dict[int, int] = field(default_factory=dict)
    bond_ops: list[dict[str, Any]] = field(default_factory=list)
    element_delta: dict[str, int] = field(default_factory=dict)
    radical_delta: int = 0
    status: Literal["enabled", "refused", "irreversible"] = "enabled"
    status_reason: str = ""
    orientation: Literal["as_generated", "reversed"] = "as_generated"
    thermo_provenance: dict[str, Any] = field(default_factory=dict)
    provenance: dict[str, Any] = field(default_factory=dict)

    # F1 fields added by the offline compiler.  They deliberately have defaults
    # so rewrite_record/0.1 files produced by the atom-map fixture remain
    # readable while the interface is still 0.x.
    site_type: str = ""
    rate_source: dict[str, Any] = field(default_factory=dict)
    coproducts: list[dict[str, Any]] = field(default_factory=list)
    reverse_of: Optional[str] = None
    inventory_class: Optional[str] = None
    implicit_h_delta: int = 0
    formula_delta: dict[str, int] = field(default_factory=dict)
    mass_delta: int = 0
    reactant_pair_convention: str = "ordered-unimolecular"
    resonance_derived: bool = False
    reactant_graphs: list[str] = field(default_factory=list)
    product_graphs: list[str] = field(default_factory=list)

    # Strand-specific placeholders (null until compiler defines them)
    cut_offset: Optional[int] = None
    feature_ops: Optional[list[dict[str, Any]]] = None
    inheritance: Optional[dict[str, Any]] = None
    junction_ops: Optional[list[dict[str, Any]]] = None

    _KNOWN_FIELDS: set[str] = field(
        default_factory=lambda: {f.name for f in fields(EventRecord)},
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self):
        if not self.event_id:
            object.__setattr__(self, "event_id", self._compute_event_id())

    def _compute_event_id(self) -> str:
        payload = {
            item.name: _canonical_link_handles(getattr(self, item.name), item.name)
            for item in fields(self)
            if item.name not in {"event_id", "canonical_index", "_KNOWN_FIELDS"}
        }
        canonical = json.dumps(
            payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("utf-8")
        return f"evt_{sha256(canonical).hexdigest()}"

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d.pop("_KNOWN_FIELDS", None)
        return d

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, sort_keys=True)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> EventRecord:
        version = data.get("schema_version")
        if version != "rewrite_record/0.1":
            raise ValueError(f"Unknown schema version: {version!r}")

        unknown = set(data.keys()) - _KNOWN_FIELDS_CLASS
        if unknown:
            raise ValueError(f"Unknown fields in record: {sorted(unknown)}")

        missing = _LEGACY_REQUIRED_FIELDS - set(data.keys())
        if missing:
            raise ValueError(f"Missing required fields: {sorted(missing)}")

        # JSON serializes dict keys as strings; convert atom_map keys back to int
        if "atom_map" in data and isinstance(data["atom_map"], dict):
            data = dict(data)
            data["atom_map"] = {int(k): v for k, v in data["atom_map"].items()}

        return cls(**data)

    @classmethod
    def from_json(cls, s: str) -> EventRecord:
        return cls.from_dict(json.loads(s))

    def validate(self) -> None:
        """Validate internal consistency."""
        if self.status not in ("enabled", "refused", "irreversible"):
            raise ValueError(f"Invalid status: {self.status!r}")
        if self.orientation not in ("as_generated", "reversed"):
            raise ValueError(f"Invalid orientation: {self.orientation!r}")
        if self.k_table is not None:
            if "T" not in self.k_table or "k" not in self.k_table:
                raise ValueError("k_table must have 'T' and 'k' keys")
            if len(self.k_table["T"]) != len(self.k_table["k"]):
                raise ValueError("k_table T and k arrays must have same length")
            if self.k_table.get("interpolation") != "linear-ln-k":
                raise ValueError("k_table interpolation must be 'linear-ln-k'")
            if self.k_table.get("extrapolation") != "refuse":
                raise ValueError("k_table extrapolation must be 'refuse'")
            if any(
                not math.isfinite(float(value)) or float(value) <= 0
                for value in self.k_table["k"]
            ):
                raise ValueError("k_table rates must be finite and positive")
            temperatures = [float(value) for value in self.k_table["T"]]
            if any(
                not math.isfinite(value) for value in temperatures
            ) or temperatures != sorted(set(temperatures)):
                raise ValueError("k_table temperatures must be finite and increasing")
        if len(self.participant_site_types) != len(self.reactant_multiplicities):
            raise ValueError("participant site types and multiplicities must align")
        if self.arity != sum(self.reactant_multiplicities):
            raise ValueError("arity must equal the sum of participant multiplicities")
        if self.raw_path_degeneracy <= 0 or self.ssa_multiplier <= 0:
            raise ValueError("degeneracy and SSA multiplier must be positive")
        if len(self.atom_map) != len(set(self.atom_map.values())):
            raise ValueError("atom_map must be a bijection")
        if self.status != "enabled" and not self.status_reason:
            raise ValueError("non-enabled records require a status reason")
        if not _is_full_event_id(self.event_id):
            raise ValueError("event_id must contain a full 256-bit SHA-256 digest")
        if self.reverse_of is not None and not _is_full_event_id(self.reverse_of):
            raise ValueError("reverse_of must contain a full event id")
        expected = self._compute_event_id()
        if self.event_id != expected:
            raise ValueError("event_id does not match canonical semantic content")


def _is_full_event_id(value: Any) -> bool:
    if not isinstance(value, str) or len(value) != 68 or not value.startswith("evt_"):
        return False
    return all(character in "0123456789abcdef" for character in value[4:])


def _canonical_link_handles(value: Any, key: str | None = None) -> Any:
    """Break reverse-record hash cycles while retaining link semantics.

    The concrete full-width partner handle is checked by artifact validation.
    Its position and presence remain semantic in each record digest, represented
    by a stable relation token so reciprocal record IDs do not require finding a
    cryptographic fixed point.
    """
    if key in {"reverse_of", "reverse_event_handle", "forward_event_id"}:
        return None if value is None else "__linked_reverse_event_id__"
    if isinstance(value, dict):
        return {
            str(item_key): _canonical_link_handles(item_value, str(item_key))
            for item_key, item_value in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [_canonical_link_handles(item) for item in value]
    return value


def ssa_multiplier_for(reactant_species: list) -> float:
    """
    Return the SSA multiplier for a bimolecular reaction's reactant species.

    RMG's merged degeneracy already includes the 1/2 factor for A+A reactions.
    The stochastic simulator counts pairs as N(N-1)/2 for A+A and N_A*N_B for A+B.
    To match RMG's rate law at large N, A+A needs a multiplier of 2.0, A+B needs 1.0.

    Args:
        reactant_species: List of reactant Species objects (length 1 or 2).

    Returns:
        2.0 for A+A (same species), 1.0 for A+B (different species) or unimolecular.
    """
    if len(reactant_species) == 1:
        return 1.0
    if len(reactant_species) == 2:
        sp1, sp2 = reactant_species[0], reactant_species[1]
        if sp1 is sp2 or (hasattr(sp1, "is_isomorphic") and sp1.is_isomorphic(sp2)):
            return 2.0
        return 1.0
    return 1.0


_KNOWN_FIELDS_CLASS = {f.name for f in fields(EventRecord) if f.name != "_KNOWN_FIELDS"}
_LEGACY_REQUIRED_FIELDS = _KNOWN_FIELDS_CLASS - {
    "canonical_index",
    "site_type",
    "rate_source",
    "coproducts",
    "reverse_of",
    "inventory_class",
    "implicit_h_delta",
    "formula_delta",
    "mass_delta",
    "reactant_pair_convention",
    "resonance_derived",
    "degeneracy",
    "reactant_graphs",
    "product_graphs",
}
