from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field, fields
from typing import Any, Literal, Optional
from hashlib import sha256


@dataclass
class EventRecord:
    """
    Draft schema rewrite_record/0.1 for a compiled kMC event.

    Strand-specific fields (cut_offset, feature_ops, inheritance, junction_ops)
    are documented null placeholders; the compiler will define them later.
    """
    schema_version: str = "rewrite_record/0.1"
    event_id: str = ""
    family: str = ""
    template: str = ""
    arity: int = 0
    participant_site_types: list[str] = field(default_factory=list)
    reactant_multiplicities: list[int] = field(default_factory=list)
    raw_path_degeneracy: float = 1.0
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
            self.event_id = self._compute_event_id()

    def _compute_event_id(self) -> str:
        parts = [
            self.family,
            self.template,
            str(self.arity),
            ",".join(self.participant_site_types),
            ",".join(str(m) for m in self.reactant_multiplicities),
            json.dumps(self.atom_map, sort_keys=True),
            json.dumps(self.bond_ops, sort_keys=True),
        ]
        h = sha256("|".join(parts).encode()).hexdigest()[:16]
        return f"evt_{h}"

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

        missing = _KNOWN_FIELDS_CLASS - set(data.keys())
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