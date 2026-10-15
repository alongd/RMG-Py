"""Recertified closed windows, chain-end contexts and the short-chain ladder
(I-137).  Skeleton: the public API the tests target; bodies follow."""

from dataclasses import dataclass, field
from typing import Any, Callable, Optional

RING_LIMIT = 8
PORT_DEPTH_LIMIT = 6
FIXED_POINT_RETRIES = 3
LOOP_STATUSES = ("fixed_point", "no_fixed_point", "ceiling_exceeded", "certificate_unstable")


class WindowError(ValueError):
    """Bad input."""


@dataclass(frozen=True)
class Support:
    atoms: frozenset
    reacting: frozenset = frozenset()
    sources: tuple = ()
    guards: tuple = ()
    searches: tuple = ()
    refusal: Optional[str] = None


@dataclass(frozen=True)
class Window:
    atoms: frozenset
    support: Support
    added: tuple = ()
    refusal: Optional[str] = None


@dataclass(frozen=True)
class Observation:
    support: Support
    viewer: Callable = field(compare=False)
    certificate: Any = field(default=None, compare=False)
    refusal: Optional[str] = None

    def view(self, inside):
        return self.viewer(inside)


@dataclass(frozen=True)
class FamilyConstants:
    family: str
    radius: int
    reacting_distance: int
    free_label_search: bool
    bimolecular: bool = False


def fixed_point(witness, observe, providers=None, *, ceiling=None, port_limit=PORT_DEPTH_LIMIT,
                retries=FIXED_POINT_RETRIES):
    raise NotImplementedError


def close(witness, support, *, ring_limit=RING_LIMIT, sides=("reactants", "products")):
    raise NotImplementedError


def context_assignments(witness, near_ports=None, providers=None):
    raise NotImplementedError


@dataclass(frozen=True)
class ContextResult:
    assignment: tuple
    status: str
    reason: str = ""
    witness: Any = field(default=None, compare=False)
    window: Any = field(default=None, compare=False)
    observation: Any = field(default=None, compare=False)


def merge_contexts(results):
    raise NotImplementedError


def span_bound(constants):
    raise NotImplementedError


def ladder_plan(constants, bound=None):
    raise NotImplementedError


def chain_backbone_atoms(spec):
    raise NotImplementedError
