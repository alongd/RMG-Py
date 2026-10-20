"""The initial melt of a linear statistical copolymer (schema v2).

A ``PolymerSpec`` declares the chain-length ensemble, the composition, the
sequence statistics, the sampler and its seed (``initial_ensemble``); this
module turns them into one residue sequence per chain, the realized
composition, and the masses the melt kMC reads in place of a repeat unit's.
Nothing here names a polymer.

**Sampler** ``markov_chain_v1``.  The chain is a first-order Markov chain:

- ``ideal_random``: states are the present residues (composition > 0) in
  name order, pi is the composition and every row of T is pi;
- ``markov1``: the same states and pi, T the declared transition rows;
- ``alternating``: states are the positions of the declared cycle, pi is
  uniform over them and T steps to the next position; a state emits the
  residue at its position.

Chains are numbered in materialization order (DP ascending, then count).
Draw ``(chain, step)`` takes u = the first 8 bytes (big-endian) of
sha256(``"markov_chain_v1:<seed>:<chain>:<step>"``) over 2^64 and returns
the first state, in state order, whose cumulative weight exceeds u (exact
rationals).  Step 0 draws from pi, step k from the row of the state drawn
at step k - 1.  Same declaration, same sequences, on every platform.

**Composition** (Q4(a)).  The realized composition is the residue count
over the whole melt's residue count, exact.  A residue declared at x > 0
is refused (``refused:composition_tolerance``) when
``|realized - x| / x`` exceeds ``composition_relative_tolerance``; a
declaration that declares no tolerance cannot be sampled.

**Masses.**  ``M_i`` is residue i's molar mass (every hydrogen counted),
``Mbar0 = sum x_i M_i`` the mean residue molar mass.  The MET length of a
component of molar mass ``M_c`` is ``max(1, floor(M_c / Mbar0))``.  The
representative masses, volume and scale factor are the v1 ones
(``LinearHomopolymerSpec.initialization``) over the sampled chains."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, Mapping, Tuple

from rmgpy import constants
from rmgpy.kmc.melt_state import MeltStateError, _element_mass, _path_keys

SAMPLER = "markov_chain_v1"
COMPOSITION_REFUSAL = "refused:composition_tolerance"
_SCALE = 2 ** 64


def _kernel(ensemble):
    """(states, pi, rows, residue of each state) of the declared statistics."""
    composition = dict(ensemble.composition)
    statistics = ensemble.statistics
    if statistics.model == "alternating":
        cycle = statistics.cycle
        states = tuple(range(len(cycle)))
        pi = {k: Fraction(1, len(cycle)) for k in states}
        rows = {k: {(k + 1) % len(cycle): Fraction(1)} for k in states}
        return states, pi, rows, dict(enumerate(cycle))
    states = tuple(sorted(r for r, x in composition.items() if x > 0))
    pi = {r: composition[r] for r in states}
    if statistics.model == "ideal_random":
        rows = {i: dict(pi) for i in states}
    elif statistics.model == "markov1":
        rows = {i: dict(row) for i, row in statistics.transition}
    else:
        raise MeltStateError(f"no sampler for the sequence model {statistics.model!r}")
    return states, pi, rows, {r: r for r in states}


def _draw(seed: int, chain: int, step: int, states, weights: Mapping) -> object:
    digest = hashlib.sha256(f"{SAMPLER}:{seed}:{chain}:{step}".encode()).digest()
    u = Fraction(int.from_bytes(digest[:8], "big"), _SCALE)
    total = Fraction(0)
    for state in states:
        total += weights.get(state, Fraction(0))
        if u < total:
            return state
    raise MeltStateError(f"the {SAMPLER} weights {dict(weights)} do not sum to one")


def sample_sequences(spec) -> Tuple[Tuple[str, ...], ...]:
    """One residue sequence per chain, in materialization order (module doc)."""
    ensemble = spec.initial_ensemble
    if ensemble.sampler != SAMPLER:
        raise MeltStateError(f"no sampler {ensemble.sampler!r}; the melt kMC samples {SAMPLER}")
    states, pi, rows, residue = _kernel(ensemble)
    out = []
    for dp, count in sorted(ensemble.chains.counts.items()):
        if dp < 1:
            raise MeltStateError(f"chains need DP >= 1, got {dp}")
        for _ in range(count):
            chain = len(out)
            state = _draw(ensemble.seed, chain, 0, states, pi)
            sequence = [residue[state]]
            for step in range(1, dp):
                state = _draw(ensemble.seed, chain, step, states, rows[state])
                sequence.append(residue[state])
            out.append(tuple(sequence))
    return tuple(out)


def realized_composition(spec, sequences=None) -> Dict[str, Fraction]:
    """Per declared residue, its count over the melt's residue count (exact)."""
    sequences = sample_sequences(spec) if sequences is None else sequences
    total = sum(len(s) for s in sequences)
    counts = {name: 0 for name, _residue in spec.residues}
    for sequence in sequences:
        for name in sequence:
            counts[name] += 1
    return {name: Fraction(n, total) if total else Fraction(0) for name, n in sorted(counts.items())}


def check_composition(spec, realized: Mapping[str, Fraction]) -> None:
    """Refuse a sample whose composition is off the declared one by more than
    the declared relative tolerance, or a declaration without one (module doc)."""
    ensemble = spec.initial_ensemble
    tolerance = ensemble.composition_tolerance
    if tolerance is None:
        raise MeltStateError(
            f"{COMPOSITION_REFUSAL}: initial_ensemble declares no composition_relative_tolerance; the melt "
            f"kMC samples residue sequences ({SAMPLER}, seed {ensemble.seed}) and accepts a sample only "
            "within a declared relative tolerance of the declared composition")
    declared = dict(ensemble.composition)
    worst = max((abs(realized[r] - x) / x for r, x in declared.items() if x > 0), default=Fraction(0))
    if worst > tolerance:
        def show(values):
            return {r: str(values[r]) for r in sorted(values)}
        raise MeltStateError(
            f"{COMPOSITION_REFUSAL}: the sampled melt ({SAMPLER}, seed {ensemble.seed}) has composition "
            f"{show(realized)} against the declared {show(declared)}: relative deviation {worst} "
            f"({float(worst):.6g}) exceeds composition_relative_tolerance {tolerance}")


def _graph_mass(graph) -> Fraction:
    """g/mol of a declared graph, every hydrogen counted (exact)."""
    total = Fraction(0)
    for atom in graph.atoms:
        total += _element_mass(atom.element, None if atom.isotope == -1 else atom.isotope)
        total += atom.hydrogens * _element_mass("H", None)
    return total * 1000


def residue_molar_mass(spec, name: str) -> Fraction:
    return _graph_mass(dict(spec.residues)[name].graph)


def mean_residue_molar_mass(spec) -> Fraction:
    """``Mbar0 = sum x_i M_i`` (g/mol, exact)."""
    return sum((x * residue_molar_mass(spec, r) for r, x in spec.initial_ensemble.composition), Fraction(0))


def residue_sites(spec, name: str):
    """(entry site, exit site) of residue ``name``."""
    from rmgpy.polymer_spec import _entry_and_exit
    entries, exits = _entry_and_exit(dict(spec.residues)[name], dict(spec.links))
    return entries[0][1], exits[0][1]


def residue_backbone_keys(spec, name: str) -> frozenset:
    """Atom keys on residue ``name``'s shortest entry-site -> exit-site path."""
    entry, exit_ = residue_sites(spec, name)
    return _path_keys(dict(spec.residues)[name].graph, entry.atom, exit_.atom)


@dataclass(frozen=True)
class SequenceInitialization:
    """The sampled melt's representative masses, volume and scale factor
    (the v1 ``SharedInitialization`` fields the melt kMC reads)."""

    sequences: Tuple[Tuple[str, ...], ...]
    representative_polymer_mass_kg: float
    representative_mobile_mass_kg: float
    volume_m3: float
    macroscopic_polymer_mass_kg: float
    scale_factor: float


def initialization(spec, sequences=None) -> SequenceInitialization:
    sequences = sample_sequences(spec) if sequences is None else sequences
    masses = {name: residue_molar_mass(spec, name) for name, _residue in spec.residues}
    ends = _graph_mass(spec.head_end_group) + _graph_mass(spec.tail_end_group)
    polymer = sum((sum((masses[r] for r in s), Fraction(0)) + ends for s in sequences), Fraction(0))
    small = {name: graph for name, graph, _disposition in spec.small_species}
    mobile = sum((count * _graph_mass(small[name]) for name, count in spec.initial_mobile_counts), Fraction(0))
    polymer_kg = float(polymer / 1000) / constants.Na
    mobile_kg = float(mobile / 1000) / constants.Na
    total_kg = float((polymer + mobile) / 1000) / constants.Na
    return SequenceInitialization(
        sequences, polymer_kg, mobile_kg, total_kg / spec.density_kg_m3, spec.initial_mass_kg,
        spec.initial_mass_kg / polymer_kg if polymer_kg else 0.0)
