"""Observables of a network melt (``NetworkPolymerSpec``, S7).

Every quantity is read from the atoms of a ``melt_state.MeltState`` (a
molecule is a connected component), in representative molar-mass units
(g/mol summed over the representative melt; times scale / N_A for grams):

* **volatile**: a component isomorphic to a declared species
  (``small_species``, the monomer products, the built-ins H and H2:
  ``melt_engine.SpeciesIdentity``) has left the condensed melt;
* **condensed**: every other component; the **gel** is the component with the
  most backbone-origin atoms (the least id on a tie), as the engine reads it
  (``MeltEngine.gel_component``) -- a finite-ensemble convention: the largest
  molecule is reported as the gel whatever its size, so the sol's largest
  molecule (``sol_molecule_sizes``) is reported beside it; the **sol** is the
  rest of the condensed melt;
* **inert**: each declared inert phase's mass, ``f_i / (1 - sum f)`` times the
  initial resin mass (the declaration's fractions are of the whole
  composite), constant: it takes no part in the chemistry;
* ``gel_mass_fraction`` / ``sol_mass_fraction``: of the condensed resin mass;
* ``residue_mass_fraction``: (condensed + inert) over the initial composite
  mass (resin at t = 0, polymer and mobile, plus inert);
  ``mass_loss_fraction`` is one minus it, so a filler is residue and lowers
  the composite's loss by its share; mass that leaves the state (an atom
  removed) is loss.

``observed_run`` records these at given times with the engine's leak per
class (``MeltEngine.leak_by_class``).  Nothing here names a polymer."""

from __future__ import annotations

import math
from collections import Counter
from fractions import Fraction
from typing import Dict, Iterable, List, Optional

NETWORK_OBSERVABLES_REFUSAL = "refused:network_observables"


def _require_network(spec):
    from rmgpy.polymer_spec import NetworkPolymerSpec
    if not isinstance(spec, NetworkPolymerSpec):
        raise TypeError(f"{NETWORK_OBSERVABLES_REFUSAL}: network observables read a NetworkPolymerSpec, got "
                        f"{type(spec).__name__}")


def initial_resin_mass(spec) -> Fraction:
    """The representative initial resin mass (g/mol summed): every residue
    with its vacant sites capped, plus the initial mobile species."""
    from rmgpy.kmc.melt_branched import mean_residue_molar_mass
    from rmgpy.kmc.melt_sequence import _graph_mass
    residues = sum(dict(spec.initial_ensemble.residue_counts).values())
    small = {name: graph for name, graph, _d in spec.small_species}
    mobile = sum((count * _graph_mass(small[name]) for name, count in spec.initial_mobile_counts), Fraction(0))
    return mean_residue_molar_mass(spec) * residues + mobile


def inert_masses(spec) -> Dict[str, Fraction]:
    """Each inert phase's representative mass (g/mol summed), exact."""
    resin = initial_resin_mass(spec)
    rest = 1 - sum((p.mass_fraction for p in spec.inert_phases), Fraction(0))
    return {p.name: p.mass_fraction / rest * resin for p in spec.inert_phases}


def _units(state, cid) -> int:
    from rmgpy.kmc.melt_state import segment_residue
    units = set()
    for a in state.components.atoms(cid):
        origin = state.origin[a]
        if origin.kind == "chain" and state.atoms[a].element != "H" and segment_residue(origin.segment):
            units.add((origin.ordinal, origin.segment))
    return len(units)


def gel_component(state) -> int:
    comps = state.components
    return min(comps.ids(), key=lambda c: (-comps.backbone_count(c), c))


def observe(spec, state, identity=None) -> dict:
    """The network observables of ``state`` (module doc)."""
    _require_network(spec)
    if identity is None:
        from rmgpy.kmc.melt_engine import SpeciesIdentity
        identity = SpeciesIdentity(spec)
    comps = state.components
    gel = gel_component(state)
    volatile, condensed, gel_mass = Counter(), [], 0.0
    sol_sizes = Counter()
    for cid in comps.ids():
        mass = float(comps.molar_mass_g_mol(cid))
        name = identity.identify(state, cid)
        if name is not None:
            volatile[name] += mass
            continue
        condensed.append(mass)
        if cid == gel:
            gel_mass = mass
        else:
            sol_sizes[_units(state, cid)] += 1
    condensed_mass = math.fsum(condensed)
    inert = {name: float(m) for name, m in inert_masses(spec).items()}
    inert_total = math.fsum(inert.values())
    total0 = float(initial_resin_mass(spec)) + inert_total
    residue = (condensed_mass + inert_total) / total0
    return {
        "gel_component": gel,
        "gel_residues": _units(state, gel),
        "gel_mass_fraction": gel_mass / condensed_mass if condensed_mass else 0.0,
        "sol_mass_fraction": (condensed_mass - gel_mass) / condensed_mass if condensed_mass else 0.0,
        "sol_molecule_sizes": dict(sorted(sol_sizes.items())),
        "condensed_mass_g_mol": condensed_mass,
        "volatile_mass_g_mol": math.fsum(volatile.values()),
        "volatile_by_species_g_mol": dict(sorted(volatile.items())),
        "inert_mass_g_mol": inert_total,
        "inert_by_phase_g_mol": dict(sorted(inert.items())),
        "inert_composition": {p.name: dict(p.elements) for p in spec.inert_phases},
        "initial_resin_mass_g_mol": float(initial_resin_mass(spec)),
        "residue_mass_fraction": residue,
        "mass_loss_fraction": 1.0 - residue,
    }


def observed_run(engine, times: Iterable[float]) -> List[dict]:
    """Run ``engine`` (a network ``MeltEngine``) through ``times`` (ascending)
    and record the observables and the leak per class at each."""
    _require_network(engine.spec)
    rows = []
    for t in sorted(float(x) for x in times):
        if t > engine.time:
            engine.run(until=t)
        row = {"time": t, "events": engine.events_fired, "leak_num": engine.leak_num,
               "leak_by_class": dict(engine.leak_by_class)}
        row.update(observe(engine.spec, engine.state, engine.identity))
        rows.append(row)
    return rows
