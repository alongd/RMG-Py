"""S5 observables of the atom-native melt, comparable with the method of moments.

Every quantity is computed from the atoms of a ``melt_state.MeltState`` (a
molecule is a connected component) with the definitions the method-of-moments
(MoM) side uses (design I-127 v3 §7), exactly where the two can agree and with
a stated bound where their conventions differ.  Nothing here names a polymer:
``b`` (backbone atoms per repeat unit), the repeat-unit and end-group
compositions and the scale factor are read from the declaration.

Per component c (exact, from atoms)
-----------------------------------
* ``comp_c``: integer atom counts per (element, isotope);
* ``M_c = Composition(comp_c).molar_mass_g_mol()`` (``rmgpy.polymer_mass``);
* ``B_c``: atoms whose immutable origin is backbone (the repeat unit's
  head-port -> tail-port path; an initially mobile monomer gets it through its
  unique map onto the repeat unit, ``melt_state.monomer_backbone_keys``);
* ``n_c = floor(B_c / b)``, an integer;
* ``caps_c``: the atoms whose origin is a declared end group (head or tail);
  ``E_c = mass(caps_c) / e`` is the component's physical end-group inventory
  in end pairs (``e`` = head + tail end-group molar mass).

Classes (run-time rule, ``SPECIES_RULE``)
-----------------------------------------
A component isomorphic (under resonance, all H explicit) to a declared species
-- ``spec.small_species``, ``monomer_product``, the built-ins H and H2 -- IS
that species, whatever its ancestry: a fragment cut from a chain that matches
a declared species leaves the polymer class.  Any other component is
*polymer-derived* when ``n_c >= 1`` (it carries at least one whole repeat
unit's backbone: MoM's chains have DP >= 1, explicit ladder 1..xs, tail
>= xs+1) and an *undeclared light product* when ``n_c = 0``.  A declaration
whose declared species is isomorphic to an initial chain, or whose small
species label equals a built-in label, is refused (``melt_engine``).  The
release report lists every undeclared component outside the tail
(``n_c <= xs``) by canonical key (the smallest SMILES over its resonance
forms); nothing is refused at run time.

Tiers (design §7)
-----------------
1. Total mass per class and the number of polymer-derived molecules.  The
   ``mom_definition`` values are an accounting identity on the kMC partition,
   not a comparison with a MoM run: the MoM mass law
   (``EvolvedPoolMass``: repeat x ``mu1`` + end-group inventory ``E`` +
   signed defect inventory) and the MoM molecule count
   (``mu0_tail + amount(explicit)``) are rebuilt from the kMC components and
   must equal the kMC sums (to rounding).  The comparison against MoM output
   is ``side_by_side`` with the MoM run's numbers.
2. Mass moments over the polymer-derived molecules, exact
   (``sum amount M_c^k``), beside the MoM convention
   (``observables_from_moments`` of the integer DP moments, which assigns
   every chain exactly ``a n_c + e``).  ``delta_c = M_c - (a n_c + e)``; the
   report gives the exact differences and the bounds ``sum |delta_c|`` (M1)
   and ``sum |delta_c| (2 (a n_c + e) + |delta_c|)`` (M2).
3. Tail / explicit partition with the integer ``n_c`` (tail iff ``n_c > xs``,
   so the tail support invariant ``mu1 >= (xs + 1) mu0`` holds); diagnostic
   only.

MoM mass law on the exported state (``mom_mass_law``)
-----------------------------------------------------
The exported MoM state is the tail moments ``mu``, the PHYSICAL end-group
inventory of the tail ``E = sum_tail E_c`` (from atom origins: a bridge
scission creates no end group, so it leaves ``E`` unchanged) and the explicit
amounts per DP.  MoM's mass law on it is
``a mu1 + e E + sum_n amount_n (a n + e)`` (explicit DP-n species at the
capped DP-n molar mass).  Its gap to the kMC polymer-derived mass is exact
and signed: ``sum_tail mass(comp_c - n_c repeat - caps_c) + sum_explicit
mass(comp_c - n_c repeat - ends)`` (zero after any bridge scission or
recombination of whole units; nonzero for a split inside a repeat unit, a
radical end lacking a cap H in the explicit ladder, a transferred H), bounded
by the same sum of absolute values (``bound_g``).

Units: ``amount = scale N / N_A`` (mol) with ``scale`` the declaration's
``initialization().scale_factor``; the volume is constant
(``volume_m3 = scale x`` the representative volume).

Exact invariants (``InvariantChecker``): after every event, ``sum comp`` and
``sum B`` are constant; the molecule count changes by ``|products| -
|reactants|``; a bridge scission adds ``-2 M_a M_b`` to ``sum M^2`` and a
recombination ``+2 M_a M_b``; components outside the event are unchanged.
Given the other checks the ``sum M^2`` increment of the index's stored masses
is implied algebraically, so the expected increment is computed from masses
recomputed from each component's composition (exact rationals): the line then
catches stored-mass drift between the event's components.

Export: ``write_time_series`` writes one JSON line per requested time (the
engine's horizon API), ``mom_reported_quantities`` maps one observation to the
quantities the MoM solver reports, ``side_by_side`` compares them with a MoM
run's own numbers.
"""

from __future__ import annotations

import json
import math
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, Iterable, List, Mapping, Optional, Tuple

from rmgpy import constants
from rmgpy.kmc.melt_engine import (BUILTIN_SPECIES, MeltEngine, SpeciesIdentity,
                                   refuse_initial_species_overlap)
from rmgpy.kmc.melt_state import _Builder, _backbone_keys, _element_mass
from rmgpy.polymer_mass import (Composition, EvolvedPoolMass, ExplicitInventory, InventoryTerm,
                                observables_from_moments)

OBSERVABLES_SCHEMA = "kmc-melt-observables/2"
UNITS = {"amount": "mol (scale N / N_A)", "mass": "g", "molar_mass": "g/mol", "M2": "g^2/mol",
         "dp_moments": "mol, mol, mol (sum of amount n^k, k = 0, 1, 2)", "time": "s", "volume": "m^3",
         "end_inventory": "mol of end-group pairs (head + tail)"}
CONVENTION_DEFINITION = (
    "delta_c = M_c - (a n_c + e) = molar mass of comp_c - n_c repeat - ends (signed); "
    "M1_difference = sum amount delta_c, M1_bound = sum amount |delta_c|; "
    "M2_difference = sum amount delta_c (2 (a n_c + e) + delta_c), "
    "M2_bound = sum amount |delta_c| (2 (a n_c + e) + |delta_c|); exact differences, not tolerances")
MASS_LAW_DEFINITION = (
    "value_g = a mu1_tail + e E + sum_n amount_n (a n + e), E = physical end-group inventory of the tail "
    "(mass of end-group-origin atoms / e); difference_g = kMC polymer-derived mass - value_g = "
    "amount x [sum_tail mass(comp_c - n_c repeat - caps_c) + sum_explicit mass(comp_c - n_c repeat - ends)]; "
    "bound_g = the same sum of absolute values")
TIER1_NOTE = (
    "mom_definition is an accounting identity on the kMC partition (the MoM mass law and molecule count "
    "rebuilt from the kMC components), not a comparison with a MoM run; compare MoM output with side_by_side")
SPECIES_RULE = (
    "a component isomorphic under resonance (all H explicit) to a declared species (small species, "
    "monomer_product, built-in H and H2) is that declared species whatever its ancestry; any other component "
    "is polymer-derived iff n_c >= 1, else an undeclared light product; a declared species isomorphic to an "
    "initial chain, or a small-species label equal to a built-in label, is refused at run start")
END_SEGMENTS = ("head_end", "tail_end")


class ObservablesError(ValueError):
    """An observable cannot be formed (bad cutoff, time, declaration)."""


class InvariantError(ObservablesError):
    """An exact per-event invariant of the melt is broken."""


# --------------------------------------------------------------------------
# compositions
# --------------------------------------------------------------------------
def _key(counts: Mapping) -> Tuple:
    """Canonical composition tuple ((element, isotope or None), count), sorted, nonzero."""
    return tuple(sorted(((k, int(v)) for k, v in counts.items() if v), key=lambda kv: (kv[0][0], kv[0][1] or 0)))


def _from_declared(composition: Composition) -> Counter:
    return Counter({(symbol, None if isotope == -1 else isotope): count
                    for symbol, isotope, count in composition.entries})


def _minus(counts: Mapping, *parts: Tuple[int, Mapping]) -> Tuple:
    rest = Counter(counts)
    for times, part in parts:
        for k, v in part.items():
            rest[k] -= times * v
    return _key(rest)


_MASSES: Dict[Tuple, float] = {}


def molar_mass_g_mol(composition: Tuple) -> float:
    """``Composition(comp).molar_mass_g_mol()`` (the shared mass routine), cached;
    signed compositions give signed masses."""
    if composition not in _MASSES:
        _MASSES[composition] = (Composition.from_mapping(dict(composition)).molar_mass_g_mol()
                                if composition else 0.0)
    return _MASSES[composition]


def backbone_atoms_per_unit(spec) -> int:
    """``b``: atoms on the declared repeat unit's head-port -> tail-port path."""
    return len(_backbone_keys(spec.repeat_unit))


def canonical_key(state, cid) -> str:
    """The smallest SMILES over the component's resonance forms (a canonical,
    representation-independent key for reporting)."""
    from rmgpy.kmc import window_matcher as wm
    from rmgpy.kmc.rmg_moves import molecule_from_rows
    keys = set()
    for form in wm.component_forms(state, cid):
        keys.add(molecule_from_rows(wm.form_rows(state, cid, form)).to_smiles())
    return min(keys)


@dataclass(frozen=True)
class ComponentRecord:
    """Exact observables of one component."""

    cid: int
    composition: Tuple        # ((element, isotope or None), count), sorted
    molar_mass_g_mol: float   # Composition(composition).molar_mass_g_mol()
    backbone: int             # B_c
    n: int                    # floor(B_c / b)
    species: Optional[str]    # the declared species it is, or None
    category: str             # "species", "polymer" or "light"
    caps: Tuple               # composition of the end-group-origin atoms
    defect: Tuple             # polymer: comp - n repeat - ends (signed); else ()
    law_defect: Tuple         # polymer: comp - n repeat - caps (signed); else ()


def _fsum(values) -> float:
    return math.fsum(values)


def _rel(a: float, b: float) -> float:
    scale = max(abs(a), abs(b))
    return abs(a - b) / scale if scale else 0.0


def declared_species_names(spec, identity: SpeciesIdentity) -> List[str]:
    """Every name the identity can return: small species, ``monomer_product``
    (when it is not also a small species) and the built-ins."""
    names = [name for name, _g, _d in spec.small_species]
    monomer = getattr(spec, "monomer_product", None)
    if monomer is not None:
        builder = _Builder()
        builder.add_segment(monomer, "mobile", 0, "monomer_product")
        alone = builder.state()
        if identity.identify(alone, alone.components.ids()[0]) == "monomer_product":
            names.append("monomer_product")
    return names + list(BUILTIN_SPECIES)


# --------------------------------------------------------------------------
# the observer
# --------------------------------------------------------------------------
class MeltObserver:
    """Tier 1-3 observables of a melt declared by ``spec`` with the MoM cutoff
    ``xs`` (chains with ``n_c <= xs`` are the explicit ladder, ``n_c > xs`` the
    moment tail).  The declaration is checked like at run start: a declared
    species isomorphic to an initial chain, or a small-species label equal to
    a built-in label, is refused (``DeclarationError``)."""

    def __init__(self, spec, *, xs: int, identity: Optional[SpeciesIdentity] = None):
        if isinstance(xs, bool) or not isinstance(xs, int) or xs < 0:
            raise ObservablesError(f"xs must be a non-negative integer, got {xs!r}")
        self.spec = spec
        self.xs = xs
        self.b = backbone_atoms_per_unit(spec)
        self.repeat = _from_declared(spec.repeat_unit.composition)
        self.ends = _from_declared(spec.end_composition)
        self.a = spec.repeat_unit.composition.molar_mass_g_mol()
        self.e = spec.end_composition.molar_mass_g_mol()
        init = spec.initialization()
        self.scale = init.scale_factor
        if not (math.isfinite(self.scale) and self.scale > 0.0):
            raise ObservablesError(f"the declaration's scale factor {self.scale!r} is not positive")
        self.unit = self.scale / constants.Na
        self.representative_volume_m3 = init.volume_m3
        self.initial_polymer_mass_g = spec.initial_mass_kg * 1000.0
        self.identity = identity if identity is not None else SpeciesIdentity(spec)
        refuse_initial_species_overlap(spec, self.identity)
        self.declared = declared_species_names(spec, self.identity)
        declared = {name: graph.composition for name, graph, _d in spec.small_species}
        if "monomer_product" in self.declared:
            declared["monomer_product"] = spec.monomer_product.composition
        declared["H"] = Composition.from_mapping({"H": 1})
        declared["H2"] = Composition.from_mapping({"H": 2})
        self.declared_molar_mass = {k: declared[k].molar_mass_g_mol() for k in self.declared}

    # -- per component ----------------------------------------------------
    def components(self, state, engine: Optional[MeltEngine] = None) -> List[ComponentRecord]:
        comps = state.components
        origin = state.origin
        nodes = state.atoms
        synced = (engine is not None and engine.state is state
                  and engine._revision == state.revision and set(engine.entries) == set(comps.ids()))
        out = []
        for cid in comps.ids():
            counts = comps.composition(cid)
            composition = _key(counts)
            backbone = comps.backbone_count(cid)
            n = backbone // self.b
            caps = Counter()
            for x in comps.atoms(cid):
                o = origin[x]
                if o.kind == "chain" and o.segment in END_SEGMENTS:
                    caps[(nodes[x].element, nodes[x].isotope)] += 1
            caps = _key(caps)
            species = engine.entries[cid].species if synced else self.identity.identify(state, cid)
            if species is not None:
                category, defect, law_defect = "species", (), ()
            elif n >= 1:
                category = "polymer"
                defect = _minus(counts, (n, self.repeat), (1, self.ends))
                law_defect = _minus(counts, (n, self.repeat), (1, dict(caps)))
            else:
                category, defect, law_defect = "light", (), ()
            out.append(ComponentRecord(cid, composition, molar_mass_g_mol(composition), backbone, n,
                                       species, category, caps, defect, law_defect))
        return out

    # -- one observation --------------------------------------------------------
    def observe(self, target) -> dict:
        """One observation (a JSON-ready dict) of an engine or a bare state."""
        engine = target if isinstance(target, MeltEngine) else None
        state = engine.state if engine is not None else target
        records = self.components(state, engine)
        u = self.unit
        a, e = self.a, self.e
        polymer = [r for r in records if r.category == "polymer"]
        light = [r for r in records if r.category == "light"]
        tail = [r for r in polymer if r.n > self.xs]
        explicit = [r for r in polymer if r.n <= self.xs]

        def moments(rs):
            return [len(rs), sum(r.n for r in rs), sum(r.n * r.n for r in rs)]

        rep_all, rep_tail = moments(polymer), moments(tail)
        mol_all = [u * x for x in rep_all]
        mol_tail = [u * x for x in rep_tail]
        end_pairs = _fsum(molar_mass_g_mol(r.caps) for r in tail) / e if e else 0.0
        end_mol = u * end_pairs

        # tier 1: accounting identity (MoM mass law and count rebuilt from the kMC partition)
        defects = Counter(r.law_defect for r in tail if r.law_defect)
        pool = EvolvedPoolMass(
            "polymer", Composition.from_mapping(dict(self.repeat)), tuple(mol_tail),
            end_inventory=((InventoryTerm("ends", Composition.from_mapping(dict(self.ends)), end_mol, True),)
                           if end_mol else ()),
            defects=tuple(InventoryTerm(f"defect {i}", Composition.from_mapping(dict(d)), u * k)
                          for i, (d, k) in enumerate(sorted(defects.items()))))
        by_comp = Counter(r.composition for r in explicit)
        explicit_inventory = [ExplicitInventory(f"explicit {i}", Composition.from_mapping(dict(c)), u * k)
                              for i, (c, k) in enumerate(sorted(by_comp.items()))]
        mom_mass = _fsum([pool.mass_g] + [x.mass_g for x in explicit_inventory])
        mom_amount = _fsum([mol_tail[0]] + [x.amount_mol for x in explicit_inventory])
        species = {}
        for name in self.declared:
            rs = [r for r in records if r.species == name]
            amount = u * len(rs)
            species[name] = {"molecules": len(rs), "amount_mol": amount,
                             "mass_g": u * _fsum(r.molar_mass_g_mol for r in rs),
                             "mom_definition": {"amount_mol": amount,
                                                "mass_g": amount * self.declared_molar_mass[name]}}
        polymer_mass = u * _fsum(r.molar_mass_g_mol for r in polymer)
        tier1 = {
            "mom_definition_note": TIER1_NOTE,
            "polymer_derived": {"molecules": len(polymer), "amount_mol": u * len(polymer),
                                "mass_g": polymer_mass,
                                "mom_definition": {"amount_mol": mom_amount, "mass_g": mom_mass}},
            "species": species,
            "undeclared_light": {"molecules": len(light), "amount_mol": u * len(light),
                                 "mass_g": u * _fsum(r.molar_mass_g_mol for r in light)},
            "total_mass_g": u * _fsum(r.molar_mass_g_mol for r in records),
        }

        # MoM mass law on the exported state (mu_tail, physical E, explicit amounts)
        ladder = Counter(r.n for r in explicit)
        ladder_mass = Counter()
        for r in explicit:
            ladder_mass[r.n] += r.molar_mass_g_mol
        value = _fsum([a * mol_tail[1], e * end_mol] + [u * k * (a * n + e) for n, k in ladder.items()])
        gaps = ([molar_mass_g_mol(r.law_defect) for r in tail]
                + [molar_mass_g_mol(r.defect) for r in explicit])
        mass_law = {"definition": MASS_LAW_DEFINITION, "end_inventory_mol": end_mol,
                    "tail_repeat_mass_g": a * mol_tail[1], "tail_end_mass_g": e * end_mol,
                    "explicit_capped_mass_g": _fsum(u * k * (a * n + e) for n, k in ladder.items()),
                    "value_g": value, "kmc_mass_g": polymer_mass,
                    "difference_g": u * _fsum(gaps), "bound_g": u * _fsum(abs(g) for g in gaps)}

        # tier 2: mass moments over the polymer-derived molecules
        m = [r.molar_mass_g_mol for r in polymer]
        m1, m2 = _fsum(m), _fsum(x * x for x in m)
        conv = observables_from_moments(mol_all, a, e)
        deltas = [(molar_mass_g_mol(r.defect), a * r.n + e) for r in polymer]
        tier2 = {
            "M0_mol": u * len(polymer), "M1_g": u * m1, "M2_g2_mol": u * m2,
            "Mn_g_mol": m1 / len(polymer) if polymer else 0.0, "Mw_g_mol": m2 / m1 if m1 else 0.0,
            "dp_moments_representative": rep_all, "dp_moments_mol": mol_all,
            "mom_convention": {"mass_g": conv.mass_g, "Mn_g_mol": conv.mn_g_mol, "Mw_g_mol": conv.mw_g_mol,
                               "M2_g2_mol": conv.mw_g_mol * conv.mass_g},
            "end_mass_convention": {
                "definition": CONVENTION_DEFINITION,
                "delta_min_g_mol": min((d for d, _b in deltas), default=0.0),
                "delta_max_g_mol": max((d for d, _b in deltas), default=0.0),
                "M1_difference_g": u * _fsum(d for d, _b in deltas),
                "M1_bound_g": u * _fsum(abs(d) for d, _b in deltas),
                "M2_difference_g2_mol": u * _fsum(d * (2.0 * b + d) for d, b in deltas),
                "M2_bound_g2_mol": u * _fsum(abs(d) * (2.0 * b + abs(d)) for d, b in deltas),
            },
        }

        # tier 3: tail / explicit partition (diagnostic)
        tier3 = {
            "xs": self.xs, "diagnostic_only": True,
            "tail": {"molecules": len(tail), "dp_moments_representative": rep_tail, "dp_moments_mol": mol_tail,
                     "end_inventory_mol": end_mol,
                     "support_invariant_holds": rep_tail[1] >= (self.xs + 1) * rep_tail[0]},
            "explicit": {"dp_molecules": {str(n): k for n, k in sorted(ladder.items())},
                         "dp_amount_mol": {str(n): u * k for n, k in sorted(ladder.items())},
                         "dp_mass_g": {str(n): u * ladder_mass[n] for n in sorted(ladder)}},
        }

        # release report: undeclared, not in the tail, by canonical key
        undeclared: Dict[Tuple, dict] = {}
        for r in records:
            if r.species is None and r.n <= self.xs:
                key = canonical_key(state, r.cid)
                slot = undeclared.setdefault((r.category, r.n, key), {
                    "key": key, "molecules": 0, "n": r.n, "category": r.category,
                    "molar_mass_g_mol": r.molar_mass_g_mol})
                slot["molecules"] += 1
        release = {"undeclared": [undeclared[k] for k in sorted(undeclared)]}

        total_comp = Counter()
        for r in records:
            total_comp.update(dict(r.composition))
        counts = {"events": 0, "null_proposals": 0, "family": {}, "template": {}, "met": {}}
        leak = None
        time_s = 0.0
        temperature = None
        rep_volume = self.representative_volume_m3
        if engine is not None:
            counts = {"events": engine.events_fired, "null_proposals": engine.null_proposals,
                      "family": dict(sorted(engine.family_counts.items())),
                      "template": dict(sorted(engine.template_counts.items())),
                      "met": dict(sorted(engine.met_counts.items()))}
            leak = {k: v for k, v in engine.leak_estimate().items() if k != "definition"}
            time_s, temperature, rep_volume = engine.time, engine.temperature, engine.volume
        return {
            "schema": OBSERVABLES_SCHEMA, "time_s": time_s, "engine_time_s": time_s,
            "temperature_K": temperature, "units": UNITS,
            "scale_factor": self.scale, "avogadro": constants.Na, "xs": self.xs,
            "backbone_atoms_per_unit": self.b, "repeat_unit_molar_mass_g_mol": a,
            "end_group_molar_mass_g_mol": e, "initial_polymer_mass_g": self.initial_polymer_mass_g,
            "representative_volume_m3": rep_volume, "volume_m3": self.scale * rep_volume,
            "declared_species": list(self.declared), "species_rule": SPECIES_RULE,
            "tier1": tier1, "mom_mass_law": mass_law, "tier2": tier2, "tier3": tier3, "release": release,
            "counts": counts, "leak_estimate": leak,
            "invariants": {"components": len(records), "backbone_atoms": sum(r.backbone for r in records),
                           "composition": {f"{el}{iso or ''}": k for (el, iso), k in _key(total_comp)}},
        }


def observe(target, spec=None, *, xs: int) -> dict:
    """One observation of an engine (or of a bare state with its ``spec``)."""
    if spec is None:
        if not isinstance(target, MeltEngine):
            raise ObservablesError("a bare state needs its declaration")
        spec = target.spec
    return MeltObserver(spec, xs=xs).observe(target)


# --------------------------------------------------------------------------
# exact invariants
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Totals:
    """Exact melt totals over the component index."""

    count: int
    mass: Fraction        # sum M_c (g/mol, exact)
    mass2: Fraction       # sum M_c^2
    composition: Tuple    # sum comp_c
    backbone: int         # sum B_c


def _snapshot(state) -> Dict[int, Tuple]:
    comps = state.components
    return {cid: (comps.molar_mass_g_mol(cid), _key(comps.composition(cid)), comps.backbone_count(cid))
            for cid in comps.ids()}


def _totals(snapshot) -> Totals:
    comp = Counter()
    mass = mass2 = Fraction(0)
    backbone = 0
    for m, c, b in snapshot.values():
        comp.update(dict(c))
        mass += m
        mass2 += m * m
        backbone += b
    return Totals(len(snapshot), mass, mass2, _key(comp), backbone)


def _composition_mass(composition: Tuple) -> Fraction:
    """Exact molar mass (g/mol) recomputed from a composition (the element
    masses the state's index uses)."""
    return sum((_element_mass(el, iso) * k for (el, iso), k in composition), Fraction(0)) * 1000


def exact_totals(state) -> Totals:
    return _totals(_snapshot(state))


class InvariantChecker:
    """Checks the exact invariants after every event; usable as
    ``engine.run(check=checker)`` (called as ``checker(engine, event)``) or
    directly as ``checker.check(state, event)`` (``event=None``: totals only).
    ``counts`` tallies split / merge / other events."""

    def __init__(self, state):
        self._prev = _snapshot(state)
        self._tot = _totals(self._prev)
        self.checks = 0
        self.counts: Counter = Counter()

    def __call__(self, engine, event):
        self.check(engine.state, event)

    def check(self, state, event=None):
        snap = _snapshot(state)
        tot = _totals(snap)
        prev, ptot = self._prev, self._tot
        if tot.composition != ptot.composition:
            raise InvariantError(f"sum of compositions changed: {ptot.composition} -> {tot.composition}")
        if tot.backbone != ptot.backbone:
            raise InvariantError(f"sum of backbone counts changed: {ptot.backbone} -> {tot.backbone}")
        if tot.mass != ptot.mass:
            raise InvariantError(f"sum of masses changed: {ptot.mass} -> {tot.mass}")
        if event is not None:
            reactants, products = tuple(event.components), tuple(event.products)
            if any(c not in prev for c in reactants) or any(c not in snap for c in products):
                raise InvariantError(f"event components {reactants} -> {products} are not in the index")
            outside_before = {c: v for c, v in prev.items() if c not in reactants}
            outside_after = {c: v for c, v in snap.items() if c not in products}
            if outside_before != outside_after:
                raise InvariantError("components outside the event changed")
            expected = len(products) - len(reactants)
            if tot.count - ptot.count != expected:
                raise InvariantError(f"molecule count changed by {tot.count - ptot.count}, the event "
                                     f"implies {expected}")
            # expected increment from masses recomputed from each component's
            # composition (not the index's stored mass): with sum comp conserved
            # and the components outside the event unchanged, the increment of
            # the stored masses is implied algebraically, so only this
            # independent source makes the M^2 line a check (index mass drift)
            m_r = [_composition_mass(prev[c][1]) for c in reactants]
            m_p = [_composition_mass(snap[c][1]) for c in products]
            if len(reactants) == 1 and len(products) == 2:
                kind, delta = "split", -2 * m_p[0] * m_p[1]
            elif len(reactants) == 2 and len(products) == 1:
                kind, delta = "merge", 2 * m_r[0] * m_r[1]
            else:
                kind, delta = "other", sum(m * m for m in m_p) - sum(m * m for m in m_r)
            if tot.mass2 - ptot.mass2 != delta:
                raise InvariantError(f"sum M^2 changed by {tot.mass2 - ptot.mass2}, a {kind} implies {delta}")
            self.counts[kind] += 1
        self._prev, self._tot = snap, tot
        self.checks += 1


# --------------------------------------------------------------------------
# export
# --------------------------------------------------------------------------
def write_time_series(engine, times: Iterable[float], path, *, xs: int, check: bool = True,
                      observer: Optional[MeltObserver] = None, append: bool = False) -> List[dict]:
    """Advance ``engine`` to each requested time (``run(until=t)``, the horizon
    API: observing never changes the trajectory) and write one observation per
    time as a JSON line; with ``check`` the exact invariants are asserted after
    every event (``invariant_checks`` counts them).  Returns the observations."""
    observer = observer if observer is not None else MeltObserver(engine.spec, xs=xs)
    checker = InvariantChecker(engine.state) if check else None
    times = [float(t) for t in times]
    if any(not math.isfinite(t) for t in times) or times != sorted(times) or (times and times[0] < engine.time):
        raise ObservablesError(f"observation times must be finite, sorted and >= the engine time {engine.time}")
    out = []
    with open(path, "a" if append else "w") as handle:
        for t in times:
            if t > engine.time:
                engine.run(until=t, check=checker)
            obs = observer.observe(engine)
            obs["time_s"] = t
            obs["invariant_checks"] = checker.checks if checker is not None else None
            handle.write(json.dumps(obs, sort_keys=True, allow_nan=False) + "\n")
            handle.flush()
            out.append(obs)
    return out


def mom_reported_quantities(obs: Mapping) -> dict:
    """The kMC state of one observation expressed in the quantities the MoM
    solver reports, by the MoM definitions: tail moments (pool moment slots),
    the physical end-group inventory of the tail, explicit-ladder amounts per
    DP, the condensed polymer mass by MoM's mass law on those
    (``mom_mass_law.value_g``), Mn / Mw of the tail moments
    (``observables_from_moments``), polymer conversion, the amount of every
    declared species (0 when absent) and the polymer-derived molecules."""
    t1, t3, law = obs["tier1"], obs["tier3"], obs["mom_mass_law"]
    mu = list(t3["tail"]["dp_moments_mol"])
    tail = observables_from_moments(mu, obs["repeat_unit_molar_mass_g_mol"], obs["end_group_molar_mass_g_mol"])
    return {
        "mu0_mol": mu[0], "mu1_mol": mu[1], "mu2_mol": mu[2],
        "end_inventory_mol": law["end_inventory_mol"],
        "explicit_dp_mol": dict(t3["explicit"]["dp_amount_mol"]),
        "polymer_condensed_mass_g": law["value_g"],
        "Mn_g_mol": tail.mn_g_mol, "Mw_g_mol": tail.mw_g_mol,
        "polymer_conversion": 1.0 - law["value_g"] / obs["initial_polymer_mass_g"],
        "species_mol": {k: v["amount_mol"] for k, v in t1["species"].items()},
        "polymer_derived_molecules_mol": t1["polymer_derived"]["mom_definition"]["amount_mol"],
    }


_TIER1 = ("polymer_condensed_mass_g", "polymer_derived_mass_exact_g", "polymer_derived_molecules_mol",
          "polymer_conversion")
_ZERO_WHEN_ABSENT = ("explicit_dp_mol.", "species_mol.")
EQUAL_RELATIVE = 1.0e-12


def _flatten(quantities: Mapping) -> Dict[str, Optional[float]]:
    out = {}
    for key, value in quantities.items():
        if isinstance(value, Mapping):
            for sub, v in value.items():
                out[f"{key}.{sub}"] = None if v is None else float(v)
        else:
            out[key] = None if value is None else float(value)
    return out


def side_by_side(obs: Mapping, mom: Optional[Mapping] = None) -> List[dict]:
    """Rows comparing one kMC observation with a MoM run's own numbers
    (``mom``: a dict keyed like ``mom_reported_quantities``; nested
    ``explicit_dp_mol`` / ``species_mol`` mappings).

    One row per MoM quantity -- every key of ``mom_reported_quantities``,
    every declared species, every key ``mom`` supplies (a DP or species the
    kMC run never produced has kMC value 0) -- plus
    ``polymer_derived_mass_exact_g``: the exact kMC polymer-derived mass
    against the MoM condensed mass, with the mass-law gap bound.  Every
    ``difference`` (kMC - MoM) and ``agreement`` is computed from the
    supplied MoM value: ``equal`` (relative difference <= 1e-12), ``within_bound``
    (|difference| <= the stated bound), ``differs``, ``missing_in_mom``,
    ``missing_in_kmc``; without ``mom`` every row is ``no_mom_values`` and no
    difference is claimed."""
    kmc = _flatten(mom_reported_quantities(obs))
    for name in obs.get("declared_species", ()):
        kmc.setdefault(f"species_mol.{name}", 0.0)
    kmc["polymer_derived_mass_exact_g"] = float(obs["tier1"]["polymer_derived"]["mass_g"])
    bounds = {"polymer_derived_mass_exact_g": float(obs["mom_mass_law"]["bound_g"])}
    given = _flatten(mom) if mom is not None else {}
    if mom is not None and "polymer_condensed_mass_g" in given:
        given["polymer_derived_mass_exact_g"] = given["polymer_condensed_mass_g"]
    keys = list(kmc) + [k for k in given if k not in kmc]
    rows = []
    for key in keys:
        k_value = kmc.get(key)
        if k_value is None and key.startswith(_ZERO_WHEN_ABSENT):
            k_value = 0.0
        m_value = given.get(key)
        bound = bounds.get(key)
        difference = relative = None
        if mom is None:
            agreement = "no_mom_values"
        elif m_value is None:
            agreement = "missing_in_mom"
        elif k_value is None:
            agreement = "missing_in_kmc"
        else:
            difference = k_value - m_value
            relative = _rel(k_value, m_value)
            if relative <= EQUAL_RELATIVE:
                agreement = "equal"
            elif bound is not None and abs(difference) <= bound:
                agreement = "within_bound"
            else:
                agreement = "differs"
        tier = 1 if key in _TIER1 or key.startswith("species_mol.") else 3
        rows.append({"tier": tier, "quantity": key, "kmc": k_value, "mom": m_value, "difference": difference,
                     "relative_difference": relative, "bound": bound, "agreement": agreement})
    return rows


__all__ = [
    "CONVENTION_DEFINITION", "ComponentRecord", "InvariantChecker", "InvariantError", "MASS_LAW_DEFINITION",
    "MeltObserver", "OBSERVABLES_SCHEMA", "ObservablesError", "SPECIES_RULE", "TIER1_NOTE", "Totals",
    "backbone_atoms_per_unit", "canonical_key", "declared_species_names", "exact_totals",
    "mom_reported_quantities", "molar_mass_g_mol", "observe", "side_by_side", "write_time_series",
]
