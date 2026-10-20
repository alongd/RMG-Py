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
by the same sum of absolute values (``bound_g``).  The law carries no
side-group depletion: the solver's depletion inventory Z is a model variable
(stationary tail mu1, tail only) that a kMC state does not determine.

Side-group loss (``side_group_channels``)
-----------------------------------------
Compared through physical observables only, per channel of the MoM deck's
``side_group_homolysis`` declaration (validated with the solver's rules; no
declaration, no channel): the cumulative amount of the channel's gas product
released, counted from events (``SideGroupReleaseTracker``: an X atom that
becomes a lone X radical, attributed to the site class it held on the current
graph), and the remaining inventory of X at the channel's site class in
polymer-derived molecules, by current-graph site matching, split tail /
explicit.  Nothing is reconstructed from provenance.  ``side_by_side``
compares the released amounts and the condensed mass (the exact kMC
polymer-derived mass) directly with the MoM run's.  Losses outside every
declared channel are ``mom_mass_law.undeclared_loss``.

Units: ``amount = scale N / N_A`` (mol) with ``scale`` the declaration's
``initialization().scale_factor``; the volume is constant
(``volume_m3 = scale x`` the representative volume).

Linear statistical copolymers (``PolymerSpec``)
-----------------------------------------------
No one repeat unit exists, so ``n_c`` is mass-based, ``floor(M_c / Mbar0)``
with ``Mbar0 = sum x_i M_i`` the mean residue molar mass (the MET length,
``melt_engine.TransportLaw``), and ``a`` in the reported fields is
``Mbar0``.  The scale and volume are those of the sampled melt
(``melt_sequence.initialization``).  The MoM-convention blocks (tier 1
``mom_definition``, ``mom_mass_law``, tier 2 ``mom_convention`` and
``end_mass_convention``, tier 3) and the side-group channels compare with
a MoM run, which reads v1 declarations only: they are None, and
``mom_reported_quantities`` and side-group channels are refused
(``refused:schema_v2_mom``).  ``composition`` reports, per residue, the
declared mole fraction, the realized one of the sampled melt (from the
atoms' origin, constant in time) and that of the polymer-derived molecules
now: a residue's atoms (its declared atom keys, by origin) in polymer-derived
molecules over its atom count, normalized over the residues.

Exact invariants (``InvariantChecker``): after every event, ``sum comp`` and
``sum B`` are constant; the molecule count changes by ``|products| -
|reactants|``; a bridge scission adds ``-2 M_a M_b`` to ``sum M^2`` and a
recombination ``+2 M_a M_b``; components outside the event are unchanged.
Given the other checks the ``sum M^2`` increment of the index's stored masses
is implied algebraically, so the expected increment is computed from masses
recomputed from each component's composition (exact rationals): the line then
catches stored-mass drift between the event's components.

Export: ``write_time_series`` writes one JSON line per requested time (the
engine's horizon API) and, last, one typed stop record (``"record": "stop"``,
schema ``STOP_SCHEMA``): why the call stopped (``STOP_REASONS``) and its event
count; observation lines carry no ``record`` key.  ``mom_reported_quantities`` maps one observation to the
quantities the MoM solver reports, ``side_by_side`` compares them with a MoM
run's own numbers.
"""

from __future__ import annotations

import json
import math
import weakref
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, Iterable, List, Mapping, Optional, Tuple

from rmgpy import constants
from rmgpy.kmc.melt_engine import (BUILTIN_SPECIES, MeltEngine, SpeciesIdentity, declared_initialization,
                                   monomer_product_species, refuse_initial_species_overlap)
from rmgpy.kmc.melt_state import _Builder, _backbone_keys, _element_mass, segment_residue
from rmgpy.polymer_mass import (Composition, EvolvedPoolMass, ExplicitInventory, InventoryTerm,
                                observables_from_moments)
from rmgpy.polymer_spec import PolymerSpec

OBSERVABLES_SCHEMA = "kmc-melt-observables/6"
STOP_SCHEMA = "kmc-melt-stop/1"
# The reasons ``write_time_series`` can stop with: the last requested time was
# reached ("horizon"), or the engine had no enabled channel when it was
# ("absorbing"), or the call's event cap was reached first ("max_events").
STOP_REASONS = ("horizon", "absorbing", "max_events")
UNITS = {"amount": "mol (scale N / N_A)", "mass": "g", "molar_mass": "g/mol", "M2": "g^2/mol",
         "dp_moments": "mol, mol, mol (sum of amount n^k, k = 0, 1, 2)", "time": "s", "volume": "m^3",
         "end_inventory": "mol of end-group pairs (head + tail)"}
CONVENTION_DEFINITION = (
    "delta_c = M_c - (a n_c + e) = molar mass of comp_c - n_c repeat - ends (signed); "
    "M1_difference = sum amount delta_c, M1_bound = sum amount |delta_c|; "
    "M2_difference = sum amount delta_c (2 (a n_c + e) + delta_c), "
    "M2_bound = sum amount |delta_c| (2 (a n_c + e) + |delta_c|); exact differences, not tolerances")
UNDECLARED_LOSS_DEFINITION = (
    "chain-origin atoms (repeat units and end groups) that are in no polymer-derived molecule, except the X "
    "atoms the release tracker attributed to a declared channel's release from a polymer-derived molecule "
    "(wherever they are now outside the polymer: free, or in a gas species formed later); a release from a "
    "site with no declared class stays in this loss; without a continuous tracker nothing is exempted "
    "(release_attribution_tracked = False); reported with mass, never folded into the MoM projection")
MASS_LAW_DEFINITION = (
    "value_g = a mu1_tail + e E + sum_n amount_n (a n + e): MoM's mass law WITHOUT side-group depletion "
    "(the solver's depletion inventory Z is a model variable that a kMC state does not determine; a loss of "
    "declared X atoms is visible through side_group_channels' released amounts and the exact condensed "
    "mass); E = physical end-group inventory of the tail (mass of end-group-origin atoms / e); "
    "difference_g = kMC polymer-derived mass - value_g = amount x [sum_tail mass(comp_c - n_c repeat - "
    "caps_c) + sum_explicit mass(comp_c - n_c repeat - ends)]; bound_g = the same with absolute values; a "
    "convention gap, never a tolerance for comparing with a MoM run")
SIDE_GROUP_DEFINITION = (
    "per declared side_group_homolysis channel, from a SideGroupReleaseTracker that follows the run (events "
    "or rescans, never provenance): released_count = GROSS count of releases from polymer-derived molecules "
    "(an X atom that becomes a free gas product, exactly the declared one: element, isotope, charge, one "
    "radical electron, lone pairs), attributed to the site class the X held, on the current graph, when last "
    "seen attached (RMG aromaticity perception, resonance-invariant); an X that re-enters a polymer-derived "
    "molecule after a release is a recapture (recaptured_count), released_net_count = released - recaptured; "
    "MoM's cumulative release Z(0) - Z is a gross channel flux (Z is never replenished); "
    "released_from_other_count = releases from non-polymer molecules (mobile or small species, light "
    "fragments) at the channel's site class, never booked on the polymer channel; released values are None "
    "unless the tracker can prove it saw every event (released_tracked); free_now = free gas-product atoms "
    "now; inventory = X now at the channel's site class in polymer-derived molecules, tail (n_c > xs) / "
    "explicit (n_c <= xs)")
BALANCE_DEFINITION = (
    "per channel element: initial_count (at tracker start; without a tracker the current total) = "
    "condensed_at_declared_sites (in polymer-derived molecules at a declared channel's site class) + "
    "condensed_elsewhere (in polymer-derived molecules otherwise) + free (free gas-product atoms) + other "
    "(every other molecule: mobile or small species, gas species such as HX, light fragments)")
TIER1_NOTE = (
    "mom_definition is an accounting identity on the kMC partition (the MoM mass law and molecule count "
    "rebuilt from the kMC components), not a comparison with a MoM run; compare MoM output with side_by_side")
SPECIES_RULE = (
    "a component isomorphic under resonance (all H explicit) to a declared species (small species, "
    "monomer_product, built-in H and H2) is that declared species whatever its ancestry; any other component "
    "is polymer-derived iff n_c >= 1, else an undeclared light product; a declared species isomorphic to an "
    "initial chain, or a small-species label equal to a built-in label, is refused at run start")
END_SEGMENTS = ("head_end", "tail_end")
MOM_V2_REFUSAL = ("refused:schema_v2_mom: the method-of-moments comparison reads v1 declarations only; a "
                  "linear copolymer (schema v2) has no MoM run to compare with yet")


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


def backbone_atoms_per_unit(spec) -> Optional[int]:
    """``b``: atoms on the declared repeat unit's head-port -> tail-port path
    (None for a copolymer, whose length is mass-based)."""
    if isinstance(spec, PolymerSpec):
        return None
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
    n: int                    # floor(B_c / b); a copolymer: floor(M_c / Mbar0)
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


SIDE_GROUP_CHANNEL_KEYS = ("label", "A", "n", "Ea", "site_selector", "sites_per_unit", "gas_product")
SIDE_GROUP_SITE_SELECTORS = ("aryl", "benzylic", "aliphatic")


def _sanitize(label: str) -> str:
    return "".join(c if (c.isalnum() or c == "_") else "_" for c in label)


def side_group_site_indices(monomer, element_symbol, site_selector) -> Tuple[int, ...]:
    """Python mirror of the MoM solver's ``side_group_site_atom_indices``
    (``rmgpy/solver/polymer.pyx``): indices in ``monomer.atoms`` of the
    removable X atoms (element match, unlabeled, closed-shell, terminal,
    single-bonded, on a carbon) whose carbon environment is the selector's:
    ``aryl`` (the carbon is in an aromatic ring), ``benzylic`` (not aromatic,
    bonded to an aromatic carbon), ``aliphatic`` (neither)."""
    aromatic_atoms = set()
    for ring in monomer.get_aromatic_rings()[0]:
        aromatic_atoms.update(ring)
    matches = []
    for idx, atom in enumerate(monomer.atoms):
        if (atom.symbol != element_symbol or atom.label
                or atom.radical_electrons != 0 or len(atom.bonds) != 1
                or not next(iter(atom.bonds.values())).is_single()):
            continue
        neighbor = next(iter(atom.bonds.keys()))
        if not neighbor.is_carbon():
            continue
        if neighbor in aromatic_atoms:
            cls = "aryl"
        elif any(nb in aromatic_atoms and nb.is_carbon() for nb in neighbor.bonds.keys()):
            cls = "benzylic"
        else:
            cls = "aliphatic"
        if cls == site_selector:
            matches.append(idx)
    return tuple(matches)


def validate_side_group_channels(spec, channels) -> List[dict]:
    """Validate a ``side_group_homolysis`` channel list (the polymer() deck
    structure the MoM solver reads) against the declaration, with the
    solver's rules (``validate_side_group_homolysis`` mirrored in Python:
    exact keys, label, selector vocabulary, A > 0, Ea >= 0, sites_per_unit
    > 0, a monoatomic mono-radical gas_product, unique sanitized labels;
    structurally against the repeat unit's legacy molecule: the selector
    matches >= 1 site, sites_per_unit equals the match count, no two
    channels on the same atom set; plus the shared-spec isotope check).
    Returns one normalized row per channel: label, site_selector,
    sites_per_unit (int), gas_product, element, isotope, molar_mass_g_mol
    (the solver's M_X), site_keys (repeat-unit origin keys of the X atoms) and
    neighbor (site key -> declared neighbor key)."""
    from rmgpy.molecule import Molecule
    if isinstance(channels, Mapping) or not isinstance(channels, (list, tuple)):
        raise ObservablesError("side_group_homolysis must be a list of channel dicts")
    if not channels:
        raise ObservablesError("side_group_homolysis must not be an empty list (omit it to declare no channel)")
    rows, seen = [], {}
    for pos, ch in enumerate(channels):
        if not isinstance(ch, Mapping):
            raise ObservablesError(f"side_group_homolysis entry {pos} must be a channel dict")
        missing = [k for k in SIDE_GROUP_CHANNEL_KEYS if k not in ch]
        extra = sorted(set(ch) - set(SIDE_GROUP_CHANNEL_KEYS))
        if missing or extra:
            raise ObservablesError(f"side_group_homolysis entry {pos} must have exactly the keys "
                                   f"{SIDE_GROUP_CHANNEL_KEYS}; missing {missing}, unknown {extra}")
        label = ch["label"]
        if not isinstance(label, str) or not label.strip():
            raise ObservablesError(f"side_group_homolysis entry {pos} label must be a non-empty string")
        selector = ch["site_selector"]
        if selector not in SIDE_GROUP_SITE_SELECTORS:
            raise ObservablesError(f"channel {label!r} site_selector={selector!r} must be one of "
                                   f"{SIDE_GROUP_SITE_SELECTORS}")
        values = {}
        for key in ("A", "n", "Ea", "sites_per_unit"):
            v = ch[key]
            if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(float(v)):
                raise ObservablesError(f"channel {label!r} {key}={v!r} must be a finite number")
            values[key] = float(v)
        if values["A"] <= 0.0 or values["Ea"] < 0.0 or values["sites_per_unit"] <= 0.0:
            raise ObservablesError(f"channel {label!r}: A > 0, Ea >= 0 and sites_per_unit > 0 are required")
        gas = ch["gas_product"]
        if not isinstance(gas, str) or not gas.strip():
            raise ObservablesError(f"channel {label!r} gas_product must be a SMILES string")
        try:
            gmol = Molecule().from_smiles(gas)
        except Exception as error:
            raise ObservablesError(f"channel {label!r} gas_product={gas!r} does not parse ({error})") from None
        if gmol.get_radical_count() != 1:
            raise ObservablesError(f"channel {label!r} gas_product={gas!r} must be a mono-radical")
        if len(gmol.atoms) != 1:
            raise ObservablesError(f"channel {label!r} gas_product={gas!r} has {len(gmol.atoms)} atoms; only "
                                   "monoatomic radical gas products can be represented (the solver's v1 rule)")
        sanitized = _sanitize(label)
        if sanitized in seen:
            raise ObservablesError(f"side_group_homolysis has duplicate channel labels {seen[sanitized]!r} "
                                   f"and {label!r}")
        seen[sanitized] = label
        gas_atom = gmol.atoms[0]
        isotope = None if gas_atom.element.isotope == -1 else gas_atom.element.isotope
        rows.append({"label": label, "site_selector": selector, "sites_per_unit": values["sites_per_unit"],
                     "gas_product": gas, "element": gas_atom.symbol, "isotope": isotope,
                     "signature": (gas_atom.symbol, isotope, int(gas_atom.charge),
                                   int(gas_atom.radical_electrons), int(gas_atom.lone_pairs)),
                     "molar_mass_g_mol": gmol.get_molecular_weight() * 1000.0})
    # structural layer, on the legacy repeat-unit molecule the solver validates against
    repeat_molecule, _ends = spec.legacy_molecules()
    sites_seen = {}
    for row in rows:
        idxs = side_group_site_indices(repeat_molecule, row["element"], row["site_selector"])
        if not idxs:
            raise ObservablesError(f"channel {row['label']!r} site_selector={row['site_selector']!r} matches NO "
                                   f"removable side-group {row['element']} atom in the repeat unit")
        if float(len(idxs)) != row["sites_per_unit"]:
            raise ObservablesError(f"channel {row['label']!r} sites_per_unit={row['sites_per_unit']:g} contradicts "
                                   f"the structural match count {len(idxs)}")
        isotopes = {repeat_molecule.atoms[i].element.isotope for i in idxs}
        if isotopes != {-1 if row["isotope"] is None else row["isotope"]}:
            raise ObservablesError(f"channel {row['label']!r} source isotope(s) {sorted(isotopes)} disagree with "
                                   "its gas_product")
        key = (row["element"], frozenset(idxs))
        if key in sites_seen:
            raise ObservablesError(f"channels {sites_seen[key]!r} and {row['label']!r} resolve to the SAME "
                                   f"{row['element']} atom set")
        sites_seen[key] = row["label"]
        _symbol, _iso, charge, radical, lone_pairs = row["signature"]
        for i in idxs:
            site = repeat_molecule.atoms[i]
            if (int(site.charge), int(site.radical_electrons) + 1, int(site.lone_pairs)) != (charge, radical, lone_pairs):
                raise ObservablesError(
                    f"channel {row['label']!r} gas_product={row['gas_product']!r} (charge {charge}, {radical} "
                    f"radical electron(s), {lone_pairs} lone pairs) cannot be matched: homolysis of its site "
                    f"gives charge {site.charge}, {site.radical_electrons + 1} radical electron(s), "
                    f"{site.lone_pairs} lone pairs")
        row["sites_per_unit"] = len(idxs)
    return rows


PERCEPTION_STATS = {"calls": 0, "atoms": 0, "ring_system_hits": 0, "whole_fallbacks": 0}
_RING_SYSTEM_CACHE: Dict[Tuple, frozenset] = {}
_COMPONENT_MEMO = weakref.WeakKeyDictionary()   # state -> (revision, {cid: atoms})
_CAPS = {2.0: ("C", 2), 3.0: ("C", 1)}          # a capped multiple bond keeps its order


def clear_aromaticity_cache():
    _RING_SYSTEM_CACHE.clear()
    _COMPONENT_MEMO.clear()


def _perceive(rows) -> frozenset:
    """Atom keys of ``rows`` on RMG-aromatic rings (``get_aromatic_rings``)."""
    from rmgpy.kmc.rmg_moves import molecule_from_rows
    PERCEPTION_STATS["calls"] += 1
    PERCEPTION_STATS["atoms"] += len(rows)
    rings, _bonds = molecule_from_rows(rows).get_aromatic_rings()
    return frozenset(atom.id for ring in rings for atom in ring)


def _ring_atoms_and_bridges(state, atoms):
    """(atoms on a cycle, bridge edges) of one component (iterative Tarjan)."""
    nodes = state.atoms
    disc, low, bridges, counter = {}, {}, set(), 0
    for root in atoms:
        if root in disc:
            continue
        disc[root] = low[root] = counter
        counter += 1
        stack = [(root, -1, iter(sorted(nodes[root].edges)))]
        while stack:
            v, parent, it = stack[-1]
            w = next(it, None)
            if w is None:
                stack.pop()
                if stack:
                    u = stack[-1][0]
                    low[u] = min(low[u], low[v])
                    if low[v] > disc[u]:
                        bridges.add((min(u, v), max(u, v)))
                continue
            if w == parent:
                continue
            if w in disc:
                low[v] = min(low[v], disc[w])
            else:
                disc[w] = low[w] = counter
                counter += 1
                stack.append((w, v, iter(sorted(nodes[w].edges))))
    ring = {a for a in atoms if any((min(a, b), max(a, b)) not in bridges for b in nodes[a].edges)}
    return ring, bridges


def _ring_system_rows(state, group):
    """Rows of one ring system (connected ring atoms) and its first shell
    (direct neighbours), each first-shell bond leaving the system replaced by
    a cap of the same order (H for single, CH2= for double, CH# for triple):
    the atoms whose element, charge, electrons and bonds RDKit's aromaticity
    model reads.  Ids are ranks (sorted atom id); returns (rows, rank -> atom)
    or None when a leaving bond cannot be capped (falls back to the whole
    component).

    The second shell never decides the answer: RDKit counts a ring atom's
    pi electrons from that atom's own bonds, charge and electrons and the
    element of an exocyclic double-bond partner (first shell; C=O versus
    C=CH2 can flip it), and RMG keeps only six-membered carbon rings, all
    inside the ring system.  The cap's order matters only for validity: it
    keeps every first-shell atom's bond-order sum, so RMG's atom typing
    (``Molecule.update``) and the RDKit conversion accept the graph."""
    nodes = state.atoms
    shell = sorted({b for a in group for b in nodes[a].edges} - set(group))
    included = sorted(group) + shell
    rank = {a: k for k, a in enumerate(sorted(included))}
    rows, extra = {}, len(rank)
    for a in sorted(included):
        n = nodes[a]
        edges = []
        for b, order in n.edges.items():
            if b in rank:
                edges.append((rank[b], float(order)))
                continue
            order = float(order)
            if order == 1.0:
                rows[extra] = ("H", None, 0, 0, 0, ((rank[a], 1.0),))
                edges.append((extra, 1.0))
                extra += 1
            elif order in _CAPS:
                element, hydrogens = _CAPS[order]
                cap = extra
                cap_edges = [(rank[a], order)]
                extra += 1
                for _h in range(hydrogens):
                    rows[extra] = ("H", None, 0, 0, 0, ((cap, 1.0),))
                    cap_edges.append((extra, 1.0))
                    extra += 1
                rows[cap] = (element, None, 0, 0, 0, tuple(sorted(cap_edges)))
                edges.append((cap, order))
            else:
                return None
        rows[rank[a]] = (n.element, n.isotope, n.charge, n.radical, n.lone_pairs, tuple(sorted(edges)))
    return rows, {k: a for a, k in rank.items()}


def aromatic_atoms(state, cid) -> frozenset:
    """Atoms of component ``cid`` on an aromatic ring by RMG's perception
    (``Molecule.get_aromatic_rings``: six-membered carbon rings, RDKit
    aromaticity, independent of the Kekule form stored), the perception the
    solver's site selector uses.  Perceived per ring system (connected ring
    atoms plus their capped first shell; aromaticity is a property of the
    ring system and its direct substituents), cached by the ring system's
    labelled graph across components and events, and memoized per state
    revision; equal to ``aromatic_atoms_whole``."""
    memo = None
    if getattr(state, "_tx", None) is None:            # no open transaction: revision pins the graph
        held = _COMPONENT_MEMO.get(state)
        if held is None or held[0] != state.revision:
            held = _COMPONENT_MEMO[state] = (state.revision, {})
        memo = held[1]
        if cid in memo:
            return memo[cid]
    atoms = state.components.atoms(cid)
    ring, _bridges = _ring_atoms_and_bridges(state, atoms)
    out = set()
    seen = set()
    for start in sorted(ring):
        if start in seen:
            continue
        group, stack = [], [start]
        seen.add(start)
        while stack:
            a = stack.pop()
            group.append(a)
            for b in state.atoms[a].edges:
                if b in ring and b not in seen:
                    seen.add(b)
                    stack.append(b)
        if len(group) < 6:
            continue
        built = _ring_system_rows(state, group)
        if built is None:
            PERCEPTION_STATS["whole_fallbacks"] += 1
            out = set(aromatic_atoms_whole(state, cid))
            break
        rows, back = built
        graph_key = tuple(sorted(rows.items()))
        found = _RING_SYSTEM_CACHE.get(graph_key)
        if found is None:
            found = _RING_SYSTEM_CACHE[graph_key] = _perceive(rows)
            if len(_RING_SYSTEM_CACHE) > 20000:
                _RING_SYSTEM_CACHE.clear()
        else:
            PERCEPTION_STATS["ring_system_hits"] += 1
        out.update(back[k] for k in found if k in back)
    result = frozenset(out)
    if memo is not None:
        memo[cid] = result
    return result


def aromatic_atoms_whole(state, cid) -> frozenset:
    """Whole-component RMG perception (the reference ``aromatic_atoms`` must equal)."""
    rows = {}
    for a in state.components.atoms(cid):
        n = state.atoms[a]
        rows[a] = (n.element, n.isotope, n.charge, n.radical, n.lone_pairs,
                   tuple(sorted((b, float(o)) for b, o in n.edges.items())))
    return _perceive(rows)


def site_class(state, x, aromatic: Optional[frozenset] = None) -> Optional[str]:
    """The solver's site predicate on the current graph: ``x`` is a terminal,
    single-bonded, closed-shell atom on a carbon; its class is ``aryl`` (the
    carbon is on an aromatic ring), ``benzylic`` (bonded to an aromatic-ring
    carbon) or ``aliphatic``; None otherwise.  ``aromatic``: the component's
    ``aromatic_atoms`` (computed when not given)."""
    node = state.atoms[x]
    if node.radical or len(node.edges) != 1:
        return None
    (c, order), = node.edges.items()
    if order != 1.0 or state.atoms[c].element != "C":
        return None
    if aromatic is None:
        aromatic = aromatic_atoms(state, state.components.component_of(x))
    if c in aromatic:
        return "aryl"
    if any(b in aromatic and state.atoms[b].element == "C" for b in state.atoms[c].edges if b != x):
        return "benzylic"
    return "aliphatic"


def declared_species_names(spec, identity: SpeciesIdentity) -> List[str]:
    """Every name the identity can return: small species, ``monomer_product``
    (when it is not also a small species) and the built-ins."""
    names = [name for name, _g, _d in spec.small_species]
    for label, monomer in monomer_product_species(spec):
        builder = _Builder()
        builder.add_segment(monomer, "mobile", 0, label)
        alone = builder.state()
        if identity.identify(alone, alone.components.ids()[0]) == label:
            names.append(label)
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

    def __init__(self, spec, *, xs: int, identity: Optional[SpeciesIdentity] = None,
                 side_group_homolysis=None):
        if isinstance(xs, bool) or not isinstance(xs, int) or xs < 0:
            raise ObservablesError(f"xs must be a non-negative integer, got {xs!r}")
        self.spec = spec
        self.xs = xs
        self.b = backbone_atoms_per_unit(spec)
        self.sequence = isinstance(spec, PolymerSpec)
        if self.sequence:
            from rmgpy.kmc import melt_sequence as msq
            if side_group_homolysis is not None:
                raise ObservablesError(f"{MOM_V2_REFUSAL} (side_group_homolysis channels)")
            self.mean_residue = msq.mean_residue_molar_mass(spec)
            self.repeat = None
            ends = spec.head_end_group.composition.plus(spec.tail_end_group.composition)
            self.ends = _from_declared(ends)
            self.a = float(self.mean_residue)
            self.e = ends.molar_mass_g_mol()
            self.residue_keys = {name: frozenset(atom.key for atom in residue.graph.atoms)
                                 for name, residue in spec.residues}
        else:
            self.repeat = _from_declared(spec.repeat_unit.composition)
            self.ends = _from_declared(spec.end_composition)
            self.a = spec.repeat_unit.composition.molar_mass_g_mol()
            self.e = spec.end_composition.molar_mass_g_mol()
        init = declared_initialization(spec)
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
        # side-group-homolysis channels: exactly the MoM declaration's (none: no depletion rows)
        self.channels = ([] if side_group_homolysis is None
                         else validate_side_group_channels(spec, side_group_homolysis))
        self.channel_by_site = {(row["element"], row["isotope"], row["site_selector"]): i
                                for i, row in enumerate(self.channels)}
        self.channel_elements = {(row["element"], row["isotope"]) for row in self.channels}
        self.gas_signatures = {row["signature"] for row in self.channels}
        self.declaration = frozenset(channel_identity(row) for row in self.channels)

    def category(self, state, cid) -> str:
        """"species", "polymer" or "light": the module's classification of one component."""
        if self.identity.identify(state, cid) is not None:
            return "species"
        return "polymer" if self.length(state, cid) >= 1 else "light"

    def length(self, state, cid, backbone: Optional[int] = None) -> int:
        """``n_c``: ``floor(B_c / b)``, or ``floor(M_c / Mbar0)`` for a copolymer."""
        comps = state.components
        if self.sequence:
            return math.floor(comps.molar_mass_g_mol(cid) / self.mean_residue)
        return (comps.backbone_count(cid) if backbone is None else backbone) // self.b

    def channel_of(self, state, x, aromatic: Optional[frozenset] = None) -> Optional[int]:
        """The declared channel whose X atom ``x`` is, by current-graph site matching."""
        node = state.atoms[x]
        if (node.element, node.isotope) not in self.channel_elements:
            return None
        cls = site_class(state, x, aromatic)
        return None if cls is None else self.channel_by_site.get((node.element, node.isotope, cls))

    def is_free_gas_product(self, state, x) -> bool:
        """``x`` is a lone atom exactly matching a declared gas product
        (element, isotope, charge, radical electrons, lone pairs)."""
        node = state.atoms[x]
        return (not node.edges and (node.element, node.isotope, node.charge, node.radical, node.lone_pairs)
                in self.gas_signatures)

    def _x_atoms(self, state, cid) -> List[int]:
        return [x for x in state.components.atoms(cid)
                if (state.atoms[x].element, state.atoms[x].isotope) in self.channel_elements]

    def site_inventory(self, state, records=None) -> Dict[int, Dict[str, int]]:
        """Per channel index: X atoms at the channel's site class in polymer-derived
        molecules now, {"tail": n, "explicit": n} (representative counts)."""
        out = {i: {"tail": 0, "explicit": 0} for i in range(len(self.channels))}
        if not self.channels:
            return out
        records = self.components(state) if records is None else records
        for r in records:
            if r.category != "polymer":
                continue
            xs_here = self._x_atoms(state, r.cid)
            if not xs_here:
                continue
            aromatic = aromatic_atoms(state, r.cid)
            part = "tail" if r.n > self.xs else "explicit"
            for x in xs_here:
                i = self.channel_of(state, x, aromatic)
                if i is not None:
                    out[i][part] += 1
        return out

    def x_balance(self, state, tracker: Optional["SideGroupReleaseTracker"] = None, records=None) -> dict:
        """Per channel element (representative counts): initial = condensed at
        declared sites + condensed elsewhere + free + other (``BALANCE_DEFINITION``)."""
        out = {}
        for element in sorted(self.channel_elements, key=str):
            out[element] = {"condensed_at_declared_sites_count": 0, "condensed_elsewhere_count": 0,
                            "free_count": 0, "other_count": 0}
        if not self.channels:
            return {}
        records = {r.cid: r for r in (self.components(state) if records is None else records)}
        for cid, r in records.items():
            xs_here = self._x_atoms(state, cid)
            if not xs_here:
                continue
            aromatic = aromatic_atoms(state, cid) if r.category == "polymer" else None
            for x in xs_here:
                slot = out[(state.atoms[x].element, state.atoms[x].isotope)]
                if r.category == "polymer":
                    at_site = self.channel_of(state, x, aromatic) is not None
                    slot["condensed_at_declared_sites_count" if at_site else "condensed_elsewhere_count"] += 1
                elif self.is_free_gas_product(state, x):
                    slot["free_count"] += 1
                else:
                    slot["other_count"] += 1
        for element, slot in out.items():
            total = sum(slot.values())
            slot["initial_count"] = (tracker.initial[element] if tracker is not None and tracker.continuous
                                     else total)
        return {f"{el}{iso or ''}": slot for (el, iso), slot in out.items()}

    def release_tracker(self, target) -> "SideGroupReleaseTracker":
        """The tracker that follows ``target``: for an engine, the run's one
        tracker for this channel declaration (registered, shared with
        ``write_time_series``); for a bare state, a new tracker updated by
        hand."""
        if isinstance(target, MeltEngine):
            return _run_tracker(target, self)
        return SideGroupReleaseTracker(self, target)

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
            n = self.length(state, cid, backbone)
            caps = Counter()
            for x in comps.atoms(cid):
                o = origin[x]
                if o.kind == "chain" and o.segment in END_SEGMENTS:
                    caps[(nodes[x].element, nodes[x].isotope)] += 1
            caps = _key(caps)
            species = engine.entries[cid].species if synced else self.identity.identify(state, cid)
            if species is not None:
                category, defect, law_defect = "species", (), ()
            elif n >= 1 and self.sequence:
                category, defect, law_defect = "polymer", (), ()
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
    def observe(self, target, releases: Optional["SideGroupReleaseTracker"] = None) -> dict:
        """One observation (a JSON-ready dict) of an engine or a bare state;
        ``releases``: the tracker that followed the run (released amounts are
        None without one)."""
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
        mom_mass = mom_amount = None
        if not self.sequence:
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
                             "mom_definition": None if self.sequence else {
                                 "amount_mol": amount, "mass_g": amount * self.declared_molar_mass[name]}}
        polymer_mass = u * _fsum(r.molar_mass_g_mol for r in polymer)
        tier1 = {
            "mom_definition_note": TIER1_NOTE,
            "polymer_derived": {"molecules": len(polymer), "amount_mol": u * len(polymer),
                                "mass_g": polymer_mass,
                                "mom_definition": (None if self.sequence else
                                                   {"amount_mol": mom_amount, "mass_g": mom_mass})},
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
        chain_defect = 0.0
        explicit_capped = _fsum(u * k * (a * n + e) for n, k in ladder.items())
        value = _fsum([a * mol_tail[1], e * end_mol, -mol_tail[0] * chain_defect, explicit_capped])
        # exact gap: per tail component comp - n repeat - caps, per explicit chain comp - n repeat - ends
        gaps = ([molar_mass_g_mol(r.law_defect) for r in tail]
                + [molar_mass_g_mol(r.defect) for r in explicit])
        # undeclared loss: chain-origin atoms outside every polymer-derived molecule, except free
        # gas-product atoms of declared channels (those are the channels' releases)
        if releases is not None and engine is not None:
            releases.sync(engine)
        tracked = (releases is not None and releases.continuous
                   and releases.declaration == self.declaration)
        exempt = releases.released_from if tracked else {}
        lost_atoms = Counter()
        for r in records:
            if r.category == "polymer":
                continue
            for x in state.components.atoms(r.cid):
                if state.origin[x].kind != "chain" or x in exempt:
                    continue
                lost_atoms[(state.atoms[x].element, state.atoms[x].isotope)] += 1
        lost_atoms = _key(lost_atoms)
        # side-group channels: released (from events) and current-graph inventory
        inventory = self.site_inventory(state, records)
        free_now = Counter()
        for x, node in state.atoms.items():
            if not node.edges and self.channels and self.is_free_gas_product(state, x):
                free_now[(node.element, node.isotope, node.charge, node.radical, node.lone_pairs)] += 1
        channels = []
        for i, row in enumerate(self.channels):
            m_x = row["molar_mass_g_mol"]

            ident = channel_identity(row)
            gross = releases.counts[ident] if tracked else None
            recaptured = releases.recaptured[ident] if tracked else None
            channels.append({"label": row["label"], "site_selector": row["site_selector"],
                             "sites_per_unit": row["sites_per_unit"], "gas_product": row["gas_product"],
                             "molar_mass_g_mol": m_x, "released_tracked": tracked,
                             "released_count": gross, "recaptured_count": recaptured,
                             "released_net_count": None if not tracked else gross - recaptured,
                             "released_from_other_count": releases.other_counts[ident] if tracked else None,
                             "released_mol": u * gross if tracked else None,
                             "released_mass_g": u * gross * m_x if tracked else None,
                             "released_net_mol": None if not tracked else u * (gross - recaptured),
                             "free_now_mol": u * free_now[row["signature"]],
                             "inventory_tail_mol": u * inventory[i]["tail"],
                             "inventory_explicit_mol": u * inventory[i]["explicit"]})
        balance = {}
        for element, slot in self.x_balance(state, releases if tracked else None, records).items():
            balance[element] = dict(slot, **{k[:-6] + "_mol": u * v for k, v in slot.items()})
        mass_law = {"definition": MASS_LAW_DEFINITION, "end_inventory_mol": end_mol,
                    "tail_repeat_mass_g": a * mol_tail[1], "tail_end_mass_g": e * end_mol,
                    "chain_mass_defect_g_mol": chain_defect,
                    "explicit_capped_mass_g": explicit_capped,
                    "value_g": value, "kmc_mass_g": polymer_mass,
                    "difference_g": u * _fsum(gaps),
                    "bound_g": u * _fsum(abs(g) for g in gaps),
                    "undeclared_loss": {"definition": UNDECLARED_LOSS_DEFINITION,
                                        "release_attribution_tracked": tracked,
                                        "atoms": {f"{el}{iso or ''}": k for (el, iso), k in lost_atoms},
                                        "mass_g": u * molar_mass_g_mol(lost_atoms)}}

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

        if self.sequence:       # the MoM-convention blocks (module doc)
            mass_law = tier3 = None
            tier2["mom_convention"] = tier2["end_mass_convention"] = None

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
        out = {
            "schema": OBSERVABLES_SCHEMA, "time_s": time_s, "engine_time_s": time_s,
            "temperature_K": temperature, "units": UNITS,
            "scale_factor": self.scale, "avogadro": constants.Na, "xs": self.xs,
            "backbone_atoms_per_unit": self.b, "repeat_unit_molar_mass_g_mol": a,
            "end_group_molar_mass_g_mol": e, "initial_polymer_mass_g": self.initial_polymer_mass_g,
            "representative_volume_m3": rep_volume, "volume_m3": self.scale * rep_volume,
            "declared_species": list(self.declared), "species_rule": SPECIES_RULE,
            "tier1": tier1, "mom_mass_law": mass_law, "tier2": tier2, "tier3": tier3, "release": release,
            "side_group_definition": SIDE_GROUP_DEFINITION, "side_group_channels": channels,
            "side_group_unattributed_releases": dict(releases.unattributed) if tracked else None,
            "side_group_balance_definition": BALANCE_DEFINITION, "side_group_balance": balance,
            "counts": counts, "leak_estimate": leak,
            "invariants": {"components": len(records), "backbone_atoms": sum(r.backbone for r in records),
                           "composition": {f"{el}{iso or ''}": k for (el, iso), k in _key(total_comp)}},
        }
        if self.sequence:
            out.update({"length_rule": "mean_residue_molar_mass",
                        "mean_residue_molar_mass_g_mol": self.a,
                        "composition": self.composition(state, records)})
        return out

    def composition(self, state, records=None) -> dict:
        """Per residue of a copolymer: declared, realized (initial, from the
        atoms' origin) and polymer-derived mole fractions, and the polymer-derived
        residue amounts (module doc)."""
        records = self.components(state) if records is None else records
        category = {r.cid: r.category for r in records}
        names = [name for name, _residue in self.spec.residues]
        units = {name: set() for name in names}
        atoms = {name: 0 for name in names}
        for i, origin in state.origin.items():
            if origin.kind != "chain":
                continue
            name = segment_residue(origin.segment)
            if name is None or origin.key not in self.residue_keys[name]:
                continue
            units[name].add((origin.ordinal, origin.segment))
            if category[state.components.component_of(i)] == "polymer":
                atoms[name] += 1
        size = {name: len(keys) for name, keys in self.residue_keys.items()}
        equivalents = {name: Fraction(atoms[name], size[name]) for name in names}

        def fractions(counts):
            total = sum(counts.values(), Fraction(0))
            return {name: float(Fraction(counts[name]) / total) if total else 0.0 for name in names}

        tolerance = self.spec.initial_ensemble.composition_tolerance
        return {"declared": {name: float(x) for name, x in self.spec.initial_ensemble.composition},
                "realized_initial": fractions({name: len(units[name]) for name in names}),
                "polymer_derived": fractions(equivalents),
                "polymer_derived_residue_mol": {name: self.unit * float(x) for name, x in equivalents.items()},
                "relative_tolerance": None if tolerance is None else float(tolerance)}


class SideGroupReleaseTracker:
    """Follows one run and counts side-group releases from what happens to
    the graph, never from provenance.

    An X atom (element of a declared channel) that becomes a free gas
    product (exactly the declared one) is one release, booked by where it was
    when last seen attached: from a polymer-derived molecule at a declared
    channel's site class -> that channel (gross); from a polymer-derived
    molecule elsewhere -> ``unattributed["polymer"]`` (its mass stays in the
    undeclared loss); from any other molecule -> the channel's
    ``other_counts`` or ``unattributed["other"]``.  An X released from a
    polymer through a channel that re-enters a polymer-derived molecule is a
    recapture of that channel.

    Feeding: ``engine.run(check=tracker)`` (one call per fired event), or
    ``update(state)`` after hand edits (a rescan of every component).
    Continuity: the tracker counts the engine's fired events; ``sync(engine)``
    (called by ``observe``) rescans when exactly one event was missed (one
    event cannot release and recapture the same atom, so the rescan is
    exact, e.g. after a post-commit index failure) and drops ``continuous``
    when more were missed or when the tracker started on an engine that had
    already fired events."""

    _UNSEEN = object()

    def __init__(self, observer: MeltObserver, target):
        engine = target if isinstance(target, MeltEngine) else None
        state = engine.state if engine is not None else target
        self.observer = observer
        self.idents = [channel_identity(row) for row in observer.channels]
        self.declaration = frozenset(self.idents)
        self.engine = weakref.ref(engine) if engine is not None else None
        self.events_seen = engine.events_fired if engine is not None else 0
        self.continuous = engine is None or engine.events_fired == 0
        self.counts = Counter({k: 0 for k in self.idents})          # keyed by channel identity
        self.recaptured = Counter({k: 0 for k in self.idents})
        self.other_counts = Counter({k: 0 for k in self.idents})
        self.unattributed = {"polymer": 0, "other": 0}
        self.attached: Dict[int, Tuple[Optional[Tuple], str]] = {}  # X -> (channel identity or None, source)
        self.released_from: Dict[int, Tuple] = {}                   # X released from a polymer via a channel
        self.updates = 0
        self.initial = {}
        for element in observer.channel_elements:
            self.initial[element] = sum(1 for node in state.atoms.values()
                                        if (node.element, node.isotope) == element)
        self.update(state)
        self.updates = 0

    def __call__(self, engine, event):
        """One fired event.  Idempotent in ``events_fired``: an event already
        caught up by ``sync`` is skipped, not treated as missed."""
        if self.engine is None or self.engine() is not engine:
            self.continuous = False
        new = engine.events_fired - self.events_seen
        if new <= 0:
            return
        if new == 1:
            self.update(engine.state, event)
        else:
            self.update(engine.state)
            self.continuous = False
        self.events_seen = engine.events_fired

    def sync(self, engine):
        """Catch up with events fired without the tracker (see the class doc)."""
        if self.engine is None or self.engine() is not engine:
            self.continuous = False
            return
        missed = engine.events_fired - self.events_seen
        if missed:
            self.update(engine.state)
            if missed > 1:
                self.continuous = False
            self.events_seen = engine.events_fired

    def update(self, state, event=None):
        obs = self.observer
        if not obs.channels:
            return
        comps = state.components
        cids = comps.ids() if event is None else [c for c in event.products if c in comps._rec]
        for cid in cids:
            xs_here = obs._x_atoms(state, cid)
            if not xs_here:
                continue
            category = aromatic = None
            for x in xs_here:
                if obs.is_free_gas_product(state, x):
                    previous = self.attached.pop(x, self._UNSEEN)
                    if previous is self._UNSEEN:
                        continue
                    channel, source = previous           # channel identity or None
                    if source == "polymer":
                        if channel is None:
                            self.unattributed["polymer"] += 1
                            self.released_from.pop(x, None)
                        else:
                            self.counts[channel] += 1
                            self.released_from[x] = channel
                    elif channel is None:
                        self.unattributed["other"] += 1
                    else:
                        self.other_counts[channel] += 1
                    continue
                if category is None:
                    category = obs.category(state, cid)
                    aromatic = aromatic_atoms(state, cid)
                source = "polymer" if category == "polymer" else "other"
                if source == "polymer" and x in self.released_from:
                    self.recaptured[self.released_from.pop(x)] += 1
                index = obs.channel_of(state, x, aromatic)
                self.attached[x] = (None if index is None else self.idents[index], source)
        self.updates += 1


class _Checks:
    def __init__(self, *checks):
        self.checks = [c for c in checks if c is not None]

    def __call__(self, engine, event):
        for check in self.checks:
            check(engine, event)


def observe(target, spec=None, *, xs: int, side_group_homolysis=None) -> dict:
    """One observation of an engine (or of a bare state with its ``spec``)."""
    if spec is None:
        if not isinstance(target, MeltEngine):
            raise ObservablesError("a bare state needs its declaration")
        spec = target.spec
    return MeltObserver(spec, xs=xs, side_group_homolysis=side_group_homolysis).observe(target)


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
_RUN_TRACKERS = weakref.WeakKeyDictionary()     # engine -> {channel declaration: tracker}


def channel_identity(row) -> Tuple:
    """A channel's identity: (label, site selector, gas product), never its position."""
    return (row["label"], row["site_selector"], row["gas_product"])


def _run_tracker(engine, observer) -> SideGroupReleaseTracker:
    """The release tracker of this engine run for these channels: created
    once (continuous only if the engine has fired no event yet) and reused by
    every later export, so released amounts are cumulative across appends."""
    key = observer.declaration
    trackers = _RUN_TRACKERS.setdefault(engine, {})
    tracker = trackers.get(key)
    if tracker is None:
        tracker = trackers[key] = SideGroupReleaseTracker(observer, engine)
    return tracker


class TimeSeries(list):
    """The observations of one ``write_time_series`` call, and why it stopped:
    ``stop_reason`` (one of ``STOP_REASONS``), ``events`` (fired in this call),
    ``total_events`` (the engine's count) and ``time_s`` (engine time)."""

    def __init__(self, observations=(), *, stop_reason="horizon", events=0, total_events=0, time_s=0.0):
        super().__init__(observations)
        self.stop_reason = stop_reason
        self.events = events
        self.total_events = total_events
        self.time_s = time_s


def write_time_series(engine, times: Iterable[float], path, *, xs: int, check: bool = True,
                      observer: Optional[MeltObserver] = None, append: bool = False,
                      side_group_homolysis=None, max_events: Optional[int] = None) -> TimeSeries:
    """Advance ``engine`` to each requested time (``run(until=t)``, the horizon
    API: observing never changes the trajectory) and write one observation per
    time as a JSON line; with ``check`` the exact invariants are asserted after
    every event (``invariant_checks`` counts them).

    ``max_events`` caps the events fired by this call (all horizons together).
    When the cap is reached before a requested time, the call stops there and
    that observation is not written; a time the engine has reached (the
    capping event fired exactly at it) is written.  The stop reason is
    "max_events" exactly when a requested time was not reached because of the
    cap; otherwise it is the last advance's (``RunResult.stop_reason``):
    "horizon", or "absorbing" (no enabled channel) at the last requested time,
    and "horizon" when no advance was needed or when the cap and the last
    requested time are reached by the same event.  It ends the file as a stop
    record and is returned with the observations (``TimeSeries``)."""
    if max_events is not None and (isinstance(max_events, bool) or not isinstance(max_events, int)
                                   or max_events < 0):
        raise ObservablesError(f"max_events must be None or a non-negative integer, got {max_events!r}")
    observer = observer if observer is not None else MeltObserver(
        engine.spec, xs=xs, identity=engine.identity, side_group_homolysis=side_group_homolysis)
    checker = InvariantChecker(engine.state) if check else None
    tracker = _run_tracker(engine, observer)
    hook = _Checks(checker, tracker)
    times = [float(t) for t in times]
    if any(not math.isfinite(t) for t in times) or times != sorted(times) or (times and times[0] < engine.time):
        raise ObservablesError(f"observation times must be finite, sorted and >= the engine time {engine.time}")
    out = TimeSeries()
    start = engine.events_fired
    stop = "horizon"
    with open(path, "a" if append else "w") as handle:
        for t in times:
            tracker.sync(engine)                 # catch up before advancing
            if t > engine.time:
                left = None if max_events is None else max_events - (engine.events_fired - start)
                stop = engine.run(until=t, max_events=left, check=hook).stop_reason
                if stop not in STOP_REASONS:
                    raise ObservablesError(f"engine stopped with unknown reason {stop!r}")
                if stop == "max_events":
                    if engine.time < t:
                        break                    # t not reached: no observation
                    # the capping event fired exactly at t: t is reached and
                    # observed; the series stops for the cap only if a later
                    # requested time remains (the next advance has no events left)
                    stop = "horizon"
            obs = observer.observe(engine, releases=tracker)
            obs["time_s"] = t
            obs["invariant_checks"] = checker.checks if checker is not None else None
            handle.write(json.dumps(obs, sort_keys=True, allow_nan=False) + "\n")
            handle.flush()
            out.append(obs)
        out.stop_reason, out.events = stop, engine.events_fired - start
        out.total_events, out.time_s = engine.events_fired, engine.time
        record = {"record": "stop", "schema": STOP_SCHEMA, "stop_reason": stop, "events": out.events,
                  "total_events": out.total_events, "observations": len(out), "time_s": out.time_s}
        handle.write(json.dumps(record, sort_keys=True, allow_nan=False) + "\n")
        handle.flush()
    return out


def mom_reported_quantities(obs: Mapping) -> dict:
    """The kMC state of one observation in the quantities the MoM solver
    reports: tail moments (pool moment slots), the physical end-group
    inventory of the tail, explicit-ladder amounts per DP, Mn / Mw of the
    tail moments (``observables_from_moments``), the amount of every declared
    species (0 when absent), the polymer-derived molecules, and the physical
    observables compared directly: the condensed polymer mass (the exact kMC
    polymer-derived mass; MoM's mass law on the kMC moments is
    ``mom_mass_law.value_g``, its gap the ``convention_gap_g`` row), the
    polymer conversion from it, and per declared side-group channel the
    released gas-product amount (None when untracked)."""
    if obs.get("length_rule") == "mean_residue_molar_mass":
        raise ObservablesError(MOM_V2_REFUSAL)
    t1, t3, law = obs["tier1"], obs["tier3"], obs["mom_mass_law"]
    mu = list(t3["tail"]["dp_moments_mol"])
    tail = observables_from_moments(mu, obs["repeat_unit_molar_mass_g_mol"], obs["end_group_molar_mass_g_mol"])
    return {
        "mu0_mol": mu[0], "mu1_mol": mu[1], "mu2_mol": mu[2],
        "end_inventory_mol": law["end_inventory_mol"],
        "explicit_dp_mol": dict(t3["explicit"]["dp_amount_mol"]),
        "polymer_condensed_mass_g": t1["polymer_derived"]["mass_g"],
        "Mn_g_mol": tail.mn_g_mol, "Mw_g_mol": tail.mw_g_mol,
        "polymer_conversion": 1.0 - t1["polymer_derived"]["mass_g"] / obs["initial_polymer_mass_g"],
        "released_mol": {row["label"]: row["released_mol"] for row in obs.get("side_group_channels", ())},
        "species_mol": {k: v["amount_mol"] for k, v in t1["species"].items()},
        "polymer_derived_molecules_mol": t1["polymer_derived"]["mom_definition"]["amount_mol"],
    }


_TIER1 = ("polymer_condensed_mass_g", "polymer_derived_molecules_mol", "polymer_conversion", "convention_gap_g")
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
    ``explicit_dp_mol`` / ``species_mol`` / ``side_group_Z_mol`` mappings).

    One row per MoM quantity -- every key of ``mom_reported_quantities``,
    every declared species, every key ``mom`` supplies (a DP, species or side
    group the kMC run never produced has kMC value 0).  MoM numbers are
    compared ONLY with the kMC value in the same convention (the kMC state
    expressed by the MoM definitions, e.g. the condensed mass by MoM's mass
    law), and a compared row is ``equal`` (relative difference <= 1e-12) or
    ``differs``; ``missing_in_mom`` / ``missing_in_kmc`` when a side has no
    value; without ``mom`` every row is ``no_mom_values``.  No convention
    bound ever turns a kMC-vs-MoM difference into agreement.  The kMC
    convention gap (exact polymer-derived mass minus its own mass-law value)
    is its own row, ``convention_gap_g`` with its bound, labelled
    ``not_a_comparison``."""
    kmc = _flatten(mom_reported_quantities(obs))
    for name in obs.get("declared_species", ()):
        kmc.setdefault(f"species_mol.{name}", 0.0)
    given = _flatten(mom) if mom is not None else {}
    keys = list(kmc) + [k for k in given if k not in kmc and k != "convention_gap_g"]
    rows = []
    for key in keys:
        k_value = kmc.get(key)
        if k_value is None and key.startswith(_ZERO_WHEN_ABSENT):
            k_value = 0.0
        m_value = given.get(key)
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
            agreement = "equal" if relative <= EQUAL_RELATIVE else "differs"
        tier = 1 if key in _TIER1 or key.startswith(("species_mol.", "released_mol.")) else 3
        rows.append({"tier": tier, "quantity": key, "kmc": k_value, "mom": m_value, "difference": difference,
                     "relative_difference": relative, "bound": None, "agreement": agreement})
    law = obs["mom_mass_law"]
    rows.append({"tier": 1, "quantity": "convention_gap_g", "kmc": float(law["difference_g"]), "mom": None,
                 "difference": None, "relative_difference": None, "bound": float(law["bound_g"]),
                 "agreement": "not_a_comparison"})
    return rows


__all__ = [
    "CONVENTION_DEFINITION", "ComponentRecord", "InvariantChecker", "InvariantError", "MASS_LAW_DEFINITION",
    "SIDE_GROUP_DEFINITION", "BALANCE_DEFINITION", "SideGroupReleaseTracker", "site_class", "aromatic_atoms",
    "aromatic_atoms_whole", "channel_identity", "clear_aromaticity_cache", "PERCEPTION_STATS",
    "MeltObserver", "OBSERVABLES_SCHEMA", "ObservablesError", "SPECIES_RULE", "STOP_REASONS", "STOP_SCHEMA",
    "TimeSeries", "TIER1_NOTE", "Totals",
    "backbone_atoms_per_unit", "canonical_key", "declared_species_names", "exact_totals",
    "side_group_site_indices", "validate_side_group_channels",
    "mom_reported_quantities", "molar_mass_g_mol", "observe", "side_by_side", "write_time_series",
]
