"""S4 melt engine: the loaded template catalog fired on an atom-explicit melt.

The engine couples four layers that already exist and adds nothing chemical:

* the state is ``melt_state.MeltState`` (atoms are the only authority, a
  molecule is a connected component, id = minimum atom id);
* the chemistry is the ``window_matcher.Template`` set a catalog loads at the
  requested temperature (``template_emission.loaded_templates``); the engine
  never edits a weight, and temperature enters only through that call;
* transport is a declared law (``TransportLaw``) read from the polymer
  declaration's ``transport_model["met"]``; nothing here names a polymer;
* time is the observation-independent clock: a waiting time that crosses a
  horizon is kept and reused while the state, the propensities and the
  environment are unchanged, and a rejected (null) proposal consumes its draw.

Index (two layers)
------------------
Structural matches are recomputed over every touched component after each
event (no dirty radius): per unimolecular template the deduplicated
transitions of the component, split into allowed and leak; per bimolecular
template and participant the edit-role-distinct sites of the component,
grouped by their *branching kind* (below).  Per component the index also keeps
the MET length ``n_c = max(1, floor(B_c / b))`` (``B_c`` backbone-origin atoms,
``b`` backbone atoms per repeat unit), the radical count and the declared
species the component is (if any).

Global aggregates ``S_p = sum_c n_p(c)`` and ``X_pq = sum_c n_p(c) n_q(c)`` are
kept per population ``p = (template, participant, kind)``; a component's
contribution is removed and re-added whenever it is touched, so a split
enables (and a merge disables) remote pairs through the ``sum_c`` term.

Propensities
------------
* unimolecular: ``a = w * |transitions|``;
* bimolecular: ``a = w * P / (N_A V)`` with ``P = n0 n1 - sum_c n0(c) n1(c)``
  (ordered pairs on different components), divided by 2 under ``role_swap``
  (and by the order of any participant-preserving edit symmetry);
* MET (radical-decreasing bimolecular classes): per component pair the
  activation ``k_act(c, c') = sum over classes of w T`` (T: the class's
  transitions between c and c', branch-creating ones included) is summed
  FIRST and Collins-Kimball is applied once:
  ``k = k_act * phi``, ``phi = k_D / (k_D + k_act)``, ``k_D`` from the MET
  lengths of c and c'; the allowed part ``phi * k_allowed / (N_A V)`` fires,
  the branch-creating part accrues as leak.  Within the pair a class is drawn
  proportional to ``w T`` and a transition uniformly.  The reverse of
  every MET class (any radical-increasing class: homolysis, reverse
  disproportionation) fires by exact thinning: it is proposed at its
  intrinsic weight, applied tentatively, ``phi`` of the PRODUCT pair is
  computed from the product components' own index entries, and the event is
  accepted with probability ``phi`` (else rolled back: a null proposal), so
  ``K = k_act / k_diss = Kc`` is preserved.  A pair one of whose components is
  a declared small species (``spec.small_species``, ``monomer_product``, the
  built-ins H and H2) bypasses MET in both directions (``phi = 1``) and is
  counted.

Exact pair sampler
------------------
``ComponentPairSampler``: the left component is drawn with weight
``n0(c) (n1 - n1(c))`` by a Fenwick descent over component slots (slot = the
component id, i.e. an atom id), the left site uniformly in it, the right site
uniformly among the ``n1 - n1(c)`` sites outside it.

Branching leak
--------------
A transition whose formed bond would give a backbone-origin atom a third
backbone-origin neighbour is refused at the transition level.  Unimolecular
transitions are checked one by one.  A bimolecular site carries a kind (per
cross-participant formed bond: whether its atom is backbone and its backbone
degree after the participant's own edits); a kind pair is allowed or leak, so
the leak and allowed pair counts are exact.  Leak propensity accrues to the
leak integral (``leak_num``) and never fires.

Matching cost: structural matches use the matcher's own embedding, finish
and evaluate steps on a component-local view, with one optimization: hydrogen
leaves that no edit names are enumerated in increasing atom order only
(``participant_embeddings``); this leaves every deduplicated transition and
site unchanged (tested against the full enumeration).

Transport declaration (``transport_model["met"]``, schema ``MET_SCHEMA``)::

    schema: "kmc-met-transport/1"
    segment_diffusivity: {law: wlf, D_ref_m2_s, T_ref_K, C1, C2_K, T_wlf_K}
                       | {law: arrhenius, D_ref_m2_s, T_ref_K, Ea_J_mol}
    chain_scaling: {law: rouse} | {law: crossover, entanglement_units}
    mean_square_end_to_end_per_molar_mass_m2_mol_per_g: <R^2>/M
    contact_radius_m: sigma
    spin_factor: s in [0.25, 1]

``M0`` and ``b`` come from the declared repeat unit.  A declaration without
the block, or with a missing or unknown field, is refused at run start.

Not here: observables tiers, the MoM comparison, captures/J_ring (R0 only),
evaporation, volume change.
"""

from __future__ import annotations

import math
import weakref
from collections import Counter
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple

from rmgpy import constants
from rmgpy.kmc import window_matcher as wm
from rmgpy.kmc.melt_state import MeltState, SetBond, SetRadical, materialize, _Builder, _backbone_keys

MET_SCHEMA = "kmc-met-transport/1"
SPIN_FACTOR_RANGE = (0.25, 1.0)
BUILTIN_SPECIES = ("H", "H2")


class EngineError(ValueError):
    """The engine cannot proceed (an inconsistent index, a broken invariant)."""


class DeclarationError(EngineError):
    """The polymer declaration is refused at run start."""


# --------------------------------------------------------------------------
# transport (declared, generic)
# --------------------------------------------------------------------------
def _number(value, label, *, positive=True):
    try:
        out = float(value)
    except (TypeError, ValueError):
        raise DeclarationError(f"transport_model.met: {label} must be a number, got {value!r}") from None
    if not math.isfinite(out) or (positive and out <= 0.0):
        raise DeclarationError(f"transport_model.met: {label} must be finite and positive, got {value!r}")
    return out


@dataclass(frozen=True)
class TransportLaw:
    """The declared MET transport law of one linear polymer.

    ``k_D(n1, n2) = 4 pi N_A s (D(n1) + D(n2)) R``,
    ``R = max(sigma, 2 Rg(min(n1, n2)))``, ``Rg(n)^2 = C M0 n / 6``,
    ``D(n) = D0(T) / n`` (Rouse) or ``D0 N_e / n^2`` above ``N_e`` (crossover).
    ``M0`` (g/mol) is the declared repeat unit's molar mass and ``b`` its
    backbone atom count (head-port to tail-port path); every other number is
    read from ``transport_model["met"]`` (schema ``MET_SCHEMA``)."""

    segment_diffusivity: Tuple
    chain_scaling: Tuple
    r2_per_molar_mass: float
    contact_radius_m: float
    spin_factor: float
    repeat_unit_molar_mass_g_mol: float
    backbone_atoms_per_unit: int

    @classmethod
    def from_spec(cls, spec) -> "TransportLaw":
        model = getattr(spec, "transport_model", None)
        if not isinstance(model, Mapping) or "met" not in model:
            raise DeclarationError(
                "the polymer declaration has no MET transport parameters: "
                "transport_model['met'] is required (schema " + MET_SCHEMA + ": "
                "segment_diffusivity, chain_scaling, "
                "mean_square_end_to_end_per_molar_mass_m2_mol_per_g, contact_radius_m, "
                "spin_factor); the engine has no built-in polymer constants")
        met = model["met"]
        if not isinstance(met, Mapping):
            raise DeclarationError("transport_model['met'] must be a mapping")
        if met.get("schema") != MET_SCHEMA:
            raise DeclarationError(f"transport_model['met'].schema must be {MET_SCHEMA!r}, "
                                   f"got {met.get('schema')!r}")
        required = {"schema", "segment_diffusivity", "chain_scaling",
                    "mean_square_end_to_end_per_molar_mass_m2_mol_per_g",
                    "contact_radius_m", "spin_factor"}
        missing = required - set(met)
        unknown = set(met) - required
        if missing:
            raise DeclarationError(f"transport_model['met'] lacks {sorted(missing)}")
        if unknown:
            raise DeclarationError(f"transport_model['met'] has unknown fields {sorted(unknown)}")
        d0 = met["segment_diffusivity"]
        if not isinstance(d0, Mapping):
            raise DeclarationError("segment_diffusivity must be a mapping")
        law = d0.get("law")
        if law == "wlf":
            keys = ("D_ref_m2_s", "T_ref_K", "C1", "C2_K", "T_wlf_K")
        elif law == "arrhenius":
            keys = ("D_ref_m2_s", "T_ref_K", "Ea_J_mol")
        else:
            raise DeclarationError(f"segment_diffusivity.law must be wlf or arrhenius, got {law!r}")
        if set(d0) != {"law", *keys}:
            raise DeclarationError(f"segment_diffusivity ({law}) needs exactly {sorted(keys)}")
        d0_values = tuple([law] + [_number(d0[k], f"segment_diffusivity.{k}",
                                           positive=(k != "Ea_J_mol")) for k in keys])
        scaling = met["chain_scaling"]
        if not isinstance(scaling, Mapping):
            raise DeclarationError("chain_scaling must be a mapping")
        if scaling.get("law") == "rouse" and set(scaling) == {"law"}:
            scaling_values = ("rouse",)
        elif scaling.get("law") == "crossover" and set(scaling) == {"law", "entanglement_units"}:
            scaling_values = ("crossover", _number(scaling["entanglement_units"],
                                                   "chain_scaling.entanglement_units"))
        else:
            raise DeclarationError("chain_scaling must be {law: rouse} or "
                                   "{law: crossover, entanglement_units}")
        spin = _number(met["spin_factor"], "spin_factor")
        if not SPIN_FACTOR_RANGE[0] <= spin <= SPIN_FACTOR_RANGE[1]:
            raise DeclarationError(f"spin_factor must lie in {SPIN_FACTOR_RANGE}, got {spin}")
        repeat = spec.repeat_unit
        m0 = float(repeat.composition.molar_mass_g_mol())
        b = len(_backbone_keys(repeat))
        return cls(d0_values, scaling_values,
                   _number(met["mean_square_end_to_end_per_molar_mass_m2_mol_per_g"],
                           "mean_square_end_to_end_per_molar_mass_m2_mol_per_g"),
                   _number(met["contact_radius_m"], "contact_radius_m"), spin, m0, b)

    def d0(self, temperature: float) -> float:
        t = float(temperature)
        law = self.segment_diffusivity[0]
        if law == "wlf":
            _l, d_ref, t_ref, c1, c2, t_wlf = self.segment_diffusivity

            def shift(x):
                return -c1 * (x - t_wlf) / (c2 + x - t_wlf)
            return d_ref * 10.0 ** (shift(t_ref) - shift(t))
        _l, d_ref, t_ref, ea = self.segment_diffusivity
        return d_ref * math.exp(-ea / constants.R * (1.0 / t - 1.0 / t_ref))

    def diffusivity(self, temperature: float, units: int) -> float:
        d0 = self.d0(temperature)
        if self.chain_scaling[0] == "rouse" or units <= self.chain_scaling[1]:
            return d0 / units
        return d0 * self.chain_scaling[1] / units ** 2

    def k_diffusion(self, temperature: float, n1: int, n2: int) -> float:
        """Encounter rate constant (m^3 mol^-1 s^-1) of two MET lengths."""
        rg = math.sqrt(self.r2_per_molar_mass * self.repeat_unit_molar_mass_g_mol * min(n1, n2) / 6.0)
        radius = max(self.contact_radius_m, 2.0 * rg)
        return (4.0 * math.pi * constants.Na * self.spin_factor
                * (self.diffusivity(temperature, n1) + self.diffusivity(temperature, n2)) * radius)

    def met_length(self, backbone_atoms: int) -> int:
        return max(1, backbone_atoms // self.backbone_atoms_per_unit)


def escape_factor(k_act: float, k_d: float) -> float:
    """phi = k_D / (k_D + k_act): the Collins-Kimball factor of a summed activation."""
    return k_d / (k_d + k_act)


# --------------------------------------------------------------------------
# catalog
# --------------------------------------------------------------------------
class PayloadCatalog:
    """Emitted template payloads; the templates at T come only from
    ``template_emission.loaded_templates(payloads, T)``."""

    def __init__(self, payloads):
        self.payloads = list(payloads)

    def templates(self, temperature: float):
        from rmgpy.kmc.template_emission import loaded_templates
        return loaded_templates(self.payloads, float(temperature))


@dataclass(frozen=True)
class _Class:
    """Static classification of one loaded template."""

    index: int
    name: str
    family: str
    bimolecular: bool
    radical_delta: int
    met_forward: bool
    met_reverse: bool
    divisor: int
    cross_bonds: Tuple[Tuple[int, int], ...]          # formed bonds across participants
    internal: Tuple[Tuple[Tuple[int, int], float, float], ...]  # (pair, before, after) inside a participant


def _classify(index: int, t: wm.Template) -> _Class:
    before = sum(a.radical for a in t.atoms)
    after = sum(t.atom_final.get(i, {}).get("radical", a.radical) for i, a in enumerate(t.atoms))
    delta = after - before
    bi = len(t.participants) == 2
    cross, internal = [], []
    for (i, j), order in sorted(t.bond_final.items()):
        old = t.bond.get((i, j), 0.0)
        if t.atoms[i].participant != t.atoms[j].participant:
            if old == 0.0 and order > 0.0:
                cross.append((i, j) if t.atoms[i].participant == 0 else (j, i))
        else:
            internal.append(((i, j), old, order))
    divisor = len(t.edit_symmetry) * (2 if t.role_swap else 1) if bi else 1
    return _Class(index, t.name, t.family, bi, delta, bi and delta < 0, delta > 0, divisor,
                  tuple(cross), tuple(internal))


# --------------------------------------------------------------------------
# component-local matching view
# --------------------------------------------------------------------------
class _Node:
    __slots__ = ("element", "isotope", "charge", "radical", "lone_pairs", "edges")

    def __init__(self, node):
        self.element = node.element
        self.isotope = node.isotope
        self.charge = node.charge
        self.radical = node.radical
        self.lone_pairs = node.lone_pairs
        self.edges = dict(node.edges)


class _LocalView(wm._View):
    """``window_matcher._View`` restricted to whole components (components are
    closed under adjacency, so every matcher read stays inside); nodes carry
    plain-dict edges.  Matching on it equals matching on the full state
    restricted to these components."""

    def __init__(self, state, cids):
        self._state = weakref.ref(state)
        self.revision = state.revision
        comps = state.components
        nodes = state.atoms
        members = {c: comps.atoms(c) for c in sorted(set(cids))}
        self.ids = tuple(sorted(a for atoms in members.values() for a in atoms))
        self.node = {a: _Node(nodes[a]) for a in self.ids}
        self.adj = {a: tuple(sorted(self.node[a].edges)) for a in self.ids}
        self.sig, self.free = {}, {}
        by_sig = {}
        for a in self.ids:
            n = self.node[a]
            sig = (n.element, n.isotope, n.charge, n.lone_pairs, len(n.edges))
            self.sig[a] = sig
            self.free[a] = wm._valence(n.element) - n.charge - 2 * n.lone_pairs - len(n.edges)
            by_sig.setdefault(sig, []).append(a)
        self.by_sig = {k: tuple(v) for k, v in by_sig.items()}
        self.in_ring = wm._ring_atoms(self.adj, self.ids)
        self.comp = {a: c for c, atoms in members.items() for a in atoms}
        self.subs, self.sub_of, self.comp_forms = [], {}, {}
        for cid, atoms in members.items():
            subs, total = wm._component_subs(atoms, self._value)
            self.comp_forms[cid] = total
            for sub in subs:
                for a in sub.atoms:
                    self.sub_of[a] = len(self.subs)
                self.subs.append(sub)
        self._base_rows, self._dist, self._hash = {}, {}, None


# --------------------------------------------------------------------------
# H-leaf quotient enumeration (a matcher optimization only)
# --------------------------------------------------------------------------
_H_ORDER: Dict[int, Tuple] = {}


def _h_order(template: wm.Template) -> Tuple:
    """Per participant piece, per plan position: the plan positions of the
    earlier interchangeable hydrogen leaves (same parent, not edited) whose
    image must be smaller.  Embeddings that differ only by permuting such
    leaves give the same edits, site, resonance choices and outcome, so
    keeping the increasing assignment alone leaves every deduplicated
    transition and every site unchanged (the embedding COUNT is not kept;
    the engine never reads it).  Equivalence to full enumeration is tested."""
    key = id(template)
    cached = _H_ORDER.get(key)
    if cached is not None and cached[0] is template:
        return cached[1]
    edited = set(template.edited_all)
    out = []
    for pieces in template.pieces:
        per_piece = []
        for piece in pieces:
            seen: Dict[int, List[int]] = {}
            constraint = []
            for t, parent, _checks in piece.plan:
                atom = template.atoms[t]
                if (atom.element == "H" and parent is not None and t not in edited
                        and atom.kind == "closed" and len(template.neighbours[t]) == 1):
                    constraint.append(tuple(seen.get(parent, ())))
                    seen.setdefault(parent, []).append(t)
                else:
                    constraint.append(())
            per_piece.append(tuple(constraint))
        out.append(tuple(per_piece))
    out = tuple(out)
    _H_ORDER[key] = (template, out)
    return out


def _piece_embeddings_quotient(view, template, piece, order):
    plan = piece.plan
    root = plan[0][0]
    results = []
    mapping = {}
    used = set()
    nodes = view.node

    def visit(pos):
        if pos == len(plan):
            results.append(tuple((t, mapping[t]) for t, _p, _c in plan))
            return
        t, parent, checks = plan[pos]
        before = order[pos]
        sig = template.sig[t]
        for a in view.adj[mapping[parent]]:
            if a in used or view.sig[a] != sig:
                continue
            if before and any(mapping[u] > a for u in before):
                continue
            edges = nodes[a].edges
            if any((mapping[u] in edges) != bonded for u, bonded in checks):
                continue
            mapping[t] = a
            used.add(a)
            visit(pos + 1)
            used.discard(a)
            del mapping[t]

    for a in view.by_sig.get(template.sig[root], ()):
        mapping[root] = a
        used.add(a)
        visit(1)
        used.discard(a)
        del mapping[root]
    return results


def participant_embeddings(view, template, participant, *, quotient=True):
    """``window_matcher._participant_embeddings`` with the hydrogen-leaf
    quotient (``quotient=False``: the matcher's own full enumeration)."""
    if not quotient:
        return wm._participant_embeddings(view, template, participant)
    orders = _h_order(template)[participant]
    per_piece = [_piece_embeddings_quotient(view, template, piece, order)
                 for piece, order in zip(template.pieces[participant], orders)]
    if any(not p for p in per_piece):
        return []
    by_comp = []
    for embeddings in per_piece:
        grouped = {}
        for emb in embeddings:
            grouped.setdefault(view.comp[emb[0][1]], []).append(emb)
        by_comp.append(grouped)
    out = []
    for cid in sorted(by_comp[0]):
        if any(cid not in g for g in by_comp):
            continue
        for combo in wm._piece_products([g[cid] for g in by_comp], view):
            mapping = {}
            for emb in combo:
                mapping.update(emb)
            embedding = wm._finish(view, template, cid, mapping)
            if embedding is not None:
                out.append(embedding)
    return out


# --------------------------------------------------------------------------
# transitions and index entries
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class EngineTransition:
    """One transition ready to apply (atom-level; independent of any view)."""

    template: str
    family: str
    class_index: int
    edit_key: Tuple[int, ...]
    components: Tuple[int, ...]
    edits: Tuple
    form_bonds: Tuple
    form_radicals: Tuple
    delta: Tuple
    leak: bool = False


@dataclass
class _Entry:
    """Index entry of one component."""

    cid: int
    atoms: Tuple[int, ...]
    uni: Dict[int, Tuple[Tuple[EngineTransition, ...], int]]   # class -> (allowed, leak count)
    sites: Dict[Tuple[int, int, Tuple], Tuple[Tuple[int, ...], ...]]  # (class, participant, kind) -> sites
    backbone: int
    met_length: int
    radicals: int
    species: Optional[str]

    def counts(self) -> Dict[Tuple[int, int, Tuple], int]:
        return {k: len(v) for k, v in self.sites.items()}

    def canonical(self):
        return (self.atoms,
                tuple(sorted((k, tuple(t.edit_key for t in v[0]), v[1]) for k, v in self.uni.items())),
                tuple(sorted(self.sites.items())), self.backbone, self.met_length,
                self.radicals, self.species)


def _bb_degree(state, a, removed=(), added=0):
    origin = state.origin
    edges = state.atoms[a].edges
    return sum(1 for b in edges if origin[b].is_backbone and b not in removed) + added


def _creates_branch(state, edits) -> bool:
    """A formed bond leaves some backbone-origin atom with >= 3 backbone-origin
    neighbours (edits are concrete ``SetBond``s on state atoms)."""
    origin = state.origin
    formed, broken = {}, {}
    for e in edits:
        if isinstance(e, SetBond) and (e.old == 0) != (e.new == 0):
            target = formed if e.new != 0 else broken
            target.setdefault(e.a, set()).add(e.b)
            target.setdefault(e.b, set()).add(e.a)
    for a, partners in formed.items():
        if not origin[a].is_backbone:
            continue
        gained = sum(1 for b in partners if origin[b].is_backbone)
        if gained and _bb_degree(state, a, broken.get(a, ()), gained) >= 3:
            return True
    return False


def _site_kind(state, cls: _Class, t: wm.Template, participant: int, site) -> Tuple:
    """Branching kind of a bimolecular site: (internal leak, ((template atom,
    backbone?, backbone degree after this participant's own edits), ...))."""
    pos = {ti: a for ti, a in zip(t.edited[participant], site)}
    removed, added = {}, {}
    internal_formed = []
    for (i, j), old, new in cls.internal:
        if t.atoms[i].participant != participant:
            continue
        a, b = pos[i], pos[j]
        if old != 0 and new == 0:
            removed.setdefault(a, set()).add(b)
            removed.setdefault(b, set()).add(a)
        elif old == 0 and new != 0:
            internal_formed.append((a, b))
            for x, y in ((a, b), (b, a)):
                if state.origin[y].is_backbone:
                    added[x] = added.get(x, 0) + 1
    origin = state.origin

    def degree(a):
        return _bb_degree(state, a, removed.get(a, ()), added.get(a, 0))

    internal_leak = any(origin[x].is_backbone and degree(x) >= 3
                        for a, b in internal_formed for x in (a, b))
    cross = []
    for i0, i1 in cls.cross_bonds:
        ti = i0 if participant == 0 else i1
        a = pos[ti]
        cross.append((ti, bool(origin[a].is_backbone), degree(a)))
    return (internal_leak, tuple(sorted(set(cross))))


def _pair_leak(cls: _Class, kind0: Tuple, kind1: Tuple) -> bool:
    if kind0[0] or kind1[0]:
        return True
    info = {}
    for ti, bb, deg in kind0[1] + kind1[1]:
        info[ti] = (bb, deg)
    gained = Counter()
    for i0, i1 in cls.cross_bonds:
        if info[i1][0]:
            gained[i0] += 1
        if info[i0][0]:
            gained[i1] += 1
    return any(info[ti][0] and info[ti][1] + g >= 3 for ti, g in gained.items())


# --------------------------------------------------------------------------
# Fenwick trees and the exact pair sampler
# --------------------------------------------------------------------------
class _Fenwick:
    __slots__ = ("n", "tree", "top")

    def __init__(self, n: int):
        self.n = n
        self.tree = [0] * (n + 1)
        self.top = 1
        while self.top * 2 <= n:
            self.top *= 2

    def add(self, k: int, delta):
        i = k + 1
        while i <= self.n:
            self.tree[i] += delta
            i += i & -i

    def prefix(self, k: int):
        """Sum over slots [0, k)."""
        total, i = 0, k
        while i > 0:
            total += self.tree[i]
            i -= i & -i
        return total

    def find(self, r):
        """Least slot whose inclusive prefix exceeds r (0 <= r < total)."""
        pos, step = 0, self.top
        while step:
            nxt = pos + step
            if nxt <= self.n and self.tree[nxt] <= r:
                r -= self.tree[nxt]
                pos = nxt
            step //= 2
        return pos, r


class ComponentPairSampler:
    """Exact uniform sampler of ordered (participant-0 site, participant-1
    site) pairs on different components, over a dynamic component set.

    Slots are component ids (atom ids, ``0 <= cid < size``).  The left
    component is drawn with weight ``n0(c) (n1 - n1(c))`` by a Fenwick descent
    over two trees (``n0(c)`` and ``n0(c) n1(c)``), then a left site uniformly
    in it and a right site uniformly among the ``n1 - n1(c)`` sites outside
    (a third tree, ``n1(c)``)."""

    def __init__(self, size: int):
        self.size = int(size)
        self._a = _Fenwick(self.size)
        self._b = _Fenwick(self.size)
        self._c = _Fenwick(self.size)
        self._s0: Dict[int, Tuple] = {}
        self._s1: Dict[int, Tuple] = {}
        self.n0 = 0
        self.n1 = 0
        self.x = 0

    def set(self, cid: int, sites0=(), sites1=()):
        old0, old1 = len(self._s0.get(cid, ())), len(self._s1.get(cid, ()))
        sites0, sites1 = tuple(sites0), tuple(sites1)
        new0, new1 = len(sites0), len(sites1)
        if sites0:
            self._s0[cid] = sites0
        else:
            self._s0.pop(cid, None)
        if sites1:
            self._s1[cid] = sites1
        else:
            self._s1.pop(cid, None)
        if new0 != old0:
            self._a.add(cid, new0 - old0)
        if new0 * new1 != old0 * old1:
            self._b.add(cid, new0 * new1 - old0 * old1)
        if new1 != old1:
            self._c.add(cid, new1 - old1)
        self.n0 += new0 - old0
        self.n1 += new1 - old1
        self.x += new0 * new1 - old0 * old1

    @property
    def ordered_total(self) -> int:
        return self.n0 * self.n1 - self.x

    def _left_component(self, rng) -> int:
        r = int(rng.integers(self.ordered_total))
        pos, step = 0, self._a.top
        a, b, n1 = self._a.tree, self._b.tree, self.n1
        while step:
            nxt = pos + step
            if nxt <= self.size:
                w = n1 * a[nxt] - b[nxt]
                if w <= r:
                    r -= w
                    pos = nxt
            step //= 2
        return pos

    def draw(self, rng):
        """(left component, left site, right component, right site)."""
        if self.ordered_total <= 0:
            raise EngineError("no pair to draw")
        c0 = self._left_component(rng)
        s0 = self._s0[c0]
        left = s0[int(rng.integers(len(s0)))]
        n1c = len(self._s1.get(c0, ()))
        g = int(rng.integers(self.n1 - n1c))
        if g >= self._c.prefix(c0):
            g += n1c
        c1, g = self._c.find(g)
        return c0, left, c1, self._s1[c1][g]


# --------------------------------------------------------------------------
# declared species identity
# --------------------------------------------------------------------------
def _graph(rows):
    import networkx as nx
    g = nx.Graph()
    for a, row in rows.items():
        g.add_node(a, label=row[:5])
        for b, o in row[5]:
            if a < b:
                g.add_edge(a, b, order=float(o))
    return g


class SpeciesIdentity:
    """Exact isomorphism (under resonance) against the declared species."""

    def __init__(self, spec):
        declared = [(name, graph) for name, graph, _d in spec.small_species]
        if getattr(spec, "monomer_product", None) is not None:
            declared.append(("monomer_product", spec.monomer_product))
        self.forms = []
        for name, graph in declared:
            builder = _Builder()
            builder.add_segment(graph, "mobile", 0, name)
            self._add(name, builder.state())
        from rmgpy.kmc.melt_state import Node, Origin
        h1 = {0: Node("H", None, 0, 1, 0, {})}
        h2 = {0: Node("H", None, 0, 0, 0, {1: 1.0}), 1: Node("H", None, 0, 0, 0, {0: 1.0})}
        self._add("H", MeltState(h1, {0: Origin("mobile", 0, "H", "H", False)}))
        self._add("H2", MeltState(h2, {i: Origin("mobile", 0, "H2", f"H{i}", False) for i in (0, 1)}))
        self.max_atoms = max(len(rows) for _n, _c, rows in self.forms)
        self._cache: Dict[Tuple, Optional[str]] = {}

    def _add(self, name, state):
        cid = state.components.ids()[0]
        comp = tuple(sorted(state.components.composition(cid).items(), key=str))
        for form in wm.component_forms(state, cid):
            self.forms.append((name, comp, wm.form_rows(state, cid, form)))

    def identify(self, state, cid) -> Optional[str]:
        atoms = state.components.atoms(cid)
        if len(atoms) > self.max_atoms:
            return None
        comp = tuple(sorted(state.components.composition(cid).items(), key=str))
        nodes = state.atoms
        rows = {a: (nodes[a].element, nodes[a].isotope, nodes[a].charge, nodes[a].radical,
                    nodes[a].lone_pairs, tuple(sorted(nodes[a].edges.items()))) for a in atoms}
        key = (comp, tuple(sorted((r[:5], tuple(sorted(o for _b, o in r[5]))) for r in rows.values())))
        if key in self._cache:
            return self._cache[key]
        import networkx as nx
        found = None
        g = _graph(rows)
        for name, c, declared in self.forms:
            if c == comp and nx.is_isomorphic(
                    g, _graph(declared), node_match=lambda x, y: x["label"] == y["label"],
                    edge_match=lambda x, y: x["order"] == y["order"]):
                found = name
                break
        self._cache[key] = found
        return found


# --------------------------------------------------------------------------
# events and results
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Event:
    time: float
    template: str
    family: str
    channel: str          # "uni", "bi" or "met"
    components: Tuple[int, ...]
    products: Tuple[int, ...]
    revision: int


@dataclass(frozen=True)
class Propensities:
    enabled: float
    leak: float
    channels: Tuple[Tuple[Tuple, float], ...]

    @property
    def total(self) -> float:
        return self.enabled


@dataclass(frozen=True)
class PairRate:
    k_act: float          # summed activation over every MET class and site pair (m3/mol/s)
    k_act_allowed: float
    k_d: float
    phi: float
    bypass: bool
    allowed: float        # propensity (1/s) of the allowed MET transitions of the pair
    leak: float           # propensity (1/s) of the refused (branching) ones
    classes: Tuple[Tuple[int, float, float], ...]   # (class, T allowed, T leak)


@dataclass(frozen=True)
class RunResult:
    time: float
    events: int
    total_events: int
    stop_reason: str
    leak_num: float
    leak_den: float
    leak: float
    null_proposals: int
    family_counts: Mapping[str, int]
    template_counts: Mapping[str, int]
    met: Mapping[str, int]


@dataclass(frozen=True)
class _Pending:
    elapsed: float
    absorbing: bool
    enabled: float
    leak: float
    horizon_time: float
    token: Tuple


# --------------------------------------------------------------------------
# the engine
# --------------------------------------------------------------------------
class MeltEngine:
    """Isothermal kMC of a declared linear-polymer melt over a loaded catalog.

    ``spec``: the ``LinearHomopolymerSpec`` (transport is read from it and the
    run is refused without it); ``catalog``: anything with
    ``templates(temperature)`` (``PayloadCatalog``); ``rng``: a
    ``numpy.random.Generator``; ``state``: a ``MeltState`` (default: the
    declaration materialized, declared contents only); ``volume_m3``: the
    constant volume (default: the declaration's representative volume)."""

    def __init__(self, spec, catalog, temperature: float, rng, *, state=None, volume_m3=None,
                 quotient: bool = True):
        from rmgpy.polymer_spec import LinearHomopolymerSpec
        if not isinstance(spec, LinearHomopolymerSpec):
            raise DeclarationError(
                f"only LinearHomopolymerSpec declarations are supported, got {type(spec).__name__}")
        self.spec = spec
        self.quotient = bool(quotient)
        self.transport = TransportLaw.from_spec(spec)
        self.catalog = catalog
        self.rng = rng
        self.state = materialize(spec) if state is None else state
        isotopes = sorted({n.isotope for n in self.state.atoms.values() if n.isotope is not None})
        if isotopes:
            raise DeclarationError(f"non-natural isotopes are refused in v1 (found mass numbers {isotopes})")
        self.volume = float(spec.initialization().volume_m3 if volume_m3 is None else volume_m3)
        if not math.isfinite(self.volume) or self.volume <= 0.0:
            raise DeclarationError("the volume must be finite and positive")
        self.nav = constants.Na * self.volume
        self.identity = SpeciesIdentity(spec)
        self.temperature = None
        self.templates: Tuple[wm.Template, ...] = ()
        self.classes: Tuple[_Class, ...] = ()
        self._load(temperature)
        # clock
        self.time = 0.0
        self.leak_num = self.leak_den = 0.0
        self._segment_start = 0.0
        self._leak_num_c = self._leak_den_c = 0.0
        self._pending: Optional[_Pending] = None
        self._stop = "horizon"
        # counters
        self.events: List[Event] = []
        self.events_fired = 0
        self.null_proposals = 0
        self.family_counts: Counter = Counter()
        self.template_counts: Counter = Counter()
        self.met_counts: Counter = Counter()
        self.conservation_checks = 0
        self._index_revision = 0
        self.rebuild()

    # -- catalog and temperature ----------------------------------------
    def _load(self, temperature):
        templates = tuple(self.catalog.templates(float(temperature)))
        names = [t.name for t in templates]
        if len(set(names)) != len(names):
            raise EngineError("template names must be unique")
        if self.templates:
            if names != [t.name for t in self.templates]:
                raise EngineError("the catalog loaded a different template set at the new temperature")
            for old, new in zip(self.templates, templates):
                if (old.atoms, old.bond, old.edits) != (new.atoms, new.bond, new.edits):
                    raise EngineError(f"template {new.name} changed structure with temperature")
        self.templates = templates
        self.classes = tuple(_classify(i, t) for i, t in enumerate(templates))
        self.weights = tuple(t.weight for t in templates)
        self.temperature = float(temperature)

    def set_temperature(self, temperature: float):
        """Reload the catalog's weights at ``temperature`` (structure unchanged)
        and refresh every rate; the pending draw is invalidated."""
        self._load(temperature)
        self._refresh_met()

    # -- index -------------------------------------------------------------
    def _scan(self, cids, state=None) -> Dict[int, _Entry]:
        state = self.state if state is None else state
        quotient = self.quotient
        out = {}
        for cid in cids:
            view = _LocalView(state, [cid])
            atoms = state.components.atoms(cid)
            uni, sites = {}, {}
            for cls, t in zip(self.classes, self.templates):
                if not cls.bimolecular:
                    found, cache = {}, {}
                    for emb in participant_embeddings(view, t, 0, quotient=quotient):
                        edits = wm._concrete_edits(t, emb.mapping)
                        for choice in emb.choices:
                            key = (edits, choice)
                            if key in cache:
                                continue
                            edit_key, identity, edits_, delta, fb, fr = wm._evaluate(
                                view, t, cls.index, emb.mapping, (cid,), choice)
                            cache[key] = identity
                            prev = found.get(identity)
                            if prev is None or edit_key < prev.edit_key:
                                found[identity] = EngineTransition(
                                    t.name, t.family, cls.index, edit_key, (cid,), edits_, fb, fr,
                                    delta, _creates_branch(state, edits_))
                    if found:
                        trs = sorted(found.values(), key=lambda x: x.edit_key)
                        allowed = tuple(x for x in trs if not x.leak)
                        uni[cls.index] = (allowed, len(trs) - len(allowed))
                else:
                    for p in (0, 1):
                        by = wm._sites(t, p, participant_embeddings(view, t, p, quotient=quotient)).get(cid, ())
                        grouped = {}
                        for site in by:
                            grouped.setdefault(_site_kind(state, cls, t, p, site), []).append(site)
                        for kind, ss in grouped.items():
                            sites[(cls.index, p, kind)] = tuple(sorted(ss))
            backbone = state.components.backbone_count(cid)
            out[cid] = _Entry(cid, atoms, uni, sites, backbone, self.transport.met_length(backbone),
                              state.components.radical_count(cid), self.identity.identify(state, cid))
        return out

    def rebuild(self):
        """Full recompute of every index layer from the current state."""
        size = max(self.state.atoms) + 1
        self._size = size
        self.entries: Dict[int, _Entry] = {}
        self.S: Counter = Counter()
        self.X: Counter = Counter()
        self.kinds: Dict[Tuple[int, int], set] = {}
        self.channels: Dict[Tuple[int, Tuple, Tuple], ComponentPairSampler] = {}
        self.uni_trees = {c.index: _Fenwick(size) for c in self.classes if not c.bimolecular}
        self.uni_allowed = Counter()
        self.uni_leak = Counter()
        self.met_comps: Dict[int, Dict] = {}
        self.pairs: Dict[Tuple[int, int], PairRate] = {}
        for entry in self._scan(self.state.components.ids()).values():
            self._add(entry)
        self._refresh_met()
        self._revision = self.state.revision
        self._index_revision += 1
        self._pending = None

    def _ensure_kind(self, ci, p, kind):
        known = self.kinds.setdefault((ci, p), set())
        if kind in known:
            return
        known.add(kind)
        if self.classes[ci].met_forward:
            return
        for other in sorted(self.kinds.get((ci, 1 - p), ()), key=repr):
            k0, k1 = (kind, other) if p == 0 else (other, kind)
            sampler = ComponentPairSampler(self._size)
            self.channels[(ci, k0, k1)] = sampler
            for cid, entry in self.entries.items():
                s0 = entry.sites.get((ci, 0, k0), ())
                s1 = entry.sites.get((ci, 1, k1), ())
                if s0 or s1:
                    sampler.set(cid, s0, s1)

    def _x_terms(self, entry: _Entry, sign: int):
        for (ci, p, kind), sites in entry.sites.items():
            if p != 0:
                continue
            for (cj, q, other), s1 in entry.sites.items():
                if cj == ci and q == 1:
                    self.X[(ci, kind, other)] += sign * len(sites) * len(s1)

    def _add(self, entry: _Entry):
        cid = entry.cid
        self.entries[cid] = entry
        for ci, (allowed, leak) in entry.uni.items():
            if allowed:
                self.uni_trees[ci].add(cid, len(allowed))
            self.uni_allowed[ci] += len(allowed)
            self.uni_leak[ci] += leak
        for (ci, p, kind), sites in entry.sites.items():
            self.S[(ci, p, kind)] += len(sites)
            self._ensure_kind(ci, p, kind)
        self._x_terms(entry, +1)
        for (ci, k0, k1), sampler in self.channels.items():
            s0 = entry.sites.get((ci, 0, k0), ())
            s1 = entry.sites.get((ci, 1, k1), ())
            if s0 or s1:
                sampler.set(cid, s0, s1)
        met = {k: len(v) for k, v in entry.sites.items() if self.classes[k[0]].met_forward}
        if met:
            self.met_comps[cid] = met

    def _remove(self, cid: int):
        entry = self.entries.pop(cid)
        for ci, (allowed, leak) in entry.uni.items():
            if allowed:
                self.uni_trees[ci].add(cid, -len(allowed))
            self.uni_allowed[ci] -= len(allowed)
            self.uni_leak[ci] -= leak
        for key, sites in entry.sites.items():
            self.S[key] -= len(sites)
        self._x_terms(entry, -1)
        for (ci, k0, k1), sampler in self.channels.items():
            if entry.sites.get((ci, 0, k0)) or entry.sites.get((ci, 1, k1)):
                sampler.set(cid, (), ())
        self.met_comps.pop(cid, None)
        for pair in [p for p in self.pairs if cid in p]:
            del self.pairs[pair]

    # -- MET ------------------------------------------------------------------
    def pair_rate(self, e0: _Entry, e1: _Entry) -> Optional[PairRate]:
        """MET pair rate of two components (None when no MET class pairs them)."""
        c0 = {k: len(v) for k, v in e0.sites.items() if self.classes[k[0]].met_forward}
        c1 = {k: len(v) for k, v in e1.sites.items() if self.classes[k[0]].met_forward}
        if not c0 or not c1:
            return None
        per = {}
        for x, y in ((c0, c1), (c1, c0)):
            for (ci, p, k0), n0 in x.items():
                if p != 0:
                    continue
                cls = self.classes[ci]
                for (cj, q, k1), n1 in y.items():
                    if cj != ci or q != 1:
                        continue
                    slot = per.setdefault(ci, [0, 0])
                    slot[1 if _pair_leak(cls, k0, k1) else 0] += n0 * n1
        if not per:
            return None
        classes = tuple((ci, a / self.classes[ci].divisor, l / self.classes[ci].divisor)
                        for ci, (a, l) in sorted(per.items()))
        k_allowed = math.fsum(self.weights[ci] * a for ci, a, _l in classes)
        k_leak = math.fsum(self.weights[ci] * l for ci, _a, l in classes)
        k_act = k_allowed + k_leak
        if k_act <= 0.0:
            return None
        k_d = self.transport.k_diffusion(self.temperature, e0.met_length, e1.met_length)
        bypass = e0.species is not None or e1.species is not None
        phi = 1.0 if bypass else escape_factor(k_act, k_d)
        return PairRate(k_act, k_allowed, k_d, phi, bypass, phi * k_allowed / self.nav,
                        phi * k_leak / self.nav, classes)

    def _refresh_met(self, changed: Optional[Iterable[int]] = None):
        comps = sorted(self.met_comps)
        if changed is None:
            self.pairs = {}
            todo = [(a, b) for i, a in enumerate(comps) for b in comps[i + 1:]]
        else:
            changed = sorted(set(changed) & set(comps))
            todo = sorted({tuple(sorted((a, b))) for a in changed for b in comps if a != b})
        for a, b in todo:
            rate = self.pair_rate(self.entries[a], self.entries[b])
            if rate is not None:
                self.pairs[(a, b)] = rate
            else:
                self.pairs.pop((a, b), None)

    def met_pair_rates(self) -> Dict[Tuple[int, int], PairRate]:
        return dict(self.pairs)

    # -- propensities ---------------------------------------------------------
    def propensities(self) -> Propensities:
        channels, leak = [], []
        for cls in self.classes:
            w = self.weights[cls.index]
            if not cls.bimolecular:
                channels.append((("uni", cls.index), w * self.uni_allowed[cls.index]))
                leak.append(w * self.uni_leak[cls.index])
            elif not cls.met_forward:
                for (ci, k0, k1) in sorted((k for k in self.channels if k[0] == cls.index), key=repr):
                    pairs = (self.S[(ci, 0, k0)] * self.S[(ci, 1, k1)] - self.X[(ci, k0, k1)])
                    a = w * pairs / (cls.divisor * self.nav)
                    if _pair_leak(cls, k0, k1):
                        leak.append(a)
                    else:
                        channels.append((("bi", ci, k0, k1), a))
        met = math.fsum(r.allowed for r in self.pairs.values())
        channels.append((("met",), met))
        leak.append(math.fsum(r.leak for r in self.pairs.values()))
        return Propensities(math.fsum(a for _k, a in channels), math.fsum(leak), tuple(channels))

    # -- snapshots (incremental == full recompute) ---------------------------
    def snapshot(self):
        """Canonical image of every index layer (entries, aggregates, sampler
        totals, MET pair table, propensities)."""
        entries = tuple(self.entries[c].canonical() for c in sorted(self.entries))
        s = tuple(sorted((k, v) for k, v in self.S.items() if v))
        x = tuple(sorted((k, v) for k, v in self.X.items() if v))
        samplers = tuple(sorted((k, (v.n0, v.n1, v.x, v.ordered_total))
                                for k, v in self.channels.items() if v.n0 and v.n1))
        uni = tuple(sorted((k, self.uni_allowed[k], self.uni_leak[k], self.uni_trees[k].prefix(self._size))
                           for k in self.uni_trees))
        pairs = tuple(sorted((k, (v.classes, v.k_act, v.k_d, v.phi, v.bypass, v.allowed, v.leak))
                             for k, v in self.pairs.items()))
        prop = self.propensities()
        return {"entries": entries, "S": s, "X": x, "samplers": samplers, "uni": uni,
                "pairs": pairs, "enabled": prop.enabled, "leak": prop.leak}

    def fresh_snapshot(self):
        """The snapshot a full recompute gives on the current state (the
        engine's own index is left untouched)."""
        other = object.__new__(MeltEngine)
        other.__dict__.update({k: v for k, v in self.__dict__.items()
                               if k in ("spec", "transport", "catalog", "rng", "state", "volume", "nav",
                                        "identity", "temperature", "templates", "classes", "weights",
                                        "quotient")})
        other._index_revision = 0
        other.rebuild()
        return other.snapshot()

    def check_consistency(self):
        """Sampler totals equal the S/X aggregates; unimolecular trees equal counts."""
        for (ci, k0, k1), sampler in self.channels.items():
            pairs = self.S[(ci, 0, k0)] * self.S[(ci, 1, k1)] - self.X[(ci, k0, k1)]
            if sampler.ordered_total != pairs:
                raise EngineError(f"channel {(ci, k0, k1)}: sampler {sampler.ordered_total} != aggregates {pairs}")
        for ci, tree in self.uni_trees.items():
            if tree.prefix(self._size) != self.uni_allowed[ci]:
                raise EngineError(f"template {ci}: unimolecular tree != count")

    # -- building and applying transitions -----------------------------------
    def _pair_transition(self, ci, c0, s0, c1, s1) -> Tuple[EngineTransition, Any]:
        t = self.templates[ci]
        view = _LocalView(self.state, sorted({c0, c1}))
        e0 = [e for e in participant_embeddings(view, t, 0, quotient=self.quotient)
              if e.cid == c0 and tuple(e.mapping[x] for x in t.edited[0]) == s0]
        e1 = [e for e in participant_embeddings(view, t, 1, quotient=self.quotient)
              if e.cid == c1 and tuple(e.mapping[x] for x in t.edited[1]) == s1]
        comps = tuple(sorted((c0, c1)))
        found = {}
        for a in e0:
            for b in e1:
                mapping = tuple(x if x is not None else y for x, y in zip(a.mapping, b.mapping))
                for ca in a.choices:
                    for cb in b.choices:
                        edit_key, identity, edits, delta, fb, fr = wm._evaluate(
                            view, t, ci, mapping, comps, ca + cb)
                        prev = found.get(identity)
                        if prev is None or edit_key < prev.edit_key:
                            found[identity] = EngineTransition(t.name, t.family, ci, edit_key, comps,
                                                               edits, fb, fr, delta,
                                                               _creates_branch(self.state, edits))
        if len(found) != 1:
            raise EngineError(f"{t.name}: site pair {s0} x {s1} gives {len(found)} outcomes, expected 1")
        return next(iter(found.values())), view

    def _begin_apply(self, tr: EngineTransition, view):
        """Apply ``tr`` inside an open transaction (returned): matched form,
        edits, first-form normalization, product check."""
        state = self.state
        predicted = {a: view.base_row(a) for c in tr.components for a in state.components.atoms(c)}
        predicted.update(dict(tr.delta))
        tx = state.begin()
        try:
            switch = []
            for a, b, o in tr.form_bonds:
                cur = float(state.atoms[a].edges[b])
                if cur != o:
                    switch.append(SetBond(a, b, cur, o))
            for a, r in tr.form_radicals:
                cur = state.atoms[a].radical
                if cur != r:
                    switch.append(SetRadical(a, cur, r))
            if switch:
                tx.apply(switch)
            tx.apply(tr.edits)
            cids = sorted({state.components.component_of(a) for a in predicted})
            normalize = [e for c in cids for e in wm._normalizing_edits(state, c)]
            if normalize:
                tx.apply(normalize)
            for c in cids:
                if wm._normalizing_edits(state, c):
                    raise EngineError(f"component {c} is not at its normalization fixed point")
            for a in sorted(predicted):
                n = state.atoms[a]
                row = (n.element, n.isotope, n.charge, n.radical, n.lone_pairs,
                       tuple(sorted((b, float(o)) for b, o in n.edges.items())))
                if row != predicted[a]:
                    raise EngineError(f"atom {a}: applied product differs from the matched one")
            return tx, cids
        except BaseException:
            if tx.open:
                tx.rollback()
            raise

    def _conserved(self, journal):
        old_comp, new_comp = Counter(), Counter()
        old_bb = new_bb = 0
        for rec in journal._old_records.values():
            old_comp.update(dict(rec.composition))
            old_bb += rec.backbone
        comps = self.state.components
        for cid, _atoms in journal.new_components:
            new_comp.update(comps.composition(cid))
            new_bb += comps.backbone_count(cid)
        if old_comp != new_comp or old_bb != new_bb:
            raise EngineError(f"conservation broken: {old_comp} -> {new_comp}, backbone {old_bb} -> {new_bb}")
        self.conservation_checks += 1

    def _fire(self, tr: EngineTransition, view, channel: str, thinning: bool) -> Optional[Event]:
        if tr.leak:
            raise EngineError(f"{tr.template}: a branching (leak) transition reached the fire step")
        tx, products = self._begin_apply(tr, view)
        new = None
        if thinning:
            if len(products) == 2:
                new = self._scan(products)
                rate = self.pair_rate(new[products[0]], new[products[1]])
                if rate is None:
                    phi, bypass = 1.0, (new[products[0]].species is not None
                                        or new[products[1]].species is not None)
                else:
                    phi, bypass = rate.phi, rate.bypass
                if bypass:
                    self.met_counts["reverse_bypass"] += 1
                if phi < 1.0 and self.rng.random() >= phi:
                    tx.rollback()
                    self.met_counts["reverse_rejected"] += 1
                    return None
                self.met_counts["reverse_accepted"] += 1
            else:
                self.met_counts["reverse_no_pair"] += 1
        journal = tx.commit()
        self._conserved(journal)
        old = [c for c, _a in journal.old_components]
        added = [c for c, _a in journal.new_components]
        self._update(old, added, new)
        event = Event(self.time, tr.template, tr.family, channel, tr.components, tuple(added),
                      self.state.revision)
        self.events.append(event)
        self.events_fired += 1
        self.family_counts[tr.family] += 1
        self.template_counts[tr.template] += 1
        return event

    def _update(self, old, added, scanned=None):
        for c in old:
            if c in self.entries:
                self._remove(c)
        scanned = dict(scanned or {})
        todo = [c for c in added if c not in scanned]
        scanned.update(self._scan(todo))
        for c in added:
            self._add(scanned[c])
        self._refresh_met(added)
        self._revision = self.state.revision
        self._index_revision += 1

    # -- channel draws ----------------------------------------------------------
    def _draw_uni(self, ci) -> Tuple[EngineTransition, Any]:
        tree = self.uni_trees[ci]
        cid, r = tree.find(int(self.rng.integers(self.uni_allowed[ci])))
        allowed = self.entries[cid].uni[ci][0]
        tr = allowed[r]
        return tr, _LocalView(self.state, [cid])

    def _draw_met(self, total: float) -> Tuple[EngineTransition, Any, PairRate]:
        r = self.rng.random() * total
        pairs = sorted(self.pairs.items())
        chosen = None
        for key, rate in pairs:
            if rate.allowed <= 0.0:
                continue
            chosen = (key, rate)
            if r < rate.allowed:
                break
            r -= rate.allowed
        (a, b), rate = chosen
        weights = [(ci, self.weights[ci] * t_a) for ci, t_a, _l in rate.classes if t_a > 0]
        r = self.rng.random() * math.fsum(w for _c, w in weights)
        ci = weights[-1][0]
        for c, w in weights:
            if r < w:
                ci = c
                break
            r -= w
        cls = self.classes[ci]
        ea, eb = self.entries[a], self.entries[b]
        options = []
        for x, y in ((ea, eb), (eb, ea)):
            for (cj, p, k0), s0 in x.sites.items():
                if cj != ci or p != 0:
                    continue
                for (ck, q, k1), s1 in y.sites.items():
                    if ck == ci and q == 1 and not _pair_leak(cls, k0, k1):
                        options.append((x.cid, s0, y.cid, s1))
        options.sort(key=lambda o: (o[0], o[2], o[1], o[3]))
        total_pairs = sum(len(s0) * len(s1) for _x, s0, _y, s1 in options)
        g = int(self.rng.integers(total_pairs))
        for x, s0, y, s1 in options:
            n = len(s0) * len(s1)
            if g < n:
                left, right = s0[g // len(s1)], s1[g % len(s1)]
                tr, view = self._pair_transition(ci, x, left, y, right)
                return tr, view, rate
            g -= n
        raise EngineError("MET draw fell through")

    # -- clock ----------------------------------------------------------------
    def _token(self):
        return (id(self.state), self.state.revision, self._index_revision, self.temperature,
                self.volume, id(self.catalog), id(self.rng))

    def _sync(self):
        if self.state.revision != self._revision:
            self.rebuild()

    def _integrate(self, elapsed, enabled, leak):
        self.leak_num = self._leak_num_c + leak * elapsed
        self.leak_den = self._leak_den_c + (enabled + leak) * elapsed
        self.time = self._segment_start + elapsed
        self._commit()

    def _commit(self):
        self._segment_start = self.time
        self._leak_num_c = self.leak_num
        self._leak_den_c = self.leak_den
        self._pending = None

    def _view_to_horizon(self, until, enabled, leak):
        elapsed = until - self._segment_start
        self.leak_num = self._leak_num_c + leak * elapsed
        self.leak_den = self._leak_den_c + (enabled + leak) * elapsed
        self.time = self._segment_start + elapsed

    def step(self, *, until: Optional[float] = None) -> Optional[Event]:
        """Advance through null proposals until one event fires or a horizon ends."""
        self._sync()
        while until is None or self.time < until:
            prop = self.propensities()
            enabled, leak = prop.enabled, prop.leak
            pending = self._pending
            if pending is not None and not (
                    self.time == pending.horizon_time
                    and (pending.enabled, pending.leak) == (enabled, leak)
                    and pending.token == self._token()):
                self._commit()
                pending = None
            if pending is None:
                if enabled <= 0.0:
                    pending = _Pending(math.inf, True, enabled, leak, self.time, ())
                else:
                    pending = _Pending(float(self.rng.exponential(1.0 / enabled)), False,
                                       enabled, leak, self.time, ())
            if pending.absorbing or (until is not None and self._segment_start + pending.elapsed > until):
                self._stop = "absorbing" if pending.absorbing else "horizon"
                if pending.absorbing and until is None:
                    return None
                self._view_to_horizon(until, enabled, leak)
                self._pending = _Pending(pending.elapsed, pending.absorbing, enabled, leak,
                                         self.time, self._token())
                return None
            self._integrate(pending.elapsed, enabled, leak)
            r = self.rng.random() * enabled
            key = None
            for k, a in prop.channels:
                if a <= 0.0:
                    continue
                key = k
                if r < a:
                    break
                r -= a
            event = self._fire_channel(key, dict(prop.channels)[key])
            if event is None:
                self.null_proposals += 1
                continue
            return event
        return None

    def _fire_channel(self, key, total) -> Optional[Event]:
        kind = key[0]
        if kind == "uni":
            ci = key[1]
            tr, view = self._draw_uni(ci)
            return self._fire(tr, view, "uni", self.classes[ci].met_reverse)
        if kind == "bi":
            _k, ci, k0, k1 = key
            c0, s0, c1, s1 = self.channels[(ci, k0, k1)].draw(self.rng)
            tr, view = self._pair_transition(ci, c0, s0, c1, s1)
            return self._fire(tr, view, "bi", self.classes[ci].met_reverse)
        tr, view, rate = self._draw_met(total)
        event = self._fire(tr, view, "met", False)
        self.met_counts["forward_bypass" if rate.bypass else "forward"] += 1
        return event

    def run(self, *, until: Optional[float] = None, max_events: Optional[int] = None,
            check=None) -> RunResult:
        """Run to a horizon and/or an event cap.  ``check(engine, event)`` is
        called after every fired event."""
        if until is None and max_events is None:
            raise ValueError("run needs a time or event limit")
        start = self.events_fired
        stop = "horizon"
        while (until is None or self.time < until) and (
                max_events is None or self.events_fired - start < max_events):
            event = self.step(until=until)
            if event is None:
                stop = self._stop
                break
            if check is not None:
                check(self, event)
        else:
            if max_events is not None and self.events_fired - start >= max_events:
                stop = "max_events"
        return RunResult(self.time, self.events_fired - start, self.events_fired, stop,
                         self.leak_num, self.leak_den,
                         self.leak_num / self.leak_den if self.leak_den else 0.0,
                         self.null_proposals, dict(self.family_counts), dict(self.template_counts),
                         dict(self.met_counts))

    # -- report -----------------------------------------------------------------
    def report(self) -> dict:
        comps = self.state.components
        total = Counter()
        backbone = 0
        for cid in comps.ids():
            total.update(comps.composition(cid))
            backbone += comps.backbone_count(cid)
        return {
            "time": self.time, "temperature": self.temperature, "volume_m3": self.volume,
            "events": self.events_fired, "null_proposals": self.null_proposals,
            "leak_num": self.leak_num, "leak_den": self.leak_den,
            "leak_fraction": self.leak_num / self.leak_den if self.leak_den else 0.0,
            "family_counts": dict(self.family_counts), "template_counts": dict(self.template_counts),
            "met": dict(self.met_counts), "conservation_checks": self.conservation_checks,
            "components": len(comps), "composition": {f"{k[0]}{k[1] or ''}": v for k, v in sorted(total.items(), key=str)},
            "backbone_atoms": backbone,
        }


__all__ = [
    "ComponentPairSampler", "DeclarationError", "EngineError", "EngineTransition", "Event",
    "MET_SCHEMA", "MeltEngine", "PairRate", "PayloadCatalog", "Propensities", "RunResult",
    "SpeciesIdentity", "TransportLaw", "escape_factor",
]
