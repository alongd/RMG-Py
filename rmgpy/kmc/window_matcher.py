"""Window matcher: hand-written templates + an atom-native MeltState -> transitions.

A template is an executable window: atoms (element, isotope, charge, radical
electrons, lone pairs), the induced bonds among them, and for each atom whether
it is ``closed`` (its full degree lies inside the window) or a ``port`` (its
other neighbours are described by an ``exterior`` multiset of (bond order,
element, in_ring)).  Edits are final assignments on closed atoms only.  Atoms
are split into one or two participants (molecules); a participant's window may
be disconnected, and path predicates bound the shortest heavy-atom path between
two of its atoms, so intramolecular H shifts that run through atoms outside the
window are expressible.  ``role_swap`` (two participants exchangeable by an
automorphism that fixes the edit set) is computed here, never supplied.

Matching is pure and deterministic: ``match(templates, state)``.

* Embeddings are injective, induced and signature-equal; closed atoms have
  exact degree; port exteriors are equal as multisets.  Every embedding is
  enumerated (no one-per-anchor shortcut), across every supported resonance
  form of the component.
* Resonance forms are enumerated per complete component, factorized over its
  conjugated subsystems (atoms with free valence, joined by bonds).  Forms move
  radicals and pi bonds inside a subsystem (allyl shifts, alternating-cycle
  flips); charges and lone pairs never move.  A six-cycle whose pi bonds
  alternate is written as an aromatic sextet (all 1.5).  At most
  ``RESONANCE_FORM_LIMIT`` forms per component (checked while exploring) and no
  conjugated ring larger than ``CONJUGATED_RING_LIMIT``; either overflow raises.
  A window touches only the subsystems that contain its atoms; the forms of the
  other subsystems multiply the embedding count but cannot change the outcome,
  so they are counted, not re-evaluated.
* A transition's identity is (template class, the touched components, the
  complete mapped product of those components with every conjugated subsystem
  written in its first form).  It is stored as the difference from the
  reactant components written in their first forms, which is equivalent
  because both sides are fixed points of the same normalization.  Embeddings
  with equal identity in one class are one transition; equal outcomes from
  different classes stay separate pathways.
* Two-participant classes are counted by the participant quotient:
  ``n_i(c)`` edit-role-distinct sites of participant i in component c,
  ``P = n0*n1 - sum_c n0(c)*n1(c)``, halved under ``role_swap`` (and divided
  by the order of any participant-preserving symmetry that moves edited
  atoms; trivial for ordinary windows).
  ``PairSampler`` draws such pairs exactly (Fenwick tree over components).
* ``apply_transition`` moves the touched subsystems to the matched form,
  applies the edits and normalizes the product components to their first form
  inside one MeltState transaction; any failure rolls the whole transaction
  back.

Rates are not part of this module: every template carries a caller-given
weight.
"""

import math
import weakref
from dataclasses import dataclass
from typing import Dict, Mapping, Optional, Sequence, Tuple

from rmgpy.kmc.melt_state import BOND_ORDERS, SetBond, SetCharge, SetLonePairs, SetRadical
from rmgpy.molecule.element import PeriodicSystem

RESONANCE_FORM_LIMIT = 64
CONJUGATED_RING_LIMIT = 8
RAW_PLACEMENT_LIMIT = 4096  # integer pi placements explored per subsystem


class WindowMatcherError(ValueError):
    """Base class of every refusal of this module."""


class TemplateError(WindowMatcherError):
    """A template is malformed."""


class ResonanceLimitError(WindowMatcherError):
    """A component has more resonance forms than the bound allows."""

    def __init__(self, message, count=None):
        super().__init__(message)
        self.count = count


class ConjugatedRingError(ResonanceLimitError):
    """A conjugated subsystem contains a ring larger than the supported size."""


class StaleTransitionError(WindowMatcherError):
    """The state changed since the transition was matched."""


# --------------------------------------------------------------------------
# chemistry helpers
# --------------------------------------------------------------------------
def _valence(element):
    try:
        return PeriodicSystem.valence_electrons[element]
    except KeyError:
        raise WindowMatcherError(f"element {element!r} has no valence data") from None


def _key(a, b):
    return (a, b) if a < b else (b, a)


# --------------------------------------------------------------------------
# template
# --------------------------------------------------------------------------
def _int(value, name, minimum=None):
    if isinstance(value, bool) or not isinstance(value, int):
        raise TemplateError(f"{name} must be an integer, got {value!r}")
    if minimum is not None and value < minimum:
        raise TemplateError(f"{name} must be >= {minimum}, got {value!r}")
    return value


@dataclass(frozen=True)
class WindowAtom:
    """One window atom.  ``exterior`` (ports only) is the multiset of
    (bond order, element, in_ring) of the neighbours outside the window."""

    element: str
    participant: int = 0
    kind: str = "closed"
    exterior: Tuple = ()
    isotope: Optional[int] = None
    charge: int = 0
    radical: int = 0
    lone_pairs: int = 0

    def __post_init__(self):
        if not isinstance(self.element, str) or not self.element:
            raise TemplateError("a window atom needs an element symbol")
        if self.kind not in ("closed", "port"):
            raise TemplateError(f"kind must be closed or port, got {self.kind!r}")
        if isinstance(self.participant, bool) or self.participant not in (0, 1):
            raise TemplateError(f"participant must be 0 or 1, got {self.participant!r}")
        _int(self.charge, "charge")
        _int(self.radical, "radical", 0)
        _int(self.lone_pairs, "lone_pairs", 0)
        if self.isotope is not None:
            _int(self.isotope, "isotope", 1)
        try:
            exterior = tuple(sorted((float(o), str(e), bool(r)) for o, e, r in (self.exterior or ())))
        except (TypeError, ValueError):
            raise TemplateError(f"exterior must hold (order, element, in_ring) triples, got "
                                f"{self.exterior!r}") from None
        for order, _e, _r in exterior:
            if order not in BOND_ORDERS:
                raise TemplateError(f"exterior bond order {order!r} is not one of {BOND_ORDERS}")
        object.__setattr__(self, "exterior", exterior)


@dataclass(frozen=True)
class BondEdit:
    """Final bond order between two closed window atoms (0 = no bond)."""

    a: int
    b: int
    order: float


@dataclass(frozen=True)
class AtomEdit:
    """Final radical / charge / lone pairs of a closed window atom (None = unchanged)."""

    atom: int
    radical: Optional[int] = None
    charge: Optional[int] = None
    lone_pairs: Optional[int] = None


@dataclass(frozen=True)
class PathPredicate:
    """Atoms ``a`` and ``b`` of one participant lie in one component and their
    shortest heavy-atom path has between ``minimum`` and ``maximum`` bonds."""

    a: int
    b: int
    minimum: int
    maximum: int


@dataclass(frozen=True)
class _Piece:
    plan: Tuple  # ((template atom, parent or None, ((earlier atom, bonded?), ...)), ...)


class Template:
    """A hand-written executable template (see the module docstring).

    ``edited`` holds, per participant, the sorted indices of the atoms named by
    an edit; a participant's site is the image of those atoms.
    """

    def __init__(self, name, family, atoms, bonds, edits, *, weight, paths=()):
        if not isinstance(name, str) or not name:
            raise TemplateError("a template needs a name")
        self.name = name
        self.family = str(family)
        self.atoms = tuple(atoms)
        if not self.atoms or not all(isinstance(a, WindowAtom) for a in self.atoms):
            raise TemplateError(f"{name}: atoms must be a non-empty sequence of WindowAtom")
        try:
            weight = float(weight)
        except (TypeError, ValueError):
            raise TemplateError(f"{name}: weight must be a number") from None
        if not math.isfinite(weight) or weight <= 0:
            raise TemplateError(f"{name}: weight must be finite and positive, got {weight!r}")
        self.weight = weight
        n = len(self.atoms)
        self.participants = tuple(sorted({a.participant for a in self.atoms}))
        if self.participants not in ((0,), (0, 1)):
            raise TemplateError(f"{name}: participants must be 0 or 0 and 1")
        for i, atom in enumerate(self.atoms):
            if atom.kind == "port" and not atom.exterior:
                raise TemplateError(f"{name}: port {i} needs a non-empty exterior")
            if atom.kind == "closed" and atom.exterior:
                raise TemplateError(f"{name}: closed atom {i} cannot have an exterior")
            _valence(atom.element)
        # bonds
        self.bond: Dict[Tuple[int, int], float] = {}
        for item in bonds:
            i, j, order = item
            self._index(i), self._index(j)
            if i == j:
                raise TemplateError(f"{name}: bond {item!r} needs two atoms")
            if float(order) not in BOND_ORDERS:
                raise TemplateError(f"{name}: bond order {order!r} is not one of {BOND_ORDERS}")
            if self.atoms[i].participant != self.atoms[j].participant:
                raise TemplateError(f"{name}: bond {item!r} joins two participants")
            if _key(i, j) in self.bond:
                raise TemplateError(f"{name}: bond {item!r} is given twice")
            self.bond[_key(i, j)] = float(order)
        self.bond = dict(sorted(self.bond.items()))
        self.bonds = tuple((i, j, o) for (i, j), o in self.bond.items())
        neighbours = [[] for _ in range(n)]
        for i, j in self.bond:
            neighbours[i].append(j)
            neighbours[j].append(i)
        self.neighbours = tuple(tuple(sorted(x)) for x in neighbours)
        # edits
        self.edits = tuple(edits)
        atom_final, bond_final, seen = {}, {}, set()
        for edit in self.edits:
            if isinstance(edit, AtomEdit):
                i = self._index(edit.atom)
                self._closed(i, edit)
                if ("atom", i) in seen:
                    raise TemplateError(f"{name}: atom {i} is edited twice")
                seen.add(("atom", i))
                final = {}
                for field in ("radical", "charge", "lone_pairs"):
                    value = getattr(edit, field)
                    if value is None:
                        continue
                    _int(value, field, None if field == "charge" else 0)
                    if value == getattr(self.atoms[i], field):
                        raise TemplateError(f"{name}: {edit} is a no-op edit")
                    final[field] = value
                if not final:
                    raise TemplateError(f"{name}: {edit} changes nothing (no-op edit)")
                atom_final[i] = final
            elif isinstance(edit, BondEdit):
                i, j = self._index(edit.a), self._index(edit.b)
                if i == j:
                    raise TemplateError(f"{name}: {edit} needs two atoms")
                self._closed(i, edit)
                self._closed(j, edit)
                order = float(edit.order)
                if order != 0 and order not in BOND_ORDERS:
                    raise TemplateError(f"{name}: {edit} order is not 0 or one of {BOND_ORDERS}")
                if ("bond", _key(i, j)) in seen:
                    raise TemplateError(f"{name}: bond {i}-{j} is edited twice")
                seen.add(("bond", _key(i, j)))
                if order == self.bond.get(_key(i, j), 0.0):
                    raise TemplateError(f"{name}: {edit} is a no-op edit")
                bond_final[_key(i, j)] = order
            else:
                raise TemplateError(f"{name}: unsupported edit {edit!r}")
        if not any((self.bond.get(k, 0.0) == 0) != (o == 0) for k, o in bond_final.items()):
            raise TemplateError(f"{name}: the connectivity change is empty (a pure pi move)")
        self.atom_final = atom_final
        self.bond_final = bond_final
        self._check_valence()
        edited = sorted(set(atom_final) | {i for k in bond_final for i in k})
        self.edited_all = tuple(edited)
        self.edited = tuple(tuple(i for i in edited if self.atoms[i].participant == p)
                            for p in self.participants)
        if any(not e for e in self.edited):
            raise TemplateError(f"{name}: every participant needs an edited atom")
        # paths
        self.paths = tuple(paths)
        for path in self.paths:
            if not isinstance(path, PathPredicate):
                raise TemplateError(f"{name}: paths must be PathPredicate")
            a, b = self._index(path.a), self._index(path.b)
            _int(path.minimum, "path minimum", 1)
            _int(path.maximum, "path maximum", path.minimum)
            if a == b or self.atoms[a].participant != self.atoms[b].participant:
                raise TemplateError(f"{name}: a path joins two atoms of one participant")
        # derived matching data
        self.sig = tuple((a.element, a.isotope, a.charge, a.lone_pairs,
                          len(self.neighbours[i]) + len(a.exterior))
                         for i, a in enumerate(self.atoms))
        self.ext_sigma = tuple(tuple(sorted((e, r) for _o, e, r in a.exterior)) for a in self.atoms)
        self.pieces = tuple(self._pieces(p) for p in self.participants)
        # symmetry of the edited atoms (site tuple = edited[0] + edited[1])
        self.site_order = tuple(t for e in self.edited for t in e)
        position = {t: k for k, t in enumerate(self.site_order)}
        preserving = self._restrictions(swap=False)
        swapping = self._restrictions(swap=True) if len(self.participants) == 2 else ()
        self.edit_symmetry = tuple(tuple(position[t] for t in r) for r in preserving)
        self.role_swap = bool(swapping)
        self.swap_positions = (tuple(position[t] for t in swapping[0])
                               if swapping else None)

    def __repr__(self):
        return f"Template({self.name!r}, {self.family!r}, {len(self.atoms)} atoms)"

    # -- validation helpers ------------------------------------------------
    def _index(self, i):
        if isinstance(i, bool) or not isinstance(i, int) or not 0 <= i < len(self.atoms):
            raise TemplateError(f"{self.name}: atom index {i!r} is out of range")
        return i

    def _closed(self, i, edit):
        if self.atoms[i].kind != "closed":
            raise TemplateError(f"{self.name}: {edit} touches port {i}; edits never touch a port")

    def _check_valence(self):
        for stage in ("before", "after"):
            orders = dict(self.bond)
            if stage == "after":
                for k, o in self.bond_final.items():
                    if o == 0:
                        orders.pop(k, None)
                    else:
                        orders[k] = o
            for i, atom in enumerate(self.atoms):
                final = self.atom_final.get(i, {}) if stage == "after" else {}
                radical = final.get("radical", atom.radical)
                charge = final.get("charge", atom.charge)
                lone_pairs = final.get("lone_pairs", atom.lone_pairs)
                used = (sum(o for k, o in orders.items() if i in k)
                        + sum(o for o, _e, _r in atom.exterior) + radical + 2 * lone_pairs)
                if used != _valence(atom.element) - charge:
                    raise TemplateError(
                        f"{self.name}: atom {i} ({atom.element}) violates valence {stage} the edits")

    def _pieces(self, participant):
        members = [i for i, a in enumerate(self.atoms) if a.participant == participant]
        edited = set(self.edited_all)
        seen, pieces = set(), []
        for start in members:
            if start in seen:
                continue
            component, stack = [], [start]
            seen.add(start)
            while stack:
                i = stack.pop()
                component.append(i)
                for j in self.neighbours[i]:
                    if j not in seen:
                        seen.add(j)
                        stack.append(j)
            component.sort()
            root = min(component, key=lambda i: (
                i not in edited, self.atoms[i].element == "H", self.atoms[i].kind != "closed",
                -len(self.neighbours[i]), i))
            order, parent, placed = [root], {root: None}, {root}
            k = 0
            while k < len(order):
                for j in self.neighbours[order[k]]:
                    if j not in placed:
                        placed.add(j)
                        parent[j] = order[k]
                        order.append(j)
                k += 1
            plan = tuple((t, parent[t], tuple((u, _key(t, u) in self.bond) for u in order[:pos]))
                         for pos, t in enumerate(order))
            pieces.append(_Piece(plan))
        return tuple(pieces)

    def _restrictions(self, swap):
        """Images of the edited atoms (in ``site_order``) under every template
        automorphism that fixes atom labels, bonds, the edit set and the path
        predicates, and keeps (``swap=False``) or exchanges (``swap=True``) the
        participants.  Sorted; the identity comes first when ``swap=False``."""
        n = len(self.atoms)
        label = []
        for i, a in enumerate(self.atoms):
            final = self.atom_final.get(i, {})
            label.append((a.element, a.isotope, a.charge, a.radical, a.lone_pairs, a.kind,
                          a.exterior, tuple(sorted(final.items()))))
        pair = {}
        for k, o in self.bond.items():
            pair.setdefault(k, [None, None, None])[0] = o
        for k, o in self.bond_final.items():
            pair.setdefault(k, [None, None, None])[1] = o
        for p in self.paths:
            pair.setdefault(_key(p.a, p.b), [None, None, None])[2] = (p.minimum, p.maximum)
        pair = {k: tuple(v) for k, v in pair.items()}
        edited = set(self.site_order)
        order = list(self.site_order) + [i for i in range(n) if i not in edited]
        depth = len(self.site_order)
        perm, used, found = {}, set(), set()

        def fits(i, j):
            if j in used or label[j] != label[i]:
                return False
            if (self.atoms[j].participant != self.atoms[i].participant) != swap:
                return False
            return all(pair.get(_key(i, u)) == pair.get(_key(j, v)) for u, v in perm.items())

        def complete(k):
            if k == n:
                return True
            i = order[k]
            for j in range(n):
                if fits(i, j):
                    perm[i] = j
                    used.add(j)
                    done = complete(k + 1)
                    del perm[i]
                    used.discard(j)
                    if done:
                        return True
            return False

        def assign(k):
            if k == depth:
                image = tuple(perm[t] for t in self.site_order)
                if image not in found and complete(k):
                    found.add(image)
                return
            i = order[k]
            for j in range(n):
                if fits(i, j):
                    perm[i] = j
                    used.add(j)
                    assign(k + 1)
                    del perm[i]
                    used.discard(j)

        assign(0)
        return tuple(sorted(found))


# --------------------------------------------------------------------------
# resonance forms of one conjugated subsystem
# --------------------------------------------------------------------------
_SHAPE_CACHE: Dict[Tuple, Tuple] = {}


def _ring_check_and_cycles(n, bonds, incident):
    """Refuse a conjugated ring larger than the limit; return the simple cycles
    of length <= the limit as tuples of bond indices in ring order."""
    for k, (i, j) in enumerate(bonds):
        dist = {i: 0}
        frontier = [i]
        while frontier and j not in dist:
            following = []
            for v in frontier:
                for kk, w in incident[v]:
                    if kk != k and w not in dist:
                        dist[w] = dist[v] + 1
                        following.append(w)
            frontier = following
        if j in dist and dist[j] + 1 > CONJUGATED_RING_LIMIT:
            raise ConjugatedRingError(
                f"conjugated ring of {dist[j] + 1} atoms exceeds {CONJUGATED_RING_LIMIT}")
    cycles = []
    for s in range(n):
        stack = [(s, (s,), ())]
        while stack:
            v, path, pbonds = stack.pop()
            for k, w in incident[v]:
                if w == s and len(path) >= 3 and k != pbonds[-1]:
                    if path[1] < path[-1]:
                        cycles.append(pbonds + (k,))
                elif w > s and w not in path and len(path) < CONJUGATED_RING_LIMIT:
                    stack.append((w, path + (w,), pbonds + (k,)))
    cycles.sort()
    return tuple(cycles)


def _kekulize(free, bonds, orders, radicals, incident):
    for o in orders:
        if o not in (1.0, 1.5, 2.0, 3.0):
            raise WindowMatcherError(f"bond order {o!r} inside a conjugated subsystem")
    pi = [None if o == 1.5 else int(o) - 1 for o in orders]
    n = len(free)
    need = [free[a] - radicals[a] - sum(pi[k] for k, _ in incident[a] if pi[k] is not None)
            for a in range(n)]
    left = [sum(1 for k, _ in incident[a] if pi[k] is None) for a in range(n)]
    unknown = [k for k, p in enumerate(pi) if p is None]
    if any(x < 0 for x in need):
        raise WindowMatcherError("the stored form violates the free valence of a subsystem")

    def visit(pos):
        if pos == len(unknown):
            return all(x == 0 for x in need)
        k = unknown[pos]
        a, b = bonds[k]
        left[a] -= 1
        left[b] -= 1
        for value in (1, 0):
            if need[a] < value or need[b] < value:
                continue
            need[a] -= value
            need[b] -= value
            if need[a] <= left[a] and need[b] <= left[b] and visit(pos + 1):
                pi[k] = value
                return True
            need[a] += value
            need[b] += value
        left[a] += 1
        left[b] += 1
        return False

    if not visit(0):
        raise WindowMatcherError("aromatic bonds of a subsystem have no Kekule assignment")
    return tuple(pi)


def _render(pi, radicals, six_cycles):
    orders = [1.0 + p for p in pi]
    for cycle in six_cycles:
        values = [pi[k] for k in cycle]
        if all(v in (0, 1) for v in values) and all(values[i] != values[i - 1] for i in range(6)):
            for k in cycle:
                orders[k] = 1.5
    return tuple(orders), tuple(radicals)


def _form_order(form):
    """First form: most aromatic sextets, then id-ordered bond orders, radicals."""
    orders, radicals = form
    return (-sum(1 for o in orders if o == 1.5), orders, radicals)


def _moves(pi, rad, incident, even_cycles):
    """Allyl shifts of a radical and flips of alternating even cycles."""
    for a in range(len(rad)):
        if rad[a] == 0:
            continue
        for k, b in incident[a]:
            if pi[k] >= 2:
                continue
            for k2, c in incident[b]:
                if k2 == k or pi[k2] == 0:
                    continue
                p = list(pi)
                r = list(rad)
                p[k] += 1
                p[k2] -= 1
                r[a] -= 1
                r[c] += 1
                yield tuple(p), tuple(r)
    for cycle in even_cycles:
        values = [pi[k] for k in cycle]
        if all(v in (0, 1) for v in values) and all(values[i] != values[i - 1]
                                                    for i in range(len(values))):
            p = list(pi)
            for k in cycle:
                p[k] = 1 - p[k]
            yield tuple(p), rad


def _shape_forms(shape, budget):
    """Rendered forms of one subsystem shape, first form first.  Raises
    ResonanceLimitError as soon as more than ``budget`` forms are found."""
    cached = _SHAPE_CACHE.get(shape)
    if cached is not None:
        if len(cached) > budget:
            raise ResonanceLimitError(f"more than {budget} resonance forms", budget + 1)
        return cached
    free, bonds, orders, radicals = shape
    n = len(free)
    incident = [[] for _ in range(n)]
    for k, (i, j) in enumerate(bonds):
        incident[i].append((k, j))
        incident[j].append((k, i))
    cycles = _ring_check_and_cycles(n, bonds, incident)
    six = tuple(c for c in cycles if len(c) == 6)
    even = tuple(c for c in cycles if len(c) % 2 == 0)
    start = (_kekulize(free, bonds, orders, radicals, incident), tuple(radicals))
    seen = {start}
    queue = [start]
    rendered = {}
    head = 0
    while head < len(queue):
        pi, rad = queue[head]
        head += 1
        form = _render(pi, rad, six)
        if form not in rendered:
            rendered[form] = len(rendered)
            if len(rendered) > budget:
                raise ResonanceLimitError(f"more than {budget} resonance forms", budget + 1)
        for nxt in _moves(pi, rad, incident, even):
            if nxt not in seen:
                if len(seen) >= RAW_PLACEMENT_LIMIT:
                    raise ResonanceLimitError(
                        f"more than {RAW_PLACEMENT_LIMIT} pi placements in one subsystem")
                seen.add(nxt)
                queue.append(nxt)
    forms = tuple(sorted(rendered, key=_form_order))
    if len(_SHAPE_CACHE) > 100000:
        _SHAPE_CACHE.clear()
    _SHAPE_CACHE[shape] = forms
    return forms


@dataclass(frozen=True)
class _Sub:
    atoms: Tuple[int, ...]
    bonds: Tuple[Tuple[int, int], ...]
    apos: Mapping
    bpos: Mapping
    forms: Tuple


def _subsystems(members, value, outside=None):
    """Conjugated subsystems among ``members`` (sorted atom ids).

    ``value(a)`` -> (element, charge, lone pairs, radical, {neighbour: order});
    ``outside(b)`` gives the free valence of a neighbour outside ``members``
    (it must be zero: a subsystem never leaves the region)."""
    free = {}
    for a in members:
        element, charge, lone_pairs, _radical, edges = value(a)
        free[a] = _valence(element) - charge - 2 * lone_pairs - len(edges)
        if free[a] < 0:
            raise WindowMatcherError(f"atom {a} has negative free valence")
    out, seen = [], set()
    for start in members:
        if free[start] <= 0 or start in seen:
            continue
        seen.add(start)
        stack, group = [start], [start]
        while stack:
            a = stack.pop()
            for b in sorted(value(a)[4]):
                if b not in free:
                    if outside is None or outside(b) > 0:
                        raise WindowMatcherError(f"conjugation leaves the region at atom {b}")
                    continue
                if free[b] > 0 and b not in seen:
                    seen.add(b)
                    stack.append(b)
                    group.append(b)
        group.sort()
        out.append((tuple(group), free))
    return out


def _sub_shape(group, free, value):
    index = {a: k for k, a in enumerate(group)}
    bonds = tuple(sorted((a, b) for a in group for b in value(a)[4] if b in index and a < b))
    shape = (tuple(free[a] for a in group),
             tuple((index[a], index[b]) for a, b in bonds),
             tuple(float(value(a)[4][b]) for a, b in bonds),
             tuple(value(a)[3] for a in group))
    return bonds, shape


def _component_subs(members, value, limit=RESONANCE_FORM_LIMIT):
    """Subsystems of one component with their forms, enforcing the
    per-component bound while exploring."""
    subs, product = [], 1
    for group, free in _subsystems(members, value):
        bonds, shape = _sub_shape(group, free, value)
        budget = limit // product
        try:
            forms = _shape_forms(shape, budget)
        except ConjugatedRingError:
            raise
        except ResonanceLimitError as error:
            count = product * (budget + 1) if error.count is not None else None
            raise ResonanceLimitError(
                f"component at atom {members[0]}: more than {limit} resonance forms "
                f"(stopped at form {count})", count) from None
        product *= len(forms)
        subs.append(_Sub(group, bonds, {a: k for k, a in enumerate(group)},
                         {b: k for k, b in enumerate(bonds)}, forms))
    return subs, product


# --------------------------------------------------------------------------
# state view
# --------------------------------------------------------------------------
def _ring_atoms(adj, ids):
    """Atoms on a cycle (incident to a non-bridge edge)."""
    disc, low, bridges = {}, {}, set()
    counter = 0
    for root in ids:
        if root in disc:
            continue
        disc[root] = low[root] = counter
        counter += 1
        stack = [(root, -1, 0)]
        while stack:
            v, parent, pos = stack[-1]
            neighbours = adj[v]
            if pos < len(neighbours):
                stack[-1] = (v, parent, pos + 1)
                w = neighbours[pos]
                if w == parent:
                    continue
                if w in disc:
                    low[v] = min(low[v], disc[w])
                else:
                    disc[w] = low[w] = counter
                    counter += 1
                    stack.append((w, v, 0))
            else:
                stack.pop()
                if stack:
                    u = stack[-1][0]
                    low[u] = min(low[u], low[v])
                    if low[v] > disc[u]:
                        bridges.add(_key(u, v))
    return frozenset(a for a in ids if any(_key(a, b) not in bridges for b in adj[a]))


class _View:
    """Read-only matching data of one state revision."""

    def __init__(self, state):
        self.state = state
        self.revision = state.revision
        nodes = state.atoms
        self.ids = tuple(sorted(nodes))
        self.node = {a: nodes[a] for a in self.ids}
        self.adj = {a: tuple(sorted(nodes[a].edges)) for a in self.ids}
        self.sig, self.free = {}, {}
        by_sig = {}
        for a in self.ids:
            n = self.node[a]
            sig = (n.element, n.isotope, n.charge, n.lone_pairs, len(n.edges))
            self.sig[a] = sig
            self.free[a] = _valence(n.element) - n.charge - 2 * n.lone_pairs - len(n.edges)
            by_sig.setdefault(sig, []).append(a)
        self.by_sig = {k: tuple(v) for k, v in by_sig.items()}
        self.in_ring = _ring_atoms(self.adj, self.ids)
        self.comp = {a: state.components.component_of(a) for a in self.ids}
        self.subs = []
        self.sub_of = {}
        self.comp_forms = {}
        for cid in state.components.ids():
            subs, total = _component_subs(state.components.atoms(cid), self._value)
            self.comp_forms[cid] = total
            for sub in subs:
                for a in sub.atoms:
                    self.sub_of[a] = len(self.subs)
                self.subs.append(sub)
        self._base_rows = {}
        self._dist = {}
        self._hash = None

    def _value(self, a):
        n = self.node[a]
        return n.element, n.charge, n.lone_pairs, n.radical, n.edges

    def full_hash(self):
        if self._hash is None:
            self._hash = self.state.full_hash()
        return self._hash

    def order(self, a, b):
        return float(self.node[a].edges[b])

    def base_row(self, a):
        """Row of atom ``a`` with its subsystem in the first form."""
        row = self._base_rows.get(a)
        if row is None:
            n = self.node[a]
            radical = n.radical
            edges = dict(n.edges)
            s = self.sub_of.get(a)
            if s is not None:
                sub = self.subs[s]
                orders, rads = sub.forms[0]
                radical = rads[sub.apos[a]]
                for b in edges:
                    k = sub.bpos.get(_key(a, b))
                    if k is not None:
                        edges[b] = orders[k]
            row = (n.element, n.isotope, n.charge, radical, n.lone_pairs,
                   tuple(sorted((b, float(o)) for b, o in edges.items())))
            self._base_rows[a] = row
        return row

    def heavy_distance(self, a, b, maximum):
        """Shortest path from a to b through heavy atoms (None beyond maximum)."""
        cached = self._dist.get(a)
        if cached is None or cached[0] < maximum:
            dist = {a: 0} if self.node[a].element != "H" else {}
            frontier = list(dist)
            depth = 0
            while frontier and depth < maximum:
                depth += 1
                following = []
                for v in frontier:
                    for w in self.adj[v]:
                        if w not in dist and self.node[w].element != "H":
                            dist[w] = depth
                            following.append(w)
                frontier = following
            cached = (maximum, dist)
            self._dist[a] = cached
        d = cached[1].get(b)
        return d if d is not None and d <= maximum else None


_VIEWS = weakref.WeakKeyDictionary()


def _view_for(state):
    if getattr(state, "_tx", None) is not None:
        return _View(state)  # an open transaction: never cache a half-applied state
    cached = _VIEWS.get(state)
    if cached is None or cached.revision != state.revision:
        cached = _View(state)
        _VIEWS[state] = cached
    return cached


# --------------------------------------------------------------------------
# public form helpers
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Form:
    """One resonance form of a component: the orders of its subsystem bonds and
    the radicals of its subsystem atoms."""

    bonds: Tuple[Tuple[int, int, float], ...]
    radicals: Tuple[Tuple[int, int], ...]


def component_forms(state, cid) -> Tuple[Form, ...]:
    """Every supported resonance form of a component, first form first.  The
    list does not depend on the form the state stores."""
    members = state.components.atoms(cid)
    nodes = state.atoms

    def value(a):
        n = nodes[a]
        return n.element, n.charge, n.lone_pairs, n.radical, n.edges

    subs, _total = _component_subs(members, value)
    forms = [((), ())]
    for sub in subs:
        forms = [(bonds + tuple((a, b, o) for (a, b), o in zip(sub.bonds, orders)),
                  rads + tuple(zip(sub.atoms, r)))
                 for bonds, rads in forms for orders, r in sub.forms]
    return tuple(Form(tuple(sorted(b)), tuple(sorted(r))) for b, r in forms)


def form_rows(state, cid, form) -> Dict[int, Tuple]:
    """atom -> (element, isotope, charge, radical, lone pairs, ((neighbour,
    order), ...)) for every atom of the component written in ``form``."""
    radicals = dict(form.radicals)
    orders = {}
    for a, b, o in form.bonds:
        orders[_key(a, b)] = o
    rows = {}
    for a in state.components.atoms(cid):
        n = state.atoms[a]
        rows[a] = (n.element, n.isotope, n.charge, radicals.get(a, n.radical), n.lone_pairs,
                   tuple(sorted((b, float(orders.get(_key(a, b), o))) for b, o in n.edges.items())))
    return rows


# --------------------------------------------------------------------------
# transitions
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Transition:
    """One physically distinct transition of one template class."""

    template: str
    family: str
    class_index: int
    edit_key: Tuple[int, ...]
    components: Tuple[int, ...]
    mapping: Tuple[int, ...]
    edits: Tuple
    form_bonds: Tuple[Tuple[int, int, float], ...]
    form_radicals: Tuple[Tuple[int, int], ...]
    identity: Tuple
    delta: Tuple
    embeddings: int
    weight: float
    state_revision: int
    state_hash: str


@dataclass(frozen=True)
class ClassStats:
    """Per-class counts.  ``embeddings`` counts embeddings over every full
    resonance form of the touched components; ``evaluated`` counts the
    (embedding, touched-subsystem form) outcomes computed before dedup;
    ``sites`` holds per participant {component: sorted edit-role-distinct
    sites}; ``pair_count`` is the participant quotient P."""

    embeddings: int
    evaluated: int
    transitions: Optional[int]
    sites: Tuple
    pair_count: Optional[int]
    role_swap: bool


@dataclass(frozen=True)
class _Embedding:
    cid: int
    mapping: Tuple  # state atom per template atom (None outside the participant)
    choices: Tuple  # tuples of (subsystem, form index) over the touched subsystems
    per_choice: int  # full-form embeddings behind each choice


class MatchResult:
    """Transitions (ordered by edit key) and per-class statistics."""

    def __init__(self, templates, transitions, classes):
        self.templates = templates
        self.transitions = transitions
        self.classes = classes

    def sampler(self, name):
        template = next(t for t in self.templates if t.name == name)
        stats = self.classes[name]
        if stats.pair_count is None:
            raise WindowMatcherError(f"{name} is not a two-participant class")
        return PairSampler(stats.sites[0], stats.sites[1], unordered=template.role_swap,
                           swap=template.swap_positions, symmetry=template.edit_symmetry)


def _piece_embeddings(view, template, piece):
    """Every induced, signature-equal embedding of one connected piece, as
    tuples of (template atom, state atom) in plan order."""
    plan = piece.plan
    root = plan[0][0]
    results = []
    mapping = {}
    used = set()

    def visit(pos):
        if pos == len(plan):
            results.append(tuple((t, mapping[t]) for t, _p, _c in plan))
            return
        t, parent, checks = plan[pos]
        for a in view.adj[mapping[parent]]:
            if a in used or view.sig[a] != template.sig[t]:
                continue
            if any((mapping[u] in view.node[a].edges) != bonded for u, bonded in checks):
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


def _piece_products(lists, view):
    """Disjoint, mutually non-adjacent combinations of piece embeddings."""
    def visit(k, chosen, used):
        if k == len(lists):
            yield tuple(chosen)
            return
        for emb in lists[k]:
            atoms = [a for _t, a in emb]
            if any(a in used for a in atoms):
                continue
            if any(b in used for a in atoms for b in view.adj[a]):
                continue
            chosen.append(emb)
            yield from visit(k + 1, chosen, used | frozenset(atoms))
            chosen.pop()

    yield from visit(0, [], frozenset())


def _participant_embeddings(view, template, participant):
    """Embeddings of one participant's window, with their consistent forms."""
    per_piece = [_piece_embeddings(view, template, piece)
                 for piece in template.pieces[participant]]
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
        for combo in _piece_products([g[cid] for g in by_comp], view):
            mapping = {}
            for emb in combo:
                mapping.update(emb)
            embedding = _finish(view, template, cid, mapping)
            if embedding is not None:
                out.append(embedding)
    return out


def _finish(view, template, cid, mapping):
    """Exterior, path and form checks of a full participant mapping."""
    image = frozenset(mapping.values())
    atoms = template.atoms
    for t in sorted(mapping):
        if atoms[t].kind == "port":
            a = mapping[t]
            outside = tuple(sorted((view.node[b].element, b in view.in_ring)
                                   for b in view.adj[a] if b not in image))
            if outside != template.ext_sigma[t]:
                return None
    for path in template.paths:
        if path.a in mapping:
            d = view.heavy_distance(mapping[path.a], mapping[path.b], path.maximum)
            if d is None or not path.minimum <= d <= path.maximum:
                return None
    # group the pi-level constraints by subsystem (None: fixed by the state)
    groups = {}
    for t in sorted(mapping):
        groups.setdefault(view.sub_of.get(mapping[t]), []).append(("atom", t))
    for (i, j) in template.bond:
        if i in mapping:
            a, b = mapping[i], mapping[j]
            s = view.sub_of.get(a)
            if s is None or _key(a, b) not in view.subs[s].bpos:
                s = None
            groups.setdefault(s, []).append(("bond", i, j))
    static = groups.pop(None, [])
    if not _satisfied(view, template, mapping, image, static, None, None):
        return None
    touched = sorted(groups)
    consistent = []
    for s in touched:
        sub = view.subs[s]
        ok = [k for k, form in enumerate(sub.forms)
              if _satisfied(view, template, mapping, image, groups[s], sub, form)]
        if not ok:
            return None
        consistent.append([(s, k) for k in ok])
    rest = view.comp_forms[cid]
    for s in touched:
        rest //= len(view.subs[s].forms)
    choices = [()]
    for options in consistent:
        choices = [c + (o,) for c in choices for o in options]
    full = tuple(mapping.get(i) for i in range(len(atoms)))
    return _Embedding(cid, full, tuple(choices), rest)


def _satisfied(view, template, mapping, image, constraints, sub, form):
    orders, rads = form if form is not None else (None, None)

    def order(a, b):
        if sub is not None:
            k = sub.bpos.get(_key(a, b))
            if k is not None:
                return orders[k]
        return view.order(a, b)

    for c in constraints:
        if c[0] == "atom":
            t = c[1]
            a = mapping[t]
            atom = template.atoms[t]
            radical = rads[sub.apos[a]] if sub is not None else view.node[a].radical
            if radical != atom.radical:
                return False
            if atom.kind == "port":
                exterior = tuple(sorted((order(a, b), view.node[b].element, b in view.in_ring)
                                        for b in view.adj[a] if b not in image))
                if exterior != atom.exterior:
                    return False
        else:
            _tag, i, j = c
            if order(mapping[i], mapping[j]) != template.bond[(i, j)]:
                return False
    return True


def _concrete_edits(template, mapping):
    edits = []
    for edit in template.edits:
        if isinstance(edit, BondEdit):
            old = template.bond.get(_key(edit.a, edit.b), 0.0)
            edits.append(SetBond(mapping[edit.a], mapping[edit.b], old, float(edit.order)))
        else:
            atom = template.atoms[edit.atom]
            a = mapping[edit.atom]
            if edit.radical is not None:
                edits.append(SetRadical(a, atom.radical, edit.radical))
            if edit.charge is not None:
                edits.append(SetCharge(a, atom.charge, edit.charge))
            if edit.lone_pairs is not None:
                edits.append(SetLonePairs(a, atom.lone_pairs, edit.lone_pairs))
    return tuple(edits)


def _product_delta(view, edits, choice):
    """Product rows (first form) that differ from the reactant rows (first
    form), over the touched subsystems and the edited atoms."""
    region = set()
    for s, _k in choice:
        region.update(view.subs[s].atoms)
    for e in edits:
        if isinstance(e, SetBond):
            region.update((e.a, e.b))
        else:
            region.add(e.atom)
    region = sorted(region)
    vals = {}
    for a in region:
        n = view.node[a]
        vals[a] = [n.radical, n.charge, n.lone_pairs, dict(n.edges)]
    form_bonds, form_radicals = [], []
    for s, k in choice:
        sub = view.subs[s]
        orders, rads = sub.forms[k]
        for a, r in zip(sub.atoms, rads):
            vals[a][0] = r
            form_radicals.append((a, r))
        for (a, b), o in zip(sub.bonds, orders):
            vals[a][3][b] = o
            vals[b][3][a] = o
            form_bonds.append((a, b, o))
    for e in edits:
        if isinstance(e, SetBond):
            if e.new == 0:
                del vals[e.a][3][e.b]
                del vals[e.b][3][e.a]
            else:
                vals[e.a][3][e.b] = e.new
                vals[e.b][3][e.a] = e.new
        elif isinstance(e, SetRadical):
            vals[e.atom][0] = e.new
        elif isinstance(e, SetCharge):
            vals[e.atom][1] = e.new
        else:
            vals[e.atom][2] = e.new

    def value(a):
        v = vals[a]
        return view.node[a].element, v[1], v[2], v[0], v[3]

    for group, free in _subsystems(region, value, outside=view.free.__getitem__):
        bonds, shape = _sub_shape(group, free, value)
        orders, rads = _shape_forms(shape, RESONANCE_FORM_LIMIT)[0]
        for a, r in zip(group, rads):
            vals[a][0] = r
        for (a, b), o in zip(bonds, orders):
            vals[a][3][b] = o
            vals[b][3][a] = o
    delta = []
    for a in region:
        n = view.node[a]
        v = vals[a]
        row = (n.element, n.isotope, v[1], v[0], v[2],
               tuple(sorted((b, float(o)) for b, o in v[3].items())))
        if row != view.base_row(a):
            delta.append((a, row))
    return tuple(delta), tuple(sorted(form_bonds)), tuple(sorted(form_radicals))


def _evaluate(view, template, index, mapping, comps, choice):
    edits = _concrete_edits(template, mapping)
    delta, form_bonds, form_radicals = _product_delta(view, edits, choice)
    identity = (template.name, comps, delta)
    edit_key = (index,) + tuple(mapping[t] for t in template.edited_all)
    return edit_key, identity, edits, delta, form_bonds, form_radicals


def _sites(template, participant, embeddings):
    by_comp = {}
    for emb in embeddings:
        site = tuple(emb.mapping[t] for t in template.edited[participant])
        by_comp.setdefault(emb.cid, set()).add(site)
    return {cid: tuple(sorted(by_comp[cid])) for cid in sorted(by_comp)}


def _pair_count(sites0, sites1, template):
    """P = n0*n1 - sum_c n0(c)*n1(c), halved under role_swap.  When the
    template has participant-preserving symmetries that move edited atoms
    (``edit_symmetry`` larger than the identity), several site pairs are the
    same transition and P is divided by that group's order too."""
    n0 = sum(len(v) for v in sites0.values())
    n1 = sum(len(v) for v in sites1.values())
    pairs = n0 * n1 - sum(len(v) * len(sites1.get(c, ())) for c, v in sites0.items())
    divisor = len(template.edit_symmetry) * (2 if template.role_swap else 1)
    if pairs % divisor:
        raise WindowMatcherError(f"{template.name}: pair count {pairs} is not a multiple of "
                                 f"its symmetry {divisor}")
    return pairs // divisor


def _joint(embeddings):
    for a in embeddings[0]:
        for b in embeddings[1]:
            if a.cid == b.cid:
                continue
            mapping = tuple(x if x is not None else y for x, y in zip(a.mapping, b.mapping))
            comps = tuple(sorted((a.cid, b.cid)))
            choices = tuple(ca + cb for ca in a.choices for cb in b.choices)
            yield mapping, comps, choices, a.per_choice * b.per_choice


def match(templates: Sequence[Template], state, *, joint: bool = True) -> MatchResult:
    """All transitions of ``templates`` on ``state`` (pure and deterministic).

    With ``joint=False`` two-participant classes are only counted (sites and
    the participant quotient), not enumerated pair by pair.
    """
    templates = tuple(templates)
    names = [t.name for t in templates]
    if len(set(names)) != len(names):
        raise TemplateError("template names must be unique")
    view = _view_for(state)
    state_hash = view.full_hash()
    transitions, classes = [], {}
    for index, template in enumerate(templates):
        embeddings = [_participant_embeddings(view, template, p) for p in template.participants]
        sites = tuple(_sites(template, p, embeddings[p]) for p in template.participants)
        two = len(template.participants) == 2
        found = {}
        total = evaluated = 0
        if not two:
            candidates = ((e.mapping, (e.cid,), e.choices, e.per_choice) for e in embeddings[0])
        elif joint:
            candidates = _joint(embeddings)
        else:
            candidates = ()
            total = sum(a.per_choice * len(a.choices) * b.per_choice * len(b.choices)
                        for a in embeddings[0] for b in embeddings[1] if a.cid != b.cid)
        for mapping, comps, choices, count in candidates:
            for choice in choices:
                evaluated += 1
                total += count
                edit_key, identity, edits, delta, fb, fr = _evaluate(
                    view, template, index, mapping, comps, choice)
                previous = found.get(identity)
                if previous is None or edit_key < previous[0]:
                    merged = count + (previous[1] if previous else 0)
                    found[identity] = (edit_key, merged, mapping, comps, edits, delta, fb, fr)
                else:
                    found[identity] = previous[:1] + (previous[1] + count,) + previous[2:]
        mine = []
        for identity, (edit_key, count, mapping, comps, edits, delta, fb, fr) in found.items():
            mine.append(Transition(template.name, template.family, index, edit_key, comps, mapping,
                                   edits, fb, fr, identity, delta, count, template.weight,
                                   view.revision, state_hash))
        mine.sort(key=lambda t: t.edit_key)
        for a, b in zip(mine, mine[1:]):
            if a.edit_key == b.edit_key:
                raise WindowMatcherError(
                    f"{template.name}: two outcomes share edit key {a.edit_key}")
        transitions.extend(mine)
        pair_count = _pair_count(sites[0], sites[1], template) if two else None
        classes[template.name] = ClassStats(total, evaluated,
                                            len(mine) if (joint or not two) else None,
                                            sites, pair_count, template.role_swap)
    transitions.sort(key=lambda t: t.edit_key)
    return MatchResult(templates, tuple(transitions), classes)


def product_rows(state, transition) -> Dict[int, Tuple]:
    """The complete mapped product of the touched components (every
    conjugated subsystem in its first form), atom -> row."""
    view = _view_for(state)
    _check_fresh(view, transition)
    rows = {}
    for cid in transition.components:
        for a in state.components.atoms(cid):
            rows[a] = view.base_row(a)
    rows.update(dict(transition.delta))
    return rows


def _check_fresh(view, transition):
    if view.revision != transition.state_revision or view.full_hash() != transition.state_hash:
        raise StaleTransitionError(
            f"transition matched at revision {transition.state_revision}, state is at "
            f"{view.revision}")


# --------------------------------------------------------------------------
# apply
# --------------------------------------------------------------------------
def _normalizing_edits(state, cid):
    first = component_forms(state, cid)[0]
    edits = []
    for a, b, o in first.bonds:
        current = float(state.atoms[a].edges[b])
        if current != o:
            edits.append(SetBond(a, b, current, o))
    for a, r in first.radicals:
        current = state.atoms[a].radical
        if current != r:
            edits.append(SetRadical(a, current, r))
    return edits


def apply_transition(state, transition, *, fault=None):
    """Apply a transition in one transaction: switch the touched subsystems to
    the matched form, apply the edits, normalize the product components to
    their first form (a fixed point) and check the result equals the matched
    product.  Any failure rolls everything back.  ``fault(stage)`` is a test
    seam called at ``after_form``, ``after_edits`` and ``after_normalize``."""
    view = _view_for(state)
    _check_fresh(view, transition)
    if not any(isinstance(e, SetBond) and (e.old == 0) != (e.new == 0) for e in transition.edits):
        raise WindowMatcherError("the connectivity change is empty (a pure pi move); refused")
    predicted = product_rows(state, transition)
    tx = state.begin()
    try:
        switch = []
        for a, b, o in transition.form_bonds:
            current = float(state.atoms[a].edges[b])
            if current != o:
                switch.append(SetBond(a, b, current, o))
        for a, r in transition.form_radicals:
            current = state.atoms[a].radical
            if current != r:
                switch.append(SetRadical(a, current, r))
        if switch:
            tx.apply(switch)
        if fault:
            fault("after_form")
        tx.apply(transition.edits)
        if fault:
            fault("after_edits")
        cids = sorted({state.components.component_of(a) for a in predicted})
        normalize = [e for cid in cids for e in _normalizing_edits(state, cid)]
        if normalize:
            tx.apply(normalize)
        for cid in cids:
            if _normalizing_edits(state, cid):
                raise WindowMatcherError(f"component {cid} is not at its normalization fixed point")
        if fault:
            fault("after_normalize")
        for a in sorted(predicted):
            n = state.atoms[a]
            actual = (n.element, n.isotope, n.charge, n.radical, n.lone_pairs,
                      tuple(sorted((b, float(o)) for b, o in n.edges.items())))
            if actual != predicted[a]:
                raise WindowMatcherError(f"atom {a}: applied product differs from the matched one")
        return tx.commit()
    except BaseException:
        if tx.open:
            tx.rollback()
        raise


# --------------------------------------------------------------------------
# exact pair sampler
# --------------------------------------------------------------------------
class PairSampler:
    """Exact uniform sampler of participant pairs on different components.

    Component c is drawn with weight n0(c)*(n1 - n1(c)) by a Fenwick search
    over two trees (n0(c) and n0(c)*n1(c)), then a left site uniformly in c and
    a right site uniformly among the n1 - n1(c) sites outside c.  With
    ``unordered`` (role_swap) and ``symmetry`` (``Template.edit_symmetry``)
    the ordered draw is mapped to the canonical representative of its
    transition (the least image of left+right under the symmetry group and,
    if unordered, the exchange ``swap`` = ``Template.swap_positions``, or
    plain exchange when ``swap`` is None); every transition then has
    probability exactly 1/P.
    """

    def __init__(self, sites0, sites1, *, unordered=False, swap=None, symmetry=None):
        self.unordered = bool(unordered)
        self.swap = swap
        self.symmetry = tuple(symmetry) if symmetry else ()
        self._divisor = max(1, len(self.symmetry)) * (2 if self.unordered else 1)
        self._cids = tuple(sorted(set(sites0) | set(sites1)))
        self._pos = {c: k for k, c in enumerate(self._cids)}
        n = len(self._cids)
        self._s0 = [tuple(sites0.get(c, ())) for c in self._cids]
        self._s1 = [tuple(sites1.get(c, ())) for c in self._cids]
        self._a = [0] * (n + 1)
        self._b = [0] * (n + 1)
        self._c = [0] * (n + 1)
        self._n1 = 0
        self._top = 1
        while self._top * 2 <= n:
            self._top *= 2
        for k in range(n):
            n0c, n1c = len(self._s0[k]), len(self._s1[k])
            self._add(k, n0c, n0c * n1c, n1c)
            self._n1 += n1c

    def _add(self, k, da, db, dc):
        i = k + 1
        while i < len(self._a):
            self._a[i] += da
            self._b[i] += db
            self._c[i] += dc
            i += i & -i

    def _prefix(self, tree, k):
        total, i = 0, k
        while i > 0:
            total += tree[i]
            i -= i & -i
        return total

    @property
    def ordered_total(self):
        n = len(self._cids)
        return self._n1 * self._prefix(self._a, n) - self._prefix(self._b, n)

    @property
    def total(self):
        return self.ordered_total // self._divisor

    def _canonical(self, left, right):
        if self._divisor == 1:
            return left, right
        whole = tuple(left) + tuple(right)
        images = [whole]
        if self.unordered:
            images.append(tuple(whole[q] for q in self.swap) if self.swap is not None
                          else tuple(right) + tuple(left))
        if self.symmetry:
            images = [tuple(image[q] for q in g) for image in images for g in self.symmetry]
        best = min(images)
        return best[:len(left)], best[len(left):]

    def update(self, cid, sites0, sites1):
        """Replace the sites of one (existing) component."""
        k = self._pos[cid]
        old0, old1 = len(self._s0[k]), len(self._s1[k])
        self._s0[k], self._s1[k] = tuple(sites0), tuple(sites1)
        new0, new1 = len(self._s0[k]), len(self._s1[k])
        self._add(k, new0 - old0, new0 * new1 - old0 * old1, new1 - old1)
        self._n1 += new1 - old1

    def _left_component(self, rng):
        """Component index drawn with weight n0(c)*(n1 - n1(c))."""
        r = rng.randrange(self.ordered_total)
        pos, step, n = 0, self._top, len(self._cids)
        while step:
            nxt = pos + step
            if nxt <= n:
                w = self._n1 * self._a[nxt] - self._b[nxt]
                if w <= r:
                    r -= w
                    pos = nxt
            step //= 2
        return pos

    def _right(self, k, rng):
        """A participant-1 site uniformly outside component index k."""
        n1c = len(self._s1[k])
        g = rng.randrange(self._n1 - n1c)
        if g >= self._prefix(self._c, k):
            g += n1c
        pos, step, n = 0, self._top, len(self._cids)
        while step:
            nxt = pos + step
            if nxt <= n and self._c[nxt] <= g:
                g -= self._c[nxt]
                pos = nxt
            step //= 2
        return self._s1[pos][g]

    def draw(self, rng):
        """One pair; ``rng`` needs ``randrange``."""
        if self.ordered_total <= 0:
            raise WindowMatcherError("no pair to draw")
        k = self._left_component(rng)
        left = self._s0[k][rng.randrange(len(self._s0[k]))]
        right = self._right(k, rng)
        return self._canonical(left, right)
