"""Recertified closed windows, chain-end contexts and the short-chain ladder
for the declared v1 PS domain (I-137).

For one record member (a channel member: one mapped outcome, with every RMG
mapping of it) on a boundary witness (``rmgpy.kmc.boundaries``), this module
builds a CLOSED WINDOW that provably contains everything RMG read to select
the member's rate, stable on every extension of the witness, or refuses by
name.  It enumerates the chain-end contexts of the member and the finite
short-chain ladder per family.  No templates or weights are cut here.

Support (item 1)
----------------
The support of a member is the union, over every mapping of the member (A+A
template slot orders and reactant resonance variants included, none
skipped), of

* the reacting atoms: labelled atoms (H-only centres seed with their labelled
  H), every atom of ``Move.mapped`` and every atom of the recipe's edits;
* for every predicate RMG tried in every descent (root checks, matched AND
  failed children, logic components in RMG's evaluation order): its anchors,
  the candidates of its free labels (when bounded, below), its ball (the
  predicate's radius around anchors and candidates, ``certificates``), the
  images of its embeddings and the labelled atoms RMG ignored;
* for a reverse record nothing else: the certificate is built on the family
  template side (the record's products), which shares the witness's atom-id
  namespace, so the support is transported by the identity map (``identity``
  is checked: every support atom is a witness atom).

Radius unit: ordinary bonds as in the certificate (``RADIUS_UNIT``); explicit
hydrogens are graph atoms, so a ball reaches the H of its outer atoms only
when the radius does; closure adds every H of an included heavy atom.  The
boundary signature (full bond lists of the ball rows) names the bonds that
leave a ball: in the declared domain (C and H only, no charge or lone pair)
RMG's atom type of a ball atom is a function of its own electrons and bond
orders, never of the far atom's element, so the far atom is not read; the
bond's ORDER is, and the closure makes it a window property (any bond of
order > 1 lies in a pi system or a ring the closure includes whole; checked:
``boundary_bond_problems``).

Not atom sets (each with an explicit disposition, ``Support.dispositions``):

* free-label searches (R_Recombination's two ``*``): RMG scans the whole
  structure.  ``bounded_by_domain`` when no atom an extension can add (an
  interior repeat unit of the provider, its atom types computed by RMG) is a
  specific case of the free group atom (so no extension adds a candidate:
  structural, not observed) AND every candidate is a reacting atom; a
  candidate outside the reacting atoms (a remote radical site) makes the
  selection depend on the whole species: ``whole_species`` -> refused; an
  extension able to add candidates -> refused (``unbounded_search``);
* group multiplicity and radical-count guards (whole molecule): a context
  condition, the species' unpaired-electron count under the declared
  maximum-multiplicity convention; extension invariant because every provider
  fragment is closed shell (checked per fragment); varied only by the radical
  inventory of the contexts, which ``context_assignments`` always enumerates;
  the guard values are part of the certificate view;
* element-count guards: monotone under extension (the net element delta of a
  provider step is non-negative, checked): a satisfied minimum stays
  satisfied; a failed one with zero delta stays failed; a failed one with a
  positive delta could flip at some depth -> refused;
* ``inRing``: closure includes every ring <= ``RING_LIMIT`` through a window
  atom;
* the depository search (whole species): a depository match makes the
  certificate ``whole_species`` -> refused (no window); a rule estimate needs
  a non-refused depository exclusion;
* the multiplicity convention: declared, checked on every species by
  ``domain.check``;
* the channel kind (``kept`` | ``synthetic:unchanged_species``) and the
  identical-reactant pairing: whole-species identity verdicts (an extension of
  one participant changes them), not reads of the selection; both kinds are
  rated by the family's estimate for the template, so the kind is reported
  as a disposition and left out of the certificate view.

Closure (item 2)
----------------
On BOTH reaction sides (the window is the union): every ring <= ``RING_LIMIT``
through a window atom, whole (a reacting atom in a larger ring -> refused; a
larger ring touching the window without a reacting atom is refused too: the
declared domain admits six-membered phenyl rings only, so such a ring is a
domain error, never silently truncated); the whole conjugated pi / radical
system (atoms with an unpaired electron or a bond of order > 1, joined by
bonds) touching the window; every H of an included heavy atom; repeated to a
fixed point.

Recertified fixed point (item 3)
--------------------------------
``fixed_point``: observe (certify, support, close); while an artificial port
(port atom or a declared cap) lies in the window, extend every such port by
one repeat (``boundaries.extend``), recertify and re-close.  Then the
independent check: for EVERY artificial port alone, +1 repeat, recertify and
close; the fixed point needs every check to leave the window's marked graph
(identity correspondence: surviving atoms keep their ids) and the
certificate view (``certificate_view``) unchanged.  Failed checks extend
their ports by one and the round repeats.  Pins: ports in vector order
(participant, port id); ``PORT_DEPTH_LIMIT`` repeats per port; at most
``FIXED_POINT_RETRIES`` re-extensions after failed checks; the heavy-atom
ceiling ``PROXY_PADDING_MAX_HEAVY_ATOMS``.  Statuses: ``fixed_point``,
``no_fixed_point`` (a port is still in the window at its limit, or windows
keep changing past the retries), ``ceiling_exceeded`` (never truncated),
``certificate_unstable`` (the certificate keeps changing past the retries),
and ``refused:<reason>`` from the observation.

Span bound (item 4)
-------------------
Per family, from the declared domain and the family's tree constants
(``family_constants``): ``r`` the largest predicate radius of the tree,
``D`` the largest bond distance between labelled heavy atoms of one
connected template structure (reactants, and products after the recipe).
Every support atom lies within ``r`` bonds of a reacting atom, so the
support's backbone extent is at most ``D + 2r`` bonds.  Closure extends the
backbone by at most ``CLOSURE_EXTENSION = 1`` bond per side: rings are
pendant phenyls (any ring atom's backbone carbon is nearer the centre than
the ring atom), radical pi systems are a backbone carbon plus its pendant ring
(allylic radicals are outside the domain), and the only backbone pi system
is an end alkene (two carbons).  Hence a window spans at most
``A = D + 2(r + 1) + 1`` backbone atoms of one chain, i.e. at most
``ceil((A + 1) / 2)`` repeat units (two backbone carbons per repeat, any
phase).  A window exceeding its bound fails the run (``SpanBoundError``).

Ladder (item 6)
---------------
A chain of DP ``n`` has ``2n - h + t`` backbone carbons (``h``/``t``: the
head/tail starts/ends on the other repeat phase, ``PsChainSpec.head_shift``
/``tail_shift``), at least ``2n - 1``.  A window can hold both physical ends
of a chain iff the chain's backbone fits its span: ``L = floor((A + 1) / 2)``
(the largest n with 2n - 1 <= A; in-phase chains need 2n <= A, so the
repeat phase and the end groups decide which chains of DP L reach).  DP
``L + 1`` (``2L + 1 > A`` backbone carbons in every phase) is the
non-overlap control.

Contexts (item 5)
-----------------
Per artificial port: artificial, or physical with a declared end chemistry
built by the provider (``boundaries.reconstruct_end``: caps deleted, fresh
end atoms, never a relabelled cap); the event is rematched on the context
by ``rmg_moves`` (same labels, same mapped outcome), recertified and closed.
Two contexts merge only if their closed windows are equal as marked graphs
(both sides, participant roles, reacting atoms, physical ends, oriented
exterior bonds) and their certificate views are equal under that
correspondence.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import combinations, product
from typing import Any, Callable, Optional, Sequence

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import domain as dom

RING_LIMIT = 8
PORT_DEPTH_LIMIT = 6
FIXED_POINT_RETRIES = 3
CLOSURE_EXTENSION = 1
ISOMORPHISM_LIMIT = 256
LOOP_STATUSES = ("fixed_point", "no_fixed_point", "ceiling_exceeded", "certificate_unstable")
DISPOSITIONS = ("bounded_by_domain", "context_condition", "monotone_satisfied",
                "invariant_failed", "closure", "declared", "whole_species",
                "whole_species_identity", "refused")


class WindowError(ValueError):
    """Bad input."""


class SpanBoundError(AssertionError):
    """A window exceeds its family's derived span bound (a bug or a domain
    error): fails the run."""


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Support:
    """``atoms``: the support; ``reacting``: the centres (labels, mapped atoms,
    edits); ``sources``: (node, result, atoms) per predicate; ``dispositions``:
    (item, disposition, argument) for every non-atom-set read; ``refusal``:
    why no window can be built (None: it can)."""

    atoms: frozenset
    reacting: frozenset = frozenset()
    sources: tuple = ()
    dispositions: tuple = ()
    refusal: Optional[str] = None


@dataclass(frozen=True)
class Window:
    atoms: frozenset
    support: Support
    added: tuple = ()
    refusal: Optional[str] = None


@dataclass(frozen=True)
class Observation:
    """One observation of a witness: support, the certificate view (``viewer``:
    ``(inside, transport=None) -> hashable``), the certificate itself and
    the recertification against the base witness."""

    support: Support
    viewer: Callable = field(compare=False)
    certificate: Any = field(default=None, compare=False)
    refusal: Optional[str] = None
    recertification: Any = field(default=None, compare=False)

    def view(self, inside, transport=None):
        if transport is None:
            return self.viewer(inside)
        return self.viewer(inside, transport)


@dataclass(frozen=True)
class FamilyConstants:
    family: str
    radius: int
    reacting_distance: int
    free_label_search: bool
    bimolecular: bool = False
    directions: tuple = ()


@dataclass(frozen=True)
class FixedPointResult:
    status: str
    reason: str
    depths: tuple = ()
    witness: Any = field(default=None, compare=False)
    window: Optional[Window] = field(default=None, compare=False)
    observation: Optional[Observation] = field(default=None, compare=False)
    checks: tuple = field(default=(), compare=False)
    retries: int = 0
    history: tuple = field(default=(), compare=False)
    flips: tuple = field(default=(), compare=False)

    @property
    def reached(self) -> bool:
        return self.status == "fixed_point"


@dataclass(frozen=True)
class ContextResult:
    assignment: tuple
    status: str
    reason: str = ""
    witness: Any = field(default=None, compare=False)
    window: Any = field(default=None, compare=False)
    observation: Any = field(default=None, compare=False)
    fixed_point: Any = field(default=None, compare=False)


# ---------------------------------------------------------------------------
# graph helpers on witness rows
# ---------------------------------------------------------------------------
def _side_rows(witness, side):
    return witness.reactant_rows() if side == "reactants" else witness.product_rows()


def _neighbours(rows, atom):
    return [n for n, _o in rows[atom][5]]


def _pi_atoms(rows) -> set:
    return {a for a, row in rows.items()
            if row[3] > 0 or any(float(o) > 1.0 for _n, o in row[5])}


def _pi_component(rows, pi, start) -> set:
    seen, stack = {start}, [start]
    while stack:
        a = stack.pop()
        for n in _neighbours(rows, a):
            if n in pi and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def _cyclic_atoms(rows) -> set:
    """Atoms on a cycle (incident to a non-bridge bond)."""
    import networkx as nx
    graph = nx.Graph()
    graph.add_nodes_from(rows)
    for a, row in rows.items():
        for n, _o in row[5]:
            if a < n:
                graph.add_edge(a, n)
    bridges = {tuple(sorted(e)) for e in nx.bridges(graph)}
    return {a for a in rows for n in _neighbours(rows, a) if tuple(sorted((a, n))) not in bridges}


def smallest_ring(rows, atom, cyclic=None, limit: int = 64) -> Optional[tuple]:
    """Atoms of a smallest ring through ``atom`` (None if it is on none, a
    tuple longer than ``limit`` atoms is never searched: ``()`` then)."""
    cyclic = _cyclic_atoms(rows) if cyclic is None else cyclic
    if atom not in cyclic:
        return None
    best = None
    for first in _neighbours(rows, atom):
        if first not in cyclic:
            continue
        # shortest path first -> atom avoiding the bond atom-first
        parent = {first: None}
        frontier = [first]
        depth = 0
        found = False
        while frontier and depth < limit and not found:
            depth += 1
            nxt = []
            for a in frontier:
                for n in _neighbours(rows, a):
                    if a == first and n == atom:
                        continue
                    if n == atom:
                        path = [a]
                        while parent[path[-1]] is not None:
                            path.append(parent[path[-1]])
                        ring = tuple(sorted(set(path) | {atom}))
                        if best is None or len(ring) < len(best):
                            best = ring
                        found = True
                        break
                    if n not in parent:
                        parent[n] = a
                        nxt.append(n)
                if found:
                    break
            frontier = nxt
    return best if best is not None else ()


# ---------------------------------------------------------------------------
# closure (item 2)
# ---------------------------------------------------------------------------
def close(witness, support, *, ring_limit=RING_LIMIT, sides=("reactants", "products")):
    """The closed window of ``support`` (a ``Support`` or an atom set) on
    ``witness`` (see the module docstring); a refusal is returned as a
    ``Window`` with ``refusal`` set and no atoms beyond the support."""
    if not isinstance(support, Support):
        support = Support(atoms=frozenset(support))
    rows_of = {side: _side_rows(witness, side) for side in sides}
    every = set(witness.reactant_rows())
    unknown = sorted(set(support.atoms) - every)
    if unknown:
        raise WindowError(f"support atoms {unknown[:5]} are not witness atoms")
    window = set(support.atoms)
    centres = set(support.reacting) or set(support.atoms)
    info = {side: (rows, _cyclic_atoms(rows), _pi_atoms(rows)) for side, rows in rows_of.items()}
    ring_cache = {}
    added = []
    changed = True
    while changed:
        changed = False
        for side, (rows, cyclic, pi) in info.items():
            for atom in sorted(window):
                if rows[atom][0] == "H":
                    continue
                key = (side, atom)
                if key not in ring_cache:
                    ring_cache[key] = smallest_ring(rows, atom, cyclic, limit=4 * ring_limit)
                ring = ring_cache[key]
                if ring is not None:
                    if ring == () or len(ring) > ring_limit:
                        size = len(ring) if ring else f"> {4 * ring_limit}"
                        if atom in centres:
                            return Window(frozenset(support.atoms), support, tuple(added),
                                          f"refused:ring: reacting atom {atom} is in a "
                                          f"{size}-membered ring (> {ring_limit}) on the "
                                          f"{side} side")
                        return Window(frozenset(support.atoms), support, tuple(added),
                                      f"refused:large_ring: window atom {atom} is in a "
                                      f"{size}-membered ring (> {ring_limit}) on the {side} "
                                      "side (the declared domain admits phenyl rings only)")
                    new = set(ring) - window
                    if new:
                        window |= new
                        added.append(("ring", side, tuple(sorted(new))))
                        changed = True
                if atom in pi:
                    new = _pi_component(rows, pi, atom) - window
                    if new:
                        window |= new
                        added.append(("pi", side, tuple(sorted(new))))
                        changed = True
            hydrogens = {n for a in window if rows[a][0] != "H" for n in _neighbours(rows, a)
                         if rows[n][0] == "H"} - window
            if hydrogens:
                window |= hydrogens
                added.append(("H", side, tuple(sorted(hydrogens))))
                changed = True
    return Window(frozenset(window), support, tuple(added), None)


def boundary_bond_problems(witness, window) -> list:
    """Bonds leaving the window whose order is not 1 on some side (the order
    of a boundary bond must be a window property: pi systems and rings are
    closed)."""
    out = []
    atoms = window.atoms if isinstance(window, Window) else set(window)
    for side in ("reactants", "products"):
        rows = _side_rows(witness, side)
        for a in sorted(atoms):
            for n, o in rows[a][5]:
                if n not in atoms and float(o) != 1.0:
                    out.append(f"{side}: bond {a}-{n} of order {o} leaves the window")
    return out


def port_inside(witness, port_key, atoms) -> bool:
    port = witness.port(port_key)
    return port.atom_id in atoms or any(c in atoms for c in port.cap_atom_ids)


# ---------------------------------------------------------------------------
# window marked graphs
# ---------------------------------------------------------------------------
def _port_atoms(witness) -> dict:
    out = {}
    for participant in witness.participants:
        for port in participant.metadata.ports:
            out[port.atom_id] = (participant.index, port.port_id)
            for cap in port.cap_atom_ids:
                out[cap] = (participant.index, port.port_id)
    return out


def _roles(witness) -> dict:
    return {a: p.index for p in witness.participants for a in p.atom_ids}


def window_labels(witness, atoms) -> tuple:
    """(node labels {atom: label}, edge labels {(a, b): label}) of the window's
    marked graph: attributes on both sides, participant role, reacting mark
    (the edit support), physical-end mark, and per exterior bond its orders
    and the ports its far side reaches (orientation)."""
    atoms = set(atoms)
    before, after = witness.reactant_rows(), witness.product_rows()
    support = b.edit_support(witness)
    roles = _roles(witness)
    ports = _port_atoms(witness)
    physical = {}
    for participant in witness.participants:
        for port in participant.metadata.ports:
            if not port.artificial:
                physical[port.atom_id] = (participant.index, port.port_id, port.tag)
    bonds_before, bonds_after = witness.reactant_bonds(), witness.product_bonds()
    adjacency = defaultdict(set)
    for pair in set(bonds_before) | set(bonds_after):
        adjacency[pair[0]].add(pair[1])
        adjacency[pair[1]].add(pair[0])
    reach_cache = {}

    def reach(start):
        # port keys reachable from ``start`` outside the window (either side)
        if start in reach_cache:
            return reach_cache[start]
        seen, stack, found = {start}, [start], set()
        while stack:
            a = stack.pop()
            if a in ports:
                found.add(ports[a])
            for n in adjacency[a]:
                if n not in atoms and n not in seen:
                    seen.add(n)
                    stack.append(n)
        value = tuple(sorted(found))
        for a in seen:
            reach_cache[a] = value
        return value

    nodes = {}
    for a in sorted(atoms):
        exterior = []
        for n in sorted(adjacency[a]):
            if n in atoms:
                continue
            pair = (min(a, n), max(a, n))
            exterior.append((bonds_before.get(pair, 0.0), bonds_after.get(pair, 0.0), reach(n)))
        nodes[a] = (b._attrs(before[a]), b._attrs(after[a]), roles.get(a, -1),
                    a in support, physical.get(a), tuple(sorted(exterior, key=repr)))
    edges = {}
    for a in atoms:
        for n in adjacency[a]:
            if n in atoms and a < n:
                edges[(a, n)] = (bonds_before.get((a, n), 0.0), bonds_after.get((a, n), 0.0))
    return nodes, edges


def same_window(first, atoms_first, second, atoms_second) -> bool:
    """Identity correspondence (surviving atoms keep their ids)."""
    if set(atoms_first) != set(atoms_second):
        return False
    return window_labels(first, atoms_first) == window_labels(second, atoms_second)


def _collapsed_graph(witness, atoms):
    """networkx graph of the window: heavy atoms and reacting H as nodes,
    the other H folded into their heavy neighbour's label; reacting atoms
    carry their id (marks fix them)."""
    import networkx as nx
    nodes, edges = window_labels(witness, atoms)
    rows = witness.reactant_rows()
    support = b.edit_support(witness)
    kept = {a for a in nodes if rows[a][0] != "H" or a in support}
    folded = defaultdict(list)
    for a in nodes:
        if a in kept:
            continue
        heavy = [n for n, _o in rows[a][5] if n in kept]
        folded[heavy[0] if heavy else None].append(repr(nodes[a]))
    graph = nx.Graph()
    for a in kept:
        mark = f"r{a}" if a in support else ""
        graph.add_node(a, label=repr((nodes[a], mark, tuple(sorted(folded.get(a, ()))))))
    graph.graph["orphans"] = tuple(sorted(folded.get(None, ())))
    for (x, y), label in edges.items():
        if x in kept and y in kept:
            graph.add_edge(x, y, label=repr(label))
    return graph


def window_isomorphisms(first, atoms_first, second, atoms_second, limit=ISOMORPHISM_LIMIT):
    """Mark-preserving isomorphisms (second -> first ids over heavy atoms and
    reacting H) of two windows' marked graphs (at most ``limit``)."""
    from networkx.algorithms import isomorphism
    g1, g2 = _collapsed_graph(first, atoms_first), _collapsed_graph(second, atoms_second)
    if g1.graph["orphans"] != g2.graph["orphans"] or len(g1) != len(g2):
        return []
    matcher = isomorphism.GraphMatcher(
        g2, g1, node_match=lambda x, y: x["label"] == y["label"],
        edge_match=lambda x, y: x["label"] == y["label"])
    out = []
    for mapping in matcher.isomorphisms_iter():
        out.append(dict(mapping))
        if len(out) >= limit:
            break
    return out


# ---------------------------------------------------------------------------
# certificate support and view (item 1)
# ---------------------------------------------------------------------------
def _label_set(mapping) -> frozenset:
    return frozenset((lab, a) for _s, lab, a in mapping["labels"])


def member_entries(cert, wanted) -> list:
    """(channel, member) of every certified member one of whose mappings has
    a label set in ``wanted``."""
    out = []
    for ch in cert.get("channels", []):
        for mem in ch.get("members", []):
            if any(_label_set(mp) in wanted for mp in mem["mappings"]):
                out.append((ch, mem))
    return out


def mapping_atoms(mapping) -> set:
    atoms = {a for _s, _l, a in mapping["labels"]}
    atoms |= {a for _s, _i, _j, a in mapping.get("mapped", [])}
    edits = mapping.get("edits") or {}
    atoms |= {a for a, _b, _c in edits.get("atoms", [])}
    for a, c, _o1, _o2 in edits.get("bonds", []):
        atoms |= {a, c}
    return atoms


class SearchClassifier:
    """Whether an extension can add a candidate to a free-label search:
    tests the RMG group atoms of the label against every atom of an interior
    repeat unit of each provider (RMG atom types)."""

    def __init__(self, family, providers: Sequence[str] = ("polystyrene", "linear_alkane")):
        self.family = family
        self.fresh = []
        for name in providers:
            self.fresh.extend(_interior_unit_atoms(name))
        self._cache = {}

    def extensible(self, node, label) -> bool:
        key = (node, label)
        if key not in self._cache:
            from rmgpy.kmc import channels as chm
            entry = self.family.groups.entries.get(node)
            if entry is None:
                self._cache[key] = True    # unknown semantics: never claimed bounded
            else:
                hit = False
                for group in chm._group_atoms(self.family.groups, entry):
                    centres = group.get_all_labeled_atoms().get(label)
                    if centres is None:
                        continue
                    gatoms = centres if isinstance(centres, list) else [centres]
                    if any(atom.is_specific_case_of(g) for atom in self.fresh for g in gatoms):
                        hit = True
                        break
                self._cache[key] = hit
        return self._cache[key]


_INTERIOR = {}


def _interior_unit_atoms(provider: str) -> list:
    """RMG atoms (atom types perceived) of an interior repeat unit: unit 1 of a
    pristine PS DP3, carbon 2 of a C5 alkane, and their H."""
    if provider not in _INTERIOR:
        if provider == "polystyrene":
            build = b.PROVIDERS["polystyrene"].build("interior", b.PsChainSpec(3))
            coords = {a: (u, k) for a, u, k in build.metadata.coordinates}
            heavy = {a for a, (u, _k) in coords.items() if u == 1}
        elif provider == "linear_alkane":
            build = b.PROVIDERS["linear_alkane"].build(5)
            heavy = {2}
        else:
            _INTERIOR[provider] = []
            return []
        mol = build.molecule
        atoms = [a for a in mol.atoms if a.id in heavy]
        atoms += [h for a in atoms for h in a.bonds if h.is_hydrogen()]
        _INTERIOR[provider] = atoms
    return _INTERIOR[provider]


def _fragment_facts(providers=None) -> dict:
    """Per provider: whether every extension fragment is closed shell and the
    net element delta of one step (fresh atoms minus the deleted H cap)."""
    out = {}
    providers = b.PROVIDERS if providers is None else providers
    for name, provider in providers.items():
        facts = {}
        for tag, coordinate in (("backbone_tail", (0, 1)), ("backbone_head", (0, 0))):
            port = b.Port("p", 0, (1,), tag, "artificial")
            try:
                fragment = provider.fragment(port, coordinate)
            except b.BoundaryError:
                continue
            delta = defaultdict(int)
            for atom in fragment.atoms:
                delta[atom[0]] += 1
            delta["H"] -= 1
            facts[tag] = {"closed_shell": all(atom[3] == 0 for atom in fragment.atoms),
                          "delta": dict(delta)}
        out[name] = facts
    return out


def certificate_support(cert, wanted, *, classifier: Optional[SearchClassifier] = None,
                        fragment_facts=None, include_failed: bool = True) -> Support:
    """The support of the member of ``cert`` selected by ``wanted`` (label
    sets), with the disposition of every non-atom-set read (see the module
    docstring)."""
    entries = member_entries(cert, wanted)
    if not entries:
        return Support(frozenset(), refusal="refused:member: no certified member carries the "
                                            "wanted labels")
    if len({id(ch) for ch, _m in entries}) > 1:
        return Support(frozenset(), refusal=f"refused:member: the member is in "
                                            f"{len(entries)} channels")
    facts = fragment_facts if fragment_facts is not None else _fragment_facts()
    closed_shell = all(f["closed_shell"] for p in facts.values() for f in p.values())
    deltas = defaultdict(set)
    for p in facts.values():
        for f in p.values():
            for el, d in f["delta"].items():
                deltas[el].add(d)
    reacting, atoms, sources = set(), set(), []
    dispositions = {}
    refusal = None
    for ch, mem in entries:
        for mp in mem["mappings"]:
            reacting |= mapping_atoms(mp)
    for ch, mem in entries:
        for mp in mem["mappings"]:
            if mp.get("descent") is None:
                refusal = refusal or f"refused:incomplete: mapping without descent: {mp.get('gap')}"
                continue
            for d in cert["descents"][mp["descent"]]:
                if d.get("skipped"):
                    continue
                for p in _iter_predicates(d):
                    if p.get("logic"):
                        continue
                    if p["result"] == "failed" and not include_failed:
                        continue
                    pa = set((p.get("anchors") or {}).values())
                    pa |= {i for i, _d in p.get("ball", [])}
                    emb = p.get("embeddings") or {}
                    pa |= set(emb.get("images", []))
                    pa |= set(p.get("ignored", []))
                    for lab in p.get("free_labels", []):
                        cands = set((p.get("candidates") or {}).get(lab, []))
                        key = ("search", p["node"], lab)
                        if classifier is None or classifier.extensible(p["node"], lab):
                            dispositions[key] = ("refused", "an extension can add candidates "
                                                 "of this free label (unbounded search)")
                            refusal = refusal or (f"refused:unbounded_search: {p['node']} "
                                                  f"label {lab}")
                        elif not cands <= reacting:
                            dispositions[key] = ("whole_species", "a candidate outside the "
                                                 f"reacting atoms: {sorted(cands - reacting)}")
                            refusal = refusal or (f"refused:whole_species_search: {p['node']} "
                                                  f"label {lab} candidates "
                                                  f"{sorted(cands - reacting)[:4]}")
                        else:
                            dispositions[key] = ("bounded_by_domain", "no interior repeat-unit "
                                                 "atom matches the free group atom (no extension "
                                                 "adds a candidate); every candidate is a "
                                                 "reacting atom")
                    for name, g in sorted((p.get("guards") or {}).items()):
                        if not g.get("read"):
                            continue
                        if name in ("multiplicity", "radical_count"):
                            if not closed_shell:
                                dispositions[("guard", name)] = ("refused", "a provider "
                                                                 "fragment is not closed shell")
                                refusal = refusal or f"refused:guard: {name}"
                            else:
                                dispositions[("guard", name)] = (
                                    "context_condition", "whole-species unpaired-electron count "
                                    "(declared maximum multiplicity); extension invariant "
                                    "(closed-shell fragments); varied by the contexts' radical "
                                    "inventory; its value is in the certificate view")
                        elif name == "element_count":
                            minimum = g.get("group_minimum") or {}
                            have = g.get("structure") or {}
                            for el, need in sorted(minimum.items()):
                                count = have.get(el, 0)
                                ds = deltas.get(el, {0})
                                key = ("guard", f"element_count:{el}")
                                if any(x < 0 for x in ds):
                                    dispositions[key] = ("refused", "an extension removes "
                                                         f"{el} atoms")
                                    refusal = refusal or f"refused:guard: element_count {el}"
                                elif count >= need:
                                    dispositions[key] = ("monotone_satisfied", "satisfied; "
                                                         "extension never lowers the count")
                                elif ds == {0}:
                                    dispositions.setdefault(key, (
                                        "invariant_failed", "failed; extension never adds "
                                        f"{el}"))
                                else:
                                    dispositions[key] = ("refused", f"failed ({count} < {need}); "
                                                         "an extension could satisfy it")
                                    refusal = refusal or f"refused:guard: element_count {el}"
                    sources.append((p["node"], p["result"], frozenset(pa)))
                    atoms |= pa
    atoms |= reacting
    for ch, _mem in entries:
        dispositions[("global", "channel_kind")] = (
            "whole_species_identity", f"{ch.get('kind')}: kept vs synthetic:unchanged_species "
            "(products equal to the reactant species) and the identical-reactant pairing are "
            "decided by whole-species identity; both kinds are rated by the family's own "
            "estimate for the template (channels.channel_kinetics), which the view carries")
    for guard in cert.get("global_guards", []):
        name = guard["guard"]
        if name == "multiplicity_convention":
            dispositions[("global", name)] = ("declared", "domain.MULTIPLICITY_CONVENTION, "
                                              "checked on every species by domain.check")
        elif name == "inRing":
            dispositions[("global", name)] = ("closure", "every ring <= RING_LIMIT through a "
                                              "window atom is in the window")
        elif name == "depository_whole_species":
            for ch, _mem in entries:
                if ch.get("whole_species"):
                    dispositions[("global", name)] = ("whole_species", "a depository entry "
                                                      "matched the whole species")
                    refusal = refusal or "refused:whole_species: depository match"
                else:
                    verdict = (ch.get("depository_exclusion") or {}).get("verdict", "")
                    if verdict.startswith("refused") or not verdict:
                        dispositions[("global", name)] = ("refused", f"depository exclusion "
                                                          f"{verdict or 'missing'}")
                        refusal = refusal or f"refused:depository_exclusion: {verdict}"
                    else:
                        dispositions[("global", name)] = ("bounded_by_domain",
                                                          f"depository exclusion: {verdict}")
    return Support(frozenset(atoms), frozenset(reacting), tuple(sources),
                   tuple(sorted((k, v[0], v[1]) for k, v in dispositions.items())), refusal)


def _iter_predicates(descent):
    from rmgpy.kmc import certificates as cm
    return cm.iter_descent_predicates(descent)


def certificate_refusal(cert) -> Optional[str]:
    """Why a certificate cannot carry a window (None: it can)."""
    if cert is None:
        return "refused:certificate: none"
    if not cert.get("complete"):
        return f"refused:incomplete: {(cert.get('gaps') or ['?'])[0]}"[:300]
    if not cert.get("eligible"):
        why = cert.get("input_refusals") or [
            (ch.get("eligibility") or {}).get("disposition") for ch in cert["channels"]
            if (ch.get("eligibility") or {}).get("disposition") != "eligible"]
        return f"refused:domain: {why[:2]}"[:300]
    return None


def member_refusal(cert, wanted) -> Optional[str]:
    for ch, _mem in member_entries(cert, wanted):
        if ch.get("support") != "supported":
            return f"refused:{ch.get('support')}"
        if not ch.get("publishable", True):
            return f"refused:unpublishable: {(ch.get('eligibility') or {}).get('disposition')}"
    return None


def certificate_view(cert, wanted, inside, transport=None):
    """The member's certificate, correspondence-aware, restricted to the
    window ``inside``: per mapping (identity with ids through the map) its
    channel selection and status, descent signature (every predicate tried,
    verdicts, choices, leaves) and every predicate's local support (anchors,
    candidates, ball, embeddings, ignored atoms) and guard verdicts, and the
    atom-table rows of window atoms.  Non-reacting H are named by their heavy
    neighbour; atoms outside the window are ``out``."""
    from rmgpy.kmc import certificates as cm
    inside = set(inside)
    entries = member_entries(cert, wanted)
    reacting = set()
    for _ch, mem in entries:
        for mp in mem["mappings"]:
            reacting |= mapping_atoms(mp)
    witness_atoms = {row[0]: row for row in cert["witness"]["atoms"]}
    base = (lambda a: a) if transport is None else (lambda a: transport.get(a, ("unmapped", a)))

    def t(a):
        if a in reacting:
            return ("r", a)
        if a not in inside:
            return "out"
        row = witness_atoms.get(a)
        if row is not None and row[1] == "H" and row[3]:
            return ("H", t(row[3][0]))
        return base(a)

    def tl(values):
        return tuple(sorted((t(v) for v in values), key=repr))

    out = []
    for ch, mem in entries:
        # the channel kind (kept | synthetic:unchanged_species) is a
        # whole-species identity verdict, not a read of the selection: it is
        # a disposition (``certificate_support``), not part of the view
        selection = list(cm._channel_selection(ch))
        selection[1] = None
        sel = tuple(selection) + (
            ch.get("support"), (ch.get("eligibility") or {}).get("disposition"),
            ch.get("whole_species"), (ch.get("depository_exclusion") or {}).get("verdict"))
        for mp in mem["mappings"]:
            descents = cert["descents"].get(mp.get("descent"), []) if mp.get("descent") else []
            preds = []
            rows = []
            for d in descents:
                if d.get("skipped"):
                    continue
                for p in _iter_predicates(d):
                    if p.get("logic"):
                        preds.append((p["node"], p["result"], "logic"))
                        continue
                    emb = p.get("embeddings") or {}
                    guards = []
                    for name, g in sorted((p.get("guards") or {}).items()):
                        if name == "element_count":
                            have = g.get("structure") or {}
                            guards.append((name, g.get("read"), tuple(sorted(
                                (el, have.get(el, 0) >= need)
                                for el, need in (g.get("group_minimum") or {}).items()))))
                        else:
                            guards.append((name, g.get("read"), g.get("structure")))
                    preds.append((p["node"], p["result"], p.get("anchoring"),
                                  tuple(sorted((lab, t(a)) for lab, a in
                                               (p.get("anchors") or {}).items())),
                                  tuple(sorted((lab, tl(ids)) for lab, ids in
                                               (p.get("candidates") or {}).items())),
                                  tuple(sorted(((t(i), dd) for i, dd in p.get("ball", [])),
                                               key=repr)),
                                  (emb.get("count"),
                                   tuple(sorted((int(g), tl(ids))
                                                for g, ids in emb.get("by_group", []))),
                                   tl(emb.get("images", []))),
                                  tuple(guards), tl(p.get("ignored", []))))
                for row in d.get("atoms", []):
                    if row[0] in inside:
                        rows.append((t(row[0]),) + tuple(row[1:8]) + (
                            tuple(sorted(((t(n), o) for n, o in row[8]), key=repr)),))
            key = cm.mapping_key(mp, t)
            out.append((key, sel, cm.descent_signature(descents), tuple(preds),
                        tuple(sorted(rows, key=repr))))
    return tuple(sorted(out, key=repr))


# ---------------------------------------------------------------------------
# observers
# ---------------------------------------------------------------------------
def connectivity_change(family, labels) -> frozenset:
    """The bonds the FORWARD (family-template) recipe makes or breaks, as atom
    pairs, for a template-side label set ``{(label, atom id), ...}`` (a label
    may name two atoms: R_Recombination's two ``*``).  Two mappings with the
    same mapped outcome make and break the same bonds (the outcome fixes the
    product connectivity; resonance moves only bond orders and electrons)."""
    atoms = defaultdict(list)
    for label, atom in labels:
        atoms[label].append(atom)
    out = set()
    for action in family.forward_recipe.actions:
        if action[0] not in ("FORM_BOND", "BREAK_BOND"):
            continue
        first, second = action[1], action[3]
        if first == second:
            pairs = [tuple(sorted(atoms[first]))] if len(atoms[first]) == 2 else []
        else:
            pairs = [tuple(sorted((x, y))) for x in atoms[first] for y in atoms[second]]
        out.update((action[0], pair) for pair in pairs)
    return frozenset(out)


def _candidate_labels(candidate) -> frozenset:
    return frozenset((lab, a) for _s, lab, a in candidate.labels)


class ScopedCertifier:
    """A ``certificates.Certifier`` whose forward channel sets keep only the
    mappings that make and break the member's bonds (``connectivity_change``):
    every mapping of the member (all slot orders, all resonance variants) and
    every mapping of any other pathway with the same mapped outcome (pathway
    ambiguity), none of the unrelated ones.  The member's mappings, descents,
    selection and channel status are those of the full certificate (checked:
    ``windowsTest``); the channel's other members (other bonds) are other
    windows.  Shares every other cache with the full certifier."""

    def __new__(cls, base, family, change):
        from rmgpy.kmc import certificates as cm

        class _Scoped(cm.Certifier):
            def channel_set(self, family_, reactants):
                from rmgpy.kmc import channels as chm
                from rmgpy.kmc import rmg_moves as rm
                if family_.label != self._scope_family:
                    return cm.Certifier.channel_set(self, family_, reactants)
                roots = [cm.aromatic_root(m) for m in reactants]
                key = ("scoped", family_.label, tuple(sorted(self._scope_change)),
                       tuple(rm._signature(rm.molecule_rows(m)) for m in roots))
                cs = self._channel_sets.get(key)
                if cs is None:
                    cs = chm.build_channels(
                        family_, roots, "forward", resonance=self.resonance,
                        select=lambda c: connectivity_change(
                            family_, _candidate_labels(c)) == self._scope_change)
                    self._channel_sets[key] = cs
                return cs

        scoped = _Scoped.__new__(_Scoped)
        scoped.__dict__.update(base.__dict__)
        scoped._scope_family = family.label
        scoped._scope_change = frozenset(change)
        return scoped


def scoped_certifier(certifier, family, wanted):
    """A ``ScopedCertifier`` for the member whose template-side label sets are
    ``wanted`` (all must share one connectivity change)."""
    changes = {connectivity_change(family, labels) for labels in wanted}
    if len(changes) != 1:
        raise WindowError(f"the member's mappings make different bonds: {changes}")
    return ScopedCertifier(certifier, family, changes.pop())


def template_molecules(witness, direction: str):
    """RMG molecules of the family-template side of ``witness``: its reactants
    for a forward record, its products for a reverse record (components by
    least atom id; atom ids kept)."""
    from rmgpy.kmc import rmg_moves as rm
    rows = witness.reactant_rows() if direction == "forward" else witness.product_rows()
    return [rm.molecule_from_rows(part) for part in rm.split_rows(dict(rows))]


class MemberObserver:
    """Observe a member on a witness: certify the family-template side (the
    first observation is the base; later ones are ``recertify_certificate``
    against it through the identity map of surviving atoms), then support,
    dispositions and the certificate view."""

    def __init__(self, certifier, family, direction, wanted, *, base_certificate=None,
                 classifier=None, include_failed=True, scoped=True):
        self.wanted = frozenset(wanted)
        self.certifier = scoped_certifier(certifier, family, self.wanted) if scoped else \
            certifier
        self.family = family
        self.direction = direction
        self.base = base_certificate
        self.base_atoms = None if base_certificate is None else \
            {row[0] for row in base_certificate["witness"]["atoms"]}
        self.classifier = classifier or SearchClassifier(family)
        self.facts = _fragment_facts()
        self.include_failed = include_failed
        self.calls = 0
        self.flips = []

    def __call__(self, witness) -> Observation:
        self.calls += 1
        mols = template_molecules(witness, self.direction)
        recert = None
        if self.base is None:
            cert = self.certifier.certify(self.family, mols, member_labels=self.wanted)
            self.base = cert
            self.base_atoms = {row[0] for row in cert["witness"]["atoms"]}
        else:
            present = {a.id for m in mols for a in m.atoms}
            amap = {a: a for a in self.base_atoms if a in present}
            recert = self.certifier.recertify_certificate(self.family, self.base, mols, amap)
            cert = recert["certificate"]
            for f in recert.get("flips", []):
                self.flips.append((dict(witness.extension_counts), f["node"], f["old"], f["new"]))
            if recert.get("refusals"):
                return Observation(Support(frozenset()), lambda inside, t=None: None, cert,
                                   "refused:recertification: " + recert["refusals"][0][:240],
                                   recert)
        refusal = certificate_refusal(cert) or member_refusal(cert, self.wanted)
        support = certificate_support(cert, self.wanted, classifier=self.classifier,
                                      fragment_facts=self.facts,
                                      include_failed=self.include_failed)
        refusal = refusal or support.refusal
        missing = sorted(set(support.atoms) - set(witness.reactant_rows()))
        if missing and not refusal:
            refusal = f"refused:transport: support atoms {missing[:4]} are not witness atoms"
        wanted = self.wanted

        def viewer(inside, transport=None, cert=cert):
            return certificate_view(cert, wanted, inside, transport)
        return Observation(support, viewer, cert, refusal, recert)


# ---------------------------------------------------------------------------
# fixed point (item 3)
# ---------------------------------------------------------------------------
def _key(witness):
    return tuple(sorted((k, v) for k, v in witness.extension_counts if v))


def fixed_point(witness, observe, providers=None, *, ceiling=None, port_limit=PORT_DEPTH_LIMIT,
                retries=FIXED_POINT_RETRIES, closure=None, check_ports: bool = True,
                joint_check: bool = False) -> FixedPointResult:
    """The recertified fixed point of ``witness`` under ``observe`` (see the
    module docstring).  ``check_ports``/``joint_check`` exist for the
    mutation tests only."""
    closure = closure or close
    limit = b._ceiling(ceiling)
    ports = b.vector_ports(witness)
    cache = {}
    history = []

    def look(w):
        key = _key(w)
        if key not in cache:
            obs = observe(w)
            win = None if obs.refusal else closure(w, obs.support)
            cache[key] = (obs, win)
            history.append((dict(key), None if win is None else len(win.atoms),
                            obs.refusal or (win.refusal if win else None)))
        return cache[key]

    def depths(w):
        return tuple((p, w.depth(p)) for p in ports)

    def stop(status, reason, w=None, win=None, obs=None, checks=(), used=0):
        return FixedPointResult(status, reason, depths(w or witness), w, win, obs,
                                tuple(checks), used, tuple(history),
                                tuple(getattr(observe, "flips", ())))

    def grow(w, keys):
        for p in keys:
            if w.depth(p) >= port_limit:
                return None, stop("no_fixed_point", f"port {p} reached the per-port limit "
                                  f"{port_limit} without a fixed point", w)
        try:
            for p in keys:
                w, _record = b.extend(w, p, providers, ceiling=limit)
        except b.CeilingExceeded as error:
            return None, stop("ceiling_exceeded", f"{error} (never truncated)", w)
        except b.ExtensionRefused as error:
            if "reacting atom" in str(error):
                return None, stop("refused:reacting_port", f"a reacting atom sits on an "
                                  f"artificial port (a cut-point event: its long-chain and "
                                  f"physical-end analogues are other witnesses): {error}", w)
            return None, stop("refused:extension", f"an artificial port in the window cannot "
                              f"be extended: {error}", w)
        return w, None

    current = witness
    obs, win = look(current)
    if obs.refusal or win.refusal:
        return stop(obs.refusal or win.refusal, "refused", current, win, obs)
    rounds = []
    used = 0
    while True:
        inside = [p for p in ports if port_inside(current, p, win.atoms)]
        while inside:
            current, failure = grow(current, inside)
            if failure:
                return failure
            obs, win = look(current)
            if obs.refusal or win.refusal:
                return stop(obs.refusal or win.refusal, "refused on extension", current, win,
                            obs, rounds, used)
            inside = [p for p in ports if port_inside(current, p, win.atoms)]
        view = obs.view(win.atoms)
        checks = []
        groups = [list(ports)] if joint_check and ports else [[p] for p in ports]
        if not check_ports:
            groups = groups[:1]
        for group in groups:
            other, failure = grow(current, group)
            if failure:
                if failure.status == "ceiling_exceeded":
                    return failure
                return stop(failure.status, failure.reason, current, None, obs, rounds, used)
            o2, w2 = look(other)
            if o2.refusal or w2.refusal:
                return stop("certificate_unstable", f"+1 at {group}: "
                            f"{o2.refusal or w2.refusal}", current, None, obs, rounds, used)
            same_w = same_window(current, win.atoms, other, w2.atoms)
            same_c = view == o2.view(w2.atoms)
            for p in group:
                checks.append({"port": p, "same_window": same_w, "same_certificate": same_c,
                               "window_atoms": len(w2.atoms)})
        rounds.append(tuple(checks))
        failed = [c for c in checks if not (c["same_window"] and c["same_certificate"])]
        if not failed:
            return stop("fixed_point", "", current, win, obs, rounds, used)
        if used >= retries:
            unstable = any(not c["same_certificate"] for c in failed)
            return stop("certificate_unstable" if unstable else "no_fixed_point",
                        f"checks still fail after {used} re-extensions: "
                        f"{[c['port'] for c in failed]}", current, None, obs, rounds, used)
        used += 1
        current, failure = grow(current, sorted({c["port"] for c in failed}))
        if failure:
            return failure
        obs, win = look(current)
        if obs.refusal or win.refusal:
            return stop(obs.refusal or win.refusal, "refused on re-extension", current, win, obs,
                        rounds, used)


# ---------------------------------------------------------------------------
# contexts (item 5)
# ---------------------------------------------------------------------------
def context_assignments(witness, near_ports=None, providers=None, *, cross: bool = True) -> tuple:
    """Assignments ((port, end chemistry or None), ...) over the artificial
    ports, in this order: all artificial; each port alone with each declared
    end; every combination of radical ends over two or more ports (the
    radical inventory is read whole-species, across a merged structure too);
    per participant, the full product over its ``near_ports`` (default: every
    artificial port); with ``cross``, the full product over the near ports of
    all participants (needed when one predicate reads several participants:
    a merged-structure descent)."""
    ports = list(b.vector_ports(witness))
    chem = {p: b.end_chemistries(witness, p, providers) for p in ports}
    near = ports if near_ports is None else [p for p in ports if p in set(near_ports)]
    seen, out = set(), []

    def add(choice):
        assignment = tuple((p, choice.get(p)) for p in ports)
        if assignment not in seen:
            seen.add(assignment)
            out.append(assignment)

    def full(group):
        for size in range(2, len(group) + 1):
            for subset in combinations(group, size):
                for combo in product(*(chem[p] for p in subset)):
                    add(dict(zip(subset, combo)))

    add({})
    for p in ports:
        for e in chem[p]:
            if p in near or e in b.RADICAL_ENDS:
                add({p: e})
    radical = {p: [e for e in chem[p] if e in b.RADICAL_ENDS] for p in ports}
    for size in range(2, len(ports) + 1):
        for subset in combinations(ports, size):
            for combo in product(*(radical[p] for p in subset)):
                add(dict(zip(subset, combo)))
    for index in sorted({p[0] for p in near}):
        full([p for p in near if p[0] == index])
    if cross:
        full(near)
    return tuple(out)


def reconstruct(witness, assignment, providers=None):
    """The witness with every assigned end built (``boundaries.reconstruct_end``)."""
    for port, end in assignment:
        if end is not None:
            witness, _record = b.reconstruct_end(witness, port, end, providers)
    return witness


def contexts_merge(first: ContextResult, second: ContextResult) -> bool:
    if first.status != "fixed_point" or second.status != "fixed_point":
        return False
    w1, w2 = first.window, second.window
    if set(w1.atoms) == set(w2.atoms) and \
            window_labels(first.witness, w1.atoms) == window_labels(second.witness, w2.atoms):
        if first.observation.view(w1.atoms) == second.observation.view(w2.atoms):
            return True
    for mapping in window_isomorphisms(first.witness, w1.atoms, second.witness, w2.atoms):
        if first.observation.view(w1.atoms, {a: a for a in w1.atoms}) == \
                second.observation.view(w2.atoms, mapping):
            return True
    return False


def merge_contexts(results) -> tuple:
    """(representative index, merged index) pairs: a context merges into the
    first earlier class representative whose closed window is equal as a
    marked graph and whose certificate view is equal under that
    correspondence."""
    representatives, merges = [], []
    for j, result in enumerate(results):
        if result.status != "fixed_point":
            continue
        for i in representatives:
            if contexts_merge(results[i], result):
                merges.append((i, j))
                break
        else:
            representatives.append(j)
    return tuple(merges)


# ---------------------------------------------------------------------------
# span bound (item 4) and ladder (item 6)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class SpanBound:
    family: str
    radius: int
    reacting_distance: int
    reach: int
    backbone_atoms: int
    repeats: int
    L: int
    formula: str
    argument: str


SPAN_ARGUMENT = (
    "every support atom lies within r bonds of a reacting atom (balls of the predicate radius "
    "around anchors and bounded candidates, which are reacting atoms; embeddings stay in the "
    "ball); reacting atoms of one chain span at most D bonds; closure adds whole pendant "
    "phenyl rings (no backbone extent), radical pi systems (a backbone carbon and its pendant "
    "ring: allylic radicals are outside the domain) and end alkenes (one backbone bond beyond "
    "the support), so reach = r + 1 per side")


def span_bound(constants: FamilyConstants) -> SpanBound:
    reach = int(constants.radius) + CLOSURE_EXTENSION
    atoms = int(constants.reacting_distance) + 2 * reach + 1
    repeats = math.ceil((atoms + 1) / 2)
    L = (atoms + 1) // 2
    formula = (f"reach = r + {CLOSURE_EXTENSION} = {reach}; backbone atoms A = D + 2*reach + 1 "
               f"= {constants.reacting_distance} + {2 * reach} + 1 = {atoms}; repeats = "
               f"ceil((A + 1)/2) = {repeats}; L = floor((A + 1)/2) = {L} (largest n with "
               f"2n - 1 <= A: one end on each repeat phase; in-phase chains need 2n <= A)")
    return SpanBound(constants.family, int(constants.radius), int(constants.reacting_distance),
                     reach, atoms, repeats, L, formula, SPAN_ARGUMENT)


def chain_backbone_atoms(spec) -> int:
    """Backbone carbons of a PS chain spec (any phase)."""
    if getattr(spec, "monomer", False):
        return 2
    return 2 * spec.units - int(bool(getattr(spec, "head_shift", False))) + \
        int(bool(getattr(spec, "tail_shift", False)))


def family_constants(family) -> FamilyConstants:
    """``r`` (largest predicate radius over the family's tree nodes) and ``D``
    (largest bond distance between labelled heavy atoms of one connected
    template structure, reactants and products after the recipe), read from
    the loaded RMG family."""
    import networkx as nx
    from rmgpy.data.base import LogicNode
    from rmgpy.kmc import channels as chm
    radius = 0
    for entry in family.groups.entries.values():
        if isinstance(entry.item, LogicNode):
            continue
        for group in chm._group_atoms(family.groups, entry):
            r = chm._group_radius(group)
            if r is not None:
                radius = max(radius, r)
    slots = [chm._group_atoms(family.groups, e) for e in family.forward_template.reactants]
    actions = family.forward_recipe.actions
    distance = 0
    free = False
    for combo in product(*slots):
        graph = nx.Graph()
        labels = defaultdict(list)
        heavy = set()
        for s, group in enumerate(combo):
            index = {a: k for k, a in enumerate(group.atoms)}
            for a in group.atoms:
                node = (s, index[a])
                graph.add_node(node)
                if a.label:
                    labels[a.label].append(node)
                    if not all(t.label == "H" for t in a.atomtype):
                        heavy.add(node)
                for c in a.bonds:
                    graph.add_edge(node, (s, index[c]))
        if any(len(v) > 1 for v in labels.values()):
            free = True

        def ends(l1, l2):
            a, c = labels.get(l1, []), labels.get(l2, [])
            if l1 == l2:
                return [tuple(a[:2])] if len(a) >= 2 else []
            return [(x, y) for x in a for y in c]

        product_graph = graph.copy()
        for action in actions:
            if action[0] == "FORM_BOND":
                for x, y in ends(action[1], action[3]):
                    product_graph.add_edge(x, y)
            elif action[0] == "BREAK_BOND":
                for x, y in ends(action[1], action[3]):
                    if product_graph.has_edge(x, y):
                        product_graph.remove_edge(x, y)
        for g in (graph, product_graph):
            for comp in nx.connected_components(g):
                nodes = [n for n in comp if n in heavy]
                sub = g.subgraph(comp)
                for x, y in combinations(nodes, 2):
                    distance = max(distance, nx.shortest_path_length(sub, x, y))
    directions = ("forward",) + (("reverse",) if family.reverse_template is not None else ())
    molecularity = len(family.forward_template.reactants)
    if family.auto_generated:
        root = slots[0][0] if slots and slots[0] else None
        molecularity = 1 if root is None else nx.number_connected_components(
            _group_graph(root))
    reverse_mol = len(family.forward_template.products)
    bimolecular = molecularity == 2 or (family.reverse_template is not None and reverse_mol == 2)
    return FamilyConstants(family.label, radius, distance, free, bimolecular,
                           tuple((d, molecularity if d == "forward" else reverse_mol)
                                 for d in directions))


def _group_graph(group):
    import networkx as nx
    graph = nx.Graph()
    index = {a: k for k, a in enumerate(group.atoms)}
    for a in group.atoms:
        graph.add_node(index[a])
        for c in a.bonds:
            graph.add_edge(index[a], index[c])
    return graph


def backbone_index(provider: str, unit: int, position: int) -> Optional[int]:
    if provider == "polystyrene":
        return 2 * unit + position if position in (0, 1) else None
    if provider == "linear_alkane":
        return unit
    return None


def window_span(witness, atoms) -> dict:
    """Per participant: the window's backbone extent (atoms) and repeat units,
    and the physical ends (port keys) inside the window."""
    atoms = set(atoms)
    out = {}
    for participant in witness.participants:
        meta = participant.metadata
        idx, units = [], []
        for a, u, k in meta.coordinates:
            if a in atoms:
                i = backbone_index(meta.provider, u, k)
                if i is not None:
                    idx.append(i)
                    units.append(u)
        ends = tuple(sorted(port.port_id for port in meta.ports
                            if not port.artificial and port.atom_id in atoms))
        all_ends = tuple(sorted(port.port_id for port in meta.ports if not port.artificial))
        out[participant.index] = {
            "backbone_atoms": (max(idx) - min(idx) + 1) if idx else 0,
            "units": (max(units) - min(units) + 1) if units else 0,
            "physical_ends_inside": ends, "physical_ends": all_ends,
            "both_ends": len(all_ends) == 2 and len(ends) == 2}
    return out


def check_span(witness, atoms, bound: SpanBound) -> dict:
    """``window_span`` with the bound verdict; raises ``SpanBoundError`` when
    a participant's window exceeds the bound."""
    span = window_span(witness, atoms)
    for index, s in span.items():
        if s["backbone_atoms"] > bound.backbone_atoms or s["units"] > bound.repeats:
            raise SpanBoundError(f"{bound.family}: participant {index} window spans "
                                 f"{s['backbone_atoms']} backbone atoms / {s['units']} units "
                                 f"> bound {bound.backbone_atoms} / {bound.repeats}")
    return span


@dataclass(frozen=True)
class LadderItem:
    label: str
    kinds: tuple
    dp: int
    specs: tuple
    control: bool = False
    long_index: Optional[int] = None   # the artificial (long) participant of a pair


@dataclass(frozen=True)
class LadderPlan:
    family: str
    L: int
    items: tuple
    excluded: tuple

    @property
    def excluded_keys(self) -> set:
        return {key for key, _reason in self.excluded}


LONG_UNITS = 4


def _ends_label(spec):
    head = spec.head_end or "cut"
    tail = spec.tail_end or "cut"
    extra = []
    if spec.head_shift:
        extra.append("hs")
    if spec.tail_shift:
        extra.append("ts")
    for name, flag in (("bz", spec.radical_units), ("sec", spec.secondary_radical_units),
                       ("pl", spec.pendant_loss_units)):
        if flag:
            extra.append(f"{name}{','.join(map(str, flag))}")
    return f"DP{spec.units}[{head}|{tail}{'|' + ','.join(extra) if extra else ''}]"


def ladder_plan(constants: FamilyConstants, bound: Optional[SpanBound] = None) -> LadderPlan:
    """The short-chain ladder of a family up to L, plus the DP L + 1 control:
    intact in-phase chains in every declared end chemistry; scission
    fragments of both backbone cut phases (phase A: in-phase fragments with
    the cut end radical or alkene; phase B: shifted fragments); pendant loss
    and C-H cleavage products (one per interior unit position, saturated
    ends); for a bimolecular family, short/short and short/long pairs (one
    representative per local environment: the reacting end chemistry and the
    DP).  Species outside the declared domain are excluded with the
    domain's reason; species built twice are merged (their kinds joined)."""
    bound = bound or span_bound(constants)
    L = bound.L
    P = b.PROVIDERS["polystyrene"]
    items, excluded = {}, []
    by_species = {}

    def species_key(spec):
        from rmgpy.kmc import channels as chm
        build = P.build("ladder", spec)
        verdict = dom.check(build.molecule)
        return build, verdict, chm.SpeciesKeys()(build.molecule)

    def add(key, kind, spec, dp):
        try:
            _build, verdict, skey = species_key(spec)
        except b.BoundaryError as error:
            excluded.append((key, f"not constructible: {error}"))
            return
        if not verdict.accepted:
            excluded.append((key, f"outside the declared domain: {verdict.reason}"))
            return
        if skey in by_species:
            label = by_species[skey]
            item = items[label]
            if kind not in item.kinds:
                items[label] = LadderItem(item.label, item.kinds + (kind,), item.dp, item.specs,
                                          item.control)
            if key[0] != item.label:
                excluded.append((key, f"duplicate species of {label}"))
            return
        label = _ends_label(spec)
        by_species[skey] = label
        items[label] = LadderItem(label, (kind,), dp, (spec,), dp > L)

    for dp in range(1, L + 2):
        for head in b.HEAD_ENDS:
            for tail in b.TAIL_ENDS:
                spec = b.PsChainSpec(dp, head_end=head, tail_end=tail)
                add((dp, head, tail), "intact", spec, dp)
        for head in b.HEAD_ENDS:     # phase A: a cut between repeats
            for tail in ("benzylic_radical", "styryl"):
                add(("A", dp, head, tail), "scission:phase_A",
                    b.PsChainSpec(dp, head_end=head, tail_end=tail), dp)
        for head in ("primary_radical", "vinylidene"):
            for tail in b.TAIL_ENDS:
                add(("A", dp, head, tail), "scission:phase_A",
                    b.PsChainSpec(dp, head_end=head, tail_end=tail), dp)
        for head in b.HEAD_ENDS:     # phase B: a cut inside a repeat
            for tail in b.HEAD_ENDS:
                add(("B", dp, head, tail), "scission:phase_B",
                    b.PsChainSpec(dp, head_end=head, tail_end=tail, tail_shift=True), dp)
        for head in b.TAIL_ENDS:
            for tail in b.TAIL_ENDS:
                add(("B", dp, head, tail), "scission:phase_B",
                    b.PsChainSpec(dp, head_end=head, tail_end=tail, head_shift=True), dp)
        for u in range(dp):
            if 2 * u + 1 < 2 * dp - 1:
                add(("pl", dp, u), "pendant_loss",
                    b.PsChainSpec(dp, head_end="methyl", tail_end="benzyl",
                                  pendant_loss_units=(u,)), dp)
                add(("bz", dp, u), "ch_cleavage:benzylic",
                    b.PsChainSpec(dp, head_end="methyl", tail_end="benzyl", radical_units=(u,)),
                    dp)
            if u >= 1:
                add(("sec", dp, u), "ch_cleavage:secondary",
                    b.PsChainSpec(dp, head_end="methyl", tail_end="benzyl",
                                  secondary_radical_units=(u,)), dp)
    out = list(items.values())
    if constants.bimolecular:
        radical = [("benzylic_radical", lambda d: b.PsChainSpec(d, head_end="methyl",
                                                                tail_end="benzylic_radical")),
                   ("primary_radical", lambda d: b.PsChainSpec(d, head_end="primary_radical",
                                                               tail_end="benzyl"))]
        partners = [("closed", lambda d: b.PsChainSpec(d, head_end="methyl", tail_end="benzyl")),
                    ("vinylidene", lambda d: b.PsChainSpec(d, head_end="vinylidene",
                                                           tail_end="benzyl")),
                    ("styryl", lambda d: b.PsChainSpec(d, head_end="methyl", tail_end="styryl"))]
        partners += radical
        long_closed = b.PsChainSpec(LONG_UNITS)
        long_radical = b.PsChainSpec(LONG_UNITS, radical_units=(LONG_UNITS // 2,))
        for dp in sorted({1, L, L + 1}):
            for rname, rspec in radical:
                for pname, pspec in partners:
                    a, c = rspec(dp), pspec(dp)
                    if not (dom.check(P.build("x", a).molecule).accepted and
                            dom.check(P.build("x", c).molecule).accepted):
                        excluded.append((("pair", dp, rname, pname), "outside the declared "
                                         "domain"))
                        continue
                    out.append(LadderItem(f"pair DP{dp} {rname}+{pname}", ("pair:short_short",),
                                          dp, (a, c), dp > L))
        for dp in range(1, L + 2):
            for rname, rspec in radical:
                a = rspec(dp)
                if not dom.check(P.build("x", a).molecule).accepted:
                    excluded.append((("long", dp, rname), "outside the declared domain"))
                    continue
                out.append(LadderItem(f"pair DP{dp} {rname}+long closed", ("pair:short_long",),
                                      dp, (a, long_closed), dp > L, 1))
                out.append(LadderItem(f"pair DP{dp} {rname}+long radical",
                                      ("pair:short_long",), dp, (a, long_radical), dp > L, 1))
            closed = b.PsChainSpec(dp, head_end="methyl", tail_end="benzyl")
            out.append(LadderItem(f"pair DP{dp} closed+long radical", ("pair:short_long",), dp,
                                  (closed, long_radical), dp > L, 1))
    return LadderPlan(constants.family, L, tuple(out), tuple(excluded))


# ---------------------------------------------------------------------------
# spin and resonance diagnostics (item 7)
# ---------------------------------------------------------------------------
def diagnostics(witness, atoms) -> dict:
    """Radical sites in the window (both sides), the separated-radical verdict
    of every species touched and the matcher's resonance form count per whole
    component used (or the overflow cutoff, never reported as a count)."""
    from rmgpy.kmc import rmg_moves as rm
    from rmgpy.kmc import window_matcher as wm
    atoms = set(atoms)
    out = {"radical_sites": {}, "spin": [], "resonance": []}
    for side in ("reactants", "products"):
        rows = dict(_side_rows(witness, side))
        out["radical_sites"][side] = sorted([a, rows[a][3]] for a in atoms if rows[a][3])
        for part in rm.split_rows(rows):
            if not set(part) & atoms:
                continue
            mol = rm.molecule_from_rows(part)
            verdict = dom.spin_verdict(mol)
            out["spin"].append({"side": side, "first_atom": min(part), "status": verdict.status,
                                "reason": verdict.reason,
                                "radical_electrons": verdict.radical_electrons})

            def value(a, part=part):
                row = part[a]
                return row[0], row[2], row[4], row[3], dict(row[5])
            try:
                _subs, total = wm._component_subs(sorted(part), value)
                out["resonance"].append({"side": side, "first_atom": min(part), "forms": total})
            except wm.ConjugatedRingError as error:
                out["resonance"].append({"side": side, "first_atom": min(part),
                                         "cutoff": f"conjugated ring > "
                                                   f"{wm.CONJUGATED_RING_LIMIT}: {error}"})
            except wm.ResonanceLimitError:
                out["resonance"].append({"side": side, "first_atom": min(part),
                                         "cutoff": f"> {wm.RESONANCE_FORM_LIMIT} forms"})
    return out


# ---------------------------------------------------------------------------
# record members on provider witnesses
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class MemberCase:
    """One member of a channel on a provider witness: ``witness`` the
    representative mapping (record direction, unextended), ``moves`` every
    mapping of the member, ``wanted`` their family-template label sets."""

    label: str
    family: str
    direction: str
    pathway: tuple
    witness: Any = field(compare=False)
    moves: tuple = field(compare=False)
    wanted: frozenset = field(compare=False)
    support: str = "supported"
    species: tuple = ()
    reactants: tuple = field(default=(), compare=False)
    forms: tuple = field(default=(), compare=False)
    kind: str = "kept"


def template_labels(move, direction) -> frozenset:
    if direction == "forward":
        return frozenset((lab, a) for _s, lab, a in move.labels)
    return frozenset(tuple(x) for x in (move.product_labels or ()))


def channel_moves(certifier, family, molecules, direction, resonance="all"):
    """(reactant molecules the moves are on, their forms, channel set)."""
    from rmgpy.kmc import channels as chm
    from rmgpy.kmc import rmg_moves as rm
    if direction == "forward":
        cs = certifier.channel_set(family, molecules)
        forms = tuple(certifier._forms(cs, k) for k in range(len(cs.reactants)))
        return tuple(cs.reactants), forms, cs
    cs = chm.build_channels(family, molecules, "reverse", resonance=resonance)
    forms = tuple(rm.reactant_forms(m, resonance) for m in cs.reactants)
    return tuple(cs.reactants), forms, cs


def _member_moves(cs, channel):
    from rmgpy.kmc import channels as chm
    from rmgpy.kmc import rmg_moves as rm
    by_outcome = defaultdict(list)
    for m in cs.moves:
        if m.status in chm.ALLOWED_STATUSES and m.template is not None and \
                tuple(m.template) == channel.pathway and m.products is not None:
            by_outcome[rm.product_key(m)].append(m)
    return {o: sorted(v, key=rm.Move.sort_key) for o, v in by_outcome.items()}


def ring_atoms(molecules) -> set:
    out = set()
    for m in molecules:
        out |= {a.id for a in m.atoms if m.is_atom_in_cycle(a)}
    return out


def ring_change_select(family, molecules):
    """``select`` refusing a mapping whose bond changes touch a ring atom (a
    ring carbon gains or loses a bond: the declared domain's rings are
    monosubstituted phenyls, so a species of that channel is outside the
    domain: aryl radical, cyclohexadienyl, a second ring substituent or an
    opened ring)."""
    ring = ring_atoms(molecules)

    def select(candidate):
        change = connectivity_change(family, _candidate_labels(candidate))
        return not any(a in ring for _kind, pair in change for a in pair)
    return select


def member_cases(certifier, family, builds, direction, *, label="", molecules=None,
                 participants=None, select=None, skip_ring_changes: bool = True, stats=None):
    """Every member of every channel of ``family`` (``direction``) on the
    provider builds (or given molecules + participants).  With
    ``skip_ring_changes`` mappings whose bond changes touch a ring atom are
    not run through the recipe (``ring_change_select``: refused by the
    domain); ``stats`` counts them."""
    from rmgpy.kmc import certificates as cm
    from rmgpy.kmc import channels as chm
    from rmgpy.kmc import rmg_moves as rm
    if molecules is None:
        molecules, participants = b.witness_reactants(builds)
    roots = [cm.aromatic_root(m) for m in molecules] if direction == "forward" else \
        list(molecules)
    gate = ring_change_select(family, roots) if skip_ring_changes else None
    cs = chm.build_channels(family, roots, direction, resonance=certifier.resonance,
                            select=gate)
    reactants = tuple(cs.reactants)
    forms = tuple(rm.reactant_forms(m, certifier.resonance) for m in reactants)
    if stats is not None:
        stats["ring_change_mappings"] = stats.get("ring_change_mappings", 0) + sum(
            1 for m in cs.moves if m.status == "excluded")
        stats["mappings"] = stats.get("mappings", 0) + len(cs.moves)
        stats["excluded_outcomes"] = stats.get("excluded_outcomes", 0) + len(cs.excluded)
    out = []
    for ch in cs.channels:
        moves_of = _member_moves(cs, ch)
        for k, outcome in enumerate(ch.members):
            moves = tuple(moves_of.get(outcome, ()))
            if not moves:
                continue
            if select is not None and not select(ch, moves):
                continue
            witness = b.witness_from_move(moves[0], reactants, participants, forms=list(forms))
            species = tuple(ch.reactants) + tuple(ch.products)
            out.append(MemberCase(f"{label}|{family.label}|{direction}|{';'.join(ch.pathway)}|"
                                  f"{','.join(ch.products)}|m{k}", family.label, direction,
                                  tuple(ch.pathway), witness, moves,
                                  frozenset(template_labels(m, direction) for m in moves),
                                  ch.support, species, reactants, forms, ch.kind))
    return out


def species_domain_refusal(witness) -> Optional[str]:
    """``domain.check`` of every species on both sides (None: all accepted)."""
    from rmgpy.kmc import rmg_moves as rm
    for side in ("reactants", "products"):
        for part in rm.split_rows(dict(_side_rows(witness, side))):
            verdict = dom.check(rm.molecule_from_rows(part))
            if not verdict.accepted:
                return f"refused:domain: {side} species at atom {min(part)}: {verdict.reason}"
    return None


def mapping_problems(case, window_atoms, fp_witness=None) -> list:
    """Every mapping of the member: reacting atoms (labels, mapped atoms) and
    its edits on BOTH sides inside the window."""
    problems = []
    atoms = set(window_atoms)
    for move in case.moves:
        mw = b.witness_from_move(move, case.reactants, case.witness.participants,
                                 forms=list(case.forms))
        edited = set(b.edit_support(mw))
        labelled = {a for _s, _l, a in move.labels} | {a for _s, _i, _j, a in move.mapped}
        if move.product_labels:
            labelled |= {a for _l, a in move.product_labels}
        outside = sorted((edited | labelled) - atoms)
        if outside:
            problems.append(f"mapping {move.order}/{move.forms}: atoms {outside[:6]} outside "
                            "the window")
    return problems


def run_context(certifier, family, case: MemberCase, assignment, *, bound=None,
                providers=None, start_depths=None, **loop) -> ContextResult:
    """A physical context: ends built by the provider on the ORIGINAL witness
    (the end sits where the witness declares it), the event rematched by
    ``rmg_moves`` (same labels, same mapped outcome), the other artificial
    ports pre-extended to ``start_depths`` (the base fixed point's), then
    recertified and closed to its own fixed point."""
    from rmgpy.kmc import rmg_moves as rm
    try:
        built = reconstruct(case.witness, assignment, providers)
    except b.BoundaryError as error:
        return ContextResult(assignment, "refused:end", str(error))
    why = species_domain_refusal(built)
    if why:
        return ContextResult(assignment, "refused:domain", why, built)
    rows = dict(built.reactant_rows())
    molecules = [rm.molecule_from_rows({a: rows[a] for a in p.atom_ids})
                 for p in built.participants]
    rep = case.moves[0]
    try:
        if case.direction == "forward":
            scoped = scoped_certifier(certifier, family, case.wanted)
            reactants, forms, cs = channel_moves(scoped, family, molecules, "forward")
        else:
            from rmgpy.kmc import channels as chm
            cs = chm.build_channels(family, molecules, "reverse",
                                    select=lambda c: c.labels == rep.labels)
            reactants = tuple(cs.reactants)
            forms = tuple(rm.reactant_forms(m) for m in reactants)
    except Exception as error:   # reported, never hidden
        return ContextResult(assignment, "refused:rematch", f"{type(error).__name__}: {error}",
                             built)
    target = rm.product_key(dict(built.product_rows()))
    found = None
    for ch in cs.channels:
        for outcome, moves in _member_moves(cs, ch).items():
            if outcome == target and any(m.labels == rep.labels for m in moves):
                found = (ch, moves)
                break
        if found:
            break
    if found is None:
        same = [m for m in cs.moves if m.labels == rep.labels]
        return ContextResult(assignment, "refused:rematch",
                             f"no move with the labels and the transported outcome "
                             f"({len(same)} moves with the labels: "
                             f"{sorted({m.status for m in same})})", built)
    _ch, moves = found
    move = next(m for m in moves if m.labels == rep.labels)
    witness = b.witness_from_move(move, reactants, built.participants, forms=list(forms))
    wanted = frozenset(template_labels(m, case.direction) for m in moves)
    for port, depth in (start_depths or ()):
        if witness.port(port).artificial:
            for _k in range(depth):
                try:
                    witness, _record = b.extend(witness, port, providers)
                except b.BoundaryError:
                    break
    observer = MemberObserver(certifier, family, case.direction, wanted)
    fp = fixed_point(witness, observer, providers, **loop)
    if fp.reached and bound is not None:
        check_span(fp.witness, fp.window.atoms, bound)
    return ContextResult(assignment, fp.status, fp.reason, fp.witness, fp.window,
                         fp.observation, fp)


def near_ports(witness, atoms, bound) -> tuple:
    """Artificial ports within ``bound.reach + 1`` backbone bonds of the window."""
    rows = witness.reactant_rows()
    heavy = {a for a, row in rows.items() if row[0] != "H"}
    dist = {a: 0 for a in atoms if a in heavy}
    frontier = list(dist)
    depth = 0
    while frontier and depth < bound.reach + 1:
        depth += 1
        nxt = []
        for a in frontier:
            for n, _o in rows[a][5]:
                if n in heavy and n not in dist:
                    dist[n] = depth
                    nxt.append(n)
        frontier = nxt
    return tuple(p for p in b.vector_ports(witness) if witness.port(p).atom_id in dist)


def run_member(certifier, family, case: MemberCase, *, bound=None, contexts: bool = True,
               near_only: bool = False, providers=None, **loop) -> dict:
    """Items 1-5 and 7 for one member: refusal by status, or the fixed point
    (window, span against ``bound``, per-mapping edit checks, diagnostics,
    recertification flips) and the context list with merges."""
    import time
    start = time.perf_counter()
    out = {"label": case.label, "family": case.family, "direction": case.direction,
           "pathway": list(case.pathway), "mappings": len(case.moves),
           "slot_orders": sorted({m.order for m in case.moves}),
           "resonance_variants": sorted({m.forms for m in case.moves})}
    if case.support != "supported":
        out.update(status=f"refused:{case.support}", reason="channel support")
        return out
    why = species_domain_refusal(case.witness)
    if why:
        out.update(status="refused:domain", reason=why)
        return out
    observer = MemberObserver(certifier, family, case.direction, case.wanted)
    fp = fixed_point(case.witness, observer, providers, **loop)
    out.update(status=fp.status, reason=fp.reason, depths=[[list(p), d] for p, d in fp.depths],
               retries=fp.retries, observations=observer.calls,
               flips=sorted({(n, o, w) for _k, n, o, w in fp.flips}),
               check_rounds=len(fp.checks))
    if not fp.reached:
        if contexts and fp.status == "refused:reacting_port":
            _contexts(out, certifier, family, case, fp, None, bound, providers, loop)
        out["seconds"] = round(time.perf_counter() - start, 1)
        return out
    window = fp.window
    out["window_atoms"] = len(window.atoms)
    out["window_heavy"] = sum(1 for a in window.atoms if fp.witness.reactant_rows()[a][0] != "H")
    out["boundary_bond_problems"] = boundary_bond_problems(fp.witness, window)
    out["mapping_problems"] = mapping_problems(case, window.atoms)
    span = window_span(fp.witness, window.atoms)
    out["span"] = {str(k): v for k, v in span.items()}
    if bound is not None:
        check_span(fp.witness, window.atoms, bound)
        out["span_bound_held"] = True
    out["dispositions"] = sorted({(str(k[0]), d) for k, d, _a in window.support.dispositions})
    out["diagnostics"] = diagnostics(fp.witness, window.atoms)
    if contexts:
        near = near_ports(case.witness, window.atoms, bound) if (near_only and bound) else None
        _contexts(out, certifier, family, case, fp, near, bound, providers, loop)
    out["seconds"] = round(time.perf_counter() - start, 1)
    return out


def merged_structure(certificate) -> bool:
    """Whether some descent reads the merged reactants (one predicate can see
    several participants)."""
    if not certificate:
        return True
    return any(d.get("structure") == "merged reactants"
               for ds in certificate.get("descents", {}).values() for d in ds)


def _contexts(out, certifier, family, case, fp, near, bound, providers, loop):
    results = [ContextResult(tuple((p, None) for p in b.vector_ports(case.witness)),
                             fp.status, fp.reason, fp.witness, fp.window, fp.observation, fp)]
    cross = merged_structure(fp.observation.certificate if fp.observation else None)
    out["context_policy"] = {"near_ports": None if near is None else [list(p) for p in near],
                             "cross_participant": cross}
    for assignment in context_assignments(case.witness, near, cross=cross)[1:]:
        results.append(run_context(certifier, family, case, assignment, bound=bound,
                                   providers=providers, start_depths=fp.depths, **loop))
    merges = merge_contexts(results)
    merged = {j for _i, j in merges}
    rows = []
    for k, r in enumerate(results):
        span = window_span(r.witness, r.window.atoms) if r.window and r.status == "fixed_point" \
            else None
        rows.append({"assignment": [[list(p), e] for p, e in r.assignment if e],
                     "status": r.status, "reason": r.reason[:200],
                     "window_atoms": len(r.window.atoms) if span else None,
                     "both_ends": any(s["both_ends"] for s in span.values()) if span else None,
                     "physical_ends_inside": {str(i): s["physical_ends_inside"]
                                              for i, s in span.items()} if span else None,
                     "merged_into": next((i for i, j in merges if j == k), None)})
    out["contexts"] = rows
    out["context_classes"] = len([k for k, r in enumerate(results)
                                  if r.status == "fixed_point" and k not in merged])
    out["merges"] = [list(x) for x in merges]
