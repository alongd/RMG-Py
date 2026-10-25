"""The initial melt of a branched declaration (``BranchedPolymerSpec``):
``random_bond_v1`` and the molecules it builds.

**random_bond_v1.**  The residue instances are the declared representative
counts, in residue-name then index order.  Each link of each type (link name
order, then index) draws a free site of its source class, then one of its
target class, each uniformly among the free sites of that class (instance,
then site name, order): the draw is ``u = int(sha256(random_bond_v1:<seed>:
<link>:<index>:<attempt>:<end>)[:8]) / 2**64``.  A draw is refused and drawn
again (the attempt counter advances) when it links an instance to itself, adds
a second link between one pair of instances, or closes a cycle of fewer than
``minimum_cycle_bonds`` bonds (the shortest path between its two site atoms,
over the residues' heavy-atom paths between their site atoms and the links
already drawn, plus the new bond); ``refused:random_bond`` after
``MAX_ATTEMPTS`` refused draws of one link.  The realized link counts are the
declared ones by construction.  The declaration has already refused counts at
or above the gel point (``polymer_spec.gel_criterion``).

**Molecules** are the connected components, in order of their least instance;
a molecule's units are its instances in instance order, so a unit's segment is
``unit <k> <residue>`` (k from 1, as for a linear v2 chain) and the molecule is
a ``chain`` origin ordinal.  A vacant site is capped by one H; ``is_backbone``
flags the atoms on the shortest paths between a residue's site atoms.

``realized`` reports the realized melt against the declaration: the link
counts, the branch points (instances with three or more linked sites) against
their exact expectation under uniform site occupancy (each class's link ends
spread uniformly over its sites without replacement, classes independent), the
molecule sizes in residues and the spectral radius of the half-edge branching
matrix (below 1: below the gel point)."""

from __future__ import annotations

import hashlib
import heapq
import weakref
from collections import Counter
from dataclasses import dataclass, replace
from fractions import Fraction
from math import comb
from typing import Dict, Mapping, Sequence, Tuple

from rmgpy import constants
from rmgpy.kmc.melt_state import MeltStateError, _Builder, _path_keys, segment_residue

RANDOM_BOND_REFUSAL = "refused:random_bond"
MAX_ATTEMPTS = 10 ** 4


@dataclass(frozen=True)
class RandomBond:
    """``instances``: residue name per instance; ``links``: (instance, site,
    instance, site) per drawn link, source end first, with ``types`` the link
    name of each; ``refused``: draws refused and drawn again; ``molecules``:
    the instances of each connected component."""

    instances: Tuple[str, ...]
    links: Tuple[Tuple[int, str, int, str], ...]
    types: Tuple[str, ...]
    refused: int
    molecules: Tuple[Tuple[int, ...], ...]


def _draw(seed, link, index, attempt, end, n) -> int:
    digest = hashlib.sha256(f"random_bond_v1:{seed}:{link}:{index}:{attempt}:{end}".encode()).digest()
    return (int.from_bytes(digest[:8], "big") * n) >> 64


def _site_distances(residue) -> Dict[Tuple[str, str], int]:
    """Heavy-atom bond distances between the site atoms of a residue."""
    graph = residue.graph
    heavy = {a.key for a in graph.atoms if a.element != "H"}
    adjacency = {k: set() for k in heavy}
    for bond in graph.bonds:
        p, q = bond.atoms
        if p in heavy and q in heavy:
            adjacency[p].add(q)
            adjacency[q].add(p)
    atoms = sorted({site.atom for _s, site in residue.sites})
    out = {}
    for a in atoms:
        distance, frontier = {a: 0}, [a]
        while frontier:
            following = []
            for x in frontier:
                for y in adjacency[x]:
                    if y not in distance:
                        distance[y] = distance[x] + 1
                        following.append(y)
            frontier = following
        for b in atoms:
            out[(a, b)] = distance[b]
    return out


def _path_bonds(start, goal, adjacency) -> int:
    """Fewest bonds from ``start`` to ``goal`` (None: not connected)."""
    best, heap = {start: 0}, [(0, start)]
    while heap:
        d, node = heapq.heappop(heap)
        if node == goal:
            return d
        if d > best[node]:
            continue
        for nxt, w in adjacency.get(node, ()):
            if d + w < best.get(nxt, d + w + 1):
                best[nxt] = d + w
                heapq.heappush(heap, (d + w, nxt))
    return None


def random_bond(spec) -> RandomBond:
    """``random_bond_v1`` on the declared counts (module doc)."""
    residues = dict(spec.residues)
    ensemble = spec.initial_ensemble
    counts = dict(ensemble.residue_counts)
    instances = tuple(name for name in sorted(residues) for _ in range(counts[name]))
    sites = {name: dict(residue.sites) for name, residue in residues.items()}
    distances = {name: _site_distances(residue) for name, residue in residues.items()}
    free = {}
    for i, name in enumerate(instances):
        for s, site in sorted(sites[name].items()):
            free.setdefault(site.site_class, []).append((i, s))
    adjacency = {}
    for i, name in enumerate(instances):
        atoms = sorted({site.atom for site in sites[name].values()})
        for a in atoms:
            adjacency[(i, a)] = [((i, b), distances[name][(a, b)]) for b in atoms if b != a]
    parent = list(range(len(instances)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    linked, links, types, refused = set(), [], [], 0
    for lname, link in spec.links:
        for index in range(dict(ensemble.link_counts)[lname]):
            for attempt in range(MAX_ATTEMPTS + 1):
                if attempt == MAX_ATTEMPTS:
                    raise MeltStateError(
                        f"{RANDOM_BOND_REFUSAL}: link {lname} #{index} was refused {MAX_ATTEMPTS} times (a "
                        "self-link, a second link between one pair of residues or a cycle shorter than "
                        f"{spec.minimum_cycle_bonds} bonds): the declared counts admit no such melt")
                ends = []
                for end, cls in enumerate((link.source, link.target)):
                    pool = [x for x in free.get(cls, ()) if x not in ends]
                    if not pool:
                        raise MeltStateError(f"{RANDOM_BOND_REFUSAL}: no free site of class {cls!r} is left "
                                             f"for link {lname} #{index}")
                    ends.append(pool[_draw(ensemble.seed, lname, index, attempt, end, len(pool))])
                (i, s), (j, t) = ends
                a, b = (i, sites[instances[i]][s].atom), (j, sites[instances[j]][t].atom)
                if i == j or frozenset((i, j)) in linked:
                    refused += 1
                    continue
                if root(i) == root(j):
                    path = _path_bonds(a, b, adjacency)
                    if path is not None and path + 1 < spec.minimum_cycle_bonds:
                        refused += 1
                        continue
                break
            for x in ends:
                free[sites[instances[x[0]]][x[1]].site_class].remove(x)
            linked.add(frozenset((i, j)))
            adjacency[a].append((b, 1))
            adjacency[b].append((a, 1))
            parent[root(i)] = root(j)
            links.append((i, s, j, t))
            types.append(lname)
    groups = {}
    for i in range(len(instances)):
        groups.setdefault(root(i), []).append(i)
    molecules = tuple(sorted(tuple(sorted(g)) for g in groups.values()))
    return RandomBond(instances, tuple(links), tuple(types), refused, molecules)


def backbone_keys(residue) -> frozenset:
    """Atom keys on the shortest paths between the residue's site atoms."""
    atoms = sorted({site.atom for _s, site in residue.sites})
    keys = set(atoms)
    for x in atoms:
        for y in atoms:
            if x < y:
                keys |= _path_keys(residue.graph, x, y)
    return frozenset(keys)


def _capped(residue, linked_sites) -> object:
    """The residue graph with one H on each site not in ``linked_sites``."""
    extra = Counter(site.atom for s, site in residue.sites if s not in linked_sites)
    atoms = tuple(replace(a, hydrogens=a.hydrogens + extra[a.key]) for a in residue.graph.atoms)
    return replace(residue.graph, atoms=atoms)


def _add_molecules(spec, molecules, builder, first=0) -> int:
    residues = dict(spec.residues)
    for ordinal, (names, links) in enumerate(molecules, start=first):
        linked = {}
        for i, s, j, t in links:
            for unit, site in ((i, s), (j, t)):
                if site not in dict(residues[names[unit]].sites):
                    raise MeltStateError(f"residue {names[unit]} has no site {site!r}")
                if site in linked.setdefault(unit, set()):
                    raise MeltStateError(f"site {names[unit]}.{site} of unit {unit} is linked twice")
                linked[unit].add(site)
        ids = []
        for unit, name in enumerate(names):
            residue = residues[name]
            vacant = [s for s, site in residue.sites if s not in linked.get(unit, ()) and site.vacancy != "H"]
            if vacant:
                raise MeltStateError(f"the forbidden-vacant sites {vacant} of {name} (unit {unit}) are not linked")
            ids.append(builder.add_segment(_capped(residue, linked.get(unit, ())), "chain", ordinal,
                                           f"unit {unit + 1} {name}", backbone_keys(residue)))
        for i, s, j, t in links:
            si, sj = dict(residues[names[i]].sites)[s], dict(residues[names[j]].sites)[t]
            builder.bond(ids[i][si.atom], ids[j][sj.atom], si.bond_order)
    return first + len(molecules)


def molecule_words(result: RandomBond):
    """(residue names, links in molecule-local units) of each molecule of a
    ``random_bond`` result."""
    out = []
    for members in result.molecules:
        local = {i: k for k, i in enumerate(members)}
        links = [(local[i], s, local[j], t) for i, s, j, t in result.links if i in local]
        out.append(([result.instances[i] for i in members], links))
    return out


def chains(spec, result: RandomBond = None):
    """(``_Builder``, molecule count) of the ``random_bond_v1`` melt."""
    result = random_bond(spec) if result is None else result
    builder = _Builder()
    return builder, _add_molecules(spec, molecule_words(result), builder)


def build_molecules(spec, molecules: Sequence[Tuple[Sequence[str], Sequence[Tuple]]]):
    """The ``MeltState`` of given molecules: (residue names, links (unit, site,
    unit, site)) each, built as ``materialize`` builds the random melt."""
    builder = _Builder()
    _add_molecules(spec, molecules, builder)
    return builder.state()


def residue_atoms(state) -> Dict[Tuple[int, int], Dict[str, int]]:
    """(molecule ordinal, unit from 0) -> declaration atom key -> atom id, for
    the heavy atoms of every residue unit."""
    out = {}
    for a, origin in state.origin.items():
        if origin.kind == "chain" and segment_residue(origin.segment) and state.atoms[a].element != "H":
            unit = int(origin.segment.split(" ")[1]) - 1
            out.setdefault((origin.ordinal, unit), {})[origin.key] = a
    return out


def _hypergeometric(k, sites, ends, j) -> Fraction:
    return Fraction(comb(k, j) * comb(sites - k, ends - j), comb(sites, ends))


def expected_branch_points(spec) -> Fraction:
    """Instances with three or more linked sites expected under uniform site
    occupancy (module doc), exact."""
    from rmgpy.polymer_spec import _class_tables
    counts, link_counts = dict(spec.initial_ensemble.residue_counts), dict(spec.initial_ensemble.link_counts)
    sites, ends, _pairs = _class_tables(spec.residues, spec.links, counts, link_counts)
    total = Fraction(0)
    for name, residue in spec.residues:
        per = Counter(site.site_class for _s, site in residue.sites)
        dist = {0: Fraction(1)}
        for cls, k in sorted(per.items()):
            nxt = {}
            for j in range(k + 1):
                p = _hypergeometric(k, sites[cls], ends.get(cls, 0), j)
                for have, q in dist.items():
                    nxt[have + j] = nxt.get(have + j, Fraction(0)) + p * q
            dist = nxt
        total += counts[name] * sum((p for linked, p in dist.items() if linked >= 3), Fraction(0))
    return total


def realized(spec, result: RandomBond = None) -> dict:
    """The realized melt against the declaration (module doc)."""
    import numpy as np
    from rmgpy.polymer_spec import gel_criterion
    result = random_bond(spec) if result is None else result
    linked = Counter()
    for i, _s, j, _t in result.links:
        linked[i] += 1
        linked[j] += 1
    counts, link_counts = dict(spec.initial_ensemble.residue_counts), dict(spec.initial_ensemble.link_counts)
    _below, matrix = gel_criterion(spec.residues, spec.links, counts, link_counts)
    radius = max(abs(np.linalg.eigvals(np.array(matrix, dtype=float)))) if matrix else 0.0
    return {
        "generator": spec.initial_ensemble.generator,
        "seed": spec.initial_ensemble.seed,
        "link_counts": {name: Counter(result.types)[name] for name in sorted(link_counts)},
        "declared_link_counts": dict(sorted(link_counts.items())),
        "refused_draws": result.refused,
        "branch_points": sum(1 for n in linked.values() if n >= 3),
        "expected_branch_points": expected_branch_points(spec),
        "molecule_sizes": dict(sorted(Counter(len(m) for m in result.molecules).items())),
        "gel_branching_ratio": float(radius),
    }


def _residue_mass(residue) -> Fraction:
    from rmgpy.kmc.melt_sequence import _graph_mass
    return _graph_mass(residue.graph)


def mean_residue_molar_mass(spec) -> Fraction:
    """The melt's exact mean residue molar mass (g/mol): every residue with its
    sites occupied, plus one H per vacant site (sites minus twice the links),
    over the residue count."""
    from rmgpy.kmc.melt_sequence import _graph_mass
    from rmgpy.polymer_spec import _vacancy_cap
    counts = dict(spec.initial_ensemble.residue_counts)
    residues = sum(counts.values())
    if not residues:
        raise MeltStateError("a branched declaration with no residues has no mean residue mass")
    mass = sum((counts[name] * _residue_mass(residue) for name, residue in spec.residues), Fraction(0))
    vacant = sum(counts[name] * len(residue.sites) for name, residue in spec.residues) \
        - 2 * sum(dict(spec.initial_ensemble.link_counts).values())
    return (mass + vacant * _graph_mass(_vacancy_cap())) / residues


@dataclass(frozen=True)
class BranchedInitialization:
    """The random melt's representative masses, volume and scale factor (the
    v1 ``SharedInitialization`` fields the melt kMC reads)."""

    molecules: Tuple[Tuple[int, ...], ...]
    representative_polymer_mass_kg: float
    representative_mobile_mass_kg: float
    volume_m3: float
    macroscopic_polymer_mass_kg: float
    scale_factor: float


def initialization(spec, result: RandomBond = None) -> BranchedInitialization:
    from rmgpy.kmc.melt_sequence import _graph_mass
    result = random_bond(spec) if result is None else result
    polymer = mean_residue_molar_mass(spec) * len(result.instances)
    small = {name: graph for name, graph, _disposition in spec.small_species}
    mobile = sum((count * _graph_mass(small[name]) for name, count in spec.initial_mobile_counts), Fraction(0))
    polymer_kg = float(polymer / 1000) / constants.Na
    mobile_kg = float(mobile / 1000) / constants.Na
    total_kg = float((polymer + mobile) / 1000) / constants.Na
    return BranchedInitialization(
        result.molecules, polymer_kg, mobile_kg, total_kg / spec.density_kg_m3, spec.initial_mass_kg,
        spec.initial_mass_kg / polymer_kg if polymer_kg else 0.0)


def link_classes(spec) -> frozenset:
    """The declared linkable class pairs, both orders."""
    return frozenset(pair for _n, link in spec.links
                     for pair in ((link.source, link.target), (link.target, link.source)))


def residue_site_table(spec) -> Mapping[str, Mapping[str, Tuple[str, ...]]]:
    """residue -> site atom key -> the classes of its sites (one per site)."""
    out = {}
    for name, residue in spec.residues:
        for _s, site in residue.sites:
            out.setdefault(name, {}).setdefault(site.atom, []).append(site.site_class)
    return {name: {atom: tuple(sorted(c)) for atom, c in table.items()} for name, table in out.items()}


# --------------------------------------------------------------------------
# declared linearity (the melt engine's leak rule for a branched declaration)
# --------------------------------------------------------------------------
DECLARED = "declared"


def _occurrence(state, a):
    """(ordinal, segment, residue) of a heavy residue atom, else None."""
    origin = state.origin[a]
    if origin.kind != "chain" or state.atoms[a].element == "H":
        return None
    residue = segment_residue(origin.segment)
    return None if residue is None else (origin.ordinal, origin.segment, residue)


class DeclaredLinearity:
    """The declared linearity rule of a ``BranchedPolymerSpec`` (it replaces
    the linear guard of ``melt_engine`` for that declaration only).  A formed
    bond between heavy atoms of two residue units is allowed iff
    1. each atom is a site atom with a site left vacant (its inter-residue
       bonds, after the transition's broken and formed bonds, at most its
       sites) and a class of one pairs with a class of the other in a
       declared link;
    2. each unit keeps at most its site count of inter-residue bonds;
    3. it closes no cycle of fewer than ``minimum_cycle_bonds`` bonds (the
       shortest heavy-atom path between its atoms, the transition's other
       bonds applied, plus the bond).
    Every other inter-residue bond is leak.  A bond at a heavy atom of no
    residue unit (a small species, a released fragment) is leak when it joins
    two pieces that each carry a residue atom (an undeclared link), and is
    allowed otherwise; a bond inside one unit is allowed.

    A bimolecular site's kind is ``(internal leak, items)``: per cross-bond
    atom ``(template atom, DECLARED, unit label, vacant classes, linkable
    partner classes, free sites, spare unit sites, carries a residue)`` read
    after the participant's own edits, the unit label the least template atom
    of that unit among the site's cross atoms; two participants are on two
    components, so a pair closes no cycle (``declared_pair_leak``)."""

    def __init__(self, spec):
        self.sites = residue_site_table(spec)
        self.capacity = {name: len(residue.sites) for name, residue in spec.residues}
        self.links = link_classes(spec)
        self.minimum = int(spec.minimum_cycle_bonds)
        self._index = weakref.WeakKeyDictionary()

    def _units(self, state):
        index = self._index.get(state)
        if index is None:
            index = {}
            for a in state.atoms:
                occ = _occurrence(state, a)
                if occ is not None and state.origin[a].key in self.sites[occ[2]]:
                    index.setdefault(occ, []).append(a)
            index = self._index[state] = {k: tuple(v) for k, v in index.items()}
        return index

    def _inter(self, state, a, removed, added=()):
        """Inter-residue bonds of heavy atom ``a`` after ``removed`` (atom ->
        partners) and ``added`` (partners)."""
        occ = _occurrence(state, a)
        cut = removed.get(a, ())
        partners = [b for b in state.atoms[a].edges if b not in cut] + list(added)
        return sum(1 for b in partners if _occurrence(state, b) not in (None, occ))

    def _unit_links(self, state, occ, removed, formed):
        return sum(self._inter(state, a, removed, formed.get(a, ())) for a in self._units(state).get(occ, ()))

    def _carries(self, state, a, removed) -> bool:
        seen, stack = {a}, [a]
        while stack:
            x = stack.pop()
            if _occurrence(state, x) is not None:
                return True
            cut = removed.get(x, ())
            for y in state.atoms[x].edges:
                if y not in seen and y not in cut and state.atoms[y].element != "H":
                    seen.add(y)
                    stack.append(y)
        return False

    def _path_within(self, state, a, b, removed, extra, limit) -> bool:
        """Whether a heavy-atom path of at most ``limit`` bonds joins a to b."""
        distance, frontier = {a: 0}, [a]
        while frontier:
            following = []
            for x in frontier:
                if distance[x] >= limit:
                    continue
                cut = removed.get(x, ())
                for y in list(state.atoms[x].edges) + list(extra.get(x, ())):
                    if y in cut or y in distance or state.atoms[y].element == "H":
                        continue
                    if y == b:
                        return True
                    distance[y] = distance[x] + 1
                    following.append(y)
            frontier = following
        return False

    def _site_ok(self, state, a, removed, formed) -> Tuple[tuple, int]:
        """(classes of ``a``'s sites, free sites after the edits)."""
        occ = _occurrence(state, a)
        classes = self.sites[occ[2]].get(state.origin[a].key, ())
        return classes, len(classes) - self._inter(state, a, removed, formed.get(a, ()))

    def bond_leaks(self, state, a, b, removed, formed) -> bool:
        if "H" in (state.atoms[a].element, state.atoms[b].element):
            return False
        oa, ob = _occurrence(state, a), _occurrence(state, b)
        if oa is not None and oa == ob:
            return False
        if oa is None or ob is None:
            extra = {x: {y for y in ys if {x, y} != {a, b}} for x, ys in formed.items()}
            return not self._path_within(state, a, b, removed, extra, len(state.atoms)) \
                and self._carries(state, a, removed) and self._carries(state, b, removed)
        (ca, fa), (cb, fb) = self._site_ok(state, a, removed, formed), self._site_ok(state, b, removed, formed)
        if not ca or not cb or fa < 0 or fb < 0:
            return True
        if not any((x, y) in self.links for x in ca for y in cb):
            return True
        if any(self._unit_links(state, o, removed, formed) > self.capacity[o[2]] for o in (oa, ob)):
            return True
        extra = {x: {y for y in ys if {x, y} != {a, b}} for x, ys in formed.items()}
        return self._path_within(state, a, b, removed, extra, self.minimum - 2)

    def leaks(self, state, edits) -> bool:
        """Whether a transition (concrete ``SetBond`` edits) forms a bond the
        declaration does not allow (class doc)."""
        from rmgpy.kmc.melt_state import SetBond
        formed, broken = {}, {}
        for e in edits:
            if isinstance(e, SetBond) and (e.old == 0) != (e.new == 0):
                target = formed if e.new != 0 else broken
                target.setdefault(e.a, set()).add(e.b)
                target.setdefault(e.b, set()).add(e.a)
        return any(a < b and self.bond_leaks(state, a, b, broken, formed)
                   for a, partners in formed.items() for b in partners)

    def site_kind(self, state, cls, t, participant: int, site) -> Tuple:
        pos = {ti: a for ti, a in zip(t.edited[participant], site)}
        removed, formed = {}, {}
        for (i, j), old, new in cls.internal:
            if t.atoms[i].participant != participant:
                continue
            a, b = pos[i], pos[j]
            if (old == 0) != (new == 0):
                target = formed if new != 0 else removed
                target.setdefault(a, set()).add(b)
                target.setdefault(b, set()).add(a)
        internal_leak = any(a < b and self.bond_leaks(state, a, b, removed, formed)
                            for a, partners in formed.items() for b in partners)
        cross = [i0 if participant == 0 else i1 for i0, i1 in cls.cross_bonds]
        label = {}
        for ti in sorted(cross):
            label.setdefault(_occurrence(state, pos[ti]), ti)
        items = []
        for ti in cross:
            a = pos[ti]
            occ = _occurrence(state, a)
            if occ is None:
                items.append((ti, DECLARED, -1, (), (), 0, 0, self._carries(state, a, removed)))
                continue
            classes, free = self._site_ok(state, a, removed, formed)
            vacant = classes if free > 0 else ()
            partners = tuple(sorted({y for x in vacant for (p, y) in self.links if p == x}))
            spare = self.capacity[occ[2]] - self._unit_links(state, occ, removed, formed)
            items.append((ti, DECLARED, label[occ], vacant, partners, free, spare, True))
        return (internal_leak, tuple(sorted(set(items))))


def is_declared_kind(kind) -> bool:
    return bool(kind[1]) and kind[1][0][1] == DECLARED


def declared_pair_leak(cls, kind0, kind1) -> bool:
    """The pair verdict of two declared kinds (``DeclaredLinearity``)."""
    if kind0[0] or kind1[0]:
        return True
    info = {}
    for side, kind in ((0, kind0), (1, kind1)):
        for ti, _tag, unit, vacant, partners, free, spare, carries in kind[1]:
            info[ti] = (side, unit, vacant, partners, free, spare, carries)
    atom_gain, unit_gain = Counter(), Counter()
    for i0, i1 in cls.cross_bonds:
        x, y = info[i0], info[i1]
        if x[1] < 0 or y[1] < 0:
            # a bond at an atom of no residue unit: leak when both pieces carry a residue atom
            if x[6] and y[6]:
                return True
            continue
        if not set(x[3]) & set(y[2]):
            return True
        for ti, item in ((i0, x), (i1, y)):
            atom_gain[ti] += 1
            unit_gain[(item[0], item[1])] += 1
    for ti, gain in atom_gain.items():
        if gain > info[ti][4]:
            return True
    spare = {(item[0], item[1]): item[5] for item in info.values() if item[1] >= 0}
    return any(gain > spare[key] for key, gain in unit_gain.items())
