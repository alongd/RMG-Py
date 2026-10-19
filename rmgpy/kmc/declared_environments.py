"""The environment set of a linear copolymer declaration (schema v2).

**Why one period suffices.**
- A site-class environment (``site_class.extract_site_class`` at the
  declared locality radius r, plus its r + 1 boundary shell) sees 2(r + 1) + 1
  consecutive backbone atoms.
- Those atoms touch at most k = 1 + ceil(2(r + 1) / n_min) residues, where
  n_min is the shortest residue backbone (atoms on the shortest path from the
  entry site's atom to the exit site's atom).
- So one periodic chain whose period holds every k-residue word the sequence
  support allows shows every environment any chain of the declaration can
  have.
- The period is a closed walk over the support's k-words
  (``covering_period``).
- The set depends on the support alone, never on composition or statistics,
  just as the chemistry digest does.

**Bond orders are the declared ones** (aromatic 1.5), as in
``declared_chain._atom_classes``. A Kekulé chain would split the classes of
aromatic atoms that the declared orders keep equivalent.

**Keys compare across declarations and periods.** Every key comes from one
site-class policy with a fixed ``max_heavy_atoms`` (``MAX_HEAVY_ATOMS``).

**Each key is computed on a ball, not the whole chain.**
- The ball has heavy-atom radius r + 3 around the atom.
- A cut single bond gets a hydrogen cap; at a multiple bond the ball grows
  across it instead.
- Reason: ``Molecule.update`` is super-linear in ring count, so whole chains
  with many rings are too slow to extract.
- ``whole_chain=True`` extracts on the whole chain; tests use it to show the
  two keys are equal.
"""

from collections import deque
from fractions import Fraction
from math import ceil
from typing import Dict, List, Mapping, Sequence, Tuple

from rmgpy.kmc import declared_chain as dc

MAX_HEAVY_ATOMS = 10 ** 6
BALL_EXTRA = 3


def word_order(radius: int, n_min: int) -> int:
    """k: the most residues one environment of radius ``radius`` (with its
    boundary shell) can touch when the shortest residue backbone has
    ``n_min`` atoms."""
    return 1 + ceil(2 * (int(radius) + 1) / int(n_min))


def _ends(junctions: Mapping) -> Tuple[Dict[str, str], Dict[str, str]]:
    """(entry atom, exit atom) of every residue named by ``junctions``."""
    entry, exit_ = {}, {}
    for (a, b), (left, right, _order) in junctions.items():
        exit_[a] = left
        entry[b] = right
    return entry, exit_


def _backbone_atoms(residue, entry: str, exit_: str) -> int:
    atoms = {a.key: a for a in residue.graph.atoms if a.element != "H"}
    adjacency = {k: set() for k in atoms}
    for bond in residue.graph.bonds:
        x, y = bond.atoms
        if x in atoms and y in atoms:
            adjacency[x].add(y)
            adjacency[y].add(x)
    distance, queue = {entry: 0}, deque([entry])
    while queue:
        key = queue.popleft()
        for nxt in adjacency[key]:
            if nxt not in distance:
                distance[nxt] = distance[key] + 1
                queue.append(nxt)
    return distance[exit_] + 1


def spec_word_order(spec) -> int:
    """k of a ``PolymerSpec``: its locality radius and shortest residue backbone."""
    entry, exit_ = _ends(spec.junctions)
    n_min = min(_backbone_atoms(residue, entry[name], exit_[name]) for name, residue in spec.residues)
    return word_order(spec.locality_policy["radius"], n_min)


def _de_bruijn(n: int, k: int) -> List[int]:
    """The lexicographically least de Bruijn sequence B(n, k) (FKM)."""
    a, out = [0] * n * k, []

    def extend(t, p):
        if t > k:
            if k % p == 0:
                out.extend(a[1:p + 1])
            return
        a[t] = a[t - p]
        extend(t + 1, p)
        for j in range(a[t - p] + 1, n):
            a[t] = j
            extend(t + 1, t)
    extend(1, 1)
    return out


def covering_period(support: Mapping[str, Sequence[str]], k: int) -> Tuple[str, ...]:
    """A cyclic residue sequence whose cyclic k-windows are exactly the
    k-words the support allows (every one at least once).

    - **Complete support** (every residue may follow every residue): the
      lexicographically least de Bruijn sequence B(n, k). Every k-word
      appears exactly once.
    - **Otherwise:** a closed walk over the graph whose nodes are
      (k-1)-words and whose edges are k-words. The walk takes an unused edge
      when one leaves the current node, and otherwise moves to the nearest
      node that has one. The walk exists when the support is strongly
      connected, which ``PolymerSpec`` requires."""
    if int(k) < 2:
        raise ValueError(f"a covering period needs k >= 2, got {k}")
    names = sorted(support)
    successors = {a: tuple(sorted(support[a])) for a in names}
    if len(names) == 1 and successors[names[0]] == tuple(names):
        return tuple(names)
    if all(successors[a] == tuple(names) for a in names):
        return tuple(names[i] for i in _de_bruijn(len(names), int(k)))
    nodes = [(a,) for a in names]
    for _ in range(k - 2):
        nodes = [w + (b,) for w in nodes for b in successors[w[-1]]]
    edges = {node: [node[1:] + (b,) for b in successors[node[-1]]] for node in nodes}
    unused = {node: list(targets) for node, targets in edges.items()}
    remaining = sum(len(targets) for targets in unused.values())

    def path_to(origin, goal):
        previous, queue = {origin: None}, deque([origin])
        while queue:
            node = queue.popleft()
            if node != origin and goal(node):
                path = [node]
                while previous[path[-1]] is not None:
                    path.append(previous[path[-1]])
                return path[::-1]
            for nxt in edges[node]:
                if nxt not in previous:
                    previous[nxt] = node
                    queue.append(nxt)
        raise ValueError("the sequence support is not strongly connected: no single chain period "
                         "holds every word")

    start = current = nodes[0]
    walk = [start]
    while remaining:
        if unused[current]:
            current = unused[current].pop(0)
            remaining -= 1
            walk.append(current)
        else:
            path = path_to(current, lambda node: bool(unused[node]))
            walk.extend(path[1:])
            current = path[-1]
    if current != start:
        walk.extend(path_to(current, lambda node: node == start)[1:])
    return tuple(node[0] for node in walk[:-1])


def _chain(residues, junctions, period, radius):
    """The periodic chain of ``period``, 2q + 1 copies long with H-capped
    ends, as ``declared_chain`` rows: (rows, residue occurrences [(index,
    name, key -> row id)], row id -> occurrence index, q)."""
    entry, exit_ = _ends(junctions)
    period = tuple(period)
    for a, b in zip(period, period[1:] + period[:1]):
        if (a, b) not in junctions:
            raise ValueError(f"the period step {a}->{b} is not in the sequence support")
    width = sum(_backbone_atoms(residues[name], entry[name], exit_[name]) for name in period)
    q = ceil((int(radius) + 2) / width)
    asm = dc._Assembly()
    placed = []
    for index in range(len(period) * (2 * q + 1)):
        name = period[index % len(period)]
        graph = residues[name].graph
        atoms = {a.key: a for a in graph.atoms}
        keys = tuple(sorted(atoms))
        handles = asm.add_period(atoms, keys, {k: i for i, k in enumerate(keys)},
                                 dc._bonds_of(graph), index, set(keys))
        if placed:
            left, right, order = junctions[(placed[-1][1], name)]
            asm.bond(placed[-1][2][left], handles[right], order)
        placed.append((index, name, handles))
    _left, _right, order = junctions[(period[-1], period[0])]
    for _ in range(int(order)):
        asm.add_h(placed[0][2][entry[period[0]]], slot=4)
        asm.add_h(placed[-1][2][exit_[period[-1]]], slot=4)
    rows, ids = asm.rows()
    occurrences = [(index, name, {key: ids[h] for key, h in handles.items()})
                   for index, name, handles in placed]
    occurrence_of = {row: index for index, _name, table in occurrences for row in table.values()}
    return rows, occurrences, occurrence_of, q


def _ball(rows, heavy, root, reach):
    """Heavy-atom distances from ``root`` up to ``reach``, grown across
    multiple bonds, and the rows of that ball (its hydrogens kept, every
    other cut bond capped by a hydrogen)."""
    distance, queue = {root: 0}, deque([root])
    while queue:
        i = queue.popleft()
        if distance[i] == reach:
            continue
        for j, _order in rows[i][5]:
            if j in heavy and j not in distance:
                distance[j] = distance[i] + 1
                queue.append(j)
    grown = True
    while grown:
        grown = False
        for i in list(distance):
            for j, order in rows[i][5]:
                if j in heavy and j not in distance and float(order) != 1.0:
                    distance[j] = distance[i] + 1
                    grown = True
    keep = set(distance) | {j for i in distance for j, _order in rows[i][5] if j not in heavy}
    ball, cap = {}, max(rows) + 1
    for i in keep:
        edges = []
        for j, order in rows[i][5]:
            if j in keep:
                edges.append((j, order))
            else:
                ball[cap] = ("H", None, 0, 0, 0, ((i, 1.0),))
                edges.append((cap, 1.0))
                cap += 1
        ball[i] = rows[i][:5] + (tuple(sorted(edges)),)
    return distance, ball


def _occurrences(residues: Mapping, junctions: Mapping, locality_policy: Mapping,
                 period: Sequence[str], whole_chain: bool = False) -> List[tuple]:
    """``(site key, word, position, residue, atom key, start)`` for every
    heavy atom of the middle copy of the periodic chain of ``period``:
    ``word`` the residues within the environment's reach, ``position`` the
    atom's residue in it, ``start`` the period index of the word's first
    residue (any integer; read modulo ``len(period)``)."""
    from rmgpy.kmc.rmg_moves import molecule_from_rows
    from rmgpy.kmc.site_class import extract_site_class

    radius = int(locality_policy["radius"])
    period = tuple(period)
    rows, occurrences, occurrence_of, q = _chain(residues, junctions, period, radius)
    heavy = {i for i, row in rows.items() if row[0] != "H"}
    policy = dc.site_class_policy(locality_policy, MAX_HEAVY_ATOMS)
    whole = molecule_from_rows(rows) if whole_chain else None
    whole_atoms = {atom.id: atom for atom in whole.atoms} if whole_chain else None
    out = []
    for index, name, table in occurrences[q * len(period):(q + 1) * len(period)]:
        for key, root in sorted(table.items()):
            if root not in heavy:
                continue
            distance, ball = _ball(rows, heavy, root, radius + BALL_EXTRA)
            reach = sorted({occurrence_of[j] for j, d in distance.items() if d <= radius + 1})
            word = tuple(period[i % len(period)] for i in reach)
            if whole_chain:
                site, _mapping, _pattern = extract_site_class(whole, [("site", whole_atoms[root])], policy)
            else:
                molecule = molecule_from_rows(ball)
                atom = next(a for a in molecule.atoms if a.id == root)
                site, _mapping, _pattern = extract_site_class(molecule, [("site", atom)], policy)
            out.append((site, word, reach.index(index), name, key, reach[0]))
    return out


def period_environments(residues: Mapping, junctions: Mapping, locality_policy: Mapping,
                        period: Sequence[str], whole_chain: bool = False) -> Dict[str, List[tuple]]:
    """Site-class key -> sorted occurrences ``(word, residue, atom key)`` over
    the heavy atoms of the middle copy of the periodic chain of ``period``.

    ``word`` lists, in chain order, the residues whose atoms lie within the
    environment's reach (radius + 1 bonds). ``whole_chain`` extracts every
    key on the whole chain instead of a ball."""
    found: Dict[str, set] = {}
    for site, word, _position, name, key, _start in _occurrences(
            residues, junctions, locality_policy, period, whole_chain):
        found.setdefault(site, set()).add((word, name, key))
    return {site: sorted(items) for site, items in sorted(found.items())}


def sequence_environments(spec) -> Dict[str, List[tuple]]:
    """Every site-class environment a chain of the ``PolymerSpec`` can have:
    ``period_environments`` over the support's covering period at
    ``spec_word_order``."""
    period = covering_period(dict(spec.sequence_support), spec_word_order(spec))
    return period_environments(dict(spec.residues), spec.junctions, spec.locality_policy, period)


def environment_witnesses(spec) -> Tuple[Tuple[str, ...], Dict[str, List[tuple]]]:
    """``(period, witnesses)``: the covering period of ``sequence_environments``
    and, per site-class key, its occurrences ``(word, position, atom key,
    start)`` sorted shortest word first (``position``: the atom's residue in
    the word; ``start``: the period index of the word's first residue, read
    modulo the period, for the context on either side)."""
    period = covering_period(dict(spec.sequence_support), spec_word_order(spec))
    found: Dict[str, set] = {}
    for site, word, position, _name, key, start in _occurrences(
            dict(spec.residues), spec.junctions, spec.locality_policy, period):
        found.setdefault(site, set()).add((word, position, key, start % len(period)))
    witnesses = {site: sorted(items, key=lambda item: (len(item[0]),) + item)
                 for site, items in sorted(found.items())}
    return period, witnesses


def word_frequencies(spec, k: int) -> Dict[Tuple[str, ...], Fraction]:
    """Exact frequency of every k-residue word of a long chain drawn from
    the declaration's composition and sequence statistics.

    The frequency of a word w is π(w1) · Π T(w_i, w_{i+1}), with π the
    composition and T the transition matrix. An alternating cycle instead
    counts each of its rotations once. Only positive frequencies are kept."""
    composition = dict(spec.initial_ensemble.composition)
    statistics = spec.initial_ensemble.statistics
    if statistics.model == "alternating":
        cycle, out = statistics.cycle, {}
        for start in range(len(cycle)):
            word = tuple(cycle[(start + i) % len(cycle)] for i in range(k))
            out[word] = out.get(word, Fraction(0)) + Fraction(1, len(cycle))
        return out
    present = [r for r, x in sorted(composition.items()) if x > 0]
    if statistics.model == "ideal_random":
        transition = {i: {j: composition[j] for j in present} for i in present}
    else:
        transition = {i: dict(row) for i, row in statistics.transition}
    words = {(i,): composition[i] for i in present}
    for _ in range(int(k) - 1):
        words = {word + (j,): p * t for word, p in words.items()
                 for j, t in transition[word[-1]].items() if t}
    return words
