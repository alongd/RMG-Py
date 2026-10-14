"""Independent transition oracle for windowMatcherTest.py (not collected).

It enumerates elementary moves per family straight from chemistry rules over
the raw ``MeltState`` graph and builds complete mapped product states.  It
shares NO code with ``rmgpy.kmc.window_matcher``: no template matching, no
outcome key, no resonance enumeration.  Resonance equivalence is decided by a
brute-force pi-placement enumeration over integer (Kekule) bond orders, using
networkx for graph bookkeeping; the matcher instead enumerates by backtracking
and writes aromatic sextets as 1.5.

A product state is a mapping ``atom -> (element, isotope or 0, charge, radical,
lone_pairs, ((neighbour, order), ...))`` over every atom of the touched
reactant components.  ``canonical(product)`` returns a hashable key that is
equal for two products iff they are the same labelled sigma graph with the
same charges and lone pairs and resonance-equivalent pi/radical placement
(same placement set per conjugated region).

Families (scope is stated in chemistry terms and is the same for every state):

- ``h_abstraction_benzylic``: a benzylic carbon radical X (radical carbon not
  in a ring, bonded to a carbon of a six-membered carbon ring) takes an H
  bonded to an sp3 carbon of ANOTHER molecule.  Ordered (X, H).
- ``beta_scission_midchain``: a mid-chain carbon radical X (not in a ring,
  at least two non-ring carbon neighbours) breaks the C(alpha)-C(beta) bond
  where alpha is a non-ring sp3 carbon neighbour of X and beta a non-ring
  sp3 carbon neighbour of alpha (not X); X=alpha forms, beta becomes the
  radical.
- ``recombination_chain_end``: two chain-end carbon radicals (not in a ring,
  exactly one non-ring heavy neighbour) on different molecules bond.  Unordered.
- ``disproportionation``: a benzylic radical X takes an H from an sp3 non-ring
  carbon Y bonded to a carbon radical Z (Z not in a ring) of another molecule;
  Z=Y forms.  Ordered (acceptor, donor).
- ``h_shift``: intramolecular 1,n H shift: a benzylic radical X takes an H
  bonded to an sp3 carbon Y of the same molecule whose shortest heavy-atom
  path to X has n-1 bonds.
- ``h_atom_recombination``: two free H atoms (no bonds, one radical) bond.
"""

import itertools
from collections import Counter

import networkx as nx

VALENCE = {"H": 1, "C": 4, "N": 5, "O": 6, "F": 7, "Cl": 7, "Br": 7, "I": 7}
MAX_REGION_BONDS = 16

FAMILIES = ("h_abstraction_benzylic", "beta_scission_midchain", "recombination_chain_end",
            "disproportionation", "h_shift", "h_atom_recombination")


class OracleError(AssertionError):
    """The oracle met an input it cannot judge (never a silent pass)."""


# --------------------------------------------------------------------------
# brute-force resonance (integer pi placements per conjugated region)
# --------------------------------------------------------------------------
_REGION_CACHE = {}


def _free(rec):
    element, _iso, charge, _rad, lone_pairs, edges = rec
    return VALENCE[element] - charge - 2 * lone_pairs - len(edges)


def _regions(product):
    free = {a: _free(rec) for a, rec in product.items()}
    for a, value in free.items():
        if value < 0:
            raise OracleError(f"atom {a} has negative free valence")
    graph = nx.Graph()
    pi_atoms = [a for a in product if free[a] > 0]
    graph.add_nodes_from(pi_atoms)
    for a in pi_atoms:
        for b, _order in product[a][5]:
            if b in free and free[b] > 0:
                graph.add_edge(a, b)
    return [tuple(sorted(c)) for c in nx.connected_components(graph)], free


def _placements(region, product, free):
    """All integer pi placements of a region with its radical count fixed:
    (bonds, [(pi per bond, radical per atom)])."""
    members = set(region)
    bonds = tuple(sorted((a, b) for a in region for b, _o in product[a][5] if b > a and b in members))
    radicals = sum(product[a][3] for a in region)
    key = (region, bonds, tuple(free[a] for a in region), radicals)
    if key in _REGION_CACHE:
        return bonds, _REGION_CACHE[key]
    if len(bonds) > MAX_REGION_BONDS:
        raise OracleError(f"region of {len(bonds)} bonds is too large for brute force")
    ranges = [range(min(free[a], free[b], 2) + 1) for a, b in bonds]
    found = []
    for combo in itertools.product(*ranges):
        used = Counter()
        for (a, b), p in zip(bonds, combo):
            used[a] += p
            used[b] += p
        if any(used[a] > free[a] for a in region):
            continue
        rad = tuple(free[a] - used[a] for a in region)
        if sum(rad) != radicals:
            continue
        found.append((combo, rad))
    found.sort()
    if len(_REGION_CACHE) > 200000:
        _REGION_CACHE.clear()
    _REGION_CACHE[key] = found
    return bonds, found


def _consistent(region, bonds, placement, product):
    """Is the given (possibly aromatic-notation) form a rendering of this
    integer placement?  1.5 bonds must form alternating six-cycles."""
    combo, rad = placement
    for a, r in zip(region, rad):
        if product[a][3] != r:  # the radical field of the record
            return False
    arom = []
    for (a, b), p in zip(bonds, combo):
        order = dict(product[a][5])[b]
        if order == 1.5:
            if p > 1:
                return False
            arom.append((a, b, p))
        elif order != 1 + p:
            return False
    if not arom:
        return True
    per_atom = Counter()
    pi_on = Counter()
    graph = nx.Graph()
    for a, b, p in arom:
        per_atom[a] += 1
        per_atom[b] += 1
        pi_on[a] += p
        pi_on[b] += p
        graph.add_edge(a, b)
    if any(per_atom[a] != 2 or pi_on[a] != 1 for a in per_atom):
        return False
    return all(len(c) == 6 for c in nx.connected_components(graph))


_CONSISTENT_CACHE = {}


def _region_matches(region, bonds, placements, product):
    """Is the form of this region (as written in ``product``) a rendering of
    one of its integer placements?  Cached on the region's written values."""
    key = (region, tuple(product[a][3] for a in region),
           tuple(dict(product[a][5])[b] for a, b in bonds))
    if key not in _CONSISTENT_CACHE:
        if len(_CONSISTENT_CACHE) > 200000:
            _CONSISTENT_CACHE.clear()
        _CONSISTENT_CACHE[key] = any(_consistent(region, bonds, pl, product) for pl in placements)
    return _CONSISTENT_CACHE[key]


def canonical(product):
    """Resonance-invariant key of a complete mapped product state."""
    for a, rec in product.items():
        element, _iso, charge, rad, lp, edges = rec
        total = sum(o for _, o in edges) + rad + 2 * lp
        if total != VALENCE[element] - charge:
            raise OracleError(f"atom {a} violates valence in the product")
        for b, o in edges:
            if b not in product or dict(product[b][5]).get(a) != o:
                raise OracleError(f"asymmetric or dangling bond {a}-{b}")
    regions, free = _regions(product)
    pi_order = {}
    radical = {a: product[a][3] for a in product}
    for region in regions:
        bonds, placements = _placements(region, product, free)
        if not _region_matches(region, bonds, placements, product):
            raise OracleError(f"form of region {region[:4]}... is not a valid resonance form")
        combo, rad = placements[0]
        for (a, b), p in zip(bonds, combo):
            pi_order[(a, b)] = 1.0 + p
        for a, r in zip(region, rad):
            radical[a] = r
    rows = []
    for a in sorted(product):
        element, iso, charge, _rad, lp, edges = product[a]
        canon_edges = tuple((b, pi_order[(min(a, b), max(a, b))]) if (min(a, b), max(a, b)) in pi_order
                            else (b, float(o)) for b, o in sorted(edges))
        rows.append((a, element, iso or 0, charge, lp, radical[a], canon_edges))
    return tuple(rows)


# --------------------------------------------------------------------------
# the oracle over one state
# --------------------------------------------------------------------------
class Oracle:
    """Elementary moves of a MeltState, enumerated from chemistry rules."""

    def __init__(self, state):
        self.records = {}
        graph = nx.Graph()
        for i, node in state.atoms.items():
            self.records[i] = (node.element, node.isotope or 0, node.charge, node.radical,
                               node.lone_pairs, tuple(sorted((j, float(o)) for j, o in node.edges.items())))
            graph.add_node(i)
        for i, node in state.atoms.items():
            for j in node.edges:
                graph.add_edge(i, j)
        self.graph = graph
        self.molecule = {}
        for comp in nx.connected_components(graph):
            members = tuple(sorted(comp))
            for a in members:
                self.molecule[a] = members
        self.in_ring = set()
        self.six_ring_carbon = set()
        for cycle in nx.cycle_basis(graph):
            self.in_ring.update(cycle)
            if len(cycle) == 6 and all(self.element(a) == "C" for a in cycle):
                self.six_ring_carbon.update(cycle)
        self.heavy = graph.subgraph([a for a in graph if self.element(a) != "H"])
        self._sites = None

    # -- atom predicates -----------------------------------------------------
    def element(self, a):
        return self.records[a][0]

    def neighbours(self, a):
        return [b for b, _ in self.records[a][5]]

    def sp3_carbon(self, a):
        return self.element(a) == "C" and len(self.records[a][5]) == 4

    def radical_sites(self):
        """Atoms carrying an unpaired electron in some resonance placement."""
        if self._sites is None:
            sites = set()
            done = set()
            for a in sorted(self.records):
                mol = self.molecule[a]
                if mol in done:
                    continue
                done.add(mol)
                product = {x: self.records[x] for x in mol}
                regions, free = _regions(product)
                for region in regions:
                    _bonds, placements = _placements(region, product, free)
                    for _combo, rad in placements:
                        sites.update(x for x, r in zip(region, rad) if r > 0)
            self._sites = sorted(sites)
        return self._sites

    def benzylic(self, x):
        return (self.element(x) == "C" and x not in self.in_ring
                and any(b in self.six_ring_carbon for b in self.neighbours(x)))

    def non_ring_heavy_neighbours(self, x):
        return [b for b in self.neighbours(x) if self.element(b) != "H" and b not in self.in_ring]

    # -- product construction -------------------------------------------------
    def _product(self, *atoms):
        mols = sorted({self.molecule[a] for a in atoms})
        out = {}
        for mol in mols:
            for x in mol:
                el, iso, chg, rad, lp, edges = self.records[x]
                out[x] = [el, iso, chg, rad, lp, dict(edges)]
        return out

    @staticmethod
    def _localize(product, x):
        """Rewrite x's conjugated region (integer orders) so x carries a radical."""
        frozen = {a: (r[0], r[1], r[2], r[3], r[4], tuple(sorted(r[5].items()))) for a, r in product.items()}
        regions, free = _regions(frozen)
        for region in regions:
            if x not in region:
                continue
            bonds, placements = _placements(region, frozen, free)
            for combo, rad in placements:
                if rad[region.index(x)] > 0:
                    for (a, b), p in zip(bonds, combo):
                        product[a][5][b] = product[b][5][a] = 1.0 + p
                    for a, r in zip(region, rad):
                        product[a][3] = r
                    return
        raise OracleError(f"atom {x} carries no radical in any placement")

    @staticmethod
    def _bond(product, a, b, order):
        if order == 0:
            del product[a][5][b]
            del product[b][5][a]
        else:
            product[a][5][b] = product[b][5][a] = float(order)

    @staticmethod
    def _freeze(product):
        return {a: (r[0], r[1], r[2], r[3], r[4], tuple(sorted(r[5].items()))) for a, r in product.items()}

    def _key(self, product):
        return canonical(self._freeze(product))

    # -- families ------------------------------------------------------------
    def moves(self, family, n=None):
        """dict pathway key (family, product key) -> description, plus the raw
        move count as ``self.last_raw``."""
        generator = getattr(self, "_" + family)
        out = {}
        raw = 0
        for description, product in (generator(n) if family == "h_shift" else generator()):
            raw += 1
            key = (family, self._key(product))
            out.setdefault(key, description)
        self.last_raw = raw
        return out

    def _h_abstraction_benzylic(self):
        for x in self.radical_sites():
            if not self.benzylic(x):
                continue
            for y in sorted(self.records):
                if self.molecule[y] == self.molecule[x] or not self.sp3_carbon(y):
                    continue
                for h in self.neighbours(y):
                    if self.element(h) != "H":
                        continue
                    p = self._product(x, y)
                    self._localize(p, x)
                    p[x][3] -= 1
                    self._bond(p, x, h, 1)
                    self._bond(p, y, h, 0)
                    p[y][3] += 1
                    yield ("HA", x, y, h), p

    def _beta_scission_midchain(self):
        for x in self.radical_sites():
            if self.element(x) != "C" or x in self.in_ring:
                continue
            carbons = [b for b in self.non_ring_heavy_neighbours(x) if self.element(b) == "C"]
            if len(carbons) < 2:
                continue
            for alpha in carbons:
                if not self.sp3_carbon(alpha):
                    continue
                for beta in self.neighbours(alpha):
                    if beta == x or not self.sp3_carbon(beta) or beta in self.in_ring:
                        continue
                    p = self._product(x)
                    self._localize(p, x)
                    if p[x][5][alpha] != 1.0:
                        raise OracleError("alpha bond is not single")
                    p[x][3] -= 1
                    self._bond(p, x, alpha, 2)
                    self._bond(p, alpha, beta, 0)
                    p[beta][3] += 1
                    yield ("BS", x, alpha, beta), p

    def _chain_end(self, x):
        return (self.element(x) == "C" and x not in self.in_ring
                and len(self.non_ring_heavy_neighbours(x)) == 1)

    def _recombination_chain_end(self):
        ends = [x for x in self.radical_sites() if self._chain_end(x)]
        for i, x in enumerate(ends):
            for y in ends[i + 1:]:
                if self.molecule[x] == self.molecule[y]:
                    continue
                p = self._product(x, y)
                self._localize(p, x)
                self._localize(p, y)
                p[x][3] -= 1
                p[y][3] -= 1
                self._bond(p, x, y, 1)
                yield ("RC", x, y), p

    def _disproportionation(self):
        sites = self.radical_sites()
        for x in sites:
            if not self.benzylic(x):
                continue
            for z in sites:
                if self.element(z) != "C" or z in self.in_ring or self.molecule[z] == self.molecule[x]:
                    continue
                for y in self.neighbours(z):
                    if not self.sp3_carbon(y) or y in self.in_ring:
                        continue
                    for h in self.neighbours(y):
                        if self.element(h) != "H":
                            continue
                        p = self._product(x, z)
                        self._localize(p, x)
                        self._localize(p, z)
                        p[x][3] -= 1
                        self._bond(p, x, h, 1)
                        self._bond(p, y, h, 0)
                        self._bond(p, z, y, p[z][5][y] + 1)
                        p[z][3] -= 1
                        yield ("DS", x, z, y, h), p

    def _h_shift(self, n):
        if n is None or n < 3:
            raise OracleError("an H shift needs n >= 3")
        for x in self.radical_sites():
            if not self.benzylic(x):
                continue
            dist = nx.single_source_shortest_path_length(self.heavy, x, cutoff=n - 1)
            for y in sorted(dist):
                if dist[y] != n - 1 or not self.sp3_carbon(y):
                    continue
                for h in self.neighbours(y):
                    if self.element(h) != "H":
                        continue
                    p = self._product(x)
                    self._localize(p, x)
                    p[x][3] -= 1
                    self._bond(p, x, h, 1)
                    self._bond(p, y, h, 0)
                    p[y][3] += 1
                    yield ("HS", n, x, y, h), p

    def _h_atom_recombination(self):
        free = [a for a in sorted(self.records)
                if self.element(a) == "H" and not self.records[a][5] and self.records[a][3] == 1]
        for i, a in enumerate(free):
            for b in free[i + 1:]:
                p = self._product(a, b)
                p[a][3] = p[b][3] = 0
                self._bond(p, a, b, 1)
                yield ("HH", a, b), p
