"""The environment set of a network declaration (schema v2,
``architecture: network``; ``rmgpy.polymer_spec.NetworkPolymerSpec``).

**Trees.**  Above the gel point almost every residue is linked on several
sites, so an environment is read on a neighbourhood of the network, not on a
chain word.  ``minimum_cycle_bonds`` >= 2(r + 1) + 2 keeps every
neighbourhood an environment reads a tree (a shorter cycle would fit inside
one environment, and the generator and the runtime both refuse one), so the
environment set is the set of rooted residue trees over the occupied, vacant
and branch alphabet:

- the root is one residue; each of its sites is vacant (``vacancy: H``, an H
  cap) or linked to a site of any residue whose class a declared link pairs
  with it;
- a linked residue's other sites are expanded the same way while the site's
  atom lies within ``site_reach`` = r + 1 bonds of the root residue (an
  environment reads r bonds whole and the heavy bond orders of the atoms
  one bond further, so a site at r + 1 bonds changes a key by being occupied);
- beyond the reach a site is vacant when it may be, and a forbidden-vacant
  one is continued by its least partner (residue, then site name) with every
  other site vacant or so continued: what lies there is out of every reach.

Every heavy atom of the root residue is an environment root; its key is the
``declared_environments`` ball key (declared bond orders, the size-independent
policy), so keys compare with the S6 lift's (``declared_branches``).

**The compile reuses the lift.**  The lift's words (chains of residue
variants with stub arms) hold every neighbourhood of the declaration
(``network_environment_gap`` is empty for a declaration whose lift covers its
trees), so a network compiles as its lift; the compile refuses a network with
a gap by name (``refused:network_environment``).  Nothing here names a
polymer."""

from __future__ import annotations

import hashlib
import itertools
from typing import Dict, Optional, Tuple

from rmgpy.kmc import declared_chain as dc
from rmgpy.kmc import declared_environments as de

GAP_REFUSAL = "refused:network_environment"


def site_reach(spec) -> int:
    """Bonds from the root residue within which a site is expanded: r + 1."""
    return int(spec.locality_policy["radius"]) + 1


def _site_distances(graph) -> Dict[Tuple[str, str], int]:
    heavy = {a.key for a in graph.atoms if a.element != "H"}
    adjacency = {k: set() for k in heavy}
    for bond in graph.bonds:
        p, q = bond.atoms
        if p in heavy and q in heavy:
            adjacency[p].add(q)
            adjacency[q].add(p)
    out = {}
    for a in heavy:
        distance, frontier = {a: 0}, [a]
        while frontier:
            following = []
            for x in frontier:
                for y in adjacency[x]:
                    if y not in distance:
                        distance[y] = distance[x] + 1
                        following.append(y)
            frontier = following
        for b, d in distance.items():
            out[(a, b)] = d
    return out


class _Trees:
    def __init__(self, spec, reach):
        self.spec = spec
        self.reach = int(reach)
        self.graphs = {name: residue.graph for name, residue in spec.residues}
        self.sites = {name: dict(residue.sites) for name, residue in spec.residues}
        pairs = set()
        for _n, link in spec.links:
            pairs |= {(link.source, link.target), (link.target, link.source)}
        self.pairs = pairs
        self.distances = {name: _site_distances(graph) for name, graph in self.graphs.items()}
        self.radius = int(spec.locality_policy["radius"])
        self.policy = dc.site_class_policy(spec.locality_policy, de.MAX_HEAVY_ATOMS)

    def partners(self, cls):
        return [(r, s) for r in sorted(self.sites) for s, site in sorted(self.sites[r].items())
                if (cls, site.site_class) in self.pairs]

    def _leaf(self, name, entered):
        nodes, links = [name], []
        for s, site in sorted(self.sites[name].items()):
            if s == entered or site.vacancy == "H":
                continue
            partner, psite = self.partners(site.site_class)[0]
            sub_nodes, sub_links = self._leaf(partner, psite)
            base = len(nodes)
            nodes += sub_nodes
            links += [(i + base, a, j + base, b) for i, a, j, b in sub_links]
            links.append((0, s, base, psite))
        return nodes, links

    def trees(self, name, entered=None, at=0):
        """(nodes, links) of every tree rooted at ``name`` entered through
        site ``entered`` whose atom is ``at`` bonds from the root residue."""
        sites = [s for s in sorted(self.sites[name]) if s != entered]
        entry_atom = None if entered is None else self.sites[name][entered].atom
        choices = []
        for s in sites:
            site = self.sites[name][s]
            d = at if entry_atom is None else at + self.distances[name][(entry_atom, site.atom)]
            if d > self.reach:
                if site.vacancy == "H":
                    choices.append([None])
                else:
                    partner, psite = self.partners(site.site_class)[0]
                    choices.append([(psite, self._leaf(partner, psite))])
                continue
            options = [None] if site.vacancy == "H" else []
            for partner, psite in self.partners(site.site_class):
                options += [(psite, tree) for tree in self.trees(partner, psite, d + 1)]
            choices.append(options)
        for combo in itertools.product(*choices):
            nodes, links = [name], []
            for s, option in zip(sites, combo):
                if option is None:
                    continue
                psite, (sub_nodes, sub_links) = option
                base = len(nodes)
                nodes += sub_nodes
                links += [(i + base, a, j + base, b) for i, a, j, b in sub_links]
                links.append((0, s, base, psite))
            yield nodes, links

    def rows(self, nodes, links):
        asm = dc._Assembly()
        handles, occupied = [], set()
        for i, name in enumerate(nodes):
            graph = self.graphs[name]
            atoms = {a.key: a for a in graph.atoms}
            keys = tuple(sorted(atoms))
            handles.append(asm.add_period(atoms, keys, {k: j for j, k in enumerate(keys)},
                                          dc._bonds_of(graph), i, set(keys)))
        for i, s, j, t in links:
            si, sj = self.sites[nodes[i]][s], self.sites[nodes[j]][t]
            asm.bond(handles[i][si.atom], handles[j][sj.atom], si.bond_order)
            occupied |= {(i, s), (j, t)}
        for i, name in enumerate(nodes):
            for s, site in sorted(self.sites[name].items()):
                if (i, s) not in occupied:
                    asm.add_h(handles[i][site.atom], 4)
        rows, ids = asm.rows()
        return rows, {(i, k): ids[h] for i, table in enumerate(handles) for k, h in table.items()}

    def environments(self):
        from rmgpy.kmc.rmg_moves import molecule_from_rows
        from rmgpy.kmc.site_class import extract_site_class
        found = {}
        for centre in sorted(self.sites):
            keys = sorted(a.key for a in self.graphs[centre].atoms if a.element != "H")
            for nodes, links in self.trees(centre):
                rows, table = self.rows(nodes, links)
                heavy = {i for i, row in rows.items() if row[0] != "H"}
                for key in keys:
                    root = table[(0, key)]
                    _distance, ball = de._ball(rows, heavy, root, self.radius + de.BALL_EXTRA)
                    molecule = molecule_from_rows(ball)
                    atom = next(a for a in molecule.atoms if a.id == root)
                    site = extract_site_class(molecule, [("site", atom)], self.policy)[0]
                    found.setdefault(site, set()).add((centre, key))
        return {site: found[site] for site in sorted(found)}


_TREES: Dict[tuple, Dict[str, set]] = {}


def tree_environments(spec, reach: Optional[int] = None) -> Dict[str, set]:
    """site-class key -> {(residue, atom key)} over every rooted tree of the
    declaration (module doc); ``reach`` defaults to ``site_reach``."""
    from rmgpy.polymer_spec import BranchedPolymerSpec, normalized_polymer_spec
    spec = normalized_polymer_spec(spec)
    if not isinstance(spec, BranchedPolymerSpec):
        raise TypeError(f"tree environments take a branched or network declaration, got {type(spec).__name__}")
    reach = site_reach(spec) if reach is None else int(reach)
    key = (hashlib.sha256(spec.normalized_json).hexdigest(), reach)
    if key not in _TREES:
        _TREES[key] = _Trees(spec, reach).environments()
    return {site: set(found) for site, found in _TREES[key].items()}


def network_environment_gap(spec) -> Tuple[str, ...]:
    """The tree environments the S6 lift's environment set lacks (sorted)."""
    from rmgpy.kmc.declared_branches import branched_environments
    lift = branched_environments(spec)
    return tuple(sorted(set(tree_environments(spec)) - set(lift)))
