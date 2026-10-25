"""Shared builders for the branched (schema v2, ``architecture: branched``)
tests: the branched novolac test fixture (``fixtures/decks/novolac_v2_deck.py``),
mutations of it, and an oracle for its environment set that reads the
declaration mapping alone (no lift, no word): every neighbourhood of every
residue with every site choice up to ``depth`` residue hops, vacant leaves
beyond, a forbidden site at a leaf continued by its least partner."""

import copy
import itertools
from pathlib import Path

from polymer_spec_v2_support import polymer_keyword

DECKS = Path(__file__).parent / "fixtures" / "decks"
BRANCHED_DECK = DECKS / "novolac_v2_deck.py"


def branched_mapping():
    """The branched fixture's declaration, as a fresh mutable mapping."""
    return copy.deepcopy(polymer_keyword(BRANCHED_DECK, "polymer_spec"))


def branched_spec(mutate=None):
    from rmgpy.polymer_spec import normalized_polymer_spec
    mapping = branched_mapping()
    if mutate is not None:
        mutate(mapping)
    return normalized_polymer_spec(mapping)


def counts(residues=None, links=None, seed=None):
    def mutate(mapping):
        ensemble = mapping["initial_ensemble"]
        if residues is not None:
            ensemble["representative_residue_counts"] = dict(residues)
        if links is not None:
            ensemble["link_counts"] = dict(links)
        if seed is not None:
            ensemble["seed"] = seed
    return mutate


class Oracle:
    """The declaration's residues, sites and links read from the mapping."""

    def __init__(self, mapping):
        from rmgpy.kmc import declared_chain as dc
        from rmgpy.kmc import declared_environments as de
        from rmgpy.polymer_spec import GraphDeclaration

        self.graphs = {name: GraphDeclaration.from_mapping(residue["graph"])
                       for name, residue in mapping["residues"].items()}
        self.sites = {name: {s: (v["atom"], v["class"], v["vacancy"]) for s, v in residue["sites"].items()}
                      for name, residue in mapping["residues"].items()}
        self.links = {(link["from"], link["to"]) for link in mapping["links"].values()}
        self.locality = mapping["locality_policy"]
        self.radius = int(self.locality["radius"])
        self.policy = dc.site_class_policy(self.locality, de.MAX_HEAVY_ATOMS)

    def linkable(self, c1, c2):
        return (c1, c2) in self.links or (c2, c1) in self.links

    def partners(self, cls):
        return [(r, s) for r in sorted(self.sites) for s, (_a, c, _v) in sorted(self.sites[r].items())
                if self.linkable(cls, c)]

    def build(self, nodes, links):
        """``nodes`` residue names, ``links`` (node, site, node, site): rows
        and (node, atom key) -> row id, every vacant site H-capped."""
        from rmgpy.kmc import declared_chain as dc
        asm = dc._Assembly()
        handles, occupied = [], set()
        for i, name in enumerate(nodes):
            graph = self.graphs[name]
            atoms = {a.key: a for a in graph.atoms}
            keys = tuple(sorted(atoms))
            handles.append(asm.add_period(atoms, keys, {k: j for j, k in enumerate(keys)},
                                          dc._bonds_of(graph), i, set(keys)))
        for i, s, j, t in links:
            asm.bond(handles[i][self.sites[nodes[i]][s][0]], handles[j][self.sites[nodes[j]][t][0]], 1.0)
            occupied |= {(i, s), (j, t)}
        for i, name in enumerate(nodes):
            for s, (atom, _c, vacancy) in sorted(self.sites[name].items()):
                if (i, s) not in occupied:
                    assert vacancy == "H", f"forbidden site {name}.{s} left vacant"
                    asm.add_h(handles[i][atom], 4)
        asm.check_valence("oracle")
        rows, ids = asm.rows()
        return rows, {(i, k): ids[h] for i, table in enumerate(handles) for k, h in table.items()}

    def keys(self, rows, roots):
        """site-class key of every ``(label, row id)`` root (the environment ball)."""
        from rmgpy.kmc import declared_environments as de
        from rmgpy.kmc.rmg_moves import molecule_from_rows
        from rmgpy.kmc.site_class import extract_site_class
        heavy = {i for i, row in rows.items() if row[0] != "H"}
        out = {}
        for label, root in roots:
            _distance, ball = de._ball(rows, heavy, root, self.radius + de.BALL_EXTRA)
            molecule = molecule_from_rows(ball)
            atom = next(a for a in molecule.atoms if a.id == root)
            out[label] = extract_site_class(molecule, [("site", atom)], self.policy)[0]
        return out

    def _leaf(self, name, entered):
        nodes, links = [name], []
        for s in sorted(self.sites[name]):
            _atom, cls, vacancy = self.sites[name][s]
            if s == entered or vacancy == "H":
                continue
            partner, site = self.partners(cls)[0]
            sub_nodes, sub_links = self._leaf(partner, site)
            base = len(nodes)
            nodes += sub_nodes
            links += [(i + base, a, j + base, b) for i, a, j, b in sub_links]
            links.append((0, s, base, site))
        return nodes, links

    def _trees(self, name, entered, depth):
        sites = [s for s in sorted(self.sites[name]) if s != entered]
        choices = []
        for s in sites:
            _atom, cls, vacancy = self.sites[name][s]
            if depth <= 0:
                choices.append([None] if vacancy == "H" else [(t, self._leaf(p, t))
                                                              for p, t in self.partners(cls)[:1]])
                continue
            options = [None] if vacancy == "H" else []
            for partner, site in self.partners(cls):
                options += [(site, tree) for tree in self._trees(partner, site, depth - 1)]
            choices.append(options)
        for combo in itertools.product(*choices):
            nodes, links = [name], []
            for s, option in zip(sites, combo):
                if option is None:
                    continue
                site, (sub_nodes, sub_links) = option
                base = len(nodes)
                nodes += sub_nodes
                links += [(i + base, a, j + base, b) for i, a, j, b in sub_links]
                links.append((0, s, base, site))
            yield nodes, links

    def environments(self, depth):
        """site-class key -> {(residue, atom key)} over every neighbourhood."""
        found = {}
        for centre in sorted(self.sites):
            keys = sorted(a.key for a in self.graphs[centre].atoms if a.element != "H")
            for nodes, links in self._trees(centre, None, depth):
                rows, table = self.build(nodes, links)
                for key, site in self.keys(rows, [(k, table[(0, k)]) for k in keys]).items():
                    found.setdefault(site, set()).add((centre, key))
        return found


def state_keys(state, oracle):
    """site-class key -> {(residue, atom key)} of every heavy residue atom of a
    melt state (``Origin.segment`` ``unit <k> <residue>``)."""
    from rmgpy.kmc.melt_state import segment_residue
    rows = {a: (n.element, n.isotope, n.charge, n.radical, n.lone_pairs, tuple(sorted(n.edges.items())))
            for a, n in state.atoms.items()}
    roots = [((segment_residue(o.segment), o.key), a) for a, o in state.origin.items()
             if state.atoms[a].element != "H" and segment_residue(o.segment)]
    found = {}
    for (label, a), key in zip(roots, oracle.keys(rows, [(i, a) for i, (_l, a) in enumerate(roots)]).values()):
        found.setdefault(key, set()).add(label)
    return found
