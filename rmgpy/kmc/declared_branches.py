"""The environment set and compile chains of a branched declaration (schema
v2, ``architecture: branched``; ``rmgpy.polymer_spec.BranchedPolymerSpec``).

**The lift.**  A branched declaration is lifted to a linear alphabet of
residue *variants* (``branched_lift``), so the linear-v2 machinery (the
covering period, the sequence chain provider, the compile structures and
the roster) runs on it unchanged.  A variant is one residue with
- one *entry* and one *exit* site: the two sites the chain runs through
  (occupied in chain);
- every other site either *vacant* (``H``: an H cap) or a *branch*: the site
  linked to a stub arm.  A forbidden-vacant site is always a branch.

A **stub arm** is a word of residues hanging off the branch site: the first
``stub_depth`` = ceil((r + 1) / n_min) residues with every entry site the
links allow enumerated (each continuing through its least other site, or its
forbidden one), then canonical residues (the least partner) until one can end
the arm (every other site vacant-capable).  n_min is the fewest backbone
atoms (entry -> exit path) of any residue.  Stub atoms carry keys
``~<site>.<i>.<key>`` (``ARM_PREFIX``); they are never environment roots and
never take a radical edit (``refused:stub_site``): their environments are the
other words'.  Variants are deduplicated by isomorphism of (graph, entry
atom, exit atom); the variant name is the residue tree (``variant_name``,
``parse_variant``), e.g. ``A(x>y;z=[B(u>v)/A(x>y)/B(u>v)/A(x)])``.

The lift's support joins an exit to an entry whose classes a declared link
joins, in either direction; its links are ``out:X -> in:Y`` per direction, so
``polymer_spec.linear_junctions`` validates every junction (rules 2 and 2b).

**Word order.**  An environment reads r + 1 bonds on each side of its atom,
so it touches at most rho residues beyond its own on one side, rho = 1 + the
most residues of one support path whose backbone atoms sum to at most r; a
centred window of k = 2 rho + 1 residues holds every environment (the S1
window argument made exact over the support).

**The environment set** (``BranchedLift.environments``) is the union of three
kinds of word, every key computed on the ball of
``declared_environments`` around a non-stub heavy atom, once per touched
residue word (the residues within r + 1 bonds, in chain order, with the real
ends it touches):
- ``interior``: every position of the covering period at k, rho + 1 residues
  of context on each side;
- ``end``: an H-capped terminal residue (its entry site vacant), positions 0
  to rho of every word of k + 1 residues inward;
- ``short``: every position of every word of 1 to k residues whose two end
  sites are vacant (an isolated residue, a dimer, ...).
A terminal residue has one occupied site and an in-chain one two, so the
interior words alone miss ends; molecules shorter than 2 rho + 1 residues
show their two ends to one atom.

**Compile** (``BranchedChainProvider``): the lift's sequence chain provider
under ``branched:<sha256[:16]>`` of the declaration.  Its depository growth
units are the residues' tree increments (a residue added anywhere in a tree
adds its formula with every site occupied plus (sites - 2) H), so the size
and formula exclusion arguments hold for any tree growth; it names no growth
sites, so a decision that needs the analogues built is ``refused:unknown``
(branch growth by construction is not built).  Nothing here names a polymer."""

from __future__ import annotations

import hashlib
import itertools
import re
from collections import deque
from dataclasses import dataclass
from math import ceil
from typing import Dict, Mapping, Optional, Sequence, Tuple

from rmgpy.kmc import declared_chain as dc
from rmgpy.kmc import declared_environments as de
from rmgpy.kmc.boundaries import PROVIDERS
from rmgpy.kmc.declared_sequence import SequenceChainProvider

ARM_PREFIX = "~"
STUB_REFUSAL = "refused:stub_site"
STUB_ARM_REFUSAL = "refused:stub_arm"
KINDS = ("interior", "end", "short")


@dataclass(frozen=True)
class Variant:
    """One residue of the lift: ``roles`` ((site, "H" | arm), ...) sorted by
    site for every site but the entry and the exit; an arm is ((residue,
    entry site, exit site or None), ...) outward."""

    residue: str
    entry: str
    exit: str
    roles: tuple


def variant_name(variant: Variant) -> str:
    def arm_text(arm):
        return "[" + "/".join(f"{r}({e}>{x})" if x is not None else f"{r}({e})" for r, e, x in arm) + "]"
    roles = ";".join(f"{site}=" + ("H" if role == "H" else arm_text(role)) for site, role in variant.roles)
    return f"{variant.residue}({variant.entry}>{variant.exit}" + (f";{roles}" if roles else "") + ")"


_STEP = re.compile(r"^([^()\[\]/;=>]+)\(([^()\[\]/;=>]+)(?:>([^()\[\]/;=>]+))?\)$")


def parse_variant(name: str) -> Variant:
    """The residue tree a variant name spells (``variant_name``'s inverse)."""
    head, _, rest = name.partition("(")
    if not rest.endswith(")"):
        raise ValueError(f"{name!r} is not a variant name")
    body = rest[:-1]
    parts, depth, current = [], 0, ""
    for ch in body:
        if ch == ";" and depth == 0:
            parts.append(current)
            current = ""
            continue
        depth += ch == "["
        depth -= ch == "]"
        current += ch
    parts.append(current)
    entry, _, exit_ = parts[0].partition(">")
    roles = []
    for part in parts[1:]:
        site, _, role = part.partition("=")
        if role == "H":
            roles.append((site, "H"))
            continue
        if not (role.startswith("[") and role.endswith("]")):
            raise ValueError(f"{name!r}: role {part!r} is neither H nor an arm")
        arm = []
        for step in role[1:-1].split("/"):
            match = _STEP.match(step)
            if match is None:
                raise ValueError(f"{name!r}: arm step {step!r}")
            arm.append((match.group(1), match.group(2), match.group(3)))
        roles.append((site, tuple(arm)))
    variant = Variant(head, entry, exit_, tuple(roles))
    if variant_name(variant) != name:
        raise ValueError(f"{name!r} is not a canonical variant name")
    return variant


@dataclass(frozen=True)
class Context:
    word: tuple
    position: int


def _lifted_atom(atom, key, extra_h):
    mapping = atom.as_mapping()
    mapping["key"] = key
    mapping["hydrogens"] = atom.hydrogens + extra_h
    return mapping


class BranchedLift:
    """The lift of one ``BranchedPolymerSpec`` (module doc).  Duck-types the
    parts of a linear ``PolymerSpec`` that the sequence chain provider, the
    environment code and the compile builder read: ``residues``,
    ``junctions``, ``sequence_support``, the end groups (the vacancy cap),
    ``monomer_products``, ``small_species`` and ``locality_policy``."""

    def __init__(self, spec, stub_depth: Optional[int] = None):
        from rmgpy.polymer_spec import LinkDeclaration, ResidueDeclaration, linear_junctions

        self.base = spec
        self.label = spec.label
        self.locality_policy = spec.locality_policy
        self.radius = int(spec.locality_policy["radius"])
        self.monomer_products = spec.monomer_products
        self.small_species = spec.small_species
        self.head_end_group = self.tail_end_group = spec.vacancy_cap
        self.residue_graphs = {name: residue.graph for name, residue in spec.residues}
        self.sites = {name: dict(residue.sites) for name, residue in spec.residues}
        self.link_classes = {(link.source, link.target) for _n, link in spec.links}
        for name, sites in self.sites.items():
            if len(sites) < 2:
                raise ValueError(f"{STUB_ARM_REFUSAL}: residue {name} has one site: it is never in a chain, "
                                 "so the lift has no word for it")
        self.n_min = min(self._path_atoms(name, e, x) for name in self.sites
                         for e, x in itertools.permutations(sorted(self.sites[name]), 2))
        self.stub_depth = ceil((self.radius + 1) / self.n_min) if stub_depth is None else int(stub_depth)
        variants = {}
        for name in sorted(self.sites):
            sites = sorted(self.sites[name])
            for entry, exit_ in itertools.permutations(sites, 2):
                others = [s for s in sites if s not in (entry, exit_)]
                options = []
                for s in others:
                    site = self.sites[name][s]
                    choices = ["H"] if site.vacancy == "H" else []
                    choices += self._arms(site.site_class)
                    options.append(choices)
                for combo in itertools.product(*options):
                    variant = Variant(name, entry, exit_, tuple(zip(others, combo)))
                    variants[variant_name(variant)] = variant
        graphs = {n: self._graph(v) for n, v in variants.items()}
        kept = self._distinct(graphs, variants)
        self.variants = {n: variants[n] for n in kept}
        self.residues = {}
        for n in kept:
            v = self.variants[n]
            entry, exit_ = self.sites[v.residue][v.entry], self.sites[v.residue][v.exit]
            atoms, bonds = graphs[n]
            self.residues[n] = ResidueDeclaration.from_mapping({
                "graph": {"atoms": atoms, "bonds": bonds, "electronic_state": {"multiplicity": 1}},
                "sites": {"entry": {"atom": entry.atom, "bond_order": 1, "class": "in:" + entry.site_class,
                                    "vacancy": "forbidden"},
                          "exit": {"atom": exit_.atom, "bond_order": 1, "class": "out:" + exit_.site_class,
                                   "vacancy": "forbidden"}}}, n)
        links = {}
        for x, y in sorted(self.link_classes):
            for a, b in ((x, y), (y, x)):
                links[f"{a}>{b}"] = LinkDeclaration("out:" + a, "in:" + b, 1.0)
        self.links = tuple(sorted(links.items()))
        support = {a: tuple(sorted(b for b in self.residues
                                   if self.linkable(self._class(a, "exit"), self._class(b, "entry"))))
                   for a in self.residues}
        self.sequence_support = tuple(sorted(support.items()))
        self.support = support
        self.junction_table = tuple(sorted(linear_junctions(self.residues, links, support).items()))
        self.entry, self.exit = de._ends(self.junctions)
        self.backbone = {n: de._backbone_atoms(r, self.entry[n], self.exit[n]) for n, r in self.residues.items()}
        self.rho = 1 + self._reach_residues()
        self.k = 2 * self.rho + 1
        self.policy = dc.site_class_policy(self.locality_policy, de.MAX_HEAVY_ATOMS)
        self._memo = {}
        self._environments = {}
        self.residues = dict(sorted(self.residues.items()))

    # -- declaration reads --------------------------------------------------------
    @property
    def junctions(self) -> Dict[Tuple[str, str], Tuple[str, str, float]]:
        return dict(self.junction_table)

    @property
    def normalized_json(self) -> bytes:
        return self.base.normalized_json

    def linkable(self, c1, c2) -> bool:
        return (c1, c2) in self.link_classes or (c2, c1) in self.link_classes

    def _class(self, name, role) -> str:
        v = self.variants[name]
        return self.sites[v.residue][v.entry if role == "entry" else v.exit].site_class

    def vacant_entry(self, name) -> bool:
        v = self.variants[name]
        return self.sites[v.residue][v.entry].vacancy == "H"

    def vacant_exit(self, name) -> bool:
        v = self.variants[name]
        return self.sites[v.residue][v.exit].vacancy == "H"

    def stub_keys(self, name) -> frozenset:
        return frozenset(a.key for a in self.residues[name].graph.atoms if a.key.startswith(ARM_PREFIX))

    def _path_atoms(self, name, e, x) -> int:
        graph = self.residue_graphs[name]
        a, b = self.sites[name][e].atom, self.sites[name][x].atom
        heavy = {atom.key for atom in graph.atoms if atom.element != "H"}
        adjacency = {k: set() for k in heavy}
        for bond in graph.bonds:
            p, q = bond.atoms
            if p in heavy and q in heavy:
                adjacency[p].add(q)
                adjacency[q].add(p)
        distance, queue = {a: 0}, deque([a])
        while queue:
            node = queue.popleft()
            for nxt in adjacency[node]:
                if nxt not in distance:
                    distance[nxt] = distance[node] + 1
                    queue.append(nxt)
        return distance[b] + 1

    def _partners(self, cls):
        return [(r, s) for r in sorted(self.sites) for s, site in sorted(self.sites[r].items())
                if self.linkable(cls, site.site_class)]

    def _arms(self, cls) -> list:
        out = []

        def grow(prefix, cls, i):
            options = self._partners(cls) if i < self.stub_depth else self._partners(cls)[:1]
            for residue, entry in options:
                others = sorted(s for s in self.sites[residue] if s != entry)
                forbidden = [s for s in others if self.sites[residue][s].vacancy != "H"]
                if len(forbidden) > 1:
                    raise ValueError(f"{STUB_ARM_REFUSAL}: residue {residue} entered at {entry} keeps "
                                     f"{len(forbidden)} forbidden-vacant sites; an arm continues through one")
                if i + 1 >= self.stub_depth and not forbidden:
                    out.append(tuple(prefix + [(residue, entry, None)]))
                    continue
                exit_ = forbidden[0] if forbidden else others[0]
                grow(prefix + [(residue, entry, exit_)], self.sites[residue][exit_].site_class, i + 1)
        grow([], cls, 0)
        return out

    def _graph(self, variant):
        """(atom mappings, bond mappings) of a variant: the residue with its
        vacant sites H-capped and its arms attached (module doc)."""
        graph = self.residue_graphs[variant.residue]
        extra = {}
        bonds = [b.as_mapping() for b in graph.bonds]
        arms = []
        for site, role in variant.roles:
            atom = self.sites[variant.residue][site].atom
            if role == "H":
                extra[atom] = extra.get(atom, 0) + 1
            else:
                arms.append((site, atom, role))
        atoms = {a.key: _lifted_atom(a, a.key, extra.get(a.key, 0)) for a in graph.atoms}
        for site, atom, arm in arms:
            previous = atom
            for i, (residue, entry, exit_) in enumerate(arm):
                prefix = f"{ARM_PREFIX}{site}.{i}."
                sub = self.residue_graphs[residue]
                caps = {}
                for s, item in self.sites[residue].items():
                    if s not in (entry, exit_):
                        caps[item.atom] = caps.get(item.atom, 0) + 1
                for a in sub.atoms:
                    atoms[prefix + a.key] = _lifted_atom(a, prefix + a.key, caps.get(a.key, 0))
                bonds += [{"atoms": [prefix + b.atoms[0], prefix + b.atoms[1]], "order": b.as_mapping()["order"]}
                          for b in sub.bonds]
                bonds.append({"atoms": [previous, prefix + self.sites[residue][entry].atom], "order": 1})
                if exit_ is not None:
                    previous = prefix + self.sites[residue][exit_].atom
        return [atoms[k] for k in sorted(atoms)], bonds

    def _distinct(self, graphs, variants) -> list:
        import networkx as nx
        from networkx.algorithms.isomorphism import GraphMatcher

        built = {}
        for name, (atoms, bonds) in graphs.items():
            v = variants[name]
            roles = {self.sites[v.residue][v.entry].atom: "entry"}
            roles[self.sites[v.residue][v.exit].atom] = roles.get(self.sites[v.residue][v.exit].atom, "") + "exit"
            g = nx.Graph()
            for a in atoms:
                g.add_node(a["key"], label=(a["element"], a.get("isotope"), a["charge"],
                                            a["electronic_state"]["radical_electrons"], a["lone_pairs"],
                                            a["hydrogens"], roles.get(a["key"]),
                                            a["key"].startswith(ARM_PREFIX)))
            for b in bonds:
                g.add_edge(*b["atoms"], order=float(b["order"]))
            built[name] = g
        kept = []
        for name in sorted(graphs, key=lambda n: (len(n), n)):
            if not any(GraphMatcher(built[name], built[other],
                                    node_match=lambda x, y: x["label"] == y["label"],
                                    edge_match=lambda x, y: x["order"] == y["order"]).is_isomorphic()
                       for other in kept):
                kept.append(name)
        return sorted(kept)

    def _reach_residues(self) -> int:
        """The most residues of one support path whose backbone atoms sum to
        at most r (module doc, word order)."""
        states = {(n, self.backbone[n]) for n in self.residues if self.backbone[n] <= self.radius}
        best = 0
        while states:
            best += 1
            states = {(b, total + self.backbone[b]) for a, total in states for b in self.support[a]
                      if total + self.backbone[b] <= self.radius}
        return best

    # -- environments -------------------------------------------------------------
    def _rows(self, word):
        asm = dc._Assembly()
        placed = []
        for i, name in enumerate(word):
            graph = self.residues[name].graph
            atoms = {a.key: a for a in graph.atoms}
            keys = tuple(sorted(atoms))
            handles = asm.add_period(atoms, keys, {k: j for j, k in enumerate(keys)},
                                     dc._bonds_of(graph), i, set(keys))
            if placed:
                left, right, order = self.junctions[(word[i - 1], name)]
                asm.bond(placed[-1][left], handles[right], order)
            placed.append(handles)
        asm.add_h(placed[0][self.entry[word[0]]], 4)
        asm.add_h(placed[-1][self.exit[word[-1]]], 4)
        rows, ids = asm.rows()
        return rows, [{k: ids[h] for k, h in table.items()} for table in placed]

    def _main_keys(self, name):
        return sorted(a.key for a in self.residues[name].graph.atoms
                      if a.element != "H" and not a.key.startswith(ARM_PREFIX))

    def _word_occurrences(self, word, positions, head_real=False, tail_real=False):
        """(touched word, position in it, touched indices, atom key, site key)
        of every non-stub heavy atom of ``word`` at ``positions``; the caps at
        a real end are vacant sites, the others out of every reach."""
        from rmgpy.kmc.rmg_moves import molecule_from_rows
        from rmgpy.kmc.site_class import extract_site_class

        rows, tables = self._rows(word)
        heavy = {i for i, row in rows.items() if row[0] != "H"}
        unit = {row: i for i, table in enumerate(tables) for row in table.values()}
        out = []
        for p in positions:
            for key in self._main_keys(word[p]):
                root = tables[p][key]
                distance, ball = de._ball(rows, heavy, root, self.radius + de.BALL_EXTRA)
                touched = sorted({unit[x] for x, d in distance.items() if d <= self.radius + 1})
                memo = (tuple(word[i] for i in touched), touched.index(p), key,
                        head_real and touched[0] == 0, tail_real and touched[-1] == len(word) - 1)
                if memo not in self._memo:
                    molecule = molecule_from_rows(ball)
                    atom = next(a for a in molecule.atoms if a.id == root)
                    self._memo[memo] = extract_site_class(molecule, [("site", atom)], self.policy)[0]
                out.append((memo[0], memo[1], touched, key, self._memo[memo]))
        return out

    def word_keys(self, word, position, head_real=False, tail_real=False) -> Dict[str, str]:
        """atom key -> site-class key of the non-stub heavy atoms of ``word[position]``."""
        return {key: site for _w, _p, _t, key, site in
                self._word_occurrences(tuple(word), [position], head_real, tail_real)}

    def context(self, name) -> Context:
        """``name`` with rho + 1 arm-free residues of the support on each side
        (the least such residue at each step)."""
        plain = [n for n in sorted(self.residues) if all(r == "H" for _s, r in self.variants[n].roles)]
        preceding = {n: [a for a in plain if n in self.support[a]] for n in self.residues}
        left, right = [], []
        for _ in range(self.rho + 1):
            first = left[0] if left else name
            left.insert(0, preceding[first][0])
            last = right[-1] if right else name
            right.append(next(b for b in self.support[last] if b in plain))
        return Context(tuple(left) + (name,) + tuple(right), self.rho + 1)

    def _words(self, length, first_ok):
        words = [(n,) for n in sorted(self.residues) if first_ok(n)]
        for _ in range(length - 1):
            words = [w + (b,) for w in words for b in self.support[w[-1]]]
        return words

    def environments(self, kinds: Sequence[str] = KINDS) -> Dict[str, list]:
        """site-class key -> sorted occurrences (kind, word, position, start,
        atom key) of the kinds of word asked for (module doc).  ``word`` is
        the touched residue word of an interior occurrence (``start`` its
        period index) and the word from the real end to its last touched
        residue otherwise (``start`` None)."""
        kinds = tuple(k for k in KINDS if k in kinds)
        if kinds in self._environments:
            return self._environments[kinds]
        found = {}
        if "interior" in kinds:
            period = self.period
            pad = self.rho + 1
            size = len(period)
            for start in range(size):
                word = tuple(period[(start - pad + i) % size] for i in range(2 * pad + 1))
                for touched_word, position, touched, key, site in self._word_occurrences(word, [pad]):
                    found.setdefault(site, set()).add(
                        ("interior", touched_word, position, (start - pad + touched[0]) % size, key))
        if "end" in kinds:
            for word in self._words(self.k + 1, self.vacant_entry):
                for _w, position, touched, key, site in self._word_occurrences(word, range(self.rho + 1), True):
                    if touched[0] == 0:
                        found.setdefault(site, set()).add(
                            ("end", word[:touched[-1] + 1], touched[position], None, key))
        if "short" in kinds:
            for length in range(1, self.k + 1):
                for word in self._words(length, self.vacant_entry):
                    if not self.vacant_exit(word[-1]):
                        continue
                    for _w, position, touched, key, site in self._word_occurrences(
                            word, range(length), True, True):
                        found.setdefault(site, set()).add(("short", word, touched[position], None, key))
        out = {site: sorted(items, key=lambda o: (len(o[1]), o[0], o[1], o[2], str(o[3]), o[4]))
               for site, items in sorted(found.items())}
        self._environments[kinds] = out
        return out

    @property
    def period(self) -> tuple:
        if "_period" not in self.__dict__:
            self.__dict__["_period"] = de.covering_period(self.support, self.k)
        return self.__dict__["_period"]

    def environment_witnesses(self):
        """``(period, witnesses)`` of the interior words, as
        ``declared_environments.environment_witnesses`` gives them for a
        linear declaration: per key, (word, position, atom key, start),
        shortest word first."""
        witnesses = {}
        for site, occurrences in self.environments(("interior",)).items():
            items = {(word, position, key, start) for _kind, word, position, start, key in occurrences}
            witnesses[site] = sorted(items, key=lambda item: (len(item[0]),) + item)
        return self.period, witnesses

    def h_bearing(self, kind, word, position, key) -> bool:
        """Whether atom ``key`` of ``word[position]`` carries a hydrogen in
        that word (a vacant end site's cap counts at a real end)."""
        name = word[position]
        atom = next(a for a in self.residues[name].graph.atoms if a.key == key)
        if atom.element == "H":
            return False
        if atom.hydrogens:
            return True
        if kind in ("end", "short") and position == 0 and key == self.entry[name]:
            return True
        return kind == "short" and position == len(word) - 1 and key == self.exit[name]


def branched_lift(spec, stub_depth: Optional[int] = None) -> BranchedLift:
    """The lift of a ``BranchedPolymerSpec`` (module doc), cached per
    declaration and stub depth."""
    from rmgpy.polymer_spec import BranchedPolymerSpec, normalized_polymer_spec

    spec = normalized_polymer_spec(spec)
    if not isinstance(spec, BranchedPolymerSpec):
        raise TypeError(f"a branched lift takes a BranchedPolymerSpec, got {type(spec).__name__}")
    key = (hashlib.sha256(spec.normalized_json).hexdigest(), stub_depth)
    if key not in _LIFTS:
        _LIFTS[key] = BranchedLift(spec, stub_depth)
    return _LIFTS[key]


_LIFTS: Dict[tuple, BranchedLift] = {}


def branched_environments(spec, stub_depth: Optional[int] = None,
                          kinds: Sequence[str] = KINDS) -> Dict[str, list]:
    """Every site-class environment of the branched declaration (module doc)."""
    return branched_lift(spec, stub_depth).environments(kinds)


def h_bearing_environments(spec, environments: Mapping[str, list]) -> set:
    """The keys of ``environments`` one of whose occurrences carries a hydrogen."""
    lift = branched_lift(spec)
    return {site for site, occurrences in environments.items()
            if any(lift.h_bearing(kind, word, position, key) for kind, word, position, _s, key in occurrences)}


def _formula(graph, extra_h=0) -> tuple:
    counts = {}
    for atom in graph.atoms:
        counts[atom.element] = counts.get(atom.element, 0) + 1
        if atom.hydrogens:
            counts["H"] = counts.get("H", 0) + atom.hydrogens
    if extra_h:
        counts["H"] = counts.get("H", 0) + extra_h
    return tuple(sorted(counts.items()))


class BranchedChainProvider(SequenceChainProvider):
    """The sequence chain provider of a branched declaration's lift (module
    doc), registered under ``branched:<sha256[:16]>`` of the declaration."""

    def __init__(self, spec):
        self.branched_spec = spec
        self.lift = branched_lift(spec)
        super().__init__(self.lift, name="branched:" + hashlib.sha256(spec.normalized_json).hexdigest()[:16])

    def _apply(self, asm, edit, chain, sites, port_handles, cap_handles):
        if 0 <= edit.unit < len(chain.residues) and edit.key.startswith(ARM_PREFIX):
            raise dc.ChainSpecError(f"{STUB_REFUSAL}: {edit}: ({edit.unit}, {edit.key!r}) is a stub atom; "
                                    "its environments belong to other words and it is never a radical site")
        return super()._apply(asm, edit, chain, sites, port_handles, cap_handles)

    def unit_formula(self) -> tuple:
        """The tree increments of the residues (module doc): formula with every
        site occupied plus (sites - 2) H, the set sorted."""
        return tuple(sorted({_formula(residue.graph, len(residue.sites) - 2)
                             for _name, residue in self.branched_spec.residues}))

    def growth_sites(self, molecule):
        """None: branch growth by construction is not built (module doc)."""
        return None


def branched_provider(spec) -> BranchedChainProvider:
    """The registered provider of this exact normalized branched declaration."""
    name = "branched:" + hashlib.sha256(spec.normalized_json).hexdigest()[:16]
    found = PROVIDERS.get(name)
    if isinstance(found, BranchedChainProvider) and \
            found.branched_spec.normalized_json == spec.normalized_json:
        return found
    return BranchedChainProvider(spec)

