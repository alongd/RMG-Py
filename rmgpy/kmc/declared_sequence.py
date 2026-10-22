"""Chain witnesses of a linear copolymer declaration (schema v2).

``SequenceChainProvider`` builds the chains of one ``PolymerSpec``
(``rmgpy/polymer_spec.py``) from the declared residues, links and end
groups; nothing here names a polymer.  It is the copolymer counterpart of
``declared_chain.DeclaredChainProvider`` and is registered in
``boundaries.PROVIDERS`` under ``sequence:<sha256[:16]>`` of the
normalized declaration.

Builds (``SequenceChainSpec``):

- ``residues``: one declared residue name per unit, in chain order; every
  step ``a -> b`` must be in the sequence support, and residue ``b`` is
  joined to ``a`` by the declared junction (exit atom of ``a`` to entry atom
  of ``b``, ``PolymerSpec.junctions``).
- ``head`` / ``tail``: ``artificial`` (H caps on the port, port
  artificial), ``declared`` (the declared end group attached, port
  physical) or ``radical`` (port physical, the port atom keeps the cap's
  electrons), as in ``declared_chain``.  The head port is the entry atom of
  the first residue, the tail port the exit atom of the last.
- ``edits``: ``declared_chain.EditOp``, addressing ``(unit, residue atom
  key)``; a declared end group's heavy atom is ``(-1, key)`` at the head
  and ``(units, key)`` at the tail.

Coordinates are ``(unit, position)`` with ``position`` the index of the key
in that residue's sorted keys; ``BoundaryMetadata.residues`` carries the
residue name of every unit.  Aromatic (1.5) residues are built in one
deterministic Kekule form per residue, as ``declared_chain`` builds its
periods.

**Extension names its residue**: the residue beyond an artificial port is
not fixed by the declaration, so one extension step adds a NAMED residue,
one the sequence support allows there (``next_residues``; ``fragment`` with
``residues``), and a step without one is refused
(``refused:sequence_extension``), as is an end edited onto a chain
(``end_edit``).  The window fixed point tries every allowed residue
(``windows.fixed_point``, residue branches).

**The growth unit** (``unit_formula``) of the depository exclusion is the
set of residue formulas, one per residue: a long-chain analogue gains any
non-negative integer combination of them (``certificates.growth_units``).
The exclusion is decided by construction at the growth sites
(``growth_sites``): every capped terminal residue of a witness species, read
as every residue and side it matches from which the whole species parses
under the sequence support (``_parses``: residue by residue on the heavy
skeleton, ending on a whole residue, the end group or a part of a residue a
reaction cut; never a join checked by element alone), on either side of the
reaction (``certificates.reaction_growth_sites``); at each reading every residue word
the sequence support allows next to it is an insertion
(``growth_options``), built whole by ``growth_fragment``; an extension step
is one such insertion of one residue at a port."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Mapping

from rmgpy.kmc import declared_chain as dc
from rmgpy.kmc.boundaries import (
    PROVIDERS,
    BoundaryMetadata,
    BoundaryProvider,
    ExtensionRefused,
    Fragment,
    MetadataError,
    Port,
    ProxyBuild,
)

SEQUENCE_EXTENSION_REFUSAL = "refused:sequence_extension"


@dataclass(frozen=True)
class SequenceChainSpec:
    """A chain of declared residues (module doc)."""

    residues: tuple
    head: str = "artificial"
    tail: str = "artificial"
    edits: tuple = ()

    def __post_init__(self):
        residues = tuple(self.residues)
        if not residues or not all(isinstance(r, str) and r for r in residues):
            raise dc.ChainSpecError("a sequence chain needs one or more residue names")
        object.__setattr__(self, "residues", residues)
        for end in (self.head, self.tail):
            if end not in dc.END_KINDS:
                raise dc.ChainSpecError(f"end {end!r} is not one of {dc.END_KINDS}")
        edits = tuple(dc.EditOp.from_mapping(e) if isinstance(e, Mapping) else e
                      for e in self.edits)
        if any(not isinstance(e, dc.EditOp) for e in edits):
            raise dc.ChainSpecError("edits must be EditOp entries")
        object.__setattr__(self, "edits", edits)

    @classmethod
    def from_mapping(cls, value: Mapping) -> "SequenceChainSpec":
        return dc._from_mapping(cls, value)


def _formula(graph) -> dict:
    counts = {}
    for atom in graph.atoms:
        counts[atom.element] = counts.get(atom.element, 0) + 1
        if atom.hydrogens:
            counts["H"] = counts.get("H", 0) + atom.hydrogens
    return counts


def _embed(elements, adjacent, root, inward, t, j, target, heavy, available):
    """Every map of a residue skeleton (key -> element, key -> neighbour
    keys) into the heavy atoms ``available`` of ``target`` rows (``heavy``:
    atom -> heavy neighbours), ``root`` on ``t``, keeping elements and bonds,
    whose image atoms have no heavy neighbour but the images of their
    skeleton neighbours, ``j`` at the root and at most one more at
    ``inward``: yields (image, that atom or None)."""
    if target[t][0] != elements[root]:
        return
    order, parent = [root], {}
    for k in order:
        for r in sorted(adjacent[k]):
            if r not in parent and r != root:
                parent[r] = k
                order.append(r)
    if len(order) != len(elements):
        return
    image, used = {root: t}, {t}

    def finish():
        beyond = None
        for k, a in image.items():
            extra = set(heavy[a]) - {image[r] for r in adjacent[k]}
            if k == root:
                extra.discard(j)
            if k == inward and len(extra) <= 1:
                beyond = next(iter(extra), None)
            elif extra:
                return
        yield dict(image), beyond

    def extend(i):
        if i == len(order):
            yield from finish()
            return
        k = order[i]
        for a in sorted(heavy[image[parent[k]]]):
            if a in used or a not in available or target[a][0] != elements[k] or \
                    any(r in image and image[r] not in heavy[a] for r in adjacent[k]):
                continue
            image[k] = a
            used.add(a)
            yield from extend(i + 1)
            del image[k]
            used.discard(a)

    yield from extend(1)


def _cut(elements, adjacent, root, t, target, heavy, available) -> bool:
    """Whether the heavy atoms ``available`` (connected through ``t``) are a
    connected part of a residue skeleton with ``t`` on its ``root``: elements
    kept, every bond among them a skeleton bond."""
    if target[t][0] != elements[root]:
        return False
    order, parent = [t], {}
    for a in order:
        for n in sorted(heavy[a]):
            if n in available and n not in parent and n != t:
                parent[n] = a
                order.append(n)
    image, used = {t: root}, {root}

    def extend(i):
        if i == len(order):
            return True
        a = order[i]
        for k in sorted(adjacent[image[parent[a]]]):
            if k in used or elements[k] != target[a][0] or \
                    any(n in image and image[n] not in adjacent[k] for n in heavy[a]):
                continue
            image[a] = k
            used.add(k)
            if extend(i + 1):
                return True
            del image[a]
            used.discard(k)
        return False

    return extend(1)


class SequenceChainProvider(BoundaryProvider):
    """Chains of one linear copolymer declaration (module doc)."""

    RADICAL_END_CHEMISTRIES = ("radical",)

    def __init__(self, spec, name=None):
        from rmgpy.kmc.declared_environments import _backbone_atoms, _ends

        self.spec = spec
        self.name = name or "sequence:" + hashlib.sha256(spec.normalized_json).hexdigest()[:16]
        self.residues = dict(spec.residues)
        self.junctions = spec.junctions
        self.support = {a: tuple(successors) for a, successors in spec.sequence_support}
        self.entry, self.exit = _ends(self.junctions)
        self.keys = {name: tuple(sorted(atom.key for atom in residue.graph.atoms))
                     for name, residue in self.residues.items()}
        self.bonds = {}
        for name, residue in self.residues.items():
            atoms = {atom.key: atom for atom in residue.graph.atoms}
            sites = {site.atom: site.bond_order for _s, site in residue.sites}
            self.bonds[name] = dc._kekule(atoms, dc._bonds_of(residue.graph), sites)
        self.backbone = {name: _backbone_atoms(residue, self.entry[name], self.exit[name])
                         for name, residue in self.residues.items()}
        self.coordinate_backbone_atoms = min(self.backbone.values())
        self.END_CHEMISTRIES = {"backbone_head": ("declared", "radical"),
                                "backbone_tail": ("declared", "radical")}
        PROVIDERS[self.name] = self

    # -- builds ----------------------------------------------------------------
    def _check(self, chain) -> None:
        if not isinstance(chain, SequenceChainSpec):
            raise dc.ChainSpecError(f"a sequence chain needs a SequenceChainSpec, got {chain!r}")
        for name in chain.residues:
            if name not in self.residues:
                raise dc.ChainSpecError(f"{name!r} is not a declared residue "
                                        f"({sorted(self.residues)})")
        for a, b in zip(chain.residues, chain.residues[1:]):
            if b not in self.support[a]:
                raise dc.ChainSpecError(f"the step {a}->{b} is not in the sequence support")

    def _order(self, name, atom) -> float:
        return next(site.bond_order for _s, site in self.residues[name].sites if site.atom == atom)

    def build(self, name: str, chain: SequenceChainSpec) -> ProxyBuild:
        """The chain of ``chain`` (module doc); RMG's default multiplicity."""
        from rmgpy.kmc.rmg_moves import molecule_from_rows

        self._check(chain)
        asm = dc._Assembly()
        sites, placed = {}, []
        for unit, residue in enumerate(chain.residues):
            keys = self.keys[residue]
            atoms = {atom.key: atom for atom in self.residues[residue].graph.atoms}
            handles = asm.add_period(atoms, keys, {k: i for i, k in enumerate(keys)},
                                     self.bonds[residue], unit, set(keys))
            if placed:
                left, right, order = self.junctions[(chain.residues[unit - 1], residue)]
                asm.bond(placed[-1][left], handles[right], order)
            placed.append(handles)
            for key, handle in handles.items():
                sites[(unit, key)] = handle
        units = len(chain.residues)
        first, last = chain.residues[0], chain.residues[-1]
        ports = []
        for side, end, handle, order, tag in (
                ("head", chain.head, placed[0][self.entry[first]],
                 self._order(first, self.entry[first]), "backbone_head"),
                ("tail", chain.tail, placed[-1][self.exit[last]],
                 self._order(last, self.exit[last]), "backbone_tail")):
            caps = []
            if end == "artificial":
                caps = [asm.add_h(handle, 4) for _ in range(int(order))]
            elif end == "radical":
                asm.attrs[handle][3] += int(order)
            else:
                group_unit = -1 if side == "head" else units
                # the declared end group, attached as declared_chain attaches it
                for key, atom in dc.DeclaredChainProvider._end_group(
                        self, asm, side, handle, order, group_unit).items():
                    sites[(group_unit, key)] = atom
            ports.append((side, handle, caps, tag,
                          "artificial" if end == "artificial" else "physical"))
        port_handles = {h for _s, h, _c, _t, _p in ports}
        cap_handles = {c for _s, _h, caps, _t, _p in ports for c in caps}
        for edit in chain.edits:
            self._apply(asm, edit, chain, sites, port_handles, cap_handles)
        asm.check_valence(name)
        asm.check_connected(name)
        rows, ids = asm.rows()
        built_ports = tuple(Port(side, ids[h], tuple(ids[c] for c in caps), tag, physicality)
                            for side, h, caps, tag, physicality in ports)
        coordinates = tuple((ids[h], u, k) for h, (u, k) in asm.coord.items())
        metadata = BoundaryMetadata(self.name, built_ports, coordinates, residues=chain.residues)
        metadata.check(rows, name)
        return ProxyBuild(name, molecule_from_rows(rows), metadata, tuple(range(len(rows))))

    def _apply(self, asm, edit, chain, sites, port_handles, cap_handles):
        units = len(chain.residues)

        def locate(unit, key):
            if unit in (-1, units) and (unit, key) in sites:
                return sites[(unit, key)]
            if not 0 <= unit < units:
                raise dc.EditError(f"{edit}: unit {unit} is outside the chain (0..{units - 1})")
            residue = chain.residues[unit]
            if key not in self.keys[residue]:
                raise dc.EditError(f"{edit}: key {key!r} is not an atom of residue {residue} "
                                   f"{list(self.keys[residue])}")
            handle = sites.get((unit, key))
            if handle is None or handle not in asm.attrs:
                raise dc.EditError(f"{edit}: atom ({unit}, {key!r}) is not in this chain")
            return handle

        handle = locate(edit.unit, edit.key)
        if edit.kind == "set_radical":
            asm.attrs[handle][3] = int(edit.value)
        elif edit.kind == "add_h":
            for _ in range(edit.count):
                asm.add_h(handle, 2)
        elif edit.kind == "remove_h":
            free = [h for h in asm.hydrogens(handle) if h not in cap_handles]
            if len(free) < edit.count:
                raise dc.EditError(f"{edit}: atom ({edit.unit}, {edit.key!r}) has {len(free)} "
                                   "removable hydrogens")
            for h in free[len(free) - edit.count:]:
                asm.remove(h)
        elif edit.kind == "set_bond_order":
            other = locate(*edit.other)
            if other not in asm.adj[handle]:
                raise dc.EditError(f"{edit}: ({edit.unit}, {edit.key!r}) and {edit.other} "
                                   "are not bonded")
            asm.adj[handle][other] = asm.adj[other][handle] = float(edit.value)
        elif edit.kind == "delete_atom":
            if handle in port_handles:
                raise dc.EditError(f"{edit}: ({edit.unit}, {edit.key!r}) carries a port")
            for h in [h for h in asm.adj[handle] if asm.attrs[h][0] == "H"
                      and len(asm.adj[h]) == 1]:
                asm.remove(h)
            asm.remove(handle)

    def position(self, residue: str, key: str) -> int:
        """The coordinate position of ``key`` in ``residue``."""
        return self.keys[residue].index(key)

    # -- what the windows read ---------------------------------------------------
    def declared_elements(self) -> frozenset:
        """Every element the declaration names: residues, end groups,
        monomer products and small species."""
        spec = self.spec
        graphs = [residue.graph for residue in self.residues.values()]
        graphs += [spec.head_end_group, spec.tail_end_group]
        graphs += [graph for _name, graph in spec.monomer_products]
        graphs += [graph for _name, graph, _disposition in spec.small_species]
        return frozenset(atom.element for graph in graphs for atom in graph.atoms)

    def unit_formula(self) -> tuple:
        """The depository growth units (module doc): one formula per residue,
        each ((element, count), ...) sorted, the set sorted."""
        return tuple(sorted({tuple(sorted(_formula(residue.graph).items()))
                             for residue in self.residues.values()}))

    # -- depository growth (module doc) ------------------------------------------
    def _preceding(self) -> dict:
        return {name: tuple(a for a in sorted(self.support) if name in self.support[a])
                for name in self.residues}

    def _growth_patterns(self) -> tuple:
        """(side, residue, port, pattern rows, exact atoms, join neighbour,
        joined atom) of every capped terminal residue: the residue at the
        artificial ``side`` of a two-residue chain with each neighbour the
        support allows, and of the residue alone; its atoms with every
        hydrogen exact, the join neighbour matched by element and bond order
        (``declared_chain``'s end patterns), the rest of the species parsed
        from it (``_parses``).  The residue alone has no join neighbour
        (None, None): its exact atoms are the whole species."""
        cached = self.__dict__.get("_growth")
        if cached is None:
            from rmgpy.kmc.rmg_moves import molecule_rows
            preceding = self._preceding()
            cached = []
            for name in sorted(self.residues):
                words = [("tail", (q, name)) for q in preceding[name]]
                words += [("head", (name, q)) for q in self.support[name]]
                words += [("head", (name,)), ("tail", (name,))]
                for side, word in words:
                    build = self.build("growth", SequenceChainSpec(word))
                    rows = molecule_rows(build.molecule)
                    unit = 0 if side == "head" else len(word) - 1
                    atoms = {a for a, u, _k in build.metadata.coordinates if u == unit}
                    exact = atoms | {h for a in atoms for h, _o in rows[a][5] if rows[h][0] == "H"}
                    stubs = {n for a in exact for n, _o in rows[a][5] if n not in exact}
                    stub = min(stubs) if stubs else None
                    inner = None if stub is None else \
                        next(n for n, _o in rows[stub][5] if n in atoms)
                    cached.append((side, name, build.metadata.port(side),
                                   {a: rows[a] for a in exact | stubs}, frozenset(exact),
                                   stub, inner))
            cached = tuple(cached)
            self.__dict__["_growth"] = cached
        return cached

    def growth_sites(self, molecule) -> tuple:
        """Where residues insert in ``molecule`` (atoms carry ids): every atom
        that is the capped end of a whole terminal residue (``_growth_patterns``)
        from which the species parses under the sequence support
        (``_parses``), one site per end atom; its coordinate lists every
        reading ``(side, residue, caps)`` it matches, the site's port caps are
        all their caps.  A radical end, a vinyl end, a species without a whole
        residue or an end whose species does not parse is no site."""
        from rmgpy.kmc.rmg_moves import molecule_rows
        target = molecule_rows(molecule)
        heavy = {a: frozenset(n for n, _o in row[5] if target[n][0] != "H")
                 for a, row in target.items() if row[0] != "H"}
        memo = {}
        readings = {}
        for side, name, port, pattern, exact, stub, inner in self._growth_patterns():
            for anchor in sorted(target):
                if (side, name) in readings.get(anchor, {}):
                    continue
                image = dc._match_end(pattern, exact, port.atom_id, target, anchor)
                if image is not None and stub is not None and not self._parses(
                        target, heavy, image[stub], image[inner], name, side, memo):
                    image = None
                if image is not None:
                    readings.setdefault(anchor, {})[(side, name)] = \
                        tuple(sorted(image[c] for c in port.cap_atom_ids))
        return tuple(self._site(anchor, readings[anchor]) for anchor in sorted(readings))

    @staticmethod
    def _site(anchor, readings) -> tuple:
        """(Port, coordinate) of a growth site from its readings
        {(side, residue): caps}: the coordinate every reading sorted, the
        port's caps all their caps."""
        coordinate = tuple((side, name, caps) for (side, name), caps in sorted(readings.items()))
        caps = tuple(sorted({c for _s, _n, cs in coordinate for c in cs}))
        return (Port(f"growth:{anchor}", anchor, caps, f"backbone_{coordinate[0][0]}",
                     "artificial"), coordinate)

    def _skeletons(self) -> dict:
        """The heavy-atom skeleton (key -> element, key -> heavy neighbour
        keys) of every residue, and of the head and tail end groups
        (``("end", side)``, with their port atom)."""
        cached = self.__dict__.get("_skeleton")
        if cached is None:
            cached = {}
            graphs = [(name, residue.graph, None) for name, residue in self.residues.items()]
            graphs += [(("end", side), graph, graph.ports[0][1].atom) for side, graph in
                       (("head", self.spec.head_end_group), ("tail", self.spec.tail_end_group))]
            for name, graph, port in graphs:
                elements = {atom.key: atom.element for atom in graph.atoms if atom.element != "H"}
                adjacent = {key: set() for key in elements}
                for bond in graph.bonds:
                    x, y = bond.atoms
                    if x in elements and y in elements:
                        adjacent[x].add(y)
                        adjacent[y].add(x)
                cached[name] = (elements, {k: frozenset(v) for k, v in adjacent.items()}, port)
            self.__dict__["_skeleton"] = cached
        return cached

    def _parses(self, target, heavy, t, j, prev, side, memo) -> bool:
        """Whether the heavy atoms beyond the bond ``j``-``t`` (``j`` in
        residue ``prev``, read at the ``side`` end) parse as the rest of a
        chain of the declaration: residue after residue inward (before
        ``prev`` toward the head from a tail end, after it toward the tail from
        a head end), each one the sequence support allows next to the last,
        joined at its declared entry and exit atoms, matched on the heavy
        skeleton (elements and bonds: the reaction may change hydrogens,
        radicals and bond orders); the parse ends where the atoms do, on a
        whole residue, on the declared end group, or on a connected part of
        one residue (a residue cut by the reaction).  Never by element alone:
        an end whose neighbours are not the residue the support allows there
        is no residue end."""
        key = (t, j, prev, side)
        if key in memo:
            return memo[key]
        memo[key] = False
        available, frontier = {t}, [t]
        while frontier:
            a = frontier.pop()
            for n in heavy[a]:
                if n != j and n not in available:
                    available.add(n)
                    frontier.append(n)
        skeletons = self._skeletons()
        found = False
        for q in (self._preceding()[prev] if side == "tail" else self.support[prev]):
            elements, adjacent, _port = skeletons[q]
            root, inward = (self.exit[q], self.entry[q]) if side == "tail" else \
                (self.entry[q], self.exit[q])
            for image, beyond in _embed(elements, adjacent, root, inward, t, j, target, heavy,
                                        available):
                if beyond is None or self._parses(target, heavy, beyond, image[inward], q,
                                                  side, memo):
                    found = True
                    break
            if found or (len(available) < len(elements) and
                         _cut(elements, adjacent, root, t, target, heavy, available)):
                found = True
                break
        if not found:
            elements, adjacent, port = skeletons[("end", "head" if side == "tail" else "tail")]
            found = port in elements and len(elements) == len(available) and any(
                beyond is None for _image, beyond in
                _embed(elements, adjacent, port, None, t, j, target, heavy, available))
        memo[key] = found
        return found

    def merged_site(self, first, second) -> tuple:
        """One growth site read on both sides of a reaction
        (``certificates.reaction_growth_sites``): every reading of either,
        a reading on both keeping ``first``'s caps."""
        readings = {(side, name): caps for side, name, caps in second[1]}
        readings.update({(side, name): caps for side, name, caps in first[1]})
        return self._site(first[0].atom_id, readings)

    def _site_order(self, name, side) -> float:
        return self._order(name, (self.entry if side == "head" else self.exit)[name])

    def growth_options(self, coordinate, budget):
        """Every insertion at a growth site (``growth_sites``): per reading
        (side, residue, caps) of the site, every non-empty word of residues
        the sequence support allows next to that residue on that side
        (appended after it at a tail, prepended before it at a head; chain
        order), adding at most ``budget`` heavy atoms.  Yields ((reading,
        word), formula added, heavy atoms added), shorter words first; a
        word whose junction order is not the reading's cap count is no
        insertion."""
        preceding = self._preceding()
        formulas = {name: _formula(residue.graph) for name, residue in self.residues.items()}
        heavy = {name: sum(n for e, n in f.items() if e != "H") for name, f in formulas.items()}
        for reading in coordinate:
            side, name, caps = reading
            frontier = [((), 0)]
            while frontier:
                grown = []
                for word, size in frontier:
                    edge = name if not word else word[-1] if side == "tail" else word[0]
                    for nxt in (self.support[edge] if side == "tail" else preceding[edge]):
                        if size + heavy[nxt] > budget:
                            continue
                        new = word + (nxt,) if side == "tail" else (nxt,) + word
                        grown.append((new, size + heavy[nxt]))
                for word, size in grown:
                    first, last = (name, word[0]) if side == "tail" else (word[-1], name)
                    if int(self.junctions[(first, last)][2]) != len(caps):
                        continue
                    added = {}
                    for r in word:
                        for e, n in formulas[r].items():
                            added[e] = added.get(e, 0) + n
                    outer = word[-1] if side == "tail" else word[0]
                    caps_added = int(self._site_order(outer, side)) - len(caps)
                    if caps_added:
                        added["H"] = added.get("H", 0) + caps_added
                    yield (reading, word), added, size
                frontier = grown

    def fragment_facts(self) -> dict:
        """The facts of one growth step the window guards read
        (``windows._fragment_facts``): per side, residue end and residue the
        support allows next to it, whether that residue's fragment
        (``growth_fragment``) is closed shell and its net element delta
        (fresh atoms minus the caps it replaces)."""
        preceding = self._preceding()
        out = {}
        for name in sorted(self.residues):
            for side, nexts in (("head", preceding[name]), ("tail", self.support[name])):
                caps = tuple(range(int(self._site_order(name, side))))
                for nxt in nexts:
                    fragment = self.growth_fragment(((side, name, caps), (nxt,)))
                    delta = {}
                    for atom in fragment.atoms:
                        delta[atom[0]] = delta.get(atom[0], 0) + 1
                    delta["H"] = delta.get("H", 0) - len(caps)
                    out[f"backbone_{side}:{name}:{nxt}"] = {
                        "closed_shell": all(atom[3] == 0 for atom in fragment.atoms),
                        "delta": {e: n for e, n in delta.items() if n}}
        return out

    def growth_fragment(self, option) -> Fragment:
        """The word of a ``growth_options`` insertion as one fragment: its
        residues joined as a chain, the end facing the site bonded to the
        site's atom by the declared junction, the far end capped."""
        (side, name, caps), word = option
        key = (side, name, len(caps), word)
        cache = self.__dict__.setdefault("_fragments", {})
        if key not in cache:
            from rmgpy.kmc.rmg_moves import molecule_rows
            build = self.build("growth", SequenceChainSpec(word))
            rows = molecule_rows(build.molecule)
            attach = build.metadata.port("head" if side == "tail" else "tail")
            new = build.metadata.port(side)
            first, last = (name, word[0]) if side == "tail" else (word[-1], name)
            order = float(self.junctions[(first, last)][2])
            kept = sorted(a for a in rows if a not in attach.cap_atom_ids)
            local = {a: i for i, a in enumerate(kept)}
            cache[key] = Fragment(
                tuple(tuple(rows[a][:5]) for a in kept),
                tuple(sorted((local[a], local[b], float(o)) for a in kept for b, o in rows[a][5]
                             if b in local and local[a] < local[b])),
                ((local[attach.atom_id], order),), local[new.atom_id],
                tuple(local[c] for c in new.cap_atom_ids),
                tuple((local[a], u, k) for a, u, k in build.metadata.coordinates if a in local))
        return cache[key]

    def _neighbours(self, name) -> tuple:
        before = min(a for a in sorted(self.support) if name in self.support[a])
        return before, min(self.support[name])

    def interior_atoms_rmg(self) -> list:
        """Perceived RMG atoms of every residue in the middle of a
        three-residue artificial chain (its declared neighbours), with their
        hydrogens."""
        out = []
        for name in sorted(self.residues):
            before, after = self._neighbours(name)
            build = self.build("interior", SequenceChainSpec((before, name, after)))
            molecule = build.molecule
            molecule.update()
            middle = {a for a, u, _k in build.metadata.coordinates if u == 1}
            atoms = [a for a in molecule.atoms if a.id in middle]
            out += atoms + [h for a in atoms for h in a.bonds if h.is_hydrogen()]
        return out

    def probe_ports(self) -> tuple:
        """(port, coordinate) of the head and tail ports of a two-residue
        artificial chain."""
        name = min(self.residues)
        build = self.build("probe", SequenceChainSpec((name, self._neighbours(name)[1])))
        return tuple((port, build.metadata.coordinate(port.atom_id))
                     for port in build.metadata.ports)

    def witness_backbone_indices(self, rows, metadata: BoundaryMetadata) -> dict:
        """Atom id -> backbone index for the heavy atoms on a shortest
        head-port -> tail-port path of the participant (the distance from the
        head port, counted from the backbone atoms of the residues before it)."""
        port_ids = {port.port_id for port in metadata.ports}
        if not {"head", "tail"} <= port_ids or metadata.residues is None:
            raise MetadataError(f"refused:backbone_ports: a participant with ports "
                                f"{sorted(port_ids)} has no head-to-tail backbone")
        from_head, from_tail, span = dc.DeclaredChainProvider._path_over_rows(rows, metadata)
        unit, _position = metadata.coordinate(metadata.port("head").atom_id)
        base = sum(self.backbone[name] for name in metadata.residues[:max(unit, 0)])
        return {a: base + d for a, d in from_head.items()
                if a in from_tail and d + from_tail[a] == span}

    def backbone_index(self, unit: int, position: int, build=None):
        if build is None:
            raise MetadataError("refused:sequence_coordinates: a sequence chain has a backbone "
                                "index only on a built chain")
        from rmgpy.kmc.boundaries import rows_of_molecule

        indices = self.witness_backbone_indices(dict(rows_of_molecule(build.molecule)),
                                                build.metadata)
        for a, u, k in build.metadata.coordinates:
            if (u, k) == (unit, position):
                return indices.get(a)
        return None

    def next_residues(self, residue, tag, caps) -> tuple:
        """The residues one extension step may add at a port on ``residue``
        (module doc): after it at a tail, before it at a head, every residue
        the sequence support allows there whose junction order is the port's
        cap count."""
        if tag == "backbone_tail":
            pairs = [((residue, nxt), nxt) for nxt in self.support[residue]]
        elif tag == "backbone_head":
            pairs = [((nxt, residue), nxt) for nxt in self._preceding()[residue]]
        else:
            return ()
        return tuple(sorted(nxt for pair, nxt in pairs if int(self.junctions[pair][2]) == caps))

    def next_words(self, residue, tag, caps, length) -> tuple:
        """Every word of ``length`` residues that ``length`` extension steps
        may add at a port on ``residue`` (``next_residues`` step by step, each
        new port's cap count its residue's site order), innermost first."""
        side = "tail" if tag == "backbone_tail" else "head"
        words = [((), residue, caps)]
        for _step in range(length):
            words = [(word + (nxt,), nxt, int(self._site_order(nxt, side)))
                     for word, last, count in words
                     for nxt in self.next_residues(last, tag, count)]
        return tuple(word for word, _last, _count in words)

    def fragment(self, port, coordinate, residues=None):
        """One residue at an artificial port: ``residues`` is (the residue at
        the port, the residue added), an allowed step (``next_residues``)
        joined by its declared junction, the new residue one unit outward.
        Without ``residues`` the step is refused: the declaration does not fix
        the residue beyond the port."""
        if residues is None:
            raise ExtensionRefused(f"{SEQUENCE_EXTENSION_REFUSAL}: {self.name}: the residue "
                                   "beyond an artificial port of a copolymer chain is not fixed "
                                   "by the declaration")
        name, nxt = residues
        if nxt not in self.residues:
            raise ExtensionRefused(f"{self.name}: {nxt!r} is not a declared residue "
                                   f"({sorted(self.residues)})")
        if nxt not in self.next_residues(name, port.tag, len(port.cap_atom_ids)):
            raise ExtensionRefused(f"{self.name}: {nxt} is not allowed next to {name} at a "
                                   f"{port.tag} port with {len(port.cap_atom_ids)} caps")
        side = "tail" if port.tag == "backbone_tail" else "head"
        step = self.growth_fragment(((side, name, tuple(port.cap_atom_ids)), (nxt,)))
        unit = coordinate[0] + (1 if side == "tail" else -1)
        return Fragment(step.atoms, step.bonds, step.attachments, step.new_port, step.new_caps,
                        tuple((i, unit, k) for i, _u, k in step.coordinates))

    def end_edit(self, metadata, port, end, rows, protected):
        raise ExtensionRefused(f"{SEQUENCE_EXTENSION_REFUSAL}: {self.name}: an end of a "
                               "copolymer chain is built, never edited onto an extended chain")


def sequence_provider(spec) -> SequenceChainProvider:
    """The registered provider of this exact normalized v2 declaration (a
    branched one's is its lift's, ``declared_branches.branched_provider``)."""
    from rmgpy.polymer_spec import BranchedPolymerSpec, normalized_polymer_spec

    spec = normalized_polymer_spec(spec)
    if isinstance(spec, BranchedPolymerSpec):
        from rmgpy.kmc.declared_branches import branched_provider
        return branched_provider(spec)
    name = "sequence:" + hashlib.sha256(spec.normalized_json).hexdigest()[:16]
    found = PROVIDERS.get(name)
    if isinstance(found, SequenceChainProvider) and \
            found.spec.normalized_json == spec.normalized_json:
        return found
    return SequenceChainProvider(spec)
