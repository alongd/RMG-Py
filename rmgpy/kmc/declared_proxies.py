"""The compile proxies of a declared linear homopolymer.

``declared_proxy_set(spec, units)`` builds the centred finite structures the
compiler generates family reactions on from the declaration alone, with the
builders the emit roster uses (``emit.reference_systems``): the declaration's
own chain provider (``emit.declared_provider``) through the vocabulary
adapter of the provider whose repeat unit it declares (``emit.match_provider``,
``spec_for_declared_provider``), and the monomer the declaration names or
derives (``reference_monomer``).  Nothing here names a polymer.
``declared_structures`` gives the same structures with their witness specs;
the roster of a declared artifact is built from it
(``emit.declared_reference_systems``).

The structures, at ``units`` repeat units with artificial (capped) ends:

- ``pristine``: the closed chain;
- one radical per class of hydrogen-bearing heavy atoms of the primitive
  period (``radical_sites``: the orbits of ``ChainGeometry.atom_classes``):
  ``interior_radical:<i>`` for a backbone class, on the centre unit
  (``emit.centred_radical_unit``), and ``pendant_radical:<i>`` for a
  side-group class, on the first unit; ``i`` counts the classes of that
  kind in order of their least key;
- ``doubly_featured``: radicals on the backbone of the first two units;
- ``head_end_radical`` / ``tail_end_radical``: a radical end at the head or
  the tail (the other end artificial);
- ``head_group_radical:<i>`` / ``tail_group_radical:<i>``: the declared end
  group attached at that end (the other end artificial) with a radical on
  its ``i``-th hydrogen-bearing heavy atom (none for a hydrogen cap);
- ``monomer``: the declared or derived monomer (no ports).

A structure that is an earlier one (same molecule, same port kinds; a chain
read backwards is the same structure) is dropped, so each class of
equivalent sites gives one structure.

Proxies: every structure but the monomer alone; every pair of single-radical
structures, itself included; every single-radical structure with
``pristine``; every single-radical structure with the monomer, at full
length (``<radical>+monomer``) and at ``units - 1``
(``<radical>_short+monomer``).  Each proxy is repeated at ``units + 2`` as
``<site type>@<units + 2>`` (the full contexts of the radius-one
completeness bound).  No proxy names family candidates: the compile tries
the requested families on all.

The pristine chain of every primitive-period count below ``units`` repeat
units is a proxy of its own too (``pristine~<k>``, k periods; no ``@``
copy).  A template's window closes on a side only where its witness ends
inside it, so an event whose fragments both end inside the window (PE: the
cut of n-butane into two ethyl radicals and its reverse) is witnessed only on
a chain that short.  Its ends are artificial like every other proxy's: an
event whose reacting atom is an end atom stays a cut point.

A schema-v2 (linear copolymer) declaration is refused by the functions
above; its proxies are ``sequence_proxy_set`` and its roster structures
``sequence_declared_structures`` (the section at the end of this module).
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import emit
from rmgpy.kmc import windows as win
from rmgpy.kmc.compiler import SiteProxy, construction_ports
from rmgpy.kmc.declared_chain import ChainSpec, EditOp
from rmgpy.kmc.proxy_padding import RepeatUnitGraph

DEFAULT_PROXY_UNITS = 3
# the extra units of the full-context copies
COVERAGE_EXTRA_UNITS = 2
SHORT = "_short"
# ``pristine~<k>``: the pristine chain of k primitive periods, below the proxies
SHORT_CHAIN = "~"
SINGLE_RADICALS = ("interior_radical", "pendant_radical", "head_end_radical",
                   "tail_end_radical", "head_group_radical", "tail_group_radical")
CONTEXT_CLASS = {"pristine": "pristine", "interior_radical": "featured",
                 "doubly_featured": "featured", "head_end_radical": "end-proximal",
                 "tail_end_radical": "end-proximal", "head_group_radical": "end-proximal",
                 "tail_group_radical": "end-proximal", "pendant_radical": "pendant",
                 "head_terminal_radical": "end-proximal", "short_radical": "end-proximal",
                 "monomer": "monomer"}


def radical_sites(declared) -> tuple:
    """``((kind, keys), ...)``: one entry per class of hydrogen-bearing heavy
    atoms of the declared primitive period (``ChainGeometry.atom_classes``:
    orbits under the locality policy's site classes); ``kind`` is
    ``interior`` for a backbone class and ``pendant`` otherwise; interior
    classes first, each kind in order of its least key."""
    g = declared.geometry
    bonded_h = {key for left, right in g.bonds for key, other in ((left, right), (right, left))
                if g.atoms[other].element == "H"}
    out = []
    for _site, positions in g.atom_classes:
        keys = tuple(sorted(g.keys[p] for p in positions))
        atom = g.atoms[keys[0]]
        if atom.element == "H" or not (atom.hydrogens or keys[0] in bonded_h):
            continue
        out.append(("interior" if keys[0] in g.backbone_offsets else "pendant", keys))
    return tuple(sorted(out, key=lambda item: (item[0] != "interior", item[1][0])))


def declared_proxy_set(spec, units: int = DEFAULT_PROXY_UNITS) -> tuple:
    """The compile proxies of a linear-homopolymer declaration at ``units``
    repeat units, with their ``units + 2`` full-context copies and the
    shorter pristine chains (module doc)."""
    builder = _Builder.of(spec, units)
    wide = units + COVERAGE_EXTRA_UNITS
    extended = tuple(replace(proxy, site_type=f"{proxy.site_type}@{wide}",
                             metadata={**proxy.metadata, "coverage_site_type": proxy.site_type})
                     for proxy in _proxies(builder, wide))
    short = builder.short_chains(units)
    return _proxies(builder, units) + extended + _proxies(builder, units, short,
                                                          [(name, (name,)) for name in short])


def declared_structures(spec, units: int = DEFAULT_PROXY_UNITS) -> tuple:
    """``(structures, proxies)`` of ``declared_proxy_set(spec, units)``:
    ``structures`` maps a participant name (``<name>@<depth>``) to
    ``(witness spec or None, build)`` (None: the portless monomer, an
    explicit build); ``proxies`` maps each site type to its participant
    names, in order."""
    builder = _Builder.of(spec, units)
    wide = units + COVERAGE_EXTRA_UNITS
    structures, proxies = {}, {}
    for depth, suffix in ((units, ""), (wide, f"@{wide}")):
        builds = builder.builds(depth)
        for name, item in builds.items():
            structures[f"{name}@{depth}"] = item
        for site_type, participants in _declarations(builds):
            proxies[site_type + suffix] = tuple(f"{name}@{depth}" for name in participants)
    for name, item in builder.short_chains(units).items():
        structures[name] = item
        proxies[name] = (name,)
    return structures, proxies


def _is_single_radical(name) -> bool:
    return not name.endswith(SHORT) and name.split(":")[0] in SINGLE_RADICALS


def _declarations(builds) -> list:
    """(site type, participant names) of every proxy over ``builds``."""
    radicals = [name for name in builds if _is_single_radical(name)]
    declarations = [(name, (name,)) for name in builds
                    if name != "monomer" and not name.endswith(SHORT)]
    for i, first in enumerate(radicals):
        for second in radicals[i:]:
            declarations.append((f"{first}+{second}", (first, second)))
    for radical in radicals:
        declarations.append((f"{radical}+pristine", (radical, "pristine")))
    if "monomer" in builds:
        for radical in radicals:
            declarations.append((f"{radical}+monomer", (radical, "monomer")))
            if radical + SHORT in builds:
                declarations.append((f"{radical}{SHORT}+monomer", (radical + SHORT, "monomer")))
    return declarations


def _proxies(builder, units: int, structures=None, declarations=None) -> tuple:
    """The proxies of ``declarations`` over ``structures`` (default: every
    structure at ``units`` and ``_declarations`` over them)."""
    from rmgpy.species import Species

    if structures is None:
        structures = builder.builds(units)
    builds = {name: build for name, (_spec, build) in structures.items()}
    species = {name: Species(molecule=[build.molecule]) for name, build in builds.items()}
    proxies = []
    for site_type, participants in (_declarations(builds) if declarations is None else declarations):
        boundaries = tuple(
            port
            for index, name in enumerate(participants)
            for port in construction_ports(name, builds[name].molecule, index,
                                           builds[name].metadata))
        site_types = tuple(name[:-len(SHORT)] if name.endswith(SHORT) else name.split(SHORT_CHAIN)[0]
                           for name in participants)
        proxies.append(SiteProxy(
            site_type,
            tuple(species[name] for name in participants),
            participant_site_types=site_types,
            artificial_boundaries=boundaries,
            repeat_unit=builder.repeat_unit if boundaries else None,
            boundary_builds=tuple(builds[name] for name in participants),
            metadata={"centred": True,
                      "context_class": CONTEXT_CLASS[site_types[0].split(":")[0]],
                      "proxy_units": units, "declaration": builder.declared.name},
        ))
    return tuple(proxies)


class _Builder:
    """The roster's builders for one declaration."""

    def __init__(self, spec):
        self.spec = spec
        self.declared = emit.declared_provider(spec)
        self.provider = emit.match_provider(spec.repeat_unit) or self.declared
        self.repeat_unit = self._repeat_unit()
        self.sites = radical_sites(self.declared)
        self._builds = {}

    @classmethod
    def of(cls, spec, units):
        from rmgpy.polymer_spec import normalized_linear_homopolymer_spec, refuse_schema_v2

        refuse_schema_v2(spec, "declared_proxies")
        if isinstance(units, bool) or not isinstance(units, int) or units < 2:
            raise ValueError(f"declared proxies need at least two units (a short variant has "
                             f"units - 1), got {units!r}")
        return cls(normalized_linear_homopolymer_spec(spec))

    def chain(self, units, head=None, tail=None, radical_unit=None):
        return self.provider.spec_for_declared_provider(self.declared, units, head, tail,
                                                        radical_unit)

    def builds(self, units: int) -> dict:
        """name -> (witness spec or None, build) of every structure at
        ``units``, with the ``_short`` single radicals at ``units - 1``."""
        if units not in self._builds:
            builds = self.structures(units)
            radicals = [name for name in builds if _is_single_radical(name)]
            if "monomer" in builds:
                short = self.structures(units - 1, only=radicals)
                builds.update({name + SHORT: item for name, item in short.items()})
            self._builds[units] = builds
        return self._builds[units]

    def short_chains(self, units: int) -> dict:
        """``pristine~<k>`` -> (witness spec, build) of the pristine chain of
        every k primitive periods below ``units`` repeat units (the chains
        the declared chain builds)."""
        out = {}
        for k in range(1, units * self.declared.periods_per_repeat):
            spec = win.DeclaredSpec(self.declared.name, ChainSpec(k))
            try:
                out[f"pristine{SHORT_CHAIN}{k}"] = (spec, win.build_spec(spec, "pristine"))
            except (b.BoundaryError, ValueError):
                continue
        return out

    def structures(self, units: int, only=None) -> dict:
        """name -> (witness spec, build) of every structure buildable at
        ``units`` (the monomer's spec is None)."""
        recipes = [("pristine", lambda: self.chain(units))]
        recipes += self._site_recipes(units, "interior")
        recipes.append(("doubly_featured", lambda: self._joined(
            self.chain(units, radical_unit=0), self.chain(units, radical_unit=1))))
        recipes += [("head_end_radical", lambda: self.chain(units, head="radical")),
                    ("tail_end_radical", lambda: self.chain(units, tail="radical"))]
        recipes += self._site_recipes(units, "pendant")
        recipes += self._group_recipes(units)
        out = {}
        for name, recipe in recipes:
            if only is not None and name not in only:
                continue
            try:
                spec = recipe()
                build = win.build_spec(spec, name)
            except (b.BoundaryError, ValueError):
                continue
            if not any(_same_structure(build, other) for _spec, other in out.values()):
                out[name] = (spec, build)
        if only is None:
            monomer = self._monomer()
            if monomer is not None:
                out["monomer"] = (None, monomer)
        return out

    def _site_recipes(self, units, kind) -> list:
        """One radical per class of ``kind`` (``radical_sites``): a backbone
        class on the centre unit, a side-group class on the first."""
        unit = emit.centred_radical_unit(units) * self.declared.periods_per_repeat \
            if kind == "interior" else 0
        classes = [keys for k, keys in self.sites if k == kind]
        return [(f"{kind}_radical:{i}",
                 lambda keys=keys: self._edited(self.chain(units), unit, keys[0]))
                for i, keys in enumerate(classes)]

    def _group_recipes(self, units) -> list:
        """One radical per hydrogen-bearing heavy atom of each declared end
        group, the group attached at its end (equivalent atoms give one
        structure: an equal one is dropped)."""
        out = []
        for side, graph in (("head", self.spec.head_end_group),
                            ("tail", self.spec.tail_end_group)):
            hydrogens = {atom.key for atom in graph.atoms if atom.element == "H"}
            bonded_h = {key for bond in graph.bonds for key in bond.atoms
                        if key not in hydrogens and set(bond.atoms) & hydrogens}
            keys = sorted(atom.key for atom in graph.atoms if atom.element != "H"
                          and (atom.hydrogens or atom.key in bonded_h))
            for i, key in enumerate(keys):
                def recipe(side=side, key=key):
                    spec = self.chain(units, **{side: "declared"})
                    return self._edited(spec, -1 if side == "head" else spec.chain.units, key)
                out.append((f"{side}_group_radical:{i}", recipe))
        return out

    @staticmethod
    def _edited(spec, unit, key):
        """``spec`` with one radical on (``unit``, ``key``)."""
        return replace(spec, chain=replace(spec.chain, edits=spec.chain.edits + (
            EditOp("remove_h", unit, key), EditOp("set_radical", unit, key, 1))))

    @staticmethod
    def _joined(first, second):
        """One declared spec carrying the edits of both (same chain otherwise)."""
        if replace(first.chain, edits=()) != replace(second.chain, edits=()):
            raise ValueError("the two radical specs are not one chain")
        return replace(first, chain=replace(first.chain,
                                            edits=first.chain.edits + second.chain.edits))

    def _monomer(self):
        molecule = self.provider.reference_monomer(self.spec)
        if molecule is None:
            return None
        if self.spec.monomer_product is not None:
            return self.declared.monomer("monomer")
        return b.explicit_species_build(
            "monomer", adjacency_list=molecule.to_adjacency_list(remove_h=False))

    def _repeat_unit(self) -> RepeatUnitGraph:
        """One primitive period of the declaration, its ports capped with
        hydrogen: the unit the padding extends by (``check_padding_declarations``)."""
        build = self.declared.build("repeat_unit", ChainSpec(1))
        heavy = [atom.id for atom in build.molecule.atoms if atom.element.number != 1]
        ports = {port.tag: port.atom_id for port in build.metadata.ports}
        return RepeatUnitGraph(self.declared.name, build.molecule,
                               heavy.index(ports["backbone_head"]),
                               heavy.index(ports["backbone_tail"]))


def _same_structure(first, second) -> bool:
    """The same molecule with the same port kinds, port atom onto port atom
    of the same kind (a chain read backwards is the same structure)."""
    kinds = lambda build: sorted((p.artificial, len(p.cap_atom_ids))
                                 for p in build.metadata.ports)
    if kinds(first) != kinds(second):
        return False
    a, c = first.molecule, second.molecule
    if not a.is_isomorphic(c, strict=False):
        return False
    a_by_id = {atom.id: atom for atom in a.atoms}
    c_by_id = {atom.id: atom for atom in c.atoms}
    by_kind = {}
    for port in second.metadata.ports:
        by_kind.setdefault((port.artificial, len(port.cap_atom_ids)), []).append(port.atom_id)
    pins = {}
    for port in first.metadata.ports:
        candidates = by_kind[(port.artificial, len(port.cap_atom_ids))]
        if len(candidates) == 1:
            pins[a_by_id[port.atom_id]] = c_by_id[candidates[0]]
    # whether one pinned isomorphism exists: the first decides
    return a.find_first_isomorphism(c, initial_map=pins, strict=False) is not None if pins else True


# ---------------------------------------------------------------------------
# linear copolymer declarations (schema v2)
# ---------------------------------------------------------------------------
# The compile structures of a ``PolymerSpec``, built by its residue-sequence
# chain provider (``declared_sequence``) on the witness words of its
# environment set (``declared_environments.environment_witnesses``):
#
# - ``interior_radical:<i>``: one per class of hydrogen-bearing heavy atoms
#   (``i`` counts the classes in key order), on the class's shortest witness
#   word with one more residue of the covering period on each side (the key
#   reads the heavy neighbours of its boundary shell, so a cap must not sit
#   on the shell), widened from the period until it holds ``units`` residues
#   with the radical's residue in the centre.  A class no single-residue
#   chain shows is a junction class (``SequenceStructure.junction``); those
#   come out of the words, not out of a rule of their own;
# - ``pristine:<i>``: the distinct chains of those radicals, unedited;
# - ``head_end_radical:<i>`` / ``tail_end_radical:<i>``: one per environment
#   of the radical end atom over every k-residue word at that end, on the
#   least such word (the other end artificial).  A shorter word can show the
#   same key only through its cap (a cap and a side group standing in for
#   the next residue), which any growth of the chain undoes, so the
#   witnesses are k-residue words;
# - ``head_group_radical:<i>`` / ``tail_group_radical:<i>``: likewise, per
#   (hydrogen-bearing heavy atom of the declared end group x environment);
# - ``monomer:<residue>``: one per ``monomer_products`` entry (no ports).
#
# Each structure is repeated at ``units + 2`` with one more residue at each
# artificial end (the full-context copies: the same names, no new ones).  Proxies: every structure but
# the monomers alone, and the pairs over family-tree classes
# (``pair_classes.class_pairs``) of each depth, computed from the loaded
# families by the compile and recorded with the artifact; the roster takes
# the recorded pairs (``sequence_declared_structures``).


@dataclass(frozen=True)
class SequenceStructure:
    """One compile structure of a v2 declaration: its witness spec
    (``windows.DeclaredSpec`` of the sequence provider; None for a monomer),
    its build, the environment key of the class it stands for (radicals) and
    the witness word (interior radicals)."""

    witness: object
    build: object
    environment: object = None
    junction: bool = False
    word: tuple = ()


PAIR_RULE = "family_tree_classes"


def _h_bearing(graph, key) -> bool:
    atoms = {atom.key: atom for atom in graph.atoms}
    hydrogens = {k for k, atom in atoms.items() if atom.element == "H"}
    bonded = {k for bond in graph.bonds for k in bond.atoms if set(bond.atoms) & hydrogens}
    return atoms[key].element != "H" and bool(atoms[key].hydrogens or key in bonded - hydrogens)


class _SequenceBuilder:
    """The structures of one v2 declaration (comment above)."""

    def __init__(self, spec):
        from rmgpy.kmc import declared_environments as de
        from rmgpy.kmc.declared_chain import site_class_policy
        from rmgpy.kmc.declared_sequence import sequence_provider

        self.spec = spec
        self.provider = sequence_provider(spec)
        self.period, self.witnesses = de.environment_witnesses(spec)
        self.k = de.spec_word_order(spec)
        self.policy = site_class_policy(spec.locality_policy, de.MAX_HEAVY_ATOMS)
        self.homopolymer = set()
        for name, _residue in spec.residues:
            self.homopolymer |= set(de.period_environments(
                dict(spec.residues), spec.junctions, spec.locality_policy, (name,)))
        self._built = {}

    # -- chains ----------------------------------------------------------------
    def _from_period(self, word, position, start, units):
        """``word`` (at period index ``start``) with one more residue of the
        covering period on each side, widened from the period until it has
        ``units`` residues and ``position`` is central."""
        size = len(self.period)
        left, right = start - 1, start + len(word)
        residues = [self.period[left % size]] + list(word) + [self.period[right % size]]
        position += 1
        while True:
            before, after = position, len(residues) - 1 - position
            if len(residues) >= units and abs(before - after) <= 1:
                return tuple(residues), position
            if before < after or (before == after and len(residues) % 2 == 0):
                left -= 1
                residues.insert(0, self.period[left % size])
                position += 1
            else:
                right += 1
                residues.append(self.period[right % size])

    def _support_paths(self, length):
        paths = [(name,) for name in sorted(self.provider.support)]
        for _ in range(length - 1):
            paths = [p + (b,) for p in paths for b in self.provider.support[p[-1]]]
        return paths

    def _padded(self, residues, side, units):
        """``residues`` grown at the end opposite ``side`` (the artificial
        one) along the support until it has ``units`` residues."""
        residues = list(residues)
        support = self.provider.support
        while len(residues) < units:
            if side == "head":
                residues.append(min(support[residues[-1]]))
            else:
                residues.insert(0, min(a for a in sorted(support) if residues[0] in support[a]))
        return tuple(residues)

    def _key(self, build, unit, key):
        from rmgpy.kmc.site_class import extract_site_class

        if unit in (-1, len(build.metadata.residues)):
            graph = self.spec.head_end_group if unit == -1 else self.spec.tail_end_group
            position = sorted(atom.key for atom in graph.atoms).index(key)
        else:
            position = self.provider.position(build.metadata.residues[unit], key)
        atom_id = next(a for a, u, p in build.metadata.coordinates if (u, p) == (unit, position))
        atom = next(atom for atom in build.molecule.atoms if atom.id == atom_id)
        return extract_site_class(build.molecule, [("site", atom)], self.policy)[0]

    @staticmethod
    def _radical(chain, unit, key):
        from rmgpy.kmc.declared_chain import EditOp
        return replace(chain, edits=chain.edits + (EditOp("remove_h", unit, key),
                                                   EditOp("set_radical", unit, key, 1)))

    def _interior(self, units) -> list:
        from rmgpy.kmc.declared_sequence import SequenceChainSpec

        residues = dict(self.spec.residues)
        out, index = [], 0
        for key, occurrences in self.witnesses.items():
            word, position, atom, start = occurrences[0]
            if not _h_bearing(residues[word[position]].graph, atom):
                continue
            chain, unit = self._from_period(word, position, start, units)
            out.append((f"interior_radical:{index}",
                        self._radical(SequenceChainSpec(chain), unit, atom),
                        dict(environment=key, junction=key not in self.homopolymer, word=word)))
            index += 1
        return out

    def _ends(self, units) -> list:
        """End radicals and end-group radicals (comment above)."""
        from rmgpy.kmc.declared_sequence import SequenceChainSpec

        out = []
        for side in ("head", "tail"):
            group = self.spec.head_end_group if side == "head" else self.spec.tail_end_group
            sites = [("end_radical", None)] + [
                ("group_radical", key) for key in sorted(atom.key for atom in group.atoms)
                if _h_bearing(group, key)]
            for kind, key in sites:
                found = {}
                for path in self._support_paths(self.k):
                    found.setdefault(self._key(*self._end_probe(path, side, key)), []).append(path)
                index = len([n for n, _c, _i in out if n.startswith(f"{side}_{kind}:")])
                for environment in sorted(found):
                    path = min(found[environment])
                    chain = self._end_chain(SequenceChainSpec(self._padded(path, side, units)),
                                            side, key)
                    out.append((f"{side}_{kind}:{index}", chain, dict(environment=environment)))
                    index += 1
        return out

    def _end_probe(self, path, side, key):
        from rmgpy.kmc.declared_sequence import SequenceChainSpec
        chain = self._end_chain(SequenceChainSpec(path), side, key)
        unit, atom = self._end_site(chain, side, key)
        return self.provider.build("end", chain), unit, atom

    def _end_site(self, chain, side, key):
        """(unit, key) of the radical of an end structure."""
        if key is None:
            residue = chain.residues[0] if side == "head" else chain.residues[-1]
            return ((0, self.provider.entry[residue]) if side == "head"
                    else (len(chain.residues) - 1, self.provider.exit[residue]))
        return (-1 if side == "head" else len(chain.residues)), key

    def _end_chain(self, chain, side, key):
        if key is None:
            return replace(chain, **{side: "radical"})
        chain = replace(chain, **{side: "declared"})
        unit, atom = self._end_site(chain, side, key)
        return self._radical(chain, unit, atom)

    def _wide(self, chain):
        """``chain`` with one more residue at each artificial end, from the
        support (every edit off the head group moves with its unit)."""
        from rmgpy.kmc.declared_sequence import SequenceChainSpec

        residues = list(chain.residues)
        support = self.provider.support
        edits = chain.edits
        if chain.head == "artificial":
            residues.insert(0, min(a for a in sorted(support) if residues[0] in support[a]))
            edits = tuple(replace(e, unit=e.unit + 1) if e.unit >= 0 else e for e in edits)
        if chain.tail == "artificial":
            residues.append(min(support[residues[-1]]))
        return SequenceChainSpec(tuple(residues), chain.head, chain.tail, edits)

    def _monomers(self) -> list:
        out = []
        for name, graph in self.spec.monomer_products:
            molecule = graph.to_legacy_molecule({})
            molecule.multiplicity = graph.multiplicity
            out.append((f"monomer:{name}", b.explicit_species_build(
                f"monomer:{name}", adjacency_list=molecule.to_adjacency_list(remove_h=False))))
        return out

    def structures(self, units: int, wide: bool = False) -> dict:
        """name -> ``SequenceStructure`` at ``units`` (``wide``: the
        full-context copy of each of those structures, under its name)."""
        if (units, wide) in self._built:
            return self._built[(units, wide)]
        if wide:
            out = {}
            for name, s in self.structures(units).items():
                if s.witness is None:
                    out[name] = s
                    continue
                chain = self._wide(s.witness.chain)
                out[name] = replace(s, witness=win.DeclaredSpec(self.provider.name, chain),
                                    build=self.provider.build(name, chain))
            self._built[(units, wide)] = out
            return out
        radicals = self._interior(units)
        recipes = radicals + [("pristine", replace(chain, edits=()), {})
                              for _name, chain, _info in radicals] + self._ends(units)
        out = {}
        for name, chain, info in recipes:
            build = self.provider.build(name, chain)
            if any(_same_structure(build, other.build) for other in out.values()):
                if name.startswith("interior_radical:"):
                    raise ValueError(f"{name}: two environment classes gave one structure")
                continue
            if name == "pristine":
                name = f"pristine:{sum(n.startswith('pristine:') for n in out)}"
            out[name] = SequenceStructure(win.DeclaredSpec(self.provider.name, chain), build,
                                          **info)
        for name, build in self._monomers():
            out[name] = SequenceStructure(None, build)
        self._built[(units, wide)] = out
        return out


# The compile structures of a ``BranchedPolymerSpec`` are those of its lift
# (``declared_branches``: a linear alphabet of residue variants, built by the
# lift's chain provider under ``branched:<sha16>``), with its ends read from
# the end and short words of its environment set:
#
# - ``interior_radical:<i>`` / ``pristine:<i>``: as above, on the interior
#   witnesses of the lift (stub atoms are never witnesses);
# - ``head_end_radical:<i>`` / ``tail_end_radical:<i>``: a radical at a
#   forbidden-vacant end site (a cut link: no H cap can stand there), one per
#   environment over the support words of rho + 2 variants from that end;
# - ``head_terminal_radical:<i>``: one per environment of a hydrogen-bearing
#   atom first seen near a real end (an H-capped vacant entry site) and not on
#   an interior radical, on its end word with one more variant inward;
# - ``short_radical:<i>``: likewise for a molecule both of whose ends the
#   environment reads (a word of at most k variants, both ends capped).
# A radical on a vacant end site's own atom takes the cap's place
# (``head``/``tail`` ``radical``); any other is an H removed.


class _BranchedBuilder(_SequenceBuilder):
    """The structures of one branched declaration (comment above)."""

    def __init__(self, spec):
        from rmgpy.kmc.declared_branches import branched_lift, branched_provider
        from rmgpy.kmc.declared_chain import site_class_policy
        from rmgpy.kmc import declared_environments as de

        self.declaration = spec
        self.lift = branched_lift(spec)
        self.spec = self.lift
        self.provider = branched_provider(spec)
        self.period, self.witnesses = self.lift.environment_witnesses()
        self.k = self.lift.k
        self.policy = site_class_policy(spec.locality_policy, de.MAX_HEAVY_ATOMS)
        self.homopolymer = set()
        self._built = {}

    def _checked(self, name, chain, environment, unit, key):
        """``chain`` after checking that its pristine atom ``(unit, key)``
        reads ``environment`` in the lift's own key space (the declared
        bonds; a built chain is Kekule, whose keys differ on an aromatic
        residue), so its word, caps and padding are the environment's."""
        found = self.lift.word_keys(chain.residues, unit, True, chain.tail != "artificial")[key]
        if found != environment:
            raise ValueError(f"{name}: its pristine atom ({unit}, {key!r}) reads {found!r}, not the "
                             f"environment {environment!r} it stands for")
        return chain

    def _terminal(self, name, word, position, key, both):
        from rmgpy.kmc.declared_sequence import SequenceChainSpec

        lift = self.lift
        residues = tuple(word)
        if not both:
            residues = residues + (min(lift.support[residues[-1]]),)
        chain = SequenceChainSpec(residues, "declared", "declared" if both else "artificial")
        atom = next(a for a in lift.residues[residues[position]].graph.atoms if a.key == key)
        if not atom.hydrogens and position == 0 and key == lift.entry[residues[0]]:
            return replace(chain, head="radical")
        if not atom.hydrogens and both and position == len(residues) - 1 and key == lift.exit[residues[-1]]:
            return replace(chain, tail="radical")
        return self._radical(chain, position, key)

    def _ends(self, units) -> list:
        """Open-end, terminal and short-molecule radicals (comment above)."""
        from rmgpy.kmc.declared_sequence import SequenceChainSpec

        lift, out = self.lift, []
        for side in ("head", "tail"):
            found = {}
            for path in self._support_paths(lift.rho + 2):
                end = path[0] if side == "head" else path[-1]
                if (lift.vacant_entry(end) if side == "head" else lift.vacant_exit(end)):
                    continue
                found.setdefault(self._key(*self._end_probe(path, side, None)), []).append(path)
            for index, environment in enumerate(sorted(found)):
                path = min(found[environment])
                chain = self._end_chain(SequenceChainSpec(self._padded(path, side, units)), side, None)
                out.append((f"{side}_end_radical:{index}", chain, dict(environment=environment)))
        covered = {key for key, occurrences in self.witnesses.items()
                   if _h_bearing(lift.residues[occurrences[0][0][occurrences[0][1]]].graph, occurrences[0][2])}
        environments = lift.environments()
        counters = {"end": 0, "short": 0}
        for environment, occurrences in environments.items():
            if environment in covered:
                continue
            for kind, word, position, _start, key in occurrences:
                if kind == "interior" or not lift.h_bearing(kind, word, position, key):
                    continue
                prefix = "head_terminal_radical" if kind == "end" else "short_radical"
                name = f"{prefix}:{counters[kind]}"
                chain = self._terminal(name, word, position, key, kind == "short")
                out.append((name, self._checked(name, chain, environment, position, key),
                            dict(environment=environment)))
                counters[kind] += 1
                break
        return out


def _sequence_builder(spec, units) -> _SequenceBuilder:
    from rmgpy.polymer_spec import BranchedPolymerSpec, PolymerSpec

    if not isinstance(spec, (PolymerSpec, BranchedPolymerSpec)):
        raise TypeError(f"a sequence proxy set takes a PolymerSpec or a BranchedPolymerSpec (schema "
                        f"v2), got {type(spec).__name__}")
    if isinstance(units, bool) or not isinstance(units, int) or units < 1:
        raise ValueError(f"sequence proxies need at least one residue per chain, got {units!r}")
    return _BranchedBuilder(spec) if isinstance(spec, BranchedPolymerSpec) else _SequenceBuilder(spec)


def sequence_structures(spec, units: int = DEFAULT_PROXY_UNITS, wide: bool = False) -> dict:
    """name -> ``SequenceStructure`` of a ``PolymerSpec`` at ``units``
    residues (``wide``: their ``units + 2`` full-context copies)."""
    return _sequence_builder(spec, units).structures(units, wide)


def _sequence_declarations(structures, pairs) -> list:
    """(site type, participant names) of every proxy over ``structures``."""
    declarations = [(name, (name,)) for name in structures if not name.startswith("monomer")]
    for first, second in pairs:
        if first not in structures or second not in structures:
            raise ValueError(f"the pair {first}+{second} names a structure that is not built")
        declarations.append((f"{first}+{second}", (first, second)))
    return declarations


def _depths(units) -> tuple:
    return ((units, False, ""), (units + COVERAGE_EXTRA_UNITS, True,
                                 f"@{units + COVERAGE_EXTRA_UNITS}"))


def _check_pairs(pairs, units) -> dict:
    expected = {str(depth) for depth, _wide, _suffix in _depths(units)}
    if not isinstance(pairs, dict) or set(pairs) != expected:
        raise ValueError(f"the recorded pairs must give the pairs of depths {sorted(expected)}")
    return {depth: [tuple(pair) for pair in value] for depth, value in pairs.items()}


def sequence_proxy_set(spec, units: int, families) -> tuple:
    """``(proxies, record)`` of a ``PolymerSpec``: every structure but the
    monomers alone and the pairs over family-tree classes of ``families``
    (name -> loaded family) at ``units`` and at its full-context depth;
    ``record`` holds the pairs per depth (``rule``: ``family_tree_classes``)
    for the artifact and its roster."""
    from rmgpy.kmc.pair_classes import class_pairs

    builder = _sequence_builder(spec, units)
    pairs, classes = {}, {}
    for depth, wide, _suffix in _depths(units):
        structures = builder.structures(units, wide)
        molecules = {name: s.build.molecule for name, s in structures.items()}
        found, members = class_pairs(families, molecules)
        pairs[str(depth)] = [list(pair) for pair in found]
        classes[str(depth)] = [list(group) for group in members]
    record = {"rule": PAIR_RULE, "families": sorted(families), "pairs": pairs, "classes": classes}
    return _sequence_proxies(builder, units, pairs), record


def _sequence_proxies(builder, units, pairs) -> tuple:
    from rmgpy.species import Species

    pairs = _check_pairs(pairs, units)
    proxies = []
    for depth, wide, suffix in _depths(units):
        structures = builder.structures(units, wide)
        builds = {name: s.build for name, s in structures.items()}
        species = {name: Species(molecule=[build.molecule]) for name, build in builds.items()}
        for site_type, participants in _sequence_declarations(structures, pairs[str(depth)]):
            boundaries = tuple(
                port for index, name in enumerate(participants)
                for port in construction_ports(name, builds[name].molecule, index,
                                               builds[name].metadata))
            metadata = {"centred": True,
                        "context_class": CONTEXT_CLASS[participants[0].split(":")[0]],
                        "proxy_units": depth, "declaration": builder.provider.name}
            if wide:
                metadata["coverage_site_type"] = site_type
            proxies.append(SiteProxy(
                site_type + suffix,
                tuple(species[name] for name in participants),
                participant_site_types=tuple(participants),
                artificial_boundaries=boundaries,
                repeat_unit=None,
                boundary_builds=tuple(builds[name] for name in participants),
                metadata=metadata,
            ))
    return tuple(proxies)


def sequence_declared_structures(spec, units: int, pairs) -> tuple:
    """``(structures, proxies)`` of the compile of a ``PolymerSpec`` at
    ``units`` with the recorded ``pairs``: ``structures`` maps a participant
    name (``<name>@<depth>``) to ``(witness spec or None, build)`` (None: a
    monomer, an explicit build); ``proxies`` maps each site type to its
    participant names, in order (``declared_structures`` of a v1
    declaration)."""
    builder = _sequence_builder(spec, units)
    pairs = _check_pairs(pairs, units)
    structures, proxies = {}, {}
    for depth, wide, suffix in _depths(units):
        found = builder.structures(units, wide)
        for name, s in found.items():
            structures[f"{name}@{depth}"] = (s.witness, s.build)
        for site_type, participants in _sequence_declarations(found, pairs[str(depth)]):
            proxies[site_type + suffix] = tuple(f"{name}@{depth}" for name in participants)
    return structures, proxies
