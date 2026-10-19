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
"""

from __future__ import annotations

from dataclasses import replace

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
SINGLE_RADICALS = ("interior_radical", "pendant_radical", "head_end_radical",
                   "tail_end_radical", "head_group_radical", "tail_group_radical")
CONTEXT_CLASS = {"pristine": "pristine", "interior_radical": "featured",
                 "doubly_featured": "featured", "head_end_radical": "end-proximal",
                 "tail_end_radical": "end-proximal", "head_group_radical": "end-proximal",
                 "tail_group_radical": "end-proximal", "pendant_radical": "pendant",
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
    repeat units, with their ``units + 2`` full-context copies (module doc)."""
    builder = _Builder.of(spec, units)
    wide = units + COVERAGE_EXTRA_UNITS
    extended = tuple(replace(proxy, site_type=f"{proxy.site_type}@{wide}",
                             metadata={**proxy.metadata, "coverage_site_type": proxy.site_type})
                     for proxy in _proxies(builder, wide))
    return _proxies(builder, units) + extended


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


def _proxies(builder, units: int) -> tuple:
    from rmgpy.species import Species

    builds = {name: build for name, (_spec, build) in builder.builds(units).items()}
    species = {name: Species(molecule=[build.molecule]) for name, build in builds.items()}
    proxies = []
    for site_type, participants in _declarations(builds):
        boundaries = tuple(
            port
            for index, name in enumerate(participants)
            for port in construction_ports(name, builds[name].molecule, index,
                                           builds[name].metadata))
        site_types = tuple(name[:-len(SHORT)] if name.endswith(SHORT) else name
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
    return bool(a.find_isomorphism(c, initial_map=pins, strict=False)) if pins else True
