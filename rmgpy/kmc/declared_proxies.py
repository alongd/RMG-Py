"""The compile proxies of a declared linear homopolymer.

``declared_proxy_set(spec, units)`` builds the centred finite structures the
compiler generates family reactions on from the declaration alone, with the
builders the emit roster uses (``emit.reference_systems``): the declaration's
own chain provider (``emit.declared_provider``) through the vocabulary
adapter of the provider whose repeat unit it declares (``emit.match_provider``,
``spec_for_declared_provider``), and the monomer the declaration names or
derives (``reference_monomer``).  Nothing here names a polymer.

The structures, at ``units`` repeat units with artificial (capped) ends:

- ``pristine``: the closed chain;
- ``interior_radical``: one radical on the centre unit's backbone;
- ``doubly_featured``: radicals on the backbone of the first two units;
- ``head_end_radical`` / ``tail_end_radical``: a radical end at the head or
  the tail (the other end artificial); the tail one is dropped when it is the
  head one read backwards (same molecule, same port kinds);
- ``pendant_radical``: one radical on the first unit's side-group atom
  farthest from the backbone (the least key when several are equally far;
  none for a chain without a hydrogen-bearing side-group atom);
- ``monomer``: the declared or derived monomer (no ports).

Proxies: every structure but the monomer alone; every pair of single-radical
structures (``interior_radical``, the end radicals, ``pendant_radical``),
itself included; every single-radical structure with ``pristine``; every
single-radical structure at ``units - 1`` (its ``_short`` variant) with the
monomer.  Each proxy is repeated at ``units + 2`` as ``<site type>@<units +
2>`` (the full contexts of the radius-one completeness bound).  No proxy
names family candidates: the compile tries the requested families on all.
"""

from __future__ import annotations

from collections import deque
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
SINGLE_RADICALS = ("interior_radical", "head_end_radical", "tail_end_radical",
                   "pendant_radical")
CONTEXT_CLASS = {"pristine": "pristine", "interior_radical": "featured",
                 "doubly_featured": "featured", "head_end_radical": "end-proximal",
                 "tail_end_radical": "end-proximal", "pendant_radical": "pendant",
                 "monomer": "monomer"}


def declared_proxy_set(spec, units: int = DEFAULT_PROXY_UNITS) -> tuple:
    """The compile proxies of a linear-homopolymer declaration at ``units``
    repeat units, with their ``units + 2`` full-context copies (module doc)."""
    from rmgpy.polymer_spec import normalized_linear_homopolymer_spec

    if isinstance(units, bool) or not isinstance(units, int) or units < 2:
        raise ValueError(f"declared proxies need at least two units (a short variant has "
                         f"units - 1), got {units!r}")
    spec = normalized_linear_homopolymer_spec(spec)
    wide = units + COVERAGE_EXTRA_UNITS
    extended = tuple(replace(proxy, site_type=f"{proxy.site_type}@{wide}",
                             metadata={**proxy.metadata, "coverage_site_type": proxy.site_type})
                     for proxy in _proxies(spec, wide))
    return _proxies(spec, units) + extended


def _proxies(spec, units: int) -> tuple:
    from rmgpy.species import Species

    builder = _Builder(spec)
    builds = builder.structures(units)
    radicals = [name for name in SINGLE_RADICALS if name in builds]
    short = builder.structures(units - 1, only=radicals) if "monomer" in builds else {}
    builds.update({f"{name}_short": build for name, build in short.items()})
    species = {name: Species(molecule=[build.molecule]) for name, build in builds.items()}

    declarations = [(name, (name,)) for name in builds
                    if name != "monomer" and not name.endswith("_short")]
    for i, first in enumerate(radicals):
        for second in radicals[i:]:
            declarations.append((f"{first}+{second}", (first, second)))
    for radical in radicals:
        declarations.append((f"{radical}+pristine", (radical, "pristine")))
    for radical in radicals:
        if f"{radical}_short" in builds:
            declarations.append((f"{radical}+monomer", (f"{radical}_short", "monomer")))

    proxies = []
    for site_type, participants in declarations:
        boundaries = tuple(
            port
            for index, name in enumerate(participants)
            for port in construction_ports(name, builds[name].molecule, index,
                                           builds[name].metadata))
        site_types = tuple(name[:-len("_short")] if name.endswith("_short") else name
                           for name in participants)
        proxies.append(SiteProxy(
            site_type,
            tuple(species[name] for name in participants),
            participant_site_types=site_types,
            artificial_boundaries=boundaries,
            repeat_unit=builder.repeat_unit if boundaries else None,
            boundary_builds=tuple(builds[name] for name in participants),
            metadata={"centred": True, "context_class": CONTEXT_CLASS[site_types[0]],
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

    def chain(self, units, head=None, tail=None, radical_unit=None):
        return self.provider.spec_for_declared_provider(self.declared, units, head, tail,
                                                        radical_unit)

    @staticmethod
    def build(spec, name):
        return win.build_spec(spec, name)

    def structures(self, units: int, only=None) -> dict:
        """name -> ``ProxyBuild`` of every structure buildable at ``units``."""
        recipes = {
            "pristine": lambda: [self.chain(units)],
            "interior_radical": lambda: [self.chain(units, radical_unit=units // 2)],
            "doubly_featured": lambda: [self._joined(self.chain(units, radical_unit=0),
                                                     self.chain(units, radical_unit=1))]
            if units > 1 else [],
            "head_end_radical": lambda: [self.chain(units, head="radical")],
            "tail_end_radical": lambda: [self.chain(units, tail="radical")],
            "pendant_radical": lambda: self._pendant(units),
        }
        out = {}
        for name, recipe in recipes.items():
            if only is not None and name not in only:
                continue
            try:
                specs = recipe()
            except (b.BoundaryError, ValueError):
                continue
            for spec in specs:
                build = self.build(spec, name)
                if not any(_same_structure(build, other) for other in out.values()):
                    out[name] = build
                    break
        if only is None:
            monomer = self._monomer()
            if monomer is not None:
                out["monomer"] = monomer
        return out

    @staticmethod
    def _joined(first, second):
        """One declared spec carrying the edits of both (same chain otherwise)."""
        if replace(first.chain, edits=()) != replace(second.chain, edits=()):
            raise ValueError("the two radical specs are not one chain")
        return replace(first, chain=replace(first.chain,
                                            edits=first.chain.edits + second.chain.edits))

    def _pendant(self, units) -> list:
        """The pristine chain with one radical on a side-group atom of its first
        unit, per side-group atom farthest from the backbone."""
        g = self.declared.geometry
        backbone = set(g.backbone_offsets)
        neighbours = {key: set() for key in g.keys}
        for left, right in g.bonds:
            neighbours[left].add(right)
            neighbours[right].add(left)
        distance = {key: 0 for key in backbone}
        queue = deque(sorted(backbone))
        while queue:
            key = queue.popleft()
            for other in sorted(neighbours[key]):
                if other not in distance:
                    distance[other] = distance[key] + 1
                    queue.append(other)
        sites = [key for key in g.keys if key not in backbone and key in distance
                 and g.atoms[key].element != "H" and g.atoms[key].hydrogens]
        if not sites:
            return []
        far = max(distance[key] for key in sites)
        pristine = self.chain(units)
        return [replace(pristine, chain=replace(pristine.chain, edits=pristine.chain.edits + (
            EditOp("remove_h", 0, key), EditOp("set_radical", 0, key, 1))))
            for key in sorted(sites) if distance[key] == far]

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
