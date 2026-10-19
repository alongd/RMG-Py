"""Template catalog emission (I-155): a compiled event-set artifact plus a roster
of declared reference systems becomes one versioned, provenance-stamped,
gate-checked ``kmc-template-catalog/1`` file the melt engine can load.

Two steps, both from the command line::

    python -m rmgpy.kmc.emit ARTIFACT ROSTER --database DB --out rows.jsonl
        [--chunk i/n] [--family F]... [--limit N] [--no-contexts] [--legacy-artifact]
    python -m rmgpy.kmc.emit assemble ROWS... --out catalog.json [--database DB]
    python -m rmgpy.kmc.emit catalog-set CATALOG... --out SET.json

A one-family artifact (``compile --family F``) emits and assembles into a
complete catalog; ``--family`` naming every family of the artifact is the same
full selection.  ``catalog-set`` ties per-family catalogs into one
``kmc-catalog-set/1`` manifest (``rmgpy.kmc.catalog_set``) that the runner
accepts wherever it accepts catalogs.

**Emit.** For each record (or the shard ``--chunk i/n``: records ``i::n`` of the
event-id-sorted list): certify the records, match the reference cases of the
roster's systems, ``template_emission.emit_record`` with the production rates
and the reference Kc.  One JSON line per record: the payload plus status,
refusal and timing.  An error is a row (``status: error``), never silent.  The
first line of a rows file is its header (artifact, database, code, roster and
shard identity).

**Roster.** A JSON file ``{"schema": "kmc-reference-roster/2", "polymer":
<LinearHomopolymerSpec mapping>, "systems": {name: [participant, ...]}}``.  A
participant is ``{"chain": {"owner", "dp", "head", "tail", "radical_unit"}}``
(a chain of the declared repeat unit; ``owner`` names the exact registered provider whose
``END_CHEMISTRIES`` vocabulary owns ``head``/``tail``, ``null`` the artificial H-capped end),
``{"smiles": ...}`` or ``{"adjacency_list": ...}``.  The declared repeat unit is
rebuilt by a declaration provider registered from the embedded polymer; each chain's vocabulary
owner must be registered and declare that repeat unit, or it is refused (``refused:owner``).
A ``kmc-reference-roster/3`` roster holds a declared artifact's compile structures
(``{"declared_chain": ChainSpec}`` participants); a ``kmc-reference-roster/4`` roster those
of an artifact compiled from a schema-v2 (copolymer) declaration, with the pairs its compile
recorded (``{"declared_sequence": SequenceChainSpec}`` participants).

**Binding.** An artifact compiled from a deck carries the declaration's
``chemistry_sha256`` and ``declaration_sha256`` (``rmgpy.kmc.run``) in
``provenance.proxy_source``; the roster's embedded declaration must have the
same ones, or the run is refused (``refused:declaration``) before anything is
built, whatever the flags.  An artifact carrying neither (a compile without a
deck) is refused (``refused:unbound_artifact``) unless ``--legacy-artifact``.

**Assemble.** Merge the shards, then ``resolve_collisions``, ``derive_weights``
and ``check_cycles`` over the whole set.  The catalog is written even when a gate
fails; ``complete`` is true only when every gate passes and no record errored.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import traceback
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from itertools import combinations
from pathlib import Path
from typing import Any, Callable, Optional

from rmgpy.kmc import boundaries as b
from rmgpy.kmc.catalog_set import CatalogSetError, build_catalog_set  # noqa: F401
from rmgpy.kmc import template_emission as te
from rmgpy.kmc import windows as win

SCHEMA = "kmc-template-catalog/1"
ROWS_SCHEMA = "kmc-emit-rows/1"
ROSTER_SCHEMA = "kmc-reference-roster/2"
# a declared artifact's roster: its compile structures (``declared_reference_systems``)
ROSTER_SCHEMA_V3 = "kmc-reference-roster/3"
# the roster of an artifact compiled from a schema-v2 (copolymer) declaration:
# its compile structures and the pairs the compile recorded
ROSTER_SCHEMA_V4 = "kmc-reference-roster/4"
# Kept equal to the acceptance harness' labels so that rows of the two agree
# (neither string selects data: the Kc label tags the reference thermo, the
# rate label tags the production-rate provenance)
REFERENCE_THERMO_LABEL = "i147-reference"
REFERENCE_RATES_LABEL = "i138-reference"
PAYLOAD_FIELDS = ("key", "name", "family", "status", "direction", "executable", "weight", "derived")
HEADER_MATCH = ("schema", "artifact_id", "artifact_sha256", "artifact_records", "families",
                "database_path", "database_head", "database_content_sha256", "code_sha",
                "code_sources_sha256", "code_extensions_sha256", "roster_sha256", "temperatures",
                "selection")
CYCLE_REFUSED_FRACTION = 0.05
MIDPOINT_THRESHOLD = 0.01     # the harness' gate on the primary-weight midpoint error


# ---------------------------------------------------------------------------
# reference systems from a declaration
# ---------------------------------------------------------------------------
class RosterRefusal(ValueError):
    """A roster participant that cannot be built (``status`` names why)."""

    def __init__(self, status, reason):
        super().__init__(f"{status}: {reason}")
        self.status = status
        self.reason = reason


@dataclass(frozen=True)
class Roster:
    polymer: Any
    systems: dict
    sha256: str
    path: str
    schema: str = ROSTER_SCHEMA
    # v3: the proxy units of the compile whose structures the roster holds
    proxy_units: Optional[int] = None
    # v4: the pair proxies the compile recorded, per depth
    pairs: Optional[dict] = None

    def digests(self) -> dict:
        """The digests of the declaration the roster carries (``run``'s
        ``chemistry_sha256`` and ``declaration_sha256``; a deck of the same
        declaration has the same ones)."""
        from rmgpy.kmc.run import chemistry_sha256, declaration_sha256
        return {"chemistry_sha256": chemistry_sha256(self.polymer),
                "declaration_sha256": declaration_sha256(self.polymer)}


@dataclass
class Participant:
    """One reference participant: the spec the windows machinery reads and the
    builder of its fresh proxy build (a build is made per use, as the
    acceptance harness does)."""

    label: str
    spec: Any
    make: Callable

    def build(self, name):
        return self.make(name)


ROSTER_FIELDS = {ROSTER_SCHEMA: {"schema", "polymer", "systems"},
                 ROSTER_SCHEMA_V3: {"schema", "polymer", "proxy_units", "systems"},
                 ROSTER_SCHEMA_V4: {"schema", "polymer", "proxy_units", "pairs", "systems"}}


def load_roster(path) -> Roster:
    """A roster file.  A v4 roster carries a schema-v2 declaration and every
    other schema a v1 one; the other combination is refused by name."""
    from rmgpy.polymer_spec import (
        is_schema_v2,
        normalized_linear_homopolymer_spec,
        normalized_polymer_spec,
        refuse_schema_v2,
    )
    raw = Path(path).read_bytes()
    data = json.loads(raw)
    schema = data.get("schema") if isinstance(data, dict) else None
    if schema not in ROSTER_FIELDS:
        raise RosterRefusal("refused:roster", f"roster schema must be one of "
                                              f"{sorted(ROSTER_FIELDS)}, got {schema!r}")
    if set(data) != ROSTER_FIELDS[schema]:
        raise RosterRefusal("refused:roster", f"a {schema} roster has exactly "
                                              f"{', '.join(sorted(ROSTER_FIELDS[schema]))}")
    units = data.get("proxy_units")
    if schema == ROSTER_SCHEMA_V3 and (isinstance(units, bool) or not isinstance(units, int)
                                       or units < 2):
        raise RosterRefusal("refused:roster", f"proxy_units must be an integer >= 2, got {units!r}")
    if schema == ROSTER_SCHEMA_V4:
        if isinstance(units, bool) or not isinstance(units, int) or units < 1:
            raise RosterRefusal("refused:roster", f"proxy_units must be an integer >= 1, got "
                                                  f"{units!r}")
        if not is_schema_v2(data["polymer"]):
            raise RosterRefusal("refused:roster", f"a {ROSTER_SCHEMA_V4} roster carries a "
                                                  "schema-v2 declaration; this one is schema v1")
        pairs = data["pairs"]
        if not isinstance(pairs, dict) or not all(
                isinstance(v, list) and all(isinstance(p, list) and len(p) == 2 for p in v)
                for v in pairs.values()):
            raise RosterRefusal("refused:roster", "pairs must map a depth to a list of "
                                                  "participant-name pairs")
        polymer = normalized_polymer_spec(data["polymer"])
    else:
        polymer = normalized_linear_homopolymer_spec(refuse_schema_v2(data["polymer"],
                                                                      "emit.load_roster"))
    systems = data["systems"]
    if not isinstance(systems, dict) or not all(isinstance(v, list) for v in systems.values()):
        raise RosterRefusal("refused:roster", "systems must map a name to a participant list")
    declared_provider(polymer)
    return Roster(polymer, systems, hashlib.sha256(raw).hexdigest(), str(path), schema, units,
                  data.get("pairs"))


# the declaration digests an artifact compiled from a deck carries in
# ``provenance.proxy_source`` and a deck or roster must match
BINDING_KEYS = ("chemistry_sha256", "declaration_sha256")
LEGACY_ARTIFACT_FLAG = "--legacy-artifact"


class ArtifactBindingRefusal(ValueError):
    """An artifact that cannot be combined with the given declaration
    (``status`` names why)."""

    def __init__(self, status, reason):
        super().__init__(f"{status}: {reason}")
        self.status = status
        self.reason = reason


def artifact_declaration(artifact: dict) -> dict:
    """The declaration digests the artifact's proxies were built from (empty
    for a compile without a deck: the built-in proxy set or no ``proxy_source``)."""
    provenance = artifact.get("provenance")
    source = provenance.get("proxy_source") if isinstance(provenance, dict) else None
    source = source if isinstance(source, dict) else {}
    return {key: source[key] for key in BINDING_KEYS if key in source}


def bind_artifact(artifact: dict, digests: dict, against: str, legacy_artifact: bool = False) -> tuple:
    """The declaration digests the artifact carries and ``digests`` (those of
    ``against``, a deck or roster) match.  A carried digest that differs is
    refused whatever the flags; an artifact carrying none (a legacy compile)
    is refused unless ``legacy_artifact``, and then binds nothing."""
    carried = artifact_declaration(artifact)
    differ = [f"{key} {str(value)[:16]}... against {str(digests.get(key))[:16]}..."
              for key, value in carried.items() if value != digests.get(key)]
    if differ:
        raise ArtifactBindingRefusal(
            "refused:declaration",
            f"the artifact's declaration digests differ from {against}'s ({'; '.join(differ)}): "
            "the artifact was compiled from another declaration")
    if not carried and not legacy_artifact:
        raise ArtifactBindingRefusal(
            "refused:unbound_artifact",
            f"the artifact carries no declaration digest ({' / '.join(BINDING_KEYS)} in "
            f"provenance.proxy_source), so nothing ties it to {against}'s declaration; pass "
            f"{LEGACY_ARTIFACT_FLAG} to use it unbound")
    return tuple(carried)


def certification_rules(artifact: dict, spec) -> dict:
    """The ``certificates.Certifier`` rules an artifact is certified under:
    for an artifact compiled from a deck (it carries the declaration's
    digests, ``artifact_declaration``), the declaration's own species domain
    (``windows.declared_interim`` over its elements, the domain its windows
    are evaluated under) and its depository growth unit; for a legacy
    artifact none (the v1 rules).  ``spec`` is the declaration bound to the
    artifact (``bind_artifact``)."""
    if not artifact_declaration(artifact):
        return {}
    if spec is None:
        raise ArtifactBindingRefusal("refused:declaration",
                                     "an artifact compiled from a deck is certified under its "
                                     "declaration, and none was given")
    declared = declared_provider(spec)
    policy = win.declared_interim(declared.declared_elements())
    return {"domain": win.certificate_domain(policy), "repeat": declared.unit_formula()}


def make_certifier(fam, rules: dict):
    """The record certifier of emission and of the mapping check."""
    from rmgpy.kmc import certificates as cm
    return cm.Certifier(fam, max_channel_sets=8, max_entries=1024, **rules)


def roster_bytes(roster: Roster) -> bytes:
    """Canonical bytes of a reference roster, used for both writing and its
    provenance hash."""
    document = {"schema": roster.schema, "polymer": roster.polymer.as_mapping(),
                "systems": roster.systems}
    if roster.schema in (ROSTER_SCHEMA_V3, ROSTER_SCHEMA_V4):
        document["proxy_units"] = roster.proxy_units
    if roster.schema == ROSTER_SCHEMA_V4:
        document["pairs"] = roster.pairs
    return (json.dumps(document, sort_keys=True, separators=(",", ":")) + "\n").encode()


def declared_roster_units(artifact: dict) -> Optional[int]:
    """The proxy units of an artifact compiled from a deck (its roster is
    built from those compile structures); None for a legacy artifact."""
    if not artifact_declaration(artifact):
        return None
    units = artifact["provenance"]["proxy_source"].get("proxy_units")
    if isinstance(units, bool) or not isinstance(units, int):
        raise RosterRefusal("refused:roster", f"an artifact compiled from a deck carries no "
                                              f"proxy_units in its proxy_source ({units!r})")
    return units


def declared_roster_pairs(artifact: dict) -> Optional[dict]:
    """The pair proxies, per depth, an artifact compiled from a schema-v2
    deck recorded (None: an artifact of a v1 deck or a legacy one)."""
    if not artifact_declaration(artifact):
        return None
    return artifact["provenance"]["proxy_source"].get("pairs")


def check_roster_for_artifact(roster: Roster, artifact: dict) -> None:
    """A declared artifact is emitted against the roster of its own compile
    structures (same proxy units): v3 for a v1 deck, v4 (the same recorded
    pairs too) for a schema-v2 deck; a legacy artifact against a v2 roster;
    refused by name otherwise."""
    units = declared_roster_units(artifact)
    if units is None and roster.schema != ROSTER_SCHEMA:
        raise RosterRefusal("refused:roster", f"a legacy artifact takes a {ROSTER_SCHEMA} roster, "
                                              f"got {roster.schema}")
    pairs = declared_roster_pairs(artifact)
    schema = ROSTER_SCHEMA_V3 if pairs is None else ROSTER_SCHEMA_V4
    if units is not None and (roster.schema != schema or roster.proxy_units != units):
        raise RosterRefusal("refused:roster", f"an artifact compiled from a deck takes the "
                                              f"{schema} roster of its compile structures "
                                              f"at proxy_units {units}, got {roster.schema} at "
                                              f"{roster.proxy_units}")
    if units is not None and pairs is not None and roster.pairs != pairs:
        raise RosterRefusal("refused:roster", "the roster's pairs are not the pairs the artifact's "
                                              "compile recorded")


def declared_reference_systems(spec, units: int, pairs: Optional[dict] = None) -> Roster:
    """The roster of a declared artifact: one system per compile proxy (its
    participants; both orders for two different structures), every chain
    participant its exact witness spec, a monomer an explicit species.  A v1
    declaration gives the v3 roster of
    ``declared_proxies.declared_structures(spec, units)`` (``declared_chain``
    participants); a schema-v2 one the v4 roster of
    ``declared_proxies.sequence_declared_structures(spec, units, pairs)``
    (``declared_sequence`` participants), ``pairs`` the pair proxies its
    compile recorded (``declared_roster_pairs``)."""
    from dataclasses import asdict
    from rmgpy.kmc.declared_proxies import declared_structures, sequence_declared_structures
    from rmgpy.polymer_spec import (
        is_schema_v2,
        normalized_linear_homopolymer_spec,
        normalized_polymer_spec,
    )

    if is_schema_v2(spec):
        if pairs is None:
            raise RosterRefusal("refused:roster", "the roster of a schema-v2 artifact holds the "
                                                  "pairs its compile recorded; none were given")
        spec = normalized_polymer_spec(spec)
        structures, proxies = sequence_declared_structures(spec, units, pairs)
        kind, schema = "declared_sequence", ROSTER_SCHEMA_V4
    else:
        if pairs is not None:
            raise RosterRefusal("refused:roster", "a v1 declaration's roster takes no recorded "
                                                  "pairs")
        spec = normalized_linear_homopolymer_spec(spec)
        structures, proxies = declared_structures(spec, units)
        kind, schema = "declared_chain", ROSTER_SCHEMA_V3
    declared = declared_provider(spec)
    entries = {}
    for name, (witness, build) in structures.items():
        if witness is None:
            entries[name] = _adjacency_entry(build.molecule)
        elif isinstance(witness, win.DeclaredSpec) and witness.provider == declared.name:
            entries[name] = {kind: asdict(witness.chain)}
        else:
            raise RosterRefusal("refused:roster", f"{name}: a compile structure that is not a "
                                                  f"chain of the declaration ({witness!r})")
    systems = {}
    for site_type, names in proxies.items():
        if len(names) == 1:
            systems[f"single:{site_type}"] = [entries[names[0]]]
        else:
            systems[f"pair:{site_type}:left"] = [entries[n] for n in names]
            if len(set(names)) > 1:
                systems[f"pair:{site_type}:right"] = [entries[n] for n in reversed(names)]
    provisional = Roster(spec, systems, "", "", schema, units,
                         None if pairs is None else {k: [list(p) for p in v]
                                                     for k, v in pairs.items()})
    return replace(provisional, sha256=hashlib.sha256(roster_bytes(provisional)).hexdigest())


def write_roster(roster: Roster, path) -> Roster:
    """Write ``roster`` in the emit CLI's input format and return the loaded
    file-backed roster (whose hash is therefore the exact file hash)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(roster_bytes(roster))
    return load_roster(path)


def _adjacency_entry(molecule) -> dict:
    return {"adjacency_list": molecule.to_adjacency_list()}


def centred_radical_unit(dp: int) -> int:
    """The unit of a ``dp``-unit chain the declared compile proxies put
    their interior radical on (``declared_proxies``); the roster offers it too."""
    return dp // 2


def _chain_entry(owner, dp, head=None, tail=None, radical_unit=None) -> dict:
    return {"chain": {"owner": owner, "dp": dp, "head": head, "tail": tail,
                      "radical_unit": radical_unit}}


def _reference_chain_entries(provider, dps, phases, declared=None) -> list:
    """``(label, roster entry, radical)`` for the provider's capped,
    physical-end and interior-radical reference chains."""
    out = []

    def add(label, dp, head=None, tail=None, radical_unit=None):
        entry = _chain_entry(provider.name, dp, head, tail, radical_unit)
        spec = provider.spec_for_declared_provider(
            declared, dp, head, tail, radical_unit)
        participant = Participant(label, spec, lambda name, s=spec: win.build_spec(s, name))
        built = participant.build("reference_check")
        out.append((label, entry, bool(built.molecule.get_radical_count())))

    def add_interior_radical(dp):
        first = None
        for unit in range(dp):
            entry = _chain_entry(provider.name, dp, radical_unit=unit)
            try:
                spec = provider.spec_for_declared_provider(
                    declared, dp, interior_radical=unit)
                built = win.build_spec(spec, "reference_check")
            except (b.BoundaryError, ValueError):
                continue
            rows = dict(b.rows_of_molecule(built.molecule))
            indices = declared.witness_backbone_indices(rows, built.metadata)
            radicals = [atom for atom, row in rows.items()
                        if row[0] != "H" and row[3] and atom in indices]
            if any({indices[atom] - 1, indices[atom] + 1} <=
                   {indices[n] for n, _order in rows[atom][5] if n in indices}
                   for atom in radicals):
                out.append((f"chain:dp{dp}:interior-radical", entry, True))
                first = unit
                break
        # the compile proxies' own interior radical (centred), when it differs
        centred = centred_radical_unit(dp)
        if centred != first:
            try:
                add(f"chain:dp{dp}:interior-radical:centred", dp, radical_unit=centred)
            except (b.BoundaryError, ValueError):
                pass

    end_table = getattr(provider, "END_CHEMISTRIES", {})
    head_ends = tuple(end_table.get("backbone_head", ()))
    tail_ends = tuple(end_table.get("backbone_tail", ()))
    radical_ends = set(getattr(provider, "RADICAL_END_CHEMISTRIES", ()))
    alkene_ends = set(getattr(provider, "ALKENE_END_CHEMISTRIES", ()))
    closed_head = next((end for end in head_ends
                        if end not in radical_ends | alkene_ends), None)
    closed_tail = next((end for end in tail_ends
                        if end not in radical_ends | alkene_ends), None)
    for dp in dps:
        add(f"chain:dp{dp}:capped", dp)
        if phases and closed_head is not None and closed_tail is not None:
            add(f"chain:dp{dp}:declared", dp, "declared", "declared")
            for end in head_ends:
                if end in radical_ends:
                    add(f"chain:dp{dp}:end-radical:head:{end}", dp, end, closed_tail)
                elif end in alkene_ends:
                    add(f"chain:dp{dp}:alkene:head:{end}", dp, end, closed_tail)
            for end in tail_ends:
                if end in radical_ends:
                    add(f"chain:dp{dp}:end-radical:tail:{end}", dp, closed_head, end)
                elif end in alkene_ends:
                    add(f"chain:dp{dp}:alkene:tail:{end}", dp, closed_head, end)
        add_interior_radical(dp)
    if phases:
        for label, dp, head, tail, radical_unit in provider.reference_special_declarations():
            add(f"chain:special:{label}", dp, head, tail, radical_unit)
    return out


def _reference_small_entries(spec, small) -> list:
    """``(label, roster entry, radical, declared)`` for built-ins and the
    declaration's physical small species."""
    if small is False:
        return []
    selected = None if small is True else set(small)
    declared = {name: graph for name, graph, _disposition in spec.small_species}
    collisions = sorted(set(declared) & {"H", "H2"})
    if collisions:
        raise RosterRefusal("refused:small_name_collision",
                            f"declared small species shadow built-ins: {collisions}")
    out = []
    builtins = {
        "H": "1 H u1 p0 c0\n",
        "H2": "1 H u0 p0 c0 {2,S}\n2 H u0 p0 c0 {1,S}\n",
    }
    for name, adjacency in builtins.items():
        if selected is None or name in selected:
            out.append((f"small:{name}", {"adjacency_list": adjacency}, name == "H", False))
    for name, graph in declared.items():
        if selected is not None and name not in selected:
            continue
        molecule = graph.to_legacy_molecule({})
        molecule.multiplicity = graph.multiplicity
        out.append((f"small:{name}", _adjacency_entry(molecule),
                    bool(molecule.get_radical_count()), True))
    return out


def reference_systems(spec, provider, *, dps=(1, 2, 3, 8), phases=True,
                      small=True, monomer=True, molecularity=frozenset((1, 2))) -> Roster:
    """Generate an emit roster solely from a chain declaration and provider.

    Singles cover every requested chain length and end role.  Bimolecular
    systems contain every ordered radical pair (including two copies for
    identical species), both orders of chain plus H or a declared small
    species, and both orders of chain plus the provider-built monomer.
    """
    from rmgpy.polymer_spec import normalized_linear_homopolymer_spec, refuse_schema_v2
    spec = normalized_linear_homopolymer_spec(refuse_schema_v2(spec, "emit.reference_systems"))
    declared = declared_provider(spec)
    provider = declared if provider is None else provider
    dps = tuple(dps)
    if not dps or any(isinstance(dp, bool) or not isinstance(dp, int) or dp < 1 for dp in dps):
        raise RosterRefusal("refused:roster", f"reference DPs must be positive integers: {dps!r}")
    molecularity = frozenset(molecularity)
    if not molecularity <= {1, 2} or not molecularity:
        raise RosterRefusal("refused:roster", "molecularity must contain 1, 2, or both")

    chains = _reference_chain_entries(provider, dps, bool(phases), declared)
    small_entries = _reference_small_entries(spec, small)
    participants = [(label, entry, radical, "chain") for label, entry, radical in chains]
    participants += [(label, entry, radical, "declared-small" if declared else "small")
                     for label, entry, radical, declared in small_entries]

    monomer_item = None
    if monomer:
        try:
            molecule = provider.reference_monomer(spec)
        except b.BoundaryError as error:
            raise RosterRefusal("refused:monomer", str(error)) from error
        if molecule is not None:
            monomer_item = ("monomer", _adjacency_entry(molecule),
                            bool(molecule.get_radical_count()), "monomer")
            participants.append(monomer_item)

    systems = {}
    if 1 in molecularity:
        for label, entry, _radical, _kind in participants:
            systems[f"single:{label}"] = [entry]
    if 2 in molecularity:
        chains = [item for item in participants if item[3] == "chain"]
        radicals = [item for item in chains if item[2]]
        for radical in radicals:
            for chain in chains:
                prefix = f"pair:radical:{radical[0]}:{chain[0]}"
                systems[f"{prefix}:left"] = [radical[1], chain[1]]
                systems[f"{prefix}:right"] = [chain[1], radical[1]]
        h_item = next((item for item in participants if item[0] == "small:H"), None)
        for chain in [item for item in participants if item[3] == "chain"]:
            if h_item is not None:
                systems[f"pair:H:{chain[0]}:left"] = [chain[1], h_item[1]]
                systems[f"pair:H:{chain[0]}:right"] = [h_item[1], chain[1]]
            for species in [item for item in participants if item[3] == "declared-small"]:
                systems[f"pair:{species[0]}:{chain[0]}:left"] = [chain[1], species[1]]
                systems[f"pair:{species[0]}:{chain[0]}:right"] = [species[1], chain[1]]
            if monomer_item is not None:
                systems[f"pair:monomer:{chain[0]}:left"] = [chain[1], monomer_item[1]]
                systems[f"pair:monomer:{chain[0]}:right"] = [monomer_item[1], chain[1]]

    provisional = Roster(spec, systems, "", "")
    raw = roster_bytes(provisional)
    return Roster(spec, systems, hashlib.sha256(raw).hexdigest(), "")


def _normal_graph(graph):
    """(atom label by key, bond order by key pair, port atom by name) of a graph
    declaration with explicit single-bonded terminal H atoms folded into their
    neighbour's hydrogen count."""
    ports = {name: port.atom for name, port in graph.ports}
    edges = defaultdict(dict)
    for bond in graph.bonds:
        left, right = bond.atoms
        edges[left][right] = bond.order
        edges[right][left] = bond.order
    hydrogens = {atom.key: atom.hydrogens for atom in graph.atoms}
    folded = set()
    for atom in graph.atoms:
        if (atom.element == "H" and atom.isotope == -1 and atom.key not in ports.values()
                and len(edges[atom.key]) == 1 and atom.radical_electrons == 0):
            (other, order), = edges[atom.key].items()
            if order == 1.0:
                hydrogens[other] += 1
                folded.add(atom.key)
    labels = {a.key: (a.element, a.isotope, a.charge, a.radical_electrons, hydrogens[a.key])
              for a in graph.atoms if a.key not in folded}
    bonds = {k: {o: v for o, v in edges[k].items() if o not in folded} for k in labels}
    return labels, bonds, ports


def isomorphic(first, second) -> bool:
    """Whether two graph declarations are the same repeat unit: a bijection of
    atoms preserving element, isotope, charge, radicals, hydrogen count, bond
    orders, and the head and tail ports."""
    (la, ba, pa), (lb, bb, pb) = _normal_graph(first), _normal_graph(second)
    if set(pa) != set(pb) or len(la) != len(lb) or sorted(la.values()) != sorted(lb.values()):
        return False
    if any(first_port.bond_order != second_port.bond_order for (_n, first_port), (_m, second_port)
           in zip(first.ports, second.ports)):
        return False
    order, seen = [], set()
    stack = [pa["head"]]
    while stack:
        key = stack.pop()
        if key not in seen:
            seen.add(key)
            order.append(key)
            stack.extend(sorted(set(ba[key]) - seen, reverse=True))
    if len(order) != len(la):
        return False

    def extend(index, mapping, used):
        if index == len(order):
            return True
        key = order[index]
        for candidate in sorted(lb):
            if candidate in used or lb[candidate] != la[key]:
                continue
            if index == 0 and candidate != pb["head"]:
                continue
            if len(ba[key]) != len(bb[candidate]):
                continue
            if any(mapping[other] not in bb[candidate] or bb[candidate][mapping[other]] != order_
                   for other, order_ in ba[key].items() if other in mapping):
                continue
            if sum(1 for other in ba[key] if other in mapping) != \
                    sum(1 for other in bb[candidate] if other in used):
                continue
            mapping[key] = candidate
            used.add(candidate)
            if extend(index + 1, mapping, used):
                return True
            del mapping[key]
            used.discard(candidate)
        return False
    mapping = {}
    if not extend(0, mapping, set()):
        return False
    return mapping[pa["tail"]] == pb["tail"]


def match_provider(graph, declaration=None):
    """The registered provider whose ``REPEAT_UNIT`` is the declared repeat unit,
    preferring the provider of this exact declaration when supplied, or None."""
    from rmgpy.polymer_spec import GraphDeclaration
    if declaration is not None:
        exact = [provider for provider in b.PROVIDERS.values()
                 if getattr(provider, "spec", None) == declaration]
        if exact:
            return min(exact, key=lambda provider: provider.name)
    for name in sorted(b.PROVIDERS, key=lambda value: (value.startswith("declared:"), value)):
        unit = getattr(b.PROVIDERS[name], "REPEAT_UNIT", None)
        if unit is not None and isomorphic(
                graph, GraphDeclaration.from_mapping(unit, required_ports=("head", "tail"))):
            return b.PROVIDERS[name]
    return None


def declared_provider(spec):
    """The chain provider owned by this exact normalized declaration: the
    residue-sequence provider of a schema-v2 declaration
    (``declared_sequence``)."""
    from rmgpy.polymer_spec import is_schema_v2
    if is_schema_v2(spec):
        from rmgpy.kmc.declared_sequence import sequence_provider
        return sequence_provider(spec)
    exact = match_provider(spec.repeat_unit, spec)
    if exact is not None and getattr(exact, "spec", None) == spec:
        return exact
    from rmgpy.kmc.declared_chain import DeclaredChainProvider
    return DeclaredChainProvider(spec)


def explicit_build(name, smiles=None, adjacency_list=None):
    """A provider-free build of an explicit small molecule: no ports, one
    repeat coordinate per heavy atom."""
    return b.explicit_species_build(name, smiles=smiles, adjacency_list=adjacency_list)


CHAIN_KEYS = {"owner", "dp", "head", "tail", "radical_unit"}


def _owns_declaration(provider, graph):
    exact = getattr(provider, "spec", None)
    if exact is not None:
        return exact == graph
    unit = getattr(provider, "REPEAT_UNIT", None)
    if unit is None:
        return False
    from rmgpy.polymer_spec import GraphDeclaration
    return isomorphic(graph.repeat_unit,
                      GraphDeclaration.from_mapping(unit, required_ports=("head", "tail")))


def _participant(provider, graph, entry, label):
    if not isinstance(entry, dict) or len(entry) != 1:
        raise RosterRefusal("refused:roster", f"{label}: exactly one of chain, declared_chain, "
                                              "smiles, adjacency_list")
    (kind, value), = entry.items()
    if kind in ("declared_chain", "declared_sequence"):
        # a v3 (v4) participant: the exact chain of the roster's declaration
        from rmgpy.kmc.declared_chain import ChainSpec
        from rmgpy.kmc.declared_sequence import SequenceChainProvider, SequenceChainSpec
        chain_type = SequenceChainSpec if isinstance(provider, SequenceChainProvider) else ChainSpec
        if kind != ("declared_sequence" if chain_type is SequenceChainSpec else "declared_chain"):
            raise RosterRefusal("refused:roster", f"{label}: a {kind} participant is not a chain "
                                                  "of this roster's declaration")
        try:
            spec = win.DeclaredSpec(provider.name, chain_type.from_mapping(value))
            participant = Participant(label, spec, lambda name, s=spec: win.build_spec(s, name))
            participant.build("roster_check")
        except Exception as error:   # reported as a refusal, never dropped
            raise RosterRefusal("refused:variant", f"{label}: {type(error).__name__}: {error}")
        return participant
    if kind == "chain":
        if not isinstance(value, dict) or not {"owner", "dp"} <= set(value) or \
                set(value) - CHAIN_KEYS:
            raise RosterRefusal("refused:roster", f"{label}: a chain declares owner, dp and "
                                                  f"optionally head, tail, radical_unit")
        owner_name = value["owner"]
        if not isinstance(owner_name, str) or not owner_name:
            raise RosterRefusal("refused:owner", f"{label}: owner must be a non-empty name")
        owner = b.PROVIDERS.get(owner_name)
        if owner is None:
            raise RosterRefusal("refused:owner", f"{label}: vocabulary owner {owner_name!r} "
                                                 "is not registered")
        if not _owns_declaration(owner, graph):
            raise RosterRefusal("refused:owner", f"{label}: vocabulary owner {owner_name!r} "
                                                 "does not declare this roster's repeat unit")
        try:
            spec = owner.spec_for_declared_provider(
                provider, value["dp"], value.get("head"), value.get("tail"),
                value.get("radical_unit"))
            participant = Participant(label, spec, lambda name, s=spec: win.build_spec(s, name))
            participant.build("roster_check")
        except RosterRefusal:
            raise
        except Exception as error:   # reported as a refusal, never dropped
            raise RosterRefusal("refused:variant", f"{label}: {type(error).__name__}: {error}")
        return participant
    if kind in ("smiles", "adjacency_list"):
        try:
            def make(name, kind=kind, value=value):
                return explicit_build(name, **{kind: value})
            make("roster_check")
        except Exception as error:
            raise RosterRefusal("refused:variant", f"{label}: {type(error).__name__}: {error}")
        return Participant(label, ("small", kind, value), make)
    raise RosterRefusal("refused:roster", f"{label}: unknown participant kind {kind!r}")


def build_reference_systems(roster: Roster):
    """(systems, refusals): ``systems`` maps a name to its participants; a
    system with a refused participant is left out and each refusal is a row
    (``kind: roster_refusal``, system, participant, status, reason)."""
    from rmgpy.polymer_spec import GraphDeclaration  # noqa: F401  (documents the declaration)
    provider = declared_provider(roster.polymer)
    systems, refusals = {}, []
    for name, entries in roster.systems.items():
        built, bad = [], []
        for index, entry in enumerate(entries):
            try:
                built.append(_participant(provider, roster.polymer, entry,
                                          f"{name}[{index}]"))
            except RosterRefusal as error:
                bad.append({"kind": "roster_refusal", "system": name, "participant": index,
                            "status": error.status, "reason": error.reason})
        if bad:
            refusals.extend(bad)
        elif built:
            systems[name] = built
        else:
            refusals.append({"kind": "roster_refusal", "system": name, "participant": None,
                             "status": "refused:roster", "reason": "an empty system"})
    return systems, refusals


# ---------------------------------------------------------------------------
# emission (the acceptance harness' path, on declared systems)
# ---------------------------------------------------------------------------
def load_database(families, path):
    """The compile's own load (``rmgpy.kmc.compile.compile_artifact``)."""
    from rmgpy.data.rmg import RMGDatabase
    from rmgpy.kmc.compiler import prepare_rate_rules
    path = Path(path)
    db = RMGDatabase()
    db.load_kinetics(str(path / "input/kinetics"), reaction_libraries=[], seed_mechanisms=None,
                     kinetics_families=list(families), kinetics_depositories=["training"])
    db.load_thermo(str(path / "input/thermo"), thermo_libraries=["primaryThermoLibrary"],
                   depository=True)
    prepare_rate_rules(db.kinetics, db.thermo, verbose=True)
    return db


def rss_mb():
    try:
        with open("/proc/self/status") as handle:
            for line in handle:
                if line.startswith("VmHWM:"):
                    return int(line.split()[1]) // 1024
    except OSError:
        pass
    return -1


def molecularity(family, direction):
    from rmgpy.kmc import template_emission as te
    return te._molecularity(family, direction)


def reference_cases(certifier, fam, systems, species):
    """Every member case of every family direction on every subset of every
    system's participants of the direction's molecularity, keyed by (family,
    template-side species, template, other-side species)."""
    out = defaultdict(list)
    stats = Counter()
    for name, family in sorted(fam.items()):
        directions = ["forward"] + (["reverse"] if family.reverse_template is not None else [])
        for direction in directions:
            n = molecularity(family, direction)
            for participants in systems.values():
                for subset in sorted(set(combinations(range(len(participants)), n))):
                    chosen = [participants[i] for i in subset]
                    builds = [p.build(f"p{k}") for k, p in enumerate(chosen)]
                    specs = [p.spec for p in chosen]
                    cases = win.member_cases(certifier, family, builds, direction,
                                             label=f"{name}|{direction}|{subset}")
                    stats[f"{name}|{direction}"] += len(cases)
                    for case in cases:
                        mols = win.template_molecules(case.witness, case.direction)
                        other = win.template_molecules(
                            case.witness, "reverse" if case.direction == "forward" else "forward")
                        key = (name, tuple(sorted(species(m) for m in mols)),
                               ";".join(case.pathway), tuple(sorted(species(m) for m in other)))
                        out[key].append((case, specs))
    return out, dict(stats)


def certified_records(records, fam, certifier):
    """The artifact's records as the certified compile publishes them
    (``certify_records`` with ``required``: a required record that does not
    certify is republished refused by status under a new event id).  Returns
    (record dicts in input order, certificate store, replacements old id ->
    new id)."""
    from rmgpy.kmc import certificates as cm
    from rmgpy.kmc.event_record import EventRecord
    records = [EventRecord.from_dict(r) for r in records]
    certified, store, replacements = cm.certify_records(records, fam, certifier, required=True)
    return ([r.to_dict() if hasattr(r, "to_dict") else r for r in certified], store,
            replacements)


def record_certificate(record, store):
    """The stored certificate a certified record references (None: none)."""
    return store.get((record.get("rate_certificate") or {}).get("sha256"))


def record_key(certifier, record, species):
    """(family, template-side species, template, other-side species): the
    record's channel (members of other channels on the same species and
    template are other sites)."""
    _side, mols, others = certifier.template_side(record)
    return (record["family"], tuple(sorted(species(m) for m in mols)), record["template"],
            tuple(sorted(species(m) for m in others)))


def bind(certifier, record, rc, case, owner):
    """``recertify`` the record on the member's template side through an RMG
    isomorphism of the template sides (``windows_acceptance.bind_record``,
    any direction)."""
    from rmgpy.kmc import certificates as cm
    from rmgpy.kmc import windows as win
    _side, mols, _others = certifier.template_side(record)
    proxy = win.template_molecules(case.witness, case.direction)
    roots_r = [cm.aromatic_root(m) for m in mols]
    roots_p = [cm.aromatic_root(m) for m in proxy]
    amap, order = {}, []
    for mr in roots_r:
        for k, mp in enumerate(roots_p):
            if k in order or len(mp.atoms) != len(mr.atoms):
                continue
            maps = mr.find_isomorphism(mp, strict=False)
            if maps:
                order.append(k)
                amap.update({a.id: c.id for a, c in maps[0].items()})
                break
        else:
            return {"same_selection": False, "reason": "no isomorphic member fragment"}
    # the member's molecules in the record's reactant order: a certificate
    # mapping's identity carries the reactant index, so an unequal order
    # (the member's components go by least atom id) swaps the roles
    proxy = [proxy[k] for k in order]
    res = certifier.recertify(record, proxy, amap, certificate=rc, rate_owner=owner)
    return {"same_selection": res["same_selection"], "same_support": res["same_support"],
            "refusals": res["refusals"][:2]}


def select_records(records, chunk=None, limit=None):
    """The event-id-sorted records, then shard ``chunk`` ("i/n": ``i::n``), then
    at most ``limit``."""
    picked = sorted(records, key=lambda r: r["event_id"])
    if chunk:
        k, n = parse_chunk(chunk)
        picked = picked[k::n]
    if limit:
        picked = picked[:limit]
    return picked


def parse_chunk(chunk):
    if chunk is None:
        return 0, 1
    k, n = (int(x) for x in str(chunk).split("/"))
    if not 0 <= k < n:
        raise ValueError(f"bad chunk {chunk!r}: need 0 <= i < n")
    return k, n


def git_head(path) -> Optional[str]:
    try:
        out = subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"], capture_output=True,
                             text=True, timeout=60)
        return out.stdout.strip() or None if out.returncode == 0 else None
    except (OSError, subprocess.SubprocessError):
        return None


def code_sha() -> Optional[str]:
    """HEAD of the code repository, ``-dirty`` when tracked files differ."""
    here = Path(__file__).resolve().parent
    head = git_head(here)
    if head is None:
        return None
    try:
        dirty = subprocess.run(["git", "-C", str(here), "status", "--porcelain", "-uno"],
                               capture_output=True, text=True, timeout=60).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        dirty = ""
    return head + ("-dirty" if dirty else "")


CODE_SOURCE_SUFFIXES = (".py", ".pyx", ".pxd")


def code_sources_sha256(root=None) -> str:
    """sha256 over the sorted relative paths and bytes of every ``.py``,
    ``.pyx`` and ``.pxd`` file under ``root`` (default: the ``rmgpy``
    package, whose code beyond ``rmgpy/kmc`` also decides the chemistry:
    reactions, families, thermo).  The code's identity where git is not
    available (``code_sha`` null).  About 230 files, 8 MB: well under a
    second per shard header or assemble."""
    root = Path(root) if root is not None else Path(__file__).resolve().parent.parent
    digest = hashlib.sha256()
    for path in sorted(p for p in root.rglob("*") if p.suffix in CODE_SOURCE_SUFFIXES
                       and "__pycache__" not in p.parts and p.is_file()):
        digest.update(path.relative_to(root).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def code_extensions_sha256(root=None, modules=None) -> Optional[str]:
    """sha256 over the sorted relative paths and bytes of every compiled
    extension module file under ``root`` (default: the ``rmgpy`` package):
    the code that actually runs where sources were edited next to binaries
    built from another revision.  None (unknown) when a loaded ``rmgpy``
    extension module (``modules``, default ``sys.modules``) comes from a file
    outside that set.  About 50 files, 18 MB: well under a second."""
    import importlib.machinery
    import sys
    root = Path(root) if root is not None else Path(__file__).resolve().parent.parent
    suffixes = tuple(importlib.machinery.EXTENSION_SUFFIXES)
    files = sorted(p for p in root.rglob("*") if p.name.endswith(suffixes)
                   and "__pycache__" not in p.parts and p.is_file())
    known = {p.resolve() for p in files}
    for module in list(modules if modules is not None else sys.modules.values()):
        name, path = getattr(module, "__name__", ""), getattr(module, "__file__", None)
        if (name == "rmgpy" or str(name).startswith("rmgpy.")) and path and \
                str(path).endswith(suffixes) and Path(path).resolve() not in known:
            return None
    digest = hashlib.sha256()
    for path in files:
        digest.update(path.relative_to(root).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def database_content_sha256(path) -> Optional[str]:
    """The database's content digest (``database_provenance``), known where
    git is not (``database_head`` null); None if the tree cannot be read."""
    from rmgpy.kmc.database_provenance import database_content_digest
    try:
        return database_content_digest(path)[0]
    except (OSError, ValueError):
        return None


class RatesByGrid:
    """``ProductionRates`` per temperature grid (cached): the midpoint check
    regenerates a record's own kinetics on the grid midpoints."""

    def __init__(self, database, families, database_path):
        self.database = database
        self.families = list(families)
        self.database_path = database_path
        self._cache = {}

    def __call__(self, grid):
        grid = tuple(grid)
        if grid not in self._cache:
            self._cache[grid] = te.ProductionRates(
                self.database, self.families, database_path=self.database_path,
                database_sha=REFERENCE_RATES_LABEL, temperature_grid=grid)
        return self._cache[grid]


def midpoint_entries(record, payload, rates_for) -> dict:
    """Per emitted primary template with a weight table, the worst error of the
    runtime interpolation at the grid midpoints against the record's OWN
    kinetics there (the harness' ``midpoint_table``, one record): the weight
    table interpolated linearly in 1/T at (T_i + T_i+1) / 2 against ``ssa * k /
    T_ch`` with ``k`` regenerated at the midpoints for the emitted witness's
    channel.  {key: {status: ok, worst, worst_T, n} or {status: unavailable,
    reason}}; a derived template is not entered (its weight follows from its
    primary's)."""
    out = {}
    for t in payload.get("templates", []):
        if t.get("status") != "emitted" or not t.get("weight") or t.get("derived"):
            continue
        wit = ((t.get("maps") or {}).get("witnesses") or {}).get("emitted")
        lrow = next((r for r in (t.get("locality") or {}).get("rows", [])
                     if r.get("label") == "emitted" and r.get("status") == "ok"), None)
        if record is None or wit is None or lrow is None:
            out[t["key"]] = {"status": "unavailable",
                             "reason": "no emitted witness or record to regenerate on"}
            continue
        if not isinstance(t.get("member_reverse"), bool):
            out[t["key"]] = {"status": "unavailable",
                             "reason": "the template names no generation side for its witness"}
            continue
        temps = t["weight"]["T"]
        mids = tuple((a + b) / 2.0 for a, b in zip(temps, temps[1:]))
        r_rows, p_rows = te.rows_from_json(wit["r"]), te.rows_from_json(wit["p"])
        # the witness is stored in template direction: regenerate on the side
        # the template's own member was generated on
        gen, prods = (p_rows, r_rows) if t["member_reverse"] else (r_rows, p_rows)
        # the regenerated k is the EMITTED witness's, so its own count (the one the
        # weight was calibrated with) is the denominator, not the record's
        tch = (lrow.get("counts") or {}).get("matcher")
        try:
            result = rates_for(mids)(te.RateRequest(record["family"], record["template"],
                                                    record["orientation"], gen, prods))
            if result.problems or tuple(result.temperatures) != mids:
                raise te.EmissionError("; ".join(result.problems) or "grid mismatch")
            truth = te.weight_table(mids, result.k, result.ssa_multiplier, tch)
        except Exception as error:   # reported, never hidden
            out[t["key"]] = {"status": "unavailable",
                             "reason": f"{type(error).__name__}: {str(error)[:120]}"}
            continue
        worst, at = 0.0, None
        for T, w in zip(mids, truth):
            err = abs(te._interpolate(temps, t["weight"]["w"], T) / w - 1.0)
            if err > worst:
                worst, at = err, T
        out[t["key"]] = {"status": "ok", "worst": worst, "worst_T": at, "n": len(mids)}
    return out


class Emitter:
    """The emission session of one artifact against a roster's systems."""

    def __init__(self, artifact, database, systems, families=None, contexts=True, polymer=None):
        self.artifact = Path(artifact)
        self.database = Path(database)
        self.systems = systems
        # the roster's declaration: a declared artifact is certified under it
        self.polymer = polymer
        self.only = list(families) if families else None
        self.contexts = contexts
        self.state = {}

    def init(self):
        from rmgpy.kmc import certificates as cm
        from rmgpy.kmc import channels as chm
        from rmgpy.kmc.reference_thermo import GasPhaseRMGReferenceThermo
        S = self.state
        art = json.loads(self.artifact.read_text())
        S["artifact_records"] = len(art["records"])
        if self.only:
            # one family's records (certified among themselves: a reversed
            # record's rate owner is in its own family)
            art["records"] = [r for r in art["records"] if r["family"] in self.only]
        families = sorted({r["family"] for r in art["records"]})
        t = time.time()
        db = load_database(families, self.database)
        fam = db.kinetics.families
        S.update(db=db, fam=fam, families=families, load_s=round(time.time() - t, 1))
        certifier = make_certifier(fam, certification_rules(art, self.polymer))
        t = time.time()
        certified, store, _replacements = certified_records(art["records"], fam, certifier)
        S.update(certify_s=round(time.time() - t, 1), records=certified, store=store,
                 by_id={r["event_id"]: r for r in certified}, certifier=certifier)
        species = chm.SpeciesKeys()
        t = time.time()
        cases, stats = reference_cases(certifier, fam, self.systems, species)
        S.update(cases=cases, case_stats=stats, cases_s=round(time.time() - t, 1), species=species,
                 bounds={n: win.span_bound(win.family_constants(f)) for n, f in fam.items()},
                 counter=te.RmgCounter(fam),
                 kc=te.ThermoKc(GasPhaseRMGReferenceThermo(db.thermo, REFERENCE_THERMO_LABEL)),
                 rates=te.ProductionRates(db, families, database_path=self.database,
                                          database_sha=REFERENCE_RATES_LABEL),
                 rates_by_grid=RatesByGrid(db, families, self.database))

    @property
    def records(self):
        return self.state["records"]

    def run_record(self, record):
        from rmgpy.kmc import certificates as cm
        S = self.state
        contexts = self.contexts
        certifier, fam = S["certifier"], S["fam"]
        start, cpu = time.time(), time.process_time()
        rc = record_certificate(record, S["store"])
        row = {"kind": "row", "event_id": record["event_id"][:24], "family": record["family"],
               "orientation": record["orientation"], "template": record["template"][:100],
               "record_status": record["status"], "rate_kind": (record.get("rate_source") or {})
               .get("kind"), "certified": (rc or {}).get("certified")}
        family = fam[record["family"]]
        why = te.upstream_refusal(record, rc, S["by_id"])
        members, member_rows = [], []
        try:
            if why is None:
                key = record_key(certifier, record, S["species"])
                found = S["cases"].get(key, [])
                row["members"] = len(found)
                if found:
                    owner = cm._rate_owner(record, S["by_id"])
                    row["bind"] = bind(certifier, record, rc, found[0][0], owner)
                    if not row["bind"].get("same_selection"):
                        why = f"refused:recertify: {row['bind']}"[:300]
                if why is None:
                    for case, specs in found:
                        out = []
                        keep = []
                        res = win.run_member(certifier, family, case,
                                             bound=S["bounds"][record["family"]],
                                             contexts=contexts, specs=specs, keep=keep,
                                             contexts_out=out)
                        if out:
                            ctx = out[0]
                        elif keep:
                            fp = keep[0]
                            ctx = ((win.ContextResult(tuple((p, None) for p in
                                                            win.b.vector_ports(case.witness)),
                                                      fp.status, fp.reason, fp.witness, fp.window,
                                                      fp.observation, fp),), ())
                        else:
                            ctx = ((), ())
                        members.append((case, res, ctx))
                        member_rows.append({k: res.get(k) for k in (
                            "label", "status", "reason", "window_atoms", "context_classes",
                            "covered", "seconds")})
            if why is not None:
                payload = te.finish({"schema": te.SCHEMA, "status": "refused",
                                     "record_event_id": record["event_id"],
                                     "refusal": {"by": "upstream", "status": te.status_of(why),
                                                 "reason": why}, "templates": [],
                                     "contexts": []})
            else:
                payload = te.emit_record(record, family=family, members=members,
                                         counter=S["counter"], rates=S["rates"],
                                         record_certificate=rc, by_id=S["by_id"],
                                         check_upstream=False, kc=S["kc"])
            row["payload"] = payload
            row["status"] = payload["status"]
            row["refusal"] = (payload.get("refusal") or {}).get("status")
            row["midpoint"] = midpoint_entries(record, payload, S["rates_by_grid"])
        except te.EmissionError as error:
            row.update(status="error", refusal="error:emission", reason=str(error)[:400])
        except Exception:   # reported, never hidden
            row.update(status="error", refusal="error:exception",
                       reason=traceback.format_exc()[-1500:])
        row["member_rows"] = member_rows
        row["seconds"] = round(time.time() - start, 1)
        row["cpu"] = round(time.process_time() - cpu, 1)
        row["rss_mb"] = rss_mb()
        return row


def shard_header(emitter, roster, refusals, picked, chunk, partial):
    from rmgpy.kmc.compiler import DEFAULT_T_GRID
    S = emitter.state
    artifact = emitter.artifact
    return {"kind": "header", "schema": ROWS_SCHEMA, "artifact_id": artifact.stem,
            "artifact_sha256": hashlib.sha256(artifact.read_bytes()).hexdigest(),
            "artifact_records": S["artifact_records"], "families": S["families"],
            "database_path": str(emitter.database), "database_head": git_head(emitter.database),
            "database_content_sha256": database_content_sha256(emitter.database),
            "code_sha": code_sha(), "code_sources_sha256": code_sources_sha256(),
            "code_extensions_sha256": code_extensions_sha256(),
            "roster_sha256": roster.sha256,
            "temperatures": [float(t) for t in DEFAULT_T_GRID], "chunk": chunk,
            "selection": "partial" if partial else "all",
            "assigned": [r["event_id"] for r in picked], "roster_refusals": refusals,
            "systems": {n: [p.label for p in ps] for n, ps in emitter.systems.items()}}


def is_partial_selection(only, limit, artifact_records, selected_records) -> bool:
    """Whether the emitted records are less than the artifact: ``--limit``, or
    ``--family`` selecting fewer records than the artifact holds (a one-family
    artifact, or ``--family`` naming all its families, is the whole artifact)."""
    if limit:
        return True
    return bool(only) and selected_records != artifact_records


def run_emit(args):
    roster = load_roster(args.roster)
    # the roster's declaration must be the artifact's, before anything is built
    try:
        artifact = json.loads(Path(args.artifact).read_text())
        bound_by = bind_artifact(artifact, roster.digests(), f"the roster {args.roster}",
                                 args.legacy_artifact)
        check_roster_for_artifact(roster, artifact)
    except (ArtifactBindingRefusal, RosterRefusal) as error:
        print("REFUSED", error, flush=True)
        return 1
    print("BINDING", {"bound_by": list(bound_by), "legacy_artifact": not bound_by}, flush=True)
    systems, refusals = build_reference_systems(roster)
    for refusal in refusals:
        print("ROSTER", json.dumps(refusal), flush=True)
    emitter = Emitter(args.artifact, args.database, systems, families=args.family,
                      contexts=not args.no_contexts, polymer=roster.polymer)
    t0 = time.time()
    emitter.init()
    picked = select_records(emitter.records, args.chunk, args.limit)
    header = shard_header(emitter, roster, refusals, picked, args.chunk,
                          partial=is_partial_selection(args.family, args.limit,
                                                       emitter.state["artifact_records"],
                                                       len(emitter.records)))
    print("INIT", {k: emitter.state.get(k) for k in ("load_s", "certify_s", "cases_s", "case_stats")},
          len(emitter.records), "records, selected", len(picked), flush=True)
    errors = 0
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as log:
        log.write(json.dumps(header, default=str) + "\n")
        for refusal in refusals:
            log.write(json.dumps(refusal) + "\n")
        for k, record in enumerate(picked):
            row = emitter.run_record(record)
            errors += row["status"] == "error"
            log.write(json.dumps(row, default=str) + "\n")
            log.flush()
            print(f"[{time.time() - t0:6.0f}s] {k + 1}/{len(picked)} {row['event_id'][:16]} "
                  f"{row['family']} {row['orientation']} -> {row['status']} "
                  f"{row.get('refusal') or ''} t={row['seconds']}s rss={row['rss_mb']}MB",
                  flush=True)
    return 1 if errors else 0


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------
@dataclass
class Services:
    """What the balance functions need from the database."""

    counter: Any
    kc: Any
    # the content digest of the database the counter and thermo come from
    database_content_sha256: Optional[str] = None


def load_services(database, families):
    from rmgpy.kmc.reference_thermo import GasPhaseRMGReferenceThermo
    db = load_database(families, database)
    return Services(te.RmgCounter(db.kinetics.families),
                    te.ThermoKc(GasPhaseRMGReferenceThermo(db.thermo, REFERENCE_THERMO_LABEL)),
                    database_content_sha256(database))


def dumps(catalog) -> str:
    return json.dumps(catalog, sort_keys=True, indent=1, allow_nan=False) + "\n"


def _canon(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"))


def read_rows(paths):
    """([(path, header or None, rows)], failures): rows files in a stable order."""
    files, failures = [], []
    for path in paths:
        header, rows, refusals = None, [], []
        for number, line in enumerate(Path(path).read_text().splitlines(), 1):
            if not line.strip():
                continue
            try:
                item = json.loads(line)
            except ValueError as error:
                failures.append(f"{path}:{number}: unreadable row ({error})")
                continue
            kind = item.get("kind")
            if kind == "header":
                if header is not None:
                    failures.append(f"{path}: more than one header")
                header = item
            elif kind == "roster_refusal":
                refusals.append(item)
            else:
                rows.append(item)
        files.append({"path": str(path), "header": header, "rows": rows})
        if header is None:
            failures.append(f"{path}: no shard header (artifact, database, roster and shard "
                            f"identity unknown)")

    def order(f):
        h = f["header"] or {}
        return (parse_chunk(h.get("chunk")) if h else (10 ** 9, 1), _canon(h), f["path"])
    files.sort(key=order)
    return files, failures


def stripped_payloads(payloads):
    """The fields ``loaded_templates`` reads, emitted templates only."""
    out = []
    for p in payloads:
        templates = [{k: t[k] for k in PAYLOAD_FIELDS if k in t}
                     for t in p.get("templates", []) if t.get("status") == "emitted"]
        if templates:
            out.append({"record_event_id": p["record_event_id"], "family": p.get("family"),
                        "reverse_of": p.get("reverse_of"), "templates": templates})
    return sorted(out, key=lambda p: p["record_event_id"])


def _refusal_reason(row):
    payload = row.get("payload") or {}
    if row.get("reason"):
        return row["reason"]
    if (payload.get("refusal") or {}).get("reason"):
        return payload["refusal"]["reason"]
    refused = Counter((t.get("refusal") or {}).get("status") for t in payload.get("templates", [])
                      if t.get("status") == "refused")
    return ", ".join(f"{k} x{v}" for k, v in sorted(refused.items(), key=str)) or None


def midpoint_gate(rows, payloads):
    """Aggregate the rows' midpoint entries over the final emitted templates:
    the worst error, the templates above the threshold, and every template
    whose check is not evaluable (a gate that evaluates nothing must not pass).
    Derived templates are counted ``unavailable`` and not gated, as in the
    harness.  Returns the gate table plus ``failures``."""
    by_record = {(r.get("payload") or {}).get("record_event_id"): r for r in rows if r.get("payload")}
    worst, at, who, n = 0.0, None, None, 0
    unavailable, unverified, failing, failures = Counter(), [], [], []
    for p in payloads:
        rid = p.get("record_event_id") or ""
        entries = (by_record.get(rid) or {}).get("midpoint") or {}
        for t in p.get("templates", []):
            if t.get("status") != "emitted" or not t.get("weight"):
                continue
            name = f"{rid[:16]}:{t['key'][:12]}"
            if t.get("derived"):
                unavailable["derived"] += 1
                continue
            entry = entries.get(t["key"])
            if entry is None:
                unverified.append([name, "no midpoint entry in the row (not computed)"])
            elif entry.get("status") != "ok":
                unverified.append([name, entry.get("reason", "unavailable")])
            else:
                n += 1
                if entry["worst"] > worst:
                    worst, at, who = entry["worst"], entry["worst_T"], name
                if entry["worst"] > MIDPOINT_THRESHOLD:
                    failing.append([name, entry["worst"], entry["worst_T"]])
    for name, error, T in failing:
        failures.append(f"primary-weight midpoint error {error:.3g} > "
                        f"{MIDPOINT_THRESHOLD:.0%} at {T} K ({name})")
    for name, reason in unverified:
        failures.append(f"primary-weight midpoint check not evaluable for {name}: {reason}")
    return {"evaluated": not unverified, "threshold": MIDPOINT_THRESHOLD, "templates": n,
            "worst": worst, "worst_T": at, "worst_template": who, "failing": failing,
            "unverified": unverified, "unavailable": dict(unavailable), "failures": failures}


class AssembleRefused(RuntimeError):
    """Rows that this assemble cannot vouch for: emitted with another code or
    database than the one deriving and checking them here."""


IDENTITY_FIELDS = ("database_content_sha256", "code_sources_sha256", "code_extensions_sha256")


def assembler_identity(database) -> dict:
    """This process's own database content digest and code digest."""
    return {"database_content_sha256": database_content_sha256(database),
            "code_sources_sha256": code_sources_sha256(),
            "code_extensions_sha256": code_extensions_sha256()}


def assemble(paths, *, database=None, services=None) -> dict:
    """The catalog of the rows files ``paths`` (see the module docstring).
    Refused (``AssembleRefused``) unless every shard header carries the
    database content digest and the code digest of this process: the
    derivation and the cycle check run here, with this code and database."""
    failures = []
    files, read_failures = read_rows(paths)
    failures += read_failures
    headers = [f["header"] for f in files if f["header"]]
    first = headers[0] if headers else {}
    own = assembler_identity(database or first.get("database_path"))
    for field in IDENTITY_FIELDS:
        theirs = sorted({str(h.get(field)) for h in headers})
        if own[field] is None or any(h.get(field) != own[field] for h in headers):
            raise AssembleRefused(
                f"{field}: the shards were emitted with {theirs}, this assemble runs with "
                f"{own[field]!r}; derivation and the balance check would use other code or data "
                "than the rows record (emit and assemble with the same code and database)")
    if services is not None and \
            getattr(services, "database_content_sha256", None) != own["database_content_sha256"]:
        raise AssembleRefused(
            f"database_content_sha256: the supplied services were built from "
            f"{getattr(services, 'database_content_sha256', None)!r}, this assemble checks "
            f"{own['database_content_sha256']!r}; the balance check would run on another database "
            "than the provenance names")
    for field in HEADER_MATCH:
        values = {_canon(h.get(field)) for h in headers}
        if len(values) > 1:
            failures.append(f"shard headers disagree on {field}: {sorted(values)}")
    # coverage: every shard once, every assigned record a row, the whole artifact assigned
    ns = {parse_chunk(h.get("chunk"))[1] for h in headers}
    complete_coverage = bool(headers)
    if len(ns) > 1:
        failures.append(f"shard counts disagree: {sorted(ns)}")
        complete_coverage = False
    elif ns:
        (n,) = ns
        seen = Counter(parse_chunk(h.get("chunk"))[0] for h in headers)
        missing = [k for k in range(n) if not seen[k]]
        if len(missing) <= 8:
            failures += [f"missing shard {k}/{n}" for k in missing]
        elif missing:
            failures.append(f"missing shards {missing[0]}/{n} .. {missing[-1]}/{n} "
                            f"({len(missing)} of {n})")
        complete_coverage = complete_coverage and not missing
        for k in range(n):
            if seen[k] > 1:
                failures.append(f"duplicate shard {k}/{n} ({seen[k]} files)")
                complete_coverage = False
    if first.get("selection", "all") != "all":
        failures.append("the rows are a partial selection (--family or --limit)")
        complete_coverage = False
    by_id = {}
    for f in files:
        for row in f["rows"]:
            by_id.setdefault(row["event_id"], row)
    assigned = {i[:24] for h in headers for i in h.get("assigned", [])}
    if headers and first.get("selection", "all") == "all" and \
            len(assigned) != first.get("artifact_records"):
        failures.append(f"the shards assign {len(assigned)} of the artifact's "
                        f"{first.get('artifact_records')} records")
        complete_coverage = False
    if headers:
        absent = sorted(assigned - set(by_id))
        if absent:
            failures.append(f"{len(absent)} assigned record(s) have no row, e.g. {absent[0]}")
            complete_coverage = False
        extra = sorted(set(by_id) - assigned)
        if extra:
            failures.append(f"{len(extra)} row(s) were not assigned to any shard, e.g. {extra[0]}")
    roster_refusals = {}
    for f in files:
        for h in ([f["header"]] if f["header"] else []):
            for item in h.get("roster_refusals", []):
                roster_refusals[_canon(item)] = item
    for item in sorted(roster_refusals.values(), key=_canon):
        failures.append(f"roster {item.get('status')}: system {item.get('system')} "
                        f"participant {item.get('participant')}: {item.get('reason')}")
    rows = [by_id[k] for k in sorted(by_id)]
    for row in rows:
        if row.get("status") == "error":
            failures.append(f"error {row['event_id']}: {str(row.get('reason', ''))[:300]}")
    if rows and not any(r.get("status") in ("emitted", "partial") for r in rows):
        failures.append("every record was refused (a vacuous catalog)")

    have = [k for k, r in enumerate(rows) if r.get("payload")]
    payloads = [rows[k]["payload"] for k in have]
    # every key that enters collision resolution (refused copies included): a
    # catalog set compares them across its members (catalog_set)
    entry_keys = sorted({t["key"] for p in payloads for t in p.get("templates", [])
                         if t.get("status") == "emitted"})
    collisions = []
    try:
        merged, collisions = te.resolve_collisions(payloads)
    except Exception as error:    # a split reciprocal pair, a bad weight
        # the catalog is refused: unmerged copies would load as independent
        # templates and fail (or fake) every gate after this one
        failures.append(f"catalog consistency: {type(error).__name__}: {error}"[:600])
        have, merged = [], []
    refused_collisions = [c for c in collisions if c["disposition"] != "merged"]
    for c in refused_collisions:
        failures.append(f"{c['disposition']} {c['key'][:16]}: {c['copies']} copies, "
                        f"log10 diameter {c['log10_diameter']} (records {c['records']})")
    payloads = merged
    cycle_table = []
    if payloads:
        try:
            services = services or load_services(database or first.get("database_path"),
                                                 first.get("families", []))
            payloads = te.derive_weights(payloads, counter=services.counter)
            payloads, cycle_table = te.check_cycles(payloads, counter=services.counter,
                                                    kc=services.kc, complete=complete_coverage)
        except Exception as error:   # reported, never hidden
            failures.append(f"balance check could not run: {type(error).__name__}: "
                            f"{str(error)[:500]}")
    for k, payload in zip(have, payloads):
        rows[k]["payload"] = payload
        rows[k]["status"] = payload["status"]
        rows[k]["refusal"] = (payload.get("refusal") or {}).get("status")
    refused_pairs = [r for r in cycle_table if r.get("status") == "refused:cycle"]
    failures += te.cycle_gates(cycle_table, payloads, complete_coverage)
    if cycle_table and len(refused_pairs) > CYCLE_REFUSED_FRACTION * len(cycle_table):
        failures.append(f"more than {CYCLE_REFUSED_FRACTION:.0%} of the loaded reciprocal pairs are "
                        f"refused by the cycle check: {len(refused_pairs)} of {len(cycle_table)}")
    calibration = []
    for row in rows:
        p = row.get("payload") or {}
        for where, cal in [("record", p.get("calibration"))] + \
                [(f"template {t.get('key', '')[:12]}", t.get("calibration"))
                 for t in p.get("templates", []) if t.get("status") == "emitted"]:
            for c in cal or []:
                calibration.append(c)
                if not c.get("ok"):
                    failures.append(f"calibration {row['event_id']} {where} at {c.get('T')} K: {c}")
    checked = [r for r in cycle_table if r.get("status") != "not_evaluable" and r.get("rows")]
    balance = {
        "calibration_rows": len(calibration),
        "calibration_max_relative_error": max([c["relative_error"] for c in calibration
                                               if c.get("relative_error") is not None], default=None),
        "calibration_ok": all(c.get("ok") for c in calibration),
        "worst_ratio_error": max([r["max_ratio_error"] for r in checked], default=None),
        "worst_grid_error": max([r.get("max_grid_error", 0.0) for r in cycle_table], default=0.0),
        "worst_offgrid_error": max([r.get("max_offgrid_error", 0.0) for r in cycle_table],
                                   default=0.0),
        "self_inverse_templates": sum(1 for p in payloads for t in p.get("templates", [])
                                      if t.get("status") == "emitted" and t.get("self_inverse"))}
    midpoint = midpoint_gate(rows, payloads)
    failures += midpoint.pop("failures")
    out_payloads = stripped_payloads(payloads)
    expected = sum(1 for p in payloads for t in p.get("templates", [])
                   if t.get("status") == "emitted")
    temperatures = [float(t) for t in first.get("temperatures", [])]
    round_trip = {"temperatures": len(temperatures), "templates": expected, "ok": True}
    for T in temperatures:
        try:
            loaded = len(te.loaded_templates(out_payloads, T))
        except Exception as error:
            round_trip["ok"] = False
            failures.append(f"round trip at {T} K: {type(error).__name__}: {str(error)[:300]}")
            break
        if loaded != expected:
            round_trip["ok"] = False
            failures.append(f"round trip at {T} K loads {loaded} templates, rows emitted {expected}")
            break
    if not temperatures:
        round_trip["ok"] = False
        failures.append("round trip: no provenance temperatures")
    completeness = defaultdict(Counter)
    refusals = []
    for row in rows:
        completeness[f"{row.get('family')}|{row.get('orientation')}"][row.get("status")] += 1
        if row.get("status") != "emitted":
            full = (row.get("payload") or {}).get("record_event_id") or row["event_id"]
            refusals.append({"record_event_id": full, "family": row.get("family"),
                             "orientation": row.get("orientation"), "status": row.get("status"),
                             "refusal": row.get("refusal"), "reason": _refusal_reason(row)})
    completeness = {k: {s: v.get(s, 0) for s in ("emitted", "partial", "refused", "error")}
                    for k, v in sorted(completeness.items())}
    failures = list(dict.fromkeys(failures))
    return {
        "schema": SCHEMA,
        "provenance": {k: first.get(k) for k in (
            "artifact_id", "artifact_sha256", "artifact_records", "database_path", "database_head",
            "database_content_sha256", "code_sha", "code_sources_sha256",
            "code_extensions_sha256", "roster_sha256",
            "families", "temperatures", "systems")},
        "payloads": out_payloads,
        "refusals": refusals,
        "gates": {
            "collisions": {"keys": len(collisions),
                           "merged": len(collisions) - len(refused_collisions),
                           "refused": len(refused_collisions),
                           "copies": sum(c["copies"] for c in collisions),
                           "tolerance_log10": te.LOCALITY_TOLERANCE_LOG10,
                           "entry_keys": entry_keys},
            "cycles": {"pairs": len(cycle_table), "refused": len(refused_pairs),
                       "not_evaluable": len(te.cycle_not_evaluable(cycle_table)),
                       "by_status": dict(Counter(r.get("status") for r in cycle_table))},
            "balance": balance,
            "round_trip": round_trip,
            "coverage": {"complete": complete_coverage, "rows": len(rows),
                         "assigned": len(assigned)},
            "midpoint": midpoint,
            "failures": failures},
        "completeness": completeness,
        "complete": not failures,
    }


def run_assemble(args):
    services = None
    try:
        catalog = assemble(args.rows, database=args.database, services=services)
    except AssembleRefused as error:
        print("REFUSED", error, flush=True)
        return 1
    Path(args.out).write_text(dumps(catalog))
    print("CATALOG", args.out, "complete" if catalog["complete"] else "INCOMPLETE",
          len(catalog["payloads"]), "payloads", flush=True)
    for failure in catalog["gates"]["failures"]:
        print("GATE", failure, flush=True)
    return 0 if catalog["complete"] else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "assemble":
        parser = argparse.ArgumentParser(prog="python -m rmgpy.kmc.emit assemble")
        parser.add_argument("rows", nargs="+", help="rows files (shards) of one artifact")
        parser.add_argument("--out", required=True)
        parser.add_argument("--database", default=None,
                            help="RMG-database path for the balance check (default: the "
                                 "shard headers')")
        return run_assemble(parser.parse_args(argv[1:]))
    if argv and argv[0] == "catalog-set":
        parser = argparse.ArgumentParser(prog="python -m rmgpy.kmc.emit catalog-set")
        parser.add_argument("catalogs", nargs="+", help="one kmc-template-catalog/1 file per family")
        parser.add_argument("--out", required=True)
        args = parser.parse_args(argv[1:])
        try:
            doc = build_catalog_set(args.catalogs, args.out)
        except CatalogSetError as error:
            print("REFUSED", error, flush=True)
            return 1
        print("CATALOG-SET", args.out, len(doc["catalogs"]), "catalogs,",
              len(doc["dependencies"]), "dependencies", flush=True)
        return 0
    parser = argparse.ArgumentParser(prog="python -m rmgpy.kmc.emit")
    parser.add_argument("artifact")
    parser.add_argument("roster")
    parser.add_argument("--database", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--chunk", default=None, help="i/n: records i::n of the sorted list")
    parser.add_argument("--family", action="append", default=None,
                        help="only these families' records (repeatable; a partial catalog unless "
                             "they are all the artifact's records)")
    parser.add_argument("--limit", type=int, default=None, help="at most N records (a partial catalog)")
    parser.add_argument("--no-contexts", action="store_true")
    parser.add_argument(LEGACY_ARTIFACT_FLAG, action="store_true",
                        help="emit an artifact that carries no declaration digest (a compile "
                             "without --deck) against the roster's declaration, unbound")
    return run_emit(parser.parse_args(argv[1:] if argv and argv[0] == "run" else argv))


if __name__ == "__main__":
    sys.exit(main())
