"""I-125 probe: can the kMC runtime run a declared melt of long PS chains?

Not collected by pytest (no ``Test.py`` suffix); ``longChainRuntimeTest.py``
imports the helpers below.  Run as a script to print the evidence for the
claims C1..C7:

    python long_chain_runtime_probe.py report [C1 C2 ...]

and to cold-compile the reduced artifact (five families, PS end-radical proxies):

    python long_chain_runtime_probe.py build DATABASE OUTPUT_DIR

Environment: I125_DATABASE (RMG-database snapshot), I125_ARTIFACT_DIR
(where the reduced artifact lives; built there when absent).

Nothing here changes product code: the declared melt is materialized into
``KMCState`` strands by ``materialize`` below because the repository has no
such step (itself the finding behind C1/C3).
"""

import copy
import hashlib
import itertools
import json
import os
import subprocess
import sys
import time
import uuid
from collections import Counter
from pathlib import Path

import networkx as nx
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
DATABASE = os.environ.get("I125_DATABASE", "/home/alon/runs/phase2b/database")
ARTIFACT_DIR = Path(os.environ.get("I125_ARTIFACT_DIR", "/home/alon/runs/i125/art"))
FAMILIES = (
    "Disproportionation",
    "H_Abstraction",
    "R_Addition_MultipleBond",
    "R_Recombination",
    "intra_H_migration",
)
PROXIES = (
    "pristine",
    "end_radical",
    "benzylic_end_radical",
    "end_radical+styrene",
    "benzylic_end_radical+styrene",
    "interior_radical",
)
DP = 20


# --------------------------------------------------------------------------
# reduced artifact
# --------------------------------------------------------------------------
def _build_main(database, output):
    """Cold-compile the reduced artifact (recipe of the I-122 probe driver)."""
    from rmgpy.kmc import compile as compile_module
    from rmgpy.kmc.compile import CompileRequest, compile_artifact

    original = compile_module.ps_proxy_set
    wanted = set(PROXIES)

    def proxy_set(*options, **keywords):
        return tuple(
            proxy
            for proxy in original(*options, **keywords)
            if proxy.site_type in wanted
        )

    compile_module.ps_proxy_set = proxy_set
    start = time.time()
    path, _ = compile_artifact(
        CompileRequest(
            database_path=database,
            output_dir=output,
            families=FAMILIES,
            family_universe=FAMILIES,
            database_sha="determinism-probe",
        )
    )
    digest = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    print(f"RESULT {digest} {time.time() - start:.1f}s", flush=True)


def build_artifact(directory=None, database=None):
    """Return the reduced artifact path, building it in a subprocess if absent."""
    directory = Path(directory or ARTIFACT_DIR)
    existing = sorted(directory.glob("*.json"))
    if existing:
        return existing[0]
    directory.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ)
    environment.update(
        PYTHONHASHSEED="0",
        PYTHONPATH=str(REPO_ROOT),
        MPLCONFIGDIR=str(directory / "matplotlib"),
    )
    with (directory / "build.stdout.log").open("w") as out, (
        directory / "build.stderr.log"
    ).open("w") as err:
        subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "build",
             str(database or DATABASE), str(directory)],
            cwd=REPO_ROOT, env=environment, stdout=out, stderr=err,
            check=True, timeout=900,
        )
    return sorted(directory.glob("*.json"))[0]


def load_records(path=None):
    from rmgpy.kmc.event_record import EventRecord

    path = Path(path) if path else build_artifact()
    artifact = json.loads(path.read_bytes())
    return tuple(EventRecord.from_dict(data) for data in artifact["records"])


# --------------------------------------------------------------------------
# declared PS chain -> KMCState strands
# --------------------------------------------------------------------------
def _atom(key, element="C", hydrogens=0):
    return {"key": key, "element": element, "hydrogens": hydrogens,
            "electronic_state": {"radical_electrons": 0}}


def ps_spec(dp=DP, counts=None):
    """Declared PS: C8H8 repeat unit (head c1 / tail c2), H end groups."""
    from rmgpy.polymer_spec import LinearHomopolymerSpec

    atoms = [_atom("c1", hydrogens=2), _atom("c2", hydrogens=1)]
    atoms += [_atom(f"p{i}", hydrogens=0 if i == 1 else 1) for i in range(1, 7)]
    bonds = [(("c1", "c2"), 1), (("c2", "p1"), 1)]
    bonds += [((f"p{i}", f"p{i + 1}"), 1.5) for i in range(1, 6)]
    bonds.append((("p1", "p6"), 1.5))
    repeat = {
        "atoms": atoms,
        "bonds": [{"atoms": list(pair), "order": order} for pair, order in bonds],
        "ports": {"head": {"atom": "c1", "bond_order": 1},
                  "tail": {"atom": "c2", "bond_order": 1}},
        "electronic_state": {"multiplicity": 1},
    }

    def cap(port):
        return {"atoms": [_atom("h", "H")],
                "bonds": [], "ports": {port: {"atom": "h", "bond_order": 1}},
                "electronic_state": {"multiplicity": 1}}

    return LinearHomopolymerSpec.from_mapping({
        "schema_version": 1, "label": "PS", "repeat_unit": repeat,
        "end_groups": {"head": cap("cap"), "tail": cap("cap")},
        "small_species": {},
        "initial_ensemble": {"dp_counts": dict(counts or {dp: 1}),
                             "rounding": "largest_remainder_dp_ascending_v1"},
        "initial_mass_kg": 1.25, "density_kg_m3": 1050.0,
        "locality_policy": {"radius": 2, "bounds": {
            "elements": ["C", "H"], "isotopes": [-1], "formal_charges": [0],
            "radical_electrons": [0, 1], "lone_pairs": [0],
            "bond_orders": [1, 1.5, 2, 3], "spin_multiplicities": [1, 2, 3],
            "maximum_heavy_atom_degree": 4,
            "maximum_explicit_hydrogens_per_heavy_atom": 4,
            "maximum_total_hydrogens_on_free_small_species": 12,
            "maximum_family_arity": 2, "maximum_root_tuple_size": 4,
            "maximum_free_small_species_heavy_atoms": 8,
            "attachment_classes": [], "topology_metadata_schema": {}}},
        "family_manifest": {}, "library_manifest": {},
        "transport_model": {}, "diagnostics": {},
    })


class ChainGraph:
    """Explicit-H graph of one declared chain with its backbone order."""

    def __init__(self, dp, name, radicals=()):
        from rmgpy.kmc.state import _parse_adjacency
        from rmgpy.polymer_shared_support import _authoritative_capped_graph

        self.dp, self.name = dp, name
        molecule = _authoritative_capped_graph(ps_spec(dp), dp)
        self.nodes = _parse_adjacency(molecule.to_adjacency_list())
        self.backbone = self._ordered_backbone()
        self.removed = []
        for kind in radicals:
            self.abstract(self.site(kind))

    def _ordered_backbone(self):
        nodes = self.nodes
        backbone = {i for i, n in nodes.items() if n["element"] == "C"
                    and all(float(o) == 1.0 for o in n["edges"].values())}
        ends = sorted(i for i in backbone if len(set(nodes[i]["edges"]) & backbone) == 1)
        assert len(ends) == 2, ends
        # start at the CH2 end (head); the other end is the benzylic carbon
        start = min(ends, key=lambda i: self.is_benzylic(i, backbone))
        order = [start]
        while len(order) < len(backbone):
            nxt = (set(nodes[order[-1]]["edges"]) & backbone) - set(order)
            assert len(nxt) == 1
            order.append(nxt.pop())
        return order

    def is_benzylic(self, index, backbone=None):
        backbone = backbone if backbone is not None else set(self.backbone)
        return any(self.nodes[j]["element"] == "C" and j not in backbone
                   for j in self.nodes[index]["edges"])

    def site(self, kind):
        """Backbone atom index for a radical site kind."""
        if kind == "end":  # primary CH2 at the head terminus
            return self.backbone[0]
        if kind == "benzylic_end":  # benzylic CH2Ph at the tail terminus
            return self.backbone[-1]
        if kind == "interior":  # tertiary benzylic carbon at mid chain
            benzylic = [i for i in self.backbone if self.is_benzylic(i)]
            return benzylic[len(benzylic) // 2]
        raise ValueError(kind)

    def abstract(self, carbon):
        hydrogen = max(j for j in self.nodes[carbon]["edges"]
                       if self.nodes[j]["element"] == "H")
        del self.nodes[carbon]["edges"][hydrogen]
        del self.nodes[hydrogen]
        self.nodes[carbon]["radical"] = 1
        self.removed.append(hydrogen)

    def position_of(self, index):
        position = {atom: p for p, atom in enumerate(self.backbone)}
        if index in position:
            return position[index]
        seen, frontier = {index}, {index}
        while True:
            hit = frontier & set(position)
            if hit:
                return min(position[a] for a in hit)
            frontier = {j for a in frontier for j in self.nodes[a]["edges"]} - seen
            assert frontier
            seen |= frontier


def materialize(chain, repeat_unit_formula=None, labeler=None):
    """ChainGraph -> state.Strand (uuids default to the deterministic name:index)."""
    from rmgpy.kmc.state import AtomRef, Strand

    strand_id = chain.name
    refs, graph, formula = {}, {}, Counter()
    backbone = set(chain.backbone)
    uuids = {i: (labeler(chain.name, i) if labeler else f"{chain.name}:{i}")
             for i in chain.nodes}
    for key, (i, node) in enumerate(sorted(chain.nodes.items())):
        refs[key] = AtomRef(uuids[i], strand_id, chain.position_of(i),
                            "backbone" if i in backbone else "feature")
        graph[uuids[i]] = {
            "element": node["element"], "radical": int(node["radical"]),
            "charge": 0, "lone_pairs": int(node["lone_pairs"]),
            "implicit_hydrogens": 0,
            "edges": {uuids[j]: order for j, order in node["edges"].items()},
        }
        formula[node["element"]] += 1
    radicals = sum(n["radical"] for n in chain.nodes.values())
    return Strand(strand_id, len(chain.backbone), {}, dict(formula), radicals,
                  refs, graph, dict(repeat_unit_formula or {}))


def styrene_strand(name="styrene", labeler=None):
    """A declared mobile small species as a length-1 strand (no backbone)."""
    from rmgpy.kmc.state import AtomRef, Strand, _parse_adjacency

    text = """multiplicity 1
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,S} {6,S}
3 C u0 p0 c0 {1,S} {7,B} {11,B}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 C u0 p0 c0 {3,B} {8,B} {12,S}
8 C u0 p0 c0 {7,B} {9,B} {13,S}
9 C u0 p0 c0 {8,B} {10,B} {14,S}
10 C u0 p0 c0 {9,B} {11,B} {15,S}
11 C u0 p0 c0 {10,B} {3,B} {16,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
"""
    nodes = _parse_adjacency(text)
    refs, graph, formula = {}, {}, Counter()
    for i, node in nodes.items():
        ident = labeler(name, i) if labeler else f"{name}:{i}"
        refs[i] = AtomRef(ident, name, 0, "feature")
        graph[ident] = {**{k: node[k] for k in ("element", "radical", "charge",
                                                "lone_pairs", "implicit_hydrogens")},
                        "edges": {(labeler(name, j) if labeler else f"{name}:{j}"): o
                                  for j, o in node["edges"].items()}}
        formula[node["element"]] += 1
    return Strand(name, 1, {}, dict(formula), 0, refs, graph)


def seeded_labeler(seed):
    """uuid4-shaped labels from a seeded stream (what uuid.uuid4 would give)."""
    import random

    rng = random.Random(seed)
    issued = {}

    def label(name, index):
        key = (name, index)
        if key not in issued:
            issued[key] = str(uuid.UUID(int=rng.getrandbits(128), version=4))
        return issued[key]

    return label


def build_state(chains, with_styrene=True, repeat_unit_formula=None, labeler=None):
    from rmgpy.kmc.state import KMCState

    strands = [materialize(c, repeat_unit_formula, labeler) for c in chains]
    if with_styrene:
        strands.append(styrene_strand(labeler=labeler))
    return KMCState(strands)


def scenario_chains(kind):
    """Scenario -> list of ChainGraph (the target chain first, a closed partner second)."""
    targets = {
        "closed": (),
        "interior": ("interior",),
        "end": ("end",),
        "benzylic_end": ("benzylic_end",),
        "interior+end": ("interior", "end"),
    }[kind]
    return [ChainGraph(DP, "target", targets), ChainGraph(DP, "partner")]


SCENARIOS = ("closed", "interior", "end", "benzylic_end", "interior+end")


# --------------------------------------------------------------------------
# C1: long-chain matching
# --------------------------------------------------------------------------
def _nx_chain(state):
    """networkx view of every atom graph in the state (uuid -> node attrs)."""
    graph = nx.Graph()
    for strand in state.strands.values():
        for ident, node in strand.atom_graph.items():
            graph.add_node(ident, element=node["element"], radical=node["radical"],
                           degree=len(node["edges"]), strand=strand.strand_id)
    for strand in state.strands.values():
        for ident, node in strand.atom_graph.items():
            for other, order in node["edges"].items():
                graph.add_edge(ident, other, order=float(order))
    return graph


def expected_sites(record, state_graph):
    """Chemically valid sites: radius-1 reaction-centre environments.

    Independent of SiteIndex: per reactant graph, the atoms the record edits
    (bond_ops / set_* atoms) are the centre; the pattern is the centre plus its
    neighbours (leaf atoms, unconstrained beyond element/radical), centre atoms
    must have exactly the pattern degree.  Returns, per participant, the set of
    distinct centre images (tuples of uuids) found in the state.
    """
    from rmgpy.kmc.state import _parse_adjacency

    graphs = [_parse_adjacency(text) for text in record.reactant_graphs]
    offsets = list(itertools.accumulate([0] + [len(g) for g in graphs]))
    edited = set()
    for op in record.bond_ops or ():
        if "atoms" in op:
            edited.update(op["atoms"])
        if "atom" in op:
            edited.add(op["atom"])
    per_participant = []
    for graph, offset in zip(graphs, offsets):
        centre = sorted(i for i in graph if (offset + i) in edited)
        pattern = nx.Graph()
        for i in centre:
            pattern.add_node(("c", i), element=graph[i]["element"],
                             radical=graph[i]["radical"], degree=len(graph[i]["edges"]),
                             centre=True)
        for i in centre:
            for j, order in graph[i]["edges"].items():
                key = ("c", j) if j in centre else ("n", i, j)
                if key not in pattern:
                    pattern.add_node(key, element=graph[j]["element"],
                                     radical=graph[j]["radical"], degree=None,
                                     centre=False)
                pattern.add_edge(("c", i), key, order=float(order))

        def node_match(target, wanted):
            if target["element"] != wanted["element"] or target["radical"] != wanted["radical"]:
                return False
            return not wanted["centre"] or target["degree"] == wanted["degree"]

        def edge_match(target, wanted):
            return target["order"] == wanted["order"]

        matcher = nx.algorithms.isomorphism.GraphMatcher(
            state_graph, pattern, node_match=node_match, edge_match=edge_match)
        images = set()
        for mapping in matcher.subgraph_monomorphisms_iter():
            inverse = {v: k for k, v in mapping.items()}
            images.add(tuple(inverse[("c", i)] for i in centre))
        per_participant.append(images)
    return per_participant


def location_class(site_graph, length):
    """'end' when the matched backbone atoms touch a terminus, else 'interior'."""
    positions = [node["position"] for node in site_graph.values()
                 if node.get("role") == "backbone"]
    if not positions:
        return "no-backbone"
    if min(positions) <= 1 or max(positions) >= length - 2:
        return "end"
    return "interior"


def c1_scenario(kind, records):
    """Per-record SiteIndex candidate counts split by location, plus expected."""
    from rmgpy.kmc.ssa import SiteIndex

    state = build_state(scenario_chains(kind))
    index = SiteIndex(state, records)
    state_graph = _nx_chain(state)
    out = {}
    for record in records:
        candidates = index.candidates(record.event_id)
        classes = Counter()
        for combination in candidates:
            for site in combination:
                if site.strand_id in ("target", "partner"):
                    length = state.strands[site.strand_id].length
                    classes[(site.strand_id, location_class(site.graph, length))] += 1
        out[record.event_id] = {"n": len(candidates), "classes": classes}
    return state, index, out, state_graph


def build_state_new_atom(chains, with_styrene=True):
    """Build through the only production uuid source, ``KMCState.new_atom`` (uuid4)."""
    from rmgpy.kmc.state import KMCState, Strand

    strands = []
    for chain in chains:
        formula = Counter(n["element"] for n in chain.nodes.values())
        radicals = sum(n["radical"] for n in chain.nodes.values())
        strands.append(Strand(chain.name, len(chain.backbone), {}, dict(formula), radicals))
    if with_styrene:
        strands.append(Strand("styrene", 1, {}, dict(Counter("C" if i < 8 else "H" for i in range(16))), 0))
    state = KMCState(strands)
    backbone_of = {c.name: set(c.backbone) for c in chains}
    for chain in chains:
        mapping = {}
        for i in sorted(chain.nodes):
            role = "backbone" if i in backbone_of[chain.name] else "feature"
            mapping[i] = state.new_atom(chain.name, chain.position_of(i), role).uuid
        strand = state.strands[chain.name]
        for i, node in chain.nodes.items():
            strand.atom_graph[mapping[i]] = {
                "element": node["element"], "radical": int(node["radical"]),
                "charge": 0, "lone_pairs": int(node["lone_pairs"]),
                "implicit_hydrogens": 0,
                "edges": {mapping[j]: o for j, o in node["edges"].items()}}
    if with_styrene:
        template = styrene_strand()
        ident = {}
        for key, ref in sorted(template.atom_refs.items()):
            ident[ref.uuid] = state.new_atom("styrene", 0, "feature").uuid
        for old, node in template.atom_graph.items():
            state.strands["styrene"].atom_graph[ident[old]] = {
                **{k: v for k, v in node.items() if k != "edges"},
                "edges": {ident[j]: o for j, o in node["edges"].items()}}
        state.strands["styrene"].formula = dict(template.formula)
    state._refresh_components()
    state.assert_uuid_uniqueness()
    return state


def end_vs_interior(state, uuid_tuple, window=3):
    """Classify a centre image by its backbone positions (end = within ``window``)."""
    refs = [state._ref(u) for u in uuid_tuple]
    length = state.strands[refs[0].strand_id].length
    positions = [r.position for r in refs]
    return "end" if min(positions) < window or max(positions) >= length - window else "interior"


def c1_report(records, scenarios=SCENARIOS, with_expected=True):
    """Per scenario and family: SiteIndex matches vs radius-1 expected sites."""
    rows = {}
    for kind in scenarios:
        state, index, out, state_graph = c1_scenario(kind, records)
        table = {}
        for record in records:
            key = (record.family, record.arity)
            row = table.setdefault(key, Counter())
            row["records"] += 1
            matched = out[record.event_id]["n"]
            row["siteindex_candidates"] += matched
            row["siteindex_records_matching"] += bool(matched)
            for (strand, cls), count in out[record.event_id]["classes"].items():
                row[f"siteindex_{cls}"] += count
            if with_expected:
                images = expected_sites(record, state_graph)
                if record.arity == 1:
                    chosen = images[0]
                    row["expected_sites"] += len(chosen)
                    row["expected_records_with_sites"] += bool(chosen)
                    for image in chosen:
                        row[f"expected_{end_vs_interior(state, image)}"] += 1
                else:
                    bound = min(len(i) for i in images)
                    row["expected_pair_bound"] += bound
                    row["expected_records_with_sites"] += bool(bound)
        rows[kind] = table
    return rows


def cap_diagnostic(records):
    """Why the DP3 graph cannot embed: its cap carbons need terminal hydrogens.

    For the pristine homolysis record and the chain, the hydrogen count on the
    backbone carbons that carry a cap (degree-1 in the backbone) versus every
    backbone carbon of the DP20 chain by position class.
    """
    from rmgpy.kmc.state import _parse_adjacency

    record = next(r for r in records if r.family == "R_Recombination"
                  and r.arity == 1 and r.participant_site_types == ["pristine"])
    graph = _parse_adjacency(record.reactant_graphs[0])
    backbone = {i for i, n in graph.items() if n["element"] == "C"
                and all(float(o) == 1.0 for o in n["edges"].values())}
    caps = {}
    for i in backbone:
        if len(set(graph[i]["edges"]) & backbone) == 1:
            caps[i] = sum(graph[j]["element"] == "H" for j in graph[i]["edges"])
    chain = ChainGraph(DP, "diag")
    per_position = {}
    for p, atom in enumerate(chain.backbone):
        per_position[p] = sum(chain.nodes[j]["element"] == "H" for j in chain.nodes[atom]["edges"])
    interior = Counter(per_position[p] for p in range(2, len(chain.backbone) - 2))
    return {"record_cap_carbon_H_counts": sorted(caps.values()),
            "chain_terminal_carbon_H_counts": [per_position[0], per_position[len(chain.backbone) - 1]],
            "chain_interior_backbone_H_counts": dict(interior),
            "record_graph_atoms": len(graph)}


def h_census(chain):
    """Hydrogens of a closed chain by environment (the H_Abstraction expectation)."""
    backbone = set(chain.backbone)
    census = Counter()
    for atom, node in chain.nodes.items():
        if node["element"] != "H":
            continue
        heavy = next(iter(node["edges"]))
        if heavy not in backbone:
            census["aromatic"] += 1
            continue
        neighbours = len(set(chain.nodes[heavy]["edges"]) & backbone)
        if neighbours == 1:
            census["terminal carbon"] += 1
        elif chain.is_benzylic(heavy):
            census["tertiary benzylic"] += 1
        else:
            census["secondary backbone"] += 1
    return dict(census)


# --------------------------------------------------------------------------
# C2: scission ownership
# --------------------------------------------------------------------------
def c2_dp3_control(records, limit=6):
    """Backbone-cutting events DO apply on a DP3 chain (window = whole strand)."""
    from rmgpy.kmc.ssa import SiteIndex

    chains = [ChainGraph(3, "target", ("end",)), ChainGraph(3, "partner", ("benzylic_end",))]
    state = build_state(chains)
    index = SiteIndex(state, records)
    outcomes = []
    tried = Counter()
    for record in records:
        candidates = index.candidates(record.event_id)
        if not candidates:
            continue
        for ordinal in range(min(len(candidates), 4)):
            copy_state = copy.deepcopy(state)
            copy_index = SiteIndex(copy_state, [record])
            site = copy_index.candidates(record.event_id)[ordinal]
            try:
                entry = copy_index.apply(record, site)
                if not entry.derived_cut_offsets:
                    tried["applied, no backbone cut"] += 1
                    continue
                lengths = sorted(s.length for s in copy_state.strands.values())
                outcomes.append((record.event_id, record.family, "ok",
                                 entry.derived_cut_offsets, lengths))
                tried["applied, backbone cut"] += 1
            except Exception as error:  # noqa: BLE001 - the probe reports every failure
                tried[f"{type(error).__name__}"] += 1
        if len(outcomes) >= limit:
            break
    return outcomes, dict(tried)


def _window(chain_state, strand_id, first, second, radius=3):
    """Atoms within ``radius`` bonds of the two cut atoms (a reaction window)."""
    graph = chain_state.strands[strand_id].atom_graph
    frontier = {first, second}
    seen = set(frontier)
    for _ in range(radius):
        frontier = {j for a in frontier for j in graph[a]["edges"]} - seen
        seen |= frontier
    return seen


def c2_dp20_split(window_radius=3):
    """Call ``_split_strand`` exactly as ``apply`` does, on a DP20 radical chain.

    Returns {variant: outcome}; ``apply`` passes the mapped window's elements,
    so the window is the atoms within ``window_radius`` of the cut bond.
    """
    from rmgpy.kmc.state import LedgerError

    results = {}
    variants = {
        "repeat_unit_formula={} (as the DP3 melt builder)": {},
        "repeat_unit_formula=C4H4 per backbone position": {"C": 4, "H": 4},
    }
    for label, unit in variants.items():
        for scope in ("window", "whole-strand"):
            chain = ChainGraph(DP, "t", ("interior", "end"))
            state = build_state([chain], with_styrene=False, repeat_unit_formula=unit)
            strand = state.strands["t"]
            by_position = {r.position: r.uuid for r in strand.atom_refs.values() if r.role == "backbone"}
            cut = len(chain.backbone) // 2 - 1
            first, second = by_position[cut], by_position[cut + 1]
            members = (_window(state, "t", first, second, window_radius) if scope == "window"
                       else set(strand.atom_graph))
            elements = {u: strand.atom_graph[u]["element"] for u in members}
            before = dict(strand.formula)
            try:
                right_id, offset = state._split_strand(first, second, elements)
                left, right = state.strands["t"], state.strands[right_id]
                outcome = ("ok", dict(left.formula), dict(right.formula),
                           _graph_formula(left), _graph_formula(right))
            except LedgerError as error:
                outcome = ("LedgerError", str(error))
            results[(label, scope)] = {"before": before, "mapped_atoms": len(members), "outcome": outcome}
    return results


def _graph_formula(strand):
    return dict(sorted(Counter(n["element"] for n in strand.atom_graph.values()).items()))


# --------------------------------------------------------------------------
# C3: initiation on a closed melt
# --------------------------------------------------------------------------
def radical_free_inventory(records, temperature=700.0):
    """Records whose reactants carry no radical: the only initiation candidates."""
    from rmgpy.kmc.ssa import record_rate
    from rmgpy.kmc.state import _parse_adjacency

    table = Counter()
    rows = []
    for record in records:
        graphs = [_parse_adjacency(text) for text in record.reactant_graphs]
        if any(n["radical"] for g in graphs for n in g.values()):
            continue
        rate = record_rate(record, temperature)
        table[(record.family, record.arity, record.orientation)] += 1
        rows.append((record.event_id, record.family, record.arity, tuple(record.participant_site_types), rate))
    return table, rows


def closed_melt_chains(counts, prefix="m"):
    chains = []
    for dp, number in sorted(counts.items()):
        for k in range(number):
            chains.append(ChainGraph(dp, f"{prefix}{dp}-{k}"))
    return chains


def c3_run(records, counts, temperature=700.0, volume=1e-24, max_events=5, seed=1):
    from rmgpy.kmc.ssa import IsothermalSSA

    state = build_state(closed_melt_chains(counts))
    engine = IsothermalSSA(state, records, temperature=temperature, volume=volume,
                           rng=np.random.default_rng(seed),
                           inventory="R1")
    report = engine.channel_propensities()
    before = {"chains": len(state.strands) - 1, "radicals": state.total_radicals}
    positive = {k: v for k, v in report.propensities.items() if v > 0}
    active = sum(bool(engine.index.candidates(r.event_id)) for r in records)
    result = engine.run(max_events=max_events)
    return {**before, "records_with_candidates": active, "positive_channels": len(positive),
            "total_enabled": report.total_enabled, "events_fired_of_requested_5": len(result.events),
            "engine_time_s": result.time}


# --------------------------------------------------------------------------
# C4: coproduct atoms
# --------------------------------------------------------------------------
def c4_coproduct(records, formula_wanted, clear_reverse_of):
    """Apply one coproduct-producing event on a DP3 radical melt and inspect it."""
    from rmgpy.kmc.event_record import EventRecord
    from rmgpy.kmc.ssa import SiteIndex

    chains = [ChainGraph(3, "target", ("end",)), ChainGraph(3, "partner", ("benzylic_end",)),
              ChainGraph(3, "third")]
    state = build_state(chains, with_styrene=False)
    pool = [r for r in records if r.coproducts and r.coproducts[0].get("formula") == formula_wanted
            and r.arity == 1]
    index = SiteIndex(state, pool)
    record = next(r for r in pool if index.candidates(r.event_id))
    if clear_reverse_of:
        data = record.to_dict()
        data["reverse_of"] = ""
        applied_record = EventRecord.from_dict(data)
        # the index keys candidates by event id, so rebuild it for the edited record
        index = SiteIndex(state, [applied_record])
    else:
        applied_record = record
    site = index.candidates(applied_record.event_id)[0]
    graph_before = {s.strand_id: _graph_formula(s) for s in state.strands.values()}
    entry = index.apply(applied_record, site)
    wanted = Counter(formula_wanted)
    coproduct_atoms = None
    for component in entry.product_graphs:
        if Counter(row[1] for row in component) == wanted:
            coproduct_atoms = {row[0] for row in component}
            break
    in_graph = {u for s in state.strands.values() for u in s.atom_graph}
    # all candidate sites of the full record set after the event, to see whether
    # coproduct atoms can still be matched into a site
    full = SiteIndex(state, records)
    touching = sum(
        1 for r in records for combo in full.candidates(r.event_id)
        for site in combo for node in site.graph.values()
        if node["atom_ref"].uuid in (coproduct_atoms or ()))
    mismatch = {s.strand_id: (dict(s.formula), _graph_formula(s)) for s in state.strands.values()
                if dict(s.formula) != _graph_formula(s)}
    return {
        "record": applied_record.event_id[:12], "family": applied_record.family,
        "reverse_of_set": bool(applied_record.reverse_of), "coproduct_formula": formula_wanted,
        "sink": dict(state.sink), "sink_radicals": state.sink_radicals,
        "coproduct_atoms": len(coproduct_atoms or ()),
        "coproduct_atoms_still_in_a_strand_graph": len((coproduct_atoms or set()) & in_graph),
        "sites_touching_coproduct_atoms": touching,
        "strands_where_formula_ne_graph": mismatch,
        "total_ledger": state.total_ledger,
    }


# --------------------------------------------------------------------------
# C5: moment basis
# --------------------------------------------------------------------------
def c5_report(counts=None):
    counts = counts or {10: 3, 20: 4, 30: 3}
    spec = ps_spec(max(counts), counts)
    realized = spec.initial_ensemble.counts
    chains = closed_melt_chains(realized)
    state = build_state(chains, with_styrene=False)
    dp = []
    for strand in state.strands.values():
        dp.append(strand.formula.get("C", 0) // 8)  # independent: carbons / C8 repeat unit
    mu = (len(dp), sum(dp), sum(d * d for d in dp))
    declared = (sum(realized.values()), sum(k * v for k, v in realized.items()),
                sum(k * k * v for k, v in realized.items()))
    masses = [s.mass for s in state.strands.values()]
    return {"declared_realized_counts": realized, "declared_dp_moments": declared,
            "independent_dp_moments_over_strands": mu, "state_moments": state.moments,
            "state_mass_moments_recomputed": (len(masses), sum(masses), sum(m * m for m in masses)),
            "strand_lengths": sorted(s.length for s in state.strands.values()),
            "spec_macroscopic_moments_mol": spec.macroscopic_moments_mol,
            "spec_mn_g_mol": spec.ensemble_observables.mn_g_mol}


# --------------------------------------------------------------------------
# C6 / C7: trajectories
# --------------------------------------------------------------------------
def _melt_for_trajectory(labeler=None, new_atom_path=False):
    chains = [ChainGraph(3, f"c{i}") for i in range(4)] + [
        ChainGraph(3, "r0", ("end",)), ChainGraph(3, "r1", ("benzylic_end",))]
    if new_atom_path:
        return build_state_new_atom(chains)
    return build_state(chains, labeler=labeler)


def _signature(state):
    strands = sorted((s.length, tuple(sorted(s.formula.items())), s.radical_count)
                     for s in state.strands.values())
    return {"strands": strands, "sink": sorted(state.sink.items()), "moments": list(state.moments)}


def run_trajectory(records, mode, seed, schedule="single", horizon=5.0e-5, label_seed=None):
    """Event sequence of the DP3 control melt.  mode: named | uuid4 | seeded."""
    from rmgpy.kmc.ssa import IsothermalSSA

    if mode == "uuid4":
        state = _melt_for_trajectory(new_atom_path=True)
    elif mode == "seeded":
        state = _melt_for_trajectory(labeler=seeded_labeler(label_seed))
    else:
        state = _melt_for_trajectory()
    engine = IsothermalSSA(state, records, temperature=700.0, volume=1e-24,
                           rng=np.random.default_rng(seed),
                           inventory="R1")
    if schedule == "single":
        engine.run(until=horizon)
    else:
        for i in range(1, 11):
            engine.run(until=horizon * i / 10)
    return {"events": [(repr(e.time), e.event_id[4:12]) for e in engine.events],
            "applied": [(list(a.participants), list(a.derived_cut_offsets)) for a in state.event_log],
            "time": repr(engine.time), "signature": _signature(state)}


def c6_report(records, seeds=(1, 2, 3)):
    out = []
    for seed in seeds:
        single = run_trajectory(records, "named", seed, "single")
        piecewise = run_trajectory(records, "named", seed, "ten")
        same_ids = [e[1] for e in single["events"]] == [e[1] for e in piecewise["events"]]
        out.append({"seed": seed, "events_single": len(single["events"]),
                    "events_ten_calls": len(piecewise["events"]),
                    "same_event_id_sequence": same_ids,
                    "same_times": single["events"] == piecewise["events"],
                    "same_final_signature": single["signature"] == piecewise["signature"]})
    return out


def trace_in_subprocess(mode, seed, extra=()):
    environment = dict(os.environ)
    environment.update(PYTHONHASHSEED="0", PYTHONPATH=str(REPO_ROOT),
                       MPLCONFIGDIR=environment.get("MPLCONFIGDIR", "/tmp/i125-mpl"))
    completed = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), "trace", mode, str(seed), *extra],
        cwd=REPO_ROOT, env=environment, capture_output=True, text=True, check=True, timeout=600)
    return json.loads(completed.stdout.strip().splitlines()[-1])


def c7_report(records, seed=1):
    first = trace_in_subprocess("uuid4", seed)
    second = trace_in_subprocess("uuid4", seed)
    named_a = trace_in_subprocess("named", seed)
    named_b = trace_in_subprocess("named", seed)
    relabel = [run_trajectory(records, "seeded", seed, label_seed=k) for k in range(1, 6)]
    base = relabel[0]

    def ids(trace):
        return [e[1] for e in trace["events"]]

    return {
        "uuid4_process_a_events": len(first["events"]), "uuid4_process_b_events": len(second["events"]),
        "uuid4_same_event_ids": ids(first) == ids(second),
        "uuid4_same_times": first["events"] == second["events"],
        "uuid4_same_final_signature": first["signature"] == second["signature"],
        "uuid4_same_applied_sites": first["applied"] == second["applied"],
        "named_same_trace": named_a == named_b,
        "relabel_same_event_ids": [ids(r) == ids(base) for r in relabel[1:]],
        "relabel_same_times": [r["events"] == base["events"] for r in relabel[1:]],
        "relabel_same_final_signature": [r["signature"] == base["signature"] for r in relabel[1:]],
        "relabel_same_applied_sites": [r["applied"] == base["applied"] for r in relabel[1:]],
    }


# --------------------------------------------------------------------------
# command line
# --------------------------------------------------------------------------
def _print(title, value):
    print(f"--- {title}")
    if isinstance(value, dict):
        for key, item in value.items():
            print(f"  {key}: {item}")
    elif isinstance(value, (list, tuple)):
        for item in value:
            print(f"  {item}")
    else:
        print(f"  {value}")


def report(which):
    records = load_records()
    print(f"artifact records: {len(records)}; families {sorted({r.family for r in records})}")
    if "C1" in which:
        _print("C1 cap diagnostic", cap_diagnostic(records))
        _print("C1 H census of the closed DP20 chain", h_census(ChainGraph(DP, "census")))
        for dp, rad in ((3, ()), (3, ("end",)), (4, ("end",)), (5, ("end",))):
            state = build_state([ChainGraph(dp, "t", rad), ChainGraph(dp, "p")])
            from rmgpy.kmc.ssa import SiteIndex
            index = SiteIndex(state, records)
            n = sum(len(index.candidates(r.event_id)) for r in records)
            m = sum(bool(index.candidates(r.event_id)) for r in records)
            print(f"  control DP{dp} radicals={rad}: records matching={m} candidates={n}")
        for kind, table in c1_report(records).items():
            print(f"  scenario {kind}")
            for key, row in sorted(table.items()):
                print(f"    {key}: {dict(sorted(row.items()))}")
    if "C2" in which:
        outcomes, tried = c2_dp3_control(records)
        _print("C2 DP3 control (backbone-cut events)", {"tried": tried, "first": outcomes[:3]})
        _print("C2 DP20 _split_strand", {str(k): v for k, v in c2_dp20_split().items()})
    if "C3" in which:
        table, rows = radical_free_inventory(records)
        _print("C3 radical-free-substrate records", dict(table))
        _print("C3 DP3 closed melt", c3_run(records, {3: 5}))
        _print("C3 DP20 closed melt", c3_run(records, {10: 2, 20: 2, 30: 1}))
    if "C4" in which:
        for formula in ({"H": 1}, {"C": 6, "H": 5}):
            for clear in (False, True):
                _print(f"C4 coproduct {formula} reverse_of cleared={clear}",
                       c4_coproduct(records, formula, clear))
    if "C5" in which:
        _print("C5 moment basis", c5_report())
    if "C6" in which:
        _print("C6 snapshot schedule", c6_report(records))
    if "C7" in which:
        _print("C7 uuid4", c7_report(records))


def main(argv):
    command = argv[1] if len(argv) > 1 else "report"
    if command == "build":
        _build_main(argv[2], argv[3])
    elif command == "trace":
        records = load_records()
        result = run_trajectory(records, argv[2], int(argv[3]),
                                label_seed=int(argv[4]) if len(argv) > 4 else None)
        print(json.dumps(result))
    elif command == "report":
        report(set(argv[2:]) or {"C1", "C2", "C3", "C4", "C5", "C6", "C7"})
    else:
        raise SystemExit(f"unknown command {command}")


if __name__ == "__main__":
    main(sys.argv)
