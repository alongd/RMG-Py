"""Run: PYTHONPATH=$PWD python test/rmgpy/kmc/fixtures/i049_probe/runtime.py DIR

Seed true backbone positions, then execute the small subset written by artifact.py.
No artifact load, RMG family generation, or event compilation.
"""

from collections import Counter, deque
import copy
import json
from pathlib import Path
import sys

if __package__:
    from .common import benzylic_tail, describe, molecule
else:
    from common import benzylic_tail, describe, molecule
from rmgpy.kmc.ssa import SiteIndex, radical_site_class
from rmgpy.kmc.state import AtomRef, KMCState, Strand, _parse_adjacency
from rmgpy.molecule.molecule import Molecule


def seed(record):
    strands = []
    for ordinal, text in enumerate(record["reactant_graphs"]):
        graph = _parse_adjacency(text)
        mol = molecule(text)
        backbone = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and not a.IsInRing()}
        adjacency = {i: sorted(j for j in graph[i]["edges"] if j in backbone) for i in backbone}
        endpoints = sorted(i for i in backbone if len(adjacency[i]) <= 1)
        assert endpoints and all(len(v) <= 2 for v in adjacency.values())
        path, previous, current = [], None, endpoints[0]
        while current is not None:
            path.append(current)
            remaining = [i for i in adjacency[current] if i != previous]
            previous, current = current, remaining[0] if remaining else None
        assert set(path) == backbone
        positions = {i: position for position, i in enumerate(path)}
        for local in graph:
            if local in positions:
                continue
            seen, queue = {local}, deque([local])
            while queue:
                node = queue.popleft()
                if node in backbone:
                    positions[local] = positions[node]
                    break
                for neighbor in graph[node]["edges"]:
                    if neighbor not in seen:
                        seen.add(neighbor)
                        queue.append(neighbor)
            assert local in positions
        strand_id = f"chain-{ordinal}"
        refs = {i: AtomRef(f"atom-{ordinal}-{i}", strand_id, positions[i],
                           "backbone" if i in backbone else "feature") for i in graph}
        owned = {}
        features = {}
        for local, node in graph.items():
            owned[refs[local].uuid] = {**copy.deepcopy(node), "edges": {
                refs[other].uuid: order for other, order in node["edges"].items()}}
            if node["radical"]:
                features.setdefault(positions[local], {}).setdefault("atoms", {})[
                    refs[local].uuid] = {"radical": node["radical"]}
        strands.append(Strand(strand_id, len(path), features,
                              dict(Counter(n["element"] for n in graph.values())),
                              sum(n["radical"] for n in graph.values()), refs, owned))
    return KMCState(strands)


def snapshot(state):
    # Reconstruct chemistry from persistent state-owned nodes, independently of labels.
    nodes = {u: n for strand in state.strands.values() for u, n in strand.atom_graph.items()}
    components = []
    unseen = set(nodes)
    while unseen:
        root = min(unseen)
        members, queue = set(), [root]
        while queue:
            node = queue.pop()
            if node in members:
                continue
            members.add(node)
            queue.extend(set(nodes[node]["edges"]) - members)
        unseen -= members
        indices = {u: i + 1 for i, u in enumerate(sorted(members))}
        lines = []
        orders = {"1.0": "S", "2.0": "D", "3.0": "T", "1.5": "B"}
        for uuid in sorted(members):
            node = nodes[uuid]
            edges = " ".join(f"{{{indices[v]},{orders[str(float(o))]}}}" for v, o in node["edges"].items())
            lines.append(f"{indices[uuid]} {node['element']} u{node['radical']} p{node['lone_pairs']} c{node['charge']} {edges}")
        components.append(describe("\n".join(lines)))
    return {"components": components, "strand_lengths": sorted(s.length for s in state.strands.values()),
            "ledger_radicals": state.total_radicals,
            "graph_radicals": sum(n["radical"] for n in nodes.values()),
            "sink_radicals": state.sink_radicals}


def main():
    output = Path(sys.argv[1])
    data = json.loads((output / "runtime-records.json").read_text())
    result = {"seeds": [], "events": []}
    for units in (1, 2, 3):
        adj = Molecule(smiles=benzylic_tail(units)).to_adjacency_list(remove_h=False)
        record = {"event_id": f"seed-{units}", "arity": 1, "reactant_graphs": [adj],
                  "participant_site_types": ["benzylic_end_radical"], "status": "enabled"}
        state = seed(record)
        candidates = SiteIndex(state, [record]).candidates(record["event_id"])
        assert candidates and radical_site_class(candidates[0][0], state) == "end"
        snap = snapshot(state)
        assert snap["ledger_radicals"] == snap["graph_radicals"] == 1
        assert sum(r["terminal_benzylic"] for c in snap["components"] for r in c["radicals"]) == 1
        result["seeds"].append({"units": units, "class": "end", "matches": len(candidates), **snap})
    for family, record in data["producers"].items():
        state = seed(record)
        before = snapshot(state)
        candidates = SiteIndex(state, [record]).candidates(record["event_id"])
        assert candidates
        applied = state.apply(record, candidates[0])
        after = snapshot(state)
        assert after["graph_radicals"] == after["ledger_radicals"]
        assert after["ledger_radicals"] == before["ledger_radicals"] + record["radical_delta"]
        assert any(r["terminal_benzylic"] for c in after["components"] for r in c["radicals"])
        inverse = data["inverse_records"][family]
        inverse_candidates = SiteIndex(state, [inverse]).candidates(inverse["event_id"])
        assert inverse_candidates, family
        state.apply(inverse, inverse_candidates[0])
        restored = snapshot(state)
        assert sorted(c["smiles"] for c in restored["components"]) == sorted(c["smiles"] for c in before["components"])
        assert restored["ledger_radicals"] == before["ledger_radicals"]
        row = {"channel": family, "family": record["family"], "event_id": record["event_id"], "orientation": record["orientation"],
               "participant_site_types": record["participant_site_types"],
               "derived_cut_offsets": applied.derived_cut_offsets,
               "before": before, "after": after, "inverse_matches": len(inverse_candidates),
               "round_trip": True}
        result["events"].append(row)
        print(json.dumps(row, indent=2), flush=True)
    (output / "runtime-results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"PASS: three benzylic-end seeds represented/indexed; {len(result['events'])} real producer/inverse round trips")


if __name__ == "__main__":
    main()
