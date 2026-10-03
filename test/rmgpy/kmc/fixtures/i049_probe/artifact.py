"""Run: PYTHONPATH=$PWD python test/rmgpy/kmc/fixtures/i049_probe/artifact.py --output DIR

One JSON load of the dispatch artifact; cache only distinct molecule graphs.
Write a small record subset for runtime.py so it never reloads the artifact.
"""

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import resource
import time

from common import benzylic_tail, canonical_smiles, describe, terminal_benzylic

ARTIFACT = Path("/home/alon/runs/i046-rules-from-training/cache/artifact/"
                "7491ed418f3ce2f0633688a9487cd12da278704b8f57e4f595176ddf80fb0296.json")


def brief(record):
    return {key: record[key] for key in ("event_id", "family", "site_type", "template",
                                        "arity", "orientation", "reverse_of")} | {
        "reactants": [describe(g)["smiles"] for g in record["reactant_graphs"]],
        "products": [describe(g)["smiles"] for g in record["product_graphs"]]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    start = time.monotonic()
    with ARTIFACT.open() as stream:
        data = json.load(stream)
    assert data["provenance"]["rmg_database_sha"] == "4a12d36fcdc193ede82c8d1ab5c1653495d445bc"
    records = data["records"]
    assert len(records) == 14998
    by_id = {r["event_id"]: r for r in records}
    targets = {n: canonical_smiles(benzylic_tail(n)) for n in (1, 2, 3, 5)}
    reaction_counts, site_counts = Counter(), Counter()
    producers, consumers = Counter(), Counter()
    standard = {n: {"reactant_occurrences": 0, "product_occurrences": 0} for n in targets}
    examples = defaultdict(list)
    selected = {}
    benz_styrene = []
    for r in records:
        reaction_counts[r["family"]] += 1
        site_counts[r["site_type"].split("@")[0]] += 1
        rb = any(terminal_benzylic(g) for g in r["reactant_graphs"])
        pb = any(terminal_benzylic(g) for g in r["product_graphs"])
        if rb:
            consumers[r["family"]] += 1
        if pb:
            producers[r["family"]] += 1
            if len(examples[r["family"]]) < 2:
                examples[r["family"]].append(brief(r))
            # Prefer the smallest example, retaining all fields for real execution.
            key = r["family"]
            eligible = True
            if key == "R_Addition_MultipleBond":
                key += ":addition" if r["arity"] == 2 else ":scission"
            if r["family"] == "R_Recombination":
                # Probe backbone homolysis rather than H-atom loss from a ring.
                eligible = len(r["product_graphs"]) == 2 and all(
                    describe(g)["formula"].get("C", 0) >= 2 for g in r["product_graphs"])
            previous = selected.get(key)
            if eligible and (previous is None or sum(map(len, r["reactant_graphs"])) < sum(map(len, previous["reactant_graphs"]))):
                selected[key] = r
        for n, target in targets.items():
            standard[n]["reactant_occurrences"] += sum(describe(g)["smiles"] == target for g in r["reactant_graphs"])
            standard[n]["product_occurrences"] += sum(describe(g)["smiles"] == target for g in r["product_graphs"])
        if r["family"] == "R_Addition_MultipleBond":
            rs = any(describe(g)["smiles"] == "C=Cc1ccccc1" for g in r["reactant_graphs"])
            ps = any(describe(g)["smiles"] == "C=Cc1ccccc1" for g in r["product_graphs"])
            if (rb or pb) and (rs or ps):
                benz_styrene.append(brief(r))
    head_to_tail = [r for r in records if r["family"] == "R_Addition_MultipleBond"
                    and r["arity"] == 2
                    and any(describe(g)["smiles"] == "C=Cc1ccccc1" for g in r["reactant_graphs"])
                    and any(terminal_benzylic(g) for g in r["reactant_graphs"])]
    assert not head_to_tail
    assert set(selected) == {"R_Recombination", "intra_H_migration",
                             "R_Addition_MultipleBond:addition", "R_Addition_MultipleBond:scission"}
    ceiling = []
    for pair in data["ps_ceiling_pairs"]:
        prop = by_id[pair["propagation_event_id"]]
        dep = by_id[pair["depropagation_event_id"]]
        radical_reactant = next(describe(g)["radicals"] for g in prop["reactant_graphs"] if describe(g)["radicals"])
        radical_product = describe(prop["product_graphs"][0])["radicals"]
        assert len(radical_reactant) == len(radical_product) == 1
        assert radical_reactant[0]["H"] == radical_product[0]["H"] == 2
        assert not radical_reactant[0]["benzylic"] and not radical_product[0]["benzylic"]
        assert dep["reverse_of"] == prop["event_id"]
        ceiling.append({"pair": pair, "propagation": brief(prop), "depropagation": brief(dep)})
    result = {"artifact": str(ARTIFACT), "bytes": ARTIFACT.stat().st_size,
              "records": len(records), "discovery": len(data["discovery"]),
              "families": dict(reaction_counts), "site_types": dict(site_counts),
              "by_context": dict(Counter(r["site_type"] for r in records)),
              "terminal_benzylic_producers": dict(producers), "terminal_benzylic_consumers": dict(consumers),
              "standard_tail_occurrences": standard, "producer_examples": dict(examples),
              "terminal_benzylic_and_styrene": benz_styrene, "head_to_tail_propagations": 0,
              "ceiling": ceiling,
              "discovery_by_context": dict(Counter(d["site_type"] for d in data["discovery"])),
              "elapsed_s": time.monotonic() - start,
              "max_rss_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    subset = {"producers": selected,
              "inverse_records": {family: by_id[r["reverse_of"]] for family, r in selected.items()}}
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "artifact-results.json").write_text(json.dumps(result, indent=2) + "\n")
    (args.output / "runtime-records.json").write_text(json.dumps(subset, indent=2) + "\n")
    for key in ("records", "discovery", "families", "site_types", "terminal_benzylic_producers",
                "terminal_benzylic_consumers", "standard_tail_occurrences", "head_to_tail_propagations"):
        print(key, result[key], flush=True)
    for row in benz_styrene:
        print("benzylic+styrene", row)
    print("PASS: zero benzylic-end+styrene additions; existing ceiling pairs use primary ends")


if __name__ == "__main__":
    main()
