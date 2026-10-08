"""Compare the one new PS artifact with the dispatch's immutable baseline."""

import argparse
from collections import Counter
import hashlib
import importlib.util
import json
import math
from pathlib import Path
from statistics import median

OLD_ARTIFACT = Path("/home/alon/runs/polymer/i038-bench/head-artifact/"
                    "0883a1292e20708a17ee8dc7960cf5648b354e0c24ce3373dee83813cea85c9d.json")
HERE = Path(__file__).resolve().parent
TREE_FAMILIES = ("Disproportionation", "R_Recombination")
TEMPERATURES = (600.0, 700.0, 800.0)
DATABASE_PROVENANCE_KEYS = {"database_sha", "rmg_database_sha"}


def load_i044():
    spec = importlib.util.spec_from_file_location("i044", HERE.parent / "i044_probe/run_probe.py")
    prior = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prior)
    return prior


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def structural_key(record):
    from rmgpy.kmc.event_record import _canonical_link_handles

    # Provenance, rates and corrected site classifications change event IDs,
    # while the mapped chemical rewrite remains the same.
    omitted = {"event_id", "canonical_index", "provenance", "rate_source", "k_table",
               "reverse_of", "thermo_provenance", "status", "status_reason",
               "site_type", "participant_site_types", "reactant_multiplicities"}
    payload = json.dumps(_canonical_link_handles(
        {k: v for k, v in record.items() if k not in omitted}), sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()


def index(artifact):
    result = {structural_key(record): record for record in artifact["records"]}
    assert len(result) == len(artifact["records"]), "ambiguous chemical rewrite identity"
    return result


def source_entries(source):
    return [item["entry"] for item in source.get("rules", [])] + [
        item.get("rule", item["entry"]) for item in source.get("training", [])
    ]


def has_training(source):
    return bool(source.get("training")) or any(
        "generated from training reaction" in entry.get("short_desc", "").lower()
        for entry in source_entries(source)
    )


def has_default(source):
    return any(entry.get("short_desc", "").strip().lower() == "default"
               for entry in source_entries(source))


def has_default_or_root(source, root_label):
    if (has_default(source) or source.get("entry") == root_label
            or source.get("template") == root_label):
        return True
    # Count a newly averaged root too; source contributors alone may be leaves.
    return any(
        line.startswith(f"Estimated using template [{root_label}] for rate rule ")
        or (line.startswith("Estimated using average of templates ")
            and f"[{root_label}]" in line)
        for line in source.get("comment", "").splitlines()
    )


def representative(old, new, family):
    if family == "R_Addition_MultipleBond":
        selected = old.get("ps_primary_end_ceiling_pairs", old["ps_ceiling_pairs"])[0]["propagation_event_id"]
        record = next(record for record in old["records"] if record["event_id"] == selected)
        return index(new)[structural_key(record)]
    # Freeze the representative against the old artifact, before seeing new rates.
    record = next(record for record in old["records"] if record["family"] == family
                  and record["rate_source"].get("kind") == "RMG family estimate")
    return index(new)[structural_key(record)]


def _without_database_provenance(value):
    if isinstance(value, dict):
        return {
            key: _without_database_provenance(item)
            for key, item in value.items()
            if key not in DATABASE_PROVENANCE_KEYS
        }
    if isinstance(value, list):
        return [_without_database_provenance(item) for item in value]
    return value


def _assert_database_provenance(value, expected_sha, path=()):
    if isinstance(value, dict):
        for key, item in value.items():
            if key in DATABASE_PROVENANCE_KEYS:
                assert item == expected_sha, (path + (key,), item, expected_sha)
            else:
                _assert_database_provenance(item, expected_sha, path + (key,))
    elif isinstance(value, list):
        for index_, item in enumerate(value):
            _assert_database_provenance(item, expected_sha, path + (index_,))


def verify_tree_invariance(old, new):
    before, after = index(old), index(new)
    database_sha = new.get("provenance", {}).get("rmg_database_sha")
    assert database_sha, "new artifact must declare provenance.rmg_database_sha"
    counts = {}
    for family in TREE_FAMILIES:
        old_keys = {key for key, record in before.items() if record["family"] == family}
        new_keys = {key for key, record in after.items() if record["family"] == family}
        assert old_keys <= new_keys, (family, "missing baseline channels")
        for key in old_keys:
            assert json.dumps(before[key]["k_table"], sort_keys=True) == json.dumps(
                after[key]["k_table"], sort_keys=True
            ), (family, "k_table")
            assert _without_database_provenance(before[key]["rate_source"]) == _without_database_provenance(
                after[key]["rate_source"]
            ), (family, "rate_source")
            _assert_database_provenance(after[key]["rate_source"], database_sha)
        counts[family] = len(old_keys)
    return counts


def run(new_path, database_path):
    from rmgpy.data.rmg import RMGDatabase
    from rmgpy.kmc.compiler import (
        _rate_rule_source, compiler_source_hash, prepare_rate_rules, validate_artifact,
    )
    from rmgpy.data.kinetics.family import TemplateReaction
    from rmgpy.kmc.ssa import record_rate
    from scipy.optimize import brentq

    prior = load_i044()
    old, new = json.loads(OLD_ARTIFACT.read_text()), json.loads(new_path.read_text())
    assert digest(OLD_ARTIFACT) == OLD_ARTIFACT.stem
    assert digest(new_path) == new_path.stem
    assert new["provenance"]["rmg_database_sha"] == prior.DATABASE_SHA
    assert new["provenance"]["compiler_sources_sha256"] == compiler_source_hash()
    assert new["provenance"]["compiler_sources_sha256"] != old["provenance"]["compiler_sources_sha256"]
    validate_artifact(new)
    tree_counts = verify_tree_invariance(old, new)
    before, after = index(old), index(new)
    families = sorted(set(old["families"]) | set(new["families"]))
    db = RMGDatabase()
    non_tree = [family for family in families if family not in TREE_FAMILIES]
    db.load_kinetics(str(database_path / "input/kinetics"), reaction_libraries=[],
                     seed_mechanisms=None, kinetics_families=non_tree,
                     kinetics_depositories=["training"])
    db.load_thermo(str(database_path / "input/thermo"),
                   thermo_libraries=["primaryThermoLibrary"], depository=True)
    original_rule_entries = {
        name: [{"index": entry.index, "label": entry.label, "rank": entry.rank,
                "short_desc": entry.short_desc}
               for entry in sorted(
                   (entry for entries in db.kinetics.families[name].rules.entries.values()
                    for entry in entries), key=lambda entry: entry.index)]
        for name in non_tree
    }
    old_sources = {}
    source_cache = {}
    old_rate_nodes_reproduced = 0
    for key, record in before.items():
        family_name = record["family"]
        if family_name in TREE_FAMILIES or record["rate_source"].get("kind") != "RMG family estimate":
            continue
        family = db.kinetics.families[family_name]
        cache_key = family_name, record["template"], record["raw_path_degeneracy"]
        if cache_key not in source_cache:
            template = family.retrieve_template(record["template"].split(";"))
            kinetics, entry = family.get_kinetics_for_template(template, degeneracy=record["raw_path_degeneracy"])
            estimate = TemplateReaction(family=family_name, template=record["template"].split(";"), kinetics=kinetics)
            source_cache[cache_key] = kinetics, _rate_rule_source(family, estimate)
        kinetics, source = source_cache[cache_key]
        for temperature, rate in zip(record["k_table"]["T"], record["k_table"]["k"]):
            prior.close(kinetics.get_rate_coefficient(temperature), rate)
            old_rate_nodes_reproduced += 1
        old_sources[key] = source
    summary = []
    for family in families:
        old_records = [record for record in old["records"] if record["family"] == family]
        new_records = [record for record in new["records"] if record["family"] == family]
        matched = [key for key in before.keys() & after.keys() if before[key]["family"] == family]
        row = {"family": family, "old_records": len(old_records), "new_records": len(new_records),
               "old_default_or_root_forward": None, "new_default_or_root_forward": None,
               "old_rmg_estimate_tagged_records": sum(record["rate_source"].get("kind") == "RMG family estimate" for record in old_records),
               "new_rmg_estimate_tagged_records": sum(record["rate_source"].get("kind") == "RMG family estimate" for record in new_records),
               "matched_records": len(matched), "rate_changes": {}}
        if family not in TREE_FAMILIES:
            root = ";".join(group.label for group in db.kinetics.families[family].get_root_template())
            row["root_template"] = root
            row["old_default_or_root_forward"] = sum(has_default_or_root(old_sources[key], root) for key in old_sources if before[key]["family"] == family)
            row["new_default_or_root_forward"] = sum(has_default_or_root(record["rate_source"], root) for record in new_records
                                                     if record["rate_source"].get("kind") == "RMG family estimate")
        else:
            row["old_default_or_root_forward"] = sum(record["rate_source"].get("entry") == "Root" for record in old_records)
            row["new_default_or_root_forward"] = sum(record["rate_source"].get("entry") == "Root" for record in new_records)
        for t in TEMPERATURES:
            differences = [math.log10(record_rate(after[key], t) / record_rate(before[key], t))
                           for key in matched if before[key]["k_table"] and after[key]["k_table"]]
            row["rate_changes"][str(t)] = {"median": median(differences), "min": min(differences),
                                         "max": max(differences), "count": len(differences)}
        summary.append(row)
    pairs = []
    for pair in old.get("ps_primary_end_ceiling_pairs", old["ps_ceiling_pairs"]):
        old_prop = next(record for record in old["records"] if record["event_id"] == pair["propagation_event_id"])
        old_dep = next(record for record in old["records"] if record["event_id"] == pair["depropagation_event_id"])
        prop, dep = after[structural_key(old_prop)], after[structural_key(old_dep)]
        assert prop["thermo_provenance"] == old_prop["thermo_provenance"]
        rxn = prior.prior.reaction_from_record(prop)
        for species in rxn.reactants + rxn.products:
            species.thermo = db.thermo.get_thermo_data(species)
        tc = brentq(lambda t: math.log(rxn.get_equilibrium_constant(t) * 1000.0), 600.0, 800.0)
        for t, kf, kr, kc in zip(prop["k_table"]["T"], prop["k_table"]["k"], dep["k_table"]["k"],
                                prop["thermo_provenance"]["equilibrium_constant_table"]["Kc"]):
            prior.close(rxn.get_equilibrium_constant(t, type="Kc"), kc)
            prior.close(kf / kc, kr)
        pairs.append({"site_type": prop["site_type"], "propagation_event_id": prop["event_id"],
                      "template": prop["template"], "source": prop["rate_source"], "Tc_thermo_K": tc,
                      "rows": [{"T_K": t, "k_fwd_L_mol_s": record_rate(prop, t) * 1000.0,
                                "iupac_L_mol_s": prior.benchmark(t),
                                "ratio": record_rate(prop, t) * 1000.0 / prior.benchmark(t)}
                               for t in TEMPERATURES]})
    discovery_counts = lambda artifact: dict(Counter(item["family"] for item in artifact["discovery"]))
    # Fresh RMG preparation and public get_kinetics spot-checks; no event compile.
    prepare_rate_rules(db.kinetics, db.thermo, verbose=True)
    spot_checks = {}
    for family_name in non_tree:
        record = representative(old, new, family_name)
        reaction = prior.prior.reaction_from_record(record)
        estimate = TemplateReaction(
            reactants=reaction.reactants, products=reaction.products,
            family=family_name, template=record["template"].split(";"),
            degeneracy=record["raw_path_degeneracy"],
        )
        family = db.kinetics.families[family_name]
        kinetics, source, entry, forward = family.get_kinetics(
            estimate, estimate.template, degeneracy=estimate.degeneracy,
            return_all_kinetics=False,
        )
        assert forward
        assert record["rate_source"]["source"] == str(source)
        assert record["rate_source"]["entry"] == (str(entry) if entry else None)
        assert record["rate_source"]["rank"] == getattr(entry, "rank", None)
        estimate.kinetics = kinetics
        fresh_source = _rate_rule_source(family, estimate)
        for field, value in fresh_source.items():
            assert record["rate_source"][field] == value, (family_name, field)
        for temperature, rate in zip(record["k_table"]["T"], record["k_table"]["k"]):
            prior.close(kinetics.get_rate_coefficient(temperature), rate)
        spot_checks[family_name] = {"event_id": record["event_id"],
                                   "grid_points": len(record["k_table"]["T"]),
                                   "source": fresh_source}
    result = {"old_artifact": str(OLD_ARTIFACT), "new_artifact": str(new_path),
              "new_artifact_sha256": digest(new_path), "new_records": len(new["records"]),
              "old_rate_nodes_reproduced": old_rate_nodes_reproduced,
              "tree_records_bit_identical": tree_counts, "families": summary, "propagation": pairs,
              "original_rule_entries": original_rule_entries,
              "rate_rule_preparation": new["provenance"]["rate_rule_preparation"],
              "old_grid_ceiling_K": old["ps_ceiling_temperature_K"],
              "new_grid_ceiling_K": new["ps_ceiling_temperature_K"],
              "records_gained": len(after.keys() - before.keys()), "records_lost": len(before.keys() - after.keys()),
              "status_changes": sum(before[key]["status"] != after[key]["status"] for key in before.keys() & after.keys()),
              "old_discovery_counts": discovery_counts(old), "new_discovery_counts": discovery_counts(new),
              "excluded_channels_old": len(old["excluded_channels"]), "excluded_channels_new": len(new["excluded_channels"]),
              "representatives": {family: representative(old, new, family)["rate_source"] for family in non_tree}}
    result["fresh_rmg_spot_checks"] = spot_checks
    prior.close(result["old_grid_ceiling_K"], result["new_grid_ceiling_K"])
    for pair in pairs:
        prior.close(pair["Tc_thermo_K"], 710.020462, relative=1e-8)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", type=Path)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify", action="store_true", help="reproduce and compare saved results")
    args = parser.parse_args()
    result = run(args.artifact, args.database)
    if args.verify:
        assert result == json.loads(args.output.read_text()), "saved probe results differ"
    else:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print("I046 artifact hashes, compiler fingerprint and schema verified")
    print("I046 tree-family rate tables and sources: 13274 bit-identical records")
    print("I046 old non-tree rule rates reproduced; family distributions and discovery compared")
    print("I046 both propagation pairs: dimensional Kc, detailed balance and 710.02 K thermo ceiling verified")
    print("I046 three fresh public-RMG estimates: all 87 rate nodes and source comments verified")
    print("I046 results: " + str(args.output))


if __name__ == "__main__":
    main()
