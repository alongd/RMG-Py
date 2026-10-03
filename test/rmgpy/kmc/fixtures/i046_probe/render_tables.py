"""Render or verify every numerical block in the I046 report."""

import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPORT = HERE.parent / "I046_rules_from_training.md"
spec = importlib.util.spec_from_file_location("i044_tables", HERE.parent / "i044_probe/render_tables.py")
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)


def blocks(result):
    number, table = prior.number, prior.table
    output = {}
    output["families"] = table(
        ["Family", "Old records", "New records", "Old Default/root estimates", "New Default/root estimates"],
        [[row["family"], row["old_records"], row["new_records"],
          row["old_default_or_root_forward"], row["new_default_or_root_forward"]]
         for row in result["families"]],
    )
    output["rate_changes"] = table(
        ["Family", "T (K)", "Matched records", "Median log10(new/old)", "Minimum", "Maximum"],
        [[row["family"], t, values["count"], number(values["median"]),
          number(values["min"]), number(values["max"])]
         for row in result["families"] for t, values in row["rate_changes"].items()],
    )
    output["preparation"] = table(
        ["Family", "Rules loaded", "Original Default entries", "After training", "After averaging"],
        [[family, values["rules_before"],
          sum(entry["short_desc"].strip().lower() == "default"
              for entry in result["original_rule_entries"][family]),
          values["rules_after_training"], values["rules_after_averaging"]]
         for family, values in result["rate_rule_preparation"]["families"].items()],
    )
    output["propagation"] = table(
        ["Proxy", "T (K)", "k_fwd (L/mol/s)", "IUPAC extrapolation (L/mol/s)", "Ratio"],
        [[pair["site_type"], row["T_K"], number(row["k_fwd_L_mol_s"]),
          number(row["iupac_L_mol_s"]), number(row["ratio"])]
         for pair in result["propagation"] for row in pair["rows"]],
    )
    output["provenance"] = table(
        ["Quantity", "Reproduced value"],
        [["New artifact SHA-256", result["new_artifact_sha256"]],
         ["New record count", result["new_records"]],
         ["Tree records with identical rates/sources", sum(result["tree_records_bit_identical"].values())],
         ["Old non-tree rate nodes reproduced", result["old_rate_nodes_reproduced"]],
         ["Old artifact-grid ceiling (K)", number(result["old_grid_ceiling_K"])],
         ["New artifact-grid ceiling (K)", number(result["new_grid_ceiling_K"])],
         *[[pair["site_type"] + " continuous thermo ceiling (K)", f'{pair["Tc_thermo_K"]:.6f}']
           for pair in result["propagation"]],
         ["Chemical records gained", result["records_gained"]],
         ["Chemical records lost", result["records_lost"]],
         ["Status changes", result["status_changes"]],
         ["Old excluded channels", result["excluded_channels_old"]],
         ["New excluded channels", result["excluded_channels_new"]]],
    )
    output["discovery"] = table(
        ["Family", "Old discovery records", "New discovery records"],
        [[row["family"], result["old_discovery_counts"].get(row["family"], 0),
          result["new_discovery_counts"].get(row["family"], 0)] for row in result["families"]],
    )
    output["sources"] = "\n\n".join(
        "Proxy `" + pair["site_type"] + "`, template `" + pair["template"] + "`:\n\n```text\n"
        + pair["source"]["comment"] + "\n```\n\n```json\n"
        + json.dumps({key: pair["source"].get(key) for key in ("entry", "rank", "exact", "rules", "training")},
                     indent=2, sort_keys=True) + "\n```"
        for pair in result["propagation"]
    )
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("--update", action="store_true")
    args = parser.parse_args()
    content = REPORT.read_text()
    for name, value in blocks(json.loads(args.results.read_text())).items():
        begin, end = f"<!-- BEGIN I046:{name} -->", f"<!-- END I046:{name} -->"
        left, rest = content.split(begin, 1)
        stored, right = rest.split(end, 1)
        expected = "\n" + value + "\n"
        if args.update:
            content = left + begin + expected + end + right
        else:
            assert stored == expected, name
    if args.update:
        REPORT.write_text(content)
    print("I046 all seven report blocks verified" if not args.update else "I046 report blocks updated")


if __name__ == "__main__":
    main()
