"""Check the canonical JSON emitted by run_probe.py against the report claims."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def close_list(actual: list[float], expected: list[float]) -> bool:
    return len(actual) == len(expected) and all(
        math.isclose(left, right, rel_tol=1.0e-12, abs_tol=0.0)
        for left, right in zip(actual, expected)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    args = parser.parse_args()
    result = json.loads(args.results.read_text())
    report = Path(__file__).resolve().parent.parent / "I032_initiation_probe.md"
    report_text = report.read_text()

    assert result["provenance"] == {
        "database_sha": "4a12d36fcdc193ede82c8d1ab5c1653495d445bc",
        "repository_base_sha": "0a56f94a68a463adf438a365d9f932df83acf326",
    }
    artifact = result["artifact"]
    assert artifact["families"] == [
        "Disproportionation",
        "R_Addition_MultipleBond",
        "R_Recombination",
        "intra_H_migration",
    ]
    assert artifact["record_count"] == 314
    assert artifact["records_by_family"] == {
        "Disproportionation": 181,
        "R_Addition_MultipleBond": 22,
        "R_Recombination": 7,
        "intra_H_migration": 104,
    }
    assert artifact["records_by_inventory_class"] == {
        "R1:J_ring": 6,
        "R1:quinoid_disproportionation": 180,
        "null": 128,
    }
    assert artifact["records_by_status"] == {"enabled": 6, "irreversible": 308}
    assert artifact["non_ring_records"] == 308
    assert artifact["non_ring_irreversible_records"] == 308
    assert artifact["non_ring_linked_reverse_records"] == 0
    assert artifact["linked_reverse_records"] == 6
    assert artifact["no_ring_zero_radical_initiators"] == 0
    assert artifact["ring_zero_radical_initiators"] == 3
    assert artifact["ps_ceiling_pair_count"] == 1
    assert math.isclose(
        artifact["ps_ceiling_temperature_K"],
        710.2486726349888,
        rel_tol=1.0e-12,
    )

    assert all(
        "H_Abstraction" not in proxy["scheduled_families"]
        for proxy in result["proxies"].values()
    )
    assert result["proxies"]["pristine"]["scheduled_families"] == []
    h_abstraction = result["h_abstraction"]
    assert h_abstraction["forced_counts_on_declared_inputs"] == {
        "doubly_featured": 0,
        "end_radical": 0,
        "end_radical+end_radical": 14,
        "end_radical+styrene": 5,
        "interior_radical": 0,
        "junction_radical": 0,
        "junction_radical+end_radical": 28,
        "pristine": 0,
    }
    assert h_abstraction["forced_total_on_declared_inputs"] == 47
    assert h_abstraction["augmented_radical_plus_pristine_counts"] == {
        "end_radical+pristine": 14,
        "interior_radical+pristine": 44,
        "junction_radical+pristine": 14,
    }
    assert close_list(
        h_abstraction["Kc_forward"],
        [3.0277701967396972e-05, 1.25829988581265e-04, 3.650914241632471e-04],
    )
    assert close_list(
        h_abstraction["k_reverse_m3_mol_s"],
        [1.504827616776421, 1.1999724279753776, 1.0158114743285866],
    )

    homolysis = result["homolysis"]
    assert homolysis["forced_pristine_counts_all_candidate_families"] == {
        "Disproportionation": 0,
        "H_Abstraction": 0,
        "R_Addition_MultipleBond": 0,
        "R_Recombination": 23,
        "intra_H_migration": 0,
    }
    assert homolysis["forced_pristine_dissociations"] == 23
    assert homolysis["orientation"] == "reversed"
    assert close_list(
        homolysis["Kc_association"],
        [1.1689496351602934e18, 1.63104476045059e14, 2.1751754270151132e11],
    )
    assert close_list(
        homolysis["k_homolysis_s-1"],
        [1.8759785971763749e-13, 2.024922054745522e-09, 2.043769064160407e-06],
    )
    for temperature, rate in zip((600, 700, 800), homolysis["k_homolysis_s-1"]):
        bracket = homolysis["expected_A_1e15_s-1_bracket"][str(temperature)]
        assert bracket["Ea_320_kJ_mol"] <= rate <= bracket["Ea_290_kJ_mol"]

    assert result["thermo_coverage"] == {"failures": [], "unique_species": 195}
    assert result["solvation_solute_coverage"] == {
        "failures": [],
        "unique_species": 195,
    }
    assert result["frontier_overwrite"] == {
        "admitted_by_met_filter": True,
        "all_production_proxies_frontier_false": True,
        "synthetic_frontier_reason": "condensed-phase reference thermochemistry unavailable",
        "synthetic_frontier_status": "irreversible",
    }
    assert result["acceptance_tests"] == {
        "c3_compares_only_site_family_template": True,
        "c3_compiled_source": "artifact keys excluding site_type J_ring",
        "c3_required_source": "oracle keys generated only from each proxy's scheduled_families",
        "c9_checks_irreversible_list_nonempty": True,
        "c9_has_literature_comparison": False,
        "c9_recomputes_compiled_ceiling": True,
        "faithful_c3_generated_reactions_by_site": {
            "doubly_featured": 115,
            "end_radical": 54,
            "end_radical+end_radical": 546,
            "end_radical+styrene": 261,
            "interior_radical": 89,
            "junction_radical": 46,
            "junction_radical+end_radical": 1045,
            "pristine": 39,
        },
        "faithful_c3_generated_reactions_by_site_and_family": {
            "doubly_featured:R_Addition_MultipleBond": 8,
            "doubly_featured:R_Recombination": 35,
            "doubly_featured:intra_H_migration": 72,
            "end_radical+end_radical:Disproportionation": 501,
            "end_radical+end_radical:H_Abstraction": 24,
            "end_radical+end_radical:R_Addition_MultipleBond": 20,
            "end_radical+end_radical:R_Recombination": 1,
            "end_radical+styrene:Disproportionation": 250,
            "end_radical+styrene:H_Abstraction": 5,
            "end_radical+styrene:R_Addition_MultipleBond": 6,
            "end_radical:R_Addition_MultipleBond": 3,
            "end_radical:R_Recombination": 37,
            "end_radical:intra_H_migration": 14,
            "interior_radical:R_Addition_MultipleBond": 5,
            "interior_radical:R_Recombination": 37,
            "interior_radical:intra_H_migration": 47,
            "junction_radical+end_radical:Disproportionation": 957,
            "junction_radical+end_radical:H_Abstraction": 48,
            "junction_radical+end_radical:R_Addition_MultipleBond": 39,
            "junction_radical+end_radical:R_Recombination": 1,
            "junction_radical:R_Addition_MultipleBond": 1,
            "junction_radical:R_Recombination": 38,
            "junction_radical:intra_H_migration": 7,
            "pristine:R_Recombination": 39,
        },
        "faithful_c3_missing_unique_key_count": 128,
        "faithful_c3_missing_unique_keys_by_site_and_family": {
            "doubly_featured:R_Addition_MultipleBond": 1,
            "doubly_featured:R_Recombination": 6,
            "doubly_featured:intra_H_migration": 24,
            "end_radical+end_radical:H_Abstraction": 5,
            "end_radical+end_radical:R_Addition_MultipleBond": 4,
            "end_radical+styrene:Disproportionation": 9,
            "end_radical+styrene:H_Abstraction": 3,
            "end_radical:R_Recombination": 7,
            "end_radical:intra_H_migration": 4,
            "interior_radical:R_Recombination": 7,
            "interior_radical:intra_H_migration": 13,
            "junction_radical+end_radical:Disproportionation": 5,
            "junction_radical+end_radical:H_Abstraction": 11,
            "junction_radical+end_radical:R_Addition_MultipleBond": 8,
            "junction_radical:R_Addition_MultipleBond": 1,
            "junction_radical:R_Recombination": 4,
            "junction_radical:intra_H_migration": 7,
            "pristine:R_Recombination": 9,
        },
        "faithful_c3_polymer_units": 5,
        "faithful_c3_span_radius": 1,
    }
    assert (
        sum(
            result["acceptance_tests"][
                "faithful_c3_generated_reactions_by_site"
            ].values()
        )
        == 2195
    )
    assert result["option_irreversible_share"] == {
        "a_gas_phase_Kc": {
            "compiled_pairs": 311,
            "irreversible_pairs": 0,
        },
        "b_solvated_Kc_if_reference_available": {
            "compiled_pairs": 311,
            "irreversible_pairs": 0,
        },
        "c_status_quo": {
            "compiled_pairs": 311,
            "irreversible_pairs": 308,
        },
    }
    for claim in (
        "| **Total** |  | **6** | **308** | **0** | **314** |",
        "produce 47 H-abstractions",
        "produces 44, 14, and 14 reactions",
        "produces 23 `R_Recombination`",
        "**710.2487 K**",
        "**0/311**",
        "**308/311 (99.04%)**",
        "128 missing unique keys",
        "2,195 reactions",
        "1.87598×10⁻¹³",
        "2.02492×10⁻⁹",
        "2.04377×10⁻⁶",
        "1.50483",
        "1.19997",
        "1.01581",
    ):
        assert claim in report_text, f"report is missing reproduced claim: {claim}"
    print("I-032 report claims reproduced from run_probe.py output")


if __name__ == "__main__":
    main()
