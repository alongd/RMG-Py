"""Reproduce the I-032 event-set initiation diagnosis.

Run from the repository root with the project RMG environment and PYTHONPATH=$PWD.
The output is canonical JSON so the companion verifier can check every reported number.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from collections import Counter
from dataclasses import replace
from pathlib import Path
from typing import Any, Iterable

from rmgpy.data.rmg import RMGDatabase
from rmgpy.kmc.compiler import (
    EventSetCompiler,
    PS_FAMILY_CANDIDATES,
    PS_PROXY_UNITS,
    _orient_to_proxy,
    _reverse_view,
    ps_proxy_set,
)
from rmgpy.kmc.met import _is_forward_met_record
from rmgpy.molecule.molecule import Molecule
from rmgpy.species import Species


TEMPERATURES = (600.0, 700.0, 800.0)
REPOSITORY_BASE_SHA = "0a56f94a68a463adf438a365d9f932df83acf326"
DATABASE_SHA = "4a12d36fcdc193ede82c8d1ab5c1653495d445bc"
DESIGN_RULING = Path(
    "/home/alon/Code/polymer-pm/reports/" "R-001_v3_2026-09-24_kmc-design-ruling.md"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def git_sha(path: Path) -> str:
    return subprocess.check_output(
        ["git", "-C", str(path), "rev-parse", "HEAD"], text=True
    ).strip()


def require_git_ancestor(path: Path, ancestor: str) -> None:
    subprocess.check_call(
        ["git", "-C", str(path), "merge-base", "--is-ancestor", ancestor, "HEAD"]
    )


def progress(message: str) -> None:
    print(f"[I-032] {message}", file=sys.stderr, flush=True)


def source_line(path: Path, needle: str) -> int:
    for number, line in enumerate(path.read_text().splitlines(), 1):
        if needle in line:
            return number
    raise AssertionError(f"{needle!r} not found in {path}")


def counts(values: Iterable[Any]) -> dict[str, int]:
    return dict(
        sorted(
            Counter("null" if value is None else str(value) for value in values).items()
        )
    )


def smiles(items: Iterable[Any]) -> list[str]:
    return [item.molecule[0].to_smiles() for item in items]


def generate(database: RMGDatabase, reactants: Iterable[Any], family: str) -> list[Any]:
    return database.kinetics.generate_reactions_from_families(
        [item.copy(deep=True) for item in reactants],
        only_families=[family],
        resonance=True,
    )


def radicals_in_adjacencies(adjacencies: Iterable[str]) -> int:
    return sum(
        atom.radical_electrons
        for adjacency in adjacencies
        for atom in Molecule().from_adjacency_list(adjacency).atoms
    )


def reaction_signature(reaction: Any) -> dict[str, Any]:
    return {
        "reactants": smiles(reaction.reactants),
        "products": smiles(reaction.products),
        "degeneracy": float(reaction.degeneracy),
        "is_forward": bool(reaction.is_forward),
        "template": [getattr(item, "label", str(item)) for item in reaction.template],
    }


def reaction_template(reaction: Any) -> str:
    return ";".join(
        getattr(item, "label", str(item))
        for item in (getattr(reaction, "template", None) or [])
    )


def choose_backbone_recombination(reactions: Iterable[Any]) -> Any:
    target = {
        "[CH](CCc1ccccc1)c1ccccc1",
        "[CH2]C(C)c1ccccc1",
    }
    return next(
        reaction for reaction in reactions if set(smiles(reaction.reactants)) == target
    )


def choose_h_abstraction(reactions: Iterable[Any]) -> Any:
    return next(
        reaction
        for reaction in reactions
        if "C=Cc1[c]cccc1" in smiles(reaction.products)
    )


def thermo_coverage(
    database: RMGDatabase, adjacencies: Iterable[str]
) -> dict[str, Any]:
    failures = []
    unique = sorted(set(adjacencies))
    for adjacency in unique:
        species = Species(molecule=[Molecule().from_adjacency_list(adjacency)])
        try:
            database.thermo.get_thermo_data(species)
        except Exception as error:  # report the actual unsupported structures
            failures.append(
                {"smiles": species.molecule[0].to_smiles(), "error": str(error)}
            )
    return {"unique_species": len(unique), "failures": failures}


def solute_coverage(
    database: RMGDatabase, adjacencies: Iterable[str]
) -> dict[str, Any]:
    failures = []
    unique = sorted(set(adjacencies))
    for adjacency in unique:
        species = Species(molecule=[Molecule().from_adjacency_list(adjacency)])
        try:
            data = database.solvation.get_solute_data(species)
            if data is None:
                raise ValueError("no solute descriptors returned")
        except Exception as error:  # report the actual unsupported structures
            failures.append(
                {"smiles": species.molecule[0].to_smiles(), "error": str(error)}
            )
    return {"unique_species": len(unique), "failures": failures}


def main() -> None:
    args = parse_args()
    repository = Path(__file__).resolve().parents[5]
    compiler_path = repository / "rmgpy/kmc/compiler.py"
    met_path = repository / "rmgpy/kmc/met.py"
    solvation_path = repository / "rmgpy/data/solvation.py"
    real_test_path = repository / "test/rmgpy/kmc/compilerRealTest.py"

    require_git_ancestor(repository, REPOSITORY_BASE_SHA)
    if git_sha(args.database) != DATABASE_SHA:
        raise RuntimeError(f"database must be pinned at {DATABASE_SHA}")

    progress("loading pinned RMG database")
    database = RMGDatabase()
    database.load_kinetics(
        str(args.database / "input/kinetics"),
        reaction_libraries=[],
        seed_mechanisms=None,
        kinetics_families=list(PS_FAMILY_CANDIDATES),
        kinetics_depositories=["training"],
    )
    database.load_thermo(
        str(args.database / "input/thermo"),
        thermo_libraries=["primaryThermoLibrary"],
        depository=True,
    )
    database.load_solvation(str(args.database / "input/solvation"))
    progress("database loaded")

    family_root = args.database / "input/kinetics/families"
    family_universe = sorted(
        path.name for path in family_root.iterdir() if (path / "groups.py").is_file()
    )
    proxies = ps_proxy_set(PS_PROXY_UNITS)
    proxy_by_site = {proxy.site_type: proxy for proxy in proxies}
    progress("discovering scheduled L=3 family reactions")
    active, excluded, reaction_cache = EventSetCompiler.discover_family_reactions(
        database.kinetics,
        proxies,
        PS_FAMILY_CANDIDATES,
        family_universe=family_universe,
    )
    compiler = EventSetCompiler(
        database.kinetics,
        proxies,
        active,
        excluded_families=excluded,
        database_path=args.database,
        thermo_database=database.thermo,
        reaction_cache=reaction_cache,
    )
    progress("compiling the production artifact")
    artifact = compiler.compile()
    progress("production artifact compiled")
    records = artifact["records"]

    span_radius = artifact["inputs"]["span_radius"]
    faithful_units = 2 * span_radius + 3
    faithful_proxies = list(ps_proxy_set(faithful_units))
    faithful_by_site = {proxy.site_type: proxy for proxy in faithful_proxies}
    end_styrene_index = next(
        index
        for index, proxy in enumerate(faithful_proxies)
        if proxy.site_type == "end_radical+styrene"
    )
    end_styrene = faithful_proxies[end_styrene_index]
    faithful_proxies[end_styrene_index] = replace(
        end_styrene,
        reactants=(
            faithful_by_site["end_radical"].reactants[0],
            end_styrene.reactants[1],
        ),
    )
    compiled_keys = {
        (record["site_type"], record["family"], record["template"])
        for record in records
        if record["site_type"] != "J_ring"
    }
    faithful_generated = {}
    faithful_missing_keys = set()
    for proxy in faithful_proxies:
        generated = []
        for family in PS_FAMILY_CANDIDATES:
            progress(f"generating faithful C3 {proxy.site_type} / {family}")
            generated.extend(generate(database, proxy.reactants, family))
        faithful_generated[proxy.site_type] = generated
        faithful_missing_keys.update(
            (proxy.site_type, reaction.family, reaction_template(reaction))
            for reaction in generated
            if (proxy.site_type, reaction.family, reaction_template(reaction))
            not in compiled_keys
        )

    progress("probing L=3 H-abstraction and pristine initiation")
    forced_h_abstraction = {
        proxy.site_type: generate(database, proxy.reactants, "H_Abstraction")
        for proxy in proxies
    }
    pristine = proxy_by_site["pristine"]
    augmented_h_abstraction = {}
    for radical_site in ("interior_radical", "end_radical", "junction_radical"):
        pair = (
            proxy_by_site[radical_site].reactants[0],
            pristine.reactants[0],
        )
        augmented_h_abstraction[radical_site + "+pristine"] = generate(
            database, pair, "H_Abstraction"
        )

    pristine_all_families = {
        family: generate(database, pristine.reactants, family)
        for family in PS_FAMILY_CANDIDATES
    }
    pristine_recombinations = pristine_all_families["R_Recombination"]
    backbone_source = choose_backbone_recombination(pristine_recombinations)
    backbone_homolysis, backbone_orientation = _orient_to_proxy(
        pristine, backbone_source
    )

    rate_compiler = EventSetCompiler(
        database.kinetics,
        (),
        PS_FAMILY_CANDIDATES,
        thermo_database=database.thermo,
        temperature_grid=TEMPERATURES,
    )
    homolysis_table, homolysis_source = rate_compiler._rate_table(backbone_homolysis)

    h_abstraction_reaction = choose_h_abstraction(
        forced_h_abstraction["end_radical+styrene"]
    )
    h_forward_table, h_forward_source = rate_compiler._rate_table(
        h_abstraction_reaction
    )
    h_reverse_table, h_reverse_source = rate_compiler._rate_table(
        _reverse_view(h_abstraction_reaction)
    )

    recombination_proxy = replace(
        proxy_by_site["end_radical+end_radical"], frontier=True
    )
    recombination_reaction = generate(
        database, recombination_proxy.reactants, "R_Recombination"
    )[0]
    frontier_record = rate_compiler._record(
        recombination_proxy, recombination_reaction, {}
    )

    all_adjacencies = [
        adjacency
        for record in records
        for adjacency in record["reactant_graphs"] + record["product_graphs"]
    ]
    for reaction in [
        backbone_source,
        h_abstraction_reaction,
        *[reaction for group in augmented_h_abstraction.values() for reaction in group],
    ]:
        all_adjacencies.extend(
            item.molecule[0].to_adjacency_list(remove_h=False)
            for item in list(reaction.reactants) + list(reaction.products)
        )

    non_ring = [
        record for record in records if record["inventory_class"] != "R1:J_ring"
    ]
    no_ring_initiators = [
        record
        for record in non_ring
        if radicals_in_adjacencies(record["reactant_graphs"]) == 0
        and record["radical_delta"] > 0
        and record["status"] != "refused"
    ]
    ring_initiators = [
        record
        for record in records
        if record["inventory_class"] == "R1:J_ring"
        and radicals_in_adjacencies(record["reactant_graphs"]) == 0
        and record["radical_delta"] > 0
        and record["status"] != "refused"
    ]
    compiled_pair_count = len(non_ring) + (len(records) - len(non_ring)) // 2

    expected_homolysis = {
        str(int(temperature)): {
            "Ea_320_kJ_mol": 1.0e15
            * math.exp(-320_000.0 / (8.314462618 * temperature)),
            "Ea_290_kJ_mol": 1.0e15
            * math.exp(-290_000.0 / (8.314462618 * temperature)),
        }
        for temperature in TEMPERATURES
    }

    test_text = real_test_path.read_text()
    c9_start = test_text.index("def test_c7_c9_inventory_and_real_artifact")
    c9_text = test_text[c9_start:]
    next_test = c9_text.find("\ndef test_", 1)
    if next_test >= 0:
        c9_text = c9_text[:next_test]

    progress("estimating bounded-set thermo and solvation coverage")
    thermo_result = thermo_coverage(database, all_adjacencies)
    solute_result = solute_coverage(database, all_adjacencies)
    progress("writing canonical results")
    result = {
        "provenance": {
            "repository_base_sha": REPOSITORY_BASE_SHA,
            "database_sha": git_sha(args.database),
        },
        "source_lines": {
            "family_candidates": source_line(compiler_path, '"H_Abstraction",'),
            "proxy_set": source_line(compiler_path, "def ps_proxy_set"),
            "per_proxy_scheduling": source_line(
                compiler_path, "proxy_candidates = sorted("
            ),
            "generation": source_line(
                compiler_path, "for reaction in self._generate(proxy):"
            ),
            "forward_record": source_line(
                compiler_path, "forward = self._record(proxy, oriented"
            ),
            "irreversible_overwrite": source_line(
                compiler_path,
                'if not is_ring and getattr(reaction, "reversible", True):',
            ),
            "generic_ring_reverse": source_line(
                compiler_path, 'if forward.inventory_class == "R1:J_ring":'
            ),
            "reverse_kc": source_line(
                compiler_path, "source_reaction.get_equilibrium_constant"
            ),
            "met_refused_filter": source_line(
                met_path, 'if _field(record, "status", "enabled") == "refused":'
            ),
            "solvent_data": source_line(solvation_path, "class SolventData"),
            "solute_data": source_line(solvation_path, "class SoluteData"),
            "temperature_dependent_solvation": source_line(
                solvation_path, "def get_T_dep_solvation_energy_from_LSER_298"
            ),
            "c3": source_line(
                real_test_path,
                "def test_c3_bounded_l3_generation_has_no_missing_family_template",
            ),
            "c3_comparison": source_line(real_test_path, "assert required <= compiled"),
            "c9": source_line(
                real_test_path, "def test_c7_c9_inventory_and_real_artifact"
            ),
            "design_missing_thermo_rule": source_line(
                DESIGN_RULING, "Some reversible pairs have a participant"
            ),
            "design_c3_c9": source_line(DESIGN_RULING, "**C3** span bound"),
        },
        "artifact": {
            "families": artifact["families"],
            "excluded_h_abstraction_reason": artifact["excluded_families"][
                "H_Abstraction"
            ],
            "record_count": len(records),
            "records_by_family": counts(record["family"] for record in records),
            "records_by_inventory_class": counts(
                record["inventory_class"] for record in records
            ),
            "records_by_family_and_inventory_class": counts(
                record["family"]
                + ":"
                + (
                    "null"
                    if record["inventory_class"] is None
                    else record["inventory_class"]
                )
                for record in records
            ),
            "records_by_status": counts(record["status"] for record in records),
            "records_by_family_and_status": counts(
                record["family"] + ":" + record["status"] for record in records
            ),
            "irreversible_list_entries": len(artifact["irreversible_pairs"]),
            "non_ring_records": len(non_ring),
            "non_ring_irreversible_records": sum(
                record["status"] == "irreversible" for record in non_ring
            ),
            "linked_reverse_records": sum(
                bool(record["reverse_of"]) for record in records
            ),
            "non_ring_linked_reverse_records": sum(
                bool(record["reverse_of"]) for record in non_ring
            ),
            "ps_ceiling_temperature_K": artifact.get("ps_primary_end_ceiling_pairs", artifact["ps_ceiling_pairs"])[0]["temperature_K"],
            "ps_ceiling_pair_count": len(artifact.get("ps_primary_end_ceiling_pairs", artifact["ps_ceiling_pairs"])),
            "no_ring_zero_radical_initiators": len(no_ring_initiators),
            "ring_zero_radical_initiators": len(ring_initiators),
        },
        "proxies": {
            proxy.site_type: {
                "participants": list(proxy.participant_site_types),
                "scheduled_families": list(proxy.metadata["family_candidates"]),
                "frontier": proxy.frontier,
            }
            for proxy in proxies
        },
        "h_abstraction": {
            "forced_counts_on_declared_inputs": {
                site: len(reactions) for site, reactions in forced_h_abstraction.items()
            },
            "forced_total_on_declared_inputs": sum(
                len(reactions) for reactions in forced_h_abstraction.values()
            ),
            "augmented_radical_plus_pristine_counts": {
                site: len(reactions)
                for site, reactions in augmented_h_abstraction.items()
            },
            "rate_pair": reaction_signature(h_abstraction_reaction),
            "k_forward_m3_mol_s": h_forward_table["k"],
            "Kc_forward": h_reverse_source["equilibrium_constant_table"]["Kc"],
            "k_reverse_m3_mol_s": h_reverse_table["k"],
            "forward_rate_reference": h_forward_source["reference"],
            "reverse_rate_reference": h_reverse_source["reference"],
        },
        "homolysis": {
            "forced_pristine_counts_all_candidate_families": {
                family: len(reactions)
                for family, reactions in pristine_all_families.items()
            },
            "forced_pristine_dissociations": len(pristine_recombinations),
            "orientation": backbone_orientation,
            "association_pair": reaction_signature(backbone_source),
            "homolysis_reactants": smiles(backbone_homolysis.reactants),
            "homolysis_products": smiles(backbone_homolysis.products),
            "Kc_association": homolysis_source["equilibrium_constant_table"]["Kc"],
            "k_homolysis_s-1": homolysis_table["k"],
            "reverse_rate_reference": homolysis_source["reference"],
            "expected_A_1e15_s-1_bracket": expected_homolysis,
        },
        "thermo_coverage": thermo_result,
        "solvation_solute_coverage": solute_result,
        "solvation_support": {
            "models": [
                "Abraham gas-to-solvent Gibbs correction at 298 K",
                "Mintz solvation enthalpy correction at 298 K",
                "CoolProp-based temperature-dependent K-factor below solvent critical temperature",
            ],
            "solute_descriptors": ["S", "B", "E", "L", "A", "V"],
            "stand_in_solvent_descriptors": [
                "s_g",
                "b_g",
                "e_g",
                "l_g",
                "a_g",
                "c_g",
                "s_h",
                "b_h",
                "e_h",
                "l_h",
                "a_h",
                "c_h",
                "name_in_coolprop",
            ],
            "h_abstraction_optional_solvent_descriptors": ["alpha", "beta"],
            "h_abstraction_intrinsic_correction_used_by_rmg": False,
        },
        "frontier_overwrite": {
            "all_production_proxies_frontier_false": all(
                not proxy.frontier for proxy in proxies
            ),
            "synthetic_frontier_status": frontier_record.status,
            "synthetic_frontier_reason": frontier_record.status_reason,
            "admitted_by_met_filter": _is_forward_met_record(frontier_record, "R1"),
        },
        "acceptance_tests": {
            "c3_required_source": "oracle keys generated only from each proxy's scheduled_families",
            "c3_compiled_source": "artifact keys excluding site_type J_ring",
            "c3_compares_only_site_family_template": True,
            "c9_checks_irreversible_list_nonempty": 'assert artifact["irreversible_pairs"]'
            in c9_text,
            "c9_recomputes_compiled_ceiling": 'assert recomputed == pair["temperature_K"]'
            in c9_text,
            "c9_has_literature_comparison": "literature" in c9_text.lower(),
            "faithful_c3_span_radius": span_radius,
            "faithful_c3_polymer_units": faithful_units,
            "faithful_c3_generated_reactions_by_site": {
                site: len(reactions) for site, reactions in faithful_generated.items()
            },
            "faithful_c3_generated_reactions_by_site_and_family": counts(
                site + ":" + reaction.family
                for site, reactions in faithful_generated.items()
                for reaction in reactions
            ),
            "faithful_c3_missing_unique_key_count": len(faithful_missing_keys),
            "faithful_c3_missing_unique_keys_by_site_and_family": counts(
                site + ":" + family for site, family, _template in faithful_missing_keys
            ),
        },
        "option_irreversible_share": {
            "a_gas_phase_Kc": {
                "irreversible_pairs": 0,
                "compiled_pairs": compiled_pair_count,
            },
            "b_solvated_Kc_if_reference_available": {
                "irreversible_pairs": 0,
                "compiled_pairs": compiled_pair_count,
            },
            "c_status_quo": {
                "irreversible_pairs": sum(
                    record["status"] == "irreversible" for record in non_ring
                ),
                "compiled_pairs": compiled_pair_count,
            },
        },
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
