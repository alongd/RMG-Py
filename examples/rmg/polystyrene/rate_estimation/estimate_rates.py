"""Fit the polystyrene MOM channels from a compact kMC event-set artifact.

Records are selected by family, template, and canonical reactant/product
structures. Event IDs are output provenance only and never selectors.
"""

import argparse
import csv
import hashlib
import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np

from rmgpy import settings
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule


PS_DENSITY_KG_M3 = 1050.0
PS_REPEAT_SMILES = "[CH2][CH]c1ccccc1"
QSSA_R_GAS = 8.314

END_RADICAL = "CC(CC(C[CH]c1ccccc1)c1ccccc1)c1ccccc1"
PRISTINE = "CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1"
SATURATED_END = PRISTINE

RECOMBINATION_TEMPLATE = (
    "Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_"
    "Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R"
)
INITIATION_RECOMBINATION_TEMPLATE = (
    "Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_"
    "Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R"
)
DISPROPORTIONATION_TEMPLATE = (
    "Root_Ext-1R!H-R_N-4R->O_N-Sp-5R!H=1R!H_Ext-4CHNS-R_"
    "N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClClFFFIIINNNOOOPPPSiSiSi#4C_"
    "6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_2R!H->C"
)

SELECTORS = [
    {
        "label": "initiation",
        "family": "R_Recombination",
        "template": INITIATION_RECOMBINATION_TEMPLATE,
        "reactants": [PRISTINE],
        "products": ["[CH](CCc1ccccc1)c1ccccc1", "[CH2]C(C)c1ccccc1"],
    },
    {
        "label": "depropagation",
        "family": "R_Addition_MultipleBond",
        "template": "Cds-HH_Cds-CbH;CsJ-CbCsH",
        "reactants": [END_RADICAL],
        "products": ["CC(C[CH]c1ccccc1)c1ccccc1", "C=Cc1ccccc1"],
    },
    {
        "label": "termination_recombination",
        "family": "R_Recombination",
        "template": RECOMBINATION_TEMPLATE,
        "reactants": [END_RADICAL, END_RADICAL],
        "products": [
            "CC(CC(CC(c1ccccc1)C(CC(CC(C)c1ccccc1)c1ccccc1)c1ccccc1)"
            "c1ccccc1)c1ccccc1"
        ],
    },
    {
        "label": "termination_disproportionation",
        "family": "Disproportionation",
        "template": DISPROPORTIONATION_TEMPLATE,
        "reactants": [END_RADICAL, END_RADICAL],
        "products": [
            SATURATED_END,
            "CC(CC(C=Cc1ccccc1)c1ccccc1)c1ccccc1",
        ],
    },
]

TRANSFER_PRODUCTS = [
    ("C/H/Cs2Cb;C_rad/H/CbCs", "CC(CC(CCc1ccccc1)=C1C=C[CH]C=C1)c1ccccc1"),
    ("Cb_H;C_rad/H/CbCs", "CC(CC(CCc1[c]cccc1)c1ccccc1)c1ccccc1"),
    ("Cb_H;C_rad/H/CbCs", "CC(CC(CCc1ccccc1)c1ccccc1)c1c[c]ccc1"),
    ("Cb_H;C_rad/H/CbCs", "CC(CC(CCc1ccccc1)c1cc[c]cc1)c1ccccc1"),
    ("Cb_H;C_rad/H/CbCs", "CC(CC(CCc1ccccc1)c1[c]cccc1)c1ccccc1"),
    ("C/H2/NonDeC;C_rad/H/CbCs", "CC([CH]C(CCc1ccccc1)c1ccccc1)c1ccccc1"),
    ("Cb_H;C_rad/H/CbCs", "CC(CC(CCc1ccccc1)c1ccccc1)c1[c]cccc1"),
    ("C/H3/Cs;C_rad/H/CbCs", "[CH2]C(CC(CCc1ccccc1)c1ccccc1)c1ccccc1"),
    ("C/H2/NonDeC;C_rad/H/CbCs", "CC(CC([CH]Cc1ccccc1)c1ccccc1)c1ccccc1"),
    ("C/H/Cs2Cb;C_rad/H/CbCs", "CC(CC(CCc1ccccc1)c1ccccc1)=C1C=C[CH]C=C1"),
    ("Cb_H;C_rad/H/CbCs", "CC(CC(CCc1ccccc1)c1ccccc1)c1cc[c]cc1"),
    ("Cb_H;C_rad/H/CbCs", "CC(CC(CCc1ccccc1)c1c[c]ccc1)c1ccccc1"),
    ("Cb_H;C_rad/H/CbCs", "CC(CC(CCc1cc[c]cc1)c1ccccc1)c1ccccc1"),
    ("Cb_H;C_rad/H/CbCs", "CC(CC(CCc1c[c]ccc1)c1ccccc1)c1ccccc1"),
]

for index, (template, radical_product) in enumerate(TRANSFER_PRODUCTS, start=1):
    SELECTORS.append(
        {
            "label": f"transfer_{index:02d}",
            "family": "H_Abstraction",
            "template": template,
            "reactants": [END_RADICAL, PRISTINE],
            "products": [SATURATED_END, radical_product],
        }
    )


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "artifact",
        type=Path,
        help="kmc_event_set/0.2 JSON artifact to fit",
    )
    parser.add_argument(
        "--database",
        type=Path,
        help=(
            "RMG-database input directory used by the artifact; defaults to "
            "database.directory from rmgrc (validated for provenance only)"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="directory for rate_points.csv (default: script directory)",
    )
    return parser.parse_args()


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@lru_cache(maxsize=None)
def graph_smiles(adjacency_list):
    """Return a canonical, hydrogen-suppressed structural representation."""
    return Molecule().from_adjacency_list(adjacency_list).to_smiles()


def graph_side(graphs):
    return tuple(sorted(graph_smiles(graph) for graph in graphs))


def target_side(smiles):
    # These literals are canonical SMILES emitted from the reference reaction
    # graphs. Keep them literal: parsing aromatic radical SMILES and writing
    # them back can select a Kekule representation instead.
    return tuple(sorted(smiles))


def select_records(records):
    """Return one uniquely matched artifact record for every selector."""
    matches = {}
    for selector in SELECTORS:
        candidates = [
            record
            for record in records
            if record.get("family") == selector["family"]
            and record.get("template") == selector["template"]
            and len(record.get("reactant_graphs", [])) == len(selector["reactants"])
            and len(record.get("product_graphs", [])) == len(selector["products"])
        ]
        expected_reactants = target_side(tuple(selector["reactants"]))
        expected_products = target_side(tuple(selector["products"]))
        structural = [
            record
            for record in candidates
            if graph_side(tuple(record["reactant_graphs"])) == expected_reactants
            and graph_side(tuple(record["product_graphs"])) == expected_products
        ]
        if len(structural) != 1:
            ids = [record.get("event_id", "<missing>") for record in structural]
            raise SystemExit(
                f"Structural selector {selector['label']!r} matched "
                f"{len(structural)} records (event IDs: {ids}); expected exactly one"
            )
        matches[selector["label"]] = structural[0]
    return matches


def rate_table(record, temperatures):
    table = record.get("k_table", {})
    record_temperatures = np.asarray(table.get("T", []), dtype=float)
    if not np.array_equal(record_temperatures, temperatures):
        raise SystemExit(f"Temperature grid mismatch in {record['event_id']}")
    rates = np.asarray(table.get("k", []), dtype=float)
    if rates.shape != temperatures.shape or np.any(rates <= 0.0):
        raise SystemExit(f"Invalid rate table in {record['event_id']}")
    return rates


def fit_channel(name, temperatures, rates, units, event_ids):
    fitted = Arrhenius().fit_to_data(temperatures, rates, units, T0=1.0)
    unconstrained_ea = float(fitted.Ea.value_si)
    constraint = "unconstrained"
    if unconstrained_ea < 0.0:
        design = np.column_stack([np.ones_like(temperatures), np.log(temperatures)])
        log_a, exponent = np.linalg.lstsq(design, np.log(rates), rcond=None)[0]
        fitted = Arrhenius(
            A=(math.exp(log_a), units),
            n=float(exponent),
            Ea=(0.0, "J/mol"),
            T0=(1.0, "K"),
        )
        constraint = f"Ea fixed at zero (unconstrained Ea={unconstrained_ea:.12g})"

    rmg_rates = np.array(
        [fitted.get_rate_coefficient(temperature) for temperature in temperatures]
    )
    qssa_rates = (
        float(fitted.A.value_si)
        * temperatures ** float(fitted.n.value_si)
        * np.exp(-float(fitted.Ea.value_si) / (QSSA_R_GAS * temperatures))
    )
    return {
        "channel": name,
        "A": float(fitted.A.value_si),
        "n": float(fitted.n.value_si),
        "Ea_J_mol": float(fitted.Ea.value_si),
        "units": units,
        "event_ids": sorted(event_ids),
        "rates": rates,
        "rmg_rates": rmg_rates,
        "qssa_rates": qssa_rates,
        "max_rmg_error": float(np.max(np.abs(rmg_rates / rates - 1.0))),
        "max_qssa_error": float(np.max(np.abs(qssa_rates / rates - 1.0))),
        "constraint": constraint,
    }


def write_rate_points(path, artifact_sha, temperatures, channels):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "channel",
                "T_K",
                "k_event_set",
                "k_fit_rmg",
                "k_fit_qssa",
                "relative_error_rmg",
                "relative_error_qssa",
                "event_ids",
                "artifact_sha256",
            ]
        )
        for channel in channels:
            event_ids = ";".join(channel["event_ids"])
            for temperature, observed, rmg_rate, qssa_rate in zip(
                temperatures,
                channel["rates"],
                channel["rmg_rates"],
                channel["qssa_rates"],
            ):
                writer.writerow(
                    [
                        channel["channel"],
                        f"{temperature:.1f}",
                        f"{observed:.12e}",
                        f"{rmg_rate:.12e}",
                        f"{qssa_rate:.12e}",
                        f"{rmg_rate / observed - 1.0:.12e}",
                        f"{qssa_rate / observed - 1.0:.12e}",
                        event_ids,
                        artifact_sha,
                    ]
                )


def main():
    args = parse_arguments()
    artifact_path = args.artifact.expanduser().resolve()
    database_path = (
        args.database or Path(settings["database.directory"])
    ).expanduser().resolve()
    output_directory = args.output_dir.expanduser().resolve()
    if not artifact_path.is_file():
        raise SystemExit(f"Artifact does not exist: {artifact_path}")
    if not (database_path / "thermo").is_dir() or not (
        database_path / "kinetics"
    ).is_dir():
        raise SystemExit(f"Not an RMG-database input directory: {database_path}")

    artifact_sha = file_sha256(artifact_path)
    with artifact_path.open() as handle:
        artifact = json.load(handle)
    if artifact.get("schema_version") != "kmc_event_set/0.2":
        raise SystemExit(
            "Expected kmc_event_set/0.2, got "
            f"{artifact.get('schema_version')!r}"
        )

    selected = select_records(artifact.get("records", []))
    initiation_record = selected["initiation"]
    temperatures = np.asarray(initiation_record["k_table"]["T"], dtype=float)
    if 1000.0 not in temperatures:
        raise SystemExit("Artifact temperature grid does not include deck temperature 1000 K")

    repeat_mw_kg_mol = Molecule().from_smiles(
        PS_REPEAT_SMILES
    ).get_molecular_weight()
    repeat_concentration = PS_DENSITY_KG_M3 / repeat_mw_kg_mol

    transfer_records = [
        selected[f"transfer_{index:02d}"]
        for index in range(1, len(TRANSFER_PRODUCTS) + 1)
    ]
    channels = [
        fit_channel(
            "initiation",
            temperatures,
            2.0 * rate_table(initiation_record, temperatures),
            "s^-1",
            [initiation_record["event_id"]],
        ),
        fit_channel(
            "depropagation",
            temperatures,
            rate_table(selected["depropagation"], temperatures),
            "s^-1",
            [selected["depropagation"]["event_id"]],
        ),
        fit_channel(
            "termination",
            temperatures,
            rate_table(selected["termination_recombination"], temperatures)
            + rate_table(selected["termination_disproportionation"], temperatures),
            "m^3/(mol*s)",
            [
                selected["termination_recombination"]["event_id"],
                selected["termination_disproportionation"]["event_id"],
            ],
        ),
        fit_channel(
            "transfer",
            temperatures,
            sum(rate_table(record, temperatures) for record in transfer_records)
            / 3.0
            * repeat_concentration,
            "s^-1",
            [record["event_id"] for record in transfer_records],
        ),
    ]

    output_directory.mkdir(parents=True, exist_ok=True)
    output_path = output_directory / "rate_points.csv"
    write_rate_points(output_path, artifact_sha, temperatures, channels)

    print("STRUCTURAL_MATCHES")
    print("| selector | family | template | event ID |")
    print("|---|---|---|---|")
    for selector in SELECTORS:
        record = selected[selector["label"]]
        print(
            f"| {selector['label']} | {record['family']} | "
            f"{record['template']} | {record['event_id']} |"
        )
    print("RATE_ESTIMATION_COMPLETE")
    print(f"artifact={artifact_path}")
    print(f"artifact_sha256={artifact_sha}")
    print(f"database={database_path}")
    print(f"output={output_path}")
    print(f"temperature_points={len(temperatures)}")
    print(f"repeat_mw_kg_mol={repeat_mw_kg_mol:.15g}")
    print(f"repeat_concentration_mol_m3={repeat_concentration:.15g}")
    for channel in channels:
        index_1000 = int(np.where(temperatures == 1000.0)[0][0])
        print(
            f"{channel['channel']}: A={channel['A']:.12e} "
            f"n={channel['n']:.12g} Ea={channel['Ea_J_mol']:.12g} J/mol "
            f"k1000_qssa={channel['qssa_rates'][index_1000]:.12e} "
            f"max_rmg_error={100.0 * channel['max_rmg_error']:.6f}% "
            f"max_qssa_error={100.0 * channel['max_qssa_error']:.6f}% "
            f"constraint={channel['constraint']} "
            f"event_ids={','.join(channel['event_ids'])}"
        )


if __name__ == "__main__":
    main()
