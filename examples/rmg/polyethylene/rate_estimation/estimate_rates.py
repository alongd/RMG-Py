"""Reproduce the polyethylene method-of-moments channel fits.

All chemistry is evaluated from an RMG-database input directory supplied by
``--database`` or configured as ``database.directory`` in ``rmgrc``.
"""

import argparse
from pathlib import Path

from rmgpy import settings
from examples.rmg.polymer_rate_estimation import (
    ProductionRateEstimator,
    fit_channel,
    print_summary,
    reaction_thermo_sources,
    repeat_properties,
    write_artifacts,
)


PE_DENSITY_KG_M3 = 950.0
REPEAT_MW_KG_MOL, REPEAT_CONCENTRATION_MOL_M3 = repeat_properties(
    "[CH2][CH2]", PE_DENSITY_KG_M3
)
FAMILIES = [
    "R_Recombination",
    "Disproportionation",
    "R_Addition_MultipleBond",
    "H_Abstraction",
    "intra_H_migration",
]


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--database",
        type=Path,
        help=(
            "RMG-database input directory; defaults to database.directory " "from rmgrc"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help=(
            "directory for generated CSV, JSON, and Markdown "
            "(default: script directory)"
        ),
    )
    return parser.parse_args()


args = parse_arguments()
DATABASE = (
    (args.database or Path(settings["database.directory"])).expanduser().resolve()
)
OUT = args.output_dir.expanduser().resolve()
if not (DATABASE / "thermo").is_dir() or not (DATABASE / "kinetics").is_dir():
    raise SystemExit(f"Not an RMG-database input directory: {DATABASE}")
OUT.mkdir(parents=True, exist_ok=True)


estimator = ProductionRateEstimator(DATABASE, FAMILIES)

# One central C--C bond in n-hexane: n-C3H7 + n-C3H7 -> n-C6H14.
# Reverse the R_Recombination estimate with RMG thermo to obtain a per-bond
# homolysis frequency. The solver's moment convention contributes two PE
# backbone bonds per ethylene repeat.
(initiation_rxn,) = estimator.generate(
    ["[CH2]CC", "[CH2]CC"], "R_Recombination", ["CCCCCC"]
)
initiation_per_bond = initiation_rxn.requested_reverse_rates()
initiation = fit_channel(
    "initiation",
    2.0 * initiation_per_bond,
    "s^-1",
    {
        "family": "R_Recombination (thermodynamic reverse)",
        "surrogate": "n-hexane -> 2 n-propyl radicals (central C--C bond)",
        "kinetics_source": initiation_rxn.selected_kinetics,
        "thermo_sources": reaction_thermo_sources(initiation_rxn),
        "normalization": (
            "one central-bond homolysis rate multiplied by two PE backbone "
            "C--C bonds per ethylene repeat, because the solver applies "
            "initiation to mu1-mu0 repeat-bond units"
        ),
        "per_bond_rates_s^-1": [float(value) for value in initiation_per_bond],
        "backbone_bonds_per_repeat": 2.0,
    },
)

# Chain-end beta-scission: reverse n-butyl addition to ethylene.
(deprop_rxn,) = estimator.generate(
    ["[CH2]CCC", "C=C"],
    "R_Addition_MultipleBond",
    ["[CH2]CCCCC"],
)
depropagation = fit_channel(
    "depropagation",
    deprop_rxn.requested_reverse_rates(),
    "s^-1",
    {
        "family": "R_Addition_MultipleBond (thermodynamic reverse)",
        "surrogate": "1-hexyl radical -> ethylene + 1-butyl radical",
        "kinetics_source": deprop_rxn.selected_kinetics,
        "thermo_sources": reaction_thermo_sources(deprop_rxn),
        "normalization": "one radical chain end; one ethylene per event",
    },
)

# Chain-end termination: parallel recombination and disproportionation of two
# primary 1-hexyl radicals. The QSSA solver accepts their summed kt.
(recomb_rxn,) = estimator.generate(["[CH2]CCCCC", "[CH2]CCCCC"], "R_Recombination")
(disp_rxn,) = estimator.generate(["[CH2]CCCCC", "[CH2]CCCCC"], "Disproportionation")
termination = fit_channel(
    "termination",
    recomb_rxn.requested_forward_rates() + disp_rxn.requested_forward_rates(),
    "m^3/(mol*s)",
    {
        "families": ["R_Recombination", "Disproportionation"],
        "surrogates": [
            "2 1-hexyl radicals -> n-dodecane",
            "2 1-hexyl radicals -> n-hexane + 1-hexene",
        ],
        "kinetics_sources": [
            recomb_rxn.selected_kinetics,
            disp_rxn.selected_kinetics,
        ],
        "normalization": "sum of the two bimolecular disappearance channels",
    },
    require_nonnegative_ea=True,
)

# Intermolecular transfer: retain the three secondary-octyl abstraction
# products. Together they span three ethylene-repeat equivalents; convert the
# summed bimolecular rate to the solver's pseudo-first-order convention.
inter_rxns = estimator.generate(
    ["[CH2]CCCCC", "CCCCCCCC"],
    "H_Abstraction",
    predicate=lambda reaction: reaction.degeneracy == 4.0,
)
assert len(inter_rxns) == 3
inter_bimolecular = sum(reaction.requested_forward_rates() for reaction in inter_rxns)
inter_pseudo_first = inter_bimolecular / 3.0 * REPEAT_CONCENTRATION_MOL_M3

# Intramolecular transfer: retain all six distinct 1-octyl H-migration paths.
intra_rxns = estimator.generate(["[CH2]CCCCCCC"], "intra_H_migration")
assert len(intra_rxns) == 6
intra_first = sum(reaction.requested_forward_rates() for reaction in intra_rxns)
transfer = fit_channel(
    "transfer",
    intra_first + inter_pseudo_first,
    "s^-1",
    {
        "families": ["H_Abstraction", "intra_H_migration"],
        "surrogates": [
            "1-hexyl + n-octane -> n-hexane + secondary 2/3/4-octyl",
            "1-octyl -> secondary 2/3/4-octyl by all six distinct H-migration paths",
        ],
        "kinetics_sources": [
            reaction.selected_kinetics for reaction in inter_rxns + intra_rxns
        ],
        "normalization": (
            "intermolecular sum / 3 interior ethylene-repeat equivalents * "
            f"{REPEAT_CONCENTRATION_MOL_M3:.12g} mol/m^3 repeat concentration; "
            "then add the unimolecular intra_H_migration sum"
        ),
        "repeat_concentration_mol_m3": REPEAT_CONCENTRATION_MOL_M3,
        "density_kg_m3": PE_DENSITY_KG_M3,
        "repeat_mw_kg_mol": REPEAT_MW_KG_MOL,
        "intermolecular_rates_s^-1": [float(value) for value in inter_pseudo_first],
        "intramolecular_rates_s^-1": [float(value) for value in intra_first],
    },
)

channels = [initiation, depropagation, termination, transfer]
write_artifacts(OUT, DATABASE, channels)
print_summary(DATABASE, OUT, REPEAT_CONCENTRATION_MOL_M3, channels)
