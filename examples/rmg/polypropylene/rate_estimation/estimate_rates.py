"""Reproduce the polypropylene method-of-moments channel fits.

Chemistry comes from an RMG-database input directory supplied by
``--database`` or ``database.directory`` in ``rmgrc``.
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


PP_DENSITY_KG_M3 = 900.0
REPEAT_MW_KG_MOL, REPEAT_CONCENTRATION_MOL_M3 = repeat_properties(
    "[CH2][CH]C", PP_DENSITY_KG_M3
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

# One secondary--tertiary backbone C--C bond in 2,4-dimethylpentane:
# isobutyl + isopropyl -> 2,4-dimethylpentane. Use the reverse of the
# R_Recombination estimate and thermo for per-bond homolysis
# frequency. The solver's moment convention contributes two PP backbone
# bonds per propylene repeat.
(initiation_rxn,) = estimator.generate(
    ["[CH2]C(C)C", "C[CH]C"],
    "R_Recombination",
    ["CC(C)CC(C)C"],
)
initiation_per_bond = initiation_rxn.requested_reverse_rates()
initiation = fit_channel(
    "initiation",
    2.0 * initiation_per_bond,
    "s^-1",
    {
        "family": "R_Recombination (thermodynamic reverse)",
        "surrogate": (
            "2,4-dimethylpentane -> isobutyl + isopropyl radicals "
            "(one secondary--tertiary backbone C--C bond)"
        ),
        "kinetics_source": initiation_rxn.selected_kinetics,
        "thermo_sources": reaction_thermo_sources(initiation_rxn),
        "normalization": (
            "one secondary--tertiary bond homolysis rate multiplied by two "
            "PP backbone C--C bonds per propylene repeat, because the "
            "solver applies initiation to mu1-mu0 repeat-bond units"
        ),
        "per_bond_rates_s^-1": [float(value) for value in initiation_per_bond],
        "backbone_bonds_per_repeat": 2.0,
    },
)

# Chain-end beta-scission: reverse isopropyl addition to propylene at
# terminal CH2. This head-to-tail route produces a secondary PP
# chain-end radical and releases propylene in reverse.
(deprop_rxn,) = estimator.generate(
    ["C[CH]C", "C=CC"],
    "R_Addition_MultipleBond",
    ["C[CH]CC(C)C"],
)
depropagation = fit_channel(
    "depropagation",
    deprop_rxn.requested_reverse_rates(),
    "s^-1",
    {
        "family": "R_Addition_MultipleBond (thermodynamic reverse)",
        "surrogate": ("4-methyl-2-pentyl radical -> propylene + isopropyl radical"),
        "kinetics_source": deprop_rxn.selected_kinetics,
        "thermo_sources": reaction_thermo_sources(deprop_rxn),
        "normalization": ("one secondary radical chain end; one propylene per event"),
    },
)

# Chain-end termination: parallel recombination and both distinct
# disproportionation products of two 4-methyl-2-pentyl radicals.
# The QSSA solver accepts their summed kt.
(recomb_rxn,) = estimator.generate(["C[CH]CC(C)C", "C[CH]CC(C)C"], "R_Recombination")
disp_rxns = estimator.generate(["C[CH]CC(C)C", "C[CH]CC(C)C"], "Disproportionation")
assert len(disp_rxns) == 2
termination = fit_channel(
    "termination",
    recomb_rxn.requested_forward_rates()
    + sum(reaction.requested_forward_rates() for reaction in disp_rxns),
    "m^3/(mol*s)",
    {
        "families": ["R_Recombination", "Disproportionation"],
        "surrogates": [
            "2 4-methyl-2-pentyl radicals -> recombination product",
            (
                "2 4-methyl-2-pentyl radicals -> 2-methylpentane + "
                "the two distinct hexene products"
            ),
        ],
        "kinetics_sources": [
            recomb_rxn.selected_kinetics,
            *[reaction.selected_kinetics for reaction in disp_rxns],
        ],
        "normalization": (
            "sum of recombination and both disproportionation rows, each with "
            "the reaction-path degeneracy assigned by RMG"
        ),
    },
    require_nonnegative_ea=True,
)


# Intermolecular transfer: retain tertiary C--H abstraction from the
# three PP-like sites in 2,4,6-trimethylheptane. Convert the sum
# bimolecular coefficient to the solver's pseudo-first-order convention.
def has_tertiary_radical(molecule):
    return any(
        atom.radical_electrons == 1
        and sum(neighbor.element.symbol == "C" for neighbor in atom.edges) == 3
        for atom in molecule.atoms
    )


def has_tertiary_radical_product(reaction, carbon_count=None):
    for species in reaction.products:
        molecule = species.molecule[0]
        molecule_carbon_count = sum(
            atom.element.symbol == "C" for atom in molecule.atoms
        )
        if carbon_count is not None and molecule_carbon_count != carbon_count:
            continue
        if has_tertiary_radical(molecule):
            return True
    return False


def is_tertiary_c10_radical(reaction):
    return has_tertiary_radical_product(reaction, carbon_count=10)


inter_rxns = estimator.generate(
    ["C[CH]CC(C)C", "CC(C)CC(C)CC(C)C"],
    "H_Abstraction",
    predicate=is_tertiary_c10_radical,
)
assert len(inter_rxns) == 2
assert sum(reaction.requested_degeneracy for reaction in inter_rxns) == 3.0
inter_bimolecular = sum(reaction.requested_forward_rates() for reaction in inter_rxns)
inter_pseudo_first = inter_bimolecular / 3.0 * REPEAT_CONCENTRATION_MOL_M3


# Intramolecular transfer: retain the tertiary-C--H 1,5-shift
# from a PP-trimer-length secondary chain-end radical.
def is_tertiary_r5h_migration(reaction):
    return reaction.template_labels[0].startswith(
        "R5H"
    ) and has_tertiary_radical_product(reaction)


intra_rxns = estimator.generate(
    ["C[CH]CC(C)CC(C)C"],
    "intra_H_migration",
    predicate=is_tertiary_r5h_migration,
)
assert len(intra_rxns) == 1
intra_first = sum(reaction.requested_forward_rates() for reaction in intra_rxns)
transfer = fit_channel(
    "transfer",
    intra_first + inter_pseudo_first,
    "s^-1",
    {
        "families": ["H_Abstraction", "intra_H_migration"],
        "surrogates": [
            (
                "4-methyl-2-pentyl + 2,4,6-trimethylheptane -> "
                "2-methylpentane + tertiary PP-like radical"
            ),
            (
                "PP-trimer-length secondary chain-end radical -> tertiary "
                "mid-chain radical by intramolecular 1,5 H shift"
            ),
        ],
        "kinetics_sources": [
            reaction.selected_kinetics for reaction in inter_rxns + intra_rxns
        ],
        "normalization": (
            "tertiary-C--H intermolecular sum / 3 propylene-repeat "
            f"equivalents * {REPEAT_CONCENTRATION_MOL_M3:.12g} mol/m^3 "
            "repeat concentration; "
            "then add the unimolecular tertiary 1,5 intra_H_migration row"
        ),
        "repeat_concentration_mol_m3": REPEAT_CONCENTRATION_MOL_M3,
        "density_kg_m3": PP_DENSITY_KG_M3,
        "repeat_mw_kg_mol": REPEAT_MW_KG_MOL,
        "intermolecular_rates_s^-1": [float(value) for value in inter_pseudo_first],
        "intramolecular_rates_s^-1": [float(value) for value in intra_first],
    },
)

channels = [initiation, depropagation, termination, transfer]
write_artifacts(OUT, DATABASE, channels)
print_summary(DATABASE, OUT, REPEAT_CONCENTRATION_MOL_M3, channels)
