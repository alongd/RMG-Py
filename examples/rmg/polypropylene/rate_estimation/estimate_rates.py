"""Reproduce the polypropylene method-of-moments channel fits.

Chemistry comes from an RMG-database input directory supplied by
``--database`` or ``database.directory`` in ``rmgrc``.
"""

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.kinetics import Arrhenius, ArrheniusBM, ArrheniusEP
from rmgpy.molecule import Molecule


TEMPERATURES = np.arange(300.0, 1000.1, 50.0)
PP_DENSITY_KG_M3 = 900.0
REPEAT_MW_KG_MOL = 0.04207974
REPEAT_CONCENTRATION_MOL_M3 = PP_DENSITY_KG_M3 / REPEAT_MW_KG_MOL
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
            "RMG-database input directory; defaults to database.directory "
            "from rmgrc"
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
    args.database or Path(settings["database.directory"])
).expanduser().resolve()
OUT = args.output_dir.expanduser().resolve()
if not (DATABASE / "thermo").is_dir() or not (DATABASE / "kinetics").is_dir():
    raise SystemExit(f"Not an RMG-database input directory: {DATABASE}")
OUT.mkdir(parents=True, exist_ok=True)


def mol(smiles):
    return Molecule().from_smiles(smiles)


def add_thermo(reaction):
    for species in reaction.reactants + reaction.products:
        species.generate_resonance_structures()
        species.thermo = db.thermo.get_thermo_data(species)


def select_model_generation_kinetics(reaction):
    add_thermo(reaction)
    family = db.kinetics.families[reaction.family]
    kinetics, source, entry, is_forward = family.get_kinetics(
        reaction,
        template_labels=reaction.template,
        degeneracy=reaction.degeneracy,
        estimator="rate rules",
        return_all_kinetics=False,
    )
    assert is_forward, (
        f"selected kinetics for {reaction} are defined in the reverse "
        "direction"
    )
    original_type = type(kinetics).__name__
    dHrxn298 = reaction.get_enthalpy_of_reaction(298.0)
    if isinstance(kinetics, (ArrheniusBM, ArrheniusEP)):
        kinetics = kinetics.to_arrhenius(dHrxn298)
    reaction.kinetics = kinetics
    reaction.selected_kinetics = {
        "source": getattr(source, "label", source),
        "entry_label": getattr(entry, "label", None),
        "entry_index": getattr(entry, "index", None),
        "is_forward": is_forward,
        "original_type": original_type,
        "reaction_degeneracy": float(reaction.degeneracy),
        "dHrxn298_J_mol": float(dHrxn298),
        "comment": kinetics.comment,
    }
    return reaction


def generate_with_model_kinetics(reactants, family, products=None):
    reactions = db.kinetics.generate_reactions_from_families(
        [mol(item) for item in reactants],
        None if products is None else [mol(item) for item in products],
        only_families=[family],
    )
    return [
        select_model_generation_kinetics(reaction)
        for reaction in reactions
    ]


def reverse_rates(reaction):
    add_thermo(reaction)
    return np.array(
        [
            reaction.kinetics.get_rate_coefficient(T)
            / reaction.get_equilibrium_constant(T)
            for T in TEMPERATURES
        ]
    )


def forward_rates(reaction):
    return np.array(
        [reaction.kinetics.get_rate_coefficient(T) for T in TEMPERATURES]
    )


def fit(name, rates, units, provenance, require_nonnegative_ea=False):
    fitted = Arrhenius().fit_to_data(TEMPERATURES, rates, units, T0=1.0)
    unconstrained_ea = float(fitted.Ea.value_si)
    if require_nonnegative_ea and unconstrained_ea < 0.0:
        design = np.column_stack(
            [np.ones_like(TEMPERATURES), np.log(TEMPERATURES)]
        )
        log_a, n = np.linalg.lstsq(design, np.log(rates), rcond=None)[0]
        fitted = Arrhenius(
            A=(math.exp(log_a), units),
            n=float(n),
            Ea=(0.0, "J/mol"),
            T0=(1.0, "K"),
        )
    fit_rates = np.array(
        [fitted.get_rate_coefficient(T) for T in TEMPERATURES]
    )
    rel = np.abs(fit_rates / rates - 1.0)
    return {
        "channel": name,
        "A": float(fitted.A.value_si),
        "n": float(fitted.n.value_si),
        "Ea_J_mol": float(fitted.Ea.value_si),
        "A_units": units,
        "temperature_range_K": [300.0, 1000.0],
        "temperature_step_K": 50.0,
        "max_relative_fit_error": float(rel.max()),
        "rms_log_error": float(
            np.sqrt(np.mean(np.log(fit_rates / rates) ** 2))
        ),
        "fit_constraint": (
            "Ea fixed at 0 J/mol because the unconstrained fit was negative "
            f"({unconstrained_ea:.12g} J/mol) and the solver requires Ea >= 0"
            if require_nonnegative_ea and unconstrained_ea < 0.0
            else "unconstrained"
        ),
        "provenance": provenance,
        "rates": [float(value) for value in rates],
        "fit_rates": [float(value) for value in fit_rates],
    }


db = RMGDatabase()
db.load_thermo(
    str(DATABASE / "thermo"),
    [
        "primaryThermoLibrary",
        "thermo_DFT_CCSDTF12_BAC",
        "DFT_QCI_thermo",
        "CBS_QB3_1dHR",
    ],
    depository=False,
    surface=False,
)
db.load_kinetics(
    str(DATABASE / "kinetics"),
    reaction_libraries=[],
    seed_mechanisms=[],
    kinetics_families=FAMILIES,
    kinetics_depositories=["training"],
)
db.load_forbidden_structures(str(DATABASE / "forbiddenStructures.py"))
for family in db.kinetics.families.values():
    if not family.auto_generated:
        family.add_rules_from_training(thermo_database=db.thermo)
        family.fill_rules_by_averaging_up(verbose=True)

# One secondary--tertiary backbone C--C bond in 2,4-dimethylpentane:
# isobutyl + isopropyl -> 2,4-dimethylpentane. Reverse the
# Use the R_Recombination estimate and thermo for per-bond homolysis
# frequency. The solver's moment convention contributes two PP backbone
# bonds per propylene repeat.
initiation_rxn, = generate_with_model_kinetics(
    ["[CH2]C(C)C", "C[CH]C"],
    "R_Recombination",
    ["CC(C)CC(C)C"],
)
initiation_per_bond = reverse_rates(initiation_rxn)
initiation = fit(
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
        "thermo_sources": {
            species.molecule[0].to_smiles(): species.thermo.comment
            for species in initiation_rxn.reactants + initiation_rxn.products
        },
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
deprop_rxn, = generate_with_model_kinetics(
    ["C[CH]C", "C=CC"],
    "R_Addition_MultipleBond",
    ["C[CH]CC(C)C"],
)
depropagation = fit(
    "depropagation",
    reverse_rates(deprop_rxn),
    "s^-1",
    {
        "family": "R_Addition_MultipleBond (thermodynamic reverse)",
        "surrogate": (
            "4-methyl-2-pentyl radical -> propylene + isopropyl radical"
        ),
        "kinetics_source": deprop_rxn.selected_kinetics,
        "thermo_sources": {
            species.molecule[0].to_smiles(): species.thermo.comment
            for species in deprop_rxn.reactants + deprop_rxn.products
        },
        "normalization": (
            "one secondary radical chain end; one propylene per event"
        ),
    },
)

# Chain-end termination: parallel recombination and both distinct
# disproportionation products of two 4-methyl-2-pentyl radicals.
# The QSSA solver accepts their summed kt.
recomb_rxn, = generate_with_model_kinetics(
    ["C[CH]CC(C)C", "C[CH]CC(C)C"], "R_Recombination"
)
disp_rxns = generate_with_model_kinetics(
    ["C[CH]CC(C)C", "C[CH]CC(C)C"], "Disproportionation"
)
assert len(disp_rxns) == 2
termination = fit(
    "termination",
    forward_rates(recomb_rxn)
    + sum((forward_rates(reaction) for reaction in disp_rxns)),
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
inter_rxns_all = generate_with_model_kinetics(
    ["C[CH]CC(C)C", "CC(C)CC(C)CC(C)C"], "H_Abstraction"
)


def has_tertiary_radical(molecule):
    return any(
        atom.radical_electrons == 1
        and sum(
            neighbor.element.symbol == "C" for neighbor in atom.edges
        ) == 3
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


inter_rxns = list(filter(is_tertiary_c10_radical, inter_rxns_all))
assert len(inter_rxns) == 2
assert sum(reaction.degeneracy for reaction in inter_rxns) == 3.0
inter_bimolecular = sum((forward_rates(reaction) for reaction in inter_rxns))
inter_pseudo_first = inter_bimolecular / 3.0 * REPEAT_CONCENTRATION_MOL_M3

# Intramolecular transfer: retain the tertiary-C--H 1,5-shift
# from a PP-trimer-length secondary chain-end radical.
intra_rxns_all = generate_with_model_kinetics(
    ["C[CH]CC(C)CC(C)C"], "intra_H_migration"
)


def is_tertiary_r5h_migration(reaction):
    template = reaction.template[0]
    template_label = getattr(template, "label", str(template))
    return (
        template_label.startswith("R5H")
        and has_tertiary_radical_product(reaction)
    )


intra_rxns = list(filter(is_tertiary_r5h_migration, intra_rxns_all))
assert len(intra_rxns) == 1
intra_first = sum((forward_rates(reaction) for reaction in intra_rxns))
transfer = fit(
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
        "intermolecular_rates_s^-1": [
            float(value) for value in inter_pseudo_first
        ],
        "intramolecular_rates_s^-1": [float(value) for value in intra_first],
    },
)

channels = [initiation, depropagation, termination, transfer]
payload = {
    "database": str(DATABASE),
    "temperatures_K": [float(value) for value in TEMPERATURES],
    "channels": channels,
}
(OUT / "rate_mapping.json").write_text(json.dumps(payload, indent=2) + "\n")

with (OUT / "rate_points.csv").open("w", newline="") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(["channel", "T_K", "k_RMG", "k_fit", "relative_error"])
    for channel in channels:
        for T, observed, fitted in zip(
            TEMPERATURES, channel["rates"], channel["fit_rates"]
        ):
            writer.writerow(
                [
                    channel["channel"],
                    f"{T:.1f}",
                    f"{observed:.12e}",
                    f"{fitted:.12e}",
                    f"{fitted / observed - 1.0:.12e}",
                ]
            )

lines = [
    (
        "| Channel | Family / surrogate | A (SI) | n | Ea (J/mol) | "
        "max fit error |"
    ),
    "|---|---|---:|---:|---:|---:|",
]
for channel in channels:
    provenance = channel["provenance"]
    family = (
        provenance["family"]
        if "family" in provenance
        else " + ".join(provenance["families"])
    )
    surrogate = (
        provenance["surrogate"]
        if "surrogate" in provenance
        else " + ".join(provenance["surrogates"])
    )
    lines.append(
        f"| {channel['channel']} | {family}; {surrogate} | "
        f"{channel['A']:.9e} {channel['A_units']} | {channel['n']:.9g} | "
        f"{channel['Ea_J_mol']:.9f} | "
        f"{100.0 * channel['max_relative_fit_error']:.3f}% |"
    )
(OUT / "rate_table.md").write_text("\n".join(lines) + "\n")

print("RATE_ESTIMATION_COMPLETE")
print(f"database={DATABASE}")
print(f"output_dir={OUT}")
print(f"temperature_points={len(TEMPERATURES)} (300--1000 K, 50 K step)")
print(f"repeat_concentration_mol_m3={REPEAT_CONCENTRATION_MOL_M3:.12g}")
for channel in channels:
    print(
        f"{channel['channel']}: A={channel['A']:.12e} "
        f"n={channel['n']:.12g} Ea={channel['Ea_J_mol']:.12g} J/mol "
        f"max_error={100.0 * channel['max_relative_fit_error']:.6f}%"
    )
