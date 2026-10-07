"""Shared production-path helpers for polymer MOM example rate fits."""

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from rmgpy.data.rmg import RMGDatabase
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.rmg.model import CoreEdgeReactionModel


TEMPERATURES = np.arange(300.0, 1000.1, 50.0)
QSSA_R_GAS = 8.314


def molecule(smiles):
    """Return an RMG molecule parsed from ``smiles``."""
    return Molecule().from_smiles(smiles)


def repeat_properties(repeat_smiles, density_kg_m3):
    """Return exact RMG repeat MW and the initial melt concentration."""
    repeat_mw_kg_mol = molecule(repeat_smiles).get_molecular_weight()
    return repeat_mw_kg_mol, density_kg_m3 / repeat_mw_kg_mol


def _side_matches_molecules(species, molecules):
    if len(species) != len(molecules):
        return False

    def match(remaining_species, remaining_molecules):
        if not remaining_species:
            return True
        candidate = remaining_species[0]
        for index, target in enumerate(remaining_molecules):
            if candidate.is_isomorphic(target):
                if match(
                    remaining_species[1:],
                    remaining_molecules[:index] + remaining_molecules[index + 1 :],
                ):
                    return True
        return False

    return match(list(species), list(molecules))


def _template_labels(reaction, requested_forward):
    template = reaction.template
    if not requested_forward and getattr(reaction, "reverse", None):
        template = reaction.reverse.template
    return [getattr(item, "label", str(item)) for item in template]


def _requested_view(reaction, requested_reactants):
    if _side_matches_molecules(reaction.reactants, requested_reactants):
        return RequestedReactionView(
            reactants=reaction.reactants,
            products=reaction.products,
            degeneracy=float(reaction.degeneracy),
            template_labels=_template_labels(reaction, True),
        )
    if _side_matches_molecules(reaction.products, requested_reactants):
        reverse = getattr(reaction, "reverse", None)
        degeneracy = reverse.degeneracy if reverse else reaction.degeneracy
        return RequestedReactionView(
            reactants=reaction.products,
            products=reaction.reactants,
            degeneracy=float(degeneracy),
            template_labels=_template_labels(reaction, False),
        )
    raise ValueError(
        "Generated reaction cannot be related to the requested physical "
        f"reactants: {reaction}"
    )


class _RecordingReactionModel(CoreEdgeReactionModel):
    """Production model that records the selection made by its own seam."""

    def __init__(self):
        super().__init__()
        self.last_selection = None

    def generate_kinetics(self, reaction):
        selection = super().generate_kinetics(reaction)
        kinetics, source, entry, is_forward = selection
        self.last_selection = {
            "source": getattr(source, "label", str(source)),
            "entry_label": getattr(entry, "label", None),
            "entry_index": getattr(entry, "index", None),
            "entry_rank": getattr(entry, "rank", None),
            "selection_is_forward": bool(is_forward),
            "raw_kinetics_type": type(kinetics).__name__,
            "raw_kinetics_comment": kinetics.comment,
        }
        return selection


@dataclass(frozen=True)
class RequestedReactionView:
    """Pre-production chemistry oriented as requested by the caller."""

    reactants: list
    products: list
    degeneracy: float
    template_labels: list


@dataclass
class ProductionReaction:
    """A production-selected reaction plus its immutable physical direction."""

    reaction: object
    requested_reactant_smiles: list
    requested_product_smiles: list
    requested_forward_in_final: bool
    selected_kinetics: dict

    @property
    def requested_degeneracy(self):
        """Return the reaction-path degeneracy before canonicalization."""
        return float(self.selected_kinetics["requested_degeneracy"])

    def requested_forward_rates(self, temperatures=TEMPERATURES):
        rates = self._final_forward_rates(temperatures)
        if self.requested_forward_in_final:
            return rates
        return rates / self._equilibrium_constants(temperatures)

    def requested_reverse_rates(self, temperatures=TEMPERATURES):
        rates = self._final_forward_rates(temperatures)
        if self.requested_forward_in_final:
            return rates / self._equilibrium_constants(temperatures)
        return rates

    def _final_forward_rates(self, temperatures):
        return np.array(
            [self.reaction.kinetics.get_rate_coefficient(T) for T in temperatures]
        )

    def _equilibrium_constants(self, temperatures):
        return np.array(
            [self.reaction.get_equilibrium_constant(T, type="Kc") for T in temperatures]
        )


class ProductionRateEstimator:
    """Deep example module for production reaction selection and rates."""

    def __init__(self, database_path, families):
        self.database_path = Path(database_path).expanduser().resolve()
        if not (self.database_path / "thermo").is_dir():
            raise ValueError(f"Missing thermo database: {self.database_path}")
        if not (self.database_path / "kinetics").is_dir():
            raise ValueError(f"Missing kinetics database: {self.database_path}")

        self.database = RMGDatabase()
        self.database.load_thermo(
            str(self.database_path / "thermo"),
            [
                "primaryThermoLibrary",
                "thermo_DFT_CCSDTF12_BAC",
                "DFT_QCI_thermo",
                "CBS_QB3_1dHR",
            ],
            depository=False,
            surface=False,
        )
        self.database.load_kinetics(
            str(self.database_path / "kinetics"),
            reaction_libraries=[],
            seed_mechanisms=[],
            kinetics_families=families,
            kinetics_depositories=["training"],
        )
        self.database.load_forbidden_structures(
            str(self.database_path / "forbiddenStructures.py")
        )
        for family in self.database.kinetics.families.values():
            if not family.auto_generated:
                family.add_rules_from_training(thermo_database=self.database.thermo)
                family.fill_rules_by_averaging_up(verbose=True)

        self.model = _RecordingReactionModel()
        self.model.kinetics_estimator = "rate rules"
        self.model.pressure_dependence = True

    def generate(
        self,
        reactant_smiles,
        family,
        product_smiles=None,
        predicate=None,
    ):
        """Generate and production-process reactions in the requested direction."""
        requested_reactants = [molecule(item) for item in reactant_smiles]
        requested_products = (
            None
            if product_smiles is None
            else [molecule(item) for item in product_smiles]
        )
        raw_reactions = self.database.kinetics.generate_reactions_from_families(
            requested_reactants,
            requested_products,
            only_families=[family],
        )

        selected = []
        for raw_reaction in raw_reactions:
            view = _requested_view(raw_reaction, requested_reactants)
            if predicate is not None and not predicate(view):
                continue
            pre_production = {
                "requested_degeneracy": view.degeneracy,
                "requested_template_labels": view.template_labels,
            }
            reaction, is_new = self.model.make_new_reaction(
                raw_reaction,
                check_existing=False,
                generate_thermo=True,
                generate_kinetics=True,
                perform_cut=False,
            )
            if reaction is None or not is_new:
                raise RuntimeError(f"Production rejected surrogate: {raw_reaction}")

            if _side_matches_molecules(reaction.reactants, requested_reactants):
                requested_forward = True
                physical_products = reaction.products
            elif _side_matches_molecules(reaction.products, requested_reactants):
                requested_forward = False
                physical_products = reaction.reactants
            else:
                raise RuntimeError(
                    "Production reaction lost the requested physical direction: "
                    f"{reaction}"
                )

            if requested_products is not None and not _side_matches_molecules(
                physical_products, requested_products
            ):
                raise RuntimeError(
                    "Production reaction does not match requested products: "
                    f"{reaction}"
                )

            selection = dict(self.model.last_selection)
            selection.update(pre_production)
            selection.update(
                {
                    "requested_forward_in_final": requested_forward,
                    "final_degeneracy": float(reaction.degeneracy),
                    "final_template_labels": [
                        getattr(item, "label", str(item)) for item in reaction.template
                    ],
                    "final_kinetics_type": type(reaction.kinetics).__name__,
                    "final_kinetics_comment": reaction.kinetics.comment,
                    "final_reaction": str(reaction),
                }
            )
            selected.append(
                ProductionReaction(
                    reaction=reaction,
                    requested_reactant_smiles=list(reactant_smiles),
                    requested_product_smiles=[
                        species.molecule[0].to_smiles() for species in physical_products
                    ],
                    requested_forward_in_final=requested_forward,
                    selected_kinetics=selection,
                )
            )
        return selected


def fit_channel(
    name,
    rates,
    units,
    provenance,
    require_nonnegative_ea=False,
):
    """Fit one channel and also evaluate the solver's pinned gas constant."""
    fitted = Arrhenius().fit_to_data(TEMPERATURES, rates, units, T0=1.0)
    unconstrained_ea = float(fitted.Ea.value_si)
    if require_nonnegative_ea and unconstrained_ea < 0.0:
        design = np.column_stack([np.ones_like(TEMPERATURES), np.log(TEMPERATURES)])
        log_a, exponent = np.linalg.lstsq(design, np.log(rates), rcond=None)[0]
        fitted = Arrhenius(
            A=(math.exp(log_a), units),
            n=float(exponent),
            Ea=(0.0, "J/mol"),
            T0=(1.0, "K"),
        )

    fit_rates = np.array([fitted.get_rate_coefficient(T) for T in TEMPERATURES])
    qssa_fit_rates = (
        float(fitted.A.value_si)
        * TEMPERATURES ** float(fitted.n.value_si)
        * np.exp(-float(fitted.Ea.value_si) / (QSSA_R_GAS * TEMPERATURES))
    )
    relative_fit_error = fit_rates / rates - 1.0
    relative_qssa_error = qssa_fit_rates / rates - 1.0
    return {
        "channel": name,
        "A": float(fitted.A.value_si),
        "n": float(fitted.n.value_si),
        "Ea_J_mol": float(fitted.Ea.value_si),
        "A_units": units,
        "temperature_range_K": [300.0, 1000.0],
        "temperature_step_K": 50.0,
        "max_relative_fit_error": float(np.abs(relative_fit_error).max()),
        "max_relative_qssa_law_error": float(np.abs(relative_qssa_error).max()),
        "qssa_error_at_1000K": float(relative_qssa_error[-1]),
        "rms_log_error": float(np.sqrt(np.mean(np.log(fit_rates / rates) ** 2))),
        "fit_constraint": (
            "Ea fixed at 0 J/mol because the unconstrained fit was negative "
            f"({unconstrained_ea:.12g} J/mol) and the solver requires Ea >= 0"
            if require_nonnegative_ea and unconstrained_ea < 0.0
            else "unconstrained"
        ),
        "provenance": provenance,
        "rates": [float(value) for value in rates],
        "fit_rates": [float(value) for value in fit_rates],
        "qssa_fit_rates": [float(value) for value in qssa_fit_rates],
    }


def reaction_thermo_sources(production_reaction):
    """Return production-processed thermo comments for one surrogate."""
    reaction = production_reaction.reaction
    return {
        species.molecule[0].to_smiles(): species.thermo.comment
        for species in reaction.reactants + reaction.products
    }


def write_artifacts(output_directory, database_path, channels):
    """Write the common JSON, CSV, and Markdown fit artifacts."""
    output_directory = Path(output_directory).expanduser().resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    payload = {
        "database": str(Path(database_path).expanduser().resolve()),
        "database_revision": (
            "declared cd86d4e1c; unavailable from snapshot (no .git metadata)"
        ),
        "kinetics_path": "CoreEdgeReactionModel.make_new_reaction",
        "temperatures_K": [float(value) for value in TEMPERATURES],
        "channels": channels,
    }
    (output_directory / "rate_mapping.json").write_text(
        json.dumps(payload, indent=2) + "\n"
    )

    with (output_directory / "rate_points.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "channel",
                "T_K",
                "k_RMG",
                "k_fit",
                "relative_error",
                "k_QSSA",
                "relative_QSSA_error",
            ]
        )
        for channel in channels:
            for temperature, observed, fitted, qssa_fitted in zip(
                TEMPERATURES,
                channel["rates"],
                channel["fit_rates"],
                channel["qssa_fit_rates"],
            ):
                writer.writerow(
                    [
                        channel["channel"],
                        f"{temperature:.1f}",
                        f"{observed:.12e}",
                        f"{fitted:.12e}",
                        f"{fitted / observed - 1.0:.12e}",
                        f"{qssa_fitted:.12e}",
                        f"{qssa_fitted / observed - 1.0:.12e}",
                    ]
                )

    lines = [
        (
            "| Channel | Family / surrogate | A (SI) | n | Ea (J/mol) | "
            "max RMG fit error | max solver-law error |"
        ),
        "|---|---|---:|---:|---:|---:|---:|",
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
            f"{channel['A']:.9e} {channel['A_units']} | "
            f"{channel['n']:.9g} | {channel['Ea_J_mol']:.9f} | "
            f"{100.0 * channel['max_relative_fit_error']:.3f}% | "
            f"{100.0 * channel['max_relative_qssa_law_error']:.3f}% |"
        )
    (output_directory / "rate_table.md").write_text("\n".join(lines) + "\n")


def print_summary(database_path, output_directory, repeat_concentration, channels):
    """Print the common completion summary."""
    print("RATE_ESTIMATION_COMPLETE")
    print(f"database={Path(database_path).expanduser().resolve()}")
    print(f"output_dir={Path(output_directory).expanduser().resolve()}")
    print(f"temperature_points={len(TEMPERATURES)} (300--1000 K, 50 K step)")
    print(f"repeat_concentration_mol_m3={repeat_concentration:.12g}")
    for channel in channels:
        print(
            f"{channel['channel']}: A={channel['A']:.12e} "
            f"n={channel['n']:.12g} Ea={channel['Ea_J_mol']:.12g} J/mol "
            f"max_error={100.0 * channel['max_relative_fit_error']:.6f}% "
            "max_qssa_law_error="
            f"{100.0 * channel['max_relative_qssa_law_error']:.6f}%"
        )
