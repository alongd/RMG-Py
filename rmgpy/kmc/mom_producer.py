"""Method-of-moments producer for the melt kMC runner's ``--mom`` file (I-169).

Command line::

    python -m rmgpy.kmc.mom_producer DECK --mom-checkout PATH --temperature T
        --end E --snapshots N --families F [F ...] --catalog FILE [FILE ...] --out DIR
        [--run | --write-only | --convert] [--allow-family-mismatch]
        [--database DIR] [--species-constraint KEY=VALUE ...]

The same deck the kMC runner reads (``rmgpy.kmc.run.read_deck``: one
``polymer(..., linear_homopolymer_spec=...)``) is written as an RMG
method-of-moments (MoM) input for a polymer-mom checkout, whose ``polymer()``
takes legacy arguments only.  The translation, and every MoM setting the deck
does not determine, is a named convention in ``<out>/conventions.json``
(``CONVENTIONS``, written in every mode).

**Chemistry scope.**  Only the RMG families given (``database(kineticsFamilies=
[...])``; no reaction libraries, no seed mechanisms) and none of the MoM moment
kernels (``k_scission``, ``k_unzip``, ``k_homolysis``, ``k_depropagation``,
``radical_qssa_unzip``, ``side_group_homolysis``, ``explicit_dp``).  The
families must equal the union of the kMC catalogs' ``provenance.families``;
otherwise the producer refuses unless ``--allow-family-mismatch``, and then the
mismatch is recorded in ``conventions.json`` and ``producer.json``.

**Modes.**  ``--write-only`` (default) writes ``<out>/rmg/input.py``,
``<out>/rmg/rmgrc`` (the kinetics database: ``--database`` or the catalogs'
common ``provenance.database_path``), ``conventions.json`` and
``producer.json``.  ``--run`` also runs RMG from the checkout in a subprocess
(``PYTHONPATH`` = the checkout, which must be the ``rmgpy`` imported; both
streams to ``<out>/rmg/stdout.log`` / ``stderr.log``) and converts.
``--convert`` converts an existing ``<out>/rmg`` run.

**Conversion.**  The final-core solver profile ``rmg/solver/simulation_1_<N>.csv``
(highest ``N``; columns are mol) and the pool table
``rmg/chemkin/polymer_pools.json`` give, at exactly the deck's observation
times (``rmgpy.kmc.run.observation_times``), the quantities of
``QUANTITY_SOURCES`` keyed like ``melt_observables.mom_reported_quantities``.
Each profile column is interpolated linearly in time between the two profile
rows that bracket the requested time (a time on a row takes the row); derived
quantities are computed from the interpolated columns; a requested time
outside the profile is refused.  ``NOT_PRODUCED`` keys are absent, with the
reason.  Output: ``<out>/mom.json`` (the ``--mom`` file, a JSON list of
``{"time_s", "quantities"}``) and ``<out>/mom_diagnostics.json`` (the MoM's own
code-convention numbers at the same times: per-pool moments, its end-group-free
condensed mass and Mn / Mw).
"""

import argparse
import ast
import bisect
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from rmgpy import constants
from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import melt_observables as mo
from rmgpy.kmc import run as kr
from rmgpy.polymer_mass import observables_from_moments

SCHEMA = "kmc-mom-producer/1"
CONVENTIONS_SCHEMA = "kmc-mom-conventions/1"
BATH_GAS_TOTAL_MOL = 1.0
DEFAULT_PRESSURE = (1.0, "bar")
DEFAULT_THERMO_LIBRARIES = ("primaryThermoLibrary",)
LABEL_MAX = 10
MASS_RELATIVE = 1.0e-9
BUILTIN_STRUCTURES = {"H": "[H]", "H2": "[H][H]"}

QUANTITY_SOURCES = {
    "mu0_mol": "sum over pools p of profile column <p>_mu0 (mol of chains)",
    "mu1_mol": "sum over pools p of profile column <p>_mu1 (mol of repeat units)",
    "mu2_mol": "sum over pools p of profile column <p>_mu2 (mol x DP^2)",
    "Mn_g_mol": "observables_from_moments(sum mu, a, e).mn_g_mol (every chain a n + e)",
    "Mw_g_mol": "observables_from_moments(sum mu, a, e).mw_g_mol (every chain a n + e)",
    "polymer_condensed_mass_g": "sum_p (M_u,p mu1_p - defect_p mu0_p) + e sum_p mu0_p",
    "polymer_conversion": "1 - polymer_condensed_mass_g(t) / polymer_condensed_mass_g(profile t = 0)",
    "polymer_derived_molecules_mol": "sum over pools p of <p>_mu0 (no explicit-DP species)",
    "species_mol": "profile column of each declared species (declared as a MoM input species)",
}
NOT_PRODUCED = {
    "end_inventory_mol": (
        "MoM's end-group ledger L_p = mu0_p(0) - mu0_p - H_p (polymer.pyx end_group_ledger) nets "
        "chains minted without end groups (scission) against chains removed with them "
        "(recombination, release to gas); the tail's physical end-group inventory is not "
        "determined, and the ledger is written to no file"),
    "explicit_dp_mol": (
        "no explicit-DP species: explicit_dp=False (the deck requests none; MoM's explicit_dp is a "
        "single DP == cutoff species, the kMC ladder 1..xs has no MoM counterpart)"),
    "released_mol": (
        "per-channel side-group release is in no MoM output (the depletion inventory Z is an "
        "auxiliary state outside the profile; gas products are shared by channels) and the "
        "side-group kernel is off (chemistry from the families only)"),
}


class ProducerError(ValueError):
    """An input of the producer is refused, or the MoM run cannot be converted."""


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


# --------------------------------------------------------------------------
# the runner's quantity keys
# --------------------------------------------------------------------------
def runner_quantity_keys() -> frozenset:
    """The top-level keys ``melt_observables.mom_reported_quantities`` returns
    (evaluated on a minimal observation, so the runner's own definition is
    the reference)."""
    obs = {"tier1": {"polymer_derived": {"mass_g": 1.0, "mom_definition": {"amount_mol": 1.0}},
                     "species": {}},
           "tier3": {"tail": {"dp_moments_mol": [1.0, 2.0, 4.0]}, "explicit": {"dp_amount_mol": {}}},
           "mom_mass_law": {"end_inventory_mol": 1.0},
           "repeat_unit_molar_mass_g_mol": 1.0, "end_group_molar_mass_g_mol": 1.0,
           "initial_polymer_mass_g": 1.0}
    return frozenset(mo.mom_reported_quantities(obs))


def check_quantity_keys(quantities: dict, declared: Sequence[str]) -> None:
    """Refuse a quantity the runner does not report, or a species it does not declare."""
    known = runner_quantity_keys()
    unknown = sorted(set(quantities) - known)
    if unknown:
        raise ProducerError(f"quantities {unknown} are not keys of melt_observables.mom_reported_quantities "
                            f"({sorted(known)})")
    clash = sorted(set(quantities) & set(NOT_PRODUCED))
    if clash:
        raise ProducerError(f"quantities {clash} are declared not produced (NOT_PRODUCED)")
    stray = sorted(set(quantities.get("species_mol", {})) - set(declared))
    if stray:
        raise ProducerError(f"species {stray} are not declared species of the deck ({list(declared)})")


# --------------------------------------------------------------------------
# names
# --------------------------------------------------------------------------
def mom_label(label: str) -> str:
    """The MoM pool label: the deck label's letters and digits, at most
    ``LABEL_MAX`` characters, starting with a letter ("P" prefixed when not;
    "POLY" when nothing is left)."""
    text = re.sub(r"[^A-Za-z0-9]", "", label)
    if not text:
        text = "POLY"
    if not text[0].isalpha():
        text = "P" + text
    return text[:LABEL_MAX]


def moment_column(label: str) -> str:
    """The profile column of an index-less, unreactive polymer moment dummy
    (``rmgpy.chemkin.get_species_identifier`` of the MoM checkout): the label
    up to 16 characters, else ``<stem>-<md5[:6]><_muN>`` of exactly 16."""
    if 0 < len(label) <= 10:
        return label
    if re.search(r"[^A-Za-z0-9\-_,\(\)\*#.:\[\]]+", label):
        raise ProducerError(f"moment dummy label {label!r} has characters the profile writer replaces")
    if len(label) <= 16:
        return label
    digest = hashlib.md5(label.encode()).hexdigest()[:6]
    match = re.search(r"(_mu[0-9]+)$", label)
    suffix = match.group(1) if match else ""
    stem_len = 16 - len(suffix) - 7
    if stem_len < 1:
        return f"{label[:9]}-{digest}"
    return f"{label[:stem_len]}-{digest}{suffix}"


def species_column(header: Sequence[str], label: str) -> int:
    """Index of a species' profile column: ``label`` or ``label(<index>)``."""
    pattern = re.compile(re.escape(label) + r"(\(\d+\))?$")
    found = [k for k, h in enumerate(header) if pattern.fullmatch(h)]
    if len(found) != 1:
        raise ProducerError(f"species {label!r}: {len(found)} profile columns match (header {list(header)})")
    return found[0]


# --------------------------------------------------------------------------
# the deck's other statements
# --------------------------------------------------------------------------
def _module_calls(tree, name):
    return [stmt.value for stmt in tree.body if isinstance(stmt, ast.Expr) and isinstance(stmt.value, ast.Call)
            and isinstance(stmt.value.func, ast.Name) and stmt.value.func.id == name]


def _keyword_literal(call, key, path):
    for keyword in call.keywords:
        if keyword.arg == key:
            try:
                return ast.literal_eval(keyword.value)
            except (ValueError, TypeError, SyntaxError) as error:
                raise ProducerError(f"deck {path}: {key} (line {keyword.value.lineno}) is not a literal") from error
    return None


def deck_context(path) -> dict:
    """What the producer reads from the deck besides the declaration:
    ``database(thermoLibraries=...)``, the module-level ``species(...)``
    statements (source text, by label) and ``simpleReactor``'s pressure and
    initial mole fractions (the bath)."""
    path = Path(path)
    text = path.read_text()
    tree = ast.parse(text)
    out = {"thermo_libraries": list(DEFAULT_THERMO_LIBRARIES), "species": {}, "pressure": None, "bath": {}}
    databases = _module_calls(tree, "database")
    if len(databases) > 1:
        raise ProducerError(f"deck {path}: {len(databases)} database() calls")
    if databases:
        libraries = _keyword_literal(databases[0], "thermoLibraries", path)
        if libraries is not None:
            if not isinstance(libraries, (list, tuple)) or not all(isinstance(x, str) for x in libraries):
                raise ProducerError(f"deck {path}: thermoLibraries must be a list of names")
            out["thermo_libraries"] = list(libraries)
    for call in _module_calls(tree, "species"):
        label = _keyword_literal(call, "label", path)
        if not isinstance(label, str):
            raise ProducerError(f"deck {path}: species() at line {call.lineno} has no literal label")
        if label in out["species"]:
            raise ProducerError(f"deck {path}: species {label!r} declared twice")
        out["species"][label] = ast.get_source_segment(text, call)
    reactors = _module_calls(tree, "simpleReactor")
    if len(reactors) > 1:
        raise ProducerError(f"deck {path}: {len(reactors)} simpleReactor() calls; the producer reads one")
    if reactors:
        pressure = _keyword_literal(reactors[0], "pressure", path)
        if pressure is not None:
            if (not isinstance(pressure, tuple) or len(pressure) != 2
                    or not isinstance(pressure[0], (int, float)) or not isinstance(pressure[1], str)):
                raise ProducerError(f"deck {path}: simpleReactor pressure must be (value, 'units')")
            out["pressure"] = (float(pressure[0]), pressure[1])
        bath = _keyword_literal(reactors[0], "initialMoleFractions", path) or {}
        if not isinstance(bath, dict):
            raise ProducerError(f"deck {path}: initialMoleFractions must be a mapping")
        for label, x in bath.items():
            if label not in out["species"]:
                raise ProducerError(f"deck {path}: bath species {label!r} is not a deck species()")
            if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) or x < 0:
                raise ProducerError(f"deck {path}: bath mole fraction of {label!r} must be finite and >= 0")
        total = math.fsum(bath.values())
        out["bath"] = {k: float(v) / total for k, v in sorted(bath.items())} if total > 0 else {}
    return out


# --------------------------------------------------------------------------
# families and database
# --------------------------------------------------------------------------
def catalog_provenance(paths: Sequence) -> Tuple[List[str], List[Optional[str]]]:
    """(sorted union of the catalogs' provenance.families, their database paths)."""
    families, databases = set(), []
    for p in paths:
        try:
            doc = json.loads(Path(p).read_text())
        except (OSError, ValueError) as error:
            raise ProducerError(f"catalog {p}: {error}") from error
        provenance = doc.get("provenance") if isinstance(doc, dict) else None
        if not isinstance(provenance, dict):
            raise ProducerError(f"catalog {p}: no provenance")
        fams = provenance.get("families")
        if not isinstance(fams, list) or not fams or not all(isinstance(f, str) for f in fams):
            raise ProducerError(f"catalog {p}: provenance.families must be a non-empty list of family names")
        families.update(fams)
        databases.append(provenance.get("database_path"))
    return sorted(families), databases


def family_check(catalog_families: Sequence[str], families: Sequence[str], allow: bool) -> dict:
    """The family comparison; refused on a mismatch unless ``allow``."""
    given, catalog = sorted(set(families)), sorted(set(catalog_families))
    record = {"families": given, "catalog_families": catalog, "match": given == catalog,
              "only_in_families": sorted(set(given) - set(catalog)),
              "only_in_catalogs": sorted(set(catalog) - set(given))}
    if not record["match"]:
        if not allow:
            raise ProducerError(
                f"--families {given} differ from the catalogs' provenance families {catalog} "
                f"(only in --families: {record['only_in_families']}, only in the catalogs: "
                f"{record['only_in_catalogs']}); pass --allow-family-mismatch to run anyway, recorded")
        record["allowed_mismatch"] = True
    return record


def database_directory(given: Optional[str], catalog_databases: Sequence[Optional[str]]) -> Tuple[str, str]:
    """(RMG ``database.directory``, its source).  ``given`` (a database root or
    its ``input``) wins; else the catalogs' common ``database_path``."""
    if given is not None:
        root, source = Path(given), "--database"
    else:
        paths = sorted({p for p in catalog_databases if p})
        if len(paths) != 1 or any(p is None for p in catalog_databases):
            raise ProducerError(f"the catalogs name {paths or 'no'} database path(s); pass --database")
        root, source = Path(paths[0]), "catalog provenance.database_path"
    directory = root if (root / "kinetics").is_dir() else root / "input"
    if not (directory / "kinetics" / "families").is_dir():
        raise ProducerError(f"database {root}: no kinetics/families under {directory}")
    return str(directory), source


# --------------------------------------------------------------------------
# the MoM input
# --------------------------------------------------------------------------
def _species_statement(label: str, adjacency: Optional[str] = None, smiles: Optional[str] = None) -> str:
    if smiles is not None:
        return f"species(label={label!r}, reactive=True, structure=SMILES({smiles!r}))"
    return f"species(label={label!r}, reactive=True, structure=adjacencyList(\"\"\"{adjacency}\"\"\"))"


def declared_species(spec) -> Dict[str, Tuple[Optional[str], Optional[str]]]:
    """Every name the kMC observables report under ``species_mol`` ->
    (adjacency list, SMILES) of its structure."""
    identity = me.SpeciesIdentity(spec)
    out = {}
    graphs = {name: graph for name, graph, _d in spec.small_species}
    for name in mo.declared_species_names(spec, identity):
        if name in BUILTIN_STRUCTURES:
            out[name] = (None, BUILTIN_STRUCTURES[name])
        elif name == "monomer_product":
            out[name] = (spec.monomer_product.to_legacy_molecule({}).to_adjacency_list(), None)
        else:
            out[name] = (graphs[name].to_legacy_molecule({}).to_adjacency_list(), None)
    return out


def _parse_constraints(items: Sequence[str]) -> Dict[str, Any]:
    out = {}
    for item in items or ():
        if "=" not in item:
            raise ProducerError(f"--species-constraint {item!r}: expected KEY=VALUE")
        key, raw = item.split("=", 1)
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9]*", key):
            raise ProducerError(f"--species-constraint {item!r}: bad key")
        try:
            value = ast.literal_eval(raw)
        except (ValueError, SyntaxError) as error:
            raise ProducerError(f"--species-constraint {item!r}: value is not a literal") from error
        out[key] = value
    return out


def build_input(deck, context: dict, *, temperature: float, end: float, families: Sequence[str],
                constraints: Optional[Dict[str, Any]] = None) -> Tuple[str, dict]:
    """The MoM input text and the facts the conversion needs."""
    spec = deck.spec
    if deck.xs < 2:
        raise ProducerError(f"deck cutoff {deck.xs}: the MoM polymer() needs cutoff >= 2")
    label = mom_label(deck.label)
    repeat, (head, tail) = spec.legacy_molecules()
    for name, group in (("head", head), ("tail", tail)):
        if group.get_radical_count() != 1:
            raise ProducerError(f"{name} end group: the MoM polymer() needs a mono-radical end group "
                                f"(the port), got {group.get_radical_count()} radical electrons")
    a = spec.repeat_unit.composition.molar_mass_g_mol()
    e = spec.end_composition.molar_mass_g_mol()
    mu = [float(x) for x in spec.macroscopic_moments_mol]
    scale = spec.initialization().scale_factor
    species = declared_species(spec)
    clash = sorted((set(species) | {label}) & set(context["species"]))
    if clash:
        raise ProducerError(f"deck species {clash} collide with the declared species / pool label")
    mobile = {name: count * scale / constants.Na for name, count in spec.initial_mobile_counts if count}
    graphs = {name: graph for name, graph, _d in spec.small_species}
    mobile_kg = {name: amount * graphs[name].composition.molar_mass_kg_mol() for name, amount in mobile.items()}
    phase = [label] + (sorted(mobile) if spec.density_basis == "whole_initial_mixture" else [])
    initial = {label: mu[0]}
    initial.update({name: amount for name, amount in sorted(mobile.items())})
    initial.update({name: x * BATH_GAS_TOTAL_MOL for name, x in context["bath"].items()})
    pressure = context["pressure"] or DEFAULT_PRESSURE
    constraints = dict(constraints or {})
    if constraints:
        constraints.setdefault("allowed", ["input species"])   # the declared species stay admissible
    lines = [
        "# RMG method-of-moments input written by rmgpy.kmc.mom_producer from the kMC deck",
        f"# {deck.path} (sha256 {deck.sha256}); conventions: ../conventions.json",
        "database(",
        f"    thermoLibraries={list(context['thermo_libraries'])!r},",
        "    reactionLibraries=[],",
        "    seedMechanisms=[],",
        "    kineticsDepositories='default',",
        f"    kineticsFamilies={sorted(set(families))!r},",
        "    kineticsEstimator='rate rules',",
        ")",
    ]
    lines += [context["species"][k] for k in sorted(context["species"])]
    for name in sorted(species):
        adjacency, smiles = species[name]
        lines.append(_species_statement(name, adjacency, smiles))
    lines += [
        "polymer(",
        f"    label={label!r},",
        f"    monomer=\"\"\"{repeat.to_adjacency_list()}\"\"\",",
        f"    end_groups=[{head.to_smiles()!r}, {tail.to_smiles()!r}],",
        f"    cutoff={deck.xs},",
        f"    moments=[{mu[0]!r}, {mu[1]!r}, {mu[2]!r}],",
        f"    initial_mass={a * mu[1] / 1000.0!r},",
        ")",
        f"phase = polymer_phase(label={label + '_melt'!r}, species={phase!r}, solvent={label!r}, "
        f"density=({spec.density_kg_m3!r}, 'kg/m^3'))",
        "hybridPolymerReactor(",
        f"    temperature=({float(temperature)!r}, 'K'),",
        f"    pressure=({pressure[0]!r}, {pressure[1]!r}),",
        f"    initialMoles={initial!r},",
        "    polymerPhase=phase,",
        f"    terminationTime=({float(end)!r}, 's'),",
        "    constant_gas_volume=False,",
        ")",
        "model(toleranceMoveToCore=0.1, toleranceInterruptSimulation=0.1)",
        "simulator(atol=1e-16, rtol=1e-08)",
        "options(units='si', generateOutputHTML=False, generatePlots=False, saveSimulationProfiles=True, "
        "saveEdgeSpecies=False)",
    ]
    if constraints:
        args = ", ".join(f"{k}={v!r}" for k, v in sorted(constraints.items()))
        lines.append(f"generatedSpeciesConstraints({args})")
    facts = {"label": label, "a_g_mol": a, "e_g_mol": e, "moments_mol": mu, "species": sorted(species),
             "phase_species": phase, "initial_moles": initial, "pressure": list(pressure),
             "constraints": constraints, "families": sorted(set(families)),
             "volume_ratio_mom_over_kmc": (a * mu[1] / 1000.0 + math.fsum(
                 mobile_kg[n] for n in phase[1:])) / (spec.initial_mass_kg + math.fsum(mobile_kg.values()))}
    return "\n".join(lines) + "\n", facts


def conventions(deck, facts: dict, *, temperature: float, family: dict, database: Tuple[str, str]) -> dict:
    """``conventions.json``: every MoM setting or definition that is not the
    deck's own, named, with its value."""
    spec = deck.spec
    rows = [
        ("label", facts["label"], f"MoM pool label from {deck.label!r}: letters and digits, <= {LABEL_MAX} "
                                  "characters (so its moment columns are the label verbatim)"),
        ("legacy_translation", "LinearHomopolymerSpec.legacy_molecules",
         "repeat unit as an adjacency list with the head port *2 and the tail port *1 (radical at each "
         "port); end groups as mono-radical SMILES (radical at the port)"),
        ("initial_moments", facts["moments_mol"],
         "macroscopic_moments_mol of the declared ensemble (exact DP moments, mol); the reactor takes "
         "mu0 from initialMoles and mu1, mu2 from the moment ratios"),
        ("mass_basis", "MoM excludes end groups",
         "polymer(initial_mass = a mu1 / 1000 kg); MoM's own condensed mass is "
         "sum_p (M_u mu1 - defect mu0) without end groups"),
        ("capped_chains", facts["e_g_mol"],
         "Mn, Mw and polymer_condensed_mass_g assign every MoM chain a n + e (e = head + tail end-group "
         "molar mass, g/mol), as the kMC's Mn / Mw do; kMC's condensed mass is exact, so after a bridge "
         "scission they differ by e per chain minted without end groups"),
        ("moments_all_dp", True,
         "MoM pool moments cover every DP of the pool (no explicit ladder); kMC's are the tail n > xs"),
        ("all_pools_summed", True, "moments, chain count and mass sum every pool in polymer_pools.json, "
                                   "spawned pools included"),
        ("volume", facts["volume_ratio_mom_over_kmc"],
         "MoM V_poly = (a mu1 + phase species mass) / rho (end groups excluded), constant; value = "
         "MoM V_poly / kMC volume (polymer with end groups + mobile species, same density)"),
        ("phase_species", facts["phase_species"],
         "density_basis " + spec.density_basis + ": species whose mass enters V_poly"),
        ("declared_species", facts["species"],
         "every name kMC reports under species_mol is a reactive MoM input species (gas phase in MoM, "
         "melt in kMC); initial amount count x scale / N_A from initial_mobile_counts, else 0"),
        ("bath_gas", facts["initial_moles"], f"the deck's simpleReactor mole fractions x "
                                             f"{BATH_GAS_TOTAL_MOL} mol (MoM gas phase; kMC has none)"),
        ("pressure", facts["pressure"], "the deck's simpleReactor pressure, else 1 bar"),
        ("temperature_K", float(temperature), "isothermal"),
        ("families", facts["families"], "database(kineticsFamilies=...); reactionLibraries = "
                                        "seedMechanisms = []; kineticsEstimator 'rate rules'; "
                                        "kineticsDepositories 'default'; no pressureDependence"),
        ("family_check", family, "the catalogs' provenance.families against --families"),
        ("database", database[0], "RMG database.directory (rmgrc), from " + database[1]),
        ("kernels_off", ["k_scission=0", "k_unzip=0", "k_homolysis", "k_depropagation",
                         "radical_qssa_unzip", "side_group_homolysis", "explicit_dp=False"],
         "no MoM moment kernel; the deck's side_group_homolysis (" +
         ("none" if deck.side_group_homolysis is None else "present") + ") is not passed"),
        ("no_cage", True, "MoM rate coefficients are intrinsic: the deck's transport_model (MET, "
                          "Collins-Kimball phi) has no MoM counterpart"),
        ("model_settings", {"toleranceMoveToCore": 0.1, "toleranceInterruptSimulation": 0.1,
                            "atol": 1e-16, "rtol": 1e-08, "constant_gas_volume": False},
         "RMG model / simulator settings the deck does not set"),
        ("species_constraints", facts["constraints"], "generatedSpeciesConstraints as given "
                                                      "(--species-constraint; allowed=['input species'] "
                                                      "unless given, so the declared input species stay "
                                                      "admissible), else none"),
        ("profile", "rmg/solver/simulation_1_<N>.csv, highest N",
         "the final-core solver profile (rows at the integrator's steps, amounts in mol)"),
        ("time_interpolation", "linear in time",
         "each profile column interpolated between the two rows bracketing the requested time (a row "
         "time takes the row); derived quantities from the interpolated columns; outside the profile "
         "refused"),
        ("polymer_conversion_basis", "profile t = 0", "1 - mass(t) / mass at the profile's first row"),
    ]
    return {"schema": CONVENTIONS_SCHEMA,
            "conventions": [{"name": n, "value": v, "statement": s} for n, v, s in rows],
            "quantity_sources": dict(QUANTITY_SOURCES), "not_produced": dict(NOT_PRODUCED)}


# --------------------------------------------------------------------------
# the profile
# --------------------------------------------------------------------------
def find_profile(rmg_dir) -> Path:
    found = []
    for p in (Path(rmg_dir) / "solver").glob("simulation_1_*.csv"):
        m = re.fullmatch(r"simulation_1_(\d+)\.csv", p.name)
        if m:
            found.append((int(m.group(1)), p))
    if not found:
        raise ProducerError(f"no solver profile simulation_1_<N>.csv under {Path(rmg_dir) / 'solver'}")
    return max(found)[1]


def read_profile(path) -> Tuple[List[str], List[List[float]]]:
    with open(path, newline="") as handle:
        rows = list(csv.reader(handle))
    if not rows or rows[0][:2] != ["Time (s)", "Volume (m^3)"]:
        raise ProducerError(f"profile {path}: header must start with 'Time (s)', 'Volume (m^3)'")
    header, data = rows[0], []
    for k, row in enumerate(rows[1:], 2):
        if len(row) != len(header):
            raise ProducerError(f"profile {path}: row {k} has {len(row)} fields, header {len(header)}")
        try:
            values = [float(x) for x in row]
        except ValueError as error:
            raise ProducerError(f"profile {path}: row {k}: {error}") from error
        if not all(math.isfinite(v) for v in values):
            raise ProducerError(f"profile {path}: row {k} is not finite")
        data.append(values)
    if not data:
        raise ProducerError(f"profile {path}: no rows")
    times = [r[0] for r in data]
    if any(b < a for a, b in zip(times, times[1:])):
        raise ProducerError(f"profile {path}: time is not non-decreasing")
    return header, data


def read_pools(path) -> List[dict]:
    try:
        doc = json.loads(Path(path).read_text())
    except (OSError, ValueError) as error:
        raise ProducerError(f"pool table {path}: {error}") from error
    pools = doc.get("pools") if isinstance(doc, dict) else None
    if not isinstance(pools, list) or not pools:
        raise ProducerError(f"pool table {path}: no pools")
    for p in pools:
        if not isinstance(p, dict) or not isinstance(p.get("label"), str) or not _positive(p.get("monomer_mw_g_mol")):
            raise ProducerError(f"pool table {path}: each pool needs a label and monomer_mw_g_mol")
    return pools


def _positive(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x) and x > 0


def interpolate(data: List[List[float]], t: float) -> List[float]:
    """The profile columns at time ``t``: a row at ``t`` exactly, else linear in
    time between the last row before ``t`` and the first row after it."""
    times = [r[0] for r in data]
    if not times[0] <= t <= times[-1]:
        raise ProducerError(f"time {t!r} s is outside the profile [{times[0]!r}, {times[-1]!r}] s")
    at = [r for r in data if r[0] == t]
    if at:
        if any(r != at[0] for r in at[1:]):
            raise ProducerError(f"time {t!r} s: several profile rows with different values")
        return list(at[0])
    right = bisect.bisect_right(times, t)
    lo, hi = data[right - 1], data[right]
    w = (t - lo[0]) / (hi[0] - lo[0])
    return [x + (y - x) * w for x, y in zip(lo, hi)]


def convert(rmg_dir, times: Sequence[float], facts: dict) -> Tuple[List[dict], List[dict], dict]:
    """(the --mom entries, the diagnostics entries, provenance of the inputs read)."""
    rmg_dir = Path(rmg_dir)
    profile = find_profile(rmg_dir)
    pools_path = rmg_dir / "chemkin" / "polymer_pools.json"
    header, data = read_profile(profile)
    pools = read_pools(pools_path)
    a, e = facts["a_g_mol"], facts["e_g_mol"]
    by_label = {p["label"]: p for p in pools}
    if facts["label"] not in by_label:
        raise ProducerError(f"pool table has no pool {facts['label']!r} (pools {sorted(by_label)})")
    m_u = by_label[facts["label"]]["monomer_mw_g_mol"]
    if abs(m_u - a) > MASS_RELATIVE * a:
        raise ProducerError(f"MoM repeat-unit mass {m_u!r} differs from the deck's {a!r} g/mol")
    columns = {}
    for p in pools:
        cols = []
        for k in range(3):
            name = moment_column(f"{p['label']}_mu{k}")
            if name not in header:
                raise ProducerError(f"profile {profile}: no column {name!r} for pool {p['label']!r}")
            cols.append(header.index(name))
        columns[p["label"]] = cols
    used = {i for cols in columns.values() for i in cols}
    stray = [h for i, h in enumerate(header) if re.search(r"_mu[0-9]+$", h) and i not in used]
    if stray:
        raise ProducerError(f"profile {profile}: moment columns {stray} belong to no pool of {pools_path}")
    species = {name: species_column(header, name) for name in facts["species"]}

    def state(row):
        mu = [math.fsum(row[columns[p["label"]][k]] for p in pools) for k in range(3)]
        code = math.fsum(row[columns[p["label"]][1]] * p["monomer_mw_g_mol"]
                         - row[columns[p["label"]][0]] * float(p.get("chain_mass_defect_g_mol") or 0.0)
                         for p in pools)
        return mu, code

    mu_initial, code_initial = state(data[0])
    if data[0][0] != 0.0:
        raise ProducerError(f"profile {profile}: the first row is at t = {data[0][0]!r} s, not 0")
    mass_initial = code_initial + e * mu_initial[0]
    entries, diagnostics = [], []
    for t in times:
        row = interpolate(data, t)
        mu, code = state(row)
        mass = code + e * mu[0]
        obs = observables_from_moments(mu, a, e)
        quantities = {"mu0_mol": mu[0], "mu1_mol": mu[1], "mu2_mol": mu[2],
                      "Mn_g_mol": obs.mn_g_mol, "Mw_g_mol": obs.mw_g_mol,
                      "polymer_condensed_mass_g": mass,
                      "polymer_conversion": 1.0 - mass / mass_initial,
                      "polymer_derived_molecules_mol": mu[0],
                      "species_mol": {name: row[i] for name, i in species.items()}}
        check_quantity_keys(quantities, facts["species"])
        entries.append({"time_s": t, "quantities": quantities})
        diagnostics.append({"time_s": t, "volume_m3": row[1], "code_condensed_mass_g": code,
                            "native_Mn_g_mol": a * mu[1] / mu[0] if mu[0] else None,
                            "native_Mw_g_mol": a * mu[2] / mu[1] if mu[1] else None,
                            "pools": {p["label"]: [row[i] for i in columns[p["label"]]] for p in pools}})
    read = {"profile": str(profile), "profile_sha256": _sha256(profile.read_bytes()),
            "profile_rows": len(data), "profile_time_range_s": [data[0][0], data[-1][0]],
            "pools": str(pools_path), "pools_sha256": _sha256(pools_path.read_bytes()),
            "pool_labels": [p["label"] for p in pools]}
    return entries, diagnostics, read


# --------------------------------------------------------------------------
# running RMG
# --------------------------------------------------------------------------
BOOT = (
    "import os, runpy, sys\n"
    "checkout = os.path.realpath(sys.argv[1])\n"
    "import rmgpy\n"
    "where = os.path.realpath(rmgpy.__file__)\n"
    "open('rmgpy_path.txt', 'w').write(where + '\\n')\n"
    "if not where.startswith(checkout + os.sep):\n"
    "    sys.exit('rmgpy is imported from %s, not the MoM checkout %s' % (where, checkout))\n"
    "sys.argv = ['rmg.py', 'input.py']\n"
    "runpy.run_path(os.path.join(checkout, 'rmg.py'), run_name='__main__')\n")


def checkout_identity(path) -> dict:
    def git(*args):
        try:
            done = subprocess.run(["git", "-C", str(path), *args], capture_output=True, text=True, timeout=60)
        except (OSError, subprocess.SubprocessError):
            return None
        return done.stdout.strip() if done.returncode == 0 else None
    status = git("status", "--porcelain", "--untracked-files=no")
    return {"path": str(path), "sha": git("rev-parse", "HEAD"), "dirty": None if status is None else bool(status)}


def run_rmg(rmg_dir, checkout, python=None) -> int:
    rmg_dir, checkout = Path(rmg_dir), Path(checkout)
    if not (checkout / "rmg.py").is_file() or not (checkout / "rmgpy").is_dir():
        raise ProducerError(f"{checkout} is not an RMG-Py checkout (no rmg.py / rmgpy)")
    env = dict(os.environ, PYTHONPATH=str(checkout), PYTHONHASHSEED="0")
    with open(rmg_dir / "stdout.log", "ab") as out, open(rmg_dir / "stderr.log", "ab") as err:
        done = subprocess.run([python or sys.executable, "-c", BOOT, str(checkout)], cwd=str(rmg_dir), env=env,
                              stdout=out, stderr=err)
    if done.returncode != 0:
        tail = (rmg_dir / "stderr.log").read_text(errors="replace").splitlines()[-15:]
        raise ProducerError(f"RMG exited with {done.returncode} (see {rmg_dir / 'stderr.log'}):\n" + "\n".join(tail))
    return done.returncode


# --------------------------------------------------------------------------
# command line
# --------------------------------------------------------------------------
def produce(deck_path, *, checkout, temperature: float, end: float, snapshots: int, families: Sequence[str],
            catalogs: Sequence, out, mode: str = "write", allow_family_mismatch: bool = False,
            database: Optional[str] = None, species_constraints: Sequence[str] = ()) -> dict:
    """Write (and with ``mode`` "run" / "convert", run and convert); returns ``producer.json``."""
    if mode not in ("write", "run", "convert"):
        raise ProducerError(f"mode {mode!r}")
    out = Path(out)
    rmg_dir = out / "rmg"
    if mode == "convert":
        if not (rmg_dir / "input.py").is_file():
            raise ProducerError(f"--convert needs an existing run under {rmg_dir}")
    elif out.exists() and any(out.iterdir()):
        raise ProducerError(f"output directory {out} exists and is not empty")
    temperature = float(temperature)
    if not math.isfinite(temperature) or temperature <= 0:
        raise ProducerError(f"temperature {temperature!r} must be positive")
    times = kr.observation_times(end=end, snapshots=snapshots)
    if not families:
        raise ProducerError("--families needs at least one family")
    try:
        deck = kr.read_deck(deck_path)
    except kr.RunError as error:
        raise ProducerError(str(error)) from error
    context = deck_context(deck_path)
    catalog_families, catalog_databases = catalog_provenance(catalogs)
    family = family_check(catalog_families, families, allow_family_mismatch)
    db = database_directory(database, catalog_databases)
    text, facts = build_input(deck, context, temperature=temperature, end=times[-1], families=families,
                              constraints=_parse_constraints(species_constraints))
    conv = conventions(deck, facts, temperature=temperature, family=family, database=db)
    manifest = {"schema": SCHEMA, "status": "written", "mode": mode,
                "deck": {"path": deck.path, "sha256": deck.sha256, **deck.digests()},
                "catalogs": [{"path": str(p), "sha256": _sha256(Path(p).read_bytes())} for p in catalogs],
                "family_check": family, "temperature_K": temperature, "times_s": times,
                "mom_checkout": checkout_identity(checkout), "input_sha256": _sha256(text.encode()),
                "facts": facts}
    rmg_dir.mkdir(parents=True, exist_ok=True)
    if mode == "convert":
        if (rmg_dir / "input.py").read_text() != text:
            raise ProducerError(f"{rmg_dir / 'input.py'} is not the input these arguments write")
    else:
        (rmg_dir / "input.py").write_text(text)
        (rmg_dir / "rmgrc").write_text(f"database.directory : {db[0]}\n")
    (out / "conventions.json").write_text(json.dumps(conv, indent=1, sort_keys=True, allow_nan=False) + "\n")
    if mode == "run":
        run_rmg(rmg_dir, checkout)
        manifest["status"] = "ran"
    if mode in ("run", "convert"):
        entries, diagnostics, read = convert(rmg_dir, times, facts)
        (out / "mom.json").write_text(json.dumps(entries, indent=1, allow_nan=False) + "\n")
        (out / "mom_diagnostics.json").write_text(json.dumps(diagnostics, indent=1, allow_nan=False) + "\n")
        rmgpy_path = rmg_dir / "rmgpy_path.txt"
        manifest.update({"status": "converted", "read": read, "mom_file": str(out / "mom.json"),
                         "mom_sha256": _sha256((out / "mom.json").read_bytes()),
                         "rmgpy_imported": rmgpy_path.read_text().strip() if rmgpy_path.exists() else None})
    (out / "producer.json").write_text(json.dumps(manifest, indent=1, sort_keys=True, allow_nan=False) + "\n")
    return manifest


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m rmgpy.kmc.mom_producer", description=__doc__.split("\n\n")[0])
    parser.add_argument("deck")
    parser.add_argument("--mom-checkout", required=True)
    parser.add_argument("--temperature", type=float, required=True)
    parser.add_argument("--end", type=float, required=True)
    parser.add_argument("--snapshots", type=int, required=True)
    parser.add_argument("--families", nargs="+", required=True)
    parser.add_argument("--catalog", nargs="+", required=True)
    parser.add_argument("--out", required=True)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--run", action="store_const", dest="mode", const="run")
    group.add_argument("--write-only", action="store_const", dest="mode", const="write")
    group.add_argument("--convert", action="store_const", dest="mode", const="convert")
    parser.add_argument("--allow-family-mismatch", action="store_true")
    parser.add_argument("--database")
    parser.add_argument("--species-constraint", action="append", default=[])
    args = parser.parse_args(argv)
    try:
        manifest = produce(args.deck, checkout=args.mom_checkout, temperature=args.temperature, end=args.end,
                           snapshots=args.snapshots, families=args.families, catalogs=args.catalog,
                           out=args.out, mode=args.mode or "write",
                           allow_family_mismatch=args.allow_family_mismatch, database=args.database,
                           species_constraints=args.species_constraint)
    except (ProducerError, kr.RunError) as error:
        print(f"refused: {error}", file=sys.stderr)
        return 2
    print(f"{manifest['status']}: {args.out}" + (f" -> {manifest['mom_file']}" if "mom_file" in manifest else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
