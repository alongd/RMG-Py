"""Method-of-moments producer for the melt kMC runner's ``--mom`` file (I-169).

Command line::

    python -m rmgpy.kmc.mom_producer DECK --mom-checkout PATH --temperature T
        --end E --snapshots N --families F [F ...] --catalog FILE [FILE ...] --out DIR
        [--run | --write-only | --convert] [--allow-family-mismatch] [--allow-inactive]
        [--database DIR] [--species-constraint KEY=VALUE ...]

The same deck the kMC runner reads (``rmgpy.kmc.run.read_deck``: one
``polymer(..., linear_homopolymer_spec=...)``) is written as an RMG
method-of-moments (MoM) input for a polymer-mom checkout, whose ``polymer()``
takes legacy arguments only.  The translation, and every MoM setting the deck
does not determine, is a named convention in ``<out>/conventions.json``
(written in every mode).

**Chemistry scope.**  Only the RMG families given (``database(kineticsFamilies=
[...])``; no reaction libraries, no seed mechanisms) and none of the MoM moment
kernels (``k_scission``, ``k_unzip``, ``k_homolysis``, ``k_depropagation``,
``radical_qssa_unzip``, ``side_group_homolysis``, ``explicit_dp``).  The
families must equal the union of the kMC catalogs' ``provenance.families``;
otherwise the producer refuses unless ``--allow-family-mismatch``, and then the
mismatch is recorded in ``conventions.json`` and ``producer.json``.

**Kernel mode** (``--kernel-rates TABLE --kernel ID ... [--phi PHI]``, no
families, no catalogs): the chemistry is the kernel-rates table's homolysis and
end-recombination kernels (``kernel_input``) as polymer-mom's ``kernel_homolysis``
on the chain pool and ``kernel_recombination`` on an explicit end-radical pool
``<label>_re``, its two-participant transfer kernels (end + end and mid-chain
disproportionation, H transfer, reverse disproportionation) as
``kernel_transfer`` between those pools, an explicit terminal-alkene pool
``<label>_ut`` and explicit interior-radical and internal-alkene pools
``<label>_rm`` / ``<label>_ui``, and its intramolecular H transfer
(backbiting) as ``kernel_move`` (templates it cannot map listed in the
record; ``--kernel-far-law`` maps a law exact only from a longer length from
there), at the table's temperature; the run is active when its profile moves
(``kernel_census``).

**Modes.**  ``--write-only`` (default) writes ``<out>/rmg/input.py``,
``<out>/rmg/rmgrc`` (the kinetics database: ``--database`` or the catalogs'
common ``provenance.database_path``), ``conventions.json`` and
``producer.json``.  ``--run`` also runs RMG from the checkout in a subprocess
(``PYTHONPATH`` = the checkout, which must be the ``rmgpy`` imported; both
streams to ``<out>/rmg/stdout.log`` / ``stderr.log``) and converts.
``--convert`` converts an existing ``<out>/rmg`` run.

**Only kMC's definitions.**  A key of ``melt_observables.mom_reported_quantities``
is exported only when the MoM run gives that key's kMC definition
(``QUANTITY_SOURCES``); every other key is ``NOT_PRODUCED``, with the reason.
In particular kMC's ``mu*_mol`` / ``Mn_g_mol`` / ``Mw_g_mol`` are tail
quantities (chains with n > xs), which three moments of a MoM pool holding
every DP do not determine; the whole-ensemble keys ``whole_mu*_mol`` /
``whole_Mn_g_mol`` / ``whole_Mw_g_mol`` (``melt_observables.
WHOLE_ENSEMBLE_DEFINITION``) are on the mass DP (M - e) / a on both sides:
the pools' moments shifted to it plus every polymer-class MoM species at its
mass DP, withheld where MoM's pool moments leave that coordinate
(``pool_coordinate_breaks``, ``whole_moments``) or are not realizable within
MoM's own floors (``realizable``).  The condensed mass is MoM's own ledger booking
(``polymer.pyx`` ``physical_elemental_inventory``: per pool
``M_u mu1 - defect mu0 + e (mu0 + L)``, ``L = mu0(0) - mu0 - H`` with ``H = 0``
as no explicit-DP species exists) plus every other MoM species that kMC's
own classifier (``MeltObserver.category`` on the species' structure, the
backbone-origin count bracketed by element) puts in the polymer class, for the
mass and the molecule count alike.  The mass is exported only where the kMC
inventory (pools + every non-bath species) and the bath each close to
``MASS_CLOSURE_RELATIVE`` of their initial values, and the polymer rows only
where no present species is unclassifiable (``KmcClasses``).

**Refusals and flags.**  Conversion refuses: a first profile row not at
t = 0; conflicting rows at one time anywhere in the profile; a pool repeat
mass off the deck's; an initial inventory (summed moments, booked mass, every
species column and the bath) off the deck's by more than ``INITIAL_RELATIVE``;
a requested time outside the profile.  The pool-coupled reactions of the final
mechanism are counted per family (``census``); when no row of the requested
families is live (``ACTIVITY_RULE``) the run is refused (``status = "refused:inactive_chemistry"``, no ``--mom``
file) unless ``--allow-inactive``, and then flagged.  The flags are written to
``conventions.json``, ``producer.json`` and the ``--mom`` file's header (the
runner copies the header into ``run.json`` and the flag names onto every
comparison row): ``inactive_chemistry`` and ``kinetics_not_comparable`` (a declared
species' concentration per amount in MoM, at its residence, over kMC's melt at
t = 0 differs from 1 by more than the deck's ``observable_relative_tolerance``).

**Residence.**  ``polymer_phase(residence=...)`` carries where each species
lives (``species_residence``): the declared species, H and H2 at the
declaration's residence (melt), ``'*'`` (every other reactive species) melt,
the bath gas.  MoM condenses the melt species in its gas mask under kMC's
single convention (gas thermo, Kc at C0 = P0/RT, melt concentrations).

**Conversion.**  The final-core solver profile ``rmg/solver/simulation_1_<N>.csv``
(highest ``N``; columns are mol), the pool table ``rmg/chemkin/polymer_pools.json``
and the species dictionary ``rmg/chemkin/species_dictionary.txt`` (columns are
matched to structures, so RMG's formula-named identifiers resolve) give the
quantities at exactly the deck's observation times
(``rmgpy.kmc.run.observation_times``).  Each profile column is interpolated
linearly in time between the two profile rows bracketing the requested time
(a time on a row takes the row); derived quantities come from the interpolated
columns.  Output: ``<out>/mom.json`` (the ``--mom`` file, header form:
``{"header": {"flags", "conventions_sha256", "producer"}, "entries": [...]}``) and
``<out>/mom_diagnostics.json`` (MoM's own numbers at the same times:
whole-ensemble moments per pool and summed, its end-group-free code mass, its
native Mn / Mw, the mass-closure residual).
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
from rmgpy.kmc import melt_state as ms
from rmgpy.kmc import rmg_moves as rr
from rmgpy.kmc import run as kr
from rmgpy.kmc.event_record import same_species
from rmgpy.molecule import Molecule
from rmgpy import species_identity as _species_identity
from rmgpy.molecule.molecule import Bond
from rmgpy.polymer_spec import DEFAULT_RESIDENCE
from rmgpy.quantity import Quantity

SCHEMA = "kmc-mom-producer/1"
CONVENTIONS_SCHEMA = "kmc-mom-conventions/1"
BATH_GAS_TOTAL_MOL = 1.0
DEFAULT_PRESSURE = (1.0, "bar")
DEFAULT_THERMO_LIBRARIES = ("primaryThermoLibrary",)
LABEL_MAX = 10
MASS_RELATIVE = 1.0e-9           # pool repeat mass against the deck's
INITIAL_RELATIVE = 1.0e-9        # t = 0 moments and mass against the declared ensemble
MASS_CLOSURE_RELATIVE = 1.0e-6   # booked total mass against its t = 0 value
SIMULATOR_ATOL = 1.0e-16         # simulator(atol=...) the producer writes
SIMULATOR_RTOL = 1.0e-8          # simulator(rtol=...) the producer writes
# the profile moves (ACTIVITY_RULE) when a physical quantity changes by more than both floors: a relative
# one, ten times the integrator's relative tolerance, and an absolute one in grams, the integrator's relative
# tolerance on the deck's initial polymer mass (the error it allows on the condensed mass)
ACTIVITY_RELATIVE_FLOOR = 10.0 * SIMULATOR_RTOL
ACTIVITY_MASS_FLOOR_FRACTION = SIMULATOR_RTOL
MOM_SMALL_EPS = 1.0e-30          # rmgpy/solver/polymer.pyx:73 SMALL_EPS
MOM_EXHAUSTION_FLOOR_K = 100.0   # rmgpy/solver/polymer.pyx:153 EXHAUSTION_FLOOR_K
# MoM's per-moment accepted-state floor f = max(SMALL_EPS, EXHAUSTION_FLOOR_K atol) (polymer.pyx:4288)
MOM_MOMENT_FLOOR = max(MOM_SMALL_EPS, MOM_EXHAUSTION_FLOOR_K * SIMULATOR_ATOL)
STATUS_INACTIVE = "refused:inactive_chemistry"
BUILTIN_STRUCTURES = {"H": "[H]", "H2": "[H][H]"}
KERNEL_RATES_SCHEMAS = ("kmc-kernel-rates/2", "kmc-kernel-rates/3")
RADICAL_POOL_SUFFIX = "_re"      # the explicit end-radical pool of kernel mode: <label>_re
ALKENE_POOL_SUFFIX = "_ut"       # the explicit terminal-alkene pool of kernel mode: <label>_ut
INTERIOR_RADICAL_POOL_SUFFIX = "_rm"   # the explicit interior-radical pool of kernel mode: <label>_rm
INTERIOR_ALKENE_POOL_SUFFIX = "_ui"    # the explicit internal-alkene pool of kernel mode: <label>_ui

QUANTITY_SOURCES = {
    "polymer_condensed_mass_g": (
        "kMC's polymer-derived mass: MoM's ledger booking of its pools, sum over pools p of "
        "M_u,p mu1_p - defect_p mu0_p + e_p (mu0_p + L_p), L_p = mu0_p(0) - mu0_p (H_p = 0: no "
        "explicit-DP species), plus every other MoM species kMC's classifier puts in the polymer "
        "class (amount x molar mass of its structure); exported only where the kMC inventory closes "
        "and the bath closes (MASS_CLOSURE_RELATIVE) and no present species is unclassifiable"),
    "polymer_conversion": "1 - polymer_condensed_mass_g / (1000 initial_mass_kg); with the mass only",
    "polymer_derived_molecules_mol": (
        "kMC's count of polymer-derived molecules (tail mu0 + explicit ladder: every component of class "
        "polymer): sum over pools p of <p>_mu0 plus the amount of every other MoM species in the polymer "
        "class; not exported where a present species is unclassifiable"),
    "species_mol": ("profile column of each declared species, matched to its structure in the species "
                    "dictionary by event_record.same_species (resonance- and spin-aware)"),
    "whole_mu0_mol": ("sum over pools of <p>_mu0 plus the amount of every MoM species kMC's classifier "
                      "puts in the polymer class (melt_observables.WHOLE_ENSEMBLE_DEFINITION)"),
    "whole_mu1_mol": ("sum of mass DPs n_m = (M - e) / a (deck a, e): (MoM's booked pool chain mass - e "
                      "sum mu0_p) / a plus every polymer-class species' amount x n_m of its structure"),
    "whole_mu2_mol": ("sum of n_m^2: each pool's mu2 shifted to the mass DP by (e_p - defect_p - e) / a, plus "
                      "polymer-class species' amount x n_m^2; withheld where the end-group ledger books a "
                      "deficit (whole_moments)"),
    "whole_Mn_g_mol": "(a mu1 + e mu0) / mu0 of the whole moments (melt_observables.whole_averages)",
    "whole_Mw_g_mol": ("(a^2 mu2 + 2 a e mu1 + e^2 mu0) / (a mu1 + e mu0) of the whole moments, after MoM's "
                       "accepted-state floors (realizable)"),
}
_TAIL = ("kMC's {} is a tail quantity (chains with n > xs, xs = the deck cutoff >= 2); a MoM pool's three "
         "moments cover every DP and do not determine the tail")
NOT_PRODUCED = {
    "mu0_mol": _TAIL.format("mu0_mol"),
    "mu1_mol": _TAIL.format("mu1_mol"),
    "mu2_mol": _TAIL.format("mu2_mol"),
    "Mn_g_mol": _TAIL.format("Mn_g_mol (of the tail moments)"),
    "Mw_g_mol": _TAIL.format("Mw_g_mol (of the tail moments)"),
    "end_inventory_mol": (
        "kMC's is the physical end-group inventory of the tail; MoM's ledger nets minted against "
        "removed chains over the whole pool and is written to no file"),
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


def _rel(x: float, y: float) -> float:
    scale = max(abs(x), abs(y))
    return abs(x - y) / scale if scale else 0.0


# --------------------------------------------------------------------------
# the runner's quantity keys
# --------------------------------------------------------------------------
def runner_quantity_keys() -> frozenset:
    """The top-level keys ``melt_observables.mom_reported_quantities`` returns
    (evaluated on a minimal observation, so the runner's own definition is
    the reference)."""
    obs = {"tier1": {"polymer_derived": {"mass_g": 1.0, "mom_definition": {"amount_mol": 1.0}},
                     "species": {}},
           "tier2": {"mass_dp_moments_mol": [1.0, 2.0, 4.0]},
           "tier3": {"tail": {"dp_moments_mol": [1.0, 2.0, 4.0]}, "explicit": {"dp_amount_mol": {}}},
           "mom_mass_law": {"end_inventory_mol": 1.0},
           "repeat_unit_molar_mass_g_mol": 1.0, "end_group_molar_mass_g_mol": 1.0,
           "initial_polymer_mass_g": 1.0}
    return frozenset(mo.mom_reported_quantities(obs))


def check_quantity_keys(quantities: dict, declared: Sequence[str]) -> None:
    """Refuse a quantity the runner does not report, one declared not
    produced, or a species the deck does not declare."""
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


def read_species_dictionary(path) -> Dict[str, Molecule]:
    """``chemkin/species_dictionary.txt``: identifier -> molecule (blocks of an
    identifier line and an adjacency list, separated by blank lines)."""
    try:
        text = Path(path).read_text()
    except OSError as error:
        raise ProducerError(f"species dictionary {path}: {error}") from error
    out = {}
    for block in re.split(r"\n\s*\n", text.strip()):
        lines = block.strip().splitlines()
        if not lines:
            continue
        name, adjacency = lines[0].strip(), "\n".join(lines[1:])
        try:
            out[name] = Molecule().from_adjacency_list(adjacency)
        except Exception as error:     # the adjacency-list parser raises several types
            raise ProducerError(f"species dictionary {path}: {name!r}: {error}") from error
    return out


def resolve_species(header: Sequence[str], dictionary: Dict[str, Molecule],
                    structures: Dict[str, Molecule]) -> Dict[str, int]:
    """Declared name -> profile column, by structure: the column whose
    dictionary molecule is the declared species (``event_record.same_species``:
    one resonance set at one multiplicity, on copies; RMG names a column
    ``label(index)`` or, past 16 characters, ``formula(index)``)."""
    out = {}
    for name, molecule in structures.items():
        found = [k for k, h in enumerate(header) if h in dictionary and same_species(dictionary[h], molecule)]
        if len(found) != 1:
            raise ProducerError(f"species {name!r}: {len(found)} profile columns hold its structure "
                                f"({[header[k] for k in found]})")
        out[name] = found[0]
    return out


# --------------------------------------------------------------------------
# kMC's classes for MoM species
# --------------------------------------------------------------------------
def _element(decl) -> Tuple[str, Optional[int]]:
    return decl.element, (None if decl.isotope in (-1, None) else int(decl.isotope))


def ancestry_elements(spec) -> Tuple[frozenset, frozenset]:
    """(elements every atom of which is backbone, elements some atom of which
    is backbone), as (element, isotope) over the atoms of the declared melt:
    the repeat unit's head-tail path (``melt_state._backbone_keys``) and an
    initial monomer's mapped atoms (``melt_state.monomer_backbone_keys``) are
    backbone; every other repeat atom, every end-group atom, every other
    initial mobile atom and every hydrogen are not."""
    backbone, other = set(), set()

    def add(graph, flagged):
        for decl in graph.atoms:
            (backbone if decl.key in flagged else other).add(_element(decl))
            if decl.hydrogens:
                other.add(("H", None))

    add(spec.repeat_unit, ms._backbone_keys(spec.repeat_unit))
    add(spec.head_end_group, frozenset())
    add(spec.tail_end_group, frozenset())
    graphs = {name: graph for name, graph, _d in spec.small_species}
    for name, count in spec.initial_mobile_counts:
        if count:
            add(graphs[name], ms.monomer_backbone_keys(spec, graphs[name]) or frozenset())
    return frozenset(backbone - other), frozenset(backbone)


class KmcClasses:
    """kMC's class of a MoM species: ``melt_observables.MeltObserver.category``
    (declared species by identity, else polymer when n_c = floor(B_c / b) >= 1,
    else light) on the species' structure.  MoM keeps no ancestry, so B_c is
    bracketed: the atoms of elements that are only ever backbone (lower) and
    of elements that are ever backbone (upper); the class is kMC's when both
    give the same, else "undetermined" (also for a species the kMC melt
    cannot represent, such as RMG's inert Ar / He / Ne: it blocks the polymer
    rows only where present)."""

    def __init__(self, spec, xs: int):
        try:
            self.observer = mo.MeltObserver(spec, xs=xs)
        except (mo.ObservablesError, me.DeclarationError, ms.MeltStateError) as error:
            raise ProducerError(f"kMC observer of the deck: {error}") from error
        self.only, self.possible = ancestry_elements(spec)

    def classify(self, molecule) -> str:
        mol = molecule.copy(deep=True)
        for i, atom in enumerate(mol.atoms):
            atom.id = i
        rows = rr.molecule_rows(mol)
        nodes = {i: ms.Node(el, iso, charge, radical, lone, dict(edges))
                 for i, (el, iso, charge, radical, lone, edges) in rows.items()}
        classes = set()
        for allowed in (self.only, self.possible):
            origin = {i: ms.Origin("mobile", 0, "mom", str(i), (row[0], row[1]) in allowed)
                      for i, row in rows.items()}
            try:
                state = ms.MeltState(nodes, origin)
            except ms.MeltStateError:      # an atom the kMC melt cannot hold (RMG's inert Ar, He, Ne)
                return "undetermined"
            ids = state.components.ids()
            if len(ids) != 1:
                raise ProducerError(f"species {molecule.to_smiles()} is {len(ids)} molecules")
            classes.add(self.observer.category(state, ids[0]))
        return classes.pop() if len(classes) == 1 else "undetermined"


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


def species_structures(facts: dict) -> Dict[str, Molecule]:
    out = {}
    for name, (adjacency, smiles) in facts["species_structures"].items():
        out[name] = Molecule(smiles=smiles) if smiles is not None else Molecule().from_adjacency_list(adjacency)
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


def phase_kinetics(spec, species: Sequence[str], initial: Dict[str, float], gas: Sequence[str],
                   temperature: float, pressure: Tuple[float, str], *, v_poly: float,
                   residence: Dict[str, str]) -> dict:
    """Per declared species, its concentration per amount in MoM over kMC's
    melt at t = 0: ``V_kMC / V_poly`` for a melt-resident species (MoM: c =
    n / V_poly, its polymer volume), ``V_kMC / V_gas`` for a gas one (c = n /
    V_gas, the ideal-gas volume of the gas moles, 1 m^3 when there are none);
    kMC: c = n / V, V the melt volume of the declaration."""
    init = spec.initialization()
    v_kmc = init.volume_m3 * init.scale_factor
    n_gas = math.fsum(initial[k] for k in gas if k in initial)
    p_pa = Quantity(pressure[0], pressure[1]).value_si
    v_gas = constants.R * float(temperature) * n_gas / p_pa if n_gas > 0 else 1.0
    tolerance = spec.observable_relative_tolerance
    ratios = {name: v_kmc / (v_poly if residence.get(name, residence["*"]) == "melt" else v_gas)
              for name in species}
    return {"ratio_by_species": ratios, "V_kmc_m3": v_kmc, "V_poly_m3": v_poly, "V_gas_m3": v_gas,
            "gas_mol": n_gas, "tolerance": tolerance,
            "not_comparable": sorted(n for n, r in ratios.items() if abs(r - 1.0) > tolerance)}


def species_residence(spec, species: Sequence[str], bath: Sequence[str]) -> Dict[str, str]:
    """The ``polymer_phase(residence=...)`` MoM takes: every declared species
    (H and H2 included) at its declared residence, ``'*'`` (every other
    reactive species, the undeclared lights and generated ones) at the
    declaration's default, the bath gas (kMC has no gas phase; the bath only
    sets MoM's gas volume)."""
    out = {"*": DEFAULT_RESIDENCE}
    out.update({name: spec.residence(name) for name in species})
    out.update({name: "gas" for name in bath})
    return dict(sorted(out.items()))


def _uniform_cut(law: dict, a: float, radical: Optional[str], where: str) -> Tuple[int, str]:
    """(min_fragment e, the radical class) of a homolysis law whose exact
    outcomes are, on every n it was measured on, two radicals of one class,
    k and n - k for every k in [e, n - e], a sites per position; refused
    otherwise."""
    fragments = [p for _n, outs in law["exact"].items() for _m, prods in outs for p in prods]
    if not fragments:
        raise ProducerError(f"{where}: the law cuts no chain up to n = {law['exact_to_n']}")
    e = min(int(p[1]) for p in fragments)
    classes = {p[0] for p in fragments}
    if len(classes) != 1 or (radical is not None and classes != {radical}):
        raise ProducerError(f"{where}: fragments of classes {sorted(classes)}, not one end-radical pool "
                            f"{radical or ''}")
    for n_text, outs in sorted(law["exact"].items(), key=lambda x: int(x[0])):
        n, got, want = int(n_text), {}, {}
        for m, prods in outs:
            pair = tuple(sorted(int(p[1]) for p in prods))
            if len(prods) != 2 or sum(pair) != n:
                raise ProducerError(f"{where}: at n = {n} an outcome {prods} is not two fragments of n")
            got[pair] = got.get(pair, 0) + m
        for k in range(e, n - e + 1):
            pair = tuple(sorted((k, n - k)))
            want[pair] = want.get(pair, 0) + a
        if got != want:
            raise ProducerError(f"{where}: at n = {n} the outcomes {got} are not the uniform cut over "
                                f"[{e}, {n - e}] with {a!r} sites per position ({want})")
    return e, classes.pop()


class _Unmapped(Exception):
    """A template (group) of a pair-transfer kernel with no exact moment form: listed, never dropped."""


SITE_LAW_RELATIVE = 1.0e-9       # a kernel's moment law against the table's site sums on the checked lengths
SITE_LAW_MARGIN = 4              # lengths checked past every exact law and chain edge (the far rule's own)


def _label_of(families: dict):
    """The table's class label of (family, n), None below the family's min_n
    (kernel_reference.KernelModel.label_of)."""
    def label_of(family, n):
        f = families.get(family)
        if f is None or n < int(f["min_n"]):
            return None
        return f.get("labels", {}).get(str(n), f["label_far"])
    return label_of


def _sites_at(law: dict, n: int, nbb: int, label_of, family_of: Dict[str, str]) -> Dict[str, float]:
    """{product family: sites} of one participant law on a molecule of n
    periods, as the kernel reference reads it (kernel_reference.KernelModel.
    _outcomes): the exact law up to exact_to_n, the far rule beyond, an
    interior law's far outcomes wherever their products exist; none where n
    is not the law's class.  Every outcome must be one product of the same
    length (_Unmapped otherwise)."""
    if label_of(law["family"], n) != law["class"]:
        return {}
    out: Dict[str, float] = {}

    def add(m, family):
        out[family] = out.get(family, 0.0) + float(m)

    def relative(m, prods, interior):
        if len(prods) != 1 or prods[0][0] == "join" or list(prods[0][1:]) != [1, 0, 0]:
            raise _Unmapped(f"its outcome {prods} is not one product of the same length")
        if label_of(prods[0][0], n) is None:
            if interior:
                return                    # no position places it: no site (the reference's own reading)
            raise _Unmapped(f"its product {prods[0][0]!r} is no state at {n} periods")
        add(m, prods[0][0])

    anchor = law["anchor"]
    if anchor in ("chain", "end") and n <= int(law["exact_to_n"]):
        for m, prods in law["exact"].get(str(n), []):
            if len(prods) != 1 or prods[0][0] == "join" or int(prods[0][1]) != n:
                raise _Unmapped(f"at {n} periods an outcome {prods} is not one product of the same length")
            if prods[0][0] not in family_of:
                raise _Unmapped(f"product class {prods[0][0]!r} is no family of the table")
            add(m, family_of[prods[0][0]])
        return out
    far = law["far"]
    if anchor == "chain":
        L, edge = n * nbb, int(far["edge"])
        for k in range(L):
            rows = (far["first"][str(k)] if k < edge else far["last"][str(L - 1 - k)] if L - 1 - k < edge
                    else far["bulk"])
            for m, prods in rows:
                relative(m, prods, False)
    else:
        for m, prods in far["outcomes"]:
            relative(m, prods, anchor == "interior")
    return out


def _law_top(laws: Sequence[dict], floor: int, nbb: int) -> int:
    """The longest length a site-law check reads: past every exact law and
    chain edge by SITE_LAW_MARGIN periods, where only far rules act."""
    top = floor
    for law in laws:
        if law["anchor"] in ("chain", "end"):
            top = max(top, int(law["exact_to_n"]))
        if law["anchor"] == "chain":
            top = max(top, -(-2 * int(law["far"]["edge"]) // nbb))
    return top + SITE_LAW_MARGIN


def _bilinear(F: Dict[Tuple[int, int], float], floor: int, top: int) -> Tuple[Optional[int], Optional[list], float]:
    """(n0, [c00, c10, c01, c11], worst) with F(n, m) = c00 + c10 n + c01 m
    + c11 n m on every n, m in [n0, top], n0 the smallest from ``floor``
    (None if none), and worst the largest relative difference on the lengths
    below n0."""
    scale = max(abs(v) for v in F.values())
    tol = SITE_LAW_RELATIVE * scale
    for n0 in range(floor, top - 1):
        f = lambda a, b: F[(a, b)]
        raw11 = f(n0 + 1, n0 + 1) - f(n0, n0 + 1) - f(n0 + 1, n0) + f(n0, n0)
        raw10 = f(n0 + 1, n0) - f(n0, n0) - raw11 * n0
        raw01 = f(n0, n0 + 1) - f(n0, n0) - raw11 * n0
        # a slope worth less than the fit tolerance over the checked lengths is rounding in the table's sums:
        # the law without it is tried first, so a constant stays a constant
        for c10, c01, c11 in ((0.0 if abs(raw10) * top <= tol else raw10, 0.0 if abs(raw01) * top <= tol else raw01,
                               0.0 if abs(raw11) * top * top <= tol else raw11), (raw10, raw01, raw11)):
            c00 = f(n0, n0) - c10 * n0 - c01 * n0 - c11 * n0 * n0
            law = lambda a, b: c00 + c10 * a + c01 * b + c11 * a * b
            if all(abs(law(a, b) - v) <= tol for (a, b), v in F.items() if a >= n0 and b >= n0):
                worst = max([abs(law(a, b) - v) / max(abs(law(a, b)), abs(v)) for (a, b), v in F.items()
                             if (a < n0 or b < n0) and max(abs(law(a, b)), abs(v)) > 0.0], default=0.0)
                return n0, [c00, c10, c01, c11], worst
    return None, None, 0.0


def _affine(r: Dict[int, float], floor: int, top: int) -> Tuple[Optional[int], Optional[list], float]:
    """(n0, [B, A], worst) with r(n) = A n + B on every n in [n0, top]
    (_bilinear's rule in one length)."""
    scale = max(abs(v) for v in r.values())
    tol = SITE_LAW_RELATIVE * scale
    for n0 in range(floor, top):
        raw = r[n0 + 1] - r[n0]
        for A in ((0.0 if abs(raw) * top <= tol else raw), raw):
            B = r[n0] - A * n0
            if all(abs(A * n + B - v) <= tol for n, v in r.items() if n >= n0):
                worst = max([abs(A * n + B - v) / max(abs(A * n + B), abs(v)) for n, v in r.items()
                             if n < n0 and max(abs(A * n + B), abs(v)) > 0.0], default=0.0)
                return n0, [B, A], worst
    return None, None, 0.0


def _exact_down(n0: int, good) -> int:
    """The shortest length L <= n0 with good(n) on every n in [L, n0): the
    scan down from a law's fit start, stopping at its first miscount."""
    length = n0
    while length > 1 and good(length - 1):
        length -= 1
    return length


def _site_or_none(law, n, nbb, label_of, family_of) -> Optional[Dict[str, float]]:
    try:
        return _sites_at(law, n, nbb, label_of, family_of)
    except _Unmapped:
        return None


def _pair_exact_from(members, laws, fam0, fam1, c, n0, top, F, nbb, label_of, family_of) -> List[int]:
    """Per participant, the shortest length from which a pair term's law
    c00 + c10 n + c01 m + c11 n m is the table's sum_T (w/divisor) s0(n)
    s1(m) for every length of the other participant from the fit start n0
    to the checked top (the scan goes below the reachable floor, down to the
    first miscount or a length the table has no state for)."""
    tol = SITE_LAW_RELATIVE * max(abs(v) for v in F.values())
    cache = {}

    def site(name, slot, n):
        if (name, slot, n) not in cache:
            cache[(name, slot, n)] = _site_or_none(laws[name][slot], n, nbb, label_of, family_of)
        return cache[(name, slot, n)]

    def ok(n, m):
        terms = []
        for name, k, _s in members:
            s0, s1 = site(name, 0, n), site(name, 1, m)
            if s0 is None or s1 is None:
                return False
            terms.append(k * s0.get(fam0, 0.0) * s1.get(fam1, 0.0))
        return abs(c[0] + c[1] * n + c[2] * m + c[3] * n * m - math.fsum(terms)) <= tol
    return [_exact_down(n0, lambda n: all(ok(n, m) for m in range(n0, top + 1))),
            _exact_down(n0, lambda m: all(ok(n, m) for n in range(n0, top + 1)))]


def _move_exact_at(members, laws, fam, c, n, r, nbb, label_of, family_of) -> bool:
    """A move term's law c0 + c1 n is the table's sum_T w s(n) at n."""
    terms = []
    for name, w, _s in members:
        s = _site_or_none(laws[name], n, nbb, label_of, family_of)
        if s is None:
            return False
        terms.append(w * s.get(fam, 0.0))
    return abs(c[0] + c[1] * n - math.fsum(terms)) <= SITE_LAW_RELATIVE * max(abs(v) for v in r.values())


def kernel_input(table: dict, spec, kernels: Sequence[str], *, temperature: float, phi: float, label: str,
                 radical_label: str, volume_ratio: float, alkene_label: Optional[str] = None,
                 interior_radical_label: Optional[str] = None, interior_alkene_label: Optional[str] = None,
                 far_law: bool = False) -> dict:
    """polymer-mom's ``kernel_homolysis``, ``kernel_recombination``,
    ``kernel_transfer`` and ``kernel_move`` terms of the requested kernels of a kernel-rates
    table, and the record of the checks that make their moment laws exact on
    this deck.

    Homolysis (``homolysis``, a chain class): per template the table's weight
    w, its sites a*n + b (n in the table's period, ``periods_per_unit`` per
    declared unit) and ``min_fragment`` e from its exact law, which must be
    the uniform cut (``_uniform_cut``) with b = a*(1 - 2e); phi (``phi_thinned``
    templates) the given cage-escape factor; the fragments into the end-radical
    pool.  The sites are the count on every chain the H + C system holds only
    if no initial chain is shorter than 2e - 1 periods (refused otherwise;
    recombination products are at least 2e).
    Recombination (``recombination``, two radicals of the fragments' class
    joining into one chain, gamma 0): k_pair = (w / divisor) phi^thinned
    V_poly / V, kMC's per-pair coefficient in MoM's concentration n / V_poly;
    sites_per_chain the table's per_feature, whose profile (zero sites at
    short radicals) must lie below e, the shortest radical homolysis makes.
    Pair transfer (any other two-participant kernel, e.g. end + end
    disproportionation and its reverse, the mid-chain ones, H transfer):
    each participant changes family, length unchanged, into a pool kernel
    mode carries (``chain`` the chain pool, ``radical_end`` the end-radical
    pool, ``alkene_terminal`` the terminal-alkene pool ``alkene_label``,
    ``radical_interior`` / ``alkene_internal`` the interior-radical and
    internal-alkene pools ``interior_radical_label`` /
    ``interior_alkene_label``); the term is per kernel and product pair, not
    per template, because a table splits one site set over templates by
    context (a short donor's own template).  Its sites are the kernel
    reference's reading of each law (``_sites_at``), and its templates map
    together when F(n, m) = sum_T (w/divisor)_T s_T0(n) s_T1(m) is one
    bilinear law c00 + c10 n + c01 m + c11 n m on every pair of lengths from
    the reachable floor (the shortest length the requested kernels and the
    deck make: homolysis' min_fragment, else the shortest initial chain) up
    to past every exact law (``_law_top``); a constant law is the I-221 term
    (k_pair = K V_poly / V), any other carries ``sites`` (k_pair the largest
    coefficient times V_poly / V); phi on phi_thinned templates (one value
    per group).  Pool move (a one-participant intramolecular H transfer,
    backbiting): per product family sum_T w_T s_T(n) = A n + B the same way,
    ``kernel_move`` k*(a n + b) with k the larger of |A|, |B|; one that keeps
    its pool (a mid-chain radical's own shift) changes no moment and is
    listed under ``identity``.  A law exact only from a longer length n0 is
    listed under ``unmapped`` unless ``far_law``, which maps it from n0 and
    lists the shorter lengths and their largest relative difference under
    ``short_lengths``.  A template (group) whose law fails, or whose products
    have no pool, is listed under ``unmapped`` with its reason; one the table
    gives no law (``law_inert``, unmodelled) under ``inert`` (rate 0 by the
    table's own law, as in the kernel reference); never dropped silently; a
    kernel whose participant classes have no pool is refused.
    Every term declares its validity length ``valid_from_periods``, the
    shortest chain its site law counts correctly (homolysis 2e - 1, where the
    cut count a*n + b turns negative below; recombination one past its
    zero-site profile; a pair transfer or move the length its law is exact
    from), and the record's ``validity`` lists them with the participant
    families, for the kernel reference to report the fraction of each
    family's chains below (``kernel_reference.law_domain_fractions``).
    A depository override of a mapped template (RMG's own rate for the
    reaction at exact molecules, which kMC fires at its own weight) has no
    moment form, since moments cannot single out sizes: MoM's kernel carries
    the template law at every size, and the record lists each override under
    ``depository_departures`` (kernel, template, override, depository entry,
    reactant and product states, multiplicity, ``w_kmc`` the override's
    weight, ``w_mom`` the template's, both before divisor and phi, which apply
    alike, and ``kmc_over_mom``).
    Refused: another schema, temperature or chemistry digest, kernels not in
    the table or of another type, inert or unmodelled homolysis or
    recombination laws."""
    from rmgpy.kmc import kernel_rates as K
    from rmgpy.kmc.melt_engine import declared_initialization
    if table.get("schema") not in KERNEL_RATES_SCHEMAS:
        raise ProducerError(f"not a kernel-rates table of {KERNEL_RATES_SCHEMAS} (got {table.get('schema')!r})")
    if float(table["temperature_K"]) != float(temperature):
        raise ProducerError(f"the table is at temperature {table['temperature_K']!r} K, the run at {temperature!r} "
                            "K: kernel mode is isothermal at the table's temperature")
    bound = (table.get("catalog") or {}).get("chemistry_sha256")
    if bound != kr.chemistry_sha256(spec):
        raise ProducerError(f"the table's chemistry digest {bound!r} is not the deck's ({kr.chemistry_sha256(spec)})")
    if isinstance(phi, bool) or not isinstance(phi, (int, float)) or not 0.0 <= phi <= 1.0:
        raise ProducerError(f"phi {phi!r} must be in [0, 1]")
    unit = table["moment_unit"]
    ppu, nbb = int(unit["periods_per_declared_unit"]), int(unit["backbone_atoms_per_period"])
    a_unit = spec.repeat_unit.composition.molar_mass_g_mol()
    if abs(float(unit["molar_mass_g_mol"]) * ppu - a_unit) > MASS_RELATIVE * a_unit:
        raise ProducerError(f"{ppu} periods of {unit['molar_mass_g_mol']!r} g/mol are not the declared unit's "
                            f"{a_unit!r} g/mol")
    unknown = sorted(set(kernels) - set(table["kernels"]))
    if unknown or not kernels:
        raise ProducerError(f"kernels {unknown or 'none'} not in the table (kernels {sorted(table['kernels'])})")
    rows = {r["template"]: r for r in table["templates"]}
    n_min = min(dp for dp, _c in declared_initialization(spec).dp_counts) * ppu
    by_type = {"homolysis": [], "recombination": [], "transfer": [], "move": []}
    for kid in sorted(set(kernels)):
        kernel = table["kernels"][kid]
        kind = kernel["type"] if kernel["type"] in ("homolysis", "recombination") else (
            "transfer" if len(kernel["classes"]) == 2 else
            "move" if kernel["type"] == "intra_h_transfer" and len(kernel["classes"]) == 1 else None)
        if kind is None:
            raise ProducerError(f"{kid}: no MoM kernel carries type {kernel['type']!r} (homolysis, "
                                "recombination, two-participant pair transfer and intramolecular H transfer "
                                "only)")
        by_type[kind].append((kid, kernel))
    radical = None
    for kid, kernel in by_type["recombination"]:
        if len(set(kernel["classes"])) != 1 or not kernel["classes"][0].startswith("radical_end:"):
            raise ProducerError(f"{kid}: no MoM kernel joins classes {kernel['classes']} (one end-radical class)")
        radical = kernel["classes"][0] if radical is None else radical
        if kernel["classes"][0] != radical:
            raise ProducerError(f"{kid}: a second end-radical class {kernel['classes'][0]!r} ({radical!r})")
    homolysis, recombination, checks, e_min, validity, laws_out = [], [], {}, None, [], []
    for kid, kernel in by_type["homolysis"]:
        if kernel["classes"] != ["chain"]:
            raise ProducerError(f"{kid}: no MoM kernel cuts class {kernel['classes']} (the chain class only)")
        for term in kernel["templates"]:
            name, where = term["template"], f"{kid} template {term['template']}"
            if term.get("law") is None or term.get("law_inert"):
                raise ProducerError(f"{where}: an inert or unmodelled law has no moment form")
            sites = rows[name]["participants"][0]["chain"]
            a, b = float(sites["a"]), float(sites["b"])
            e, radical = _uniform_cut(term["law"], a, radical, where)
            if b != a * (1 - 2 * e):
                raise ProducerError(f"{where}: sites a*n + b with b = {b!r}, not b = a*(1 - 2*min_fragment) = "
                                    f"{a * (1 - 2 * e)!r} of its cut law")
            if n_min < 2 * e - 1:
                raise ProducerError(f"{where}: initial chains of {n_min} periods are shorter than {2 * e - 1}, "
                                    f"where a*n + b = {a * n_min + b!r} is not their site count")
            homolysis.append({"template": name, "w": float(rows[name]["weight"]),
                              "phi": float(phi) if rows[name]["phi_thinned"] else 1.0, "a": a, "b": b,
                              "min_fragment": e, "periods_per_unit": ppu, "daughter": radical_label,
                              "valid_from_periods": 2 * e - 1})
            # a n + b = a (n - 2e + 1) counts the cuts leaving both fragments at least e long from 2e - 1
            # (zero there); shorter chains have no cut while the law goes negative
            validity.append({"kernel": kid, "template": name, "families": ["chain"],
                             "valid_from_periods": 2 * e - 1, "fit_from_periods": 2 * e - 1,
                             "participants": [{"slot": 0, "family": "chain", "exact_from_periods": 2 * e - 1}]})
            w_h = float(rows[name]["weight"]) * (float(phi) if rows[name]["phi_thinned"] else 1.0)
            laws_out.append({"kernel": kid, "template": name, "templates": [name], "kind": "uni",
                             "families": ["chain"], "products": None, "coefficients": [w_h * b, w_h * a]})
            checks[name] = {"min_fragment": e, "exact_to_n": term["law"]["exact_to_n"],
                            "exact_from_n": sites.get("exact_from_n"), "initial_min_periods": n_min}
            e_min = e if e_min is None else min(e_min, e)
    for kid, kernel in by_type["recombination"]:
        for term in kernel["templates"]:
            name, where = term["template"], f"{kid} template {term['template']}"
            if term.get("law") is None or term.get("law_inert"):
                raise ProducerError(f"{where}: an inert or unmodelled law has no moment form")
            row = rows[name]
            sites = [p[radical] for p in row["participants"]]
            c = float(sites[0]["per_feature"])
            zero = set()
            for s in sites:
                if float(s["per_feature"]) != c:
                    raise ProducerError(f"{where}: participants with {s['per_feature']!r} and {c!r} sites")
                for key, counts in s.get("profile", {}).items():
                    m = re.fullmatch(r"fragment_backbone_atoms=(\d+)", key)
                    if not m or list(counts) != [0] or int(m.group(1)) % nbb:
                        raise ProducerError(f"{where}: a profile entry {key}: {counts} has no moment form")
                    zero.add(int(m.group(1)) // nbb)
            if e_min is not None and zero and max(zero) >= e_min:
                raise ProducerError(f"{where}: its profile gives no site at {sorted(zero)} periods, where "
                                    f"homolysis makes radicals (from {e_min})")
            join = [c, [["join", "chain", 0]]]
            for slot in ("S0", "S1"):
                law = term["law"][slot]
                if law["far"]["outcomes"] != [join] or any(
                        outs != ([] if int(n) in zero else [join]) for n, outs in law["exact"].items()):
                    raise ProducerError(f"{where}: its {slot} law is not one chain of n0 + n1 periods per site")
            k = float(row["weight"]) / float(row["divisor"]) * (float(phi) if row["phi_thinned"] else 1.0)
            valid = max(zero) + 1 if zero else 1
            recombination.append({"template": name, "k_pair": k * float(volume_ratio), "sites_per_chain": c,
                                  "periods_per_unit": ppu, "product": label, "valid_from_periods": valid})
            validity.append({"kernel": kid, "template": name, "families": ["radical_end"],
                             "valid_from_periods": valid, "fit_from_periods": valid,
                             "participants": [{"slot": p, "family": "radical_end", "exact_from_periods": valid}
                                              for p in (0, 1)]})
            laws_out.append({"kernel": kid, "template": name, "templates": [name], "kind": "pair",
                             "families": ["radical_end", "radical_end"], "products": None,
                             "coefficients": [k * c * c, 0.0, 0.0, 0.0]})
            checks[name] = {"w_over_divisor": k, "zero_site_periods": sorted(zero)}
    floor = n_min if e_min is None else min(n_min, e_min)
    families = table["families"]
    family_of = {lab: fam for fam, f in families.items() for lab in [f["label_far"], *f.get("labels", {}).values()]}
    label_of = _label_of(families)
    pools = {"chain": label, "radical_end": radical_label, "alkene_terminal": alkene_label,
             "radical_interior": interior_radical_label, "alkene_internal": interior_alkene_label}
    transfer, move, unmapped, inert, identity, short = {}, {}, [], [], [], []

    def owner_families(kid, kernel):
        nonlocal radical
        owners = []
        for cls in kernel["classes"]:
            fam = family_of.get(cls)
            if fam not in pools or pools[fam] is None or cls != families[fam]["label_far"]:
                raise ProducerError(f"{kid}: no MoM pool carries its participant class {cls!r} (pools: "
                                    f"{sorted(f for f, p in pools.items() if p)} at their far classes)")
            if fam == "radical_end":
                radical = cls if radical is None else radical
                if cls != radical:
                    raise ProducerError(f"{kid}: a second end-radical class {cls!r} ({radical!r})")
            owners.append(fam)
        return owners

    def usable(kid, term):
        """False (listed) for a template the table gives no law."""
        name = term["template"]
        if term.get("law") is None or term.get("law_inert"):
            inert.append({"kernel": kid, "templates": [name],
                          "reason": term.get("law_inert") or term.get("law_refused") or "no product law"})
            return False
        return True

    def short_lengths(kid, names, n0, worst, what):
        """The far-law record (and the strict refusal) of a law exact only from n0 > floor."""
        if far_law:
            short.append({"kernel": kid, "templates": names, "from_periods": floor, "law_from_periods": n0,
                          "worst_relative": worst})
            return True
        unmapped.append({"kernel": kid, "templates": names, "reason": (
            f"{what} on lengths from {floor} periods (it is only from {n0} periods; its difference there up to "
            f"{worst:.3g} relative)")})
        return False

    for kid, kernel in by_type["transfer"]:
        owners = owner_families(kid, kernel)
        laws = {}
        for term in kernel["templates"]:
            if usable(kid, term):
                laws[term["template"]] = (term["law"]["S0"], term["law"]["S1"])
        top = _law_top([law for pair in laws.values() for law in pair], floor, nbb)
        groups = {}
        for name, (law0, law1) in laws.items():
            try:
                sites = [{n: _sites_at(law, n, nbb, label_of, family_of) for n in range(floor, top + 1)}
                         for law in (law0, law1)]
                fams = [sorted({f for s in slot.values() for f in s}) for slot in sites]
                for slot_fams, owner in zip(fams, owners):
                    for fam in slot_fams:
                        if pools.get(fam) is None:
                            raise _Unmapped(f"products of family {fam!r}: no MoM pool carries it")
                        if fam == owner:
                            raise _Unmapped(f"a participant stays in its family {fam!r}")
            except _Unmapped as why:
                unmapped.append({"kernel": kid, "templates": [name], "reason": str(why)})
                continue
            row = rows[name]
            k = float(row["weight"]) / float(row["divisor"])
            for fam0 in fams[0] or [None]:
                for fam1 in fams[1] or [None]:
                    groups.setdefault((fam0, fam1, bool(row["phi_thinned"])), []).append((name, k, sites))
        for (fam0, fam1, thinned), members in sorted(groups.items(), key=lambda x: [str(v) for v in x[0]]):
            names = sorted({n for n, _k, _s in members})
            F = {(n, m): math.fsum(k * s[0][n].get(fam0, 0.0) * s[1][m].get(fam1, 0.0) for _n, k, s in members)
                 for n in range(floor, top + 1) for m in range(floor, top + 1)}
            if fam0 is None or fam1 is None or max(F.values()) <= 0.0:
                inert.append({"kernel": kid, "templates": names,
                              "reason": f"no site on lengths from {floor} periods (inert on this deck)"})
                continue
            n0, c, worst = _bilinear(F, floor, top)
            if n0 is None:
                lo, hi = min(F, key=F.get), max(F, key=F.get)
                unmapped.append({"kernel": kid, "templates": names, "reason": (
                    f"sum_T (w/divisor) s0(n) s1(m) is not one constant, nor bilinear in n and m, on lengths from "
                    f"{floor} periods: {F[lo]!r} at (n, m) = {lo}, {F[hi]!r} at {hi}")})
                continue
            if n0 > floor and not short_lengths(
                    kid, names, n0, worst, "sum_T (w/divisor) s0(n) s1(m) is not one constant, nor bilinear in n "
                                           "and m,"):
                continue
            joined = "+".join(names)
            entry = {"template": joined, "phi": float(phi) if thinned else 1.0, "partner": pools[owners[1]],
                     "product": pools[fam0], "partner_product": pools[fam1], "valid_from_periods": n0}
            exact = _pair_exact_from(members, laws, fam0, fam1, c, n0, top, F, nbb, label_of, family_of)
            validity.append({"kernel": kid, "template": joined, "families": sorted(set(owners)),
                             "valid_from_periods": n0, "fit_from_periods": n0,
                             "participants": [{"slot": p, "family": owners[p], "exact_from_periods": exact[p]}
                                              for p in (0, 1)]})
            laws_out.append({"kernel": kid, "template": joined, "templates": names, "kind": "pair",
                             "families": list(owners), "products": [fam0, fam1],
                             "coefficients": [entry["phi"] * x for x in c]})
            check = {"kernel": kid, "templates": names, "reachable_from_periods": floor, "law_from_periods": n0,
                     "checked_to_periods": top}
            if not any(c[1:]):
                entry["k_pair"] = c[0] * float(volume_ratio)
                check["w_over_divisor_sites"] = c[0]
            else:
                scale = max(abs(x) for x in c)
                entry.update({"k_pair": scale * float(volume_ratio), "sites": [x / scale for x in c],
                              "periods_per_unit": ppu})
                check["sites_law"] = c
            transfer.setdefault(pools[owners[0]], []).append(
                {key: entry[key] for key in ("template", "k_pair", "phi", "partner", "product", "partner_product",
                                             "sites", "periods_per_unit", "valid_from_periods") if key in entry})
            checks[joined] = check
    for kid, kernel in by_type["move"]:
        (owner,) = owner_families(kid, kernel)
        laws = {term["template"]: term["law"] for term in kernel["templates"] if usable(kid, term)}
        top = _law_top(list(laws.values()), floor, nbb)
        groups = {}
        for name, law in laws.items():
            try:
                sites = {n: _sites_at(law, n, nbb, label_of, family_of) for n in range(floor, top + 1)}
                fams = sorted({f for s in sites.values() for f in s})
                for fam in fams:
                    if pools.get(fam) is None:
                        raise _Unmapped(f"products of family {fam!r}: no MoM pool carries it")
            except _Unmapped as why:
                unmapped.append({"kernel": kid, "templates": [name], "reason": str(why)})
                continue
            row = rows[name]
            w = float(row["weight"]) * (float(phi) if row["phi_thinned"] else 1.0)
            for fam in fams:
                groups.setdefault(fam, []).append((name, w, sites))
        for fam, members in sorted(groups.items()):
            names = sorted({n for n, _w, _s in members})
            if fam == owner:
                identity.append({"kernel": kid, "templates": names,
                                 "reason": f"its products stay in the {fam!r} pool, lengths unchanged"})
                continue
            r = {n: math.fsum(w * s[n].get(fam, 0.0) for _n, w, s in members) for n in range(floor, top + 1)}
            if max(r.values()) <= 0.0:
                inert.append({"kernel": kid, "templates": names,
                              "reason": f"no site on lengths from {floor} periods (inert on this deck)"})
                continue
            n0, c, worst = _affine(r, floor, top)
            if n0 is None:
                unmapped.append({"kernel": kid, "templates": names, "reason": (
                    f"sum_T w s(n) is not affine in n on lengths from {floor} periods")})
                continue
            if n0 > floor and not short_lengths(kid, names, n0, worst, "sum_T w s(n) is not affine in n"):
                continue
            joined = "+".join(names)
            scale = max(abs(x) for x in c)
            move.setdefault(pools[owner], []).append({"template": joined, "k": scale, "a": c[1] / scale,
                                                      "b": c[0] / scale, "periods_per_unit": ppu,
                                                      "product": pools[fam], "valid_from_periods": n0})
            exact = _exact_down(n0, lambda n: _move_exact_at(members, laws, fam, c, n, r, nbb, label_of, family_of))
            validity.append({"kernel": kid, "template": joined, "families": [owner], "valid_from_periods": n0,
                             "fit_from_periods": n0,
                             "participants": [{"slot": 0, "family": owner, "exact_from_periods": exact}]})
            laws_out.append({"kernel": kid, "template": joined, "templates": names, "kind": "uni",
                             "families": [owner], "products": [fam], "coefficients": [c[0], c[1]]})
            checks[joined] = {"kernel": kid, "templates": names, "reachable_from_periods": floor,
                              "law_from_periods": n0, "checked_to_periods": top, "rate_law": [c[0], c[1]]}
    # the overrides of every mapped template: kMC fires them at their weight, MoM's kernel at the template's
    departures = []
    for kid, name in dict.fromkeys((law["kernel"], n) for law in laws_out for n in law["templates"]):
        (term,) = [t for t in table["kernels"][kid]["templates"] if t["template"] == name]
        w_mom = float(rows[name]["weight"])
        for o in term.get("overrides") or ():
            departures.append({"kernel": kid, "template": name, "override": o["name"], "entry": o.get("entry"),
                               "reactants": o["reactants"], "products": o["products"],
                               "multiplicity": o["multiplicity"], "w_kmc": o["weight"], "w_mom": w_mom,
                               "kmc_over_mom": o["weight"] / w_mom})
    record = {"kernels": sorted(set(kernels)), "excluded_kernels": sorted(set(table["kernels"]) - set(kernels)),
              "temperature_K": float(table["temperature_K"]), "phi": float(phi), "moment_unit": unit,
              "volume_ratio_mom_over_kmc": float(volume_ratio), "radical_class": radical, "checks": checks,
              "far_law": bool(far_law), "unmapped": unmapped, "inert": inert, "identity": identity,
              "short_lengths": short, "validity": validity, "laws": laws_out,
              "depository_departures": departures, "chemistry_sha256": bound,
              "table_templates_sha256": K.template_rows_sha256(table, [n for law in laws_out
                                                                       for n in law["templates"]])}
    return {"homolysis": homolysis, "recombination": recombination, "transfer": transfer, "move": move,
            "record": record}


def build_input(deck, context: dict, *, temperature: float, end: float, families: Sequence[str],
                constraints: Optional[Dict[str, Any]] = None,
                kernels: Optional[dict] = None) -> Tuple[str, dict]:
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
    radical_label = label + RADICAL_POOL_SUFFIX if kernels is not None else None
    alkene_label = label + ALKENE_POOL_SUFFIX if kernels is not None else None
    mid_label = label + INTERIOR_RADICAL_POOL_SUFFIX if kernels is not None else None
    internal_label = label + INTERIOR_ALKENE_POOL_SUFFIX if kernels is not None else None
    pools = [label] + [p for p in (radical_label, alkene_label, mid_label, internal_label) if p]
    generated = {f"{label}_melt"} | {p + s for p in pools for s in ("", "_mu0", "_mu1", "_mu2")}
    clash = sorted((set(species) & generated) | ((set(species) | generated) & set(context["species"])))
    if clash:
        raise ProducerError(f"species labels {clash} collide with the MoM pool label {label!r}, its moment "
                            "or phase labels, or another species of the deck")
    mobile = {name: count * scale / constants.Na for name, count in spec.initial_mobile_counts if count}
    graphs = {name: graph for name, graph, _d in spec.small_species}
    mobile_kg = {name: amount * graphs[name].composition.molar_mass_kg_mol() for name, amount in mobile.items()}
    phase_mobile = sorted(mobile) if spec.density_basis == "whole_initial_mixture" else []
    phase = pools + phase_mobile
    initial = {label: mu[0]}
    initial.update({name: amount for name, amount in sorted(mobile.items())})
    initial.update({name: x * BATH_GAS_TOTAL_MOL for name, x in context["bath"].items()})
    pressure = context["pressure"] or DEFAULT_PRESSURE
    residence = species_residence(spec, sorted(species), sorted(context["bath"]))
    v_poly = (a * mu[1] / 1000.0 + math.fsum(mobile_kg[n] for n in phase_mobile)) / spec.density_kg_m3
    volume_ratio = (a * mu[1] / 1000.0 + math.fsum(mobile_kg[n] for n in phase_mobile)) / (
        spec.initial_mass_kg + math.fsum(mobile_kg.values()))
    terms = None
    if kernels is not None:
        terms = kernel_input(kernels["table"], spec, kernels["ids"], temperature=temperature, phi=kernels["phi"],
                             label=label, radical_label=radical_label, alkene_label=alkene_label,
                             interior_radical_label=mid_label, interior_alkene_label=internal_label,
                             volume_ratio=volume_ratio, far_law=kernels.get("far_law", False))
        terms["record"].update(kernels["source"])
        # the terminal-alkene and mid-chain pools exist only when a transfer or move term moves chains into
        # or out of them
        touched = {p for owner, ts in terms["transfer"].items() for t in ts
                   for p in (owner, t["partner"], t["product"], t["partner_product"])}
        touched |= {p for owner, ts in terms["move"].items() for t in ts for p in (owner, t["product"])}
        for name in (alkene_label, mid_label, internal_label):
            if name not in touched:
                pools.remove(name)
                phase.remove(name)
        alkene_label, mid_label, internal_label = (p if p in touched else None
                                                   for p in (alkene_label, mid_label, internal_label))
        # Each pool's chain end mass, from its kernel-rates family (a chain of n periods weighs n periods
        # plus the family's offset): MoM's polymer() takes only mono-radical end groups, so every pool is
        # declared with the deck's pair, and the conversion reads the pool's true end mass here.
        families = kernels["table"]["families"]
        pool_family = {pool: fam for fam, pool in (("chain", label), ("radical_end", radical_label),
                                                   ("alkene_terminal", alkene_label),
                                                   ("radical_interior", mid_label),
                                                   ("alkene_internal", internal_label)) if pool is not None}
        ends = {pool: float(families[fam]["molar_mass_offset_g_mol"]) for pool, fam in pool_family.items()}
        if abs(ends[label] - e) > MASS_RELATIVE * e:
            raise ProducerError(f"the table's chain family carries {ends[label]!r} g/mol of end groups, the deck "
                                f"{e!r} g/mol")
        terms["record"]["pool_end_mass_g_mol"] = ends
        # kMC's identity of each family's states (the table's probes): a pool chain that kMC classes as a
        # declared species, by its length in periods
        terms["record"]["pool_family"] = pool_family
        terms["record"]["declared_species_by_pool"] = {
            pool: sorted([name, int(key.rsplit("@", 1)[1])]
                         for key, name in (families[fam].get("declared_species") or {}).items())
            for pool, fam in pool_family.items()}

    def transfer_line(pool):
        return ([f"    kernel_transfer={terms['transfer'][pool]!r},"] if terms and pool in terms["transfer"] else []) \
            + ([f"    kernel_move={terms['move'][pool]!r},"] if terms and pool in terms["move"] else [])

    def born_empty(pool, declaration):
        return [
            "polymer(",
            f"    label={pool!r},",
            f"    monomer=\"\"\"{repeat.to_adjacency_list()}\"\"\",",
            f"    end_groups=[{head.to_smiles()!r}, {tail.to_smiles()!r}],",
            f"    cutoff={deck.xs},",
            "    moments=[0.0, 0.0, 0.0],",
            "    initial_mass=0.0,",
            f"    {declaration},",
        ]
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
        f"    kineticsFamilies={'none' if terms is not None else sorted(set(families))!r},",
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
    ]
    if terms is not None and terms["homolysis"]:
        lines.append(f"    kernel_homolysis={terms['homolysis']!r},")
    lines += transfer_line(label) + [")"]
    if terms is not None:
        lines += [
            "polymer(",
            f"    label={radical_label!r},",
            f"    monomer=\"\"\"{repeat.to_adjacency_list()}\"\"\",",
            f"    end_groups=[{head.to_smiles()!r}, {tail.to_smiles()!r}],",
            f"    cutoff={deck.xs},",
            "    moments=[0.0, 0.0, 0.0],",
            "    initial_mass=0.0,",
            "    end_radical_site='primary',",
        ] + ([f"    kernel_recombination={terms['recombination']!r},"] if terms["recombination"] else []) + \
            transfer_line(radical_label) + [")"]
    for pool, declaration in ((alkene_label, "end_alkene_site='primary'"),
                              (mid_label, "interior_feature='radical'"),
                              (internal_label, "interior_feature='alkene'")):
        if pool is not None:
            lines += born_empty(pool, declaration) + transfer_line(pool) + [")"]
    lines += [
        f"phase = polymer_phase(label={label + '_melt'!r}, species={phase!r}, solvent={label!r}, "
        f"density=({spec.density_kg_m3!r}, 'kg/m^3'), residence={residence!r})",
        "hybridPolymerReactor(",
        f"    temperature=({float(temperature)!r}, 'K'),",
        f"    pressure=({pressure[0]!r}, {pressure[1]!r}),",
        f"    initialMoles={initial!r},",
        "    polymerPhase=phase,",
        f"    terminationTime=({float(end)!r}, 's'),",
        "    constant_gas_volume=False,",
        ")",
        "model(toleranceMoveToCore=0.1, toleranceInterruptSimulation=0.1)",
        f"simulator(atol={SIMULATOR_ATOL!r}, rtol={SIMULATOR_RTOL!r})",
        "options(units='si', generateOutputHTML=False, generatePlots=False, saveSimulationProfiles=True, "
        "saveEdgeSpecies=False)",
    ]
    if constraints:
        args = ", ".join(f"{k}={v!r}" for k, v in sorted(constraints.items()))
        lines.append(f"generatedSpeciesConstraints({args})")
    facts = {"label": label, "a_g_mol": a, "e_g_mol": e, "moments_mol": mu,
             "initial_mass_g": spec.initial_mass_kg * 1000.0, "species": sorted(species),
             "species_structures": {k: list(v) for k, v in sorted(species.items())},
             "end_groups": [head.to_smiles(), tail.to_smiles()],
             "phase_species": phase, "initial_moles": initial, "bath": sorted(context["bath"]),
             "pressure": list(pressure),
             "constraints": constraints, "families": sorted(set(families)),
             "volume_ratio_mom_over_kmc": volume_ratio,
             "residence": residence,
             "phase_kinetics": phase_kinetics(spec, sorted(species), initial, sorted(context["bath"]), temperature,
                                              pressure, v_poly=v_poly, residence=residence)}
    if terms is not None:
        facts.update({"kernels": terms, "radical_pool": radical_label, "alkene_pool": alkene_label,
                      "interior_radical_pool": mid_label, "interior_alkene_pool": internal_label})
    return "\n".join(lines) + "\n", facts


def conventions(deck, facts: dict, *, temperature: float, family: dict, database: Tuple[str, str],
                flags: dict) -> dict:
    """``conventions.json``: every MoM setting or definition that is not the
    deck's own, named, with its value, and the run's flags."""
    spec = deck.spec
    rows = [
        ("label", facts["label"], f"MoM pool label from {deck.label!r}: letters and digits, <= {LABEL_MAX} "
                                  "characters (so its moment columns are the label verbatim); a declared "
                                  "species label equal to it, its moment or phase labels is refused"),
        ("legacy_translation", "LinearHomopolymerSpec.legacy_molecules",
         "repeat unit as an adjacency list with the head port *2 and the tail port *1 (radical at each "
         "port); end groups as mono-radical SMILES (radical at the port)"),
        ("initial_moments", facts["moments_mol"],
         "macroscopic_moments_mol of the declared ensemble (exact DP moments, mol); the reactor takes "
         "mu0 from initialMoles and mu1, mu2 from the moment ratios; polymer(initial_mass = a mu1 / 1000 "
         "kg), MoM's end-group-free basis"),
        ("initial_inventory", INITIAL_RELATIVE,
         "at the profile's t = 0 row the summed pool moments, the booked mass, every declared species "
         "(its initial amount, else 0), the bath (its initial amount) and every other species column (0) "
         "must equal the deck's to this relative tolerance, else refused"),
        ("tail_keys", sorted(k for k in NOT_PRODUCED if k.startswith(("mu", "M"))),
         "not produced: kMC's are tail quantities (n > xs); MoM's pool moments cover every DP "
         "(whole-ensemble moments are in mom_diagnostics.json, not compared)"),
        ("mass_booking", "MoM physical_elemental_inventory ledger",
         "condensed mass = sum_p M_u,p mu1_p - defect_p mu0_p + e_p (mu0_p + L_p), L_p = mu0_p(0) - mu0_p "
         "(H_p = 0, no explicit-DP species); e_p from the pool's end-group SMILES"),
        ("mass_closure", MASS_CLOSURE_RELATIVE,
         "two closures, each within this relative tolerance of its own t = 0 value: the kMC inventory "
         "(booked pools + every non-bath, non-bookkeeping species by its structure's molar mass) and the "
         "bath separately; where either fails, the mass and the conversion are not exported (residuals "
         "in mom_diagnostics.json)"),
        ("whole_ensemble", "mass DP n_m = (M - e) / a",
         mo.WHOLE_ENSEMBLE_DEFINITION + ".  MoM side: whole_moments (pools shifted by (e_p - defect_p - e) / a, "
         "mu1 from the booked chain mass, mu2 withheld on a nonzero end-group ledger), realizable (MoM's "
         "accepted-state floors, max(1e-30, 100 atol)), pool_coordinate_breaks (mechanisms that leave the "
         "coordinate withhold every whole key after t = 0); mu2 under MoM's closures (mu3 log-Lagrange, "
         "gamma) is MoM's model, compared as such"),
        ("pool_chain_classes", "kMC classes before polymer",
         "a pool chain kMC would class as a declared species (pool_species_candidates: chain mass at a whole DP "
         "equal to the species' mass, structure-checked for pools without defect) is booked as that species "
         "only where the pool has no defect, its moments pin that DP and the end-group ledger is zero; a "
         "defect pool keeps no structure, so a candidate in it is never booked by mass alone; any other "
         "non-empty candidate pool withholds the polymer rows, at every row including t = 0; a species row "
         "is exported only when every pool that may hold the species is booked or provably empty (all "
         "moments within f), never as a partial count; a pool receiving volatile ejections of positive "
         "moment units (eject_units corrected for the pools' defects, polymer.pyx 4110-4113) may hold chains "
         "below one unit's backbone (kMC class light): its polymer rows are withheld after t = 0; pool "
         "moments below -f, units without chains (mu0 = 0) or a variance below MoM's census tolerance are "
         "corruption (pool_state) and withhold mass, molecules and the whole keys"),
        ("kmc_classes", "melt_observables.MeltObserver.category",
         "every MoM species outside the pools and the bath is put in kMC's class (declared species, "
         "polymer-derived when floor(B / b) >= 1, else light) by kMC's own classifier on its structure; "
         "B, the backbone-origin atom count MoM does not keep, is bracketed by the elements only ever "
         "backbone and those ever backbone in the declared melt; a polymer-class species adds its mass and "
         "amount to the polymer rows; where a species whose bracket gives two classes is present, the "
         "polymer rows are not exported"),
        ("volume", facts["volume_ratio_mom_over_kmc"],
         "MoM V_poly = (a mu1 + phase species mass) / rho (end groups excluded), constant; value = "
         "MoM V_poly / kMC volume (polymer with end groups + mobile species, same density)"),
        ("phase_species", facts["phase_species"],
         "density_basis " + spec.density_basis + ": species whose mass enters V_poly"),
        ("declared_species", facts["species"],
         "every name kMC reports under species_mol is a reactive MoM input species (at its residence in "
         "MoM, melt in kMC); initial amount count x scale / N_A from initial_mobile_counts, else 0"),
        ("residence", facts["residence"],
         "polymer_phase(residence=...): each declared species (H, H2 included) at the declaration's "
         "residence (LinearHomopolymerSpec.residence: melt), '*' (every other reactive species, generated "
         "ones included) melt, the bath gas; MoM condenses melt species in its gas mask (n / V_poly, gas "
         "thermo, Kc at C0 = P0/RT: kMC's single convention) and, the convention declared, passes wholly "
         "melt rows through its reference-state tripwire"),
        ("phase_kinetics", facts["phase_kinetics"],
         "concentration per amount of each declared species, MoM (melt: n / V_poly; gas: n / V_gas) over "
         "kMC melt (n / V), at t = 0; kinetics_not_comparable when a ratio differs from 1 by more than the "
         "deck's observable_relative_tolerance"),
        ("bath_gas", facts["initial_moles"], f"the deck's simpleReactor mole fractions x "
                                             f"{BATH_GAS_TOTAL_MOL} mol (MoM gas phase; kMC has none): "
                                             "outside every polymer row and closed on its own"),
        ("pressure", facts["pressure"], "the deck's simpleReactor pressure, else 1 bar"),
        ("temperature_K", float(temperature), "isothermal"),
        ("families", facts["families"], "database(kineticsFamilies=...); reactionLibraries = "
                                        "seedMechanisms = []; kineticsEstimator 'rate rules'; "
                                        "kineticsDepositories 'default'; no pressureDependence"),
        ("family_check", family, "the catalogs' provenance.families against --families"),
        ("inactive_chemistry", STATUS_INACTIVE,
         "pool-coupled reactions of the final mechanism counted per family (polymer_pools.json "
         "reactions[].refused; family from the serialized id of a Cantera-null row, else from the "
         "chem_annotated.inp entries with the same species, a fall-off collider (+M) or (+<name>) dropped; "
         "entries without one counted by equation as unattributed requested-family rows); " + ACTIVITY_RULE + "; no live row: refused, unless "
         "--allow-inactive (then flagged)"),
        ("database", database[0], "RMG database.directory (rmgrc), from " + database[1]),
        ("kernels_off", ["k_scission=0", "k_unzip=0", "k_homolysis", "k_depropagation",
                         "radical_qssa_unzip", "side_group_homolysis", "explicit_dp=False"],
         "no MoM moment kernel; the deck's side_group_homolysis (" +
         ("none" if deck.side_group_homolysis is None else "present") + ") is not passed"),
        ("no_cage", True, "MoM rate coefficients are intrinsic: the deck's transport_model (MET, "
                          "Collins-Kimball phi) has no MoM counterpart"),
        ("model_settings", {"toleranceMoveToCore": 0.1, "toleranceInterruptSimulation": 0.1,
                            "atol": SIMULATOR_ATOL, "rtol": SIMULATOR_RTOL, "constant_gas_volume": False},
         "RMG model / simulator settings the deck does not set"),
        ("species_constraints", facts["constraints"], "generatedSpeciesConstraints as given "
                                                      "(--species-constraint; allowed=['input species'] "
                                                      "unless given, so the declared input species stay "
                                                      "admissible), else none"),
        ("profile", "rmg/solver/simulation_1_<N>.csv, highest N",
         "the final-core solver profile (rows at the integrator's steps, amounts in mol); columns matched "
         "to species by structure (chemkin/species_dictionary.txt); a first row not at t = 0 or "
         "conflicting rows at one time are refused"),
        ("time_interpolation", "linear in time",
         "each profile column interpolated between the two rows bracketing the requested time (a row "
         "time takes the row); derived quantities from the interpolated columns; outside the profile "
         "refused"),
    ]
    if "kernels" in facts:
        record = facts["kernels"]["record"]
        rows += [
            ("kernels", {"homolysis": facts["kernels"]["homolysis"],
                         "recombination": facts["kernels"]["recombination"],
                         "transfer": facts["kernels"]["transfer"], "move": facts["kernels"]["move"],
                         "record": record},
             "kernel mode: kineticsFamilies 'none'; the kernel-rates table's homolysis kernel on the chain "
             f"pool and end recombination on the explicit end-radical pool {facts['radical_pool']!r} "
             "(polymer(end_radical_site='primary'), born empty, the chain pool's end groups), per template: "
             "w, a, b (sites a*n + b, n in the table's period, periods_per_unit per declared unit) and "
             "min_fragment from its exact law, which must be the uniform cut; recombination per_feature "
             "sites, joined into one chain of n0 + n1 periods; pair transfer (end + end disproportionation, "
             "its reverse) per kernel and product pair, lengths unchanged, k_pair the constant sum_T "
             "(w/divisor) s0 s1 on every reachable length pair, into the chain, end-radical and "
             f"terminal-alkene pool {facts.get('alkene_pool')!r} (polymer(end_alkene_site='primary'), born "
             "empty, the chain pool's end groups), the interior-radical and internal-alkene pools "
             f"{facts.get('interior_radical_pool')!r}, {facts.get('interior_alkene_pool')!r} "
             "(polymer(interior_feature='radical' / 'alkene'), born empty, the chain pool's end groups); sites "
             "that grow with a chain (H transfer from its interior, end-mid reverse disproportionation) as a "
             "bilinear law c00 + c10 n + c01 m + c11 n m of the pair's lengths, exact on every reachable length "
             "pair; intramolecular H transfer (backbiting) as a pool move k(a n + b) per chain, lengths "
             "unchanged; with record.far_law a law exact only from a longer length is mapped from there and "
             "the shorter lengths recorded (record.short_lengths), else listed; templates it cannot map are "
             "listed in record.unmapped, templates the table gives no law in record.inert (rate 0 in the "
             "table's own law), moves that keep their pool in record.identity; a depository override at "
             "exact sizes (kMC fires it at its own weight) has no moment form: the kernel carries the "
             "template law at every size and each override is listed in record.depository_departures with "
             "its sizes and kmc_over_mom; refusals in "
             "mom_producer.kernel_input; isothermal at the table's temperature (no k(T))"),
            ("kernel_phi", record["phi"],
             "cage-escape factor on the radical-creating (phi_thinned) templates; 1 is the well-mixed bound"),
            ("kernel_volume", record["volume_ratio_mom_over_kmc"],
             "k_pair = (w / divisor) phi^thinned V_poly / V: kMC's pair rate (w / divisor) S0 S1 / V on MoM's "
             "concentrations n / V_poly, V_poly / V this ratio (end groups excluded from V_poly)"),
        ]
    return {"schema": CONVENTIONS_SCHEMA,
            "conventions": [{"name": n, "value": v, "statement": s} for n, v, s in rows],
            "quantity_sources": dict(QUANTITY_SOURCES), "not_produced": dict(NOT_PRODUCED),
            "flags": flags,
            "runner_note": "the flags travel in the --mom file's header: run.json mom.header and the "
                           "mom_flags column of comparison.csv"}


# --------------------------------------------------------------------------
# the run's outputs
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
    for x, y in zip(data, data[1:]):
        if x[0] == y[0] and x != y:
            raise ProducerError(f"profile {path}: duplicate rows at t = {x[0]!r} s with different values")
    if data[0][0] != 0.0:
        raise ProducerError(f"profile {path}: the first row is at t = {data[0][0]!r} s, not 0")
    return header, data


def read_pools(path) -> dict:
    try:
        doc = json.loads(Path(path).read_text())
    except (OSError, ValueError) as error:
        raise ProducerError(f"pool table {path}: {error}") from error
    if isinstance(doc, dict) and doc.get("refused"):
        reason = doc.get("reason", "unspecified")
        species = doc.get("species") or ()
        named = f" for {', '.join(map(str, species))}" if species else ""
        raise ProducerError(f"refused:{reason}{named}")
    pools = doc.get("pools") if isinstance(doc, dict) else None
    if not isinstance(pools, list) or not pools:
        raise ProducerError(f"pool table {path}: no pools")
    for p in pools:
        if not isinstance(p, dict) or not isinstance(p.get("label"), str) or not _positive(p.get("monomer_mw_g_mol")):
            raise ProducerError(f"pool table {path}: each pool needs a label and monomer_mw_g_mol")
        if p.get("explicit_dp"):
            raise ProducerError(f"pool table {path}: pool {p['label']!r} has explicit-DP species; the ledger "
                                "handoff H is not in the profile")
    return doc


def _positive(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x) and x > 0


def without_state(pool: dict) -> bool:
    """Whether a pool-table row is a pool the solver holds no state for: MoM
    lists every registered pool, also one whose proxy and moment dummies never
    reached the core (``phase_species`` and ``bookkeeping_species`` both
    present and empty), and the profile has no column for it.  Such a pool
    must have zero moments (refused otherwise: amount with no state).  A row
    without the keys (an artifact before them) is a pool with state."""
    if "bookkeeping_species" not in pool or "phase_species" not in pool:
        return False
    if pool["bookkeeping_species"] or pool["phase_species"]:
        return False
    moments = pool.get("moments") or ()
    if any(x != 0.0 for x in moments):
        raise ProducerError(f"pool {pool['label']!r} has no core species but moments {list(moments)}")
    return True


def species_key(label: str) -> Tuple:
    """A species' identity across RMG's outputs: its index (``label(N)`` in
    the sidecar and chem.yaml, ``label(N)`` / ``formula(N)`` in Chemkin), else
    its bare label (index-less moment dummies)."""
    m = re.fullmatch(r"(.*)\((\d+)\)", label)
    return ("index", int(m.group(2))) if m else ("label", label)


def _strip_collider(side: str) -> str:
    """A Chemkin side without its trailing fall-off collider ``(+M)`` or
    ``(+<name>)``, where ``<name>`` may itself hold parentheses
    (``P(3)(+H2(2))``)."""
    side = side.strip()
    if side.endswith(")"):
        depth = 0
        for i in range(len(side) - 1, -1, -1):
            if side[i] == ")":
                depth += 1
            elif side[i] == "(":
                depth -= 1
                if depth == 0:
                    return side[:i] if side.startswith("(+", i) else side
    return side


def _split_side(side: str) -> List[str]:
    side = _strip_collider(side)
    if side.endswith("+M"):
        side = side[:-2]
    parts, depth, current = [], 0, ""
    for ch in side:
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        if ch == "+" and depth == 0:
            parts.append(current)
            current = ""
        else:
            current += ch
    parts.append(current)
    return [x.strip() for x in parts if x.strip()]


def _is_number(token: str) -> bool:
    try:
        float(token)
    except ValueError:
        return False
    return True


def read_annotated(path) -> List[dict]:
    """``chem_annotated.inp`` reactions as RMG writes them: consecutive
    entries of one RMG number (an expanded MultiArrhenius writes several) are
    one reaction, with its Chemkin numbers, family (``Template reaction``
    line; library / pdep prefixed; None when the writer gave none: an
    expanded MultiArrhenius), reactant and product keys (a fall-off collider
    dropped) and reversibility, in file order.  The rate coefficients are
    not read here: RMG's own reader gives them (``rmg_reactions``)."""
    try:
        text = Path(path).read_text()
    except OSError as error:
        raise ProducerError(f"annotated mechanism {path}: {error}") from error
    out, pending = [], None
    for line in text.splitlines():
        m = re.match(r"! Reaction index: Chemkin #(\d+); RMG #(-?\d+)", line)
        if m:
            pending = {"chemkin": int(m.group(1)), "rmg": int(m.group(2)), "family": None}
            continue
        if pending is None or not line.strip():
            continue
        if line.startswith("!"):
            m = re.match(r"! (Template|Library|PDep) reaction: (\S+)", line)
            if m and pending["family"] is None:
                pending["family"] = (m.group(2) if m.group(1) == "Template"
                                     else f"{m.group(1).lower()}:{m.group(2)}")
            continue
        tokens = line.split()
        trailing = 0
        while tokens and _is_number(tokens[-1]):
            tokens.pop()
            trailing += 1
        if trailing != 3:
            raise ProducerError(f"annotated mechanism {path}: Chemkin #{pending['chemkin']}: no reaction line")
        equation = "".join(tokens)
        arrow = "<=>" if "<=>" in equation else ("=>" if "=>" in equation else "=")
        if arrow not in equation:
            raise ProducerError(f"annotated mechanism {path}: Chemkin #{pending['chemkin']}: no reaction line")
        left, right = equation.split(arrow, 1)
        entry = {"reactants": [species_key(x) for x in _split_side(left)],
                 "products": [species_key(x) for x in _split_side(right)], "reversible": arrow != "=>"}
        last = out[-1] if out else None
        if (last is not None and last["rmg"] == pending["rmg"] and pending["family"] is None
                and _reaction_key(last["reactants"], last["products"])
                == _reaction_key(entry["reactants"], entry["products"])):
            last["chemkin"].append(pending["chemkin"])
        else:
            out.append({"chemkin": [pending["chemkin"]], "rmg": pending["rmg"], "family": pending["family"],
                        **entry})
        pending = None
    return out


def _reaction_key(reactants, products) -> Tuple:
    return tuple(sorted(reactants)), tuple(sorted(products))


def rmg_reactions(annotated, dictionary, reactions: List[dict]) -> Tuple[Dict[int, Any], List[Any]]:
    """Each Chemkin entry's reaction (by Chemkin number; its kinetics, and its
    species' thermo from the file's THERM block) and the species, as RMG's own
    reader builds them (``rmgpy.chemkin.load_chemkin_file`` with the species
    dictionary beside the file; duplicates kept apart, so an expanded
    MultiArrhenius stays one entry per term).  The reader returns the entries
    in file order: the k-th is ``reactions``' k-th Chemkin number, checked by
    its species."""
    from rmgpy.chemkin import get_species_identifier, load_chemkin_file
    try:
        species, loaded = load_chemkin_file(str(annotated), str(dictionary), check_duplicates=False,
                                             generate_resonance_structures=False)
    except Exception as error:
        raise ProducerError(f"RMG's Chemkin reader refuses {annotated} with {dictionary}: "
                            f"{type(error).__name__}: {error}") from error
    entries = [(c, _reaction_key(r["reactants"], r["products"])) for r in reactions for c in r["chemkin"]]
    if len(loaded) != len(entries):
        raise ProducerError(f"RMG's Chemkin reader gives {len(loaded)} reactions of {annotated}, its annotations "
                            f"{len(entries)}")
    out = {}
    for (number, key), reaction in zip(entries, loaded):
        mine = _reaction_key([species_key(get_species_identifier(x)) for x in reaction.reactants],
                             [species_key(get_species_identifier(x)) for x in reaction.products])
        if mine != key:
            raise ProducerError(f"{annotated}: Chemkin #{number} is {reaction} for RMG's reader, not the annotated "
                                "reaction's species")
        out[number] = reaction
    return out, species


def reaction_rate(terms: Sequence[Any], temperature: float, pressure: float) -> Tuple[Optional[float], Optional[str]]:
    """The rate coefficient of one reaction (the sum of an expanded
    MultiArrhenius' terms, as RMG's MultiArrhenius / MultiPDepArrhenius) by
    RMG's ``get_rate_coefficient(T, P)``, and why it is not usable (None when
    it is): an error RMG raises, a value that is not finite, or a negative
    one."""
    from rmgpy.kinetics import Arrhenius, MultiArrhenius, MultiPDepArrhenius, PDepArrhenius
    try:
        if len(terms) == 1:
            kinetics = terms[0]
        elif all(isinstance(k, Arrhenius) for k in terms):
            kinetics = MultiArrhenius(arrhenius=list(terms))
        elif all(isinstance(k, PDepArrhenius) for k in terms):
            kinetics = MultiPDepArrhenius(arrhenius=list(terms))
        else:
            raise ValueError(f"duplicate entries of mixed kinetics {[type(k).__name__ for k in terms]}")
        k = kinetics.get_rate_coefficient(temperature, pressure)
    except Exception as error:
        return None, f"RMG cannot evaluate the rate coefficient ({type(error).__name__}: {error})"
    if not math.isfinite(k):
        return k, f"the rate coefficient {k!r} is not finite"
    if k < 0.0:
        return k, f"the rate coefficient {k!r} is negative"
    return k, None


def dead(k: Optional[float]) -> bool:
    """A rate coefficient that moves nothing: zero, negative, not finite, or not evaluable."""
    return k is None or not math.isfinite(k) or k <= 0.0


AMBIGUOUS_MAX = 10           # gas-mask ambiguity units whose 2^n gas sets the census still enumerates


def _effective_pressures(efficiencies, gas: Sequence[int], pressure: float,
                         compositions: Sequence[Sequence[Optional[float]]]) -> List[float]:
    """MoM's effective pressure at each composition in time order over the
    core positions ``gas`` (polymer.pyx 6546-6564): P sum(eff y) / sum(y); a
    gas total at or below 1e-20 mol skips the update and keeps the previous
    value, the total pressure of generate_rate_coefficients (4744) before any
    state."""
    out, current = [], pressure
    for y in compositions:
        amounts = [y[j] for j in gas]
        total = math.fsum(amounts)
        if total > 1e-20:
            current = pressure * math.fsum(efficiencies[j] * v for j, v in zip(gas, amounts)) / total
        out.append(current)
    return out


def _units(part) -> frozenset:
    """The ambiguity units in one half of a rule: a unit string, a set of
    them, or none (None, True, "gas", "condensed")."""
    if isinstance(part, str):
        return frozenset() if part in ("gas", "condensed") else frozenset([part])
    if isinstance(part, (set, frozenset, tuple, list)):
        return frozenset(part)
    return frozenset()


def _simplified(rule) -> Tuple:
    """A rule with what it cannot depend on removed: a species a configured
    pool's override condenses, or the first stage condenses outright, is
    condensed under every combination and carries no unit."""
    configured, stage1 = rule
    if configured is True or stage1 == "condensed":
        return True, "condensed"
    return _units(configured) or None, _units(stage1) or "gas"


def _merged(a, b) -> Tuple:
    """Two pools' rules for one species, as MoM applies them: the overrides
    only condense (polymer.pyx 3185-3261) and the first stage is one mask,
    so either condensation stands whatever the pool order."""
    (ca, sa), (cb, sb) = a, b
    configured = True if True in (ca, cb) else (_units(ca) | _units(cb)) or None
    stage1 = "condensed" if "condensed" in (sa, sb) else (_units(sa) | _units(sb)) or "gas"
    return configured, stage1


def _condensed(rule, chosen) -> bool:
    """Whether a core species is condensed in MoM's gas mask under the
    ambiguity units ``chosen`` (``core_phases``): its rule is (configured,
    stage1), configured None (no pool override), True (a configured pool's
    override) or units "configured:<pool>" (any chosen condenses); stage1
    "gas", "condensed" or units "stage1:<species>" (any chosen condenses)."""
    configured, stage1 = rule
    if configured is True or _units(configured) & set(chosen):
        return True
    return stage1 == "condensed" or bool(_units(stage1) & set(chosen))


def mom_rates(group: Sequence[Any], reversible: bool, temperature: float, pressure: float,
              core: Sequence[Tuple[Optional[int], Any, Tuple]], compositions: Sequence[Sequence[float]]) -> dict:
    """Whether one reaction (its Chemkin entries ``group``) moves forward and
    backward as MoM's solver computes its rates (polymer.pyx at b1285f08f).
    kf: ``get_rate_coefficient(T, P)``, a non-finite value taken as 0
    (``generate_rate_coefficients``, 4744-4748); a pressure-dependent
    reaction with collider efficiencies at MoM's effective pressure, MoM
    recomputing it at every state: here at every composition in
    ``compositions`` (profile rows, in time order; ``_effective_pressures``).
    The efficiencies are assigned over the FULL core in MoM's order
    (``get_effective_collider_efficiencies(core_species)``, polymer_input.py
    428-434: each efficiency goes to the first isomorphic core species,
    model.pyx 415-421), and only then the gas mask weighs the amounts
    (``core``, from ``core_phases``, in MoM's core order: with the solver's
    recorded mask every rule is settled and one gas set is evaluated).  Every
    ambiguity unit of a reconstructed mask (a spawned pool configured or not;
    a species whose first-stage phase the sidecar cannot settle) that decides
    the phase of a species holding an amount at some composition, counted on
    the simplified rules (``_simplified``), is taken both ways, every combination (at
    most ``AMBIGUOUS_MAX`` units, else dead with the reason): where the
    forward rate's liveness differs it is dead with the reason "gas mask
    ambiguous for ...".  Compositions are the profile's saved rows (and the
    horizon), samples of the solver's states: a dead verdict and every reason
    say "at the N sampled compositions".  kb = kf / Keq with Keq from the
    species' thermo, and kb = 0 where Keq is not finite or not positive or kb
    is not finite (4749-4760), under every combination likewise: where the
    reverse liveness differs the reverse is dead and the reason a rate
    error.  A direction
    is live when its coefficient is positive and finite at some composition."""
    terms = [r.kinetics for r in group]
    kinetics = terms[0]

    def forward_at(pressures):
        values, why = [], None
        for p in dict.fromkeys(pressures):
            k, reason = reaction_rate(terms, temperature, p)
            values.append(k)
            why = why or reason
        return values, why

    sampled = f"{len(compositions)} sampled compositions"
    out = {"reverse": False, "reverse_reason": None}
    if len(terms) == 1 and kinetics.is_pressure_dependent() and getattr(kinetics, "efficiencies", None):
        efficiencies = kinetics.get_effective_collider_efficiencies([c[1] for c in core])
        rows = [[None if c[0] is None else y[c[0]] for c in core] for y in compositions]
        rows = [[0.0 if v is None else v for v in y] for y in rows]      # a core species without a column holds 0
        # the ambiguity units that decide the phase of a core species holding an amount at some composition,
        # counted on the simplified rules (a species condensed whatever its units brings none)
        rules = [_simplified(c[2]) for c in core]
        present = [j for j in range(len(core)) if any(y[j] != 0.0 for y in rows)]
        units = sorted({u for j in present for part in rules[j] for u in _units(part)})
        if len(units) > AMBIGUOUS_MAX:
            out.update(forward=False, error=f"gas mask ambiguous for {units} (more than {AMBIGUOUS_MAX} to "
                                            f"enumerate) at the {sampled}")
            return out
        variants = []
        for mask in range(2 ** len(units)):
            chosen = {u for b, u in enumerate(units) if mask >> b & 1}
            members = [j for j in range(len(core)) if not _condensed(rules[j], chosen)]
            values, why = forward_at(_effective_pressures(efficiencies, members, pressure, rows))
            variants.append((any(not dead(k) for k in values), values, why))
        if len({v[0] for v in variants}) > 1:
            out.update(forward=False, error=f"gas mask ambiguous for {units}: the rate is live with some of them "
                                            f"in MoM's gas and dead with others at the {sampled}")
            return out
    else:
        values, why = forward_at([pressure])
        units, variants = [], [(any(not dead(k) for k in values), values, why)]
        sampled += " (this rate does not depend on composition)"
    live, forward, error = variants[0]
    out.update(forward=live, error=error,
               dead_reason=None if live else (f"{error} " if error else "") + f"kf is {forward!r} at the {sampled}")
    if reversible:
        try:
            keq = float(group[0].get_equilibrium_constant(temperature))
        except Exception as e:
            out["reverse_reason"] = f"RMG cannot evaluate Keq ({type(e).__name__}: {e})"
            return out
        if not math.isfinite(keq) or keq <= 0.0:
            out["reverse_reason"] = f"Keq = {keq!r}: MoM sets kb = 0 at the {sampled}"
            return out
        backward = []
        for _live, values, _why in variants:            # kb under every gas mask, as the forward rate
            kb = [(k / keq) if (k is not None and math.isfinite(k)) else 0.0 for k in values]
            backward.append((any(not dead(x) for x in kb), kb))
        if len({b[0] for b in backward}) > 1:
            reason = (f"gas mask ambiguous for {units}: the reverse rate kb = kf / Keq is live with some of them "
                      f"in MoM's gas and dead (not finite, MoM's kb = 0) with others at the {sampled}")
            out["reverse_reason"] = reason
            out["error"] = reason if error is None else f"{error}; {reason}"
            return out
        out["reverse"] = backward[0][0]
        if not out["reverse"]:
            out["reverse_reason"] = f"kb = kf / Keq is {backward[0][1]!r} at the {sampled}"
    return out


ACTIVITY_RULE = (
    "a pool-coupled row is live when it is admitted, its rate coefficient (evaluated for admitted rows only, "
    "by RMG's own kinetics objects from rmgpy.chemkin.load_chemkin_file on chem_annotated.inp and the "
    "species dictionary, an expanded MultiArrhenius as the sum of its terms, get_rate_coefficient at the run "
    "temperature and the deck's pressure, a pressure-dependent reaction with collider efficiencies at MoM's "
    "effective pressure at every gas composition the profile visits up to the horizon) is positive and finite "
    "at some such composition (zero, negative, non-finite or unevaluable: dead, each error a stated reason), "
    "and every reactant is present at t = 0 (a profile column > 0, or a pool proxy whose pool has mu0 > 0) or "
    "is made by another live row (admitted sidecar rows and every Chemkin reaction that is no sidecar row, each "
    "only with a live rate coefficient); a reversible row also runs from its products when its reverse "
    "coefficient, kb = kf / Keq with Keq from the species' thermo and kb = 0 where Keq is not finite or not "
    "positive (as MoM's solver), is live; the chemistry is active when a live row belongs to the requested "
    "families AND the profile moves: some polymer pool's mu0, mu1 or mu2, or the amount of a declared species, "
    "interpolated as conversion exports it at a requested observation time (up to the horizon), differs from "
    f"its t = 0 value by more than ACTIVITY_RELATIVE_FLOOR = {ACTIVITY_RELATIVE_FLOOR!r} of that value AND by "
    f"more than an absolute floor of ACTIVITY_MASS_FLOOR_FRACTION = {ACTIVITY_MASS_FLOOR_FRACTION!r} of the "
    "deck's initial polymer mass in grams (a pool moment mu_k weighed by the deck's grams per unit of mu_k, "
    "its initial mass over its initial mu_k: the floor is then that fraction of each moment's own deck total; "
    "a species by its molar mass); bookkeeping proxy columns and other species do not count")


def profile_motion(header, data, pools, facts: dict, dictionary, times: Sequence[float]) -> Optional[dict]:
    """The first physical change of the profile beyond both activity floors
    (``ACTIVITY_RULE``) at the requested observation times, interpolated as
    conversion exports them, or None when nothing moves."""
    floor_g = ACTIVITY_MASS_FLOOR_FRACTION * float(facts["initial_mass_g"])
    per_unit = [float(facts["initial_mass_g"]) / float(m) for m in facts["moments_mol"]]   # g per unit of mu_k
    quantities = []
    for p in pools:
        for k, grams in ((0, per_unit[0]), (1, per_unit[1]), (2, per_unit[2])):
            column = moment_column(f"{p['label']}_mu{k}")
            if column in header:
                quantities.append((f"{p['label']} mu{k}", header.index(column), grams))
    for name, i in sorted(resolve_species(header, dictionary, species_structures(facts)).items()):
        quantities.append((name, i, dictionary[header[i]].get_molecular_weight() * 1000.0))
    first = data[0]
    for t in times:
        row = interpolate(data, t)
        for name, i, grams in quantities:
            change = abs(row[i] - first[i])
            if change > ACTIVITY_RELATIVE_FLOOR * abs(first[i]) and change * grams > floor_g:
                return {"quantity": name, "time_s": t, "change": row[i] - first[i], "change_g": change * grams}
    return None


def _h_loss_radical_daughter(molecule, proxies) -> bool:
    """MoM's structural H-loss radical-daughter test (rmgpy/polymer.py
    ``is_h_loss_radical_daughter`` at b1285f08f, conditions (ii)-(iv), (vi)):
    radical-bearing and neutral, the heavy-atom composition of a proxy, H
    short by its radical count, and, saturated, isomorphic to that proxy under
    any resonance form of either; any structure error: False."""
    try:
        n_rad = molecule.get_radical_count()
        comp = molecule.get_element_count()
        if n_rad < 1 or molecule.get_net_charge() != 0:
            return False
    except Exception:
        return False
    heavy = {el: n for el, n in comp.items() if el != "H"}
    saturated = None
    for proxy in proxies:
        try:
            pcomp = proxy.get_element_count()
        except Exception:
            continue
        if {el: n for el, n in pcomp.items() if el != "H"} != heavy or pcomp.get("H", 0) - comp.get("H", 0) != n_rad:
            continue
        try:
            if saturated is None:
                sat = molecule.copy(deep=True)
                sat.saturate_radicals()
                sat.update()
                saturated = [sat] + list(sat.generate_resonance_structures() or [])
            forms = [proxy] + list(proxy.generate_resonance_structures() or [])
            if any(x.is_isomorphic(y) for x in saturated for y in forms):
                return True
        except Exception:
            return False
    return False


# one implementation, shared with polymer-mom's serializer (rmgpy.polymer)
_same_structure = _species_identity.same_structure


def recorded_molecule(structure) -> Optional[Molecule]:
    """A recorded structure ({"adjacency_list", "multiplicity"}) as RMG's
    molecule, its multiplicity the recorded one; None when it cannot be
    read (absent, malformed, or its multiplicity disagreeing)."""
    if not isinstance(structure, dict) or not isinstance(structure.get("adjacency_list"), str) \
            or not isinstance(structure.get("multiplicity"), int):
        return None
    try:
        molecule = Molecule().from_adjacency_list(structure["adjacency_list"])
    except Exception:
        return None
    return molecule if molecule.multiplicity == structure["multiplicity"] else None


def _recorded_mask(conventions: dict, profile) -> Optional[List[Tuple[str, Optional[str], bool]]]:
    """The solver's gas mask MoM recorded for the profile read
    (conventions.solver_gas_mask: per polymer reactor, its 1-based number,
    the profile it wrote and the core species in solver order, each with its
    structure and gas flag), as (species, structure, gas), or None when the
    sidecar records no mask.  The producer reads reactor 1's profile: a
    record that cannot be bound, or none for reactor 1 naming that profile,
    is refused."""
    record = conventions.get("solver_gas_mask")
    if record is None:
        return None
    reactors = record.get("reactors") if isinstance(record, dict) else None

    def entry_ok(e):
        return isinstance(e, dict) and isinstance(e.get("species"), str) and isinstance(e.get("gas"), bool) \
            and "structure" in e
    if (not isinstance(record, dict) or record.get("source") != "solver" or not isinstance(reactors, list)
            or not all(isinstance(r, dict) and isinstance(r.get("reactor"), int) and isinstance(r.get("profile"), str)
                       and isinstance(r.get("core"), list) and all(entry_ok(e) for e in r["core"])
                       for r in reactors)):
        raise ProducerError("conventions.solver_gas_mask is not a solver mask record bound to its reactor, profile "
                            f"and species structures: {record!r}")
    wanted = f"solver/{Path(profile).name}"
    mine = [r for r in reactors if r["reactor"] == 1 and r["profile"] == wanted]
    if not mine:
        raise ProducerError(f"conventions.solver_gas_mask records reactors and profiles "
                            f"{[(r['reactor'], r['profile']) for r in reactors]}; the producer reads reactor 1's "
                            f"profile {wanted}, which no recorded mask names")
    return [(e["species"], e["structure"], e["gas"]) for e in mine[0]["core"]]


def species_labels(annotated) -> Dict[Tuple, str]:
    """Each species' own RMG label(index), by species key, where the Chemkin
    file's SPECIES block gives it: RMG's verbose writer comments each line
    with it (``<identifier> ! <label(index)>``, chemkin.pyx save_chemkin_file;
    read as RMG's reader reads an alias: one identifier, one comment token)."""
    out, inside = {}, False
    for line in Path(annotated).read_text().splitlines():
        head = line.split("!", 1)[0].split()
        if not inside:
            inside = bool(head) and head[0].upper() in ("SPECIES", "SPEC")
            continue
        if head and head[0].upper() == "END":
            break
        parts = line.split("!")
        if len(parts) == 2 and len(parts[0].split()) == 1 and len(parts[1].split()) == 1:
            out[species_key(parts[0].strip())] = parts[1].strip()
    return out


def core_phases(header: Sequence[str], doc: dict, chemkin_species: Sequence[Any],
                profile, own_labels: Optional[Dict[Tuple, str]] = None) -> Tuple[List[Tuple[Optional[int], Any, Tuple]], str]:
    """MoM's core in its order as (profile column or None, RMG species, rule),
    and the phase source: the order is the profile's columns (the solver
    writes its core in promotion order, rmg/model.py; RMG's Chemkin writer
    sorts the SPECIES block by index, chemkin.pyx ~2153, so the block is not
    it), then any Chemkin species the profile lacks (amount 0).

    ``"solver_mask"``: the sidecar records MoM's own gas mask
    (conventions.solver_gas_mask) for reactor 1's profile read here
    (``_recorded_mask``).  It is the only phase source: each rule is the
    recorded flag, no unit, nothing reconstructed.  Each recorded species
    be the Chemkin species it names by RMG's isomorphism (``_same_structure``
    of ``recorded_molecule``; refused otherwise, and a recorded structure
    that cannot be read equals nothing).  The mask is static
    over the simulation that wrote the profile, whose columns are its core in
    its order, so the record must name the profile's species columns in that
    order (refused otherwise).

    ``"reconstructed"`` (older sidecars): the rule (``_condensed``) rebuilds
    the solver's gas mask from what the sidecar records:

    * the second stage, per pool (``_apply_pool_phase_overrides``, polymer.pyx
      3185-3243 at b1285f08f), condenses the FIRST core species whose base
      label (index stripped) is the pool's label, then breaks, plus the pool's
      mu dummies and explicit-DP species; a configured (root) pool always, a
      spawned pool (conventions.spawned_pools) only if the solver configured
      it, which the sidecar does not record: the unit "configured:<pool>";
    * the first stage (``get_gas_mask``, polymer_input.py 709-781): a root
      pool's mu dummies are condensed, a spawned pool's gas; section D
      (787-791) condenses EVERY core species whose base label is that of an
      H-loss radical daughter of a root pool's proxy (880-990), which the
      sidecar cannot settle (the durable gas veto and Polymer-ness are not
      recorded): each candidate daughter (a pool member or any other core
      species, not a mu dummy, its base no root pool's label) is the unit
      "stage1:<species>", and a pool member's first stage is condensed when
      any candidate of its base is.  The base is the RMG label's, as MoM's:
      a pool member's from the sidecar's name, another species' from its
      own label(index) in ``own_labels`` (the SPECIES block's comments,
      ``species_labels``), never from a Chemkin identifier, which may be a
      formula; a species without either is its own unit and reaches no
      pool's member;
    * a species in several pools takes every pool's rule (``_merged``: any
      condensation stands, whatever the pool order);
    * a species no pool lists: condensed iff conventions.condensed_species
      names it (derived from the solver's own final-core mask; the
      serializer's spawned-pool closure only adds pool members).

    Names and columns are matched by identity, never by raw string: a moment
    by ``moment_column`` of ``<pool>_mu<k>``, a species by ``species_key``.
    A profile column that is no species of the Chemkin file is refused."""
    from rmgpy.chemkin import get_species_identifier
    conventions = doc.get("conventions") if isinstance(doc.get("conventions"), dict) else {}
    condensed = conventions.get("condensed_species")
    if not isinstance(condensed, list):
        raise ProducerError("the pool table has no conventions.condensed_species: MoM's gas set is not recorded")
    labels = {p["label"] for p in doc["pools"]}
    spawned = {str(x) for x in conventions.get("spawned_pools") or ()}

    def identity(name: str):
        m = re.fullmatch(r"(.+)_mu([0-2])", str(name))
        if m and m.group(1) in labels:
            return species_key(moment_column(str(name)))
        return species_key(str(name))

    def base(name: str) -> str:
        m = re.fullmatch(r"(.*)\((\d+)\)", str(name))
        return m.group(1) if m else str(name)

    by_key = {}
    for x in chemkin_species:
        # a long moment dummy: RMG's reader labels it by the SPECIES comment (its full label), the
        # profile column is its 16-character identifier
        key = identity(x.label) if x.index == -1 else species_key(get_species_identifier(x))
        by_key.setdefault(key, x)
    order, unmapped = [], []
    for i, h in enumerate(header):            # MoM's core order: the profile's columns
        if i < 2:
            continue
        if species_key(h) in by_key:
            order.append((i, species_key(h)))
        else:
            unmapped.append(h)
    if unmapped:
        raise ProducerError(f"profile {profile}: columns {unmapped} are neither a pool moment nor a species of the "
                            "Chemkin file")
    placed = {key for _i, key in order}
    order += [(None, key) for key in by_key if key not in placed]    # no profile column: last, amount 0

    recorded = _recorded_mask(conventions, profile)
    if recorded is not None:
        keys = [identity(x) for x, _s, _gas in recorded]
        columns = [species_key(h) for h in header[2:]]
        if keys != columns:
            raise ProducerError(f"conventions.solver_gas_mask names {[x for x, _s, _g in recorded]}, the profile "
                                f"{profile} {list(header[2:])}: the recorded mask is not the mask of the profile's "
                                "core in its order")
        for key, (x, structure, _g) in zip(keys, recorded):
            if key not in by_key:
                continue
            molecule = recorded_molecule(structure)
            if molecule is None:
                raise ProducerError(f"conventions.solver_gas_mask: the structure recorded for {x} cannot be read "
                                    f"({structure!r}); it equals no species")
            verdict = _same_structure([molecule], [m for m in (by_key[key].molecule or []) if m is not None])
            if verdict is _species_identity.UNDECIDABLE:
                raise ProducerError(f"conventions.solver_gas_mask: {x} cannot be compared to the Chemkin species "
                                    "because resonance expansion is unavailable")
            if verdict is _species_identity.MALFORMED:
                raise ProducerError(f"conventions.solver_gas_mask: {x} has a malformed mixed-spin molecule list")
            if verdict is not True:
                raise ProducerError(f"conventions.solver_gas_mask: {x} has another structure in the Chemkin file "
                                    f"than the solver's mask was for (not RMG-isomorphic, resonance and "
                                    f"multiplicity considered)")
        gas = dict(zip(keys, (g for _x, _s, g in recorded)))
        return [(i, by_key[key], (None, "gas" if gas.get(key, True) else "condensed")) for i, key in order], \
            "solver_mask"

    position = {key: j for j, (_i, key) in enumerate(order)}
    for key, species in by_key.items():
        molecules = [m for m in (species.molecule or []) if m is not None]
        verdict = _same_structure(molecules, molecules)
        if verdict is _species_identity.UNDECIDABLE:
            raise ProducerError(f"reconstructed phase: {get_species_identifier(species)} cannot be compared because resonance expansion is unavailable")
        if verdict is _species_identity.MALFORMED:
            raise ProducerError(f"reconstructed phase: {get_species_identifier(species)} has a malformed mixed-spin molecule list")
    # pool members: the mu dummies, the species of the pool's base label, the explicit-DP species
    member_of, names_of = {}, {}
    root_proxies = []
    for p in doc["pools"]:
        label = p["label"]
        mus = {identity(f"{label}_mu{k}") for k in range(3)}
        names = list(p.get("phase_species") or ()) + list(p.get("bookkeeping_species") or ())
        labelled = {identity(x) for x in names if base(x) == label and identity(x) not in mus}
        explicit = {identity(x) for x in names if identity(x) not in mus and identity(x) not in labelled}
        member_of[label] = (mus, labelled, explicit)
        for x in names:
            names_of.setdefault(identity(x), str(x))
        if label not in spawned:
            root_proxies += [by_key[k].molecule[0] for k in labelled if k in by_key and by_key[k].molecule]
    all_mus = {k for mus, _l, _e in member_of.values() for k in mus}
    root_labels = labels - spawned

    # section D: candidate H-loss daughters of a root proxy, by base label (a pool member by its sidecar name)
    daughters: Dict[str, set] = {}
    for key, _j in sorted(position.items(), key=lambda kv: kv[1]):
        x = by_key[key]
        name = names_of.get(key) or (own_labels or {}).get(key)
        mol = x.molecule[0] if x.molecule else None
        if key in all_mus or mol is None or (name is not None and base(name) in root_labels):
            continue
        if name is None:                  # no RMG label: its own unit, which reaches no pool member
            continue
        if _h_loss_radical_daughter(mol, root_proxies):
            daughters.setdefault(base(name), set()).add(f"stage1:{get_species_identifier(x)}")

    rules = {}
    for p in doc["pools"]:
        label = p["label"]
        mus, labelled, explicit = member_of[label]
        configured = True if label not in spawned else f"configured:{label}"
        first = min((k for k in labelled if k in position), key=position.get, default=None)
        for key in mus | labelled | explicit:
            if key not in by_key:
                continue
            if key in mus:
                stage1 = "condensed" if label not in spawned else "gas"
            else:
                stage1 = frozenset(daughters.get(base(names_of.get(key, "")), ())) or "gas"
            condensable = key in mus or key in explicit or key == first
            rule = (configured if condensable else None, stage1)
            rules[key] = _merged(rules[key], rule) if key in rules else rule
    kept = {identity(x) for x in condensed}
    return [(i, by_key[key], rules.get(key, (None, "condensed" if key in kept else "gas"))) for i, key in order], \
        "reconstructed"


def kernel_census(rmg_dir, facts: dict, times: Sequence[float]) -> dict:
    """Kernel mode's activity: the kernels written (no reaction rows: families
    are off) and the profile's first change beyond both floors
    (``profile_motion``), or none."""
    rmg_dir = Path(rmg_dir)
    header, data = read_profile(find_profile(rmg_dir))
    doc = read_pools(rmg_dir / "chemkin" / "polymer_pools.json")
    dictionary = read_species_dictionary(rmg_dir / "chemkin" / "species_dictionary.txt")
    motion = profile_motion(header, data, doc["pools"], facts, dictionary, times)
    return {"mode": "kernels", "kernels": facts["kernels"]["record"]["kernels"],
            "profile_changes": motion is not None, "motion": motion}


def census(rmg_dir, facts: dict, temperature: float, times: Sequence[float]) -> dict:
    """Pool-coupled reactions of the final mechanism per family (rows and
    admitted rows), the live rows and whether the profile moves
    (``ACTIVITY_RULE``).  Each sidecar row is matched one-to-one to its own
    Chemkin reaction: the sidecar and the Chemkin file both list the final
    core reactions in order, so the k-th sidecar row with given reactants and
    products (species by index) is the k-th Chemkin reaction with them; a
    different count is refused.  Family: for a row without a Cantera index,
    the serializer's own id ``d<family>:...`` (refused if its Chemkin
    reaction names another); otherwise its Chemkin reaction's ``Template
    reaction``.  A reaction the Chemkin writer gave no family (an expanded
    MultiArrhenius: plain Reaction entries) is counted by equation under
    ``unattributed``; the input admits no reaction source but the requested
    families (no libraries, no seed, no pdep), so it counts as a
    requested-family row.  Rate coefficients: its Chemkin reaction's, by RMG's
    own reader and kinetics as MoM's solver evaluates them (``rmg_reactions``,
    ``mom_rates``), only for admitted rows; motion: ``profile_motion`` at
    ``times``."""
    rmg_dir = Path(rmg_dir)
    temperature = float(temperature)
    families = facts["families"]
    pressure = Quantity(*facts["pressure"]).value_si
    dictionary_path = rmg_dir / "chemkin" / "species_dictionary.txt"
    doc = read_pools(rmg_dir / "chemkin" / "polymer_pools.json")
    annotated = rmg_dir / "chemkin" / "chem_annotated.inp"
    reactions = read_annotated(annotated)
    rmg, rmg_species = rmg_reactions(annotated, dictionary_path, reactions)
    horizon = float(times[-1])
    header, data = read_profile(find_profile(rmg_dir))
    dictionary = read_species_dictionary(dictionary_path)
    core, phase_source = core_phases(header, doc, rmg_species, find_profile(rmg_dir), species_labels(annotated))
    visited = [r for r in data if r[0] < horizon] + [interpolate(data, horizon)]
    by_key: Dict[Tuple, List[dict]] = {}
    for reaction in reactions:
        by_key.setdefault(_reaction_key(reaction["reactants"], reaction["products"]), []).append(reaction)
    sidecar = list(doc.get("reactions") or ())
    keyed: Dict[Tuple, List[int]] = {}
    for j, row in enumerate(sidecar):
        key = _reaction_key([species_key(x) for x in row.get("reactants") or ()],
                            [species_key(x) for x in row.get("products") or ()])
        keyed.setdefault(key, []).append(j)
    own = {}
    for key, js in keyed.items():
        group = by_key.get(key, [])
        if len(group) != len(js):
            name = (sidecar[js[0]].get("cantera") or {}).get("equation") or sidecar[js[0]].get("id")
            raise ProducerError(f"{len(js)} pool reaction row(s) and {len(group)} reaction(s) of {annotated} "
                                f"have the species of {name}: no one-to-one match")
        own.update(zip(js, group))
    requested = set(families)
    counts: Dict[str, Dict[str, int]] = {}
    unattributed: Dict[str, Dict[str, int]] = {}
    refused_reasons: Dict[str, int] = {}
    rate_errors: Dict[str, str] = {}
    reverse_disabled: Dict[str, str] = {}
    dead_reasons: Dict[str, str] = {}
    rows = []
    for j, row in enumerate(sidecar):
        reaction = own[j]
        cantera = row.get("cantera")
        name = (cantera or {}).get("equation") or str(row.get("id"))
        family = reaction["family"]
        if cantera is None:
            m = re.match(r"d([^:]+):", str(row.get("id")))
            if not m:
                raise ProducerError(f"pool reaction {row.get('id')!r}: no Cantera index and no serialized family")
            serialized = None if m.group(1) == "rxn" else m.group(1)
            if family is not None and serialized is not None and family != serialized:
                raise ProducerError(f"pool reaction {row.get('id')!r}: serialized family {serialized!r}, Chemkin "
                                    f"#{reaction['chemkin'][0]} family {family!r}: no one-to-one match")
            family = serialized if serialized is not None else family
        admitted = not row.get("refused")
        slot = (counts if family is not None else unattributed).setdefault(
            family if family is not None else name, {"rows": 0, "admitted": 0})
        slot["rows"] += 1
        slot["admitted"] += admitted
        if not admitted:
            reason = str(row.get("refused_reason"))
            refused_reasons[reason] = refused_reasons.get(reason, 0) + 1
        kinetics = row.get("kinetics") or {}
        reversible = kinetics.get("reversible")
        if reversible is None:
            reversible = reaction["reversible"]
        zero, backward = None, False
        if admitted:                      # a refused row's rate is never evaluated
            rates = mom_rates([rmg[c] for c in reaction["chemkin"]], bool(reversible), temperature, pressure,
                              core, visited)
            zero, backward = not rates["forward"], rates["reverse"]
            if rates["error"] is not None:
                rate_errors[str(row.get("id"))] = rates["error"]
            if zero:
                dead_reasons[str(row.get("id"))] = rates.get("dead_reason") or rates["error"]
            if reversible and not backward:
                reverse_disabled[str(row.get("id"))] = rates["reverse_reason"]
        rows.append({"id": row.get("id"), "admitted": admitted, "requested": family is None or
                     set(family.split("|")) <= requested, "zero_rate": zero,
                     "reactants": {species_key(x) for x in row.get("reactants") or ()},
                     "products": {species_key(x) for x in row.get("products") or ()},
                     "reversible": bool(reversible), "backward": backward,
                     "proxies": list(row.get("proxy_reactants") or ()) + list(row.get("proxy_products") or ())})
    # presence at t = 0, and whether the profile moves
    first = data[0]
    present = {species_key(h) for h, v in zip(header[2:], first[2:]) if v > 0}
    motion = profile_motion(header, data, doc["pools"], facts, dictionary, times)
    mu0 = {}
    for p in doc["pools"]:
        column = moment_column(f"{p['label']}_mu0")
        mu0[p["label"]] = first[header.index(column)] if column in header else 0.0
    for r in rows:
        for proxy in r["proxies"]:
            if mu0.get(re.sub(r"\(\d+\)$", "", proxy), 0.0) > 0:
                present.add(species_key(proxy))
    own_reactions = {id(x) for x in own.values()}
    movers = [(r, r["reactants"], r["products"], r["backward"]) for r in rows
              if r["admitted"] and not r["zero_rate"]]
    for x in reactions:
        if id(x) in own_reactions:
            continue
        rates = mom_rates([rmg[c] for c in x["chemkin"]], x["reversible"], temperature, pressure, core, visited)
        if rates["error"] is not None:
            rate_errors[f"chemkin #{x['chemkin'][0]}"] = rates["error"]
        if rates["forward"]:
            movers.append((None, set(x["reactants"]), set(x["products"]), rates["reverse"]))
    fired, changed = set(), True
    while changed:
        changed = False
        for j, (_r, left, right, reversible) in enumerate(movers):
            if j not in fired and (left <= present or (reversible and right <= present)):
                fired.add(j)
                present |= left | right
                changed = True
    live = [m[0]["id"] for j, m in enumerate(movers) if j in fired and m[0] is not None and m[0]["requested"]]
    admitted_requested = sum(r["admitted"] and r["requested"] for r in rows)
    return {"families": dict(sorted(counts.items())), "unattributed": dict(sorted(unattributed.items())),
            "requested_families": sorted(requested), "admitted_in_requested_families": admitted_requested,
            "zero_rate_rows": [r["id"] for r in rows if r["admitted"] and r["zero_rate"]],
            "live_rows": live, "profile_changes": motion is not None, "profile_motion": motion,
            "activity_floors": {"relative": ACTIVITY_RELATIVE_FLOOR,
                                "absolute_g": ACTIVITY_MASS_FLOOR_FRACTION * float(facts["initial_mass_g"])},
            "rate_errors": rate_errors, "reverse_disabled": reverse_disabled, "dead_reasons": dead_reasons,
            "temperature_K": temperature, "pressure_Pa": pressure,
            "activity_rule": ACTIVITY_RULE, "refused_reasons": refused_reasons, "phase_source": phase_source,
            "stale_topology": bool((doc.get("conventions") or {}).get("stale_topology"))}


def interpolate(data: List[List[float]], t: float) -> List[float]:
    """The profile columns at time ``t``: a row at ``t`` exactly, else linear in
    time between the last row before ``t`` and the first row after it."""
    times = [r[0] for r in data]
    if not times[0] <= t <= times[-1]:
        raise ProducerError(f"time {t!r} s is outside the profile [{times[0]!r}, {times[-1]!r}] s")
    at = [r for r in data if r[0] == t]
    if at:
        return list(at[0])
    right = bisect.bisect_right(times, t)
    lo, hi = data[right - 1], data[right]
    w = (t - lo[0]) / (hi[0] - lo[0])
    return [x + (y - x) * w for x, y in zip(lo, hi)]


def _end_mass(smiles: Sequence[str]) -> float:
    return math.fsum(Molecule(smiles=s).get_molecular_weight() * 1000.0 for s in smiles)


WHOLE_KEYS = mo.WHOLE_ENSEMBLE_KEYS


def whole_moments(row, first, pools, columns, ends, shift, polymer_cols, n_mass, *, booked: float, e: float,
                  a: float) -> Tuple[List[Optional[float]], Optional[str]]:
    """Whole-ensemble mass-DP moments (``melt_observables.WHOLE_ENSEMBLE_DEFINITION``)
    of one profile row: every pool chain and every polymer-class species.
    mu0 = sum of pool mu0 + species amounts; mu1 from MoM's booked chain mass
    (sum of n_m = (booked - e mu0_pools) / a, exact whatever the ledger);
    mu2 = sum of pool mu2 shifted by each pool's offset (e_p - defect_p - e) / a
    (n_m = pool DP + offset for a chain with its pool's end groups) -- exact
    only when the ledger books no end-group deficit or excess (sum of
    e_p L_p = 0 within MoM's floors); otherwise None and the reason."""
    mu0_pools = math.fsum(row[columns[p["label"]][0]] for p in pools)
    w0 = mu0_pools + math.fsum(row[i] for i in polymer_cols)
    w1 = (booked - e * mu0_pools) / a + math.fsum(row[i] * n_mass[i] for i in polymer_cols)
    ledger, zero = whole_ledger(row, first, pools, columns, ends)
    if not zero:
        return [w0, w1, None], (f"MoM's end-group ledger books {ledger!r} g/mol x mol of end groups minted "
                                "or removed (sum of e_p L_p): how that deficit spreads over chains, hence "
                                "mu2 and Mw, is not in the profile")
    w2 = 0.0
    for p in pools:
        m0, m1, m2 = (row[i] for i in columns[p["label"]])
        d = shift[p["label"]]
        w2 += m2 + 2.0 * d * m1 + d * d * m0
    w2 += math.fsum(row[i] * n_mass[i] ** 2 for i in polymer_cols)
    return [w0, w1, w2], None


def whole_ledger(row, first, pools, columns, ends) -> Tuple[float, bool]:
    """(sum of e_p L_p, whether it is zero within what MoM resolves of each
    pool's mu0 at both times: its floor f plus the simulator's relative
    tolerance, e_p (2 f + rtol (|mu0(0)| + |mu0(t)|)) summed over the pools).
    A conserving kernel set leaves only the integrator's noise (~1e-13 g on
    a run of ~3 mol of chains); a real end-group deficit is events x
    end masses."""
    ledger = math.fsum(ends[p["label"]] * (first[columns[p["label"]][0]] - row[columns[p["label"]][0]])
                       for p in pools)
    return ledger, abs(ledger) <= ledger_threshold(row, first, pools, columns, ends)


def ledger_threshold(row, first, pools, columns, ends) -> float:
    """The largest end-group ledger (g) whole_ledger calls zero at this row:
    e_p (2 f + rtol (|mu0(0)| + |mu0(t)|)) summed over the pools."""
    return math.fsum(ends[p["label"]] * (2.0 * MOM_MOMENT_FLOOR + SIMULATOR_RTOL * (
        abs(first[columns[p["label"]][0]]) + abs(row[columns[p["label"]][0]]))) for p in pools)


def pool_state(row, pools, columns, floor: float = MOM_MOMENT_FLOOR) -> Tuple[Dict[str, List[float]], List[str]]:
    """Each pool's moments under MoM's accepted-state check
    (``_assert_pool_moments_accepted``): a moment below ``-f`` is corruption,
    within ``[-f, 0)`` it is clamped to 0.  After the clamp a pool is corrupt
    when it holds units without chains (mu0 = 0 exactly while mu1 or mu2 is
    beyond ``f``) or when its variance mu0 mu2 - mu1^2 is below minus MoM's
    variance-census tolerance f |mu2| + f |mu0| + 2 f |mu1| + 2 f^2
    (polymer.pyx 5570-5580 at b1285f08f, all three floors f).  A trace
    population whose mu0 is under ``f`` (DP 30 at (5e-15, 1.5e-13, 4.5e-12))
    is a valid state."""
    out, corrupt = {}, []
    for p in pools:
        label = p["label"]
        m = [row[i] for i in columns[label]]
        for k in range(3):
            if m[k] < -floor:
                corrupt.append(f"pool {label} mu{k} = {m[k]!r} mol is below MoM's accepted-state floor -{floor!r} "
                               "(corrupt)")
            elif m[k] < 0.0:
                m[k] = 0.0
        if m[0] == 0.0 and (abs(m[1]) > floor or abs(m[2]) > floor):
            corrupt.append(f"pool {label} holds units without chains (mu0 = {m[0]!r}, mu1 = {m[1]!r}, "
                           f"mu2 = {m[2]!r} mol; corrupt)")
        elif m[0] * m[2] - m[1] * m[1] < -(floor * (abs(m[2]) + abs(m[0]) + 2.0 * abs(m[1])) + 2.0 * floor * floor):
            corrupt.append(f"pool {label} variance mu0 mu2 - mu1^2 = {m[0] * m[2] - m[1] * m[1]!r} is below MoM's "
                           f"variance-census tolerance (mu = {m!r} mol; corrupt)")
        out[label] = m
    return out, corrupt


def pinned_at(m: Sequence[float], n: int, floor: float = MOM_MOMENT_FLOOR) -> bool:
    """Every chain of a pool with moments ``m`` is at DP ``n`` (mu1 = n mu0,
    mu2 = n^2 mu0 within the moments' floors)."""
    return abs(m[1] - n * m[0]) <= floor * (n + 1) and abs(m[2] - n * n * m[0]) <= floor * (n * n + 1)


def _capped_chain(repeat, heads: Sequence[str], n: int) -> Optional[Molecule]:
    """The chain head - n repeat units - tail of a pool's declared structure
    (repeat unit: head port *2, tail port *1, radical at each; end groups:
    mono-radical SMILES), or None when it cannot be built."""
    try:
        chain = Molecule(smiles=heads[0])
        open_atom = [a for a in chain.atoms if a.radical_electrons == 1]
        if len(open_atom) != 1:
            return None
        open_atom = open_atom[0]
        for _ in range(n):
            unit = repeat.copy(deep=True)
            head = next(a for a in unit.atoms if a.label == "*2")
            tail = next(a for a in unit.atoms if a.label == "*1")
            chain = chain.merge(unit)
            chain.add_bond(Bond(open_atom, head, 1))
            open_atom.decrement_radical()
            head.decrement_radical()
            open_atom = tail
        cap = Molecule(smiles=heads[1])
        port = [a for a in cap.atoms if a.radical_electrons == 1]
        if len(port) != 1:
            return None
        chain = chain.merge(cap)
        chain.add_bond(Bond(open_atom, port[0], 1))
        open_atom.decrement_radical()
        port[0].decrement_radical()
        for atom in chain.atoms:
            atom.label = ""
        chain.update_multiplicity()
        chain.update()
        return chain
    except Exception:          # RMG raises several types for an unbuildable structure
        return None


def pool_species_candidates(pools, structures: Dict[str, Molecule], masses: Dict[str, float],
                            deck_repeat, *, kernel_record: Optional[dict] = None
                            ) -> Dict[str, List[Tuple[str, int, float]]]:
    """Per pool, the declared species a chain of it may be (kMC's class
    "species" comes before "polymer"): a species whose molar mass is the
    pool's chain mass a_p n + e_p - defect_p at a whole DP n >= 1.  For a pool
    without defect the chain's structure is built from the pool's declared
    repeat unit (else the deck's) and end groups and must be the species
    (``same_species``); a defect pool keeps no structure, so its formula
    match stands as a candidate that can never be booked (``convert``).
    In kernel mode (``kernel_record`` with ``pool_family``) a pool's chains
    are its family's states, whose end mass is the family's, not the declared
    caps': the candidates are the states kMC identifies as declared species
    (``declared_species_by_pool``, periods), at a whole number of units, each
    checked against the pool's chain mass and refused by name off it."""
    out = {}
    if kernel_record and kernel_record.get("pool_family"):
        ppu = int(kernel_record["moment_unit"]["periods_per_declared_unit"])
        period = float(kernel_record["moment_unit"]["molar_mass_g_mol"])
        for p in pools:
            label = p["label"]
            found = []
            for name, periods in kernel_record["declared_species_by_pool"].get(label, []):
                expected = period * periods + float(kernel_record["pool_end_mass_g_mol"][label])
                if abs(masses[name] - expected) > 1e-6 * expected:
                    raise ProducerError(f"declared species {name!r} is {periods} periods of pool {label!r}'s "
                                        f"family, {expected!r} g/mol, but weighs {masses[name]!r} g/mol")
                if periods % ppu == 0:            # MoM's pools hold whole units only
                    found.append((name, periods // ppu, masses[name]))
            out[label] = found
        return out
    for p in pools:
        a_p = p["monomer_mw_g_mol"]
        e_p = _end_mass(p["end_groups"])
        defect = float(p.get("chain_mass_defect_g_mol") or 0.0)
        repeat = deck_repeat
        if p.get("monomer_adj_list"):
            try:
                repeat = Molecule().from_adjacency_list(p["monomer_adj_list"])
            except Exception:
                repeat = deck_repeat
        found = []
        for name, molecule in structures.items():
            x = (masses[name] - e_p + defect) / a_p
            n = round(x)
            if n < 1 or abs(x - n) > 1e-6:
                continue
            if defect == 0.0:
                chain = _capped_chain(repeat, p["end_groups"], n)
                if chain is not None and not same_species(chain, molecule):
                    continue
            found.append((name, n, masses[name]))
        out[p["label"]] = found
    return out


BREAKING_ARCHETYPES = {
    "discrete_chip/1": ("a chip debits round(chip MW / a) units while the gas gains chip MW: mass (and a bare "
                        "cap at a = 0) is created per event (solver polymer.pyx DISCRETE_CHIP arm)"),
    "legacy_mu1/1": "the legacy flux moves exactly one unit (mass a) per event whatever mass moves (UNRESOLVED)",
    "moment_credit_conduit/1": ("a conduit credit sets n = (M + d) / a: the credited chain's end groups are "
                                "counted as units"),
}
BREAKING_POOL_BLOCKS = {
    "homolysis_initiation": "the k_homolysis kernel mints end-radical chains declaring both caps",
    "side_group_homolysis": "side-group loss leaves n unchanged while mass leaves (Z is not in the profile)",
    "explicit_dp": "the explicit-DP handshake mints capped oligomers outside the mu columns",
}


def pool_coordinate_breaks(doc: dict, label_mass: Dict[str, float], a: float) -> Dict[str, str]:
    """The reason the MoM pool moments of this run stop being moments of the
    mass DP (M - e) / a after t = 0 (pool-coordinate audit of polymer-mom
    b1285f08f), or {} when no such mechanism is live: an admitted row of a
    breaking archetype (``BREAKING_ARCHETYPES``); an admitted migration
    between pools whose per-chain offsets e_p - defect_p differ (no solver
    guard); a pool kernel with a nonzero rate (``channels``: k_scission,
    k_unzip, ...) or a breaking kernel block (``BREAKING_POOL_BLOCKS``); a
    motif pool (``spawn_event_metadata.triggering_dp``: its a is the parent's
    placeholder).  Mechanisms that keep the coordinate: same-pool and
    equal-offset migration rows, volatile ejection (mass-exact, fractional
    shift), scission fragments (exact in mu1 through the ledger; mu2 is
    withheld on a nonzero ledger), defect and daughter pools (offset
    (e_p - defect_p - e) / a)."""
    reasons = []
    pools = {p["label"]: p for p in doc["pools"]}
    for row in doc.get("reactions") or ():
        if row.get("refused"):
            continue
        archetype = row.get("archetype")
        if archetype in BREAKING_ARCHETYPES:
            reasons.append(f"admitted {archetype} row {row.get('id')}: {BREAKING_ARCHETYPES[archetype]}")
        elif archetype == "migration/1":
            src, dst = row.get("src_pool"), row.get("dst_pool")
            if src in label_mass and dst in label_mass and abs(label_mass[src] - label_mass[dst]) > MASS_RELATIVE * a:
                reasons.append(f"admitted migration row {row.get('id')} between pools {src} and {dst} whose "
                               "per-chain offsets e_p - defect_p differ")
    for label, p in pools.items():
        for name, channel in sorted((p.get("channels") or {}).items()):
            rate = channel.get("A") if isinstance(channel, dict) else None
            if rate is None or rate != 0.0:
                reasons.append(f"pool {label} kernel {name} has a nonzero rate ({rate!r})")
        for block, why in BREAKING_POOL_BLOCKS.items():
            if p.get(block):
                reasons.append(f"pool {label} {block}: {why}")
        if (p.get("spawn_event_metadata") or {}).get("triggering_dp") is not None:
            reasons.append(f"pool {label} is a motif pool whose repeat-unit mass is the parent's placeholder")
    if not reasons:
        return {}
    text = "MoM's pool moments leave the mass-DP coordinate: " + "; ".join(reasons)
    return {k: text for k in WHOLE_KEYS}


def realizable(whole: Sequence[float], floor: float = MOM_MOMENT_FLOOR, *, a: float = 1.0,
               e: float = 0.0) -> Tuple[List[float], Optional[str]]:
    """Whole-ensemble moments under MoM's own accepted-state tolerances
    (rmgpy/solver/polymer.pyx ``_assert_pool_moments_accepted``): a moment
    below ``-f`` is corruption, within ``[-f, 0)`` it is clamped to 0; the
    variance ``mu0 mu2 - mu1^2`` below ``-(f |mu2| + f |mu0| + 2 f |mu1| + 2 f^2)``
    (the variance census tolerance) leaves the realizable set, within it
    ``mu2`` is clamped to ``mu1^2 / mu0`` (zero variance).  Returns the
    clamped moments and None, or the moments and the reason Mn / Mw cannot
    be formed.  mu1 is neither required nor clamped non-negative (a mass DP
    (M - e) / a is negative for a molecule lighter than the end pair); the
    cone ``mu1 >= mu0`` is not required; Mn = (a mu1 + e mu0) / mu0 and the Mw
    numerator must be positive; mu0 = 0 with mu1 or mu2 beyond f is
    refused."""
    w = [float(x) for x in whole]
    for k in (0, 2):                     # mu1 is a sum of mass DPs, which may be negative
        if w[k] < -floor:
            return w, f"whole mu{k} = {w[k]!r} mol is below MoM's accepted-state floor -{floor!r}"
        if w[k] < 0.0:
            w[k] = 0.0
    if w[0] == 0.0:
        if abs(w[1]) > floor or w[2] > floor:
            return w, f"whole moments {w} hold mass DP without molecules"
        return w, None
    variance = w[0] * w[2] - w[1] * w[1]
    tolerance = floor * (abs(w[2]) + abs(w[0]) + 2.0 * abs(w[1])) + 2.0 * floor * floor
    if variance < -tolerance:
        return w, (f"whole mu0 mu2 - mu1^2 = {variance!r} is below MoM's accepted-state variance tolerance "
                   f"-{tolerance!r}: no distribution has these moments")
    if variance < 0.0:
        w[2] = w[1] * w[1] / w[0]
    mass = a * w[1] + e * w[0]
    if mass <= 0.0 or a * a * w[2] + 2.0 * a * e * w[1] + e * e * w[0] <= 0.0:
        return w, f"whole Mn / Mw of {w} are not positive (a mu1 + e mu0 = {mass!r})"
    return w, None


def convert(rmg_dir, times: Sequence[float], facts: dict, classes: KmcClasses) -> Tuple[List[dict], List[dict], dict]:
    """(the --mom entries, the diagnostics entries, provenance of the inputs read)."""
    rmg_dir = Path(rmg_dir)
    profile = find_profile(rmg_dir)
    pools_path = rmg_dir / "chemkin" / "polymer_pools.json"
    header, data = read_profile(profile)
    pool_doc = read_pools(pools_path)
    stateless = [p["label"] for p in pool_doc["pools"] if without_state(p)]
    pools = [p for p in pool_doc["pools"] if not without_state(p)]
    dictionary = read_species_dictionary(rmg_dir / "chemkin" / "species_dictionary.txt")
    a = facts["a_g_mol"]
    by_label = {p["label"]: p for p in pools}
    if facts["label"] not in by_label:
        raise ProducerError(f"pool table has no pool {facts['label']!r} (pools {sorted(by_label)})")
    m_u = by_label[facts["label"]]["monomer_mw_g_mol"]
    if abs(m_u - a) > MASS_RELATIVE * a:
        raise ProducerError(f"MoM repeat-unit mass {m_u!r} differs from the deck's {a!r} g/mol")
    columns, ends = {}, {}
    # kernel mode: the pools' end masses from their families (build_input), not their declared end groups
    pool_ends = ((facts.get("kernels") or {}).get("record") or {}).get("pool_end_mass_g_mol") or {}
    for p in pools:
        cols = []
        for k in range(3):
            name = moment_column(f"{p['label']}_mu{k}")
            if name not in header:
                raise ProducerError(f"profile {profile}: no column {name!r} for pool {p['label']!r}")
            cols.append(header.index(name))
        columns[p["label"]] = cols
        groups = p.get("end_groups")
        if not isinstance(groups, list) or len(groups) != 2:
            raise ProducerError(f"pool {p['label']!r}: no end groups in {pools_path}")
        ends[p["label"]] = pool_ends.get(p["label"], _end_mass(groups))
    used = {i for cols in columns.values() for i in cols}
    stray = [h for i, h in enumerate(header) if re.search(r"_mu[0-9]+$", h) and i not in used]
    if stray:
        raise ProducerError(f"profile {profile}: moment columns {stray} belong to no pool of {pools_path}")
    bookkeeping = set(used) | {header.index(s) for p in pools for s in p.get("bookkeeping_species") or ()
                               if s in header}
    species = resolve_species(header, dictionary, species_structures(facts))
    bath = {}
    for label in facts["bath"]:
        found = [i for i, h in enumerate(header) if i >= 2 and re.sub(r"\(\d+\)$", "", h) == label]
        if len(found) != 1:
            raise ProducerError(f"profile {profile}: {len(found)} columns for the bath species {label!r}")
        bath[label] = found[0]
    system = [i for i in range(2, len(header)) if i not in bookkeeping and i not in bath.values()]
    molar = {i: (dictionary[header[i]].get_molecular_weight() * 1000.0 if header[i] in dictionary else None)
             for i in system + list(bath.values())}
    kind = {i: (classes.classify(dictionary[header[i]]) if header[i] in dictionary else "undetermined")
            for i in system}
    polymer_cols = [i for i in system if kind[i] == "polymer"]
    a_deck, e_deck = facts["a_g_mol"], facts["e_g_mol"]
    n_mass = {i: (molar[i] - e_deck) / a_deck for i in polymer_cols}       # mass DP of a polymer-class species
    shift = {p["label"]: (ends[p["label"]] - float(p.get("chain_mass_defect_g_mol") or 0.0) - e_deck) / a_deck
             for p in pools}                                                # pool DP -> mass DP offset
    first = data[0]

    def booked(row):
        total = 0.0
        for p in pools:
            mu0, mu1 = row[columns[p["label"]][0]], row[columns[p["label"]][1]]
            ledger = first[columns[p["label"]][0]] - mu0
            total += (mu1 * p["monomer_mw_g_mol"] - mu0 * float(p.get("chain_mass_defect_g_mol") or 0.0)
                      + ends[p["label"]] * (mu0 + ledger))
        return total

    def mass_of(row, cols):
        missing = [header[i] for i in cols if row[i] != 0.0 and molar[i] is None]
        return (None, missing) if missing else (math.fsum(row[i] * molar[i] for i in cols if row[i] != 0.0), [])

    # the initial inventory: pools, every species column and the bath against the deck
    mu_first = [math.fsum(first[columns[p["label"]][k]] for p in pools) for k in range(3)]
    mass_first = booked(first)
    expected = {i: 0.0 for i in system}
    expected.update({i: facts["initial_moles"].get(name, 0.0) for name, i in species.items()})
    expected.update({i: facts["initial_moles"][label] for label, i in bath.items()})
    off = sorted(header[i] for i, x in expected.items() if _rel(first[i], x) > INITIAL_RELATIVE)
    bad = [(k, x, y) for k, (x, y) in enumerate(zip(mu_first, facts["moments_mol"])) if _rel(x, y) > INITIAL_RELATIVE]
    if bad or off or _rel(mass_first, facts["initial_mass_g"]) > INITIAL_RELATIVE:
        raise ProducerError(f"profile {profile}: the initial inventory (summed moments {mu_first}, booked mass "
                            f"{mass_first!r} g, species columns off: {off}) is not the declared ensemble's "
                            f"({facts['moments_mol']}, {facts['initial_mass_g']!r} g, initial amounts "
                            f"{facts['initial_moles']}) to relative {INITIAL_RELATIVE}")
    gas_first, missing_first = mass_of(first, system)
    system_first = None if gas_first is None else mass_first + gas_first
    bath_first, _m = mass_of(first, list(bath.values()))
    species_mass = {name: m.get_molecular_weight() * 1000.0 for name, m in species_structures(facts).items()}
    candidates = pool_species_candidates(pools, species_structures(facts), species_mass,
                                         classes.observer.spec.legacy_molecules()[0],
                                         kernel_record=(facts.get("kernels") or {}).get("record"))
    light_pools = {}
    for r in pool_doc.get("reactions") or ():
        units = (r.get("params") or {}).get("eject_units")
        if (r.get("refused") or r.get("archetype") != "volatile_ejection/1"
                or not isinstance(units, (int, float)) or r.get("dst_pool") not in columns):
            continue
        # MoM moves eject_units corrected for the pools' defects between source and destination
        # (polymer.pyx 4110-4113 at b1285f08f: a_mom = a_mass - (defect_dst - defect_src) / mw_src, zero
        # within 1e-9 max(1, |a_mass|)); a source the table does not resolve keeps the raw units.
        a_mom = float(units)
        src, dst = by_label.get(r.get("src_pool")), by_label[r["dst_pool"]]
        if src is not None and src is not dst and (src.get("monomer_mw_g_mol") or 0.0) > 0.0:
            a_mom = units - (float(dst.get("chain_mass_defect_g_mol") or 0.0)
                             - float(src.get("chain_mass_defect_g_mol") or 0.0)) / src["monomer_mw_g_mol"]
            if abs(a_mom) <= 1.0e-9 * max(1.0, abs(units)):
                a_mom = 0.0
        if a_mom > 0.0:
            light_pools[r["dst_pool"]] = max(a_mom, light_pools.get(r["dst_pool"], 0.0))
    pool_breaks = pool_coordinate_breaks(
        pool_doc, {p["label"]: ends[p["label"]] - float(p.get("chain_mass_defect_g_mol") or 0.0) for p in pools},
        a_deck)
    entries, diagnostics = [], []
    for t in times:
        row = interpolate(data, t)
        mu = [math.fsum(row[columns[p["label"]][k]] for p in pools) for k in range(3)]
        mass = booked(row)
        gas, missing = mass_of(row, system)
        bath_mass, _m = mass_of(row, list(bath.values()))
        residual = None if (gas is None or system_first is None) else mass + gas - system_first
        bath_residual = None if (bath_mass is None or bath_first is None) else bath_mass - bath_first
        unclassified = sorted(header[i] for i in system if kind[i] == "undetermined" and row[i] != 0.0)
        oligomer_mass, _m = mass_of(row, polymer_cols)
        quantities = {"species_mol": {name: row[i] for name, i in species.items()}}
        closes = (residual is not None and abs(residual) <= MASS_CLOSURE_RELATIVE * system_first
                  and bath_residual is not None and abs(bath_residual) <= MASS_CLOSURE_RELATIVE * bath_first)
        withheld = {}
        pm, corrupt = pool_state(row, pools, columns)
        booked_as, blocked = {}, []
        corrupt_pools = {label for label in pm if any(c.startswith(f"pool {label} ") for c in corrupt)}
        unresolved = {}                       # pool label -> why its chains are not booked
        end_ledger, ledger_zero = whole_ledger(row, first, pools, columns, ends)
        for label, cands in candidates.items():  # classified at every row, t = 0 included
            m = pm[label]
            if not cands or (label not in corrupt_pools and all(abs(x) <= MOM_MOMENT_FLOOR for x in m)):
                continue                      # provably empty
            names = sorted({c[0] for c in cands})
            dps = sorted({c[1] for c in cands})
            pinned = [(name, n) for name, n, _M in cands if pinned_at(m, n)]
            if label in corrupt_pools:
                why = f"pool {label} may hold chains of declared species {names} and its moments are corrupt"
            elif float(by_label[label].get("chain_mass_defect_g_mol") or 0.0) != 0.0:
                why = (f"pool {label} may hold chains of the formula of declared species {names} (DP {dps}); a "
                       "defect pool keeps no structure, so they are not booked by mass alone")
            elif ledger_zero and len(pinned) == 1:
                booked_as[label] = pinned[0]
                continue
            else:
                why = (f"pool {label} may hold chains that kMC classes as declared species {names} (DP {dps}) and "
                       "its moments do not pin them" + ("" if ledger_zero else " (end-group ledger nonzero)"))
            unresolved[label] = why
            blocked.append(why)
        if t > data[0][0]:
            blocked += [f"pool {label} receives volatile ejections of {units} moment units: chains below one "
                        "repeat unit's backbone (kMC class light) may be in it"
                        for label, units in sorted(light_pools.items()) if pm[label][0] > MOM_MOMENT_FLOOR]
        # a species row stands only when every pool that may hold it is booked or provably empty
        species_withheld = {}
        for name in sorted({c[0] for cands in candidates.values() for c in cands}):
            open_pools = [label for label, cands in candidates.items()
                          if label in unresolved and any(c[0] == name for c in cands)]
            if open_pools:
                species_withheld[name] = "; ".join(unresolved[label] for label in open_pools)
                quantities["species_mol"].pop(name, None)
            else:
                into = [pm[label][0] for label, (s_name, _n) in booked_as.items() if s_name == name]
                if into:
                    quantities["species_mol"][name] = quantities["species_mol"].get(name, 0.0) + math.fsum(into)
        reasons = []
        if unclassified:
            reasons.append(f"unclassifiable species present: {unclassified}")
        reasons += corrupt + [b for b in blocked if b not in corrupt]
        polymer_withheld = "; ".join(reasons) or None
        if polymer_withheld is not None:
            withheld.update({k: polymer_withheld for k in WHOLE_KEYS})
        else:
            clamped_row = list(row)
            for label in pm:
                for k, i in enumerate(columns[label]):
                    clamped_row[i] = pm[label][k]
            booked_mu0 = math.fsum(pm[label][0] for label in booked_as)
            booked_mass = math.fsum(pm[label][0] * species_mass[name] for label, (name, _n) in booked_as.items())
            quantities["polymer_derived_molecules_mol"] = (math.fsum(v[0] for v in pm.values()) - booked_mu0
                                                           + math.fsum(row[i] for i in polymer_cols))
            whole, reason = whole_moments(clamped_row, first, pools, columns, ends, shift, polymer_cols, n_mass,
                                          booked=mass, e=e_deck, a=a_deck)
            for label, (name, _n) in booked_as.items():
                n_s = (species_mass[name] - e_deck) / a_deck
                whole[0] -= pm[label][0]
                whole[1] -= pm[label][0] * n_s
                if whole[2] is not None:
                    whole[2] -= pm[label][0] * n_s * n_s
            if not closes:
                for k in ("whole_mu1_mol", "whole_mu2_mol", "whole_Mn_g_mol", "whole_Mw_g_mol"):
                    withheld[k] = "the kMC inventory or the bath does not close: MoM's booked chain mass is not trusted"
            if reason is not None:
                for k in ("whole_mu2_mol", "whole_Mw_g_mol"):
                    withheld.setdefault(k, reason)
            if t > data[0][0]:                    # nothing has moved the pools at the first row
                for k in pool_breaks:
                    withheld.setdefault(k, pool_breaks[k])
            if whole[2] is None:              # Mn needs mu0 and mu1 only: a zero-variance stand-in for mu2
                whole[2] = whole[1] * whole[1] / whole[0] if whole[0] > 0.0 else 0.0
            clamped, bad = realizable(whole, a=a_deck, e=e_deck)
            if bad is not None:
                for k in ("whole_Mn_g_mol", "whole_Mw_g_mol"):
                    withheld.setdefault(k, bad)
            mn, mw = mo.whole_averages(clamped, a_deck, e_deck) if bad is None else (None, None)
            values = {"whole_mu0_mol": clamped[0], "whole_mu1_mol": clamped[1], "whole_mu2_mol": clamped[2],
                      "whole_Mn_g_mol": mn, "whole_Mw_g_mol": mw}
            quantities.update({k: v for k, v in values.items() if k not in withheld})
            if closes and oligomer_mass is not None:
                polymer = mass + oligomer_mass - booked_mass
                quantities["polymer_condensed_mass_g"] = polymer
                quantities["polymer_conversion"] = 1.0 - polymer / facts["initial_mass_g"]
        check_quantity_keys(quantities, facts["species"])
        entries.append({"time_s": t, "quantities": quantities})
        code = math.fsum(row[columns[p["label"]][1]] * p["monomer_mw_g_mol"]
                         - row[columns[p["label"]][0]] * float(p.get("chain_mass_defect_g_mol") or 0.0)
                         for p in pools)
        diagnostics.append({"time_s": t, "volume_m3": row[1], "whole_ensemble_moments_mol": mu,
                            "end_group_ledger_g": end_ledger,
                            "end_group_ledger_threshold_g": ledger_threshold(row, first, pools, columns, ends),
                            "booked_pool_mass_g": mass, "polymer_class_species_mass_g": oligomer_mass,
                            "code_condensed_mass_g": code, "kmc_inventory_species_mass_g": gas,
                            "species_without_structure": missing or missing_first,
                            "unclassified_species": unclassified, "whole_withheld": withheld,
                            "polymer_withheld": polymer_withheld,
                            "pools_booked_as_species": {k: list(v) for k, v in booked_as.items()},
                            "species_withheld": species_withheld,
                            "mass_closure_residual_g": residual, "bath_closure_residual_g": bath_residual,
                            "mass_exported": "polymer_condensed_mass_g" in quantities,
                            "native_Mn_g_mol": a * mu[1] / mu[0] if mu[0] else None,
                            "native_Mw_g_mol": a * mu[2] / mu[1] if mu[1] else None,
                            "pools": {p["label"]: [row[i] for i in columns[p["label"]]] for p in pools}})
    read = {"profile": str(profile), "profile_sha256": _sha256(profile.read_bytes()),
            "profile_rows": len(data), "profile_time_range_s": [data[0][0], data[-1][0]],
            "pools": str(pools_path), "pools_sha256": _sha256(pools_path.read_bytes()),
            "pool_labels": [p["label"] for p in pools], "pools_without_state": stateless,
            "species_columns": {name: header[i] for name, i in species.items()},
            "bath_columns": {label: header[i] for label, i in bath.items()},
            "species_classes": {header[i]: kind[i] for i in system},
            "polymer_species_mass_dp": {header[i]: n_mass[i] for i in polymer_cols},
            "pool_mass_dp_shift": shift,
            "pool_species_candidates": {k: [list(c) for c in v] for k, v in candidates.items()}}
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
def _write_json(path, doc):
    Path(path).write_text(json.dumps(doc, indent=1, sort_keys=True, allow_nan=False) + "\n")


def produce(deck_path, *, checkout, temperature: float, end: float, snapshots: int, families: Sequence[str],
            catalogs: Sequence, out, mode: str = "write", allow_family_mismatch: bool = False,
            allow_inactive: bool = False, database: Optional[str] = None,
            species_constraints: Sequence[str] = (), kernel_rates: Optional[str] = None,
            kernels: Sequence[str] = (), phi: float = 1.0, far_law: bool = False) -> dict:
    """Write (and with ``mode`` "run" / "convert", run and convert); returns ``producer.json``.

    Kernel mode (``kernel_rates``, a kernel-rates table, and ``kernels``, its
    kernel ids): no families and no catalogs; ``kernel_input``'s terms at the
    table's temperature, ``phi`` on radical-creating templates, ``far_law``
    its far-law option; the run is active when its profile moves (no family
    census)."""
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
    kernel_source = None
    if kernel_rates is not None:
        if families or catalogs:
            raise ProducerError("kernel mode writes no families and reads no catalogs (the kernels are its "
                                "chemistry)")
        raw = Path(kernel_rates).read_bytes()
        kernel_source = {"table": json.loads(raw), "ids": list(kernels), "phi": phi, "far_law": bool(far_law),
                         "source": {"table": str(Path(kernel_rates).resolve()), "table_sha256": _sha256(raw)}}
    elif not families:
        raise ProducerError("--families needs at least one family")
    try:
        deck = kr.read_deck(deck_path)
    except kr.RunError as error:
        raise ProducerError(str(error)) from error
    context = deck_context(deck_path)
    if kernel_source is None:
        catalog_families, catalog_databases = catalog_provenance(catalogs)
        family = family_check(catalog_families, families, allow_family_mismatch)
    else:
        catalog_databases, family = [], {"mode": "kernels", "families": []}
    db = database_directory(database, catalog_databases)
    text, facts = build_input(deck, context, temperature=temperature, end=times[-1], families=families,
                              constraints=_parse_constraints(species_constraints), kernels=kernel_source)
    flags = {}
    kinetics = facts["phase_kinetics"]
    if kinetics["not_comparable"]:
        flags["kinetics_not_comparable"] = {
            "ratio_by_species": {n: kinetics["ratio_by_species"][n] for n in kinetics["not_comparable"]},
            "tolerance": kinetics["tolerance"]}
    manifest = {"schema": SCHEMA, "status": "written", "mode": mode,
                "deck": {"path": deck.path, "sha256": deck.sha256, **deck.digests()},
                "catalogs": [{"path": str(p), "sha256": _sha256(Path(p).read_bytes())} for p in catalogs],
                "family_check": family, "temperature_K": temperature, "times_s": times,
                "mom_checkout": checkout_identity(checkout), "input_sha256": _sha256(text.encode()),
                "facts": facts, "flags": flags}

    def finish(status=None):
        if status is not None:
            manifest["status"] = status
        _write_json(out / "conventions.json",
                    conventions(deck, facts, temperature=temperature, family=family, database=db, flags=flags))
        _write_json(out / "producer.json", manifest)

    rmg_dir.mkdir(parents=True, exist_ok=True)
    if mode == "convert":
        if (rmg_dir / "input.py").read_text() != text:
            raise ProducerError(f"{rmg_dir / 'input.py'} is not the input these arguments write")
    else:
        (rmg_dir / "input.py").write_text(text)
        (rmg_dir / "rmgrc").write_text(f"database.directory : {db[0]}\n")
    finish()
    if mode == "write":
        return manifest
    if mode == "run":
        run_rmg(rmg_dir, checkout)
        manifest["status"] = "ran"
    rmgpy_path = rmg_dir / "rmgpy_path.txt"
    manifest["rmgpy_imported"] = rmgpy_path.read_text().strip() if rmgpy_path.exists() else None
    for stale in ("mom.json", "mom_diagnostics.json"):
        if (out / stale).exists():
            (out / stale).unlink()
    # the run's integrity first (its refusals take precedence), then its activity
    entries, diagnostics, read = convert(rmg_dir, times, facts, KmcClasses(deck.spec, deck.xs))
    if kernel_source is not None:
        manifest["census"] = kernel_census(rmg_dir, facts, times)
        if not manifest["census"]["profile_changes"]:
            flags["inactive_chemistry"] = manifest["census"]
            if not allow_inactive:
                finish(STATUS_INACTIVE)
                raise ProducerError(f"{STATUS_INACTIVE}: the profile of the kernels "
                                    f"{facts['kernels']['record']['kernels']} does not move ({ACTIVITY_RULE}); "
                                    "pass --allow-inactive to convert it anyway, flagged")
    else:
        manifest["census"] = census(rmg_dir, facts, temperature, times)
    if kernel_source is None and (not manifest["census"]["live_rows"] or not manifest["census"]["profile_changes"]):
        flags["inactive_chemistry"] = manifest["census"]
        if not allow_inactive:
            finish(STATUS_INACTIVE)
            raise ProducerError(f"{STATUS_INACTIVE}: no pool-coupled reaction of the families "
                                f"{facts['families']} is live in the MoM mechanism, or its profile "
                                f"does not move (live rows {manifest['census']['live_rows']}, profile "
                                f"changes {manifest['census']['profile_changes']}) "
                                f"({manifest['census']['admitted_in_requested_families']} admitted; census "
                                f"{manifest['census']['families']}, unattributed "
                                f"{manifest['census']['unattributed']}; {ACTIVITY_RULE}); pass "
                                "--allow-inactive to convert it anyway, flagged")
    finish()                                   # conventions.json as the header digests it
    raw_diagnostics = (json.dumps(diagnostics, indent=1, allow_nan=False) + "\n").encode()
    header = {"flags": flags, "conventions_sha256": _sha256((out / "conventions.json").read_bytes()),
              "diagnostics_sha256": _sha256(raw_diagnostics),
              "run": {"temperature_K": temperature, "pressure": list(facts["pressure"]), "horizon_s": times[-1]},
              "producer": {"schema": SCHEMA, "deck_sha256": deck.sha256, "input_sha256": manifest["input_sha256"],
                           "mom_checkout": manifest["mom_checkout"]}}
    (out / "mom.json").write_text(json.dumps({"header": header, "entries": entries}, indent=1, allow_nan=False)
                                  + "\n")
    (out / "mom_diagnostics.json").write_bytes(raw_diagnostics)
    manifest.update({"read": read, "mom_file": str(out / "mom.json"),
                     "mom_sha256": _sha256((out / "mom.json").read_bytes())})
    finish("converted")
    return manifest


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m rmgpy.kmc.mom_producer", description=__doc__.split("\n\n")[0])
    parser.add_argument("deck")
    parser.add_argument("--mom-checkout", required=True)
    parser.add_argument("--temperature", type=float, required=True)
    parser.add_argument("--end", type=float, required=True)
    parser.add_argument("--snapshots", type=int, required=True)
    parser.add_argument("--families", nargs="+", default=[])
    parser.add_argument("--catalog", nargs="+", default=[])
    parser.add_argument("--kernel-rates", help="kernel mode: a kernel-rates table (no --families, no --catalog)")
    parser.add_argument("--kernel", action="append", default=[], help="kernel mode: a kernel id (repeatable)")
    parser.add_argument("--phi", type=float, default=1.0,
                        help="kernel mode: cage-escape factor of radical-creating templates (1: well mixed)")
    parser.add_argument("--kernel-far-law", action="store_true",
                        help="kernel mode: map a law exact only from a longer length from there, its shorter "
                             "lengths recorded (default: listed unmapped)")
    parser.add_argument("--out", required=True)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--run", action="store_const", dest="mode", const="run")
    group.add_argument("--write-only", action="store_const", dest="mode", const="write")
    group.add_argument("--convert", action="store_const", dest="mode", const="convert")
    parser.add_argument("--allow-family-mismatch", action="store_true")
    parser.add_argument("--allow-inactive", action="store_true")
    parser.add_argument("--database")
    parser.add_argument("--species-constraint", action="append", default=[])
    args = parser.parse_args(argv)
    try:
        manifest = produce(args.deck, checkout=args.mom_checkout, temperature=args.temperature, end=args.end,
                           snapshots=args.snapshots, families=args.families, catalogs=args.catalog,
                           out=args.out, mode=args.mode or "write",
                           allow_family_mismatch=args.allow_family_mismatch, allow_inactive=args.allow_inactive,
                           database=args.database, species_constraints=args.species_constraint,
                           kernel_rates=args.kernel_rates, kernels=args.kernel, phi=args.phi,
                           far_law=args.kernel_far_law)
    except (ProducerError, kr.RunError) as error:
        print(f"refused: {error}", file=sys.stderr)
        return 2
    print(f"{manifest['status']}: {args.out}" + (f" -> {manifest['mom_file']}" if "mom_file" in manifest else "")
          + (f" flags {sorted(manifest['flags'])}" if manifest["flags"] else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
