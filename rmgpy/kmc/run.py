"""Run the melt kMC from one RMG deck and emitted template catalogs (I-124).

Command line::

    python -m rmgpy.kmc.run DECK CATALOG [CATALOG ...] --temperature T
        (--times t0 t1 ... | --end T_END --snapshots N) [--max-events N]
        --seeds S [S ...] --out DIR [--allow-incomplete] [--allow-unbound-catalog]
        [--allow-unmatched-catalog] [--mom FILE] [--jobs N]

**Deck.** An RMG input file, read without executing it: it is parsed with
``ast`` and must contain exactly one ``polymer(...)`` call, as a module-level
statement; the name ``polymer`` may appear nowhere else (no alias, no call
inside a conditional, loop, function or expression) and nothing may bind it
(def / class names, assignment or loop targets, imports including ``*``,
global / nonlocal, parameters, ``except ... as``).  Its arguments bind like
``rmgpy.rmg.input.polymer``'s (``POLYMER_PARAMETERS``), ``*args`` / ``**kw``
are refused, and the arguments read are ``ast.literal_eval``-ed: a non-literal
one is refused.  ``linear_homopolymer_spec`` (a literal mapping or JSON
string) is normalized by ``rmgpy.polymer_spec.normalized_linear_homopolymer_spec``
(the input file's own reader), a schema-2 ``polymer_spec`` by
``normalized_polymer_spec``; a deck carries exactly one of the two.  A missing
declaration, both, a label that differs from the declaration's, or legacy
chemistry arguments beside it are refused.  A ``polymer_spec`` deck (linear
statistical copolymer) runs on the melt sampled from its declared sequence
statistics (``melt_sequence``); ``run.json`` carries its declared and
realized composition (``composition``).  It has no MoM counterpart yet: a
``--mom`` file is refused and no comparison is written
(``refused:schema_v2_mom``, ``comparison_refused``).
From the same call the run reads ``cutoff`` (the MoM hybrid threshold ``xs``,
default as in ``polymer()``) and ``side_group_homolysis`` (the MoM side-group
channels), so the observables are partitioned exactly as the MoM run of the
same deck.  The rest of the deck is never evaluated.

**Declaration digests.**  ``declaration_sha256`` = sha256 of the
declaration's ``normalized_json`` (``LinearHomopolymerSpec`` or
``PolymerSpec``: the bytes of the normalized declaration).  ``chemistry_sha256`` = sha256 of the compact sorted-key JSON
(``separators=(",", ":")``) of the normalized declaration MINUS the fields in
``CHEMISTRY_EXCLUDED`` (``chemistry_digest``); every other field, including
any field the declaration grows later, enters it.  Excluded, each because it
does not enter catalog emission (``rmgpy.kmc.emit`` reads only the repeat
unit and the roster's systems) and is a size, scale or runtime block of one
particular run:

* ``label``: a display string; no code reads chemistry from it;
* ``initial_ensemble`` (dp_counts; for a ``PolymerSpec`` also composition,
  sequence statistics, sampler and seed), ``initial_mobile_counts``:
  population sizes and the initial state; templates are local windows and
  emission builds its reference systems from the roster's ``systems``, not
  from the ensemble;
* ``initial_mass_kg``: the macroscopic scale (``scale_factor``) only;
* ``Mn_g_mol``, ``Mw_g_mol``: summary statistics of the declared ensemble
  (validated against it);
* ``density_kg_m3``, ``density_basis``: the volume, which enters the engine
  as ``N_A V`` in bimolecular propensities, never a template weight;
* ``transport_model``: the runtime MET block the engine reads from the deck
  (emission rosters carry ``{}``);
* ``observable_relative_tolerance``: a tolerance of the observables checks;
* ``diagnostics``: a free-form block that no emission, engine or
  observables code reads.

Included therefore: schema_version, repeat_unit, end_groups, small_species
and monomer_product (declared species change species identity and MET
bypass, hence rates), locality_policy, family_manifest, library_manifest;
for a ``PolymerSpec`` also architecture, residues, links, sequence_support
(it decides which local environments exist) and monomer_products.
Both digests are written to ``run.json`` (``deck``) and are the values an
emitter writes into a catalog's ``provenance`` under the same names.

**Catalogs.** Each file is ``kmc-template-catalog/1`` (``rmgpy.kmc.emit``)
with ``provenance``, ``payloads``, ``gates.failures`` and the ``complete``
flag; a flag that contradicts the gate failures is refused.  Refused unless
flagged, and then labelled in every output file:

* ``complete: false`` -> ``--allow-incomplete``; label ``catalog_complete``
  (and ``catalog_gate_failures``);
* a catalog whose provenance carries no declaration digest ->
  ``--allow-unbound-catalog``; label ``catalog_bound``.  A digest that
  differs from the deck's is refused whatever the flags (``bound_by`` lists
  the digests that matched);
* no loaded template matches any site of the initial state (no unimolecular
  transition, no bimolecular participant site) ->
  ``--allow-unmatched-catalog``; label ``catalog_matches_initial_state``
  (``initial_match`` counts the matched templates and sites).

**Combining catalogs.**  One emission owns each reaction family (the
emitter's coverage gate makes the templates within one catalog disjoint), so
several catalogs are combined only when every catalog's
``provenance.families`` is a non-empty list and the lists are pairwise
disjoint; shared families are refused by name, and a catalog without
``provenance.families`` cannot be combined with another (alone it runs).
Nested windows of one family (a window and the same window grown by a
shell) would otherwise both fire the same reactions.  Two further guards:

* the initial-state site check (more than one catalog): on the initial
  melt, no reaction site may be matched by templates of two catalogs
  (``cross_catalog_site``; the first such site is reported).  A
  unimolecular site is the transition's concrete edit set; a bimolecular
  one is the pair of participant sites with each participant's concrete
  change (its atom edits, intra-participant bond edits and the order of the
  bonds it forms to the partner), matched for both participants (in either
  role).  Evidence on the initial state, not a proof for later states;

Always refused: loaded template names that collide, and two loaded
templates with the same executable, within a catalog or across catalogs (it
would fire twice).  Executables are compared as the ``window_matcher.Template``
the engine loads, by labelled-graph isomorphism (``template_identity``):
atoms labelled by element, kind, sorted exterior, isotope, charge, radicals,
lone pairs and final atom edits; every bonded, edited or same-participant
pair labelled by (same participant, bond order before, after).  Invariant to
the order of atoms, bonds and edits and to exchanging the two participants.  Each
catalog must load at the run temperature (no extrapolation).

**Run.** Isothermal.  Per seed ``n`` (sequentially, in one process, unless
``--jobs``): ``MeltEngine(spec, catalog, T, numpy.random.default_rng(n))`` and
``melt_observables.write_time_series`` with the invariant checker on, into
``<out>/seed-<n>/observables.jsonl``; every observation line carries, as it is
written, ``seed`` and the labels above, and the file ends with the time
series' stop record (``melt_observables.STOP_SCHEMA``).  The optional event cap
is run-wide per seed and is passed to ``write_time_series(max_events=...)``:
when it is reached the seed stops and a requested time the engine has not
reached is not observed (``observations`` counts the observation lines).  The
seed's stop reason is the one ``write_time_series`` returns (``max_events``
when the cap left a requested time unreached, else ``horizon`` or
``absorbing``; ``melt_observables.STOP_REASONS``), recorded as returned;
an exception in a seed (engine, observation or report) stops it with stop
reason ``error`` (``error`` names it; the traceback goes to stderr), the other
seeds still run and the run's status is ``failed``.

**Manifest** ``<out>/run.json`` (``kmc-melt-run/1``)::

    schema, status ("running" | "complete" | "failed"[, error])
    deck: {path, sha256, label, declaration (normalized mapping),
           declaration_sha256, chemistry_sha256, xs, side_group_homolysis}
    catalogs: [{path, sha256, schema, complete, gate_failures, bound_by,
                provenance, templates (count at T), template_names_sha256}]
    catalog_complete, catalog_gate_failures: [{catalog, path, sha256, failures}]
    catalog_bound, catalog_matches_initial_state, initial_match: {templates, sites}
    temperature_K, isothermal, times_s, max_events, seeds, jobs, observables_schema
    code: {sha (git HEAD of the checkout), dirty, dirty_paths}
    environment: {python, numpy, PYTHONHASHSEED}
    runs: [{seed, path, observables_sha256, observations, stop_reason, [error,]
            engine_time_s, events, null_proposals, wall_s,
            report (engine report without its leak entries, with its ``thinning``
            decomposition: proposals, accepted, phi_sum and escape probability
            per radical-increasing class; null for a failed seed)}]
    censoring: {seeds_requested, censored_times_s,
                stopped_early: [{seed, stop_reason, observations,
                                 last_observed_time_s, engine_time_s}]}
    leak: {estimator, definition, note, excluded_seeds (seeds that failed),
           per_seed: [{seed, leak_num, leak_den, non_thinned_integral, phi_sum,
                       thinned_proposals, open_thinned, open_phi_source}],
           pooled (melt_engine.pooled_leak_fraction over per_seed)}
    mom: null | {path, sha256, times, header (the --mom file's header, null for the list forms)}
    comparison: "comparison.csv"
    rows: "rows.json" (or rows_refused: why the rows could not be mapped)
    wall_s

The run's leak number is ``leak.pooled.fraction`` only: per seed the
manifest carries the additive parts, never a single-run fraction.  Rerunning
with the same inputs gives byte-identical observation files and a manifest
that differs only in the ``wall_s`` entries.

**Comparison** ``<out>/comparison.csv`` (``kmc-melt-comparison/2``; a first
``#`` line carries the schema and the labels): one row per observation time
and ``melt_observables.side_by_side`` quantity.  ``seeds_requested``;
``seeds`` = the seeds that reached the time with a value; ``censored`` =
fewer seeds reached the time than were requested.  An uncensored row gives
``kmc_mean``, ``kmc_stderr`` (sample standard deviation / sqrt(seeds); empty
for one seed), the MoM value, ``difference = kmc_mean - mom`` and
``relative_difference`` = difference / |mom|.  A censored row leaves
``kmc_mean``, ``kmc_stderr``, ``difference`` and ``relative_difference``
empty and gives the statistics of the seeds that reached the time only as
``censored_kmc_mean`` / ``censored_kmc_stderr``.  ``explicit_dp_mol.*`` and
``species_mol.*`` absent from a seed count as 0 (as in ``side_by_side``).
``mom_flags``: the names of the MoM run's flags (the header form's
``header.flags``), sorted and ";"-joined, on every row.  Without ``--mom`` the
MoM columns are empty.

**Rows** ``<out>/rows.json``: the run's comparison rows
(``comparison_rows.kmc_rows``, ``kmc-comparison-rows/1``: element inventory
and residual, mass, declared species, the ladder and length distribution in
primitive periods, whole-ensemble moments), joined with a kernel reference's
or a MoM file's rows by ``python -m rmgpy.kmc.comparison_rows join``.

``--mom FILE``: JSON (a list) or JSON lines of ``{"time_s": t, "quantities":
{...}}`` with ``quantities`` keyed like ``mom_reported_quantities`` (the
format ``side_by_side`` takes): every value a finite JSON number (not a
string, boolean or null) or a mapping of such numbers.  A MoM time matches the
observation time within relative ``MOM_TIME_RELATIVE`` (1e-12, float
round-off only) of it; a finite time matching none, or more than one, is
refused, and so is a non-finite one.  Header form: one JSON object
``{"header": {"flags": {name: details}, ...}, "entries": [...]}`` (what
``rmgpy.kmc.mom_producer`` writes), refused unless ``header.producer.deck_sha256``
is the sha256 of the deck being run, ``header.conventions_sha256`` that of
the ``conventions.json`` beside the file, and ``header.run`` the run's
temperature and last observation time (``horizon_s``); the entries as above, the header copied
whole into ``run.json`` (``mom.header``) and its flag names onto every
comparison row.

Nothing here names a polymer: the chemistry is the deck's declaration and the
catalogs' templates.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import io
import json
import math
import os
import platform
import statistics
import subprocess
import sys
import time
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy

from rmgpy.kmc import catalog_set as cs
from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import melt_observables as mo
from rmgpy.polymer_spec import (LinearHomopolymerSpec, PolymerSpec, normalized_linear_homopolymer_spec,
                                normalized_polymer_spec)

RUN_SCHEMA = "kmc-melt-run/1"
COMPARISON_SCHEMA = "kmc-melt-comparison/2"
CATALOG_SCHEMA = "kmc-template-catalog/1"     # rmgpy.kmc.emit.SCHEMA
CATALOG_KEYS = ("schema", "provenance", "payloads", "gates", "complete")
# rmgpy.rmg.input.polymer's parameters, in order (positional binding), and its cutoff default
POLYMER_PARAMETERS = (
    "label", "monomer", "end_groups", "cutoff", "Mn", "Mw", "moments", "initial_mass", "k_unzip",
    "k_scission", "discrete_dp_threshold", "monomer_product", "radical_qssa_unzip", "explicit_dp",
    "k_homolysis", "k_depropagation", "side_group_homolysis", "thermal_analysis_inputs",
    "linear_homopolymer_spec", "polymer_spec")
POLYMER_CUTOFF_DEFAULT = 4
LEGACY_CHEMISTRY = ("monomer", "end_groups", "moments", "Mn", "Mw", "monomer_product")
# the declaration fields that do NOT enter chemistry_sha256 (module docstring: why each one)
CHEMISTRY_EXCLUDED = ("label", "initial_ensemble", "initial_mobile_counts", "initial_mass_kg", "Mn_g_mol",
                      "Mw_g_mol", "density_kg_m3", "density_basis", "transport_model",
                      "observable_relative_tolerance", "diagnostics")
BINDING_KEYS = ("chemistry_sha256", "declaration_sha256")
LEAK_NOTE = ("the run's leak number is pooled.fraction (sum of leak_num over seeds / sum of leak_den); "
             "per seed only the additive parts are reported; seeds that failed are excluded")
LEAK_PARTS = ("leak_num", "leak_den", "non_thinned_integral", "phi_sum", "thinned_proposals",
              "open_thinned", "open_phi_source")
COMPARISON_COLUMNS = ("time_s", "quantity", "tier", "seeds_requested", "seeds", "censored", "kmc_mean",
                      "kmc_stderr", "censored_kmc_mean", "censored_kmc_stderr", "mom", "difference",
                      "relative_difference", "note", "mom_flags")
MOM_TIME_RELATIVE = 1.0e-12


class RunError(ValueError):
    """An input of the run is refused (deck, catalog, times, seeds, MoM file, output directory)."""


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _plain(value):
    """``value`` as plain JSON data (refused if it is not JSON-representable)."""
    try:
        return json.loads(json.dumps(value, sort_keys=True, allow_nan=False))
    except (TypeError, ValueError) as error:
        raise RunError(f"not JSON-representable: {error}") from error


def declaration_sha256(spec: Union[LinearHomopolymerSpec, PolymerSpec]) -> str:
    return _sha256(spec.normalized_json)


def chemistry_digest(normalized: Dict[str, Any]) -> str:
    """sha256 of a normalized declaration mapping minus ``CHEMISTRY_EXCLUDED``."""
    return _sha256(json.dumps({k: v for k, v in normalized.items() if k not in CHEMISTRY_EXCLUDED},
                              sort_keys=True, separators=(",", ":")).encode())


def chemistry_sha256(spec: Union[LinearHomopolymerSpec, PolymerSpec]) -> str:
    return chemistry_digest(json.loads(spec.normalized_json))


# --------------------------------------------------------------------------
# the deck
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Deck:
    path: str
    sha256: str
    spec: Union[LinearHomopolymerSpec, PolymerSpec]
    label: str
    xs: int
    side_group_homolysis: Optional[list]

    def digests(self) -> Dict[str, str]:
        return {"chemistry_sha256": chemistry_sha256(self.spec),
                "declaration_sha256": declaration_sha256(self.spec)}

    def manifest(self) -> dict:
        return dict({"path": self.path, "sha256": self.sha256, "label": self.label,
                     "declaration": json.loads(self.spec.normalized_json),
                     "xs": self.xs, "side_group_homolysis": self.side_group_homolysis}, **self.digests())


def _binds_polymer(node) -> bool:
    """Whether ``node`` binds (or, for ``import *``, may bind) the name ``polymer``."""
    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
        return node.name == "polymer"
    if isinstance(node, ast.Name):
        return node.id == "polymer" and isinstance(node.ctx, (ast.Store, ast.Del))
    if isinstance(node, (ast.Import, ast.ImportFrom)):
        return any(alias.name == "*" or (alias.asname or alias.name.split(".")[0]) == "polymer"
                   for alias in node.names)
    if isinstance(node, (ast.Global, ast.Nonlocal)):
        return "polymer" in node.names
    if isinstance(node, ast.arg):
        return node.arg == "polymer"
    if isinstance(node, ast.ExceptHandler):
        return node.name == "polymer"
    for kind, field in (("MatchAs", "name"), ("MatchStar", "name"), ("MatchMapping", "rest")):
        cls = getattr(ast, kind, None)
        if cls is not None and isinstance(node, cls):
            return getattr(node, field, None) == "polymer"
    return False


def _is_none(node) -> bool:
    return isinstance(node, ast.Constant) and node.value is None


def read_deck(path) -> Deck:
    """The one ``linear_homopolymer_spec`` declaration of an RMG deck, read
    from its syntax tree without executing anything (module docstring)."""
    path = Path(path)
    try:
        raw = path.read_bytes()
        tree = ast.parse(raw.decode("utf-8"), filename=str(path))
    except (OSError, UnicodeDecodeError, SyntaxError, ValueError) as error:
        raise RunError(f"deck {path} could not be read: {type(error).__name__}: {error}") from error
    bindings = sorted({getattr(node, "lineno", 0) for node in ast.walk(tree) if _binds_polymer(node)})
    if bindings:
        raise RunError(f"deck {path}: the name 'polymer' is bound or may be bound (line {bindings}: a def, "
                       "class, assignment, import, global/nonlocal or parameter); the deck is read "
                       "without executing it and polymer must be the input file's own")
    statements = [stmt.value for stmt in tree.body
                  if isinstance(stmt, ast.Expr) and isinstance(stmt.value, ast.Call)
                  and isinstance(stmt.value.func, ast.Name) and stmt.value.func.id == "polymer"]
    allowed = {id(call.func) for call in statements}
    stray = [node for node in ast.walk(tree)
             if isinstance(node, ast.Name) and node.id == "polymer" and id(node) not in allowed]
    if stray:
        lines = sorted({node.lineno for node in stray})
        raise RunError(f"deck {path}: 'polymer' is used other than as one module-level polymer(...) "
                       f"statement (line {lines}); the deck is read without executing it")
    if not statements:
        raise RunError(f"deck {path} declares no polymer(); the melt kMC needs one polymer() call "
                       "carrying linear_homopolymer_spec")
    if len(statements) > 1:
        raise RunError(f"deck {path} declares {len(statements)} polymer() calls (lines "
                       f"{[c.lineno for c in statements]}); the melt kMC runs exactly one "
                       "linear_homopolymer_spec declaration")
    call = statements[0]
    if any(isinstance(a, ast.Starred) for a in call.args) or any(k.arg is None for k in call.keywords):
        raise RunError(f"deck {path}: polymer() with *args or **kwargs is refused (arguments must be literal)")
    if len(call.args) > len(POLYMER_PARAMETERS):
        raise RunError(f"deck {path}: polymer() takes at most {len(POLYMER_PARAMETERS)} positional arguments")
    nodes = dict(zip(POLYMER_PARAMETERS, call.args))
    for keyword in call.keywords:
        if keyword.arg not in POLYMER_PARAMETERS:
            raise RunError(f"deck {path}: polymer() got an unexpected keyword argument {keyword.arg!r}")
        if keyword.arg in nodes:
            raise RunError(f"deck {path}: polymer() got multiple values for argument {keyword.arg!r}")
        nodes[keyword.arg] = keyword.value

    def literal(name, default=None):
        node = nodes.get(name)
        if node is None:
            return default
        try:
            return ast.literal_eval(node)
        except (ValueError, TypeError, SyntaxError, MemoryError, RecursionError) as error:
            raise RunError(f"deck {path}: polymer() argument {name} (line {node.lineno}) is not a literal; "
                           "the deck is read without executing it") from error

    if "label" not in nodes:
        raise RunError(f"deck {path}: polymer() is missing its label")
    label = literal("label")
    declared = {name: literal(name) for name in ("linear_homopolymer_spec", "polymer_spec")}
    declared = {name: value for name, value in declared.items() if value is not None}
    if not declared:
        raise RunError(f"deck {path}: polymer({label!r}) carries no linear_homopolymer_spec or polymer_spec")
    if len(declared) > 1:
        raise RunError(f"deck {path}: polymer({label!r}) carries both linear_homopolymer_spec and "
                       "polymer_spec; a deck declares exactly one")
    (keyword, value), = declared.items()
    legacy = [k for k in LEGACY_CHEMISTRY if k in nodes and not _is_none(nodes[k])]
    if legacy:
        raise RunError(f"deck {path}: polymer({label!r}) combines {keyword} "
                       f"with legacy arguments {legacy}")
    normalize = normalized_polymer_spec if keyword == "polymer_spec" else normalized_linear_homopolymer_spec
    try:
        spec = normalize(value)
    except (TypeError, ValueError, KeyError) as error:
        raise RunError(f"deck {path}: {keyword} refused: {error}") from error
    if label != spec.label:
        raise RunError(f"deck {path}: polymer label {label!r} differs from the declaration "
                       f"label {spec.label!r}")
    xs = literal("cutoff", POLYMER_CUTOFF_DEFAULT)
    if isinstance(xs, bool) or not isinstance(xs, int) or xs < 0:
        raise RunError(f"deck {path}: polymer cutoff must be a non-negative integer, got {xs!r}")
    sgh = literal("side_group_homolysis")
    if sgh is not None:
        try:
            sgh = _plain(sgh)
        except RunError as error:
            raise RunError(f"deck {path}: side_group_homolysis {error}") from error
    return Deck(str(path), _sha256(raw), spec, spec.label, xs, sgh)


# --------------------------------------------------------------------------
# catalogs
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class CatalogRecord:
    path: str
    sha256: str
    provenance: Any
    complete: bool
    gate_failures: Tuple[str, ...]
    template_names: Tuple[str, ...]
    bound_by: Tuple[str, ...] = ()

    def manifest(self) -> dict:
        return {"path": self.path, "sha256": self.sha256, "schema": CATALOG_SCHEMA,
                "complete": self.complete, "gate_failures": list(self.gate_failures),
                "bound_by": list(self.bound_by),
                "provenance": self.provenance, "templates": len(self.template_names),
                "template_names_sha256": _sha256("\n".join(self.template_names).encode())}


def _bound_by(path, provenance, deck: Deck, allow_unbound: bool) -> Tuple[str, ...]:
    """The declaration digests the catalog carries and the deck matches;
    a digest that differs from the deck's is refused."""
    mine = deck.digests()
    carried = provenance if isinstance(provenance, dict) else {}
    out = []
    for key in BINDING_KEYS:
        if key not in carried:
            continue
        if carried[key] != mine[key]:
            raise RunError(f"catalog {path}: provenance {key} {str(carried[key])[:16]}... differs from the "
                           f"deck's {key} {mine[key][:16]}...: the catalog was emitted for another declaration")
        out.append(key)
    if not out and not allow_unbound:
        raise RunError(f"catalog {path} carries no declaration digest ({' / '.join(BINDING_KEYS)} in its "
                       "provenance), so nothing ties it to the deck's polymer; pass --allow-unbound-catalog "
                       "to run it, labelled catalog_bound=false")
    return tuple(out)


def expand_catalog_sets(paths: Sequence) -> List[Tuple[Path, Optional[bytes]]]:
    """(path, bytes) for each catalog to load: every ``kmc-catalog-set/1``
    file is replaced by its member catalogs with the very bytes whose sha256,
    shared inputs, collision-entry keys and dependencies were verified (a set
    loads exactly as its members given one by one, and a member replaced on
    disk after the check is never read).  A plain catalog comes with the bytes
    read here, or None when unreadable (``load_catalogs`` reports it)."""
    out: List[Tuple[Path, Optional[bytes]]] = []
    for path in paths:
        path = Path(path)
        try:
            raw = path.read_bytes()
            doc = json.loads(raw)
        except (OSError, ValueError):
            out.append((path, None))             # unreadable: load_catalogs reports it
            continue
        if not cs.is_catalog_set(doc):
            out.append((path, raw))
            continue
        try:
            out.extend(cs.read_set_members(path, doc))
        except cs.CatalogSetError as error:
            raise RunError(str(error)) from error
    return out


def load_catalogs(paths: Sequence, temperature: float, *, allow_incomplete: bool = False,
                  deck: Optional[Deck] = None, allow_unbound: bool = False
                  ) -> Tuple[me.PayloadCatalog, List[CatalogRecord]]:
    """The union of the catalogs as one ``PayloadCatalog`` and one record per
    file.  Refused: a wrong schema or missing field, a ``complete`` flag that
    contradicts the gate failures, an incomplete catalog without
    ``allow_incomplete``, a catalog that does not load at ``temperature`` or
    loads no template, template names or executables shared by two emitted
    templates; with ``deck``, a declaration digest that differs from the
    deck's, and a catalog without one unless ``allow_unbound``.  Several
    catalogs (plain, from one set or from several) must together pass the
    set checks (``catalog_set.check_union``: shared inputs, owned families
    and records, collision-entry keys, links), whatever ``allow_incomplete``."""
    if not paths:
        raise RunError("at least one catalog file is needed")
    sources = expand_catalog_sets(paths)
    temperature = float(temperature)
    records: List[CatalogRecord] = []
    payloads: List[dict] = []
    owner: Dict[str, str] = {}
    executables: Dict[str, List[Tuple[Any, str]]] = {}
    verified: List[Tuple[Path, bytes]] = []
    for path, raw in sources:
        path = Path(path)
        try:
            if raw is None:
                raw = path.read_bytes()
            doc = json.loads(raw)
        except (OSError, ValueError) as error:
            raise RunError(f"catalog {path}: unreadable ({error})") from error
        if not isinstance(doc, dict) or doc.get("schema") != CATALOG_SCHEMA:
            schema = doc.get("schema") if isinstance(doc, dict) else None
            raise RunError(f"catalog {path}: schema {schema!r} is not {CATALOG_SCHEMA!r}")
        missing = [k for k in CATALOG_KEYS if k not in doc]
        if missing:
            raise RunError(f"catalog {path}: missing fields {missing} (a catalog carries its "
                           "provenance, gates and complete flag)")
        complete = doc["complete"]
        gates = doc["gates"]
        failures = gates.get("failures") if isinstance(gates, dict) else None
        if not isinstance(complete, bool):
            raise RunError(f"catalog {path}: the complete flag must be true or false, got {complete!r}")
        if not isinstance(failures, list) or not all(isinstance(f, str) for f in failures):
            raise RunError(f"catalog {path}: gates.failures must be a list of strings")
        if complete == bool(failures):
            raise RunError(f"catalog {path}: complete={complete} contradicts its {len(failures)} gate "
                           "failures")
        if not complete and not allow_incomplete:
            raise RunError(f"catalog {path} is incomplete (complete: false; gate failures: {failures}); "
                           "pass --allow-incomplete to run it, labelled incomplete")
        bound_by = _bound_by(path, doc["provenance"], deck, allow_unbound) if deck is not None else ()
        if not isinstance(doc["payloads"], list):
            raise RunError(f"catalog {path}: payloads must be a list")
        try:
            templates = me.PayloadCatalog(doc["payloads"]).templates(temperature)
            names = [t.name for t in templates]
        except Exception as error:   # EmissionError, malformed payloads
            raise RunError(f"catalog {path} does not load at {temperature} K: {error}") from error
        if not names:
            raise RunError(f"catalog {path} loads no template at {temperature} K")
        for name in names:
            if name in owner:
                raise RunError(f"catalogs {owner[name]} and {path} collide on template name {name!r}")
            owner[name] = str(path)
        for t in templates:
            graph = _template_graph(t)
            bucket = executables.setdefault(_graph_hash(graph), [])
            for other_graph, where in bucket:
                if _same_graph(graph, other_graph):
                    raise RunError(f"two loaded templates have the same executable (up to the order of atoms, "
                                   f"bonds and edits and participant exchange): {where} and {t.name!r} of "
                                   f"{path}; loading both would fire it twice")
            bucket.append((graph, f"{t.name!r} of {path}"))
        payloads += doc["payloads"]
        verified.append((path, raw))
        records.append(CatalogRecord(str(path), _sha256(raw), _plain(doc["provenance"]), complete,
                                     tuple(failures), tuple(names), bound_by))
    if len(records) > 1:
        families = []
        for r in records:
            f = r.provenance.get("families") if isinstance(r.provenance, dict) else None
            if not isinstance(f, list) or not f or not all(isinstance(x, str) and x for x in f):
                raise RunError(f"catalog {r.path} carries no provenance.families (a non-empty list of family "
                               "names), so it cannot be combined with another catalog: catalogs combine only "
                               "when their families are pairwise disjoint")
            families.append(set(f))
        for i in range(len(records)):
            for j in range(i + 1, len(records)):
                shared = sorted(families[i] & families[j])
                if shared:
                    raise RunError(f"catalogs {records[i].path} and {records[j].path} share the families "
                                   f"{shared}; one emission owns each family, so catalogs combine only when "
                                   "their provenance.families are pairwise disjoint")
        # the set checks hold for the union of everything loaded together
        # (plain catalogs, one set or several sets); allow_incomplete relaxes
        # completeness only, never these
        try:
            cs.check_union([cs.member_from_bytes(path, raw) for path, raw in verified])
        except cs.CatalogSetError as error:
            raise RunError(f"the catalogs loaded together are not one consistent set: {error}") from error
    return me.PayloadCatalog(payloads), records


def _template_graph(t):
    """The labelled graph of a loaded ``window_matcher.Template`` (module
    docstring): equal executables give isomorphic graphs whatever their order."""
    import networkx as nx
    final_atom, final_bond = {}, {}
    for edit in t.edits:
        if hasattr(edit, "order"):
            final_bond[(min(edit.a, edit.b), max(edit.a, edit.b))] = float(edit.order)
        else:
            final_atom[edit.atom] = (edit.radical, edit.charge, edit.lone_pairs)
    g = nx.Graph()
    for i, a in enumerate(t.atoms):
        g.add_node(i, lab=repr((a.element, a.kind, tuple(sorted(a.exterior)), a.isotope, a.charge, a.radical,
                                a.lone_pairs, final_atom.get(i))))
    n = len(t.atoms)
    for i in range(n):
        for j in range(i + 1, n):
            same = t.atoms[i].participant == t.atoms[j].participant
            before = float(t.bond.get((i, j), 0.0))
            after = final_bond.get((i, j), before)
            if same or before or after:
                g.add_edge(i, j, lab=repr((same, before, after)))
    return g


def _graph_hash(g) -> str:
    import networkx as nx
    return nx.weisfeiler_lehman_graph_hash(g, node_attr="lab", edge_attr="lab", iterations=4)


def _same_graph(g, h) -> bool:
    from networkx.algorithms import isomorphism as iso
    return iso.GraphMatcher(g, h, node_match=lambda x, y: x["lab"] == y["lab"],
                            edge_match=lambda x, y: x["lab"] == y["lab"]).is_isomorphic()


def template_identity(t) -> str:
    """A hash of a loaded template's executable (equal for equal executables;
    a collision is resolved by isomorphism in ``load_catalogs``)."""
    return _graph_hash(_template_graph(t))


def _edit_key(edit) -> Tuple:
    if hasattr(edit, "a") and hasattr(edit, "b"):
        a, b = sorted((edit.a, edit.b))
        return (type(edit).__name__, a, b, float(edit.old), float(edit.new))
    return (type(edit).__name__,) + tuple(getattr(edit, k) for k in ("atom", "old", "new"))


def _participant_change(t, participant: int, site) -> frozenset:
    """The concrete change a bimolecular template makes at one participant site."""
    where = dict(zip(t.edited[participant], site))
    out = set()
    for edit in t.edits:
        if hasattr(edit, "order"):
            before = float(t.bond.get((min(edit.a, edit.b), max(edit.a, edit.b)), 0.0))
            ends = [x for x in (edit.a, edit.b) if x in where]
            if len(ends) == 2:
                a, b = sorted((where[edit.a], where[edit.b]))
                out.add(("bond", a, b, before, float(edit.order)))
            elif ends:
                out.add(("partner", where[ends[0]], before, float(edit.order)))
        elif edit.atom in where:
            atom = t.atoms[edit.atom]
            for field in ("radical", "charge", "lone_pairs"):
                new = getattr(edit, field)
                if new is not None:
                    out.add((field, where[edit.atom], getattr(atom, field), new))
    return frozenset(out)


def cross_catalog_site(engine: me.MeltEngine, owner: Dict[str, int]) -> Optional[str]:
    """The first reaction site of the engine's current state matched by
    templates of two catalogs (``owner``: template name -> catalog index), or
    None (module docstring)."""
    uni: Dict[frozenset, Tuple[int, str]] = {}
    pair_sites: Dict[int, List[set]] = {}
    for cid in sorted(engine.entries):
        entry = engine.entries[cid]
        per_class = [(ci, entry.uni[ci][0]) for ci in sorted(entry.uni)] + \
            [(ci, entry.uni_o[(ci, o)][0]) for ci, o in sorted(entry.uni_o)]
        for ci, transitions in per_class:
            name = engine.templates[ci].name
            for tr in transitions:
                key = frozenset(_edit_key(e) for e in tr.edits)
                seen = uni.setdefault(key, (owner[name], name))
                if seen[0] != owner[name]:
                    atoms = sorted({x for e in tr.edits for x in (getattr(e, "a", None), getattr(e, "b", None),
                                                                     getattr(e, "atom", None)) if x is not None})
                    return (f"atoms {atoms} (one transition, edits {sorted(key)}) are matched by template "
                            f"{seen[1]!r} of catalog {seen[0]} and template {name!r} of catalog {owner[name]}")
        for (ci, participant, _kind), sites in sorted(entry.sites.items(), key=lambda kv: repr(kv[0])):
            t = engine.templates[ci]
            slot = pair_sites.setdefault(ci, [set(), set()])
            for site in sites:
                slot[participant].add(_participant_change(t, participant, site))
    names = sorted(pair_sites, key=lambda ci: engine.templates[ci].name)
    for x, ci in enumerate(names):
        for cj in names[x + 1:]:
            a, b = engine.templates[ci].name, engine.templates[cj].name
            if owner[a] == owner[b]:
                continue
            s, u = pair_sites[ci], pair_sites[cj]
            for p in (0, 1):                 # u's participants in the same or the exchanged role
                first, second = s[0] & u[p], s[1] & u[1 - p]
                if first and second:
                    site = sorted(min(first, key=repr))
                    return (f"a bimolecular site (participant change {site}, with a matching partner) is "
                            f"matched by template {a!r} of catalog {owner[a]} and template {b!r} of catalog "
                            f"{owner[b]}")
    return None


def initial_match(engine: me.MeltEngine) -> dict:
    """The templates (classes) matching some site of the engine's current
    state: a unimolecular transition (allowed or leak) or a bimolecular
    participant site; ``sites`` counts both."""
    classes, sites = set(), 0
    for entry in engine.entries.values():
        for (ci, *_o), (allowed, leak) in [((ci,), v) for ci, v in entry.uni.items()] + \
                list(entry.uni_o.items()):
            if allowed or leak:
                classes.add(ci)
                sites += len(allowed) + leak
        for (ci, _participant, _kind), found in entry.sites.items():
            if found:
                classes.add(ci)
                sites += len(found)
    return {"templates": len(classes), "sites": sites}


# --------------------------------------------------------------------------
# times, seeds, MoM file
# --------------------------------------------------------------------------
def observation_times(times: Optional[Iterable[float]] = None, end: Optional[float] = None,
                      snapshots: Optional[int] = None) -> List[float]:
    """The observation times: given, or ``snapshots`` evenly spaced times from
    0 to ``end`` inclusive.  Finite, non-negative and strictly increasing."""
    if times is not None and (end is not None or snapshots is not None):
        raise RunError("give the observation times or an end time with a snapshot count, not both")
    if times is None:
        if end is None or snapshots is None:
            raise RunError("observation times are needed: --times, or --end with --snapshots")
        end = float(end)
        if isinstance(snapshots, bool) or not isinstance(snapshots, int) or snapshots < 2:
            raise RunError(f"the snapshot count must be an integer >= 2, got {snapshots!r}")
        if not math.isfinite(end) or end <= 0.0:
            raise RunError(f"the end time must be finite and positive, got {end!r}")
        out = [end * k / (snapshots - 1) for k in range(snapshots)]
        out[-1] = end
    else:
        out = [float(t) for t in times]
    if not out:
        raise RunError("at least one observation time is needed")
    if any(not math.isfinite(t) or t < 0.0 for t in out):
        raise RunError(f"observation times must be finite and non-negative: {out}")
    if any(b <= a for a, b in zip(out, out[1:])):
        raise RunError(f"observation times must be strictly increasing: {out}")
    return out


def _seeds(seeds) -> List[int]:
    out = list(seeds)
    if not out:
        raise RunError("at least one seed is needed")
    for s in out:
        if isinstance(s, bool) or not isinstance(s, int) or s < 0:
            raise RunError(f"seeds must be non-negative integers, got {s!r}")
    if len(set(out)) != len(out):
        raise RunError(f"duplicate seed in {out}: each seed is one independent run")
    return out


def _finite_number(value) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _quantities_problem(quantities: dict) -> Optional[str]:
    """Why MoM quantities are not finite JSON numbers or mappings of them (None: they are)."""
    for key, value in quantities.items():
        if isinstance(value, dict):
            for sub, v in value.items():
                if not _finite_number(v):
                    return f"{key}.{sub} = {v!r}"
        elif not _finite_number(value):
            return f"{key} = {value!r}"
    return None


def _non_finite(value, where: str = "header") -> List[str]:
    """The paths of every number in a decoded JSON value that is not finite
    (JSON decoding accepts NaN and Infinity)."""
    if isinstance(value, float) and not math.isfinite(value):
        return [f"{where} = {value!r}"]
    if isinstance(value, dict):
        return [x for k, v in value.items() for x in _non_finite(v, f"{where}.{k}")]
    if isinstance(value, list):
        return [x for i, v in enumerate(value) for x in _non_finite(v, f"{where}[{i}]")]
    return []


def load_mom(path, times: Sequence[float], *, deck_sha256: Optional[str] = None,
             temperature_K: Optional[float] = None) -> Tuple[List[Optional[dict]], dict]:
    """Per observation time the MoM quantities (None where the file has none);
    matching in the module docstring.  A header-form file is bound to the
    deck being run (``deck_sha256``, the sha256 of its bytes: refused when
    absent or different from ``header.producer.deck_sha256``) and to the
    ``conventions.json`` beside it (refused when missing or when its sha256
    is not ``header.conventions_sha256``) and to the run's conditions
    (``header.run``: ``temperature_K`` equal to ``temperature_K``,
    ``horizon_s`` the last observation time within ``MOM_TIME_RELATIVE``;
    the pressure is the deck's, bound by its digest)."""
    path = Path(path)
    try:
        raw = path.read_bytes()
        text = raw.decode("utf-8")
    except (OSError, UnicodeDecodeError) as error:
        raise RunError(f"MoM file {path}: {error}") from error
    try:
        doc = json.loads(text)
    except ValueError:
        try:
            doc = [json.loads(line) for line in text.splitlines() if line.strip()]
        except ValueError as error:
            raise RunError(f"MoM file {path}: neither JSON nor JSON lines ({error})") from error
    header = None
    if isinstance(doc, dict) and set(doc) == {"header", "entries"}:
        header, doc = doc["header"], doc["entries"]
        if not isinstance(header, dict) or not isinstance(header.get("flags"), dict):
            raise RunError(f"MoM file {path}: the header must be an object with an object of flags "
                           "({\"flags\": {name: details}, ...})")
        if not isinstance(doc, list):
            raise RunError(f"MoM file {path}: entries must be a list of {{time_s, quantities}} entries")
        non_finite = _non_finite(header)
        if non_finite:
            raise RunError(f"MoM file {path}: header values {non_finite} are not finite numbers")
        producer = header.get("producer")
        bound = producer.get("deck_sha256") if isinstance(producer, dict) else None
        if deck_sha256 is None or bound != deck_sha256:
            raise RunError(f"MoM file {path}: the header's deck digest {bound!r} is not the deck being run "
                           f"({deck_sha256!r}): the MoM run is not of this deck")
        conventions = path.parent / "conventions.json"
        try:
            digest = _sha256(conventions.read_bytes())
        except OSError as error:
            raise RunError(f"MoM file {path}: no conventions.json beside it to bind the header to "
                           f"({error})") from error
        if header.get("conventions_sha256") != digest:
            raise RunError(f"MoM file {path}: the header's conventions digest {header.get('conventions_sha256')!r} "
                           f"is not the sha256 of {conventions} ({digest})")
        conditions = header.get("run") if isinstance(header.get("run"), dict) else {}
        mom_t = conditions.get("temperature_K")
        if temperature_K is None or isinstance(mom_t, bool) or not isinstance(mom_t, (int, float)) \
                or float(mom_t) != float(temperature_K):
            raise RunError(f"MoM file {path}: the MoM run's temperature {mom_t!r} K is not the run's "
                           f"({temperature_K!r} K)")
        horizon = conditions.get("horizon_s")
        last = float(times[-1]) if len(times) else None
        if (last is None or isinstance(horizon, bool) or not isinstance(horizon, (int, float))
                or abs(float(horizon) - last) > MOM_TIME_RELATIVE * max(abs(float(horizon)), abs(last))):
            raise RunError(f"MoM file {path}: the MoM run's horizon {horizon!r} s is not the run's last "
                           f"observation time ({last!r} s)")
    elif isinstance(doc, dict):
        doc = [doc]
    if not isinstance(doc, list):
        raise RunError(f"MoM file {path}: expected a list of {{time_s, quantities}} entries")
    out: List[Optional[dict]] = [None] * len(times)
    for entry in doc:
        if (not isinstance(entry, dict) or set(entry) != {"time_s", "quantities"}
                or not isinstance(entry["quantities"], dict)
                or isinstance(entry["time_s"], bool) or not isinstance(entry["time_s"], (int, float))):
            raise RunError(f"MoM file {path}: each entry is {{\"time_s\": t, \"quantities\": {{...}}}}")
        t = float(entry["time_s"])
        if not math.isfinite(t):
            raise RunError(f"MoM file {path}: time {t} is not finite")
        problem = _quantities_problem(entry["quantities"])
        if problem is not None:
            raise RunError(f"MoM file {path}: at time {t} s, {problem} is not a finite JSON number "
                           "(strings, booleans, null and lists are refused)")
        match = [k for k, ti in enumerate(times) if abs(t - ti) <= MOM_TIME_RELATIVE * max(abs(t), abs(ti))]
        if not match:
            raise RunError(f"MoM file {path}: time {t!r} s is not an observation time of the run {list(times)} "
                           f"(relative tolerance {MOM_TIME_RELATIVE})")
        if len(match) > 1:
            raise RunError(f"MoM file {path}: time {t!r} s is ambiguous: observation times "
                           f"{[times[k] for k in match]} are all within relative {MOM_TIME_RELATIVE}")
        k = match[0]
        if out[k] is not None:
            raise RunError(f"MoM file {path}: two entries at time {times[k]} s")
        out[k] = entry["quantities"]
    return out, {"path": str(path), "sha256": _sha256(raw), "times": sum(q is not None for q in out),
                 "header": header}


# --------------------------------------------------------------------------
# code identity
# --------------------------------------------------------------------------
def code_identity() -> dict:
    """``git rev-parse HEAD`` of the checkout this module runs from, and
    whether ``rmgpy/`` differs from it (tracked changes or untracked files)."""
    package = Path(__file__).resolve().parents[1]

    def git(*args):
        try:
            done = subprocess.run(["git", "-C", str(package), *args], capture_output=True, text=True,
                                  timeout=60)
        except (OSError, subprocess.SubprocessError):
            return None
        return done.stdout if done.returncode == 0 else None

    head = git("rev-parse", "HEAD")
    status = git("status", "--porcelain", "--untracked-files=normal", "--", ".")
    paths = sorted(line[3:] for line in status.splitlines() if line.strip()) if status is not None else None
    return {"sha": head.strip() if head else None,
            "dirty": None if paths is None else bool(paths),
            "dirty_paths": None if paths is None else paths[:50]}


# --------------------------------------------------------------------------
# one seed
# --------------------------------------------------------------------------
class _LabelledObserver(mo.MeltObserver):
    """The observer ``write_time_series`` would build, adding the run's labels
    to every observation as it is written."""

    def __init__(self, *args, label: Dict[str, Any], **kwargs):
        super().__init__(*args, **kwargs)
        self.label = dict(label)

    def observe(self, target, releases=None) -> dict:
        obs = super().observe(target, releases=releases)
        obs.update(self.label)
        return obs


@dataclass(frozen=True)
class _Context:
    spec: Union[LinearHomopolymerSpec, PolymerSpec]
    catalog: me.PayloadCatalog
    temperature: float
    times: Tuple[float, ...]
    xs: int
    side_group_homolysis: Optional[list]
    max_events: Optional[int]
    out: Path
    label: Tuple[Tuple[str, Any], ...]


def _observation_lines(text: str) -> List[str]:
    """The observation lines of a time series (typed records such as the stop
    record carry a ``record`` key and are not observations)."""
    return [line for line in text.splitlines() if "record" not in json.loads(line)]


def _run_seed(context: _Context, seed: int) -> dict:
    start = time.perf_counter()
    directory = context.out / f"seed-{seed}"
    directory.mkdir()
    path = directory / "observables.jsonl"
    error = None
    engine = None
    report = leak = None
    try:
        engine = me.MeltEngine(context.spec, context.catalog, context.temperature,
                               numpy.random.default_rng(seed))
        observer = _LabelledObserver(context.spec, xs=context.xs, identity=engine.identity,
                                     side_group_homolysis=context.side_group_homolysis,
                                     label=dict(context.label, seed=seed))
        series = mo.write_time_series(engine, context.times, path, xs=context.xs, check=True,
                                      observer=observer, max_events=context.max_events)
        stop = series.stop_reason
        report = engine.report()
        estimate = report.pop("leak_estimate")
        for key in ("leak_fraction", "leak_definition", "leak_num", "leak_den"):
            report.pop(key, None)
        leak = {"seed": seed}
        leak.update({k: estimate[k] for k in LEAK_PARTS})
        report, leak = _plain(report), _plain(leak)
    except Exception as exc:   # recorded as the seed's stop; the run's status becomes failed
        stop = "error"
        error = f"{type(exc).__name__}: {exc}"
        report = leak = None
        print(f"seed {seed} failed:", file=sys.stderr, flush=True)
        traceback.print_exc(file=sys.stderr)
        sys.stderr.flush()
    data = path.read_bytes() if path.exists() else b""
    rows = _observation_lines(data.decode())
    result = {"seed": seed, "path": f"seed-{seed}/observables.jsonl", "observables_sha256": _sha256(data),
              "observations": len(rows), "stop_reason": stop,
              "last_observed_time_s": json.loads(rows[-1])["time_s"] if rows else None,
              "engine_time_s": engine.time if engine is not None else None,
              "events": engine.events_fired if engine is not None else 0,
              "null_proposals": engine.null_proposals if engine is not None else 0,
              "report": report, "leak": leak}
    if error is not None:
        result["error"] = error
    result["wall_s"] = time.perf_counter() - start
    return result


_CONTEXT: Optional[_Context] = None


def _install(context: _Context):
    global _CONTEXT
    _CONTEXT = context


def _seed_task(seed: int) -> dict:
    return _run_seed(_CONTEXT, seed)


def _executor(jobs: int, context: _Context):
    """A pool of ``jobs`` forked workers that hold ``context``."""
    import concurrent.futures
    import multiprocessing
    return concurrent.futures.ProcessPoolExecutor(
        max_workers=jobs, mp_context=multiprocessing.get_context("fork"),
        initializer=_install, initargs=(context,))


# --------------------------------------------------------------------------
# side by side
# --------------------------------------------------------------------------
def comparison_rows(observations: Dict[int, List[dict]], times: Sequence[float],
                    mom: Sequence[Optional[dict]], seeds_requested: Optional[int] = None,
                    mom_flags: Sequence[str] = ()) -> List[dict]:
    """Per time and ``side_by_side`` quantity: the seed mean and standard
    error of the kMC value, the MoM value and their difference; a time that
    fewer than ``seeds_requested`` seeds reached gives censored rows.  Every
    row carries the MoM run's flag names (``mom_flags``, ";"-joined, sorted)."""
    requested = len(observations) if seeds_requested is None else seeds_requested
    flags = ";".join(sorted(mom_flags))
    rows = []
    for k, t in enumerate(times):
        seeds = [s for s, obs in observations.items() if len(obs) > k]
        if not seeds:
            continue
        censored = len(seeds) < requested
        slots: Dict[str, dict] = {}
        for s in seeds:
            for r in mo.side_by_side(observations[s][k], mom[k]):
                slot = slots.setdefault(r["quantity"], {"tier": r["tier"], "mom": r["mom"], "values": {}})
                slot["values"][s] = r["kmc"]
        for quantity, slot in slots.items():
            values = []
            for s in seeds:
                v = slot["values"].get(s)
                if v is None and quantity.startswith(mo._ZERO_WHEN_ABSENT):
                    v = 0.0
                if v is not None:
                    values.append(float(v))
            n = len(values)
            mean = math.fsum(values) / n if n else None
            stderr = statistics.stdev(values) / math.sqrt(n) if n >= 2 else None
            m = slot["mom"]
            row = {"time_s": t, "quantity": quantity, "tier": slot["tier"], "seeds_requested": requested,
                   "seeds": n, "censored": censored, "kmc_mean": None, "kmc_stderr": None,
                   "censored_kmc_mean": None, "censored_kmc_stderr": None, "mom": m, "difference": None,
                   "relative_difference": None,
                   "note": "not_a_comparison" if quantity == "convention_gap_g" else "", "mom_flags": flags}
            if censored:
                row["censored_kmc_mean"], row["censored_kmc_stderr"] = mean, stderr
            else:
                difference = mean - m if (mean is not None and m is not None) else None
                row.update({"kmc_mean": mean, "kmc_stderr": stderr, "difference": difference,
                            "relative_difference": (difference / abs(m)
                                                    if (difference is not None and m != 0.0) else None)})
            rows.append(row)
    return rows


def _cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return repr(value)
    return str(value)


def write_comparison(path, rows: List[dict], labels: Dict[str, Any]):
    buffer = io.StringIO()
    head = " ".join(f"{k}={json.dumps(v, sort_keys=True)}" for k, v in labels.items())
    buffer.write(f"# {COMPARISON_SCHEMA} {head}\n")
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(COMPARISON_COLUMNS)
    for row in rows:
        writer.writerow([_cell(row[c]) for c in COMPARISON_COLUMNS])
    Path(path).write_text(buffer.getvalue())


# --------------------------------------------------------------------------
# the run
# --------------------------------------------------------------------------
def sequence_composition(spec: PolymerSpec) -> dict:
    """The declared and realized composition of a copolymer deck's sampled
    melt (exact rationals as strings) and the declared tolerance."""
    from rmgpy.kmc import melt_sequence as msq
    realized = msq.realized_composition(spec)
    declared = dict(spec.initial_ensemble.composition)
    worst = max((abs(realized[r] - x) / x for r, x in declared.items() if x > 0), default=0)
    tolerance = spec.initial_ensemble.composition_tolerance
    return {"sampler": spec.initial_ensemble.sampler, "seed": spec.initial_ensemble.seed,
            "declared": {r: str(x) for r, x in sorted(declared.items())},
            "realized": {r: str(x) for r, x in sorted(realized.items())},
            "max_relative_deviation": float(worst),
            "relative_tolerance": None if tolerance is None else str(tolerance)}


def _write_manifest(out: Path, manifest: dict):
    tmp = out / "run.json.tmp"
    tmp.write_text(json.dumps(manifest, sort_keys=True, indent=1, allow_nan=False) + "\n")
    tmp.replace(out / "run.json")


def run(deck, catalogs: Sequence, temperature: float, times: Sequence[float], seeds: Sequence[int], out, *,
        max_events: Optional[int] = None, allow_incomplete: bool = False, allow_unbound_catalog: bool = False,
        allow_unmatched_catalog: bool = False, mom=None, jobs: int = 1) -> dict:
    """Run every seed and write the outputs (module docstring); returns the
    manifest.  Inputs are refused (``RunError``) before anything is written; a
    seed that fails is recorded and the status is ``failed``."""
    start = time.perf_counter()
    deck = read_deck(deck)
    sequence = isinstance(deck.spec, PolymerSpec)
    if sequence and mom is not None:
        raise RunError(f"{mo.MOM_V2_REFUSAL} (--mom {mom})")
    temperature = float(temperature)
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise RunError(f"the temperature must be finite and positive, got {temperature!r}")
    times = observation_times(times=times)
    seeds = _seeds(seeds)
    if max_events is not None and (isinstance(max_events, bool) or not isinstance(max_events, int)
                                   or max_events < 1):
        raise RunError(f"the event cap must be a positive integer, got {max_events!r}")
    if isinstance(jobs, bool) or not isinstance(jobs, int) or jobs < 1:
        raise RunError(f"jobs must be a positive integer, got {jobs!r}")
    catalog, records = load_catalogs(catalogs, temperature, allow_incomplete=allow_incomplete, deck=deck,
                                     allow_unbound=allow_unbound_catalog)
    mom_by_time, mom_record = (load_mom(mom, times, deck_sha256=deck.sha256, temperature_K=temperature)
                               if mom is not None
                               else ([None] * len(times), None))
    try:   # the declaration as the engine and the observables read it, before any output
        me.TransportLaw.from_spec(deck.spec).check_range(temperature, temperature)
        mo.MeltObserver(deck.spec, xs=deck.xs, side_group_homolysis=deck.side_group_homolysis)
        initial = me.MeltEngine(deck.spec, catalog, temperature, numpy.random.default_rng(0))
        matched = initial_match(initial)
    except ValueError as error:
        raise RunError(f"deck {deck.path}: {error}") from error
    if len(records) > 1:
        owner = {name: i for i, r in enumerate(records) for name in r.template_names}
        overlap = cross_catalog_site(initial, owner)
        if overlap is not None:
            raise RunError(f"a reaction site of the initial state is matched by templates of two catalogs "
                           f"({overlap}); combined, it would fire twice")
    del initial
    if not matched["templates"] and not allow_unmatched_catalog:
        raise RunError("no catalog template matches any site of the initial state (the catalogs may be "
                       "another polymer's); pass --allow-unmatched-catalog to run it, labelled "
                       "catalog_matches_initial_state=false")
    out = Path(out)
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise RunError(f"output directory {out} exists and is not empty")
    out.mkdir(parents=True, exist_ok=True)

    complete = all(r.complete for r in records)
    failures = [{"catalog": i, "path": r.path, "sha256": r.sha256, "failures": list(r.gate_failures)}
                for i, r in enumerate(records) if not r.complete]
    labels = {"catalog_complete": complete, "catalog_gate_failures": failures,
              "catalog_bound": all(r.bound_by for r in records),
              "catalog_matches_initial_state": matched["templates"] > 0}
    context = _Context(deck.spec, catalog, temperature, tuple(times), deck.xs, deck.side_group_homolysis,
                       max_events, out, tuple(sorted(labels.items())))
    manifest = dict({
        "schema": RUN_SCHEMA, "status": "running",
        "deck": deck.manifest(),
        "catalogs": [r.manifest() for r in records],
        "initial_match": matched,
        "temperature_K": temperature, "isothermal": True, "times_s": times, "max_events": max_events,
        "seeds": seeds, "jobs": jobs, "observables_schema": mo.OBSERVABLES_SCHEMA,
        "code": code_identity(),
        "environment": {"python": platform.python_version(), "numpy": numpy.__version__,
                        "PYTHONHASHSEED": os.environ.get("PYTHONHASHSEED")},
        "runs": [],
        "censoring": {"seeds_requested": len(seeds), "stopped_early": [], "censored_times_s": []},
        "leak": {"estimator": me.LEAK_ESTIMATOR + "/pooled", "definition": me.LEAK_DEFINITION,
                 "note": LEAK_NOTE, "per_seed": [], "excluded_seeds": [],
                 "pooled": me.pooled_leak_fraction([])},
        "mom": mom_record, "comparison": None if sequence else "comparison.csv", "wall_s": None,
    }, **labels)
    if sequence:
        manifest["comparison_refused"] = mo.MOM_V2_REFUSAL
        manifest["composition"] = sequence_composition(deck.spec)
    _write_manifest(out, manifest)

    def record(result):
        leak = result.pop("leak")
        last = result.pop("last_observed_time_s")
        manifest["runs"].append(result)
        if leak is None:
            manifest["leak"]["excluded_seeds"].append(result["seed"])
        else:
            manifest["leak"]["per_seed"].append(leak)
        manifest["leak"]["pooled"] = me.pooled_leak_fraction(manifest["leak"]["per_seed"])
        if result["observations"] < len(times):
            manifest["censoring"]["stopped_early"].append(
                {"seed": result["seed"], "stop_reason": result["stop_reason"],
                 "observations": result["observations"], "last_observed_time_s": last,
                 "engine_time_s": result["engine_time_s"]})
            manifest["censoring"]["censored_times_s"] = [
                t for k, t in enumerate(times) if any(r["observations"] <= k for r in manifest["runs"])]
        manifest["wall_s"] = time.perf_counter() - start
        _write_manifest(out, manifest)

    _install(context)
    try:
        if jobs == 1:
            for seed in seeds:
                record(_run_seed(context, seed))
        else:
            with _executor(jobs, context) as pool:
                for result in pool.map(_seed_task, seeds):
                    record(result)
        if not sequence:
            observations = {r["seed"]: [json.loads(x) for x in _observation_lines((out / r["path"]).read_text())]
                            if (out / r["path"]).exists() else [] for r in manifest["runs"]}
            write_comparison(out / "comparison.csv",
                             comparison_rows(observations, times, mom_by_time, len(seeds),
                                             sorted(((mom_record or {}).get("header") or {}).get("flags", {}))),
                             labels)
            from rmgpy.kmc import comparison_rows as cr
            try:
                rows = cr.kmc_rows(out, deck=deck.path)
            except cr.ComparisonRowsError as error:
                manifest["rows_refused"] = str(error)
            else:
                (out / "rows.json").write_text(json.dumps(rows, sort_keys=True, indent=1, allow_nan=False) + "\n")
                manifest["rows"] = "rows.json"
    except BaseException as error:
        manifest["status"] = "failed"
        manifest["error"] = f"{type(error).__name__}: {error}"
        manifest["wall_s"] = time.perf_counter() - start
        _write_manifest(out, manifest)
        raise
    manifest["status"] = "failed" if any(r["stop_reason"] == "error" for r in manifest["runs"]) else "complete"
    manifest["wall_s"] = time.perf_counter() - start
    _write_manifest(out, manifest)
    return json.loads(json.dumps(manifest, sort_keys=True))


# --------------------------------------------------------------------------
# command line
# --------------------------------------------------------------------------
def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m rmgpy.kmc.run",
        description="Isothermal melt kMC of the deck's linear_homopolymer_spec or polymer_spec over emitted catalogs.")
    parser.add_argument("deck", help="RMG input file with one polymer(..., linear_homopolymer_spec=...)")
    parser.add_argument("catalogs", nargs="+",
                        help="kmc-template-catalog/1 files, or kmc-catalog-set/1 manifests (their members)")
    parser.add_argument("--temperature", type=float, required=True, help="K (isothermal)")
    parser.add_argument("--times", type=float, nargs="+", default=None, help="observation times (s)")
    parser.add_argument("--end", type=float, default=None, help="end time (s), with --snapshots")
    parser.add_argument("--snapshots", type=int, default=None, help="evenly spaced observations, 0..end")
    parser.add_argument("--max-events", type=int, default=None, help="event cap per seed")
    parser.add_argument("--seeds", type=int, nargs="+", required=True)
    parser.add_argument("--out", required=True, help="output directory (new or empty)")
    parser.add_argument("--allow-incomplete", action="store_true",
                        help="run catalogs with complete: false, labelled incomplete")
    parser.add_argument("--allow-unbound-catalog", action="store_true",
                        help="run catalogs that carry no declaration digest, labelled catalog_bound=false")
    parser.add_argument("--allow-unmatched-catalog", action="store_true",
                        help="run when no template matches the initial state, labelled")
    parser.add_argument("--mom", default=None, help="MoM reported quantities at the observation times")
    parser.add_argument("--jobs", type=int, default=1, help="worker processes (default 1: sequential)")
    args = parser.parse_args(argv)
    try:
        times = observation_times(args.times, args.end, args.snapshots)
        manifest = run(args.deck, args.catalogs, args.temperature, times, args.seeds, args.out,
                       max_events=args.max_events, allow_incomplete=args.allow_incomplete,
                       allow_unbound_catalog=args.allow_unbound_catalog,
                       allow_unmatched_catalog=args.allow_unmatched_catalog, mom=args.mom, jobs=args.jobs)
    except RunError as error:
        print(f"refused: {error}", file=sys.stderr, flush=True)
        return 2
    for r in manifest["runs"]:
        when = "-" if r["engine_time_s"] is None else f"{r['engine_time_s']:.6g} s"
        print(f"seed {r['seed']}: {r['stop_reason']} at {when}, {r['events']} events, "
              f"{r['observations']} observations, {r['wall_s']:.1f} s wall"
              + (f" ({r['error']})" if "error" in r else ""), flush=True)
    pooled = manifest["leak"]["pooled"]
    print(f"pooled leak fraction {pooled['fraction']:.6g} over {pooled['runs']} seeds; status "
          f"{manifest['status']} (catalog_complete={manifest['catalog_complete']}, "
          f"catalog_bound={manifest['catalog_bound']}, "
          f"catalog_matches_initial_state={manifest['catalog_matches_initial_state']})", flush=True)
    return 1 if manifest["status"] == "failed" else 0


if __name__ == "__main__":
    sys.exit(main())
