"""Run the melt kMC from one RMG deck and emitted template catalogs (I-124).

Command line::

    python -m rmgpy.kmc.run DECK CATALOG [CATALOG ...] --temperature T
        (--times t0 t1 ... | --end T_END --snapshots N) [--max-events N]
        --seeds S [S ...] --out DIR [--allow-incomplete] [--allow-unbound-catalog]
        [--allow-unmatched-catalog] [--mom FILE] [--jobs N]

**Deck.** An RMG input file, read without executing it: it is parsed with
``ast`` and must contain exactly one ``polymer(...)`` call, as a module-level
statement; the name ``polymer`` may appear nowhere else (no alias, no call
inside a conditional, loop, function or expression).  Its arguments bind like
``rmgpy.rmg.input.polymer``'s (``POLYMER_PARAMETERS``), ``*args`` / ``**kw``
are refused, and the arguments read are ``ast.literal_eval``-ed: a non-literal
one is refused.  ``linear_homopolymer_spec`` (a literal mapping or JSON
string) is normalized by ``rmgpy.polymer_spec.normalized_linear_homopolymer_spec``
(the input file's own reader); a missing declaration, a label that differs
from the declaration's, or legacy chemistry arguments beside it are refused.
From the same call the run reads ``cutoff`` (the MoM hybrid threshold ``xs``,
default as in ``polymer()``) and ``side_group_homolysis`` (the MoM side-group
channels), so the observables are partitioned exactly as the MoM run of the
same deck.  The rest of the deck is never evaluated.

**Declaration digests.**  ``declaration_sha256`` = sha256 of
``LinearHomopolymerSpec.normalized_json`` (the bytes of the normalized
declaration); ``chemistry_sha256`` = sha256 of the compact sorted-key JSON
(``separators=(",", ":")``) of the normalized declaration restricted to
``CHEMISTRY_FIELDS`` (schema_version, repeat_unit, end_groups,
locality_policy, family_manifest, library_manifest): the chemistry a catalog
is emitted for, independent of the ensemble, mass, density, transport and
declared species of a particular run.  Both are written to ``run.json``
(``deck``) and are the values an emitter writes into a catalog's
``provenance`` under the same names.

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

Always refused: loaded template names that collide, and an executable
(``sha256_json`` of the template's ``executable``) that two emitted templates
share, within a catalog or across catalogs (it would fire twice).  Each
catalog must load at the run temperature (no extrapolation).

**Run.** Isothermal.  Per seed ``n`` (sequentially, in one process, unless
``--jobs``): ``MeltEngine(spec, catalog, T, numpy.random.default_rng(n))`` and
``melt_observables.write_time_series`` with the invariant checker on, into
``<out>/seed-<n>/observables.jsonl``; every line carries, as it is written,
``seed`` and the labels above.  The optional event cap is run-wide per seed:
when it is reached the seed stops, the observation in progress is not written
(``observations`` counts the lines), and the stop reason is ``max_events``;
an exception in a seed stops it with stop reason ``error`` (``error`` names
it; the traceback goes to stderr), the other seeds still run and the run's
status is ``failed``; otherwise the stop reason is that of the last horizon
(``horizon`` or ``absorbing``).

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
            report (engine report without its leak entries; null for a failed seed)}]
    censoring: {seeds_requested, censored_times_s,
                stopped_early: [{seed, stop_reason, observations,
                                 last_observed_time_s, engine_time_s}]}
    leak: {estimator, definition, note, excluded_seeds (seeds that failed),
           per_seed: [{seed, leak_num, leak_den, non_thinned_integral, phi_sum,
                       thinned_proposals, open_thinned, open_phi_source}],
           pooled (melt_engine.pooled_leak_fraction over per_seed)}
    mom: null | {path, sha256, times}
    comparison: "comparison.csv"
    wall_s

The run's leak number is ``leak.pooled.fraction`` only: per seed the
manifest carries the additive parts, never a single-run fraction.  Rerunning
with the same inputs gives byte-identical observation files and a manifest
that differs only in the ``wall_s`` entries.

**Comparison** ``<out>/comparison.csv`` (``kmc-melt-comparison/1``; a first
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
Without ``--mom`` the MoM columns are empty.

``--mom FILE``: JSON (a list) or JSON lines of ``{"time_s": t, "quantities":
{...}}`` with ``quantities`` keyed like ``mom_reported_quantities`` (the
format ``side_by_side`` takes), every number finite.  A MoM time matches the
observation time within relative ``MOM_TIME_RELATIVE`` (1e-12, float
round-off only) of it; a finite time matching none, or more than one, is
refused, and so is a non-finite one.

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
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy

from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import melt_observables as mo
from rmgpy.kmc.template_emission import sha256_json
from rmgpy.polymer_spec import LinearHomopolymerSpec, normalized_linear_homopolymer_spec

RUN_SCHEMA = "kmc-melt-run/1"
COMPARISON_SCHEMA = "kmc-melt-comparison/1"
CATALOG_SCHEMA = "kmc-template-catalog/1"     # rmgpy.kmc.emit.SCHEMA
CATALOG_KEYS = ("schema", "provenance", "payloads", "gates", "complete")
# rmgpy.rmg.input.polymer's parameters, in order (positional binding), and its cutoff default
POLYMER_PARAMETERS = (
    "label", "monomer", "end_groups", "cutoff", "Mn", "Mw", "moments", "initial_mass", "k_unzip",
    "k_scission", "discrete_dp_threshold", "monomer_product", "radical_qssa_unzip", "explicit_dp",
    "k_homolysis", "k_depropagation", "side_group_homolysis", "thermal_analysis_inputs",
    "linear_homopolymer_spec")
POLYMER_CUTOFF_DEFAULT = 4
LEGACY_CHEMISTRY = ("monomer", "end_groups", "moments", "Mn", "Mw", "monomer_product")
# the declaration fields that fix the chemistry a catalog was emitted for (chemistry_sha256)
CHEMISTRY_FIELDS = ("schema_version", "repeat_unit", "end_groups", "locality_policy", "family_manifest",
                    "library_manifest")
BINDING_KEYS = ("chemistry_sha256", "declaration_sha256")
LEAK_NOTE = ("the run's leak number is pooled.fraction (sum of leak_num over seeds / sum of leak_den); "
             "per seed only the additive parts are reported; seeds that failed are excluded")
LEAK_PARTS = ("leak_num", "leak_den", "non_thinned_integral", "phi_sum", "thinned_proposals",
              "open_thinned", "open_phi_source")
COMPARISON_COLUMNS = ("time_s", "quantity", "tier", "seeds_requested", "seeds", "censored", "kmc_mean",
                      "kmc_stderr", "censored_kmc_mean", "censored_kmc_stderr", "mom", "difference",
                      "relative_difference", "note")
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


def declaration_sha256(spec: LinearHomopolymerSpec) -> str:
    return _sha256(spec.normalized_json)


def chemistry_sha256(spec: LinearHomopolymerSpec) -> str:
    normalized = json.loads(spec.normalized_json)
    return _sha256(json.dumps({k: normalized[k] for k in CHEMISTRY_FIELDS}, sort_keys=True,
                              separators=(",", ":")).encode())


# --------------------------------------------------------------------------
# the deck
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Deck:
    path: str
    sha256: str
    spec: LinearHomopolymerSpec
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
    value = literal("linear_homopolymer_spec")
    if value is None:
        raise RunError(f"deck {path}: polymer({label!r}) carries no linear_homopolymer_spec")
    legacy = [k for k in LEGACY_CHEMISTRY if k in nodes and not _is_none(nodes[k])]
    if legacy:
        raise RunError(f"deck {path}: polymer({label!r}) combines linear_homopolymer_spec "
                       f"with legacy arguments {legacy}")
    try:
        spec = normalized_linear_homopolymer_spec(value)
    except (TypeError, ValueError, KeyError) as error:
        raise RunError(f"deck {path}: linear_homopolymer_spec refused: {error}") from error
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


def load_catalogs(paths: Sequence, temperature: float, *, allow_incomplete: bool = False,
                  deck: Optional[Deck] = None, allow_unbound: bool = False
                  ) -> Tuple[me.PayloadCatalog, List[CatalogRecord]]:
    """The union of the catalogs as one ``PayloadCatalog`` and one record per
    file.  Refused: a wrong schema or missing field, a ``complete`` flag that
    contradicts the gate failures, an incomplete catalog without
    ``allow_incomplete``, a catalog that does not load at ``temperature`` or
    loads no template, template names or executables shared by two emitted
    templates; with ``deck``, a declaration digest that differs from the
    deck's, and a catalog without one unless ``allow_unbound``."""
    if not paths:
        raise RunError("at least one catalog file is needed")
    temperature = float(temperature)
    records: List[CatalogRecord] = []
    payloads: List[dict] = []
    owner: Dict[str, str] = {}
    executables: Dict[str, str] = {}
    for path in paths:
        path = Path(path)
        try:
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
            names = [t.name for t in me.PayloadCatalog(doc["payloads"]).templates(temperature)]
            emitted = [(p.get("record_event_id"), t.get("key"), sha256_json(t["executable"]))
                       for p in doc["payloads"] for t in p.get("templates", ())
                       if t.get("status") == "emitted"]
        except Exception as error:   # EmissionError, malformed payloads
            raise RunError(f"catalog {path} does not load at {temperature} K: {error}") from error
        if not names:
            raise RunError(f"catalog {path} loads no template at {temperature} K")
        for name in names:
            if name in owner:
                raise RunError(f"catalogs {owner[name]} and {path} collide on template name {name!r}")
            owner[name] = str(path)
        for record, key, digest in emitted:
            where = f"{path} (record {str(record)[:16]}, key {str(key)[:12]})"
            if digest in executables:
                raise RunError(f"two emitted templates share the executable {digest[:16]}: "
                               f"{executables[digest]} and {where}; loading both would fire it twice")
            executables[digest] = where
        payloads += doc["payloads"]
        records.append(CatalogRecord(str(path), _sha256(raw), _plain(doc["provenance"]), complete,
                                     tuple(failures), tuple(names), bound_by))
    return me.PayloadCatalog(payloads), records


def initial_match(engine: me.MeltEngine) -> dict:
    """The templates (classes) matching some site of the engine's current
    state: a unimolecular transition (allowed or leak) or a bimolecular
    participant site; ``sites`` counts both."""
    classes, sites = set(), 0
    for entry in engine.entries.values():
        for ci, (allowed, leak) in entry.uni.items():
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


def _finite_numbers(value) -> bool:
    if isinstance(value, bool) or value is None or isinstance(value, str):
        return True
    if isinstance(value, (int, float)):
        return math.isfinite(value)
    if isinstance(value, dict):
        return all(_finite_numbers(v) for v in value.values())
    if isinstance(value, list):
        return all(_finite_numbers(v) for v in value)
    return True


def load_mom(path, times: Sequence[float]) -> Tuple[List[Optional[dict]], dict]:
    """Per observation time the MoM quantities (None where the file has none);
    matching in the module docstring."""
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
    if isinstance(doc, dict):
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
        if not _finite_numbers(entry["quantities"]):
            raise RunError(f"MoM file {path}: a quantity at time {t} s is not finite")
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
    return out, {"path": str(path), "sha256": _sha256(raw), "times": sum(q is not None for q in out)}


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
class _EventCapReached(Exception):
    pass


class _SeedEngine(me.MeltEngine):
    """The engine with a run-wide event cap across ``write_time_series``'
    horizon calls; records each call's stop reason.  Chemistry, clock and
    index are the engine's own."""

    def __init__(self, *args, event_cap: Optional[int] = None, **kwargs):
        super().__init__(*args, **kwargs)
        self.event_cap = event_cap
        self.stop_reasons: List[str] = []

    def run(self, *, until=None, max_events=None, check=None):
        if self.event_cap is not None:
            left = max(0, self.event_cap - self.events_fired)
            max_events = left if max_events is None else min(max_events, left)
        result = super().run(until=until, max_events=max_events, check=check)
        self.stop_reasons.append(result.stop_reason)
        if (result.stop_reason == "max_events" and self.event_cap is not None
                and self.events_fired >= self.event_cap):
            raise _EventCapReached()
        return result


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
    spec: LinearHomopolymerSpec
    catalog: me.PayloadCatalog
    temperature: float
    times: Tuple[float, ...]
    xs: int
    side_group_homolysis: Optional[list]
    max_events: Optional[int]
    out: Path
    label: Tuple[Tuple[str, Any], ...]


def _run_seed(context: _Context, seed: int) -> dict:
    start = time.perf_counter()
    directory = context.out / f"seed-{seed}"
    directory.mkdir()
    path = directory / "observables.jsonl"
    error = None
    engine = None
    try:
        engine = _SeedEngine(context.spec, context.catalog, context.temperature, numpy.random.default_rng(seed),
                             event_cap=context.max_events)
        observer = _LabelledObserver(context.spec, xs=context.xs, identity=engine.identity,
                                     side_group_homolysis=context.side_group_homolysis,
                                     label=dict(context.label, seed=seed))
        mo.write_time_series(engine, context.times, path, xs=context.xs, check=True, observer=observer)
        stop = engine.stop_reasons[-1] if engine.stop_reasons else "horizon"
    except _EventCapReached:
        stop = "max_events"
    except Exception as exc:   # recorded as the seed's stop; the run's status becomes failed
        stop = "error"
        error = f"{type(exc).__name__}: {exc}"
        print(f"seed {seed} failed:", file=sys.stderr, flush=True)
        traceback.print_exc(file=sys.stderr)
        sys.stderr.flush()
    data = path.read_bytes() if path.exists() else b""
    rows = data.decode().splitlines()
    result = {"seed": seed, "path": f"seed-{seed}/observables.jsonl", "observables_sha256": _sha256(data),
              "observations": len(rows), "stop_reason": stop,
              "last_observed_time_s": json.loads(rows[-1])["time_s"] if rows else None,
              "engine_time_s": engine.time if engine is not None else None,
              "events": engine.events_fired if engine is not None else 0,
              "null_proposals": engine.null_proposals if engine is not None else 0,
              "report": None, "leak": None}
    if error is not None:
        result["error"] = error
    if engine is not None and error is None:
        report = engine.report()
        estimate = report.pop("leak_estimate")
        for key in ("leak_fraction", "leak_definition", "leak_num", "leak_den"):
            report.pop(key, None)
        leak = {"seed": seed}
        leak.update({k: estimate[k] for k in LEAK_PARTS})
        result["report"], result["leak"] = _plain(report), _plain(leak)
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
                    mom: Sequence[Optional[dict]], seeds_requested: Optional[int] = None) -> List[dict]:
    """Per time and ``side_by_side`` quantity: the seed mean and standard
    error of the kMC value, the MoM value and their difference; a time that
    fewer than ``seeds_requested`` seeds reached gives censored rows."""
    requested = len(observations) if seeds_requested is None else seeds_requested
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
                   "note": "not_a_comparison" if quantity == "convention_gap_g" else ""}
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
    mom_by_time, mom_record = load_mom(mom, times) if mom is not None else ([None] * len(times), None)
    try:   # the declaration as the engine and the observables read it, before any output
        me.TransportLaw.from_spec(deck.spec).check_range(temperature, temperature)
        mo.MeltObserver(deck.spec, xs=deck.xs, side_group_homolysis=deck.side_group_homolysis)
        matched = initial_match(me.MeltEngine(deck.spec, catalog, temperature, numpy.random.default_rng(0)))
    except ValueError as error:
        raise RunError(f"deck {deck.path}: {error}") from error
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
        "mom": mom_record, "comparison": "comparison.csv", "wall_s": None,
    }, **labels)
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
        observations = {r["seed"]: [json.loads(x) for x in (out / r["path"]).read_text().splitlines()]
                        if (out / r["path"]).exists() else [] for r in manifest["runs"]}
        write_comparison(out / "comparison.csv",
                         comparison_rows(observations, times, mom_by_time, len(seeds)), labels)
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
        description="Isothermal melt kMC of the deck's linear_homopolymer_spec over emitted catalogs.")
    parser.add_argument("deck", help="RMG input file with one polymer(..., linear_homopolymer_spec=...)")
    parser.add_argument("catalogs", nargs="+", help="kmc-template-catalog/1 files")
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
