"""Run the melt kMC from one RMG deck and emitted template catalogs (I-124).

Command line::

    python -m rmgpy.kmc.run DECK CATALOG [CATALOG ...] --temperature T
        (--times t0 t1 ... | --end T_END --snapshots N) [--max-events N]
        --seeds S [S ...] --out DIR [--allow-incomplete] [--mom FILE] [--jobs N]

**Deck.** An RMG input file.  Only its ``polymer(...)`` call is read: the
deck is executed with every other input-file function (``database``,
``species``, reactors, ``model``, ...) replaced by an inert stub, so nothing
is loaded and no reaction model is set up.  Exactly one ``polymer()`` call is
accepted and it must carry ``linear_homopolymer_spec`` (normalized by
``rmgpy.polymer_spec.normalized_linear_homopolymer_spec``, the input file's
own reader); no call, several calls (with or without a declaration), a label
that differs from the declaration's, or legacy chemistry arguments beside the
declaration are refused.  From the same call the run also reads ``cutoff``
(the MoM hybrid threshold ``xs``, default as in ``polymer()``) and
``side_group_homolysis`` (the MoM side-group channels), so the observables
are partitioned exactly as the MoM run of the same deck.

**Catalogs.** Each file is ``kmc-template-catalog/1`` (``rmgpy.kmc.emit``)
with ``provenance``, ``payloads``, ``gates.failures`` and the ``complete``
flag; a flag that contradicts the gate failures is refused.  A catalog with
``complete: false`` is refused unless ``--allow-incomplete``; then every
output file carries ``catalog_complete: false`` and the gate failures
(``catalog_gate_failures``).  Catalogs whose loaded template names collide are
refused.  Each catalog must load at the run temperature (no extrapolation).

**Run.** Isothermal.  Per seed ``n`` (sequentially, in one process, unless
``--jobs``): ``MeltEngine(spec, catalog, T, numpy.random.default_rng(n))`` and
``melt_observables.write_time_series`` with the invariant checker on, into
``<out>/seed-<n>/observables.jsonl``; each line is then labelled with
``seed``, ``catalog_complete`` and ``catalog_gate_failures``.  The optional
event cap is run-wide per seed: when it is reached the seed stops, the
observation in progress is not written (``observations`` counts the lines),
and the stop reason is ``max_events``; otherwise the stop reason is that of
the last horizon (``horizon`` or ``absorbing``).

**Manifest** ``<out>/run.json`` (``kmc-melt-run/1``)::

    schema, status ("running" | "complete" | "failed"[, error])
    deck: {path, sha256, label, declaration (normalized mapping),
           declaration_sha256, xs, side_group_homolysis}
    catalogs: [{path, sha256, schema, complete, gate_failures, provenance,
                templates (count at T), template_names_sha256}]
    catalog_complete, catalog_gate_failures: [{catalog, path, sha256, failures}]
    temperature_K, isothermal, times_s, max_events, seeds, jobs, observables_schema
    code: {sha (git HEAD of the checkout), dirty, dirty_paths}
    environment: {python, numpy, PYTHONHASHSEED}
    runs: [{seed, path, observables_sha256, observations, stop_reason,
            engine_time_s, events, null_proposals, wall_s,
            report (engine report without its leak entries)}]
    leak: {estimator, definition, per_seed: [{seed, leak_num, leak_den,
           non_thinned_integral, phi_sum, thinned_proposals, open_thinned,
           open_phi_source}], pooled (melt_engine.pooled_leak_fraction)}
    mom: null | {path, sha256, times}
    comparison: "comparison.csv"
    wall_s

The run's leak number is ``leak.pooled.fraction`` only: per seed the
manifest carries the additive parts, never a single-run fraction.  Rerunning
with the same inputs gives byte-identical observation files and a manifest
that differs only in the ``wall_s`` entries.

**Comparison** ``<out>/comparison.csv`` (``kmc-melt-comparison/1``; a first
``#`` line carries the schema and the catalog label): one row per
observation time and ``melt_observables.side_by_side`` quantity, with
``seeds`` (the seeds that reached the time), ``kmc_mean``, ``kmc_stderr``
(sample standard deviation / sqrt(seeds); empty for one seed), the MoM value
and ``difference = kmc_mean - mom`` (``relative_difference`` = difference /
|mom|).  ``explicit_dp_mol.*`` and ``species_mol.*`` absent from a seed count
as 0 (as in ``side_by_side``).  Without ``--mom`` the MoM columns are empty.
``--mom FILE``: JSON (a list) or JSON lines of ``{"time_s": t, "quantities":
{...}}`` with ``quantities`` keyed like ``mom_reported_quantities`` (the
format ``side_by_side`` takes); every MoM time must be an observation time.

Nothing here names a polymer: the chemistry is the deck's declaration and the
catalogs' templates.
"""

from __future__ import annotations

import argparse
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
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy

from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import melt_observables as mo
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
# the names rmgpy.rmg.input.read_input_file exposes to a deck (besides True / False)
DECK_VOCABULARY = (
    "database", "catalystProperties", "species", "polymer", "forbidden", "SMARTS", "fragment_adj",
    "fragment_SMILES", "SMILES", "InChI", "adjacencyList", "adjacencyListGroup", "react",
    "simpleReactor", "constantVIdealGasReactor", "constantTPIdealGasReactor", "liquidSurfaceReactor",
    "constantTVLiquidReactor", "liquidReactor", "surfaceReactor", "mbsampledReactor",
    "hybridPolymerReactor", "polymer_phase", "PolymerPool", "MassTransfer", "simulator", "solvation",
    "SolventData", "liquidVolumetricMassTransferCoefficientPowerLaw", "model", "quantumMechanics",
    "mlEstimator", "pressureDependence", "options", "generatedSpeciesConstraints",
    "generatePolymerConstraints", "thermoCentralDatabase", "uncertainty", "restartFromSeed")
LEAK_NOTE = ("the run's leak number is pooled.fraction (sum of leak_num over seeds / sum of leak_den); "
             "per seed only the additive parts are reported")
LEAK_PARTS = ("leak_num", "leak_den", "non_thinned_integral", "phi_sum", "thinned_proposals",
              "open_thinned", "open_phi_source")
COMPARISON_COLUMNS = ("time_s", "quantity", "tier", "seeds", "kmc_mean", "kmc_stderr", "mom",
                      "difference", "relative_difference", "note")
MOM_TIME_RELATIVE = 1.0e-9


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


# --------------------------------------------------------------------------
# the deck
# --------------------------------------------------------------------------
class _Inert:
    """Stands for every input-file function except ``polymer``: accepts any
    call, attribute or index and does nothing."""

    __slots__ = ()

    def __call__(self, *args, **kwargs):
        return self

    def __getattr__(self, name):
        return self

    def __getitem__(self, key):
        return self

    def __repr__(self):
        return "<inert>"


_INERT = _Inert()


@dataclass(frozen=True)
class Deck:
    path: str
    sha256: str
    spec: LinearHomopolymerSpec
    label: str
    xs: int
    side_group_homolysis: Optional[list]

    def manifest(self) -> dict:
        return {"path": self.path, "sha256": self.sha256, "label": self.label,
                "declaration": json.loads(self.spec.normalized_json),
                "declaration_sha256": _sha256(self.spec.normalized_json),
                "xs": self.xs, "side_group_homolysis": self.side_group_homolysis}


def _bind_polymer(args, kwargs) -> dict:
    if len(args) > len(POLYMER_PARAMETERS):
        raise TypeError(f"polymer() takes at most {len(POLYMER_PARAMETERS)} positional arguments")
    bound = dict(zip(POLYMER_PARAMETERS, args))
    for key, value in kwargs.items():
        if key not in POLYMER_PARAMETERS:
            raise TypeError(f"polymer() got an unexpected keyword argument {key!r}")
        if key in bound:
            raise TypeError(f"polymer() got multiple values for argument {key!r}")
        bound[key] = value
    if "label" not in bound:
        raise TypeError("polymer() missing required argument 'label'")
    return bound


def read_deck(path) -> Deck:
    """The one ``linear_homopolymer_spec`` declaration of an RMG deck (see the
    module docstring); nothing else of the deck is evaluated for effect."""
    path = Path(path)
    try:
        raw = path.read_bytes()
    except OSError as error:
        raise RunError(f"deck {path}: {error}") from error
    calls: List[dict] = []

    def polymer(*args, **kwargs):
        calls.append(_bind_polymer(args, kwargs))

    local_context: Dict[str, Any] = {name: _INERT for name in DECK_VOCABULARY}
    local_context.update({"__builtins__": None, "True": True, "False": False, "polymer": polymer})
    try:
        exec(compile(raw.decode("utf-8"), str(path), "exec"), {"__builtins__": None}, local_context)
    except Exception as error:   # the deck is user input: any failure is a refusal
        raise RunError(f"deck {path} could not be read: {type(error).__name__}: {error}") from error
    if not calls:
        raise RunError(f"deck {path} declares no polymer(); the melt kMC needs one polymer() call "
                       "carrying linear_homopolymer_spec")
    if len(calls) > 1:
        labels = [c["label"] for c in calls]
        raise RunError(f"deck {path} declares {len(calls)} polymer() calls ({labels}); the melt kMC "
                       "runs exactly one linear_homopolymer_spec declaration")
    call = calls[0]
    value = call.get("linear_homopolymer_spec")
    if value is None:
        raise RunError(f"deck {path}: polymer({call['label']!r}) carries no linear_homopolymer_spec")
    legacy = [k for k in LEGACY_CHEMISTRY if call.get(k) is not None]
    if legacy:
        raise RunError(f"deck {path}: polymer({call['label']!r}) combines linear_homopolymer_spec "
                       f"with legacy arguments {legacy}")
    try:
        spec = normalized_linear_homopolymer_spec(value)
    except (TypeError, ValueError, KeyError) as error:
        raise RunError(f"deck {path}: linear_homopolymer_spec refused: {error}") from error
    if call["label"] != spec.label:
        raise RunError(f"deck {path}: polymer label {call['label']!r} differs from the declaration "
                       f"label {spec.label!r}")
    xs = call.get("cutoff", POLYMER_CUTOFF_DEFAULT)
    if isinstance(xs, bool) or not isinstance(xs, int) or xs < 0:
        raise RunError(f"deck {path}: polymer cutoff must be a non-negative integer, got {xs!r}")
    sgh = call.get("side_group_homolysis")
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

    def manifest(self) -> dict:
        return {"path": self.path, "sha256": self.sha256, "schema": CATALOG_SCHEMA,
                "complete": self.complete, "gate_failures": list(self.gate_failures),
                "provenance": self.provenance, "templates": len(self.template_names),
                "template_names_sha256": _sha256("\n".join(self.template_names).encode())}


def load_catalogs(paths: Sequence, temperature: float, *, allow_incomplete: bool = False
                  ) -> Tuple[me.PayloadCatalog, List[CatalogRecord]]:
    """The union of the catalogs as one ``PayloadCatalog`` and one record per
    file.  Refused: a wrong schema or missing field, a ``complete`` flag that
    contradicts the gate failures, an incomplete catalog without
    ``allow_incomplete``, a catalog that does not load at ``temperature`` or
    loads no template, template names that collide between catalogs."""
    if not paths:
        raise RunError("at least one catalog file is needed")
    temperature = float(temperature)
    records: List[CatalogRecord] = []
    payloads: List[dict] = []
    owner: Dict[str, str] = {}
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
        if not isinstance(doc["payloads"], list):
            raise RunError(f"catalog {path}: payloads must be a list")
        try:
            names = [t.name for t in me.PayloadCatalog(doc["payloads"]).templates(temperature)]
        except Exception as error:   # EmissionError, malformed payloads
            raise RunError(f"catalog {path} does not load at {temperature} K: {error}") from error
        if not names:
            raise RunError(f"catalog {path} loads no template at {temperature} K")
        for name in names:
            if name in owner:
                raise RunError(f"catalogs {owner[name]} and {path} collide on template name {name!r}")
            owner[name] = str(path)
        payloads += doc["payloads"]
        records.append(CatalogRecord(str(path), _sha256(raw), _plain(doc["provenance"]), complete,
                                     tuple(failures), tuple(names)))
    return me.PayloadCatalog(payloads), records


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


def load_mom(path, times: Sequence[float]) -> Tuple[List[Optional[dict]], dict]:
    """Per observation time the MoM quantities (None where the file has none)."""
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
        match = [k for k, ti in enumerate(times)
                 if t == ti or abs(t - ti) <= MOM_TIME_RELATIVE * max(abs(t), abs(ti))]
        if not match:
            raise RunError(f"MoM file {path}: time {t} s is not an observation time of the run {list(times)}")
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
    engine = _SeedEngine(context.spec, context.catalog, context.temperature, numpy.random.default_rng(seed),
                         event_cap=context.max_events)
    try:
        mo.write_time_series(engine, context.times, path, xs=context.xs, check=True,
                             side_group_homolysis=context.side_group_homolysis)
        stop = engine.stop_reasons[-1] if engine.stop_reasons else "horizon"
    except _EventCapReached:
        stop = "max_events"
    label = dict(context.label, seed=seed)
    labelled = []
    for line in path.read_text().splitlines():
        row = json.loads(line)
        row.update(label)
        labelled.append(json.dumps(row, sort_keys=True, allow_nan=False) + "\n")
    data = "".join(labelled).encode()
    path.write_bytes(data)
    report = engine.report()
    estimate = report.pop("leak_estimate")
    for key in ("leak_fraction", "leak_definition", "leak_num", "leak_den"):
        report.pop(key, None)
    leak = {"seed": seed}
    leak.update({k: estimate[k] for k in LEAK_PARTS})
    return {"seed": seed, "path": f"seed-{seed}/observables.jsonl", "observables_sha256": _sha256(data),
            "observations": len(labelled), "stop_reason": stop, "engine_time_s": engine.time,
            "events": engine.events_fired, "null_proposals": engine.null_proposals,
            "report": _plain(report), "leak": _plain(leak), "wall_s": time.perf_counter() - start}


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
                    mom: Sequence[Optional[dict]]) -> List[dict]:
    """Per time and ``side_by_side`` quantity: the seed mean and standard
    error of the kMC value, the MoM value and their difference."""
    rows = []
    for k, t in enumerate(times):
        seeds = [s for s, obs in observations.items() if len(obs) > k]
        if not seeds:
            continue
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
            difference = mean - m if (mean is not None and m is not None) else None
            relative = difference / abs(m) if (difference is not None and m != 0.0) else None
            rows.append({"time_s": t, "quantity": quantity, "tier": slot["tier"], "seeds": n,
                         "kmc_mean": mean, "kmc_stderr": stderr, "mom": m, "difference": difference,
                         "relative_difference": relative,
                         "note": "not_a_comparison" if quantity == "convention_gap_g" else ""})
    return rows


def _cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return repr(value)
    return str(value)


def write_comparison(path, rows: List[dict], complete: bool, failures: list):
    buffer = io.StringIO()
    buffer.write(f"# {COMPARISON_SCHEMA} catalog_complete={'true' if complete else 'false'} "
                 f"catalog_gate_failures={json.dumps(failures, sort_keys=True)}\n")
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
        max_events: Optional[int] = None, allow_incomplete: bool = False, mom=None, jobs: int = 1) -> dict:
    """Run every seed and write the outputs (module docstring); returns the manifest."""
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
    catalog, records = load_catalogs(catalogs, temperature, allow_incomplete=allow_incomplete)
    mom_by_time, mom_record = load_mom(mom, times) if mom is not None else ([None] * len(times), None)
    try:   # the declaration as the engine and the observables read it, before any output
        me.TransportLaw.from_spec(deck.spec).check_range(temperature, temperature)
        mo.MeltObserver(deck.spec, xs=deck.xs, side_group_homolysis=deck.side_group_homolysis)
    except ValueError as error:
        raise RunError(f"deck {deck.path}: {error}") from error
    out = Path(out)
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise RunError(f"output directory {out} exists and is not empty")
    out.mkdir(parents=True, exist_ok=True)

    complete = all(r.complete for r in records)
    failures = [{"catalog": i, "path": r.path, "sha256": r.sha256, "failures": list(r.gate_failures)}
                for i, r in enumerate(records) if not r.complete]
    context = _Context(deck.spec, catalog, temperature, tuple(times), deck.xs, deck.side_group_homolysis,
                       max_events, out, (("catalog_complete", complete), ("catalog_gate_failures", failures)))
    manifest = {
        "schema": RUN_SCHEMA, "status": "running",
        "deck": deck.manifest(),
        "catalogs": [r.manifest() for r in records],
        "catalog_complete": complete, "catalog_gate_failures": failures,
        "temperature_K": temperature, "isothermal": True, "times_s": times, "max_events": max_events,
        "seeds": seeds, "jobs": jobs, "observables_schema": mo.OBSERVABLES_SCHEMA,
        "code": code_identity(),
        "environment": {"python": platform.python_version(), "numpy": numpy.__version__,
                        "PYTHONHASHSEED": os.environ.get("PYTHONHASHSEED")},
        "runs": [],
        "leak": {"estimator": me.LEAK_ESTIMATOR + "/pooled", "definition": me.LEAK_DEFINITION,
                 "note": LEAK_NOTE, "per_seed": [], "pooled": me.pooled_leak_fraction([])},
        "mom": mom_record, "comparison": "comparison.csv", "wall_s": None,
    }
    _write_manifest(out, manifest)

    def record(result):
        leak = result.pop("leak")
        manifest["runs"].append(result)
        manifest["leak"]["per_seed"].append(leak)
        manifest["leak"]["pooled"] = me.pooled_leak_fraction(manifest["leak"]["per_seed"])
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
                        for r in manifest["runs"]}
        write_comparison(out / "comparison.csv", comparison_rows(observations, times, mom_by_time),
                         complete, failures)
    except BaseException as error:
        manifest["status"] = "failed"
        manifest["error"] = f"{type(error).__name__}: {error}"
        manifest["wall_s"] = time.perf_counter() - start
        _write_manifest(out, manifest)
        raise
    manifest["status"] = "complete"
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
    parser.add_argument("--mom", default=None, help="MoM reported quantities at the observation times")
    parser.add_argument("--jobs", type=int, default=1, help="worker processes (default 1: sequential)")
    args = parser.parse_args(argv)
    try:
        times = observation_times(args.times, args.end, args.snapshots)
        manifest = run(args.deck, args.catalogs, args.temperature, times, args.seeds, args.out,
                       max_events=args.max_events, allow_incomplete=args.allow_incomplete, mom=args.mom,
                       jobs=args.jobs)
    except RunError as error:
        print(f"refused: {error}", file=sys.stderr, flush=True)
        return 2
    for r in manifest["runs"]:
        print(f"seed {r['seed']}: {r['stop_reason']} at {r['engine_time_s']:.6g} s, {r['events']} events, "
              f"{r['observations']} observations, {r['wall_s']:.1f} s wall", flush=True)
    pooled = manifest["leak"]["pooled"]
    print(f"pooled leak fraction {pooled['fraction']:.6g} over {pooled['runs']} seeds "
          f"(catalog_complete={manifest['catalog_complete']})", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
