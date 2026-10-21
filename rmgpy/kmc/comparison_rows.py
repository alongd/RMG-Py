"""Comparison rows (I-217, slice S7 of the shared-chemistry design): the
exact inventory rows of the method comparison, one schema for every
producer -- a kMC run (``run``), a kernel reference (``kernel_reference``:
the population balance or the kernel SSA) and a MoM ``--mom`` file -- and
the join of two row files.

**A row document** (``SCHEMA``) holds one pass of one producer:
``producer`` (``kind``: ``kmc``, ``kernel_pbe``, ``kernel_ssa`` or ``mom``,
and its source), ``binding`` (the deck's sha256 and chemistry digest, the
temperature, the observation times, the pass: ``declared`` or
``well_mixed_bound``), ``geometry`` (periods per declared unit, backbone
atoms per period, the ladder cutoff in periods), ``withheld_families``
(a row family the producer cannot give at all, with the reason) and
``rows``.  A row is ``{time_s, row, producer, kind, value, stderr, lower,
upper, seeds, basis, reason}`` with ``kind`` one of ``KINDS``:

* ``exact``: ``value`` (a mean over seeds with its standard error where the
  producer is stochastic);
* ``interval``: ``[lower, upper]``, the value is not determined (``value``
  None);
* ``withheld``: the producer cannot give the row here (``reason``);
* ``diagnostic``: a number about the producer itself, never compared.

**The rows** (``DEFINITIONS``): element inventory and its residual against
t = 0; total mass; each declared species; the mass outside them; the
explicit ladder under the cutoff and the length distribution above it, in
primitive periods; the whole-ensemble moments.  A molecule of B backbone
atoms has the length ``[floor(B / nbb), ceil(B / nbb)]`` periods
(``length_interval``); a bin's amount is the interval from the molecules
surely in it to those possibly in it, exact when every molecule is on the
period grid.  The kMC run's own convention (``n_c = floor(B / b)``) lies in
that interval.

**The join** (``join``) refuses two documents whose deck, temperature, time
grid or pass differ, naming the field; per time and row it gives the
difference (and its standard error) only where both sides are exact, and
otherwise says which side withholds (and why), which side lacks the row,
or, for an interval, whether the other side's value lies in it.
"""

import argparse
import csv
import hashlib
import io
import json
import math
import statistics
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

SCHEMA = "kmc-comparison-rows/1"
JOIN_SCHEMA = "kmc-comparison-join/1"
KINDS = ("exact", "interval", "withheld", "diagnostic")
PASS_DECLARED = "declared"
PASS_WELL_MIXED = "well_mixed_bound"
DEFINITIONS = {
    "element_mol.<El>": "atoms of the element in every molecule of the system (mol)",
    "element_residual_mol.<El>": "element_mol at t minus element_mol at t = 0 of the same run (mol); 0 for a "
                                 "conserving producer",
    "element_drift_mol.<El>": "diagnostic: the population balance's integrated inventory minus its initial one "
                              "(round-off of the integrator; the balance's inventory is the per-term invariant)",
    "total_mass_g": "mass of every molecule (g)",
    "species_mol.<name>": "amount of the declared species (small species, monomer product, H, H2) (mol)",
    "mass_outside_species_g": "total mass minus the mass of every declared species (g)",
    "ladder_mol.<p>": "polymer-derived molecules of length p primitive periods, p up to the cutoff "
                      "(cutoff x periods per declared unit); a molecule of B backbone atoms spans "
                      "[floor(B / nbb), ceil(B / nbb)] (mol)",
    "mwd_mol.<p>": "the same above the cutoff: the length distribution (mol)",
    "mwd_overflow_mol": "polymer-derived molecules longer than a reference's length grid (mol)",
    "whole_<key>": "melt_observables.WHOLE_ENSEMBLE_DEFINITION",
}
ROW_FIELDS = ("time_s", "row", "producer", "kind", "value", "stderr", "lower", "upper", "seeds", "basis",
              "reason")
JOIN_COLUMNS = ("time_s", "row", "status", "a", "stderr_a", "lower_a", "upper_a", "b", "stderr_b", "lower_b",
                "upper_b", "difference", "difference_stderr", "z", "inside", "reason", "source_a", "source_b")
_WHOLE = ("whole_mu0_mol", "whole_mu1_mol", "whole_mu2_mol", "whole_Mn_g_mol", "whole_Mw_g_mol")
_CENSORED = "censored: {} of {} seeds reached this time"


class ComparisonRowsError(ValueError):
    """A producer output cannot be mapped, or two row documents cannot be joined."""


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def length_interval(backbone: int, nbb: int) -> Tuple[int, int]:
    """[floor, ceil] of a molecule's length in primitive periods."""
    return backbone // nbb, -(-backbone // nbb)


def geometry(spec) -> Tuple[int, int]:
    """(periods per declared unit, backbone atoms per period) of a declaration
    (``declared_chain.ChainGeometry``, the kernel tables' moment unit)."""
    from rmgpy.kmc.declared_chain import ChainGeometry
    g = ChainGeometry.of(spec)
    ppu = g.d_rep // g.period
    if ppu < 1 or g.d_rep % g.period or len(g.backbone_keys) % ppu:
        raise ComparisonRowsError("the declared unit is not a whole number of primitive periods")
    return ppu, len(g.backbone_keys) // ppu


def _length_bins(molecules: Counter, cut: int, unit: float, full: bool) -> Dict[str, Tuple[float, float]]:
    """{row: (lower, upper)} of molecules by period interval {(lo, hi): count}:
    the ladder bins 1..cut always, the bins above when ``full``."""
    lower, upper = Counter(), Counter()
    for (lo, hi), k in molecules.items():
        for p in range(lo, hi + 1):
            upper[p] += k
        if lo == hi:
            lower[lo] += k
    out = {f"ladder_mol.{p}": (unit * lower[p], unit * upper[p]) for p in range(1, cut + 1)}
    if full:
        for p in sorted(upper):
            if p > cut:
                out[f"mwd_mol.{p}"] = (unit * lower[p], unit * upper[p])
    return out


def _row(time_s, name, producer, kind, value=None, stderr=None, lower=None, upper=None, seeds=None,
         basis=None, reason=None) -> dict:
    return {"time_s": time_s, "row": name, "producer": producer, "kind": kind, "value": value, "stderr": stderr,
            "lower": lower, "upper": upper, "seeds": seeds, "basis": basis, "reason": reason}


def _aggregate(time_s, producer, per_seed: List[dict], basis: Dict[str, str]) -> List[dict]:
    """Rows of one time from per-seed values {row: (lower, upper)} (or
    ("diagnostic", value)): a seed mean (with its standard error) where
    every seed is exact, the mean bounds where one is not."""
    names = []
    for values in per_seed:
        names += [n for n in values if n not in names]
    rows = []
    n = len(per_seed)
    for name in names:
        zero = name.startswith(("mwd_mol.", "ladder_mol."))
        pairs = [v.get(name, (0.0, 0.0) if zero else None) for v in per_seed]
        if any(p is None for p in pairs):
            raise ComparisonRowsError(f"row {name} at {time_s} s is missing from a seed")
        family = name.split(".", 1)[0]
        if pairs[0][0] == "diagnostic":
            values = [p[1] for p in pairs]
            rows.append(_row(time_s, name, producer, "diagnostic", value=math.fsum(values) / n, seeds=n,
                             basis=basis.get(family)))
            continue
        se = (lambda xs: statistics.stdev(xs) / math.sqrt(n) if n >= 2 else None)
        if all(lo == hi for lo, hi in pairs):
            values = [lo for lo, _hi in pairs]
            rows.append(_row(time_s, name, producer, "exact", value=math.fsum(values) / n, stderr=se(values),
                             seeds=n, basis=basis.get(family)))
        else:
            rows.append(_row(time_s, name, producer, "interval", lower=math.fsum(p[0] for p in pairs) / n,
                             upper=math.fsum(p[1] for p in pairs) / n, seeds=n, basis=basis.get(family)))
    return rows


def _withheld_time(time_s, producer, names, reason) -> List[dict]:
    return [_row(time_s, name, producer, "withheld", reason=reason) for name in names]


def _document(kind, source, binding, geometry_, withheld_families, rows, **extra) -> dict:
    doc = {"schema": SCHEMA, "producer": dict({"kind": kind}, **source), "binding": binding,
           "geometry": geometry_, "definitions": DEFINITIONS, "withheld_families": withheld_families,
           "rows": rows}
    doc.update(extra)
    return doc


# --------------------------------------------------------------------------
# (a) a kMC run
# --------------------------------------------------------------------------
def _observation_lines(path: Path) -> List[dict]:
    return [o for o in (json.loads(x) for x in path.read_text().splitlines() if x.strip()) if "record" not in o]


def _kmc_values(obs: dict, first: dict, ppu: int, nbb: int, cut: int) -> Tuple[dict, dict, dict]:
    """({row: (lower, upper)}, {family: basis}, {family: withheld reason}) of
    one observation; ``first`` is the same seed's t = 0 observation."""
    from rmgpy.kmc.melt_observables import mom_reported_quantities
    u = obs["scale_factor"] / obs["avogadro"]
    t1 = obs["tier1"]

    def atoms(o):
        blocks = o["tier1"]
        if "composition" in blocks["polymer_derived"]:
            parts = [blocks["polymer_derived"]["composition"], blocks["undeclared_light"]["composition"]]
            parts += [s["composition"] for _name, s in sorted(blocks["species"].items())]
            total = Counter()
            for part in parts:
                total.update(part)
            return total, "sum over the classes (polymer-derived, undeclared light, each declared species)"
        return Counter(o["invariants"]["composition"]), \
            "invariants.composition (the observation has no per-class compositions)"

    now, basis_atoms = atoms(obs)
    then, _b = atoms(first)
    values, basis, withheld = {}, {"element_mol": basis_atoms, "element_residual_mol": basis_atoms}, {}
    for el in sorted(set(now) | set(then)):
        values[f"element_mol.{el}"] = (u * now[el],) * 2
        values[f"element_residual_mol.{el}"] = (u * (now[el] - then[el]),) * 2
    values["total_mass_g"] = (t1["total_mass_g"],) * 2
    species_mass = []
    for name in obs["declared_species"]:
        values[f"species_mol.{name}"] = (t1["species"][name]["amount_mol"],) * 2
        species_mass.append(t1["species"][name]["mass_g"])
    values["mass_outside_species_g"] = (t1["total_mass_g"] - math.fsum(species_mass),) * 2
    molecules = Counter()
    distribution = obs["tier2"].get("backbone_distribution")
    if distribution is not None:
        for b, k in distribution.items():
            molecules[length_interval(int(b), nbb)] += k
        basis["ladder_mol"] = basis["mwd_mol"] = "tier2.backbone_distribution (B per molecule)"
    else:
        for d, k in obs["tier3"]["explicit"]["dp_molecules"].items():
            molecules[(int(d) * ppu, (int(d) + 1) * ppu)] += k
        basis["ladder_mol"] = ("tier3 explicit ladder in declared units: n_c = floor(B / b) = d spans "
                               "[d ppu, (d + 1) ppu] periods")
        withheld["mwd_mol"] = "the observation carries no length distribution (observables schema before 7)"
    values.update(_length_bins(molecules, cut, u, distribution is not None))
    quantities = mom_reported_quantities(obs)
    for key in _WHOLE:
        values[key] = (quantities[key],) * 2
    return values, basis, withheld


def kmc_rows(run_dir, *, pass_name: str = PASS_DECLARED, deck=None) -> dict:
    """The rows of one pass of a kMC run directory (``run.json`` and the
    seeds' observables): ``declared``, or ``well_mixed_bound`` when the
    manifest carries a ``well_mixed`` pass.  The deck (``run.json``'s path
    unless given) must be the run's (sha256)."""
    from rmgpy.kmc import run as kr
    run_dir = Path(run_dir)
    manifest = json.loads((run_dir / "run.json").read_text())
    if pass_name == PASS_DECLARED:
        section = manifest
    elif pass_name == PASS_WELL_MIXED:
        section = manifest.get("well_mixed")
        if not isinstance(section, dict):
            raise ComparisonRowsError(f"{run_dir}: run.json has no well-mixed pass")
    else:
        raise ComparisonRowsError(f"unknown pass {pass_name!r} ({PASS_DECLARED}, {PASS_WELL_MIXED})")
    deck_path = Path(deck) if deck is not None else Path(manifest["deck"]["path"])
    try:
        declared = kr.read_deck(deck_path)
    except (OSError, kr.RunError) as error:
        raise ComparisonRowsError(f"deck {deck_path}: {error}") from error
    if declared.sha256 != manifest["deck"]["sha256"]:
        raise ComparisonRowsError(f"deck {deck_path} is not the run's deck (sha256 {declared.sha256}, "
                                  f"run {manifest['deck']['sha256']})")
    ppu, nbb = geometry(declared.spec)
    times = [float(t) for t in manifest["times_s"]]
    requested = len(manifest["seeds"])
    series = {r["seed"]: _observation_lines(run_dir / r["path"]) if (run_dir / r["path"]).exists() else []
              for r in section["runs"]}
    cut = None
    rows, withheld_families, names = [], {}, []
    for k, t in enumerate(times):
        reached = [s for s, obs in sorted(series.items()) if len(obs) > k]
        if cut is None and reached:
            cut = series[reached[0]][0]["xs"] * ppu
        if len(reached) < requested:
            rows += _withheld_time(t, "kmc", names, _CENSORED.format(len(reached), requested))
            continue
        per_seed, basis = [], {}
        for s in reached:
            values, basis, withheld = _kmc_values(series[s][k], series[s][0], ppu, nbb, cut)
            withheld_families.update(withheld)
            per_seed.append(values)
        block = _aggregate(t, "kmc", per_seed, basis)
        names = names or [r["row"] for r in block if not r["row"].startswith("mwd_mol.")]
        rows += block
    binding = {"deck_sha256": manifest["deck"]["sha256"],
               "chemistry_sha256": manifest["deck"].get("chemistry_sha256"),
               "temperature_K": float(manifest["temperature_K"]), "times_s": times, "pass": pass_name}
    source = {"source": str(run_dir.resolve()),
              "observables_sha256": {str(r["seed"]): r.get("observables_sha256") for r in section["runs"]},
              "seeds_requested": requested}
    return _document("kmc", source, binding, {"periods_per_declared_unit": ppu, "backbone_atoms_per_period": nbb,
                                              "cutoff_periods": cut},
                     withheld_families, rows, mwd_support_periods=[1, None])


# --------------------------------------------------------------------------
# (b) a kernel reference
# --------------------------------------------------------------------------
def _reference_values(obs: dict, first: dict, doc: dict, producer: str, cut: int) -> Tuple[dict, dict]:
    u = doc["unit_mol"]
    unbalanced = doc.get("unbalanced_terms") or []
    values, basis = {}, {}
    elements = sorted(set(obs["elements"]) | set(first["elements"]))
    if producer == "pbe" and not unbalanced:
        basis["element_mol"] = basis["element_residual_mol"] = basis["total_mass_g"] = (
            "model invariant: every term balances each element exactly (checked at model build); the "
            "integrated value's round-off is element_drift_mol")
        for el in elements:
            values[f"element_mol.{el}"] = (u * first["elements"].get(el, 0),) * 2
            values[f"element_residual_mol.{el}"] = (0.0, 0.0)
            values[f"element_drift_mol.{el}"] = ("diagnostic", u * (obs["elements"].get(el, 0)
                                                                    - first["elements"].get(el, 0)))
        values["total_mass_g"] = (u * first["mass_g"],) * 2
    else:
        note = ("integer state counts" if producer == "ssa" else "integrated state counts")
        if unbalanced:
            note += "; terms that do not balance: " + ", ".join(x["template"] for x in unbalanced)
        basis["element_mol"] = basis["element_residual_mol"] = basis["total_mass_g"] = note
        for el in elements:
            now, then = obs["elements"].get(el, 0), first["elements"].get(el, 0)
            values[f"element_mol.{el}"] = (u * now,) * 2
            values[f"element_residual_mol.{el}"] = (u * (now - then),) * 2
        values["total_mass_g"] = (u * obs["mass_g"],) * 2
    for name in doc["declared_species"]:
        values[f"species_mol.{name}"] = (u * obs["species"].get(name, 0.0),) * 2
    values["mass_outside_species_g"] = (u * (obs["mass_g"] - obs["species_mass_g"]),) * 2
    molecules = Counter({(int(p), int(p)): k for p, k in obs["lengths"].items()})
    values.update(_length_bins(molecules, cut, u, True))
    values["mwd_overflow_mol"] = (u * obs["overflow_molecules"],) * 2
    basis["ladder_mol"] = basis["mwd_mol"] = "states on the period grid (n periods)"
    for key in _WHOLE:
        values[key] = (obs["whole"][key],) * 2
    return values, basis


def reference_rows(path, producer: str) -> dict:
    """The rows of a ``kernel_reference`` output: its population balance
    (``pbe``) or its kernel SSA (``ssa``, seed means)."""
    from rmgpy.kmc import kernel_reference as R
    path = Path(path)
    raw = path.read_bytes()
    doc = json.loads(raw)
    if doc.get("schema") != R.SCHEMA:
        raise ComparisonRowsError(f"{path}: not a {R.SCHEMA} output (got {doc.get('schema')!r})")
    if producer not in ("pbe", "ssa"):
        raise ComparisonRowsError(f"producer must be pbe or ssa, got {producer!r}")
    geometry_ = doc["geometry"]
    cut = doc["deck"]["cutoff"] * geometry_["periods_per_declared_unit"]
    times = [float(t) for t in doc["times_s"]]
    observed = doc["row_observations"]
    runs = [observed["pbe"]] if producer == "pbe" else observed["ssa"]["runs"]
    rows = []
    for k, t in enumerate(times):
        per_seed, basis = [], {}
        for series in runs:
            values, basis = _reference_values(series[k], series[0], doc, producer, cut)
            per_seed.append(values)
        rows += _aggregate(t, f"kernel_{producer}", per_seed, basis)
    binding = {"deck_sha256": doc["deck"]["sha256"], "chemistry_sha256": doc["deck"]["chemistry_sha256"],
               "temperature_K": float(doc["temperature_K"]), "times_s": times, "pass": doc["pass"]}
    source = {"source": str(path.resolve()), "sha256": _sha256(raw), "table_sha256": doc["table_sha256"],
              "phi": doc["phi"], "unbalanced_terms": doc.get("unbalanced_terms") or [],
              "not_modelled_terms": doc.get("not_modelled_terms") or [],
              "seeds": None if producer == "pbe" else observed["ssa"]["seeds"]}
    return _document(f"kernel_{producer}", source, binding, dict(geometry_, cutoff_periods=cut), {}, rows,
                     mwd_support_periods=[1, doc["n_max"]])


# --------------------------------------------------------------------------
# (c) a MoM --mom file
# --------------------------------------------------------------------------
MOM_NO_INVENTORY = ("the --mom file carries no element or total inventory (only the quantities of "
                    "melt_observables.mom_reported_quantities)")
MOM_MOMENTS_ONLY = "MoM reports moments: its pools hold no length distribution"


def mom_rows(path, deck) -> dict:
    """The rows of a header-form ``--mom`` file of ``deck`` (its sha256 must
    be the header's ``producer.deck_sha256``).  MoM's rate coefficients are
    intrinsic (no cage): its pass is the well-mixed bound."""
    from rmgpy.kmc import run as kr
    from rmgpy.kmc.melt_engine import SpeciesIdentity
    from rmgpy.kmc.melt_observables import declared_species_names
    path = Path(path)
    raw = path.read_bytes()
    doc = json.loads(raw)
    if not (isinstance(doc, dict) and set(doc) == {"header", "entries"}):
        raise ComparisonRowsError(f"{path}: not a header-form --mom file")
    header, entries = doc["header"], doc["entries"]
    declared = kr.read_deck(deck)
    bound = (header.get("producer") or {}).get("deck_sha256")
    if bound != declared.sha256:
        raise ComparisonRowsError(f"{path}: the header's deck digest {bound!r} is not the deck {deck} "
                                  f"({declared.sha256})")
    conventions_path = path.parent / "conventions.json"
    not_produced = {}
    if conventions_path.exists():
        conventions = json.loads(conventions_path.read_text())
        if header.get("conventions_sha256") == _sha256(conventions_path.read_bytes()):
            not_produced = conventions.get("not_produced") or {}
    ppu, nbb = geometry(declared.spec)
    cut = declared.xs * ppu
    names = declared_species_names(declared.spec, SpeciesIdentity(declared.spec))
    ladder = not_produced.get("explicit_dp_mol") or MOM_MOMENTS_ONLY
    rows = []
    times = sorted(float(e["time_s"]) for e in entries)
    for entry in sorted(entries, key=lambda e: float(e["time_s"])):
        t, q = float(entry["time_s"]), entry["quantities"]
        species = q.get("species_mol") or {}
        for name in names:
            if name in species:
                rows.append(_row(t, f"species_mol.{name}", "mom", "exact", value=float(species[name])))
            else:
                rows.append(_row(t, f"species_mol.{name}", "mom", "withheld",
                                 reason=not_produced.get("species_mol") or "not exported by the MoM producer"))
        for key in _WHOLE:
            if key in q:
                rows.append(_row(t, key, "mom", "exact", value=float(q[key])))
            else:
                rows.append(_row(t, key, "mom", "withheld",
                                 reason=not_produced.get(key) or "not exported by the MoM producer"))
        for name in ("total_mass_g", "mass_outside_species_g"):
            rows.append(_row(t, name, "mom", "withheld", reason=MOM_NO_INVENTORY))
        rows += [_row(t, f"ladder_mol.{p}", "mom", "withheld", reason=ladder) for p in range(1, cut + 1)]
    binding = {"deck_sha256": declared.sha256, "chemistry_sha256": kr.chemistry_sha256(declared.spec),
               "temperature_K": float((header.get("run") or {}).get("temperature_K")), "times_s": times,
               "pass": PASS_WELL_MIXED}
    source = {"source": str(path.resolve()), "sha256": _sha256(raw), "flags": sorted(header.get("flags", {}))}
    return _document("mom", source, binding, {"periods_per_declared_unit": ppu, "backbone_atoms_per_period": nbb,
                                              "cutoff_periods": cut},
                     {"element_mol": MOM_NO_INVENTORY, "element_residual_mol": MOM_NO_INVENTORY,
                      "mwd_mol": MOM_MOMENTS_ONLY}, rows, mwd_support_periods=None)


# --------------------------------------------------------------------------
# the join
# --------------------------------------------------------------------------
def _check_binding(a: dict, b: dict):
    from rmgpy.kmc.run import MOM_TIME_RELATIVE
    x, y = a["binding"], b["binding"]
    if x["deck_sha256"] != y["deck_sha256"]:
        raise ComparisonRowsError(f"the deck differs: {x['deck_sha256']} against {y['deck_sha256']}")
    if float(x["temperature_K"]) != float(y["temperature_K"]):
        raise ComparisonRowsError(f"the temperature differs: {x['temperature_K']} K against {y['temperature_K']} K")
    tx, ty = [float(t) for t in x["times_s"]], [float(t) for t in y["times_s"]]
    if len(tx) != len(ty) or any(abs(p - q) > MOM_TIME_RELATIVE * max(abs(p), abs(q)) for p, q in zip(tx, ty)):
        raise ComparisonRowsError(f"the time grid differs: {tx} against {ty}")
    if x["pass"] != y["pass"]:
        raise ComparisonRowsError(f"the pass differs: {x['pass']} against {y['pass']}")


def _side(doc: dict, table: dict, k: int, name: str) -> Optional[dict]:
    row = table.get((k, name))
    if row is not None:
        return row
    family = name.split(".", 1)[0]
    reason = doc.get("withheld_families", {}).get(family)
    if reason is not None:
        return {"kind": "withheld", "reason": reason}
    support = doc.get("mwd_support_periods")
    if family == "mwd_mol" and support is not None:
        p = int(name.split(".", 1)[1])
        if p >= support[0] and (support[1] is None or p <= support[1]):
            return {"kind": "exact", "value": 0.0, "stderr": None}
        return {"kind": "withheld", "reason": f"beyond the producer's length grid ({support[1]} periods)"}
    return None


def join(a: dict, b: dict) -> List[dict]:
    """Per time and row of two row documents of one deck, temperature, time
    grid and pass (``ComparisonRowsError`` naming the field otherwise)."""
    for doc in (a, b):
        if doc.get("schema") != SCHEMA:
            raise ComparisonRowsError(f"not a {SCHEMA} document (got {doc.get('schema')!r})")
    _check_binding(a, b)
    times = [float(t) for t in a["binding"]["times_s"]]
    index = []
    for doc in (a, b):
        grid = [float(t) for t in doc["binding"]["times_s"]]
        index.append({(grid.index(float(r["time_s"])), r["row"]): r for r in doc["rows"]})
    out = []
    for k, t in enumerate(times):
        names = []
        for table in index:
            names += [name for (kk, name) in table if kk == k and name not in names]
        for name in names:
            ra, rb = _side(a, index[0], k, name), _side(b, index[1], k, name)
            row = {c: None for c in JOIN_COLUMNS}
            row.update(time_s=t, row=name, source_a=a["producer"]["kind"], source_b=b["producer"]["kind"])
            for tag, r in (("a", ra), ("b", rb)):
                if r is not None:
                    row[tag], row[f"stderr_{tag}"] = r.get("value"), r.get("stderr")
                    row[f"lower_{tag}"], row[f"upper_{tag}"] = r.get("lower"), r.get("upper")
            kinds = (ra and ra["kind"], rb and rb["kind"])
            if "diagnostic" in kinds:
                row["status"] = "not_a_comparison"
            elif ra is None or rb is None:
                row["status"] = "missing_a" if ra is None else "missing_b"
            elif "withheld" in kinds:
                row["status"] = {("withheld", "withheld"): "withheld_both"}.get(
                    kinds, "withheld_a" if kinds[0] == "withheld" else "withheld_b")
                row["reason"] = "; ".join(r["reason"] for r in (ra, rb) if r["kind"] == "withheld")
            elif "interval" in kinds:
                lo = [r["value"] if r["kind"] == "exact" else r["lower"] for r in (ra, rb)]
                hi = [r["value"] if r["kind"] == "exact" else r["upper"] for r in (ra, rb)]
                row["status"] = "interval"
                row["inside"] = max(lo) <= min(hi)
            else:
                row["status"] = "compared"
                row["difference"] = ra["value"] - rb["value"]
                errors = [r["stderr"] for r in (ra, rb) if r.get("stderr") is not None]
                if errors:
                    se = math.sqrt(math.fsum(e * e for e in errors))
                    row["difference_stderr"] = se
                    row["z"] = row["difference"] / se if se > 0.0 else None
            out.append(row)
    return out


def _cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return repr(value)
    return str(value)


def join_text(a: dict, b: dict, rows: List[dict]) -> str:
    buffer = io.StringIO()
    head = {"a": a["producer"].get("source"), "b": b["producer"].get("source"),
            "deck_sha256": a["binding"]["deck_sha256"], "temperature_K": a["binding"]["temperature_K"],
            "pass": a["binding"]["pass"]}
    buffer.write(f"# {JOIN_SCHEMA} " + " ".join(f"{k}={json.dumps(v)}" for k, v in head.items()) + "\n")
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(JOIN_COLUMNS)
    for row in rows:
        writer.writerow([_cell(row[c]) for c in JOIN_COLUMNS])
    return buffer.getvalue()


# --------------------------------------------------------------------------
# the command line
# --------------------------------------------------------------------------
def _write(doc: dict, path: Path):
    path.write_text(json.dumps(doc, sort_keys=True, indent=1, allow_nan=False) + "\n")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m rmgpy.kmc.comparison_rows",
                                     description="comparison rows of kMC, kernel-reference and MoM outputs, and "
                                                 "their join")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("kmc", help="rows of a kMC run directory")
    p.add_argument("run_dir")
    p.add_argument("--pass", dest="pass_name", default=PASS_DECLARED, choices=(PASS_DECLARED, PASS_WELL_MIXED))
    p.add_argument("--deck", default=None, help="the run's deck, when run.json's path moved")
    p.add_argument("--out", default=None, help="default: where the runner writes them, <run_dir>/rows.json or "
                   "<run_dir>/well_mixed/rows.json")
    p = sub.add_parser("reference", help="rows of a kernel_reference output")
    p.add_argument("reference")
    p.add_argument("--producer", required=True, choices=("pbe", "ssa"))
    p.add_argument("--out", default=None, help="default: <reference stem>.<producer>.rows.json")
    p = sub.add_parser("mom", help="rows of a header-form --mom file")
    p.add_argument("mom")
    p.add_argument("--deck", required=True)
    p.add_argument("--out", default=None, help="default: <mom stem>.rows.json")
    p = sub.add_parser("join", help="the per-row, per-time comparison of two row files")
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("--out", default=None, help="CSV (default: standard output)")
    args = parser.parse_args(argv)
    try:
        if args.command == "kmc":
            doc = kmc_rows(args.run_dir, pass_name=args.pass_name, deck=args.deck)
            where = "rows.json" if args.pass_name == PASS_DECLARED else "well_mixed/rows.json"
            out = Path(args.out) if args.out else Path(args.run_dir) / where
        elif args.command == "reference":
            doc = reference_rows(args.reference, args.producer)
            out = Path(args.out) if args.out else Path(args.reference).with_suffix(f".{args.producer}.rows.json")
        elif args.command == "mom":
            doc = mom_rows(args.mom, args.deck)
            out = Path(args.out) if args.out else Path(args.mom).with_suffix(".rows.json")
        else:
            a, b = (json.loads(Path(x).read_text()) for x in (args.a, args.b))
            text = join_text(a, b, join(a, b))
            if args.out:
                Path(args.out).write_text(text)
            else:
                sys.stdout.write(text)
            return 0
    except (ComparisonRowsError, OSError, ValueError, KeyError) as error:
        print(f"refused: {error}", file=sys.stderr, flush=True)
        return 2
    _write(doc, out)
    print(f"written: {out} ({len(doc['rows'])} rows)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
