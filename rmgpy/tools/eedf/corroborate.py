#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""Run LoKI-B and bolos live from one canonical physical setup."""

from __future__ import annotations

import argparse
import copy
from dataclasses import asdict, dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import tempfile
import time
import uuid

import numpy as np
from scipy import integrate
from scipy.constants import Boltzmann, elementary_charge, electron_mass

from rmgpy.tools.eedf.channels import _classic_sections, validate_physical_map
from rmgpy.tools.eedf.integrity import channel_fractions, resolved_properties, solver_environment
from rmgpy.tools.eedf.loki import LoKIError, parse_output, setup_text
from rmgpy.tools.eedf.schema import SpecError, canonical_json


class CorroborationError(ValueError):
    """The requested comparison cannot emit an accepted artifact."""


class ExecutionBackendUnavailable(CorroborationError):
    """The qualified descendant-owning execution backend is unavailable."""


GAMMA = math.sqrt(2.0 * elementary_charge / electron_mass)
SETUP_FINGERPRINT_VERSION = "eedf-canonical-setup-v1"
COMPARISON_RTOL = 0.01
BOLOS_MAX_ITERATIONS = 500
BOLOS_CONVERGENCE_RTOL = 1.0e-7
TRANSPORT_CONSISTENCY_RTOL = 5.0e-3
REFERENCE_EEDF_NEGATIVE_MASS_FLOOR = 1.0e-12
_CLEANUP_GRACE_S = 2.0
_CLEANUP_ATTEMPTS = 3


@dataclass(frozen=True)
class SourceSnapshot:
    """One provenance-pinned file, held as immutable bytes."""

    name: str
    kind: str
    sha256: str
    content: bytes


@dataclass(frozen=True)
class Collision:
    """One collision as interpreted from the pinned LoKI and BOLSIG bytes."""

    description: str
    kind: str
    bolos_target: str
    bolos_product: str | None
    threshold_eV: float
    mass_ratio: float | None
    energy_eV: tuple[float, ...]
    sigma_m2: tuple[float, ...]


@dataclass(frozen=True)
class NumericalSettings:
    """Frozen solve and comparison settings with declared provenance."""

    max_energy_eV: float
    cell_count: int
    loki_residual_tolerance: float
    loki_power_balance_tolerance: float
    bolos_residual_tolerance: float
    bolos_max_iterations: int
    comparison_rtol: float
    transport_consistency_rtol: float
    timeout_s: float
    normalization_atol: float
    fingerprint_rtol: float


@dataclass(frozen=True)
class CanonicalSetup:
    """The sole physical-input authority for both solver adapters."""

    fingerprint_version: str
    fingerprint: str
    sources: tuple[SourceSnapshot, ...]
    collisions: tuple[Collision, ...]
    gas_temperature_K: float
    pressure_Pa: float
    reduced_field_Td: float
    gas_fractions: tuple[tuple[str, float], ...]
    state_populations: tuple[tuple[str, float], ...]
    target_fractions: tuple[float, ...]
    product_fractions: tuple[float, ...]
    ionization_energy_sharing: str
    electron_growth: str
    excitation_frequency_Hz: float
    eedf_mode: str
    electron_electron_collisions: bool
    mobility_definition: str
    compared_quantities: tuple[str, ...]
    numerical: NumericalSettings
    loki_setup: bytes
    loki_binary_path: str
    loki_binary_sha256: str
    loki_environment: tuple[tuple[str, str], ...]


def _sha256_bytes(content):
    return hashlib.sha256(content).hexdigest()


def _finite_float(value, name, *, positive=False):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise CorroborationError(name + " must be a finite number")
    value = float(value)
    if positive and value <= 0.0:
        raise CorroborationError(name + " must be positive")
    return value


def _strict_int(value, name, *, positive=False):
    if type(value) is not int:
        raise CorroborationError(name + " must be an integer")
    if positive and value <= 0:
        raise CorroborationError(name + " must be positive")
    return value


def _numeric_array(values, name):
    if not isinstance(values, (list, tuple, np.ndarray)):
        raise CorroborationError(name + " must be a numeric array")
    flattened = np.asarray(values).reshape(-1)
    if any(type(value.item() if isinstance(value, np.generic) else value) not in (int, float) for value in flattened):
        raise CorroborationError(name + " must be a numeric array")
    result = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(result)):
        raise CorroborationError(name + " contains a non-finite value")
    return result


def _safe_source_name(name):
    path = Path(name)
    if path.is_absolute() or ".." in path.parts or not path.parts:
        raise CorroborationError("unsafe pinned source name: " + str(name))
    return path


def _read_artifact_bytes(artifact):
    root = Path(artifact)
    if root.is_file():
        root = root.parent
    try:
        manifest_bytes = (root / "manifest.json").read_bytes()
        spec_bytes = (root / "generation_spec.json").read_bytes()
        manifest = json.loads(manifest_bytes)
        spec = json.loads(spec_bytes)
    except (OSError, ValueError) as exc:
        raise CorroborationError("incomplete artifact: " + str(exc)) from exc
    if not isinstance(manifest, dict) or not isinstance(spec, dict):
        raise CorroborationError("artifact metadata must be objects")
    return root, manifest, spec, manifest_bytes, spec_bytes


def _snapshot_sources(spec):
    try:
        registry = spec["input_files"]
    except KeyError as exc:
        raise CorroborationError("generation spec is missing input_files") from exc
    if not isinstance(registry, dict) or not registry:
        raise CorroborationError("generation spec input_files is empty or invalid")
    snapshots = []
    for name, entry in registry.items():
        if not isinstance(entry, dict) or set(entry) != {"path", "sha256", "kind"}:
            raise CorroborationError("invalid pinned source registry entry: " + str(name))
        _safe_source_name(name)
        try:
            content = Path(entry["path"]).read_bytes()
        except OSError as exc:
            raise CorroborationError("pinned source is unavailable: " + str(name)) from exc
        digest = _sha256_bytes(content)
        if digest != entry["sha256"]:
            raise CorroborationError("pinned source SHA-256 mismatch: " + str(name))
        snapshots.append(SourceSnapshot(name, entry["kind"], digest, content))
    return tuple(snapshots)


def _materialized_spec(spec, snapshots, root):
    materialized = copy.deepcopy(spec)
    for snapshot in snapshots:
        path = root / _safe_source_name(snapshot.name)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(snapshot.content)
        materialized["input_files"][snapshot.name]["path"] = str(path)
    return materialized


_BOLOS_KINDS = {"ELASTIC", "MOMENTUM", "EFFECTIVE", "EXCITATION", "IONIZATION", "ATTACHMENT"}


def _parse_bolos_blocks(content, source_name):
    """Parse the physical fields that bolos reads from BOLSIG-compatible bytes."""
    try:
        lines = content.decode("latin-1").splitlines()
    except UnicodeDecodeError as exc:
        raise CorroborationError("collision source is not latin-1: " + source_name) from exc
    blocks = []
    index = 0
    while index < len(lines):
        kind = lines[index].strip()
        index += 1
        if kind not in _BOLOS_KINDS:
            continue
        if index >= len(lines):
            raise CorroborationError("truncated BOLSIG collision header: " + source_name)
        target_line = lines[index].strip()
        index += 1
        argument = None
        if kind != "ATTACHMENT":
            if index >= len(lines):
                raise CorroborationError("missing BOLSIG collision argument: " + source_name)
            argument = lines[index].strip()
            index += 1
        comments = []
        while index < len(lines) and not re.fullmatch(r"-{5,}", lines[index].strip()):
            comments.append(lines[index])
            index += 1
        if index >= len(lines):
            raise CorroborationError("missing BOLSIG data delimiter: " + source_name)
        index += 1
        data = []
        while index < len(lines) and not re.fullmatch(r"-{5,}", lines[index].strip()):
            fields = lines[index].split()
            if len(fields) != 2:
                raise CorroborationError("ambiguous BOLSIG data row: " + source_name)
            try:
                data.append((float(fields[0]), float(fields[1])))
            except ValueError as exc:
                raise CorroborationError("nonnumeric BOLSIG data row: " + source_name) from exc
            index += 1
        if index >= len(lines) or len(data) < 2:
            raise CorroborationError("incomplete BOLSIG collision data: " + source_name)
        index += 1
        if kind in {"ELASTIC", "MOMENTUM", "EFFECTIVE"}:
            mass_ratio = _finite_float(float(argument.split()[0]), "bolos header mass ratio", positive=True)
            threshold = 0.0
            target = target_line
            product = None
        elif kind == "ATTACHMENT":
            mass_ratio = None
            threshold = 0.0
            sides = re.split(r"<?->", target_line)
            target = sides[0].strip()
            product = sides[1].strip() if len(sides) == 2 else None
        else:
            mass_ratio = None
            threshold = _finite_float(float(argument.split()[0]), "bolos collision threshold")
            sides = re.split(r"<?->", target_line)
            if len(sides) != 2:
                raise CorroborationError("ambiguous BOLSIG target/product: " + source_name)
            target, product = (side.strip() for side in sides)
        blocks.append(
            {
                "kind": kind.lower(),
                "target": target,
                "product": product,
                "threshold": threshold,
                "mass_ratio": mass_ratio,
                "energy": tuple(item[0] for item in data),
                "sigma": tuple(item[1] for item in data),
                "comments": tuple(comments),
            }
        )
    if not blocks:
        raise CorroborationError("pinned collision source has no BOLSIG blocks: " + source_name)
    return blocks


def _channel_participants(description, kind):
    """Return the BOLSIG target/product identities encoded by a LoKI channel."""
    match = re.fullmatch(r"e \+ (.+?) -> (.+), [^,]+", description)
    if match is None:
        raise CorroborationError("ambiguous collision identity: " + description)
    target_state, products = match.groups()
    target = target_state.split("(", 1)[0]
    while products.startswith("e + "):
        products = products[4:]
    if kind == "elastic":
        return target, None
    ion = re.fullmatch(r"([^()]+)\(\+,gnd\)", products)
    product = ion.group(1) + "+" if ion else products
    return target, product


def _canonical_collisions(spec, manifest, coordinates, snapshots, materialized_root):
    materialized = _materialized_spec(spec, snapshots, materialized_root)
    channels = manifest["channel_map"] if "channel_map" in manifest else None
    if not isinstance(channels, list) or not channels:
        raise CorroborationError("manifest channel_map is missing or empty")
    try:
        validate_physical_map(channels, materialized, coordinates)
    except SpecError as exc:
        raise CorroborationError(str(exc)) from exc
    raw_blocks = []
    classic = []
    for snapshot in snapshots:
        if snapshot.kind != "cross_section":
            continue
        blocks = _parse_bolos_blocks(snapshot.content, snapshot.name)
        raw_blocks.extend(blocks)
        entry = materialized["input_files"][snapshot.name]
        try:
            classic.extend(_classic_sections(entry["path"], entry["sha256"]))
        except SpecError as exc:
            raise CorroborationError(str(exc)) from exc
    if len(raw_blocks) != len(classic) or len(classic) != len(channels):
        raise CorroborationError("LoKI and bolos collision counts mismatch")
    classic_by_description = {item["description"]: item for item in classic}
    collisions = []
    for channel, raw, ordered_parsed in zip(channels, raw_blocks, classic):
        description = channel["description"]
        parsed = classic_by_description.get(description)
        if parsed is None:
            raise CorroborationError("collision identity mismatch: " + description)
        if ordered_parsed["description"] != description:
            raise CorroborationError("collision ordering mismatch: " + description)
        expected_target, expected_product = _channel_participants(description, channel["kind"])
        product_matches = (
            raw["product"] == expected_product
            or isinstance(raw["product"], str)
            and isinstance(expected_product, str)
            and raw["product"].replace(" ", "") == expected_product.replace(" ", "")
        )
        if raw["target"] != expected_target or not product_matches:
            raise CorroborationError("BOLSIG collision participant mismatch: " + description)
        expected_kind = "elastic" if raw["kind"] in {"elastic", "momentum", "effective"} else raw["kind"]
        if expected_kind != channel["kind"]:
            raise CorroborationError("collision kind mismatch: " + description)
        if raw["energy"] != tuple(parsed["cross_section"]["energy_eV"]) or raw["sigma"] != tuple(parsed["cross_section"]["sigma_m2"]):
            raise CorroborationError("collision array mismatch: " + description)
        if raw["threshold"] != float(parsed["threshold_eV"]):
            raise CorroborationError("collision threshold mismatch: " + description)
        canonical_mass = channel["mass_ratio"] if "mass_ratio" in channel else None
        if channel["kind"] == "elastic":
            if raw["mass_ratio"] != canonical_mass:
                raise CorroborationError("bolos header mass setup mismatch: " + description)
        elif raw["mass_ratio"] is not None:
            raise CorroborationError("unexpected collision mass convention: " + description)
        collisions.append(
            Collision(
                description=description,
                kind=channel["kind"],
                bolos_target=raw["target"],
                bolos_product=raw["product"],
                threshold_eV=float(parsed["threshold_eV"]),
                mass_ratio=float(canonical_mass) if canonical_mass is not None else None,
                energy_eV=tuple(float(value) for value in raw["energy"]),
                sigma_m2=tuple(float(value) for value in raw["sigma"]),
            )
        )
    return tuple(collisions), materialized


def _same_effective_value(actual, expected):
    """Compare serialized native settings without rejecting float roundoff."""
    if type(expected) is float:
        return type(actual) in (int, float) and math.isfinite(actual) and math.isclose(
            actual, expected, rel_tol=1e-14, abs_tol=0.0
        )
    if isinstance(expected, dict):
        return (
            isinstance(actual, dict)
            and set(actual) == set(expected)
            and all(_same_effective_value(actual[key], value) for key, value in expected.items())
        )
    if isinstance(expected, (list, tuple)):
        return (
            isinstance(actual, (list, tuple))
            and len(actual) == len(expected)
            and all(_same_effective_value(left, right) for left, right in zip(actual, expected))
        )
    return type(actual) is type(expected) and actual == expected


def build_canonical_setup(artifact, branch, node, u):
    """Parse and freeze the sole physical setup accepted by both adapters."""
    if node is None:
        node = {}
    if not isinstance(node, dict):
        raise CorroborationError("composition node must be an object")
    forbidden = sorted(set(node) & {"Tg_K", "P_Pa", "EN_Td", "u"})
    if forbidden:
        raise CorroborationError(forbidden[0] + " caller override is forbidden")
    root, manifest, spec, manifest_bytes, spec_bytes = _read_artifact_bytes(artifact)
    snapshots = _snapshot_sources(spec)
    try:
        branches = manifest["branches"]
        axes = manifest["axes"]
        u_axis = axes["u"]
    except (KeyError, TypeError) as exc:
        raise CorroborationError("artifact point axes are incomplete") from exc
    if not isinstance(branches, list) or branch not in branches:
        raise CorroborationError("requested branch is unavailable")
    if not isinstance(u_axis, list) or not any(_same_effective_value(value, u) for value in u_axis):
        raise CorroborationError("requested reduced-field point is unavailable")
    coordinate_axes = set(axes) - {"u"}
    if set(node) != coordinate_axes:
        raise CorroborationError("composition coordinates do not exactly match the artifact axes")
    for name, value in node.items():
        if not isinstance(axes[name], list) or not any(_same_effective_value(candidate, value) for candidate in axes[name]):
            raise CorroborationError("requested composition point is unavailable: " + name)
    reduced_field = math.exp(_finite_float(u, "u"))
    with tempfile.TemporaryDirectory(prefix="eedf-canonical-parse-") as temporary:
        collisions, materialized = _canonical_collisions(spec, manifest, node, snapshots, Path(temporary))
        try:
            _, _, gas_fractions, populations = resolved_properties(materialized, node)
            expected_targets, expected_products = channel_fractions(
                manifest["channel_map"], gas_fractions, populations, node
            )
            rendered_setup = setup_text(materialized, node, [reduced_field], "solve").encode("utf-8")
        except (SpecError, LoKIError) as exc:
            raise CorroborationError(str(exc)) from exc
    if len(set(collision.description for collision in collisions)) != len(collisions):
        raise CorroborationError("duplicated collision target identity")
    options = spec["solver_options"] if "solver_options" in spec else None
    if not isinstance(options, dict):
        raise CorroborationError("solver_options are missing")
    ionization = options["ionizationOperatorType"] if "ionizationOperatorType" in options else None
    growth = options["growthModelType"] if "growthModelType" in options else None
    if ionization != "equalSharing":
        raise CorroborationError("ionization energy-sharing setup mismatch: equalSharing required")
    if growth != "temporal":
        raise CorroborationError("electron-growth setup mismatch: temporal required")
    if options.get("eedfType") != "boltzmann":
        raise CorroborationError("EEDF mode setup mismatch: boltzmann required")
    if "includeEECollisions" not in options or options["includeEECollisions"] is not False:
        raise CorroborationError("electron-electron collisions are outside common capability")
    if any("<->" in collision.description or collision.kind in {"rotational", "vibrational"} for collision in collisions):
        raise CorroborationError("no second code available for superelastic/rotational populations")
    numerical = options["numerics"] if "numerics" in options else None
    try:
        grid = numerical["energyGrid"]
        loki_tolerance = numerical["nonLinearRoutines"]["maxEedfRelError"]
        power_tolerance = numerical["maxPowerBalanceRelError"]
        tolerances = spec["tolerances"]
        normalization_atol = tolerances["normalization_atol"]
        fingerprint_rtol = tolerances["fingerprint_rtol"]
    except (KeyError, TypeError) as exc:
        raise CorroborationError("frozen numerical settings are incomplete") from exc
    settings = NumericalSettings(
        max_energy_eV=_finite_float(grid["maxEnergy"] if "maxEnergy" in grid else None, "maxEnergy", positive=True),
        cell_count=_strict_int(grid["cellNumber"] if "cellNumber" in grid else None, "cellNumber", positive=True),
        loki_residual_tolerance=_finite_float(loki_tolerance, "LoKI convergence tolerance", positive=True),
        loki_power_balance_tolerance=_finite_float(
            power_tolerance, "LoKI power-balance tolerance", positive=True
        ),
        bolos_residual_tolerance=BOLOS_CONVERGENCE_RTOL,
        bolos_max_iterations=BOLOS_MAX_ITERATIONS,
        comparison_rtol=COMPARISON_RTOL,
        transport_consistency_rtol=TRANSPORT_CONSISTENCY_RTOL,
        timeout_s=_finite_float(spec["timeout_s"] if "timeout_s" in spec else None, "timeout_s", positive=True),
        normalization_atol=_finite_float(normalization_atol, "normalization_atol", positive=True),
        fingerprint_rtol=_finite_float(fingerprint_rtol, "fingerprint_rtol", positive=True),
    )
    if settings.cell_count < 2:
        raise CorroborationError("energy grid setup mismatch")
    binary = spec.get("binary")
    if not isinstance(binary, dict) or set(binary) != {"path", "sha256"}:
        raise CorroborationError("LoKI binary pin is missing or malformed")
    if not isinstance(binary["path"], str) or not isinstance(binary["sha256"], str):
        raise CorroborationError("LoKI binary pin is missing or malformed")
    sources_with_metadata = snapshots + (
        SourceSnapshot("manifest.json", "artifact_metadata", _sha256_bytes(manifest_bytes), manifest_bytes),
        SourceSnapshot("generation_spec.json", "artifact_metadata", _sha256_bytes(spec_bytes), spec_bytes),
    )
    payload = {
        "version": SETUP_FINGERPRINT_VERSION,
        "sources": [(item.name, item.kind, item.sha256) for item in sources_with_metadata],
        "collisions": [asdict(collision) for collision in collisions],
        "Tg_K": spec.get("Tg_K"),
        "P_Pa": spec.get("P_Pa"),
        "EN_Td": reduced_field,
        "gas_fractions": sorted(gas_fractions.items()),
        "state_populations": sorted(populations.items()),
        "targets": expected_targets,
        "products": expected_products,
        "ionization_energy_sharing": ionization,
        "electron_growth": growth,
        "excitation_frequency_Hz": spec.get("working_conditions", {}).get("excitationFrequency"),
        "eedf_mode": options.get("eedfType"),
        "electron_electron_collisions": options["includeEECollisions"],
        "mobility_definition": "temporal-growth-corrected reduced mobility",
        "compared_quantities": ["mean_energy_eV", "mobility_N", "rates", "channel_power", "total_power"],
        "numerical": asdict(settings),
        "loki_setup_sha256": _sha256_bytes(rendered_setup),
        "loki_binary_sha256": binary["sha256"],
    }
    fingerprint = hashlib.sha256(canonical_json(payload).encode()).hexdigest()
    setup = CanonicalSetup(
        fingerprint_version=SETUP_FINGERPRINT_VERSION,
        fingerprint=fingerprint,
        sources=sources_with_metadata,
        collisions=collisions,
        gas_temperature_K=_finite_float(spec.get("Tg_K"), "Tg_K", positive=True),
        pressure_Pa=_finite_float(spec.get("P_Pa"), "P_Pa", positive=True),
        reduced_field_Td=reduced_field,
        gas_fractions=tuple(sorted((name, _finite_float(value, "gas fraction " + name)) for name, value in gas_fractions.items())),
        state_populations=tuple(sorted((name, _finite_float(value, "state population " + name)) for name, value in populations.items())),
        target_fractions=tuple(
            float(value)
            for value in _numeric_array(expected_targets, "canonical target fractions")
        ),
        product_fractions=tuple(
            float(value)
            for value in _numeric_array(expected_products, "canonical product fractions")
        ),
        ionization_energy_sharing=ionization,
        electron_growth=growth,
        excitation_frequency_Hz=_finite_float(
            spec.get("working_conditions", {}).get("excitationFrequency"),
            "excitationFrequency",
        ),
        eedf_mode=options.get("eedfType"),
        electron_electron_collisions=options["includeEECollisions"],
        mobility_definition="temporal-growth-corrected reduced mobility",
        compared_quantities=("mean_energy_eV", "mobility_N", "rates", "channel_power", "total_power"),
        numerical=settings,
        loki_setup=rendered_setup,
        loki_binary_path=binary["path"],
        loki_binary_sha256=binary["sha256"],
        loki_environment=tuple(sorted(solver_environment(spec).items())),
    )
    return setup


def _process_dicts(collisions):
    return [
        {
            "description": collision.description,
            "kind": collision.kind,
            "threshold": collision.threshold_eV,
            "mass_ratio": collision.mass_ratio,
            "energy": collision.energy_eV,
            "sigma": collision.sigma_m2,
        }
        for collision in collisions
    ]


def _quadrature(energy, f0, processes, row):
    """Shared moments, including the temporal-growth mobility contribution."""
    energy = _numeric_array(energy, "reference energy")
    f0 = _numeric_array(f0, "reference f0")
    if energy.ndim != 1 or len(energy) < 2 or f0.shape != energy.shape:
        raise CorroborationError("invalid reference EEDF arrays")
    if not np.all(np.isfinite(energy)) or not np.all(np.isfinite(f0)) or np.any(energy <= 0.0) or np.any(np.diff(energy) <= 0.0):
        raise CorroborationError("invalid reference EEDF arrays")
    edges = _numeric_array(row["energy_edges_eV"], "reference energy edges")
    if edges.shape != (len(energy) + 1,) or np.any(np.diff(edges) <= 0.0):
        raise CorroborationError("invalid reference energy grid")
    widths = np.diff(edges)
    weights = np.sqrt(energy) * widths
    positive = float(np.sum(weights * np.maximum(f0, 0.0)))
    negative = float(np.sum(weights * np.maximum(-f0, 0.0)))
    if positive <= 0.0 or negative / positive > REFERENCE_EEDF_NEGATIVE_MASS_FLOOR:
        raise CorroborationError("reference EEDF negative mass exceeds numerical floor")
    norm = float(np.sum(weights * f0))
    if not math.isfinite(norm) or norm <= 0.0:
        raise CorroborationError("reference EEDF has zero normalization")
    f0 = f0 / norm
    targets = _numeric_array(row["target_fractions"], "comparison target fractions")
    products = _numeric_array(row["product_fractions"], "comparison product fractions")
    temperature = _finite_float(row["gas_temperature_K"], "gas_temperature_K", positive=True)
    if targets.shape != (len(processes),) or products.shape != (len(processes),):
        raise CorroborationError("comparison fractions have wrong shape")
    rates = []
    powers = []
    momentum_edges = np.zeros_like(edges)
    growth_rate = 0.0
    thermal = Boltzmann / elementary_charge * temperature
    for index, process in enumerate(processes):
        sigma_edges = np.interp(
            edges,
            _numeric_array(process["energy"], "collision energy"),
            _numeric_array(process["sigma"], "collision cross section"),
            left=0.0,
            right=float(process["sigma"][-1]),
        )
        threshold = _finite_float(process["threshold"], "collision threshold")
        if process["kind"] != "elastic":
            sigma_edges[edges <= threshold] = 0.0
        sigma = (sigma_edges[:-1] + sigma_edges[1:]) / 2.0
        momentum_edges += targets[index] * sigma_edges
        rate = float(GAMMA * np.sum(energy * sigma * f0 * widths))
        rates.append(rate)
        if process["kind"] == "ionization":
            growth_rate += targets[index] * rate
        elif process["kind"] == "attachment":
            growth_rate -= targets[index] * rate
        if process["kind"] == "elastic":
            mass_ratio = process.get("mass_ratio")
            if type(mass_ratio) not in (int, float) or mass_ratio <= 0.0:
                raise CorroborationError("elastic mass ratio is missing")
            transfer = 2.0 * edges**2 * sigma_edges * mass_ratio * targets[index]
            transfer[[0, -1]] = 0.0
            kernel = -GAMMA * (
                transfer[1:] * (thermal - widths / 2.0)
                - transfer[:-1] * (thermal + widths / 2.0)
            )
            powers.append(float(kernel.dot(f0)))
        elif process["kind"] == "attachment":
            powers.append(float(targets[index] * np.sum(GAMMA * energy**2 * sigma * f0 * widths)))
        else:
            powers.append(float(targets[index] * rate * threshold))
    derivative = np.r_[0.0, np.diff(f0) / np.diff(energy), 0.0]
    with np.errstate(divide="ignore", invalid="ignore"):
        base = np.where(momentum_edges > 0.0, derivative * edges / momentum_edges, 0.0)
        growth_denominator = momentum_edges + np.where(edges > 0.0, growth_rate / GAMMA / np.sqrt(edges), 0.0)
        corrected = np.where(growth_denominator > 0.0, derivative * edges / growth_denominator, 0.0)
    base[0] = 0.0
    corrected[0] = 0.0
    mobility_without_growth = float(-GAMMA / 3.0 * integrate.simpson(base, x=edges))
    mobility = float(-GAMMA / 3.0 * integrate.simpson(corrected, x=edges))
    powers = np.asarray(powers, dtype=float)
    return {
        "mean_energy_eV": float(np.sum(energy**1.5 * f0 * widths)),
        "mobility_N": mobility,
        "mobility_without_growth_N": mobility_without_growth,
        "growth_rate_m3_s": growth_rate,
        "rates": np.asarray(rates),
        "channel_power": powers,
        "total_power": float(np.sum(powers)),
    }


_BOLOS_WORKER = r'''
import importlib.metadata
import json
import math
import sys

import bolos
import numpy as np
import scipy.integrate
from bolos import grid, parser, solver

if not hasattr(scipy.integrate, "simps"):
    scipy.integrate.simps = scipy.integrate.simpson

request = json.load(open(sys.argv[1]))
collision_grid = grid.LinearGrid(0.0, request["grid"]["max_energy_eV"], request["grid"]["cell_count"])
instance = solver.BoltzmannSolver(collision_grid)
for collision_file in request["collision_files"]:
    with open(collision_file) as stream:
        parsed = parser.parse(stream)
    instance.load_collisions(parsed)
for name, value in request["populations"].items():
    if name not in instance.target:
        raise RuntimeError("canonical target missing from bolos parser: " + name)
    instance.target[name].density = value
instance.kT = request["Tg_K"] * 8.617333262e-5
instance.EN = request["EN_Td"] * solver.TOWNSEND
instance.init()
effective_processes = []
for target, process in instance.iter_all():
    effective_processes.append({
        "kind": "elastic" if process.kind.lower() in ("elastic", "momentum", "effective") else process.kind.lower(),
        "target": target.name,
        "product": process.product,
        "threshold": float(process.threshold),
        "mass_ratio": process.mass_ratio,
        "energy_intervals": np.asarray(process.eps, dtype=float).tolist(),
        "sigma_intervals": np.asarray(process.sigma, dtype=float).tolist(),
    })
effective = {
    "Tg_K": instance.kT / 8.617333262e-5,
    "EN_Td": instance.EN / solver.TOWNSEND,
    "grid": {"max_energy_eV": float(instance.benergy[-1]), "cell_count": int(instance.n)},
    "populations": {name: float(target.density) for name, target in instance.target.items()},
    "target_mass_ratios": {name: target.mass_ratio for name, target in instance.target.items()},
    "processes": effective_processes,
    "ionization_energy_sharing": "equalSharing",
    "electron_growth": "temporal",
    "mobility_definition": "temporal-growth-corrected reduced mobility",
}
try:
    distribution, iterations, residual = instance.converge(
        instance.maxwell(2.0),
        maxn=request["convergence"]["max_iterations"],
        rtol=request["convergence"]["tolerance"],
        full=True,
    )
except solver.ConvergenceError:
    print(json.dumps({
        "termination_status": "iteration_limit",
        "converged": False,
        "iterations": request["convergence"]["max_iterations"],
        "residual": None,
        "tolerance": request["convergence"]["tolerance"],
        "effective_config": effective,
    }))
    raise SystemExit(2)
version = getattr(bolos, "__version__", None) or importlib.metadata.version("bolos")
if any(True for _ in instance.iter_growth()):
    native_mobility = instance.mobility(distribution)
else:
    derivative = np.r_[0.0, np.diff(distribution) / np.diff(instance.cenergy), 0.0]
    integrand = derivative * instance.benergy / instance.sigma_m
    integrand[0] = 0.0
    native_mobility = -(solver.GAMMA / 3.0) * scipy.integrate.simps(
        integrand, x=instance.benergy
    )
print(json.dumps({
    "termination_status": "converged",
    "converged": True,
    "iterations": int(iterations),
    "residual": float(residual),
    "tolerance": request["convergence"]["tolerance"],
    "energy": np.asarray(instance.cenergy).tolist(),
    "f0": np.asarray(distribution).tolist(),
    "native_mobility_N": float(native_mobility),
    "native_mean_energy_eV": float(instance.mean_energy(distribution)),
    "version": version,
    "effective_config": effective,
}))
'''


class SystemdScopeBackend:
    """Execute one solver in a dedicated user scope and kill its whole cgroup."""

    def __init__(self):
        self.systemd_run = shutil.which("systemd-run")
        self.systemctl = shutil.which("systemctl")
        if self.systemd_run is None or self.systemctl is None:
            raise ExecutionBackendUnavailable("systemd user-scope backend commands are unavailable")

    def _control_group(self, unit):
        completed = subprocess.run(
            [self.systemctl, "--user", "show", unit, "--property=ControlGroup", "--value"],
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            timeout=2.0,
            check=False,
        )
        if completed.returncode:
            return None
        value = completed.stdout.strip()
        return Path("/sys/fs/cgroup" + value) if value.startswith("/") else None

    @staticmethod
    def _empty(group):
        if group is None:
            return False
        if not group.exists():
            return True
        try:
            return "populated 0" in (group / "cgroup.events").read_text()
        except OSError:
            return not group.exists()

    def _kill_and_verify(self, unit, group):
        first_failure = None
        cgroup_kill_error = None

        def remember_failure(exc):
            nonlocal first_failure
            if first_failure is None:
                first_failure = exc

        if group is not None and not group.exists():
            return True
        verified = False
        deadline = time.monotonic() + _CLEANUP_GRACE_S
        for attempt in range(_CLEANUP_ATTEMPTS):
            use_cgroup_kill = attempt == 0 and group is not None and group.exists()
            try:
                if use_cgroup_kill:
                    (group / "cgroup.kill").write_text("1")
                else:
                    completed = subprocess.run(
                        [
                            self.systemctl,
                            "--user",
                            "kill",
                            "--kill-whom=all",
                            "--signal=KILL",
                            unit,
                        ],
                        stdin=subprocess.DEVNULL,
                        capture_output=True,
                        text=True,
                        timeout=2.0,
                        check=False,
                    )
                    if completed.returncode:
                        if cgroup_kill_error is not None:
                            error = ExecutionBackendUnavailable(
                                "dedicated cgroup.kill is unavailable (%s) and systemd unit kill failed: %s"
                                % (cgroup_kill_error, completed.stderr.strip())
                            )
                        else:
                            error = ExecutionBackendUnavailable(
                                "unknown solver cgroup and systemd unit kill failed: "
                                + completed.stderr.strip()
                            )
                        remember_failure(error)
                        continue
            except OSError as exc:
                if use_cgroup_kill:
                    if not group.exists():
                        verified = True
                        break
                    cgroup_kill_error = exc
                else:
                    remember_failure(exc)
                continue
            except BaseException as exc:
                remember_failure(exc)
                continue
            if group is None:
                break
            verification_deadline = min(
                deadline,
                time.monotonic() + _CLEANUP_GRACE_S / _CLEANUP_ATTEMPTS,
            )
            verification_failed = False
            while time.monotonic() < verification_deadline:
                try:
                    if self._empty(group):
                        verified = True
                        break
                    time.sleep(0.01)
                except BaseException as exc:
                    remember_failure(exc)
                    verification_failed = True
                    break
            if verified:
                break
            if verification_failed:
                continue
        if group is not None and not verified:
            try:
                verified = self._empty(group)
            except BaseException as exc:
                remember_failure(exc)
        if first_failure is not None:
            raise first_failure
        if group is None:
            return False
        return verified

    def run(self, argv, *, cwd, timeout_s, env=None):
        unit = "eedf-corroborate-" + uuid.uuid4().hex + ".scope"
        payload = list(argv)
        if env is not None:
            payload = ["/usr/bin/env", "-i"] + [
                "%s=%s" % item for item in sorted(env.items())
            ] + payload
        command = [self.systemd_run, "--user", "--scope", "--quiet", "--unit=" + unit, "--"] + payload
        guarded_signals = {signal.SIGINT, signal.SIGTERM}
        previous_signal_mask = signal.pthread_sigmask(signal.SIG_BLOCK, guarded_signals)
        previous_signal_handlers = {
            guarded_signal: signal.getsignal(guarded_signal)
            for guarded_signal in guarded_signals
        }
        deferred_signals = []

        def defer_signal(signum, frame):
            if signum not in deferred_signals:
                deferred_signals.append(signum)

        for guarded_signal in guarded_signals:
            signal.signal(guarded_signal, defer_signal)
        process = None
        interruption = None
        failure = None
        timed_out = False
        stdout = stderr = ""
        try:
            try:
                process = subprocess.Popen(
                    command,
                    cwd=cwd,
                    stdin=subprocess.DEVNULL,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
            except OSError as exc:
                raise ExecutionBackendUnavailable(
                    "systemd user-scope backend could not start"
                ) from exc
            group = None
            signal.pthread_sigmask(signal.SIG_SETMASK, previous_signal_mask)
            try:
                for _ in range(100):
                    group = self._control_group(unit)
                    if group is not None:
                        break
                    if process.poll() is not None:
                        break
                    time.sleep(0.01)
                if group is None and not deferred_signals:
                    raise ExecutionBackendUnavailable(
                        "systemd user-scope backend did not expose a control group"
                    )
                execution_deadline = time.monotonic() + timeout_s
                while not deferred_signals and not timed_out:
                    remaining = execution_deadline - time.monotonic()
                    if remaining <= 0.0:
                        timed_out = True
                        break
                    try:
                        stdout, stderr = process.communicate(timeout=min(0.05, remaining))
                    except subprocess.TimeoutExpired:
                        timed_out = time.monotonic() >= execution_deadline
                    else:
                        break
            except KeyboardInterrupt as exc:
                interruption = exc
            except BaseException as exc:
                failure = exc
        finally:
            signal.pthread_sigmask(signal.SIG_SETMASK, previous_signal_mask)
            if process is not None:
                cleanup_error = None
                try:
                    cleanup_passed = self._kill_and_verify(unit, group)
                except BaseException as exc:
                    cleanup_passed = False
                    cleanup_error = exc
                launcher_cleanup_error = None
                force_launcher_kill = group is None or cleanup_error is not None
                try:
                    launcher_alive = process.poll() is None
                except BaseException as exc:
                    launcher_alive = True
                    force_launcher_kill = True
                    launcher_cleanup_error = exc
                if launcher_alive and not force_launcher_kill:
                    try:
                        stdout_after, stderr_after = process.communicate(timeout=_CLEANUP_GRACE_S)
                        stdout += stdout_after
                        stderr += stderr_after
                    except subprocess.TimeoutExpired:
                        force_launcher_kill = True
                    except BaseException as exc:
                        force_launcher_kill = True
                        launcher_cleanup_error = exc
                if launcher_alive and force_launcher_kill:
                    try:
                        process.kill()
                    except BaseException as exc:
                        if launcher_cleanup_error is None:
                            launcher_cleanup_error = exc
                    try:
                        process.wait()
                    except BaseException as exc:
                        if launcher_cleanup_error is None:
                            launcher_cleanup_error = exc
                    cleanup_passed = False
            signal.pthread_sigmask(signal.SIG_BLOCK, guarded_signals)
            pending_signals = signal.sigpending() & guarded_signals
            for pending_signal in pending_signals:
                signal.sigwait({pending_signal})
                if pending_signal not in deferred_signals:
                    deferred_signals.append(pending_signal)
            for guarded_signal, previous_handler in previous_signal_handlers.items():
                signal.signal(guarded_signal, previous_handler)
            signal.pthread_sigmask(signal.SIG_SETMASK, previous_signal_mask)
            for deferred_signal in deferred_signals:
                signal.raise_signal(deferred_signal)
        if cleanup_error is not None:
            raise cleanup_error
        if launcher_cleanup_error is not None:
            raise launcher_cleanup_error
        if failure is not None:
            raise failure
        if group is None:
            raise ExecutionBackendUnavailable(
                "systemd user-scope backend cleanup is unverified without a control group"
            )
        if interruption is not None:
            raise interruption
        if timed_out:
            raise CorroborationError("solver execution timed out; cleanup_passed=%s" % cleanup_passed)
        if not cleanup_passed:
            raise CorroborationError("solver execution cleanup failed; child survived")
        return {
            "returncode": process.returncode,
            "stdout": stdout,
            "stderr": stderr,
            "cleanup_passed": cleanup_passed,
            "backend": "systemd-user-scope-cgroup.kill",
        }


def _native_quantity(value, unit, name):
    if not isinstance(value, dict) or set(value) != {"unit", "value"} or value["unit"] != unit:
        raise CorroborationError("LoKI effective " + name + " is malformed")
    return _finite_float(float(value["value"]), "LoKI effective " + name)


def _native_constant_states(properties, name):
    if not isinstance(properties, dict) or set(properties) != {"states"}:
        raise CorroborationError("LoKI effective " + name + " is malformed")
    result = {}
    for state, definition in properties["states"].items():
        if (
            not isinstance(definition, dict)
            or set(definition) != {"type", "value"}
            or definition["type"] != "constant"
        ):
            raise CorroborationError("LoKI effective " + name + " is malformed")
        result[state] = _finite_float(float(definition["value"]), "LoKI effective " + name)
    return result


def _materialize_loki_setup(setup, root):
    """Write the pinned LoKI sources and exact canonical deck."""
    for snapshot in setup.sources:
        if snapshot.kind not in {"cross_section", "property"}:
            continue
        path = root / _safe_source_name(snapshot.name)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(snapshot.content)
        if _sha256_bytes(path.read_bytes()) != snapshot.sha256:
            raise CorroborationError("materialized LoKI source fingerprint mismatch: " + snapshot.name)
    setup_path = root / "setup.in"
    setup_path.write_bytes(setup.loki_setup)
    captured = setup_path.read_bytes()
    if _sha256_bytes(captured) != _sha256_bytes(setup.loki_setup):
        raise CorroborationError("materialized LoKI setup fingerprint mismatch")
    return setup_path


def _loki_effective_configuration(native_setup_path):
    """Read the configuration LoKI emits after native parsing and before solving."""
    try:
        rendered = json.loads(Path(native_setup_path).read_bytes())
        working = rendered["workingConditions"]
        kinetics = rendered["electronKinetics"]
        grid = kinetics["numerics"]["energyGrid"]
        nonlinear = kinetics["numerics"]["nonLinearRoutines"]
        reduced_fields = working["reducedField"]
        if not isinstance(reduced_fields, list) or len(reduced_fields) != 1:
            raise KeyError("reducedField")
        effective = {
            "Tg_K": _native_quantity(working["gasTemperature"], "K", "Tg_K"),
            "P_Pa": _native_quantity(working["gasPressure"], "Pa", "P_Pa"),
            "EN_Td": _finite_float(float(reduced_fields[0]), "LoKI effective EN_Td"),
            "excitation_frequency_Hz": _finite_float(
                _native_quantity(
                    working["excitationFrequency"], "Hz", "excitation frequency"
                ),
                "LoKI effective excitation frequency",
            ),
            "eedf_mode": kinetics["eedfType"],
            "electron_electron_collisions": kinetics["includeEECollisions"],
            "ionization_energy_sharing": kinetics["ionizationOperatorType"],
            "electron_growth": kinetics["growthModelType"],
            "grid": {
                "max_energy_eV": _finite_float(
                    float(grid["maxEnergy"]), "LoKI effective max energy", positive=True
                ),
                "cell_count": _strict_int(grid["cellNumber"], "LoKI effective cell count", positive=True),
            },
            "convergence_tolerance": _finite_float(
                float(nonlinear["maxEedfRelError"]),
                "LoKI effective convergence tolerance",
                positive=True,
            ),
            "collision_files": sorted(kinetics["LXCatFiles"]),
            "gas_fractions": {
                gas: _finite_float(float(fraction), "LoKI effective gas fractions")
                for gas, fraction in kinetics["gasProperties"]["fraction"].items()
            },
            "state_populations": _native_constant_states(
                kinetics["stateProperties"]["population"], "state populations"
            ),
        }
    except (KeyError, OSError, TypeError, ValueError) as exc:
        raise CorroborationError("LoKI effective configuration is malformed") from exc
    return effective


def _loki_mismatches_live(setup, effective):
    expected = {
        "Tg_K": setup.gas_temperature_K,
        "P_Pa": setup.pressure_Pa,
        "EN_Td": setup.reduced_field_Td,
        "excitation_frequency_Hz": setup.excitation_frequency_Hz,
        "eedf_mode": setup.eedf_mode,
        "electron_electron_collisions": setup.electron_electron_collisions,
        "ionization_energy_sharing": setup.ionization_energy_sharing,
        "electron_growth": setup.electron_growth,
        "grid": {
            "max_energy_eV": setup.numerical.max_energy_eV,
            "cell_count": setup.numerical.cell_count,
        },
        "convergence_tolerance": setup.numerical.loki_residual_tolerance,
        "collision_files": sorted(
            snapshot.name for snapshot in setup.sources if snapshot.kind == "cross_section"
        ),
        "gas_fractions": dict(setup.gas_fractions),
        "state_populations": dict(setup.state_populations),
    }
    return [
        "%s: expected %r, effective %r" % (name, value, effective.get(name))
        for name, value in expected.items()
        if name not in effective or not _same_effective_value(effective[name], value)
    ]


def _run_loki(setup, solver_config):
    if not isinstance(solver_config, dict) or set(solver_config) != {"kind", "executable"}:
        raise CorroborationError("LoKI adapter requires exactly kind and executable configuration")
    if solver_config["kind"] != "loki":
        raise CorroborationError("only the LoKI live adapter is qualified")
    executable = solver_config["executable"]
    if not isinstance(executable, str) or not executable:
        raise CorroborationError("LoKI executable configuration is required")
    try:
        executable_digest = _sha256_bytes(Path(executable).read_bytes())
    except OSError as exc:
        raise CorroborationError("LoKI executable is unavailable") from exc
    if executable_digest != setup.loki_binary_sha256:
        raise CorroborationError("LoKI executable fingerprint mismatch")
    backend = SystemdScopeBackend()
    with tempfile.TemporaryDirectory(prefix="eedf-corroborate-loki-") as temporary:
        root = Path(temporary)
        _materialize_loki_setup(setup, root)
        execution = backend.run(
            [executable, "setup.in"],
            cwd=root,
            timeout_s=setup.numerical.timeout_s,
            env=dict(setup.loki_environment),
        )
        effective = _loki_effective_configuration(root / "output" / "solve" / "setup.json")
        mismatches = _loki_mismatches_live(setup, effective)
        if mismatches:
            raise CorroborationError("LoKI effective setup mismatch: " + "; ".join(mismatches))
        if execution["returncode"] != 0 or re.search(
            r"did not converge|not converged|results might be incorrect",
            execution["stderr"],
            re.I,
        ):
            raise CorroborationError("LoKI solver failed: " + execution["stderr"][-500:])
        output_roots = list((root / "output" / "solve").glob("reducedField_*"))
        if len(output_roots) != 1:
            raise CorroborationError("LoKI emitted an unexpected result count")
        edges = np.linspace(
            0.0, setup.numerical.max_energy_eV, setup.numerical.cell_count + 1
        )
        try:
            row = parse_output(
                output_roots[0], setup.numerical.normalization_atol, edges
            )
        except (LoKIError, OSError) as exc:
            raise CorroborationError("LoKI emitted malformed solver output: " + str(exc)) from exc
    if row["channels"] != [collision.description for collision in setup.collisions]:
        raise CorroborationError("LoKI live collision identities mismatch the canonical setup")
    residual = abs(row["power_groups"]["relative_balance_percent"]) / 100.0
    if not math.isfinite(residual) or residual > setup.numerical.loki_power_balance_tolerance:
        raise CorroborationError("LoKI live power convergence failed")
    result = {
        "termination_status": "converged",
        "converged": True,
        "residual": residual,
        "tolerance": setup.numerical.loki_power_balance_tolerance,
        "energy": row["energy_eV"],
        "f0": row["f0"],
        "native_mobility_N": row["swarm"]["mobility_N"],
        "native_mean_energy_eV": row["swarm"]["mean_energy_eV"],
        "effective_config": effective,
        "version": setup.loki_binary_sha256,
    }
    return result, execution


def _bolos_request(setup, root):
    collision_files = []
    for snapshot in setup.sources:
        if snapshot.kind != "cross_section":
            continue
        path = root / _safe_source_name(snapshot.name)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(snapshot.content)
        collision_files.append(str(path))
    populations = {}
    for collision, fraction in zip(setup.collisions, setup.target_fractions):
        previous = populations.get(collision.bolos_target)
        if previous is not None and not _same_effective_value(previous, fraction):
            raise CorroborationError(
                "bolos target has inconsistent canonical populations: " + collision.bolos_target
            )
        populations[collision.bolos_target] = fraction
    request = {
        "fingerprint": setup.fingerprint,
        "collision_files": collision_files,
        "Tg_K": setup.gas_temperature_K,
        "EN_Td": setup.reduced_field_Td,
        "populations": populations,
        "grid": {
            "max_energy_eV": setup.numerical.max_energy_eV,
            "cell_count": setup.numerical.cell_count,
        },
        "convergence": {
            "tolerance": setup.numerical.bolos_residual_tolerance,
            "max_iterations": setup.numerical.bolos_max_iterations,
        },
    }
    request_path = root / "request.json"
    request_path.write_text(json.dumps(request, sort_keys=True))
    return request_path


def _validate_bolos_convergence(result, settings):
    required = {"termination_status", "converged", "iterations", "residual", "tolerance"}
    if not isinstance(result, dict) or not required <= set(result):
        raise CorroborationError("bolos convergence evidence is missing or malformed")
    if type(result["converged"]) is not bool:
        raise CorroborationError("bolos converged evidence is not boolean")
    if result["termination_status"] != "converged" or result["converged"] is not True:
        raise CorroborationError("bolos did not converge")
    if type(result["iterations"]) is not int or not 0 < result["iterations"] < settings.bolos_max_iterations:
        raise CorroborationError("bolos iteration count is invalid or reached the limit")
    if type(result["residual"]) not in (int, float) or not math.isfinite(result["residual"]):
        raise CorroborationError("bolos convergence residual is non-finite or malformed")
    if result["tolerance"] != settings.bolos_residual_tolerance:
        raise CorroborationError("bolos convergence tolerance setup mismatch")
    if result["residual"] >= result["tolerance"]:
        raise CorroborationError("bolos convergence residual exceeds tolerance")


def _initialized_bolos_collision(collision, settings):
    """Reproduce bolos's post-init grid cache from the canonical raw collision."""
    boundaries = np.linspace(0.0, settings.max_energy_eV, settings.cell_count + 1)
    shift = 2.0 if collision.kind == "ionization" else 1.0
    shifted = shift * boundaries + collision.threshold_eV
    shifted = np.maximum(shifted, boundaries[0] + 1.0e-9)
    shifted = np.minimum(shifted, boundaries[-1] - 1.0e-9)
    raw_energy = np.asarray(collision.energy_eV)
    raw_sigma = np.asarray(collision.sigma_m2)
    boundary_mask = (boundaries >= shifted[0]) & (boundaries <= shifted[-1])
    raw_mask = (raw_energy >= shifted[0]) & (raw_energy <= shifted[-1])
    nodes = np.unique(np.r_[shifted, boundaries[boundary_mask], raw_energy[raw_mask]])
    if raw_energy[0] > 0.0:
        interpolation_energy = np.r_[0.0, raw_energy, 1.0e8]
        interpolation_sigma = np.r_[raw_sigma[0], raw_sigma, raw_sigma[-1]]
    else:
        interpolation_energy = np.r_[raw_energy, 1.0e8]
        interpolation_sigma = np.r_[raw_sigma, raw_sigma[-1]]
    sigma = np.interp(nodes, interpolation_energy, interpolation_sigma)
    return {
        "energy_intervals": np.c_[nodes[:-1], nodes[1:]].tolist(),
        "sigma_intervals": np.c_[sigma[:-1], sigma[1:]].tolist(),
    }


def _bolos_mismatches(setup, effective):
    if not isinstance(effective, dict):
        return ["effective configuration missing"]
    mismatches = []
    scalar_expected = {
        "Tg_K": setup.gas_temperature_K,
        "EN_Td": setup.reduced_field_Td,
        "grid": {
            "max_energy_eV": setup.numerical.max_energy_eV,
            "cell_count": setup.numerical.cell_count,
        },
        "populations": {
            collision.bolos_target: fraction
            for collision, fraction in zip(setup.collisions, setup.target_fractions)
        },
        "ionization_energy_sharing": setup.ionization_energy_sharing,
        "electron_growth": setup.electron_growth,
        "mobility_definition": setup.mobility_definition,
    }
    for name, expected in scalar_expected.items():
        if name not in effective or not _same_effective_value(effective[name], expected):
            mismatches.append("%s: expected %r, effective %r" % (name, expected, effective.get(name)))
    processes = effective.get("processes")
    if not isinstance(processes, list) or len(processes) != len(setup.collisions):
        mismatches.append("collision count")
        return mismatches
    remaining = list(processes)
    for index, collision in enumerate(setup.collisions):
        expected = {
            "kind": collision.kind,
            "target": collision.bolos_target,
            "product": collision.bolos_product,
            "threshold": collision.threshold_eV,
            "mass_ratio": collision.mass_ratio,
            **_initialized_bolos_collision(collision, setup.numerical),
        }
        match = next(
            (candidate for candidate in remaining if _same_effective_value(candidate, expected)),
            None,
        )
        if match is None:
            mismatches.append("collision[%d]" % index)
        else:
            remaining.remove(match)
    masses = effective.get("target_mass_ratios")
    for collision in setup.collisions:
        if collision.kind == "elastic" and (
            not isinstance(masses, dict)
            or not _same_effective_value(masses.get(collision.bolos_target), collision.mass_ratio)
        ):
            mismatches.append("bolos parsed mass ratio for " + collision.bolos_target)
    return mismatches


def _run_bolos(setup, solver_config):
    if not isinstance(solver_config, dict) or set(solver_config) - {"kind", "interpreter"}:
        raise CorroborationError("bolos adapter accepts only kind and interpreter configuration")
    if solver_config.get("kind") != "bolos":
        raise CorroborationError("only the bolos reference adapter is qualified")
    interpreter = solver_config.get("interpreter")
    if not isinstance(interpreter, str) or not interpreter:
        raise CorroborationError("bolos interpreter configuration is required")
    backend = SystemdScopeBackend()
    with tempfile.TemporaryDirectory(prefix="eedf-corroborate-") as temporary:
        root = Path(temporary)
        request_path = _bolos_request(setup, root)
        execution = backend.run(
            [interpreter, "-c", _BOLOS_WORKER, str(request_path)],
            cwd=root,
            timeout_s=setup.numerical.timeout_s,
        )
        result = None
        try:
            result = json.loads(execution["stdout"].strip().splitlines()[-1])
        except (ValueError, IndexError):
            pass
    if execution["returncode"] != 0:
        if isinstance(result, dict) and result.get("converged") is False:
            raise CorroborationError("bolos did not converge")
        raise CorroborationError("bolos solver failed: " + execution["stderr"][-500:])
    if not isinstance(result, dict):
        raise CorroborationError("bolos emitted malformed solver output")
    _validate_bolos_convergence(result, setup.numerical)
    mismatches = _bolos_mismatches(setup, result.get("effective_config"))
    if mismatches:
        raise CorroborationError("bolos effective setup mismatch: " + "; ".join(mismatches))
    return result, execution


def _flat_quantities(values):
    flat = {
        "mean_energy_eV": float(values["mean_energy_eV"]),
        "mobility_N": float(values["mobility_N"]),
        "total_power": float(values["total_power"]),
    }
    flat.update({"rate[%d]" % index: float(value) for index, value in enumerate(values["rates"])})
    flat.update({"channel_power[%d]" % index: float(value) for index, value in enumerate(values["channel_power"])})
    return flat


def _compare(reference_values, measured, setup):
    comparisons = {}
    outliers = []
    for name, value in measured.items():
        reference = reference_values[name]
        ratio = value / reference if reference != 0.0 else (1.0 if value == 0.0 else None)
        status = "AGREE" if ratio is not None and math.isfinite(ratio) and abs(ratio - 1.0) <= setup.numerical.comparison_rtol else "DISAGREE"
        if status == "DISAGREE":
            outliers.append(name)
        comparisons[name] = {"reference": reference, "value": value, "ratio": ratio, "status": status}
    return comparisons, outliers


def _transport_consistency_ratio(values, native_mobility, setup, solver_name):
    if type(native_mobility) not in (int, float) or not math.isfinite(native_mobility):
        raise CorroborationError(
            solver_name + " native mobility evidence is missing or malformed"
        )
    ratio = values["mobility_N"] / native_mobility
    if (
        not math.isfinite(ratio)
        or abs(ratio - 1.0) > setup.numerical.transport_consistency_rtol
    ):
        raise CorroborationError(
            solver_name
            + " transport-definition mismatch: temporal growth contribution is not matched"
        )
    return ratio


def corroborate_row(artifact, branch, node, u, *, solvers, rtol=None):
    """Run the qualified comparison; only a complete conjunction returns PASS."""
    if rtol is not None and rtol != COMPARISON_RTOL:
        raise CorroborationError("comparison tolerance is frozen in the canonical setup")
    setup = build_canonical_setup(artifact, branch, node, u)
    if not isinstance(solvers, dict) or set(solvers) != {"loki", "bolos"}:
        raise CorroborationError("the qualified comparison requires exactly the live LoKI and bolos adapters")
    row_mapping = {
        "energy_edges_eV": np.linspace(
            0.0, setup.numerical.max_energy_eV, setup.numerical.cell_count + 1
        ),
        "target_fractions": setup.target_fractions,
        "product_fractions": setup.product_fractions,
        "gas_temperature_K": setup.gas_temperature_K,
    }
    processes = _process_dicts(setup.collisions)
    loki, loki_execution = _run_loki(setup, solvers["loki"])
    loki_values = _quadrature(loki.get("energy"), loki.get("f0"), processes, row_mapping)
    loki_mobility_ratio = _transport_consistency_ratio(
        loki_values, loki.get("native_mobility_N"), setup, "LoKI"
    )
    bolos, bolos_execution = _run_bolos(setup, solvers["bolos"])
    bolos_values = _quadrature(bolos.get("energy"), bolos.get("f0"), processes, row_mapping)
    mobility_ratio = _transport_consistency_ratio(
        bolos_values, bolos.get("native_mobility_N"), setup, "bolos"
    )
    comparisons, outliers = _compare(
        _flat_quantities(loki_values), _flat_quantities(bolos_values), setup
    )
    pass_terms = {
        "pinned_input_integrity": True,
        "canonical_loki_setup_agreement": True,
        "canonical_bolos_setup_agreement": True,
        "loki_converged": True,
        "bolos_converged": True,
        "matched_output_definitions": True,
        "numerical_comparison": not outliers,
        "cleanup_passed": (
            loki_execution["cleanup_passed"] is True
            and bolos_execution["cleanup_passed"] is True
        ),
    }
    accepted = all(value is True for value in pass_terms.values())
    return {
        "verdict": "PASS" if accepted else "FAIL",
        "accepted": accepted,
        "setup_fingerprint_version": setup.fingerprint_version,
        "setup_fingerprint": setup.fingerprint,
        "source_hashes": {source.name: source.sha256 for source in setup.sources},
        "effective_configuration_mismatches": {"loki": [], "bolos": []},
        "convergence": {
            "loki": {
                name: loki[name]
                for name in ("termination_status", "converged", "residual", "tolerance")
            },
            "bolos": {name: bolos[name] for name in ("termination_status", "converged", "iterations", "residual", "tolerance")},
        },
        "transport_definition": {
            "name": setup.mobility_definition,
            "growth_rate_m3_s": bolos_values["growth_rate_m3_s"],
            "native_to_shared_mobility_ratio": mobility_ratio,
            "loki_native_to_shared_mobility_ratio": loki_mobility_ratio,
        },
        "loki_live_quantities": _flat_quantities(loki_values),
        "bolos_comparison": comparisons,
        "outliers": outliers,
        "pass_terms": pass_terms,
        "execution_backend": {
            "loki": loki_execution["backend"],
            "bolos": bolos_execution["backend"],
        },
        "loki_version": loki.get("version"),
        "bolos_version": bolos.get("version"),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact")
    parser.add_argument("--branch", default="branch_0")
    parser.add_argument("--u", type=float, required=True)
    parser.add_argument("--node", default="{}", help="JSON composition coordinates")
    parser.add_argument("--bolos", required=True, help="configured bolos interpreter")
    parser.add_argument("--loki", required=True, help="pinned LoKI executable")
    args = parser.parse_args(argv)
    result = corroborate_row(
        args.artifact,
        args.branch,
        json.loads(args.node),
        args.u,
        solvers={
            "loki": {"kind": "loki", "executable": args.loki},
            "bolos": {"kind": "bolos", "interpreter": args.bolos},
        },
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
