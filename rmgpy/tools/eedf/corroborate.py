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

"""Corroborate a stored LoKI-B EEDF from one canonical physical setup."""

from __future__ import annotations

import argparse
import copy
from dataclasses import asdict, dataclass
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time
import uuid

import h5py
import numpy as np
from scipy import integrate
from scipy.constants import Boltzmann, elementary_charge, electron_mass

from rmgpy.tools.eedf.channels import _classic_sections, validate_physical_map
from rmgpy.tools.eedf.integrity import channel_fractions, resolved_properties
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
MASS_RATIO_CONSTANT_RTOL = 5.0e-8
_CLEANUP_GRACE_S = 2.0


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
    bolos_residual_tolerance: float
    bolos_max_iterations: int
    comparison_rtol: float
    transport_consistency_rtol: float
    timeout_s: float


@dataclass(frozen=True)
class StoredRow:
    """The complete stored LoKI result needed for corroboration."""

    energy_eV: tuple[float, ...]
    energy_edges_eV: tuple[float, ...]
    f0: tuple[float, ...]
    target_fractions: tuple[float, ...]
    product_fractions: tuple[float, ...]
    rates: tuple[float, ...]
    channel_power: tuple[float, ...]
    mean_energy_eV: float
    mobility_N: float
    rate_floors: tuple[float, ...]
    field_power: float
    termination_status: str
    converged: bool
    convergence_residual: float
    iteration_count: int
    convergence_tolerance: float
    channels: tuple[str, ...]


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
    mobility_definition: str
    compared_quantities: tuple[str, ...]
    numerical: NumericalSettings
    row: StoredRow
    loki_effective_config: bytes
    loki_effective_sha256: str
    loki_convergence_evidence: bytes
    loki_convergence_sha256: str
    loki_convergence_log: bytes
    loki_convergence_log_sha256: str


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
        table_bytes = (root / "table.h5").read_bytes()
        manifest = json.loads(manifest_bytes)
        spec = json.loads(spec_bytes)
    except (OSError, ValueError) as exc:
        raise CorroborationError("incomplete artifact: " + str(exc)) from exc
    if not isinstance(manifest, dict) or not isinstance(spec, dict):
        raise CorroborationError("artifact metadata must be objects")
    return root, manifest, spec, manifest_bytes, spec_bytes, table_bytes


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


def _h5_scalar(group, name, index):
    if name not in group:
        raise CorroborationError("stored LoKI convergence evidence is missing: " + name)
    value = group[name][index]
    return value.item() if isinstance(value, np.generic) else value


def _stored_row(table_bytes, manifest, branch, coordinates, u):
    if "branches" not in manifest or branch not in manifest["branches"]:
        raise CorroborationError("unknown branch: " + str(branch))
    try:
        axis = _numeric_array(manifest["axes"]["u"], "stored u axis")
    except (KeyError, TypeError, ValueError) as exc:
        raise CorroborationError("invalid stored u axis") from exc
    matches = np.flatnonzero(np.isclose(axis, u, rtol=0.0, atol=1.0e-12))
    if len(matches) != 1:
        raise CorroborationError("requested u is not one unique stored row")
    indices = [int(matches[0])]
    for name, values in manifest["axes"].items():
        if name == "u":
            continue
        if name in {"Tg_K", "P_Pa"}:
            raise CorroborationError(name + " axis cannot override a canonical setup")
        if name not in coordinates:
            raise CorroborationError("missing composition coordinate: " + name)
        values = _numeric_array(values, "composition axis " + name)
        coordinate = _finite_float(coordinates[name], "composition coordinate " + name)
        selected = np.flatnonzero(np.isclose(values, coordinate, rtol=0.0, atol=1.0e-12))
        if len(selected) != 1:
            raise CorroborationError("composition coordinate is not one unique stored value: " + name)
        indices.append(int(selected[0]))
    index = tuple(indices)
    try:
        h5 = h5py.File(io.BytesIO(table_bytes), "r")
        group = h5["branches"][branch]
        energy = tuple(_numeric_array(manifest["energy_eV"], "stored energy grid"))
        edges = tuple(_numeric_array(manifest["energy_edges_eV"], "stored energy edges"))
        channels_raw = group["channels"][index]
        channels = tuple(
            value.decode() if isinstance(value, bytes) else str(value)
            for value in np.atleast_1d(channels_raw)
        )
        converged_raw = _h5_scalar(group, "converged", index)
        if type(converged_raw) not in (bool, np.bool_):
            raise CorroborationError("stored LoKI converged evidence is not boolean")
        termination_raw = _h5_scalar(group, "termination_status", index)
        termination = termination_raw.decode() if isinstance(termination_raw, bytes) else str(termination_raw)
        row = StoredRow(
            energy_eV=energy,
            energy_edges_eV=edges,
            f0=tuple(_numeric_array(group["f0"][index], "stored LoKI f0")),
            target_fractions=tuple(_numeric_array(group["target_fractions"][index], "stored target fractions")),
            product_fractions=tuple(_numeric_array(group["product_fractions"][index], "stored product fractions")),
            rates=tuple(_numeric_array(group["k_ine"][index], "stored rates")),
            channel_power=tuple(_numeric_array(group["channel_power"][index], "stored channel powers")),
            mean_energy_eV=_finite_float(
                _h5_scalar(group["swarm"], "mean_energy_eV", index), "stored mean energy"
            ),
            mobility_N=_finite_float(
                _h5_scalar(group["swarm"], "mobility_N", index), "stored reduced mobility"
            ),
            rate_floors=tuple(_numeric_array(group["rate_floors"][index], "stored rate floors")),
            field_power=_finite_float(
                _h5_scalar(group["power_groups"], "field", index), "stored field power"
            ),
            termination_status=termination,
            converged=converged_raw,
            convergence_residual=_finite_float(
                _h5_scalar(group, "convergence_residual", index), "stored LoKI convergence residual"
            ),
            iteration_count=_strict_int(
                _h5_scalar(group, "iteration_count", index), "stored LoKI iteration count", positive=True
            ),
            convergence_tolerance=_finite_float(
                _h5_scalar(group, "convergence_tolerance", index),
                "stored LoKI convergence tolerance",
                positive=True,
            ),
            channels=channels,
        )
        h5.close()
    except CorroborationError:
        raise
    except (KeyError, OSError, TypeError, ValueError) as exc:
        raise CorroborationError("stored LoKI row is incomplete or malformed: " + str(exc)) from exc
    return row, float(math.exp(axis[indices[0]]))


def _validate_loki_convergence(row, settings):
    if row.termination_status != "converged":
        raise CorroborationError("stored LoKI termination status is not converged")
    if row.converged is not True:
        raise CorroborationError("stored LoKI converged evidence is false")
    if not math.isfinite(row.convergence_residual):
        raise CorroborationError("stored LoKI convergence residual is non-finite")
    if row.iteration_count <= 0:
        raise CorroborationError("stored LoKI iteration count is not positive")
    if row.convergence_tolerance != settings.loki_residual_tolerance:
        raise CorroborationError("stored LoKI convergence tolerance setup mismatch")
    if row.convergence_residual >= row.convergence_tolerance:
        raise CorroborationError("stored LoKI convergence residual exceeds tolerance")


def _effective_quantity(value, unit):
    if not isinstance(value, dict) or set(value) != {"unit", "value"} or value["unit"] != unit:
        return None
    return value["value"]


def _effective_population(value):
    if not isinstance(value, dict) or set(value) != {"type", "value"} or value["type"] != "constant":
        return None
    return value["value"]


def _same_effective_value(actual, expected):
    """Compare serialized native settings without rejecting float roundoff."""
    if type(expected) is float:
        return type(actual) in (int, float) and math.isfinite(actual) and math.isclose(
            actual, expected, rel_tol=1e-14, abs_tol=1e-14
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


def _read_loki_effective(root, manifest):
    entry = manifest.get("loki_effective_config")
    if not isinstance(entry, dict) or set(entry) != {"path", "sha256"}:
        raise CorroborationError("LoKI native effective configuration pin is missing")
    path = root / _safe_source_name(entry["path"])
    try:
        content = path.read_bytes()
    except OSError as exc:
        raise CorroborationError("LoKI native effective configuration is unavailable") from exc
    digest = _sha256_bytes(content)
    if digest != entry["sha256"]:
        raise CorroborationError("LoKI native effective configuration fingerprint mismatch")
    return content, digest


def _read_loki_convergence_evidence(root, manifest):
    entry = manifest.get("loki_convergence_evidence")
    if not isinstance(entry, dict) or set(entry) != {"path", "sha256"}:
        raise CorroborationError("LoKI native convergence evidence pin is missing")
    path = root / _safe_source_name(entry["path"])
    try:
        content = path.read_bytes()
        evidence = json.loads(content)
    except (OSError, ValueError) as exc:
        raise CorroborationError("LoKI native convergence evidence is unavailable or malformed") from exc
    digest = _sha256_bytes(content)
    if digest != entry["sha256"]:
        raise CorroborationError("LoKI native convergence evidence fingerprint mismatch")
    log_entry = manifest.get("loki_convergence_log")
    if not isinstance(log_entry, dict) or set(log_entry) != {"path", "sha256"}:
        raise CorroborationError("LoKI native convergence log pin is missing")
    log_path = root / _safe_source_name(log_entry["path"])
    try:
        log_content = log_path.read_bytes()
    except OSError as exc:
        raise CorroborationError("LoKI native convergence log is unavailable") from exc
    log_digest = _sha256_bytes(log_content)
    if log_digest != log_entry["sha256"]:
        raise CorroborationError("LoKI native convergence log fingerprint mismatch")
    return content, evidence, digest, log_content, log_digest


def _validate_loki_native_evidence(row, evidence, native_log, native_log_digest, spec):
    expected = {
        "termination_status": row.termination_status,
        "converged": row.converged,
        "residual": row.convergence_residual,
        "iteration_count": row.iteration_count,
        "tolerance": row.convergence_tolerance,
    }
    if not isinstance(evidence, dict) or set(evidence) != set(expected) | {"mode", "provenance"}:
        raise CorroborationError("LoKI native convergence evidence is malformed")
    if type(evidence["converged"]) is not bool:
        raise CorroborationError("LoKI native converged evidence is not boolean")
    if not isinstance(evidence["mode"], str) or not evidence["mode"]:
        raise CorroborationError("LoKI native convergence mode is malformed")
    provenance_keys = {
        "base_solver_commit",
        "base_solver_sha256",
        "evidence_log_sha256",
        "instrumentation_source_sha256",
        "instrumented_solver_sha256",
    }
    provenance = evidence["provenance"]
    if not isinstance(provenance, dict) or set(provenance) != provenance_keys:
        raise CorroborationError("LoKI native convergence provenance is malformed")
    if provenance["base_solver_commit"] != spec.get("loki_commit"):
        raise CorroborationError("LoKI native convergence solver commit mismatch")
    binary = spec.get("binary")
    if not isinstance(binary, dict) or provenance["base_solver_sha256"] != binary.get("sha256"):
        raise CorroborationError("LoKI native convergence solver fingerprint mismatch")
    for name in ("instrumentation_source_sha256", "instrumented_solver_sha256"):
        if not isinstance(provenance[name], str) or re.fullmatch(r"[0-9a-f]{64}", provenance[name]) is None:
            raise CorroborationError("LoKI native convergence provenance is malformed")
    if provenance["evidence_log_sha256"] != native_log_digest:
        raise CorroborationError("LoKI native convergence evidence log mismatch")
    pattern = re.compile(
        rb"(?:^|/)T9_CONVERGENCE mode=(\S+) converged=([01]) iterations=(\d+) "
        rb"residual=(\S+) tolerance=(\S+)$",
        re.MULTILINE,
    )
    records = []
    for match in pattern.finditer(native_log):
        records.append(
            {
                "mode": match.group(1).decode("ascii"),
                "termination_status": "converged" if match.group(2) == b"1" else "not_converged",
                "converged": match.group(2) == b"1",
                "iteration_count": int(match.group(3)),
                "residual": _finite_float(float(match.group(4)), "LoKI native log residual"),
                "tolerance": _finite_float(float(match.group(5)), "LoKI native log tolerance", positive=True),
            }
        )
    native = [record for record in records if record["mode"] == evidence["mode"]]
    if len(native) != 1:
        raise CorroborationError("LoKI native convergence log is missing or ambiguous")
    for name in expected:
        if not _same_effective_value(native[0][name], evidence[name]):
            raise CorroborationError("LoKI convergence sidecar differs from native log: " + name)
    mismatches = [name for name, value in expected.items() if not _same_effective_value(evidence[name], value)]
    if mismatches:
        raise CorroborationError("stored LoKI convergence differs from native evidence: " + ", ".join(mismatches))


def _loki_mismatches(setup):
    try:
        effective = json.loads(setup.loki_effective_config)
        working = effective["workingConditions"]
        kinetics = effective["electronKinetics"]
        grid = kinetics["numerics"]["energyGrid"]
        native_gas = effective["nativeGasProperties"]
        native_collisions = tuple(effective["nativeCollisionDescriptions"])
        native_populations = kinetics["stateProperties"]["population"]["states"]
    except (KeyError, TypeError, ValueError) as exc:
        return ["native effective configuration malformed: " + str(exc)]
    expected = {
        "Tg_K": setup.gas_temperature_K,
        "P_Pa": setup.pressure_Pa,
        "EN_Td": setup.reduced_field_Td,
        "ionization_energy_sharing": setup.ionization_energy_sharing,
        "electron_growth": setup.electron_growth,
        "max_energy_eV": setup.numerical.max_energy_eV,
        "cell_count": setup.numerical.cell_count,
        "collision_files": sorted(snapshot.name for snapshot in setup.sources if snapshot.kind == "cross_section"),
        "gas_fractions": dict(setup.gas_fractions),
        "state_populations": dict(setup.state_populations),
        "collision_identities": tuple(collision.description for collision in setup.collisions),
    }
    actual = {
        "Tg_K": _effective_quantity(working.get("gasTemperature"), "K"),
        "P_Pa": _effective_quantity(working.get("gasPressure"), "Pa"),
        "EN_Td": _effective_quantity(working.get("reducedField"), "Td"),
        "ionization_energy_sharing": kinetics.get("ionizationOperatorType"),
        "electron_growth": kinetics.get("growthModelType"),
        "max_energy_eV": grid.get("maxEnergy"),
        "cell_count": grid.get("cellNumber"),
        "collision_files": sorted(kinetics.get("LXCatFiles", [])),
        "gas_fractions": native_gas.get("fraction"),
        "state_populations": {
            name: _effective_population(value) for name, value in native_populations.items()
        },
        "collision_identities": native_collisions,
    }
    mismatches = [name for name in expected if not _same_effective_value(actual[name], expected[name])]
    masses = native_gas.get("mass")
    if not isinstance(masses, dict) or "e" not in masses:
        mismatches.append("elastic mass ratios")
    else:
        for collision in setup.collisions:
            if collision.kind != "elastic":
                continue
            species = collision.description.split(" ", 3)[2].split("(", 1)[0]
            if species not in masses or not math.isclose(
                masses["e"] / masses[species],
                collision.mass_ratio,
                rel_tol=MASS_RATIO_CONSTANT_RTOL,
                abs_tol=0.0,
            ):
                mismatches.append("elastic mass ratio for " + species)
    if setup.row.channels != expected["collision_identities"]:
        mismatches.append("stored collision identities")
    return ["%s: expected %r, effective %r" % (name, expected.get(name), actual.get(name)) for name in mismatches]


def build_canonical_setup(artifact, branch, node, u):
    """Parse and freeze the sole physical setup accepted by both adapters."""
    if node is None:
        node = {}
    if not isinstance(node, dict):
        raise CorroborationError("composition node must be an object")
    forbidden = sorted(set(node) & {"Tg_K", "P_Pa", "EN_Td", "u"})
    if forbidden:
        raise CorroborationError(forbidden[0] + " caller override is forbidden")
    root, manifest, spec, manifest_bytes, spec_bytes, table_bytes = _read_artifact_bytes(artifact)
    snapshots = _snapshot_sources(spec)
    row, reduced_field = _stored_row(table_bytes, manifest, branch, node, u)
    with tempfile.TemporaryDirectory(prefix="eedf-canonical-parse-") as temporary:
        collisions, materialized = _canonical_collisions(spec, manifest, node, snapshots, Path(temporary))
        try:
            _, _, gas_fractions, populations = resolved_properties(materialized, node)
            expected_targets, expected_products = channel_fractions(
                manifest["channel_map"], gas_fractions, populations, node
            )
        except SpecError as exc:
            raise CorroborationError(str(exc)) from exc
    if len(set(collision.description for collision in collisions)) != len(collisions):
        raise CorroborationError("duplicated collision target identity")
    if tuple(expected_targets) != row.target_fractions or tuple(expected_products) != row.product_fractions:
        raise CorroborationError("stored fractions mismatch pinned populations")
    options = spec["solver_options"] if "solver_options" in spec else None
    if not isinstance(options, dict):
        raise CorroborationError("solver_options are missing")
    ionization = options["ionizationOperatorType"] if "ionizationOperatorType" in options else None
    growth = options["growthModelType"] if "growthModelType" in options else None
    if ionization != "equalSharing":
        raise CorroborationError("ionization energy-sharing setup mismatch: equalSharing required")
    if growth != "temporal":
        raise CorroborationError("electron-growth setup mismatch: temporal required")
    if "includeEECollisions" not in options or options["includeEECollisions"] is not False:
        raise CorroborationError("electron-electron collisions are outside common capability")
    if any("<->" in collision.description or collision.kind in {"rotational", "vibrational"} for collision in collisions):
        raise CorroborationError("no second code available for superelastic/rotational populations")
    numerical = options["numerics"] if "numerics" in options else None
    try:
        grid = numerical["energyGrid"]
        loki_tolerance = numerical["nonLinearRoutines"]["maxEedfRelError"]
    except (KeyError, TypeError) as exc:
        raise CorroborationError("frozen numerical settings are incomplete") from exc
    settings = NumericalSettings(
        max_energy_eV=_finite_float(grid["maxEnergy"] if "maxEnergy" in grid else None, "maxEnergy", positive=True),
        cell_count=_strict_int(grid["cellNumber"] if "cellNumber" in grid else None, "cellNumber", positive=True),
        loki_residual_tolerance=_finite_float(loki_tolerance, "LoKI convergence tolerance", positive=True),
        bolos_residual_tolerance=BOLOS_CONVERGENCE_RTOL,
        bolos_max_iterations=BOLOS_MAX_ITERATIONS,
        comparison_rtol=COMPARISON_RTOL,
        transport_consistency_rtol=TRANSPORT_CONSISTENCY_RTOL,
        timeout_s=_finite_float(spec["timeout_s"] if "timeout_s" in spec else None, "timeout_s", positive=True),
    )
    if settings.cell_count < 2 or len(row.energy_eV) != settings.cell_count:
        raise CorroborationError("stored energy grid setup mismatch")
    effective, effective_digest = _read_loki_effective(root, manifest)
    convergence_content, convergence_evidence, convergence_digest, convergence_log, convergence_log_digest = (
        _read_loki_convergence_evidence(root, manifest)
    )
    sources_with_metadata = snapshots + (
        SourceSnapshot("manifest.json", "artifact_metadata", _sha256_bytes(manifest_bytes), manifest_bytes),
        SourceSnapshot("generation_spec.json", "artifact_metadata", _sha256_bytes(spec_bytes), spec_bytes),
        SourceSnapshot("table.h5", "stored_loki_result", _sha256_bytes(table_bytes), table_bytes),
        SourceSnapshot(
            "loki-effective.json",
            "native_effective_config",
            effective_digest,
            effective,
        ),
        SourceSnapshot(
            "loki-convergence.json",
            "native_convergence_evidence",
            convergence_digest,
            convergence_content,
        ),
        SourceSnapshot(
            "loki-convergence.log",
            "native_convergence_log",
            convergence_log_digest,
            convergence_log,
        ),
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
        "mobility_definition": "temporal-growth-corrected reduced mobility",
        "compared_quantities": ["mean_energy_eV", "mobility_N", "rates", "channel_power", "total_power"],
        "numerical": asdict(settings),
        "loki_effective_sha256": effective_digest,
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
        target_fractions=tuple(_numeric_array(expected_targets, "canonical target fractions")),
        product_fractions=tuple(_numeric_array(expected_products, "canonical product fractions")),
        ionization_energy_sharing=ionization,
        electron_growth=growth,
        mobility_definition="temporal-growth-corrected reduced mobility",
        compared_quantities=("mean_energy_eV", "mobility_N", "rates", "channel_power", "total_power"),
        numerical=settings,
        row=row,
        loki_effective_config=effective,
        loki_effective_sha256=effective_digest,
        loki_convergence_evidence=convergence_content,
        loki_convergence_sha256=convergence_digest,
        loki_convergence_log=convergence_log,
        loki_convergence_log_sha256=convergence_log_digest,
    )
    _validate_loki_convergence(setup.row, setup.numerical)
    _validate_loki_native_evidence(
        setup.row, convergence_evidence, setup.loki_convergence_log,
        setup.loki_convergence_log_sha256, spec,
    )
    mismatches = _loki_mismatches(setup)
    if mismatches:
        raise CorroborationError("LoKI effective setup mismatch: " + "; ".join(mismatches))
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
native_processes = []
for collision_file in request["collision_files"]:
    with open(collision_file) as stream:
        parsed = parser.parse(stream)
    native_processes.extend(parsed)
    instance.load_collisions(parsed)
for name, value in request["gas_fractions"].items():
    if name not in instance.target:
        raise RuntimeError("canonical target missing from bolos parser: " + name)
    instance.target[name].density = value
instance.kT = request["Tg_K"] * 8.617333262e-5
instance.EN = request["EN_Td"] * solver.TOWNSEND
instance.init()
effective_processes = []
for process in native_processes:
    effective_processes.append({
        "kind": "elastic" if process["kind"].lower() in ("elastic", "momentum", "effective") else process["kind"].lower(),
        "target": process.get("target"),
        "product": process.get("product"),
        "threshold": float(process.get("threshold", 0.0)),
        "mass_ratio": process.get("mass_ratio"),
        "energy": [float(row[0]) for row in process["data"]],
        "sigma": [float(row[1]) for row in process["data"]],
    })
effective = {
    "Tg_K": instance.kT / 8.617333262e-5,
    "EN_Td": instance.EN / solver.TOWNSEND,
    "grid": {"max_energy_eV": float(instance.benergy[-1]), "cell_count": int(instance.n)},
    "gas_fractions": {name: float(target.density) for name, target in instance.target.items()},
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
        if group is None or not group.exists():
            return True
        try:
            return "populated 0" in (group / "cgroup.events").read_text()
        except OSError:
            return not group.exists()

    def _kill_and_verify(self, group):
        if group is not None and group.exists():
            try:
                (group / "cgroup.kill").write_text("1")
            except OSError as exc:
                raise ExecutionBackendUnavailable("dedicated cgroup.kill is unavailable: " + str(exc)) from exc
        deadline = time.monotonic() + _CLEANUP_GRACE_S
        while time.monotonic() < deadline:
            if self._empty(group):
                return True
            time.sleep(0.01)
        return self._empty(group)

    def run(self, argv, *, cwd, timeout_s):
        unit = "eedf-corroborate-" + uuid.uuid4().hex + ".scope"
        command = [self.systemd_run, "--user", "--scope", "--quiet", "--unit=" + unit, "--"] + list(argv)
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
            raise ExecutionBackendUnavailable("systemd user-scope backend could not start") from exc
        group = None
        for _ in range(100):
            group = self._control_group(unit)
            if group is not None:
                break
            if process.poll() is not None:
                break
            time.sleep(0.01)
        if group is None and process.poll() is None:
            process.kill()
            process.wait()
            raise ExecutionBackendUnavailable("systemd user-scope backend did not expose a control group")
        interrupted = False
        timed_out = False
        try:
            stdout, stderr = process.communicate(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            timed_out = True
            stdout = stderr = ""
        except KeyboardInterrupt:
            interrupted = True
            stdout = stderr = ""
        cleanup_passed = self._kill_and_verify(group)
        if process.poll() is None:
            try:
                stdout_after, stderr_after = process.communicate(timeout=_CLEANUP_GRACE_S)
                stdout += stdout_after
                stderr += stderr_after
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
                cleanup_passed = False
        if interrupted:
            raise KeyboardInterrupt
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


def _bolos_request(setup, root):
    collision_files = []
    for snapshot in setup.sources:
        if snapshot.kind != "cross_section":
            continue
        path = root / _safe_source_name(snapshot.name)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(snapshot.content)
        collision_files.append(str(path))
    request = {
        "fingerprint": setup.fingerprint,
        "collision_files": collision_files,
        "Tg_K": setup.gas_temperature_K,
        "EN_Td": setup.reduced_field_Td,
        "gas_fractions": dict(setup.gas_fractions),
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
        "gas_fractions": dict(setup.gas_fractions),
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
    for index, (collision, actual) in enumerate(zip(setup.collisions, processes)):
        expected = {
            "kind": collision.kind,
            "target": collision.bolos_target,
            "product": collision.bolos_product,
            "threshold": collision.threshold_eV,
            "mass_ratio": collision.mass_ratio,
            "energy": list(collision.energy_eV),
            "sigma": list(collision.sigma_m2),
        }
        if not _same_effective_value(actual, expected):
            mismatches.append("collision[%d]" % index)
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


def _stored_quantities(row):
    values = {
        "mean_energy_eV": row.mean_energy_eV,
        "mobility_N": row.mobility_N,
        "total_power": float(sum(row.channel_power)),
    }
    values.update({"rate[%d]" % index: value for index, value in enumerate(row.rates)})
    values.update({"channel_power[%d]" % index: value for index, value in enumerate(row.channel_power)})
    if not all(math.isfinite(value) for value in values.values()):
        raise CorroborationError("stored comparison quantity is non-finite")
    return values


def _compare(stored, measured, setup):
    comparisons = {}
    outliers = []
    power_floor = abs(setup.row.field_power) * 1.0e-4
    for name, value in measured.items():
        reference = stored[name]
        if name.startswith("rate["):
            floor = setup.row.rate_floors[int(name[5:-1])]
        elif name.startswith("channel_power[") or name == "total_power":
            floor = power_floor
        else:
            floor = None
        below = floor is not None and abs(reference) < floor and abs(value) < floor
        if below:
            ratio = None
            status = "BELOW_FLOOR"
        else:
            ratio = value / reference if reference != 0.0 else (1.0 if value == 0.0 else None)
            status = "AGREE" if ratio is not None and math.isfinite(ratio) and abs(ratio - 1.0) <= setup.numerical.comparison_rtol else "DISAGREE"
        if status == "DISAGREE":
            outliers.append(name)
        comparisons[name] = {"reference": reference, "value": value, "ratio": ratio, "status": status, "floor": floor}
    return comparisons, outliers


def corroborate_row(artifact, branch, node, u, *, solvers, rtol=None):
    """Run the qualified comparison; only a complete conjunction returns PASS."""
    if rtol is not None and rtol != COMPARISON_RTOL:
        raise CorroborationError("comparison tolerance is frozen in the canonical setup")
    setup = build_canonical_setup(artifact, branch, node, u)
    if not isinstance(solvers, dict) or set(solvers) != {"bolos"}:
        raise CorroborationError("the qualified comparison requires exactly the bolos adapter")
    row_mapping = {
        "energy_edges_eV": setup.row.energy_edges_eV,
        "target_fractions": setup.target_fractions,
        "product_fractions": setup.product_fractions,
        "gas_temperature_K": setup.gas_temperature_K,
    }
    processes = _process_dicts(setup.collisions)
    loki = _quadrature(setup.row.energy_eV, setup.row.f0, processes, row_mapping)
    stored = _stored_quantities(setup.row)
    loki_comparison, loki_outliers = _compare(stored, _flat_quantities(loki), setup)
    if loki_outliers:
        raise CorroborationError("stored LoKI quantities mismatch canonical definitions: " + ", ".join(loki_outliers))
    bolos, execution = _run_bolos(setup, solvers["bolos"])
    bolos_values = _quadrature(bolos.get("energy"), bolos.get("f0"), processes, row_mapping)
    native_mobility = bolos.get("native_mobility_N")
    if type(native_mobility) not in (int, float) or not math.isfinite(native_mobility):
        raise CorroborationError("bolos native mobility evidence is missing or malformed")
    mobility_ratio = bolos_values["mobility_N"] / native_mobility
    if not math.isfinite(mobility_ratio) or abs(mobility_ratio - 1.0) > setup.numerical.transport_consistency_rtol:
        raise CorroborationError("transport-definition mismatch: temporal growth contribution is not matched")
    comparisons, outliers = _compare(stored, _flat_quantities(bolos_values), setup)
    pass_terms = {
        "pinned_input_integrity": True,
        "canonical_loki_setup_agreement": True,
        "canonical_bolos_setup_agreement": True,
        "loki_converged": True,
        "bolos_converged": True,
        "matched_output_definitions": True,
        "numerical_comparison": not outliers,
        "cleanup_passed": execution["cleanup_passed"] is True,
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
                "termination_status": setup.row.termination_status,
                "converged": setup.row.converged,
                "iterations": setup.row.iteration_count,
                "residual": setup.row.convergence_residual,
                "tolerance": setup.row.convergence_tolerance,
            },
            "bolos": {name: bolos[name] for name in ("termination_status", "converged", "iterations", "residual", "tolerance")},
        },
        "transport_definition": {
            "name": setup.mobility_definition,
            "growth_rate_m3_s": bolos_values["growth_rate_m3_s"],
            "native_to_shared_mobility_ratio": mobility_ratio,
        },
        "loki_definition_check": loki_comparison,
        "bolos_comparison": comparisons,
        "outliers": outliers,
        "pass_terms": pass_terms,
        "execution_backend": execution["backend"],
        "bolos_version": bolos.get("version"),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact")
    parser.add_argument("--branch", default="branch_0")
    parser.add_argument("--u", type=float, required=True)
    parser.add_argument("--node", default="{}", help="JSON composition coordinates")
    parser.add_argument("--bolos", required=True, help="configured bolos interpreter")
    args = parser.parse_args(argv)
    result = corroborate_row(
        args.artifact,
        args.branch,
        json.loads(args.node),
        args.u,
        solvers={"bolos": {"kind": "bolos", "interpreter": args.bolos}},
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
