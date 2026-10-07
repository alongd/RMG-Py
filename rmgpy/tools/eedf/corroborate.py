#!/usr/bin/env python3
"""Independent two-term corroboration of a stored EEDF table row.

Reference solvers are deliberately subprocesses.  The solver interpreters and
their package names are configuration, never repository inputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import tempfile

import h5py
import numpy as np
from scipy.constants import Boltzmann, elementary_charge, electron_mass

from rmgpy.tools.eedf.channels import _classic_sections, validate_physical_map
from rmgpy.tools.eedf.integrity import resolved_properties
from rmgpy.tools.eedf.moments import distribution_moments
from rmgpy.tools.eedf.schema import SpecError


class CorroborationError(ValueError):
    """The artifact or requested comparison is not admissible."""


class Unavailable(CorroborationError):
    """A reference solver cannot represent or converge the requested row."""


GAMMA = np.sqrt(2 * elementary_charge / electron_mass)
REFERENCE_EEDF_NEGATIVE_MASS_FLOOR = 1e-12
_TIMEOUT_CLEANUP_GRACE_S = 0.25
_WORKER = r"""
import importlib.metadata
import json
import sys

import bolos
import numpy as np
import scipy.integrate
from bolos import grid, parser, solver

cfg = json.load(open(sys.argv[1]))
kind = cfg["kind"]
if kind == "bolos":
    if not hasattr(scipy.integrate, "simps"):
        scipy.integrate.simps = scipy.integrate.simpson

    collision_grid = grid.LinearGrid(0., cfg["emax"], cfg["ncells"])
    boltzmann_solver = solver.BoltzmannSolver(collision_grid)
    for collision_file in cfg["collision_files"]:
        with open(collision_file) as collision_stream:
            boltzmann_solver.load_collisions(parser.parse(collision_stream))
    boltzmann_solver.target["Ar"].density = 1.
    boltzmann_solver.kT = cfg["Tg"] * 8.617333262e-5
    boltzmann_solver.EN = cfg["EN"] * solver.TOWNSEND
    boltzmann_solver.init()
    try:
        distribution = boltzmann_solver.converge(
            boltzmann_solver.maxwell(2.), maxn=500, rtol=1e-7
        )
    except Exception as exc:
        print(json.dumps({"status": "unconverged", "reason": str(exc)}))
        raise SystemExit
    version = getattr(bolos, "__version__", None) or importlib.metadata.version("bolos")
    print(json.dumps({
        "status": "converged",
        "energy": np.asarray(boltzmann_solver.cenergy).tolist(),
        "f0": np.asarray(distribution).tolist(),
        "version": version,
    }))
"""


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _artifact_files(artifact):
    root = Path(artifact)
    if root.is_file():
        root = root.parent
    try:
        manifest = json.loads((root / "manifest.json").read_text())
        spec = json.loads((root / "generation_spec.json").read_text())
    except OSError as exc:
        raise CorroborationError("incomplete artifact: " + str(exc)) from exc
    try:
        input_files = spec["input_files"]
    except KeyError as exc:
        raise CorroborationError("generation spec is missing input_files") from exc
    if not isinstance(input_files, dict) or not input_files:
        raise CorroborationError("generation spec input_files is empty or invalid")
    for name, entry in input_files.items():
        path = Path(entry["path"])
        if not path.exists() or _sha256(path) != entry["sha256"]:
            raise CorroborationError("pinned input SHA-256 mismatch: " + name)
    return root, manifest, spec


def _row(root, manifest, branch, node, u):
    if branch not in manifest.get("branches", []):
        raise CorroborationError("unknown branch: " + str(branch))
    axis = np.asarray(manifest["axes"]["u"], float)
    index = int(np.argmin(abs(axis - u)))
    if not np.isclose(axis[index], u, rtol=0, atol=1e-8):
        raise CorroborationError("requested u is not a stored row")
    names = [name for name in manifest["axes"] if name != "u"]
    node = node or {}
    indices = [index]
    for name in names:
        values = np.asarray(manifest["axes"][name], float)
        if name not in node:
            raise CorroborationError("missing composition node: " + name)
        j = int(np.argmin(abs(values - node[name])))
        if not np.isclose(values[j], node[name], rtol=0, atol=1e-8):
            raise CorroborationError("node is not stored")
        indices.append(j)
    with h5py.File(root / "table.h5", "r") as h5:
        group = h5["branches"][branch]
        out = {
            "energy_eV": np.asarray(manifest["energy_eV"]),
            "energy_edges_eV": np.asarray(manifest["energy_edges_eV"]),
        }
        for key, value in group.items():
            if isinstance(value, h5py.Group):
                out[key] = {k: np.asarray(v[tuple(indices)]) for k, v in value.items()}
            elif key not in ("energy_eV", "energy_edges_eV"):
                out[key] = np.asarray(value[tuple(indices)])
        out["EN_Td"] = float(np.exp(axis[index]))
    return out


def _processes(spec, manifest, node):
    """Rebuild collision arrays from every pinned file, never from the manifest."""
    channels = manifest.get("channel_map") or []
    try:
        validate_physical_map(channels, spec, node)
    except SpecError as exc:
        raise CorroborationError(str(exc)) from exc

    physical = {}
    for entry in spec["input_files"].values():
        if entry.get("kind") != "cross_section":
            continue
        try:
            sections = _classic_sections(entry["path"], entry["sha256"])
        except SpecError as exc:
            raise CorroborationError(str(exc)) from exc
        for section in sections:
            if section["description"] in physical:
                raise CorroborationError("duplicate physical collision")
            physical[section["description"]] = section
    if set(physical) != {channel["description"] for channel in channels}:
        raise CorroborationError(
            "channel map does not cover physical collision identities"
        )

    processes = []
    for channel in channels:
        section = physical[channel["description"]]
        description = section["description"]
        processes.append(
            {
                "description": description,
                "kind": section["kind"],
                "product": description.split(" -> ")[-1]
                .rsplit(", ", 1)[0]
                .replace("e + ", "")
                .strip(),
                "threshold": section["threshold_eV"],
                "energy": section["cross_section"]["energy_eV"],
                "sigma": section["cross_section"]["sigma_m2"],
                "mass_ratio": channel.get("mass_ratio"),
                "opb_eV": channel.get("opb_eV"),
                "statistical_weight_ratio": channel.get("statistical_weight_ratio"),
            }
        )
    return processes


def _unavailable_reason(kind, spec, node, processes):
    """Name any physical or numerical setup an adapter would substitute."""
    if kind == "cantera":
        return "grid not controllable"
    if kind != "bolos":
        return "unknown reference solver adapter"
    if any(
        "<->" in process["description"]
        or process["kind"] in ("rotational", "vibrational")
        for process in processes
    ):
        return "no second code available for superelastic/rotational populations"
    arm = spec.get("arm", {})
    if (
        arm.get("gases") != ["Ar"]
        or arm.get("additive") is not None
        or arm.get("feed_partial_pressure_Pa") != 0
    ):
        return "reference adapter supports pure Ar only"
    try:
        _, states, fractions, populations = resolved_properties(spec, node)
    except SpecError as exc:
        return "physical setup cannot be resolved exactly: " + str(exc)
    if fractions != {"Ar": 1.0}:
        return "reference adapter supports pure Ar only"
    if populations.get("Ar(1S0)") != 1.0 or any(
        value != 0.0 for state, value in populations.items() if state != "Ar(1S0)"
    ):
        return "excited-state population is not representable exactly"
    state_energies = {
        entry.split(" = ", 1)[0]: float(entry.split(" = ", 1)[1])
        for entry in states.get("energy", [])
    }
    if state_energies and state_energies != {"Ar(1S0)": 0.0}:
        return "nonzero or excited-state energies are not representable exactly"
    options = spec.get("solver_options", {})
    if options.get("includeEECollisions") is not False:
        return "electron-electron collision setting is not representable exactly"
    if spec.get("working_conditions", {}).get("excitationFrequency") != 0:
        return "non-DC excitation is not representable exactly"
    if options.get("eedfType") != "boltzmann":
        return "EEDF operator is not the two-term Boltzmann operator"
    if options.get("growthModelType") != "temporal":
        return "growth operator is not representable exactly"
    if not any(process["kind"] == "ionization" for process in processes):
        return None
    ionization_operator = options.get("ionizationOperatorType")
    if ionization_operator == "equalSharing":
        return None
    if ionization_operator == "usingSDCS":
        return "ionization operator usingSDCS is not representable exactly"
    return "ionization operator is not a validated exact match: " + str(
        ionization_operator
    )


def _physical_condition(spec, manifest, node, name):
    """Resolve a row condition exactly as LoKI's setup builder does."""
    state = {
        key: envelope.get("reference")
        for key, envelope in spec.get("envelopes", {}).items()
    }
    state.update(node)
    value = state.get(name, spec.get(name, manifest.get("row_inputs", {}).get(name)))
    if value is None or not np.isfinite(value):
        raise CorroborationError("row is missing physical condition: " + name)
    return float(value)


def _quadrature(energy, f0, processes, row):
    """Shared cell quadrature; elastic power is net transfer including thermal gain."""
    energy, f0 = np.asarray(energy, float), np.asarray(f0, float)
    if (
        energy.ndim != 1
        or len(energy) < 2
        or f0.shape != energy.shape
        or not np.all(np.isfinite(energy))
        or not np.all(np.isfinite(f0))
        or np.any(energy < 0)
        or np.any(np.diff(energy) <= 0)
    ):
        raise CorroborationError("invalid reference EEDF arrays")
    stored_energy = np.asarray(row["energy_eV"], float)
    if len(energy) == len(stored_energy) and np.allclose(
        energy, stored_energy, rtol=0, atol=1e-12
    ):
        edges = np.asarray(row["energy_edges_eV"], float)
    else:
        edges = np.r_[
            energy[0] - (energy[1] - energy[0]) / 2,
            (energy[:-1] + energy[1:]) / 2,
            energy[-1] + (energy[-1] - energy[-2]) / 2,
        ]
    widths = np.diff(edges)
    if np.any(widths <= 0):
        raise CorroborationError("invalid reference energy grid")
    mass_weights = np.sqrt(energy) * widths
    positive_mass = float(np.sum(mass_weights * np.maximum(f0, 0.0)))
    negative_mass = float(np.sum(mass_weights * np.maximum(-f0, 0.0)))
    negative_fraction = negative_mass / positive_mass if positive_mass > 0 else np.inf
    if negative_fraction > REFERENCE_EEDF_NEGATIVE_MASS_FLOOR:
        raise CorroborationError(
            "reference EEDF negative mass %.6g exceeds numerical floor %.6g"
            % (negative_fraction, REFERENCE_EEDF_NEGATIVE_MASS_FLOOR)
        )
    norm = np.sum(mass_weights * f0)
    if not np.isfinite(norm) or norm <= 0:
        raise CorroborationError("reference EEDF has zero normalization")
    f0 /= norm
    mean = float(np.sum(energy**1.5 * f0 * widths))
    try:
        targets = np.asarray(row["target_fractions"], float)
        products = np.asarray(row["product_fractions"], float)
        gas_temperature = float(row["gas_temperature_K"])
    except KeyError as exc:
        raise CorroborationError(
            "stored row is missing required comparison quantity: " + str(exc.args[0])
        ) from exc
    if targets.shape != (len(processes),) or not np.all(np.isfinite(targets)):
        raise CorroborationError("invalid stored target_fractions")
    if products.shape != (len(processes),) or not np.all(np.isfinite(products)):
        raise CorroborationError("invalid stored product_fractions")
    if not np.isfinite(gas_temperature):
        raise CorroborationError("invalid stored gas_temperature_K")
    rates, powers, sigtot = [], [], np.zeros_like(energy)
    thermal = Boltzmann / elementary_charge * gas_temperature
    for index, process in enumerate(processes):
        try:
            section_energy = process["energy"]
            section_sigma = process["sigma"]
            threshold = process["threshold"]
        except KeyError as exc:
            raise CorroborationError(
                "collision process is missing required comparison quantity: "
                + str(exc.args[0])
            ) from exc
        sigma_edges = np.interp(
            edges, section_energy, section_sigma, left=0.0, right=0.0
        )
        if process["kind"] != "elastic":
            sigma_edges[edges <= threshold] = 0.0
        sigma = (sigma_edges[:-1] + sigma_edges[1:]) / 2
        sigtot += targets[index] * sigma
        rate = GAMMA * np.sum(energy * sigma * f0 * widths)
        rates.append(rate)
        if process["kind"] == "elastic":
            mass_ratio = process.get("mass_ratio")
            if not mass_ratio:
                raise CorroborationError("elastic mass ratio is missing")
            transfer = 2 * edges**2 * sigma_edges * mass_ratio * targets[index]
            transfer[[0, -1]] = 0.0
            kernel = -GAMMA * (
                transfer[1:] * (thermal - widths / 2)
                - transfer[:-1] * (thermal + widths / 2)
            )
            powers.append(float(kernel.dot(f0)))
        elif process["kind"] == "attachment":
            powers.append(
                float(targets[index] * np.sum(GAMMA * energy**2 * sigma * f0 * widths))
            )
        else:
            lower = min(int(threshold / widths[0]), len(energy))
            effective_threshold = edges[lower] if lower < len(energy) else threshold
            powers.append(
                float(
                    (targets[index] * rate - products[index] * 0.0)
                    * effective_threshold
                )
            )
    dfde = np.gradient(f0, energy)
    with np.errstate(divide="ignore", invalid="ignore"):
        mobility = (
            -GAMMA
            / 3
            * np.trapz(np.where(sigtot > 0, energy / sigtot * dfde, 0.0), energy)
        )
    powers = np.asarray(powers)
    return {
        "mean_energy_eV": mean,
        "mobility_N": float(mobility),
        "rates": np.asarray(rates),
        "channel_power": powers,
        "total_power": float(powers.sum()),
    }


def _loki_quadrature(row, processes, ionization_operator="usingSDCS"):
    """Recompute LoKI moments with LoKI's recorded ionization operator."""
    channels = []
    for process in processes:
        channel = dict(process)
        channel["threshold_eV"] = process["threshold"]
        channel["cross_section"] = {
            "energy_eV": process["energy"],
            "sigma_m2": process["sigma"],
        }
        channels.append(channel)
    moments = distribution_moments(row, channels)
    transport = _quadrature(row["energy_eV"], row["f0"], processes, row)
    if ionization_operator == "equalSharing":
        return transport
    powers = np.asarray(moments["channel_power"], float)
    return {
        "mean_energy_eV": moments["mean_energy_eV"],
        "mobility_N": transport["mobility_N"],
        "rates": np.asarray(moments["k_ine"], float),
        "channel_power": powers,
        "total_power": float(powers.sum()),
    }


def _solver_spec(name, value):
    if isinstance(value, str):
        return {"interpreter": value, "kind": name}
    result = dict(value)
    result.setdefault("kind", name)
    return result


def _linux_descendants(pid):
    """Return live descendants while Linux still exposes their parent links."""
    children_file = Path("/proc") / str(pid) / "task" / str(pid) / "children"
    if not children_file.exists():
        return None
    descendants, pending = set(), [pid]
    while pending:
        parent = pending.pop()
        path = Path("/proc") / str(parent) / "task" / str(parent) / "children"
        try:
            children = [int(value) for value in path.read_text().split()]
        except (FileNotFoundError, ProcessLookupError):
            continue
        new_children = [child for child in children if child not in descendants]
        descendants.update(new_children)
        pending.extend(new_children)
    return descendants


def _terminate_reference_tree(process):
    """Best-effort bounded cleanup, including Linux descendants in new sessions."""
    try:
        os.killpg(process.pid, signal.SIGSTOP)
    except ProcessLookupError:
        pass
    descendants = _linux_descendants(process.pid)
    if descendants is not None:
        for pid in descendants:
            try:
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    for stream in (process.stdout, process.stderr):
        if stream is not None:
            stream.close()
    try:
        process.wait(timeout=_TIMEOUT_CLEANUP_GRACE_S)
    except subprocess.TimeoutExpired:
        pass
    return descendants is not None


def _run_reference(config, solver):
    for key in ("continuation_from", "continuation_steps"):
        if key in solver:
            config[key] = solver[key]
    config.update({"kind": solver.get("kind", "cantera")})
    with tempfile.TemporaryDirectory(prefix="eedf-corroborate-") as tmp:
        path = Path(tmp) / "request.json"
        path.write_text(json.dumps(config))
        interpreter = solver.get("interpreter") or solver.get("python")
        if not interpreter:
            raise Unavailable("solver interpreter is required")
        timeout = float(solver.get("timeout_s", config.get("timeout_s", 60.0)))
        try:
            process = subprocess.Popen(
                [interpreter, "-c", _WORKER + "\n", str(path)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                start_new_session=True,
            )
        except OSError as exc:
            raise Unavailable("reference solver could not start: " + str(exc)) from exc
        try:
            stdout, stderr = process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            enumerated_descendants = _terminate_reference_tree(process)
            limitation = (
                ""
                if enumerated_descendants
                else "; escaped descendants cannot be enumerated on this platform"
            )
            raise Unavailable(
                "reference solver timeout after %.6g s%s" % (timeout, limitation)
            ) from exc
        if process.returncode:
            raise Unavailable("reference solver failed: " + stderr[-500:])
        try:
            result = json.loads(stdout.strip().splitlines()[-1])
        except (ValueError, IndexError) as exc:
            raise Unavailable("invalid solver output") from exc
    if not isinstance(result, dict):
        raise Unavailable("invalid solver output")
    if result.get("status") != "converged":
        raise Unavailable("unconverged: " + result.get("reason", "unknown"))
    return result


def _flat_quantities(values):
    flat = {
        "mean_energy_eV": float(values["mean_energy_eV"]),
        "mobility_N": float(values["mobility_N"]),
        "total_power": float(values["total_power"]),
    }
    flat.update(
        {"rate[%d]" % i: float(value) for i, value in enumerate(values["rates"])}
    )
    flat.update(
        {
            "channel_power[%d]" % i: float(value)
            for i, value in enumerate(values["channel_power"])
        }
    )
    return flat


def _stored_quantities(row, count):
    try:
        stored = {
            "mean_energy_eV": float(row["swarm"]["mean_energy_eV"]),
            "mobility_N": float(row["swarm"]["mobility_N"]),
            "total_power": float(np.sum(row["channel_power"])),
        }
        rates = np.asarray(row["k_ine"], float)
        powers = np.asarray(row["channel_power"], float)
    except (KeyError, TypeError, ValueError) as exc:
        raise CorroborationError(
            "stored row is missing a corroboration quantity"
        ) from exc
    if rates.shape != (count,) or powers.shape != (count,):
        raise CorroborationError("stored row corroboration quantity shape mismatch")
    stored.update({"rate[%d]" % i: float(value) for i, value in enumerate(rates)})
    stored.update(
        {"channel_power[%d]" % i: float(value) for i, value in enumerate(powers)}
    )
    if not all(np.isfinite(value) for value in stored.values()):
        raise CorroborationError("stored row has nonfinite corroboration quantity")
    return stored


def _floors(row, manifest, count):
    try:
        rate_floors = np.asarray(row["rate_floors"], float)
        field_power = abs(float(row["power_groups"]["field"]))
        absolute_power_share = float(manifest["floors"]["absolute_power_share"])
    except KeyError as exc:
        raise CorroborationError(
            "stored comparison metadata is missing: " + str(exc.args[0])
        ) from exc
    if rate_floors.shape != (count,) or not np.all(np.isfinite(rate_floors)):
        raise CorroborationError("stored row rate-floor shape/value mismatch")
    if not np.isfinite(field_power) or not np.isfinite(absolute_power_share):
        raise CorroborationError("stored row has nonfinite power-floor metadata")
    if np.any(rate_floors < 0) or absolute_power_share < 0:
        raise CorroborationError("stored row has negative comparison floor")
    power_floor = absolute_power_share * field_power
    floors = {"rate[%d]" % i: float(value) for i, value in enumerate(rate_floors)}
    floors.update({"channel_power[%d]" % i: power_floor for i in range(count)})
    floors["total_power"] = power_floor
    return floors


def _compare_quantities(stored, measured, floors, below_floor, rtol):
    comparisons, ratios, outliers = {}, {}, []
    for key, value in measured.items():
        reference = stored[key]
        floor = floors.get(key)
        if below_floor[key]:
            ratio, status = None, "BELOW_FLOOR"
        else:
            ratio = (
                float(value / reference) if reference else (1.0 if value == 0 else None)
            )
            status = (
                "AGREE"
                if ratio is not None and np.isfinite(ratio) and abs(ratio - 1) <= rtol
                else "DISAGREE"
            )
            if status == "DISAGREE":
                outliers.append(key)
        ratios[key] = ratio
        comparisons[key] = {
            "status": status,
            "ratio": ratio,
            "reference": reference,
            "value": float(value),
            "floor": floor,
        }
    return comparisons, ratios, outliers


def corroborate_row(artifact, branch, node, u, *, solvers, rtol=0.01):
    """Corroborate one stored row; ``rtol`` is proposed, not frozen."""
    root, manifest, spec = _artifact_files(artifact)
    row = _row(root, manifest, branch, node, u)
    if not bool(row.get("converged", False)):
        raise CorroborationError("stored LoKI row is not converged")
    node = node or {}
    processes = _processes(spec, manifest, node)
    row["gas_temperature_K"] = _physical_condition(spec, manifest, node, "Tg_K")
    pressure = _physical_condition(spec, manifest, node, "P_Pa")
    stored = _stored_quantities(row, len(processes))
    floors = _floors(row, manifest, len(processes))
    loki = _flat_quantities(
        _loki_quadrature(
            row, processes, spec["solver_options"]["ionizationOperatorType"]
        )
    )
    loki_below = {
        key: (
            key in floors
            and abs(stored[key]) < floors[key]
            and abs(loki[key]) < floors[key]
        )
        for key in stored
    }
    consistency, consistency_ratios, consistency_bad = _compare_quantities(
        stored, loki, floors, loki_below, rtol
    )
    result = {
        "rtol": rtol,
        "rtol_status": "proposed, not frozen",
        "EN_Td": row["EN_Td"],
        "elastic_power_definition": "net elastic transfer including thermal gain",
        "loki_consistency": {
            "verdict": "DISAGREE" if consistency_bad else "AGREE",
            "ratios": consistency_ratios,
            "quantities": consistency,
            "outliers": consistency_bad,
        },
        "solvers": {},
    }
    collision_files = [
        entry["path"]
        for entry in spec["input_files"].values()
        if entry.get("kind") == "cross_section"
    ]
    base = {
        "EN": row["EN_Td"],
        "Tg": row["gas_temperature_K"],
        "P": pressure,
        "processes": processes,
        "emax": float(row["energy_edges_eV"][-1]),
        "ncells": len(row["energy_eV"]),
        "collision_files": collision_files,
        "timeout_s": float(spec.get("timeout_s", 60.0)),
    }
    successful = {}
    solver_specs = (
        solvers.items()
        if isinstance(solvers, dict)
        else ((Path(str(value)).stem, value) for value in solvers)
    )
    for name, value in solver_specs:
        solver = _solver_spec(name, value)
        try:
            reason = _unavailable_reason(
                solver.get("kind", name), spec, node, processes
            )
            if reason:
                raise Unavailable(reason)
            independent = _run_reference(dict(base), solver)
            energy = np.asarray(independent.get("energy"), float)
            expected = np.asarray(row["energy_eV"], float)
            if energy.shape != expected.shape or not np.allclose(
                energy, expected, rtol=0, atol=1e-12
            ):
                raise Unavailable("returned grid does not match stored grid")
            successful[name] = {
                "quantities": _flat_quantities(
                    _quadrature(energy, independent.get("f0"), processes, row)
                ),
                "version": independent.get("version", "unknown"),
            }
        except (Unavailable, CorroborationError, OSError, TypeError, ValueError) as exc:
            result["solvers"][name] = {"verdict": "UNAVAILABLE", "reason": str(exc)}

    all_below = {}
    for key in stored:
        floor = floors.get(key)
        all_below[key] = (
            floor is not None
            and abs(stored[key]) < floor
            and all(
                abs(item["quantities"][key]) < floor for item in successful.values()
            )
        )
    for name, item in successful.items():
        comparisons, ratios, bad = _compare_quantities(
            stored, item["quantities"], floors, all_below, rtol
        )
        result["solvers"][name] = {
            "verdict": "DISAGREE" if bad else "AGREE",
            "ratios": ratios,
            "quantities": comparisons,
            "outliers": bad,
            "version": item["version"],
        }
    verdicts = [item["verdict"] for item in result["solvers"].values()]
    result["overall_verdict"] = (
        "DISAGREE"
        if "DISAGREE" in verdicts
        else "AGREE" if "AGREE" in verdicts else "UNAVAILABLE"
    )
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact")
    parser.add_argument("--branch", default="branch_0")
    parser.add_argument("--u", type=float, required=True)
    parser.add_argument("--node", default="{}", help="JSON composition node")
    parser.add_argument("--rtol", type=float, default=0.01)
    parser.add_argument("--cantera")
    parser.add_argument("--bolos")
    parser.add_argument("--cantera-continuation-from", type=float)
    parser.add_argument("--cantera-continuation-steps", type=int, default=60)
    args = parser.parse_args(argv)
    solvers = {}
    if args.cantera:
        solvers["cantera"] = {"interpreter": args.cantera, "kind": "cantera"}
        if args.cantera_continuation_from is not None:
            solvers["cantera"].update(
                {
                    "continuation_from": args.cantera_continuation_from,
                    "continuation_steps": args.cantera_continuation_steps,
                }
            )
    if args.bolos:
        solvers["bolos"] = {"interpreter": args.bolos, "kind": "bolos"}
    print(
        json.dumps(
            corroborate_row(
                args.artifact,
                args.branch,
                json.loads(args.node),
                args.u,
                solvers=solvers,
                rtol=args.rtol,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
