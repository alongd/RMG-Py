#!/home/alon/anaconda3/envs/rmg_env/bin/python
"""Independent Rouse first-passage reference and tube asymptote report.

Default: simulate and write results.json/results.md. --verify-again reproduces
every deterministic number in results.json into a separate directory and checks
all displayed numbers in pack.md. No RMG package import or database is needed: only the
explicitly named met.py is loaded to map inputs and evaluate the candidate.
"""

import os

# One thread even if the calling environment permits more.
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.util
import json
import math
import multiprocessing
from pathlib import Path
import resource
import sys
import time

import numpy as np
from scipy.fft import rfft, irfft

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
MET_PATH = next(parent / "rmgpy/kmc/met.py" for parent in HERE.parents
                if (parent / "rmgpy/kmc/met.py").is_file())
RUN_DIRECTORY = Path("/home/alon/runs/i030-met-kernel-reference/fixture-verification/reference")


def load_met():
    source_hash = hashlib.sha256(MET_PATH.read_bytes()).hexdigest()
    spec = importlib.util.spec_from_file_location("met_reference_target", MET_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    if hashlib.sha256(MET_PATH.read_bytes()).hexdigest() != source_hash:
        raise AssertionError("target met.py changed during import")
    module._reference_source_sha256 = source_hash
    return module


def mapped_parameters(cfg, met):
    n = cfg["beads_per_chain"]
    temperature = cfg["temperature_K"]
    arm = met.TRANSPORT_ARMS[cfg["rouse_arm"]]
    rg2 = met.PS_C_R2 * met.PS_M0 * n / 6.0
    # Exact finite discrete-chain Rg, NOT the continuum N->infinity identity.
    b2 = rg2 * 6.0 * n / (n * n - 1.0)
    dcm = arm.chain_diffusivity(temperature, n)
    rg = math.sqrt(rg2)
    t0 = rg2 / dcm
    p = np.arange(1, n)
    lam = 4.0 * np.sin(math.pi * p / (2.0 * n)) ** 2
    omega = 3.0 * n / (b2 / rg2) * lam
    return {
        "N": n, "temperature_K": temperature, "Rg_m": rg,
        "bond_rms_m": math.sqrt(b2), "sigma0_m": met.SIGMA_CONTACT,
        "sigma_effective_m": max(met.SIGMA_CONTACT, 2.0 * rg),
        "D0_m2_s": arm.d0(temperature), "D_CM_m2_s": dcm,
        "bead_friction_kg_s": met.K_B * temperature / (n * dcm),
        "spring_constant_N_m": 3.0 * met.K_B * temperature / b2,
        "time_unit_s": t0, "tau_R_reduced": float(1.0 / omega[0]),
        "tau_R_s": float(t0 / omega[0]), "capture_reduced": met.SIGMA_CONTACT / rg,
        "rate_unit_m3_s": dcm * rg, "rate_unit_m3_mol_s": met.N_A * dcm * rg,
        "bond_squared_reduced": b2 / rg2,
        "entanglement_threshold_units": met.TRANSPORT_ARMS[cfg["entangled_arm"]].n_e,
        "mid_bead_index_zero_based": n // 2 - 1,
    }


def mode_coefficients(n, bond_squared, pair_class):
    p = np.arange(1, n)
    lam = 4.0 * np.sin(math.pi * p / (2.0 * n)) ** 2
    omega = 3.0 * n / bond_squared * lam
    end = np.sqrt(2.0 / n) * np.cos(math.pi * p * 0.5 / n)
    mid = np.sqrt(2.0 / n) * np.cos(math.pi * p * (n // 2 - 0.5) / n)
    sites = {"end/end": (end, end), "end/mid": (end, mid), "mid/mid": (mid, mid)}
    left, right = sites[pair_class]
    variances = (left * left + right * right) * bond_squared / (3.0 * lam)
    # Relative COM D=2; the other modes supply 2(N-1), hence local D=2N.
    assert abs(float(np.sum(variances * omega)) + 2.0 - 2.0 * n) < 1e-10
    return omega, variances


def gaussian_spectrum(steps, dt, omega, variances):
    length = 1 << (2 * steps).bit_length()
    lag = np.arange(length // 2 + 1, dtype=float) * dt
    cov = np.zeros_like(lag)
    for rate, variance in zip(omega, variances):
        cov += variance * np.exp(-rate * lag)
    first_column = np.concatenate((cov, cov[-2:0:-1]))
    spectrum = rfft(first_column, workers=1).real
    if np.min(spectrum) < -1e-10 * cov[0]:
        raise AssertionError("circulant embedding is not positive")
    return np.sqrt(np.maximum(spectrum, 0.0)), length


def survival_and_rates(hit, volume, edges):
    count = len(hit)
    alive = [int(np.sum(hit > t)) for t in edges]
    survival = np.asarray(alive, dtype=float) / count
    se = np.sqrt(survival * (1.0 - survival) / count)
    bins = []
    for i, (t1, t2) in enumerate(zip(edges[:-1], edges[1:])):
        start, end = alive[i:i+2]
        if end == 0:
            raise AssertionError("insufficient survivors")
        k = volume / (t2 - t1) * math.log(start / end)
        kse = volume / (t2 - t1) * math.sqrt((start - end) / (start * end))
        bins.append({"t1": t1, "t2": t2, "survivors_start": start,
                     "survivors_end": end, "k_reduced": k, "se_reduced": kse})
    return {"survival": survival.tolist(), "survival_se": se.tolist(), "bins": bins}


def corrected_long(hit, run, cfg, lower=None, upper=None):
    lo, hi = cfg["long_window_reduced"]
    lo = lo if lower is None else lower
    hi = hi if upper is None else upper
    start, end = int(np.sum(hit > lo)), int(np.sum(hit > hi))
    if end == 0 or start == end:
        raise AssertionError("long-time reference needs nonzero events and survivors")
    volume = run["box_over_Rg"] ** 3
    kbox = volume / (hi - lo) * math.log(start / end)
    error = volume / (hi - lo) * math.sqrt((start - end) / (start * end))
    # Leading periodic Green-function correction, relative COM diffusivity=2.
    correction = cfg["finite_box_green_constant"] / (4.0 * math.pi * 2.0 * run["box_over_Rg"])
    k = kbox / (1.0 + correction * kbox)
    se = error / (1.0 + correction * kbox) ** 2
    return {"window": [lo, hi], "events": start - end, "k_box_reduced": kbox,
            "se_box_reduced": error, "k_infinite_reduced": k, "se_infinite_reduced": se,
            "effective_radius_over_Rg": k / (8.0 * math.pi)}


def simulate(run, cfg, mapped, pair_class, class_index):
    dt = run["dt_reduced"]
    tmax = run.get("tmax_reduced", cfg["time_edges_reduced"][-1])
    steps = int(round(tmax / dt))
    n = mapped["N"]
    omega, variances = mode_coefficients(n, mapped["bond_squared_reduced"], pair_class)
    root_spectrum, length = gaussian_spectrum(steps, dt, omega, variances)
    rng = np.random.Generator(np.random.PCG64(run["seed"] + class_index))
    # FFT length determines batch size; all paths stay under 256 MiB of arrays.
    batch_size = max(1, min(32, int(2**20 // length)))
    a = mapped["capture_reduced"]
    shifted_a = a + cfg["boundary_shift_constant"] * math.sqrt(4.0 * n * dt)
    side = run["box_over_Rg"]
    hit = np.full(run["replicas"], np.inf)
    for offset in range(0, len(hit), batch_size):
        batch = min(batch_size, len(hit) - offset)
        initial = rng.uniform(-side / 2.0, side / 2.0, (batch, 3))
        inside = np.sum(initial * initial, axis=1) <= a * a
        while np.any(inside):
            initial[inside] = rng.uniform(-side / 2.0, side / 2.0, (int(inside.sum()), 3))
            inside = np.sum(initial * initial, axis=1) <= a * a
        squared_distance = np.zeros((batch, steps + 1))
        for axis in range(3):
            white = rng.standard_normal((batch, length))
            projected = irfft(rfft(white, axis=1, workers=1) * root_spectrum,
                              n=length, axis=1, workers=1)[:, :steps+1].copy()
            del white
            projected -= projected[:, :1].copy()
            increments = rng.normal(0.0, math.sqrt(4.0 * dt), (batch, steps))
            projected[:, 1:] += np.cumsum(increments, axis=1)
            del increments
            projected += initial[:, axis, None]
            projected -= side * np.floor(projected / side + 0.5)
            squared_distance += projected * projected
            del projected
        contact = squared_distance <= shifted_a * shifted_a
        # Initial pairs were conditioned outside the physical microscopic sink.
        contact[:, 0] = False
        first = np.argmax(contact, axis=1)
        occurred = np.any(contact, axis=1)
        hit[offset:offset+batch][occurred] = first[occurred] * dt
    edges = [edge for edge in cfg["time_edges_reduced"] if edge <= tmax]
    summary = survival_and_rates(hit, side**3, edges)
    summary.update({"name": run["name"], "pair_class": pair_class,
                    "seed": run["seed"] + class_index, "replicas": run["replicas"],
                    "dt_reduced": dt, "box_over_Rg": side, "numerical_capture_reduced": shifted_a,
                    "concentration_mol_m3": 1.0 / (mapped["Rg_m"]**3 * side**3 * 6.02214076e23),
                    "hit_indices_sha256": hashlib.sha256(
                        np.where(np.isfinite(hit), np.rint(hit / dt), -1).astype("<i8").tobytes()).hexdigest()})
    if tmax >= cfg["long_window_reduced"][1]:
        summary["long"] = corrected_long(hit, run, cfg)
        summary["late_halves"] = [corrected_long(hit, run, cfg, 4.0, 6.0),
                                  corrected_long(hit, run, cfg, 6.0, 8.0)]
    return summary


def gaussian_checks(mapped):
    """Check exact covariance sampling against an independently built bond Laplacian."""
    n = mapped["N"]
    lap = np.diag([1.0] + [2.0] * (n-2) + [1.0])
    lap += np.diag([-1.0] * (n-1), 1) + np.diag([-1.0] * (n-1), -1)
    eig, vec = np.linalg.eigh(lap)
    rg2 = float(np.sum(mapped["bond_squared_reduced"] / eig[1:]) / n)
    assert abs(rg2 - 1.0) < 1e-12
    errors = []
    for pair_class in ("end/end", "end/mid", "mid/mid"):
        omega, var = mode_coefficients(n, mapped["bond_squared_reduced"], pair_class)
        indices = {"end/end": (0, 0), "end/mid": (0, n//2-1), "mid/mid": (n//2-1, n//2-1)}
        left, right = indices[pair_class]
        direct = ((vec[left, 1:]**2 + vec[right, 1:]**2)
                  * mapped["bond_squared_reduced"] / (3.0 * eig[1:]))
        errors.append(float(np.max(np.abs(direct - var))))
        assert errors[-1] < 1e-12
        root, length = gaussian_spectrum(1000, 1e-4, omega, var)
        recovered = irfft(root * root, n=length)
        target = np.array([np.sum(var * np.exp(-omega * lag * 1e-4)) for lag in range(10)])
        assert np.max(np.abs(recovered[:10] - target)) < 1e-12
    return {"Rg_squared_reduced": rg2, "max_mode_variance_error": max(errors),
            "FFT_covariance_check": "PASS", "local_relative_D_reduced": 2.0*n}


def compare(rows, cfg, kernel):
    output = []
    z = cfg["statistical_sigmas"]
    for pair_class in cfg["pair_classes"]:
        own = {row["name"]: row for row in rows if row["pair_class"] == pair_class}
        baseline = own["base"]["long"]
        checks = {}
        for label, allowance in (("step", cfg["step_relative_allowance"]),
                                 ("box", cfg["box_relative_allowance"])):
            alternative = own[label]["long"]
            delta = abs(alternative["k_infinite_reduced"] - baseline["k_infinite_reduced"])
            band = z * math.hypot(alternative["se_infinite_reduced"], baseline["se_infinite_reduced"])
            band += allowance * baseline["k_infinite_reduced"]
            checks[label] = {"absolute_difference": delta, "allowed_difference": band, "pass": delta <= band}
        first, second = own["base"]["late_halves"]
        delta = abs(first["k_infinite_reduced"] - second["k_infinite_reduced"])
        band = z * math.hypot(first["se_infinite_reduced"], second["se_infinite_reduced"])
        band += cfg["plateau_relative_allowance"] * baseline["k_infinite_reduced"]
        checks["plateau"] = {"absolute_difference": delta, "allowed_difference": band, "pass": delta <= band}
        reference = own["step"]["long"]
        kref = reference["k_infinite_reduced"]
        se = reference["se_infinite_reduced"]
        log_band = math.log(cfg["long_model_factor"]) + z * se / kref
        discrepancy = abs(math.log(kernel[pair_class] / kref))
        bins = []
        for row in own["transient"]["bins"]:
            if row["t2"] <= 0.2:
                predicted = kernel[pair_class]
                band = cfg["transient_relative_allowance"] * row["k_reduced"] + z * row["se_reduced"]
                bins.append({"t1": row["t1"], "t2": row["t2"],
                             "kernel_over_reference": predicted / row["k_reduced"] if row["k_reduced"] else None,
                             "diagnostic_agreement": abs(predicted - row["k_reduced"]) <= band})
        output.append({"pair_class": pair_class, "numerical_checks": checks,
                       "reference_reduced": kref, "reference_se_reduced": se,
                       "kernel_reduced": kernel[pair_class], "kernel_over_reference": kernel[pair_class] / kref,
                       "absolute_log_discrepancy": discrepancy, "allowed_log_discrepancy": log_band,
                       "long_scale_pass": discrepancy <= log_band,
                       "reference_valid": all(check["pass"] for check in checks.values()),
                       "transient_reporting_only": bins})
    return output


def tube_asymptotes(cfg, met):
    """Literature exponents; convention-dependent amplitudes are NOT simulated."""
    n = cfg["entangled_units"]
    arm = met.TRANSPORT_ARMS[cfg["entangled_arm"]]
    ne = arm.n_e
    b2 = met.PS_C_R2 * met.PS_M0
    d0 = arm.d0(cfg["temperature_K"])
    ta = b2 / d0  # scaling time convention, not a measured local collision time
    te, tr, td = ne**2 * ta, n**2 * ta, n**3 / ne * ta
    rg = math.sqrt(b2 * n / 6.0)
    dcm = arm.chain_diffusivity(cfg["temperature_K"], n)
    assert ta < te < tr < td
    regimes = [
        {"name": "pre-tube Rouse", "window_s": [ta, te], "rms_exponent": 0.25, "k_exponent": -0.25},
        {"name": "tube breathing", "window_s": [te, tr], "rms_exponent": 0.125, "k_exponent": -0.625},
        {"name": "coherent curvilinear reptation", "window_s": [tr, td], "rms_exponent": 0.25, "k_exponent": -0.25},
        {"name": "relaxed COM", "window_s": [td, math.inf], "rms_exponent": 0.5, "k_exponent": 0.0},
    ]
    regimes[-1]["window_s"][1] = None
    rows = []
    for pair_class in cfg["pair_classes"]:
        rows.append({"pair_class": pair_class, "transient_regimes": regimes,
                     "k_long_m3_s": "C_pair * D_CM * Rg; C_pair not identified",
                     "C_pair_from_candidate": met.diffusion_rate(arm, cfg["temperature_K"], n, n,
                           pair_class, spin_factor=1.0) / (met.N_A * dcm * rg),
                     "numeric_reference_error": None, "comparison_policy": "report-only asymptotic scaling"})
    return {"units": n, "Ne_units": ne, "D_CM_m2_s": dcm, "Rg_m": rg,
            "bond_rms_m": math.sqrt(b2), "ta_s": ta, "te_s": te, "tau_R_s": tr,
            "tau_d_s": td, "rows": rows,
            "assumption": "extend representative-monomer tube exponents to all three pairs; no site prefactor claim"}


def markdown(result):
    mapped = result["mapped"]
    unit = mapped["rate_unit_m3_mol_s"]
    lines = ["<!-- BEGIN REPRODUCED NUMBERS -->", "",
             "All ± entries below are one Monte Carlo standard error; systematic/model errors are separate.", "",
             "| Quantity | Value |", "|---|---:|"]
    for key, value in mapped.items():
        lines.append(f"| `{key}` | {value:.12g} |")
    lines += ["", "| Run | Pair | seed | microscopic numerical radius / Rg | L / Rg | replicas | dt / t0 | late events | raw k_box / (D_CM Rg) | corrected k∞ / (D_CM Rg) |",
              "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for row in result["simulation"]:
        if "long" not in row:
            continue
        long = row["long"]
        lines.append(f"| {row['name']} | {row['pair_class']} | {row['seed']} | {row['numerical_capture_reduced']:.8g} | {row['box_over_Rg']:g} | {row['replicas']} | {row['dt_reduced']:g} | {long['events']} | {long['k_box_reduced']:.8g} ± {long['se_box_reduced']:.8g} | {long['k_infinite_reduced']:.8g} ± {long['se_infinite_reduced']:.8g} |")
    lines += ["", "Final long-time oracle uses the step-refined run, as declared before simulation.", "",
              "| Pair | Reference k∞ (m³ mol⁻¹ s⁻¹) | Candidate k (m³ mol⁻¹ s⁻¹) | Candidate / reference | numerical audits | factor-two scale comparison |",
              "|---|---:|---:|---:|---|---|"]
    for row in result["comparison"]:
        lines.append(f"| {row['pair_class']} | {row['reference_reduced']*unit:.8g} ± {row['reference_se_reduced']*unit:.8g} | {row['kernel_reduced']*unit:.8g} | {row['kernel_over_reference']:.8g} | {'PASS' if row['reference_valid'] else 'FAIL'} | {'PASS' if row['long_scale_pass'] else 'FAIL'} |")
    for row in result["simulation"]:
        if row["name"] not in ("step", "transient"):
            continue
        label = "Step-refined" if row["name"] == "step" else "Dedicated transient"
        lines += ["", f"{label} {row['pair_class']}; L/Rg={row['box_over_Rg']:g}, replicas={row['replicas']}, "
                  f"seed={row['seed']}, dt/t0={row['dt_reduced']:g}. Time is in t0; bin k is a finite-box hazard coefficient.", "",
                  "| t1 | t2 | S(t2) ± SE | k_bin / (D_CM Rg) ± SE | k_bin (m³ mol⁻¹ s⁻¹) ± SE |",
                  "|---:|---:|---:|---:|---:|"]
        for index, item in enumerate(row["bins"]):
            lines.append(f"| {item['t1']:g} | {item['t2']:g} | {row['survival'][index+1]:.8g} ± {row['survival_se'][index+1]:.8g} | {item['k_reduced']:.8g} ± {item['se_reduced']:.8g} | {item['k_reduced']*unit:.8g} ± {item['se_reduced']*unit:.8g} |")
    lines += ["", "| Pair | Transient t1 | t2 | Candidate / finite-box reference | 20% + 4 SE diagnostic |",
              "|---|---:|---:|---:|---|"]
    for row in result["comparison"]:
        for item in row["transient_reporting_only"]:
            ratio = item["kernel_over_reference"]
            ratio_text = f"{ratio:.8g}" if ratio is not None else "undefined (no events)"
            lines.append(f"| {row['pair_class']} | {item['t1']:g} | {item['t2']:g} | {ratio_text} | {'within band' if item['diagnostic_agreement'] else 'outside band'} (report-only) |")
    lines += ["", "| Entangled theory parameter | Value |", "|---|---:|"]
    for key in ("units", "Ne_units", "D_CM_m2_s", "Rg_m", "bond_rms_m", "ta_s", "te_s", "tau_R_s", "tau_d_s"):
        lines.append(f"| `{key}` | {result['tube'][key]:.12g} |")
    lines += ["", "Entangled entries are scaling predictions with unidentified prefactors, not measured rates or error bars.",
              "", "<!-- END REPRODUCED NUMBERS -->"]
    return "\n".join(lines) + "\n"


def run_one(run, cfg, mapped, pair_class, index):
    begin = time.monotonic()
    print(f"RUN {run['name']} {pair_class} seed={run['seed']+index}", flush=True)
    row = simulate(run, cfg, mapped, pair_class, index)
    if "long" in row:
        value = f"k∞={row['long']['k_infinite_reduced']:.8g} ± {row['long']['se_infinite_reduced']:.8g} (D_CM Rg)"
    else:
        value = f"S(0.2t0)={row['survival'][-1]:.8g}; short-time bins measured"
    print(f"MEASURED {run['name']} {pair_class}: {value}; elapsed={time.monotonic()-begin:.1f}s", flush=True)
    return row


def execute_jobs(jobs, workers):
    if workers == 1:
        rows = [run_one(*job) for job in jobs]
    else:
        with ProcessPoolExecutor(max_workers=workers,
                                 mp_context=multiprocessing.get_context("spawn")) as pool:
            futures = [pool.submit(run_one, *job) for job in jobs]
            # Stable order and separate per-job seeds preserve the serial results exactly.
            rows = [future.result() for future in futures]
    return rows


def run_all(cfg, workers=3):
    met = load_met()
    mapped = mapped_parameters(cfg, met)
    checks = gaussian_checks(mapped)
    jobs = [(run, cfg, mapped, pair_class, index) for run in cfg["runs"]
            for index, pair_class in enumerate(cfg["pair_classes"])]
    rows = execute_jobs(jobs, workers)
    kernel = {pair: met.diffusion_rate(met.TRANSPORT_ARMS[cfg["rouse_arm"]], cfg["temperature_K"],
                                    mapped["N"], mapped["N"], pair, spin_factor=cfg["spin_factor"])
              / mapped["rate_unit_m3_mol_s"] for pair in cfg["pair_classes"]}
    result = {"schema": 1, "parameters_sha256": hashlib.sha256((HERE/"parameters.json").read_bytes()).hexdigest(),
              "met_source_sha256": met._reference_source_sha256,
              "mapped": mapped, "independent_checks": checks, "simulation": rows,
              "comparison": compare(rows, cfg, kernel), "tube": tube_asymptotes(cfg, met)}
    if hashlib.sha256(MET_PATH.read_bytes()).hexdigest() != met._reference_source_sha256:
        raise AssertionError("target met.py changed during the run; results cannot be labelled with one source hash")
    for row in result["comparison"]:
        print(f"CANDIDATE {row['pair_class']}: numerical_audits={'PASS' if row['reference_valid'] else 'FAIL'}, "
              f"long_scale={'PASS' if row['long_scale_pass'] else 'FAIL'}, "
              f"ratio={row['kernel_over_reference']:.8g}; transient=REPORT-ONLY", flush=True)
    return result


def supplement_transient(cfg, saved, workers):
    """Retain completed long-time measurements; add the predeclared precision extension.

    The complete --verify-again run does not use this shortcut: it repeats all
    twelve ensembles and checks every value and contact-index hash together.
    """
    initial = json.loads((HERE / "parameters_initial.json").read_text())
    for key, value in initial.items():
        if key != "runs":
            assert cfg[key] == value, f"original policy/model changed: {key}"
    assert cfg["runs"][:3] == initial["runs"]
    assert {row["name"] for row in saved["simulation"]} == {"base", "step", "box"}
    met = load_met()
    mapped = mapped_parameters(cfg, met)
    assert saved["mapped"] == mapped
    assert saved["met_source_sha256"] == hashlib.sha256(MET_PATH.read_bytes()).hexdigest()
    run = next(row for row in cfg["runs"] if row["name"] == "transient")
    jobs = [(run, cfg, mapped, pair, index) for index, pair in enumerate(cfg["pair_classes"])]
    saved["simulation"].extend(execute_jobs(jobs, workers))
    kernel = {row["pair_class"]: row["kernel_reduced"] for row in saved["comparison"]}
    saved["comparison"] = compare(saved["simulation"], cfg, kernel)
    saved["parameters_sha256"] = hashlib.sha256((HERE / "parameters.json").read_bytes()).hexdigest()
    print("SUPPLEMENT COMPLETE: added three dedicated transient ensembles; original model and tolerances unchanged.")
    return saved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=RUN_DIRECTORY)
    parser.add_argument("--verify-again", action="store_true")
    parser.add_argument("--supplement-transient", action="store_true",
                        help="append the declared short-time precision extension to an initial nine-run result")
    parser.add_argument("--workers", type=int, choices=(1, 2, 3), default=3,
                        help="independent one-thread ensembles; hard cap is three cores")
    args = parser.parse_args()
    if args.verify_again and args.supplement_transient:
        parser.error("verification always reruns every ensemble")
    cfg = json.loads((HERE / "parameters.json").read_text())
    previous = json.loads((HERE / "results/results.json").read_text()) if args.verify_again else None
    if args.supplement_transient:
        saved = json.loads((args.output / "results.json").read_text())
        result = supplement_transient(cfg, saved, args.workers)
    else:
        result = run_all(cfg, args.workers)
    destination = args.output / "reproduced" if args.verify_again else args.output
    destination.mkdir(parents=True, exist_ok=True)
    rendered = markdown(result)
    (destination / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    (destination / "results.md").write_text(rendered)
    if args.verify_again:
        if previous != result:
            raise AssertionError("full rerun differs from the saved deterministic reference")
        pack = (HERE.parent / "pack.md").read_text()
        if rendered.rstrip() not in pack:
            raise AssertionError("pack numbers do not match the freshly reproduced table")
        print(f"VERIFIER PASS: all {len(result['simulation'])} first-passage runs reproduced exactly; every generated pack number matches.")
    else:
        print(f"RESULTS written to {destination}")
    self_peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    child_peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024.0
    print(f"PARENT_PEAK_RSS_MiB={self_peak:.1f}; CHILD_MAX_RSS_MiB={child_peak:.1f}; "
          f"workers={args.workers}; FFT/BLAS_threads_per_worker=1", flush=True)


if __name__ == "__main__":
    main()
