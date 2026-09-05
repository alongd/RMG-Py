#!/usr/bin/env python3
"""I-099 two-temperature polymer mass-loss probe.

Standalone measurement script (NOT a pytest test -- it defines no `test_`
functions and no `Test*` classes, so pytest collection skips it).

Loads each of the four i099 demonstration decks through RMG's own input
reader, initializes the HybridPolymerReactor/HybridPolymerSystem against the
model core exactly as `rmgpy/rmg/main.py` does before simulating, then drives
the solver by hand over an explicit log-spaced time grid (the pattern proven
in `rmgpy/tools/polymer_moments_runner.py`'s module docstring and its
`run_segments()` loop: for a single, already-isothermal segment,
`initialize_model()` performs the full `set_initial_conditions()` /
`generate_rate_coefficients()` / `set_initial_derivative()` /
`initialize_solver()` sequence internally, so no manual restart dance is
needed -- just call `solver.advance(t)` once per grid point).

For each grid point it records the time, the condensed-phase polymer mass,
and the molar amount of every core species (which for these decks is exactly
the gas species plus each polymer pool's mu0/mu1/mu2 moment coordinates).
Condensed mass is computed with a local transcription of the NORMATIVE
formula documented at `rmgpy/tools/polymer_moments_runner.py:649-707` and
implemented as `PolymerPoolConfig.condensed_mass_g()` at
`rmgpy/solver/polymer.pyx:1171-1180`:

    condensed_mass_g = mu1_mol * monomer_mw_g_mol - mu0_mol * chain_mass_defect_g_mol

summed over every pool in `solver.polymer_pools`. `chain_mass_defect_g_mol`
defaults to 0.0 on an ordinary pool (neither i099 deck declares a feature
pool), so this reduces to `mu1 * monomer_mw_g_mol` here -- but the full
two-term form is transcribed so the script stays correct if that ever
changes.

Per-pool `monomer_mw_g_mol` / `chain_mass_defect_g_mol` are read directly off
the LIVE `PolymerPoolConfig` objects the solver itself built
(`solver.polymer_pools`), not off the `<run_dir>/chemkin/polymer_pools.json`
sidecar the brief mentions as the normal source -- that sidecar is a
serialization of the exact same values and does not exist for
`arrhenius_T700` (confirmed: only `scalar_T700`, `scalar_T900` and
`arrhenius_T900` have one on disk), and reading the live pool config is
guaranteed consistent with what the residual itself used, so this script
does that uniformly for all four decks and says so once here rather than
falling back deck-by-deck.

Each deck runs in its own subprocess. The RMG thermo database is a
process-global singleton (see MEMORY: I-189, `project_i189_db_singleton_isolation.md`);
re-using one process across decks with different pool parameterisations is
exactly the kind of cross-deck leakage that singleton has bitten before, and
correctness comes first per the brief, so this script always pays the
database-load cost once per deck rather than risk it.

Usage:
    python test/rmgpy/solver/i099_temperature_response.py [--decks DIR [DIR ...]]
    python test/rmgpy/solver/i099_temperature_response.py --help

Run from `/home/alon/Code/RMG-Py-i099-arrhenius-empty` with `rmg_env` active:
    source /home/alon/anaconda3/etc/profile.d/conda.sh && conda activate rmg_env
"""

import argparse
import contextlib
import io
import json
import os
import subprocess
import sys
import time

REPO_ROOT = "/home/alon/Code/RMG-Py-i099-arrhenius-empty"

DEFAULT_DECKS = [
    "/home/alon/runs/RMG/i099_arrhenius/scalar_T700",
    "/home/alon/runs/RMG/i099_arrhenius/scalar_T900",
    "/home/alon/runs/RMG/i099_arrhenius/arrhenius_T700",
    "/home/alon/runs/RMG/i099_arrhenius/arrhenius_T900",
]

RESULT_MARKER = "I099_RESULT_JSON:"

DEFAULT_BUDGET_S = 12 * 60.0
DEFAULT_N_POINTS = 40
DEFAULT_ATOL = 1e-16
DEFAULT_RTOL = 1e-8


def _check_environment():
    """Gate section: verify we are running the tree we think we are, and
    print the compiled solver backend's MRO so every run shows
    `pydas.daspk.DASPK` is really what is being driven."""
    import rmgpy

    real = os.path.realpath(rmgpy.__file__)
    assert real.startswith(REPO_ROOT), (
        f"rmgpy resolved to {real!r}, not under {REPO_ROOT!r} -- refusing to run "
        "against the wrong checkout."
    )
    import rmgpy.solver.base as b

    mro = [c.__module__ + "." + c.__name__ for c in b.ReactionSystem.__mro__]
    print("ReactionSystem MRO:", mro, flush=True)
    assert any(m == "pydas.daspk.DASPK" for m in mro), (
        f"pydas.daspk.DASPK not found in ReactionSystem.__mro__ ({mro}); "
        "this is not the compiled DASPK backend the brief requires."
    )
    return mro


def time_grid(n_points=DEFAULT_N_POINTS):
    """t = 0 plus n_points log-spaced points from 1e-6 s to 1.0 s."""
    import numpy as np

    return np.concatenate(([0.0], np.logspace(-6.0, 0.0, n_points)))


def condensed_mass_g(mu0_mol, mu1_mol, monomer_mw_g_mol, chain_mass_defect_g_mol):
    """Local transcription of the NORMATIVE formula, verbatim from
    rmgpy/tools/polymer_moments_runner.py:649-707 and
    rmgpy/solver/polymer.pyx:1171-1180:

        condensed_mass_g = mu1*monomer_mw_g_mol - mu0*chain_mass_defect_g_mol

    `chain_mass_defect_g_mol` is 0.0 on an ordinary (non-feature) pool, which
    is every pool in the four i099 decks, so this term is a no-op here -- but
    the full two-term form is kept so the transcription matches the source
    exactly rather than a simplified special case.
    """
    return mu1_mol * monomer_mw_g_mol - mu0_mol * chain_mass_defect_g_mol


def _total_condensed_mass_g(pools, y):
    total = 0.0
    for pool in pools:
        mu0_idx, mu1_idx, _mu2_idx = pool.mu_indices
        total += condensed_mass_g(
            float(y[mu0_idx]), float(y[mu1_idx]),
            float(pool.monomer_mw_g_mol), float(pool.chain_mass_defect_g_mol),
        )
    return total


def _write_csv(path, rows, species_labels):
    import csv

    fieldnames = ["t_s", "condensed_mass_g"] + list(species_labels)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def run_one_deck(deck_dir, budget_s=DEFAULT_BUDGET_S, n_points=DEFAULT_N_POINTS,
                  atol=DEFAULT_ATOL, rtol=DEFAULT_RTOL):
    """Worker body: runs entirely inside its own subprocess. Returns a JSON-
    serializable result dict; never raises for an ordinary integration
    failure (that is reported as data), only for setup problems that make
    the deck impossible to even attempt."""
    import tempfile

    mro = _check_environment()

    from rmgpy.rmg.main import RMG

    deck_dir = os.path.abspath(deck_dir)
    input_file = os.path.join(deck_dir, "input.py")
    if not os.path.isfile(input_file):
        raise FileNotFoundError(input_file)

    scratch = tempfile.mkdtemp(prefix="i099_scratch_")
    print(f"[{deck_dir}] loading deck via RMG(input_file=...).initialize() "
          f"(scratch={scratch}) ...", flush=True)

    t_setup_start = time.monotonic()
    rmg = RMG(input_file=input_file, output_directory=scratch)
    rmg.initialize()
    print(f"[{deck_dir}] initialize() done in {time.monotonic() - t_setup_start:.1f}s; "
          f"core species={len(rmg.reaction_model.core.species)}, "
          f"core reactions={len(rmg.reaction_model.core.reactions)}", flush=True)

    rs = rmg.reaction_systems[0]
    core = rmg.reaction_model.core
    edge = rmg.reaction_model.edge

    rs.initialize_model(
        core_species=core.species,
        core_reactions=core.reactions,
        edge_species=edge.species,
        edge_reactions=edge.reactions,
        atol=atol,
        rtol=rtol,
    )
    solver = rs.solver  # wrapper does NOT proxy ODE state -- go through .solver
    species_labels = [s.label for s in core.species]
    pools = list(solver.polymer_pools)
    print(f"[{deck_dir}] pools: " + ", ".join(
        f"{p.label} (monomer_mw_g_mol={p.monomer_mw_g_mol:g}, "
        f"chain_mass_defect_g_mol={p.chain_mass_defect_g_mol:g}, "
        f"mu_indices={p.mu_indices})" for p in pools), flush=True)

    grid = time_grid(n_points)
    rows = []
    furthest_t_s = 0.0
    error = None
    t_integrate_start = time.monotonic()
    for t in grid:
        t = float(t)
        if t == 0.0:
            y = solver.y0
        else:
            elapsed = time.monotonic() - t_integrate_start
            if elapsed > budget_s:
                error = (f"wall-clock budget of {budget_s:g}s exceeded "
                         f"(elapsed={elapsed:.1f}s) before reaching t={t:g}s; "
                         f"furthest time actually reached = {furthest_t_s:g}s")
                print(f"[{deck_dir}] {error}", flush=True)
                break
            try:
                solver.advance(t)
            except Exception as exc:  # noqa: BLE001 -- report, do not swallow
                error = f"{type(exc).__name__} at t={t:g}s: {exc}"
                print(f"[{deck_dir}] INTEGRATION FAILED: {error}", flush=True)
                break
            y = solver.y
        mass = _total_condensed_mass_g(pools, y)
        row = {"t_s": t, "condensed_mass_g": mass}
        for i, lab in enumerate(species_labels):
            row[lab] = float(y[i])
        rows.append(row)
        furthest_t_s = t

    wall_s = time.monotonic() - t_integrate_start
    csv_path = os.path.join(deck_dir, "i099_trajectory.csv")
    _write_csv(csv_path, rows, species_labels)
    print(f"[{deck_dir}] wrote {csv_path} ({len(rows)} rows, "
          f"furthest_t_s={furthest_t_s:g}, wall_s={wall_s:.1f})", flush=True)

    return {
        "deck": deck_dir,
        "csv": csv_path,
        "mro_ok": True,
        "n_points_requested": len(grid),
        "n_rows_written": len(rows),
        "furthest_t_s": furthest_t_s,
        "error": error,
        "wall_s": wall_s,
        "grid": [float(t) for t in grid],
        "condensed_mass_g": [r["condensed_mass_g"] for r in rows],
        "t_s": [r["t_s"] for r in rows],
    }


def _run_worker_subprocess(deck_dir, budget_s, n_points, atol, rtol):
    cmd = [
        sys.executable, os.path.abspath(__file__),
        "--worker", deck_dir,
        "--budget-s", str(budget_s),
        "--n-points", str(n_points),
        "--atol", str(atol),
        "--rtol", str(rtol),
    ]
    print(f"=== {deck_dir}: launching worker subprocess ===", flush=True)
    try:
        proc = subprocess.run(
            cmd, cwd=REPO_ROOT, capture_output=True, text=True,
            timeout=budget_s + 900.0,  # generous allowance for DB load, on top
                                        # of the in-script 12-minute integration
                                        # budget the worker enforces itself
        )
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout or ""
        stderr = exc.stderr or ""
        sys.stdout.write(stdout)
        sys.stderr.write(stderr)
        return {
            "deck": deck_dir, "error": f"subprocess hard-timeout: {exc}",
            "furthest_t_s": None, "condensed_mass_g": [], "t_s": [], "grid": [],
        }

    sys.stdout.write(proc.stdout)
    sys.stderr.write(proc.stderr)

    for line in reversed(proc.stdout.splitlines()):
        if line.startswith(RESULT_MARKER):
            return json.loads(line[len(RESULT_MARKER):])

    return {
        "deck": deck_dir,
        "error": (f"worker subprocess exited {proc.returncode} without printing a "
                  f"result marker; stderr tail: {proc.stderr[-2000:]}"),
        "furthest_t_s": None, "condensed_mass_g": [], "t_s": [], "grid": [],
    }


def _print_summary(results):
    print("\n=== i099 summary: condensed mass [g] vs. time, by parameterisation ===")
    by_arm = {}
    for r in results:
        base = os.path.basename(r["deck"].rstrip("/"))
        if "_T" not in base:
            continue
        arm, t_tag = base.rsplit("_T", 1)
        by_arm.setdefault(arm, {})[t_tag] = r

    any_failure = any(r.get("error") for r in results)

    for arm, by_temp in sorted(by_arm.items()):
        temps = sorted(by_temp.keys(), key=lambda s: float(s))
        if len(temps) < 2:
            print(f"\n-- arm '{arm}': only {len(temps)} deck(s) present, no comparison")
            continue
        t_lo, t_hi = temps[0], temps[1]
        r_lo, r_hi = by_temp[t_lo], by_temp[t_hi]
        print(f"\n-- arm '{arm}': T={t_lo}K vs T={t_hi}K")
        if r_lo.get("error"):
            print(f"   T={t_lo}K FAILED to integrate: {r_lo['error']} "
                  f"(furthest_t_s={r_lo.get('furthest_t_s')})")
        if r_hi.get("error"):
            print(f"   T={t_hi}K FAILED to integrate: {r_hi['error']} "
                  f"(furthest_t_s={r_hi.get('furthest_t_s')})")

        grid_lo = r_lo.get("t_s") or []
        grid_hi = r_hi.get("t_s") or []
        common_t = sorted(set(round(t, 12) for t in grid_lo) &
                           set(round(t, 12) for t in grid_hi))
        if not common_t:
            print("   no common time points recorded on both sides -- nothing to compare")
            continue

        idx_lo = {round(t, 12): i for i, t in enumerate(grid_lo)}
        idx_hi = {round(t, 12): i for i, t in enumerate(grid_hi)}

        # Print a representative subset: every ~5th common point.
        stride = max(1, len(common_t) // 8)
        sample_t = common_t[::stride]
        if common_t[-1] not in sample_t:
            sample_t.append(common_t[-1])

        print(f"   {'t_s':>14}  {t_lo+'K':>16}  {t_hi+'K':>16}  {'rel_diff':>12}")
        max_abs_diff = 0.0
        max_rel_diff = 0.0
        for t in sample_t:
            m_lo = r_lo["condensed_mass_g"][idx_lo[t]]
            m_hi = r_hi["condensed_mass_g"][idx_hi[t]]
            abs_diff = abs(m_hi - m_lo)
            denom = max(abs(m_lo), abs(m_hi), 1e-300)
            rel_diff = abs_diff / denom
            print(f"   {t:14.6e}  {m_lo:16.8e}  {m_hi:16.8e}  {rel_diff:12.6e}")

        for t in common_t:
            m_lo = r_lo["condensed_mass_g"][idx_lo[t]]
            m_hi = r_hi["condensed_mass_g"][idx_hi[t]]
            abs_diff = abs(m_hi - m_lo)
            denom = max(abs(m_lo), abs(m_hi), 1e-300)
            max_abs_diff = max(max_abs_diff, abs_diff)
            max_rel_diff = max(max_rel_diff, abs_diff / denom)

        if max_abs_diff == 0.0:
            print(f"   -> bit-identical condensed-mass columns across all "
                  f"{len(common_t)} common points (max abs diff = 0.0)")
        else:
            print(f"   -> max abs diff over {len(common_t)} common points = "
                  f"{max_abs_diff:.6e} g; max rel diff = {max_rel_diff:.6e}")

    print()
    for r in results:
        if r.get("error"):
            print(f"DECK FAILURE: {r['deck']}: {r['error']} "
                  f"(furthest_t_s={r.get('furthest_t_s')})")

    return 0 if not any_failure else 0  # failures are data, per the brief; report, don't fail the runner


def build_arg_parser():
    p = argparse.ArgumentParser(
        description="I-099: integrate the four scalar/arrhenius x 700K/900K "
                     "polymer decks on an explicit time grid and compare "
                     "condensed-phase mass loss between temperatures within "
                     "each parameterisation.")
    p.add_argument("--decks", nargs="+", default=DEFAULT_DECKS,
                    help="deck directories to process (default: the four "
                         "canonical i099 decks)")
    p.add_argument("--budget-s", type=float, default=DEFAULT_BUDGET_S,
                    help="wall-clock integration budget per deck, seconds "
                         "(default: 720 = 12 minutes)")
    p.add_argument("--n-points", type=int, default=DEFAULT_N_POINTS,
                    help="number of log-spaced grid points from 1e-6s to 1.0s "
                         "(t=0 is always prepended)")
    p.add_argument("--atol", type=float, default=DEFAULT_ATOL)
    p.add_argument("--rtol", type=float, default=DEFAULT_RTOL)
    p.add_argument("--worker", metavar="DECK_DIR", help=argparse.SUPPRESS)
    return p


def main(argv=None):
    args = build_arg_parser().parse_args(argv)

    if args.worker:
        # Internal subprocess entry point: run exactly one deck in this
        # process and print its result as the last stdout line.
        result = run_one_deck(args.worker, budget_s=args.budget_s,
                               n_points=args.n_points, atol=args.atol, rtol=args.rtol)
        print(RESULT_MARKER + json.dumps(result))
        return 0

    _check_environment()
    results = []
    for deck in args.decks:
        result = _run_worker_subprocess(deck, args.budget_s, args.n_points,
                                         args.atol, args.rtol)
        results.append(result)

    return _print_summary(results)


if __name__ == "__main__":
    sys.exit(main())
