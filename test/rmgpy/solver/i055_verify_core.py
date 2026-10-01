#!/usr/bin/env python3
"""Post-fix numeric checks: cone invariance, mass ledger, cone detector.

Run under rmg_env with THIS worktree first on PYTHONPATH. Exits 0 only if
every check passes. Every number is recomputed here; nothing is restated.
"""

import contextlib
import io
import logging
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from unzip_realizability_probe import (  # noqa: E402
    MU0, MU1, MU2, MONOMER, build_system,
)

T_END = 100.0
K_SWEEP = (0.05, 1.0, 100.0)
# Coupled arm mirrors the poly_104 deck shape: scission feeds mu0 while unzip
# drains mu1. k_unzip = 100 sits far ABOVE the measured k_scission/4 = 0.25
# negative-mu1 threshold, so this arm is discriminating rather than vacuous.
K_SCISSION_ARM = 1.0
IC_ARM = (0.099, 0.7378, 8.0)


def sweep(k_unzip, k_scission, moments, t_end=T_END, atol=1e-14):
    rs = build_system(k_unzip, moments=moments, k_scission=k_scission)
    y0 = rs.y.copy()
    y0[MU0], y0[MU1], y0[MU2] = moments
    zeros = np.zeros_like(y0)

    def rhs(t, y):
        return rs.residual(t, y, zeros)[0]

    sol = solve_ivp(rhs, (0.0, t_end), y0, method="LSODA", rtol=1e-10,
                    atol=atol, max_step=t_end / 2000.0)
    if sol.status != 0:
        raise RuntimeError(f"integration failed: {sol.message}")
    mu0, mu1 = sol.y[MU0], sol.y[MU1]
    drained = moments[1] - mu1[-1]
    released = sol.y[MONOMER][-1] - sol.y[MONOMER][0]
    return {
        "min_gap": float((mu1 - mu0).min()),
        "min_mu0": float(mu0.min()),
        "min_mu1": float(mu1.min()),
        "drained": float(drained),
        "released": float(released),
        "defect": float(drained - released),
        "scale": float(max(abs(moments[1]), abs(drained))),
    }


ATOL = 1e-14
# LSODA controls LOCAL error per step at atol; error accumulated over a
# ~100 s integration with thousands of steps is a small multiple of that.
# K = 10 is that multiple. It is not a free knob: check_error_scaling below
# refutes the alternative reading -- that a residual excursion is a genuine
# cone exit -- by showing the excursion tracks atol linearly rather than
# converging to a fixed physical value. This is a tolerance on THIS probe's
# scipy integration, not on any solver detector, floor or product tolerance.
K_ATOL = 10.0


def check_error_scaling():
    """The worst case must be integration error, not a real cone exit.

    Discriminator: integration error shrinks with atol; a genuine excursion
    out of the cone converges to a fixed nonzero value. Runs the tightest
    case (scission arm at k_unzip = 100, which fully exhausts the pool) at
    three atol values three decades apart."""
    print("  atol refinement on the worst case (with scission, k_unzip=100):")
    vals = []
    for atol in (1e-12, 1e-14, 1e-16):
        r = sweep(100.0, K_SCISSION_ARM, IC_ARM, atol=atol)
        vals.append((atol, r["min_gap"], r["min_mu0"]))
        print(f"    atol={atol:<8g} min(mu1-mu0)={r['min_gap']:+.6e}  "
              f"min(mu0)={r['min_mu0']:+.6e}  "
              f"|gap|/atol={abs(r['min_gap']) / atol:.2f}")
    worst = max(abs(g) for _, g, _ in vals)
    # If the excursion were physical it would be ~constant across three
    # decades of atol; require instead that every excursion sits within
    # K_ATOL of its OWN atol, i.e. that it tracks the tolerance.
    tracks = all(abs(g) <= K_ATOL * a and abs(m) <= K_ATOL * a
                 for a, g, m in vals)
    if tracks:
        print(f"    OK: every excursion is within {K_ATOL:g}*atol of its own "
              f"atol -- it tracks the integration tolerance across three "
              f"decades, so it is integration error, not a cone exit "
              f"(a real exit would be atol-independent).")
    else:
        print(f"    FAIL: excursions do not track atol (worst {worst:.3e}); "
              f"this may be a genuine cone exit, NOT integration error.")
    return tracks


def check_cone_and_mass():
    """Verifier items 2 and 3."""
    ok = True
    tol = K_ATOL * ATOL
    for label, k_s, ic in (("unzip only ", 0.0, (1.0, 5.0, 30.0)),
                           ("with scission", K_SCISSION_ARM, IC_ARM)):
        for k in K_SWEEP:
            r = sweep(k, k_s, ic)
            gap_ok = r["min_gap"] >= -tol
            mu0_ok = r["min_mu0"] >= -tol
            rel_defect = abs(r["defect"]) / max(r["scale"], 1e-300)
            mass_ok = rel_defect <= 1e-9
            ok &= gap_ok and mu0_ok and mass_ok
            print(f"  [{label}] k_unzip={k:<8g} t_end={T_END:g}s  "
                  f"min(mu1-mu0)={r['min_gap']:+.6e}  "
                  f"min(mu0)={r['min_mu0']:+.6e}  "
                  f"{'OK' if gap_ok and mu0_ok else 'FAIL'}")
            print(f"  {'':>{len(label) + 4}}mass: drained={r['drained']:.12e}  "
                  f"released={r['released']:.12e}  "
                  f"rel defect={rel_defect:.3e}  "
                  f"{'OK' if mass_ok else 'FAIL'}")
    print(f"  (cone/positivity tolerance {tol:g} = {K_ATOL:g} * solver atol "
          f"{ATOL:g}; mass tolerance 1e-9 relative)")
    ok &= check_error_scaling()
    return ok


def check_cone_detector():
    """The cone census must fire on a cone-exiting ACCEPTED state that is
    NOT negative -- exactly the sub-threshold blind spot -- and must stay
    silent on a healthy state. It must never raise; the r81 negative check
    keeps sole ownership of raising."""
    ok = True
    rs = build_system(1.0, moments=(1.0, 5.0, 30.0))

    class _Grab(logging.Handler):
        def __init__(self):
            super().__init__()
            self.msgs = []

        def emit(self, record):
            self.msgs.append(record.getMessage())

    def run_state(mu0, mu1, mu2, fresh):
        if fresh:
            rs._realizability_warned = set()
        rs.y[MU0], rs.y[MU1], rs.y[MU2] = mu0, mu1, mu2
        h = _Grab()
        root = logging.getLogger()
        root.addHandler(h)
        old = root.level
        root.setLevel(logging.WARNING)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                rs._assert_pool_moments_accepted()
        finally:
            root.removeHandler(h)
            root.setLevel(old)
        return [m for m in h.msgs if "CONE CENSUS" in m]

    # Healthy, well inside the cone: silent.
    hits = run_state(1.0, 5.0, 30.0, fresh=True)
    print(f"  healthy state (mu0=1, mu1=5): {len(hits)} cone-census line(s) "
          f"-- {'OK' if not hits else 'FAIL (false positive)'}")
    ok &= not hits

    # Cone exit with mu1 STILL POSITIVE: the sub-threshold blind spot. The
    # r81 negative check cannot see this state; the cone census must.
    hits = run_state(1.0, 0.5, 30.0, fresh=True)
    print(f"  cone exit, mu1=0.5 > 0 (r81 negative check is blind here): "
          f"{len(hits)} cone-census line(s) -- "
          f"{'OK' if len(hits) == 1 else 'FAIL'}")
    if hits:
        print(f"    {hits[0]}")
    ok &= len(hits) == 1

    # Warn-once: same pool, same violation, no second line.
    hits2 = run_state(1.0, 0.4, 30.0, fresh=False)
    print(f"  repeat violation, same pool: {len(hits2)} cone-census line(s) "
          f"-- {'OK (warn-once)' if not hits2 else 'FAIL (log spam)'}")
    ok &= not hits2

    # The r81 negative raise is untouched.
    rs._realizability_warned = set()
    rs.y[MU0], rs.y[MU1], rs.y[MU2] = -1.0, -2.0, -3.0
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            rs._assert_pool_moments_accepted()
    except ValueError as e:
        good = "beyond the exhaustion floor" in str(e)
        print(f"  r81 negative raise still fires: "
              f"{'OK' if good else 'FAIL (wrong message)'}")
        ok &= good
    else:
        print("  r81 negative raise still fires: FAIL (did not raise)")
        ok = False
    return ok


def main():
    print("--- item 2/3: cone invariance and mass ledger, post-fix ---")
    ok = check_cone_and_mass()
    print("--- added detector: accepted-state cone census ---")
    ok &= check_cone_detector()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
