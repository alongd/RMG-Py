#!/usr/bin/env python3
"""Measure, on the PRE-FIX kernel, the two distinct k_unzip crossings.

They are not the same crossing and only one of them has a threshold:

  * CONE EXIT (mu1 < mu0) -- leaving the realizable set. At the boundary
    mu1 = mu0 the vector field reads dmu0/dt = k_s*(mu1 - mu0) = 0 and
    dmu1/dt = -k_u*mu0 < 0, so the boundary is crossed outward by ANY
    trajectory that reaches it, at every k_unzip > 0. No threshold.
  * NEGATIVE mu1 -- what _assert_pool_moments_accepted actually fires on.
    The (mu0, mu1) subsystem is linear with matrix [[-k_s, k_s], [-k_u, 0]];
    its discriminant is k_s*(k_s - 4*k_u), so the eigenvalues are complex --
    the trajectory spirals through mu1 = 0 -- iff k_u > k_s/4, and are real
    and negative (no zero crossing) below that. Threshold at k_s/4.

The consequence is the detector gap: below k_s/4 the pool leaves the cone
and nothing fires.

Bisection below measures the threshold rather than assuming it.
"""

import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from unzip_realizability_probe import MU0, MU1, MU2, build_system  # noqa: E402


K_SCISSION = 1.0
IC = (0.099, 0.7378, 8.0)
T_END = 400.0


def run(k_unzip, t_end=T_END):
    rs = build_system(k_unzip, moments=IC, k_scission=K_SCISSION)
    y0 = rs.y.copy()
    y0[MU0], y0[MU1], y0[MU2] = IC
    zeros = np.zeros_like(y0)

    def rhs(t, y):
        return rs.residual(t, y, zeros)[0]

    def cone(t, y):
        return y[MU1] - y[MU0]

    def neg(t, y):
        return y[MU1]

    cone.terminal = False
    cone.direction = -1
    neg.terminal = True
    neg.direction = -1
    sol = solve_ivp(rhs, (0.0, t_end), y0, method="LSODA", rtol=1e-11,
                    atol=1e-16, events=(cone, neg), max_step=t_end / 4000.0)
    t_cone = float(sol.t_events[0][0]) if len(sol.t_events[0]) else None
    t_neg = float(sol.t_events[1][0]) if len(sol.t_events[1]) else None
    return t_cone, t_neg


def main():
    print("=" * 78)
    print(f"PRE-FIX threshold measurement: k_scission={K_SCISSION}, "
          f"(mu0, mu1, mu2)={IC}")
    print("=" * 78)

    print("cone exit (mu1 < mu0) vs negative mu1, swept over 5 decades:")
    cone_at_all = True
    for k in (1e-3, 1e-2, 0.05, 0.1, 0.2, 0.25, 0.3, 1.0, 10.0, 100.0):
        t_cone, t_neg = run(k)
        cone_at_all &= t_cone is not None
        print(f"  k_unzip={k:<8g} k_u/k_s={k / K_SCISSION:<8g} "
              f"cone t={('%.5g' % t_cone) if t_cone is not None else 'NONE':<10} "
              f"mu1<0 t={('%.5g' % t_neg) if t_neg is not None else 'NONE'}")

    def bisect(t_end):
        """Bisect on the PRESENCE of a mu1 = 0 crossing within t_end."""
        lo, hi = 0.05, 1.0
        assert run(lo, t_end)[1] is None, "lower bracket already crosses zero"
        assert run(hi, t_end)[1] is not None, "upper bracket does not cross zero"
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if run(mid, t_end)[1] is None:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    predicted = K_SCISSION / 4.0

    # Just above the threshold the mu1 = 0 crossing is TANGENTIAL: the first
    # zero of the spiral occurs at a time that diverges as k_u -> k_s/4 from
    # above. A finite horizon therefore MISSES crossings in a shrinking band
    # above the threshold, biasing the bisected value HIGH by an amount that
    # must shrink as the horizon grows. Measuring at two horizons turns that
    # from an unexplained discrepancy into a checked prediction -- rather
    # than a tolerance chosen to make one number pass.
    print("-" * 78)
    print("bisected threshold vs integration horizon "
          "(bias must be positive and shrinking):")
    measurements = []
    for t_end in (100.0, 400.0, 1600.0):
        m = bisect(t_end)
        measurements.append((t_end, m))
        print(f"  horizon={t_end:<8g} measured k_unzip threshold = {m:.10g}   "
              f"bias = {m - predicted:+.3e}  ({(m - predicted) / predicted:+.3e} rel)")
    measured = measurements[-1][1]

    print("-" * 78)
    print(f"MEASURED mu1<0 threshold : k_unzip = {measured:.10g}  "
          f"(longest horizon, {measurements[-1][0]:g} s)")
    print(f"PREDICTED (k_scission/4) : k_unzip = {predicted:.10g}")
    print(f"relative difference      : {abs(measured - predicted) / predicted:.3e}")
    print(f"cone exit occurred at every swept k_unzip: {cone_at_all}")
    print("-" * 78)

    ok = True
    biases = [m - predicted for _, m in measurements]
    if not all(b > 0 for b in biases):
        print(f"FAIL: bias is not uniformly positive ({biases}) -- the "
              f"finite-horizon explanation does not hold")
        ok = False
    elif not all(a > b for a, b in zip(biases, biases[1:])):
        print(f"FAIL: bias does not shrink with horizon ({biases}) -- the "
              f"finite-horizon explanation does not hold")
        ok = False
    else:
        print("OK: bias is positive and shrinks monotonically with the "
              "integration horizon, which is the signature of a tangential "
              "crossing being missed -- not a disagreement with k_scission/4.")
    if abs(measured - predicted) / predicted > 1e-3:
        print("FAIL: measured threshold disagrees with k_scission/4 by >1e-3")
        ok = False
    else:
        print(f"OK: measured threshold {measured:.8g} agrees with "
              f"k_scission/4 = {predicted:g} to "
              f"{abs(measured - predicted) / predicted:.2e} relative")
    if not cone_at_all:
        print("FAIL: cone exit was NOT structural across the sweep")
        ok = False
    else:
        print("OK: cone exit is STRUCTURAL -- present at every k_unzip > 0, "
              "including far below the k_scission/4 negative-mu1 threshold, "
              "where the accepted-state negative check never fires.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
