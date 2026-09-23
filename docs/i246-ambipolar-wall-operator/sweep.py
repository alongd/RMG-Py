#!/usr/bin/env python3
"""
I-246 envelope sweep.

Classification rules, conditioning threshold and observables were frozen in
``envelope.md`` and committed BEFORE this script was written or run; see
`git log --diff-filter=A -- docs/i246-ambipolar-wall-operator/envelope.md`.

Every frequency below is normalised PER CHARGED PARTICLE, so nu_ion, nu_wall and
nu_RR are directly comparable. Mixing a per-electron rate with a per-neutral one
is a real error made earlier in this project.

Writes sweep.csv and prints the summary.
"""

import csv
import os
import sys

import numpy as np

import rmgpy.constants as constants

import argon_wall_model as M

KB = constants.kB
HERE = os.path.dirname(os.path.abspath(__file__))

# --- predeclared thresholds, copied from envelope.md section 3 --------------
NEAR_THRESHOLD_RELATIVE = 1.0e-2
CONDITIONING_LIMIT = 1.0e2
ALPHA_CEILING = 1.0e-3


def alpha_rr(te_ev):
    """Radiative recombination Ar+ + e- -> Ar + h*nu, m^3/s.

    ORDER OF MAGNITUDE ONLY (Seaton-type scaling for a singly-charged ion). Its
    role here is solely to establish whether volume recombination can compete
    with the wall; it is not used in any steady state and no conclusion rests on
    its precision."""
    return 2.7e-19 * te_ev ** -0.75


def alpha_3b(te_ev):
    """Collisional (three-body) recombination Ar+ + 2e- -> Ar + e-, m^6/s.
    ORDER OF MAGNITUDE ONLY, same caveat."""
    return 1.0e-38 * te_ev ** -4.5


def wall_constant(te_ev, lam, mu0=M.MU0_AR_IN_AR):
    """K such that nu_wall = K / n_neutral.  K = mu_0 N_0 (k T_e/e) / Lambda^2."""
    from rmgpy.solver.plasma import PLASMA_LOSCHMIDT
    return mu0 * PLASMA_LOSCHMIDT * (constants.R / constants.Na) * (te_ev * M.EV_TO_K) \
        / constants.e / (lam * lam)


def grid_point(te_ev, pressure_pa, radius, length=M.L_NOMINAL, mu0=M.MU0_AR_IN_AR,
               tgas=M.TGAS):
    """Everything envelope.md asks to be reported at one grid point."""
    lam = M.diffusion_length(radius, length)
    n_ar0 = pressure_pa / (KB * tgas)             # un-depleted neutral density
    k_iz = M.k_ionisation(te_ev)
    K = wall_constant(te_ev, lam, mu0)

    nu_ion = k_iz * n_ar0                          # per electron
    nu_wall = K / n_ar0                            # per charged particle
    dnu = nu_ion - nu_wall
    rel = abs(dnu) / max(nu_ion, nu_wall)

    row = dict(te_ev=te_ev, pressure_torr=pressure_pa / M.TORR_TO_PA, radius_m=radius,
               lambda_m=lam, n_ar0=n_ar0, k_iz=k_iz,
               nu_ion=nu_ion, nu_wall=nu_wall, dnu=dnu, rel_dnu=rel,
               p_lambda_torr_cm=(pressure_pa / M.TORR_TO_PA) * lam * 100.0)

    # --- the two branches ---------------------------------------------------
    # (a) sub-threshold, source-sustained: S + nu_ion*n_e = nu_wall*n_e.
    #     Exists only where nu_wall > nu_ion. Reported as the INTERVAL the
    #     cosmic-ray source is known to, never narrowed.
    if dnu < 0.0:
        s_lo = M.cosmic_ray_source(pressure_pa, tgas, high=False)
        s_hi = M.cosmic_ray_source(pressure_pa, tgas, high=True)
        row['ne_sub_lo'] = s_lo / (-dnu)
        row['ne_sub_hi'] = s_hi / (-dnu)
        row['alpha_sub_hi'] = row['ne_sub_hi'] / n_ar0
    else:
        row['ne_sub_lo'] = row['ne_sub_hi'] = row['alpha_sub_hi'] = float('nan')

    # (b) above threshold, the constant-pressure stationary state. With the
    #     composition-dependent mobility nu_wall = K/n_Ar, the balance
    #     k_iz*n_Ar = K/n_Ar gives n_Ar* = sqrt(K/k_iz) -- a SQUARE ROOT, not the
    #     linear n_Ar* = nu_wall/k_iz a constant-mobility model would give. The
    #     constant-pressure EOS then fixes n_e as a DIFFERENCE of large numbers.
    n_ar_star = np.sqrt(K / k_iz) if k_iz > 0 else np.inf
    row['n_ar_star'] = n_ar_star
    f_ion = 1.0 - n_ar_star / n_ar0                # ionisation fraction at the state
    row['f_star'] = f_ion
    te_k = te_ev * M.EV_TO_K
    row['ne_star'] = ((pressure_pa / KB) - n_ar_star * tgas) / (tgas + te_k) \
        if 0.0 < f_ion < 1.0 else float('nan')
    row['alpha_star'] = row['ne_star'] / n_ar_star if 0.0 < f_ion < 1.0 else float('nan')

    # --- conditioning, on the branch, never across the boundary -------------
    # C_wall = |dln n_e/dln K| = (1-f)/(2f); the 1/2 is the mobility's own
    # 1/n_neutral dependence partly self-correcting, and is absent from a
    # constant-mobility model. C_Te carries the extra |1 - dln k_iz/dln Te|.
    if 0.0 < f_ion < 1.0:
        row['C_wall'] = (1.0 - f_ion) / (2.0 * f_ion)
        h = 1e-4
        dlnk_dlnte = (np.log(M.k_ionisation(te_ev * (1 + h)))
                      - np.log(M.k_ionisation(te_ev * (1 - h)))) / (2 * h)
        row['dlnk_dlnTe'] = dlnk_dlnte
        row['C_Te'] = abs(row['C_wall'] * (1.0 - dlnk_dlnte))
        row['C_Ar'] = row['C_wall'] * 2.0     # |dln n_e/dln n_Ar| has no 1/2
    else:
        row['C_wall'] = row['C_Te'] = row['C_Ar'] = float('nan')
        row['dlnk_dlnTe'] = float('nan')

    # --- volume recombination, for comparison only --------------------------
    # nu_RR depends on n_e, which is exactly the ill-determined quantity, so it
    # is quoted at the branch density where one exists and flagged otherwise.
    ne_ref = row['ne_sub_hi'] if dnu < 0.0 else row['ne_star']
    row['ne_ref'] = ne_ref
    if np.isfinite(ne_ref) and ne_ref > 0:
        row['nu_RR'] = alpha_rr(te_ev) * ne_ref
        row['nu_3b'] = alpha_3b(te_ev) * ne_ref ** 2
    else:
        row['nu_RR'] = row['nu_3b'] = float('nan')

    # --- classification, by the predeclared rules ---------------------------
    ill = (rel < NEAR_THRESHOLD_RELATIVE) or any(
        np.isfinite(row[c]) and row[c] > CONDITIONING_LIMIT
        for c in ('C_wall', 'C_Te', 'C_Ar'))
    if np.isfinite(row['alpha_star']) and row['alpha_star'] > ALPHA_CEILING and dnu > 0:
        row['class'] = 'volume-growth (steady state OUT OF SUPPORT)'
    elif ill:
        row['class'] = 'near-threshold/ill-conditioned'
    elif dnu > 0:
        row['class'] = 'volume-growth'
    else:
        row['class'] = 'extinction'
    return row


def main():
    print("plasma module: {0}".format(M.assert_provenance()))
    te_grid = [0.5, 0.7, 0.8, 0.85, 0.9, 0.92, 0.93, 0.95, 1.0, 1.2, 1.5, 2.0, 3.0, 4.0]
    p_grid_torr = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
    r_grid = [0.005, 0.01, 0.02, 0.05, 0.10, 0.20]

    rows = []
    for te in te_grid:
        for p in p_grid_torr:
            for r in r_grid:
                rows.append(grid_point(te, p * M.TORR_TO_PA, r))

    fields = list(rows[0].keys())
    with open(os.path.join(HERE, 'sweep.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print("wrote sweep.csv: {0} grid points x {1} columns".format(len(rows), len(fields)))

    # ---------------------------------------------------------------- summary
    print("\n" + "=" * 100)
    print("CLASSIFICATION CENSUS (rules frozen in envelope.md before this ran)")
    print("=" * 100)
    from collections import Counter
    census = Counter(r['class'] for r in rows)
    for k, v in sorted(census.items(), key=lambda kv: -kv[1]):
        print("  {0:6d}  {1:.1f}%   {2}".format(v, 100.0 * v / len(rows), k))
    both_sides = sum(1 for r in rows if r['dnu'] > 0), sum(1 for r in rows if r['dnu'] < 0)
    print("\n  BOTH SIDES SAMPLED: {0} points with Delta_nu > 0, {1} with Delta_nu < 0".format(
        *both_sides))
    assert both_sides[0] > 0 and both_sides[1] > 0, "the sweep sampled only one side"

    print("\n" + "=" * 100)
    print("THE SUSTAINMENT BOUNDARY  (Delta_nu = 0), located by bisection in Te")
    print("=" * 100)
    print("{0:>10} {1:>9} {2:>11} {3:>12} {4:>13} {5:>13}".format(
        "p (torr)", "R (m)", "Lambda(cm)", "p*Lambda", "Te_thr (eV)", "nu at thr"))
    print("-" * 100)
    boundary = []
    for p in p_grid_torr:
        for r in r_grid:
            lo, hi = 0.3, 6.0
            if grid_point(lo, p * M.TORR_TO_PA, r)['dnu'] > 0:
                continue
            if grid_point(hi, p * M.TORR_TO_PA, r)['dnu'] < 0:
                continue
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                if grid_point(mid, p * M.TORR_TO_PA, r)['dnu'] < 0:
                    lo = mid
                else:
                    hi = mid
            g = grid_point(lo, p * M.TORR_TO_PA, r)
            boundary.append((p, r, g))
            print("{0:>10.2f} {1:>9.3f} {2:>11.4f} {3:>12.4f} {4:>13.5f} {5:>13.4e}".format(
                p, r, g['lambda_m'] * 100, g['p_lambda_torr_cm'], lo, g['nu_wall']))

    print("\nSIMILARITY: nu_ion ~ p and nu_wall ~ 1/(p*Lambda^2), so the boundary should")
    print("collapse onto a single curve in the product p*Lambda. Sorted by p*Lambda:")
    print("{0:>14} {1:>14} {2:>10} {3:>9}".format("p*Lambda(torr cm)", "Te_thr (eV)", "p(torr)", "R(m)"))
    print("-" * 100)
    for p, r, g in sorted(boundary, key=lambda b: b[2]['p_lambda_torr_cm'])[:12]:
        lo, hi = 0.3, 6.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if grid_point(mid, p * M.TORR_TO_PA, r)['dnu'] < 0:
                lo = mid
            else:
                hi = mid
        print("{0:>14.4f} {1:>14.5f} {2:>10.2f} {3:>9.3f}".format(
            g['p_lambda_torr_cm'], lo, p, r))

    print("\n" + "=" * 100)
    print("VOLUME RECOMBINATION AGAINST THE WALL (order of magnitude)")
    print("=" * 100)
    finite = [r for r in rows if np.isfinite(r['nu_RR'])]
    if finite:
        ratios = [r['nu_RR'] / r['nu_wall'] for r in finite]
        ratios3 = [r['nu_3b'] / r['nu_wall'] for r in finite if np.isfinite(r['nu_3b'])]
        print("  nu_RR / nu_wall  over {0} points: min {1:.3e}  median {2:.3e}  max {3:.3e}".format(
            len(ratios), min(ratios), float(np.median(ratios)), max(ratios)))
        print("  nu_3b / nu_wall  over {0} points: min {1:.3e}  median {2:.3e}  max {3:.3e}".format(
            len(ratios3), min(ratios3), float(np.median(ratios3)), max(ratios3)))
        print("  -> the wall is the dominant charged-particle sink wherever it dominates;")
        print("     that is what makes a wall operator the right object to build.")

    print("\n" + "=" * 100)
    print("CONDITIONING ALONG THE ABOVE-THRESHOLD BRANCH, at 5 torr, R = 0.05 m")
    print("=" * 100)
    print("{0:>8} {1:>12} {2:>13} {3:>12} {4:>12} {5:>12}  {6}".format(
        "Te(eV)", "f_star", "alpha_star", "C_wall", "C_Te", "C_Ar", "class"))
    print("-" * 110)
    for te in [0.9, 0.92, 0.93, 0.95, 1.0, 1.1, 1.2, 1.5, 2.0, 3.0]:
        g = grid_point(te, 5.0 * M.TORR_TO_PA, 0.05)
        print("{0:>8.2f} {1:>12.6f} {2:>13.4e} {3:>12.4e} {4:>12.4e} {5:>12.4e}  {6}".format(
            te, g['f_star'], g['alpha_star'], g['C_wall'], g['C_Te'], g['C_Ar'], g['class']))

    print("\n" + "=" * 100)
    print("THE SUB-THRESHOLD SOURCE-SUSTAINED BRANCH, at 5 torr, R = 0.05 m")
    print("=" * 100)
    print("It exists only where nu_wall > nu_ion, and its width is the cosmic-ray")
    print("source interval, which is NOT narrowed:")
    print("{0:>8} {1:>13} {2:>13} {3:>13} {4:>13} {5:>11}".format(
        "Te(eV)", "nu_ion", "nu_wall", "n_e low", "n_e high", "alpha high"))
    print("-" * 100)
    for te in [0.5, 0.7, 0.8, 0.85, 0.9, 0.92]:
        g = grid_point(te, 5.0 * M.TORR_TO_PA, 0.05)
        if np.isfinite(g['ne_sub_hi']):
            print("{0:>8.2f} {1:>13.4e} {2:>13.4e} {3:>13.4e} {4:>13.4e} {5:>11.3e}".format(
                te, g['nu_ion'], g['nu_wall'], g['ne_sub_lo'], g['ne_sub_hi'],
                g['alpha_sub_hi']))
    print("\nThe ratio n_e(high)/n_e(low) is exactly the source interval's width, 5.0,")
    print("at every point: the absolute sub-threshold density is known no better than")
    print("the background ionisation rate is.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
