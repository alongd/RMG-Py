###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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

"""
Biphasic hybrid reactor solver for polymer pyrolysis/degradation in RMG-Py.

Phase Policy: "Produce-then-Transfer"
1. Reactions are strictly phase-pure.
   - Gas reactants -> Gas products.
   - Polymer reactants -> Polymer products.
   - Cross-phase core products disqualify a reaction (rate=0).
   - Edge products (not in solver state) allowed for diagnostics.
2. Mass Transfer is explicit via kLa driving force.

Hybrid Representation:
- Explicit Oligomers: Species tracking chains with degree of polymerization (DP) <= xs.
- Moment Tail: Statistical moments (mu0, mu1, mu2) for chains with DP > xs.
- Hybrid Handshake: Dynamic flux from Tail -> Explicit at boundary xs.

Units:
- Time: s
- Residual: mol/s (dn/dt - dydt)
- Moments (Input/State): Moles of Moments (Mu_k = mu_k * V_poly)
- Tail kinetics interface: tail_kinetics receives μk/V_poly [mol/m^3] and returns d(μk/V_poly)/dt [mol/m^3/s]
  (and small-species sources in mol/m^3/s). The solver multiplies by V_poly when accumulating into dn_dt [mol/s].
- Diagnostics: Core rates in mol/m^3_phase/s (or just 'rate' units)
"""

from __future__ import annotations

import copy
import itertools
import logging
import math
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

cimport cython
from libc.math cimport exp as _c_exp
from libc.math cimport fabs
import numpy as np
cimport numpy as np

import rmgpy.constants as constants
from rmgpy.quantity import Quantity
from rmgpy.solver.base import ReactionSystem


SMALL_EPS = 1e-30


# M18.2d (round-40 adjudicated one-shot): C1 extension of the max(0, y)
# state clamp into NEGATIVE trial-state territory ONLY.
#   v >= 0:  returns v EXACTLY (bit-for-bit the old clamp -- all valid
#            nonnegative states, and every place a pool can PARK, see the
#            identity law untouched; no epsilon, no multiplier);
#   v <  0:  v * exp(v / lam) -- slope 1 at 0- (C1 joint with the v > 0
#            side), bounded extremum -lam/e at v = -lam, decaying to the
#            old clamped 0 for deeply negative v.
# lam is the r81 floor scale max(SMALL_EPS, 100*atol), used purely as the
# DECAY LENGTH in v < 0 territory (the region only BDF predictor
# overshoots and Newton trial states visit).
#
# MECHANISM (M18.2c/d forensics, rounds 38-41): every pathologically
# parked pool (poly_102 mod_5 at 0.2-5 floors) sits within one BDF
# predictor overshoot of the {y_i = 0} kink planes of its own moment
# slots, and measured Newton trial states cross them (trial
# mu0 = -1.7e-10 against accepted ~2e-11). At a negative trial base the
# hard clamp flattened every row fed by that slot, so one-sided FD
# Jacobian columns were EXACTLY ZERO -- a bistable Newton matrix that
# flipped between zero and full columns across trial sequences, at every
# park. Gas radicals in QSS (H at 1e-9..1e-12 mol) contributed the same
# structure through the concentration clamp. The C1 extension gives the
# corrector a bounded restoring derivative instead of a blind zero
# column. Battery evidence (round-40): both standing IDID=-7 deaths
# eliminated (exact-crash rtol=1e-6 at t=24.639; from-deck rtol=1e-4 at
# t=14.2445, first-ever full 0-100 s traversal) with ~3.4x wall speedup
# on the deck-tolerance crash window; K2 accuracy fixtures bit-for-bit
# unchanged (the kink field was never their mechanism).
#
# SITE INVENTORY (of 48 direct max(0, y[...]) state reads): 41 VALUE
# sites carry the softclamp -- the vectorized C_gas/C_poly concentration
# reads (the widest field: one kink plane per species), the dispatch
# site scalings and throttle pairs (forward/reverse/chip/VE a>0/a<0/
# legacy), scission/channel/QSSA kernel moment reads, the limiter/cone
# machinery state reads, and their get_reaction_rates mirrors (accepted
# states can park in-band-negative, and reported rates must never
# diverge from the residual's law). 7 sites keep the hard clamp:
# bookkeeping (_explicit_moment_contributions, condensed-mass, SGH
# census/tripwire reads on accepted states) and branch-guard
# COMPARISONS (equivalent branch structure either way, and piecewise-
# constant selectors are not differentiated). Two pre-existing normative
# structures deliberately re-flatten parts of the extension downstream:
# pool-moment reads inside the limiter/cone machinery hit its internal
# <= 0 guards (rounds 29/30 law), and the section-4 moment dispatch is
# gated on rf > 0 / rr > 0, so MOMENT rows ignore reversed flux while
# GAS-SPECIES rows carry the full smooth extension -- the battery shows
# the gas-row field is the load-bearing one.
cdef inline double _softclamp(double v, double lam):
    if v >= 0.0:
        return v
    return v * _c_exp(v / lam)
# Attribution trust floor multiplier (item #14a, spec 2026-06-11 §2/§3):
# trust = max(SMALL_EPS, ATTRIBUTION_TRUST_K * atol_mu0). K = 100 buys ~two
# decades of separation between the integrator's own error budget (the atol
# on the pool's mu0 slot) and what the spawn-gate snapshot will TRUST as an
# E[n] denominator — a margin decision on record (like the thermo bounds),
# not a default. Tolerance-anchored on purpose: any absolute mole constant
# is wrong across deck scales (a ~10 mg TGA sample at Mn ~1e4-1e5 g/mol puts
# a whole physical deck's mu0 at ~1e-10-1e-9 mol). SMALL_EPS itself is
# untouched and keeps its job as the solver's REALIZABILITY clamp; ONLY the
# gate snapshot's attribution distrusts the band between them — the gate
# deliberately distrusts what the solver still legally computes (two
# constants, two jobs).
ATTRIBUTION_TRUST_K = 100.0
# Exhaustion-tail census/tripwire floor multiplier (adjudicated round 81,
# the PP run-7 IDID=-7 wall; re-adjudicated round 86: DIAGNOSTIC ONLY --
# the r81 RHS read-projection is reverted, the RHS always computes from the
# raw state): per-state floor = max(SMALL_EPS, K * atol[state]).
# A pool is counted EXHAUSTED for the census only when |mu0|, |mu1| AND
# |mu2| ALL sit below their floors (never mu0 alone -- tiny mu0 with
# nontrivial mu1 is a few long chains, not an empty pool). K = 100 is the
# same two-decade tolerance-anchored margin as ATTRIBUTION_TRUST_K (a
# separate constant on purpose: the spawn-gate trust band and the solver
# exhaustion band are two jobs that may diverge). Any moment more negative
# than -floor is a HARD error (integrator corruption), never
# max(..., SMALL_EPS). Generic solver infrastructure OUTSIDE the 2.8 kernel
# contract -- the sidecar recipe strings are deliberately untouched.
EXHAUSTION_FLOOR_K = 100.0
# I-244 ACCEPTED-STATE RELATIVE FLOOR. The r81 hard error (an ACCEPTED moment
# more negative than its floor is integrator corruption) and the
# pool_exhausted event read, per pool and moment order k,
#     F_k = max(EXHAUSTION_FLOOR_K*atol[slot], EXHAUSTION_RELATIVE_K*rtol[slot]
#               * peak_k),   peak_k = max |mu_k| over the run's accepted states
# (mu_k(0) for a pool that only drains). The absolute term alone (1e-14 mol at
# the default atol) refused a pool that depolymerizes completely on round-off:
# the i061 error-weight floor makes 100*atol the moment slots' own DASPK atol,
# so the floor sat at ONE error weight, and the RMS error norm lets one slot of
# an N-slot state carry up to sqrt(N) of them (mu2 = -2e-14 mol from
# mu2(0) = 2e3 mol on an accepted state). The relative term is DASPK's own
# error budget at the pool's scale: while the pool held ~peak_k, every accepted
# step was allowed a local error of rtol*peak_k in that slot, and a drain to
# zero does not take that error back, so a residue inside rtol*peak_k cannot
# be told from zero. K = 1: exactly the integrator's budget, no margin -- a
# truly negative state (a sign error drains O(peak_k) past zero) crosses it
# after a fraction rtol of its excursion. Nothing is fitted. The absolute floors
# (_pool_mu_floors) and every consumer of them -- the trial-state census, the
# cone/variance censuses, the bundle limiter, the softclamp -- are untouched,
# so F_k only ever relaxes the hard error, never tightens it.
EXHAUSTION_RELATIVE_K = 1.0
# i061 MOMENT ERROR-WEIGHT FLOOR (DASPK step-collapse fix). The moment-
# coordinate state slots (a pool's mu0/mu1/mu2 positions and any is_moment_dummy
# core species) are BOOKKEEPING coordinates carrying dmu/dt, not molar species
# amounts. The chemistry already declares any moment at |mu_k| <= EXHAUSTION_
# FLOOR_K*atol to be AT THE FLOOR -- indistinguishable from zero for the r81
# negative-moment tripwire, the exhaustion census, the softclamp and the
# attribution trust band. Yet with the deck atol (e.g. 1e-14) DASPK's error
# weight ewt = rtol*|y| + atol still demands the integrator resolve those slots
# to atol, i.e. TWO DECADES BELOW the model's own accepted-state floor. At the
# all-monomer moment boundary ~43/70 slots sit there at once and DDASPK
# collapses its step to ~1e-9 s chasing noise the chemistry discards (i061:
# ~40-55 days to reach terminationTime). We therefore floor the DASPK atol of
# the MOMENT slots ONLY to the model's own accepted-state floor, so error
# control stops one decade before the chemistry does. Anchored to EXHAUSTION_
# FLOOR_K's rationale but a SEPARATE constant per the idiom above (error control
# and chemistry exhaustion are two jobs that may diverge). NEVER applied to a
# physical species -- see the include-mask scoping in initialize_model.
MOMENT_EWT_FLOOR_K = 100.0
# Near-exhaustion bundle limiter band (P1-A, round-27 re-adjudication of the
# P1-1 hard-min): the cross-pool bundle cap S_cap is TAIL-ONLY. The band is
# measured by the debited pool's ACCEPTED-STATE FLOOR DISTANCE
#     E = min(mu0/f0, mu1/f1, mu2/f2),   f_k = self._pool_mu_floors[p, k]
#                                            = max(SMALL_EPS, 100*atol[state])
# i.e. E is DIMENSIONLESS: how many r81 accepted-state floors the pool's
# thinnest debited moment sits above zero. f_k is the SAME per-state floor
# the r81 negative-moment tripwire and the exhaustion census use, so the
# limiter is deck-independent (no absolute mole constant) and scale-aware
# (tracks the integrator's own atol basis), and it is guaranteed to be fully
# active BEFORE any moment can approach the r81 floor (-f_k) -- no cliff at
# the floor. Band edges (dimensionless floor distances):
#   E <= E_LO (100 floors  = two decades above the census band |mu_k|<=f_k):
#             the soft cap is fully active, S_free = softmin_p(S_base,
#             S_cap terms) -- every debited moment's drain vanishes
#             linearly near exhaustion (softmin_p <= every term);
#   E >= E_HI (1e4 floors = four decades above the floor): S_free == S_base
#             EXACTLY (early return -- bit-for-bit the pre-hard-min bulk
#             law; healthy pools never feel the limiter);
#   between:  C1 smoothstep blend w = 3e^2 - 2e^3 on the normalized
#             e = (E - E_LO)/(E_HI - E_LO),
#             S_free = w*S_base + (1-w)*softmin_p(S_base, S_cap terms).
# S_free then passes through the INDEPENDENT cone-margin drain gate
# (round-30 N1 below) to give the direction's S_eff.
# With the default atol = 1e-16 the band spans ~[1e-12, 1e-10] mol on the
# thinnest moment -- far below any healthy pool, exactly the depletion band
# where the regen-#2 crash lives. (Mirrored in the numpy oracle consumer --
# keep in sync.)
#
# Round-29/30 adjudication (regen-#3 IDID=-7 post-mortem): the
# near-exhaustion bundle limiter (tail-only smoothstep; C1 soft-min,
# cone-margin drain guard) adds two pieces. N2 lives strictly inside the
# tail-only exhaustion structure; N1 is a STRUCTURALLY INDEPENDENT second
# gate with its own band (round-30 re-adjudication: folding S_cone into
# the E blend diluted the guard to w*S_base in-band and IGNORED it
# entirely for out-of-cone pools with bulk-scale moments, E >= E_hi).
# Bulk healthy IN-CONE pools stay bitwise S_base through both gates:
#
# (N1) CONE-MARGIN DRAIN GATE (S_cone), independent of E: a length-biased
#      bundle debit (b1 = mu2/mu1; mu1/mu0 for end-group picks) moves the
#      realizability cone margin Q10 = mu1 - mu0 by -(b1 - b0) per event;
#      when b1 > b0 the drain pushes the pool toward/deeper into the
#      unrealizable mu1 < mu0 cone, where b1 self-accelerates
#      (d(b1)/d(mu1) = -mu2/mu1^2) and the mu3 closure freezes (regen-#3
#      forensics: daughter parked at mu = [1e-3, 1.5e-7, 1.8e-3] carried
#      b1 ~ 1.2e4 and d(b1)/d(mu1) ~ 8e10 into the RHS). For a
#      cone-shrinking debit (b1 > b0):
#        * Q10 <= 0 (pool AT or OUTSIDE the cone): the debited direction
#          returns 0 REGARDLESS of E -- the cone-deepening drain is off
#          even at bulk-scale moments;
#        * Q10 > 0: the gate runs on its own dimensionless margin
#          distance M = Q10 / max(f0, f1) (r81 floor units -- the same
#          deck-independent, scale-aware anchoring as E; degree-1
#          homogeneous in the moments over fixed floors), with band
#          edges M_LO/M_HI mirroring E_LO/E_HI: M >= M_HI returns the
#          exhaustion-limited S_free EXACTLY (margin safely bulk). Below
#          M_HI, form the narrowed completion
#          C = u*softmin_p(S_free, S_cone), where S_cone =
#          Q10/(V_poly*(b1 - b0)) is the event-site density that would
#          spend the whole margin and u is I-090's reverse smoothstep on
#          the unresolved b1 neighbourhood (zero outside it). M <= M_INNER
#          returns C; a C1 smoothstep reaches the original soft-min cap at
#          M_LO, above which the original band blends it into S_free.
#      The gate throttles the EVENT RATE only (per-event moment ratios
#      and mass bookkeeping untouched -- the bundle pick stays
#      length-biased) and is direction-aware exactly like the exhaustion
#      cap (forward debits src, reverse debits dst).
#
# (N2) C1 SOFT-MIN: every in-band hard min() -- the moment-min in E, the
#      per-moment cap min, and min(S_base, S_cap) -- is a RECIPROCAL
#      P-NORM soft-min over positive terms,
#          softmin_p(x_i) = (sum_i x_i^(-p))^(-1/p),  p = SOFTMIN_P
#      (positive, scale-homogeneous of degree 1, C^inf where all terms are
#      positive, softmin_p <= min <= every term; NOT log-sum-exp, which
#      needs a dimensional temperature and can go negative). Quantitative
#      tie bias at p = 8: a two-term tie returns 2^(-1/8) = 0.917*min, a
#      three-term tie 3^(-1/8) = 0.872*min -- a numerical regularization
#      of near-exhaustion / near-cone states, NOT physical kinetics (ties
#      only bind inside the bands; bulk is exact). The hard min()s put C0
#      kinks INSIDE the band wherever an argmin switches, which is what
#      fed DASPK's quasi-Newton the regen-#3 intermittent IDID=-7. The
#      band edges stay C1: w' = 0 at both smoothstep ends and the soft
#      term carries the (1-w) -> 0 factor at E_HI, so the exact early
#      return introduces no derivative cliff. (The single accepted C0
#      kink is at Q10 = 0 exactly, where the N1 hard zero meets the
#      softmin_p(S_free, S_cone) -> 0 limit CONTINUOUSLY with a bounded
#      one-sided slope -- adjudicated round-30, unlike the unbounded
#      hard-min kinks this change removes.)
BUNDLE_LIMITER_E_LO = 1.0e2
BUNDLE_LIMITER_E_HI = 1.0e4
# Cone-margin gate band (round-30 N1): dimensionless margin distance
# M = Q10 / max(f0, f1) = how many r81 accepted-state floors of cone
# margin remain. Same numeric edges as the exhaustion band on purpose
# (both measure "distance to a singular manifold" in floor units), but
# SEPARATE constants -- the two gates are independent jobs and their
# bands may diverge under future adjudication. (Mirrored in the numpy
# oracle consumer -- keep in sync.)
CONE_MARGIN_M_LO = 1.0e2
# I-072 numerical regularization widens I-069's smoothstep from W=1
# (M=99..100) to W=50 accepted-state floors (M=50..100). These are distances
# measured IN floor units, not sub-floor noise. A cubic smoothstep's maximum
# normalized slope is 1.5/W and curvature is 6/W^2, so W=50 is 50x gentler
# in slope and 2500x gentler in curvature than W=1. Instrumentation found a
# daughter parked at M=99.02: W=1 accumulated 173 convergence failures and
# first-order 1.8e-5 s steps; widths 1 and 10 missed the replay budget, while
# W=50 completed it. PE/PS final moments stayed bit-identical. The depletion-
# handoff reversible pin near solverPolymerTest.py's directional-gate test
# deliberately changes from both directions off to its M=93 side live and
# M=46 side off. W=100 (M_inner=0) also completed, but removes the retained
# exact-zero dead band entirely, so it is not interchangeable with W=50.
# M >= M_LO remains bit-for-bit unchanged.
CONE_MARGIN_M_INNER = 5.0e1
CONE_MARGIN_M_HI = 1.0e4
# Soft-min sharpness p (round-29 N2), moderate by design. Chosen by the
# Lipschitz slope-jump sweeps (the 3001-point C1 band sweep, the regen-#3
# broken-state RHS sweep, and a targeted in-band cone-crossing sweep that
# drives mu1 through mu0 inside the band -- the state family whose argmin
# switches carried the old hard-min's C0 kinks): the old hard-min law
# jumps its finite-difference slope by ~0.72 of the slope scale at the
# crossing (an O(1) kink); the soft law measures ~1.8e-3 for p in
# {6, 8, 12} (grid-curvature floor, statistically flat in p) rising to
# ~2.5e-3 at p = 16 as the corner re-sharpens toward the hard min. p = 8
# is the adjudicated moderate choice: on the smooth plateau, while staying
# within ~8.3% of the hard min when two terms tie (2^(-1/8) = 0.917).
# Smaller p over-throttles states with several comparable terms.
BUNDLE_LIMITER_SOFTMIN_P = 8.0
LN_EXP_OVERFLOW_GUARD = 700.0
# N5b mu3-closure boundary band (round-62 DASSL-hang root cause):
# _safe_mu3_from_mu012's realizability guard (mu1 >= mu0) is exact, but
# right at the all-monomer boundary mu1 -> mu0 the two moments are
# catastrophic-cancellation noise of each other (agreeing to 12+ digits,
# per the round-62 investigation), so a solver step landing a rounding
# ULP on the "wrong" side of mu1 == mu0 must not see the guard's full
# 0.0 branch -- the physical limit AT the boundary is mu3 -> mu1 (every
# chain has length 1, so sum(n^3) == sum(n) == mu1), not zero. sqrt(machine
# epsilon) is the standard heuristic scale for "this difference is
# rounding noise, not signal" (same scale already used for FD step
# sizing at line ~5736); band = MU3_CLOSURE_BOUNDARY_REL * mu0 is
# therefore the width within which mu0 - mu1 > 0 is treated as noise
# astride the boundary. Genuine deep violations (mu1 < mu0 - band)
# still return 0.0 unchanged.
MU3_CLOSURE_BOUNDARY_REL = float(np.sqrt(np.finfo(np.float64).eps))
CONE_B1_NOISE_REL_FLOOR = MU3_CLOSURE_BOUNDARY_REL
# Cone-margin dead-band NARROWING (I-090). See _bundle_limited_site's
# stage-2 block for the law; the reasoning about its WIDTH is here.
#
# Inside the dead band the round-62 law returns 0.0 for EVERY b1c. That is
# the right answer where the stage-2 bound S_cone = Q10/(V*(b1c - 1)) is a
# real bound, and the wrong one right at b1c == 1, where S_cone DIVERGES
# and therefore caps nothing: the hard zero there discards a rate this
# gate has no business capping, and does it with a jump (S_free on the
# b1c <= 1 side, 0.0 on the other) that DASPK's corrector cannot step
# across -- the round-62/I-077b never-returns signature.
#
# So the completion is kept on a neighbourhood of the b1c == 1 surface and
# handed back to the exact hard zero outside it. The width of that
# neighbourhood is the scale below which b1c - 1 = (mu2 - mu1)/mu1 is not a
# difference the accepted state ever RESOLVED -- which is just the
# propagated uncertainty of the two moments the difference is taken
# between, in exactly the error-weight form the integrator itself enforces:
#
#     ewt(mu_k) = rtol*|mu_k| + f_k          (DASPK's own per-slot weight,
#                                             f_k the r81 accepted-state
#                                             floor already read here)
#     band      = (ewt(mu1) + ewt(mu2)) / mu1
#               = (rtol*(mu1 + mu2) + f1 + f2) / mu1
#
# The two uncertainties ADD because the numerator is their difference;
# dividing by mu1 converts the result from moment units into b1c units.
# Every input is a scale this function already carries for other reasons --
# the deck's declared rtol, and the SAME r81 per-moment floors the dead
# band itself rests on -- so nothing here is tunable against a test.
# CONE_B1_NOISE_REL_FLOOR is the degenerate-case floor UNDER both terms,
# not the band: sqrt(machine eps), this file's own scale for "this
# difference is rounding noise, not signal" (MU3_CLOSURE_BOUNDARY_REL, used
# identically one moment down at _safe_mu3_from_mu012's mu0 - mu1 guard,
# line ~1176, and for FD step sizing, line ~5736). Nothing is resolvable
# below it however tight the tolerances are asked to be.
#
# The band must be RELATIVE (a width in b1c - 1), not absolute (a width in
# mu2 - mu1), and that is not a stylistic choice -- it is what the divergent
# -bound argument is about. The obvious absolute alternative is seductive:
# measure mu2 - mu1 against the SAME CONE_MARGIN_M_LO floors-count this gate
# applies to Q10 = mu1 - mu0 one line earlier, which introduces no constant
# at all. It is REFUTED, by measurement rather than taste:
# test_bundle_limiter_two_regime_unit_pins drives the TAIL case at
# mu = (1.9, 26, 610)*1e-15 -- mu2 - mu1 = 58 floors, so "noise" by that
# standard, yet b1c - 1 = 22.46, where S_cone is a perfectly finite 1.07e-15
# and the completion has nothing to repair. The absolute form returns a
# nonzero site there and the pin correctly refuses it. An absolute width
# says "b1c is near 1" only while mu1 is far from its floor; the relative
# width says it always.
TAIL_CONC_MIN = 1e-9  # Minimum concentration (mol/m^3) to actuate handshake

# Gas constant for the radical_qssa_unzip Arrhenius evaluations, J/(mol K).
# Pinned to 8.314 by the M1/M2 channel contract (Ea is stored in J/mol by
# convention, not dimensionally enforced) -- deliberately NOT constants.R
# so the channel's k(T) is bit-reproducible against the documented law.
QSSA_R_GAS = 8.314

# +inf sentinel for the RHS runtime degenerate-rate guard: `x < QSSA_INF` is
# a single comparison that rejects BOTH inf and NaN (NaN fails every
# comparison), keeping the per-RHS-call hot path to a plain compare chain.
QSSA_INF = float("inf")


def _raise_degenerate_qssa_rates(pool_label, T, ki, kdp, kt, ktr=None,
                                 kia=None, ktrec=None, ktdisp=None):
    """Cold path for the radical_qssa_unzip runtime degenerate-rate guard.

    M1 validates config FINITENESS, but Arrhenius EVALUATION at the solver
    temperature can still degenerate: kt(T) can underflow to exactly 0.0
    (R_ss = sqrt(f*ki*B/kt) then divides by zero into inf) and a large A or n
    can overflow A*T**n to inf (NaN downstream). Fail loud, naming the pool,
    the temperature and every offending constant -- never silently zero or
    clamp the rate. Only called off the hot path, when the cheap comparison
    chain in the RHS has already found a degenerate value."""
    problems = []
    for name, val in (("initiation ki", ki), ("depropagation kdp", kdp),
                      ("termination kt", kt), ("transfer ktr", ktr),
                      ("initiation_allyl kia", kia),
                      ("termination_recombination ktrec", ktrec),
                      ("termination_disproportionation ktdisp", ktdisp)):
        if val is None:
            continue
        if not math.isfinite(val):
            problems.append(f"{name}(T)={val!r} is non-finite "
                            f"(Arrhenius overflow at evaluation)")
        elif name == "termination kt" and val <= 0.0:
            problems.append(
                f"{name}(T)={val!r} underflowed to <= 0 "
                f"(R_ss = sqrt(f*ki*B/kt) would divide by zero)")
    raise ValueError(
        f"Pool {pool_label}: radical_qssa_unzip rate evaluation is "
        f"degenerate at T={T:g} K: " + "; ".join(problems) + ". Config-time "
        f"validation (M1) pins finiteness of A/n/Ea but cannot catch "
        f"runtime overflow/underflow of k(T) = A*T**n*exp(-Ea/(R*T)); fix "
        f"the offending Arrhenius block for the reactor temperature range.")

# Pool moment-flux archetypes. Mirror of rmgpy.polymer.PolymerFluxArchetype
# (not imported to avoid a solver->polymer module cycle); equality is pinned
# by test_flux_archetype_constants_match_enum.
FLUX_NONE = 0
FLUX_SAME_POOL = 1
FLUX_MIGRATION = 2
FLUX_SCISSION_FRAGMENT = 3
FLUX_UNRESOLVED = 4
FLUX_DISCRETE_CHIP = 5
FLUX_VOLATILE_EJECTION = 6
FLUX_MOMENT_CREDIT_CONDUIT = 7


# ======================================================================================
# DATA CONFIGURATION
# ======================================================================================

# Pinned basis string for the radical QSSA unzip channel (forward-compat pin):
# the M2 rate law will count initiation sites from the pool's backbone-bond
# count (mu1 - mu0). Any other basis is rejected until a rate law that
# consumes it actually ships, so a deck cannot silently request semantics
# that do not exist yet.
RADICAL_QSSA_UNZIP_BASIS = "backbone_bonds_mu1_minus_mu0"

_RADICAL_QSSA_ARRHENIUS_KEYS = ("A", "n", "Ea")
_RADICAL_QSSA_MANDATORY_BLOCKS = ("initiation", "depropagation", "termination")

# Weak-link allyl/U-state vocabulary (schema-2.2 milestone i: config +
# validation ONLY -- no RHS reads it yet, and the sidecar schema bump is a
# later milestone). ALL-OR-NOTHING as a group, and mutually exclusive with
# the legacy SUMMED 'termination' block: the U-state (unsaturated tail ends)
# is PRODUCED by the disproportionation branch specifically, so a summed kt
# cannot source U -- the split blocks are structurally required.
_RADICAL_QSSA_WEAKLINK_TRIPLET_BLOCKS = (
    "initiation_allyl", "termination_recombination",
    "termination_disproportionation")
_RADICAL_QSSA_WEAKLINK_KEYS = _RADICAL_QSSA_WEAKLINK_TRIPLET_BLOCKS + (
    "unsaturated_tail_ends_initial",)

_RADICAL_QSSA_ALLOWED_KEYS = frozenset(
    _RADICAL_QSSA_MANDATORY_BLOCKS
    + ("transfer", "efficiency", "monomer_yield", "basis")
    + _RADICAL_QSSA_WEAKLINK_KEYS)


def _validate_qssa_arrhenius_triplet(pool_label, block_name, triplet):
    """Validate one radical_qssa_unzip Arrhenius block {A, n, Ea}; return it
    normalized to plain floats.

    Rules (from day one, unlike the legacy k_unzip scalar): every value must
    be a real number and FINITE (NaN/inf rejected explicitly), A > 0,
    Ea >= 0, n any finite float. Raises ValueError naming the pool and the
    offending block/key.
    """
    if not isinstance(triplet, dict):
        raise ValueError(
            f"Pool {pool_label}: radical_qssa_unzip block '{block_name}' must "
            f"be a dict {{A, n, Ea}}, got {type(triplet).__name__}.")
    missing = [k for k in _RADICAL_QSSA_ARRHENIUS_KEYS if k not in triplet]
    extra = sorted(set(triplet) - set(_RADICAL_QSSA_ARRHENIUS_KEYS))
    if missing or extra:
        raise ValueError(
            f"Pool {pool_label}: radical_qssa_unzip block '{block_name}' must "
            f"have exactly the keys {{A, n, Ea}}; missing {missing or 'none'}, "
            f"unknown {extra or 'none'}.")
    out = {}
    for key in _RADICAL_QSSA_ARRHENIUS_KEYS:
        val = triplet[key]
        if isinstance(val, bool) or not isinstance(val, (int, float)):
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip {block_name}.{key}="
                f"{val!r} must be a number.")
        val = float(val)
        if not math.isfinite(val):
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip {block_name}.{key}="
                f"{val!r} is not finite (NaN/inf are rejected).")
        out[key] = val
    if out["A"] <= 0.0:
        raise ValueError(
            f"Pool {pool_label}: radical_qssa_unzip {block_name}.A={out['A']:g} "
            f"must be > 0 (a non-positive pre-exponential is not a valid rate "
            f"constant).")
    if out["Ea"] < 0.0:
        raise ValueError(
            f"Pool {pool_label}: radical_qssa_unzip {block_name}.Ea="
            f"{out['Ea']:g} must be >= 0 [J/mol].")
    return out


def validate_radical_qssa_unzip(pool_label, channel):
    """Validate a radical_qssa_unzip channel config dict and return it
    normalized (defaults filled, all numerics coerced to float).

    This is the single source of truth for the channel's FIELD rules, shared
    by the deck helper (rmgpy/rmg/input.py, re-raised as InputError),
    PolymerPool.to_config (rmgpy/rmg/polymer_input.py) and the solver's own
    validate_configuration (last line of defense: a directly-constructed
    PolymerPoolConfig bypasses both upstream layers). The per-layer CROSS
    invariants (resolvable monomer routing; mutual exclusion with k_unzip > 0)
    stay with the callers, whose field names differ.

    Channel contract (M1; the QSSA rate law lands in M2 -- until then the
    stored config is inert, nothing in the residual reads it):

    - ``initiation``, ``depropagation``, ``termination``: mandatory Arrhenius
      triplets ``{A, n, Ea}``. SI units BY CONVENTION, not dimensionally
      enforced: A [s^-1] for the unimolecular blocks (initiation,
      depropagation), A [m^3 mol^-1 s^-1] for the bimolecular termination
      block; Ea [J/mol]; n dimensionless. All values finite; A > 0; Ea >= 0.
    - ``transfer``: optional Arrhenius triplet (default None), same
      finite/positivity rules. Units A [s^-1]: the M2 rate law is
      PSEUDO-FIRST-ORDER in the active end R (ktr multiplies R directly in
      the balance ktr*R + 2*kt*R^2 = G_R), so a literature bimolecular
      k_tr [L mol^-1 s^-1 / m^3 mol^-1 s^-1] must be premultiplied by the
      relevant substrate concentration [mol/m^3, SI] BEFORE entering this
      config -- do not drop that concentration factor silently.
    - ``efficiency`` (f_i) and ``monomer_yield`` (y_m): floats in (0, 1],
      default 1.0.
    - ``basis``: must equal RADICAL_QSSA_UNZIP_BASIS (default; forward-compat
      pin).

    Weak-link allyl/U-state vocabulary (schema-2.2 milestone i; OPTIONAL as
    a GROUP, all-or-nothing within it):

    - ``initiation_allyl``: Arrhenius triplet (A [s^-1]), same rules as
      ``initiation`` -- the weak-link (allylic) initiation channel driven by
      the pool's unsaturated-tail-ends state U.
    - ``termination_recombination`` and ``termination_disproportionation``:
      Arrhenius triplets (A [m^3 mol^-1 s^-1]), same rules as the legacy
      summed ``termination``.
    - ``unsaturated_tail_ends_initial``: finite float >= 0, the initial U
      amount [mol] -- the SAME amount basis as mu0 (consumers divide by
      V_poly).
    - MUTUAL EXCLUSION: if ANY weak-link key is present, ALL four must be
      present and the legacy summed ``termination`` must be ABSENT. U is
      produced by the disproportionation branch specifically, so a summed
      kt cannot source U. A channel with NONE of the weak-link keys
      normalizes EXACTLY as before this vocabulary existed (legacy freeze).

    Raises ValueError naming the pool on any violation.
    """
    if not isinstance(channel, dict):
        raise ValueError(
            f"Pool {pool_label}: radical_qssa_unzip must be a dict, got "
            f"{type(channel).__name__}.")
    unknown = sorted(set(channel) - _RADICAL_QSSA_ALLOWED_KEYS)
    if unknown:
        raise ValueError(
            f"Pool {pool_label}: radical_qssa_unzip has unknown key(s) "
            f"{unknown}; allowed keys are {sorted(_RADICAL_QSSA_ALLOWED_KEYS)}.")

    weaklink_present = [k for k in _RADICAL_QSSA_WEAKLINK_KEYS
                        if k in channel]
    if weaklink_present:
        if "termination" in channel:
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip carries the legacy "
                f"SUMMED 'termination' block together with weak-link key(s) "
                f"{weaklink_present}. They are mutually exclusive: the "
                f"weak-link U-state channel requires the SPLIT termination "
                f"blocks (termination_recombination + "
                f"termination_disproportionation) because U production is "
                f"sourced by the disproportionation branch specifically -- "
                f"a summed termination cannot source U. Remove 'termination' "
                f"(split it into the two branches) or remove all weak-link "
                f"keys.")
        missing_weaklink = [k for k in _RADICAL_QSSA_WEAKLINK_KEYS
                            if k not in channel]
        if missing_weaklink:
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip weak-link channel is "
                f"incomplete: key(s) {weaklink_present} present but "
                f"{missing_weaklink} missing. The weak-link vocabulary is "
                f"all-or-nothing: initiation_allyl, "
                f"termination_recombination, termination_disproportionation "
                f"and unsaturated_tail_ends_initial must ALL be present "
                f"(and the legacy summed 'termination' absent).")
        mandatory_blocks = (("initiation", "depropagation")
                            + _RADICAL_QSSA_WEAKLINK_TRIPLET_BLOCKS)
    else:
        mandatory_blocks = _RADICAL_QSSA_MANDATORY_BLOCKS

    normalized = {}
    for block in mandatory_blocks:
        if block not in channel:
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip is missing the "
                f"mandatory Arrhenius block '{block}' (required blocks: "
                f"{list(mandatory_blocks)}).")
        normalized[block] = _validate_qssa_arrhenius_triplet(
            pool_label, block, channel[block])

    if weaklink_present:
        u0 = channel["unsaturated_tail_ends_initial"]
        if isinstance(u0, bool) or not isinstance(u0, (int, float)):
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip "
                f"unsaturated_tail_ends_initial={u0!r} must be a number "
                f"[mol; same amount basis as mu0].")
        u0 = float(u0)
        if not math.isfinite(u0):
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip "
                f"unsaturated_tail_ends_initial={u0!r} is not finite "
                f"(NaN/inf are rejected).")
        if u0 < 0.0:
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip "
                f"unsaturated_tail_ends_initial={u0:g} must be >= 0 "
                f"[mol; same amount basis as mu0].")
        normalized["unsaturated_tail_ends_initial"] = u0

    transfer = channel.get("transfer", None)
    normalized["transfer"] = (
        None if transfer is None
        else _validate_qssa_arrhenius_triplet(pool_label, "transfer", transfer))

    for name in ("efficiency", "monomer_yield"):
        val = channel.get(name, 1.0)
        if isinstance(val, bool) or not isinstance(val, (int, float)):
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip {name}={val!r} must "
                f"be a number.")
        val = float(val)
        if not math.isfinite(val) or not (0.0 < val <= 1.0):
            raise ValueError(
                f"Pool {pool_label}: radical_qssa_unzip {name}={val!r} must "
                f"be a finite float in (0, 1].")
        normalized[name] = val

    basis = channel.get("basis", RADICAL_QSSA_UNZIP_BASIS)
    if basis != RADICAL_QSSA_UNZIP_BASIS:
        raise ValueError(
            f"Pool {pool_label}: radical_qssa_unzip basis={basis!r} is not "
            f"supported; the only allowed basis is "
            f"'{RADICAL_QSSA_UNZIP_BASIS}' (forward-compat pin: other bases "
            f"are rejected until a rate law that consumes them ships).")
    normalized["basis"] = basis

    return normalized


# Ratified label convention (radical-homolysis conduit, adjudicated round 66;
# wording per round 67 ruling a) for the two end-radical daughter pools
# spawned by a k_homolysis parent. The contract is POSITIONAL: 'primary' in
# the suffix means the open-*1 end radical, 'secondary' the open-*2 end
# radical of the backbone C-C cut. The names are NOT a claim of
# primary/secondary radical character -- that reading holds only for PP under
# its head-to-tail repeat convention, is orientation-dependent, and is false
# for other polymers.
K_HOMOLYSIS_DAUGHTER_SUFFIXES = ("_rad_primary_end", "_rad_secondary_end")


def _validate_kernel_arrhenius_triplet(pool_label, param_name, triplet):
    """Shared field validator for pool-level kernel Arrhenius triplets
    {A, n, Ea} (k_homolysis, k_depropagation); returns the triplet
    normalized to plain floats. Error text stays byte-compatible with the
    historical k_homolysis messages (parameter name interpolated)."""
    if not isinstance(triplet, dict):
        raise ValueError(
            f"Pool {pool_label}: {param_name} must be a dict {{A, n, Ea}} "
            f"(Arrhenius triplet, SI convention: A [s^-1], Ea [J/mol]), got "
            f"{type(triplet).__name__}. A bare scalar is rejected by design "
            f"(round 66): the kernel must evaluate k(T) at the solver's "
            f"runtime temperature.")
    missing = [k for k in _RADICAL_QSSA_ARRHENIUS_KEYS if k not in triplet]
    extra = sorted(set(triplet) - set(_RADICAL_QSSA_ARRHENIUS_KEYS))
    if missing or extra:
        raise ValueError(
            f"Pool {pool_label}: {param_name} must have exactly the keys "
            f"{{A, n, Ea}}; missing {missing or 'none'}, unknown "
            f"{extra or 'none'}.")
    out = {}
    for key in _RADICAL_QSSA_ARRHENIUS_KEYS:
        val = triplet[key]
        if isinstance(val, bool) or not isinstance(val, (int, float)):
            raise ValueError(
                f"Pool {pool_label}: {param_name} {key}={val!r} must be a "
                f"number.")
        val = float(val)
        if not math.isfinite(val):
            raise ValueError(
                f"Pool {pool_label}: {param_name} {key}={val!r} is not finite "
                f"(NaN/inf are rejected).")
        out[key] = val
    if out["A"] <= 0.0:
        raise ValueError(
            f"Pool {pool_label}: {param_name} A={out['A']:g} must be > 0 (a "
            f"non-positive pre-exponential is not a valid rate constant).")
    if out["Ea"] < 0.0:
        raise ValueError(
            f"Pool {pool_label}: {param_name} Ea={out['Ea']:g} must be >= 0 "
            f"[J/mol].")
    return out


def validate_k_homolysis(pool_label, triplet):
    """Validate a k_homolysis Arrhenius triplet {A, n, Ea} and return it
    normalized to plain floats (single source of truth for the field rules,
    shared by the deck helper rmgpy/rmg/input.py -- re-raised as InputError
    -- PolymerPool.to_config and the solver's own validate_configuration).

    Same rules and SI convention as the radical_qssa_unzip Arrhenius blocks
    (deliberately NOT a bare scalar -- round 66: 'a scalar precomputed at
    1100 K will poison any ramp/TA replay'; the solver evaluates
    k(T) = A*T^n*exp(-Ea/(R_gas*T)) at its runtime T): exactly the keys
    {A, n, Ea}; every value a finite real number; A > 0 [s^-1]; Ea >= 0
    [J/mol]; n dimensionless. Units are convention, not dimensionally
    enforced (matches the QSSA blocks' posture).
    """
    return _validate_kernel_arrhenius_triplet(pool_label, "k_homolysis",
                                              triplet)


def validate_k_depropagation(pool_label, triplet):
    """Validate a k_depropagation Arrhenius triplet {A, n, Ea} and return it
    normalized to plain floats (single source of truth for the field rules,
    shared by the deck helper rmgpy/rmg/input.py -- re-raised as InputError
    -- derive_daughter_pool_configs and the solver's own
    validate_configuration).

    End-radical DEPROPAGATION kernel (adjudicated round 74 SS2, the run-6
    no-outlet wall fix): per-chain-end unzip frequency
    k_dep(T) = A*T^n*exp(-Ea/(R_gas*T)) evaluated at the solver's RUNTIME
    temperature (same not-a-bare-scalar posture as k_homolysis, round 66).
    A [s^-1], Ea [J/mol], n dimensionless; units are convention, not
    dimensionally enforced."""
    return _validate_kernel_arrhenius_triplet(pool_label, "k_depropagation",
                                              triplet)


# Catalog-driven pool kernels (slice S5 a + b of the shared-chemistry design,
# /home/alon/runs/i204/design.md §2.2): their coefficients are a kMC kernel-
# rates table's, per template, at one temperature. Lengths in a term count
# the term's period, `periods_per_unit` of them per monomer unit of the pool
# (the pools keep their monomer unit; the mu3 closure is homogeneous in the
# length scale, so working in periods is exact).
KERNEL_HOMOLYSIS_KEYS = ("template", "w", "phi", "a", "b", "min_fragment",
                         "periods_per_unit", "daughter")
# I-278: a homolysis term whose two fragments land in different pools, the
# chain losing `residue` periods to the cut, carries all three.
KERNEL_HOMOLYSIS_PARTNER_KEYS = ("partner_daughter", "partner_min_fragment", "residue")
KERNEL_RECOMBINATION_KEYS = ("template", "k_pair", "sites_per_chain",
                             "periods_per_unit", "product")
# I-266: a recombination term joining one chain of this pool and one of
# another pool (two end-radical classes; an addition of a radical to an
# alkene pool) carries both.
KERNEL_RECOMBINATION_PARTNER_KEYS = ("partner", "partner_sites_per_chain")
# I-266: one chain of the pool adds one monomer unit, its monomer_product
# consumed.
KERNEL_PROPAGATION_KEYS = ("template", "k_pair")
KERNEL_TRANSFER_KEYS = ("template", "k_pair", "phi", "partner", "product",
                        "partner_product")
# S5 d: a transfer term whose sites grow with a chain carries both.
KERNEL_TRANSFER_SITE_KEYS = ("sites", "periods_per_unit")
KERNEL_MOVE_KEYS = ("template", "k", "a", "b", "periods_per_unit", "product")
# I-273: every chain of the pool breaks at k per chain into fragments, each
# product pool receiving fixed coefficients of the parent's moments.
KERNEL_INTERIOR_SCISSION_KEYS = ("template", "k", "periods_per_unit", "products")
# I-278: the products' lengths sum to n - residue periods per event.
KERNEL_INTERIOR_SCISSION_RESIDUE_KEYS = ("residue",)
# Any term may also declare the shortest chain (periods) its site law counts
# correctly; homolysis' is 2*min_fragment - 1 when undeclared (the cut count
# a*n + b turns negative below it). Reported, never acted on.
KERNEL_VALIDITY_KEY = "valid_from_periods"
CATALOG_KERNELS = ("kernel_homolysis", "kernel_recombination",
                   "kernel_transfer", "kernel_move", "kernel_propagation",
                   "kernel_interior_scission")


def catalog_kernel_connections(pool):
    """Every chain move of `pool`'s (validated) catalog kernel terms, as
    (kernel, template, source pool label, target pool label): homolysis
    pool -> daughter, recombination pool -> product (and partner ->
    product for a join of two pools), transfer pool -> product and partner
    -> partner_product, move pool -> product, interior scission pool ->
    each product. Propagation keeps its chains in the pool."""
    out = [("kernel_interior_scission", term["template"], pool.label, product["pool"])
           for term in getattr(pool, "kernel_interior_scission", None) or ()
           for product in term["products"]]
    for kernel, moves in (("kernel_homolysis", (("", "daughter"), ("", "partner_daughter"))),
                          ("kernel_recombination", (("", "product"),
                                                    ("partner", "product"))),
                          ("kernel_transfer", (("", "product"),
                                               ("partner", "partner_product"))),
                          ("kernel_move", (("", "product"),))):
        for term in getattr(pool, kernel, None) or ():
            for src, dst in moves:
                if (src and src not in term) or dst not in term:
                    continue
                out.append((kernel, term["template"],
                            term[src] if src else pool.label, term[dst]))
    return out


def _range_power_sums(lo, shift):
    """Exact coefficients, in powers of n, of S_j(n) = sum_{k=lo}^{n-shift} k^j
    for j = 1, 2: ([c0, c1, c2], [d0, d1, d2, d3]) with S_1 = sum c_i n^i and
    S_2 = sum d_i n^i, valid wherever n - shift >= lo - 1 (an empty range is 0).
    Faulhaber's sums at the upper end M = n - shift less the same at lo - 1."""
    from fractions import Fraction

    def at(m):          # coefficients in M of sum_{k=1}^{M} k and k^2
        return ([Fraction(0), Fraction(1, 2), Fraction(1, 2)],
                [Fraction(0), Fraction(1, 6), Fraction(1, 2), Fraction(1, 3)])

    def shifted(poly):  # poly(M) with M = n - shift, as a polynomial in n
        out = [Fraction(0)] * len(poly)
        for power, coeff in enumerate(poly):
            for i in range(power + 1):
                out[i] += coeff * math.comb(power, i) * Fraction(-shift) ** (power - i)
        return out

    def at_value(poly, m):
        return sum(c * Fraction(m) ** i for i, c in enumerate(poly))

    s1, s2 = at(0)
    out = []
    for poly in (s1, s2):
        top = shifted(poly)
        top[0] -= at_value(poly, lo - 1)
        out.append([float(c) for c in top])
    return out[0], out[1]


def _kernel_number(pool_label, kernel, template, key, value, *, positive=False,
                   unit_interval=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) \
            or not math.isfinite(value):
        raise ValueError(f"Pool {pool_label}: {kernel} term {template!r}: "
                         f"{key}={value!r} must be a finite number.")
    value = float(value)
    if positive and not value > 0.0:
        raise ValueError(f"Pool {pool_label}: {kernel} term {template!r}: "
                         f"{key}={value!r} must be > 0.")
    if unit_interval and not 0.0 <= value <= 1.0:
        raise ValueError(f"Pool {pool_label}: {kernel} term {template!r}: "
                         f"{key}={value!r} must be in [0, 1].")
    return value


def _kernel_integer(pool_label, kernel, template, key, value):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"Pool {pool_label}: {kernel} term {template!r}: "
                         f"{key}={value!r} must be an integer >= 1.")
    return int(value)


def _kernel_terms(pool_label, kernel, terms, keys, target_key, optional=()):
    """The shared shape rules of a kernel term list; yields (template, term).
    `optional` keys come all together or not at all."""
    if not isinstance(terms, (list, tuple)) or not terms:
        raise ValueError(f"Pool {pool_label}: {kernel} must be a non-empty "
                         f"list of term dicts, got {terms!r}.")
    allowed = [set(keys)] + ([set(keys) | set(optional)] if optional else [])
    seen = set()
    for term in terms:
        if not isinstance(term, dict) or set(term) - {KERNEL_VALIDITY_KEY} not in allowed:
            extra = (f" (or with all of {list(optional)} too)"
                     if optional else "")
            extra += f", and optionally {KERNEL_VALIDITY_KEY!r}"
            raise ValueError(f"Pool {pool_label}: {kernel} term {term!r} must "
                             f"have exactly the keys {list(keys)}{extra}.")
        template = term["template"]
        if not isinstance(template, str) or not template:
            raise ValueError(f"Pool {pool_label}: {kernel} template "
                             f"{template!r} must be a non-empty string.")
        if template in seen:
            raise ValueError(f"Pool {pool_label}: {kernel} has a duplicate "
                             f"template {template!r}.")
        seen.add(template)
        if target_key is not None:
            target = term[target_key]
            if not isinstance(target, str) or not target or target == pool_label:
                raise ValueError(f"Pool {pool_label}: {kernel} term {template!r}: "
                                 f"{target_key}={target!r} must name another pool.")
        yield template, term


def _kernel_validity(pool_label, kernel, template, term, floor=None,
                     floor_name=None):
    """The term's declared validity length {KERNEL_VALIDITY_KEY: n} (an
    integer >= 1, and >= the law's own `floor` where it has one), or {}."""
    if KERNEL_VALIDITY_KEY not in term:
        return {}
    n = _kernel_integer(pool_label, kernel, template, KERNEL_VALIDITY_KEY,
                        term[KERNEL_VALIDITY_KEY])
    if floor is not None and n < floor:
        raise ValueError(
            f"Pool {pool_label}: {kernel} term {template!r}: "
            f"{KERNEL_VALIDITY_KEY}={n} is below {floor_name} = {floor}, "
            f"where the site law already counts negative sites.")
    return {KERNEL_VALIDITY_KEY: n}


def validate_kernel_homolysis(pool_label, terms):
    """Validate a kernel_homolysis term list and return it normalized (single
    source of truth, shared by the deck helper -- re-raised as InputError --
    PolymerPool.to_config and the solver's validate_configuration).

    Each term is one template of a kernel-rates table: a chain of n periods
    has a*n + b sites (a > 0), each cut at a position k with
    min_fragment <= k <= n - min_fragment, so b = a*(1 - 2*min_fragment)
    (refused otherwise: the moment law below is that cut law); the event rate
    per site is w*phi [1/s], w > 0 the table's weight, phi in [0, 1] the
    template's cage-escape factor. Both fragments land in the `daughter`
    pool. The law is exact in moments for chains of a*n + b >= 0 sites; the
    producer refuses tables whose states leave that domain.

    I-278, two daughters: a term with `partner_daughter` (a pool; it may be
    the daughter), `partner_min_fragment` e' (>= 1) and `residue` r (>= 0),
    all three or none, has b = a*(1 - r - e - e'): a cut leaves k periods in
    `daughter` and n - r - k in `partner_daughter`, k uniform over
    [e, n - r - e']."""
    out = []
    for template, term in _kernel_terms(pool_label, "kernel_homolysis", terms,
                                        KERNEL_HOMOLYSIS_KEYS, "daughter",
                                        optional=KERNEL_HOMOLYSIS_PARTNER_KEYS):
        num = lambda key, **kw: _kernel_number(pool_label, "kernel_homolysis",
                                               template, key, term[key], **kw)
        w = num("w", positive=True)
        phi = num("phi", unit_interval=True)
        a = num("a", positive=True)
        b = num("b")
        e = _kernel_integer(pool_label, "kernel_homolysis", template,
                            "min_fragment", term["min_fragment"])
        c = _kernel_integer(pool_label, "kernel_homolysis", template,
                            "periods_per_unit", term["periods_per_unit"])
        two = {}
        floor, floor_name = 2 * e - 1, "2*min_fragment - 1"
        if "partner_daughter" in term:
            partner = term["partner_daughter"]
            if not isinstance(partner, str) or not partner or partner == pool_label:
                raise ValueError(
                    f"Pool {pool_label}: kernel_homolysis term {template!r}: "
                    f"partner_daughter={partner!r} must name another pool "
                    f"(it may equal the daughter).")
            e2 = _kernel_integer(pool_label, "kernel_homolysis", template,
                                 "partner_min_fragment",
                                 term["partner_min_fragment"])
            residue = term["residue"]
            if isinstance(residue, bool) or not isinstance(residue, int) \
                    or residue < 0:
                raise ValueError(
                    f"Pool {pool_label}: kernel_homolysis term {template!r}: "
                    f"residue={residue!r} must be an integer >= 0.")
            want = a * (1 - residue - e - e2)
            if not math.isclose(b, want, rel_tol=1e-12, abs_tol=0.0):
                raise ValueError(
                    f"Pool {pool_label}: kernel_homolysis term {template!r}: "
                    f"b={b!r} is not b = a*(1 - residue - min_fragment - "
                    f"partner_min_fragment) = {want!r}; the sites a*n + b "
                    f"are not the cut positions min_fragment..n - residue - "
                    f"partner_min_fragment the moment law integrates.")
            two = {"partner_daughter": partner, "partner_min_fragment": e2,
                   "residue": int(residue)}
            floor = e + e2 + residue - 1
            floor_name = "min_fragment + partner_min_fragment + residue - 1"
        elif not math.isclose(b, a * (1 - 2 * e), rel_tol=1e-12, abs_tol=0.0):
            raise ValueError(
                f"Pool {pool_label}: kernel_homolysis term {template!r}: "
                f"b={b!r} is not b = a*(1 - 2*min_fragment) = "
                f"{a * (1 - 2 * e)!r}; the sites a*n + b are not the cut "
                f"positions min_fragment..n - min_fragment the moment law "
                f"integrates.")
        out.append({"template": template, "w": w, "phi": phi, "a": a, "b": b,
                    "min_fragment": e, "periods_per_unit": c,
                    "daughter": term["daughter"], **two,
                    **_kernel_validity(pool_label, "kernel_homolysis", template,
                                       term, floor, floor_name)})
    return out


def validate_kernel_recombination(pool_label, terms):
    """Validate a kernel_recombination term list and return it normalized
    (shared like validate_kernel_homolysis).

    Each term is one template joining two chains of this pool, each with
    `sites_per_chain` sites (> 0, the same on every chain), into one chain of
    n0 + n1 periods in the `product` pool at k_pair*S**2 per volume, S the
    pool's site concentration [mol/m^3] and k_pair [m^3/(mol s)] > 0 the
    rate coefficient per ordered site pair in MoM's concentration
    convention (n / V_poly).

    I-266, a join of two pools: a term with `partner` (another pool) and
    `partner_sites_per_chain` (> 0), both or neither, joins one chain of
    this pool and one of the partner into one chain of n + m periods in
    `product` at k_pair*S*S' per volume, S' the partner's site
    concentration."""
    out = []
    for template, term in _kernel_terms(pool_label, "kernel_recombination",
                                        terms, KERNEL_RECOMBINATION_KEYS,
                                        "product",
                                        optional=KERNEL_RECOMBINATION_PARTNER_KEYS):
        num = lambda key, **kw: _kernel_number(pool_label,
                                               "kernel_recombination",
                                               template, key, term[key], **kw)
        partner = {}
        if "partner" in term:
            if not isinstance(term["partner"], str) or not term["partner"] \
                    or term["partner"] == pool_label:
                raise ValueError(
                    f"Pool {pool_label}: kernel_recombination term "
                    f"{template!r}: partner={term['partner']!r} must name "
                    f"another pool (a join within this pool has no partner).")
            partner = {"partner": term["partner"],
                       "partner_sites_per_chain": num("partner_sites_per_chain",
                                                      positive=True)}
        out.append({"template": template, "k_pair": num("k_pair", positive=True),
                    "sites_per_chain": num("sites_per_chain", positive=True),
                    "periods_per_unit": _kernel_integer(
                        pool_label, "kernel_recombination", template,
                        "periods_per_unit", term["periods_per_unit"]),
                    "product": term["product"], **partner,
                    **_kernel_validity(pool_label, "kernel_recombination",
                                       template, term)})
    return out


def validate_kernel_propagation(pool_label, terms):
    """Validate a kernel_propagation term list and return it normalized
    (shared like validate_kernel_homolysis).

    Each term (I-266) is one kernel of a kernel-rates table whose events add
    one monomer unit to a chain of this pool, consuming one of the pool's
    monomer_product, at k_pair*c_M per chain, k_pair [m^3/(mol s)] > 0 (the
    sites folded in) and c_M the monomer's amount over the polymer volume
    (kMC's convention: every small species at its melt concentration)."""
    out = []
    for template, term in _kernel_terms(pool_label, "kernel_propagation",
                                        terms, KERNEL_PROPAGATION_KEYS, None):
        out.append({"template": template,
                    "k_pair": _kernel_number(pool_label, "kernel_propagation",
                                             template, "k_pair",
                                             term["k_pair"], positive=True),
                    **_kernel_validity(pool_label, "kernel_propagation",
                                       template, term)})
    return out


def validate_kernel_transfer(pool_label, terms):
    """Validate a kernel_transfer term list and return it normalized (shared
    like validate_kernel_homolysis).

    Each term (S5 c + e) is one kernel of a kernel-rates table whose events
    take one chain of this pool and one chain of the `partner` pool (this
    pool too, e.g. end + end disproportionation) and move them, lengths
    unchanged, into the `product` and `partner_product` pools at
    k_pair*phi*c*c_partner per volume, c the pools' chain concentrations
    [mol/m^3], k_pair [m^3/(mol s)] > 0 the coefficient per ordered chain
    pair in MoM's concentration convention (the sites folded in: the
    producer maps a kernel only where they are the same on every chain),
    phi in [0, 1] the cage-escape factor of a radical-creating kernel.

    S5 d: a term may also carry `sites` [c00, c10, c01, c11] (finite, not all
    zero) and `periods_per_unit`, both or neither: the rate per ordered pair
    of a chain of n periods here and one of m periods in the partner pool is
    then k_pair*phi*(c00 + c10*n + c01*m + c11*n*m) (sites that grow with a
    chain, e.g. a dead chain's interior H; the chain is picked in proportion
    to them). Without them the rate is k_pair*phi per pair, as before."""
    out = []
    for template, term in _kernel_terms(pool_label, "kernel_transfer", terms,
                                        KERNEL_TRANSFER_KEYS, "product",
                                        optional=KERNEL_TRANSFER_SITE_KEYS):
        num = lambda key, **kw: _kernel_number(pool_label, "kernel_transfer",
                                               template, key, term[key], **kw)
        k_pair = num("k_pair", positive=True)
        phi = num("phi", unit_interval=True)
        site_law = {}
        if "sites" in term:
            sites = term["sites"]
            if not isinstance(sites, (list, tuple)) or len(sites) != 4 or any(
                    isinstance(c, bool) or not isinstance(c, (int, float))
                    or not math.isfinite(c) for c in sites) \
                    or not any(sites):
                raise ValueError(f"Pool {pool_label}: kernel_transfer term "
                                 f"{template!r}: sites={sites!r} must be four "
                                 f"finite numbers [c00, c10, c01, c11], not "
                                 f"all zero.")
            site_law = {"sites": [float(c) for c in sites],
                        "periods_per_unit": _kernel_integer(
                            pool_label, "kernel_transfer", template,
                            "periods_per_unit", term["periods_per_unit"])}
        partner, partner_product = term["partner"], term["partner_product"]
        if not isinstance(partner, str) or not partner:
            raise ValueError(f"Pool {pool_label}: kernel_transfer term "
                             f"{template!r}: partner={partner!r} must name a "
                             f"pool (this one included).")
        if not isinstance(partner_product, str) or not partner_product \
                or partner_product == partner:
            raise ValueError(f"Pool {pool_label}: kernel_transfer term "
                             f"{template!r}: partner_product="
                             f"{partner_product!r} must name a pool other "
                             f"than the partner {partner!r}.")
        out.append({"template": template, "k_pair": k_pair, "phi": phi,
                    "partner": partner, "product": term["product"],
                    "partner_product": partner_product, **site_law,
                    **_kernel_validity(pool_label, "kernel_transfer", template,
                                       term)})
    return out


def validate_kernel_move(pool_label, terms):
    """Validate a kernel_move term list and return it normalized (shared like
    validate_kernel_homolysis).

    Each term (S5 d) is one kernel of a kernel-rates table whose events take
    one chain of this pool into the `product` pool, length unchanged
    (backbiting: an end radical becomes a mid-chain radical, and its
    reverse), at k*(a*n + b) per chain [1/s], n in periods
    (`periods_per_unit` per monomer unit), k > 0, a and b finite, not both
    zero."""
    out = []
    for template, term in _kernel_terms(pool_label, "kernel_move", terms,
                                        KERNEL_MOVE_KEYS, "product"):
        num = lambda key, **kw: _kernel_number(pool_label, "kernel_move",
                                               template, key, term[key], **kw)
        k, a, b = num("k", positive=True), num("a"), num("b")
        if a == 0.0 and b == 0.0:
            raise ValueError(f"Pool {pool_label}: kernel_move term "
                             f"{template!r}: a and b are both zero (no site "
                             f"on any chain).")
        out.append({"template": template, "k": k, "a": a, "b": b,
                    "periods_per_unit": _kernel_integer(
                        pool_label, "kernel_move", template,
                        "periods_per_unit", term["periods_per_unit"]),
                    "product": term["product"],
                    **_kernel_validity(pool_label, "kernel_move", template,
                                       term)})
    return out


def validate_kernel_interior_scission(pool_label, terms):
    """Validate a kernel_interior_scission term list and return it normalized
    (shared like validate_kernel_homolysis).

    Each term (I-273) is one kernel of a kernel-rates table on one position
    class of mid-chain radicals: every chain of this pool breaks at k [1/s]
    per chain (k > 0) into fragments, and each product {pool, moments}
    receives per event the fragments whose lengths L (periods,
    `periods_per_unit` per monomer unit) have E[sum L^j | parent of n
    periods] = sum_{i<=j} q_ji n^i, `moments` the lower-triangular
    [[q00], [q10, q11], [q20, q21, q22]] (finite, q00 >= 0). The products
    must conserve the parent's periods: sum q10 = 0 and sum q11 = 1."""
    out = []
    kernel = "kernel_interior_scission"
    for template, term in _kernel_terms(pool_label, kernel, terms,
                                        KERNEL_INTERIOR_SCISSION_KEYS, None,
                                        optional=KERNEL_INTERIOR_SCISSION_RESIDUE_KEYS):
        residue = term.get("residue", 0)
        if isinstance(residue, bool) or not isinstance(residue, int) or residue < 0:
            raise ValueError(f"Pool {pool_label}: {kernel} term {template!r}: "
                             f"residue={residue!r} must be an integer >= 0.")
        products = term["products"]
        if not isinstance(products, (list, tuple)) or not products:
            raise ValueError(f"Pool {pool_label}: {kernel} term {template!r}: "
                             f"products={products!r} must be a non-empty list "
                             f"of {{pool, moments}}.")
        norm = []
        for product in products:
            if not isinstance(product, dict) or set(product) != {"pool", "moments"}:
                raise ValueError(f"Pool {pool_label}: {kernel} term {template!r}: "
                                 f"product {product!r} must have exactly the "
                                 f"keys ['pool', 'moments'].")
            target = product["pool"]
            if not isinstance(target, str) or not target or target == pool_label:
                raise ValueError(f"Pool {pool_label}: {kernel} term {template!r}: "
                                 f"pool={target!r} must name another pool.")
            q = product["moments"]
            if not isinstance(q, (list, tuple)) or len(q) != 3 or any(
                    not isinstance(row, (list, tuple)) or len(row) != j + 1
                    for j, row in enumerate(q)):
                raise ValueError(f"Pool {pool_label}: {kernel} term {template!r}: "
                                 f"moments={q!r} must be [[q00], [q10, q11], "
                                 f"[q20, q21, q22]].")
            q = [[_kernel_number(pool_label, kernel, template, f"moments[{j}][{i}]", v)
                  for i, v in enumerate(row)] for j, row in enumerate(q)]
            if q[0][0] < 0.0:
                raise ValueError(f"Pool {pool_label}: {kernel} term {template!r}: "
                                 f"q00={q[0][0]!r} fragments per event must be >= 0.")
            norm.append({"pool": target, "moments": q})
        q10 = math.fsum(p["moments"][1][0] for p in norm)
        q11 = math.fsum(p["moments"][1][1] for p in norm)
        scale = max([1.0] + [abs(p["moments"][1][0]) for p in norm])
        if abs(q10 + residue) > 1e-12 * scale or abs(q11 - 1.0) > 1e-12:
            raise ValueError(f"Pool {pool_label}: {kernel} term {template!r}: the "
                             f"products' lengths sum to {q11!r} n + {q10!r} periods "
                             f"per event, not the parent's n - {residue} (residue "
                             f"{residue}): they must conserve its periods but "
                             f"the residue.")
        out.append({"template": template,
                    "k": _kernel_number(pool_label, kernel, template, "k", term["k"],
                                        positive=True),
                    "periods_per_unit": _kernel_integer(
                        pool_label, kernel, template, "periods_per_unit",
                        term["periods_per_unit"]),
                    "products": norm,
                    **({"residue": int(residue)} if "residue" in term else {}),
                    **_kernel_validity(pool_label, kernel, template, term)})
    return out


# Side-group homolysis initiation kernel (FR1-K1, adjudicated adversarial
# round 70): pool-level side-group X-loss (e.g. aliphatic C-Br -> Br.(gas) +
# mid-chain backbone radical; the chain length is UNCHANGED -- this is NOT a
# chain cut, so it is deliberately a distinct kernel from k_homolysis, and
# the S2 H-loss conduit is NOT reused: the ruling requires its own machinery
# and fingerprint class). SGH kernel-v2 spawns NO feature pool -- Br leaves
# through the auxiliary Z-inventory depletion, not a moment-debited daughter.
# The v1 X-loss-feature-pool label convention '{parent}_sidegrp_{sanitized
# channel label}' is RETAINED only for diagnostic tagging / legacy-artifact
# detection (never spawned in v2).
SIDE_GROUP_DAUGHTER_INFIX = "_sidegrp_"

_SIDE_GROUP_CHANNEL_KEYS = ("label", "A", "n", "Ea", "site_selector",
                            "sites_per_unit", "gas_product")

# Structural site-selector vocabulary (round-72 P1 adjudication): every
# side_group_homolysis channel MUST name the environment of the atom its
# leaving group hangs off -- the kinetics label alone must never pick the
# site. Implemented as explicit in-repo atom-environment predicates over
# RMG's own molecular graph (adjudication-equivalent alternative to
# site_smarts: full SMARTS would require an RMG<->RDKit atom-index
# round-trip whose ordering is not contract-stable -- a misordered mapping
# would silently select the wrong atom, the exact bug class being fixed).
# The classes PARTITION the side-group sites (mutually exclusive, jointly
# exhaustive): 'aryl', 'benzylic' and 'aliphatic' for a group on a carbon,
# 'hetero' for a group on any other atom (I-257: an ester's methyl leaves
# from its oxygen). Two channels can only resolve to the SAME sites by
# duplicating a (selector, leaving group) pair -- which the validator
# hard-fails as a double-carry.
SIDE_GROUP_SITE_SELECTORS = ("aryl", "benzylic", "aliphatic", "hetero")

_leaving_group_species = {}


def _leaving_group(leaving):
    """The leaving-group radical as a Species with its resonance structures
    (cached): ``leaving`` is a gas_product SMILES, a Molecule, or the element
    symbol of a monoatomic X (the round-72 API, e.g. 'Br')."""
    from rmgpy.molecule import Molecule
    from rmgpy.molecule.element import get_element
    from rmgpy.species import Species
    if not isinstance(leaving, str):
        spc = Species(molecule=[leaving.copy(deep=True)])
        spc.generate_resonance_structures()
        return spc
    spc = _leaving_group_species.get(leaving)
    if spc is None:
        try:
            get_element(leaving)
            smiles = f"[{leaving}]"
        except Exception:
            smiles = leaving
        spc = Species(molecule=[Molecule().from_smiles(smiles)])
        spc.generate_resonance_structures()
        _leaving_group_species[leaving] = spc
    return spc


def _far_side(monomer, attach, root):
    """Indices of the atoms reached from ``root`` without crossing its bond
    to ``attach``."""
    seen = {root}
    stack = [root]
    while stack:
        atom = stack.pop()
        for nb in atom.bonds:
            if nb is attach and atom is root:
                continue
            if nb not in seen:
                seen.add(nb)
                stack.append(nb)
    return frozenset(monomer.atoms.index(a) for a in seen)


def _side_fragment(monomer, side, root_idx):
    """The atoms ``side`` of the repeat unit cut off as a radical at the
    root."""
    frag = monomer.copy(deep=True)
    keep = [frag.atoms[i] for i in sorted(side)]
    root = frag.atoms[root_idx]
    for atom in [a for a in frag.atoms if a not in keep]:
        frag.remove_atom(atom)
    root.increment_radical()
    frag.clear_labeled_atoms()
    frag.update()
    # a graph edited in place can carry stale connectivity values that
    # false-negative VF2; its adjacency-list round trip does not
    from rmgpy.molecule import Molecule
    return Molecule().from_adjacency_list(frag.to_adjacency_list())


def _same_radical(group, fragment):
    """True if the fragment is the leaving group in any resonance (or
    Kekule) drawing of either: a ring cut off a unit keeps the unit's bond
    representation, the SMILES parse its own."""
    from rmgpy.species import Species
    frag = Species(molecule=[fragment])
    frag.generate_resonance_structures(keep_isomorphic=False)
    return any(g.is_isomorphic(f) for g in group.molecule for f in frag.molecule)


def side_group_sites(monomer, leaving, site_selector):
    """Deterministic structural site selection for the side_group_homolysis
    kernel (round-72 P1; I-257 for polyatomic leaving groups): one
    ``(attach index, root index, side atom indices)`` per side group of the
    repeat unit that is the leaving group and hangs off an attach atom of
    class ``site_selector``, in ``monomer.atoms`` order of the root (the
    canonical deterministic order; the producer removes the FIRST match).

    A side group is the far side of a single, acyclic bond attach-root of
    the unit: the atoms reached from the root without crossing that bond,
    none of them a stitch-labeled atom (so no backbone atom), the root
    closed-shell. It is the leaving group when that side, with one radical
    electron on the root, is isomorphic to the leaving-group radical. For a
    monoatomic X the side is X alone, which is the round-72 removable-X
    predicate (terminal, single-bonded, closed-shell, unlabeled).

    The selector classifies the attach atom, the group's own atoms aside:

    - 'aryl':      a carbon in an aromatic ring (Molecule.get_aromatic_rings
                   -- robust to kekulized vs aromatic bond representations);
    - 'benzylic':  a carbon NOT aromatic but bonded to an aromatic-ring
                   carbon outside the group;
    - 'aliphatic': a carbon neither aromatic nor benzylic;
    - 'hetero':    any other element.
    """
    group = _leaving_group(leaving)
    n_group = len(group.molecule[0].atoms)
    aromatic_atoms = set()
    for ring in monomer.get_aromatic_rings()[0]:
        aromatic_atoms.update(ring)
    matches = []
    for root_idx, root in enumerate(monomer.atoms):
        if root.label or root.radical_electrons != 0:
            continue
        for attach, bond in root.bonds.items():
            if not bond.is_single() or monomer.is_bond_in_cycle(bond):
                continue
            side = _far_side(monomer, attach, root)
            if len(side) != n_group or any(monomer.atoms[i].label for i in side):
                continue
            if not _same_radical(group, _side_fragment(monomer, side, root_idx)):
                continue
            side_atoms = {monomer.atoms[i] for i in side}
            if not attach.is_carbon():
                cls = "hetero"
            elif attach in aromatic_atoms:
                cls = "aryl"
            elif any(nb in aromatic_atoms and nb.is_carbon() and nb not in side_atoms
                     for nb in attach.bonds.keys()):
                cls = "benzylic"
            else:
                cls = "aliphatic"
            if cls == site_selector:
                matches.append((monomer.atoms.index(attach), root_idx, side))
    return tuple(matches)


def side_group_site_atom_indices(monomer, leaving, site_selector):
    """The root atom indices of :func:`side_group_sites` (the X atoms
    themselves for a monoatomic X)."""
    return tuple(root for _attach, root, _side in
                 side_group_sites(monomer, leaving, site_selector))


def sanitize_side_group_channel_label(channel_label):
    """Deterministic sanitizer for a side_group_homolysis channel label used
    inside pool labels and fingerprints: every character outside [A-Za-z0-9_]
    becomes '_' (e.g. 'aliphatic_C-Br' -> 'aliphatic_C_Br'). Single source of
    truth shared by the validator (sanitized-collision duplicate check), the
    producer (rmgpy.polymer.Polymer.generate_side_loss_daughters) and the
    flattener's destination-pool resolution."""
    return "".join(c if (c.isalnum() or c == "_") else "_"
                   for c in channel_label)


def side_group_daughter_pool_label(pool_label, channel_label):
    """The ratified v1 X-loss-feature-pool label convention for one
    (parent, channel) pair: '{parent}_sidegrp_{sanitized channel label}'.
    SGH kernel-v2 spawns no feature pool; this label is retained only for
    diagnostic tagging / legacy-artifact detection."""
    return (f"{pool_label}{SIDE_GROUP_DAUGHTER_INFIX}"
            f"{sanitize_side_group_channel_label(channel_label)}")


def _side_group_gas_mw_g_mol(gas_product):
    """Molar mass [g/mol] of a channel's gas_product SMILES (M_X of the
    ejected X radical). Lazy Molecule import (cycle-avoidance idiom)."""
    from rmgpy.molecule import Molecule
    return Molecule().from_smiles(gas_product).get_molecular_weight() * 1000.0


# How a carrier books its side-group losses (I-257). The inventory
# accounting seeds Z(0) = sites*mu1(0) once and needs the carrier's mu1
# stationary; beside a kernel that moves the carrier's chains it cannot hold,
# so such a carrier moves each reacting chain, intact, into the channel's
# X-loss feature daughter, which carries the exact per-chain defect M_X.
SIDE_GROUP_INVENTORY = "inventory"
SIDE_GROUP_DAUGHTER = "daughter"
# The X-loss daughter itself (I-260): it keeps its chains and counts their
# further losses on an inventory of the channel's intact sites, which its
# parent's moves and its own releases also move (side_loss_parent).
SIDE_GROUP_IN_POOL = "in_pool"
SIDE_GROUP_CHAIN_MOVERS = ("radical_qssa_unzip", "k_unzip", "k_homolysis",
                           "k_depropagation")


def side_group_accounting(pool):
    """SIDE_GROUP_IN_POOL for an X-loss daughter (``side_loss_parent`` set),
    else SIDE_GROUP_DAUGHTER when a kernel of ``pool`` (a PolymerPoolConfig
    or a Polymer: both carry the kernel fields) moves its chains
    (SIDE_GROUP_CHAIN_MOVERS), else SIDE_GROUP_INVENTORY."""
    if getattr(pool, "side_loss_parent", None):
        return SIDE_GROUP_IN_POOL
    for name in SIDE_GROUP_CHAIN_MOVERS:
        value = getattr(pool, name, None)
        if name == "k_unzip":
            if float(value or 0.0) > 0.0:
                return SIDE_GROUP_DAUGHTER
        elif value is not None:
            return SIDE_GROUP_DAUGHTER
    return SIDE_GROUP_INVENTORY


def validate_side_group_homolysis(pool_label, channels, monomer=None):
    """Validate a side_group_homolysis channel LIST and return it normalized
    (single source of truth for the field rules, shared by the deck helper
    rmgpy/rmg/input.py -- re-raised as InputError -- PolymerPool.to_config,
    the producer and the solver's validate_configuration; mirrors
    validate_k_homolysis).

    Each channel is a dict with EXACTLY the keys
    {label, A, n, Ea, site_selector, sites_per_unit, gas_product}:
    - label: non-empty str naming the bond class (e.g. 'aliphatic_C-Br').
      Channels stay SEPARATE, never lumped (round-70 ruling); duplicate
      labels -- raw or after sanitization -- are rejected.
    - A > 0 [s^-1 per site], n dimensionless, Ea >= 0 [J/mol]: same SI
      convention and finite-rejection posture as validate_k_homolysis (the
      solver evaluates k(T) = A*T^n*exp(-Ea/(R_gas*T)) at its runtime T).
    - site_selector: REQUIRED structural site selector (round-72 P1), one
      of SIDE_GROUP_SITE_SELECTORS ('aryl' | 'benzylic' | 'aliphatic' |
      'hetero') -- names the environment of the atom the channel's leaving
      group hangs off. The kinetics label alone must NEVER pick the site:
      on a mixed-site monomer a channel labeled 'aliphatic_C-Br' would
      otherwise remove the first X in atom order (an ARYL one), minting the
      wrong -- or the same -- defect structure under a different rate.
    - sites_per_unit > 0: side-group sites per monomer repeat unit
      (e.g. aryl C-Br = 4, benzylic C-Br = 2). CHECKED, not trusted, when
      ``monomer`` is given (below).
    - gas_product: SMILES of the ejected leaving-group radical, a NEUTRAL
      MONO-RADICAL of any size (I-257: '[Br]', the ester radical 'CO[C]=O',
      phenyl); the group is read from the repeat unit's structure
      (side_group_sites).

    STRUCTURAL layer (round-72 P1): when ``monomer`` -- the pool's parsed
    repeat-unit Molecule -- is given, each channel's selector is resolved
    against it via side_group_site_atom_indices and the law is enforced:
    (a) the selector must match >= 1 side group that is the channel's
    leaving group (no match -> hard fail naming the channel);
    (b) sites_per_unit must EQUAL the selector's structural match count
    (mismatch -> hard fail with both numbers); (c) no two channels may
    resolve to the SAME atom set (that is a double-carry the distinct
    labels were hiding -> hard fail naming both channels). The layers that
    hold the monomer structure (deck helper, PolymerPool.to_config, the
    producer) pass it; the solver's directly-constructed PolymerPoolConfig
    backstop carries no structure and validates shape + vocabulary only.
    """
    if isinstance(channels, dict) or not isinstance(channels, (list, tuple)):
        raise ValueError(
            f"Pool {pool_label}: side_group_homolysis must be a list of "
            f"channel dicts (one entry per bond class, e.g. "
            f"[{{'label': 'aliphatic_C-Br', 'A': ..., 'n': ..., 'Ea': ..., "
            f"'sites_per_unit': ..., 'gas_product': '[Br]'}}]), got "
            f"{type(channels).__name__}.")
    if len(channels) == 0:
        raise ValueError(
            f"Pool {pool_label}: side_group_homolysis must not be an empty "
            f"list -- omit the parameter entirely to disable the kernel "
            f"(an empty list would be a silently inert channel).")
    normalized = []
    seen_sanitized = {}
    for pos, ch in enumerate(channels):
        if not isinstance(ch, dict):
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis entry {pos} must "
                f"be a channel dict, got {type(ch).__name__}.")
        missing = [k for k in _SIDE_GROUP_CHANNEL_KEYS if k not in ch]
        extra = sorted(set(ch) - set(_SIDE_GROUP_CHANNEL_KEYS))
        if missing or extra:
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis entry {pos} "
                f"(label={ch.get('label', '<unset>')!r}) must have exactly "
                f"the keys {{label, A, n, Ea, site_selector, sites_per_unit, "
                f"gas_product}}; missing {missing or 'none'}, unknown "
                f"{extra or 'none'}.")
        label = ch["label"]
        if not isinstance(label, str) or not label.strip():
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis entry {pos} label "
                f"must be a non-empty string naming the bond class, got "
                f"{label!r}.")
        out = {"label": label}
        sel = ch["site_selector"]
        if not isinstance(sel, str) or sel not in SIDE_GROUP_SITE_SELECTORS:
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis channel '{label}' "
                f"site_selector={sel!r} must be one of "
                f"{SIDE_GROUP_SITE_SELECTORS} -- the REQUIRED structural "
                f"site selector pins WHICH X atom the channel removes "
                f"(round-72: a kinetics label alone must never pick the "
                f"site).")
        out["site_selector"] = sel
        for key in ("A", "n", "Ea", "sites_per_unit"):
            val = ch[key]
            if isinstance(val, bool) or not isinstance(val, (int, float)):
                raise ValueError(
                    f"Pool {pool_label}: side_group_homolysis channel "
                    f"'{label}' {key}={val!r} must be a number.")
            val = float(val)
            if not math.isfinite(val):
                raise ValueError(
                    f"Pool {pool_label}: side_group_homolysis channel "
                    f"'{label}' {key}={val!r} is not finite (NaN/inf are "
                    f"rejected).")
            out[key] = val
        if out["A"] <= 0.0:
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis channel '{label}' "
                f"A={out['A']:g} must be > 0 (a non-positive pre-exponential "
                f"is not a valid rate constant).")
        if out["Ea"] < 0.0:
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis channel '{label}' "
                f"Ea={out['Ea']:g} must be >= 0 [J/mol].")
        if out["sites_per_unit"] <= 0.0:
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis channel '{label}' "
                f"sites_per_unit={out['sites_per_unit']:g} must be > 0 "
                f"(side-group X sites per monomer repeat unit).")
        gp = ch["gas_product"]
        if not isinstance(gp, str) or not gp.strip():
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis channel '{label}' "
                f"gas_product must be a SMILES string, got {gp!r}.")
        from rmgpy.molecule import Molecule
        try:
            gmol = Molecule().from_smiles(gp)
        except Exception as e:
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis channel '{label}' "
                f"gas_product={gp!r} does not parse as SMILES ({e}).")
        if gmol.get_radical_count() != 1:
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis channel '{label}' "
                f"gas_product={gp!r} must be a mono-radical (the homolysis "
                f"co-product is the ejected X radical); got radical count "
                f"{gmol.get_radical_count()}.")
        if gmol.get_net_charge() != 0:
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis channel '{label}' "
                f"gas_product={gp!r} carries net charge "
                f"{gmol.get_net_charge()}; a homolysis leaves a neutral "
                f"radical.")
        out["gas_product"] = gp
        san = sanitize_side_group_channel_label(label)
        if san in seen_sanitized:
            raise ValueError(
                f"Pool {pool_label}: side_group_homolysis has duplicate "
                f"channel labels: {seen_sanitized[san]!r} and {label!r} "
                f"collide (sanitized '{san}'). Two decks of the same channel "
                f"would double-carry the same bond class; channels must be "
                f"unique.")
        seen_sanitized[san] = label
        normalized.append(out)

    # STRUCTURAL layer (round-72 P1) -- only when the caller holds the
    # parsed repeat-unit structure; see the docstring for the three laws.
    if monomer is not None:
        from rmgpy.molecule import Molecule
        seen_sites = {}
        for out in normalized:
            label = out["label"]
            gas = Molecule().from_smiles(out["gas_product"])
            sym = (gas.atoms[0].symbol if len(gas.atoms) == 1
                   else out["gas_product"])
            idxs = side_group_site_atom_indices(monomer, gas,
                                                out["site_selector"])
            if not idxs:
                raise ValueError(
                    f"Pool {pool_label}: side_group_homolysis channel "
                    f"'{label}' site_selector='{out['site_selector']}' "
                    f"matches NO removable side-group {sym} (the far side "
                    f"of a single acyclic bond, holding no stitch atom, "
                    f"that is the leaving group with its radical at the "
                    f"cut, on a matching attach atom) in the repeat unit "
                    f"-- the selector must structurally pin which {sym} "
                    f"the channel removes; a non-matching selector would "
                    f"mint the wrong defect (round-72 P1).")
            if float(len(idxs)) != out["sites_per_unit"]:
                raise ValueError(
                    f"Pool {pool_label}: side_group_homolysis channel "
                    f"'{label}' sites_per_unit={out['sites_per_unit']:g} "
                    f"contradicts the structural match count: "
                    f"site_selector='{out['site_selector']}' matches "
                    f"{len(idxs)} {sym} site(s) in the repeat unit. "
                    f"sites_per_unit is CHECKED against the monomer, never "
                    f"trusted (round-72 P1); set it to {len(idxs)} or fix "
                    f"the selector.")
            key = frozenset(idxs)
            if key in seen_sites:
                raise ValueError(
                    f"Pool {pool_label}: side_group_homolysis channels "
                    f"'{seen_sites[key]}' and '{label}' resolve to the SAME "
                    f"{sym} atom set {sorted(idxs)} in the repeat unit -- "
                    f"two rate channels on one structural site double-carry "
                    f"the loss (the distinct labels were hiding it; "
                    f"round-72 P1). Merge them or pick disjoint selectors.")
            seen_sites[key] = label
    return normalized


@dataclass(frozen=True)
class PolymerPoolConfig:
    """
    Configuration for one polymer pool.
    Indices must refer to valid CORE species indices.
    """
    label: str
    xs: int
    explicit_dp_to_species_index: Dict[int, int]
    mu_indices: Tuple[int, int, int]
    monomer_poly_index: Optional[int] = None

    # Monomer (repeat-unit) molecular weight [g/mol]. Consumed by
    # spawn_gate_flux_snapshot (E[n]*MW event-mass calibration, spec
    # 2026-06-10-mass-flux-spawn-gate-design.md §3). Default 0.0 zeroes the
    # pool's snapshot contribution (honest degradation: the spawn gate
    # defers; nothing is fabricated).
    monomer_mw_g_mol: float = 0.0

    # Monomer (repeat-unit) heavy-atom (non-H) count. The SECOND axis of the
    # r89 dual-axis polymer-sized melt gate in _reference_state_tripwire
    # (mirror of rmgpy.polymer._discrete_is_polymer_sized: mass AND heavy
    # structure, both computable, both at/above
    # _IMPOSTOR_DISCRETE_MONOMER_UNITS monomer-equivalents against >= 1
    # pool). Populated by every config producer from the monomer STRUCTURE
    # (rmgpy/rmg/polymer_input.py to_pool_config /
    # derive_daughter_pool_configs; consumer world reconstructs it from the
    # artifact's monomer_adj_list) -- per the r89 data-flow constraint the
    # heavy axis is never approximated from mass. Default 0 = axis
    # UNCOMPUTABLE (honest degradation: tag-branch candidates classify
    # conservative-gas and the tripwire emits an axis-undecidable census
    # warning; nothing is classified blind).
    monomer_heavy_atoms: int = 0

    # Kinetics Parameters
    k_scission: float = 0.0

    # Recession Rate [1/s]
    # REQUIRED for the Hybrid Handshake (Tail -> Explicit flux).
    # This parameter sets the timescale for both the default recession chemistry
    # (n -> n-1) AND the physical handshake flux across the boundary xs.
    k_unzip: float = 0.0

    # Radical QSSA unzip channel (M1: config + validation ONLY -- nothing in
    # the residual reads this yet; the QSSA rate law lands in M2). Normalized
    # dict per validate_radical_qssa_unzip: {initiation, depropagation,
    # termination: {A, n, Ea}, transfer: {A, n, Ea} | None, efficiency,
    # monomer_yield: float in (0, 1], basis: RADICAL_QSSA_UNZIP_BASIS}.
    # Units BY CONVENTION (documented, not dimensionally enforced):
    # initiation A [s^-1], depropagation A [s^-1], termination A
    # [m^3 mol^-1 s^-1] (the only bimolecular block), transfer A [s^-1]
    # (PSEUDO-first-order: ktr multiplies R directly in the M2 rate law --
    # premultiply a literature bimolecular k_tr by the substrate
    # concentration [mol/m^3] before configuring it); Ea [J/mol].
    # Mutually exclusive with k_unzip > 0
    # (double-count guard) and requires monomer_poly_index (the channel
    # reuses the pool's existing monomer routing -- no new routing field).
    radical_qssa_unzip: Optional[Dict[str, object]] = None

    # Radical-homolysis initiation kernel (Stage 1, adjudicated round 66).
    # Normalized Arrhenius triplet {A, n, Ea} per validate_k_homolysis
    # (SI convention: A [s^-1], Ea [J/mol]; NOT a bare scalar -- the RHS
    # evaluates k(T) = A*T^n*exp(-Ea/(R_gas*T)) at the runtime T). Event
    # rate R = k(T)*max(mu1-mu0, 0) [backbone bonds]; the parent pool is
    # debited the bond-weighted breaking chain and BOTH end-radical daughter
    # pools ({label}_rad_primary_end / {label}_rad_secondary_end, created at
    # model setup by the producer) are credited the uniform-cut fragments.
    # One-way (no reverse leg); releases NO gas species. Mutually exclusive
    # with k_scission > 0 (same random backbone-break event; double-count),
    # with radical_qssa_unzip (QSSA initiation IS random backbone
    # homolysis; double-count), and (round 67) with k_unzip > 0 (legacy
    # closed-chain monomer-loss channel vs radical-end pools feeding
    # explicit beta-scission/unzip chemistry; double-carried
    # depolymerization).
    k_homolysis: Optional[Dict[str, float]] = None

    # End-radical DEPROPAGATION kernel (adjudicated round 74 SS2, the run-6
    # no-outlet wall fix): normalized Arrhenius triplet {A, n, Ea} per
    # validate_k_depropagation (SI convention: A [s^-1], Ea [J/mol]; NOT a
    # bare scalar -- the RHS evaluates k_dep(T) = A*T^n*exp(-Ea/(R_gas*T))
    # at the runtime T) or None. Consumption channel for RADICAL-END pools
    # (one active radical end per chain): per unzip event one monomer
    # volatile is released to monomer_poly_index (gas amount basis).
    #   R    = k_dep(T) * mu0 * g     [g = smooth exhaustion gate, == 1 in
    #                                  the realizable region mean DP >= 1]
    #   dmu1 = -R;  gas = +R (the SAME float -- exact mass invariant)
    #   dmu2 = -k_dep*(2*mu1 - mu0)   [smooth positive-part on the excess]
    #   dmu0 = -k_dep*N1, N1 = mu0 * P(DP=1) from the EXISTING
    #          discrete/gamma closure (never a permanent dmu0 = 0 -- r74:
    #          that stalls at a one-repeat-per-chain residue)
    # Mutually exclusive on one pool with k_unzip > 0 (legacy
    # phenomenological form of the SAME chain-end monomer-release event),
    # with radical_qssa_unzip (its depropagation block IS lumped chain-end
    # depropagation) and with k_homolysis (multi-generation homolysis of
    # radical-ended chains is DEFERRED, r74 SS3; and a closed-chain
    # initiation pool has no radical end to depropagate). Requires a
    # resolvable monomer_poly_index (released units must land somewhere).
    k_depropagation: Optional[Dict[str, float]] = None

    # Side-group homolysis initiation kernel (FR1-K1, adjudicated round 70):
    # normalized channel LIST per validate_side_group_homolysis (each entry
    # exactly {label, A, n, Ea, site_selector, sites_per_unit, gas_product};
    # site_selector is the round-72 REQUIRED structural site selector) or
    # None. Per
    # None. SGH kernel-v2 (explicit Br-inventory depletion): per channel with
    # s = sites_per_unit and k(T) = A*T^n*exp(-Ea/(R_gas*T)) per site, the
    # reacting chain stays INTACT (no chain cut) and the removable X leaves
    # via a solver-owned auxiliary inventory Z (Z(0) = s*mu1(0), dZ/dt = -k*Z)
    # rather than a moment-debited feature pool. The carrier's mu0/mu1/mu2 are
    # STATIONARY (no debit, no feature-pool credit); the gas X radical is
    # credited +R and the carrier's condensed mass is booked the lost-X mass
    # directly (mu1*MW - Sum_row max(0, s*mu1 - Z)*M_X). v2 spawns NO feature
    # pool. Inventory-depletion SGH requires a STATIONARY carrier mu1 (the §1h
    # build guard refuses explicit_dp, tail_kinetics and resolved row-flux
    # touching the carrier; k_scission, verified mu1-preserving, MAY coexist).
    # A carrier with k_unzip > 0, radical_qssa_unzip, k_homolysis or
    # k_depropagation instead moves each reacting chain, intact, into the
    # channel's X-loss feature daughter '{label}_sidegrp_{channel}'
    # (side_group_accounting, I-257; the v1 law, the daughter carrying the
    # exact per-chain defect M_X).
    side_group_homolysis: Optional[List[Dict[str, object]]] = None
    # Resolved core index of each channel's gas_product species, aligned
    # with side_group_homolysis (PolymerPool.to_config resolves it; a
    # kernel-enabled pool without resolved indices is a hard config error --
    # the emitted X radical would silently vanish).
    side_group_gas_indices: Optional[Tuple[int, ...]] = None

    # Catalog-driven pool kernels (slice S5 a + b): normalized term lists per
    # validate_kernel_homolysis (this pool's chains cut into the `daughter`
    # end-radical pool) and validate_kernel_recombination (this pool's
    # chains join into the `product` pool), or None.
    kernel_homolysis: Optional[List[Dict[str, object]]] = None
    kernel_recombination: Optional[List[Dict[str, object]]] = None
    # Pair transfer (S5 c + e) per validate_kernel_transfer: one chain of
    # this pool and one of its partner change pool, lengths unchanged.
    kernel_transfer: Optional[List[Dict[str, object]]] = None
    # Pool move (S5 d) per validate_kernel_move: one chain of this pool
    # changes pool, length unchanged (backbiting and its reverse).
    kernel_move: Optional[List[Dict[str, object]]] = None
    # Propagation (I-266) per validate_kernel_propagation: one chain of this
    # pool adds one monomer unit, its monomer_poly_index species consumed.
    kernel_propagation: Optional[List[Dict[str, object]]] = None
    # Interior scission (I-273) per validate_kernel_interior_scission: every
    # chain of this pool breaks into fragments in its product pools.
    kernel_interior_scission: Optional[List[Dict[str, object]]] = None

    # Exact per-chain mass defect [g/mol] of an X-loss feature pool (FR1-K1
    # mass contract, the round-70 #1 P1 trap): the feature pool keeps the
    # parent's monomer_mw_g_mol (the chain transferred INTACT), but every
    # feature chain lost exactly ONE X atom (v1: no multi-loss cascade), so
    # the pool's EXACT condensed mass is
    #     mass [g] = mu1*monomer_mw_g_mol - mu0*chain_mass_defect_g_mol
    # (see condensed_mass_g). The producer pins it to M_X of the spawning
    # channel's gas_product; _flatten_side_group_state hard-errors on any
    # kernel destination whose defect does not pin M_X exactly. 0.0 on
    # ordinary pools (the accessor then reduces to mu1*MW).
    chain_mass_defect_g_mol: float = 0.0

    # An X-loss feature daughter's parent pool label (I-260), or None. Such a
    # pool carries its parent's backbone kernels and its spawning channel as
    # side_group_homolysis, booked in-pool (side_group_accounting): its Z
    # inventory holds the channel's intact sites on its chains, the first
    # loss per chain being its chain_mass_defect_g_mol and the rest
    # s*mu1 - Z - mu0. Its radical_qssa_unzip shares one radical population
    # with the parent's. side_loss_monomer_index is the core index of the
    # modified monomer (the unit less its group) its defect units release,
    # required when it carries a release kernel.
    side_loss_parent: Optional[str] = None
    side_loss_monomer_index: Optional[int] = None

    # Chain end-group bookkeeping (I-141). ``repeat_unit_atoms`` is the
    # element count of one repeat unit; ``end_group_atoms`` the element count
    # of the TWO declared end groups of one chain (head + tail), both as
    # ((symbol, count), ...). The moment closure books condensed atoms as
    # mu1*unit only; the physical inventory is mu1*unit + mu0*end groups, and
    # the atoms of chains that vanish (or are minted by scission) are carried
    # by one non-feeding "orphaned end groups" accumulator per pool (an
    # appended aux slot, see _book_end_group_ledger). A pool with
    # ``end_group_atoms`` None has no ledger; mixing declared and undeclared
    # pools in one system is refused at setup.
    repeat_unit_atoms: Optional[Tuple[Tuple[str, int], ...]] = None
    end_group_atoms: Optional[Tuple[Tuple[str, int], ...]] = None

    # True for a pool whose config is derived from an engine-spawned daughter
    # (rmgpy.polymer.is_engine_spawned_pool_daughter, set by
    # derive_daughter_pool_configs). Such a pool is born empty by design --
    # zero moments are its initial condition -- so set_initial_conditions
    # does not warn that it has none. Deck pools keep False.
    engine_spawned: bool = False

    # Custom Kinetics Hook
    # f(T, P, mu0, mu1, mu2, mu3) -> (dmu0_dt, dmu1_dt, dmu2_dt, small_species_sources)
    # Units:
    #   mu* are intensive moments [mol/m^3_poly]
    #   dmu*/dt are [mol/m^3_poly/s]
    #   small_species_sources: Dict[core_index -> mol/m^3_poly/s]
    # The solver multiplies these volumetric rates by V_poly to update extensive state [mol] / [mol/s].
    # Contract: If provided, this handles "Chemical Event" derivatives.
    # It does NOT handle the "Physical Handshake" flux.
    tail_kinetics: Optional[Callable[[float, float, float, float, float, float],
                                     Tuple[float, float, float, Dict[int, float]]]] = None

    def condensed_mass_g(self, mu0_mol, mu1_mol):
        """EXACT condensed mass [g] of this pool at the given extensive
        moments [mol] (FR1-K1 mass contract): mu1*monomer_mw_g_mol minus the
        per-chain X-loss defect mu0*chain_mass_defect_g_mol. For ordinary
        pools (defect 0) this is the plain moment-derived mass; for X-loss
        feature pools it carries the one-lost-X-per-chain correction so that
        d(condensed mass)/dt + d(gas X mass)/dt = 0 holds exactly on the
        side-group kernel's contributions (pinned by test). Also works on
        rates: feeding d(mu)/dt returns d(mass)/dt."""
        return (mu1_mol * self.monomer_mw_g_mol
                - mu0_mol * self.chain_mass_defect_g_mol)

    def __repr__(self):
        return f"<PolymerPool '{self.label}' xs={self.xs}>"


@dataclass(frozen=True)
class MassTransferConfig:
    """
    Explicit mass transfer between a gas-phase species and its dissolved counterpart.
    """
    gas_index: int
    poly_index: int
    K: float    # K = C_poly_eq / C_gas_eq
    kLa: float  # [1/s]


# ======================================================================================
# HELPERS
# ======================================================================================

def _safe_mu3_from_mu012(mu0: float, mu1: float, mu2: float) -> float:
    """Log-Lagrange closure μ3 = μ0·(μ2/μ1)³, with a realizability guard.

    A valid chain-length distribution (k ≥ 1) always satisfies μ1 ≥ μ0. If
    solver noise or a bad source term pushes the state out of the realizable
    cone, the (μ2/μ1)³ extrapolation can explode into a DASSL singularity, so
    we return 0.0 (no closure contribution) rather than amplifying garbage.

    N5b boundary-consistent closure (round-62 DASSL-hang root cause): the
    guard above is exact, but AT the all-monomer boundary μ1 = μ0 (every
    chain has length 1) the physical limit is μ3 → μ1 (Σn³ = Σn = μ1 when
    every n = 1), not 0.0 -- and solver noise routinely lands a μ1 a
    rounding ULP below μ0 right at that boundary (μ0, μ1 agreeing to 12+
    digits per the round-62 investigation). Flatly returning 0.0 for any
    μ1 < μ0, however small the violation, injects an O(k_scission·μ1)
    spurious jump into dμ2/dt against an error weight orders of magnitude
    smaller, which is what stalled DASSL's corrector. Within
    MU3_CLOSURE_BOUNDARY_REL of the boundary we instead return the
    boundary-consistent limit μ1; genuine deep violations (beyond the
    band) still return 0.0 unchanged.
    """
    if mu0 <= SMALL_EPS or mu1 <= SMALL_EPS or mu2 <= SMALL_EPS:
        return 0.0
    if mu1 < mu0:  # unrealizable: more chains than monomer units
        if mu0 - mu1 <= MU3_CLOSURE_BOUNDARY_REL * mu0:
            return mu1  # noise-scale violation astride the boundary
        return 0.0      # genuine deep realizability violation
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        try:
            ln_mu3 = 3.0 * np.log(mu2) - 3.0 * np.log(mu1) + np.log(mu0)
            if ln_mu3 > LN_EXP_OVERFLOW_GUARD:
                return float("inf")
            return float(np.exp(ln_mu3))
        except (FloatingPointError, ValueError):
            return 0.0

def _explicit_moment_contributions(V_poly: float,
                                   y: np.ndarray,
                                   explicit_dp_to_species_index: Dict[int, int],
                                   ) -> Tuple[float, float, float]:
    """Calculate explicit oligomer contributions to moments."""
    if V_poly <= 0: return 0.0, 0.0, 0.0
    mu0, mu1, mu2 = 0.0, 0.0, 0.0

    for n, idx in explicit_dp_to_species_index.items():
        if n <= 0: continue
        if not (0 <= idx < len(y)):
            raise IndexError(f"Explicit species index {idx} out of bounds.")

        Nn = max(0.0, y[idx]) / V_poly
        mu0 += Nn
        mu1 += n * Nn
        mu2 += (n * n) * Nn
    return mu0, mu1, mu2


def _gamma_params_from_mu012(mu0: float, mu1: float, mu2: float) -> Optional[Tuple[float, float]]:
    if mu0 <= SMALL_EPS or mu1 <= SMALL_EPS or mu2 <= SMALL_EPS:
        return None
    mean = mu1 / mu0
    if mean <= SMALL_EPS: return None
    pdi = (mu2 * mu0) / (mu1 * mu1)

    # Monodisperse fallback trigger
    if not np.isfinite(pdi) or pdi <= 1.0 + 1e-6:
        return None

    k = 1.0 / (pdi - 1.0)
    theta = mean / k
    if k <= 0.0 or theta <= 0.0 or not np.isfinite(k) or not np.isfinite(theta):
        return None
    return float(k), float(theta)


def _gamma_prob_conditional_hybrid(m: int, xs: int, k: float, theta: float) -> float:
    if m <= 0:
        return 0.0
    try:
        from scipy.special import gammainc
        def cdf(x):
            return float(gammainc(k, x/theta))
        F_b = cdf(m + 0.5)
        F_a = cdf(max(0.0, m - 0.5))
        F_cut = cdf(xs + 0.5)

        bin_prob = max(0.0, F_b - F_a)
        tail_prob = max(0.0, 1.0 - F_cut)
        if tail_prob <= 1e-12: return 0.0
        return bin_prob / tail_prob
    except ImportError:
        return _discrete_gamma_fallback(m, xs, k, theta)


def _discrete_gamma_fallback(target: int, xs: int, k: float, theta: float) -> float:
    if target <= xs:
        return 0.0
    def logpdf(n):
        x = float(n)
        if x <= 0: return -1e9
        return (k - 1.0) * math.log(x) - (x / theta) - k * math.log(theta) - math.lgamma(k)

    mean = k * theta
    sigma = math.sqrt(k * theta * theta)
    nmax = int(max(xs + 1, mean + 6.0 * sigma + 20.0))
    nmax = min(nmax, xs + 2000)

    logs = []
    target_idx = -1
    for i, n in enumerate(range(xs + 1, nmax + 1)):
        val = logpdf(n)
        logs.append(val)
        if n == target:
            target_idx = i

    if not logs:
        return 0.0
    max_log = max(logs)
    denom = sum(math.exp(v - max_log) for v in logs)
    if denom <= 0.0 or target_idx == -1:
        return 0.0
    return math.exp(logs[target_idx] - max_log) / denom


# End-radical DEPROPAGATION kernel (adjudicated round 74 SS2) smooth-gate
# width, dimensionless in mean DP. sp(x) = x^3/(x^2 + W^2) is a C2 smooth
# positive-part: exactly 0 for x <= 0, asymptotically x - W^2/x for x >> W.
# The gate g = 1 - sp(1 - mean) is therefore EXACTLY 1 in the realizable
# region mean >= 1 (the healthy law is exact, no perturbation) and rolls
# off C2-smoothly below it (r74 SS5: no hard max(...,0) cliff at
# exhaustion -- the DASPK grind / IDID=-7 class). This gate alone left a
# residual monomer release of k_dep*mu0*W^2 at the pathological mu1 = 0,
# mu0 > 0 noise state (it bottoms out at W^2/(1 + W^2), not at zero);
# since I-065 that residual is closed to EXACTLY zero by the additional
# _release_units_gate factor, which is itself C2 and exactly 1.0 for
# mean >= 1 so it costs the smoothness nothing. The chain count still
# drains at the full -k_dep*mu0 there, so the state cannot become a stiff
# no-outlet grind.
KDEP_GATE_WIDTH = 1.0e-2


def _smooth_pos(x: float, w: float) -> float:
    """C2 smooth positive-part x^3/(x^2 + w^2): 0 for x <= 0, ~x for
    x >> w, with value/first/second derivative all continuous at 0
    (strictly smoother than the max(x, 0) clamps elsewhere)."""
    if x <= 0.0:
        return 0.0
    return x * x * x / (x * x + w * w)


def _release_units_gate(mu0: float, mu1: float) -> float:
    """Availability gate on chain-end monomer RELEASE: the fraction of the
    chain-end event rate the pool actually holds repeat units for.

    Every chain-end release channel sets its event rate from the CHAIN count
    (one active end per chain), never from the repeat units on offer. Inside
    the realizable cone mu1 >= mu0 >= 0 that is right -- every chain carries
    at least one unit, so every event has one to give. OUTSIDE the cone it
    fabricates: at mu1 = 0, mu0 > 0 the rate is still k*mu0 > 0, which
    debits mu1 further negative AND emits that much gas monomer with no
    repeat unit behind it. Defence in depth for I-055: the realizable-cone
    invariant should not be the only thing between the solver and mass it
    never had.

    Returns a C2 quintic smootherstep in t = mu1/mu0:
        t <= 0  -> 0.0 EXACTLY  (no units, no release, no gas)
        t >= 1  -> 1.0 EXACTLY  (on/inside the cone: bit-for-bit no-op,
                                 x*1.0 == x in IEEE-754)
        else    -> 10t^3 - 15t^4 + 6t^5, with value, first AND second
                   derivative continuous at both ends.
    The quintic is not decoration. A hard min(mu0, mu1) limiter is exactly
    zero at mu1 = 0 too, but introduces a first-derivative kink at mean
    DP = 1 -- the boundary a fully-unzipping pool legitimately ENDS on --
    and that kink is a real cliff: it fails the existing
    test_deprop_smooth_exhaustion_gate_no_cliff pin (release rate must move
    by < 1e-4 relative across mean = 1; a linear rolloff moves it by 1e-3)
    and reopens the DASPK grind / IDID=-7 class the r74 SS5 gate was
    designed against. The quintic is C2 there and exactly 1.0 above it.

    mu0 <= 0 returns 0.0: no chains means no chain ends, hence no events.
    On the legacy channel that replaces a NEGATIVE release rate (k*mu0 < 0,
    i.e. gas emitted as a negative amount) with none at all.

    Applies to the release only, never to the I-055 chain-termination debit
    on mu0: that debit must keep draining at full strength off the cone, or
    a mean-DP < 1 noise state has no way back (r74: dmu0 = 0 stalls the
    residue instead of terminating it)."""
    cdef double t
    if mu0 <= 0.0:
        return 0.0
    if mu1 >= mu0:
        return 1.0
    if mu1 <= 0.0:
        return 0.0
    t = mu1 / mu0
    return t * t * t * (10.0 + t * (6.0 * t - 15.0))


def _deprop_dp1_fraction(mu0: float, mu1: float, mu2: float) -> float:
    """DP=1 chain fraction p1 (N1 = mu0 * p1) for the end-radical
    depropagation kernel (r74 SS2: 'Estimate N1 using the existing
    discrete/gamma-style moment closure').

    Gamma leg -- the EXISTING closure machinery: shape/scale from
    _gamma_params_from_mu012 (k = 1/(PDI-1), theta = mean/k), DP=1 bin
    probability from _gamma_prob_conditional_hybrid(1, 0, k, theta), i.e.
    [F(1.5) - F(0.5)] / [1 - F(0.5)] on the half-integer-bin discretization
    conditioned on n >= 1 (the same convention the hybrid handshake uses).

    Smooth terminal boundary floor (disclosed per r74): the gamma leg is
    unavailable for degenerate/monodisperse shapes (PDI <= 1 + 1e-6 returns
    no params), so p1 is floored by a C1 smoothstep 1 -> 0 over mean DP in
    [1, 2]. The floor is what guarantees NO stall: any state at or below
    one repeat unit per chain has p1 = 1, so dmu0 = -k_dep*mu0 terminates
    the residue instead of freezing it (r74: no permanent dmu0 = 0).

    Returns a value in [0, 1]; smooth (max of two continuous legs) in the
    moment state."""
    if mu0 <= SMALL_EPS:
        return 0.0
    mean = mu1 / mu0
    t = mean - 1.0
    if t <= 0.0:
        p_floor = 1.0
    elif t >= 1.0:
        p_floor = 0.0
    else:
        p_floor = 1.0 - (3.0 * t * t - 2.0 * t * t * t)
    p_gamma = 0.0
    params = _gamma_params_from_mu012(mu0, mu1, mu2)
    if params:
        k_shape, theta = params
        p_gamma = _gamma_prob_conditional_hybrid(1, 0, k_shape, theta)
    p1 = p_gamma if p_gamma > p_floor else p_floor
    if p1 < 0.0:
        return 0.0
    if p1 > 1.0:
        return 1.0
    return p1


# ======================================================================================
# THERMO REFERENCE-STATE TRIPWIRE (spec 2026-06-11)
# ======================================================================================

# Standard-state constants for the unpaired reference-state magnitude U.
# P0 is the NASA-polynomial gas standard pressure; C0 = 1 mol/m^3 is the
# standard concentration that makes log10(P0/(R*T*C0)) dimensionless -- at
# 1000 K the term is 1.08 decades. Omit C0 in a re-derivation in other units
# and the formula merely LOOKS wrong (spec §5.1).
REFERENCE_STATE_P0 = 1.0e5      # Pa
REFERENCE_STATE_C0 = 1.0        # mol/m^3
# Bounds (spec §6 -- chemistry decisions on record, not defaults): census
# above the measured benign ceiling (0.33 decades + small-margin items);
# refuse with a wide margin below the pathological floor (~10 decades),
# grounded in bimodality plus the §5.2 sign argument. The 0.5-3.0 band is
# deliberately wide and deliberately EMPTY on today's measurements; anything
# that ever lands in it is news AND the trigger to compute rotation into U
# properly (spec §5.2 escalation, on record in multi_pool_design.md §5.2).
REFERENCE_STATE_CENSUS_DECADES = 0.5
REFERENCE_STATE_REFUSE_DECADES = 3.0
# Chain-scale window slack on top of the largest pool monomer MW. ONE window
# (max pool monomer MW + slack) serves BOTH the physically-melt class
# tag-branch gate (spec §5.1, C3-amended) and the §5.3 provenance
# counterparty predicate -- one definition, two uses, so the class and the
# sensor cannot drift apart.
REFERENCE_STATE_MW_SLACK_G_MOL = 10.0


def _sackur_tetrode_s_trans(mw_kg_mol, T):
    """Exact Sackur-Tetrode translational entropy [J/mol/K] of an ideal gas
    of molar mass ``mw_kg_mol`` [kg/mol] at ``T`` [K] and P = REFERENCE_STATE_P0."""
    m = mw_kg_mol / constants.Na
    return constants.R * (
        math.log((2.0 * constants.pi * m * constants.kB * T / constants.h ** 2) ** 1.5
                 * constants.kB * T / REFERENCE_STATE_P0) + 2.5)


def _unpaired_reference_decades(reactant_melt_mws, product_melt_mws, T):
    """Unpaired reference-state magnitude U [decades] for one reversible
    reaction (spec 2026-06-11 §5.1):

        U = |sum_prod,melt S_trans - sum_react,melt S_trans| / (R ln10)
            + |dn_melt| * log10(P0 / (R*T*C0))

    Inputs are the molar masses [kg/mol] of the PHYSICALLY-MELT participants
    (condensed-mask OR is_polymer_proxy-tagged with chain-scale MW -- spec
    §5.1, C3-amended) on each side, multiplicity
    included; gas participants are excluded because their gas reference
    state is simply CORRECT. Paired same-mass chains cancel in the signed
    sum exactly as they do in the real thermo -- no pairing optimization is
    needed. Translational-only BY DESIGN: the measured U distribution is
    bimodal with a four-decade gap (benign <= 0.33, pathological >= ~10), so
    the omitted bounded rotational term (+1.6-3.1 decades, same sign where
    it matters) cannot flip a decision against the §6 thresholds. Escalation
    trigger (spec §5.2): a census entry in the 1-5 decade band means rotation
    must be computed into U properly."""
    s_sum = 0.0
    for mw in product_melt_mws:
        s_sum += _sackur_tetrode_s_trans(mw, T)
    for mw in reactant_melt_mws:
        s_sum -= _sackur_tetrode_s_trans(mw, T)
    dn_melt = len(product_melt_mws) - len(reactant_melt_mws)
    return (abs(s_sum) / (constants.R * math.log(10.0))
            + abs(dn_melt) * math.log10(
                REFERENCE_STATE_P0 / (constants.R * T * REFERENCE_STATE_C0)))


def _dual_axis_polymer_sized(mw_g_mol, heavy_count, pool_axes, units):
    """r89 dual-axis polymer-sized test at the solver seam -- the EXACT
    mirror of rmgpy.polymer._discrete_is_polymer_sized measured against the
    configured pools, on pre-plumbed numbers instead of structures (the
    solver's PolymerPoolConfig sidecars and consumer-world species carry no
    Molecule):

    * mass axis:  ``mw_g_mol`` [g/mol] vs ``units * pool monomer MW``;
    * heavy axis: ``heavy_count`` (non-H atoms; -1 = unknown, e.g. a
      consumer-world species whose artifact carried no composition) vs
      ``units * pool monomer_heavy_atoms`` (0 = unknown denominator).

    Returns ``(sized, missing)``. ``sized`` is True iff at least one pool
    decides True on BOTH computable axes (any-pool, matching
    _discrete_is_chain_scale_proxy_derived). ``missing`` is None when every
    consulted pool DECIDED (True somewhere, or computable evidence decided
    False everywhere); otherwise it names the uncomputable axis/axes of a
    pool whose computable evidence did not already decide False -- the r85
    P2(a) undecidable case: an uncomputable axis NEVER defers to the other
    (a single-axis melt call is exactly the run-9 false positive in reverse),
    the answer degrades to conservative-gas, and the caller must announce it
    (mirroring rmgpy.polymer._warn_impostor_axis_undecidable semantics)."""
    # r95 absolute chain-scale floor: mirror rmgpy.polymer._discrete_is_polymer_sized
    # EXACTLY (the drift guard test_solver_dual_axis_gate_agrees_with_polymer_module_classifier
    # depends on it). ONE source of truth for the two constants -- imported from
    # rmgpy.polymer, never re-hardcoded here (lazy import, matching the tripwire's
    # lazy import of _IMPOSTOR_DISCRETE_MONOMER_UNITS; avoids the circular
    # rmgpy.polymer <-> rmgpy.solver.polymer top-level import).
    from rmgpy.polymer import ABS_CHAIN_SCALE_MW, ABS_CHAIN_SCALE_HEAVY
    missing = None
    for monomer_mw, monomer_heavy in pool_axes:
        mass_computable = monomer_mw > 0.0 and mw_g_mol > 0.0
        mass_sized = mass_computable and mw_g_mol >= max(units * monomer_mw,
                                                         ABS_CHAIN_SCALE_MW)
        heavy_computable = monomer_heavy > 0 and heavy_count >= 0
        heavy_sized = heavy_computable and heavy_count >= max(units * monomer_heavy,
                                                              ABS_CHAIN_SCALE_HEAVY)
        if mass_computable and heavy_computable:
            if mass_sized and heavy_sized:
                return True, None
            continue  # both axes computable: decided False for this pool
        # At least one axis uncomputable for this pool: undecidable UNLESS a
        # computable axis already decided False (then it is a decided False).
        if not any([mass_computable and not mass_sized,
                    heavy_computable and not heavy_sized]):
            missing = ("mass and structure"
                       if not mass_computable and not heavy_computable
                       else ("mass" if not mass_computable else "structure"))
    return False, missing


def _assert_chain_scale_melt_member(label, mw_kg_mol, heavy_count,
                                    gas_classified, pool_axes, units):
    """Cannot-happen leak guard (spec §5.1 C3 amendment 2026-06-11, r89
    dual-axis form), called for every species entering the melt sum. A
    gas-classified species can be physically-melt ONLY via the tag branch,
    whose dual-axis polymer-sized conjunct (_dual_axis_polymer_sized: mass
    AND heavy axes, both computable, both >= units monomer-equivalents
    against >= 1 pool) is part of the CLASS DEFINITION -- so a gas-classified
    member below the threshold (or with an undecidable axis) cannot reach
    the sum through the amended gate. A tagged below-threshold species being
    EXCLUDED by the gate is expected and silent (the family.py:1657
    over-tagging fingerprint, H2 on every proxy-touching reaction; the PP
    run-9 DP-2 hexadiene at 2.0 monomer-equivalents); the raise here fires
    only if such a species REACHES the sum (a future refactor recomputing
    membership without the conjunct), converting the mistake into a loud
    CLASSIFICATION error instead of a silently-large U and a misattributed
    reference-state refusal. Condensed-branch members (gas_classified=False;
    pool-configured by input) are exempt from the gate and from this
    guard."""
    if not gas_classified:
        return
    sized, _missing = _dual_axis_polymer_sized(
        mw_kg_mol * 1000.0, heavy_count, pool_axes, units)
    if not sized:
        raise ValueError(
            "THERMO REFERENCE-STATE TRIPWIRE (classification leak): the "
            "is_polymer_proxy tag includes a non-chain species in the melt "
            f"sum ({label}, MW = {mw_kg_mol * 1000.0:.2f} g/mol, "
            f"heavy atoms = {heavy_count}, not polymer-sized at "
            f"{units} monomer-equivalents on both axes against any pool); "
            "physically-melt class definition violated -- this is a "
            "classification leak, NOT a thermo problem; do not respond by "
            "touching reference states. See the proxy-tag propagation chain "
            "(family.py) and the invariant section of "
            "docs/multi_pool_design.md.")


def _thermo_provenance(spc):
    """Classify a species' thermo source from its comment: 'library', 'gav',
    or None (no thermo / unrecognized -- never warned on). Substrings pinned
    against rmgpy/data/thermo.py: 'Thermo library: <label>' (:1818) and
    'Thermo group additivity estimation: ...' (:232/:2237; gav_keywords
    :2845). 'Thermo library' is checked FIRST: an HBI radical estimated from
    a library parent ('Thermo library: X + radical(Y)') classifies 'library',
    which is correct -- a library-parent HBI resolves through the same parent
    and PRESERVES the structural cancellation (spec §2)."""
    comment = getattr(getattr(spc, "thermo", None), "comment", "") or ""
    if "Thermo library" in comment:
        return "library"
    if "group additivity" in comment:
        return "gav"
    return None


# ======================================================================================
# SOLVER CLASS
# ======================================================================================

class HybridPolymerSystem(ReactionSystem):
    """
    Biphasic Hybrid Reactor Solver.

    Residual returns mol/s.
    Enforces 'Produce-then-Transfer' policy for interphase transport.
    """

    def __init__(
        self,
        T,
        P,
        initial_mole_fractions: Dict,
        V_poly: float,
        polymer_pools: List[PolymerPoolConfig],
        mass_transfer: Optional[List[MassTransferConfig]] = None,
        polymer_species_labels=None,
        gas_species_mask: Optional[np.ndarray] = None,
        prospective_gas_mask: Optional[np.ndarray] = None,
        prospective_classifier=None,
        prospective_condensed_edge_daughter_classifier=None,
        allow_default_prospective_edge: bool = False,
        constant_gas_volume: bool = False,
        V_gas0: Optional[float] = None,
        initial_polymer_moments: Optional[Dict[str, Tuple[float, float, float]]] = None,
        initial_explicit_species: Optional[Dict[str, Dict[int, float]]] = None,
        termination=None,
        sensitive_species=None,
        sensitivity_threshold: float = 1e-3,
        sens_conditions=None,
        const_spc_names: Optional[List[str]] = None,
        pdep_collision_reaction_indices: Optional[np.ndarray] = None,
        pdep_collider_kinetics: Optional[List] = None,
        collider_efficiencies: Optional[np.ndarray] = None,
        debug_check_realizability: bool = False,
        allow_unpaired_reference_state: bool = False,
        declared_melt_convention: bool = False,
        allow_unstamped_proxy_rows: bool = False,
        qssa_kout_floor_s: float = 0.0,
        polymer_scoped_jacobian: bool = False,
    ):
        super().__init__(termination=termination,
                         sensitive_species=sensitive_species,
                         sensitivity_threshold=sensitivity_threshold)
        cdef np.ndarray species_to_pool_indices

        self.T = Quantity(T)
        self.P = Quantity(P)
        # Deck-tied floor for the census-only QSSA k_out threshold policy
        # (round-20 increment 7): callers may set it from the RMG filter
        # threshold; the rebuild census uses
        # max(qssa_kout_floor_s, 1/terminationTime).
        self.qssa_kout_floor_s = float(qssa_kout_floor_s or 0.0)
        # M18.2d softclamp decay length (finalized from atol in
        # initialize_model; SMALL_EPS default keeps the extension
        # numerically OFF until tolerances are known).
        self._softclamp_lam = SMALL_EPS
        self.initial_mole_fractions = initial_mole_fractions

        self.V_poly = float(V_poly)
        if self.V_poly <= 0.0:
            raise ValueError(f"V_poly must be > 0, got {self.V_poly}")

        self.constant_gas_volume = bool(constant_gas_volume)
        self.V_gas0 = float(V_gas0) if V_gas0 is not None else None
        self.V_gas = 0.0
        self.V = 0.0

        self.polymer_pools = list(polymer_pools)
        self.mass_transfer = list(mass_transfer) if mass_transfer else []
        self.gas_species_mask = gas_species_mask
        # Item 17 (spec 2026-06-12 SS3(a)/(b)): stage-1 seed for the
        # prospective mask over chain(core, edge), normally computed by
        # to_solver_object with the SAME config-keyed classifier as
        # gas_species_mask (one classifier, longer list). Stored as a seed;
        # initialize_model rebuilds self.prospective_gas_mask on EVERY call.
        # RIDER R3: prospective_gas_mask is GATE-INPUT ONLY -- nothing but
        # the residual product gates (and riders R1/R2, which verify and
        # report them) may read it; gas_species_mask stays the sole source
        # of truth for every other phase behavior (concentrations, volumes,
        # validation, the thermo tripwire's melt class, the artifact mask
        # harvest).
        self._prospective_gas_mask_seed = prospective_gas_mask
        self.prospective_gas_mask = None
        # Item 17 A5-2 (spec 2026-06-12 SS3(a)/(d)): the §3(a) staging fix.
        # _prospective_classifier is a CALLABLE f(species_list) -> bool array
        # (the bound polymerPhase.get_gas_mask), plumbed onto the engine by
        # to_solver_object so initialize_model can re-run stage-1 over the
        # LIVE chain(core, edge) on EVERY call (base.pyx:simulate already
        # re-runs initialize_model with the live edge, so no reactor-level
        # trigger is needed). Used for CLASSIFICATION ONLY -- R3-clean, never a
        # phase verdict for any core behavior. None on direct-test/runner
        # construction. _allow_default_prospective_edge is the opt-in flag that
        # marks a build as legitimately fallback-permitted (direct-test/runner
        # with no blueprint phase object); R1-EDGE raises iff a build is NOT so
        # flagged AND its prospective edge suffix was default-filled.
        # _prospective_edge_provenance is the per-build int8 marker over the
        # edge entries (1 = stage-1-classified, 0 = default-filled).
        self._prospective_classifier = prospective_classifier
        # Spec 2026-06-29: callable f(combined_species) -> set of qualifying
        # daughter BASE labels (the bound PolymerPhase.get_condensed_edge_daughter_bases),
        # plumbed by to_solver_object. Re-run over the LIVE chain(core, edge) in
        # initialize_model on EVERY call (never a frozen set -> cannot go stale
        # on the engine-reuse path). None on direct-test/runner construction.
        self.prospective_condensed_edge_daughter_classifier = \
            prospective_condensed_edge_daughter_classifier
        self._allow_default_prospective_edge = bool(allow_default_prospective_edge)
        # r71 FIX 4 escape hatch (test/runner-only, mirroring
        # allow_default_prospective_edge): permits the legacy mu1-only
        # NONE -> UNRESOLVED demotion for unstamped LIVE proxy rows instead
        # of the production hard-fail at initialize_model. A production
        # build (rmgpy.rmg construction) leaves this default-False.
        self._allow_unstamped_proxy_rows = bool(allow_unstamped_proxy_rows)
        self._prospective_edge_provenance = None
        self.const_spc_names = const_spc_names or []
        self.const_spc_indices = None
        self.sens_conditions = sens_conditions

        self.pdep_collision_reaction_indices = pdep_collision_reaction_indices if pdep_collision_reaction_indices is not None else np.array([], int)
        self.pdep_collider_kinetics = pdep_collider_kinetics if pdep_collider_kinetics is not None else []
        self.collider_efficiencies = collider_efficiencies if collider_efficiencies is not None else np.array([[]], float)
        self.initial_polymer_moments = initial_polymer_moments or {}

        # Opt-in moment-realizability diagnostic (off by default; zero cost when
        # off). When on, the residual logs ONCE per pool if its moment state
        # leaves the realizable cone (μ1≥μ0≥0, μ0·μ2≥μ1²) — a fast way to localize
        # a bad moment source term. It only logs (never raises): the μ3 closure
        # already degrades gracefully on out-of-cone states.
        self.debug_check_realizability = bool(debug_check_realizability)
        self._realizability_warned = set()

        # Thermo reference-state tripwire (spec 2026-06-11 §7). The override
        # is the deck author's CONSCIOUS assertion that their thermo handles
        # unpaired reference states -- the eventual per-deck switch-on point
        # for full melt consistency. Census/diagnostics are repopulated by
        # _reference_state_tripwire on every initialize_model rebuild:
        # reference_state_census = [(str(rxn), U)] for U > census bound;
        # reference_state_max_decades = max U over reversible core reactions.
        self.allow_unpaired_reference_state = bool(allow_unpaired_reference_state)
        self.reference_state_census = []
        self.reference_state_max_decades = 0.0
        # Declared single convention (polymer_phase residence, shared-chemistry
        # slice S4): the deck states that every melt-resident (condensed)
        # species takes the gas reference at its melt concentration, as the
        # kinetic Monte Carlo does. A reversible row wholly condensed is then
        # listed in reference_state_declared_convention [(rxn, U)] instead of
        # refused; rows across the phase boundary are refused as before.
        self.declared_melt_convention = bool(declared_melt_convention)
        self.reference_state_declared_convention = []
        self.initial_explicit_species = initial_explicit_species or {}
        self.polymer_species_labels = set(polymer_species_labels) if polymer_species_labels else set()

        self.jacobian_matrix = None

        # 18.2e P5 (option (d), spar round-44): scoped-RHS user-side FD
        # Jacobian -- feature-gated, DEFAULT OFF. When OFF, no `jacobian`
        # attribute is EVER set on this instance, so pydas' attribute-
        # presence detection (checked inside DASPK.initialize, proven by
        # probe) sees nothing and DASPK keeps its internal dense FD: the
        # shipped default path is bit-identical to today. When requested,
        # _setup_scoped_jacobian() (called from initialize_model just
        # before initialize_solver) validates the layout assumptions and
        # either arms self.jacobian = self._scoped_fd_jacobian or falls
        # back LOUDLY to no-jacobian (guard 5, one warning).
        self._scoped_jacobian_requested = bool(polymer_scoped_jacobian)
        self._scoped_jacobian_active = False
        # J-scope flag read by residual(): True ONLY inside
        # _scoped_fd_jacobian's FD sweep. Every real (DASPK-originated)
        # residual call runs with False -> the RHS path is untouched.
        self._jac_scope = False
        # Guard-2 tripwire: True while the LAST residual-kernel activity
        # was a scoped J-sweep (edge_*/census diagnostics stale from a
        # perturbed, edge-skipped call); every real residual call clears
        # it. DDASPK always follows a Jacobian (re)form with corrector
        # residual calls, so after any step()/advance() this must read
        # False (pinned by test) -- i.e. the diagnostics the enlargement
        # loop reads after step() always come from a REAL residual call.
        self._scoped_jac_diag_stale = False
        # DDAWTS weight vectors captured at _setup_scoped_jacobian time
        # (pydas' own atol/rtol arrays are cdef-private).
        self._jac_wt_atol = None
        self._jac_wt_rtol = None
        # i061: pre-floor copy of atol_array (see MOMENT_EWT_FLOOR_K). The
        # moment error-weight floor raises atol_array on the mu-slots for DASPK
        # error control; chemistry noise-floor consumers read this preserved
        # array so their behaviour is bitwise unchanged. Set in initialize_model.
        self._chem_atol_array = None

        self._scratch_C_gas = None
        self._scratch_C_poly = None
        self._scratch_dn_dt = None
        # Weak-link U-state (milestones ii+iii): number of appended per-pool
        # U slots and their dU/dt scratch. 0/None until initiate_tolerances /
        # validate_configuration lay the state out; legacy-only systems keep
        # num_qssa_u == 0 and the pre-milestone neq == n_core layout.
        self.num_qssa_u = 0
        self._scratch_du_dt = None
        # SGH kernel-v2 auxiliary Br-inventory depletion state (Z slots): the
        # number of appended per-(pool, SGH-channel-row) Z slots and their
        # dZ/dt scratch. Z = moles of remaining removable X sites for the row;
        # Z(0) = sgh_sites[row]*mu1(0). 0/None until initiate_tolerances /
        # validate_configuration lay the state out; legacy/non-SGH systems keep
        # num_sgh_z == 0 and the residual path is bit-for-bit unchanged.
        self.num_sgh_z = 0
        self._scratch_dz_dt = None
        # I-141 orphaned-end-group ledger. L_p = mol of chains' worth of pool-p
        # end groups that left the moment accounting (net of chains minted):
        #     L_p = mu0_p(0) - mu0_p(t) - H_p
        # where H_p = chains handed to the pool's explicit-oligomer species
        # (which carries its own ends, so those orphan nothing). mu0 is
        # already a state, so L needs NO state of its own; only H_p does, and
        # only on pools with an explicit-DP handshake: one trailing aux slot
        # (state order [core | U | Z | H]), output-only, NEVER read into a
        # rate. (A passive slot on every pool was tried and rejected: any
        # extra slot changes DASPK's step sequence, measured 2e-9 on the
        # poly102 crash-window replay.) end_ledger_rel[p] is the slot offset
        # inside the H block (-1: none).
        self.num_end_ledger = 0
        self.end_ledger_rel = np.full(len(polymer_pools), -1, dtype=np.int32)
        self._scratch_dl_dt = None
        self._core_species_atoms = None

        self.pool_mu0_indices = np.full(len(self.polymer_pools), -1, dtype=np.int32)
        self.pool_mu1_indices = np.full(len(polymer_pools), -1, dtype=np.int32)

        # Exhaustion-tail conditioning state (r81 B). None/empty until
        # initialize_model lays out the state vector and the atol-derived
        # per-moment floors; the census set is per-rebuild (cleared there).
        self._pool_exhausted = None
        self._pool_worst_trial_excursion = None
        self._pool_mu_floors = None
        self._exhaustion_census_emitted = set()
        # I-244 relative accepted floor and pool_exhausted events (see
        # EXHAUSTION_RELATIVE_K), laid out per rebuild by initialize_model.
        self._pool_mu_peak = None
        self._pool_mu_rtol = None
        self._pool_ended = None
        self._pool_exhaustion_residue = None
        self._pool_exhaustion_z = {}
        self.pool_exhaustion_events = []
        # I-090 cone-margin dead-band narrowing: the rtol half of the b1c - 1
        # noise band, rebound from the live tolerances by
        # initiate_tolerances. The base-class default stands in so the law is
        # well-defined on a system whose tolerances were never initiated.
        self._cone_b1_rtol = 1.0e-8

    def initiate_tolerances(self, atol=1e-16, rtol=1e-8, sensitivity=False,
                            sens_atol=1e-6, sens_rtol=1e-4):
        """Extend the base ODE layout with one U slot per weak-link pool
        (milestones ii+iii).

        Layout contract: y = [core species amounts (mol) | U slots (mol)],
        one appended slot per pool whose radical_qssa_unzip channel carries
        the weak-link vocabulary, in pool order. Pools with a legacy-only
        channel (or none) get NO slot: the legacy state layout is UNCHANGED
        (neq == n_core, pinned bitwise by the golden RHS test). Slot
        positions are bound in _flatten_radical_qssa_state (qssa_u_slot);
        the count here keys on the same "initiation_allyl" key presence, and
        the flattening pass cross-checks the two counts.

        The DASPK sensitivity layout interleaves per-species blocks
        (neq = n_core*(n_rxn+n_core+1)); a trailing U slot does not fit it,
        so sensitivity + weak-link is refused loudly rather than silently
        corrupting the sensitivity state vector.
        """
        ReactionSystem.initiate_tolerances(self, atol, rtol, sensitivity,
                                           sens_atol, sens_rtol)
        # I-090: bind the rtol half of the cone-margin dead-band b1c noise
        # band to the tolerance actually in force. Taken from the scalar rtol
        # here for the same reason _softclamp_lam takes the scalar atol --
        # these are scalar-tolerance decks; the per-state generalization is
        # deferred with its sibling's. The floor half of the band comes from
        # the pool's own r81 floors, read at the call site.
        self._cone_b1_rtol = float(rtol)
        n_weak = 0
        for pool in self.polymer_pools:
            q = pool.radical_qssa_unzip
            if isinstance(q, dict) and "initiation_allyl" in q:
                n_weak += 1
        self.num_qssa_u = n_weak
        # SGH kernel-v2 (§1c): one trailing Z aux-state slot per (pool, SGH
        # channel row). Counted here on the RAW configs (before validation),
        # mirroring the U-slot count posture above; _flatten_side_group_state
        # cross-checks this count against the slots it binds. State-vector
        # order is [core | U slots | Z slots]; the Z block is concatenated
        # AFTER the U block so the appended-slot ordering matches the absolute
        # slot indices bound in flatten (num_core_species + num_qssa_u +
        # z_count).
        n_sgh_z = 0
        for pool in self.polymer_pools:
            channels = pool.side_group_homolysis
            if channels and side_group_accounting(pool) != SIDE_GROUP_DAUGHTER:
                n_sgh_z += len(channels)
        self.num_sgh_z = n_sgh_z
        # I-141: declared-and-undeclared end groups in one system are refused
        # rather than guessed at; a handoff-counter slot H is allocated only
        # on declared pools that own an explicit-DP species, AFTER U and Z.
        declared = [p for p in self.polymer_pools if p.end_group_atoms]
        self.end_ledger_rel = np.full(len(self.polymer_pools), -1, dtype=np.int32)
        if declared and len(declared) != len(self.polymer_pools):
            # Not fatal at setup: hand-built pool configs predate the end-group
            # declaration (every deck/artifact producer now sets it). The
            # reporters refuse by name (_require_end_group_atoms) instead of
            # guessing a composition for the undeclared pools.
            logging.warning(
                "Pool(s) %s declare no end_group_atoms while other pools do; "
                "end_group_ledger/orphaned_end_atoms/physical_elemental_inventory "
                "will refuse until every pool declares them.",
                [p.label for p in self.polymer_pools if not p.end_group_atoms])
        n_led = 0
        for i_p, pool in enumerate(self.polymer_pools):
            if pool.end_group_atoms and pool.explicit_dp_to_species_index:
                self.end_ledger_rel[i_p] = n_led
                n_led += 1
        self.num_end_ledger = n_led
        if n_weak == 0 and n_sgh_z == 0 and n_led == 0:
            return
        if sensitivity:
            raise ValueError(
                f"{n_weak} weak-link allyl/U-state pool(s) and {n_sgh_z} "
                f"side-group-homolysis Z aux-state slot(s) configured with "
                f"sensitivity analysis enabled: the DASPK sensitivity state "
                f"layout interleaves per-species blocks and has no room for "
                f"the appended trailing U/Z slots. Disable sensitivity or "
                f"remove the weak-link / side_group_homolysis channel(s).")
        if n_weak > 0:
            self.neq += n_weak
            self.atol_array = np.concatenate(
                [self.atol_array, np.full(n_weak, atol, dtype=float)])
            self.rtol_array = np.concatenate(
                [self.rtol_array, np.full(n_weak, rtol, dtype=float)])
        if n_sgh_z > 0:
            self.neq += n_sgh_z
            self.atol_array = np.concatenate(
                [self.atol_array, np.full(n_sgh_z, atol, dtype=float)])
            self.rtol_array = np.concatenate(
                [self.rtol_array, np.full(n_sgh_z, rtol, dtype=float)])
        if n_led > 0:
            self.neq += n_led
            self.atol_array = np.concatenate(
                [self.atol_array, np.full(n_led, atol, dtype=float)])
            self.rtol_array = np.concatenate(
                [self.rtol_array, np.full(n_led, rtol, dtype=float)])

    def validate_configuration(self):
        """Strict validation of indices and masks."""
        n_core = self.num_core_species

        if self.gas_species_mask is None:
            raise ValueError("gas_species_mask must be set before validation.")

        if len(self.gas_species_mask) != n_core:
            raise ValueError(f"Gas mask len {len(self.gas_species_mask)} != num_core {n_core}")

        moment_indices = set()

        for pool in self.polymer_pools:
            for idx in pool.mu_indices:
                if not (0 <= idx < n_core):
                    raise ValueError(f"Pool {pool.label} moment index {idx} out of range.")
                if self.gas_species_mask[idx]:
                    raise ValueError(f"Pool {pool.label} moment index {idx} masked as GAS.")
                moment_indices.add(idx)

            for dp, idx in pool.explicit_dp_to_species_index.items():
                if not (0 <= idx < n_core):
                    raise ValueError(f"Pool {pool.label} explicit DP={dp} index {idx} out of range.")
                if self.gas_species_mask[idx]:
                    raise ValueError(f"Pool {pool.label} explicit DP={dp} index {idx} masked as GAS.")

            if pool.monomer_poly_index is not None:
                if not (0 <= pool.monomer_poly_index < n_core):
                    raise ValueError(f"Pool {pool.label} monomer index out of range.")
                # INVERTED contract (incident 2026-07-03, design B-prime): the
                # released-monomer target is the deck's gas volatile and the
                # unzip/QSSA release deposits into the GAS amount basis. A
                # condensed-masked target would re-create the reference-state
                # conflation the tripwire refuses (one Species carrying two
                # phase residences).
                if not self.gas_species_mask[pool.monomer_poly_index]:
                    raise ValueError(
                        f"Pool {pool.label} monomer index masked as CONDENSED. "
                        f"The unzip/QSSA release target (monomer_product) is "
                        f"emitted to the GAS phase; it must be classified gas.")

            # Unzip-channel invariants. This is the LAST line of defense: the
            # deck helper (rmg/input.py), PolymerPool.to_config
            # (rmg/polymer_input.py) and the artifact runner
            # (tools/polymer_moments_runner.py) all guard the same shapes, but
            # a directly-constructed PolymerPoolConfig bypasses every one of
            # them. The residual drains the condensed moments unconditionally
            # when k_unzip > 0 (dmu1_dt -= k_unzip*mu0; dmu2_dt -=
            # k_unzip*(2*mu1 - mu0)) while the released-monomer emission is
            # gated on monomer_poly_index is not None -- so k_unzip > 0 with
            # no monomer index makes mass leave the condensed phase
            # un-conserved (the drained repeat units go nowhere).
            # Finiteness first: NaN passes BOTH the `< 0` and `> 0` gates as
            # False, so a non-finite k_unzip/k_scission would make the
            # channel SILENTLY INERT (a laundered no-op) while inf poisons
            # the residual. Mirror the QSSA triplet validator's posture.
            for _rate_name, _rate_val in (("k_unzip", pool.k_unzip),
                                          ("k_scission", pool.k_scission)):
                if not math.isfinite(_rate_val):
                    raise ValueError(
                        f"Pool {pool.label}: {_rate_name}={_rate_val!r} is not "
                        f"finite -- NaN/inf is not a valid rate constant (NaN "
                        f"would silently disable the channel; inf would poison "
                        f"the residual). Set {_rate_name} to a finite value "
                        f">= 0.")
            if pool.k_unzip < 0.0:
                raise ValueError(
                    f"Pool {pool.label}: k_unzip={pool.k_unzip:g} is negative -- "
                    f"a negative k_unzip is not a valid rate constant (every unzip "
                    f"consumer is gated on k_unzip > 0, so it would silently "
                    f"masquerade as a frozen pool). Set k_unzip >= 0.")
            if pool.k_unzip > 0.0 and pool.monomer_poly_index is None:
                raise ValueError(
                    f"Pool {pool.label}: k_unzip={pool.k_unzip:g} > 0 but "
                    f"monomer_poly_index is None. The unzip channel drains this "
                    f"pool's condensed moments unconditionally, and without a "
                    f"released-monomer emission target the drained mass leaves "
                    f"the condensed phase silently un-conserved. Wire "
                    f"monomer_poly_index (the released monomer's core index) or "
                    f"set k_unzip=0.")

            # Radical-QSSA unzip channel invariants (M1: the channel is inert
            # -- nothing in the residual reads it until the M2 rate law -- but
            # the solver is the LAST line of defense, same rationale as the
            # k_unzip guards above: a directly-constructed PolymerPoolConfig
            # bypasses the deck helper and PolymerPool.to_config, so the
            # solver re-validates the channel shape itself).
            if pool.radical_qssa_unzip is not None:
                # CONSUME the validator's normalized return and store it back
                # (review round 21, finding 1): a directly-constructed
                # PolymerPoolConfig with a minimal-but-valid channel would
                # otherwise pass validation but stay UN-normalized in storage
                # (missing efficiency/monomer_yield/transfer/basis keys), so
                # any future q["efficiency"] read would KeyError. The config
                # is a frozen dataclass, so the store-back uses
                # object.__setattr__ (the standard frozen-dataclass idiom; no
                # prior precedent in this codebase). deepcopy makes the stored
                # structure defensively independent of every caller-provided
                # sub-dict (finding 2's aliasing half): post-hoc mutation of
                # the caller's dict cannot reach the stored config.
                normalized_qssa = validate_radical_qssa_unzip(
                    pool.label, pool.radical_qssa_unzip)
                object.__setattr__(pool, 'radical_qssa_unzip',
                                   copy.deepcopy(normalized_qssa))
                # (The milestone-i anti-silent-no-op guard that refused
                # weak-link configs here was REMOVED by milestones ii+iii:
                # the U-state slot and the allyl initiation RHS below now
                # consume the vocabulary for real.)
                if pool.monomer_poly_index is None:
                    raise ValueError(
                        f"Pool {pool.label}: radical_qssa_unzip is configured "
                        f"but monomer_poly_index is None. The QSSA unzip "
                        f"channel releases monomer through the pool's existing "
                        f"monomer routing; without an emission target the "
                        f"depropagated repeat units would leave the condensed "
                        f"phase silently un-conserved. Wire monomer_poly_index "
                        f"(the released monomer's core index) or remove "
                        f"radical_qssa_unzip.")
                if pool.k_unzip > 0.0:
                    raise ValueError(
                        f"Pool {pool.label}: radical_qssa_unzip is configured "
                        f"AND k_unzip={pool.k_unzip:g} > 0. These are two "
                        f"representations of the SAME chain-end depropagation "
                        f"channel and are mutually exclusive on a pool "
                        f"(enabling both would double-count the unzip flux). "
                        f"Set k_unzip=0 or remove radical_qssa_unzip.")

            # Catalog-driven pool kernels (S5 a + b): normalized and stored
            # back like the channels above. kernel_homolysis IS backbone
            # homolysis, so k_homolysis, k_scission > 0 or the QSSA channel
            # (whose initiation is homolysis) on the same pool would
            # double-count the event.
            if pool.kernel_homolysis is not None:
                object.__setattr__(pool, 'kernel_homolysis',
                                   validate_kernel_homolysis(
                                       pool.label, pool.kernel_homolysis))
                for other, on in (("k_homolysis", pool.k_homolysis is not None),
                                  ("k_scission", pool.k_scission > 0.0),
                                  ("radical_qssa_unzip",
                                   pool.radical_qssa_unzip is not None)):
                    if on:
                        raise ValueError(
                            f"Pool {pool.label}: kernel_homolysis and {other} "
                            f"are mutually exclusive -- both carry the "
                            f"backbone homolysis of this pool's chains.")
            if pool.kernel_recombination is not None:
                object.__setattr__(pool, 'kernel_recombination',
                                   validate_kernel_recombination(
                                       pool.label, pool.kernel_recombination))
            if pool.kernel_transfer is not None:
                object.__setattr__(pool, 'kernel_transfer',
                                   validate_kernel_transfer(
                                       pool.label, pool.kernel_transfer))
            if pool.kernel_move is not None:
                object.__setattr__(pool, 'kernel_move',
                                   validate_kernel_move(
                                       pool.label, pool.kernel_move))
            if pool.kernel_propagation is not None:
                object.__setattr__(pool, 'kernel_propagation',
                                   validate_kernel_propagation(
                                       pool.label, pool.kernel_propagation))
                if pool.monomer_poly_index is None:
                    raise ValueError(
                        f"Pool {pool.label}: kernel_propagation is configured "
                        f"but monomer_poly_index is None: the monomer each "
                        f"event consumes has no core species, so the added "
                        f"units would come from nowhere.")
            if pool.kernel_interior_scission is not None:
                object.__setattr__(pool, 'kernel_interior_scission',
                                   validate_kernel_interior_scission(
                                       pool.label, pool.kernel_interior_scission))

            # Radical-homolysis kernel invariants (Stage 1, adjudicated round
            # 66). LAST line of defense, same rationale as above: the deck
            # helper and PolymerPool.to_config guard these shapes, but a
            # directly-constructed PolymerPoolConfig bypasses both. The
            # normalized triplet is stored back (frozen-dataclass store-back
            # idiom, mirroring the QSSA channel above) so the flattener can
            # trust its shape.
            if pool.k_homolysis is not None:
                normalized_khom = validate_k_homolysis(
                    pool.label, pool.k_homolysis)
                object.__setattr__(pool, 'k_homolysis',
                                   copy.deepcopy(normalized_khom))
                if pool.k_scission > 0.0:
                    raise ValueError(
                        f"Pool {pool.label}: k_homolysis is configured AND "
                        f"k_scission={pool.k_scission:g} > 0. Both "
                        f"parameterize the SAME random backbone-break event "
                        f"(homolysis IS random scission, with the products "
                        f"routed to the end-radical pools) and are mutually "
                        f"exclusive on a pool -- enabling both would "
                        f"double-count chain initiation. Set k_scission=0 or "
                        f"remove k_homolysis.")
                if pool.radical_qssa_unzip is not None:
                    raise ValueError(
                        f"Pool {pool.label}: k_homolysis is configured AND "
                        f"radical_qssa_unzip is configured. The QSSA "
                        f"channel's initiation block IS random backbone "
                        f"homolysis, so the two are mutually exclusive on a "
                        f"pool -- enabling both would double-count "
                        f"initiation. Remove one of them.")
                if pool.k_unzip > 0.0:
                    raise ValueError(
                        f"Pool {pool.label}: k_homolysis is configured AND "
                        f"k_unzip={pool.k_unzip:g} > 0. Legacy k_unzip is a "
                        f"phenomenological closed-chain monomer-loss "
                        f"channel, while k_homolysis creates radical-end "
                        f"pools that feed explicit beta-scission/unzip "
                        f"chemistry; the two are mutually exclusive on a "
                        f"pool -- enabling both would double-carry "
                        f"depolymerization. Set k_unzip=0 or remove "
                        f"k_homolysis.")

            # End-radical depropagation kernel invariants (adjudicated round
            # 74 SS2). LAST line of defense, same rationale as the kernels
            # above: the deck helper (parent-declared, daughter-applied) and
            # derive_daughter_pool_configs guard these shapes, but a
            # directly-constructed PolymerPoolConfig bypasses both. The
            # normalized triplet is stored back (frozen-dataclass store-back
            # idiom) so the flattener can trust its shape.
            if pool.k_depropagation is not None:
                normalized_kdep = validate_k_depropagation(
                    pool.label, pool.k_depropagation)
                object.__setattr__(pool, 'k_depropagation',
                                   copy.deepcopy(normalized_kdep))
                if pool.monomer_poly_index is None:
                    raise ValueError(
                        f"Pool {pool.label}: k_depropagation is configured "
                        f"but monomer_poly_index is None. The kernel "
                        f"releases ONE monomer volatile per unzip event; "
                        f"without a resolvable gas emission target the "
                        f"released units would leave the condensed phase "
                        f"silently un-conserved. Wire monomer_poly_index "
                        f"(the released monomer's core index) or remove "
                        f"k_depropagation.")
                if pool.k_unzip > 0.0:
                    raise ValueError(
                        f"Pool {pool.label}: k_depropagation is configured "
                        f"AND k_unzip={pool.k_unzip:g} > 0. Legacy k_unzip "
                        f"is the phenomenological scalar form of the SAME "
                        f"chain-end monomer-release event (and its law "
                        f"stalls at a permanent dmu0 = 0), so the two are "
                        f"mutually exclusive on a pool -- enabling both "
                        f"would double-carry depropagation. Set k_unzip=0 "
                        f"or remove k_depropagation.")
                if pool.radical_qssa_unzip is not None:
                    raise ValueError(
                        f"Pool {pool.label}: k_depropagation is configured "
                        f"AND radical_qssa_unzip is configured. The QSSA "
                        f"channel's depropagation block IS lumped chain-end "
                        f"depropagation, so the two are mutually exclusive "
                        f"on a pool -- enabling both would double-count the "
                        f"unzip flux. Remove one of them.")
                if pool.k_homolysis is not None:
                    raise ValueError(
                        f"Pool {pool.label}: k_depropagation is configured "
                        f"AND k_homolysis is configured. Multi-generation "
                        f"homolysis of radical-ended chains is DEFERRED "
                        f"(adjudicated round 74 SS3: internal homolysis of a "
                        f"radical-ended chain creates radical-state pools "
                        f"this producer cannot generate), and a closed-chain "
                        f"initiation pool has no radical end to "
                        f"depropagate; the two kernels are mutually "
                        f"exclusive on one pool. Remove one of them.")

            # Side-group homolysis kernel invariants (FR1-K1, adjudicated
            # round 70). LAST line of defense, same rationale as the two
            # kernels above: deck helper and PolymerPool.to_config guard
            # these shapes, but a directly-constructed PolymerPoolConfig
            # bypasses both. NOTE: a carrier whose chains another kernel moves
            # books its losses through X-loss feature daughters
            # (side_group_accounting, I-257); an inventory carrier (kernel-v2)
            # MAY coexist only with k_scission, which acts on backbone C-C
            # bonds and is verified mu1-preserving.
            if (not math.isfinite(pool.chain_mass_defect_g_mol)
                    or pool.chain_mass_defect_g_mol < 0.0):
                raise ValueError(
                    f"Pool {pool.label}: chain_mass_defect_g_mol="
                    f"{pool.chain_mass_defect_g_mol!r} must be a finite "
                    f"value >= 0 [g/mol] (the X-loss feature-pool mass "
                    f"contract; garbage here silently mints/destroys "
                    f"condensed mass).")
            if pool.side_group_homolysis is not None:
                # Shape + selector-vocabulary backstop only: the config
                # carries no monomer structure, so the round-72 STRUCTURAL
                # selector layer (match >= 1, sites_per_unit == match
                # count, disjoint atom sets) runs where the parsed monomer
                # lives -- deck helper, PolymerPool.to_config, producer.
                normalized_sgh = validate_side_group_homolysis(
                    pool.label, pool.side_group_homolysis)
                object.__setattr__(pool, 'side_group_homolysis',
                                   copy.deepcopy(normalized_sgh))
                # Beside k_unzip, radical_qssa_unzip, k_homolysis or
                # k_depropagation the carrier books its losses through X-loss
                # feature daughters (side_group_accounting, I-257): its
                # radicals stay on the daughter's chains and feed none of
                # those kernels, so nothing is carried twice.
                gi = pool.side_group_gas_indices
                if gi is None or len(gi) != len(normalized_sgh):
                    raise ValueError(
                        f"Pool {pool.label}: side_group_homolysis is "
                        f"configured but the gas-product core indices are "
                        f"not resolved (side_group_gas_indices={gi!r} for "
                        f"{len(normalized_sgh)} channel(s)). The kernel "
                        f"emits one X radical per event; without a resolved "
                        f"gas destination the ejected X would silently "
                        f"vanish (un-conserved mass).")
                for ci, g_idx in enumerate(gi):
                    if not (0 <= g_idx < n_core):
                        raise ValueError(
                            f"Pool {pool.label}: side_group_homolysis "
                            f"channel '{normalized_sgh[ci]['label']}' gas "
                            f"index {g_idx} out of range.")
                    if not self.gas_species_mask[g_idx]:
                        raise ValueError(
                            f"Pool {pool.label}: side_group_homolysis "
                            f"channel '{normalized_sgh[ci]['label']}' "
                            f"gas_product index {g_idx} is NOT masked as "
                            f"GAS -- the ejected X radical is a gas-phase "
                            f"species; emitting it into a condensed slot "
                            f"would corrupt the phase bookkeeping.")

        for mt in self.mass_transfer:
            if not (0 <= mt.poly_index < n_core):
                raise ValueError(f"Mass transfer poly_index {mt.poly_index} out of range.")
            if not (0 <= mt.gas_index < n_core):
                raise ValueError(f"Mass transfer gas_index {mt.gas_index} out of range.")

            if self.gas_species_mask[mt.poly_index]:
                raise ValueError(f"Mass transfer poly_index {mt.poly_index} is masked as GAS.")
            if not self.gas_species_mask[mt.gas_index]:
                raise ValueError(f"Mass transfer gas_index {mt.gas_index} is NOT masked as GAS.")

        # Moment isolation check
        ir = self.reactant_indices
        ip = self.product_indices
        all_rxn_indices = {i for i in itertools.chain(ir.flatten(), ip.flatten()) if i >= 0}
        overlap = moment_indices.intersection(all_rxn_indices)
        if overlap:
            raise ValueError(f"Configuration Error: Moment indices {overlap} appear in reaction stoichiometry. "
                             f"Moments must evolve only via tail_kinetics.")

        # Flatten the (now validated + normalized) radical_qssa_unzip channels
        # into solver-owned per-pool flat arrays. Runs AFTER the per-pool loop
        # above so every stored channel dict is guaranteed normalized.
        self._flatten_radical_qssa_state()

        # Flatten the (now validated + normalized) k_homolysis kernels and
        # resolve each enabled parent's end-radical daughter pools. Same
        # timing contract as the QSSA flattener above.
        self._flatten_homolysis_state()

        # Flatten the (now validated + normalized) side_group_homolysis
        # channels and seed each enabled carrier's per-channel Z inventory +
        # gas destinations (FR1-K1, kernel-v2). v2 spawns NO feature pool --
        # Br leaves through the Z depletion state, not a moment-debited
        # daughter. Same timing contract.
        self._flatten_side_group_state()

        # Flatten the (now validated + normalized) k_depropagation kernels
        # and resolve each enabled pool's released-monomer gas destination
        # (r74 SS2). Same timing contract.
        self._flatten_depropagation_state()

        # Flatten the (now validated + normalized) catalog-driven kernel
        # terms and resolve their target pools (S5 a + b). Same timing.
        self._flatten_catalog_kernel_state()

    def _flatten_catalog_kernel_state(self):
        """Flatten every pool's kernel_homolysis / kernel_recombination terms
        into per-pool lists of (coefficients, target moment slots); the RHS
        reads only these. HARD ERROR when a term's daughter or product pool
        is not configured: its fragments or products would vanish."""
        label_to_pool = {p.label: p for p in self.polymer_pools}
        # One repeat unit across every chain move: the moments count units,
        # so a move between pools of different units keeps the count and
        # books the wrong atoms and mass (a C2 -> C1 move loses carbon). A
        # unit declared on one side only cannot be checked: refused too.
        for pool in self.polymer_pools:
            for kernel, template, src, dst in catalog_kernel_connections(pool):
                a, b = label_to_pool.get(src), label_to_pool.get(dst)
                if a is None or b is None:
                    continue                    # named below
                atoms = [None if p.repeat_unit_atoms is None
                         else sorted(p.repeat_unit_atoms) for p in (a, b)]
                mass = [float(p.monomer_mw_g_mol or 0.0) or None for p in (a, b)]
                for what, (x, y), same in (
                        ("atoms", atoms, atoms[0] == atoms[1]),
                        ("molar mass [g/mol]", mass,
                         None not in mass and math.isclose(
                             mass[0], mass[1], rel_tol=1e-9, abs_tol=0.0))):
                    if (x is None and y is None) or same:
                        continue
                    raise ValueError(
                        f"Pool {pool.label}: {kernel} term {template!r} "
                        f"moves chains across a repeat unit: pool {src!r} "
                        f"has {what} {x!r}, pool {dst!r} {y!r} "
                        f"(the moments count units, not atoms; undeclared "
                        f"on one side cannot be checked).")
        # Each term's validity length, the shortest chain (periods) its site
        # law counts correctly: declared, or homolysis' own 2e - 1. Below it
        # the law miscounts (a chain of 4 periods under sites n - 7 counts
        # -3); that is the method's, reported here for the comparison side
        # (which knows the distribution) to weigh -- never corrected.
        self.kernel_validity = []
        for pool in self.polymer_pools:
            for kernel in CATALOG_KERNELS:
                for term in getattr(pool, kernel, None) or ():
                    n = term.get(KERNEL_VALIDITY_KEY)
                    source = "declared"
                    if n is None and kernel == "kernel_homolysis":
                        n, source = (term["min_fragment"]
                                     + term.get("partner_min_fragment", term["min_fragment"])
                                     + term.get("residue", 0) - 1), "derived"
                    elif n is None:
                        source = "undeclared"
                    self.kernel_validity.append(
                        {"pool": pool.label, "kernel": kernel,
                         "template": term["template"],
                         "valid_from_periods": n, "source": source})
        self._kernel_guard_warned = set()
        self.kernel_homolysis_terms = [[] for _ in self.polymer_pools]
        # I-278: two-daughter homolysis terms, (term, daughter slots, partner
        # slots, S1 and S2 coefficients of each daughter's length sums).
        self.kernel_homolysis_two_terms = [[] for _ in self.polymer_pools]
        self.kernel_recombination_terms = [[] for _ in self.polymer_pools]
        for i, pool in enumerate(self.polymer_pools):
            for kernel, key, terms, out in (
                    ("kernel_homolysis", "daughter", pool.kernel_homolysis,
                     self.kernel_homolysis_terms[i]),
                    ("kernel_recombination", "product",
                     pool.kernel_recombination,
                     self.kernel_recombination_terms[i])):
                for term in terms or ():
                    slots = []
                    for k in (key, "partner", "partner_daughter"):
                        if k not in term:
                            continue
                        target = label_to_pool.get(term[k])
                        if target is None:
                            raise ValueError(
                                f"Pool {pool.label}: {kernel} term "
                                f"{term['template']!r} names the {k} pool "
                                f"{term[k]!r}, which is not among the "
                                f"configured solver pools "
                                f"{sorted(label_to_pool)}.")
                        slots.append(tuple(target.mu_indices))
                    if kernel == "kernel_recombination":
                        # (term, product slots, partner slots or None)
                        out.append((dict(term), slots[0],
                                    slots[1] if len(slots) > 1 else None))
                    elif "partner_daughter" in term:
                        e1, e2 = term["min_fragment"], term["partner_min_fragment"]
                        r = term["residue"]
                        self.kernel_homolysis_two_terms[i].append(
                            (dict(term), slots[0], slots[1],
                             _range_power_sums(e1, r + e2),
                             _range_power_sums(e2, r + e1)))
                    else:
                        out.append((dict(term), slots[0]))
        # kernel_transfer: (term, partner, product, partner_product) moment
        # slots; every named pool must be configured or its chains vanish.
        self.kernel_transfer_terms = [[] for _ in self.polymer_pools]
        for i, pool in enumerate(self.polymer_pools):
            for term in pool.kernel_transfer or ():
                slots = []
                for key in ("partner", "product", "partner_product"):
                    target = label_to_pool.get(term[key])
                    if target is None:
                        raise ValueError(
                            f"Pool {pool.label}: kernel_transfer term "
                            f"{term['template']!r} names the {key} pool "
                            f"{term[key]!r}, which is not among the "
                            f"configured solver pools "
                            f"{sorted(label_to_pool)}.")
                    slots.append(tuple(target.mu_indices))
                self.kernel_transfer_terms[i].append((dict(term), *slots))
        # kernel_move: (term, product) moment slots, the same refusal.
        self.kernel_move_terms = [[] for _ in self.polymer_pools]
        for i, pool in enumerate(self.polymer_pools):
            for term in pool.kernel_move or ():
                target = label_to_pool.get(term["product"])
                if target is None:
                    raise ValueError(
                        f"Pool {pool.label}: kernel_move term "
                        f"{term['template']!r} names the product pool "
                        f"{term['product']!r}, which is not among the "
                        f"configured solver pools {sorted(label_to_pool)}.")
                self.kernel_move_terms[i].append(
                    (dict(term), tuple(target.mu_indices)))
        # kernel_propagation (I-266): the terms; the monomer is the pool's
        # monomer_poly_index (validate_configuration requires it).
        self.kernel_propagation_terms = [
            [dict(term) for term in pool.kernel_propagation or ()]
            for pool in self.polymer_pools]
        # kernel_interior_scission (I-273): (term, [(q, product moment
        # slots)]); every product pool must be configured or its fragments
        # vanish.
        self.kernel_interior_scission_terms = [[] for _ in self.polymer_pools]
        for i, pool in enumerate(self.polymer_pools):
            for term in getattr(pool, "kernel_interior_scission", None) or ():
                products = []
                for product in term["products"]:
                    target = label_to_pool.get(product["pool"])
                    if target is None:
                        raise ValueError(
                            f"Pool {pool.label}: kernel_interior_scission term "
                            f"{term['template']!r} names the product pool "
                            f"{product['pool']!r}, which is not among the "
                            f"configured solver pools {sorted(label_to_pool)}.")
                    products.append(([list(row) for row in product["moments"]],
                                     tuple(target.mu_indices)))
                self.kernel_interior_scission_terms[i].append((dict(term), products))
        self._assert_catalog_kernels_count_once()

    def _row_kernel_channel(self, combined_rxns, int r):
        """(pool index, carried channel, kernels) when a kernel on the pool
        row ``r`` decomposes carries its channel (I-253 channel identity,
        rmgpy.polymer.kernels_carrying_channel), else None. The pool is the
        configured pool of the row's one-species side; the channel is the
        row's generation stamp, or its species' reading; the row moves
        repeat units when it is a scission or chip row, or a volatile
        ejection with a nonzero moment shift (reaction_eject_units_moment)."""
        from rmgpy.polymer import (kernel_carried_channel, kernels_carrying_channel,
                                   reaction_channel)
        cdef int n_core = self.num_core_species
        ir = [int(i) for i in self.reactant_indices[r] if i != -1]
        ip = [int(i) for i in self.product_indices[r] if i != -1]
        if len(ir) == 1 and len(ip) == 2:
            whole, pieces = ir[0], ip
        elif len(ir) == 2 and len(ip) == 1:
            whole, pieces = ip[0], ir
        else:
            return None
        if whole >= n_core or self.species_to_pool_indices[whole] < 0:
            return None
        p = int(self.species_to_pool_indices[whole])
        arch = self.reaction_flux_archetype[r]
        moves = (arch == FLUX_SCISSION_FRAGMENT or arch == FLUX_DISCRETE_CHIP
                 or (arch == FLUX_VOLATILE_EJECTION
                     and self.reaction_eject_units_moment[r] != 0.0))
        carried = kernel_carried_channel(reaction_channel(combined_rxns[r]), moves)
        released = [self._channel_species[i] for i in pieces
                    if not (i < n_core and self.species_to_pool_indices[i] >= 0)]
        kernels = kernels_carrying_channel(self.polymer_pools[p], carried, released)
        if not kernels:
            return None
        return p, carried, kernels

    def _refuse_unbalanced_proxy_rows(self, combined_rxns):
        """Refuse by name every live pool-coupled row (a core pool proxy among
        its species) whose species do not balance
        (rmgpy.polymer.proxy_balance_defect: the net charge changes, or the
        species read another channel than the row's generation stamp) --
        archetype NONE, zero flux -- and record it in ``proxy_balance_census``
        ({reaction, reason, charge, channel}); warn-once per row. A row with
        one species on both sides (a chip, a chain-end fold) books a change
        that species does not show, and is not read for radicals."""
        from rmgpy.polymer import proxy_balance_defect, warn_once_proxy_balance
        cdef int n_core = self.num_core_species
        self.proxy_balance_census = []
        for r in range(len(combined_rxns)):
            if self.reaction_refused[r]:
                continue
            ir = {int(i) for i in self.reactant_indices[r] if i != -1}
            ip = {int(i) for i in self.product_indices[r] if i != -1}
            if not any(i < n_core and self.species_to_pool_indices[i] >= 0 for i in ir | ip):
                continue
            entry = proxy_balance_defect(combined_rxns[r], folded=bool(ir & ip))
            if entry is None:
                continue
            self.reaction_refused[r] = 1
            self.reaction_flux_archetype[r] = FLUX_NONE
            self.proxy_balance_census.append(entry)
            warn_once_proxy_balance(entry)

    def _refuse_kernel_superseded_rows(self, combined_rxns):
        """Refuse by name every live row whose channel a kernel on the pool
        it decomposes carries (archetype NONE, zero flux) and record it in
        ``kernel_supersession_census`` ({reaction, pool, channel, kernels,
        reason}); one warn-once census line per pool. Rows whose daughter
        endpoint the enlarge boundary held back (already refused at pool
        resolution, ``_held_superseded``) are recorded first, by the same
        rule read on their species."""
        from rmgpy.polymer import _reaction_census_label, warn_once_kernel_supersession
        self.kernel_supersession_census = []
        by_pool = {}
        for r, rxn_r, (pool_label, carried, kernels) in getattr(self, "_held_superseded", ()):
            label = _reaction_census_label(rxn_r)
            self.kernel_supersession_census.append(
                {"reaction": label, "pool": pool_label, "channel": carried,
                 "kernels": list(kernels), "reason": "superseded:" + carried})
            by_pool.setdefault(pool_label, []).append((label, carried, tuple(kernels)))
        for r in range(len(combined_rxns)):
            if self.reaction_refused[r]:
                continue
            hit = self._row_kernel_channel(combined_rxns, r)
            if hit is None:
                continue
            p, carried, kernels = hit
            self.reaction_refused[r] = 1
            self.reaction_flux_archetype[r] = FLUX_NONE
            label = _reaction_census_label(combined_rxns[r])
            self.kernel_supersession_census.append(
                {"reaction": label, "pool": self.polymer_pools[p].label,
                 "channel": carried, "kernels": list(kernels),
                 "reason": "superseded:" + carried})
            by_pool.setdefault(self.polymer_pools[p].label, []).append(
                (label, carried, tuple(kernels)))
        for pool_label, rows in by_pool.items():
            warn_once_kernel_supersession({"pool": pool_label, "superseded_rows": rows})

    def _census_kernel_channel_double_counts(self, combined_rxns):
        """The double-count census, by the same rule: every LIVE row whose
        channel a kernel on the pool it decomposes carries
        (``kernel_double_count_census``: {reaction, pool, channel, kernels},
        warn-once per row). The supersession above leaves none, so an entry
        means a row went live past it. Diagnostic-only."""
        from rmgpy.polymer import _reaction_census_label, warn_once_kernel_double_count
        self.kernel_double_count_census = []
        for r in range(len(combined_rxns)):
            if self.reaction_refused[r]:
                continue
            hit = self._row_kernel_channel(combined_rxns, r)
            if hit is None:
                continue
            p, carried, kernels = hit
            entry = {"reaction": _reaction_census_label(combined_rxns[r]),
                     "pool": self.polymer_pools[p].label,
                     "channel": carried, "kernels": list(kernels)}
            self.kernel_double_count_census.append(entry)
            warn_once_kernel_double_count(entry)

    def _assert_catalog_kernels_count_once(self):
        """One event, one count: refuse a catalog template listed on two
        pools (each listing fires), and a live core reaction row that moves
        chains along a catalog term's connection (source -> target pool;
        both directions of a reversible row) -- the generated row and the
        kernel would both integrate the same events."""
        owner = {}
        for pool in self.polymer_pools:
            for kernel in CATALOG_KERNELS:
                for term in getattr(pool, kernel, None) or ():
                    first = owner.setdefault((kernel, term["template"]),
                                             pool.label)
                    if first != pool.label:
                        raise ValueError(
                            f"{kernel} template {term['template']!r} is "
                            f"listed on pools {first!r} and {pool.label!r}: "
                            f"each listing fires, so its events would count "
                            f"twice.")
        index = {p.label: i for i, p in enumerate(self.polymer_pools)}
        moves = {}
        for pool in self.polymer_pools:
            for kernel, template, src, dst in catalog_kernel_connections(pool):
                moves.setdefault((index[src], index[dst]), []).append(
                    (kernel, template, src, dst))
        arch = getattr(self, "reaction_flux_archetype", None)
        if not moves or arch is None:
            return
        refused = getattr(self, "reaction_refused", None)
        reversible = getattr(self, "_core_row_reversible", None)
        for r in range(min(arch.shape[0], self.num_core_reactions)):
            if (refused is not None and refused[r]) or arch[r] not in (
                    FLUX_MIGRATION, FLUX_SCISSION_FRAGMENT,
                    FLUX_VOLATILE_EJECTION):
                continue
            s, d = int(self.reaction_src_pool[r]), int(self.reaction_dst_pool[r])
            legs = [(s, d)]
            if reversible is None or reversible[r]:
                legs.append((d, s))
            for leg in legs:
                for kernel, template, src, dst in moves.get(leg, ()):
                    raise ValueError(
                        f"{kernel} term {template!r} moves chains {src!r} "
                        f"-> {dst!r}, and so does the live core reaction "
                        f"row {r} (flux archetype {int(arch[r])}): both "
                        f"would integrate the same events. Keep one of "
                        f"them.")

    def _flatten_depropagation_state(self):
        """Flatten each pool's validated+normalized k_depropagation triplet
        into solver-owned per-pool flat arrays (r74 SS2 end-radical
        depropagation kernel).

        Contract (mirrors _flatten_homolysis_state):
        - Populated ONLY here, on the initialize_model ->
          validate_configuration path, from the normalized dict just stored.
        - The RHS reads ONLY these arrays, never the dict.

        Layout (index = pool position in self.polymer_pools):
        - kdep_enabled[i] (int8): 1 iff the pool has k_depropagation.
        - kdep_A/n/Ea[i]: the Arrhenius triplet (A [s^-1], Ea [J/mol]).
        - kdep_gas[i] (int32): core index of the released monomer volatile
          (== pool.monomer_poly_index, guaranteed non-None and gas-masked by
          validate_configuration); -1 on kernel-free pools.
        """
        n_pools = len(self.polymer_pools)
        self.kdep_enabled = np.zeros(n_pools, dtype=np.int8)
        self.kdep_A = np.zeros(n_pools, dtype=float)
        self.kdep_n = np.zeros(n_pools, dtype=float)
        self.kdep_Ea = np.zeros(n_pools, dtype=float)
        self.kdep_gas = np.full(n_pools, -1, dtype=np.int32)
        for i, pool in enumerate(self.polymer_pools):
            trip = pool.k_depropagation
            if trip is None:
                continue
            self.kdep_enabled[i] = 1
            self.kdep_A[i] = trip["A"]
            self.kdep_n[i] = trip["n"]
            self.kdep_Ea[i] = trip["Ea"]
            self.kdep_gas[i] = pool.monomer_poly_index

    def _flatten_homolysis_state(self):
        """Flatten each pool's validated+normalized k_homolysis triplet into
        solver-owned per-pool flat arrays, and resolve the two end-radical
        daughter pools' moment slots (Stage 1, adjudicated round 66).

        Contract (mirrors _flatten_radical_qssa_state):
        - Populated ONLY here, on the initialize_model ->
          validate_configuration path, from the normalized dict just stored.
        - The RHS reads ONLY these arrays, never the dict.

        Layout (index = pool position in self.polymer_pools):
        - khom_enabled[i] (int8): 1 iff the pool has k_homolysis configured.
        - khom_A/n/Ea[i]: the Arrhenius triplet (A [s^-1], Ea [J/mol]).
        - khom_prim_mu0/mu1/mu2[i] and khom_sec_mu0/mu1/mu2[i] (int32): the
          absolute core indices of the {label}_rad_primary_end /
          {label}_rad_secondary_end daughter pools' moment slots; -1 on
          kernel-free pools.

        HARD ERROR when a kernel-enabled pool's daughter pools are missing
        from the configured pool list: the producer creates them at model
        setup, so an unresolvable daughter here means the kernel's fragment
        credits would silently vanish (mass laundered out of the melt).
        """
        n_pools = len(self.polymer_pools)
        # warn-once registry for the runtime moment-shape guard (mirrors
        # self._realizability_warned): one loud line per pool per rebuild.
        self._khom_guard_warned = set()
        self.khom_enabled = np.zeros(n_pools, dtype=np.int8)
        self.khom_A = np.zeros(n_pools, dtype=float)
        self.khom_n = np.zeros(n_pools, dtype=float)
        self.khom_Ea = np.zeros(n_pools, dtype=float)
        self.khom_prim_mu0 = np.full(n_pools, -1, dtype=np.int32)
        self.khom_prim_mu1 = np.full(n_pools, -1, dtype=np.int32)
        self.khom_prim_mu2 = np.full(n_pools, -1, dtype=np.int32)
        self.khom_sec_mu0 = np.full(n_pools, -1, dtype=np.int32)
        self.khom_sec_mu1 = np.full(n_pools, -1, dtype=np.int32)
        self.khom_sec_mu2 = np.full(n_pools, -1, dtype=np.int32)

        label_to_pool = {p.label: p for p in self.polymer_pools}
        for i, pool in enumerate(self.polymer_pools):
            trip = pool.k_homolysis
            if trip is None:
                continue
            missing = [f"{pool.label}{suffix}"
                       for suffix in K_HOMOLYSIS_DAUGHTER_SUFFIXES
                       if f"{pool.label}{suffix}" not in label_to_pool]
            if missing:
                raise ValueError(
                    f"Pool {pool.label}: k_homolysis is enabled but the "
                    f"end-radical daughter pool(s) {missing} are not among "
                    f"the configured solver pools (expected "
                    f"'{pool.label}{K_HOMOLYSIS_DAUGHTER_SUFFIXES[0]}' and "
                    f"'{pool.label}{K_HOMOLYSIS_DAUGHTER_SUFFIXES[1]}'). The "
                    f"producer spawns both at model setup; without them the "
                    f"kernel's fragment credits would silently vanish.")
            self.khom_enabled[i] = 1
            self.khom_A[i] = trip["A"]
            self.khom_n[i] = trip["n"]
            self.khom_Ea[i] = trip["Ea"]
            prim = label_to_pool[
                f"{pool.label}{K_HOMOLYSIS_DAUGHTER_SUFFIXES[0]}"]
            sec = label_to_pool[
                f"{pool.label}{K_HOMOLYSIS_DAUGHTER_SUFFIXES[1]}"]
            (self.khom_prim_mu0[i], self.khom_prim_mu1[i],
             self.khom_prim_mu2[i]) = prim.mu_indices
            (self.khom_sec_mu0[i], self.khom_sec_mu1[i],
             self.khom_sec_mu2[i]) = sec.mu_indices

    def _flatten_side_group_state(self):
        """Flatten each pool's validated+normalized side_group_homolysis
        channel list into solver-owned per-channel row arrays (SGH kernel-v2,
        explicit Br-inventory depletion; adjudicated spar r103 / plan r104).

        Contract (mirrors _flatten_homolysis_state):
        - Populated ONLY here, on the initialize_model ->
          validate_configuration path, from the normalized channel lists.
        - The RHS reads ONLY these arrays, never the pool dicts.

        Layout (one row per (pool, channel)):
        - sgh_enabled[i] (int8, per pool): 1 iff pool i has channels.
        - sgh_row_start[i]/sgh_row_end[i] (int32, per pool): this pool's
          contiguous slice of the row arrays below (one row per channel).
        - sgh_A/n/Ea/sites[row]: the channel's Arrhenius triplet
          (A [s^-1 per site], Ea [J/mol]) and sites_per_unit.
        - sgh_gas[row] (int32): core index of the ejected X radical.
        - sgh_mw_X[row] (float): molar mass M_X [g/mol] of the ejected X
          radical (the channel's gas_product); the condensed-mass tally uses
          it to book Br leaving via Z as lost*M_X.
        - sgh_z_slot[row] (int32): ABSOLUTE ODE index of this row's Z
          aux-state slot (num_core_species + num_qssa_u + z_count); Z = moles
          of remaining removable X sites, Z(0) = sites*mu1(0), dZ/dt = -R.
          -1 on a daughter-accounting carrier (side_group_accounting, I-257),
          whose rows instead credit
        - sgh_dst_mu0/mu1/mu2[row] (int32): the X-loss feature daughter's
          moment slots (-1 on an inventory row).
        - sgh_channel_labels[row]: channel label (guard warn-once keys and
          error messages).

        v2 deletes the SGH X-loss FEATURE pool entirely: Br leaves the melt
        through the Z depletion state, not a moment-debited daughter pool, so
        the feature-pool resolution / mass-defect-pin / monomer_mw-divergence
        guards are gone (Br mass is now carried by sgh_mw_X on the Z path).

        HARD BUILD ERRORS (§1h + P1-B, never silent): an SGH carrier pool
        that ALSO carries a mu1-MOVING kernel or has any live resolved
        row-flux touching its mu1 is refused -- inventory-depletion SGH
        requires a STATIONARY carrier mu1 (Z(0) = sites*mu1(0) is seeded once
        and depletes independently; a moving mu1 would strand or double-count
        the Br inventory).
        """
        n_pools = len(self.polymer_pools)
        # warn-once registries keyed (pool label, channel label); one loud
        # line per key per rebuild. _sgh_guard_warned is retained for
        # sibling-kernel API parity; _sgh_inventory_census_emitted backs the
        # accepted-state Z-inventory in-band census (§1g).
        self._sgh_guard_warned = set()
        self._sgh_inventory_census_emitted = set()
        self.sgh_enabled = np.zeros(n_pools, dtype=np.int8)
        self.sgh_row_start = np.zeros(n_pools, dtype=np.int32)
        self.sgh_row_end = np.zeros(n_pools, dtype=np.int32)
        rows_A = []
        rows_n = []
        rows_Ea = []
        rows_sites = []
        rows_gas = []
        rows_mw_x = []
        rows_z_slot = []
        rows_dst = ([], [], [])
        self.sgh_channel_labels = []
        pool_index = {p.label: j for j, p in enumerate(self.polymer_pools)}

        # Absolute Z-slot base: state vector is [core | U slots | Z slots], so
        # the first Z slot follows the core block and every U slot. Mirrors
        # qssa_u_slot's absolute-index binding (num_core_species + u_count).
        z_base = self.num_core_species + int(self.num_qssa_u)
        z_count = 0
        # The Z slot of each X-loss daughter's in-pool row, keyed (daughter,
        # channel), bound in the same order as below: its parent's move
        # credits the sites a moved chain brings (I-260).
        in_pool_slot = {}
        for pool in self.polymer_pools:
            if (pool.side_group_homolysis
                    and side_group_accounting(pool) != SIDE_GROUP_DAUGHTER):
                for ch in pool.side_group_homolysis:
                    if side_group_accounting(pool) == SIDE_GROUP_IN_POOL:
                        in_pool_slot[(pool.label, ch["label"])] = (
                            z_base + len(in_pool_slot) + z_count)
                    else:
                        z_count += 1
        z_count = 0
        rows_in_pool = []
        rows_dst_z = []
        self.slr_row = np.full(n_pools, -1, dtype=np.int32)
        self.slr_mod_idx = np.full(n_pools, -1, dtype=np.int32)
        for i, pool in enumerate(self.polymer_pools):
            self.sgh_row_start[i] = len(rows_A)
            self.sgh_row_end[i] = len(rows_A)
            channels = pool.side_group_homolysis
            if not channels:
                continue
            accounting = side_group_accounting(pool)
            inventory = accounting != SIDE_GROUP_DAUGHTER
            if accounting == SIDE_GROUP_INVENTORY:
                # §1h mu1-mover exclusion (SOLVER-BUILD guard, deck-kernel
                # arm): inventory-depletion SGH requires a stationary carrier
                # mu1.
                self._assert_sgh_carrier_mu1_stationary(i, pool)
            elif accounting == SIDE_GROUP_IN_POOL:
                self._bind_side_loss_release(i, pool, channels, len(rows_A))
            for ci, ch in enumerate(channels):
                rows_A.append(ch["A"])
                rows_n.append(ch["n"])
                rows_Ea.append(ch["Ea"])
                rows_sites.append(ch["sites_per_unit"])
                rows_gas.append(int(pool.side_group_gas_indices[ci]))
                # M_X of the ejected X radical -- the value v1 threw away; v2
                # carries Br mass out of the melt through the Z inventory.
                m_x = _side_group_gas_mw_g_mol(ch["gas_product"])
                rows_mw_x.append(m_x)
                rows_in_pool.append(int(accounting == SIDE_GROUP_IN_POOL))
                if inventory:
                    rows_z_slot.append(z_base + z_count)
                    z_count += 1
                    for k in range(3):
                        rows_dst[k].append(-1)
                    rows_dst_z.append(-1)
                else:
                    rows_z_slot.append(-1)
                    dst = self._side_group_daughter(pool, ch["label"], m_x,
                                                    pool_index)
                    for k in range(3):
                        rows_dst[k].append(int(dst.mu_indices[k]))
                    rows_dst_z.append(in_pool_slot.get((dst.label, ch["label"]), -1))
                self.sgh_channel_labels.append(ch["label"])
            self.sgh_enabled[i] = 1
            self.sgh_row_end[i] = len(rows_A)

        self.sgh_A = np.asarray(rows_A, dtype=float)
        self.sgh_n = np.asarray(rows_n, dtype=float)
        self.sgh_Ea = np.asarray(rows_Ea, dtype=float)
        self.sgh_sites = np.asarray(rows_sites, dtype=float)
        self.sgh_gas = np.asarray(rows_gas, dtype=np.int32)
        self.sgh_mw_X = np.asarray(rows_mw_x, dtype=float)
        self.sgh_z_slot = np.asarray(rows_z_slot, dtype=np.int32)
        self.sgh_dst_mu0 = np.asarray(rows_dst[0], dtype=np.int32)
        self.sgh_dst_mu1 = np.asarray(rows_dst[1], dtype=np.int32)
        self.sgh_dst_mu2 = np.asarray(rows_dst[2], dtype=np.int32)
        self.sgh_in_pool = np.asarray(rows_in_pool, dtype=np.int8)
        self.sgh_dst_z = np.asarray(rows_dst_z, dtype=np.int32)

        # Layout cross-check (mirrors the U-slot drift guard): initiate_
        # tolerances allocated num_sgh_z Z slots from the raw configs; the
        # flattening above bound z_count of them. A mismatch means the pool
        # configs changed between ODE layout and validation -- the Z slots
        # would alias core/U state. Fail loud, never integrate a misaligned
        # vector.
        if z_count != self.num_sgh_z:
            raise ValueError(
                f"SGH Z-state layout drift: initiate_tolerances allocated "
                f"{self.num_sgh_z} Z slot(s) but validate_configuration bound "
                f"{z_count} side_group_homolysis channel row(s). The pool "
                f"configs changed between ODE layout and validation; rebuild "
                f"the reactor (initialize_model) instead of mutating pool "
                f"configs in place.")

    def _bind_side_loss_release(self, int pool_i, pool, channels, int row):
        """Bind an X-loss daughter's release split (I-260): when it carries
        a release kernel (radical_qssa_unzip, k_unzip, k_depropagation), a
        released unit is a monomer while its group is intact and the
        modified monomer (side_loss_monomer_index) once lost, read from the
        channel's intact-site inventory, which needs one channel of one site
        per unit (several sites per unit would release several modified
        monomers)."""
        if not any(p.label == pool.side_loss_parent for p in self.polymer_pools):
            raise ValueError(
                f"Pool {pool.label}: an X-loss daughter of "
                f"'{pool.side_loss_parent}', which is not configured; its "
                f"inventory is credited by that parent's moves.")
        if not (pool.radical_qssa_unzip is not None or pool.k_unzip > 0.0
                or pool.k_depropagation is not None):
            return
        if len(channels) != 1 or channels[0]["sites_per_unit"] != 1.0:
            raise ValueError(
                f"Pool {pool.label}: an X-loss daughter with a release kernel "
                f"needs one side_group_homolysis channel of one site per unit "
                f"(got {[(c['label'], c['sites_per_unit']) for c in channels]}): "
                f"its released units are a monomer or the one modified "
                f"monomer.")
        m = pool.side_loss_monomer_index
        if m is None or not (0 <= m < self.num_core_species) or not self.gas_species_mask[m]:
            raise ValueError(
                f"Pool {pool.label}: an X-loss daughter with a release kernel "
                f"releases its defect units as the modified monomer, but "
                f"side_loss_monomer_index={m!r} is not a gas core species; "
                f"those units would leave the melt un-conserved.")
        self.slr_row[pool_i] = row
        self.slr_mod_idx[pool_i] = m

    def _split_side_loss_release(self, int pool_i, int mono, double r_release,
                                 double mu1_mol, double V_poly, y,
                                 dict small_src, dz_dt):
        """An X-loss daughter's release r_release [mol/(m^3 s)], already
        booked to its monomer index ``mono``, split by its intact-site
        fraction phi = Z/mu1 (one site per unit, _bind_side_loss_release):
        the units whose group is gone, r*(1 - phi), are the modified
        monomer, and the intact sites leave with their units,
        dZ/dt -= r*phi*V_poly. No-op on any other pool."""
        cdef int row = self.slr_row[pool_i]
        cdef int z_slot, mod
        cdef double z, r_mod
        cdef double phi = 1.0
        if row < 0 or r_release == 0.0:
            return
        z_slot = self.sgh_z_slot[row]
        z = y[z_slot]
        if z < 0.0:
            z = 0.0
        if mu1_mol > 0.0:
            phi = z / mu1_mol
            if phi > 1.0:
                phi = 1.0
        r_mod = r_release * (1.0 - phi)
        mod = self.slr_mod_idx[pool_i]
        small_src[mono] = small_src.get(mono, 0.0) - r_mod
        small_src[mod] = small_src.get(mod, 0.0) + r_mod
        dz_dt[z_slot - self.num_core_species - self.num_qssa_u] -= r_release * phi * V_poly

    def _side_group_daughter(self, pool, channel_label, m_x, pool_index):
        """The configured X-loss feature daughter of a daughter-accounting
        carrier's channel (side_group_daughter_pool_label), whose per-chain
        defect must be the carrier's plus M_X and whose repeat unit mass the
        carrier's: the chain moves intact and its mass leaves as M_X."""
        label = side_group_daughter_pool_label(pool.label, channel_label)
        j = pool_index.get(label)
        if j is None:
            raise ValueError(
                f"Pool {pool.label}: side_group_homolysis channel "
                f"'{channel_label}' books its losses into the X-loss feature "
                f"daughter '{label}' (the carrier's chains are moved by "
                f"another kernel), but no such pool is configured; the "
                f"deck spawns it at setup (Polymer.generate_side_loss_"
                f"daughters).")
        dst = self.polymer_pools[j]
        want = float(pool.chain_mass_defect_g_mol) + m_x
        if not math.isclose(dst.chain_mass_defect_g_mol, want, rel_tol=1e-9):
            raise ValueError(
                f"Pool {dst.label}: chain_mass_defect_g_mol="
                f"{dst.chain_mass_defect_g_mol!r} is not the carrier "
                f"{pool.label}'s defect plus M_X = {want!r}; the X-loss "
                f"daughter would mint or destroy condensed mass.")
        if not math.isclose(dst.monomer_mw_g_mol, pool.monomer_mw_g_mol,
                            rel_tol=1e-9):
            raise ValueError(
                f"Pool {dst.label}: monomer_mw_g_mol={dst.monomer_mw_g_mol!r} "
                f"differs from its carrier {pool.label}'s "
                f"{pool.monomer_mw_g_mol!r}; the chain moves intact.")
        return dst

    def _assert_sgh_carrier_mu1_stationary(self, int pool_i, pool):
        """§1h + P1-B (SGH kernel-v2): refuse any SGH carrier pool whose mu1
        is NOT stationary. Inventory-depletion SGH seeds Z(0) = sites*mu1(0)
        ONCE and lets Z deplete on its own (dZ/dt = -R); the carrier's
        mu0/mu1/mu2 are UNCHANGED by the SGH kernel. If some OTHER kernel or
        resolved reaction row moves the carrier's mu1, the seeded Br inventory
        would strand (mu1 grows) or double-count (mu1 shrinks), so it is a
        hard build error.

        Two arms:
        - DECK-KERNEL arm: the carrier pool config also declares a mu1-moving
          kernel -- k_homolysis (dmu1 -= ...), k_depropagation (dmu1 -= ...),
          tail_kinetics (arbitrary dmu1), k_unzip (dmu1 -= k_unzip*mu0),
          radical_qssa_unzip (its depropagation drains mu1). k_scission is
          ALLOWED: verified mu1-PRESERVING (RHS touches only dmu0/dmu2, never
          dmu1, and routes nothing to another pool). k_unzip / radical_qssa_
          unzip are already mutually exclusive with SGH at the deck level;
          the refusal here is a defensive backstop.
        - RUNTIME-ROW-FLUX arm (P1-B, r104): a live resolved pool moment-flux
          archetype row (MIGRATION / SCISSION_FRAGMENT / DISCRETE_CHIP /
          VOLATILE_EJECTION, or the legacy UNRESOLVED mu1 transfer) whose
          source or destination pool is this carrier -- every such row writes
          the pool's mu_indices[1]. Deck-kernel enumeration alone is
          insufficient (a resolved apportionment/migration row can move mu1
          with no carrier-side kernel declared). Also catches an
          end-radical-daughter mu1 CREDIT if the carrier were ever spawned as
          a k_homolysis daughter.
        """
        label = pool.label
        # --- Deck-kernel arm (read pool config directly; robust to the
        # per-kernel flatten ordering) ---
        if pool.k_homolysis is not None:
            self._raise_sgh_mu1_mover(label, "k_homolysis")
        if pool.k_depropagation is not None:
            self._raise_sgh_mu1_mover(label, "k_depropagation")
        if pool.tail_kinetics is not None:
            self._raise_sgh_mu1_mover(label, "tail_kinetics")
        if float(pool.k_unzip or 0.0) > 0.0:
            self._raise_sgh_mu1_mover(label, "k_unzip")
        if pool.radical_qssa_unzip is not None:
            self._raise_sgh_mu1_mover(label, "radical_qssa_unzip")
        # P1-B (solver-build, covers BOTH deck-built and directly-built /
        # artifact-replay configs): an SGH carrier that ALSO carries explicit
        # DP is a silent mass-corruption hole. get_total_polymer_condensed_
        # mass_g counts the explicit-DP repeat-unit mass in mu1_mol, but the
        # SGH `lost` term and Z(0) use the tail-only mu1 -- so the Br inventory
        # is tracked only for the tail population while explicit-DP repeat
        # units are silently reported as still brominated. The deck refuses
        # this pairing, but a direct PolymerPoolConfig or an artifact replay
        # can still build it, so refuse it here at the solver-build layer.
        if pool.explicit_dp_to_species_index:
            raise ValueError(
                f"SGH-v2 carrier pool '{label}' also carries a non-empty "
                f"explicit_dp_to_species_index ("
                f"{dict(pool.explicit_dp_to_species_index)}). Inventory-"
                f"depletion side_group_homolysis (kernel-v2) seeds and "
                f"depletes the Z inventory on the TAIL mu1 basis "
                f"(Z(0) = sites*mu1, and the lost-X term uses tail-only mu1) "
                f"and cannot coexist with explicit-DP on the same carrier: "
                f"get_total_polymer_condensed_mass_g counts the explicit-DP "
                f"repeat-unit mass in mu1_mol, so the condensed-mass tally "
                f"would mix tail-only Z with explicit-augmented mu1 and "
                f"silently corrupt Br mass/speciation (explicit-DP repeat "
                f"units would be reported as still brominated). Move SGH to a "
                f"pool without explicit DP.")
        # k_scission is intentionally NOT listed: it is mu1-PRESERVING
        # (dmu0_dt += k_s*(mu1-mu0); dmu2_dt += k_s*(mu1-mu3)/3; mu1 untouched,
        # no cross-pool routing). Verified against the RHS.

        # Catalog pool kernels: every term moves whole chains (their mu1)
        # off or onto a pool, and none carries the carrier's Z with them --
        # whichever pool owns the term.
        for other in self.polymer_pools:
            for kernel, template, src, dst in catalog_kernel_connections(other):
                if label in (src, dst):
                    self._raise_sgh_mu1_mover(
                        label, f"{kernel} term {template!r} of pool "
                               f"{other.label!r} (moves chains {src!r} -> "
                               f"{dst!r})")

        mu1_idx = int(pool.mu_indices[1])

        # --- End-radical-daughter mu1 CREDIT arm: if any pool's k_homolysis
        # credits a fragment mu1 into THIS carrier's mu1 slot (only possible
        # if the carrier were spawned as a homolysis daughter). khom arrays
        # are flattened before this pass. ---
        prim_mu1 = getattr(self, "khom_prim_mu1", None)
        sec_mu1 = getattr(self, "khom_sec_mu1", None)
        khom_en = getattr(self, "khom_enabled", None)
        if khom_en is not None and prim_mu1 is not None and sec_mu1 is not None:
            for j in range(len(khom_en)):
                if not khom_en[j]:
                    continue
                if int(prim_mu1[j]) == mu1_idx or int(sec_mu1[j]) == mu1_idx:
                    self._raise_sgh_mu1_mover(
                        label,
                        f"k_homolysis end-radical fragment credit from pool "
                        f"'{self.polymer_pools[j].label}'")

        # --- Runtime-row-flux arm (P1-B): scan resolved flux archetypes.
        # reaction_flux_archetype / reaction_src_pool / reaction_dst_pool are
        # resolved (initialize_model) BEFORE validate_configuration calls this
        # flattener; if absent (e.g. a direct validate_configuration in a unit
        # test with no reaction rows) there is nothing to scan. ---
        arch = getattr(self, "reaction_flux_archetype", None)
        if arch is None:
            return
        src = getattr(self, "reaction_src_pool", None)
        dst = getattr(self, "reaction_dst_pool", None)
        stp = getattr(self, "species_to_pool_indices", None)
        proxy = getattr(self, "is_pool_proxy", None)
        ir = getattr(self, "reactant_indices", None)
        ip = getattr(self, "product_indices", None)
        refused = getattr(self, "reaction_refused", None)
        n_core = self.num_core_species
        for r in range(arch.shape[0]):
            # Refused / conduit-deferred ("stamp-but-keep") rows hit the RHS
            # `continue` guard before ANY reactant/product/moment write, so
            # they apply zero mu1 flux and cannot violate carrier-mu1
            # stationarity. Excluding them keeps this guard aligned with what
            # the RHS actually integrates (a non-refused row on the carrier
            # still raises below).
            if refused is not None and refused[r]:
                continue
            a = arch[r]
            if a == FLUX_NONE or a == FLUX_SAME_POOL:
                continue
            if a == FLUX_UNRESOLVED:
                # Legacy mu1-only transfer: the moved pools are the ROW's
                # proxy species' pools (src/dst stay -1). Map via the proxy
                # membership arrays.
                if (proxy is None or stp is None
                        or ir is None or ip is None):
                    continue
                touched = False
                for slot in range(3):
                    for tbl in (ir, ip):
                        sp = tbl[r, slot]
                        if (sp != -1 and sp < n_core and proxy[sp]
                                and int(stp[sp]) == pool_i):
                            touched = True
                            break
                    if touched:
                        break
                if touched:
                    self._raise_sgh_mu1_mover(
                        label, f"live FLUX_UNRESOLVED mu1 transfer (row {r})")
                continue
            # MIGRATION / SCISSION_FRAGMENT / DISCRETE_CHIP / VOLATILE_EJECTION
            # all write a src (and possibly dst) pool's mu_indices[1].
            if src is not None and int(src[r]) == pool_i:
                self._raise_sgh_mu1_mover(
                    label,
                    f"live resolved row flux (archetype={int(a)}, row {r}, "
                    f"source pool)")
            if dst is not None and int(dst[r]) == pool_i:
                self._raise_sgh_mu1_mover(
                    label,
                    f"live resolved row flux (archetype={int(a)}, row {r}, "
                    f"destination pool)")

    def _raise_sgh_mu1_mover(self, label, competitor):
        raise ValueError(
            f"SGH-v2 carrier pool '{label}' requires a stationary mu1, but "
            f"it also carries a competing mu1-moving kernel/row: "
            f"{competitor}. Inventory-depletion side_group_homolysis "
            f"(kernel-v2) seeds the removable-X inventory once as "
            f"Z(0) = sites*mu1(0) and depletes it independently of the "
            f"carrier moments; a moving carrier mu1 would strand or "
            f"double-count the Br inventory. Remove the competing "
            f"kernel/flux from the SGH carrier (or move SGH to a stationary "
            f"pool).")

    def _flatten_radical_qssa_state(self):
        """Flatten each pool's validated+normalized radical_qssa_unzip channel
        into solver-owned per-pool flat numpy arrays (review round 21,
        finding 1 -- the M2-prep half; finding 2's mutation-bypass closure).

        Contract:
        - Populated ONLY here, on the initialize_model -> validate_configuration
          path, from the normalized dict validate_configuration just stored.
        - The M2 QSSA rate law will read ONLY these arrays, NEVER the dict.
          Once flattening has run, mutating the (mutable) nested dict behind
          the frozen PolymerPoolConfig is harmless to solver dynamics: the
          arrays are the solver's sole source of truth for the channel.
        - M1: NOTHING in the residual/RHS reads these yet (pinned by the
          bitwise residual-equality test); they are inert state.

        Layout (index = pool position in self.polymer_pools):
        - qssa_enabled[i] (int8): 1 iff the pool has a channel configured.
        - qssa_ki_A/n/Ea[i]: initiation Arrhenius triplet (A [s^-1]).
        - qssa_kdp_A/n/Ea[i]: depropagation triplet (A [s^-1]).
        - qssa_kt_A/n/Ea[i]: termination triplet (A [m^3 mol^-1 s^-1]; the
          only bimolecular block).
        - qssa_efficiency[i], qssa_monomer_yield[i]: f_i, y_m in (0, 1]
          (default 1.0 on channel-absent pools -- the inert value).
        - qssa_has_transfer[i] (int8) + qssa_ktr_A/n/Ea[i]: optional transfer
          triplet (A [s^-1], pseudo-first-order: ktr multiplies R directly
          in the M2 rate law); zeros when absent.
        Disabled rows keep zero Arrhenius slots: every M2 consumer must gate
        on qssa_enabled, never on A != 0.

        Weak-link U-state extension (milestones ii+iii):
        - qssa_weaklink[i] (int8): 1 iff the channel carries the weak-link
          vocabulary. On weak-link rows the legacy summed qssa_kt_* slots
          STAY ZERO (there is no summed 'termination' block); the RHS must
          gate on qssa_weaklink, never on kt_A != 0.
        - qssa_kia_A/n/Ea[i]: initiation_allyl triplet (A [s^-1]).
        - qssa_ktrec_A/n/Ea[i] / qssa_ktdisp_A/n/Ea[i]: split termination
          triplets (A [m^3 mol^-1 s^-1]); the RHS uses kt_total =
          kt_rec(T) + kt_disp(T) wherever the legacy summed kt appeared,
          and ONLY kt_disp sources dU/dt.
        - qssa_u0[i]: unsaturated_tail_ends_initial [mol] (state IC, set
          into y0 by set_initial_conditions after the census trap there).
        - qssa_u_slot[i] (int32): absolute ODE index of the pool's U slot
          (n_core + running weak-link ordinal, matching the count
          initiate_tolerances used to extend neq), or -1.
        """
        n_pools = len(self.polymer_pools)
        self.qssa_enabled = np.zeros(n_pools, dtype=np.int8)
        self.qssa_ki_A = np.zeros(n_pools, dtype=float)
        self.qssa_ki_n = np.zeros(n_pools, dtype=float)
        self.qssa_ki_Ea = np.zeros(n_pools, dtype=float)
        self.qssa_kdp_A = np.zeros(n_pools, dtype=float)
        self.qssa_kdp_n = np.zeros(n_pools, dtype=float)
        self.qssa_kdp_Ea = np.zeros(n_pools, dtype=float)
        self.qssa_kt_A = np.zeros(n_pools, dtype=float)
        self.qssa_kt_n = np.zeros(n_pools, dtype=float)
        self.qssa_kt_Ea = np.zeros(n_pools, dtype=float)
        self.qssa_efficiency = np.ones(n_pools, dtype=float)
        self.qssa_monomer_yield = np.ones(n_pools, dtype=float)
        self.qssa_has_transfer = np.zeros(n_pools, dtype=np.int8)
        self.qssa_ktr_A = np.zeros(n_pools, dtype=float)
        self.qssa_ktr_n = np.zeros(n_pools, dtype=float)
        self.qssa_ktr_Ea = np.zeros(n_pools, dtype=float)
        self.qssa_weaklink = np.zeros(n_pools, dtype=np.int8)
        self.qssa_kia_A = np.zeros(n_pools, dtype=float)
        self.qssa_kia_n = np.zeros(n_pools, dtype=float)
        self.qssa_kia_Ea = np.zeros(n_pools, dtype=float)
        self.qssa_ktrec_A = np.zeros(n_pools, dtype=float)
        self.qssa_ktrec_n = np.zeros(n_pools, dtype=float)
        self.qssa_ktrec_Ea = np.zeros(n_pools, dtype=float)
        self.qssa_ktdisp_A = np.zeros(n_pools, dtype=float)
        self.qssa_ktdisp_n = np.zeros(n_pools, dtype=float)
        self.qssa_ktdisp_Ea = np.zeros(n_pools, dtype=float)
        self.qssa_u0 = np.zeros(n_pools, dtype=float)
        self.qssa_u_slot = np.full(n_pools, -1, dtype=np.int32)

        u_count = 0
        for i, pool in enumerate(self.polymer_pools):
            q = pool.radical_qssa_unzip
            if q is None:
                continue
            self.qssa_enabled[i] = 1
            self.qssa_ki_A[i] = q["initiation"]["A"]
            self.qssa_ki_n[i] = q["initiation"]["n"]
            self.qssa_ki_Ea[i] = q["initiation"]["Ea"]
            self.qssa_kdp_A[i] = q["depropagation"]["A"]
            self.qssa_kdp_n[i] = q["depropagation"]["n"]
            self.qssa_kdp_Ea[i] = q["depropagation"]["Ea"]
            if "initiation_allyl" in q:
                # Weak-link channel: split termination blocks (validator
                # invariant: no legacy summed 'termination' key exists on
                # this channel, so qssa_kt_* stays 0 on this row).
                self.qssa_weaklink[i] = 1
                self.qssa_kia_A[i] = q["initiation_allyl"]["A"]
                self.qssa_kia_n[i] = q["initiation_allyl"]["n"]
                self.qssa_kia_Ea[i] = q["initiation_allyl"]["Ea"]
                self.qssa_ktrec_A[i] = q["termination_recombination"]["A"]
                self.qssa_ktrec_n[i] = q["termination_recombination"]["n"]
                self.qssa_ktrec_Ea[i] = q["termination_recombination"]["Ea"]
                self.qssa_ktdisp_A[i] = q["termination_disproportionation"]["A"]
                self.qssa_ktdisp_n[i] = q["termination_disproportionation"]["n"]
                self.qssa_ktdisp_Ea[i] = q["termination_disproportionation"]["Ea"]
                self.qssa_u0[i] = q["unsaturated_tail_ends_initial"]
                self.qssa_u_slot[i] = self.num_core_species + u_count
                u_count += 1
            else:
                self.qssa_kt_A[i] = q["termination"]["A"]
                self.qssa_kt_n[i] = q["termination"]["n"]
                self.qssa_kt_Ea[i] = q["termination"]["Ea"]
            self.qssa_efficiency[i] = q["efficiency"]
            self.qssa_monomer_yield[i] = q["monomer_yield"]
            if q["transfer"] is not None:
                self.qssa_has_transfer[i] = 1
                self.qssa_ktr_A[i] = q["transfer"]["A"]
                self.qssa_ktr_n[i] = q["transfer"]["n"]
                self.qssa_ktr_Ea[i] = q["transfer"]["Ea"]

        # Layout cross-check: initiate_tolerances extended neq by ITS count
        # of weak-link pools (key presence on the raw configs, before
        # validation); the flattening above counted the normalized ones. A
        # mismatch means the pool configs changed between the two passes --
        # the U slots would alias core species. Fail loud, never integrate
        # a misaligned state vector.
        if u_count != self.num_qssa_u:
            raise ValueError(
                f"Weak-link U-state layout drift: initiate_tolerances "
                f"allocated {self.num_qssa_u} U slot(s) but "
                f"validate_configuration found {u_count} weak-link pool(s). "
                f"The pool configs changed between ODE layout and "
                f"validation; rebuild the reactor (initialize_model) instead "
                f"of mutating pool configs in place.")

        # QSSA lineages (I-260): an X-loss daughter's chains are its
        # parent's melt, so their radicals terminate with the parent's and
        # the two share one radical population (qssa_root: the lineage's
        # parent pool, qssa_lineage_n: its member count, qssa_lineage_members:
        # the pools of lineages with more than one member). The daughter
        # carries a copy of the parent's channel; anything else is refused,
        # and so is a weak-link channel, whose U state is per pool.
        self.qssa_root = np.arange(n_pools, dtype=np.int32)
        self.qssa_lineage_n = np.ones(n_pools, dtype=np.int32)
        pool_index = {p.label: j for j, p in enumerate(self.polymer_pools)}
        for i, pool in enumerate(self.polymer_pools):
            parent = getattr(pool, "side_loss_parent", None)
            if not (parent and self.qssa_enabled[i]):
                continue
            j = pool_index.get(parent)
            if j is None or not self.qssa_enabled[j]:
                continue
            if self.qssa_weaklink[i] or self.qssa_weaklink[j]:
                raise ValueError(
                    f"Pool {pool.label}: an X-loss daughter of {parent} "
                    f"carrying a weak-link radical_qssa_unzip channel; the "
                    f"U state is per pool and has no shared form.")
            if pool.radical_qssa_unzip != self.polymer_pools[j].radical_qssa_unzip:
                raise ValueError(
                    f"Pool {pool.label}: its radical_qssa_unzip channel "
                    f"differs from its parent {parent}'s; one radical "
                    f"population has one set of constants.")
            self.qssa_root[i] = j
            self.qssa_lineage_n[j] += 1
        self.qssa_lineage_members = [
            i for i in range(n_pools)
            if self.qssa_enabled[i] and self.qssa_lineage_n[self.qssa_root[i]] > 1]

    def _reference_state_tripwire(self, core_species, core_reactions):
        """Build-time thermo reference-state tripwire (spec 2026-06-11).

        Invariant (docs/multi_pool_design.md §5.2): per species, per reaction
        side, the thermo reference state must match the phase residence.
        Today every species carries the gas reference UNIFORMLY, so
        boundary-crossing reactions between same-mass chains cancel exactly
        -- and the cancellation is STRUCTURAL (HBI saturates the radical onto
        the same GAV parent the proxy runs on). The danger is the PARTIAL
        fix: a melt correction applied to the condensed set alone injects up
        to ~11.6 decades into every boundary-crossing Keq. This pass
        measures, per reversible CORE reaction, the unpaired reference-state
        magnitude U over the physically-melt participants (condensed
        gas_species_mask OR is_polymer_proxy-tagged, unvetoed AND dual-axis
        polymer-sized -- spec §5.1 C3-amended, r89 dual-axis), logs a census
        above
        REFERENCE_STATE_CENSUS_DECADES, refuses above
        REFERENCE_STATE_REFUSE_DECADES (unless
        allow_unpaired_reference_state -- the census still logs), and warns
        on mixed library-vs-GAV thermo provenance among chain-scale
        counterparties (the spec-§5.3 decoupling fingerprint that U, being
        MW-only, is structurally blind to). Cost: ~|core reactions|
        Sackur-Tetrode evaluations per rebuild -- nil.
        """
        cdef int n_core = self.num_core_species
        T = self.T.value_si
        self.reference_state_census = []
        self.reference_state_max_decades = 0.0
        self.reference_state_declared_convention = []

        mask = self.gas_species_mask

        # Chain-scale window (spec §5.3): largest configured pool monomer MW
        # + slack, in kg/mol. Since the r89 dual-axis amendment this window
        # serves ONLY the provenance counterparty predicate (i) below -- the
        # physically-melt class gate (tag branch) uses the dual-axis
        # polymer-sized test instead (the window form melt-classified DP-2
        # gas volatiles: PP run-9 1,5-hexadiene, 82.15 g/mol > propene
        # window 52.1 but only 2.0 monomer-equivalents). Predicate (ii)
        # (sharing the proxy's saturated parent) was REJECTED as the
        # cheap-at-init test: it needs graph saturation + isomorphism per
        # gas species per rebuild; (i) is a float compare on data already in
        # hand. FR1-K2 mass audit note (round-72 P2): deliberately NOT
        # defect-aware. This is a PER-REPEAT-UNIT MW scale window, not a
        # condensed-mass computation; chain_mass_defect_g_mol is a
        # PER-CHAIN correction (one lost X per chain), and X-loss feature
        # pools keep the parent's monomer_mw by design (intact-chain
        # transfer), so the per-unit window is exactly right as is.
        chain_window_kg = (max((float(p.monomer_mw_g_mol) for p in self.polymer_pools),
                               default=0.0)
                           + REFERENCE_STATE_MW_SLACK_G_MOL) / 1000.0

        # r89 dual-axis gate inputs. The threshold is imported at call time
        # (documented solver <-> polymer module cycle, same posture as
        # rmgpy.polymer's function-local solver imports) so the gate and
        # rmgpy.polymer._discrete_is_polymer_sized share ONE constant; the
        # per-pool axes come from the flattened configs (monomer_mw_g_mol +
        # monomer_heavy_atoms), plumbed producer-side from the monomer
        # STRUCTURE -- the heavy axis is never approximated from mass.
        from rmgpy.polymer import _IMPOSTOR_DISCRETE_MONOMER_UNITS as _units
        pool_axes = [(float(getattr(p, "monomer_mw_g_mol", 0.0) or 0.0),
                      int(getattr(p, "monomer_heavy_atoms", 0) or 0))
                     for p in self.polymer_pools]

        is_melt = [False] * n_core
        mws = [0.0] * n_core
        heavies = [-1] * n_core   # -1 = unknown (no structure, no stamp)
        # Backstop census: polymer-sized (dual-axis) members whose melt
        # classification the durable gas veto SUPPRESSED. create_reacted_copy
        # returns None both for a genuine gas volatile and for a wing-match
        # FAILURE on a real chain-scale fragment; the veto trusts that None as
        # authoritative "gas", so a mis-vetoed genuine chain would be silently
        # dropped from the melt sum instead of loudly refused. Recording +
        # warning here keeps that case visible for a human without regressing
        # the alpha-methylstyrene fix (the build no longer refuses). r89: the
        # census keys on the SAME dual-axis chain-scale notion as the gate --
        # a genuinely-volatile vetoed species (AMS, 1.13 monomer-equivalents)
        # is no longer announced; only chain-scale suppressions are.
        veto_suppressed_chain_scale = []
        # r89 undecidable-axis census (mirror of rmgpy.polymer
        # ._warn_impostor_axis_undecidable): a tag-branch candidate whose
        # polymer-sized verdict could not be established on BOTH axes is
        # NEVER classified blind -- it degrades to conservative-gas -- but
        # the degradation is announced, per species, so a silent axis-data
        # gap cannot quietly unpair a genuine chain.
        axis_undecidable = []
        for i in range(n_core):
            spc = core_species[i]
            # Input contract: consumer-world species carry no structure
            # (molecule == []); MW arrives via the species-level
            # molecular_weight quantity, populated by the runner from the
            # artifact's composition. Prefer it (value_si is kg/molecule);
            # the structure branch below is a defensive fallback (normally
            # unreachable: the molecular_weight property lazily derives from
            # the structure when one exists -- species.py:278-282 -- so the
            # quantity path above short-circuits it).
            mol_list = getattr(spc, "molecule", None)
            mw_q = getattr(spc, "molecular_weight", None)
            if mw_q is not None and mw_q.value_si > 0.0:
                mws[i] = mw_q.value_si * constants.Na
            elif mol_list:
                mws[i] = mol_list[0].get_molecular_weight()
            else:
                mws[i] = 0.0
            # Heavy-atom (non-H) count for the r89 dual-axis gate: from the
            # structure when one exists; otherwise from the species-level
            # props stamp the consumer runner populates from the artifact's
            # composition (real composition data both ways -- the r89
            # data-flow constraint forbids deriving this axis from mass).
            # -1 = unknown: the heavy axis is then uncomputable and the
            # dual-axis verdict degrades to conservative-gas + census.
            _props = getattr(spc, "props", None)
            if mol_list:
                _m0 = mol_list[0]
                heavies[i] = _m0.get_num_atoms() - _m0.get_num_atoms('H')
            elif isinstance(_props, dict):
                # HARDCODED literal must equal rmgpy.polymer
                # .POLYMER_HEAVY_ATOM_COUNT_KEY -- pinned by
                # test_heavy_atom_key_literal_matches_solver_gate
                # (rename both).
                _hv = _props.get("polymer_heavy_atom_count")
                if _hv is not None:
                    heavies[i] = int(_hv)
            # Physically-melt CLASS (spec §5.1 C3-AMENDED, r89 dual-axis):
            # the condensed branch (pool-configured by input)
            # unconditionally; the tag branch only when DUAL-AXIS
            # POLYMER-SIZED (mass AND heavy axes both computable and both
            # >= _IMPOSTOR_DISCRETE_MONOMER_UNITS monomer-equivalents
            # against >= 1 pool -- _dual_axis_polymer_sized, the exact
            # mirror of rmgpy.polymer._discrete_is_polymer_sized). The size
            # conjunct is part of the class DEFINITION, not bolted onto the
            # provenance set: family.py:1657 blanket-tags every structure
            # of a proxy-touching reaction (including H2), so a raw
            # tag-read would be correct only by spawn-pass ordering -- the
            # structural gate cannot be broken by a lifecycle reordering
            # the way a tag-read can. A stale tag on a below-threshold
            # species simply FAILS the conjunct and is excluded: expected
            # and silent (its gas reference state is CORRECT -- PP run-9's
            # 1,5-hexadiene, a genuine DP-2 gas volatile at 2.0
            # monomer-equivalents, is exactly this class; the pre-r89 MW
            # window wrongly melt-classified it and refused an all-gas
            # allyl recombination at U ~ 11 decades).
            #
            # DURABLE GAS-VOLATILE VETO (rmgpy.polymer): the tag branch keys off
            # is_polymer_proxy, which is a monotonic multi-writer sticky cache
            # with no authoritative "gas" clear point -- a genuine discrete gas
            # volatile (e.g. alpha-methylstyrene) that got proxy-contaminated
            # and happens to clear the MW window would be wrongly counted as a
            # melt participant (the reference-state-tripwire leak). The polymer
            # handshake / chip demotion stamps such a species with a POSITIVE,
            # durable veto in props (POLYMER_REFERENCE_STATE_GAS_VETO_KEY) that
            # the proxy stamping machinery never touches; honor it here. C3
            # mask-lagged chains never receive the veto, so they stay melt.
            _gas_veto = False
            if isinstance(_props, dict):
                # HARDCODED literal must equal rmgpy.polymer
                # .POLYMER_REFERENCE_STATE_GAS_VETO_KEY -- pinned by
                # test_gas_veto_key_literal_matches_solver_gate (rename both).
                _gas_veto = bool(_props.get("polymer_reference_state_gas_veto", False))
            _proxy_i = bool(getattr(spc, "is_polymer_proxy", False))
            # r89 gate order: veto precedence is unchanged (a vetoed species
            # NEVER classifies melt via the tag branch, whatever its size);
            # the dual-axis verdict is evaluated for every gas-masked
            # tagged candidate because both censuses key on it.
            _sized = False
            _missing = None
            if mask[i] and _proxy_i:
                _sized, _missing = _dual_axis_polymer_sized(
                    mws[i] * 1000.0, heavies[i], pool_axes, _units)
                if not _sized and _missing is not None:
                    axis_undecidable.append(
                        (getattr(spc, "label", "?"), _missing))
            is_melt[i] = ((not mask[i])
                          or (_proxy_i and not _gas_veto and _sized))
            # The veto is only a backstop concern when it suppressed an
            # otherwise-melt POLYMER-SIZED tag-branch member (gas-masked,
            # proxy, dual-axis sized). A below-threshold vetoed species
            # would fail the size conjunct anyway -- its exclusion is
            # correct and silent (r89: alpha-methylstyrene, 1.13
            # monomer-equivalents, is no longer announced).
            if _gas_veto and _proxy_i and mask[i] and _sized:
                veto_suppressed_chain_scale.append(
                    (getattr(spc, "label", "?"), mws[i]))

        # Preserve the reference-state backstop for the ambiguous handshake
        # false-None case (code-review IMPORTANT #1): a genuine chain-scale
        # fragment the wing-matcher missed is now gas-vetoed rather than
        # refused, so surface it for a human. NOT a refusal -- census only.
        self.gas_veto_census = list(veto_suppressed_chain_scale)
        if veto_suppressed_chain_scale:
            logging.warning(
                "THERMO REFERENCE-STATE GAS VETO: %d chain-scale "
                "(dual-axis polymer-sized at >= %.1f monomer-equivalents "
                "of mass AND heavy atoms) product(s) were excluded from the "
                "melt reference-state by the durable gas-volatile veto: %s. "
                "If any is actually a polymer chain, it is a handshake "
                "create_reacted_copy false-None (wing-match failure), NOT a "
                "thermo problem -- investigate the handshake, do not touch "
                "reference states.",
                len(veto_suppressed_chain_scale), _units,
                ", ".join("%s (%.1f g/mol)" % (lbl, mw * 1000.0)
                          for lbl, mw in veto_suppressed_chain_scale))

        # r89 undecidable-axis census: conservative-gas degradations are
        # announced, never silent (mirror of rmgpy.polymer
        # ._warn_impostor_axis_undecidable -- never classify blind; the
        # solver tripwire's condensed branch and the r71 hard-fail stay the
        # loud backstops).
        self.reference_state_axis_undecidable = list(axis_undecidable)
        if axis_undecidable:
            logging.warning(
                "THERMO REFERENCE-STATE AXIS UNDECIDABLE (r89 census): %d "
                "proxy-tagged gas-masked candidate(s) had an uncomputable "
                "polymer-sized axis, so the dual-axis melt conjunct could "
                "not be established on BOTH axes; classifying "
                "conservative-gas (never melt-classify blind): %s. If one "
                "is a genuine chain, plumb its missing axis data (species "
                "composition / pool monomer structure), do not touch "
                "reference states.",
                len(axis_undecidable),
                ", ".join("%s (missing: %s)" % (lbl, miss)
                          for lbl, miss in axis_undecidable))

        offenders = []
        for rxn in core_reactions:
            if not getattr(rxn, "reversible", False):
                continue
            if getattr(rxn, "polymer_refused", False):
                # Refused rows (item-18 stamp-but-keep + the PP v1
                # gas-association refusal, adjudicated round 63) have their
                # WHOLE flux suppressed in the residual (reaction_refused),
                # so their reference-state pairing is moot: no Keq of theirs
                # ever moves a mole. Skip them here -- keyed on the refused
                # stamp ONLY, no other relaxation (the SAME row unrefused
                # still trips; pinned by
                # test_refused_row_with_unpaired_reference_state_passes_initialize).
                continue
            j = self.reaction_index[rxn]
            r_idx = [int(k) for k in self.reactant_indices[j] if 0 <= k < n_core]
            p_idx = [int(k) for k in self.product_indices[j] if 0 <= k < n_core]
            melt_r = [k for k in r_idx if is_melt[k]]
            melt_p = [k for k in p_idx if is_melt[k]]
            if not melt_r and not melt_p:
                continue  # all-gas reaction: gas reference uniformly correct

            # Leak self-assertion (spec §5.1 C3 amendment, r89 dual-axis):
            # cannot-happen guard on every species entering the melt sum. A
            # gas-classified member can only be here via the tag branch,
            # whose dual-axis polymer-sized conjunct is enforced in the gate
            # above; if a below-threshold (or axis-undecidable) species ever
            # reaches this point (a future refactor recomputing membership
            # without the conjunct), raise the CLASSIFICATION error loudly
            # instead of computing a large U and misattributing it to thermo.
            for k in melt_r + melt_p:
                # Log-domain guard (input contract): a melt participant with
                # no resolvable molar mass would send mw <= 0 into
                # _sackur_tetrode_s_trans's math.log. A melt chain with no
                # structure is the tripwire's MAIN case (consumer world), not
                # an edge case -- never skip it silently, and never let a raw
                # 'math domain error' escape.
                if mws[k] <= 0.0:
                    raise ValueError(
                        "THERMO REFERENCE-STATE TRIPWIRE (input contract): no "
                        f"molar mass available for melt participant "
                        f"'{core_species[k].label}' (molecule list empty and "
                        "species-level molecular_weight unset); consumer-world "
                        "species must carry molecular_weight populated from "
                        "the artifact's composition/mw fields -- this is an "
                        "input-contract violation, NOT a thermo problem; do "
                        "not respond by touching reference states.")
                _assert_chain_scale_melt_member(
                    core_species[k].label, mws[k], heavies[k], bool(mask[k]),
                    pool_axes, _units)

            u = _unpaired_reference_decades(
                [mws[k] for k in melt_r], [mws[k] for k in melt_p], T)
            if u > self.reference_state_max_decades:
                self.reference_state_max_decades = u
            if u > REFERENCE_STATE_CENSUS_DECADES:
                self.reference_state_census.append((str(rxn), u))
                logging.warning(
                    "THERMO REFERENCE-STATE CENSUS: reaction %s: U = %.2f "
                    "decades at T = %.1f K (above the %.1f-decade census "
                    "bound; see the invariant section of "
                    "docs/multi_pool_design.md).",
                    rxn, u, T, REFERENCE_STATE_CENSUS_DECADES)
            if u > REFERENCE_STATE_REFUSE_DECADES:
                # Declared single convention (polymer_phase residence): every
                # participant condensed in the mask means every participant
                # takes the gas reference at its melt concentration, the
                # convention the deck declared -- nothing is unpaired against
                # another phase, so the row is listed, not refused. A row
                # with a gas-masked participant is checked as before.
                if self.declared_melt_convention and not any(mask[k] for k in r_idx + p_idx):
                    self.reference_state_declared_convention.append((rxn, u))
                else:
                    offenders.append((rxn, u))

            # Provenance sensor (spec §5.3): counterparty set = melt
            # participants + gas participants inside the shared chain-scale
            # MW window of a melt participant. Deliberately NARROW: small
            # co-reactants (H2, CH4, every abstraction partner) legitimately
            # take library thermo while chains take GAV -- sweeping them in
            # would warn on every healthy deck (alarm fatigue re-arming the
            # exact landmine this sensor guards).
            melt_set = set(melt_r) | set(melt_p)
            counterparties = set(melt_set)
            for k in set(r_idx) | set(p_idx):
                if k in melt_set:
                    continue
                if any(abs(mws[k] - mws[m]) <= chain_window_kg for m in melt_set):
                    counterparties.add(k)
            examples = {}
            for k in sorted(counterparties):
                prov = _thermo_provenance(core_species[k])
                if prov is not None and prov not in examples:
                    examples[prov] = core_species[k].label
            if "library" in examples and "gav" in examples:
                logging.warning(
                    "THERMO REFERENCE-STATE PROVENANCE: mixed thermo "
                    "provenance among chain-scale counterparties in reaction "
                    "%s (library: %s; group additivity: %s) -- the structural "
                    "cancellation that keeps gas-reference thermo safe on "
                    "melt chains may be broken for this pair; see the "
                    "invariant section of docs/multi_pool_design.md.",
                    rxn, examples["library"], examples["gav"])

        if offenders:
            if self.allow_unpaired_reference_state:
                logging.warning(
                    "THERMO REFERENCE-STATE TRIPWIRE: "
                    "allow_unpaired_reference_state=True: bypassing the "
                    "reference-state refusal for %d reaction(s); max U = "
                    "%.2f decades. The deck author asserts the thermo "
                    "handles the melt reference state.",
                    len(offenders), self.reference_state_max_decades)
            else:
                lines = "\n".join(
                    f"  {rxn}: U = {u:.2f} decades" for rxn, u in offenders)
                raise ValueError(
                    "THERMO REFERENCE-STATE TRIPWIRE: this deck has "
                    f"{len(offenders)} reversible core reaction(s) with an "
                    "unpaired reference-state term (U > "
                    f"{REFERENCE_STATE_REFUSE_DECADES} decades) at "
                    f"T = {T:.1f} K:\n{lines}\n"
                    "Do NOT fix this by applying a melt correction to the "
                    "condensed set alone -- that injects the full mismatch "
                    "into every boundary-crossing Keq (decapping kb by up to "
                    "~12 orders of magnitude). See the thermo reference-state "
                    "invariant section of docs/multi_pool_design.md. If the "
                    "deck's thermo genuinely handles the melt reference "
                    "state, set allow_unpaired_reference_state=True on the "
                    "reactor.")

    def get_const_spc_indices(self, core_species):
        """Identify indices of constant species."""
        if not self.const_spc_names:
            return
        self.const_spc_indices = list()
        for name in self.const_spc_names:
            for spc in core_species:
                if spc.label == name:
                    self.const_spc_indices.append(core_species.index(spc))
                    break

    def convert_initial_keys_to_species_objects(self, species_dict):
        """Convert initial mole fraction keys to species objects."""
        initial_mole_fractions = dict()
        for label, moleFrac in self.initial_mole_fractions.items():
            initial_mole_fractions[species_dict[label]] = moleFrac
        self.initial_mole_fractions = initial_mole_fractions

        conditions = dict()
        if self.sens_conditions is not None:
            for label, value in self.sens_conditions.items():
                if label in ("T", "P"):
                    conditions[label] = value
                else:
                    conditions[species_dict[label]] = value
        self.sens_conditions = conditions

    def _apply_pool_phase_overrides(self, mask_arr, species_list,
                                    record_indices):
        """Stage 2 of the two-stage phase classifier (item 17, spec
        2026-06-12 SS3(a)): the per-pool config-label override pass, factored
        so ONE code path serves both masks. Run with
        (gas_species_mask, core_species, record_indices=True) -- behavior
        bit-identical to the historical inline pass, including the index
        bookkeeping (species_to_pool_indices / is_pool_proxy /
        pool_mu1_indices / pool_mu0_indices and the mu1 fallback) -- and
        with (prospective_gas_mask, chain(core, edge),
        record_indices=False), where ONLY mask writes happen (the index
        arrays are core-sized; an edge match must never write them)."""
        # ONE base-label convention (function-local import: documented
        # solver<->polymer cycle note at the module head): trailing-index
        # strip ONLY -- the former first-'(' truncation could alias a
        # SMILES-labelled species (structural parentheses) onto a pool label.
        from rmgpy.polymer import strip_rmg_index_suffix
        n_core = self.num_core_species
        n_total = len(species_list)
        for pool_i, pool in enumerate(self.polymer_pools):
            for i in range(n_total):
                base_label = strip_rmg_index_suffix(species_list[i].label)
                if base_label == pool.label:
                    mask_arr[i] = False
                    if record_indices:
                        self.species_to_pool_indices[i] = pool_i
                        self.is_pool_proxy[i] = 1
                    break
            mu1_target_label = f"{pool.label}_mu1"
            for i in range(n_total):
                # Handle RMG renaming: "PS_mu1(2)" -> "PS_mu1"
                base_label = strip_rmg_index_suffix(species_list[i].label)
                if base_label == mu1_target_label:
                    if record_indices:
                        self.pool_mu1_indices[pool_i] = i
                    mask_arr[i] = False  # Ensure it's poly phase

            mu0_target_label = f"{pool.label}_mu0"
            for i in range(n_total):
                base_label = strip_rmg_index_suffix(species_list[i].label)
                if base_label == mu0_target_label:
                    if record_indices:
                        self.pool_mu0_indices[pool_i] = i
                    mask_arr[i] = False

            # map explicit oligomers
            for dp, idx in pool.explicit_dp_to_species_index.items():
                if 0 <= idx < n_core:
                    if record_indices:
                        self.species_to_pool_indices[idx] = pool_i
                    mask_arr[idx] = False

            # map the moment indices
            for idx in pool.mu_indices:
                if 0 <= idx < n_core:
                    if record_indices:
                        self.species_to_pool_indices[idx] = pool_i
                    mask_arr[idx] = False

            # monomer_poly_index (the unzip/QSSA release target) is
            # deliberately NOT mapped or condensed here (incident 2026-07-03,
            # design B-prime): the release target is the mechanism's gas
            # volatile (e.g. styrene for PS), and force-condensing it made
            # every reversible gas core reaction producing it carry an
            # ~11-decade unpaired reference-state term (correctly refused by
            # the tripwire). The release deposits into the GAS amount basis
            # (small_src -> dn_dt, mol/s); validate_configuration now REQUIRES
            # the monomer index to be masked GAS.

            if record_indices and self.pool_mu1_indices[pool_i] == -1:
                # Fallback: Try to use the config index if label lookup failed
                # (This catches cases where renaming didn't happen as expected)
                cfg_idx = pool.mu_indices[1]
                if 0 <= cfg_idx < n_core:
                    self.pool_mu1_indices[pool_i] = cfg_idx
                else:
                    print(f"WARNING: Could not locate mu1 species for pool {pool.label}. Polymer chemistry will be inert.")

    def initialize_model(self, core_species, core_reactions, edge_species, edge_reactions,
                      surface_species=None, surface_reactions=None, pdep_networks=None,
                      atol=1e-16, rtol=1e-8, sensitivity=False, sens_atol=1e-6, sens_rtol=1e-4,
                      filter_reactions=False, conditions=None, **kwargs):
        """
        Initialize the polymer hybrid reactor model.
        """
        ReactionSystem.initialize_model(
            self, core_species, core_reactions, edge_species, edge_reactions,
            surface_species, surface_reactions, pdep_networks,
            atol, rtol, sensitivity, sens_atol, sens_rtol, filter_reactions, conditions
        )

        # Round-27 P1-C: engine rebuild signature for the sidecar
        # stale_topology predicate -- stable (label, index) identity keys of
        # the core this engine was ACTUALLY built against. save_everything
        # compares it against the current model core at sidecar write; any
        # mismatch (including a same-count species/reaction swap the old
        # count predicate missed) marks the emission stale in-band.
        # Function-local import: avoids a solver->polymer module cycle.
        from rmgpy.polymer import core_topology_signature
        self.core_topology_signature = core_topology_signature(
            core_species, core_reactions)

        cdef int n_core = self.num_core_species
        cdef int n_rxn = len(core_reactions) + len(edge_reactions)
        # Flag end-group (terminal) reactions so their proxy rate scales by
        # chain-end density (mu0) instead of monomer-unit density (mu1). Iterate
        # in the SAME order that builds kf/kb and the ir/ip arrays
        # (generate_rate_coefficients below) so the index matches r_idx.
        self.is_end_group_reaction = np.zeros(n_rxn, dtype=np.int8)
        for i, rxn in enumerate(itertools.chain(core_reactions, edge_reactions)):
            if getattr(rxn, "is_end_group_reaction", False):
                self.is_end_group_reaction[i] = 1

        # r71 FIX 2 (PP run-5 stall) -- rebuild RESTAMPING at the last honest
        # chokepoint. The generation-time refusal stamps are ad-hoc object
        # attributes that have PROVEN losable between stamping and this
        # rebuild (canonical dedup in check_for_existing_reaction; the
        # __reduce__ serialization contract -- rmgpy/reaction.py documents
        # "the solver re-derives them at initialize_model", and THIS is that
        # re-derivation). Idempotently re-run the shape predicate over every
        # row (core AND edge) BEFORE reaction_refused is captured and BEFORE
        # the thermo reference-state tripwire scans core rows: this loop sees
        # the final post-promotion/post-dedup list regardless of construction
        # path. stamp_gas_association_refusal early-returns on already-refused
        # rows, so an upstream qssa-invalid census reason is never overwritten
        # (precedence preserved) and stamped rows are unchanged. The other
        # refusal sibling (stamp_polymer_flux_archetype's UNRESOLVED
        # FEATURE-radical branch) is NOT re-runnable here -- it needs the
        # generation-time handshake context (resolved polymer_reactants, chip
        # surgery) -- which is exactly why FIX 4 below hard-fails any
        # remaining unclassified live proxy row instead of trusting it.
        # Function-local import: avoids a solver->polymer module cycle.
        # pool_registry (r87 shape B): the reference-state-split
        # isomerization row has no Polymer participant, so the restamp
        # supplies the registered pools collected from the SPECIES lists --
        # real Polymer objects with monomer structures (the r85 dual-axis
        # size gate stays computable), NOT the structureless
        # PolymerPoolConfig sidecars. Collected once per rebuild.
        from rmgpy.polymer import collect_polymer_pool_registry, stamp_gas_association_refusal
        restamp_pool_registry = collect_polymer_pool_registry(
            core_species or [], edge_species or [])
        for rxn in itertools.chain(core_reactions, edge_reactions):
            stamp_gas_association_refusal(rxn, pool_registry=restamp_pool_registry)

        # item 18: capture the Task-3 refusal stamp in the SAME chain(core, edge)
        # order so self.reaction_refused[r_idx] aligns with the residual's r_idx
        # exactly like is_end_group_reaction. A refused FEATURE-radical reaction's
        # whole flux is suppressed in the residual so it fabricates no MW-weighted
        # backbone mass (stamp-but-keep: it stays in the edge/model).
        self.reaction_refused = np.zeros(n_rxn, dtype=np.int8)
        for i, rxn in enumerate(itertools.chain(core_reactions, edge_reactions)):
            if getattr(rxn, "polymer_refused", False):
                self.reaction_refused[i] = 1

        # item 18 (T3): correct-but-loud census of refused FEATURE-radical
        # reactions, tagged by radical class (eliminating vs accumulating) and
        # refuse reason (conduit-deferred vs qssa-invalid /
        # qssa-unassessable). Built once per rebuild over the SAME
        # chain(core, edge) order as the capture above so
        # reaction_refused[i] aligns. Changes NO flux/state -- it only reports
        # (mirrors the thermo reference_state_census posture). Helpers are
        # function-local imports to avoid a solver->polymer module cycle.
        #
        # Round-20 increment 7 (rate-derived QSSA diagnostic, CENSUS-ONLY):
        # for each ACCUMULATING refusal the rebuild computes the lost
        # radical's directional pseudo-first-order outflux k_out [1/s] from
        # the rebuilt model's own kinetics (kf forward / kf/Keq reverse per
        # consuming orientation; co-participant concentrations from the
        # deck's initial state). The census entry carries k_out, the
        # explicit threshold policy and the reactor T. When the radical or
        # its consumers are NOT visible in the rebuilt model the reason is
        # the distinct 'qssa-unassessable' -- missing evidence is never
        # called 'slow'. The computed reason is written back onto
        # rxn.polymer_refused_reason: the artifact emitter
        # (compile_polymer_reaction_entries) reads that same attr, so
        # census and artifact stay single-sourced. reaction_refused is
        # UNTOUCHED -- no flux changes here.
        from rmgpy.polymer import (_warn_once_refused, _reaction_census_label,
                                   assess_refused_qssa_kout)
        # Threshold policy (explicit, logged, deck-tied; never a bare
        # constant): at least 1/terminationTime, raised by the optional
        # deck-tied floor self.qssa_kout_floor_s (callers may set it from
        # the RMG filter threshold). None when the deck provides neither.
        t_term = None
        for _term in (self.termination or []):
            if type(_term).__name__ == "TerminationTime":
                _tv = float(_term.time.value_si)
                t_term = _tv if t_term is None else min(t_term, _tv)
        _kout_floor = float(getattr(self, "qssa_kout_floor_s", 0.0) or 0.0)
        kout_threshold_s = None
        if t_term is not None and t_term > 0.0:
            kout_threshold_s = max(_kout_floor, 1.0 / t_term)
        elif _kout_floor > 0.0:
            kout_threshold_s = _kout_floor
        _T_census = float(self.T.value_si)
        _P_census = float(self.P.value_si)

        def _initial_concentration(spc):
            """Deck-tied initial concentration of a census co-participant
            [mol/m^3]; None when unknown (the direction is then counted
            visible but unquantified). Pool proxies use the pool's initial
            mu1 site density; gas species use the ideal-gas initial state."""
            lbl = getattr(spc, "label", None)
            if lbl and lbl in (self.initial_polymer_moments or {}):
                mom = self.initial_polymer_moments[lbl]
                return max(0.0, float(mom[1])) / self.V_poly
            x = None
            try:
                x = self.initial_mole_fractions.get(spc)
                if x is None and lbl is not None:
                    x = self.initial_mole_fractions.get(lbl)
            except Exception:
                x = None
            if x is None:
                return None
            return float(x) * _P_census / (constants.R * _T_census)

        all_census_rxns = list(itertools.chain(core_reactions, edge_reactions))
        all_census_spcs = list(itertools.chain(core_species or [],
                                               edge_species or []))
        _kout_policy_logged = False
        self.refused_reaction_census = []
        for i, rxn in enumerate(all_census_rxns):
            if not self.reaction_refused[i]:
                continue
            accumulating = bool(getattr(rxn, "polymer_refused_accumulating", False))
            reason = "qssa-invalid" if accumulating else "conduit-deferred"
            # Single-source stamp for EVERY refused row (a stale reason from
            # an earlier rebuild -- e.g. across an r92 flip-restamp that
            # changed the accumulating bit -- must never survive into the
            # artifact emitter). The accumulating branch below may overwrite
            # it with 'qssa-unassessable'.
            rxn.polymer_refused_reason = reason
            entry = {
                "reaction": _reaction_census_label(rxn),
                "radical_class": "accumulating" if accumulating else "eliminating",
                "reason": reason,
            }
            if accumulating:
                if not _kout_policy_logged:
                    _kout_policy_logged = True
                    logging.warning(
                        "QSSA k_out census policy (round-20, census-only): "
                        "threshold = max(deck floor qssa_kout_floor_s=%.6g "
                        "1/s, 1/terminationTime=%s 1/s) = %s 1/s at reactor "
                        "T = %.2f K. k_out below/without evidence is "
                        "reported, never converted into flux.",
                        _kout_floor,
                        ("%.6g" % (1.0 / t_term)) if t_term else "n/a",
                        ("%.6g" % kout_threshold_s)
                        if kout_threshold_s is not None else "n/a",
                        _T_census)
                diag = assess_refused_qssa_kout(
                    rxn, all_census_rxns, all_census_spcs,
                    _T_census, _P_census,
                    concentration_of=_initial_concentration)
                if not diag["visible"]:
                    # the radical / its consumers are not visible in the
                    # rebuilt model: distinct reason, never 'slow'.
                    reason = "qssa-unassessable"
                    entry["reason"] = reason
                k_out = diag["k_out_s"]
                if not diag["visible"] or k_out is None or kout_threshold_s is None:
                    verdict = "unassessable"
                elif k_out >= kout_threshold_s:
                    # a lower bound >= threshold is safely 'fast'
                    verdict = "fast"
                elif diag["k_out_is_lower_bound"]:
                    # under-threshold LOWER BOUND is not evidence of
                    # slowness (unquantified directions remain)
                    verdict = "unassessable"
                else:
                    verdict = "slow"
                entry["k_out_s"] = k_out
                entry["k_out_threshold_s"] = kout_threshold_s
                entry["k_out_verdict"] = verdict
                entry["reactor_T_K"] = _T_census
                entry["radical_label"] = diag["radical_label"]
                # Single-source stamp: the artifact emitter reads this attr.
                rxn.polymer_refused_reason = reason
            self.refused_reaction_census.append(entry)
            _warn_once_refused(entry)

        # Per-reaction pool moment-flux archetype (spec 2026-06-09). Same
        # chain(core, edge) order as is_end_group_reaction so indices match
        # r_idx in the residual.
        self.reaction_flux_archetype = np.zeros(n_rxn, dtype=np.int8)
        self.reaction_src_pool = np.full(n_rxn, -1, dtype=np.int32)
        self.reaction_dst_pool = np.full(n_rxn, -1, dtype=np.int32)
        # Stamped chip repeat-unit counts (spec 2026-06-10 §4.3); same
        # chain(core, edge) order so the index matches r_idx in the residual.
        self.reaction_chip_units = np.zeros(n_rxn, dtype=np.int32)
        # Stamped volatile-ejection unit debit a = Sigma(volatile MW) /
        # source-pool monomer MW (fractional, e.g. ~1.135 for
        # alpha-methylstyrene off a styrene pool). Same chain(core, edge)
        # order so the index matches r_idx in the residual.
        self.reaction_eject_units = np.zeros(n_rxn, dtype=np.float64)
        # M18.4 conduit-admission (arch 7 = FLUX_MOMENT_CREDIT_CONDUIT) row
        # params, flattened from the object-side stamp
        # (rmgpy/polymer.py admit arm: polymer_conduit_dst_pool +
        # polymer_conduit_params). dst is stamped as a pool LABEL string and
        # is the AUTHORITATIVE destination (resolved to an index here, NOT via
        # the product-side pool scan below). u = chain_units is the only field
        # the residual credit arm consumes; conduit_gas_mw + candidate_key ride
        # along for the (later) accepted-step grams booking -- mw+key are
        # sufficient, so NO gas species index is stored (admission stamps the
        # gas token as "label(index)", which a bare-label core map would
        # mis-resolve for a real admitted row). Same chain(core, edge) order so
        # the index matches r_idx in the residual.
        self.conduit_dst_pool_index = np.full(n_rxn, -1, dtype=np.int32)
        self.conduit_chain_units_u = np.zeros(n_rxn, dtype=np.float64)
        self.conduit_gas_units_a = np.zeros(n_rxn, dtype=np.float64)
        self.conduit_gas_mw = np.zeros(n_rxn, dtype=np.float64)
        self.conduit_candidate_key = [None] * n_rxn
        # M18.4 increment 2 (adversarial-review rebuild): per-row TRAPEZOID
        # endpoint memory F_prev_accepted [mol/s] -- the arch-7 net event flux
        # recomputed side-effect-freely from the PREVIOUS accepted state. The
        # accepted-step writer books 0.5*(F_prev + F_curr)*dt*gas_mw, so this
        # holds F at the previous accepted snapshot (0 until first primed). It
        # is NOT written by the residual (the residual leaves state at a
        # Jacobian-FD-perturbed / rejected-trial point, never the accepted
        # solution); _book_conduit_flux_accepted recomputes F from self.y.
        # Re-created (nulled) here on every model rebuild together with the
        # prev-time below, so no trapezoid ever spans a rebuild boundary. Same
        # chain(core, edge) index order as the sibling conduit arrays.
        self.conduit_prev_accepted_flux = np.zeros(n_rxn, dtype=np.float64)
        # Previous ACCEPTED-step time for the conduit grams delta (dt = t_now -
        # t_prev). None => not yet primed (first accepted snapshot books dt=0);
        # re-nulled on every model rebuild so a fresh simulate() never books a
        # stale/negative dt across the boundary.
        self._conduit_prev_accepted_t = None
        _pool_label_to_index = {p.label: j
                                for j, p in enumerate(self.polymer_pools)}
        for i, rxn in enumerate(itertools.chain(core_reactions, edge_reactions)):
            self.reaction_flux_archetype[i] = int(getattr(rxn, "polymer_flux_archetype", 0))
            self.reaction_chip_units[i] = int(getattr(rxn, "polymer_chip_units", 0))
            self.reaction_eject_units[i] = float(getattr(rxn, "polymer_eject_units", 0.0))
            if self.reaction_flux_archetype[i] == FLUX_MOMENT_CREDIT_CONDUIT:
                # Validate the admitted arch-7 stamp BEFORE trusting any field:
                # a missing/degenerate chain_units, absent gas product, or a
                # REVERSIBLE stamp (signed r_mol_s = rf-rr has no reverse
                # moment law -- admitted rows are supposed to be
                # irreversible-export rewritten) all signal mis-stamping.
                # Refuse the whole row LOUDLY (flux-dead), mirroring the
                # dst-less refusal below; never default a bad field to a
                # zero-u / phantom credit. Fail-closed on any parse error (G7).
                params = getattr(rxn, "polymer_conduit_params", None) or {}
                refuse_reason = None
                try:
                    cu = params.get("chain_units", None)
                    gas_products = params.get("gas_products") or []
                    gp0 = gas_products[0] if gas_products else {}
                    gas_mw = gp0.get("mw_g_mol", None)
                    if bool(getattr(rxn, "reversible", False)):
                        refuse_reason = (
                            "it is reversible (admitted conduit rows must be "
                            "irreversible-export rewritten; signed r_mol_s has "
                            "no reverse moment law)")
                    elif (cu is None or not math.isfinite(float(cu))
                            or float(cu) < 1.0):
                        refuse_reason = (
                            "chain_units %r is missing, non-finite, or < 1.0"
                            % (cu,))
                    elif not gas_products or not gp0.get("species"):
                        refuse_reason = (
                            "gas_products is empty or has no species token")
                    elif (gas_mw is None or not math.isfinite(float(gas_mw))
                            or float(gas_mw) <= 0.0):
                        refuse_reason = (
                            "gas product mw_g_mol %r is missing, non-finite, "
                            "or <= 0" % (gas_mw,))
                    elif not params.get("candidate_key"):
                        refuse_reason = "candidate_key is missing"
                except Exception as exc:
                    refuse_reason = ("malformed conduit params (%s: %s)"
                                     % (type(exc).__name__, exc))
                if refuse_reason is not None:
                    logging.warning(
                        "MOMENT-CREDIT CONDUIT: row %d (%s) stamped arch 7 but "
                        "%s; refusing the row (zero flux).",
                        i, str(rxn), refuse_reason)
                    self.reaction_refused[i] = 1
                    self.reaction_flux_archetype[i] = FLUX_NONE
                else:
                    dst_label = getattr(rxn, "polymer_conduit_dst_pool", None)
                    self.conduit_dst_pool_index[i] = _pool_label_to_index.get(
                        dst_label, -1)
                    self.conduit_chain_units_u[i] = float(cu)
                    self.conduit_gas_units_a[i] = float(
                        params.get("gas_units", 0.0))
                    self.conduit_candidate_key[i] = params.get("candidate_key")
                    self.conduit_gas_mw[i] = float(gas_mw)
        # Item 17 (spec 2026-06-12 SS3(e)): generation-time stamps, captured
        # AFTER the stamp-read loop and BEFORE the init demotion pass
        # (:NONE->UNRESOLVED + unresolvable stamped shapes) mutates the
        # array in place. For reactions kinetics-FLIPPED at GENERATION time
        # (restamp_flipped_polymer_archetype, r92: restamp-or-refuse mutates
        # the object itself), "pre-demotion" deliberately means
        # pre-SOLVER-demotion: the captured value is the restamped
        # flipped-direction archetype (or UNRESOLVED for flip-REFUSED rows,
        # which are flux-dead via reaction_refused) and the census reports
        # it as such.
        self.reaction_pre_demotion_archetype = self.reaction_flux_archetype.copy()

        if n_core <= 0:
            raise ValueError(f"Solver received an empty core species list (n_core={n_core}).")

        if len(core_species) != n_core:
            raise ValueError(
                f"Core species length mismatch: len(core_species)={len(core_species)} != num_core_species={n_core}. "
                "Refusing to proceed with inconsistent allocations."
            )

        if self.gas_species_mask is None:
            mask = np.ones(n_core, dtype=bool)
            poly_labels = getattr(self, "polymer_species_labels", set())
            for i in range(n_core):
                if core_species[i].label in poly_labels:
                    mask[i] = False
            self.gas_species_mask = mask
        else:
            if len(self.gas_species_mask) != n_core:
                raise ValueError(
                    f"Provided gas_species_mask length {len(self.gas_species_mask)} != n_core {n_core}. "
                    "Mask must match core size exactly."
                )
            self.gas_species_mask = np.asarray(self.gas_species_mask, dtype=bool)

        # Build mapping: core species index -> polymer pool index (or -1)
        self.species_to_pool_indices = np.full(n_core, -1, dtype=np.int32)
        self.is_pool_proxy = np.zeros(n_core, dtype=np.int8)
        self.pool_mu1_indices.fill(-1)
        self._apply_pool_phase_overrides(self.gas_species_mask, core_species,
                                         record_indices=True)

        # Characteristic-flux include-mask (fix #2a, 2026-06-27). The base
        # enlargement char_rate is an L2 norm over core_species_rates; the
        # moment-coordinate positions (PS_mu0/_mu1/_mu2) carry moment
        # derivatives (the residual writes dn_dt[mu_indices] = dmu*_dt), not
        # molar species fluxes, and contaminate that norm: the lumped k_unzip
        # channel scales char_rate ~linearly in k_unzip and buries family
        # chemistry below tol_move_to_core, leaving the core empty.
        #
        # Authority = the pool mu_indices, NOT is_moment_dummy alone. Species
        # carry an is_moment_dummy flag (set in model.py when the dummies are
        # created), but Species.copy() (species.py:784) does NOT preserve it,
        # so a copied core dummy reads False and an is_moment_dummy-only mask
        # would be silently live-inert. mu_indices are integer core positions
        # the solver binds itself from the pool configs -- copy-proof and
        # exactly the positions the residual treats as moment coordinates. We
        # union with the flag too (cheap; covers any flagged dummy outside an
        # active pool's mu_indices). Everything real -- the proxy, the solvent,
        # the monomer-release routing species -- stays True: we drop bookkeeping
        # coordinates, NOT real chemistry (2a, not 2b). Rebuilt every
        # initialize_model (per RMG iteration) so it tracks the live core.
        include_mask = np.array(
            [not getattr(core_species[i], "is_moment_dummy", False)
             for i in range(n_core)], dtype=bool)
        for pool in self.polymer_pools:
            for idx in pool.mu_indices:
                if 0 <= idx < n_core:
                    include_mask[idx] = False
        self._char_rate_include_mask = include_mask

        # --- prospective_gas_mask (item 17, spec 2026-06-12 SS3(a)) -------
        # A SECOND array over chain(core, edge), never a resize of
        # gas_species_mask (SS3(b): the core size is load-bearing -- six hard
        # raises, ~10 length-coupled consumers and two silent fallbacks key
        # on it). Stage 1 = the blueprint classifier seed (combined-list
        # get_gas_mask, passed by to_solver_object); stage 2 = the SAME
        # per-pool override pass run on the combined list (label-keyed
        # writes scan core+edge; index-keyed writes land identically in the
        # core prefix). NEVER polymer-identity shortcuts: probed INCONSISTENT
        # with the post-promotion mask (spec SS2).
        n_edge_spc = len(edge_species) if edge_species is not None else 0
        combined_species = list(core_species) + (list(edge_species)
                                                 if edge_species else [])
        seed = self._prospective_gas_mask_seed
        # A5-2 (spec 2026-06-12 SS3(a) REDESIGNED): stage-1 precedence, in
        # order. The edge-provenance marker records HOW the edge suffix was
        # produced (1 = stage-1-classified, 0 = default-filled) so R1-EDGE can
        # raise on a production fallback R1's core-prefix check is blind to.
        edge_provenance = np.ones(n_edge_spc, dtype=np.int8)
        if self._prospective_classifier is not None:
            # (1) PRODUCTION PATH: rebuild stage 1 against the LIVE
            # chain(core, edge) via the plumbed classifier handle
            # (polymerPhase.get_gas_mask over the CURRENT lists -- no frozen
            # seed to go stale). The edge suffix is genuinely classified, so
            # provenance stays all-1.
            stage1 = np.asarray(
                self._prospective_classifier(combined_species), dtype=bool)
            if stage1.shape[0] != n_core + n_edge_spc:
                raise ValueError(
                    "PROSPECTIVE-MASK CLASSIFIER: prospective_classifier "
                    "returned length %d for a chain(core, edge) of length %d."
                    % (stage1.shape[0], n_core + n_edge_spc))
            self.prospective_gas_mask = stage1.copy()
        elif seed is not None and len(seed) == n_core + n_edge_spc:
            # (2) DIRECT-TEST SEED PATH: a full doctored seed supplied at
            # construction (e.g. R1's T4). The seed IS a stage-1 product by
            # contract, so provenance is stage-1-classified.
            self.prospective_gas_mask = np.asarray(seed, dtype=bool).copy()
        else:
            # (3) FALLBACK -- genuine last resort, direct-test/runner only.
            # Edge defaults prospectively-GAS; the edge suffix provenance is
            # default-filled (0), so R1-EDGE raises unless the build is flagged
            # allow_default_prospective_edge. The PROSPECTIVE-MASK SEED STALE
            # warning is retained ONLY for back-compat when a stale seed was
            # present (a production build never reaches here -- classifier
            # present => branch 1).
            if seed is not None:
                logging.warning(
                    "PROSPECTIVE-MASK SEED STALE: stage-1 seed length %d != "
                    "n_core + n_edge = %d; rebuilding prospective mask from "
                    "the fallback (core mask + edge defaults GAS).",
                    len(seed), n_core + n_edge_spc)
            # Stated premise, not an implication (amendment A3): edge species
            # default prospectively-GAS here. Fixtures express a
            # prospectively-condensed EDGE species through configured pool
            # labels (stage 2) -- exactly how production daughters become
            # condensed under item 16.
            self.prospective_gas_mask = np.concatenate([
                np.asarray(self.gas_species_mask, dtype=bool).copy(),
                np.ones(n_edge_spc, dtype=bool)])
            edge_provenance = np.zeros(n_edge_spc, dtype=np.int8)
        self._prospective_edge_provenance = edge_provenance
        self._apply_pool_phase_overrides(self.prospective_gas_mask,
                                         combined_species,
                                         record_indices=False)

        # Edge-daughter prospective condensed-mask (spec 2026-06-29). After the
        # stage-2 override, BEFORE R1: flip qualifying EDGE daughters CONDENSED
        # so Gate B (is_poly_event AND not has_condensed_prod) does not zero
        # their real scission flux. The classifier is re-run over the LIVE
        # combined list here (callable, never a frozen set). EDGE slots only --
        # R1 (next) re-verifies the core prefix is untouched. LOUD on failure:
        # a raising classifier must NOT silently default the daughter to GAS
        # (that would re-hide Gate B).
        if self.prospective_condensed_edge_daughter_classifier is not None:
            try:
                condensed_bases = \
                    self.prospective_condensed_edge_daughter_classifier(combined_species)
            except Exception:
                logging.error(
                    "EDGE-DAUGHTER CONDENSED-MASK: classifier raised while "
                    "building prospective_gas_mask (n_core=%d, n_edge=%d); not "
                    "defaulting to GAS.", n_core, n_edge_spc)
                raise
            # LOCKSTEP with the classifier's base-label convention
            # (rmgpy.polymer.strip_rmg_index_suffix, via polymer_input's
            # _base_label): trailing-index strip ONLY. The former first-'('
            # truncation aliased SMILES-labelled species onto shared bases
            # (PP run-5 tripwire crash).
            from rmgpy.polymer import strip_rmg_index_suffix
            for i in range(n_core, n_core + n_edge_spc):
                if strip_rmg_index_suffix(combined_species[i].label) in condensed_bases:
                    self.prospective_gas_mask[i] = False

        # RIDER R1 -- core-prefix parity tripwire (spec SS3(d)). The
        # architecture's central claim ("the prospective mask is the real
        # mask evaluated early") made self-verifying: the gates may use
        # prospective[p] uniformly for all p precisely BECAUSE this prefix
        # is proven equal, every build. Raise, never warn. Re-runs on every
        # rebuild, so spawned daughters and config changes are re-checked.
        if not np.array_equal(
                np.asarray(self.prospective_gas_mask[:n_core], dtype=bool),
                np.asarray(self.gas_species_mask, dtype=bool)):
            diverging = [
                f"index {i} ({core_species[i].label}): core mask says "
                f"{'GAS' if self.gas_species_mask[i] else 'CONDENSED'}, "
                f"prospective says "
                f"{'GAS' if self.prospective_gas_mask[i] else 'CONDENSED'}"
                for i in range(n_core)
                if bool(self.prospective_gas_mask[i])
                != bool(self.gas_species_mask[i])]
            raise ValueError(
                "PROSPECTIVE-MASK TRIPWIRE: prospective_gas_mask core prefix "
                "diverges from gas_species_mask: " + "; ".join(diverging)
                + ". Known suspects: (a) duplicate-label fallback "
                "disablement in the combined get_gas_mask call (a duplicate "
                "label present only in the edge disables label_fallback_safe "
                "for the combined list -- PolymerPhase.get_gas_mask, "
                "polymer_input.py); (b) any base-label producer/consumer "
                "drifting from rmgpy.polymer.strip_rmg_index_suffix (the PP "
                "run-5 crash was first-'(' truncation aliasing SMILES "
                "labels onto shared bases -- fixed 2026-07-05).")

        # RIDER R1-EDGE -- edge-suffix provenance guard (item 17 A5-2, spec
        # SS3(d); DISTINCT from R1 by construction). R1's core-prefix check
        # CANNOT see the §3(a) staging defect: when a build falls through to
        # edge-defaults-GAS, the core PREFIX still matches gas_species_mask
        # verbatim (the fallback copies it), so R1 passes green while the edge
        # SUFFIX is silently default-filled. A PRODUCTION build (not flagged
        # allow_default_prospective_edge) whose edge suffix came from defaults
        # means the live-edge stage-1 rebuild never fired -- the latent no-op
        # A5-2 exists to kill. Raise, never warn. The legitimate last-resort
        # fallback (direct-test/runner, no blueprint phase) is opt-in via the
        # flag. Runs at the same point R1 does, every rebuild.
        if (not self._allow_default_prospective_edge
                and self._prospective_edge_provenance is not None
                and self._prospective_edge_provenance.shape[0] > 0
                and not np.all(self._prospective_edge_provenance)):
            n_default = int(np.count_nonzero(
                self._prospective_edge_provenance == 0))
            raise ValueError(
                "PROSPECTIVE-MASK PROVENANCE: fallback fired during a "
                "production enlargement: %d edge entries, %d default-filled; "
                "a production build must classify the live edge via stage-1 "
                "(prospective_classifier over chain(core, edge)) -- the "
                "rebuild trigger missed an edge-list change, or no classifier "
                "handle was plumbed onto the engine. Direct-test/runner "
                "construction that legitimately needs the fallback must set "
                "allow_default_prospective_edge=True."
                % (self._prospective_edge_provenance.shape[0], n_default))

        # Resolve per-reaction source/target pools from species indices (no
        # label matching). r71 FIX 4 (PP run-5 stall): an unstamped
        # proxy-touching row (archetype NONE after the FIX 2 restamp) that is
        # NOT refused has no classification the solver can conserve mass
        # under -- the former silent NONE -> UNRESOLVED remap (legacy
        # mu1-only flux) is exactly what ran live against run-5's exhausted
        # pool. Such rows now HARD-FAIL unless the build explicitly opts into
        # the legacy fallback (allow_unstamped_proxy_rows=True --
        # direct-test/runner construction only, mirroring
        # allow_default_prospective_edge). Refused rows stay archetype NONE:
        # they are flux-dead via reaction_refused, and routing them onto the
        # legacy path would misreport an adjudicated row as legacy-live.
        ir_arr = self.reactant_indices
        ip_arr = self.product_indices
        n_unstamped = 0
        n_unresolvable = 0
        unstamped_live_rows = []
        combined_rxns = list(itertools.chain(core_reactions, edge_reactions))
        # r91 (PP run-10, spar-adjudicated): spawned-pool demotion refusal
        # bookkeeping. A stamped pool-coupled row whose required pool
        # endpoint resolves to -1 BECAUSE the endpoint is a mid-run-SPAWNED,
        # config-less pool (a Polymer species present in the species lists
        # but absent from self.polymer_pools) is REFUSED conduit-deferred
        # (reaction_refused -> zero whole-row flux, the same mechanics as
        # every other refused row) instead of demoted to the legacy mu1-only
        # transfer: the legacy path applied an unapportioned mu1-only drain
        # (run-10 drove pool polypropylene_mod to mu0=+0.0047 / mu1=-5.6e-5
        # at simulate-leg 21; the r81 negative-beyond-floor tripwire killed
        # the run). The refusal is re-derived EVERY rebuild and deliberately
        # NOT stamped on the reaction object: once item 16 configures the
        # spawned pool, the SAME row resolves and runs fully apportioned.
        # Unresolved endpoints NOT attributable to a spawned config-less
        # Polymer keep the legacy demotion unchanged, censused separately as
        # a loud anomaly (should be impossible for configured pools; the
        # constructible routes are a non-Polymer daughter species, an
        # edge-only proxy, or a duplicate core base label).
        from rmgpy.polymer import Polymer as _Polymer
        from rmgpy.polymer import strip_rmg_index_suffix as _strip_idx
        from rmgpy.polymer import generated_row_superseding_kernels
        _configured_pool_labels = {p.label for p in self.polymer_pools}
        _arch_names = {FLUX_MIGRATION: "MIGRATION",
                       FLUX_SCISSION_FRAGMENT: "SCISSION_FRAGMENT",
                       FLUX_DISCRETE_CHIP: "DISCRETE_CHIP",
                       FLUX_VOLATILE_EJECTION: "VOLATILE_EJECTION",
                       FLUX_MOMENT_CREDIT_CONDUIT: "MOMENT_CREDIT_CONDUIT"}
        spawned_refusals = []      # (row_idx, rxn, archetype name, labels)
        held_superseded = []       # (row_idx, rxn, (pool, channel, kernels))
        demotion_anomalies = []    # (row_idx, rxn, archetype name)
        self.spawned_pool_refusal_census = []
        for i in range(n_rxn):
            src = -1
            for slot in range(3):
                ridx = ir_arr[i, slot]
                if ridx != -1 and ridx < n_core and self.species_to_pool_indices[ridx] != -1:
                    src = self.species_to_pool_indices[ridx]
                    break
            dst = -1
            for slot in range(3):
                pidx = ip_arr[i, slot]
                if pidx != -1 and pidx < n_core:
                    pool_j = self.species_to_pool_indices[pidx]
                    if pool_j != -1:
                        if pool_j != src:
                            dst = pool_j  # prefer the cross-pool product
                            break
                        # pool_j == src here: keep the same-pool fold-back only
                        # as a fallback while no cross-pool product was found.
                        if dst == -1:
                            dst = pool_j
            self.reaction_src_pool[i] = src
            self.reaction_dst_pool[i] = dst
            if self.reaction_flux_archetype[i] == FLUX_MOMENT_CREDIT_CONDUIT:
                # Arch 7 (credit-only conduit): the destination is
                # AUTHORITATIVE from the stamped pool label
                # (conduit_dst_pool_index), NOT the product-side pool scan;
                # src stays NULL (no source-pool debit). Overwrite so
                # reaction_dst_pool -- and hence the SGH destination-pool arm
                # and the residual credit arm -- see the stamped destination.
                dst = int(self.conduit_dst_pool_index[i])
                self.reaction_dst_pool[i] = dst
                if dst == -1:
                    # dst-less conduit credit row: its destination pool label
                    # did not resolve to a configured pool, so it can credit
                    # NOTHING. Refuse it LOUDLY (flux-dead) rather than leaving
                    # a live arch-7 stamp that silently injects zero moment
                    # flux.
                    logging.warning(
                        "MOMENT-CREDIT CONDUIT: row %d (%s) stamped arch 7 but "
                        "its destination pool %r did not resolve to a "
                        "configured pool; refusing the row (zero flux).",
                        i, str(combined_rxns[i]),
                        getattr(combined_rxns[i], "polymer_conduit_dst_pool",
                                None))
                    self.reaction_refused[i] = 1
                    self.reaction_flux_archetype[i] = FLUX_NONE
            if self.reaction_flux_archetype[i] == FLUX_NONE and (src != -1 or dst != -1):
                if self.reaction_refused[i]:
                    # Refused row (adjudicated, flux-dead via
                    # reaction_refused): stays archetype NONE -- never
                    # routed onto the legacy mu1-only path (r71 FIX 4).
                    pass
                elif self._allow_unstamped_proxy_rows:
                    self.reaction_flux_archetype[i] = FLUX_UNRESOLVED
                    n_unstamped += 1
                else:
                    unstamped_live_rows.append(
                        "%s row %d: %s" % (
                            "core" if i < len(core_reactions) else "edge",
                            i, str(combined_rxns[i])))
            arch_i = self.reaction_flux_archetype[i]
            if ((arch_i in (FLUX_MIGRATION, FLUX_SCISSION_FRAGMENT, FLUX_VOLATILE_EJECTION)
                    and (src == -1 or dst == -1))
                    or (arch_i == FLUX_DISCRETE_CHIP
                        and src == -1)):
                # A stamped archetype needs its pool(s) resolved in the
                # solver: MIGRATION/SCISSION_FRAGMENT/VOLATILE_EJECTION need
                # BOTH src and dst
                # (e.g. scission daughters are registered as core Polymer
                # species but have no pool config yet); DISCRETE_CHIP needs
                # only src (no dst: complement folds back to the same pool
                # and the chip is a plain gas species).
                # src == dst is deliberately NOT touched: for that shape
                # (fold-back proxy product + non-pool daughter) the dispatch
                # skip and the legacy transfer produce identical pool flux
                # (reactant -r and fold-back +r cancel on the same mu1).
                #
                # r91: attribute each MISSING required endpoint. If every
                # missing endpoint is explained by a spawned config-less
                # Polymer participant on that side, REFUSE the row
                # (conduit-deferred, zero flux -- archetype cleared to NONE
                # so no dispatch/census path sees a live stamped shape; the
                # generation stamp survives in
                # reaction_pre_demotion_archetype). Otherwise keep the
                # legacy mu1-only demotion (parent drain never silently
                # zeroed) and census the anomaly loudly.
                rxn_i = combined_rxns[i]
                spawned_labels = set()
                attributable = True
                sides = []
                if src == -1:
                    sides.append(rxn_i.reactants or [])
                if dst == -1 and arch_i != FLUX_DISCRETE_CHIP:
                    sides.append(rxn_i.products or [])
                for side in sides:
                    side_labels = [
                        _strip_idx(s.label) for s in side
                        if isinstance(s, _Polymer)
                        and _strip_idx(s.label) not in _configured_pool_labels]
                    if side_labels:
                        spawned_labels.update(side_labels)
                    else:
                        attributable = False
                if attributable and spawned_labels:
                    self.reaction_refused[i] = 1
                    self.reaction_flux_archetype[i] = FLUX_NONE
                    # A daughter the enlarge boundary held back because its
                    # only row is superseded (I-253): the row's name is the
                    # supersession, read by the same species-level rule.
                    superseded_i = generated_row_superseding_kernels(rxn_i)
                    if superseded_i:
                        held_superseded.append((i, rxn_i, superseded_i))
                    else:
                        spawned_refusals.append(
                            (i, rxn_i, _arch_names[arch_i],
                             sorted(spawned_labels)))
                else:
                    self.reaction_flux_archetype[i] = FLUX_UNRESOLVED
                    n_unresolvable += 1
                    demotion_anomalies.append((i, rxn_i, _arch_names[arch_i]))
        if spawned_refusals:
            self.spawned_pool_refusal_census = [
                {"reaction": _reaction_census_label(r_rxn),
                 "archetype": r_arch,
                 "pools": r_pools,
                 "reason": "conduit-deferred"}
                for (_, r_rxn, r_arch, r_pools) in spawned_refusals]
            logging.warning(
                "SPAWNED-POOL DEMOTION REFUSAL: %d stamped pool-coupled "
                "row(s) target spawned polymer pool(s) with no solver pool "
                "config (the item-16 enlarge-boundary promotion did not "
                "reach them: a daughter without spawn markers, or one still "
                "edge-resident); refused conduit-deferred instead of demoted "
                "to legacy mu1-only; archetypes=%s, pools=%s, first_rows=%s",
                len(spawned_refusals),
                ",".join(sorted({a for (_, _, a, _) in spawned_refusals})),
                ",".join(sorted({lbl for (_, _, _, ls) in spawned_refusals
                                 for lbl in ls})),
                "; ".join("%s row %d: %s" % (
                    "core" if r_i < len(core_reactions) else "edge",
                    r_i, str(r_rxn))
                    for (r_i, r_rxn, _, _) in spawned_refusals[:3]))
        if n_unresolvable:
            logging.warning(
                "%d reactions stamped MIGRATION/SCISSION_FRAGMENT/DISCRETE_CHIP "
                "could not resolve their solver pool(s); demoted to legacy "
                "mu1-only moment flux (UNRESOLVED).", n_unresolvable)
            logging.warning(
                "CONFIGURED-POOL UNRESOLVED DEMOTION ANOMALY: %d stamped "
                "pool-coupled row(s) had an unresolved required pool "
                "endpoint NOT attributable to an unconfigured spawned pool "
                "(should be impossible for configured pools; constructible "
                "routes: non-Polymer daughter species, edge-only proxy, "
                "duplicate core base label); legacy mu1-only demotion kept "
                "unchanged. rows=%s",
                len(demotion_anomalies),
                "; ".join("%s row %d [%s]: %s" % (
                    "core" if a_i < len(core_reactions) else "edge",
                    a_i, a_arch, str(a_rxn))
                    for (a_i, a_rxn, a_arch) in demotion_anomalies[:5]))
        if unstamped_live_rows:
            raise ValueError(
                "UNSTAMPED PROXY ROW(S): %d proxy-touching reaction(s) "
                "arrived at the solver rebuild with no polymer_flux_archetype "
                "stamp and no polymer_refused adjudication -- the legacy "
                "mu1-only fallback would apply silent unclassified flux to "
                "the pool (the PP run-5 stall channel, adjudicated r71). "
                "Offending row(s): %s. A production build must classify "
                "every live proxy row at generation "
                "(stamp_polymer_flux_archetype) or refuse it "
                "(stamp_gas_association_refusal). Direct-test/runner "
                "construction that deliberately exercises the legacy "
                "mu1-only fallback must set allow_unstamped_proxy_rows=True."
                % (len(unstamped_live_rows), "; ".join(unstamped_live_rows)))
        if n_unstamped:
            logging.warning(
                "%d proxy-touching reactions arrived without a polymer_flux_archetype "
                "stamp; applying legacy mu1-only pool moment flux for them "
                "(allow_unstamped_proxy_rows escape hatch -- direct-test/"
                "runner construction only).",
                n_unstamped)

        # LIVE-UNRESOLVED census (adjudicated): enumerate every non-refused
        # FLUX_UNRESOLVED proxy-touching row that will apply the legacy
        # mu1-only fallback at the RHS. Sources: generation-time UNRESOLVED
        # classification (stamp_polymer_flux_archetype) OR the solver-side
        # non-attributable pool-demotion above. Pure additive diagnostic --
        # ONE prominent warning at initialize_model, no raise, no flux
        # change: mu0's sign is unknowable for these shapes
        # (survivor-vs-consumed lineage is unstamped), so the fallback stays
        # explicit-and-loud rather than silently "corrected". Proxy-touch
        # membership deliberately reuses the SAME check as the
        # unstamped_live_rows gate / SGH mover scan (is_pool_proxy over
        # reactant_indices/product_indices core slots).
        live_unresolved_rows = []
        for i in range(n_rxn):
            if self.reaction_flux_archetype[i] != FLUX_UNRESOLVED:
                continue
            if self.reaction_refused[i]:
                continue
            touches_proxy = False
            for slot in range(3):
                r_sp = ir_arr[i, slot]
                if r_sp != -1 and r_sp < n_core and self.is_pool_proxy[r_sp]:
                    touches_proxy = True
                    break
                p_sp = ip_arr[i, slot]
                if p_sp != -1 and p_sp < n_core and self.is_pool_proxy[p_sp]:
                    touches_proxy = True
                    break
            if touches_proxy:
                live_unresolved_rows.append(i)
        if live_unresolved_rows:
            logging.warning(
                "LIVE UNRESOLVED MOMENT-FLUX CENSUS: %d non-refused "
                "FLUX_UNRESOLVED proxy-touching row(s) will apply legacy "
                "mu1-only flux at the RHS (updates mu1 only; CANNOT conserve "
                "mu0/mu2 for these shapes). Sources: generation-time "
                "UNRESOLVED classification OR solver-side non-attributable "
                "pool-demotion. first rows: %s",
                len(live_unresolved_rows),
                "; ".join("%s row %d: %s" % (
                    "core" if r_i < len(core_reactions) else "edge",
                    r_i, str(combined_rxns[r_i]))
                    for r_i in live_unresolved_rows[:5]))

        # P1-2 (atom-transfer mass defect): per-row MOMENT-space eject units.
        # reaction_eject_units stores the MASS defect a_mass = (net
        # non-polymer mass)/monomer_mw. For an atom-transfer conduit
        # (H + parent <=> H2 + daughter) that mass is carried by the dst
        # pool's per-chain chain_mass_defect_g_mol, NOT by a repeat-unit
        # decrement: booking the moment shift with a_mass lands monomeric
        # chains at mu1/mu0 = 1 - a, below the realizability cone mu1 >= mu0
        # BY CONSTRUCTION (regen-#2 forensics: parent ground to DP~1 by
        # scission, every forward event seeds the daughter unrealizable, the
        # mu3 closure freezes, the reverse mu1 drain loses its self-limit).
        # The moment shift must instead use
        #     a_moment = a_mass - (defect_dst - defect_src)/monomer_mw,
        # which is exactly 0 for H-loss conduits (chains transfer with
        # UNCHANGED mu1; the gas still gains the H mass through the ordinary
        # species stoichiometry, and the condensed-mass closure
        # mu1*MW - mu0*defect books the loss) and exactly a_mass for true
        # monomer/chip ejection rows (defect_delta = 0 -- byte-identical).
        # Placement is load-bearing: AFTER the demotion loop, so src/dst are
        # FINAL. Exact-cancellation roundoff (structural defect vs
        # reaction-mass a, both g/mol sums differing by summation order) is
        # snapped to a true zero shift.
        from rmgpy.polymer import volatile_ejection_moment_units
        self.reaction_eject_units_moment = self.reaction_eject_units.copy()
        for r_ve in range(n_rxn):
            if self.reaction_flux_archetype[r_ve] != FLUX_VOLATILE_EJECTION:
                continue
            src_ve = self.reaction_src_pool[r_ve]
            dst_ve = self.reaction_dst_pool[r_ve]
            if src_ve == -1 or dst_ve == -1 or src_ve == dst_ve:
                continue
            defect_src_ve = float(getattr(
                self.polymer_pools[src_ve], "chain_mass_defect_g_mol", 0.0) or 0.0)
            defect_dst_ve = float(getattr(
                self.polymer_pools[dst_ve], "chain_mass_defect_g_mol", 0.0) or 0.0)
            if defect_dst_ve == defect_src_ve:
                continue
            mw_src_ve = float(getattr(
                self.polymer_pools[src_ve], "monomer_mw_g_mol", 0.0) or 0.0)
            if mw_src_ve <= 0.0:
                raise ValueError(
                    f"VOLATILE_EJECTION row {r_ve} connects pools with "
                    f"differing chain_mass_defect_g_mol "
                    f"({defect_src_ve!r} -> {defect_dst_ve!r}) but source "
                    f"pool '{self.polymer_pools[src_ve].label}' has no "
                    f"positive monomer_mw_g_mol; cannot convert the mass "
                    f"defect to a moment shift without fabricating mass.")
            self.reaction_eject_units_moment[r_ve] = volatile_ejection_moment_units(
                float(self.reaction_eject_units[r_ve]), defect_src_ve,
                defect_dst_ve, mw_src_ve)

        # Kernel supersession by channel identity (I-253 ruling; it replaces
        # the I-250 row-class rule and the item-18 tripwire): a live generated
        # row whose reaction channel a kernel on the pool it decomposes
        # carries (rmgpy.polymer.kernels_carrying_channel) counts that
        # kernel's events again. It is refused by name (archetype NONE, zero
        # flux; non-sticky: re-derived per rebuild, the reaction object is
        # not stamped) and censused. Every other row stays live and is never
        # counted as a double count, whatever its archetype. Placement is
        # load-bearing: AFTER the demotion loop and the moment-shift loop
        # above (FINAL archetypes, pools and moment shifts: they decide
        # whether a row moves repeat units) and BEFORE the double-count
        # census, the catalog count-once guard and the stationary-carrier
        # guard (validate_configuration), so none of them sees a superseded
        # row.
        self._channel_species = combined_species
        self._held_superseded = held_superseded
        # Proxy balance (I-255): a pool-coupled row whose proxies change the
        # net charge, or no longer read the channel stamped from its raw
        # molecules, is refused by name first, so a defective proxy is named
        # for its defect rather than superseded or integrated.
        self._refuse_unbalanced_proxy_rows(combined_rxns)
        self._refuse_kernel_superseded_rows(combined_rxns)
        self._census_kernel_channel_double_counts(combined_rxns)

        # Enforce the moment-isolation invariant and pool/mass-transfer index
        # sanity now that gas_species_mask and the reaction index tables are
        # populated. Moments must evolve only via the tail block, never through
        # generic reaction stoichiometry.
        # RIDER R2 static half (spec SS3(e)): once per rebuild by
        # construction -- the enumeration runs exactly once per
        # initialize_model; no keying set needed.
        self._static_phase_gate_census(core_species, core_reactions)

        # read by the catalog-kernel overlap check (a reversible row moves
        # chains both ways)
        self._core_row_reversible = [bool(getattr(rxn, "reversible", False))
                                     for rxn in core_reactions]
        self.validate_configuration()

        # radical_qssa_unzip double-count census (M2), kernel against kernel:
        # QSSA initiation is random backbone homolysis, as is the pool's own
        # k_scission. A generated row beside the QSSA kernel is the channel
        # rule's (superseded above when the kernel carries its channel, live
        # and uncounted otherwise). WARN-ONCE, NEVER refuse; the k_unzip
        # co-presence stays a hard error (M1 mutual exclusion). Placement:
        # AFTER validate_configuration (qssa_enabled -- the flattened
        # solver-owned gate, the ONLY channel signal the solver trusts --
        # exists here).
        from rmgpy.polymer import warn_once_qssa_double_count
        self.qssa_double_count_census = []
        for p_idx, pool in enumerate(self.polymer_pools):
            if not self.qssa_enabled[p_idx]:
                continue
            # QSSA + the pool's own k_scission: BOTH are random backbone
            # homolysis -- the most direct initiation double-count of all
            # (unlike generated scission, k_scission cannot claim to cover
            # different bonds by construction; only the user knows whether
            # the two were parameterized for disjoint physics). Same
            # warn-once helper, distinct census key.
            if pool.k_scission > 0.0:
                entry = {"pool": pool.label, "overlap": "k_scission",
                         "k_scission": pool.k_scission}
                self.qssa_double_count_census.append(entry)
                warn_once_qssa_double_count(entry)

        # k_homolysis supersession census (Stage 1, adjudicated round 66):
        # refused explicit homolysis/association rows co-existing with the
        # kernel are FINE -- they carry zero flux (reaction_refused) -- but
        # the census must state the pairing EXPLICITLY: a warn-level line
        # names the refused conduit-deferred rows as superseded by the
        # kernel. Placement mirrors the QSSA census above (after
        # validate_configuration: khom_enabled -- the flattened solver-owned
        # gate -- exists here; after the refused census: superseded rows are
        # final). Diagnostic-only: NO flux/state change, warn NEVER refuse.
        from rmgpy.polymer import warn_once_homolysis_supersession
        self.homolysis_supersession_census = []
        deferred_rows = [c["reaction"] for c in self.refused_reaction_census
                         if c.get("reason") == "conduit-deferred"]
        if deferred_rows and np.any(self.khom_enabled):
            for p_idx, pool in enumerate(self.polymer_pools):
                if not self.khom_enabled[p_idx]:
                    continue
                entry = {"pool": pool.label,
                         "superseded_rows": list(deferred_rows)}
                self.homolysis_supersession_census.append(entry)
                warn_once_homolysis_supersession(entry)

        # side_group_homolysis supersession census (FR1-K1): refused
        # explicit gas-radical<->condensed rows REMAIN refused next to the
        # live kernel (they carry zero flux -- reaction_refused), and the
        # census pairs them explicitly, same contract as the k_homolysis
        # census above. Diagnostic-only: NO flux/state change, warn NEVER
        # refuse. Same placement rationale (after validate_configuration:
        # sgh_enabled exists; after the refused census: rows are final).
        from rmgpy.polymer import warn_once_side_group_supersession
        self.side_group_supersession_census = []
        if deferred_rows and np.any(self.sgh_enabled):
            for p_idx, pool in enumerate(self.polymer_pools):
                if not self.sgh_enabled[p_idx]:
                    continue
                entry = {"pool": pool.label,
                         "superseded_rows": list(deferred_rows)}
                self.side_group_supersession_census.append(entry)
                warn_once_side_group_supersession(entry)

        # Thermo reference-state tripwire (spec 2026-06-11 §7): runs AFTER
        # the archetype demotion pass and validate_configuration (masks,
        # pool membership and archetypes are final here) and BEFORE
        # generate_rate_coefficients computes any kb from Keq -- refusal
        # gates SOLVING. Re-runs on every rebuild, so spawned daughters are
        # checked automatically.
        self._reference_state_tripwire(core_species, core_reactions)

        self._scratch_C_gas = np.zeros(n_core, float)
        self._scratch_C_poly = np.zeros(n_core, float)
        self._scratch_dn_dt = np.zeros(n_core, float)
        self._scratch_du_dt = np.zeros(self.num_qssa_u, float)
        self._scratch_dz_dt = np.zeros(self.num_sgh_z, float)
        self._scratch_dl_dt = np.zeros(self.num_end_ledger, float)
        self._core_species_atoms = self._collect_species_atoms(core_species)
        self._core_labels = [getattr(sp, 'label', '') for sp in core_species]
        self._scratch_proxy_activity = np.zeros(n_core, float)

        # RIDER R2 dynamic half (item 17, spec 2026-06-12 SS3(e)): the edge
        # counterfactual -- what enlargement WOULD have seen absent the
        # gates. Ungated rows mirror their real writes; gate-zeroed edge
        # rows carry their counterfactual here and ONLY here. Warn-once
        # keying (gate_code, edge reaction index) clears per engine rebuild
        # = per RMG iteration: a persistent gated channel re-announces once
        # per iteration, deliberately (correct-but-loud; measured ~14
        # lines/iteration on the reference EPDM deck).
        self.edge_reaction_gate_code = np.zeros(self.num_edge_reactions,
                                                dtype=np.int8)
        self.edge_reaction_rates_ungated = np.zeros(self.num_edge_reactions,
                                                    float)
        self.edge_species_rates_ungated = np.zeros(self.num_edge_species,
                                                   float)
        self._phase_gate_census_emitted = set()

        # Exhaustion-tail conditioning (r81 B): per-pool per-moment floors
        # derived from the solver ABSOLUTE tolerances -- floor_k =
        # max(SMALL_EPS, EXHAUSTION_FLOOR_K * atol[state]) on the pool's own
        # mu0/mu1/mu2 state slots (mole basis, same basis atol governs).
        # Census keying is per-rebuild, mirroring the sibling warn-once sets.
        n_pools_exh = len(self.polymer_pools)
        self._pool_exhausted = np.zeros(n_pools_exh, dtype=np.uint8)
        # P1-3: most-negative moment value seen per pool on TRIAL states
        # since the last rebuild (diagnostic census; residual never raises).
        self._pool_worst_trial_excursion = np.zeros(n_pools_exh, dtype=float)
        self._exhaustion_census_emitted = set()
        self._pool_mu_floors = np.full((n_pools_exh, 3), SMALL_EPS,
                                       dtype=float)
        atol_arr_exh = getattr(self, "atol_array", None)
        for p in range(n_pools_exh):
            mu_idx_exh = self.polymer_pools[p].mu_indices
            for k in range(3):
                s_idx = int(mu_idx_exh[k])
                if atol_arr_exh is None:
                    a_exh = 0.0
                elif np.ndim(atol_arr_exh) == 0:
                    a_exh = float(atol_arr_exh)
                elif s_idx < len(atol_arr_exh):
                    a_exh = float(atol_arr_exh[s_idx])
                else:
                    a_exh = float(np.max(atol_arr_exh))
                self._pool_mu_floors[p, k] = max(
                    SMALL_EPS, EXHAUSTION_FLOOR_K * a_exh)
        # Realizability census (record_realizability_census), per rebuild.
        self.realizability_census = {
            pool.label: {"accepted_steps": 0, "steps_outside": 0,
                         "first_outside_t": None, "last_outside_t": None,
                         "worst_variance_deficit": 0.0,
                         "worst_cone_gap": 0.0, "excursions": []}
            for pool in self.polymer_pools}
        self._realizability_census_t = None

        # M18.2d: softclamp decay length = the r81 floor scale (scalar-atol
        # decks; uniform across slots -- per-state generalization deferred,
        # see the m18d design doc).
        self._softclamp_lam = max(SMALL_EPS,
                                  EXHAUSTION_FLOOR_K * float(atol))

        self.get_const_spc_indices(core_species)
        self.set_initial_conditions()

        # I-244: the relative accepted floor's scale (the peak |mu_k| over
        # accepted states, from the initial state on) and rtol per moment
        # slot; the exhaustion events and residue ledger, per rebuild.
        rtol_arr_exh = getattr(self, "rtol_array", None)
        self._pool_mu_peak = np.zeros((n_pools_exh, 3), dtype=float)
        self._pool_mu_rtol = np.zeros((n_pools_exh, 3), dtype=float)
        for p in range(n_pools_exh):
            mu_idx_exh = self.polymer_pools[p].mu_indices
            for k in range(3):
                s_idx = int(mu_idx_exh[k])
                self._pool_mu_peak[p, k] = abs(float(self.y0[s_idx]))
                if rtol_arr_exh is None:
                    r_exh = 0.0
                elif np.ndim(rtol_arr_exh) == 0:
                    r_exh = float(rtol_arr_exh)
                elif s_idx < len(rtol_arr_exh):
                    r_exh = float(rtol_arr_exh[s_idx])
                else:
                    r_exh = float(np.max(rtol_arr_exh))
                self._pool_mu_rtol[p, k] = r_exh
        self._pool_ended = np.zeros(n_pools_exh, dtype=np.uint8)
        self._pool_exhaustion_residue = np.zeros((n_pools_exh, 3), dtype=float)
        self._pool_exhaustion_z = {}
        self.pool_exhaustion_events = []

        if filter_reactions:
            ReactionSystem.set_initial_reaction_thresholds(self)
        self.generate_rate_coefficients(core_reactions, edge_reactions)
        ReactionSystem.compute_network_variables(self, pdep_networks)

        ReactionSystem.set_initial_derivative(self)
        # 18.2e P5: (re)arm or clear the pydas `jacobian` capability
        # attribute BEFORE initialize_solver -- pydas detects the user
        # Jacobian by attribute presence inside DASPK.initialize (probe-
        # proven; detection does NOT re-run per step). Runs on every
        # rebuild so the layout validation tracks the live model.
        self._setup_scoped_jacobian()

        # i061 moment error-weight floor (see MOMENT_EWT_FLOOR_K). Applied HERE,
        # AFTER _pool_mu_floors (4060), _softclamp_lam (4082) and the scoped
        # Jacobian's _jac_wt_atol have all been computed from the ORIGINAL
        # atol_array, and immediately BEFORE initialize_solver hands the array to
        # DASPK -- so the only consumer that sees the floored value is DASPK's
        # own error weight. Chemistry noise-floor consumers are left bitwise
        # unchanged: _chem_atol_array preserves the pre-floor array for the one
        # LIVE reader (attribution trust, get_polymer_pool_stats). Scope is the
        # complement of _char_rate_include_mask -- exactly the pool mu_indices
        # and is_moment_dummy core positions (bookkeeping coordinates), NEVER a
        # physical species -- a claim that i075 turned into a GUARD rather than
        # a comment: _assert_moment_slots_carry_no_molar_amount refuses the
        # floor outright if any slot in the mask complement looks like real
        # chemistry. Inert when there are no moment slots.
        self._chem_atol_array = np.array(self.atol_array, dtype=float)
        if self._char_rate_include_mask is not None:
            _floor = MOMENT_EWT_FLOOR_K * atol
            _moment_slots = [i for i in range(self.num_core_species)
                             if not self._char_rate_include_mask[i]]
            self._assert_moment_slots_carry_no_molar_amount(
                _moment_slots, core_species, _floor)
            _floored = []
            for i in _moment_slots:
                if self.atol_array[i] < _floor:
                    self.atol_array[i] = _floor
                    _floored.append(i)
            if _floored:
                _labels = ", ".join(
                    "%d:%s" % (i, getattr(core_species[i], "label", "?"))
                    for i in _floored)
                logging.info(
                    "i061 moment error-weight floor: raised DASPK atol to "
                    "%.3e on %d moment-coordinate slot(s) [%s]; physical "
                    "species untouched.", _floor, len(_floored), _labels)

        ReactionSystem.initialize_solver(self)

        self.diagnose_polymer_mapping(core_species)

    def _assert_moment_slots_carry_no_molar_amount(self, moment_slots,
                                                   core_species, floor):
        """i075 d3: refuse the moment error-weight floor on anything that could
        be a physical species.

        The floor's scope is the complement of _char_rate_include_mask. That
        complement is NOT provably free of physical species from the mask's own
        construction:

          * the flag arm (`is_moment_dummy`) is a plain mutable Species
            attribute, set at exactly one site (rmgpy/rmg/model.py) but never
            re-validated afterwards; and
          * the authoritative arm, a pool's `mu_indices`, is resolved by LABEL
            -- polymer_input.to_config maps `mu_species` through the core
            spc_map, and those mu_species come from
            `species_dict["{proxy}_mu{k}"]`. `_register_polymer` deliberately
            SKIPS creating the dummy when a species of that label "already
            exists (e.g., from the input file)", and no code path checks
            `is_moment_dummy` on the species it reuses. A deck that names a
            real species `PS_mu0` therefore binds it as a moment coordinate.

        So the claim is enforced here instead of assumed. Three properties of a
        bookkeeping coordinate, each independent of the mask that selected the
        slot and each a direct statement that the slot carries a REAL MOLAR
        AMOUNT when it fails:

          1. it takes part in no reaction anywhere in the network -- a moment
             coordinate carries dmu/dt, never a chemical flux;
          2. the deck declares no mole fraction for it -- a moment coordinate's
             t=0 value comes from initial_polymer_moments, never from
             initial_mole_fractions;
          3. it is not classified GAS -- moment coordinates are condensed-phase
             bookkeeping, and a gas slot is by definition a molar amount in the
             gas volume.

        Any of the three failing means relaxing DASPK's error control -- two
        decades, silently -- on something that is not a bookkeeping coordinate,
        which is exactly the failure this guard exists to make impossible. Hard
        error, never a warning: the alternative is a quietly under-resolved
        species.

        NOT checked, deliberately: `reactive`. Moment dummies are created
        reactive=False by rmgpy/rmg/model.py, but nothing in the solver contract
        requires it, and every synthetic pool fixture in the solver test suite
        builds its mu-species with the Species default (reactive=True). Keying
        on it would refuse legal models. LIMIT of what is checked: a species
        that is genuinely inert, carries no deck loading and is condensed --
        an unreacting condensed diluent -- would still pass all three arms if it
        were mis-bound into a mu slot. That case is unreachable through
        polymer_input (a mu slot is bound only from `{proxy}_mu{k}`), but it is
        not excluded by this guard.
        """
        cdef int i
        if not moment_slots:
            return
        _network_species = None
        _mf = getattr(self, "initial_mole_fractions", None) or {}
        _mf_labels = set()
        for _k in _mf:
            _lbl = getattr(_k, "label", None)
            _mf_labels.add(_lbl if _lbl is not None else _k)
        gas_mask = self.gas_species_mask
        for i in moment_slots:
            spc = core_species[i] if i < len(core_species) else None
            label = getattr(spc, "label", "?")
            reason = None
            if _network_species is None:
                _network_species = set()
                for arr in (self.reactant_indices, self.product_indices):
                    if arr is None:
                        continue
                    flat = np.asarray(arr).ravel()
                    _network_species.update(int(v) for v in flat[flat >= 0])
            if i in _network_species:
                reason = ("it appears in the reaction network "
                          "(reactant_indices/product_indices), so it carries "
                          "chemical flux")
            elif spc is not None and (spc in _mf or label in _mf_labels):
                reason = ("the deck declares an initial MOLE FRACTION for it; "
                          "a moment coordinate is loaded through "
                          "initial_polymer_moments")
            elif (gas_mask is not None and i < len(gas_mask)
                    and bool(gas_mask[i])):
                reason = ("it is classified GAS, so its state slot is a molar "
                          "amount in the gas volume, not condensed-phase "
                          "bookkeeping")
            if reason is not None:
                raise ValueError(
                    "i061 moment error-weight floor: core slot %d (%r) is in "
                    "the complement of _char_rate_include_mask -- the solver "
                    "is treating it as a moment COORDINATE and would relax "
                    "its DASPK atol to %.3e -- but %s. A physical species must "
                    "never receive the moment floor. Fix the pool's "
                    "mu_indices / is_moment_dummy binding rather than the "
                    "floor: a species labelled like a moment dummy but "
                    "carrying real molar amounts is already being integrated "
                    "as dmu/dt."
                    % (i, label, float(floor), reason))

    def diagnose_polymer_mapping(self, core_species):
        w = 90
        print(f"\n{'=' * w}")
        print(f"{'POLYMER SOLVER DIAGNOSTIC':^{w}}")
        print(f"{'=' * w}")

        # --- Per-pool sections ---
        for pool_i, pool in enumerate(self.polymer_pools):
            # k_homolysis kernel + end-radical daughter pools must be
            # traceable in run logs (round 66: radical inventories); the
            # daughters themselves enumerate here as ordinary pools, with
            # their role tagged so a reader can find their mu0/mu1 growth.
            khom = pool.k_homolysis
            if khom is not None:
                khom_str = (f"A={khom['A']:.2e},n={khom['n']:g},"
                            f"Ea={khom['Ea']:.4g}J/mol")
            else:
                khom_str = "off"
            # side_group_homolysis channels (FR1-K1) must be traceable in
            # run logs too: the parent names every channel; feature pools
            # are tagged; the v1 saturation limitation is disclosed.
            sgh_channels = pool.side_group_homolysis
            if sgh_channels:
                sgh_str = ("[" + ",".join(ch["label"] for ch in sgh_channels)
                           + "]")
            else:
                sgh_str = "off"
            # k_depropagation kernel (r74 SS2): the radical-end consumption
            # channel must be traceable in run logs too.
            kdep = pool.k_depropagation
            if kdep is not None:
                kdep_str = (f"A={kdep['A']:.2e},n={kdep['n']:g},"
                            f"Ea={kdep['Ea']:.4g}J/mol")
            else:
                kdep_str = "off"
            role_tag = ""
            if pool.label.endswith(K_HOMOLYSIS_DAUGHTER_SUFFIXES):
                role_tag = "  [end-radical daughter]"
            elif SIDE_GROUP_DAUGHTER_INFIX in pool.label:
                role_tag = "  [side-group X-loss daughter]"
            print(f"\n--- Pool {pool_i}: '{pool.label}'  (xs={pool.xs}, "
                  f"k_scission={pool.k_scission:.2e}, "
                  f"k_unzip={pool.k_unzip:.2e}, "
                  f"k_homolysis={khom_str}, "
                  f"k_depropagation={kdep_str}, "
                  f"side_group_homolysis={sgh_str}){role_tag} ---")
            if sgh_channels:
                # SGH kernel-v2: Br leaves the melt through an explicit
                # per-channel Z inventory (Z = moles of removable X sites,
                # Z(0) = sites*mu1(0), dZ/dt = -R). Report each row's seeded
                # inventory and Z slot so a reader can trace the depletion.
                rs = int(self.sgh_row_start[pool_i])
                re_ = int(self.sgh_row_end[pool_i])
                for sgh_row in range(rs, re_):
                    z_slot = int(self.sgh_z_slot[sgh_row])
                    if z_slot < 0:
                        print(f"  side_group_homolysis channel "
                              f"'{self.sgh_channel_labels[sgh_row]}': X-loss "
                              f"feature daughter (mu0 slot "
                              f"{int(self.sgh_dst_mu0[sgh_row])}), M_X="
                              f"{self.sgh_mw_X[sgh_row]:g} g/mol -> gas index "
                              f"{int(self.sgh_gas[sgh_row])}.")
                        continue
                    z0 = self.y0[z_slot] if 0 <= z_slot < len(self.y0) \
                        else float('nan')
                    print(f"  side_group_homolysis[v2] channel "
                          f"'{self.sgh_channel_labels[sgh_row]}': "
                          f"Z inventory slot {z_slot}, Z(0)={z0:.4e} mol "
                          f"(= sites*mu1(0)), M_X={self.sgh_mw_X[sgh_row]:g} "
                          f"g/mol -> gas index {int(self.sgh_gas[sgh_row])}.")
            print(f"  {'Role':<16} {'Species Label':<30} {'Index':<7} {'y0':>13}")
            print(f"  {'-'*16} {'-'*30} {'-'*7} {'-'*13}")

            # Moments
            for k, mu_label in enumerate(["mu0", "mu1", "mu2"]):
                idx = pool.mu_indices[k]
                lbl = core_species[idx].label if 0 <= idx < len(core_species) else "???"
                val = self.y0[idx] if 0 <= idx < len(self.y0) else float('nan')
                print(f"  {'Moment ('+mu_label+')':<16} {lbl:<30} {idx:<7} {val:>13.4e}")

            # Proxy species
            for i in range(len(core_species)):
                if self.species_to_pool_indices[i] == pool_i and self.is_pool_proxy[i] and i not in pool.mu_indices:
                    print(f"  {'Proxy':<16} {core_species[i].label:<30} {i:<7} {self.y0[i]:>13.4e}")

            # Explicit oligomers
            if pool.explicit_dp_to_species_index:
                for dp in sorted(pool.explicit_dp_to_species_index):
                    idx = pool.explicit_dp_to_species_index[dp]
                    lbl = core_species[idx].label if 0 <= idx < len(core_species) else "???"
                    val = self.y0[idx] if 0 <= idx < len(self.y0) else float('nan')
                    print(f"  {'Explicit DP='+str(dp):<16} {lbl:<30} {idx:<7} {val:>13.4e}")

            # Monomer
            if pool.monomer_poly_index is not None:
                idx = pool.monomer_poly_index
                lbl = core_species[idx].label if 0 <= idx < len(core_species) else "???"
                val = self.y0[idx] if 0 <= idx < len(self.y0) else float('nan')
                print(f"  {'Monomer(poly)':<16} {lbl:<30} {idx:<7} {val:>13.4e}")

            # Custom kinetics
            if pool.tail_kinetics is not None:
                print(f"  Custom tail_kinetics: YES")

        # --- Unassigned non-gas species (catch anything missed) ---
        unassigned = []
        for i in range(len(core_species)):
            if self.species_to_pool_indices[i] == -1 and not self.gas_species_mask[i]:
                unassigned.append(i)
        if unassigned:
            print(f"\n--- Unassigned non-gas species (not in any pool) ---")
            print(f"  {'Species Label':<30} {'Index':<7} {'y0':>13}")
            print(f"  {'-'*30} {'-'*7} {'-'*13}")
            for i in unassigned:
                print(f"  {core_species[i].label:<30} {i:<7} {self.y0[i]:>13.4e}")

        # --- Phase / Volume summary (DEBUG) ---
        n_gas_species = int(np.sum(self.gas_species_mask))
        n_poly_species = len(core_species) - n_gas_species
        n_gas_moles = float(np.sum(self.y0[:len(core_species)][self.gas_species_mask]))
        n_poly_moles = float(np.sum(self.y0[:len(core_species)][~self.gas_species_mask]))
        n_gas_via_residual = float(np.sum(self.y0[:len(core_species)] * self.gas_species_mask))
        print(f"\n--- Phase / Volume summary (DEBUG) ---")
        print(f"  T.value_si = {self.T.value_si}  P.value_si = {self.P.value_si}  "
              f"R = {constants.R}  constant_gas_volume = {self.constant_gas_volume}  V_gas0 = {self.V_gas0}")
        print(f"  n_gas (mask*y0 sum) = {n_gas_via_residual:.6e} mol   "
              f"ideal V_gas = {constants.R*self.T.value_si*n_gas_via_residual/self.P.value_si:.6e} m^3")
        print(f"  V_gas = {self.V_gas:.6e} m^3   V_poly = {self.V_poly:.6e} m^3   V = {self.V:.6e} m^3")
        print(f"  gas species: {n_gas_species} (sum y0 = {n_gas_moles:.6e} mol)   "
              f"poly species: {n_poly_species} (sum y0 = {n_poly_moles:.6e} mol)")
        print(f"  {'GAS Species':<30} {'Index':<7} {'y0':>13}")
        for i in range(len(core_species)):
            if self.gas_species_mask[i] and self.y0[i] != 0.0:
                print(f"  {core_species[i].label:<30} {i:<7} {self.y0[i]:>13.4e}")
        # largest-magnitude initial gas concentration (collapse check)
        if self.V_gas > 0:
            cmax = 0.0
            cmax_lbl = "-"
            for i in range(len(core_species)):
                if self.gas_species_mask[i]:
                    c = abs(self.y0[i]) / self.V_gas
                    if c > cmax:
                        cmax = c
                        cmax_lbl = core_species[i].label
            print(f"  max |C_gas0| = {cmax:.4e} mol/m^3  ({cmax_lbl})")

        # --- Mass transfer ---
        if self.mass_transfer:
            print(f"\n--- Mass Transfer ---")
            for mt in self.mass_transfer:
                gas_lbl = core_species[mt.gas_index].label if 0 <= mt.gas_index < len(core_species) else "???"
                poly_lbl = core_species[mt.poly_index].label if 0 <= mt.poly_index < len(core_species) else "???"
                print(f"  {gas_lbl} (gas, idx={mt.gas_index}) <-> {poly_lbl} (poly, idx={mt.poly_index})  "
                      f"kLa={mt.kLa:.2e}")

        print(f"\n{'=' * w}\n")







    def set_initial_conditions(self):
        """
        Sets initial state y0.

        WARNING: The dictionary `initial_mole_fractions` is named by the base class contract,
        but for this reactor type, the values are interpreted as MOLES of gas species.
        If `constant_gas_volume=False`, `V_gas` is inferred from the total gas moles via the
        ideal gas law at (T, P).
        """
        cdef int n_core = self.num_core_species

        ReactionSystem.set_initial_conditions(self)

        if self.gas_species_mask is None or self.gas_species_mask.shape[0] != n_core:
            raise ValueError(f"CRITICAL DIMENSION MISMATCH: gas_species_mask size "
                             f"({None if self.gas_species_mask is None else self.gas_species_mask.shape[0]}) "
                             f"does not match num_core_species ({n_core}).")

        # 1. Gas Species (Interpret values as MOLES)
        for spec, val_moles in self.initial_mole_fractions.items():
            i = self.get_species_index(spec)
            self.y0[i] = val_moles

        # 2. Polymer Explicit Species
        for pool in self.polymer_pools:
            if pool.label in self.initial_explicit_species:
                dp_map = self.initial_explicit_species[pool.label]
                for dp, moles in dp_map.items():
                    if dp in pool.explicit_dp_to_species_index:
                        idx = pool.explicit_dp_to_species_index[dp]
                        self.y0[idx] = max(0.0, moles)
            elif (pool.label not in self.initial_polymer_moments
                    and not getattr(pool, "engine_spawned", False)):
                # An engine-spawned daughter is born empty: its moment slots
                # stay at zero, its declared initial condition.
                logging.warning(f"Polymer pool '{pool.label}' has no initial conditions provided (Explicit or Moments).")

        # 3. Polymer Moments
        for pool in self.polymer_pools:
            if pool.label in self.initial_polymer_moments:
                moms = self.initial_polymer_moments[pool.label]
                for k in range(3):
                    idx = pool.mu_indices[k]
                    self.y0[idx] = max(0.0, moms[k])

        # 4. Initialize Gas Volume
        mask = self.gas_species_mask
        if self.constant_gas_volume:
            if self.V_gas0 is None or self.V_gas0 <= 0.0:
                raise ValueError("constant_gas_volume=True requires V_gas0 > 0")
            self.V_gas = self.V_gas0
        else:
            # Sum moles of gas species only
            n_gas0 = float(np.sum(self.y0[:n_core][mask]))
            self.V_gas = constants.R * self.T.value_si * n_gas0 / self.P.value_si if n_gas0 > 0 else 1.0
        self.V = self.V_gas + self.V_poly

        # 5. Set Concentrations
        for j in range(n_core):
            if mask[j]:
                self.core_species_concentrations[j] = self.y0[j] / self.V_gas
            else:
                self.core_species_concentrations[j] = self.y0[j] / self.V_poly

        # 6. Moment Consistency Check
        tol = 1e-12
        for pool in self.polymer_pools:
            idx_mu0, idx_mu1, idx_mu2 = pool.mu_indices

            mu0_tot = max(0.0, self.y0[idx_mu0]) / self.V_poly
            mu1_tot = max(0.0, self.y0[idx_mu1]) / self.V_poly
            mu2_tot = max(0.0, self.y0[idx_mu2]) / self.V_poly

            mu0_exp, mu1_exp, mu2_exp = _explicit_moment_contributions(
                V_poly=self.V_poly,
                y=self.y0,
                explicit_dp_to_species_index=pool.explicit_dp_to_species_index)

            r0 = mu0_tot - mu0_exp
            r1 = mu1_tot - mu1_exp
            r2 = mu2_tot - mu2_exp

            if r1 < -tol:
                logging.warning(f"Polymer pool {pool.label}: Explicit mass exceeds Total. Clamping tail to 0.")

            self.y0[idx_mu0] = max(0.0, r0) * self.V_poly
            self.y0[idx_mu1] = max(0.0, r1) * self.V_poly
            self.y0[idx_mu2] = max(0.0, r2) * self.V_poly

        # 7. Weak-link U-state slots (milestones ii+iii): U(0) =
        # unsaturated_tail_ends_initial [mol -- SAME amount basis as mu0],
        # behind the CENSUS TRAP (adversarial-review requirement): each
        # chain carries at most 2 tail ends, so U0 must fit the TAIL-
        # DISTRIBUTION chain-end capacity 2*mu0_tail (the tail mu0 written
        # back in step 6). TAIL-ONLY basis (review r37): U is a tail-
        # distribution state -- B and the RHS capacity throttle count the
        # tail moments only -- so the census must use the SAME basis;
        # explicit-species ends must not back a U0 that would then evolve
        # against tail-only capacity (mixed semantics). A mol vs mol/L (or
        # other per-volume) typo would otherwise become a silent hidden
        # initiation source. U0 = 0 always passes (mu0 may be 0 for an
        # empty pool).
        u_slots = getattr(self, "qssa_u_slot", None)
        if u_slots is not None:
            for i, pool in enumerate(self.polymer_pools):
                slot = int(u_slots[i])
                if slot < 0:
                    continue
                u0 = float(self.qssa_u0[i])
                mu0_amount = max(0.0, self.y0[pool.mu_indices[0]])
                capacity = 2.0 * mu0_amount
                if u0 > capacity:
                    raise ValueError(
                        f"Pool {pool.label}: radical_qssa_unzip "
                        f"unsaturated_tail_ends_initial={u0:g} mol exceeds "
                        f"the pool's tail-distribution chain-end capacity "
                        f"2*mu0 = {capacity:g} mol (initial TAIL mu0 = "
                        f"{mu0_amount:g} mol; each chain carries at most 2 "
                        f"tail ends, and U is a TAIL-distribution state -- "
                        f"explicit-species chain ends do not back U, "
                        f"matching the RHS capacity throttle basis). U is "
                        f"on the SAME amount basis as mu0 [mol] -- a mol "
                        f"vs mol/L (or other per-volume) typo here would "
                        f"silently become a hidden initiation source. Fix "
                        f"the amount basis or the value.")
                self.y0[slot] = u0

        # 8. SGH kernel-v2 Z aux-state slots (§1c): Z(0) = sites * mu1(0)
        # [mol], the moles of removable X (Br) sites initially present on the
        # carrier's tail chains. AMOUNT basis, NOT per-volume: mu1(0) is the
        # extensive moment slot value written in step 3/clamped in step 6
        # (the SAME tail-mu1 the condensed-mass tally and the RHS Z-bound
        # guard read back). The §1h build guard has already refused any
        # carrier whose mu1 moves, so this seed stays consistent with
        # sites*mu1(t) for all t. A mol-vs-per-volume typo here would silently
        # mis-seed the entire Br budget, so the amount basis is explicit.
        z_slots = getattr(self, "sgh_z_slot", None)
        if z_slots is not None and len(z_slots) > 0:
            for i, pool in enumerate(self.polymer_pools):
                if not self.sgh_enabled[i]:
                    continue
                mu1_amount = max(0.0, self.y0[pool.mu_indices[1]])
                for sgh_row in range(int(self.sgh_row_start[i]),
                                     int(self.sgh_row_end[i])):
                    z_slot = int(z_slots[sgh_row])
                    if z_slot < 0:
                        continue
                    self.y0[z_slot] = self.sgh_sites[sgh_row] * mu1_amount
                    if self.sgh_in_pool[sgh_row]:
                        # An X-loss daughter's chains each lost one site
                        # already (I-260).
                        self.y0[z_slot] = max(
                            0.0, self.y0[z_slot]
                            - max(0.0, self.y0[pool.mu_indices[0]]))

    def generate_rate_coefficients(self, core_reactions, edge_reactions):
        for rxn in itertools.chain(core_reactions, edge_reactions):
            j = self.reaction_index[rxn]
            kf = rxn.get_rate_coefficient(self.T.value_si, self.P.value_si)
            # Guard against non-finite forward coefficients (e.g. extreme PLOG/Troe
            # extrapolation) so a single bad reaction cannot poison the residual.
            self.kf[j] = kf if np.isfinite(kf) else 0.0
            if rxn.reversible:
                Keq = rxn.get_equilibrium_constant(self.T.value_si)
                # When Keq under/overflows (thermo at extreme T for large fragments),
                # kb = kf/Keq -> inf, and inf * 0 (zero product conc) = NaN, which
                # crashes DASPK ("nans in moles"). Treat such degenerate reverse
                # directions as irreversible rather than emitting a non-finite kb.
                if (not np.isfinite(Keq)) or Keq <= 0.0:
                    self.Keq[j] = np.inf
                    self.kb[j] = 0.0
                else:
                    kb = self.kf[j] / Keq
                    self.Keq[j] = Keq
                    self.kb[j] = kb if np.isfinite(kb) else 0.0
            else:
                self.kb[j] = 0.0
                self.Keq[j] = np.inf

    def get_threshold_rate_constants(self, model_settings):
        """
        Get the theoretical maximum rate constants for reaction filtering.
        These are used to normalize edge reaction fluxes.
        """
        # Unimolecular limit: kB * T / h (Transition State Theory limit)
        unimolecular_threshold_rate_constant = (constants.kB / constants.h) * self.T.value_si

        # Bimolecular limit: User-defined (usually 1e8 - 1e12 m^3/mol/s)
        bimolecular_threshold_rate_constant = model_settings.filter_threshold

        # Trimolecular limit: Adjusted for volume units
        trimolecular_threshold_rate_constant = model_settings.filter_threshold / 1e3

        return (unimolecular_threshold_rate_constant,
                bimolecular_threshold_rate_constant,
                trimolecular_threshold_rate_constant)

    def _recompute_conduit_accepted_flux(self, r_idx):
        """SIDE-EFFECT-FREE recompute of the arch-7 net event flux
        F = r_mol_s [mol/s] at the ACCEPTED state ``self.y``.

        The residual() CANNOT be trusted as the flux source at the accepted
        step: it is a DASPK callback re-invoked on Newton/corrector iterations
        and by the finite-difference Jacobian helper
        :meth:`_scoped_fd_jacobian` (which sets ``_jac_scope=True`` and calls
        perturbed residuals), so after :meth:`step` returns the accepted step
        the last residual evaluation reflects a Jacobian-FD-perturbed or
        rejected-trial state, NOT the accepted solution -- and base.pyx's
        simulate() does not refresh rates before the accepted-step hook, so
        core_reaction_rates is stale for the same reason. We therefore
        RECOMPUTE F from ``self.y`` here, mirroring residual()'s rate
        primitives EXACTLY but writing ONLY to locals (no dn_dt, no
        diagnostics, no shared arrays).

        Mirrors residual() §2-§6 for one row, read-only:
          * V_gas from the accepted gas moles (ideal-gas, residual §2), local;
          * concentrations = softclamp(y)/V with the SAME negative-only
            softclamp extension residual §3 applies, into a LOCAL y_head;
          * rf = kf * prod(_C(reactant slots)) with residual's _C read rules
            (pool proxy -> 1, edge -> 0, else C_gas/C_poly);
          * V_rxn == V_poly: arch-7 admission REQUIRES a condensed chain
            reactant (polymer_conduit.evaluate_conduit_admission G5 chain
            phase, ~L908-917), so the residual's phase selection
            (is_poly_event -> V_rxn = V_poly, ~L5795/5799) is always poly for
            an arch-7 row; and irreversible export (kr == 0) makes
            r_mol_s = rf * V_rxn = rf * V_poly = F. Returns F [mol/s]."""
        n_core = self.num_core_species
        y = self.y
        mask = self.gas_species_mask
        # V_gas at the ACCEPTED state (residual §2), into a LOCAL -- never
        # touch self.V_gas / self.V (those belong to the residual).
        if self.constant_gas_volume:
            V_gas = self.V_gas0
        else:
            n_gas = float(np.sum(y[:n_core][mask]))
            V_gas = (constants.R * self.T.value_si * n_gas / self.P.value_si
                     if n_gas > 0.0 else 1.0)
        V_poly = self.V_poly
        # softclamp exactly as residual §3 (negative entries only), into a
        # LOCAL copy so self.y is never mutated.
        y_head = y[:n_core]
        neg = y_head < 0.0
        if neg.any():
            y_head = y_head.copy()
            yn = y_head[neg]
            y_head[neg] = yn * np.exp(yn / self._softclamp_lam)
        s2p = self.species_to_pool_indices
        ir = self.reactant_indices
        rf = self.kf[r_idx]
        for slot in range(3):
            idx = ir[r_idx, slot]
            if idx == -1:
                continue
            if idx >= n_core:
                # Edge reactant: no state in y -> concentration 0 (residual
                # _C parity), rf vanishes.
                rf = 0.0
                break
            if s2p[idx] != -1:
                continue                 # pool proxy: unit activity (C == 1)
            rf *= (y_head[idx] / V_gas) if mask[idx] else (y_head[idx] / V_poly)
        return rf * V_poly

    def _book_conduit_flux_accepted(self, double dt):
        """M18.4 increment 2 WRITER (DESIGN §4.4), adversarial-review rebuild:
        accumulate, per admitted arch-7 conduit candidate_key, the GRAMS of
        gas emitted over ONE ACCEPTED integrator step of length ``dt`` [s],
        into the run-level ``rmgpy.polymer_conduit._CONDUIT_FLUX_TOTALS`` store
        ``{key: {"grams": float, "revoked": bool}}``.

        Called ONLY from :meth:`_phase_gate_flux_census`, itself invoked once
        per ACCEPTED snapshot from base.pyx's simulate() loop (never from the
        residual). F is recomputed side-effect-freely from the ACCEPTED
        ``self.y`` (:meth:`_recompute_conduit_accepted_flux`) -- NOT read from
        any residual scratch, which after step() reflects a Jacobian-perturbed
        / rejected-trial state.

        TRAPEZOID rule (fix 2): grams += 0.5*(F_prev + F_curr) * dt * gas_mw,
        where F_prev is the accepted-state F of the PREVIOUS accepted step
        (``conduit_prev_accepted_flux``; 0 at the first primed snapshot). The
        right-endpoint rectangle F_curr*dt it replaces is badly biased for
        coarse steps.

        REVOKED handling (DESIGN: "a REVOKED row's accumulated mass stays
        counted -- it happened"): a key already flagged ``revoked`` is NOT
        accumulated further, but its existing grams are left untouched.

        Fix 5 -- CORE rows only: the arch-7 credit state is meaningful for core
        rows; an edge conduit row is neither designed nor tested here and is
        skipped explicitly.

        Fix 4 -- GUARDS: a non-finite/negative F, dt, or grams is skipped
        (logged), never allowed to poison the totals or F_prev.

        Inert in normal runs: with CONDUIT_ADMISSION_ENABLED False no reaction
        is stamped archetype 7, so the loop finds no live conduit row and books
        nothing. Never raises (best-effort diagnostic accumulator)."""
        keys = self.conduit_candidate_key
        if keys is None:
            return
        prev_flux = self.conduit_prev_accepted_flux
        n_core_rxn = self.num_core_reactions
        # dt SEMANTICS: dt == 0 is the DESIGNED prime (first accepted snapshot /
        # rebuild-boundary re-null) -- book nothing but still advance F_prev in
        # the loop below. A non-finite or NEGATIVE dt is anomalous (clock went
        # backward / NaN): log it and return WITHOUT booking AND without
        # advancing the trapezoid endpoints, so a garbage delta neither drops
        # mass silently nor seeds the next step's trapezoid with a poisoned
        # F_prev. Reserve the silent no-book path strictly for dt == 0.
        # (fix 4 / round 84)
        if (not np.isfinite(dt)) or dt < 0.0:
            logging.warning(
                "MOMENT-CREDIT CONDUIT (M18.4 increment 2): anomalous accepted-"
                "step dt %r; skipping booking and NOT advancing trapezoid "
                "endpoints this step.", dt)
            return
        book = dt > 0.0
        totals = None
        if book:
            try:
                from rmgpy.polymer_conduit import _CONDUIT_FLUX_TOTALS as totals
            except Exception as exc:
                # The run-level census is load-bearing: surface the failure
                # loudly rather than silently dropping this step's grams.
                logging.warning(
                    "MOMENT-CREDIT CONDUIT (M18.4 increment 2): could not "
                    "import _CONDUIT_FLUX_TOTALS (%s: %s); accepted-step grams "
                    "NOT booked this step.", type(exc).__name__, exc)
                return
        for r_idx in range(len(keys)):
            if self.reaction_flux_archetype[r_idx] != FLUX_MOMENT_CREDIT_CONDUIT:
                continue
            if r_idx >= n_core_rxn:
                continue  # fix 5: edge arch-7 rows are out of scope here
            F_curr = self._recompute_conduit_accepted_flux(r_idx)
            if not np.isfinite(F_curr) or F_curr < 0.0:
                logging.warning(
                    "MOMENT-CREDIT CONDUIT (M18.4 increment 2): recomputed "
                    "accepted-state flux for core row %d is non-finite/negative "
                    "(%r); skipping (no booking, F_prev zeroed).",
                    r_idx, F_curr)
                prev_flux[r_idx] = 0.0
                continue
            F_prev = prev_flux[r_idx]
            # Advance the trapezoid endpoint for the NEXT accepted step even on
            # a non-booking (dt<=0) priming snapshot.
            prev_flux[r_idx] = F_curr
            if not book:
                continue
            key = keys[r_idx]
            if not key:
                continue
            entry = totals.get(key)
            # REVOKED: keep what was already booked, add nothing more.
            if entry is not None and entry.get("revoked"):
                continue
            grams = 0.5 * (F_prev + F_curr) * dt * self.conduit_gas_mw[r_idx]
            if not np.isfinite(grams) or grams < 0.0:
                logging.warning(
                    "MOMENT-CREDIT CONDUIT (M18.4 increment 2): computed grams "
                    "for key %r is non-finite/negative (%r); skipping to avoid "
                    "poisoning the totals.", key, grams)
                continue
            if entry is None:
                totals[key] = {"grams": grams, "revoked": False}
            else:
                # Guard the PRE-EXISTING accumulator too (round 84): a total
                # already poisoned (non-finite/negative, from any source) must
                # not be silently propagated -- the artifact emitter floats it
                # blindly (rmgpy/polymer.py). Validate BOTH prior and post-add;
                # leave the entry untouched and log if either is bad.
                prior = entry.get("grams", 0.0)
                new_total = prior + grams
                if (not np.isfinite(prior) or prior < 0.0
                        or not np.isfinite(new_total) or new_total < 0.0):
                    logging.warning(
                        "MOMENT-CREDIT CONDUIT (M18.4 increment 2): pre-existing "
                        "or updated total for key %r is non-finite/negative "
                        "(prior=%r, new=%r); leaving totals unchanged.",
                        key, prior, new_total)
                    continue
                entry["grams"] = new_total

    def _phase_gate_flux_census(self, core_species, edge_species,
                                edge_reactions, char_rate, tol_move_to_core):
        """RIDER R2 dynamic half (item 17, spec 2026-06-12 SS3(e)): emit one
        census line per (gate, edge reaction) per engine rebuild when the
        UNGATED ratio would have cleared the enlargement bar -- the
        species-level quantity base.pyx:1046 actually tests. Lives in
        simulate() via the base hook because tol_move_to_core is a simulate
        local and char_rate exists only per accepted snapshot; reads the
        most recent residual evaluation's arrays with EXACTLY the staleness
        of the enlargement read it audits (amendment A2 -- a feature).
        String assembly happens here (python level, once per emission),
        never in the residual."""
        # SGH kernel-v2 ACCEPTED-state Z-inventory invariant (§1g, P1-A):
        # base.pyx invokes this override once per ACCEPTED snapshot with
        # self.y holding the integrator's accepted solution (the stale-array
        # note above governs the ungated-rate ARRAYS this census audits, not
        # self.y), so this is the right context for the HARD physical-band
        # check the residual guard (trial state) deliberately does not make.
        # Self-guards on num_sgh_z == 0, so non-SGH runs pay ~nothing.
        self._assert_sgh_inventory_accepted()
        # P1-3: pool-moment negative tripwire, same accepted-state hook --
        # the r81 hard-error semantics live here now (the residual's
        # trial-state census never raises).
        self._assert_pool_moments_accepted()
        # M18.4 increment 2: book admitted-conduit gas GRAMS for THIS accepted
        # step (DESIGN §4.4). Placed on the accepted-step hook (never the
        # residual) so REJECTED Newton trials cannot overcount, and BEFORE the
        # char_rate==0 return so a no-enlargement step still books its flux.
        # dt is the accepted-step length t_now - t_prev; the first primed
        # snapshot books dt=0, and a rebuild-boundary re-null makes dt<=0 (no
        # booking) rather than a stale/negative delta.
        conduit_dt = 0.0
        if self._conduit_prev_accepted_t is not None:
            conduit_dt = self.t - self._conduit_prev_accepted_t
        self._conduit_prev_accepted_t = self.t
        self._book_conduit_flux_accepted(conduit_dt)
        if char_rate == 0.0:
            return  # the base.pyx singularity path owns the no-flux case
        gate_codes = getattr(self, "edge_reaction_gate_code", None)
        if gate_codes is None:
            return
        n_core = self.num_core_species
        n_core_rxn = self.num_core_reactions
        ip = self.product_indices
        for j in range(gate_codes.shape[0]):
            code = int(gate_codes[j])
            if code == 0:
                continue
            key = (code, j)
            if key in self._phase_gate_census_emitted:
                continue
            best_ratio = 0.0
            best_label = "<none>"
            condensed_labels = []
            for slot in range(3):
                p = ip[n_core_rxn + j, slot]
                if p == -1:
                    continue
                if not self.prospective_gas_mask[p]:
                    condensed_labels.append(
                        edge_species[p - n_core].label if p >= n_core
                        else core_species[p].label)
                if p >= n_core:
                    ratio = abs(self.edge_species_rates_ungated[p - n_core]) \
                        / char_rate
                    if ratio > best_ratio:
                        best_ratio = ratio
                        best_label = edge_species[p - n_core].label
            if best_ratio <= tol_move_to_core:
                continue
            self._phase_gate_census_emitted.add(key)
            rxn = edge_reactions[j]
            idx = n_core_rxn + j
            decisive = (", ".join(condensed_labels) if code == 1
                        else "no prospectively-condensed product")
            fam = (getattr(rxn, "family", None)
                   or getattr(rxn, "label", "") or type(rxn).__name__)
            logging.warning(
                "PHASE-GATE FLUX CENSUS: reaction %s gate=%s "
                "ungated_ratio=%.6e via edge species %s "
                "(tol_move_to_core=%.3e); decisive=%s; archetype "
                "pre-demotion=%d post-demotion=%d family=%s",
                rxn, "A" if code == 1 else "B", best_ratio, best_label,
                tol_move_to_core, decisive,
                int(self.reaction_pre_demotion_archetype[idx]),
                int(self.reaction_flux_archetype[idx]), fam)

    def _static_phase_gate_census(self, core_species, core_reactions):
        """RIDER R2 static half (item 17, spec 2026-06-12 SS3(e), amendment
        A1): enumerate CORE reactions whose phase-gate verdict is zero and
        announce each once per rebuild. Covers reactions that arrive
        core-gated without crossing the edge on their own flux: the third
        route (independent species promotion -- each participant promoted on
        other channels' flux) and legacy/restart cores. Gate verdicts for
        core rows are STATIC -- masks + event type only, no rates -- so this
        runs at init with zero residual cost; the kinetics loop's bare
        continue for core-gated rows survives untouched. The gate
        classification below MUST mirror residual section 6 (the
        prospective-mask gates) -- keep in sync."""
        n_core = self.num_core_species
        ir = self.reactant_indices
        ip = self.product_indices
        mask = self.gas_species_mask
        pmask = self.prospective_gas_mask
        for i in range(len(core_reactions)):
            r0, r1, r2 = ir[i, 0], ir[i, 1], ir[i, 2]
            # Defensive parity with the residual's reactant skips (core
            # reactions have all-core participants by construction).
            if r0 == -1 or r0 >= n_core:
                continue
            if (r1 != -1 and r1 >= n_core) or (r2 != -1 and r2 >= n_core):
                continue
            is_poly_event = ((not mask[r0])
                             or (r1 != -1 and not mask[r1])
                             or (r2 != -1 and not mask[r2]))
            has_condensed_prod = False
            condensed_labels = []
            for slot in range(3):
                p = ip[i, slot]
                if p == -1:
                    continue
                if not pmask[p]:
                    has_condensed_prod = True
                    if p < n_core:
                        condensed_labels.append(core_species[p].label)
            if (not is_poly_event) and has_condensed_prod:
                code = 1
            elif is_poly_event and not has_condensed_prod:
                code = 2
            else:
                continue
            rxn = core_reactions[i]
            decisive = (", ".join(condensed_labels) if code == 1
                        else "no prospectively-condensed product")
            fam = (getattr(rxn, "family", None)
                   or getattr(rxn, "label", "") or type(rxn).__name__)
            logging.warning(
                "PHASE-GATE FLUX CENSUS: reaction %s gate=%s static (core, "
                "init-time); decisive=%s; archetype pre-demotion=%d "
                "post-demotion=%d family=%s",
                rxn, "A" if code == 1 else "B", decisive,
                int(self.reaction_pre_demotion_archetype[i]),
                int(self.reaction_flux_archetype[i]), fam)

    def spawn_gate_flux_snapshot(self, motif_counts_by_pool=None):
        """Engine half of the mass-flux spawn gate.

        Spec: docs/superpowers/specs/2026-06-10-mass-flux-spawn-gate-design.md
        §3/§4.1 (AMENDED). The engine cannot attribute representatives (it
        has no ledger), so the split is: the engine reads arrays, the
        python-side gate does ledger-dependent attribution. Returns a
        3-tuple ``(gross, pool_stats, proxy_event_mass_total)``:

        1. ``gross``: dict core-species label ->
           ``max(0, core_species_production_rates[i])`` for ALL core species
           (ordinary species have real entries since change (a), spec §4.6).
           Labels are ``Species.label`` verbatim (the same labels the ledger
           records at absorption). Duplicate core labels RAISE ValueError —
           uniqueness is NOT enforced on the standard non-RMS path
           (model.py's dedup is gated on edge.phase_system), and silent
           overwrite would misattribute gate flux; the main.py caller turns
           the raise into warn + None snapshot -> the gate defers.
        2. ``pool_stats``: dict pool label ->
           ``(E_n, monomer_mw_g_mol, chain_mass_defect_g_mol)``
           with ``E_n = y[mu1]/y[mu0] if y[mu0] > trust else 0.0`` where
           ``trust = max(SMALL_EPS, ATTRIBUTION_TRUST_K * atol_mu0)`` is the
           ATTRIBUTION TRUST FLOOR (item #14a): mu0 at or below the
           integrator's own noise floor cannot substantiate a mass
           attribution, so the band (SMALL_EPS, trust] zeroes E_n (the band
           defers) and logs a ``SPAWN-GATE ATTRIBUTION CENSUS`` warning once
           per snapshot per pool. atol_mu0 is the pool's mu0-slot entry when
           the tolerance is a vector, the scalar atol otherwise. mu0
           exhausted (<= SMALL_EPS) still defers SILENTLY (unchanged). The
           solver's own E[n] consumers and the residual keep SMALL_EPS
           realizability semantics — errs toward deferral either way.
           ``motif_counts_by_pool`` (optional dict pool label -> int) lets
           the main.py stash report in the census line how many ledger
           motifs currently attribute to the zeroed pool (the engine has no
           ledger; the caller does).
        3. ``proxy_event_mass_total``: float — sum of
           ``gross_j * max(0, condensed_mass_g(1, E_n[pool(j)]))``
           = ``gross_j * max(0, E_n*mw - chain_mass_defect)`` over
           CANONICAL pool proxies
           (``species_to_pool_indices[j] != -1 and is_pool_proxy[j]``
           — the engine CAN attribute those). The defect term is the
           FR1-K2 mass-consumer audit (round-72 P2): an X-loss feature
           pool's chains each lost exactly one X, so the EXACT per-chain
           event mass is the pool's normative condensed_mass_g at
           (mu0=1, mu1=E_n) — the raw E_n*mw would overstate every event
           by M_X. Ordinary pools (defect 0) reduce to the legacy
           product bit-identically; the max(0, .) clamp errs toward
           deferral.

        GROSS production, never net dn_dt: canonical proxies have
        dn_dt ~= 0 BY DESIGN (the archetype apportionment reroutes their
        flux to pool moments) and ordinary species net to ~0 at steady
        state. E[n] is read LIVE from the state vector (never
        recorded-and-stale). All terms are polymer-phase volumetric, so
        V_poly cancels in any fraction of these numbers. Pure read of
        already-maintained state — no bookkeeping here beyond change (a)'s
        residual writes.
        """
        gross = {}
        pool_stats = {}
        proxy_event_mass_total = 0.0
        stp = getattr(self, "species_to_pool_indices", None)
        prod = getattr(self, "core_species_production_rates", None)
        y = getattr(self, "y", None)
        if stp is None or prod is None or y is None:
            return gross, pool_stats, proxy_event_mass_total
        n_core = self.num_core_species
        # Live per-pool stats, gated by the ATTRIBUTION TRUST FLOOR
        # (item #14a, spec 2026-06-11 §2/§3): per-pool mu0-slot atol when
        # the tolerance is a vector (the floor is PER-POOL — a deck mixing
        # large and trace pools gets per-pool trust, a feature), scalar atol
        # otherwise. atol_array unavailable (snapshot before
        # initialize_model) -> atol_mu0 = 0 -> the floor degenerates to
        # SMALL_EPS (pre-floor behavior, honest).
        n_pools = len(self.polymer_pools)
        e_n_by_pool = [0.0] * n_pools
        # i061: the attribution trust band is a CHEMISTRY noise floor and must
        # stay anchored to the deck atol, so read the PRE-floor copy. The moment
        # error-weight floor (initialize_model) raised atol_array on the mu-slots
        # for DASPK error control only. Falls back to atol_array when no floor
        # was applied (e.g. a snapshot before initialize_model, atol_mu0->0.0,
        # the pre-existing honest degeneracy).
        atol_arr = getattr(self, "_chem_atol_array", None)
        if atol_arr is None:
            atol_arr = getattr(self, "atol_array", None)
        for p in range(n_pools):
            i0 = self.pool_mu0_indices[p]
            i1 = self.pool_mu1_indices[p]
            if i0 < 0 or i1 < 0:
                continue
            mu0 = y[i0]
            mu1 = y[i1]
            if atol_arr is None:
                atol_mu0 = 0.0
            elif np.ndim(atol_arr) == 0:
                atol_mu0 = float(atol_arr)
            elif i0 < len(atol_arr):
                atol_mu0 = float(atol_arr[i0])
            else:
                atol_mu0 = float(np.max(atol_arr))
            trust = max(SMALL_EPS, ATTRIBUTION_TRUST_K * atol_mu0)
            # r81 (B) negative-moment rule: the pre-r81 mu0 <= SMALL_EPS
            # branch zeroed SILENTLY and would swallow a negative. Now a
            # mu0 beyond -floor is a hard error (integrator corruption);
            # inside [-floor, SMALL_EPS] it projects to E[n] = 0 WITH the
            # exhaustion census below. Floor: the accepted-state floor the
            # r81 check reads (I-244: relative to the pool's scale, see
            # EXHAUSTION_RELATIVE_K) when available (post-initialize_model);
            # SMALL_EPS otherwise (pre-init snapshot -> honest degenerate,
            # mirroring atol_mu0 = 0).
            floors_exh = self.pool_moment_accepted_floors()
            floor0 = (float(floors_exh[p, 0]) if floors_exh is not None
                      and p < len(floors_exh) else SMALL_EPS)
            if mu0 < -floor0:
                raise ValueError(
                    f"Polymer pool '{self.polymer_pools[p].label}' mu0="
                    f"{mu0!r} mol is negative beyond the exhaustion floor "
                    f"-{floor0:.6e} mol (= -max(SMALL_EPS, "
                    f"{EXHAUSTION_FLOOR_K:.0f}*atol[state], "
                    f"{EXHAUSTION_RELATIVE_K:g}*rtol*peak|mu0|)) at the "
                    f"spawn-gate snapshot. A negative moment beyond the "
                    f"integrator's own error budget is integrator "
                    f"corruption, not exhaustion -- refusing to attribute "
                    f"from it (r81 negative-moment rule: hard error, never "
                    f"a silent 0.0).")
            if mu0 > trust:
                e_n_by_pool[p] = mu1 / mu0
            elif mu0 > SMALL_EPS:
                # The distrust band (SMALL_EPS, trust]: integrator-noise-
                # scale mu0. Zero-out is the honest measurement ("this pool
                # cannot substantiate any mass attribution this iteration"),
                # drags the window statistic down, and preserves the
                # subsystem asymmetry: under-attribution defers, never
                # falsely spawns. The census line is the greppable reason a
                # deck's spawns mysteriously defer AND the standing sensor
                # for the floor ever biting legitimate chemistry (the
                # rejected-absolute-constant failure mode, made observable).
                n_motifs = 0
                if motif_counts_by_pool:
                    n_motifs = int(motif_counts_by_pool.get(
                        self.polymer_pools[p].label, 0))
                logging.warning(
                    "SPAWN-GATE ATTRIBUTION CENSUS: pool %s mu0=%.6e is in "
                    "the distrust band (SMALL_EPS=%.0e < mu0 <= trust "
                    "floor=%.6e = max(SMALL_EPS, %.0f * atol_mu0=%.6e)); "
                    "E[n] attribution zeroed for this snapshot; ledger "
                    "motifs attributed to this pool: %d.",
                    self.polymer_pools[p].label, mu0, SMALL_EPS, trust,
                    ATTRIBUTION_TRUST_K, atol_mu0, n_motifs)
            else:
                # mu0 in [-floor0, SMALL_EPS]: exhausted/empty -> E[n] = 0.0
                # WITH the exhaustion census (r81 B5, once per pool per
                # rebuild) -- never silently: a swallowed in-band negative
                # would hide the onset of integrator corruption.
                key_exh = (self.polymer_pools[p].label, "spawn-gate")
                if key_exh not in self._exhaustion_census_emitted:
                    self._exhaustion_census_emitted.add(key_exh)
                    logging.warning(
                        "POOL EXHAUSTION CENSUS: pool %s spawn-gate "
                        "mu0=%.6e mol is within the exhaustion band "
                        "[-%.6e, %.0e]; E[n] attribution projected to 0.0 "
                        "for this snapshot (floor = max(SMALL_EPS, "
                        "%.0f*atol[mu0]); negative beyond -floor is a "
                        "hard error).",
                        self.polymer_pools[p].label, mu0, floor0, SMALL_EPS,
                        EXHAUSTION_FLOOR_K)
        for p in range(n_pools):
            pool = self.polymer_pools[p]
            mw = float(getattr(pool, "monomer_mw_g_mol", 0.0) or 0.0)
            defect = float(getattr(pool, "chain_mass_defect_g_mol", 0.0)
                           or 0.0)
            pool_stats[pool.label] = (e_n_by_pool[p], mw, defect)
        # index -> label for CORE species (species_index covers core+edge).
        labels = {}
        for spc, idx in self.species_index.items():
            if idx < n_core:
                labels[idx] = spc.label
        for i in range(n_core):
            label = labels.get(i)
            if label is None:
                continue
            g = max(0.0, float(prod[i]))
            if label in gross:
                # Label uniqueness is NOT an enforced invariant on the
                # standard (non-RMS) path: model.py's dedup loop is gated on
                # edge.phase_system ("!!! Not maintained when
                # require_rms=False?"). A duplicate key would silently
                # misattribute gate flux. Raise instead: the main.py stash
                # wraps this call, logs a warning and leaves the snapshot
                # None -> the gate defers honestly (spec §4.5).
                raise ValueError(
                    f"duplicate core species label {label!r} in "
                    f"spawn_gate_flux_snapshot(); spawn-gate attribution "
                    f"would be unreliable")
            gross[label] = g
            p = stp[i]
            if p >= 0 and self.is_pool_proxy[i]:
                pool = self.polymer_pools[p]
                e_n = pool_stats[pool.label][0]
                # EXACT per-chain event mass (FR1-K2 mass audit, round-72
                # P2): the pool's normative condensed_mass_g at
                # (mu0=1, mu1=E_n) = E_n*mw - chain_mass_defect. Ordinary
                # pools (defect 0) reduce to the legacy E_n*mw
                # bit-identically; max(0, .) errs toward deferral.
                # SGH kernel-v2 (§1e.2): an SGH carrier's Br leaves through
                # the Z inventory, NOT the chain_mass_defect (which v2 pins to
                # 0.0), so this per-chain CHAIN-population event mass is inert
                # to SGH and reduces to E_n*mw. Assert that invariant here so
                # a future defect-bearing SGH carrier can't silently
                # double-book Br into both Z and the spawn-gate snapshot;
                # tied to the §1h stationary-mu1 contract (the FR1 carrier is
                # not a spawn-gate proxy on the reference deck, so this branch
                # is not even reached for it).
                sgh_en2 = getattr(self, "sgh_enabled", None)
                # An X-loss daughter (I-260) is not such a carrier: its
                # defect is the first loss per chain, its Z the rest.
                if (sgh_en2 is not None and p < len(sgh_en2) and sgh_en2[p]
                        and not self.sgh_in_pool[self.sgh_row_start[p]]
                        and abs(float(pool.chain_mass_defect_g_mol)) > 0.0):
                    raise ValueError(
                        f"SGH-v2 carrier pool '{pool.label}' is a spawn-gate "
                        f"proxy with a nonzero chain_mass_defect_g_mol="
                        f"{pool.chain_mass_defect_g_mol!r}; v2 carries Br "
                        f"through the Z inventory and pins the carrier defect "
                        f"to 0.0, so a nonzero defect would double-book Br "
                        f"into both Z and the spawn-gate event mass.")
                proxy_event_mass_total += g * max(
                    0.0, pool.condensed_mass_g(1.0, e_n))
        return gross, pool_stats, proxy_event_mass_total

    def _update_pool_exhaustion(self, y):
        """Centralized pool-exhaustion CENSUS/TRIPWIRE driver (r81 B, the
        PP run-7 IDID=-7 wall; re-adjudicated r86: the r81 RHS
        read-projection is reverted -- ``self._pool_exhausted`` is kept
        for DIAGNOSTIC LOGGING ONLY and has ZERO RHS effect; the RHS
        always computes from the raw state vector). Called once per
        residual evaluation; refreshes ``self._pool_exhausted`` from the
        RAW state vector ``y``:

        * a pool is counted EXHAUSTED (census once per pool per rebuild)
          only when ``abs(mu0)``, ``abs(mu1)`` AND ``abs(mu2)`` are ALL
          below their per-state floors
          ``max(SMALL_EPS, EXHAUSTION_FLOOR_K * atol[state])`` (never mu0
          alone: tiny mu0 with nontrivial mu1 is a few long chains);
        * any moment more negative than ``-floor_k`` on this TRIAL state
          is CENSUSED (worst excursion recorded per pool, warn-once) --
          never raised here (P1-3): residual() runs on raw Newton TRIAL
          states and pydas cannot propagate the exception (SystemError ->
          "ITERATION MATRIX IS SINGULAR" -> dead run). The r81 hard-error
          tripwire (integrator corruption, never ``max(..., SMALL_EPS)``)
          moved verbatim to ``_assert_pool_moments_accepted()`` on the
          ACCEPTED state (the SGH accepted-state-hook precedent);
        * a negative inside ``[-floor_k, +floor_k]`` is ANNOUNCED through
          the exhaustion census instead of silence (the RHS's
          pre-existing local ``max(0.0, .)`` reads are what clamp it).

        ``y`` itself is never mutated. mu3 < mu2 OUTSIDE the exhaustion
        band deliberately stays a real guard event (the existing warn-once
        kernel guards are untouched)."""
        cdef int p, n_pools
        cdef double f0, f1, f2, raw0, raw1, raw2
        n_pools = len(self.polymer_pools)
        for p in range(n_pools):
            idx0, idx1, idx2 = self.polymer_pools[p].mu_indices
            f0 = self._pool_mu_floors[p, 0]
            f1 = self._pool_mu_floors[p, 1]
            f2 = self._pool_mu_floors[p, 2]
            raw0 = y[idx0]
            raw1 = y[idx1]
            raw2 = y[idx2]
            if raw0 < -f0 or raw1 < -f1 or raw2 < -f2:
                # TRIAL-STATE-SAFE (P1-3, regen-#2 forensics): residual() is
                # evaluated on RAW Newton TRIAL states -- pydas cannot
                # propagate a Python exception out of the corrector loop, so
                # a raise here surfaces as SystemError -> "ITERATION MATRIX
                # IS SINGULAR" -> dead run even when the trial would simply
                # have been rejected. Census only: record the worst
                # excursion per pool and warn once per rebuild. The r81
                # hard-error semantics (negative moment beyond
                # -max(SMALL_EPS, 100*atol) on an ACCEPTED state =
                # integrator corruption, hard stop, never
                # max(..., SMALL_EPS)) live verbatim in
                # _assert_pool_moments_accepted(), invoked per accepted
                # snapshot from _phase_gate_flux_census beside the SGH
                # accepted-state hook (the SGH precedent).
                worst = min(raw0, raw1, raw2)
                if (self._pool_worst_trial_excursion is not None
                        and worst < self._pool_worst_trial_excursion[p]):
                    self._pool_worst_trial_excursion[p] = worst
                self._emit_pool_exhaustion_census(
                    p, raw0, raw1, raw2,
                    "beyond-floor negative moment on a TRIAL state (census "
                    "only -- accepted-state hard check owns the r81 raise)")
                self._pool_exhausted[p] = 0
                continue
            if (abs(raw0) <= f0 and abs(raw1) <= f1 and abs(raw2) <= f2):
                self._pool_exhausted[p] = 1
                # A pool born empty (initial moments exactly zero) that has
                # not filled past its floor has not been drained: not
                # announced. A pool that starts inside the band is.
                if self._pool_held(p) or not self._pool_born_empty(p):
                    self._emit_pool_exhaustion_census(
                        p, raw0, raw1, raw2,
                        "exhausted (census only -- r86: RHS computes live "
                        "from the raw state)")
            else:
                self._pool_exhausted[p] = 0
                if raw0 < 0.0 or raw1 < 0.0 or raw2 < 0.0:
                    self._emit_pool_exhaustion_census(
                        p, raw0, raw1, raw2,
                        "in-band negative moment (census only; pool NOT "
                        "exhausted)")

    def _emit_pool_exhaustion_census(self, int p, double raw0, double raw1,
                                     double raw2, reason):
        """Once-per-pool-per-rebuild exhaustion census (r81 B3,
        re-adjudicated r86: log-only diagnostics -- nothing is suppressed):
        label, raw moments, floors, and which kernels / pool-coupled
        reaction rows draw on the pool while it sits in the census band.
        Cold path -- string and row-scan work happen at most once per
        (pool, reason) per rebuild."""
        pool = self.polymer_pools[p]
        key = (pool.label, reason)
        if key in self._exhaustion_census_emitted:
            return
        self._exhaustion_census_emitted.add(key)
        f0 = self._pool_mu_floors[p, 0]
        f1 = self._pool_mu_floors[p, 1]
        f2 = self._pool_mu_floors[p, 2]
        kernels = []
        if getattr(pool, "tail_kinetics", None):
            kernels.append("tail_kinetics")
        if float(getattr(pool, "k_scission", 0.0) or 0.0) > 0.0:
            kernels.append("k_scission")
        if float(getattr(pool, "k_unzip", 0.0) or 0.0) > 0.0:
            kernels.append("k_unzip")
        khom = getattr(self, "khom_enabled", None)
        if khom is not None and len(khom) > p and khom[p]:
            kernels.append("k_homolysis")
        sgh = getattr(self, "sgh_enabled", None)
        if sgh is not None and len(sgh) > p and sgh[p]:
            kernels.append("side_group_homolysis")
        qssa = getattr(self, "qssa_enabled", None)
        if qssa is not None and len(qssa) > p and qssa[p]:
            kernels.append("radical_qssa_unzip")
        kdep = getattr(self, "kdep_enabled", None)
        if kdep is not None and len(kdep) > p and kdep[p]:
            kernels.append("k_depropagation")
        coupled_rows = []
        n_rows = 0
        ir = getattr(self, "reactant_indices", None)
        ip = getattr(self, "product_indices", None)
        stp = getattr(self, "species_to_pool_indices", None)
        if ir is not None and ip is not None and stp is not None:
            n_core = self.num_core_species
            n_rows = ir.shape[0]
            for j in range(n_rows):
                hit = False
                for slot in range(3):
                    a = ir[j, slot]
                    b = ip[j, slot]
                    if ((a != -1 and a < n_core and stp[a] == p)
                            or (b != -1 and b < n_core and stp[b] == p)):
                        hit = True
                        break
                if hit:
                    coupled_rows.append(j)
        logging.warning(
            "POOL EXHAUSTION CENSUS: pool %s %s (raw mu0=%.6e, mu1=%.6e, "
            "mu2=%.6e mol; floors %.6e/%.6e/%.6e mol = max(SMALL_EPS, "
            "%.0f*atol[state])); pool kernels drawing on this state: %s; "
            "pool-coupled reaction rows: %d of %d (core+edge)%s.",
            pool.label, reason, raw0, raw1, raw2, f0, f1, f2,
            EXHAUSTION_FLOOR_K, kernels if kernels else "none",
            len(coupled_rows), n_rows,
            " (first rows: %s)" % coupled_rows[:8]
            if coupled_rows else "")

    def _assert_pool_moments_accepted(self):
        """Pool-moment ACCEPTED-state negative tripwire (P1-3): the r81
        hard-error semantics moved verbatim off the residual's TRIAL
        states. Reads the trial-unsafe ``self.y`` -- base.pyx invokes
        ``_phase_gate_flux_census`` once per ACCEPTED snapshot with
        ``self.y`` holding the integrator's accepted solution (the SGH
        ``_assert_sgh_inventory_accepted`` precedent). Any moment more
        negative than its accepted floor
        ``-max(SMALL_EPS, EXHAUSTION_FLOOR_K*atol[state],
        EXHAUSTION_RELATIVE_K*rtol[state]*peak|mu_k|)`` (I-244,
        ``pool_moment_accepted_floors``) on an ACCEPTED state is integrator
        corruption, not exhaustion: HARD stop, never ``max(..., SMALL_EPS)``.
        A pool that empties lands inside the floor and is ended by
        ``end_exhausted_pools`` first."""
        cdef int p, n_pools
        cdef double f0, f1, f2, raw0, raw1, raw2
        if self.polymer_pools is None or self._pool_mu_floors is None:
            return
        y = self.y
        accepted = self.pool_moment_accepted_floors()
        n_pools = len(self.polymer_pools)
        for p in range(n_pools):
            idx0, idx1, idx2 = self.polymer_pools[p].mu_indices
            f0 = self._pool_mu_floors[p, 0]
            f1 = self._pool_mu_floors[p, 1]
            f2 = self._pool_mu_floors[p, 2]
            raw0 = y[idx0]
            raw1 = y[idx1]
            raw2 = y[idx2]
            # CONE census on the ACCEPTED state. The negative check below
            # fires on mu1 < 0, which is DOWNSTREAM of the actual violation:
            # the realizable set is mu1 >= mu0 >= 0 (a k>=1 distribution has
            # at least one repeat unit per chain), and a pool exits it at
            # mu1 = mu0, strictly before mu1 = 0. Measured on the legacy
            # k_unzip kernel against a scission-fed pool, the two are ~10 ms
            # apart at k_unzip/k_scission = 100 -- but the mu1 = 0 crossing
            # exists ONLY above k_unzip = k_scission/4 (there the (mu0, mu1)
            # subsystem's eigenvalues go complex and the trajectory spirals
            # through zero), whereas the cone exit happens at EVERY positive
            # k_unzip. Below that ratio a pool could leave the realizable
            # cone and no accepted-state check would ever fire.
            # Deliberately a warn-once CENSUS, not a raise: leaving the cone
            # says a moment SOURCE TERM is wrong, which is the diagnosis the
            # (default-off) debug_check_realizability warning already
            # carries, whereas the r81 raise below is scoped to integrator
            # corruption. A pool that fully unzips legitimately ENDS at the
            # boundary mu1 = mu0, so promoting this to a hard error would
            # kill runs on boundary ULP noise; that is a policy change to
            # adjudicate on its own, not a rider on a kernel fix. The
            # tolerance is the SUM of the two moments' own exhaustion
            # floors -- each is trusted to within its floor, so their
            # difference is trusted to within both. No existing floor,
            # tolerance or check is altered.
            if (raw1 < raw0 - (f0 + f1)
                    and self.polymer_pools[p].label
                    not in self._realizability_warned):
                pool = self.polymer_pools[p]
                self._realizability_warned.add(pool.label)
                logging.warning(
                    "POOL CONE CENSUS: pool %s left the realizable cone on an "
                    "ACCEPTED state: mu1=%.6e < mu0=%.6e mol (gap %.6e, "
                    "tolerance %.6e = f_mu0 + f_mu1). Require mu1 >= mu0 >= 0 "
                    "-- every chain carries at least one repeat unit. This is "
                    "a moment SOURCE TERM defect, not exhaustion, and it "
                    "precedes any negative moment; census only, the r81 "
                    "negative check below is unchanged.",
                    pool.label, raw1, raw0, raw1 - raw0, f0 + f1)

            # VARIANCE census on the ACCEPTED state -- the OTHER half of
            # three-moment realizability, which nothing checked (I-065
            # defect 2). mu1 >= mu0 >= 0 above is necessary but not
            # sufficient: any distribution also satisfies Cauchy-Schwarz,
            #     mu0*mu2 >= mu1^2   (equivalently Var[n] = mu2/mu0 -
            #     (mu1/mu0)^2 >= 0, i.e. a NON-NEGATIVE chain-length
            #     variance), with equality only for a monodisperse pool.
            # A state that violates it has an imaginary spread: PDI < 1,
            # _gamma_params_from_mu012 returns no params, and the mu3
            # log-Lagrange closure mu0*(mu2/mu1)^3 stops describing any
            # distribution at all -- the same class of corruption the cone
            # exit is, and just as invisible to the mu1 < 0 raise below.
            # Tolerance is first-order propagation of the moments' own
            # exhaustion floors through the product, plus the two
            # second-order terms so an all-at-floor state cannot trip it:
            #     d(mu0*mu2 - mu1^2) <= f0*|mu2| + f2*|mu0| + 2*f1*|mu1|
            #                           + f0*f2 + f1*f1
            # No existing floor, tolerance or check is altered, and the
            # same warn-once/never-raise policy as the cone census applies
            # -- promoting either to a hard error is a policy change to
            # adjudicate on its own, not a rider on a kernel fix. Keyed
            # separately in the warn set so a pool that already reported a
            # cone exit can still report a variance violation.
            if (raw0 * raw2 < raw1 * raw1 - (f0 * fabs(raw2)
                                             + f2 * fabs(raw0)
                                             + 2.0 * f1 * fabs(raw1)
                                             + f0 * f2 + f1 * f1)
                    and (self.polymer_pools[p].label, "variance")
                    not in self._realizability_warned):
                pool = self.polymer_pools[p]
                self._realizability_warned.add((pool.label, "variance"))
                logging.warning(
                    "POOL VARIANCE CENSUS: pool %s violated moment "
                    "realizability on an ACCEPTED state: mu0*mu2=%.6e < "
                    "mu1^2=%.6e (mu0=%.6e, mu1=%.6e, mu2=%.6e mol, deficit "
                    "%.6e, tolerance %.6e). Cauchy-Schwarz requires "
                    "mu0*mu2 >= mu1^2 for ANY chain-length distribution; "
                    "below it the chain-length variance is negative, the "
                    "gamma closure has no parameters and the mu3 "
                    "log-Lagrange closure describes nothing. This is a "
                    "moment SOURCE TERM defect; census only, the r81 "
                    "negative check below is unchanged.",
                    pool.label, raw0 * raw2, raw1 * raw1, raw0, raw1, raw2,
                    raw0 * raw2 - raw1 * raw1,
                    (f0 * fabs(raw2) + f2 * fabs(raw0) + 2.0 * f1 * fabs(raw1)
                     + f0 * f2 + f1 * f1))

            # I-244: the hard error reads the accepted floor, relative to the
            # pool's scale (EXHAUSTION_RELATIVE_K); the censuses above keep
            # the absolute floors.
            f0 = accepted[p, 0]
            f1 = accepted[p, 1]
            f2 = accepted[p, 2]
            if raw0 < -f0 or raw1 < -f1 or raw2 < -f2:
                pool = self.polymer_pools[p]
                raise ValueError(
                    f"Polymer pool '{pool.label}' moment state went "
                    f"negative beyond the exhaustion floor on the "
                    f"ACCEPTED state: mu0={raw0!r}, "
                    f"mu1={raw1!r}, mu2={raw2!r} mol vs floors "
                    f"-{f0:.6e}/-{f1:.6e}/-{f2:.6e} mol "
                    f"(= -max(SMALL_EPS, {EXHAUSTION_FLOOR_K:.0f}*"
                    f"atol[state], {EXHAUSTION_RELATIVE_K:g}*rtol*peak"
                    f"|mu_k|)). A negative moment beyond the "
                    f"integrator's own error budget is integrator "
                    f"corruption, not exhaustion -- refusing to integrate "
                    f"it (r81 negative-moment rule: hard error, never "
                    f"max(..., SMALL_EPS)).")

    def pool_moment_accepted_floors(self):
        """(n_pools, 3) accepted-state floors [mol] (I-244, see
        EXHAUSTION_RELATIVE_K): per pool and moment order k,
        ``max(_pool_mu_floors[p, k], EXHAUSTION_RELATIVE_K*rtol[slot]*
        peak|mu_k|)``, the peak over the accepted states so far (the initial
        state's |mu_k| for a pool that only drains; the absolute floor alone
        for a pool born empty until it fills). None before
        initialize_model."""
        if self._pool_mu_floors is None:
            return None
        if self._pool_mu_peak is None:
            return np.array(self._pool_mu_floors, dtype=float)
        return np.maximum(self._pool_mu_floors, EXHAUSTION_RELATIVE_K
                          * self._pool_mu_rtol * self._pool_mu_peak)

    def _pool_held(self, p):
        """True once pool ``p`` has held more than its absolute floor in some
        moment on an accepted state, the initial state included (the I-244
        peak; before the peak is laid out, the initial state itself). A pool
        born empty that never filled has not been drained: it is neither
        ended by end_exhausted_pools nor announced exhausted."""
        floors = self._pool_mu_floors
        peak = self._pool_mu_peak
        if peak is not None:
            return any(peak[p, k] > floors[p, k] for k in range(3))
        y0 = getattr(self, "y0", None)
        if y0 is None or floors is None:
            return True
        idx = self.polymer_pools[p].mu_indices
        return any(abs(float(y0[int(idx[k])])) > floors[p, k] for k in range(3))

    def _pool_born_empty(self, p):
        """True if pool ``p``'s initial moments are exactly zero (a spawned
        daughter, or a deck pool given none)."""
        y0 = getattr(self, "y0", None)
        if y0 is None:
            return False
        idx = self.polymer_pools[p].mu_indices
        return all(float(y0[int(idx[k])]) == 0.0 for k in range(3))

    def end_exhausted_pools(self):
        """End every pool that has emptied, on the ACCEPTED state ``self.y``
        (I-244; simulate's step() calls it once per accepted step, before the
        accepted-state check). A pool is exhausted when every moment lies
        inside its accepted floor, ``|mu_k| <= F_k``
        (``pool_moment_accepted_floors``), and it ever held more than its
        absolute floor (a pool born empty that never filled is not ended).
        It is then ended by one named event per episode:

        * ``pool_exhaustion_events`` gains ``{"event": "pool_exhausted",
          "pool", "time_s", "residual_moments_mol", "floors_mol",
          "residual_mass_g"}`` (the residue's condensed_mass_g) and a
          ``POOL EXHAUSTED`` warning is logged;
        * the pool's mu0/mu1/mu2 are set to exactly 0.0 and the residue is
          added to ``_pool_exhaustion_residue``, so the balances stay exact:
          ``code_convention_inventory() + exhausted_residue_atoms()`` and
          ``physical_elemental_inventory() + exhausted_residue_atoms()`` are
          what they were before the clamp (the residue chains' end groups go
          to the end-group ledger, like any chain that leaves the moment
          accounting), and ``exhausted_residue_mass_g()`` names the mass;
        * DASPK is restarted from the clamped state at the same time (its
          history holds the unclamped one), and the reactor continues.

        While the pool stays inside its floors it is not ended again;
        leaving them (refilled) starts a new episode. Returns the pools
        ended now."""
        if (self._pool_ended is None or self._pool_mu_floors is None
                or not self.polymer_pools):
            return []
        y = self.y
        floors = self.pool_moment_accepted_floors()
        ended = []
        clamp = np.array(y, dtype=float, copy=True)
        for p, pool in enumerate(self.polymer_pools):
            idx = [int(i) for i in pool.mu_indices]
            mu = [float(y[i]) for i in idx]
            inside = all(abs(mu[k]) <= floors[p, k] for k in range(3))
            if not inside:
                self._pool_ended[p] = 0
                continue
            if self._pool_ended[p] or not self._pool_held(p):
                continue
            self._pool_ended[p] = 1
            for k in range(3):
                self._pool_exhaustion_residue[p, k] += mu[k]
                clamp[idx[k]] = 0.0
            # An X-loss daughter's intact-site inventory leaves with its
            # chains (I-260): booked beside the residue moments.
            if getattr(self, "sgh_enabled", None) is not None and self.sgh_enabled[p]:
                for row in range(int(self.sgh_row_start[p]), int(self.sgh_row_end[p])):
                    if self.sgh_in_pool[row]:
                        z_slot = int(self.sgh_z_slot[row])
                        self._pool_exhaustion_z[row] = (
                            self._pool_exhaustion_z.get(row, 0.0) + float(y[z_slot]))
                        clamp[z_slot] = 0.0
            event = {"event": "pool_exhausted", "pool": pool.label,
                     "time_s": float(self.t), "residual_moments_mol": mu,
                     "floors_mol": [float(f) for f in floors[p]],
                     "residual_mass_g": float(pool.condensed_mass_g(mu[0],
                                                                    mu[1]))}
            self.pool_exhaustion_events.append(event)
            ended.append(pool.label)
            logging.warning(
                "POOL EXHAUSTED: pool %s at t = %.6g s: every moment inside "
                "its accepted floor (mu0=%.6e, mu1=%.6e, mu2=%.6e mol; "
                "floors %.6e/%.6e/%.6e mol = max(SMALL_EPS, %.0f*atol, "
                "%g*rtol*peak|mu_k|)). Moments set to 0; residue %.6e g "
                "booked (exhausted_residue_mass_g); the reactor continues.",
                pool.label, self.t, mu[0], mu[1], mu[2], floors[p, 0],
                floors[p, 1], floors[p, 2], EXHAUSTION_FLOOR_K,
                EXHAUSTION_RELATIVE_K, event["residual_mass_g"])
        if ended and np.any(clamp != np.asarray(y, dtype=float)):
            dydt = self.residual(self.t, clamp, np.zeros(self.neq, float),
                                 self.senpar)[0]
            self.initialize(self.t, clamp, dydt, self.senpar,
                            self.atol_array, self.rtol_array)
        return ended

    def _update_pool_mu_peak(self):
        """Raise each pool's peak |mu_k| to the ACCEPTED state ``self.y``
        (the scale of the relative accepted floor)."""
        if self._pool_mu_peak is None:
            return
        y = self.y
        for p, pool in enumerate(self.polymer_pools):
            for k in range(3):
                v = abs(float(y[int(pool.mu_indices[k])]))
                if v > self._pool_mu_peak[p, k]:
                    self._pool_mu_peak[p, k] = v

    def exhausted_residue_atoms(self, defect_atoms=None):
        """{element: mol} the pool_exhausted events took from the moment
        state, in the code convention (repeat units minus per-chain defect
        atoms; the residue chains' ends are on the end-group ledger):
        sum_pools R1*repeat_unit_atoms - R0*defect atoms, R the residue
        booked so far; ``defect_atoms`` as code_convention_inventory takes
        it. Empty when no pool was exhausted."""
        defect_atoms = defect_atoms or {}
        out = {}
        residue = self._pool_exhaustion_residue
        if residue is None:
            return out
        for p, pool in enumerate(self.polymer_pools):
            if not (residue[p, 0] or residue[p, 1]):
                continue
            if not pool.repeat_unit_atoms:
                raise ValueError(
                    f"Pool {pool.label}: repeat_unit_atoms is not declared; "
                    f"the exhausted residue cannot be reported in elements.")
            for el, n in pool.repeat_unit_atoms:
                out[el] = out.get(el, 0.0) + float(residue[p, 1]) * n
            if pool.chain_mass_defect_g_mol > 0.0:
                if pool.label not in defect_atoms:
                    raise ValueError(
                        f"Pool {pool.label}: chain_mass_defect_g_mol="
                        f"{pool.chain_mass_defect_g_mol:g} > 0 but no "
                        f"defect_atoms[{pool.label!r}] (elements lost per "
                        f"chain) was supplied.")
                for el, n in defect_atoms[pool.label].items():
                    out[el] = out.get(el, 0.0) - float(residue[p, 0]) * n
            for row, lost in self._exhausted_further(p):
                for el, n in self._core_species_atoms[int(self.sgh_gas[row])].items():
                    out[el] = out.get(el, 0.0) - lost * n
        return out

    def _exhausted_further(self, int p):
        """The further losses (I-260) of an exhausted X-loss daughter's
        residue, [(row, mol)]: s*R1 - Z_residue - R0 per in-pool row."""
        out = []
        residue = self._pool_exhaustion_residue
        for row, z in self._pool_exhaustion_z.items():
            if self.sgh_row_start[p] <= row < self.sgh_row_end[p]:
                out.append((row, self.sgh_sites[row] * float(residue[p, 1]) - z
                            - float(residue[p, 0])))
        return out

    def exhausted_residue_mass_g(self):
        """{pool label: g} the pool_exhausted events took from each pool
        (condensed_mass_g of the booked residue), for the pools ended so
        far."""
        out = {}
        residue = self._pool_exhaustion_residue
        if residue is None:
            return out
        for p, pool in enumerate(self.polymer_pools):
            if any(e["pool"] == pool.label for e in self.pool_exhaustion_events):
                out[pool.label] = float(pool.condensed_mass_g(residue[p, 0],
                                                              residue[p, 1]))
                for row, lost in self._exhausted_further(p):
                    out[pool.label] -= lost * float(self.sgh_mw_X[row])
        return out

    def _update_sgh_inventory_guard(self, y):
        """SGH kernel-v2 Z-inventory census on a TRIAL state (§1g, P1-A).
        Called from residual() on the RAW Newton trial state -- which is NOT
        an accepted solution -- so this method is TRIAL-STATE-SAFE: it NEVER
        raises and NEVER mutates ``y``. Z decays stiffly (dZ/dt = -kZ, FR1 k
        huge), so a recoverable Newton trial can momentarily undershoot the
        tiny -floor as Z -> 0; hard-raising here would turn that transient
        into a fatal integration error (the trial-state-hard-raise bug class).
        The RHS's local ``max(0.0, y[z_slot])`` read already sanitizes a
        negative trial Z, so this method only CENSUSES (warn-once per
        (pool, channel)) any Z sitting outside the physical band
        ``[0, sites*mu1]`` -- BOTH the beyond-floor case and the in-band case
        emit ``_emit_sgh_inventory_census``.

        The HARD physical-band invariant is enforced on the ACCEPTED state
        ONLY, by ``_assert_sgh_inventory_accepted`` (self.y), where an
        out-of-band Z is genuine corruption rather than a recoverable Newton
        trial excursion. This mirrors the TA-side resolution (clamp/sanitize
        the trial state, hard-check the accepted state); the residual is not
        an accepted-state context. Physical band ``0 <= Z <= sites*mu1``: mu1
        is pinned STATIONARY by the §1h build guard, so the upper bound equals
        the seeded Z(0). Cheap no-op when ``num_sgh_z == 0`` (legacy/non-SGH
        runs untouched)."""
        cdef int p, sgh_row, z_slot
        cdef double z_val, upper, floor_z, mu1_amount, a_z
        if self.num_sgh_z == 0:
            return
        sgh_en = getattr(self, "sgh_enabled", None)
        if sgh_en is None:
            return
        atol_arr = getattr(self, "atol_array", None)
        for p in range(len(self.polymer_pools)):
            if not sgh_en[p]:
                continue
            pool = self.polymer_pools[p]
            mu1_amount = max(0.0, y[int(pool.mu_indices[1])])
            for sgh_row in range(int(self.sgh_row_start[p]),
                                 int(self.sgh_row_end[p])):
                z_slot = int(self.sgh_z_slot[sgh_row])
                if z_slot < 0:
                    continue
                z_val = y[z_slot]
                upper = self.sgh_sites[sgh_row] * mu1_amount
                # atol-derived floor on the Z slot (mole basis; same basis
                # atol governs -- mirrors the moment-floor derivation). Used
                # only to size the census-reported band; never a raise gate
                # here (trial-state-safe).
                if atol_arr is None:
                    a_z = 0.0
                elif np.ndim(atol_arr) == 0:
                    a_z = float(atol_arr)
                elif z_slot < len(atol_arr):
                    a_z = float(atol_arr[z_slot])
                else:
                    a_z = float(np.max(atol_arr))
                floor_z = max(SMALL_EPS, EXHAUSTION_FLOOR_K * a_z)
                if z_val < 0.0 or z_val > upper:
                    self._emit_sgh_inventory_census(
                        p, sgh_row, z_val, upper, floor_z)

    def _assert_sgh_inventory_accepted(self):
        """SGH kernel-v2 ACCEPTED-state Z-inventory invariant (§1g, P1-A): the
        HARD physical-band check the residual guard is deliberately too
        trial-unsafe to make. Reads ``self.y`` -- the ACCEPTED integrator
        solution (see the _phase_gate_flux_census call site: base.pyx invokes
        that override once per accepted snapshot, with self.y holding the
        accepted state) -- so a Z that has left ``[0, sites*max(0, mu1)]``
        beyond the atol-derived floor is genuine integrator corruption, not a
        recoverable Newton trial excursion. HARD-RAISES with the SAME
        physical-band condition the residual guard used to enforce:
        ``z < -floor_z or z > sites*max(0, mu1) + floor_z``,
        ``floor_z = max(SMALL_EPS, EXHAUSTION_FLOOR_K*atol[z_slot])``. Mirrors
        the moment driver's accepted-state posture (hard error, never
        max(..., SMALL_EPS)). Cheap no-op when ``num_sgh_z == 0``."""
        cdef int p, sgh_row, z_slot
        cdef double z_val, upper, floor_z, mu1_amount, a_z
        if self.num_sgh_z == 0:
            return
        sgh_en = getattr(self, "sgh_enabled", None)
        if sgh_en is None:
            return
        y = self.y
        atol_arr = getattr(self, "atol_array", None)
        for p in range(len(self.polymer_pools)):
            if not sgh_en[p]:
                continue
            pool = self.polymer_pools[p]
            mu1_amount = max(0.0, y[int(pool.mu_indices[1])])
            for sgh_row in range(int(self.sgh_row_start[p]),
                                 int(self.sgh_row_end[p])):
                z_slot = int(self.sgh_z_slot[sgh_row])
                if z_slot < 0:
                    continue
                z_val = y[z_slot]
                upper = self.sgh_sites[sgh_row] * mu1_amount
                if atol_arr is None:
                    a_z = 0.0
                elif np.ndim(atol_arr) == 0:
                    a_z = float(atol_arr)
                elif z_slot < len(atol_arr):
                    a_z = float(atol_arr[z_slot])
                else:
                    a_z = float(np.max(atol_arr))
                floor_z = max(SMALL_EPS, EXHAUSTION_FLOOR_K * a_z)
                if z_val < -floor_z or z_val > upper + floor_z:
                    raise ValueError(
                        f"Polymer pool '{pool.label}' side_group_homolysis "
                        f"channel '{self.sgh_channel_labels[sgh_row]}' Z "
                        f"inventory left the physical band beyond the floor on "
                        f"the ACCEPTED state: Z={z_val!r} mol vs "
                        f"[0, sites*mu1={upper:.6e}] +/- {floor_z:.6e} mol (= "
                        f"max(SMALL_EPS, {EXHAUSTION_FLOOR_K:.0f}*"
                        f"atol[z_slot])). A Z outside [0, sites*mu1] beyond the "
                        f"integrator's own error budget on an accepted state "
                        f"is integrator corruption, not depletion -- refusing "
                        f"to continue (SGH-v2 Z-bound rule: hard error, never "
                        f"max(..., SMALL_EPS)).")

    def _emit_sgh_inventory_census(self, int p, int sgh_row, double z_val,
                                   double upper, double floor_z):
        """Once-per-(pool, channel)-per-rebuild in-band Z census (mirrors
        _emit_pool_exhaustion_census): log-only, nothing suppressed. Cold
        path -- string work happens at most once per key per rebuild."""
        pool = self.polymer_pools[p]
        ch = self.sgh_channel_labels[sgh_row]
        key = (pool.label, ch)
        emitted = getattr(self, "_sgh_inventory_census_emitted", None)
        if emitted is None:
            emitted = set()
            self._sgh_inventory_census_emitted = emitted
        if key in emitted:
            return
        emitted.add(key)
        logging.warning(
            "SGH INVENTORY CENSUS: pool %s channel %s Z=%.6e mol sits "
            "in-band outside [0, sites*mu1=%.6e] within +/- floor %.6e mol "
            "(census only; the RHS max(0.0, .) read clamps it, y unchanged).",
            pool.label, ch, z_val, upper, floor_z)

    def _side_loss_further(self, int p, y):
        """An X-loss daughter's further losses beyond the first per chain
        (I-260), per in-pool row: ([(row, mol)], g), mol = s*mu1 - Z - mu0
        over the moment and inventory slots of state ``y``."""
        pool = self.polymer_pools[p]
        mu0 = float(y[pool.mu_indices[0]])
        mu1 = float(y[pool.mu_indices[1]])
        rows, grams = [], 0.0
        for row in range(int(self.sgh_row_start[p]), int(self.sgh_row_end[p])):
            if not self.sgh_in_pool[row]:
                continue
            n = (self.sgh_sites[row] * mu1 - float(y[int(self.sgh_z_slot[row])])
                 - mu0)
            rows.append((row, n))
            grams += n * self.sgh_mw_X[row]
        return rows, grams

    def get_total_polymer_condensed_mass_g(self, y=None):
        """r86 terminationPolymerConversion metric (defect-adjusted):

            M_pool = max(0, mu1_p*monomer_mw_g_mol_p
                            - mu0_p*chain_mass_defect_g_mol_p)
            M_poly = sum over ALL solver polymer pools of M_pool

        summed over EVERY pool the solver carries at this rebuild --
        configured, spawned, homolysis-daughter and deprop-daughter pools
        apply the schema-2.7 mass-defect formula through
        PolymerPoolConfig.condensed_mass_g. SGH kernel-v2 carriers instead
        book their Br loss through the explicit Z inventory (there is no
        side-group feature pool in v2): for an sgh_enabled carrier the mass
        is mu1*MW - Σ_row lost_X, lost_X = max(0, sites*tail_mu1 - Z)*M_X, so
        an FR1 Br-loss carrier never falsely retains the Br that has left as
        gas X. Explicit-DP
        oligomer carriers hold genuine polymer repeat-unit mass OUTSIDE
        the mu state slots (the init consistency check subtracts their
        share from mu1, r86 probe 3), so their contribution is added back
        via _explicit_moment_contributions. Only each pool's FINAL mass
        contribution is clamped at zero (numerical fuzz); individual
        moments are never clamped. Returns grams; reads self.y when no
        state vector is passed."""
        cdef int p, i0, i1
        cdef double total, mu0_mol, mu1_mol, V_poly
        if y is None:
            y = self.y
        if self.polymer_pools is None or len(self.polymer_pools) == 0:
            return None
        V_poly = self.V_poly
        total = 0.0
        for p in range(len(self.polymer_pools)):
            pool = self.polymer_pools[p]
            i0 = self.pool_mu0_indices[p] if self.pool_mu0_indices is not None else -1
            i1 = self.pool_mu1_indices[p] if self.pool_mu1_indices is not None else -1
            if i0 < 0 or i1 < 0:
                # Pre-initialize_model fallback: the config-time layout.
                i0, i1 = pool.mu_indices[0], pool.mu_indices[1]
            mu0_mol = y[i0]
            mu1_mol = y[i1]
            if pool.explicit_dp_to_species_index:
                e0, e1, _e2 = _explicit_moment_contributions(
                    V_poly, y, pool.explicit_dp_to_species_index)
                mu0_mol += e0 * V_poly
                mu1_mol += e1 * V_poly
            sgh_en = getattr(self, "sgh_enabled", None)
            if (sgh_en is not None and p < len(sgh_en) and sgh_en[p]
                    and self.sgh_in_pool[self.sgh_row_start[p]]):
                # An X-loss daughter (I-260): the first group per chain is
                # its defect, the further ones s*mu1 - Z - mu0 (a function
                # of state, never clamped: every kernel's terms then close).
                total += max(0.0, pool.condensed_mass_g(mu0_mol, mu1_mol)
                             - self._side_loss_further(p, y)[1])
            elif (sgh_en is not None and p < len(sgh_en) and sgh_en[p]
                    and self.sgh_z_slot[self.sgh_row_start[p]] >= 0):
                # SGH kernel-v2 (§1e): Br leaves the melt through the
                # per-channel Z depletion state, NOT a moment debit or a
                # feature pool. Br lost so far on each row is the seeded
                # inventory sites*mu1(0) minus the remaining Z(t). The §1h
                # guard pins the carrier mu1 STATIONARY, so sites*tail_mu1(t)
                # == Z(0) at all times and lost = max(0, sites*tail_mu1 - Z)
                # tracks depleted Br EXACTLY: 0 at t0 (Z == Z(0)) and
                # sites*mu1(0)*M_X at full depletion (Z -> 0 == ALL the Br).
                # A v2 carrier's chain_mass_defect_g_mol is 0, so the intact
                # mass reduces to mu1*MW; subtract the Br that has left.
                # (tail_mu1 matches the Z-seed basis -- the moment slot, not
                # the explicit-DP-augmented mu1_mol.)
                tail_mu1 = max(0.0, y[i1])
                lost_x = 0.0
                for sgh_row in range(int(self.sgh_row_start[p]),
                                     int(self.sgh_row_end[p])):
                    z_val = y[int(self.sgh_z_slot[sgh_row])]
                    lost_x += max(0.0, self.sgh_sites[sgh_row] * tail_mu1
                                  - z_val) * self.sgh_mw_X[sgh_row]
                total += max(0.0, mu1_mol * pool.monomer_mw_g_mol - lost_x)
            else:
                total += max(0.0, pool.condensed_mass_g(mu0_mol, mu1_mol))
        return total

    def _chain_bundle(self, int pool_idx, y, double V_poly, bint end_group):
        """
        Expected (chains, units, units^2) carried by ONE picked chain of pool
        ``pool_idx``: (1, E[k], E[k^2]). end_group=True -> uniform chain pick
        (rate ~ mu0): (1, mu1/mu0, mu2/mu0). Otherwise length-biased pick
        (rate ~ mu1): (1, mu2/mu1, mu3/mu1) with the guarded mu3 closure.
        Returns (b0, b1, b2, mu2_ok); b0 == 0.0 means the pool is too empty to
        move a chain (denominator below SMALL_EPS) -- caller skips the term.
        mu2_ok False means apply b0/b1 but skip the mu2 component (mu3 = inf).
        """
        idx0, idx1, idx2 = self.polymer_pools[pool_idx].mu_indices
        mu0 = _softclamp(y[idx0], self._softclamp_lam) / V_poly
        mu1 = _softclamp(y[idx1], self._softclamp_lam) / V_poly
        mu2 = _softclamp(y[idx2], self._softclamp_lam) / V_poly
        if end_group:
            if mu0 <= SMALL_EPS:
                return 0.0, 0.0, 0.0, False
            return 1.0, mu1 / mu0, mu2 / mu0, True
        if mu1 <= SMALL_EPS:
            return 0.0, 0.0, 0.0, False
        mu3 = _safe_mu3_from_mu012(mu0, mu1, mu2)
        if np.isfinite(mu3):
            return 1.0, mu2 / mu1, mu3 / mu1, True
        return 1.0, mu2 / mu1, 0.0, False

    def _bundle_availability_cap(self, int pool_idx, y, double V_poly,
                                 bint end_group):
        """P1-1 directional bundle throttle cap (adjudicated regen-#2 fix):
        the largest event-site density [mol/m^3_poly] the DEBITED pool can
        sustain so that EVERY moment's drain vanishes linearly near
        exhaustion:

            S_cap = min over POSITIVE bundle terms b_k of mu_k/(V_poly*b_k)

        One cross-pool VE/MIGRATION event debits the FULL per-chain bundle
        (b0, b1, b2) from the debited pool, so a leg whose event rate scales
        with mu1 drains mu1 at rate ~ k*C*mu1*(mu2/mu1) = k*C*mu2 -- the
        per-event mu2/mu1 debit CANCELS the linear mu1-site self-limit the
        pre-P1-1 code relied on (that round-22 P2b TODO's premise is FALSE),
        and once the out-of-cone mu3 closure freezes (returns 0) the mu2
        debit stops while the mu1 drain stays constant-rate: the regen-#2
        mu1-negative crash. Scaling the direction's event rate by
        min(S_base, S_cap) bounds every moment's drain by ~k*C*mu_k.
        For end-group bundles at a consistent state all three caps equal
        mu0/V_poly == S_base (no-op). Returns 0.0 for an empty bundle: no
        chains to move, so the direction's GAS flux is throttled off WITH
        its moment flux (they must never diverge).

        Round-27 P1-A re-adjudication: for healthy polydisperse pools the
        mu1-scaled cap mu1/b1 = mu1^2/mu2 < mu1 binds IN BULK whenever
        PDI > 1, rewriting healthy cross-pool kinetics -- so this cap is
        applied TAIL-ONLY via _bundle_limited_site (near-exhaustion bundle
        limiter, tail-only smoothstep; C1 soft-min, cone-margin drain
        guard), never as a global hard-min. Callers must go through
        _bundle_limited_site.

        Round-29 N2 (regen-#3 IDID=-7 adjudication; see the
        BUNDLE_LIMITER_SOFTMIN_P block comment): the min over cap terms
        is the reciprocal p-norm soft-min
        softmin_p(x_i) = (sum x_i^(-p))^(-1/p) (C1 in-band; evaluated in
        the overflow-safe factored form m*(sum (m/x_i)^p)^(-1/p), m the
        hard min, algebraically identical). The round-30 N1 cone-margin
        drain gate is NOT part of this cap: it is a structurally
        independent second gate inside _bundle_limited_site (folding it
        in here would dilute it through the E blend).
        (mirrored in the numpy oracle consumer -- keep in sync)"""
        cdef double m, acc, p, t0, t1, t2
        b0, b1, b2, mu2_ok = self._chain_bundle(pool_idx, y, V_poly,
                                                end_group)
        if b0 <= 0.0:
            return 0.0
        mu_idx = self.polymer_pools[pool_idx].mu_indices
        t0 = _softclamp(y[mu_idx[0]], self._softclamp_lam) / (V_poly * b0)
        t1 = -1.0
        t2 = -1.0
        if b1 > 0.0:
            t1 = _softclamp(y[mu_idx[1]], self._softclamp_lam) / (V_poly * b1)
        if b2 > 0.0:
            t2 = _softclamp(y[mu_idx[2]], self._softclamp_lam) / (V_poly * b2)
        # softmin_p over the included (non-negative) terms; any zero term
        # zeroes the whole cap (softmin_p <= min).
        m = t0
        if t1 >= 0.0 and t1 < m:
            m = t1
        if t2 >= 0.0 and t2 < m:
            m = t2
        if m <= 0.0:
            return 0.0
        p = BUNDLE_LIMITER_SOFTMIN_P
        acc = (m / t0) ** p
        if t1 >= 0.0:
            acc += (m / t1) ** p
        if t2 >= 0.0:
            acc += (m / t2) ** p
        return m * acc ** (-1.0 / p)

    def _pool_floor_distance(self, int pool_idx, y):
        """E: the debited pool's accepted-state floor distance,

            E = softmin_p(mu0/f0, mu1/f1, mu2/f2)

        with f_k = self._pool_mu_floors[pool_idx, k]
                 = max(SMALL_EPS, EXHAUSTION_FLOOR_K*atol[state]) -- the SAME
        per-state r81 accepted-state floors the negative-moment tripwire
        uses, so E is a dimensionless "how many r81 floors above exhaustion
        is the thinnest moment" measure (deck-independent, scale-aware; see
        the BUNDLE_LIMITER_E_LO/_E_HI block comment). Round-29 N2: the
        moment-min is the reciprocal p-norm SOFT-min (softmin_p <= min;
        within a factor 3^(-1/p) ~ 0.87 of it when all three terms tie), so
        E is C1 in y and the smoothstep weight w(E) carries no argmin-switch
        kink into the band. Moments are read from raw y with the same
        max(0, .) clamp as every site factor. Falls back to per-moment
        floors of max(SMALL_EPS, EXHAUSTION_FLOOR_K*atol) semantics via
        _pool_mu_floors, which initialize_model always lays out for
        configured pools."""
        cdef double m, acc, p, e_k
        cdef int k
        mu_idx = self.polymer_pools[pool_idx].mu_indices
        floors = self._pool_mu_floors
        cdef double e0 = _softclamp(y[mu_idx[0]], self._softclamp_lam) / floors[pool_idx, 0]
        cdef double e1 = _softclamp(y[mu_idx[1]], self._softclamp_lam) / floors[pool_idx, 1]
        cdef double e2 = _softclamp(y[mu_idx[2]], self._softclamp_lam) / floors[pool_idx, 2]
        m = e0
        if e1 < m:
            m = e1
        if e2 < m:
            m = e2
        if m <= 0.0:
            return 0.0
        p = BUNDLE_LIMITER_SOFTMIN_P
        acc = (m / e0) ** p + (m / e1) ** p + (m / e2) ** p
        return m * acc ** (-1.0 / p)

    def _bundle_limited_site(self, int pool_idx, y, double V_poly,
                             bint end_group, double s_base):
        """Near-exhaustion bundle limiter (tail-only smoothstep; C1
        soft-min, cone-margin drain guard) -- round-27 P1-A directional
        bundle throttle, round-29 N2 soft-min, round-30 N1 independent
        cone gate. TWO structurally independent gates in series:

        Stage 1, exhaustion tail limiter (E band):
            S_cap  = softmin_p over POSITIVE bundle terms mu_k/(V_poly*b_k)
            E      = softmin_p(mu0/f0, mu1/f1, mu2/f2)  (accepted-state
                                                  floor distance, unitless)
            w      = smoothstep(E; E_lo, E_hi) = 3e^2 - 2e^3,
                     e = clamp((E - E_lo)/(E_hi - E_lo), 0, 1)
            S_free = w*S_base + (1-w)*softmin_p(S_base, S_cap)
                     (== s_base EXACTLY when E >= E_hi -- early return,
                     bit-for-bit the pre-hard-min bulk law)

        Stage 2, cone-margin drain gate (M band, round-30: INDEPENDENT of
        E -- an out-of-cone pool with bulk-scale moments must not slip
        through the E early return, and an in-band out-of-cone pool must
        not keep draining at w*S_base). Applies ONLY to non-end-group
        (length-biased, b1 = mu2/mu1) cone-shrinking debits (b1 > b0;
        b0 = 1):
            Q10    = mu1 - mu0            (realizability cone margin)
            Q10 <= 0            -> S_eff = 0 REGARDLESS of E
            M      = Q10 / max(f0, f1)    (margin distance, floor units)
            M >= M_hi           -> S_eff = S_free EXACTLY (early return)
            S_cone = Q10 / (V_poly*(b1 - b0))
            M <= M_inner        -> narrowed completion: 0 EXACTLY for a
                                 resolved b1-1, u*softmin_p(S_free, S_cone)
                                 inside I-090's unresolved b1 neighbourhood
            M_inner < M < M_lo  -> C1 blend from the narrowed completion
                                 to softmin_p(S_free, S_cone)
            M_lo <= M < M_hi    -> original v-smoothstep blend of that
                                 soft-min cap into S_free
        Non-cone-shrinking debits (b1 <= b0) pass S_free through
        untouched, as do empty pools (stage 1 already throttled those).
        EXCEPTION, end-group rows (round-62 N5b adjudicated fix): the
        gate above is skipped entirely for end-group rows, which return
        S_free unconditionally. End-group uniform-pick drains debit
        (dmu0, dmu1) = (1, mu1/mu0)*rate, so Q10 = mu1 - mu0 decays
        multiplicatively and never crosses zero -- the "out of cone"
        state this gate exists to catch cannot occur on this surface.
        The stage-2 cap also simplifies EXACTLY here (b1c - 1 =
        q10/y0c, so S_cone == y0c/V_poly identically, which is always
        >= S_base >= S_free), i.e. the cap is structurally non-binding
        for end-group rows, so evaluating the q10/b1c cancellation
        surface at all would only inject noise (up to ~1.08% measured
        near the all-monomer boundary) for a gate that can never fire.
        A prior revision instead forced a hard zero at end-group Q10 <=
        0 by analogy with the non-end-group path; that analogy does not
        hold here (see the in-code comment at the end_group branch) and
        broke a correct HARD PIN regression test -- see
        test_rtol_1e6_near_floor_trajectory_faithful.

        softmin_p(x_i) = (sum x_i^(-p))^(-1/p), the reciprocal p-norm
        soft-min (p = BUNDLE_LIMITER_SOFTMIN_P; see its block comment and
        the CONE_MARGIN_M_LO/_M_HI comment for the adjudication and the
        band design).

        Healthy IN-CONE bulk pools take both early returns: s_base comes
        back bitwise, so byte-pins and bulk reversible-row detailed
        balance (C_G* = Keq*S_base(A)/S_base(B)) are untouched. N5b
        (round-62 DASSL-hang root cause) keeps resolved b1 shapes at exact
        zero through M_inner=50 accepted-state floors. I-072's M=50..100
        interval is numerical regularization, not a claim that those states
        are sub-floor noise: widening the cubic handoff reduces its maximum
        blend slope and curvature by 50x and 2500x versus W=1. It reaches the
        original band cap at M_lo with zero endpoint slopes. The law for
        M > M_lo is bit-for-bit unchanged; equality necessarily takes the
        continuous edge value. See the constant block for measured replay,
        deck, reversible-pair, and W=100 counterfactual effects.
        s_base is
        the direction's adjudicated site law (mu1/V_poly or mu0/V_poly
        per row scaling, including the pre-existing a>0/a<0 VE
        min(mu0, mu1/|a|) throttle), computed by the caller from the pool
        this direction DEBITS (forward debits src, reverse debits dst).
        (mirrored in get_reaction_rates and in the numpy oracle consumer --
        keep in sync)"""
        cdef double e_dist, e_n, w, cap, narrowed, m, acc, p, s_free
        cdef double y0c, y1c, y2c, b1c, q10, m_dist, v_n, v, s_cone, f_c
        cdef double z_n, z
        # ---- stage 1: near-exhaustion tail limiter --------------------
        e_dist = self._pool_floor_distance(pool_idx, y)
        if e_dist >= BUNDLE_LIMITER_E_HI:
            s_free = s_base
        else:
            cap = self._bundle_availability_cap(pool_idx, y, V_poly,
                                                end_group)
            # softmin_p(s_base, cap): the reciprocal p-norm is
            # associative, so folding s_base into the cap's own softmin
            # IS the flat softmin over all terms. A zero cap (empty
            # bundle) zeroes the whole soft cap.
            if cap <= 0.0 or s_base <= 0.0:
                cap = 0.0
            else:
                m = cap if cap < s_base else s_base
                p = BUNDLE_LIMITER_SOFTMIN_P
                acc = (m / s_base) ** p + (m / cap) ** p
                cap = m * acc ** (-1.0 / p)
            if e_dist <= BUNDLE_LIMITER_E_LO:
                s_free = cap
            else:
                e_n = ((e_dist - BUNDLE_LIMITER_E_LO)
                       / (BUNDLE_LIMITER_E_HI - BUNDLE_LIMITER_E_LO))
                w = e_n * e_n * (3.0 - 2.0 * e_n)
                s_free = w * s_base + (1.0 - w) * cap
        # ---- stage 2: cone-margin drain gate (independent of E) -------
        mu_idx = self.polymer_pools[pool_idx].mu_indices
        if end_group:
            # End-group uniform-pick drains debit (dmu0, dmu1) =
            # (1, mu1/mu0)*rate, so d(Q10)/dt = -(rate/mu0)*Q10 -- Q10
            # decays multiplicatively and never crosses zero (every
            # length-1 chain still has its own end group; "nothing left
            # to drain at b1c==1" is a non-end-group/internal-unit
            # notion, it does not apply here). The stage-2 cone cap
            # simplifies EXACTLY on this surface: b1c - 1 = q10/y0c, so
            # s_cone == y0c/V_poly identically, which is >= s_base >=
            # s_free always -- the cap never binds for end-group rows.
            # The gate is therefore a structural no-op here: evaluating
            # the q10/b1c cancellation surface at all only adds error
            # (up to ~1.08% measured near the all-monomer boundary) for
            # a cap that can never bind. Return s_free directly --
            # byte-exact s_base in the bulk, and correct (not a hard
            # zero) right at the b1c==1 boundary. (A prior revision of
            # this branch forced a hard zero at q10<=0 by analogy with
            # the non-end-group path below; that analogy does not hold
            # for end-group rows and broke a correct HARD PIN regression
            # test -- see test_rtol_1e6_near_floor_trajectory_faithful.)
            return s_free
        y0c = _softclamp(y[mu_idx[0]], self._softclamp_lam)
        y1c = _softclamp(y[mu_idx[1]], self._softclamp_lam)
        y2c = _softclamp(y[mu_idx[2]], self._softclamp_lam)
        if y1c / V_poly <= SMALL_EPS:
            return s_free
        b1c = y2c / y1c         # length-biased pick E[k]
        if b1c <= 1.0:          # b0 = 1: debit does not shrink Q10
            return s_free
        q10 = y1c - y0c
        if q10 <= 0.0:
            return 0.0          # out of cone: drain OFF, whatever E
        floors = self._pool_mu_floors
        f_c = max(floors[pool_idx, 0], floors[pool_idx, 1])
        m_dist = q10 / f_c
        if m_dist >= CONE_MARGIN_M_HI:
            return s_free           # margin safely bulk: gate inactive
        # N5b cone-gate dead-band fix: below M_LO, q10 itself is
        # sub-floor-scale catastrophic-cancellation noise (the same
        # noise regime that motivated the b1c/q10 consistency fix
        # above). At this scale the row's true rate is bounded by
        # roughly f_c/(V_poly*(b1c-1)) ~ 1e2*atol/(V_poly*(b1c-1)),
        # i.e. of order the accepted-state floor itself -- about 15
        # orders of magnitude below any deck-observable rate, squarely
        # inside the model's own error budget (the same adjudicated
        # "below-floor is noise, not signal" logic as the r81 floors).
        # Returning cap here trusted that noise's magnitude AND sign, so the
        # deep dead band returns the exact hard zero instead.
        #
        # I-090 NARROWING. That argument has a domain: it holds wherever
        # S_cone = q10/(V_poly*(b1c - 1)) is a bound at all. It is VOID on
        # the b1c == 1 surface, where S_cone DIVERGES and therefore caps
        # nothing -- there the hard zero is not "declining to trust noise",
        # it is discarding a rate this gate has no business capping, and
        # doing it across a jump (s_free on the b1c <= 1 side of the
        # surface, 0.0 on the other) that the corrector cannot step over.
        # That jump is the measured never-returns signature: the poly_102
        # mod_2..mod_5 pools ride |b1c - 1| <= 5e-5 through the crash
        # window and DASPK's step size collapses against it indefinitely.
        #
        # So the dead band now runs the divergence-repairing completion
        # softmin_p(s_free, S_cone) on a noise-scale neighbourhood of
        # b1c == 1 ONLY, and hands back to the SAME exact hard zero outside
        # it. Both handoffs are C1, so the narrowing gives back nothing the
        # completion bought and introduces no new edge:
        #   * b1c -> 1+ : S_cone -> +inf so the completion -> s_free, and
        #     u -> 1 with u' = 0, so the law meets the b1c <= 1 branch's
        #     s_free continuously and with matching slope.
        #   * b1c - 1 -> b1_band : u -> 0 with u' = 0, so the law meets the
        #     N5b hard zero continuously, and BIT-FOR-BIT 0.0 beyond it --
        #     not merely small. At O(1) b1c - 1 the round-62 hard zero remains
        #     bit-for-bit unchanged through M_INNER; M_INNER < M < M_LO ramps
        #     to the original band cap. That smoothstep has
        # zero slope at both ends, removing the pre-existing M_LO jump without
        # changing any value above M_LO. The value exactly at M_LO necessarily
        # becomes the continuous edge cap.
        # The neighbourhood half-width is derived, not fitted.  I-090 used
        # only the propagated error weight of the two moments b1c is built
        # from.  Its absolute floor terms divided by y1c, however, diverge as
        # an entire bundle becomes sub-floor: a fixed, resolved distribution
        # shape is then misclassified as arbitrarily uncertain and its drain
        # reopens.  I-067 conditions that band with the scale-free fraction
        # of the complete three-moment inventory carried by the cone margin:
        #
        #   B_err = (rtol*(y1+y2) + f1+f2) / y1
        #   B_Q   = softmin_p(B_err, Q10/(y0+y1+y2))
        #   h     = smoothstep(E/E_lo)
        #   B     = h*B_err + (1-h)*B_Q              (E < E_lo)
        #
        # Q10/(y0+y1+y2) is dimensionless, scale invariant at fixed shape,
        # and contains no fitted constant.  Thus a resolved shape remains
        # resolved below the inventory floor, and the band follows the
        # remaining cone reserve when the absolute error term loses meaning.
        # h and h' meet (1, 0) at E_lo, so the E-axis blend is continuous,
        # joins the old law there with matching slope, and is exactly the old
        # law throughout the bulk.  The hard floor applied to B below can
        # still introduce a derivative kink, so the complete floored law is
        # not generally C1.  The b1 completion itself has zero endpoint
        # slopes at both of its existing handoffs.  See the I-090 width
        # derivation near CONE_B1_NOISE_REL_FLOOR.  This law is mirrored in
        # the numpy consumer and _s_eff; keep all three in sync.
        b1_band = ((self._cone_b1_rtol * (y1c + y2c)
                    + floors[pool_idx, 1] + floors[pool_idx, 2]) / y1c)
        if e_dist < BUNDLE_LIMITER_E_LO:
            e_n = e_dist / BUNDLE_LIMITER_E_LO
            w = e_n * e_n * (3.0 - 2.0 * e_n)
            cap = q10 / (y0c + y1c + y2c)
            m = cap if cap < b1_band else b1_band
            p = BUNDLE_LIMITER_SOFTMIN_P
            acc = (m / b1_band) ** p + (m / cap) ** p
            cap = m * acc ** (-1.0 / p)
            b1_band = w * b1_band + (1.0 - w) * cap
        if b1_band < CONE_B1_NOISE_REL_FLOOR:
            b1_band = CONE_B1_NOISE_REL_FLOOR
        b1_n = (b1c - 1.0) / b1_band
        s_cone = q10 / (V_poly * (b1c - 1.0))
        if s_free <= 0.0:
            return s_free
        m = s_cone if s_cone < s_free else s_free
        p = BUNDLE_LIMITER_SOFTMIN_P
        acc = (m / s_free) ** p + (m / s_cone) ** p
        cap = m * acc ** (-1.0 / p)
        if b1_n >= 1.0:
            narrowed = 0.0
        else:
            u = 1.0 - b1_n * b1_n * (3.0 - 2.0 * b1_n)
            narrowed = u * cap
        if m_dist <= CONE_MARGIN_M_INNER:
            return narrowed
        if m_dist < CONE_MARGIN_M_LO:
            z_n = ((m_dist - CONE_MARGIN_M_INNER)
                   / (CONE_MARGIN_M_LO - CONE_MARGIN_M_INNER))
            z = z_n * z_n * (3.0 - 2.0 * z_n)
            return z * cap + (1.0 - z) * narrowed
        v_n = ((m_dist - CONE_MARGIN_M_LO)
               / (CONE_MARGIN_M_HI - CONE_MARGIN_M_LO))
        v = v_n * v_n * (3.0 - 2.0 * v_n)
        return v * s_free + (1.0 - v) * cap

    # ------------------------------------------------------------------
    # 18.2e P5 (option (d), spar round-44): scoped-RHS user-side FD
    # Jacobian. PROTOTYPE pending spar review; feature-gated, default OFF
    # (constructor kwarg polymer_scoped_jacobian=False). The Jacobian
    # path is purely ADDITIVE: no residual law, softclamp, r81 floor, or
    # §4 dispatch semantics change anywhere.
    # ------------------------------------------------------------------
    def advance(self, tout):
        """Advance and leave scoped-Jacobian diagnostics at the accepted state.

        DASPK normally evaluates a full residual after a user Jacobian, but
        may return at ``tout`` immediately after the Jacobian sweep. Refresh
        only in that exceptional ordering; the discarded residual restores
        rate/census diagnostics without changing the accepted state.
        """
        result = super().advance(tout)
        if self._scoped_jac_diag_stale:
            self.residual(self.t, np.array(self.y, copy=True),
                          np.array(self.dydt, copy=True))
        self.record_realizability_census()
        return result

    def step(self, tout):
        """One integrator step (simulate's loop), then on the accepted state:
        the pools' peak moments, the end of any pool that has emptied
        (I-244, end_exhausted_pools), and the realizability census."""
        result = super().step(tout)
        self._update_pool_mu_peak()
        self.end_exhausted_pools()
        self.record_realizability_census()
        return result

    def record_realizability_census(self):
        """Census of the accepted state ``self.y`` (once per accepted time):
        per pool, whether its moments lie outside the realizable set of a
        chain-length distribution -- the cone mu1 >= mu0 >= 0 or
        mu0*mu2 >= mu1^2, each within the moments' exhaustion floors (the
        tolerances of the warn-once cone and variance censuses).  Pure
        moments can leave it (a kernel's daughter flux need not be a
        distribution's); that is a result to report, so nothing is changed:
        ``realizability_census[label]`` counts the accepted steps and those
        outside, the first and last time outside, the worst relative
        variance deficit (mu0*mu2 - mu1^2) / mu1^2 and cone gap
        (mu1 - mu0) / mu0, and the excursions [t_enter, t_leave or None];
        entering and leaving are logged."""
        census = getattr(self, "realizability_census", None)
        if not census or self._pool_mu_floors is None \
                or self._realizability_census_t == self.t:
            return
        self._realizability_census_t = self.t
        y = self.y
        for p, pool in enumerate(self.polymer_pools):
            c = census[pool.label]
            i0, i1, i2 = pool.mu_indices
            m0, m1, m2 = float(y[i0]), float(y[i1]), float(y[i2])
            f0, f1, f2 = (float(self._pool_mu_floors[p, k]) for k in range(3))
            cone = m1 < m0 - (f0 + f1)
            deficit = m0 * m2 - m1 * m1
            variance = deficit < -(f0 * abs(m2) + f2 * abs(m0)
                                   + 2.0 * f1 * abs(m1) + f0 * f2 + f1 * f1)
            c["accepted_steps"] += 1
            outside = cone or variance
            if outside:
                c["steps_outside"] += 1
                if c["first_outside_t"] is None:
                    c["first_outside_t"] = self.t
                c["last_outside_t"] = self.t
                if variance and m1 > 0.0:
                    c["worst_variance_deficit"] = min(
                        c["worst_variance_deficit"], deficit / (m1 * m1))
                if cone and m0 > 0.0:
                    c["worst_cone_gap"] = min(c["worst_cone_gap"],
                                              (m1 - m0) / m0)
            excursions = c["excursions"]
            if outside and (not excursions or excursions[-1][1] is not None):
                excursions.append([self.t, None])
                logging.warning(
                    "POOL REALIZABILITY CENSUS: pool %s left the realizable "
                    "set at t = %.6g s (mu0=%.6e, mu1=%.6e, mu2=%.6e mol; "
                    "%s). Recorded, not corrected.", pool.label, self.t, m0,
                    m1, m2, "mu1 < mu0" if cone else "mu0*mu2 < mu1^2")
            elif not outside and excursions and excursions[-1][1] is None:
                excursions[-1][1] = self.t
                logging.warning(
                    "POOL REALIZABILITY CENSUS: pool %s re-entered the "
                    "realizable set at t = %.6g s.", pool.label, self.t)

    def realizability_summary(self):
        """The census in one JSON-ready summary: per pool the accepted steps
        and those outside the realizable set, and the pools ever outside."""
        census = getattr(self, "realizability_census", None) or {}
        steps = {label: c["accepted_steps"] for label, c in census.items()}
        return {"accepted_steps": max(steps.values(), default=0),
                "steps_outside": {label: c["steps_outside"]
                                  for label, c in census.items()},
                "pools_outside": sorted(label for label, c in census.items()
                                        if c["steps_outside"]),
                "first_outside_t": {label: c["first_outside_t"]
                                    for label, c in census.items()
                                    if c["steps_outside"]}}

    def request_scoped_jacobian(self, enable=True):
        """Toggle the scoped-Jacobian request on an already-constructed
        system and re-run the arming validation immediately (the caller
        must still (re)run initialize_solver()/initialize() afterwards --
        pydas detects the capability inside DASPK.initialize)."""
        self._scoped_jacobian_requested = bool(enable)
        self._setup_scoped_jacobian()
        return self._scoped_jacobian_active

    def _setup_scoped_jacobian(self):
        """Arm or clear the pydas `jacobian` capability attribute.

        Guard 5 (fail-safe fallback): every layout assumption the scoped
        FD sweep relies on is validated here against the LIVE model; any
        mismatch clears the capability (DASPK falls back to its internal
        dense FD -- today's behavior) with ONE loud warning naming the
        failed assumption. When the feature is not requested this method
        only guarantees the attribute is absent (default-off path:
        hasattr(self, 'jacobian') stays False, bit-identical to today).
        """
        cdef int n_core, neq_expected
        self._scoped_jacobian_active = False
        self._jac_scope = False
        # Default OFF: make sure pydas can see no capability attribute.
        self.__dict__.pop("jacobian", None)
        if not getattr(self, "_scoped_jacobian_requested", False):
            return

        problems = []
        n_core = self.num_core_species
        neq_expected = n_core + self.num_qssa_u + self.num_sgh_z + self.num_end_ledger
        # Guard 3 baseline: appended U/Z slots must exactly close the
        # state vector (initiate_tolerances extended neq by ITS count).
        if self.neq != neq_expected:
            problems.append(
                f"neq {self.neq} != n_core+U+Z {neq_expected}")
        # DDAWTS weights need the full tolerance arrays.
        if (self.atol_array is None or self.rtol_array is None
                or len(self.atol_array) != self.neq
                or len(self.rtol_array) != self.neq):
            problems.append("atol/rtol arrays missing or wrong length")
        elif not np.all(np.asarray(self.atol_array) > 0.0):
            # wt_j > 0 for every component is what makes the DMATD
            # increment nonzero at y_j = dydt_j = 0.
            problems.append("non-positive atol entries (FD del could be 0)")
        # Edge rows must be the contiguous row suffix the scope drops.
        if (self.reactant_indices is None
                or self.reactant_indices.shape[0]
                != self.num_core_reactions + self.num_edge_reactions):
            problems.append("reactant_indices rows != core+edge count")
        # Sensitivity mode drives extra residual semantics pydas-side;
        # the scoped sweep has no sensitivity arm -- refuse.
        if getattr(self, "sensitivity", False):
            problems.append("sensitivity mode enabled")
        if self.gas_species_mask is None \
                or self.gas_species_mask.shape[0] != n_core:
            problems.append("gas_species_mask missing/mismatched")

        if problems:
            logging.warning(
                "18.2e scoped Jacobian REQUESTED but NOT armed -- "
                "falling back to DASPK internal FD (no user jacobian): %s",
                "; ".join(problems))
            return

        self._jac_wt_atol = np.array(self.atol_array, dtype=float)
        self._jac_wt_rtol = np.array(self.rtol_array, dtype=float)
        # Arm: instance attribute only (pydas dispatches by attribute
        # presence; the class itself never grows a `jacobian` member).
        self.jacobian = self._scoped_fd_jacobian
        self._scoped_jacobian_active = True
        logging.info(
            "18.2e scoped-RHS FD Jacobian ARMED (neq=%d, core rows=%d, "
            "edge rows skipped in J-sweeps=%d)",
            self.neq, self.num_core_reactions, self.num_edge_reactions)

    def _scoped_fd_jacobian(self, t, y, dydt, cj,
                            senpar=np.zeros(1, float)):
        """Dense one-sided FD iteration matrix J = dG/dy + cj*dG/dy',
        mirroring DDASPK's internal DMATD column increments exactly:

            del_j = sqrt(uround) * max(|y_j|, |h*y'_j|, |wt_j|)
            del_j = sign(del_j, h*y'_j);  del_j = (y_j + del_j) - y_j
            perturb y_j += del_j AND y'_j += cj*del_j, one column per
            residual call, column = (G_pert - G_base)/del_j

        with h approximated by 1/cj (DMATD receives the true step h; the
        callback signature only carries cj = alpha_s/h, alpha_s ~ O(1) --
        an increment-size choice, not a law change) and wt_j the DDAWTS
        weight rtol_j*|y_j| + atol_j.

        Each of the neq+1 residual calls runs the SAME residual() code
        path with _jac_scope=True, which prunes ONLY (i) edge rows and
        (ii) the two diagnostic censuses -- both proven residual-
        invisible (P4 probe: returned delta bitwise identical; pinned by
        the golden-J guard tests). In particular EVERY column -- gas,
        moment, proxy, U/Z -- differences the full core RHS including the
        V_gas = f(sum gas y) global-volume coupling and the pdep
        collider-efficiency recomputation (guard 1: no coupling column is
        omitted or approximated).
        """
        cdef int j, neq, n_core
        cdef double h, ysave, ypsave, hyp, delj

        neq = len(y)
        n_core = self.num_core_species
        # Guard 3 (shape assert, loud): the scoped J shape must track the
        # appended-U/Z state layout EXACTLY; a drifted layout must never
        # produce a silently wrong-shaped/misaligned J.
        if (neq != self.neq
                or neq != n_core + self.num_qssa_u + self.num_sgh_z + self.num_end_ledger
                or len(dydt) != neq
                or self._jac_wt_atol is None
                or len(self._jac_wt_atol) != neq
                or len(self._jac_wt_rtol) != neq):
            raise RuntimeError(
                "18.2e scoped Jacobian: state layout drifted from armed "
                f"assumptions (len(y)={neq}, self.neq={self.neq}, "
                f"n_core+U+Z+L={n_core + self.num_qssa_u + self.num_sgh_z + self.num_end_ledger})"
                " -- refusing to return a misaligned J.")
        if cj == 0.0:
            raise RuntimeError("18.2e scoped Jacobian: cj == 0")

        h = 1.0 / cj
        squr = np.sqrt(np.finfo(np.float64).eps)
        wt = self._jac_wt_rtol * np.abs(y) + self._jac_wt_atol

        ycol = np.array(y, dtype=float, copy=True)
        ypcol = np.array(dydt, dtype=float, copy=True)
        pd = np.empty((neq, neq), dtype=float)
        self._jac_scope = True
        try:
            # base residual (scoped): residual() returns a FRESH array
            # (dn_dt - dydt allocates), safe to hold across the sweep.
            base = np.asarray(
                self.residual(t, ycol, ypcol, senpar)[0])
            for j in range(neq):
                ysave = ycol[j]
                ypsave = ypcol[j]
                hyp = h * ypsave
                delj = squr * max(abs(ysave), abs(hyp), abs(wt[j]))
                # Fortran SIGN(del, hyp): |del| if hyp >= 0 else -|del|
                if hyp < 0.0:
                    delj = -delj
                # representable increment, exactly as DMATD forms it
                delj = (ysave + delj) - ysave
                if delj == 0.0:
                    # unreachable while wt_j > 0 (validated at arming);
                    # defensive: never divide by zero below
                    delj = squr * wt[j] if wt[j] > 0.0 else squr
                ycol[j] = ysave + delj
                ypcol[j] = ypsave + cj * delj
                col = np.asarray(
                    self.residual(t, ycol, ypcol, senpar)[0])
                pd[:, j] = (col - base) / delj
                ycol[j] = ysave
                ypcol[j] = ypsave
        finally:
            self._jac_scope = False
            # Guard-2 tripwire: diagnostics (edge_* arrays, censuses) are
            # stale until the next REAL residual call clears this.
            self._scoped_jac_diag_stale = True
        return pd

    @cython.boundscheck(False)
    def residual(self, double t, np.ndarray[np.float64_t, ndim=1] y,
                 np.ndarray[np.float64_t, ndim=1] dydt,
                 np.ndarray[np.float64_t, ndim=1] senpar = np.zeros(1, float)):
        """
        Compute the residual (dn/dt - dydt) at time t for state y.
        """
        cdef int i
        cdef int n_core = len(self.core_species_rates)
        cdef int n_rxn = len(self.core_reaction_rates)
        cdef int p0_idx, p1_idx, p2_idx, prod_p_idx, p_slot, p_idx_tmp
        cdef int n_rows_scope
        cdef bint jac_scope
        cdef double mu1_0, mu1_1, rf_hijack

        # 18.2e P5: J-scope flag -- True ONLY inside _scoped_fd_jacobian's
        # FD sweep. It prunes work PROVEN residual-invisible (edge rows:
        # every dn_dt/du_dt/dz_dt write below is core_rxn-gated; the two
        # exhaustion/SGH censuses: diagnostic-only per r86/§1g) and
        # changes NOTHING else -- same code path, same laws. Every real
        # (DASPK-originated) call runs with the flag False.
        jac_scope = self._jac_scope

        if self.gas_species_mask.shape[0] != n_core:
            raise ValueError(f"State/Mask mismatch: y={n_core}, mask={self.gas_species_mask.shape[0]}")

        # Disable caching
        self.jacobian_matrix = None

        ir = self.reactant_indices
        ip = self.product_indices
        inet = self.network_indices
        knet = self.network_leak_coefficients

        n_core = len(self.core_species_rates)
        n_rxn = len(self.core_reaction_rates)
        n_edge = len(self.edge_species_rates)

        # Item 17 (spec 2026-06-12 SS3(b)): combined-length raise mirroring
        # the gas_species_mask raise above -- the product gates below index
        # prospective_gas_mask with raw product slots (core AND edge).
        if (self.prospective_gas_mask is None
                or self.prospective_gas_mask.shape[0] != n_core + n_edge):
            raise ValueError(
                f"State/Prospective-mask mismatch: n_core+n_edge="
                f"{n_core + n_edge}, prospective="
                f"{None if self.prospective_gas_mask is None else self.prospective_gas_mask.shape[0]}")

        # 1. Lazy Allocation (Safety check)
        if self._scratch_C_gas is None or len(self._scratch_C_gas) != n_core:
            self._scratch_C_gas = np.zeros(n_core, float)
            self._scratch_C_poly = np.zeros(n_core, float)
            self._scratch_dn_dt = np.zeros(n_core, float)
        if (self._scratch_du_dt is None
                or len(self._scratch_du_dt) != self.num_qssa_u):
            self._scratch_du_dt = np.zeros(self.num_qssa_u, float)
        if (self._scratch_dz_dt is None
                or len(self._scratch_dz_dt) != self.num_sgh_z):
            self._scratch_dz_dt = np.zeros(self.num_sgh_z, float)
        if (self._scratch_dl_dt is None
                or len(self._scratch_dl_dt) != self.num_end_ledger):
            self._scratch_dl_dt = np.zeros(self.num_end_ledger, float)

        # 2. Update Volumes
        if self.constant_gas_volume:
            V_gas = self.V_gas0
        else:
            n_gas = float(np.sum(y[:n_core][self.gas_species_mask]))
            V_gas = constants.R * self.T.value_si * n_gas / self.P.value_si if n_gas > 0 else 1.0
        self.V_gas = V_gas
        V_poly = self.V_poly
        self.V = self.V_gas + self.V_poly

        # 3. Concentrations (Vectorized)
        C_gas = self._scratch_C_gas
        C_poly = self._scratch_C_poly
        C_gas[:] = 0.0
        C_poly[:] = 0.0

        mask = self.gas_species_mask
        pmask = self.prospective_gas_mask  # gate-input only (rider R3)
        # NOTE the [:n_core] slice: y may carry trailing weak-link U slots
        # beyond the core-species block; the mask covers only the prefix.
        # M18.2d (round-40 one-shot): the concentration clamp is the widest
        # kink field the corrector differentiates (one {y_i = 0} kink plane
        # per species). softclamp extends it C1 into NEGATIVE trial-state
        # territory only; for y >= 0 (all valid states) the fast path below
        # is bit-for-bit the old np.maximum law (positives untouched).
        _y_head = y[:n_core]
        _neg = _y_head < 0.0
        if _neg.any():
            # evaluate the extension ONLY on the negative entries (no exp
            # overflow on large positives; positives stay untouched views)
            _y_head = _y_head.copy()
            _yn = _y_head[_neg]
            _y_head[_neg] = _yn * np.exp(_yn / self._softclamp_lam)
        C_gas[mask] = _y_head[mask] / V_gas
        C_poly[~mask] = _y_head[~mask] / V_poly

        # Sync diagnostic concentrations
        self.core_species_concentrations[mask] = C_gas[mask]
        self.core_species_concentrations[~mask] = C_poly[~mask]

        # 4. Clear Accumulators
        dn_dt = self._scratch_dn_dt
        dn_dt[:] = 0.0

        du_dt = self._scratch_du_dt
        du_dt[:] = 0.0

        dz_dt = self._scratch_dz_dt
        dz_dt[:] = 0.0

        dl_dt = self._scratch_dl_dt
        dl_dt[:] = 0.0

        proxy_activity = self._scratch_proxy_activity
        proxy_activity[:] = 0.0

        # r81 (B) exhaustion census/tripwire, re-adjudicated r86: refresh
        # the per-pool exhausted flags from the RAW state ONCE per residual
        # evaluation (hard-errors on any moment beyond -floor; census once
        # per pool per rebuild). The flags are DIAGNOSTIC ONLY -- no RHS
        # read below consults them (the r81 read-projection is reverted);
        # site factors, bundles and kernels always compute from raw y.
        # 18.2e P5: skipped under jac_scope (diagnostic-only; guard-2
        # tripwire below pins that a REAL call always follows a J-sweep
        # before any consumer reads these diagnostics).
        if not jac_scope:
            self._update_pool_exhaustion(y)

        # SGH kernel-v2 Z-inventory census (§1g, P1-A): TRIAL-STATE-SAFE --
        # residual() is called on the raw Newton trial state, so this NEVER
        # raises and NEVER mutates y. It only census/logs (warn-once) any Z
        # sitting outside the physical band [0, sites*mu1]; the RHS's local
        # max(0.0, .) reads sanitize a negative trial Z. The HARD physical-
        # band invariant is enforced on the ACCEPTED state only, by
        # _assert_sgh_inventory_accepted (self.y), invoked once per accepted
        # snapshot from _phase_gate_flux_census. Cheap no-op when
        # num_sgh_z == 0.
        # 18.2e P5: skipped under jac_scope (census-only); real calls also
        # clear the guard-2 staleness tripwire here -- after any step()/
        # advance() the diagnostics were refreshed by a REAL call.
        if not jac_scope:
            self._update_sgh_inventory_guard(y)
            self._scoped_jac_diag_stale = False

        self.core_reaction_rates[:] = 0.0
        self.edge_reaction_rates[:] = 0.0
        self.edge_species_rates[:] = 0.0
        self.edge_reaction_gate_code[:] = 0
        self.edge_reaction_rates_ungated[:] = 0.0
        self.edge_species_rates_ungated[:] = 0.0
        self.core_species_consumption_rates[:] = 0.0
        self.core_species_production_rates[:] = 0.0
        self.network_leak_rates[:] = 0.0

        # 5. Recalculate P-dependent rates in the GAS phase
        if self.pdep_collision_reaction_indices.shape[0] != 0:
            T_val = self.T.value_si
            P_val = self.P.value_si

            y_gas = y[:n_core] * mask
            sum_y_gas = np.sum(y_gas)

            if sum_y_gas > 1e-20:
                for i in range(self.pdep_collision_reaction_indices.shape[0]):
                    j = self.pdep_collision_reaction_indices[i]

                    Peff = P_val * np.sum(self.collider_efficiencies[i] * y_gas / sum_y_gas)

                    kf_j = self.pdep_collider_kinetics[i].get_rate_coefficient(T_val, Peff)
                    self.kf[j] = kf_j if np.isfinite(kf_j) else 0.0
                    if self.Keq[j] != np.inf and self.Keq[j] > 0.0:
                        kb_j = self.kf[j] / self.Keq[j]
                        self.kb[j] = kb_j if np.isfinite(kb_j) else 0.0
            else:
                pass

        # 6. Standard Kinetics (Optimized Loop)
        # 18.2e P5: under jac_scope only the n_rxn CORE rows run. Edge
        # rows (r_idx >= n_rxn; contiguous suffix by construction in
        # generate_reactant_product_indices) write ONLY edge_*/proxy-
        # activity diagnostics, never dn_dt/du_dt/dz_dt (every state
        # write is core_rxn-gated) -- P4 probe: the returned delta is
        # BITWISE identical with edge rows skipped, and the golden-J
        # guard tests pin it.
        n_rows_scope = n_rxn if jac_scope else ir.shape[0]
        for r_idx in range(n_rows_scope):
            r0, r1, r2 = ir[r_idx, 0], ir[r_idx, 1], ir[r_idx, 2]

            if r0 == -1:
                continue

            # simple.pyx parity (:443-459, :509-515): a row with an EDGE
            # REACTANT is still evaluated -- the edge species has no state in
            # y, so its concentration is ZERO and rf vanishes; the discovery
            # flux for reverse-stored decomposition channels (daughter as
            # edge reactant, rr = kb*[core products]) flows through the
            # reactant-side edge accumulator in section 4b. Phase verdict for
            # an edge reactant comes from prospective_gas_mask (the only mask
            # covering edge slots; its core prefix is bit-identical to
            # gas_species_mask by rider R1).
            r0_is_gas = mask[r0] if r0 < n_core else pmask[r0]
            r1_is_gas = True
            r2_is_gas = True

            if r1 != -1:
                r1_is_gas = mask[r1] if r1 < n_core else pmask[r1]

            if r2 != -1:
                r2_is_gas = mask[r2] if r2 < n_core else pmask[r2]

            is_poly_event = (not r0_is_gas) or (r1 != -1 and not r1_is_gas) or (r2 != -1 and not r2_is_gas)

            if is_poly_event:
                phase = 'poly'
                V_rxn = V_poly
            else:
                phase = 'gas'
                V_rxn = V_gas

            # Product Checks
            p0, p1, p2 = ip[r_idx, 0], ip[r_idx, 1], ip[r_idx, 2]

            prods_phase_ok = True
            has_edge_prod = False
            has_any_prod = False
            has_condensed_prod = False
            gate_code = 0  # 0 = open, 1 = Gate A, 2 = Gate B

            # PROSPECTIVE-MASK GATES (item 17, spec 2026-06-12 SS3(c)). The
            # phase verdict for EVERY product, core or edge alike, is
            # prospective_gas_mask[p]; for p < n_core it is bit-identical to
            # mask[p] by rider R1 (core-prefix parity raise at init), so
            # CORE-product behavior is unchanged by construction. The old
            # has_edge_prod PHASE bypass (maskless edge product -> exempt)
            # is dead: it is what let enlargement promote on flux the
            # post-promotion model zeroes (the umbrella invariant's mask
            # projection). has_edge_prod survives ONLY for the reverse-rate
            # concentration-availability hole below (:rr block) -- an edge
            # product has no state in y, so rr is UNCOMPUTABLE there; that
            # is simple.pyx parity (Z6), not a phase verdict.
            def _check_prod(p):
                nonlocal has_edge_prod, has_any_prod, has_condensed_prod
                if p == -1:
                    return
                has_any_prod = True
                if p >= n_core:
                    has_edge_prod = True
                if not pmask[p]:
                    has_condensed_prod = True

            _check_prod(p0)
            _check_prod(p1)
            _check_prod(p2)

            if (not is_poly_event) and has_condensed_prod:
                # Gate A: a gas event with ANY prospectively-condensed
                # product (core or edge) is zeroed.
                prods_phase_ok = False
                gate_code = 1
            elif is_poly_event and (not has_condensed_prod):
                # Gate B: a poly event with NO prospectively-condensed
                # product (core or edge) is zeroed.
                prods_phase_ok = False
                gate_code = 2

            gated = False
            if not prods_phase_ok:
                if r_idx >= n_rxn:
                    # RIDER R2 dynamic half: record the gate verdict and fall
                    # through the EXISTING rate computation (pool mapping,
                    # _C, site scaling, throttle; the rr hole applies as
                    # usual) with every REAL write suppressed -- the zeros in
                    # edge_reaction_rates/edge_species_rates ARE the
                    # consistency point; the counterfactual lands only in
                    # the *_ungated arrays.
                    self.edge_reaction_gate_code[r_idx - n_rxn] = gate_code
                    gated = True
                else:
                    # Core-gated rows keep the bare continue: zero residual
                    # cost. Their loudness is the STATIC census at
                    # initialize_model (spec SS3(e) static half).
                    continue

            # 1. Map Reactants to Polymer Pools  (MOVED UP before _C)
            #    Edge reactants map to -1: species_to_pool_indices covers the
            #    core prefix only, and an edge species is never a pool proxy.
            p0_pool_idx = -1 if r0 >= n_core else self.species_to_pool_indices[r0]
            p1_pool_idx = -1 if (r1 == -1 or r1 >= n_core) else self.species_to_pool_indices[r1]
            p2_pool_idx = -1 if (r2 == -1 or r2 >= n_core) else self.species_to_pool_indices[r2]  # optional but consistent

            # Rate Calculation
            kf = self.kf[r_idx]
            kb = self.kb[r_idx]

            def _C(idx):
                if idx == -1:
                    return 1.0
                if idx >= n_core:
                    # Edge species: no state in y -> zero concentration, so
                    # rf vanishes for edge-reactant rows (simple.pyx parity;
                    # NOT a guard -- definitional, no silent zeroing).
                    return 0.0
                if self.species_to_pool_indices[idx] != -1:
                    return 1.0
                return C_gas[idx] if mask[idx] else C_poly[idx]

            if r1 == -1:
                rf = kf * _C(r0)
            elif r2 == -1:
                rf = kf * _C(r0) * _C(r1)
            else:
                rf = kf * _C(r0) * _C(r1) * _C(r2)

            # Reverse Rate (Unrolled, Hole-Safe)
            rr = 0.0
            if has_any_prod and not has_edge_prod:
                rr = kb
                if p0 != -1: rr *= _C(p0)
                if p1 != -1: rr *= _C(p1)
                if p2 != -1: rr *= _C(p2)

            # 2. Polymer scaling
            if p0_pool_idx != -1 or p1_pool_idx != -1 or p2_pool_idx != -1:
                if p0_pool_idx != -1 or p1_pool_idx != -1 or p2_pool_idx != -1:
                    # Determine the target pool index
                    target_pool_idx = p0_pool_idx
                    if target_pool_idx == -1: target_pool_idx = p1_pool_idx
                    if target_pool_idx == -1: target_pool_idx = p2_pool_idx

                    # Default to Mu1 (Site Density)
                    moment_idx = self.polymer_pools[target_pool_idx].mu_indices[1]

                    # Check for specific "End Group" physics (Tier 1 / Initiation)
                    if self.is_end_group_reaction[r_idx]:
                        moment_idx = self.polymer_pools[target_pool_idx].mu_indices[0]  # Scale by Mu0 (Chain Density)

                    # Calculate the site/chain concentration
                    # Note: We use the *mapped* index from the pool, ensuring we grab the real state variable
                    site = _softclamp(y[moment_idx], self._softclamp_lam) / V_poly

                    # Exhaustion throttle (spec 2026-06-10 s5 AMENDMENT): a
                    # mu0-scaled DISCRETE_CHIP drains mu1 but never mu0, so
                    # its rate would not self-limit -- mu1 would run linearly
                    # negative past exhaustion while chips keep being created.
                    # Counting inequality (not a heuristic): every chain with
                    # n >= a holds at least a units, so donatable ends are
                    # bounded by mu1/a always. Throttling the SITE (not
                    # write-gating the mu1 leg, which would duplicate mass)
                    # gives dmu1/dt >= -kf*mu1 globally; it multiplies BOTH
                    # directions by design. a = 0 chips (drain nothing) and
                    # mu1-scaled chips (already self-limiting) are exempt.
                    # (mirrored in get_reaction_rates' hijack block -- keep in sync)
                    if (self.reaction_flux_archetype[r_idx] == FLUX_DISCRETE_CHIP
                            and self.is_end_group_reaction[r_idx]
                            and self.reaction_chip_units[r_idx] > 0):
                        mu_idx = self.polymer_pools[target_pool_idx].mu_indices
                        site = min(
                            _softclamp(y[mu_idx[0]], self._softclamp_lam),
                            _softclamp(y[mu_idx[1]], self._softclamp_lam) / float(self.reaction_chip_units[r_idx]),
                        ) / V_poly
                    elif (self.reaction_flux_archetype[r_idx] == FLUX_VOLATILE_EJECTION
                            and self.is_end_group_reaction[r_idx]
                            and self.reaction_src_pool[r_idx] != -1
                            and self.reaction_eject_units_moment[r_idx] > 0.0):
                        # a>0 VE shares the DISCRETE_CHIP exhaustion
                        # structure: mu0-scaled drain of mu1 that never touches
                        # mu0, so it would run mu1 linearly negative past
                        # exhaustion. Throttle site = min(mu0, mu1/a). Guard
                        # a>0: a<0 GROWS the chain (no exhaustion) and mu1/a<0
                        # would be a spurious negative site. Ruling round 20
                        # item C extends this from same-pool to CROSS-POOL
                        # rows too (the forward leg debits src and sheds a
                        # units/event it must be able to source; previously
                        # only the b0==0 skip guarded it). INTERIOR
                        # (mu1-scaled) rows stay unthrottled, relying on the
                        # linear self-limit of the mu1 site -- TODO (Codex
                        # round-22 P2b): add a near-exhaustion stress-pin
                        # test for that reliance. (mirrored in
                        # get_reaction_rates' hijack block -- keep in sync)
                        mu_idx = self.polymer_pools[target_pool_idx].mu_indices
                        site = min(
                            _softclamp(y[mu_idx[0]], self._softclamp_lam),
                            _softclamp(y[mu_idx[1]], self._softclamp_lam) / float(self.reaction_eject_units_moment[r_idx]),
                        ) / V_poly

                    # P1-1 DIRECTIONAL BUNDLE THROTTLE, round-27 P1-A
                    # tail-only form (cross-pool VE/MIGRATION legs; closes
                    # the round-22 P2b TODO -- its "linear self-limit of the
                    # mu1 site" premise is FALSE for the reverse leg, whose
                    # per-event mu2/mu1 debit cancels the mu1 scaling and
                    # drains mu1 at constant rate once the frozen
                    # out-of-cone closure stops debiting mu2; regen-#2
                    # forensics). Each debited direction's event rate runs
                    # at the near-exhaustion bundle limiter (tail-only
                    # smoothstep; C1 soft-min, cone-margin drain guard)
                    # site -- two independent gates in series:
                    #   S_free = w*S_base + (1-w)*softmin_p(S_base, S_cap),
                    #   S_cap  = softmin_p(mu0/b0, mu1/b1, mu2/b2)/V_poly
                    #            over the POSITIVE bundle terms of the
                    #            POOL BEING DEBITED,
                    #   then the round-30 cone-margin drain gate on the
                    #   margin Q10 = mu1 - mu0 (only for b1 > b0 debits;
                    #   INDEPENDENT of E): Q10 <= 0 -> 0 regardless of E;
                    #   otherwise its own M = Q10/max(f0,f1) gate ramps from
                    #   I-090's narrowed completion to the original soft-min
                    #   cap immediately below M_LO, then blends into S_free,
                    #   S_cone = Q10/(V_poly*(b1 - b0)),
                    #   w = smoothstep of the debited pool's accepted-state
                    #       floor distance E (see _bundle_limited_site) --
                    # S_eff == S_base EXACTLY in bulk (E >= E_hi; healthy
                    # pools keep the adjudicated site law bit-for-bit), and
                    # the hard cap is fully active near exhaustion
                    # (E <= E_lo). Applied BEFORE rate = rf - rr so species
                    # writes, core_reaction_rates, edge promotion and the
                    # section-5 moment dispatch all see the SAME limited
                    # directional rate -- gas flux and moment flux must
                    # never diverge. Forward debits src; the reverse limiter
                    # is applied in the dst block below. (mirrored in
                    # get_reaction_rates -- keep in sync)
                    cross_bundle_throttle = False
                    if self.reaction_flux_archetype[r_idx] in (
                            FLUX_MIGRATION, FLUX_VOLATILE_EJECTION):
                        src_bt = self.reaction_src_pool[r_idx]
                        dst_bt = self.reaction_dst_pool[r_idx]
                        if src_bt != -1 and dst_bt != -1 and src_bt != dst_bt:
                            cross_bundle_throttle = True
                            site = self._bundle_limited_site(
                                src_bt, y, V_poly,
                                self.is_end_group_reaction[r_idx], site)

                    rf *= site
                    # Direction-specific source availability (run-5 DASPK
                    # IDID=-7 forensics, adjudicated Part C): each direction's
                    # event flux must scale with the moments of the pool
                    # ACTUALLY BEING DEBITED in that direction. The REVERSE
                    # leg of a cross-pool exchange debits the dst pool
                    # (section-5 legs run (rf, src->dst), (rr, dst->src)), so
                    # its site factor comes from the dst pool's OWN moments --
                    # with the forward-reactant site a near-empty spawned pool
                    # was drained at a rate set by the healthy source pool
                    # (mu2 negative in ~1e-25 s, corrector divergence,
                    # h -> ~4e-21, IDID=-7). Same moment order as the forward
                    # site (mu1; mu0 for end-group rows), so the outflux -> 0
                    # CONTINUOUSLY (linearly) as the debited pool empties: no
                    # dt-dependent clamp, no activation threshold. Same-pool
                    # rows (incl. the throttled DISCRETE_CHIP / same-pool-VE
                    # shapes above, whose throttle multiplies BOTH directions
                    # by design) and rows with no product-side pool (dst ==
                    # -1: the reverse leg debits gas/explicit species only,
                    # whose availability is already in _C) keep the shared
                    # reactant-pool site. (mirrored in get_reaction_rates'
                    # hijack block -- keep in sync)
                    dst_pool_idx = self.reaction_dst_pool[r_idx]
                    if dst_pool_idx != -1 and dst_pool_idx != target_pool_idx:
                        moment_idx = self.polymer_pools[dst_pool_idx].mu_indices[1]
                        if self.is_end_group_reaction[r_idx]:
                            moment_idx = self.polymer_pools[dst_pool_idx].mu_indices[0]
                        if (self.reaction_flux_archetype[r_idx] == FLUX_VOLATILE_EJECTION
                                and self.is_end_group_reaction[r_idx]
                                and self.reaction_eject_units_moment[r_idx] < 0.0):
                            # Ruling round 20 item C, mirror direction: for
                            # a<0 the REVERSE leg of a cross-pool VE sheds
                            # |a| units/event from the dst pool it debits
                            # (sa = +a < 0 shrinks the chain landing back
                            # in src), so its availability carries the same
                            # min(mu0, mu1/|a|) exhaustion throttle as the
                            # a>0 forward site above. The forward leg (chain
                            # GROWTH landing in dst) stays exempt. (mirrored
                            # in get_reaction_rates -- keep in sync)
                            mu_idx_dst = self.polymer_pools[dst_pool_idx].mu_indices
                            site_rev = min(
                                _softclamp(y[mu_idx_dst[0]], self._softclamp_lam),
                                _softclamp(y[mu_idx_dst[1]], self._softclamp_lam) / (-float(self.reaction_eject_units_moment[r_idx])),
                            ) / V_poly
                        else:
                            site_rev = _softclamp(y[moment_idx], self._softclamp_lam) / V_poly
                        # P1-1 mirror direction (round-27 P1-A tail-only
                        # form): the reverse leg debits the dst pool's FULL
                        # per-chain bundle, so it carries the same
                        # near-exhaustion bundle limiter (tail-only
                        # smoothstep; C1 soft-min, cone-margin drain guard)
                        # S_eff as the forward leg -- computed from the dst
                        # pool's OWN moments, floor distance and cone
                        # margin. (mirrored in get_reaction_rates --
                        # keep in sync)
                        if cross_bundle_throttle:
                            site_rev = self._bundle_limited_site(
                                dst_pool_idx, y, V_poly,
                                self.is_end_group_reaction[r_idx], site_rev)
                        rr *= site_rev
                    else:
                        rr *= site
            elif has_any_prod and not has_edge_prod:
                # PRODUCT-side proxy (review round 51): a reverse-stored
                # association row (e.g. edge daughter + edge daughter <=>
                # pool proxy, or gas radical + gas radical <=> pool proxy)
                # carries the proxy on the PRODUCT side, so rr's _C(proxy)
                # placeholder (1.0) must be replaced by the pool site
                # density, mirroring the reactant-side hijack above --
                # otherwise reverse homolysis discovery flux is unscaled.
                # Products are all core here (has_edge_prod is False), so
                # species_to_pool_indices indexing is safe. rf is untouched:
                # no reactant is a proxy on this branch (elif), so the
                # forward direction has no proxy concentration to replace.
                prod_pool_idx = -1
                for prod_p_idx in (p0, p1, p2):
                    if prod_p_idx != -1 and self.species_to_pool_indices[prod_p_idx] != -1:
                        prod_pool_idx = self.species_to_pool_indices[prod_p_idx]
                        break
                if prod_pool_idx != -1:
                    # Default to Mu1 (site density); end-group physics
                    # scales by Mu0 (chain density) -- same selection as
                    # the reactant-side block.
                    moment_idx = self.polymer_pools[prod_pool_idx].mu_indices[1]
                    if self.is_end_group_reaction[r_idx]:
                        moment_idx = self.polymer_pools[prod_pool_idx].mu_indices[0]
                    rr *= _softclamp(y[moment_idx], self._softclamp_lam) / V_poly

            # Net rate (volumetric, in the phase volume chosen earlier)
            rate = rf - rr
            r_mol_s = rate * V_rxn
            abs_flux = abs(r_mol_s)

            if r_idx < n_rxn:
                self.core_reaction_rates[r_idx] = rate
            elif gated or self.reaction_refused[r_idx]:
                # counterfactual only -- the real entry stays 0 (the
                # consistency point enlargement reads). r71 FIX 3: refused
                # edge rows are routed here too -- their real
                # edge_reaction_rates entry must be zero (enlargement
                # input), the ungated entry stays a strictly-diagnostic
                # counterfactual (census only, never flux).
                self.edge_reaction_rates_ungated[r_idx - n_rxn] = rate
            else:
                self.edge_reaction_rates[r_idx - n_rxn] = rate
                self.edge_reaction_rates_ungated[r_idx - n_rxn] = rate

            core_rxn = r_idx < n_rxn

            if self.reaction_refused[r_idx]:
                # item 18 + r71 FIX 3: a refused reaction must contribute NO
                # flux ANYWHERE (stamp-but-keep): skip all reactant/product/
                # moment writes so no MW-weighted backbone mass is fabricated
                # in core -- AND no edge promotion flux is delivered. The
                # pre-r71 core-only gate left refused EDGE rows feeding
                # edge_species_rates / proxy_activity, which is exactly how
                # PP run-5's dead chain-scale radicals were promoted to core
                # (forensics: run5/RMG.log:910-955), where the adjudication
                # was lost and legacy mu1 flux collapsed DASPK. Refused rows
                # must not move species to core; the per-reaction ungated
                # counterfactual above is the only surviving record
                # (diagnostic, never flux).
                continue

            # 3. Apply Fluxes to Reactants (core reactions only -- edge
            #    reactions are diagnostic-only, matching simple.pyx).
            #    proxy_activity stays UNGATED on purpose: it feeds proxy
            #    similarity/spawn diagnostics, which want edge flux too --
            #    except gate-zeroed edge rows, suppressed just below (T9).
            #    Pool MOMENT flux is handled once per reaction in section 5.
            if r0 >= n_core:
                # Edge reactant: no core state, never a proxy; its
                # reactant-side edge flux is accumulated in section 4b
                # (diagnostic-only).
                pass
            elif self.is_pool_proxy[r0]:
                if not gated:
                    # Ghost-flux suppression (spec SS3(e)/T9): a gate-zeroed
                    # edge row's |flux| is judged unphysical under core
                    # semantics and must not feed spawn/similarity
                    # diagnostics.
                    proxy_activity[r0] += abs_flux
                if core_rxn:
                    self.core_species_consumption_rates[r0] += rf
                    self.core_species_production_rates[r0] += rr
            elif core_rxn:
                dn_dt[r0] -= r_mol_s
                # Change (a) (spec 2026-06-10 mass-flux-spawn-gate s3.1/s4.6):
                # gross bookkeeping for ORDINARY core species too (simple.pyx
                # parity), so absorbed-representative production has a real
                # record for spawn_gate_flux_snapshot(). Cost-gated on the
                # EPDM deck (<= ~5% wall-clock; fallback: ledger-tracked-only
                # writes via a flag array).
                self.core_species_consumption_rates[r0] += rf
                self.core_species_production_rates[r0] += rr

            if r1 != -1:
                if r1 >= n_core:
                    pass  # edge reactant: see the r0 branch above
                elif self.is_pool_proxy[r1]:
                    if not gated:
                        proxy_activity[r1] += abs_flux
                    if core_rxn:
                        self.core_species_consumption_rates[r1] += rf
                        self.core_species_production_rates[r1] += rr
                elif core_rxn:
                    dn_dt[r1] -= r_mol_s
                    self.core_species_consumption_rates[r1] += rf
                    self.core_species_production_rates[r1] += rr

            if r2 != -1 and core_rxn:
                # r2 is always ordinary today (mirror the proxy logic above
                # if polymer is ever allowed in r2); change (a) gross writes
                # apply as for any ordinary reactant.
                dn_dt[r2] -= r_mol_s
                self.core_species_consumption_rates[r2] += rf
                self.core_species_production_rates[r2] += rr

            # 4. Apply Fluxes to Products (same core_rxn gate as section 3)
            for p_slot in range(3):
                p_idx_tmp = ip[r_idx, p_slot]
                if p_idx_tmp == -1:
                    continue

                if p_idx_tmp < n_core:
                    if self.is_pool_proxy[p_idx_tmp]:
                        if not gated:
                            proxy_activity[p_idx_tmp] += abs_flux
                        if core_rxn:
                            self.core_species_production_rates[p_idx_tmp] += rf
                            self.core_species_consumption_rates[p_idx_tmp] += rr
                    elif core_rxn:
                        # feature polymer species (e.g., PS_rad, scission oligomer, etc.) stays explicit
                        dn_dt[p_idx_tmp] += r_mol_s
                        # Change (a): gross bookkeeping for ordinary products
                        # (simple.pyx parity; see the reactant-side comment).
                        self.core_species_production_rates[p_idx_tmp] += rf
                        self.core_species_consumption_rates[p_idx_tmp] += rr
                else:
                    self.edge_species_rates_ungated[p_idx_tmp - n_core] += rate
                    if not gated:
                        self.edge_species_rates[p_idx_tmp - n_core] += rate

            # 4b. Reactant-side edge accumulation (simple.pyx :509-515 sign
            #     parity): edge reactants get -= rate, so a reverse-stored
            #     decomposition row (rf = 0, net rate = -rr) registers
            #     POSITIVE daughter production for enlargement. Diagnostic-
            #     only: no dn_dt / core bookkeeping (edge species have no
            #     state), exactly like the product-side branch above.
            for p_slot in range(3):
                p_idx_tmp = ir[r_idx, p_slot]
                if p_idx_tmp == -1:
                    continue
                if p_idx_tmp >= n_core:
                    self.edge_species_rates_ungated[p_idx_tmp - n_core] -= rate
                    if not gated:
                        self.edge_species_rates[p_idx_tmp - n_core] -= rate

            # 5. Pool moment flux -- archetype dispatch (core reactions only).
            #    Spec: docs/superpowers/specs/2026-06-09-proxy-moment-flux-
            #    apportionment-design.md. SAME_POOL and NONE apply nothing:
            #    fold-back flux is net-zero by construction, so skipping avoids
            #    roundoff and closure calls.
            if core_rxn:
                arch = self.reaction_flux_archetype[r_idx]
                if arch == FLUX_MIGRATION:
                    src = self.reaction_src_pool[r_idx]
                    dst = self.reaction_dst_pool[r_idx]
                    # -1 cannot reach here (init demotes unresolved pools to
                    # UNRESOLVED); the checks are defensive only.
                    if src != -1 and dst != -1 and src != dst:
                        # Per-direction bundles: forward moves A-statistics
                        # chains A->B, reverse moves B-statistics chains B->A.
                        # Each direction is guarded by its OWN source pool.
                        for ev_rate, from_pool, to_pool in (
                                (rf, src, dst), (rr, dst, src)):
                            if ev_rate <= 0.0:
                                continue
                            ev_mol = ev_rate * V_rxn
                            b0, b1, b2, mu2_ok = self._chain_bundle(
                                from_pool, y, V_poly,
                                self.is_end_group_reaction[r_idx])
                            if b0 == 0.0:
                                continue
                            f_idx = self.polymer_pools[from_pool].mu_indices
                            t_idx = self.polymer_pools[to_pool].mu_indices
                            dn_dt[f_idx[0]] -= ev_mol * b0
                            dn_dt[f_idx[1]] -= ev_mol * b1
                            dn_dt[t_idx[0]] += ev_mol * b0
                            dn_dt[t_idx[1]] += ev_mol * b1
                            if mu2_ok:
                                dn_dt[f_idx[2]] -= ev_mol * b2
                                dn_dt[t_idx[2]] += ev_mol * b2
                elif arch == FLUX_VOLATILE_EJECTION:
                    # Volatile ejection (spec 2026-06-2x): a MIGRATION mirror
                    # A(pool) <=> volatile_gas + B(other pool). The moved chain
                    # loses `a` = reaction_eject_units backbone units to the
                    # discrete gas volatile as it lands in the destination
                    # pool, so the to-leg carries the SAME b0 chains but with
                    # E[n] shifted by sa (sa = -a forward: chain shrinks
                    # landing in dst; +a reverse: chain re-grows landing in
                    # src). from-leg loses the FULL bundle (identical to
                    # MIGRATION). mu2 shift is the exact quadratic
                    # (n + sa)^2 - n^2 = 2*sa*n + sa^2; sa*sa == a*a in both
                    # directions -> (b2 - 2a*b1 + a^2*b0) forward,
                    # (b2 + 2a*b1 + a^2*b0) reverse. The gas volatile itself is
                    # credited by the standard section-4 net-rate product path
                    # (NO gas moles added inside this branch).
                    src = self.reaction_src_pool[r_idx]
                    dst = self.reaction_dst_pool[r_idx]
                    a = float(self.reaction_eject_units_moment[r_idx])
                    # -1 cannot reach here (init demotes unresolved pools);
                    # the checks are defensive only (mirror MIGRATION).
                    if src != -1 and dst != -1 and src != dst:
                        for ev_rate, from_pool, to_pool, sa in (
                                (rf, src, dst, -a), (rr, dst, src, +a)):
                            if ev_rate <= 0.0:
                                continue
                            ev_mol = ev_rate * V_rxn
                            b0, b1, b2, mu2_ok = self._chain_bundle(
                                from_pool, y, V_poly,
                                self.is_end_group_reaction[r_idx])
                            if b0 == 0.0:
                                continue
                            f_idx = self.polymer_pools[from_pool].mu_indices
                            t_idx = self.polymer_pools[to_pool].mu_indices
                            # from-leg loses the FULL bundle
                            dn_dt[f_idx[0]] -= ev_mol * b0
                            dn_dt[f_idx[1]] -= ev_mol * b1
                            # to-leg gains the a-SHIFTED bundle
                            dn_dt[t_idx[0]] += ev_mol * b0
                            dn_dt[t_idx[1]] += ev_mol * (b1 + sa * b0)
                            if mu2_ok:
                                dn_dt[f_idx[2]] -= ev_mol * b2
                                dn_dt[t_idx[2]] += ev_mol * (
                                    b2 + 2.0 * sa * b1 + a * a * b0)
                    elif src != -1 and src == dst:
                        # SAME-POOL volatile ejection (spec round-13): unzip /
                        # depropagation (a>0 sheds mass) or monomer/radical
                        # addition (a<0 gains mass) that lands back in the SAME
                        # pool. A dedicated single-pool chip-style write, NOT
                        # the cross-pool per-direction cancellation above: the
                        # cancellation can push mu2 the wrong way near
                        # exhaustion, so we write mu2 with the DISCRETE_CHIP
                        # clamp (>0) instead. `a` (signed) was resolved above.
                        # mu0 is untouched (no chain created/destroyed); the gas
                        # volatile flows through the standard section-4 path.
                        b0, b1, _b2, _mu2_ok = self._chain_bundle(
                            src, y, V_poly, self.is_end_group_reaction[r_idx])
                        if b0 != 0.0:
                            s_idx = self.polymer_pools[src].mu_indices
                            e_n = b1
                            if rf > 0.0:
                                # Forward: Delta n = -a. Signed a: for a<0 this
                                # ADDS mass (chain growth) -- correct. Clamp the
                                # mu2 decrement at >0 (identical to
                                # DISCRETE_CHIP): 2a*E[n]-a^2 < 0 is unphysical
                                # per-chain for a>0 but reachable in expectation
                                # near chip-size mean length; for a<0 it is
                                # always <0, so the mu2 growth term is dropped
                                # here (the exact +extension lives on reverse).
                                rf_mol = rf * V_rxn
                                dn_dt[s_idx[1]] -= rf_mol * a
                                mu2_dec = 2.0 * a * e_n - a * a
                                if mu2_dec > 0.0:
                                    dn_dt[s_idx[2]] -= rf_mol * mu2_dec
                            if rr > 0.0:
                                # Reverse: EXACT extension form (n+a)^2 - n^2 =
                                # +(2a*E[n] + a^2); unconditionally positive for
                                # a>0 (no clamp), mirrors DISCRETE_CHIP reverse.
                                rr_mol = rr * V_rxn
                                dn_dt[s_idx[1]] += rr_mol * a
                                dn_dt[s_idx[2]] += rr_mol * (
                                    2.0 * a * e_n + a * a)
                elif arch == FLUX_SCISSION_FRAGMENT:
                    src = self.reaction_src_pool[r_idx]
                    dst = self.reaction_dst_pool[r_idx]
                    # -1 cannot reach here (init demotes unresolved pools).
                    # src == dst (fold-back + non-pool daughter) is skipped:
                    # outcome-identical to the legacy transfer (see the init
                    # demotion comment).
                    if src != -1 and dst != -1 and src != dst:
                        s_idx = self.polymer_pools[src].mu_indices
                        d_idx = self.polymer_pools[dst].mu_indices
                        # Length-biased parent-chain statistics; mirrors
                        # _chain_bundle's non-end-group branch (kept inline
                        # because the scission factors 1/2, 2/3, 1/3 differ
                        # per moment and per side).
                        mu0_p = _softclamp(y[s_idx[0]], self._softclamp_lam) / V_poly
                        mu1_p = _softclamp(y[s_idx[1]], self._softclamp_lam) / V_poly
                        mu2_p = _softclamp(y[s_idx[2]], self._softclamp_lam) / V_poly
                        ok = mu1_p > SMALL_EPS
                        if ok and r_mol_s < 0.0:
                            # Net reverse = coupling bookkeeping; it depletes
                            # the DAUGHTER, so guard its moments too. (Stated
                            # approximation: parent statistics, sign-flipped.)
                            if (max(0.0, y[d_idx[0]]) / V_poly <= SMALL_EPS or
                                    max(0.0, y[d_idx[1]]) / V_poly <= SMALL_EPS):
                                ok = False
                        if ok:
                            # Complement-stays-in-parent accounting: parent
                            # mu0 net 0 (chain broke, complement remains);
                            # fragment (uniform cut of a length-biased chain:
                            # E[a] = E[n]/2, E[a^2] = E[n^2]/3) moves to the
                            # daughter. mu1 conserves exactly per reaction.
                            e_n = mu2_p / mu1_p
                            dn_dt[s_idx[1]] -= r_mol_s * e_n / 2.0
                            dn_dt[d_idx[0]] += r_mol_s
                            dn_dt[d_idx[1]] += r_mol_s * e_n / 2.0
                            mu3_p = _safe_mu3_from_mu012(mu0_p, mu1_p, mu2_p)
                            if np.isfinite(mu3_p):
                                e_n2 = mu3_p / mu1_p
                                dn_dt[s_idx[2]] -= r_mol_s * (2.0 / 3.0) * e_n2
                                dn_dt[d_idx[2]] += r_mol_s * e_n2 / 3.0
                elif arch == FLUX_DISCRETE_CHIP:
                    # Discrete chip (spec 2026-06-10 §5): an end-anchored cut
                    # ejects a stamped a-unit chip to the gas phase; the
                    # complement folds back into the SAME pool, so mu0 is
                    # unchanged and there is no dst pool. Closure-free in
                    # effect: only E[n] (an exact ratio of tracked moments)
                    # is consumed -- b2/mu2_ok are ignored here, so no mu2_ok
                    # branch is needed. (For mu1-scaled chips _chain_bundle's
                    # length-biased path still computes the mu3 closure
                    # internally; its output is unused by this branch.) The
                    # bundle and the rate scaling key on the same stored flag;
                    # in the throttled exhaustion regime the rate measure is
                    # donor-limited while the pick stays population-based --
                    # accepted approximation, see spec s5 amendment. The chip
                    # species itself gains/loses moles through the standard
                    # gas dn_dt path (section 4) -- no special handling here.
                    src = self.reaction_src_pool[r_idx]
                    # src == -1 cannot reach here (init demotes); defensive.
                    if src != -1:
                        chip_a = float(self.reaction_chip_units[r_idx])
                        b0, b1, _b2, _mu2_ok = self._chain_bundle(
                            src, y, V_poly, self.is_end_group_reaction[r_idx])
                        if b0 != 0.0:
                            s_idx = self.polymer_pools[src].mu_indices
                            chip_e_n = b1
                            if rf > 0.0:
                                # Forward (ejection): Delta n = -a,
                                # Delta n^2 = -(2aE[n] - a^2). Clamp the mu2
                                # decrement at >= 0: 2aE[n] - a^2 < 0 is
                                # impossible per-chain (n >= a always) but
                                # reachable in expectation when the pool mean
                                # length has decayed toward chip size; by
                                # then the moment description is marginal.
                                rf_mol = rf * V_rxn
                                dn_dt[s_idx[1]] -= rf_mol * chip_a
                                chip_dmu2 = 2.0 * chip_a * chip_e_n - chip_a * chip_a
                                if chip_dmu2 > 0.0:
                                    dn_dt[s_idx[2]] -= rf_mol * chip_dmu2
                            if rr > 0.0:
                                # Reverse (re-attachment) -- EXACT extension
                                # form, NOT the forward sign-flip:
                                # (n+a)^2 - n^2 = +(2aE[n] + a^2). E[n] on the
                                # same single pool (the re-formed chain
                                # extends a parent-statistics chain by a).
                                # Unconditionally positive: no clamp.
                                rr_mol = rr * V_rxn
                                dn_dt[s_idx[1]] += rr_mol * chip_a
                                dn_dt[s_idx[2]] += rr_mol * (
                                    2.0 * chip_a * chip_e_n + chip_a * chip_a)
                elif arch == FLUX_MOMENT_CREDIT_CONDUIT:
                    # Moment-credit conduit (M18.4, DESIGN §2): an admitted
                    # refused-class row consumes a real chain-scale discrete
                    # and CREDITS the destination pool a Dirac point-mass
                    # bundle (1, u, u^2) per event -- credit-ONLY, so there is
                    # NO source-pool debit (src is NULL) and NO extra gas
                    # dn/dt write here (the reactant discrete and the single
                    # gas product are already moved by the standard section-4
                    # net-rate stoichiometry). u = chain_units is precomputed
                    # at admission and flattened into conduit_chain_units_u.
                    # F is the net event flux r_mol_s (mol/s), the SAME
                    # variable the sibling archetype arms consume; irreversible
                    # export (OQ-2) means rr = 0, so r_mol_s = rf*V_rxn = F.
                    # NOTE: the residual deliberately does NOT stash F for the
                    # accepted-step grams booking. The residual is a DASPK
                    # callback re-run on Newton trials and by the FD Jacobian
                    # (_scoped_fd_jacobian, perturbed state), so any scratch it
                    # left would reflect a perturbed/rejected point, not the
                    # accepted solution. _book_conduit_flux_accepted instead
                    # RECOMPUTES F from the accepted self.y off the accepted-
                    # step hook (side-effect-free mirror of this arm's rate).
                    dst = self.reaction_dst_pool[r_idx]
                    # dst == -1 cannot reach here (init refuses dst-less arch-7
                    # rows, flux-dead); the check is defensive only.
                    if dst != -1:
                        F = r_mol_s
                        u = self.conduit_chain_units_u[r_idx]
                        d_idx = self.polymer_pools[dst].mu_indices
                        dn_dt[d_idx[0]] += F
                        dn_dt[d_idx[1]] += F * u
                        dn_dt[d_idx[2]] += F * u * u
                elif arch == FLUX_UNRESOLVED:
                    # Legacy mu1-only transfer (pre-apportionment behavior),
                    # replicated exactly: -r per reactant proxy, +r per
                    # product proxy. NOTE: mu0-scaled shapes here share the
                    # chip exhaustion structure but are deliberately NOT
                    # throttled (bit-exact legacy contract; see
                    # docs/multi_pool_design.md limitation 14 / spec A4).
                    if self.is_pool_proxy[r0]:
                        dn_dt[self.polymer_pools[p0_pool_idx].mu_indices[1]] -= r_mol_s
                    if r1 != -1 and self.is_pool_proxy[r1]:
                        dn_dt[self.polymer_pools[p1_pool_idx].mu_indices[1]] -= r_mol_s
                    for p_slot in range(3):
                        p_idx_tmp = ip[r_idx, p_slot]
                        if (p_idx_tmp != -1 and p_idx_tmp < n_core
                                and self.is_pool_proxy[p_idx_tmp]):
                            pool_idx = self.species_to_pool_indices[p_idx_tmp]
                            dn_dt[self.polymer_pools[pool_idx].mu_indices[1]] += r_mol_s

        # 7. Network Leaks
        for j in range(inet.shape[0]):
            s0 = inet[j, 0]
            if s0 == -1 or s0 >= n_core or not mask[s0]:
                continue

            s1, s2 = inet[j, 1], inet[j, 2]
            if (s1 != -1 and not mask[s1]) or (s2 != -1 and not mask[s2]):
                continue

            k = knet[j]
            rate = k * C_gas[s0]
            if s1 != -1: rate *= C_gas[s1]
            if s2 != -1: rate *= C_gas[s2]

            # Record only: the leak rate drives model enlargement (a network
            # whose leak flux is large gets expanded) and is NOT a sink. This
            # matches simple.pyx, which never debits the core balance; a
            # debit here removed core gas species and broke carbon closure.
            self.network_leak_rates[j] = rate

        # The backbone bonds of each QSSA lineage (a parent and its X-loss
        # daughters, I-260), from the same moments the loop below reads.
        qssa_lineage_B = None
        if getattr(self, "qssa_lineage_members", None):
            qssa_lineage_B = np.zeros(len(self.polymer_pools), dtype=float)
            for p_lin in self.qssa_lineage_members:
                idx_lin = self.polymer_pools[p_lin].mu_indices
                qssa_lineage_B[self.qssa_root[p_lin]] += max(
                    _softclamp(y[idx_lin[1]], self._softclamp_lam) / V_poly
                    - _softclamp(y[idx_lin[0]], self._softclamp_lam) / V_poly,
                    0.0)

        # 8. Polymer Tail & Handshake
        for pool_i, pool in enumerate(self.polymer_pools):
            idx_mu0, idx_mu1, idx_mu2 = pool.mu_indices
            xs = pool.xs

            mu0_mol = _softclamp(y[idx_mu0], self._softclamp_lam)
            mu1_mol = _softclamp(y[idx_mu1], self._softclamp_lam)
            mu2_mol = _softclamp(y[idx_mu2], self._softclamp_lam)

            mu0 = mu0_mol / V_poly
            mu1 = mu1_mol / V_poly
            mu2 = mu2_mol / V_poly
            mu3 = _safe_mu3_from_mu012(mu0, mu1, mu2)

            if self.debug_check_realizability and pool.label not in self._realizability_warned:
                # Use raw (pre-clamp) moles so negatives are caught too. All terms
                # are extensive, so the inequalities are scale-invariant in V_poly.
                raw0, raw1, raw2 = y[idx_mu0], y[idx_mu1], y[idx_mu2]
                if (raw0 < -1e-9 or raw1 < -1e-9 or raw2 < -1e-9
                        or raw1 + 1e-9 < raw0
                        or raw0 * raw2 + 1e-9 < raw1 * raw1):
                    self._realizability_warned.add(pool.label)
                    logging.warning(
                        "Polymer pool '%s' moment state left the realizable cone "
                        "(mu0=%.6g, mu1=%.6g, mu2=%.6g; require mu1>=mu0>=0 and "
                        "mu0*mu2>=mu1^2). The mu3 closure is guarded, but a moment "
                        "source term is likely wrong.", pool.label, raw0, raw1, raw2)

            dmu0_dt = 0.0
            dmu1_dt = 0.0
            dmu2_dt = 0.0
            small_src = dict()

            if pool.tail_kinetics:
                dmu0_dt, dmu1_dt, dmu2_dt, small_src = pool.tail_kinetics(
                    self.T.value_si, self.P.value_si, mu0, mu1, mu2, mu3)
            else:
                if pool.k_scission > 0:
                    # Random backbone scission, discrete-bond (Ziff-McGrady)
                    # convention: a chain of length k has (k-1) breakable bonds,
                    # so the distribution's breakable-bond count is (mu1 - mu0):
                    #   dμ0/dt = k_s·(μ1 − μ0)     (one new chain per bond broken)
                    #   dμ1/dt = 0                 (monomer units / mass conserved)
                    #   dμ2/dt = (k_s/3)·(μ1 − μ3)
                    # See docs/multi_pool_design.md §5. The μ0 term MUST carry the
                    # −μ0 depletion part: with the bare +k_s·μ1 form μ0 grows past
                    # μ1, the state leaves the realizable cone (μ1 ≥ μ0 always for a
                    # k≥1 distribution), and the μ3 closure blows up to a DASSL
                    # singularity. The (μ1 − μ0) form self-limits (rate → 0 as the
                    # pool reaches all-monomer) and structurally keeps μ0 ≤ μ1.
                    dmu0_dt += pool.k_scission * (mu1 - mu0)
                    if np.isfinite(mu3):
                        dmu2_dt += pool.k_scission * (mu1 - mu3) / 3.0

                if pool.k_unzip > 0:
                    # Chain-end monomer release: one active end per chain, so
                    # events occur at k_unzip*mu0 and each removes ONE repeat
                    # unit from its chain.
                    #   dmu1 = -k_u*mu0
                    #   dmu2 = -k_u*(2*mu1 - mu0)   (a DP=n chain loses 2n-1)
                    #   dmu0 = -k_u*N1              (see below)
                    # The dmu0 term is what the legacy form never carried. At
                    # the all-monomer boundary mu1 = mu0 every chain has
                    # length 1, and unzipping its last unit removes a CHAIN,
                    # not just a unit. Without that debit the drain -k_u*mu0
                    # stays at full strength as mu1 falls to mu0, so mu1 is
                    # pushed THROUGH mu0, out of the realizable cone
                    # (mu1 >= mu0 >= 0 for any k>=1 distribution) and then
                    # through zero -- the accepted-state r81 raise, not a
                    # trial-state excursion. With it the boundary is an
                    # INVARIANT set of the vector field:
                    #   d(mu1 - mu0)/dt = -k_u*mu0*(1 - p1) -> 0 as p1 -> 1,
                    # which is the same self-limiting structure the scission
                    # kernel above gets from its (mu1 - mu0) factor.
                    # N1 = mu0*p1 comes from the SAME closure the
                    # k_depropagation sibling uses (_deprop_dp1_fraction:
                    # gamma leg plus a smooth terminal floor with p1 == 1 for
                    # mean DP <= 1). That floor is also what keeps this a
                    # TERMINATION rather than a stall: the last repeat unit
                    # per chain still drains, at dmu0 = -k_u*mu0, instead of
                    # freezing the residue (r74).
                    # Deliberately NOT gated: the unit drain, the mu2 drain
                    # and the released-monomer flux are unchanged, bit for
                    # bit, so this fixes realizability without touching how
                    # much mass the unzip channel removes.
                    #
                    # The debit is written as (release rate)*N1 because that
                    # is the CHANNEL-INDEPENDENT form of the law: releases
                    # occur one per active chain end, so a fraction N1 of
                    # them land on a DP=1 chain and remove it. The QSSA
                    # channel below carries the SAME term against its own
                    # release rate r_qssa -- the two are the same chain-end
                    # event, and the solver PINS them equivalent, so the
                    # debit must be applied to both or to neither.
                    # NOT gated on the explicit-tail handshake. The handshake
                    # removes chains at DP=xs and hands them to an explicit
                    # species; this removes chains at DP=1 that unzipped
                    # their last unit -- different populations, both real. A
                    # gate on the handshake's map entry was tried and is
                    # WRONG: it makes the debit differ between two configs
                    # that differ ONLY in whether the tail is explicit, which
                    # breaks the auto-gen path's ON/OFF differencing pin.
                    #
                    # I-065 defect 1: the EVENT rate is set by the chain
                    # count, so off the cone (mu1 = 0, mu0 > 0) it drained
                    # mu1 further negative and emitted k_unzip*mu0 of gas
                    # monomer with no repeat unit behind it -- the solver
                    # reporting mass it never had. _release_units_gate is
                    # EXACTLY 1.0 for mu1 >= mu0, so on every state inside
                    # the realizable cone this is a bit-for-bit no-op
                    # (x*1.0 == x), and it is exactly 0.0 at mu1 <= 0. It
                    # scales the RELEASE -- the mu1 drain, the mu2 drain and
                    # the gas emission, which must move together or the
                    # ledger fabricates/destroys mass -- and deliberately
                    # NOT the I-055 chain-termination debit on mu0 below,
                    # which keeps draining at full strength so a mean-DP < 1
                    # state heals back toward the cone instead of stalling.
                    r_events = pool.k_unzip * mu0
                    r_release = r_events * _release_units_gate(mu0, mu1)
                    dmu1_dt -= r_release
                    dmu2_dt -= (pool.k_unzip * (2.0 * mu1 - mu0)
                                * _release_units_gate(mu0, mu1))
                    dmu0_dt -= r_events * _deprop_dp1_fraction(mu0, mu1, mu2)
                    if pool.monomer_poly_index is not None:
                        # Released monomer is emitted to the GAS species
                        # amount basis (incident 2026-07-03, design B-prime):
                        # small_src [mol/(m^3_poly s)] lands as
                        # dn_dt += r*V_poly [mol/s] on the gas-masked
                        # monomer_poly_index. Mass conservation: one gas
                        # monomer mole per drained mu1 repeat unit.
                        small_src[pool.monomer_poly_index] = r_release
                        self._split_side_loss_release(
                            pool_i, pool.monomer_poly_index, r_release,
                            mu1_mol, V_poly, y, small_src, dz_dt)

            # Radical-homolysis initiation kernel (Stage 1, adjudicated round
            # 66). Independent of tail_kinetics (a custom tail closure does
            # not describe chain initiation); mutually exclusive with
            # k_scission > 0 and with the QSSA channel (validate_configuration
            # hard-errors on both shapes). Reads ONLY the flattened khom_*
            # arrays from _flatten_homolysis_state -- never the pool dict.
            #
            # Event: random backbone C-C homolysis at the RUNTIME temperature
            # (round 66: 'a scalar precomputed at 1100 K will poison any
            # ramp/TA replay'):
            #   k(T) = A * T^n * exp(-Ea/(R_gas*T))
            #   R    = k(T) * max(mu1 - mu0, 0)   [a chain of length n has
            #                                      n-1 breakable bonds]
            # Bond-weighted breaking-chain moments (E[n] = B1, E[n^2] = B2):
            #   B1 = (mu2 - mu1)/(mu1 - mu0);  B2 = (mu3 - mu2)/(mu1 - mu0)
            # with mu3 from the established log-Lagrange closure (computed
            # once per pool above). Parent debit:
            #   dmu0 -= R;  dmu1 -= R*B1;  dmu2 -= R*B2
            # EACH end-radical daughter pool (uniform cut at bond j of an
            # n-chain: E[j] = n/2, E[j^2] = n(2n-1)/6):
            #   dmu0 += R;  dmu1 += R*B1/2;  dmu2 += R*(2*B2 - B1)/6
            # STABLE FORMS (adjudicated round 67, P2): dividing by
            # (mu1 - mu0) and re-multiplying by R cancels catastrophically
            # near DP -> 1 exhaustion (mu1 ~= mu0), so every product-side
            # term is computed in its algebraically identical direct form:
            #   R*B1         = k*(mu2 - mu1)
            #   R*B2         = k*(mu3 - mu2)
            #   R*(2B2-B1)/6 = k*(2*mu3 - 3*mu2 + mu1)/6
            # (R = k*max(mu1-mu0, 0) itself stays the mu0 event rate).
            # Totals reproduce the legacy Ziff-McGrady random-scission moment
            # sources EXACTLY: dmu0_tot = +R = k*(mu1-mu0); dmu1_tot = 0
            # (mass conserved, machine precision); dmu2_tot = -R*(B1+B2)/3
            # = k*(mu1-mu3)/3. NO gas species term (homolysis releases no
            # volatiles) and NO reverse leg (one-way; round 66
            # §Reversibility -- recombination arrives via the discovered
            # chemistry conduit, not this kernel).
            if self.khom_enabled[pool_i]:
                B_hom = mu1 - mu0
                if B_hom > 0.0 and mu0 > 0.0:
                    T_hom = self.T.value_si
                    k_hom = (self.khom_A[pool_i]
                             * T_hom ** self.khom_n[pool_i]
                             * math.exp(-self.khom_Ea[pool_i]
                                        / (QSSA_R_GAS * T_hom)))
                    if not (0.0 < k_hom < QSSA_INF):
                        raise ValueError(
                            f"Pool {pool.label}: k_homolysis(T={T_hom:g} K) "
                            f"evaluated to a degenerate rate constant "
                            f"({k_hom!r}) -- A={self.khom_A[pool_i]:g}, "
                            f"n={self.khom_n[pool_i]:g}, "
                            f"Ea={self.khom_Ea[pool_i]:g} J/mol. Refusing "
                            f"to integrate a poisoned kernel.")
                    mu3_ok = np.isfinite(mu3)
                    # Direct (stable) moment differences: d1 = R*B1/k,
                    # d2 = R*B2/k, c2 = 6*(daughter mu2 credit)/k.
                    d1_hom = mu2 - mu1
                    d2_hom = (mu3 - mu2) if mu3_ok else -1.0
                    c2_hom = (2.0 * mu3 - 3.0 * mu2 + mu1) if mu3_ok else -1.0
                    if ((not mu3_ok) or d1_hom < 0.0 or d2_hom < 0.0
                            or c2_hom < 0.0):
                        # Loud guard (warn once per pool per rebuild,
                        # mirroring the realizability idiom above): the
                        # moment state is outside the kernel's domain
                        # (mu2 < mu1, a degenerate/overflowed mu3 closure,
                        # or -- round 67 P1 -- a nonrealizable combination
                        # where B1, B2 >= 0 but the daughter mu2 credit
                        # factor 2*B2 - B1 = c2/(mu1-mu0) is NEGATIVE, which
                        # would credit negative mu2 to both end-radical
                        # daughters). Contribute ZERO flux this step rather
                        # than integrate garbage fragment moments.
                        if pool.label not in self._khom_guard_warned:
                            self._khom_guard_warned.add(pool.label)
                            logging.warning(
                                "Polymer pool '%s': k_homolysis kernel "
                                "skipped -- nonfinite/negative moment "
                                "combination (mu0=%.6g, mu1=%.6g, mu2=%.6g, "
                                "mu3=%.6g; mu2-mu1=%.6g, mu3-mu2=%.6g, "
                                "daughter-mu2 factor 2*mu3-3*mu2+mu1=%.6g). "
                                "The kernel contributes zero flux until the "
                                "state re-enters its realizable domain.",
                                pool.label, mu0, mu1, mu2, mu3,
                                d1_hom, d2_hom, c2_hom)
                    else:
                        r_hom = k_hom * B_hom
                        dmu0_dt -= r_hom
                        dmu1_dt -= k_hom * d1_hom
                        dmu2_dt -= k_hom * d2_hom
                        # Fragment credits go straight to the daughter
                        # pools' moment slots (they are core species; the
                        # daughter's own accumulator flush is additive).
                        frag_mu0 = r_hom * V_poly
                        frag_mu1 = (k_hom * d1_hom / 2.0) * V_poly
                        frag_mu2 = (k_hom * c2_hom / 6.0) * V_poly
                        dn_dt[self.khom_prim_mu0[pool_i]] += frag_mu0
                        dn_dt[self.khom_prim_mu1[pool_i]] += frag_mu1
                        dn_dt[self.khom_prim_mu2[pool_i]] += frag_mu2
                        dn_dt[self.khom_sec_mu0[pool_i]] += frag_mu0
                        dn_dt[self.khom_sec_mu1[pool_i]] += frag_mu1
                        dn_dt[self.khom_sec_mu2[pool_i]] += frag_mu2

            # Catalog-driven kernels (S5 a + b), per template, in the term's
            # period: m_j = c^j mu_j (c = periods_per_unit); every move is
            # converted back by c^-j. Read only the flattened term lists.
            #
            # kernel_homolysis: a chain of n periods has a*n + b sites, a cut
            # at k in [e, n - e] (e = min_fragment) leaves k and n - k:
            #   events  r0 = w*phi*(a*m1 + b*m0)
            #   parent  dm_j -= w*phi*(a*m_{j+1} + b*m_j)        (j = 0, 1, 2)
            #   daughter dm0 += 2*r0;  dm1 += w*phi*(a*m2 + b*m1)
            #           dm2 += w*phi*2a*sum_n N_n sum_{k=e}^{n-e} k^2
            #              = w*phi*(a/3)*(2*m3 + (3 - 6e)*m2
            #                 + (6e^2 - 6e + 1)*m1 + (-4e^3 + 6e^2 - 2e)*m0)
            # exact in moments but for m3, the log-Lagrange closure (exact on
            # a monodisperse pool). Periods are conserved to rounding.
            if mu0 > 0.0 and mu1 > 0.0:
                for term, (d0, d1, d2) in self.kernel_homolysis_terms[pool_i]:
                    wphi = term["w"] * term["phi"]
                    if wphi == 0.0:
                        continue
                    a_k, b_k = term["a"], term["b"]
                    e_k = float(term["min_fragment"])
                    c_k = float(term["periods_per_unit"])
                    m1_k, m2_k = c_k * mu1, c_k * c_k * mu2
                    m3_k = c_k * c_k * c_k * mu3
                    r0 = wphi * (a_k * m1_k + b_k * mu0)
                    r1 = wphi * (a_k * m2_k + b_k * m1_k)
                    r2 = wphi * (a_k * m3_k + b_k * m2_k)
                    q2 = wphi * (a_k / 3.0) * (
                        2.0 * m3_k + (3.0 - 6.0 * e_k) * m2_k
                        + (6.0 * e_k * e_k - 6.0 * e_k + 1.0) * m1_k
                        + (-4.0 * e_k ** 3 + 6.0 * e_k * e_k - 2.0 * e_k) * mu0)
                    if not (np.isfinite(m3_k) and r0 >= 0.0 and r1 >= 0.0
                            and r2 >= 0.0 and q2 >= 0.0):
                        # outside the law's domain (chains shorter than the
                        # sites allow, or a degenerate closure): no flux,
                        # loudly, rather than negative chains
                        if pool.label not in self._kernel_guard_warned:
                            self._kernel_guard_warned.add(pool.label)
                            logging.warning(
                                "Polymer pool '%s': kernel_homolysis term %s "
                                "skipped -- moments (%.6g, %.6g, %.6g, mu3 "
                                "%.6g) give rates (%.6g, %.6g, %.6g, %.6g) "
                                "outside the cut law's domain.", pool.label,
                                term["template"], mu0, mu1, mu2, mu3, r0, r1,
                                r2, q2)
                        continue
                    dmu0_dt -= r0
                    dmu1_dt -= r1 / c_k
                    dmu2_dt -= r2 / (c_k * c_k)
                    dn_dt[d0] += 2.0 * r0 * V_poly
                    dn_dt[d1] += (r1 / c_k) * V_poly
                    dn_dt[d2] += (q2 / (c_k * c_k)) * V_poly

            # kernel_homolysis, two daughters (I-278): a chain of n periods
            # has a*n + b sites, b = a*(1 - r - e - e'); a cut at k leaves k
            # periods in the daughter (k from e to n - r - e') and n - r - k
            # in the partner daughter (their range is e' to n - r - e):
            #   events  r0 = w*phi*(a*m1 + b*m0);  parent as above
            #   daughter dm0 += r0; dm_j += w*phi*a*sum_n N_n S_j(n), j = 1, 2
            # with S_j(n) the exact polynomial sum_{k=lo}^{n-shift} k^j
            # (_range_power_sums), so dm1 reads m0..m2 and dm2 m0..m3 (the
            # same closure for m3).
            if mu0 > 0.0 and mu1 > 0.0:
                for term, (d0, d1, d2), (q0, q1, q2_), sums_d, sums_q in \
                        self.kernel_homolysis_two_terms[pool_i]:
                    wphi = term["w"] * term["phi"]
                    if wphi == 0.0:
                        continue
                    a_k, b_k = term["a"], term["b"]
                    c_k = float(term["periods_per_unit"])
                    m_k = (mu0, c_k * mu1, c_k * c_k * mu2,
                           c_k * c_k * c_k * mu3)
                    r0 = wphi * (a_k * m_k[1] + b_k * m_k[0])
                    r1 = wphi * (a_k * m_k[2] + b_k * m_k[1])
                    r2 = wphi * (a_k * m_k[3] + b_k * m_k[2])
                    credits = []
                    for (s1, s2) in (sums_d, sums_q):
                        credits.append((
                            wphi * a_k * sum(c * m_k[i] for i, c in enumerate(s1)),
                            wphi * a_k * sum(c * m_k[i] for i, c in enumerate(s2))))
                    if not (np.isfinite(m_k[3]) and r0 >= 0.0 and r1 >= 0.0
                            and r2 >= 0.0
                            and all(np.isfinite(v) and v >= 0.0
                                    for credit in credits for v in credit)):
                        if pool.label not in self._kernel_guard_warned:
                            self._kernel_guard_warned.add(pool.label)
                            logging.warning(
                                "Polymer pool '%s': kernel_homolysis term %s "
                                "skipped -- moments (%.6g, %.6g, %.6g, mu3 "
                                "%.6g) give rates (%.6g, %.6g, %.6g) and "
                                "credits %s outside the cut law's domain.",
                                pool.label, term["template"], mu0, mu1, mu2,
                                mu3, r0, r1, r2, credits)
                        continue
                    dmu0_dt -= r0
                    dmu1_dt -= r1 / c_k
                    dmu2_dt -= r2 / (c_k * c_k)
                    for (cr1, cr2), (s0, s1_, s2_) in zip(
                            credits, ((d0, d1, d2), (q0, q1, q2_))):
                        dn_dt[s0] += r0 * V_poly
                        dn_dt[s1_] += (cr1 / c_k) * V_poly
                        dn_dt[s2_] += (cr2 / (c_k * c_k)) * V_poly

            # kernel_recombination: two chains of this pool, s sites each,
            # join into one chain of n0 + n1 periods in the product pool at
            # k_pair*(s*m0)^2 per volume (coagulation, closed: no m3):
            #   this    dm0 -= 2*r;  dm1 -= 2*k s^2 m0 m1;  dm2 -= 2*k s^2 m0 m2
            #   product dm0 += r;    dm1 += 2*k s^2 m0 m1
            #           dm2 += k s^2 (2*m0*m2 + 2*m1^2)
            # (mu2 gains 2*(m1/m0)^2 per event beyond the transfer).
            # I-266, a join of two pools (a term with a partner): one chain of
            # this pool (s sites) and one of the partner (s' sites, moments
            # nu_j) join into one chain of n + m in the product pool at
            # K = k_pair*s*s' per pair (a convolution, closed: no m3; the
            # same in any length unit):
            #   this    dmu_j -= K mu_j nu0;   partner dnu_j -= K mu0 nu_j
            #   product dmu0 += K mu0 nu0;     dmu1 += K (mu1 nu0 + mu0 nu1)
            #           dmu2 += K (mu2 nu0 + 2 mu1 nu1 + mu0 nu2)
            if mu0 > 0.0:
                for term, (p0, p1, p2), partner in self.kernel_recombination_terms[pool_i]:
                    if partner is not None:
                        q0, q1, q2 = partner
                        nu0 = _softclamp(y[q0], self._softclamp_lam) / V_poly
                        if not nu0 > 0.0:
                            continue
                        nu1 = _softclamp(y[q1], self._softclamp_lam) / V_poly
                        nu2 = _softclamp(y[q2], self._softclamp_lam) / V_poly
                        kss = (term["k_pair"] * term["sites_per_chain"]
                               * term["partner_sites_per_chain"])
                        dmu0_dt -= kss * mu0 * nu0
                        dmu1_dt -= kss * mu1 * nu0
                        dmu2_dt -= kss * mu2 * nu0
                        dn_dt[q0] -= kss * mu0 * nu0 * V_poly
                        dn_dt[q1] -= kss * mu0 * nu1 * V_poly
                        dn_dt[q2] -= kss * mu0 * nu2 * V_poly
                        dn_dt[p0] += kss * mu0 * nu0 * V_poly
                        dn_dt[p1] += kss * (mu1 * nu0 + mu0 * nu1) * V_poly
                        dn_dt[p2] += kss * (mu2 * nu0 + 2.0 * mu1 * nu1
                                            + mu0 * nu2) * V_poly
                        continue
                    c_k = float(term["periods_per_unit"])
                    ks2 = (term["k_pair"] * term["sites_per_chain"]
                           * term["sites_per_chain"])
                    m1_k, m2_k = c_k * mu1, c_k * c_k * mu2
                    r_ev = ks2 * mu0 * mu0
                    g1 = 2.0 * ks2 * mu0 * m1_k
                    dmu0_dt -= 2.0 * r_ev
                    dmu1_dt -= g1 / c_k
                    dmu2_dt -= 2.0 * ks2 * mu0 * m2_k / (c_k * c_k)
                    dn_dt[p0] += r_ev * V_poly
                    dn_dt[p1] += (g1 / c_k) * V_poly
                    dn_dt[p2] += (ks2 * (2.0 * mu0 * m2_k + 2.0 * m1_k * m1_k)
                                  / (c_k * c_k)) * V_poly

            # kernel_transfer (S5 c + e): one chain of this pool and one of
            # the partner pool change pool, lengths unchanged, at
            # k*phi*c*c_partner per volume (no closure, any length unit):
            #   this    dmu_j -= k phi mu_j mu0'   -> product
            #   partner dmu_j' -= k phi mu0 mu_j'  -> partner_product
            # (partner == this pool: both debits land on it, ordered pairs).
            # S5 d, a term with `sites` (c00, c10, c01, c11) in its period
            # (m_j = c^j mu_j here, n_j = c^j nu_j in the partner): a pair of
            # n and m periods reacts at k phi (c00 + c10 n + c01 m + c11 n m),
            # so each chain is picked in proportion to its sites:
            #   this    dm_j -= k phi (c00 m_j n0 + c10 m_{j+1} n0
            #                          + c01 m_j n1 + c11 m_{j+1} n1)
            #   partner dn_j -= k phi (c00 m0 n_j + c10 m1 n_j
            #                          + c01 m0 n_{j+1} + c11 m1 n_{j+1})
            # exact but for m3 / n3 on a side whose sites grow (the
            # log-Lagrange closure); a flow outside the law's domain (a
            # negative or unbounded moment) is skipped, loudly.
            if mu0 > 0.0:
                for term, (q0, q1, q2), (a0, a1, a2), (b0, b1, b2) in \
                        self.kernel_transfer_terms[pool_i]:
                    kphi = term["k_pair"] * term["phi"]
                    if kphi == 0.0:
                        continue
                    nu0 = _softclamp(y[q0], self._softclamp_lam) / V_poly
                    if not nu0 > 0.0:
                        continue
                    nu1 = _softclamp(y[q1], self._softclamp_lam) / V_poly
                    nu2 = _softclamp(y[q2], self._softclamp_lam) / V_poly
                    sites = term.get("sites")
                    if sites is None:
                        own = (kphi * mu0 * nu0 * V_poly,
                               kphi * mu1 * nu0 * V_poly,
                               kphi * mu2 * nu0 * V_poly)
                        other = (kphi * mu0 * nu0 * V_poly,
                                 kphi * mu0 * nu1 * V_poly,
                                 kphi * mu0 * nu2 * V_poly)
                    else:
                        c00, c10, c01, c11 = sites
                        c_k = float(term["periods_per_unit"])
                        m = (mu0, c_k * mu1, c_k * c_k * mu2,
                             c_k ** 3 * mu3)
                        n = (nu0, c_k * nu1, c_k * c_k * nu2,
                             c_k ** 3 * _safe_mu3_from_mu012(nu0, nu1, nu2))
                        own = tuple([kphi * (c00 * m[j] * n[0]
                                             + c10 * m[j + 1] * n[0]
                                             + c01 * m[j] * n[1]
                                             + c11 * m[j + 1] * n[1])
                                     / c_k ** j * V_poly for j in range(3)])
                        other = tuple([kphi * (c00 * m[0] * n[j]
                                               + c10 * m[1] * n[j]
                                               + c01 * m[0] * n[j + 1]
                                               + c11 * m[1] * n[j + 1])
                                       / c_k ** j * V_poly for j in range(3)])
                        if not all([np.isfinite(v) and v >= 0.0
                                    for v in own + other]):
                            if pool.label not in self._kernel_guard_warned:
                                self._kernel_guard_warned.add(pool.label)
                                logging.warning(
                                    "Polymer pool '%s': kernel_transfer term "
                                    "%s skipped -- flows %s / %s outside its "
                                    "site law's domain.", pool.label,
                                    term["template"], own, other)
                            continue
                    dmu0_dt -= own[0] / V_poly
                    dmu1_dt -= own[1] / V_poly
                    dmu2_dt -= own[2] / V_poly
                    dn_dt[a0] += own[0]
                    dn_dt[a1] += own[1]
                    dn_dt[a2] += own[2]
                    dn_dt[q0] -= other[0]
                    dn_dt[q1] -= other[1]
                    dn_dt[q2] -= other[2]
                    dn_dt[b0] += other[0]
                    dn_dt[b1] += other[1]
                    dn_dt[b2] += other[2]

            # kernel_move (S5 d): a chain of n periods leaves for the product
            # pool, length unchanged, at k*(a n + b) (m_j = c^j mu_j):
            #   this dm_j -= k (a m_{j+1} + b m_j) -> product
            # exact but for m3 when a != 0 (the log-Lagrange closure); a flow
            # outside the law's domain is skipped, loudly.
            if mu0 > 0.0:
                for term, (p0, p1, p2) in self.kernel_move_terms[pool_i]:
                    c_k = float(term["periods_per_unit"])
                    m = (mu0, c_k * mu1, c_k * c_k * mu2, c_k ** 3 * mu3)
                    flow = tuple([term["k"] * (term["a"] * m[j + 1]
                                               + term["b"] * m[j])
                                  / c_k ** j for j in range(3)])
                    if not all([np.isfinite(v) and v >= 0.0 for v in flow]):
                        if pool.label not in self._kernel_guard_warned:
                            self._kernel_guard_warned.add(pool.label)
                            logging.warning(
                                "Polymer pool '%s': kernel_move term %s "
                                "skipped -- flows %s outside its site "
                                "law's domain.", pool.label,
                                term["template"], flow)
                        continue
                    dmu0_dt -= flow[0]
                    dmu1_dt -= flow[1]
                    dmu2_dt -= flow[2]
                    dn_dt[p0] += flow[0] * V_poly
                    dn_dt[p1] += flow[1] * V_poly
                    dn_dt[p2] += flow[2] * V_poly

            # kernel_propagation (I-266): a chain of n units adds one monomer
            # unit at k_pair*c_M per chain, c_M = n_M / V_poly (kMC's
            # convention: every small species at its melt concentration),
            # the monomer consumed (small_src -> dn_dt*V_poly):
            #   dmu0 = 0;  dmu1 = k c_M mu0;  dmu2 = k c_M (2 mu1 + mu0)
            if mu0 > 0.0 and self.kernel_propagation_terms[pool_i]:
                m_idx = pool.monomer_poly_index
                c_mon = _softclamp(y[m_idx], self._softclamp_lam) / V_poly
                if c_mon > 0.0:
                    for term in self.kernel_propagation_terms[pool_i]:
                        kc = term["k_pair"] * c_mon
                        dmu1_dt += kc * mu0
                        dmu2_dt += kc * (2.0 * mu1 + mu0)
                        small_src[m_idx] = small_src.get(m_idx, 0.0) - kc * mu0

            # kernel_interior_scission (I-273): every chain of this pool (one
            # position class of mid-chain radicals) breaks at k per chain;
            # each product pool receives fixed coefficients of this pool's
            # moments, in the term's period (m_j = c^j mu_j; no closure):
            #   this    dmu_j -= k mu_j
            #   product dm_j  += k sum_{i<=j} q_ji m_i
            # q_ji the coefficients of E[sum of its fragment lengths^j | n].
            # A credit outside the law's domain (a negative or unbounded
            # moment) skips the term, loudly.
            if mu0 > 0.0:
                for term, products in self.kernel_interior_scission_terms[pool_i]:
                    k_s = term["k"]
                    c_k = float(term["periods_per_unit"])
                    m = (mu0, c_k * mu1, c_k * c_k * mu2)
                    credits = [[k_s * math.fsum([q[j][i] * m[i] for i in range(j + 1)])
                                / c_k ** j for j in range(3)] for q, _slots in products]
                    if not all([np.isfinite(v) and v >= 0.0
                                for credit in credits for v in credit]):
                        if pool.label not in self._kernel_guard_warned:
                            self._kernel_guard_warned.add(pool.label)
                            logging.warning(
                                "Polymer pool '%s': kernel_interior_scission "
                                "term %s skipped -- credits %s outside its "
                                "fragment law's domain.", pool.label,
                                term["template"], credits)
                        continue
                    dmu0_dt -= k_s * mu0
                    dmu1_dt -= k_s * mu1
                    dmu2_dt -= k_s * mu2
                    for credit, (_q, (p0, p1, p2)) in zip(credits, products):
                        dn_dt[p0] += credit[0] * V_poly
                        dn_dt[p1] += credit[1] * V_poly
                        dn_dt[p2] += credit[2] * V_poly

            # Side-group homolysis initiation kernel (SGH kernel-v2, explicit
            # Br-inventory depletion; adjudicated spar r103 / plan r104).
            # Reads ONLY the flattened sgh_* row arrays from
            # _flatten_side_group_state -- never the pool dicts. Kernel-v2 is
            # mutually exclusive with k_homolysis (both move/route the carrier
            # moments): the §1h build guard refuses any carrier whose mu1
            # moves, including the k_homolysis kernel above. It MAY coexist
            # only with k_scission (verified mu1-preserving backbone
            # scission).
            #
            # v2 carries an EXPLICIT auxiliary extensive state Z per (pool,
            # channel row): Z = moles of remaining removable X (Br) sites,
            # seeded Z(0) = sites*mu1(0) [set_initial_conditions step 8] and
            # depleting as the X atoms leave the melt. Per row at the RUNTIME
            # T:
            #   k(T)   = A*T^n*exp(-Ea/(R_gas*T))   [s^-1 per site]
            #   z_conc = Z / V_poly                 [mol/m^3]
            #   R      = k(T) * z_conc              [mol/(m^3 s)]
            # There is NO extra `sites` factor in R -- `sites` was FOLDED INTO
            # Z(0)=sites*mu1(0), so at t0 this reproduces v1's k*sites*mu1
            # rate exactly (Z(0)=sites*mu1). Any residual `*sites` here would
            # double-count (§6.1). Then:
            #   dZ/dt -= R*V_poly                    [mol/s]  (inventory drains)
            #   gas X += R                            through the SAME
            #            small_src -> dn_dt*V_poly path as the unzip/QSSA
            #            releases (one X radical per event; additive across
            #            channels sharing a gas species).
            # The carrier's mu0/mu1/mu2 are UNCHANGED (no moment debit, no
            # feature-pool credit): Br mass leaves through Z, and
            # get_total_polymer_condensed_mass_g books it as lost*M_X so that
            # d(condensed mass)/dt + d(gas X mass)/dt = 0 exactly. The RHS
            # NEVER raises for Z (only sanitizes trial fuzz with max(0.0,.));
            # the accepted-state guard owns the hard [0, sites*mu1] bounds.
            if self.sgh_enabled[pool_i]:
                T_sgh = self.T.value_si
                for sgh_row in range(self.sgh_row_start[pool_i],
                                     self.sgh_row_end[pool_i]):
                    k_sgh = (self.sgh_A[sgh_row]
                             * T_sgh ** self.sgh_n[sgh_row]
                             * math.exp(-self.sgh_Ea[sgh_row]
                                        / (QSSA_R_GAS * T_sgh)))
                    if not (0.0 < k_sgh < QSSA_INF):
                        # Degenerate k is SETUP-invariant (A/n/Ea + T), not a
                        # trial-state artifact, so it stays a hard raise (v1
                        # contract preserved).
                        raise ValueError(
                            f"Pool {pool.label}: side_group_homolysis "
                            f"channel "
                            f"'{self.sgh_channel_labels[sgh_row]}' "
                            f"k(T={T_sgh:g} K) evaluated to a degenerate "
                            f"rate constant ({k_sgh!r}) -- "
                            f"A={self.sgh_A[sgh_row]:g}, "
                            f"n={self.sgh_n[sgh_row]:g}, "
                            f"Ea={self.sgh_Ea[sgh_row]:g} J/mol. Refusing "
                            f"to integrate a poisoned kernel.")
                    if self.sgh_z_slot[sgh_row] < 0:
                        # X-loss feature daughter (I-257; the v1 law): every
                        # repeat unit carries s sites, so a chain of length n
                        # reacts at k*s*n and moves, intact, to the daughter:
                        #   parent  dmu_j -= k*s*mu_{j+1}   (j = 0, 1, 2; mu3
                        #           from the closure above, direct products)
                        #   daughter gains exactly that, gas X gains k*s*mu1;
                        # the daughter's defect M_X books the mass that left.
                        # Zero flux outside the realizable domain (mu2 < mu1,
                        # mu3 nonfinite or < mu2) and on an empty carrier.
                        if (mu0 <= 0.0 or mu1 <= 0.0 or mu2 < mu1
                                or not np.isfinite(mu3) or mu3 < mu2):
                            continue
                        ks_sgh = k_sgh * self.sgh_sites[sgh_row]
                        r0_sgh = ks_sgh * mu1
                        r1_sgh = ks_sgh * mu2
                        r2_sgh = ks_sgh * mu3
                        dmu0_dt -= r0_sgh
                        dmu1_dt -= r1_sgh
                        dmu2_dt -= r2_sgh
                        dn_dt[self.sgh_dst_mu0[sgh_row]] += r0_sgh * V_poly
                        dn_dt[self.sgh_dst_mu1[sgh_row]] += r1_sgh * V_poly
                        dn_dt[self.sgh_dst_mu2[sgh_row]] += r2_sgh * V_poly
                        g_idx_sgh = int(self.sgh_gas[sgh_row])
                        small_src[g_idx_sgh] = (
                            small_src.get(g_idx_sgh, 0.0) + r0_sgh)
                        # A moved chain of length n brings s*n sites, one of
                        # them lost (the daughter's per-chain defect): its
                        # intact-site inventory gains s*k*s*mu2 - k*s*mu1
                        # (I-260).
                        if self.sgh_dst_z[sgh_row] >= 0:
                            dz_dt[self.sgh_dst_z[sgh_row] - n_core
                                  - self.num_qssa_u] += (
                                (self.sgh_sites[sgh_row] * r1_sgh - r0_sgh)
                                * V_poly)
                        continue
                    # Trial-state fuzz sanitize ONLY (never raise for Z here).
                    z_amount = y[int(self.sgh_z_slot[sgh_row])]
                    if z_amount < 0.0:
                        z_amount = 0.0
                    r_sgh = k_sgh * z_amount / V_poly     # mol/(m^3 s)
                    # Deplete the Br inventory. z_rel is this row's index in
                    # the [Z] block (state order core | U | Z):
                    # z_slot - n_core - num_qssa_u.
                    z_rel = (int(self.sgh_z_slot[sgh_row])
                             - n_core - self.num_qssa_u)
                    dz_dt[z_rel] -= r_sgh * V_poly
                    # Gas X radical: +R (same small_src path as v1 ~5167).
                    g_idx_sgh = int(self.sgh_gas[sgh_row])
                    small_src[g_idx_sgh] = (
                        small_src.get(g_idx_sgh, 0.0) + r_sgh)

            # radical_qssa_unzip channel (M2 rate law). Reads ONLY the
            # flattened solver-owned qssa_* arrays from M1 -- NEVER the pool
            # dict -- and gates on qssa_enabled, never on A != 0. Mutually
            # exclusive with k_unzip > 0 (M1 hard error), independent of
            # tail_kinetics (a custom tail closure does not describe
            # radical depropagation).
            #
            # Rate law (isothermal QSSA on the active chain-END radical R):
            #   B    = max(mu1 - mu0, 0)         breakable backbone bonds
            #                                    (mol/m^3 of condensed volume:
            #                                    a chain of length k has k-1
            #                                    bonds, so sum = mu1 - mu0)
            #   G_R  = 2 f ki(T) B               radical generation (homolysis
            #                                    makes TWO end radicals/bond)
            #   no transfer:  0 = G_R - 2 kt R^2          -> R = sqrt(f ki B / kt)
            #   transfer:     0 = G_R - ktr R - 2 kt R^2  -> R = (sqrt(ktr^2
            #                     + 8 kt G_R) - ktr) / (4 kt)
            #     (transfer = first-order active-END loss: H-abstraction turns
            #     an unzipping end radical into a mid-chain radical)
            #   r_mono = y_m kdp(T) R            monomer release, mol/(m^3 s)
            # Arrhenius k(T) = A T^n exp(-Ea/(R_gas T)), Ea in J/mol (M1 SI
            # pin). Units: ki/kdp in s^-1; kt in m^3/(mol s) (bimolecular);
            # ktr in s^-1 -- PSEUDO-FIRST-ORDER, because ktr multiplies R
            # DIRECTLY in the balance above (a literature bimolecular k_tr
            # [L/(mol s) or m^3/(mol s)] must be premultiplied by the
            # relevant substrate concentration [mol/m^3, SI] BEFORE entering
            # the config). The mol/m^3 concentration basis established at
            # C_poly = y/V_poly above means NO unit conversion is applied
            # here.
            r_qssa = 0.0
            if self.qssa_enabled[pool_i]:
                B_qssa = max(mu1 - mu0, 0.0)
                if self.qssa_weaklink[pool_i]:
                    # Weak-link allyl/U-state channel (milestones ii+iii;
                    # rate law per review r36). Extensions over the legacy
                    # algebra (everything else identical -- the bridge pin
                    # proves it):
                    #   kt_total = kt_rec(T) + kt_disp(T)   replaces kt
                    #   u_active = min(max(U,0)/V_poly, B)  [mol/m^3]
                    #     ACTIVE-SITE CLAMP (r36 P1-2): the weak channel
                    #     fissions a backbone bond ALLYLIC to the
                    #     unsaturated end; ends beyond the remaining
                    #     backbone-bond count have nothing to fission. At
                    #     B = 0 the whole channel is INERT (gate below) --
                    #     no monomer drain that could push mu1 < mu0 past
                    #     the DP->1 self-termination floor, no U sink, no
                    #     production.
                    #   G_R = 2 f ki(T) B + 1 f ki_allyl(T) u_active
                    #     (nu = 1: ONE unzipping radical per weak-link
                    #     fission; the allylic co-fragment does not unzip)
                    #   dU/dt = + kt_disp R_ss^2 * max(0, 1 - U/(2*mu0))
                    #             * V_poly
                    #           - f ki_allyl(T) u_active V_poly
                    #     Production: halved radical-disappearance
                    #     convention (kt_disp*R^2 IS the disproportionation
                    #     EVENT rate, 1 unsaturated end per event) -- ONLY
                    #     disproportionation sources U; recombination and
                    #     random initiation do not. NO f on production
                    #     (r36): R_ss already carries the escape
                    #     efficiency; termination is not a caged event.
                    #     CAPACITY THROTTLE (r36 P1-3): the linear factor
                    #     max(0, 1 - U/(2*mu0_amount)) is EXACTLY zero at
                    #     the chain-end capacity 2*mu0 (each chain has at
                    #     most 2 tail ends), so U(t) cannot integrate past
                    #     it -- the t=0 census guards the IC, this guards
                    #     the trajectory. mu0 here is the TAIL mu0 amount
                    #     [mol], the same distribution B is counted on; at
                    #     mu0 = 0 the throttle is 0 (no ends, no
                    #     capacity).
                    #     Sink: f-SYMMETRIC with the G_R term (r36 P1-1) --
                    #     f is the escaped-radical-pair efficiency, a
                    #     caged recombination restores the allylic bond,
                    #     so the caged fraction consumes no U.
                    #   Volume convention matches the moment ODEs below
                    #   exactly: concentration-rate * V_poly -> amount-rate
                    #   (dn_dt[idx_mu1] += dmu1_dt * V_poly).
                    u_amount = _softclamp(y[self.qssa_u_slot[pool_i]], self._softclamp_lam)
                    if B_qssa > 0.0:
                        u_active = u_amount / V_poly
                        if u_active > B_qssa:
                            u_active = B_qssa
                        T_qssa = self.T.value_si
                        RT_qssa = QSSA_R_GAS * T_qssa
                        ki_qssa = (self.qssa_ki_A[pool_i]
                                   * T_qssa ** self.qssa_ki_n[pool_i]
                                   * math.exp(-self.qssa_ki_Ea[pool_i] / RT_qssa))
                        kdp_qssa = (self.qssa_kdp_A[pool_i]
                                    * T_qssa ** self.qssa_kdp_n[pool_i]
                                    * math.exp(-self.qssa_kdp_Ea[pool_i] / RT_qssa))
                        kia_qssa = (self.qssa_kia_A[pool_i]
                                    * T_qssa ** self.qssa_kia_n[pool_i]
                                    * math.exp(-self.qssa_kia_Ea[pool_i] / RT_qssa))
                        ktrec_qssa = (self.qssa_ktrec_A[pool_i]
                                      * T_qssa ** self.qssa_ktrec_n[pool_i]
                                      * math.exp(-self.qssa_ktrec_Ea[pool_i] / RT_qssa))
                        ktdisp_qssa = (self.qssa_ktdisp_A[pool_i]
                                       * T_qssa ** self.qssa_ktdisp_n[pool_i]
                                       * math.exp(-self.qssa_ktdisp_Ea[pool_i] / RT_qssa))
                        kt_qssa = ktrec_qssa + ktdisp_qssa
                        # Same degenerate-rate posture as the legacy branch:
                        # kt_total underflowing to 0 divides R_ss by zero;
                        # any inf/NaN evaluation poisons the residual.
                        if not (ki_qssa < QSSA_INF and kdp_qssa < QSSA_INF
                                and kia_qssa < QSSA_INF
                                and ktrec_qssa < QSSA_INF
                                and ktdisp_qssa < QSSA_INF
                                and 0.0 < kt_qssa < QSSA_INF):
                            _raise_degenerate_qssa_rates(
                                pool.label, T_qssa, ki_qssa, kdp_qssa,
                                kt_qssa, kia=kia_qssa, ktrec=ktrec_qssa,
                                ktdisp=ktdisp_qssa)
                        fkiB = self.qssa_efficiency[pool_i] * ki_qssa * B_qssa
                        G_R_qssa = (2.0 * fkiB
                                    + self.qssa_efficiency[pool_i]
                                    * kia_qssa * u_active)
                        if self.qssa_has_transfer[pool_i]:
                            ktr_qssa = (self.qssa_ktr_A[pool_i]
                                        * T_qssa ** self.qssa_ktr_n[pool_i]
                                        * math.exp(-self.qssa_ktr_Ea[pool_i] / RT_qssa))
                            if not ktr_qssa < QSSA_INF:
                                _raise_degenerate_qssa_rates(
                                    pool.label, T_qssa, ki_qssa, kdp_qssa,
                                    kt_qssa, ktr_qssa, kia=kia_qssa,
                                    ktrec=ktrec_qssa, ktdisp=ktdisp_qssa)
                            R_ss = ((math.sqrt(ktr_qssa * ktr_qssa
                                               + 8.0 * kt_qssa * G_R_qssa)
                                     - ktr_qssa) / (4.0 * kt_qssa))
                        else:
                            R_ss = math.sqrt(G_R_qssa / (2.0 * kt_qssa))
                        r_qssa = (self.qssa_monomer_yield[pool_i]
                                  * kdp_qssa * R_ss)
                        # dU/dt (amount basis, mol/s): capacity-throttled
                        # production by the disproportionation branch ONLY
                        # (no f -- R_ss carries the escape efficiency);
                        # f-symmetric sink by allyl fission of the ACTIVE
                        # unsaturated ends. See the law block above.
                        cap_qssa = 2.0 * mu0_mol  # tail chain-end capacity
                        if cap_qssa > 0.0:
                            throttle_qssa = 1.0 - u_amount / cap_qssa
                            if throttle_qssa < 0.0:
                                throttle_qssa = 0.0
                        else:
                            throttle_qssa = 0.0
                        du_dt[self.qssa_u_slot[pool_i] - n_core] += (
                            ktdisp_qssa * R_ss * R_ss * throttle_qssa
                            - self.qssa_efficiency[pool_i]
                            * kia_qssa * u_active) * V_poly
                elif B_qssa > 0.0:
                    T_qssa = self.T.value_si
                    RT_qssa = QSSA_R_GAS * T_qssa
                    ki_qssa = (self.qssa_ki_A[pool_i]
                               * T_qssa ** self.qssa_ki_n[pool_i]
                               * math.exp(-self.qssa_ki_Ea[pool_i] / RT_qssa))
                    kdp_qssa = (self.qssa_kdp_A[pool_i]
                                * T_qssa ** self.qssa_kdp_n[pool_i]
                                * math.exp(-self.qssa_kdp_Ea[pool_i] / RT_qssa))
                    kt_qssa = (self.qssa_kt_A[pool_i]
                               * T_qssa ** self.qssa_kt_n[pool_i]
                               * math.exp(-self.qssa_kt_Ea[pool_i] / RT_qssa))
                    # Runtime degenerate-rate guard (fail-loud, cheap hot
                    # path). M1 pins config finiteness, but EVALUATION can
                    # still degenerate at solver T: kt underflows to 0.0
                    # (sqrt(f*ki*B/kt) divides by zero -> inf) or a large
                    # A/n overflows to inf/NaN. `x < QSSA_INF` rejects inf
                    # AND NaN in one compare; ki/kdp underflow to 0.0 is
                    # fine (the rate smoothly reaches 0), kt is not.
                    if not (ki_qssa < QSSA_INF and kdp_qssa < QSSA_INF
                            and 0.0 < kt_qssa < QSSA_INF):
                        _raise_degenerate_qssa_rates(
                            pool.label, T_qssa, ki_qssa, kdp_qssa, kt_qssa)
                    # A parent and its X-loss daughters are one melt with
                    # one radical population (I-260): R from the lineage's
                    # backbone bonds, released in proportion to the bonds
                    # (the radicals) this pool holds. One pool: unchanged.
                    B_eval_qssa = B_qssa
                    if qssa_lineage_B is not None and \
                            self.qssa_lineage_n[self.qssa_root[pool_i]] > 1:
                        B_eval_qssa = qssa_lineage_B[self.qssa_root[pool_i]]
                    fkiB = self.qssa_efficiency[pool_i] * ki_qssa * B_eval_qssa
                    G_R_qssa = 2.0 * fkiB  # two end radicals per homolysis
                    if self.qssa_has_transfer[pool_i]:
                        ktr_qssa = (self.qssa_ktr_A[pool_i]
                                    * T_qssa ** self.qssa_ktr_n[pool_i]
                                    * math.exp(-self.qssa_ktr_Ea[pool_i] / RT_qssa))
                        if not ktr_qssa < QSSA_INF:
                            _raise_degenerate_qssa_rates(
                                pool.label, T_qssa, ki_qssa, kdp_qssa,
                                kt_qssa, ktr_qssa)
                        R_ss = ((math.sqrt(ktr_qssa * ktr_qssa
                                           + 8.0 * kt_qssa * G_R_qssa)
                                 - ktr_qssa) / (4.0 * kt_qssa))
                    else:
                        # == sqrt(G_R / (2 kt)), simplified
                        R_ss = math.sqrt(fkiB / kt_qssa)
                    r_qssa = self.qssa_monomer_yield[pool_i] * kdp_qssa * R_ss
                    if B_eval_qssa != B_qssa:
                        r_qssa = r_qssa * (B_qssa / B_eval_qssa)
                if r_qssa > 0.0:
                    # Chain-END monomer release signature: mu1 drains one
                    # unit per release, mu2 drains (2 E[n] - 1) per release
                    # with the same-pool-VE clamp (>0 only: the drain must
                    # never make mu2 increase; mu0 ~ 0 guarded by the eps
                    # clamp).
                    # monomer_yield already scales r_qssa, so the moment
                    # drain and the gas emission below scale TOGETHER --
                    # scaling only one side would fabricate/destroy mass.
                    # I-055 (see the k_unzip kernel above for the derivation).
                    # This channel used to leave mu0 untouched ("no
                    # chain created/destroyed"), which is wrong for the same
                    # reason it was wrong on the legacy k_unzip channel
                    # above -- a release that consumes a DP=1 chain's LAST
                    # repeat unit destroys the CHAIN, not just the unit.
                    # Without the debit mu1 is driven down through mu0, out
                    # of the realizable cone mu1 >= mu0 >= 0, and then
                    # negative: the accepted-state r81 raise that killed run
                    # poly_104. A fraction N1 of releases land on a DP=1
                    # chain, so the debit is r_qssa*N1 -- the SAME
                    # (release rate)*N1 law the legacy channel carries.
                    # It must be applied to both channels or to neither:
                    # the solver pins them equivalent
                    # (test_qssa_handshake_equivalence_with_k_unzip), and
                    # fixing only one breaks that pin.
                    # I-065 defect 1, sibling channel. The availability gate
                    # carried by the legacy channel above is applied here
                    # too, because the two are pinned equivalent and a
                    # release law must hold on both or on neither. On THIS
                    # channel it is provably dead code today: the release
                    # rate is r_qssa ~ sqrt(B_qssa) with
                    # B_qssa = max(mu1 - mu0, 0) and the whole block is
                    # gated `elif B_qssa > 0.0`, so the channel is already
                    # identically zero everywhere the gate is < 1, and
                    # exactly 1.0 (bit-for-bit no-op) everywhere it fires.
                    # It is kept because that safety is B_qssa's accident,
                    # not this kernel's contract: loosen the initiation
                    # gate and the fabrication returns.
                    r_release_qssa = r_qssa * _release_units_gate(mu0, mu1)
                    dmu1_dt -= r_release_qssa
                    dmu0_dt -= r_qssa * _deprop_dp1_fraction(mu0, mu1, mu2)
                    qssa_mu2_dec = 2.0 * (mu1 / max(mu0, SMALL_EPS)) - 1.0
                    if qssa_mu2_dec > 0.0:
                        dmu2_dt -= r_release_qssa * qssa_mu2_dec
                    # monomer_poly_index is non-None whenever enabled (M1
                    # invariant); emission flows through the SAME small_src
                    # -> dn_dt * V_poly path as the k_unzip channel, i.e.
                    # to the GAS species amount basis (incident 2026-07-03,
                    # design B-prime).
                    small_src[pool.monomer_poly_index] = (
                        small_src.get(pool.monomer_poly_index, 0.0)
                        + r_release_qssa)
                    self._split_side_loss_release(
                        pool_i, pool.monomer_poly_index, r_release_qssa,
                        mu1_mol, V_poly, y, small_src, dz_dt)

            # End-radical DEPROPAGATION kernel (adjudicated round 74 SS2, the
            # run-6 no-outlet wall fix). Reads ONLY the flattened kdep_*
            # arrays from _flatten_depropagation_state -- never the pool
            # dict. Gates on kdep_enabled; mutually exclusive with
            # k_unzip > 0 / radical_qssa_unzip / k_homolysis on the same
            # pool (validate_configuration hard-errors), so at most one
            # chain-end release arm fires per pool.
            #
            # Law, per radical-end pool (ONE active radical end per chain),
            # k_dep(T) = A*T^n*exp(-Ea/(R_gas*T)) at the RUNTIME T:
            #   R    = k_dep * mu0 * g   unzip events == monomer release
            #        (g carries BOTH the r74 exhaustion gate and the I-065
            #         availability gate; both are exactly 1 for mean >= 1)
            #   gas  = +R at kdep_gas    (the SAME float as the mu1 drain:
            #                             d(condensed) + d(gas monomer) = 0
            #                             EXACTLY under MW multiplication)
            #   dmu1 = -R
            #   dmu2 = -k_dep*(2*mu1 - mu0), computed as
            #          -k_dep*mu0*(g + 2*sp(mean - 1)): a DP=n chain loses
            #          (2n - 1) from mu2 per event; the smooth positive-part
            #          sp keeps the excess term nonnegative without a
            #          max(...,0) cliff (identical to the nominal law up to
            #          O(W^2/(mean-1)) in the realizable region).
            #   dmu0 = -k_dep * N1, N1 = mu0 * p1 from the EXISTING
            #          discrete/gamma closure + smooth terminal boundary
            #          floor (_deprop_dp1_fraction): DP=1 chains release
            #          their last unit and TERMINATE -- never a permanent
            #          dmu0 = 0 (r74: that stalls at a one-repeat-per-chain
            #          residue).
            # SMOOTH exhaustion gate (r74 SS5, designed in from the start):
            # g = 1 - sp(1 - mean) is EXACTLY 1 for mean >= 1 and rolls off
            # C2-smoothly below; dmu0 is deliberately NOT gated by g, so a
            # noise state mean < 1 drains chains FASTER than units
            # (self-healing back toward the realizable cone, no stiff
            # no-outlet grind).
            if self.kdep_enabled[pool_i] and mu0 > 0.0:
                T_kdep = self.T.value_si
                k_dep = (self.kdep_A[pool_i]
                         * T_kdep ** self.kdep_n[pool_i]
                         * math.exp(-self.kdep_Ea[pool_i]
                                    / (QSSA_R_GAS * T_kdep)))
                if not (0.0 < k_dep < QSSA_INF):
                    raise ValueError(
                        f"Pool {pool.label}: k_depropagation(T={T_kdep:g} K) "
                        f"evaluated to a degenerate rate constant "
                        f"({k_dep!r}) -- A={self.kdep_A[pool_i]:g}, "
                        f"n={self.kdep_n[pool_i]:g}, "
                        f"Ea={self.kdep_Ea[pool_i]:g} J/mol. Refusing to "
                        f"integrate a poisoned kernel.")
                mean_kdep = mu1 / mu0
                # I-065 defect 1, sibling channel. The r74 SS5 exhaustion
                # gate g is exactly 1 for mean >= 1 but bottoms out at
                # W^2/(1 + W^2) = 1e-4, so at the mu1 = 0, mu0 > 0 noise
                # state this channel still emitted k_dep*mu0*1e-4 of gas
                # with no repeat unit behind it -- bounded fabrication, but
                # fabrication. _release_units_gate takes it to EXACTLY zero
                # there while staying EXACTLY 1.0 (bit-for-bit no-op) for
                # mean >= 1, and is C2 at the boundary, so the no-cliff
                # property g was designed for survives
                # (test_deprop_smooth_exhaustion_gate_no_cliff).
                # dmu0 below stays UNGATED, as r74 requires: chains keep
                # draining at -k_dep*mu0 so a mean < 1 state heals back
                # toward the cone rather than grinding.
                g_kdep = ((1.0 - _smooth_pos(1.0 - mean_kdep, KDEP_GATE_WIDTH))
                          * _release_units_gate(mu0, mu1))
                r_kdep = k_dep * mu0 * g_kdep
                if r_kdep > 0.0:
                    dmu1_dt -= r_kdep
                    dmu2_dt -= k_dep * mu0 * (
                        g_kdep + 2.0 * _smooth_pos(mean_kdep - 1.0,
                                                   KDEP_GATE_WIDTH))
                    g_idx_kdep = int(self.kdep_gas[pool_i])
                    small_src[g_idx_kdep] = (
                        small_src.get(g_idx_kdep, 0.0) + r_kdep)
                    self._split_side_loss_release(
                        pool_i, g_idx_kdep, r_kdep, mu1_mol, V_poly, y,
                        small_src, dz_dt)
                dmu0_dt -= k_dep * mu0 * _deprop_dp1_fraction(mu0, mu1, mu2)

            # Hybrid Handshake
            tail_mean = mu1 / mu0 if mu0 > SMALL_EPS else 0.0
            # I-060. `mu0 > TAIL_CONC_MIN` is the half that does real work
            # and it is kept unchanged: below 1e-9 mol/m^3 the tail is
            # numerically empty, every quantity the handshake computes comes
            # through the division mu1/mu0, and depositing that quotient into
            # a REAL core species would seed an oligomer out of round-off.
            #
            # The other half, `tail_mean > xs + 1e-9`, asked a real question --
            # can the pool afford the chain the handshake is about to remove?
            # -- and answered it with an average, as a hard boolean, on the
            # wrong side of the flux's own maximum. Measured on a scission-fed
            # k_unzip pool (xs=3, k_unzip/k_scission = 5, PDI 1.5): the flux is
            # LARGEST just above the threshold and the boolean took it to zero
            # there -- F(mean = xs + 1e-6) = 3.8390e-02, F(mean = xs) = 0.0,
            # a finite JUMP in the residual on a surface the trajectory
            # crosses, with dmu0 flipping sign across it
            # (-2.1377e-02 -> +1.7013e-02). That is a zeroth-order
            # discontinuity, strictly worse than the first-derivative kink the
            # _release_units_gate quintic exists to avoid, and downstream of
            # it the handshake -- the tail's ONLY outlet into the explicit
            # ladder -- stayed off for the rest of the run: the DP=xs species
            # froze and the whole sub-cutoff residue left as gas monomer.
            #
            # The boolean is therefore removed here and replaced, at the flux
            # below, by explicit per-state BUDGETS on the boundary population
            # -- one from each realizability inequality, derived there. An
            # average over the whole tail was never the right instrument: the
            # question is how many chains sit at DP = xs, and mean DP > xs is
            # neither necessary nor sufficient for that.
            #
            # Not to be confused with the tail's REPRESENTATION invariant,
            # which is a different inequality: the tail holds chains with
            # n > xs (module docstring 'Moment Tail: ... chains with DP > xs';
            # polymer_input.PolymerPool xs doc), so its support starts at xs+1
            # and it satisfies mu1 >= (xs+1)*mu0. The handshake is the only
            # term that pushes THAT quantity back up -- dmu0 -= F and
            # dmu1 -= xs*F give d(mu1 - (xs+1)*mu0)/dt += +F exactly -- which
            # is why losing it below the boundary is a defect and not merely
            # a branch not taken.
            #
            # Removing the condition and putting nothing in its place was
            # tried and is WRONG: at full strength below the boundary the
            # handshake removes chains it has no units for, and the pool is
            # driven clean out of the mu1 >= mu0 cone (min(mu1 - mu0) =
            # -1.51e-01 on the trajectory above, against +1.85e-20 unfixed).
            # So was a smooth availability FRACTION -- the budgets below are
            # absolute chain counts, and a fraction of the wrong base is not
            # a bound (see the derivation at the flux clamp).
            valid_tail = mu0 > TAIL_CONC_MIN

            # Per-chain unzip frequency feeding the handshake: the legacy
            # k_unzip IS that frequency; the QSSA equivalent is
            # k = r_mono / mu0 (release events per chain per second), so
            # QSSA pools do not strand low-DP condensed residue. The two
            # channels are mutually exclusive (M1), so at most one arm fires.
            k_chain_handshake = 0.0
            if pool.k_unzip > 0:
                k_chain_handshake = pool.k_unzip
            elif self.qssa_enabled[pool_i] and r_qssa > 0.0:
                k_chain_handshake = r_qssa / max(mu0, SMALL_EPS)

            if valid_tail and k_chain_handshake > 0.0:
                explicit_xs_idx = pool.explicit_dp_to_species_index.get(xs, None)
                if explicit_xs_idx is not None:
                    params = _gamma_params_from_mu012(mu0, mu1, mu2)
                    p_cond = 0.0

                    if params:
                        k_shape, theta = params
                        p_cond = _gamma_prob_conditional_hybrid(xs + 1, xs, k_shape, theta)
                    else:
                        # I-098 MONODISPERSE FALLBACK.
                        #
                        # The closure has declined -- any moment <=
                        # SMALL_EPS, PDI <= 1 + 1e-6, or non-finite params --
                        # so there is no fitted distribution here and the
                        # distributional assumption has to be STATED. It is
                        # the substance of this branch, not a footnote.
                        #
                        # What p_cond means. The gamma leg above computes
                        # P(DP = xs+1 | DP > xs): the fraction of tail chains
                        # in the BOUNDARY BIN. The tail's support starts at
                        # xs + 1 (module docstring; the representation
                        # invariant mu1 >= (xs+1)*mu0 is spelled out at the
                        # valid_tail note above), so the boundary bin IS the
                        # support's minimum. Hence tail_mean = xs + 1 is the
                        # mean equalling the minimum, which on a non-negative
                        # lattice forces EVERY chain onto that minimum:
                        # p_cond = 1 exactly, with no assumption at all.
                        #
                        # The triangle that stood here returned 0.0 there and
                        # peaked at xs + 1.5, i.e. it ran monotone in the
                        # wrong direction across the entire lower half of its
                        # support. Measured against this branch's own inputs:
                        # 0.0000 at xs+1 (truth 1.0), 2.0e-12 as
                        # tail_mean -> (xs+1)+, 1.0 at xs+1.5 (truth 0.5).
                        #
                        # THE ASSUMPTION. The closure refuses precisely
                        # because the variance collapsed, so take the
                        # MINIMUM-VARIANCE distribution on the tail's integer
                        # support {xs+1, xs+2, ...} with this mean. On a
                        # lattice that distribution is unique -- all mass on
                        # the two sites bracketing the mean, split so the
                        # mean is reproduced -- and its mass at xs + 1 is
                        # 1 - (tail_mean - (xs+1)) on [xs+1, xs+2], 0 above,
                        # 1 below (a mean under the support's minimum is not
                        # realizable; the nearest state that is puts
                        # everything on the minimum). That is the single
                        # clamped expression below. Its only constants are
                        # the lattice offsets 1 and 2, in units of one repeat
                        # unit; nothing is fitted, and nothing needs a scale
                        # this function does not already carry.
                        #
                        # Why this is not merely one defensible choice among
                        # many: the trigger is itself a realizability filter.
                        # A lattice distribution with mean n + f has variance
                        # >= f*(1-f), so PDI - 1 >= f*(1-f)/mean^2, and
                        # reaching PDI <= 1 + 1e-6 needs f <= 1.6e-5 at
                        # xs = 3 (4.0e-6 at xs = 1, 2.6e-3 at xs = 50). Every
                        # state that can actually arrive here therefore sits
                        # within ~1e-5 of an INTEGER mean, where the answer is
                        # not an interpolation at all: 1 at xs+1, 0 at xs+2
                        # and beyond. The old triangle returned 0 at BOTH of
                        # those. The interior is lattice-unrealizable, and the
                        # linear form is the unique CONTINUOUS interpolant of
                        # the two realizable endpoints -- continuity is not
                        # optional, since a step would put a zeroth-order jump
                        # in the residual on a surface the trajectory crosses,
                        # which is the defect the I-060 note above exists to
                        # record.
                        #
                        # Meeting the gamma leg. As PDI -> 1+ the gamma leg
                        # tends to nearest-lattice-site ROUNDING -- its
                        # half-integer bins give 1.0 below xs+1.5, 0.0 above,
                        # 0.500133 at it -- and this bracket is exactly that
                        # step's continuous interpolant, so the two agree at
                        # frac = 0, 0.5 and 1. Measured across the trigger at
                        # fixed mean (PDI = 1+1e-6 vs 1+1e-6+1e-9), over the
                        # lattice-realizable states that can reach it,
                        # sup|fallback - gamma| falls from 1.000000 to
                        # 1.0e-05. Over the whole interval including the
                        # unrealizable interior it falls from 1.000000 to
                        # 0.475000; that residue is the gamma leg's own
                        # rounding step at frac = 0.5, which a continuous law
                        # cannot reproduce and should not try to.
                        #
                        # Monotone non-increasing, as the truth is: the mass
                        # on the lowest site can only fall as the mean leaves
                        # it. The triangle was not.
                        p_cond = max(0.0, min(1.0, (xs + 2.0) - tail_mean))

                    p_cond = min(1.0, max(0.0, p_cond))

                    # Population (conc) at boundary
                    N_boundary = mu0 * p_cond

                    # Flux Clamp
                    N_boundary = min(N_boundary, mu0)
                    if xs > 0:
                        N_boundary = min(N_boundary, mu1 / xs)
                        N_boundary = min(N_boundary, mu2 / (xs * xs))
                    # I-060 REALIZABILITY BUDGETS.
                    #
                    # p_cond above is the gamma closure's ESTIMATE of what
                    # fraction of the tail sits at the boundary. The three
                    # clamps above are the elementary bounds on the same
                    # population: mu0 = SUM c_n >= c_xs, mu1 = SUM n c_n >=
                    # xs*c_xs, mu2 = SUM n^2 c_n >= xs^2*c_xs. The two below
                    # are the bounds the two REALIZABILITY inequalities give,
                    # and neither is implied by those three. Each is an
                    # absolute upper bound on c_xs; each is ATTAINED by an
                    # explicit realizable distribution, so neither is slack
                    # that could be dropped; and a handshake respecting all
                    # five can never remove chains the pool provably is not
                    # holding.
                    #
                    # (a) From mu1 >= mu0 -- every chain carries at least one
                    #     repeat unit:
                    #         mu1 - mu0 = SUM (n - 1) c_n >= (xs - 1)*c_xs
                    #     so c_xs <= (mu1 - mu0)/(xs - 1). One event spends
                    #     (xs - 1) of that excess, so this IS the excess
                    #     budget, and it makes mu1 = mu0 an invariant set:
                    #     d(mu1 - mu0)/dt carries -(xs - 1)*F with F now
                    #     bounded by k*(mu1 - mu0), which vanishes there. The
                    #     I-055 chain-termination debit cannot supply that --
                    #     its own release term vanishes at the edge (p1 -> 1)
                    #     while F did not. xs <= 1 spends no excess per event,
                    #     so no budget is owed and none is applied.
                    #
                    # (b) From mu0*mu2 >= mu1^2 (Cauchy-Schwarz). Take moments
                    #     about xs, m_k = SUM (n - xs)^k c_n. Then
                    #     Q := mu0*mu2 - mu1^2 = m_0*m_2 - m_1^2, since Q is
                    #     mu0^2 * Var and so translation invariant. A
                    #     handshake event takes m_0 -> m_0 - F and leaves m_1
                    #     and m_2 untouched, because it removes chains at
                    #     exactly n = xs, where (n - xs) = 0. Hence
                    #         dQ/dt = -F * m_2 = -F * SUM (n - xs)^2 c_n <= 0
                    #     for EVERY realizable state -- unconditional, not a
                    #     corner case. That sign is not itself the defect:
                    #     deleting a monodisperse sub-population is supposed
                    #     to cut the variance, and what is left of a real
                    #     distribution is still a real distribution. The
                    #     defect is running the drain past its budget. Split
                    #     m_0 = c_xs + r over the n = xs atom and the rest;
                    #     m_1 and m_2 see only the rest, so
                    #         Q = c_xs*m_2 + (r*m_2 - m_1^2) >= c_xs*m_2
                    #     by Cauchy-Schwarz on the restricted measure, giving
                    #         c_xs <= Q / m_2.
                    #     TIGHT: for an atom c_xs at xs plus an atom at
                    #     xs + d it returns exactly c_xs. Capping here makes
                    #     Q = 0 an invariant set (dQ/dt >= -k*Q) without
                    #     clamping Q itself, which would hide the drain
                    #     rather than stop it. m_2 = 0 means the pool is
                    #     monodisperse at exactly xs: Q = 0, dQ/dt = 0 for any
                    #     F, no budget is owed and the whole pool may
                    #     legitimately cross -- so the bound is skipped rather
                    #     than evaluated as 0/0.
                    #
                    # On the BASIS question these replace. Round 1 wrote (a) as
                    # _release_units_gate((xs-1)*mu0, mu1-mu0), i.e. as a
                    # FRACTION with mu0 in the denominator, and review asked
                    # whether the denominator should have been N_boundary --
                    # the population the rate is actually drawn from -- rather
                    # than all tail chains. By _release_units_gate's own
                    # semantics (units on offer / units the event rate implies)
                    # review is right that mu0 is the wrong base: the demand is
                    # (xs-1)*N_boundary. But N_boundary is not the fix either.
                    # A gate of the form N*h(B/N) is self-referential -- it
                    # DECREASES as the attempt grows, going to zero like
                    # 10*B^3/N^2 -- and it still is not a bound, since
                    # h(t) > t for t > 0.5. The error was not the choice of
                    # base; it was treating a BUDGET as a fraction at all. A
                    # budget is an absolute count of chains, so the operator
                    # that respects it is min(attempt, budget), which needs no
                    # denominator and is what is written here.
                    #
                    # Hard min rather than the C2 quintic _release_units_gate
                    # carries for the 1-unit release channels, and the choice
                    # is forced: a smooth multiplier that is exactly 1.0 above
                    # the boundary must exceed the budget just below it (the
                    # quintic has h(t) > t for t > 0.5), so it cannot be both
                    # inert above and a bound below. Round 1 shipped the
                    # quintic and was a bound only by measurement, not by
                    # construction. Inertness is not given up: (a) is exactly
                    # slack for mean >= xs, where (mu1 - mu0)/(xs - 1) >= mu0
                    # >= N_boundary. (b) is NOT slack there and does change
                    # the flux above the cutoff -- deliberately, since the
                    # over-budget removal is a pre-existing defect that the
                    # old boolean never addressed (measured: 8 of 77 swept
                    # states on the unfixed build, all at mean >= xs + 0.5).
                    # The kink a binding min leaves is first-order, against
                    # the zeroth-order jump the old boolean put on a surface
                    # the trajectory crosses; the three clamps above are
                    # already min() in this same expression; and the no-cliff
                    # pin (test_deprop_smooth_exhaustion_gate_no_cliff) is
                    # about the release rate at mean = 1, not this population.
                    if xs > 1:
                        N_boundary = min(N_boundary, (mu1 - mu0) / (xs - 1))
                    m2_about_xs = mu2 - 2.0 * xs * mu1 + xs * xs * mu0
                    q_cone = mu0 * mu2 - mu1 * mu1
                    if m2_about_xs > 0.0:
                        N_boundary = min(N_boundary, q_cone / m2_about_xs)
                    elif m2_about_xs < 0.0 or q_cone < 0.0:
                        # m2_about_xs is SUM (n - xs)^2 c_n, a sum of squares
                        # against non-negative concentrations: it cannot be
                        # negative for ANY realizable tail. Reaching here means
                        # the triple is already off the cone, so no budget
                        # above describes a realizable population -- and the
                        # one that would have caught it was skipped by the
                        # m2 <= 0 test itself. Do NOT read the other clamps as
                        # protection here: at xs = 2, (mu0, mu1, mu2) =
                        # (1, 3.5, 9) has q_cone = -3.25 and m2_about_xs = -1,
                        # every other budget is positive, and before this
                        # branch existed the handshake ran there at its full
                        # unclamped rate (measured F = 1.000000e-01).
                        #
                        # m2_about_xs == 0 with q_cone == 0 is the one
                        # degenerate case that IS realizable -- the whole tail
                        # sitting exactly at n = xs -- and it falls through to
                        # the clamps above, which are correct for it.
                        # m2_about_xs == 0 forces q_cone = -m1^2 <= 0, so
                        # q_cone > 0 cannot occur in this branch at all.
                        #
                        # Refusing the flux is the conservative reading, not a
                        # diagnosis. Whether an accepted off-cone state should
                        # be fatal rather than warned about belongs to the
                        # variance census, which owns that call.
                        N_boundary = 0.0
                    # A negative N_boundary would run the handshake BACKWARDS,
                    # inventing tail chains out of the explicit species. It
                    # arises where a budget IS computed and comes out negative:
                    # q_cone < 0 with m2_about_xs > 0, or mu1 < mu0.
                    if N_boundary < 0.0:
                        N_boundary = 0.0

                    F = k_chain_handshake * N_boundary

                    if F > 0.0:
                        dn_dt[explicit_xs_idx] += F * V_poly
                        dmu0_dt -= F
                        # I-141: the explicit oligomer species carries its own
                        # end groups, so this chain loss orphans nothing; the
                        # H counter lets the ledger subtract it from the
                        # mu0 decrease (L = mu0(0) - mu0 - H).
                        if self.end_ledger_rel[pool_i] >= 0:
                            dl_dt[self.end_ledger_rel[pool_i]] += F * V_poly
                        dmu1_dt -= xs * F
                        dmu2_dt -= (xs * xs) * F

            dn_dt[idx_mu0] += dmu0_dt * V_poly
            dn_dt[idx_mu1] += dmu1_dt * V_poly
            dn_dt[idx_mu2] += dmu2_dt * V_poly

            for spc_idx, r_vol in small_src.items():
                if 0 <= spc_idx < n_core:
                    dn_dt[spc_idx] += r_vol * V_poly

        # 9. Mass Transfer
        for mt in self.mass_transfer:
            ig = mt.gas_index
            ipoly = mt.poly_index

            Cg = C_gas[ig]
            Cp = C_poly[ipoly]

            J = mt.kLa * (Cp - mt.K * Cg)
            dn = J * V_poly

            dn_dt[ig] += dn
            dn_dt[ipoly] -= dn

        # 10. Constants
        if self.const_spc_indices is not None:
            for i in self.const_spc_indices:
                dn_dt[i] = 0.0

        # 11. Diagnostics
        self.core_species_rates[mask] = dn_dt[mask] / V_gas
        self.core_species_rates[~mask] = dn_dt[~mask] / V_poly

        for i in range(n_core):
            if self.is_pool_proxy[i]:
                self.core_species_rates[i] = proxy_activity[i] / V_poly

        # --- DEBUG: report first non-finite / blow-up in residual ---
        if not getattr(self, '_dbg_reported', False):
            if (not np.isfinite(V_gas)) or (not np.all(np.isfinite(dn_dt))):
                self._dbg_reported = True
                bad = np.where(~np.isfinite(dn_dt))[0]
                print(f"\n[POLY-DBG] NON-FINITE in residual at t={t:.3e}")
                print(f"[POLY-DBG]   V_gas={V_gas:.6e} V_poly={V_poly:.6e} "
                      f"n_gas={float(np.sum(y[:n_core][mask])):.6e}")
                cg = C_gas[mask]
                print(f"[POLY-DBG]   max|C_gas|={float(np.max(np.abs(cg))) if cg.size else 0.0:.4e} "
                      f"max|dn_dt finite|={float(np.max(np.abs(dn_dt[np.isfinite(dn_dt)]))) if np.any(np.isfinite(dn_dt)) else 0.0:.4e}")
                print(f"[POLY-DBG]   non-finite dn_dt indices (first 10): {list(bad[:10])}")
                for bi in bad[:10]:
                    print(f"[POLY-DBG]     idx {bi}: y={y[bi]:.4e} gas={bool(mask[bi])}")

        if self.num_qssa_u > 0 or self.num_sgh_z > 0 or self.num_end_ledger > 0:
            # Appended-state layout: y = [core species | U slots | Z slots].
            # The residual vector follows the SAME layout, so the blocks are
            # concatenated in slot-index order (core | U | Z) to match the
            # absolute slot assignments in flatten. Zero-length blocks are
            # harmless (concatenate drops them), but the legacy
            # (num_qssa_u == 0 AND num_sgh_z == 0) return below stays the
            # bit-for-bit pre-milestone expression -- a golden bitwise-RHS
            # test depends on it.
            return (np.concatenate((dn_dt, du_dt, dz_dt, dl_dt)) - dydt), 1
        return (dn_dt - dydt), 1

    # ------------------------------------------------------------------
    # I-141 chain end-group bookkeeping: reporters (read-only; none of this
    # feeds a rate).
    # ------------------------------------------------------------------
    @staticmethod
    def _collect_species_atoms(core_species):
        """Element counts per core species ({} when the species carries no
        structure, e.g. a label-only consumer-world species)."""
        out = []
        for sp in core_species:
            try:
                out.append(dict(sp.molecule[0].get_element_count()))
            except Exception:
                out.append({})
        return out

    def _ledger_state(self, y):
        if y is None:
            y = self.y if getattr(self, "y", None) is not None else self.y0
        return np.asarray(y, dtype=float)

    def _require_end_group_atoms(self):
        missing = [p.label for p in self.polymer_pools if not p.end_group_atoms]
        if missing:
            raise ValueError(
                f"Pool(s) {missing} declare no end_group_atoms: the chain "
                f"end-group ledger / physical inventory will not guess a "
                f"composition. Declare end_group_atoms (and repeat_unit_atoms) "
                f"on every pool.")

    def end_group_ledger(self, y=None):
        """Named output: {pool label: mol of chains' worth of that pool's end
        groups that left the moment accounting (net of chains minted)}:
        L_p = mu0_p(0) - mu0_p(y) - H_p. Every mu0 change by any channel
        (scission/homolysis mint chains: L falls; recombination and release
        to gas remove them: L rises; cross-pool transfers move L between
        pools) is booked, except chains handed to an explicit-oligomer
        species (H_p), whose ends travel with that species."""
        self._require_end_group_atoms()
        y = self._ledger_state(y)
        y0 = np.asarray(self.y0, dtype=float)
        base = self.num_core_species + self.num_qssa_u + self.num_sgh_z
        out = {}
        for i, pool in enumerate(self.polymer_pools):
            m0 = int(pool.mu_indices[0])
            h = 0.0
            if self.end_ledger_rel[i] >= 0:
                h = float(y[base + int(self.end_ledger_rel[i])])
            out[pool.label] = float(y0[m0] - y[m0] - h)
        return out

    def orphaned_end_atoms(self, y=None):
        """Named output "orphaned end atoms": {element: mol} = sum over pools
        of ledger_p * declared end-group composition of pool p. Compare with
        the small products of an atom-explicit simulation."""
        led = self.end_group_ledger(y)
        out = {}
        for pool in self.polymer_pools:
            if pool.label in led:
                for el, n in pool.end_group_atoms:
                    out[el] = out.get(el, 0.0) + led[pool.label] * n
        return out

    def code_convention_inventory(self, y=None, species_atoms=None,
                                  defect_atoms=None):
        """Elemental inventory under the moment code's own convention:
        non-moment core species (gas and any condensed explicit species) plus
        sum_pools mu1*repeat_unit - mu0*per-chain defect atoms. Chain end
        groups are NOT included. ``species_atoms`` {label: {el: n}} overrides
        / supplies compositions; ``defect_atoms`` {pool label: {el: n}} gives
        the elements lost per chain by an X-loss feature pool (required when
        its chain_mass_defect_g_mol > 0)."""
        y = self._ledger_state(y)
        species_atoms = species_atoms or {}
        defect_atoms = defect_atoms or {}
        n_core = self.num_core_species
        moment = set()
        for pool in self.polymer_pools:
            moment.update(int(i) for i in pool.mu_indices)
        out = {}
        for i in range(n_core):
            if i in moment or self.is_pool_proxy[i] or y[i] == 0.0:
                continue
            comp = species_atoms.get(self._core_labels[i]) \
                if getattr(self, "_core_labels", None) else None
            if comp is None:
                comp = self._core_species_atoms[i]
            if not comp:
                raise ValueError(
                    f"Core species index {i} ({getattr(self, '_core_labels', ['?'] * n_core)[i]}) holds "
                    f"{y[i]:g} mol but has no element composition; pass "
                    f"species_atoms={{label: {{element: count}}}}.")
            for el, n in comp.items():
                out[el] = out.get(el, 0.0) + y[i] * n
        for pool in self.polymer_pools:
            if not pool.repeat_unit_atoms:
                raise ValueError(
                    f"Pool {pool.label}: repeat_unit_atoms is not declared; "
                    f"the elemental inventory cannot be reported.")
            mu0 = float(y[pool.mu_indices[0]])
            mu1 = float(y[pool.mu_indices[1]])
            for el, n in pool.repeat_unit_atoms:
                out[el] = out.get(el, 0.0) + mu1 * n
            if pool.chain_mass_defect_g_mol > 0.0:
                if pool.label not in defect_atoms:
                    raise ValueError(
                        f"Pool {pool.label}: chain_mass_defect_g_mol="
                        f"{pool.chain_mass_defect_g_mol:g} > 0 but no "
                        f"defect_atoms[{pool.label!r}] (elements lost per "
                        f"chain) was supplied.")
                for el, n in defect_atoms[pool.label].items():
                    out[el] = out.get(el, 0.0) - mu0 * n
        # An inventory-accounting side-group carrier keeps its moments while
        # its groups leave: the groups lost so far, max(0, sites*mu1 - Z) per
        # row, are off its chains (they are in the gas).
        sgh_en = getattr(self, "sgh_enabled", None)
        for p, pool in enumerate(self.polymer_pools):
            if sgh_en is None or p >= len(sgh_en) or not sgh_en[p]:
                continue
            if self.sgh_in_pool[self.sgh_row_start[p]]:
                # An X-loss daughter's further losses, s*mu1 - Z - mu0 per
                # row (I-260; the first per chain is its defect, above).
                for row, lost in self._side_loss_further(p, y)[0]:
                    for el, n in self._core_species_atoms[int(self.sgh_gas[row])].items():
                        out[el] = out.get(el, 0.0) - lost * n
                continue
            mu1 = max(0.0, float(y[pool.mu_indices[1]]))
            for row in range(int(self.sgh_row_start[p]), int(self.sgh_row_end[p])):
                z_slot = int(self.sgh_z_slot[row])
                if z_slot < 0:
                    continue
                lost = max(0.0, self.sgh_sites[row] * mu1 - float(y[z_slot]))
                for el, n in self._core_species_atoms[int(self.sgh_gas[row])].items():
                    out[el] = out.get(el, 0.0) - lost * n
        return out

    def physical_elemental_inventory(self, y=None, species_atoms=None,
                                     defect_atoms=None):
        """Physical elemental inventory = code convention + sum_pools
        mu0*(declared end groups) + orphaned end atoms (the ledger). With the
        ledger, this is conserved by the moment model for any declared
        linear polymer; the ledger itself is orphaned_end_atoms()."""
        y = self._ledger_state(y)
        out = dict(self.code_convention_inventory(y, species_atoms, defect_atoms))
        for pool in self.polymer_pools:
            if not pool.end_group_atoms:
                raise ValueError(
                    f"Pool {pool.label}: end_group_atoms is not declared; the "
                    f"physical inventory cannot be reported.")
            mu0 = float(y[pool.mu_indices[0]])
            for el, n in pool.end_group_atoms:
                out[el] = out.get(el, 0.0) + mu0 * n
        for el, v in self.orphaned_end_atoms(y).items():
            out[el] = out.get(el, 0.0) + v
        return out

    def get_reaction_rates(self, y_in):
        """
        Return the net rates of all core and edge reactions at state y.
        """
        # 1. Cast Input to Memoryview
        cdef np.ndarray[np.float64_t, ndim=1] y_arr = np.asarray(y_in, dtype=np.float64)
        cdef double[:] y = y_arr

        cdef int n_core = self.num_core_species
        cdef int n_rxn = self.num_core_reactions + self.num_edge_reactions

        # Output Array
        cdef np.ndarray[np.float64_t, ndim=1] rates_arr = np.zeros(n_rxn, float)
        cdef double[:] reaction_rates = rates_arr

        # Scratch Arrays
        cdef np.ndarray[np.float64_t, ndim=1] Cg_arr = np.zeros(n_core, float)
        cdef double[:] C_gas = Cg_arr

        cdef np.ndarray[np.float64_t, ndim=1] Cp_arr = np.zeros(n_core, float)
        cdef double[:] C_poly = Cp_arr

        # Declarations
        cdef int r_idx, r0, r1, r2, p0_pool_idx
        cdef int i  # <--- Moved declaration to top
        cdef int p0, p1, p2
        cdef double kf, kb, rf, rr, site, V_rxn
        cdef double V_gas, V_poly
        cdef double n_gas

        # 2. Volumes
        V_poly = self.V_poly
        if self.constant_gas_volume:
            V_gas = self.V_gas0
        else:
            n_gas = 0.0
            for i in range(n_core):
                if self.gas_species_mask[i]:
                    n_gas += y[i]

            V_gas = constants.R * self.T.value_si * n_gas / self.P.value_si if n_gas > 0 else 1.0

        # 3. Concentrations
        for i in range(n_core):
            if self.gas_species_mask[i]:
                C_gas[i] = _softclamp(y[i], self._softclamp_lam) / V_gas
            else:
                C_poly[i] = _softclamp(y[i], self._softclamp_lam) / V_poly

        ir = self.reactant_indices
        ip = self.product_indices

        # 4. Kinetics Loop
        for r_idx in range(n_rxn):
            r0 = ir[r_idx, 0]
            r1 = ir[r_idx, 1]
            r2 = ir[r_idx, 2]

            # Determine Phase for Volume Basis. Edge participants (review
            # round 51: edge-reactant rows are evaluated, not skipped) are
            # judged by prospective_gas_mask -- the only mask covering edge
            # slots; its core prefix is bit-identical to gas_species_mask by
            # rider R1. Raw indexing of the core-sized arrays with an edge
            # slot would be an IndexError.
            is_poly = False
            if r0 != -1 and not (self.gas_species_mask[r0] if r0 < n_core
                                 else self.prospective_gas_mask[r0]):
                is_poly = True
            elif r1 != -1 and not (self.gas_species_mask[r1] if r1 < n_core
                                   else self.prospective_gas_mask[r1]):
                is_poly = True
            elif r2 != -1 and not (self.gas_species_mask[r2] if r2 < n_core
                                   else self.prospective_gas_mask[r2]):
                is_poly = True

            V_rxn = V_poly if is_poly else V_gas

            kf = self.kf[r_idx]
            kb = self.kb[r_idx]

            # Forward Rate. An EDGE reactant has no state in y -> zero
            # concentration, so rf vanishes (simple.pyx :443-450 parity).
            rf = kf
            if r0 >= n_core:
                rf *= 0.0
            elif self.species_to_pool_indices[r0] != -1:
                rf *= 1.0
            elif self.gas_species_mask[r0]:
                rf *= C_gas[r0]
            else:
                rf *= C_poly[r0]

            if r1 != -1:
                if r1 >= n_core:
                    rf *= 0.0
                elif self.species_to_pool_indices[r1] != -1:
                    rf *= 1.0
                elif self.gas_species_mask[r1]:
                    rf *= C_gas[r1]
                else:
                    rf *= C_poly[r1]

            if r2 != -1:
                if r2 >= n_core:
                    rf *= 0.0
                elif self.species_to_pool_indices[r2] != -1:
                    rf *= 1.0
                elif self.gas_species_mask[r2]:
                    rf *= C_gas[r2]
                else:
                    rf *= C_poly[r2]

            # Reverse Rate. An EDGE product has no state in y, so rr is
            # UNCOMPUTABLE and stays 0 (simple.pyx :452-453 parity -- the
            # concentration-availability hole, same as the residual).
            p0 = ip[r_idx, 0]
            p1 = ip[r_idx, 1]
            p2 = ip[r_idx, 2]
            rr = 0.0
            if kb > 0 and p0 < n_core and p1 < n_core and p2 < n_core:
                rr = kb

                if p0 != -1:
                    if self.species_to_pool_indices[p0] != -1:
                        rr *= 1.0
                    elif self.gas_species_mask[p0]:
                        rr *= C_gas[p0]
                    else:
                        rr *= C_poly[p0]

                if p1 != -1:
                    if self.species_to_pool_indices[p1] != -1:
                        rr *= 1.0
                    elif self.gas_species_mask[p1]:
                        rr *= C_gas[p1]
                    else:
                        rr *= C_poly[p1]

                if p2 != -1:
                    if self.species_to_pool_indices[p2] != -1:
                        rr *= 1.0
                    elif self.gas_species_mask[p2]:
                        rr *= C_gas[p2]
                    else:
                        rr *= C_poly[p2]

            # [THE HIJACK] Scale by Moment if Reactant is Proxy (edge
            # reactants are never proxies: pool proxies are core species)
            p0_pool_idx = self.species_to_pool_indices[r0] if (r0 != -1 and r0 < n_core) else -1
            if p0_pool_idx != -1:
                moment_idx = self.polymer_pools[p0_pool_idx].mu_indices[1]

                if self.pool_mu1_indices[p0_pool_idx] != -1:
                    moment_idx = self.pool_mu1_indices[p0_pool_idx]

                if self.is_end_group_reaction[r_idx]:
                    moment_idx = self.pool_mu0_indices[p0_pool_idx]

                site = _softclamp(y[moment_idx], self._softclamp_lam) / V_poly

                # Exhaustion throttle -- parity with the residual's exhaustion
                # throttle (spec 2026-06-10 s5 AMENDMENT): the diagnostic rate
                # must not report a flux the residual refuses to deliver for a
                # mu0-scaled DISCRETE_CHIP past unit exhaustion. Same condition
                # and formula as the residual's section-2 site scaling; here
                # moment_idx is already the resolved mu0 index (end-group
                # branch above), and the mu1 index is resolved with the same
                # pool_mu1_indices override as the default branch.
                if (self.reaction_flux_archetype[r_idx] == FLUX_DISCRETE_CHIP
                        and self.is_end_group_reaction[r_idx]
                        and self.reaction_chip_units[r_idx] > 0):
                    mu1_idx = self.polymer_pools[p0_pool_idx].mu_indices[1]
                    if self.pool_mu1_indices[p0_pool_idx] != -1:
                        mu1_idx = self.pool_mu1_indices[p0_pool_idx]
                    site = min(
                        _softclamp(y[moment_idx], self._softclamp_lam),
                        _softclamp(y[mu1_idx], self._softclamp_lam) / float(self.reaction_chip_units[r_idx]),
                    ) / V_poly
                elif (self.reaction_flux_archetype[r_idx] == FLUX_VOLATILE_EJECTION
                        and self.is_end_group_reaction[r_idx]
                        and self.reaction_src_pool[r_idx] != -1
                        and self.reaction_eject_units_moment[r_idx] > 0.0):
                    # a>0 VE exhaustion throttle -- parity with the
                    # residual's section-2 site scaling (keep in sync). site =
                    # min(mu0, mu1/a); guard a>0 (a<0 grows, no throttle).
                    # Ruling round 20 item C: covers same-pool AND cross-pool
                    # rows (the forward leg debits src either way).
                    mu1_idx = self.polymer_pools[p0_pool_idx].mu_indices[1]
                    if self.pool_mu1_indices[p0_pool_idx] != -1:
                        mu1_idx = self.pool_mu1_indices[p0_pool_idx]
                    site = min(
                        _softclamp(y[moment_idx], self._softclamp_lam),
                        _softclamp(y[mu1_idx], self._softclamp_lam) / float(self.reaction_eject_units_moment[r_idx]),
                    ) / V_poly

                # P1-1 directional bundle throttle, round-27 P1-A tail-only
                # form -- parity with the residual's section-2 scaling
                # (keep in sync): cross-pool VE/MIGRATION legs run each
                # debited direction at the near-exhaustion bundle limiter
                # (tail-only smoothstep; C1 soft-min, cone-margin drain
                # guard): the E-band S_free = w*S_base +
                # (1-w)*softmin_p(S_base, S_cap) of the debited pool
                # (forward debits src) followed by the independent
                # cone-margin drain gate (Q10 <= 0 -> 0 regardless of E);
                # S_eff == S_base exactly for in-cone bulk pools.
                cross_bundle_throttle = False
                if self.reaction_flux_archetype[r_idx] in (
                        FLUX_MIGRATION, FLUX_VOLATILE_EJECTION):
                    src_bt = self.reaction_src_pool[r_idx]
                    dst_bt = self.reaction_dst_pool[r_idx]
                    if src_bt != -1 and dst_bt != -1 and src_bt != dst_bt:
                        cross_bundle_throttle = True
                        site = self._bundle_limited_site(
                            src_bt, y_arr, V_poly,
                            self.is_end_group_reaction[r_idx], site)

                rf *= site
                # Direction-specific source availability -- keep in sync with
                # the residual's section-2 scaling (adjudicated Part C): the
                # reverse leg of a cross-pool exchange debits the dst pool,
                # so rr's site comes from the dst pool's own moments (mu1;
                # mu0 for end-group rows; same pool_mu1_indices override as
                # the forward site). Same-pool and no-dst rows keep the
                # shared (possibly throttled) reactant-pool site.
                dst_pool_idx = self.reaction_dst_pool[r_idx]
                if dst_pool_idx != -1 and dst_pool_idx != p0_pool_idx:
                    moment_idx = self.polymer_pools[dst_pool_idx].mu_indices[1]
                    if self.pool_mu1_indices[dst_pool_idx] != -1:
                        moment_idx = self.pool_mu1_indices[dst_pool_idx]
                    if self.is_end_group_reaction[r_idx]:
                        moment_idx = self.pool_mu0_indices[dst_pool_idx]
                    if (self.reaction_flux_archetype[r_idx] == FLUX_VOLATILE_EJECTION
                            and self.is_end_group_reaction[r_idx]
                            and self.reaction_eject_units_moment[r_idx] < 0.0):
                        # Ruling round 20 item C mirror direction -- parity
                        # with the residual (keep in sync): a<0 reverse leg
                        # sheds |a|/event from the dst pool it debits, so
                        # its availability is min(mu0, mu1/|a|).
                        dst_mu1_idx = self.polymer_pools[dst_pool_idx].mu_indices[1]
                        if self.pool_mu1_indices[dst_pool_idx] != -1:
                            dst_mu1_idx = self.pool_mu1_indices[dst_pool_idx]
                        site_rev = min(
                            _softclamp(y[self.pool_mu0_indices[dst_pool_idx]], self._softclamp_lam),
                            _softclamp(y[dst_mu1_idx], self._softclamp_lam) / (-float(self.reaction_eject_units_moment[r_idx])),
                        ) / V_poly
                    else:
                        site_rev = _softclamp(y[moment_idx], self._softclamp_lam) / V_poly
                    # P1-1 mirror direction, round-27 P1-A tail-only form
                    # -- parity with the residual (keep in sync): reverse
                    # leg debits the dst pool's full bundle, same tail-only
                    # smoothstep (C1 soft-min, cone-margin drain guard)
                    # S_eff from dst's own moments, floor distance and cone
                    # margin.
                    if cross_bundle_throttle:
                        site_rev = self._bundle_limited_site(
                            dst_pool_idx, y_arr, V_poly,
                            self.is_end_group_reaction[r_idx], site_rev)
                    rr *= site_rev
                else:
                    rr *= site
            elif rr != 0.0:
                # PRODUCT-side proxy mirror (review round 51; keep in sync
                # with the residual's product-side scaling): a reverse-stored
                # association row carries the proxy on the PRODUCT side, so
                # rr's 1.0 placeholder is replaced by the pool site density.
                # rr != 0.0 implies all products are core (the rr hole above),
                # so species_to_pool_indices indexing is safe. rf untouched:
                # no reactant is a proxy on this branch (elif).
                prod_pool_idx = -1
                if p0 != -1 and self.species_to_pool_indices[p0] != -1:
                    prod_pool_idx = self.species_to_pool_indices[p0]
                elif p1 != -1 and self.species_to_pool_indices[p1] != -1:
                    prod_pool_idx = self.species_to_pool_indices[p1]
                elif p2 != -1 and self.species_to_pool_indices[p2] != -1:
                    prod_pool_idx = self.species_to_pool_indices[p2]
                if prod_pool_idx != -1:
                    moment_idx = self.polymer_pools[prod_pool_idx].mu_indices[1]
                    if self.pool_mu1_indices[prod_pool_idx] != -1:
                        moment_idx = self.pool_mu1_indices[prod_pool_idx]
                    if self.is_end_group_reaction[r_idx]:
                        moment_idx = self.pool_mu0_indices[prod_pool_idx]
                    rr *= _softclamp(y[moment_idx], self._softclamp_lam) / V_poly

            reaction_rates[r_idx] = (rf - rr)

        return rates_arr

    def get_edge_reaction_rates(self, core_y, edge_rates):
        """
        Calculate ONLY edge reaction rates using the hijacked logic.
        """
        # Call full calculation
        full_rates = self.get_reaction_rates(core_y)

        # Copy to output buffer
        cdef int n_core_rxn = self.num_core_reactions
        cdef int n_edge_rxn = self.num_edge_reactions
        cdef double[:] edge_view = edge_rates  # Cast input to view
        cdef double[:] full_view = full_rates

        cdef int i
        for i in range(n_edge_rxn):
            edge_view[i] = full_view[n_core_rxn + i]
