#!/usr/bin/env bash
# Verifier for I-060: the hybrid handshake, and the two realizability budgets
# that bound it.
#
#   bash test/rmgpy/solver/i060_verify.sh
#
# Exits 0 only if every check passes. Every number below is RECOMPUTED here;
# nothing is restated from a report. Nothing below relaxes any solver detector,
# floor, tolerance or test.
#
# TWO ARMS, and both are pinned to their source rather than to a timestamp:
#   FIXED   = the tree this script is shipped in. Derived from the script's own
#             location, never hardcoded -- a verifier that names its birth
#             worktree keeps measuring that worktree after the change merges,
#             and dies once the worktree is deleted. It is REBUILT here, and
#             the rebuild must not change the .so hash: that, and not an mtime,
#             is the proof the binary matches the source. (The build is
#             byte-reproducible on this toolchain; if it ever stops being, this
#             check fails loudly rather than passing quietly.)
#   UNFIXED = /home/alon/Code/RMG-Py-i065-mainline @ 9e7e0c4d5, EXACTLY this
#             branch's base with none of the fix. Pinned to a sha ON PURPOSE (a
#             "before" build must not follow mainline), asserted clean, and READ
#             ONLY -- so it is not rebuilt; its binary is pinned by hash below.
#
# Every python invocation runs from /tmp: `python -m pytest` and `python -c` put
# the cwd FIRST on sys.path, ahead of PYTHONPATH, so running from the worktree
# would silently import the FIXED extension into the unfixed arm and make every
# before/after comparison vacuous by construction.
set -uo pipefail

FIXED=${I060_FIXED_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd -P)}
UNFIXED=${I060_UNFIXED_DIR:-/home/alon/Code/RMG-Py-i065-mainline}   # READ ONLY
UNFIXED_SHA=9e7e0c4d5e9572b2a67c1b1d8bb8880287a793e1
UNFIXED_SO_SHA=435bc01c205f3369099627e790f9dab326463adda482224bec5107a5b9d1c051
XS=3                      # the probe fixture's cutoff (i060_tail_probe.XS)
T=${TMPDIR:-/tmp}/i060_verify.$$
mkdir -p "$T"
FAIL=0

step()    { printf '\n=== %s ===\n' "$*"; }
verdict() { if [ "$1" -eq 0 ]; then printf 'PASS: %s\n' "$2"
            else printf 'FAIL: %s (rc=%s)\n' "$2" "$1"; FAIL=1; fi; }
export PYTHONDONTWRITEBYTECODE=1

# shellcheck disable=SC1091
source /home/alon/anaconda3/etc/profile.d/conda.sh
conda activate rmg_env

################################################################################
step "1. the extension under test is BUILT FROM the source under test"
################################################################################
printf '  FIXED   tree: %s\n' "$FIXED"
printf '  UNFIXED tree: %s\n' "$UNFIXED"
SO_REL=rmgpy/solver/polymer.cpython-39-x86_64-linux-gnu.so
SO_F="$FIXED/$SO_REL"
SO_U="$UNFIXED/$SO_REL"
BEFORE=""
[ -f "$SO_F" ] && BEFORE=$(sha256sum "$SO_F" | cut -d' ' -f1)
printf '  forcing a rebuild of the tree under test (~30 s)...\n'
touch "$FIXED/rmgpy/solver/polymer.pyx"
( cd "$FIXED" && python setup.py build_ext --inplace -j 8 ) \
  >"$T/build.log" 2>"$T/build.err"
BRC=$?
printf '  build rc=%s (log %s)\n' "$BRC" "$T/build.log"
[ "$BRC" -eq 0 ] || { printf 'FAIL: the tree under test does not build\n'; tail -5 "$T/build.err"; FAIL=1; }
AFTER=$(sha256sum "$SO_F" | cut -d' ' -f1)
printf '  fixed   .so: %s\n' "$SO_F"
printf '    sha256 before rebuild: %s\n' "${BEFORE:-<absent>}"
printf '    sha256 after  rebuild: %s\n' "$AFTER"
if [ -z "$BEFORE" ]; then
  printf 'NOTE: no prior binary to compare -- built fresh, so it matches by construction\n'
elif [ "$BEFORE" = "$AFTER" ]; then
  verdict 0 "the binary under test is the compile of the source under test (rebuild was a no-op)"
else
  verdict 1 "the binary under test was STALE -- it has just been rebuilt and every number below would have been read off the wrong .so"
fi
# the baseline: pinned commit, clean, and pinned BINARY (it is read-only, so it
# is not rebuilt; the hash is what stands in for that).
H=$(git -C "$UNFIXED" rev-parse HEAD)
TRACKED=$(git -C "$UNFIXED" status --porcelain | grep -v '^??' || true)
UNTRACKED=$(git -C "$UNFIXED" status --porcelain | grep '^??' || true)
U_SHA=$(sha256sum "$SO_U" | cut -d' ' -f1)
printf '  unfixed HEAD=%s\n  unfixed .so sha256=%s\n' "$H" "$U_SHA"
[ "$H" = "$UNFIXED_SHA" ]        || { printf 'FAIL: unfixed arm is not at %s\n' "$UNFIXED_SHA"; FAIL=1; }
[ "$U_SHA" = "$UNFIXED_SO_SHA" ] || { printf 'FAIL: unfixed .so is not the pinned baseline binary\n'; FAIL=1; }
# The baseline is a baseline iff its TRACKED source is the pinned commit and
# its BINARY is the pinned hash. Untracked files cannot affect either, and step
# 4 below writes some itself -- running the suite against this tree drops
# chemkin output next to the package, not in the cwd. They are listed rather
# than ignored, but they are not a reason to call the baseline unusable; the
# .so hash pin above is the load-bearing check and is far stronger than a
# dirty-tree test.
if [ -n "$TRACKED" ]; then
  printf 'FAIL: unfixed arm has MODIFIED TRACKED files -- not a usable baseline:\n'
  printf '%s\n' "$TRACKED" | sed 's/^/    /'
  FAIL=1
else
  printf '  unfixed tracked tree: clean (0 modified)\n'
fi
if [ -n "$UNTRACKED" ]; then
  printf '  unfixed untracked scratch (test output; source and binary unaffected):\n'
  printf '%s\n' "$UNTRACKED" | sed 's/^/    /'
fi
if [ "$AFTER" != "$U_SHA" ]; then verdict 0 "the two arms are genuinely different binaries"
else verdict 1 "the two arms are the SAME binary -- every comparison below is vacuous"; fi
for arm in fixed unfixed; do
  case $arm in fixed) P=$FIXED; W=$SO_F ;; unfixed) P=$UNFIXED; W=$SO_U ;; esac
  M=$(cd /tmp && PYTHONPATH="$P" python -c 'import rmgpy.solver.polymer as p; print(p.__file__)')
  printf '  %-8s imports %s\n' "$arm" "$M"
  [ "$M" = "$W" ] || { printf 'FAIL: %s arm imported the wrong build\n' "$arm"; FAIL=1; }
done

################################################################################
step "2. the probe, on both arms, state-matched"
################################################################################
for arm in unfixed fixed; do
  case $arm in fixed) P=$FIXED ;; unfixed) P=$UNFIXED ;; esac
  (cd /tmp && PYTHONPATH="$P" python \
     "$FIXED/test/rmgpy/solver/i060_tail_probe.py" "$T/probe.$arm.json") \
     >"$T/probe.$arm.out" 2>"$T/probe.$arm.err"
  [ -s "$T/probe.$arm.json" ] || { printf 'FAIL: %s probe produced nothing\n' "$arm"; FAIL=1; }
done
# ...and the UNFIXED build replayed on the FIXED build's own accepted states,
# so the per-step comparison is one build against another AT THE SAME STATE,
# not two different trajectories talked about as if they were one.
python -c "
import json
d=json.load(open('$T/probe.fixed.json'))
json.dump(d['traj']['states'], open('$T/states.json','w'))"
(cd /tmp && PYTHONPATH="$UNFIXED" python \
   "$FIXED/test/rmgpy/solver/i060_tail_probe.py" "$T/probe.replay.json" \
   "$T/states.json") >"$T/probe.replay.out" 2>"$T/probe.replay.err"

python - "$T/probe.unfixed.json" "$T/probe.fixed.json" "$T/probe.replay.json" "$XS" <<'PY'
import json, sys
u = json.load(open(sys.argv[1]))
f = json.load(open(sys.argv[2]))
r = json.load(open(sys.argv[3]))
XS = int(sys.argv[4])
K = u['k_unzip']                       # k_chain_handshake on this fixture
hx = float.fromhex
ok = True
print(f"  unfixed build: {u['build']}")
print(f"  fixed   build: {f['build']}")
if u['build'] == f['build']:
    print("  FAIL: both arms loaded the SAME extension"); sys.exit(1)

# ---------------------------------------------------------------- round 2 ---
print("\n  [A] THE BLOCKER: mu0*mu2 >= mu1^2. A handshake event removes one")
print("      chain at exactly n = xs, so with moments about xs it moves m_0")
print("      alone and dQ/dt = -F*m_2 <= 0 unconditionally. The budget that")
print("      keeps Q >= 0 is c_xs <= Q/m_2; capping there is exactly")
print("      -dQ/dt <= k*Q. Swept over the corner adversarial for Q (LOW PDI,")
print("      where Q = mu0^2*Var is small) -- a different corner from the one")
print("      adversarial for mu1 - mu0.")
for name, d in (("unfixed", u), ("fixed", f)):
    q = d['q_field']
    over = [x for x in q if x['N'] > x['B_Q'] * (1.0 + 1e-12)]
    above = [x for x in over if x['mean'] > XS + 1e-9]
    worst = max((x['F'] * x['m2'] / x['Q']) for x in q if x['Q'] > 0.0)
    print(f"    {name:<8} states {len(q):<5} over the Q budget {len(over):<5} "
          f"(of those, above the old cutoff: {len(above)})")
    print(f"             max drain -dQ/dt per unit Q = {worst:.6e}   "
          f"(the budget permits k = {K:g})")
    d['_over'], d['_above'], d['_worst'] = len(over), len(above), worst
if u['_over'] <= 0:
    ok = False; print("    FAIL: the defect did not reproduce on the unfixed arm")
if u['_above'] != u['_over']:
    print(f"    NOTE: {u['_over'] - u['_above']} of the unfixed over-budget "
          f"states are at or below the cutoff")
else:
    print(f"    ALL {u['_over']} unfixed over-budget states are ABOVE the old "
          f"cutoff, where the removed boolean was TRUE -- so this is a "
          f"PRE-EXISTING defect, not one round 1 introduced")
if f['_over'] != 0:
    ok = False; print("    FAIL: the fixed arm still exceeds the Q budget")
if not (f['_worst'] <= K * (1 + 1e-9)):
    ok = False; print("    FAIL: the fixed arm drains Q faster than k*Q")
else:
    print(f"    fixed arm's worst drain is {f['_worst']:.6e} against k = {K:g}: "
          f"the bound is saturated and never exceeded, so Q = 0 is an "
          f"invariant set")
print(f"    unfixed/fixed worst-drain ratio: {u['_worst'] / f['_worst']:.4g}x")

print("\n  [B] the same thing on a trajectory, EVERY step, both arms.")
print("      Start state: the worst realizable state a seeded search named,")
print("      mean DP = xs + 1.5 -- ABOVE the cutoff -- with k_scission = 0 so")
print("      the handshake is the only channel that can drain Q.")
for name, d in (("unfixed", u), ("fixed", f)):
    t = d['q_traj']['rows']
    neg = [x for x in t if x['Q'] < 0.0]
    print(f"    {name:<8} steps {len(t)}  Q(0) = {t[0]['Q']:.6e}  "
          f"min Q = {min(x['Q'] for x in t):.6e}  steps with Q < 0: {len(neg)}")
    d['_neg'] = len(neg)
if u['_neg'] <= 0:
    ok = False
    print("    FAIL: the unfixed arm did not leave the variance cone -- the "
          "trajectory does not exercise the defect")
if f['_neg'] != 0:
    ok = False
    print("    FAIL: the fixed arm leaves the variance cone")
tq = f['q_traj']['rows']
bad = [x for x in tq if x['mu1'] - x['mu0'] < 0.0]
print(f"    fixed arm also holds mu1 >= mu0 on {len(tq) - len(bad)}/{len(tq)} "
      f"steps, min {min(x['mu1'] - x['mu0'] for x in tq):.6e}")
if bad:
    ok = False; print("    FAIL: fixed arm violated mu1 >= mu0 here")

# ---------------------------------------------------------------- round 1 ---
print("\n  [C] the round-1 defect still fixed: the outlet no longer switches")
print("      off at the cutoff, and the residual no longer jumps doing it.")
print(f"  {'tail mean DP':>14}{'F unfixed':>15}{'F fixed':>15}"
      f"{'d(mu1-(xs+1)mu0)/dt unf':>26}{'  fixed':>14}")
for a, b in zip(u['field'], f['field']):
    print(f"  {a['mean']:>14.6f}{hx(a['handshake']):>15.6e}"
          f"{hx(b['handshake']):>15.6e}{hx(a['de_tail']):>26.6e}"
          f"{hx(b['de_tail']):>14.6e}")
hi = [a for a in u['field'] if abs(a['mean'] - (XS + 1e-6)) < 1e-12][0]
at = [a for a in u['field'] if a['mean'] == float(XS)][0]
bhi = [a for a in f['field'] if abs(a['mean'] - (XS + 1e-6)) < 1e-12][0]
bat = [a for a in f['field'] if a['mean'] == float(XS)][0]
jhi, jat = hx(hi['handshake']), hx(at['handshake'])
khi, kat = hx(bhi['handshake']), hx(bat['handshake'])
print(f"\n  unfixed: F(xs+1e-6) = {jhi:.6e} -> F(xs) = {jat:.6e}   "
      f"relative step {abs(jhi - jat) / jhi:.6f}")
print(f"  fixed  : F(xs+1e-6) = {khi:.6e} -> F(xs) = {kat:.6e}   "
      f"relative step {abs(khi - kat) / khi:.3e}")
if not (jhi > 0.0 and jat == 0.0 and abs(jhi - jat) / jhi > 0.5):
    ok = False; print("  FAIL: the unfixed jump did not reproduce")
if abs(khi - kat) / khi > 1e-4:
    ok = False; print("  FAIL: the fixed arm still jumps across the cutoff")
print(f"  unfixed: dmu0(xs+1e-6) = {hx(hi['dmu0']):.6e} -> dmu0(xs) = "
      f"{hx(at['dmu0']):.6e}  (sign flip across the surface)")
if not (hx(hi['dmu0']) < 0.0 < hx(at['dmu0'])):
    ok = False; print("  FAIL: expected the chain-count derivative to flip sign")
below = [b for b in f['field'] if 1.0 < b['mean'] < XS]
if not all(hx(b['handshake']) > 0.0 for b in below):
    ok = False; print("  FAIL: fixed arm still has zero flux below the cutoff")
else:
    print(f"  fixed  : flux strictly positive at all {len(below)} swept states "
          f"with 1 < mean < xs (unfixed: exactly 0.0 at every one)")
edge = [b for b in f['field'] if b['mean'] == 1.0]
if not all(hx(b['handshake']) == 0.0 for b in edge):
    ok = False; print("  FAIL: flux is nonzero at the cone edge mu1 == mu0")
else:
    print("  fixed  : EXACTLY 0.0 at mu1 == mu0, so that edge is invariant too")

print("\n  [D] two-arm, per-step, ON THE SAME STATES (unfixed replayed on the")
print("      fixed trajectory's own accepted states).")
tf = f['traj']['rows'][::10]
assert len(tf) == len(r['replay']), (len(tf), len(r['replay']))
n_worse = n_shut = 0
worst = 0.0
for a, b in zip(tf, r['replay']):
    df = a['dmu1'] - (XS + 1) * a['dmu0']
    du = hx(b['dmu1']) - (XS + 1) * hx(b['dmu0'])
    Ffix, Funf = a['handshake'], hx(b['handshake'])
    worst = max(worst, abs((df - du) - (Ffix - Funf)))
    if Funf == 0.0 and Ffix > 0.0:
        n_shut += 1
    if df < du:
        n_worse += 1
print(f"  {len(tf)} accepted states replayed. On every one,")
print( "    d(mu1-(xs+1)mu0)/dt_fixed - _unfixed == F_fixed - F_unfixed")
print(f"    exactly; largest residue {worst:.3e} (so the whole difference")
print( "    between the builds is the handshake flux and nothing else)")
if worst > 1e-12:
    ok = False; print("  FAIL: the difference is not accounted for by the flux alone")
print(f"  the unfixed outlet is shut (F == 0.0 exactly) on {n_shut} of those "
      f"states; drift worse on {n_worse}")
if n_worse:
    ok = False; print("  FAIL: the fix made the tail-support drift worse somewhere")

print("\n  [E] WHAT CHANGED where the removed boolean was TRUE. The round-1")
print("      claim was 'bit-for-bit inert above the cutoff'. That is no longer")
print("      the whole truth and is not asserted as such: the Q budget binds")
print("      above the cutoff too, deliberately, because that is where the")
print("      pre-existing over-removal lives. What IS asserted, per state:")
print("        - where NEITHER budget binds: bit-for-bit identical;")
print("        - where one binds: the flux STRICTLY FALLS, never rises.")
n_inert = n_bind = n_bad = 0
for key in ("field", "field_mono", "random"):
    for a, b in zip(u[key], f[key]):
        if a['mean'] <= XS + 1e-9:
            continue                      # not the old boolean's region
        mu0, mu1, mu2 = a['mu0'], a['mu1'], a['mu2']
        m2 = mu2 - 2.0 * XS * mu1 + XS * XS * mu0
        q = mu0 * mu2 - mu1 * mu1
        b_exc = (mu1 - mu0) / (XS - 1) if XS > 1 else float('inf')
        b_q = q / m2 if m2 > 0.0 else float('inf')
        n_att = hx(a['handshake']) / K
        binds = n_att > min(b_exc, b_q) * (1.0 + 1e-12)
        same = all(a[fld] == b[fld]
                   for fld in ("handshake", "dmu0", "dmu1", "dmu2", "gas"))
        if binds:
            n_bind += 1
            if hx(b['handshake']) > hx(a['handshake']) + 1e-18:
                n_bad += 1
                ok = False
                print(f"    FAIL {key} mean={a['mean']!r}: flux ROSE "
                      f"{hx(a['handshake'])!r} -> {hx(b['handshake'])!r}")
        else:
            n_inert += 1
            if not same:
                n_bad += 1
                ok = False
                print(f"    FAIL {key} mean={a['mean']!r}: no budget binds yet "
                      f"the residual differs")
print(f"    states above the cutoff: {n_inert + n_bind}")
print(f"      no budget binds: {n_inert}  -- all bit-for-bit identical on all "
      f"five residual components")
print(f"      a budget binds : {n_bind}  -- flux strictly reduced on every one")
print(f"    violations: {n_bad}")
cone = [a for a in u['random'] if a['tag'] == 'cone']
print(f"    ({len(cone)} of these states are randomized, seed 20260904: mu0 "
      f"over 8 decades, mean DP from xs+1e-9 to xs+1e3, PDI in (1, 6])")

print("\n  [F] THE INVARIANTS, at EVERY accepted step of the scission-fed")
print("      trajectory -- not at the endpoints.")
rows = f['traj']['rows']
u_rows = u['traj']['rows']
held0 = rows[0]['mu1']
bad_cone = [x for x in rows if x['mu1'] - x['mu0'] < 0.0]
bad_var = [x for x in rows if x['mu0'] * x['mu2'] - x['mu1'] ** 2 < 0.0]
bad_mass = [x for x in rows
            if abs(XS * x['pxs'] + x['gas'] + x['mu1'] - held0) > 1e-9 * held0]
print(f"    steps checked: {len(rows)}")
print(f"    mu1 >= mu0        : violated on {len(bad_cone)}; "
      f"min {min(x['mu1'] - x['mu0'] for x in rows):.6e} "
      f"(unfixed {min(x['mu1'] - x['mu0'] for x in u_rows):.6e})")
print(f"    mu0*mu2 >= mu1^2  : violated on {len(bad_var)}; "
      f"min {min(x['mu0'] * x['mu2'] - x['mu1'] ** 2 for x in rows):.6e} "
      f"(unfixed {min(x['mu0'] * x['mu2'] - x['mu1'] ** 2 for x in u_rows):.6e})")
print(f"    repeat-unit ledger: xs*explicit + gas + mu1 == mu1(0) = {held0:g} "
      f"to 1e-9 rel on {len(rows) - len(bad_mass)}/{len(rows)} steps")
if bad_cone or bad_var or bad_mass:
    ok = False; print("    FAIL: an invariant above is violated on the fixed arm")
mf = min(x['e_tail'] for x in rows)
mu_ = min(x['e_tail'] for x in u_rows)
print(f"    FINDING, not a pass/fail: the tail's own SUPPORT invariant "
      f"mu1 >= (xs+1)*mu0 is violated on BOTH arms (min fixed {mf:.6e}, "
      f"unfixed {mu_:.6e}). The scission kernel manufactures sub-cutoff chains "
      f"inside the tail and nothing routes them out; the handshake can only "
      f"slow it. Pre-existing, and better fixed than unfixed by {mf - mu_:.6e}.")

print("\n  the consequence, end to end (same fixture, same rate constants):")
for nm, rr in (("unfixed", u_rows), ("fixed", rows)):
    last = rr[-1]
    print(f"    {nm:<8} explicit DP={XS}: {last['pxs']:.6f} mol   "
          f"gas monomer: {last['gas']:.6f} mol   "
          f"ledger {XS * last['pxs'] + last['gas'] + last['mu1']:.9f}")
if not rows[-1]['pxs'] > 1.5 * u_rows[-1]['pxs']:
    ok = False
    print("  FAIL: the fixed arm did not route materially more mass into the "
          "explicit ladder")
sys.exit(0 if ok else 1)
PY
verdict "$?" "both cone halves: defects reproduce unfixed, are gone fixed, budgets bound every state"

################################################################################
step "3. the added tests: the EXACT failing set on the unfixed build, by name"
################################################################################
# Not `len(red) >= 1`: that is satisfied by an import error or a missing symbol
# on the unfixed arm, i.e. by the code being ABSENT rather than by the behaviour
# being wrong. Each expected failure is named, and each must fail for its own
# reason (regex over the assertion message), so a collection error cannot pass
# for a reproduction.
CLS=test/rmgpy/solver/solverPolymerTest.py::TestTailExhaustionHandshake
for arm in unfixed fixed; do
  case $arm in fixed) P=$FIXED ;; unfixed) P=$UNFIXED ;; esac
  (cd /tmp && PYTHONPATH="$P" python -m pytest "$FIXED/$CLS" \
     --no-cov -q --tb=long -p no:cacheprovider \
     --junitxml="$T/new.$arm.xml" >"$T/new.$arm.log" 2>&1)
  printf '  %-8s %s\n' "$arm" "$(grep -E 'passed|failed' "$T/new.$arm.log" | tail -1)"
done
python - "$T/new.unfixed.xml" "$T/new.fixed.xml" <<'PY'
import re, sys, xml.etree.ElementTree as ET

# name -> (must fail unfixed, regex its failure message must match)
EXPECTED = {
    "test_handshake_survives_the_boundary_crossing":
        r"outlet into the explicit ladder is shut",
    "test_handshake_flux_is_continuous_across_the_cutoff":
        r"residual jumps across the cutoff",
    "test_flux_never_rises_where_the_old_boolean_was_true":
        r"no state above the cutoff had a budget bind",
    "test_the_cutoff_slack_is_a_declared_change":
        r"the old boolean's own boundary",
    "test_scission_fed_pool_keeps_its_outlet_through_exhaustion":
        r"outlet shut and the residue left as gas",
    "test_handshake_q_drain_is_bounded_by_k_times_q":
        r"drains Q at up to .* against the k = .* the budget permits",
    "test_handshake_cannot_drive_the_variance_cone_negative":
        r"mu0\*mu2 - mu1\^2 went NEGATIVE on \d+ of \d+ steps",
    "test_degenerate_cutoffs_still_respect_the_variance_cone":
        r"drains Q at .* per unit Q against the k",
}


def load(p):
    out = {}
    for tc in ET.parse(p).getroot().iter("testcase"):
        state, msg = "passed", ""
        for ch in tc:
            if ch.tag in ("failure", "error"):
                state = ch.tag
                msg = (ch.get("message") or "") + (ch.text or "")
            elif ch.tag == "skipped":
                state = "skipped"
        out[tc.get("name")] = (state, msg)
    return out


pre, post = load(sys.argv[1]), load(sys.argv[2])
ok = True
if not pre:
    print("  FAIL: the unfixed arm collected NO tests -- an absent-code result, "
          "not a behavioural one")
    sys.exit(1)
red = {k for k, (s, _) in pre.items() if s != "passed"}
print(f"  unfixed: {len(pre)} collected, {len(red)} failing")
missing = sorted(set(EXPECTED) - red)
extra = sorted(red - set(EXPECTED))
if missing:
    ok = False
    print(f"  FAIL: expected to fail on the unfixed build but did not: {missing}")
if extra:
    ok = False
    print(f"  FAIL: failed on the unfixed build but was not expected to: {extra}")
print("  each expected failure, and the reason it actually failed for:")
for name in sorted(EXPECTED):
    state, msg = pre.get(name, ("<absent>", ""))
    pat = EXPECTED[name]
    hit = bool(re.search(pat, msg, re.S))
    m = " ".join(msg.split())
    m = (m[-110:] if len(m) > 110 else m)
    print(f"    [{'ok ' if hit and state == 'failure' else 'BAD'}] {name}")
    print(f"          state={state}  matches /{pat}/: {hit}")
    print(f"          ...{m}")
    if not (hit and state == "failure"):
        ok = False
        print(f"    FAIL: {name} did not fail for its stated reason "
              f"(an import or collection error would land here)")
grn = sorted(k for k, (s, _) in post.items() if s != "passed")
print(f"  fixed: {len(post)} collected, still failing: {grn or 'none'}")
if grn:
    ok = False
print("  green on BOTH arms (property pins, not discriminators, and reported "
      "as such rather than counted as evidence):")
for k in sorted(k for k in pre if pre[k][0] == "passed" and post[k][0] == "passed"):
    print(f"    {k}")
sys.exit(0 if ok else 1)
PY
verdict "$?" "the exact expected set fails before, each for its own reason, and all pass after"

################################################################################
step "4. the full existing suite, both arms, pass/fail sets named"
################################################################################
# --no-cov: the repo config enables --cov by default and the suite does not
#           finish in reasonable time with it. The coverage CONFIG is not edited.
# -k "not <SLOW>": the poly_102 crash-window replay does not terminate in
#           reasonable time on EITHER build. Deselected, not skipped, and
#           disclosed here -- it is non-terminating, not failing.
SLOW=test_v2_v3_crash_window_twin_and_episode_invariance
SUITE=(test/rmgpy/solver/solverPolymerTest.py
       test/rmgpy/solver/solverPolymerJacobianTest.py
       test/rmgpy/solver/solverPolymerConduitTest.py
       test/rmgpy/polymerTest.py
       test/rmgpy/rmg/inputTest.py)
for arm in unfixed fixed; do
  case $arm in fixed) P=$FIXED ;; unfixed) P=$UNFIXED ;; esac
  PATHS=(); for s in "${SUITE[@]}"; do PATHS+=("$FIXED/$s"); done
  (cd /tmp && PYTHONPATH="$P" python -m pytest "${PATHS[@]}" \
     --no-cov -q --tb=no -p no:cacheprovider -k "not $SLOW" \
     --junitxml="$T/suite.$arm.xml" >"$T/suite.$arm.log" 2>&1)
  printf '  %-8s %s\n' "$arm" "$(tail -1 "$T/suite.$arm.log")"
done
python - "$T/suite.unfixed.xml" "$T/suite.fixed.xml" <<'PY'
import sys, xml.etree.ElementTree as ET
def load(p):
    r = {}
    for tc in ET.parse(p).getroot().iter("testcase"):
        k = "passed"
        for ch in tc:
            if ch.tag in ("failure", "error"): k = ch.tag
            elif ch.tag == "skipped": k = "skipped"
        r[f'{tc.get("classname")}::{tc.get("name")}'] = k
    return r
pre, post = load(sys.argv[1]), load(sys.argv[2])
c = lambda r, k: sum(1 for v in r.values() if v == k)
for n, r in (("unfixed", pre), ("fixed", post)):
    print(f"  {n:<8} total={len(r):<5} passed={c(r,'passed'):<5} "
          f"failed={c(r,'failure') + c(r,'error'):<4} skipped={c(r,'skipped')}")
ok = True
if len(post) < 800:
    ok = False; print(f"  FAIL: only {len(post)} tests collected; expected >= 800")
prf = {k for k, v in pre.items() if v in ("failure", "error")}
fin = {k for k, v in post.items() if v in ("failure", "error")}
own = lambda s: {k for k in s if "TestTailExhaustionHandshake" in k}
print(f"  failing on the UNFIXED arm, EXCLUDING the added class "
      f"({len(prf - own(prf))}) -- all pre-existing:")
for k in sorted(prf - own(prf)): print(f"    {k}")
print(f"  failing on the FIXED arm ({len(fin)}):")
for k in sorted(fin): print(f"    {k}")
changed = sorted(k for k in set(pre) | set(post) if pre.get(k) != post.get(k))
print(f"  tests whose state CHANGED between the arms: {len(changed)}")
for k in changed: print(f"    {k}: {pre.get(k)} -> {post.get(k)}")
expected = own(set(changed))
if len(changed) != len(expected):
    ok = False
    print(f"  FAIL: {len(changed) - len(expected)} state change(s) OUTSIDE the "
          f"added class")
else:
    print(f"  every one of the {len(changed)} change(s) is in "
          f"TestTailExhaustionHandshake, and the count matches step 3")
if c(post, "passed") < c(pre, "passed"):
    ok = False
    print(f"  FAIL: passed count fell from {c(pre,'passed')} to {c(post,'passed')}")
sys.exit(0 if ok else 1)
PY
verdict "$?" "no regression; every state change is the added class and the count matches"

################################################################################
step "5. nothing out of scope changed"
################################################################################
printf '  files changed on this branch vs its base:\n'
git -C "$FIXED" diff --stat "$UNFIXED_SHA" -- | sed 's/^/    /'
printf '  worktree status:\n'
git -C "$FIXED" status --porcelain | sed 's/^/    /'
printf '    (docs/contracts and build artifacts are git-ignored)\n'
printf '  baseline tree left untouched: %s modified/untracked entries\n' \
  "$(git -C "$UNFIXED" status --porcelain | wc -l)"

################################################################################
printf '\n================================================================\n'
if [ "$FAIL" -eq 0 ]; then printf 'I-060 VERIFIER: ALL CHECKS PASS\n'
else printf 'I-060 VERIFIER: FAILURES ABOVE\n'; fi
printf '================================================================\n'
rm -rf "$T"
exit "$FAIL"
