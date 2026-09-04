#!/usr/bin/env bash
# Verifier for I-060: the hybrid handshake across the tail-exhaustion boundary.
#
#   bash test/rmgpy/solver/i060_verify.sh
#
# Exits 0 only if every check passes. Every number below is RECOMPUTED here;
# nothing is restated from a report. Nothing below relaxes any solver detector,
# floor, tolerance or test.
#
# TWO ARMS, and they are provably different builds:
#   FIXED   = the tree this script is shipped in, rebuilt in place. Derived
#             from the script's own location, never hardcoded -- a verifier
#             that names its birth worktree keeps measuring that worktree
#             after the change merges, and dies once the worktree is deleted.
#   UNFIXED = /home/alon/Code/RMG-Py-i065-mainline @ 9e7e0c4d5, which is
#             EXACTLY this branch's base with none of the fix. Pinned to a sha
#             ON PURPOSE (a "before" build must not follow mainline), asserted
#             clean, and READ ONLY throughout.
#
# Every python invocation runs from /tmp: `python -m pytest` and `python -c`
# put the cwd FIRST on sys.path, ahead of PYTHONPATH, so running from the
# worktree would silently import the FIXED extension into the unfixed arm and
# make every before/after comparison vacuous by construction.
set -uo pipefail

FIXED=${I060_FIXED_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd -P)}
UNFIXED=${I060_UNFIXED_DIR:-/home/alon/Code/RMG-Py-i065-mainline}   # READ ONLY
UNFIXED_SHA=9e7e0c4d5e9572b2a67c1b1d8bb8880287a793e1
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
step "1. the extension under test is the one built in THIS tree"
################################################################################
printf '  FIXED   tree: %s\n' "$FIXED"
printf '  UNFIXED tree: %s\n' "$UNFIXED"
SO_F=$(ls "$FIXED"/rmgpy/solver/polymer.cpython-*.so 2>/dev/null | head -1)
SO_U=$(ls "$UNFIXED"/rmgpy/solver/polymer.cpython-*.so 2>/dev/null | head -1)
[ -n "$SO_F" ] || { printf 'FAIL: no built polymer extension in %s\n' "$FIXED"; exit 1; }
[ -n "$SO_U" ] || { printf 'FAIL: no built polymer extension in %s\n' "$UNFIXED"; exit 1; }
printf '  fixed   .so: %s\n    sha256 %s\n' "$SO_F" "$(sha256sum "$SO_F" | cut -d' ' -f1)"
printf '  unfixed .so: %s\n    sha256 %s\n' "$SO_U" "$(sha256sum "$SO_U" | cut -d' ' -f1)"
case "$SO_F" in "$FIXED"/*) ;; *) printf 'FAIL: fixed .so is outside the tree under test\n'; FAIL=1 ;; esac
# fresh, not stale
PYX="$FIXED/rmgpy/solver/polymer.pyx"
printf '  pyx mtime: %s\n  so  mtime: %s\n' "$(date -r "$PYX" -Is)" "$(date -r "$SO_F" -Is)"
if [ "$SO_F" -nt "$PYX" ]; then verdict 0 "compiled extension is newer than its .pyx"
else verdict 1 "compiled extension is STALE relative to the .pyx"; fi
# the baseline is the pinned commit, clean, and a genuinely different binary
H=$(git -C "$UNFIXED" rev-parse HEAD); D=$(git -C "$UNFIXED" status --porcelain)
printf '  unfixed HEAD=%s dirty=%s\n' "$H" "${D:-<clean>}"
[ "$H" = "$UNFIXED_SHA" ] || { printf 'FAIL: unfixed arm is not at %s\n' "$UNFIXED_SHA"; FAIL=1; }
[ -z "$D" ] || { printf 'FAIL: unfixed arm is dirty -- not a usable baseline\n'; FAIL=1; }
if [ "$(sha256sum "$SO_F" | cut -c1-16)" != "$(sha256sum "$SO_U" | cut -c1-16)" ]
then verdict 0 "the two arms are genuinely different binaries"
else verdict 1 "the two arms are the SAME binary -- every comparison below is vacuous"; fi
# and the module each arm actually imports is the one we just hashed
for arm in fixed unfixed; do
  case $arm in fixed) P=$FIXED; W=$SO_F ;; unfixed) P=$UNFIXED; W=$SO_U ;; esac
  M=$(cd /tmp && PYTHONPATH="$P" python -c 'import rmgpy.solver.polymer as p; print(p.__file__)')
  printf '  %-8s imports %s\n' "$arm" "$M"
  [ "$M" = "$W" ] || { printf 'FAIL: %s arm imported the wrong build\n' "$arm"; FAIL=1; }
done

################################################################################
step "2/3/4. the probe, on both arms, state-matched"
################################################################################
for arm in unfixed fixed; do
  case $arm in fixed) P=$FIXED ;; unfixed) P=$UNFIXED ;; esac
  (cd /tmp && PYTHONPATH="$P" python \
     "$FIXED/test/rmgpy/solver/i060_tail_probe.py" "$T/probe.$arm.json") \
     >"$T/probe.$arm.out" 2>"$T/probe.$arm.err"
  [ -s "$T/probe.$arm.json" ] || { printf 'FAIL: %s probe produced nothing\n' "$arm"; FAIL=1; }
done
# ...and the UNFIXED build replayed on the FIXED build's own accepted states,
# so the per-step comparison below is one build against another AT THE SAME
# STATE, not two different trajectories talked about as if they were one.
python -c "
import json,sys
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
hx = float.fromhex
ok = True
print(f"  unfixed build: {u['build']}")
print(f"  fixed   build: {f['build']}")
if u['build'] == f['build']:
    print("  FAIL: both arms loaded the SAME extension"); sys.exit(1)

print("\n  [2] THE DEFECT -- the outlet switches off at the cutoff, and the")
print("      residual JUMPS doing it. F is the handshake flux into the")
print("      explicit DP=xs species, which nothing else in the fixture feeds.")
print(f"  {'tail mean DP':>14}{'F unfixed':>15}{'F fixed':>15}"
      f"{'d(mu1-(xs+1)mu0)/dt unf':>26}{'  fixed':>14}")
for a, b in zip(u['field'], f['field']):
    print(f"  {a['mean']:>14.6f}{hx(a['handshake']):>15.6e}"
          f"{hx(b['handshake']):>15.6e}{hx(a['de_tail']):>26.6e}"
          f"{hx(b['de_tail']):>14.6e}")
# the jump, unfixed, across mean = xs
hi = [a for a in u['field'] if abs(a['mean'] - (XS + 1e-6)) < 1e-12][0]
at = [a for a in u['field'] if a['mean'] == float(XS)][0]
jhi, jat = hx(hi['handshake']), hx(at['handshake'])
print(f"\n  unfixed: F(mean = xs + 1e-6) = {jhi:.6e}, F(mean = xs) = {jat:.6e}"
      f"  -> relative step {abs(jhi - jat) / jhi:.6f}")
if not (jhi > 0.0 and jat == 0.0):
    ok = False; print("  FAIL: the unfixed arm does not show the switch-off")
if abs(jhi - jat) / jhi < 0.5:
    ok = False; print("  FAIL: the unfixed step is not a jump")
# ...and the sign flip of the chain-count derivative across the same surface
print(f"  unfixed: dmu0(mean = xs + 1e-6) = {hx(hi['dmu0']):.6e}, "
      f"dmu0(mean = xs) = {hx(at['dmu0']):.6e}")
if not (hx(hi['dmu0']) < 0.0 < hx(at['dmu0'])):
    ok = False; print("  FAIL: expected the chain-count derivative to flip sign")

print("\n  [3] FIXED, same states: continuous, and the flux survives")
bhi = [a for a in f['field'] if abs(a['mean'] - (XS + 1e-6)) < 1e-12][0]
bat = [a for a in f['field'] if a['mean'] == float(XS)][0]
khi, kat = hx(bhi['handshake']), hx(bat['handshake'])
print(f"  fixed  : F(mean = xs + 1e-6) = {khi:.6e}, F(mean = xs) = {kat:.6e}"
      f"  -> relative step {abs(khi - kat) / khi:.3e}")
if abs(khi - kat) / khi > 1e-4:
    ok = False; print("  FAIL: the fixed arm still jumps across the cutoff")
below = [b for b in f['field'] if 1.0 < b['mean'] < XS]
if not all(hx(b['handshake']) > 0.0 for b in below):
    ok = False; print("  FAIL: fixed arm still has zero flux below the cutoff")
else:
    print(f"  fixed  : flux is strictly positive at all {len(below)} swept "
          f"states with 1 < mean < xs (unfixed: exactly 0.0 at every one)")
edge = [b for b in f['field'] if b['mean'] == 1.0]
if not all(hx(b['handshake']) == 0.0 for b in edge):
    ok = False; print("  FAIL: flux is nonzero at the cone edge mu1 == mu0")
else:
    print("  fixed  : EXACTLY 0.0 at mu1 == mu0, so the cone edge stays an "
          "invariant set")

print("\n  [4] two-arm, per-step, ON THE SAME STATES (unfixed replayed on the")
print("      fixed trajectory's own accepted states).")
tf = f['traj']['rows'][::10]
assert len(tf) == len(r['replay']), (len(tf), len(r['replay']))
n_worse = n_better = n_shut = 0
worst = 0.0
for a, b in zip(tf, r['replay']):
    df = a['dmu1'] - (XS + 1) * a['dmu0']          # fixed, at this state
    du = hx(b['dmu1']) - (XS + 1) * hx(b['dmu0'])  # unfixed, SAME state
    Ffix, Funf = a['handshake'], hx(b['handshake'])
    # The identity the fix rests on: the handshake contributes EXACTLY +F to
    # the drift of the tail's support invariant (dmu0 -= F, dmu1 -= xs*F), and
    # the two builds differ in nothing else. So the whole per-state difference
    # in that drift must be the difference of the two fluxes -- no residue.
    worst = max(worst, abs((df - du) - (Ffix - Funf)))
    if Funf == 0.0 and Ffix > 0.0: n_shut += 1
    if df > du: n_better += 1
    elif df < du: n_worse += 1
print(f"  {len(tf)} accepted states replayed. On every one,")
print(f"    d(mu1-(xs+1)mu0)/dt_fixed - _unfixed == F_fixed - F_unfixed")
print(f"    exactly; largest residue over all states {worst:.3e}")
if worst > 1e-12:
    ok = False; print("  FAIL: the drift difference is not accounted for by the flux alone")
print(f"  the unfixed arm's outlet is shut (F == 0.0 exactly) on {n_shut} of "
      f"those states, where the fixed arm's is open")
print(f"  drift strictly better on {n_better} states, worse on {n_worse}")
if n_worse:
    ok = False; print("  FAIL: the fix made the tail-support drift worse somewhere")

print("\n  [5] NO KINETICS CHANGE where the removed boolean was already true.")
print("      Hex-float EQUALITY on all five residual components; approximate")
print("      agreement would hide exactly the regression at issue.")
n = same = 0
for key in ("field", "field_mono", "random"):
    for a, b in zip(u[key], f[key]):
        if a['mean'] <= XS + 1e-9:      # not the no-op region
            continue
        for fld in ("handshake", "dmu0", "dmu1", "dmu2", "gas"):
            n += 1
            if a[fld] == b[fld]:
                same += 1
            else:
                ok = False
                print(f"  DIFFERS {key} mean={a['mean']!r} {fld}: "
                      f"unfixed={hx(a[fld])!r} fixed={hx(b[fld])!r}")
cone = [a for a in u['random'] if a['tag'] == 'cone']
print(f"  compared {n} residual components bit for bit across "
      f"{n // 5} states above the cutoff")
print(f"  ({len(cone)} of them randomized, seed 20260904: mu0 over 8 decades, "
      f"mean DP from xs+1e-9 to xs+1e3, PDI in (1, 6])")
print("  every one bit-for-bit identical" if same == n else "  NOT identical")

print("\n  [6] THE INVARIANTS, asserted at EVERY accepted step of the fixed")
print("      trajectory -- not at the endpoints.")
rows = f['traj']['rows']
u_rows = u['traj']['rows']
held0 = rows[0]['mu1']
bad_cone = [x for x in rows if x['mu1'] - x['mu0'] < 0.0]
bad_var = [x for x in rows if x['mu0'] * x['mu2'] - x['mu1'] ** 2 < 0.0]
bad_mass = [x for x in rows
            if abs(XS * x['pxs'] + x['gas'] + x['mu1'] - held0)
            > 1e-9 * held0]
print(f"  steps checked: {len(rows)}")
print(f"  mu1 >= mu0        : violated on {len(bad_cone)} step(s); "
      f"min(mu1 - mu0) = {min(x['mu1'] - x['mu0'] for x in rows):.6e} "
      f"(unfixed {min(x['mu1'] - x['mu0'] for x in u_rows):.6e})")
print(f"  mu0*mu2 >= mu1^2  : violated on {len(bad_var)} step(s); "
      f"min = {min(x['mu0'] * x['mu2'] - x['mu1'] ** 2 for x in rows):.6e}")
print(f"  repeat-unit ledger: xs*explicit + gas + mu1 == mu1(0) = {held0:g} "
      f"to 1e-9 rel on {len(rows) - len(bad_mass)}/{len(rows)} steps")
if bad_cone or bad_var or bad_mass:
    ok = False
    print("  FAIL: an invariant asserted above is violated on the fixed arm")
# The tail's REPRESENTATION invariant is reported, not asserted, and both arms
# are shown -- it is violated on BOTH, by the scission kernel, which
# manufactures sub-cutoff chains inside the tail and is not this ticket.
mf = min(x['e_tail'] for x in rows)
mu_ = min(x['e_tail'] for x in u_rows)
print(f"  FINDING, not a pass/fail: the tail's own support invariant "
      f"mu1 >= (xs+1)*mu0 is violated on BOTH arms "
      f"(min fixed {mf:.6e}, min unfixed {mu_:.6e}). The scission kernel "
      f"manufactures sub-cutoff chains inside the tail and nothing routes "
      f"them out; the handshake can only slow it. Pre-existing, and better "
      f"fixed than unfixed by {mf - mu_:.6e}.")

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
verdict "$?" "defect reproduces unfixed, is fixed, exact no-op above the cutoff, invariants hold every step"

################################################################################
step "5. the added tests: RED on the unfixed build, GREEN on the fixed one"
################################################################################
CLS=test/rmgpy/solver/solverPolymerTest.py::TestTailExhaustionHandshake
for arm in unfixed fixed; do
  case $arm in fixed) P=$FIXED ;; unfixed) P=$UNFIXED ;; esac
  (cd /tmp && PYTHONPATH="$P" python -m pytest "$FIXED/$CLS" \
     --no-cov -q --tb=no -p no:cacheprovider \
     --junitxml="$T/new.$arm.xml" >"$T/new.$arm.log" 2>&1)
  printf '  %-8s %s\n' "$arm" "$(grep -E '^[0-9]+ (passed|failed)|passed|failed' "$T/new.$arm.log" | tail -1)"
done
python - "$T/new.unfixed.xml" "$T/new.fixed.xml" <<'PY'
import sys, xml.etree.ElementTree as ET
def load(p):
    r = {}
    for tc in ET.parse(p).getroot().iter("testcase"):
        k = "passed"
        for ch in tc:
            if ch.tag in ("failure", "error"): k = ch.tag
            elif ch.tag == "skipped": k = "skipped"
        r[tc.get("name")] = k
    return r
pre, post = load(sys.argv[1]), load(sys.argv[2])
red = sorted(k for k, v in pre.items() if v != "passed")
grn = sorted(k for k, v in post.items() if v != "passed")
print(f"  RED before ({len(red)}):")
for k in red: print(f"    {k}")
print(f"  still failing after ({len(grn)}): {grn or 'none'}")
both = sorted(k for k in pre if pre[k] == "passed" and post[k] == "passed")
print(f"  green on BOTH arms ({len(both)}) -- these pin properties rather than")
print( "  discriminate, and are reported as such, not counted as evidence:")
for k in both: print(f"    {k}")
ok = len(red) >= 1 and not grn
sys.exit(0 if ok else 1)
PY
verdict "$?" "the added tests fail before and pass after"

################################################################################
step "6. the full existing suite, both arms, pass/fail sets named"
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
print(f"  failing on the UNFIXED arm ({len(prf)}) -- all pre-existing:")
for k in sorted(prf): print(f"    {k}")
print(f"  failing on the FIXED arm ({len(fin)}):")
for k in sorted(fin): print(f"    {k}")
changed = sorted(k for k in set(pre) | set(post) if pre.get(k) != post.get(k))
print(f"  tests whose state CHANGED between the arms: {len(changed)}")
for k in changed: print(f"    {k}: {pre.get(k)} -> {post.get(k)}")
# The added class is collected on both arms (it lives in the FIXED tree, which
# supplies the test file to both), so the state changes it produces are the
# ones named in step 5 and nowhere else. Anything else changing state is a
# regression.
expected = {k for k in changed if "TestTailExhaustionHandshake" in k}
if len(changed) != len(expected):
    ok = False
    print(f"  FAIL: {len(changed) - len(expected)} state change(s) outside the "
          f"added class")
else:
    print(f"  every one of the {len(changed)} change(s) is in "
          f"TestTailExhaustionHandshake, and the count matches step 5")
if c(post, "passed") < c(pre, "passed"):
    ok = False
    print(f"  FAIL: passed count fell from {c(pre,'passed')} to {c(post,'passed')}")
sys.exit(0 if ok else 1)
PY
verdict "$?" "no regression; every state change is the added class and the count matches"

################################################################################
step "7. nothing out of scope changed"
################################################################################
printf '  files changed on this branch vs its base:\n'
git -C "$FIXED" diff --stat "$UNFIXED_SHA" -- | sed 's/^/    /'
printf '  worktree status:\n'
git -C "$FIXED" status --porcelain | sed 's/^/    /'
printf '    (docs/contracts and build artifacts are git-ignored)\n'
printf '  baseline tree left untouched: %s\n' \
  "$(git -C "$UNFIXED" status --porcelain | wc -l) modified/untracked entries"

################################################################################
printf '\n================================================================\n'
if [ "$FAIL" -eq 0 ]; then printf 'I-060 VERIFIER: ALL CHECKS PASS\n'
else printf 'I-060 VERIFIER: FAILURES ABOVE\n'; fi
printf '================================================================\n'
rm -rf "$T"
exit "$FAIL"
