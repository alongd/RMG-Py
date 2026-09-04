#!/usr/bin/env bash
# Verifier for the four I-065 merge blockers.
#
#   bash test/rmgpy/solver/i065_verify.sh
#
# Exits 0 only if every check passes. Every number is RECOMPUTED here; nothing
# is restated from a report. Nothing below relaxes any solver detector, floor,
# tolerance or test.
#
# TWO ARMS, and they are provably different builds:
#   FIXED   = this worktree (i065-merge-blockers), rebuilt in place
#   UNFIXED = /home/alon/Code/RMG-Py-i055-unzip-floor @ b34c787ec, which is
#             EXACTLY this branch's base with none of the four fixes. Its HEAD
#             and cleanliness are asserted, and it is READ ONLY throughout.
# For defects 3/4 the unfixed arm is /home/alon/Code/RMG-Py-i057-mass-reconcile
# @ 4a9abd3be, likewise read-only, since that is where those two live.
#
# Every python invocation runs from /tmp: `python -m pytest` and `python -c`
# put the cwd FIRST on sys.path, ahead of PYTHONPATH, so running from the
# worktree would silently import the FIXED extension into the unfixed arm and
# make every before/after comparison vacuous by construction (the trap that
# made I-055's first no-regression measurement meaningless).
set -uo pipefail

# FIXED is the tree UNDER TEST, and it is derived from this script's own
# location -- never hardcoded. A verifier that names its birth worktree keeps
# testing that worktree after the change is merged, so it reports the branch's
# result while claiming to report mainline's, and it dies outright once the
# branch worktree is deleted. Deriving it means this script always measures
# whatever tree it is shipped in, which is the only reading that stays true
# after this lands. Override with I065_FIXED_DIR only to test a foreign tree.
FIXED=${I065_FIXED_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd -P)}
# The two UNFIXED arms are the opposite case: they are pinned to specific
# commits ON PURPOSE, because they are the "before" builds the comparison is
# against, and they must NOT follow mainline. They are read-only siblings; if
# either worktree is ever deleted, step 0 fails loudly rather than silently
# comparing a build against itself.
UNFIXED=/home/alon/Code/RMG-Py-i055-unzip-floor          # READ ONLY
UNFIXED_SHA=b34c787ecc05ae16740a3bd9b8dce55a8be5687a
UNFIXED57=/home/alon/Code/RMG-Py-i057-mass-reconcile     # READ ONLY
UNFIXED57_SHA=4a9abd3be7ed9b3ca9b74beab8db8ecb8d23b0d9
DECKS=${I065_DECK_DIR:-}    # optional: dir holding the RMG production runs
# Session reference: the commit time of this branch's base, b34c787ec.
SESSION_START_TS=1788434439
T=${TMPDIR:-/tmp}/i065_verify.$$
mkdir -p "$T"
FAIL=0

step() { printf '\n=== %s ===\n' "$*"; }
verdict() {  # verdict <rc> <label>
  if [ "$1" -eq 0 ]; then printf 'PASS: %s\n' "$2"
  else printf 'FAIL: %s (rc=%s)\n' "$2" "$1"; FAIL=1; fi
}
export PYTHONDONTWRITEBYTECODE=1

# shellcheck disable=SC1091
source /home/alon/anaconda3/etc/profile.d/conda.sh
conda activate rmg_env

################################################################################
step "0. the two arms are different builds, and the unfixed one is untouched"
################################################################################
for pair in "$UNFIXED $UNFIXED_SHA" "$UNFIXED57 $UNFIXED57_SHA"; do
  set -- $pair
  H=$(git -C "$1" rev-parse HEAD)
  D=$(git -C "$1" status --porcelain)
  printf '  %s\n    HEAD=%s\n    dirty=%s\n' "$1" "$H" "${D:-<clean>}"
  [ "$H" = "$2" ] || { printf 'FAIL: %s is not at %s\n' "$1" "$2"; FAIL=1; }
  [ -z "$D" ] || { printf 'FAIL: %s is dirty -- not a usable baseline\n' "$1"; FAIL=1; }
done
SO_F=$(ls "$FIXED"/rmgpy/solver/polymer.cpython-*.so | head -1)
SO_U=$(ls "$UNFIXED"/rmgpy/solver/polymer.cpython-*.so | head -1)
H_F=$(sha256sum "$SO_F" | cut -c1-16)
H_U=$(sha256sum "$SO_U" | cut -c1-16)
printf '  fixed   .so sha256[:16] = %s\n  unfixed .so sha256[:16] = %s\n' "$H_F" "$H_U"
if [ "$H_F" != "$H_U" ]; then verdict 0 "the two arms are genuinely different binaries"
else verdict 1 "the two arms are the SAME binary -- every comparison below is vacuous"; fi

################################################################################
step "1. defect 1: the fabrication reproduces UNFIXED and is gone FIXED"
################################################################################
for arm in fixed unfixed; do
  case $arm in fixed) P=$FIXED ;; unfixed) P=$UNFIXED ;; esac
  (cd /tmp && PYTHONPATH="$P" python \
     "$FIXED/test/rmgpy/solver/i065_channel_probe.py") \
     >"$T/probe.$arm.json" 2>"$T/probe.$arm.err"
  [ -s "$T/probe.$arm.json" ] || { printf 'FAIL: %s probe produced nothing\n' "$arm"; FAIL=1; }
done
python - "$T/probe.unfixed.json" "$T/probe.fixed.json" <<'PY1'
import json, sys
u = json.load(open(sys.argv[1])); f = json.load(open(sys.argv[2]))
print(f"  unfixed build: {u['build']}")
print(f"  fixed   build: {f['build']}")
if u['build'] == f['build']:
    print("  FAIL: both arms loaded the SAME extension"); sys.exit(1)
key = lambda r: (r['channel'], tuple(r['mu']))
U = {key(r): r for r in u['rows']}
F = {key(r): r for r in f['rows']}
off = [tuple(s) for s in f['off_cone_states']]
ok, fabricated = True, 0.0
print("  OFF the realizable cone -- gas emitted with no repeat units behind it")
print(f"  {'channel':<20}{'mu':<28}{'gas UNFIXED':>16}{'gas FIXED':>18}")
for ch in ("legacy_k_unzip", "k_depropagation", "radical_qssa_unzip"):
    for mu in off:
        if mu[1] > 0.0 and mu[1] >= mu[0]:
            continue
        gu = float.fromhex(U[(ch, mu)]['gas'])
        gf = float.fromhex(F[(ch, mu)]['gas'])
        print(f"  {ch:<20}{str(mu):<28}{gu:>16.6e}{gf:>18.6e}")
        if mu[1] == 0.0:                     # zero units: MUST be exactly zero
            fabricated += gu
            if gf != 0.0:
                ok = False
                print(f"    FAIL: fixed build still emits {gf!r} at mu1 = 0")
        elif gf > gu:
            ok = False
            print("    FAIL: fixed build releases MORE than the unfixed one")
print(f"  total fabricated gas at mu1 = 0 (unfixed): {fabricated:.6e} mol/s")
if fabricated <= 0.0:
    ok = False
    print("  FAIL: the unfixed arm fabricated nothing -- the defect did not reproduce")
sys.exit(0 if ok else 1)
PY1
verdict "$?" "defect 1 reproduces unfixed, is exactly zero fixed"

################################################################################
step "2. defect 1 is a bit-for-bit no-op INSIDE the realizable cone"
################################################################################
python - "$T/probe.unfixed.json" "$T/probe.fixed.json" <<'PY2'
import json, sys
u = json.load(open(sys.argv[1])); f = json.load(open(sys.argv[2]))
key = lambda r: (r['channel'], tuple(r['mu']))
U = {key(r): r for r in u['rows']}
F = {key(r): r for r in f['rows']}
cone = [tuple(s) for s in f['cone_states']]
ok, n = True, 0
for ch in ("legacy_k_unzip", "k_depropagation", "radical_qssa_unzip"):
    for mu in cone:
        for fld in ("dmu0", "dmu1", "dmu2", "gas"):
            a, b = U[(ch, mu)][fld], F[(ch, mu)][fld]
            n += 1
            if a != b:
                ok = False
                da, db = float.fromhex(a), float.fromhex(b)
                print(f"  DIFFERS {ch} mu={mu} {fld}: unfixed={da!r} fixed={db!r} "
                      f"delta={db - da!r}")
print(f"  compared {n} derivatives across {len(cone)} fully-realizable states "
      f"x 3 channels")
print(f"  hand-picked: {[tuple(s) for s in cone[:5]]}")
print(f"  plus {len(cone) - 5} randomized states (seed 20260903), mu0 over 8 "
      f"decades, mean DP 1..1e4, PDI > 1, all satisfying mu1 >= mu0 > 0 AND "
      f"mu0*mu2 >= mu1^2")
print("  every one bit-for-bit identical" if ok else "  NOT identical (above)")
sys.exit(0 if ok else 1)
PY2
verdict "$?" "no kinetics change on any valid state (hex-float equality, not approx)"

################################################################################
step "3. the legacy/QSSA equivalence pin still passes, by name"
################################################################################
PIN=test/rmgpy/solver/solverPolymerTest.py::TestHybridPolymerReactor::test_qssa_handshake_equivalence_with_k_unzip
printf '  %s\n' "$PIN"
(cd /tmp && PYTHONPATH="$FIXED" python -m pytest "$FIXED/$PIN" \
   --no-cov -q --tb=short -p no:cacheprovider >"$T/pin.log" 2>&1)
RC=$?
tail -2 "$T/pin.log" | sed 's/^/  /'
verdict "$RC" "test_qssa_handshake_equivalence_with_k_unzip"
# ...and the other half of the pin: the added tests that assert the release
# law is carried by BOTH channels and is a no-op on the QSSA one.
(cd /tmp && PYTHONPATH="$FIXED" python -m pytest \
   "$FIXED/test/rmgpy/solver/solverPolymerTest.py::TestReleaseAvailabilityGate" \
   --no-cov -q --tb=short -p no:cacheprovider >"$T/gate.log" 2>&1)
RC=$?
tail -2 "$T/gate.log" | sed 's/^/  /'
verdict "$RC" "TestReleaseAvailabilityGate (all three release channels)"

################################################################################
step "4a. defect 2: the variance census fires on a constructed invalid state"
################################################################################
(cd /tmp && PYTHONPATH="$FIXED" python -m pytest \
   "$FIXED/test/rmgpy/solver/solverPolymerTest.py::TestAcceptedStateVarianceCensus" \
   --no-cov -q --tb=short -p no:cacheprovider >"$T/var.log" 2>&1)
RC=$?
tail -2 "$T/var.log" | sed 's/^/  /'
verdict "$RC" "TestAcceptedStateVarianceCensus (fires on mu0*mu2 < mu1^2, silent at the floors, never raises)"

################################################################################
step "4b. defect 2 over real DASPK trajectories, both arms"
################################################################################
# The census hook run per ACCEPTED snapshot, exactly as base.pyx's simulate()
# calls it, on the production decks' own pool parameters and t=0 moments plus a
# sweep around them. Reported, never suppressed: a nonzero count here is a
# FINDING, and the arms are compared so that a violation cannot be blamed on
# this fix when the unfixed build reaches the same state.
for arm in fixed unfixed; do
  case $arm in fixed) P=$FIXED ;; unfixed) P=$UNFIXED ;; esac
  (cd /tmp && PYTHONPATH="$P" python \
     "$FIXED/test/rmgpy/solver/i065_census_replay.py" "$T/replay.$arm.json") \
     >/dev/null 2>"$T/replay.$arm.err"
done
python - "$T/replay.unfixed.json" "$T/replay.fixed.json" <<'PY4'
import json, sys
u = {r['arm']: r for r in json.load(open(sys.argv[1]))}
f = {r['arm']: r for r in json.load(open(sys.argv[2]))}
print(f"  {'arm':<38}{'acc':>5}{'cone':>6}{'var':>5}"
      f"{'min(mu0mu2-mu1^2)':>20}{'min(mu1-mu0)':>15}")
ok = True
for k, r in f.items():
    print(f"  {k:<38}{r['accepted_states']:>5}{r['cone_hits']:>6}"
          f"{r['variance_hits']:>5}{r['min_mu0mu2_minus_mu1sq']:>20.4e}"
          f"{r['min_mu1_minus_mu0']:>15.4e}")
    if r['error']:
        print(f"      integrator: {r['error']}")
    # A violation the UNFIXED build also reaches is pre-existing, not caused
    # here. Compare the geometry, which both builds compute, rather than the
    # hit count, which only the fixed build can produce.
    du = u[k]['min_mu0mu2_minus_mu1sq']
    if r['variance_hits'] and du >= 0.0:
        ok = False
        print(f"      FAIL: variance violated on the FIXED build only "
              f"(unfixed min = {du:.4e}) -- this fix would be the cause")
    if r['variance_hits'] and du < 0.0:
        print(f"      FINDING: violation present on the UNFIXED build too "
              f"(min = {du:.4e}) -- pre-existing, now merely VISIBLE")
sys.exit(0 if ok else 1)
PY4
verdict "$?" "no variance violation is introduced by this change (violations found are pre-existing)"

################################################################################
step "4c. defect 2 over the RMG production deck runs"
################################################################################
if [ -z "$DECKS" ] || [ ! -d "$DECKS" ]; then
  printf '  NOT RUN. An RMG generation is hours, not seconds, so this verifier\n'
  printf '  scans deck runs rather than performing them, and it found none.\n'
  printf '  WHAT THIS LEAVES UNVERIFIED: whether the variance census fires\n'
  printf '  during a full RMG generation on a production deck. Nothing else in\n'
  printf '  this script covers that; step 4b covers the moment subsystem on the\n'
  printf '  decks parameters, not the decks generated mechanism.\n'
  printf '  To close it, run each deck and re-run with the directory named:\n'
  printf '    bash %s/test/rmgpy/solver/i065_run_deck.sh <run-dir> <wall-seconds>\n' "$FIXED"
  printf '    I065_DECK_DIR=<parent-of-run-dirs> bash %s/test/rmgpy/solver/i065_verify.sh\n' "$FIXED"
  verdict 1 "production-deck scan NOT RUN -- see the gap stated above"
else
  RC=0
  for d in "$DECKS"/*/; do
    [ -f "$d/RMG.log" ] || continue
    HITS=$(cat "$d"/RMG.log "$d"/stdout.log "$d"/stderr.log 2>/dev/null \
             | grep -c 'VARIANCE CENSUS')
    CONE=$(cat "$d"/RMG.log "$d"/stdout.log "$d"/stderr.log 2>/dev/null \
             | grep -c 'CONE CENSUS')
    ITERS=$(grep -c 'Model Enlargement Summary\|Generating initial reactions' \
              "$d/RMG.log" 2>/dev/null)
    printf '  %-28s variance=%s cone=%s  (rc=%s, %s)\n' \
      "$(basename "$d")" "$HITS" "$CONE" \
      "$(sed -n 's/^exit_rc=//p' "$d/meta.txt" 2>/dev/null | tail -1)" \
      "$(sed -n 's/^finished_utc=//p' "$d/meta.txt" 2>/dev/null | tail -1)"
    if [ "$HITS" -gt 0 ]; then
      printf '    FINDING -- offending states:\n'
      grep -h 'VARIANCE CENSUS' "$d"/RMG.log "$d"/stdout.log "$d"/stderr.log \
        2>/dev/null | sed 's/^/      /' | head -5
      RC=1
    fi
  done
  if [ "$RC" -eq 0 ]; then verdict 0 "the variance census never fired on a production deck run"
  else printf 'REPORTED (not a failure of the fix): the census fired on a production deck; states above\n'; fi
fi

################################################################################
step "5. defect 3: the Mn == 0.0 path does what the comment says"
################################################################################
(cd /tmp && PYTHONPATH="$FIXED" python -m pytest \
   "$FIXED/test/rmgpy/polymerTest.py" -k "mn_zero or zero_declared_mass or reconciled_initial_mass" \
   --no-cov -q --tb=short -p no:cacheprovider >"$T/d3.log" 2>&1)
RC=$?
tail -2 "$T/d3.log" | sed 's/^/  /'
verdict "$RC" "Mn == 0.0 skips the reconciliation and warns; zero declared mass is a real disagreement"
# and the same tests RED on the unfixed i057 arm
(cd /tmp && PYTHONPATH="$UNFIXED57" python -m pytest \
   "$FIXED/test/rmgpy/polymerTest.py" -k "mn_zero or zero_declared_mass" \
   --no-cov -q --tb=no -p no:cacheprovider >"$T/d3red.log" 2>&1)
LAST=$(tail -1 "$T/d3red.log")
printf '  on %s: %s\n' "$(basename "$UNFIXED57")" "$LAST"
case "$LAST" in
  *failed*) printf 'PASS: the defect-3 tests are RED on the unfixed arm\n' ;;
  *) printf 'FAIL: the defect-3 tests are not red unfixed -- they prove nothing\n'; FAIL=1 ;;
esac

################################################################################
step "6. the ordering constraint holds: the disagreement warning is still alive"
################################################################################
(cd /tmp && PYTHONPATH="$FIXED" I065_FIXED="$FIXED" python - <<'PY6'
import logging, os, sys
# Derived from FIXED, not hardcoded: see the note at the top of this script.
sys.path.insert(0, os.path.join(os.environ["I065_FIXED"], "test", "rmgpy"))
from polymerTest import _build_compile_inputs
from rmgpy.rmg.polymer_input import compile_polymer_phase

class Grab(logging.Handler):
    def __init__(self):
        super().__init__(); self.msgs = []
    def emit(self, r): self.msgs.append(r.getMessage())

# A deck whose two declarations disagree: initial_mass=1 kg with Mn=5000 implies
# mu0 = 0.2 chains, initialMoles says 0.01. Reconciling BEFORE the detector
# would equalise them by construction and silence it forever.
bp, im, sd, poly = _build_compile_inputs(moles=0.01)
before = poly.initial_mass_g
h = Grab(); root = logging.getLogger(); root.addHandler(h); root.setLevel(logging.WARNING)
compile_polymer_phase(bp, im, sd)
root.removeHandler(h)
hits = [m for m in h.msgs if "initialMoles" in m]
print(f"  declared initial_mass_g before = {before!r}")
print(f"  reconciled initial_mass_g after = {poly.initial_mass_g!r} (= mu0*Mn = {0.01*5000.0!r})")
print(f"  disagreement warnings fired: {len(hits)}")
for m in hits: print(f"    {m[:150]}")
ok = len(hits) == 1 and poly.initial_mass_g == 0.01 * 5000.0
# and the detector must NOT fire on a consistent deck (no false positive)
bp, im, sd, poly2 = _build_compile_inputs(moles=0.2)
h2 = Grab(); root.addHandler(h2)
compile_polymer_phase(bp, im, sd)
root.removeHandler(h2)
quiet = [m for m in h2.msgs if "initialMoles" in m]
print(f"  consistent deck warnings: {len(quiet)} (must be 0)")
ok = ok and not quiet
sys.exit(0 if ok else 1)
PY6
)
verdict "$?" "the detector still fires on a disagreeing deck and stays quiet on a consistent one"

################################################################################
step "7. no regression: the suite, both arms, from provably different builds"
################################################################################
# Deselected by NAME, not by --deselect: with the suite given as absolute
# paths and the invocation dir /tmp, pytest's --deselect nodeid does not match
# (measured: 16 collected, 0 deselected) and the non-terminating test runs.
SLOW=test_v2_v3_crash_window_twin_and_episode_invariance
SUITE=(test/rmgpy/solver/solverPolymerTest.py
       test/rmgpy/solver/solverPolymerJacobianTest.py
       test/rmgpy/solver/solverPolymerConduitTest.py
       test/rmgpy/polymerTest.py
       test/rmgpy/rmg/inputTest.py)
# --no-cov: the repo config enables --cov by default and the suite does not
#           finish in 45 min with it. The coverage CONFIG is not edited.
# --deselect ONE test, disclosed: the poly_102 crash-window replay does not
#           terminate in reasonable time on EITHER build. Not skipped to get
#           green -- it is non-terminating, not failing, on both arms alike.
for arm in unfixed fixed; do
  case $arm in fixed) P=$FIXED ;; unfixed) P=$UNFIXED ;; esac
  PATHS=(); for s in "${SUITE[@]}"; do PATHS+=("$FIXED/$s"); done
  (cd /tmp && PYTHONPATH="$P" python -m pytest "${PATHS[@]}" \
     --no-cov -q --tb=no -p no:cacheprovider -k "not $SLOW" \
     --junitxml="$T/suite.$arm.xml" >"$T/suite.$arm.log" 2>&1)
  MOD=$(cd /tmp && PYTHONPATH="$P" python -c \
    'import rmgpy.solver.polymer as p; print(p.__file__)')
  printf '  %-8s extension under test: %s\n' "$arm" "$MOD"
  case "$MOD" in "$P"/*) ;; *) printf 'FAIL: %s arm imported the wrong build\n' "$arm"; FAIL=1 ;; esac
done
python - "$T/suite.unfixed.xml" "$T/suite.fixed.xml" <<'PY7'
import sys, xml.etree.ElementTree as ET
def load(p):
    res = {}
    for tc in ET.parse(p).getroot().iter("testcase"):
        k = "passed"
        for ch in tc:
            if ch.tag in ("failure", "error"): k = ch.tag
            elif ch.tag == "skipped": k = "skipped"
        res[f'{tc.get("classname")}::{tc.get("name")}'] = k
    return res
pre, post = load(sys.argv[1]), load(sys.argv[2])
c = lambda r, k: sum(1 for v in r.values() if v == k)
for n, r in (("unfixed", pre), ("fixed", post)):
    print(f"  {n:<8} total={len(r):<5} passed={c(r,'passed'):<5} "
          f"failed={c(r,'failure') + c(r,'error'):<4} skipped={c(r,'skipped')}")
ok = c(post, "passed") >= c(pre, "passed")
if not ok:
    print(f"  FAIL: passed count fell from {c(pre,'passed')} to {c(post,'passed')}")
fin = {k for k, v in post.items() if v in ("failure", "error")}
prf = {k for k, v in pre.items() if v in ("failure", "error")}
# The 10 + 4 tests added by this work do not exist on the unfixed arm as
# PASSING tests, so they are expected to be red there; a NEW red on the fixed
# arm is a regression.
new = sorted(fin - prf)
if new:
    ok = False
    print(f"  FAIL: {len(new)} NEW failure(s) on the fixed build:")
    for k in new: print(f"    {k}")
else:
    print(f"  OK: {len(fin)} failure(s) on the fixed build, every one of them "
          f"also failing unfixed -- zero regressions")
cured = sorted(prf - fin)
print(f"  {len(cured)} test(s) red unfixed and green fixed (the added coverage):")
for k in cured: print(f"    {k}")
if len(post) < 800:
    ok = False
    print(f"  FAIL: only {len(post)} tests collected; expected >= 800")
sys.exit(0 if ok else 1)
PY7
verdict "$?" "no regression; passed count did not fall; every remaining failure is pre-existing"

################################################################################
step "8. the build is fresh -- not a stale binary"
################################################################################
PYX="$FIXED/rmgpy/solver/polymer.pyx"
printf '  pyx: %s  %s\n' "$(date -r "$PYX" -Is)" "$PYX"
printf '  so : %s  %s\n' "$(date -r "$SO_F" -Is)" "$SO_F"
if [ "$SO_F" -nt "$PYX" ]; then verdict 0 "compiled extension is newer than its .pyx"
else verdict 1 "compiled extension is STALE relative to the .pyx"; fi
MOD=$(cd /tmp && PYTHONPATH="$FIXED" python -c \
  'import rmgpy.solver.polymer as p; print(p.__file__)')
printf '  module imported under test: %s\n' "$MOD"
[ "$MOD" = "$SO_F" ] || { printf 'FAIL: the imported module is not the rebuilt .so\n'; FAIL=1; }

################################################################################
step "9. nothing out of scope changed"
################################################################################
for R in /home/alon/Code/CKMG /home/alon/Code/TA; do
  # Both carry long-standing UNTRACKED scratch from other work, and TA also
  # carries a tracked .gitignore edit from 2026-07-22. The brief's bar is
  # "nothing ATTRIBUTABLE to you", so the check is: no tracked file modified
  # AND touched since this session's reference time. Anything older is listed
  # with its date rather than passed silently -- the check is not weakened, it
  # is made attributable, since a dirty-tree equality test would fail forever
  # on other people's edits and say nothing about whether this work touched
  # the repo.
  printf '  %s\n' "$R"
  N=0
  while read -r _ f; do
    [ -n "$f" ] || continue
    TS=$(stat -c %Y "$R/$f" 2>/dev/null || echo 0)
    if [ "$TS" -gt "$SESSION_START_TS" ]; then
      printf '    MODIFIED THIS SESSION: %s (%s)\n' "$f" "$(date -d @"$TS" -Is)"
      N=$((N + 1))
    else
      printf '    pre-existing, untouched here: %s (%s)\n' "$f" "$(date -d @"$TS" -Is)"
    fi
  done < <(git -C "$R" status --porcelain | grep -v '^??' || true)
  if [ "$N" -eq 0 ]; then printf '    no tracked file modified since %s\n' \
      "$(date -d @"$SESSION_START_TS" -Is)"
  else printf 'FAIL: %s has tracked files modified during this session\n' "$R"; FAIL=1; fi
done
# /home/alon/runs: this session never wrote there (its deck copies live in the
# session scratchpad). Filesystem mtimes cannot attribute a write to a process,
# so what is asserted is the stronger, checkable thing: the files this session
# READ are byte-identical to their pinned digests, and anything else newer is
# listed rather than silently passed.
while read -r want path; do
  [ -e "$path" ] || { printf 'FAIL: %s vanished\n' "$path"; FAIL=1; continue; }
  got=$(sha256sum "$path" | cut -d' ' -f1)
  if [ "$got" = "$want" ]; then printf '  unchanged: %s\n' "$path"
  else printf 'FAIL: %s CHANGED (%s != %s)\n' "$path" "$got" "$want"; FAIL=1; fi
done <<'RUNS'
0b9b55440ac3c869c2264badba6c72b9b778a0ea046e2b1c72b2eb5674cd842b /home/alon/runs/RMG/poly_105/input.py
RUNS
printf '  files under /home/alon/runs newer than this session start:\n'
find /home/alon/runs -newermt "@$SESSION_START_TS" -type f 2>/dev/null \
  | sed 's/^/    /' | head -20
printf '    (attribution by mtime is not possible; any entry above belongs to a\n'
printf '     concurrent session -- this one only READ /home/alon/runs)\n'
printf '  this worktree:\n'
git -C "$FIXED" status --porcelain | sed 's/^/    /'
printf '    (empty = everything committed; docs/contracts is git-ignored)\n'

################################################################################
printf '\n================================================================\n'
if [ "$FAIL" -eq 0 ]; then printf 'I-065 VERIFIER: ALL CHECKS PASS\n'
else printf 'I-065 VERIFIER: FAILURES ABOVE\n'; fi
printf '================================================================\n'
rm -rf "$T"
exit "$FAIL"
