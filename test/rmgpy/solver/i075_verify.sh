#!/usr/bin/env bash
# Verifier for I-075: landing the i061 DASPK moment-slot error-weight floor.
#
#   bash test/rmgpy/solver/i075_verify.sh
#
# Exits 0 only if every check passes. Every number is recomputed; nothing here
# restates a report. No test, tolerance, assertion or detector is relaxed.
#
# LOCATION-INDEPENDENT BY CONSTRUCTION. The tree under test is derived from
# THIS SCRIPT'S OWN PATH, and every commit it compares against is pinned by
# sha, never by a branch name. It therefore keeps working when run from the
# merged mainline with the feature branch deleted -- which is where it is
# meant to be run.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WT="$(cd "$HERE/../../.." && pwd)"

# --- pinned SHAs (never branch names) ---------------------------------------
# The tip of `polymer` immediately BEFORE this ticket's merge. The fix must be
# absent here and present at the branch's current tip.
POLYMER_PRE_MERGE=9e7e0c4d5e9572b2a67c1b1d8bb8880287a793e1
# The two commits this ticket lands.
DIAG_COMMIT=ee73b7890c79354c549c742bb0208cce4d797940
FIX_COMMIT=74df5e84df1782ebce5b80490aeeeab9715762c0
# The polymer-suite test the i061 fix turned from red to green.
FLOOR_TEST='test/rmgpy/solver/solverPolymerTest.py::TestR81ExhaustionTailConditioning::test_floor_crossing_pool_wake_up_integrates_smoothly'
# Non-terminating on every build measured (documented by the I-055 verifier and
# re-confirmed here): deselected from the suite arms and disclosed, never
# skipped to make a failure disappear -- it does not fail, it does not finish.
NONTERM='test/rmgpy/solver/solverPolymerJacobianTest.py::TestScopedJacobianPoly102::test_v2_v3_crash_window_twin_and_episode_invariance'

T="${TMPDIR:-/tmp}/i075_verify.$$"
mkdir -p "$T"
FAIL=0

step()    { printf '\n=== %s ===\n' "$*"; }
verdict() { if [ "$1" -eq 0 ]; then printf 'PASS: %s\n' "$2"
            else printf 'FAIL: %s (rc=%s)\n' "$2" "$1"; FAIL=1; fi; }

# shellcheck disable=SC1091
source /home/alon/anaconda3/etc/profile.d/conda.sh
conda activate rmg_env

printf 'tree under test : %s\n' "$WT"
printf 'scratch         : %s\n' "$T"

################################################################################
step "1. the extension under test is the one built in THIS tree"
################################################################################
# Every number below is void without this. A green test on a stale .so is the
# most convincing wrong answer available in this codebase.
for MOD in polymer base; do
  PYX="$WT/rmgpy/solver/$MOD.pyx"
  SO=$(ls "$WT"/rmgpy/solver/$MOD.cpython-*.so 2>/dev/null | head -1)
  if [ -z "$SO" ]; then
    printf 'FAIL: no compiled %s extension in %s/rmgpy/solver\n' "$MOD" "$WT"
    printf '      build it first:\n'
    printf '        source /home/alon/anaconda3/etc/profile.d/conda.sh && conda activate rmg_env\n'
    printf '        cd %s && python utilities.py check-pydas && python setup.py build_ext --inplace -j 8\n' "$WT"
    FAIL=1
    continue
  fi
  printf '  %-8s pyx %s  %s\n' "$MOD" "$(date -r "$PYX" -Is)" "$PYX"
  printf '  %-8s so  %s  %s\n' "$MOD" "$(date -r "$SO" -Is)" "$SO"
  printf '  %-8s sha256 %s\n' "$MOD" "$(sha256sum "$SO" | cut -d' ' -f1)"
  if [ "$SO" -nt "$PYX" ]; then
    verdict 0 "$MOD: compiled extension is newer than its .pyx"
  else
    verdict 1 "$MOD: compiled extension is STALE relative to its .pyx -- rebuild before believing any number"
  fi
  IMPORTED=$(cd /tmp && PYTHONPATH="$WT" python -c \
    "import rmgpy.solver.$MOD as m; print(m.__file__)")
  printf '  %-8s imported under test: %s\n' "$MOD" "$IMPORTED"
  if [ "$IMPORTED" = "$SO" ]; then
    verdict 0 "$MOD: the imported module IS the extension in this tree"
  else
    verdict 1 "$MOD: the imported module is NOT this tree's extension"
  fi
done

################################################################################
step "2. the mainline did not carry the fix before the merge, and does after"
################################################################################
# Asserted against the RECORDED pre-merge sha, not against wherever the branch
# points when this runs.
printf '  pre-merge polymer tip (pinned): %s\n' "$POLYMER_PRE_MERGE"
if git -C "$WT" cat-file -e "$POLYMER_PRE_MERGE^{commit}" 2>/dev/null; then
  BEFORE=$(git -C "$WT" grep -c _chem_atol_array "$POLYMER_PRE_MERGE" -- \
             rmgpy/solver/polymer.pyx 2>/dev/null | awk -F: '{s+=$NF} END{print s+0}')
  printf '  _chem_atol_array occurrences at %s: %s\n' \
         "${POLYMER_PRE_MERGE:0:9}" "$BEFORE"
  if [ "$BEFORE" -eq 0 ]; then
    verdict 0 "the pinned pre-merge mainline does NOT carry the fix"
  else
    verdict 1 "the pinned pre-merge mainline ALREADY carried the fix -- this ticket was void"
  fi
else
  verdict 1 "pinned pre-merge sha $POLYMER_PRE_MERGE is not in this repository"
fi

AFTER=$(git -C "$WT" grep -c _chem_atol_array HEAD -- rmgpy/solver/polymer.pyx \
          2>/dev/null | awk -F: '{s+=$NF} END{print s+0}')
printf '  _chem_atol_array occurrences at HEAD (%s): %s\n' \
       "$(git -C "$WT" rev-parse --short HEAD)" "$AFTER"
if [ "$AFTER" -ge 1 ]; then
  verdict 0 "the tree under test carries the fix"
else
  verdict 1 "the tree under test does NOT carry the fix"
fi

for C in "$DIAG_COMMIT" "$FIX_COMMIT"; do
  if git -C "$WT" merge-base --is-ancestor "$C" HEAD 2>/dev/null; then
    verdict 0 "${C:0:9} ($(git -C "$WT" log -1 --format=%s "$C")) is an ancestor of HEAD"
  else
    verdict 1 "${C:0:9} is NOT an ancestor of HEAD"
  fi
done
# This one is the LANDING check: it is red on the un-merged feature branch by
# construction (the branch was cut before the mainline moved on) and green only
# once the work has actually landed on top of the recorded mainline tip. Red
# here means "not landed yet", not "broken".
if git -C "$WT" merge-base --is-ancestor "$POLYMER_PRE_MERGE" HEAD 2>/dev/null; then
  verdict 0 "the pinned pre-merge mainline is an ancestor of HEAD -- the work is LANDED on top of it, nothing rewritten away"
else
  verdict 1 "the pinned pre-merge mainline is NOT an ancestor of HEAD -- run this from the MERGED tree (on the un-merged branch this is expected)"
fi

################################################################################
step "3+4+5+7. gate, descriptors, floor behaviour, scope guard (recomputed)"
################################################################################
(cd /tmp && PYTHONPATH="$WT" python "$HERE/i075_verify_core.py")
verdict "$?" "verifier core: env gate, bounded descriptors, floor unchanged, scope guard fires"

################################################################################
step "6. the previously-failing floor-crossing test passes"
################################################################################
printf '  %s\n' "$FLOOR_TEST"
(cd "$WT" && PYTHONPATH="$WT" python -m pytest "$FLOOR_TEST" \
   --no-cov -q --tb=short -p no:cacheprovider) >"$T/floor.log" 2>&1
RC=$?
tail -3 "$T/floor.log" | sed 's/^/    /'
verdict "$RC" "test_floor_crossing_pool_wake_up_integrates_smoothly passes (it was DASPK IDID=-6 before i061)"

################################################################################
step "8. the full unit suite, before vs after"
################################################################################
# The BEFORE arm is a RECORDED measurement: the same worktree, built from
# 74df5e84d (the i061 pair, before the i075 edits), full unit suite, same
# deselect. It cannot be recomputed here without destroying the build under
# test, so it is pinned below by name -- and the AFTER arm is run live and
# diffed against it. Any test that changes state must appear in EXPECTED_CHANGE.
BEFORE_TXT="$HERE/i075_suite_before.txt"
if [ ! -f "$BEFORE_TXT" ]; then
  verdict 1 "recorded BEFORE suite result $BEFORE_TXT is missing"
else
  printf '  recorded BEFORE: %s (see its header for the exact invocation)\n' "$BEFORE_TXT"
  printf '  running the AFTER arm with the SAME flags (~4 min)...\n'
  (cd "$WT" && PYTHONPATH="$WT" python -m pytest \
      -m "not functional and not database" \
      -p no:cacheprovider --no-cov -q --tb=no \
      --deselect "$NONTERM" \
      --junitxml="$T/suite_after.xml") >"$T/suite_after.log" 2>&1
  tail -2 "$T/suite_after.log" | sed 's/^/    /'
  python - "$BEFORE_TXT" "$T/suite_after.xml" <<'PY_CMP'
import sys, xml.etree.ElementTree as ET

# The ONLY state changes this work is allowed to cause, each with the reason it
# is expected. Anything else -- in either direction -- fails the check.
_WHY = ("re-pinned to the post-i061 slot partition: these two asserted "
        "np.all(atol_array == deck_atol), i.e. the pre-fix uniformity the "
        "floor exists to break. Now pinned per slot class AND on the "
        "chemistry copy -- strictly more than before (they are RED with the "
        "floor disabled, which the originals were not).")
EXPECTED_CHANGE = {
    "test.rmgpy.tools.polymerMomentsRunnerTest.TestAtolReplayParity::"
    "test_atol_reaches_solver_floors_and_numpy_consumer": _WHY,
    "test.rmgpy.tools.polymerMomentsRunnerTest.TestAtolReplayParity::"
    "test_default_tolerances_unchanged": _WHY,
}

def load_txt(path):
    res = {}
    for ln in open(path):
        if ln.startswith("#") or not ln.strip():
            continue
        outcome, nodeid = ln.rstrip("\n").split("\t", 1)
        res[nodeid] = outcome
    return res

def load_xml(path):
    res = {}
    for tc in ET.parse(path).getroot().iter("testcase"):
        k = "passed"
        for ch in tc:
            if ch.tag in ("failure", "error"):
                k = "failed"
            elif ch.tag == "skipped":
                k = "skipped"
        res["%s::%s" % (tc.get("classname"), tc.get("name"))] = k
    return res

before, after = load_txt(sys.argv[1]), load_xml(sys.argv[2])
def tally(d):
    return {k: sum(1 for v in d.values() if v == k)
            for k in ("passed", "failed", "skipped")}
print("  BEFORE total=%d %s" % (len(before), tally(before)))
print("  AFTER  total=%d %s" % (len(after), tally(after)))

ok = True
only_before = sorted(set(before) - set(after))
only_after = sorted(set(after) - set(before))
if only_before or only_after:
    ok = False
    print("  FAIL: the collected test set changed")
    for k in only_before[:20]:
        print("    disappeared: %s" % k)
    for k in only_after[:20]:
        print("    appeared:    %s" % k)

changed = {k: (before[k], after[k]) for k in set(before) & set(after)
           if before[k] != after[k]}
print("  state changes: %d" % len(changed))
for k, (b, a) in sorted(changed.items()):
    why = EXPECTED_CHANGE.get(k)
    print("    %s: %s -> %s%s" % (k, b, a, ("  [%s]" % why) if why else
                                  "  [UNEXPLAINED]"))
    if why is None:
        ok = False
if len(changed) != len(EXPECTED_CHANGE):
    ok = False
    print("  FAIL: %d state change(s) measured, %d explained"
          % (len(changed), len(EXPECTED_CHANGE)))
else:
    print("  OK: %d state change(s) measured, all %d explained"
          % (len(changed), len(EXPECTED_CHANGE)))

failed_after = sorted(k for k, v in after.items() if v == "failed")
print("  FAILING AFTER (%d) -- every one of them also failing BEFORE:"
      % len(failed_after))
for k in failed_after:
    print("    %s%s" % (k, "" if before.get(k) == "failed" else "   <-- NEW"))
if len(after) < 3000:
    ok = False
    print("  FAIL: only %d tests collected; a truncated run must not pass"
          % len(after))
sys.exit(0 if ok else 1)
PY_CMP
  verdict "$?" "full unit suite: collected set identical, every state change explained"
  printf '  NOTE: %s is DESELECTED in both arms.\n' "${NONTERM##*::}"
  printf '        It does not fail -- it does not terminate, on this build and on the\n'
  printf '        pre-i061 build alike (independently recorded by the I-055 verifier).\n'
fi

################################################################################
printf '\n================================================================\n'
if [ "$FAIL" -eq 0 ]; then printf 'I-075 VERIFIER: ALL CHECKS PASSED\n'
else printf 'I-075 VERIFIER: FAILURES PRESENT\n'; fi
printf '================================================================\n'
exit "$FAIL"
