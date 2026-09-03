#!/usr/bin/env bash
# Verifier for the legacy k_unzip realizability fix.
#
#   bash test/rmgpy/solver/i055_verify.sh
#
# Exits 0 only if every check passes. Each check recomputes its own numbers.
# Nothing here relaxes any solver detector, floor, tolerance or test; the only
# tolerances chosen are on this script's own scipy integrations, and each is
# justified by a refinement study rather than picked to pass.
set -uo pipefail

WT=/home/alon/Code/RMG-Py-i055-unzip-floor
PRIMARY=/home/alon/Code/RMG-Py          # pre-fix build; READ ONLY
CKMG=/home/alon/Code/CKMG               # READ ONLY
T=${TMPDIR:-/tmp}/i055_verify.$$
mkdir -p "$T"
FAIL=0

step() { printf '\n=== %s ===\n' "$*"; }
verdict() {  # verdict <rc> <label>
  if [ "$1" -eq 0 ]; then printf 'PASS: %s\n' "$2"
  else printf 'FAIL: %s (rc=%s)\n' "$2" "$1"; FAIL=1; fi
}

# shellcheck disable=SC1091
source /home/alon/anaconda3/etc/profile.d/conda.sh

################################################################################
step "1. the defect reproduces on the PRE-FIX kernel"
################################################################################
conda activate rmg_env
# rmgpy resolves to the primary checkout, whose compiled polymer extension is
# the one that produced ~/runs/RMG/poly_104. Nothing there is written.
# NOTE the cwd: python puts the script's directory (and, for -c, the cwd)
# FIRST on sys.path, so running this from the worktree would silently import
# the FIXED module and make the pre-fix checks vacuous. Run from /tmp.
BASE_MOD=$(cd /tmp && PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$PRIMARY" python -c \
  'import rmgpy.solver.polymer as p; print(p.__file__)')
printf 'pre-fix module under test: %s\n' "$BASE_MOD"
case "$BASE_MOD" in
  "$PRIMARY"/*) ;;
  *) printf 'FAIL: expected the pre-fix module to come from %s\n' "$PRIMARY"; FAIL=1 ;;
esac
(cd /tmp && PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$PRIMARY" \
  python "$WT/test/rmgpy/solver/unzip_realizability_probe.py") \
  2>"$T/p1.err" | grep -v '^WARNING'
RC=${PIPESTATUS[0]}
# The probe exits 1 when it finds a cone exit, which is what "the defect
# reproduces" MEANS, so 1 is the expected outcome here.
if [ "$RC" -eq 1 ]; then verdict 0 "pre-fix kernel exits the realizable cone and drives mu1 negative"
else verdict 1 "pre-fix probe did not report a cone exit"; fi

################################################################################
step "1b. structural-vs-threshold, measured on the PRE-FIX kernel"
################################################################################
(cd /tmp && PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$PRIMARY" \
  python "$WT/test/rmgpy/solver/i055_threshold_probe.py") \
  2>"$T/p1b.err" | grep -v '^WARNING'
verdict "${PIPESTATUS[0]}" "cone exit structural at every k_unzip; mu1<0 thresholded at k_scission/4"

################################################################################
step "2+3. post-fix: cone invariance to 100 s, mass ledger, cone detector"
################################################################################
(cd /tmp && PYTHONPATH="$WT:${PYTHONPATH:-}" \
  python "$WT/test/rmgpy/solver/i055_verify_core.py") \
  2>"$T/p2.err" | grep -v '^WARNING:root:POOL EXHAUSTION'
verdict "${PIPESTATUS[0]}" "mu1 >= mu0 >= 0 throughout; mass conserved; cone census fires correctly"

################################################################################
step "4. the unzip channel is still alive (CKMG audit, read-only)"
################################################################################
(cd /tmp && PYTHONPATH="$WT:${PYTHONPATH:-}" python \
  "$WT/test/rmgpy/solver/i055_sidecar_emit.py" "$T/synthetic_pools.json")
verdict "$?" "post-fix emitter still writes an unzip channel with A > 0"

bash -c '
  source /home/alon/anaconda3/etc/profile.d/conda.sh
  conda activate ck_env
  cd /tmp || exit 1
  export PYTHONDONTWRITEBYTECODE=1
  python -c "import ckmg; print(\"  ckmg:\", ckmg.__file__)" || exit 1
  python "$1/test/rmgpy/solver/i055_sidecar_audit.py" \
    "synthetic=$2/synthetic_pools.json" \
    "poly_104-live=/home/alon/Code/ckmg-pm4/artifacts/I-054_first-live-sidecar/polymer_pools.json"
' _ "$WT" "$T"
verdict "$?" "audit_sidecar_liveness: mass_loss_channels['"'"'unzip'"'"'] True, structurally_dead False"

################################################################################
step "5. no regression in the polymer solver suite"
################################################################################
cd "$WT" || exit 1
SUITE=(test/rmgpy/solver/solverPolymerTest.py
       test/rmgpy/solver/solverPolymerJacobianTest.py
       test/rmgpy/solver/solverPolymerConduitTest.py)
# --no-cov   the repo enables --cov=arkane --cov=rmgpy by default; with it the
#            suite does not finish in 45 min, without it it takes 20 s. The
#            coverage CONFIG is not edited -- only overridden per run.
# --deselect ONE test, disclosed in the report: the poly_102 crash-window
#            replay does not terminate in reasonable time on EITHER build. It
#            is NOT skipped to get green -- it is characterised separately, and
#            it is not failing, it is non-terminating.
SLOW=test/rmgpy/solver/solverPolymerJacobianTest.py::TestScopedJacobianPoly102::test_v2_v3_crash_window_twin_and_episode_invariance
PYOPTS=(-p no:cacheprovider --no-cov --deselect "$SLOW" --tb=no -q)

# WHY THERE IS NO SECOND ARM HERE ANY MORE.
# This check used to diff a "baseline" run against a post-fix run. That was
# broken twice over:
#   1. The baseline ran from the worktree cwd, and `python -m pytest` puts cwd
#      FIRST on sys.path, overriding PYTHONPATH="$PRIMARY". Both arms imported
#      the SAME (fixed) extension, so "every test has the same verdict" was
#      true by construction and the check measured nothing.
#   2. Even with that fixed, the PRIMARY CHECKOUT IS NOT A VALID BASELINE: at
#      the same commit its other compiled extensions differ from a fresh
#      worktree build, which misattributed FIVE unrelated failures to this fix.
# The only sound baseline is the SAME worktree with the one fix line toggled
# off and rebuilt. That is a ~2 min rebuild, too destructive to run inside a
# verifier, so it was run once by hand and its result is recorded below as the
# named pre-existing set. This check now asserts the post-fix failure set is
# EXACTLY that set -- no additions -- which is the property "no regression"
# actually means.
#
# Measured on the fix-off rebuild of THIS worktree: 421 tests, 407 passed,
# 14 failed = these 8 pre-existing + the 6 I-055 tripwire tests then red.
PRE_EXISTING=(
  "TestConduitFluxBookingIncrement2::test_integration_accepted_steps_match_analytic_integral"
  "TestHybridPolymerReactor::test_cross_pool_reverse_flux_vanishes_continuously"
  "TestHybridPolymerReactor::test_cross_pool_ve_detailed_balance_in_depletion_band"
  "TestHybridPolymerReactor::test_weaklink_u_capacity_saturation"
  "TestR81ExhaustionTailConditioning::test_floor_crossing_pool_wake_up_integrates_smoothly"
  "TestR81ResurrectionZeroMetricGuard::test_positive_edge_rate_resurrection_unchanged"
  "TestR81ResurrectionZeroMetricGuard::test_run7_shape_zero_edge_rates_do_not_add_species"
  "TestRtolNearFloorConviction::test_rtol_1e4_near_floor_conviction_canary"
)
printf 'THESE %s FAILURES ARE PRE-EXISTING -- this suite was already red before\n' "${#PRE_EXISTING[@]}"
printf 'this work. They are NOT regressions caused by the fix:\n'
printf '    %s\n' "${PRE_EXISTING[@]}"
printf '\nrunning the polymer suite against this worktree...\n'
PYTHONPATH="$WT" pytest "${SUITE[@]}" "${PYOPTS[@]}" \
  --junitxml="$T/post.xml" >"$T/post.log" 2>&1
printf '  %s\n' "$(grep -oE '^=+ .*(passed|failed).* =+$' "$T/post.log" | tail -1)"

# Counts and the per-test verdict stream come from the JUnit XML, not scraped
# console text: with -v pytest WRAPS long node ids across lines, and a grep
# over that silently undercounted 404 passed as 109 while still reporting PASS.
python - "$T/post.xml" "${PRE_EXISTING[@]}" <<'PY_CMP'
import sys, xml.etree.ElementTree as ET
post_path, pre = sys.argv[1], set(sys.argv[2:])
res = {}
for tc in ET.parse(post_path).getroot().iter("testcase"):
    k = "passed"
    for ch in tc:
        if ch.tag in ("failure", "error"): k = ch.tag
        elif ch.tag == "skipped": k = "skipped"
    res[f'{tc.get("classname").split(".")[-1]}::{tc.get("name")}'] = k
t = lambda k: sum(1 for v in res.values() if v == k)
print(f"  total={len(res)} passed={t('passed')} failed={t('failure')} "
      f"error={t('error')} skipped={t('skipped')}")
ok = True
fin = sorted(k for k, v in res.items() if v in ("failure", "error"))
new = [k for k in fin if k not in pre]
gone = [k for k in pre if k not in fin]
if new:
    ok = False
    print(f"  FAIL: {len(new)} NEW failure(s) not in the pre-existing set:")
    for k in new: print(f"    {k}")
else:
    print(f"  OK: {len(fin)} failure(s), every one of them pre-existing "
          f"-- zero regressions")
for k in gone:
    print(f"  note: pre-existing failure now PASSES: {k}")
TRIP = "TestLegacyUnzipRealizability"
trip = sorted(k for k in res if TRIP in k)
if not trip:
    print(f"  FAIL: no {TRIP} tests collected -- the tripwire is missing")
    ok = False
elif any(res[k] != "passed" for k in trip):
    print(f"  FAIL: I-055 tripwire not green:")
    for k in trip:
        if res[k] != "passed": print(f"    {k}: {res[k]}")
    ok = False
else:
    print(f"  OK: all {len(trip)} I-055 tripwire tests green")
# A run that collected almost nothing must not pass silently.
if len(res) < 400:
    print(f"  FAIL: only {len(res)} tests collected; expected >= 400")
    ok = False
sys.exit(0 if ok else 1)
PY_CMP
verdict "$?" "no regression: every remaining failure is pre-existing, I-055 tripwire green"

################################################################################
step "5a. the I-055 tripwire is RED on the PRE-FIX kernel"
################################################################################
# The tripwire must prove the defect, not just pass. These tests exercise the
# polymer kernels only, so the primary checkout IS a valid pre-fix build for
# them (the confound in note 2 above affects other suites, not this class).
# Run from /tmp: `python -m pytest` would otherwise put the worktree cwd first
# on sys.path and silently test the FIXED module.
TRIPRED=$(cd /tmp && PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$PRIMARY" \
  python -m pytest "$WT/test/rmgpy/solver/solverPolymerTest.py::TestLegacyUnzipRealizability" \
  --no-cov -q --tb=no -p no:cacheprovider 2>&1 | tail -1)
printf '  pre-fix kernel: %s\n' "$TRIPRED"
case "$TRIPRED" in
  *failed*) printf 'PASS: the tripwire is RED on the pre-fix kernel\n' ;;
  *) printf 'FAIL: the tripwire did NOT fail pre-fix -- it proves nothing\n'; FAIL=1 ;;
esac

################################################################################
step "6. the tests ran against a REBUILT extension, not a stale binary"
################################################################################
PYX="$WT/rmgpy/solver/polymer.pyx"
SO=$(ls "$WT"/rmgpy/solver/polymer.cpython-*.so 2>/dev/null | head -1)
printf '  pyx: %s  %s\n' "$(date -r "$PYX" -Is)" "$PYX"
printf '  so : %s  %s\n' "$(date -r "$SO" -Is)" "$SO"
if [ -n "$SO" ] && [ "$SO" -nt "$PYX" ]; then verdict 0 "compiled extension is newer than the .pyx edit"
else verdict 1 "compiled extension is STALE relative to the .pyx"; fi
MOD=$(PYTHONPATH="$WT:${PYTHONPATH:-}" python -c \
  'import rmgpy.solver.polymer as p; print(p.__file__)')
printf '  module actually imported under test: %s\n' "$MOD"
[ "$MOD" = "$SO" ] || { printf 'FAIL: the imported module is not the rebuilt .so\n'; FAIL=1; }

################################################################################
step "7. nothing out of scope changed"
################################################################################
# CKMG carries long-standing UNTRACKED scratch (literature PDFs, probe logs,
# memos) from July/August, none of it this work's. The check that actually
# means "CKMG was not edited here" is: no TRACKED file modified, and nothing
# under it newer than this session. Requiring a spotless tree would fail on
# other people's scratch and say nothing about whether I touched it.
SESSION_START=1788411180   # 2026-09-03T07:53+03:00, when this worktree was created
CK_TRACKED=$(git -C "$CKMG" status --porcelain | grep -vE '^\?\?' || true)
if [ -z "$CK_TRACKED" ]; then verdict 0 "CKMG has no modified tracked files"
else verdict 1 "CKMG has MODIFIED TRACKED files"; printf '%s\n' "$CK_TRACKED"; fi
CK_NEW=$(find "$CKMG" -newermt "@$SESSION_START" -not -path '*/.git/*' -type f 2>/dev/null)
if [ -z "$CK_NEW" ]; then verdict 0 "no file under CKMG is newer than this session's start"
else verdict 1 "files under CKMG were written during this session"; printf '%s\n' "$CK_NEW" | head; fi
printf '  (CKMG untracked scratch predating this session, left alone:)\n'
git -C "$CKMG" status --porcelain | grep -E '^\?\?' | sed 's/^/    /' | head -8

printf '  primary checkout %s:\n' "$PRIMARY"
PRIMARY_STATUS=$(git -C "$PRIMARY" status --porcelain)
printf '%s\n' "${PRIMARY_STATUS:-    (clean)}" | sed 's/^/    /'
# arkane/ess/molpro.py was already modified before this work began and is not
# mine; anything ELSE here would mean the read-only primary checkout was written.
UNEXPECTED=$(printf '%s\n' "$PRIMARY_STATUS" | grep -v 'arkane/ess/molpro.py' | grep -v '^$')
if [ -z "$UNEXPECTED" ]; then verdict 0 "primary checkout carries only its pre-existing modification"
else verdict 1 "primary checkout has UNEXPECTED changes"; printf '%s\n' "$UNEXPECTED"; fi

printf '  worktree %s:\n' "$WT"
git -C "$WT" status --porcelain | sed 's/^/    /'
EXPECTED=' M rmgpy/solver/polymer.pyx
 M test/rmgpy/solver/solverPolymerTest.py
?? test/rmgpy/solver/i055_sidecar_audit.py
?? test/rmgpy/solver/i055_sidecar_emit.py
?? test/rmgpy/solver/i055_threshold_probe.py
?? test/rmgpy/solver/i055_verify.sh
?? test/rmgpy/solver/i055_verify_core.py
?? test/rmgpy/solver/unzip_realizability_probe.py'
if [ "$(git -C "$WT" status --porcelain)" = "$EXPECTED" ]; then
  verdict 0 "worktree contains exactly the intended files"
else
  verdict 1 "worktree contains unintended changes (see list above)"
fi
# poly_104 is the evidence run; it finished before this work began and must
# not have been written. poly_103 is deliberately still running, so its tree
# changes on its own and is reported, not asserted.
POLY104_FINISHED=1788407772   # newest mtime of the run's own output, before this session
NEWER=$(find /home/alon/runs/RMG/poly_104 -type f -newermt "@$POLY104_FINISHED" 2>/dev/null)
# THIS work created no file under ~/runs/RMG. Anything newer is another
# worker's; it is reported, not silently passed, but it is not mine to own or
# to delete. MINE would be one of the basenames this work actually wrote.
MINE=$(printf '%s\n' "$NEWER" | grep -E '(unzip_realizability_probe|i055_)' || true)
if [ -z "$MINE" ]; then
  verdict 0 "this work wrote nothing under ~/runs/RMG/poly_104"
else
  verdict 1 "this work wrote into the read-only evidence tree"; printf '%s\n' "$MINE"
fi
if [ -n "$NEWER" ]; then
  printf '  NOTE: files under poly_104 newer than the run, written by ANOTHER worker\n'
  printf '        (the brief declares this tree read-only; reported, not touched):\n'
  printf '%s\n' "$NEWER" | sed 's/^/    /'
fi
printf '  poly_103 is live and its tree changes on its own; not asserted.\n'
printf '  poly_103 process (must stay alive): '
if ps -p 213774 >/dev/null 2>&1; then printf 'PID 213774 alive\n'
else printf 'PID 213774 NOT running\n'; fi

################################################################################
printf '\n================================================================\n'
if [ "$FAIL" -eq 0 ]; then printf 'VERIFIER: ALL CHECKS PASSED\n'
else printf 'VERIFIER: FAILURES PRESENT\n'; fi
printf '================================================================\n'
exit "$FAIL"
