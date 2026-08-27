#!/usr/bin/env bash
# Serial per-branch green/red sweep for the I-123 third pass.
# One stepping worktree, checked out detached at each contributing branch tip in turn,
# check-pydas'd, built in place, and run with the same unit-suite command.
# Serial on purpose: concurrent pytest runs across worktrees share a fixed QM scratch path
# and produce false failures.

set -u
source ~/anaconda3/etc/profile.d/conda.sh
conda activate rmg_env

WT=/home/alon/Code/RMG-Py-i123c-step
L=/tmp/claude-1000/-home-alon-Code-RMG-Py-i123c-audit/098e215d-fcb0-4d98-86b2-f99bf77793a7/scratchpad/logs
mkdir -p "$L"

BASELINE_DB=/home/alon/Code/RMG-database-i123c-baseline/input
UNION_DB=/home/alon/Code/RMG-database-i123c-audit/input
FIRSTPASS_DB=/home/alon/Code/RMG-database-i123c-step1/input

# branch  db
BRANCHES=(
  "i110-make-guard              $BASELINE_DB"
  "i119-rr-registry             $BASELINE_DB"
  "i115-preflight-deck          $BASELINE_DB"
  "i112-marcus-work-terms       $BASELINE_DB"
  "i102-quarantine              $BASELINE_DB"
  "i126-chemkin-electrons       $FIRSTPASS_DB"
  "i134-duplicate-drops-sink    $UNION_DB"
  "i135-tdep-roundtrip          $UNION_DB"
)

cd "$WT" || exit 1
export PYTHONPATH="$WT:${PYTHONPATH:-}"

for entry in "${BRANCHES[@]}"; do
  set -- $entry
  branch=$1
  db=$2
  out="$L/branch-$branch"
  echo "################ $branch  (db=$db) ################"

  git checkout -- . 2>/dev/null
  git clean -fdq -e '*.so' -e '*.c' -e 'build' 2>/dev/null
  git checkout --detach "$branch" > "$out.checkout.log" 2>&1
  if [ $? -ne 0 ]; then echo "CHECKOUT FAILED for $branch"; continue; fi
  echo "  HEAD: $(git log --oneline -1)"

  printf 'database.directory = %s\n' "$db" > rmgrc

  python utilities.py check-pydas > "$out.checkpydas.log" 2>&1
  echo "  check-pydas exit: $?"

  make build > "$out.build.log" 2>&1
  bexit=$?
  echo "  make build exit: $bexit"
  if [ $bexit -ne 0 ]; then
    echo "  BUILD FAILED -- skipping suite for $branch"
    continue
  fi

  python -c "import rmgpy, rmgpy.molecule.molecule as m; from rmgpy import settings; \
print('rmgpy      :', rmgpy.__file__); print('compiled   :', m.__file__); \
print('database   :', settings['database.directory'])" > "$out.paths.log" 2>&1
  cat "$out.paths.log"

  python -m pytest -m "not functional and not database" > "$out.suite.log" 2>&1
  echo "  suite exit: $?"
  tail -1 "$out.suite.log" | sed 's/^/  /'
  grep -E '^FAILED|^ERROR' "$out.suite.log" | sed 's/^/  /'
  echo
done

echo "################ SWEEP DONE ################"
