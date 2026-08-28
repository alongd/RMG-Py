#!/usr/bin/env bash
# Serial per-branch green/red sweep for the I-123 fourth pass.
#
# One stepping worktree, checked out detached at the authoritative local shared branch and
# then at each contributing branch tip in turn, check-pydas'd, built in place, and run with
# the same unit-suite command as the union. The shared branch runs FIRST and through the very
# same worktree and command, so its count is a baseline for the union's rather than a number
# from somewhere else.
#
# Serial on purpose: concurrent pytest runs across worktrees share a fixed scratch path
# outside every worktree and produce false failures.
#
# The build is incremental (`make build`, no `make clean`) exactly as the third pass's sweep
# was, so the counts are comparable with that pass's. git checkout refreshes the mtime of
# every file it changes, which is what setup.py's dependency check reads.

set -u
source ~/anaconda3/etc/profile.d/conda.sh
conda activate rmg_env

WT=/home/alon/Code/RMG-Py-i123d-step
L=/home/alon/Code/RMG-Py-i123d-audit/docs/i123d-audit/evidence/sweep
mkdir -p "$L"

BASELINE_DB=/home/alon/Code/RMG-database-i123d-baseline/input   # detached at RMG-database `plasma`
UNION_DB=/home/alon/Code/RMG-database-i123d-audit/input         # this pass's database union
FIRSTPASS_DB=/home/alon/Code/RMG-database-i123d-step1/input     # detached at i123-integration-db

# ref  database-it-is-paired-with
BRANCHES=(
  "plasma                       $BASELINE_DB"
  "i110-make-guard              $BASELINE_DB"
  "i119-rr-registry             $BASELINE_DB"
  "i115-preflight-deck          $BASELINE_DB"
  "i112-marcus-work-terms       $BASELINE_DB"
  "i102-quarantine              $BASELINE_DB"
  "i126-chemkin-electrons       $FIRSTPASS_DB"
  "i134-duplicate-drops-sink    $UNION_DB"
  "i135-tdep-roundtrip          $UNION_DB"
  "i145-units-laundering        $UNION_DB"
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
  cat "$out.paths.log" | sed 's/^/  /'

  python -m pytest -m "not functional and not database" > "$out.suite.log" 2>&1
  echo "  suite exit: $?"
  tail -1 "$out.suite.log" | sed 's/^/  /'
  grep -E '^FAILED|^ERROR' "$out.suite.log" | sed 's/^/  /'
  echo
done

echo "################ SWEEP DONE ################"
