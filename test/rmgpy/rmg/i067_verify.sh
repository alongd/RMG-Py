#!/usr/bin/env bash
# Verifier for the polymer-tier generation bound.
#
#   bash test/rmgpy/rmg/i067_verify.sh [runs-root] [snapshot-dir] [pytest-after-report]
#
# Exits 0 only if every check passes. Every number it prints is recomputed here,
# from artifacts on disk, at the moment it is printed -- no claim made elsewhere
# is restated. Nothing here relaxes any tolerance, detector or test.
#
# The three arguments default to this session's scratchpad, which is where the
# comparison runs were made. Point them elsewhere to re-verify a different set.
set -uo pipefail

WT=/home/alon/Code/RMG-Py-i067-enlarge-bound
DEFAULT_SCRATCH=/tmp/claude-1000/-home-alon-Code-RMG-Py-i067-enlarge-bound/fe3e9b88-e327-4961-8ce0-1e2b0b8899cc/scratchpad

RUNS=${1:-$DEFAULT_SCRATCH/runs}
SNAP=${2:-$DEFAULT_SCRATCH}
PYTEST_AFTER=${3:-$DEFAULT_SCRATCH/pytest_after.txt}

# shellcheck disable=SC1091
source /home/alon/anaconda3/etc/profile.d/conda.sh
conda activate rmg_env

printf 'worktree:   %s\n' "$WT"
printf 'HEAD:       %s\n' "$(git -C "$WT" rev-parse HEAD)"
printf 'dirty:      %s\n' "$(git -C "$WT" status --porcelain | tr '\n' ';')"
printf 'python:     %s\n' "$(command -v python)"

# cwd matters: python puts the cwd first on sys.path for `-c`/`-m`, and
# PYTHONPATH here carries /home/alon/Code/RMG-Py (the shared checkout). Run from
# the worktree and prepend it explicitly, or the verifier measures the wrong engine.
cd "$WT" || exit 2
PYTHONPATH="$WT:${PYTHONPATH:-}" python "$WT/test/rmgpy/rmg/i067_verify.py" \
    "$RUNS" "$SNAP" "$PYTEST_AFTER"
rc=$?
printf '\ni067_verify.sh exiting %d\n' "$rc"
exit $rc
