#!/usr/bin/env bash
# Run one production deck against the I-065 build, both streams persisted.
#   run_deck.sh <deck-dir> <wall-seconds>
set -uo pipefail
DIR="$1"; WALL="$2"
# Derived from this script's own location, never hardcoded: a deck runner that
# names its birth worktree keeps running the branch build after the change is
# merged, and dies once that worktree is deleted.
WT=${I065_FIXED_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd -P)}

source /home/alon/anaconda3/etc/profile.d/conda.sh
conda activate rmg_env

# Fresh-build refusal, same rule as the production launcher: a stale .so makes
# a green run meaningless.
PYX="$WT/rmgpy/solver/polymer.pyx"
SO=$(ls -1 "$WT"/rmgpy/solver/polymer*.so | head -1)
[ "$SO" -nt "$PYX" ] || { echo "REFUSE: stale build" >&2; exit 1; }

cd "$DIR" || exit 1
{
  echo "deck=$DIR"
  echo "engine=$WT @ $(git -C "$WT" rev-parse HEAD)"
  echo "so=$SO ($(date -r "$SO" '+%F %T'))"
  echo "pyx=$(date -r "$PYX" '+%F %T')"
  echo "wall_cap_s=$WALL"
  echo "started_utc=$(date -u '+%FT%TZ')"
} > "$DIR/meta.txt"

timeout --signal=TERM "$WALL" python -c "
import sys
sys.path.insert(0, '$WT')
sys.argv = ['rmg.py', '$DIR/input.py']
from rmgpy import __main__
__main__.main()
" > >(tee -a "$DIR/stdout.log") 2> >(tee -a "$DIR/stderr.log" >&2)
rc=$?
{ echo "exit_rc=$rc"; echo "finished_utc=$(date -u '+%FT%TZ')"; } >> "$DIR/meta.txt"
exit "$rc"
