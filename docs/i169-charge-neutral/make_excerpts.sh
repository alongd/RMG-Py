#!/bin/bash
# Regenerate the RMG.log.charge-excerpt files from the full RMG.log of each run.
#
# The full logs are NOT committed: they are 500+ lines of generated model-building
# output, of which about ten lines bear on this ticket. These excerpts carry those
# lines with their ORIGINAL line numbers, plus the two censuses the Verifier claims
# rest on (every Warning: line in the run, and every charge/neutral mention across
# all three streams), so a reviewer can check the claims without re-running anything.
#
# Run from docs/i169-charge-neutral/ in a directory tree that still holds the full
# logs, i.e. immediately after the runs. Re-running the decks reproduces them.

set -u

emit() {
    local dir="$1" title="$2"
    local out="${dir}/RMG.log.charge-excerpt"
    {
        echo "# ${title}"
        echo "#"
        echo "# EXCERPT of ${dir}/RMG.log -- NOT the whole log. Line numbers are the"
        echo "# original ones. Regenerate with ../make_excerpts.sh after re-running the deck."
        echo "# Full log: $(wc -l < "${dir}/RMG.log") lines."
        echo
        echo "## Every 'Warning:' line in the entire run"
        grep -n "^Warning" "${dir}/RMG.log" || echo "  (none)"
        echo
        echo "## Every charge/neutrality line RMG itself emitted"
        echo "## (lines belonging to the echoed input deck -- the log replays the deck"
        echo "##  verbatim, comments and all -- and lines that only match on the worktree"
        echo "##  path 'i169-charge-neutral' are excluded; that is what makes this a census"
        echo "##  of RMG's own output rather than of the deck's prose.)"
        grep -vn 'i169-charge-neutral' "${dir}/RMG.log" \
            | grep -v '^[0-9]*:[[:space:]]*#' \
            | grep -i 'charge\|neutral' \
            || echo "  (none -- RMG emitted no message mentioning charge or neutrality)"
        echo
        echo "## Same census over stdout.log and stderr.log"
        for stream in stdout.log stderr.log; do
            if [ -f "${dir}/${stream}" ]; then
                printf '%s: %s RMG-emitted charge/neutral mentions\n' "${stream}" \
                    "$(grep -v 'i169-charge-neutral' "${dir}/${stream}" \
                        | grep -v '^[[:space:]]*#' \
                        | grep -ic 'charge\|neutral')"
            fi
        done
        echo
        echo "## Context: the reaction system, and how the run ended"
        sed -n '/^plasmaReactor(/,/^)/p' "${dir}/RMG.log" | sed 's/^/    /'
        echo
        tail -6 "${dir}/RMG.log" | sed 's/^/    /'
    } > "${out}"
    echo "wrote ${out} ($(wc -l < "${out}") lines)"
}

emit before-nonneutral \
    "BEFORE: electronDensity, no cation, unmodified plasma@4f18dc389 -- the silence"
emit after-nonneutral \
    "AFTER: the same deck, new code -- the warning naming the net charge"
emit after-neutral \
    "AFTER: the same deck plus chargeBalanceSpecies='Arp' -- neutral, and silent"
