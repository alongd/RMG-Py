#!/usr/bin/env python3
"""Round 106 HIGH 1 probe -- the departing test must be PER SPECIES, not on the aggregate.

Reproduces the reviewer's sequence: two sub-tolerance neutrals, one ('a') RISING toward its
own transient (log-log slope 5e-8 -> 1e-7 -> 1.5e-7 -> ...), the other ('b') larger and
FALLING (6.77e-7 -> 6.01e-7 -> 5.26e-7 -> ...). The aggregate maximum is 'b' at every step,
and it is falling, so a departing test that reads only the aggregate maximum stops treating
the system as departing and lets the external (electron) arm through -- while a live species
'a' is still climbing toward its timescale. The system is NOT settled and must not arm.

Asserts the fix (per-species trend) holds and EXITS NON-ZERO on failure, so a regression is a
red CI signal, not a line of console output nobody reads. Run against the round-100 built
module this asserts False (the aggregate armed the criterion); against the round-106 module it
passes.
"""
import sys

import numpy as np

sys.path.insert(0, '.')
from rmgpy.solver.termination import TerminationSteadyState

ATOL = 1e-16
BIG = 1e12


def drive(term, a_slopes, b_slopes):
    """Two live sub-tolerance species over a constant BIG neutral, one e-fold per step."""
    a = b = 1.0
    for k in range(len(a_slopes)):
        t_prev, t_now = float(np.exp(k)), float(np.exp(k + 1))       # dlnt = 1
        a_next, b_next = a * np.exp(a_slopes[k]), b * np.exp(b_slopes[k])
        term.update(np.array([BIG, a_next, b_next]), t_now,
                    np.array([BIG, a, b]), t_prev, ATOL,
                    labels=['big', 'a', 'b'], external_armed=True)
        a, b = a_next, b_next


failures = []

# Hidden-rising: 'a' climbs, 'b' (the aggregate max) falls. Must NOT arm.
term = TerminationSteadyState(tolerance=1e-6, window=3)
drive(term,
      [5e-8, 1e-7, 1.5e-7, 2e-7, 2.5e-7, 3e-7, 3.5e-7, 4e-7],
      [6.77e-7, 6.01e-7, 5.26e-7, 5.0e-7, 4.8e-7, 4.7e-7, 4.6e-7, 4.55e-7])
print('hidden-rising:   worst_label={0!r}  armed_external={1}  armed={2}'.format(
    term.worst_label, term.armed_external, term.armed))
if term.worst_label != 'b':
    failures.append('the aggregate maximum should be the falling species b')
if term.armed_external or term.armed:
    failures.append('a rising hidden species must keep the criterion un-armed (HIGH 1 present)')

# Control: both species genuinely settle -> the external arm IS granted.
term2 = TerminationSteadyState(tolerance=1e-6, window=3)
drive(term2,
      [3e-7, 2e-7, 1e-7, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9],
      [6e-7, 4e-7, 2e-7, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9])
print('genuinely settled: armed_external={0}'.format(term2.armed_external))
if not term2.armed_external:
    failures.append('a genuinely settled system must still arm (fix must not over-block)')

if failures:
    print('\nFAIL:')
    for f in failures:
        print('  -', f)
    sys.exit(1)
print('\nPASS: per-species departing blocks the arm the aggregate maximum hid.')
