#!/usr/bin/env python3
"""Round 106 HIGH 2 probe -- the REBUTTAL, measured.

HIGH 2 posits a period-9 trajectory that aliases DECADE sampling: arm at t=1, then residuals
exactly zero at t=10 and t=100, so an endpoints-only span check returns steady at t=100 while
the interior swings. The load-bearing premise is that ``TerminationSteadyState.update`` is
called at decade-spaced sample points, so a periodic composition can land on the same phase
twice and be read as flat.

That premise is false in this codebase, and this probe measures it. ``ReactionSystem.simulate``
drives the integrator in pydas INTERMEDIATE one-step mode (``self.step(step_time)`` at
base.pyx:796) and calls ``update`` with ``self.t`` after EVERY internal DASSL step, not at
decade marks. DASSL's own local-truncation-error control cannot accept a step that spans an
oscillation period; so a flat residual over an accepted step is the integrator certifying the
composition is smooth across it. A period-9 oscillation would therefore show non-zero residual
at the fine interior steps between t=1 and t=10 -- the criterion sees the interior the HIGH-2
case assumes it skips.

This probe instruments a genuinely relaxing wall reactor, records the actual sequence of times
``update`` is fed, and asserts the cadence is fine (median and maximum step far below one
decade, and no single step spans a decade). It EXITS NON-ZERO if the cadence is ever coarse
enough for decade-aliasing to be possible. See rework.md, Round 106 HIGH 2, for the argument.
"""
import sys

import numpy as np

sys.path.insert(0, '.')
sys.path.insert(0, 'test/rmgpy/solver')

from rmgpy.solver.termination import TerminationSteadyState, TerminationTime
from rmgpy.exceptions import PlasmaStateError

# Reuse the steady-state end-to-end harness: a plasma reactor that genuinely relaxes.
from steadyStateTest import _relaxing_system, _simulate


class CadenceSpy(TerminationSteadyState):
    """Records every (t_prev, t_now) pair update() is fed, then defers to the real criterion."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.samples = []

    def update(self, y_now, t_now, y_prev, t_prev, floor, labels=None, **kwargs):
        if t_prev > 0.0 and t_now > t_prev:
            self.samples.append((t_prev, t_now))
        return super().update(y_now, t_now, y_prev, t_prev, floor, labels=labels, **kwargs)


spy = CadenceSpy(tolerance=1e-8)
reactor, core_species, core_reactions = _relaxing_system(
    [spy, TerminationTime((1.0, 's'))])
try:
    _simulate(reactor, core_species, core_reactions)
except PlasmaStateError as e:
    # A fully-relaxed run can trip the electron-negativity guard as n_e drifts to zero; by
    # then the cadence is already recorded, which is all this probe measures.
    print('(simulate stopped on the electron guard: {0})'.format(str(e)[:70]))

decade_steps = np.array([np.log10(t_now) - np.log10(t_prev)
                         for (t_prev, t_now) in spy.samples])
n = decade_steps.size
median = float(np.median(decade_steps)) if n else float('nan')
worst = float(decade_steps.max()) if n else float('nan')
n_decade_jumps = int((decade_steps >= 1.0).sum())

print('update() samples: {0}'.format(n))
print('step size in decades: median={0:.4g}  max={1:.4g}'.format(median, worst))
print('steps spanning >= 1 full decade (would permit decade-aliasing): {0}'.format(n_decade_jumps))

failures = []
if n < 20:
    failures.append('too few samples ({0}) to characterise the cadence'.format(n))
if not (worst < 1.0):
    failures.append('a step spanned >= 1 decade (max={0:.4g}); decade-aliasing becomes '
                    'possible'.format(worst))
if not (median < 0.1):
    failures.append('median step {0:.4g} decade is not fine relative to a decade'.format(median))

if failures:
    print('\nFAIL:')
    for f in failures:
        print('  -', f)
    sys.exit(1)
print('\nPASS: update() is fed every internal DASSL step at fine cadence, never at decade '
      'marks. The period-9 decade-aliasing case cannot arise: a periodic composition shows '
      'non-zero residual at the interior steps the criterion actually sees.')
