#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import numpy as np

from rmgpy.quantity import Quantity
"""
Contains classes for termination criteria
"""
class TerminationTime:
    """
    Represent a time at which the simulation should be terminated. This class
    has one attribute: the termination `time` in seconds.
    """

    def __init__(self, time=(0.0, 's')):
        self.time = Quantity(time)


################################################################################

class TerminationConversion:
    """
    Represent a conversion at which the simulation should be terminated. This
    class has two attributes: the `species` to monitor and the fractional
    `conversion` at which to terminate.
    """

    def __init__(self, spec=None, conv=0.0):
        self.species = spec
        self.conversion = conv


class TerminationRateRatio:
    """
    Represent a fraction of the maximum characteristic rate of the simulation
    at which the simulation should be terminated.  This class has one attribute
    the ratio between the current and maximum characteristic rates at which
    to terminate
    """

    def __init__(self, ratio=0.01):
        self.ratio = ratio


################################################################################

class TerminationSteadyState:
    """
    Represent a *stationary state* at which the simulation should be terminated: the
    endpoint of a low-pressure discharge, where the composition has stopped changing
    while nothing is being consumed net and no rate ratio is decaying. None of the other
    three criteria can express that.

    THE CRITERION IS A PROPOSAL. What follows is argued from a measured trajectory (see
    ``docs/i170-steady-state/FINDING.md``), not from first principles, but the quantity,
    the tolerance and the window are physics choices that have not yet been ratified.

    **Residual.** For every core species that is *live* -- resolved by the integrator,
    i.e. holding more moles than the simulator's own ``atol`` -- take the discrete
    log-log slope of its mole fraction between two consecutive solver steps, and take
    the largest:

    .. math::
        R = \\max_i \\left| \\frac{\\ln x_i(t_k) - \\ln x_i(t_{k-1})}
                                 {\\ln t_k - \\ln t_{k-1}} \\right|

    ``R`` is dimensionless, independent of the time unit, and independent of the
    solver's step size -- which matters more than it looks. It also has a direct
    physical reading: since :math:`|d \\ln x_i / dt| = 1 / \\tau_i` is the inverse of
    species *i*'s relaxation time, :math:`R = t / \\tau_{\\mathrm{min}}` is **the elapsed
    integration time expressed in units of the fastest chemistry still running**.

    **Why the residual alone is not enough, and why this class latches.** Because
    ``R = t / tau_min``, a small ``R`` has *two* causes: the system has converged, or the
    integration has not yet reached its own chemistry. Both are real. Measured on the
    lithium design case, a criterion without the latch below terminates at
    ``t = 3.0e-15 s`` -- step 2 of 497 -- on a system whose true endpoint is
    ``1.84e-4 s``, too early by a factor of 6e10. No tolerance fixes this: to exclude the
    early plateau you would need a tolerance below 2.8e-10, which is beneath the
    integrator's own ``rtol`` and therefore beneath the noise. So the criterion **arms**:
    it may only fire once ``R`` has reached 1, i.e. once the integration has actually run
    for longer than the fastest relaxation time in the system.

    **Termination.** Both of:

    1. armed -- ``R >= 1`` has held at least once; and
    2. flat -- ``R < tolerance`` on ``window`` consecutive steps.

    **A system that carries no flux at all is handled by the solver, and it is NOT a steady
    state.** Its composition cannot change again, so the integration has to stop -- but a
    system that never started has converged to nothing, and reporting that as a satisfied
    criterion is precisely the failure this criterion exists to remove. The latch is what
    tells the two apart, with no extra machinery: a run that went through a transient and
    then froze arrives here armed, with a residual of exactly zero that clears any
    tolerance and fills the window like any other flat tail, and it counts. A run that
    never armed is terminated by the solver with ``steady_state_reached`` left False and a
    warning saying no steady state was demonstrated.

    Attributes:

    `tolerance`     the residual below which the composition counts as no longer changing
    `window`        how many consecutive steps must satisfy it
    `armed`         whether the integration has passed its fastest relaxation time
    `streak`        consecutive steps satisfied so far
    `residual`      the most recently evaluated residual (nan before the second step)
    `worst_label`   the species carrying that residual, for the log
    """

    def __init__(self, tolerance=1e-6, window=3):
        tolerance = float(tolerance)
        if not np.isfinite(tolerance) or tolerance <= 0.0:
            raise ValueError('terminationSteadyState tolerance must be finite and strictly '
                             'positive; got {0!r}.'.format(tolerance))
        window = int(window)
        if window < 1:
            raise ValueError('terminationSteadyState window must be at least 1 step; got '
                             '{0!r}.'.format(window))
        self.tolerance = tolerance
        self.window = window
        self.reset()

    def reset(self):
        """
        Forget the evaluation state. Every call to ``simulate`` restarts the integration
        from t=0, so the latch must not carry over from the previous one -- a termination
        object outlives one simulation and is reused across the whole model-generation run.
        """
        self.armed = False
        self.streak = 0
        self.residual = float('nan')
        self.worst_label = None

    def update(self, y_now, t_now, y_prev, t_prev, floor, labels=None):
        """
        Fold one solver step into the criterion and report whether it is now satisfied.

        `y_now`/`y_prev` are the core-species mole vectors at `t_now`/`t_prev`; `floor`
        is the moles below which a species is beneath the integrator's resolution (pass
        the simulator's ``atol``). The residual itself is taken on **mole fractions** --
        the DSL's own currency, and immune to the drift in total moles that a
        two-temperature equation of state produces without any chemistry happening.

        Returns True only if armed and the flat condition has held for `window` steps.
        """
        r, self.worst_label = self.compute_residual(y_now, t_now, y_prev, t_prev, floor,
                                                    labels=labels)
        self.residual = r
        if not np.isfinite(r):
            # Not evaluable (first step, zero/degenerate time interval, or no live
            # species yet): no information either way, so do not advance the streak and
            # do not arm. `inf` lands here too, which is correct -- it is emitted when a
            # species crosses up through the floor, which is the opposite of steady.
            self.streak = 0
            return False
        if r >= 1.0:
            self.armed = True
        if r < self.tolerance:
            self.streak += 1
        else:
            self.streak = 0
        return self.armed and self.streak >= self.window

    @staticmethod
    def compute_residual(y_now, t_now, y_prev, t_prev, floor, labels=None):
        """
        The residual R for one step, and the label (or index) of the species carrying it.

        Returns ``(nan, None)`` when the step carries no information: the interval is not
        a positive interval in log time, or no species is live at both ends. Returns
        ``(inf, label)`` when a species has crossed *up* through `floor`, i.e. appeared --
        a system growing a new species is emphatically not stationary.
        """
        y_now = np.asarray(y_now, dtype=np.float64)
        y_prev = np.asarray(y_prev, dtype=np.float64)
        if t_prev <= 0.0 or t_now <= t_prev:
            return float('nan'), None
        dlnt = np.log(t_now) - np.log(t_prev)
        if not np.isfinite(dlnt) or dlnt <= 0.0:
            return float('nan'), None

        # A negative population is a solver failure, not a settled composition. The live
        # mask below tests ``y > floor``, which would silently DROP a large-negative species
        # -- neither counting it in the residual nor reporting it -- so the criterion could
        # declare steady state while a species sits large and negative. A moles value below
        # ``-floor`` is negative beyond the integrator's own resolution, i.e. not round-off
        # noise around zero. Poison the residual instead: such a species is emphatically not
        # stationary, so return inf (which resets the streak and can never arm) and surface it
        # in the label, mirroring how an *appearing* species is handled below.
        negative = (y_now < -floor) | (y_prev < -floor)
        if negative.any():
            idx = int(np.argmax(negative))
            return float('inf'), _label_of(labels, idx)

        total_now = y_now.sum()
        total_prev = y_prev.sum()
        if not (np.isfinite(total_now) and np.isfinite(total_prev)) or total_now <= 0.0 or total_prev <= 0.0:
            return float('nan'), None

        # The floor is applied to MOLES (that is where the integrator's absolute tolerance
        # lives), while the residual is taken on mole fractions.
        live_now = y_now > floor
        live_prev = y_prev > floor

        appeared = live_now & ~live_prev
        if appeared.any():
            idx = int(np.argmax(appeared))
            return float('inf'), _label_of(labels, idx)

        both = live_now & live_prev
        if not both.any():
            return float('nan'), None

        x_now = y_now[both] / total_now
        x_prev = y_prev[both] / total_prev
        slope = np.abs(np.log(x_now) - np.log(x_prev)) / dlnt
        if not np.isfinite(slope).any():
            return float('nan'), None
        slope = np.where(np.isfinite(slope), slope, -1.0)
        k = int(np.argmax(slope))
        idx = int(np.flatnonzero(both)[k])
        return float(slope[k]), _label_of(labels, idx)


def _label_of(labels, index):
    """The caller's name for core species `index`, or the bare index if it gave none."""
    if labels is None:
        return index
    try:
        return labels[index]
    except (IndexError, TypeError):
        return index