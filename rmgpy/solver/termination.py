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

    1. armed -- ``R >= 1`` has held at least once (or, for a sub-floor channel the reactor
       tracks externally, its own relaxation time has passed while no generic channel is
       still departing -- see :meth:`update`); and
    2. flat -- ``R < tolerance`` and the flat run has persisted for at least one SYSTEM
       RELAXATION TIME. Persistence is a physical span, not a step count: counting accepted
       solver steps would make the verdict depend on the integrator's step controller, so
       the same plateau would terminate under fine stepping and not under coarse. The
       reactor supplies the relaxation time (``1/nu_wall`` for a wall-bounded discharge);
       an ordinary reactor that knows no such time falls back to one e-fold of absolute
       time. The step ``window`` survives only as a cheap fluke guard -- a flat interval
       must be at least two accepted samples -- and is never sufficient on its own.

    **A system that carries no flux at all is handled by the solver, and it is NOT a steady
    state.** Its composition cannot change again, so the integration has to stop -- but a
    system that never started has converged to nothing, and reporting that as a satisfied
    criterion is precisely the failure this criterion exists to remove. The arm is what
    tells the two apart: a run that went through a transient and then froze arrives here
    armed and flat and counts (after one relaxation time confirms it, like any other tail);
    a run that never armed is terminated by the solver with ``steady_state_reached`` left
    False and a warning saying no steady state was demonstrated. There is no special case
    for a residual of exactly zero -- equal endpoints do not prove a frozen structure, so a
    zero residual earns the same one-relaxation-time confirmation as any other flat tail.

    Attributes:

    `tolerance`     the residual below which the composition counts as no longer changing
    `window`        the fluke-guard floor on how many accepted samples a flat interval spans
    `armed`         whether the integration has passed its fastest relaxation time
    `streak`        consecutive flat steps satisfied so far
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
        self.armed_generic = False
        self.armed_external = False
        self._r_gen_prev = float('nan')
        self._steps_since_generic_rise = 0
        self.streak = 0
        self._t_flat_start = float('nan')
        self.residual = float('nan')
        self.worst_label = None

    def update(self, y_now, t_now, y_prev, t_prev, floor, labels=None,
               external_residual=float('nan'), external_armed=False,
               relaxation_time=float('nan')):
        """
        Fold one solver step into the criterion and report whether it is now satisfied.

        `y_now`/`y_prev` are the core-species mole vectors at `t_now`/`t_prev`; `floor`
        is the moles below which a species is beneath the integrator's resolution (pass
        the simulator's ``atol``). The residual itself is taken on **mole fractions** --
        the DSL's own currency, and immune to the drift in total moles that a
        two-temperature equation of state produces without any chemistry happening.

        `external_residual`/`external_armed` carry a channel this criterion cannot see on
        its own: a state variable the reactor resolves BELOW `floor`. A plasma electron
        seeded from zero by an external source saturates far under ``atol``, so
        :meth:`compute_residual` (which excludes sub-floor species) reads only the inert
        neutrals and would fire at ignition instead of at saturation. The reactor passes the
        electron's own **mole-fraction** slope as `external_residual` -- the same intensive
        quantity every other residual measures -- folded into the residual by MAX so firing
        waits for it to go flat, and `external_armed` to arm the criterion from the reactor's
        known relaxation time, since a saturating-from-zero slope is bounded by 1 and never
        reaches the generic ``R >= 1`` arm. Both default to nan/False, so an ordinary reactor
        is byte-for-byte unaffected: with no external channel, arming reduces to the generic
        ``R >= 1`` latch exactly as before.

        `external_armed` vouches ONLY for the electron. It may arm the whole criterion only
        while the generic (neutral) residual is not still DEPARTING -- and departing is
        judged over a sequence (no step-over-step rise for `window` consecutive samples),
        not from the two most recent values, so one flat or noisy sample cannot license it.
        That arm is re-evaluated every step and does not latch: a neutral mid-transit --
        below tolerance only because it has not yet reached its timescale -- keeps the system
        unsteady even after the electron has saturated, and a settled neutral that resumes
        moving un-arms it. See the per-channel arming block below.

        `relaxation_time` is the system timescale the flat run must persist across before the
        state counts as steady (the reactor's ``1/nu_wall`` for a wall-bounded discharge).
        Anchoring persistence to a physical time rather than to a count of accepted steps is
        what makes the verdict independent of the integrator's step controller. Defaults to
        nan, in which case persistence falls back to one e-fold of absolute time.

        Returns True only if armed and the flat condition has held for `window` steps.
        """
        r_gen, self.worst_label = self.compute_residual(y_now, t_now, y_prev, t_prev, floor,
                                                        labels=labels)
        # Fold in the reactor's sub-floor channel by MAX for the FLAT test: the composition
        # is settled only when the slowest of everything -- neutrals AND the invisible
        # electron -- is flat.
        r = r_gen
        if np.isfinite(external_residual):
            if not np.isfinite(r) or external_residual > r:
                r = external_residual
                self.worst_label = '<external channel>'
        self.residual = r

        # Arming is PER CHANNEL -- the thing that arms must be the thing declared steady.
        #  - the GENERIC channel arms on its OWN residual crossing 1 (a neutral has run
        #    through its relaxation time). That is a HISTORICAL FACT about the trajectory --
        #    it either happened or it did not -- so this latch is permanent and sound.
        #  - the EXTERNAL (electron) channel is licensed to arm the whole criterion only
        #    while NO generic channel is still DEPARTING (climbing toward its own arm). A
        #    saturating-from-zero electron slope is bounded by 1 and never reaches the
        #    ``R >= 1`` arm, so the reactor vouches for it from its known relaxation time;
        #    but that vouch says nothing about a neutral still in transit, so it may license
        #    the system only while the generic channel is not departing. "Not departing" is
        #    a LIVE condition, not a historical fact -- a neutral can start moving again --
        #    so this arm is RE-EVALUATED every step and NOT latched. An irreversible arm on a
        #    reversible condition is exactly the defect this replaces.
        #
        # "Departing" is a property of a SEQUENCE, not of the two most recent samples. The
        # old test compared r_gen to the single preceding value, so one flat or noisy sample
        # -- or the very first step, whose predecessor is nan -- read as "not rising" and
        # armed the channel forever. Here the generic residual counts as still climbing until
        # it has failed to set a step-over-step rise for `window` CONSECUTIVE finite samples:
        # a single flat/noisy sample inside a genuine climb is undone by the next rise and
        # can no longer license the arm, while a residual that has truly stopped rising for a
        # full window releases it. Threshold-free -- it reads the sign of the trend over the
        # window, never a magnitude floor.
        if np.isfinite(r_gen):
            if np.isfinite(self._r_gen_prev) and r_gen > self._r_gen_prev:
                self._steps_since_generic_rise = 0    # rose vs the previous sample: departing
            else:
                self._steps_since_generic_rise += 1   # flat or falling
            self._r_gen_prev = r_gen
        if np.isfinite(r_gen) and r_gen >= 1.0:
            self.armed_generic = True
        generic_departing = self._steps_since_generic_rise < self.window
        self.armed_external = external_armed and not generic_departing
        self.armed = self.armed_generic or self.armed_external

        if not np.isfinite(r):
            # Not evaluable (first step, zero/degenerate time interval, or no live
            # species yet): no information either way, so do not advance the streak and
            # do not arm. `inf` lands here too, which is correct -- it is emitted when a
            # species crosses up through the floor, which is the opposite of steady.
            self.streak = 0
            return False
        if r < self.tolerance:
            if self.streak == 0:
                # Anchor the flat interval at its FIRST endpoint (t_prev of the first flat
                # step), not its current endpoint: the residual at this step describes the
                # interval [t_prev, t_now], so the flat run began at t_prev. Setting the
                # start to t_now instead threw away that first interval and made a late tail
                # need an extra span measured from where the window happened to open.
                self._t_flat_start = t_prev
            self.streak += 1
        else:
            self.streak = 0
            self._t_flat_start = float('nan')

        # Persistence is a PHYSICAL span, not a step count. Counting accepted solver steps
        # makes the verdict a function of the integrator's controller: the same plateau
        # terminates under fine stepping (many steps in the window) and not under coarse
        # stepping (too few), though the physics is identical. So the binding requirement is
        # that the flat run has persisted for at least one SYSTEM RELAXATION TIME -- the
        # timescale on which this system settles, which the reactor supplies
        # (`relaxation_time`, e.g. ``1/nu_wall`` for a wall-bounded discharge). That is
        # step-controller-independent AND anchored to the physics rather than to absolute
        # elapsed time: a late-converging tail needs one more tau to confirm, never a fixed
        # multiple of whenever the window opened. Where the reactor knows no system time (an
        # ordinary reactor passes nan), fall back to one e-fold of absolute time -- the
        # native scale of a ``d/d ln t`` criterion.
        #
        # The step `window` is kept ONLY as a cheap fluke guard -- a flat interval must be at
        # least two accepted samples, so a single flat step that happens to span a whole tau
        # cannot alone terminate -- and is deliberately NOT sufficient: the physical span
        # must also hold. Two is the floor because a span needs two endpoints; requiring more
        # would put back the step-count dependence this removes. There is no exact-zero
        # waiver: equal endpoints do not prove a frozen structure (they alias an oscillation
        # or a stop-restart), so a residual of zero earns the same one-relaxation-time
        # confirmation as any other flat tail. Anchoring to tau rather than to an absolute
        # e-fold keeps that confirmation short enough that a fully-pumped (gamma=0) discharge
        # settling to n_e -> 0 is recognised before the electron drifts past the wall guard.
        if not (self.streak >= 2 and np.isfinite(self._t_flat_start)
                and self._t_flat_start > 0.0 and t_now > 0.0):
            return False
        if np.isfinite(relaxation_time) and relaxation_time > 0.0:
            span_ok = (t_now - self._t_flat_start) >= relaxation_time
        else:
            span_ok = (np.log(t_now) - np.log(self._t_flat_start)) >= 1.0
        return self.armed and span_ok

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