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

    **What ``window`` buys, and why one number governs two things.** ``window`` sets both (a)
    the flat-streak floor -- how many consecutive flat samples a tail must span before it can
    terminate -- and (b) the departing test's patience -- how many consecutive samples a species
    must fail to rise before it stops counting as still climbing (see :meth:`update`). These are
    NOT opposed knobs pulling in different directions: raising ``window`` makes BOTH strictly
    more conservative. A longer streak is harder to satisfy, so persistence tightens; and a
    species is held as "departing" for more samples, so the external channel is licensed to arm
    later. Both effects push termination LATER, never earlier -- a larger ``window`` can only
    delay a verdict, monotonically (pinned by
    ``test_raising_window_makes_the_criterion_no_less_conservative_in_both_roles``). So the
    single number reads coherently as one quantity: *how much evidence to require before trusting
    a "no longer changing" verdict*. Splitting it into two knobs would let a user buy a stricter
    streak while weakening the departing test, which is not a combination that means anything
    physical here.

    Attributes:

    `tolerance`     the residual below which the composition counts as no longer changing
    `window`        evidence floor (samples): both the flat-streak length AND the departing
                    test's consecutive-non-rise count; larger is strictly more conservative
    `armed`         whether the integration has passed its fastest relaxation time
    `streak`        consecutive flat steps satisfied so far
    `residual`      the most recently evaluated residual (nan before the second step)
    `worst_label`   the species carrying that residual, for the log
    """

    def __init__(self, tolerance=1e-6, window=2):
        tolerance = float(tolerance)
        if not np.isfinite(tolerance) or tolerance <= 0.0:
            raise ValueError('terminationSteadyState tolerance must be finite and strictly '
                             'positive; got {0!r}.'.format(tolerance))
        window = int(window)
        if window < 2:
            raise ValueError('terminationSteadyState window must be at least 2 accepted '
                             'samples: a flat interval needs two endpoints, and a single '
                             'step -- however long it spans -- is not an interval. Got '
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
        # HIGH 1 (round 106): the departing test is PER SPECIES, not on the aggregate
        # maximum. A dict of the previous log-log slope for each live species, and a dict
        # counting how many consecutive samples each has failed to rise. A species is still
        # departing until it has not risen for `window` consecutive samples; the generic
        # channel is departing while ANY live species still is.
        self._slope_prev = {}
        self._steps_since_rise = {}
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
        analysis = self._slope_analysis(y_now, t_now, y_prev, t_prev, floor, labels=labels)
        r_gen = analysis['r']
        self.worst_label = analysis['label']
        # Fold in the reactor's sub-floor channel by MAX for the FLAT test: the composition
        # is settled only when the slowest of everything -- neutrals AND the invisible
        # electron -- is flat.
        r = r_gen
        # A non-finite external residual must be closed as a CLASS, not member by member
        # (round 111 HIGH 1). The dispatch is ``np.isfinite`` -- whose negation is True for
        # ``nan`` BY CONSTRUCTION -- and never a comparison or a maximum, which ``nan`` slips
        # through silently (``nan > r``, ``nan < r``, ``nan == r`` are all False). Round 109
        # made the fold ignore a non-finite value; round 110 poisoned +inf/-inf but its own
        # underflow-addendum rerouted the subnormal-fraction case from +inf onto ``nan``, which
        # then took the benign "no information" branch and let a flat generic channel fire on an
        # armed-but-unusable electron channel. So the branch below cannot be written by listing
        # the values checked; it must make NO non-finite value reach the flat test as a usable
        # number.
        if np.isfinite(external_residual):
            # A usable finite external residual: MAX-combine for the flat test. Substitute only
            # when the generic channel carried NO information (``nan``) or is the smaller of the
            # two. A generic ``inf`` is preserved, since ``external_residual > inf`` is False.
            if np.isnan(r) or external_residual > r:
                r = external_residual
                self.worst_label = '<external channel>'
        elif bool(external_armed) or not np.isnan(external_residual):
            # A NON-FINITE external residual (nan, +inf or -inf) from a channel that is IN PLAY:
            # it has ARMED -- a live discharge whose electron passed its relaxation time
            # ``t*nu_wall >= 1`` -- or it is an actual +/-inf, which only an active channel can
            # compute (an ordinary reactor emits only the default ``nan``). The channel is active
            # yet cannot certify a steady electron, so it must never authorise termination through
            # ANY channel: not the external arm, and not a flat generic channel left to fire on
            # its own. POISON the folded residual with ``inf`` so the non-finite guard below resets
            # the streak and it can never arm, and NAME the invalid value in diagnostics. Testing
            # ``not np.isnan`` here (True only for the infinities) keeps the +/-inf poison
            # UNCONDITIONAL -- as in round 110 -- while ``external_armed`` adds the armed-``nan``
            # case round 110 left open.
            r = float('inf')
            self.worst_label = '<external channel: non-finite residual {0!r}>'.format(external_residual)
        else:
            # A ``nan`` from a channel that is NOT in play: the ordinary reactor (``external_armed``
            # False, ``external_residual`` the default ``nan``), or a plasma channel that has not
            # yet armed and merely reported no number. Genuinely NO information -- keep the generic
            # residual and judge the generic channel on its own. Byte-for-byte the pre-plasma path.
            pass
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
        # The generic R>=1 arm is a HISTORICAL FACT on the aggregate maximum -- did the
        # fastest chemistry ever run through its own relaxation time. The maximum is exactly
        # right for that ("did ANY species reach R=1"), and it either happened or it did not,
        # so this latch is permanent and sound.
        if np.isfinite(r_gen) and r_gen >= 1.0:
            self.armed_generic = True

        # "Departing" is PER SPECIES, and a property of each species' own SEQUENCE (HIGH 1,
        # round 106). The previous test asked whether the aggregate MAXIMUM residual rose;
        # that hides a species climbing toward its own transient whenever a different,
        # decaying species holds the maximum, and it compares the slopes of two DIFFERENT
        # species across a step in which the largest changed. "Is anything still moving" is a
        # per-species question. Track each live species' log-slope against its OWN previous
        # value: it counts as still climbing until it has failed to rise for `window`
        # CONSECUTIVE samples. A single flat/noisy sample inside a genuine climb is undone by
        # the next rise; a species that has truly stopped rising for a full window releases
        # its hold. The generic channel is departing while ANY live species still is.
        # Threshold-free -- it reads the sign of each species' trend, never a magnitude floor.
        # `nan`, `inf` and `ok` are THREE different answers, and the counters do a different
        # thing on each -- decided here explicitly (round 109 HIGH, second half):
        if analysis['status'] == 'ok':
            # A real per-step slope for each live species: count each species' trend against
            # its OWN previous sample, and prune species no longer live so a stale, frozen
            # counter cannot vouch that the present composition has stopped moving.
            live_keys = set(analysis['slopes'])
            for key, slope in analysis['slopes'].items():
                prev = self._slope_prev.get(key)
                if prev is not None and slope > prev:
                    self._steps_since_rise[key] = 0                       # rose vs its own last sample
                else:
                    self._steps_since_rise[key] = self._steps_since_rise.get(key, 0) + 1
                self._slope_prev[key] = slope
            for key in [k for k in self._steps_since_rise if k not in live_keys]:
                self._steps_since_rise.pop(key, None)
                self._slope_prev.pop(key, None)
        elif analysis['status'] == 'inf':
            # STRONG NEGATIVE information: a species appeared from nothing or went negative --
            # a STRUCTURAL change to the composition, emphatically not steady. The per-species
            # trend history described a composition that no longer holds; leaving it in place
            # would let a stale "settled" count from before the change keep `generic_departing`
            # False, so the very step that proves the system is not steady could arm and
            # terminate on a substituted external residual (that was the round 109 defect).
            # DISCARD the history: with no counters the generic channel reads as departing
            # (below), `armed_external` drops, and the appearance cannot be mistaken for a
            # settled tail. A genuine re-settling afterward rebuilds the counters over a fresh
            # window -- the confirmation the appearance interrupted is paid again, correctly.
            self._steps_since_rise.clear()
            self._slope_prev.clear()
        # else `status == 'nan'`: NO information (degenerate interval, dead totals, no live
        # species). The step says nothing about whether anything is moving, so the per-species
        # counters are left EXACTLY as they are -- neither advanced nor discarded. (The flat
        # STREAK is still broken by the non-finite guard below; only the trend memory persists,
        # because we did not observe any species change.)
        if not self._steps_since_rise:
            generic_departing = True     # nothing measured yet -- conservatively still departing
        else:
            generic_departing = any(c < self.window for c in self._steps_since_rise.values())
        # The external arm and the external RESIDUAL must travel together (round 110 HIGH 1).
        # Round 109 made the FOLD above ignore a non-finite `external_residual`, but the boolean
        # `external_armed` arrives on this same call and was read on its own -- so a channel that
        # had just reported NO usable number (`nan`) or a value that blew up (`inf`, e.g. a
        # subnormal electron fraction underflowing to zero in the hook) still granted permission,
        # and a flat generic channel then terminated on it. Ignoring a value is not the same as
        # distrusting its source: `external_armed` vouches for the electron, but if the electron
        # channel produced no usable number it cannot vouch for anything, so the arm is licensed
        # only while `external_residual` is finite. An ordinary reactor passes `external_armed`
        # False and `external_residual` nan, so this is a no-op there and arming reduces to the
        # generic R>=1 latch exactly as before.
        #
        # ``bool(...)`` is not cosmetic: ``np.isfinite`` returns a ``numpy.bool_``, and Python's
        # short-circuiting ``and`` returns that scalar when it is the deciding operand -- so
        # without the coercion ``armed_external`` (and, through the ``or``, ``armed``) silently
        # became a ``numpy.bool_``, which is truthy/falsy but fails an ``is False`` identity test
        # a caller may make (round 110: the same "almost right on a contract that says exactly
        # right" class as HIGH 1). These attributes cross into the compiled solver, whose declared
        # type here is a Python ``bool``; keep them that.
        self.armed_external = bool(external_armed and np.isfinite(external_residual)
                                   and not generic_departing)
        self.armed = bool(self.armed_generic or self.armed_external)

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
        # The step `window` is a fluke-guard floor on the sample count -- a flat interval
        # must span at least `window` accepted samples -- and it is HONOURED as declared
        # (round 106 MEDIUM): the previous code hard-coded a floor of two regardless of
        # `window`, so a knob the user set was silently ignored. It defaults to two, the
        # irreducible floor (a span needs two endpoints; a single flat step, however long,
        # is not an interval), and at the default the verdict is step-controller-independent
        # because any positive span yields at least two samples. Raising `window` above two
        # trades that independence for extra insurance against a fluke, at the user's explicit
        # choice -- but it is never SUFFICIENT: the physical span below must also hold, so the
        # physics, not the sample count, remains the control. There is no exact-zero waiver:
        # equal endpoints do not prove a frozen structure (they alias an oscillation or a
        # stop-restart), so a residual of zero earns the same one-relaxation-time confirmation
        # as any other flat tail. Anchoring to tau rather than to an absolute e-fold keeps
        # that confirmation short enough that a fully-pumped (gamma=0) discharge settling to
        # n_e -> 0 is recognised before the electron drifts past the wall guard.
        if not (self.streak >= self.window and np.isfinite(self._t_flat_start)
                and self._t_flat_start > 0.0 and t_now > 0.0):
            return False
        if np.isfinite(relaxation_time) and relaxation_time > 0.0:
            span_ok = (t_now - self._t_flat_start) >= relaxation_time
        else:
            span_ok = (np.log(t_now) - np.log(self._t_flat_start)) >= 1.0
        # ``bool(...)`` on the final verdict, not just the arm attributes (round 111 MEDIUM):
        # the e-fold fallback computes ``span_ok`` from ``np.log`` -- a numpy.float64 comparison
        # yielding a ``numpy.bool_`` -- so ``self.armed and span_ok`` returns that scalar when
        # armed is True, and update() would leak a ``numpy.bool_`` an ``is True`` caller fails.
        # This is the single chokepoint every return-True path flows through.
        return bool(self.armed and span_ok)

    @staticmethod
    def _slope_analysis(y_now, t_now, y_prev, t_prev, floor, labels=None):
        """
        The full per-step slope analysis, from which both the aggregate residual and the
        per-species departing test are derived (HIGH 1, round 106).

        Returns a dict with:

        - ``status``: ``'nan'`` (no information -- degenerate interval, dead totals, no
          live species), ``'inf'`` (a species went negative beyond the floor or appeared
          from nothing -- emphatically not steady), or ``'ok'`` (a real aggregate slope).
        - ``r``: the aggregate residual -- ``nan``/``inf``/``max_i |slope_i|`` matching
          ``status``.
        - ``label``: the caller's name (or bare index) for the species carrying ``r``, or
          ``None`` when there is none.
        - ``slopes``: ``{int core-species index -> float log-log slope}`` for every species
          live at BOTH endpoints with a finite slope, present only when ``status == 'ok'``
          (empty otherwise). Keyed by INTEGER INDEX -- guaranteed hashable, and stable
          across steps as the set of live species changes -- so :meth:`update` can track each
          species' own trend rather than the aggregate maximum. A species that HELD the
          maximum on one step and yields it on the next is one falling entry here; a species
          climbing toward its own transient while a different, decaying species holds the
          maximum is one rising entry the aggregate could never show.
        """
        empty = {}
        y_now = np.asarray(y_now, dtype=np.float64)
        y_prev = np.asarray(y_prev, dtype=np.float64)
        if t_prev <= 0.0 or t_now <= t_prev:
            return {'status': 'nan', 'r': float('nan'), 'label': None, 'slopes': empty}
        dlnt = np.log(t_now) - np.log(t_prev)
        if not np.isfinite(dlnt) or dlnt <= 0.0:
            return {'status': 'nan', 'r': float('nan'), 'label': None, 'slopes': empty}

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
            return {'status': 'inf', 'r': float('inf'), 'label': _label_of(labels, idx), 'slopes': empty}

        total_now = y_now.sum()
        total_prev = y_prev.sum()
        if not (np.isfinite(total_now) and np.isfinite(total_prev)) or total_now <= 0.0 or total_prev <= 0.0:
            return {'status': 'nan', 'r': float('nan'), 'label': None, 'slopes': empty}

        # The floor is applied to MOLES (that is where the integrator's absolute tolerance
        # lives), while the residual is taken on mole fractions.
        live_now = y_now > floor
        live_prev = y_prev > floor

        appeared = live_now & ~live_prev
        if appeared.any():
            idx = int(np.argmax(appeared))
            return {'status': 'inf', 'r': float('inf'), 'label': _label_of(labels, idx), 'slopes': empty}

        both = live_now & live_prev
        if not both.any():
            return {'status': 'nan', 'r': float('nan'), 'label': None, 'slopes': empty}

        both_indices = np.flatnonzero(both)
        x_now = y_now[both] / total_now
        x_prev = y_prev[both] / total_prev
        slope = np.abs(np.log(x_now) - np.log(x_prev)) / dlnt
        finite = np.isfinite(slope)
        if not finite.any():
            return {'status': 'nan', 'r': float('nan'), 'label': None, 'slopes': empty}
        # Per-species slopes, keyed by the core-species integer index (hashable and stable),
        # for every both-live species with a finite slope. This is what the departing test
        # tracks per species; the aggregate max below is only what the generic R>=1 arm reads.
        slopes = {int(both_indices[j]): float(slope[j]) for j in range(slope.shape[0]) if finite[j]}
        slope_masked = np.where(finite, slope, -1.0)
        k = int(np.argmax(slope_masked))
        idx = int(both_indices[k])
        return {'status': 'ok', 'r': float(slope_masked[k]), 'label': _label_of(labels, idx), 'slopes': slopes}

    @staticmethod
    def compute_residual(y_now, t_now, y_prev, t_prev, floor, labels=None):
        """
        The residual R for one step, and the label (or index) of the species carrying it.

        A thin wrapper over :meth:`_slope_analysis` that keeps the historical 2-tuple
        contract (``(r, label)``) for the many callers and tests that read only the
        aggregate. Returns ``(nan, None)`` when the step carries no information: the
        interval is not a positive interval in log time, or no species is live at both
        ends. Returns ``(inf, label)`` when a species has crossed *up* through `floor`,
        i.e. appeared -- a system growing a new species is emphatically not stationary.
        """
        analysis = TerminationSteadyState._slope_analysis(
            y_now, t_now, y_prev, t_prev, floor, labels=labels)
        return analysis['r'], analysis['label']


def _label_of(labels, index):
    """The caller's name for core species `index`, or the bare index if it gave none."""
    if labels is None:
        return index
    try:
        return labels[index]
    except (IndexError, TypeError):
        return index