#!/usr/bin/env python3
"""
I-246 PROBE D -- can the solver RMG actually links against integrate an
index-1 DAE row of the shape the quasineutral-electron substitution needs?

The plan is to replace the electron's ODE row

    delta[e] = res[e] - dydt[e]

with a pure algebraic charge-conservation row

    delta[e] = sum_j z_j * y_j          (no dydt term at all)

Whether that converges is a property of the pydas wrapper and the compiled
DASPK/DASSL underneath it, not something to be assumed. Three things are
checked here, each against the exact same DASx class rmgpy.solver.base imports:

  D1  a plain ODE pair integrates (positive control -- the harness is sound)
  D2  the same system with row 1 made algebraic integrates, and agrees with D1
  D3  what happens when the algebraic row's initial state is INCONSISTENT
      (y1 != y0 at t0), since that is the failure mode a non-neutral deck
      would produce

Also recorded: RMG's set_initial_derivative computes dydt0 = -residual(y0, 0),
and every RMG reactor writes delta = res - dydt, so dydt0 comes out as MINUS
the true derivative. That is pre-existing and project-wide; a DAE is the place
it could start to matter, so D4 measures whether it does.
"""

import numpy as np

import rmgpy
from rmgpy.solver.base import ReactionSystem   # noqa: F401  (forces the same build)

# the exact solver class rmgpy/solver/base.pyx selects, via the same settings.pxi
try:
    from pydas.daspk import DASPK as DASx
    from pydas.daspk import DASPKError as DASxError
    BACKEND = 'DASPK'
except ImportError:                                            # pragma: no cover
    from pydas.dassl import DASSL as DASx
    from pydas.dassl import DASSLError as DASxError
    BACKEND = 'DASSL'

print("rmgpy from : {0}".format(rmgpy.__file__))
print("backend    : {0}".format(BACKEND))
print()

K = 3.0
TEND = 1.0


class Ode(DASx):
    """dy0/dt = -K*y0 ;  dy1/dt = -K*y1.  Both rows differential."""
    algebraic = False

    def residual(self, t, y, dydt, senpar=np.zeros(1, float)):
        res = np.array([-K * y[0], -K * y[1]])
        delta = res - dydt
        if self.algebraic:
            # row 1 becomes the constraint y1 - y0 = 0, with NO dydt term
            delta[1] = y[1] - y[0]
        return delta, 1

    def jacobian(self, t, y, dydt, cj, senpar=np.zeros(1, float)):
        pd = -cj * np.identity(2, float)
        pd[0, 0] += -K
        if self.algebraic:
            pd[1, 0] = -1.0
            pd[1, 1] = +1.0          # note: no -cj on an algebraic row
        else:
            pd[1, 1] += -K
        return pd


def run(algebraic, y0, dydt0_sign=+1, label=""):
    s = Ode()
    s.algebraic = algebraic
    y0 = np.array(y0, float)
    # mirror rmgpy.solver.base.ReactionSystem.set_initial_derivative exactly:
    #     dydt0 = -residual(t0, y0, zeros)[0]
    dydt0 = -s.residual(0.0, y0, np.zeros(2, float))[0] * dydt0_sign
    atol = np.ones(2, float) * 1e-16
    rtol = np.ones(2, float) * 1e-8
    try:
        s.initialize(0.0, y0, dydt0, np.zeros(1, float), atol, rtol)
        s.advance(TEND)
        exact = y0[0] * np.exp(-K * TEND)
        err0 = abs(s.y[0] - exact) / exact
        err1 = abs(s.y[1] - s.y[0]) / abs(s.y[0])
        print("  {0:<46} OK   y = [{1:.12e}, {2:.12e}]".format(label, s.y[0], s.y[1]))
        print("  {0:<46}      dydt0 used = {1}".format("", dydt0))
        print("  {0:<46}      rel err vs exp(-Kt) = {1:.3e}; |y1-y0|/y0 = {2:.3e}".format(
            "", err0, err1))
        return True, s.y.copy()
    except DASxError as e:
        print("  {0:<46} FAIL {1}".format(label, e))
        return False, None
    except Exception as e:                                     # pragma: no cover
        print("  {0:<46} FAIL {1}: {2}".format(label, type(e).__name__, e))
        return False, None


print("D1  positive control: two differential rows")
ok1, y1 = run(False, [1.0, 1.0], label="ODE, consistent y0=[1,1]")
print()

print("D2  the substitution: row 1 algebraic (y1 - y0 = 0)")
ok2, y2 = run(True, [1.0, 1.0], label="DAE, consistent y0=[1,1]")
print()

print("D3  inconsistent initial state on the algebraic row (y1 != y0)")
ok3, y3 = run(True, [1.0, 1.5], label="DAE, INCONSISTENT y0=[1,1.5]")
print()

print("D4  the dydt0 sign: RMG's -residual(y0,0) against the true derivative")
ok4, y4 = run(True, [1.0, 1.0], dydt0_sign=-1, label="DAE, dydt0 sign flipped")
print()

print("=" * 78)
print("D1 differential control integrates : {0}".format(ok1))
print("D2 ALGEBRAIC ROW INTEGRATES        : {0}".format(ok2))
if ok1 and ok2:
    print("   D1 vs D2 agreement on y0        : {0:.3e} relative".format(
        abs(y2[0] - y1[0]) / y1[0]))
print("D3 inconsistent IC                 : {0} {1}".format(
    ok3, "(accepted -- so a non-neutral deck must be refused by US, not by the solver)"
    if ok3 else "(refused by the solver)"))
print("D4 dydt0 sign sensitivity          : {0} {1}".format(
    ok4, "(insensitive)" if ok4 else "(SENSITIVE -- the pre-existing sign matters here)"))
print("=" * 78)
