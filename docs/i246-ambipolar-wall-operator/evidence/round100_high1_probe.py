import numpy as np
import sys
sys.path.insert(0, '.')
from rmgpy.solver.termination import TerminationSteadyState

ATOL = 1e-16

def feed(term, r_gen_seq, ext_armed_seq):
    """Drive update() with a synthetic generic residual sequence and an external-armed
    flag per step. One constant species (1e12) dominates; the second at exp(r) sets the
    generic residual over an e-fold. external_residual is left nan (electron already flat)."""
    fired = None
    for k, (rg, ea) in enumerate(zip(r_gen_seq, ext_armed_seq)):
        j = getattr(term, '_pf', 0)
        big = 1e12
        y_prev = np.array([big, float(np.exp(0.0))])
        y_now  = np.array([big, float(np.exp(rg))])   # dln x2 = rg over one e-fold -> R=rg
        t_prev = float(np.exp(j)); t_now = float(np.exp(j+1))
        term._pf = j + 1
        got = term.update(y_now, t_now, y_prev, t_prev, ATOL,
                          external_residual=float('nan'), external_armed=ea)
        if got and fired is None:
            fired = k
    return fired

# A generic (neutral) channel RISING THROUGHOUT but sub-tolerance, with ONE flat sample
# at step 3. The electron has passed its relaxation time the whole time (external_armed=True).
# The neutral is still departing -> the criterion MUST NOT arm/terminate.
rising = [1e-9, 2e-9, 3e-9, 3e-9, 4e-9, 5e-9, 6e-9, 7e-9, 8e-9, 9e-9]
ext    = [True]*len(rising)

term = TerminationSteadyState(tolerance=1e-6, window=3)
fired = feed(term, rising, ext)
print("CURRENT built module:")
print("  armed_external =", term.armed_external, " armed =", term.armed, " fired_at =", fired)
print("  -> HIGH 1 present if armed_external is True (a still-rising neutral armed the criterion)")
