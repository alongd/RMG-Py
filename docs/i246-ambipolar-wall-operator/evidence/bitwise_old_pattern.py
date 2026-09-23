#!/usr/bin/env python3
"""The OLD bitwise test's logic: build two WALL-LESS reactors and compare their
residual and Jacobian. This never evaluates the wall path, so a corruption of
_apply_wall_terms leaves it green -- which is exactly why it could not fail. Run
this under a corrupted build alongside the rewritten pytest to show the contrast.
"""
import numpy as np

import sys
sys.path.insert(0, 'test/rmgpy/solver')
from plasmaWallTest import _build_reactor, _state_at  # noqa: E402


def main():
    ra, _, _ = _build_reactor(wall=False, with_chemistry=True)
    rb, _, _ = _build_reactor(wall=False, with_chemistry=True)
    y = _state_at(ra, 1.0e-6)
    dydt = np.zeros(ra.num_core_species, float)
    da, ia = ra.residual(0.0, y.copy(), dydt.copy())
    db, ib = rb.residual(0.0, y.copy(), dydt.copy())
    pa = np.array(ra.jacobian(0.0, y.copy(), dydt.copy(), 1.0), float)
    pb = np.array(rb.jacobian(0.0, y.copy(), dydt.copy(), 1.0), float)
    ok = np.array_equal(da, db) and ia == ib and np.array_equal(pa, pb)
    print('OLD two-wall-less pattern:', 'GREEN (vacuous)' if ok else 'RED')


if __name__ == '__main__':
    main()
