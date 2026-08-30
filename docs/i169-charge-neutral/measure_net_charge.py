#!/usr/bin/env python
"""
Read a plasma deck through rmgpy.rmg.input.read_input_file and report, for every
PlasmaReactor it defines, the initial composition with per-species charge, the net
charge per mole, and the total charge magnitude the net is measured against.

Usage:  python measure_net_charge.py <input.py>
Run from a directory holding an rmgrc; PYTHONPATH must point at this worktree.
"""
import sys

from rmgpy.rmg.main import RMG
from rmgpy.rmg.input import read_input_file
from rmgpy.solver.plasma import PlasmaReactor


def charge_of(spc):
    try:
        return spc.get_net_charge()
    except Exception:
        return None


def main(path):
    rmg = RMG(input_file=path, output_directory='.')
    read_input_file(path, rmg)
    for system in rmg.reaction_systems:
        if not isinstance(system, PlasmaReactor):
            continue
        print('reactor: PlasmaReactor T {0} P {1} Te {2}'.format(
            system.T.value_si, system.P.value_si, system.Te.value_si))
        net = 0.0
        scale = 0.0
        total_x = 0.0
        for spc, x in system.initial_mole_fractions.items():
            z = charge_of(spc)
            print('  {0:<8s} x={1!r:<26s} charge={2}'.format(str(spc), repr(x), z))
            total_x += x
            if z:
                net += x * z
                scale += abs(x * z)
        print('  sum x             = {0!r}'.format(total_x))
        print('  net charge/mol    = {0!r}'.format(net))
        print('  charge magnitude  = {0!r}'.format(scale))
        if scale > 0.0:
            print('  |net| / magnitude = {0!r}'.format(abs(net) / scale))


if __name__ == '__main__':
    main(sys.argv[1])
