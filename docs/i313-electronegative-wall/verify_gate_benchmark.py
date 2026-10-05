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

"""Retained 120-pair integration regression against a compiled base solver."""
import argparse
import importlib.util
import json
import logging
import os
from pathlib import Path
import statistics
import sys
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline-so', required=True)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    os.sched_setaffinity(0, {max(os.sched_getaffinity(0)) - 1})
    root = Path.cwd()
    sys.path[:0] = [str(root), str(root / 'test/rmgpy/solver')]
    import rmgpy
    import rmgpy.solver.plasma as tip
    assert root in Path(rmgpy.__file__).resolve().parents
    assert root in Path(tip.__file__).resolve().parents
    import plasmaEnergyBalanceTest as fixture
    import rmgpy.data.rmg as data
    data.database = None
    spec = importlib.util.spec_from_file_location('rmgpy.solver.plasma', args.baseline_so)
    base = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base)
    assert tip.PlasmaReactor is not base.PlasmaReactor
    logging.getLogger().setLevel(logging.ERROR)
    logger = logging.getLogger('i313_gate_benchmark')
    logger.setLevel(logging.INFO)
    logger.addHandler(logging.StreamHandler(sys.stdout))
    logger.propagate = False
    logger.info('IMPORTS %s %s', tip.__file__, base.__file__)

    def cls_for(parent):
        class Reactor(parent):
            def _declared_entry_kinetics(self, library, index):
                return fixture.ToyLibraryReactor.toy_entries[index]

            def __init__(self, *args, **kwargs):
                scalar = kwargs.pop('ion_reduced_mobility')
                kwargs['ion_reduced_mobilities'] = {'Ar+': scalar}
                super().__init__(*args, **kwargs)
        return Reactor

    classes = {'base': cls_for(base.PlasmaReactor), 'tip': cls_for(tip.PlasmaReactor)}
    results = {}
    for mode in ('batch', 'five_endpoints'):
        times = {'base': [], 'tip': []}
        ratios = []
        for i in range(124):
            states, pair = {}, {}
            order = ('base', 'tip') if i % 2 == 0 else ('tip', 'base')
            for name in order:
                reactor, _, _ = fixture._build(cls=classes[name], elastic=False)
                start = time.perf_counter_ns()
                if mode == 'batch':
                    reactor.advance(1e-6)
                else:
                    for endpoint in (1e-10, 1e-9, 1e-8, 1e-7, 1e-6):
                        reactor.advance(endpoint)
                pair[name] = (time.perf_counter_ns() - start) / 1e9
                states[name] = (reactor.y.tobytes(), reactor.wall_flux.tobytes(),
                                reactor.nu_wall_latched, reactor.energy_budget)
            assert states['base'] == states['tip']
            if i >= 4:
                for name in pair:
                    times[name].append(pair[name])
                ratios.append(pair['tip'] / pair['base'])
        results[mode] = {'times_s': times, 'ratios': ratios}
        logger.info('PAIRED %s base_ms %s tip_ms %s ratio_median %s ratio_quartiles %s',
                    mode, statistics.median(times['base']) * 1000,
                    statistics.median(times['tip']) * 1000, statistics.median(ratios),
                    statistics.quantiles(ratios, n=4))
    args.output.write_text(json.dumps(results, indent=2) + '\n')
    for mode, result in results.items():
        assert statistics.median(result['ratios']) <= 1.05, (
            mode, statistics.median(result['ratios']))


if __name__ == '__main__':
    main()
