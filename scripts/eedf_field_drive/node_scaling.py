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

"""Measured node budget gate; no generation runs until this gate passes."""
import copy
import json
import os
from pathlib import Path
import re
import time
import uuid

import numpy as np

from job_loki import prepare_spec, install_refusal_capture, install_row_timing
from rmgpy.tools.eedf.channels import _load_channel_map
from rmgpy.tools.eedf.generate import _SolvePool, _worker_count
from rmgpy.tools.eedf.loki import LoKIDriver
from rmgpy.tools.eedf.schema import _validate_spec

# Recorded complete local rows and actual round-0 call counts. No science is cut.
LOCAL_SERIAL_SIX_SECONDS = 41.34543256700272
LOCAL_SCREEN_SECONDS = 4 * 106.652 + 1510.3972046090057
LOCAL_ROUND0_SECONDS = (
    1108 * (106.652 / 222)
    + 6012 * (1533.093495823996 / 222)
    + 2060 * (1510.3972046090057 / 222))
BUILD_RESERVE_SECONDS = 1309
MAX_PROJECTED_SECONDS = 16.8 * 3600
POINTS = [(f, field) for f in (1e-5, 1e-4)
          for field in (0.002447633990682018, 0.3, 30.)]


def projection(serial_seconds, parallel_seconds, elapsed_seconds):
    """Replace the aggregate node proxy with the measured six-point node ratio."""
    values = (serial_seconds, parallel_seconds, elapsed_seconds)
    if not all(np.isfinite(x) for x in values) or min(values[:2]) <= 0 or elapsed_seconds < 0:
        raise ValueError('invalid node timing')
    factor = serial_seconds / LOCAL_SERIAL_SIX_SECONDS
    speedup = serial_seconds / parallel_seconds
    seconds = (elapsed_seconds + BUILD_RESERVE_SECONDS
               + LOCAL_SCREEN_SECONDS * factor
               + LOCAL_ROUND0_SECONDS * factor / speedup)
    return {'node_local_factor': factor, 'speedup': speedup,
            'projected_total_seconds': seconds,
            'projected_total_hours': seconds / 3600,
            'limit_hours': 16.8, 'continue': seconds <= MAX_PROJECTED_SECONDS,
            'assumptions': 'six-point scaling represents full scans/cold/held/T10; '
            'unmeasured fraction bands use measured 1e-5/1e-4 full-row averages; '
            'round 0 only, refinement extent unknown; remaining build reserve 1309s'}


def numeric_outputs(work):
    """Read every reported EEDF/anisotropy, swarm, power and rate numeric value."""
    folders = list((Path(work) / 'output/solve').glob('reducedField_*'))
    if len(folders) != 1:
        raise ValueError('probe must produce exactly one field')
    folder = folders[0]
    rates = []
    for line in (folder / 'rateCoefficients.txt').read_text().splitlines()[1:]:
        text = line.strip()
        if not text or set(text) == {'-'} or text == '* Extra Rate Coefficients *':
            continue
        fields = text.split()
        rates.append([float(fields[0]), float(fields[1])])
    if not rates:
        raise ValueError('empty rate output')
    result = {'eedf': np.loadtxt(folder / 'eedf.txt', skiprows=1),
              'rates': np.asarray(rates)}
    for name, filename in [('swarm', 'swarmParameters.txt'), ('power', 'powerBalance.txt')]:
        values = {}
        for line in (folder / filename).read_text().splitlines():
            match = re.match(r'\s*(.*?)\s*=\s*([+\-0-9.eE]+)', line)
            if match:
                label, value = match.groups()
                # This percentage is derived from the near-zero dimensional residual.
                if label == 'Relative Power Balance':
                    continue
                values[label] = float(value)
        if not values:
            raise ValueError('empty numeric probe output: ' + name)
        result[name] = values
    if 'Power Balance' not in result['power']:
        raise ValueError('missing dimensional cancellation residual')
    return result


def compare_values(reference, candidate):
    """Apply the PM's execution-only bounds, never table qualification tolerances."""
    if set(reference) != set(candidate):
        raise ValueError('different numeric output keys')
    maxima = {}
    for group in reference:
        a, b = reference[group], candidate[group]
        if isinstance(a, dict):
            if set(a) != set(b):
                raise ValueError('different ' + group + ' values')
            pairs = [(key, a[key], b[key]) for key in a]
        else:
            pairs = [('', a, b)]
        for key, x, y in pairs:
            x, y = np.asarray(x), np.asarray(y)
            if x.shape != y.shape or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
                raise ValueError('invalid numeric output: ' + group + '.' + key)
            if group == 'power' and key == 'Power Balance':
                if max(float(abs(x)), float(abs(y))) > 1e-20:
                    raise ValueError('dimensional cancellation residual exceeds 1e-20')
            elif not np.allclose(y, x, rtol=1e-6, atol=1e-20):
                raise ValueError('numeric equivalence failed: ' + group + '.' + key)
            maxima[group + '.' + key] = float(np.max(np.abs(y - x)))
    return maxima


def run_probe(root, spec):
    results = root / 'results'
    record = {'status': 'running', 'continue': False,
              'criteria': {'rtol': 1e-6, 'atol': 1e-20,
                           'power_balance': 'absolute residual <=1e-20 in both runs'},
              'points': POINTS, 'measurements': {}}
    target = results / 'node-scaling.json'

    def save():
        target.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')

    save()
    try:
        if _worker_count(spec) != 3:
            raise ValueError('node allocation cannot support three five-thread solvers')
        prefix = 'scaling_' + uuid.uuid4().hex
        mapping = _load_channel_map(spec['channel_map'])
        outputs = {}
        for mode, count in [('serial', 1), ('parallel', 3)]:
            candidate = copy.deepcopy(spec)
            candidate['max_parallel'] = count
            labels = [prefix + '_' + mode + '_' + str(i) for i in range(len(POINTS))]
            start = time.monotonic()
            with _SolvePool(LoKIDriver(candidate), mapping, candidate) as pool:
                pending = [pool.submit({'electron_density': 5.2e15, 'f_Ar4s': f},
                                       [float(np.log(field))], label)
                           for (f, field), label in zip(POINTS, labels)]
                for future in pending:
                    if len(future.result()) != 1:
                        raise ValueError('probe row count')
            record['measurements'][mode] = {'wall_seconds': time.monotonic() - start}
            outputs[mode] = [numeric_outputs(Path(spec['scratch_root']) / label)
                             for label in labels]
            save()
        record['numeric_equivalence'] = [
            compare_values(a, b) for a, b in zip(outputs['serial'], outputs['parallel'])]
        elapsed = time.time() - float((results / 'job-started-epoch').read_text())
        record.update(projection(record['measurements']['serial']['wall_seconds'],
                                 record['measurements']['parallel']['wall_seconds'], elapsed))
        record['status'] = 'passed' if record['continue'] else 'budget_refused'
    except Exception as error:
        record.update(status='refused',
                      exception=type(error).__name__, message=str(error))
        record['continue'] = False
    finally:
        save()
    return record


def main():
    root = Path(__file__).resolve().parent
    spec = prepare_spec(root)
    _validate_spec(spec)
    install_refusal_capture(root)
    install_row_timing(root)
    record = run_probe(root, spec)
    print(json.dumps(record, indent=2), flush=True)


if __name__ == '__main__':
    main()
