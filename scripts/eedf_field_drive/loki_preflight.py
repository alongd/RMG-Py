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

"""One frozen screen-reference solve before compiling RMG native extensions."""
import json
from pathlib import Path
import uuid
import numpy as np

from job_loki import prepare_spec, install_refusal_capture, capture_logs
from rmgpy.tools.eedf.loki import LoKIDriver
from rmgpy.tools.eedf.schema import _validate_spec

root = Path(__file__).resolve().parent
spec = prepare_spec(root)
_validate_spec(spec)
install_refusal_capture(root)
# Index18 is the reproduced refusal; this selects a diagnostic witness,
# without editing the frozen grid or numerical settings.
index = 18
coordinate = {item['name']: item['reference'] for item in spec['screen']['candidates']}
field = float(np.exp(spec['axes']['u'][index]))
name = 'node_preflight_' + uuid.uuid4().hex
record = {'job_name': name, 'u_index': index, 'u': spec['axes']['u'][index],
          'EN_Td': field, 'coordinate': coordinate}
(root / 'results/loki-preflight.json').write_text(json.dumps(record, indent=2) + '\n')
print('LOKI_PREFLIGHT_STARTED ' + json.dumps(record), flush=True)
try:
    rows = LoKIDriver(spec).run(coordinate, [field], name)
    if len(rows) != 1:
        raise RuntimeError('single-coordinate preflight did not return one row')
except Exception as error:
    record.update(status='refused', exception=type(error).__name__, message=str(error))
    raise
else:
    record['status'] = 'passed'
    print('LOKI_PREFLIGHT_PASS', flush=True)
finally:
    (root / 'results/loki-preflight.json').write_text(json.dumps(record, indent=2) + '\n')
    capture_logs(root, Path(spec['scratch_root']) / name, 'loki-preflight', record['status'])
