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

import os
import hashlib
import json
import logging
from pathlib import Path
import re
import subprocess
import sys
import time
import tarfile
import numpy as np
from job_loki import prepare_spec, install_row_timing, install_refusal_capture
import scipy
import scipy.constants
import rmgpy
from rmgpy.solver import plasma
from rmgpy.solver.eedf import EEDFTable
from rmgpy.tools.eedf.generate import generate, model_inputs
from rmgpy.tools.eedf.integrity import _solver_environment
from rmgpy.tools.eedf.schema import _file_hash, _validate_spec
if sys.argv[1:] not in ([], ['--preflight-only']):
    raise SystemExit('usage: run_table.py [--preflight-only]')
root = Path(__file__).resolve().parent
results = root / 'results'
assert Path(rmgpy.__file__).resolve().is_relative_to(root / 'rmg')
assert Path(plasma.__file__).resolve().is_relative_to(root / 'native')
spec = prepare_spec(root)
install_refusal_capture(root)
install_row_timing(root)
# Persist numerical identities even when the physical gate refuses this runtime.
from rmgpy.tools.eedf.channels import _gas_values, _load_channel_map, electron_mass
from rmgpy.tools.eedf.integrity import _physical_properties
gas_properties, _ = _physical_properties(spec)
mass_pin = spec['input_files'][gas_properties['mass']]
masses = _gas_values(mass_pin['path'], mass_pin['sha256'])
identity = {'python': sys.version, 'numpy': np.__version__, 'scipy': scipy.__version__,
            'scipy_path': scipy.__file__, 'electron_mass': electron_mass,
            'electron_mass_provenance': 'CODATA 2018, pinned by rmgpy.tools.eedf.constants',
            'scipy_electron_mass': scipy.constants.electron_mass,
            'mass_pin': mass_pin, 'channel_map_pin': spec['channel_map'], 'elastic': []}
for channel in _load_channel_map(spec['channel_map']):
    if channel['kind'] == 'elastic':
        species = channel['description'].split(' ', 3)[2].split('(', 1)[0]
        mass = masses[species]
        identity['elastic'].append({'description': channel['description'], 'mass': mass,
                                   'declared_ratio': channel['mass_ratio'],
                                   'computed_ratio': electron_mass / mass if mass > 0 else None})
(results / 'physical-identity.json').write_text(json.dumps(identity, indent=2) + '\n')
_validate_spec(spec)
inputs = model_inputs(spec)  # Full identity/physical closure before any solve.
(results / 'model-inputs.json').write_text(json.dumps(inputs, indent=2) + '\n')
if sys.argv[1:] == ['--preflight-only']:
    print('PHYSICAL_PREFLIGHT_PASS', flush=True)
    raise SystemExit(0)
# Preserve the exact paths/bytes needed for later model_inputs and normal admission.
# Native Python extensions and raw solver scratch stay on Zeus.
with tarfile.open(results / 'provenance.tar.gz', 'w:gz', compresslevel=9) as archive:
    for name in ('inputs', 'loki-build/app/loki', 'loki-build/CMakeCache.txt',
                 'loki-build/runtime', 'loki-build/lib', 'loki-build/growth-mixing.patch',
                 'loki-build/growth-mixing-provenance.json',
                 'loki-build/growth-mixing-build.json'):
        archive.add(root / name, arcname=name)

(results / 'remote-spec.json').write_text(json.dumps(spec, indent=2) + '\n')
logging.basicConfig(level=logging.INFO)
print('LOW_FIELD_GENERATION_STARTED', len(spec['axes']['u']), spec['axes']['f_Ar4s'], flush=True)
path = generate(spec, command_line=['python', 'run_table.py'])
manifest = json.loads((path / 'manifest.json').read_text())
# Keep the published artifact byte-for-byte, including rejected-table evidence.
# Its verbose verdicts compress well; full build and solver trees never enter results/.
with tarfile.open(results / 'table.tar.gz', 'w:gz', compresslevel=9) as archive:
    archive.add(path, arcname=str(path.relative_to(root)))
(results / 'model-inputs.json').write_text(json.dumps(model_inputs(spec), indent=2) + '\n')
result = {'artifact': str(path), 'artifact_sha256': manifest['artifact_sha256'],
          'accepted': manifest['accepted'], 'T10': manifest['T10 interpolation qualification'],
          'source': manifest['A6b-source'], 'rmgpy': str(rmgpy.__file__),
          'plasma': str(plasma.__file__), 'python': sys.version}
(results / 'qualification-result.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result), flush=True)
if not manifest['accepted']:
    raise SystemExit(2)
EEDFTable.load(path, model_inputs(spec), artifact_sha256=manifest['artifact_sha256'])
result['loader_admitted'] = True
(results / 'qualification-result.json').write_text(json.dumps(result, indent=2) + '\n')
print('LOW_FIELD_TABLE_QUALIFIED', path, flush=True)
