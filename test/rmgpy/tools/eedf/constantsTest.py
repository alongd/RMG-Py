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

"""Keep frozen EEDF physical identities independent of the installed SciPy."""

import copy
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
import scipy.constants

from rmgpy.tools.eedf.constants import (
    AVOGADRO, BOLTZMANN, ELECTRON_MASS, ELEMENTARY_CHARGE,
)
from rmgpy.tools.eedf.generate import model_inputs
from rmgpy.tools.eedf.schema import SpecError, _file_hash, load_spec


def test_constants_are_codata_2018_and_exact_si():
    assert ELECTRON_MASS == 9.1093837015e-31
    assert BOLTZMANN == 1.380649e-23
    assert ELEMENTARY_CHARGE == 1.602176634e-19
    assert AVOGADRO == 6.02214076e23


def test_frozen_model_inputs_ignore_scipy_mass_update(tmp_path):
    path = os.environ.get('EEDF_REAL_SPEC')
    if path is None:
        pytest.skip('EEDF_REAL_SPEC supplies the complete pinned fixture')
    # Reloading in the pytest process replaces exception classes retained by
    # other collected tests. Keep all strict checks in a separate interpreter.
    result = subprocess.run(
        [sys.executable, '-c',
         "import runpy, sys; from pathlib import Path; "
         "runpy.run_path(sys.argv[1])['_check_frozen_model_inputs']"
         "(sys.argv[2], Path(sys.argv[3]))",
         str(Path(__file__).resolve()), path, str(tmp_path)],
        capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def _check_frozen_model_inputs(path, tmp_path):
    spec = load_spec(path)
    original_map = Path(spec['channel_map']['path']).read_bytes()
    mapping = json.loads(original_map)
    elastic = next(channel for channel in mapping if channel['kind'] == 'elastic')
    original = model_inputs(spec)
    modules = [importlib.import_module('rmgpy.tools.eedf.' + name) for name in (
        'moments', 'channels', 'loki', 'channel_map', 'vibrational',
        'maxwellian', 'corroborate',
    )]
    try:
        with pytest.MonkeyPatch.context() as patch:
            patch.setattr(scipy.constants, 'electron_mass', 9.1093837139e-31)
            # Reload after the patch: already cached SciPy values must not mask it.
            for module in modules:
                importlib.reload(module)
                if hasattr(module, 'electron_mass'):
                    assert module.electron_mass == ELECTRON_MASS
            assert model_inputs(spec) == original
            assert Path(spec['channel_map']['path']).read_bytes() == original_map

            # The gate remains exact: even a one-ULP change in a separately
            # hashed test map must fail, not pass under a numeric tolerance.
            elastic['mass_ratio'] = float(np.nextafter(elastic['mass_ratio'], np.inf))
            changed_map = tmp_path / 'changed-map.json'
            changed_map.write_text(json.dumps(mapping))
            changed = copy.deepcopy(spec)
            changed['channel_map'] = {'path': str(changed_map),
                                      'sha256': _file_hash(changed_map)}
            with pytest.raises(SpecError, match='mass_ratio differs from hashed LoKI mass'):
                model_inputs(changed)
    finally:
        for module in modules:
            importlib.reload(module)
