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

"""Public and native lifecycle regressions in isolated Python owners."""
import os
from pathlib import Path
import subprocess
import sys

import pytest

PROBE = Path(__file__).with_name('_native_protocol_probe.py')


@pytest.mark.parametrize('mode', [
    'blocking-hook', 'gil', 'healthy', 'completion', 'race', 'expiry',
    'bootstrap-finalizer', 'external-reap', 'reap-cache', 'accept-before-reap',
    'exited-owner-stall', 'lost-popen', 'missing-handshake', 'create-error',
    'foreign-thread', 'bad-envelope',
    'unsupported-spawn', 'unsupported-forkserver', 'unsupported-extension',
    'artifact-lease', 'artifact-platform', 'artifact-startup',
])
def test_native_worker_protocol(mode):
    environment = dict(os.environ, PATH='')
    environment.pop('PYTHONPATH', None)
    result = subprocess.run([sys.executable, str(PROBE), mode], env=environment,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=15)
    assert result.returncode == 0, result.stdout + result.stderr
    assert ('NATIVE_PROTOCOL ' + mode in result.stdout
            or 'PUBLIC_ARTIFACT ' + mode + ' PASS' in result.stdout)


@pytest.mark.parametrize('mode', ['cleanup-before', 'cleanup-after', 'bookkeeping', 'unreadable',
                                  'exit-immediately', 'storm-unbounded'])
def test_native_interrupt_boundaries(mode):
    environment = dict(os.environ, PATH='')
    environment.pop('PYTHONPATH', None)
    probe = PROBE.with_name('_native_interrupt_probe.py')
    result = subprocess.run([sys.executable, str(probe), mode], env=environment,
                            capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'NATIVE_INTERRUPT' in result.stdout


@pytest.mark.parametrize('mode', ['entry', 'exception', 'return', 'early-fork', 'after-cleanup', 'report'])
def test_native_guard_boundaries(mode):
    environment = dict(os.environ, PATH='')
    environment.pop('PYTHONPATH', None)
    result = subprocess.run([sys.executable, str(PROBE.with_name('_native_guard_probe.py')), mode],
                            env=environment, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'NATIVE_GUARD' in result.stdout


def test_native_cleanup_statement_matrix():
    environment = dict(os.environ, PATH='')
    environment.pop('PYTHONPATH', None)
    result = subprocess.run([sys.executable, str(PROBE.with_name('_native_cleanup_matrix.py')), 'green'],
                            env=environment, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'CLEANUP_MATRIX' in result.stdout


@pytest.mark.parametrize('mode', ['reference-release'])
def test_native_guard_ownership_and_references(mode):
    environment = dict(os.environ, PATH='')
    environment.pop('PYTHONPATH', None)
    result = subprocess.run([sys.executable, str(PROBE.with_name('_native_guard_probe.py')), mode],
                            env=environment, capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'NATIVE_GUARD' in result.stdout


@pytest.mark.native_slow
@pytest.mark.skipif(os.environ.get('RMG_NATIVE_SLOW') != '1', reason='explicit 1000-fork qualification gate')
def test_native_fork_storm_1000(tmp_path):
    environment = dict(os.environ, PATH='')
    environment.pop('PYTHONPATH', None)
    result = subprocess.run([sys.executable, str(PROBE.with_name('_native_storm_probe.py')), str(tmp_path)],
                            env=environment, capture_output=True, text=True, timeout=240)
    print(result.stdout)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'NATIVE_STORM_RESULT' in result.stdout
    print('NATIVE_STORM_ARTIFACT', str(tmp_path))
