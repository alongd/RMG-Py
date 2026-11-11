#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'), #
# to deal in the Software without restriction, including without limitation  #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:       #
#                                                                             #
# The above copyright notice and this permission notice shall be included in #
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

"""
CoolProp must load before any openbabel format plugin.

CoolProp's extension and openbabel 3.2's JSON format plugins each embed rapidjson, whose static
tables are STB_GNU_UNIQUE symbols bound process-wide to the first copy loaded. If an openbabel
plugin loads first, CoolProp's static initialiser spins forever inside ``import CoolProp``.

Each test runs a child process: this pytest process has already loaded openbabel through
``test/conftest.py``, so an in-process ``import CoolProp`` could hang the whole run. The timeout
is the assertion.
"""

import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
PROBE = REPO_ROOT / "test" / "rmgpy" / "test_data" / "coolprop_load_order" / "solvation_import_probe.py"
TIMEOUT_S = 60

needs_coolprop = pytest.mark.skipif(
    importlib.util.find_spec("CoolProp") is None,  # finds the package without loading the extension
    reason="CoolProp is not installed",
)


def _text(data):
    if data is None:
        return ""
    return data.decode(errors="replace") if isinstance(data, bytes) else data


def _tail(text, n=20):
    return "\n".join(text.splitlines()[-n:])


def _run_child(args):
    """Run ``python <args>`` from the repo root with this tree on PYTHONPATH.

    Returns ``(returncode or None on timeout, stdout, stderr)``.
    """
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join([str(REPO_ROOT)] + [p for p in [env.get("PYTHONPATH")] if p])
    env["PYTHONUNBUFFERED"] = "1"
    try:
        proc = subprocess.run(
            [sys.executable] + args,
            cwd=REPO_ROOT,
            env=env,
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            timeout=TIMEOUT_S,
        )
    except subprocess.TimeoutExpired as exc:
        return None, _text(exc.stdout), _text(exc.stderr)
    return proc.returncode, proc.stdout, proc.stderr


def _report(stdout, stderr):
    return f"\n--- last 20 lines of stdout ---\n{_tail(stdout)}\n--- last 20 lines of stderr ---\n{_tail(stderr)}"


@needs_coolprop
def test_coolprop_import_after_openbabel_conversion():
    """A program that imports rmgpy and converts a molecule can then import CoolProp."""
    code = (
        "import rmgpy\n"
        "from rmgpy.molecule import Molecule\n"
        "Molecule().from_smiles('CCO').to_inchi()\n"
        "import CoolProp\n"
        "print('OK')\n"
    )
    rc, out, err = _run_child(["-c", code])
    if rc is None:
        pytest.fail(
            f"child did not finish in {TIMEOUT_S} s: import CoolProp hung (an openbabel format plugin "
            f"loaded before CoolProp's extension){_report(out, err)}"
        )
    assert rc == 0, f"child exited with {rc}{_report(out, err)}"
    assert "OK" in out.splitlines(), f"child did not reach the end{_report(out, err)}"


@needs_coolprop
def test_plain_pytest_through_conftest_imports_solvation():
    """A plain ``python -m pytest`` in the tree (so test/conftest.py applies) can import solvation."""
    rc, out, err = _run_child(["-m", "pytest", "-p", "no:cacheprovider", "--no-cov", "-q", str(PROBE)])
    lines = out.splitlines()
    last = lines[-1].strip() if lines else ""
    if rc is None:
        pytest.fail(
            f"child pytest did not finish in {TIMEOUT_S} s: it hung (last stdout line {last!r}; "
            f"test/conftest.py loaded openbabel before CoolProp){_report(out, err)}"
        )
    assert last != "collecting ...", f"child pytest stopped while collecting{_report(out, err)}"
    assert rc == 0, f"child pytest exited with {rc}{_report(out, err)}"
    assert "1 passed" in out, f"probe test did not pass{_report(out, err)}"


def test_rmgpy_imports_without_coolprop():
    """With CoolProp absent, ``import rmgpy`` behaves as it does without the load-order import."""
    code = (
        "import sys\n"
        "class _NoCoolProp:\n"
        "    def find_spec(self, name, path=None, target=None):\n"
        "        if name == 'CoolProp' or name.startswith('CoolProp.'):\n"
        "            raise ImportError('CoolProp blocked for this test')\n"
        "        return None\n"
        "sys.meta_path.insert(0, _NoCoolProp())\n"
        "import rmgpy\n"
        "assert rmgpy.settings is not None\n"
        "assert 'CoolProp' not in sys.modules, 'CoolProp was loaded'\n"
        "print('OK')\n"
    )
    rc, out, err = _run_child(["-c", code])
    assert rc == 0, f"import rmgpy failed without CoolProp (rc={rc}){_report(out, err)}"
    assert "OK" in out.splitlines(), f"child did not reach the end{_report(out, err)}"
