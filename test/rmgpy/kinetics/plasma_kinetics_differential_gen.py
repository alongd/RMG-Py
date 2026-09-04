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

"""
Generator for the ``origin/99`` reference fixture consumed by
``plasmaKineticsConsolidationTest.py``.

The fixture pins the numeric output of the four plasma kinetics classes
(``TwoTemperaturePlasma``, ``ElectronCollisionPlasma``, ``BadnellRRArrhenius``,
``VoronovEIArrhenius``) as computed by the *older* ``origin/99`` branch, so that
the consolidated ``plasma`` implementation can be shown to reproduce it. It is a
differential test: the fixture is the reference side, the live import is the
side under test.

To regenerate the reference side against ``origin/99`` (not needed for the test
itself, only to refresh the committed fixture):

    git worktree add --detach ../RMG-Py-99-ref origin/99
    cp rmgrc ../RMG-Py-99-ref/rmgrc        # database.directory -> plasma DB
    ( cd ../RMG-Py-99-ref && python utilities.py check-pydas \\
        && PYTHONPATH=$PWD python setup.py build_ext --inplace -j 8 )
    PYTHONPATH=../RMG-Py-99-ref python test/rmgpy/kinetics/plasma_kinetics_differential_gen.py out.json

The same script, run against the current tree, produces the "under test" side;
the test does that in-process via :func:`evaluate` rather than shelling out.

The construction parameters below are the single source of truth shared by the
fixture and the live comparison; the test module imports them from here.
"""
import json
import math
import sys

import numpy as np

import rmgpy.constants as constants
from rmgpy.kinetics.arrhenius import (
    TwoTemperaturePlasma,
    ElectronCollisionPlasma,
    BadnellRRArrhenius,
    VoronovEIArrhenius,
)

Na = constants.Na

# --- temperature spreads (log-spaced, several decades) -----------------------
SPREAD_T = [round(x, 6) for x in np.logspace(math.log10(300), math.log10(30000), 15)]
SPREAD_TE = [round(x, 6) for x in np.logspace(math.log10(1000), math.log10(50000), 15)]
GRID_T = [300.0, 1000.0, 3000.0]
GRID_TE = [1000.0, 5000.0, 20000.0]

# --- shared construction parameters ------------------------------------------
TWOTEMP_KW = dict(A=(2.0e-13, "cm^3/(molecule*s)"), n=0.5,
                  Ea_g=(50000.0, "J/mol"), Ea_e=(30000.0, "J/mol"), T0=(300.0, "K"))
ECP_ENERGIES = ([0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0], "eV/molecule")
ECP_SIGMA = ([0.0, 1.0e-21, 5.0e-21, 8.0e-21, 9.0e-21, 9.5e-21, 1.0e-20], "m^2")
BADNELL_KW = dict(A=(1.0e-12, "cm^3/(molecule*s)"), B=0.7, T0=(100.0, "K"),
                  T1=(1.0e7, "K"), C=0.3, T2=(1.0e5, "K"))
VORONOV_KW = dict(A=(2.91e-8, "cm^3/(molecule*s)"), P=0.0, X=0.232, K=0.39, dE=(13.6, "eV"))
BADNELL_ZN = [(1, 0), (2, 1), (6, 3)]
VORONOV_ZN = [(1, 1), (2, 2), (6, 6)]
A_SETTER_SPELLINGS = ["cm^3/(molecule*s)", "cm^3/(mol*s)", "m^3/(molecule*s)",
                      "m^3/(mol*s)", "cm^3/s", "m^3/s"]


def safe(fn, *a):
    try:
        return fn(*a)
    except Exception as e:  # noqa: BLE001 - surface the failure verbatim in the dump
        return "ERR:" + type(e).__name__ + ":" + str(e)[:80]


def call_opt(obj, method, *a):
    fn = getattr(obj, method, None)
    if fn is None:
        return "NO_METHOD"
    return safe(fn, *a)


def evaluate(include_yaml=True):
    """Return the full case dictionary for whichever build is imported.

    ``include_yaml`` gates the (Z,N) table-driven cases, which require the plasma
    RMG-database to be configured.
    """
    out = {"meta": {"spread_T": SPREAD_T, "spread_Te": SPREAD_TE,
                    "grid_T": GRID_T, "grid_Te": GRID_TE, "Na": Na}}

    ttp = TwoTemperaturePlasma(**TWOTEMP_KW)
    out["TwoTemp_grc"] = {str(T): safe(ttp.get_rate_coefficient, T) for T in SPREAD_T}
    out["TwoTemp_two_temp"] = {f"{T},{Te}": safe(ttp.get_rate_coefficient_two_temp, T, Te)
                               for T in GRID_T for Te in GRID_TE}
    out["TwoTemp_A_si"] = ttp._A.value_si

    ecp = ElectronCollisionPlasma(energies=ECP_ENERGIES, sigma=ECP_SIGMA)
    out["ECP_grc"] = {str(Te): safe(ecp.get_rate_coefficient, Te) for Te in SPREAD_TE}
    out["ECP_electron_temp"] = {str(Te): call_opt(ecp, "get_rate_coefficient_electron_temp", Te)
                                for Te in SPREAD_TE}

    brr = BadnellRRArrhenius(**BADNELL_KW)
    out["Badnell_grc"] = {str(T): safe(brr.get_rate_coefficient, T) for T in SPREAD_T}
    out["Badnell_electron_temp"] = {str(Te): call_opt(brr, "get_rate_coefficient_electron_temp", Te)
                                    for Te in SPREAD_T}
    out["Badnell_A_si"] = brr._A.value_si

    vei = VoronovEIArrhenius(**VORONOV_KW)
    out["Voronov_grc"] = {str(Te): safe(vei.get_rate_coefficient, Te) for Te in SPREAD_TE}
    out["Voronov_electron_temp"] = {str(Te): call_opt(vei, "get_rate_coefficient_electron_temp", Te)
                                    for Te in SPREAD_TE}

    setter = {}
    for spelling in A_SETTER_SPELLINGS:
        try:
            v = VoronovEIArrhenius(A=(1.0e-8, spelling), P=0.0, X=0.3, K=0.25, dE=(13.6, "eV"))
            setter[spelling] = {"A_si": v._A.value_si, "units": v._A.units,
                                "k_10000K": safe(v.get_rate_coefficient, 10000.0)}
        except Exception as e:  # noqa: BLE001
            setter[spelling] = "ERR:" + type(e).__name__ + ":" + str(e)[:120]
    out["Voronov_A_setter"] = setter

    if include_yaml:
        for (Z, N) in BADNELL_ZN:
            try:
                b = BadnellRRArrhenius(Z=Z, N=N)
                out[f"Badnell_yaml_{Z}_{N}"] = {str(T): safe(b.get_rate_coefficient, T) for T in SPREAD_T}
            except Exception as e:  # noqa: BLE001
                out[f"Badnell_yaml_{Z}_{N}"] = "CONSTRUCT_ERR:" + type(e).__name__ + ":" + str(e)[:120]
        for (Z, N) in VORONOV_ZN:
            try:
                v = VoronovEIArrhenius(Z=Z, N=N)
                out[f"Voronov_yaml_{Z}_{N}"] = {str(Te): safe(v.get_rate_coefficient, Te) for Te in SPREAD_TE}
            except Exception as e:  # noqa: BLE001
                out[f"Voronov_yaml_{Z}_{N}"] = "CONSTRUCT_ERR:" + type(e).__name__ + ":" + str(e)[:120]

    return out


if __name__ == "__main__":
    dst = sys.argv[1] if len(sys.argv) > 1 else "plasma_kinetics_dump.json"
    with open(dst, "w") as fh:
        json.dump(evaluate(), fh, indent=1, sort_keys=True)
    import rmgpy.kinetics.arrhenius as _arr
    print("WROTE", dst, "from", _arr.__file__)
