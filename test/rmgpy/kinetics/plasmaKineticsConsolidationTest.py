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
Differential numeric test: the four plasma kinetics classes on the consolidated
``plasma`` branch reproduce, across a spread of gas and electron temperatures,
the numbers computed by their older ``origin/99`` counterparts.

This exists so that "the engine has a single consolidated source of truth for
plasma kinetics" rests on a measurement rather than a presence check. The four
classes are declared on both branches; declaring the same class is not evidence
of computing the same number. Here the *fixture*
(``plasma_kinetics_origin99_reference.json``) is the number ``origin/99``
actually computes, captured from its compiled module at commit
``a9c3e8d148``; the *live* import is the side under test.

The rate bodies of these classes are byte-identical source compiled by the same
toolchain, so agreement is expected to machine precision. The tolerance is
``rtol = 1e-12`` — a few ulps of double precision, room only for benign
floating-point reassociation, not for a real difference. The observed maximum
relative deviation over every spread point is ~3.3e-16 (about 1.5 ulp).

The one genuine divergence is the ``VoronovEIArrhenius`` A-factor unit setter for
the ``m^3/(molecule*s)`` spelling: ``origin/99`` relabels metre as centimetre and
understates the rate by exactly 10**6, while ``plasma`` normalizes it correctly
(multiply by Avogadro's number). :meth:`test_voronov_A_setter_divergence` pins
that divergence with ``plasma`` on the dimensionally-correct side.

The table-driven (Z,N) cases require the plasma RMG-database and are marked
``database``; they check that the Badnell/Voronov data tables — one physical file
read by both branches through ``settings['database.directory']`` — parse to the
same rate on both.
"""
import hashlib
import importlib.util
import json
import math
import os

import pytest

import rmgpy.constants as constants
import rmgpy.kinetics.arrhenius as arrhenius_module
from rmgpy.kinetics.arrhenius import BadnellRRArrhenius, VoronovEIArrhenius

_HERE = os.path.dirname(os.path.abspath(__file__))

# The generator/shared-parameters module sits beside this test under ``test/`` and
# is not part of the installed ``rmgpy`` package, so load it by path.
_gen_spec = importlib.util.spec_from_file_location(
    "plasma_kinetics_differential_gen",
    os.path.join(_HERE, "plasma_kinetics_differential_gen.py"),
)
gen = importlib.util.module_from_spec(_gen_spec)
_gen_spec.loader.exec_module(gen)

_FIXTURE_PATH = os.path.join(_HERE, "plasma_kinetics_origin99_reference.json")

with open(_FIXTURE_PATH) as _fh:
    REFERENCE = json.load(_fh)

# Live values from the build under test, computed once. YAML cases are excluded
# here so the core comparison needs no database; they are exercised separately
# under the ``database`` marker.
LIVE = gen.evaluate(include_yaml=False)

RTOL = 1e-12

# Rate-coefficient cases that are numeric on BOTH branches across the whole spread.
RATE_CASES = [
    "TwoTemp_grc",
    "TwoTemp_two_temp",
    "ECP_grc",
    "ECP_electron_temp",
    "Badnell_grc",
    "Voronov_grc",
]

# Expected sha256 of the data tables the fixture's YAML cases were generated from.
_BADNELL_SHA = REFERENCE["meta"]["provenance"]["badnell_yaml_sha256"]
_VORONOV_SHA = REFERENCE["meta"]["provenance"]["voronov_yaml_sha256"]


def _relerr(a, b):
    if a == b:
        return 0.0
    denom = max(abs(a), abs(b))
    return abs(a - b) / denom if denom else 0.0


class PlasmaKineticsConsolidationTest:
    """``plasma`` reproduces ``origin/99`` for the four plasma kinetics classes."""

    def test_reference_and_live_share_the_spread(self):
        """The two sides sample the same temperatures (a guard on the comparison)."""
        for key in ("spread_T", "spread_Te", "grid_T", "grid_Te"):
            assert LIVE["meta"][key] == REFERENCE["meta"][key], key
        assert LIVE["meta"]["Na"] == REFERENCE["meta"]["Na"]

    @pytest.mark.parametrize("case", RATE_CASES)
    def test_rate_paths_reproduce_origin99(self, case):
        """k(T)/k(Te)/k(T,Te) match origin/99 to machine precision across the spread."""
        live, ref = LIVE[case], REFERENCE[case]
        keys = set(live) & set(ref)
        assert keys, f"no shared points for {case}"
        worst = 0.0
        for k in keys:
            lv, rv = live[k], ref[k]
            assert not isinstance(lv, str), f"{case}[{k}] live errored: {lv}"
            assert not isinstance(rv, str), f"{case}[{k}] ref errored: {rv}"
            worst = max(worst, _relerr(lv, rv))
        assert worst <= RTOL, f"{case}: max relerr {worst:.3e} exceeds {RTOL:.0e}"

    def test_scalar_A_normalization_matches(self):
        """The SI pre-exponential agrees for the shipped cm^3/(molecule*s) spelling."""
        for key in ("TwoTemp_A_si", "Badnell_A_si"):
            assert math.isclose(LIVE[key], REFERENCE[key], rel_tol=RTOL), key

    def test_voronov_A_setter_divergence(self):
        """The single genuine divergence, resolved in favour of the consolidated branch.

        For every A-factor spelling except ``m^3/(molecule*s)`` the two branches
        already agree. For ``m^3/(molecule*s)`` origin/99 understates by exactly
        10**6 (it relabels metre as centimetre); the consolidated branch
        normalizes correctly to ``A * Avogadro``. This test would go red if the
        divergence were papered over on either side.
        """
        live = LIVE["Voronov_A_setter"]
        ref = REFERENCE["Voronov_A_setter"]

        # Every spelling except the divergent one agrees between the branches.
        for spelling in ["cm^3/(molecule*s)", "cm^3/(mol*s)", "m^3/(mol*s)",
                         "cm^3/s", "m^3/s"]:
            lv, rv = live[spelling]["A_si"], ref[spelling]["A_si"]
            assert math.isclose(lv, rv, rel_tol=RTOL), \
                f"{spelling}: live {lv!r} vs ref {rv!r} unexpectedly differ"

        # The divergence itself.
        spelling = "m^3/(molecule*s)"
        live_a = live[spelling]["A_si"]
        ref_a = ref[spelling]["A_si"]
        # Consolidated branch is dimensionally correct: 1e-8 m^3/(molecule*s) * Na.
        assert math.isclose(live_a, 1.0e-8 * constants.Na, rel_tol=RTOL), \
            f"consolidated A_si {live_a!r} is not 1e-8 * Na"
        # origin/99 understated by exactly 10**6.
        assert math.isclose(ref_a * 1.0e6, live_a, rel_tol=1e-9), \
            f"expected origin/99 to understate by 1e6: ref {ref_a!r}, live {live_a!r}"
        assert live_a / ref_a > 1.0e5, "divergence collapsed — the branches now agree here"

    def test_electron_temp_delegation(self):
        """The consolidated-branch-only ``get_rate_coefficient_electron_temp``
        delegates exactly to ``get_rate_coefficient`` for Badnell and Voronov."""
        brr = BadnellRRArrhenius(**gen.BADNELL_KW)
        for Te in gen.SPREAD_T:
            assert math.isclose(brr.get_rate_coefficient_electron_temp(Te),
                                brr.get_rate_coefficient(Te), rel_tol=RTOL), Te
        vei = VoronovEIArrhenius(**gen.VORONOV_KW)
        for Te in gen.SPREAD_TE:
            assert math.isclose(vei.get_rate_coefficient_electron_temp(Te),
                                vei.get_rate_coefficient(Te), rel_tol=RTOL), Te

    @pytest.mark.database
    def test_table_driven_rates_reproduce_origin99(self):
        """Badnell/Voronov built from the (Z,N) YAML tables match origin/99.

        Both branches read one physical table file through
        ``settings['database.directory']``; this confirms table + parser produce
        the same rate. Skipped unless the plasma tables (by sha256) are present.
        """
        from rmgpy.data.base import Database  # noqa: F401 - ensures settings importable
        from rmgpy import settings

        kinetics_dir = os.path.join(settings["database.directory"], "kinetics")
        badnell = os.path.join(kinetics_dir, "badnell.yaml")
        voronov = os.path.join(kinetics_dir, "voronov.yaml")
        for path, sha in ((badnell, _BADNELL_SHA), (voronov, _VORONOV_SHA)):
            if not os.path.exists(path):
                pytest.skip(f"plasma data table not present: {path}")
            got = hashlib.sha256(open(path, "rb").read()).hexdigest()
            if got != sha:
                pytest.skip(f"data table differs from fixture provenance: {path}")

        for (Z, N) in gen.BADNELL_ZN:
            ref = REFERENCE[f"Badnell_yaml_{Z}_{N}"]
            b = BadnellRRArrhenius(Z=Z, N=N)
            for T_str, rv in ref.items():
                lv = b.get_rate_coefficient(float(T_str))
                assert math.isclose(lv, rv, rel_tol=RTOL), f"Badnell ({Z},{N}) @ {T_str}"

        for (Z, N) in gen.VORONOV_ZN:
            ref = REFERENCE[f"Voronov_yaml_{Z}_{N}"]
            v = VoronovEIArrhenius(Z=Z, N=N)
            for Te_str, rv in ref.items():
                lv = v.get_rate_coefficient(float(Te_str))
                assert math.isclose(lv, rv, rel_tol=RTOL), f"Voronov ({Z},{N}) @ {Te_str}"

    def test_module_is_the_build_under_test(self):
        """A live guard that the compiled arrhenius module was actually imported."""
        assert arrhenius_module.__file__.endswith(".so") or \
            arrhenius_module.__file__.endswith(".pyx"), arrhenius_module.__file__
