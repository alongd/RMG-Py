#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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
The ``Marcus`` work terms ``wr`` and ``wp``.

``Marcus`` has always declared ``wr`` (work to assemble the precursor complex) and ``wp`` (work to
assemble the successor complex), parsed them out of the database, carried them through ``__repr__``,
``__reduce__``, rule averaging, and the RMS export -- and never let either one touch a rate.
``get_gibbs_activation_energy`` computed the bare quadratic ``lmbd/4 * (1 + dGrxn/lmbd)^2``, so an
author who wrote ``wr=(50,'kJ/mol')`` on an entry got a rate bit-identical to ``wr=0``. Nothing
warned, nothing failed: the stored value and its absence were indistinguishable.

That mattered here because the family that motivated this ticket models an overall rate as an
association step times an electron-transfer step, and the association contribution is precisely
what ``wr``/``wp`` carry. A silent zero there is a missing physical term wearing the costume of a
measured one.

The tests below pin three things:

* a supplied work term now moves the barrier by the amount Marcus theory says it should
  (:class:`TestWorkTermsEnterTheBarrier`);
* the terms cannot leak into ReactionMechanismSimulator, whose ``Marcus`` struct declares them and
  whose rate body never reads them, without a loud refusal
  (:class:`TestWorkTermsRefusedOnRmsExport`);
* every entry that leaves the work terms at zero -- which is every entry in the database today --
  is numerically untouched (:class:`TestZeroWorkTermsAreANoOp`).
"""

import numpy as np
import pytest

from rmgpy.exceptions import MechanismWriterError
from rmgpy.kinetics.arrhenius import Marcus, check_marcus_work_terms_exportable

#: Total reorganization energy, J/mol. Carried entirely on the inner sphere so that ``lmbd`` is
#: temperature-independent and the arithmetic below stays checkable by hand.
LAMBDA = 1.0e5

#: Gas constant, J/(mol*K). Local so a failure here cannot be masked by a constants-module change.
R = 8.314472


def marcus(wr=0.0, wp=0.0):
    """A temperature-independent ``Marcus`` model with total reorganization energy ``LAMBDA``."""
    return Marcus(
        A=(1.0e10, "m^3/(mol*s)"),
        n=0.0,
        lmbd_i_coefs=[LAMBDA, 0.0, 0.0, 0.0],
        beta=(1.2e10, "1/m"),
        wr=(wr, "J/mol"),
        wp=(wp, "J/mol"),
        lmbd_o=(0.0, "J/mol"),
    )


def expected_barrier(dGrxn, wr=0.0, wp=0.0):
    """The Marcus barrier with work terms, written out independently of the implementation."""
    return wr + LAMBDA / 4.0 * (1.0 + (dGrxn + wp - wr) / LAMBDA) ** 2


class TestWorkTermsEnterTheBarrier:
    """
    RED-before: every assertion in this class passed trivially with ``ratio == 1.0`` when the work
    terms were discarded, because the two models it compares were numerically the same model.
    """

    def test_nonzero_wr_raises_the_barrier(self):
        """``wr`` is work spent before the barrier, so it adds to the barrier."""
        T, dGrxn, wr = 1000.0, -2.0e4, 5.0e4
        bare = marcus().get_gibbs_activation_energy(T, dGrxn)
        worked = marcus(wr=wr).get_gibbs_activation_energy(T, dGrxn)
        assert worked > bare
        assert worked == pytest.approx(expected_barrier(dGrxn, wr=wr), rel=1e-12)

    def test_nonzero_wr_changes_the_rate(self):
        """The defect stated as a number: this ratio was exactly 1.0 before the fix."""
        T, dGrxn, wr = 1000.0, -2.0e4, 5.0e4
        k_bare = marcus().get_rate_coefficient(T, dGrxn)
        k_worked = marcus(wr=wr).get_rate_coefficient(T, dGrxn)
        assert k_worked != k_bare
        expected_ratio = np.exp(
            -(expected_barrier(dGrxn, wr=wr) - expected_barrier(dGrxn)) / (R * T)
        )
        assert k_worked / k_bare == pytest.approx(expected_ratio, rel=1e-9)

    def test_wp_shifts_the_driving_force(self):
        """
        ``wp`` enters only through ``dG' = dGrxn + wp - wr``, so a reaction with ``wp`` set is the
        same rate as one whose driving force is that much less favourable and whose ``wp`` is zero.
        """
        T, dGrxn, wp = 1000.0, -2.0e4, 3.0e4
        with_wp = marcus(wp=wp).get_rate_coefficient(T, dGrxn)
        folded_in = marcus().get_rate_coefficient(T, dGrxn + wp)
        assert with_wp == pytest.approx(folded_in, rel=1e-12)

    def test_equal_work_terms_do_not_cancel(self):
        """
        ``wr`` and ``wp`` cancel inside the quadratic but ``wr`` also sits outside it, so setting
        both equal is not a no-op -- it raises the barrier by exactly ``wr``. Getting this wrong is
        the natural way to mis-wire the terms, and it would leave the rate unchanged again.
        """
        T, dGrxn, w = 1000.0, -2.0e4, 4.0e4
        bare = marcus().get_gibbs_activation_energy(T, dGrxn)
        both = marcus(wr=w, wp=w).get_gibbs_activation_energy(T, dGrxn)
        assert both - bare == pytest.approx(w, rel=1e-12)

    def test_work_terms_break_identity(self):
        """
        Two models differing only in their work terms are now different rate laws, so
        ``is_identical_to`` must say so. It could not before: it compared an ``electrons``
        attribute ``Marcus`` has never defined and raised ``AttributeError`` on every pair.
        """
        assert marcus().is_identical_to(marcus())
        assert not marcus().is_identical_to(marcus(wr=5.0e4))
        assert not marcus().is_identical_to(marcus(wp=5.0e4))


class TestWorkTermsRefusedOnRmsExport:
    """
    ReactionMechanismSimulator's ``Marcus`` struct carries ``wr::K = 0.0`` and ``wp::K = 0.0`` and
    its rate body reads neither. Now that RMG's barrier does read them, exporting a non-zero work
    term would make the two runtimes compute different rates for the same reaction with nothing
    reporting the disagreement. Refuse instead.
    """

    def test_nonzero_wr_is_refused(self):
        with pytest.raises(MechanismWriterError) as exc:
            check_marcus_work_terms_exportable(marcus(wr=5.0e4), "an RMS YAML file", "reaction R1")
        message = str(exc.value)
        assert "reaction R1" in message
        assert "an RMS YAML file" in message
        assert "wr=50000" in message
        assert "never reads them" in message

    def test_nonzero_wp_is_refused(self):
        with pytest.raises(MechanismWriterError):
            check_marcus_work_terms_exportable(marcus(wp=1.0), "an RMS YAML file", "reaction R1")

    def test_zero_work_terms_export(self):
        """Zero is the one value both runtimes provably agree on, and it must pass silently."""
        assert check_marcus_work_terms_exportable(marcus(), "an RMS YAML file", "reaction R1") is None

    def test_non_marcus_kinetics_are_not_the_guard_s_business(self):
        assert check_marcus_work_terms_exportable(None, "an RMS YAML file", "reaction R1") is None

    def test_yaml_rms_writer_refuses(self):
        """The guard is actually reached from the RMS YAML writer, not merely importable."""
        from rmgpy.yaml_rms import obj_to_dict

        assert obj_to_dict(marcus(), [])["type"] == "Marcus"
        with pytest.raises(MechanismWriterError):
            obj_to_dict(marcus(wr=5.0e4), [])


class TestZeroWorkTermsAreANoOp:
    """
    The negative control. Every ``Marcus`` entry in the database carries ``wr=0`` and ``wp=0``, so
    the change must not move a single one of their rates.
    """

    @pytest.mark.parametrize("T", [300.0, 800.0, 1500.0, 2500.0])
    @pytest.mark.parametrize("dGrxn", [-2.0e5, -5.0e4, 0.0, 5.0e4])
    def test_barrier_reduces_to_the_bare_quadratic(self, T, dGrxn):
        kin = marcus()
        lmbd = kin.get_lmbd_i(T) + kin.lmbd_o.value_si
        bare = lmbd / 4.0 * (1.0 + dGrxn / lmbd) ** 2
        assert kin.get_gibbs_activation_energy(T, dGrxn) == bare

    @pytest.mark.parametrize("T", [300.0, 800.0, 1500.0, 2500.0])
    @pytest.mark.parametrize("dGrxn", [-2.0e5, -5.0e4, 0.0, 5.0e4])
    def test_rate_reduces_to_the_bare_expression(self, T, dGrxn):
        kin = marcus()
        lmbd = kin.get_lmbd_i(T) + kin.lmbd_o.value_si
        bare = lmbd / 4.0 * (1.0 + dGrxn / lmbd) ** 2
        expected = kin.A.value_si * T ** kin.n.value_si * np.exp(-bare / (R * T))
        assert kin.get_rate_coefficient(T, dGrxn) == pytest.approx(expected, rel=1e-12)

    def test_defaulted_work_terms_are_zero(self):
        """An entry that names no work terms gets zeros, not an unanswered question."""
        kin = Marcus(A=(1.0e10, "m^3/(mol*s)"), n=0.0, lmbd_i_coefs=[LAMBDA, 0.0, 0.0, 0.0])
        assert kin.wr.value_si == 0.0
        assert kin.wp.value_si == 0.0

    def test_round_trip_preserves_work_terms(self):
        """``repr`` round-trip: the terms survive serialisation, as they always did."""
        kin = marcus(wr=5.0e4, wp=3.0e4)
        clone = eval(repr(kin), {"Marcus": Marcus, "array": np.array, "np": np})
        assert clone.wr.value_si == pytest.approx(5.0e4)
        assert clone.wp.value_si == pytest.approx(3.0e4)
        assert clone.get_rate_coefficient(1000.0, -2.0e4) == pytest.approx(
            kin.get_rate_coefficient(1000.0, -2.0e4), rel=1e-12
        )
