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
I-135 -- a failed export no longer reports success.

RMG's end-of-run Chemkin-to-Cantera translation caught every exception, logged
it, and carried on. ``execute`` then returned normally and the process exited 0,
so a run that left ``cantera_from_ck/`` empty told its caller that everything
had worked -- the same failure shape that has already let green verifications
sit on top of a broken pipeline in this project.

The failure is now recorded as it is caught, reported in a banner at the end of
the run, and raised once everything that did work is safely on disk.
"""

import pytest

from rmgpy.exceptions import MechanismWriterError


class TestAFailedExportIsNotSilent:
    """The second defect: the end-of-run Chemkin-to-Cantera translation caught
    every exception, logged it, and let the process exit 0 -- so a run that
    produced no Cantera file at all reported success to its caller."""

    def _rmg(self):
        from rmgpy.rmg.main import RMG
        return RMG(input_file=None, output_directory=None)

    def test_a_clean_run_records_no_failure_and_does_not_raise(self):
        rmg = self._rmg()
        assert rmg.export_failures == []
        assert rmg.report_export_failures() is None

    def test_a_recorded_failure_fails_the_run(self):
        rmg = self._rmg()
        rmg.export_failures.append(('Chemkin-to-Cantera translation (cantera_from_ck/)',
                                    ValueError('could not convert string to float')))

        with pytest.raises(MechanismWriterError) as exc:
            rmg.report_export_failures()
        assert 'cantera_from_ck' in str(exc.value)

    def test_execute_consults_the_record(self):
        """The wiring, not just the helper. Without this call the helper is
        dead code and the run still exits 0, which is the whole defect."""
        import inspect

        from rmgpy.rmg.main import RMG

        source = inspect.getsource(RMG.execute)
        assert 'self.report_export_failures()' in source

    def test_the_translation_step_records_rather_than_swallows(self):
        """The other half of the wiring: the handlers that used to only log."""
        import inspect

        from rmgpy.rmg.main import RMG

        source = inspect.getsource(RMG.execute)
        assert source.count('self.export_failures.append(') == 2, source.count('self.export_failures.append(')
