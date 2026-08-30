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

    def test_execute_still_runs_the_export_step(self):
        """
        The export steps live in :meth:`RMG.generate_end_of_run_cantera_files`,
        which ``execute`` calls. Without that call the handlers below are dead
        code and the run exits 0 regardless of what they record.
        """
        import inspect

        from rmgpy.rmg.main import RMG

        source = inspect.getsource(RMG.execute)
        assert 'self.generate_end_of_run_cantera_files()' in source

    def test_the_export_steps_record_rather_than_swallow(self):
        """
        The other half of the wiring: the handlers that used to only log.

        Counting appends alone would pass a block that grew a third ``except``
        which only logged, so this checks the stronger property -- that every
        exception handler in the export step records what it caught.
        """
        import inspect
        import re

        from rmgpy.rmg.main import RMG

        source = inspect.getsource(RMG.generate_end_of_run_cantera_files)
        handlers = re.findall(r'^        except .*?(?=^        \S|\Z)', source,
                              flags=re.MULTILINE | re.DOTALL)
        assert handlers, 'the export step has no exception handlers at all'
        for handler in handlers:
            assert 'self.export_failures.append(' in handler, (
                'an export-step exception handler logs without recording:\n' + handler)


class TestCompletionBannerFollowsTheExportCheck:
    """The third defect: ``execute`` logged ``MODEL GENERATION COMPLETED`` *before*
    ``report_export_failures`` raised. A run that failed export therefore printed
    success and then a traceback, and any tool that greps for that banner saw only
    success. The banner must never be emitted for a run whose export failed."""

    def test_report_precedes_the_banner_in_execute(self):
        """``execute`` is too heavy to drive without a full model and database, so
        this checks the property directly on its source -- the same idiom the
        wiring tests above use. ``report_export_failures`` raises iff a failure was
        recorded, so placing its call before the banner makes the banner
        unreachable on the failure path while leaving the healthy path unchanged."""
        import inspect

        from rmgpy.rmg.main import RMG

        source = inspect.getsource(RMG.execute)
        assert 'self.report_export_failures()' in source
        assert 'MODEL GENERATION COMPLETED' in source
        assert source.index('self.report_export_failures()') < source.index('MODEL GENERATION COMPLETED'), (
            'the completion banner is emitted before the export-failure check that '
            'can invalidate it')

    def test_a_recorded_failure_raises_before_any_completion_banner(self, caplog):
        """Behavioral guard on the ordering: replaying ``execute``'s tail
        (report, then banner) with a failure recorded must raise before the banner
        is ever logged."""
        import logging

        from rmgpy.rmg.main import RMG

        rmg = RMG(input_file=None, output_directory=None)
        rmg.export_failures.append((
            'Cantera export for plasma mechanism (cantera2/)',
            MechanismWriterError('No Cantera artifact was produced')))

        with caplog.at_level(logging.INFO):
            with pytest.raises(MechanismWriterError):
                # execute()'s fixed tail order: the check that can invalidate the
                # banner runs first and raises, so the banner is never reached.
                rmg.report_export_failures()
                logging.info('MODEL GENERATION COMPLETED')
        assert 'MODEL GENERATION COMPLETED' not in caplog.text, (
            'the completion banner was logged for a run whose export failed')
