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
I-198 -- a mis-spelled or mis-formatted rmgrc must not silently load the default database.

``Settings.require_database_directory()`` exists to stop a run proceeding against a database the
user never named. It closed two doors -- no rmgrc found at all, and a ``database.directory`` naming
a non-existent directory -- but left a third open: an rmgrc that *was* found yet never set
``database.directory``, because the key was mis-spelled, written as an INI ``[section]``, or absent.
The loader's substring match finds nothing, the value stays the compiled-in default (a real
existing directory on most checkouts), and both prior checks pass. The run then measures nothing
against the wrong database -- the exact failure the guard exists to prevent.

These tests pin the guard directly. The INI-section and typo'd-key cases are the hole (RED before
the fix); the two existing refusals and the correctly-spelled happy path guard against the fix
degrading into a blanket refusal.
"""

import os

import pytest

from rmgpy import Settings, SettingsError


def _settings_from_rmgrc(tmp_path, contents):
    """Write ``contents`` to an rmgrc in ``tmp_path`` and load a Settings from it explicitly."""
    rmgrc = tmp_path / "rmgrc"
    rmgrc.write_text(contents)
    return Settings(path=str(rmgrc))


class I198RmgrcGuardTest:
    """require_database_directory() must refuse a found-but-unset rmgrc, and only that."""

    def test_ini_section_spelling_refused(self, tmp_path):
        """An INI ``[database]`` / ``directory =`` spelling sets nothing -- must raise, not default."""
        s = _settings_from_rmgrc(tmp_path, "[database]\ndirectory = /nonexistent/path/input\n")
        assert s.filename is not None  # the file WAS found
        assert s.sources["database.directory"].startswith("Default")  # yet nothing was set
        with pytest.raises(SettingsError):
            s.require_database_directory()

    def test_typoed_key_refused(self, tmp_path):
        """A typo in the key (``databse.directory``) sets nothing -- must raise, not default."""
        s = _settings_from_rmgrc(tmp_path, "databse.directory = /nonexistent/path/input\n")
        assert s.sources["database.directory"].startswith("Default")
        with pytest.raises(SettingsError):
            s.require_database_directory()

    def test_unrelated_key_refused(self, tmp_path):
        """An rmgrc with only unrelated content sets nothing -- must raise, not default."""
        s = _settings_from_rmgrc(tmp_path, "totally.bogus = 1\n")
        assert s.sources["database.directory"].startswith("Default")
        with pytest.raises(SettingsError):
            s.require_database_directory()

    def test_missing_file_still_refused(self, tmp_path):
        """Existing refusal 1: no rmgrc found at all (filename is None) must still raise."""
        s = _settings_from_rmgrc(tmp_path, "totally.bogus = 1\n")
        s.filename = None  # simulate the "searched everywhere, found nothing" state
        with pytest.raises(SettingsError):
            s.require_database_directory()

    def test_nonexistent_directory_still_refused(self, tmp_path):
        """Existing refusal 2: a database.directory naming a non-existent dir must still raise."""
        s = _settings_from_rmgrc(tmp_path, "database.directory = /this/path/does/not/exist/input\n")
        assert s.sources["database.directory"].startswith("from ")  # it WAS set from the file
        with pytest.raises(SettingsError):
            s.require_database_directory()

    def test_correctly_spelled_resolves(self, tmp_path):
        """A correctly-spelled rmgrc naming an existing dir must resolve normally, not raise."""
        db = tmp_path / "db_input"
        db.mkdir()
        s = _settings_from_rmgrc(tmp_path, "database.directory = {0}\n".format(db))
        assert s.sources["database.directory"].startswith("from ")
        resolved = s.require_database_directory()
        assert resolved == os.path.abspath(str(db))

    def test_programmatic_set_not_refused(self, tmp_path):
        """A value set programmatically is not the default and must not be refused.

        The discriminator is 'value IS the compiled-in default' (source == DEFAULT_SOURCE),
        not 'value came from a file'. ``__setitem__`` writes source '-', so keying on
        provenance would false-refuse an explicitly-set, correct value -- a false positive
        worse than the hole this guard closes.
        """
        db = tmp_path / "db_input"
        db.mkdir()
        # Start from an rmgrc that sets nothing (would otherwise be refused), then override
        # database.directory programmatically -- the writer path with source '-'.
        s = _settings_from_rmgrc(tmp_path, "totally.bogus = 1\n")
        s["database.directory"] = str(db)
        assert s.sources["database.directory"] == "-"
        assert s.sources["database.directory"] != s.DEFAULT_SOURCE
        resolved = s.require_database_directory()
        assert resolved == os.path.abspath(str(db))
