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
This is the rmg module.
"""

import os
import os.path

from rmgpy.version import __version__
from rmgpy.exceptions import SettingsError


################################################################################


class Settings(dict):
    """
    A dictionary-like object containing global settings for RMG jobs. These
    settings are generally loaded from a file named ``rmgrc`` located at one
    of several possible places on disk. The ``rmgrc`` file used is stored in
    the `filename` attribute. This class inherits from the built-in dict class
    and adds methods for loading and resetting the settings, as well as a
    custom :meth:`__setitem__` method for processing the setting values before
    adding them to the dictionary.

    In general you should be working with the module-level variable
    ``settings`` in this module, which is an instance of this class.
    """

    #: Source marker written by :meth:`reset` for values left at their compiled-in
    #: default. ``load`` overwrites the marker with ``"from {filename}"`` on a match and
    #: :meth:`__setitem__` writes ``"-"`` for a programmatic assignment, so equality
    #: against this constant is the one reliable test for "still the compiled-in default".
    #: :meth:`require_database_directory` keys its refusal on that, not on provenance.
    DEFAULT_SOURCE = "Default, relative to RMG-Py source code"

    def __init__(self, path=None):
        super(Settings, self).__init__()
        self.filename = None
        self.sources = dict()
        self.load(path)

    def __setitem__(self, key, value):
        if key == "database.directory":
            value = os.path.abspath(os.path.expandvars(value))
        elif key == "test_data.directory":
            value = os.path.abspath(os.path.expandvars(value))
        else:
            raise SettingsError('Unexpecting setting "{0}" encountered.'.format(key))
        self.sources[key] = "-"
        super(Settings, self).__setitem__(key, value)

    def report(self):
        """
        Returns a string saying what is set and where things came from, suitable for logging
        """
        lines = ["Global RMG Settings:"]
        for key in self.keys():
            lines.append("   {0:20s} = {1:20s} ({2})".format(key, self[key], self.sources[key]))
        return "\n".join(lines)

    def load(self, path=None):
        """
        Load settings from a file on disk. If an explicit file is not specified,
        the following locations will be searched for a settings file, and the
        first one found will be loaded:

        * An rmgrc file in the current working directory

        * An rmgrc file in the user's $HOME/.rmg directory

        * An rmgrc file in the same directory as this package

        If none of these can be found, a SettingsError is raised.
        """
        # First set all settings to their default values
        self.reset()

        if path:
            # The user specified an explicit file to use for the settings
            # Make sure that it exists, fail if it does not
            if not os.path.exists(path):
                raise SettingsError('Specified RMG settings file "{0}" does not exist.'.format(path))
            else:
                self.filename = path
        else:
            # The user did not specify an explicit file to use for the settings
            # Load one of the default settings files instead
            working_dir = os.path.abspath(os.path.dirname(__file__))
            if os.path.exists("rmgrc"):
                self.filename = "rmgrc"
            elif os.path.exists(os.path.expanduser("~/.rmg/rmgrc")):
                self.filename = os.path.expanduser("~/.rmg/rmgrc")
            elif os.path.exists(os.path.join(working_dir, "rmgrc")):
                self.filename = os.path.join(working_dir, "rmgrc")
            else:
                return  # fail silently, instead of raising an error

        # From here on we assume that we have identified the appropriate
        # settings file to load

        with open(self.filename, "r") as f:
            for line in f:
                # Remove any comments from the line
                index = line.find("#")
                if index != -1:
                    line = line[:index]
                # Is there a key-value pair remaining?
                if line.find("database.directory") != -1:
                    value = line.split()[-1]  # Get the last token from this line
                    value = value.strip()
                    self["database.directory"] = value
                    self.sources["database.directory"] = "from {0}".format(self.filename)

                elif line.find("test_data.directory") != -1:
                    value = line.split()[-1]  # Get the last token from this line
                    value = value.strip()
                    self["test_data.directory"] = value
                    self.sources["test_data.directory"] = "from {0}".format(self.filename)

    def require_database_directory(self):
        """
        Resolve and validate ``database.directory`` for an actual run, returning the
        resolved absolute path. Two conditions are refused **loudly** rather than
        papered over with a silently-chosen default, because a silent default has
        repeatedly turned out to be the *wrong* database and produced a green test
        run that measured nothing:

        * No ``rmgrc`` configuration file was found at all (``self.filename is None``).
          RMG will not guess a database location; the run stops naming the file it
          wanted and how to create it.
        * An ``rmgrc`` *was* found but ``database.directory`` is still the compiled-in
          default, i.e. ``sources["database.directory"] == DEFAULT_SOURCE``. The rule the
          refusal enforces is exact: **no line in the file matched the loader's
          ``database.directory`` substring scan**, so nothing overrode the default. That
          covers an absent key, an INI ``[section]`` spelling, and a mis-spelling that does
          not contain the substring (e.g. ``databse.directory``). It does **not** cover a
          mis-spelling that still contains it -- ``database.directory_old = /x`` matches the
          scan, sets the source to ``from <file>``, and passes this refusal; such a line is
          the loose-substring-match weakness left out of scope, caught only by the existing
          directory check if the mis-parsed value is not an existing directory. This is the
          likeliest of the three refusals to be walked through: the default is usually a real
          directory, so the check below would pass and the run would proceed against a
          database the user never named.

          Keying on ``DEFAULT_SOURCE`` rather than on "did the value come from a file" is
          deliberate: a value set programmatically (``settings["database.directory"] = ...``,
          source ``"-"``) is not the default and must not be refused.
        * ``database.directory`` names a path that is not an existing directory.

        Raises :class:`SettingsError` in each case. This is intentionally *not*
        wired into settings loading (importing :mod:`rmgpy` must always succeed); it
        is called from the run path, where a missing or wrong database is fatal.

        Deliberate compatibility break: an ``rmgrc`` that sets only ``test_data.directory``
        and relied on the compiled-in ``database.directory`` default now hard-fails here.
        That reliance on a silently-chosen default is exactly what this guard exists to
        stop, so the break is intended, not incidental.
        """
        package_dir = os.path.abspath(os.path.dirname(__file__))
        repo_template = os.path.abspath(os.path.join(package_dir, "..", "rmgrc.template"))
        packaged_template = os.path.join(package_dir, "rmgrc_template")
        if self.filename is None:
            raise SettingsError(
                "No RMG configuration file (rmgrc) was found, so database.directory is unset.\n"
                "RMG refuses to guess a database location: a silently-chosen default has\n"
                "repeatedly been the wrong database and produced a meaningless-green run.\n"
                "\n"
                "Fix, from the repository root of this checkout:\n"
                "    cp rmgrc.template rmgrc\n"
                "then edit ./rmgrc so that database.directory points at your RMG-database\n"
                "input directory. Searched (in order): ./rmgrc, ~/.rmg/rmgrc, {packaged}\n"
                "Templates: {repo_template} or {packaged_template}".format(
                    packaged=os.path.join(package_dir, "rmgrc"),
                    repo_template=repo_template,
                    packaged_template=packaged_template,
                )
            )
        if self.sources["database.directory"] == self.DEFAULT_SOURCE:
            # A file was found (filename is not None) but no line in it matched the loader's
            # ``database.directory`` substring scan, so the value is still the compiled-in
            # default set by reset(). Equality against DEFAULT_SOURCE -- not
            # ``startswith("from ")`` -- is the discriminator, so a value assigned
            # programmatically (source "-") is correctly left alone rather than refused.
            raise SettingsError(
                "An RMG configuration file was found:\n"
                "    {filename}\n"
                "but no 'database.directory' line in it took effect, so the value is still\n"
                "the compiled-in default:\n"
                "    {path}\n"
                "RMG refuses to fall back to that default: a silently-chosen default has\n"
                "repeatedly been the wrong database and produced a meaningless-green run.\n"
                "An absent key, an INI-style [section] header, or a mis-spelled key is the\n"
                "usual cause.\n"
                "\n"
                "Fix: edit {filename} so it contains one line spelled exactly\n"
                "    database.directory = /path/to/your/RMG-database/input\n"
                "(a single 'key = value' line; no [section] headers). Template: {repo_template}\n"
                "or {packaged_template}.".format(
                    filename=self.filename,
                    path=self["database.directory"],
                    repo_template=repo_template,
                    packaged_template=packaged_template,
                )
            )
        path = self["database.directory"]
        if not os.path.isdir(path):
            raise SettingsError(
                "database.directory resolves to:\n"
                "    {path}\n"
                "which is not an existing directory. It was read from {source} ({filename}).\n"
                "RMG refuses to fall back to a default location.\n"
                "Fix: edit {filename} so database.directory names your RMG-database input\n"
                "directory, or copy {repo_template} and set it there.".format(
                    path=path,
                    source=self.sources["database.directory"],
                    filename=self.filename,
                    repo_template=repo_template,
                )
            )
        return path

    def reset(self):
        """
        Reset all settings to their default values.
        """
        self.filename = None
        rmgpy_module_dir = os.path.abspath(os.path.dirname(__file__))
        self["database.directory"] = os.path.realpath(os.path.join(rmgpy_module_dir, "..", "..", "RMG-database", "input"))
        self.sources["database.directory"] = self.DEFAULT_SOURCE
        self["test_data.directory"] = os.path.realpath(os.path.join(rmgpy_module_dir, "..", "test", "rmgpy", "test_data"))
        self.sources["test_data.directory"] = self.DEFAULT_SOURCE


# The global settings object
settings = Settings(path=None)

################################################################################


def get_path():
    """
    Return the directory that this file is found in on disk.
    """
    return os.path.abspath(os.path.dirname(__file__))
