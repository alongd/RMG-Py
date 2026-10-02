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


import os


class QMFileNameError(ValueError):
    """A molecule's QM files exceed a target directory's filename byte limit."""


def _check_file_names(species, file_names):
    """Check the longest planned filename per directory before any QM writes.

    ``file_names`` contains (directory, filename) pairs. A directory that does
    not exist yet, including a future temporary directory represented by None,
    uses NAME_MAX=255 until its actual filesystem limit can be queried.
    """
    longest = {}
    for directory, filename in file_names:
        byte_count = len(os.fsencode(filename))
        if directory not in longest or byte_count > longest[directory][1]:
            longest[directory] = (filename, byte_count)
    for directory, (filename, byte_count) in longest.items():
        name_max = os.pathconf(directory, 'PC_NAME_MAX') if directory and os.path.isdir(directory) else 255
        if 0 <= name_max < byte_count:
            target = directory if directory is not None else 'temporary QM directory (not yet created)'
            raise QMFileNameError(
                "QM species {}: longest filename {!r} is {} bytes; NAME_MAX={} in target directory {!r}".format(
                    species, filename, byte_count, name_max, target))
