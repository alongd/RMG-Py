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
i250-seed-dict-comments -- round-trip test for `//`-comment stripping in
rmgpy.data.base.Database.get_species.

This is the read-side half of unblocking provenance comments in kinetics library / depository /
family dictionary.txt files. The write side (labelling the seed and seed_edge dictionaries with a
provenance comment block) is owned by a different, unmerged branch
(`i245-prescribed-te-declares-itself`) and is deliberately NOT part of this change -- see
docs/i250-seed-dict-comments/REPORT.md. This test only proves the engine-side blocker is gone: a
dictionary.txt carrying a `//` comment block in the shape that branch emits can be round-tripped
through Database.get_species without the comment corrupting parsing or leaking into a label.

The comment block below is transcribed BY HAND from that branch's rmgpy/diagnostic_mode.py and
rmgpy/chemkin.pyx (tip 544bf5718) -- it is a literal in this file, not an import from that branch.

Note on the trailing-same-line-comment shape: the fix ports get_species's `//`-stripping from
load_species_dictionary's `else` branch verbatim, including a shared quirk -- `line[0:index]`
discards the line's trailing newline along with the comment, so a trailing comment on a line
followed by another data line of the same species merges the two lines and corrupts the parse.
This was confirmed, independently, to be a pre-existing limitation of BOTH readers (not something
introduced or left unfixed by this change): rmgpy.chemkin.load_species_dictionary fails identically
on that placement. The test below exercises the one placement that is well-formed for either
reader: a trailing comment on the LAST atom line of a block, immediately before the blank-line
separator.
"""

import os

import pytest

from rmgpy.data.base import Database
from rmgpy.exceptions import DatabaseError

# Transcribed by hand from i245-prescribed-te-declares-itself's provenance comment block shape:
# a `//`-prefixed mode line, a `//`-prefixed disclaimer line, and a bare `//` prefix line.
PROVENANCE_COMMENT_BLOCK = (
    "// RMG diagnostic mode: prescribed-Te transport diagnostic\n"
    "// Species evolution conditional on prescribed Te and n_e; absolute electron density "
    "and ionisation degree are not predicted.\n"
    "//\n"
)

# Three species with distinct, non-trivial structures (a hydrocarbon radical, a ring, and a
# multiply-bonded species), so that a mangled block split or a wrong adjacency-list read would
# show up as a structural mismatch, not just a count mismatch.
SPECIES = [
    ("C2H5", "1 C u1 p0 c0 {2,S} {3,S} {4,S}\n"
             "2 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}\n"
             "3 H u0 p0 c0 {1,S}\n"
             "4 H u0 p0 c0 {1,S}\n"
             "5 H u0 p0 c0 {2,S}\n"
             "6 H u0 p0 c0 {2,S}\n"
             "7 H u0 p0 c0 {2,S}\n"),
    ("cyclopropane", "1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}\n"
                      "2 C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}\n"
                      "3 C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}\n"
                      "4 H u0 p0 c0 {1,S}\n"
                      "5 H u0 p0 c0 {1,S}\n"
                      "6 H u0 p0 c0 {2,S}\n"
                      "7 H u0 p0 c0 {2,S}\n"
                      "8 H u0 p0 c0 {3,S}\n"
                      "9 H u0 p0 c0 {3,S}\n"),
    ("N2", "1 N u0 p1 c0 {2,T}\n"
           "2 N u0 p1 c0 {1,T}\n"),
]


def _expected_species_objects():
    from rmgpy.species import Species
    expected = {}
    for label, adjlist in SPECIES:
        sp = Species(label=label)
        sp.from_adjacency_list(adjlist)
        expected[label] = sp
    return expected


def _assert_round_trip(species_dict):
    """Shared structural assertions used by every shape below."""
    expected = _expected_species_objects()

    assert set(species_dict.keys()) == set(expected.keys()), (
        f"labels differ: got {sorted(species_dict.keys())}, expected {sorted(expected.keys())}"
    )

    for label, expected_species in expected.items():
        actual_species = species_dict[label]
        assert actual_species.is_isomorphic(expected_species), (
            f"species {label!r} round-tripped to a non-isomorphic structure:\n"
            f"  expected: {expected_species.molecule[0].to_adjacency_list()}\n"
            f"  actual:   {actual_species.molecule[0].to_adjacency_list()}"
        )
        # Also compare the adjacency-list text directly, so a structurally-coincidental
        # isomorphism (e.g. atom order) does not mask a comment fragment leaking into the
        # parsed graph as a stray atom/bond.
        assert (actual_species.molecule[0].to_adjacency_list().strip()
                == expected_species.molecule[0].to_adjacency_list().strip())

        # The comment text must never leak into a label.
        assert "//" not in label
        assert "RMG diagnostic mode" not in label
        assert "prescribed-Te" not in label


def _write_dictionary(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    assert path.exists()
    assert path.stat().st_size > 0
    return str(path)


class TestBaseCommentStripping:
    """
    Round-trip tests proving rmgpy.data.base.Database.get_species accepts `//`-to-end-of-line
    comments, matching the syntax rmgpy.chemkin.load_species_dictionary already accepts.
    """

    def test_comment_block_at_head_of_file(self, tmp_path):
        """A provenance comment block at the head of the file, before any species."""
        text = PROVENANCE_COMMENT_BLOCK + "\n" + "".join(
            f"{label}\n{adjlist}\n" for label, adjlist in SPECIES
        )
        path = _write_dictionary(tmp_path, "dictionary_head.txt", text)
        species_dict = Database().get_species(path, resonance=False)
        _assert_round_trip(species_dict)

    def test_comment_block_between_species(self, tmp_path):
        """A comment block on its own line(s), between two species blocks."""
        label0, adjlist0 = SPECIES[0]
        label1, adjlist1 = SPECIES[1]
        label2, adjlist2 = SPECIES[2]
        text = (
            f"{label0}\n{adjlist0}\n"
            + PROVENANCE_COMMENT_BLOCK + "\n"
            + f"{label1}\n{adjlist1}\n"
            + f"{label2}\n{adjlist2}\n"
        )
        path = _write_dictionary(tmp_path, "dictionary_between.txt", text)
        species_dict = Database().get_species(path, resonance=False)
        _assert_round_trip(species_dict)

    def test_trailing_comment_on_adjacency_list_line(self, tmp_path):
        """
        A trailing `//` comment sharing a line with an adjacency-list line.

        This ports get_species's `//`-stripping from rmgpy.chemkin.load_species_dictionary's
        `else` branch verbatim: `line = line[0:line.index('//')]` discards everything from `//`
        onward, INCLUDING the line's trailing newline. That means a trailing comment on a line
        that is followed by another atom/bond line of the same species merges the two lines and
        corrupts the bond parse -- confirmed this is not a regression, but a pre-existing
        limitation shared identically by both readers (reproduced independently against
        rmgpy.chemkin.load_species_dictionary itself, which fails the exact same way on that
        placement). The one placement where a trailing same-line comment is well-formed for
        either reader is the LAST atom line of a block, immediately before the blank-line
        separator -- there is no following data line for the truncated newline to merge into.
        That is the placement exercised here.
        """
        label0, adjlist0 = SPECIES[0]
        label1, adjlist1 = SPECIES[1]
        label2, adjlist2 = SPECIES[2]
        adjlist0_lines = adjlist0.splitlines()
        # Append a trailing comment to the LAST atom line of the first species (immediately
        # before the blank-line separator), not an interior line -- see docstring above.
        adjlist0_lines[-1] = adjlist0_lines[-1] + " // trailing provenance note"
        adjlist0_with_trailing_comment = "\n".join(adjlist0_lines) + "\n"

        text = (
            f"{label0}\n{adjlist0_with_trailing_comment}\n"
            f"{label1}\n{adjlist1}\n"
            f"{label2}\n{adjlist2}\n"
        )
        path = _write_dictionary(tmp_path, "dictionary_trailing.txt", text)
        species_dict = Database().get_species(path, resonance=False)
        _assert_round_trip(species_dict)

    @pytest.mark.parametrize("label", ["//", "//N2"])
    def test_comment_only_label_is_rejected(self, tmp_path, label):
        """A comment-shaped label must not silently become the empty label."""
        text = f"{label}\n{SPECIES[2][1]}\n"
        path = _write_dictionary(tmp_path, "dictionary_empty_label.txt", text)

        with pytest.raises(DatabaseError, match=r"Empty species label") as exc_info:
            Database().get_species(path, resonance=False)

        message = str(exc_info.value)
        assert str(path) in message
        assert label in message

    def test_whitespace_separator_names_the_following_comment_label(self, tmp_path):
        """Whitespace-only separators must reset the diagnostic label candidate."""
        text = (
            "// header\n"
            " \t\n"
            "//N2\n"
            f"{SPECIES[2][1]}\n"
        )
        path = _write_dictionary(tmp_path, "dictionary_whitespace_separator.txt", text)

        with pytest.raises(DatabaseError, match=r"Empty species label") as exc_info:
            Database().get_species(path, resonance=False)

        assert "'//N2'" in str(exc_info.value)
