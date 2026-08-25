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
No plasma deck in this repository reaches ``Cation_R_Recombination``.

``Cation_R_Recombination`` is a lithium-ion-battery SEI family: Marcus electron transfer at a
Li(110) electrode in liquid ethylene carbonate at 298.15 K under declared electrode potentials.
Read as plasma kinetics its rates land 30-230 orders of magnitude below anything physical, and a
binding ruling excludes it from every plasma configuration while *preserving* it for the
electrochemical domain, where it is correct.

The exclusion rests on declared configuration alone on this branch, for two reasons that are both
true here:

* the kinetics database carries no ``plasma`` family set (``kinetics/families/recommended.py``
  defines ``default``, ``electrochem``, ``surface``, ... but no plasma set), so a deck that names
  the family individually reaches it regardless of any curated selection; and
* this repository carries no quarantine loader that would refuse the family at load time -- that
  hard-fail gate lives on a branch that has not merged.

So a plasma deck naming the family is a live path, not a stale reference, and these tests are what
holds the exclusion in place until the runtime gate lands.

Assertion styles used here, stated per test so the strength of each is visible:

* *parsed-configuration* -- the deck is executed by RMG's own input reader
  (:func:`rmgpy.rmg.input.read_input_file`) and the resolved ``kinetics_families`` list is
  inspected. Stronger than a grep over deck text: it sees the list RMG actually acts on, including
  any list a deck computes rather than spells out.
* *behavioural* -- the deck's own declared families and its own declared species are pushed
  through the real kinetics database, and the generated reactions are counted by family. This is
  the strongest form available without running a full RMG job: a deck that loads and produces zero
  reactions from the family cannot introduce electrode kinetics into a plasma mechanism.
* *AST-literal* -- the deck's ``database(kineticsFamilies=[...])`` argument is read statically.
  Used only for the electrochemical negative-control deck, which declares a ``liquidSurfaceReactor``
  and therefore cannot be executed in an environment without ReactionMechanismSimulator/Julia.
"""

import ast
import itertools
import os
from collections import Counter

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

#: The reclassified battery-SEI family. Preserved in the database as provenance evidence and
#: correct for electrochemistry; excluded from every plasma configuration.
QUARANTINED_FAMILY = "Cation_R_Recombination"

#: Deck-DSL call that makes a deck a plasma deck.
PLASMA_REACTOR_CALL = "plasmaReactor"

#: Directories that hold no RMG decks and are expensive or meaningless to walk.
SKIPPED_DIRS = {".git", "build", "dist", "external", "testing", ".pytest_cache", "__pycache__"}

#: The plasma decks known to exist when this test was written. Asserted to still be discovered so
#: that a rename or a broken discovery sweep cannot make the exclusion pass vacuously.
KNOWN_PLASMA_DECKS = {
    "docs/m7-preflight/input.py",
    "docs/m7-preflight/input_placement.py",
    "docs/m7-preflight/input_secondary.py",
}

#: The electrochemical deck that legitimately keeps the family (Li(110) electrode, liquid ACN,
#: 298.15 K, declared liquid/surface potentials). The negative control.
SEI_DECK = "examples/rmg/SEI_pure_ACN/input.py"


def _top_level_calls(path):
    """
    Return the set of names called at module top level in the file at `path`, or None if the file
    is not parseable Python.

    RMG input files are flat sequences of top-level DSL calls, so this identifies a deck and its
    reactor type without importing anything. Source and test modules that merely *define* or
    *invoke inside a function* the same names are not matched.
    """
    try:
        with open(path, "r") as f:
            tree = ast.parse(f.read(), filename=path)
    except (SyntaxError, ValueError, UnicodeDecodeError):
        return None
    names = set()
    for node in tree.body:
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
            func = node.value.func
            if isinstance(func, ast.Name):
                names.add(func.id)
    return names


def _declared_families_literal(path):
    """
    Return the literal `kineticsFamilies` argument of the top-level `database(...)` call in the
    deck at `path`, or None if it is absent or not a literal.
    """
    with open(path, "r") as f:
        tree = ast.parse(f.read(), filename=path)
    for node in tree.body:
        if not (isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)):
            continue
        func = node.value.func
        if not (isinstance(func, ast.Name) and func.id == "database"):
            continue
        for kw in node.value.keywords:
            if kw.arg == "kineticsFamilies":
                try:
                    return ast.literal_eval(kw.value)
                except ValueError:
                    return None
    return None


def _discover_decks():
    """
    Walk the repository and return every RMG input deck as a repo-relative path.

    A deck is a parseable Python file whose top-level calls include `database` and at least one
    call whose name contains `Reactor`.
    """
    decks = []
    for dirpath, dirnames, filenames in os.walk(REPO_ROOT):
        dirnames[:] = [d for d in dirnames if d not in SKIPPED_DIRS]
        for filename in filenames:
            if not filename.endswith(".py"):
                continue
            full = os.path.join(dirpath, filename)
            calls = _top_level_calls(full)
            if calls is None:
                continue
            if "database" in calls and any("Reactor" in name for name in calls):
                decks.append((os.path.relpath(full, REPO_ROOT), calls))
    return sorted(decks)


ALL_DECKS = _discover_decks()
PLASMA_DECKS = [relpath for relpath, calls in ALL_DECKS if PLASMA_REACTOR_CALL in calls]
NON_PLASMA_DECKS = [relpath for relpath, calls in ALL_DECKS if PLASMA_REACTOR_CALL not in calls]


def _read_deck(relpath, tmp_path):
    """Parse a deck with RMG's own input reader and return the populated RMG object."""
    from rmgpy.rmg.main import RMG
    from rmgpy.rmg.input import read_input_file

    full = os.path.join(REPO_ROOT, relpath)
    rmg = RMG(input_file=full, output_directory=str(tmp_path))
    read_input_file(full, rmg)
    return rmg


class PlasmaDeckDiscoveryTest:
    """
    Guard the sweep itself. An exclusion test that found no decks would pass for the wrong reason.
    """

    def test_plasma_decks_are_found(self):
        """At least one plasma deck exists, so the exclusion assertions are not vacuous."""
        assert PLASMA_DECKS, (
            "no plasma deck discovered under {0} -- the exclusion assertions below would pass "
            "vacuously. Either deck discovery broke or every plasma deck was removed.".format(REPO_ROOT)
        )

    def test_known_plasma_decks_still_discovered(self):
        """The decks this test was written against are still classified as plasma decks."""
        missing = KNOWN_PLASMA_DECKS - set(PLASMA_DECKS)
        assert not missing, (
            "these plasma decks are no longer discovered: {0}. If they were renamed or deleted, "
            "update KNOWN_PLASMA_DECKS; do not leave the sweep unable to see them.".format(sorted(missing))
        )

    def test_sei_deck_is_not_classified_as_plasma(self):
        """The electrochemical deck is discovered, and it is NOT a plasma deck."""
        relpaths = [relpath for relpath, _ in ALL_DECKS]
        assert SEI_DECK in relpaths, "the electrochemical negative-control deck was not discovered"
        assert SEI_DECK in NON_PLASMA_DECKS, (
            "{0} declares a plasmaReactor -- it is not the electrochemical deck this test assumes".format(SEI_DECK)
        )


class PlasmaDeckFamilyExclusionTest:
    """No plasma deck declares the battery-SEI family."""

    @pytest.mark.parametrize("relpath", PLASMA_DECKS)
    def test_plasma_deck_does_not_declare_family(self, relpath, tmp_path):
        """
        Assertion style: *parsed-configuration*. RMG's own input reader executes the deck and the
        resolved family list is inspected -- not the deck's text.
        """
        rmg = _read_deck(relpath, tmp_path)
        families = rmg.kinetics_families
        assert isinstance(families, list), (
            "{0} selects families with the sentinel {1!r} rather than an explicit list; this test "
            "cannot tell what that resolves to without loading the database".format(relpath, families)
        )
        assert QUARANTINED_FAMILY not in families, (
            "plasma deck {0} declares {1}, a lithium-ion-battery SEI family excluded from every "
            "plasma configuration. Declared families: {2}".format(relpath, QUARANTINED_FAMILY, families)
        )

    @pytest.mark.database
    @pytest.mark.parametrize("relpath", PLASMA_DECKS)
    def test_plasma_deck_generates_no_reaction_from_family(self, relpath, tmp_path):
        """
        Assertion style: *behavioural*. The deck's own declared families and its own declared
        reactive species are pushed through the real kinetics database over every unimolecular and
        bimolecular reactant set, and the generated reactions are counted by family. Zero reactions
        from the family is the assertion, backed by the family being absent from what the deck's
        declaration actually loaded (which also covers a deck naming a family *set*, or ``all``).
        """
        from rmgpy import settings
        from rmgpy.data.kinetics.database import KineticsDatabase

        rmg = _read_deck(relpath, tmp_path)
        families = rmg.kinetics_families
        assert isinstance(families, list), "{0} does not declare an explicit family list".format(relpath)

        depositories = rmg.kinetics_depositories
        kinetics_database = KineticsDatabase()
        kinetics_database.load(
            os.path.join(settings["database.directory"], "kinetics"),
            families=families,
            depositories=depositories if depositories is not None else [],
        )

        # A deck may name a family *set* rather than a family -- `kineticsFamilies=['electrochem']`
        # is legal and resolves to a group that contains this family, so the declaration check above
        # cannot see it. Assert on what the database actually loaded for this deck.
        assert QUARANTINED_FAMILY not in kinetics_database.families, (
            "plasma deck {0} resolves to a family set that loads {1}. Declared: {2}".format(
                relpath, QUARANTINED_FAMILY, families
            )
        )

        # generate_reactions_from_families takes ONE OR TWO reactants -- an exact reactant set, not
        # a species pool. Handing it the whole declared species list generates the reactions of that
        # single N-body combination, which for N > 2 is nothing at all; that is a test that passes
        # for the wrong reason. Enumerate the unimolecular and bimolecular combinations (self-pairs
        # included) that RMG itself would explore from the deck's reactive species.
        reactive = [spc for spc in rmg.initial_species if spc.reactive]
        assert reactive, "{0} declares no reactive species; nothing would be generated".format(relpath)
        reactant_sets = [[spc] for spc in reactive]
        reactant_sets += [[a, b] for a, b in itertools.combinations_with_replacement(reactive, 2)]

        reactions = []
        for reactant_set in reactant_sets:
            reactions.extend(
                kinetics_database.generate_reactions_from_families(reactant_set, only_families=None)
            )
        by_family = Counter(reaction.family for reaction in reactions)
        assert by_family[QUARANTINED_FAMILY] == 0, (
            "plasma deck {0} generated {1} reaction(s) from {2} out of its own declared feed. "
            "Generated by family: {3}".format(relpath, by_family[QUARANTINED_FAMILY], QUARANTINED_FAMILY, dict(by_family))
        )


class ElectrochemicalFamilyPreservedTest:
    """
    Negative control. The exclusion must be surgical: the family is still correct for the
    electrochemical/SEI domain and a binding ruling requires it be preserved there.
    """

    def test_sei_deck_still_declares_family(self):
        """
        Assertion style: *AST-literal*. The SEI deck declares a ``liquidSurfaceReactor`` and cannot
        be executed without ReactionMechanismSimulator/Julia, so its family list is read statically.
        """
        families = _declared_families_literal(os.path.join(REPO_ROOT, SEI_DECK))
        assert families is not None, "could not read kineticsFamilies from {0}".format(SEI_DECK)
        assert QUARANTINED_FAMILY in families, (
            "{0} no longer declares {1}. The family is correct for electrochemistry and must be "
            "preserved there; this ticket excludes it from plasma decks only.".format(SEI_DECK, QUARANTINED_FAMILY)
        )

    @pytest.mark.database
    def test_family_still_generates_the_sei_reaction(self):
        """
        Assertion style: *behavioural*. The family still loads from the database and still produces
        its reaction from an organolithium, so the plasma-side exclusion removed configuration and
        not the family.
        """
        from rmgpy import settings
        from rmgpy.data.kinetics.database import KineticsDatabase
        from rmgpy.species import Species

        kinetics_database = KineticsDatabase()
        kinetics_database.load(
            os.path.join(settings["database.directory"], "kinetics"),
            families=[QUARANTINED_FAMILY],
            depositories=["training"],
        )
        assert QUARANTINED_FAMILY in kinetics_database.families, (
            "{0} did not load from the database -- it must be preserved, not deleted".format(QUARANTINED_FAMILY)
        )

        methyllithium = Species(label="CH3Li").from_smiles("C[Li]")
        methyllithium.generate_resonance_structures()
        reactions = kinetics_database.generate_reactions_from_families([methyllithium], only_families=None)
        from_family = [reaction for reaction in reactions if reaction.family == QUARANTINED_FAMILY]
        assert from_family, (
            "{0} generated no reaction from CH3Li. The electrochemical path can no longer reach "
            "the family, which is a wider break than this ticket intends.".format(QUARANTINED_FAMILY)
        )
