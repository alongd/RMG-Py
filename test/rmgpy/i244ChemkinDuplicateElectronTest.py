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
The Chemkin writer must not stamp ``DUPLICATE`` on two reactions whose electron
placements differ.

Chemkin requires a duplicate pair to share a stoichiometry. RMG keeps a charged
reaction's free electrons out of ``Reaction.reactants``/``products`` and in the
scalar ``Reaction.electrons``, so ``mark_duplicate_reaction``'s reference
comparison sees only the heavy species -- and two reactions that differ only in
where their electrons sit compared equal there and were both marked. The deck
that came out contained, for one ionisation channel written by two owners:

    [Li] + e-  =>  [Lip] + e- + e-     DUPLICATE     2 reactants, 3 products
    [Li]       =>  [Lip] + e-          DUPLICATE     1 reactant,  2 products

Both carry ``electrons = +1``, so the NET scalar cannot tell them apart; only the
per-side placement pair -- ``(1, 2)`` against ``(0, 1)`` -- can. The warning is
what hid it: it prints the electron-free canonical form, in which the two really
are identical.

**Every assertion here is on the text of the emitted deck**, not on
``Reaction.duplicate``. The flag is an implementation detail; the deck is what
Chemkin reads and what the defect corrupts.

The two *un*-marking branches of ``mark_duplicate_reaction`` (mixed pressure
dependence, opposite direction and irreversible) are pinned here too, with
placement-mismatched pairs, because the refinement was deliberately NOT applied
to them: narrowing their condition would leave wrongly-marked duplicates marked,
which is the opposite of the repair.
"""

import os

import pytest

from rmgpy.chemkin import save_chemkin_file
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.electron_balance import get_electron_placement_counts
from rmgpy.kinetics import Arrhenius, Chebyshev
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial

#: An owner that declares the two-sided placement ``(1, 2)`` in
#: :data:`rmgpy.electron_placement.FAMILY_ELECTRON_PLACEMENT` -- one electron
#: incident, two produced. It is the only shape on this branch whose incident
#: count is not derivable from the net count, which is exactly what makes it the
#: shape a net-based comparison cannot separate from ``(0, 1)``.
DECLARED_TWO_SIDED_OWNER = "PlasmaElectronImpactIonization"

#: An owner that declares the one-sided placement ``(1, 0)``.
DECLARED_ONE_SIDED_OWNER = "PlasmaRadiativeRecombination"

#: An owner absent from the placement table, so its reactions fall back on the
#: net-derived rule: ``electrons = +1`` becomes ``(0, 1)``. A user's own kinetics
#: library is one; so is any container that has lost the declaration in transit.
UNDECLARED_OWNER = "AnUndeclaredKineticsLibrary"


def _make_species(label, index, molecule):
    """An RMG Species carrying the thermo the Chemkin writer needs."""
    spc = Species(label=label, molecule=[molecule])
    spc.index = index
    coeffs = [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, -11.7246]
    spc.thermo = NASA(
        polynomials=[
            NASAPolynomial(coeffs=coeffs, Tmin=(200, "K"), Tmax=(1000, "K")),
            NASAPolynomial(coeffs=coeffs, Tmin=(1000, "K"), Tmax=(6000, "K")),
        ],
        Tmin=(200, "K"),
        Tmax=(6000, "K"),
    )
    return spc


@pytest.fixture()
def charged_species():
    """``e-``, neutral ``[Li]`` and the cation ``[Lip]``, as one shared trio."""
    electron = _make_species("e-", 1, Molecule().from_adjacency_list("1 e u1 p0 c-1"))
    li = _make_species("Li", 2, Molecule().from_adjacency_list("multiplicity 2\n1 Li u1 p0 c0\n"))
    lip = _make_species("Lip", 3, Molecule().from_adjacency_list("1 Li u0 p0 c+1\n"))
    return electron, li, lip


@pytest.fixture()
def neutral_species():
    """A neutral trio for the control cases, which must not move at all."""
    ethane = _make_species("C2H6", 1, Molecule(smiles="CC"))
    methyl = _make_species("CH3", 2, Molecule(smiles="[CH3]"))
    ethyl = _make_species("C2H5", 3, Molecule(smiles="C[CH2]"))
    return ethane, methyl, ethyl


def _library_reaction(reactants, products, owner, electrons=0, kinetics=None,
                      reversible=False, duplicate=False):
    """
    A reaction built the way ``KineticsLibrary.load`` builds one: a
    :class:`LibraryReaction` whose ``library`` (and hence ``family``) names its
    owner, with ``electrons`` assigned after construction -- which is the line
    ``rmgpy/data/kinetics/library.py`` runs for a rate law that carries a count.
    """
    rxn = LibraryReaction(
        reactants=list(reactants),
        products=list(products),
        library=owner,
        kinetics=kinetics,
        reversible=reversible,
        duplicate=duplicate,
    )
    rxn.electrons = electrons
    return rxn


def _write_deck(tmp_path, species, reactions, name="chem.inp"):
    """Write the mechanism through the real Chemkin export path and read it back."""
    path = os.path.join(str(tmp_path), name)
    save_chemkin_file(path, species, reactions, verbose=False)
    with open(path, "r") as f:
        return f.read()


def _deck_entries(text):
    """
    Return ``[(equation, is_marked_duplicate), ...]`` parsed out of a deck's
    REACTIONS block.

    An entry starts at a line in column 0 that carries a reaction arrow; the
    auxiliary lines that follow it are indented, except ``DUPLICATE``, which
    Chemkin writes flush left on its own line and which therefore has to be
    recognised before the arrow test.
    """
    entries = []
    in_reactions = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("REACTIONS"):
            in_reactions = True
            continue
        if not in_reactions:
            continue
        if stripped == "END":
            break
        if stripped == "DUPLICATE":
            assert entries, "a DUPLICATE line appeared before any reaction"
            equation, _ = entries[-1]
            entries[-1] = (equation, True)
            continue
        if line[:1].strip() and ("=>" in line or "<=>" in line or "=" in line):
            entries.append((line.split()[0], False))
    return entries


def _stoichiometry(equation):
    """``(n_reactants, n_products)`` of a Chemkin equation string."""
    for arrow in ("<=>", "=>", "="):
        if arrow in equation:
            left, right = equation.split(arrow, 1)
            return len(left.split("+")), len(right.split("+"))
    raise AssertionError("no arrow in equation {0!r}".format(equation))


class TestElectronPlacementSeparatesDuplicates:
    """The repair: a placement mismatch must keep a pair out of DUPLICATE."""

    def test_two_owners_of_one_ionisation_channel_are_not_duplicates(self, tmp_path, charged_species):
        """
        The declared ``(1, 2)`` ionisation and an undeclared ``(0, 1)`` writing of
        the same channel must reach the deck as two independent reactions.

        They share their heavy species and their net electron count, and differ
        only in incident order -- which is precisely the difference Chemkin cares
        about, because it is the difference in stoichiometry.
        """
        electron, li, lip = charged_species

        declared = _library_reaction(
            [li], [lip], DECLARED_TWO_SIDED_OWNER, electrons=1,
            kinetics=Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        undeclared = _library_reaction(
            [li], [lip], UNDECLARED_OWNER, electrons=1,
            kinetics=Arrhenius(A=(1.0e10, "s^-1"), n=0.0, Ea=(0.0, "kcal/mol")),
        )

        # The premise the repair rests on, asserted rather than assumed: the NET
        # count is equal, so only the per-side pair can separate them.
        assert declared.electrons == undeclared.electrons == 1
        assert get_electron_placement_counts(declared) == (1, 2)
        assert get_electron_placement_counts(undeclared) == (0, 1)

        text = _write_deck(tmp_path, [electron, li, lip], [declared, undeclared])
        entries = _deck_entries(text)

        assert len(entries) == 2, "expected both reactions in the deck, got {0!r}".format(entries)
        stoichiometries = {_stoichiometry(eq) for eq, _ in entries}
        assert stoichiometries == {(2, 3), (1, 2)}, (
            "the two writings of the channel should differ in stoichiometry, got {0!r}".format(entries)
        )
        assert not any(marked for _, marked in entries), (
            "two reactions of different stoichiometry were marked DUPLICATE, which Chemkin "
            "forbids:\n{0}".format(text)
        )
        assert "DUPLICATE" not in text

    def test_matching_placement_is_still_a_duplicate(self, tmp_path, charged_species):
        """
        Two writings of the SAME channel, with the same placement, are a genuine
        duplicate pair and must still be marked.

        This is the other half of the refinement: it may only ever withdraw a
        mark, and only where the placements differ.
        """
        electron, li, lip = charged_species

        first = _library_reaction(
            [li], [lip], DECLARED_TWO_SIDED_OWNER, electrons=1,
            kinetics=Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        second = _library_reaction(
            [li], [lip], DECLARED_TWO_SIDED_OWNER, electrons=1,
            kinetics=Arrhenius(A=(3.0e11, "cm^3/(mol*s)"), n=0.5, Ea=(2.0, "kcal/mol")),
        )
        assert get_electron_placement_counts(first) == get_electron_placement_counts(second) == (1, 2)

        text = _write_deck(tmp_path, [electron, li, lip], [first, second])
        entries = _deck_entries(text)

        assert len(entries) == 2
        assert all(_stoichiometry(eq) == (2, 3) for eq, _ in entries)
        assert all(marked for _, marked in entries), (
            "a genuine duplicate pair lost its DUPLICATE lines:\n{0}".format(text)
        )
        assert text.count("DUPLICATE") == 2


class TestThreeBodyRecombinationShape:
    """
    The pair that started this ticket, in the form it takes in the real world.

    Radiative recombination and three-body recombination consume the same cation
    and make the same neutral, in the same direction, and both are irreversible
    -- so the heavy-species comparison matches and the writer marked them. Their
    stoichiometries differ by an electron on each side::

        A+ + e-        =>  A            placement (1, 0)
        A+ + e- + e-   =>  A + e-       placement (2, 1)

    The three-body channel is carried by branches that are not merged here, and
    **it cannot be stood in for** by reading this branch's ``(1, 2)`` declaration
    in the reversed orientation, which is what three-body recombination chemically
    is. The two halves of ``rmgpy.electron_balance`` do not agree about reversed
    orientation, and the disagreement is measured below: the duplicate check would
    separate the pair, but the export boundary refuses the reversed reaction
    before any deck is written. So the ``(2, 1)`` side has to come from an owner
    that declares ``(2, 1)`` forward, and this branch has none.
    """

    def test_the_duplicate_check_separates_the_pair(self, charged_species):
        """
        At the identity level the repair covers the three-body case: the two
        placements differ, so the refined comparison declines to mark them.

        This is the one assertion in this file that is not on an emitted deck,
        because on this branch there is no deck to assert on -- see the next test.
        """
        from rmgpy.chemkin import mark_duplicate_reactions

        _, li, lip = charged_species
        radiative = _library_reaction(
            [lip], [li], DECLARED_ONE_SIDED_OWNER, electrons=-1,
            kinetics=Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        three_body = _library_reaction(
            [lip], [li], DECLARED_TWO_SIDED_OWNER, electrons=-1,
            kinetics=Arrhenius(A=(1.0e18, "cm^6/(mol^2*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        assert get_electron_placement_counts(radiative) == (1, 0)
        assert get_electron_placement_counts(three_body) == (2, 1)
        assert radiative.electrons == three_body.electrons == -1, (
            "the net count is the same for both, which is why it cannot be the discriminator"
        )

        mark_duplicate_reactions([radiative, three_body])
        assert not radiative.duplicate and not three_body.duplicate

    def test_the_export_boundary_refuses_a_reversed_declared_reaction(self, tmp_path, charged_species):
        """
        ``get_electron_placement_counts`` reads a declaration in either
        orientation; ``expand_electrons`` accepts only the forward one and raises
        otherwise. That asymmetry is deliberate on the export side -- writing a
        wrong equation is unrecoverable -- but it means a three-body
        recombination cannot reach a deck by borrowing the ionisation family's
        declaration backwards. Pinned so that the limit is a measurement and not
        an assumption when the three-body branches land.
        """
        from rmgpy.exceptions import MechanismWriterError

        electron, li, lip = charged_species
        three_body = _library_reaction(
            [lip], [li], DECLARED_TWO_SIDED_OWNER, electrons=-1,
            kinetics=Arrhenius(A=(1.0e18, "cm^6/(mol^2*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        with pytest.raises(MechanismWriterError, match="declares the electron placement"):
            _write_deck(tmp_path, [electron, li, lip], [three_body])

    def test_an_undeclared_three_body_owner_is_still_folded_onto_the_radiative_channel(
            self, tmp_path, charged_species):
        """
        **The limit of this repair, pinned rather than left to be discovered.**

        The separation comes from the per-side placement, and a three-body
        recombination whose owner declares no placement does not have one: the
        net-derived rule sees ``electrons = -1`` and returns ``(1, 0)``, the same
        pair radiative recombination gets. The two then still compare equal and
        are still marked. That is not a hole in the duplicate check -- the writer
        would also *export* such a reaction as second order, losing the incident
        electron -- but it does mean the repair is only as good as the owner's
        declaration, so a three-body library must declare ``(2, 1)`` to be
        separated from the radiative channel.
        """
        electron, li, lip = charged_species

        radiative = _library_reaction(
            [lip], [li], DECLARED_ONE_SIDED_OWNER, electrons=-1,
            kinetics=Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        undeclared_three_body = _library_reaction(
            [lip], [li], UNDECLARED_OWNER, electrons=-1,
            kinetics=Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        assert get_electron_placement_counts(undeclared_three_body) == (1, 0)

        text = _write_deck(tmp_path, [electron, li, lip], [radiative, undeclared_three_body])
        entries = _deck_entries(text)

        assert len(entries) == 2
        assert {_stoichiometry(eq) for eq, _ in entries} == {(2, 1)}, (
            "both are exported second order, which is the deeper consequence of the missing "
            "declaration:\n{0}".format(text)
        )
        assert all(marked for _, marked in entries), (
            "this pair IS identical as exported, so it is correctly still a duplicate; if this "
            "ever stops holding, the rule has started guessing:\n{0}".format(text)
        )


class TestOrdinaryDuplicateHandlingIsUntouched:
    """Neutral chemistry, and both un-marking branches, must not move."""

    def test_neutral_duplicate_pair_is_still_marked(self, tmp_path, neutral_species):
        """
        The commonest case in RMG: two neutral reactions over the same species.
        Both placements are ``(0, 0)``, so the refinement is a no-op and the pair
        must still be marked.
        """
        ethane, methyl, ethyl = neutral_species
        kwargs = dict(reversible=True)
        first = _library_reaction(
            [ethane], [methyl, methyl], "SomeNeutralLibrary",
            kinetics=Arrhenius(A=(1.0e16, "s^-1"), n=0.0, Ea=(80.0, "kcal/mol")), **kwargs
        )
        second = _library_reaction(
            [ethane], [methyl, methyl], "SomeNeutralLibrary",
            kinetics=Arrhenius(A=(5.0e15, "s^-1"), n=0.2, Ea=(75.0, "kcal/mol")), **kwargs
        )
        assert get_electron_placement_counts(first) == get_electron_placement_counts(second) == (0, 0)

        text = _write_deck(tmp_path, [ethane, methyl, ethyl], [first, second])
        entries = _deck_entries(text)

        assert len(entries) == 2
        assert all(marked for _, marked in entries), (
            "a neutral duplicate pair lost its DUPLICATE lines:\n{0}".format(text)
        )
        assert text.count("DUPLICATE") == 2

    def test_opposite_direction_irreversible_branch_still_unmarks(self, tmp_path, charged_species):
        """
        The opposite-direction irreversible branch must keep unmarking a pair that
        arrives already flagged -- INCLUDING one whose placements differ.

        This is why the refinement was applied only to the branch that marks. Had
        the match flags themselves been narrowed, this pair would have stopped
        matching at all, the branch would never have run, and two reactions that
        are not duplicates would have kept their DUPLICATE lines: a mis-marked
        deck produced by the very change meant to prevent one.
        """
        electron, li, lip = charged_species

        ionisation = _library_reaction(
            [li], [lip], DECLARED_TWO_SIDED_OWNER, electrons=1, duplicate=True,
            kinetics=Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        recombination = _library_reaction(
            [lip], [li], DECLARED_ONE_SIDED_OWNER, electrons=-1, duplicate=True,
            kinetics=Arrhenius(A=(2.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        assert get_electron_placement_counts(ionisation) == (1, 2)
        assert get_electron_placement_counts(recombination) == (1, 0)

        text = _write_deck(tmp_path, [electron, li, lip], [ionisation, recombination])
        entries = _deck_entries(text)

        assert len(entries) == 2
        assert not any(marked for _, marked in entries), (
            "the opposite-direction irreversible branch stopped unmarking:\n{0}".format(text)
        )
        assert "DUPLICATE" not in text

    def test_mixed_pressure_dependence_branch_still_unmarks(self, tmp_path, charged_species):
        """
        The mixed-pressure-dependence branch must keep unmarking an already-flagged
        pair, with mismatched placements as well. Same reasoning as above.
        """
        electron, li, lip = charged_species

        declared = _library_reaction(
            [li], [lip], DECLARED_TWO_SIDED_OWNER, electrons=1, duplicate=True,
            kinetics=Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        undeclared = _library_reaction(
            [li], [lip], UNDECLARED_OWNER, electrons=1, duplicate=True,
            kinetics=Chebyshev(
                coeffs=[[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                kunits="s^-1",
                Tmin=(300, "K"), Tmax=(2000, "K"),
                Pmin=(0.01, "bar"), Pmax=(100, "bar"),
            ),
        )
        assert declared.kinetics.is_pressure_dependent() != undeclared.kinetics.is_pressure_dependent()

        text = _write_deck(tmp_path, [electron, li, lip], [declared, undeclared])
        entries = _deck_entries(text)

        assert len(entries) == 2
        assert not any(marked for _, marked in entries), (
            "the mixed-pressure-dependence branch stopped unmarking:\n{0}".format(text)
        )
        assert "DUPLICATE" not in text


class TestPlacementTruthTable:
    """
    The rule, enumerated: a pair may be marked only when its per-side electron
    placements agree in the orientation its heavy species agreed in.

    Every row is driven through the real writer and read off the emitted deck.
    The neutral row is the one that bounds the blast radius -- ``(0, 0)`` against
    ``(0, 0)`` for all of RMG outside the charged families and plasma libraries.
    """

    @pytest.mark.parametrize(
        "electrons_a, electrons_b, owner_a, owner_b, expect_marked, why",
        [
            (0, 0, "NeutralLibA", "NeutralLibB", True,
             "(0, 0) vs (0, 0) -- all neutral chemistry, unchanged"),
            (1, 1, DECLARED_TWO_SIDED_OWNER, DECLARED_TWO_SIDED_OWNER, True,
             "(1, 2) vs (1, 2) -- a genuine duplicate of a declared channel"),
            (1, 1, UNDECLARED_OWNER, "AnotherUndeclaredLibrary", True,
             "(0, 1) vs (0, 1) -- a genuine duplicate under the net-derived rule"),
            (1, 1, DECLARED_TWO_SIDED_OWNER, UNDECLARED_OWNER, False,
             "(1, 2) vs (0, 1) -- equal net count, different incident order"),
            (1, 1, UNDECLARED_OWNER, DECLARED_TWO_SIDED_OWNER, False,
             "(0, 1) vs (1, 2) -- the same row with the arguments swapped"),
        ],
    )
    def test_same_direction_rows(self, tmp_path, charged_species, neutral_species,
                                 electrons_a, electrons_b, owner_a, owner_b,
                                 expect_marked, why):
        electron, li, lip = charged_species
        ethane, methyl, _ = neutral_species

        if electrons_a == 0:
            species = [ethane, methyl]
            reactants, products = [ethane], [methyl, methyl]
            units_a = units_b = "s^-1"
        else:
            species = [electron, li, lip]
            reactants, products = [li], [lip]
            # The exported reactant side is one species longer wherever the owner
            # declares an incident electron, and the A-factor units have to match
            # the equation the writer emits.
            units_a = "cm^3/(mol*s)" if owner_a == DECLARED_TWO_SIDED_OWNER else "s^-1"
            units_b = "cm^3/(mol*s)" if owner_b == DECLARED_TWO_SIDED_OWNER else "s^-1"

        first = _library_reaction(
            reactants, products, owner_a, electrons=electrons_a,
            kinetics=Arrhenius(A=(1.0e12, units_a), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        second = _library_reaction(
            reactants, products, owner_b, electrons=electrons_b,
            kinetics=Arrhenius(A=(2.0e12, units_b), n=0.0, Ea=(0.0, "kcal/mol")),
        )

        text = _write_deck(tmp_path, species, [first, second])
        entries = _deck_entries(text)

        assert len(entries) == 2, text
        marked = [flag for _, flag in entries]
        assert marked == [expect_marked, expect_marked], (
            "{0}: expected marked={1}, deck says {2}\n{3}".format(why, expect_marked, marked, text)
        )
