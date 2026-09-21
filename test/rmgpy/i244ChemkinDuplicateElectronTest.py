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

**Assertions about a deck are made on the text of the emitted deck**, not on
``Reaction.duplicate``. The flag is an implementation detail; the deck is what
Chemkin reads and what the defect corrupts.

``duplicate`` is nevertheless a *group* predicate -- Chemkin requires every entry
whose equation another entry also writes to carry ``DUPLICATE``, and an entry
carrying it alone is an error too -- while ``mark_duplicate_reaction`` answers a
*pairwise* question. That mismatch is not closable by refining the pairwise
predicate, so ``chemkin_duplicate_flags`` recomputes every flag from
``chemkin_duplicate_group_key`` over a whole list, and each Chemkin renderer calls
it for the list it is about to serialize. ``TestDuplicateIsAGroupPredicate`` holds
that line. Not *every* writer consults it -- :mod:`arkane.pdep` and
:mod:`rmgpy.yaml_cantera2` do not -- which is why the answer is passed to the
writer as an argument and never stored back on the reactions;
``TestTheAnswerBelongsToTheDeck`` holds that line.

The two *un*-marking branches of the pairwise ``mark_duplicate_reaction`` (mixed
pressure dependence, opposite direction and irreversible) survive for its
incremental caller in :mod:`rmgpy.rmg.model`, and are pinned here by driving that
function directly and asserting the branch-specific warning. A deck-level
assertion cannot pin them any more: the group recompute would produce the same
deck with either branch deleted.
"""

import logging
import os
from types import SimpleNamespace

import pytest

from rmgpy.chemkin import (mark_duplicate_reaction, mark_duplicate_reactions, save_chemkin,
                           save_chemkin_file)
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


@pytest.fixture()
def permutation_species():
    """
    ``H``, ``OH``, ``H2`` and ``O``: a neutral quartet whose reaction
    ``H + OH => H2 + O`` is element-balanced and has two DISTINCT species on each
    side, which is what a permutation case needs. Cantera checks the balance
    before it checks the duplicates, so an unbalanced stand-in is rejected for a
    reason that has nothing to do with this ticket.
    """
    return (_make_species("H", 1, Molecule(smiles="[H]")),
            _make_species("OH", 2, Molecule(smiles="[OH]")),
            _make_species("H2", 3, Molecule(smiles="[H][H]")),
            _make_species("O", 4, Molecule(smiles="[O]")))


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

    def test_opposite_direction_irreversible_branch_still_unmarks(self, caplog, charged_species):
        """
        The opposite-direction irreversible branch of the *pairwise*
        ``mark_duplicate_reaction`` must keep unmarking a pair that arrives already
        flagged.

        This drives that function directly rather than writing a deck, because the
        deck no longer exercises it: ``mark_duplicate_reactions`` recomputes every
        flag from the group key, and would clear this pair whether the branch fired
        or not. An outcome-only assertion through the writer would therefore pass
        with the branch deleted. The distinctive warning the branch logs is what
        pins the branch itself -- no other code emits that sentence.
        """
        electron, li, lip = charged_species
        del electron

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

        with caplog.at_level(logging.WARNING):
            mark_duplicate_reaction(ionisation, [recombination])

        assert any("irreversible in opposite directions" in record.getMessage()
                   for record in caplog.records), (
            "the opposite-direction irreversible branch did not run; messages were {0!r}".format(
                [record.getMessage() for record in caplog.records])
        )
        assert (ionisation.duplicate, recombination.duplicate) == (False, False)

    def test_mixed_pressure_dependence_branch_still_unmarks(self, caplog, charged_species):
        """
        The mixed-pressure-dependence branch of the pairwise
        ``mark_duplicate_reaction`` must keep unmarking an already-flagged pair.
        Driven and pinned the same way, and for the same reason, as the test above.
        """
        electron, li, lip = charged_species
        del electron

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

        with caplog.at_level(logging.WARNING):
            mark_duplicate_reaction(declared, [undeclared])

        assert any("mixed pressure dependence" in record.getMessage()
                   for record in caplog.records), (
            "the mixed-pressure-dependence branch did not run; messages were {0!r}".format(
                [record.getMessage() for record in caplog.records])
        )
        assert (declared.duplicate, undeclared.duplicate) == (False, False)


class TestDuplicateIsAGroupPredicate:
    """
    ``Reaction.duplicate`` records membership in a duplicate *group*, and Chemkin's
    rule is a group rule: every entry whose equation another entry also writes must
    carry ``DUPLICATE``, and an entry carrying it alone is an error in its own
    right. The pairwise ``mark_duplicate_reaction`` decides that group question one
    pair at a time, off a boolean that cannot name *which* group it means, so it
    cannot answer it -- no refinement of the pairwise predicate can.

    These are the cases that mismatch produces. Every one of them puts
    ``duplicate = True`` on at least one input before calling, because that is the
    state the defect lives in: a pair that arrives unmarked was already handled.
    """

    def test_a_pre_marked_placement_mismatched_pair_is_cleared(self, tmp_path, charged_species):
        """
        The pair this ticket exists to separate -- ``(1, 2)`` against ``(0, 1)`` --
        arriving with both flags already set.

        The electron-placement refinement lives in the branch that marks, which is
        the ``else`` of ``if reaction1.duplicate and reaction2.duplicate``. A
        pre-marked pair never reaches it, so before the group-level recompute these
        two shipped as duplicates of each other despite writing equations of
        different stoichiometry.
        """
        electron, li, lip = charged_species

        declared = _library_reaction(
            [li], [lip], DECLARED_TWO_SIDED_OWNER, electrons=1, duplicate=True,
            kinetics=Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        undeclared = _library_reaction(
            [li], [lip], UNDECLARED_OWNER, electrons=1, duplicate=True,
            kinetics=Arrhenius(A=(2.0e12, "s^-1"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        assert (declared.duplicate, undeclared.duplicate) == (True, True)
        assert get_electron_placement_counts(declared) == (1, 2)
        assert get_electron_placement_counts(undeclared) == (0, 1)

        text = _write_deck(tmp_path, [electron, li, lip], [declared, undeclared])
        entries = _deck_entries(text)

        assert len(entries) == 2
        assert {_stoichiometry(equation) for equation, _ in entries} == {(2, 3), (1, 2)}, (
            "the two channels should write different stoichiometries:\n{0}".format(text)
        )
        assert not any(marked for _, marked in entries), (
            "a pre-marked pair with different electron placements kept its DUPLICATE "
            "lines:\n{0}".format(text)
        )

    def test_a_cross_group_comparison_does_not_clear_another_groups_flags(self, tmp_path,
                                                                          charged_species):
        """
        Four reactions in two groups of two: a ``(1, 2)`` pair that is not pressure
        dependent, and a ``(0, 1)`` pair that is.

        Comparing a member of one group against a member of the other reaches the
        mixed-pressure-dependence un-marking branch, which clears *both* flags --
        and one of them belonged to the group that was never in question. The deck
        then carried two identical equations with no ``DUPLICATE`` line at all,
        which Chemkin rejects for the opposite reason to the case above.

        **One of the four arrives unflagged on purpose.** Written with all four
        pre-marked and all four expected marked, this test passed against an
        authority that did nothing at all -- every input already held the value
        every assertion wanted, so it could only catch an authority that *cleared*,
        never one that was absent. ``undeclared_b`` therefore starts ``False`` and
        has to be *set* by the recompute, which is a state no no-op reaches. The
        defect under test is unchanged by that: the cross-group contamination runs
        between ``declared_a``/``declared_b`` and ``undeclared_a``, all three of
        which still arrive flagged.
        """
        electron, li, lip = charged_species

        def chebyshev():
            return Chebyshev(
                coeffs=[[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                kunits="s^-1",
                Tmin=(300, "K"), Tmax=(2000, "K"),
                Pmin=(0.01, "bar"), Pmax=(100, "bar"),
            )

        declared_a = _library_reaction(
            [li], [lip], DECLARED_TWO_SIDED_OWNER, electrons=1, duplicate=True,
            kinetics=Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        declared_b = _library_reaction(
            [li], [lip], DECLARED_TWO_SIDED_OWNER, electrons=1, duplicate=True,
            kinetics=Arrhenius(A=(2.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        undeclared_a = _library_reaction(
            [li], [lip], UNDECLARED_OWNER, electrons=1, duplicate=True, kinetics=chebyshev(),
        )
        undeclared_b = _library_reaction(
            [li], [lip], "AnotherUndeclaredKineticsLibrary", electrons=1, duplicate=False,
            kinetics=chebyshev(),
        )
        reactions = [declared_a, declared_b, undeclared_a, undeclared_b]
        assert [rxn.duplicate for rxn in reactions] == [True, True, True, False], (
            "the incoming state must not already be the expected one, or a no-op passes"
        )

        text = _write_deck(tmp_path, [electron, li, lip], reactions)
        entries = _deck_entries(text)

        assert len(entries) == 4
        by_equation = {}
        for equation, marked in entries:
            by_equation.setdefault(equation, []).append(marked)
        assert len(by_equation) == 2, "expected two distinct equations:\n{0}".format(text)
        for equation, marks in by_equation.items():
            assert marks == [True, True], (
                "{0!r} is written twice and must carry DUPLICATE on both entries, "
                "got {1}:\n{2}".format(equation, marks, text)
            )

    def test_a_lone_flag_that_names_no_mate_is_cleared(self, tmp_path, charged_species):
        """
        A single reaction arriving flagged, with nothing in the deck that writes its
        equation. A lone ``DUPLICATE`` line is a Chemkin error, and the pairwise
        function can never clear it because it never sees the whole deck.
        """
        electron, li, lip = charged_species

        lonely = _library_reaction(
            [li], [lip], DECLARED_TWO_SIDED_OWNER, electrons=1, duplicate=True,
            kinetics=Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        assert lonely.duplicate is True

        text = _write_deck(tmp_path, [electron, li, lip], [lonely])
        entries = _deck_entries(text)

        assert len(entries) == 1
        assert entries[0][1] is False, "a lone entry kept its DUPLICATE line:\n{0}".format(text)

    def test_a_pre_marked_mixed_reversibility_pair_is_cleared(self, tmp_path, neutral_species):
        """
        The same defect class with no electrons in it at all: ``A<=>B`` and ``A=>B``
        write different equations, so they are not duplicates -- and neither
        un-marking branch covers mixed reversibility, so an already-flagged pair of
        this shape kept two lone ``DUPLICATE`` lines.

        Included because it shows the repair is about the group predicate rather
        than about electrons; the electron placement is simply one more thing the
        group key has to carry.
        """
        ethane, methyl, ethyl = neutral_species

        reversible = _library_reaction(
            [ethane], [methyl, methyl], "SomeNeutralLibrary", duplicate=True, reversible=True,
            kinetics=Arrhenius(A=(1.0e16, "s^-1"), n=0.0, Ea=(80.0, "kcal/mol")),
        )
        irreversible = _library_reaction(
            [ethane], [methyl, methyl], "SomeNeutralLibrary", duplicate=True, reversible=False,
            kinetics=Arrhenius(A=(5.0e15, "s^-1"), n=0.2, Ea=(75.0, "kcal/mol")),
        )
        assert (reversible.duplicate, irreversible.duplicate) == (True, True)

        text = _write_deck(tmp_path, [ethane, methyl, ethyl], [reversible, irreversible])
        entries = _deck_entries(text)

        assert len(entries) == 2
        assert len({equation for equation, _ in entries}) == 2, (
            "the two entries should write different equations:\n{0}".format(text)
        )
        assert not any(marked for _, marked in entries), (
            "a mixed-reversibility pair kept its DUPLICATE lines:\n{0}".format(text)
        )

    def test_the_verdict_does_not_depend_on_the_flags_it_arrives_with(self, charged_species):
        """
        The recompute reads the reactions, not their flags. The same list must come
        out the same way whether every flag arrives set or every flag arrives clear,
        and a second pass must change nothing.

        The pairwise sweep could not have this property: arriving all-set, it had no
        branch that would clear the lone reaction below, and arriving all-clear it
        had none that would set it. The two runs disagreed.
        """
        electron, li, lip = charged_species
        del electron

        def build(flag):
            mate_a = _library_reaction(
                [li], [lip], DECLARED_TWO_SIDED_OWNER, electrons=1, duplicate=flag,
                kinetics=Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
            )
            mate_b = _library_reaction(
                [li], [lip], DECLARED_TWO_SIDED_OWNER, electrons=1, duplicate=flag,
                kinetics=Arrhenius(A=(2.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
            )
            lonely = _library_reaction(
                [li], [lip], UNDECLARED_OWNER, electrons=1, duplicate=flag,
                kinetics=Arrhenius(A=(3.0e12, "s^-1"), n=0.0, Ea=(0.0, "kcal/mol")),
            )
            return [mate_a, mate_b, lonely]

        from_set = build(True)
        mark_duplicate_reactions(from_set)
        from_clear = build(False)
        mark_duplicate_reactions(from_clear)

        expected = [True, True, False]
        assert [rxn.duplicate for rxn in from_set] == expected
        assert [rxn.duplicate for rxn in from_clear] == expected

        mark_duplicate_reactions(from_set)
        assert [rxn.duplicate for rxn in from_set] == expected, "the recompute is not idempotent"

    def test_the_production_save_path_recomputes_the_flags(self, tmp_path, charged_species):
        """
        ``save_chemkin`` is the function an RMG run writes its deck with. The repair
        reaches a real run's deck only if the recompute happens somewhere on that
        path; without it, the flags a run ships are whatever ``rmgpy.rmg.model``'s
        incremental pairwise calls left behind, and this pair ships marked.

        The recompute is per render call, over that render's own reaction list --
        ``save_chemkin`` itself no longer does one. ``TestTheAnswerBelongsToTheDeck``
        covers why.
        """
        electron, li, lip = charged_species

        declared = _library_reaction(
            [li], [lip], DECLARED_TWO_SIDED_OWNER, electrons=1, duplicate=True,
            kinetics=Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        undeclared = _library_reaction(
            [li], [lip], UNDECLARED_OWNER, electrons=1, duplicate=True,
            kinetics=Arrhenius(A=(2.0e12, "s^-1"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        model = SimpleNamespace(
            core=SimpleNamespace(species=[electron, li, lip], reactions=[declared, undeclared]),
            edge=SimpleNamespace(species=[], reactions=[]),
            output_species_list=[],
            output_reaction_list=[],
            surface_site_density=None,
        )

        path = os.path.join(str(tmp_path), "chem_annotated.inp")
        verbose_path = os.path.join(str(tmp_path), "chem_annotated_verbose.inp")
        save_chemkin(model, path, verbose_path)

        with open(path, "r") as f:
            text = f.read()
        entries = _deck_entries(text)

        assert len(entries) == 2
        assert not any(marked for _, marked in entries), (
            "save_chemkin shipped a placement-mismatched pair as duplicates:\n{0}".format(text)
        )


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


class TestTheAnswerBelongsToTheDeck:
    """
    ``DUPLICATE`` is a statement about one deck, not about one reaction.

    The group recompute gained the power to *clear* flags, which is what makes an
    RMG edge deck loadable at all -- but the same power, pointed at the wrong key
    or left lying on the wrong object, breaks decks that used to be fine. The
    three cases below are the three ways that happened, each of them measured on
    both engines before it was repaired (see
    ``docs/i244-chemkin-duplicate-electron-aware/round73_repro.py``).

    **Every duplicate claim here is asserted by loading the mechanism**, with
    ``cantera.Solution``. Neither cheaper check sees these defects: the two
    entries' rendered text genuinely differs (``CH3(2)+C2H5(3)=>...`` against
    ``C2H5(3)+CH3(2)=>...``), so a text comparison finds nothing, and
    ``ck2yaml.convert_mech`` writes the YAML without ever running
    ``Kinetics::checkDuplicates``, so it reports success on an invalid deck.
    """

    @staticmethod
    def _load_with_cantera(tmp_path, text, name="chem.inp"):
        """
        Convert `text` and then LOAD it. Returns ``None`` when Cantera accepts the
        mechanism, or a one-line reason when it does not.
        """
        import cantera as ct
        import cantera.ck2yaml as ck2yaml

        inp = os.path.join(str(tmp_path), name)
        out = os.path.join(str(tmp_path), name + ".yaml")
        with open(inp, "w") as f:
            f.write(text)
        try:
            ck2yaml.convert_mech(input_file=inp, out_name=out, quiet=True)
        except Exception as exc:
            return "ck2yaml.convert_mech rejected it: {0}".format(
                str(exc).strip().splitlines()[0])
        try:
            ct.Solution(out)
        except Exception as exc:
            # The whole message, not a line grepped out of it. Cantera echoes the
            # offending YAML into the error, so a filter for "duplicate" reliably
            # matches the echo of `duplicate: true` and reports a duplicate
            # complaint no matter what the actual objection was -- which is how an
            # unbalanced test reaction first read as a duplicate failure here.
            return "cantera.Solution rejected it:\n{0}".format(exc)
        return None

    @staticmethod
    def _permuted_pair(species):
        """
        Two entries with the same participants on each side, whose reactant LISTS
        differ only in order, both arriving already flagged. Chemkin and Cantera
        normalise the reactant multiset, so these are one equation written twice
        and both entries must keep ``DUPLICATE``.
        """
        h, oh, h2, o = species
        first = _library_reaction(
            [h, oh], [h2, o], "NeutralLibA", duplicate=True,
            kinetics=Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        second = _library_reaction(
            [oh, h], [h2, o], "NeutralLibB", duplicate=True,
            kinetics=Arrhenius(A=(2.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        return first, second

    def test_a_permuted_pair_is_not_cleared_into_a_deck_cantera_rejects(self, tmp_path,
                                                                        permutation_species):
        """
        The key keys each side as a multiset. Under an ORDERED key these two land
        in different groups, each a singleton, both flags are cleared, and the deck
        Cantera gets has one equation written twice with no ``DUPLICATE`` line --
        which ``Kinetics::checkDuplicates`` refuses.
        """
        first, second = self._permuted_pair(permutation_species)
        text = _write_deck(tmp_path, list(permutation_species), [first, second])
        entries = _deck_entries(text)

        assert len(entries) == 2, text
        assert [marked for _, marked in entries] == [True, True], (
            "a pair differing only in reactant order lost its DUPLICATE lines:\n{0}".format(text)
        )
        equations = [equation for equation, _ in entries]
        assert equations[0] != equations[1], (
            "this case only bites while the two entries render differently; if they "
            "render alike the test has stopped covering it:\n{0}".format(text)
        )
        rejection = self._load_with_cantera(tmp_path, text)
        assert rejection is None, rejection

    def test_the_lone_flag_case_still_clears(self, tmp_path, permutation_species):
        """
        The guard against over-correcting the case above: widening the key must not
        cost the clearing. A single flagged entry with no mate still loses its line,
        and Cantera still accepts the result.
        """
        h, oh, h2, o = permutation_species
        lonely = _library_reaction(
            [h, oh], [h2, o], "NeutralLibA", duplicate=True,
            kinetics=Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        text = _write_deck(tmp_path, list(permutation_species), [lonely])
        entries = _deck_entries(text)

        assert len(entries) == 1, text
        assert entries[0][1] is False, "a lone entry kept its DUPLICATE line:\n{0}".format(text)
        rejection = self._load_with_cantera(tmp_path, text)
        assert rejection is None, rejection

    def test_a_render_leaves_every_reactions_flag_exactly_as_it_found_it(self, tmp_path,
                                                                        permutation_species):
        """
        The answer is per list, so it is passed to the writer and not stored. A
        render that wrote its answer back would hand it to whichever writer ran
        next over a different list.
        """
        h, oh, h2, o = permutation_species
        first, second = self._permuted_pair(permutation_species)
        lonely = _library_reaction(
            [h2], [h, h], "NeutralLibC", duplicate=True,
            kinetics=Arrhenius(A=(3.0e12, "s^-1"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        before = [first.duplicate, second.duplicate, lonely.duplicate]
        _write_deck(tmp_path, list(permutation_species), [first, second, lonely])
        after = [first.duplicate, second.duplicate, lonely.duplicate]

        assert after == before == [True, True, True], (
            "rendering a deck rewrote the reactions' flags: {0} -> {1}".format(before, after)
        )

    def test_a_core_plus_edge_save_does_not_mark_a_core_only_cantera_entry(self, tmp_path,
                                                                          permutation_species):
        """
        An RMG run saves the core deck, then the core+edge deck, over the same
        reaction objects, and the Cantera YAML writer runs afterwards on the core
        alone. A core reaction whose only mate lives on the edge is a duplicate in
        the second deck and not in the first. If the second save stores its answer,
        the Cantera writer reads it and emits a lone ``duplicate: true``, which
        Cantera rejects for the same reason a lone ``DUPLICATE`` line is a Chemkin
        error.
        """
        from rmgpy.yaml_cantera2 import reaction_to_dict_list

        h, oh, h2, o = permutation_species
        core = _library_reaction(
            [h, oh], [h2, o], "NeutralLibA",
            kinetics=Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        edge = _library_reaction(
            [h, oh], [h2, o], "NeutralLibB",
            kinetics=Arrhenius(A=(2.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
        )
        model = SimpleNamespace(
            core=SimpleNamespace(species=list(permutation_species), reactions=[core]),
            edge=SimpleNamespace(species=[], reactions=[edge]),
            output_species_list=[],
            output_reaction_list=[],
            surface_site_density=None,
        )

        core_path = os.path.join(str(tmp_path), "chem.inp")
        edge_path = os.path.join(str(tmp_path), "chem_edge.inp")
        save_chemkin(model, core_path, os.path.join(str(tmp_path), "chem_annotated.inp"),
                     save_edge_species=False)
        save_chemkin(model, edge_path, os.path.join(str(tmp_path), "chem_edge_annotated.inp"),
                     save_edge_species=True)

        with open(core_path) as f:
            core_text = f.read()
        with open(edge_path) as f:
            edge_text = f.read()

        assert [marked for _, marked in _deck_entries(core_text)] == [False], (
            "the core deck's lone entry was marked:\n{0}".format(core_text))
        assert [marked for _, marked in _deck_entries(edge_text)] == [True, True], (
            "the core+edge deck's genuine pair was not marked:\n{0}".format(edge_text))

        entry = reaction_to_dict_list(core, list(permutation_species))[0]
        assert not entry.get("duplicate", False), (
            "the core+edge save leaked its answer to the Cantera writer, which emitted "
            "a lone duplicate:true for {0}".format(entry["equation"])
        )

    def test_a_generator_of_reactions_is_written_not_consumed(self, tmp_path,
                                                              permutation_species):
        """
        ``mark_duplicate_reaction`` has promised for years that its reaction list
        "can be any iterator". Keying the groups reads the argument once; writing
        the entries reads it again. A generator passed in was therefore exhausted
        by the first read and the writer emitted an empty mechanism -- and returned
        success, so the loss was silent.
        """
        h, oh, h2, o = permutation_species

        def one():
            return _library_reaction(
                [h, oh], [h2, o], "NeutralLibA",
                kinetics=Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol")),
            )

        from_list = _write_deck(tmp_path, list(permutation_species), [one()], name="list.inp")
        from_generator = _write_deck(tmp_path, list(permutation_species),
                                     (rxn for rxn in [one()]), name="gen.inp")

        assert len(_deck_entries(from_list)) == 1, from_list
        assert len(_deck_entries(from_generator)) == len(_deck_entries(from_list)), (
            "a generator argument lost the mechanism:\n{0}".format(from_generator)
        )
