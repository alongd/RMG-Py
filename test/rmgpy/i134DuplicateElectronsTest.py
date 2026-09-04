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
I-134 -- reaction identity includes the electron, per side.

The defect these tests pin: RMG decided whether a proposed reaction was already
in the model by comparing reactant and product *references* and *labels*, and
nothing else. The canonical representation keeps the electron out of the
participant lists and in the scalar ``Reaction.electrons``, so two reactions over
the same heavy species were indistinguishable to that comparison. The two shipped
lithium plasma channels are exactly that shape --

    electron-impact ionisation   Li  + e-  =>  Li+ + 2 e-    placement (1, 2)
    radiative recombination      Li+ + e-  =>  Li  + hv      placement (1, 0)

-- so the model kept whichever was offered first and discarded the other as a
duplicate, silently, at ``INFO`` level, with the run exiting 0. A deck asking for
both wrote a mechanism with a cation source and no cation loss channel.

**Why the net scalar is not the fix, which is the subtle part.** The two carry
``electrons = +1`` and ``-1``: exactly equal and opposite, which is the relation a
genuine reverse pair has. A comparison built on the net count therefore also calls
them the same reaction -- and ``Reaction.is_isomorphic`` did, measured, before this
change. They are not reverses: reversing the ionisation gives three-body
recombination, ``(2, 1)``, third order, not the radiative channel's ``(1, 0)``,
second order. Only the per-side placement separates them, which is why the repair
runs through :func:`rmgpy.electron_balance.get_electron_placement_counts` and not
through ``reaction.electrons``.

**The blast radius, which is why the repair is in the identity predicates and not
in the plasma code.** The blind comparison is shared by all four return sites of
``check_for_existing_reaction`` -- two of them the family branches -- and by
``Reaction.is_isomorphic``, which ``ReactionModel.merge``,
``KineticsLibrary.check_for_duplicates`` and the pressure-dependence code all
consume. ``TestChargedFamilyReactionsAreNotConfused`` is the proof that this was
never a lithium special case: it is a real family-generated pair, out of the
database, with no lithium in it.

The classes below pin, in order: the contract for the overwhelming majority of
reactions, which have no electron placement declaration and whose verdicts must
not move at all; that genuine duplicates still collapse, which is what the check
is for; that a real charged family's reactions are no longer confused with each
other; that the real lithium mechanism reaches the model with both channels; and
that a kinetics library carrying both channels can be loaded, which it could not
be before.
"""

import os

import pytest

import rmgpy.data.rmg
from rmgpy import settings
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import KineticsLibrary, LibraryReaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.electron_balance import get_electron_placement_counts, get_placement_declaration
from rmgpy.electron_placement import FAMILY_ELECTRON_PLACEMENT
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import CoreEdgeReactionModel, are_identical_species_references
from rmgpy.species import Species

IONISATION = 'PlasmaElectronImpactIonization'
RECOMBINATION = 'PlasmaRadiativeRecombination'
ATTACHMENT = 'Plasma_Electron_Attachment'

#: Triplet dioxygen, the reactant the attachment family is trained on.
O2_ADJACENCY_LIST = """
multiplicity 3
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""


def _net_rule_counts(electrons):
    """The placement the net-derived rule produces: the rule
    :func:`rmgpy.electron_balance.expand_electrons` falls back to, restated here
    independently so the tests below compare against a second statement of it
    rather than against the implementation under test."""
    if electrons < 0:
        return -electrons, 0
    return 0, electrons


def _species(smiles):
    return Species(molecule=[Molecule().from_smiles(smiles)])


class _GlobalKineticsDatabase:
    """Register a kinetics database as the module-level singleton, and put back
    whatever was there before.

    ``check_for_existing_reaction`` resolves a reaction's owner through
    ``get_family_library_object``, which reads ``rmgpy.data.rmg.database`` and
    not the model's own ``kinetics_database`` attribute -- so a test driving the
    real model path has to register one. Constructing an :class:`RMGDatabase`
    does that as a side effect of ``__init__``, which means it also clobbers the
    singleton for every test that runs afterwards in the same session; hence the
    restore.
    """

    def __init__(self, kinetics_database):
        self.kinetics_database = kinetics_database
        self.saved = None

    def __enter__(self):
        self.saved = rmgpy.data.rmg.database
        database = RMGDatabase()
        database.kinetics = self.kinetics_database
        rmgpy.data.rmg.database = database
        return database

    def __exit__(self, *exc_info):
        rmgpy.data.rmg.database = self.saved
        return False


def _reaction(reactants, products, electrons=0, family=None):
    """A reaction over shared Species objects.

    ``family`` left at ``None`` gives a plain :class:`Reaction`, which is an
    undeclared owner -- and, being a ``cdef`` class, cannot be given a ``family``
    attribute at all, so an owner has to be expressed by using a subclass that has
    one. ``LibraryReaction`` sets ``family = library``, which is how a kinetics
    library declares its electron placement on the same terms a family does.
    """
    if family is None:
        return Reaction(reactants=list(reactants), products=list(products),
                        electrons=electrons)
    return LibraryReaction(reactants=list(reactants), products=list(products),
                           electrons=electrons, library=family, reversible=False)


class TestTheUndeclaredOwnerContract:
    """The contract for every reaction whose owner has never heard of electrons,
    which is all of RMG outside the plasma families and libraries.

    This is stated as a test rather than left as an inference from the suite
    passing, because it is what bounds what the repair can move. The claim is not
    "nothing broke", and it is deliberately different for the two predicates:

    * ``Reaction.is_isomorphic`` already compared the net scalar, and the
      per-side counts *reduce* to exactly that comparison for an undeclared
      owner, so its verdicts are unchanged by construction over the whole space
      of electron counts;
    * ``are_identical_species_references`` compared no electron information at
      all, so for a CHARGED reaction its verdict genuinely changes -- that is the
      repair. What bounds it is that a neutral reaction's counts are ``(0, 0)``,
      so neutral chemistry, which is all of RMG outside the charged families and
      the plasma libraries, is untouched.
    """

    ELECTRON_COUNTS = [-3, -2, -1, 0, 1, 2, 3]

    @pytest.mark.parametrize('electrons', ELECTRON_COUNTS)
    def test_counts_are_the_net_rule(self, electrons):
        rxn = _reaction([_species('[OH]')], [_species('O')], electrons=electrons)
        assert get_electron_placement_counts(rxn) == _net_rule_counts(electrons)

    def test_a_reaction_with_no_family_attribute_at_all(self):
        """A plain :class:`Reaction` -- the reactor's placement view is one -- has
        no ``family`` attribute, so the declaration lookup must not raise."""
        rxn = Reaction(reactants=[_species('[OH]')], products=[_species('O')], electrons=-1)
        assert not hasattr(rxn, 'family') or rxn.family is None
        assert get_electron_placement_counts(rxn) == (1, 0)

    def test_an_unrecognised_owner_is_undeclared(self):
        rxn = _reaction([_species('[OH]')], [_species('O')], electrons=2,
                        family='H_Abstraction')
        assert get_placement_declaration(rxn) is None
        assert get_electron_placement_counts(rxn) == (0, 2)

    @pytest.mark.parametrize('electrons1', ELECTRON_COUNTS)
    @pytest.mark.parametrize('electrons2', ELECTRON_COUNTS)
    def test_is_isomorphic_forward_verdict_is_exactly_its_old_net_comparison(
            self, electrons1, electrons2):
        """``is_isomorphic`` already compared ``self.electrons == other.electrons``
        forward. Over the whole space of pairs, the per-side comparison agrees with
        it exactly, which is what makes the change to that predicate a
        restatement for undeclared owners rather than a tightening."""
        a, b = _species('[OH]'), _species('O')
        rxn1 = _reaction([a], [b], electrons=electrons1)
        rxn2 = _reaction([a], [b], electrons=electrons2)
        assert rxn1.is_isomorphic(rxn2, either_direction=False) is (electrons1 == electrons2)

    @pytest.mark.parametrize('electrons1', ELECTRON_COUNTS)
    @pytest.mark.parametrize('electrons2', ELECTRON_COUNTS)
    def test_is_isomorphic_reverse_verdict_is_exactly_its_old_net_comparison(
            self, electrons1, electrons2):
        """Reverse identity used to be ``self.electrons == -other.electrons``.
        Same agreement, sides swapped."""
        a, b = _species('[OH]'), _species('O')
        rxn1 = _reaction([a], [b], electrons=electrons1)
        rxn2 = _reaction([b], [a], electrons=electrons2)
        assert rxn1.is_isomorphic(rxn2, either_direction=True) is (electrons1 == -electrons2)

    @pytest.mark.parametrize('electrons1', ELECTRON_COUNTS)
    @pytest.mark.parametrize('electrons2', ELECTRON_COUNTS)
    def test_the_model_predicate_is_tightened_but_only_for_charged_reactions(
            self, electrons1, electrons2):
        """``are_identical_species_references`` compared NO electron information
        before, so unlike ``is_isomorphic`` this IS a change of verdict -- for
        charged reactions. The bound on it is that two neutral reactions both have
        the counts ``(0, 0)``, so nothing about neutral chemistry moves; the
        neutral corner of this grid is the old behaviour, and the rest is the
        repair."""
        a, b = _species('[OH]'), _species('O')
        rxn1 = _reaction([a], [b], electrons=electrons1)
        rxn2 = _reaction([b], [a], electrons=electrons2)
        expected = electrons1 == -electrons2
        assert are_identical_species_references(rxn1, rxn2) is expected
        if electrons1 == 0 and electrons2 == 0:
            assert expected is True, 'two neutral reactions must still collapse'

    @pytest.mark.parametrize('owner', sorted(o for o, d in FAMILY_ELECTRON_PLACEMENT.items()
                                             if 0 in d))
    def test_one_sided_declared_owners_also_sit_on_the_net_rule(self, owner):
        """A one-sided declaration places exactly where the net rule would, so
        these owners' verdicts do not move either. Only a two-sided declaration
        can change an answer."""
        reactant_count, product_count = FAMILY_ELECTRON_PLACEMENT[owner]
        electrons = product_count - reactant_count
        rxn = _reaction([_species('[Li]')], [_species('[Li+]')],
                        electrons=electrons, family=owner)
        assert get_electron_placement_counts(rxn) == _net_rule_counts(electrons)

    def test_the_two_sided_owner_is_the_only_one_that_moves(self):
        two_sided = sorted(o for o, d in FAMILY_ELECTRON_PLACEMENT.items() if 0 not in d)
        # I-206 added the Plasma_Electron_Impact_Ionization FAMILY at (1, 2), a
        # second two-sided owner mirroring the IONISATION library's pair for the
        # same chemistry. Both "move" -- their placement differs from the net
        # rule -- which is what this test checks via the IONISATION reaction below.
        assert two_sided == [IONISATION, 'Plasma_Electron_Impact_Ionization'], two_sided
        rxn = _reaction([_species('[Li]')], [_species('[Li+]')],
                        electrons=1, family=IONISATION)
        assert get_electron_placement_counts(rxn) == (1, 2)
        assert get_electron_placement_counts(rxn) != _net_rule_counts(1)

    def test_a_reversed_declared_reaction_reports_swapped_counts(self):
        """A declaration is stated in its owner's forward orientation. A reaction
        stored the other way round carries the negated net count, and that is what
        identifies the orientation."""
        rxn = _reaction([_species('[Li+]')], [_species('[Li]')],
                        electrons=-1, family=IONISATION)
        assert get_electron_placement_counts(rxn) == (2, 1)

    def test_a_reaction_contradicting_its_owner_falls_back_rather_than_raising(self):
        """The export boundary refuses such a reaction, loudly, and should. The
        duplicate check is consulted on every comparison of every run, so it reads
        the reaction's own account of itself instead of turning a data
        inconsistency into a traceback from an unrelated place."""
        rxn = _reaction([_species('[Li]')], [_species('[Li+]')],
                        electrons=7, family=IONISATION)
        assert get_electron_placement_counts(rxn) == (0, 7)


class TestGenuineDuplicatesStillCollapse:
    """The duplicate check exists for a real reason, and the repair must not cost
    it. Two reactions that are the same reaction still have to collapse."""

    def test_identical_neutral_reactions_collapse(self):
        a, b = _species('[OH]'), _species('O')
        assert are_identical_species_references(_reaction([a], [b]), _reaction([a], [b]))

    def test_a_neutral_reaction_and_its_reverse_collapse(self):
        a, b = _species('[OH]'), _species('O')
        assert are_identical_species_references(_reaction([a], [b]), _reaction([b], [a]))

    def test_identical_charged_reactions_collapse(self):
        """Same heavy species, same owner, same placement: the same reaction."""
        a, b = _species('[Li]'), _species('[Li+]')
        rxn1 = _reaction([a], [b], electrons=1, family=IONISATION)
        rxn2 = _reaction([a], [b], electrons=1, family=IONISATION)
        assert are_identical_species_references(rxn1, rxn2)
        assert rxn1.is_isomorphic(rxn2, either_direction=True)

    def test_a_charged_reaction_and_its_true_reverse_collapse(self):
        """``A + e- => B`` and ``B => A + e-`` are one reaction seen twice, and
        RMG represents such a pair with a single object. Collapsing them is
        correct and must survive the repair."""
        a, b = _species('[Li+]'), _species('[Li]')
        forward = _reaction([a], [b], electrons=-1, family=RECOMBINATION)
        reverse = _reaction([b], [a], electrons=1, family=RECOMBINATION)
        assert get_electron_placement_counts(forward) == (1, 0)
        assert get_electron_placement_counts(reverse) == (0, 1)
        assert are_identical_species_references(forward, reverse)
        assert forward.is_isomorphic(reverse, either_direction=True)

    def test_a_different_collider_still_separates_otherwise_identical_reactions(self):
        a, b, m = _species('[OH]'), _species('O'), _species('[He]')
        rxn1 = _reaction([a], [b])
        rxn2 = _reaction([a], [b])
        rxn2.specific_collider = m
        assert not are_identical_species_references(rxn1, rxn2)

    def test_a_repeated_reaction_from_a_declared_owner_still_collapses(self):
        """A declared owner does not make two copies of one reaction distinct."""
        a, b = _species('[Li]'), _species('[Li+]')
        first = _reaction([a], [b], electrons=1, family=IONISATION)
        second = _reaction([a], [b], electrons=1, family=IONISATION)
        assert get_electron_placement_counts(first) == get_electron_placement_counts(second)
        assert are_identical_species_references(first, second)


@pytest.mark.database
class TestChargedFamilyReactionsAreNotConfused:
    """A real family-generated pair, out of the database, with no lithium in it.

    This is the artifact that shows the defect was never a property of two plasma
    libraries meeting: it is a property of RMG's notion of reaction identity, and
    it reaches every charged reaction a family produces. The reactions here come
    from ``Plasma_Electron_Attachment`` applied to triplet dioxygen, which is what
    the family is trained on.
    """

    @classmethod
    def setup_class(cls):
        families_path = os.path.join(settings['database.directory'], 'kinetics', 'families')
        cls.database = KineticsDatabase()
        cls.database.load_recommended_families(os.path.join(families_path, 'recommended.py'))
        cls.database.load_families(families_path, families=[ATTACHMENT])
        cls.family = cls.database.families[ATTACHMENT]

    def _attachment(self):
        """Regenerated per test. ``make_new_reaction`` replaces a reaction's lists
        with the model's own Species objects and the identity predicate compares
        by reference, so a shared reaction would carry one test's objects into the
        next and report a verdict that is an artefact of the test."""
        o2 = Species(molecule=[Molecule().from_adjacency_list(O2_ADJACENCY_LIST)])
        o2.generate_resonance_structures()
        reactions = self.family.generate_reactions([o2.molecule])
        assert reactions, 'the attachment family generated nothing from O2'
        return reactions[0]

    def _same_direction(self, rxn, electrons):
        return TemplateReaction(reactants=list(rxn.reactants), products=list(rxn.products),
                                electrons=electrons, family=rxn.family,
                                template=getattr(rxn, 'template', None),
                                degeneracy=rxn.degeneracy)

    def test_the_family_is_charged_and_its_reactions_carry_the_declaration(self):
        attach = self._attachment()
        assert self.family.electrons == -1
        assert attach.electrons == -1
        assert get_electron_placement_counts(attach) == (1, 0)

    def test_its_true_reverse_still_collapses(self):
        """The control. A reaction and its reverse are one reaction, and the
        check must keep saying so."""
        attach = self._attachment()
        reverse = TemplateReaction(reactants=list(attach.products),
                                   products=list(attach.reactants),
                                   electrons=-attach.electrons, family=attach.family,
                                   template=getattr(attach, 'template', None),
                                   degeneracy=attach.degeneracy)
        assert are_identical_species_references(attach, reverse)

    def test_the_same_direction_charge_reversal_survives(self):
        """``O2 + e- => O2-`` and ``O2 => O2- + e-`` are different reactions --
        different molecularity, opposite charge transfer -- and used to collapse."""
        attach = self._attachment()
        flipped = self._same_direction(attach, electrons=-attach.electrons)
        assert get_electron_placement_counts(attach) == (1, 0)
        assert get_electron_placement_counts(flipped) == (0, 1)
        assert not are_identical_species_references(attach, flipped)

    def test_the_neutral_transformation_over_the_same_heavy_species_survives(self):
        """The shape that needs no plasma library at all: a neutral reaction whose
        heavy species happen to coincide with a charged one's."""
        attach = self._attachment()
        neutral = self._same_direction(attach, electrons=0)
        assert get_electron_placement_counts(neutral) == (0, 0)
        assert not are_identical_species_references(attach, neutral)

    def test_the_model_keeps_both_through_the_real_path(self):
        with _GlobalKineticsDatabase(self.database):
            model = CoreEdgeReactionModel()
            model.kinetics_database = self.database
            _, first_new = model.make_new_reaction(self._attachment(), generate_thermo=False,
                                                   generate_kinetics=False)
            assert first_new is True
            flipped = self._same_direction(self._attachment(), electrons=1)
            _, second_new = model.make_new_reaction(flipped, generate_thermo=False,
                                                    generate_kinetics=False)
        assert second_new is True, ('the same-direction charge reversal was discarded as a '
                                    'duplicate of the attachment')


@pytest.mark.database
class TestLithiumChargeNetworkReachesTheModel:
    """The real mechanism, through the real model-building path.

    ``add_reaction_library_to_edge`` is the step that dropped the sink, and
    ``CoreEdgeReactionModel`` is the model the assertion is about: the requirement
    is that the recombination is *in the model*, not merely that it loads.
    """

    @classmethod
    def setup_class(cls):
        cls.libraries_path = os.path.join(settings['database.directory'],
                                          'kinetics', 'libraries')

    def _database(self):
        database = KineticsDatabase()
        database.load_libraries(self.libraries_path, libraries=[IONISATION, RECOMBINATION])
        return database

    def _sole_reaction(self, database, label):
        reactions = database.libraries[label].get_library_reactions()
        assert len(reactions) == 1, '{0} has {1} entries, expected 1'.format(label,
                                                                            len(reactions))
        return reactions[0]

    def test_the_two_channels_are_not_a_reverse_pair(self):
        """The measurement the whole repair turns on. Their net counts are equal
        and opposite -- the relation a reverse pair has -- but their placements
        are not mirror images, because the reverse of the ionisation is the
        three-body channel and not the radiative one."""
        database = self._database()
        ionisation = self._sole_reaction(database, IONISATION)
        recombination = self._sole_reaction(database, RECOMBINATION)

        assert ionisation.electrons == -recombination.electrons

        ion_counts = get_electron_placement_counts(ionisation)
        rec_counts = get_electron_placement_counts(recombination)
        assert ion_counts == (1, 2)
        assert rec_counts == (1, 0)
        assert (ion_counts[1], ion_counts[0]) != rec_counts

    @pytest.mark.parametrize('source_first', [True, False])
    def test_both_channels_enter_the_model(self, source_first):
        """Both orders, because the drop was order-dependent: whichever was
        offered second was the one that disappeared."""
        database = self._database()
        order = [IONISATION, RECOMBINATION] if source_first else [RECOMBINATION, IONISATION]
        verdicts = {}
        with _GlobalKineticsDatabase(database):
            model = CoreEdgeReactionModel()
            model.kinetics_database = database
            for label in order:
                _, is_new = model.make_new_reaction(self._sole_reaction(database, label),
                                                    generate_thermo=False,
                                                    generate_kinetics=False)
                verdicts[label] = is_new

        assert verdicts[IONISATION] is True, 'the ionisation source did not enter the model'
        assert verdicts[RECOMBINATION] is True, 'the recombination sink did not enter the model'

    def test_the_cation_has_a_loss_channel_in_the_model(self):
        """The deliverable, stated as chemistry rather than as a count: the cation
        appears on the reactant side of some reaction in the model."""
        database = self._database()
        with _GlobalKineticsDatabase(database):
            model = CoreEdgeReactionModel()
            model.kinetics_database = database
            for label in (IONISATION, RECOMBINATION):
                model.make_new_reaction(self._sole_reaction(database, label),
                                        generate_thermo=False, generate_kinetics=False)

        cation_labels = set()
        for rxn in model.new_reaction_list:
            for product in rxn.products:
                if product.molecule and product.molecule[0].get_net_charge() > 0:
                    cation_labels.add(product.label)
        assert cation_labels, 'no cation is produced at all'

        consumed = set()
        for rxn in model.new_reaction_list:
            for reactant in rxn.reactants:
                if reactant.label in cation_labels:
                    consumed.add(reactant.label)
        assert consumed == cation_labels, (
            'these cations are produced but never consumed: {0}'.format(
                sorted(cation_labels - consumed))
        )


@pytest.mark.database
class TestTheLibraryLoadLandmineIsNotReachedFromHere:
    """A landmine adjacent to this repair, measured and NOT fixed by it.

    ``KineticsLibrary.check_for_duplicates`` runs inside ``load`` and raises
    ``DatabaseError`` on any unmarked duplicate, deciding duplication with
    ``Reaction.is_isomorphic``. So it looked as though repairing that predicate
    would also let one library carry both channels. Measured, it does not, and the
    reason is upstream of anything this change touches: the check runs over
    ``entry.item``, which ``KineticsLibrary.load_entry`` builds as a plain
    :class:`Reaction` **with no owner recorded at all**. No owner means no
    placement declaration, which means the net-derived rule -- under which the two
    channels are ``(0, 1)`` and ``(1, 0)``, exact mirrors, and therefore still
    duplicates.

    The owner only appears later, in ``get_library_reactions``, which wraps each
    entry in a :class:`LibraryReaction` whose ``family`` is the library label.
    Everything downstream of that point -- the model builder, the reactor, both
    writers -- sees the declaration and is repaired. The library loader is the one
    consumer that sees the reaction before its owner is attached.

    The verdict is also **unchanged** by this commit rather than broken by it: the
    old code compared ``electrons == -electrons`` (``+1`` against ``-(-1)``, true)
    and the new code compares ``(0, 1)`` against the mirror of ``(1, 0)`` (equal,
    true). Same answer, same reason.

    Fixing it means giving library entries their owner at load, and it needs a
    design decision this ticket has no mandate for: the placement registry is
    keyed per owner, so one library cannot declare two different placements even
    once the owner is recorded. That is a separate ticket, and the strict xfail
    below is what will tell whoever takes it that they have succeeded.

    The library here is not a fixture: its entries are the shipped
    ``PlasmaElectronImpactIonization`` and ``PlasmaRadiativeRecombination``
    entries, loaded out of the database by RMG's own loader.
    """

    @classmethod
    def setup_class(cls):
        cls.libraries_path = os.path.join(settings['database.directory'],
                                          'kinetics', 'libraries')

    def _database(self):
        database = KineticsDatabase()
        database.load_libraries(self.libraries_path, libraries=[IONISATION, RECOMBINATION])
        return database

    def _merged_library(self, database):
        merged = KineticsLibrary(label='PlasmaBothChannels')
        merged.entries = {}
        index = 0
        for label in (IONISATION, RECOMBINATION):
            for entry in database.libraries[label].entries.values():
                index += 1
                entry.index = index
                merged.entries[index] = entry
        assert len(merged.entries) == 2
        return merged

    def test_a_loaded_entry_carries_no_owner(self):
        """The measurement the whole limitation rests on."""
        database = self._database()
        for label in (IONISATION, RECOMBINATION):
            for entry in database.libraries[label].entries.values():
                assert getattr(entry.item, 'family', None) is None
                assert get_placement_declaration(entry.item) is None

    def test_the_owner_appears_when_the_entry_becomes_a_library_reaction(self):
        """And from there on, the placement is available and this repair applies."""
        database = self._database()
        placements = {}
        for label in (IONISATION, RECOMBINATION):
            reactions = database.libraries[label].get_library_reactions()
            assert len(reactions) == 1
            assert reactions[0].family == label
            placements[label] = get_electron_placement_counts(reactions[0])
        assert placements == {IONISATION: (1, 2), RECOMBINATION: (1, 0)}

    @pytest.mark.xfail(strict=True, reason=(
        'Not reached by this repair: check_for_duplicates runs over entry.item, which '
        'carries no owner, so both channels fall back to the net rule and stay mirrors. '
        'Fixing it requires recording the owner at load, and a decision about a registry '
        'that is keyed one placement per owner. When this starts passing, that ticket has '
        'landed and this marker should be removed.'))
    def test_one_library_carrying_both_channels_can_be_loaded(self):
        self._merged_library(self._database()).check_for_duplicates()
