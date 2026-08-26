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
I-116: the electron-impact ionisation path, from library load to reactor
acceptance.

The production diff this file exists for is one line --
``'PlasmaElectronImpactIonization': (1, 2)`` in
:data:`rmgpy.electron_placement.FAMILY_ELECTRON_PLACEMENT` -- and the reason it
takes a file to defend is that the line is the *only* thing standing between a
plasma mechanism that can create its own charge and one that cannot. Every stage
below fails independently, and three of them fail quietly.

Six stages, asserted separately:

1. the ``PlasmaElectronImpactIonization`` kinetics library **loads**, picking its
   Voronov fit up out of ``voronov.yaml`` by ``(Z, N)``;
2. the canonical entry ``Li => Li+`` **balances**, the 0 -> +1 charge gap closed
   by the ``electrons = +1`` the rate law propagates onto the reaction;
3. it **resolves** through :func:`resolve_electron_placement` into
   ``Li + e- => Li+ + 2 e-``, leaving the canonical reaction untouched;
4. the declaration's net change is **validated against** ``reaction.electrons``
   -- and a perturbed scalar is refused;
5. the **rate-order cross-check resolves** to a real integer and compares it
   against the view's reactant count. ``get_plasma_rate_order`` answers ``None``
   for any rate law outside its units table, and a ``None`` that read as
   agreement would make this whole file decorative, so the resolvability is
   asserted on its own;
6. :class:`~rmgpy.solver.plasma.PlasmaReactor` **accepts** it through its own
   ``initialize_model`` validation path -- not a stand-in -- and prices it at
   the Voronov coefficient evaluated at ``Te``, which is the assertion that
   catches a view accepted at the wrong reaction order (accepting it is not
   enough: a wrong-order view does not raise).

Two ways this could pass without being true are guarded explicitly:
``TestRateOrderCrossCheckResolved`` (a silent ``None`` must not read as
agreement) and ``TestOrderComesFromDeclarationNotScalar`` (the reactant count
must come from the declaration, never from anything derived from the reaction's
own electron count).

The three negative controls are ``TestNeighbouringDeclarationsUnperturbed``,
``TestUncoveredSpeciesGetsNothing`` and ``TestWrongDeclarationIsRefused``.

DATABASE DEPENDENCE. These tests need a database carrying the
``PlasmaElectronImpactIonization`` library. If ``settings['database.directory']``
does not have it, every test here SKIPS with that reason named -- it does not
inject the library, and it does not pass. A green run of this file is therefore
evidence about the runtime's database as well as about RMG-Py.
"""

import os

import pytest

from rmgpy import settings
from rmgpy.data.kinetics import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.electron_balance import get_plasma_rate_order
from rmgpy.electron_placement import (
    FAMILY_ELECTRON_PLACEMENT,
    RATE_ORDER_AGREES,
    RATE_ORDER_UNRESOLVABLE,
    _place_declared_electrons,
    resolve_electron_placement,
)
from rmgpy.exceptions import ElectronPlacementError
from rmgpy.kinetics import Arrhenius
from rmgpy.kinetics.arrhenius import VoronovEIArrhenius
from rmgpy.reaction import Reaction
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.species import Species

################################################################################

#: The kinetics library that owns electron-impact ionisation, and -- because
#: ``LibraryReaction.__init__`` sets ``self.family = library`` -- also its
#: placement key.
LIBRARY = 'PlasmaElectronImpactIonization'

#: The declaration this ticket ships: one electron incident, two liberated.
DECLARATION = (1, 2)

T_GAS = 500.0        # K
T_E = 11604.5        # K, ~1 eV, the bottom of the Voronov fit's validity range
P0 = 1.0e5           # Pa
Y_E0 = 1.0e-4        # mol of seeded electrons

#: Attachment and cation-recombination rate coefficients as measured on the base
#: commit, before this ticket's line. Negative control 1 compares against these
#: literals rather than against a recomputation.
BASE_ATTACHMENT_K_1000 = 1000000.0000000001
BASE_RECOMBINATION_K_1000 = 10000.000000000002


def _libraries_path():
    return os.path.join(settings['database.directory'], 'kinetics', 'libraries')


def _electron():
    return Species(label='e').from_adjacency_list('1 e u1 p0 c-1')


def _n_electrons(participants):
    return sum(1 for spc in participants if spc.is_electron())


@pytest.fixture(scope='module')
def library():
    """The real ionisation library, loaded from the configured database.

    Skips -- rather than injecting anything -- when the configured database has
    no such library, so a green run of this file always means the path was
    exercised against real data.
    """
    path = os.path.join(_libraries_path(), LIBRARY)
    if not os.path.isdir(path):
        pytest.skip(
            'kinetics library {0!r} is not present in the configured database '
            '({1}); this suite exercises the real library and does not inject '
            'a substitute.'.format(LIBRARY, settings['database.directory']))
    database = KineticsDatabase()
    database.load_libraries(_libraries_path(), libraries=[LIBRARY])
    return database.libraries[LIBRARY]


@pytest.fixture
def canonical(library):
    """A fresh canonical ``Li => Li+`` off the library.

    Fresh per test on purpose: the reactor indexes species by identity, and the
    negative controls mutate the registry around this object.
    """
    reactions = library.get_library_reactions()
    assert len(reactions) == 1
    return reactions[0]


@pytest.fixture
def declare_wrongly():
    """Temporarily replace the shipped ionisation declaration, then restore it.

    Restoration is unconditional (``try/finally`` inside the fixture teardown
    contract), which is what lets a test prove a wrong declaration is refused
    without leaving the registry wrong for every test that follows.
    """
    saved = FAMILY_ELECTRON_PLACEMENT[LIBRARY]

    def _declare(value):
        FAMILY_ELECTRON_PLACEMENT[LIBRARY] = value

    yield _declare
    FAMILY_ELECTRON_PLACEMENT[LIBRARY] = saved


################################################################################
# The premise the one-line approach rests on
################################################################################

class TestLibraryLabelIsAPlacementKey:
    """``LibraryReaction`` sets ``family = library``, which is the only reason a
    LIBRARY name can be a key in a mapping called ``FAMILY_ELECTRON_PLACEMENT``.
    Asserted rather than taken on trust: if it stopped holding, this ticket's
    whole approach would need revisiting rather than patching."""

    def test_library_reaction_reports_its_library_as_its_family(self):
        reaction = LibraryReaction(reactants=[], products=[], library=LIBRARY)
        assert reaction.family == LIBRARY

    def test_the_shipped_declaration_is_keyed_by_that_label(self):
        assert FAMILY_ELECTRON_PLACEMENT[LIBRARY] == DECLARATION


################################################################################
# Stage 1 -- the library loads
################################################################################

@pytest.mark.database
class TestStage1LibraryLoads:

    def test_library_carries_exactly_the_one_covered_reaction(self, library):
        reactions = library.get_library_reactions()
        assert len(reactions) == 1
        assert str(reactions[0]) == '[Li] => [Lip]'

    def test_entry_carries_the_real_voronov_fit(self, canonical):
        """Read out of ``voronov.yaml`` by ``(Z, N)``, not transcribed. ``Z``
        and ``N`` are constructor arguments rather than stored attributes on the
        cdef class, so the stage is identified from the comment the loader
        writes and from the coefficients it produced."""
        kinetics = canonical.kinetics
        assert isinstance(kinetics, VoronovEIArrhenius)
        assert 'Z=3, N=3' in kinetics.comment
        assert kinetics.dE == pytest.approx(5.4, rel=1e-6)   # eV, the Li I threshold
        assert kinetics.uses_electron_temperature
        assert kinetics.uses_electron_density
        assert kinetics.get_rate_coefficient_electron_temp(10000.0) > 0.0

    def test_the_loader_propagates_the_rate_law_electron_count(self, canonical):
        """``load_entry`` takes no electron argument, so the count reaches the
        reaction only because ``KineticsLibrary.load`` copies it off the rate
        law. Without that, stage 2 fails."""
        assert canonical.electrons == 1

    def test_canonical_form_carries_no_explicit_electron(self, canonical):
        assert _n_electrons(canonical.reactants) == 0
        assert _n_electrons(canonical.products) == 0

    def test_the_reaction_is_attributed_to_the_library(self, canonical):
        assert canonical.family == LIBRARY


################################################################################
# Stage 2 -- the reaction balances
################################################################################

@pytest.mark.database
class TestStage2ReactionBalances:

    def test_canonical_reaction_balances(self, canonical):
        """``Li => Li+`` closes its 0 -> +1 charge gap through the metadata
        electron count, which is what makes the library loadable at all --
        ``KineticsLibrary.load`` raises ``DatabaseError`` on an unbalanced
        entry."""
        assert canonical.is_balanced() is True


################################################################################
# Stage 3 -- it resolves through electron placement
################################################################################

@pytest.mark.database
class TestStage3ResolvesThroughPlacement:

    def test_view_has_one_reactant_and_two_product_electrons(self, canonical):
        view = resolve_electron_placement(canonical, [_electron()])
        assert _n_electrons(view.reactants) == 1
        assert _n_electrons(view.products) == 2

    def test_view_is_the_reactor_facing_form(self, canonical):
        view = resolve_electron_placement(canonical, [_electron()])
        assert view.electrons == 0
        assert len(view.reactants) == 2      # Li + e-, i.e. order 2
        assert len(view.products) == 3       # Li+ + 2 e-
        assert view.kinetics is canonical.kinetics
        assert view.reversible is canonical.reversible

    def test_resolution_leaves_the_canonical_reaction_untouched(self, canonical):
        reactants_before = list(canonical.reactants)
        products_before = list(canonical.products)

        resolve_electron_placement(canonical, [_electron()])

        assert canonical.electrons == 1
        assert canonical.reactants == reactants_before
        assert canonical.products == products_before
        assert _n_electrons(canonical.reactants) == 0
        assert _n_electrons(canonical.products) == 0

    def test_view_balances_in_the_electron_pseudo_element(self, canonical):
        """Step 10 of the resolver would have refused otherwise; asserting the
        outcome directly means this file does not depend on the resolver's own
        internal check to notice a wrong-sided placement."""
        from rmgpy.electron_balance import get_species_electron_count
        view = resolve_electron_placement(canonical, [_electron()])
        left = sum(get_species_electron_count(s) for s in view.reactants)
        right = sum(get_species_electron_count(s) for s in view.products)
        assert left == right


################################################################################
# Stage 4 -- the declared count validates against the reaction's own scalar
################################################################################

@pytest.mark.database
class TestStage4ScalarValidatesAgainstDeclaration:

    def test_scalar_equals_the_declared_net_change(self, canonical):
        reactant_count, product_count = FAMILY_ELECTRON_PLACEMENT[LIBRARY]
        assert canonical.electrons == product_count - reactant_count

    @pytest.mark.parametrize('wrong_scalar,shape', [
        (-1, 'attachment-shaped'),
        (0, 'excitation-shaped'),
        (2, 'ionization-shaped'),
    ])
    def test_a_perturbed_scalar_is_refused_by_name(self, canonical, wrong_scalar, shape):
        """The declaration is untouched; only the reaction's own count moves.
        The refusal names the shape it saw and the net the declaration wants."""
        canonical.electrons = wrong_scalar
        with pytest.raises(ElectronPlacementError) as excinfo:
            resolve_electron_placement(canonical, [_electron()])
        message = str(excinfo.value)
        assert shape in message
        assert LIBRARY in message
        assert 'declares net 1' in message

    def test_a_disagreeing_rate_law_electron_count_is_refused(self, canonical):
        """``VoronovEIArrhenius`` carries its own net count, and it is a THIRD
        source that must agree -- not a second one that must be absent."""
        canonical.kinetics = VoronovEIArrhenius(Z=3, N=3, electrons=2)
        with pytest.raises(ElectronPlacementError,
                           match='declaring a net electron change of 2'):
            resolve_electron_placement(canonical, [_electron()])


################################################################################
# Stage 5 -- the rate-order cross-check RESOLVED (guard: a None is not agreement)
################################################################################

@pytest.mark.database
class TestRateOrderCrossCheckResolved:
    """``get_plasma_rate_order`` returns ``None`` for any rate law whose units
    are not in its lookup table, and a silent ``None`` reported as success is
    worse than no check at all. So the resolvability is asserted on its own,
    before anything is concluded from the agreement."""

    def test_order_resolves_to_a_real_integer_for_this_kinetics_class(self, canonical):
        order = get_plasma_rate_order(canonical.kinetics)
        assert order is not None, (
            'the cross-check abstained for {0}; a None must never read as '
            'agreement'.format(type(canonical.kinetics).__name__))
        assert isinstance(order, int) and not isinstance(order, bool)
        assert order == 2

    def test_the_check_compared_two_numbers(self, canonical):
        """Both sides of the comparison, computed independently here."""
        order = get_plasma_rate_order(canonical.kinetics)
        view = resolve_electron_placement(canonical, [_electron()])
        assert order == 2
        assert len(view.reactants) == 2
        assert order == len(view.reactants)

    def test_the_view_records_agreement_and_not_abstention(self, canonical):
        view = resolve_electron_placement(canonical, [_electron()])
        assert RATE_ORDER_AGREES in view.comment
        assert RATE_ORDER_UNRESOLVABLE not in view.comment
        assert '(order 2)' in view.comment

    def test_an_unresolvable_order_is_reported_as_its_own_outcome(self, canonical):
        """The contrast case, so the assertion above is known to discriminate:
        a plain Arrhenius is outside the units table, and the view says so
        rather than claiming agreement. (Reached with a declaration of
        ``(1, 0)`` because a thermal rate law carries no electron count of its
        own; the point is the OUTCOME STRING, not the chemistry.)"""
        thermal = TemplateReaction(
            reactants=[Species(label='O2').from_smiles('[O][O]')],
            products=[Species(label='O2-').from_smiles('[O][O-]')],
            kinetics=Arrhenius(A=(1.0e12, 'cm^3/(mol*s)'), n=0.0, Ea=(0.0, 'kJ/mol')),
            family='Plasma_Electron_Attachment',
            reversible=False, is_forward=True, electrons=-1)
        assert get_plasma_rate_order(thermal.kinetics) is None
        view = resolve_electron_placement(thermal, [_electron()])
        assert RATE_ORDER_UNRESOLVABLE in view.comment
        assert RATE_ORDER_AGREES not in view.comment


################################################################################
# Guard -- the order comes from the DECLARATION, never from the scalar
################################################################################

@pytest.mark.database
class TestOrderComesFromDeclarationNotScalar:
    """The scalar is validated AGAINST the declaration and read for nothing
    else. If the reactant count were derived from ``reaction.electrons``, a
    ``(1, 2)`` and a ``(0, 1)`` declaration would be indistinguishable -- both
    are net +1 -- and the ``(0, 1)`` reading is exactly the factor-of-electron-
    density error the two-sided declaration exists to prevent."""

    def test_the_placement_primitive_never_reads_the_scalar(self, canonical):
        """Called with the scalar deliberately zeroed, it still places 1 and 2."""
        canonical.electrons = 0
        reactants, products = _place_declared_electrons(canonical, _electron(), 1, 2)
        assert _n_electrons(reactants) == 1
        assert _n_electrons(products) == 2

    @pytest.mark.parametrize('declaration,expected_reactants', [
        ((0, 1), 1),    # same net +1, one reactant -- the silent-factor error
        ((2, 3), 3),    # same net +1, three reactants
    ])
    def test_the_reactant_count_tracks_the_declaration_with_the_scalar_fixed(
            self, canonical, declare_wrongly, declaration, expected_reactants):
        """``reaction.electrons`` stays at +1 throughout, and is the CORRECT net
        for every declaration below. Only the declaration moves, and the view's
        reactant count moves with it -- reported by the cross-check's own
        refusal, which is what catches these."""
        declare_wrongly(declaration)
        assert canonical.electrons == 1
        with pytest.raises(ElectronPlacementError) as excinfo:
            resolve_electron_placement(canonical, [_electron()])
        assert 'built a view with {0:d} reactant(s)'.format(expected_reactants) \
            in str(excinfo.value)


################################################################################
# Stage 6 -- the plasma reactor accepts it
################################################################################

@pytest.mark.database
class TestStage6ReactorAccepts:
    """Through :meth:`PlasmaReactor.initialize_model`, which is the reactor's
    own validation path and the one production uses -- never a stand-in. The
    reactor resolves the placement itself at its generation-to-reactor handoff,
    so these tests hand it the CANONICAL reaction and never call the resolver."""

    @staticmethod
    def _reactor(canonical, electron):
        imf = {electron: Y_E0, canonical.reactants[0]: 1.0}
        return PlasmaReactor(T_GAS, P0, imf, (T_E, 'K'), n_sims=1, termination=[])

    @staticmethod
    def _core(canonical, electron):
        return [electron, canonical.reactants[0], canonical.products[0]]

    def test_reactor_accepts_the_canonical_reaction(self, canonical):
        electron = _electron()
        reactor = self._reactor(canonical, electron)
        reactor.initialize_model(self._core(canonical, electron), [canonical], [], [])
        assert reactor.kf[0] > 0.0

    def test_the_accepted_rate_is_the_voronov_coefficient_at_te(self, canonical):
        """Acceptance alone would not be worth asserting: a view accepted at the
        WRONG reaction order does not raise. Comparing the packed kf against the
        rate law's own electron-temperature evaluator is what makes the order
        visible, and it also proves the kinetics survived the handoff rather
        than being collapsed to a one-temperature k(T, P)."""
        electron = _electron()
        reactor = self._reactor(canonical, electron)
        reactor.initialize_model(self._core(canonical, electron), [canonical], [], [])
        expected = canonical.kinetics.get_rate_coefficient_electron_temp(T_E)
        assert reactor.kf[0] == pytest.approx(expected, rel=1e-12)

    def test_no_reverse_rate_is_manufactured(self, canonical):
        electron = _electron()
        reactor = self._reactor(canonical, electron)
        reactor.initialize_model(self._core(canonical, electron), [canonical], [], [])
        assert reactor.kb[0] == 0.0

    def test_the_reactor_leaves_the_canonical_reaction_untouched(self, canonical):
        electron = _electron()
        reactants_before = list(canonical.reactants)
        reactor = self._reactor(canonical, electron)
        reactor.initialize_model(self._core(canonical, electron), [canonical], [], [])
        assert canonical.electrons == 1
        assert canonical.reactants == reactants_before
        assert _n_electrons(canonical.reactants) == 0

    def test_the_reactors_own_validation_is_live_on_this_reaction(self, canonical):
        """The same reaction, marked reversible, is refused by name inside
        ``_validate_reactions``. Without this, "the reactor accepted it" could
        mean the reactor never looked."""
        from rmgpy.exceptions import NonEquilibriumReverseRateError
        electron = _electron()
        canonical.reversible = True
        reactor = self._reactor(canonical, electron)
        with pytest.raises(NonEquilibriumReverseRateError,
                           match='reversal is undefined'):
            reactor.initialize_model(self._core(canonical, electron), [canonical], [], [])


################################################################################
# Negative control 1 -- the neighbouring declarations are unperturbed
################################################################################

class TestNeighbouringDeclarationsUnperturbed:
    """Adding a key to a dict should not perturb the others. Asserted rather
    than assumed, because "the shape I was working on still works" is the
    cheapest thing to check and the easiest to skip.

    Both shapes are resolved twice -- once with the ionisation declaration
    present and once with it removed -- and compared participant by
    participant, so the control does not depend on remembering what the base
    commit produced. The rate coefficients are ALSO compared against literals
    measured on the base commit, which does."""

    @staticmethod
    def _attachment():
        return TemplateReaction(
            reactants=[Species(label='O2').from_smiles('[O][O]')],
            products=[Species(label='O2-').from_smiles('[O][O-]')],
            family='Plasma_Electron_Attachment', electrons=-1,
            reversible=False, is_forward=True,
            kinetics=Arrhenius(A=(1.0e12, 'cm^3/(mol*s)'), n=0.0, Ea=(0.0, 'kJ/mol')))

    @staticmethod
    def _recombination():
        return TemplateReaction(
            reactants=[Species(label='Li+').from_adjacency_list('1 Li u0 p0 c+1'),
                       Species(label='CH3').from_smiles('[CH3]')],
            products=[Species(label='CH3Li').from_smiles('C[Li]')],
            family='Cation_R_Recombination', electrons=-1,
            reversible=False, is_forward=True,
            kinetics=Arrhenius(A=(1.0e10, 'cm^3/(mol*s)'), n=0.0, Ea=(0.0, 'kJ/mol')))

    @staticmethod
    def _fingerprint(builder):
        view = resolve_electron_placement(builder(), [_electron()])
        return (
            [str(s) for s in view.reactants],
            [str(s) for s in view.products],
            _n_electrons(view.reactants),
            _n_electrons(view.products),
            view.electrons,
            view.reversible,
            view.comment,
        )

    @pytest.mark.parametrize('builder_name', ['_attachment', '_recombination'])
    def test_view_is_identical_with_and_without_the_new_entry(self, builder_name):
        builder = getattr(self, builder_name)
        with_entry = self._fingerprint(builder)
        removed = FAMILY_ELECTRON_PLACEMENT.pop(LIBRARY)
        try:
            without_entry = self._fingerprint(builder)
        finally:
            FAMILY_ELECTRON_PLACEMENT[LIBRARY] = removed
        assert with_entry == without_entry

    def test_attachment_still_places_one_reactant_electron_and_none_produced(self):
        view = resolve_electron_placement(self._attachment(), [_electron()])
        assert _n_electrons(view.reactants) == 1
        assert _n_electrons(view.products) == 0
        assert view.kinetics.get_rate_coefficient(1000.0) == pytest.approx(
            BASE_ATTACHMENT_K_1000, rel=1e-12)

    def test_recombination_still_places_one_reactant_electron_and_none_produced(self):
        view = resolve_electron_placement(self._recombination(), [_electron()])
        assert _n_electrons(view.reactants) == 1
        assert _n_electrons(view.products) == 0
        assert view.kinetics.get_rate_coefficient(1000.0) == pytest.approx(
            BASE_RECOMBINATION_K_1000, rel=1e-12)


################################################################################
# Negative control 2 -- an uncovered species gets nothing
################################################################################

@pytest.mark.database
class TestUncoveredSpeciesGetsNothing:
    """No match, no estimate, no fabricated rate. A library has no matcher, so
    this holds structurally -- but "it holds structurally" is exactly the claim
    that is never checked, so it is checked here."""

    def test_the_library_exposes_no_group_or_template_matcher(self, library):
        assert not hasattr(type(library), 'groups')
        assert not hasattr(type(library), 'forward_template')
        assert not hasattr(type(library), 'generate_reactions')

    @pytest.mark.parametrize('label,smiles', [('Ar', '[Ar]'), ('He', '[He]')])
    def test_an_uncovered_species_gets_no_reaction(self, library, label, smiles):
        database = KineticsDatabase()
        species = Species(label=label).from_smiles(smiles)
        species.generate_resonance_structures()
        assert database.generate_reactions_from_library(library, [species]) == []

    def test_the_covered_species_does_get_one_through_the_same_call(self, library):
        """Without this, the emptiness above would be equally consistent with a
        broken call. (Note the SINGULAR ``generate_reactions_from_library``: the
        plural iterates ``KineticsDatabase.library_order``, which
        ``load_libraries(path, libraries=[...])`` never populates, so it returns
        ``[]`` for every species including this one.)"""
        database = KineticsDatabase()
        species = Species(label='Li').from_adjacency_list('1 Li u1 p0 c0')
        found = database.generate_reactions_from_library(library, [species])
        assert len(found) == 1

    def test_the_lookup_route_does_not_resolve_but_fails_loudly(self, library):
        """A finding, pinned so it cannot regress into silence.

        ``KineticsDatabase.generate_reactions_from_library`` constructs
        ``LibraryReaction(library=library)`` -- the OBJECT -- while
        ``KineticsLibrary.get_library_reactions`` passes ``self.label``. Only
        the label is a registry key, so a reaction obtained through the lookup
        route carries a ``KineticsLibrary`` instance as its ``family`` and does
        NOT resolve. The route that matters -- the one
        ``CoreEdgeReactionModel.add_reaction_library_to_edge`` uses to bring a
        ``kineticsLibraries`` entry into a model -- is
        ``get_library_reactions``, which is why stage 3 passes.

        What is asserted here is that the other route fails by NAME rather than
        by silently resolving to something wrong. (It also loses the electron
        count, since only ``KineticsLibrary.load`` propagates it, so the shape
        would be wrong twice over.)"""
        database = KineticsDatabase()
        species = Species(label='Li').from_adjacency_list('1 Li u1 p0 c0')
        via_object = database.generate_reactions_from_library(library, [species])[0]
        assert via_object.family is library          # the object, not the label
        assert not isinstance(via_object.family, str)
        with pytest.raises(ElectronPlacementError,
                           match='has no electron-placement declaration'):
            resolve_electron_placement(via_object, [_electron()])


################################################################################
# Negative control 3 -- a deliberately wrong declaration is refused
################################################################################

@pytest.mark.database
class TestWrongDeclarationIsRefused:
    """What proves the validation is live rather than decorative: the
    declaration is temporarily made wrong, the resolver is shown refusing by
    name against the reaction's ACTUAL electron count, and the shipped value is
    then shown to be back and still working."""

    @pytest.mark.parametrize('wrong', [(1, 1), (2, 2)])
    def test_a_wrong_declaration_is_refused_against_the_actual_count(
            self, canonical, declare_wrongly, wrong):
        declare_wrongly(wrong)
        with pytest.raises(ElectronPlacementError) as excinfo:
            resolve_electron_placement(canonical, [_electron()])
        message = str(excinfo.value)
        assert LIBRARY in message
        assert 'carries electrons=1' in message
        assert 'ionization-shaped (net electron production)' in message
        assert 'declares net 0' in message

    def test_the_shipped_declaration_is_back_and_still_resolves(self, canonical):
        """Runs after the parametrized refusals above and depends on their
        teardown having restored the registry -- which is the point."""
        assert FAMILY_ELECTRON_PLACEMENT[LIBRARY] == DECLARATION
        view = resolve_electron_placement(canonical, [_electron()])
        assert _n_electrons(view.reactants) == 1
        assert _n_electrons(view.products) == 2
