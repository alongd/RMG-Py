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
I-113: the two-sided electron-placement declaration.

The declaration used to be ``(side, count)``, which works only because
attachment is the degenerate case where the electron's contribution to reaction
ORDER and its NET stoichiometric change are the same number. Ionisation
separates them: ``Li + e- -> Li+ + 2e-`` has incident order 1 -- which is what
makes the ``cm^3/(molecule*s)`` Voronov coefficient dimensionally meaningful --
and net production +1. A bare ``electrons = +1`` cannot tell that apart from
``Li -> Li+ + e-``, which is first order and would need ``s^-1``; get it wrong
and the rate is off by a factor of the electron density while the file looks
well formed.

The declaration is now ``(reactant_count, product_count)``. Incident order is
``reactant_count``, declared rather than derived; net is
``product_count - reactant_count``, which is the only thing
``Reaction.electrons`` is ever checked against.

These tests pin, on real ``voronov.yaml`` / ``badnell.yaml`` kinetics wherever a
rate law is load-bearing:

* the ionisation view resolving at all (it could not, before);
* the scalar NOT being the order source -- one net count, two declarations, two
  different views;
* the rate-order cross-check refusing the wrong-order case by name, and
  reporting "not resolvable" as its OWN outcome rather than as agreement;
* step 10, the structural ``E`` pseudo-element verification, still passing for
  attachment AND ionisation and still biting a wrong-sided declaration;
* the new two-sided list-append primitive agreeing with ``expand_electrons``
  everywhere ``expand_electrons`` can speak -- which is how the placement path
  and the export path stay in step;
* attachment and cation recombination coming out numerically unchanged.

No ionisation family is added to :data:`FAMILY_ELECTRON_PLACEMENT` here: which
family (or library) should own Li ionisation is a data question and a separate
ticket. The declaration is injected under a synthetic label, as the I-108
measurement did. (That separate ticket has since answered: I-116 declared the
``PlasmaElectronImpactIonization`` kinetics library at ``(1, 2)``. This file
keeps its synthetic label and its injection, so what it tests stays the code
path rather than the shipped data.)
"""

import pytest

from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.electron_balance import expand_electrons, get_species_electron_count
from rmgpy.electron_placement import (
    FAMILY_ELECTRON_PLACEMENT,
    RATE_ORDER_AGREES,
    RATE_ORDER_UNRESOLVABLE,
    _place_declared_electrons,
    resolve_electron_placement,
)
from rmgpy.exceptions import AtomTypeError, ElectronPlacementError
from rmgpy.kinetics.arrhenius import (Arrhenius, BadnellRRArrhenius,
                                      TwoTemperaturePlasma, VoronovEIArrhenius)
from rmgpy.reaction import Reaction
from rmgpy.species import Species

################################################################################

#: Same spellings the plasma database's own families use.
ADJACENCY_LISTS = {
    'e': '1 e u1 p0 c-1',
    'Li': '1 Li u1 p0 c0',
    'Liplus': '1 Li u0 p0 c+1',
}

#: A synthetic ionisation family label, deliberately NOT shipped in
#: ``FAMILY_ELECTRON_PLACEMENT``: this ticket delivers the code path, not the
#: data. Injected per-test so the shipped registry stays closed.
IONISATION_FAMILY = 'Li_Electron_Impact_Ionization'

ATTACHMENT_FAMILY = 'Plasma_Electron_Attachment'
RECOMBINATION_FAMILY = 'Cation_R_Recombination'

#: Attachment and recombination rate coefficients as measured on the base
#: commit (98b0e70b4), before the declaration was widened. The negative control
#: compares against these literals, not against a recomputation.
BASE_ATTACHMENT_K_1000 = 1000000.0000000001
BASE_RECOMBINATION_K_1000 = 10.000000000000002


def _species():
    """A fresh label -> Species map. Fresh per test: species are indexed by
    identity, and a shape that appears twice must not share objects."""
    return {label: Species(label=label).from_adjacency_list(adjlist)
            for label, adjlist in ADJACENCY_LISTS.items()}


def _voronov_li(electrons=1):
    """Real electron-impact ionisation kinetics for Li I -> Li II, read from the
    plasma database's ``voronov.yaml`` (Z=3, N=3). Order 2 by its own units."""
    return VoronovEIArrhenius(Z=3, N=3, electrons=electrons)


def _badnell_li(electrons=-1):
    """Real radiative recombination kinetics for Li II -> Li I, from
    ``badnell.yaml`` (Z=3, N=2). Order 2 by its own units."""
    return BadnellRRArrhenius(Z=3, N=2, electrons=electrons)


def _thermal_bimolecular():
    """An ordinary Arrhenius rate law. ``get_plasma_rate_order`` returns ``None``
    for it, which is exactly the silent gap the cross-check must not read as
    agreement."""
    return Arrhenius(A=(1.0e12, 'cm^3/(mol*s)'), n=0.0, Ea=(0.0, 'kJ/mol'))


def _ionisation(declaration_electrons=1, kinetics=None, family=IONISATION_FAMILY,
                **overrides):
    """``Li => Li+`` in canonical (fully collapsed metadata) form: the electrons
    live only in the scalar. This is the shape the loader accepts and the
    reactor refuses; the placement view is what stands between them."""
    spc = _species()
    kwargs = dict(
        reactants=[spc['Li']],
        products=[spc['Liplus']],
        family=family,
        electrons=declaration_electrons,
        reversible=False,
        is_forward=True,
        kinetics=_voronov_li() if kinetics is None else kinetics,
    )
    kwargs.update(overrides)
    return TemplateReaction(**kwargs)


@pytest.fixture
def declare(monkeypatch):
    """Inject a placement declaration for the duration of one test, then remove
    it. The shipped table is a closed, hand-maintained registry."""
    def _declare(label, value):
        monkeypatch.setitem(FAMILY_ELECTRON_PLACEMENT, label, value)
    return _declare


################################################################################
# The headline: an ionisation-shaped reaction resolving through placement
################################################################################

class TestIonisationPlacementView:

    def test_real_voronov_ionisation_resolves(self, declare):
        """``Li + e- -> Li+ + 2e-`` with the real ``voronov.yaml`` Li I -> Li II
        coefficient. Before the widening this was unreachable at any family:
        ``(side, count)`` cannot say "one in, two out"."""
        declare(IONISATION_FAMILY, (1, 2))
        reaction = _ionisation()
        electron = _species()['e']

        reactants_before = list(reaction.reactants)
        products_before = list(reaction.products)

        view = resolve_electron_placement(
            reaction, [electron] + reactants_before + products_before)

        assert isinstance(view, Reaction)
        assert view.electrons == 0
        assert sum(1 for spc in view.reactants if spc.is_electron()) == 1
        assert sum(1 for spc in view.products if spc.is_electron()) == 2
        assert len(view.reactants) == 2 and len(view.products) == 3

        # The canonical reaction is untouched: same scalar, same lists.
        assert reaction.electrons == 1
        assert reaction.reactants == reactants_before
        assert reaction.products == products_before

        # The rate law is shared by reference, not copied or refitted.
        assert view.kinetics is reaction.kinetics

    def test_incident_order_is_the_declared_reactant_count(self, declare):
        """Second order overall -- one Li, one electron -- which is what makes
        the Voronov coefficient's ``cm^3/(mol*s)`` dimensionally meaningful."""
        declare(IONISATION_FAMILY, (1, 2))
        view = resolve_electron_placement(_ionisation(), [_species()['e']])
        assert len(view.reactants) == 2

    def test_same_net_count_two_declarations_two_views(self, declare):
        """The central invariant, shown rather than asserted: ``electrons = +1``
        is held fixed and only the DECLARATION moves. ``(1, 2)`` gives the
        second-order electron-impact view; ``(0, 1)`` gives the first-order
        spontaneous-ionisation view. If the scalar were the order source both
        would come out the same."""
        declare(IONISATION_FAMILY, (1, 2))
        impact = resolve_electron_placement(
            _ionisation(kinetics=_thermal_bimolecular()), [_species()['e']])

        declare(IONISATION_FAMILY, (0, 1))
        spontaneous = resolve_electron_placement(
            _ionisation(kinetics=_thermal_bimolecular()), [_species()['e']])

        assert len(impact.reactants) == 2 and len(impact.products) == 3
        assert len(spontaneous.reactants) == 1 and len(spontaneous.products) == 2
        assert sum(1 for s in impact.reactants if s.is_electron()) == 1
        assert sum(1 for s in spontaneous.reactants if s.is_electron()) == 0

    def test_scalar_that_contradicts_the_declaration_is_refused_by_name(self, declare):
        """Step 5 generalises to ``electrons == product_count - reactant_count``
        with no change of character: the scalar is still only ever validated."""
        declare(IONISATION_FAMILY, (1, 2))
        with pytest.raises(ElectronPlacementError, match='declares net') as excinfo:
            resolve_electron_placement(_ionisation(declaration_electrons=2),
                                       [_species()['e']])
        assert IONISATION_FAMILY in str(excinfo.value)

    def test_second_charge_state_is_blocked_upstream_of_placement(self):
        """``voronov.yaml`` carries Li II -> Li III as its own (Z=3, N=2) entry,
        and a per-charge-state declaration would cover it one step at a time.
        It is unreachable for a reason that has nothing to do with placement:
        RMG's atomtype table has no doubly-charged lithium, so the PRODUCT
        cannot be built. Pinned here so the next reader does not spend the
        widening looking for the wall in this module."""
        with pytest.raises(AtomTypeError, match=r'Li\+\+'):
            Species(label='Li2plus').from_adjacency_list('1 Li u0 p0 c+2')


################################################################################
# Step 10: the structural verification, unedited
################################################################################

class TestStructuralVerificationStillHolds:

    def test_ionisation_view_balances_in_the_e_pseudo_element(self, declare):
        """Measured, not assumed: under the writers' charge rule the finished
        ionisation view is E=(1, 1). Step 10 therefore passes it with no
        weakening -- it is the same check attachment passes."""
        declare(IONISATION_FAMILY, (1, 2))
        view = resolve_electron_placement(_ionisation(), [_species()['e']])
        left = sum(get_species_electron_count(s) for s in view.reactants)
        right = sum(get_species_electron_count(s) for s in view.products)
        assert (left, right) == (1, 1)

    def test_attachment_view_balances_in_the_e_pseudo_element(self):
        o2 = Species(label='O2').from_smiles('[O][O]')
        o2m = Species(label='O2-').from_smiles('[O][O-]')
        reaction = TemplateReaction(
            reactants=[o2], products=[o2m], family=ATTACHMENT_FAMILY,
            electrons=-1, reversible=False, is_forward=True,
            kinetics=_thermal_bimolecular())
        view = resolve_electron_placement(reaction, [_species()['e'], o2, o2m])
        left = sum(get_species_electron_count(s) for s in view.reactants)
        right = sum(get_species_electron_count(s) for s in view.products)
        assert (left, right) == (1, 1)

    def test_wrong_sided_two_sided_declaration_is_caught_by_step_ten(self, declare):
        """A declaration of ``(2, 1)`` on a net -1 reaction passes the scalar
        validation -- ``1 - 2 == -1`` -- and is still refused, because the
        finished view does not balance. Step 10 checks the participants
        themselves, which no producer controls, and it keeps biting once the
        declaration has two sides."""
        declare(IONISATION_FAMILY, (2, 1))
        reaction = _ionisation(declaration_electrons=-1,
                               kinetics=_thermal_bimolecular())
        with pytest.raises(ElectronPlacementError,
                           match='does not balance in the E pseudo-element'):
            resolve_electron_placement(reaction, [_species()['e']])


################################################################################
# The rate-order cross-check -- a check, never an order source
################################################################################

class TestRateOrderCrossCheck:

    def test_agreement_is_recorded_on_the_view(self, declare):
        declare(IONISATION_FAMILY, (1, 2))
        view = resolve_electron_placement(_ionisation(), [_species()['e']])
        assert RATE_ORDER_AGREES in view.comment

    def test_first_order_declaration_with_bimolecular_kinetics_is_refused(self, declare):
        """The failure the cross-check exists for. ``(0, 1)`` declares
        ``Li -> Li+ + e-`` -- first order, which would need ``s^-1`` -- while the
        Voronov coefficient is second order. Accepting it would make the rate
        wrong by a factor of the electron density with the file looking well
        formed. This is now a placement-time refusal instead of an export-time
        one."""
        declare(IONISATION_FAMILY, (0, 1))
        with pytest.raises(ElectronPlacementError,
                           match='rate coefficient of order') as excinfo:
            resolve_electron_placement(_ionisation(), [_species()['e']])
        message = str(excinfo.value)
        assert IONISATION_FAMILY in message
        assert 'VoronovEIArrhenius' in message

    def test_unresolvable_order_is_its_own_named_outcome(self, declare):
        """``get_plasma_rate_order`` returns ``None`` for any rate law outside
        its units table -- a silent gap. A silent ``None`` must not read as
        agreement, so it gets its own name on the view."""
        declare(IONISATION_FAMILY, (1, 2))
        view = resolve_electron_placement(
            _ionisation(kinetics=_thermal_bimolecular()), [_species()['e']])
        assert RATE_ORDER_UNRESOLVABLE in view.comment
        assert RATE_ORDER_AGREES not in view.comment

    def test_the_two_outcomes_are_distinguishable(self):
        assert RATE_ORDER_AGREES != RATE_ORDER_UNRESOLVABLE

    def test_cross_check_is_not_an_order_source(self, declare):
        """Order 2 kinetics with a ``(0, 1)`` declaration is REFUSED, not
        silently corrected to ``(1, 2)``. Placement never comes from the units."""
        declare(IONISATION_FAMILY, (0, 1))
        with pytest.raises(ElectronPlacementError):
            resolve_electron_placement(_ionisation(), [_species()['e']])


################################################################################
# The list-append primitive: a second one, kept in step with the export path
################################################################################

class TestPlacementPrimitive:

    @pytest.mark.parametrize('count', [1, 2, 3])
    def test_matches_expand_electrons_on_consumption(self, count):
        """``expand_electrons`` is net-derived and shared with the export path,
        so it cannot express a two-sided placement. The placement primitive is a
        separate function -- and on every shape ``expand_electrons`` CAN express,
        the two agree exactly. That equality is what keeps the reactor boundary
        and the export boundary from drifting apart."""
        spc = _species()
        reaction = Reaction(reactants=[spc['Li']], products=[spc['Liplus']],
                            electrons=-count)
        expected = expand_electrons(reaction, [spc['e']])
        got = _place_declared_electrons(reaction, spc['e'], count, 0)
        assert [str(s) for s in got[0]] == [str(s) for s in expected[0]]
        assert [str(s) for s in got[1]] == [str(s) for s in expected[1]]

    @pytest.mark.parametrize('count', [1, 2, 3])
    def test_matches_expand_electrons_on_production(self, count):
        spc = _species()
        reaction = Reaction(reactants=[spc['Li']], products=[spc['Liplus']],
                            electrons=count)
        expected = expand_electrons(reaction, [spc['e']])
        got = _place_declared_electrons(reaction, spc['e'], 0, count)
        assert [str(s) for s in got[0]] == [str(s) for s in expected[0]]
        assert [str(s) for s in got[1]] == [str(s) for s in expected[1]]

    def test_does_not_mutate_the_reaction(self):
        spc = _species()
        reaction = Reaction(reactants=[spc['Li']], products=[spc['Liplus']],
                            electrons=1)
        reactants_before, products_before = reaction.reactants, reaction.products
        _place_declared_electrons(reaction, spc['e'], 1, 2)
        assert reaction.reactants is reactants_before
        assert reaction.products is products_before
        assert len(reaction.reactants) == 1 and len(reaction.products) == 1

    def test_places_on_both_sides(self):
        spc = _species()
        reaction = Reaction(reactants=[spc['Li']], products=[spc['Liplus']],
                            electrons=1)
        reactants, products = _place_declared_electrons(reaction, spc['e'], 1, 2)
        assert sum(1 for s in reactants if s.is_electron()) == 1
        assert sum(1 for s in products if s.is_electron()) == 2


################################################################################
# The declaration schema itself
################################################################################

class TestDeclarationSchema:

    def test_shipped_registry_is_still_closed_and_two_sided(self):
        """The table is a closed, hand-maintained list. Two families and, since
        I-119, one kinetics library at ``(1, 0)`` -- one electron in, none out --
        plus, since I-116, the ``PlasmaElectronImpactIonization`` library at
        ``(1, 2)``, plus, since I-154, three families at ``(0, 1)`` -- none in,
        one out. Widening the schema in I-113 must not have added chemistry to
        the table by itself, and it did not: every later entry arrived by a
        decision, and any EIGHTH one still fails here.

        Note what the widened schema bought and what it did not. Only
        ``PlasmaElectronImpactIonization`` needs two differing NON-ZERO numbers;
        the recombination entry I-119 added is attachment-shaped and would have
        fitted the old ``(side, count)`` spelling, and so would the three I-154
        added, as ``('products', 1)``. The schema is wider than six of its
        seven users need, which is correct -- the one that needs it could not be
        expressed at all before."""
        assert FAMILY_ELECTRON_PLACEMENT == {
            'Plasma_Electron_Attachment': (1, 0),
            'Cation_R_Recombination': (1, 0),
            'PlasmaRadiativeRecombination': (1, 0),
            'PlasmaElectronImpactIonization': (1, 2),
            'Plasma_Associative_Ionization_Alkali_Alkali': (0, 1),
            'Plasma_Associative_Ionization_Alkali_Alkaline': (0, 1),
            'Plasma_Associative_Ionization_Alkaline_Alkaline': (0, 1),
        }

    def test_no_electron_impact_ionisation_family_is_shipped(self):
        """I-113 delivered the code path and left the owner open. I-116 answered
        it with a kinetics LIBRARY, ``PlasmaElectronImpactIonization``, and I-154
        confirmed that ruling by declining to carry the earlier development
        line's ``Plasma_Electron_Impact_Ionization`` FAMILY, which would have
        been a second, disagreeing owner for the same chemistry.

        Read this test's name precisely, because I-154 DID ship ionisation
        families -- three of them, at ``(0, 1)``. Those are ASSOCIATIVE
        ionisation, where no electron is incident.
        What is still unshipped as a family is ELECTRON-IMPACT ionisation,
        ``A + e- -> A+ + 2 e-``, which is this file's subject -- so its synthetic
        family label is still undeclared and every test here still injects its
        own declaration. The test was renamed for that distinction; the old name
        asserted something that has become false."""
        assert IONISATION_FAMILY not in FAMILY_ELECTRON_PLACEMENT
        assert not any(key.startswith('Plasma_Electron_Impact_Ionization')
                       for key in FAMILY_ELECTRON_PLACEMENT)

    @pytest.mark.parametrize('bad', [
        ('reactants', 1),   # the old two-field shape, now meaningless
        (1,),               # too short
        (1, 0, 0),          # too long
        (-1, 2),            # negative count
        (0, 0),             # places no electron at all
        (1.5, 0),           # not integral
        None,
    ])
    def test_malformed_declaration_is_refused_by_name(self, declare, bad):
        declare(IONISATION_FAMILY, bad)
        with pytest.raises(ElectronPlacementError, match='declaration'):
            resolve_electron_placement(_ionisation(), [_species()['e']])


################################################################################
# Kinetics-level electron counts
################################################################################

class TestKineticsLevelElectronCount:

    def test_plasma_rate_law_net_count_is_validated_against_the_declaration(self, declare):
        """``VoronovEIArrhenius.electrons`` and ``BadnellRRArrhenius.electrons``
        are documented to be the SAME signed net quantity as
        ``Reaction.electrons``. They are therefore validated against the
        declaration, exactly as the reaction scalar is -- a third source that
        must agree, not a second one that must be absent."""
        declare(IONISATION_FAMILY, (1, 2))
        view = resolve_electron_placement(_ionisation(kinetics=_voronov_li(1)),
                                          [_species()['e']])
        assert view.electrons == 0

    def test_plasma_rate_law_net_count_that_disagrees_is_refused(self, declare):
        declare(IONISATION_FAMILY, (1, 2))
        with pytest.raises(ElectronPlacementError,
                           match='disagrees with the family declaration'):
            resolve_electron_placement(_ionisation(kinetics=_voronov_li(2)),
                                       [_species()['e']])

    def test_badnell_recombination_net_count_is_validated(self):
        """The mirror case on the shipped ``(1, 0)`` declaration, with real
        ``badnell.yaml`` radiative-recombination kinetics."""
        spc = _species()
        reaction = TemplateReaction(
            reactants=[spc['Liplus']], products=[spc['Li']],
            family=RECOMBINATION_FAMILY, electrons=-1, reversible=False,
            is_forward=True, kinetics=_badnell_li(-1))
        view = resolve_electron_placement(reaction, [spc['e']])
        assert sum(1 for s in view.reactants if s.is_electron()) == 1
        assert len(view.reactants) == 2          # order 2, matching badnell.yaml
        assert RATE_ORDER_AGREES in view.comment

    def test_other_kinetics_electron_field_is_still_refused(self):
        """Charge-transfer rate laws carry ``electrons`` meaning "electrons
        transferred", which is not the same quantity. Only the two plasma laws
        that document the net convention are reconciled; everything else stays a
        second placement source and stays fatal."""
        class DeclaredElectronKinetics(Arrhenius):
            electrons = -1

        o2 = Species(label='O2').from_smiles('[O][O]')
        o2m = Species(label='O2-').from_smiles('[O][O-]')
        reaction = TemplateReaction(
            reactants=[o2], products=[o2m], family=ATTACHMENT_FAMILY,
            electrons=-1, reversible=False, is_forward=True,
            kinetics=DeclaredElectronKinetics(A=(1.0e12, 'cm^3/(mol*s)'), n=0.0,
                                              Ea=(0.0, 'kJ/mol')))
        with pytest.raises(ElectronPlacementError, match='its own electron count'):
            resolve_electron_placement(reaction, [_species()['e'], o2, o2m])


################################################################################
# Negative control: the shapes that already worked must not have moved
################################################################################

class TestDeclaredFamiliesNumericallyUnchanged:

    def test_attachment_view_is_unchanged(self):
        """Participants, order and rate coefficient identical to the values
        measured on the base commit before the widening."""
        o2 = Species(label='O2').from_smiles('[O][O]')
        o2m = Species(label='O2-').from_smiles('[O][O-]')
        reaction = TemplateReaction(
            reactants=[o2], products=[o2m], family=ATTACHMENT_FAMILY,
            electrons=-1, reversible=False, is_forward=True,
            kinetics=_thermal_bimolecular())
        view = resolve_electron_placement(reaction, [_species()['e'], o2, o2m])

        assert [str(s) for s in view.reactants] == ['O2', 'e']
        assert [str(s) for s in view.products] == ['O2-']
        assert view.electrons == 0
        assert view.kinetics.get_rate_coefficient(1000.0) == pytest.approx(
            BASE_ATTACHMENT_K_1000, rel=1e-12)

    def test_recombination_view_is_unchanged(self):
        spc = _species()
        r = Species(label='CH3').from_smiles('[CH3]')
        rli = Species(label='LiCH3').from_smiles('[Li]C')
        reaction = TemplateReaction(
            reactants=[spc['Liplus'], r], products=[rli],
            family=RECOMBINATION_FAMILY, electrons=-1, reversible=False,
            is_forward=True,
            kinetics=TwoTemperaturePlasma(A=(1.0e4, 'm^6/(mol^2*s)'), n=-1.0,
                                          Ea_g=(0.0, 'J/mol'),
                                          Ea_e=(0.0, 'J/mol')))
        view = resolve_electron_placement(reaction, [spc['e'], spc['Liplus'], r, rli])

        assert [str(s) for s in view.reactants] == ['Liplus', 'CH3', 'e']
        assert [str(s) for s in view.products] == ['LiCH3']
        assert view.electrons == 0
        assert view.kinetics.get_rate_coefficient(1000.0) == pytest.approx(
            BASE_RECOMBINATION_K_1000, rel=1e-12)
