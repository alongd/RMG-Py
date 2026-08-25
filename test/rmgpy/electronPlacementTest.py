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
Tests for :mod:`rmgpy.electron_placement` — the electron-representation
boundary between the database representation of electron attachment (electron
carried as ``Reaction.electrons`` metadata) and the reactor representation
(electron as an explicit reactant).

The resolver under test must place the electron from the FAMILY's placement
declaration, never from the net scalar; must never mutate the canonical
reaction; and must fail BY NAME for everything outside single-electron
attachment in the forward direction of ``Plasma_Electron_Attachment``.
"""

import os.path

import numpy as np
import pytest

import rmgpy.electron_placement as electron_placement
from rmgpy import settings
from rmgpy.data.kinetics.common import get_molecularity
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.electron_balance import check_electron_balance
from rmgpy.electron_placement import FAMILY_ELECTRON_PLACEMENT, resolve_electron_placement
from rmgpy.exceptions import (
    ElectronPlacementError,
    NonEquilibriumReverseRateError,
    PlasmaStateError,
)
from rmgpy.kinetics import Arrhenius
from rmgpy.kinetics.arrhenius import TwoTemperaturePlasma
from rmgpy.reaction import Reaction
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.species import Species

FAMILY = "Plasma_Electron_Attachment"
# An electron-bearing family that also declares electrons = -1 but carries NO
# placement declaration; proves the resolver did not silently generalize to
# every family that happens to consume an electron. (Until I-088 this was
# Cation_R_Recombination, which is now declared -- a surface electrochemistry
# family is used instead, deliberately far from the plasma families.)
OTHER_ELECTRON_FAMILY = "Surface_Proton_Electron_Reduction_Alpha"

NEUTRAL_O2 = "[O][O]"
ANION_O2 = "[O][O-]"

T_GAS = 500.0     # K
T_E = 11604.5     # K (~1 eV)
P0 = 1.0e5        # Pa
Y_E0 = 1.0e-4     # mol, seeded electrons


def _electron():
    return Species(label="e-").from_adjacency_list("1 e u1 p0 c-1")


def _o2():
    return Species(label="O2").from_smiles(NEUTRAL_O2)


def _o2_anion():
    return Species(label="O2-").from_smiles(ANION_O2)


def _arrhenius():
    return Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kJ/mol"))


def _attachment(**overrides):
    """A synthetic canonical attachment reaction, family-attributed, with the
    electron carried only as ``electrons = -1`` metadata."""
    kwargs = dict(
        reactants=[_o2()],
        products=[_o2_anion()],
        family=FAMILY,
        electrons=-1,
        reversible=False,
        is_forward=True,
        kinetics=_arrhenius(),
    )
    kwargs.update(overrides)
    return TemplateReaction(**kwargs)


class TestElectronPlacementResolver:
    """Unit tests: the resolver's acceptance conditions and named refusals,
    on synthetic reactions (no database required)."""

    def test_synthetic_attachment_resolves_nonmutating(self):
        """Acceptance: the view carries the electron explicitly with
        electrons=0, while the canonical reaction is untouched."""
        reaction = _attachment()
        electron = _electron()
        species_list = [electron, reaction.reactants[0], reaction.products[0]]

        reactants_before = list(reaction.reactants)
        products_before = list(reaction.products)
        reactant_list_object = reaction.reactants
        product_list_object = reaction.products

        view = resolve_electron_placement(reaction, species_list)

        # The view: explicit electron on the reactant side, zero metadata count.
        assert isinstance(view, Reaction)
        assert view.electrons == 0
        assert sum(1 for spc in view.reactants if spc.is_electron()) == 1
        assert any(spc is electron for spc in view.reactants)
        assert len(view.reactants) == len(reactants_before) + 1
        assert not any(spc.is_electron() for spc in view.products)
        # Copied lists, shared species and kinetics objects.
        assert view.reactants is not reaction.reactants
        assert view.products is not reaction.products
        assert view.reactants[0] is reactants_before[0]
        assert view.products[0] is products_before[0]
        assert view.kinetics is reaction.kinetics
        assert view.reversible is False

        # The canonical reaction, inspected directly after resolution: same
        # metadata count, same participant list objects, same contents, no
        # explicit electron.
        assert reaction.electrons == -1
        assert reaction.reactants is reactant_list_object
        assert reaction.products is product_list_object
        assert reaction.reactants == reactants_before
        assert reaction.products == products_before
        assert not any(spc.is_electron() for spc in reaction.reactants)
        assert not any(spc.is_electron() for spc in reaction.products)

    def test_molecularity_reconciliation(self):
        """get_molecularity's net-derived path and the placement view agree by
        construction: the canonical form counts len(reactants)=1 plus one
        metadata electron; the view counts the electron as an ordinary
        reactant with zero metadata. Both are bimolecular."""
        reaction = _attachment()
        electron = _electron()
        view = resolve_electron_placement(
            reaction, [electron, reaction.reactants[0], reaction.products[0]])
        assert get_molecularity(reaction) == 2
        assert get_molecularity(view) == 2

    def test_non_attachment_electron_family_fails_by_name(self):
        """A different family that also declares electrons = -1 must be
        refused by name: this increment declares placement for attachment
        only, and there is no general mechanism to inherit."""
        reaction = _attachment(family=OTHER_ELECTRON_FAMILY)
        with pytest.raises(ElectronPlacementError,
                           match="no electron-placement declaration") as excinfo:
            resolve_electron_placement(reaction, [_electron()])
        assert OTHER_ELECTRON_FAMILY in str(excinfo.value)

    def test_family_free_reaction_fails_by_name(self):
        """A plain Reaction carries no family attribution, so no declaration
        can be named; the net count must not be used instead."""
        reaction = Reaction(reactants=[_o2()], products=[_o2_anion()],
                            electrons=-1, reversible=False, kinetics=_arrhenius())
        with pytest.raises(ElectronPlacementError, match="no family attribution"):
            resolve_electron_placement(reaction, [_electron()])

    def test_ionization_shaped_fails_by_name(self):
        """Net +1 (electron production) handed to the attachment resolver is a
        named refusal, not a product-side placement."""
        reaction = _attachment(electrons=1)
        with pytest.raises(ElectronPlacementError, match="ionization-shaped") as excinfo:
            resolve_electron_placement(reaction, [_electron()])
        assert FAMILY in str(excinfo.value)

    def test_excitation_shaped_fails_by_name(self):
        """Net 0 (no net electron change) is a named refusal, not a no-op."""
        reaction = _attachment(electrons=0)
        with pytest.raises(ElectronPlacementError, match="excitation-shaped"):
            resolve_electron_placement(reaction, [_electron()])

    def test_double_representation_fails_by_name(self):
        """An explicit electron participant AND a nonzero metadata count would
        double-count electron stoichiometry or rate order; fatal by name."""
        electron = _electron()
        reaction = _attachment(reactants=[_o2(), electron])
        with pytest.raises(ElectronPlacementError,
                           match="represents the electron twice"):
            resolve_electron_placement(reaction, [electron])

    def test_already_explicit_fails_by_name(self):
        """An explicit electron with electrons=0 is already in reactor form;
        placing another would double it."""
        electron = _electron()
        reaction = _attachment(reactants=[_o2(), electron], electrons=0)
        with pytest.raises(ElectronPlacementError, match="already carries"):
            resolve_electron_placement(reaction, [electron])

    def test_reverse_generated_reaction_resolves(self):
        """I-088 replaces the old ``is_forward=False`` refusal.

        ``is_forward`` records how a reaction was FOUND, not how its participant
        lists are ordered: ``KineticsFamily._create_reaction`` stores
        family-forward molecular orientation in BOTH generation directions. So a
        reverse-generated reaction's electron is still a reactant, and refusing
        on ``is_forward`` excluded exactly the chemistry the plasma decks
        generate. It resolves, to the same view a forward-found one gives."""
        reaction = _attachment(is_forward=False)
        view = resolve_electron_placement(reaction, [_electron()])
        assert sum(1 for spc in view.reactants if spc.is_electron()) == 1
        assert not any(spc.is_electron() for spc in view.products)
        assert view.electrons == 0

    def test_unknown_direction_resolves(self):
        """Same reasoning for ``is_forward=None``. The old BLOCKER-1 check
        refused a direction-unknown reaction on the grounds that it "would
        silently manufacture a forward-direction view from ambiguous input" --
        but the view is not manufactured from the direction at all, and since
        I-088 the side is verified structurally instead (see
        ``test_wrong_side_placement_refused_by_E_balance``, which is what
        actually catches a bad orientation and which ``is_forward`` never
        could).

        A duck-typed stand-in is used because the compiled ``Reaction`` stores
        ``is_forward`` as a C ``bint`` that coerces ``None`` to ``False``;
        ``resolve_electron_placement`` is pure Python and sees the real
        ``None``."""

        class _UnknownDirectionReaction:
            family = FAMILY
            electrons = -1
            reversible = False
            is_forward = None  # never established -- pure-Python default
            reactants = [_o2()]
            products = [_o2_anion()]
            kinetics = _arrhenius()
            label = ''
            duplicate = False
            degeneracy = 1

            def __str__(self):
                return 'O2 + e- -> O2- (direction unknown)'

        view = resolve_electron_placement(_UnknownDirectionReaction(), [_electron()])
        assert sum(1 for spc in view.reactants if spc.is_electron()) == 1
        assert not any(spc.is_electron() for spc in view.products)

    def test_wrong_side_placement_refused_by_E_balance(self):
        """The check that REPLACED the ``is_forward`` refusals, and is strictly
        stronger than either of them.

        This reaction is genuinely reverse-ORIENTED -- its reactant side is the
        family's product side (the anion) -- while its metadata still claims the
        family-forward ``electrons = -1``. That is the I-086 defect class: the
        producer's orientation and its electron sign disagree. Every metadata
        check passes (declared family, correct net count, no explicit electron),
        and ``is_forward=True`` here, so NEITHER old direction check would have
        fired. Placing the electron on the declared reactant side gives E=2 on
        the left against E=0 on the right, and the structural check refuses it."""
        reaction = _attachment(reactants=[_o2_anion()], products=[_o2()],
                               is_forward=True)
        with pytest.raises(ElectronPlacementError,
                           match="does not balance in the E pseudo-element"):
            resolve_electron_placement(reaction, [_electron()])

    def test_reversible_resolves_and_stays_reversible(self):
        """I-088 replaces the old ``reversible=True`` refusal.

        Reversibility does not make the electron's side ambiguous -- the stored
        orientation still has a definite reactant side. What it makes ambiguous
        is the reverse RATE, and that is the consumer's policy: the resolver
        passes ``reversible`` through untouched so ``PlasmaReactor`` can refuse
        it by name (locked by
        ``TestCationRecombinationPlacement.test_05_...`` against the real
        reactor)."""
        reaction = _attachment(reversible=True)
        view = resolve_electron_placement(reaction, [_electron()])
        assert sum(1 for spc in view.reactants if spc.is_electron()) == 1
        assert view.reversible is True

    def test_missing_kinetics_fails_by_name(self):
        reaction = _attachment(kinetics=None)
        with pytest.raises(ElectronPlacementError, match="has no kinetics"):
            resolve_electron_placement(reaction, [_electron()])

    def test_pressure_dependent_kinetics_fails_by_name(self):
        class DeclaredPdepKinetics(Arrhenius):
            def is_pressure_dependent(self):
                return True

        reaction = _attachment(
            kinetics=DeclaredPdepKinetics(A=(1.0e12, "cm^3/(mol*s)"), n=0.0,
                                          Ea=(0.0, "kJ/mol")))
        with pytest.raises(ElectronPlacementError, match="pressure-dependent"):
            resolve_electron_placement(reaction, [_electron()])

    def test_kinetics_level_electron_field_fails_by_name(self):
        """Charge-transfer-style kinetics carry their own electron count; that
        is a second placement source and must be refused, not reconciled."""
        class DeclaredElectronKinetics(Arrhenius):
            electrons = -1

        reaction = _attachment(
            kinetics=DeclaredElectronKinetics(A=(1.0e12, "cm^3/(mol*s)"), n=0.0,
                                              Ea=(0.0, "kJ/mol")))
        with pytest.raises(ElectronPlacementError, match="its own electron count"):
            resolve_electron_placement(reaction, [_electron()])

    def test_missing_electron_species_fails_by_name(self):
        reaction = _attachment()
        with pytest.raises(ElectronPlacementError,
                           match="No electron species is resolvable"):
            resolve_electron_placement(
                reaction, [reaction.reactants[0], reaction.products[0]])

    def test_ambiguous_electron_species_fails_by_name(self):
        reaction = _attachment()
        with pytest.raises(ElectronPlacementError, match="must be unique"):
            resolve_electron_placement(reaction, [_electron(), _electron()])

    def test_declaration_registry_is_explicit_and_closed(self):
        """The registry is a closed, hand-maintained list. Two families are
        declared, both in the single-electron-on-the-reactant-side shape --
        spelled ``(reactant_count, product_count)`` since I-113 widened the
        declaration; any other entry appearing here means placement semantics
        generalized without someone deciding they should."""
        assert FAMILY_ELECTRON_PLACEMENT == {
            "Plasma_Electron_Attachment": (1, 0),
            "Cation_R_Recombination": (1, 0),
        }


@pytest.mark.database
class TestAttachmentPlacementDatabase:
    """The resolver against REAL family generation (not a library lookup),
    and the reactor-facing half of the boundary: the reactor rejects the
    canonical form and accepts the resolved view."""

    @classmethod
    def setup_class(cls):
        cls.families_path = os.path.join(settings["database.directory"],
                                         "kinetics", "families")
        database = KineticsDatabase()
        database.load_recommended_families(
            os.path.join(cls.families_path, "recommended.py"))
        database.load_families(cls.families_path, families=[FAMILY])

        species = Species().from_smiles(NEUTRAL_O2)
        species.generate_resonance_structures()
        reactions = database.generate_reactions_from_families(
            [species], only_families=[FAMILY])
        assert len(reactions) == 1
        cls.reaction = reactions[0]
        assert isinstance(cls.reaction, TemplateReaction)
        assert cls.reaction.family == FAMILY

        # Give the generated reaction its kinetics from the family's own rate
        # rules -- the same step the model builder performs. A fresh family is
        # loaded because add_rules_from_training/fill_rules_by_averaging_up
        # mutate the family's rule table.
        kinetics_database = KineticsDatabase()
        kinetics_database.load_families(cls.families_path, families=[FAMILY])
        family = kinetics_database.families[FAMILY]
        family.add_rules_from_training(thermo_database=None)
        family.fill_rules_by_averaging_up()
        template = family.retrieve_template(cls.reaction.template)
        cls.reaction.kinetics = family.get_kinetics_for_template(
            template, degeneracy=cls.reaction.degeneracy)[0]

        cls.electron = _electron()
        cls.species_list = [cls.electron] + list(cls.reaction.reactants) \
            + list(cls.reaction.products)

    def _resolve(self):
        return resolve_electron_placement(self.reaction, self.species_list)

    def _reactor(self):
        imf = {self.electron: Y_E0, self.reaction.reactants[0]: 1.0}
        return PlasmaReactor(T_GAS, P0, imf, (T_E, "K"), n_sims=1, termination=[])

    def _core_species(self):
        return [self.electron, self.reaction.reactants[0], self.reaction.products[0]]

    def test_generated_reaction_resolves_nondestructively(self):
        """The family-generated O2 + e- -> O2- reaction resolves to a view
        with the electron among the reactants and electrons=0, while the
        canonical reaction object -- inspected directly afterwards -- still
        carries electrons=-1 and its original participant lists."""
        reactants_before = list(self.reaction.reactants)
        products_before = list(self.reaction.products)

        view = self._resolve()

        assert view.electrons == 0
        assert sum(1 for spc in view.reactants if spc.is_electron()) == 1
        assert any(spc is self.electron for spc in view.reactants)
        assert len(view.reactants) == 2
        assert view.reactants[0] is reactants_before[0]
        assert view.products == products_before
        assert view.kinetics is self.reaction.kinetics

        # The non-destructive claim, asserted directly on the canonical object.
        assert self.reaction.electrons == -1
        assert self.reaction.reactants == reactants_before
        assert self.reaction.products == products_before
        assert not any(spc.is_electron() for spc in self.reaction.reactants)
        assert not any(spc.is_electron() for spc in self.reaction.products)

    def test_molecularity_reconciliation_on_generated_reaction(self):
        """Both representations of the generated reaction are bimolecular:
        net-derived on the canonical form, participant-derived on the view."""
        view = self._resolve()
        assert get_molecularity(self.reaction) == 2
        assert get_molecularity(view) == 2

    def test_reactor_autoresolves_canonical_form(self):
        """After the electron-representation boundary was wired into
        PlasmaReactor.initialize_model, the SAME reactor that once refused the
        metadata-only representation now accepts it automatically: the reactor
        invokes the resolver at its generation-to-reactor handoff, so the
        canonical reaction reaches the solver as a resolved view with a
        positive forward rate. (Previously this raised PlasmaStateError
        "metadata-only electron count"; that rejection now lives only on the
        resolver's contract, exercised for shapes the family cannot place.)"""
        reactor = self._reactor()
        reactor.initialize_model(self._core_species(), [self.reaction], [], [])
        assert reactor.kf[0] > 0.0
        # And the automatic resolution left the canonical object untouched.
        assert self.reaction.electrons == -1
        assert not any(spc.is_electron() for spc in self.reaction.reactants)

    def test_reactor_accepts_resolved_view(self):
        """The same reactor configuration that rejects the canonical form
        accepts the resolved view, and gives it a positive forward rate."""
        view = self._resolve()
        reactor = self._reactor()
        reactor.initialize_model(self._core_species(), [view], [], [])
        assert reactor.kf[0] > 0.0
        # And resolution for the reactor still left the canonical object alone.
        assert self.reaction.electrons == -1
        assert len(self.reaction.reactants) == 1


def _load_family_and_generate_attachment():
    """Load the real ``Plasma_Electron_Attachment`` family from the database and
    generate the ``O2 + e- -> O2-`` reaction from it (NOT a library lookup),
    with kinetics filled in from the family's own rate rules -- the same steps
    the model builder performs. Returns the canonical :class:`TemplateReaction`
    (electron carried only as ``electrons = -1`` metadata). Loading a fresh
    database each call is deliberate: it lets a test prove the reactor path is
    exercised anew, never from a cached or special-cased result."""
    families_path = os.path.join(settings["database.directory"], "kinetics", "families")
    database = KineticsDatabase()
    database.load_recommended_families(os.path.join(families_path, "recommended.py"))
    database.load_families(families_path, families=[FAMILY])

    species = Species().from_smiles(NEUTRAL_O2)
    species.generate_resonance_structures()
    reactions = database.generate_reactions_from_families([species], only_families=[FAMILY])
    assert len(reactions) == 1
    reaction = reactions[0]
    assert isinstance(reaction, TemplateReaction)
    assert reaction.family == FAMILY
    assert reaction.electrons == -1
    assert not any(spc.is_electron() for spc in reaction.reactants)

    kinetics_database = KineticsDatabase()
    kinetics_database.load_families(families_path, families=[FAMILY])
    family = kinetics_database.families[FAMILY]
    family.add_rules_from_training(thermo_database=None)
    family.fill_rules_by_averaging_up()
    template = family.retrieve_template(reaction.template)
    reaction.kinetics = family.get_kinetics_for_template(
        template, degeneracy=reaction.degeneracy)[0]
    return reaction


def _net_charge(spc):
    """Net electrical charge of a reactor species; the electron is -1."""
    if spc.is_electron():
        return -1
    return spc.molecule[0].get_net_charge()


class TestElectronPlacementReactorIntegration:
    """The electron-representation boundary, end to end and AUTOMATIC.

    A real family-generated attachment reaction (electron carried only as
    ``electrons = -1`` metadata) is handed, through the ordinary
    direct-construction reactor handoff, to :class:`PlasmaReactor`. The reactor
    must resolve the electron placement ITSELF, at its single production call
    site, without the test ever invoking :func:`resolve_electron_placement`.

    None of the tests in this class call the resolver directly; the resolver is
    reached only through ``PlasmaReactor.initialize_model``.
    """

    @classmethod
    def setup_class(cls):
        cls.reaction = _load_family_and_generate_attachment()
        cls.electron = _electron()
        cls.o2 = cls.reaction.reactants[0]
        cls.o2_anion = cls.reaction.products[0]

    def _core_species(self):
        return [self.electron, self.o2, self.o2_anion]

    def _reactor(self):
        imf = {self.electron: Y_E0, self.o2: 1.0}
        return PlasmaReactor(T_GAS, P0, imf, (T_E, "K"), n_sims=1, termination=[])

    def _initialized_reactor(self):
        """Build a fresh reactor and run the FULL production pipeline: the
        canonical reaction goes in, automatic resolution happens inside
        initialize_model. The test never touches the resolver."""
        reactor = self._reactor()
        reactor.initialize_model(self._core_species(), [self.reaction], [], [])
        return reactor

    def _packed_reaction(self, reactor):
        """The single reaction object the reactor exposes through reaction_index.
        After the identity fix this is the CANONICAL (model) reaction, not the
        internal resolved view -- callers must be able to bridge reaction_index
        back to model reactions. The view's electron placement is inspected
        instead through the solver's packed index arrays (see _electron_reactant_
        count)."""
        reactions = list(reactor.reaction_index.keys())
        assert len(reactions) == 1
        return reactions[0]

    def _electron_reactant_count(self, reactor, reaction):
        """How many times the electron index appears on the reactant side of the
        packed (view-derived) row for `reaction` -- i.e. the incident-electron
        rate order the solver will actually evaluate."""
        j = reactor.reaction_index[reaction]
        ie = reactor.electron_index
        return sum(1 for m in range(reactor.reactant_indices.shape[1])
                   if reactor.reactant_indices[j, m] == ie)

    def _electron_product_count(self, reactor, reaction):
        j = reactor.reaction_index[reaction]
        ie = reactor.electron_index
        return sum(1 for m in range(reactor.product_indices.shape[1])
                   if reactor.product_indices[j, m] == ie)

    # --- Verifier assertions 1-16 -------------------------------------------

    def test_01_family_generated_reaction_reaches_reactor_automatically(self):
        """1. The real family-generated attachment reaction reaches the reactor
        automatically (no manual resolution, no direct resolver call)."""
        reactor = self._initialized_reactor()
        assert reactor.kf[0] > 0.0

    def test_02_canonical_reaction_unchanged_after_pipeline(self):
        """2. The canonical stored/generated reaction is unchanged after the
        full pipeline runs."""
        reactants_before = list(self.reaction.reactants)
        products_before = list(self.reaction.products)
        self._initialized_reactor()
        assert self.reaction.electrons == -1
        assert self.reaction.reactants == reactants_before
        assert self.reaction.products == products_before
        assert not any(spc.is_electron() for spc in self.reaction.reactants)
        assert not any(spc.is_electron() for spc in self.reaction.products)

    def test_03_view_has_exactly_one_incident_electron_on_reactant_side(self):
        """3. The reactor-facing view contains exactly one incident electron on
        the forward (reactant) side and none on the product side. The view is
        not exposed publicly (reaction_index maps the canonical reaction), so
        its placement is read from the solver's packed index arrays, which are
        built from the view."""
        reactor = self._initialized_reactor()
        assert self._electron_reactant_count(reactor, self.reaction) == 1
        assert self._electron_product_count(reactor, self.reaction) == 0

    def test_04_electron_enters_forward_rate_exactly_once(self):
        """4. Electron concentration enters the forward rate expression exactly
        once: the packed reactant-index row references the electron once."""
        reactor = self._initialized_reactor()
        j = reactor.reaction_index[self._packed_reaction(reactor)]
        ie = reactor.electron_index
        occurrences = sum(1 for m in range(reactor.reactant_indices.shape[1])
                          if reactor.reactant_indices[j, m] == ie)
        assert occurrences == 1

    def test_05_first_order_scaling_in_electron_concentration(self):
        """5. Changing the electron concentration produces first-order scaling
        of this attachment reaction's rate (rate = k*y_e*y_O2/V)."""
        rates = {}
        volumes = {}
        for y_e in (Y_E0, 2.0 * Y_E0):
            reactor = self._reactor()
            reactor.initial_mole_fractions = {self.electron: y_e, self.o2: 1.0}
            reactor.initialize_model(self._core_species(), [self.reaction], [], [])
            ncs = reactor.num_core_species
            y0 = np.array(reactor.y0[:ncs], float)
            res, _ = reactor.residual(0.0, y0.copy(), np.zeros(ncs, float))
            ie = reactor.electron_index
            rates[y_e] = -res[ie]  # electron consumption rate, mol/s
            volumes[y_e] = reactor.V
        expected_ratio = 2.0 * volumes[Y_E0] / volumes[2.0 * Y_E0]
        actual_ratio = rates[2.0 * Y_E0] / rates[Y_E0]
        assert abs(actual_ratio - expected_ratio) / expected_ratio < 1e-10

    def test_06_electron_derivative_negative(self):
        """6. The electron derivative dN_e/dt is negative (electron consumed)."""
        reactor = self._initialized_reactor()
        ncs = reactor.num_core_species
        y0 = np.array(reactor.y0[:ncs], float)
        res, _ = reactor.residual(0.0, y0.copy(), np.zeros(ncs, float))
        assert res[reactor.electron_index] < 0.0

    def test_07_o2_derivative_negative(self):
        """7. The O2 derivative is negative (reactant consumed)."""
        reactor = self._initialized_reactor()
        ncs = reactor.num_core_species
        y0 = np.array(reactor.y0[:ncs], float)
        res, _ = reactor.residual(0.0, y0.copy(), np.zeros(ncs, float))
        assert res[reactor.get_species_index(self.o2)] < 0.0

    def test_08_o2_anion_derivative_positive(self):
        """8. The O2- derivative is positive (product formed)."""
        reactor = self._initialized_reactor()
        ncs = reactor.num_core_species
        y0 = np.array(reactor.y0[:ncs], float)
        res, _ = reactor.residual(0.0, y0.copy(), np.zeros(ncs, float))
        assert res[reactor.get_species_index(self.o2_anion)] > 0.0

    def test_09_net_charge_conserved(self):
        """9. Net electrical charge is conserved across the reaction:
        sum_i charge_i * dN_i/dt = 0."""
        reactor = self._initialized_reactor()
        ncs = reactor.num_core_species
        y0 = np.array(reactor.y0[:ncs], float)
        res, _ = reactor.residual(0.0, y0.copy(), np.zeros(ncs, float))
        charge_rate = 0.0
        for spc in self._core_species():
            charge_rate += _net_charge(spc) * res[reactor.get_species_index(spc)]
        scale = max(abs(res[reactor.electron_index]), 1e-300)
        assert abs(charge_rate) < 1e-10 * scale

    def test_10_electron_not_double_counted(self):
        """10. The electron is not counted twice in molecularity, rate order, or
        stoichiometry: the bimolecular view has one electron reactant, and
        doubling C_e doubles (not quadruples) the forward rate."""
        reactor = self._initialized_reactor()
        packed = self._packed_reaction(reactor)
        assert get_molecularity(packed) == 2
        j = reactor.reaction_index[packed]
        ie = reactor.electron_index
        ncs = reactor.num_core_species
        y0 = np.array(reactor.y0[:ncs], float)
        reactor.residual(0.0, y0.copy(), np.zeros(ncs, float))
        conc = np.array(reactor.core_species_concentrations, float)
        kf = float(reactor.kf[j])

        def forward_rate(c_e):
            cx = conc.copy()
            cx[ie] = c_e
            rate = kf
            for m in range(reactor.reactant_indices.shape[1]):
                s = reactor.reactant_indices[j, m]
                if s != -1:
                    rate *= cx[s]
            return rate

        r1 = forward_rate(conc[ie])
        r2 = forward_rate(2.0 * conc[ie])
        assert abs(r2 / r1 - 2.0) < 1e-12  # first order: 2x, not 4x

    def test_11_residual_and_jacobian_share_coefficients(self):
        """11. Residual and Jacobian consume the same resolved coefficients and
        the same EOS from the same state -- not two independently derived
        sets."""
        reactor = self._initialized_reactor()
        ncs = reactor.num_core_species
        y0 = np.array(reactor.y0[:ncs], float)
        zeros = np.zeros(ncs, float)
        reactor.residual(0.0, y0.copy(), zeros.copy())
        v_residual = reactor.V
        reactor.jacobian(0.0, y0.copy(), zeros.copy(), 0.0)
        v_jacobian = reactor.compute_volume(y0)
        assert v_residual == v_jacobian
        # No reconstructed reverse rate for the irreversible attachment.
        assert float(reactor.kb[0]) == 0.0

    def test_12_analytic_electron_jacobian_column_matches_fd(self):
        """12. The analytic electron-state Jacobian column agrees with a
        finite-difference approximation of the same column."""
        reactor = self._initialized_reactor()
        ncs = reactor.num_core_species
        y0 = np.array(reactor.y0[:ncs], float)
        zeros = np.zeros(ncs, float)
        jac = reactor.jacobian(0.0, y0.copy(), zeros.copy(), 0.0)
        z = reactor.electron_index
        h = max(1.0e-7 * abs(y0[z]), 1.0e-9)
        yp = y0.copy(); yp[z] += h
        ym = y0.copy(); ym[z] -= h
        rp, _ = reactor.residual(0.0, yp, zeros.copy())
        rm, _ = reactor.residual(0.0, ym, zeros.copy())
        fd_col = (rp - rm) / (2.0 * h)
        assert np.abs(fd_col).max() > 0.0
        assert np.allclose(jac[:, z], fd_col, rtol=1e-4, atol=1e-8 * np.abs(fd_col).max())

    def test_13_running_pipeline_twice_does_not_accumulate_electrons(self):
        """13. Reconstructing the reactor does not accumulate additional
        electron participants: running the pipeline twice keeps exactly one
        incident electron on the packed reactant side, and the canonical is
        untouched. (See test_regression_reinit_same_reactor_does_not_accumulate
        for the same-object re-init variant.)"""
        for _ in range(2):
            reactor = self._initialized_reactor()
            assert self._electron_reactant_count(reactor, self.reaction) == 1
            # canonical never grows an explicit electron
            assert not any(spc.is_electron() for spc in self.reaction.reactants)
            assert self.reaction.electrons == -1

    def test_14_thermal_reactions_pass_through_unchanged(self):
        """14. Ordinary thermal reactions are untouched by the boundary: a
        reaction with no metadata electron is returned by identity, never
        wrapped, copied, or mutated."""
        thermal = Reaction(reactants=[self.o2], products=[self.o2_anion],
                           reversible=False, kinetics=_arrhenius())
        assert getattr(thermal, "electrons", 0) == 0
        reactor = self._reactor()
        resolved = reactor._resolve_electron_placements(
            [thermal], self._core_species(), [])
        assert resolved[0] is thermal  # identity: thermal path unchanged

    def test_15_reload_and_regenerate_uses_same_automatic_path(self):
        """15. Reloading the database and regenerating the reaction from scratch
        still exercises the same automatic production path -- a fresh canonical
        object reaches the reactor and is resolved there, not from a cache."""
        fresh_reaction = _load_family_and_generate_attachment()
        assert fresh_reaction is not self.reaction
        electron = _electron()
        core_species = [electron, fresh_reaction.reactants[0], fresh_reaction.products[0]]
        imf = {electron: Y_E0, fresh_reaction.reactants[0]: 1.0}
        reactor = PlasmaReactor(T_GAS, P0, imf, (T_E, "K"), n_sims=1, termination=[])
        reactor.initialize_model(core_species, [fresh_reaction], [], [])
        # The fresh canonical reaction bridges into reaction_index by identity
        # (not a cache, not a view) and was resolved anew AT the reactor: the
        # solver packed an incident electron on its reactant side, with kf > 0.
        assert fresh_reaction in reactor.reaction_index
        ie = reactor.electron_index
        j = reactor.reaction_index[fresh_reaction]
        assert any(reactor.reactant_indices[j, m] == ie
                   for m in range(reactor.reactant_indices.shape[1]))
        assert reactor.kf[0] > 0.0
        assert fresh_reaction.electrons == -1

    def test_16_generated_attachment_remains_irreversible(self):
        """16. The generated attachment reaction remains explicitly
        irreversible: the view is irreversible and the reactor builds no
        reverse rate for it."""
        reactor = self._initialized_reactor()
        packed = self._packed_reaction(reactor)
        assert packed.reversible is False
        assert float(reactor.kb[0]) == 0.0
        assert np.isinf(reactor.Keq[0])

    # --- Production-path assertion ------------------------------------------

    def test_production_path_invokes_resolver_not_the_test(self, monkeypatch):
        """The resolver is invoked by PRODUCTION code (initialize_model), not by
        this test. A spy on the module-level resolver records that the reactor
        called it exactly once, with the CANONICAL reaction (production made the
        view, the test did not)."""
        calls = []
        original = electron_placement.resolve_electron_placement

        def spy(reaction, species_list):
            calls.append(reaction)
            return original(reaction, species_list)

        monkeypatch.setattr(electron_placement, "resolve_electron_placement", spy)
        reactor = self._reactor()
        reactor.initialize_model(self._core_species(), [self.reaction], [], [])

        assert len(calls) == 1
        assert calls[0] is self.reaction  # production handed the canonical in
        assert reactor.kf[0] > 0.0

    # --- Negative controls --------------------------------------------------

    def test_negative_control_missing_electron_species(self):
        """Missing canonical electron species: the attachment reaction is in the
        model but no reactor electron species exists. The production path must
        FAIL with the resolver's named representation error, not silently omit
        the reaction from the reactor's reaction set."""
        core_species = [self.o2, self.o2_anion]  # no electron
        imf = {self.o2: 1.0}
        reactor = PlasmaReactor(T_GAS, P0, imf, (T_E, "K"), n_sims=1, termination=[])
        with pytest.raises(ElectronPlacementError, match="[Nn]o electron species"):
            reactor.initialize_model(core_species, [self.reaction], [], [])

    def test_negative_control_duplicate_explicit_representation(self):
        """Duplicate explicit representation: a reaction that already carries an
        explicit electron reactant AND a nonzero metadata electron count. The
        production path must fail by name (double representation), never
        double-count the electron."""
        duplicated = TemplateReaction(
            reactants=[self.o2, self.electron], products=[self.o2_anion],
            family=FAMILY, electrons=-1, reversible=False, is_forward=True,
            kinetics=_arrhenius())
        reactor = self._reactor()
        with pytest.raises(ElectronPlacementError, match="twice|double"):
            reactor.initialize_model(self._core_species(), [duplicated], [], [])

    def test_negative_control_unsupported_ionization_shaped(self):
        """Unsupported ionisation-shaped case: net electron PRODUCTION
        (electrons=+1) does not determine incident-electron order. The
        production path must HARD-FAIL through the resolver's family-declaration
        contract -- proving it does NOT take the shortcut
        `incident order = abs(Reaction.electrons)`, which would (wrongly)
        succeed here because abs(+1) == the declared count of 1."""
        ionization_shaped = TemplateReaction(
            reactants=[self.o2], products=[self.o2_anion],
            family=FAMILY, electrons=1, reversible=False, is_forward=True,
            kinetics=_arrhenius())
        reactor = self._reactor()
        with pytest.raises(ElectronPlacementError, match="ionization-shaped"):
            reactor.initialize_model(self._core_species(), [ionization_shaped], [], [])

    # --- Merge-gate regression coverage -------------------------------------

    def test_regression_plasma_kinetics_without_incident_electron_refused(self):
        """BLOCKER 2: a plasma-shaped (electron-temperature dependent) rate law
        with NO explicit incident electron and a zero net electron count skips
        the metadata gate and the resolver entirely, so the reactor itself must
        refuse it -- otherwise it would be evaluated at the wrong reaction order
        (missing its electron-density factor)."""
        plasma_no_electron = Reaction(
            reactants=[self.o2], products=[self.o2_anion], reversible=False,
            kinetics=TwoTemperaturePlasma(A=(1.0e-3, "m^3/(mol*s)"), n=0.5,
                                          Ea_g=(0.0, "J/mol"), Ea_e=(0.0, "J/mol")))
        assert getattr(plasma_no_electron, "electrons", 0) == 0  # skips the gate
        assert not any(spc.is_electron() for spc in plasma_no_electron.reactants)
        reactor = self._reactor()
        with pytest.raises(PlasmaStateError, match="plasma-shaped kinetics"):
            reactor.initialize_model(self._core_species(), [plasma_no_electron], [], [])

    def test_regression_reinit_same_reactor_does_not_accumulate(self):
        """BLOCKER 3: re-initializing the SAME reactor with the same canonical
        metadata reaction must not accumulate stale resolved-view keys in
        reaction_index. (Assertion 13 missed this by building a fresh reactor
        each loop.)"""
        reactor = self._initialized_reactor()
        assert len(reactor.reaction_index) == 1
        # Re-initialize the same reactor object with the same canonical reaction.
        reactor.initialize_model(self._core_species(), [self.reaction], [], [])
        assert len(reactor.reaction_index) == 1
        reactor.initialize_model(self._core_species(), [self.reaction], [], [])
        assert len(reactor.reaction_index) == 1

    def test_regression_reaction_index_maps_canonical_not_view(self):
        """MAJOR: the public reaction_index must bridge back to the canonical
        (model) reaction a caller passed in, not the reactor-internal resolved
        view. Covers edge promotion: a metadata reaction supplied as an edge
        reaction, then promoted to core, keeps model identity at the correct
        positional index without leaking a view."""
        # (a) core: canonical reaction is the key, no view leaks in.
        reactor = self._initialized_reactor()
        keys = list(reactor.reaction_index.keys())
        assert self.reaction in reactor.reaction_index  # identity, not a view
        assert reactor.reaction_index[self.reaction] == 0
        assert all(k is self.reaction for k in keys)
        assert all(k.electrons == -1 for k in keys)  # canonical, un-resolved

        # (b) edge promotion: same canonical reaction as an EDGE reaction beside
        # a thermal core reaction, then promoted to core. Identity is preserved
        # and the index tracks the position, with no accumulation.
        thermal = Reaction(reactants=[self.o2], products=[self.o2_anion],
                           reversible=False, kinetics=_arrhenius())
        reactor.initialize_model(self._core_species(), [thermal], [], [self.reaction])
        assert reactor.reaction_index[thermal] == 0
        assert reactor.reaction_index[self.reaction] == 1  # edge, after core
        assert len(reactor.reaction_index) == 2

        reactor.initialize_model(self._core_species(), [thermal, self.reaction], [], [])
        assert reactor.reaction_index[self.reaction] == 1  # promoted into core
        assert len(reactor.reaction_index) == 2  # no stale view left behind


# ---------------------------------------------------------------------------
# I-088: the reverse-generated cation recombination.
#
# The chemistry this project actually needs is generated in the family's
# REVERSE direction: RMG has only the neutral CH3Li in the model, matches it
# against Cation_R_Recombination's product template, and reconstructs the
# reactant side Li+ + CH3. `_create_reaction` stores that back in family-forward
# molecular orientation (reactants=[Li+, CH3], products=[CH3Li]) while flagging
# is_forward=False to record how it was FOUND (family.py:1753-1770, I-086). The
# family is reversible=True with an auto-derived reverse.
#
# So the reaction that must resolve carries is_forward=False AND reversible=True,
# and its electron is physically a REACTANT (electrons=-1). These tests assert
# the two quantities the ticket names -- the side the electron lands on and how
# many land there -- plus the structural proof that the side is right (the E
# pseudo-element and the net charge both balance across the view), and the proof
# that placement did not launder the reverse-rate problem: the view is still
# reversible, and the reactor still refuses it for that reason, by name.
# ---------------------------------------------------------------------------

CATION_FAMILY = "Cation_R_Recombination"
NEUTRAL_CH3LI = "C[Li]"


def _generate_reverse_cation_recombination():
    """Load the real ``Cation_R_Recombination`` family and reach the lithium
    cation the way the model builder does: feed the NEUTRAL product ``CH3Li``
    and let the family generate in its reverse direction. Returns the canonical
    :class:`TemplateReaction` with kinetics filled in from the family's own rate
    rules. No library lookup, no hand-built reaction."""
    families_path = os.path.join(settings["database.directory"], "kinetics", "families")
    database = KineticsDatabase()
    database.load_recommended_families(os.path.join(families_path, "recommended.py"))
    database.load_families(families_path, families=[CATION_FAMILY])

    species = Species().from_smiles(NEUTRAL_CH3LI)
    species.generate_resonance_structures()
    reactions = database.generate_reactions_from_families(
        [species], only_families=[CATION_FAMILY])

    # Li+ + CH3 <=> CH3Li: two reactants, one of them the lithium cation, and a
    # single neutral product. There must be exactly one.
    matches = [rxn for rxn in reactions
               if len(rxn.reactants) == 2 and len(rxn.products) == 1
               and any(_net_charge(spc) == 1 for spc in rxn.reactants)]
    assert len(matches) == 1, (
        "expected exactly one reverse Li+ + CH3 <=> CH3Li reaction, got "
        "{0}: {1}".format(len(matches), [str(r) for r in reactions]))
    reaction = matches[0]

    # Kinetics from the family's own rate rules, as the model builder does. A
    # fresh family is loaded because add_rules_from_training and
    # fill_rules_by_averaging_up mutate the family's rule table.
    kinetics_database = KineticsDatabase()
    kinetics_database.load_families(families_path, families=[CATION_FAMILY])
    family = kinetics_database.families[CATION_FAMILY]
    family.add_rules_from_training(thermo_database=None)
    family.fill_rules_by_averaging_up()
    template = family.retrieve_template(reaction.template)
    reaction.kinetics = family.get_kinetics_for_template(
        template, degeneracy=reaction.degeneracy)[0]
    return reaction


@pytest.mark.database
class TestCationRecombinationPlacement:
    """I-088: placement resolves for the reverse-generated, reversible cation
    recombination, and puts the electron on the reactant side exactly once."""

    @classmethod
    def setup_class(cls):
        cls.reaction = _generate_reverse_cation_recombination()
        cls.electron = _electron()
        cls.species_list = ([cls.electron] + list(cls.reaction.reactants)
                            + list(cls.reaction.products))

    def _resolve(self):
        return resolve_electron_placement(self.reaction, self.species_list)

    # -- preconditions: this really is the shape the three refusals rejected ---

    def test_00_reaction_is_reverse_generated_and_reversible(self):
        """The reaction under test carries ALL THREE properties the resolver
        used to refuse: an undeclared family, is_forward=False, reversible=True.
        If any of these drifts, the rest of this class stops testing I-088."""
        rxn = self.reaction
        assert rxn.family == CATION_FAMILY
        assert rxn.is_forward is False
        assert rxn.reversible is True
        assert rxn.electrons == -1
        # Family-forward orientation: the CATION is on the reactant side, so the
        # electron must be a reactant to make up the missing negative charge.
        assert sorted(_net_charge(spc) for spc in rxn.reactants) == [0, 1]
        assert [_net_charge(spc) for spc in rxn.products] == [0]

    # -- the two quantities the ticket names ----------------------------------

    def test_01_electron_placed_on_reactant_side_exactly_once(self):
        """THE assertion: the electron lands on the REACTANT side, and exactly
        one of it does. A product-side placement or a count of 0 or 2 fails."""
        view = self._resolve()
        assert sum(1 for spc in view.reactants if spc.is_electron()) == 1
        assert sum(1 for spc in view.products if spc.is_electron()) == 0
        assert any(spc is self.electron for spc in view.reactants)
        assert view.electrons == 0
        assert len(view.reactants) == len(self.reaction.reactants) + 1
        assert len(view.products) == len(self.reaction.products)
        assert get_molecularity(view) == 3

    def test_02_view_balances_in_charge_and_E_pseudo_element(self):
        """Structural proof that the side is right, independent of any metadata:
        with the electron on the reactant side both the net electrical charge and
        the ``E`` pseudo-element balance across the view. Placing it on the
        product side would leave charge -1 vs 0 and E -1 vs +1."""
        view = self._resolve()
        assert (sum(_net_charge(spc) for spc in view.reactants)
                == sum(_net_charge(spc) for spc in view.products) == 0)
        # The Chemkin writer's own guard, applied to the view's two sides. It
        # raises MechanismWriterError on an E imbalance; silence is the assertion.
        check_electron_balance(view, view.reactants, view.products, str(view))

    def test_03_canonical_reaction_is_not_mutated(self):
        """The database representation survives resolution untouched."""
        reactants_before = list(self.reaction.reactants)
        products_before = list(self.reaction.products)
        self._resolve()
        assert self.reaction.electrons == -1
        assert self.reaction.reactants == reactants_before
        assert self.reaction.products == products_before
        assert not any(spc.is_electron() for spc in self.reaction.reactants)
        assert not any(spc.is_electron() for spc in self.reaction.products)

    # -- placement must not launder the reverse-rate problem -------------------

    def test_04_view_stays_reversible(self):
        """Placement converts a REPRESENTATION, it does not decide a rate
        policy. The view carries reversible=True through unchanged, so the
        reactor's reverse-rate guard still sees what it needs to see."""
        assert self._resolve().reversible is True

    def test_05_reactor_still_refuses_the_reversible_view_by_name(self):
        """The protection the resolver's own ``reversible`` refusal provided is
        NOT lost: PlasmaReactor._validate_reactions refuses a reversible
        electron-containing reaction by name, because kr = kf/Keq(Tgas) would
        price the electron's thermochemistry at the gas temperature. This is the
        tripwire for I-088's removal of that check from the resolver -- if the
        reactor ever stops refusing, this test goes red."""
        view = self._resolve()
        core_species = [self.electron] + list(view.reactants[:-1]) + list(view.products)
        imf = {self.electron: Y_E0, view.reactants[0]: 1.0}
        reactor = PlasmaReactor(T_GAS, P0, imf, (T_E, "K"), n_sims=1, termination=[])
        with pytest.raises(NonEquilibriumReverseRateError,
                           match="electron-containing reaction"):
            reactor.initialize_model(core_species, [view], [], [])
