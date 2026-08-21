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

import pytest

from rmgpy import settings
from rmgpy.data.kinetics.common import get_molecularity
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.electron_placement import FAMILY_ELECTRON_PLACEMENT, resolve_electron_placement
from rmgpy.exceptions import ElectronPlacementError, PlasmaStateError
from rmgpy.kinetics import Arrhenius
from rmgpy.reaction import Reaction
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.species import Species

FAMILY = "Plasma_Electron_Attachment"
# Any of the other electron-bearing families that also declare electrons = -1;
# proves the resolver did not silently generalize beyond attachment.
OTHER_ELECTRON_FAMILY = "Cation_R_Recombination"

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

    def test_reverse_direction_fails_by_name(self):
        """The placement view is directional; a reverse-generated reaction is
        refused rather than given a forward-side electron."""
        reaction = _attachment(is_forward=False)
        with pytest.raises(ElectronPlacementError, match="reverse direction"):
            resolve_electron_placement(reaction, [_electron()])

    def test_reversible_fails_by_name(self):
        """The attachment declaration covers the irreversible forward form
        only; a reversible reaction's electron side is direction-ambiguous."""
        reaction = _attachment(reversible=True)
        with pytest.raises(ElectronPlacementError, match="reversible"):
            resolve_electron_placement(reaction, [_electron()])

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

    def test_declaration_registry_is_attachment_only(self):
        """This increment declares exactly one family; a broader registry
        would mean attachment semantics silently generalized."""
        assert FAMILY_ELECTRON_PLACEMENT == {FAMILY: ("reactants", 1)}


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

    def test_reactor_rejects_canonical_form(self):
        """The reactor refuses the metadata-only representation by name: the
        exact class of net-derived form the resolver exists to replace."""
        reactor = self._reactor()
        with pytest.raises(PlasmaStateError, match="metadata-only electron count"):
            reactor.initialize_model(self._core_species(), [self.reaction], [], [])
        # The refusal must not have mutated the canonical reaction either.
        assert self.reaction.electrons == -1

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
