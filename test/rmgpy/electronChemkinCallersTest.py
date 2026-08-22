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
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL     #
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER  #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
Tests for the caller sites where ``Reaction.to_chemkin`` was invoked without a usable species
list and would crash on a charged reaction once the missing/wrong-list holes were closed.

Two of the six sites have GENUINE regression locks that fail when their production line is
reverted (verified by actual revert, RED evidence in the review record):

- rmgpy/solver/base.pyx:1280 (sensitivity-CSV header) ->
  test/rmgpy/solver/simpleTest.py::SimpleReactorTest::test_sensitivity_header_places_electron_for_charged_reaction
- rmgpy/rmg/output.py (edge HTML render) ->
  test/rmgpy/rmg/outputTest.py::TestOutput::test_save_output_html_edge_places_electron_for_charged_reaction

The three sites below cannot be driven end-to-end in this test environment without
disproportionate or unavailable scaffolding; each test says so and binds to the real production
object/attribute where that object is constructible. What each proves, and why it cannot
fail-on-revert, is stated in its docstring -- a reasoned exception, not silence.
"""

from rmgpy.chemkin import get_species_identifier
from rmgpy.exceptions import MechanismWriterError
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species

import pytest


class TestElectronSurvivesChemkinCallerSites:
    """Charged-reaction coverage for the three caller sites whose enclosing method is not
    drivable here (see module docstring for the two that are)."""

    def setup_class(self):
        self.neutral = Species(label="OH", molecule=[Molecule().from_smiles("[OH]")])
        self.anion = Species(label="OHm", molecule=[Molecule().from_smiles("[OH-]")])
        self.electron = Species(label="e", molecule=[Molecule().from_smiles("e")])
        self.electron_id = get_species_identifier(self.electron)
        # Attachment OH + e- -> OH-: electrons=-1 folds the electron onto the reactant side and
        # the anion carries the charge, so it balances in the E pseudo-element.
        self.charged = Reaction(reactants=[self.neutral], products=[self.anion], electrons=-1,
                                kinetics=Arrhenius(A=(1e13, "cm^3/(mol*s)"), n=0, Ea=(10, "kJ/mol")))
        self.species_list = [self.neutral, self.anion, self.electron]

    def test_local_uncertainty_kinetics_intermediate(self):
        """rmgpy/tools/uncertainty.py:1532 -- Uncertainty.local_analysis_intermediate builds
        self.all_kinetics_intermediates from
        self.reaction_list[i].to_chemkin(species_list=self.species_list, kinetics=False).

        Driving local_analysis_intermediate needs full local-sensitivity results (a reaction
        system, DASPK sensitivity output, MUQ-free but still disproportionate scaffolding), so
        this binds to the REAL production attributes instead: it constructs the actual Uncertainty
        object and reads self.species_list / self.reaction_list -- the exact attributes the fixed
        line reads, paired the same way. It does not call the method, so it does not fail if that
        line is reverted; what it proves is that the paired lists the object holds serialize a
        charged reaction with its electron rather than dropping it."""
        from rmgpy.tools.uncertainty import Uncertainty

        u = Uncertainty(species_list=self.species_list, reaction_list=[self.charged])
        label = 'k' + str(u.reaction_list[0].index) + ': ' + u.reaction_list[0].to_chemkin(
            species_list=u.species_list, kinetics=False)
        assert self.electron_id in label, label

    def test_collision_rate_violator_report(self):
        """rmgpy/rmg/main.py:1676 -- RMG.check_model's collision-rate-violators report, now
        violator[0].to_chemkin(self.reaction_model.core.species) unconditionally (kinetics=True).

        Reasoned exception: driving check_model needs a fully populated RMG object AND a charged
        *bimolecular* reaction whose forward rate exceeds the collision limit so it is collected
        as a violator (a unimolecular attachment is never collision-limited); RMG itself imports
        h5py, which is binary-incompatible with the pinned numpy in this env. That is
        disproportionate scaffolding for one report line. This proves the fixed call form
        -- to_chemkin(core.species) with kinetics=True -- places the electron and survives
        write_kinetics_entry, and that the previous no-list form raises."""
        entry = self.charged.to_chemkin(self.species_list)
        assert self.electron_id in entry, entry
        with pytest.raises(MechanismWriterError, match="electron"):
            self.charged.to_chemkin()

    def test_global_uncertainty_reaction_description(self):
        """rmgpy/tools/globaluncertainty.py:471 -- ReactorPCEFactory.analyze_results builds a
        rate-sensitivity description from
        cantera.reaction_list[d].to_chemkin(species_list=cantera.species_list, kinetics=False).

        Reasoned exception: globaluncertainty.py imports muq at module top level, and its Cantera
        object imports the cantera package; neither muq nor cantera is installed in this env, so
        neither the module nor its real object can be imported at all -- the enclosing method
        cannot be reached or its attribute read here. This proves the with-list call form the
        fixed line uses places the electron; the review independently confirmed at source that
        cantera.species_list is the RMG Species list paired with cantera.reaction_list."""
        description = 'dln[{0}]/dln[{1}]'.format(
            'OH', self.charged.to_chemkin(species_list=self.species_list, kinetics=False))
        assert self.electron_id in description, description
        with pytest.raises(MechanismWriterError, match="electron"):
            self.charged.to_chemkin(kinetics=False)
