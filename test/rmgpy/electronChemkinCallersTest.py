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
Regression tests for the four caller sites where ``Reaction.to_chemkin`` was
invoked with no ``species_list`` and only worked while ``write_reaction_string``
tolerated a missing list. Once that tolerance became a hard
``MechanismWriterError`` for charged reactions, each site had to pass the
species list already in its scope. These tests pin, per site, that the list the
fix passes is the load-bearing thing: the fixed call form places the metadata
electron, and the previous no-list form raises.

Each enclosing method (ReactionSystem.simulate, RMG.check_model,
ReactorPCEFactory.analyze_results, Uncertainty.local_analysis_intermediate) needs
heavy scaffolding (DASPK integration, a full RMG object, MUQ), so these exercise
the exact ``to_chemkin`` call expression each site now uses rather than driving
the method end-to-end -- the construction-level test the review explicitly
sanctioned for the sensitivity-header path, applied to all four.
"""

import pytest

from rmgpy.chemkin import get_species_identifier
from rmgpy.exceptions import MechanismWriterError
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species


class TestElectronSurvivesChemkinCallerSites:
    """One charged-reaction test per caller site fixed in review round 14."""

    def setup_class(self):
        self.neutral = Species(label="OH", molecule=[Molecule().from_smiles("[OH]")])
        self.anion = Species(label="OHm", molecule=[Molecule().from_smiles("[OH-]")])
        self.electron = Species(label="e", molecule=[Molecule().from_smiles("e")])
        self.electron_id = get_species_identifier(self.electron)
        # Attachment OH + e- -> OH-: electrons=-1 folds the electron onto the reactant
        # side, keeping the written reaction bimolecular (cm^3/(mol*s) A-factor).
        self.charged = Reaction(reactants=[self.neutral], products=[self.anion], electrons=-1,
                                kinetics=Arrhenius(A=(1e13, "cm^3/(mol*s)"), n=0, Ea=(10, "kJ/mol")))
        self.species_list = [self.neutral, self.anion, self.electron]

    def test_base_solver_sensitivity_header(self):
        """rmgpy/solver/base.pyx:1280 -- sensitivity-CSV header in ReactionSystem.simulate.
        The header format is 'dln[..]/dln[k{j}]: {rxn}' with the reaction rendered via
        core_reactions[j].to_chemkin(species_list=core_species, kinetics=False)."""
        header = 'dln[{0}]/dln[k{1}]: {2}'.format(
            'OH', 1, self.charged.to_chemkin(species_list=self.species_list, kinetics=False))
        assert self.electron_id in header, header
        with pytest.raises(MechanismWriterError, match="electron"):
            self.charged.to_chemkin(kinetics=False)

    def test_collision_rate_violator_report(self):
        """rmgpy/rmg/main.py:1676 -- collision-rate-violators report in RMG.check_model,
        now violator[0].to_chemkin(self.reaction_model.core.species) unconditionally
        (kinetics=True, so the electron must survive write_kinetics_entry too)."""
        entry = self.charged.to_chemkin(self.species_list)
        assert self.electron_id in entry, entry
        with pytest.raises(MechanismWriterError, match="electron"):
            self.charged.to_chemkin()

    def test_global_uncertainty_reaction_description(self):
        """rmgpy/tools/globaluncertainty.py:471 -- rate-sensitivity description in
        ReactorPCEFactory.analyze_results, now
        cantera.reaction_list[d].to_chemkin(species_list=cantera.species_list, kinetics=False)."""
        description = 'dln[{0}]/dln[{1}]'.format(
            'OH', self.charged.to_chemkin(species_list=self.species_list, kinetics=False))
        assert self.electron_id in description, description
        with pytest.raises(MechanismWriterError, match="electron"):
            self.charged.to_chemkin(kinetics=False)

    def test_local_uncertainty_kinetics_intermediate(self):
        """rmgpy/tools/uncertainty.py:1532 -- kinetics intermediate label in
        Uncertainty.local_analysis_intermediate, now
        reaction_list[i].to_chemkin(species_list=self.species_list, kinetics=False)."""
        label = 'k' + str(self.charged.index) + ': ' + self.charged.to_chemkin(
            species_list=self.species_list, kinetics=False)
        assert self.electron_id in label, label
        with pytest.raises(MechanismWriterError, match="electron"):
            self.charged.to_chemkin(kinetics=False)
