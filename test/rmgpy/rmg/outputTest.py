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

import os
import shutil


from rmgpy.rmg.model import CoreEdgeReactionModel, ReactionModel
from rmgpy.rmg.output import save_output_html
from rmgpy.chemkin import load_chemkin_file


class TestOutput:
    def test_save_output_html(self):
        """
        This example is to test if an HTML file can be generated
        for the provided chemkin model.
        """
        folder = os.path.join(os.path.dirname(__file__), "..", "test_data", "saveOutputHTML")

        chemkin_path = os.path.join(folder, "eg6", "chem_annotated.inp")
        dictionary_path = os.path.join(folder, "eg6", "species_dictionary.txt")

        # load_chemkin_file
        species, reactions = load_chemkin_file(chemkin_path, dictionary_path)

        # convert it into a reaction model:
        core = ReactionModel(species, reactions)
        cerm = CoreEdgeReactionModel(core)

        out = os.path.join(folder, "output.html")
        save_output_html(out, cerm)

        assert os.path.isfile(out)
        os.remove(out)
        shutil.rmtree(os.path.join(folder, "species"))

    def test_save_output_html_edge_places_electron_for_charged_reaction(self):
        """
        Regression lock for the edge render in save_output_html (rmgpy/rmg/output.py).

        The edge page serializes reactions against a list that must include the electron, which is
        core reactor state absent from edge.species. Before the fix the template rendered
        rxn.to_chemkin against the displayed edge `species` list (edge.species + output), so a
        charged EDGE reaction crashed in expand_electrons ('no electron species'). The fix keeps
        that displayed list untouched and serializes against a complete `serialization_species`
        (core + edge + output). Reverting that single line makes this test fail.
        """
        from copy import deepcopy
        from rmgpy.species import Species
        from rmgpy.molecule import Molecule
        from rmgpy.data.kinetics.library import LibraryReaction
        from rmgpy.kinetics import Arrhenius
        from rmgpy.chemkin import get_species_identifier

        folder = os.path.join(os.path.dirname(__file__), "..", "test_data", "saveOutputHTML")
        chemkin_path = os.path.join(folder, "eg6", "chem_annotated.inp")
        dictionary_path = os.path.join(folder, "eg6", "species_dictionary.txt")
        species, reactions = load_chemkin_file(chemkin_path, dictionary_path)
        cerm = CoreEdgeReactionModel(ReactionModel(species, reactions))

        # A genuine attachment OH + e- -> OH-, which balances in the E pseudo-element (the anion
        # carries the charge). Borrow real thermo from a loaded species so the template's
        # enthalpy/entropy rendering works; the electron lives in the core, never in edge.species.
        neutral = Species(label="OHrad", molecule=[Molecule().from_smiles("[OH]")])
        anion = Species(label="OHminus", molecule=[Molecule().from_smiles("[OH-]")])
        neutral.thermo = deepcopy(species[0].thermo)
        anion.thermo = deepcopy(species[0].thermo)
        electron = Species(label="eNEG", molecule=[Molecule().from_adjacency_list("1 e u1 p0 c-1")])
        cerm.core.species.append(electron)
        cerm.edge.species.extend([neutral, anion])

        # electrons=-1 folds the electron onto the reactant side, so the written reaction is
        # bimolecular (cm^3/(mol*s) A-factor).
        charged_edge_rxn = LibraryReaction(
            index=9001, library="electron_test", reversible=False,
            reactants=[neutral], products=[anion], electrons=-1,
            kinetics=Arrhenius(A=(1e13, "cm^3/(mol*s)"), n=0, Ea=(0, "kJ/mol")))
        cerm.edge.reactions.append(charged_edge_rxn)

        out = os.path.join(folder, "output_edge.html")
        try:
            save_output_html(out, cerm, part_core_edge="edge")
            assert os.path.isfile(out)
            with open(out) as f:
                html = f.read()
            assert get_species_identifier(electron) in html
        finally:
            if os.path.isfile(out):
                os.remove(out)
            sp_dir = os.path.join(folder, "species")
            if os.path.isdir(sp_dir):
                shutil.rmtree(sp_dir)
