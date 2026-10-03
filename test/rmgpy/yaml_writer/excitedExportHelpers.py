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

"""Shared export fixtures and flux witness, loaded by exact file path."""

import numpy as np

from rmgpy.data.base import Entry
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species
from rmgpy.statmech import Conformer
from rmgpy.thermo import NASA, NASAPolynomial


def nitrogen(label='N2', level=-1):
    mol = Molecule(smiles='N#N')
    mol.vibrational_level = level
    return Species(label=label, molecule=[mol], conformer=Conformer(), thermo=NASA(polynomials=[
        NASAPolynomial(coeffs=[3.5, 0, 0, 0, 0, 0, 2], Tmin=(200, 'K'), Tmax=(1000, 'K')),
        NASAPolynomial(coeffs=[3.5, 0, 0, 0, 0, 0, 2], Tmin=(1000, 'K'), Tmax=(3000, 'K')),
    ], Tmin=(200, 'K'), Tmax=(3000, 'K')))


def database(cls, level, rate):
    obj = cls(name='probe')
    obj.entries[1] = Entry(index=1, label='A <=> B',
        item=Reaction(reactants=[nitrogen('A')], products=[nitrogen('B', level)]),
        data=Arrhenius(A=(rate, 's^-1'), Ea=(0, 'J/mol')))
    return obj


class FluxDiagramWitness:
    def test_review_7_flux_nodes_keep_identity(self, tmp_path, monkeypatch):
        import pydot
        from rmgpy.tools import fluxdiagram
        ground, excited = nitrogen(), nitrogen(level=1)
        rxn = Reaction(reactants=[ground], products=[excited], pairs=[(ground, excited)])
        model = CoreEdgeReactionModel()
        model.core.species = [ground, excited]; model.core.reactions = [rxn]
        monkeypatch.setattr(fluxdiagram, 'video_fps', 1)
        monkeypatch.setattr(fluxdiagram, 'initial_padding', 1)
        monkeypatch.setattr(fluxdiagram, 'final_padding', 1)
        fluxdiagram.generate_flux_diagram(model, np.array([1.]), np.array([[1., 2.]]), np.array([[1.]]), str(tmp_path))
        graph = pydot.graph_from_dot_file(str(tmp_path / 'flux_diagram_0001.dot'))[0]
        nodes = {n.get_name().strip('"') for n in graph.get_nodes()} - {'node', 'graph', 'edge', '\\n'}
        assert nodes == {'N2', 'N2(v1)'}
        assert {(e.get_source().strip('"'), e.get_destination().strip('"')) for e in graph.get_edges()} == {('N2', 'N2(v1)')}
