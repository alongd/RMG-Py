#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu),
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module contains methods for generation of unidirectional transitions of molecules.
These transitions are unsymmetrical in the sense that input molecule A may return molecule B, but not vice versa.
For example, if [O.][:S.]=O is formed during model generation, the methods here will identify it as O=S=O, but O=S=O
will not return the open shell [O.][:S.]=O structure.

Currently supported transitions:
    - ``generate_birad_multiple_bond_unidirectional_transitions``: two adjacent radicals shift with a multiple bond (unidirectional)
    - ``generate_OS_unidirectional_transitions``: transitions of excited O=O (or its sulfur-containing isoelectronic structures) into the ground [O][O] state (unidirectional)
"""

import cython

from .graph import Vertex
from .molecule import Atom, Bond, Molecule
import rmgpy.molecule.pathfinder as pathfinder
from rmgpy.exceptions import AtomTypeError
import rmgpy.molecule.filtration as filtration


def generate_birad_multiple_bond_unidirectional_transitions(mol):
    """
    Generate all structures formed by relaxation of an excited two adjacent radicals state to a
    multiple bond state. This is a unidirectional transition (i.e., if the [O.][S.]=O structure is present in the system
    it should be marked as `.reactive=False` and the reactive structure O=S=O is added to the molecule list. However, if
    O=S=O is present is the system, this doesn't necessarily mean that the excited [O.][S.]=O structure is generated.
    This function should NOT convert ground state [O.][O.] into an excited O=O state (or its isoelectronic structures
    [S.][S.] and [S.][O.]).
    This should target the following transitions:
    - [:N.]=[:N.] => N#N
    - [O.][S.](=O)=O => O=S(=O)=O
    - [O.][:S.]=O => O=S=O
    - C=[S+][C.]([O.])[S-] => C=[S+]C(=O)[S-]
    - [CH2.][CH2.] => C=C
    (where ':' denotes a lone pair, '.' denotes a radical, '-' not in [] denotes a single bond, '-'/'+' denote charge)
    * This is treated as a unidirectional transition and not a family since: (a) otherwise the excited species and all
    of its resonance structures will still be reactive and contribute to the combinatorial complexity of the system;
    (b) each model that involves [O][O] will have (at least in the edge) the excited O=O state due to the reverse
    reaction in the family.
    """
    cython.declare(structures=list, paths=list, index=cython.int, structure=Molecule)
    cython.declare(atom=Atom, atom1=Atom, atom2=Atom, bond12=Bond)
    cython.declare(v1=Vertex, v2=Vertex)

    # from rmgpy.rmg.input import get_species_constaints

    structures = []
    # if mol.isRadical() and pathfinder.is_OS(mol) == 0 and not get_species_constaints('allowAdjacentRadicals'):
    if mol.isRadical() and filtration.is_OS(mol) == 0:
        # This is neither O2, S2, SO, and the user did not specify "allowAdjacentRadicals=True" in the input file
        for atom in mol.vertices:
            paths = pathfinder.find_birad_multiple_bond_delocalization_paths(atom)
            for atom1, atom2, bond12 in paths:
                atom1.decrementRadical()
                atom2.decrementRadical()
                bond12.incrementOrder()
                atom1.updateCharge()
                atom2.updateCharge()
                # Make a copy of structure
                structure = mol.copy(deep=True)
                # Also copy the connectivity values, since they are the same for all structures
                for index in xrange(len(mol.vertices)):
                    v1 = mol.vertices[index]
                    v2 = structure.vertices[index]
                    v2.connectivity1 = v1.connectivity1
                    v2.connectivity2 = v1.connectivity2
                    v2.connectivity3 = v1.connectivity3
                    v2.sortingLabel = v1.sortingLabel
                # Restore current structure
                atom1.incrementRadical()
                atom2.incrementRadical()
                bond12.decrementOrder()
                atom1.updateCharge()
                atom2.updateCharge()
                try:
                    # update both atomTypes and multiplicity
                    structure.update(log_species_while_updating_atom_types=False)
                except AtomTypeError:
                    pass  # Don't append structure if it creates an undefined atomType
                else:
                    structures.append(structure)
    return structures


def generate_OS_unidirectional_transitions(mol):
    """
    This is function deals with the opposite case for generate_birad_multiple_bond_unidirectional_transitions(),
    addressing three specific isoelectronic structures: O=O, S=S, and S=O. Here the transition is from a double bond
    structure into a single bond structure with two adjacent radicals. This is also a unidirectional transition (i.e.,
    if the O=O structure is present in the system it should be marked as `.reactive=False` and the reactive structure
    [O.][O.] is added to the molecule list.
    This function assumes that mol consists of two atoms only, which should be checked by the calling function.
    This transition function is disabled if the user chooses to allow singlet oxygen in the input file (in which case,
    the user is responsible to rates for the transition in the form of a family or a library).
    """
    cython.declare(index=cython.int, structure=Molecule)
    cython.declare(atom=Atom, atom1=Atom, atom2=Atom, bond12=Bond)
    cython.declare(v1=Vertex, v2=Vertex)

    # from rmgpy.rmg.input import get_species_constaints

    # if pathfinder.is_OS(mol) == 2 and not get_species_constaints('allowSingletO2'):  # check if mol is the excited state, e.g., O=O
    if filtration.is_OS(mol) == 2:  # check if mol is the excited state, e.g., O=O
        # This is either O=O, S=S, or S=O, and the user did not specify "allowSingletO2=True" in the input file
        mol.vertices[0].incrementRadical()
        mol.vertices[1].incrementRadical()
        mol.vertices[0].bonds[mol.vertices[1]].decrementOrder()
        # Make a copy of structure
        structure = mol.copy(deep=True)
        # Also copy the connectivity values, since they are the same for all structures
        for index in xrange(len(mol.vertices)):
            v1 = mol.vertices[index]
            v2 = structure.vertices[index]
            v2.connectivity1 = v1.connectivity1
            v2.connectivity2 = v1.connectivity2
            v2.connectivity3 = v1.connectivity3
            v2.sortingLabel = v1.sortingLabel
        # Restore current structure
        mol.vertices[0].decrementRadical()
        mol.vertices[1].decrementRadical()
        mol.vertices[0].bonds[mol.vertices[1]].incrementOrder()
        try:
            structure.update(log_species_while_updating_atom_types=False)  # update both atomTypes and multiplicity
        except AtomTypeError:
            pass  # Don't append structure if it creates an undefined atomType
        else:
            return [structure]
    return []
