#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),
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
This module contains methods for filtering a list of Molecules representing a single species,
keeping only the representative structures. Relevant for filtration of resonance structures.
"""

import logging
from .molecule import Molecule


def filter_structures(mol_list, mark_unreactive=True):
    """
    We often get too many resonance structures from the combination of all rules, particularly for species containing
    lone pairs. This method filters them out by minimizing the number of C/N/O/S atoms without a full octet.
    """

    # Get an octet deviation list
    octet_deviation_list = get_octet_deviation_list(mol_list)

    # Filter mol_list using the octet rule and the respective octet deviation list
    filtered_list, charge_span_list = octet_filtration(mol_list, octet_deviation_list)

    # Filter by charge span
    filtered_list = charge_span_filtration(filtered_list, charge_span_list)

    if mark_unreactive:
        # Mark selected unreactive structures if OS and/or adjacent birad unidirectional transitions were used
        mark_unreactive_structures(filtered_list, mol_list)

    # Check that there's at least one reactive structure in the list
    check_reactive(filtered_list)

    return filtered_list


def get_octet_deviation_list(mol_list):
    """
    Returns the a list of octet deviations for a respective list of :class:Molecule objects
    """

    octet_deviation_list = []
    for mol in mol_list:
        octet_deviation_list.append(get_octet_deviation(mol))

    return octet_deviation_list


def get_octet_deviation(mol):
    """
    Returns the octet deviation for a :class:Molecule object
    """

    assert isinstance(mol, Molecule), "Octet deviation could only be determined for Molecule objects."

    octet_deviation = 0  # This is the overall "score" for the molecule, summed across all C/N/O/S atoms
    for atom in mol.vertices:
        val_electrons = 2 * (int(atom.getBondOrdersForAtom()) + atom.lonePairs) + atom.radicalElectrons
        if atom.isCarbon():
            octet_deviation += abs(8 - val_electrons)  # expecting C to be near octet
        elif atom.isNitrogen():
            if atom.lonePairs:
                octet_deviation += abs(8 - val_electrons)  # expecting N p1/2/3 to be near octet
            else:
                octet_deviation += min(abs(10 - val_electrons), abs(8 - val_electrons))  # N p0 could be near octet or
                # dectet. N p0 could be closer to an octet rather than a dectet such as in O=[N+][O-]
            if val_electrons > 8:
                octet_deviation += 1  # penalty for N with valance greater than 8 (as in O=[N.]=O,
                # [NH2.]=[:NH.], N#[N.]O, CCN=N#N)
        elif atom.isOxygen():
            octet_deviation += abs(8 - val_electrons)  # expecting O to be near octet
            if atom.atomType.label in ['O4sc', 'O4dc', 'O4tc']:
                octet_deviation += 1  # penalty for O p1 c+1
                # as in [N-2][N+]#[O+], [O-]S#[O+], OS(S)([O-])#[O+], [OH+]=S(O)(=O)[O-], [OH.+][S-]=O.
                # [C-]#[O+] and [O-][O+]=O which are correct structures also get penalized here, but that's OK
                # since they are still eventually selected as representative structures according to the rules here.
        elif atom.isSulfur():
            if atom.lonePairs == 0:
                octet_deviation += abs(12 - val_electrons)  # duodectet on S p0, eg O=S(=O)(O)O val 12, O[S](=O)=O val 11
            elif atom.lonePairs == 1:
                octet_deviation += min(abs(8 - val_electrons), abs(10 - val_electrons))  # octet/dectet on S p1,
                # eg [O-][S+]=O val 8, O[S]=O val 9, OS([O])=O val 10
            elif atom.lonePairs == 2:
                octet_deviation += min(abs(8 - val_electrons), abs(10 - val_electrons))  # octet/dectet on S p2,
                # eg [S][S] val 7, OS[O] val 8, [NH+]#[N+][S-][O-] val 9, O[S-](O)[N+]#N val 10
            elif atom.lonePairs == 3:
                octet_deviation += abs(8 - val_electrons)  # octet on S p3, eg [S-][O+]=O
            for atom2, bond in atom.bonds.iteritems():
                if atom2.isSulfur() and bond.isTriple():
                    octet_deviation += 0.5  # penalty for S#S substructures. Often times sulfur can have a triple
                    # bond to another sulfur in a structure that obeys the octet rule, but probably shouldn't be a
                    # correct resonance structure. This adds to the combinatorial effect of resonance structures
                    # when generating reactions, yet probably isn't too important for reactivity. The penalty value
                    # is 0.5 since S#S substructures are captured twice (once for each S atom).
                    # Examples: CS(=O)SC <=> CS(=O)#SC;
                    # [O.]OSS[O.] <=> [O.]OS#S[O.] <=> [O.]OS#[S.]=O; N#[N+]SS[O-] <=> N#[N+]C#S[O-]
        if is_OS(mol) and atom.radicalElectrons >= 2:
            octet_deviation += atom.radicalElectrons + 1  # This helps to distinguish between a birad site and two
            # adjacent radicals on S2 or SO structures, e.g., [::O.][::S.] vs. [::O]=[:S..]

    return octet_deviation


def octet_filtration(mol_list, octet_deviation_list):
    """
    Returns the a filtered list based on the octet_deviation_list. Also computes and returns a charge_span_list.
    Filtering using the octet deviation criterion rules out most unrepresentative structures. However, since some
    charge-strained species are still kept (e.g., [NH]N=S=O <-> [NH+]#[N+][S-][O-]), we also generate during the same
    loop a charge_span_list to keep track of the charge spans. This is used for further filtering.
    """

    filtered_list = []
    charge_span_list = []
    for index in xrange(len(mol_list)):
        if octet_deviation_list[index] == min(octet_deviation_list):
            filtered_list.append(mol_list[index])
            charge_span_list.append(mol_list[index].getChargeSpan())

    return filtered_list, charge_span_list


def get_charge_span_list(mol_list):
    """
    Returns the a list of charge spans for a respective list of :class:Molecule objects
    This is also calculated in the octet_filtration() method along with the octet filtration process
    """

    charge_span_list = []
    for mol in mol_list:
        charge_span_list.append(mol.getChargeSpan())

    return charge_span_list


def charge_span_filtration(filtered_list, charge_span_list):
    """
    Returns the a new filtered_list, filtered based on charge_span_list
    If the species is a radical we first check whether keeping an extra charge span separation might be important for
    reactivity by relocating the radical site. If so, wee keep these structures.
    For example:
    - Both of NO2's resonance structures will be kept: [O]N=O <=> O=[N+.][O-]
    - NCO will only have two resonance structures [N.]=C=O <=> N#C[O.], and will loose the third structure which has
      the same octet deviation, has a charge separation, but the radical site has already been considered: [N+.]#C[O-]
    - CH2NO keeps all three structures, since a new radical site is introduced: [CH2.]N=O <=> C=N[O.] <=> C=[N+.][O-]
    However, if the species is not a radical we only keep the structures with the minimal charge span.
    For example:
    - NSH will only keep N#S and not [N-]=[SH+]
    - The following species will loose two thirds of its resonance structures, which are charged: CS(=O)SC <=>
      CS(=O)#SC <=> C[S+]([O-]SC <=> CS([O-])=[S+]C <=> C[S+]([O-])#SC <=> C[S+](=O)=[S-]C
    - The azide structure is know to have three resonance structures: [NH-][N+]#N <=> N=[N+]=[N-] <=> [NH+]#[N+][N-2];
      here we'll lose the third one, which is theoretically "true", but doesn't contribute to reactivity.
    """

    min_charge_span = min(charge_span_list)
    if len(set(charge_span_list)) > 1:
        if filtered_list[0].isRadical():
            # Proceed only if there are structures with different charge spans and the species is a radical.
            charged_list = [filtered_list[index] for index in xrange(len(filtered_list)) if
                            charge_span_list[index] == min_charge_span + 1]  # save the 2nd charge span layer in a temp list
            filtered_list = [filtered_list[index] for index in xrange(len(filtered_list)) if
                            charge_span_list[index] == min_charge_span]  # keep at least one charge span layer
            # Find the radical sites in all filtered_list structures:
            sorting_list = []
            for mol in filtered_list:
                for atom in mol.vertices:
                    if atom.radicalElectrons:
                        sorting_list.append(int(atom.sortingLabel))
            # Find unique radical sites in charged_list and append these structures to filtered_list:
            for mol in charged_list:
                for atom in mol.vertices:
                    if atom.radicalElectrons and int(atom.sortingLabel) not in sorting_list:
                        filtered_list.append(mol)
        else:
            filtered_list = [filtered_list[index] for index in xrange(len(filtered_list)) if
                            charge_span_list[index] == min_charge_span]  # keep one charge span layer for non-radicals

    return filtered_list


def mark_unreactive_structures(filtered_list, mol_list):
    """
    Mark selected structures in filtered_list with the Molecule.reactive flag set to `False` (it is `True` by default)
    if OS and/or adjacent birad unidirectional transitions were used
    Changes the filtered_list object, and does not return anything
    """

    # Special treatment for O=O, S=S, S=O, where the correct ground state structure is being filtered out due to a
    # higher octet deviation
    if is_OS(filtered_list[0]):
        # TODO: should depend on allowSingletO2 option from the input file, in which case the function should work for
        # all OS structures but the singlet
        for mol in filtered_list:
            if is_OS(mol) != 1:  # This is not the triplet ground state [O][O], [S][S], or [S][O]
                mol.reactive = False
        for mol in filtered_list:
            if is_OS(mol) == 1:  # Check whether filtered_list already contains the ground state
                break
        else:  # the ground state was filtered out, append it to filtered_list
            for mol in mol_list:
                if is_OS(mol) == 1:  # Check if mol is the ground state structure, e.g., [O][O] etc.
                    filtered_list.append(mol)   # If so, append mol (as the reactive structure)
                    break
            else:  # could not find the ground state in mol_list
                for mol in mol_list:
                    logging.info('\n')
                    logging.info(mol.toSMILES())
                    logging.info(mol.toAdjacencyList())
                    logging.info('reactive = {0}'.format(mol.reactive))
                raise ValueError('The ground state for species {0} could not be found in its structure'
                                 ' list.'.format(filtered_list[0].toSMILES()))

    # Special treatment for generate_birad_multiple_bond_resonance_structures, where excited states should be marked
    # as unreactive
    else:
        # TODO: should depend on input file option
        min_num_rads = max_num_rads = filtered_list[0].getRadicalCount()
        for mol in filtered_list:
            if mol.getRadicalCount() < min_num_rads:
                min_num_rads = mol.getRadicalCount()
            elif mol.getRadicalCount() > max_num_rads:
                max_num_rads = mol.getRadicalCount()
        if max_num_rads != min_num_rads:  # mark as unreactive only if the species contains structures with
            # different number of radicals
            for mol in filtered_list:
                if mol.getRadicalCount() > min_num_rads:
                    mol.reactive = False
                    logging.debug('marked structure {0} as unreactive'.format(mol.toSMILES()))

    # sort all structures in filtered_list so that the reactive ones are first
    filtered_list.sort(key=lambda mol: mol.reactive, reverse=True)

    # Make sure that the first original structure is always first in the list (unless it was filtered out).
    # Important whenever Species.molecule[0] is expected to be used (e.g., training reactions) after generating the
    # resonance structures. If it was filtered out, it should be appended to the list.
    for index in xrange(len(filtered_list)):
        if filtered_list[index].isIsomorphic(mol_list[0]):
            filtered_list.insert(0, filtered_list.pop(index))
            break
    else:
        # Append the original structure to mol_list and set `reactive` to `False`.
        # This structure may very well deviate from the octet rule or have other attributes by which it should have
        # been filtered out. However, for processing reactions (e.g., degeneracy calculations) it should be kept
        # (e.g., [::N]O <=> [::N][::O.] + [H.], where [::N][::O.] should be recognized as [:N.]=[::O]).
        mol = mol_list[0]
        mol.reactive = False
        filtered_list.append(mol)


def check_reactive(filtered_list):
    """
    Check that there's at least one reactive structure in the returned list.
    If not, raise an error (does not return anything)
    """

    if not any([mol.reactive for mol in filtered_list]):
        logging.info('\n\n')
        logging.error('No reactive structures were attributed to species {0}'.format(filtered_list[0].toSMILES()))
        for mol in filtered_list:
            logging.info('Structure: {0}\n{1}Reactive: {2}'.format(mol.toSMILES(),mol.toAdjacencyList(),mol.reactive))
            logging.info('\n')
        raise AssertionError('Each species must have at least one reactive structure. Something went wrong'
                             ' when processing species {0}'.format(filtered_list[0].toSMILES()))


def is_OS(mol):
    """
    This is a helper function that decides whether a molecule is either O2, S2, or SO,
    and whether it is the ground ([O.][O.], [S.][S.], or [S.][O.]) or the excited (O=O, S=S, or S=O) state.
    returns an integer:
    0 - neither O2, S2, or SO
    1 - triplet ground state ([O.][O.], [S.][S.], or [S.][O.])
    2 - singlet excited state (O=O, S=S, or S=O)
    3 - a O2/S2/SO structure which is neither case `1` or `2` (e.g., [:S..][:S]), and will eventually be filtered out
    """
    if (len(mol.vertices) == 2 and
            ((mol.vertices[0].isOxygen() or mol.vertices[0].isSulfur()) and
            (mol.vertices[1].isOxygen() or mol.vertices[1].isSulfur()))):
        # This is O2/S2/SO
        if (mol.vertices[0].bonds[mol.vertices[1]].isSingle() and
                mol.vertices[0].radicalElectrons == mol.vertices[1].radicalElectrons == 1):
            # This is the ground state
            return 1
        if mol.vertices[0].bonds[mol.vertices[1]].isDouble() and mol.multiplicity == 1:
            # This is the excited state
            return 2
        # this is another O2/S2/SO structure, which will eventually be filtered out
        return 3
    return 0
