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
This module contains classes and functions that are used by multiple modules
in this subpackage.
"""
import itertools
import logging
import json
import re
import os

from rmgpy.data.base import LogicNode
from rmgpy.exceptions import DatabaseError, SpeciesIdentityError
from rmgpy.molecule import Group, Molecule
from rmgpy.molecule.fragment import Fragment
from rmgpy.reaction import Reaction
from rmgpy.species import Species


################################################################################


_EXTERNAL_PROVENANCE_START = '[RMG external library provenance v1] '
_EXTERNAL_PROVENANCE_END = ' [/RMG external library provenance v1]'


def format_external_library_provenance(label, source):
    """Encode a reserved full-line restart record with a semantic label and source."""
    return (_EXTERNAL_PROVENANCE_START
            + json.dumps({'label': label, 'source': os.path.realpath(source)}, sort_keys=True)
            + _EXTERNAL_PROVENANCE_END)


def parse_external_library_provenance(description):
    """Read only reserved full-line records, refusing malformed or conflicting ones."""
    def unique_keys(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise DatabaseError('Conflicting external library provenance JSON keys.')
            result[key] = value
        return result

    provenance = None
    for line in description.splitlines():
        line = line.strip()
        if not line.startswith(_EXTERNAL_PROVENANCE_START.rstrip()):
            continue
        if not line.startswith(_EXTERNAL_PROVENANCE_START) or not line.endswith(_EXTERNAL_PROVENANCE_END):
            raise DatabaseError('Malformed external library provenance record.')
        try:
            record = json.loads(line[len(_EXTERNAL_PROVENANCE_START):-len(_EXTERNAL_PROVENANCE_END)],
                                object_pairs_hook=unique_keys)
        except (TypeError, ValueError) as error:
            raise DatabaseError('Malformed external library provenance JSON.') from error
        if (not isinstance(record, dict) or set(record) != {'label', 'source'}
                or not all(isinstance(value, str) and value for value in record.values())
                or not os.path.isabs(record['source'])):
            raise DatabaseError('Invalid external library provenance fields.')
        record['source'] = os.path.realpath(record['source'])
        if provenance is not None and provenance != record:
            raise DatabaseError('Conflicting external library provenance records.')
        provenance = record
    return provenance


def check_smiles_keyed_efficiencies(efficiencies):
    """Refuse resolved colliders before a library writer reduces keys to SMILES."""
    from rmgpy.export import has_export_state
    for collider in efficiencies:
        if isinstance(collider, Molecule) and has_export_state(collider):
            raise SpeciesIdentityError(
                'Cannot serialize resolved collider as a SMILES-keyed kinetics-library '
                'efficiency. Use a state-aware Chemkin or RMS export instead:\n{0}'.format(
                    collider.to_adjacency_list()))


def validate_generic_collider_dictionary(collider, species_dict):
    """Refuse generic M syntax when it could hide a resolved dictionary collider."""
    if collider.upper().strip() != "(+M)":
        return
    for label in ('M', 'm'):
        species = species_dict.get(label)
        if species is not None and any(mol.has_resolved_state() for mol in species.molecule):
            raise SpeciesIdentityError(
                'Cannot reload resolved collider "{0}" through generic third-body syntax; '
                'M and m are reserved generic collider labels.'.format(label))


def validate_participant_labels(species):
    """Refuse resolved labels that the library reaction grammar splits at plus."""
    for spc in species:
        if (isinstance(spc, Species) and '+' in spc.label
                and any(mol.has_resolved_state() for mol in spc.molecule)):
            raise SpeciesIdentityError(
                'Cannot serialize or reload resolved participant "{0}"; '
                '+ is an ambiguous reaction separator in library labels.'.format(spc.label))


def validate_dictionary_participant_labels(reactants, products, species_dict):
    """Check whole dictionary labels before tokenization can discard a state."""
    import re
    referenced = [spc for label, spc in species_dict.items() if '+' in label
                  and any(mol.has_resolved_state() for mol in spc.molecule)
                  and any(re.search(r'(?:^|\+)\s*' + re.escape(label) + r'\s*(?:\+|$)', side)
                          for side in (reactants, products))]
    validate_participant_labels(referenced)


def library_reaction_equation(entry, declarations=None):
    """Validate an equation against its emitted dictionary, or refuse missing context."""
    from collections import Counter
    from rmgpy.export import SpeciesReferences, describe_species, refuse_resolved_species, resolve_species_reference
    from rmgpy.util import get_reaction_collider, parse_reaction_equation
    reaction = entry.item
    references = reaction.reactants + reaction.products
    if reaction.specific_collider is not None:
        references = references + [reaction.specific_collider]
    references += list(getattr(entry.data, 'coverage_dependence', None) or {})
    if declarations is None:
        refuse_resolved_species(references, 'kinetics equation without emitted dictionary declarations')
        declarations = SpeciesReferences(list(dict.fromkeys(references)), context='kinetics library',
                                         allow_ground_collisions=True)
    elif not isinstance(declarations, SpeciesReferences):
        declarations = SpeciesReferences(declarations, context='kinetics library', allow_ground_collisions=True)
    reactants = [resolve_species_reference(spc, declarations) for spc in reaction.reactants]
    products = [resolve_species_reference(spc, declarations) for spc in reaction.products]
    for spc in (getattr(entry.data, 'coverage_dependence', None) or {}):
        resolve_species_reference(spc, declarations)
    collider = resolve_species_reference(reaction.specific_collider, declarations) if reaction.specific_collider else None
    suffix = ' (+{0})'.format(collider) if collider is not None else ''
    exact_tokens = list(zip(reactants, reaction.reactants)) + list(zip(products, reaction.products))
    if collider is not None:
        exact_tokens.append((collider, reaction.specific_collider))
    tokens = exact_tokens
    if entry.label:
        left, right, reversible = parse_reaction_equation(entry.label)
        left_collider, right_collider = get_reaction_collider(left), get_reaction_collider(right)
        generic = left_collider is not None and left_collider.upper() == '(+M)'
        if generic and collider is None and hasattr(entry.data, 'efficiencies'):
            suffix = ' (+M)'
        expected_collider = suffix.strip() or None
        if left_collider:
            left = left.replace(left_collider, '', 1)
        if right_collider:
            right = right.replace(right_collider, '', 1)
        # Reaction.__str__ includes assigned RMG indices. They are display
        # annotations, not different library species; only the object's actual
        # index may annotate a resolved declaration name during validation.
        indexed_reactants = [name + ('({0})'.format(spc.index) if spc.index >= 0 else '')
                             for name, spc in zip(reactants, reaction.reactants)]
        indexed_products = [name + ('({0})'.format(spc.index) if spc.index >= 0 else '')
                            for name, spc in zip(products, reaction.products)]
        left_names = Counter(token.strip() for token in left.split('+'))
        right_names = Counter(token.strip() for token in right.split('+'))
        if (reversible != reaction.reversible or left_collider != expected_collider
                or right_collider != expected_collider
                or left_names not in (Counter(reactants), Counter(indexed_reactants))
                or right_names not in (Counter(products), Counter(indexed_products))):
            raise SpeciesIdentityError('Entry label equation {0!r} disagrees with its reaction identities.'.format(entry.label))
        # Counters permit reordered terms; bind tokens to the corresponding
        # spelling of the typed participants, rather than their written order.
        tokens = []
        for species, exact, indexed, names in (
                (reaction.reactants, reactants, indexed_reactants, left_names),
                (reaction.products, products, indexed_products, right_names)):
            tokens.extend(zip(exact if names == Counter(exact) else indexed, species))
        if collider is not None:
            tokens.append((collider, reaction.specific_collider))

    canonical = (' + '.join(reactants) + suffix + (' <=> ' if reaction.reversible else ' => ')
                 + ' + '.join(products) + suffix)
    if not any(mol.has_resolved_state() for spc in declarations for mol in spc.molecule):
        return entry.label or canonical

    # Use the reader's exact-key-first rule against the dictionary we emit.
    # A display annotation is safe only when it reloads as the intended full
    # identity; canonical names are checked by the same rule before writing.
    for equation, bindings in ((entry.label or canonical, tokens), (canonical, exact_tokens)):
        dictionary = {}
        for name, species in zip(declarations.names, declarations):
            dictionary.setdefault(name, species)  # render_dictionary writes the first record for each name.
        add_dictionary_index_aliases(equation, dictionary)
        mismatches = [(token, intended, dictionary.get(token)) for token, intended in bindings
                      if token not in dictionary or not dictionary[token].is_isomorphic(intended)]
        if not mismatches:
            return equation
    token, intended, selected = mismatches[0]
    raise SpeciesIdentityError(
        'Kinetics equation token "{0}" cannot unambiguously reload {1}; dictionary selects {2}.'.format(
            token, describe_species(intended),
            describe_species(selected) if selected is not None else 'no identity'))


class _SerializedCoverageReference:
    """An already-resolved key for the kinetics classes' repr protocol."""

    def __init__(self, identifier):
        self.identifier = identifier

    def to_chemkin(self):
        return self.identifier


def library_serializable_kinetics(kinetics, declarations):
    """Resolve repr keys on a copy, leaving live kinetics untouched."""
    import copy
    from rmgpy.export import SpeciesReferences, resolve_species_reference
    result = copy.copy(kinetics)
    if hasattr(result, 'arrhenius'):
        result.arrhenius = [library_serializable_kinetics(rate, declarations) for rate in result.arrhenius]
    if hasattr(result, 'efficiencies'):
        check_smiles_keyed_efficiencies(result.efficiencies)
        efficiencies = {}
        for molecule, value in result.efficiencies.items():
            molecule = Molecule(smiles=molecule) if isinstance(molecule, str) else molecule
            # The SMILES field explicitly declares this unresolved molecular
            # identity inline; it need not be a reaction dictionary participant.
            inline = Species(molecule=[molecule])
            inline_declarations = SpeciesReferences([inline], lambda spc: spc.molecule[0].to_smiles(),
                                                   context='library efficiency')
            efficiencies[resolve_species_reference(molecule, inline_declarations)] = value
        result.efficiencies = dict(sorted(efficiencies.items()))
    coverage = getattr(result, 'coverage_dependence', None) or {}
    if coverage:
        from rmgpy.export import has_export_state
        resolved = any(has_export_state(species) for species in declarations)
        # Ground-only repr uses its legacy indexed declaration names. References
        # still validate against the reaction inventory and use the same resolver.
        coverage_names = declarations if resolved else SpeciesReferences(
            coverage, lambda species: species.to_chemkin(), context='ground library coverage')
        serialized = {}
        for spc, parameters in coverage.items():
            identifier = resolve_species_reference(spc, declarations)
            if not resolved:
                identifier = resolve_species_reference(spc, coverage_names)
            serialized[_SerializedCoverageReference(identifier)] = parameters
        result.coverage_dependence = serialized
    return result


def save_entry(f, entry, declarations=None):
    """
    Save an `entry` in the kinetics database by writing a string to
    the given file object `f`.
    """

    if isinstance(entry.item, Reaction):
        entry.item.check_resolved_species_reversibility(kinetics=entry.data)
        validate_participant_labels(entry.item.reactants + entry.item.products)

    def sort_efficiencies(efficiencies0):
        efficiencies = {}
        for mol, eff in efficiencies0.items():
            if isinstance(mol, str):
                # already in SMILES string format
                smiles = mol
            else:
                smiles = mol.to_smiles()

    collider = getattr(entry.item, 'specific_collider', None)
    if (collider is not None and collider.label.strip().upper() == 'M'
            and any(mol.has_resolved_state() for mol in collider.molecule)):
        raise SpeciesIdentityError(
            'Cannot serialize resolved named collider "{0}"; '
            'M and m are reserved generic collider labels.'.format(collider.label))

    if hasattr(entry.data, 'efficiencies'):
        check_smiles_keyed_efficiencies(entry.data.efficiencies)

    if declarations is None:
        from rmgpy.export import kinetics_references, refuse_resolved_species
        refuse_resolved_species(kinetics_references(entry.data),
                                'kinetics entry writer without emitted dictionary declarations',
                                reactions=[entry.item] if isinstance(entry.item, Reaction) else [])

    label = entry.label
    if (isinstance(entry.item, Reaction) and entry.item.reactants
            and all(isinstance(spc, Species) for spc in entry.item.reactants + entry.item.products)):
        label = library_reaction_equation(entry, declarations)

    f.write('entry(\n')
    f.write('    index = {0:d},\n'.format(entry.index))
    if label != '':
        f.write('    label = "{0}",\n'.format(label))

    # Entries for kinetic rules, libraries, training reactions
    # and depositories will have a Reaction object for its item
    if isinstance(entry.item, Reaction):
        # Write out additional data if depository or library
        # kinetic rules would have a Group object for its reactants instead of Species
        if isinstance(entry.item.reactants[0], Species):
            # Add degeneracy if the reaction is coming from a depository or kinetics library
            f.write('    degeneracy = {0:.1f},\n'.format(entry.item.degeneracy))
            if entry.item.duplicate:
                f.write('    duplicate = {0!r},\n'.format(entry.item.duplicate))
            if not entry.item.reversible:
                f.write('    reversible = {0!r},\n'.format(entry.item.reversible))
            if entry.item.allow_pdep_route:
                f.write('    allow_pdep_route = {0!r},\n'.format(entry.item.allow_pdep_route))
            if entry.item.elementary_high_p:
                f.write('    elementary_high_p = {0!r},\n'.format(entry.item.elementary_high_p))
            if entry.item.allow_max_rate_violation:
                f.write('    allow_max_rate_violation = {0!r},\n'.format(entry.item.allow_max_rate_violation))
    # Entries for groups with have a group or logicNode for its item
    elif isinstance(entry.item, Group):
        f.write('    group = \n')
        f.write('"""\n')
        f.write(entry.item.to_adjacency_list())
        f.write('""",\n')
    elif isinstance(entry.item, LogicNode):
        f.write('    group = "{0}",\n'.format(entry.item))
    else:
        raise DatabaseError("Encountered unexpected item of type {0} while "
                            "saving database.".format(entry.item.__class__))

    # Write kinetics
    if isinstance(entry.data, str):
        f.write('    kinetics = "{0}",\n'.format(entry.data))
    elif entry.data is not None:
        if declarations is None:
            from rmgpy.export import SpeciesReferences
            declarations = SpeciesReferences(
                [spc for spc in getattr(entry.item, 'reactants', []) + getattr(entry.item, 'products', [])
                 if isinstance(spc, Species)], context='kinetics library', allow_ground_collisions=True)
        serialized = library_serializable_kinetics(entry.data, declarations)
        kinetics = repr(serialized)
        kinetics = '    kinetics = {0},\n'.format(kinetics.replace('\n', '\n    '))
        f.write(kinetics)
    else:
        f.write('    kinetics = None,\n')

    # Write reference
    if entry.reference is not None:
        reference = entry.reference.to_pretty_repr()
        lines = reference.splitlines()
        f.write('    reference = {0}\n'.format(lines[0]))
        for line in lines[1:-1]:
            f.write('    {0}\n'.format(line))
        f.write('    ),\n'.format(lines[0]))

    if entry.reference_type != "":
        f.write('    referenceType = "{0}",\n'.format(entry.reference_type))
    if entry.rank is not None:
        f.write('    rank = {0},\n'.format(entry.rank))

    if entry.short_desc.strip() != '':
        f.write(f'    shortDesc = """{entry.short_desc.strip()}""",\n')
    if entry.long_desc.strip() != '':
        if parse_external_library_provenance(entry.long_desc):
            # repr preserves JSON escapes through the Python library-file reader.
            f.write('    longDesc = {!r},\n'.format(entry.long_desc.strip()))
        else:
            f.write(f'    longDesc = \n"""\n{entry.long_desc.strip()}\n""",\n')

    # write metal attributes
    if entry.metal:
        f.write('    metal = "{0}",\n'.format(entry.metal))
    if entry.facet:
        f.write('    facet = "{0}",\n'.format(entry.facet))
    if entry.site:
        f.write('    site = "{0}",\n'.format(entry.site))

    f.write(')\n\n')


def get_molecularity(reaction):
    """
    Return the molecularity of `reaction`, i.e. the number of particles that collide to form the
    transition state, which is the order of its rate coefficient and therefore determines the
    dimensions of its A factor.

    This is not simply ``len(reaction.reactants)``. Electrons transferred by a reaction are genuine
    reactant particles, yet they are never listed among the reactants: they are carried by the
    ``electrons`` attribute, which reaction families declare in their ``groups.py`` and which is
    negative when electrons are consumed. Non-dissociative electron attachment, ``A + e- => A-``,
    is stored with a single reactant but is bimolecular, and its rate is expressed in
    cm^3/(mol*s) accordingly. Electrons released by a reaction are products and do not contribute.

    Returns an ``int``.
    """
    electrons = getattr(reaction, 'electrons', 0) or 0
    # Only electrons on the reactant side (a negative count) add to the molecularity.
    return len(reaction.reactants) + max(-electrons, 0)


def ensure_species(input_list, resonance=False, keep_isomorphic=False):
    """
    The input list of :class:`Species` or :class:`Molecule` objects is modified
    in place to only have :class:`Species` objects. Returns None.
    """
    for index, item in enumerate(input_list):
        if isinstance(item, Molecule) or isinstance(item, Fragment):
            new_item = Species(molecule=[item])
        elif isinstance(item, Species):
            new_item = item
        else:
            raise TypeError('Only Molecule or Species objects can be handled.')
        if resonance:
            if not any([mol.reactive for mol in new_item.molecule]):
                # if generating a reaction containing a Molecule with a reactive=False flag (e.g., for degeneracy
                # calculations), that was now converted into a Species, first mark as reactive=True
                new_item.molecule[0].reactive = True
            new_item.generate_resonance_structures(keep_isomorphic=keep_isomorphic)
        input_list[index] = new_item


def generate_molecule_combos(input_species):
    """
    Generate combinations of molecules from the given species objects.
    """
    if len(input_species) == 1:
        combos = [(mol,) for mol in input_species[0].molecule]
    elif len(input_species) == 2:
        combos = itertools.product(input_species[0].molecule, input_species[1].molecule)
    elif len(input_species) == 3:
        combos = itertools.product(input_species[0].molecule, input_species[1].molecule, input_species[2].molecule)
    else:
        raise ValueError('Reaction generation can be done for 1, 2, or 3 species, not {0}.'.format(len(input_species)))

    return combos


def ensure_independent_atom_ids(input_species, resonance=True):
    """
    Given a list or tuple of :class:`Species` or :class:`Molecule` objects,
    ensure that atom ids are independent.
    The `resonance` argument can be set to False to not generate
    resonance structures.

    Modifies the list in place (replacing :class:`Molecule` with :class:`Species`).
    Returns None.
    """
    ensure_species(input_species)  # do not generate resonance structures since we do so below

    # Inline check that all atom IDs across all species' first molecule are unique.
    # Building a list then converting to set is faster than incremental set.add()
    # because set(list) has a vectorized C path.
    ids = []
    for spcs in input_species:
        atoms = spcs.molecule[0].atoms
        ids.extend([atom.id for atom in atoms])

    if len(set(ids)) != len(ids):
        # Collision: reassign IDs and remake resonance structures
        for species in input_species:
            reactive_mols = []
            unreactive_mols = []
            for m in species.molecule:
                (reactive_mols if m.reactive else unreactive_mols).append(m)
            mol = reactive_mols[0]  # Choose first reactive molecule
            mol.assign_atom_ids()
            species.molecule = [mol]
            # Remake resonance structures with new labels
            if resonance:
                species.generate_resonance_structures(keep_isomorphic=True)
            if unreactive_mols:
                species.molecule.extend(unreactive_mols)
    elif resonance:
        # IDs are already independent, generate resonance structures if needed
        for species in input_species:
            species.generate_resonance_structures(keep_isomorphic=True)

def check_for_same_reactants(reactants):
    """
    Given a list reactants, check if the reactants are the same.
    If they refer to the same memory address, then make a deep copy so they can be manipulated independently.
    
    Returns a tuple containing the modified reactants list, and an integer containing the number of identical reactants in the reactants list. 
    
    """

    same_reactants = 0
    if len(reactants) == 2:
        if reactants[0] is reactants[1]:
            reactants[1] = reactants[1].copy(deep=True)
            same_reactants = 2
        elif reactants[0].is_isomorphic(reactants[1]):
            same_reactants = 2
    elif len(reactants) == 3:
        same_01 = reactants[0] is reactants[1]
        same_02 = reactants[0] is reactants[2]
        if same_01 and same_02:
            same_reactants = 3
            reactants[1] = reactants[1].copy(deep=True)
            reactants[2] = reactants[2].copy(deep=True)
        elif same_01:
            same_reactants = 2
            reactants[1] = reactants[1].copy(deep=True)
        elif same_02:
            same_reactants = 2
            reactants[2] = reactants[2].copy(deep=True)
        elif reactants[1] is reactants[2]:
            same_reactants = 2
            reactants[2] = reactants[2].copy(deep=True)
        else:
            same_01 = reactants[0].is_isomorphic(reactants[1])
            same_02 = reactants[0].is_isomorphic(reactants[2])
            if same_01 and same_02:
                same_reactants = 3
            elif same_01 or same_02:
                same_reactants = 2
            elif reactants[1].is_isomorphic(reactants[2]):
                same_reactants = 2
    elif len(reactants) > 3:
        raise ValueError('Cannot check for duplicate reactants if provided number of reactants is greater than 3. ' 
                         'Got: {} reactants'.format(len(reactants))) 
        
    return reactants, same_reactants

def find_degenerate_reactions(rxn_list, same_reactants=None, template=None, kinetics_database=None,
                              kinetics_family=None, save_order=False, resonance=True):
    """
    Given a list of Reaction objects, this method combines degenerate
    reactions and increments the reaction degeneracy value. For multiple
    transition states, this method keeps them as duplicate reactions.

    If a template is specified, then the reaction list will be filtered
    to leave only reactions which match the specified template, then the
    degeneracy will be calculated as usual.

    A KineticsDatabase or KineticsFamily instance can also be provided to
    calculate the degeneracy for reactions generated in the reverse direction.
    If not provided, then it will be retrieved from the global database.

    This algorithm used to exist in family._generate_reactions, but was moved
    here so it could operate across reaction families.

    This method returns an updated list with degenerate reactions removed.

    Args:
        rxn_list (list):                                reactions to be analyzed
        same_reactants (bool, optional):                indicate whether the reactants are identical
        template (list, optional):                      specify a specific template to filter by
        kinetics_database (KineticsDatabase, optional): provide a KineticsDatabase instance for calculating degeneracy
        kinetics_family (KineticsFamily, optional):     provide a KineticsFamily instance for calculating degeneracy
        save_order (bool, optional):                    reset atom order after performing atom isomorphism
        resonance (bool, optional):                     whether to consider resonance when computing degeneracy 

    Returns:
        Reaction list with degenerate reactions combined with proper degeneracy values
    """
    # If a specific reaction template is requested, filter by that template
    if template is not None:
        selected_rxns = []
        template = frozenset(template)
        for rxn in rxn_list:
            if template == frozenset(rxn.template):
                selected_rxns.append(rxn)
        if not selected_rxns:
            # Only log a warning here. If a non-empty output is expected, then the caller should raise an exception
            logging.warning('No reactions matched the specified template, {0}'.format(template))
            return []
    else:
        selected_rxns = rxn_list

    # We want to sort all the reactions into sublists composed of isomorphic reactions
    # with degenerate transition states
    sorted_rxns = []
    for rxn0 in selected_rxns:
        rxn0.ensure_species(save_order=save_order)
        if len(sorted_rxns) == 0:
            # This is the first reaction, so create a new sublist
            sorted_rxns.append([rxn0])
        else:
            # Loop through each sublist, which represents a unique reaction
            for sub_list in sorted_rxns:
                # Try to determine if the current rxn0 is identical or isomorphic to any reactions in the sublist
                isomorphic = False
                identical = False
                same_template = True
                for rxn in sub_list:
                    isomorphic = rxn0.is_same_reaction(rxn, check_identical=False, strict=False,
                                                    check_template_rxn_products=True, save_order=save_order)
                    if isomorphic:
                        identical = rxn0.is_same_reaction(rxn, check_identical=True, strict=False,
                                                       check_template_rxn_products=True, save_order=save_order)
                        if identical:
                            # An exact copy of rxn0 is already in our list, so we can move on
                            break
                        same_template = frozenset(rxn.template) == frozenset(rxn0.template)
                    else:
                        # This sublist contains a different product
                        break

                # Process the reaction depending on the results of the comparisons
                if identical:
                    # This reaction does not contribute to degeneracy
                    break
                elif isomorphic:
                    if same_template:
                        # We found the right sublist, and there is no identical reaction
                        # We should add rxn0 to the sublist as a degenerate rxn, and move on to the next rxn
                        sub_list.append(rxn0)
                        break
                    else:
                        # We found an isomorphic sublist, but the reaction templates are different
                        # We need to mark this as a duplicate and continue searching the remaining sublists
                        rxn0.duplicate = True
                        sub_list[0].duplicate = True
                        continue
                else:
                    # This is not an isomorphic sublist, so we need to continue searching the remaining sublists
                    # Note: This else statement is not technically necessary but is included for clarity
                    continue
            else:
                # We did not break, which means that there was no isomorphic sublist, so create a new one
                sorted_rxns.append([rxn0])

    rxn_list = []
    for sub_list in sorted_rxns:
        # Collapse our sorted reaction list by taking one reaction from each sublist
        rxn = sub_list[0]
        # The degeneracy of each reaction is the number of reactions that were in the sublist
        rxn.degeneracy = sum([reaction0.degeneracy for reaction0 in sub_list])
        rxn_list.append(rxn)

    for rxn in rxn_list:
        if rxn.is_forward:
            reduce_same_reactant_degeneracy(rxn, same_reactants)
        else:
            # fix the degeneracy of (not ownReverse) reactions found in the backwards direction
            try:
                family = kinetics_family or kinetics_database.families[rxn.family]
            except AttributeError:
                from rmgpy.data.rmg import get_db
                family = get_db('kinetics').families[rxn.family]
            if not family.own_reverse:
                rxn.degeneracy = family.calculate_degeneracy(rxn, resonance=resonance)

    return rxn_list


def reduce_same_reactant_degeneracy(reaction, same_reactants=None):
    """
    This method reduces the degeneracy of reactions with identical reactants,
    since translational component of the transition states are already taken
    into account (so swapping the same reactant is not valid)

    same_reactants can be None or an integer. If it is None, then isomorphism
    checks will be done to determine if the reactions are the same. If it is an
    integer, that integer denotes the number of reactants that are isomorphic.

    This comes from work by Bishop and Laidler in 1965
    """
    if not (same_reactants == 0 or same_reactants == 1):
        if len(reaction.reactants) == 2:
            if ((reaction.is_forward and same_reactants == 2) or
                    reaction.reactants[0].is_isomorphic(reaction.reactants[1])):
                reaction.degeneracy *= 0.5
                logging.debug(
                    'Degeneracy of reaction {} was decreased by 50% to {} since the reactants are identical'.format(
                        reaction, reaction.degeneracy)
                )
        elif len(reaction.reactants) == 3:
            if reaction.is_forward:
                if same_reactants == 3:
                    reaction.degeneracy /= 6.0
                    logging.debug(
                        'Degeneracy of reaction {} was divided by 6 to give {} since all of the reactants '
                        'are identical'.format(reaction, reaction.degeneracy)
                    )
                elif same_reactants == 2:
                    reaction.degeneracy *= 0.5
                    logging.debug(
                        'Degeneracy of reaction {} was decreased by 50% to {} since two of the reactants '
                        'are identical'.format(reaction, reaction.degeneracy)
                    )
            else:
                same_01 = reaction.reactants[0].is_isomorphic(reaction.reactants[1])
                same_02 = reaction.reactants[0].is_isomorphic(reaction.reactants[2])
                if same_01 and same_02:
                    reaction.degeneracy /= 6.0
                    logging.debug(
                        'Degeneracy of reaction {} was divided by 6 to give {} since all of the reactants '
                        'are identical'.format(reaction, reaction.degeneracy)
                    )
                elif same_01 or same_02:
                    reaction.degeneracy *= 0.5
                    logging.debug(
                        'Degeneracy of reaction {} was decreased by 50% to {} since two of the reactants '
                        'are identical'.format(reaction, reaction.degeneracy)
                    )
                elif reaction.reactants[1].is_isomorphic(reaction.reactants[2]):
                    reaction.degeneracy *= 0.5
                    logging.debug(
                        'Degeneracy of reaction {} was decreased by 50% to {} since two of the reactants '
                        'are identical'.format(reaction, reaction.degeneracy)
                    )


def add_dictionary_index_aliases(equation, species_dict):
    """Add missing display-index aliases only for resolved dictionary species."""
    from rmgpy.util import get_reaction_collider, parse_reaction_equation
    left, right, _ = parse_reaction_equation(equation)
    names = []
    for side in (left, right):
        collider = get_reaction_collider(side)
        if collider:
            if collider.upper() != '(+M)':
                names.append(collider[2:-1])
            side = side.replace(collider, '', 1)
        names.extend(token.strip() for token in side.split('+'))
    for name in names:
        if name in species_dict:
            continue
        match = re.fullmatch(r'(.+)\(\d+\)', name)
        if match and match.group(1) in species_dict:
            species = species_dict[match.group(1)]
            if any(mol.has_resolved_state() for mol in species.molecule):
                species_dict[name] = species
