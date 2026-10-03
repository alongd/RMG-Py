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

"""State-aware declarations and reference resolution for mechanism writers."""

from rmgpy.exceptions import SpeciesIdentityError
from rmgpy.molecule import Molecule
from rmgpy.species import Species


def export_molecule(reference, molecule=None):
    """Represent a declared manifold member as explicit v=0 without mutating it."""
    if isinstance(reference, Species):
        reference._molecule_state_key()
    if molecule is None:
        molecule = reference.molecule[0] if isinstance(reference, Species) else reference
    declaration = (getattr(reference, 'props', {}).get('vibrational_manifold')
                   or molecule.props.get('vibrational_manifold'))
    if not declaration:
        return molecule
    if molecule.vibrational_level not in (-1, 0):
        raise SpeciesIdentityError('A vibrationalManifold declaration cannot export a nonzero fixed level.')
    molecule = molecule.copy(deep=True)
    molecule.props.pop('vibrational_manifold', None)
    molecule.vibrational_level = 0
    return molecule


def has_export_state(reference):
    """Include the declared v=0 member at state-aware or refusing export boundaries."""
    molecules = reference.molecule if isinstance(reference, Species) else [reference]
    return any(export_molecule(reference, molecule).has_resolved_state() for molecule in molecules)


class SpeciesReferences(list):
    """A writer's declared species and their emitted names.

    The list interface keeps existing kinetics APIs usable. Names are allocated
    by the writer once, and lookup and identity checks are shared by all writers.
    Treat this declaration set as immutable while serializing a mechanism.
    ``allow_ground_collisions`` preserves writers' historical ground-only
    behavior for conflicting labels; resolved collisions still refuse.
    """

    def __init__(self, species, identifiers=None, error_type=SpeciesIdentityError, context='export',
                 collision_type=None, allow_ground_collisions=True):
        super().__init__(species)
        self.error_type = error_type
        self.context = context
        self.identifier = identifiers if callable(identifiers) else None
        self.names = tuple(identifiers(spc) for spc in self) if callable(identifiers) else tuple(
            identifiers if identifiers is not None else (spc.label for spc in self))
        if len(self.names) != len(self):
            raise ValueError('One emitted identifier is required per declared species.')
        self.positions = {id(spc): index for index, spc in enumerate(self)}
        self.by_name = {}
        for index, (spc, name) in enumerate(zip(self, self.names)):
            previous = self.by_name.get(name)
            if previous is not None and self[previous] is not spc:
                resolved = any(has_export_state(species) for species in (self[previous], spc))
                if not allow_ground_collisions or resolved:
                    same = any(export_molecule(self[previous], first).is_isomorphic(export_molecule(spc, second))
                               for first in self[previous].molecule for second in spc.molecule)
                    if not same:
                        raise (collision_type or error_type)(
                            '{0} identifier "{1}" collides under a shared label between {2} and {3}.'.format(
                                context, name, describe_species(self[previous]), describe_species(spc)))
            self.by_name[name] = index


def describe_species(reference):
    """Describe identity-bearing metadata without using SMILES as an identity."""
    molecules = reference.molecule if isinstance(reference, Species) else [reference]
    states = [export_molecule(reference, mol).state_suffix() for mol in molecules]
    return '"{0}" (index={1}, state={2!r})'.format(
        getattr(reference, 'label', 'molecular reference'), getattr(reference, 'index', -1), states)


def resolve_species_reference(reference, declarations, identifiers=None, allow_missing_efficiency=False):
    """Return the declared writer identifier for a full molecular identity.

    Accepts Species or Molecule references. Object membership uses the same
    declaration as the output record; other objects must match molecular identity,
    including electronic and vibrational state. Missing identities and colliding
    declarations raise a named error instead of falling back to a matching label.
    ``declarations`` can be a SpeciesReferences context or a species sequence.
    For an efficiency only, ``allow_missing_efficiency`` retains the legacy
    omission of an absent ground-state collider. An absent resolved collider
    still raises; every emitted efficiency must resolve to a declaration.
    """
    if not isinstance(declarations, SpeciesReferences):
        declarations = SpeciesReferences(declarations, identifiers)
    position = declarations.positions.get(id(reference))
    if position is not None:
        return declarations.names[position]
    if not isinstance(reference, (Species, Molecule)):
        raise declarations.error_type('{0} reference has no molecular identity: {1!r}.'.format(
            declarations.context, reference))
    for index, spc in enumerate(declarations):
        if spc.molecule and any(export_molecule(spc, mol).is_isomorphic(export_molecule(reference, other))
                                for mol in spc.molecule
                                for other in (reference.molecule if isinstance(reference, Species) else [reference])):
            return declarations.names[index]
    molecules = reference.molecule if isinstance(reference, Species) else [reference]
    if allow_missing_efficiency and not has_export_state(reference):
        return None
    name = declarations.identifier(reference) if (
        isinstance(reference, Species) and declarations.identifier is not None) else getattr(reference, 'label', '')
    collision = ' (collides with a declared identifier)' if name in declarations.by_name else ''
    raise declarations.error_type(
        '{0} references undeclared identity {1}{2}.'.format(
            declarations.context, describe_species(reference), collision))


def kinetics_references(kinetics, with_kind=False):
    """Yield every efficiency and coverage reference, including wrapped rates."""
    pending = [kinetics]
    while pending:
        rate = pending.pop()
        if rate is None:
            continue
        pending.extend(getattr(rate, 'arrhenius', None) or [])
        pending.extend([getattr(rate, 'arrheniusLow', None), getattr(rate, 'arrheniusHigh', None)])
        for reference in getattr(rate, 'efficiencies', {}):
            reference = Molecule(smiles=reference) if isinstance(reference, str) else reference
            yield (reference, True) if with_kind else reference
        for reference in (getattr(rate, 'coverage_dependence', None) or {}):
            yield (reference, False) if with_kind else reference


def validate_reaction_references(reactions, declarations):
    """Resolve every reference before native conversion or writer filtering."""
    for reaction in reactions:
        for reference in reaction.reactants + reaction.products:
            resolve_species_reference(reference, declarations)
        if reaction.specific_collider is not None:
            resolve_species_reference(reaction.specific_collider, declarations)
        for reference, efficiency in kinetics_references(reaction.kinetics, with_kind=True):
            resolve_species_reference(reference, declarations, allow_missing_efficiency=efficiency)


def refuse_resolved_species(species, context, reactions=(), allow_manifold=False):
    """Refuse resolved identity at a boundary that only carries ground labels."""
    references = list(species)
    for reaction in reactions:
        references.extend(reaction.reactants + reaction.products)
        if reaction.specific_collider is not None:
            references.append(reaction.specific_collider)
        references.extend(kinetics_references(reaction.kinetics))
    for reference in references:
        if isinstance(reference, Species):
            molecules = list(reference.molecule or [])
        else:
            molecules = [reference] if isinstance(reference, Molecule) else []
        for adjacency in (getattr(getattr(reference, 'thermo', None), 'thermo_coverage_dependence', None) or {}):
            molecules = molecules + [Molecule().from_adjacency_list(adjacency)]
        if (any(mol.has_resolved_state() for mol in molecules)
                or (not allow_manifold and (getattr(reference, 'props', {}).get('vibrational_manifold')
                    or any(mol.props.get('vibrational_manifold') for mol in molecules)))):
            raise SpeciesIdentityError(
                'Resolved electronic or vibrational states are not supported by {0}: {1}.'.format(
                    context, describe_species(reference)))


def refuse_resolved_job(job, context):
    """Inspect Arkane species, reaction, network and explorer jobs before output."""
    species = [job.species] if hasattr(job, 'species') else []
    reactions = [job.reaction] if hasattr(job, 'reaction') else []
    network = getattr(job, 'network', None)
    if network is not None:
        species.extend(network.get_all_species())
        reactions.extend(network.path_reactions + network.net_reactions)
    species.extend(getattr(job, 'source', None) or [])
    refuse_resolved_species(species, context, reactions)
    if getattr(job, 'pdepjob', None) is not None:
        refuse_resolved_job(job.pdepjob, context)
