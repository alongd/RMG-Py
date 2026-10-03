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
Contains classes for providing pressure-dependent kinetics estimation
functionality to RMG.
"""

import logging
from functools import lru_cache, wraps
from types import GetSetDescriptorType, MemberDescriptorType
import rmgpy.kinetics as kinetics
from rmgpy.kinetics.uncertainties import RateUncertainty
from rmgpy.data.solvation import SoluteData, SoluteTSData
import os.path
import shutil

import mpmath as mp
import numpy as np
import scipy.optimize as opt

import rmgpy.pdep.network
import rmgpy.reaction
from rmgpy.constants import R
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.data.base import Database, Entry
from rmgpy.electron_balance import get_placement_owner
from rmgpy.exceptions import PressureDependenceError, NetworkError, InvalidMicrocanonicalRateError
from rmgpy.pdep import Configuration
from rmgpy.rmg.react import react_species
from rmgpy.statmech import (Conformer, Rotation, LinearRotor, NonlinearRotor, KRotor,
    SphericalTopRotor, Translation, IdealGasTranslation, Torsion, HinderedRotor,
    Vibration, HarmonicOscillator)
from rmgpy.statmech.mode import Mode
from rmgpy.statmech.torsion import FreeRotor
from rmgpy.species import Species, TransitionState
from rmgpy.molecule import Molecule, Group, Atom
from rmgpy.quantity import ScalarQuantity, ArrayQuantity
from rmgpy.thermo import ThermoData, Wilhoit, NASA, NASAPolynomial
from rmgpy.thermo.model import HeatCapacityModel
from rmgpy.transport import TransportData
from rmgpy.pdep.collision import SingleExponentialDown
from rmgpy.molecule.molecule import Bond
from rmgpy.molecule.graph import Graph, Vertex, Edge
from rmgpy.molecule.element import Element
from rmgpy.molecule.atomtype import AtomType
from rmgpy.data.vaporLiquidMassTransfer import (HenryLawConstantData,
    LiquidVolumetricMassTransferCoefficientData)


################################################################################

class _UnclassifiedSpeciesReference(ValueError):
    """An unsupported type or material position cannot prove electron absence."""


# Membership is deliberate: adding a shipped kinetics class requires reviewing
# its species-bearing fields. The completeness regression audits this list.
_KINETICS_TYPES = frozenset((
    kinetics.KineticsModel, kinetics.PDepKineticsModel, kinetics.TunnelingModel,
    kinetics.Arrhenius, kinetics.ArrheniusEP, kinetics.ArrheniusBM,
    kinetics.ArrheniusChargeTransfer, kinetics.ArrheniusChargeTransferBM,
    kinetics.Marcus, kinetics.TwoTemperaturePlasma, kinetics.ElectronCollisionPlasma,
    kinetics.BadnellRRArrhenius, kinetics.VoronovEIArrhenius,
    kinetics.PDepArrhenius, kinetics.MultiArrhenius, kinetics.MultiPDepArrhenius,
    kinetics.Chebyshev, kinetics.ThirdBody, kinetics.Lindemann, kinetics.Troe,
    kinetics.KineticsData, kinetics.PDepKineticsData, kinetics.Wigner, kinetics.Eckart,
    kinetics.SurfaceArrhenius, kinetics.SurfaceArrheniusBEP,
    kinetics.StickingCoefficient, kinetics.StickingCoefficientBEP,
    kinetics.SurfaceChargeTransfer, kinetics.SurfaceChargeTransferBEP,
    RateUncertainty,
))
_KINETICS_REFERENCE_FIELDS = {
    cls: tuple(name for name in (
        'efficiencies', '_coverage_dependence', 'coverage_dependence', 'highPlimit', 'arrhenius',
        'arrheniusLow', 'arrheniusHigh', 'solute', 'uncertainty',
    ) if hasattr(cls, name) and not (name == 'coverage_dependence'
        and hasattr(cls, '_coverage_dependence')))
    for cls in _KINETICS_TYPES
}
# These descriptors belong only to the explicitly listed native classes. Their
# non-material values still require an exact supported numeric/record type.
def _native_data_fields(cls, reference_fields=()):
    fields = {
        name for base in cls.__mro__ for name, descriptor in vars(base).items()
        if type(descriptor) in (GetSetDescriptorType, MemberDescriptorType)
        and not name.startswith('__') and name.lstrip('_') not in {field.lstrip('_') for field in reference_fields}
    }
    # A public quantity accessor and its native backing descriptor expose
    # the same value; check the backing field once.
    return tuple(name for name in sorted(fields) if '_' + name not in fields)


_KINETICS_DATA_FIELDS = {
    cls: _native_data_fields(cls, _KINETICS_REFERENCE_FIELDS[cls])
    for cls in _KINETICS_TYPES
}
# Only declared object-valued physical slots can carry an unknown type.
# Native int/bint/double slots and derived rotational constants cannot.
_NATIVE_RECORD_FIELDS = {
    TransitionState: ('label', '_frequency', 'conformer', 'tunneling'),
    Conformer: ('_E0', 'modes', '_number', '_mass', '_coordinates'),
    Mode: (), Rotation: (), Translation: (), Torsion: (), Vibration: (),
    LinearRotor: ('_inertia',), NonlinearRotor: ('_inertia',),
    KRotor: ('_inertia',), SphericalTopRotor: ('_inertia',),
    IdealGasTranslation: ('_mass',),
    HinderedRotor: ('_inertia', '_fourier', '_barrier', 'energies'),
    FreeRotor: ('_inertia',), HarmonicOscillator: ('_frequencies',),
}
# Channel-owned records share one complete slot schema in both checkers.
_NATIVE_RECORD_FIELDS.update({
    Species: ('label', 'thermo', 'conformer', 'transport_data', 'molecule',
        '_molecular_weight', 'energy_transfer_model', 'props', '_aug_inchi',
        '_state_cache_key', 'liquid_volumetric_mass_transfer_coefficient_data',
        'henry_law_constant_data', '_fingerprint', '_inchi', '_smiles'),
    Molecule: ('vertices', 'ordered_vertices', 'props', '_symm_sssr', '_sssr',
        'metal', 'facet', '_electronic_state', '_fingerprint', '_inchi', '_smiles'),
    Graph: ('vertices', 'ordered_vertices'),
    Atom: ('edges', 'mapping', 'element', 'label', 'atomtype', 'coords',
        'site', 'morphology', 'props'),
    Vertex: ('edges', 'mapping'), Bond: ('vertex1', 'vertex2'),
    Edge: ('vertex1', 'vertex2'),
    Element: ('name', 'symbol', 'chemkin_name'),
    AtomType: ('label', 'generic', 'specific', 'increment_bond', 'decrement_bond',
        'form_bond', 'break_bond', 'increment_radical', 'decrement_radical',
        'increment_lone_pair', 'decrement_lone_pair', 'increment_charge',
        'decrement_charge', 'single', 'all_double', 'r_double', 'o_double',
        's_double', 'triple', 'quadruple', 'benzene', 'lone_pairs', 'charge'),
    HeatCapacityModel: ('_Tmin', '_Tmax', '_E0', '_Cp0', '_CpInf', 'comment', 'label'),
    NASAPolynomial: ('_Tmin', '_Tmax', '_E0', '_Cp0', '_CpInf', 'comment', 'label'),
    NASA: ('_Tmin', '_Tmax', '_E0', '_Cp0', '_CpInf', 'comment', 'label',
        'poly1', 'poly2', 'poly3', '_thermo_coverage_dependence'),
    ThermoData: ('_Tmin', '_Tmax', '_E0', '_Cp0', '_CpInf', 'comment', 'label',
        '_H298', '_S298', '_Tdata', '_Cpdata', '_thermo_coverage_dependence'),
    Wilhoit: ('_Tmin', '_Tmax', '_E0', '_Cp0', '_CpInf', 'comment', 'label',
        '_B', '_H0', '_S0', '_thermo_coverage_dependence'),
    SingleExponentialDown: ('_alpha0', '_t0'),
})
_NATIVE_RECORD_READERS = {
    Species: rmgpy.reaction._native_species_record_values,
    Molecule: rmgpy.reaction._native_molecule_record_values,
    Atom: rmgpy.reaction._native_atom_record_values,
    AtomType: rmgpy.reaction._native_atomtype_record_values,
    Element: rmgpy.reaction._native_element_record_values,
}
_FAST_TYPED_RECORD_CHECKERS = {AtomType: '_native_atomtype_lists_valid'}
_PYTHON_RECORD_FIELDS = {
    TransportData: ('shapeIndex', 'epsilon', 'sigma', 'dipoleMoment',
        'polarizability', 'rotrelaxcollnum', 'comment'),
    HenryLawConstantData: ('Ts', 'kHs'),
    LiquidVolumetricMassTransferCoefficientData: ('Ts', 'kLAs'),
    SoluteData: ('S', 'B', 'E', 'L', 'A', 'V', 'comment'),
    SoluteTSData: ('Sg_g', 'Bg_g', 'Eg_g', 'Lg_g', 'Ag_g', 'Cg_g', 'Sh_g',
        'Bh_g', 'Eh_g', 'Lh_g', 'Ah_g', 'Ch_g', 'K_g', 'Sg_h', 'Bg_h',
        'Eg_h', 'Lg_h', 'Ag_h', 'Cg_h', 'Sh_h', 'Bh_h', 'Eh_h', 'Lh_h',
        'Ah_h', 'Ch_h', 'K_h', 'comment'),
}
_NUMERIC_RECORD_FIELDS = {
    ScalarQuantity: ('units', 'uncertainty_type'),
    ArrayQuantity: ('units', 'uncertainty_type', 'value_si', 'uncertainty_si'),
}
_NUMERIC_SLOT_ALIASES = {'_uncertainty_type': 'uncertainty_type'}
_SIMPLE_RECORD_FIELDS = dict(_NATIVE_RECORD_FIELDS)
_SIMPLE_RECORD_FIELDS.update(_PYTHON_RECORD_FIELDS)
_SIMPLE_RECORD_FIELDS.update({
    cls: _KINETICS_DATA_FIELDS[cls]
    for cls in (kinetics.TunnelingModel, kinetics.Wigner, kinetics.Eckart)
})
_CONTAINERS = frozenset((list, tuple, dict, set, frozenset, np.ndarray))
_NUMBERS = frozenset((str, bytes, bool, int, float, complex, type(None),
                     ScalarQuantity, ArrayQuantity)) | frozenset(
    cls for cls in np.sctypeDict.values() if isinstance(cls, type)
    and issubclass(cls, np.number)
)
def _is_numeric_data(value):
    cls = type(value)
    if type(cls) is not type or cls not in _NUMBERS:
        return False
    if cls in (ScalarQuantity, ArrayQuantity):
        if type(value.units) is not str or type(value.uncertainty_type) is not str:
            return False
    if cls is ArrayQuantity:
        for array in (value.value_si, value.uncertainty_si):
            if array is not None and (type(array) is not np.ndarray or array.dtype.hasobject):
                return False
    return True


_REACTION_MATERIAL_FIELDS = ('reactants', 'products', 'specific_collider', 'pairs')
_REACTION_RATE_FIELDS = ('kinetics', 'network_kinetics', 'SurfaceArrhenius',
                         'SurfaceChargeTransfer', 'reverse')
_REACTION_PROVENANCE_FIELDS = ('entry', 'depository', 'template', 'network',
                              'library', 'family')
_SOLUTE_DATA_FIELDS = frozenset(('S', 'B', 'E', 'L', 'A', 'V', 'comment',
    'Sg_g', 'Bg_g', 'Eg_g', 'Lg_g', 'Ag_g', 'Cg_g', 'Sh_g', 'Bh_g', 'Eh_g',
    'Lh_g', 'Ah_g', 'Ch_g', 'K_g', 'Sg_h', 'Bg_h', 'Eg_h', 'Lg_h', 'Ag_h',
    'Cg_h', 'Sh_h', 'Bh_h', 'Eh_h', 'Lh_h', 'Ah_h', 'Ch_h', 'K_h'))
_REACTION_DATA_FIELDS = frozenset((
    'index', 'label', 'reversible', 'transition_state', 'duplicate', 'degeneracy',
    '_degeneracy', 'electrons', 'protons', '_protons', 'allow_pdep_route',
    'elementary_high_p', 'comment', 'k_effective_cache', 'is_forward',
    'allow_max_rate_violation', 'rank', 'estimator',
))

_REACTION_FIELDS = (_REACTION_DATA_FIELDS | frozenset(_REACTION_MATERIAL_FIELDS)
                    | frozenset(_REACTION_RATE_FIELDS) | frozenset(_REACTION_PROVENANCE_FIELDS)
                    | frozenset(('labeled_atoms',)))


@lru_cache(maxsize=1)
def _reaction_types():
    # Keep family/depository imports out of this module's initialization cycle.
    from rmgpy.data.kinetics.family import TemplateReaction
    from rmgpy.data.kinetics.depository import DepositoryReaction
    return (rmgpy.reaction.Reaction, LibraryReaction, TemplateReaction,
            DepositoryReaction, PDepReaction)


@lru_cache(maxsize=1)
def _provenance_types():
    from rmgpy.data.kinetics.depository import KineticsDepository
    from rmgpy.data.kinetics.library import KineticsLibrary
    from rmgpy.data.kinetics.family import KineticsFamily
    from rmgpy.data.base import LogicOr, LogicAnd
    return (Database, Entry, KineticsLibrary, KineticsDepository, KineticsFamily,
            Group, LogicOr, LogicAnd, rmgpy.pdep.network.Network, PDepNetwork)


def _check_template_atom_labels(reaction, labels, position):
    """Labels may only alias native atoms in this reaction's checked molecules."""
    def refuse():
        raise _UnclassifiedSpeciesReference('Unclassified atom labels at ' + position)

    if type(labels) is not dict:
        refuse()
    for side, mapping in labels.items():
        if type(side) is not str or side not in ('reactants', 'products') or type(mapping) is not dict:
            refuse()
        if not mapping:
            continue
        participants = getattr(reaction, side)
        if type(participants) is not list:
            refuse()
        owned = set()
        for participant in participants:
            if type(participant) is Species and type(participant.molecule) is list:
                molecules = participant.molecule
            elif type(participant) is Molecule:
                molecules = [participant]
            else:
                refuse()
            for molecule in molecules:
                if type(molecule) is not Molecule or type(molecule.atoms) is not list:
                    refuse()
                if any(type(atom) is not Atom for atom in molecule.atoms):
                    refuse()
                owned.update(id(atom) for atom in molecule.atoms)
        for label, reference in mapping.items():
            if type(label) is not str:
                refuse()
            atoms = reference if type(reference) is list else [reference]
            if any(type(atom) is not Atom or id(atom) not in owned for atom in atoms):
                refuse()


def _reaction_species_references(reaction, seen=None, atomtype_seen=None):
    """Traverse current declared schemas in the compiled general walker."""
    yield from rmgpy.reaction._native_reaction_species_references(
        reaction, seen, atomtype_seen)


def _has_electron_participant(reaction, seen=None, atomtype_seen=None):
    """Recognize known material/metadata electrons; unknown schemas refuse."""
    if type(type(reaction)) is not type or type(reaction) not in _reaction_types():
        return True
    if reaction.electrons:
        return True
    if type(reaction) is rmgpy.reaction.Reaction:
        simple = rmgpy.reaction._simple_native_channel_verdict(
            reaction, _KINETICS_REFERENCE_FIELDS[kinetics.Arrhenius],
            _KINETICS_DATA_FIELDS[kinetics.Arrhenius], _NUMBERS, _SIMPLE_RECORD_FIELDS,
            _NATIVE_RECORD_READERS, _PYTHON_RECORD_FIELDS, seen, atomtype_seen)
        if simple >= 0:
            return bool(simple)
    try:
        if any(type(reference) is Molecule and reference.is_electron()
               for reference in _reaction_species_references(reaction, seen, atomtype_seen)):
            return True
    except _UnclassifiedSpeciesReference:
        return True
    return (type(reaction) is not rmgpy.reaction.Reaction
            and get_placement_owner(reaction) is not None)


def _check_electron_channel_routing(reaction, library=None, seen=None, atomtype_seen=None):
    """Refuse every electron reaction before network state changes."""
    if _has_electron_participant(reaction, seen, atomtype_seen):
        if type(type(reaction)) is not type or type(reaction) not in _reaction_types():
            source = library if type(library) is str and library else 'Unclassified'
            raise NetworkError('Electron reactions cannot enter pressure-dependent networks: '
                               '<unclassified reaction> (library: {0}). Unsupported exact reaction type.'.format(source))
        census_error = ''
        try:
            for reference in _reaction_species_references(reaction):
                pass
        except _UnclassifiedSpeciesReference as error:
            census_error = str(error) + '. '
        candidates = (library, getattr(reaction, 'library', None),
                      getattr(reaction, 'family', None))
        names = []
        for candidate in candidates:
            if type(candidate) is str:
                names.append(candidate)
            elif type(type(candidate)) is type and type(candidate) in _provenance_types()[:5]:
                names.append(candidate.label)
        source = next((name for name in names if type(name) is str and name), 'Unclassified')
        # The diagnostic must not execute a refused participant's properties
        # or __str__, and must not initialize any labels on the candidate.
        def describe(participant):
            if type(participant) is Species:
                label = participant.label
                if type(label) is not str or not label:
                    molecules = participant.molecule
                    label = (molecules[0].to_smiles()
                             if type(molecules) is list and molecules
                             and type(molecules[0]) is Molecule else '<species>')
                if type(participant.index) is int and participant.index != -1:
                    return '{0}({1:d})'.format(label, participant.index)
                return label
            if type(participant) is Molecule:
                return str(participant)
            return '<unsupported participant type>'

        def side(participants):
            if type(participants) not in (list, tuple):
                return '<unsupported participant-list type>'
            return ' + '.join(describe(participant) for participant in participants)

        reactants = side(reaction.reactants)
        products = side(reaction.products)
        collider = (' (+{0})'.format(describe(reaction.specific_collider))
                    if reaction.specific_collider is not None else '')
        arrow = ' <=> ' if reaction.reversible else ' => '
        equation = reactants + collider + arrow + products + collider
        raise NetworkError(
            'Electron reactions cannot enter pressure-dependent networks: '
            '{0} (library: {1}). {2}Keep the rate explicit; disable '
            'elementary_high_p/allow_pdep_route and cached network kinetics, '
            'or disable pressure dependence.'.format(equation, source, census_error))


class PDepReaction(rmgpy.reaction.Reaction):

    def __init__(self,
                 index=-1,
                 label='',
                 reactants=None,
                 products=None,
                 specific_collider=None,
                 network=None,
                 kinetics=None,
                 network_kinetics=None,
                 reversible=True,
                 transition_state=None,
                 duplicate=False,
                 degeneracy=1,
                 pairs=None,
                 electrons=0
                 ):
        rmgpy.reaction.Reaction.__init__(self,
                                         index=index,
                                         label=label,
                                         reactants=reactants,
                                         products=products,
                                         specific_collider=specific_collider,
                                         kinetics=kinetics,
                                         network_kinetics=network_kinetics,
                                         reversible=reversible,
                                         transition_state=transition_state,
                                         duplicate=duplicate,
                                         degeneracy=degeneracy,
                                         pairs=pairs,
                                         electrons=electrons
                                         )
        self.network = network

    def __reduce__(self):
        """
        A helper function used when pickling an object.

        No constructor arguments: every field travels in the state dict, discovered from
        the object. The fourteen-item list this replaces dropped `comment`, `rank`,
        `is_forward`, `elementary_high_p`, `allow_pdep_route` and
        `allow_max_rate_violation` -- measured at `13e3227b2`, a round trip turned
        ``allow_max_rate_violation=True`` into ``False`` and ``rank=3`` into ``None``. A
        pressure-dependent reaction is the one shape for which `allow_pdep_route` and
        `elementary_high_p` are most likely to be set, and this is the reducer that lost
        them.

        Imported inside the method rather than at the top of the module, to keep
        `rmgpy.rmg` out of the import path of `rmgpy.data.kinetics`; `depository.py`
        states the mirror image of this note.
        """
        from rmgpy.data.kinetics.family import _NOT_REPRODUCED, reaction_state, state_fields
        return (PDepReaction, (), reaction_state(self, _NOT_REPRODUCED,
                                                 state_fields(self)))

    def __setstate__(self, state):
        """
        Restore what `__reduce__` sent; see `TemplateReaction.__setstate__` for why the
        default unpickling path will not do.
        """
        from rmgpy.data.kinetics.family import apply_reaction_state
        apply_reaction_state(self, state)

    def copy(self):
        """
        Create a deep copy of this reaction, as a `PDepReaction`.

        Inherited, `Reaction.copy` builds a *base* `Reaction`, so the copy lost `network`
        -- the only thing that says which pressure-dependent network this reaction belongs
        to, and what `get_source` below reports. `rmgpy/tools/isotopes.py:394` copies
        whatever `Reaction` it is handed, and a core model that ran pressure dependence is
        full of these.

        The `network` back-pointer is carried by reference, deliberately: it holds this
        reaction in turn, so deepening it would reproduce the whole network around one of
        its own members. That is written down in `_COPIED_BY_REFERENCE` beside the reason,
        which is what keeps it a decision rather than an omission.
        """
        from rmgpy.data.kinetics.family import _NOT_COPIED_BY_REFERENCE, copy_reaction
        return copy_reaction(self, _NOT_COPIED_BY_REFERENCE)

    def __deepcopy__(self, memo):
        """
        The same override, for the same reason; see `TemplateReaction.__deepcopy__`.
        """
        memo[id(self)] = other = self.copy()
        return other

    def get_source(self):
        """
        Get the source of this PDepReaction
        """
        return str(self.network)


################################################################################

def _check_network_reactions(network):
    """Validate each current channel and shared physical record once per call.

    The visited sets are local to this operation. Admission never lends a
    verdict to a later computation, so subsequent mutations are inspected.
    Plain Network is supported for direct numerical and Arkane callers.
    """
    if type(network) not in (PDepNetwork, rmgpy.pdep.network.Network):
        raise NetworkError('Electron reactions cannot enter pressure-dependent networks: '
                           '<unsupported network type> (library: Unclassified).')
    seen, atomtype_seen, channels = set(), set(), set()
    for reactions in (network.path_reactions, network.net_reactions):
        if type(reactions) not in (list, tuple):
            raise NetworkError('Electron reactions cannot enter pressure-dependent networks: '
                               '<unsupported reaction-list type> (library: Unclassified).')
        for reaction in reactions:
            if id(reaction) not in channels:
                channels.add(id(reaction))
                _check_electron_channel_routing(reaction, seen=seen, atomtype_seen=atomtype_seen)


def _checked_network_entry(function):
    """Check current channels once at an explicit computation/routing entry."""
    @wraps(function)
    def checked(self, *args, **kwargs):
        if function.__name__ in ('add_path_reaction', 'add_net_reaction'):
            if type(self) is not PDepNetwork:
                raise NetworkError('Electron reactions cannot enter pressure-dependent networks: '
                                   '<unsupported network type> (library: Unclassified).')
        else:
            PDepNetwork._check_reactions(self)
        return function(self, *args, **kwargs)
    checked._electron_channel_entry = True
    return checked


_NETWORK_ACCESSOR_EXCLUSIONS = {
    'invalidate': 'Only clears the validity flag; it neither reads nor routes any channel.',
    'get_all_species': 'Enumerates existing configurations and bath gas without evaluating or registering reactions.',
    'log_summary': 'Formats diagnostic state for logging; no rate, grid, routing or network computation is performed.',
    'cleanup': 'Releases existing numerical arrays and configuration caches; it does not compute or route channels.',
}


def _guard_network_entries(cls):
    """Cover the whole public network API; only named pure accessors are exempt."""
    import inspect
    for name, function in inspect.getmembers(cls, inspect.isroutine):
        if (not name.startswith('_') and name not in _NETWORK_ACCESSOR_EXCLUSIONS
                and not getattr(function, '_electron_channel_entry', False)):
            setattr(cls, name, _checked_network_entry(function))
    return cls


@_guard_network_entries
class PDepNetwork(rmgpy.pdep.network.Network):
    """
    A representation of a *partial* unimolecular reaction network. Each partial
    network has a single `source` isomer or reactant channel, and is responsible
    only for :math:`k(T,P)` values for net reactions with source as the
    reactant. Multiple partial networks can have the same source, but networks
    with the same source and any explored isomers must be combined.

    =================== ======================= ================================
    Attribute           Type                    Description
    =================== ======================= ================================
    `source`            ``list``                The isomer or reactant channel that acts as the source
    `explored`          ``list``                A list of the unimolecular isomers whose reactions have been fully explored
    =================== ======================= ================================

    """

    def __init__(self, index=-1, source=None):
        rmgpy.pdep.network.Network.__init__(self, label="PDepNetwork #{0}".format(index))
        self.index = index
        self.source = source
        self.energy_correction = None
        self.explored = []
        self.products_cache = []

    def __str__(self):
        return "PDepNetwork #{0}".format(self.index)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (PDepNetwork, (self.index, self.source), self.__dict__)

    def __setstate__(self, state):
        if type(state) is not dict:
            raise NetworkError('Electron reactions cannot enter pressure-dependent networks: '
                               '<unsupported network state> (library: Unclassified).')
        for name in ('path_reactions', 'net_reactions'):
            reactions = state.get(name, [])
            if type(reactions) not in (list, tuple):
                raise NetworkError('Electron reactions cannot enter pressure-dependent networks: '
                                   '<unsupported reaction-list type> (library: Unclassified).')
            for reaction in reactions:
                _check_electron_channel_routing(reaction)
        self.__dict__.update(state)

    def _check_reactions(self):
        """Validate current path/net channels at routing/computation boundaries.

        A reaction in a cyclic pickle may still be a placeholder in __setstate__.
        Entry checks run on the complete graph, and never cache a verdict.
        Attribute access, printing and inspection do not compute or route it.
        """
        if type(self) is not PDepNetwork:
            raise NetworkError('Electron reactions cannot enter pressure-dependent networks: '
                               '<unsupported network type> (library: Unclassified).')
        _check_network_reactions(self)

    def add_net_reaction(self, reaction):
        """Validate only the incoming channel before registration."""
        _check_electron_channel_routing(reaction)
        self.net_reactions.append(reaction)

    def cleanup(self):
        """
        Delete intermedate arrays used to compute k(T,P) values.
        """
        for isomer in self.isomers:
            isomer.cleanup()
        for reactant in self.reactants:
            reactant.cleanup()
        for product in self.products:
            product.cleanup()

        self.e_list = None
        self.j_list = None
        self.dens_states = None
        self.coll_freq = None
        self.Mcoll = None
        self.Kij = None
        self.Fim = None
        self.Gnj = None
        self.E0 = None
        self.n_grains = 0
        self.n_j = 0

        self.K = None
        self.p0 = None

    @_checked_network_entry
    def get_leak_coefficient(self, T, P):
        """
        Return the pressure-dependent rate coefficient :math:`k(T,P)` describing
        the total rate of "leak" from this network. This is defined as the sum
        of the :math:`k(T,P)` values for all net reactions to nonexplored
        unimolecular isomers.
        """
        k = 0.0
        if len(self.net_reactions) == 0 and len(self.path_reactions) == 1:
            # The network is of the form A + B -> C* (with C* nonincluded)
            # For this special case we use the high-pressure limit k(T) to
            # ensure that we're estimating the total leak flux
            rxn = self.path_reactions[0]
            if rxn.kinetics is None:
                if rxn.reverse.kinetics is not None:
                    rxn = rxn.reverse
                else:
                    raise PressureDependenceError('Path reaction {0} with no high-pressure-limit kinetics encountered '
                                                  'in PDepNetwork #{1:d} while evaluating leak flux.'.format(rxn, self.index))
            if rxn.products is self.source:
                rxn.check_resolved_species_reversibility(reversible=True)
                k = rxn.get_rate_coefficient(T, P) / rxn.get_equilibrium_constant(T)
            else:
                k = rxn.get_rate_coefficient(T, P)
        else:
            # The network has at least one included isomer, so we can calculate
            # the leak flux normally
            for rxn in self.net_reactions:
                if len(rxn.products) == 1 and rxn.products[0] not in self.explored:
                    k += rxn.get_rate_coefficient(T, P)
        return k

    @_checked_network_entry
    def get_maximum_leak_species(self, T, P):
        """
        Get the unexplored (unimolecular) isomer with the maximum leak flux.
        Note that the leak rate coefficients vary with temperature and
        pressure, so you must provide these in order to get a meaningful result.
        """
        # Choose species with maximum leak flux
        max_k = 0.0
        max_species = None
        if len(self.net_reactions) == 0 and len(self.path_reactions) == 1:
            max_k = self.get_leak_coefficient(T, P)
            rxn = self.path_reactions[0]
            if rxn.products == self.source:
                assert len(rxn.reactants) == 1
                max_species = rxn.reactants[0]
            else:
                assert len(rxn.products) == 1
                max_species = rxn.products[0]
        else:
            for rxn in self.net_reactions:
                if len(rxn.products) == 1 and rxn.products[0] not in self.explored:
                    k = rxn.get_rate_coefficient(T, P)
                    if max_species is None or k > max_k:
                        max_species = rxn.products[0]
                        max_k = k

        # Make sure we've identified a species
        if max_species is None:
            raise NetworkError('No unimolecular isomers left to explore!')
        # Return the species
        return max_species

    @_checked_network_entry
    def get_leak_branching_ratios(self, T, P):
        """
        Return a dict with the unexplored isomers in the partial network as the
        keys and the fraction of the total leak coefficient as the values.
        """
        ratios = {}
        if len(self.net_reactions) == 0 and len(self.path_reactions) == 1:
            rxn = self.path_reactions[0]
            assert rxn.reactants == self.source or rxn.products == self.source
            if rxn.products == self.source:
                assert len(rxn.reactants) == 1
                ratios[rxn.reactants[0]] = 1.0
            else:
                assert len(rxn.products) == 1
                ratios[rxn.products[0]] = 1.0
        else:
            for rxn in self.net_reactions:
                if len(rxn.products) == 1 and rxn.products[0] not in self.explored:
                    ratios[rxn.products[0]] = rxn.get_rate_coefficient(T, P)

        kleak = sum(ratios.values())
        for spec in ratios:
            ratios[spec] /= kleak

        return ratios

    @_checked_network_entry
    def explore_isomer(self, isomer):
        """
        Explore a previously-unexplored unimolecular `isomer` in this partial
        network using the provided core-edge reaction model `reaction_model`,
        returning the new reactions and new species.
        """

        if isomer in self.explored:
            logging.warning('Already explored isomer {0} in pressure-dependent network #{1:d}'.format(isomer,
                                                                                                      self.index))
            return []

        assert isomer not in self.source, "Attempted to explore isomer {0}, but that is the source configuration for this network.".format(isomer)

        for product in self.products:
            if product.species == [isomer]:
                break
        else:
            raise Exception('Attempted to explore isomer {0}, but that species not found in product channels.'.format(isomer))

        logging.info('Exploring isomer {0} in pressure-dependent network #{1:d}'.format(isomer, self.index))

        for mol in isomer.molecule:
            mol.update()

        # Find reactions involving the found species as unimolecular
        # reactant or product (e.g. A <---> products)

        # Don't find reactions involving the new species as bimolecular
        # reactants or products with itself (e.g. A + A <---> products)
        # Don't find reactions involving the new species as bimolecular
        # reactants or products with other core species (e.g. A + B <---> products)

        new_reactions = react_species((isomer,))
        for reaction in new_reactions:
            _check_electron_channel_routing(reaction)
        self.explored.append(isomer)
        self.isomers.append(product)
        self.products.remove(product)

        return new_reactions

    def add_path_reaction(self, newReaction):
        """
        Add a path reaction to the network. If the path reaction already exists,
        no action is taken.
        """
        _check_electron_channel_routing(newReaction)
        newReaction.check_resolved_species_reversibility()
        if newReaction.network_kinetics is not None:
            newReaction.check_resolved_species_reversibility(kinetics=newReaction.network_kinetics)
        # Add this reaction to that network if not already present
        found = False
        for rxn in self.path_reactions:
            if newReaction.reactants == rxn.reactants and newReaction.products == rxn.products:
                found = True
                break
            elif (newReaction.allows_reverse_match(rxn)
                  and newReaction.products == rxn.reactants and newReaction.reactants == rxn.products):
                found = True
                break
        if not found:
            self.path_reactions.append(newReaction)
            self.invalidate()

    @_checked_network_entry
    def get_energy_filtered_reactions(self, T, tol):
        """
        Returns a list of products and isomers that are greater in Free Energy
        than a*R*T + Gfsource(T)
        """
        dE = tol * R * T
        for conf in self.isomers + self.products + self.reactants:
            if len(conf.species) == len(self.source):
                if len(self.source) == 1:
                    if self.source[0].is_isomorphic(conf.species[0]):
                        E0source = conf.E0
                        break
                elif len(self.source) == 2:
                    boo00 = self.source[0].is_isomorphic(conf.species[0])
                    boo01 = self.source[0].is_isomorphic(conf.species[1])
                    if boo00 or boo01:  # if we found source[0]
                        boo10 = self.source[1].is_isomorphic(conf.species[0])
                        boo11 = self.source[1].is_isomorphic(conf.species[1])
                        if (boo00 and boo11) or (boo01 and boo10):
                            E0source = conf.E0
                            break
        else:
            raise ValueError('No isomer, product or reactant channel is isomorphic to the source')

        filtered_rxns = []
        for rxn in self.path_reactions:
            E0 = rxn.transition_state.conformer.E0.value_si
            if E0 - E0source > dE:
                filtered_rxns.append(rxn)

        return filtered_rxns

    @_checked_network_entry
    def get_rate_filtered_products(self, T, P, tol):
        """
        determines the set of path_reactions that have fluxes less than
        tol at steady state where all A => B + C reactions are irreversible
        and there is a constant flux from/to the source configuration of 1.0
        """
        c = self.solve_ss_network(T, P)
        isomer_spcs = [iso.species[0] for iso in self.isomers]
        filtered_prod = []
        if c is not None:
            for rxn in self.net_reactions:
                val = 0.0
                val2 = 0.0
                if rxn.reactants[0] in isomer_spcs:
                    ind = isomer_spcs.index(rxn.reactants[0])
                    kf = rxn.get_rate_coefficient(T, P)
                    val = kf * c[ind]
                if rxn.products[0] in isomer_spcs:
                    ind2 = isomer_spcs.index(rxn.products[0])
                    rxn.check_resolved_species_reversibility(reversible=True)
                    kr = rxn.get_rate_coefficient(T, P) / rxn.get_equilibrium_constant(T)
                    val2 = kr * c[ind2]

                if max(val, val2) < tol:
                    filtered_prod.append(rxn.products)

            return filtered_prod

        else:
            logging.warning("Falling back flux reduction from Steady State analysis to rate coefficient analysis")
            ks = np.array([rxn.get_rate_coefficient(T, P) for rxn in self.net_reactions])
            frs = ks / ks.sum()
            inds = [i for i in range(len(frs)) if frs[i] < tol]
            filtered_prod = [self.net_reactions[i].products for i in inds]
            return filtered_prod

    @_checked_network_entry
    def solve_ss_network(self, T, P):
        """
        calculates the steady state concentrations if all A => B + C
        reactions are irreversible and the flux from/to the source
        configuration is 1.0
        """
        A = np.zeros((len(self.isomers), len(self.isomers)))
        b = np.zeros(len(self.isomers))
        bimolecular = len(self.source) > 1

        isomer_spcs = [iso.species[0] for iso in self.isomers]

        for rxn in self.net_reactions:
            if rxn.reactants[0] in isomer_spcs:
                ind = isomer_spcs.index(rxn.reactants[0])
                kf = rxn.get_rate_coefficient(T, P)
                A[ind, ind] -= kf
            else:
                ind = None
            if rxn.products[0] in isomer_spcs:
                ind2 = isomer_spcs.index(rxn.products[0])
                rxn.check_resolved_species_reversibility(reversible=True)
                kr = rxn.get_rate_coefficient(T, P) / rxn.get_equilibrium_constant(T)
                A[ind2, ind2] -= kr
            else:
                ind2 = None

            if ind is not None and ind2 is not None:
                A[ind, ind2] += kr
                A[ind2, ind] += kf

            if bimolecular:
                if rxn.reactants[0] == self.source:
                    kf = rxn.get_rate_coefficient(T, P)
                    b[ind2] += kf
                elif rxn.products[0] == self.source:
                    rxn.check_resolved_species_reversibility(reversible=True)
                    kr = rxn.get_rate_coefficient(T, P) / rxn.get_equilibrium_constant(T)
                    b[ind] += kr

        if not bimolecular:
            ind = isomer_spcs.index(self.source[0])
            b[ind] = -1.0  # flux at source
        else:
            total_source_flux = b.sum()
            if total_source_flux == 0:
                return None
            b = -b / total_source_flux  # 1.0 flux from source

        if len(b) == 1:
            if A[0, 0] == 0:
                return None
            return np.array([b[0] / A[0, 0]])

        con = np.linalg.cond(A)

        if np.log10(con) < 15:
            c = np.linalg.solve(A, b)
        else:
            logging.warning("Matrix Ill-conditioned, attempting to use Arbitrary Precision Arithmetic")
            mp.dps = 30 + int(np.log10(con))
            Amp = mp.matrix(A.tolist())
            bmp = mp.matrix(b.tolist())

            try:
                c = mp.qr_solve(Amp, bmp)

                c = np.array(list(c[0]))

                if any(c <= 0.0):
                    c, rnorm = opt.nnls(A, b)

                c = c.astype(float)
            except:  # fall back to raw flux analysis rather than solve steady state problem
                return None
        
        if np.isnan(c).any():
            return None
        
        return c

    def remove_disconnected_reactions(self):
        """
        gets rid of reactions/isomers/products not connected to the source by a reaction sequence
        """
        kept_reactions = []
        kept_products = [self.source]
        incomplete = True
        while incomplete:
            s = len(kept_reactions)
            for rxn in self.path_reactions:
                if not rxn in kept_reactions:
                    if rxn.reactants in kept_products:
                        kept_products.append(rxn.products)
                        kept_reactions.append(rxn)
                    elif rxn.products in kept_products:
                        kept_products.append(rxn.reactants)
                        kept_reactions.append(rxn)

            incomplete = s != len(kept_reactions)

        logging.info('Removing disconnected items')
        for rxn in self.path_reactions:
            if rxn not in kept_reactions:
                logging.info('Removing rxn: {}'.format(rxn))
                self.path_reactions.remove(rxn)

        nrxns = []
        for nrxn in self.net_reactions:
            if nrxn.products not in kept_products or nrxn.reactants not in kept_products:
                logging.info('Removing net rxn: {}'.format(nrxn))
            else:
                logging.info('Keeping net rxn: {}'.format(nrxn))
                nrxns.append(nrxn)
        self.net_reactions = nrxns

        prods = []
        for prod in self.products:
            if prod.species not in kept_products:
                logging.info('Removing product: {}'.format(prod))
            else:
                logging.info("Keeping product: {}".format(prod))
                prods.append(prod)

        self.products = prods

        rcts = []
        for rct in self.reactants:
            if rct.species not in kept_products:
                logging.info('Removing product: {}'.format(rct))
            else:
                logging.info("Keeping product: {}".format(rct))
                rcts.append(rct)
        self.reactants = rcts

        isos = []
        for iso in self.isomers:
            if iso.species not in kept_products:
                logging.info('Removing isomer: {}'.format(iso))
            else:
                logging.info("Keeping isomer: {}".format(iso))
                isos.append(iso)

        self.isomers = isos
        self.explored = [iso.species[0] for iso in isos]

        self.n_isom = len(self.isomers)
        self.n_reac = len(self.reactants)
        self.n_prod = len(self.products)

    def remove_reactions(self, reaction_model, networks, rxns=None, prods=None):
        """
        removes a list of reactions from the network and all reactions/products
        left disconnected by removing those reactions
        """
        if rxns:
            for rxn in rxns:
                self.path_reactions.remove(rxn)

        if prods:
            isomers = [x.species[0] for x in self.isomers]

            for prod in prods:
                prod = [x for x in prod]
                if prod[0] in isomers:  # skip isomers
                    continue
                for rxn in self.path_reactions:
                    if rxn.products == prod or rxn.reactants == prod:
                        self.path_reactions.remove(rxn)

            prodspc = [x[0] for x in prods]
            for prod in prods:
                prod = [x for x in prod]
                if prod[0] in isomers:  # deal with isomers
                    for rxn in self.path_reactions:
                        if rxn.reactants == prod and rxn.products[0] not in isomers and rxn.products[0] not in prodspc:
                            break
                        if rxn.products == prod and rxn.reactants[0] not in isomers and rxn.reactants not in prodspc:
                            break
                    else:
                        for rxn in self.path_reactions:
                            if rxn.reactants == prod or rxn.products == prod:
                                self.path_reactions.remove(rxn)

        self.remove_disconnected_reactions()

        self.cleanup()

        self.invalidate()

        assert self.path_reactions != [], 'Reduction process removed all reactions, cannot update network with no reactions'

        reaction_model.update_unimolecular_reaction_networks()

        if reaction_model.pressure_dependence.output_file:
            path0 = os.path.join(reaction_model.pressure_dependence.output_file, 'pdep')
            path = os.path.join(reaction_model.pressure_dependence.output_file, 'pdep','final')
            if not os.path.exists(path):
                os.mkdir(path)
            for name in os.listdir(path0):
                if name.endswith('.py') and '_' in name:
                    s1,s2 = name.split('_')
                    index = int(s1[7:])
                    N_isomers = int(s2.split('.')[0]) 
                    if index == self.index and N_isomers == len(self.isomers):
                        shutil.copy(os.path.join(path0, name),
                                  os.path.join(path, 'network{}_reduced.py'.format(networks.index(self))))


    def merge(self, other):
        """
        Merge the partial network `other` into this network.
        """
        PDepNetwork._check_reactions(other)
        # Make sure the two partial networks have the same source configuration
        assert self.source == other.source

        # Merge isomers
        for isomer in other.isomers:
            if isomer not in self.isomers:
                self.isomers.append(isomer)
        # Merge explored
        for isomer in other.explored:
            if isomer not in self.explored:
                self.explored.append(isomer)
        # Merge reactants
        for reactants in other.reactants:
            if reactants not in self.reactants:
                self.reactants.append(reactants)
        # Merge products
        for products in other.products:
            if products not in self.products:
                self.products.append(products)

        # However, products that have been explored are actually isomers
        # These should be removed from the list of products!
        products_to_remove = []
        for products in self.products:
            if len(products.species) == 1 and products.species[0] in self.isomers:
                products_to_remove.append(products)
        for products in products_to_remove:
            self.products.remove(products)

        # Merge path reactions
        for reaction in other.path_reactions:
            found = False
            for rxn in self.path_reactions:
                if reaction.reactants == rxn.reactants and reaction.products == rxn.products:
                    # NB the isEquivalent() method that used to be on the previous line also checked reverse direction.
                    # I am not sure which is appropriate 
                    found = True
                    break
            if not found:
                self.path_reactions.append(reaction)

        # Also merge net reactions (so that when we update the network in the
        # future, we update the existing net reactions rather than making new ones)
        # Q: What to do when a net reaction exists in both networks being merged?
        for reaction in other.net_reactions:
            found = False
            for rxn in self.net_reactions:
                if reaction.reactants == rxn.reactants and reaction.products == rxn.products:
                    # NB the isEquivalent() method that used to be on the previous line also checked reverse direction.
                    # I am not sure which is appropriate 
                    found = True
                    break
            if not found:
                self.add_net_reaction(reaction)

        # Mark this network as invalid
        self.valid = False

    @_checked_network_entry
    def update_configurations(self, reaction_model):
        """
        Sort the reactants and products of each of the network's path reactions
        into isomers, reactant channels, and product channels. You must pass 
        the current `reaction_model` because some decisions on sorting are made
        based on which species are in the model core. 
        """
        reactants = []
        products = []

        # All explored species are isomers
        isomers = self.explored[:]

        # The source configuration is an isomer (if unimolecular) or a reactant channel (if bimolecular)
        if len(self.source) == 1:
            # The source is a unimolecular isomer
            if self.source[0] not in isomers: isomers.insert(0, self.source[0])
        else:
            # The source is a bimolecular reactant channel
            self.source.sort()
            reactants.append(self.source)

        # Iterate over path reactions and make sure each set of reactants and products is classified
        for rxn in self.path_reactions:
            # Sort bimolecular configurations so that we always encounter them in the
            # same order
            # The actual order doesn't matter, as long as it is consistent
            rxn.reactants.sort()
            rxn.products.sort()
            # Reactants of the path reaction
            if len(rxn.reactants) == 1 and rxn.reactants[0] not in isomers and rxn.reactants not in products:
                # We've encountered a unimolecular reactant that is not classified
                # These are always product channels (since they would be in source or explored otherwise)
                products.append(rxn.reactants)
            elif len(rxn.reactants) > 1 and rxn.reactants not in reactants and rxn.reactants not in products:
                # We've encountered bimolecular reactants that are not classified
                if all([reactant in reaction_model.core.species for reactant in rxn.reactants]):
                    # Both reactants are in the core, so treat as reactant channel
                    reactants.append(rxn.reactants)
                else:
                    # One or more reactants is an edge species, so treat as product channel
                    products.append(rxn.reactants)
            # Products of the path reaction
            if len(rxn.products) == 1 and rxn.products[0] not in isomers and rxn.products not in products:
                # We've encountered a unimolecular product that is not classified
                # These are always product channels (since they would be in source or explored otherwise)
                products.append(rxn.products)
            elif len(rxn.products) > 1 and rxn.products not in reactants and rxn.products not in products:
                # We've encountered bimolecular products that are not classified
                if all([product in reaction_model.core.species for product in rxn.products]):
                    # Both products are in the core, so treat as reactant channel
                    reactants.append(rxn.products)
                else:
                    # One or more reactants is an edge species, so treat as product channel
                    products.append(rxn.products)

        # Clear existing configurations
        self.isomers = []
        self.reactants = []
        self.products = []

        # Make a configuration object for each
        for isomer in isomers:
            self.isomers.append(Configuration(isomer))
        for reactant in reactants:
            self.reactants.append(Configuration(*reactant))
        for product in products:
            self.products.append(Configuration(*product))
        if self.energy_correction:
            for spec in self.reactants + self.products + self.isomers:
                spec.energy_correction = self.energy_correction

    def add_products_to_reactants(self):
        self.products_cache = self.products
        self.products = []
        self.reactants = self.reactants + self.products_cache

    def remove_products_from_reactants(self):
        if self.products_cache != []:
            for prod in self.products_cache:
                self.reactants.remove(prod)
            self.products = self.products_cache

    def update(self, reaction_model, pdep_settings, requires_rms=False):
        """
        Regenerate the :math:`k(T,P)` values for this partial network if the
        network is marked as invalid.
        """
        from rmgpy.thermo.state import require_network_thermo_allowed
        require_network_thermo_allowed(self)
        from rmgpy.kinetics import Arrhenius, KineticsData, MultiArrhenius

        # Path rates can have been assigned or their participants resolved after
        # admission. A pressure-dependent fit must not erase Te dependence.
        for reaction in self.path_reactions:
            reaction.check_resolved_species_reversibility(reversible=True)
            if reaction.network_kinetics is not None:
                reaction.check_resolved_species_reversibility(
                    kinetics=reaction.network_kinetics, reversible=True)

        # Get the parameters for the pressure dependence calculation
        job = pdep_settings
        job.network = self
        output_directory = pdep_settings.output_file

        Tmin = job.Tmin.value_si
        Tmax = job.Tmax.value_si
        Pmin = job.Pmin.value_si
        Pmax = job.Pmax.value_si
        Tlist = job.Tlist.value_si
        Plist = job.Plist.value_si
        maximum_grain_size = job.maximum_grain_size.value_si if job.maximum_grain_size is not None else 0.0
        minimum_grain_count = job.minimum_grain_count
        method = job.method
        interpolation_model = job.interpolation_model
        active_j_rotor = job.active_j_rotor
        active_k_rotor = job.active_k_rotor
        rmgmode = job.rmgmode

        # Figure out which configurations are isomers, reactant channels, and product channels
        self.update_configurations(reaction_model)

        if "simulation least squares" in method or method == "chemically-significant eigenvalues georgievskii":
            self.add_products_to_reactants()

        # Make sure we have high-P kinetics for all path reactions
        for rxn in self.path_reactions:
            if rxn.kinetics is None and rxn.reverse.kinetics is None:
                raise PressureDependenceError('Path reaction {0} with no high-pressure-limit kinetics encountered in '
                                              'PDepNetwork #{1:d}.'.format(rxn, self.index))
            elif rxn.kinetics is not None and rxn.kinetics.is_pressure_dependent() and rxn.network_kinetics is None:
                raise PressureDependenceError('Pressure-dependent kinetics encountered for path reaction {0} in '
                                              'PDepNetwork #{1:d}.'.format(rxn, self.index))

        # Do nothing if the network is already valid
        if self.valid:
            self.remove_products_from_reactants()
            return
        # Do nothing if there are no explored wells
        if len(self.explored) == 0 and len(self.source) > 1:
            self.remove_products_from_reactants()
            return
        # Log the network being updated
        logging.info("Updating {0!s}".format(self))

        # Generate states data for unimolecular isomers and reactants if necessary
        for isomer in self.isomers:
            spec = isomer.species[0]
            if not spec.has_statmech():
                spec.generate_statmech()
        for reactants in self.reactants:
            for spec in reactants.species:
                if not spec.has_statmech():
                    spec.generate_statmech()
        # Also generate states data for any path reaction reactants, so we can
        # always apply the ILT method in the direction the kinetics are known
        for reaction in self.path_reactions:
            for spec in reaction.reactants:
                if not spec.has_statmech():
                    spec.generate_statmech()
        # While we don't need the frequencies for product channels, we do need
        # the E0, so create a conformer object with the E0 for the product
        # channel species if necessary
        for products in self.products:
            for spec in products.species:
                if spec.conformer is None:
                    spec.conformer = Conformer(E0=spec.get_thermo_data().E0)

        # Use the lowest E0 as the reference energy (`energy_correction`) for the network
        # The `energy_correction` will be added to the free energies and enthalpies for each
        # configuration in the network.
        energy_correction = -min(sum(spec.conformer.E0.value_si for spec in stationary_point.species)
                                 for stationary_point in self.reactants + self.isomers + self.products)
        for spec in self.reactants + self.products + self.isomers:
            spec.energy_correction = energy_correction
        self.energy_correction = energy_correction

        # Determine transition state energies on potential energy surface
        # In the absence of any better information, we simply set it to
        # be the reactant ground-state energy + the activation energy
        # Note that we need Arrhenius kinetics in order to do this
        for rxn in self.path_reactions:
            if rxn.kinetics is None:
                raise Exception('Path reaction "{0}" in PDepNetwork #{1:d} has no kinetics!'.format(rxn, self.index))
            elif isinstance(rxn.kinetics, KineticsData):
                if len(rxn.reactants) == 1:
                    kunits = 's^-1'
                elif len(rxn.reactants) == 2:
                    kunits = 'm^3/(mol*s)'
                elif len(rxn.reactants) == 3:
                    kunits = 'm^6/(mol^2*s)'
                else:
                    kunits = ''
                rxn.kinetics = Arrhenius().fit_to_data(Tlist=rxn.kinetics.Tdata.value_si,
                                                       klist=rxn.kinetics.kdata.value_si, kunits=kunits)
            elif isinstance(rxn.kinetics, MultiArrhenius):
                logging.info('Converting multiple kinetics to a single Arrhenius expression for reaction {rxn}'.format(
                    rxn=rxn))
                rxn.kinetics = rxn.kinetics.to_arrhenius(Tmin=Tmin, Tmax=Tmax)
            elif not isinstance(rxn.kinetics, Arrhenius) and rxn.network_kinetics is None:
                raise Exception('Path reaction "{0}" in PDepNetwork #{1:d} has invalid kinetics '
                                'type "{2!s}".'.format(rxn, self.index, rxn.kinetics.__class__))
            rxn.fix_barrier_height(force_positive=True)
            if rxn.network_kinetics is None:
                E0 = sum([spec.conformer.E0.value_si for spec in rxn.reactants]) + rxn.kinetics.Ea.value_si + energy_correction
            else:
                E0 = sum([spec.conformer.E0.value_si for spec in rxn.reactants]) + rxn.network_kinetics.Ea.value_si + energy_correction
            rxn.transition_state = rmgpy.species.TransitionState(conformer=Conformer(E0=(E0 * 0.001, "kJ/mol")))

        # Set collision model
        bath_gas = [spec for spec in reaction_model.core.species if not spec.reactive]
        assert len(bath_gas) > 0, 'No unreactive species to identify as bath gas'

        self.bath_gas = {}
        for spec in bath_gas:
            # is this really the only/best way to weight them?
            self.bath_gas[spec] = 1.0 / len(bath_gas)

        # Save input file
        if not self.label:
            self.label = str(self.index)

        if output_directory:
            job.save_input_file(
                os.path.join(output_directory, 'pdep', 'network{0:d}_{1:d}.py'.format(self.index, len(self.isomers))))
            if getattr(pdep_settings, 'generate_PES_diagrams', False):
                job.draw(os.path.join(output_directory, 'pdep'), filename_stem=f'network{self.index:d}_{len(self.isomers):d}', file_format='pdf')

        # Calculate the rate coefficients
        self.initialize(Tmin, Tmax, Pmin, Pmax, maximum_grain_size, minimum_grain_count, active_j_rotor, active_k_rotor,
                        rmgmode)
        try:
            K = self.calculate_rate_coefficients(Tlist, Plist, method)
        except InvalidMicrocanonicalRateError:
            if output_directory:
                filename_stem = f'network{self.index:d}_{len(self.isomers):d}'
                job.draw(output_directory, filename_stem=filename_stem, file_format='pdf')
                logging.info(f"Network {self.index} has been drawn and saved as {filename_stem}.pdf in {output_directory} to aid debugging.")
            raise

        # Generate PDepReaction objects
        configurations = []
        configurations.extend([isom.species[:] for isom in self.isomers])
        configurations.extend([reactant.species[:] for reactant in self.reactants])
        configurations.extend([product.species[:] for product in self.products])
        j = configurations.index(self.source)

        for i in range(K.shape[2]):
            if i != j:
                # Find the path reaction
                net_reaction = None
                for r in self.net_reactions:
                    if r.has_template(configurations[j], configurations[i]):
                        net_reaction = r
                # If net reaction does not already exist, make a new one
                if net_reaction is None:
                    net_reaction = PDepReaction(
                        reactants=configurations[j],
                        products=configurations[i],
                        network=self,
                        kinetics=None,
                    )
                    net_reaction = reaction_model.make_new_pdep_reaction(net_reaction)
                    self.add_net_reaction(net_reaction)

                    # Place the net reaction in the core or edge if necessary
                    # Note that leak reactions are not placed in the edge
                    if all([s in reaction_model.core.species for s in net_reaction.reactants]) \
                            and all([s in reaction_model.core.species for s in net_reaction.products]):
                        # Check whether netReaction already exists in the core as a LibraryReaction
                        for rxn in reaction_model.core.reactions:
                            if isinstance(rxn, LibraryReaction) \
                                    and rxn.is_same_reaction(net_reaction, either_direction=True) \
                                    and not rxn.allow_pdep_route \
                                    and (rxn.kinetics.is_pressure_dependent() or not rxn.elementary_high_p):
                                logging.info(f'Network reaction {net_reaction} matched an existing core reaction {rxn} '
                                             f'from the {rxn.library} library, and was not added to the model')
                                break
                        else:
                            reaction_model.add_reaction_to_core(net_reaction, requires_rms=requires_rms)
                    else:
                        # Check whether netReaction already exists in the edge as a LibraryReaction
                        for rxn in reaction_model.edge.reactions:
                            if isinstance(rxn, LibraryReaction) \
                                    and rxn.is_same_reaction(net_reaction, either_direction=True) \
                                    and not rxn.allow_pdep_route \
                                    and (rxn.kinetics.is_pressure_dependent() or not rxn.elementary_high_p):
                                logging.info(f'Network reaction {net_reaction} matched an existing edge reaction {rxn} '
                                             f'from the {rxn.library} library, and was not added to the model')
                                break
                        else:
                            reaction_model.add_reaction_to_edge(net_reaction, requires_rms=requires_rms)

                # Set/update the net reaction kinetics using interpolation model
                kdata = K[:, :, i, j].copy()
                order = len(net_reaction.reactants)
                kdata *= 1e6 ** (order - 1)
                kunits = {1: 's^-1', 2: 'cm^3/(mol*s)', 3: 'cm^6/(mol^2*s)'}[order]
                net_reaction.kinetics = job.fit_interpolation_model(Tlist, Plist, kdata, kunits)
                net_reaction.check_resolved_species_reversibility()

                # Check: For each net reaction that has a path reaction, make
                # sure the k(T,P) values for the net reaction do not exceed
                # the k(T) values of the path reaction
                # Only check the k(T,P) value at the highest P and lowest T,
                # as this is the one most likely to be in the high-pressure 
                # limit
                t = 0
                p = len(Plist) - 1
                for pathReaction in self.path_reactions:
                    if pathReaction.is_isomerization():
                        # Don't check isomerization reactions, since their
                        # k(T,P) values potentially contain both direct and
                        # well-skipping contributions, and therefore could be
                        # significantly larger than the direct k(T) value
                        # (This can also happen for association/dissociation
                        # reactions, but the effect is generally not too large)
                        continue
                    if pathReaction.reactants == net_reaction.reactants and pathReaction.products == net_reaction.products:
                        if pathReaction.network_kinetics is not None:
                            kinf = pathReaction.network_kinetics.get_rate_coefficient(Tlist[t])
                        else:
                            kinf = pathReaction.kinetics.get_rate_coefficient(Tlist[t])
                        if K[t, p, i, j] > 2 * kinf:  # To allow for a small discretization error
                            logging.warning('k(T,P) for net reaction {0} exceeds high-P k(T) by {1:g} at {2:g} K, '
                                            '{3:g} bar'.format(net_reaction, K[t, p, i, j] / kinf, Tlist[t], Plist[p] / 1e5))
                            logging.info('    k(T,P) = {0:9.2e}    k(T) = {1:9.2e}'.format(K[t, p, i, j], kinf))
                        break
                    elif pathReaction.products == net_reaction.reactants and pathReaction.reactants == net_reaction.products:
                        pathReaction.check_resolved_species_reversibility(reversible=True)
                        if pathReaction.network_kinetics is not None:
                            pathReaction.check_resolved_species_reversibility(
                                kinetics=pathReaction.network_kinetics, reversible=True)
                            kinf = pathReaction.network_kinetics.get_rate_coefficient(
                                Tlist[t]) / pathReaction.get_equilibrium_constant(Tlist[t])
                        else:
                            kinf = pathReaction.kinetics.get_rate_coefficient(
                                Tlist[t]) / pathReaction.get_equilibrium_constant(Tlist[t])
                        if K[t, p, i, j] > 2 * kinf:  # To allow for a small discretization error
                            logging.warning('k(T,P) for net reaction {0} exceeds high-P k(T) by {1:g} at {2:g} K, '
                                            '{3:g} bar'.format(net_reaction, K[t, p, i, j] / kinf, Tlist[t], Plist[p] / 1e5))
                            logging.info('    k(T,P) = {0:9.2e}    k(T) = {1:9.2e}'.format(K[t, p, i, j], kinf))
                        break

        self.log_summary(level=logging.INFO)

        # Delete intermediate arrays to conserve memory
        self.cleanup()

        self.remove_products_from_reactants()

        # We're done processing this network, so mark it as valid
        self.valid = True
