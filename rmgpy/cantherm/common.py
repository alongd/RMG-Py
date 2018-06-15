#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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

import numpy
import os.path
import logging
import time
import yaml
import re

import rmgpy.constants as constants
from rmgpy.molecule.translator import toInChI, toInChIKey
from rmgpy.thermo import NASA, NASAPolynomial, Wilhoit, ThermoData
from rmgpy.statmech.conformer import Conformer
from rmgpy.statmech import IdealGasTranslation, LinearRotor, NonlinearRotor, HarmonicOscillator
from rmgpy.statmech.torsion import HinderedRotor, FreeRotor
from rmgpy.pdep.collision import SingleExponentialDown
from rmgpy.thermo.nasa import NASA
from rmgpy.thermo.wilhoit import Wilhoit
from rmgpy.transport import TransportData
from rmgpy.quantity import Quantity
from rmgpy.molecule.element import getElement
try:
    from yaml import CDumper as Dumper, CLoader as Loader, CSafeLoader as SafeLoader
except ImportError:
    from yaml import Dumper, Loader, SafeLoader


################################################################################


def check_conformer_energy(Vlist,path):
    """
    Check to see that the starting energy of the species in the potential energy scan calculation
    is not 0.5 kcal/mol (or more) higher than any other energies in the scan. If so, print and 
    log a warning message.  
    """    
    Vlist = numpy.array(Vlist, numpy.float64)
    Vdiff = (Vlist[0] - numpy.min(Vlist))*constants.E_h*constants.Na/1000
    if Vdiff >= 2: #we choose 2 kJ/mol to be the critical energy
        logging.warning('the species corresponding to ' + str(os.path.basename(path)) +
                        ' is different in energy from the lowest energy conformer by ' + "%0.2f" % Vdiff +
                        ' kJ/mol. This can cause significant errors in your computed rate constants. ')


def determine_molecular_weight(species):
    """
    Determine Species.molecularWeight
    """
    if species.molecule is None or len(species.molecule) == 0:
        raise ValueError('Cannot determine molecularWeight for structureless species')
    return Quantity(sum([getElement(int(atom.number)).mass for atom in species.molecule[0].atoms]), 'kg/mol')


################################################################################


def yaml_registers():
    yaml_representers()
    yaml_constructors()
    yaml_implicit_resolvers()


def yaml_representers():
    """
    Register YAML representers
    (Representers convert objects to tagged scalar nodes)
    """
    yaml.add_representer(Conformer, conformer_representer)
    yaml.add_representer(SingleExponentialDown, single_exponential_down_representer)
    yaml.add_representer(NASA, nasa_representer)
    yaml.add_representer(Wilhoit, wilhoit_representer)
    yaml.add_representer(ThermoData, thermodata_representer)
    yaml.add_representer(TransportData, transport_data_representer)


def yaml_constructors():
    """
    Register YAML constructors
    (Constructors rebuild objects)
    """
    yaml.SafeLoader.add_constructor(u'!conformer', conformer_constructor)
    yaml.SafeLoader.add_constructor(u'!single_exponential_down', single_exponential_down_constructor)
    yaml.SafeLoader.add_constructor(u'!nasa', nasa_constructor)
    yaml.SafeLoader.add_constructor(u'!wilhoit', wilhoit_constructor)
    yaml.SafeLoader.add_constructor(u'!thermodata', thermodata_constructor)
    yaml.SafeLoader.add_constructor(u'!transport_data', transport_data_constructor)


def yaml_implicit_resolvers():
    """
    Define and register YAML implicit resolvers
    (Implicit resolvers detect the `!object` tags using regex without explicitly writing the tags in the .yml file)
    Two implicit resolvers are registered per tag, one for dump/load and the other for safe_load
    """
    pattern = re.compile(r'Conformer.*')
    yaml.add_implicit_resolver(u'!conformer', pattern)
    yaml.SafeLoader.add_implicit_resolver(u'!conformer', pattern, first='Conformer')
    pattern = re.compile(r'SingleExponentialDown.*')
    yaml.add_implicit_resolver(u'!single_exponential_down', pattern)
    yaml.SafeLoader.add_implicit_resolver(u'!single_exponential_down', pattern, first='SingleExponentialDown')
    pattern = re.compile(r'NASA.*')
    yaml.add_implicit_resolver(u'!nasa', pattern)
    yaml.SafeLoader.add_implicit_resolver(u'!nasa', pattern, first='NASA')
    pattern = re.compile(r'Wilhoit.*')
    yaml.add_implicit_resolver(u'!wilhoit', pattern)
    yaml.SafeLoader.add_implicit_resolver(u'!wilhoit', pattern, first='Wilhoit')
    pattern = re.compile(r'ThermoData.*')
    yaml.add_implicit_resolver(u'!thermodata', pattern)
    yaml.SafeLoader.add_implicit_resolver(u'!thermodata', pattern, first='ThermoData')
    pattern = re.compile(r'TransportData.*')
    yaml.add_implicit_resolver(u'!transport_data', pattern)
    yaml.SafeLoader.add_implicit_resolver(u'!transport_data', pattern, first='TransportData')


def conformer_representer(dumper, data):
    """
    A helper function for YAML parsing of Conformer objects
    """
    return dumper.represent_scalar(u'!conformer', data.__repr__())


def conformer_constructor(loader, node):
    """
    A helper function for YAML parsing of Conformer objects
    """
    value = loader.construct_scalar(node)
    if value.startswith('Conformer(E0=('):
        return eval(value)
    else:
        raise ValueError('Could not interpret Conformer object')


def single_exponential_down_representer(dumper, data):
    """
    A helper function for YAML parsing of SingleExponentialDown objects
    """
    return dumper.represent_scalar(u'!single_exponential_down', data.__repr__())


def single_exponential_down_constructor(loader, node):
    """
    A helper function for YAML parsing of SingleExponentialDown objects
    """
    value = loader.construct_scalar(node)
    if value.startswith('SingleExponentialDown(alpha0=('):
        return eval(value)
    else:
        raise ValueError('Could not interpret SingleExponentialDown object')


def nasa_representer(dumper, data):
    """
    A helper function for YAML parsing of NASA objects
    """
    return dumper.represent_scalar(u'!nasa', data.__repr__())


def nasa_constructor(loader, node):
    """
    A helper function for YAML parsing of NASA objects
    """
    value = loader.construct_scalar(node)
    if value.startswith('NASA(polynomials=['):
        return eval(value)
    else:
        raise ValueError('Could not interpret NASA object')


def wilhoit_representer(dumper, data):
    """
    A helper function for YAML parsing of Wilhoit objects
    """
    return dumper.represent_scalar(u'!nasa', data.__repr__())


def wilhoit_constructor(loader, node):
    """
    A helper function for YAML parsing of Wilhoit objects
    """
    value = loader.construct_scalar(node)
    if value.startswith('Wilhoit(Cp0=('):
        return eval(value)
    else:
        raise ValueError('Could not interpret Wilhoit object')


def thermodata_representer(dumper, data):
    """
    A helper function for YAML parsing of ThermoData objects
    """
    return dumper.represent_scalar(u'!thermodata', data.__repr__())


def thermodata_constructor(loader, node):
    """
    A helper function for YAML parsing of ThermoData objects
    """
    value = loader.construct_scalar(node)
    if value.startswith('ThermoData('):
        return eval(value)
    else:
        raise ValueError('Could not interpret ThermoData object')


def transport_data_representer(dumper, data):
    """
    A helper function for YAML parsing of TransportData objects
    """
    return dumper.represent_scalar(u'!transport_data', data.__repr__())


def transport_data_constructor(loader, node):
    """
    A helper function for YAML parsing of TransportData objects
    """
    value = loader.construct_scalar(node)
    if value.startswith('TransportData('):
        return eval(value)
    else:
        raise ValueError('Could not interpret TransportData object')


class CanthermSpecies(object):
    """
    A class for parsing Cantherm species .yml files
    """
    def __init__(self, species=None, author='', level_of_theory='', model_chemistry='',
                 frequency_scale_factor=None, use_hindered_rotors=None, use_atom_corrections=None,
                 use_bond_corrections=None, atom_energies='', chemkin_thermo_string=''):
        if species is None:
            raise ValueError('No species was passed to CanthermSpecies')

        d = dict()

        d['label'] = species.label
        d['datetime'] = time.strftime("%Y-%m-%d %H:%M")
        d['author'] = author
        d['level_of_theory'] = level_of_theory
        d['model_chemistry'] = model_chemistry
        d['frequency_scale_factor'] = frequency_scale_factor
        d['use_hindered_rotors'] = use_hindered_rotors
        d['use_atom_corrections'] = use_atom_corrections
        d['use_bond_corrections'] = use_bond_corrections
        d['atom_energies'] = atom_energies

        if species.molecule is not None and len(species.molecule) > 0:
            d['SMILES'] = species.molecule[0].toSMILES()
            d['adjacencyList'] = species.molecule[0].toAdjacencyList()
            try:
                inchi = toInChI(species.molecule[0], backend='try-all', aug_level=0)
            except ValueError:
                inchi = ''
            try:
                inchi_key = toInChIKey(species.molecule[0], backend='try-all', aug_level=0)
            except ValueError:
                inchi_key = ''
            d['InChI'] = inchi
            d['InChI_Key'] = inchi_key

        d['conformer'] = species.conformer

        if species.symmetryNumber != -1:
            d['symmetryNumber'] = species.symmetryNumber
        if species.transportData is not None:
            d['transport_data'] = species.transportData  # collisionModel
        if species.molecularWeight is not None:
            d['molecularWeight'] = [species.molecularWeight.value, species.molecularWeight.units]  # decode
        if species.energyTransferModel is not None:
            d['energy_transfer_model'] = species.energyTransferModel

        if species.thermo is not None:
            d['thermo'] = species.thermo
        d['chemkin_thermo_string'] = chemkin_thermo_string

        self.species_dictionary = d

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the object
        """
        result = '{0!r}'.format(self.__class__.__name__)
        result += '{'
        result += 'label: {0!r}, date: {1!r}, author: {2!r}, level_of_theory: {3!r}, ' \
                  'model_chemistry: {4!r}, frequency_scale_factor: {5!r}, use_hindered_rotors: {6!r}, ' \
                  'use_atom_corrections: {7!r}, use_bond_corrections: {8!r}, atom_energies: {9!r}, SMILES: {10!r}, ' \
                  'adjacencyList: {11!r}, InChI: {12!r}, InChI_Key: {13!r}, conformer: {14!r}, ' \
                  'transport_data: {15!r}, molecularWeight: {16!r}, energy_transfer_model: {17!r}, ' \
                  'thermo: {18!r}, chemkin_thermo_string: {19!r}'.format(
                    self.species_dictionary['label'],
                    self.species_dictionary['datetime'],
                    self.species_dictionary['author'],
                    self.species_dictionary['level_of_theory'],
                    self.species_dictionary['model_chemistry'],
                    self.species_dictionary['frequency_scale_factor'],
                    self.species_dictionary['use_hindered_rotors'],
                    self.species_dictionary['use_atom_corrections'],
                    self.species_dictionary['use_bond_corrections'],
                    self.species_dictionary['atom_energies'],
                    self.species_dictionary['SMILES'],
                    self.species_dictionary['adjacencyList'],
                    self.species_dictionary['InChI'],
                    self.species_dictionary['InChI_Key'],
                    self.species_dictionary['conformer'],
                    self.species_dictionary['transport_data'],
                    self.species_dictionary['molecularWeight'],
                    self.species_dictionary['energy_transfer_model'],
                    self.species_dictionary['thermo'],
                    self.species_dictionary['chemkin_thermo_string'],
        )
        result += '}'
        return result

    def __reduce__(self):
        """
        A helper function used when pickling the object
        """
        return (CanthermSpecies, (self.species_dictionary))

    def update(self, species):
        """
        Update the object with a new species while keeping other attributes unchanged
        """
        self.__init__(species=species,
                      author=self.species_dictionary['author'],
                      level_of_theory=self.species_dictionary['level_of_theory'],
                      model_chemistry=self.species_dictionary['model_chemistry'],
                      frequency_scale_factor=self.species_dictionary['frequency_scale_factor'],
                      use_hindered_rotors=self.species_dictionary['use_hindered_rotors'],
                      use_atom_corrections=self.species_dictionary['use_atom_corrections'],
                      use_bond_corrections=self.species_dictionary['use_bond_corrections'],
                      chemkin_thermo_string=self.species_dictionary['chemkin_thermo_string'],
                      )
