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

import os.path
import logging
import time
import string

import numpy as np
import yaml
try:
    from yaml import CDumper as Dumper, CLoader as Loader, CSafeLoader as SafeLoader
except ImportError:
    from yaml import Dumper, Loader, SafeLoader

from rmgpy.rmgpyobject import RMGPYObject
from rmgpy import __version__ as version
import rmgpy.constants as constants
from rmgpy.quantity import ScalarQuantity, ArrayQuantity
from rmgpy.molecule.element import elementList
from rmgpy.molecule.translator import toInChI, toInChIKey
from rmgpy.statmech.conformer import Conformer
from rmgpy.statmech.rotation import LinearRotor, NonlinearRotor, KRotor, SphericalTopRotor
from rmgpy.statmech.torsion import HinderedRotor, FreeRotor
from rmgpy.statmech.translation import IdealGasTranslation
from rmgpy.statmech.vibration import HarmonicOscillator
from rmgpy.pdep.collision import SingleExponentialDown
from rmgpy.transport import TransportData
from rmgpy.thermo import NASA, Wilhoit
from rmgpy.cantherm.pdep import PressureDependenceJob

################################################################################


def check_conformer_energy(Vlist,path):
    """
    Check to see that the starting energy of the species in the potential energy scan calculation
    is not 0.5 kcal/mol (or more) higher than any other energies in the scan. If so, print and 
    log a warning message.  
    """    
    Vlist = np.array(Vlist, np.float64)
    Vdiff = (Vlist[0] - np.min(Vlist))*constants.E_h*constants.Na/1000
    if Vdiff >= 2:  # we choose 2 kJ/mol to be the critical energy
        logging.warning('the species corresponding to ' + str(os.path.basename(path)) +
                        ' is different in energy from the lowest energy conformer by ' + "%0.2f" % Vdiff +
                        ' kJ/mol. This can cause significant errors in your computed rate constants. ')


def is_pdep(jobList):
    for job in jobList:
        if isinstance(job, PressureDependenceJob):
            return True
    return False


################################################################################


class CanthermSpecies(RMGPYObject):
    """
    A class for parsing Cantherm species with statMech data into .yml files
    """
    def __init__(self, species=None, conformer=None, author='', level_of_theory='', model_chemistry='',
                 frequency_scale_factor=None, use_hindered_rotors=None, use_bond_corrections=None, atom_energies='',
                 chemkin_thermo_string='', SMILES=None, adjacency_list=None, InChI=None, InChI_Key=None, xyz=None,
                 molecular_weight=None, symmetryNumber=None, transport_data=None, energy_transfer_model=None,
                 thermo=None, label=None, datetime=None, cantherm_version=None):
        if species is None and conformer is None:
            # Expecting to get a `species` when generating the object within Cantherm,
            # or a `conformer` when parsing from YAML.
            raise ValueError('No species or conformer was passed to the CanthermSpecies object')
        if conformer is not None:
            self.conformer = conformer
        if label is None and species is not None:
            self.label = species.label
        else:
            self.label = label
        self.author = author
        self.level_of_theory = level_of_theory
        self.model_chemistry = model_chemistry
        self.frequency_scale_factor = frequency_scale_factor
        self.use_hindered_rotors = use_hindered_rotors
        self.use_bond_corrections = use_bond_corrections
        self.atom_energies = atom_energies
        self.chemkin_thermo_string = chemkin_thermo_string
        self.SMILES = SMILES
        self.adjacency_list = adjacency_list
        self.InChI = InChI
        self.InChI_Key = InChI_Key
        self.xyz = xyz
        self.molecular_weight = molecular_weight
        self.symmetryNumber = symmetryNumber
        self.transport_data = transport_data
        self.energy_transfer_model = energy_transfer_model  # check pdep flag
        self.thermo = thermo
        if species is not None:
            self.update_species_attributes(species)
        self.cantherm_version = version
        self.datetime = time.strftime("%Y-%m-%d %H:%M")

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the object
        """
        result = '{0!r}'.format(self.__class__.__name__)
        result += '{'
        for key, value in self.as_dict().iteritems():
            if key != 'class':
                result += '{0!r}: {1!r}'.format(str(key), str(value))
        result += '}'
        return result

    def update_species_attributes(self, species=None):
        """
        Update the object with a new species (while keeping non-species-dependent attributes unchanged)
        """
        if species is None:
            raise ValueError('No species was passed to CanthermSpecies')
        self.label = species.label
        if species.molecule is not None and len(species.molecule) > 0:
            self.SMILES = species.molecule[0].toSMILES()
            self.adjacency_list = species.molecule[0].toAdjacencyList()
            try:
                inchi = toInChI(species.molecule[0], backend='try-all', aug_level=0)
            except ValueError:
                inchi = ''
            try:
                inchi_key = toInChIKey(species.molecule[0], backend='try-all', aug_level=0)
            except ValueError:
                inchi_key = ''
            self.InChI = inchi
            self.InChI_Key = inchi_key
            if species.conformer is not None:
                self.conformer = species.conformer
                self.xyz = self.update_xyz_string()
            self.molecular_weight = species.molecularWeight
            if species.symmetryNumber != -1:
                self.symmetryNumber = species.symmetryNumber
            if species.transportData is not None:
                self.transport_data = species.transportData  # called `collisionModel` in Cantherm
            if species.energyTransferModel is not None:
                self.energy_transfer_model = species.energyTransferModel
            if species.thermo is not None:
                self.thermo = species.thermo.as_dict()

    def update_xyz_string(self):
        if self.conformer is not None and self.conformer.number is not None:
            # generate the xyz-format string from the Conformer coordinates
            xyz_string = '{0}\n{1}'.format(len(self.conformer.number.value_si), self.label)
            for i, coorlist in enumerate(self.conformer.coordinates.value_si):
                for element in elementList:
                    if element.number == int(self.conformer.number.value_si[i]):
                        element_symbol = element.symbol
                        break
                else:
                    raise ValueError('Could not find element symbol corresponding to atom number {0}'.format(
                        self.conformer.number.value_si[i]))
                xyz_string += '\n{0} {1} {2} {3}'.format(element_symbol,
                                                         coorlist[0],
                                                         coorlist[1],
                                                         coorlist[2])
        else:
            xyz_string = ''
        return xyz_string

    def save_yaml(self, path):
        """
        Save the species with all statMech data to a .yml file
        """
        if not os.path.exists(os.path.join(os.path.abspath(path),'SpeciesDatabase', '')):
            os.mkdir(os.path.join(os.path.abspath(path),'SpeciesDatabase', ''))
        valid_chars = "-_.()<=>+ %s%s" % (string.ascii_letters, string.digits)
        filename = os.path.join('SpeciesDatabase',
                                ''.join(c for c in self.label if c in valid_chars) + '.yml')
        full_path = os.path.join(path, filename)
        with open(full_path, 'w') as f:
            yaml.dump(data=self.as_dict(), stream=f, canonical=False)
        # remove empty lines from the file (multiline strings are have excess new line brakes for some reason):
        with open(full_path, 'r') as f:
            lines = f.readlines()
        with open(full_path, 'w') as f:
            for line in lines:
                if not line.isspace():
                    f.write(line)
        logging.debug('Dumping species {0} data as {1}'.format(self.label, filename))

    def load_yaml(self, path, species, pdep=False):
        """
        Load the all statMech data from the .yml file in `path` into `species`
        `pdep` is a boolean specifying whether or not jobList includes a pressureDependentJob.
        """
        logging.info('Loading statistical mechanics parameters for {0} from .yml file...'.format(species.label))
        with open(path, 'r') as f:
            data = yaml.safe_load(stream=f)
        try:
            if species.label != data['label']:
                logging.warning('Found different labels for species: {0} in input file, and {1} in the .yml file. '
                                'Using the label "{0}" for this species.'.format(species.label, data['label']))
        except KeyError:
            # Lacking label in the YAML file is strange, but accepted
            logging.debug('Did not find label for species {0} in .yml file.'.format(species.label))
        try:
            class_name = data['class']
        except KeyError:
            raise KeyError("Can only make objects if the `class` attribute in the dictionary is known")
        if class_name != 'CanthermSpecies':
            raise KeyError("Expected a CanthermSpecies object, but got {0}".format(class_name))
        del data['class']
        class_dict = {'ScalarQuantity': ScalarQuantity,
                      'ArrayQuantity': ArrayQuantity,
                      'Conformer': Conformer,
                      'LinearRotor': LinearRotor,
                      'NonlinearRotor': NonlinearRotor,
                      'KRotor': KRotor,
                      'SphericalTopRotor': SphericalTopRotor,
                      'HinderedRotor': HinderedRotor,
                      'FreeRotor': FreeRotor,
                      'IdealGasTranslation': IdealGasTranslation,
                      'HarmonicOscillator': HarmonicOscillator,
                      'TransportData': TransportData,
                      'SingleExponentialDown': SingleExponentialDown,
                      'Wilhoit': Wilhoit,
                      'NASA': NASA,
                      }
        self.make_object(data=data, class_dict=class_dict)
        if pdep and (self.transport_data is None or self.energy_transfer_model is None):
            raise ValueError('Transport data and an energy transfer model must be given if pressure-dependent '
                             'calculations are requested. Check file {0}'.format(path))
        if pdep and self.SMILES is None and self.adjacency_list is None\
                and self.InChI is None and self.molecular_weight is None:
            raise ValueError('The molecular weight was not specified, and a structure was not given so it could '
                             'not be calculated. Specify either the molecular weight or structure if '
                             'pressure-dependent calculations are requested. Check file {0}'.format(path))
        logging.debug("Parsed all YAML objects")
