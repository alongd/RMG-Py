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
import rmgpy.constants as constants
from rmgpy.species import Species
import yaml
import re

from rmgpy.molecule.translator import toInChI, toInChIKey
from rmgpy.thermo import NASA, Wilhoit
from rmgpy.statmech.conformer import Conformer
from rmgpy.statmech import IdealGasTranslation, LinearRotor, NonlinearRotor, HarmonicOscillator
from rmgpy.statmech.torsion import HinderedRotor, FreeRotor
from rmgpy.pdep.collision import SingleExponentialDown
from rmgpy.quantity import ScalarQuantity
try:
    from yaml import CDumper as Dumper, CLoader as Loader
except ImportError:
    from yaml import Dumper, Loader


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


# def yaml_constructors():
#     yaml.add_constructor(u'scalar_quantity', scalar_quantity_constructor)
#
#
# def yaml_representers(species):
#     yaml.add_representer(Conformer, species.conformer.conformer_representer)
#     yaml.add_representer(SingleExponentialDown, species.energyTransferModel.single_exponential_down_representer)
#     yaml.add_representer(ScalarQuantity, scalar_quantity_representer)
#
#
# def yaml_resolvers():
#     conformer_pattern = re.compile('Conformer')
#     yaml.add_implicit_resolver(u'!conformer', conformer_pattern)
#     sed_pattern = re.compile('SingleExponentialDown')
#     yaml.add_implicit_resolver(u'!single_exponential_down', sed_pattern)
#     scalar_q_pattern = re.compile('!scalar_quantity')
#     yaml.add_implicit_resolver(u'!scalar_quantity', scalar_q_pattern)
#
#
# def scalar_quantity_representer(Dumper, data):
#     """
#     A helper function for YAML parsing of the ScalarQuantity object
#     """
#     return Dumper.represent_scalar(u'!scalar_quantity', data.__repr__())
#
#
# def scalar_quantity_constructor(Loader, node):
#     """
#     A helper function for YAML parsing of the ScalarQuantity object
#     """
#     uncertainty = None
#     uncertaintyType = '+|-'
#     data = Loader.construct_scalar(node).replace('(', '').replace(')', '').split(',')
#     value = float(data[0])
#     units = data[1]
#     if len(data) > 2:
#         uncertainty = float(data[2])
#     if len(data) > 3:
#         uncertaintyType = data[3]
#     return ScalarQuantity(value=value, units=units, uncertainty=uncertainty, uncertaintyType=uncertaintyType)


def yaml_constructors(cantherm_species):
    yaml.add_constructor(u'CanthermSpecies', cantherm_species.cantherm_species_constructor)



class CanthermSpecies(yaml.YAMLObject):
    """
    A class for parsing Cantherm species .yml files
    """

    yaml_tag = u'!CanthermSpecies'

    def __init__(self, species=None, author='', level_of_theory='', model_chemistry='',
                 frequency_scale_factor=None, use_hindered_rotors=None, use_atom_corrections=None,
                 use_bond_corrections=None, atom_energies='', chemkin_thermo_string=''):
        if species is None:
            raise ValueError('No species was passed to CanthermSpecies')

        dictionary = dict()

        dictionary['label'] = species.label
        dictionary['date'] = time.strftime("%Y-%m-%d %H:%M")
        dictionary['author'] = author
        dictionary['level_of_theory'] = level_of_theory
        dictionary['model_chemistry'] = model_chemistry
        dictionary['frequency_scale_factor'] = frequency_scale_factor
        dictionary['use_hindered_rotors'] = use_hindered_rotors
        dictionary['use_atom_corrections'] = use_atom_corrections
        dictionary['use_bond_corrections'] = use_bond_corrections
        dictionary['atom_energies'] = atom_energies

        inchi = ''
        inchi_key = ''
        if species.molecule is not None and len(species.molecule) > 0:
            dictionary['smiles'] = species.molecule[0].toSMILES()
            dictionary['adjacencyList'] = species.molecule[0].toAdjacencyList()
            try:
                inchi = toInChI(species.molecule[0], backend='try-all', aug_level=0)
            except ValueError:
                pass
            try:
                inchi_key = toInChIKey(species.molecule[0], backend='try-all', aug_level=0)
            except ValueError:
                pass
            dictionary['inchi'] = inchi
            dictionary['inchi_key'] = inchi_key

        species_dict = dict()

        # species.conformer
        if species.conformer is not None and species.conformer.E0 is not None:
            conformer_dict = dict()
            conformer_dict['spin_multiplicity'] = species.conformer.spinMultiplicity
            conformer_dict['optical_isomers'] = species.conformer.opticalIsomers
            conformer_dict['E0'] = [species.conformer.E0.value, species.conformer.E0.units]
            # species.conformer.modes
            modes = dict()
            for mode in species.conformer.modes:
                mode_dict = dict()
                if isinstance(mode, IdealGasTranslation):
                    mode_dict['mass'] = [mode.mass.value, mode.mass.units]
                if isinstance(mode, LinearRotor):
                    if mode.inertia is not None:
                        mode_dict['inertia'] = [[float(value) for value in mode.inertia.value], mode.inertia.units]
                    if mode.rotationalConstant is not None:
                        mode_dict['rotationalConstant'] = [[float(value) for value in mode.rotationalConstant.value],
                                                           mode.inertia.units]
                    mode_dict['symmetry'] = mode.symmetry
                    mode_dict['quantum'] = mode.quantum
                if isinstance(mode, NonlinearRotor):
                    if mode.inertia is not None:
                        mode_dict['inertia'] = [[float(value) for value in mode.inertia.value], mode.inertia.units]
                    if mode.rotationalConstant is not None:
                        mode_dict['rotationalConstant'] = [[float(value) for value in mode.rotationalConstant.value],
                                                           mode.inertia.units]
                    mode_dict['symmetry'] = mode.symmetry
                    mode_dict['quantum'] = mode.quantum
                if isinstance(mode, HarmonicOscillator):
                    if mode.frequencies is not None:
                        mode_dict['frequencies'] = [[float(value) for value in mode.frequencies.value],
                                                    mode.frequencies.units]
                    if mode.quantum is not None:
                        mode_dict['quantum'] = mode.quantum
                if isinstance(mode, HinderedRotor):
                    if mode.inertia is not None:
                        mode_dict['inertia'] = [mode.inertia.value, mode.inertia.units]
                    if mode.rotationalConstant is not None:
                        mode_dict['rotationalConstant'] = [mode.rotationalConstant.value, mode.rotationalConstant.units]
                    mode_dict['symmetry'] = mode.symmetry
                    if mode.fourier is not None:
                        mode_dict['fourier'] = [[[float(value) for value in array] for array in mode.fourier.value],
                                                mode.fourier.units]
                    if mode.barrier is not None:
                        mode_dict['barrier'] = [mode.barrier.value, mode.barrier.units]
                    mode_dict['quantum'] = mode.quantum
                    mode_dict['semiclassical'] = mode.semiclassical
                if isinstance(mode, FreeRotor):
                    if mode.inertia is not None:
                        mode_dict['inertia'] = [[float(value) for value in mode.inertia.value], mode.inertia.units]
                    if mode.rotationalConstant is not None:
                        mode_dict['rotationalConstant'] = [mode.rotationalConstant.value, mode.rotationalConstant.units]
                modes['{0}'.format(mode.__class__.__name__)] = mode_dict
            conformer_dict['modes'] = modes
            species_dict['conformer'] = conformer_dict

        # species.thermo
        if species.thermo is not None:
            thermo_dict = dict()
            if isinstance(species.thermo, NASA):
                thermo_dict['type'] = 'NASA'
                poly_list = []
                poly_list.append(dict())
                poly_list[0]['coeffs'] = [float(value) for value in species.thermo.polynomials[0].coeffs]
                poly_list[0]['Tmin'] = [species.thermo.polynomials[0].Tmin.value,
                                                    species.thermo.polynomials[0].Tmin.units]
                poly_list[0]['Tmax'] = [species.thermo.polynomials[0].Tmax.value,
                                                    species.thermo.polynomials[0].Tmax.units]
                if species.thermo.polynomials[0].E0 is not None:
                    poly_list[0]['E0'] = [species.thermo.polynomials[0].E0.value,
                                                    species.thermo.polynomials[0].E0.units]
                if species.thermo.polynomials[0].comment is not '':
                    poly_list[0]['comment'] = species.thermo.polynomials[0].comment
                poly_list.append(dict())
                poly_list[1]['coeffs'] = [float(value) for value in species.thermo.polynomials[1].coeffs]
                poly_list[1]['Tmin'] = [species.thermo.polynomials[1].Tmin.value,
                                                    species.thermo.polynomials[1].Tmin.units]
                poly_list[1]['Tmax'] = [species.thermo.polynomials[1].Tmax.value,
                                                    species.thermo.polynomials[1].Tmax.units]
                if species.thermo.polynomials[1].E0 is not None:
                    poly_list[1]['E0'] = [species.thermo.polynomials[1].E0.value,
                                                    species.thermo.polynomials[1].E0.units]
                if species.thermo.polynomials[1].comment is not '':
                    poly_list[1]['comment'] = species.thermo.polynomials[1].comment
                thermo_dict['polynomials'] = poly_list
                if species.thermo.Tmin is not None:
                    thermo_dict['Tmin'] = [species.thermo.Tmin.value, species.thermo.Tmin.units]
                if species.thermo.Tmax is not None:
                    thermo_dict['Tmax'] = [species.thermo.Tmax.value, species.thermo.Tmax.units]
                if species.thermo.E0 is not None:
                    thermo_dict['E0'] = [species.thermo.E0.value, species.thermo.E0.units]
                if species.thermo.Cp0 is not None:
                    thermo_dict['Cp0'] = [species.thermo.Cp0.value, species.thermo.Cp0.units]
                if species.thermo.CpInf is not None:
                    thermo_dict['CpInf'] = [species.thermo.CpInf.value, species.thermo.CpInf.units]
                if species.thermo.label is not '':
                    thermo_dict['label'] = species.thermo.label
                if species.thermo.comment is not '':
                    thermo_dict['comment'] = species.thermo.comment
            elif isinstance(species.thermo, Wilhoit):
                thermo_dict['type'] = 'Wilhoit'
                thermo_dict['a0'] = species.thermo.a0
                thermo_dict['a1'] = species.thermo.a1
                thermo_dict['a2'] = species.thermo.a2
                thermo_dict['a3'] = species.thermo.a3
                if species.thermo.H0 is not None:
                    thermo_dict['H0'] = species.thermo.H0
                if species.thermo.S0 is not None:
                    thermo_dict['S0'] = species.thermo.S0
                if species.thermo.B is not None:
                    thermo_dict['B'] = [species.thermo.B.value, species.thermo.B.units]
                if species.thermo.Tmin is not None:
                    thermo_dict['Tmin'] = [species.thermo.Tmin.value, species.thermo.Tmin.units]
                if species.thermo.Tmax is not None:
                    thermo_dict['Tmax'] = [species.thermo.Tmax.value, species.thermo.Tmax.units]
                if species.thermo.E0 is not None:
                    thermo_dict['E0'] = [species.thermo.E0.value, species.thermo.E0.units]
                if species.thermo.Cp0 is not None:
                    thermo_dict['Cp0'] = [species.thermo.Cp0.value, species.thermo.Cp0.units]
                if species.thermo.CpInf is not None:
                    thermo_dict['CpInf'] = [species.thermo.CpInf.value, species.thermo.CpInf.units]
                if species.thermo.label is not '':
                    thermo_dict['label'] = species.thermo.label
                if species.thermo.comment is not '':
                    thermo_dict['comment'] = species.thermo.comment
            species_dict['thermo'] = thermo_dict

        if species.transportData is not None:
            transport_data_dict = dict()
            transport_data_dict['epsilon'] = [species.transportData.epsilon.value, species.transportData.epsilon.units]
            transport_data_dict['sigma'] = [species.transportData.sigma.value, species.transportData.sigma.units]
            species_dict['transport_data'] = transport_data_dict

        if species.molecularWeight is not None:
            species_dict['molecular_weight'] = [species.molecularWeight.value, species.molecularWeight.units]

        if species.energyTransferModel is not None:
            energy_transfer_model_dict = dict()
            if isinstance(species.energyTransferModel, SingleExponentialDown):
                energy_transfer_model_dict['type'] = 'SingleExponentialDown'
                energy_transfer_model_dict['alpha0'] = [species.energyTransferModel.alpha0.value,
                                                          species.energyTransferModel.alpha0.units]
                energy_transfer_model_dict['T0'] = [species.energyTransferModel.T0.value,
                                                      species.energyTransferModel.T0.units]
                energy_transfer_model_dict['n'] = species.energyTransferModel.n
            species_dict['energy_transfer_model'] = energy_transfer_model_dict

        species_dict['symmetry_number'] = species.symmetryNumber

        dictionary['species'] = species_dict

        dictionary['chemkin_thermo_string'] = chemkin_thermo_string

        self.species_data = dictionary

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the object
        """
        result = '{0!r}(species_data{'.format(self.__class__.__name__)
        result += 'label: {0!r}, date: {1!r}, author: {2!r}, level_of_theory: {3!r}, ' \
                  'model_chemistry: {4!r}, frequency_scale_factor: {5!r}, use_hindered_rotors: {6!r}, ' \
                  'use_atom_corrections: {7!r}, use_bond_corrections: {8!r}, atom_energies: {9!r}, smiles: {10!r}, ' \
                  'adjacencyList: {11!r}, inchi: {12!r}, inchi_key: {13!r}, species: {14!r},' \
                  'chemkin_thermo_string: {15!r})'.format(
                    self.species_data['label'],
                    self.species_data['date'],
                    self.species_data['author'],
                    self.species_data['level_of_theory'],
                    self.species_data['model_chemistry'],
                    self.species_data['frequency_scale_factor'],
                    self.species_data['use_hindered_rotors'],
                    self.species_data['use_atom_corrections'],
                    self.species_data['use_bond_corrections'],
                    self.species_data['atom_energies'],
                    self.species_data['smiles'],
                    self.species_data['adjacencyList'],
                    self.species_data['inchi'],
                    self.species_data['inchi_key'],
                    self.species_data['species'],
                    self.species_data['chemkin_thermo_string'])
        result += '}'
        return result

    # def __reduce__(self):
    #     """
    #     A helper function used when pickling the object
    #     """
    #     return (CanthermSpecies, (self.label, self.date, self.author, self.level_of_theory,self.model_chemistry,
    #                 self.frequency_scale_factor, self.use_hindered_rotors, self.use_atom_corrections,
    #                 self.use_bond_corrections, self.smiles, self.adjacencyList, self.inchi, self.inchi_key,
    #                 self.species, self.chemkin_thermo_string))

    def update(self, species):
        """
        Update the object with a new species while keeping other attributes unchanged
        """
        self.__init__(species=species,
                      author=self.species_data['author'],
                      level_of_theory=self.species_data['level_of_theory'],
                      model_chemistry=self.species_data['model_chemistry'],
                      frequency_scale_factor=self.species_data['frequency_scale_factor'],
                      use_hindered_rotors=self.species_data['use_hindered_rotors'],
                      use_atom_corrections=self.species_data['use_atom_corrections'],
                      use_bond_corrections=self.species_data['use_bond_corrections'],
                      chemkin_thermo_string=self.species_data['chemkin_thermo_string'])

    def cantherm_species_constructor(self, Loader, node):
        """
        A helper function for YAML parsing of the object
        """
        value = Loader.construct_scalar(node)
        print value

        return CanthermSpecies()


