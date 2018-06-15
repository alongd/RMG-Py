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
import rmgpy.constants as constants
from rmgpy.quantity import Quantity
from rmgpy.molecule.element import getElement
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

