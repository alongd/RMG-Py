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

"""
This module contains the :class:`Sensitivity` class
"""

import os
import os.path
import sys
import logging
import argparse
import time
try:
    import matplotlib
    matplotlib.rc('mathtext', default='regular')
except ImportError:
    pass

from rmgpy.chemkin import writeElementsSection, writeThermoEntry, writeKineticsEntry

from rmgpy.cantherm.input import loadInputFile

from rmgpy.cantherm.kinetics import KineticsJob
from rmgpy.cantherm.statmech import StatMechJob
from rmgpy.cantherm.thermo import ThermoJob
from rmgpy.cantherm.pdep import PressureDependenceJob

################################################################################


class Sensitivity(object):
    """
    The :class:`Sensitivity` class represents an instance of a sensitivity analysis job.
    The attributes are:
    
    =================== =============================================================================
    Attribute           Description
    =================== =============================================================================
    `conditions`        A list of the conditions at which the sensitivity coefficients are calculated
    `inputFile`         The path of the input file defining the jobs to execute
    `outputDirectory`   The directory in which to write the output files
    `verbose`           The level of detail in the generated logging messages
    =================== =============================================================================
    """
    
    def __init__(self):
        pass


class KineticsSensitivity(object):
    """
    The :class:`Sensitivity` class represents an instance of a sensitivity analysis job
    performed for a KineticsJob. The attributes are:

    =================== =============================================================================
    Attribute           Description
    =================== =============================================================================
    `conditions`        A list of the conditions at which the sensitivity coefficients are calculated
    `jobList`           A list of the jobs to execute
    `inputFile`         The path of the input file defining the jobs to execute
    `outputDirectory`   The directory in which to write the output files
    `verbose`           The level of detail in the generated logging messages
    =================== =============================================================================

    The output directory defaults to the same directory as the input file if
    not explicitly specified.

    To use this class programmatically, create an instance and set its
    attributes using either the :meth:`__init__()` method or by directly
    accessing the attributes, and then invoke the :meth:`execute()` method.
    You can also populate the attributes from the command line using the
    :meth:`parseCommandLineArguments()` method before running :meth:`execute()`.
    """

    def __init__(self, kinetics, conditions=None):
        self.conditions = conditions if conditions is not None else [(1000, 'K')]
        self.kinetics = kinetics





    def perturb_and_run_job(self, species):
        logging.debug("Perturbing E0 of species {0} by 5%...".format(self.sa_species))
        species.conformer.E0.value *= 1.05



        species.conformer.E0.value /= 1.05  # restore E0 to its original value


class PDepSensitivity(object):
    """
    The :class:`Sensitivity` class represents an instance of a sensitivity analysis job
    performed for a PressureDependenceJob. The attributes are:

    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `jobList`           A list of the jobs to execute
    `inputFile`         The path of the input file defining the jobs to execute
    `outputDirectory`   The directory in which to write the output files
    `verbose`           The level of detail in the generated logging messages
    =================== ========================================================

    The output directory defaults to the same directory as the input file if
    not explicitly specified.

    To use this class programmatically, create an instance and set its
    attributes using either the :meth:`__init__()` method or by directly
    accessing the attributes, and then invoke the :meth:`execute()` method.
    You can also populate the attributes from the command line using the
    :meth:`parseCommandLineArguments()` method before running :meth:`execute()`.
    """

    def __init__(self, inputFile=None, outputDirectory=None, verbose=logging.INFO):
        self.conditions = conditions if conditions is not None else [(1000, 'K'),(1, 'bar')]
        self.jobList = []
        self.inputFile = inputFile
        self.outputDirectory = outputDirectory
        self.verbose = verbose

































    if sensitivity:
        time_array = []
        normSens_array = [[] for spec in self.sensitiveSpecies]
        RTP = constants.R * self.T.value_si / self.P.value_si
        # identify sensitive species indices
        sensSpeciesIndices = numpy.array([speciesIndex[spec] for spec in self.sensitiveSpecies],
                                         numpy.int)  # index within coreSpecies list of the sensitive species

    if sensitivity:
        time_array.append(self.t)
        moleSens = self.y[numCoreSpecies:]  #
        volume = self.V

        dVdk = numpy.zeros(numCoreReactions + numCoreSpecies, numpy.float64)
        if not self.constantVolume:
            for j in xrange(numCoreReactions + numCoreSpecies):
                dVdk[j] = numpy.sum(
                    moleSens[j * numCoreSpecies:(j + 1) * numCoreSpecies]) * RTP  # Contains [ dV_dk and dV_dG ]
        for i in xrange(len(self.sensitiveSpecies)):
            normSens = numpy.zeros(numCoreReactions + numCoreSpecies, numpy.float64)
            c = self.coreSpeciesConcentrations[sensSpeciesIndices[i]]
            if c != 0:
                for j in xrange(numCoreReactions):
                    normSens[j] = 1 / volume * (moleSens[j * numCoreSpecies + sensSpeciesIndices[i]] - c * dVdk[j]) * \
                                  forwardRateCoefficients[j] / c
                for j in xrange(numCoreReactions, numCoreReactions + numCoreSpecies):
                    normSens[j] = 1 / volume * (moleSens[j * numCoreSpecies + sensSpeciesIndices[i]] - c * dVdk[
                        j]) / c * 4184  # no normalization against dG, converstion to kcal/mol units
            normSens_array[i].append(normSens)

    if sensitivity:
        for i in xrange(len(self.sensitiveSpecies)):
            with open(sensWorksheet[i], 'wb') as outfile:
                worksheet = csv.writer(outfile)
                reactionsAboveThreshold = []
                for j in xrange(numCoreReactions + numCoreSpecies):
                    for k in xrange(len(time_array)):
                        if abs(normSens_array[i][k][j]) > self.sensitivityThreshold:
                            reactionsAboveThreshold.append(j)
                            break
                species_name = getSpeciesIdentifier(self.sensitiveSpecies[i])
                headers = ['Time (s)']
                headers.extend(['dln[{0}]/dln[k{1}]: {2}'.format(species_name, j + 1, coreReactions[j].toChemkin(
                    kinetics=False)) if j < numCoreReactions
                                else 'dln[{0}]/dG[{1}]'.format(species_name,
                                                               getSpeciesIdentifier(coreSpecies[j - numCoreReactions]))
                                for j in reactionsAboveThreshold])
                worksheet.writerow(headers)

                for k in xrange(len(time_array)):
                    row = [time_array[k]]
                    row.extend([normSens_array[i][k][j] for j in reactionsAboveThreshold])
                    worksheet.writerow(row)








    def parseCommandLineArguments(self):
        """
        Parse the command-line arguments being passed to CanTherm. This uses the
        :mod:`argparse` module, which ensures that the command-line arguments are
        sensible, parses them, and returns them.
        """
    
        parser = argparse.ArgumentParser(description=
        """
        CanTherm is a Python toolkit for computing chemical reaction rates and other
        properties used in detailed kinetics models using various methodologies
        and theories.
        """)
        parser.add_argument('file', metavar='FILE', type=str, nargs=1,
            help='a file describing the job to execute')
    
        # Options for controlling the amount of information printed to the console
        # By default a moderate level of information is printed; you can either
        # ask for less (quiet), more (verbose), or much more (debug)
        group = parser.add_mutually_exclusive_group()
        group.add_argument('-q', '--quiet', action='store_const', const=logging.WARNING, default=logging.INFO, dest='verbose', help='only print warnings and errors')
        group.add_argument('-v', '--verbose', action='store_const', const=logging.DEBUG, default=logging.INFO, dest='verbose', help='print more verbose output')
        group.add_argument('-d', '--debug', action='store_const', const=0, default=logging.INFO, dest='verbose', help='print debug information')
    
        # Add options for controlling what directories files are written to
        parser.add_argument('-o', '--output-directory', type=str, nargs=1, default='',
            metavar='DIR', help='use DIR as output directory')

        # Add options for controlling generation of plots
        parser.add_argument('-p', '--plot', action='store_true', default=False, help='generate plots of results')

        args = parser.parse_args()
        
        # Extract the input file
        self.inputFile = args.file[0] 
        
        # Extract the log verbosity
        self.verbose = args.verbose
        
        # Extract the plot settings
        self.plot = args.plot
        
        # Determine the output directory
        # By default the directory containing the input file is used, unless an
        # alternate directory is specified using the -o flag
        if args.output_directory and os.path.isdir(args.output_directory[0]):
            self.outputDirectory = os.path.abspath(args.output_directory[0])
        else:
            self.outputDirectory = os.path.dirname(os.path.abspath(args.file[0]))
    
    def initializeLog(self, verbose=logging.INFO, logFile=None):
        """
        Set up a logger for CanTherm to use to print output to stdout. The
        `verbose` parameter is an integer specifying the amount of log text seen
        at the console; the levels correspond to those of the :data:`logging` module.
        """
        # Create logger
        logger = logging.getLogger()
        logger.setLevel(verbose)
    
        # Use custom level names for cleaner log output
        logging.addLevelName(logging.CRITICAL, 'Critical: ')
        logging.addLevelName(logging.ERROR, 'Error: ')
        logging.addLevelName(logging.WARNING, 'Warning: ')
        logging.addLevelName(logging.INFO, '')
        logging.addLevelName(logging.DEBUG, '')
        logging.addLevelName(0, '')
    
        # Create formatter and add to handlers
        formatter = logging.Formatter('%(levelname)s%(message)s')
        
        # Remove old handlers before adding ours
        while logger.handlers:
            logger.removeHandler(logger.handlers[0])
       
        # Create console handler; send everything to stdout rather than stderr
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(verbose)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    
        # Create file handler; always be at least verbose in the file
        if logFile:
            fh = logging.FileHandler(filename=logFile)
            fh.setLevel(min(logging.DEBUG,verbose))
            fh.setFormatter(formatter)
            logger.addHandler(fh)
    
    def logHeader(self, level=logging.INFO):
        """
        Output a header containing identifying information about CanTherm to the log.
        """
        from rmgpy import __version__
        logging.log(level, 'CanTherm execution initiated at {0}'.format(time.asctime()))
        logging.log(level, '')
    
        logging.log(level, '###############################################################')
        logging.log(level, '#                                                             #')
        logging.log(level, '#                          CanTherm                           #')
        logging.log(level, '#                                                             #')
        logging.log(level, '#   Version: {0:48s} #'.format(__version__))
        logging.log(level, '#   Authors: RMG Developers (rmg_dev@mit.edu)                 #')
        logging.log(level, '#   P.I.s:   William H. Green (whgreen@mit.edu)               #')
        logging.log(level, '#            Richard H. West (r.west@neu.edu)                 #')
        logging.log(level, '#   Website: http://reactionmechanismgenerator.github.io/     #')
        logging.log(level, '#                                                             #')
        logging.log(level, '###############################################################')
        logging.log(level, '')
    
    def logFooter(self, level=logging.INFO):
        """
        Output a footer to the log.
        """
        logging.log(level, '')
        logging.log(level, 'CanTherm execution terminated at {0}'.format(time.asctime()))

    def loadInputFile(self, inputFile):
        """
        Load a set of jobs from the given `inputFile` on disk. Returns the
        loaded set of jobs as a list.
        """
        self.inputFile = inputFile
        self.jobList = loadInputFile(self.inputFile)
        logging.info('')
        return self.jobList
        
    def execute(self):
        """
        Execute, in order, the jobs found in input file specified by the
        `inputFile` attribute.
        """
        
        # Initialize the logging system (both to the console and to a file in the
        # output directory)
        self.initializeLog(self.verbose, os.path.join(self.outputDirectory, 'cantherm.log'))
        
        # Print some information to the beginning of the log
        self.logHeader()
        
        # Load the input file for the job
        self.jobList = self.loadInputFile(self.inputFile)
        logging.info('')
        
        # Initialize (and clear!) the output files for the job
        if self.outputDirectory is None:
            self.outputDirectory = os.path.dirname(os.path.abspath(self.inputFile))
        outputFile = os.path.join(self.outputDirectory, 'output.py')
        with open(outputFile, 'w') as f:
            pass
        chemkinFile = os.path.join(self.outputDirectory, 'chem.inp')

        # write the chemkin files and run the thermo and then kinetics jobs
        with open(chemkinFile, 'w') as f:
            writeElementsSection(f)
            
            f.write('SPECIES\n\n')

            # write each species in species block
            for job in self.jobList:
                if isinstance(job,ThermoJob):
                    f.write(job.species.toChemkin())
                    f.write('\n')

            f.write('\nEND\n\n\n\n')
            f.write('THERM ALL\n')
            f.write('    300.000  1000.000  5000.000\n\n')

        # run thermo jobs (printing out thermo stuff)
        for job in self.jobList:
            if isinstance(job,ThermoJob) or isinstance(job, StatMechJob):
                job.execute(outputFile=outputFile, plot=self.plot)

        with open(chemkinFile, 'a') as f:
            f.write('\n')
            f.write('END\n\n\n\n')
            f.write('REACTIONS    KCAL/MOLE   MOLES\n\n')

        # run kinetics and pdep jobs (and outputing chemkin stuff)
        for job in self.jobList:
            if isinstance(job,KineticsJob) or isinstance(job, PressureDependenceJob):
                job.execute(outputFile=outputFile, plot=self.plot)

        with open(chemkinFile, 'a') as f:
            f.write('END\n\n')

        # Print some information to the end of the log
        self.logFooter()
