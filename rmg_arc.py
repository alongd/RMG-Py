#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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
RMG is an automatic chemical mechanism generator. It is awesomely awesome.
ARC is an automatic electronic structure computation scheduler. It is magnificently magnificent.
Here we iteratively execute RMG and ARC to generate and refine a chemical kinetic model.
"""

import os
import shutil
import logging
import pandas as pd

import yaml

import rmgpy
from rmgpy.rmg.main import RMG, initializeLog, processProfileStats, makeProfileGraph
from rmgpy.chemkin import loadChemkinFile
from rmgpy import settings
from rmgpy.data.rmg import getDB
from rmgpy.data.thermo import ThermoLibrary
from rmgpy.thermo import NASAPolynomial, NASA, ThermoData, Wilhoit

try:
    from arc import ARC
    from arc.settings import valid_chars
except ImportError:
    is_arc_available = False
else:
    is_arc_available = True

################################################################################


def main(args):
    """
    Run the master modeler module, calling RMG and ARC as needed
    """
    input_file = args.file
    output_directory = args.output_directory
    kwargs = {
        'restart': args.restart,
        'walltime': args.walltime,
        'kineticsdatastore': args.kineticsdatastore
    }
    unconverged_species = list()
    # extract ARC's part of the input file:
    rmg_input_file, arc_arguments = parse_input_file(input_file)
    additional_calcs_required = True
    max_iterations = 100
    i = 0
    while additional_calcs_required or i == max_iterations:
        run_directory = os.path.join(output_directory, 'iteration_{0}'.format(i))
        if not os.path.exists(run_directory):
            os.mkdir(run_directory)
        # run RNG
        new_rmg_input_path = os.path.join(run_directory, 'rmg_input.py')
        with open(new_rmg_input_path, 'w') as f:
            f.write(rmg_input_file)
        print('\n\n\n    Running RMG, iteration {0}:\n\n\n'.format(i))
        run_rmg(input_file=new_rmg_input_path, output_directory=run_directory, args=args, kwargs=kwargs)

        # determine species to calculate
        species_to_calc, additional_calcs_required = determine_species_to_calculate(run_directory, arc_arguments,
                                                                                    unconverged_species)
        if additional_calcs_required:
            # run ARC
            print('\n\n\n  Running ARC (iteration {0}) for the following species:'.format(i))
            for spc in species_to_calc:
                print('{0}: {1}'.format(spc.label, spc.molecule[0].toSMILES()))
            print('\n\n\n')
            run_arc(arc_arguments['kwargs'], run_directory, species_to_calc)
            # add the calculated RMG libraries to the database and input file
            unconverged_species.extend(get_unconverged_species(run_directory))
            rmg_input_file = add_rmg_libraries(run_directory, rmg_input_file, i)
        i += 1
    # clear the RMG log files in the parent directory - they don't contain actual info
    for item in os.listdir(output_directory):
        if item.endswith('.log'):
            os.remove(os.path.join(output_directory, item))


def run_rmg(input_file, output_directory, args, kwargs):
    # Parse the command-line arguments (requires the argparse module)

    if args.postprocess:
        print("Postprocessing the profiler statistics (will be appended to RMG.log)")
    else:
        # Initialize the logging system (resets the RMG.log file)
        level = logging.INFO
        if args.debug:
            level = 0
        elif args.verbose:
            level = logging.DEBUG
        elif args.quiet:
            level = logging.WARNING
        initializeLog(level, os.path.join(output_directory, 'RMG.log'))

    logging.info(rmgpy.settings.report())

    if args.profile:
        import cProfile
        global_vars = {}
        local_vars = {
            'inputFile': args.file,
            'output_dir': output_directory,
            'kwargs': kwargs,
            'RMG': RMG
        }

        command = """rmg = RMG(inputFile=inputFile, outputDirectory=output_dir); rmg.execute(**kwargs)"""

        stats_file = os.path.join(output_directory, 'RMG.profile')
        print("Running under cProfile")
        if not args.postprocess:
            # actually run the program!
            cProfile.runctx(command, global_vars, local_vars, stats_file)
        # postprocess the stats
        log_file = os.path.join(output_directory, 'RMG.log')
        processProfileStats(stats_file, log_file)
        makeProfileGraph(stats_file)

    else:
        rmg = RMG(inputFile=input_file, outputDirectory=output_directory)
        rmg.execute(**kwargs)


def run_arc(kwargs, run_directory, species_to_calc):
    """
    Run ARC with the given kwargs in `arc_arguments`
    Calculate thermo from all species in `species_to_calc`, which is a list of RMG Species objects
    ARC Project directory is `run_directory`/ARC
    """
    kwargs['project_directory'] = os.path.join(run_directory, 'ARC')
    if not os.path.exists(kwargs['project_directory']):
        os.mkdir(kwargs['project_directory'])
    i = 0
    for spc in species_to_calc:
        for char in spc.label:
            if char not in valid_chars:
                spc.label = 'unique_ARC_label_{0}'.format(i)
                i += 1
    kwargs['arc_species_list'] = species_to_calc  # ARC accepts RMG Species objects
    kwargs['project'] = 'rmg_arc'
    arc0 = ARC(**kwargs)
    arc0.execute()


# def clear_rmg_database():
#     """Set the RMG database object instant to None"""
#     rmg_db = getDB()
#     rmg_db.thermo = None
#     rmg_db.transport = None
#     rmg_db.forbiddenStructures = None
#     rmg_db.kinetics = None
#     rmg_db.statmech = None
#     rmg_db.solvation = None


def parse_input_file(input_file):
    """
    Extract ARC-related parameters from the input file,
    returns:
    `rmg_input_file` ,the legal RMG input file section, a multiline string
    `arc_arguments`, a dictionary of ARC's parameters

    The ARC section of the input file (arc_arguments) has the following structure (no leading spaces):
arc(
'SA species': 10,
'SA reactions': 10,
'SA pdep': True,
'all core species': False,
'SA observables': True
'collision violators': True
'kwargs': {dict of ARC's kwargs}
)
    The above ARC arguments default to `False` or 0 if not given.
    `kwargs` is an optional dictionary of key word arguments to pass to the ARC object (e.g., level of theory, etc.)
    """
    with open(input_file, 'r') as f:
        lines = f.readlines()
    read_rmg = True
    rmg_input_file = ''
    arc_arguments  = ''
    for line in lines:
        if 'arc(' in line:
            read_rmg = False
            continue
        if read_rmg:
            rmg_input_file += line
        else:
            if line == ')\n':
                read_rmg = True
                continue
            arc_arguments += line
    arc_arguments = yaml.safe_load(arc_arguments)
    return rmg_input_file, arc_arguments


def determine_species_to_calculate(run_directory, arc_arguments, unconverged_species):
    """
    Determine which species in the executed RMG job located in the `run_directory` path
    should be calculated by ARC using the specified criteria in the `arc_arguments` dictionary
    parsed from the arc() section of the input file.
    """
    species_to_calc = list()  # contains RMG Species objects
    spcs_labels_to_calc = list()  # contains labels only
    sa_spcs_to_calc = list()  # contains labels only
    pdep_rxns_to_explore = list()  # contains RMG Reaction objects

    chemkin_path = os.path.join(run_directory, 'chemkin', 'chem_annotated.inp')
    dictionary_path = os.path.join(run_directory, 'chemkin', 'species_dictionary.txt')
    rmg_species, rmg_reactions = loadChemkinFile(chemkin_path, dictionary_path)

    if 'all core species' in arc_arguments and arc_arguments['all core species']:
        # we'll calculate all core species, if needed based on their thermo comment
        # don't bother checking sensitivities, collision rate violators, etc.
        for spc in rmg_species:
            if calc_based_on_thermo_comment(spc) and spc.label not in unconverged_species:
                species_to_calc.append(spc)
        print_summary(run_directory, species_to_calc, header='All core species:')
    else:
        if 'SA species' in arc_arguments or 'SA reactions' in arc_arguments\
                or 'SA pdep' in arc_arguments or 'SA observables' in arc_arguments:
            sa_path = os.path.join(run_directory, 'solver')
            if not os.path.exists(sa_path):
                print("\n\nError: Could not find path to RMG's sensitivity analysis output. Not executing calculations "
                      "in ARC based on sensitivity analysis!\n\n")
            else:
                sa_files = list()
                for file1 in os.listdir(sa_path):
                    if file1.endswith(".csv"):
                        sa_files.append(file1)
                for sa_file in sa_files:
                    # iterate through all SA .csv files in the solver folder
                    df = pd.read_csv(sa_file)
                    sa_dict = dict()
                    sa_dict['rxn'], sa_dict['spc'] = dict(), dict()
                    for header in df.columns:
                        # iterate through all headers in the SA .csv file, but skip the `Time (s)` column
                        sa_type = None
                        if 'dlnk' in header and 'SA reactions' in arc_arguments and arc_arguments['SA reactions']:
                            sa_type = 'rxn'
                        elif 'dG' in header and 'SA species' in arc_arguments and arc_arguments['SA species']:
                            sa_type = 'spc'
                        if sa_type is not None:
                            # proceed only if we care about this column
                            entry = dict()
                            is_pdep = False
                            observable = header.split('[')[1].split(']')[0]
                            if observable not in spcs_labels_to_calc and observable not in unconverged_species\
                                    and 'SA observables' in arc_arguments and arc_arguments['SA observables']:
                                spcs_labels_to_calc.append(observable)
                                print_summary(run_directory, observable, header='observable')
                            if observable not in sa_dict[sa_type]:
                                sa_dict[sa_type][observable] = list()
                            # parameter extraction examples:
                            # for species get `C2H4(8)` from `dln[ethane(1)]/dG[C2H4(8)]`
                            # for reaction, get int::8 from `dln[ethane(1)]/dln[k8]: H(6)+ethane(1)=H2(12)+C2H5(5)`
                            parameter = header.split('[')[2].split(']')[0]
                            if sa_type == 'rxn':
                                parameter = int(header.split('[')[2].split(']')[0][1:])
                                if rmg_reactions[parameter].kinetics.isPressureDependent():
                                    is_pdep = True
                            entry['parameter'] = parameter  # rxn number or spc label
                            entry['is_pdep'] = is_pdep
                            entry['max_sa'] = max(df[header].max(), abs(df[header].min()))
                            sa_dict[sa_type][observable].append(entry)
                    num_rxn = arc_arguments['SA reactions'] if 'SA reactions' in arc_arguments else 0
                    num_spc = arc_arguments['SA species'] if 'SA species' in arc_arguments else 0
                    for _, sa_list in sa_dict['rxn'].iteritems():
                        # we don't use the `observable` key here, but we do need to iterate observables separately
                        sa_list_sorted = sorted(sa_list, key=lambda item: item['max_sa'], reverse=True)
                        for i in xrange(num_rxn):
                            for spc in rmg_reactions[sa_list_sorted[i]['parameter']].reactants\
                                    + rmg_reactions[sa_list_sorted[i]['parameter']].products:
                                if spc.label not in sa_spcs_to_calc and calc_based_on_thermo_comment(spc):
                                    sa_spcs_to_calc.append(spc)
                                    print_summary(run_directory, spc, header='from sensitive reactions')
                            if sa_list_sorted[i]['is_pdep']\
                                    and rmg_reactions[sa_list_sorted[i]['parameter']] not in pdep_rxns_to_explore\
                                    and 'SA pdep' in arc_arguments and arc_arguments['SA pdep']:
                                pdep_rxns_to_explore.append(rmg_reactions[sa_list_sorted[i]['parameter']])
                    for _, sa_list in sa_dict['spc'].iteritems():
                        sa_list_sorted = sorted(sa_list, key=lambda item: item['max_sa'], reverse=True)
                        for i in xrange(num_spc):
                            if sa_list_sorted[i]['parameter'] not in spcs_labels_to_calc\
                                    and sa_list_sorted[i]['parameter'] not in unconverged_species:
                                spcs_labels_to_calc.append(sa_list_sorted[i]['parameter'])
                                print_summary(run_directory, sa_list_sorted[i]['parameter'],
                                              header='from sensitive species')

        for rxn in pdep_rxns_to_explore:
            for spc in rxn.network.reactants + rxn.network.isomers + rxn.network.products:
                if spc not in sa_spcs_to_calc and calc_based_on_thermo_comment(spc)\
                        and spc.label not in unconverged_species:
                    species_to_calc.append(spc)
                    print_summary(run_directory, spc, header='from sensitive reaction networks')

        coll_violators_path = os.path.join(run_directory, 'collision_rate_violators.log')
        with open(coll_violators_path, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if '<=>' in line:
                labels = line.split()
                for label in labels:
                    if label != '<=>' and label not in spcs_labels_to_calc and label not in unconverged_species:
                        spcs_labels_to_calc.append(label)
                        print_summary(run_directory, label, header='from collision violators')

        # convert labels to RMG Species objects
        for label in spcs_labels_to_calc::
            for spc in rmg_species:
                if spc.label == label:
                    if spc not in sa_spcs_to_calc and calc_based_on_thermo_comment(spc):
                        species_to_calc.append(spc)
                    break

    additional_calcs_required = bool(len(species_to_calc))
    return species_to_calc, additional_calcs_required


def calc_based_on_thermo_comment(spc):
    """A helper function for reading the Species `spc` thermo comment and determining whether to calculate it in ARC"""
    if 'group additivity' in spc.thermo.comment or '+ radical(' in spc.thermo.comment:
        return True
    return False


def print_summary(run_directory, species, header):
    """
    Report the species to be calculated
    The report will be saved as `species_to_calculate_in_ARC.log` in `run_directory`, the RMG run folder
    `species_list` could be either a list or one instance of either RMG Species or labels corresponding to RMG Species
    `header` is the reason these species were chosen for calculation
    """
    path = os.path.join(run_directory, 'species_to_calculate_in_ARC.log')
    with open(path, 'w'):
        # create the file if it doesn't exist, overwrite the file if it exists
        pass
    if isinstance(species, list):
        species_labels = list()
        for spc in species:
            label = spc if isinstance(spc, (str, unicode)) else spc.label
            species_labels.append(label)
        with open(path, 'a') as f:
            f.write(header + '\n')
            for label in species_labels:
                f.write(label + '\n')
    else:
        label = species if isinstance(species, (str, unicode)) else species.label
        with open(path, 'a') as f:
            f.write('{0} ({1})\n'.format(label, header))


def add_rmg_libraries(run_directory, rmg_input_file, i):
    """
    Copy the RMG libraries generated by ARC (located in `run_directory`/ARC/output/RMG libraries/)
    to the RMG database, and add them to the RMG input file to be used in the next RMG run.
    A single library, `rmg_arc`, is used to which the ARC output is appended
    `i` is the RMG/ARC iteration number
    """
    new_rmg_input_file = ''
    arc_thermo_lib_path = os.path.join(run_directory, 'ARC', 'output', 'RMG libraries', 'thermo', 'rmg_arc.py')
    rmg_thermo_lib_path = os.path.join(settings['database.directory'], 'thermo', 'libraries')
    unique_library_name = False
    j = 0
    while not unique_library_name:
        if j:
            library_name = 'arc_thermo' + '_' + str(j) + '.py'
        else:
            library_name = 'arc_thermo.py'
        if not os.path.isfile(os.path.join(rmg_thermo_lib_path, library_name):
            unique_library_name = True





            
    local_context = {
        'ThermoData': ThermoData,
        'Wilhoit': Wilhoit,
        'NASAPolynomial': NASAPolynomial,
        'NASA': NASA,
    }
    if os.path.isfile(rmg_thermo_lib_path):
        # the rmg_arc thermo library already exists. load it, append, and save
        thermo_lib, arc_thermo_lib = ThermoLibrary(), ThermoLibrary()
        thermo_lib.load(path=rmg_thermo_lib_path, local_context=local_context, global_context=dict())
        arc_thermo_lib.load(path=arc_thermo_lib_path, local_context=local_context, global_context=dict())
        for entry in arc_thermo_lib.entries.itervalues():
            thermo_lib.entries[entry.label] = entry
        thermo_lib.save(path=rmg_thermo_lib_path)
        # TODO: if same label: check not isomorphic. Anyway, check not isomeprhpic with all other entries. different levels? so per project, add an index
    else:
        # the rmg_arc thermo library doesn't exist. just copy the lobrary generated by ARC
        shutil.copy(arc_thermo_lib_path, rmg_thermo_lib_path)
    for line in rmg_input_file.splitlines():
        if 'arc_thermo' in line:
            break
    else:
        # the arc_thermo library isn't present in the RMG input file, add it appropriately
        combine_lines = False
        combined_lines = list()
        for line in rmg_input_file.splitlines():
            if 'thermoLibraries' in line:
                if ']' in line.split('thermoLibraries')[1]:
                    new_rmg_input_file += add_new_library_to_line(line, 'arc_thermo') + '\n'
                else:
                    combine_lines = True
                    combined_lines.append(line)
            elif combine_lines:
                combined_lines.append(line)
                if ']' in line:
                    combine_lines = False
            else:
                new_rmg_input_file += line + '\n'

    return new_rmg_input_file or rmg_input_file


def add_new_library_to_line(line, library_name):
    """
    A helper function for adding a new thermo `library` into a single `line` that contains a list of libraries
    `line` could look like:
        thermoLibraries=['BurkeH2O2','primaryNS'],
    or:
        thermoLibraries=['BurkeH2O2','primaryNS'], reactionLibraries=['BurkeH2O2inN2'],
    or:
        seedMechanisms=[], thermoLibraries=['BurkeH2O2','primaryNS'], reactionLibraries=['BurkeH2O2inN2'],
    """
    if isinstance(line, list):
        line = ' '.join(l for l in line.splitlines())
    items = line.split('thermoLibraries')
    for i, char in enumerate(items[1]):
        if char == ']':
            break
    new_line = items[0] + "thermoLibraries" + items[1][:i] + ",'" + library_name + "'" + items[1][i:]
    return new_line


def get_unconverged_species(run_directory):
    """Get the labels of unconverged species from the present ARC iteration"""
    unconverged_species = list()
    unconverged_path = os.path.join(run_directory, 'ARC', 'output', 'unconverged_species.log')
    if os.path.isfile(unconverged_path):
        with open(unconverged_path, 'r') as f:
            for line in f:
                unconverged_species.append(line.split('\n')[0])
    return unconverged_species
