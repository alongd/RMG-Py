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

import csv
import json
import os

from rmgpy.chemkin import get_species_identifier
from rmgpy.tools.plot import SimulationPlot


def _minimum_confinement(entries):
    """Decode the manifest's positive-infinity sentinel before numeric use."""
    return min(float('inf') if item['conf'] == 'infinite' else item['conf']
               for item in entries.values())


class SimulationProfileWriter(object):
    """
    SimulationProfileWriter listens to a ReactionSystem subject
    and writes the species mole numbers as a function of the reaction time
    to a csv file.


    A new instance of the class can be appended to a subject as follows:
    
    reaction_system = ...
    listener = SimulationProfileWriter()
    reaction_system.attach(listener)

    Whenever the subject calls the .notify() method, the
    .update() method of the listener will be called.

    To stop listening to the subject, the class can be detached
    from its subject:

    reaction_system.detach(listener)

    """

    def __init__(self, output_directory, reaction_sys_index, core_species):
        super(SimulationProfileWriter, self).__init__()

        self.output_directory = output_directory
        self.reaction_sys_index = reaction_sys_index
        self.core_species = core_species

    def update(self, reaction_system):
        """
        Opens a file with filename referring to:
            - reaction system
            - number of core species

        Writes to a csv file:
            - header row with species names
            - each row with number of moles of the core species in the given reaction system.
        """

        filename = os.path.join(
            self.output_directory,
            'solver',
            'simulation_{0}_{1:d}.csv'.format(
                self.reaction_sys_index + 1, len(self.core_species)
            )
        )

        from rmgpy.export import SpeciesReferences, resolve_species_reference
        wall_model = getattr(reaction_system, 'electronegative_wall_model', None)
        declarations = SpeciesReferences(self.core_species, get_species_identifier, context='simulation profile',
                                         allow_ground_collisions=True)
        header = ['Time (s)', 'Volume (m^3)']
        for spc in self.core_species:
            header.append(resolve_species_reference(spc, declarations))

        records = None
        if wall_model is not None:
            # Each profile row uses its own accepted record, never the final
            # gate values repeated over the trajectory.
            records = [reaction_system.electronegative_wall_history[row[0]]
                       for row in reaction_system.snapshots]
            manifest = reaction_system.electronegative_wall_manifest()
            identity = '{0}; {1}'.format(manifest['closure'], manifest['geometry_arm'])
            header.extend(['Wall h [{0}]'.format(identity),
                           'Wall floating potential e/(kTe)',
                           'Wall gate A minimum conf', 'Wall gate B minimum metric',
                           'Wall gate C radial error', 'Wall gate C geometry error',
                           'Wall gate C full-profile error'])
            with open(filename[:-4] + '.electronegative-wall.json', 'w') as stream:
                json.dump(manifest, stream, indent=2, allow_nan=False)
                stream.write('\n')

        with open(filename, 'w') as csvfile:
            worksheet = csv.writer(csvfile)

            # add header row:
            worksheet.writerow(header)

            # add mole fractions:
            if records is None:
                worksheet.writerows(reaction_system.snapshots)
            else:
                for row, record in zip(reaction_system.snapshots, records):
                    gates = record['gates']
                    potential = record['floating_potential_e_over_kTe']
                    values = [record['h'], potential if potential is not None else float('nan')]
                    if isinstance(gates, dict):
                        values.extend([_minimum_confinement(gates['A']),
                                       min(gates['B']['metrics']), gates['C_radial']['error'],
                                       gates['C_geometry'].get('error', float('nan')),
                                       gates['C_full_profile'].get('error', float('nan'))])
                    else:
                        values.extend([float('nan')] * 5)
                    worksheet.writerow(list(row) + values)


class SimulationProfilePlotter(object):
    """
    SimulationProfilePlotter listens to a ReactionSystem subject
    and plots the top 10 species mole fraction profiles.

    A new instance of the class can be appended to a subject as follows:
    
    reaction_system = ...
    listener = SimulationProfilePlotter()
    reaction_system.attach(listener)

    Whenever the subject calls the .notify() method, the
    .update() method of the listener will be called.

    To stop listening to the subject, the class can be detached
    from its subject:

    reaction_system.detach(listener)
    """

    def __init__(self, output_directory, reaction_sys_index, core_species):
        super(SimulationProfilePlotter, self).__init__()

        self.output_directory = output_directory
        self.reaction_sys_index = reaction_sys_index
        self.core_species = core_species

    def update(self, reaction_system):
        """
        Saves a png with filename referring to:
            - reaction system
            - number of core species
        """

        csv_file = os.path.join(
            self.output_directory,
            'solver',
            'simulation_{0}_{1:d}.csv'.format(
                self.reaction_sys_index + 1, len(self.core_species)
            )
        )

        png_file = os.path.join(
            self.output_directory,
            'solver',
            'simulation_{0}_{1:d}.png'.format(
                self.reaction_sys_index + 1, len(self.core_species)
            )
        )

        title = ''
        if getattr(reaction_system, 'electronegative_wall_model', None) is not None:
            manifest = reaction_system.electronegative_wall_manifest()
            gates = manifest['gates']
            if isinstance(gates, dict):
                gate_labels = []
                for name, result in gates.items():
                    if name == 'A':
                        value = 'min conf={0:.4g}'.format(_minimum_confinement(result))
                    elif name == 'B':
                        value = 'min metric={0:.4g}'.format(min(result['metrics']))
                    elif 'error' in result:
                        value = '{0:.4g}/{1:.4g}'.format(result['error'], result['threshold'])
                    else:
                        value = result['status']
                    gate_labels.append(name + ': ' + value)
                gate_text = '\n'.join(gate_labels)
            else:
                gate_text = gates
            phi = manifest['floating_potential_e_over_kTe']
            title = '{0}; {1}\nt={2:.4g} s; phi={3}\n{4}'.format(
                manifest['closure'], manifest['geometry_arm'], manifest['time'],
                '{0:.4g}'.format(phi) if phi is not None else 'unavailable', gate_text)
        SimulationPlot(csv_file=csv_file, num_species=10, ylabel='Moles', title=title).plot(png_file)
