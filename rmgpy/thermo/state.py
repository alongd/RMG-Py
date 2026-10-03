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

"""Library-only thermochemistry for resolved states and declared v=0 species."""

from rmgpy.exceptions import ExcitedSpeciesThermoError, VibrationalManifoldError


def require_thermo_estimation_allowed(species_or_molecule):
    """Refuse state-blind estimation before a caller mutates or converts a graph."""
    molecules = getattr(species_or_molecule, 'molecule', [species_or_molecule])
    label = getattr(species_or_molecule, 'label', '')
    if species_or_molecule.props.get('vibrational_manifold'):
        raise VibrationalManifoldError(
            "Species {0!r} is declared by vibrationalManifold and requires a loaded "
            "thermo library entry for the same graph with vibrationallevel 0; "
            "other thermo sources are refused.".format(label or species_or_molecule.props['vibrational_manifold']))
    for molecule in molecules:
        if molecule.props.get('vibrational_manifold'):
            raise VibrationalManifoldError(
                "Species {0!r} is declared by vibrationalManifold; thermo must come "
                "from a library entry with vibrationallevel 0.".format(label or molecule.props['vibrational_manifold']))
        if molecule.has_resolved_state():
            raise ExcitedSpeciesThermoError(
                "Library-only thermochemistry is required for resolved-state species {0!r}. "
                "A loaded thermo library must contain an entry for this exact state; "
                "QM, HBI, ML and group additivity cannot estimate its thermo.\n{1}"
                .format(label or molecule.to_adjacency_list().strip(), molecule.to_adjacency_list()))


def thermo_library_species(species):
    """Return the library lookup structure without changing model identity."""
    if not species.props.get('vibrational_manifold'):
        return species
    lookup = species.copy(deep=True)
    lookup.props.pop('vibrational_manifold', None)
    for molecule in lookup.molecule:
        molecule.props.pop('vibrational_manifold', None)
        molecule.vibrational_level = 0
    return lookup


def require_vibrational_manifold_thermo(species, thermo_database):
    """Check attached v=0 thermo using the plasma library value comparator."""
    from rmgpy.solver.plasma import (
        _build_charged_thermo_library_index, _charged_species_library_thermo_match,
    )
    if thermo_database is None:
        raise VibrationalManifoldError(
            "Species {0!r} declared by vibrationalManifold requires a loaded thermo "
            "library entry with vibrationallevel 0.".format(species.label))
    outcome = _charged_species_library_thermo_match(
        species, _build_charged_thermo_library_index(thermo_database), {})
    if outcome['match'] is None:
        raise VibrationalManifoldError(
            "Thermo for species {0!r} declared by vibrationalManifold does not "
            "value-match any loaded thermo library entry for the same graph "
            "with vibrationallevel 0; other thermo is refused.".format(species.label))
