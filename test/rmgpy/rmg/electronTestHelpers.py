#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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

"""Shared fixtures and snapshots for electron channel regression tests."""

from pathlib import Path

import pytest

from rmgpy import settings
from rmgpy.data.base import ForbiddenStructures
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.rmg import RMGDatabase
import rmgpy.data.rmg as rmg_data
from rmgpy.rmg.main import RMG
import rmgpy.rmg.input as rmg_input
from rmgpy.thermo import ThermoData


LIBRARIES = ('ElectronDetachment', 'ElectronAttachment')


ORDERS = (LIBRARIES, LIBRARIES[::-1])


@pytest.fixture
def database(request):
    previous_database = rmg_data.database
    previous_rmg = rmg_input.rmg
    rmg_data.database = None
    db = RMGDatabase()
    db.kinetics = KineticsDatabase()
    db.kinetics.load_libraries(
        str(Path(settings['test_data.directory']) / 'electron_channels'),
        libraries=request.node.callspec.params['order'],
    )
    db.forbidden_structures = ForbiddenStructures()
    rmg_input.rmg = RMG()
    try:
        yield db
    finally:
        rmg_data.database = previous_database
        rmg_input.rmg = previous_rmg


def admission_state(model):
    """Snapshot registration, pending objects, indices and network containers."""
    def freeze(value):
        if isinstance(value, dict):
            return id(value), tuple((freeze(k), freeze(v)) for k, v in value.items())
        if isinstance(value, (list, tuple)):
            return id(value), tuple(freeze(v) for v in value)
        if isinstance(value, (str, int, bool, type(None))):
            return value
        return id(value)

    return (model.species_counter, model.reaction_counter, model.network_count,
            tuple(freeze(getattr(model, name)) for name in (
                'species_dict', 'reaction_dict', 'index_species_dict', 'species_cache',
                'new_species_list', 'new_reaction_list', 'network_dict', 'network_list')),
            tuple(freeze(items) for part in (model.core, model.edge)
                  for items in (part.species, part.reactions)))


def controlled_thermo(species, solvent_name=None):
    species.thermo = ThermoData(Tdata=([300, 1000], 'K'), Cpdata=([30, 30], 'J/(mol*K)'),
                               H298=(0, 'kJ/mol'), S298=(100, 'J/(mol*K)'),
                               Cp0=(30, 'J/(mol*K)'), CpInf=(30, 'J/(mol*K)'))


def registry_contents(model):
    """Compare independently built registries including membership and indices."""
    def species(s):
        return s.index, s.label, s.molecule[0].to_adjacency_list()

    def reaction(r):
        return r.index, r.family, tuple(species(s) for s in r.reactants), tuple(species(s) for s in r.products)

    return (model.species_counter, model.reaction_counter, model.network_count,
            {formula: [species(s) for s in values] for formula, values in model.species_dict.items()},
            {index: species(s) for index, s in model.index_species_dict.items()},
            [None if s is None else species(s) for s in model.species_cache],
            {family: {first: {second: [reaction(r) for r in reactions] for second, reactions in groups.items()}
                      for first, groups in values.items()} for family, values in model.reaction_dict.items()},
            [species(s) for s in model.new_species_list], [reaction(r) for r in model.new_reaction_list],
            [species(s) for s in model.edge.species], [reaction(r) for r in model.edge.reactions],
            [species(s) for s in model.core.species], [reaction(r) for r in model.core.reactions],
            model.network_dict, model.network_list)


_COLLIDER_BOUNDARIES = (
    'process', 'make_reaction', 'register', 'make_pdep', 'model_network',
    'path', 'merge', 'restore', 'configurations', 'update', 'explore',
)
