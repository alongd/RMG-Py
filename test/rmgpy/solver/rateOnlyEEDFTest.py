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

"""Exercise chemical-only columns through the real table/provider/reactor seam."""

import copy
from types import SimpleNamespace

import numpy as np
import pytest

import plasmaEEDFCoreTest as core_fixture
import rmgpy.data.rmg as data_module
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary
from rmgpy.kinetics import EEDFChannel
from rmgpy.solver.plasma import PlasmaReactor, PlasmaStateError
from rmgpy.species import Species
from rmgpy.tools.eedf.schema import _content_hash, _file_hash
from rmgpy.tools.eedf.rate_only import PROCESSES, SOURCE_SHA256, energy_schema


def fixture(tmp_path, monkeypatch, mutation=None, declared=None, duplicate=False):
    core, ionization = core_fixture.species_and_reaction()
    n2 = Species(label='N2').from_smiles('N#N')
    n = Species(label='N').from_adjacency_list('multiplicity 4\n1 N u3 p1 c0')
    d = Species(label='N(2D)').from_adjacency_list('multiplicity 2\nelectronicstate 2D\n1 N u1 p2 c0')
    for sp in (n2, n, d):
        sp.thermo = core_fixture.thermo()
    core.extend((n2, n, d))
    thermo = ThermoDatabase()
    library = ThermoLibrary(label='RateOnlyFixture')
    library.thermo_convention = 'ion'
    for i, sp in enumerate(core[1:]):
        sp.thermo.Cp0 = (20., 'J/(mol*K)')
        sp.thermo.CpInf = (20., 'J/(mol*K)')
        library.entries[i] = Entry(index=i, label=sp.label, item=sp.molecule[0].copy(deep=True),
                                  data=copy.deepcopy(sp.thermo))
    thermo.libraries = {library.label: library}
    thermo.library_order = [library.label]
    monkeypatch.setattr(data_module, 'database', SimpleNamespace(thermo=thermo, solvation=None))
    reactions = [ionization]
    chemicals = []
    for index, process, products in zip((18, 19), PROCESSES, ([n, n], [n, d])):
        marker = EEDFChannel(process=process, collision_set='test', role='rate_only')
        reaction = LibraryReaction(reactants=[core[0], n2], products=[core[0], *products],
                                   kinetics=marker, reversible=False, library='PlasmaAtomicNitrogen',
                                   entry=Entry(index=index))
        reactions.append(reaction)
        chemicals.append(dict(description=process, role='rate_only', source_sha256=SOURCE_SHA256,
                              classification='A', kind='excitation', threshold_eV=10.,
                              cross_section={'energy_eV': [10., 20.], 'sigma_m2': [0., 1e-20]},
                              reaction={'library': reaction.library, 'index': index, 'repr': repr(marker)}))
    write = core_fixture._write_artifact

    def append_chemical_columns(root, manifest, branches, held_out):
        manifest['channel_map'].extend([
            dict(description='e + N2(X,v=0) -> e + N2(sing-sup), Excitation',
                 classification='B', kind='excitation', threshold_eV=13., reaction=None),
            *chemicals])
        inputs = manifest['row_inputs']
        inputs['reactions'].extend(c['reaction'] for c in chemicals)
        inputs['rate_only'] = {'source_sha256': [SOURCE_SHA256],
                               'energy_schema': energy_schema(manifest['channel_map'])}
        for row in branches['branch_0'].flat:
            for name, values in (
                    ('k_ine', [0., 1e-16, 2e-16]), ('k_sup', [0., 0., 0.]),
                    ('channel_power', [0., 0., 0.]), ('attachment_energy_eV', [0., 0., 0.]),
                    ('target_fractions', [.05, .05, .05]), ('product_fractions', [0., 0., 0.]),
                    ('rate_floors', [1e-40, 1e-40, 1e-40])):
                row[name] = np.r_[row[name], values]
            row['below_floor'] = np.zeros((5, 2), dtype=bool)
        if mutation:
            mutation(manifest, branches)
        manifest['fingerprint'] = _content_hash(inputs)
        return write(root, manifest, branches, held_out)

    monkeypatch.setattr(core_fixture, '_write_artifact', append_chemical_columns)
    path = core_fixture.artifact(tmp_path, ionization.kinetics, composition=False, gases=['Ar', 'N2'])
    initial = dict(zip(core, [1e-6, .949998, 1e-6, 0., .05, 0., 0.]))
    energy = dict(absorbed_power=(.5, 'W'), chamber_volume=(core_fixture.VOLUME, 'm^3'),
                  sheath='floating_wall')
    if declared is not None:
        energy['electron_energies'] = declared
    reactor = PlasmaReactor(
        (core_fixture.TGAS, 'K'), (core_fixture.PRESSURE, 'Pa'), initial,
        (core_fixture.TGAS, 'K'), termination=[],
        diffusion_length=(core_fixture.LAMBDA, 'm'),
        ion_reduced_mobility=(core_fixture.MU0, 'm^2/(V*s)'),
        wall_single_bath_approximation=True,
        wall_neutralization_products={'Ar+': 'Ar'},
        electron_energy_balance=energy,
        electron_kinetics={'provider': 'loki-table', 'table': (str(path), _file_hash(path/'table.h5')),
                           'branch': 'branch_0', 'initial_reduced_field': (np.e, 'Td'),
                           'empirical_laws': {}}, thermo_source_assertions={'Ar+': 'ion'})
    if duplicate:
        reactions.append(copy.copy(reactions[-1]))
    reactor.initialize_model(core, reactions, [], [])
    return reactor, core, reactions, path


def test_rate_only_runtime_produces_atoms_without_energy_loss(tmp_path, monkeypatch):
    reactor, core, reactions, path = fixture(tmp_path, monkeypatch)
    residual = reactor.residual(0., reactor.y.copy(), np.zeros_like(reactor.y))[0]
    assert residual[reactor.species_index[core[-2]]] > 0.
    assert residual[reactor.species_index[core[-1]]] > 0.
    assert reactor.core_reaction_rates[1] > 0.
    assert reactor.core_reaction_rates[2] > 0.
    terms = reactor.electron_energy_terms
    assert np.array_equal(terms['Q_inelastic_by_reaction'][1:], [0., 0.])
    assert np.array_equal(terms['Q_inelastic_channels'][-2:], [0., 0.])
    assert np.array_equal(terms['Q_superelastic_channels'][-2:], [0., 0.])


def test_rate_only_energy_binding_refuses(tmp_path, monkeypatch):
    with pytest.raises(PlasmaStateError, match='rate-only owner cannot bind'):
        fixture(tmp_path, monkeypatch, declared={'PlasmaAtomicNitrogen:18': (10., 'eV')})


def test_duplicate_chemical_binding_refuses(tmp_path, monkeypatch):
    with pytest.raises(PlasmaStateError, match='more than one EEDF reaction'):
        fixture(tmp_path, monkeypatch, duplicate=True)


@pytest.mark.parametrize('mutation', ['metadata', 'energy', 'inverse', 'role'])
def test_bad_chemical_table_refuses(tmp_path, monkeypatch, mutation):
    def change(manifest, branches):
        if mutation == 'metadata':
            del manifest['row_inputs']['rate_only']
        elif mutation == 'role':
            manifest['channel_map'][-1]['role'] = 'collision'
        else:
            name = 'channel_power' if mutation == 'energy' else 'k_sup'
            branches['branch_0'].flat[0][name][-1] = 1e-20
    with pytest.raises((ValueError, PlasmaStateError), match='rate-only|role|energy schema'):
        fixture(tmp_path, monkeypatch, mutation=change)
