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

"""Synthetic generated ladder -> Tv table -> native plasma solver integration."""
import copy
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import rmgpy.constants as constants
import rmgpy.data.rmg as rmg_data_module
from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.kinetics import EEDFChannel
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.solver.eedf import EEDFTable, OutOfDomain
from rmgpy.species import Species
from rmgpy.tools.eedf.artifact import _write_artifact
from rmgpy.tools.eedf.schema import _content_hash, _file_hash, SpecError
from rmgpy.tools.eedf.vibrational import (
    boltzmann_populations, temperature, D6RegenerationRequired, D6CeilingExceeded,
)
from plasmaEEDFCoreTest import artifact, species_and_reaction, thermo, TGAS, PRESSURE, VOLUME, LAMBDA, MU0

_spec = importlib.util.spec_from_file_location(
    'vibrational_coupling_fixture', Path(__file__).parents[1] / 'tools' / 'eedf' / 'vibrationalTest.py')
_fixture = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_fixture)
source, small_ladder = _fixture.source, _fixture.small_ladder


@pytest.fixture(autouse=True)
def isolated_database(monkeypatch):
    monkeypatch.setattr(rmg_data_module, 'database', None)


def molecular_model(tmp_path, small_ladder, *, policy=None, wrong_endpoint=False):
    emitted, ladder = copy.deepcopy(small_ladder)
    if policy:
        ladder['D6'].update(policy)
    core, ionization = species_and_reaction()
    vib = [Species(label='renamed_' + str(s['v'])).from_adjacency_list(s['adjacency_list'])
           for s in emitted['states']]
    database, library = ThermoDatabase(), ThermoLibrary(label='SyntheticVibrational')
    library.load(emitted['states'][0]['thermo']['path'], database.local_context, {})
    argon_library = ThermoLibrary(label='SyntheticArgon', thermo_convention='ion')
    for index, species in enumerate(core[1:]):
        species.thermo.Cp0 = (20., 'J/(mol*K)')
        species.thermo.CpInf = (20., 'J/(mol*K)')
        argon_library.entries[index] = Entry(index=index, label=species.label,
                                              item=species.molecule[0].copy(deep=True),
                                              data=copy.deepcopy(species.thermo))
    database.libraries = {library.label: library, argon_library.label: argon_library}
    database.library_order = [library.label, argon_library.label]
    rmg_data_module.database = SimpleNamespace(thermo=database, solvation=None)
    for species in vib:
        species.thermo = database.get_thermo_data(species)
    core.extend(vib)
    reactions = [ionization]
    base = artifact(tmp_path / 'base', ionization.kinetics, composition=False)
    manifest = json.loads((base / 'manifest.json').read_text())
    base_table = EEDFTable.load(base, manifest['row_inputs'], artifact_sha256=_file_hash(base / 'table.h5'))
    axes = {'u': [0., 1., 2., 3.], 'Tv_N2': [400., 800., 1200.]}
    channels = manifest['channel_map']
    for pair in ladder['channels']:
        i, j = pair['initial_v'], pair['source_final_v']
        identity = None
        K = ladder['K']
        if (i < j <= K) or (i < K and j > K):
            marker = EEDFChannel(process=pair['description'], collision_set='synthetic-vibrational', side='ine')
            reaction = LibraryReaction(reactants=[core[0], vib[i]], products=[core[0], vib[min(j, K)]],
                                       kinetics=marker, reversible=False, library='SyntheticVibrational',
                                       entry=Entry(index=len(reactions) + 1))
            reactions.append(reaction)
            identity = {'library': reaction.library, 'index': reaction.entry.index, 'repr': repr(marker)}
            if j <= K:
                reverse = LibraryReaction(reactants=[core[0], vib[j]], products=[core[0], vib[i]],
                                          kinetics=EEDFChannel(process=pair['description'],
                                                               collision_set='synthetic-vibrational', side='sup'),
                                          reversible=False, library='SyntheticVibrational',
                                          entry=Entry(index=len(reactions) + 1))
                reactions.append(reverse)
        channels.append({'description': pair['description'], 'kind': 'vibrational',
                         'classification': 'A' if identity else 'B', 'reaction': identity,
                         'threshold_eV': pair['threshold_eV'],
                         'target_fraction': 'state:N2(X,v=' + str(i) + ')',
                         'product_fraction': 'state:N2(X,v=' + str(j) + ')'})
    for index, reaction in enumerate(reactions, 1):
        reaction.index = index
    if wrong_endpoint:
        reactions[1].reactants = [core[0], vib[1]]
    model = manifest['row_inputs']
    model.update(axes=axes, vibrational_ladders=[ladder],
                 arm={'gases': ['Ar', 'N2'], 'additive': 'N2', 'feed_partial_pressure_Pa': .2 * PRESSURE},
                 reactions=[c['reaction'] for c in channels if c['reaction'] is not None])
    model['engine_state_map'] = [
        {'formula': 'Ar', 'electronic_state': '', 'vibrational_level': -1,
         'multiplicity': s.molecule[0].multiplicity, 'loki_state': state, 'statistical_weight': 1.}
        for s, state in ((core[1], 'Ar(1S0)'), (core[3], 'Ar(3P2)'))]
    manifest.update(axes=axes, row_inputs=model, fingerprint=_content_hash(model))
    rows = np.empty((4, 3), object)
    for index in np.ndindex(rows.shape):
        u, tv = axes['u'][index[0]], axes['Tv_N2'][index[1]]
        row = base_table.row(u, {})._as_dict()
        row['composition'] = {'Tv_N2': tv}
        row['swarm']['mean_energy_eV'] += .001 * tv
        population = boltzmann_populations(ladder, tv)
        ine = np.array([1e-18 * (1 + u) if c['threshold_eV'] else 0. for c in ladder['channels']])
        sup = .3 * ine
        target = np.array([.2 * population[c['initial_v']] for c in ladder['channels']])
        product = np.array([.2 * population[c['source_final_v']] for c in ladder['channels']])
        threshold = np.array([c['threshold_eV'] for c in ladder['channels']])
        for name, extra in (
                ('k_ine', ine), ('k_sup', sup), ('target_fractions', target),
                ('product_fractions', product), ('channel_power', threshold * (target * ine - product * sup)),
                ('rate_floors', np.full(len(ine), 1e-40)), ('attachment_energy_eV', np.zeros(len(ine)))):
            row[name] = np.concatenate([row[name], extra])
        row['below_floor'] = np.stack([row['k_ine'] < row['rate_floors'],
                                      row['k_sup'] < row['rate_floors']], axis=-1)
        row['gas_fractions'] = {'Ar': .8, 'N2': .2}
        row['state_populations'] = {'Ar(1S0)': 1., 'Ar(3P2)': 0., 'N2(X)': 1.,
                                    **{s['loki_state']: population[s['v']] for s in ladder['levels']}}
        rows[index] = row
    path = _write_artifact(tmp_path / 'molecular', manifest, {'branch_0': rows}, [])
    # A synthetic table has no executable LoKI specification. Supply only the
    # failure-record destination; this cannot qualify a direct solver run.
    (path / 'generation_spec.json').write_text(json.dumps({'scratch_root': str(tmp_path / 'qualification')}))
    population = boltzmann_populations(ladder, 800.)
    initial = {core[0]: 1e-6, core[1]: .799998, core[2]: 1e-6, core[3]: 0.,
               **{species: .2 * population[v] for v, species in enumerate(vib)}}
    reactor = PlasmaReactor(
        (TGAS, 'K'), (PRESSURE, 'Pa'), initial, (TGAS, 'K'), termination=[],
        diffusion_length=(LAMBDA, 'm'), ion_reduced_mobility=(MU0, 'm^2/(V*s)'),
        wall_single_bath_approximation=True, wall_neutralization_products={'Ar+': 'Ar'},
        electron_energy_balance={'absorbed_power': (.5, 'W'), 'chamber_volume': (VOLUME, 'm^3'),
                                 'sheath': 'floating_wall'},
        electron_kinetics={'provider': 'loki-table', 'table': (str(path), _file_hash(path / 'table.h5')),
                          'branch': 'branch_0', 'initial_reduced_field': (np.e, 'Td'), 'empirical_laws': {}},
        thermo_source_assertions={'Ar+': 'ion'})
    reactor.initialize_model(core, reactions, [], [])
    return reactor, core, reactions, ladder


def test_generated_ladder_integrates_and_binds_both_sides(tmp_path, small_ladder):
    reactor, core, reactions, ladder = molecular_model(tmp_path, small_ladder)
    assert reactor._eedf_coordinates(reactor.y)['Tv_N2'] == pytest.approx(800.)
    for index, reaction in enumerate(reactions[1:], 1):
        column, side = reactor.eedf_reaction_map[index]
        expected = reactor.eedf_row.k_ine[column] if side == 'ine' else reactor.eedf_row.k_sup[column]
        assert reactor.kf[index] == pytest.approx(constants.Na * expected)
    initial = reactor.y.copy()
    reactor.advance(1e-10)
    assert reactor.t >= 1e-10
    assert np.isfinite(reactor.y).all()
    assert np.any(reactor.y != initial)
    terminal = reactor._eedf_qualification_state(reactor.y, reactor.t)
    assert terminal['vibrational_D6']['Tv_N2']['action'] == 'accept'
    populations = terminal['state_populations']
    assert populations['N2(X)'] == pytest.approx(1.)
    assert populations['N2(X,v=0)'] + populations['N2(X,v=1)'] == pytest.approx(1.)
    assert populations['N2(X,v=2)'] == 0.


def test_wrong_reaction_population_refuses(tmp_path, small_ladder):
    with pytest.raises(SpecError, match='exact level population'):
        molecular_model(tmp_path, small_ladder, wrong_endpoint=True)


def test_tv_axis_domain_chain_rule_and_energy_budget(tmp_path, small_ladder):
    reactor, core, reactions, ladder = molecular_model(tmp_path, small_ladder)
    provider = reactor.eedf_provider
    low, high = provider.row(1., {'Tv_N2': 400.}), provider.row(1., {'Tv_N2': 1200.})
    middle = provider.row(1., {'Tv_N2': 800.})
    assert middle.swarm['mean_energy_eV'] == pytest.approx((low.swarm['mean_energy_eV'] + high.swarm['mean_energy_eV']) / 2)
    with pytest.raises(OutOfDomain, match='Tv_N2'):
        provider.row(1., {'Tv_N2': 1201.})
    y = reactor.y.copy()
    rates = np.zeros_like(y)
    indices = reactor.eedf_axis_species['Tv_N2']['indices']
    rates[indices[0]], rates[indices[1]] = -1e-9, 1e-9
    actual = reactor._eedf_composition_energy_derivative(y, rates)
    step = 1e-4
    def mean(state):
        return provider.mean_energy_eV(provider.row(state[reactor.te_index], reactor._eedf_coordinates(state)))
    fd = (mean(y + step * rates) - mean(y - step * rates)) / (2 * step)
    fd *= y[reactor.electron_index] * constants.Na * constants.e
    assert actual == pytest.approx(fd, rel=1e-6)
    residual, status = reactor.residual(reactor.t, y, np.zeros_like(y))
    assert status == 1
    terms = reactor.electron_energy_terms
    ledger = terms['vibrational_energy']['Tv_N2']
    assert ledger['top_leak_power'] > 0.
    assert ledger['forward_power'] == pytest.approx(ledger['retained_power'] + ledger['top_leak_power'] + ledger['source_residual_power'])
    gross = sum(terms['Q_inelastic_channels'][2:])
    assert ledger['forward_power'] == pytest.approx(gross)
    # The full threshold is already in Q_inelastic, leak is diagnostic only.
    assert terms['Q_inelastic'] == pytest.approx(sum(terms['Q_inelastic_channels']) - sum(terms['Q_superelastic_channels']))
    electron_energy_rate = (
        residual[reactor.te_index] * y[reactor.electron_index] * constants.Na * constants.e * reactor.eedf_row.denergy_du
        + terms['mean_energy_eV'] * constants.Na * constants.e * terms['dNe_dt']
        + terms['composition_energy_derivative'])
    losses = sum(terms[key] for key in ('Q_inelastic', 'Q_elastic', 'Q_wall_electron', 'Q_wall_ion', 'Q_flow'))
    assert electron_energy_rate == pytest.approx(terms['P_abs'] - losses, rel=1e-12)


def test_terminal_runtime_d6_refuses_before_qualification(tmp_path, small_ladder):
    reactor, _, _, _ = molecular_model(tmp_path, small_ladder, policy={'top_fraction_limit': 0.})
    with pytest.raises(D6RegenerationRequired, match='D6_REGENERATE Tv_N2 K=6'):
        reactor._eedf_qualification_state(reactor.y, reactor.t)
    # This callback is invoked by the provider's terminal qualification transaction.
    with pytest.raises(Exception, match='D6_REGENERATE'):
        reactor.eedf_provider.qualify(reactor.y, reactor.t)
    assert reactor.eedf_provider._qualification_status == 'FAILED'


def test_runtime_ceiling_refuses_by_name(tmp_path, source):
    emitted = _fixture._generator.decoded(source, K=30)
    bindings = [{'source_sha256': p['source']['sha256'], 'source_block': p['source_block'],
                 'description': f"e + N2(X,v={p['initial_v']}) <-> e + N2(X,v={p['final_v']}), Vibrational"}
                for p in emitted['transition_audit']]
    ladder = _fixture.ladder_coupling(emitted, bindings, axis='Tv_N2', parent='N2(X)')
    reactor, _, _, _ = molecular_model(tmp_path, (emitted, ladder), policy={
        'top_fraction_limit': 0., 'top_leak_share_limit': 1.})
    with pytest.raises(D6CeilingExceeded, match='D6_REFUSE_K30'):
        reactor._eedf_qualification_state(reactor.y, reactor.t)


def test_vibrational_binding_refuses_thermal_reversal(tmp_path, small_ladder):
    from rmgpy.tools.eedf.vibrational import validate_bound_channels
    reactor, core, reactions, ladder = molecular_model(tmp_path, small_ladder)
    reactions[1].reversible = True
    with pytest.raises(SpecError, match='must be irreversible'):
        validate_bound_channels(ladder, reactor.eedf_provider.runtime_metadata['channel_map'],
                                reactions, reactor.eedf_reaction_map,
                                reactor.eedf_axis_species['Tv_N2']['indices'], core)
