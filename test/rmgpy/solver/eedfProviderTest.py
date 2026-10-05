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

import copy
import gc
import threading

import numpy as np
import pytest

from rmgpy.solver.eedf import (
    AmbiguousBranch,
    EEDFError,
    EnvelopeBreach,
    IllConditionedCoordinate,
    OutOfDomain,
)
from rmgpy.solver.eedf_provider import (
    DEVELOPMENT_UNQUALIFIED_STATUS,
    EEDFProvider,
    development_unqualified_table_route,
)
from rmgpy.tools.eedf.artifact import write_artifact
from rmgpy.tools.eedf.schema import (
    FingerprintMismatch,
    content_hash,
    file_hash,
    interpolant_identity,
)


def fixture_provider_table(tmp_path, composition=False, branches=1, decreasing=False,
                           accepted=True):
    axes = {'u': [0., 1., 2., 3.]}
    if composition:
        axes['x'] = [0., 1.]
    shape = tuple(len(axis) for axis in axes.values())

    def branch_values():
        values = np.empty(shape, dtype=object)
        for index in np.ndindex(shape):
            u_index = index[0]
            x_index = index[1] if composition else 0
            scale = (1., 2., 4., 8.)[u_index] * (1. + x_index)
            composition_values = {'x': axes['x'][x_index]} if composition else {}
            values[index] = {
                'u': axes['u'][u_index],
                'EN_Td': np.exp(axes['u'][u_index]),
                'swarm': {'mean_energy_eV': ((4. - u_index) if decreasing else (1. + u_index)) + x_index,
                          'mobility_N': 1e24 * scale,
                          'diffusion_N': 2e24 * scale},
                'k_ine': np.array([scale, 10. * scale]) * 1e-15,
                'k_sup': np.array([1.5 * scale, 15. * scale]) * 1e-15,
                'channel_power': np.array([2. * scale, 20. * scale]) * 1e-15,
                'attachment_energy_eV': np.zeros(2),
                'target_fractions': np.ones(2),
                'product_fractions': np.zeros(2),
                'rate_floors': np.full(2, 1e-40),
                'below_floor': np.zeros((2, 2), dtype=bool),
                'power_groups': {'field': scale},
                'f0': np.array([1.]),
                'energy_eV': np.array([.5]),
                'energy_edges_eV': np.array([0., 1.]),
                'gas_fractions': {'Ar': 1.},
                'state_populations': {'Ar(1S0)': 1.},
                'composition': composition_values,
                'converged': True,
                'iteration_count': -1,
                'setup': 'provider fixture',
            }
        return values

    identities = [
        {'library': 'fixture', 'index': 1, 'repr': 'first'},
        {'library': 'fixture', 'index': 2, 'repr': 'second'},
    ]
    model = {
        'axes': axes,
        'Tg_K': 298.15,
        'P_Pa': 666.61,
        'arm': {'gas': 'Ar'},
        'cross_sections': {'Ar': 'sha'},
        'channel_map_sha256': 'map',
        'reactions': identities,
        'solver_options': {'growth': 'temporal'},
        'interpolant': interpolant_identity(),
        'state_populations': {'Ar(1S0)': 1.},
    }
    manifest = {
        'accepted': accepted,
        'held_out_verdicts': [{'point': {'u': .5}, 'branch_id': 'branch_0',
                               'passed': accepted, 'checks': [{'passed': accepted}]}],
        'branch_certification': {
            f'branch_{index}': 'uncertified: unseeded scans'
            for index in range(branches)
        },
        'schema_version': 1,
        'row_inputs': model,
        'fingerprint': content_hash(model),
        'axes': axes,
        'envelopes': {
            'Tg_K': {'min': 298.15, 'max': 298.15, 'reference': 298.15},
            'metastable': {'min': 0., 'max': 1e-5, 'reference': 0.},
        },
        'tolerances': {
            'fingerprint_rtol': 1e-9,
            'u_condition': {'min_abs_denergy_du': 1e-6,
                            'min_scaled_denergy_du': 1e-6},
        },
        'floors': {'rate_absolute': 1e-40},
        'repositories': {'rmgpy': 'old', 'database': 'old'},
        'channel_map': [
            {'description': 'first', 'kind': 'excitation', 'classification': 'A',
             'reaction': identities[0]},
            {'description': 'second', 'kind': 'excitation', 'classification': 'A',
             'reaction': identities[1]},
        ],
        'energy_eV': [.5],
        'energy_edges_eV': [0., 1.],
    }
    path = write_artifact(
        tmp_path,
        manifest,
        {f'branch_{index}': branch_values() for index in range(branches)},
        [],
    )
    return path, model


def make_provider(path, model, **overrides):
    arguments = {
        'artifact_sha256': file_hash(path / 'table.h5'),
        'reaction_map': {4: (1, 'ine'), 7: (0, 'sup')},
        'branch': 'branch_0',
        'empirical_laws': {'Ar2+ recombination': {'evaluate_at': 'Te_eff'}},
    }
    arguments.update(overrides)
    return EEDFProvider(path, model, **arguments)


def test_pure_argon_and_composition_axis_use_the_same_provider_type(tmp_path):
    pure_path, pure_model = fixture_provider_table(tmp_path / 'pure')
    arm_path, arm_model = fixture_provider_table(tmp_path / 'arm', composition=True)

    pure = make_provider(pure_path, pure_model)
    arm = make_provider(arm_path, arm_model)

    assert type(pure) is EEDFProvider
    assert type(arm) is EEDFProvider
    assert pure.empirical_laws['Ar2+ recombination']['evaluate_at'] == 'Te_eff'
    assert arm.manifest['fingerprint'] == content_hash(arm_model)


def test_unqualified_table_is_development_context_only_and_carries_blockers(tmp_path):
    qualified_path, qualified_model = fixture_provider_table(
        tmp_path / 'qualified', accepted=True)
    with development_unqualified_table_route():
        with pytest.raises(EEDFError, match='already qualified'):
            make_provider(qualified_path, qualified_model)

    path, model = fixture_provider_table(tmp_path / 'unqualified', accepted=False)
    with pytest.raises(FingerprintMismatch, match='held-out qualification'):
        make_provider(path, model)

    with development_unqualified_table_route():
        provider = make_provider(path, model)

    notice = provider.development_notice
    assert notice['scientific_status'] == DEVELOPMENT_UNQUALIFIED_STATUS
    assert notice['artifact_blockers'] == [
        'artifact manifest is not accepted',
        '1 of 1 held-out verdicts failed',
        'branch_0: uncertified: unseeded scans',
    ]
    assert provider.manifest['scientific_status'] == DEVELOPMENT_UNQUALIFIED_STATUS
    with pytest.raises(EEDFError, match='refuses qualification'):
        provider.require_qualified('qualification')
    with pytest.raises(EEDFError, match='refuses export'):
        provider.require_qualified('export')

    # Leaving the context restores the production loader immediately.
    with pytest.raises(FingerprintMismatch, match='held-out qualification'):
        make_provider(path, model)


def test_complete_composition_identity_controls_the_one_row_cache(tmp_path):
    path, model = fixture_provider_table(tmp_path, composition=True)
    provider = make_provider(path, model)

    first = provider.row(1., {'x': 0.})
    first_identity = provider.cache_identity
    assert provider.row(1., {'x': 0.}) is first

    changed = provider.row(1., {'x': 1.})
    assert changed is not first
    assert provider.cache_identity != first_identity
    assert provider.reaction_rate(4, changed) == pytest.approx(2. * first.k_ine[1])
    assert provider.mean_energy_eV(changed) == pytest.approx(3.)
    assert provider.transport(changed)['mobility_N'] == pytest.approx(4e24)


def test_reaction_mapping_selects_ine_and_sup_and_refuses_bad_entries(tmp_path):
    path, model = fixture_provider_table(tmp_path)
    provider = make_provider(path, model)
    row = provider.row(1., {})

    assert provider.reaction_rate(4, row) == pytest.approx(row.k_ine[1])
    assert provider.reaction_rate(7, row) == pytest.approx(row.k_sup[0])
    with pytest.raises(KeyError, match='unmapped reaction index'):
        provider.reaction_rate(99, row)
    with pytest.raises(ValueError, match='channel column'):
        make_provider(path, model, reaction_map={0: (2, 'ine')})
    with pytest.raises(ValueError, match='channel side'):
        make_provider(path, model, reaction_map={0: (0, 'forward')})
    with pytest.raises(ValueError, match='reaction index'):
        make_provider(path, model, reaction_map={-1: (0, 'ine')})


def test_configuration_and_issued_rows_are_defensive(tmp_path):
    path, model = fixture_provider_table(tmp_path / 'one')
    other_path, other_model = fixture_provider_table(tmp_path / 'other')
    provider = make_provider(path, model)
    other = make_provider(other_path, other_model)
    row = provider.row(1., {})

    with pytest.raises(TypeError):
        provider.reaction_map[4] = (0, 'ine')
    manifest = provider.manifest
    manifest['fingerprint'] = 'mutated copy'
    assert provider.manifest['fingerprint'] != 'mutated copy'
    with pytest.raises(ValueError):
        row.k_ine[0] = 0.
    with pytest.raises(ValueError):
        row.k_ine.setflags(write=True)
    with pytest.raises(FingerprintMismatch, match='not issued'):
        other.reaction_rate(4, row)
    for u in np.linspace(0., 3., 50):
        provider.row(float(u), {})
    gc.collect()
    assert len(provider._issued_rows) <= 2


def test_cache_identity_and_row_are_published_atomically(tmp_path, monkeypatch):
    path, model = fixture_provider_table(tmp_path)
    provider = make_provider(path, model)
    first = provider.row(1., {})
    entered = threading.Event()
    release = threading.Event()
    original = provider._table.row

    def paused_row(u, composition, branch):
        if u == 2.:
            entered.set()
            assert release.wait(5.)
        return original(u, composition, branch)

    monkeypatch.setattr(provider._table, 'row', paused_row)
    writer = threading.Thread(target=lambda: provider.row(2., {}))
    writer.start()
    assert entered.wait(5.)

    reader_result = []
    reader = threading.Thread(target=lambda: reader_result.append(provider.row(1., {})))
    reader.start()
    assert reader.is_alive()
    release.set()
    writer.join(5.)
    reader.join(5.)

    assert len(reader_result) == 1
    assert reader_result[0].u == first.u
    assert provider.cache_identity == (1., ())
    assert provider.cached_row is reader_result[0]


def test_fingerprint_branch_and_named_domain_errors_are_preserved(tmp_path):
    path, model = fixture_provider_table(tmp_path, composition=True, branches=2)
    changed_model = copy.deepcopy(model)
    changed_model['Tg_K'] = 300.
    with pytest.raises(FingerprintMismatch, match='Tg_K'):
        make_provider(path, changed_model)
    with pytest.raises(AmbiguousBranch, match='declared branch missing'):
        make_provider(path, model, branch='missing')

    provider = make_provider(path, model)
    with pytest.raises(OutOfDomain, match='u'):
        provider.domain_check(-.1, {'x': .5}, y=np.zeros(1), context={'accepted': True})
    with pytest.raises(OutOfDomain, match='x'):
        provider.domain_check(1., {'x': 2.})
    with pytest.raises(EnvelopeBreach, match='metastable'):
        provider.domain_check(1., {'x': .5, 'metastable': 2e-5})


def test_provider_refuses_a_monotone_decreasing_energy_branch(tmp_path):
    path, model = fixture_provider_table(tmp_path, decreasing=True)
    with pytest.raises(IllConditionedCoordinate, match='increase strictly'):
        make_provider(path, model)
