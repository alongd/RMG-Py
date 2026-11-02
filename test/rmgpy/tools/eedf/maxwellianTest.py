"""Independent real-section Maxwellian integration and native-source integrity."""
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.constants import elementary_charge, electron_mass
from scipy.integrate import quad, simpson

from rmgpy.tools.eedf.artifact import _write_artifact
from rmgpy.tools.eedf.maxwellian import (
    _source_snapshot, _maxwellian_distribution, _maxwellian_row,
)
from rmgpy.tools.eedf.moments import _transport_moments
from rmgpy.tools.eedf.schema import _content_hash, EEDFError, FingerprintMismatch


FIXTURE = Path(__file__).parent / 'fixtures' / 'maxwellian-ar.json'


def source_data(temperatures=(1., 1.5, 2.), max_energy=30., cells=3000):
    channels = json.loads(FIXTURE.read_text())['channels']
    edges = np.linspace(0., max_energy, cells + 1)
    energy = (edges[:-1] + edges[1:]) / 2
    axes = {'u': [float(i) for i in range(len(temperatures))]}
    inputs = {'axes': axes, 'Tg_K': 298.15, 'P_Pa': 666.61184,
              'interpolant': {'version': 'historical', 'source_sha256': 'historical'},
              'solver_options': {'growthModelType': 'temporal',
                                 'ionizationOperatorType': 'usingSDCS',
                                 'includeEECollisions': False}}
    manifest = {'schema_version': 1, 'provider_kind': 'loki-table',
                'row_inputs': inputs, 'fingerprint': _content_hash(inputs),
                'axes': axes, 'energy_eV': energy.tolist(), 'energy_edges_eV': edges.tolist(),
                'channel_map': channels, 'branch_certification': {'branch_0': 'uncertified'},
                'tolerances': {'H2': {'rtol': .001, 'atol': 0}},
                'floors': {'rate_absolute': 1e-40}, 'envelopes': {}, 'accepted': False}
    rows = []
    for i, temperature in enumerate(temperatures):
        rows.append({'u': float(i), 'EN_Td': np.exp(float(i)), 'composition': {},
                     'energy_eV': energy.copy(), 'energy_edges_eV': edges.copy(),
                     'swarm': {'mean_energy_eV': 1.5 * temperature, 'mobility_N': 7., 'diffusion_N': 8.},
                     'k_ine': np.full(3, 123.), 'k_sup': np.zeros(3),
                     'channel_power': np.full(3, 456.), 'attachment_energy_eV': np.zeros(3),
                     'target_fractions': np.ones(3), 'product_fractions': np.zeros(3),
                     'rate_floors': np.full(3, 1e-40), 'below_floor': np.zeros((3, 2), dtype=bool),
                     'gas_fractions': {'Ar': 1.}, 'state_populations': {'Ar(1S0)': 1.},
                     'f0': np.ones(cells), 'power_groups': {'field': 99.},
                     'converged': True, 'iteration_count': 42})
    return manifest, rows


def independent_rate(channel, temperature, energy, edges):
    """Adaptive piecewise-linear section integral, independent of cell kernels."""
    section = channel['cross_section']
    threshold = channel['threshold_eV']
    nodes = sorted({threshold, float(edges[-1]), *[
        float(x) for x in section['energy_eV'] if threshold < x < edges[-1]]})
    integral = sum(quad(
        lambda e: e * np.interp(e, section['energy_eV'], section['sigma_m2']) * np.exp(-e / temperature),
        left, right, epsabs=1e-35, epsrel=1e-10)[0]
        for left, right in zip(nodes[:-1], nodes[1:]))
    normalization = sum(2 / np.sqrt(np.pi) * temperature**-1.5 * np.exp(-e / temperature)
                        * np.sqrt(e) * width for e, width in zip(energy, np.diff(edges)))
    return integral * np.sqrt(2 * elementary_charge / electron_mass) * 2 / np.sqrt(np.pi) * temperature**-1.5 / normalization


@pytest.mark.parametrize('temperature', [1., 2., 2.8])
def test_real_elastic_excitation_ionization_integrals(temperature):
    manifest, rows = source_data((temperature, temperature + .01))
    row = _maxwellian_row(rows[0], manifest)
    for j, channel in enumerate(manifest['channel_map']):
        expected = independent_rate(channel, temperature, row['energy_eV'], row['energy_edges_eV'])
        # 0.01-eV cells; prior grid-halving probe measured <2.8e-5 error.
        assert expected > 0
        assert row['k_ine'][j] == pytest.approx(expected, rel=5e-5, abs=0)
    assert row['swarm']['maxwellian_temperature_eV'] == pytest.approx(temperature, rel=2e-15, abs=0)
    assert row['power_groups'] == {'collision_loss': float(sum(row['channel_power']))}
    assert 'field' not in row['power_groups']
    assert not row['converged'] and row['iteration_count'] == -1
    assert row['swarm']['mobility_N'] > 1e20
    assert row['swarm']['diffusion_N'] > 1e20
    assert not np.array_equal(row['f0'], rows[0]['f0'])
    assert not np.array_equal(row['channel_power'], rows[0]['channel_power'])
    assert rows[0]['power_groups']['field'] == 99.


@pytest.mark.parametrize('temperature', [1., 2., 2.8])
def test_real_reverse_channel_detailed_balance(temperature):
    manifest, rows = source_data((temperature, temperature + .01))
    channel = manifest['channel_map'][1]
    channel['description'] = channel['description'].replace(' -> ', ' <-> ')
    channel['statistical_weight_ratio'] = 1 / 5  # Ar ground J=0, metastable J=2.
    rows[0]['product_fractions'][1] = .1
    row = _maxwellian_row(rows[0], manifest)
    forward = independent_rate(channel, temperature, row['energy_eV'], row['energy_edges_eV'])
    expected = forward * (1 / 5) * np.exp(channel['threshold_eV'] / temperature)
    # LoKI's common collision kernel snaps the threshold down to a cell face.
    # At most one cell of Boltzmann-factor error, plus measured forward error.
    tolerance = np.expm1(.01 / temperature) + 5e-5
    assert row['k_sup'][1] > 0
    assert row['k_sup'][1] == pytest.approx(expected, rel=tolerance, abs=0)
    threshold_face = int(channel['threshold_eV'] / .01) * .01
    assert row['channel_power'][1] == pytest.approx(
        (row['k_ine'][1] - .1 * row['k_sup'][1]) * threshold_face, rel=1e-12)
    assert row['channel_power'][1] < 0


def test_grid_rejects_truncated_tail_and_invalid_cells():
    manifest, rows = source_data((3., 3.1))
    with pytest.raises(EEDFError, match='unresolved Maxwellian mean'):
        _maxwellian_row(rows[0], manifest)
    energy, edges = rows[0]['energy_eV'], rows[0]['energy_edges_eV']
    bad = edges.copy()
    bad[5] += .001
    with pytest.raises(EEDFError, match='uniform source'):
        _maxwellian_distribution(energy, bad, 2., .001)
    for temperature in (0., -1., np.nan, np.inf):
        with pytest.raises(EEDFError, match='temperature'):
            _maxwellian_distribution(energy, edges, temperature, .001)


@pytest.mark.parametrize('problem', ['mass', 'weights', 'section', 'floor', 'operator'])
def test_incomplete_physical_inputs_refuse(problem):
    manifest, rows = source_data()
    if problem == 'mass':
        del manifest['channel_map'][0]['mass_ratio']
    elif problem == 'weights':
        manifest['channel_map'][1]['description'] = manifest['channel_map'][1]['description'].replace(' -> ', ' <-> ')
    elif problem == 'section':
        del manifest['channel_map'][1]['cross_section']
    elif problem == 'floor':
        rows[0]['rate_floors'][0] = 0.
    else:
        manifest['row_inputs']['solver_options']['ionizationOperatorType'] = 'conservative'
    with pytest.raises(EEDFError):
        _maxwellian_row(rows[0], manifest)


def test_transport_matches_original_mobility_formula():
    energy = np.linspace(.005, 29.995, 3000)
    edges = np.linspace(0., 30., 3001)
    f0 = np.exp(-energy / 2.)
    momentum = (1 + edges / 30.) * 1e-20
    gamma = np.sqrt(2 * elementary_charge / electron_mass)
    for growth in (0., 1e-15):
        derivative = np.r_[0., np.diff(f0) / np.diff(energy), 0.]
        base = derivative * edges / momentum
        frequency = momentum.copy()
        frequency[1:] += growth / gamma / np.sqrt(edges[1:])
        corrected = derivative * edges / frequency
        expected = -gamma / 3 * simpson(corrected, x=edges)
        transport = _transport_moments(energy, edges, f0, momentum, growth)
        assert transport['mobility_N'] == expected
        assert transport['mobility_without_growth_N'] == -gamma / 3 * simpson(base, x=edges)
        # Cell-centred diffusion integral from LoKI's D0 definition.
        cell_frequency = (frequency[:-1] + frequency[1:]) / 2
        for name, power in [('diffusion_N', 1), ('energy_diffusion_N', 2)]:
            expected_diffusion = gamma * sum(energy**power * f0 / (3 * cell_frequency) * np.diff(edges))
            assert transport[name] == pytest.approx(expected_diffusion, rel=1e-14)


def native_artifact(tmp_path):
    manifest, rows = source_data()
    path = _write_artifact(tmp_path, manifest, {'branch_0': rows}, [])
    return path, manifest, rows


def test_historical_source_uses_native_bytes_not_runtime_loader(tmp_path, monkeypatch):
    from rmgpy.solver.eedf import EEDFTable
    path, expected, rows = native_artifact(tmp_path)
    monkeypatch.setattr(EEDFTable, 'load', lambda *a, **kw: pytest.fail('runtime loader called'))
    before = (path / 'table.h5').read_bytes()
    manifest, snapshot = _source_snapshot(path, artifact_sha256=path.name, branch_id='branch_0')
    assert manifest['row_inputs']['interpolant']['version'] == 'historical'
    assert len(snapshot) == len(rows) == 3
    assert np.array_equal(snapshot[1]['k_ine'], rows[1]['k_ine'])
    assert manifest['channel_map'] == expected['channel_map']
    assert (path / 'table.h5').read_bytes() == before


@pytest.mark.parametrize('field', ['channel_map', 'axes', 'tolerances', 'row_inputs'])
def test_source_manifest_tampering_refuses(tmp_path, field):
    path, _, _ = native_artifact(tmp_path)
    manifest_path = path / 'manifest.json'
    manifest = json.loads(manifest_path.read_bytes())
    if field == 'channel_map':
        manifest[field][0]['mass_ratio'] *= 2
    elif field == 'axes':
        manifest[field]['u'][0] = -.1
    elif field == 'tolerances':
        manifest[field]['H2']['rtol'] *= 2
    else:
        manifest[field]['Tg_K'] *= 2
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(FingerprintMismatch, match='source'):
        _source_snapshot(path, artifact_sha256=path.name, branch_id='branch_0')


def test_source_hash_and_parse_use_one_read(tmp_path, monkeypatch):
    path, _, _ = native_artifact(tmp_path)
    read = Path.read_bytes
    count = []
    def replace_after_read(file):
        content = read(file)
        if file == path / 'table.h5':
            count.append(file)
            file.write_bytes(b'replaced after snapshot')
        return content
    monkeypatch.setattr(Path, 'read_bytes', replace_after_read)
    manifest, rows = _source_snapshot(path, artifact_sha256=path.name, branch_id='branch_0')
    assert len(count) == 1
    assert rows.shape == (3,)
    assert manifest['artifact_sha256'] == path.name
    with pytest.raises(FingerprintMismatch, match='SHA-256'):
        _source_snapshot(path, artifact_sha256=path.name, branch_id='branch_0')
