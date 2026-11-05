"""Validity declarations for the eliminated argon upper levels."""
import json

import pytest

from rmgpy.solver.eedf import EEDFTable
from rmgpy.tools.eedf.schema import FingerprintMismatch, _file_hash
from plasmaEEDFCoreTest import artifact, species_and_reaction


@pytest.mark.parametrize('key,value', [
    ('excited_fraction_max', 1.e-5),
    ('argon_elimination_electron_density_max_m3', 1.e18),
])
@pytest.mark.parametrize('tamper', ['change', 'remove', 'add'])
def test_table_validity_cannot_be_changed_outside_fingerprint(tmp_path, key, value, tamper):
    _, reaction = species_and_reaction()
    bounds = {} if tamper == 'add' else {key: value}
    path = artifact(tmp_path, reaction.kinetics, argon_validity=bounds)
    manifest_path = path / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    if tamper == 'remove':
        del manifest[key]
    else:
        manifest[key] = value * 2.
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(FingerprintMismatch, match=key):
        EEDFTable.load(path, manifest['row_inputs'], artifact_sha256=_file_hash(path / 'table.h5'))


from rmgpy.exceptions import PlasmaStateError
from plasmaEEDFCoreTest import build_reactor, isolate_database


@pytest.mark.parametrize('state', [None, '3P2', '3P1', '3P0', '1P1'])
def test_reactor_refuses_argon_fraction_above_declared_bound(tmp_path, state):
    # This is relative to Ar ground, even when labels are aliases.
    with pytest.raises(PlasmaStateError, match=r'excited_fraction_max.*1e-05.*observed'):
        build_reactor(tmp_path, composition=False, excited_fraction=2.e-5,
                      excited_label='metastable-alias', excited_state=state,
                      argon_validity={'excited_fraction_max': 1.e-5})


@pytest.mark.parametrize('bounds,fraction', [
    ({'excited_fraction_max': 1.e-5,
      'argon_elimination_electron_density_max_m3': 1.e18}, 5.e-6),
    ({}, 2.e-5),
])
def test_reactor_accepts_in_bound_and_undeclared_fraction(tmp_path, bounds, fraction):
    reactor, _, _ = build_reactor(tmp_path, composition=False,
                                   excited_fraction=fraction, argon_validity=bounds)
    reactor.advance(1.e-12)
    assert reactor.t >= 1.e-12


def test_reactor_refuses_density_above_elimination_limit(tmp_path):
    with pytest.raises(PlasmaStateError, match=r'eliminated-upper-level.*1e\+18.*observed'):
        build_reactor(tmp_path, composition=False, electron_fraction=1.e-5,
                      argon_validity={'argon_elimination_electron_density_max_m3': 1.e18})


def test_reactor_rechecks_elimination_density_during_advance(tmp_path):
    reactor, _, _ = build_reactor(tmp_path, composition=False, electron_fraction=6.e-6,
                                   argon_validity={'argon_elimination_electron_density_max_m3': 1.e18})
    with pytest.raises(PlasmaStateError, match='eliminated-upper-level'):
        reactor.advance(1.e-8)
    assert reactor.t > 0.


@pytest.mark.parametrize('advance', [True, False])
def test_reactor_rechecks_argon_fraction_at_accepted_steps(tmp_path, advance):
    reactor, _, _ = build_reactor(tmp_path, composition=False, excited_fraction=9.9998799e-6,
                                   argon_validity={'excited_fraction_max': 1.e-5})
    # Ionization consumes ground Ar. The fraction crosses its declared bound
    # while the excited population stays fixed and E/N remains inside the table.
    with pytest.raises(PlasmaStateError, match='excited_fraction_max'):
        if advance:
            reactor.advance(1.e-9)
        else:
            for _ in range(1000):
                reactor.step(1.e-9)
    assert reactor.t > 0.
