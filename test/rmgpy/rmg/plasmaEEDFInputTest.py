from copy import deepcopy

import pytest

import rmgpy.rmg.input as inp
from rmgpy.exceptions import InputError
from rmgpy.rmg.main import RMG


SHA256 = "A" * 64


def _electron_kinetics(**updates):
    declaration = {
        'provider': 'loki-table',
        'table': ('eedf-table.h5', SHA256),
        'branch': 'B1',
        'operatingBranch': {'id': 'reactor-000', 'path': 'branches.json'},
        'initialReducedField': (16.0, 'Td'),
        'empiricalLaws': {
            'PlasmaArgon:90': {
                'evaluate_at': 'Te_eff',
                'class': 'C',
                'basis': 'I-318 row 7',
                'sensitivity': True,
            },
        },
    }
    declaration.update(updates)
    return declaration


def _preamble():
    return (
        "species(label='e-', structure=adjacencyList(\"1 e u1 p0 c-1\"))\n"
        "species(label='Ar', structure=adjacencyList(\"1 Ar u0 p4 c0\"))\n"
    )


def _read(tmp_path, reactor_block):
    path = tmp_path / 'input.py'
    path.write_text(_preamble() + reactor_block)
    rmg = RMG()
    inp.read_input_file(str(path), rmg)
    return rmg


def _eedf_block(electron_kinetics=None, energy_balance=None, electron_temperature='',
                extra='', include_energy_balance=True):
    electron_kinetics = _electron_kinetics() if electron_kinetics is None else electron_kinetics
    energy_balance = ({'absorbedPower': (0.5, 'W'), 'sheath': 'floatingWall'}
                      if energy_balance is None else energy_balance)
    energy_line = (f"    electronEnergyBalance={energy_balance!r},\n"
                   if include_energy_balance else '')
    return (
        "plasmaReactor(\n"
        "    temperature=(298.15, 'K'),\n"
        "    pressure=(5, 'torr'),\n"
        f"{electron_temperature}"
        f"{extra}"
        "    initialMoleFractions={'Ar': 0.999999, 'e-': 1e-6},\n"
        "    chamberGeometry={'shape': 'cylinder', 'radius': (5, 'cm'), "
        "'length': (30, 'cm')},\n"
        "    ionReducedMobility=(1.535e-4, 'm^2/(V*s)'),\n"
        f"    electronKinetics={electron_kinetics!r},\n"
        f"{energy_line}"
        "    terminationTime=(1, 's'),\n"
        ")\n"
    )


def test_public_eedf_declaration_is_normalized_for_the_constructor(tmp_path, monkeypatch):
    captured = {}

    class ReactorCapture:
        def __init__(self, **kwargs):
            captured.update(kwargs)

        def convert_initial_keys_to_species_objects(self, species_dict):
            pass

    monkeypatch.setattr(inp, 'PlasmaReactor', ReactorCapture)
    _read(tmp_path, _eedf_block())

    assert captured['Te'] == (298.15, 'K')
    assert captured['electron_kinetics'] == {
        'provider': 'loki-table',
        'table': (str(tmp_path / 'eedf-table.h5'), SHA256.lower()),
        'branch': 'B1',
        'operating_branch': {'id': 'reactor-000', 'path': str(tmp_path / 'branches.json')},
        'initial_reduced_field': (16.0, 'Td'),
        'empirical_laws': {
            'PlasmaArgon:90': {
                'evaluate_at': 'Te_eff',
                'class': 'C',
                'basis': 'I-318 row 7',
                'sensitivity': True,
            },
        },
    }
    assert captured['electron_energy_balance'].keys() == {
        'absorbed_power', 'chamber_volume', 'sheath'}
    assert captured['electron_energy_balance']['sheath'] == 'floating_wall'


def test_legacy_declaration_does_not_pass_electron_kinetics(tmp_path, monkeypatch):
    captured = {}

    class ReactorCapture:
        def __init__(self, **kwargs):
            captured.update(kwargs)

        def convert_initial_keys_to_species_objects(self, species_dict):
            pass

    monkeypatch.setattr(inp, 'PlasmaReactor', ReactorCapture)
    block = (
        "plasmaReactor(temperature=(298.15, 'K'), pressure=(5, 'torr'),\n"
        "    electronTemperature=(11604.5, 'K'),\n"
        "    initialMoleFractions={'Ar': 0.999999, 'e-': 1e-6},\n"
        "    terminationTime=(1, 's'))\n"
    )
    _read(tmp_path, block)
    assert 'electron_kinetics' not in captured
    assert captured['Te'] == (11604.5, 'K')


def test_legacy_declaration_still_requires_electron_temperature(tmp_path):
    block = (
        "plasmaReactor(temperature=(298.15, 'K'), pressure=(5, 'torr'),\n"
        "    initialMoleFractions={'Ar': 0.999999, 'e-': 1e-6},\n"
        "    terminationTime=(1, 's'))\n"
    )
    with pytest.raises(TypeError, match='electronTemperature'):
        _read(tmp_path, block)


@pytest.mark.parametrize(
    'updates, match',
    [
        ({'provider': 'other'}, 'provider'),
        ({'table': ('table.h5', 'bad')}, 'sha256'),
        ({'branch': ''}, 'branch'),
        ({'operatingBranch': {'id': '', 'path': 'branches.json'}}, 'operatingBranch'),
        ({'initialReducedField': (16.0, 'V')}, 'initialReducedField'),
        ({'initialReducedField': (0.0, 'Td')}, 'initialReducedField'),
        ({'initialReducedField': (True, 'Td')}, 'initialReducedField'),
        ({'empiricalLaws': []}, 'empiricalLaws'),
    ],
)
def test_electron_kinetics_top_level_contract(updates, match):
    with pytest.raises(InputError, match=match):
        inp._plasma_electron_kinetics(_electron_kinetics(**updates))


def test_development_unqualified_route_is_not_an_input_dsl_option():
    with pytest.raises(InputError, match='unsupported key'):
        inp._plasma_electron_kinetics(
            _electron_kinetics(developmentUnqualified=True))


@pytest.mark.parametrize(
    'key, value, match',
    [
        ('evaluate_at', 'Tgas', 'evaluate_at'),
        ('class', 'D', 'class'),
        ('basis', '', 'basis'),
        ('sensitivity', 1, 'sensitivity'),
    ],
)
def test_empirical_law_contract(key, value, match):
    declaration = _electron_kinetics()
    law = deepcopy(declaration['empiricalLaws']['PlasmaArgon:90'])
    law[key] = value
    declaration['empiricalLaws'] = {'PlasmaArgon:90': law}
    with pytest.raises(InputError, match=match):
        inp._plasma_electron_kinetics(declaration)


@pytest.mark.parametrize('entry', ['', ' PlasmaArgon:90', 'PlasmaArgon:90 ', 90])
def test_empirical_law_keys_are_stable_nonempty_strings(entry):
    declaration = _electron_kinetics()
    law = declaration['empiricalLaws'].pop('PlasmaArgon:90')
    declaration['empiricalLaws'][entry] = law
    with pytest.raises(InputError, match='stable strings'):
        inp._plasma_electron_kinetics(declaration)


def test_eedf_refuses_electron_temperature(tmp_path):
    with pytest.raises(InputError, match='electronTemperature'):
        _read(tmp_path, _eedf_block(
            electron_temperature="    electronTemperature=(11604.5, 'K'),\n"))


def test_eedf_requires_energy_balance(tmp_path):
    with pytest.raises(InputError, match='electronEnergyBalance is required'):
        _read(tmp_path, _eedf_block(include_energy_balance=False))


def test_eedf_refuses_electron_density_until_initial_row_inversion_exists(tmp_path):
    with pytest.raises(InputError, match='pressure-inversion mean energy'):
        _read(tmp_path, _eedf_block(extra="    electronDensity=(1e16, 'm^-3'),\n"))


@pytest.mark.parametrize('forbidden', ['elasticCollisions'])
def test_eedf_energy_balance_refuses_table_owned_declarations(tmp_path, forbidden):
    energy = {'absorbedPower': (0.5, 'W'), 'sheath': 'floatingWall', forbidden: {}}
    with pytest.raises(InputError, match=forbidden):
        _read(tmp_path, _eedf_block(energy_balance=energy))


def test_eedf_energy_balance_passes_optional_heavy_particle_electron_energies(tmp_path, monkeypatch):
    captured = {}

    class ReactorCapture:
        def __init__(self, **kwargs):
            captured.update(kwargs)

        def convert_initial_keys_to_species_objects(self, species_dict):
            pass

    monkeypatch.setattr(inp, 'PlasmaReactor', ReactorCapture)
    energy = {'absorbedPower': (0.5, 'W'), 'sheath': 'floatingWall',
              'electronEnergies': {'PlasmaArgon:90': (-7.34, 'eV')}}
    _read(tmp_path, _eedf_block(energy_balance=energy))
    assert captured['electron_energy_balance']['electron_energies'] == {
        'PlasmaArgon:90': (-7.34, 'eV')}


def test_eedf_energy_balance_preserves_geometry_volume_conflict(tmp_path):
    energy = {
        'absorbedPower': (0.5, 'W'),
        'sheath': 'floatingWall',
        'chamberVolume': (1.0, 'm^3'),
    }
    with pytest.raises(InputError, match='two sources of truth'):
        _read(tmp_path, _eedf_block(energy_balance=energy))


def test_real_constructor_save_read_round_trip_uses_deck_relative_table(tmp_path):
    table = tmp_path / 'tables' / 'real.h5'
    table.parent.mkdir()
    table.write_bytes(b'fixture path only; provider loading occurs at initialize_model')
    (tmp_path / 'branches.json').write_text(
        '{"branches": [{"id": "reactor-000", '
        '"seed": {"u": 0.0, "n_e": 1e16}, "u": 0.0, "n_e": 1e16}]}')
    declaration = _electron_kinetics(table=('tables/real.h5', SHA256))
    input_path = tmp_path / 'input.py'
    input_path.write_text(
        "database(thermoLibraries=[], reactionLibraries=[], seedMechanisms=[], "
        "kineticsFamilies=[])\n"
        + _preamble()
        + _eedf_block(electron_kinetics=declaration)
        + "simulator(atol=1e-16, rtol=1e-8)\n"
          "model(toleranceMoveToCore=0.1, toleranceInterruptSimulation=0.1)\n")
    rmg1 = RMG()
    inp.read_input_file(str(input_path), rmg1)
    reactor1 = rmg1.reaction_systems[0]
    assert reactor1.electron_kinetics['table'] == (str(table), SHA256.lower())

    saved = tmp_path / 'saved.py'
    inp.save_input_file(str(saved), rmg1)
    text = saved.read_text()
    assert 'electronKinetics' in text
    assert 'initialReducedField' in text
    assert 'electronTemperature' not in text
    assert 'elasticCollisions' not in text
    assert 'electronEnergies' not in text

    rmg2 = RMG()
    inp.read_input_file(str(saved), rmg2)
    reactor2 = rmg2.reaction_systems[0]
    assert reactor2.electron_kinetics == reactor1.electron_kinetics
    assert set(reactor2.electron_energy_balance) == {
        'absorbed_power', 'chamber_volume', 'sheath'}
