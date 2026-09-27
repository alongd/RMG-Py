"""Engine-side acceptance tests for the PlasmaArgonDimer database library."""

import os
import shutil

import pytest

from rmgpy import constants, settings
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.thermo import ThermoDatabase


DATABASE_INPUT = settings['database.directory']
LIBRARY = 'PlasmaArgonDimer'

if not os.path.isdir(os.path.join(DATABASE_INPUT, 'kinetics', 'libraries', LIBRARY)):
    pytest.skip('kinetics library {0!r} is not present in the configured database ({1}); this '
                'suite exercises the real library and does not inject a substitute.'.format(
                    LIBRARY, DATABASE_INPUT), allow_module_level=True)


@pytest.fixture(scope='module')
def library():
    db = KineticsDatabase()
    db.load_libraries(os.path.join(DATABASE_INPUT, 'kinetics', 'libraries'), libraries=[LIBRARY])
    return db.libraries[LIBRARY]


def test_real_library_loads_with_three_balanced_entries_and_explicit_electrons(library):
    assert len(library.entries) == 3
    def key(reaction):
        return (tuple(sorted(species.label for species in reaction.reactants)),
                tuple(sorted(species.label for species in reaction.products)))
    reactions = {key(reaction): reaction for reaction in library.get_library_reactions()}
    assert set(reactions) == {
        (('Ar', 'Ar', 'Arp'), ('Ar', 'Ar2p')),
        (('Ar2p', 'e-'), ('Ar', 'Ars')),
        (('Ar2p', 'e-'), ('Ar', 'Arp', 'e-')),
    }
    for reaction in reactions.values():
        assert reaction.is_balanced()
        assert reaction.electrons == 0
    dr = reactions[(('Ar2p', 'e-'), ('Ar', 'Ars'))]
    dei = reactions[(('Ar2p', 'e-'), ('Ar', 'Arp', 'e-'))]
    assert sum(species.is_electron() for species in dr.reactants) == 1
    assert sum(species.is_electron() for species in dr.products) == 0
    assert sum(species.is_electron() for species in dei.reactants) == 1
    assert sum(species.is_electron() for species in dei.products) == 1


def test_dr_unit_gate_is_preserved_by_loaded_engine(library):
    kinetics = library.entries[1].data
    expected = {300: 9.100000e-7, 8500: 1.183422e-7,
                10444: 1.043702e-7, 13462: 8.939815e-8}
    for te, source in expected.items():
        assert kinetics.get_rate_coefficient_two_temp(298.15, te) == pytest.approx(
            source * 1e-6 * constants.Na, rel=1e-6)


def test_loader_refuses_duplicate_explicit_and_metadata_electron(tmp_path):
    source = os.path.join(DATABASE_INPUT, 'kinetics', 'libraries', LIBRARY)
    target = tmp_path / LIBRARY
    shutil.copytree(source, target)
    reaction_file = target / 'reactions.py'
    content = reaction_file.read_text()
    content = content.replace("Tmax=(8500, 'K')),", "Tmax=(8500, 'K'), electrons=-1),", 1)
    reaction_file.write_text(content)
    db = KineticsDatabase()
    with pytest.raises(Exception, match='balanced|electron'):
        db.load_libraries(str(tmp_path), libraries=[LIBRARY])


def test_termolecular_reverse_fit_can_be_evaluated_at_1000_k(library):
    reaction = next(r for r in library.get_library_reactions() if r.reversible)
    thermo = ThermoDatabase()
    thermo.load_libraries(os.path.join(DATABASE_INPUT, 'thermo', 'libraries'),
                          libraries=['PlasmaThermo', 'primaryThermoLibrary'])
    thermo.load_groups(os.path.join(DATABASE_INPUT, 'thermo', 'groups'))
    for species in reaction.reactants + reaction.products:
        species.thermo = thermo.get_thermo_data(species)
    # The collision-limit diagnostic requests a reverse fit on its own grid.  The
    # source Tmax remains 300 K; no generic range guard may turn the 1000 K call into
    # an exception.
    reverse = reaction.generate_reverse_rate_coefficient()
    assert reverse.get_rate_coefficient(1000.0) > 0.0
