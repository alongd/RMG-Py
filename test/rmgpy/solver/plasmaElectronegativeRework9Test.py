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


"""Reaction-parameter bounds, standalone temperature and scaled-energy witnesses."""
import copy
import numpy as np
import pytest
import rmgpy.data.rmg as data
from rmgpy import constants
from rmgpy.exceptions import PlasmaStateError, ElectronegativeWallRegimeError
from rmgpy.kinetics import Arrhenius, TwoTemperaturePlasma, MultiArrhenius
from rmgpy.data.kinetics.library import LibraryReaction
from plasmaElectronegativeWallTest import model, FixtureEnergyReactor, qualification
from plasmaElectronegativeRework6Test import reinitialise, decimal_operator

EVK = constants.e*constants.Na/constants.R

@pytest.fixture(autouse=True)
def isolated_database():
    previous = data.database
    data.database = None
    yield
    data.database = previous


def declared_model(threshold=0., exponent=0.):
    seed, species, reactions = model(energy=True)
    law = TwoTemperaturePlasma(A=(1e-100, 'm^3/(mol*s)'), n=exponent,
        Ea_g=(0., 'J/mol'), Ea_e=(0., 'J/mol'), T0=(seed.Te.value_si, 'K'), electrons=1)
    reactions.append(LibraryReaction(reactants=[species[0], species[1]],
        products=[species[2], species[0], species[0]], reversible=False,
        library='EN-test', kinetics=law))
    config = copy.deepcopy(seed.electron_energy_balance)
    config['electron_energies']['EN-test:3'] = (threshold, 'J/mol')
    FixtureEnergyReactor.entries = {i+1: reaction.kinetics for i, reaction in enumerate(reactions)}
    r = FixtureEnergyReactor((seed.T.value_si, 'K'), (seed.P.value_si, 'Pa'),
        dict(seed.initial_mole_fractions), (seed.Te.value_si, 'K'),
        diffusion_length=(seed.diffusion_length.value_si, 'm'),
        ion_reduced_mobilities={'Ar+': (1.535e-4, 'm^2/(V*s)')},
        anion_reduced_mobilities={'Cl-': (1.5e-4, 'm^2/(V*s)')},
        wall_diffusion_components=tuple(seed.wall_diffusion_components),
        wall_chamber_geometry=dict(seed.wall_chamber_geometry),
        electronegative_wall_model='confinedAnion',
        electronegative_wall_geometry='fullFrequency',
        electronegative_wall_qualification=qualification(),
        wall_single_bath_approximation=True, thermo_source_assertions={'Ar+': 'ion', 'Cl-': 'ion'},
        electron_energy_balance=config)
    reinitialise(r, species, reactions)
    return r, species, reactions


@pytest.mark.parametrize('threshold', (1e307, -1e307))
def test_huge_energy_is_a_named_configuration_refusal(threshold):
    with pytest.raises(PlasmaStateError, match=r'EN-test:3.*1e\+307.*1e\+08') as exc:
        declared_model(threshold, 100.)
    print('ENERGY_CONFIGURATION_REFUSAL', str(exc.value))


@pytest.mark.parametrize('exponent', (100., -100.))
def test_rate_exponent_is_a_named_configuration_refusal(exponent):
    with pytest.raises(PlasmaStateError, match=r'e- \+ Ar.*temperature exponent.*100.*50') as exc:
        declared_model(0., exponent)
    print('EXPONENT_CONFIGURATION_REFUSAL', str(exc.value))


@pytest.mark.parametrize('threshold', (-1e8, 1e8))
@pytest.mark.parametrize('exponent', (-50., 50.))
def test_new_bounds_are_inclusive(threshold, exponent):
    r, _, _ = declared_model(threshold, exponent)
    r.monitor_electronegative_wall(r.y.copy(), 7.)
    assert r.electronegative_wall_last_valid_state[0] == 7.
    print('PARAMETER_ENDPOINT', threshold, exponent)


@pytest.mark.parametrize('parameter', ('threshold', 'exponent', 'replacement', 'composite', 'refresh'))
@pytest.mark.parametrize('method', ('monitor', 'reference'))
def test_live_parameter_mutation_refuses_without_publication(parameter, method):
    r, _, reactions = declared_model()
    previous = copy.deepcopy(r.electronegative_wall_last_valid_state)
    diagnostics = copy.deepcopy(r.electronegative_wall_diagnostics)
    if parameter == 'threshold':
        r.energy_threshold[-1] = 1e307
        match = r'e- \+ Ar.*energy.*1e\+307.*1e\+08'
    elif parameter == 'exponent':
        reactions[-1].kinetics.n.value_si = 100.
        match = r'e- \+ Ar.*temperature exponent.*100.*50'
    elif parameter == 'refresh':
        index, original = r.energy_te_rate_refresh[-1]
        law = copy.deepcopy(original)
        law.n.value_si = 100.
        r.energy_te_rate_refresh[-1] = (index, law)
        match = r'e- \+ Ar.*temperature exponent.*100.*50'
    else:
        law = Arrhenius(A=(1e-100, 'm^3/(mol*s)'), n=100.)
        reactions[-1].kinetics = MultiArrhenius(arrhenius=[law]) if parameter == 'composite' else law
        match = r'e- \+ Ar.*temperature exponent.*100.*50'
    with pytest.raises(PlasmaStateError, match=match) as exc:
        if method == 'monitor': r.monitor_electronegative_wall(r.y.copy(), 7.)
        else: r.compute_reference_reaction_data(r.y.copy())
    assert previous[0] == r.electronegative_wall_last_valid_state[0]
    np.testing.assert_array_equal(previous[1], r.electronegative_wall_last_valid_state[1])
    np.testing.assert_equal(diagnostics, r.electronegative_wall_diagnostics)
    print('PARAMETER_MUTATION_REFUSAL', parameter, method, str(exc.value), 'last_valid', previous[0])


@pytest.mark.parametrize('energy', (False, True))
@pytest.mark.parametrize('temperature', (100., float('nan')))
def test_standalone_rate_checks_the_temperature_it_uses(energy, temperature):
    r, _, _ = model(energy=energy)
    previous = copy.deepcopy(r.electronegative_wall_last_valid_state)
    r.Te.value_si = temperature*EVK
    law = TwoTemperaturePlasma(A=(1., 'm^3/(mol*s)'), n=0.,
        Ea_g=(0., 'J/mol'), Ea_e=(0., 'J/mol'), T0=(30000., 'K'), electrons=1)
    with pytest.raises(ElectronegativeWallRegimeError, match=r'Te.*0.05.*50') as exc:
        r.evaluate_two_temperature_rate_coefficient(law)
    assert previous[0] == r.electronegative_wall_last_valid_state[0]
    np.testing.assert_array_equal(previous[1], r.electronegative_wall_last_valid_state[1])
    print('STANDALONE_TE_REFUSAL', energy, temperature, str(exc.value), 'last_valid', previous[0])


def test_standalone_rate_checks_exponent_before_evaluation():
    r, _, _ = model()
    law = TwoTemperaturePlasma(A=(1., 'm^3/(mol*s)'), n=100.,
        Ea_g=(0., 'J/mol'), Ea_e=(0., 'J/mol'), T0=(r.Te.value_si, 'K'), electrons=1)
    with pytest.raises(PlasmaStateError, match=r'temperature exponent.*100.*50'):
        r.evaluate_two_temperature_rate_coefficient(law)


def test_private_tiny_rate_prevents_intermediate_energy_overflow():
    # Outside the new parameter domain: probe constitutive arithmetic only.
    r, _, reactions = declared_model()
    reactions[-1].kinetics.n.value_si = 100.
    r.energy_threshold[-1] = 1e307
    y = r.y.copy()
    actual = r.jacobian(7., y, np.zeros_like(y), 0.)[-1, -1]
    assert np.isfinite(actual)
    for precision in (200, 800):
        expected = decimal_operator(r, y, full=True, precision=precision)[-1, -1]
        np.testing.assert_allclose(actual, expected, rtol=2e-10, atol=0.)
        print('SCALED_ENERGY_REFERENCE', precision, actual, expected)
    with pytest.raises(PlasmaStateError, match='reaction_energy_J_mol'):
        r.monitor_electronegative_wall(y, 7.)


@pytest.mark.parametrize('exponent', (-51., 51.))
def test_elastic_rate_exponent_is_also_bounded(exponent):
    r, _, _ = model(energy=True)
    config = copy.deepcopy(r.electron_energy_balance)
    config['elastic_collisions']['Ar'] = dict(A=(1e-20, 'm^3/s'), n=exponent, b=0., c=0.)
    with pytest.raises(PlasmaStateError, match=r'Ar.*temperature exponent.*51.*50'):
        model(energy=True, electron_energy_balance=config)


@pytest.mark.parametrize('temperature', (.05, 50.))
def test_standalone_rate_accepts_the_existing_temperature_endpoints(temperature):
    r, _, _ = model(energy=False)
    r.Te.value_si = temperature*EVK
    law = TwoTemperaturePlasma(A=(1., 'm^3/(mol*s)'), n=0.,
        Ea_g=(0., 'J/mol'), Ea_e=(0., 'J/mol'), T0=(30000., 'K'), electrons=1)
    assert r.evaluate_two_temperature_rate_coefficient(law) == 1.
