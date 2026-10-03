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


"""Cache invalidation through public reaction references and monitoring."""
import copy
from decimal import Decimal, localcontext

import numpy as np
import pytest
import rmgpy.data.rmg as data
from rmgpy.exceptions import PlasmaStateError
from rmgpy.solver.electronegative import EN_WALL_DOMAIN
from plasmaElectronegativeRework9Test import declared_model
from plasmaElectronegativeRework6Test import reinitialise


@pytest.fixture(autouse=True)
def isolated_database():
    previous = data.database
    data.database = None
    yield
    data.database = previous


def appreciable_model(prefactor=1e8):
    reactor, species, reactions = declared_model(0., 1.)
    reactions[-1].kinetics.A.value_si = prefactor
    reinitialise(reactor, species, reactions)
    return reactor, species, reactions


def coefficient(reactor, state, index):
    return next(row['coefficient_si'] for row in reactor.compute_reference_reaction_data(state)
                if row['reaction_index'] == index)


def test_same_temperature_skip_uses_the_state_law():
    reactor, _, _ = appreciable_model()
    index, law = reactor.energy_te_rate_refresh[-1]
    state = reactor.y.copy()
    state[-1] *= 2.
    reactor.Te.value_si = state[-1]
    expected = reactor.evaluate_two_temperature_rate_coefficient(law)
    actual = coefficient(reactor, state, index)
    np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=0.)
    reactor.monitor_electronegative_wall(state, 7.)
    residual = reactor.residual(7., state, np.zeros_like(state))[0].copy()
    with localcontext() as context:
        context.prec = 200
        D = lambda value: Decimal.from_float(float(value))
        expected_rate = float(D(expected)*D(state[0])*D(state[1])/D(reactor.compute_volume(state)))
    record = next(row for row in reactor.compute_reference_reaction_data(state)
                  if row['reaction_index'] == index)
    np.testing.assert_allclose(record['rate_mol_per_s'], expected_rate, rtol=1e-13)
    near = state.copy()
    near[-1] = np.nextafter(state[-1], float('inf'))
    reactor._sync_electron_temperature(near)
    reactor._sync_electron_temperature(state)
    corrected = reactor.residual(7., state, np.zeros_like(state))[0]
    np.testing.assert_array_equal(residual, corrected)
    assert reactor.electronegative_wall_last_valid_state[0] == 7.
    print('SAME_TEMPERATURE', actual, expected, 'residual_delta', corrected[0]-residual[0])


@pytest.mark.parametrize('prefactor', [float('nan'), float('inf'), 1e16], ids=['nan', 'infinity', 'cap'])
@pytest.mark.parametrize('api', ['reference', 'monitor'])
def test_cached_law_mutation_refuses_without_publication(prefactor, api):
    reactor, _, _ = appreciable_model()
    index, original = reactor.energy_te_rate_refresh[-1]
    law = copy.deepcopy(original)
    law.A.value_si = prefactor
    reactor.energy_te_rate_refresh[-1] = (index, law)
    previous = copy.deepcopy(reactor.electronegative_wall_last_valid_state)
    diagnostics = copy.deepcopy(reactor.electronegative_wall_diagnostics)
    with pytest.raises(PlasmaStateError, match='rate coefficient') as exc:
        if api == 'reference':
            reactor.compute_reference_reaction_data(reactor.y.copy())
        else:
            reactor.monitor_electronegative_wall(reactor.y.copy(), 7.)
    assert previous[0] == reactor.electronegative_wall_last_valid_state[0]
    np.testing.assert_array_equal(previous[1], reactor.electronegative_wall_last_valid_state[1])
    np.testing.assert_equal(diagnostics, reactor.electronegative_wall_diagnostics)
    print('CACHED_LAW_REFUSAL', api, prefactor, str(exc.value), 'last_valid', previous[0])


def test_same_temperature_change_rechecks_rate_cap():
    reactor, _, _ = appreciable_model(9e14)
    state = reactor.y.copy()
    state[-1] *= 2.
    reactor.Te.value_si = state[-1]
    with pytest.raises(PlasmaStateError, match='rate coefficient.*1e\\+15'):
        reactor.monitor_electronegative_wall(state, 7.)
    assert reactor.electronegative_wall_last_valid_state[0] == 0.


@pytest.mark.parametrize('parameter', ['A', 'n', 'T0', 'Ea_g', 'Ea_e', 'Tg', 'law'])
def test_cache_key_covers_live_law_inputs(parameter):
    reactor, _, reactions = appreciable_model()
    index, law = reactor.energy_te_rate_refresh[-1]
    state = reactor.y.copy()
    if parameter == 'law':
        law = copy.deepcopy(law)
        law.A.value_si *= 2.
        reactions[-1].kinetics = law
    elif parameter == 'Tg':
        law.Ea_g.value_si = 1000.
        coefficient(reactor, state, index)
        reactor.T.value_si *= 1.1
    elif parameter == 'n':
        state[-1] *= 2.
        coefficient(reactor, state, index)
        law.n.value_si = 2.
    elif parameter in ('Ea_g', 'Ea_e'):
        getattr(law, parameter).value_si = 1000.
    else:
        getattr(law, parameter).value_si *= 2.
    # The state is authoritative in energy mode, even if public Te was changed.
    reactor.Te.value_si = state[-1]
    expected = law.get_rate_coefficient_two_temp(reactor.T.value_si, state[-1])
    np.testing.assert_allclose(coefficient(reactor, state, index), expected, rtol=1e-13, atol=0.)
    print('CACHE_KEY_INPUT', parameter, expected)


def test_changed_rate_cap_is_checked_at_unchanged_temperature(monkeypatch):
    reactor, _, _ = appreciable_model()
    monkeypatch.setitem(EN_WALL_DOMAIN, 'rate_coefficient_si', (0., 1e7))
    with pytest.raises(PlasmaStateError, match='rate coefficient.*1e\\+07'):
        reactor.monitor_electronegative_wall(reactor.y.copy(), 7.)
    assert reactor.electronegative_wall_last_valid_state[0] == 0.


@pytest.mark.parametrize('parameter', ['Te', 'A', 'Tg'])
def test_prescribed_temperature_cache_inputs(parameter):
    from plasmaEnergyBalanceTest import _build
    reactor, _, reactions = _build(energy=False, elastic=False)
    state = reactor.y.copy()
    law = reactions[0].kinetics
    if parameter == 'Te':
        reactor.Te.value_si *= 2.
    elif parameter == 'Tg':
        law.Ea_g.value_si += 1000.
        coefficient(reactor, state, 0)
        reactor.T.value_si *= 1.1
    else:
        law.A.value_si *= 2.
    expected = law.get_rate_coefficient_two_temp(reactor.T.value_si, reactor.Te.value_si)
    np.testing.assert_allclose(coefficient(reactor, state, 0), expected, rtol=1e-13, atol=0.)
    reactor.residual(0., state, np.zeros_like(state))
    np.testing.assert_allclose(reactor.kf[0], expected, rtol=1e-13, atol=0.)


def test_unknown_evaluator_parameters_are_never_cached():
    reactor, _, _ = appreciable_model()
    index, original = reactor.energy_te_rate_refresh[-1]

    class MutableEvaluator:
        uses_electron_temperature = True
        factor = 1.

        def get_rate_coefficient_two_temp(self, tg, te):
            return self.factor*original.get_rate_coefficient_two_temp(tg, te)

    law = MutableEvaluator()
    reactor.energy_te_rate_refresh[-1] = (index, law)
    initial = coefficient(reactor, reactor.y.copy(), index)
    law.factor = 2.
    np.testing.assert_allclose(coefficient(reactor, reactor.y.copy(), index), 2.*initial, rtol=1e-13)



def test_live_law_replacement_after_canonical_rekey():
    reactor, _, reactions = appreciable_model()
    model_reactions = [copy.copy(reaction) for reaction in reactions]
    reactor._rekey_reaction_index_to_model(model_reactions, [])
    index, original = reactor.energy_te_rate_refresh[-1]
    law = copy.deepcopy(original)
    law.A.value_si *= 2.
    model_reactions[-1].kinetics = law
    expected = law.get_rate_coefficient_two_temp(reactor.T.value_si, reactor.Te.value_si)
    np.testing.assert_allclose(coefficient(reactor, reactor.y.copy(), index), expected,
                               rtol=1e-13, atol=0.)
