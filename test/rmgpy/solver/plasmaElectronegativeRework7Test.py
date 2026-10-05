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


"""Parameter refusal and complete electronegative residual regressions."""
import copy
import numpy as np
import pytest
from rmgpy import constants
from rmgpy.exceptions import PlasmaStateError, ElectronegativeWallRegimeError
from rmgpy.kinetics import Arrhenius, TwoTemperaturePlasma
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.solver.plasma import PlasmaReactor
from plasmaElectronegativeWallTest import model, FixtureEnergyReactor
from plasmaElectronegativeRework4Test import set_state
from plasmaElectronegativeRework6Test import decimal_operator, reinitialise, add_neutral_chemistry


def unchanged_last_valid(r, previous):
    assert r.electronegative_wall_last_valid_state[0] == previous[0]
    np.testing.assert_array_equal(r.electronegative_wall_last_valid_state[1], previous[1])


@pytest.mark.parametrize('mobility', [1.e-307, 1.e285])
@pytest.mark.parametrize('charge', ['anion', 'cation'])
def test_extreme_mobility_refused_at_configuration(mobility, charge):
    key, label = ('anion_reduced_mobilities', 'Cl-') if charge == 'anion' else ('ion_reduced_mobilities', 'Ar+')
    with pytest.raises(PlasmaStateError, match='domain:.*mobility.*expected') as error:
        model(**{key: {label: (mobility, 'm^2/(V*s)')}})
    print('PARAMETER_REFUSAL', charge, mobility, str(error.value))


@pytest.mark.parametrize('energy', [False, True])
@pytest.mark.parametrize('kind', ['ordered-triple', 'attachment', 'tiny-mobility-destruction'])
def test_r128_parameters_refuse_without_publication(energy, kind):
    r, sp, reactions = add_neutral_chemistry(energy)
    r.kf[-1] = 1.
    if kind == 'ordered-triple':
        recombination = LibraryReaction(reactants=[sp[0], sp[2]], products=[sp[1]],
            reversible=False, library='EN-test', kinetics=Arrhenius(A=(2., 'm^3/(mol*s)')))
        reactions.insert(-1, recombination)
        reinitialise(r, sp, reactions)
        r.kf[-2], r.kf[-1] = 1.e18, 1.e308
    elif kind == 'attachment':
        r.kf[1] = 1.e35
    else:
        # Constructor refusal is the primary regression; mutate the public
        # declaration as a defensive accepted-state control, too.
        r.anion_reduced_mobilities['Cl-'] = (1.e-307, 'm^2/(V*s)')
        r.kf[0] = 1.e-307
    y = set_state(r, [1.e-30, 1.e6, 1.e-6, 1.e-30, 1.e-6, 1.e-30, 0.])
    previous = copy.deepcopy(r.electronegative_wall_last_valid_state)
    with pytest.raises(PlasmaStateError, match='domain:.*expected') as error:
        r.monitor_electronegative_wall(y, 7.)
    unchanged_last_valid(r, previous)
    print('R128_REFUSAL', energy, kind, str(error.value), 'last_valid', previous[0])


@pytest.mark.parametrize('energy', [False, True])
def test_constant_rate_refused_during_configuration(energy):
    r, sp, reactions = model(energy=energy)
    reactions[1].kinetics.A.value_si = 1.e35
    with pytest.raises(PlasmaStateError, match='domain:.*Cl.*1e\+35.*expected') as error:
        reinitialise(r, sp, reactions)
    print('CONSTANT_REFUSAL', energy, str(error.value))


@pytest.mark.parametrize('prefactor', [1.e20, 1.e14])
def test_signed_energy_cost_keeps_cancellation_remainder(prefactor):
    r, sp, reactions = model(energy=True)
    reactions.append(LibraryReaction(reactants=[sp[2], sp[4]], products=[sp[1], sp[3]],
        reversible=False, kinetics=Arrhenius(A=(1.e12, 'm^3/(mol*s)')), library='EN-test'))
    law = TwoTemperaturePlasma(A=(1., 'm^3/(mol*s)'), n=-1., Ea_g=(0., 'J/mol'),
        Ea_e=(0., 'J/mol'), T0=(30000., 'K'), electrons=1)
    reactions.append(LibraryReaction(reactants=[sp[0], sp[1]], products=[sp[2], sp[0], sp[0]],
        reversible=False, library='EN-test', kinetics=law))
    reinitialise(r, sp, reactions)
    y = set_state(r, [1.e-30, 1., 1.e-6, .001, 1.e-6])
    r.kf[0] = 0.
    law.A.value_si = prefactor
    r.kf[-1] = r._evaluate_two_temperature_rate_coefficient(law)
    r.energy_threshold[-1] = -434184.05011703185
    previous = copy.deepcopy(r.electronegative_wall_last_valid_state)
    if prefactor == 1.e20:
        print('SIGNED_ENERGY_ORIGINAL', 'k', r.kf[-1],
            'actual', r.residual(7., y, np.zeros_like(y))[0][-1],
            'Decimal', decimal_operator(r, y, full=True, values=True)[-1])
        with pytest.raises(PlasmaStateError, match='domain:.*expected') as error:
            r.monitor_electronegative_wall(y, 7.)
        unchanged_last_valid(r, previous)
        print('SIGNED_ENERGY_OUTSIDE', 'k', r.kf[-1], str(error.value))
        return
    assert 0. < r.kf[-1] <= 1.e15
    expected = decimal_operator(r, y, full=True, values=True)[-1]
    actual = r.residual(7., y, np.zeros_like(y))[0][-1]
    r.monitor_electronegative_wall(y, 7.)
    print('SIGNED_ENERGY_INSIDE', 'k', r.kf[-1], 'actual', actual, 'Decimal', expected,
        'last_valid', r.electronegative_wall_last_valid_state[0])
    assert actual == pytest.approx(expected, rel=2.e-10, abs=0.)


class NonfiniteZeroAnionReactor(PlasmaReactor):
    corrupt = False
    def jacobian(self, t, y, dydt, cj):
        result = super().jacobian(t, y, dydt, cj)
        if self.corrupt:
            result[3, 1] = np.inf
        return result


@pytest.mark.parametrize('energy', [False, True])
@pytest.mark.parametrize('component', ['Jacobian', 'residual'])
def test_zero_anion_complete_operator_precedes_publication(energy, component, monkeypatch):
    r, sp, reactions = model(energy=energy, reactor_cls=NonfiniteZeroAnionReactor)
    # model's energy fixture selects its library subclass. Patch that Python
    # subclass's seam only while this test runs, retaining its real operator.
    if energy:
        original = FixtureEnergyReactor.jacobian
        def jacobian(self, t, y, dydt, cj):
            result = original(self, t, y, dydt, cj)
            if t == 7.:
                result[3, 1] = np.inf
            return result
        monkeypatch.setattr(FixtureEnergyReactor, 'jacobian', jacobian)
    else:
        r.corrupt = True
    y = set_state(r, [1.e-6, 1., 1.e-6, .001, 0.])
    previous = copy.deepcopy(r.electronegative_wall_last_valid_state)
    if component == 'residual':
        # Direct seam exercises the optional precomputed residual path.
        result = np.zeros(len(y)); result[3] = np.inf
        call = lambda: r._evaluate_electronegative_wall_regime(y, 7., result)
    else:
        call = lambda: r.monitor_electronegative_wall(y, 7.)
    with pytest.raises(ElectronegativeWallRegimeError, match=component) as error:
        call()
    unchanged_last_valid(r, previous)
    print('ZERO_ANION_REFUSAL', energy, component, str(error.value), 'last_valid', previous[0])


@pytest.mark.parametrize('energy', [False, True])
@pytest.mark.parametrize('reverse', [False, True])
def test_catalyst_has_no_residual_round_trip(energy, reverse):
    r, sp, reactions = add_neutral_chemistry(energy, reverse=reverse)
    recombination = LibraryReaction(reactants=[sp[0], sp[2]], products=[sp[1]],
        reversible=False, library='EN-test', kinetics=Arrhenius(A=(2.e9, 'm^3/(mol*s)')))
    reactions.insert(-1, recombination)
    reinitialise(r, sp, reactions)
    r.kf[-1], r.kb[-1] = (0., 1.e15) if reverse else (1.e15, 0.)
    r.wall_recycling = 0.
    y = set_state(r, [1.e-30, 1., 1.e-6, .001, 1.e-6, 1.e-5, 0.])
    expected = decimal_operator(r, y, full=True, values=True)[1]
    actual = r.residual(7., y, np.zeros_like(y))[0][1]
    r.monitor_electronegative_wall(y, 7.)
    print('CATALYST', energy, reverse, 'actual', actual, 'Decimal', expected,
        'last_valid', r.electronegative_wall_last_valid_state[0])
    assert actual == pytest.approx(expected, rel=2.e-10, abs=0.)


def test_evolving_rate_refused_after_te_refresh():
    r, sp, reactions = model(energy=True)
    law = TwoTemperaturePlasma(A=(1.e14, 'm^3/(mol*s)'), n=1., Ea_g=(0., 'J/mol'),
        Ea_e=(0., 'J/mol'), T0=(30000., 'K'), electrons=1)
    reactions.append(LibraryReaction(reactants=[sp[0], sp[1]], products=[sp[2], sp[0], sp[0]],
        reversible=False, library='EN-test', kinetics=law))
    reinitialise(r, sp, reactions)
    previous = copy.deepcopy(r.electronegative_wall_last_valid_state)
    y = set_state(r, [1.e-6, 1., 1.5e-6, .001, .5e-6], 50.)
    r.Te.value_si *= .99
    with pytest.raises(PlasmaStateError, match='domain:.*forward rate coefficient.*expected') as error:
        r.monitor_electronegative_wall(y, 7.)
    unchanged_last_valid(r, previous)
    print('EVOLVING_REFUSAL', 'k', r.kf[-1], str(error.value))


@pytest.mark.parametrize('coefficient', [np.nan, np.inf, -1., 1.e15])
def test_complete_cache_checks_zero_anions(coefficient):
    r, _, _ = model()
    y = set_state(r, [1.e-6, 1., 1.e-6, .001, 0.])
    r.kf[0] = coefficient
    previous = copy.deepcopy(r.electronegative_wall_last_valid_state)
    if coefficient == 1.e15:
        r.monitor_electronegative_wall(y, 7.)
        assert r.electronegative_wall_last_valid_state[0] == 7.
    else:
        with pytest.raises(PlasmaStateError, match='domain:.*rate coefficient.*expected'):
            r.monitor_electronegative_wall(y, 7.)
        unchanged_last_valid(r, previous)


class TinyDestructionReactor(PlasmaReactor):
    def initialize_model(self, core, reactions, edge, edge_reactions, **kwargs):
        reactions[0].kinetics.A.value_si = 1.e-307
        return super().initialize_model(core, reactions, edge, edge_reactions, **kwargs)


def test_tiny_mobility_and_destruction_refused_at_configuration():
    with pytest.raises(PlasmaStateError, match='domain:.*mobility.*1e-307.*expected') as error:
        model(reactor_cls=TinyDestructionReactor,
            anion_reduced_mobilities={'Cl-': (1.e-307, 'm^2/(V*s)')})
    print('TINY_PAIR_CONFIGURATION', str(error.value))


@pytest.mark.parametrize('energy', [False, True])
def test_zero_anion_original_infinite_chemical_operator_refuses(energy):
    r, sp, reactions = add_neutral_chemistry(energy)
    reactions[-1].reactants = [sp[5], sp[3]]
    reactions[-1].products = [sp[6]]
    reactions[-1].kinetics = Arrhenius(A=(2., 'm^3/(mol*s)'))
    reinitialise(r, sp, reactions)
    r.kf[-1] = 1.e308
    y = set_state(r, [1.e-6, 1., 1.e-6, 2., 0., 1.e-30, 0.])
    previous = copy.deepcopy(r.electronegative_wall_last_valid_state)
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        operator = r.jacobian(7., y, np.zeros_like(y), 0.)
    assert np.isinf(operator[6, 5])
    print('ZERO_ANION_ORIGINAL', energy, 'J[NCl,N]', operator[6, 5])
    with pytest.raises(PlasmaStateError, match='domain:.*rate coefficient.*expected'):
        r.monitor_electronegative_wall(y, 7.)
    unchanged_last_valid(r, previous)


def test_published_signed_inelastic_budget_keeps_cancellation_remainder():
    from decimal import Decimal, localcontext
    r, _, _ = model(energy=True)
    y = set_state(r, [1.e-30, 1., 1.e-6, .001, 1.e-6])
    r.energy_threshold[1] = -434184.05011703185
    r.monitor_electronegative_wall(y, 7.)
    with localcontext() as ctx:
        ctx.prec = 100
        D = lambda x: Decimal.from_float(float(x))
        cost = (D(r.energy_threshold[1])+D(r.energy_electrons_consumed[1])
            *D(1.5)*D(constants.R)*D(r.Te.value_si))
        expected = float(cost*D(r.core_reaction_rates[1])*D(r.compute_volume(y)))
    actual = r.electronegative_wall_diagnostics['energy_budget']['Q_inelastic']
    print('SIGNED_PUBLISHED_BUDGET', 'actual', actual, 'Decimal', expected,
        'last_valid', r.electronegative_wall_last_valid_state[0])
    assert r.electronegative_wall_last_valid_state[0] == 7.
    assert actual == pytest.approx(expected, rel=2.e-10, abs=0.)
