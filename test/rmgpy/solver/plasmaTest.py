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

"""
Tests for the two-temperature :class:`PlasmaReactor`.

The engineering fixture used throughout seeds an explicit electron population
together with an explicit balancing cation (charge-neutral initial state).
It proves reactor state handling; it does not claim that a neutral mechanism
can generate its own first cation.
"""

import logging
import pickle

import numpy as np
import pytest

import rmgpy.constants as constants
from rmgpy.exceptions import ElectronPlacementError, NonEquilibriumReverseRateError, PlasmaStateError
from rmgpy.kinetics import Arrhenius
from rmgpy.kinetics.arrhenius import TwoTemperaturePlasma
from rmgpy.reaction import Reaction
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.solver.simple import SimpleReactor
from rmgpy.species import Species
from rmgpy.thermo import ThermoData

T_GAS = 1000.0    # K
T_E = 11604.5     # K (~1 eV)
P0 = 1.0e5        # Pa

Y_E0 = 1.0e-4     # mol, seeded electrons
Y_ION0 = 1.0e-4   # mol, seeded balancing cation (charge-neutral with Y_E0)
Y_AR0 = 1.0
Y_A0 = 0.1
Y_B0 = 0.05


def _thermo(h298_kj, s298):
    return ThermoData(
        Tdata=([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0], 'K'),
        Cpdata=([20.8] * 7, 'J/(mol*K)'),
        H298=(h298_kj, 'kJ/mol'),
        S298=(s298, 'J/(mol*K)'),
    )


def _species():
    """Fresh species objects for one test (solver indexing is by identity)."""
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    ar_ion = Species(label='Ar+').from_adjacency_list('1 Ar u1 p3 c+1')
    # Two heavy species for the ordinary thermal reversible channel.
    spc_a = Species(label='A').from_adjacency_list('1 Ar u0 p4 c0')
    spc_a.thermo = _thermo(0.0, 150.0)
    spc_b = Species(label='B').from_adjacency_list('1 Ar u0 p4 c0')
    spc_b.thermo = _thermo(-15.0, 145.0)
    return electron, ar, ar_ion, spc_a, spc_b


def _ionization(electron, ar, ar_ion):
    """e + Ar -> Ar+ + e + e; incident-electron order 1, net electron +1."""
    return Reaction(
        reactants=[electron, ar], products=[ar_ion, electron, electron],
        reversible=False,
        kinetics=TwoTemperaturePlasma(A=(1.0e-3, 'm^3/(mol*s)'), n=0.5,
                                      Ea_g=(4000.0, 'J/mol'), Ea_e=(60000.0, 'J/mol')))


def _recombination(electron, ar, ar_ion):
    """e + Ar+ -> Ar; electron-consuming, incident order 1, net electron -1."""
    return Reaction(
        reactants=[electron, ar_ion], products=[ar],
        reversible=False,
        kinetics=TwoTemperaturePlasma(A=(5.0e5, 'm^3/(mol*s)'), n=-0.5,
                                      Ea_g=(0.0, 'J/mol'), Ea_e=(0.0, 'J/mol')))


def _three_body_recombination(electron, ar, ar_ion):
    """e + e + Ar+ -> Ar + e; the separately supplied reverse of ionization."""
    return Reaction(
        reactants=[electron, electron, ar_ion], products=[ar, electron],
        reversible=False,
        kinetics=TwoTemperaturePlasma(A=(1.0e4, 'm^6/(mol^2*s)'), n=-1.0,
                                      Ea_g=(0.0, 'J/mol'), Ea_e=(0.0, 'J/mol')))


def _thermal(spc_a, spc_b):
    """A <=> B, ordinary reversible thermal reaction."""
    return Reaction(
        reactants=[spc_a], products=[spc_b], reversible=True,
        kinetics=Arrhenius(A=(100.0, 's^-1'), n=0.0, Ea=(10.0, 'kJ/mol')))


def _mole_fractions(electron, ar, ar_ion, spc_a=None, spc_b=None, y_e=Y_E0):
    imf = {electron: y_e, ar: Y_AR0, ar_ion: Y_ION0}
    if spc_a is not None:
        imf[spc_a] = Y_A0
    if spc_b is not None:
        imf[spc_b] = Y_B0
    return imf


def _reactor(imf, T=T_GAS, Te=T_E, P=P0):
    return PlasmaReactor(T, P, imf, (Te, 'K') if not isinstance(Te, tuple) else Te,
                         n_sims=1, termination=[])


def _initialize(reactor, core_species, core_reactions):
    reactor.initialize_model(core_species, core_reactions, [], [])
    return reactor


def _full_system():
    """The standard fixture: electron + cation + Ar + thermal pair, all channels."""
    electron, ar, ar_ion, spc_a, spc_b = _species()
    core_species = [ar, electron, ar_ion, spc_a, spc_b]  # electron NOT first, on purpose
    core_reactions = [_ionization(electron, ar, ar_ion),
                      _recombination(electron, ar, ar_ion),
                      _thermal(spc_a, spc_b),
                      # duplicated electron reactant: exercises the
                      # repeated-species occurrence derivative in the Jacobian
                      _three_body_recombination(electron, ar, ar_ion)]
    imf = _mole_fractions(electron, ar, ar_ion, spc_a, spc_b)
    reactor = _reactor(imf)
    _initialize(reactor, core_species, core_reactions)
    return reactor, core_species, core_reactions


def _two_temp_volume(n_heavy, n_e, T=T_GAS, Te=T_E, P=P0):
    return constants.R * (n_heavy * T + n_e * Te) / P


class PlasmaReactorStateTest:
    """Ruling 2: the electron state and the two-temperature EOS."""

    def test_single_electron_state_with_positive_amount(self):
        """Post-fix test 1: exactly one electron state with N_e > 0."""
        reactor, core_species, _ = _full_system()
        electron_species = [spc for spc in core_species if spc.is_electron()]
        assert len(electron_species) == 1
        assert reactor.electron_index == core_species.index(electron_species[0]) == 1
        assert reactor.electron_species is electron_species[0]
        assert reactor.y0[reactor.electron_index] == Y_E0 > 0.0

    def test_two_temperature_eos_active(self):
        """Post-fix test 2: Te != Tgas activates the two-temperature EOS."""
        reactor, _, _ = _full_system()
        n_total = np.sum(reactor.y0[:reactor.num_core_species])
        v_one_temp = constants.R * T_GAS * n_total / P0
        v_two_temp = _two_temp_volume(n_total - Y_E0, Y_E0)
        assert abs(reactor.V - v_two_temp) / v_two_temp < 1e-12
        assert abs(reactor.V - v_one_temp) / v_one_temp > 1e-6  # genuinely two-temperature

    def test_volume_matches_settled_eos(self):
        """Post-fix test 3: V = R (N_heavy*Tgas + N_e*Te) / P."""
        reactor, _, _ = _full_system()
        n_total = np.sum(reactor.y0[:reactor.num_core_species])
        n_e = reactor.y0[reactor.electron_index]
        assert n_e == Y_E0
        expected = _two_temp_volume(n_total - n_e, n_e)
        assert abs(reactor.V - expected) / expected < 1e-12

    def test_changing_ne_changes_eos_predictably(self):
        """Post-fix test 4: changing N_e changes the EOS result predictably."""
        volumes = {}
        for y_e in (Y_E0, 5.0 * Y_E0):
            electron, ar, ar_ion, spc_a, spc_b = _species()
            reactor = _reactor(_mole_fractions(electron, ar, ar_ion, spc_a, spc_b, y_e=y_e))
            _initialize(reactor, [ar, electron, ar_ion, spc_a, spc_b], [])
            volumes[y_e] = reactor.V
        expected_delta = constants.R * T_E * (5.0 * Y_E0 - Y_E0) / P0
        actual_delta = volumes[5.0 * Y_E0] - volumes[Y_E0]
        assert abs(actual_delta - expected_delta) / expected_delta < 1e-10

    def test_changing_te_changes_only_electron_contribution(self):
        """Post-fix test 5: Te changes only the electron-pressure contribution."""
        volumes = {}
        for te in (T_E, 2.0 * T_E):
            electron, ar, ar_ion, spc_a, spc_b = _species()
            reactor = _reactor(_mole_fractions(electron, ar, ar_ion, spc_a, spc_b), Te=te)
            _initialize(reactor, [ar, electron, ar_ion, spc_a, spc_b], [])
            volumes[te] = reactor.V
        expected_delta = constants.R * Y_E0 * T_E / P0  # R*N_e*(2Te - Te)/P
        actual_delta = volumes[2.0 * T_E] - volumes[T_E]
        assert abs(actual_delta - expected_delta) / expected_delta < 1e-10

    def test_te_equal_tgas_reduces_continuously(self):
        """Post-fix test 6: Te = Tgas with N_e > 0 reduces to P V = R N_total Tgas."""
        electron, ar, ar_ion, spc_a, spc_b = _species()
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion, spc_a, spc_b), Te=T_GAS)
        _initialize(reactor, [ar, electron, ar_ion, spc_a, spc_b], [])
        n_total = np.sum(reactor.y0[:reactor.num_core_species])
        assert reactor.y0[reactor.electron_index] > 0.0  # the electron path is exercised
        v_one_temp = constants.R * T_GAS * n_total / P0
        assert abs(reactor.V - v_one_temp) / v_one_temp < 1e-12
        # Continuity: Te -> Tgas + eps stays within O(eps)
        electron, ar, ar_ion, spc_a, spc_b = _species()
        eps = 1.0e-3
        reactor_eps = _reactor(_mole_fractions(electron, ar, ar_ion, spc_a, spc_b), Te=T_GAS + eps)
        _initialize(reactor_eps, [ar, electron, ar_ion, spc_a, spc_b], [])
        assert abs(reactor_eps.V - v_one_temp) <= 2.0 * constants.R * Y_E0 * eps / P0

    def test_missing_or_zero_electrons_named_failure(self):
        """Post-fix test 7: missing/zero electrons fail loudly, no 1-T fallback."""
        # (a) no electron species at all
        _, ar, ar_ion, spc_a, spc_b = _species()
        reactor = _reactor({ar: Y_AR0, ar_ion: Y_ION0})
        with pytest.raises(PlasmaStateError, match='no electron pseudo-species'):
            reactor.initialize_model([ar, ar_ion], [], [], [])
        # (b) electron species present but zero amount
        electron, ar, ar_ion, _, _ = _species()
        reactor = _reactor({electron: 0.0, ar: Y_AR0, ar_ion: Y_ION0})
        with pytest.raises(PlasmaStateError, match='strictly positive'):
            reactor.initialize_model([ar, electron, ar_ion], [], [], [])
        # (c) electron species present but no amount supplied
        electron, ar, ar_ion, _, _ = _species()
        reactor = _reactor({ar: Y_AR0, ar_ion: Y_ION0})
        with pytest.raises(PlasmaStateError, match='no initial electron amount'):
            reactor.initialize_model([ar, electron, ar_ion], [], [], [])
        # (d) two electron pseudo-species
        electron, ar, ar_ion, _, _ = _species()
        electron2 = Species(label='e2').from_adjacency_list('1 e u1 p0 c-1')
        reactor = _reactor({electron: Y_E0, ar: Y_AR0, ar_ion: Y_ION0})
        with pytest.raises(PlasmaStateError, match='2 electron pseudo-species'):
            reactor.initialize_model([ar, electron, electron2, ar_ion], [], [], [])

    def test_missing_te_named_failure(self):
        """Post-fix test 8: missing/invalid Te fails loudly, no Tgas fallback."""
        electron, ar, ar_ion, _, _ = _species()
        imf = {electron: Y_E0, ar: Y_AR0, ar_ion: Y_ION0}
        with pytest.raises(PlasmaStateError, match='explicit electron temperature'):
            PlasmaReactor(T_GAS, P0, imf, None, n_sims=1, termination=[])
        with pytest.raises(PlasmaStateError, match='strictly positive'):
            PlasmaReactor(T_GAS, P0, imf, (0.0, 'K'), n_sims=1, termination=[])
        with pytest.raises(PlasmaStateError, match='strictly positive'):
            PlasmaReactor(T_GAS, P0, imf, (float('nan'), 'K'), n_sims=1, termination=[])

    def test_no_electron_density_api(self):
        """Ruling 2: no electron_density constructor argument exists."""
        electron, ar, ar_ion, _, _ = _species()
        imf = {electron: Y_E0, ar: Y_AR0, ar_ion: Y_ION0}
        with pytest.raises(TypeError):
            PlasmaReactor(T_GAS, P0, imf, (T_E, 'K'),
                          electron_density=(1e16, 'm^-3'), n_sims=1, termination=[])

    def test_state_packing_preserves_electron(self):
        """Post-fix test 9: packing/unpacking preserve electron amount and index."""
        electron, ar, ar_ion, spc_a, spc_b = _species()
        core_species = [spc_a, ar, electron, ar_ion, spc_b]  # electron at index 2
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion, spc_a, spc_b))
        _initialize(reactor, core_species, [_thermal(spc_a, spc_b)])
        assert reactor.electron_index == 2
        assert reactor.y0[2] == Y_E0
        # Integrate: with no electron reactions the electron slot must be preserved
        reactor.advance(1.0e-3)
        assert abs(reactor.y[2] - Y_E0) <= 1e-12 * Y_E0 + 1e-18
        # Heavy chemistry did proceed (the step was not a no-op)
        assert abs(reactor.y[0] - Y_A0) > 1e-8

    def test_jacobian_matches_finite_differences(self):
        """Post-fix test 10: analytic Jacobian vs finite differences, incl. electron column."""
        reactor, _, _ = _full_system()
        ncs = reactor.num_core_species
        y0 = np.array(reactor.y0[:ncs], float)
        zeros = np.zeros(ncs, float)

        jac = reactor.jacobian(0.0, y0.copy(), zeros.copy(), 0.0)

        fd = np.zeros((ncs, ncs), float)
        for z in range(ncs):
            h = max(1.0e-7 * abs(y0[z]), 1.0e-9)
            yp = y0.copy(); yp[z] += h
            ym = y0.copy(); ym[z] -= h
            rp, _ = reactor.residual(0.0, yp, zeros.copy())
            rm, _ = reactor.residual(0.0, ym, zeros.copy())
            fd[:, z] = (rp - rm) / (2.0 * h)

        scale = np.abs(fd).max()
        assert np.allclose(jac, fd, rtol=1e-4, atol=1e-8 * scale), (
            "analytic Jacobian disagrees with finite differences:\n"
            f"max abs diff = {np.abs(jac - fd).max():.3e}, scale = {scale:.3e}")
        # The electron-state column specifically must agree and be non-trivial
        e_col = reactor.electron_index
        assert np.abs(fd[:, e_col]).max() > 0.0
        assert np.allclose(jac[:, e_col], fd[:, e_col], rtol=1e-4, atol=1e-8 * scale)
        # And the electron-state row (d N_e / dt sensitivity to all species)
        assert np.allclose(jac[e_col, :], fd[e_col, :], rtol=1e-4, atol=1e-8 * scale)

    def test_electron_consuming_reaction_sign(self):
        """Post-fix test 11: electron-consuming reaction moves the electron derivative
        with the correct signed stoichiometry."""
        electron, ar, ar_ion, _, _ = _species()
        core_species = [ar, electron, ar_ion]
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion))
        _initialize(reactor, core_species, [_recombination(electron, ar, ar_ion)])
        ncs = reactor.num_core_species
        y0 = np.array(reactor.y0[:ncs], float)
        zeros = np.zeros(ncs, float)
        res, _ = reactor.residual(0.0, y0.copy(), zeros)
        # e + Ar+ -> Ar: dNe/dt < 0, dNAr+/dt = dNe/dt, dNAr/dt = -dNe/dt
        assert res[1] < 0.0
        assert abs(res[2] - res[1]) <= 1e-12 * abs(res[1])
        assert abs(res[0] + res[1]) <= 1e-12 * abs(res[1])

    def test_incident_electron_rate_scaling(self):
        """Post-fix test 12: one incident electron => first-order scaling in C_e."""
        rates = {}
        volumes = {}
        for y_e in (Y_E0, 2.0 * Y_E0):
            electron, ar, ar_ion, _, _ = _species()
            reactor = _reactor(_mole_fractions(electron, ar, ar_ion, y_e=y_e))
            _initialize(reactor, [ar, electron, ar_ion],
                        [_recombination(electron, ar, ar_ion)])
            ncs = reactor.num_core_species
            y0 = np.array(reactor.y0[:ncs], float)
            res, _ = reactor.residual(0.0, y0.copy(), np.zeros(ncs, float))
            rates[y_e] = -res[1]  # electron consumption rate, mol/s
            volumes[y_e] = reactor.V
        # rate(mol/s) = k * y_e * y_ion / V  =>  ratio = 2 * V1/V2
        expected_ratio = 2.0 * volumes[Y_E0] / volumes[2.0 * Y_E0]
        actual_ratio = rates[2.0 * Y_E0] / rates[Y_E0]
        assert abs(actual_ratio - expected_ratio) / expected_ratio < 1e-10

    def test_charge_and_elemental_balance_over_integration(self):
        """Post-fix test 13: charge and elemental balances invariant over a step."""
        reactor, core_species, _ = _full_system()
        # order: [Ar, e-, Ar+, A, B]; charges [0, -1, +1, 0, 0]; all carry one Ar nucleus but e-
        y_start = np.array(reactor.y0[:reactor.num_core_species], float)
        charge0 = -y_start[1] + y_start[2]
        nuclei0 = y_start[0] + y_start[2] + y_start[3] + y_start[4]
        assert abs(charge0) < 1e-16  # fixture is charge-neutral by construction

        for t in (1.0e-6, 1.0e-5, 1.0e-4):
            reactor.advance(t)
        y_end = np.array(reactor.y[:reactor.num_core_species], float)

        charge_end = -y_end[1] + y_end[2]
        nuclei_end = y_end[0] + y_end[2] + y_end[3] + y_end[4]
        assert abs(charge_end) < 1e-6 * max(y_end[1], y_end[2])
        assert abs(nuclei_end - nuclei0) < 1e-8 * nuclei0
        # The electron population actually evolved: this is a live two-temperature
        # electron state, not a stored attribute.
        assert abs(y_end[1] - y_start[1]) > 1e-12

    def test_residual_and_jacobian_consistency(self):
        """Post-fix test 14: residual and Jacobian share one reaction set, one
        electron resolution, and one reverse-rate policy."""
        reactor, _, core_reactions = _full_system()
        ncs = reactor.num_core_species
        y0 = np.array(reactor.y0[:ncs], float)
        zeros = np.zeros(ncs, float)

        # Same accepted rate-coefficient arrays feed both paths
        kb = np.array(reactor.kb, float)
        j_ion = reactor.reaction_index[core_reactions[0]]
        j_rec = reactor.reaction_index[core_reactions[1]]
        j_thermal = reactor.reaction_index[core_reactions[2]]
        assert kb[j_ion] == 0.0 and kb[j_rec] == 0.0  # no reconstructed reverse rate
        assert kb[j_thermal] > 0.0                     # ordinary thermal reverse intact

        # Same EOS from the same state in both evaluations
        reactor.residual(0.0, y0.copy(), zeros.copy())
        v_residual = reactor.V
        reactor.jacobian(0.0, y0.copy(), zeros.copy(), 0.0)
        v_jacobian = reactor.compute_volume(y0)
        assert v_residual == v_jacobian

        # And the Jacobian is the derivative of the residual (spot FD check)
        jac = reactor.jacobian(0.0, y0.copy(), zeros.copy(), 0.0)
        z = reactor.electron_index
        h = max(1.0e-7 * abs(y0[z]), 1.0e-9)
        yp = y0.copy(); yp[z] += h
        ym = y0.copy(); ym[z] -= h
        rp, _ = reactor.residual(0.0, yp, zeros.copy())
        rm, _ = reactor.residual(0.0, ym, zeros.copy())
        fd_col = (rp - rm) / (2.0 * h)
        assert np.allclose(jac[:, z], fd_col, rtol=1e-4, atol=1e-8 * np.abs(fd_col).max())


class PlasmaReverseRatePolicyTest:
    """Ruling 1: no automatic thermodynamic reversal for Te-dependent kinetics."""

    @staticmethod
    def _reversible_te_mechanism(Te=T_E):
        electron, ar, ar_ion, spc_a, spc_b = _species()
        rxn = Reaction(reactants=[spc_a], products=[spc_b], reversible=True,
                       kinetics=TwoTemperaturePlasma(A=(1.0e3, 's^-1'), n=0.0,
                                                     Ea_g=(4000.0, 'J/mol'),
                                                     Ea_e=(1000.0, 'J/mol')))
        core_species = [ar, electron, ar_ion, spc_a, spc_b]
        imf = _mole_fractions(electron, ar, ar_ion, spc_a, spc_b)
        reactor = _reactor(imf, Te=Te)
        return reactor, core_species, [rxn], rxn

    def test_reversible_te_dependent_rejected_before_integration(self):
        """RED-before target: a reversible Te-dependent reaction must be rejected
        before integration, naming the equation and the kinetics class -- not
        given a finite reverse rate through kf / Keq(Tgas)."""
        reactor, core_species, core_reactions, rxn = self._reversible_te_mechanism()
        with pytest.raises(NonEquilibriumReverseRateError) as excinfo:
            reactor.initialize_model(core_species, core_reactions, [], [])
        message = str(excinfo.value)
        assert 'A <=> B' in message
        assert 'TwoTemperaturePlasma' in message
        assert 'irreversible' in message

    def test_te_equal_tgas_does_not_bypass_guard(self):
        """Positive control 6: Te == Tgas does not bypass the declaration-time guard."""
        reactor, core_species, core_reactions, _ = self._reversible_te_mechanism(Te=T_GAS)
        with pytest.raises(NonEquilibriumReverseRateError):
            reactor.initialize_model(core_species, core_reactions, [], [])

    def test_metadata_only_electron_rejected(self):
        """A nonzero metadata-only electron count cannot express incident-electron
        order and is rejected, not guessed at. Since the electron-representation
        boundary was wired into initialize_model, a metadata-only electron
        reaction is routed to the family-declared resolver, which rejects this
        family-less reaction by name (ElectronPlacementError) -- the single
        representation-error path, not a second parallel one. The reaction is
        never silently accepted."""
        electron, ar, ar_ion, spc_a, spc_b = _species()
        rxn = Reaction(reactants=[spc_a], products=[spc_b], reversible=False,
                       electrons=1,
                       kinetics=Arrhenius(A=(100.0, 's^-1'), n=0.0, Ea=(10.0, 'kJ/mol')))
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion, spc_a, spc_b))
        with pytest.raises(ElectronPlacementError, match='no family attribution'):
            reactor.initialize_model([ar, electron, ar_ion, spc_a, spc_b], [rxn], [], [])

    def test_reversible_thermal_electron_containing_rejected(self):
        """A reversible *thermal* reaction with an explicit electron would price the
        electron's thermochemistry at Tgas; it is rejected as well."""
        electron, ar, ar_ion, spc_a, spc_b = _species()
        rxn = Reaction(reactants=[electron, ar_ion], products=[ar], reversible=True,
                       kinetics=Arrhenius(A=(1.0e6, 'm^3/(mol*s)'), n=0.0, Ea=(0.0, 'kJ/mol')))
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion, spc_a, spc_b))
        with pytest.raises(NonEquilibriumReverseRateError, match='electron-containing'):
            reactor.initialize_model([ar, electron, ar_ion, spc_a, spc_b], [rxn], [], [])

    def test_control_thermal_reverse_unchanged(self):
        """Positive control 1: an ordinary reversible thermal reaction keeps
        kr = kf / Keq(Tgas), numerically identical to SimpleReactor's."""
        reactor, _, core_reactions = _full_system()
        j_thermal = reactor.reaction_index[core_reactions[2]]
        kb_plasma = float(reactor.kb[j_thermal])
        keq_plasma = float(reactor.Keq[j_thermal])

        spc_a, spc_b = core_reactions[2].reactants[0], core_reactions[2].products[0]
        simple = SimpleReactor(T_GAS, P0, {spc_a: Y_A0, spc_b: Y_B0}, n_sims=1, termination=[])
        simple.initialize_model([spc_a, spc_b], [core_reactions[2]], [], [])
        assert kb_plasma == float(simple.kb[0])
        assert keq_plasma == float(simple.Keq[0])

    def test_control_irreversible_electron_impact(self):
        """Positive control 2: an explicitly irreversible electron-impact reaction
        evaluates at (Tgas, Te), integrates, and has exactly zero reverse rate."""
        electron, ar, ar_ion, _, _ = _species()
        rxn = _ionization(electron, ar, ar_ion)
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion))
        _initialize(reactor, [ar, electron, ar_ion], [rxn])
        j = reactor.reaction_index[rxn]
        expected_kf = rxn.kinetics.get_rate_coefficient_two_temp(T_GAS, T_E)
        assert float(reactor.kf[j]) == expected_kf
        # and NOT the one-temperature collapses
        assert float(reactor.kf[j]) != rxn.kinetics.get_rate_coefficient(T_E)   # k(Te, Te)
        assert float(reactor.kf[j]) != rxn.kinetics.get_rate_coefficient(T_GAS)  # k(T, T)
        assert float(reactor.kb[j]) == 0.0
        reactor.advance(1.0e-5)
        y = np.array(reactor.y[:reactor.num_core_species], float)
        assert np.all(np.isfinite(y))
        assert y[1] > Y_E0  # ionization produced electrons

    def test_control_forward_and_reverse_pair(self):
        """Positive control 3: a forward electron process and a separately supplied
        reverse electron process both integrate with their own kinetics."""
        electron, ar, ar_ion, _, _ = _species()
        fwd = _ionization(electron, ar, ar_ion)
        rev = _three_body_recombination(electron, ar, ar_ion)
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion))
        _initialize(reactor, [ar, electron, ar_ion], [fwd, rev])
        j_fwd = reactor.reaction_index[fwd]
        j_rev = reactor.reaction_index[rev]
        assert float(reactor.kf[j_fwd]) == fwd.kinetics.get_rate_coefficient_two_temp(T_GAS, T_E)
        assert float(reactor.kf[j_rev]) == rev.kinetics.get_rate_coefficient_two_temp(T_GAS, T_E)
        assert float(reactor.kb[j_fwd]) == 0.0
        assert float(reactor.kb[j_rev]) == 0.0
        reactor.advance(1.0e-5)
        assert np.all(np.isfinite(np.array(reactor.y[:reactor.num_core_species], float)))

    def test_control_te_changes_electron_rates_only(self):
        """Positive control 4: changing Te changes electron rates, not the
        ordinary thermal reverse-rate construction."""
        results = {}
        for te in (T_E, 2.0 * T_E):
            electron, ar, ar_ion, spc_a, spc_b = _species()
            ion = _ionization(electron, ar, ar_ion)
            thermal = _thermal(spc_a, spc_b)
            reactor = _reactor(_mole_fractions(electron, ar, ar_ion, spc_a, spc_b), Te=te)
            _initialize(reactor, [ar, electron, ar_ion, spc_a, spc_b], [ion, thermal])
            results[te] = (float(reactor.kf[reactor.reaction_index[ion]]),
                           float(reactor.kf[reactor.reaction_index[thermal]]),
                           float(reactor.kb[reactor.reaction_index[thermal]]),
                           float(reactor.Keq[reactor.reaction_index[thermal]]))
        assert results[T_E][0] != results[2.0 * T_E][0]        # electron rate moved
        assert results[T_E][1] == results[2.0 * T_E][1]        # thermal kf unchanged
        assert results[T_E][2] == results[2.0 * T_E][2]        # thermal kr unchanged
        assert results[T_E][3] == results[2.0 * T_E][3]        # Keq(Tgas) unchanged

    def test_control_tgas_changes_thermal_no_implicit_reverse(self):
        """Positive control 5: changing Tgas changes thermal thermochemistry but
        creates no implicit reverse electron rate."""
        results = {}
        for t_gas in (T_GAS, 1.2 * T_GAS):
            electron, ar, ar_ion, spc_a, spc_b = _species()
            ion = _ionization(electron, ar, ar_ion)
            thermal = _thermal(spc_a, spc_b)
            reactor = _reactor(_mole_fractions(electron, ar, ar_ion, spc_a, spc_b), T=t_gas)
            _initialize(reactor, [ar, electron, ar_ion, spc_a, spc_b], [ion, thermal])
            results[t_gas] = (float(reactor.Keq[reactor.reaction_index[thermal]]),
                              float(reactor.kb[reactor.reaction_index[ion]]))
        assert results[T_GAS][0] != results[1.2 * T_GAS][0]  # thermal Keq moved
        assert results[T_GAS][1] == 0.0
        assert results[1.2 * T_GAS][1] == 0.0                # still no reverse electron rate

    def test_control_no_reverse_flux_through_any_path(self):
        """Positive control 7: residual and Jacobian never reconstruct an implicit
        reverse rate for the irreversible electron reaction."""
        electron, ar, ar_ion, _, _ = _species()
        rxn = _ionization(electron, ar, ar_ion)
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion))
        _initialize(reactor, [ar, electron, ar_ion], [rxn])
        ncs = reactor.num_core_species
        y0 = np.array(reactor.y0[:ncs], float)
        zeros = np.zeros(ncs, float)
        reactor.residual(0.0, y0.copy(), zeros.copy())
        # Production of a pure reactant (Ar) could only come from a reverse flux
        assert float(reactor.core_species_production_rates[0]) == 0.0
        # Jacobian: no dependence of any rate on a pure product's amount except
        # through the EOS volume term; verified against finite differences
        jac = reactor.jacobian(0.0, y0.copy(), zeros.copy(), 0.0)
        z = 2  # Ar+ column (pure product)
        h = max(1.0e-7 * abs(y0[z]), 1.0e-9)
        yp = y0.copy(); yp[z] += h
        ym = y0.copy(); ym[z] -= h
        rp, _ = reactor.residual(0.0, yp, zeros.copy())
        rm, _ = reactor.residual(0.0, ym, zeros.copy())
        fd_col = (rp - rm) / (2.0 * h)
        assert np.allclose(jac[:, z], fd_col, rtol=1e-4, atol=1e-10 * max(1.0, np.abs(fd_col).max()))

    def test_control_export_reload_preserves_arrows(self):
        """Positive control 8: serialization round-trips preserve the irreversible
        arrows of the electron pair; nothing merges them into one reversible
        reaction."""
        electron, ar, ar_ion, _, _ = _species()
        fwd = _ionization(electron, ar, ar_ion)
        rev = _three_body_recombination(electron, ar, ar_ion)
        restored_fwd, restored_rev = pickle.loads(pickle.dumps((fwd, rev)))
        assert restored_fwd.reversible is False
        assert restored_rev.reversible is False
        assert isinstance(restored_fwd.kinetics, TwoTemperaturePlasma)
        assert isinstance(restored_rev.kinetics, TwoTemperaturePlasma)
        assert len(restored_fwd.reactants) == 2 and len(restored_fwd.products) == 3
        assert len(restored_rev.reactants) == 3 and len(restored_rev.products) == 2
        # The reactor's own pickle contract carries Te and the mole fractions
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion))
        cls, args = reactor.__reduce__()[0], reactor.__reduce__()[1]
        assert cls is PlasmaReactor
        assert float(args[3].value_si) == pytest.approx(T_E)


class PlasmaReactorUnsupportedTest:
    """Every unsupported configuration fails explicitly."""

    def test_sensitivity_unsupported(self):
        electron, ar, ar_ion, _, _ = _species()
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion))
        with pytest.raises(PlasmaStateError, match='sensitivity'):
            reactor.initialize_model([ar, electron, ar_ion], [], [], [], sensitivity=True)

    def test_ranged_conditions_unsupported(self):
        electron, ar, ar_ion, _, _ = _species()
        imf = {electron: Y_E0, ar: Y_AR0, ar_ion: Y_ION0}
        with pytest.raises(PlasmaStateError, match='scalar'):
            PlasmaReactor([(1000, 'K'), (1500, 'K')], P0, imf, (T_E, 'K'),
                          n_sims=1, termination=[])
        with pytest.raises(PlasmaStateError, match='scalar'):
            PlasmaReactor(T_GAS, P0, imf, [(T_E, 'K'), (2 * T_E, 'K')],
                          n_sims=1, termination=[])

    def test_pdep_kinetics_unsupported(self):
        from rmgpy.kinetics import PDepArrhenius
        electron, ar, ar_ion, spc_a, spc_b = _species()
        rxn = Reaction(reactants=[spc_a], products=[spc_b], reversible=False,
                       kinetics=PDepArrhenius(
                           pressures=([0.1, 10.0], 'bar'),
                           arrhenius=[Arrhenius(A=(100.0, 's^-1'), n=0.0, Ea=(10.0, 'kJ/mol')),
                                      Arrhenius(A=(1000.0, 's^-1'), n=0.0, Ea=(10.0, 'kJ/mol'))]))
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion, spc_a, spc_b))
        with pytest.raises(PlasmaStateError, match='pressure-dependent'):
            reactor.initialize_model([ar, electron, ar_ion, spc_a, spc_b], [rxn], [], [])

    def test_thirdbody_kinetics_unsupported(self):
        # The exact wall the 5-torr argon first-light deck hits: PlasmaAir's
        # N2 <=> N + N entry carries ThirdBody kinetics. ThirdBody is a
        # PDepKineticsModel, so is_pressure_dependent() is True and the same
        # plasma.pyx guard that refuses PDepArrhenius refuses it too. A
        # two-temperature plasma has no collider-efficiency model to give a
        # k0 * [M] rate law meaning, so the refusal is deliberate; this test
        # pins it against future weakening. Two sibling ThirdBody entries in
        # the same library, OH + e- <=> OH- and O2 + e- <=> O2-, are refused
        # identically.
        from rmgpy.kinetics.falloff import ThirdBody
        electron, ar, ar_ion, spc_a, spc_b = _species()
        rxn = Reaction(reactants=[spc_a], products=[spc_b], reversible=False,
                       kinetics=ThirdBody(
                           arrheniusLow=Arrhenius(A=(7.0e21, 'cm^3/(mol*s)'), n=-1.6,
                                                  Ea=(225000, 'cal/mol'), T0=(1, 'K'))))
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion, spc_a, spc_b))
        with pytest.raises(PlasmaStateError, match='pressure-dependent'):
            reactor.initialize_model([ar, electron, ar_ion, spc_a, spc_b], [rxn], [], [])

    def test_const_species_unsupported(self):
        electron, ar, ar_ion, _, _ = _species()
        with pytest.raises(PlasmaStateError, match='constant species'):
            PlasmaReactor(T_GAS, P0, _mole_fractions(electron, ar, ar_ion), (T_E, 'K'),
                          n_sims=1, termination=[], const_spc_names=['Ar'])

    def test_core_reaction_with_edge_species_rejected(self):
        electron, ar, ar_ion, spc_a, spc_b = _species()
        rxn = _thermal(spc_a, spc_b)  # A <=> B with B only in the edge
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion, spc_a))
        with pytest.raises(PlasmaStateError, match='non-core species'):
            reactor.initialize_model([ar, electron, ar_ion, spc_a], [rxn], [spc_b], [])

    def test_reaction_arity_out_of_range_rejected(self):
        electron, ar, ar_ion, spc_a, spc_b = _species()
        # empty product side
        rxn = Reaction(reactants=[spc_a], products=[], reversible=False,
                       kinetics=Arrhenius(A=(1.0, 's^-1'), n=0.0, Ea=(0.0, 'kJ/mol')))
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion, spc_a, spc_b))
        with pytest.raises(PlasmaStateError, match='between'):
            reactor.initialize_model([ar, electron, ar_ion, spc_a, spc_b], [rxn], [], [])
        # four products: must be the named error, not a base-class IndexError
        electron, ar, ar_ion, spc_a, spc_b = _species()
        rxn = Reaction(reactants=[spc_a], products=[electron, electron, electron, spc_b],
                       reversible=False,
                       kinetics=Arrhenius(A=(1.0, 's^-1'), n=0.0, Ea=(0.0, 'kJ/mol')))
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion, spc_a, spc_b))
        with pytest.raises(PlasmaStateError, match='between'):
            reactor.initialize_model([ar, electron, ar_ion, spc_a, spc_b], [rxn], [], [])

    def test_reaction_with_unknown_species_rejected(self):
        electron, ar, ar_ion, spc_a, spc_b = _species()
        stranger = Species(label='X').from_adjacency_list('1 Ar u0 p4 c0')
        rxn = Reaction(reactants=[spc_a], products=[stranger], reversible=False,
                       kinetics=Arrhenius(A=(1.0, 's^-1'), n=0.0, Ea=(0.0, 'kJ/mol')))
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion, spc_a, spc_b))
        with pytest.raises(PlasmaStateError, match='not in the model'):
            reactor.initialize_model([ar, electron, ar_ion, spc_a, spc_b], [rxn], [], [])

    def test_base_class_initialization_bypass_rejected(self):
        """Initializing through the base class and then calling
        generate_rate_coefficients directly must fail loudly, not silently skip
        the plasma validation."""
        from rmgpy.solver.base import ReactionSystem
        electron, ar, ar_ion, spc_a, spc_b = _species()
        rxn = _thermal(spc_a, spc_b)
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion, spc_a, spc_b))
        ReactionSystem.initialize_model(reactor,
                                        core_species=[ar, electron, ar_ion, spc_a, spc_b],
                                        core_reactions=[rxn], edge_species=[], edge_reactions=[],
                                        surface_species=[], surface_reactions=[])
        with pytest.raises(PlasmaStateError, match='validated'):
            reactor.generate_rate_coefficients([rxn], [])

    def test_electron_density_kinetics_requires_explicit_electron(self):
        """A rate law that declares uses_electron_density must carry an explicit
        incident electron; otherwise the rate silently loses its N_e dependence."""
        class DeclaredNeKinetics(Arrhenius):
            pass

        electron, ar, ar_ion, spc_a, spc_b = _species()
        kin = DeclaredNeKinetics(A=(1.0, 's^-1'), n=0.0, Ea=(0.0, 'kJ/mol'))
        kin.uses_electron_density = True
        rxn = Reaction(reactants=[spc_a], products=[spc_b], reversible=False, kinetics=kin)
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion, spc_a, spc_b))
        with pytest.raises(PlasmaStateError, match='uses_electron_density'):
            reactor.initialize_model([ar, electron, ar_ion, spc_a, spc_b], [rxn], [], [])

    def test_foreign_initial_composition_key_rejected(self):
        electron, ar, ar_ion, _, _ = _species()
        stranger = Species(label='X').from_adjacency_list('1 Ar u0 p4 c0')
        imf = _mole_fractions(electron, ar, ar_ion)
        imf[stranger] = 0.1
        reactor = _reactor(imf)
        with pytest.raises(PlasmaStateError, match='not a core species'):
            reactor.initialize_model([ar, electron, ar_ion], [], [], [])

    def test_te_kinetics_without_evaluator_unsupported(self):
        """A kinetics class that declares uses_electron_temperature but exposes
        no explicit two-temperature evaluator must be rejected, never routed
        through the standard k(T, P) interface."""
        class DeclaredTeOnlyKinetics(Arrhenius):
            pass

        electron, ar, ar_ion, spc_a, spc_b = _species()
        kin = DeclaredTeOnlyKinetics(A=(100.0, 's^-1'), n=0.0, Ea=(10.0, 'kJ/mol'))
        kin.uses_electron_temperature = True
        rxn = Reaction(reactants=[spc_a], products=[spc_b], reversible=False, kinetics=kin)
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion, spc_a, spc_b))
        with pytest.raises(PlasmaStateError, match='neither'):
            reactor.initialize_model([ar, electron, ar_ion, spc_a, spc_b], [rxn], [], [])

    def test_external_conditions_unsupported(self):
        electron, ar, ar_ion, _, _ = _species()
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion))
        with pytest.raises(PlasmaStateError, match='conditions'):
            reactor.initialize_model([ar, electron, ar_ion], [], [], [],
                                     conditions={'Te': (2 * T_E, 'K')})


@pytest.mark.database
class PlasmaDatabaseSmokeTest:
    """A minimal database-loaded plasma mechanism against the pinned database:
    H ionization (Voronov) and H+ radiative recombination (Badnell), both
    loaded from the database's YAML tables, as two separately named
    irreversible electron reactions -- the representation Ruling 1 prefers."""

    def test_database_loaded_hydrogen_plasma_smoke(self):
        import os
        from rmgpy import settings
        from rmgpy.data.thermo import ThermoDatabase
        from rmgpy.kinetics.arrhenius import BadnellRRArrhenius, VoronovEIArrhenius

        electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
        h_atom = Species(label='H').from_adjacency_list('1 H u1 p0 c0')
        h_ion = Species(label='H+').from_adjacency_list('1 H u0 p0 c+1')

        # Thermo from the pinned database (proves the database join; the two
        # irreversible reactions never consult thermo for a reverse rate)
        thermo_db = ThermoDatabase()
        thermo_db.load(os.path.join(settings['database.directory'], 'thermo'),
                       libraries=['primaryThermoLibrary'], depository=False)
        h_atom.thermo = thermo_db.get_thermo_data(h_atom)
        assert h_atom.thermo is not None

        # Kinetics from the pinned database's YAML tables
        ionization = Reaction(
            reactants=[electron, h_atom], products=[h_ion, electron, electron],
            reversible=False, kinetics=VoronovEIArrhenius(Z=1, N=1))
        recombination = Reaction(
            reactants=[electron, h_ion], products=[h_atom],
            reversible=False, kinetics=BadnellRRArrhenius(Z=1, N=0))
        assert ionization.kinetics.uses_electron_temperature
        assert recombination.kinetics.uses_electron_temperature

        core_species = [h_atom, electron, h_ion]
        reactor = _reactor({electron: Y_E0, h_atom: 1.0, h_ion: Y_ION0})
        _initialize(reactor, core_species, [ionization, recombination])

        j_ion = reactor.reaction_index[ionization]
        j_rec = reactor.reaction_index[recombination]
        assert float(reactor.kf[j_ion]) == ionization.kinetics.get_rate_coefficient_electron_temp(T_E)
        assert float(reactor.kf[j_rec]) == recombination.kinetics.get_rate_coefficient_electron_temp(T_E)
        assert float(reactor.kb[j_ion]) == 0.0
        assert float(reactor.kb[j_rec]) == 0.0

        y_start = np.array(reactor.y0[:3], float)
        for t in (1.0e-6, 1.0e-5, 1.0e-4):
            reactor.advance(t)
        y_end = np.array(reactor.y[:3], float)

        assert np.all(np.isfinite(y_end))
        assert y_end[1] > y_start[1]          # net ionization grows the electrons
        assert y_end[2] > y_start[2]          # and the cations
        charge = -y_end[1] + y_end[2]
        assert abs(charge) < 1e-6 * max(y_end[1], y_end[2])
        nuclei = y_end[0] + y_end[2]
        assert abs(nuclei - (y_start[0] + y_start[2])) < 1e-8


class PlasmaReactorDriverInterfaceTest:
    """Wall B: RMG.execute() computes Tmin/Tmax over every reaction system with
    ``x.Trange[0] if x.Trange else x.T`` (main.py:879). Every other concrete
    reactor (simple/liquid/surface) declares ``Trange`` as a public attribute
    that is None for scalar-T conditions; PlasmaReactor was the only one that
    omitted it, so the driver raised AttributeError before any chemistry ran."""

    def test_plasma_reactor_exposes_trange_as_none(self):
        """A scalar plasma reactor carries Trange == None, identical to the
        contract a scalar SimpleReactor honors -- the honest value the driver
        reads to fall back to the scalar gas temperature."""
        electron, ar, ar_ion, _, _ = _species()
        reactor = _reactor(_mole_fractions(electron, ar, ar_ion))
        # Before the fix this attribute did not exist -> AttributeError.
        assert reactor.Trange is None
        # PlasmaReactor forbids ranged conditions, so Trange is permanently None.
        assert reactor.T.value_si == T_GAS

    def test_driver_tmin_tmax_expression_reduces_to_scalar_temperature(self):
        """The exact list comprehension from RMG.execute() must run over a
        PlasmaReactor and return its scalar gas temperature, mixing cleanly with
        a ranged SimpleReactor in the same reaction_systems list."""
        electron, ar, ar_ion, spc_a, spc_b = _species()
        plasma = _reactor(_mole_fractions(electron, ar, ar_ion))
        ranged = SimpleReactor(T=[(800, 'K'), (1200, 'K')], P=(P0, 'Pa'),
                               initial_mole_fractions={spc_a: 1.0}, termination=[])
        systems = [ranged, plasma]
        # Verbatim from rmgpy/rmg/main.py:879-880; raised AttributeError pre-fix.
        tmin = min([x.Trange[0].value_si if x.Trange else x.T.value_si for x in systems])
        tmax = max([x.Trange[1].value_si if x.Trange else x.T.value_si for x in systems])
        assert tmin == 800.0            # from the ranged SimpleReactor lower bound
        assert tmax == 1200.0           # from the ranged SimpleReactor upper bound
        # And the plasma reactor alone reduces to its scalar gas temperature.
        assert min([x.Trange[0].value_si if x.Trange else x.T.value_si for x in [plasma]]) == T_GAS
        assert max([x.Trange[1].value_si if x.Trange else x.T.value_si for x in [plasma]]) == T_GAS


@pytest.mark.database
class PlasmaElectronThermoTest:
    """Wall A: the electron pseudo-species declared u1 must resolve to the
    canonical electron thermo entry. is_electron() is deliberately
    multiplicity-agnostic (an electron is an electron whether u0 or u1), but the
    structure-keyed thermo library lookup distinguished u0 from u1, so a u1
    electron missed the u0 database entry, fell through to group additivity, and
    crashed inside RDKit ('Element e not found'). The electron thermo lookup
    must be multiplicity-agnostic and resolve to the correct entry."""

    def test_u1_electron_resolves_to_canonical_electron_thermo(self):
        import os
        from rmgpy import settings
        from rmgpy.data.thermo import ThermoDatabase

        thermo_db = ThermoDatabase()
        thermo_db.load(os.path.join(settings['database.directory'], 'thermo'),
                       libraries=['electrocatThermo'], depository=False)

        electron_u1 = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
        assert electron_u1.is_electron()

        # The library lookup must HIT the (u0) electron entry despite the u1
        # declaration -- not merely avoid raising.
        hit = thermo_db.get_thermo_data_from_libraries(electron_u1)
        assert hit is not None, 'u1 electron missed the electron thermo library entry'
        thermo_data, _library, entry = hit
        assert entry.label == 'electron'
        assert abs(thermo_data.get_enthalpy(298.0)) < 1.0  # canonical electron: H298 = 0 J/mol

        # The full path must return that entry, never fall through to group
        # additivity / RDKit.
        full = thermo_db.get_thermo_data(electron_u1)
        assert abs(full.get_enthalpy(298.0)) < 1.0


class PlasmaChargeNeutralityWarningTest:
    """
    I-169: the fifth check in ``_validate_electron_state`` -- the initial
    composition's net charge.

    Before this check existed, an ``electronDensity``-seeded deck with no compensating
    cation ran to completion of model generation without one message naming the charge
    it carried: measured on the I-164 5 torr argon control deck, whose RMG.log held a
    single ``Warning:`` line, about edge-species saving.

    The check WARNS and never raises. Every test here asserts on that: a non-neutral
    composition must reach a fully initialized reactor, not an exception.
    """

    def _non_neutral_system(self, y_ion):
        """The standard fixture with the cation amount decoupled from the electron."""
        electron, ar, ar_ion, spc_a, spc_b = _species()
        core_species = [ar, electron, ar_ion, spc_a, spc_b]
        core_reactions = [_ionization(electron, ar, ar_ion),
                          _recombination(electron, ar, ar_ion),
                          _thermal(spc_a, spc_b),
                          _three_body_recombination(electron, ar, ar_ion)]
        imf = {electron: Y_E0, ar: Y_AR0, ar_ion: y_ion, spc_a: Y_A0, spc_b: Y_B0}
        return _reactor(imf), core_species, core_reactions

    def test_non_neutral_composition_warns_and_names_the_net_charge(self, caplog):
        """The deck shape the ticket is about: electrons, no compensating cation."""
        reactor, core_species, core_reactions = self._non_neutral_system(0.0)
        with caplog.at_level(logging.WARNING):
            _initialize(reactor, core_species, core_reactions)
        messages = [r.getMessage() for r in caplog.records
                    if 'NOT charge neutral' in r.getMessage()]
        assert len(messages) == 1
        # The message must carry the NUMBER, not merely the fact -- someone reading
        # RMG.log has to learn what the imbalance is without recomputing it.
        total = Y_E0 + Y_AR0 + Y_A0 + Y_B0
        expected = -Y_E0 / total
        assert repr(expected) in messages[0]
        assert 'chargeBalanceSpecies' in messages[0]     # names the way to fix it

    def test_non_neutral_composition_does_not_raise(self):
        """Non-goal made executable: warn, never raise. The reactor initializes."""
        reactor, core_species, core_reactions = self._non_neutral_system(0.0)
        _initialize(reactor, core_species, core_reactions)
        assert reactor.y0[reactor.electron_index] == Y_E0
        assert reactor.V > 0.0

    def test_neutral_composition_is_silent(self, caplog):
        """No false positives on the composition a modeller meant to write."""
        with caplog.at_level(logging.WARNING):
            _full_system()
        assert not [r for r in caplog.records if 'charge neutral' in r.getMessage()]

    def test_excess_cation_also_warns(self, caplog):
        """The check is signed-agnostic: net POSITIVE is reported too."""
        reactor, core_species, core_reactions = self._non_neutral_system(1.0e-3)
        with caplog.at_level(logging.WARNING):
            _initialize(reactor, core_species, core_reactions)
        messages = [r.getMessage() for r in caplog.records
                    if 'NOT charge neutral' in r.getMessage()]
        assert len(messages) == 1
        total = Y_E0 + Y_AR0 + 1.0e-3 + Y_A0 + Y_B0
        expected = (1.0e-3 - Y_E0) / total
        assert expected > 0.0                      # net POSITIVE, not negative
        assert repr(expected) in messages[0]

    def test_imbalance_below_tolerance_is_silent(self, caplog):
        """
        A cation off by one part in 1e9 of the electron amount is roundoff-scale, not a
        modelling error: relative imbalance 5e-10 < RTOL = 1e-6. This is the test that
        stops the check from crying wolf on arithmetic that did close.
        """
        reactor, core_species, core_reactions = self._non_neutral_system(
            Y_E0 * (1.0 + 1.0e-9))
        with caplog.at_level(logging.WARNING):
            _initialize(reactor, core_species, core_reactions)
        assert not [r for r in caplog.records if 'NOT charge neutral' in r.getMessage()]

    def test_weakly_ionized_full_imbalance_is_caught(self, caplog):
        """
        The corner an absolute-only tolerance would miss: an electron amount of 1e-10
        with no cation at all. The net charge is tiny in absolute terms but the
        composition is 100% imbalanced, and RTOL sees it.
        """
        electron, ar, ar_ion, spc_a, spc_b = _species()
        core_species = [ar, electron, ar_ion, spc_a, spc_b]
        core_reactions = [_ionization(electron, ar, ar_ion),
                          _recombination(electron, ar, ar_ion),
                          _thermal(spc_a, spc_b),
                          _three_body_recombination(electron, ar, ar_ion)]
        imf = {electron: 1.0e-10, ar: Y_AR0, ar_ion: 0.0, spc_a: Y_A0, spc_b: Y_B0}
        reactor = _reactor(imf)
        with caplog.at_level(logging.WARNING):
            _initialize(reactor, core_species, core_reactions)
        messages = [r.getMessage() for r in caplog.records
                    if 'NOT charge neutral' in r.getMessage()]
        assert len(messages) == 1
        assert 'relative imbalance 1.0' in messages[0]


class PlasmaChargeBalanceReachabilityTest:
    """
    I-185: a deck can be made charge neutral with an ion nothing can produce.

    ``chargeBalanceSpecies`` assigns the named ion a mole fraction that zeroes the net
    charge; because the composition is then arithmetically neutral,
    ``_warn_if_not_charge_neutral`` correctly stays silent. The defect this restores a
    signal for is that the balancing ion may be a phantom -- produced by no loaded
    reaction -- so the deck *looks* well posed (balanced, no warning) while the
    chemistry that would sustain the balance does not exist. The reactor now checks, at
    model initialization (the first point the reaction set exists), that the ion is
    reachable, and WARNS -- never raises -- when it is not.

    "Reachable" here == the ion's label matches a species that appears as a product of
    some reaction (core or edge), or -- for a reversible reaction -- as a participant on
    either side. Declared-but-produced-by-nothing is the phantom.

    The composition in every fixture is charge neutral by construction (cation amount ==
    electron amount), so ``_warn_if_not_charge_neutral`` is silent and the only charge
    signal that can fire is the reachability one -- which is the whole point: without
    this check, the phantom passes with no message at all.
    """

    def _neutral_imf(self):
        electron, ar, ar_ion, spc_a, spc_b = _species()
        imf = {electron: Y_E0, ar: Y_AR0, ar_ion: Y_E0, spc_a: Y_A0, spc_b: Y_B0}
        return electron, ar, ar_ion, spc_a, spc_b, imf

    def _reactor_with_balance(self, imf, label):
        return PlasmaReactor(T_GAS, P0, imf, (T_E, 'K'), n_sims=1, termination=[],
                             charge_balance_species=label)

    _UNREACHABLE = 'no reaction in the loaded model produces'

    def _reachability_msgs(self, caplog):
        return [r.getMessage() for r in caplog.records
                if self._UNREACHABLE in r.getMessage()]

    def test_defect_a_phantom_balance_is_silent_without_the_directive(self, caplog):
        """
        Verifier 1, the defect: a neutral composition whose cation is produced by NO
        reaction draws no charge signal at all when the reactor is not told the ion is a
        balance ion -- which is exactly the pre-fix behaviour of ``chargeBalanceSpecies``
        (it consumed the label into a mole fraction and told the reactor nothing). The
        neutrality warning is silent (the composition IS neutral); nothing names the
        phantom.
        """
        electron, ar, ar_ion, spc_a, spc_b, imf = self._neutral_imf()
        core_species = [ar, electron, ar_ion, spc_a, spc_b]
        core_reactions = [_thermal(spc_a, spc_b)]        # Ar+ produced by nothing
        reactor = PlasmaReactor(T_GAS, P0, imf, (T_E, 'K'), n_sims=1, termination=[])
        assert reactor.charge_balance_species is None
        with caplog.at_level(logging.WARNING):
            _initialize(reactor, core_species, core_reactions)
        assert not self._reachability_msgs(caplog)
        assert not [r for r in caplog.records if 'NOT charge neutral' in r.getMessage()]

    def test_unreachable_balance_ion_warns_and_names_it(self, caplog):
        """Verifier 2: with the directive, the phantom is diagnosed by name."""
        electron, ar, ar_ion, spc_a, spc_b, imf = self._neutral_imf()
        core_species = [ar, electron, ar_ion, spc_a, spc_b]
        core_reactions = [_thermal(spc_a, spc_b)]        # Ar+ produced by nothing
        reactor = self._reactor_with_balance(imf, 'Ar+')
        with caplog.at_level(logging.WARNING):
            _initialize(reactor, core_species, core_reactions)
        msgs = self._reachability_msgs(caplog)
        assert len(msgs) == 1
        assert "'Ar+'" in msgs[0]                        # names the species
        assert 'chargeBalanceSpecies' in msgs[0]         # says why it is suspect
        assert 'arithmetic' in msgs[0]
        # It is the SOLE charge signal: the neutrality check is (correctly) silent.
        assert not [r for r in caplog.records if 'NOT charge neutral' in r.getMessage()]

    def test_reachable_balance_ion_is_silent(self, caplog):
        """
        Verifier 3: an ion a loaded reaction produces raises nothing. This is the check
        that stops the fix becoming noise.
        """
        electron, ar, ar_ion, spc_a, spc_b, imf = self._neutral_imf()
        core_species = [ar, electron, ar_ion, spc_a, spc_b]
        core_reactions = [_ionization(electron, ar, ar_ion), _thermal(spc_a, spc_b)]
        reactor = self._reactor_with_balance(imf, 'Ar+')
        with caplog.at_level(logging.WARNING):
            _initialize(reactor, core_species, core_reactions)
        assert not self._reachability_msgs(caplog)

    def test_reversible_reaction_credits_the_reactant_side(self, caplog):
        """
        A reversible reaction that has the ion only on its reactant side still counts,
        because its reverse produces the ion; the same reaction irreversible does not.
        Exercised on the check directly to isolate the reachability rule from the
        reverse-rate policy that ``_validate_reactions`` enforces on reversible plasma
        reactions.
        """
        electron, ar, ar_ion, spc_a, spc_b, imf = self._neutral_imf()
        consumes_ion = dict(reactants=[ar_ion, electron], products=[ar],
                            kinetics=Arrhenius(A=(1.0, 's^-1'), n=0.0, Ea=(0.0, 'kJ/mol')))
        reactor = self._reactor_with_balance(imf, 'Ar+')
        with caplog.at_level(logging.WARNING):
            reactor._warn_if_balance_ion_unreachable(
                [Reaction(reversible=True, **consumes_ion)], [])
        assert not self._reachability_msgs(caplog)

        reactor2 = self._reactor_with_balance(imf, 'Ar+')
        with caplog.at_level(logging.WARNING):
            reactor2._warn_if_balance_ion_unreachable(
                [Reaction(reversible=False, **consumes_ion)], [])
        assert len(self._reachability_msgs(caplog)) == 1   # irreversible: not credited

    def test_edge_production_is_credited(self, caplog):
        """An ion produced only by an edge reaction is still loaded chemistry making it."""
        electron, ar, ar_ion, spc_a, spc_b, imf = self._neutral_imf()
        reactor = self._reactor_with_balance(imf, 'Ar+')
        with caplog.at_level(logging.WARNING):
            reactor._warn_if_balance_ion_unreachable([], [_ionization(electron, ar, ar_ion)])
        assert not self._reachability_msgs(caplog)

    def test_unreachable_balance_ion_does_not_raise(self):
        """Non-goal made executable: warn, never raise. The reactor initializes."""
        electron, ar, ar_ion, spc_a, spc_b, imf = self._neutral_imf()
        core_species = [ar, electron, ar_ion, spc_a, spc_b]
        reactor = self._reactor_with_balance(imf, 'Ar+')
        _initialize(reactor, core_species, [_thermal(spc_a, spc_b)])
        assert reactor.y0[reactor.electron_index] == Y_E0
        assert reactor.V > 0.0

    def test_reachability_warning_fires_at_most_once(self, caplog):
        """
        Latched: model generation calls initialize_model many times; the loud signal
        must not repeat every enlargement.
        """
        electron, ar, ar_ion, spc_a, spc_b, imf = self._neutral_imf()
        core_species = [ar, electron, ar_ion, spc_a, spc_b]
        core_reactions = [_thermal(spc_a, spc_b)]
        reactor = self._reactor_with_balance(imf, 'Ar+')
        with caplog.at_level(logging.WARNING):
            _initialize(reactor, core_species, core_reactions)
            _initialize(reactor, core_species, core_reactions)
        assert len(self._reachability_msgs(caplog)) == 1

    def test_charge_balance_species_label_survives_pickle(self):
        """__reduce__ carries the label, so a restored reactor still runs the check."""
        electron, ar, ar_ion, spc_a, spc_b, imf = self._neutral_imf()
        reactor = self._reactor_with_balance(imf, 'Ar+')
        restored = pickle.loads(pickle.dumps(reactor))
        assert restored.charge_balance_species == 'Ar+'
