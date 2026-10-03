#!/usr/bin/env python3
###############################################################################
# RMG - Reaction Mechanism Generator                                          #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),              #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)       #
# Permission is hereby granted, free of charge, to any person obtaining a copy #
# of this software and associated documentation files (the "Software"), to     #
# deal in the Software without restriction, including without limitation the   #
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or   #
# sell copies of the Software, and to permit persons to whom the Software is   #
# furnished to do so, subject to the following conditions:                     #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER          #
# DEALINGS IN THE SOFTWARE.                                                  #
###############################################################################
"""Passing characterisations of the pre-EEDF engine, not migration acceptance tests.

Rate and moment comparisons use archived, identified LoKI-B rows. Small
synthetic decks isolate those legacy behaviours and disclose their Ars
relaxation. The resonance case instead initializes the unchanged Report-10
example and its real database, without added channels or fitted parameters.
"""

from contextlib import ExitStack, contextmanager
import json
import logging
import hashlib
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

# Load the deck's lazy dependencies before any caller-state snapshot.
import rmgpy.electron_placement
import rmgpy.molecule.symmetry
import rmgpy.constants as constants
import rmgpy.data.rmg as rmg_data_module
from rmgpy.data.kinetics import quarantine
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.library import KineticsLibrary, LibraryReaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.thermo import ThermoDatabase
from rmgpy.exceptions import PlasmaStateError, SettingsError
from rmgpy.kinetics import Arrhenius, ElectronCollisionPlasma, TwoTemperaturePlasma
from scipy.integrate import quad
from rmgpy.solver.plasma import PlasmaReactor
import rmgpy.rmg.input as rmg_input
from rmgpy.rmg.main import RMG
from rmgpy.species import Species
from rmgpy.quantity import Quantity
from rmgpy.thermo import ThermoData

DATA = Path(__file__).parent / 'data'
AR_ION = 'e + Ar(1S0) -> e + e + Ar(+,gnd), Ionization'
AR_M1 = 'e + Ar(1S0) -> e + Ar(3P2), Excitation'
AR_M2 = 'e + Ar(1S0) -> e + Ar(3P0), Excitation'
AR_ELASTIC = 'e + Ar(1S0) -> e + Ar(1S0), Elastic'
N2_ION = 'e + N2(X) -> e + e + N2(+,X), Ionization'
TG = 298.15
PRESSURE = 5.0 * 101325.0 / 760.0
REPORT10_THERMO_LIBRARIES = ('primaryThermoLibrary', 'PlasmaThermo', 'electrocatThermo',
                            'PlasmaExcitedNeutralThermo')


def _data():
    return json.loads((DATA / 'loki_eedf_red_rows.json').read_text())


def _resonance_fixture():
    return json.loads((DATA / 'loki_ar_16td_resonance_power.json').read_text())


def _te(row):
    return (2.0 / 3.0) * row['mean_energy_eV'] * constants.e / constants.kB


def _thermo(energy_eV):
    return ThermoData(Tdata=([298, 400, 600, 800, 1000, 1500, 2000], 'K'),
                      Cpdata=([2.5 * constants.R] * 7, 'J/(mol*K)'),
                      H298=(energy_eV * constants.e * constants.Na, 'J/mol'),
                      S298=(154.8, 'J/(mol*K)'))


def _collision(data, key):
    xs = data[key]
    return ElectronCollisionPlasma(energies=(xs['energy_eV'], 'eV/molecule'),
                                   sigma=(xs['sigma_m2'], 'm^2'))


def _ashida():
    return TwoTemperaturePlasma(A=(5e-15, 'm^3/(molecule*s)'), n=.74,
                                Ea_g=(11.56, 'eV/molecule'), Ea_e=(11.56, 'eV/molecule'),
                                T0=(11604.51812, 'K'))


def _constant(rate):
    return TwoTemperaturePlasma(A=(rate, 'm^3/(molecule*s)'), n=0,
                                Ea_g=(0., 'J/mol'), Ea_e=(0., 'J/mol'), T0=(1., 'K'))


class ArchivedRowReactor(PlasmaReactor):
    """Only library-entry resolution is synthetic; the reactor is unchanged."""

    def _declared_entry_kinetics(self, library, index):
        if library != 'RedFixture':
            raise PlasmaStateError('unknown fixture library')
        return self.fixture_entries[index]


@pytest.fixture(autouse=True)
def _isolate_thermo_database():
    saved = rmg_data_module.database
    rmg_data_module.database = None
    try:
        yield
    finally:
        rmg_data_module.database = saved


def _deck(data, row=None, mixture=False, table_rates=False, drift=1., energy=True):
    row = data['ar'] if row is None else row
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    ars = Species(label='Ars').from_adjacency_list('1 Ar u2 p3 c0')
    for spc, h in ((ar, 0.), (arp, 15.76), (ars, 11.56)):
        spc.thermo = _thermo(h)
    core = [electron, ar, arp, ars]
    fractions = {electron:4e-10, arp:4e-10, ar:1.-9e-10, ars:1e-10}
    ion = (_constant(row['rates_m3_s'][AR_ION] * drift) if table_rates else
           _collision(data, 'golyatina'))
    ex = _constant(row['rates_m3_s'][AR_M1]+row['rates_m3_s'][AR_M2]) if table_rates else _ashida()
    entries = {1:ion, 2:ex}
    reactions = [LibraryReaction(reactants=[electron, ar], products=[arp, electron, electron],
                                kinetics=ion, reversible=False, library='RedFixture'),
                 LibraryReaction(reactants=[electron, ar], products=[ars, electron],
                                kinetics=ex, reversible=False, library='RedFixture')]
    energies = {'RedFixture:1':(15.76, 'eV'), 'RedFixture:2':(11.56, 'eV')}
    elastic = {'Ar':{'A':(2.336e-14, 'm^3/s'), 'n':1.609, 'b':.0618, 'c':-.1171},
               'Ars':{'ignore':'trace metastable elastic contribution omitted for this test'}}
    # RedFixture has no thermo library: these generated ThermoData objects
    # preserve the archived synthetic deck, not the database characterisation.
    assertions = {'Ar+': 'ion'}
    if mixture:
        n2 = Species(label='N2').from_smiles('N#N')
        n2p = Species(label='N2+').from_adjacency_list('multiplicity 2\n1 N u1 p0 c+1 {2,T}\n2 N u0 p1 c0 {1,T}')
        n2.thermo, n2p.thermo = _thermo(0.), _thermo(15.6)
        core.extend([n2, n2p])
        fractions[n2], fractions[n2p] = row['x_N2'], 1e-8
        fractions[electron] += 1e-8
        fractions[ar] -= row['x_N2'] + 2e-8
        entries[3] = _constant(row['rates_m3_s'][N2_ION])
        reactions.append(LibraryReaction(reactants=[electron, n2], products=[n2p, electron, electron],
                                        kinetics=entries[3], reversible=False, library='RedFixture'))
        energies['RedFixture:3'] = (15.6, 'eV')
        elastic['N2'] = {'ignore':'composition-chain/mixed-provider proxy; no N2 elastic fit'}
        assertions['N2+'] = 'ion'
    # A fast synthetic Ars relaxation keeps a nontrivial operating point;
    # its rate is not a sourced radiation or resonance model.
    reactions.append(LibraryReaction(reactants=[ars], products=[ar], reversible=False,
                     kinetics=Arrhenius(A=(1e7, 's^-1'), n=0, Ea=(0., 'J/mol'),
                                        T0=(1., 'K')), library='RedFixture'))
    if not energy:
        rec = _constant(1e-13)
        reactions.append(LibraryReaction(reactants=[arp, electron], products=[ar],
                         kinetics=rec, reversible=False, library='RedFixture'))
    kwargs = {'thermo_source_assertions':assertions}
    if energy:
        kwargs.update(diffusion_length=(.02, 'm'), ion_reduced_mobility=(1.535e-4, 'm^2/(V*s)'),
                      wall_single_bath_approximation=mixture,
                      wall_neutralization_products={'Ar+':'Ar', **({'N2+':'N2'} if mixture else {})},
                      electron_energy_balance={'absorbed_power':(.5, 'W'),
                                               'chamber_volume':(np.pi * .05**2 * .30, 'm^3'),
                                               'sheath':'floating_wall', 'electron_energies':energies,
                                               'elastic_collisions':elastic})
    if energy and mixture:
        kwargs.pop('ion_reduced_mobility')
        kwargs.pop('wall_single_bath_approximation')
        # Synthetic equal mobilities isolate electron-energy coupling, as in the wall tests.
        kwargs['ion_reduced_mobilities'] = {
            label: {'perBath': {bath: (1.535e-4, 'm^2/(V*s)') for bath in ('Ar', 'N2')}}
            for label in ('Ar+', 'N2+')}
        kwargs['wall_bath_lumping'] = {'Ars':'Ar'}
    r = ArchivedRowReactor((TG, 'K'), (PRESSURE, 'Pa'), fractions, (_te(row), 'K'), **kwargs)
    r.fixture_entries = entries
    r.initialize_model(core, reactions, [], [])
    return r, core, reactions


def _evaluate(r):
    y = np.array(r.y, float)
    residual = np.array(r.residual(r.t, y, np.zeros_like(y))[0], float)
    return y, residual


def _fit_collision_coefficient(r, te):
    """Read the declared K_m fit itself; never infer it from elastic power."""
    a, n, b, c = r.energy_elastic_params[0]
    ev = te * constants.R / constants.Na / constants.e
    lt = np.log(ev)
    return a * ev**n * np.exp(b * lt**2 + c * lt**3)


def _independent_power(r, core):
    """Extensive balance from reaction fluxes, species, and wall fluxes (W)."""
    y, rhs = _evaluate(r)
    v = r.compute_volume(y)
    te = y[r.te_index]
    ie = r.electron_index
    iar = r.species_index[core[1]]
    inelastic = float(sum(_channel_power(r)))
    mass = core[1].molecular_weight.value_si
    elastic = (3 * constants.m_e / mass * _fit_collision_coefficient(r, te)
               * y[iar] * constants.Na / v * y[ie] * constants.R * (te - TG))
    wall_e = 2 * constants.R * te * (-r.wall_loss_rates[ie])
    wall_i = float(np.dot(r.energy_ion_sheath_factor, -r.wall_loss_rates)) * constants.R * te
    stored = 1.5 * constants.R * (te * rhs[ie] + y[ie] * rhs[r.te_index])
    initial_heavy = sum(r.y0[j] for j in range(r.num_core_species) if j != ie)
    power = .5 * initial_heavy * constants.R * TG / PRESSURE / (np.pi * .05**2 * .30)
    return dict(power=power, inelastic=inelastic, elastic=elastic, wall_e=wall_e,
                wall_i=wall_i, stored=stored,
                closure=power-inelastic-elastic-wall_e-wall_i-stored)


def _run(r, end_time=40.):
    """Check finite accepted states/residuals and late-time stationarity."""
    states = []
    for t in np.logspace(-10, np.log10(end_time), 120):
        r.advance(t)
        y, rhs = _evaluate(r)
        assert np.all(np.isfinite(y))
        assert np.all(np.isfinite(rhs))
        states.append(y)
    scale = np.maximum(np.abs(states[-1]), 1e-14)
    assert np.max(np.abs(states[-1] - states[-2]) / scale) < 1e-4
    assert r.t == pytest.approx(end_time, rel=1e-12, abs=1e-12)
    return r


def _continuous_maxwellian(data, key, te):
    """Independent quadrature on the piecewise-linear sigma, scaled for tiny units."""
    xs = data[key]
    energies = np.array(xs['energy_eV'])
    sigma = np.array(xs['sigma_m2']) / 1e-20
    ev = te * constants.kB / constants.e
    integrals = [quad(lambda e: np.interp(e, energies, sigma) * e * np.exp(-e/ev),
                      a, b, epsabs=1e-13, epsrel=1e-11)[0]
                 for a, b in zip(energies[:-1], energies[1:])]
    return sum(integrals) * 1e-20 * np.sqrt(8 * constants.e/(np.pi*constants.m_e)) / ev**1.5


def _particle_rate(law, te):
    if isinstance(law, ElectronCollisionPlasma):
        return law.get_rate_coefficient_electron_temp(te) / constants.Na
    return law.get_rate_coefficient_two_temp(TG, te) / constants.Na


def test_maxwellian_shape_and_quadrature_are_separate_defects():
    """EEDF §18(a): Maxwellian central rates; superseded by same-row provider rates.

    The same IST cross section isolates EEDF shape. Integrating the identical
    piecewise-linear cross section continuously isolates the engine trapezoid.
    The observed/LoKI ratio is about 209638; only 0.714% is quadrature error.
    """
    data = _data()
    row = data['ar']
    law = _collision(data, 'ist_ionization')
    observed = _particle_rate(law, _te(row))
    assert law.get_rate_coefficient(_te(row))/constants.Na == pytest.approx(observed, rel=1e-12, abs=1e-30)
    energy = np.array(row['energy_eV'])
    sigma = np.interp(energy, data['ist_ionization']['energy_eV'], data['ist_ionization']['sigma_m2'])
    from_eedf = (np.trapz(sigma * energy * np.array(row['f0_eV_minus_3_over_2']), energy)
                 * np.sqrt(2 * constants.e/constants.m_e))
    assert from_eedf == pytest.approx(row['rates_m3_s'][AR_ION], rel=.05, abs=1e-30)
    continuous = _continuous_maxwellian(data, 'ist_ionization', _te(row))
    assert observed / row['rates_m3_s'][AR_ION] == pytest.approx(209637.6365, rel=.002, abs=1.)
    assert observed / continuous == pytest.approx(1.00713670, rel=1e-7, abs=1e-9)
    assert continuous / from_eedf > 1.9e5


@pytest.mark.parametrize('channel, factor', [('golyatina', 153325.3660), ('ashida', 60.16318387)])
def test_active_te_laws_disagree_with_the_archived_row(channel, factor):
    """EEDF §18(a): active Te laws; superseded by provider-owned channel rates.

    Golyatina also changes cross sections; its factor is not shape alone.
    At this Te its trapezoid / continuous-Maxwellian ratio is 1.05449740.
    """
    data = _data()
    law = _collision(data, 'golyatina') if channel == 'golyatina' else _ashida()
    row = data['ar']
    reference = row['rates_m3_s'][AR_ION] if channel == 'golyatina' else sum(row['rates_m3_s'][k] for k in (AR_M1, AR_M2))
    observed = _particle_rate(law, _te(row))
    assert observed == pytest.approx(reference * factor, rel=.002, abs=1e-30)
    if channel == 'golyatina':
        continuous = _continuous_maxwellian(data, 'golyatina', _te(row))
        assert observed/continuous == pytest.approx(1.05449740, rel=1e-7, abs=1e-9)


def test_one_deck_combines_inconsistent_collision_and_power_moments():
    """EEDF §18(c): inconsistent distributions; superseded by one-row moments.

    Compare K_m with LoKI's collision rate and Q_elastic with LoKI's elastic
    energy moment separately. A near power match does not imply a rate match.
    """
    data = _data()
    r, core, _ = _deck(data)
    y, _ = _evaluate(r)
    ne, nar = y[0], y[1]
    rates = [_fit_collision_coefficient(r, _te(data['ar']))]
    rates.extend(r.core_reaction_rates[j] * r.V**2 / (ne*nar*constants.Na) for j in (0, 1))
    row = data['ar']
    refs = [row['rates_m3_s'][AR_ELASTIC], row['rates_m3_s'][AR_ION],
            row['rates_m3_s'][AR_M1] + row['rates_m3_s'][AR_M2]]
    factors = [1.36785899, 153325.3660, 60.16318387]
    for rate, reference, factor in zip(rates, refs, factors):
        assert rate == pytest.approx(reference*factor, rel=.002, abs=1e-30)
    # Obtain the energy moment from the real energy residual, independent of diagnostics.
    budget = _independent_power(r, core)
    qel = budget['power']-budget['inelastic']-budget['wall_e']-budget['wall_i']-budget['stored']
    scale = ne*nar*constants.Na**2/r.V*constants.e
    assert qel/scale == pytest.approx(-row['power_eVm3_s']['Elastic collisions (net)']*.987955737,
                                    rel=.002, abs=1e-30)


def test_molecular_arm_can_use_a_maxwellian_ar_comparator():
    """EEDF §18(b): mixed provider comparison; superseded by same-provider arms.

    This is the accepted architecture at fixed archived states, not an EEDF
    derivative between the pure-Ar usingSDCS and mixed conservative rows.
    """
    data = _data()
    row = data['mixture']
    r, _, mixed_rxns = _deck(data, row=row, mixture=True)
    pure, _, pure_rxns = _deck(data)
    _evaluate(r)
    _evaluate(pure)
    ar = _particle_rate(pure_rxns[0].kinetics, _te(data['ar']))
    n2 = _particle_rate(mixed_rxns[2].kinetics, _te(row))
    assert row['x_N2'] == pytest.approx(.0384615385, rel=1e-10, abs=1e-12)
    assert n2 == pytest.approx(row['rates_m3_s'][N2_ION], rel=1e-12, abs=1e-30)
    assert n2/ar == pytest.approx(3.7066157724e-7, rel=.002, abs=1e-12)
    assert (n2/ar)/(row['rates_m3_s'][N2_ION]/data['ar']['rates_m3_s'][AR_ION]) < 7e-6


def test_drifted_row_finishes_with_finite_residuals():
    """EEDF §18(d): unqualified drift accepted; superseded by provenance refusal.

    No loader exists today: this tests acceptance of doubled numerical rates
    through real row-fed laws. It cannot test a nonexistent manifest loader.
    """
    data = _data()
    row = data['ar']
    r, _, rxns = _deck(data, table_rates=True, drift=2., energy=False)
    observed = _particle_rate(rxns[0].kinetics, _te(row))
    assert observed == pytest.approx(2*row['rates_m3_s'][AR_ION], rel=1e-12, abs=1e-30)
    _run(r)
    y, rhs = _evaluate(r)
    assert np.all(np.isfinite(rhs))
    # With no wall, the prescribed-Te particle equilibrium is k_ion/k_rec.
    assert y[0]/y[1] == pytest.approx(2*row['rates_m3_s'][AR_ION]/1e-13, rel=1e-4, abs=1e-14)


def test_maxwellian_deck_converges_with_a_temperature_coordinate():
    """C1 §8 Maxwellian rates active at convergence; superseded by u-only closure.

    This integrated legacy fixture preserves one of §8's alternative defects;
    it does not invent an E/N input that the current reactor cannot accept.
    """
    data = _data()
    r, core, rxns = _deck(data)
    _run(r)
    y, _ = _evaluate(r)
    budget = _independent_power(r, core)
    assert abs(budget['closure'])/budget['power'] < 5e-6
    assert abs(budget['stored'])/budget['power'] < 1e-4
    assert y[-1] == pytest.approx(10461.4198, rel=.002, abs=1.)
    # Recover the final electron rate from the actual converged reaction flux.
    from_flux = r.core_reaction_rates[0]*r.V**2/(y[0]*y[1]*constants.Na)
    maxwellian = _particle_rate(rxns[0].kinetics, y[-1])
    assert from_flux == pytest.approx(maxwellian, rel=1e-10, abs=1e-30)
    assert maxwellian == pytest.approx(3.54435673e-22, rel=.002, abs=1e-30)


def _unclosed_deck(data):
    """Keep one stationary legacy deck's laws, walls and state; omit only closure."""
    r, core, rxns = _deck(data)
    _run(r)
    y, _ = _evaluate(r)
    declaration = r.electron_energy_balance
    r.initial_mole_fractions = {spc: y[j]/sum(y[:r.num_core_species])
                               for j, spc in enumerate(core)}
    # Te already holds the converged value. Reinitializing the identical
    # particles and walls with just this declaration absent freezes that Te.
    r._configure_energy_balance(None)
    r.initialize_model(core, rxns, [], [])
    return r, core, rxns, declaration


def test_unclosed_electron_deck_still_converges():
    """C1 §8 energy omitted at convergence; superseded by closure refusal.

    Removing only the electron-energy declaration from an operating point
    keeps the walls, all laws and initial composition identical. The fixed-Te
    species-only run remains stationary with positive unaccounted excitation
    and wall energy demand. No flowing-reactor interface exists here.
    """
    r, core, _, _ = _unclosed_deck(_data())
    assert r.neq == r.num_core_species
    _run(r)
    y, rhs = _evaluate(r)
    assert np.all(np.isfinite(rhs))
    demand = float(np.dot(r.core_reaction_rates[:2], [15.76, 11.56])) * r.V * constants.Na * constants.e
    wall_demand = 2 * constants.R * r.Te.value_si * (-r.wall_loss_rates[r.electron_index])
    assert demand > 1.
    assert wall_demand > .01
    assert r.Te.value_si == pytest.approx(10461.4198, rel=.002, abs=1.)


def _require_input(path, directory=False, absent_ok=False):
    """Inspect every component from the filesystem root before descending."""
    import os
    import stat

    path = Path(path)
    absolute = path.absolute()  # Do not resolve away symlinks before classifying.
    for component in list(reversed(absolute.parents)) + [absolute]:
        is_directory = component != absolute or directory
        try:
            info = component.stat()
        except FileNotFoundError:
            if os.path.lexists(component):
                pytest.fail('Report-10 input {0}: dangling symlink'.format(component))
            if absent_ok:
                pytest.skip('Report-10 input is absent: {0}'.format(component))
            pytest.fail('Report-10 input {0}: lacks required content: missing path'.format(component))
        except OSError as exc:
            pytest.fail('Report-10 input {0}: is unreadable: {1}'.format(component, exc))
        correct_type = stat.S_ISDIR(info.st_mode) if is_directory else stat.S_ISREG(info.st_mode)
        if not correct_type:
            pytest.fail('Report-10 input {0}: expected {1}'.format(
                component, 'directory' if is_directory else 'regular file'))
        # Root can otherwise read mode-000 paths. Check read/search mode bits too.
        if not info.st_mode & 0o444 or (is_directory and not info.st_mode & 0o111):
            pytest.fail('Report-10 input {0}: is unreadable'.format(component))
        try:
            if is_directory:
                # Opening the directory tests access without enumerating every
                # sibling at every ancestor on every required-file validation.
                with os.scandir(component) as children:
                    next(children, None)
            else:
                with component.open('rb') as stream:
                    stream.read(1)
        except OSError as exc:
            pytest.fail('Report-10 input {0}: is unreadable: {1}'.format(component, exc))
    return path


def _canonical_species(species):
    """Canonical graph labels, independent of dictionary labels and atom order."""
    from itertools import permutations

    adjacency = []
    for molecule in species.molecule:
        molecule = molecule.copy(deep=True)
        # Exhaustive labelling is exact for these monatomic pinned species and
        # small corrupt variants. Refuse large new structures rather than guess.
        if len(molecule.atoms) > 8:
            raise ValueError('uncharacterised structure with more than eight atoms')
        candidates = []
        for order in permutations(molecule.atoms):
            molecule.vertices = list(order)
            candidates.append(molecule.to_adjacency_list())
        adjacency.append(min(candidates))
    return {'adjacency': sorted(adjacency), 'reactive': species.reactive,
            'props': _canonical_value(species.props)}


def _canonical_value(value):
    """Serialise content, never repr/identity; unknown objects fail closed."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, (float, np.number)):
        return float(value) if np.isfinite(value) else str(value)
    if isinstance(value, Species):
        return _canonical_species(value)
    if hasattr(value, 'value_si') and hasattr(value, 'units'):
        # Canonical SI magnitudes and dimensional units make unit conversions
        # neutral without discarding units or uncertainty information.
        import quantities
        units = str(quantities.Quantity(1, value.units or 'dimensionless').simplified.dimensionality)
        return {'class': type(value).__name__, 'value_si': _canonical_value(value.value_si),
                'units': units, 'uncertainty_si': _canonical_value(value.uncertainty_si),
                'uncertainty_type': value.uncertainty_type}
    if isinstance(value, np.ndarray):
        return _canonical_value(value.tolist())
    if isinstance(value, (tuple, list)):
        return [_canonical_value(item) for item in value]
    if isinstance(value, dict):
        if not all(isinstance(key, str) for key in value):
            raise TypeError('non-string contract mapping key')
        return {key: _canonical_value(item) for key, item in sorted(value.items())}
    return _canonical_object(value, neutral={'comment'})


def _canonical_object(obj, neutral=()):
    """Discover all public data descriptors, including inherited Cython fields."""
    fields = {'class': type(obj).__module__ + '.' + type(obj).__name__}
    names = set(dir(obj)) | set(getattr(obj, '__dict__', {}))
    for name in sorted(names):
        if name in neutral or name.startswith('__'):
            continue
        # Private Cython backing storage is exposed by public descriptors.
        # Python instance fields, including private ones, are also compared.
        if name.startswith('_') and name not in getattr(obj, '__dict__', {}):
            continue
        value = getattr(obj, name)
        if callable(value):
            continue
        try:
            fields[name] = _canonical_value(value)
        except Exception as exc:
            raise TypeError('field {0}: {1}'.format(name, exc)) from exc
    return fields


def _canonical_entry(entry, library_name):
    """Pin the entire entry; the neutral set is deliberately small and explicit.

    Only prose not parsed by the loader or quarantine is neutral; pin parsed
    provenance through the loader's own output even if autoGenerated is false.
    Labels are neutral only where canonical structures replace their spelling.
    For Arrhenius/TwoTemperaturePlasma, finite positive T0 is neutral when n is
    exactly zero. Every other discovered field is pinned.
    """
    fields = _canonical_object(entry, neutral={'item', 'data', 'label', 'short_desc',
                                               'long_desc', 'reference', 'reference_type'})
    reaction = _canonical_object(entry.item, neutral={'label', 'comment', 'reactants', 'products'})
    for side in ('reactants', 'products'):
        reaction[side] = sorted((_canonical_species(spc) for spc in getattr(entry.item, side)),
                                key=lambda item: json.dumps(item, sort_keys=True))
    fields['label'] = {'reactants': reaction['reactants'], 'products': reaction['products'],
                       'reversible': reaction['reversible']}
    fields['reaction'] = reaction
    fields['electron_placement'] = _canonical_value(
        rmgpy.electron_placement.FAMILY_ELECTRON_PLACEMENT.get(library_name))
    fields['kinetics'] = _canonical_object(entry.data, neutral={'comment'})
    # Parsed loader and quarantine provenance is executable, even when it
    # appears in descriptions or comments alongside neutral prose.
    fields['loader_provenance'] = _loader_provenance(entry, library_name)
    fields['authored_families'] = sorted(quarantine._families_from_provenance(
        entry.data.comment, entry.long_desc))
    if type(entry.data) in (Arrhenius, TwoTemperaturePlasma):
        if entry.data.n.value_si == 0 and np.isfinite(entry.data.T0.value_si) and entry.data.T0.value_si > 0:
            fields['kinetics']['T0']['value_si'] = 1.
    return fields


def _loader_provenance(entry, library_name):
    """Pin executable longDesc semantics using the actual provenance reader."""
    library = KineticsLibrary(label=library_name, auto_generated=True)
    library.entries = {entry.index: entry}
    comment = entry.data.comment
    try:
        reaction, = library.get_library_reactions()
        return _canonical_object(reaction, neutral={'entry', 'kinetics', 'reactants',
                                                    'products', 'label', 'comment'})
    finally:
        # The reader strips indentation in comments; preflight must not mutate.
        entry.data.comment = comment


def _canonical_library(library):
    """Discover library metadata with the same field walker as raw entries."""
    # Entries have their own contract. local_context holds constructor bindings,
    # not library input; top is retained as discovered structural metadata.
    return _canonical_object(library, neutral={'entries', 'local_context',
                                               'short_desc', 'long_desc'})


def _differing_fields(expected, actual, prefix=''):
    if isinstance(expected, dict) and isinstance(actual, dict):
        differences = []
        for key in sorted(set(expected) | set(actual)):
            path = prefix + '.' + key if prefix else key
            if key not in expected or key not in actual:
                differences.append(path)
            else:
                differences.extend(_differing_fields(expected[key], actual[key], path))
        return differences
    return [] if expected == actual else [prefix]


def _report10_libraries(path):
    """Load raw entries without ever merging away their inventory/content."""
    context = KineticsDatabase().local_context
    context['nan'] = np.nan
    libraries, problems = {}, []
    for name in ('PlasmaArgon', 'PlasmaRadiativeRecombination'):
        library = KineticsLibrary(label=name)
        library.convert_duplicates_to_multi = lambda: None
        try:
            library.load(str(path / 'kinetics/libraries' / name / 'reactions.py'), local_context=context.copy())
        except Exception as exc:
            problems.append('{0}: cannot load ({1}: {2})'.format(name, type(exc).__name__, exc))
        else:
            libraries[name] = library
    return libraries, problems


@contextmanager
def _report10_thermo(path):
    """Load the job's libraries through its thermo loader, preserving the caller.

    Forward-rate probes need library thermo, not group-additivity estimates.
    The full Report-10 job still loads the complete thermo database itself.
    """
    saved = rmg_data_module.database
    try:
        rmg_data_module.database = None
        database = RMGDatabase()
        database.thermo = ThermoDatabase()
        try:
            database.thermo.load_libraries(str(path / 'thermo/libraries'),
                                          list(REPORT10_THERMO_LIBRARIES))
        except Exception as exc:
            # The job loader adds each library only after a successful load.
            name = next(name for name in REPORT10_THERMO_LIBRARIES
                        if name not in database.thermo.libraries)
            details, cause, seen = [], exc, set()
            while cause is not None:
                if id(cause) in seen:
                    details.append('exception chain cycle detected')
                    break
                if len(seen) >= 64:
                    details.append('exception chain depth limit (64) reached')
                    break
                seen.add(id(cause))
                details.append('{0}: {1}'.format(type(cause).__name__, cause))
                cause = cause.__cause__ or cause.__context__
            raise ValueError('thermo library {0}: loading failed: {1}'.format(
                name, '; caused by '.join(details))) from exc
        yield database.thermo
    finally:
        rmg_data_module.database = saved


def _entry_reactor(entry, library_name, thermo):
    """Initialise an engine rate probe with the job's library-verified ion thermo."""
    import copy

    reaction = copy.copy(entry.item)
    reaction.kinetics = entry.data
    # Supply a complete bath/charge pool even for a purely thermal channel.
    pool = [Species(label='e-').from_adjacency_list('1 e u1 p0 c-1'),
            Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0'),
            Species(label='Arp').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')]
    if library_name == 'PlasmaRadiativeRecombination' and entry.index == 0:
        # The inactive Li+ law has no ion-convention thermo in this Ar job.
        # Its forward rate depends on the law, Tg/Te and molecularity, not
        # Li thermo. Use a verified Ar carrier of the same charge/order; the
        # raw lithium structure remains independently pinned above. This is
        # a rate probe, never an additional Report-10 chemistry channel.
        for side in ('reactants', 'products'):
            setattr(reaction, side, [(pool[2] if spc.get_net_charge() else pool[1])
                                    if spc.molecule[0].get_formula() == 'Li' else spc
                                    for spc in getattr(reaction, side)])
    core = list(dict.fromkeys(pool + reaction.reactants + reaction.products))
    if reaction.electrons or library_name in rmgpy.electron_placement.FAMILY_ELECTRON_PLACEMENT:
        # Use the same declared placement boundary as PlasmaReactor. Its owner
        # reads family on LibraryReaction, so carry that explicit provenance.
        reaction = LibraryReaction(reactants=reaction.reactants, products=reaction.products,
                                   kinetics=entry.data, library=library_name, entry=entry,
                                   reversible=reaction.reversible, electrons=reaction.electrons,
                                   specific_collider=reaction.specific_collider)
        # Reuse an explicit electron if the library already contains one.
        electron = next((spc for spc in core if spc.is_electron()), pool[0])
        core = [spc for spc in core if not spc.is_electron() or spc is electron]
        reaction = rmgpy.electron_placement.resolve_electron_placement(reaction, core)
    # Merge isomorphic bath fillers without changing the pinned participants.
    unique = []
    for spc in core:
        match = next((other for other in unique if spc.is_isomorphic(other)), None)
        if match is None:
            unique.append(spc)
        else:
            reaction.reactants = [match if item is spc else item for item in reaction.reactants]
            reaction.products = [match if item is spc else item for item in reaction.products]
    core = unique
    for spc in core:
        match = thermo.get_thermo_data_from_libraries(spc)
        if match is None and spc.get_net_charge() and not spc.is_electron():
            raise ValueError('no job thermo library entry for ion {0}'.format(spc.label))
        spc.thermo = match[0] if match is not None else _thermo(0.)
    fractions = {spc: (1.e-10 if spc.is_electron() or spc.get_net_charge() else 1.) for spc in core}
    fractions[next(spc for spc in core if spc.is_electron())] = sum(
        value * spc.get_net_charge() for spc, value in fractions.items() if not spc.is_electron())
    reactor = PlasmaReactor(T=(TG, 'K'), P=(PRESSURE, 'Pa'), Te=(11604.51812, 'K'),
                            initial_mole_fractions=fractions)
    reactor.initialize_model(core, [reaction], [], [])
    return reactor, reaction


def _engine_grid(reactor, reaction):
    """Read kf after the engine's real dispatch; never call a kinetics evaluator."""
    original_t, original_te = reactor.T, reactor.Te
    rates = []
    try:
        order = int(np.count_nonzero(reactor.reactant_indices[reactor.reaction_index[reaction]] >= 0))
        for tg in (250., TG, 600.):
            for te_ev in (.5, .9, 1., 2., 6.):
                reactor.T = Quantity(tg, 'K')
                reactor.Te = Quantity(te_ev * constants.e / constants.kB, 'K')
                reactor.generate_rate_coefficients([reaction], [])
                rates.append(float(reactor.kf[reactor.reaction_index[reaction]]) / constants.Na ** (order - 1))
    finally:
        reactor.T, reactor.Te = original_t, original_te
    return rates


def _rate_problems(expected, actual, order):
    problems = []
    units = 'm^6/s' if order == 3 else 'm^3/s' if order == 2 else 's^-1'
    grid = [(tg, te) for tg in (250., TG, 600.) for te in (.5, .9, 1., 2., 6.)]
    for (tg, te), wanted, found in zip(grid, expected, actual):
        relative = abs(found - wanted) / max(abs(wanted), 1e-300)
        if not np.isfinite(found) or relative > 1e-7:
            problems.append('rate differs at Tg={0:g} K, Te={1:g} eV ({2}): '
                            'expected={3:.17g}, actual={4:.17g}, relative difference={5:.6g}'
                            .format(tg, te, units, wanted, found, relative))
    return problems


def _report10_input_root():
    """Preserve configuration refusal and distinguish absent from corrupt paths."""
    import rmgpy

    try:
        path = Path(rmgpy.settings.require_database_directory())
    except SettingsError as exc:
        settings = rmgpy.settings
        if settings.filename is None or settings.sources.get('database.directory') == settings.DEFAULT_SOURCE:
            pytest.skip('Report-10 needs a configured plasma database: {0}'.format(exc))
        configured = settings.get('database.directory')
        if not configured:
            pytest.skip('Report-10 database input is absent: {0}'.format(exc))
        _require_input(configured, directory=True, absent_ok=True)
        pytest.fail('Report-10 configured input {0}: {1}'.format(configured, exc))
    return _require_input(path, directory=True, absent_ok=True)


def _report10_input_files(path):
    """Validate the same required deck content in preflight and scratch creation."""
    plasma_argon = _require_input(path / 'kinetics/libraries/PlasmaArgon', directory=True, absent_ok=True)
    required = []
    for name in ('PlasmaArgon', 'PlasmaRadiativeRecombination'):
        _require_input(path / 'kinetics/libraries' / name, directory=True)
        for file in ('reactions.py', 'dictionary.txt'):
            _require_input(path / 'kinetics/libraries' / name / file)
    required += ['thermo/libraries/{0}.py'.format(name) for name in REPORT10_THERMO_LIBRARIES]
    required += ['kinetics/families/Plasma_Electron_Attachment/groups.py',
                 'kinetics/families/Plasma_Electron_Attachment/rules.py']
    for file in required:
        _require_input(path / file)
    return plasma_argon


def _report10_database():
    """Require usable configured input, then compare every raw pinned entry."""
    import rmgpy

    path = _report10_input_root()
    plasma_argon = _report10_input_files(path)
    expected = json.loads((DATA / 'report10_database_contract.json').read_text())
    libraries, mismatches = _report10_libraries(path)
    with ExitStack() as stack:
        try:
            thermo = stack.enter_context(_report10_thermo(path))
        except Exception as exc:
            mismatches.append(str(exc))
            thermo = None
        for name, library in libraries.items():
            pinned = expected['libraries'][name]
            try:
                library_fields = _differing_fields(expected['library_metadata'][name],
                                                   _canonical_library(library), 'library')
                library_problem = ('library serialisation differs: ' + ', '.join(library_fields)
                                   if library_fields else None)
            except Exception as exc:
                library_problem = 'cannot serialise library: {0}: {1}'.format(type(exc).__name__, exc)
            for index in sorted(set(library.entries) | {int(key) for key in pinned}):
                identity = '{0}:{1}: '.format(name, index)
                if library_problem:
                    mismatches.append(identity + library_problem)
                if str(index) not in pinned:
                    mismatches.append(identity + 'unexpected reaction')
                    continue
                if index not in library.entries:
                    mismatches.append(identity + 'missing entry')
                    continue
                entry, wanted = library.entries[index], pinned[str(index)]
                actual_class = type(entry.data).__module__ + '.' + type(entry.data).__name__
                if actual_class != wanted['entry']['kinetics']['class']:
                    mismatches.append(identity + 'kinetics type changed (requires {0}; found {1})'.format(
                        wanted['entry']['kinetics']['class'], actual_class))
                try:
                    actual = _canonical_entry(entry, name)
                    fields = _differing_fields(wanted['entry'], actual)
                    nonfinite = [key for key in wanted['entry']['kinetics']
                                 if hasattr(getattr(entry.data, key, None), 'value_si')
                                 and not np.all(np.isfinite(getattr(entry.data, key).value_si))]
                    order = ('Tmin', 'Tmax', 'Pmin', 'Pmax')
                    nonfinite.sort(key=lambda key: (0, order.index(key)) if key in order else (1, key))
                    # New validity fields may be nonfinite even when the pin is None.
                    if nonfinite:
                        mismatches.append(identity + 'nonfinite parameters ' + ', '.join(nonfinite))
                    if thermo is None:
                        mismatches.append(identity + 'rate not evaluated because thermo failed to load')
                    elif not nonfinite:
                        try:
                            reactor, reaction = _entry_reactor(entry, name, thermo)
                            for problem in _rate_problems(wanted['rates'], _engine_grid(reactor, reaction),
                                                          int(np.count_nonzero(reactor.reactant_indices[0] >= 0))):
                                mismatches.append(identity + problem)
                        except Exception as exc:
                            mismatches.append(identity + 'rate evaluation failed: {0}: {1}'.format(type(exc).__name__, exc))
                    if fields:
                        details = []
                        if 'reaction.reactants' in fields:
                            details.append('reactant molecules differ')
                        if 'reaction.products' in fields:
                            details.append('product molecules differ')
                        for key in ('energies', 'sigma'):
                            if any(field.startswith('kinetics.' + key + '.') for field in fields):
                                details.append('cross section {0} differs'.format(key))
                        details.append('serialisation differs: ' + ', '.join(fields))
                        mismatches.append(identity + '; '.join(details))
                except Exception as exc:
                    mismatches.append(identity + 'cannot serialise: {0}: {1}'.format(type(exc).__name__, exc))
    if mismatches:
        pytest.fail('Report-10 configured database {0} has incompatible entries: {1}; '
                    'the characterisation must be regenerated'.format(path, '; '.join(mismatches)))
    return path, hashlib.sha256((plasma_argon / 'reactions.py').read_bytes()).hexdigest()


def _initialized_reaction_structure(reactor, reaction):
    """Serialize the actual solver participants with the pinned graph serializer."""
    species = {index: spc for spc, index in reactor.species_index.items()}
    row = reactor.reaction_index[reaction]
    structure = {'reversible': reaction.reversible}
    for side, indices in (('reactants', reactor.reactant_indices),
                          ('products', reactor.product_indices)):
        structure[side] = sorted((_canonical_species(species[int(index)])['adjacency']
                                  for index in indices[row] if index >= 0),
                                 key=lambda item: json.dumps(item, sort_keys=True))
    return structure


def _pinned_reaction_structure(entry):
    """Project pinned entry graphs into the declared reactor-facing electron view."""
    structure = {'reversible': entry['reaction']['reversible']}
    for side in ('reactants', 'products'):
        structure[side] = [spc['adjacency'] for spc in entry['reaction'][side]]
    placement = entry['electron_placement']
    if placement:
        electron = _canonical_species(Species().from_adjacency_list('1 e u1 p0 c-1'))['adjacency']
        if not any(electron in structure[side] for side in ('reactants', 'products')):
            for side, count in zip(('reactants', 'products'), placement):
                structure[side].extend([electron] * count)
    for side in ('reactants', 'products'):
        structure[side].sort(key=lambda item: json.dumps(item, sort_keys=True))
    return structure


def _report10_initialized_rates(reactor, reactions):
    """Compare initialized identities, actual solver structures, and engine rates."""
    contract = json.loads((DATA / 'report10_database_contract.json').read_text())
    pinned = contract['libraries']
    expected = {(name, index) for name, indices in contract['initialized_entries'].items()
                for index in indices}
    actual = {}
    for reaction in reactions:
        library = getattr(reaction, 'library', None)
        name = getattr(library, 'label', library) or '<no library>'
        entry = getattr(reaction, 'entry', None)
        index = getattr(entry, 'index', reaction.index)
        actual.setdefault((name, index), []).append(reaction)
    problems = ['{0}:{1}: missing initialized reaction'.format(*key) for key in sorted(expected - actual.keys())]
    problems += ['{0}:{1}: extra initialized reaction'.format(*key) for key in sorted(actual.keys() - expected)]
    problems += ['{0}:{1}: duplicate initialized reaction'.format(*key)
                 for key, members in sorted(actual.items()) if len(members) != 1]
    try:
        for key in sorted(expected & actual.keys()):
            name, index = key
            for reaction in actual[key]:
                identity = '{0}:{1}: '.format(name, index)
                try:
                    fields = _differing_fields(
                        _pinned_reaction_structure(pinned[name][str(index)]['entry']),
                        _initialized_reaction_structure(reactor, reaction))
                    if fields:
                        problems.append(identity + 'initialized structure differs: ' + ', '.join(fields))
                except Exception as exc:
                    problems.append(identity + 'cannot serialise initialized structure: {0}'.format(exc))
                try:
                    rates = _engine_grid(reactor, reaction)
                    problems.extend(identity + problem for problem in
                                    _rate_problems(pinned[name][str(index)]['rates'], rates,
                                                   int(np.count_nonzero(reactor.reactant_indices[reactor.reaction_index[reaction]] >= 0))))
                except Exception as exc:
                    problems.append(identity + 'rate evaluation failed: {0}'.format(exc))
    finally:
        reactor.generate_rate_coefficients(reactions, [])
    if problems:
        pytest.fail('; '.join(problems) + '; the characterisation must be regenerated')


def _report10_deck(tmp_path, monkeypatch):
    """Initialize the actual example without changing caller settings or DSL state."""
    import rmgpy

    root = Path(__file__).resolve().parents[3]
    settings = rmgpy.settings
    saved_values, saved_sources = dict(settings), dict(settings.sources)
    sources_object, filename = settings.sources, settings.filename
    # The all-module globals probe found these warning sets/disk caches change.
    # Preserve every quarantine warning/cache container, including nonempty ones.
    caches = [(name, value, value.copy()) for name, value in vars(quarantine).items()
              if name.endswith(('_WARNED', '_CACHE')) and isinstance(value, (dict, set))]
    # Inspect existing named RMG loggers without creating one in caller state.
    rmg_loggers = [logger for name, logger in logging.Logger.manager.loggerDict.copy().items()
                   if (name in ('rmgpy', 'RMG') or name.startswith('rmgpy.'))
                   and isinstance(logger, logging.Logger)]
    loggers = [(logger, logger.handlers, logger.handlers[:], logger.level, logger.propagate)
               for logger in [logging.getLogger()] + rmg_loggers]
    warning_filters = list(warnings.filters)
    with monkeypatch.context() as scoped:
        # Input loading rebinds these globals. Preserve identity, contents and
        # absence (mol_to_frag does not exist until the first input is read).
        scoped.setattr(rmg_input, 'rmg', rmg_input.rmg)
        scoped.setattr(rmg_input, 'species_dict', rmg_input.species_dict)
        scoped.setattr(rmg_input, 'mol_to_frag', {}, raising=False)
        scoped.setattr(rmg_data_module, 'database', None)
        scoped.chdir(tmp_path)
        # Numerical thermo loading also creates/updates __warningregistry__.
        # Give every loaded RMG module a private registry; undo preserves both
        # original registry identity/contents and an originally absent binding.
        for name, module in sys.modules.copy().items():
            if name == 'rmgpy' or name.startswith('rmgpy.'):
                registry = vars(module).get('__warningregistry__') or {}
                scoped.setattr(module, '__warningregistry__', registry.copy(), raising=False)
        try:
            _report10_database()
            job = RMG(input_file=str(root / 'examples/rmg/plasma_argon_energy_balance/input.py'),
                      output_directory=str(tmp_path))
            job.initialize()
            core = job.reaction_model.core.species
            rxns = job.reaction_model.core.reactions
            r = job.reaction_systems[0]
            solver_settings = job.simulator_settings_list[0]
            r.initialize_model(core, rxns, [], [], atol=solver_settings.atol, rtol=solver_settings.rtol)
            _report10_initialized_rates(r, rxns)
        finally:
            # Settings.__setitem__ changes provenance even during teardown.
            # Restore through dict primitives and keep the original sources object.
            dict.clear(settings)
            dict.update(settings, saved_values)
            settings.sources = sources_object
            sources_object.clear()
            sources_object.update(saved_sources)
            settings.filename = filename
            for name, original, contents in caches:
                setattr(quarantine, name, original)
                original.clear()
                original.update(contents)
            warnings.filters[:] = warning_filters
            for logger, handlers_object, handlers, level, propagate in loggers:
                for handler in logger.handlers[:]:
                    if handler not in handlers:
                        logger.removeHandler(handler)
                        handler.close()
                logger.handlers = handlers_object
                handlers_object[:] = handlers
                logger.setLevel(level)
                logger.propagate = propagate
    return r, core, rxns


def _channel_power(r):
    """Reconstruct each electron-power debit from real fluxes, in watts."""
    y, _ = _evaluate(r)
    return (np.asarray(r.core_reaction_rates) * r.compute_volume(y)
            * (r.energy_threshold + r.energy_electrons_consumed
               * 1.5 * constants.R * y[r.te_index]) * r.energy_participates)


@pytest.mark.database
def test_report10_omits_resonance_power_and_states_but_closes(tmp_path, monkeypatch):
    """C3 §10: real Report-10 omits resonance excitation altogether yet closes.

    LoKI-B assigns 39.8% to 1s4+1s2. Today's actual example has neither these
    states nor their direct power debit; it does not remove 39.8% into an
    energy-only sink. Entry 91's collapsed metastable mixing is a separate,
    much smaller debit. Explicit states and coupled radiation/return routes
    will supersede this characterisation.
    """
    database, fingerprint = _report10_database()
    # Provenance only: compatibility is decided by the per-entry contract.
    assert fingerprint
    r, core, rxns = _report10_deck(tmp_path, monkeypatch)
    _run(r, end_time=600.)
    y, rhs = _evaluate(r)
    powers = _channel_power(r)
    budget = _independent_power(r, core)
    # Independent flux accounting must also match the accepted solver's
    # stored-energy derivative, not only a zero-dydt residual identity.
    solver_stored = 1.5 * constants.R * (y[r.te_index] * r.dydt[r.electron_index]
                                        + y[r.electron_index] * r.dydt[r.te_index])
    assert abs(budget['closure'])/budget['power'] < 5e-6
    assert abs(budget['stored'])/budget['power'] < 1e-6
    assert abs(budget['power'] - sum(powers) - budget['elastic']
               - budget['wall_e'] - budget['wall_i'] - solver_stored)/budget['power'] < 5e-6
    row = _data()['ar']
    fixture = _resonance_fixture()
    reference = sum(ch['power_eVm3_s'] for ch in fixture['resonance_channels'])
    assert reference/row['power_eVm3_s']['Field']*100 == pytest.approx(39.8, rel=.002, abs=.05)
    resonance_species = [s for s in core if s.label in ('Ar(1s2)', 'Ar(1s4)')]
    excitation = [j for j, rxn in enumerate(rxns)
                  if any(s in resonance_species for s in rxn.products)]
    downstream = [rxn for rxn in rxns
                  if any(s in resonance_species for s in rxn.reactants)]
    q_res = sum(powers[j] for j in excitation)
    assert q_res == pytest.approx(0., rel=0., abs=1e-12)
    assert resonance_species == []
    assert downstream == []  # no radiation, superelastic, mixing or stepwise routes
    assert {s.label for s in core} == {'e-', 'Ar', 'Arp', 'Ars'}

    assert y[r.te_index] * constants.kB/constants.e == pytest.approx(.8958925, rel=.002, abs=.001)
    assert y[r.electron_index] * constants.Na/r.V == pytest.approx(1.1919978e16, rel=.003, abs=1e12)
    # Declared library-entry identities distinguish even the duplicate
    # Ars + e -> Ar + e superelastic and effective-mixing reactions.
    by_channel = {key: powers[j]/budget['power']
                  for j, key in enumerate(r.energy_reaction_keys) if key is not None}
    assert by_channel == pytest.approx({
        'PlasmaArgon:86': .0072759026,
        'PlasmaRadiativeRecombination:1': 5.2990921e-8,
        'PlasmaArgon:87': .193149541,
        'PlasmaArgon:88': .000122406111,
        'PlasmaArgon:89': -.000157999832,
        'PlasmaArgon:90': -5.05138653e-5,
        'PlasmaArgon:91': .000519357363,
    }, rel=.005, abs=1e-8)
    assert all(powers[j] == 0. for j, key in enumerate(r.energy_reaction_keys) if key is None)
    assert budget['elastic']/budget['power'] == pytest.approx(.79594136, rel=.003, abs=.001)
    assert (budget['wall_e']+budget['wall_i'])/budget['power'] == pytest.approx(.003199899, rel=.005, abs=1e-5)
    assert sum(powers) == pytest.approx(budget['inelastic'], rel=1e-12, abs=1e-8)

    # Ars is the sourced metastable reservoir, not a renamed resonance state.
    ar = next(s for s in core if s.label == 'Ar')
    ars = next(s for s in core if s.label == 'Ars')
    internal_gap = (ars.thermo.get_enthalpy(TG)-ar.thermo.get_enthalpy(TG))/constants.Na/constants.e
    assert internal_gap == pytest.approx(11.548354, rel=1e-5, abs=1e-5)
    thresholds = dict(zip(r.energy_reaction_keys, r.energy_threshold))
    assert thresholds['PlasmaArgon:87']/constants.Na/constants.e == pytest.approx(11.548354, rel=1e-6, abs=1e-8)
    assert thresholds['PlasmaArgon:91']/constants.Na/constants.e == pytest.approx(.075238, rel=1e-6, abs=1e-9)
    assert abs(rhs[r.te_index]/y[r.te_index]) < 1e-6
