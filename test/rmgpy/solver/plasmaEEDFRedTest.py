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
from rmgpy.exceptions import DatabaseError, InvalidAdjacencyListError, PlasmaStateError, SettingsError
from rmgpy.kinetics import Arrhenius, ElectronCollisionPlasma, TwoTemperaturePlasma
from scipy.integrate import quad
from rmgpy.solver.plasma import PlasmaReactor
import rmgpy.rmg.input as rmg_input
from rmgpy.rmg.main import RMG
from rmgpy.species import Species
from rmgpy.thermo import ThermoData

DATA = Path(__file__).parent / 'data'
AR_ION = 'e + Ar(1S0) -> e + e + Ar(+,gnd), Ionization'
AR_M1 = 'e + Ar(1S0) -> e + Ar(3P2), Excitation'
AR_M2 = 'e + Ar(1S0) -> e + Ar(3P0), Excitation'
AR_ELASTIC = 'e + Ar(1S0) -> e + Ar(1S0), Elastic'
N2_ION = 'e + N2(X) -> e + e + N2(+,X), Ionization'
TG = 298.15
PRESSURE = 5.0 * 101325.0 / 760.0


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
    assertions = ['Ar+']
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
        assertions.append('N2+')
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


def _report10_database():
    """Require the configured database and this example's real plasma content."""
    import rmgpy

    try:
        path = Path(rmgpy.settings.require_database_directory())
    except SettingsError as exc:
        pytest.skip('Report-10 needs a configured plasma database: {0}'.format(exc))
    plasma_argon = path / 'kinetics/libraries/PlasmaArgon'
    if not plasma_argon.is_dir():
        pytest.skip('Report-10 configured database {0} lacks the PlasmaArgon library directory'.format(path))
    required_files = [
        'thermo/libraries/{0}.py'.format(name) for name in
        ('primaryThermoLibrary', 'PlasmaThermo', 'electrocatThermo', 'PlasmaExcitedNeutralThermo')]
    required_files += [
        'kinetics/libraries/{0}/{1}'.format(name, filename)
        for name in ('PlasmaArgon', 'PlasmaRadiativeRecombination')
        for filename in ('reactions.py', 'dictionary.txt')]
    required_files += ['kinetics/families/Plasma_Electron_Attachment/groups.py',
                       'kinetics/families/Plasma_Electron_Attachment/rules.py']
    missing = [name for name in required_files if not (path / name).is_file()]
    if missing:
        pytest.fail('Report-10 configured database {0} lacks required content: {1}; '
                    'the characterisation must be regenerated'.format(path, ', '.join(missing)))
    from rmgpy.reaction import same_species_lists

    molecules = {
        'Ar': Species().from_adjacency_list('1 Ar u0 p4 c0'),
        'Arp': Species().from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1'),
        'Ars': Species().from_adjacency_list('multiplicity 3\n1 Ar u2 p3 c0'),
        'e-': Species().from_adjacency_list('1 e u1 p0 c-1'),
    }

    def two_temp(a, n, gap):
        return TwoTemperaturePlasma(A=(a, 'm^3/(molecule*s)'), n=n,
                                    Ea_g=(gap, 'eV/molecule'), Ea_e=(gap, 'eV/molecule'),
                                    T0=(11604.51812, 'K'))

    # These laws are the current database contract.  We deliberately compare
    # evaluated rates rather than storage parameters: e.g., T0 is irrelevant
    # when n is zero.  The grid spans this example's 298 K bath and its
    # approximately 0.9 eV converged electron state with conservative bounds.
    expected = {
        86: (('Ar', 'e-'), ('Arp', 'e-', 'e-'), _collision(_data(), 'golyatina')),
        87: (('Ar', 'e-'), ('Ars', 'e-'), _ashida()),
        88: (('Ars', 'e-'), ('Arp', 'e-', 'e-'), two_temp(1.187793909e-13, .3278287422, 4.401535934)),
        89: (('Ars', 'e-'), ('Ar', 'e-'), two_temp(4.3e-16, .74, 0.)),
        90: (('Ars', 'Ars'), ('Ar', 'Arp', 'e-'),
             Arrhenius(A=(6.2e-16, 'm^3/(molecule*s)'), n=0, Ea=(0., 'J/mol'), T0=(1., 'K'))),
        91: (('Ars', 'e-'), ('Ar', 'e-'), two_temp(2.e-13, 0., 0.)),
        92: (('Ar', 'Ars'), ('Ar', 'Ar'),
             Arrhenius(A=(3.e-15, 'cm^3/(molecule*s)'), n=0, Ea=(0., 'J/mol'), T0=(1., 'K'))),
        93: (('Ar', 'Ar', 'Ars'), ('Ar', 'Ar', 'Ar'),
             Arrhenius(A=(1.1e-31, 'cm^6/(molecule^2*s)'), n=0, Ea=(0., 'J/mol'), T0=(1., 'K'))),
    }
    recombination = TwoTemperaturePlasma(A=(3.77e-13, 'cm^3/(molecule*s)'), n=-.651,
                                         Ea_g=(0., 'J/mol'), Ea_e=(0., 'J/mol'),
                                         T0=(1.e4, 'K'), electrons=-1)
    context = KineticsDatabase().local_context
    mismatches = []
    for name, required in (('PlasmaArgon', expected),
                           ('PlasmaRadiativeRecombination', {1: (('Arp',), ('Ar',), recombination)})):
        library = KineticsLibrary(label=name)
        try:
            library.load(str(path / 'kinetics/libraries' / name / 'reactions.py'), local_context=context.copy())
        except (DatabaseError, InvalidAdjacencyListError, OSError, KeyError, ValueError, TypeError, SyntaxError) as exc:
            source = (path / 'kinetics/libraries' / name / 'reactions.py').read_text().lower()
            if name == 'PlasmaArgon' and 'nan' in source:
                mismatches.append('PlasmaArgon:86: nonfinite parameters sigma')
            else:
                mismatches.append('{0}: cannot load ({1})'.format(name, exc))
            continue
        if name == 'PlasmaArgon':
            for unexpected in sorted(set(library.entries) - set(required)):
                mismatches.append('PlasmaArgon:{0}: unexpected reaction'.format(unexpected))
        for index, (reactants, products, kinetics) in required.items():
            entry = library.entries.get(index)
            if entry is None:
                mismatches.append('{0}:{1}: missing entry'.format(name, index))
                continue

            def incompatible(content):
                mismatches.append('{0}:{1}: {2}'.format(name, index, content))

            for side, tokens in (('reactant', reactants), ('product', products)):
                actual = entry.item.reactants if side == 'reactant' else entry.item.products
                if not same_species_lists(actual, [molecules[token] for token in tokens], strict=True):
                    incompatible('{0} molecules (requires {1})'.format(side, ' + '.join(tokens)))
            if entry.item.reversible:
                incompatible('reaction direction (requires irreversible)')
            if not isinstance(entry.data, type(kinetics)):
                incompatible('kinetics type changed (requires {0}; found {1})'.format(
                    type(kinetics).__name__, type(entry.data).__name__))
                continue
            parameters = ('energies', 'sigma') if isinstance(kinetics, ElectronCollisionPlasma) else (
                ('A', 'n', 'T0', 'Ea_g', 'Ea_e', 'electrons') if isinstance(kinetics, TwoTemperaturePlasma)
                else ('A', 'n', 'T0', 'Ea'))
            nonfinite = [p for p in parameters if not np.all(np.isfinite(np.asarray(getattr(entry.data, p).value_si)))]
            if nonfinite:
                incompatible('nonfinite parameters ' + ', '.join(nonfinite))
                continue
            if isinstance(kinetics, ElectronCollisionPlasma):
                for parameter in ('energies', 'sigma'):
                    actual, wanted = (np.asarray(getattr(entry.data, parameter).value_si),
                                      np.asarray(getattr(kinetics, parameter).value_si))
                    if actual.shape != wanted.shape or not np.allclose(actual, wanted, rtol=1e-12, atol=0.):
                        incompatible('cross section {0} differs'.format(parameter))
                continue
            differences = []
            for tg in (250., TG, 600.):
                for te_ev in (.5, .9, 1., 2., 6.):
                    te = te_ev * constants.e / constants.kB
                    divisor = constants.Na ** (len(reactants) - 1)
                    units = 'm^6/s' if len(reactants) == 3 else 'm^3/s'
                    wanted = (kinetics.get_rate_coefficient_two_temp(tg, te) / divisor
                              if isinstance(kinetics, TwoTemperaturePlasma) else kinetics.get_rate_coefficient(tg) / divisor)
                    actual = (entry.data.get_rate_coefficient_two_temp(tg, te) / divisor
                              if isinstance(entry.data, TwoTemperaturePlasma) else entry.data.get_rate_coefficient(tg) / divisor)
                    if not np.isfinite(actual):
                        differences.append('nonfinite rate at Tg={0:g} K, Te={1:g} eV'.format(tg, te_ev))
                    elif abs(actual - wanted) / max(abs(wanted), 1e-300) > 1e-7:
                        differences.append('rate differs at Tg={0:g} K, Te={1:g} eV ({2})'.format(tg, te_ev, units))
            if differences:
                incompatible('; '.join(differences))
    if mismatches:
        pytest.fail('Report-10 configured database {0} has incompatible entries: {1}; '
                    'the characterisation must be regenerated'.format(path, '; '.join(mismatches)))
    fingerprint = hashlib.sha256((plasma_argon / 'reactions.py').read_bytes()).hexdigest()
    return path, fingerprint


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
