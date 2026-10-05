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

"""Qualification interfaces for the collisional confined-anion wall.

Kemaneci et al., arXiv:1612.07268, equations (4)-(6), (13)-(14).
The radial formula below is a source-mapping check, NEVER the independent
full-profile or finite-cylinder reference. Both references and their frozen
thresholds must be explicitly supplied by an external qualification study.
"""
import math
import copy
import importlib

import numpy as np
from scipy.special import j1, jn_zeros

import rmgpy.constants as constants
from rmgpy.exceptions import ElectronegativeWallRegimeError, PlasmaStateError

CHI = float(jn_zeros(0, 1)[0])
BETA = float(2.0 * j1(CHI) / CHI)
# Engine geometry uses 2.405 instead of 2.4048255577. This is a fixed
# numerical mapping allowance, not the unavailable Report-10 error budget.
RADIAL_MAPPING_TOLERANCE = 2.e-4
REGIME_REQUIREMENTS = ('collisional', 'unmagnetised', 'confined_anions',
                       'homogeneous_profile', 'no_boundary_ionisation',
                       'isotropic_diffusion', 'compatible_boundaries')

O2_REFERENCE_QUALIFIED_UNITY = 'o2ReferenceQualifiedUnity'
O2_QUALIFYING_REFERENCE = (
    'plasma-pm3/reports/i314-en-wall-reference/rework3/report.md; '
    'I-314 rework 3, amended 2026-10-04')
O2_QUALIFICATION_ENVELOPE = tuple(
    dict(pressure_o2_torr=pressure, absorbed_power_w=power,
         reference_status=(
             'reference-unqualified (ion-heating validity), B3-pass numerically'
             if pressure == 0.05 and power == 1.0 else 'reference-qualified'))
    for pressure in (0.05, 0.1, 0.2, 0.25, 0.5, 1.0)
    for power in (0.25, 0.5, 1.0))


def closure_scientific_metadata(model):
    """Static owner-ruling metadata carried by diagnostics and manifests."""
    if model == 'confinedAnion':
        return dict(scientific_status='FALSIFIED SIMPLIFIED CLOSURE')
    if model == O2_REFERENCE_QUALIFIED_UNITY:
        return dict(
            scientific_status='REFERENCE-QUALIFIED UNITY CATION-WALL CLOSURE',
            claim_level='B3-QUALIFIED / MODEL-CONDITIONAL',
            claim_scope=(
                'closure qualification over the recorded I-314 O2 envelope only; '
                'the engine does not infer current-run envelope membership'),
            qualifying_reference=O2_QUALIFYING_REFERENCE,
            qualification_envelope=copy.deepcopy(O2_QUALIFICATION_ENVELOPE))
    return dict(scientific_status='FINITE-CYLINDER GEOMETRIC EXTENSION')


def wall_discrepancy(simple, reference, charges):
    """Non-cancelling charge-weighted symmetric wall-source discrepancy."""
    simple = np.asarray(simple, dtype=float)
    reference = np.asarray(reference, dtype=float)
    charges = np.asarray(charges, dtype=float)
    if (simple.shape != reference.shape or simple.shape != charges.shape
            or not np.all(np.isfinite(simple)) or not np.all(np.isfinite(reference))
            or not np.all(np.isfinite(charges)) or np.any(charges == 0.)
            or np.any(simple < 0.) or np.any(reference < 0.)):
        raise ElectronegativeWallRegimeError('C: invalid cation wall-source vector')
    # Scale first so the sums do not overflow; there is no epsilon floor.
    scale = max(float(np.max(simple, initial=0.)), float(np.max(reference, initial=0.)))
    if scale == 0.:
        return 0.
    a, b = simple / scale, reference / scale
    return float(2. * np.sum(np.abs(charges) * np.abs(a-b))
                 / np.sum(np.abs(charges) * (a+b)))


def radial_source_frequency(diffusivity, temperature, mass, radius, h):
    """Complete radial eq. (6), converted from centre to volume-average density.

    The qualified attachment branch of eq. (14) sets u_BE/u_B=1.
    No axial term appears in this source expression.
    """
    bohm = math.sqrt(constants.R * temperature / (constants.Na * mass))
    edge = h / math.hypot(1., radius * bohm / (CHI * diffusivity * j1(CHI)))
    return (2. / radius) * bohm * edge / BETA


# Numerical accuracy domain for explicitly selected EN closures. Existing
# transport and scientific qualification gates further restrict acceptance.
EN_WALL_DOMAIN = {
    'concentration_mol_m3': (1.e-30, 1.e6),
    'charged_concentration_mol_m3': (1.e-30, 1.e3),
    'kz/kr': (1.e-8, 1.e8),
    'Te_eV': (0.05, 50.),
    'volume_m3': (1.e-12, 1.e3),
    'cation_ep_frequency_s_1': (1.e-12, 1.e12),
    # SI mole units for reaction order; includes attachment, detachment and
    # destruction. Zero is allowed (an absent/suppressed channel). Kemaneci
    # et al., https://arxiv.org/abs/1612.07268, Tables III-VI: binary rates
    # ~1e-18..1e-12 m^3/s per particle (~1e6..1e12 in mole units), ternary
    # rates ~1e-47..1e-38 m^6/s (~1..1e10 in mole units), radiation ~1e7/s.
    # The PM cap has at least two decades above these scales, no lower cut.
    'rate_coefficient_si': (0., 1.e15),
    # PM parameter domain: about 1000 eV/event and broad rate-law powers.
    'reaction_energy_J_mol': (-1.e8, 1.e8),
    'rate_temperature_exponent': (-50., 50.),
    # Gas-phase Ar+ ~1.7..1.8 cm^2/(V s): Romero et al., Universe 8 (2022)
    # 134, https://doi.org/10.3390/universe8020134. O2-/O4- mu*N ~5.9e19
    # (V cm s)^-1: "Drift Velocities and Interactions of Negative Ions in
    # Oxygen. II.", Phys. Rev. A 4 (1971) 1043,
    # https://doi.org/10.1103/PhysRevA.4.1043 (~2.2 cm^2/(V s) at STP).
    # [1e-8,1] m^2/(V s) leaves >3 decades either side of these values.
    'reduced_mobility_m2_V_s': (1.e-8, 1.),
    'mobility_temperature_factor': (1.e-4, 1.e4),
    'transport_reference_temperature_K': (1., 1.e6),
    'transport_temperature_exponent': (-100., 100.),
    # Loschmidt reference 2.69e25/m^3; broad normalization allowance.
    'mobility_reference_density_m_3': (1.e20, 1.e30),
    # Kemaneci Table I: R=3e-3..8e-2 m, L=.20..0.49 m; the length and
    # eigenvalue limits give >2 decades of margin, including slender tubes.
    'transport_length_m': (1.e-6, 1.e3),
    'diffusion_eigenvalue_m_2': (1.e-10, 1.e14),
    # Neutral D*p=47 cm^2 Torr/s in the Ar fixture implies D*N~1.5e20
    # /(m s) at 300 K; >2 decades of margin in either direction.
    'neutral_diffusion_density_product_m_1_s_1': (1.e14, 1.e26),
    # Probabilities and dimensionless qualification cutoffs keep their exact
    # physical ranges; decades of margin cannot apply to a probability.
    'wall_recycling': (0., 1.),
    'wall_bath_threshold': (0., .01),
    'max_ionisation_degree': (0., 1.e6),
}


def check_wall_parameter(name, value, domain_key):
    """Refuse a declared parameter without clamping it."""
    low, high = EN_WALL_DOMAIN[domain_key]
    if not math.isfinite(value) or not low <= value <= high:
        raise PlasmaStateError('domain: {0} = {1!r}; expected [{2:g}, {3:g}] ({4})'.format(
            name, value, low, high, domain_key))


def check_wall_volume(volume):
    """Validate an explicitly supplied inventory volume with the shared bounds."""
    if not math.isfinite(volume) or volume <= 0.:
        raise ElectronegativeWallRegimeError('domain: volume must be finite and positive')
    low,high = EN_WALL_DOMAIN['volume_m3']
    if not low*(1.-4.*np.finfo(float).eps) <= volume <= high*(1.+4.*np.finfo(float).eps):
        raise ElectronegativeWallRegimeError('domain: volume = {0!r} m^3; expected [{1:g}, {2:g}]'.format(volume,low,high))


def check_wall_temperature(te_ev):
    """Check the temperature actually used by a public rate evaluator."""
    low, high = EN_WALL_DOMAIN['Te_eV']
    if not math.isfinite(te_ev) or not low*(1.-4.*np.finfo(float).eps) <= te_ev <= high*(1.+4.*np.finfo(float).eps):
        raise ElectronegativeWallRegimeError('domain: Te = {0!r} eV; expected [{1:g}, {2:g}]'.format(te_ev,low,high))


def check_wall_domain(y, volume, charges, labels, te_ev, components):
    """Refuse out-of-domain states before publishing an EN result.

    Zero concentrations are permitted; finite anions still require positive
    electrons at the regime gate. Bounds are inclusive and in SI, independent
    of the arbitrary mole-inventory normalization used by the constant-P EOS.
    """
    check_wall_volume(volume)
    for value, charge, label in zip(y, charges, labels):
        c = float(value)/volume
        bounds = EN_WALL_DOMAIN['charged_concentration_mol_m3' if charge else 'concentration_mol_m3']
        if not math.isfinite(c) or (c != 0. and not bounds[0]*(1.-4.*np.finfo(float).eps) <= c <= bounds[1]*(1.+4.*np.finfo(float).eps)):
            raise ElectronegativeWallRegimeError(
                'domain: concentration of {0} = {1!r} mol/m^3; expected 0 or [{2:g}, {3:g}]'.format(label,c,*bounds))
    check_wall_temperature(te_ev)
    if components is None:
        raise ElectronegativeWallRegimeError('domain: kz/kr requires both diffusion eigenvalues')
    kr,kz = components
    ratio = kz/kr if kr > 0. else float('inf')
    low,high = EN_WALL_DOMAIN['kz/kr']
    if not math.isfinite(ratio) or not low*(1.-4.*np.finfo(float).eps) <= ratio <= high*(1.+4.*np.finfo(float).eps):
        raise ElectronegativeWallRegimeError('domain: kz/kr = {0!r}; expected [{1:g}, {2:g}]'.format(ratio,low,high))


def check_wall_frequencies(frequencies, charges, labels):
    """Bound the absolute transport timescale as well as the geometry ratio."""
    low, high = EN_WALL_DOMAIN['cation_ep_frequency_s_1']
    for frequency, charge, label in zip(frequencies, charges, labels):
        if charge > 0 and (not math.isfinite(frequency) or not
                low*(1.-4.*np.finfo(float).eps) <= frequency <= high*(1.+4.*np.finfo(float).eps)):
            raise ElectronegativeWallRegimeError(
                'domain: EP wall frequency of {0} = {1!r} s^-1; expected [{2:g}, {3:g}]'.format(label, frequency, low, high))


def electron_rate_derivative(kin, temperature, coefficient, temperature_power=0.):
    """Analytic Te derivative of the supported electron rate laws.

    With a temperature power p, return d(k*T**p)/dT / T**p.
    Shift the analytic exponent before evaluation so k*T cancellation is
    exact, including tiny activation and EOS pressure contributions.

    The cross-section law differentiates the same trapezoidal quadrature as
    its evaluator, rather than finite-differencing reaction propensities.
    """
    name = type(kin).__name__
    T = temperature
    scalar = lambda v: float(v.value_si if hasattr(v,'value_si') else v)
    if name == 'TwoTemperaturePlasma':
        return coefficient*((scalar(kin.n)+temperature_power)/T + scalar(kin.Ea_e)/constants.R/T/T)
    if name == 'ElectronCollisionPlasma':
        e = np.asarray(kin.energies.value_si)/constants.Na
        sigma = np.asarray(kin.sigma.value_si)
        terms = sigma*e*np.exp(-e/(constants.kB*T))*(e/(constants.kB*T)+(temperature_power-1.5))/T
        prefactor = math.sqrt(8./(constants.pi*constants.m_e))*(constants.kB*T)**-1.5*constants.Na
        return prefactor*np.sum(.5*(terms[:-1]+terms[1:])*np.diff(e))
    if name == 'BadnellRRArrhenius':
        t0,t1 = math.sqrt(T/scalar(kin.T0)),math.sqrt(T/scalar(kin.T1))
        C = scalar(kin.C) if kin.C is not None else 0.
        T2 = scalar(kin.T2) if kin.T2 is not None else 0.
        extra = C*math.exp(-T2/T) if T2 > 0. else 0.
        b = scalar(kin.B)+extra
        db = extra*T2/T/T
        slope = (temperature_power-.5)/T-(1.-b)*t0/(2.*T*(1.+t0))-(1.+b)*t1/(2.*T*(1.+t1))
        slope += db*(math.log1p(t0)-math.log1p(t1))
        return coefficient*slope
    if name == 'VoronovEIArrhenius':
        U = scalar(kin.dE)/(constants.kB*T/constants.e)
        if U <= 1.e-16:
            return coefficient*temperature_power/T
        p,x,k = scalar(kin.P),scalar(kin.X),scalar(kin.K)
        return coefficient*((temperature_power-k)+U+U/(x+U)-.5*p*math.sqrt(U)/(1.+p*math.sqrt(U)))/T
    raise ElectronegativeWallRegimeError('Jacobian: no analytic Te derivative for '+name)


def electron_rate_temperature_factor(kin, temperature, coefficient, temperature_power=0.):
    """Dimensionless shifted slope, with the rate prefactor kept separate.

    The solver multiplies the rate, slope and mass-action factors together
    with exponent protection. Forming d(k*T)/dT first can overflow even
    when the complete state derivative is representable.
    """
    if type(kin).__name__ == 'TwoTemperaturePlasma':
        return (float(kin.n.value_si)+temperature_power
            +float(kin.Ea_e.value_si)/constants.R/temperature)
    if type(kin).__name__ == 'ElectronCollisionPlasma':
        if coefficient == 0.:
            return temperature_power
        return (electron_rate_derivative(kin,temperature,coefficient,temperature_power)
            /coefficient*temperature)
    return electron_rate_derivative(kin,temperature,1.,temperature_power)*temperature


def closure_factor_gradient(y, electron, anions, model, geometry, components):
    """Return h, active frequency multiplier, and its analytic state gradient.

    Clipping applies only to Newton trials. At the zero-anion boundary the
    legacy operator is taken exactly, including its Jacobian.
    """
    y = np.asarray(y)
    gradient = np.zeros(len(y))
    minus = sum(max(float(y[j]), 0.) for j in anions)
    if model is None or minus == 0. or model == O2_REFERENCE_QUALIFIED_UNITY:
        return 1., 1., gradient
    ne = max(float(y[electron]), 0.)
    total = ne + minus
    h = ne / total
    weight = 1. if geometry == 'fullFrequency' else components[0] / sum(components)
    gradient[electron] = weight * (minus / total) / total if y[electron] >= 0. else 0.
    for j in anions:
        if y[j] >= 0.:
            gradient[j] = -weight * h / total
    factor = h if geometry == 'fullFrequency' else (h*components[0]+components[1])/sum(components)
    return h, factor, gradient


def check_reference(gate, reference, threshold, context, transition=False):
    """Evaluate an externally supplied reference monitor without using it in the ODE.

    The reference returns positive per-cation mol/s, domain status, provenance,
    and a declaration of independence; full-profile C additionally requires
    a declaration that its qualified model spans the profile transition.
    """
    if isinstance(reference,str):
        try:
            module, name = reference.split(':',1)
            reference = importlib.import_module(module)
            for part in name.split('.'):
                reference = getattr(reference,part)
        except Exception as exc:
            raise ElectronegativeWallRegimeError(gate + ': reference cannot be imported: ' + str(exc)) from exc
    if (not callable(reference) or isinstance(threshold, bool)
            or not isinstance(threshold, (int, float))
            or not np.isfinite(threshold) or not 0. <= threshold <= 2.):
        raise ElectronegativeWallRegimeError(gate + ': missing reference or frozen threshold')
    try:
        # The independent monitor may use mutable working arrays; it cannot
        # redefine the production vector against which its result is tested.
        simple = np.array(context['simple_cation_flux'], copy=True)
        charges = np.array(context['cation_charges'], copy=True)
        result = reference(copy.deepcopy(context))
        if (not isinstance(result, dict) or result.get('in_domain') is not True
                or result.get('independent') is not True or not result.get('provenance')
                or transition and result.get('spans_transition') is not True):
            raise ValueError('reference is outside its qualified domain or lacks independent provenance')
        error = wall_discrepancy(simple, result['cation_flux'], charges)
    except Exception as exc:
        raise ElectronegativeWallRegimeError(gate + ': reference refused: ' + str(exc)) from exc
    if error > threshold:
        raise ElectronegativeWallRegimeError('{}: error={} exceeds frozen threshold={}'.format(gate,error,threshold))
    diagnostics = {key: copy.deepcopy(value) for key, value in result.items()
                   if key not in ('cation_flux', 'in_domain', 'independent', 'provenance', 'spans_transition')}
    return dict(error=error, threshold=threshold, provenance=result['provenance'],
                reference_diagnostics=diagnostics)


def serializable_qualification(qualification):
    """Persist importable references by module:qualified-name, never repr or eval."""
    result = copy.deepcopy(qualification)
    for key in ('full_profile_reference','geometry_reference'):
        value = result.get(key)
        if callable(value):
            module, name = value.__module__, value.__qualname__
            if '<locals>' in name or '<lambda>' in name or module == '__main__':
                raise ValueError('qualification references must be importable for input-file restart')
            result[key] = module + ':' + name
    return result


def manifest_values(value):
    """Convert arrays/scalars and mathematically infinite confinement to JSON values."""
    if isinstance(value,dict):
        return {key:manifest_values(item) for key,item in value.items()}
    if isinstance(value,(list,tuple,np.ndarray)):
        return [manifest_values(item) for item in value]
    if isinstance(value,(float,np.floating)):
        if np.isnan(value):
            return 'unresolved'
        if np.isinf(value):
            return 'infinite' if value > 0. else 'negative-infinite'
        return float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value,np.integer):
        return int(value)
    return value


def qualify_envelope(samples):
    """Qualify explicitly supplied (reactor, state, time) envelope samples.

    The caller must cover the declared pressure/power/composition/geometry and
    chemistry envelope; this helper does not infer or reduce it. Every sample
    is checked through the production accepted-state monitor. Return extrema
    with their sample locations and margins to the frozen gates.
    """
    extrema = {}
    count = 0
    for count,(reactor,state,time) in enumerate(samples,1):
        # Validate metadata before the monitor can replace the last valid state.
        if not np.isfinite(time) or time < 0.:
            raise ElectronegativeWallRegimeError('envelope: time must be finite and nonnegative')
        reactor.monitor_electronegative_wall(np.asarray(state),time)
        gates = reactor.electronegative_wall_diagnostics.get('gates')
        if not isinstance(gates,dict):
            continue
        conf,label = min((row['conf'],label) for label,row in gates['A'].items())
        values = dict(min_conf=conf,min_attachment_metric=min(gates['B']['metrics']),
                      max_full_profile_error=gates['C_full_profile']['error'],
                      max_geometry_error=gates['C_geometry']['error'])
        for key,value in values.items():
            prior = extrema.get(key)
            if prior is None or (value < prior['value'] if key.startswith('min') else value > prior['value']):
                boundary = (1. if key.startswith('min') else gates[
                    'C_full_profile' if key == 'max_full_profile_error' else 'C_geometry']['threshold'])
                extrema[key] = dict(value=value,margin=(value-boundary if key.startswith('min') else boundary-value),
                    sample=count, time=time, state=np.asarray(state).tolist(),
                    pressure=reactor.P.value_si,temperature=reactor.T.value_si,
                    geometry=reactor.wall_diffusion_components,
                    anion=label if key == 'min_conf' else None)
    if count == 0:
        raise ValueError('qualification envelope must not be empty')
    return manifest_values(dict(samples=count,extrema=extrema))


def classify_charged_reaction(reactants, products, charges, electron, charged_collider=False):
    """I-314 adapter classes; unsupported physics remains visible to the reference.

    This classification never changes the engine chemistry. Directed repeated
    indices and full stoichiometry accompany it, so a reference can narrow its
    own domain further rather than trusting a label as scientific qualification.
    """
    indices = tuple(reactants) + tuple(products)
    if charged_collider:
        return 'unsupported', 'charged third body'
    if any(j != electron and abs(charges[j]) > 1 for j in indices):
        return 'unsupported', 'multiply charged ion'
    ne = reactants.count(electron)
    positive = sum(charges[j] > 0 for j in reactants if j != electron)
    negative = sum(charges[j] < 0 for j in reactants if j != electron)
    net_positive = sum(charges[j] > 0 for j in products if j != electron) - positive
    net_negative = sum(charges[j] < 0 for j in products if j != electron) - negative
    if positive == negative == 0:
        if net_positive == net_negative == 0 and products.count(electron) == ne:
            return 'neutral', None
        if ne == 1 and net_positive > 0 and net_negative == 0:
            return 'electron_ionisation', None
        if ne == 1 and net_negative > 0 and net_positive == 0:
            return 'electron_attachment', None
        return 'unsupported', 'charged production outside electron-impact classes'
    if ne == 0 and positive == 0 and negative == 1:
        return 'anion_neutral', None
    if ne == 1 and positive == 0 and negative == 1:
        return 'anion_electron', None
    if ne == 1 and positive == 1 and negative == 0:
        return 'cation_electron', None
    if ne == 0 and positive == 1 and negative == 1:
        return 'cation_anion', None
    if ne == 0 and positive == 1 and negative == 0:
        return 'cation_neutral', None
    return 'unsupported', 'charged multibody reaction outside reference classes'
