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

"""Optically thick resonance radiation with complete frequency redistribution.

Holstein, Phys. Rev. 83, 1159 (1951), Tables I/III and Eqs.6.2--6.6;
Walsh, Phys. Rev. 116, 511 (1959), Eqs.2.7--2.9. The optical length is
full gap for a slab and radius for an infinite cylinder, never a diffusion length.
The slab uses the same lifetime interpolation with Holstein's slab constants;
the mixed emission/absorption term has the leading-log geometry-independent
coefficient one. This is an asymptotic CRD closure, not finite-cylinder transport.
"""

import math


def walsh_escape_factor(optical_depth, damping, geometry):
    """Return g for Doppler line-centre depth tau and Voigt damping a.

    ``a`` is Lorentz FWHM / twice the Gaussian 1/e half-width. No optically
    thin continuation is supplied. Geometry coefficients describe the lowest
    imprisonment mode, not the 1.9/1.3 transmission fit used by some CR models.
    """
    return _walsh_terms(optical_depth, damping, geometry)[0]


def _walsh_terms(tau, a, geometry):
    if geometry not in ('slab', 'infiniteCylinder'):
        raise ValueError('radiation geometry must be slab or infiniteCylinder')
    if not math.isfinite(tau) or tau < 1.e4:
        raise ValueError('radiation requires finite Doppler optical depth >= 10000; thin/crossover out of scope')
    if not math.isfinite(a) or not 0. <= a <= .5:
        raise ValueError('radiation requires finite Voigt damping in [0, 0.5]')
    cd, cc = {'slab': (1.875, 1.150), 'infiniteCylinder': (1.60, 1.115)}[geometry]
    log_tau = math.log(tau)
    gd = cd / (tau * math.sqrt(math.pi * log_tau))
    if a == 0.:
        slope_a = (4. / (math.pi**1.5 * math.sqrt(log_tau))
                   - gd * 4. * tau / (math.pi**1.5 * cc*cc * log_tau))
        return gd, gd * (-1. - .5/log_tau) / tau, slope_a
    gc = cc * math.sqrt(a / (math.sqrt(math.pi) * tau))
    gcd = 2. * a / (math.pi * math.sqrt(log_tau))
    q = gcd / gc
    exp_q = math.exp(-q*q)
    erf_q = math.erf(q)
    value = gd*exp_q + gc*erf_q
    slope_q = exp_q * (-2.*q*gd + 2.*gc/math.sqrt(math.pi))
    slope_tau = (gd*exp_q*(-1.-.5/log_tau) - .5*gc*erf_q
                 + slope_q*q*(.5-.5/log_tau))/tau
    slope_a = (.5*gc*erf_q + .5*slope_q*q)/a
    return value, slope_tau, slope_a


def walsh_escape_derivatives(optical_depth, damping, geometry):
    """Return partial derivatives (dg/dtau, dg/da) in the admitted domain."""
    return _walsh_terms(optical_depth, damping, geometry)[1:]


def load_argon_lines(path):
    """Read the sourced vacuum lines, independently of cascade branching."""
    import json
    with open(path) as stream:
        data = json.load(stream)
    if data.get('schema_version') != 1 or set(data.get('lines', {})) != {'3P1', '1P1'}:
        raise ValueError('radiation needs the two resolved argon resonance lines')
    for state, line in data['lines'].items():
        for key in ('wavelength_m', 'A_s_1', 'g_ratio'):
            if not math.isfinite(line[key]) or line[key] <= 0.:
                raise ValueError('radiation invalid line value: '+state+' '+key)
        if line['g_ratio'] != 3. or line['lower_nist_id'] != '018001.000001':
            raise ValueError('radiation line must connect J=1 argon to ground J=0')
    return data


def argon_escape(line, ground_density, temperature, length, geometry):
    """Return (g, dg/dN_Ar, tau_D, a) from SI state values.

    NIST Spectral Line Shapes Eq.41 supplies resonance self-broadening.
    Natural width A/(2*pi) is retained. No foreign linewidth is silently
    folded into the ground-Ar opacity.
    """
    if not math.isfinite(temperature) or not 250. <= temperature <= 500.:
        raise ValueError('radiation gas temperature must be in [250, 500] K')
    if not math.isfinite(ground_density) or ground_density <= 0.:
        raise ValueError('radiation requires positive finite ground-Ar density')
    if not math.isfinite(length) or length <= 0.:
        raise ValueError('radiation geometry length must be finite and positive')
    lam, rate, ratio = line['wavelength_m'], line['A_s_1'], line['g_ratio']
    width_d = math.sqrt(2.*1.380649e-23*temperature/(39.948*1.66053906660e-27))/lam
    tau = lam**2*ground_density*ratio*rate*length/(8.*math.pi**1.5*width_d)
    angstrom = lam*1e10
    oscillator = rate*ratio*angstrom**2/6.6702e15
    width_res = (299792458./lam**2 * 1e-10 * 8.6e-30
                 * ratio**-.5 * angstrom**3 * oscillator * ground_density/1e6)
    a_res = width_res/(2.*width_d)
    a = (rate/(2.*math.pi)+width_res)/(2.*width_d)
    g, d_tau, d_a = _walsh_terms(tau,a,geometry)
    return g, (d_tau*tau+d_a*a_res)/ground_density, tau, a


def _argon_state(species):
    """Return the resolved neutral atomic Ar state (empty means ground)."""
    if not species.molecule:
        return None
    states = set()
    for mol in species.molecule:
        if len(mol.atoms) != 1 or mol.atoms[0].element.symbol != 'Ar' or mol.get_net_charge() != 0:
            return None
        if not mol.electronic_state and mol.multiplicity != 1:
            return None
        states.add(mol.electronic_state)
    if len(states) != 1:
        raise ValueError('radiation requires unambiguous argon state identity')
    return states.pop()


class ArgonRadiation:
    """Reactor-owned operator; all rates are functions of the supplied state."""

    def __init__(self, declaration, geometry, line_data):
        import copy
        if not isinstance(declaration, dict) or set(declaration) != {'model', 'geometryModel'}:
            raise ValueError('resonanceRadiation requires model and geometryModel')
        if declaration['model'] not in ('holsteinWalsh', 'complete'):
            raise ValueError('resonanceRadiation model must be holsteinWalsh or complete')
        self.declaration = copy.deepcopy(declaration)
        if not isinstance(geometry, dict):
            raise ValueError('radiation requires chamber geometry')
        mode = declaration['geometryModel']
        if mode == 'infiniteCylinder' and set(geometry) == {'shape', 'radius', 'length'} and geometry['shape'] == 'cylinder':
            radius, axial = geometry['radius'], geometry['length']
            if not all(math.isfinite(x) and x > 0 for x in (radius,axial)) or axial < 6.*radius*(1.-1e-14):
                raise ValueError('radiation infinite-cylinder geometry requires length >= 6*radius')
            self.length = radius
        elif mode == 'slab' and set(geometry) == {'shape','gap'} and geometry['shape'] == 'slab':
            self.length = geometry['gap']
            if not math.isfinite(self.length) or self.length <= 0.:
                raise ValueError('radiation slab geometry requires a positive finite full gap')
        else:
            raise ValueError('radiation geometry must match slab or infiniteCylinder; diffusionLength is insufficient')
        self.geometry = copy.deepcopy(geometry)
        self.line_data = copy.deepcopy(line_data)
        self.transitions = []
        self.ground = None
        self.neutral = []
        self.argon = []
        self.nitrogen = []
        self.charged = []
        self.size = 0

    def bind(self, core, edge):
        """Resolve by molecular identity on every core growth/restart."""
        self.size = len(core)
        self.transitions = []
        grounds = [i for i,s in enumerate(core) if _argon_state(s) == '']
        present = any(_argon_state(s) in ('3P1','1P1') for s in core+edge)
        if present and len(grounds) != 1:
            raise ValueError('radiation requires exactly one ground Ar species in core')
        self.ground = grounds[0] if grounds else None
        for state in ('3P1','1P1'):
            sites = [i for i,s in enumerate(core+edge) if _argon_state(s) == state]
            if len(sites) > 1:
                raise ValueError('radiation duplicate resolved resonance state '+state)
            if sites and sites[0] < len(core):
                self.transitions.append((sites[0],state))
        self.neutral, self.argon, self.nitrogen, self.charged = [], [], [], []
        for i, spc in enumerate(core):
            if spc.is_electron():
                continue
            if spc.get_net_charge():
                self.charged.append(i)
                continue
            self.neutral.append(i)
            if all(len(mol.atoms) == 1 and mol.atoms[0].element.symbol == 'Ar'
                   for mol in spc.molecule):
                self.argon.append(i)
            elif all(len(mol.atoms)==2 and all(atom.element.symbol=='N' for atom in mol.atoms) for mol in spc.molecule):
                self.nitrogen.append(i)

    def evaluate(self, y, volume, temperature, volume_gradient=None):
        """Return mol/s species rates, optional Jacobian, and unlatched diagnostics."""
        import numpy as np
        na = 6.02214076e23
        if not math.isfinite(volume) or volume <= 0.:
            raise ValueError('radiation requires positive finite volume')
        rate = np.zeros(self.size)
        jac = None if volume_gradient is None else np.zeros((self.size,len(volume_gradient)))
        records = []
        if not self.transitions:
            return rate, jac, dict(lines=[],photon_power_W=0.)
        density = na*y[self.ground]/volume
        if self.declaration['model'] == 'holsteinWalsh':
            neutral = sum(max(y[i],0.) for i in self.neutral)
            ar = sum(max(y[i],0.) for i in self.argon)
            nitrogen = sum(max(y[i],0.) for i in self.nitrogen)
            if neutral <= 0. or y[self.ground] < .75*neutral or neutral-ar-nitrogen > 1e-6*neutral:
                raise ValueError('radiation domain requires ground-Ar fraction >= 0.75, with N2 as the only foreign bath')
            if sum(max(y[i],0.) for i in self.charged) > 1e-3*neutral:
                raise ValueError('radiation domain requires ionization fraction <= 1e-3')
        photon = 0.
        for upper, state in self.transitions:
            line = self.line_data['lines'][state]
            if self.declaration['model'] == 'complete':
                g, dg, tau, a = 0.,0.,None,None
            else:
                if y[upper] > .003*y[self.ground]:
                    raise ValueError('radiation stimulated-emission neglect requires n_upper/(3*n_ground) <= 1e-3')
                g,dg,tau,a = argon_escape(line,density,temperature,self.length,self.declaration['geometryModel'])
            frequency = line['A_s_1']*g
            loss = frequency*y[upper]
            rate[upper] -= loss
            rate[self.ground] += loss
            if jac is not None:
                dn = -density*np.asarray(volume_gradient)/volume
                dn[self.ground] += na/volume
                gradient = y[upper]*line['A_s_1']*dg*dn
                gradient[upper] += frequency
                jac[upper] -= gradient
                jac[self.ground] += gradient
            power = loss*na*6.62607015e-34*299792458./line['wavelength_m']
            photon += power
            records.append(dict(state=state,paschen=line['paschen'],g=g,A_s_1=line['A_s_1'],
                rate_s_1=frequency,optical_depth=tau,damping=a,
                density_m_3=na*y[upper]/volume,loss_mol_s=loss,photon_power_W=power))
        return rate,jac,dict(lines=records,photon_power_W=photon,ground_density_m_3=density)


def has_argon_resonance(species):
    """Whether a molecularly resolved Ar resonance state needs an explicit model."""
    return any(_argon_state(spc) in ('3P1','1P1') for spc in species)
