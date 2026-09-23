#!/usr/bin/env python3
"""
The one argon + wall model every I-246 artifact is built on.

Envelope, parameter values and sources: see ``envelope.md`` in this directory,
which was committed BEFORE any sweep was run. Nothing here is fitted; in
particular no quantity in this file was chosen by looking at an electron density.
"""

import functools
import os

import numpy as np

import rmgpy
import rmgpy.constants as constants
from rmgpy.kinetics import VoronovEIArrhenius
from rmgpy.reaction import Reaction
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.species import Species

# --- provenance ------------------------------------------------------------
WORKTREE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def assert_provenance():
    """One shared editable install serves several worktrees; `import rmgpy`
    resolves by CWD. Refuse to measure a tree other than this one."""
    resolved = os.path.dirname(os.path.dirname(os.path.abspath(rmgpy.__file__)))
    assert resolved == WORKTREE, (
        "WRONG TREE: rmgpy resolved from {0!r}, expected {1!r}".format(resolved, WORKTREE))
    from rmgpy.solver import plasma as m
    return m.__file__


# --- envelope constants (envelope.md section 1) -----------------------------
EV_TO_K = 1.0 / 8.617333262e-5          # K per eV
TORR_TO_PA = 101325.0 / 760.0

TGAS = 298.15                            # K
P_NOMINAL = 5.0 * TORR_TO_PA             # Pa
TE_NOMINAL_EV = 3.0

R_NOMINAL = 0.05                         # m, cylinder radius
L_NOMINAL = 0.30                         # m, cylinder length

# Ellis, McDaniel & Albritton, At. Data Nucl. Data Tables 17 (1976) 177.
MU0_AR_IN_AR = 1.535e-4                  # m^2/(V s) at the Loschmidt density
MU0_AR_IN_AR_REL_UNCERTAINTY = 0.03      # the compilation's own stated accuracy

# Cosmic-ray / ambient background ionisation: 2-10 ion pairs cm^-3 s^-1 at 1 atm,
# scaled by gas density. Carried as the interval it is; NEVER narrowed.
S_EXT_PER_ATM_LOW = 2.0e6                # m^-3 s^-1 at 101325 Pa
S_EXT_PER_ATM_HIGH = 1.0e7               # m^-3 s^-1 at 101325 Pa

VORONOV_YAML = '/home/alon/Code/RMG-database-plasma/input/kinetics/voronov.yaml'


def diffusion_length(radius=R_NOMINAL, length=L_NOMINAL):
    """Lowest diffusion eigenmode of a finite cylinder: 1/L^2 = (2.405/R)^2 + (pi/L)^2."""
    return 1.0 / np.sqrt((2.405 / radius) ** 2 + (np.pi / length) ** 2)


def cosmic_ray_source(pressure_pa, gas_temperature=TGAS, high=False):
    """Background pair-production rate at this gas density, m^-3 s^-1.

    Scaled from the sea-level 1 atm value by number density. The low and high
    ends of the literature interval are both available; neither is preferred,
    and the width propagates into every sub-threshold electron density."""
    n_here = pressure_pa / (constants.kB * gas_temperature)
    n_atm = 101325.0 / (constants.kB * 273.15)
    base = S_EXT_PER_ATM_HIGH if high else S_EXT_PER_ATM_LOW
    return base * n_here / n_atm


@functools.lru_cache(maxsize=1)
def voronov_argon():
    """Argon electron-impact ionisation, from the repository's own coefficients.

    Cached because constructing it re-parses the 23 kB coefficient file and the
    sweep's bisections evaluate k_iz many thousands of times. The returned object
    is only ever read here. Callers that hand it to a Reaction get a fresh one via
    :func:`ionisation_reaction`, so no reactor shares a kinetics object with another."""
    return VoronovEIArrhenius(Z=18, N=18, yaml_path_or_obj=VORONOV_YAML)


@functools.lru_cache(maxsize=None)
def _k_ionisation_cached(te_ev):
    return voronov_argon().get_rate_coefficient(te_ev * EV_TO_K) / constants.Na


def k_ionisation(te_ev):
    """Per-particle <sigma v> for e + Ar -> Ar+ + 2e, m^3/s."""
    return _k_ionisation_cached(float(te_ev))


def nu_ion(te_ev, n_ar):
    """Ionisation frequency per electron, s^-1. Normalised PER ELECTRON:
    nu_ion = k_iz * n_Ar, so that it is directly comparable to nu_wall, which is
    also per charged particle. Comparing a per-electron rate to a per-neutral one
    is the normalisation error a previous audit in this campaign made."""
    return k_ionisation(te_ev) * n_ar


def nu_wall_closed_form(te_ev, n_neutral, lam, mu0=MU0_AR_IN_AR):
    """nu_wall = D_a/Lambda^2 with D_a = mu_i kTe/e and mu_i = mu0 N0/n_neutral.

    An INDEPENDENT re-derivation of what PlasmaReactor.compute_nu_wall computes,
    kept deliberately separate so the two can be compared rather than one being
    used to check itself."""
    from rmgpy.solver.plasma import PLASMA_LOSCHMIDT
    mu_i = mu0 * PLASMA_LOSCHMIDT / n_neutral
    d_a = mu_i * (constants.R / constants.Na) * (te_ev * EV_TO_K) / constants.e
    return d_a / (lam * lam)


def argon_species():
    """Fresh species objects (solver indexing is by identity, so never shared)."""
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    return electron, ar, arp


def ionisation_reaction(electron, ar, arp):
    """e + Ar => Ar+ + e + e. The incident electron is written out explicitly, which
    is what makes the exported equation second order and the m^3/(molecule*s)
    coefficient dimensionally meaningful."""
    return Reaction(
        reactants=[electron, ar], products=[arp, electron, electron], reversible=False,
        # a FRESH kinetics object per reaction: the cached one above is read-only
        # scaffolding for the sweep and must not be shared into a reactor.
        kinetics=VoronovEIArrhenius(Z=18, N=18, yaml_path_or_obj=VORONOV_YAML))


def build_reactor(te_ev=TE_NOMINAL_EV, pressure=P_NOMINAL, tgas=TGAS,
                  x_ion=1.0e-10, radius=R_NOMINAL, length=L_NOMINAL,
                  wall=True, gamma=1.0, source=None, quasineutral=False,
                  with_chemistry=True, max_alpha=None, lam=None,
                  mu0=MU0_AR_IN_AR, termination=None):
    """A fully initialised argon PlasmaReactor, with or without a wall.

    ``x_ion`` is the charge-neutral seed mole fraction of Ar+ (and of e-). It is
    an initial condition only; the steady states reported anywhere in this ticket
    are shown to be independent of it.
    """
    electron, ar, arp = argon_species()
    imf = {electron: x_ion, arp: x_ion, ar: 1.0 - 2.0 * x_ion}
    kwargs = {}
    if wall:
        kwargs['diffusion_length'] = (lam if lam is not None
                                      else diffusion_length(radius, length), 'm')
        kwargs['ion_reduced_mobility'] = (mu0, 'm^2/(V*s)')
        kwargs['wall_recycling'] = gamma
        if source is not None:
            kwargs['ionisation_source'] = (source, 'm^-3/s')
        if max_alpha is not None:
            kwargs['max_ionisation_degree'] = max_alpha
    reactor = PlasmaReactor(
        (tgas, 'K'), (pressure, 'Pa'), imf, (te_ev * EV_TO_K, 'K'),
        n_sims=1, termination=termination or [],
        quasineutral_electron=quasineutral, **kwargs)
    core_species = [electron, ar, arp]
    core_reactions = [ionisation_reaction(electron, ar, arp)] if with_chemistry else []
    reactor.initialize_model(core_species, core_reactions, [], [])
    return reactor, core_species, core_reactions


def state_at(reactor, alpha, n_total_mol=1.0):
    """Packed state at ionisation degree alpha = n_e/n_neutral, charge neutral."""
    y = np.zeros(reactor.num_core_species, float)
    ie = reactor.electron_index
    # y indices follow reactor.species_index; recover them by charge.
    z = reactor.species_charges
    i_ar = [j for j in range(len(z)) if z[j] == 0][0]
    i_arp = [j for j in range(len(z)) if z[j] == 1][0]
    y[i_ar] = n_total_mol
    y[i_arp] = alpha * n_total_mol
    y[ie] = alpha * n_total_mol
    return y


def indices(reactor):
    z = reactor.species_charges
    ie = reactor.electron_index
    i_ar = [j for j in range(len(z)) if z[j] == 0][0]
    i_arp = [j for j in range(len(z)) if z[j] == 1 and j != ie][0]
    return ie, i_ar, i_arp
