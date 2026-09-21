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
Contains the :class:`PlasmaReactor` class, providing a two-temperature
(gas temperature ``Tgas``, electron temperature ``Te``) homogeneous, isobaric
batch reactor with an explicit electron population in the state vector.

The equation of state is

.. math:: P V = R \\left( N_{heavy} T_{gas} + N_e T_e \\right)

where :math:`N_e` is read from the state vector at the electron species index
and :math:`N_{heavy}` is the total moles of all other (heavy) core species.
The same EOS implementation (:meth:`PlasmaReactor.compute_volume`) is used at
initialization and in every residual and Jacobian evaluation.

The reactor constructor takes no electron-density argument and never creates an
electron species on its own. It can be built directly, or from an input file via
the ``plasmaReactor(...)`` directive (see :func:`rmgpy.rmg.input.plasma_reactor`),
which converts an optional ``electronDensity`` into an electron mole fraction on
the driver side -- the constructor's ``initial_mole_fractions`` remains the single
source of the electron amount. The caller must supply exactly one
electron pseudo-species (identified structurally via ``Species.is_electron()``)
among the core species, a strictly positive initial electron amount, and a
strictly positive electron temperature. Any unsupported configuration raises
a named error (:class:`rmgpy.exceptions.PlasmaStateError` or
:class:`rmgpy.exceptions.NonEquilibriumReverseRateError`) before solver
initialization instead of degrading to a one-temperature reactor.
"""

import itertools
import logging

cimport cython
import numpy as np
cimport numpy as np

import rmgpy.constants as constants
cimport rmgpy.constants as constants
from rmgpy.exceptions import NonEquilibriumReverseRateError, PlasmaStateError
from rmgpy.quantity import Quantity
from rmgpy.quantity cimport ScalarQuantity
from rmgpy.solver.base cimport ReactionSystem


# Tolerances for the initial-composition net-charge check in
# PlasmaReactor._warn_if_not_charge_neutral. The check warns when
#
#     |net| > max(PLASMA_NET_CHARGE_ATOL, PLASMA_NET_CHARGE_RTOL * magnitude)
#
# where net = sum(x_i * z_i) per mole of mixture and magnitude = sum(|x_i * z_i|), the
# total charge the net is a cancellation of. BOTH terms are needed. A purely absolute
# test is blind to a weakly ionized deck: an argon plasma seeded at n_e = 1e16 m^-3 and
# 5 torr is 100% charge-imbalanced at a net of only -6.2e-8 per mole, which any
# absolute tolerance loose enough to survive real compositions would swallow. A purely
# relative test cries wolf in the opposite corner, where the charged fraction is itself
# near the roundoff floor and |net|/magnitude is dominated by cancellation error.
#
# ATOL is set above the double-precision accumulation bound for the sum: N terms each
# bounded by 1 accumulate at most ~N * 2.22e-16, so 1e-12 clears models up to ~4500
# species, well beyond any plasma mechanism RMG generates. RTOL at 1e-6 is far above
# that same floor in relative terms and far below any imbalance a modeller could mean.
# Neither is a physical statement: they separate "the arithmetic did not quite close"
# from "the composition is not neutral", nothing more.
PLASMA_NET_CHARGE_ATOL = 1.0e-12
PLASMA_NET_CHARGE_RTOL = 1.0e-6


# --- charged-particle wall transport ---------------------------------------
#
# Reduced (zero-field) ion mobilities in the literature are quoted at the
# Loschmidt constant: the number density of an ideal gas at 273.15 K and
# 101325 Pa. This is that constant (CODATA), and it is a UNIT CONVENTION, not a
# fitted parameter -- it is the density the tabulated mu_0 is normalised to, and
# changing it would mean reading the tabulated number as something it is not.
PLASMA_LOSCHMIDT = 2.6867811e25          # m^-3

# Default ceiling on the ionisation degree n_e / n_neutral above which the
# ion-neutral ambipolar model implemented here is outside its own assumptions
# and the reactor refuses to continue rather than extrapolate.
#
# The model takes the ion's momentum-transfer collisions to be with NEUTRALS, so
# mu_i scales as 1/n_neutral. Coulomb collisions with the other charged particles
# are neglected entirely. The ion-ion Coulomb momentum-transfer cross section
# exceeds the Ar+/Ar resonant-charge-transfer cross section (~5e-19 m^2) once the
# charged fraction rises past roughly 1e-3 at these temperatures, so above that
# the neglected channel is no longer small and the 1/n_neutral scaling is simply
# wrong -- it diverges as n_neutral -> 0, which is a singularity of the
# PARAMETRISATION, not of the physics. This ceiling is a statement about where
# the transport model stops being the right one. It is NOT chosen from, and does
# not depend on, any electron density the model is expected to produce.
PLASMA_WALL_MAX_IONISATION_DEGREE = 1.0e-3


cdef class PlasmaReactor(ReactionSystem):
    """
    A two-temperature plasma reaction system: a homogeneous, isobaric batch
    reactor at constant gas temperature ``Tgas`` and constant electron
    temperature ``Te``, with an explicit electron population carried in the
    state vector. See the module docstring for the equation of state and the
    supported-configuration rules.
    """

    cdef public ScalarQuantity T
    cdef public ScalarQuantity P
    cdef public ScalarQuantity Te
    cdef public double V
    cdef public bint constant_volume
    cdef public dict initial_mole_fractions
    cdef public list const_spc_names
    cdef public list const_spc_indices
    cdef public dict sens_conditions
    cdef public int n_sims

    # Temperature range of a ranged reactor, or None for a scalar reactor. Every
    # other concrete reactor (simple/liquid/surface) declares this; RMG.execute()
    # reads it on every reaction system to compute Tmin/Tmax. PlasmaReactor forbids
    # ranged conditions (see __init__), so this is permanently None -- the driver
    # then falls back to the scalar gas temperature self.T, exactly as it does for a
    # scalar SimpleReactor.
    cdef public list Trange

    # Electron state resolution (set during initialize_model)
    cdef public int electron_index
    cdef public object electron_species

    # True only after initialize_model has validated the composition and the
    # reaction set and resolved the electron state; guards against computing
    # rate coefficients through a base-class initialization that skipped the
    # plasma validation.
    cdef public bint _plasma_validated

    # The LABEL of the chargeBalanceSpecies the deck named (or None). The directive
    # in rmgpy/rmg/input.py assigns this ion a mole fraction so the initial
    # composition comes out neutral; carrying the label here lets the reactor check,
    # once the reaction set exists, that the ion is one the loaded chemistry can
    # actually produce -- otherwise the neutrality is arithmetic only, achieved with
    # a species nothing makes. A string, never a Species object: it survives pickling
    # untouched and is matched against reaction-participant labels at model-init.
    cdef public object charge_balance_species
    # Latched once the reachability of charge_balance_species has been settled (either
    # confirmed producible, or warned about), so the warning is emitted at most once
    # across the many initialize_model calls a model-generation run makes.
    cdef public bint _balance_reachability_settled

    # --- charged-particle wall boundary -------------------------------------
    # All of these are None/0/False on a reactor with no declared wall, and in that
    # state the residual and Jacobian take EXACTLY the paths they took before the
    # wall existed: `has_wall` gates every added term, so a wall-less reactor
    # reproduces the volume-only equations bit for bit, not approximately.
    cdef public bint has_wall
    # Lambda, the characteristic diffusion length (m). Geometry enters here and
    # ONLY here: nu_wall = D_a / Lambda**2. The plasmaReactor(...) directive
    # computes it from a named chamber shape; the reactor takes the length itself,
    # the way SurfaceReactor takes surface_volume_ratio rather than a shape.
    cdef public ScalarQuantity diffusion_length
    # Reduced zero-field ion mobility mu_0, m^2/(V*s), quoted at
    # mobility_reference_density. A measured transport property of the ion in its
    # parent gas; never fitted here.
    cdef public ScalarQuantity ion_reduced_mobility
    cdef public double mobility_reference_density        # m^-3
    # gamma: the fraction of ions neutralised at the wall whose heavy core returns
    # to the gas. 1.0 = fully recycling wall, 0.0 = fully pumping wall.
    cdef public double wall_recycling
    # Volumetric external ion-electron pair source (m^-3 s^-1): the declared
    # physical mechanism that replaces a numerical seed.
    cdef public ScalarQuantity ionisation_source
    cdef public double max_ionisation_degree
    # When True the electron stops being an integrated state and is carried by an
    # algebraic quasineutrality row instead. Independent of has_wall on purpose, so
    # the two can be A/B-ed against each other and the conditioning claim is a
    # measurement rather than an assertion.
    cdef public bint quasineutral_electron

    # --- resolved at initialize_model, all indexed by core species -----------
    cdef public np.ndarray species_charges        # int, net charge z_j
    cdef public np.ndarray neutral_heavy_mask     # uint8, 1 for uncharged non-electrons
    cdef public np.ndarray wall_recycle_target    # int, index of the neutral this
                                                  # cation returns as, or -1
    cdef public np.ndarray source_cation_target   # int, for a neutral, the index of
                                                  # its singly-charged cation, or -1
    # Diagnostics written by every residual evaluation, so a caller can integrate
    # the wall flux INDEPENDENTLY of the gas-phase deficit and compare the two.
    cdef public np.ndarray wall_loss_rates        # mol/s, per core species, signed
    cdef public double nu_wall                    # s^-1, at the last evaluation
    # Purely NUMERICAL floor on the neutral population inside compute_nu_wall, set
    # from the initial composition. It exists so that a Newton TRIAL state with no
    # neutrals left yields a large-but-finite loss frequency that the solver can
    # reject, instead of an infinity that corrupts the iteration matrix. It is not
    # a physical statement and it never binds inside the supported regime: it is
    # reached only at ionisation degrees a million times the validity ceiling.
    cdef public double wall_neutral_floor
    # Multiplier on the algebraic quasineutrality row. Scaling one row of a DAE
    # residual by a positive constant is an EXACT operation -- it changes no
    # solution, only the matrix that gets factored -- and it is necessary here.
    # The charge row's entries are the charges themselves, exactly +1 and -1,
    # while every differential row carries entries of order of the rates, 1e2 to
    # 1e5 s^-1. DASPK factors the iteration matrix with partial pivoting and no
    # equilibration, so that mismatch is a real cost: measured at alpha = 1e-8,
    # cond_2 of the iteration matrix is 5.39e5 with the electron integrated and
    # 9.41e7 with an unscaled charge row -- a 175x PENALTY from the substitution
    # that was supposed to improve the conditioning. Scaling the row to match the
    # rest of the matrix removes all of it (5.08e5, flat over three decades of
    # the multiplier). Set once, at initialization, so the residual and the
    # Jacobian can never disagree about it mid-run.
    cdef public double charge_row_scale

    def __init__(self, T, P, initial_mole_fractions, Te, n_sims=1, termination=None, sensitive_species=None,
                 sensitivity_threshold=1e-3, sens_conditions=None, const_spc_names=None,
                 charge_balance_species=None,
                 diffusion_length=None, ion_reduced_mobility=None,
                 mobility_reference_density=None, wall_recycling=1.0,
                 ionisation_source=None,
                 max_ionisation_degree=PLASMA_WALL_MAX_IONISATION_DEGREE,
                 quasineutral_electron=False):
        ReactionSystem.__init__(self, termination, sensitive_species, sensitivity_threshold)

        if isinstance(T, list) or isinstance(P, list) or isinstance(Te, list):
            raise PlasmaStateError(
                "PlasmaReactor supports only scalar T, P and Te; ranged reactor "
                "conditions are not supported in this reactor. Got T={0!r}, P={1!r}, Te={2!r}.".format(T, P, Te))
        if Te is None:
            raise PlasmaStateError(
                "PlasmaReactor requires an explicit electron temperature Te; got None. "
                "There is no fallback from Te to Tgas.")

        self.T = Quantity(T)
        self.P = Quantity(P)
        self.Te = Quantity(Te)

        for name, quantity_obj in (('Tgas', self.T), ('P', self.P), ('Te', self.Te)):
            value = quantity_obj.value_si
            if not np.isfinite(value) or value <= 0.0:
                raise PlasmaStateError(
                    "PlasmaReactor requires finite, strictly positive {0}; got {1!r} "
                    "in {2}.".format(name, value, self._identity()))

        self.initial_mole_fractions = initial_mole_fractions

        if const_spc_names:
            raise PlasmaStateError(
                "PlasmaReactor does not support constant species in this slice: "
                "freezing a species (in particular the electron) would silently "
                "break charge conservation and residual/Jacobian consistency. "
                "Got const_spc_names={0!r} in {1}.".format(const_spc_names, self._identity()))
        self.const_spc_indices = None
        self.const_spc_names = None

        self.V = 0  # will be set in initialize_model
        self.constant_volume = False

        self.electron_index = -1
        self.electron_species = None
        self._plasma_validated = False

        self.sens_conditions = sens_conditions
        self.n_sims = n_sims

        self.charge_balance_species = charge_balance_species
        self._balance_reachability_settled = False

        self._configure_wall(diffusion_length, ion_reduced_mobility,
                             mobility_reference_density, wall_recycling,
                             ionisation_source, max_ionisation_degree,
                             quasineutral_electron)

    def _configure_wall(self, diffusion_length, ion_reduced_mobility,
                        mobility_reference_density, wall_recycling,
                        ionisation_source, max_ionisation_degree,
                        quasineutral_electron):
        """
        Validate and store the charged-particle wall boundary parameters.

        A wall exists only when BOTH a diffusion length and an ion mobility are
        given: nu_wall = D_a / Lambda**2 needs both, and half of the pair is a
        deck error, not a default to be filled in. Every parameter is a declared
        geometric or transport quantity; none of them is permitted to be chosen
        by looking at an electron density.
        """
        self.wall_loss_rates = None
        self.nu_wall = 0.0
        self.species_charges = None
        self.neutral_heavy_mask = None
        self.wall_recycle_target = None
        self.source_cation_target = None

        if (diffusion_length is None) != (ion_reduced_mobility is None):
            raise PlasmaStateError(
                "a charged-particle wall needs BOTH a diffusion length and an ion "
                "reduced mobility -- nu_wall = D_a / Lambda**2 is not defined "
                "without each of them -- but only one was supplied "
                "(diffusion_length={0!r}, ion_reduced_mobility={1!r}). Supply both "
                "to declare a wall, or neither for the volume-only reactor. "
                "({2})".format(diffusion_length, ion_reduced_mobility, self._identity()))

        self.has_wall = diffusion_length is not None
        self.quasineutral_electron = bool(quasineutral_electron)

        if self.has_wall:
            self.diffusion_length = Quantity(diffusion_length)
            self.ion_reduced_mobility = Quantity(ion_reduced_mobility)
            for name, quantity_obj in (('diffusion_length', self.diffusion_length),
                                       ('ion_reduced_mobility', self.ion_reduced_mobility)):
                value = quantity_obj.value_si
                if not np.isfinite(value) or value <= 0.0:
                    raise PlasmaStateError(
                        "PlasmaReactor requires a finite, strictly positive {0}; got "
                        "{1!r}. A zero diffusion length is an infinite loss frequency "
                        "and a zero mobility is no transport at all; neither is a wall. "
                        "({2})".format(name, value, self._identity()))
        else:
            self.diffusion_length = None
            self.ion_reduced_mobility = None

        if mobility_reference_density is None:
            self.mobility_reference_density = PLASMA_LOSCHMIDT
        else:
            self.mobility_reference_density = float(Quantity(mobility_reference_density).value_si
                                                    if not isinstance(mobility_reference_density, float)
                                                    else mobility_reference_density)
        if not np.isfinite(self.mobility_reference_density) or self.mobility_reference_density <= 0.0:
            raise PlasmaStateError(
                "mobility_reference_density must be a finite, strictly positive number "
                "density (m^-3); got {0!r}. ({1})".format(
                    self.mobility_reference_density, self._identity()))

        self.wall_recycling = float(wall_recycling)
        if not np.isfinite(self.wall_recycling) or not (0.0 <= self.wall_recycling <= 1.0):
            raise PlasmaStateError(
                "wall_recycling (gamma) is the fraction of wall-neutralised ions whose "
                "heavy core returns to the gas and must lie in [0, 1]; got {0!r}. "
                "({1})".format(self.wall_recycling, self._identity()))

        if ionisation_source is None:
            self.ionisation_source = Quantity((0.0, 'm^-3/s'))
        else:
            self.ionisation_source = Quantity(ionisation_source)
        if not np.isfinite(self.ionisation_source.value_si) or self.ionisation_source.value_si < 0.0:
            raise PlasmaStateError(
                "ionisation_source is a volumetric production rate of ion-electron "
                "pairs and must be finite and non-negative (m^-3 s^-1); got {0!r}. "
                "({1})".format(self.ionisation_source.value_si, self._identity()))

        # Replaced with a commensurate value in initialize_model once the rate
        # coefficients exist; 1.0 until then, which is exact but badly scaled.
        self.charge_row_scale = 1.0

        self.max_ionisation_degree = float(max_ionisation_degree)
        if not np.isfinite(self.max_ionisation_degree) or self.max_ionisation_degree <= 0.0:
            raise PlasmaStateError(
                "max_ionisation_degree bounds the regime in which the ion-neutral "
                "ambipolar transport model is valid and must be finite and strictly "
                "positive; got {0!r}. ({1})".format(
                    self.max_ionisation_degree, self._identity()))

    def convert_initial_keys_to_species_objects(self, species_dict):
        """
        Convert the ``initial_mole_fractions`` dictionary from species labels into
        species objects, using the given ``species_dict`` (label -> Species).

        The input reader (``read_input_file``) calls this unconditionally on every
        non-``Reactor`` reaction system after the input file is parsed, once the
        species objects exist. ``plasma_reactor`` builds ``initial_mole_fractions``
        keyed by label (including the electron entry it may have inserted from an
        ``electronDensity`` directive), so this is where those labels become the
        Species objects that ``initialize_model`` validates and packs.
        """
        initial_mole_fractions = {}
        for label, mole_frac in self.initial_mole_fractions.items():
            initial_mole_fractions[species_dict[label]] = mole_frac
        self.initial_mole_fractions = initial_mole_fractions

    def _identity(self):
        """A short identity string for this reactor, used in error messages."""
        return "PlasmaReactor(label={0!r}, T={1!r}, P={2!r}, Te={3!r})".format(
            getattr(self, 'label', None), self.T, self.P, self.Te)

    def __reduce__(self):
        """
        A helper function used when pickling an object.

        Every parameter ``__init__`` accepts has to appear here, including the
        wall parameters. This is a hand-written enumeration of a signature, so
        it does not fail loudly when it falls behind: a missing entry silently
        yields a copy that reconstructs, initialises and integrates perfectly
        while having quietly lost whatever was omitted. The wall parameters were
        omitted exactly that way once already, and ``deepcopy`` inherited it
        through ``__reduce_ex__``. ``test_reduce_enumerates_every_constructor_
        parameter`` in test/rmgpy/solver/plasmaWallTest.py fails if this list
        falls behind ``__init__`` again.

        The stored ``Quantity`` objects are passed straight back:
        ``Quantity(Quantity(x))`` is idempotent in both value and units, so the
        copy reconstructs them exactly rather than round-tripping unit strings.
        ``mobility_reference_density`` is a plain float here, which is the
        branch of ``_configure_wall`` that preserves it as given rather than
        re-deriving it.
        """
        return (self.__class__,
                (self.T, self.P, self.initial_mole_fractions, self.Te, self.n_sims, self.termination,
                 self.sensitive_species, self.sensitivity_threshold, self.sens_conditions,
                 self.const_spc_names, self.charge_balance_species,
                 self.diffusion_length, self.ion_reduced_mobility,
                 self.mobility_reference_density, self.wall_recycling,
                 self.ionisation_source, self.max_ionisation_degree,
                 self.quasineutral_electron))

    cpdef initialize_model(self, list core_species, list core_reactions, list edge_species, list edge_reactions,
                          list surface_species=None, list surface_reactions=None, list pdep_networks=None,
                          atol=1e-16, rtol=1e-8, sensitivity=False, sens_atol=1e-6, sens_rtol=1e-4,
                          filter_reactions=False, dict conditions=None):
        """
        Initialize a simulation of the plasma reactor using the provided
        kinetic model. All electron-state and reverse-rate-policy validation
        happens here, before the first residual or Jacobian evaluation.
        """
        # Unsupported features fail loudly instead of being silently ignored.
        if surface_species or surface_reactions:
            raise PlasmaStateError(
                "PlasmaReactor does not support surface species or surface reactions; "
                "got {0} surface species and {1} surface reactions in {2}.".format(
                    len(surface_species or []), len(surface_reactions or []), self._identity()))
        if pdep_networks:
            raise PlasmaStateError(
                "PlasmaReactor does not support pressure-dependent networks; got {0} "
                "networks in {1}.".format(len(pdep_networks), self._identity()))
        if sensitivity or self.sensitive_species:
            raise PlasmaStateError(
                "PlasmaReactor does not support sensitivity analysis in this slice; "
                "requested in {0}.".format(self._identity()))
        if filter_reactions:
            raise PlasmaStateError(
                "PlasmaReactor does not support reaction filtering: the filtering "
                "thresholds are single-temperature quantities and are undefined for "
                "electron-impact kinetics. Requested in {0}.".format(self._identity()))
        if conditions:
            raise PlasmaStateError(
                "PlasmaReactor does not accept externally supplied ranged/override "
                "conditions {0!r}; construct the reactor directly with scalar T, P "
                "and Te. ({1})".format(conditions, self._identity()))

        # THE electron-representation boundary. This is the single production
        # call site of the family-declared electron-placement resolver: the
        # earliest point at which family-generation semantics are complete and
        # the target reactor is known, but before any stoichiometric index or
        # concentration product is built. Reactions that still carry their
        # electron as metadata (reaction.electrons != 0) are replaced, once
        # each, by the resolver's non-mutating reactor-facing view; the
        # canonical reaction objects are never touched. Because the local
        # reaction lists are rebound here, the resolved representation is the
        # ONLY one seen downstream -- by validation, the base-class packing,
        # rate-coefficient generation, residual and Jacobian evaluation --
        # so there is never a second, independently derived electron
        # representation. The canonical (model) reaction lists are retained so
        # the public reaction_index map can be re-keyed back to them at the end,
        # keeping solver-internal views out of the caller-facing identity.
        core_reactions_model = core_reactions
        edge_reactions_model = edge_reactions
        core_reactions = self._resolve_electron_placements(core_reactions, core_species, edge_species)
        edge_reactions = self._resolve_electron_placements(edge_reactions, core_species, edge_species)

        # Validate the composition and the reaction set BEFORE the base class
        # packs its stoichiometry arrays: malformed input would otherwise
        # surface as unnamed IndexError/KeyError from the packing code instead
        # of the named errors this reactor promises.
        self._validate_electron_state(core_species, edge_species)
        self._validate_reactions(core_species, edge_species, core_reactions, edge_reactions)

        # If the deck balanced its charge with a named ion, check -- now that the
        # reaction set exists -- that the loaded chemistry can actually produce that
        # ion. Warns, never raises: an unproducible balancing ion means the neutrality
        # this reactor was handed is arithmetic only, but that is the modeller's to
        # judge, not the reactor's to forbid.
        self._warn_if_balance_ion_unreachable(core_reactions, edge_reactions)

        # Now call the base class version of the method.
        # This initializes the attributes declared in the base class,
        # including the species/reaction index maps and stoichiometry arrays.
        ReactionSystem.initialize_model(self, core_species=core_species, core_reactions=core_reactions,
                                       edge_species=edge_species, edge_reactions=edge_reactions,
                                       surface_species=[], surface_reactions=[],
                                       pdep_networks=None, atol=atol, rtol=rtol, sensitivity=False,
                                       sens_atol=sens_atol, sens_rtol=sens_rtol, filter_reactions=False,
                                       conditions=None)

        # Resolve the electron's position in the packed solver state.
        self.electron_index = self.species_index[self.electron_species]
        if not (0 <= self.electron_index < self.num_core_species):
            raise PlasmaStateError(
                "electron species resolved to index {0}, outside the packed core "
                "state of size {1}. ({2})".format(
                    self.electron_index, self.num_core_species, self._identity()))
        self._plasma_validated = True

        # Resolve the charge bookkeeping the wall operator and the quasineutrality
        # row both read. Runs unconditionally (it is cheap and its results are used
        # by diagnostics either way) but only RAISES when a wall or the algebraic
        # electron is actually configured.
        self._resolve_wall_state(core_species)

        # Set initial conditions (also re-checks the packed electron amount)
        self.set_initial_conditions()

        ReactionSystem.compute_network_variables(self, None)

        # Generate forward and reverse rate coefficients using the accepted
        # reaction set and the two-temperature evaluation policy.
        self.generate_rate_coefficients(core_reactions, edge_reactions)

        # Make the algebraic charge row commensurate with the differential rows,
        # now that the rate coefficients exist and the Jacobian can be evaluated.
        # See charge_row_scale: this is an exact row scaling of the DAE residual
        # and it removes a measured 175x conditioning penalty. Taken from the
        # DIFFERENTIAL rows only, and once, so residual and Jacobian can never
        # disagree about it mid-run.
        if self.quasineutral_electron:
            self._set_charge_row_scale()

        ReactionSystem.set_initial_derivative(self)
        # Initialize the model
        ReactionSystem.initialize_solver(self)

        # Re-key the public reaction map back to the canonical (model) reactions.
        # The base class built reaction_index from the resolved views (needed by
        # generate_rate_coefficients above, which looks reactions up in it); but
        # a caller -- e.g. the model builder's remove_species, which iterates
        # reaction_system.reaction_index to find reactions touching a removed
        # species -- must see the SAME reaction objects the model holds, not
        # reactor-internal views it never created. Canonical and view lists are
        # positionally aligned, so the solver's positional indices are unchanged;
        # only the dictionary keys are restored to model identity. This runs
        # after every consumer that needs view-keyed lookups, and after
        # initialize_solver, so runtime residual/Jacobian (which index the packed
        # arrays positionally, never through reaction_index) are unaffected.
        self._rekey_reaction_index_to_model(core_reactions_model, edge_reactions_model)

    def _rekey_reaction_index_to_model(self, core_reactions, edge_reactions):
        """
        Rebuild ``self.reaction_index`` so its keys are the canonical (model)
        reactions passed to :meth:`initialize_model`, at the same positional
        indices the base class assigned to the resolved views. Reassigns the
        dict, so any views the base class registered are discarded -- which also
        means re-initializing the same reactor does not accumulate stale view
        keys, because the canonical reactions are stable across calls.
        """
        self.reaction_index = {}
        for index, rxn in enumerate(itertools.chain(core_reactions or [], edge_reactions or [])):
            self.reaction_index[rxn] = index

    def _resolve_electron_placements(self, reactions, core_species, edge_species):
        """
        Return a reactor-facing view of ``reactions`` in which every reaction
        that still carries its electron as metadata (``reaction.electrons`` is
        nonzero) has been replaced by the resolver's non-mutating placement
        view, and every already-in-reactor-form reaction (``electrons == 0``)
        is passed through unchanged, by identity.

        The order in which the electron is placed is NOT derived here from the
        net electron count: the gate below only SELECTS a reaction for
        resolution; the side and count come from the family-level declaration
        inside :func:`rmgpy.electron_placement.resolve_electron_placement`,
        which hard-fails by name (``ElectronPlacementError``) for any family,
        shape, direction, or kinetics it does not authorise -- including the
        ionisation-shaped case, whose net electron production does not
        determine incident-electron order. The resolver is looked up on the
        module at call time so a test can spy on it; the canonical reaction
        objects handed in are never mutated.

        The gate is :meth:`_needs_electron_placement`. It admits two clauses: a
        metadata-electron reaction (``electrons != 0``), which is the historical
        case; AND a reaction whose owner declares a placement, whose
        ``electrons`` scalar is 0, and which carries no explicit electron yet.
        The MOTIVATING case for the second clause is the CONSERVED-electron
        reaction -- ``AB + e- -> A + B + e-`` translated to ``AB -> A + B``, net
        zero -- which was skipped until I-208 because a net of zero closed the
        old ``if getattr(rxn, 'electrons', 0)`` gate; the resolver was never
        consulted and the reaction kept its heavy-only molecularity while its
        coefficient was second order -- a first-order evaluation of a
        second-order rate, wrong by a factor of the electron density. But net
        zero is the MOTIVATION for the clause, not its boundary: the clause keys
        on ``electrons == 0`` for ANY declared owner, so it ALSO admits a
        declared NONZERO owner (e.g. ``Plasma_Electron_Attachment`` at
        ``(1, 0)``) that arrives with a defaulted or corrupt zero ``electrons``
        scalar and no explicit electron. Such a reaction now hard-fails BY NAME
        at the resolver's net-mismatch check rather than passing through
        silently as the old one-line gate let it -- that loud failure is
        intentional, and :meth:`_needs_electron_placement` documents why the
        clause is deliberately not narrowed to exclude it. (On the full reactor
        path with plasma-flagged kinetics a wrong incident order does not reach
        the solver: ``_validate_reactions`` refuses a missing incident electron
        and RAISES. The genuinely silent form of this hazard is the EXPORT path
        with generic Arrhenius-like kinetics, whose order the writers' guard
        cannot read; that lives outside this gate and is recorded in the I-208
        report.) Everything else (an ordinary neutral with no electron and no
        declaration; a reaction already in reactor form with its electron
        explicit) is passed through by identity: a neutral must NOT reach the
        resolver, or it would fail by name.
        """
        from rmgpy import electron_placement

        if reactions is None:
            return reactions
        species_list = list(core_species or []) + list(edge_species or [])
        resolved = []
        for rxn in reactions:
            if self._needs_electron_placement(rxn):
                resolved.append(
                    electron_placement.resolve_electron_placement(rxn, species_list))
            else:
                resolved.append(rxn)
        return resolved

    def _needs_electron_placement(self, rxn):
        """True if `rxn` must be handed to the electron-placement resolver.

        The predicate that widened at I-208. A metadata-electron reaction
        (``electrons != 0``) always goes -- unchanged from before, so an
        undeclared owner still fails BY NAME inside the resolver rather than
        passing silently through, and a double representation is still caught
        there. A CONSERVED-electron reaction (``electrons == 0``) goes only when
        its owner declares a placement AND it does not already carry an explicit
        electron; that is the ``(reactant_count, product_count)`` net-zero shape
        the old ``electrons != 0`` gate dropped. Resolved on the same
        :func:`rmgpy.electron_balance.get_placement_owner` the resolver keys on,
        so the gate and the resolver cannot disagree about which owner a
        reaction belongs to.

        This clause admits ``electrons == 0`` for ANY declared owner, not only
        owners whose declared net is zero -- and that width is INTENTIONAL. A
        declared owner with a NONZERO net (e.g. ``Plasma_Electron_Attachment`` at
        ``(1, 0)``) that arrives with its ``electrons`` scalar missing or
        defaulted to 0 and no explicit electron is a CORRUPT reaction; it now
        routes to the resolver and hard-fails BY NAME at the net-mismatch check,
        where the old one-line gate passed it through silently. Loud failure of
        that corruption is the point. Narrowing the clause with a
        ``product_count == reactant_count`` condition would restore the silent
        pass-through for exactly the corrupt input this is meant to catch, so it
        is deliberately absent.

        Everything else is passed through: an ordinary neutral (no electron
        anywhere, no declaration) must not reach the resolver or it would fail
        by name, and a reaction already carrying its electron explicitly is
        already in reactor form with nothing to place.
        """
        from rmgpy import electron_placement
        from rmgpy.electron_balance import get_placement_owner

        if getattr(rxn, 'electrons', 0):
            return True
        owner = get_placement_owner(rxn)
        if owner is None or owner not in electron_placement.FAMILY_ELECTRON_PLACEMENT:
            return False
        # Use the resolver's own defensive electron test (catches AttributeError/
        # IndexError on a malformed participant) so the gate and the resolver
        # agree on what an electron is, and a malformed reaction fails INSIDE the
        # resolver with the named ElectronPlacementError rather than at the gate
        # with a raw error.
        for spc in itertools.chain(rxn.reactants or [], rxn.products or []):
            if electron_placement._is_electron(spc):
                return False  # already in reactor form; nothing to place
        return True

    def _validate_electron_state(self, list core_species, list edge_species):
        """
        Locate the single electron pseudo-species structurally (via
        ``Species.is_electron()``) and validate the initial electron amount.
        Runs before the base class packs any array. Raises
        :class:`PlasmaStateError` on any degenerate configuration; never
        falls back to a one-temperature description.
        """
        cdef list electron_core, electron_edge
        cdef set core_set

        for spc in itertools.chain(core_species, edge_species, self.initial_mole_fractions.keys()):
            if not callable(getattr(spc, 'is_electron', None)):
                raise PlasmaStateError(
                    "PlasmaReactor requires Species objects (direct construction); got "
                    "{0!r} in {1}.".format(spc, self._identity()))

        electron_core = [spc for spc in core_species if spc.is_electron()]
        electron_edge = [spc for spc in edge_species if spc.is_electron()]

        if len(electron_core) == 0:
            raise PlasmaStateError(
                "no electron pseudo-species is present in the core species; a "
                "zero-electron plasma state is a domain failure, not an implicit "
                "thermal reactor. Use SimpleReactor for ordinary thermal "
                "simulations. ({0})".format(self._identity()))
        if len(electron_core) > 1:
            raise PlasmaStateError(
                "found {0} electron pseudo-species in the core species; the electron "
                "state must be locatable unambiguously in the packed solver state. "
                "({1})".format(len(electron_core), self._identity()))
        if electron_edge:
            raise PlasmaStateError(
                "found {0} electron pseudo-species among the edge species; the "
                "electron state must be unique across the model. ({1})".format(
                    len(electron_edge), self._identity()))

        self.electron_species = electron_core[0]

        # Every initial-composition key must be a core species of this model,
        # and any electron-like key must be THE resolved electron species --
        # otherwise the electron amount could silently belong to an object
        # that is not in the packed state.
        core_set = set(core_species)
        electron_amount = None
        for spc, value in self.initial_mole_fractions.items():
            if spc not in core_set:
                raise PlasmaStateError(
                    "initial_mole_fractions key {0!r} is not a core species of "
                    "this model. ({1})".format(spc, self._identity()))
            if spc.is_electron():
                if spc is not self.electron_species:
                    raise PlasmaStateError(
                        "initial_mole_fractions contains an electron-like key "
                        "{0!r} that is not the model's electron species. "
                        "({1})".format(spc, self._identity()))
                if electron_amount is not None:
                    raise PlasmaStateError(
                        "multiple electron entries in initial_mole_fractions; the "
                        "electron amount must be single-valued. ({0})".format(self._identity()))
                electron_amount = float(value)
        if electron_amount is None:
            raise PlasmaStateError(
                "no initial electron amount was supplied in initial_mole_fractions; "
                "an explicit, finite, strictly positive electron amount is required. "
                "({0})".format(self._identity()))
        if not np.isfinite(electron_amount) or electron_amount <= 0.0:
            raise PlasmaStateError(
                "initial electron amount must be finite and strictly positive; got "
                "{0!r}. ({1})".format(electron_amount, self._identity()))

        # Fifth check, ADDED alongside the four above and replacing none of them: is
        # the initial composition charge neutral? This one WARNS and never raises. A
        # deliberately non-neutral initial condition may be perfectly legitimate and
        # this reactor does not get to decide that; what it refuses to do is stay
        # silent, which is what it did before -- a deck seeded with an electron density
        # and no compensating cation ran to completion without one message naming the
        # charge it was carrying.
        self._warn_if_not_charge_neutral()

    def _warn_if_not_charge_neutral(self):
        """
        Log a warning naming the initial composition's net charge per mole when it is
        not neutral to within
        ``max(PLASMA_NET_CHARGE_ATOL, PLASMA_NET_CHARGE_RTOL * magnitude)``.

        Never raises: non-neutrality is reported, not forbidden. A species whose charge
        cannot be determined is named in the message rather than skipped silently, so
        the reported number is never quietly incomplete.
        """
        cdef double net, magnitude, total, threshold

        net = 0.0
        magnitude = 0.0
        total = 0.0
        undetermined = []
        contributions = []
        for spc, value in self.initial_mole_fractions.items():
            amount = float(value)
            total += amount
            try:
                charge = int(spc.get_net_charge())
            except Exception:
                undetermined.append(str(spc))
                continue
            if charge:
                net += amount * charge
                magnitude += abs(amount * charge)
                contributions.append((abs(amount * charge), str(spc), charge, amount))

        # Report per mole of mixture. A directly-constructed reactor need not have been
        # handed a composition summing to one, and "net charge" is only comparable to a
        # tolerance once it is per mole of something.
        if np.isfinite(total) and total > 0.0:
            net /= total
            magnitude /= total

        if undetermined:
            logging.warning(
                "PlasmaReactor could not determine the net charge of %s in the initial "
                "composition, so the charge-neutrality check below is computed without "
                "them and may be incomplete. (%s)",
                ', '.join(sorted(undetermined)), self._identity())

        if not np.isfinite(net):
            logging.warning(
                "PlasmaReactor could not evaluate the initial composition's net charge "
                "(got %r). (%s)", net, self._identity())
            return

        threshold = max(PLASMA_NET_CHARGE_ATOL, PLASMA_NET_CHARGE_RTOL * magnitude)
        if abs(net) <= threshold:
            return

        contributions.sort(reverse=True)
        breakdown = ', '.join(
            '{0} (charge {1:+d}) x={2!r}'.format(label, charge, amount)
            for _, label, charge, amount in contributions[:10])
        if len(contributions) > 10:
            breakdown += ', ... {0} more charged species'.format(len(contributions) - 10)

        logging.warning(
            "PlasmaReactor initial composition is NOT charge neutral: net charge = %r "
            "per mole (total charge magnitude %r, relative imbalance %r), above the "
            "tolerance max(atol=%r, rtol=%r * magnitude) = %r. Charged species: %s. "
            "This is a warning, not an error -- a deliberately non-neutral initial "
            "condition is legitimate -- but a deck that seeds an electronDensity "
            "without a compensating cation lands here by accident. To have RMG compute "
            "the balancing ion's mole fraction for you, name it with "
            "chargeBalanceSpecies='<label>' in the plasmaReactor(...) block. (%s)",
            net, magnitude, (abs(net) / magnitude if magnitude > 0.0 else float('inf')),
            PLASMA_NET_CHARGE_ATOL, PLASMA_NET_CHARGE_RTOL, threshold, breakdown,
            self._identity())

    def _warn_if_balance_ion_unreachable(self, core_reactions, edge_reactions):
        """
        Warn -- once, and never raise -- when the deck balanced its charge with a
        ``chargeBalanceSpecies`` that no reaction in the loaded model produces.

        THE DEFINITION OF "reachable" this uses, stated so it can be argued with: the
        balancing ion is reachable when its label matches a species that appears as a
        PRODUCT of at least one reaction (core or edge), or -- for a reversible
        reaction, which runs in both directions -- as a participant on either side.
        A species that is merely declared, or merely consumed by irreversible
        reactions, is NOT reachable: nothing makes it. That is exactly the phantom the
        check is for -- an ion given a mole fraction to zero the net charge while no
        chemistry can generate it, so the neutral state is one the mechanism cannot
        sustain.

        WHY HERE AND NOT AT PARSE TIME. Reachability is undecidable when the input file
        is read: ``database()`` records only the NAMES of reaction libraries and seed
        mechanisms, never their contents, and the kinetics families have not run.
        ``initialize_model`` is the first point the reaction set exists, so it is the
        first point the question can be answered at all.

        WHAT THIS CANNOT CATCH, stated so the warning is not mistaken for more than it
        is:
          * It is topological, not quantitative. A production reaction that exists but
            carries negligible rate counts as reachable; the ion has a path, even if
            no meaningful flux. A deck can still be ill-posed for want of flux.
          * It is a one-hop presence check, not a graph traversal from the initial
            species. An ion produced only from other species that are themselves never
            formed still counts as reachable.
          * It matches by label. An ion a library supplies under a different label than
            the deck declared would read as unreachable (a false alarm) -- which is one
            more reason this warns rather than refuses.
          * It is settled at the first model initialization that resolves it. A
            producing reaction generated by a family only in a later enlargement is not
            retroactively credited; for a genuinely unproducible ion, which is the
            target case, this makes no difference.

        Because of these, the check restores a SIGNAL -- it names a suspect. It does not
        prove a fault, and it adds no chemistry: the gap it points at is real and stays
        unfilled.
        """
        cdef bint reachable
        cdef object label, spc, rxn

        if self.charge_balance_species is None or self._balance_reachability_settled:
            return
        label = self.charge_balance_species

        reachable = False
        for rxn in itertools.chain(core_reactions or [], edge_reactions or []):
            participants = list(rxn.products)
            if getattr(rxn, 'reversible', False):
                participants = participants + list(rxn.reactants)
            for spc in participants:
                if getattr(spc, 'label', None) == label:
                    reachable = True
                    break
            if reachable:
                break

        # Settle it either way so the (loud) warning cannot repeat across the many
        # initialize_model calls a generation run makes.
        self._balance_reachability_settled = True
        if reachable:
            return

        logging.warning(
            "PlasmaReactor: chargeBalanceSpecies %r was assigned a mole fraction to "
            "make the initial composition charge neutral, but no reaction in the "
            "loaded model produces it -- it is a product of no reaction, and a "
            "participant in no reversible one, among the %d core and %d edge "
            "reactions. The charge balance therefore rests on an ion the loaded "
            "chemistry (reaction libraries, seed mechanisms, kinetics families) cannot "
            "generate: the composition is neutral only by arithmetic, because a mole "
            "fraction was given to a species nothing makes, so the neutral state it "
            "produces is not one this mechanism can sustain. This is the "
            "charge-neutrality signal that balancing with %r silenced. Check that a "
            "library, seed mechanism or family able to produce %r is loaded, or remove "
            "the chargeBalanceSpecies directive. This is a warning, not an error, and "
            "it names a suspect rather than proving a fault. (%s)",
            label, len(core_reactions or []), len(edge_reactions or []),
            label, label, self._identity())

    def _validate_reactions(self, core_species, edge_species, core_reactions, edge_reactions):
        """
        Validate every reaction in the model against the electron
        representation rules and the reverse-rate policy. Runs before the base
        class packs any array, so malformed reactions raise named errors here
        rather than unnamed packing crashes. The result does not depend on the
        instantaneous numerical state: in particular, ``Te == Tgas`` does not
        bypass any check.
        """
        cdef set core_set, model_set
        core_set = set(core_species)
        model_set = core_set | set(edge_species)
        num_core_reactions = len(core_reactions)
        for position, rxn in enumerate(itertools.chain(core_reactions, edge_reactions)):
            kin = rxn.kinetics
            if kin is None:
                raise PlasmaStateError(
                    "reaction {0!s} has no kinetics. ({1})".format(rxn, self._identity()))

            # The packed stoichiometry arrays hold up to three participants per
            # side and encode "no species" as -1; anything outside 1..3 would be
            # read as a valid index and silently corrupt a neighbouring row.
            if not (1 <= len(rxn.reactants) <= 3) or not (1 <= len(rxn.products) <= 3):
                raise PlasmaStateError(
                    "reaction {0!s} has {1} reactants and {2} products; between "
                    "one and three of each are supported. ({3})".format(
                        rxn, len(rxn.reactants), len(rxn.products), self._identity()))

            # Every participant must belong to the model, and core reactions
            # must touch only core species: the residual/Jacobian write core
            # rows without bounds checks on that assumption.
            for spc in itertools.chain(rxn.reactants, rxn.products):
                if spc not in model_set:
                    raise PlasmaStateError(
                        "reaction {0!s} involves species {1!r} which is not in "
                        "the model. ({2})".format(rxn, spc, self._identity()))
                if position < num_core_reactions and spc not in core_set:
                    raise PlasmaStateError(
                        "core reaction {0!s} involves non-core species {1!r}. "
                        "({2})".format(rxn, spc, self._identity()))

            if kin.is_pressure_dependent():
                raise PlasmaStateError(
                    "reaction {0!s} has pressure-dependent kinetics {1}, which this "
                    "reactor does not support. ({2})".format(
                        rxn, kin.__class__.__name__, self._identity()))

            # The electron must participate through explicit reactant/product
            # entries, where its incident order and net stoichiometry are both
            # unambiguous. A nonzero metadata-only electron count
            # (reaction.electrons) carries only the net change and cannot
            # express incident-electron order, so it is rejected rather than
            # guessed at.
            if getattr(rxn, 'electrons', 0) != 0:
                raise PlasmaStateError(
                    "reaction {0!s} carries a metadata-only electron count "
                    "(electrons={1}); this representation cannot distinguish "
                    "incident-electron order from net electron production, so it is "
                    "unsupported here. Represent the electron explicitly as a "
                    "reactant/product species instead. ({2})".format(
                        rxn, rxn.electrons, self._identity()))

            n_electron_reactants = sum(1 for spc in rxn.reactants if spc.is_electron())
            n_electron_products = sum(1 for spc in rxn.products if spc.is_electron())

            # A rate law that declares itself electron-density-driven requires
            # an explicit incident electron among the reactants; without one the
            # rate would silently lose its dependence on the electron
            # population.
            if getattr(kin, 'uses_electron_density', False) and n_electron_reactants == 0:
                raise PlasmaStateError(
                    "reaction {0!s} has kinetics {1} which declares "
                    "uses_electron_density, but no explicit electron appears "
                    "among its reactants; the incident electron must be "
                    "represented explicitly. ({2})".format(
                        rxn, kin.__class__.__name__, self._identity()))

            if getattr(kin, 'uses_electron_temperature', False):
                # The kinetics must expose an unambiguous electron-temperature
                # evaluator; the standard k(T, P) interface silently collapses
                # to a one-temperature evaluation and is never used here.
                if (not hasattr(kin, 'get_rate_coefficient_two_temp')
                        and not hasattr(kin, 'get_rate_coefficient_electron_temp')):
                    raise PlasmaStateError(
                        "reaction {0!s} declares electron-temperature-dependent "
                        "kinetics {1}, but that class exposes neither "
                        "get_rate_coefficient_two_temp(T, Te) nor "
                        "get_rate_coefficient_electron_temp(Te); evaluating it "
                        "through the standard k(T, P) interface would silently "
                        "collapse to one temperature. ({2})".format(
                            rxn, kin.__class__.__name__, self._identity()))
                # RULING 1: automatic thermodynamic reversal is prohibited for
                # Te-dependent kinetics. kf(Tgas, Te) / Keq(Tgas) combines two
                # incompatible thermal closures, and Keq(Te) would price
                # heavy-species thermochemistry at the electron temperature.
                # This holds regardless of the instantaneous numerical state:
                # Te == Tgas does not bypass it.
                if rxn.reversible:
                    raise NonEquilibriumReverseRateError(
                        "automatic thermodynamic reversal is undefined for "
                        "Te-dependent reaction {0!s}, kinetics {1}; mark the "
                        "reaction irreversible or provide explicit reverse "
                        "kinetics".format(rxn, kin.__class__.__name__))
                # An electron-temperature-dependent rate law is a plasma rate
                # law: it must carry an explicit incident electron among its
                # reactants, WHATEVER the reaction's net electron count. Without
                # one, an electron-collision rate silently loses its dependence
                # on the electron population and is evaluated at the wrong
                # reaction order. The metadata guard above and the family
                # resolver together ensure a metadata electron is placed or
                # refused before here; this closes the remaining gap -- a
                # plasma-rate reaction whose electron never appears explicitly at
                # all (e.g. electrons==0, no electron reactant) -- which the
                # evaluator and reversibility checks above do not catch. Placed
                # after them so those more specific diagnostics still take
                # precedence for the reactions they name.
                if n_electron_reactants == 0:
                    raise PlasmaStateError(
                        "reaction {0!s} has plasma-shaped kinetics {1} "
                        "(electron-temperature dependent), but no explicit "
                        "electron appears among its reactants; a plasma rate law "
                        "must never be evaluated without an explicit incident "
                        "electron, whatever the reaction's net electron count. "
                        "Represent the incident electron explicitly as a "
                        "reactant. ({2})".format(
                            rxn, kin.__class__.__name__, self._identity()))
            else:
                if rxn.reversible and (n_electron_reactants or n_electron_products):
                    raise NonEquilibriumReverseRateError(
                        "automatic thermodynamic reversal is undefined for "
                        "electron-containing reaction {0!s}, kinetics {1}: Keq(Tgas) "
                        "would price the electron's thermochemistry at the gas "
                        "temperature; mark the reaction irreversible or provide "
                        "explicit reverse kinetics".format(rxn, kin.__class__.__name__))

    def evaluate_two_temperature_rate_coefficient(self, kin):
        """
        Evaluate a ``uses_electron_temperature`` kinetics object at the
        reactor's (Tgas, Te), through the evaluator the class declares.
        The standard ``get_rate_coefficient(T, P)`` interface is deliberately
        never used for these kinetics: it collapses to a one-temperature
        evaluation (k(T, Te=T)).
        """
        if hasattr(kin, 'get_rate_coefficient_two_temp'):
            return kin.get_rate_coefficient_two_temp(self.T.value_si, self.Te.value_si)
        if hasattr(kin, 'get_rate_coefficient_electron_temp'):
            return kin.get_rate_coefficient_electron_temp(self.Te.value_si)
        raise PlasmaStateError(
            "kinetics {0} declares uses_electron_temperature but exposes no "
            "two-temperature evaluator. ({1})".format(
                kin.__class__.__name__, self._identity()))

    def generate_rate_coefficients(self, core_reactions, edge_reactions):
        """
        Populates the forward rate coefficients (kf), reverse rate
        coefficients (kb) and equilibrium constants (Keq) arrays.

        Thermal kinetics retain ordinary behaviour: kf at (Tgas, P) and, for
        reversible reactions, kr(Tgas) = kf(Tgas) / Keq(Tgas).

        Electron-temperature-dependent kinetics are evaluated through their
        declared two-temperature evaluator at (Tgas, Te).
        """
        if not self._plasma_validated:
            raise PlasmaStateError(
                "generate_rate_coefficients called before this reactor's "
                "composition and reaction set were validated; initialize the "
                "reactor through PlasmaReactor.initialize_model, not through "
                "the base class. ({0})".format(self._identity()))
        for rxn in itertools.chain(core_reactions, edge_reactions):
            j = self.reaction_index[rxn]
            kin = rxn.kinetics

            if getattr(kin, 'uses_electron_temperature', False):
                if rxn.reversible:
                    # Defensive: _validate_reactions() must already have
                    # rejected this configuration before we got here.
                    raise NonEquilibriumReverseRateError(
                        "automatic thermodynamic reversal is undefined for "
                        "Te-dependent reaction {0!s}, kinetics {1}; mark the "
                        "reaction irreversible or provide explicit reverse "
                        "kinetics".format(rxn, kin.__class__.__name__))
                self.kf[j] = self.evaluate_two_temperature_rate_coefficient(kin)
                self.kb[j] = 0.0
                self.Keq[j] = np.inf
            else:
                self.kf[j] = rxn.get_rate_coefficient(self.T.value_si, self.P.value_si)
                if rxn.reversible:
                    # Defensive: _validate_reactions() must already have
                    # rejected reversible electron-containing reactions.
                    if any(spc.is_electron() for spc in itertools.chain(rxn.reactants, rxn.products)):
                        raise NonEquilibriumReverseRateError(
                            "automatic thermodynamic reversal is undefined for "
                            "electron-containing reaction {0!s}, kinetics {1}; "
                            "mark the reaction irreversible or provide explicit "
                            "reverse kinetics".format(rxn, kin.__class__.__name__))
                    # The reactor-independent refusal. Every check above is more specific
                    # and fires first, so this cannot normally be reached from here; it is
                    # kept so that the reversal stays closed in this reactor even if those
                    # plasma-specific validations are ever loosened, and so that no reactor
                    # is the odd one out. See Reaction.get_reverse_from_equilibrium_refusal.
                    rxn.check_reverse_from_equilibrium_supported()
                    self.Keq[j] = rxn.get_equilibrium_constant(self.T.value_si)
                    self.kb[j] = self.kf[j] / self.Keq[j]
                else:
                    self.kb[j] = 0.0
                    self.Keq[j] = np.inf

    def _resolve_wall_state(self, list core_species):
        """
        Build the per-core-species charge bookkeeping the wall operator needs, and
        refuse -- loudly, here, before the first residual -- any configuration the
        operator cannot represent honestly.

        Three maps come out of this:

        ``species_charges``
            the net charge z_j of every core species, read from the species itself.

        ``wall_recycle_target``
            for every positive ion, the index of the NEUTRAL core species its heavy
            core returns as when it is neutralised at the wall. Matched on element
            composition, so Ar+ -> Ar. Resolved once, here, rather than guessed per
            step.

        ``source_cation_target``
            for every neutral heavy species, the index of its SINGLY-charged cation,
            which is what the external ionisation source produces. Background
            ionising radiation makes predominantly singly-charged ions; that is a
            statement about the mechanism, not a convenience.

        An ambiguous match is refused rather than resolved by a tie-break: two
        neutral argon states in the core (ground and metastable, both Ar, both
        uncharged) would make "the neutral Ar+ returns as" a coin flip that silently
        moves population between two species with different chemistry.
        """
        cdef Py_ssize_t n = self.num_core_species
        cdef Py_ssize_t i, j

        charges = np.zeros(n, dtype=np.int_)
        neutral_mask = np.zeros(n, dtype=np.uint8)
        recycle = np.full(n, -1, dtype=np.int_)
        source_cation = np.full(n, -1, dtype=np.int_)

        formulas = [None] * n
        for i in range(n):
            spc = core_species[i]
            try:
                charges[i] = int(spc.get_net_charge())
            except Exception as exc:
                raise PlasmaStateError(
                    "PlasmaReactor could not determine the net charge of core species "
                    "{0!r}, which the wall operator and the quasineutrality row both "
                    "require: {1}. ({2})".format(spc, exc, self._identity()))
            if i == self.electron_index:
                continue
            if charges[i] == 0:
                neutral_mask[i] = 1
            try:
                formulas[i] = tuple(sorted(spc.molecule[0].get_element_count().items()))
            except Exception:
                formulas[i] = None

        # ion -> neutral counterpart (same heavy composition, zero charge). Anions
        # recycle by the same rule as cations: an anion neutralised at a wall leaves
        # the same heavy core behind, and leaving them out would make their atoms
        # vanish silently.
        for i in range(n):
            if i == self.electron_index or charges[i] == 0 or formulas[i] is None:
                continue
            matches = [j for j in range(n)
                       if neutral_mask[j] and formulas[j] is not None and formulas[j] == formulas[i]]
            if len(matches) == 1:
                recycle[i] = matches[0]
            elif len(matches) > 1 and self.has_wall and self.wall_recycling > 0.0:
                raise PlasmaStateError(
                    "the wall recycling fate of cation {0!r} is ambiguous: {1} neutral "
                    "core species share its heavy composition ({2}), so 'the neutral it "
                    "returns as' has no single answer. Recycling would silently move "
                    "population into one of them. Set wall_recycling=0.0 to make the "
                    "wall fully pumping, or remove the duplicate neutral state from the "
                    "core. ({3})".format(
                        core_species[i], len(matches),
                        ', '.join(repr(core_species[j]) for j in matches), self._identity()))

        # neutral -> its singly-charged cation (what the external source makes)
        for i in range(n):
            if not neutral_mask[i] or formulas[i] is None:
                continue
            matches = [j for j in range(n)
                       if charges[j] == 1 and j != self.electron_index
                       and formulas[j] is not None and formulas[j] == formulas[i]]
            if len(matches) == 1:
                source_cation[i] = matches[0]
            elif len(matches) > 1 and self.has_wall and self.ionisation_source.value_si > 0.0:
                raise PlasmaStateError(
                    "the external ionisation source's product is ambiguous for neutral "
                    "{0!r}: {1} singly-charged cations in the core share its heavy "
                    "composition ({2}). Remove the duplicate cation, or set "
                    "ionisation_source to zero. ({3})".format(
                        core_species[i], len(matches),
                        ', '.join(repr(core_species[j]) for j in matches), self._identity()))

        if self.has_wall:
            if not neutral_mask.any():
                raise PlasmaStateError(
                    "a charged-particle wall was declared but the core carries no "
                    "neutral heavy species; the ion mobility used here scales as "
                    "1/n_neutral and is undefined without a neutral gas to collide "
                    "with. ({0})".format(self._identity()))
            if self.wall_recycling > 0.0:
                orphans = [core_species[i] for i in range(n)
                           if charges[i] != 0 and i != self.electron_index and recycle[i] < 0]
                if orphans:
                    raise PlasmaStateError(
                        "wall_recycling={0!r} asks the wall to return the heavy core of "
                        "every neutralised ion to the gas, but {1} ion(s) have no "
                        "neutral counterpart among the core species: {2}. Their heavy "
                        "atoms would vanish. Add the neutral(s) to the model, or declare "
                        "a fully pumping wall with wall_recycling=0.0 and account for "
                        "the loss in the wall inventory. ({3})".format(
                            self.wall_recycling, len(orphans),
                            ', '.join(repr(s) for s in orphans), self._identity()))
            if self.ionisation_source.value_si > 0.0 and not (source_cation >= 0).any():
                raise PlasmaStateError(
                    "ionisation_source={0!r} m^-3 s^-1 was declared but no neutral core "
                    "species has a singly-charged cation in the core, so there is no "
                    "ion-electron pair for the source to make. Add the cation, or set "
                    "ionisation_source to zero. ({1})".format(
                        self.ionisation_source.value_si, self._identity()))

        self.species_charges = charges
        self.neutral_heavy_mask = neutral_mask
        self.wall_recycle_target = recycle
        self.source_cation_target = source_cation
        self.wall_loss_rates = np.zeros(n, float)

    cpdef double compute_nu_wall(self, np.ndarray y, double V) except -1.0:
        """
        The common ambipolar loss frequency, s^-1, shared by EVERY charged species:

        .. math::

            \\mu_i = \\mu_0 \\frac{N_{ref}}{n_{neutral}}, \\quad
            D_a \\simeq \\mu_i \\frac{k_B T_e}{e}, \\quad
            \\nu_{wall} = \\frac{D_a}{\\Lambda^2}

        One frequency for all charged species, not one lifetime each: independent
        lifetimes let the net charge drift, whereas a common nu removes charge at a
        rate proportional to the net charge itself, which is exactly zero in a
        quasineutral gas. This is what makes the zero-net-wall-current rule hold by
        construction rather than by cancellation.

        ``D_a = mu_i * k_B*T_e/e`` is the ``T_e >> T_i`` limit of
        ``D_a = D_i (1 + T_e/T_i)``; it drops a factor ``(1 + T_i/T_e)``, which at
        the ratios this reactor runs at is an underestimate of order 1%. That
        approximation is carried as part of the transport interval, not corrected
        away silently.

        The ion mobility scales as ``1/n_neutral``, so ``nu_wall`` rises without
        bound as the neutral gas is consumed. That divergence is an artefact of
        assuming ion-NEUTRAL collisions dominate, which stops being true well before
        the neutrals run out. The domain is therefore enforced -- and the run stopped
        -- by :meth:`check_wall_support`, at the ACCEPTED-STEP boundary.

        **This method never raises, and that is deliberate.** An earlier version of
        it raised ``PlasmaStateError`` from here when the ionisation degree exceeded
        the ceiling. It does not work: the solver evaluates the residual at Newton
        TRIAL states, which are not physical states and routinely sit far outside the
        domain on their way to a converged step. A Python exception raised inside the
        Fortran callback is not propagated -- DASPK carries on with a corrupted
        iteration matrix ("All zero matrix", "ITERATION MATRIX IS SINGULAR") and the
        run finally dies with an unrelated-looking ``SystemError``, which is the exact
        opposite of the loud, readable failure the domain check exists to give. The
        check belongs where Python owns the stack and the state is one the solver has
        actually accepted.
        """
        cdef double y_neutral = 0.0, n_neutral, mu_i, d_a, lam
        cdef Py_ssize_t j
        if not self.has_wall:
            return 0.0
        for j in range(self.num_core_species):
            if self.neutral_heavy_mask[j]:
                y_neutral += y[j]
        # Numerical floor only; see wall_neutral_floor. Keeps a trial state with no
        # neutrals left finite so the solver can reject it on its own terms.
        if not (y_neutral > self.wall_neutral_floor):
            y_neutral = self.wall_neutral_floor
        n_neutral = y_neutral * constants.Na / V
        mu_i = self.ion_reduced_mobility.value_si * self.mobility_reference_density / n_neutral
        # k_B*T_e/e in volts. constants.R / constants.Na is the same Boltzmann
        # constant compute_volume's EOS is built from, so the two never disagree.
        d_a = mu_i * (constants.R / constants.Na) * self.Te.value_si / constants.e
        lam = self.diffusion_length.value_si
        return d_a / (lam * lam)

    cpdef check_wall_support(self, np.ndarray y):
        """
        Refuse an ACCEPTED state that lies outside the transport model's domain.

        The ion mobility used by :meth:`compute_nu_wall` assumes the ion's
        momentum-transfer collisions are with NEUTRALS, so it scales as
        ``1/n_neutral``. Above an ionisation degree of ``max_ionisation_degree``
        that assumption fails -- Coulomb collisions with the other charged particles
        take over -- and the 1/n_neutral scaling is not merely inaccurate but the
        wrong functional form, diverging where the physics does not.

        Called on the initial composition and after every accepted solver step, never
        from inside the residual: see the note in :meth:`compute_nu_wall` for why a
        Newton trial state must not be allowed to trigger this.
        """
        cdef double y_neutral = 0.0, alpha
        cdef Py_ssize_t j
        if not self.has_wall or self.neutral_heavy_mask is None:
            return
        for j in range(self.num_core_species):
            if self.neutral_heavy_mask[j]:
                y_neutral += y[j]
        if not np.isfinite(y_neutral) or y_neutral <= 0.0:
            raise PlasmaStateError(
                "the accepted state has a neutral heavy population of {0!r} mol. The "
                "ion mobility this wall model uses scales as 1/n_neutral and is "
                "undefined without a neutral gas to collide with, so the run stops "
                "here rather than extrapolating. ({1})".format(
                    y_neutral, self._identity()))
        alpha = y[self.electron_index] / y_neutral
        if alpha > self.max_ionisation_degree:
            raise PlasmaStateError(
                "ionisation degree n_e/n_neutral = {0!r} exceeds the declared validity "
                "ceiling {1!r} of the ion-neutral ambipolar transport model: above it "
                "Coulomb collisions dominate the ion mobility and the 1/n_neutral "
                "scaling used here is the wrong functional form, not merely an "
                "inaccurate one. Refusing to extrapolate. Raise max_ionisation_degree "
                "only together with a transport model that covers the regime. "
                "({2})".format(alpha, self.max_ionisation_degree, self._identity()))

    cpdef advance(self, double tout):
        """Advance, then refuse the accepted state if it left the wall model's domain.

        ``cpdef``, matching the base class, so the check runs whether the caller
        reaches this through Python or through a C-level call: a ``def`` override of
        a ``cpdef`` method is bypassed entirely by C-level callers, which would leave
        :meth:`ReactionSystem.simulate` running unguarded.
        """
        result = ReactionSystem.advance(self, tout)
        self.check_wall_support(self.y)
        return result

    cpdef step(self, double tout):
        """Step, then refuse the accepted state if it left the wall model's domain.

        This is the one :meth:`ReactionSystem.simulate` actually calls.
        """
        result = ReactionSystem.step(self, tout)
        self.check_wall_support(self.y)
        return result

    cpdef double compute_volume(self, np.ndarray y) except -1:
        """
        The single equation-of-state implementation:

        .. math:: V = R \\left( N_{heavy} T_{gas} + N_e T_e \\right) / P

        with ``N_e = y[electron_index]`` and ``N_heavy`` the sum of all other
        core species moles, both read from the supplied state vector. Used at
        initialization and by every residual and Jacobian evaluation.
        """
        cdef double n_total, n_e, volume
        if self.electron_index < 0:
            raise PlasmaStateError(
                "compute_volume called before the electron state was resolved. "
                "({0})".format(self._identity()))
        n_total = np.sum(y[:self.num_core_species])
        n_e = y[self.electron_index]
        volume = constants.R * ((n_total - n_e) * self.T.value_si
                                + n_e * self.Te.value_si) / self.P.value_si
        if not np.isfinite(volume) or volume <= 0.0:
            raise PlasmaStateError(
                "two-temperature EOS produced a non-physical volume {0!r} m^3 "
                "from N_total={1!r} mol, N_e={2!r} mol; refusing to continue "
                "rather than integrate a corrupted state. ({3})".format(
                    volume, n_total, n_e, self._identity()))
        return volume

    def set_initial_conditions(self):
        """
        Sets the initial conditions of the rate equations that represent the
        current reactor model.

        The volume is derived from the two-temperature equation of state via
        :meth:`compute_volume`, using the user-defined pressure, gas and
        electron temperatures, and the initial species moles.
        """
        ReactionSystem.set_initial_conditions(self)

        for spec, mole_frac in self.initial_mole_fractions.items():
            i = self.get_species_index(spec)
            self.y0[i] = mole_frac

        if not np.isfinite(self.y0[self.electron_index]) or self.y0[self.electron_index] <= 0.0:
            raise PlasmaStateError(
                "packed initial state has a non-positive electron amount "
                "{0!r} at index {1}. ({2})".format(
                    self.y0[self.electron_index], self.electron_index, self._identity()))

        # The quasineutrality row is an algebraic constraint, so the packed initial
        # state must SATISFY it. A state that does not is not merely inaccurate: the
        # solver reports it as "repeated error test failures ... a singularity in the
        # solution may be present" with no mention of charge, which is unreadable.
        # Refuse here, where the composition and the number are both in hand.
        #
        # Unlike _warn_if_not_charge_neutral, which warns because a deliberately
        # non-neutral initial condition can be legitimate for an INTEGRATED electron,
        # this raises: with the electron carried algebraically there is no such thing
        # as a non-neutral state, so a non-neutral initial condition is a request the
        # equations cannot represent.
        if self.quasineutral_electron and self.species_charges is not None:
            net = self._net_charge(self.y0)
            magnitude = float(np.sum(np.abs(self.species_charges * self.y0[:self.num_core_species])))
            if abs(net) > max(PLASMA_NET_CHARGE_ATOL, PLASMA_NET_CHARGE_RTOL * magnitude):
                raise PlasmaStateError(
                    "quasineutral_electron=True carries the electron on an algebraic "
                    "charge-conservation row, but the initial composition has a net "
                    "charge of {0!r} mol (total charge magnitude {1!r}), which that row "
                    "does not admit. The solver would fail with an unrelated-looking "
                    "convergence error. Make the initial composition neutral -- "
                    "chargeBalanceSpecies='<ion>' in the plasmaReactor(...) block does "
                    "it -- or integrate the electron instead. ({2})".format(
                        net, magnitude, self._identity()))

        # Two-temperature ideal gas law, same implementation as residual/Jacobian
        self.V = self.compute_volume(self.y0)
        for j in range(self.num_core_species):
            self.core_species_concentrations[j] = self.y0[j] / self.V

        # The wall's numerical neutral floor, and the domain check on the initial
        # state. The floor is a millionth of the initial neutral population: it is
        # reached only at ionisation degrees a million times the validity ceiling,
        # so it can never bind on a state the domain check would accept.
        if self.has_wall and self.neutral_heavy_mask is not None:
            y_neutral0 = 0.0
            for j in range(self.num_core_species):
                if self.neutral_heavy_mask[j]:
                    y_neutral0 += self.y0[j]
            if y_neutral0 <= 0.0:
                raise PlasmaStateError(
                    "a charged-particle wall was declared but the initial composition "
                    "contains no neutral heavy species ({0!r} mol); the ion mobility "
                    "scales as 1/n_neutral and has nothing to scale against. "
                    "({1})".format(y_neutral0, self._identity()))
            self.wall_neutral_floor = 1.0e-6 * y_neutral0
            self.check_wall_support(self.y0)

    @cython.boundscheck(False)
    def residual(self, double t, np.ndarray[np.float64_t, ndim=1] y, np.ndarray[np.float64_t, ndim=1] dydt,
                 np.ndarray[np.float64_t, ndim=1] senpar = np.zeros(1, float)):
        """
        Return the residual function for the governing DAE system for the
        plasma reaction system.
        """
        cdef np.ndarray[np.int_t, ndim=2] ir, ip, inet
        cdef np.ndarray[np.float64_t, ndim=1] res, kf, kr, knet, delta
        cdef Py_ssize_t num_core_species, num_core_reactions, num_edge_species, num_edge_reactions
        cdef Py_ssize_t i, j, z, first, second, third
        cdef double k, V, reaction_rate, rev_reaction_rate, f_reaction_rate
        cdef np.ndarray[np.float64_t, ndim=1] core_species_concentrations, core_species_rates, core_reaction_rates
        cdef np.ndarray[np.float64_t, ndim=1] edge_species_rates, edge_reaction_rates, network_leak_rates
        cdef np.ndarray[np.float64_t, ndim=1] core_species_consumption_rates, core_species_production_rates
        cdef np.ndarray[np.float64_t, ndim=1] C, y_core_species

        ir = self.reactant_indices
        ip = self.product_indices

        num_core_species = len(self.core_species_rates)
        num_core_reactions = len(self.core_reaction_rates)
        num_edge_species = len(self.edge_species_rates)
        num_edge_reactions = len(self.edge_reaction_rates)
        kf = self.kf
        kr = self.kb

        y_core_species = y[:num_core_species]

        inet = self.network_indices
        knet = self.network_leak_coefficients

        res = np.zeros(num_core_species, float)

        core_species_concentrations = np.zeros_like(self.core_species_concentrations)
        core_species_rates = np.zeros_like(self.core_species_rates)
        core_reaction_rates = np.zeros_like(self.core_reaction_rates)
        core_species_consumption_rates = np.zeros_like(self.core_species_consumption_rates)
        core_species_production_rates = np.zeros_like(self.core_species_production_rates)
        edge_species_rates = np.zeros_like(self.edge_species_rates)
        edge_reaction_rates = np.zeros_like(self.edge_reaction_rates)
        network_leak_rates = np.zeros_like(self.network_leak_rates)

        C = np.zeros_like(self.core_species_concentrations)

        # Two-temperature ideal gas law, from the current state, through the
        # single shared EOS implementation.
        V = self.compute_volume(y)
        self.V = V

        for j in range(num_core_species):
            C[j] = y[j] / V
            core_species_concentrations[j] = C[j]

        for j in range(ir.shape[0]):
            k = kf[j]
            if ir[j, 0] >= num_core_species or ir[j, 1] >= num_core_species or ir[j, 2] >= num_core_species:
                f_reaction_rate = 0.0
            elif ir[j, 1] == -1:  # only one reactant
                f_reaction_rate = k * C[ir[j, 0]]
            elif ir[j, 2] == -1:  # only two reactants
                f_reaction_rate = k * C[ir[j, 0]] * C[ir[j, 1]]
            else:  # three reactants
                f_reaction_rate = k * C[ir[j, 0]] * C[ir[j, 1]] * C[ir[j, 2]]
            k = kr[j]
            if ip[j, 0] >= num_core_species or ip[j, 1] >= num_core_species or ip[j, 2] >= num_core_species:
                rev_reaction_rate = 0.0
            elif ip[j, 1] == -1:  # only one product
                rev_reaction_rate = k * C[ip[j, 0]]
            elif ip[j, 2] == -1:  # only two products
                rev_reaction_rate = k * C[ip[j, 0]] * C[ip[j, 1]]
            else:  # three products
                rev_reaction_rate = k * C[ip[j, 0]] * C[ip[j, 1]] * C[ip[j, 2]]

            reaction_rate = f_reaction_rate - rev_reaction_rate

            # Set the reaction and species rates
            if j < num_core_reactions:
                # The reaction is a core reaction
                core_reaction_rates[j] = reaction_rate

                # Add/substract the total reaction rate from each species rate
                # Since it's a core reaction we know that all of its reactants
                # and products are core species
                first = ir[j, 0]
                core_species_rates[first] -= reaction_rate
                core_species_consumption_rates[first] += f_reaction_rate
                core_species_production_rates[first] += rev_reaction_rate
                second = ir[j, 1]
                if second != -1:
                    core_species_rates[second] -= reaction_rate
                    core_species_consumption_rates[second] += f_reaction_rate
                    core_species_production_rates[second] += rev_reaction_rate
                    third = ir[j, 2]
                    if third != -1:
                        core_species_rates[third] -= reaction_rate
                        core_species_consumption_rates[third] += f_reaction_rate
                        core_species_production_rates[third] += rev_reaction_rate
                first = ip[j, 0]
                core_species_rates[first] += reaction_rate
                core_species_production_rates[first] += f_reaction_rate
                core_species_consumption_rates[first] += rev_reaction_rate
                second = ip[j, 1]
                if second != -1:
                    core_species_rates[second] += reaction_rate
                    core_species_production_rates[second] += f_reaction_rate
                    core_species_consumption_rates[second] += rev_reaction_rate
                    third = ip[j, 2]
                    if third != -1:
                        core_species_rates[third] += reaction_rate
                        core_species_production_rates[third] += f_reaction_rate
                        core_species_consumption_rates[third] += rev_reaction_rate

            else:
                # The reaction is an edge reaction
                edge_reaction_rates[j - num_core_reactions] = reaction_rate

                # Add/substract the total reaction rate from each species rate
                # Since it's an edge reaction its reactants and products could
                # be either core or edge species
                # We're only interested in the edge species
                first = ir[j, 0]
                if first >= num_core_species: edge_species_rates[first - num_core_species] -= reaction_rate
                second = ir[j, 1]
                if second != -1:
                    if second >= num_core_species: edge_species_rates[second - num_core_species] -= reaction_rate
                    third = ir[j, 2]
                    if third != -1:
                        if third >= num_core_species: edge_species_rates[third - num_core_species] -= reaction_rate
                first = ip[j, 0]
                if first >= num_core_species: edge_species_rates[first - num_core_species] += reaction_rate
                second = ip[j, 1]
                if second != -1:
                    if second >= num_core_species: edge_species_rates[second - num_core_species] += reaction_rate
                    third = ip[j, 2]
                    if third != -1:
                        if third >= num_core_species: edge_species_rates[third - num_core_species] += reaction_rate

        for j in range(inet.shape[0]):
            if inet[j, 0] != -1:  # all source species are in the core
                k = knet[j]
                if inet[j, 1] == -1:  # only one reactant
                    reaction_rate = k * C[inet[j, 0]]
                elif inet[j, 2] == -1:  # only two reactants
                    reaction_rate = k * C[inet[j, 0]] * C[inet[j, 1]]
                else:  # three reactants
                    reaction_rate = k * C[inet[j, 0]] * C[inet[j, 1]] * C[inet[j, 2]]
                network_leak_rates[j] = reaction_rate
            else:
                network_leak_rates[j] = 0.0

        self.core_species_concentrations = core_species_concentrations
        self.core_species_rates = core_species_rates
        self.core_species_production_rates = core_species_production_rates
        self.core_species_consumption_rates = core_species_consumption_rates
        self.core_reaction_rates = core_reaction_rates
        self.edge_species_rates = edge_species_rates
        self.edge_reaction_rates = edge_reaction_rates
        self.network_leak_rates = network_leak_rates

        res = core_species_rates * V

        # --- charged-particle wall boundary and external ionisation source ----
        #
        # Deliberately applied to `res` and NOT folded into core_species_rates:
        # that array is the gas-phase CHEMISTRY diagnostic that model enlargement
        # and the rate-ratio termination criteria compare edge fluxes against, and
        # a transport term is not a reaction flux. Keeping them apart is also what
        # constraint "separate wall transport from gas-phase mechanism chemistry"
        # asks for, in the one place it can actually be enforced.
        #
        # res is in mol/s and equals (dC/dt)*V with C = y/V, so a first-order loss
        # nu*C contributes exactly -nu*y_i: the volume cancels ALGEBRAICALLY here.
        # It is written as -nu*y_i rather than -nu*(y_i/V)*V so that the
        # cancellation is exact in floating point too, not merely to one ulp.
        if self.has_wall:
            self._apply_wall_terms(y, V, res)

        delta = res - dydt

        # --- quasineutrality replaces the electron's differential row ---------
        # delta[e] carries NO dydt term: the electron amount is not integrated, it
        # is whatever makes the packed state charge neutral. This is an index-1
        # algebraic constraint (d delta[e] / d y_e = z_e = -1, never singular), and
        # it removes from the Jacobian the stiff direction in which n_e is a small
        # difference of large ionisation and recombination fluxes.
        if self.quasineutral_electron:
            delta[self.electron_index] = self.charge_row_scale * self._net_charge(y)

        # Return DELTA, IRES.  IRES is set to 1 in order to tell DASPK to evaluate the sensitivity residuals
        return delta, 1

    def _set_charge_row_scale(self):
        """
        Choose the multiplier on the algebraic quasineutrality row so that it is
        commensurate with the differential rows around it.

        The row's own entries are the charges, exactly +/-1, while the differential
        rows carry entries of order of the rates. DASPK factors the iteration matrix
        with partial pivoting and NO equilibration, so leaving a row of order-1
        coefficients among rows of order 1e2-1e5 is a real and measurable cost, not a
        cosmetic one. The multiplier is the largest magnitude among the differential
        rows of the Jacobian at the initial state, which is what an equilibration step
        would use and which requires no chosen constant.

        Scaling one row of a DAE residual by a positive constant is exact: it changes
        the matrix that is factored and nothing about the solution. The measured
        dependence is flat over three decades of the multiplier, so the choice does
        not need to be sharp -- only commensurate.
        """
        cdef Py_ssize_t i
        pd = np.array(self.jacobian(0.0, self.y0.copy(),
                                    np.zeros(self.neq, float), 0.0), float)
        differential = np.delete(np.abs(pd), self.electron_index, axis=0)
        scale = float(np.max(differential)) if differential.size else 0.0
        if np.isfinite(scale) and scale > 0.0:
            self.charge_row_scale = scale
        else:
            self.charge_row_scale = 1.0
        logging.debug(
            'PlasmaReactor: quasineutrality row scaled by %r to match the '
            'differential rows (exact row scaling; changes no solution). (%s)',
            self.charge_row_scale, self._identity())

    cdef double _net_charge(self, np.ndarray[np.float64_t, ndim=1] y):
        """Net charge of the packed state, in mol of elementary charge."""
        cdef double net = 0.0
        cdef Py_ssize_t j
        for j in range(self.num_core_species):
            if self.species_charges[j] != 0:
                net += self.species_charges[j] * y[j]
        return net

    cdef _apply_wall_terms(self, np.ndarray[np.float64_t, ndim=1] y, double V,
                           np.ndarray[np.float64_t, ndim=1] res):
        """
        Add the wall loss, its recycled return, and the external pair source to
        ``res`` (mol/s), and record the per-species wall flux in
        ``self.wall_loss_rates`` so a caller can integrate the flux INDEPENDENTLY
        of the gas-phase deficit and compare the two.
        """
        cdef double nu, loss, y_neutral = 0.0, source_total, rate
        cdef Py_ssize_t j, target
        cdef np.ndarray[np.float64_t, ndim=1] wall

        nu = self.compute_nu_wall(y, V)
        self.nu_wall = nu
        wall = self.wall_loss_rates
        for j in range(self.num_core_species):
            wall[j] = 0.0
            if self.neutral_heavy_mask[j]:
                y_neutral += y[j]

        # One common nu for every charged species -- electron and every ion alike.
        # The net charge therefore decays at rate nu*(net charge), which is zero in
        # a neutral gas, so no charge is created or destroyed at the wall.
        for j in range(self.num_core_species):
            if self.species_charges[j] == 0:
                continue
            loss = nu * y[j]
            res[j] -= loss
            wall[j] -= loss
            if j == self.electron_index:
                continue
            target = self.wall_recycle_target[j]
            if target >= 0 and self.wall_recycling > 0.0:
                # gamma of the neutralised heavy core returns to the gas; the
                # remaining (1-gamma) stays on the wall and leaves the gas phase.
                res[target] += self.wall_recycling * loss
                wall[target] += self.wall_recycling * loss

        # External volumetric source of ion-electron pairs (the declared physical
        # mechanism that replaces the numerical seed). Its total production is
        # S_ext * V pairs per second, apportioned across the ionisable neutrals in
        # proportion to their number density -- which is what an ionising radiation
        # field does, at the stated assumption of equal per-particle ionisation
        # probability. For a single-neutral gas the apportionment is the identity.
        source_total = self.ionisation_source.value_si * V / constants.Na
        if source_total > 0.0 and y_neutral > 0.0:
            for j in range(self.num_core_species):
                target = self.source_cation_target[j]
                if not self.neutral_heavy_mask[j] or target < 0:
                    continue
                rate = source_total * y[j] / y_neutral
                res[j] -= rate
                res[target] += rate
                res[self.electron_index] += rate

    @cython.boundscheck(False)
    def jacobian(self, double t, np.ndarray[np.float64_t, ndim=1] y, np.ndarray[np.float64_t, ndim=1] dydt,
                 double cj, np.ndarray[np.float64_t, ndim=1] senpar = np.zeros(1, float)):
        """
        Return the analytical Jacobian for the reaction system.

        Uses the same reaction set (kf/kb computed once at initialization),
        the same electron resolution (electron_index), and the same EOS
        implementation (:meth:`compute_volume`) as :meth:`residual`. The
        volume derivative dV/dy_i is R*Tgas/P for heavy species and R*Te/P
        for the electron column.
        """
        cdef np.ndarray[np.int_t, ndim=2] ir, ip, rrow, prow
        cdef np.ndarray[np.float64_t, ndim=1] kf, kr, C, dVdy
        cdef np.ndarray[np.float64_t, ndim=2] pd
        cdef Py_ssize_t num_core_reactions, num_core_species, i, j, m, m2, n, direction
        cdef int z, s
        cdef double k, V, partial, prod_c, corr
        cdef list rlist, plist

        ir = self.reactant_indices
        ip = self.product_indices

        kf = self.kf
        kr = self.kb
        num_core_reactions = len(self.core_reaction_rates)
        num_core_species = len(self.core_species_concentrations)

        pd = -cj * np.identity(num_core_species, float)

        # Same EOS implementation as residual(), same state.
        V = self.compute_volume(y)

        # Analytic dV/dy_i of the same EOS, from the same live T, P, Te:
        # R*Tgas/P for heavy species, R*Te/P for the electron column.
        dVdy = np.full(num_core_species, constants.R * self.T.value_si / self.P.value_si, float)
        dVdy[self.electron_index] = constants.R * self.Te.value_si / self.P.value_si

        C = np.zeros_like(self.core_species_concentrations)
        for j in range(num_core_species):
            C[j] = y[j] / V

        for j in range(num_core_reactions):
            for direction in range(2):
                if direction == 0:
                    k = kf[j]
                    rrow = ir
                    prow = ip
                else:
                    k = kr[j]
                    rrow = ip
                    prow = ir
                if k == 0.0:
                    continue

                rlist = [rrow[j, m] for m in range(3) if rrow[j, m] != -1]
                plist = [prow[j, m] for m in range(3) if prow[j, m] != -1]
                n = len(rlist)

                # The species rate contribution of this direction is
                # +/- V * k * prod(C_r). Differentiate first through the
                # explicit y dependence (one term per reactant occurrence) ...
                for m in range(n):
                    z = rlist[m]
                    partial = k
                    for m2 in range(n):
                        if m2 != m:
                            partial *= C[rlist[m2]]
                    for s in rlist:
                        pd[s, z] -= partial
                    for s in plist:
                        pd[s, z] += partial

                # ... then through the volume: d(V^(1-n))/dy_i contributes
                # -(n-1) * k * prod(C_r) * dV/dy_i for every column i.
                if n > 1:
                    prod_c = k
                    for m in range(n):
                        prod_c *= C[rlist[m]]
                    corr = (n - 1) * prod_c
                    for i in range(num_core_species):
                        for s in rlist:
                            pd[s, i] += corr * dVdy[i]
                        for s in plist:
                            pd[s, i] -= corr * dVdy[i]

        # --- the wall operator's Jacobian, term for term with the residual ----
        #
        # nu_wall is NOT a constant: the ion mobility scales as 1/n_neutral, so
        #     nu = K * V / (Na * Y),   Y = sum of the neutral heavy moles
        # and therefore
        #     d nu / d y_k = nu * ( dVdy_k / V  -  [k is neutral heavy] / Y ).
        # This is where the brief's warning that a wall term "also perturbs V" is
        # real. It is NOT real through the C = y/V conversion, which cancels V
        # exactly; it is real through the MOBILITY, and only for a
        # composition-dependent one.
        if self.has_wall:
            self._apply_wall_jacobian(y, V, dVdy, pd)

        self.jacobian_matrix = pd + cj * np.identity(num_core_species, float)

        # The quasineutrality row is purely algebraic: d delta[e] / d dydt[e] = 0,
        # so it carries no -cj term and must be written AFTER the accumulation
        # above (which started every row at -cj*I).
        if self.quasineutral_electron:
            for i in range(num_core_species):
                pd[self.electron_index, i] = self.charge_row_scale * self.species_charges[i]
                self.jacobian_matrix[self.electron_index, i] = \
                    self.charge_row_scale * self.species_charges[i]

        return pd

    cdef _apply_wall_jacobian(self, np.ndarray[np.float64_t, ndim=1] y, double V,
                              np.ndarray[np.float64_t, ndim=1] dVdy,
                              np.ndarray[np.float64_t, ndim=2] pd):
        """Partial derivatives of exactly the terms :meth:`_apply_wall_terms` adds."""
        cdef double nu, y_neutral = 0.0, loss, dnu_rel, dloss, source_total, term, dterm
        cdef Py_ssize_t i, k, target
        cdef double gamma = self.wall_recycling
        cdef bint neutral_floored = 0

        nu = self.compute_nu_wall(y, V)
        for k in range(self.num_core_species):
            if self.neutral_heavy_mask[k]:
                y_neutral += y[k]
        # compute_nu_wall clamps the neutral moles at wall_neutral_floor, so on that
        # branch nu no longer depends on the neutral amount through 1/y_neutral --
        # only through V, which keeps its EOS dependence. Differentiating the
        # unclamped 1/y_neutral here would report a dependence the residual does not
        # have, and divide by zero once the neutrals are gone. The source term below
        # is apportioned by the UNFLOORED total, exactly as _apply_wall_terms does,
        # so the two uses cannot share one value.
        if not (y_neutral > self.wall_neutral_floor):
            neutral_floored = 1

        for i in range(self.num_core_species):
            if self.species_charges[i] == 0:
                continue
            loss = nu * y[i]
            target = self.wall_recycle_target[i] if i != self.electron_index else -1
            for k in range(self.num_core_species):
                # d(nu)/dy_k, relative to nu
                dnu_rel = dVdy[k] / V
                if self.neutral_heavy_mask[k] and not neutral_floored:
                    dnu_rel -= 1.0 / y_neutral
                dloss = loss * dnu_rel
                if k == i:
                    dloss += nu
                pd[i, k] -= dloss
                if target >= 0 and gamma > 0.0:
                    pd[target, k] += gamma * dloss

        source_total = self.ionisation_source.value_si * V / constants.Na
        if source_total > 0.0 and y_neutral > 0.0:
            for i in range(self.num_core_species):
                target = self.source_cation_target[i]
                if not self.neutral_heavy_mask[i] or target < 0:
                    continue
                term = source_total * y[i] / y_neutral
                for k in range(self.num_core_species):
                    dterm = term * dVdy[k] / V
                    if self.neutral_heavy_mask[k]:
                        dterm -= term / y_neutral
                    if k == i:
                        dterm += source_total / y_neutral
                    pd[i, k] -= dterm
                    pd[target, k] += dterm
                    pd[self.electron_index, k] += dterm
