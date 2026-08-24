# cython: embedsignature=True, cdivision=True

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

import numpy as np
cimport numpy as np
import os
from libc.math cimport exp, sqrt, log10, pow
from scipy.optimize import curve_fit, fsolve

cimport rmgpy.constants as constants
from rmgpy import settings
import rmgpy.quantity as quantity
from rmgpy.exceptions import KineticsError
from rmgpy.kinetics.uncertainties import rank_accuracy_map
from rmgpy.molecule.molecule import Bond
from rmgpy.kinetics.model import KineticsModel, PDepKineticsModel
import logging

# Prior to numpy 1.14, `numpy.linalg.lstsq` does not accept None as a value
RCOND = -1 if int(np.__version__.split('.')[1]) < 14 else None
################################################################################

cdef class Arrhenius(KineticsModel):
    """
    A kinetics model based on the (modified) Arrhenius equation. The attributes
    are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `A`             The preexponential factor
    `T0`            The reference temperature
    `n`             The temperature exponent
    `Ea`            The activation energy
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `solute`        Transition state solute data
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================

    """

    def __init__(self, A=None, n=0.0, Ea=None, T0=(1.0, "K"), Tmin=None, Tmax=None, Pmin=None, Pmax=None,
                 uncertainty=None, solute=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, uncertainty=uncertainty,
                               solute=solute, comment=comment)
        self.A = A
        self.n = n
        self.Ea = Ea
        self.T0 = T0

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Arrhenius object.
        """
        string = 'Arrhenius(A={0!r}, n={1!r}, Ea={2!r}, T0={3!r}'.format(self.A, self.n, self.Ea, self.T0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.uncertainty: string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.solute: string += ', solute={0!r}'.format(self.solute)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an Arrhenius object.
        """
        return (Arrhenius, (self.A, self.n, self.Ea, self.T0, self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                            self.uncertainty, self.solute, self.comment))

    property A:
        """The preexponential factor."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.RateCoefficient(value)

    property n:
        """The temperature exponent."""
        def __get__(self):
            return self._n
        def __set__(self, value):
            self._n = quantity.Dimensionless(value)

    property Ea:
        """The activation energy."""
        def __get__(self):
            return self._Ea
        def __set__(self, value):
            self._Ea = quantity.Energy(value)

    property T0:
        """The reference temperature."""
        def __get__(self):
            return self._T0
        def __set__(self, value):
            self._T0 = quantity.Temperature(value)

    cpdef double get_rate_coefficient(self, double T, double P=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3,
        mol, and s at temperature `T` in K.
        """
        cdef double A, n, Ea, T0
        A = self._A.value_si
        n = self._n.value_si
        Ea = self._Ea.value_si
        T0 = self._T0.value_si
        return A * (T / T0) ** n * exp(-Ea / (constants.R * T))

    cpdef change_t0(self, double T0):
        """
        Changes the reference temperature used in the exponent to `T0` in K,
        and adjusts the preexponential factor accordingly.
        """
        self._A.value_si /= (self._T0.value_si / T0) ** self._n.value_si
        self._T0.value_si = T0

    cpdef fit_to_data(self, np.ndarray Tlist, np.ndarray klist, str kunits, double T0=1, np.ndarray weights=None,
                      bint three_params=True):
        """
        Fit the Arrhenius parameters to a set of rate coefficient data `klist`
        in units of `kunits` corresponding to a set of temperatures `Tlist` in
        K. A linear least-squares fit is used, which guarantees that the
        resulting parameters provide the best possible approximation to the
        data.
        """
        import scipy.stats
        if not all(np.isfinite(klist)):
            raise  ValueError("Rates must all be finite, not inf or NaN")
        if any(klist<0):
            if not all(klist<0):
                raise ValueError("Rates must all be positive or all be negative.")
            rate_sign_multiplier = -1
            klist = -1 * klist
        else:
            rate_sign_multiplier = 1

        assert len(Tlist) == len(klist), "length of temperatures and rates must be the same"
        if len(Tlist) < 3 + three_params:
            raise KineticsError('Not enough degrees of freedom to fit this Arrhenius expression')
        if three_params:
            A = np.zeros((len(Tlist), 3), float)
            A[:, 0] = np.ones_like(Tlist)
            A[:, 1] = np.log(Tlist / T0)
            A[:, 2] = -1.0 / constants.R / Tlist
        else:
            A = np.zeros((len(Tlist), 2), float)
            A[:, 0] = np.ones_like(Tlist)
            A[:, 1] = -1.0 / constants.R / Tlist
        b = np.log(klist)
        if weights is not None:
            for n in range(b.size):
                A[n, :] *= weights[n]
                b[n] *= weights[n]
        x, residues, rank, s = np.linalg.lstsq(A, b, rcond=RCOND)

        # Determine covarianace matrix to obtain parameter uncertainties
        count = klist.size
        cov = residues[0] / (count - 3) * np.linalg.inv(np.dot(A.T, A))
        t = scipy.stats.t.ppf(0.975, count - 3)

        if not three_params:
            x = np.array([x[0], 0, x[1]])
            cov = np.array([[cov[0, 0], 0, cov[0, 1]], [0, 0, 0], [cov[1, 0], 0, cov[1, 1]]])

        self.A = (rate_sign_multiplier * exp(x[0]), kunits)
        self.n = x[1]
        self.Ea = (x[2] * 0.001, "kJ/mol")
        self.T0 = (T0, "K")
        self.Tmin = (np.min(Tlist), "K")
        self.Tmax = (np.max(Tlist), "K")
        self.solute = None
        self.comment = 'Fitted to {0:d} data points; dA = *|/ {1:g}, dn = +|- {2:g}, dEa = +|- {3:g} kJ/mol'.format(
            len(Tlist),
            exp(sqrt(cov[0, 0])),
            sqrt(cov[1, 1]),
            sqrt(cov[2, 2]) * 0.001,
        )

        return self

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Must match temperature
        and pressure range of kinetics model, as well as parameters: A, n, Ea, T0. (Shouldn't have pressure
        range if it's Arrhenius.) Otherwise returns ``False``.
        """
        if not isinstance(other_kinetics, Arrhenius):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if (not self.A.equals(other_kinetics.A) or not self.n.equals(other_kinetics.n)
                or not self.Ea.equals(other_kinetics.Ea) or not self.T0.equals(other_kinetics.T0)):
            return False

        return True

    cpdef change_rate(self, double factor):
        """
        Changes A factor in Arrhenius expression by multiplying it by a ``factor``.
        """
        self._A.value_si *= factor

    def to_cantera_kinetics(self, arrhenius_class=False):
        """
        Converts the RMG Arrhenius object to a cantera ArrheniusRate or
        the auxiliary cantera Arrhenius class (used by falloff reactions). 
        Inputs for both are (A,b,E)  where A is in units of m^3/kmol/s, b is dimensionless, and E is in J/kmol

        arrhenius_class: If ``True``, uses cantera.Arrhenius (for falloff reactions). If ``False``, uses 
        Cantera.ArrheniusRate
        """

        import cantera as ct

        rate_units_dimensionality = {'1/s': 0,
                                     's^-1': 0,
                                     'm^3/(mol*s)': 1,
                                     'm^6/(mol^2*s)': 2,
                                     'cm^3/(mol*s)': 1,
                                     'cm^6/(mol^2*s)': 2,
                                     'm^3/(molecule*s)': 1,
                                     'm^6/(molecule^2*s)': 2,
                                     'cm^3/(molecule*s)': 1,
                                     'cm^6/(molecule^2*s)': 2,
                                     }

        if self._T0.value_si != 1:
            A = self._A.value_si / (self._T0.value_si) ** self._n.value_si
        else:
            A = self._A.value_si

        try:
            A *= 1000 ** rate_units_dimensionality[self._A.units]
        except KeyError:
            raise Exception('Arrhenius A-factor units {0} not found among accepted units for converting to '
                            'Cantera Arrhenius object.'.format(self._A.units))

        b = self._n.value_si
        E = self._Ea.value_si * 1000  # convert from J/mol to J/kmol
        if arrhenius_class:
            return ct.Arrhenius(A, b, E)
        else:
            rate = ct.ArrheniusRate(A, b, E)
            if A < 0:
                rate.allow_negative_pre_exponential_factor = True
            return rate

    def set_cantera_kinetics(self, ct_reaction, species_list):
        """
        Passes in a cantera Reaction() object and sets its
        rate to a Cantera ArrheniusRate object.
        """
        import cantera as ct
        assert isinstance(ct_reaction.rate, ct.ArrheniusRate), "Must have a Cantera ArrheniusRate attribute"

        # Set the rate parameter to a cantera Arrhenius object
        ct_reaction.rate = self.to_cantera_kinetics()

    cpdef ArrheniusEP to_arrhenius_ep(self, double alpha=0.0, double dHrxn=0.0):
        """
        Converts an Arrhenius object to ArrheniusEP

        If setting alpha, you need to also input dHrxn, which must be given
        in J/mol (and vise versa).
        """

        if bool(alpha) ^ bool(dHrxn):
            raise Exception('If you set alpha or dHrxn in to_arrhenius_ep, '
                            'you need to set the other value to non-zero.')
        self.change_t0(1)
        aep = ArrheniusEP(A=self.A,
                          n=self.n,
                          alpha=alpha,
                          E0=(self.Ea.value_si - alpha * dHrxn, 'J/mol'),
                          Tmin=self.Tmin,
                          Tmax=self.Tmax,
                          Pmin=self.Pmin,
                          Pmax=self.Pmax,
                          uncertainty=self.uncertainty,
                          solute=self.solute,
                          comment=self.comment)
        return aep
################################################################################


cdef class TwoTemperaturePlasma(KineticsModel):
    """
    Two-temperature plasma kinetics, where the rate depends on both gas
    temperature T and electron temperature Te.

    The underlying functional form (Kossyi-type) is

        k(T, Te) = A * Te^n
                   * exp(-Ea_g / (R * T))
                   * exp(Ea_e * (Te - T) / (R * T * Te))

    where:
        A      : pre-exponential factor
        n      : electron temperature exponent
        Ea_g   : gas activation energy
        Ea_e   : electron activation energy

    In RMG's standard `get_rate_coefficient(T)` interface, we evaluate
    k(T, Te=T) as a reasonable fallback. A dedicated
    `get_rate_coefficient_two_temp(T, Te)` method is provided to use
    distinct gas and electron temperatures explicitly.
    """

    def __init__(self,
                 A=None,
                 n=0.0,
                 Ea_g=(0.0, "J/mol"),
                 Ea_e=(0.0, "J/mol"),
                 T0=(1.0, "K"),
                 Tmin=None, Tmax=None, Pmin=None, Pmax=None,
                 uncertainty=None, solute=None, comment=''):
        KineticsModel.__init__(
            self,
            Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax,
            uncertainty=uncertainty, solute=solute, comment=comment,
        )
        self.A = A
        self.n = n
        self.Ea_g = Ea_g
        self.Ea_e = Ea_e
        self.T0 = T0
        self.uses_electron_temperature = True

    def __repr__(self):
        string = 'TwoTemperaturePlasma(A={0!r}, n={1!r}, Ea_g={2!r}, Ea_e={3!r}'.format(self.A, self.n, self.Ea_g, self.Ea_e)
        if self.T0.value_si != 1:
            string += ', T0={0!r}'.format(self.T0)
        if self.Tmin is not None:
            string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None:
            string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None:
            string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None:
            string += ', Pmax={0!r}'.format(self.Pmax)
        if self.uncertainty:
            string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.solute:
            string += ', solute={0!r}'.format(self.solute)
        if self.comment != '':
            string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        Helper for pickling.
        """
        return (TwoTemperaturePlasma,
                (self.A, self.n, self.Ea_g, self.Ea_e, self.T0,
                 self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                 self.uncertainty, self.solute, self.comment))

    property A:
        """Pre-exponential factor (same units conventions as Arrhenius A)."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.RateCoefficient(value)

    property n:
        """Electron temperature exponent."""
        def __get__(self):
            return self._n
        def __set__(self, value):
            self._n = quantity.Dimensionless(value)

    property Ea_g:
        """Gas activation energy."""
        def __get__(self):
            return self._Ea_g
        def __set__(self, value):
            self._Ea_g = quantity.Energy(value)

    property Ea_e:
        """Electron activation energy."""
        def __get__(self):
            return self._Ea_e
        def __set__(self, value):
            self._Ea_e = quantity.Energy(value)

    property T0:
        """Reference temperature."""
        def __get__(self):
            return self._T0
        def __set__(self, value):
            self._T0 = quantity.Temperature(value)

    cpdef double get_rate_coefficient_two_temp(self, double T, double Te) except -1:
        """
        Full two-temperature rate coefficient k(T, Te) in SI molar units.
        T  : gas temperature [K]
        Te : electron temperature [K]
        """
        cdef double A_si, bval, Eag, Eae, T0

        if T <= 0.0 or Te <= 0.0:
            raise ValueError("TwoTemperaturePlasma: T and Te must be > 0 K.")

        A_si = self._A.value_si           # e.g. m^3/(mol*s)
        bval = self._n.value_si
        T0   = self._T0.value_si          # K
        Eag  = self._Ea_g.value_si        # J/mol
        Eae  = self._Ea_e.value_si        # J/mol

        return A_si * (Te / T0) ** bval * exp(-Eag / (constants.R * T)) * exp(Eae * (Te - T) / (constants.R * T * Te))

    cpdef double get_rate_coefficient(self, double T, double P=0.0) except -1:
        """
        Standard RMG interface: interpret this as k(T, Te=T).
        """
        return self.get_rate_coefficient_two_temp(T, T)

    cpdef change_t0(self, double T0):
        """
        Changes the reference temperature used in the exponent to `T0` in K,
        and adjusts the preexponential factor accordingly.
        """
        self._A.value_si /= (self._T0.value_si / T0) ** self._n.value_si
        self._T0.value_si = T0

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Check equality with another kinetics object.
        """
        if not isinstance(other_kinetics, TwoTemperaturePlasma):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if not self.A.equals(other_kinetics.A):
            return False
        if not self.n.equals(other_kinetics.n):
            return False
        if not self.Ea_g.equals(other_kinetics.Ea_g):
            return False
        if not self.Ea_e.equals(other_kinetics.Ea_e):
            return False
        # T0 scales the rate through (Te/T0)**n, so two objects that differ only in
        # their reference temperature are different rates. Compared through
        # ScalarQuantity.equals, as Arrhenius compares its own T0.
        if not self.T0.equals(other_kinetics.T0):
            return False
        return True

    cpdef change_rate(self, double factor):
        """
        Scale A by `factor`.
        """
        self._A.value_si *= factor

    def to_cantera_kinetics(self):
        """
        Convert to a Cantera TwoTempPlasmaRate object.

        Returns
        -------
        ct.TwoTempPlasmaRate
            With A in m^3/kmol/s, Ea_g and Ea_e in J/kmol.
        """
        import cantera as ct

        rate_units_dimensionality = {
            '1/s': 0,
            's^-1': 0,
            'm^3/(mol*s)': 1,
            'm^6/(mol^2*s)': 2,
            'cm^3/(mol*s)': 1,
            'cm^6/(mol^2*s)': 2,
            'm^3/(molecule*s)': 1,
            'm^6/(molecule^2*s)': 2,
            'cm^3/(molecule*s)': 1,
            'cm^6/(molecule^2*s)': 2,
        }

        A = self._A.value_si
        try:
            # convert from per-mol to per-kmol (and any cm/molecule variants)
            A *= 1000 ** rate_units_dimensionality[self._A.units]
        except KeyError:
            raise Exception(
                'TwoTemperaturePlasma A-factor units {0} not supported for '
                'conversion to Cantera TwoTempPlasmaRate.'.format(self._A.units)
            )

        b = self._n.value_si
        Ea_g = self._Ea_g.value_si * 1000.0   # J/mol -> J/kmol
        Ea_e = self._Ea_e.value_si * 1000.0   # J/mol -> J/kmol

        return ct.TwoTempPlasmaRate(A, b, Ea_g, Ea_e)

    def set_cantera_kinetics(self, ct_reaction, species_list):
        """
        Assign a Cantera Reaction's rate to a TwoTempPlasmaRate built
        from this object.
        """
        import cantera as ct

        assert isinstance(ct_reaction.rate, ct.TwoTempPlasmaRate), \
            "ct_reaction.rate must be a Cantera TwoTempPlasmaRate"
        ct_reaction.rate = self.to_cantera_kinetics()

    cpdef ArrheniusEP to_arrhenius_ep(self, double alpha=0.0, double dHrxn=0.0):
        """
        Converts an Arrhenius object to ArrheniusEP
        This function assumes you want to convert the electron activation energy Ea_e,
        not the gas activation energy Ea_g.

        If setting alpha, you need to also input dHrxn, which must be given
        in J/mol (and vise versa).
        """
        if bool(alpha) ^ bool(dHrxn):
            raise Exception('If you set alpha or dHrxn in to_arrhenius_ep, you need to set the other value to non-zero.')
        self.change_t0(1)
        aep = ArrheniusEP(A=self.A,
                          n=self.n,
                          alpha=alpha,
                          E0=(self.Ea_e.value_si - alpha * dHrxn, 'J/mol'),
                          Tmin=self.Tmin,
                          Tmax=self.Tmax,
                          Pmin=self.Pmin,
                          Pmax=self.Pmax,
                          uncertainty=self.uncertainty,
                          solute=self.solute,
                          comment=self.comment)
        return aep


################################################################################


cdef class ElectronCollisionPlasma(KineticsModel):
    """
    Te-only plasma kinetics based on a tabulated electron–collision
    cross-section σ(E) on an energy grid.

    The rate coefficient is the Maxwellian average:
        k(Te) = ⟨σ v⟩(Te)
    calculated by numerical integration over the stored cross-section.
    """

    def __init__(self,
                 energies=None,
                 sigma=None,
                 Tmin=None, Tmax=None, Pmin=None, Pmax=None,
                 uncertainty=None, solute=None, comment=''):
        KineticsModel.__init__(
            self,
            Tmin=Tmin, Tmax=Tmax,
            Pmin=Pmin, Pmax=Pmax,
            uncertainty=uncertainty,
            solute=solute,
            comment=comment,
        )

        if energies is None:
            self.energies = ([], "J/mol")
        else:
            self.energies = energies

        if sigma is None:
            self.sigma = ([], "m^2")
        else:
            self.sigma = sigma

        self.uses_electron_temperature = True

    def __repr__(self):
        string = "ElectronCollisionPlasma(energies={0!r}, sigma={1!r}".format(
            self.energies, self.sigma
        )
        if self.Tmin is not None: string += ", Tmin={0!r}".format(self.Tmin)
        if self.Tmax is not None: string += ", Tmax={0!r}".format(self.Tmax)
        if self.Pmin is not None: string += ", Pmin={0!r}".format(self.Pmin)
        if self.Pmax is not None: string += ", Pmax={0!r}".format(self.Pmax)
        if self.comment: string += ', comment="""{0}"""'.format(self.comment)
        string += ")"
        return string

    def __reduce__(self):
        return (ElectronCollisionPlasma,
                (self.energies, self.sigma,
                 self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                 None, None, self.comment))

    property energies:
        def __get__(self): return self._energies
        def __set__(self, value): self._energies = quantity.Energy(value)

    property sigma:
        def __get__(self): return self._sigma
        def __set__(self, value): self._sigma = quantity.Quantity(value)

    def __getstate__(self):
        return {
            "energies": self.energies,
            "sigma": self.sigma,
            "Tmin": self.Tmin, "Tmax": self.Tmax,
            "Pmin": self.Pmin, "Pmax": self.Pmax,
            "solute": self.solute, "comment": self.comment,
        }

    def __setstate__(self, state):
        self.__init__(**state)

    # --------------------------------------------------------------------------
    #    Core Physics (Numerical Integration)
    # --------------------------------------------------------------------------

    cpdef double integrate_rate_coefficient(self, double Te) except -1:
        """
        Calculates k(Te) = <sigma * v> by numerically integrating the stored
        cross-section over a Maxwellian electron energy distribution.
        """
        if Te <= 0.0:
            return 0.0

        # RMG stores Energy in J/mol. We need J/particle for the physics integration.
        cdef np.ndarray[np.float64_t, ndim=1] E_molar = self._energies.value_si
        cdef np.ndarray[np.float64_t, ndim=1] sigma_grid = self._sigma.value_si
        cdef int n = E_molar.shape[0]

        if n < 2:
            return 0.0

        cdef double kB = constants.kB
        cdef double Na = constants.Na
        cdef double m_e = constants.m_e  # kg

        # Pre-factor for Maxwellian flux
        cdef double prefactor = sqrt(8.0 / (constants.pi * m_e)) * pow(kB * Te, -1.5)

        cdef double integral = 0.0
        cdef double E1, E2, s1, s2, term1, term2, dE
        cdef double kTe = kB * Te
        cdef int i

        # Loop variables in Particle Units (Joules)
        cdef double E1_part, E2_part

        for i in range(n - 1):
            # Convert J/mol -> J/particle
            E1_part = E_molar[i] / Na
            E2_part = E_molar[i + 1] / Na

            s1 = sigma_grid[i]
            s2 = sigma_grid[i + 1]

            term1 = s1 * E1_part * exp(-E1_part / kTe)
            term2 = s2 * E2_part * exp(-E2_part / kTe)
            dE = E2_part - E1_part

            integral += 0.5 * (term1 + term2) * dE

        # Result is m^3/s per particle. Convert to m^3/(mol*s) for RMG output.
        return (prefactor * integral) * Na

    cpdef double get_rate_coefficient_electron_temp(self, double Te) except -1:
        """
        Returns k(Te) in m^3/(mol*s) via numerical integration.
        """
        return self.integrate_rate_coefficient(Te)

    cpdef double get_rate_coefficient(self, double T, double P=0.0) except -1:
        """
        Standard RMG interface. Interprets T as electron temperature.
        """
        return self.integrate_rate_coefficient(T)

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if the energy grid and cross-section table match (and the
        validity window agrees), ``False`` otherwise.

        The arrays are compared exactly rather than through ``ArrayQuantity.equals``:
        that helper accepts anything within 1% *or* an absolute 0.01, and a collision
        cross-section is of order 1e-20 m^2, so every table would compare equal to
        every other one. This mirrors how Chebyshev compares its coefficients.
        """
        if not isinstance(other_kinetics, ElectronCollisionPlasma):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if not np.array_equal(self._energies.value_si, other_kinetics._energies.value_si):
            return False
        if not np.array_equal(self._sigma.value_si, other_kinetics._sigma.value_si):
            return False
        return True

    cpdef change_rate(self, double factor):
        """
        Scale the overall rate by multiplying σ(E) by `factor`.
        """
        self._sigma.value_si *= factor

    # --------------------------------------------------------------------------
    #    Cantera & Legacy Conversion
    # --------------------------------------------------------------------------

    def to_cantera_kinetics(self):
        """
        Convert to a Cantera ElectronCollisionPlasmaRate object.
        Passes raw energy/cross-section arrays for internal Cantera integration.
        """
        import cantera as ct

        if hasattr(ct, "ElectronCollisionPlasmaRate"):
            # Cantera 3.0+ native support
            # Energies: Convert J/mol (RMG SI) -> eV (Cantera input)
            # Factor: Faraday constant = 96485.33212 J/mol per eV
            conversion_factor = 96485.33212

            # Use .value_si (numpy array) and divide manually.
            # Do NOT use .value("eV") as .value is not callable in Cython classes.
            energies_eV = self._energies.value_si / conversion_factor

            # Sigma: m^2 (already correct in value_si)
            return ct.ElectronCollisionPlasmaRate(
                energies_eV,
                self._sigma.value_si
            )
        else:
            # Fallback for older Cantera versions
            import warnings
            warnings.warn("Cantera.ElectronCollisionPlasmaRate not found. "
                          "Falling back to Arrhenius fit (TwoTempPlasma).")
            return self.to_two_temp_plasma().to_cantera_kinetics()

    def set_cantera_kinetics(self, ct_reaction, species_list):
        """
        Sets a Cantera Reaction's rate to an ElectronCollisionPlasmaRate
        built from this object.
        """
        ct_reaction.rate = self.to_cantera_kinetics()

    def to_arrhenius(self, double Tmin=0.0, double Tmax=0.0):
        """
        Fit the integrated cross-section data to an Arrhenius model k(Te).
        """
        if Tmin == 0.0: Tmin = 11600.0  # ~1 eV
        if Tmax == 0.0: Tmax = 116000.0  # ~10 eV

        cdef np.ndarray Te_list = 1.0 / np.linspace(1.0 / Tmax, 1.0 / Tmin, 50)
        cdef np.ndarray k_list = np.zeros_like(Te_list)

        for i, Te in enumerate(Te_list):
            k_list[i] = self.integrate_rate_coefficient(Te)

        arr = Arrhenius()
        # k_list is already molar from integrate_rate_coefficient
        arr.fit_to_data(Te_list, k_list, "m^3/(mol*s)", T0=1.0)

        arr.Tmin = (Tmin, "K")
        arr.Tmax = (Tmax, "K")
        arr.comment = f"Arrhenius fit to ElectronCollisionPlasma. {self.comment}"
        return arr

    cpdef TwoTemperaturePlasma to_two_temp_plasma(self):
        """
        Convert to TwoTemperaturePlasma by fitting k(Te).
        """
        arr = self.to_arrhenius()

        return TwoTemperaturePlasma(
            A=(arr.A.value_si, "m^3/(mol*s)"),
            n=float(arr.n.value_si),
            Ea_g=(0.0, "J/mol"),
            Ea_e=(arr.Ea.value_si, "J/mol"),
            Tmin=arr.Tmin,
            Tmax=arr.Tmax,
            comment=f"Mapped from ElectronCollisionPlasma. {self.comment}"
        )


################################################################################


cdef class BadnellRRArrhenius(KineticsModel):
    """
    Radiative recombination kinetics using the Badnell (2006) fit, evaluated at electron temperature Te.

    Rate expression (per-particle form in the paper):
        alpha_RR(Te) = A * [ sqrt(Te/T0) * (1 + sqrt(Te/T0))^(1 - B*)
                              * (1 + sqrt(Te/T1))^(1 + B*) ]^(-1)
        with B* = B + C * exp(-T2 / Te)  (use B* = B if C/T2 are not provided)

    This class returns the **per-mole** rate coefficient expected by RMG:
        k(Te) = alpha_RR(Te) * N_A        (for bimolecular RR)
    and in SI units according to the stored A units.

    Attributes (ScalarQuantity unless noted):
        A   : pre-exponential (supports 'cm^3/(molecule*s)', 'cm^3/(mol*s)', 'm^3/(mol*s)', etc.)
        B   : dimensionless Badnell parameter
        T0  : temperature in K
        T1  : temperature in K
        C   : (optional) dimensionless Badnell parameter
        T2  : (optional) temperature in K
        Z   : (optional) nuclear charge (int)
        N   : (optional) electron count before recombination (int)
        yaml_path_or_obj : (optional) source YAML file to populate this object
        comment, Tmin, Tmax, Pmin, Pmax (inherited)
    """

    def __init__(self,
                 A=None,
                 B=0.0,
                 T0=(1.0, "K"),
                 T1=(1.0, "K"),
                 C=None,
                 T2=None,
                 Z=None, N=None, yaml_path_or_obj=None,
                 Tmin=None, Tmax=None, Pmin=None, Pmax=None,
                 uncertainty=None, solute=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax,
                               uncertainty=uncertainty, solute=solute, comment=comment)
        self._Ea = None

        # Radiative recombination consumes an electron, so both channels are
        # intrinsic to this rate law regardless of how the object was built. Set
        # them here, ahead of every branch below, so that no early return can
        # leave a table-built object reporting False for either.
        self.uses_electron_temperature = True
        self.uses_electron_density = True

        # If Z/N are given, load from YAML and return early
        if Z is not None and N is not None:
            yaml_path_or_obj = yaml_path_or_obj or os.path.join(settings['database.directory'], 'kinetics', 'badnell.yaml')
            try:
                Zi = int(Z)
                Ni = int(N)
            except Exception:
                raise TypeError("Z and N must be integers.")
            self.populate_from_yaml(yaml_path_or_obj, Zi, Ni,
                                    Tmin=Tmin, Tmax=Tmax, comment=comment)
            return

        self.A = A
        self.B = B
        self.T0 = T0
        self.T1 = T1
        self.C = C
        self.T2 = T2

    property Ea:
        """
        The activation energy, estimated by fitting the kinetics to an Arrhenius form.
        Stored internally as _Ea.
        """
        def __get__(self):
            if self._Ea is None:
                self._Ea = self.to_arrhenius().Ea
            return self._Ea

    def __repr__(self):
        string = 'BadnellRRArrhenius(A={0!r}, B={1!r}, T0={2!r}, T1={3!r}'.format(self.A, self.B, self.T0, self.T1)
        if self.C is not None: string += ', C={0!r}'.format(self.C)
        if self.T2 is not None: string += ', T2={0!r}'.format(self.T2)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.uncertainty: string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.solute: string += ', solute={0!r}'.format(self.solute)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        return (BadnellRRArrhenius,
                (self.A, self.B, self.T0, self.T1, self.C, self.T2,
                 None, None, None,
                 self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                 self.uncertainty, self.solute, self.comment))

    # -------- properties --------

    property A:
        """Pre-exponential factor (Badnell A)."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            # Allow per-molecule or per-mole units; quantity.RateCoefficient handles units.
            self._A = quantity.RateCoefficient(value)

    property B:
        """Badnell B parameter (dimensionless)."""
        def __get__(self):
            return self._B
        def __set__(self, value):
            self._B = quantity.Dimensionless(value)

    property T0:
        """Badnell T0 parameter (K)."""
        def __get__(self):
            return self._T0
        def __set__(self, value):
            self._T0 = quantity.Temperature(value)

    property T1:
        """Badnell T1 parameter (K)."""
        def __get__(self):
            return self._T1
        def __set__(self, value):
            self._T1 = quantity.Temperature(value)

    property C:
        """Badnell C parameter (dimensionless, optional)."""
        def __get__(self):
            return self._C
        def __set__(self, value):
            if value is None:
                self._C = None
            else:
                self._C = quantity.Dimensionless(value)

    property T2:
        """Badnell T2 parameter (K, optional)."""
        def __get__(self):
            return self._T2
        def __set__(self, value):
            if value is None:
                self._T2 = None
            else:
                self._T2 = quantity.Temperature(value)


    cpdef double get_rate_coefficient_electron_temp(self, double Te) except -1:
        """
        Return k(Te) in SI units of m^3/(mol*s) at electron temperature Te [K].

        The Badnell radiative-recombination fit is a function of the electron
        temperature only; this method is the explicit, unambiguous entry point
        for two-temperature reactors.
        """
        return self.get_rate_coefficient(Te)

    cpdef double get_rate_coefficient(self, double T, double P=0.0) except -1:
        """
        Return k(T) in SI units of m^3/(mol*s).

        The Badnell fit is evaluated using A in SI (m^3/(mol*s)), so this
        method is agnostic to whether the user originally provided A in
        cm^3/(mol*s), cm^3/(molecule*s), or m^3/(mol*s).
        """
        cdef double A_SI, B, T0, T1, C, T2
        cdef double t0, t1, Bstar, denom

        if T <= 0.0:
            return 0.0

        # 1) Pre-exponential in SI
        A_SI = self.A.value_si  # m^3/(mol*s)

        # 2) Dimensionless / temperature parameters
        B = <double> (self.B.value_si if hasattr(self.B, "value_si") else self.B)

        C = 0.0
        if self.C is not None:
            C = <double> (self.C.value_si if hasattr(self.C, "value_si") else self.C)

        if hasattr(self.T0, "value_si"):
            T0 = <double> self.T0.value_si
        else:
            T0 = <double> self.T0

        if hasattr(self.T1, "value_si"):
            T1 = <double> self.T1.value_si
        else:
            T1 = <double> self.T1

        T2 = 0.0
        if self.T2 is not None:
            if hasattr(self.T2, "value_si"):
                T2 = <double> self.T2.value_si
            else:
                T2 = <double> self.T2

        # 3) Badnell shape factor (per Badnell 2006)
        t0 = sqrt(T / T0)
        t1 = sqrt(T / T1)

        Bstar = B
        if C != 0.0 and T2 > 0.0:
            # B* = B + C * exp(-T2 / T)
            Bstar = B + C * exp(-T2 / T)

        denom = t0 * pow(1.0 + t0, 1.0 - Bstar) * pow(1.0 + t1, 1.0 + Bstar)

        # k(T) in m^3/(mol*s)
        return A_SI / denom

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        if not isinstance(other_kinetics, BadnellRRArrhenius):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False

        # Compare all params (handle optional C/T2)
        if not self.A.equals(other_kinetics.A): return False
        if not self.B.equals(other_kinetics.B): return False
        if not self.T0.equals(other_kinetics.T0): return False
        if not self.T1.equals(other_kinetics.T1): return False

        if (self.C is None) ^ (other_kinetics.C is None): return False
        if (self.T2 is None) ^ (other_kinetics.T2 is None): return False
        if self.C is not None and not self.C.equals(other_kinetics.C): return False
        if self.T2 is not None and not self.T2.equals(other_kinetics.T2): return False
        return True

    cpdef change_rate(self, double factor):
        """
        Multiply the Badnell A factor by 'factor'.
        (This scales the entire alpha_RR, analogous to Arrhenius.change_rate.)
        """
        self._A.value_si *= factor

    def _brra__extract_row(self, obj, int Z, int N):
        """
        Internal: find a (Z,N) row in a few common YAML layouts.
        Returns a Python dict with keys A,B,T0,T1 and optional C,T2.
        Raises KeyError if not found.
        """
        cdef object units = obj.get("units", {})
        # Try list-of-blocks schema
        coeffs = obj.get("coefficients")
        if isinstance(coeffs, list):
            for blk in coeffs:
                zblk = blk.get("Z")
                if zblk is not None and int(zblk) == Z:
                    entries = blk.get("entries")
                    if isinstance(entries, list):
                        for e in entries:
                            if int(e.get("N", -999)) == N:
                                return e
                    elif isinstance(entries, dict):
                        # entries: { "N": {A:...,B:...}, ... }
                        if str(N) in entries:
                            return entries[str(N)]
                        if N in entries:
                            return entries[N]
        # Try nested map schema under "data"
        data = obj.get("data")
        if isinstance(data, dict):
            znode = data.get(str(Z)) if str(Z) in data else data.get(Z)
            if isinstance(znode, dict):
                node = znode.get(str(N)) if str(N) in znode else znode.get(N)
                if isinstance(node, dict):
                    return node
        # Try top-level Z map
        znode = obj.get(str(Z)) if str(Z) in obj else obj.get(Z)
        if isinstance(znode, dict):
            node = znode.get(str(N)) if str(N) in znode else znode.get(N)
            if isinstance(node, dict):
                return node
        raise KeyError(f"Badnell YAML: no entry for Z={Z}, N={N}")

    def _brra__units(self, obj):
        """Internal: return (a_units, t_units) with sensible defaults."""
        units = obj.get("units", {}) if isinstance(obj, dict) else {}
        a_units = units.get("A", "cm^3/(molecule*s)")
        t_units = units.get("T", "K")
        return a_units, t_units

    def _brra__compute_default_T_window(self, int Z, int N):
        """
        Badnell fit validity spans ~ z^2 * [1e1, 1e7] K with z = Z-N (initial charge).
        Use that for Tmin/Tmax if user left them unset.
        """
        cdef int z = Z - N
        cdef double Tmin = 10.0 * z * z
        cdef double Tmax = 1.0e7 * z * z
        if z < 0:
            z = 0  # defensive; but physically you shouldn't pass N>Z
        # If z==0, fall back to a broad neutral window (let RMG handle validity checks)
        if z == 0:
            Tmin, Tmax = 1.0, 3.0e9
        return Tmin, Tmax

    cpdef populate_from_yaml(self, object yaml_path_or_obj, int Z, int N,
                             bint allow_Z_gt36=False, Tmin=None, Tmax=None, comment=None):
        """
        Populate this BadnellRRArrhenius from a YAML dataset keyed by (Z,N).

        Parameters
        ----------
        yaml_path_or_obj : str | pathlib.Path | file-like | dict
            YAML filepath or already-loaded dict.
        Z : int
            Nuclear charge (1..36 supported by default).
        N : int
            Electron count BEFORE recombination (so z = Z-N).
        allow_Z_gt36 : bool
            If False, raise for Z>36.
        Tmin, Tmax : optional
            Override validity window (K). If None, uses ~z^2*[1e1,1e7] K.
        comment : optional
            Override/append comment.

        Notes
        -----
        - A is interpreted with units from YAML (default 'cm^3/(molecule*s)').
        - T0,T1,(T2) are in K.
        - Sets C,T2 only if BOTH are present; else uses B* = B.
        """
        if not allow_Z_gt36 and Z > 36:
            raise ValueError(f"Badnell YAML: Z={Z} exceeds 36 (this loader is restricted to Z<=36).")

        # Load YAML if needed
        cdef dict data
        if isinstance(yaml_path_or_obj, dict):
            data = yaml_path_or_obj
        else:
            import yaml as _yaml
            with open(yaml_path_or_obj, "r") as f:
                data = _yaml.safe_load(f)

        row = self._brra__extract_row(data, Z, N)
        a_units, t_units = self._brra__units(data)

        # Required
        A = float(row["A"])
        B = float(row["B"])
        T0 = float(row["T0"])
        T1 = float(row["T1"])

        # Optional C,T2 (use only if both provided)
        C = row.get("C", None)
        T2 = row.get("T2", None)
        if C is not None:
            C = float(C)
        if T2 is not None:
            T2 = float(T2)

        # Assign to this object
        self.A  = (A, a_units)
        self.B  = B
        self.T0 = (T0, t_units)
        self.T1 = (T1, t_units)
        if C is not None and T2 is not None:
            self.C  = C
            self.T2 = (T2, t_units)
        else:
            self.C  = None
            self.T2 = None

        # Validity window (unless user overrides)
        if Tmin is None or Tmax is None:
            dTmin, dTmax = self._brra__compute_default_T_window(Z, N)
            if Tmin is None: Tmin = dTmin
            if Tmax is None: Tmax = dTmax
        self.Tmin = (Tmin, "K")
        self.Tmax = (Tmax, "K")

        # Comment
        base = f"Badnell (2006) RR fit, Z={Z}, N={N}"
        self.comment = f"{base}; {comment}" if comment else base

        return self

    @classmethod
    def from_yaml(cls, object yaml_path_or_obj, int Z, int N,
                  bint allow_Z_gt36=False, Tmin=None, Tmax=None, comment=None):
        """
        Construct and return a new BadnellRRArrhenius populated from YAML.
        """
        obj = cls(A=(1.0e-12, "cm^3/(molecule*s)"), B=0.0, T0=(1.0, "K"), T1=(1.0, "K"))
        obj.populate_from_yaml(yaml_path_or_obj, Z, N,
                               allow_Z_gt36=allow_Z_gt36,
                               Tmin=Tmin, Tmax=Tmax, comment=comment)
        return obj

    def to_arrhenius(self, double Tmin=0.0, double Tmax=0.0):
        """
        Return an Arrhenius object that fits the Badnell kinetics over the specified temperature range.

        The fit is performed by evaluating the Badnell rate at intervals linear in 1/T
        (standard Arrhenius plot spacing) and fitting the modified Arrhenius equation to those points.

        If Tmin or Tmax are not provided (or 0.0), the method defaults to self.Tmin/self.Tmax.
        If those are also not set, it defaults to 800 K - 3000 K.
        """
        # Determine Temperature Boundaries
        if Tmin == 0.0:
            if self.Tmin is not None:
                Tmin = self.Tmin.value_si
            else:
                Tmin = 1.0e4

        if Tmax == 0.0:
            if self.Tmax is not None:
                Tmax = self.Tmax.value_si
            else:
                Tmax = 1.0e6

        if Tmin >= Tmax:
            raise ValueError(f"Tmin ({Tmin}) must be less than Tmax ({Tmax}) for Arrhenius fitting.")

        # Generate sampling points
        # Linear in 1/T space gives better weighting for Arrhenius fits
        # We use 50 points to ensure a smooth fit
        cdef np.ndarray Tlist = 1.0 / np.linspace(1.0 / Tmax, 1.0 / Tmin, 50)
        cdef np.ndarray klist = np.zeros_like(Tlist)

        cdef int i
        cdef double T_val

        for i in range(len(Tlist)):
            T_val = Tlist[i]
            klist[i] = self.get_rate_coefficient(T_val)

        # Determine units
        # get_rate_coefficient returns SI Molar: m^3/(mol*s)
        # We assume bimolecular for standard Badnell RR
        cdef str kunits = "m^3/(mol*s)"

        # Create and fit Arrhenius object
        # T0=1.0 is standard for minimizing correlation between A and n
        arr = Arrhenius()
        arr.fit_to_data(Tlist, klist, kunits, T0=1.0)

        # Carry over metadata
        arr.Tmin = (Tmin, "K")
        arr.Tmax = (Tmax, "K")
        arr.comment = f"Fitted to BadnellRRArrhenius over range {Tmin}-{Tmax} K. Original comment: {self.comment}"

        return arr

    def to_chebyshev(self, double Tmin=0.0, double Tmax=0.0, int degree_t=10):
        """
        Convert the Badnell kinetics to a Chebyshev object.

        Since standard Chemkin/Cantera formats do not support the Badnell function natively,
        the most robust, high-accuracy way to export this kinetics model is to map it
        to a Chebyshev polynomial. We use a Pressure-Independent Chebyshev fit (P_basis=1).

        Parameters:
            degree_t: Number of temperature coefficients (basis functions).
                      Badnell curves can be complex, so a higher degree (default 10)
                      is recommended compared to standard Arrhenius (typically 3-4).
        """
        from rmgpy.kinetics.chebyshev import Chebyshev

        # 1. Determine Temperature Range
        if Tmin == 0.0:
            Tmin = self.Tmin.value_si if self.Tmin is not None else 10.0
        if Tmax == 0.0:
            Tmax = self.Tmax.value_si if self.Tmax is not None else 1.0e7

        if Tmin >= Tmax:
            raise KineticsError(
                f"BadnellRRArrhenius.to_chebyshev: Tmin ({Tmin}) must be < Tmax ({Tmax})."
            )

        # 2. Define Dummy Pressure Range (Required for Chebyshev format)
        # Since Badnell is P-independent, these values don't affect the rate
        # as long as we fit with degree_p = 1.
        cdef double Pmin = 100.0        # Pa  (≈ 0.001 bar)
        cdef double Pmax = 1.0e7        # Pa  (≈ 100 bar)

        # 3. Create the Chebyshev Object
        cheb = Chebyshev(
            Tmin=(Tmin, "K"), Tmax=(Tmax, "K"),
            Pmin=(Pmin, "Pa"), Pmax=(Pmax, "Pa"),
        )

        # 4. Generate Grid and Fit
        # We need MORE grid points than polynomial degrees in both T and P.
        cdef int nT = degree_t + 1      # > degree_t
        cdef int degree_p = 1
        cdef int nP = 2                 # > degree_p

        # Use Chebyshev roots in 1/T space for optimal interpolation node spacing
        k_idx = np.arange(nT, dtype=float)

        invT_mid = 0.5 * (1.0 / Tmax + 1.0 / Tmin)
        invT_halfspan = 0.5 * (1.0 / Tmax - 1.0 / Tmin)

        T_nodes = invT_mid + invT_halfspan * np.cos(
            (2.0 * k_idx + 1.0) * np.pi / (2.0 * nT)
        )
        T_points = 1.0 / T_nodes  # Convert 1/T back to T

        # P_points: two points, but the rate is P-independent
        P_points = np.array([Pmin, Pmax], dtype=float)

        # K_data: shape (nT, nP)
        K_data = np.zeros((nT, nP), dtype=float)
        for i, T in enumerate(T_points):
            # Calculate Badnell rate at this T
            k_val = self.get_rate_coefficient(T)
            # Same k(T) for all pressures since the rate is P-independent
            K_data[i, :] = k_val

        # 5. Perform Fit
        cheb.fit_to_data(
            T_points,
            P_points,
            K_data,
            "m^3/(mol*s)",
            degree_t,
            degree_p,
            Tmin,
            Tmax,
            Pmin,
            Pmax,
        )

        # Add original kinetics description to comment for traceability
        cheb.comment = (
            f"Chebyshev fit to BadnellRR (Z={getattr(self, 'Z', '?')}, "
            f"N={getattr(self, 'N', '?')}). Original: {self}"
        )
        return cheb

    cpdef TwoTemperaturePlasma to_two_temp_plasma(self):
        """
        Map this BadnellRRArrhenius (k = k(Te) only) onto a TwoTemperaturePlasma.

        The mapping is constructed so that:
          - along T = Te, k(T, Te) matches the Badnell rate (via Arrhenius fit),
          - for general T, the rate depends only on Te (T cancels analytically when Ea_g = Ea_e = fitted Ea).
        """
        if self.Tmin is not None:
            Tlo = self.Tmin.value_si
        else:
            Tlo = 1.0e4

        if self.Tmax is not None:
            Thi = self.Tmax.value_si
        else:
            Thi = 1.0e6
        #    to_arrhenius() should already pick a physically reasonable Te window.
        arr = self.to_arrhenius(Tmin=Tlo, Tmax=Thi)

        # Convert n to a bare float if needed
        try:
            n_val = float(arr.n.value_si)
        except AttributeError:
            n_val = float(arr.n)

        # Build the canonical TwoTemperaturePlasma with Ea_g = 0, Ea_e = Ea
        plasma = TwoTemperaturePlasma(
            A=(arr.A.value_si, "m^3/(mol*s)"),
            n=n_val,
            Ea_g=(0.0, "J/mol"),
            Ea_e=(arr.Ea.value_si, "J/mol"),
            Tmin=(Tlo, "K"),
            Tmax=(Thi, "K"),
            Pmin=self.Pmin,
            Pmax=self.Pmax,
            uncertainty=self.uncertainty,
            solute=self.solute,
            comment=(f"BadnellRRArrhenius mapped to TwoTemperaturePlasma over {Tlo:g}-{Thi:g} K (Te). Original: {self.comment}"
            ),
        )

        return plasma

    def to_cantera_kinetics(self):
        """
        Export this BadnellRRArrhenius as a Cantera TwoTempPlasmaRate.

        Implementation:
          BadnellRRArrhenius → TwoTemperaturePlasma → ct.TwoTempPlasmaRate
        """
        plasma = self.to_two_temp_plasma()
        return plasma.to_cantera_kinetics()


################################################################################


cdef class VoronovEIArrhenius(KineticsModel):
    """
    Electron-impact ionization (Voronov, 1997, https://doi.org/10.1006/adnd.1997.0732) for ground-state atoms/ions.

    Per-particle fit (Table I in Voronov 1997):
        k = <sigma v>(Te_eV) = A * [ (1 + P*sqrt(U)) * U^K * exp(-U) ] / [ (X + U) ]
        where U = dE / Te_eV, Te_eV = (8.617333262e-5 eV/K) * Te_K
        Here, (Te_eV) means "function of", not multiplication.

    This class returns **molar** k in SI:
        k(Te) = <sigma v>(Te) * N_A   [m^3/(mol*s)]

    YAML schema (voronov.yaml):
      units:
        A:  "cm^3/(molecule*s)"
        dE: "eV"
        T:  "K"
        Tmin: "eV"
        Tmax: "eV"
      coefficients:
        - Z: <int>
          element: "<sym>"
          entries:
            - N: <int>
              dE, P, A, X, K, Tmin, Tmax

    Attributes:
        A : RateCoefficient (accepts per-molecule or per-mole units)
        P,X,K: Dimensionless
        dE: ionization threshold in eV (stored as double)
        Z,N: integers identifying the stage (N = electrons before ionization)
    """

    def __init__(self,
                 A=None,
                 P=0.0,
                 X=0.0,
                 K=0.0,
                 dE=None,
                 Z=None, N=None, yaml_path_or_obj=None,
                 Tmin=None, Tmax=None, Pmin=None, Pmax=None,
                 uncertainty=None, solute=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax,
                               uncertainty=uncertainty, solute=solute, comment=comment)
        self._Ea = None

        # Electron-impact ionization is driven by the electron population, so both
        # channels are intrinsic to this rate law regardless of how the object was
        # built. Set them here, ahead of every branch below, so that no early return
        # can leave a table-built object reporting False for either.
        self.uses_electron_temperature = True
        self.uses_electron_density = True

        # If Z/N are provided, load from YAML now
        if Z is not None and N is not None:
            yaml_path_or_obj = yaml_path_or_obj or os.path.join(settings['database.directory'], 'kinetics', 'voronov.yaml')
            try:
                Zi = int(Z); Ni = int(N)
            except Exception:
                raise TypeError("Z and N must be integers.")
            self.populate_from_yaml(yaml_path_or_obj, Zi, Ni, Tmin=Tmin, Tmax=Tmax, comment=comment)
            return

        # direct-parameter construction
        if dE is None:
            raise ValueError(
                "VoronovEIArrhenius requires the ionization threshold dE (in eV) when "
                "constructed from explicit parameters. The Voronov rate goes as "
                "exp(-dE / Te_eV), so a substituted threshold silently rescales every "
                "rate computed from this object; there is no defensible default. Pass "
                "dE, or build from the Voronov table by giving Z and N."
            )
        self.A = A if A is not None else (1.0e-12, "cm^3/(molecule*s)")
        self.P = P
        self.X = X
        self.K = K
        self.dE = dE

    def __repr__(self):
        # dE is emitted as a (value, units) keyword at full precision. An RMG database
        # entry is saved as repr(entry.data) and read back by evaluating it, so this
        # has to be a loadable call -- a bare "dE=13.6 eV" suffix is neither a keyword
        # argument nor valid Python -- and the threshold has to survive the trip
        # exactly, since the rate goes as exp(-dE / Te_eV). The stored value is in eV
        # by construction and the dE setter reads the value out of the tuple.
        string = "VoronovEIArrhenius(A={0!r}, P={1!r}, X={2!r}, K={3!r}, dE=({4!r}, 'eV')".format(
            self.A, self.P, self.X, self.K, self._dE_eV)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.uncertainty: string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.solute: string += ', solute={0!r}'.format(self.solute)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        return (VoronovEIArrhenius, (
            self.A, self.P, self.X, self.K, self._dE_eV,  # dE
            None, None, None,                             # Z, N, yaml_path_or_obj
            self.Tmin, self.Tmax,                         # Tmin, Tmax
            self.Pmin, self.Pmax,                         # Pmin, Pmax
            self.uncertainty, self.solute, self.comment
        ))

    # ---- properties ----
    property A:
        def __get__(self):
            return self._A
        def __set__(self, value):
            # Normalize to per-mol before creating RateCoefficient
            if isinstance(value, (tuple, list)) and len(value) >= 2:
                val = float(value[0])
                units = str(value[1])
                u = units.replace(" ", "").lower()
                if "molecule" in u or u.endswith("cm^3/s") or u == "cm^3/s":
                    # per particle -> per mole
                    val *= constants.Na
                    units = "cm^3/(mol*s)"
                value = (val, units)
            self._A = quantity.RateCoefficient(value)

    property P:
        def __get__(self):
            return self._P
        def __set__(self, value):
            self._P = quantity.Dimensionless(value)

    property X:
        def __get__(self):
            return self._X
        def __set__(self, value):
            self._X = quantity.Dimensionless(value)

    property K:
        def __get__(self):
            return self._K
        def __set__(self, value):
            self._K = quantity.Dimensionless(value)

    property dE:
        """Ionization threshold in eV (stored as a plain double in eV)."""
        def __get__(self): return self._dE_eV
        def __set__(self, value):
            if isinstance(value, tuple) or isinstance(value, list):
                # Accept (val, "eV")
                self._dE_eV = float(value[0])
            else:
                self._dE_eV = float(value)

    property dE_eV:
        def __get__(self):
            return self._dE_eV
        def __set__(self, val):
            self._dE_eV = float(val)

    # ---- core API ----
    cpdef double get_rate_coefficient_electron_temp(self, double Te) except -1:
        """
        Return bimolecular **molar** rate coefficient k(Te) in SI (m^3/(mol*s))
        at electron temperature Te [K].

        The Voronov electron-impact-ionization fit is a function of the
        electron temperature only; this method is the explicit, unambiguous
        entry point for two-temperature reactors.
        """
        return self.get_rate_coefficient(Te)

    cpdef double get_rate_coefficient(self, double T, double P=0.0) except -1:
        """
        Return bimolecular **molar** rate coefficient k(Te) in SI (m^3/(mol*s)) at electron temperature T [K].

        Uses Voronov (1997) Eq. (1):
            <σv> = A * [ (1 + P*sqrt(U)) / (X + U) ] * U^K * exp(-U)

        with U = dE / Te_eV.
        """
        cdef double A_si, Pval, Xval, Kval, Te_eV, U, numerator, denominator

        # 1. Physical Constants
        # Boltzmann constant in eV/K (approx 8.617e-5)
        cdef double kB_eVperK = constants.kB / constants.e

        if T <= 0.0:
            return 0.0

        Te_eV = kB_eVperK * T
        if Te_eV <= 0.0:
            return 0.0

        # 2. Calculate U (Dimensionless)
        # Ensure dE is in eV.
        U = self._dE_eV / Te_eV

        # Numerical stability for extremely high Te (U -> 0)
        if U < 1.0e-16:
            U = 1.0e-16

        # 3. Load Parameters
        A_si = self._A.value_si
        Pval = self._P.value_si
        Xval = self._X.value_si
        Kval = self._K.value_si

        # 4. Calculate Formula (Corrected)
        # Eq: A * [ (1 + P*sqrt(U)) / (X + U) ] * U^K * exp(-U)

        numerator = (1.0 + Pval * sqrt(U)) * pow(U, Kval) * exp(-U)
        denominator = (Xval + U)

        if denominator == 0.0:
            raise ZeroDivisionError("VoronovEIArrhenius: denominator (X + U) is zero.")

        cdef double rate_per_particle = A_si * (numerator / denominator)

        return rate_per_particle

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        if not isinstance(other_kinetics, VoronovEIArrhenius):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if not self.A.equals(other_kinetics.A): return False
        if not self.P.equals(other_kinetics.P): return False
        if not self.X.equals(other_kinetics.X): return False
        if not self.K.equals(other_kinetics.K): return False
        if self._dE_eV != other_kinetics._dE_eV: return False
        return True

    cpdef change_rate(self, double factor):
        """Scale A by `factor`."""
        self._A.value_si *= factor

    # ---- YAML helpers ----
    def _vor__extract_row(self, obj, int Z, int N):
        """
        Return the dict row for (Z,N) from a voronov.yaml-like object.
        """
        coeffs = obj.get("coefficients")
        if isinstance(coeffs, list):
            for blk in coeffs:
                zblk = blk.get("Z")
                if zblk is not None and int(zblk) == Z:
                    entries = blk.get("entries")
                    if isinstance(entries, list):
                        for e in entries:
                            if int(e.get("N", -999)) == N:
                                return e
        # Also allow top-level {Z:{N:{...}}} if ever used
        znode = obj.get(str(Z)) if isinstance(obj, dict) and str(Z) in obj else obj.get(Z) if isinstance(obj, dict) else None
        if isinstance(znode, dict):
            node = znode.get(str(N)) if str(N) in znode else znode.get(N)
            if isinstance(node, dict):
                return node
        raise KeyError(f"Voronov YAML: no entry for Z={Z}, N={N}")

    def _vor__units(self, obj):
        """
        Return (a_units, t_units, de_units, te_min_units, te_max_units).
        Defaults mirror the schema.
        """
        units = obj.get("units", {}) if isinstance(obj, dict) else {}
        a_units   = units.get("A",   "cm^3/(molecule*s)")
        t_units   = units.get("T",   "K")     # reactor Te units (K)
        de_units  = units.get("dE",  "eV")
        tmin_unit = units.get("Tmin","eV")
        tmax_unit = units.get("Tmax","eV")
        return a_units, t_units, de_units, tmin_unit, tmax_unit

    cpdef populate_from_yaml(self, object yaml_path_or_obj, int Z, int N,
                             bint allow_Z_gt28=False, Tmin=None, Tmax=None, comment=None):
        """
        Populate from voronov.yaml for the (Z,N) stage (N = electrons before ionization).
        Stores Tmin/Tmax internally in K (converted from eV in the YAML).
        """
        if not allow_Z_gt28 and Z > 28:
            raise ValueError(f"Voronov YAML: Z={Z} exceeds 28 (dataset covers H..Ni).")

        # Load YAML
        cdef dict data
        if isinstance(yaml_path_or_obj, dict):
            data = yaml_path_or_obj
        else:
            import yaml as _yaml
            with open(yaml_path_or_obj, "r") as f:
                data = _yaml.safe_load(f)

        row = self._vor__extract_row(data, Z, N)
        a_units, t_units, de_units, tmin_unit, tmax_unit = self._vor__units(data)

        # Required fields
        A  = float(row["A"])
        P  = float(row["P"])
        X  = float(row["X"])
        Kp = float(row["K"])
        dE = float(row["dE"])

        # Assign
        self.A = (A, a_units)
        self.P = P
        self.X = X
        self.K = Kp
        self.dE = dE  # eV

        # Temperature validity window: YAML gives eV; store in K
        cdef double kB_eVperK = 8.617333262e-5
        cdef double Tmin_eV = float(row.get("Tmin", 1.0))
        cdef double Tmax_eV = float(row.get("Tmax", 2.0e4))
        cdef object Tmin_K
        cdef object Tmax_K

        if Tmin is None:
            Tmin_K = (Tmin_eV / kB_eVperK, "K")
        else:
            Tmin_K = Tmin  # already something like (value, "K") from caller

        if Tmax is None:
            Tmax_K = (Tmax_eV / kB_eVperK, "K")
        else:
            Tmax_K = Tmax  # already something like (value, "K") from caller

        self.Tmin = Tmin_K
        self.Tmax = Tmax_K

        base = f"Voronov (1997) e-impact ionization fit, Z={Z}, N={N} "
        self.comment = f"{base}; {comment}" if comment else base
        return self

    @classmethod
    def from_yaml(cls, object yaml_path_or_obj, int Z, int N,
                  Tmin=None, Tmax=None, comment=None, bint allow_Z_gt28=False):
        """
        Construct a new VoronovEIArrhenius from voronov.yaml.
        """
        # Allocate without running the explicit-parameter branch of __init__. That
        # branch demands a dE, and satisfying it here would mean inventing an
        # ionization threshold -- the exact thing the branch exists to forbid --
        # for the table lookup below to overwrite. A placeholder that is always
        # overwritten is invisible until the day it isn't. So the object is built
        # empty and every parameter it ends up with comes from the table.
        obj = cls.__new__(cls)
        KineticsModel.__init__(obj, comment=comment or '')
        obj.uses_electron_temperature = True
        obj.uses_electron_density = True
        obj.populate_from_yaml(yaml_path_or_obj, Z, N,
                               Tmin=Tmin, Tmax=Tmax, comment=comment,
                               allow_Z_gt28=allow_Z_gt28)
        return obj

    def to_arrhenius(self, double Tmin=0.0, double Tmax=0.0):
        """
        Return an Arrhenius object that fits the Voronov kinetics over the specified temperature range.

        The fit is performed by evaluating the Voronov rate at intervals linear in 1/T
        (standard Arrhenius plot spacing) and fitting the modified Arrhenius equation to those points.

        If Tmin or Tmax are not provided (or 0.0), the method defaults to self.Tmin/self.Tmax.
        If those are also not set, it defaults to 10,000 K - 100,000 K (typical for ionization).
        """
        # Determine Temperature Boundaries
        if Tmin == 0.0:
            if self.Tmin is not None:
                Tmin = self.Tmin.value_si
            else:
                Tmin = 800.0

        if Tmax == 0.0:
            if self.Tmax is not None:
                Tmax = self.Tmax.value_si
            else:
                Tmax = 3000.0

        if Tmin >= Tmax:
            raise ValueError(f"Tmin ({Tmin}) must be less than Tmax ({Tmax}) for Arrhenius fitting.")

        # Generate sampling points
        # Linear in 1/T space gives better weighting for Arrhenius fits
        # We use 50 points to ensure a smooth fit
        cdef np.ndarray Tlist = 1.0 / np.linspace(1.0 / Tmax, 1.0 / Tmin, 50)
        cdef np.ndarray klist = np.zeros_like(Tlist)

        cdef int i
        cdef double T_val

        for i in range(len(Tlist)):
            T_val = Tlist[i]
            klist[i] = self.get_rate_coefficient(T_val)

        # Determine units
        # get_rate_coefficient returns SI Molar: m^3/(mol*s)
        cdef str kunits = "m^3/(mol*s)"

        # Create and fit Arrhenius object
        # T0=1.0 is standard for minimizing correlation between A and n
        arr = Arrhenius()
        arr.fit_to_data(Tlist, klist, kunits, T0=1.0)

        # Carry over metadata
        arr.Tmin = (Tmin, "K")
        arr.Tmax = (Tmax, "K")
        arr.comment = f"Fitted to VoronovEIArrhenius over range {Tmin}-{Tmax} K. Original comment: {self.comment}"

        return arr

    def to_chebyshev(self, double Tmin=0.0, double Tmax=0.0, int degree_t=10):
        """
        Convert the Voronov kinetics to a Chebyshev object.

        Uses a Pressure-Independent Chebyshev fit (P_basis=1).
        Default temperature range: 10,000 K - 100,000 K (typical for ionization).
        """
        from rmgpy.kinetics.chebyshev import Chebyshev

        # 1. Determine Temperature Range
        if Tmin == 0.0:
            Tmin = self.Tmin.value_si if self.Tmin is not None else 10000.0
        if Tmax == 0.0:
            Tmax = self.Tmax.value_si if self.Tmax is not None else 100000.0

        if Tmin >= Tmax:
            from rmgpy.exceptions import KineticsError
            raise KineticsError(
                f"VoronovEIArrhenius.to_chebyshev: Tmin ({Tmin}) must be < Tmax ({Tmax})."
            )

        # 2. Define Dummy Pressure Range (Required for Chebyshev format)
        cdef double Pmin = 100.0       # Pa  (≈ 0.001 bar)
        cdef double Pmax = 1.0e7       # Pa  (≈ 100 bar)

        # 3. Create the Chebyshev Object
        cheb = Chebyshev(
            Tmin=(Tmin, "K"), Tmax=(Tmax, "K"),
            Pmin=(Pmin, "Pa"), Pmax=(Pmax, "Pa"),
        )

        # 4. Generate Grid and Fit using Chebyshev roots
        #    IMPORTANT: need MORE grid points than polynomial degrees.
        cdef int nT = degree_t + 1      # > degree_t
        cdef int nP = 2                 # > degree_p (degree_p = 1)

        k_idx = np.arange(nT, dtype=float)

        invT_mid = 0.5 * (1.0 / Tmax + 1.0 / Tmin)
        invT_halfspan = 0.5 * (1.0 / Tmax - 1.0 / Tmin)

        # Chebyshev nodes in 1/T space
        T_nodes = invT_mid + invT_halfspan * np.cos(
            (2.0 * k_idx + 1.0) * np.pi / (2.0 * nT)
        )
        T_points = 1.0 / T_nodes

        # P_points: two points, but the rate is independent of P
        P_points = np.array([Pmin, Pmax], dtype=float)

        # K_data: shape (nT, nP)
        K_data = np.zeros((nT, nP), dtype=float)
        for i, T in enumerate(T_points):
            k_val = self.get_rate_coefficient(T)
            # Same k(T) for all pressures, since Voronov EI is P-independent
            K_data[i, :] = k_val

        # 5. Perform Fit
        # degree_t and degree_p = 1 (pressure-independent Chebyshev)
        cheb.fit_to_data(
            T_points,
            P_points,
            K_data,
            "m^3/(mol*s)",
            degree_t,
            1,          # degree_p
            Tmin,
            Tmax,
            Pmin,
            Pmax,
        )

        # Add original kinetics description to comment for traceability
        cheb.comment = (
            f"Chebyshev fit to VoronovEI (Z={getattr(self, 'Z', '?')}, "
            f"N={getattr(self, 'N', '?')}). Original: {self}"
        )
        return cheb

    cpdef TwoTemperaturePlasma to_two_temp_plasma(self):
        """
        Map this VoronovEIArrhenius (k = k(Te) only) onto a TwoTemperaturePlasma.

        The mapping is constructed so that:
          - along T = Te, k(T, Te) matches the Voronov rate (via Arrhenius fit),
          - for general T, the rate depends only on Te (T cancels analytically when Ea_g = 0).
        """
        # 1. Determine fitting window
        # Ionization is a high-energy process. Even if the object has Tmin=300K (default/room temp),
        # fitting Arrhenius from 300K includes a massive "dead zone" (rate ~ 0) that distorts
        # the curve and creates huge errors at the turn-on threshold.
        # We enforce a floor of 10,000 K (approx 1 eV) for the fit to ensure we capture the active regime.
        cdef double Tlo = 10000.0

        if self.Tmin is not None:
            # If user specified a Tmin, use it, but clamp it to at least 10,000 K to avoid artifacts.
            if self.Tmin.value_si > 10000.0:
                Tlo = self.Tmin.value_si

        cdef double Thi
        if self.Tmax is not None:
            Thi = self.Tmax.value_si
        else:
            Thi = 1.0e6  # Default to 1 MK if no max provided

        # 2. Fit Arrhenius to the electron temperature curve
        arr = self.to_arrhenius(Tmin=Tlo, Tmax=Thi)

        # Convert n to a bare float if needed
        try:
            n_val = float(arr.n.value_si)
        except AttributeError:
            n_val = float(arr.n)

        # 3. Build the canonical TwoTemperaturePlasma
        # Physics: EI is purely electron-driven. Ea_g is forced to 0.0.
        plasma = TwoTemperaturePlasma(
            A=(arr.A.value_si, "m^3/(mol*s)"),
            n=n_val,
            Ea_g=(0.0, "J/mol"),
            Ea_e=(arr.Ea.value_si, "J/mol"),
            Tmin=(Tlo, "K"),
            Tmax=(Thi, "K"),
            Pmin=self.Pmin,
            Pmax=self.Pmax,
            uncertainty=self.uncertainty,
            solute=self.solute,
            comment=(f"VoronovEIArrhenius mapped to TwoTemperaturePlasma over {Tlo:g}-{Thi:g} K (Te). Original: {self.comment}"),
        )

        return plasma

    def to_cantera_kinetics(self):
        """
        Export this VoronovEIArrhenius as a Cantera TwoTempPlasmaRate.

        Implementation:
          VoronovEIArrhenius → TwoTemperaturePlasma → ct.TwoTempPlasmaRate
        """
        plasma = self.to_two_temp_plasma()
        return plasma.to_cantera_kinetics()


################################################################################


cdef class ArrheniusEP(KineticsModel):
    """
    A kinetics model based on the (modified) Arrhenius equation, using the
    Evans-Polanyi equation to determine the activation energy. The attributes
    are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `A`             The preexponential factor
    `n`             The temperature exponent
    `alpha`         The Evans-Polanyi slope
    `E0`            The activation energy for a thermoneutral reaction
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `solute`        Transition state solute data
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================

    """

    def __init__(self, A=None, n=0.0, alpha=0.0, E0=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, uncertainty=None,
                 solute=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, uncertainty=uncertainty,
                               solute=solute, comment=comment)
        self.A = A
        self.n = n
        self.alpha = alpha
        self.E0 = E0

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        ArrheniusEP object.
        """
        string = 'ArrheniusEP(A={0!r}, n={1!r}, alpha={2!r}, E0={3!r}'.format(self.A, self.n, self.alpha, self.E0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.uncertainty is not None: string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.solute is not None: string += ', solute={0!r}'.format(self.solute)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an ArrheniusEP object.
        """
        return (ArrheniusEP, (self.A, self.n, self.alpha, self.E0, self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                              self.uncertainty, self.solute, self.comment))

    property A:
        """The preexponential factor."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.RateCoefficient(value)

    property n:
        """The temperature exponent."""
        def __get__(self):
            return self._n
        def __set__(self, value):
            self._n = quantity.Dimensionless(value)

    property alpha:
        """The Evans-Polanyi slope."""
        def __get__(self):
            return self._alpha
        def __set__(self, value):
            self._alpha = quantity.Dimensionless(value)

    property E0:
        """The activation energy for a thermoneutral reaction."""
        def __get__(self):
            return self._E0
        def __set__(self, value):
            self._E0 = quantity.Energy(value)

    cpdef double get_rate_coefficient(self, double T, double dHrxn=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3,
        mol, and s at temperature `T` in K and enthalpy of reaction `dHrxn`
        in J/mol.
        """
        cdef double A, n, Ea
        Ea = self.get_activation_energy(dHrxn)
        A = self._A.value_si
        n = self._n.value_si
        return A * T ** n * exp(-Ea / (constants.R * T))

    cpdef double get_activation_energy(self, double dHrxn) except -1:
        """
        Return the activation energy in J/mol corresponding to the given
        enthalpy of reaction `dHrxn` in J/mol.
        """
        cdef double Ea
        Ea = self._alpha.value_si * dHrxn + self._E0.value_si
        if self._E0.value_si > 0:
            if dHrxn < 0.0 and Ea < 0.0:
                Ea = 0.0
            elif dHrxn > 0.0 and Ea < dHrxn:
                Ea = dHrxn
        return Ea

    cpdef Arrhenius to_arrhenius(self, double dHrxn):
        """
        Return an :class:`Arrhenius` instance of the kinetics model using the
        given enthalpy of reaction `dHrxn` to determine the activation energy.
        """
        return Arrhenius(
            A=self.A,
            n=self.n,
            Ea=(self.get_activation_energy(dHrxn) * 0.001, "kJ/mol"),
            T0=(1, "K"),
            Tmin=self.Tmin,
            Tmax=self.Tmax,
            Pmin=self.Pmin,
            Pmax=self.Pmax,
            uncertainty=self.uncertainty,
            solute=self.solute,
            comment=self.comment,
        )

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Must match temperature
        and pressure range of kinetics model, as well as parameters: A, n, Ea, T0. (Shouldn't have pressure
        range if it's Arrhenius.) Otherwise returns ``False``.
        """
        if not isinstance(other_kinetics, ArrheniusEP):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if (not self.A.equals(other_kinetics.A) or not self.n.equals(other_kinetics.n)
                or not self.alpha.equals(other_kinetics.alpha) or not self.E0.equals(other_kinetics.E0)):
            return False

        return True

    cpdef change_rate(self, double factor):
        """
        Changes A factor by multiplying it by a ``factor``.
        """
        self._A.value_si *= factor

    def set_cantera_kinetics(self, ct_reaction, species_list):
        """
        Sets a cantera Reaction() object with the modified Arrhenius object
        converted to an Arrhenius form.
        """
        raise NotImplementedError('set_cantera_kinetics() is not implemented for ArrheniusEP class kinetics.')

################################################################################

cdef class ArrheniusBM(KineticsModel):
    """
    A kinetics model based on the (modified) Arrhenius equation, using the
    Blowers-Masel equation to determine the activation energy.
    Based on Blowers and Masel's 2000 paper Engineering Approximations for Activation
    Energies in Hydrogen Transfer Reactions.
    The attributes are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `A`             The preexponential factor
    `n`             The temperature exponent
    `w0`            The average of the bond dissociation energies of the bond formed and the bond broken
    `E0`            The activation energy for a thermoneutral reaction
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `solute`        Transition state solute data
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================

    """

    def __init__(self, A=None, n=0.0, w0=(0.0, 'J/mol'), E0=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None,
                 uncertainty=None, solute=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, uncertainty=uncertainty,
                               solute=solute, comment=comment)
        self.A = A
        self.n = n
        self.w0 = w0
        self.E0 = E0

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        ArrheniusBM object.
        """
        string = 'ArrheniusBM(A={0!r}, n={1!r}, w0={2!r}, E0={3!r}'.format(self.A, self.n, self.w0, self.E0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.uncertainty is not None: string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.solute is not None: string += ', solute={0!r}'.format(self.solute)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an ArrheniusEP object.
        """
        return (ArrheniusBM, (self.A, self.n, self.w0, self.E0, self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                              self.uncertainty, self.solute, self.comment))

    property A:
        """The preexponential factor."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.RateCoefficient(value)

    property n:
        """The temperature exponent."""
        def __get__(self):
            return self._n
        def __set__(self, value):
            self._n = quantity.Dimensionless(value)

    property w0:
        """The average of the bond dissociation energies of the bond formed and the bond broken."""
        def __get__(self):
            return self._w0
        def __set__(self, value):
            self._w0 = quantity.Energy(value)

    property E0:
        """The activation energy for a thermoneutral reaction."""
        def __get__(self):
            return self._E0
        def __set__(self, value):
            self._E0 = quantity.Energy(value)

    cpdef double get_rate_coefficient(self, double T, double dHrxn=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3,
        mol, and s at temperature `T` in K and enthalpy of reaction `dHrxn`
        in J/mol, evaluated at 298 K.
        """
        cdef double A, n, Ea
        Ea = self.get_activation_energy(dHrxn)
        A = self._A.value_si
        n = self._n.value_si
        return A * T ** n * exp(-Ea / (constants.R * T))

    cpdef double get_activation_energy(self, double dHrxn) except -1:
        """
        Return the activation energy in J/mol corresponding to the given
        enthalpy of reaction `dHrxn` in J/mol, evaluated at 298 K.
        """
        cdef double w0, E0, Ea
        E0 = self._E0.value_si
        w0 = self._w0.value_si
        if E0 < 0:
            # Negative E0 is unphysical, but could be encountered during optimization.
            Ea = E0 + max(dHrxn, 0.0)
        else:
            Vp = 2 * w0 * (w0 + E0) / (w0 - E0)
            Ea = (w0 + dHrxn / 2.0) * (Vp - 2 * w0 + dHrxn) * (Vp - 2 * w0 + dHrxn) / (Vp * Vp - (4 * w0 * w0) + dHrxn * dHrxn)
            if (Ea < 0.0) or (dHrxn < -4 * E0):
                Ea = 0.0
            elif (Ea < dHrxn) or (dHrxn > 4 * E0):
                Ea = dHrxn
        return Ea

    cpdef Arrhenius to_arrhenius(self, double dHrxn):
        """
        Return an :class:`Arrhenius` instance of the kinetics model using the
        given enthalpy of reaction `dHrxn` (in J/mol, evaluated at 298 K)
        to determine the activation energy.
        """
        return Arrhenius(
            A=self.A,
            n=self.n,
            Ea=(self.get_activation_energy(dHrxn) * 0.001, "kJ/mol"),
            T0=(1, "K"),
            Tmin=self.Tmin,
            Tmax=self.Tmax,
            uncertainty=self.uncertainty,
            solute=self.solute,
            comment=self.comment,
        )

    def fit_to_reactions(self, rxns, w0=None, recipe=None, Ts=None):
        """
        Fit an ArrheniusBM model to a list of reactions at the given temperatures,
        w0 must be either given or estimated using the family object

        WARNING: there's a lot of code duplication with ArrheniusChargeTransferBM.fit_to_reactions
                 so anything you change here you should probably change there too and vice versa!
        """
        assert w0 is not None or recipe is not None, 'either w0 or recipe must be specified'

        if Ts is None:
            Ts = np.array([300.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1500.0, 2000.0])
        if w0 is None:
            #estimate w0
            w0s = get_w0s(recipe, rxns)
            w0 = sum(w0s) / len(w0s)
        self.w0 = (w0 * 0.001, 'kJ/mol')

        if len(rxns) == 1:
            rxn = rxns[0]
            dHrxn = rxn.get_enthalpy_of_reaction(298.0)
            A = rxn.kinetics.A.value_si
            n = rxn.kinetics.n.value_si
            Ea = rxn.kinetics.Ea.value_si

            def kfcn(E0):
                Vp = 2 * w0 * (w0 + E0) / (w0 - E0)
                out = Ea - (w0 + dHrxn / 2.0) * (Vp - 2 * w0 + dHrxn) * (Vp - 2 * w0 + dHrxn) / (Vp * Vp - (4 * w0 * w0) + dHrxn * dHrxn)
                return out

            E0 = fsolve(kfcn, w0 / 10.0)[0]

            self.Tmin = rxn.kinetics.Tmin
            self.Tmax = rxn.kinetics.Tmax
            self.solute = None
            self.comment = 'Fitted to 1 reaction.'
        else:
            # define optimization function
            def kfcn(xs, lnA, n, E0):
                T = xs[:,0]
                dHrxn = xs[:,1]
                Vp = 2 * w0 * (w0 + E0) / (w0 - E0)
                Ea = (w0 + dHrxn / 2.0) * (Vp - 2 * w0 + dHrxn) * (Vp - 2 * w0 + dHrxn) / (Vp * Vp - (4 * w0 * w0) + dHrxn * dHrxn)
                Ea = np.where(dHrxn< -4.0*E0, 0.0, Ea)
                Ea = np.where(dHrxn > 4.0*E0, dHrxn, Ea)
                return lnA + np.log(T) * n + (-Ea / (constants.R * T))

            # get (T,dHrxn(T)) -> (Ln(k) mappings
            xdata = []
            ydata = []
            sigmas = []
            E0 = 0.0
            lnA = 0.0
            n = 0.0
            for rxn in rxns:
                # approximately correct the overall uncertainties to std deviations
                s = rank_accuracy_map[rxn.rank].value_si/2.0
                dHrxn = rxn.get_enthalpy_of_reaction(298.0)
                for T in Ts:
                    xdata.append([T, dHrxn])
                    ydata.append(np.log(rxn.get_rate_coefficient(T)))
                    sigmas.append(s / (constants.R * T))
                # Use BEP with alpha = 0.25 for inital guess of E0
                E0 += rxn.kinetics._Ea.value_si - 0.25 * dHrxn
                lnA += np.log(rxn.kinetics.A.value_si)
                n += rxn.kinetics.n.value_si
            # Use the averages as intial guess
            E0 /= len(rxns)
            lnA /= len(rxns)
            n /= len(rxns)
            E0 = min(E0, w0)
            w0 = max(2 * E0, w0) # Expression only works if w0>2E0, and is insensitive to w0
            self.w0 = (w0 * 0.001, 'kJ/mol')
            if E0 < 0:
                E0 = w0 / 100.0

            xdata = np.array(xdata)
            ydata = np.array(ydata)

            # fit parameters
            keep_trying = True
            xtol = 1e-8
            ftol = 1e-8
            attempts = 0
            while keep_trying:
                keep_trying = False
                try:
                    params = curve_fit(kfcn, xdata, ydata, sigma=sigmas, p0=[lnA, n, E0], xtol=xtol, ftol=ftol)
                    lnA, n, E0 = params[0].tolist()
                    if abs(E0/self.w0.value_si) > 1 and attempts < 5:
                        keep_trying = True
                        if attempts > 0:
                            self.w0.value_si *= 1.25
                            w0 = self.w0.value_si # update local variable for use in kfcn optimization function
                        attempts += 1
                        E0 = self.w0.value_si / 10.0
                except RuntimeError:
                    if xtol < 1.0:
                        keep_trying = True
                        xtol *= 10.0
                        ftol *= 10.0
                    else:
                        raise ValueError("Could not fit BM arrhenius to reactions with xtol<1.0")

            A = np.exp(lnA)

            self.Tmin = (np.min(Ts), "K")
            self.Tmax = (np.max(Ts), "K")
            self.solute = None
            self.comment = 'Fitted to {0} reactions at temperatures: {1}'.format(len(rxns), Ts)

        # fill in parameters
        A_units = ['', 's^-1', 'm^3/(mol*s)', 'm^6/(mol^2*s)']
        order = len(rxns[0].reactants)
        if order != 1 and rxn.is_surface_reaction():
            raise NotImplementedError("Units not implemented for surface reactions.")
        self.A = (A, A_units[order])

        self.n = n
        self.E0 = (E0 * 0.001, 'kJ/mol')

        return self

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Must match temperature
        and pressure range of kinetics model, as well as parameters: A, n, Ea, T0. (Shouldn't have pressure
        range if it's Arrhenius.) Otherwise returns ``False``.
        """
        if not isinstance(other_kinetics, ArrheniusBM):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if (not self.A.equals(other_kinetics.A) or not self.n.equals(other_kinetics.n)
                or not self.w0.equals(other_kinetics.w0) or not self.E0.equals(other_kinetics.E0)):
            return False

        return True

    cpdef change_rate(self, double factor):
        """
        Changes A factor by multiplying it by a ``factor``.
        """
        self._A.value_si *= factor

    def to_cantera_kinetics(self):
        """
        Converts the RMG ArrheniusBM object to a cantera BlowersMaselRate. 

        BlowersMaselRate(A, b, Ea, W)  where A is in units of m^3/kmol/s, 
        b is dimensionless, and Ea and W are in J/kmol
        """
        import cantera as ct

        rate_units_conversion = {'1/s': 1,
                                 's^-1': 1,
                                 'm^3/(mol*s)': 1000,
                                 'm^6/(mol^2*s)': 1000000,
                                 'cm^3/(mol*s)': 1000,
                                 'cm^6/(mol^2*s)': 1000000,
                                 'm^3/(molecule*s)': 1000,
                                 'm^6/(molecule^2*s)': 1000000,
                                 'cm^3/(molecule*s)': 1000,
                                 'cm^6/(molecule^2*s)': 1000000,
                                 }

        A = self._A.value_si

        try:
            A *= rate_units_conversion[self._A.units] # convert from /mol to /kmol
        except KeyError:
            raise ValueError(f'ArrheniusBM A-factor units {self._A.units} not found among accepted '
                             'units for converting to Cantera BlowersMaselRate object.')

        b = self._n.value_si
        Ea = self._E0.value_si * 1000  # convert from J/mol to J/kmol
        w = self._w0.value_si * 1000  # convert from J/mol to J/kmol

        rate = ct.BlowersMaselRate(A, b, Ea, w)
        if A < 0:
            rate.allow_negative_pre_exponential_factor = True
        return rate

    def set_cantera_kinetics(self, ct_reaction, species_list):
        """
        Accepts a cantera Reaction object and sets its rate to a Cantera BlowersMaselRate object.
        """
        import cantera as ct
        if not isinstance(ct_reaction.rate, ct.BlowersMaselRate):
            raise TypeError("ct_reaction must have a cantera BlowersMaselRate as the rate attribute")
        ct_reaction.rate = self.to_cantera_kinetics()

################################################################################

cdef class PDepArrhenius(PDepKineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient :math:`k(T,P)` where
    a set of Arrhenius kinetics are stored at a variety of pressures and
    interpolated between on a logarithmic scale. The attributes are:

    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `pressures`     The list of pressures
    `arrhenius`     The list of :class:`Arrhenius` objects at each pressure
    `Tmin`          The minimum temperature in K at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature in K at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure in bar at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure in bar at which the model is valid, or zero if unknown or undefined
    `efficiencies`  A dict associating chemical species with associated efficiencies
    `order`         The reaction order (1 = first, 2 = second, etc.)
    `comment`       Information about the model (e.g. its source)
    =============== ============================================================

    """

    def __init__(self, pressures=None, arrhenius=None, highPlimit=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None,
                 comment=''):
        PDepKineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, highPlimit=highPlimit,
                                   comment=comment)
        self.pressures = pressures
        self.arrhenius = arrhenius or []

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        PDepArrhenius object.
        """
        string = 'PDepArrhenius(pressures={0!r}, arrhenius={1!r}'.format(self.pressures, self.arrhenius)
        if self.highPlimit is not None: string += ', highPlimit={0!r}'.format(self.highPlimit)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a PDepArrhenius object.
        """
        return (PDepArrhenius, (self.pressures, self.arrhenius, self.highPlimit, self.Tmin, self.Tmax,
                                self.Pmin, self.Pmax, self.comment))

    property pressures:
        """The list of pressures."""
        def __get__(self):
            return self._pressures
        def __set__(self, value):
            self._pressures = quantity.Pressure(value)

    cdef get_adjacent_expressions(self, double P):
        """
        Returns the pressures and Arrhenius expressions for the pressures that
        most closely bound the specified pressure `P` in Pa.
        """
        cdef np.ndarray[np.float64_t, ndim=1] pressures
        cdef int i, ilow, ihigh

        pressures = self._pressures.value_si

        ilow = 0
        ihigh = -1
        for i in range(pressures.shape[0]):
            if pressures[i] <= P:
                ilow = i
            if pressures[i] >= P and ihigh == -1:
                ihigh = i
        return pressures[ilow], pressures[ihigh], self.arrhenius[ilow], self.arrhenius[ihigh]

    cpdef double get_rate_coefficient(self, double T, double P=0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3,
        mol, and s at temperature `T` in K and pressure `P` in Pa.
        """
        cdef double Plow, Phigh, klow, khigh, k
        cdef KineticsModel alow, ahigh
        cdef int j

        if P == 0:
            raise ValueError('No pressure specified to pressure-dependent PDepArrhenius.get_rate_coefficient().')

        k = 0.0
        Plow, Phigh, alow, ahigh = self.get_adjacent_expressions(P)
        if Plow == Phigh:
            k = alow.get_rate_coefficient(T)
        else:
            klow = alow.get_rate_coefficient(T)
            khigh = ahigh.get_rate_coefficient(T)
            if klow == khigh == 0.0: return 0.0
            k = klow * 10 ** (log10(P / Plow) / log10(Phigh / Plow) * log10(khigh / klow))
        return k

    cpdef fit_to_data(self, np.ndarray Tlist, np.ndarray Plist, np.ndarray K, str kunits, double T0=1):
        """
        Fit the pressure-dependent Arrhenius model to a matrix of rate
        coefficient data `K` with units of `kunits` corresponding to a set of
        temperatures `Tlist` in K and pressures `Plist` in Pa. An Arrhenius
        model is fit cpdef change_rate(self, double factor)at each pressure.
        """
        cdef int i
        self.pressures = (Plist * 1e-5, "bar")
        self.arrhenius = []
        for i in range(len(Plist)):
            arrhenius = Arrhenius().fit_to_data(Tlist, K[:, i], kunits, T0)
            self.arrhenius.append(arrhenius)
        return self

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Each duplicate
        reaction must be matched and equal to that in the other PDepArrhenius model
        in the same order.  Otherwise returns ``False``
        """
        if not isinstance(other_kinetics, PDepArrhenius):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if len(self.arrhenius) != len(other_kinetics.arrhenius):
            return False
        if not self.pressures.equals(other_kinetics.pressures):
            return False
        for index in range(len(self.arrhenius)):
            if not self.arrhenius[index].is_identical_to(other_kinetics.arrhenius[index]):
                return False
        if self.highPlimit and not self.highPlimit.equals(other_kinetics.highPlimit):
            return False

        return True

    cpdef change_rate(self, double factor):
        """
        Changes kinetics rate by a multiple ``factor``.
        """
        for kin in self.arrhenius:
            kin.change_rate(factor)
        if self.highPlimit is not None:
            self.highPlimit.change_rate(factor)

    def set_cantera_kinetics(self, ct_reaction, species_list):
        """
        Sets a Cantera PlogReaction()'s `rates` attribute with
        a list of tuples containing [(pressure in Pa, cantera arrhenius object), ...].

        A ``MultiArrhenius`` entry (from chemkin PLOG blocks with duplicate
        pressures) is expanded into one tuple per inner Arrhenius; Cantera's
        PlogRate sums duplicate-pressure entries at evaluation.
        """
        import cantera as ct
        assert isinstance(ct_reaction.rate, ct.PlogRate), "Must have a Cantera PlogRate attribute"

        rate_pairs = []
        for P, arr in zip(self._pressures.value_si, self.arrhenius):
            if isinstance(arr, MultiArrhenius):
                for sub in arr.arrhenius:
                    rate_pairs.append((P, sub.to_cantera_kinetics(arrhenius_class=True)))
            else:
                rate_pairs.append((P, arr.to_cantera_kinetics(arrhenius_class=True)))

        ct_reaction.rate = ct.PlogRate(rate_pairs)

################################################################################

cdef class MultiArrhenius(KineticsModel):
    """
    A kinetics model based on a set of (modified) Arrhenius equations, which
    are summed to obtain the overall rate. The attributes are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `arrhenius`     A list of the :class:`Arrhenius` kinetics
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================

    """

    def __init__(self, arrhenius=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.arrhenius = arrhenius

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        MultiArrhenius object.
        """
        string = 'MultiArrhenius(arrhenius={0!r}'.format(self.arrhenius)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an MultiArrhenius object.
        """
        return (MultiArrhenius, (self.arrhenius, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    cpdef double get_rate_coefficient(self, double T, double P=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3,
        mol, and s at temperature `T` in K.
        """
        cdef double k
        cdef Arrhenius arrh
        k = 0.0
        for arrh in self.arrhenius:
            k += arrh.get_rate_coefficient(T)
        return k

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Each duplicate
        reaction must be matched and equal to that in the other MultiArrhenius model
        in the same order.  Otherwise returns ``False``
        """
        if not isinstance(other_kinetics, MultiArrhenius):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if len(self.arrhenius) != len(other_kinetics.arrhenius):
            return False

        for index in range(len(self.arrhenius)):
            if not self.arrhenius[index].is_identical_to(other_kinetics.arrhenius[index]):
                return False

        return True

    cpdef Arrhenius to_arrhenius(self, double Tmin=-1, double Tmax=-1):
        """
        Return an :class:`Arrhenius` instance of the kinetics model

        Fit the Arrhenius parameters to a set of rate coefficient data generated
        from the MultiArrhenius kinetics, over the temperature range
        Tmin to Tmax, in Kelvin. If Tmin or Tmax are unspecified (or -1)
        then the MultiArrhenius's Tmin and Tmax are used.
        A linear least-squares fit is used, which guarantees that the
        resulting parameters provide the best possible approximation to the
        data.
        """
        cdef Arrhenius arrh
        cdef np.ndarray Tlist, klist
        cdef str kunits
        if Tmin == -1: Tmin = self.Tmin.value_si
        if Tmax == -1: Tmax = self.Tmax.value_si
        kunits = str(quantity.pq.Quantity(1.0, self.arrhenius[0].A.units).simplified).split()[-1]  # is this the best way to get the units returned by k??
        Tlist = np.logspace(log10(Tmin), log10(Tmax), num=25)
        klist = np.array(list(map(self.get_rate_coefficient, Tlist)), float)
        arrh = Arrhenius().fit_to_data(Tlist, klist, kunits)
        arrh.comment = "Fitted to Multiple Arrhenius kinetics over range {Tmin}-{Tmax} K. {comment}".format(
            Tmin=Tmin, Tmax=Tmax, comment=self.comment)
        return arrh

    cpdef change_rate(self, double factor):
        """
        Change kinetics rate by a multiple ``factor``.
        """
        for kin in self.arrhenius:
            kin.change_rate(factor)

    def set_cantera_kinetics(self, ct_reaction, species_list):
        """
        Sets the kinetic rates for a list of cantera `Reaction` objects
        Here, ct_reaction must be a list rather than a single cantera reaction.
        """
        if len(ct_reaction) != len(self.arrhenius):
            raise Exception('The number of Cantera Reaction objects does not match the number of Arrhenius objects')

        for i, arr in enumerate(self.arrhenius):
            arr.set_cantera_kinetics(ct_reaction[i], species_list)

################################################################################

cdef class MultiPDepArrhenius(PDepKineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient :math:`k(T,P)` where
    sets of Arrhenius kinetics are stored at a variety of pressures and
    interpolated between on a logarithmic scale. The attributes are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `arrhenius`     A list of the :class:`PDepArrhenius` kinetics at each temperature
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================

    """

    def __init__(self, arrhenius=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        PDepKineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.arrhenius = arrhenius

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        MultiPDepArrhenius object.
        """
        string = 'MultiPDepArrhenius(arrhenius={0!r}'.format(self.arrhenius)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an MultiPDepArrhenius object.
        """
        return (MultiPDepArrhenius, (self.arrhenius, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    cpdef double get_rate_coefficient(self, double T, double P=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3,
        mol, and s at temperature `T` in K and pressure `P` in Pa.
        """
        cdef double k, klow, khigh, Plow, Phigh
        cdef PDepArrhenius arrh
        cdef KineticsModel arrh_low, arrh_high
        cdef np.ndarray Plist1, Plist2
        cdef int i

        if P == 0:
            raise ValueError('No pressure specified to pressure-dependent MultiPDepArrhenius.get_rate_coefficient().')

        k = 0
        for arrh in self.arrhenius:
            k += arrh.get_rate_coefficient(T, P)

        return k

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Each duplicate
        reaction must be matched and equal to that in the other MultiArrhenius model
        in the same order.  Otherwise returns ``False``
        """
        if not isinstance(other_kinetics, MultiPDepArrhenius):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if len(self.arrhenius) != len(other_kinetics.arrhenius):
            return False

        for index in range(len(self.arrhenius)):
            if not self.arrhenius[index].is_identical_to(other_kinetics.arrhenius[index]):
                return False

        return True

    cpdef change_rate(self, double factor):
        """
        Change kinetic rate by a multiple ``factor``.
        """
        for kin in self.arrhenius:
            kin.change_rate(factor)

    def set_cantera_kinetics(self, ct_reaction, species_list):
        """
        Sets the PLOG kinetics for multiple cantera `Reaction` objects, provided in a list.
        ct_reaction is a list of cantera reaction objects.
        """
        if len(ct_reaction) != len(self.arrhenius):
            raise Exception('The number of Cantera Reaction objects does not match the number of PdepArrhenius objects')

        for i, arr in enumerate(self.arrhenius):
            arr.set_cantera_kinetics(ct_reaction[i], species_list)

################################################################################

cdef class ArrheniusChargeTransfer(KineticsModel):

    """
    A kinetics model for charge transfer reactions in the gas phase.

    It is very similar to the :class:`Arrhenius`, but the Ea is potential-dependent. `A` carries
    volumetric units (e.g. ``m^3/(mol*s)``), as used by the gas-phase cation families. The
    potential dependence is inherited from the surface charge transfer models, and it is only
    active when a reaction is evaluated away from `V0`: in the gas phase there is no electrode, so
    the reaction is evaluated at `V0` and the ``alpha * electrons * F * (V - V0)`` term contributes
    no shift. `electrons` and `alpha` are still stored and still meaningful — `electrons` records
    the electron stoichiometry of the reaction, and `alpha` remains the symmetry factor that would
    apply if an external potential were deliberately imposed.

    The attributes are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `A`             The preexponential factor
    `T0`            The reference temperature
    `n`             The temperature exponent
    `Ea`            The activation energy at the reference potential `V0`
    `electrons`     The stochiometry coeff for electrons (negative if reactant, positive if product)
    `V0`            The reference potential; the rate is evaluated here unless a potential is imposed
    `alpha`         The charge transfer (symmetry) coefficient
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `solute`        The transition state solute data
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================

    Note that the second positional argument of :meth:`get_rate_coefficient` is a potential in V,
    not a pressure in Pa. There is no pressure dependence in this model.

    """

    def __init__(self, A=None, n=0.0, Ea=None, V0=None, alpha=0.5, electrons=-1, T0=(1.0, "K"), Tmin=None, Tmax=None,
                Pmin=None, Pmax=None, solute=None, uncertainty=None, comment=''):

        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, solute=solute, uncertainty=uncertainty,
                 comment=comment)

        self.alpha = alpha
        self.A = A
        self.n = n
        self.Ea = Ea
        self.T0 = T0
        self.electrons = electrons
        self.V0 = V0

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Arrhenius object.
        """
        string = 'ArrheniusChargeTransfer(A={0!r}, n={1!r}, Ea={2!r}, V0={3!r}, alpha={4!r}, electrons={5!r}, T0={6!r}'.format(
            self.A, self.n, self.Ea, self.V0, self.alpha, self.electrons, self.T0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.solute: string += ', solute={0!r}'.format(self.solute)
        if self.uncertainty: string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a ArrheniusChargeTransfer object.
        """
        return (ArrheniusChargeTransfer, (self.A, self.n, self.Ea, self.V0, self.alpha, self.electrons, self.T0, self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                            self.solute, self.uncertainty, self.comment))

    property A:
        """The preexponential factor."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.SurfaceRateCoefficient(value)

    property n:
        """The temperature exponent."""
        def __get__(self):
            return self._n
        def __set__(self, value):
            self._n = quantity.Dimensionless(value)

    property Ea:
        """The activation energy."""
        def __get__(self):
            return self._Ea
        def __set__(self, value):
            self._Ea = quantity.Energy(value)

    property T0:
        """The reference temperature."""
        def __get__(self):
            return self._T0
        def __set__(self, value):
            self._T0 = quantity.Temperature(value)

    property V0:
        """The reference potential."""
        def __get__(self):
            return self._V0
        def __set__(self, value):
            self._V0 = quantity.Potential(value)

    property electrons:
        """The number of electrons transferred."""
        def __get__(self):
            return self._electrons
        def __set__(self, value):
            self._electrons = quantity.Dimensionless(value)

    property alpha:
        """The charge transfer coefficient."""
        def __get__(self):
            return self._alpha
        def __set__(self, value):
            self._alpha = quantity.Dimensionless(value)

    cpdef double get_activation_energy_from_potential(self, double V=0.0, bint non_negative=True):
        """
        Return the effective activation energy (in J/mol) at specificed potential (in Volts).
        """
        cdef double  electrons, alpha, Ea, V0

        electrons = self._electrons.value_si
        alpha = self._alpha.value_si
        Ea = self._Ea.value_si
        V0 = self._V0.value_si

        Ea -= alpha * electrons * constants.F * (V-V0)

        if non_negative is True:
            if Ea < 0:
                Ea = 0.0

        return Ea

    cpdef double get_rate_coefficient(self, double T, double V=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3,
        mol, and s at temperature `T` in K and potential `V` in volts.

        `V` is a potential, not a pressure. For gas-phase chemistry it should be left at the
        reference potential `V0`, which is what :meth:`rmgpy.reaction.Reaction.get_rate_coefficient`
        supplies by default.
        """
        cdef double A, n, V0, T0, Ea

        A = self._A.value_si
        n = self._n.value_si
        V0 = self._V0.value_si
        T0 = self._T0.value_si

        if V != V0:
            Ea = self.get_activation_energy_from_potential(V)
        else:
            Ea = self._Ea.value_si

        return A * (T / T0) ** n * exp(-Ea / (constants.R * T))

    cpdef change_t0(self, double T0):
        """
        Changes the reference temperature used in the exponent to `T0` in K,
        and adjusts the preexponential factor accordingly.
        """
        self._A.value_si /= (self._T0.value_si / T0) ** self._n.value_si
        self._T0.value_si = T0

    cpdef change_v0(self, double V0):
        """
        Changes the reference potential to `V0` in volts, and adjusts the
        activation energy `Ea` accordingly.
        """

        self._Ea.value_si = self.get_activation_energy_from_potential(V0)
        self._V0.value_si = V0

    cpdef fit_to_data(self, np.ndarray Tlist, np.ndarray klist, str kunits, double T0=1,
                      np.ndarray weights=None, bint three_params=False):
        """
        Fit the Arrhenius parameters to a set of rate coefficient data `klist`
        in units of `kunits` corresponding to a set of temperatures `Tlist` in
        K. A linear least-squares fit is used, which guarantees that the
        resulting parameters provide the best possible approximation to the
        data.
        """
        import scipy.stats
        if not all(np.isfinite(klist)):
            raise  ValueError("Rates must all be finite, not inf or NaN")
        if any(klist<0):
            if not all(klist<0):
                raise ValueError("Rates must all be positive or all be negative.")
            rate_sign_multiplier = -1
            klist = -1 * klist
        else:
            rate_sign_multiplier = 1

        assert len(Tlist) == len(klist), "length of temperatures and rates must be the same"
        if len(Tlist) < 3 + three_params:
            raise KineticsError('Not enough degrees of freedom to fit this Arrhenius expression')
        if three_params:
            A = np.zeros((len(Tlist), 3), np.float64)
            A[:, 0] = np.ones_like(Tlist)
            A[:, 1] = np.log(Tlist / T0)
            A[:, 2] = -1.0 / constants.R / Tlist
        else:
            A = np.zeros((len(Tlist), 2), np.float64)
            A[:, 0] = np.ones_like(Tlist)
            A[:, 1] = -1.0 / constants.R / Tlist
        b = np.log(klist)
        if weights is not None:
            for n in range(b.size):
                A[n, :] *= weights[n]
                b[n] *= weights[n]
        x, residues, rank, s = np.linalg.lstsq(A, b, rcond=RCOND)

        # Determine covarianace matrix to obtain parameter uncertainties
        count = klist.size
        cov = residues[0] / (count - 3) * np.linalg.inv(np.dot(A.T, A))
        t = scipy.stats.t.ppf(0.975, count - 3)

        if not three_params:
            x = np.array([x[0], 0, x[1]])
            cov = np.array([[cov[0, 0], 0, cov[0, 1]], [0, 0, 0], [cov[1, 0], 0, cov[1, 1]]])

        self.A = (rate_sign_multiplier * exp(x[0]), kunits)
        self.n = x[1]
        self.Ea = (x[2] * 0.001, "kJ/mol")
        self.T0 = (T0, "K")
        self.Tmin = (np.min(Tlist), "K")
        self.Tmax = (np.max(Tlist), "K")
        self.solute = None,
        self.comment = 'Fitted to {0:d} data points; dA = *|/ {1:g}, dn = +|- {2:g}, dEa = +|- {3:g} kJ/mol'.format(
            len(Tlist),
            exp(sqrt(cov[0, 0])),
            sqrt(cov[1, 1]),
            sqrt(cov[2, 2]) * 0.001,
        )

        return self

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Must match temperature
        and pressure range of kinetics model, as well as parameters: A, n, Ea, T0. (Shouldn't have pressure
        range if it's Arrhenius.) Otherwise returns ``False``.
        """
        if not isinstance(other_kinetics, ArrheniusChargeTransfer):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if (not self.A.equals(other_kinetics.A) or not self.n.equals(other_kinetics.n)
                or not self.Ea.equals(other_kinetics.Ea) or not self.T0.equals(other_kinetics.T0)
                or not self.alpha.equals(other_kinetics.alpha) or not self.electrons.equals(other_kinetics.electrons)
                or not self.V0.equals(other_kinetics.V0)):
            return False

        return True

    cpdef change_rate(self, double factor):
        """
        Changes A factor in Arrhenius expression by multiplying it by a ``factor``.
        """
        self._A.value_si *= factor

cdef class ArrheniusChargeTransferBM(KineticsModel):
    """
    A kinetics model for charge transfer reactions in the gas phase, based on the (modified)
    Arrhenius equation and using the Evans-Polanyi equation to determine the activation energy.

    This is the Blowers-Masel form of :class:`ArrheniusChargeTransfer`, used for the rate rules of
    the gas-phase cation families; `A` carries volumetric units (e.g. ``m^3/(mol*s)``). As for
    :class:`ArrheniusChargeTransfer`, there is no electrode in the gas phase, so the reaction is
    evaluated at the reference potential `V0` and the potential term contributes no shift, while
    `electrons` and `alpha` are retained on the record. A reaction generated from a family carries
    this model only until :meth:`rmgpy.reaction.Reaction.fix_barrier_height` converts it to an
    :class:`ArrheniusChargeTransfer` via :meth:`to_arrhenius_charge_transfer`; reactions from a
    kinetics library, from an imported mechanism, or built by hand are not put through that
    conversion and can reach a solver still carrying this model.

    Note that :meth:`to_arrhenius_charge_transfer` does not carry `alpha` across to the
    :class:`ArrheniusChargeTransfer` it builds, so a converted record silently takes alpha's 0.5
    default. That is invisible at `V0`, where the potential term vanishes, but not under an applied
    field.

    The attributes are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `A`             The preexponential factor
    `n`             The temperature exponent
    `w0`            The average of the bond dissociation energies of the bond formed and the bond broken
    `E0`            The activation energy for a thermoneutral reaction, at the reference potential `V0`
    `electrons`     The stochiometry coeff for electrons (negative if reactant, positive if product)
    `V0`            The reference potential; the rate is evaluated here
    `alpha`         The charge transfer (symmetry) coefficient
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `solute`        Transition state solute data
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================

    Note that the second positional argument of :meth:`get_rate_coefficient` is an enthalpy of
    reaction in J/mol -- neither a pressure nor a potential. There is no pressure dependence in
    this model.

    """

    def __init__(self, A=None, n=0.0, w0=(0.0, 'J/mol'), E0=None, V0=(0.0,'V'), alpha=0.5, electrons=-1, Tmin=None, Tmax=None,
                Pmin=None, Pmax=None, solute=None, uncertainty=None, comment=''):

        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, solute=solute, uncertainty=uncertainty,
                 comment=comment)

        self.alpha = alpha
        self.A = A
        self.n = n
        self.w0 = w0
        self.E0 = E0
        self.electrons =  electrons
        self.V0 = V0

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Arrhenius object.
        """
        string = 'ArrheniusChargeTransferBM(A={0!r}, n={1!r}, w0={2!r}, E0={3!r}, V0={4!r}, alpha={5!r}, electrons={6!r}'.format(
            self.A, self.n, self.w0, self.E0, self.V0, self.alpha, self.electrons)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.solute: string += ', solute={0!r}'.format(self.solute)
        if self.uncertainty: string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a ArrheniusChargeTransfer object.
        """
        return (ArrheniusChargeTransferBM, (self.A, self.n, self.w0, self.E0, self.V0, self.alpha, self.electrons, self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                            self.solute, self.uncertainty, self.comment))

    property A:
        """The preexponential factor."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.SurfaceRateCoefficient(value)

    property n:
        """The temperature exponent."""
        def __get__(self):
            return self._n
        def __set__(self, value):
            self._n = quantity.Dimensionless(value)

    property w0:
        """The average of the bond dissociation energies of the bond formed and the bond broken."""
        def __get__(self):
            return self._w0
        def __set__(self, value):
            self._w0 = quantity.Energy(value)

    property E0:
        """The activation energy."""
        def __get__(self):
            return self._E0
        def __set__(self, value):
            self._E0 = quantity.Energy(value)

    property V0:
        """The reference potential."""
        def __get__(self):
            return self._V0
        def __set__(self, value):
            self._V0 = quantity.Potential(value)

    property  electrons:
        """The number of electrons transferred."""
        def __get__(self):
            return self._electrons
        def __set__(self, value):
            self._electrons = quantity.Dimensionless(value)

    property alpha:
        """The charge transfer coefficient."""
        def __get__(self):
            return self._alpha
        def __set__(self, value):
            self._alpha = quantity.Dimensionless(value)

    cpdef change_v0(self, double V0):
        """
        Changes the reference potential to `V0` in volts, and adjusts the
        activation energy `E0` accordingly.
        """

        self._E0.value_si = self.get_activation_energy_from_potential(V0,0.0)
        self._V0.value_si = V0

    cpdef double get_rate_coefficient(self, double T, double dHrxn=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3,
        mol, and s at temperature `T` in K and enthalpy of reaction `dHrxn`
        in J/mol.
        """
        cdef double A, n, Ea
        Ea = self.get_activation_energy(dHrxn)
        A = self._A.value_si
        n = self._n.value_si
        return A * T ** n * exp(-Ea / (constants.R * T))

    cpdef double get_activation_energy(self, double dHrxn) except -1:
        """
        Return the activation energy in J/mol corresponding to the given
        enthalpy of reaction `dHrxn` in J/mol.
        """
        cdef double w0, E0, Ea
        E0 = self._E0.value_si
        w0 = self._w0.value_si
        if E0 < 0:
            # Negative E0 is unphysical, but could be encountered during optimization.
            Ea = E0 + max(dHrxn, 0.0)
        else:
            Vp = 2 * w0 * (w0 + E0) / (w0 - E0)
            Ea = (w0 + dHrxn / 2.0) * (Vp - 2 * w0 + dHrxn) * (Vp - 2 * w0 + dHrxn) / (Vp * Vp - (4 * w0 * w0) + dHrxn * dHrxn)
            if (Ea < 0.0) or (dHrxn < -4 * E0):
                Ea = 0.0
            elif (Ea < dHrxn) or (dHrxn > 4 * E0):
                Ea = dHrxn
        return Ea

    cpdef double get_rate_coefficient_from_potential(self, double T, double V, double dHrxn) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3,
        mol, and s at temperature `T` in K, potential `V` in volts, and
        heat of reaction `dHrxn` in J/mol.
        """
        cdef double A, n, Ea
        Ea = self.get_activation_energy_from_potential(V,dHrxn)
        Ea -= self._alpha.value_si * self._electrons.value_si * constants.F * (V-self._V0.value_si)
        A = self._A.value_si
        n = self._n.value_si
        return A * T ** n * exp(-Ea / (constants.R * T))

    def fit_to_reactions(self, rxns, w0=None, recipe=None, Ts=None):
        """
        Fit an ArrheniusChargeTransferBM model to a list of reactions at the given temperatures,
        w0 must be either given or estimated using the family object

        WARNING: there's a lot of code duplication with ArrheniusBM.fit_to_reactions
                 so anything you change here you should probably change there too and vice versa!
        """
        assert w0 is not None or recipe is not None, 'either w0 or recipe must be specified'

        for rxn in rxns:
            if rxn.kinetics._V0.value_si != 0.0:
                rxn.kinetics.change_v0(0.0)

        if Ts is None:
            Ts = np.array([300.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1500.0])
        if w0 is None:
            #estimate w0
            w0s = get_w0s(recipe, rxns)
            w0 = sum(w0s) / len(w0s)
        self.w0 = (w0 * 0.001, 'kJ/mol')

        if len(rxns) == 1:
            rxn = rxns[0]
            dHrxn = rxn.get_enthalpy_of_reaction(298.0)
            A = rxn.kinetics.A.value_si
            n = rxn.kinetics.n.value_si
            Ea = rxn.kinetics.Ea.value_si

            def kfcn(E0):
                Vp = 2 * w0 * (w0 + E0) / (w0 - E0)
                out = Ea - (w0 + dHrxn / 2.0) * (Vp - 2 * w0 + dHrxn) * (Vp - 2 * w0 + dHrxn) / (Vp * Vp - (4 * w0 * w0) + dHrxn * dHrxn)
                return out

            E0 = fsolve(kfcn, w0 / 10.0)[0]

            self.Tmin = rxn.kinetics.Tmin
            self.Tmax = rxn.kinetics.Tmax
            self.comment = 'Fitted to 1 reaction'
        else:
            # define optimization function
            def kfcn(xs, lnA, n, E0):
                T = xs[:,0]
                dHrxn = xs[:,1]
                Vp = 2 * w0 * (w0 + E0) / (w0 - E0)
                Ea = (w0 + dHrxn / 2.0) * (Vp - 2 * w0 + dHrxn) * (Vp - 2 * w0 + dHrxn) / (Vp * Vp - (4 * w0 * w0) + dHrxn * dHrxn)
                Ea = np.where(dHrxn< -4.0*E0, 0.0, Ea)
                Ea = np.where(dHrxn > 4.0*E0, dHrxn, Ea)
                return lnA + np.log(T) * n + (-Ea / (constants.R * T))

            # get (T,dHrxn(T)) -> (Ln(k) mappings
            xdata = []
            ydata = []
            sigmas = []
            E0 = 0.0
            lnA = 0.0
            n = 0.0
            for rxn in rxns:
                # approximately correct the overall uncertainties to std deviations
                s = rank_accuracy_map[rxn.rank].value_si/2.0
                dHrxn = rxn.get_enthalpy_of_reaction(298.0)
                for T in Ts:
                    xdata.append([T, dHrxn])
                    ydata.append(np.log(rxn.get_rate_coefficient(T)))
                    sigmas.append(s / (constants.R * T))
                # Use BEP with alpha = 0.25 for initial guess of E0
                E0 += rxn.kinetics._Ea.value_si - 0.25 * dHrxn
                lnA += np.log(rxn.kinetics.A.value_si)
                n += rxn.kinetics.n.value_si
            # average E0, lnA, n over all reactions
            E0 /= len(rxns)
            lnA /= len(rxns)
            n /= len(rxns)
            E0 = min(E0, w0)
            w0 = max(2 * E0, w0) # ensure w0 is at least 2*E0
            self.w0 = (w0 * 0.001, 'kJ/mol')
            if E0 < 0:
                E0 = w0 / 100.0

            xdata = np.array(xdata)
            ydata = np.array(ydata)

            # fit parameters
            keep_trying = True
            xtol = 1e-8
            ftol = 1e-8
            attempts = 0
            while keep_trying:
                keep_trying = False
                try:
                    params = curve_fit(kfcn, xdata, ydata, sigma=sigmas, p0=[lnA, n, E0], xtol=xtol, ftol=ftol)
                    lnA, n, E0 = params[0].tolist()
                    if abs(E0/self.w0.value_si) > 1.0 and attempts < 5:
                        keep_trying = True
                        if attempts > 0:
                            self.w0.value_si *= 1.25
                            w0 = self.w0.value_si # update local variable for use in kfcn optimization function
                        attempts += 1
                        E0 = self.w0.value_si / 10.0
                except RuntimeError:
                    if xtol < 1.0:
                        keep_trying = True
                        xtol *= 10.0
                        ftol *= 10.0
                    else:
                        raise ValueError("Could not fit BM arrhenius to reactions with xtol<1.0")

            A = np.exp(lnA)

            self.Tmin = (np.min(Ts), "K")
            self.Tmax = (np.max(Ts), "K")
            self.comment = 'Fitted to {0} reactions at temperatures: {1}'.format(len(rxns), Ts)

        # fill in parameters
        A_units = ['', 's^-1', 'm^3/(mol*s)', 'm^6/(mol^2*s)']
        order = len(rxns[0].reactants)
        if order != 1 and rxn.is_surface_reaction():
            raise NotImplementedError("Units not implemented for surface reactions")
        self.A = (A, A_units[order])

        self.n = n
        self.E0 = (E0 * 0.001, 'kJ/mol')
        self._V0.value_si = 0.0
        self.electrons = rxns[0].electrons

        return self

    cpdef ArrheniusChargeTransfer to_arrhenius_charge_transfer(self, double dHrxn):
        """
        Return an :class:`ArrheniusChargeTransfer` instance of the kinetics model using the
        given heat of reaction `dHrxn` to determine the activation energy.
        """
        return ArrheniusChargeTransfer(
            A=self.A,
            n=self.n,
            electrons=self.electrons,
            Ea=(self.get_activation_energy(dHrxn) * 0.001, "kJ/mol"),
            V0=self.V0,
            T0=(1, "K"),
            Tmin=self.Tmin,
            Tmax=self.Tmax,
            Pmin=self.Pmin,
            Pmax=self.Pmax,
            uncertainty=self.uncertainty,
            solute=self.solute,
            comment=self.comment,
        )

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Must match temperature
        and pressure range of kinetics model, as well as parameters: A, n, Ea, T0. (Shouldn't have pressure
        range if it's Arrhenius.) Otherwise returns ``False``.
        """
        if not isinstance(other_kinetics, ArrheniusChargeTransferBM):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if (not self.A.equals(other_kinetics.A) or not self.n.equals(other_kinetics.n)
                or not self.E0.equals(other_kinetics.E0) or not self.w0.equals(other_kinetics.w0)
                or not self.alpha.equals(other_kinetics.alpha)
                or not self.electrons.equals(other_kinetics.electrons) or not self.V0.equals(other_kinetics.V0)):
            return False

        return True

    cpdef change_rate(self, double factor):
        """
        Changes A factor by multiplying it by a ``factor``.
        """
        self._A.value_si *= factor

    def set_cantera_kinetics(self, ct_reaction, species_list):
        """
        Sets a cantera ElementaryReaction() object with the modified Arrhenius object
        converted to an Arrhenius form.
        """
        raise NotImplementedError('set_cantera_kinetics() is not implemented for ArrheniusEP class kinetics.')

cdef class Marcus(KineticsModel):
    """
    A kinetics model based on the (modified) Arrhenius equation, using the
    Evans-Polanyi equation to determine the activation energy. The attributes
    are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `A`             The preexponential factor
    `n`             The temperature exponent
    `lmbd_i_coefs`  Coefficients for inner sphere reorganization energy
    `beta`          Transmission decay coefficient
    `wr`            Work to bring reactants together
    `wp`            Work to bring products together 
    `lmbd_o`        Outer sphere reorganization energy (solvent)
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `solute`        Transition state solute data
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================

    """

    def __init__(self, A=None, n=0.0, lmbd_i_coefs=np.array([0.0,0.0,0.0,0.0]), beta=(1.2e-10,"1/m"), 
                wr=(0,"J/mol"), wp=(0,"J/mol"), lmbd_o=(0,"J/mol"), Tmin=None, Tmax=None,
                Pmin=None, Pmax=None, solute=None, uncertainty=None, comment=''):

        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, solute=solute, uncertainty=uncertainty,
                 comment=comment)

        self.A = A
        self.n = n
        self.lmbd_i_coefs = lmbd_i_coefs
        self.beta = beta 
        self.wr = wr 
        self.wp = wp 
        self.lmbd_o = lmbd_o

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Marcus object.
        """
        string = 'Marcus(A={0!r}, n={1!r}, lmbd_i_coefs={2!r}, beta={3!r}, wr={4!r}, wp={5!r}, lmbd_o={6!r}'.format(
            self.A, self.n, self.lmbd_i_coefs, self.beta, self.wr, self.wp, self.lmbd_o)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.solute: string += ', solute={0!r}'.format(self.solute)
        if self.uncertainty: string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a Marcus object.
        """
        return (Marcus, (self.A, self.n, self.lmbd_i_coefs, self.beta, self.wr, self.wp, self.lmbd_o, 
                            self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                            self.solute, self.uncertainty, self.comment))

    property A:
        """The preexponential factor."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.RateCoefficient(value)

    property n:
        """The temperature exponent."""
        def __get__(self):
            return self._n
        def __set__(self, value):
            self._n = quantity.Dimensionless(value)

    property lmbd_i_coefs:
        """Temperature polynomial coefficients for inner sphere reogranization energy"""
        def __get__(self):
            return self._lmbd_i_coefs
        def __set__(self, value):
            self._lmbd_i_coefs = quantity.Dimensionless(value)

    property beta:
        """transmission coefficient"""
        def __get__(self):
            return self._beta
        def __set__(self, value):
            self._beta = quantity.UnitType('m^-1')(value)

    property lmbd_o:
        """outer sphere reorganization energy"""
        def __get__(self):
            return self._lmbd_o
        def __set__(self, value):
            self._lmbd_o = quantity.Energy(value)

    property wr:
        """outer sphere reorganization energy"""
        def __get__(self):
            return self._wr
        def __set__(self, value):
            self._wr = quantity.Energy(value)

    property wp:
        """outer sphere reorganization energy"""
        def __get__(self):
            return self._wp
        def __set__(self, value):
            self._wp = quantity.Energy(value)

    cpdef double get_rate_coefficient(self, double T, double dGrxn=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3,
        mol, and s at temperature `T` in K and free energy of reaction `dGrxn`
        in J/mol.

        `dGrxn` is a property of the reaction, not of the rate law, so it cannot be defaulted
        here to anything meaningful: the 0.0 below means an exactly thermoneutral reaction. Callers
        holding a reaction should route through :meth:`rmgpy.reaction.Reaction.get_rate_coefficient`,
        which resolves `dGrxn` from the species thermochemistry.
        """
        cdef double A, n, dG
        dG = self.get_gibbs_activation_energy(T, dGrxn)
        A = self._A.value_si
        n = self._n.value_si
        return A * T ** n * exp(-dG / (constants.R * T))

    cpdef double get_lmbd_i(self, double T):
        """
        Return lmbd_i in J/mol
        """
        return self.lmbd_i_coefs.value_si[0]+self.lmbd_i_coefs.value_si[1]*T+self.lmbd_i_coefs.value_si[2]*T**2+self.lmbd_i_coefs.value_si[3]*T**3
    
    cpdef double get_gibbs_activation_energy(self, double T, double dGrxn) except -1:
        """
        Return the Gibbs activation energy in J/mol corresponding to the given
        free energy of reaction `dGrxn` in J/mol.
        """
        cdef double lmbd_i
        lmbd_i = self.get_lmbd_i(T)
        return (lmbd_i+self.lmbd_o.value_si)/4.0*(1.0+dGrxn/(lmbd_i+self.lmbd_o.value_si))**2

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Must match temperature
        and pressure range of kinetics model, as well as parameters: A, n, Ea, T0. (Shouldn't have pressure
        range if it's Arrhenius.) Otherwise returns ``False``.
        """
        if not isinstance(other_kinetics, Marcus):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if (not self.A.equals(other_kinetics.A) or not self.n.equals(other_kinetics.n)
                or not self.lmbd_i_coefs.equals(other_kinetics.lmbd_i_coefs) or not self.lmbd_o.equals(other_kinetics.lmbd_o)
                or not self.beta.equals(other_kinetics.beta)
                or not self.electrons.equals(other_kinetics.electrons)):
            return False

        return True

    cpdef change_rate(self, double factor):
        """
        Changes A factor by multiplying it by a ``factor``.
        """
        self._A.value_si *= factor

    def set_cantera_kinetics(self, ct_reaction, species_list):
        """
        Sets a cantera ElementaryReaction() object with the modified Arrhenius object
        converted to an Arrhenius form.
        """
        raise NotImplementedError('set_cantera_kinetics() is not implemented for Marcus class kinetics.')

def get_w0(actions, rxn):
    """
    calculates the w0 for Blower Masel kinetics by calculating wf (total bond energy of bonds formed)
    and wb (total bond energy of bonds broken) with w0 = (wf+wb)/2
    """
    mol = None
    for r in rxn.reactants:
        m = r.molecule[0]
        if mol:
            mol = mol.merge(m)
        else:
            mol = m.copy(deep=True)
    a_dict = mol.get_all_labeled_atoms()

    recipe = actions

    wb = 0.0
    wf = 0.0
    for act in recipe:

        if act[0] in ['BREAK_BOND','FORM_BOND','CHANGE_BOND']:

            if act[1] == act[3]: # the labels are the same
                atom1 = a_dict[act[1]][0]
                atom2 = a_dict[act[3]][1]
            else:
                atom1 = a_dict[act[1]]
                atom2 = a_dict[act[3]]

        if act[0] == 'BREAK_BOND':
            bd = mol.get_bond(atom1, atom2)
            wb += bd.get_bde()
        elif act[0] == 'FORM_BOND':
            bd = Bond(atom1, atom2, act[2])
            wf += bd.get_bde()
        elif act[0] == 'CHANGE_BOND':
            bd1 = mol.get_bond(atom1, atom2)

            if act[2] + bd1.order == 0.5:
                mol2 = None
                for r in rxn.products:
                    m = r.molecule[0]
                    if mol2:
                        mol2 = mol2.merge(m)
                    else:
                        mol2 = m.copy(deep=True)
                a_dict_mol2 = mol2.get_all_labeled_atoms()
                if act[1] == act[3]: # the labels are the same
                    atom1_mol2 = a_dict_mol2[act[1]][0]
                    atom2_mol2 = a_dict_mol2[act[3]][1]
                else:
                    atom1_mol2 = a_dict_mol2[act[1]]
                    atom2_mol2 = a_dict_mol2[act[3]]

                bd2 = mol2.get_bond(atom1_mol2, atom2_mol2)
            else:
                bd2 = Bond(atom1, atom2, bd1.order + act[2])

            if bd2.order == 0:
                bd2_bde = 0.0
            else:
                bd2_bde = bd2.get_bde()
            bde_diff = bd2_bde - bd1.get_bde()
            if bde_diff > 0:
                wf += abs(bde_diff)
            else:
                wb += abs(bde_diff)
    return (wf + wb) / 2.0

def get_w0s(actions, rxns):
    return [get_w0(actions, rxn) for rxn in rxns]