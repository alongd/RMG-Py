Argon resonance radiation
========================

Resolved neutral Ar ``electronicstate 3P1`` (1s4) and ``1P1`` (1s2), in either
core or edge, require an explicit reactor declaration. Labels do not identify
these states. For example, add these keywords to ``plasmaReactor``::

    chamberGeometry={'shape': 'cylinder', 'radius': (5, 'cm'), 'length': (30, 'cm')},
    resonanceRadiation={'model': 'holsteinWalsh', 'geometryModel': 'infiniteCylinder'},

``model='complete'`` selects exactly zero escape and zero radiative loss.
Both choices require geometry. ``geometryModel='slab'`` instead requires
``chamberGeometry={'shape': 'slab', 'gap': (5, 'cm')}``. The optical length is
the cylinder radius or full slab gap; ``diffusionLength`` is insufficient.
The cylinder approximation neglects axial escape and requires length/radius
at least six. This is an engineering scope limit, not a finite-cylinder error
guarantee. Geometry may be declared without charged-particle wall transport.

The operator removes each upper state at :math:`A g n_u` and returns one ground
Ar atom. It changes neither electrons nor metastables. It does not create a
mechanism reaction. Core growth rebinds the molecular identities and rejects
duplicate resonance states or missing/ambiguous ground Ar. No extra neutral
wall diffusion is inferred. An explicitly declared wall loss remains additive.

Closure and domain
------------------

The closure assumes complete frequency redistribution, a homogeneous stationary
gas, nonreflecting boundaries and the fundamental imprisonment mode. It uses
Holstein (1951), Tables I and III, and Walsh (1959), Eqs. 2.7--2.9. Define
:math:`\Delta\nu_D=\sqrt{2kT/m_{Ar}}/\lambda`,
:math:`\tau=\lambda^2 N_{Ar} (g_u/g_l) A\ell/(8\pi^{3/2}\Delta\nu_D)`, and
:math:`a=\Delta\nu_L/(2\Delta\nu_D)`. The ground-Ar density is the reactor's
own partial density, not its total gas density. Natural FWHM :math:`A/(2\pi)`
and resonance self-broadening from NIST's Spectral Line Shapes Eq. 41 supply
:math:`\Delta\nu_L`. The sourced vacuum lines are in the database's
``eedf/argon_resonance_radiation.json``:

.. list-table:: Vacuum resonance lines
   :header-rows: 1

   * - State
     - Wavelength (nm)
     - A (s\ :sup:`-1`)
   * - 1s4 / 3P1
     - 106.6659865
     - 1.320e8
   * - 1s2 / 1P1
     - 104.8219924
     - 5.32e8

The combined escape factor is

.. math::

   g_D &= C_D/[\tau\sqrt{\pi\ln\tau}],\\
   g_C &= C_C\sqrt{a/(\sqrt{\pi}\tau)},\\
   g_{CD} &= 2a/[\pi\sqrt{\ln\tau}],\\
   q &= g_{CD}/g_C,\\
   g &= g_D\exp(-q^2)+g_C\operatorname{erf}(q).

For the cylinder :math:`(C_D,C_C)=(1.60,1.115)`; for the slab these are
:math:`(1.875,1.150)`. Walsh's explicit lifetime interpolation is cylindrical;
the slab adopts that interpolation with Holstein's slab coefficients and the
leading-log mixed coefficient one. It is an asymptotic extension, not an exact
slab Voigt solution. The zero-damping limit is evaluated analytically.
Bogaerts et al. (1998), Sec. III.D, use related transmission approximations;
their 1.9/1.3 coefficients are not substituted for these lifetime coefficients.

``holsteinWalsh`` refuses :math:`\tau<10^4`, :math:`a` outside [0, 0.5], gas
T outside [250, 500] K, ground-Ar fraction below 0.75 of the neutral bath,
foreign neutrals other than N2 above a trace fraction of 1e-6, heavy charged
fraction above 1e-3, and :math:`n_u/(3n_{Ar})>10^{-3}` for either line.
These are declared engineering limits, not a theorem of transport accuracy.
There is no clamp or runtime numerical transport operator and no continuation
to thin or crossover gas. Complete trapping bypasses those physical limits
but still validates identity and geometry.

N2 admixture and the omitted linewidth
-------------------------------------

Up to 20% N2 is supported with the Ar-partial-density opacity above. N2 foreign
broadening is omitted with the following quantified approximation, not a claim
of measured VUV linewidth accuracy. Yatom et al. (2023), Eq. 15, give the impact
estimate (lambda and width in cm, density in cm\ :sup:`-3`, polarizability in
cm\ :sup:`3`, squared radius in Bohr-radius units, reduced mass in u):

.. math::

   w_W=8.18\,10^{-12}\lambda^2(\alpha\overline{R^2})^{0.4}
       (T/\mu)^{0.3}N_{N2}.

NIST CCCBDB gives mean N2 polarizability 1.710 Angstrom cubed. The NIST Coulomb
radius estimate gives upper-minus-ground squared radii 25.17 and 28.04 for these
lines. Using 30 for both, T <= 500 K and reduced mass 16.47 u bounds this
*impact estimate* by 8.3e-10 Hz cm\ :sup:`3` per N2 perturber. Compared with
self-broadening coefficients 1.07249e-9 and 4.10213e-9 Hz cm\ :sup:`3`, at
20% N2 the omitted width is at most 19.4% and 5.1% of the Ar resonance width.
At 5 torr and radius 5 cm, including this estimate changes either escape factor
by less than 10% at the worst-width temperature of 500 K. At the wider domain boundary N2/ground-Ar <= 1/3, the corresponding
linewidth bounds are 25.8% and 6.8%. Across the whole asymptotic closure the
relative escape-factor increase is bounded by the relative linewidth increase:
``g(a)/a`` decreases with positive ``a``. In the pressure limit the increase is
only ``sqrt(1 + width_ratio) - 1``.

This is a bound on an approximate impact model. It excludes uncertainty in the
Coulomb-radius estimate, molecular anisotropy and other collisional effects;
no measured Ar--N2 VUV coefficient was established here. Treat an admixture
result requiring better than this accuracy as outside the qualification of the
omission. Quenching and state mixing are separate physics and are not supplied
by this radiation operator.

Energy, diagnostics and export
------------------------------

Radiation carries :math:`A g n_u h c/\lambda` out of the gas or into the wall.
Gas temperature is prescribed; no gas thermal-energy equation receives it.
Electron excitation has already charged the electron energy budget. This
operator therefore adds no second electron cooling term. Photon power is
reported separately for the solver's EOS volume (the CSV supplies that volume).
A fixed-temperature reactor's one-mole normalization need not equal the physical
chamber volume. The residual and Jacobian use each trial state's ground
density and volume, including electron-temperature-dependent volume where
applicable. Diagnostics are latched only at accepted states.

``simulation_*.radiation.json`` accompanies the native CSV, carrying the model,
physical geometry, NIST provenance, accepted time, state identities, escape
factors, effective rates, densities and photon powers. Its trajectory aligns
with the accepted CSV snapshots. Pickle and saved input preserve the declaration.

Chemkin, Cantera and RMS mechanism files omit reactor-level operators, including
wall loss and resonance radiation, and name the omissions in the existing run
provenance header. They contain no vacuum-A reaction or frozen escape rate.
Ordinary Chemkin output remains enabled; a mechanism file alone does not
reproduce the reactor trajectory.

References
----------

* Holstein (1947), `Phys. Rev. 72, 1212 <https://doi.org/10.1103/PhysRev.72.1212>`_.
* Holstein (1951), `Phys. Rev. 83, 1159 <https://doi.org/10.1103/PhysRev.83.1159>`_.
* Walsh (1959), `Phys. Rev. 116, 511 <https://doi.org/10.1103/PhysRev.116.511>`_.
* Bogaerts, Gijbels and Vlcek (1998), `J. Appl. Phys. 84, 121 <https://doi.org/10.1063/1.368009>`_.
* NIST ASD v5.12 (2024), `doi:10.18434/T4W30F <https://doi.org/10.18434/T4W30F>`_.
  The database records the query, retrieval time, original response digest and
  per-line ASD transition reference T7574; no expanded citation is inferred.
* `NIST Spectral Line Shapes <https://physics.nist.gov/Pubs/AtSpec/node20.html>`_, Eq. 41 and Coulomb-radius estimate.
* Yatom et al. (2023), `Plasma Sources Sci. Technol. 32, 115005 <https://doi.org/10.1088/1361-6595/ad0741>`_, Eq. 15.
* `NIST CCCBDB experimental polarizabilities <https://cccbdb.nist.gov/Pollistx.Asp>`_.
