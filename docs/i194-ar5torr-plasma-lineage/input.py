# I-164 -- 5 TORR PURE-ARGON PLASMA, DIAGNOSTIC PROBE
#
# The case, as specified:
#
#   5 torr argon at room temperature, ARGON ONLY, with a seed electron density;
#   neutral at t = 0, achieved by converting part of the argon to Ar+; integrated to
#   quasi-steady state.
#
# This deck exists to be RUN, not to be believed: it is the instrument of a diagnostic
# probe that measures what stops RMG from simulating this configuration today. Nothing
# here is tuned to make the run pass.
#
# Run (from the directory holding this file, which must also hold an rmgrc pinning
# database.directory at an RMG-database-plasma checkout):
#
#   conda activate rmg_env
#   python /home/alon/Code/RMG-Py-plasma/rmg.py input.py \
#     > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)
#
# ---------------------------------------------------------------------------------
# HOW FAITHFULLY THIS EXPRESSES THE CASE
# ---------------------------------------------------------------------------------
#
# 1. CHARGE-NEUTRAL SEEDING IS NOT A DIRECTIVE; THE ARITHMETIC IS DONE BY HAND BELOW.
#    `electronDensity` inserts an electron mole fraction x_e and scales every heavy
#    species down by `heavy_scale` so the composition still sums to one
#    (rmgpy/rmg/input.py:713-732). It creates no compensating cation, so a deck that
#    uses it alone is net-negative at t = 0. To be neutral the modeller must supply
#    Ar+ explicitly AND pre-divide by heavy_scale, i.e. reproduce input.py's own
#    internal derivation outside it. With
#
#        kB = R/Na, q = n_e kB Te / P, a = q Tgas / Te
#        x_e = a / (1 - q + a),  heavy_scale = (1 - q) / (1 - q + a)
#
#    the heavy-side Ar+ fraction that lands on x_Ar+ = x_e after scaling is
#
#        r = a / (1 - q)
#
#    For P = 5 torr (666.6118421052631 Pa), Tgas = 298.15 K, Te = 34813.5 K,
#    n_e = 1e16 m^-3:
#
#        q   = 7.210383435118683e-06
#        a   = 6.175121206372917e-08
#        x_e = 6.175165350358923e-08
#        r   = 6.175165731685617e-08     <- the number typed below
#
#    This is exact, not approximate: r * heavy_scale = a/(1-q+a) = x_e identically.
#    It is nonetheless a finding, not a workaround: the API has no charge-neutral
#    seeding directive, and getting neutrality right requires re-deriving the
#    driver's internals inside the deck.
#
# 2. THERE IS NO QUASI-STEADY STATE, AND NO FINITE EXHAUSTION TIME; the stop is a
#    tolerance CUTOFF near complete ionisation (I-199). This model can only ionise
#    (Ar + e- => Arp + 2e-, irreversible; no recombination, no wall loss), so it has no
#    stationary state; mass action drives it asymptotically toward complete ionisation.
#    The reactor is constant-pressure and two-temperature (V = R(N_heavy*Tgas + Ne*Te)/P,
#    plasma.pyx:831; C = N/V, plasma.pyx:931), so in the ionised fraction f the governing
#    equation is NOT a fixed-volume logistic but a Riccati form,
#        df/dt = (k_iz*P/R) * f(1-f)/(Tgas + Te*f),   k_iz MOLAR.
#    Its initial slope is nu_iz = k_iz*P/(R*Tgas) = k_iz*n_Ar/Na = 2.23e7/s (e-folding
#    tau = 44.8 ns) at Te = 3 eV (k_iz = 8.30e7 m^3/mol/s from ElectronCollisionPlasma.-
#    integrate_rate_coefficient; n_Ar = 1.62e23 m^-3 at 5 torr / 298.15 K); its final
#    approach 1-f ~ exp(-t/5.28us) is ~117x slower (kP/(R(Tgas+Te))). n_e saturates on
#    its two-temperature plateau 1.375e21 m^-3 = n_Ar*Tgas/(Tgas+Te), NOT n_Ar.
#    terminationTime = 1e-4 s is a CUTOFF, justified as: residual neutral fraction 1-f < 1e-8
#    (reached 9.80e-5 s) AND every core quantity constant to >=4 sig figs (n_e from 2.45e-5
#    s). It is ~2.2e3 initial e-folds -- but the justification is those tolerances, not the
#    count, and it is NOT an instant of physical completion. The old 1e-3 s was ~10x past
#    this cutoff and mislabelled a "QSS substitute". Derivation/trajectory/limitations:
#    docs/i199-stoptime/report.md.
#    (terminationRateRatio would ALSO fire -- at 8.79 us, char_rate/max_char_rate = 0.0094 --
#    but NOT for the obvious reason: char_rate is in concentration units and ~ f(1-f)/
#    (Tgas+Te*f)^2, so it peaks at f* = Tgas/(2Tgas+Te) = 0.0084 and its decline is driven by
#    the two-temperature volume dilution, not reactant depletion (argon is only ~80% consumed
#    at the trigger). It stops mid-avalanche and its "reached RateRatio" log reads misleadingly
#    like convergence, so it is not used. An earlier record that it "never fires because the
#    model has no sink" is wrong -- dilution collapses char_rate regardless of any sink.)
#    `plasmaReactor` accepts only terminationConversion / terminationTime /
#    terminationRateRatio (rmgpy/rmg/input.py); none of them means "integrate until d/dt->0".
#
# 3. ROOM TEMPERATURE = 298.15 K; ELECTRON TEMPERATURE MUST BE KELVIN. The directive
#    rejects electronvolts outright, so 3 eV is typed as its kelvin equivalent
#    34813.5 K (= 3 * 11604.5). 3 eV is a standard low-pressure argon discharge Te.
#
# 4. SEED ELECTRON DENSITY = 1e16 m^-3 (1e10 cm^-3), representative of a few-torr
#    argon glow. Nothing in the case fixed it; it is stated, not derived.
#
# ---------------------------------------------------------------------------------

database(
    thermoLibraries=[
        # neutral Ar ('1 Ar u0 p4 c0')
        'primaryThermoLibrary',
        # Ar+ as a free monatomic cation ('[Arp]', '1 Ar u1 p3 c+1'), NIST-JANAF Ar-002
        'PlasmaCationThermo',
        # the electron thermo entry ('1 e u0 p0 c-1', H298 = 0). Without it the electron
        # misses every library and the reactor stops at the thermo wall.
        'electrocatThermo',
    ],
    reactionLibraries=[
        # what an argon plasma needs: Ar + e- => Ar+ + 2 e-, and Ar+ + e- => Ar + hv.
        # Both libraries are requested here as a modeller would request them; whether
        # either actually carries argon is one of the things this probe measures.
        'PlasmaElectronImpactIonization',
        'PlasmaRadiativeRecombination',
        # PlasmaArgon carries the single argon reaction the database holds,
        # Ar + e- => Arp + e- + e-, in its own file so an argon deck no longer
        # loads the whole PlasmaAir air library to reach it (I-193). PlasmaAir
        # was removed here: every wall the first-light run hit came from the air
        # ballast, never from the argon chemistry.
        'PlasmaArgon',
    ],
    seedMechanisms=[],
    kineticsDepositories=['training'],
    # An EMPTY list is falsy and KineticsDatabase.load_families reads it as "load the
    # default set" -- which would drag in Cation_R_Recombination, the quarantined
    # battery-SEI family. So an explicit list is required. Plasma_Electron_Attachment
    # is the only plasma-generic neutral family in the database and its root is a
    # LogicOr over O / OH / O2, so it cannot match argon; it is named to keep the list
    # truthy and non-default, not because it is expected to fire.
    kineticsFamilies=[
        'Plasma_Electron_Attachment',
    ],
    kineticsEstimator='rate rules',
)

# PlasmaArgon carries only the argon ionisation reaction, so the air species that forced
# a generatedSpeciesConstraints block are gone: no CH (maximumCarbeneRadicals) and no
# singlet O2s (allowSingletO2) enter the model. Both settings existed only to tolerate the
# air library and were removed with it (I-193), restoring RMG's default constraints in full.

# --- electron pseudo-species MUST be declared before plasmaReactor(...) ---
species(
    label='e-',
    reactive=True,
    structure=adjacencyList("1 e u1 p0 c-1"),
)

# --- neutral argon ---
species(
    label='Ar',
    reactive=True,
    structure=adjacencyList("""
1 Ar u0 p4 c0
"""),
)

# --- argon cation: the ionized fraction that makes the seed charge-neutral ---
species(
    label='Arp',
    reactive=True,
    structure=adjacencyList("""
multiplicity 2
1 Ar u1 p3 c+1
"""),
)

plasmaReactor(
    temperature=(298.15, 'K'),
    pressure=(5, 'torr'),
    electronTemperature=(34813.5, 'K'),     # 3 eV
    electronDensity=(1e16, 'm^-3'),         # 1e10 cm^-3
    initialMoleFractions={
        # r = a/(1-q); see note 1 in the header. After the driver's heavy_scale this
        # lands on x_Ar+ = x_e exactly, so the mixture is charge-neutral at t = 0.
        'Arp': 6.175165731685617e-08,
        'Ar': 1.0 - 6.175165731685617e-08,
    },
    terminationTime=(1e-4, 's'),            # cutoff near complete ionisation: 1-f<1e-8 and n_e stationary to >=4 sig figs; not QSS, not a finite exhaustion time (I-199, note 2)
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.5,
    toleranceInterruptSimulation=1e8,
    maximumEdgeSpecies=200,
)

options(
    units='si',
    generateOutputHTML=False,
    generatePlots=False,
    saveEdgeSpecies=True,
    saveSimulationProfiles=False,
    # The RMS writer has no branch for the plasma rate laws (VoronovEIArrhenius /
    # BadnellRRArrhenius) and would kill the run on the first save; off, as on every
    # other plasma deck in this repo.
    generateRMSYAML=False,
    generateCanteraYAML2=True,
)
