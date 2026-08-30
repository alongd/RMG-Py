restartFromSeed(path='seed')

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
# 2. "QUASI-STEADY STATE" IS NOT EXPRESSIBLE. `plasmaReactor` accepts only
#    terminationConversion / terminationTime / terminationRateRatio
#    (rmgpy/rmg/input.py:518-525). None of them means "integrate until d/dt -> 0".
#    Substituted: terminationTime = 1e-3 s, orders above the electron kinetics
#    timescale at this density, and stated here as a substitution.
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

# CONTROL RUN 2. Identical to ../input.py except that Ar+ is neither declared nor
# seeded: the deck is written the OBVIOUS way, with `electronDensity` alone. The
# directive inserts x_e and creates no compensating cation, so the mixture carries a
# net charge of -x_e at t = 0. This control measures whether anything in RMG or the
# PlasmaReactor notices.

plasmaReactor(
    temperature=(298.15, 'K'),
    pressure=(5, 'torr'),
    electronTemperature=(34813.5, 'K'),     # 3 eV
    electronDensity=(1e16, 'm^-3'),         # 1e10 cm^-3
    initialMoleFractions={
        'Ar': 1.0,
    },
    terminationTime=(1e-3, 's'),            # substitute for "quasi-steady state"
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
