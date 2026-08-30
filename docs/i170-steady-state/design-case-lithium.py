# I-170 -- LITHIUM-SEEDED ARGON PLASMA, 5 TORR, 3 eV
#
# A DESIGN CASE, not a deliverable deck. The I-164 argon deck's trajectory is exactly
# flat (the database carries no argon chemistry), so a steady-state tolerance cannot be
# argued from it. This deck keeps every I-164 condition -- 5 torr, 298.15 K gas, 3 eV
# electrons, n_e = 1e16 m^-3 -- and seeds a trace of lithium, which is the ONE element
# for which the plasma libraries actually carry both directions of the ionisation
# balance:
#
#   PlasmaElectronImpactIonization : [Li] => [Lip]   (Voronov)
#   PlasmaRadiativeRecombination   : [Lip] => [Li]   (Badnell)
#
# An alkali-seeded noble gas is a standard laboratory configuration, so this is a real
# case rather than a contrivance: the ionisation balance it relaxes to is the physical
# endpoint a steady-state criterion has to detect.
#
# Charge neutrality at t = 0 uses the same hand arithmetic as the I-164 deck: the
# neutralising cation's PRE-SCALE mole fraction is r = a/(1-q), which lands on
# x_cation = x_e after the driver's heavy_scale. Same P, Te, n_e as I-164, so the same
# r = 6.175165731685617e-08.

database(
    thermoLibraries=[
        'primaryThermoLibrary',
        'LithiumPrimaryThermo',
        'PlasmaCationThermo',
        'electrocatThermo',
    ],
    reactionLibraries=[
        'PlasmaElectronImpactIonization',
        'PlasmaRadiativeRecombination',
    ],
    seedMechanisms=[],
    kineticsDepositories=['training'],
    kineticsFamilies=[
        'Plasma_Electron_Attachment',
    ],
    kineticsEstimator='rate rules',
)

species(
    label='e-',
    reactive=True,
    structure=adjacencyList("1 e u1 p0 c-1"),
)

species(
    label='Ar',
    reactive=True,
    structure=adjacencyList("""
1 Ar u0 p4 c0
"""),
)

species(
    label='Li',
    reactive=True,
    structure=adjacencyList("""
multiplicity 2
1 Li u1 p0 c0
"""),
)

species(
    label='Lip',
    reactive=True,
    structure=adjacencyList("""
1 Li u0 p0 c+1
"""),
)

plasmaReactor(
    temperature=(298.15, 'K'),
    pressure=(5, 'torr'),
    electronTemperature=(34813.5, 'K'),     # 3 eV
    electronDensity=(1e16, 'm^-3'),         # 1e10 cm^-3
    initialMoleFractions={
        'Lip': 6.175165731685617e-08,       # r = a/(1-q); neutralises the seed electron
        'Li': 1.0e-4,
        'Ar': 1.0 - 1.0e-4 - 6.175165731685617e-08,
    },
    terminationTime=(1e2, 's'),
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
    saveSimulationProfiles=True,
    generateRMSYAML=False,
    generateCanteraYAML2=False,
)
