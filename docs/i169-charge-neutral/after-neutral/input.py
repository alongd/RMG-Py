# I-169 -- CHARGE-NEUTRAL SEEDING VIA chargeBalanceSpecies
#
# Verifier check 2. This is the SAME 5 torr argon case as
# ../before-nonneutral/input.py, with one line added: `chargeBalanceSpecies='Arp'`.
#
# THE POINT OF THIS DECK IS WHAT IT DOES *NOT* CONTAIN. There is no hand-derived
# constant anywhere in it. Compare the I-164 neutral deck, which had to type
#
#     'Arp': 6.175165731685617e-08,
#     'Ar':  1.0 - 6.175165731685617e-08,
#
# a number obtained by re-deriving rmgpy/rmg/input.py's own q / x_e / heavy_scale
# algebra outside it. Here the modeller types 'Ar': 1.0 and names the ion that carries
# the compensating charge; RMG computes its mole fraction and logs both that fraction
# and the net charge it achieved.
#
# Run (from this directory, which holds its own rmgrc):
#
#   PYTHONPATH=/home/alon/Code/RMG-Py-i169-charge-neutral \
#   python /home/alon/Code/RMG-Py-i169-charge-neutral/rmg.py input.py \
#     > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)
#
# Free parameters, stated rather than derived, carried over from I-164:
# Te = 34813.5 K (3 eV) and n_e = 1e16 m^-3 (1e10 cm^-3).

database(
    thermoLibraries=[
        # neutral Ar ('1 Ar u0 p4 c0')
        'primaryThermoLibrary',
        # Ar+ as a free monatomic cation ('[Arp]', '1 Ar u1 p3 c+1'), NIST-JANAF Ar-002
        'PlasmaCationThermo',
        # the electron thermo entry ('1 e u0 p0 c-1', H298 = 0)
        'electrocatThermo',
    ],
    reactionLibraries=[
        'PlasmaElectronImpactIonization',
        'PlasmaRadiativeRecombination',
    ],
    seedMechanisms=[],
    kineticsDepositories=['training'],
    # An EMPTY list is falsy and is read as "load the default set", which would drag in
    # the quarantined Cation_R_Recombination family; an explicit list is required.
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

# --- argon cation: named below as the species that absorbs the balancing charge.
#     It is declared but carries NO mole fraction of its own; RMG computes one. ---
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
    chargeBalanceSpecies='Arp',             # <- the whole change vs the control deck
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
    generateRMSYAML=False,
    generateCanteraYAML2=True,
)
