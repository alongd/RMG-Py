# M7-PREFLIGHT -- LITHIUM CHARGED-PATH INTEGRATION
#
# Neutral lithium-bearing feed (methyllithium) driven through a two-temperature
# plasma reactor.  The composition is NEUTRAL: no cation is seeded.  The six
# lithium-specific Cation_* families are loaded so the cation is produced by RMG
# from the neutral feed via the auto-derived reverse template of
# Cation_R_Recombination (Li+ + R.  <=>  R-Li), exercised in reverse.
#
# Run (read-only against RMG-database-plasma):
#   conda activate rmg_env
#   export PYTHONPATH=/home/alon/Code/RMG-Py-m7-preflight:$PYTHONPATH
#   python /home/alon/Code/RMG-Py-m7-preflight/rmg.py \
#       /home/alon/Code/RMG-Py-m7-preflight/docs/m7-preflight/input.py
# database.directory is taken from docs/m7-preflight/rmgrc when run from that dir,
# or set it in ~/.rmg/rmgrc -> /home/alon/Code/RMG-database-plasma/input

database(
    thermoLibraries=[
        'LithiumPrimaryThermo', 'LithiumAdditionalThermo',
        'primaryThermoLibrary',
        # electrocatThermo carries the electron thermo entry ('1 e u0 p0 c-1',
        # H298=0). I-085 note: the original preflight primary deck omitted this,
        # so NO electron thermo library was loaded and the electron missed every
        # library regardless of u0/u1 -- a deck-config gap stacked on top of the
        # actual Wall A (the u1-vs-u0 lookup mismatch). Loaded here so the driver
        # exercises the prompt's Wall A (u1 electron, electron entry present).
        'electrocatThermo',
    ],
    reactionLibraries=[],
    seedMechanisms=[],
    kineticsDepositories=['training'],
    kineticsFamilies=[
        # the six lithium-specific cation families (produce/consume Li+)
        'Cation_R_Recombination',
        'Cation_Li_Abstraction',
        'Cation_Addition_MultipleBond',
        'Cation_Addition_MultipleBond_Disprop',
        'Cation_NO_Substitution',
        'Cation_NO_Ring_Opening',
        # neutral lithium families that the neutral feed can also fire
        'Li_Abstraction',
        'Li_Addition_MultipleBond',
    ],
    kineticsEstimator='rate rules',
)

# --- electron pseudo-species MUST be declared before plasmaReactor(...) ---
species(
    label='e-',
    reactive=True,
    structure=adjacencyList("1 e u1 p0 c-1"),
)

# --- neutral lithium-bearing feed: methyllithium, CH3Li, net charge 0 ---
species(
    label='CH3Li',
    reactive=True,
    structure=SMILES("C[Li]"),
)

# --- inert bath gas ---
species(
    label='He',
    reactive=False,
    structure=SMILES("[He]"),
)

plasmaReactor(
    temperature=(1000, 'K'),
    pressure=(1, 'atm'),
    electronTemperature=(11604.5, 'K'),          # ~1 eV
    electronDensity=(1e17, 'm^-3'),              # nonzero initial electron density
    initialMoleFractions={
        'CH3Li': 0.05,
        'He': 0.95,
    },
    terminationTime=(1e-3, 's'),
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
)
