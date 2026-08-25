# M7-PREFLIGHT -- LITHIUM CHARGED-PATH INTEGRATION
#
# Neutral lithium-bearing feed (methyllithium) driven through a two-temperature
# plasma reactor.  The composition is NEUTRAL: no cation is seeded.  The
# lithium-specific Cation_* families are loaded so that any cation is produced by
# RMG from the neutral feed rather than seeded.
#
# EXCLUDED -- Cation_R_Recombination.  It has been recovered as what it actually
# is: a lithium-ion-battery SEI family, Marcus electron transfer at a Li(110)
# electrode in liquid ethylene carbonate at 298.15 K under declared electrode
# potentials.  It is not plasma chemistry; read as plasma kinetics its rates
# evaluate 30-230 orders of magnitude below anything physical, and a binding
# ruling excludes it from every plasma configuration.  The family and its data
# are preserved deliberately as provenance evidence and remain correct for the
# electrochemical domain -- see examples/rmg/SEI_pure_ACN/input.py, which keeps
# it.  Only this deck's reference to it is removed.
#
# Consequence, recorded rather than worked around: the reverse of that family's
# template was this deck's ONLY route from the neutral feed to Li+, and the only
# family here that fired at all.  With it gone the deck generates no reactions
# from CH3Li.  The charged-path exercise this deck was built for is retired, not
# relocated; restoring it needs a genuine plasma ionisation route (e.g.
# electron-impact ionisation), not this family.
#
# Pinned by test/rmgpy/preflightDeckFamilyExclusionTest.py.
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
        'primaryThermoLibrary', 'electrocatThermo',
    ],
    reactionLibraries=[],
    seedMechanisms=[],
    kineticsDepositories=['training'],
    kineticsFamilies=[
        # the lithium-specific cation families (produce/consume Li+).
        # Cation_R_Recombination is deliberately absent -- see the header.
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
    structure=adjacencyList("1 e u0 p0 c-1"),
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
