# I-123 -- THE LITHIUM CHARGE NETWORK THROUGH THE REAL INPUT-FILE PATH
#
# The companion probe (probe_lithium_charge_network.py) drives the network by calling the
# loader, the resolver and PlasmaReactor directly. This deck drives the SAME network the
# way a user would: an input file, rmg.py, the plasma reactor. Both are needed -- the probe
# proves each stage in isolation, the deck proves the input-file path that carries them.
#
# WHY THIS DECK AND NOT docs/m7-preflight/input.py. That deck feeds methyllithium and, since
# i115-preflight-deck removed Cation_R_Recombination from it, generates no reactions at all:
# the reverse of that family's template was its only route from the neutral feed to Li+, and
# its own header records that restoring the charged path "needs a genuine plasma ionisation
# route (e.g. electron-impact ionisation), not this family." That route is exactly what the
# i119 chain (RMG-Py) and the i114/i119 chain (RMG-database) deliver, and it is present for
# the first time on this integration branch. So this deck is the m7-preflight deck's own
# stated successor, not a replacement for it -- and the reason it can exist at all is the
# union under audit.
#
# The chemistry, both directions, from the shipped tables and not authored here:
#
#   source   Li  + e-  =>  Li+ + 2 e-     Voronov (1997) Z=3 N=3, via PlasmaElectronImpactIonization
#   sink     Li+ + e-  =>  Li  + hv       Badnell (2006) Z=3 N=2, via PlasmaRadiativeRecombination
#
# The composition is NEUTRAL: no cation is seeded. Li+ has to be produced from the neutral
# feed, which is the claim this deck exists to exercise.
#
# Run from the RMG-Py integration worktree root, so rmgrc resolves to the integration
# database worktree:
#
#   conda activate rmg_env
#   python rmg.py docs/i123-integration/input.py

database(
    thermoLibraries=[
        # [Li] and [Lip] as free monatomic species.
        'LithiumPrimaryThermo', 'LithiumAdditionalThermo',
        'primaryThermoLibrary',
        # carries the electron thermo entry ('1 e u0 p0 c-1', H298 = 0). Without it the
        # electron misses every library and the reactor stops at the thermo wall -- that
        # was I-085's deck-config gap, recorded in docs/m7-preflight/.
        'electrocatThermo',
    ],
    reactionLibraries=[
        'PlasmaElectronImpactIonization',
        'PlasmaRadiativeRecombination',
    ],
    seedMechanisms=[],
    kineticsDepositories=['training'],
    # The charge network under audit is carried entirely by the two libraries above. The
    # natural way to say that is `kineticsFamilies=[]`, and it is WRONG: an empty list is
    # falsy, and KineticsDatabase.load_families reads a falsy `families` as "load the
    # default set", so the deck silently gets EVERY family -- including
    # Cation_R_Recombination, the quarantined battery-SEI family that i111/i102 exist to
    # keep out of plasma decks. Measured, not assumed:
    # test_plasma_deck_generates_no_reaction_from_family fails on this deck with
    # `resolves to a family set that loads Cation_R_Recombination. Declared: []`.
    # `kineticsFamilies='none'` would load nothing, but the deck-exclusion sweep requires
    # an explicit list. So: name the two neutral lithium families the M7 deck used. Neither
    # produces a cation, so the charge network below is still library-only.
    kineticsFamilies=[
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

# --- neutral lithium feed: atomic Li, net charge 0, no cation seeded ---
species(
    label='Li',
    reactive=True,
    structure=adjacencyList("""
multiplicity 2
1 Li u1 p0 c0
"""),
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
    electronTemperature=(11604.5, 'K'),          # ~1 eV, inside the Voronov fit's range
    electronDensity=(1e17, 'm^-3'),
    initialMoleFractions={
        'Li': 0.05,
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
    # I-126. Chemkin and Cantera can both express these two rate laws -- lossily,
    # as the modified-Arrhenius reduction along T = Te, and each says so in the
    # file it writes. ReactionMechanismSimulator cannot express them at ALL:
    # `rmgpy/yaml_rms.py:252` has no branch for VoronovEIArrhenius or
    # BadnellRRArrhenius and raises
    #   ValueError: Object of type <class '...VoronovEIArrhenius'> does not have a
    #   defined conversion to ReactionMechanismSimulator format
    # That is a kinetics-coverage gap in the RMS writer, not an electron-placement
    # one -- it fires on the rate law before the equation is ever built -- and it
    # is a separate ticket. Until it lands, a deck whose whole chemistry is these
    # two rate laws must not ask for an RMS file; the writer is enabled by default
    # and would kill the run on the first save.
    generateRMSYAML=True,
    # Turned ON so the deck exercises the second writer that shares
    # `rmgpy.electron_balance` with Chemkin, rather than leaving the Cantera export
    # path proven only by unit test.
    generateCanteraYAML2=True,
)
