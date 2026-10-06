# 1. Database
database(
    thermoLibraries=['primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo', 'CBS_QB3_1dHR'],
    reactionLibraries=[],
    transportLibraries=['OneDMinN2', 'PrimaryTransportLibrary'],
    seedMechanisms=[],
    kineticsDepositories='default',
    kineticsFamilies='default',
    kineticsEstimator='rate rules',
)

# 2. Species Definitions
# N2 is the sole inert bath-gas collider and polymer-phase solvent. Keeping it
# nonreactive prevents an oxygen/nitrogen reaction library from entering the model.
species(
    label='N2',
    reactive=False,
    structure=SMILES("N#N")
)

# 3. Polymer Definition
polymer(
    label='PS',
    monomer='[CH2][CH]c1ccccc1',
    end_groups=['[CH3]', '[H]'], 
    cutoff=3,
    Mn=5000.0,
    Mw=6000.0,
    # The example charge is 50 g; at Mn=5000 g mol^-1 this is 0.01 mol of chains.
    initial_mass=0.05,
    # Release the real gas-phase monomer for every chain-end depropagation event.
    monomer_product='C=Cc1ccccc1',
    # Modified-Arrhenius fits use the complete 300--1000 K k(T) tables in the
    # PS kMC event set (RMG-database cd86d4e1c). A [s^-1] except termination
    # [m^3 mol^-1 s^-1], Ea [J/mol]. The maximum pointwise fit errors are
    # stated below so the deck preserves the table-to-kernel approximation.
    radical_qssa_unzip={
        # Backbone C--C homolysis, the reverse R_Recombination record
        # evt_1955f293aa8391769efbf18f67bf0c5116f3b33d4aeebd26b07183c835855b8e;
        # per-breakable-bond fit, max |k_fit/k_table - 1| = 2.19%.
        'initiation': {'A': 3.04130412e23, 'n': -2.45832701, 'Ea': 337516.492},
        # Chain-end beta-scission to styrene, reverse R_Addition_MultipleBond
        # evt_02fdc14b2ebcbbf62537b921e396c2db9269ef13309e2659ed0401fe4be692d2;
        # PLP-SEC propagation / RMG-thermo reverse, max fit error = 1.57%.
        'depropagation': {'A': 3.17364888e18, 'n': -1.69158135, 'Ea': 114749.955},
        # Benzylic-end termination is the sum of R_Recombination
        # evt_14528439785d2886805c6a627e725fa70e1d5f1c50f648e9fd109d8c6c5519c4
        # and Disproportionation
        # evt_ed2816b0a34741abbb5d6f63756bada20175d4007b8485bbfa926c33e8dd2698;
        # max fit error = 0.87%.
        'termination': {'A': 5.84674717e7, 'n': -0.555949338, 'Ea': 16783.0069},
        # H_Abstraction transfer from a benzylic end radical to pristine PS:
        # sum the 14 event tables below, divide by the three-unit proxy, then
        # premultiply by the reproduced initial repeat-unit concentration
        # 6484.79078956 mol/m^3 (the solver input is pseudo-first-order).
        # Max fit error = 6.27%. Records (all family H_Abstraction):
        # evt_2370b880681b61e289c8cb5146c12cd78f20d7f29aca68674e4762acec33d8bb,
        # evt_32b9290387171e0f0a06c398c0dba52bd56e4a0c684fec216e1e55ee30277e08,
        # evt_5b7f5989374e3fa65b383e8e9f46489c8c9aff5552231afc1c7058626af5479b,
        # evt_5f9fae411aca0036d4d9eaa8b423ddad3c83da9e145eec7c0c1cc94564db5ff4,
        # evt_77e872b8b5b9d29a02141768cfe865731cb7c1d1948d8222dda67893d41f3222,
        # evt_79f7b49e5c0a44b8da838dfe06f930561d6fcfdfe91f3e420629272cc4bd8eb5,
        # evt_81aea80470ae52bf633db0e2ff04daca3ea1b6a77c04bdb2f58ed4f4328ec7ff,
        # evt_89dfe1034ff2b44eb56b3cedf9cb129bb7f0391f637baba0f92bbf110d695157,
        # evt_8e052b9404ceaeb0a54c1cabc070251a7e8c56716aa4f4f2098d2449b6d1c09f,
        # evt_d112c9f8c51e6c1406c8d68e167abf4d8c7ad1281a7e0612520cad191b5c1926,
        # evt_d694ac1e75337318de35f8893e7ce18f2ac57fd814f60ccf88bcda0c1ad7bc69,
        # evt_eebf4bf7bb1c46c3b632ee8072ccbc2f1a1fe085d0559d9c6f3f4ccfc7fe339b,
        # evt_f47c72ad826f8cee629a3182530c9f5b5de9365ca2f6199f58d8928e02e9af27,
        # evt_fb3dda5907375fe173c5ee3f4f952d99e5545826b09aa456ab3d472899d459bc.
        'transfer': {'A': 6.80127073e-11, 'n': 5.72127748, 'Ea': 26917.9991},
        'efficiency': 1.0,
        'monomer_yield': 1.0,
        'basis': 'backbone_bonds_mu1_minus_mu0',
    },
)

# 4. Polymer Phase Definition
pp = polymer_phase(
    label='Polystyrene_Melt',
    species=['PS', 'N2'], 
    solvent='N2',
    density=(1050.0, 'kg/m^3'),
)

# 5. Hybrid Polymer Reactor
hybridPolymerReactor(
    temperature=(1000.0, 'K'),
    pressure=(1.0, 'bar'),
    initialMoles={
        'N2': 0.99,
        # Match the polymer declaration: 0.05 kg / 5000 g mol^-1 = 0.01 mol of chains.
        # Keeping this pool equal to initial_mass/Mn prevents the pool-consistency warning.
        'PS': 0.01,
    },
    polymerPhase=pp,
    terminationTime=(0.1, 's'),
    sensitivity=None,
    constant_gas_volume=False,
)

# 6. Model and Simulator Settings
model(
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
    filterReactions=False,
    filterThreshold=100000000.0,
    maxNumObjsPerIter=1,
    terminateAtMaxObjects=False,
)

simulator(atol=1e-16, rtol=1e-08, sens_atol=1e-06, sens_rtol=0.0001)

# 7. PDep
pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(2.0, 'kJ/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300, 2500, 'K', 10),
    pressures=(0.1, 100, 'bar', 10),
    interpolation=('Chebyshev', 6, 4),
    maximumAtoms=16,
)

options(
    name='Seed',
    generateSeedEachIteration=True,
    saveSeedToDatabase=False,
    units='si',
    generateOutputHTML=True,
    generatePlots=False,
    saveSimulationProfiles=True,
    verboseComments=False,
    saveEdgeSpecies=True,
    keepIrreversible=False,
    trimolecularProductReversible=True,
    wallTime='00:00:00:00',
    saveSeedModulus=-1,
)

generatedSpeciesConstraints(
    allowed=['input species', 'seed mechanisms', 'reaction libraries'],
    maximumCarbonAtoms=8,
    maximumOxygenAtoms=4,
    maximumNitrogenAtoms=0,
    maximumSiliconAtoms=0,
    maximumSulfurAtoms=0,
    # Styrene (the intended PS monomer product) has eight heavy atoms; five would
    # reject it as a generated species even though maximumCarbonAtoms already allows C8.
    maximumHeavyAtoms=8,
    maximumRadicalElectrons=2,
    maximumSingletCarbenes=1,
    maximumCarbeneRadicals=0,
    allowSingletO2=True,
)
