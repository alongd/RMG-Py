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
# nonreactive prevents nitrogen chemistry from entering the model.
species(
    label='N2',
    reactive=False,
    structure=SMILES('N#N'),
)

# 3. Polymer Definition
polymer(
    label='PE',
    monomer='[CH2][CH2]',
    end_groups=['[H]', '[H]'],
    cutoff=3,
    Mn=5000.0,
    Mw=6000.0,
    # The example charge is 50 g; at Mn=5000 g mol^-1 this is 0.01 mol of chains.
    initial_mass=0.05,
    # Release real gas-phase ethylene for each chain-end depropagation event.
    monomer_product='C=C',
    # Modified-Arrhenius fits use 15 RMG-computed points at 50 K intervals
    # over 300--1000 K from RMG-database polymer commit cd86d4e1c. A is
    # [s^-1] except termination [m^3 mol^-1 s^-1], and Ea is [J/mol]. No
    # product, mass-loss, MWD, or other pyrolysis data enter these fits.
    radical_qssa_unzip={
        # Backbone C--C homolysis: thermodynamic reverse
        # of family R_Recombination for n-C3H7 + n-C3H7 -> n-hexane (the
        # central n-hexane bond). The single-bond rate is multiplied by two:
        # a PE ethylene repeat contributes two backbone C--C bonds, whereas
        # the solver applies initiation to mu1-mu0 repeat-bond units. Source:
        # ArrheniusBM rate-rule node, converted with the reaction dH(298),
        # Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_
        # Ext-3R!H-R_N-Sp-3R!H=2R plus RMG thermo. Max fit error 1.775%.
        'initiation': {'A': 1.085018209146e27, 'n': -2.97626501263, 'Ea': 373605.265711},
        # Chain-end beta-scission: 1-hexyl -> ethylene + 1-butyl, the
        # thermodynamic reverse of family R_Addition_MultipleBond. Source:
        # training reaction 2905, exact rule [Cds-HH_Cds-HH;CsJ-CsHH], plus
        # RMG thermo. One ethylene is released per event. Max fit error 1.131%.
        'depropagation': {'A': 3.051203080361e9, 'n': 1.09630553255, 'Ea': 124921.404911},
        # Primary chain-end termination is the sum of (1) R_Recombination,
        # 2 1-hexyl -> n-dodecane, from exact training reaction 156,
        # and (2) Disproportionation, 2 1-hexyl -> n-hexane + 1-hexene, from
        # the database's matching ArrheniusBM Root_Ext-1R!H...Ext-6C-R tree
        # node converted with reaction dH(298). The solver requires Ea >= 0,
        # so the negligible -24.9 J/mol unconstrained result is refitted at
        # the Ea=0 boundary. The summed fit has max error 0.070%.
        'termination': {'A': 3.615275990927e6, 'n': 0.149919133130, 'Ea': 0.0},
        # The solver has one pseudo-first-order transfer sink, so this fit
        # sums both RMG-supported routes: (a) family H_Abstraction for 1-hexyl
        # + n-octane -> n-hexane + secondary 2/3/4-octyl, average rate rule
        # [C/H2/NonDeC;C_rad/H2/Cs], divided by three interior ethylene-repeat
        # equivalents and multiplied by 33864.2776785 mol/m^3 PE repeat units;
        # and (b) all six family intra_H_migration paths from 1-octyl to
        # secondary octyl, sourced by training reactions 106, 108, 110, 112,
        # 114, and 116, with model-generation source priority and every
        # ArrheniusEP/BM term converted using its reaction dH(298). Density
        # 950 kg/m^3 and repeat MW 28.05316 g/mol set the concentration
        # conversion. Max combined-fit error 2.829%.
        'transfer': {'A': 1.312764315150e-3, 'n': 4.11498243112, 'Ea': 33194.2576083},
        'efficiency': 1.0,
        'monomer_yield': 1.0,
        'basis': 'backbone_bonds_mu1_minus_mu0',
    },
)

# 4. Polymer Phase Definition
pp = polymer_phase(
    label='Polyethylene_Melt',
    species=['PE', 'N2'],
    solvent='N2',
    # Nominal example density; it also sets the repeat concentration used in
    # the pseudo-first-order transfer conversion documented above.
    density=(950.0, 'kg/m^3'),
)

# 5. Hybrid Polymer Reactor
hybridPolymerReactor(
    temperature=(1000.0, 'K'),
    pressure=(1.0, 'bar'),
    initialMoles={
        'N2': 0.99,
        # Match initial_mass/Mn: 0.05 kg / (5 kg mol^-1) = 0.01 mol chains.
        'PE': 0.01,
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

# 7. Pressure dependence
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
    maximumOxygenAtoms=0,
    maximumNitrogenAtoms=0,
    maximumSiliconAtoms=0,
    maximumSulfurAtoms=0,
    maximumHeavyAtoms=8,
    maximumRadicalElectrons=2,
    maximumSingletCarbenes=1,
    maximumCarbeneRadicals=0,
    allowSingletO2=True,
)
