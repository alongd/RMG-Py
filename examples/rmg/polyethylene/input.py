# 1. Database
database(
    thermoLibraries=[
        'primaryThermoLibrary',
        'thermo_DFT_CCSDTF12_BAC',
        'DFT_QCI_thermo',
        'CBS_QB3_1dHR',
    ],
    reactionLibraries=[],
    transportLibraries=['OneDMinN2', 'PrimaryTransportLibrary'],
    seedMechanisms=[],
    kineticsDepositories='default',
    kineticsFamilies='default',
    kineticsEstimator='rate rules',
)

# 2. Species Definitions
# N2 is the sole inert bath-gas collider. It remains in initialMoles below,
# but is deliberately absent from polymer_phase.species so its gas mass does
# not inflate the condensed volume used by polymer moment concentrations.
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
    # over 300--1000 K from the snapshot declared as RMG-database polymer
    # commit cd86d4e1c (the snapshot has no .git metadata). Each surrogate
    # is processed through CoreEdgeReactionModel.make_new_reaction, including
    # production thermo, own-reverse-family direction choice, and barrier
    # correction. A is
    # [s^-1] except termination [m^3 mol^-1 s^-1], and Ea is [J/mol]. No
    # product, mass-loss, MWD, or other pyrolysis data enter these fits.
    radical_qssa_unzip={
        # Backbone C--C homolysis: thermodynamic reverse
        # of family R_Recombination for n-C3H7 + n-C3H7 -> n-hexane (the
        # central n-hexane bond). The single-bond rate is multiplied by two:
        # a PE ethylene repeat contributes two backbone C--C bonds, whereas
        # the solver applies initiation to mu1-mu0 repeat-bond units. Thus it
        # uses 2(mu1-mu0) instead of the finite-chain count 2mu1-mu0:
        # initially one bond, or 0.2813%, fewer. Source: ArrheniusBM rate-rule
        # node selected and converted by the production path,
        # Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_
        # Ext-3R!H-R_N-Sp-3R!H=2R plus production thermo. Max fit error 1.640%.
        'initiation': {
            'A': 9.231130302962257e26,
            'n': -2.955777268894243,
            'Ea': 373487.4794598093,
        },
        # Chain-end beta-scission: 1-hexyl -> ethylene + 1-butyl, the
        # thermodynamic reverse of family R_Addition_MultipleBond. Source:
        # training reaction 2905, exact rule [Cds-HH_Cds-HH;CsJ-CsHH], plus
        # production thermo. One ethylene is released per event. Max fit error 0.909%.
        'depropagation': {
            'A': 2.611138149150429e9,
            'n': 1.1159009347437787,
            'Ea': 124806.10852008159,
        },
        # Primary chain-end termination is the sum of (1) R_Recombination,
        # 2 1-hexyl -> n-dodecane, from exact training reaction 156,
        # and (2) Disproportionation, 2 1-hexyl -> n-hexane + 1-hexene, from
        # the database's matching ArrheniusBM Root_Ext-1R!H...Ext-6C-R tree
        # node converted with reaction dH(298). The solver requires Ea >= 0,
        # so the negligible -24.9 J/mol unconstrained result is refitted at
        # the Ea=0 boundary. The summed fit has max error 0.070%.
        'termination': {'A': 3.6152759909271235e6, 'n': 0.14991913312996882, 'Ea': 0.0},
        # The solver has one pseudo-first-order transfer sink, so this fit
        # sums both RMG-supported routes: (a) family H_Abstraction for 1-hexyl
        # + n-octane -> n-hexane + secondary 2/3/4-octyl, average rate rule
        # [C/H2/NonDeC;C_rad/H2/Cs], divided by three interior ethylene-repeat
        # equivalents and multiplied by 33864.2717058 mol/m^3 initial PE
        # repeat units;
        # and (b) all six family intra_H_migration paths from 1-octyl to
        # secondary octyl, sourced by training reactions 106, 108, 110, 112,
        # 114, and 116, with model-generation source priority and every
        # production path. Density 950 kg/m^3 and exact RMG repeat MW
        # 28.0531649478 g/mol set the concentration conversion. The resulting
        # pseudo-first-order coefficient is fixed at this initial
        # concentration as mu1 falls. Max combined-fit error 3.073%.
        'transfer': {
            'A': 5.573253173648898e-9,
            'n': 5.748594368061862,
            'Ea': 26820.365227948805,
        },
        'efficiency': 1.0,
        'monomer_yield': 1.0,
        'basis': 'backbone_bonds_mu1_minus_mu0',
    },
)

# 4. Polymer Phase Definition
pp = polymer_phase(
    label='Polyethylene_Melt',
    species=['PE'],
    solvent='N2',
    # Only PE contributes to condensed volume. The solvent label supplies
    # phase metadata; N2 stays a gas species and bath-gas collider.
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
