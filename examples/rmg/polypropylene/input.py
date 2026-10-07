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
    label='PP',
    monomer='[CH2][CH]C',
    end_groups=['[H]', '[H]'],
    cutoff=3,
    Mn=5000.0,
    Mw=6000.0,
    # At Mn=5000 g mol^-1, the 50 g charge is 0.01 mol of chains.
    initial_mass=0.05,
    # Release real gas-phase propylene for each chain-end event.
    monomer_product='C=CC',
    # Modified-Arrhenius fits use 15 RMG points at 50 K intervals
    # over 300--1000 K from the snapshot declared as RMG-database polymer
    # commit cd86d4e1c (the snapshot has no .git metadata). Each surrogate
    # is processed through CoreEdgeReactionModel.make_new_reaction, including
    # production thermo, own-reverse-family direction choice, and barrier
    # correction. A is
    # [s^-1] except termination [m^3 mol^-1 s^-1], and Ea is [J/mol]. No
    # product, mass-loss, MWD, or other pyrolysis data enter these fits.
    radical_qssa_unzip={
        # Backbone secondary--tertiary C--C homolysis: reverse
        # of family R_Recombination for isobutyl + isopropyl ->
        # 2,4-dimethylpentane. The single-bond rate is doubled: a PP
        # repeat contributes two backbone C--C bonds, while the solver
        # applies initiation to mu1-mu0 repeat-bond units. Thus it uses
        # 2(mu1-mu0) instead of the finite-chain count 2mu1-mu0: initially
        # one bond, or 0.4226%, fewer. Source: ArrheniusBM rate-rule entry
        # 111, selected and converted by the production path, degeneracy 1,
        # plus production thermo. Max fit error 1.939%.
        'initiation': {
            'A': 9.026483351957035e26,
            'n': -2.636125594368046,
            'Ea': 367560.03286372335,
        },
        # Chain-end beta-scission: 4-methyl-2-pentyl -> propylene +
        # isopropyl, the reverse of R_Addition_MultipleBond. Source:
        # exact training reaction 239, propene_1 + C3H7-2 <=> C6H13-2,
        # degeneracy 1, plus RMG thermo. One propylene per event.
        # Max fit error 0.673%.
        'depropagation': {
            'A': 8.798203848271057e11,
            'n': 0.471525701306732,
            'Ea': 105603.71475817061,
        },
        # Secondary chain-end termination sums (1) R_Recombination
        # for two 4-methyl-2-pentyl radicals, ArrheniusBM entry 160
        # (degeneracy 0.5), and (2) both Disproportionation products,
        # entries 226 and 20 (degeneracies 2 and 3). Every BM term uses
        # its own reaction dH(298). The solver requires Ea >= 0, so the
        # -174.3 J/mol unconstrained fit is refitted
        # at the Ea=0 boundary. The summed fit has max error 0.511%.
        'termination': {
            'A': 1.2015531789581546e8,
            'n': -0.350299014152254,
            'Ea': 0.0,
        },
        # The solver has one pseudo-first-order transfer sink. This fit
        # sums PP's tertiary-C--H routes: (a) family H_Abstraction for a
        # 4-methyl-2-pentyl + 2,4,6-trimethylheptane, retaining two
        # tertiary products (path degeneracies 1 and 2), divided by
        # three propylene-repeat equivalents and multiplied by
        # 21387.9610774 mol/m^3 initial repeat units; and (b) the family
        # intra_H_migration R5H_CCC tertiary 1,5-shift from a trimer
        # secondary end radical (degeneracy 1). Sources are the model-
        # production-path averaged rate rules. Nominal density 900 kg/m^3
        # and exact RMG repeat MW 42.0797474217 g/mol set the conversion.
        # The resulting pseudo-first-order coefficient is fixed at this
        # initial concentration as mu1 falls. Max combined-fit error 7.756%.
        'transfer': {
            'A': 5.0957274710100526e-3,
            'n': 3.798567610993778,
            'Ea': 34881.09549034845,
        },
        'efficiency': 1.0,
        'monomer_yield': 1.0,
        'basis': 'backbone_bonds_mu1_minus_mu0',
    },
)

# 4. Polymer Phase Definition
pp = polymer_phase(
    label='Polypropylene_Melt',
    species=['PP'],
    solvent='N2',
    # Only PP contributes to condensed volume. The solvent label supplies
    # phase metadata; N2 stays a gas species and bath-gas collider.
    density=(900.0, 'kg/m^3'),
)

# 5. Hybrid Polymer Reactor
hybridPolymerReactor(
    temperature=(1000.0, 'K'),
    pressure=(1.0, 'bar'),
    initialMoles={
        'N2': 0.99,
        # Match initial_mass/Mn: 0.05 kg / 5 kg mol^-1 = 0.01 mol.
        'PP': 0.01,
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
