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
# N2 is the sole inert bath-gas collider. It remains in initialMoles below,
# but is deliberately absent from polymer_phase.species so its gas mass does
# not inflate the condensed volume used by polymer moment concentrations.
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
    # Modified-Arrhenius fits use all 29 points from 300--1000 K in compact
    # kmc_event_set/0.2 artifact
    # dfd3c856d265e5497c9b1a018ee7e330ca89cc4f96b8e5bd695102199c201860.
    # A is [s^-1] except termination [m^3 mol^-1 s^-1], and Ea is [J/mol].
    # No product, mass-loss, MWD, or other pyrolysis data enter these fits.
    radical_qssa_unzip={
        # Backbone C--C homolysis, reverse family R_Recombination, template
        # Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_
        # Ext-2R-R_N-3R!H->O_Ext-1C-R, record
        # evt_9fbefbcf21568b2e60bab8cb985e0e665dec1790d0f98a6cd410947900f3f104.
        # The per-bond table is multiplied by two because a styrene repeat has
        # two backbone C--C bonds and the solver counts mu1-mu0 sites. Artifact
        # sha256 dfd3c856d265e5497c9b1a018ee7e330ca89cc4f96b8e5bd695102199c201860;
        # maximum pointwise fit error 2.188%.
        'initiation': {'A': 6.082608239264e23, 'n': -2.45832700751, 'Ea': 320216.873028},
        # Chain-end beta-scission to styrene, reversed family
        # R_Addition_MultipleBond, template Cds-HH_Cds-CbH;CsJ-CbCsH,
        # record evt_1eab9574e11a21df23ed7c522c0753bca5e1a61eb49c98dffb213b52d6364a91.
        # Artifact sha256 dfd3c856d265e5497c9b1a018ee7e330ca89cc4f96b8e5bd695102199c201860;
        # maximum pointwise fit error 1.574%.
        'depropagation': {'A': 3.173648881804e18, 'n': -1.69158135443, 'Ea': 114750.084903},
        # Benzylic-end termination sums family R_Recombination, template
        # Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_
        # Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R, record
        # evt_c70402aba338c6984a32b76cd423e6b88c5952793ecb4fd6227d7e5a225caf7d,
        # and family Disproportionation, template Root_Ext-1R!H-R_N-4R->O_N-
        # Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClCl
        # FFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_
        # 2R!H->C, record
        # evt_79d839a45609b268524b818c1fd18c137e0046419593ac0ff74d3bcb594c045d.
        # Artifact sha256 dfd3c856d265e5497c9b1a018ee7e330ca89cc4f96b8e5bd695102199c201860;
        # maximum pointwise fit error 0.077%.
        'termination': {'A': 2.655663993173e8, 'n': -0.553042197084, 'Ea': 728.80426424},
        # H_Abstraction transfer from a benzylic end radical to pristine PS:
        # sum the 14 tables, divide by the three-repeat proxy, and premultiply
        # by mu1/V_poly = 10081.7203363 mol/m^3. Family H_Abstraction;
        # templates C/H/Cs2Cb;C_rad/H/CbCs, C/H2/NonDeC;C_rad/H/CbCs,
        # C/H3/Cs;C_rad/H/CbCs, and Cb_H;C_rad/H/CbCs. Records:
        # evt_16564f15da56c985d485d664867b311c42b922d5ece4fd864d8b25d69c9cc8d0,
        # evt_3a50995474b576fa24b539865bf3f7046d10d155b7c582ecb29644998bc15eb5,
        # evt_4bb2b692a5c44bb7e119926d01e21ae4f45bc2e8ee9a9ca4b8a2471a402c4418,
        # evt_52e491893f3f546aebf9ef3983320e8a29519977cc116e3bfa05194badb0ec36,
        # evt_6a3f8c626776f4cb009c484b017a644082fc36ae2673fc009e9e0608c1486448,
        # evt_74693b3ac0fdfef9e854bdaedc5664e41f9da55186386cd72689d5e5ee49f015,
        # evt_83b6e5d3a986c9dbd44049c4feb6bfa2e813036f44c3d2818abce2d57f2d636d,
        # evt_89660c68f37b6d0fd55f1d69770dc3c9829eec64391cf67a775c845fcbcf13e0,
        # evt_8b95b75c6f381ec3cce185ab1632d26dcfee7e4885dfcb76085aef404ac8e7b4,
        # evt_951195859cacccbc2bea426da625ed3c600a3bb2074569411e489ad3ed118d04,
        # evt_970fc40a943edee87d86564c8f8bb13a23c73c9ee6415a3dfd2f3d09f62267e6,
        # evt_b534c2dad4b2c9ed4f65343c0f2b76d38a8d8c50a5a02ac1ce582da349678247,
        # evt_c3e583adb9e4973ec543f3fa01d0a455025852877b6f80756bdf8f99d630ea8b,
        # evt_d8c1e1c2eda95f0d6f227812851eaa1b7ef5ae107f7bc6d2b00a7b7bca3b7962.
        # Artifact sha256 dfd3c856d265e5497c9b1a018ee7e330ca89cc4f96b8e5bd695102199c201860;
        # maximum pointwise fit error 6.265%.
        'transfer': {'A': 1.057374272215e-10, 'n': 5.72127747789, 'Ea': 26918.0295221},
        'efficiency': 1.0,
        'monomer_yield': 1.0,
        'basis': 'backbone_bonds_mu1_minus_mu0',
    },
)

# 4. Polymer Phase Definition
pp = polymer_phase(
    label='Polystyrene_Melt',
    species=['PS'],
    solvent='N2',
    # Only PS contributes to condensed volume. The solvent label supplies
    # phase metadata; N2 stays a gas species and bath-gas collider.
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
