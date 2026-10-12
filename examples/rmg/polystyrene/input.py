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
    # f1fa450b8627d55559b0eaa3871771c3c36de9a6fbc024a54ce654b6e732b170.
    # A is [s^-1] except termination [m^3 mol^-1 s^-1], and Ea is [J/mol].
    # No product, mass-loss, MWD, or other pyrolysis data enter these fits.
    radical_qssa_unzip={
        # Backbone C--C homolysis, reverse family R_Recombination, template
        # Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_
        # Ext-2R-R_N-3R!H->O_Ext-1C-R, record
        # evt_a09358573398f9407aa264637d54523fb6fc6451e661545044acda1c25ebdc92.
        # The per-bond table is multiplied by two because a styrene repeat has
        # two backbone C--C bonds and the solver counts mu1-mu0 sites. Artifact
        # sha256 f1fa450b8627d55559b0eaa3871771c3c36de9a6fbc024a54ce654b6e732b170;
        # maximum pointwise fit error 2.188%.
        'initiation': {'A': 6.082608239264e23, 'n': -2.45832700751, 'Ea': 320216.873028},
        # Chain-end beta-scission to styrene, reversed family
        # R_Addition_MultipleBond, template Cds-HH_Cds-CbH;CsJ-CbCsH,
        # record evt_918df20362275be38908c300a8dd81c91247366be8350cf537a0aa6d4be0648d.
        # Artifact sha256 f1fa450b8627d55559b0eaa3871771c3c36de9a6fbc024a54ce654b6e732b170;
        # maximum pointwise fit error 1.574%.
        'depropagation': {'A': 3.173648881804e18, 'n': -1.69158135443, 'Ea': 114750.084903},
        # Benzylic-end termination sums family R_Recombination, template
        # Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_
        # Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R, record
        # evt_3d1bf01d9fc851afb1e0c49315cc00143f47d4a32056941963ee8fa6680f872b,
        # and family Disproportionation, template Root_Ext-1R!H-R_N-4R->O_N-
        # Sp-5R!H=1R!H_Ext-4CHNS-R_N-6R!H->S_4CHNS->C_N-Sp-6BrBrBrCCCClClCl
        # FFFIIINNNOOOPPPSiSiSi#4C_6BrCClFINOPSi->C_N-1R!H-inRing_Ext-4C-R_
        # 2R!H->C, record
        # evt_3d96955a82837d37aa643389552079339d622160884ff39ca83ab304d79bc293.
        # Artifact sha256 f1fa450b8627d55559b0eaa3871771c3c36de9a6fbc024a54ce654b6e732b170;
        # maximum pointwise fit error 0.077%.
        'termination': {'A': 2.655663993173e8, 'n': -0.553042197084, 'Ea': 728.80426424},
        # H_Abstraction transfer from a benzylic end radical to pristine PS:
        # sum the 14 tables, divide by the three-repeat proxy, and premultiply
        # by mu1/V_poly = 10081.7203363 mol/m^3. Family H_Abstraction;
        # templates C/H/Cs2Cb;C_rad/H/CbCs, C/H2/NonDeC;C_rad/H/CbCs,
        # C/H3/Cs;C_rad/H/CbCs, and Cb_H;C_rad/H/CbCs. Records:
        # evt_122dcd3626bb9abd15b9f495953c0b141a07fc6598820357a158a485c22ceb10,
        # evt_158ab84b0b6f8777f87b92b9135d7551e715d898a54a863d7c8483d4c03645cf,
        # evt_29a619474007b769df4c189a4876b76e19a7150318935a65ac778fa8d020bde8,
        # evt_404445e9537996a29812ecb02fc644bf8ef842980872b0f9f63c910253780196,
        # evt_58903aa4397009bcb6bfd5bd956ec3cac91328ac649a53d8fa0f6d25bcb94f1b,
        # evt_5eef2bd0e85f39dd94a5d600d6a32aee8316ed242f70d2a327a3e9776df1744c,
        # evt_7b458ee57b54d3db41ad87c81e47c0ffcbdf9e6136c78bbef7812e127341408c,
        # evt_7c50c9592ba91d04b42d86cae5a52b6767d8fb14224e2616733a44e1d4a1bd9b,
        # evt_8f4c40142d4e2b2fef8719da7fcebbebf6e7f1e247bb5e467f0f944acee7c685,
        # evt_a06b983cec1bce64d2b829b9b1d323b13abffb9f73823d5ffb61b2ca8c27c2bb,
        # evt_a7a7d798656380b06b55e299fc8807b675760e419ac9201c45db0bfb88d0a287,
        # evt_c84b8a0c273076cef2e9fde028afe0609b7c67075d6e3366e19112e7786c100d,
        # evt_cb546649a572eb972dc945b263413d1bd14340ede33674362c44ee04894f5cd6,
        # evt_e222edf0d1c27595f69d1149622c9cffd365521c9bc6480491144ed5b5cbaf79.
        # Artifact sha256 f1fa450b8627d55559b0eaa3871771c3c36de9a6fbc024a54ce654b6e732b170;
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
