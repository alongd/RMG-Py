# EPDM (ethylene-propylene-diene monomer) pyrolysis -- copolymer composition example.
#
# EPDM is a RANDOM TERPOLYMER, so it is declared with `monomers=[...]` (a
# composition) rather than `monomer=` (a single repeat unit). The composition
# sets the composition-weighted repeat mass <M> = sum_i f_i M_i -- which is what
# the moments, Mn/Mw and every condensed-mass term are built on -- and the pool
# builds one dyad proxy per comonomer pair, so the mixed-neighbour backbone
# bonds (ethylene--propylene, propylene--ENB, ethylene--ENB) exist in real
# graphs and can generate chemistry. That matters here specifically: EPDM's weak
# link is the allylic C--C beside an incorporated diene unit, and that bond only
# exists where a diene unit meets a backbone unit.
#
# SCOPE LIMIT (2026-07-27, GRAPHS ONLY): this deck exercises the composition
# mass physics and the dyad graphs. The dyad ROWS are currently deferred by the
# adjudicated refusal policy (stamp-but-keep, flux zeroed), so the run completes
# with the dyad proxies present as core SPECIES but no core reaction touching
# them. Treat the mixed-neighbour chemistry as not-yet-contributing until that
# is revisited against TGA scoring.
#
# COMPOSITION OF RECORD -- Vistalon 5601 (pinned 2026-07-28):
#   Perejon, Sanchez-Jimenez, Gil-Gonzalez, Perez-Maqueda and Criado,
#   "Pyrolysis kinetics of ethylene-propylene (EPM) and ethylene-propylene-diene
#   (EPDM)", Polym Degrad Stab 98 (2013) 1571-1577, Table 1
#   (DOI 10.1016/j.polymdegradstab.2013.06.029).
#
#   Vistalon 5601 (ExxonMobil): 69.0 wt% ethylene, ENB 5.0 wt%, balance
#   propylene (26.0 wt%). Chosen among that table's four ENB grades because it
#   carries the ENB diene modeled here (not DCPD) at a typical 5 wt%, with the
#   lowest ethylene content of the ENB set -- the least crystalline, most
#   elastomeric grade, which is the character of SRM insulation gum. Swapping to
#   Vistalon 7001 (73.0/5.0), 7800 (79.0/6.0) or EPDM 70C2 (70.0/ENB 4.6) is a
#   three-number edit here and a one-line edit in the CKMG reference.
#
#   Mole fractions, converted from those weight percents:
#       x_E   = (69/28.053)  / 3.1191 = 0.7886
#       x_P   = (26/42.080)  / 3.1191 = 0.1981
#       x_ENB = (5/120.191)  / 3.1191 = 0.0133
#   giving a composition-weighted repeat mass of 32.061 g/mol.
#
#   The same paper supplies this sample's measured degradation temperatures at
#   5 K/min (T_alpha=0.1 = 705.1 K, T_alpha=0.5 = 720.2 K) and its fitted
#   kinetic triplet, which are the validation reference of record downstream.
#
# ASSUMPTIONS OF RECORD -- MATERIAL STATE:
#   Neat, UNCURED, UNFILLED gum under an INERT atmosphere. Real SRM insulation
#   is vulcanized and filled (carbon black, silica, aramid); crosslinks and
#   fillers are deliberately absent here and are a known, stated bias against
#   any measured COMPOSITE TGA curve, which must be normalized to resin mass
#   before comparison.
#
# The ENB repeat unit is the addition-polymerized unit: the strained
# 2-norbornene ring double bond is consumed by polymerization (it becomes the
# two backbone attachment sites), while the exocyclic ethylidene C=C survives as
# the pendant unsaturation -- which is exactly the site that makes EPDM
# vulcanizable and thermally weaker than EPM.

# 1. Database
database(
    thermoLibraries=['primaryThermoLibrary', 'BurkeH2O2', 'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo', 'CBS_QB3_1dHR'],
    reactionLibraries=['primaryH2O2'],
    transportLibraries=['OneDMinN2', 'PrimaryTransportLibrary', 'GRI-Mech'],
    seedMechanisms=[],
    kineticsDepositories='default',
    kineticsFamilies='default',
    kineticsEstimator='rate rules',
)

# 2. Species Definitions
# N2 is the polymer-phase solvent and third-body collider; it must be reactive
# so it picks up thermo from primaryThermoLibrary, while maximumNitrogenAtoms=0
# keeps it chemically inert.
species(
    label='N2',
    reactive=True,
    structure=SMILES("N#N")
)

# Inert bath gas: pressure dependence requires at least one nonreacting species.
species(
    label='Ar',
    reactive=False,
    structure=SMILES("[Ar]")
)

# 3. Polymer Definition -- the copolymer composition
polymer(
    label='EPDM',
    monomers=[
        # Ethylene backbone unit (C2H4).
        {'monomer': '[CH2][CH2]', 'fraction': 0.7886, 'monomer_product': 'C=C'},
        # Propylene backbone unit (C3H6); the tertiary C-H is the weak backbone
        # site in the diene-free stretches.
        {'monomer': '[CH2][CH](C)', 'fraction': 0.1981, 'monomer_product': 'C=CC'},
        # 5-ethylidene-2-norbornene (ENB) unit (C9H12), polymerized through the
        # ring double bond, pendant ethylidene retained.
        {'monomer': 'CC=C1CC2[CH][CH]C1C2', 'fraction': 0.0133},
    ],
    end_groups=['[CH3]', '[H]'],
    cutoff=3,
    Mn=5000.0,
    Mw=10000.0,
    initial_mass=1.0,
)

# 4. Polymer Phase Definition
pp = polymer_phase(
    label='EPDM_Melt',
    species=['EPDM', 'N2'],
    solvent='N2',
    density=(860.0, 'kg/m^3'),  # typical EPDM gum density
)

# 5. Hybrid Polymer Reactor
hybridPolymerReactor(
    temperature=(1000.0, 'K'),
    pressure=(1.0, 'bar'),
    initialMoles={
        'N2': 0.99,
        'EPDM': 0.01,
        'Ar': 1e-10,  # bath-gas reference for pressure dependence
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

# Two constraint tiers, and the split matters here more than for a homopolymer.
# The gas tier below bounds the volatile fragments. The dyad proxies are much
# bigger than any of them -- ENB is a C9 bicyclic unit, so the diene--diene dyad
# proxy is C28 -- and bounding them with the gas tier would either refuse the
# pool's own proxies at initialization or force a gas bound so loose that RMG
# would happily generate C28 volatiles. They are routed to the polymer tier
# instead, by heavy-atom count: the SMALLEST dyad proxy here (ethylene-propylene,
# C8H18) has 8 heavy atoms, so polymerSizeThreshold=8 puts every dyad proxy on
# the polymer tier while leaving genuine pyrolysis fragments on the gas tier.
generatedSpeciesConstraints(
    allowed=['input species', 'seed mechanisms', 'reaction libraries'],
    maximumCarbonAtoms=7,
    maximumOxygenAtoms=0,
    maximumNitrogenAtoms=0,
    maximumSiliconAtoms=0,
    maximumSulfurAtoms=0,
    maximumHeavyAtoms=7,
    maximumRadicalElectrons=2,
    maximumSingletCarbenes=1,
    maximumCarbeneRadicals=0,
    allowSingletO2=False,
)

generatePolymerConstraints(
    maximumCarbonAtoms=30,
    maximumHeavyAtoms=30,
    maximumRadicalElectrons=2,
    polymerSizeThreshold=8,
)
