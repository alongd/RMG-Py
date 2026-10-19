# Toy linear CH2-CH2 deck (the cluster PE runs' deck); its declaration digests are
# pinned by polymerSpecV2Test.  Structure and a declared toy MET block only.
database(
    thermoLibraries=['primaryThermoLibrary'],
    reactionLibraries=[],
    seedMechanisms=[],
    kineticsDepositories='default',
    kineticsFamilies='default',
    kineticsEstimator='rate rules',
)
species(label='N2', reactive=False, structure=SMILES('N#N'))
polymer('toy linear polymer (CH2-CH2)', linear_homopolymer_spec={'schema_version': 1, 'label': 'toy linear polymer (CH2-CH2)', 'repeat_unit': {'atoms': [{'key': 'c1', 'element': 'C', 'hydrogens': 2, 'electronic_state': {'radical_electrons': 0}}, {'key': 'c2', 'element': 'C', 'hydrogens': 2, 'electronic_state': {'radical_electrons': 0}}], 'bonds': [{'atoms': ['c1', 'c2'], 'order': 1}], 'ports': {'head': {'atom': 'c1', 'bond_order': 1}, 'tail': {'atom': 'c2', 'bond_order': 1}}, 'electronic_state': {'multiplicity': 1}}, 'end_groups': {'head': {'atoms': [{'key': 'h', 'element': 'H', 'hydrogens': 0, 'electronic_state': {'radical_electrons': 0}}], 'bonds': [], 'ports': {'cap': {'atom': 'h', 'bond_order': 1}}, 'electronic_state': {'multiplicity': 1}}, 'tail': {'atoms': [{'key': 'h', 'element': 'H', 'hydrogens': 0, 'electronic_state': {'radical_electrons': 0}}], 'bonds': [], 'ports': {'cap': {'atom': 'h', 'bond_order': 1}}, 'electronic_state': {'multiplicity': 1}}}, 'small_species': {}, 'initial_ensemble': {'dp_counts': {30: 100}}, 'initial_mobile_counts': {}, 'density_basis': 'dry_polymer', 'initial_mass_kg': 1.0, 'density_kg_m3': 850.0, 'locality_policy': {'radius': 2, 'bounds': {'elements': ['C', 'H'], 'isotopes': [-1], 'formal_charges': [0], 'radical_electrons': [0, 1], 'lone_pairs': [0], 'bond_orders': [1, 1.5, 2, 3], 'spin_multiplicities': [1, 2, 3], 'maximum_heavy_atom_degree': 4, 'maximum_explicit_hydrogens_per_heavy_atom': 4, 'maximum_total_hydrogens_on_free_small_species': 12, 'maximum_family_arity': 2, 'maximum_root_tuple_size': 4, 'maximum_free_small_species_heavy_atoms': 8, 'attachment_classes': [], 'topology_metadata_schema': {}}}, 'family_manifest': {}, 'library_manifest': {}, 'transport_model': {'met': {'schema': 'kmc-met-transport/1', 'segment_diffusivity': {'law': 'arrhenius', 'D_ref_m2_s': 1e-17, 'T_ref_K': 700.0, 'Ea_J_mol': 0.0}, 'chain_scaling': {'law': 'rouse'}, 'mean_square_end_to_end_per_molar_mass_m2_mol_per_g': 1e-20, 'contact_radius_m': 5e-10, 'spin_factor': 1.0}}, 'diagnostics': {}}, cutoff=4)
simpleReactor(
    temperature=(700, 'K'),
    pressure=(1.0, 'bar'),
    initialMoleFractions={'N2': 1.0},
    terminationTime=(1.0, 's'),
)
model(toleranceKeepInEdge=0.0, toleranceMoveToCore=0.1, toleranceInterruptSimulation=0.1)
options(units='si', generateOutputHTML=False)
