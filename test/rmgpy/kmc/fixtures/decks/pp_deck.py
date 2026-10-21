# Polypropylene, tacticity not declared; H caps.
# Structure and a declared MET block: no rate, yield or kinetic parameter.
# MET transport: literature values, each field cited in transport_model.sources.  The density is the deck's volume basis only.
database(
    thermoLibraries=['primaryThermoLibrary'],
    reactionLibraries=[],
    seedMechanisms=[],
    kineticsDepositories='default',
    kineticsFamilies='default',
    kineticsEstimator='rate rules',
)
species(label='N2', reactive=False, structure=SMILES('N#N'))
polymer('polypropylene (CH2-CH(CH3))', linear_homopolymer_spec={'schema_version': 1,
 'label': 'polypropylene (CH2-CH(CH3))',
 'repeat_unit': {'atoms': [{'key': 'c1',
                            'element': 'C',
                            'hydrogens': 2,
                            'electronic_state': {'radical_electrons': 0}},
                           {'key': 'c2',
                            'element': 'C',
                            'hydrogens': 1,
                            'electronic_state': {'radical_electrons': 0}},
                           {'key': 'c3',
                            'element': 'C',
                            'hydrogens': 3,
                            'electronic_state': {'radical_electrons': 0}}],
                 'bonds': [{'atoms': ['c1', 'c2'], 'order': 1},
                           {'atoms': ['c2', 'c3'], 'order': 1}],
                 'ports': {'head': {'atom': 'c1', 'bond_order': 1},
                           'tail': {'atom': 'c2', 'bond_order': 1}},
                 'electronic_state': {'multiplicity': 1}},
 'end_groups': {'head': {'atoms': [{'key': 'h',
                                    'element': 'H',
                                    'hydrogens': 0,
                                    'electronic_state': {'radical_electrons': 0}}],
                         'bonds': [],
                         'ports': {'cap': {'atom': 'h', 'bond_order': 1}},
                         'electronic_state': {'multiplicity': 1}},
                'tail': {'atoms': [{'key': 'h',
                                    'element': 'H',
                                    'hydrogens': 0,
                                    'electronic_state': {'radical_electrons': 0}}],
                         'bonds': [],
                         'ports': {'cap': {'atom': 'h', 'bond_order': 1}},
                         'electronic_state': {'multiplicity': 1}}},
 'small_species': {},
 'initial_ensemble': {'dp_counts': {30: 100}},
 'initial_mobile_counts': {},
 'density_basis': 'dry_polymer',
 'initial_mass_kg': 1.0,
 'density_kg_m3': 900.0,
 'locality_policy': {'radius': 2,
                     'bounds': {'elements': ['C', 'H'],
                                'isotopes': [-1],
                                'formal_charges': [0],
                                'radical_electrons': [0, 1],
                                'lone_pairs': [0],
                                'bond_orders': [1, 1.5, 2, 3],
                                'spin_multiplicities': [1, 2, 3],
                                'maximum_heavy_atom_degree': 4,
                                'maximum_explicit_hydrogens_per_heavy_atom': 4,
                                'maximum_total_hydrogens_on_free_small_species': 12,
                                'maximum_family_arity': 2,
                                'maximum_root_tuple_size': 4,
                                'maximum_free_small_species_heavy_atoms': 8,
                                'attachment_classes': [],
                                'topology_metadata_schema': {}}},
 'family_manifest': {},
 'library_manifest': {},
 'transport_model': {'met': {'chain_scaling': {'entanglement_units': 133.3192, 'law': 'crossover'},
          'contact_radius_m': 5.447745e-10,
          'mean_square_end_to_end_per_molar_mass_m2_mol_per_g': 6.78e-21,
          'schema': 'kmc-met-transport/1',
          'segment_diffusivity': {'D_ref_m2_s': 1.859145e-10,
                                  'Ea_J_mol': 40000.0,
                                  'T_ref_K': 362.0,
                                  'law': 'arrhenius'},
          'spin_factor': 1.0},
  'sources': {'chain_scaling': ['N_e = N_eK M_K / 42.081 = 30.59 x 183.40 / 42.081 (a-PP(463)): R. Everaers, '
                                "H. A. Karimi-Varzaneh, F. Fleck, N. Hojdis, C. Svaneborg, 'Kremer-Grest "
                                'models for commodity polymer melts: linking theory, experiment, and '
                                "simulation at the Kuhn scale', Macromolecules 53, 1901-1916 (2020), "
                                'doi:10.1021/acs.macromol.9b02428 (arXiv:1808.03509v3), Table 1'],
              'contact_radius_m': ['sigma = 2 r, r the radius of a sphere of the McGowan volume of the '
                                   "radical site C3H7 (3 C, 7 H, 9 bonds): M. H. Abraham, J. C. McGowan, 'The "
                                   'use of characteristic volumes to measure cavity terms in reversed phase '
                                   "liquid chromatography', Chromatographia 23, 243-246 (1987), Table 2 (C "
                                   '16.35, H 8.71, minus 6.56 per bond, cm3/mol), as implemented in RMG-Py '
                                   'rmgpy/data/solvation.py set_mcgowan_volume'],
              'mean_square_end_to_end_per_molar_mass_m2_mol_per_g': ['<R^2>/M = 0.678 A^2 mol/g (a-PP, '
                                                                     '298-463 K rows): R. Everaers, H. A. '
                                                                     'Karimi-Varzaneh, F. Fleck, N. Hojdis, '
                                                                     "C. Svaneborg, 'Kremer-Grest models for "
                                                                     'commodity polymer melts: linking '
                                                                     'theory, experiment, and simulation at '
                                                                     "the Kuhn scale', Macromolecules 53, "
                                                                     '1901-1916 (2020), '
                                                                     'doi:10.1021/acs.macromol.9b02428 '
                                                                     '(arXiv:1808.03509v3), Table 1'],
              'segment_diffusivity': ['anchor D0(348 K) = l_K^2 M_K / (3 pi^2 tau_K M0) = 4.628e-11 m2/s from '
                                      'the Kuhn time tau_K(348 K) = 3.99 ns (a-PP, from the experimental '
                                      'entanglement time) with l_K = 11.20 A, M_K = 183.40 g/mol (a-PP(348)): '
                                      'R. Everaers, H. A. Karimi-Varzaneh, F. Fleck, N. Hojdis, C. Svaneborg, '
                                      "'Kremer-Grest models for commodity polymer melts: linking theory, "
                                      "experiment, and simulation at the Kuhn scale', Macromolecules 53, "
                                      '1901-1916 (2020), doi:10.1021/acs.macromol.9b02428 '
                                      '(arXiv:1808.03509v3), Tables 1 and 3, Eq. 6',
                                      '348 -> 362 K (T_g + 100, T_g = 262 K) by the universal Vogel-Fulcher '
                                      'shift, C_VF = 2072 K, T_VF = T_g - 51.6 K: R. Everaers, H. A. '
                                      "Karimi-Varzaneh, F. Fleck, N. Hojdis, C. Svaneborg, 'Kremer-Grest "
                                      'models for commodity polymer melts: linking theory, experiment, and '
                                      "simulation at the Kuhn scale', Macromolecules 53, 1901-1916 (2020), "
                                      'doi:10.1021/acs.macromol.9b02428 (arXiv:1808.03509v3), Eqs. 34-35',
                                      'above 362 K Arrhenius with Ea = 40 kJ/mol, the middle of the 36-44 '
                                      'kJ/mol reported for isotactic PP melts (no atactic value given): J. A. '
                                      "Resch, 'Elastic and viscous properties of polyolefin melts with "
                                      "different molecular structures investigated in shear and elongation', "
                                      'Dr.-Ing. thesis, Universitaet Erlangen-Nuernberg (2010), p. 25 '
                                      '(review, citing Laun 1987, Gabriel 2001, Stadler et al. 2008; '
                                      'Wasserman and Graessley 1996, Eckstein et al. 1997, Auhl 2006), '
                                      'https://d-nb.info/1001254694/34'],
              'spin_factor': ['upper end of the bracket [1/4, 1] of polymer-pm report R-003 v2 (2026-09-24) '
                              'section (b): 1/4 is the singlet fraction of two uncorrelated doublets; the RMG '
                              'rates entering as k_act are thermal bulk rate constants that may already '
                              'contain spin statistics, so 1/4 on k_D could double count; 1.0 as the PS '
                              'reference arm (rmgpy/kmc/met.py); the --also-well-mixed bound brackets '
                              'transport from above']}},
 'diagnostics': {}}, cutoff=4)

simpleReactor(
    temperature=(700, 'K'),
    pressure=(1.0, 'bar'),
    initialMoleFractions={'N2': 1.0},
    terminationTime=(1.0, 's'),
)
model(toleranceKeepInEdge=0.0, toleranceMoveToCore=0.1, toleranceInterruptSimulation=0.1)
options(units='si', generateOutputHTML=False)
