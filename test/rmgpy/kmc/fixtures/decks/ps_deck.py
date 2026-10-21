# Linear atactic polystyrene, H caps; the meltStateTest PS declaration without styrene.
# Structure and a declared MET block; density_kg_m3 is the volume basis only.
# MET transport: literature values, each field cited in transport_model.sources.
database(
    thermoLibraries=['primaryThermoLibrary'],
    reactionLibraries=[],
    seedMechanisms=[],
    kineticsDepositories='default',
    kineticsFamilies='default',
    kineticsEstimator='rate rules',
)
species(label='N2', reactive=False, structure=SMILES('N#N'))
polymer('polystyrene (CH2-CH(C6H5))', linear_homopolymer_spec={'density_basis': 'dry_polymer',
 'density_kg_m3': 1050.0,
 'diagnostics': {},
 'end_groups': {'head': {'atoms': [{'electronic_state': {'radical_electrons': 0},
                                    'element': 'H',
                                    'hydrogens': 0,
                                    'key': 'h'}],
                         'bonds': [],
                         'electronic_state': {'multiplicity': 1},
                         'ports': {'cap': {'atom': 'h', 'bond_order': 1}}},
                'tail': {'atoms': [{'electronic_state': {'radical_electrons': 0},
                                    'element': 'H',
                                    'hydrogens': 0,
                                    'key': 'h'}],
                         'bonds': [],
                         'electronic_state': {'multiplicity': 1},
                         'ports': {'cap': {'atom': 'h', 'bond_order': 1}}}},
 'family_manifest': {},
 'initial_ensemble': {'dp_counts': {30: 100}},
 'initial_mass_kg': 1.0,
 'initial_mobile_counts': {},
 'label': 'polystyrene (CH2-CH(C6H5))',
 'library_manifest': {},
 'locality_policy': {'bounds': {'attachment_classes': [],
                                'bond_orders': [1, 1.5, 2, 3],
                                'elements': ['C', 'H'],
                                'formal_charges': [0],
                                'isotopes': [-1],
                                'lone_pairs': [0],
                                'maximum_explicit_hydrogens_per_heavy_atom': 4,
                                'maximum_family_arity': 2,
                                'maximum_free_small_species_heavy_atoms': 8,
                                'maximum_heavy_atom_degree': 4,
                                'maximum_root_tuple_size': 4,
                                'maximum_total_hydrogens_on_free_small_species': 12,
                                'radical_electrons': [0, 1],
                                'spin_multiplicities': [1, 2, 3],
                                'topology_metadata_schema': {}},
                     'radius': 2},
 'repeat_unit': {'atoms': [{'electronic_state': {'radical_electrons': 0},
                            'element': 'C',
                            'hydrogens': 2,
                            'key': 'c1'},
                           {'electronic_state': {'radical_electrons': 0},
                            'element': 'C',
                            'hydrogens': 1,
                            'key': 'c2'},
                           {'electronic_state': {'radical_electrons': 0},
                            'element': 'C',
                            'hydrogens': 0,
                            'key': 'p1'},
                           {'electronic_state': {'radical_electrons': 0},
                            'element': 'C',
                            'hydrogens': 1,
                            'key': 'p2'},
                           {'electronic_state': {'radical_electrons': 0},
                            'element': 'C',
                            'hydrogens': 1,
                            'key': 'p3'},
                           {'electronic_state': {'radical_electrons': 0},
                            'element': 'C',
                            'hydrogens': 1,
                            'key': 'p4'},
                           {'electronic_state': {'radical_electrons': 0},
                            'element': 'C',
                            'hydrogens': 1,
                            'key': 'p5'},
                           {'electronic_state': {'radical_electrons': 0},
                            'element': 'C',
                            'hydrogens': 1,
                            'key': 'p6'}],
                 'bonds': [{'atoms': ['c1', 'c2'], 'order': 1},
                           {'atoms': ['c2', 'p1'], 'order': 1},
                           {'atoms': ['p1', 'p2'], 'order': 1.5},
                           {'atoms': ['p2', 'p3'], 'order': 1.5},
                           {'atoms': ['p3', 'p4'], 'order': 1.5},
                           {'atoms': ['p4', 'p5'], 'order': 1.5},
                           {'atoms': ['p5', 'p6'], 'order': 1.5},
                           {'atoms': ['p1', 'p6'], 'order': 1.5}],
                 'electronic_state': {'multiplicity': 1},
                 'ports': {'head': {'atom': 'c1', 'bond_order': 1}, 'tail': {'atom': 'c2', 'bond_order': 1}}},
 'schema_version': 1,
 'small_species': {},
 'transport_model': {'met': {'chain_scaling': {'entanglement_units': 142.1, 'law': 'crossover'},
                             'contact_radius_m': 6.766081442101e-10,
                             'mean_square_end_to_end_per_molar_mass_m2_mol_per_g': 4.37e-21,
                             'schema': 'kmc-met-transport/1',
                             'segment_diffusivity': {'C1': 13.39,
                                                     'C2_K': 45.28,
                                                     'D_ref_m2_s': 5.820287308966568e-11,
                                                     'T_ref_K': 473.0,
                                                     'T_wlf_K': 378.15,
                                                     'law': 'wlf'},
                             'spin_factor': 1.0},
                     'sources': {'chain_scaling': ['N_e = M_e / M0 = 14800 / 104.15, M_e from G_N0 = 0.195 '
                                                   'MPa (4/5 convention): C. Liu, J. He, E. van Ruymbeke, R. '
                                                   "Keunings, C. Bailly, 'Evaluation of different methods "
                                                   'for the determination of the plateau modulus and the '
                                                   "entanglement molecular weight', Polymer 47, 4461-4479 "
                                                   '(2006), p. 4468'],
                                 'contact_radius_m': ['sigma = 2 r, r the radius of a sphere of the McGowan '
                                                      'volume of the radical site C8H9 (8 C, 9 H, 17 bonds): '
                                                      "M. H. Abraham, J. C. McGowan, 'The use of "
                                                      'characteristic volumes to measure cavity terms in '
                                                      "reversed phase liquid chromatography', "
                                                      'Chromatographia 23, 243-246 (1987), Table 2 (C 16.35, '
                                                      'H 8.71, minus 6.56 per bond, cm3/mol), as implemented '
                                                      'in RMG-Py rmgpy/data/solvation.py set_mcgowan_volume'],
                                 'mean_square_end_to_end_per_molar_mass_m2_mol_per_g': ['<R^2>/M = 0.437 A^2 '
                                                                                        'mol/g at 413 K '
                                                                                        '(a-PS, derived by '
                                                                                        'the authors): R. '
                                                                                        'Everaers, H. A. '
                                                                                        'Karimi-Varzaneh, F. '
                                                                                        'Fleck, N. Hojdis, '
                                                                                        'C. Svaneborg, '
                                                                                        "'Kremer-Grest "
                                                                                        'models for '
                                                                                        'commodity polymer '
                                                                                        'melts: linking '
                                                                                        'theory, experiment, '
                                                                                        'and simulation at '
                                                                                        "the Kuhn scale', "
                                                                                        'Macromolecules 53, '
                                                                                        '1901-1916 (2020), '
                                                                                        'doi:10.1021/acs.macromol.9b02428 '
                                                                                        '(arXiv:1808.03509v3), '
                                                                                        'Table 1'],
                                 'segment_diffusivity': ['D_ref = kB T / zeta0 at 473 K, log10 zeta0 = -6.95 '
                                                         '(dyn s/cm) per styrene unit at T_g + 100: J. D. '
                                                         'Ferry, Viscoelastic Properties of Polymers, 3rd '
                                                         'ed., Wiley, New York (1980), Table 12-III, pp. '
                                                         '330-331',
                                                         'WLF C1 = 13.39, C2 = 45.28 K referenced to T_g = '
                                                         '378.15 K (tan delta shift factors, 150-290 C): P. '
                                                         "Lomellini, 'Williams-Landel-Ferry versus Arrhenius "
                                                         'behaviour: polystyrene melt viscoelasticity '
                                                         "revised', Polymer 33, 4983-4989 (1992), Table 2 "
                                                         '(tan delta row), p. 4986'],
                                 'spin_factor': ['upper end of the bracket [1/4, 1] of polymer-pm report '
                                                 'R-003 v2 (2026-09-24) section (b): 1/4 is the singlet '
                                                 'fraction of two uncorrelated doublets; the RMG rates '
                                                 'entering as k_act are thermal bulk rate constants that may '
                                                 'already contain spin statistics, so 1/4 on k_D could '
                                                 'double count; 1.0 as the PS reference arm '
                                                 '(rmgpy/kmc/met.py); the --also-well-mixed bound brackets '
                                                 'transport from above']}}}, cutoff=4)
simpleReactor(
    temperature=(700, 'K'),
    pressure=(1.0, 'bar'),
    initialMoleFractions={'N2': 1.0},
    terminationTime=(1.0, 's'),
)
model(toleranceKeepInEdge=0.0, toleranceMoveToCore=0.1, toleranceInterruptSimulation=0.1)
options(units='si', generateOutputHTML=False)
