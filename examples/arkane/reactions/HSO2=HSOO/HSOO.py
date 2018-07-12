cantherm_species_dictionary = {'label': 'HSOO',
                               'timestamp': '1528897390.17',
                               'author': 'alongd',
                               'SMILES': '[O]OS',
                               'adjacencyList': 'multiplicity 2\n1 S u0 p2 c0 {2,S} {4,S}\n2 O u0 p2 c0 {1,S} {3,S}\n3 O u1 p2 c0 {2,S}\n4 H u0 p0 c0 {1,S}\n',
                               'InChI': 'InChI=1S/HO2S/c1-2-3/h3H',
                               'InChIKey': 'PWUAZNLOWMSNQX-UHFFFAOYSA-N',
                               'levelOfTheory': 'CCSD(T)-F12a/cc-pVTZ-F12//B3LYP/6-311G(2d,d,p), rotors at B3LYP/6-311G(2d,pd)',
                               'modelChemistry': 'CCSD(T)-F12/cc-pVTZ-F12',
                               'frequencyScaleFactor': 0.975,
}

CanthermSpecies(
    label = 'HSOO',
    SMILES = '[O]OS',
    adjacencyList = 'multiplicity 2\n1 S u0 p2 c0 {2,S} {4,S}\n2 O u0 p2 c0 {1,S} {3,S}\n3 O u1 p2 c0 {2,S}\n4 H u0 p0 c0 {1,S}\n',
    InChI = 'InChI=1S/HO2S/c1-2-3/h3H',
    conformer = Conformer(
        E0 = (115.143, 'kJ/mol'),
        modes = [
            IdealGasTranslation(mass=(65.0667, 'amu')),
            NonlinearRotor(
                inertia = ([91.5258, 82.869, 8.65696], 'amu*angstrom^2'),
                symmetry = 1,
            ),
            HarmonicOscillator(
                frequencies = ([266.165, 512.727, 890.356, 1175.41, 2529.35], 'cm^-1'),
            ),
            HinderedRotor(
                inertia = (4.80315, 'amu*angstrom^2'),
                symmetry = 1,
                fourier = (
                    [
                        [0.305495, -0.899654, -0.432911, -0.0818798, -0.0102349],
                        [0.000181646, 0.000219011, 0.000120536, 1.95099e-05, -1.8673e-05],
                    ],
                    'kJ/mol',
                ),
            ),
        ],
        spinMultiplicity = 2,
        opticalIsomers = 1,
        mass = ([31.9721, 1.00783, 15.9949, 15.9949], 'amu'),
        number = [16, 1, 8, 8],
        coordinates = (
            [
                [-6.95052, 3.36369, 0.00013063],
                [-7.09087, 5.91987, -0.000118519],
                [-3.50502, 3.51303, 0.000114712],
                [-2.45803, 1.33214, 0.000328714],
            ],
            'angstroms',
        ),
    ),
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [3.80904, 0.0201786, -7.63088e-05, 1.6953e-07, -1.4166e-10, 13862.2, 10.1747],
                Tmin = (10, 'K'),
                Tmax = (363.844, 'K'),
            ),
            NASAPolynomial(
                coeffs = [4.45296, 0.00805632, -5.54113e-06, 1.76698e-09, -2.12266e-13, 13848.8, 8.17796],
                Tmin = (363.844, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (115.254, 'kJ/mol'),
        Cp0 = (33.2579, 'J/(mol*K)'),
        CpInf = (78.9875, 'J/(mol*K)'),
    ),
    symmetryNumber = 1,
    collisionModel = TransportData(epsilon=(3625.11, 'J/mol'), sigma=(5.7, 'angstrom')),
    molecularWeight = (64.9697, 'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(0.5718, 'kcal/mol'), T0=(300, 'K'), n=0.85),
    chemkin_thermo_string = 'HSOO                    H   1S   1O   2     G    10.000  3000.000  363.84      1\n 4.45295619E+00 8.05631845E-03-5.54113213E-06 1.76698187E-09-2.12266359E-13    2\n 1.38487574E+04 8.17796073E+00 3.80904432E+00 2.01786111E-02-7.63088230E-05    3\n 1.69530223E-07-1.41660098E-10 1.38622319E+04 1.01747051E+01                   4\n',
)
