cantherm_species_dictionary = {'label': 'HSO2',
                               'timestamp': '1528897389.74',
                               'author': 'alongd',
                               'SMILES': 'O=[SH]=O',
                               'adjacencyList': '1 S u0 p0 c0 {2,D} {3,D} {4,S}\n2 O u0 p2 c0 {1,D}\n3 O u0 p2 c0 {1,D}\n4 H u0 p0 c0 {1,S}\n',
                               'InChI': '',
                               'InChIKey': 'LBXOVWWXGPSSJU-UHFFFAOYSA-N',
                               'levelOfTheory': 'CCSD(T)-F12a/cc-pVTZ-F12//B3LYP/6-311G(2d,d,p), rotors at B3LYP/6-311G(2d,pd)',
                               'modelChemistry': 'CCSD(T)-F12/cc-pVTZ-F12',
                               'frequencyScaleFactor': 0.975,
}

CanthermSpecies(
    label = 'HSO2',
    SMILES = 'O=[SH]=O',
    adjacencyList = '1 S u0 p0 c0 {2,D} {3,D} {4,S}\n2 O u0 p2 c0 {1,D}\n3 O u0 p2 c0 {1,D}\n4 H u0 p0 c0 {1,S}\n',
    InChI = '',
    conformer = Conformer(
        E0 = (-164.265, 'kJ/mol'),
        modes = [
            IdealGasTranslation(mass=(65.0667, 'amu')),
            NonlinearRotor(
                inertia = ([54.3242, 62.3333, 10.2918], 'amu*angstrom^2'),
                symmetry = 1,
            ),
            HarmonicOscillator(
                frequencies = ([440.169, 779.645, 954.719, 1023.18, 1218.84, 2128.76], 'cm^-1'),
            ),
        ],
        spinMultiplicity = 2,
        opticalIsomers = 1,
        mass = ([15.9949, 31.9721, 1.00783, 15.9949], 'amu'),
        number = [8, 16, 1, 8],
        coordinates = (
            [
                [-3.73877, 3.7679, 0.420763],
                [-6.33458, 3.44625, -0.480321],
                [-7.76505, 5.17143, 0.917223],
                [-7.72164, 1.05852, -0.31166],
            ],
            'angstroms',
        ),
    ),
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [4.05015, -0.00434077, 5.34925e-05, -8.86627e-08, 4.67111e-11, -19755.8, 8.27986],
                Tmin = (10, 'K'),
                Tmax = (616.212, 'K'),
            ),
            NASAPolynomial(
                coeffs = [3.04912, 0.0115323, -7.96726e-06, 2.51897e-09, -2.98402e-13, -19810.4, 11.1803],
                Tmin = (616.212, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (-164.262, 'kJ/mol'),
        Cp0 = (33.2579, 'J/(mol*K)'),
        CpInf = (83.1447, 'J/(mol*K)'),
    ),
    symmetryNumber = 1,
    collisionModel = TransportData(epsilon=(3625.11, 'J/mol'), sigma=(5.7, 'angstrom')),
    molecularWeight = (64.9697, 'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(0.5718, 'kcal/mol'), T0=(300, 'K'), n=0.85),
    chemkin_thermo_string = 'HSO2                    H   1S   1O   2     G    10.000  3000.000  616.21      1\n 3.04911514E+00 1.15322916E-02-7.96726395E-06 2.51896543E-09-2.98402330E-13    2\n-1.98104421E+04 1.11803415E+01 4.05014868E+00-4.34077417E-03 5.34924877E-05    3\n-8.86627204E-08 4.67111465E-11-1.97558179E+04 8.27985723E+00                   4\n',
)
