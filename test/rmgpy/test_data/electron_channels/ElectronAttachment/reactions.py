name = "ElectronAttachment"
shortDesc = "Dissociative attachment fixture"
entry(
    index = 1,
    label = "O2 + e- => O- + O",
    reversible = False,
    kinetics = TwoTemperaturePlasma(
        A = (6.44e14, 'cm^3/(mol*s)'), n = -1.391,
        Ea_g = (6.26, 'eV/molecule'), Ea_e = (6.26, 'eV/molecule'),
        T0 = (11604.51812, 'K'), Tmin = (5802, 'K'), Tmax = (58023, 'K'),
    ),
)
