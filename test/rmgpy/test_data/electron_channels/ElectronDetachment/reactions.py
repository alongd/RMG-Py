name = "ElectronDetachment"
shortDesc = "Associative detachment fixture"
entry(
    index = 1,
    label = "O + O- => O2 + e-",
    reversible = False,
    kinetics = Arrhenius(
        A = (1.3850923748e14, 'cm^3/(mol*s)'),
        n = 0, Ea = (0, 'J/mol'), T0 = (300, 'K'),
    ),
)
