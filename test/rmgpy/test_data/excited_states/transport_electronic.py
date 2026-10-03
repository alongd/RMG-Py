name = 'transport_electronic'
shortDesc = 'Synthetic test fixture, not physical data'
entry(index=0, label='N2A', molecule="""
electronicstate A3Su+
1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
""", transport=TransportData(shapeIndex=1, epsilon=(250,'K'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'De'), polarizability=(1.7,'angstroms^3'), rotrelaxcollnum=4))
