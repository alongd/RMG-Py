name = 'ExcitedStateFixture'
shortDesc = 'Synthetic test fixture, not physical data'

entry(index=0, label='N2', molecule="""
1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
""",
    thermo=NASA(polynomials=[NASAPolynomial(coeffs=[4.5, 0, 0, 0, 0, -1000, 3], Tmin=(200,'K'), Tmax=(6000,'K'))],
    Tmin=(200,'K'), Tmax=(6000,'K'), Cp0=(37.415081781,'J/(mol*K)'), CpInf=(37.415081781,'J/(mol*K)')))

entry(index=1, label='N2v0', molecule="""
vibrationallevel 0
1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
""",
    thermo=NASA(polynomials=[NASAPolynomial(coeffs=[3.5, 0, 0, 0, 0, -1000, 3], Tmin=(200,'K'), Tmax=(6000,'K'))],
    Tmin=(200,'K'), Tmax=(6000,'K'), Cp0=(29.100619163,'J/(mol*K)'), CpInf=(29.100619163,'J/(mol*K)')))

entry(index=2, label='N2v1', molecule="""
vibrationallevel 1
1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
""",
    thermo=NASA(polynomials=[NASAPolynomial(coeffs=[3.5, 0, 0, 0, 0, 2400, 3], Tmin=(200,'K'), Tmax=(6000,'K'))],
    Tmin=(200,'K'), Tmax=(6000,'K'), Cp0=(29.100619163,'J/(mol*K)'), CpInf=(29.100619163,'J/(mol*K)')))

entry(index=3, label='N2A', molecule="""
electronicstate A3Su+
1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
""",
    thermo=NASA(polynomials=[NASAPolynomial(coeffs=[3.5, 0, 0, 0, 0, 70000, 3], Tmin=(200,'K'), Tmax=(6000,'K'))],
    Tmin=(200,'K'), Tmax=(6000,'K'), Cp0=(29.100619163,'J/(mol*K)'), CpInf=(29.100619163,'J/(mol*K)')))
