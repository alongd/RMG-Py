name = 'MissingCpInfFixture'
shortDesc = 'Synthetic fixed-level fixture deliberately missing CpInf'

entry(index=0, label='N2v0', molecule="""
vibrationallevel 0
1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
""",
    thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'),
                     Cpdata=([29.1]*7,'J/(mol*K)'), H298=(0,'kJ/mol'),
                     S298=(190,'J/(mol*K)'), Cp0=(29.1,'J/(mol*K)')))


entry(index=1, label='N2v1', molecule="""
vibrationallevel 1
1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
""",
    thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'),
                     Cpdata=([29.1]*7,'J/(mol*K)'), H298=(0,'kJ/mol'),
                     S298=(190,'J/(mol*K)'), Cp0=(29.1,'J/(mol*K)')))
