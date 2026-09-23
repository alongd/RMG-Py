import numpy as np, sys
sys.path.insert(0,'.')
sys.argv=['x']
from rmgpy.solver.plasma import PlasmaReactor, PlasmaStateError
from rmgpy.species import Species
EV_TO_K=11604.518
def argon():
    e=Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar=Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    arp=Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    return e,ar,arp
e,ar,arp=argon()
imf={e:0.0, arp:0.0, ar:1.0}
import rmgpy.constants as c
def mk(P, source):
    r=PlasmaReactor((300.0,'K'),(P,'Pa'),imf,(2.0*EV_TO_K,'K'),n_sims=1,termination=[],
        diffusion_length=(0.03,'m'), ion_reduced_mobility=(1e-3,'m^2/(V*s)'),
        wall_recycling=0.0, ionisation_source=(source,'m^-3/s'))
    r.initialize_model([e,ar,arp],[],[],[])
    return r
# very high P -> tiny V -> source*V/Na underflows though source/Na is fine
for P,source in ((1.0e37, 1.0e-250), (1.0e40, 1.0e-245)):
    try:
        r=mk(P, source); print("P=%.0e source=%.0e -> CONSTRUCTED V=%.3e source*V/Na=%.3e"%(P,source,r.V, source*r.V/c.Na))
    except PlasmaStateError as ex:
        print("P=%.0e source=%.0e -> REFUSED: %s"%(P,source,str(ex)[:90]))

print("=== MED1 red candidates ===")
for P,source in ((1.0e40, 1.0e-260), (1.0e40, 1.0e-252)):
    import rmgpy.constants as c
    try:
        r=mk(P, source); print("P=%.0e source=%.0e -> CONSTRUCTED (source/Na=%.2e, source*V/Na=%.2e)"%(P,source,source/c.Na, source*r.V/c.Na))
    except PlasmaStateError as ex:
        print("P=%.0e source=%.0e -> REFUSED: %s"%(P,source,str(ex)[:70]))
