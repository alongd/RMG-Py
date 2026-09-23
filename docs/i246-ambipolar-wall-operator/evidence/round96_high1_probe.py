import os, sys, numpy as np, importlib.util
sys.path.insert(0, os.getcwd())
spec = importlib.util.spec_from_file_location("wt", os.path.join(os.getcwd(),"test/rmgpy/solver/plasmaWallTest.py"))
wt = importlib.util.module_from_spec(spec); spec.loader.exec_module(wt)
from rmgpy.solver.termination import TerminationSteadyState, TerminationTime
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.kinetics import Arrhenius
EV=wt.EV_TO_K

class Traced(TerminationSteadyState):
    def __init__(self,*a,**k):
        super().__init__(*a,**k); self.log=[]
    def update(self,yn,tn,yp,tp,fl,labels=None,external_residual=float('nan'),external_armed=False):
        g,gl=self.compute_residual(yn,tn,yp,tp,fl,labels=labels)
        out=super().update(yn,tn,yp,tp,fl,labels=labels,external_residual=external_residual,external_armed=external_armed)
        self.log.append((tn,g,gl,external_residual,external_armed,self.armed,self.streak,out)); return out

def build(slow_k, xseed, source=1e5, gamma=1.0, te=3.0, term=None):
    e=Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar=Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    arp=Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    arx=Species(label='Arx').from_smiles('[He]')
    imf={e:0.0,arp:0.0,ar:1.0-xseed,arx:xseed}
    r=PlasmaReactor((wt.TGAS,'K'),(wt.P_NOMINAL,'Pa'),imf,(te*EV,'K'),n_sims=1,
        termination=term or [], diffusion_length=(wt._diffusion_length(),'m'),
        ion_reduced_mobility=(wt.MU0_AR_IN_AR,'m^2/(V*s)'), wall_recycling=gamma,
        ionisation_source=(source,'m^-3/s'))
    core=[e,ar,arp,arx]
    slow=Reaction(reactants=[ar],products=[arx],reversible=False,
                  kinetics=Arrhenius(A=(slow_k,'s^-1'),n=0,Ea=(0,'J/mol')))
    r.initialize_model(core,[slow],[],[]); return r,core,[slow]

for xseed,k in [(1e-2,1e-12),(1e-2,1e-13),(1e-1,1e-11),(1e-1,1e-12)]:
    term=[Traced(tolerance=1e-8), TerminationTime((200.0,'s'))]
    r,core,rxns=build(k,xseed,term=term)
    wt._simulate(r,core,rxns)
    yv=np.asarray(r.y[:r.num_core_species],float)
    iarx=[i for i,s in enumerate(core) if s.label=='Arx'][0]
    L=term[0].log[-1]
    print(f"xseed={xseed:g} k={k:g}: steady={int(r.steady_state_reached)} t_f={r.t:.3e} "
          f"x_Arx={yv[iarx]/yv.sum():.6e} last(gen={L[1]:.2e} worst={L[2]} eA={int(L[4])} armed={int(L[5])} streak={L[6]} fired={int(L[7])})")
