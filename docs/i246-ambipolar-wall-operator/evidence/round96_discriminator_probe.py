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
        super().__init__(*a,**k); self.g=[]
    def update(self,yn,tn,yp,tp,fl,labels=None,external_residual=float('nan'),external_armed=False):
        gr,gl=self.compute_residual(yn,tn,yp,tp,fl,labels=labels)
        # electron-fraction slope (the HIGH-2-correct external residual)
        yn=np.asarray(yn,float); yp=np.asarray(yp,float)
        ie=None
        self.g.append((tn, gr, gl))
        return super().update(yn,tn,yp,tp,fl,labels=labels,external_residual=external_residual,external_armed=external_armed)

def build_high1(term):
    e=Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar=Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    arp=Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    arx=Species(label='Arx').from_smiles('[He]')
    r=PlasmaReactor((wt.TGAS,'K'),(wt.P_NOMINAL,'Pa'),{e:0.0,arp:0.0,ar:0.99,arx:0.01},(3.0*EV,'K'),
        n_sims=1,termination=term,diffusion_length=(wt._diffusion_length(),'m'),
        ion_reduced_mobility=(wt.MU0_AR_IN_AR,'m^2/(V*s)'),wall_recycling=1.0,ionisation_source=(1e5,'m^-3/s'))
    core=[e,ar,arp,arx]
    slow=Reaction(reactants=[ar],products=[arx],reversible=False,kinetics=Arrhenius(A=(1e-12,'s^-1'),n=0,Ea=(0,'J/mol')))
    r.initialize_model(core,[slow],[],[]); return r,core,[slow]

def traj(tag, r, core, rxns, term, tt):
    wt._simulate(r,core,rxns)
    g=term.g
    print(f"\n{tag}: {len(g)} steps, steady={r.steady_state_reached} t_f={r.t:.3e}")
    # sample generic-residual trajectory
    print("   t / gen_res / worst:")
    for row in [g[0]]+g[1:4]+g[len(g)//2:len(g)//2+2]+g[-4:]:
        print(f"     {row[0]:.3e}  {row[1]:.4e}  {row[2]}")
    # is gen_res (over live, finite) monotone increasing in the tail?
    tail=[x[1] for x in g[-8:] if np.isfinite(x[1])]
    print(f"   tail gen_res: {['%.2e'%v for v in tail]}")

t1=[Traced(tolerance=1e-8), TerminationTime((200.0,'s'))]
r,c,rx=build_high1(t1); traj("HIGH1 (Arx drift)", r,c,rx, t1[0], 200)

t2=[Traced(tolerance=1e-8), TerminationTime((200.0,'s'))]
r,c,rx=wt._build_reactor(wall=True, with_chemistry=False, x_ion=0.0, source=1e5, termination=t2)
traj("R93 (inert)", r,c,rx, t2[0], 200)

t3=[Traced(tolerance=1e-8), TerminationTime((50.0,'s'))]
r,c,rx=wt._build_reactor(wall=True, gamma=0.0, with_chemistry=False, x_ion=1e-4, source=1e22, termination=t3)
traj("HIGH2 (shrinking inventory)", r,c,rx, t3[0], 50)
