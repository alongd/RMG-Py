"""HIGH 2 timing probe + LOW threshold probe."""
import sys, math
import numpy as np
sys.path.insert(0, "test/rmgpy/solver")
import plasmaWallTest as T
from rmgpy.solver.termination import TerminationSteadyState, TerminationTime
import rmgpy.constants as constants

ATOL=1e-16
class Rec(TerminationSteadyState):
    def __init__(self,*a,**k):
        super().__init__(*a,**k); self.log=[]
    def update(self, y_now, t_now, y_prev, t_prev, floor, labels=None):
        out = super().update(y_now,t_now,y_prev,t_prev,floor,labels)
        # record when streak first becomes >0 and the residual/electron
        self.log.append((t_now, self.residual, self.streak, self.armed, tuple(float(v) for v in y_now)))
        return out

def run(source, tol=1e-8):
    rec = Rec(tolerance=tol)
    term=[rec, TerminationTime((200.0,'s'))]
    r,core,rxns = T._build_reactor(wall=True, with_chemistry=False, x_ion=0.0,
                                   source=source, termination=term)
    ie=r.electron_index
    T._simulate(r,core,rxns)
    V=r.compute_volume(np.asarray(r.y[:r.num_core_species],float))
    nu=r.compute_nu_wall(np.asarray(r.y[:r.num_core_species],float),V)
    n_e_ss = source/nu
    print("\n--- source=%.0e  nu=%.4f  S/nu=%.4g m^-3  n_e_final=%.6g mol (%.4g m^-3) ---" % (
        source, nu, n_e_ss, r.y[ie], r.y[ie]*constants.Na/V))
    print("reached=%s residual=%.3g armed=%s streak=%d nsteps=%d" % (
        r.steady_state_reached, r.steady_state_residual, rec.armed, rec.streak, len(rec.log)))
    # find first step where streak became >=1 (residual first < tol)
    first_flat=None
    for (t,res,streak,armed,y) in rec.log:
        if streak>=1 and first_flat is None:
            first_flat=(t,res,y[ie])
    # electron fraction of saturation when streak first started
    if first_flat:
        t,res,ne = first_flat
        # n_e at that t in density
        print("first flat step: t=%.4g  residual=%.3g  n_e=%.4g mol  saturation=%.4f" % (
            t,res,ne, ne/(r.y[ie]) if r.y[ie]>0 else float('nan')))
    # print max residual over the run and where electron slope would be
    residuals=[res for (t,res,streak,armed,y) in rec.log if np.isfinite(res)]
    print("max finite residual over run:", max(residuals) if residuals else None)
    # electron log-log slope at a few t
    for nut in (0.1,1,5,10,20,25):
        Re = nut/(math.exp(nut)-1) if nut>0 else 1.0
        print("   electron R at nu*t=%5.1f : %.3g" % (nut, Re))
    return rec

for s in (1e5, 1e12, 1e20):
    run(s)

print("\n\n====== LOW: what source magnitudes deliver a representable flux? ======")
# get_non_chemical_char_rate squares the delivered residuals. Find where the
# squared flux underflows. Delivered volumetric source S (m^-3/s) -> mol/m^3/s
# via /Na. The char rate is in mol/m^3/s. Its square underflows below sqrt(DBL_MIN).
Na=constants.Na
print("DBL_MIN normal =", sys.float_info.min, " sqrt=", math.sqrt(sys.float_info.min))
print("DBL_TRUE_MIN subnormal =", 5e-324)
for s in (1e-300, 1e-150, 1e-140, 1e-100, 1e5):
    r,core,rxns = T._build_reactor(wall=True, with_chemistry=False, x_ion=0.0, source=s)
    ncr = r.get_non_chemical_char_rate()
    # what is the raw per-species source rate in core-rate units?
    print("S=%.0e : non_chem_char_rate=%.3g  source/Na=%.3g  (source/Na)^2=%.3g" % (
        s, ncr, s/Na, (s/Na)**2))
