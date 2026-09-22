"""Round 93 premise probes against the CURRENT built module (ebe7e76c1)."""
import sys, math
import numpy as np
sys.path.insert(0, "test/rmgpy/solver")
import plasmaWallTest as T
from rmgpy.solver.plasma import PlasmaReactor, PLASMA_LOSCHMIDT, PLASMA_NET_CHARGE_RTOL
from rmgpy.solver.termination import TerminationSteadyState, TerminationTime
from rmgpy.exceptions import PlasmaStateError
import rmgpy.constants as constants

ATOL = 1e-16  # SimulatorSettings default core-species atol (base.pyx:232)

def banner(s): print("\n" + "="*70 + "\n" + s + "\n" + "="*70)

# -------------------------------------------------------------- HIGH 1
banner("HIGH 1: quasineutral ignition trajectory (net/magnitude vs inventory)")

# Subclass that records net/magnitude and applies a candidate stand-down:
# skip the net-charge guard when magnitude < ATOL (sub-resolution charged
# inventory), else enforce the relative bound. All OTHER checks skipped so the
# run proceeds and we can watch the trajectory. This is a PROTOTYPE of the fix.
class Traced(PlasmaReactor):
    trace = []
    def check_wall_support(self, y):
        if not self.has_wall or self.neutral_heavy_mask is None:
            return
        n_e = y[self.electron_index]
        n_ion = 0.0
        for j in range(self.num_core_species):
            if j != self.electron_index and self.species_charges[j] > 0:
                n_ion += y[j]
        net = n_ion - n_e
        mag = n_ion + n_e
        rel = abs(net)/mag if mag > 0 else 0.0
        Traced.trace.append((self.t, net, mag, rel))
        if mag >= ATOL and abs(net) > PLASMA_NET_CHARGE_RTOL * mag:
            raise PlasmaStateError("stand-down prototype refused: rel=%r at mag=%r" % (rel, mag))

def build_traced(**kw):
    e, ar, arp = T._argon_species()
    x = kw.get('x_ion', 0.0)
    imf = {e: x, arp: x, ar: 1.0 - 2*x}
    r = Traced((T.TGAS,'K'),(T.P_NOMINAL,'Pa'),imf,(T.TE_NOMINAL_EV*T.EV_TO_K,'K'),
               n_sims=1, termination=kw.get('termination',[]),
               quasineutral_electron=kw.get('quasineutral',True),
               diffusion_length=(T._diffusion_length(),'m'),
               ion_reduced_mobility=(T.MU0_AR_IN_AR,'m^2/(V*s)'),
               wall_recycling=1.0,
               ionisation_source=(kw['source'],'m^-3/s'))
    core=[e,ar,arp]; rxns=[]
    r.initialize_model(core,rxns,[],[])
    return r, core, rxns

# First: confirm the STOCK failure without the stand-down.
try:
    r,core,rxns = T._build_reactor(wall=True, with_chemistry=False, x_ion=0.0,
                                   source=1e5, quasineutral=True)
    r.termination=[TerminationTime((1e-4,'s'))]
    T._simulate(r,core,rxns)
    print("STOCK quasineutral run COMPLETED (unexpected) t=", r.t)
except PlasmaStateError as ex:
    print("STOCK quasineutral FAILS as reviewer reports:")
    print("  ", str(ex)[:160])

# Now the prototype stand-down: does it proceed, and does rel drop < 1e-6 once mag>ATOL?
for src in (1e5, 1e12, 1e20):
    Traced.trace = []
    try:
        r,core,rxns = build_traced(source=src, quasineutral=True)
        r.termination=[TerminationTime((1e-4,'s'))]
        T._simulate(r,core,rxns)
        tr = np.array([(t,mag,rel) for (t,net,mag,rel) in Traced.trace if mag>0])
        below = tr[tr[:,1] < ATOL]
        above = tr[tr[:,1] >= ATOL]
        print("source=%.0e: COMPLETED t=%.4g  n_e=%.4g mol" % (src, r.t, r.y[r.electron_index]))
        if len(below): print("   below-atol steps: max rel imbalance = %.3g (n=%d)" % (below[:,2].max(), len(below)))
        if len(above): print("   above-atol steps: max rel imbalance = %.3g (n=%d)  <-- must be < 1e-6" % (above[:,2].max(), len(above)))
    except PlasmaStateError as ex:
        print("source=%.0e: PROTOTYPE STILL FAILED: %s" % (src, str(ex)[:120]))

# -------------------------------------------------------------- HIGH 2
banner("HIGH 2: physical-branch steady state not recognised")
term = [TerminationSteadyState(tolerance=1e-8), TerminationTime((200.0,'s'))]
r,core,rxns = T._build_reactor(wall=True, with_chemistry=False, x_ion=0.0,
                               source=1e5, termination=term)
ie = r.electron_index
T._simulate(r,core,rxns)
V = r.compute_volume(np.asarray(r.y[:r.num_core_species],float))
nu = r.compute_nu_wall(np.asarray(r.y[:r.num_core_species],float), V)
n_e_density = r.y[ie]*constants.Na/V
print("steady_state_reached  =", r.steady_state_reached)
print("steady_state_residual =", r.steady_state_residual)
print("t_final               =", r.t)
print("nu_wall               =", nu)
print("electron moles        = %.6g   (atol=%.0e -> %s)" % (
      r.y[ie], ATOL, "ABOVE floor" if r.y[ie]>ATOL else "BELOW floor -> invisible to residual"))
print("n_e density           = %.6g m^-3" % n_e_density)
print("S/nu_wall             = %.6g m^-3" % (1e5/nu))
print("ratio n_e/(S/nu)      = %.6g" % (n_e_density/(1e5/nu)))
# is the armed latch the blocker?
for tm in term:
    if isinstance(tm, TerminationSteadyState):
        print("criterion armed       =", tm.armed, " streak=", tm.streak)

# -------------------------------------------------------------- MEDIUM 1
banner("MEDIUM 1: finite inputs -> infinite nu_wall admitted")
for label, kw in [("mu0=1e308", dict(mu0=1e308)), ("Lambda=1e-160", dict(lam=1e-160))]:
    try:
        r,_,_ = T._build_reactor(wall=True, with_chemistry=False, **kw)
        y = T._state_at(r, 1e-6); V=r.compute_volume(y)
        nu = r.compute_nu_wall(y,V)
        print("%s: ADMITTED, nu_wall=%r finite=%s" % (label, nu, np.isfinite(nu)))
        lam=r.diffusion_length.value_si
        print("    lam=%r lam^2=%r (normal-min=%.3g)" % (lam, lam*lam, sys.float_info.min))
    except (PlasmaStateError, Exception) as ex:
        print("%s: refused at construction: %s" % (label, str(ex)[:100]))

# -------------------------------------------------------------- MEDIUM 3
banner("MEDIUM 3: floor moves under (Nref*c, mu0/c) reparameterisation")
base_mu0 = T.MU0_AR_IN_AR
for c in (1.0, 1e3, 1e-3):
    e, ar, arp = T._argon_species()
    imf = {e:1e-6, arp:1e-6, ar:1-2e-6}
    r = PlasmaReactor((T.TGAS,'K'),(T.P_NOMINAL,'Pa'),imf,(T.TE_NOMINAL_EV*T.EV_TO_K,'K'),
                      n_sims=1, termination=[],
                      diffusion_length=(T._diffusion_length(),'m'),
                      ion_reduced_mobility=(base_mu0/c,'m^2/(V*s)'),
                      mobility_reference_density=(PLASMA_LOSCHMIDT*c,'m^-3'),
                      wall_recycling=1.0)
    r.initialize_model([e,ar,arp],[],[],[])
    y = T._state_at(r,1e-6); V=r.compute_volume(y)
    nu = r.compute_nu_wall(y,V)
    print("c=%.0e: floor=%.6g  nu_wall=%.6g  (mu0*Nref=%.6g)" % (
          c, r.wall_neutral_density_floor, nu, (base_mu0/c)*(PLASMA_LOSCHMIDT*c)))

# -------------------------------------------------------------- LOW
banner("LOW: subnormal source lifts the zero-electron guard, delivers nothing")
try:
    r,core,rxns = T._build_reactor(wall=True, with_chemistry=False, x_ion=0.0, source=1e-300)
    print("source=1e-300 ADMITTED zero-electron deck. source.value_si=%r" % r.ionisation_source.value_si)
    print("   normal-min double = %.3g ; is subnormal? %s" % (
          sys.float_info.min, r.ionisation_source.value_si < sys.float_info.min))
    print("   get_non_chemical_char_rate at y0 =", r.get_non_chemical_char_rate())
except PlasmaStateError as ex:
    print("source=1e-300 refused:", str(ex)[:120])

print("\nDONE")
