"""HIGH 2 red-first: a wall-only deck (no gas-phase chemistry, gamma=0) driven
through simulate() must integrate PAST t=0 -- the wall is depleting it. On the
current .so char_rate=0 reads as 'inert' and the run terminates at t=0 claiming
the composition cannot change."""
import os, sys
import numpy as np

sys.path.insert(0, os.environ['PYTHONPATH'])
import rmgpy.constants as constants
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
from rmgpy.solver.termination import TerminationSteadyState, TerminationTime

# import the test fixtures for a faithful argon wall reactor
sys.path.insert(0, os.path.join(os.environ['PYTHONPATH'], 'test', 'rmgpy', 'solver'))
import plasmaWallTest as W

term = [TerminationSteadyState(tolerance=1e-8), TerminationTime((1.0, 's'))]
r, core, rxns = W._build_reactor(wall=True, gamma=0.0, with_chemistry=False,
                                 x_ion=1.0e-4, termination=term)

res = r.simulate(core, rxns, [], [], [], [],
                 model_settings=ModelSettings(tol_keep_in_edge=0, tol_move_to_core=1e5,
                                              tol_interrupt_simulation=1e8),
                 simulator_settings=SimulatorSettings())
terminated, _res, invalid, _ss, _sr, t_final, _conv = res

ie, i_ar, i_arp = W._indices(r)
print("terminated      =", terminated, " final t =", t_final)
print("core_species_rates =", np.array(r.core_species_rates, float))
print("wall_flux       =", np.array(r.wall_flux, float), "mol/s")
print("|wall_flux|     =", float(np.sqrt(np.sum(np.array(r.wall_flux, float)**2))))
print("nu_wall_latched =", r.nu_wall_latched, "s^-1")
print("y0              =", np.array(r.y0, float))
print("y  final        =", np.array(r.y[:r.num_core_species], float))
moved = not np.array_equal(np.array(r.y[:r.num_core_species], float), np.array(r.y0, float))
print()
print("RED expectation: the wall is removing ions/electrons, so the run must")
print("integrate past t=0 and the composition must move.")
print("  final t > 0 ?      ", t_final > 0.0)
print("  composition moved ?", moved)
assert t_final > 0.0, "REPRODUCED: wall-only deck terminated at t=0 (declared inert)"
assert moved, "REPRODUCED: composition bit-identical to start (wall made invisible)"
print("\nGREEN: wall-driven deck integrated past t=0.")
