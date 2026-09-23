"""Round 90 red states, on the built module BEFORE any fix.

Each block prints RED (defect present) / GREEN (absent) so the same script is the
before/after witness. Run with PATH=rmg_env/bin, PYTHONPATH=worktree, MPLCONFIGDIR set.
"""
import os, sys, warnings
import numpy as np
sys.path.insert(0, os.environ['PYTHONPATH'])

# Reuse the test fixtures so the reactors are exactly the suite's.
sys.path.insert(0, os.path.join(os.environ['PYTHONPATH'], 'test', 'rmgpy', 'solver'))
import plasmaWallTest as T
from rmgpy.solver.plasma import _coerce_bool_flag, PlasmaStateError
from rmgpy.solver import plasma as plasma_mod
from rmgpy import constants

SimulatorSettings = T.SimulatorSettings
ModelSettings = T.ModelSettings


def _sim(reactor, core, rxns, t):
    """Drive simulate() to time t via the suite's production entry."""
    reactor.termination = [T.TerminationTime((t, 's'))]
    return T._simulate(reactor, core, rxns)


print("=== HIGH 1: external source cannot ignite neutral gas (zero e-, source declared) ===")
try:
    r, core, rxns = T._build_reactor(wall=True, with_chemistry=False, x_ion=0.0, source=1.0e5)
    # If it built, try to actually integrate from zero electrons.
    try:
        _sim(r, core, rxns, 1.0e-3)
        print("  GREEN: zero-electron deck WITH a source integrates via simulate()")
    except Exception as e:
        print("  RED(sim): built but simulate() failed:", type(e).__name__, str(e)[:90])
except Exception as e:
    print("  RED(init): zero-electron deck WITH a source is refused at init:",
          type(e).__name__, str(e)[:110])


print("\n=== HIGH 2: neutral floor is extensive (history-dependent), gating an intensive law ===")
# The reactor is isobaric, so the neutral number density is pinned at ~P/kT no matter how
# many moles remain: a wall that pumps mass out shrinks the neutral MOLES toward the old
# 1e-6-mol floor while the intensive DENSITY the transport law uses stays put. Evaluate a
# state with small neutral moles but normal density: an intensive criterion judges it on
# the density (accept, nu unclamped); the extensive moles floor judges it on its history.
r, core, rxns = T._build_reactor(wall=True, with_chemistry=False, x_ion=1.0e-9)  # ~1 mol init
ie, iar, iarp = T._indices(r)
y = np.zeros(r.num_core_species, float)
y[iar] = 1.0e-7            # neutral MOLES far below the old 1e-6-mol floor
y[iarp] = 1.0e-13
y[ie] = 1.0e-13           # alpha ~ 1e-6, far under the ceiling
V = r.compute_volume(y)
n_neutral = y[iar] * constants.Na / V
nu = r.compute_nu_wall(y, V)
lam = r.diffusion_length.value_si
mu_i = r.ion_reduced_mobility.value_si * r.mobility_reference_density / n_neutral
d_a = mu_i * (constants.R / constants.Na) * r.Te.value_si / constants.e
nu_closed = d_a / (lam * lam)               # UNCLAMPED closed form at this density
clamped = abs(nu / nu_closed - 1.0) > 1.0e-9
try:
    r.check_wall_support(y); verdict = 'accept'
except PlasmaStateError:
    verdict = 'refuse'
print(f"  neutral moles={y[iar]:.1e} mol, normal density n_neutral={n_neutral:.4e} m^-3")
print(f"  nu_wall={nu:.6e}  unclamped closed-form={nu_closed:.6e}  clamped={clamped}")
print(f"  check_wall_support verdict={verdict}")
if clamped or verdict == 'refuse':
    print("  RED: an intensive state judged by extensive moles — clamped and/or refused "
          "on its mole count")
else:
    print("  GREEN: judged on density — nu_wall unclamped and the state accepted")


print("\n=== MEDIUM 1a: non-finite nu_wall marked wall_flux 'available' ===")
try:
    r, core, rxns = T._build_reactor(wall=True, with_chemistry=False, mu0=1.0e308)
    r._latch_wall_diagnostics(r.y0, r.compute_volume(r.y0), 0.0)
    wf = np.asarray(r.wall_flux)
    avail = r.wall_energy_availability.get('wall_flux')
    finite = np.isfinite(wf).all()
    print(f"  wall_flux finite={finite}, availability={avail!r}")
    if (not finite) and avail == 'available':
        print("  RED: a non-finite wall_flux is reported 'available'")
    else:
        print("  GREEN: non-finite wall_flux is not 'available'")
except Exception as e:
    print("  note:", type(e).__name__, str(e)[:90])


print("\n=== MEDIUM 1b: sub-underflow diffusion length raises a raw ZeroDivisionError ===")
try:
    r, core, rxns = T._build_reactor(wall=True, with_chemistry=False, lam=1.0e-200)
    r.compute_nu_wall(r.y0, r.compute_volume(r.y0))
    print("  GREEN: no raw crash from a tiny diffusion length")
except PlasmaStateError as e:
    print("  GREEN(refused by name):", str(e)[:90])
except ZeroDivisionError as e:
    print("  RED: raw ZeroDivisionError:", str(e)[:90])
except Exception as e:
    print("  other:", type(e).__name__, str(e)[:90])


print("\n=== MEDIUM 2: _coerce_bool_flag still ends in bool(value) ===")
bad = [("int 2", 2), ("float 0.5", 0.5), ("float NaN", float('nan')),
       ("object()", object()), ("empty list", []), ("empty dict", {})]
red = False
for label, v in bad:
    try:
        got = _coerce_bool_flag(v, 'quasineutral_electron', 'id')
        print(f"  {label:12s} -> {got!r}  (coerced, not refused)")
        red = True
    except PlasmaStateError:
        print(f"  {label:12s} -> refused")
print("  RED: uninterpretable values are coerced, not refused" if red
      else "  GREEN: uninterpretable values are refused")


print("\n=== LOW: base.pyx char_rate=0 makes 0/0 NaN rate ratios (wall-only run) ===")
r, core, rxns = T._build_reactor(wall=True, with_chemistry=False, x_ion=1.0e-6)
try:
    with np.errstate(invalid='raise', divide='raise'):
        _sim(r, core, rxns, 1.0e-4)
    print("  GREEN: wall-only simulate() produced no invalid/divide float error")
except FloatingPointError as e:
    print("  RED: wall-only simulate() hit a float error in the rate-ratio divide:", str(e)[:80])
except Exception as e:
    print("  other:", type(e).__name__, str(e)[:100])
