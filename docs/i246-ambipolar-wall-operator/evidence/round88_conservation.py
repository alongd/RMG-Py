"""Round 88 verifier item 4: particle (heavy-atom) and charge conservation across the
wall operator, with the metastable present, after all four HIGH fixes. The wall's common
loss frequency removes every charged species at one rate, so the net charge flux must be
exactly zero on a neutral state; heavy Ar nuclei are conserved because gamma of every
neutralised ion returns as Ar and (1-gamma) leaves as accounted wall inventory."""
import os, sys
import numpy as np
sys.path.insert(0, os.environ['PYTHONPATH'])
sys.path.insert(0, os.path.join(os.environ['PYTHONPATH'], 'test', 'rmgpy', 'solver'))
import rmgpy.constants as constants
import plasmaWallTest as W

for gamma in (1.0, 0.5, 0.0):
    reactor, core = W._metastable_reactor(gamma=gamma, neutralization={'Ar+': 'Ar'})
    z = np.array(reactor.species_charges, float)
    ie = reactor.electron_index
    y = np.array(reactor.y0, float)
    delta, _ = reactor.residual(0.0, y.copy(), np.zeros_like(y))
    wall = np.array(reactor.wall_loss_rates, float)

    # net charge flux to the wall: exactly zero on the neutral seed
    net_charge_flux = float(np.sum(z * wall))
    # heavy Ar atoms: every core species except the electron carries exactly one Ar.
    # d(total Ar)/dt from the wall = sum over non-electron species of wall flux; the
    # pumped fraction (1-gamma) is the only sink, so it equals -(1-gamma)*ion loss.
    heavy_flux = float(sum(wall[j] for j in range(len(core)) if j != ie))
    ion_loss = -float(sum(wall[j] for j in range(len(core))
                          if j != ie and z[j] > 0))
    expected_heavy = -(1.0 - gamma) * ion_loss

    print(f"gamma={gamma}")
    print(f"  net charge flux to wall = {net_charge_flux:.3e}  (must be ~0)")
    print(f"  heavy-atom flux         = {heavy_flux:.6e}")
    print(f"  -(1-gamma)*ion_loss     = {expected_heavy:.6e}")
    assert abs(net_charge_flux) < 1e-25, net_charge_flux
    assert np.isclose(heavy_flux, expected_heavy, rtol=1e-12, atol=1e-30)
    print(f"  wall_neutralization_energy_flux = {reactor.wall_neutralization_energy_flux!r} W "
          f"[{reactor.wall_energy_availability['wall_neutralization_energy_flux']}]")
print("\nOK: charge conserved (zero net wall current) and heavy atoms conserved at "
      "every gamma, with the metastable present.")
