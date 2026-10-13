"""
Conservation locks for the polymer (method-of-moments) reactor:

* pdep network leak flux is a RECORDED diagnostic for model enlargement only
  (as in ``simple.pyx``); it must never be subtracted from the core balance.
* a gas radical abstracting an H from a pool chain, whose saturated partner is
  isomorphic to the pool proxy, must never fall into the legacy mu1-only flux.
"""
import numpy as np

from rmgpy.kinetics import Arrhenius
from rmgpy.reaction import Reaction
from rmgpy.solver.polymer import HybridPolymerSystem
from rmgpy.species import Species


class _FakeNetwork:
    """Minimal stand-in for a PDepNetwork: a source list and a leak coefficient."""

    def __init__(self, source, k):
        self.source = source
        self._k = k

    def get_leak_coefficient(self, T, P):
        return self._k


def test_network_leak_is_recorded_but_not_subtracted_from_core():
    A = Species().from_smiles("C")
    A.label = "A"
    B = Species().from_smiles("[CH3]")
    B.label = "B"
    k_leak = 0.7
    net = _FakeNetwork([A], k_leak)
    rs = HybridPolymerSystem(
        T=1000.0, P=1.0e5,
        initial_mole_fractions={A: 1.5, B: 0.0},
        V_poly=1.0, polymer_pools=[], mass_transfer=[],
        gas_species_mask=np.array([True, True], dtype=bool),
        constant_gas_volume=False, termination=[],
    )
    rs.initialize_model([A, B], [], [], [], pdep_networks=[net])
    dn_dt = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
    # no core reactions: the core balance must be untouched by the leak
    assert dn_dt[0] == 0.0, "network leak was subtracted from a core species"
    assert dn_dt[1] == 0.0
    # ... while the leak rate is still recorded for model enlargement
    conc = rs.y[0] / rs.V_gas
    assert np.isclose(rs.network_leak_rates[0], k_leak * conc)
    assert rs.network_leak_rates[0] > 0.0
