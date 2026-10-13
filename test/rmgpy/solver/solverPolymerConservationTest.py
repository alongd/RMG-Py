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


def test_deck_shape_carbon_closure_with_leak_and_proxy_collision_row():
    """
    Reduced-deck carbon closure: the generation-stamped proxy-collision row
    (hexyl + PE <=> PE + PE_mod) plus a pdep leak off the gas radical, run
    through the real solver residual. dC/dt over every carbon-bearing state
    (gas species plus 2 C per mu1 unit) must vanish. On the pre-fix code the
    leak debits the gas radical and the row injects mu1-only PE_mod inventory.
    """
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
    from polymerProxyCollisionTest import CASES, _stamped
    from rmgpy.solver.polymer import PolymerPoolConfig
    from rmgpy.molecule import Molecule

    _pool, stamped, _prods = _stamped(CASES["PE"])

    def sp(smi, label):
        s = Species(molecule=[Molecule(smiles=smi)])
        s.label = label
        return s

    pe, pe_mod = sp("CCCCCC", "PE"), sp("CC[CH]CCC", "PE_mod")
    hexyl = sp("[CH2]CCCCC", "hexyl")
    mu = {k: [sp("CO", k + "_mu0"), sp("C=O", k + "_mu1"), sp("C#N", k + "_mu2")]
          for k in ("PE", "PE_mod")}
    core = [pe] + mu["PE"] + [pe_mod] + mu["PE_mod"] + [hexyl]
    mask = np.array([False] * 8 + [True], dtype=bool)
    rxn = Reaction(reactants=[hexyl, pe], products=[pe, pe_mod],
                   kinetics=Arrhenius(A=(1.0e3, "m^3/(mol*s)"), n=0.0,
                                      Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
                   reversible=False)
    rxn.is_end_group_reaction = False
    rxn.polymer_flux_archetype = stamped.polymer_flux_archetype
    rxn.polymer_refused = stamped.polymer_refused
    rxn.polymer_refused_accumulating = stamped.polymer_refused_accumulating
    pools = [PolymerPoolConfig(label=l, xs=3, explicit_dp_to_species_index={},
                               mu_indices=m, monomer_poly_index=None,
                               k_scission=0.0, k_unzip=0.0, tail_kinetics=None,
                               monomer_mw_g_mol=28.05)
             for l, m in (("PE", (1, 2, 3)), ("PE_mod", (5, 6, 7)))]
    rs = HybridPolymerSystem(
        T=1000.0, P=1.0e5, initial_mole_fractions={hexyl: 2.0}, V_poly=1.0,
        polymer_pools=pools, mass_transfer=[], gas_species_mask=mask.copy(),
        constant_gas_volume=False,
        initial_polymer_moments={"PE": (1.0, 50.0, 3000.0),
                                 "PE_mod": (1.0, 50.0, 3000.0)},
        termination=[])
    rs.initialize_model(core, [rxn], [], [], pdep_networks=[_FakeNetwork([hexyl], 0.7)])
    dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
    dC = 6.0 * dn[8] + 2.0 * dn[2] + 2.0 * dn[6]
    assert abs(dC) < 1e-12, f"carbon not conserved: dC/dt = {dC:.3e} mol C/s"
    assert rs.reaction_refused[0] == 1
