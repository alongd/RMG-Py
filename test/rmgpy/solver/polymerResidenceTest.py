"""Declared residence of non-pool species (shared-chemistry slice S4).

``polymer_phase(..., residence={label: 'melt' | 'gas', '*': ...})`` says
where each non-pool species lives.  A melt-resident species is CONDENSED in
the gas mask (``PolymerPhase.get_gas_mask``, the one stage-1 classifier of
the core mask and the prospective seed): n / V_poly, with the gas thermo and
Kc at C0 = P0/RT every species already takes, i.e. the kinetic Monte Carlo's
single convention.  ``'*'`` covers every plain reactive species not named
(generated species included).  A deck that declares residence declares that
convention, so the reference-state tripwire does not refuse a reversible row
whose participants are all condensed (it is still in the census, and listed
in ``reference_state_declared_convention``); a row that crosses the phase
boundary is checked as before.  Without ``residence`` nothing changes."""

import inspect
import math

import numpy as np
import pytest

from rmgpy import constants
from rmgpy.quantity import Quantity
from rmgpy.reaction import Reaction
from rmgpy.rmg.polymer_input import (HybridPolymerReactor, PolymerPhase, PolymerPool, compile_polymer_phase,
                                     polymer_phase)

from solverPolymerTest import (_REV_KIN, _gate_pool_config, _refstate_rs, _sackur_tetrode_decades, _spc,
                               _trivial_nasa)

DECLARED = {"*": "melt", "H": "melt", "H2": "melt", "N2": "gas"}


def _phase(residence=None):
    pool = PolymerPool(label="PE", xs=2, monomer=_spc("C=C", "PE_repeat"), explicit_map={}, mu_species=[])
    return PolymerPhase(density=Quantity(850.0, "kg/m^3"), initial_moments={}, initial_explicit={},
                        pools=[pool], mass_transfer=[], residence=residence)


def _core():
    """The configured proxy, H, H2, a discrete gas-vetoed twin of the
    proxy's H-loss radical (the pe2 offender), the bath and a vetoed CH3."""
    from rmgpy.polymer import set_polymer_gas_veto
    proxy = _spc("CCCCCC", "PE(3)")
    h, h2 = _spc("[H]", "H(1)"), _spc("[H][H]", "H2(2)")
    twin = _spc("CC[CH]CCC", "CC[CH]CCC(7)")
    set_polymer_gas_veto(twin)
    n2 = _spc("N#N", "N2")
    n2.reactive = False
    ch3 = _spc("[CH3]", "CH3(17)")
    set_polymer_gas_veto(ch3)
    return [proxy, h, h2, twin, n2, ch3]


# --------------------------------------------------------------------------
# the mask
# --------------------------------------------------------------------------
def test_without_residence_the_mask_is_unchanged():
    mask = _phase().get_gas_mask(_core())
    # stage 1 alone: the configured proxy is condensed later, by the solver's pool-label pass
    assert [bool(x) for x in mask] == [True, True, True, True, True, True]


def test_melt_residence_condenses_the_named_and_the_default_species():
    mask = _phase(DECLARED).get_gas_mask(_core())
    # proxy, H, H2, the vetoed twin and CH3 condensed; the inert bath gas
    assert [bool(x) for x in mask] == [False, False, False, False, True, False]
    # named only: the rest as without residence
    mask = _phase({"H": "melt"}).get_gas_mask(_core())
    assert [bool(x) for x in mask] == [True, False, True, True, True, True]
    # an inert species is not covered by '*': it is named or it stays gas
    core = _core()
    mask = _phase({"*": "melt"}).get_gas_mask(core)
    assert bool(mask[4]) is True
    mask = _phase({"*": "melt", "N2": "melt"}).get_gas_mask(core)
    assert bool(mask[4]) is False


def test_moment_dummies_and_polymer_objects_are_left_to_the_pool_rules():
    from rmgpy.polymer import Polymer
    spawned = Polymer(label="PE_scission_tail", monomer="[CH2][CH2]", Mn=100.0, Mw=200.0)
    dummy = _spc("[Ne]", "PE_tail_mu0")
    dummy.is_moment_dummy = True
    core = [_spc("CCCCCC", "PE(3)"), spawned, dummy]
    today = [bool(x) for x in _phase().get_gas_mask(core)]
    assert today == [True, True, True]
    assert [bool(x) for x in _phase({"*": "melt"}).get_gas_mask(core)][1:] == today[1:]


def test_residence_refusals():
    with pytest.raises(ValueError, match="volatile.*later extension"):
        _phase({"H": "volatile"})
    with pytest.raises(ValueError, match="residence"):
        _phase({"H": "liquid"})
    with pytest.raises(ValueError, match="residence"):
        _phase({"": "melt"})
    with pytest.raises(ValueError, match="residence"):
        _phase(["H"])
    # a named gas the pool rules condense is a contradiction
    with pytest.raises(ValueError, match="PE.*gas"):
        _phase({"PE": "gas"}).get_gas_mask(_core())


def test_the_input_function_and_the_compiler_carry_residence():
    blueprint = polymer_phase("melt", [], "solvent", residence=dict(DECLARED))
    assert blueprint.residence == DECLARED
    with pytest.raises(ValueError, match="volatile"):
        polymer_phase("melt", [], "solvent", residence={"H": "volatile"})
    assert polymer_phase("melt", [], "solvent").residence is None
    compiled = compile_polymer_phase(blueprint, {}, {})
    assert compiled.residence == DECLARED


def _reactor(residence):
    a = _spc("CCCC", "A")
    phase = PolymerPhase(density=Quantity(1000.0, "kg/m^3"), initial_moments={}, initial_explicit={a: 1.0},
                         pools=[], mass_transfer=[], residence=residence)
    return a, HybridPolymerReactor(temperature=(800.0, "K"), pressure=(1.0e5, "Pa"), initialMoles={a: 1.0},
                                   polymerPhase=phase, terminationTime=(1.0, "s"))


def test_the_reactor_tells_the_solver_the_convention_is_declared():
    a, reactor = _reactor({"*": "melt"})
    assert reactor.to_solver_object([a], [], [], []).declared_melt_convention is True
    a, reactor = _reactor(None)
    assert reactor.to_solver_object([a], [], [], []).declared_melt_convention is False


# --------------------------------------------------------------------------
# the reference-state tripwire
# --------------------------------------------------------------------------
def _pool_core(**extra):
    """Pool proxy A (butane) at 0, its mu dummies 1-3, ``extra`` after, and a
    gas bath last (it seeds initial_mole_fractions)."""
    sp = {"A": _spc("CCCC", "A"), "A_mu0": _spc("CO", "A_mu0"), "A_mu1": _spc("C=O", "A_mu1"),
          "A_mu2": _spc("C#N", "A_mu2")}
    sp.update(extra)
    sp["bath"] = _spc("[Ar]", "bath")
    for s in sp.values():
        s.thermo = _trivial_nasa()
    return sp, list(sp.values())


def _build(core, rxns, gas, declared):
    mask = np.array([s.label in gas for s in core], dtype=bool)
    return _refstate_rs(core, rxns, mask, [_gate_pool_config()], {"A": (1.0, 5.0, 30.0)},
                        rs_kwargs={"declared_melt_convention": declared})


def test_an_all_melt_row_passes_under_the_declared_convention_and_stays_in_the_census():
    sp, core = _pool_core(G1=_spc("C", "G1"), G2=_spc("C=CC", "G2"))
    rxn = Reaction(reactants=[sp["A"]], products=[sp["G1"], sp["G2"]], **_REV_KIN)
    # liveness: the all-melt row's U is far above the refuse bound
    u = (_sackur_tetrode_decades(sp["G1"].molecule[0].get_molecular_weight(), 800.0)
         + _sackur_tetrode_decades(sp["G2"].molecule[0].get_molecular_weight(), 800.0)
         - _sackur_tetrode_decades(sp["A"].molecule[0].get_molecular_weight(), 800.0)
         + math.log10(1.0e5 / (constants.R * 800.0)))
    assert u > 3.0
    with pytest.raises(ValueError, match="unpaired reference-state"):
        _build(core, [rxn], {"bath"}, declared=False)
    rs = _build(core, [rxn], {"bath"}, declared=True)
    assert [str(r) for r, _u in rs.reference_state_declared_convention] == [str(rxn)]
    assert rs.reference_state_declared_convention[0][1] == pytest.approx(u, abs=1e-6)
    assert [name for name, _u in rs.reference_state_census] == [str(rxn)]


def test_a_row_across_the_phase_boundary_is_still_refused_under_the_declaration():
    sp, core = _pool_core(G1=_spc("C", "G1"), G2=_spc("C=CC", "G2"))
    rxn = Reaction(reactants=[sp["A"]], products=[sp["G1"], sp["G2"]], **_REV_KIN)
    with pytest.raises(ValueError, match="unpaired reference-state"):
        _build(core, [rxn], {"bath", "G1", "G2"}, declared=True)


def test_h_plus_r_to_rh_all_melt_is_not_u_zero_but_passes_under_the_declaration():
    """The design expected U = 0 once H, R and RH are all melt; it is not:
    H's translational term and dn = -1 stay unpaired (~7 decades here).  Only
    the declared convention, not the mask alone, lets the row through."""
    sp, core = _pool_core(H=_spc("[H]", "H"), R=_spc("C[CH][CH]C", "R"), RH=_spc("C[CH]CC", "RH"))
    rxn = Reaction(reactants=[sp["H"], sp["R"]], products=[sp["RH"]], **_REV_KIN)
    mw = {k: sp[k].molecule[0].get_molecular_weight() for k in ("H", "R", "RH")}
    u = abs(_sackur_tetrode_decades(mw["RH"], 800.0) - _sackur_tetrode_decades(mw["R"], 800.0)
            - _sackur_tetrode_decades(mw["H"], 800.0)) + math.log10(1.0e5 / (constants.R * 800.0))
    assert u > 3.0
    with pytest.raises(ValueError, match="unpaired reference-state"):
        _build(core, [rxn], {"bath"}, declared=False)
    rs = _build(core, [rxn], {"bath"}, declared=True)
    assert rs.reference_state_declared_convention[0][1] == pytest.approx(u, abs=1e-6)


def test_the_solver_flag_defaults_off():
    params = inspect.signature(__import__("rmgpy.solver.polymer", fromlist=["x"]).HybridPolymerSystem).parameters
    assert params["declared_melt_convention"].default is False


# --------------------------------------------------------------------------
# the volume counts what the mask condenses
# --------------------------------------------------------------------------
def _two_butanes(listed, residence=None):
    a, b = _spc("CCCC", "A"), _spc("CCCC", "B")
    phase = PolymerPhase(density=Quantity(1000.0, "kg/m^3"), initial_moments={},
                         initial_explicit={s: 1.0 for s in (a, b)[:listed]}, pools=[], mass_transfer=[],
                         residence={"*": "melt"} if residence is None else residence)
    reactor = HybridPolymerReactor(temperature=(800.0, "K"), pressure=(1.0e5, "Pa"), initialMoles={a: 1.0, b: 1.0},
                                   polymerPhase=phase, terminationTime=(1.0, "s"))
    return reactor.to_solver_object([a, b], [], [], [])


def test_a_melt_species_whose_mass_the_volume_omits_is_refused():
    """The review's case: two butanes, both declared melt with 1 mol each;
    listing one of them as a phase species halved V_poly (and doubled every
    concentration and bimolecular rate) while the gas masks agreed."""
    rs = _two_butanes(2)
    assert [bool(x) for x in rs.gas_species_mask] == [False, False]
    v_one = 1.0 * _spc("CCCC", "x").molecular_weight.value_si * constants.Na / 1000.0
    assert rs.V_poly == pytest.approx(2 * v_one, rel=1e-12)
    with pytest.raises(ValueError, match=r"'B'.*melt.*1 mol.*V_poly"):
        _two_butanes(1)


def test_a_melt_species_without_initial_moles_needs_no_listing():
    a, b = _spc("CCCC", "A"), _spc("CCCC", "B")
    phase = PolymerPhase(density=Quantity(1000.0, "kg/m^3"), initial_moments={}, initial_explicit={a: 1.0},
                         pools=[], mass_transfer=[], residence={"*": "melt"})
    reactor = HybridPolymerReactor(temperature=(800.0, "K"), pressure=(1.0e5, "Pa"), initialMoles={a: 1.0},
                                   polymerPhase=phase, terminationTime=(1.0, "s"))
    rs = reactor.to_solver_object([a, b], [], [], [])
    assert [bool(x) for x in rs.gas_species_mask] == [False, False]


def test_a_gas_species_with_initial_moles_needs_no_listing():
    rs = _two_butanes(1, residence={"A": "melt", "B": "gas"})
    assert [bool(x) for x in rs.gas_species_mask] == [False, True]


def test_listed_moles_that_differ_from_the_initial_moles_are_refused():
    a = _spc("CCCC", "A")
    phase = PolymerPhase(density=Quantity(1000.0, "kg/m^3"), initial_moments={}, initial_explicit={a: 0.5},
                         pools=[], mass_transfer=[], residence={"*": "melt"})
    reactor = HybridPolymerReactor(temperature=(800.0, "K"), pressure=(1.0e5, "Pa"), initialMoles={a: 1.0},
                                   polymerPhase=phase, terminationTime=(1.0, "s"))
    with pytest.raises(ValueError, match=r"'A'.*0.5 mol.*1 mol"):
        reactor.to_solver_object([a], [], [], [])
