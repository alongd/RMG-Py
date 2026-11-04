"""Red tests (I-302): one pool whose unit is not its proxy monomer's graph.

A linear statistical copolymer reaches the method of moments as one pool
(the kMC producer, ``rmgpy.kmc.mom_producer.copolymer_input`` on the kMC
line, seeds it from the declaration).  Its unit is one residue of mean molar
mass Mbar0 = sum x_i M_i with exact mean element counts; the proxy monomer
only names the pool.  polymer-mom's legacy ``polymer()`` takes the unit as
two numbers, ``unit_molar_mass`` (g/mol) and ``unit_atoms`` ({element:
count}, an exact rational given as "p/q"), so no declaration crosses to this
line.  Without them every pool is today's.

The numbers are the E/P test declaration's at x = (2/3, 1/3) (E = -CH2CH2-,
P = -CH2CH(CH3)-): Mbar0 = 2/3 M(C2H4) + 1/3 M(C3H6), mean residue C7/3
H14/3; its proxy unit is one E (E may follow itself), whose graph weighs
M(C2H4).
"""

from fractions import Fraction
from unittest.mock import MagicMock

import numpy as np
import pytest

from rmgpy.exceptions import InputError
from rmgpy.molecule import Molecule
from rmgpy.polymer import Polymer
from rmgpy.species import Species

PROXY = "[CH2][CH2]"
CAPS = ["[CH3]", "[CH3]"]
M_E = Molecule(smiles="C=C").get_molecular_weight() * 1000.0
M_P = Molecule(smiles="C=CC").get_molecular_weight() * 1000.0
MBAR = 2.0 / 3.0 * M_E + 1.0 / 3.0 * M_P
ATOMS = {"C": "7/3", "H": "14/3"}
EXACT = (("C", Fraction(7, 3)), ("H", Fraction(14, 3)))
MU = [1.0e-3, 8.0e-3, 7.0e-2]


def _deck_polymer(**kw):
    from rmgpy.rmg import input as rmg_input
    mock_rmg = MagicMock()
    mock_rmg.initial_species = []
    mock_rmg.reaction_model.iteration_num = 0
    mock_rmg.reaction_model.new_species_list = []
    mock_rmg.reaction_model.make_new_species.side_effect = lambda obj, **_k: (obj, True)
    old_rmg, old_sd = rmg_input.rmg, rmg_input.species_dict
    rmg_input.rmg, rmg_input.species_dict = mock_rmg, {}
    try:
        return rmg_input.polymer(**kw)
    finally:
        rmg_input.rmg, rmg_input.species_dict = old_rmg, old_sd


def _pool(**kw):
    args = dict(label="EP", monomer=PROXY, end_groups=CAPS, cutoff=4, moments=list(MU),
                initial_mass=MBAR * MU[1] / 1000.0, unit_molar_mass=MBAR, unit_atoms=ATOMS)
    args.update(kw)
    return Polymer(**args)


def test_the_deck_pool_weighs_its_declared_unit_not_its_proxy_graph():
    p = _deck_polymer(label="EP", monomer=PROXY, end_groups=CAPS, cutoff=4, moments=list(MU),
                      initial_mass=MBAR * MU[1] / 1000.0, unit_molar_mass=MBAR, unit_atoms=ATOMS)
    assert p.monomer_mw_g_mol == MBAR
    assert p.monomer_mw_g_mol != pytest.approx(M_E)
    assert p.unit_atoms == EXACT
    assert list(p.moments) == MU
    radical = _deck_polymer(label="EP_re", monomer=PROXY, end_groups=CAPS, cutoff=4, moments=[0.0, 0.0, 0.0],
                            initial_mass=0.0, unit_molar_mass=MBAR, unit_atoms=ATOMS,
                            end_radical_site="primary")
    assert list(radical.moments) == [0.0, 0.0, 0.0]
    assert radical.monomer_mw_g_mol == MBAR
    assert radical.unit_atoms == EXACT


def test_without_a_declared_unit_the_pool_is_todays():
    p = _pool(unit_molar_mass=None, unit_atoms=None)
    q = Polymer(label="EP", monomer=PROXY, end_groups=CAPS, cutoff=4, moments=list(MU),
                initial_mass=MBAR * MU[1] / 1000.0)
    assert p.monomer_mw_g_mol == q.monomer_mw_g_mol == pytest.approx(M_E)
    assert p.unit_molar_mass is None and p.unit_atoms is None
    assert p.fingerprint == q.fingerprint
    assert "_Unit-" not in q.fingerprint


def test_the_declared_unit_names_the_pool():
    """Two copolymer pools on one proxy unit (two compositions) are two pools."""
    a, b, plain = _pool(), _pool(unit_molar_mass=MBAR + 1.0), _pool(unit_molar_mass=None, unit_atoms=None)
    assert len({a.fingerprint, b.fingerprint, plain.fingerprint}) == 3


@pytest.mark.parametrize("kw, match", [
    (dict(unit_molar_mass=0.0), "unit_molar_mass"),
    (dict(unit_molar_mass=float("nan")), "unit_molar_mass"),
    (dict(unit_atoms={"C": "7/3", "H": "-1"}), "unit_atoms"),
    (dict(unit_atoms={"C": "seven"}), "unit_atoms"),
    (dict(unit_atoms={"Xx": 1}), "unit_atoms"),
    (dict(unit_atoms=None), "unit_atoms"),
    (dict(unit_molar_mass=None), "unit_molar_mass"),
])
def test_a_malformed_unit_is_refused_by_name(kw, match):
    args = dict(label="EP", monomer=PROXY, end_groups=CAPS, cutoff=4, moments=list(MU),
                initial_mass=MBAR * MU[1] / 1000.0, unit_molar_mass=MBAR, unit_atoms=ATOMS)
    args.update(kw)
    with pytest.raises(InputError, match=match):
        _deck_polymer(**args)


def test_the_phase_compiles_the_declared_unit_into_the_solver_pool():
    """The solver pool's unit is the declared one: mass and exact element
    counts, not the proxy graph's."""
    from rmgpy.rmg.polymer_input import PolymerPhaseBlueprint, compile_polymer_phase
    pool = _pool()
    dummies = [Species().from_smiles(s) for s in ("CO", "C=O", "C#N")]
    for k, d in enumerate(dummies):
        d.label = f"EP_mu{k}"
    species_dict = {pool.label: pool, **{d.label: d for d in dummies}}
    phase = compile_polymer_phase(PolymerPhaseBlueprint("EP", ["EP"], "EP"), {pool: MU[0]}, species_dict)
    assert np.allclose(np.asarray(phase.initial_moments["EP"], dtype=float), MU)
    # the melt's t = 0 volume weighs the declared unit too
    assert phase.calculate_volume() == pytest.approx(MU[1] * MBAR / 1000.0 / phase.density.value_si, rel=1e-12)
    config = phase.pools[0].to_config({s: i for i, s in enumerate([pool] + dummies)})
    assert config.monomer_mw_g_mol == MBAR
    assert config.repeat_unit_atoms == EXACT


def test_a_spawned_daughter_and_a_copy_keep_the_declared_unit():
    pool = _pool()
    for daughter in pool.generate_end_radical_daughters() + (pool.copy(deep=True),):
        assert daughter.monomer_mw_g_mol == MBAR
        assert daughter.unit_molar_mass == MBAR
        assert daughter.unit_atoms == EXACT
