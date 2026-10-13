"""
Proxy-identity collision: an H-abstraction by a gas radical from a chain whose
saturated co-product is isomorphic to the pool's own baseline proxy (hexane for
PE at cutoff 3, C9H20 for PP). The handshake folds that co-product into the
reactant pool as a BASELINE Polymer, so the product list is
[pool BASELINE, pool_mod FEATURE], which the classifier cannot place
(UNRESOLVED). It must be refused conduit-deferred at generation rather than
fall through to the legacy mu1-only flux, which fabricated PE_mod inventory.
"""
import logging

import pytest

from rmgpy.data.kinetics.family import _handshake_structures
from rmgpy.molecule import Molecule
from rmgpy.polymer import (Polymer, PolymerClass, PolymerFluxArchetype,
                           compute_h_loss_shape_evidence,
                           stamp_polymer_flux_archetype)
from rmgpy.reaction import Reaction
from rmgpy.species import Species

CASES = {
    "PE": dict(
        poly=dict(label="PE", monomer="[CH2][CH2]", end_groups=["[H]", "[H]"],
                  cutoff=3, Mn=5000.0, Mw=6000.0, initial_mass=0.05),
        radical="[CH2]CCCCC", proxy="CCCCCC", daughter="CC[CH]CCC"),
    "PP": dict(
        poly=dict(label="polypropylene", monomer="[CH2][CH](C)",
                  end_groups=["[H]", "[H]"], cutoff=3, Mn=1500.0, Mw=1800.0,
                  initial_mass=0.1485),
        radical="[CH2]CCC(C)CC(C)C", proxy="CCCC(C)CC(C)C",
        daughter="CCC[C](C)CC(C)C"),
}


def _stamped(case):
    pool = Polymer(**case["poly"])
    rad = Species(molecule=[Molecule(smiles=case["radical"])], label="rad")
    prods = [Molecule(smiles=case["proxy"]), Molecule(smiles=case["daughter"])]
    verdicts = compute_h_loss_shape_evidence([rad, pool], prods, [pool])
    _handshake_structures(prods, [pool], h_loss_verdicts=verdicts)
    rxn = Reaction(reactants=[rad, pool], products=prods)
    rxn.is_end_group_reaction = False
    stamp_polymer_flux_archetype(rxn, [rad, pool], [pool])
    return pool, rxn, prods


@pytest.mark.parametrize("name", sorted(CASES))
def test_proxy_collision_row_is_refused_conduit_deferred(name, caplog):
    with caplog.at_level(logging.WARNING):
        pool, rxn, prods = _stamped(CASES[name])
    # shape precondition: the collision really happened
    classes = [getattr(p, "_reacted_class", None) for p in prods]
    assert PolymerClass.BASELINE in classes
    assert rxn.polymer_flux_archetype == int(PolymerFluxArchetype.UNRESOLVED)
    assert rxn.polymer_refused is True
    assert rxn.polymer_refused_accumulating is False
    assert rxn.polymer_refused_detail == "proxy-identity-collision"
    assert any("proxy-identity-collision" in r.getMessage()
               for r in caplog.records)


def test_other_unresolved_shape_keeps_legacy_behaviour():
    """Control: UNRESOLVED without the BASELINE-fold-back + daughter shape is
    neither refused nor renamed."""
    case = CASES["PE"]
    pool = Polymer(**case["poly"])
    rad = Species(molecule=[Molecule(smiles=case["radical"])], label="rad")
    d1 = pool.create_reacted_copy(Molecule(smiles=case["daughter"]),
                                  h_loss_feature=True)
    d2 = pool.copy(deep=True)
    d2.label = "PE_other"
    d2._reacted_class = PolymerClass.FEATURE
    rxn = Reaction(reactants=[rad, pool], products=[d1, d2])
    rxn.is_end_group_reaction = False
    stamp_polymer_flux_archetype(rxn, [rad, pool], [pool])
    assert rxn.polymer_flux_archetype == int(PolymerFluxArchetype.UNRESOLVED)
    assert rxn.polymer_refused is False
    assert getattr(rxn, "polymer_refused_detail", None) is None
