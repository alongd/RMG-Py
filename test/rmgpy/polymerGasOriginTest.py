"""
Gas-origin provenance (I-113): a generated molecule whose heavy atoms ALL
come from genuine gas reactant molecules is gas -- never a Polymer, never a
pool daughter, never folded into a pool by any structural path.

The PE long deck stopped at 298 s on the row
``[CH2]CCC=C + PE <=> PE_scission_tail_2 + PE_mod``: gas 1-pentene (all heavy
atoms from the gas radical) was booked into a scission-tail pool by the wing
fold. These tests pin the provenance machinery that keeps it gas.
"""
import os

import pytest

from rmgpy import settings
from rmgpy.data.base import ForbiddenStructures
from rmgpy.data.rmg import RMGDatabase
from rmgpy.molecule import Molecule
from rmgpy.polymer import Polymer
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species

GAS_ORIGIN_KEY = "gas_origin"

POOLS = {
    "PE": dict(label="PE", monomer="[CH2][CH2]", end_groups=["[H]", "[H]"],
               cutoff=3, Mn=5000.0, Mw=6000.0, initial_mass=0.05),
    "PP": dict(label="polypropylene", monomer="[CH2][CH](C)",
               end_groups=["[H]", "[H]"], cutoff=3, Mn=1500.0, Mw=1800.0,
               initial_mass=0.1485),
    "PS": dict(label="PS", monomer="[CH2][CH]c1ccccc1",
               end_groups=["[H]", "[H]"], cutoff=3, Mn=5000.0, Mw=6000.0,
               initial_mass=0.05),
}

# Gas radicals whose H-abstraction product the wing fold books into a pool
# today (probed on 7ddc33280): the gas product must stay gas.
FOLD_CASES = {
    "PE": dict(radical="[CH2]CCC=C", gas_product="C=CCCC"),
    "PP": dict(radical="[CH2]CCC(C)C", gas_product="CCC(C)C"),
    "PS": dict(radical="C[CH]CCc1ccccc1", gas_product="CCCCc1ccccc1"),
}


def setup_module(module):
    db = RMGDatabase()
    path = os.path.join(settings["test_data.directory"], "testing_database")
    db.load_kinetics(os.path.join(path, "kinetics"),
                     kinetics_families=["H_Abstraction", "R_Recombination"],
                     reaction_libraries=[])
    for family in db.kinetics.families.values():
        family.forbidden = ForbiddenStructures()
    db.forbidden_structures = ForbiddenStructures()
    module.DB = db


def teardown_module(module):
    import rmgpy.data.rmg
    rmgpy.data.rmg.database = None


def _kinetics():
    return DB.kinetics


def _model(pool_name=None):
    cerm = CoreEdgeReactionModel()
    pool = None
    if pool_name is not None:
        pool, _ = cerm._register_polymer(Polymer(**POOLS[pool_name]),
                                         generate_thermo=False)
    return cerm, pool


def _gas(cerm, smiles, label):
    spc, _ = cerm.make_new_species(
        Species(molecule=[Molecule(smiles=smiles)], label=label),
        generate_thermo=False)
    return spc


def _mol(participant):
    return participant.molecule[0] if hasattr(participant, "molecule") else participant


def _all_molecules(rxn):
    out = []
    for reaction in (rxn, getattr(rxn, "reverse", None)):
        if reaction is None:
            continue
        for part in list(reaction.reactants) + list(reaction.products):
            mols = part.molecule if hasattr(part, "molecule") else [part]
            out.extend(mols)
    return out


################################################################################
# Items 1-3: genuine gas, stamp, verdict, strip
################################################################################

class TestGenuineGas:
    def test_classification(self):
        from rmgpy.gas_origin import is_genuine_gas_reactant

        gas = Species(molecule=[Molecule(smiles="[CH2]CCC=C")])
        assert is_genuine_gas_reactant(gas) is True
        assert is_genuine_gas_reactant(Polymer(**POOLS["PE"])) is False
        oligomer = Species(molecule=[Molecule(smiles="CCCCCCCCCC")])
        oligomer.props["explicit_dp_origin"] = ("PE", 5)
        assert is_genuine_gas_reactant(oligomer) is False
        dummy = Species(molecule=[Molecule(smiles="[Ne]")])
        dummy.is_moment_dummy = True
        assert is_genuine_gas_reactant(dummy) is False
        # the caller's condensed-membership verdict (solver phase mask) wins
        assert is_genuine_gas_reactant(gas, condensed=True) is False


class TestStampAndVerdict:
    def test_forward_rows_carry_aligned_verdicts_and_sources(self):
        """Gas radical at occurrence index 0 (stamp value 0: key presence,
        never truthiness), pool at index 1."""
        cerm, pe = _model("PE")
        rad = _gas(cerm, "[CH2]CCC=C", "pentenyl")
        rxns = _kinetics().generate_reactions_from_families([rad, pe])
        assert rxns, "liveness: H-abstraction rows expected"
        pentene = Molecule(smiles="C=CCCC")
        for rxn in rxns:
            assert rxn.is_forward
            assert rxn.gas_origin_side == "products"
            assert len(rxn.gas_origin_verdicts) == len(rxn.products)
            for product, verdict in zip(rxn.products, rxn.gas_origin_verdicts):
                assert verdict is _mol(product).is_isomorphic(pentene), str(rxn)
            for reactant, src in zip(rxn.reactants, rxn.gas_origin_sources):
                expected = 0 if _mol(reactant).is_isomorphic(rad.molecule[0]) else 1
                assert src == expected, str(rxn)

    def test_reverse_written_rows_put_verdicts_on_reactants(self):
        """R_Recombination is not its own reverse: a unimolecular gas
        supplied reactant yields reverse-written rows (generated radicals in
        ``reactants``)."""
        cerm, _ = _model()
        butane = _gas(cerm, "CCCC", "butane")
        rxns = _kinetics().generate_reactions_from_families([butane])
        rev = [r for r in rxns if not r.is_forward]
        assert rev, "liveness: reverse-written fission rows expected"
        for rxn in rev:
            assert rxn.gas_origin_side == "reactants"
            assert rxn.gas_origin_verdicts == tuple(True for _ in rxn.reactants)
            assert rxn.gas_origin_sources == (0,)

    def test_pool_fragments_are_not_gas(self):
        cerm, pe = _model("PE")
        rxns = _kinetics().generate_reactions_from_families([pe])
        rev = [r for r in rxns if not r.is_forward]
        assert rev, "liveness: proxy fission rows expected"
        for rxn in rev:
            assert rxn.gas_origin_side == "reactants"
            assert not any(rxn.gas_origin_verdicts)

    def test_explicit_dp_oligomer_reactant_is_not_stamped(self):
        """(c) The explicit-DP marker is what keeps the oligomer condensed:
        the same structure without it is genuine gas."""
        def octyl_verdicts(oligomer):
            cerm, _ = _model()
            methyl = _gas(cerm, "[CH3]", "methyl")
            rxns = _kinetics().generate_reactions_from_families([methyl, oligomer])
            methane = Molecule(smiles="C")
            out = set()
            for rxn in rxns:
                for product, verdict in zip(rxn.products, rxn.gas_origin_verdicts):
                    if _mol(product).is_isomorphic(methane):
                        assert verdict is True
                    else:
                        out.add(verdict)
            assert out, "liveness: octyl products expected"
            return out

        marked = Species(molecule=[Molecule(smiles="CCCCCCCC")], label="PE_dp4")
        marked.props["explicit_dp_origin"] = ("PE", 4)
        assert octyl_verdicts(marked) == {False}
        plain = Species(molecule=[Molecule(smiles="CCCCCCCC")], label="octane")
        assert octyl_verdicts(plain) == {True}

    def test_caller_condensed_flags_prevent_stamp(self):
        cerm, _ = _model()
        methyl = _gas(cerm, "[CH3]", "methyl")
        octane = _gas(cerm, "CCCCCCCC", "octane")
        rxns = _kinetics().generate_reactions_from_families(
            [methyl, octane], genuine_gas=(True, False))
        methane = Molecule(smiles="C")
        seen = False
        for rxn in rxns:
            for product, verdict in zip(rxn.products, rxn.gas_origin_verdicts):
                if not _mol(product).is_isomorphic(methane):
                    assert verdict is False
                    seen = True
        assert seen

    @pytest.mark.parametrize("supplied", ["forward", "reverse"])
    def test_no_stamp_survives_generation(self, supplied):
        """Stamps are stripped from both sides, every resonance structure and
        the reverse reaction kept by own-reverse families; the supplied
        originals are never stamped at all."""
        cerm, pe = _model("PE")
        if supplied == "forward":
            # H_Abstraction is its own reverse: every row keeps rxn.reverse
            spcs = [_gas(cerm, "C=C[CH2]", "allyl"), pe]
        else:
            spcs = [_gas(cerm, "CCCC", "butane")]
        rxns = _kinetics().generate_reactions_from_families(list(spcs))
        assert rxns
        if supplied == "forward":
            assert any(getattr(r, "reverse", None) is not None for r in rxns)
        for rxn in rxns:
            for mol in _all_molecules(rxn):
                for atom in mol.atoms:
                    assert GAS_ORIGIN_KEY not in atom.props, str(rxn)
        for spc in spcs:
            for mol in spc.molecule:
                for atom in mol.atoms:
                    assert GAS_ORIGIN_KEY not in atom.props
