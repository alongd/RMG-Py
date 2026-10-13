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
                     kinetics_families=["H_Abstraction", "R_Recombination",
                                        "Disproportionation"],
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
        rxns = _kinetics().generate_reactions_from_families(
            [rad, pe], only_families=["H_Abstraction"])
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
            # heavy-atom-free fragments ([H]) carry no provenance: False
            expected = tuple(any(not a.is_hydrogen() for a in _mol(r).atoms)
                             for r in rxn.reactants)
            assert rxn.gas_origin_verdicts == expected
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


################################################################################
# Items 4-5: origin-aware collapse, transport, census aggregation
################################################################################

def _constructed_row(products, verdicts):
    from rmgpy.data.kinetics.family import TemplateReaction

    rxn = TemplateReaction(
        reactants=[Molecule(smiles="CC[CH]CCC"), Molecule(smiles="CC[CH]CCC")],
        products=products, family="Disproportionation", template=["a", "b"],
        is_forward=True, degeneracy=1)
    rxn.gas_origin_side = "products"
    rxn.gas_origin_verdicts = tuple(verdicts)
    rxn.gas_origin_sources = (0, 1)
    return rxn


def _merge_pair(swap_verdicts, identical):
    hexene, hexane = Molecule(smiles="CCC=CCC"), Molecule(smiles="CCCCCC")
    hexene.assign_atom_ids()
    hexane.assign_atom_ids()
    first = _constructed_row([hexene, hexane], (True, False))
    if identical:
        products = [hexene.copy(deep=True), hexane.copy(deep=True)]
    else:
        products = [Molecule(smiles="CCC=CCC"), Molecule(smiles="CCCCCC")]
        for m in products:
            m.assign_atom_ids()
    second = _constructed_row(products, (False, True) if swap_verdicts else (True, False))
    return first, second


class TestOriginAwareCollapse:
    @pytest.mark.parametrize("identical", [True, False])
    def test_origin_ambiguous_merge_refuses_retained_row(self, identical):
        """(d) Two isomorphic rows whose (structure, verdict) multisets differ
        (gas hexene + pool hexane vs pool hexene + gas hexane) still merge,
        but the RETAINED representative is refused and the census counts it."""
        from rmgpy.data.kinetics.common import find_degenerate_reactions
        from rmgpy.gas_origin import ORIGIN_AMBIGUOUS_COLLAPSE, new_census

        first, second = _merge_pair(swap_verdicts=True, identical=identical)
        census = new_census()
        out = find_degenerate_reactions([first, second], same_reactants=0, census=census)
        assert len(out) == 1 and out[0] is first
        assert first.polymer_refused is True
        assert first.polymer_refused_accumulating is False
        assert first.polymer_refused_detail == ORIGIN_AMBIGUOUS_COLLAPSE
        assert census["collapse_refusals"] == 1

    @pytest.mark.parametrize("identical", [True, False])
    def test_same_provenance_merge_is_not_refused(self, identical):
        from rmgpy.data.kinetics.common import find_degenerate_reactions
        from rmgpy.gas_origin import new_census

        first, second = _merge_pair(swap_verdicts=False, identical=identical)
        census = new_census()
        out = find_degenerate_reactions([first, second], same_reactants=0, census=census)
        assert len(out) == 1
        assert out[0].polymer_refused is False
        assert out[0].polymer_refused_detail is None
        assert census["collapse_refusals"] == 0

    def _refused_row(self):
        from rmgpy.data.kinetics.common import find_degenerate_reactions

        first, second = _merge_pair(swap_verdicts=True, identical=False)
        out = find_degenerate_reactions([first, second], same_reactants=0)
        return out[0]

    def test_refusal_and_verdicts_survive_pickle(self):
        import pickle

        from rmgpy.gas_origin import ORIGIN_AMBIGUOUS_COLLAPSE

        row = self._refused_row()
        back = pickle.loads(pickle.dumps(row))
        assert back.polymer_refused is True
        assert back.polymer_refused_accumulating is False
        assert back.polymer_refused_detail == ORIGIN_AMBIGUOUS_COLLAPSE
        assert back.gas_origin_side == "products"
        assert back.gas_origin_verdicts == (True, False)
        assert back.gas_origin_sources == (0, 1)
        assert back.degeneracy == row.degeneracy

    def test_canonical_dedupe_preserves_refusal_detail_and_verdicts(self):
        """The canonical-dedupe merge (polymer.py) fills the canonical row's
        refusal detail and gas-origin metadata from the discarded candidate."""
        import pickle

        from rmgpy.data.kinetics.family import TemplateReaction
        from rmgpy.gas_origin import ORIGIN_AMBIGUOUS_COLLAPSE
        from rmgpy.polymer import merge_polymer_adjudication_stamps

        source = pickle.loads(pickle.dumps(self._refused_row()))
        canonical = TemplateReaction(reactants=list(source.reactants),
                                     products=list(source.products),
                                     family="Disproportionation", is_forward=True)
        merge_polymer_adjudication_stamps(source, canonical)
        assert canonical.polymer_refused is True
        assert canonical.polymer_refused_detail == ORIGIN_AMBIGUOUS_COLLAPSE
        assert canonical.gas_origin_side == "products"
        assert canonical.gas_origin_verdicts == (True, False)
        assert canonical.gas_origin_sources == (0, 1)

    def test_kinetics_flip_preserves_refusal_detail_and_verdicts(self, monkeypatch):
        import rmgpy.data.rmg as rmg_data
        from rmgpy.data.kinetics.family import TemplateReaction
        from rmgpy.gas_origin import ORIGIN_AMBIGUOUS_COLLAPSE
        from rmgpy.kinetics import Arrhenius

        pool = Polymer(**POOLS["PE"])
        gas = Species().from_smiles("CC")
        rxn = TemplateReaction(reactants=[pool], products=[gas],
                               pairs=[(pool, gas)], family="fake_family",
                               is_forward=True)
        rxn.polymer_refused = True
        rxn.polymer_refused_detail = ORIGIN_AMBIGUOUS_COLLAPSE
        rxn.gas_origin_side = "products"
        rxn.gas_origin_verdicts = (True,)
        rxn.gas_origin_sources = (0,)
        kin = Arrhenius(A=(1.0, "s^-1"), n=0.0, Ea=(0.0, "J/mol"), T0=(1.0, "K"))
        cerm = CoreEdgeReactionModel()
        monkeypatch.setattr(cerm, "generate_kinetics",
                            lambda reaction: (kin, "test", None, False), raising=False)

        class _FakeFamily:
            own_reverse = False

        class _FakeDb:
            families = {"fake_family": _FakeFamily()}

        monkeypatch.setattr(rmg_data, "get_db", lambda name: _FakeDb())
        cerm.apply_kinetics_to_reaction(rxn)
        assert rxn.reactants == [gas] and rxn.products == [pool]   # flipped
        assert rxn.polymer_refused is True
        assert rxn.polymer_refused_detail == ORIGIN_AMBIGUOUS_COLLAPSE
        assert rxn.gas_origin_side == "reactants"
        assert rxn.gas_origin_verdicts == (True,)
        assert rxn.gas_origin_sources == (0,)


def _pe_mod(cerm, pe):
    daughter = pe.create_reacted_copy(Molecule(smiles="CC[CH]CCC"), h_loss_feature=True)
    daughter, _ = cerm._register_polymer(daughter, generate_thermo=False)
    return daughter


class TestGeneratedCollapseAndProcnum:
    def test_generated_origin_ambiguous_collapse(self):
        """Gas 3-hexyl + the PE_mod pool (3-hexyl proxy): either side can
        donate the H, so isomorphic rows disagree on which hexane is gas."""
        from rmgpy.gas_origin import ORIGIN_AMBIGUOUS_COLLAPSE, new_census

        cerm, pe = _model("PE")
        pe_mod = _pe_mod(cerm, pe)
        hexyl = _gas(cerm, "CC[CH]CCC", "hexyl3")
        census = new_census()
        rxns = _kinetics().generate_reactions_from_families([hexyl, pe_mod], census=census)
        refused = [r for r in rxns if r.polymer_refused]
        assert census["collapse_refusals"] >= 1
        assert len(refused) >= 1
        assert all(r.polymer_refused_detail == ORIGIN_AMBIGUOUS_COLLAPSE for r in refused)

    def test_procnum_2_matches_procnum_1(self):
        """(e) Same rows, same verdicts, same census from the worker pool."""
        from rmgpy.gas_origin import new_census
        from rmgpy.rmg.react import react

        cerm, pe = _model("PE")
        pe_mod = _pe_mod(cerm, pe)
        tuples = [((_gas(cerm, "[CH2]CCC=C", "pentenyl"), pe),),
                  ((_gas(cerm, "[CH2]CCCCC", "hexyl1"), pe),),
                  ((_gas(cerm, "CC[CH]CCC", "hexyl3"), pe_mod),),
                  ((_gas(cerm, "CCCC", "butane"),),)]

        def run(procnum):
            census = new_census()
            out = react(tuples, procnum=procnum, census=census)
            rows = [sorted((str(r), r.gas_origin_side, r.gas_origin_verdicts,
                            r.gas_origin_sources, bool(r.polymer_refused),
                            r.polymer_refused_detail, r.degeneracy) for r in rxns)
                    for rxns in out]
            return rows, dict(census)

        rows1, census1 = run(1)
        rows2, census2 = run(2)
        assert rows1 == rows2
        assert census1 == census2
        assert census1["collapse_refusals"] >= 1
