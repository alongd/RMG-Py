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
import sys

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
    "PE": dict(pool="PE", radical="[CH2]CCC=C", gas_product="C=CCCC"),
    "PP": dict(pool="PP", radical="[CH2]CCC(C)C", gas_product="CCCC(C)C"),
    "PS": dict(pool="PS", radical="C[CH]CCc1ccccc1", gas_product="CCCCc1ccccc1"),
    # the proxy-identity-collision shape: gas hexane is the PE proxy's twin,
    # folded back onto PE (refused UNRESOLVED) before provenance existed
    "PE-collision": dict(pool="PE", radical="[CH2]CCCCC", gas_product="CCCCCC"),
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
        # a gas twin of the PE_mod proxy: structural registration would
        # resolve it onto the pool, so it is supplied unregistered
        hexyl = Species(molecule=[Molecule(smiles="CC[CH]CCC")], label="hexyl3")
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
                  ((Species(molecule=[Molecule(smiles="CC[CH]CCC")], label="hexyl3"), pe_mod),),
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


################################################################################
# Items 6-7: resolution, handshake bypass, spawn
################################################################################

def _veto(spc):
    from rmgpy.polymer import has_polymer_gas_veto
    return has_polymer_gas_veto(spc)


def _fold_rows(name):
    case = FOLD_CASES[name]
    cerm, pool = _model(case["pool"])
    rad = _gas(cerm, case["radical"], "rad")
    rxns = _kinetics().generate_reactions_from_families(
        [rad, pool], only_families=["H_Abstraction"])
    assert rxns
    return cerm, pool, rad, rxns


class TestGasProductResolution:
    @pytest.mark.parametrize("name", sorted(FOLD_CASES))
    def test_gas_product_never_folds_into_a_pool(self, name):
        """(a) The 298 s shape: the H-abstraction co-product made of gas
        atoms stays a gas Species; the chain radical lands in the pool's
        _mod daughter; the row is a live cross-pool VOLATILE_EJECTION whose
        eject units are the net non-polymer mass over the repeat-unit MW."""
        from rmgpy.polymer import PolymerFluxArchetype

        cerm, pool, rad, rxns = _fold_rows(name)
        gas_ref = Molecule(smiles=FOLD_CASES[name]["gas_product"])
        gas_ref_spc = Species(molecule=[gas_ref.copy(deep=True)])
        gas_ref_spc.generate_resonance_structures()
        for row in rxns:
            out, _ = cerm.make_new_reaction(row, check_existing=False, generate_thermo=False,
                                            generate_kinetics=False)
            assert out is not None
            gas = [p for p in out.products if not isinstance(p, Polymer)]
            polys = [p for p in out.products if isinstance(p, Polymer)]
            assert len(gas) == 1 and len(polys) == 1, str(out)
            assert gas_ref_spc.is_isomorphic(gas[0])
            assert _veto(gas[0]) and not gas[0].is_polymer_proxy
            assert polys[0].label.startswith(pool.label + "_mod"), polys[0].label
            assert out.polymer_flux_archetype == int(PolymerFluxArchetype.VOLATILE_EJECTION)
            assert out.polymer_refused is False
            net = gas_ref.get_molecular_weight() - rad.molecule[0].get_molecular_weight()
            assert out.polymer_eject_units == pytest.approx(
                net * 1000.0 / pool.monomer_mw_g_mol, rel=1e-9)
        assert cerm.gas_origin_census["gas_fold_bypassed"] >= 1
        assert cerm.gas_origin_census["rows_changed_class"] >= 1

    def test_supplied_participants_restored_by_source_index(self):
        """A supplied gas twin of the PE proxy (n-hexane) is restored by
        source index, never re-looked-up structurally onto the PE pool --
        even with the pool hot at the front of the MRU cache."""
        cerm = CoreEdgeReactionModel()
        twin = _gas(cerm, "CCCCCC", "hexane")
        pe, _ = cerm._register_polymer(Polymer(**POOLS["PE"]), generate_thermo=False)
        methyl = _gas(cerm, "[CH3]", "methyl")
        cerm.species_cache.insert(0, pe)
        cerm.species_cache.pop()
        rxns = _kinetics().generate_reactions_from_families(
            [methyl, twin], only_families=["H_Abstraction"])
        assert rxns
        for row in rxns:
            out, _ = cerm.make_new_reaction(row, check_existing=False, generate_thermo=False,
                                            generate_kinetics=False,
                                            supplied_species=(methyl, twin))
            assert any(r is twin for r in out.reactants), str(out)
            assert not any(isinstance(s, Polymer) for s in out.reactants + out.products)

    @pytest.mark.parametrize("order", ["pool_first", "gas_first"])
    @pytest.mark.parametrize("cache", ["cold", "warm"])
    @pytest.mark.parametrize("pool_name,twin_smiles,rad_a,rad_b", [
        ("PP", "CCCC(C)CC(C)C", "[CH2]C(C)C", "C[CH]CCC"),   # isobutyl + 2-pentyl
        ("PE", "CCCCCC", "[CH2]CC", "[CH2]CC"),              # n-propyl + n-propyl
    ])
    def test_gas_only_recombination_yields_gas_species(self, order, cache, pool_name,
                                                       twin_smiles, rad_a, rad_b):
        """(b) Two gas radicals recombine to the pool's 3-mer proxy
        structure. Both radicals are gas, so the product is a gas Species
        whatever the MRU cache holds and whichever twin was registered
        first."""
        cerm = CoreEdgeReactionModel()
        twin = _gas(cerm, twin_smiles, "twin") if order == "gas_first" else None
        pp, _ = cerm._register_polymer(Polymer(**POOLS[pool_name]), generate_thermo=False)
        ib = _gas(cerm, rad_a, "rad_a")
        pt = ib if rad_b == rad_a else _gas(cerm, rad_b, "rad_b")
        if cache == "warm":
            cerm.species_cache.insert(0, pp)
            cerm.species_cache.pop()
        else:
            cerm.species_cache = [None] * 4
        proxy = pp.molecule[0]
        rows = [r for r in _kinetics().generate_reactions_from_families(
                    [ib, pt], only_families=["R_Recombination"])
                if r.is_forward and len(r.products) == 1
                and r.products[0].is_isomorphic(proxy)]
        assert rows, "liveness: proxy-isomorphic recombination row expected"
        for row in rows:
            out, _ = cerm.make_new_reaction(row, check_existing=False, generate_thermo=False,
                                            generate_kinetics=False,
                                            supplied_species=(ib, pt))
            product = out.products[0]
            assert not isinstance(product, Polymer)
            assert _veto(product) and not product.is_polymer_proxy
            if twin is not None:
                assert product is twin
            assert out.polymer_refused is False
        if order == "pool_first" and cache == "warm":
            assert cerm.gas_origin_census["gas_lookup_unpooled"] >= 1

    def test_decay_registration_inherits_gas_verdict(self, monkeypatch):
        import rmgpy.rmg.model as model_module

        cerm = CoreEdgeReactionModel()

        def fake_decay(spc):
            return [Species(molecule=[Molecule(smiles="C=C")]),
                    Species(molecule=[Molecule(smiles="CC")])]

        monkeypatch.setattr(model_module, "decay_species", fake_decay)
        mol = Molecule(smiles="CCCC")
        mol.is_polymer_proxy = True
        out = cerm.make_new_species(mol, generate_thermo=False, check_decay=True,
                                    gas_origin=True)
        assert isinstance(out, list) and len(out) == 2
        for spc, _new in out:
            assert _veto(spc)
            assert not spc.is_polymer_proxy


class TestHandshakeBypass:
    def test_gas_verdict_bypasses_fold(self):
        from rmgpy.data.kinetics.family import _handshake_structures
        from rmgpy.polymer import compute_h_loss_shape_evidence

        pool = Polymer(**POOLS["PE"])
        rad = Species(molecule=[Molecule(smiles="[CH2]CCCCC")], label="rad")
        prods = [Molecule(smiles="CCCCCC"), Molecule(smiles="CC[CH]CCC")]
        verdicts = compute_h_loss_shape_evidence([rad, pool], prods, [pool])
        _handshake_structures(prods, [pool], h_loss_verdicts=verdicts,
                              gas_verdicts=[True, False])
        assert isinstance(prods[0], Molecule)
        assert _veto(prods[0]) and not prods[0].is_polymer_proxy
        assert isinstance(prods[1], Polymer) and prods[1].label.startswith("PE_mod")

    def test_gas_verdict_bypasses_h_loss_exemption(self):
        """An H-loss-daughter-shaped molecule made of gas atoms is gas: it is
        neither folded into PE_mod nor exempted from the durable veto."""
        from rmgpy.data.kinetics.family import _handshake_structures

        pool = Polymer(**POOLS["PE"])
        prods = [Molecule(smiles="CC[CH]CCC"), Molecule(smiles="CCCCCC")]
        _handshake_structures(prods, [pool], h_loss_verdicts=[True, False],
                              gas_verdicts=[True, False])
        assert isinstance(prods[0], Molecule)
        assert _veto(prods[0])


class TestSpawnIgnoresVetoedGas:
    def _vetoed_proxy_hexane(self):
        from rmgpy.polymer import set_polymer_gas_veto

        cand = Species(molecule=[Molecule(smiles="CCCCCC")], label="hexane")
        cand.is_polymer_proxy = True
        set_polymer_gas_veto(cand)
        return cand

    def test_multipool_classifier_skips_vetoed_candidate(self):
        """(f) polymer.py: a vetoed gas species is never a spawn/proxy candidate."""
        from rmgpy.polymer import process_polymer_candidates_multipool

        pe = Polymer(**POOLS["PE"])
        cand = self._vetoed_proxy_hexane()
        processed, intents = process_polymer_candidates_multipool(
            candidates=[cand], reaction_model=None, pool_registry=[pe])
        assert processed == [] and intents == []

    def test_model_spawn_pass_skips_vetoed_candidate(self, monkeypatch):
        """(f) model.py: the iteration-boundary spawn pass filters the veto
        before classification."""
        import rmgpy.polymer as polymer_module

        seen = []

        def record(candidates, **kwargs):
            seen.extend(candidates)
            return [], []

        monkeypatch.setattr(polymer_module, "process_polymer_candidates_multipool", record)
        cerm, _pe = _model("PE")
        cerm._apply_multipool_spawn_pass([self._vetoed_proxy_hexane()])
        assert seen == []


class TestThe298sRowSolverDeltas:
    @pytest.mark.parametrize("direction", ["forward", "reverse"])
    def test_moment_and_element_deltas(self, direction):
        """(a) pentenyl + PE -> gas 1-pentene + PE_mod through the real
        solver residual, in the generated and in the kinetics-flipped
        direction. Per event: the gas side gains (forward) or loses
        (reverse) exactly one H; the source pool loses one whole chain
        bundle and the destination gains the SAME bundle (a_moment = 0:
        the H mass rides on PE_mod's per-chain defect); C and H close."""
        import numpy as np

        from rmgpy.kinetics import Arrhenius
        from rmgpy.polymer import PolymerFluxArchetype, restamp_flipped_polymer_archetype
        from rmgpy.reaction import Reaction
        from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig

        cerm, pe, rad, rxns = _fold_rows("PE")
        out, _ = cerm.make_new_reaction(rxns[0], check_existing=False, generate_thermo=False,
                                        generate_kinetics=False, supplied_species=(rad, pe))
        gas = next(p for p in out.products if not isinstance(p, Polymer))
        mod = next(p for p in out.products if isinstance(p, Polymer))
        assert mod.chain_mass_defect_g_mol == pytest.approx(1.00794, abs=1e-3)
        reactants, products = list(out.reactants), list(out.products)
        if direction == "reverse":
            flipped = Reaction(reactants=products, products=reactants)
            flipped.polymer_flux_archetype = out.polymer_flux_archetype
            restamp_flipped_polymer_archetype(flipped)
            src = flipped
            reactants, products = products, reactants
        else:
            src = out
        assert src.polymer_flux_archetype == int(PolymerFluxArchetype.VOLATILE_EJECTION)
        assert src.polymer_refused is False

        rxn = Reaction(reactants=reactants, products=products,
                       kinetics=Arrhenius(A=(1.0e3, "m^3/(mol*s)"), n=0.0,
                                          Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
                       reversible=False)
        rxn.is_end_group_reaction = bool(src.is_end_group_reaction)
        rxn.polymer_flux_archetype = src.polymer_flux_archetype
        rxn.polymer_eject_units = src.polymer_eject_units

        def dummies(label):
            out = []
            for k, smi in enumerate(("CO", "C=O", "C#N")):
                s = Species(molecule=[Molecule(smiles=smi)], label="%s_mu%d" % (label, k))
                s.is_moment_dummy = True
                out.append(s)
            return out

        core = [pe] + dummies(pe.label) + [mod] + dummies(mod.label) + [rad, gas]
        mask = np.array([False] * 8 + [True, True], dtype=bool)
        pools = [PolymerPoolConfig(label=pe.label, xs=3, explicit_dp_to_species_index={},
                                   mu_indices=(1, 2, 3), monomer_poly_index=None,
                                   k_scission=0.0, k_unzip=0.0, tail_kinetics=None,
                                   monomer_mw_g_mol=pe.monomer_mw_g_mol),
                 PolymerPoolConfig(label=mod.label, xs=3, explicit_dp_to_species_index={},
                                   mu_indices=(5, 6, 7), monomer_poly_index=None,
                                   k_scission=0.0, k_unzip=0.0, tail_kinetics=None,
                                   monomer_mw_g_mol=pe.monomer_mw_g_mol,
                                   chain_mass_defect_g_mol=mod.chain_mass_defect_g_mol)]
        moments = {pe.label: (1.0, 50.0, 3000.0), mod.label: (1.0, 40.0, 2000.0)}
        feed = rad if direction == "forward" else gas
        rs = HybridPolymerSystem(
            T=1000.0, P=1.0e5, initial_mole_fractions={feed: 2.0}, V_poly=1.0,
            polymer_pools=pools, mass_transfer=[], gas_species_mask=mask.copy(),
            constant_gas_volume=False, initial_polymer_moments=moments, termination=[])
        rs.initialize_model(core, [rxn], [], [])
        dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
        assert rs.reaction_refused[0] == 0

        sign = 1.0 if direction == "forward" else -1.0
        r = sign * dn[9]           # events per second (mol/s) toward the alkene
        assert r > 0.0
        assert dn[8] == pytest.approx(-sign * r, rel=1e-12)
        # one length-biased chain per event (mid-chain site): the bundle is
        # (1, mu2/mu1, mu3/mu1) of the DEBITED pool, mu3 from the solver's
        # closure; it lands in the other pool unchanged (a_moment = 0)
        from rmgpy.solver.polymer import _safe_mu3_from_mu012
        assert not rxn.is_end_group_reaction
        mu0, mu1, mu2 = moments[pe.label] if direction == "forward" else moments[mod.label]
        bundle = (1.0, mu2 / mu1, _safe_mu3_from_mu012(mu0, mu1, mu2) / mu1)
        d_pe = dn[1:4]
        d_mod = dn[5:8]
        np.testing.assert_allclose(d_pe, [-sign * r * b for b in bundle], rtol=1e-9)
        np.testing.assert_allclose(d_mod, [sign * r * b for b in bundle], rtol=1e-9)
        # elemental deltas: gas C5H9 -> C5H10; chains C(2n)H(4n+2), PE_mod one H short
        c_gas = 5.0 * (dn[8] + dn[9])
        h_gas = 9.0 * dn[8] + 10.0 * dn[9]
        c_pool = 2.0 * (d_pe[1] + d_mod[1])
        h_pool = (4.0 * d_pe[1] + 2.0 * d_pe[0]) + (4.0 * d_mod[1] + 2.0 * d_mod[0] - d_mod[0])
        assert c_gas == pytest.approx(0.0, abs=1e-12 * r)
        assert c_pool == pytest.approx(0.0, abs=1e-9 * r)
        assert h_gas == pytest.approx(sign * r, rel=1e-12)
        assert h_pool == pytest.approx(-sign * r, rel=1e-9)


class TestCondensedMembershipPlumbing:
    def test_phase_mask_folds_into_genuine_gas_flags(self):
        """The solver's phase mask (False = condensed) decides for plain
        species; Polymers are never genuine gas whatever the mask says."""
        import numpy as np

        from rmgpy.rmg.main import make_condensed_membership

        cerm, pe = _model("PE")
        plain_condensed = _gas(cerm, "CC[CH]CCCCC", "pe_hloss")
        gas = _gas(cerm, "[CH3]", "methyl")

        class _Phase:
            def get_gas_mask(self, species):
                return np.array([s is not plain_condensed for s in species])

        class _System:
            polymerPhase = _Phase()

        assert make_condensed_membership([object()]) is None
        cerm.condensed_membership = make_condensed_membership([_System()])
        assert cerm.genuine_gas_flags([pe, plain_condensed, gas]) == [False, False, True]

    def test_react_all_threads_flags_per_tuple(self):
        """A species the phase mask condenses is not stamped: its fragments
        carry no gas verdict, in serial and parallel generation alike."""
        import numpy as np

        from rmgpy.gas_origin import new_census
        from rmgpy.rmg.react import react_all

        cerm, _ = _model()
        methyl = _gas(cerm, "[CH3]", "methyl")
        octane = _gas(cerm, "CCCCCCCC", "octane")
        core = [methyl, octane]
        uni = np.zeros(2, bool)
        bi = np.zeros((2, 2), bool)
        bi[0, 1] = True
        methane = Molecule(smiles="C")
        for procnum in (1, 2):
            lists, tuples = react_all(core, 2, uni, bi, procnum=procnum,
                                      genuine_gas=[True, False], census=new_census())
            seen = False
            for rxns in lists:
                for rxn in rxns:
                    for p, v in zip(getattr(rxn, rxn.gas_origin_side), rxn.gas_origin_verdicts):
                        if not _mol(p).is_isomorphic(methane):
                            assert v is False
                            seen = True
            assert seen


################################################################################
# Items 8-9: twins downstream, identity-filter counter
################################################################################

def _twin_core():
    cerm = CoreEdgeReactionModel()
    pe, _ = cerm._register_polymer(Polymer(**POOLS["PE"]), generate_thermo=False)
    twin = Species(molecule=[Molecule(smiles="CCCCCC")], label="hexane")
    twin.index = 7
    return pe, twin


class TestTwinTieBreak:
    def test_restart_reorder_keeps_gas_and_polymer_twins_apart(self):
        """(g) A restart map listing gas n-hexane before PE (both hexane
        adjacency lists) must reorder the live core by role, not by the
        first structural match."""
        from rmgpy.rmg.main import reorder_restart_core_species, restart_species_role

        pe, twin = _twin_core()
        assert restart_species_role(pe) == "polymer"
        assert restart_species_role(twin) == "gas"
        restart = [Species().from_adjacency_list(s.molecule[0].to_adjacency_list())
                   for s in (twin, pe)]
        core = [pe, twin]
        reordered = reorder_restart_core_species(restart, ["gas", "polymer"], core)
        assert reordered[0] is twin and reordered[1] is pe
        assert core == []

    def test_restart_map_round_trip_writes_roles(self, tmp_path):
        from rmgpy.rmg.main import load_restart_species_roles, write_restart_species_map

        pe, twin = _twin_core()
        write_restart_species_map(str(tmp_path), [twin, pe])
        roles = load_restart_species_roles(str(tmp_path / "species_map.yml"), 2)
        assert roles == ["gas", "polymer"]
        # restart files written before the roles existed: legacy first match
        (tmp_path / "legacy").mkdir()
        assert load_restart_species_roles(str(tmp_path / "legacy" / "species_map.yml"), 2) is None
        # a roles file inconsistent with the map is ignored
        assert load_restart_species_roles(str(tmp_path / "species_map.yml"), 3) is None

    def test_cantera_fallback_label_prefers_same_role(self):
        """(g) Cantera's structural fallback must not swap twins."""
        from rmgpy.cantera import get_label

        pe, twin = _twin_core()
        species_list = [pe, twin]
        gas_copy = Species(molecule=[Molecule(smiles="CCCCCC")], label="other")
        assert get_label(gas_copy, species_list) == "hexane(7)"
        poly_copy = pe.copy(deep=True)
        assert get_label(poly_copy, species_list).startswith(pe.label)
        species_list = [twin, pe]
        assert get_label(poly_copy, species_list).startswith(pe.label)
        assert get_label(gas_copy, species_list) == "hexane(7)"


class TestIdentityFilterCounter:
    def test_gas_hexyl_plus_pe_identity_drop_is_counted(self):
        """Item 9 (measured, not fixed): gas 1-hexyl + PE <=> gas hexane +
        chain 1-hexyl has identical participant multisets and is dropped by
        the identity filter; the census counts it as a gas/pool drop."""
        from rmgpy.gas_origin import new_census

        cerm, pe = _model("PE")
        hexyl = _gas(cerm, "[CH2]CCCCC", "hexyl1")
        census = new_census()
        _kinetics().generate_reactions_from_families(
            [hexyl, pe], only_families=["H_Abstraction"], census=census)
        assert census["identity_filter_drops_gas_mixed"] >= 1
        assert census["identity_filter_drops"] >= census["identity_filter_drops_gas_mixed"]

    def test_pure_gas_identity_drop_is_not_gas_mixed(self):
        from rmgpy.gas_origin import new_census

        cerm, _ = _model()
        census = new_census()
        _kinetics().generate_reactions_from_families(
            [_gas(cerm, "[CH3]", "methyl"), _gas(cerm, "C", "methane")],
            only_families=["H_Abstraction"], census=census)
        assert census["identity_filter_drops"] >= 1
        assert census["identity_filter_drops_gas_mixed"] == 0


class TestRestoredPoolRowsAreNeverUnstamped:
    def test_reverse_written_pool_rows_are_classified_or_refused(self):
        """Restoring a supplied pool by source index puts the Polymer on the
        products side of its reverse-written rows (PE_mod beta-scission:
        H + hexene <=> PE_mod), where no handshake or archetype stamp runs.
        The structural lookup used to hand such rows to a gas twin; now they
        reach the pool, so every one must be stamped or refused -- the
        solver's unstamped-proxy-row invariant -- and counted."""
        import rmgpy.gas_origin as gas_origin
        RESTORED_POOL_UNCLASSIFIED = getattr(gas_origin, "RESTORED_POOL_UNCLASSIFIED",
                                             "restored-pool-unclassified")

        db = RMGDatabase()
        path = os.path.join(settings["test_data.directory"], "testing_database")
        db.load_kinetics(os.path.join(path, "kinetics"),
                         kinetics_families=["R_Addition_MultipleBond"], reaction_libraries=[])
        for family in db.kinetics.families.values():
            family.forbidden = ForbiddenStructures()
        try:
            cerm, pe = _model("PE")
            pe_mod = _pe_mod(cerm, pe)
            rows = db.kinetics.generate_reactions_from_families([pe_mod])
            rev = [r for r in rows if not r.is_forward]
            assert rev, "liveness: reverse-written beta-scission rows expected"
            refused = 0
            for row in rows:
                out, _ = cerm.make_new_reaction(row, check_existing=False, generate_thermo=False,
                                                generate_kinetics=False, supplied_species=(pe_mod,))
                if out is None:
                    continue
                parts = out.reactants + out.products
                if any(isinstance(s, Polymer) for s in parts):
                    assert out.polymer_flux_archetype != 0 or out.polymer_refused, str(out)
                if out.polymer_refused_detail == RESTORED_POOL_UNCLASSIFIED:
                    refused += 1
            assert refused >= 1
            assert cerm.gas_origin_census["restored_pool_refusals"] == refused
        finally:
            setup_module(sys.modules[__name__])
