"""Fast PS chain-end regressions; no database or compiled artifact required."""

from itertools import combinations
from types import SimpleNamespace

import pytest

import rmgpy.kmc.compiler as compiler
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.kinetics.arrhenius import Arrhenius
from rmgpy.molecule.molecule import Bond, Molecule
from rmgpy.species import Species
from stateTest import ps_artifact


BENZYLIC_DECLARATIONS = {
    "benzylic_end_radical": ("benzylic_end_radical",),
    "benzylic_end_radical+styrene": ("benzylic_end_radical", "styrene"),
    "benzylic_end_radical+pristine": ("benzylic_end_radical", "pristine"),
    "benzylic_end_radical+end_radical": ("benzylic_end_radical", "end_radical"),
    "benzylic_end_radical+benzylic_end_radical": ("benzylic_end_radical",) * 2,
    "junction_radical+benzylic_end_radical": ("junction_radical", "benzylic_end_radical"),
}


@pytest.mark.parametrize("units", [1, 2, 3])
def test_both_end_orientations_have_distinct_exact_graphs(units):
    primary = Molecule(smiles=compiler.linear_ps_smiles(units, radical_end=True))
    benzylic = Molecule(smiles=compiler.benzylic_ps_end_smiles(units))
    expected = Molecule(smiles="CC(c1ccccc1)" * (units - 1) + "C[CH](c1ccccc1)")
    assert benzylic.is_isomorphic(expected)
    assert not primary.is_isomorphic(benzylic)
    assert primary.get_element_count() == benzylic.get_element_count() == {
        "C": 8 * units, "H": 8 * units + 1
    }
    radical = next(atom for atom in benzylic.atoms if atom.radical_electrons)
    assert sum(atom.element.number == 1 for atom in radical.edges) == 1
    assert sum(atom.element.number == 6 for atom in radical.edges) == 2
    assert any(benzylic.is_atom_in_cycle(atom) for atom in radical.edges)


def test_end_builder_rejects_invalid_and_incompatible_features():
    with pytest.raises(ValueError, match="repeat unit"):
        compiler.benzylic_ps_end_smiles(0)
    with pytest.raises(ValueError, match="multiple"):
        compiler.linear_ps_smiles(3, radical_unit=1, radical_units=(2,))
    with pytest.raises(ValueError, match="end"):
        compiler.linear_ps_smiles(3, radical_end=True, radical_unit=1)


@pytest.mark.parametrize("units, expected_count", [(1, 4), (2, 11), (3, 22)])
def test_catalogue_covers_every_hydrogen_deletion_graph(units, expected_count):
    catalogue = compiler.short_ps_molecule_catalogue(1)
    assert catalogue["size"] == catalogue["termination_bound"] == 37
    stored = [Molecule().from_adjacency_list(row["adjacency_list"])
              for row in catalogue["molecules"] if row["repeat_units"] == units]
    assert len(stored) == expected_count
    # Independent explicit-H deletion enumeration, no reflection quotient.
    base = Molecule(smiles="CC(c1ccccc1)" * units)
    expected = []
    for count in (0, 1, 2):
        for sites in combinations(range(2 * units), count):
            graph = base.copy(deep=True)
            backbone = [a for a in graph.atoms if a.element.number == 6
                        and not graph.is_atom_in_cycle(a)]
            for site in sites:
                atom = backbone[site]
                graph.remove_atom(next(a for a in atom.edges if a.element.number == 1))
                atom.radical_electrons += 1
            graph.update(sort_atoms=False)
            assert sum(graph.is_isomorphic(other) for other in stored) == 1
            assert not any(graph.is_isomorphic(other) for other in expected)
            assert graph.get_element_count() == {"C": 8 * units, "H": 8 * units + 2 - count}
            expected.append(graph)
    assert len(expected) == expected_count
    assert catalogue == compiler.short_ps_molecule_catalogue(1)


def test_all_six_benzylic_declarations_schedule_all_families_in_both_contexts():
    proxies = compiler.ps_proxy_set()
    assert len(proxies) == 34
    by_name = {proxy.site_type: proxy for proxy in proxies}
    for suffix, units in [("", 3), ("@5", 5)]:
        for name, labels in BENZYLIC_DECLARATIONS.items():
            proxy = by_name[name + suffix]
            assert proxy.participant_site_types == labels
            assert proxy.metadata["proxy_units"] == units
            assert set(proxy.metadata["family_candidates"]) == {
                "Disproportionation", "H_Abstraction", "R_Addition_MultipleBond",
                "R_Recombination", "intra_H_migration"
            }
            for species, label in zip(proxy.reactants, labels):
                if label == "benzylic_end_radical":
                    size = units - 1 if name.endswith("+styrene") else units
                    expected = Molecule(smiles="CC(c1ccccc1)" * (size - 1) + "C[CH](c1ccccc1)")
                    assert species.molecule[0].is_isomorphic(expected)
        # Retain the legacy primary declaration, including 5→6 styrene input.
        primary = by_name["end_radical+styrene" + suffix]
        size = 2 if units == 3 else 5
        assert primary.reactants[0].molecule[0].is_isomorphic(
            Molecule(smiles="[CH2]C(c1ccccc1)" + "CC(c1ccccc1)" * (size - 1)))


def test_benzylic_discovery_calls_each_public_candidate_pipeline():
    calls = []
    class EmptyDatabase:
        families = {family: object() for family in compiler.PS_FAMILY_CANDIDATES}

        def generate_reactions_from_families(self, reactants, only_families, resonance):
            calls.append((only_families, resonance))
            return []

    proxies = [proxy for proxy in compiler.ps_proxy_set() if "benzylic_end_radical" in proxy.site_type]
    compiler.EventSetCompiler.discover_family_reactions(EmptyDatabase(), proxies)
    assert len(calls) == 60  # six declarations × two contexts × five families
    assert all(resonance for _, resonance in calls)
    for family in compiler.PS_FAMILY_CANDIDATES:
        assert sum(selected == [family] for selected, _ in calls) == 12


@pytest.mark.parametrize("smiles, label", [
    ("[CH2]C(c1ccccc1)CC(c1ccccc1)", "end_radical"),
    ("CC(c1ccccc1)C[CH](c1ccccc1)", "benzylic_end_radical"),
    ("CC(c1ccccc1)[CH]C(c1ccccc1)", "interior_radical"),
    ("C[C](c1ccccc1)CC(c1ccccc1)", "interior_radical"),
    ("CC(c1cc[c]cc1)", "junction_radical"),
])
def test_reverse_site_labels_follow_backbone_connectivity(smiles, label):
    instance = compiler.EventSetCompiler(None, [], [])
    proxy = compiler.SiteProxy("other", [Species(molecule=[Molecule(smiles="C=Cc1ccccc1")])])
    reverse = SimpleNamespace(reactants=[Species(molecule=[Molecule(smiles=smiles)])])
    assert instance._direction_proxy(proxy, reverse).participant_site_types == (label,)


def mapped_addition(units, *, primary=False, opposite=False):
    """Known atom edits, independently specified without family generation."""
    chain = Molecule(smiles=("[CH2]C(c1ccccc1)" + "CC(c1ccccc1)" * (units - 1)
                            if primary else "CC(c1ccccc1)" * (units - 1) + "C[CH](c1ccccc1)"))
    styrene = Molecule(smiles="C=Cc1ccccc1")
    for identifier, atom in enumerate(chain.atoms + styrene.atoms, 1):
        atom.id = identifier
    radical = next(atom for atom in chain.atoms if atom.radical_electrons)
    alkene = next(bond for bond in styrene.get_all_edges()
                  if bond.is_double() and not styrene.is_atom_in_cycle(bond.atom1)
                  and not styrene.is_atom_in_cycle(bond.atom2))
    tail = next(atom for atom in (alkene.atom1, alkene.atom2)
                if sum(a.element.number == 1 for a in atom.edges) == 2)
    head = alkene.atom2 if tail is alkene.atom1 else alkene.atom1
    attacker = head if opposite else tail
    new_radical = tail if opposite else head
    # Legacy primary ceiling products preserve the primary end (head attack).
    if primary:
        attacker, new_radical = head, tail
    product = chain.copy(deep=True).merge(styrene.copy(deep=True))
    by_id = {atom.id: atom for atom in product.atoms}
    product.add_bond(Bond(by_id[radical.id], by_id[attacker.id], 1))
    product.get_bond(by_id[head.id], by_id[tail.id]).order = 1
    by_id[radical.id].radical_electrons = 0
    by_id[new_radical.id].radical_electrons = 1
    product.update(sort_atoms=False)
    reaction = TemplateReaction(
        reactants=[Species(molecule=[chain]), Species(molecule=[styrene])],
        products=[Species(molecule=[product])],
        kinetics=Arrhenius(A=(1, "m^3/(mol*s)")),
        reversible=True,
    )
    reaction.family = "R_Addition_MultipleBond"
    reaction.template = []
    reaction.is_forward = True
    return reaction, (radical.id, tail.id, head.id)


@pytest.mark.parametrize("units", [2, 3, 4])
def test_head_to_tail_maps_and_opposite_regioisomer_are_distinct(units):
    reaction, (radical, tail, head) = mapped_addition(units)
    product = reaction.products[0].molecule[0]
    assert product.is_isomorphic(Molecule(smiles="CC(c1ccccc1)" * units + "C[CH](c1ccccc1)"))
    assert not product.is_isomorphic(mapped_addition(units, opposite=True)[0].products[0].molecule[0])
    instance = compiler.EventSetCompiler(None, [], [], temperature_grid=[500, 600])
    proxy = compiler.SiteProxy("benzylic_end_radical+styrene", reaction.reactants,
                               participant_site_types=("benzylic_end_radical", "styrene"))
    record = instance._record(proxy, reaction, {})
    # Persistent map labels explicitly identify the radical attacking styrene CH2.
    atoms = [atom for species in reaction.reactants
             for atom in compiler._ordered_atoms(compiler._molecule(species))]
    labels = {atom.id: index for index, atom in enumerate(atoms)}
    formed = [op for op in record.bond_ops if op["action"] == "form"]
    assert any(set(op["atoms"]) == {labels[radical], labels[tail]} for op in formed)
    assert any(set(op["atoms"]) == {labels[head], labels[tail]}
               and float(op["order"]) == 1 for op in formed)
    assert any(op["action"] == "set_radical" and op["atom"] == labels[head]
               and op["value"] == 1 for op in record.bond_ops)


class SmallReferenceThermo:
    provenance = {"reference_thermo": "explicit fast-test Kc; no database"}

    def equilibrium_constants(self, reaction, temperatures):
        carbon_count = reaction.products[0].molecule[0].get_element_count()["C"]
        # Distinct crossings; the designated n=3 anchor is deliberately not lowest.
        midpoint = 580 if carbon_count == 24 else 520
        return [10 ** ((midpoint - temperature) / 50) for temperature in temperatures]


def small_inventory():
    proxies = compiler.ps_proxy_set()
    cache = {proxy.site_type: [] for proxy in proxies}
    for size, site, primary in [(2, "benzylic_end_radical+styrene", False),
                                (4, "benzylic_end_radical+styrene@5", False),
                                (2, "end_radical+styrene", True),
                                (5, "end_radical+styrene@5", True)]:
        cache[site] = [mapped_addition(size, primary=primary)[0]]
    # A wrong regioisomer must be discovered but never selected as a ceiling anchor.
    cache["benzylic_end_radical+styrene"].append(mapped_addition(2, opposite=True)[0])
    cache["end_radical@5"] = [mapped_addition(4, primary=True)[0]]
    return compiler.EventSetCompiler(
        None, proxies, ["R_Addition_MultipleBond"], reaction_cache=cache,
        temperature_grid=[500, 550, 600], rmgpy_sha="fast-test", rmg_database_sha="fast-test",
        reference_thermo_provider=SmallReferenceThermo(), ceiling_monomer_concentration_mol_m3=1,
    ).compile()


def test_ceiling_pairs_use_exact_benzylic_anchors_and_keep_primary_pairs():
    artifact = small_inventory()
    assert len(artifact["ps_ceiling_pairs"]) == 2
    assert len(artifact["ps_primary_end_ceiling_pairs"]) == 2
    assert artifact["ps_ceiling_anchor_event_id"] == artifact["ps_ceiling_pairs"][0]["propagation_event_id"]
    assert artifact["ps_ceiling_temperature_K"] == pytest.approx(580)
    assert artifact["ps_ceiling_pairs"][1]["temperature_K"] == pytest.approx(520)
    by_id = {record["event_id"]: record for record in artifact["records"]}
    for field in ("ps_ceiling_pairs", "ps_primary_end_ceiling_pairs"):
        for pair in artifact[field]:
            prop, dep = (by_id[pair[key]] for key in ("propagation_event_id", "depropagation_event_id"))
            assert prop["reverse_of"] == dep["event_id"]
            assert dep["reverse_of"] == prop["event_id"]
            assert prop["rate_units"] == "m^3/(mol*s)"
            assert dep["rate_units"] == "s^-1"
            assert [kf / kr for kf, kr in zip(prop["k_table"]["k"], dep["k_table"]["k"])] == pytest.approx(
                prop["thermo_provenance"]["equilibrium_constant_table"]["Kc"])
            expected_label = "benzylic_end_radical" if field == "ps_ceiling_pairs" else "end_radical"
            assert dep["participant_site_types"] == [expected_label]


@pytest.mark.parametrize("units", [1, 2, 3])
@pytest.mark.parametrize("primary", [False, True])
def test_seed_site_index_classifies_both_radicals_as_ends(units, primary):
    from fixtures.i049_probe.runtime import seed, snapshot
    from rmgpy.kmc.ssa import SiteIndex, radical_site_class

    smiles = (compiler.linear_ps_smiles(units, radical_end=True) if primary
              else compiler.benzylic_ps_end_smiles(units))
    label = "end_radical" if primary else "benzylic_end_radical"
    record = {"event_id": "seed", "arity": 1, "status": "enabled",
              "participant_site_types": [label],
              "reactant_graphs": [Molecule(smiles=smiles).to_adjacency_list(remove_h=False)]}
    state = seed(record)
    candidates = SiteIndex(state, [record]).candidates("seed")
    assert len(candidates) == 1
    assert radical_site_class(candidates[0][0], state) == "end"
    assert snapshot(state)["graph_radicals"] == state.total_radicals == 1


@pytest.mark.parametrize("units", [2, 3, 4])
def test_mapped_growth_unzip_preserves_uuids_ledger_and_alternating_backbone(units):
    from fixtures.i049_probe.runtime import seed, snapshot
    from rmgpy.kmc.ssa import SiteIndex, radical_site_class

    reaction, _ = mapped_addition(units)
    instance = compiler.EventSetCompiler(
        None, [], [], temperature_grid=[500, 550, 600],
        reference_thermo_provider=SmallReferenceThermo())
    proxy = compiler.SiteProxy("benzylic_end_radical+styrene", reaction.reactants,
                               participant_site_types=("benzylic_end_radical", "styrene"))
    records, _ = instance._linked_family_pair(proxy, reaction, {})
    prop, dep = (record.to_dict() for record in records)
    state = seed(prop)
    before, uuids = snapshot(state), set(state._uuid_owners())
    state.apply(prop, SiteIndex(state, [prop]).candidates(prop["event_id"])[0])
    after = snapshot(state)
    assert after["strand_lengths"] == [2 * (units + 1)]
    assert after["graph_radicals"] == after["ledger_radicals"] == 1
    assert after["sink_radicals"] == 0
    assert set(state._uuid_owners()) == uuids
    # Exact ordinary PS product verifies the alternating backbone as well as end chemistry.
    expected = Molecule(smiles="CC(c1ccccc1)" * units + "C[CH](c1ccccc1)")
    from fixtures.i049_probe.common import canonical_smiles
    assert [row["smiles"] for row in after["components"]] == [canonical_smiles(expected.to_smiles())]
    inverse_sites = SiteIndex(state, [dep]).candidates(dep["event_id"])
    assert inverse_sites
    assert radical_site_class(inverse_sites[0][0], state) == "end"
    state.apply(dep, inverse_sites[0])
    restored = snapshot(state)
    assert sorted(row["smiles"] for row in restored["components"]) == sorted(
        row["smiles"] for row in before["components"])
    assert restored["graph_radicals"] == restored["ledger_radicals"] == 1
    assert restored["strand_lengths"] == before["strand_lengths"]
    assert restored["sink_radicals"] == 0
    assert set(state._uuid_owners()) == uuids
    state.assert_uuid_uniqueness()


def test_independent_completeness_inputs_include_benzylic_contexts(monkeypatch):
    from completeness_oracle import full_molecule_inputs

    def forbidden(*args, **kwargs):
        raise AssertionError("oracle must not borrow compiler fragments")

    monkeypatch.setattr(compiler, "benzylic_ps_end_smiles", forbidden)
    monkeypatch.setattr(compiler, "linear_ps_smiles", forbidden)
    inputs = dict(full_molecule_inputs(1))
    assert len(inputs) == 17
    for name, labels in BENZYLIC_DECLARATIONS.items():
        for species, label in zip(inputs[name], labels):
            if label == "benzylic_end_radical":
                units = 4 if name.endswith("+styrene") else 5
                assert species.molecule[0].is_isomorphic(
                    Molecule(smiles="CC(c1ccccc1)" * (units - 1) + "C[CH](c1ccccc1)"))


def test_explicit_artifact_real_oracle_mode_cannot_launch_a_compile(monkeypatch, tmp_path):
    import compilerRealTest as real

    artifact = {"families": [], "excluded_families": {}}
    oracle = {"keys": {("benzylic_end_radical", "H_Abstraction", "mock")}}
    monkeypatch.setattr(real, "supplied_artifact", lambda: (tmp_path / "artifact.json", artifact))
    monkeypatch.setattr(real, "_load_independent_oracle", lambda kinetics: oracle)
    monkeypatch.setenv("RMG_KMC_REAL_ORACLE", "1")
    def forbidden(*args, **kwargs):
        raise AssertionError("explicit artifacts must never trigger a compile")
    monkeypatch.setattr(compiler.EventSetCompiler, "compile", forbidden)
    result = real.compilation.__wrapped__(SimpleNamespace(kinetics=None))
    assert result[1] is artifact and result[7] is oracle


@pytest.mark.phase2b
def test_compiled_benzylic_declarations_and_anchors_require_phase2b(ps_artifact):
    """Deliberate failure with the old artifact, never skip new inventory checks."""
    assert_compiled_benzylic_anchors(ps_artifact)


def assert_compiled_benzylic_anchors(artifact):
    from plpsecLibraryTest import assert_compiled_plpsec_library
    from fixtures.i049_probe.common import canonical_smiles, describe

    declared = {proxy["site_type"] for proxy in artifact["inputs"]["proxies"]}
    required = {name + suffix for name in BENZYLIC_DECLARATIONS for suffix in ("", "@5")}
    assert required <= declared, "phase-2b compile required: benzylic declarations absent"
    assert artifact["short_molecule_catalogue"]["size"] == 37
    assert len(artifact["ps_primary_end_ceiling_pairs"]) == 2
    assert {pair["product_repeat_units"] for pair in artifact["ps_ceiling_pairs"]} == {3, 5}
    by_id = {record["event_id"]: record for record in artifact["records"]}
    anchor = artifact["ps_ceiling_anchor_event_id"]
    assert anchor in {pair["propagation_event_id"] for pair in artifact["ps_ceiling_pairs"]}
    for pair in artifact["ps_ceiling_pairs"]:
        prop, dep = (by_id[pair[key]] for key in ("propagation_event_id", "depropagation_event_id"))
        units = pair["product_repeat_units"]
        assert sorted(describe(graph)["smiles"] for graph in prop["reactant_graphs"]) == sorted([
            canonical_smiles("CC(c1ccccc1)" * (units - 2) + "C[CH](c1ccccc1)"),
            canonical_smiles("C=Cc1ccccc1"),
        ])
        assert [describe(graph)["smiles"] for graph in prop["product_graphs"]] == [
            canonical_smiles("CC(c1ccccc1)" * (units - 1) + "C[CH](c1ccccc1)")]
        assert prop["reverse_of"] == dep["event_id"] and dep["reverse_of"] == prop["event_id"]
        assert dep["participant_site_types"] == ["benzylic_end_radical"]
        assert prop["rate_units"] == "m^3/(mol*s)" and dep["rate_units"] == "s^-1"
        assert [kf / kr for kf, kr in zip(prop["k_table"]["k"], dep["k_table"]["k"])] == pytest.approx(
            prop["thermo_provenance"]["equilibrium_constant_table"]["Kc"])
        if prop["event_id"] == anchor:
            assert units == 3
            assert artifact["ps_ceiling_temperature_K"] == pair["temperature_K"]

    assert assert_compiled_plpsec_library(artifact) == 2


def test_compiled_anchor_verifier_accepts_a_small_inventory():
    assert_compiled_benzylic_anchors(small_inventory())
