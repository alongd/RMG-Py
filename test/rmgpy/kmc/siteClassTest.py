import json
from pathlib import Path

import pytest

from rmgpy.kmc.site_class import extract_site_class
from rmgpy.kmc.site_class import SiteClassError
from rmgpy.molecule.element import get_element
from rmgpy.molecule.molecule import Molecule


FIXTURE_DIR = Path(__file__).with_name("fixtures") / "site_class_v1"


def policy(**changes):
    value = json.loads((FIXTURE_DIR / "policy.json").read_text())
    value.update(changes)
    return value


def heavy_atoms(molecule):
    return [atom for atom in molecule.atoms if atom.symbol != "H"]


def middle_atom(molecule):
    heavy = heavy_atoms(molecule)
    return heavy[len(heavy) // 2]


def test_instance_metadata_and_insertion_order_do_not_change_key_or_mapping():
    first = Molecule(smiles="CCCCC")
    second = first.copy(deep=True)
    first_root = middle_atom(first)
    second_root = middle_atom(second)

    for index, atom in enumerate(first.atoms):
        atom.id = index + 10
        atom.label = f"display-{index}"
        atom.props["uuid"] = f"first-{index}"
    for index, atom in enumerate(second.atoms):
        atom.id = index + 1000
        atom.label = f"renamed-{index}"
        atom.props["uuid"] = f"second-{index}"
    second.vertices.reverse()

    first_key, first_mapping, first_pattern = extract_site_class(
        first, (("center", first_root),), policy()
    )
    second_key, second_mapping, second_pattern = extract_site_class(
        second, (("center", second_root),), policy()
    )

    assert first_key == second_key
    first_position = first_pattern["ordered_root"][0]["position"]
    second_position = second_pattern["ordered_root"][0]["position"]
    assert first_mapping[first_position] is first_root
    assert second_mapping[second_position] is second_root
    assert set(first_mapping.values()).issubset(set(first.atoms))
    assert set(second_mapping.values()).issubset(set(second.atoms))


def test_root_roles_preserve_chemistry_and_symmetric_occurrences_keep_bindings():
    asymmetric = Molecule(smiles="CCCC")
    carbons = heavy_atoms(asymmetric)
    forward = extract_site_class(
        asymmetric, (("left", carbons[0]), ("right", carbons[1])), policy(radius=1)
    )[0]
    swapped = extract_site_class(
        asymmetric, (("left", carbons[1]), ("right", carbons[0])), policy(radius=1)
    )[0]
    assert forward != swapped
    renamed = extract_site_class(
        asymmetric, (("donor", carbons[0]), ("acceptor", carbons[1])), policy(radius=1)
    )[0]
    assert forward != renamed

    symmetric = Molecule(smiles="CCC")
    left, _, right = heavy_atoms(symmetric)
    left_key, left_mapping, left_pattern = extract_site_class(
        symmetric, (("end", left),), policy(radius=2)
    )
    right_key, right_mapping, right_pattern = extract_site_class(
        symmetric, (("end", right),), policy(radius=2)
    )
    assert left_key == right_key
    assert left_mapping[left_pattern["ordered_root"][0]["position"]] is left
    assert right_mapping[right_pattern["ordered_root"][0]["position"]] is right
    assert left_mapping is not right_mapping


def test_radius_is_heavy_atom_distance_and_pattern_retains_hydrogens_and_ports():
    molecule = Molecule(smiles="CCCCCCC")
    root = heavy_atoms(molecule)[3]
    _, mapping, pattern = extract_site_class(
        molecule, (("site", root),), policy(radius=1)
    )

    assert len(mapping) == 3
    assert len(pattern["atoms"]) == 3
    assert len(pattern["boundary_ports"]) == 2
    assert all(
        port["excluded"] == {"element": "C", "bond_orders": [1.0, 1.0]}
        for port in pattern["boundary_ports"]
    )
    assert sorted(len(atom["explicit_hydrogens"]) for atom in pattern["atoms"]) == [
        2,
        2,
        2,
    ]


def test_annotated_caps_are_ports_and_caps_or_reacting_atoms_cannot_be_roots():
    witness = Molecule(smiles="CCC")
    witness_carbons = heavy_atoms(witness)
    signature = {"element": "C", "bond_orders": [1.0, 1.0]}
    caps = []
    for boundary_carbon in (witness_carbons[0], witness_carbons[-1]):
        cap = next(atom for atom in boundary_carbon.bonds if atom.symbol == "H")
        cap.props["artificial_cap"] = True
        cap.props["boundary_port"] = signature
        caps.append(cap)

    interior = Molecule(smiles="CCCCCCC")
    witness_key = extract_site_class(
        witness, (("site", witness_carbons[1]),), policy(radius=1)
    )[0]
    interior_key = extract_site_class(
        interior, (("site", heavy_atoms(interior)[3]),), policy(radius=1)
    )[0]
    assert witness_key == interior_key

    with pytest.raises(SiteClassError, match="invalid_root"):
        extract_site_class(witness, (("site", caps[0]),), policy(radius=1))
    heavy_cap_witness = Molecule(smiles="CCCCC")
    heavy_cap = heavy_atoms(heavy_cap_witness)[0]
    heavy_cap.props["artificial_cap"] = True
    heavy_cap.props["boundary_port"] = signature
    with pytest.raises(SiteClassError, match="invalid_root"):
        extract_site_class(heavy_cap_witness, (("site", heavy_cap),), policy(radius=1))
    witness_carbons[1].props["reacting"] = True
    with pytest.raises(SiteClassError, match="invalid_root"):
        extract_site_class(witness, (("site", witness_carbons[1]),), policy(radius=1))


def test_every_supported_chemical_field_and_version_changes_the_key():
    def keyed(change=None, policy_change=None):
        molecule = Molecule(smiles="CCC")
        root = heavy_atoms(molecule)[1]
        if change is not None:
            change(molecule, root)
        return extract_site_class(
            molecule, (("site", root),), policy(**(policy_change or {}))
        )[0]

    def isotope(_, atom):
        atom.element = get_element("C", 13)

    def element(_, atom):
        atom.element = get_element("N")

    def charge(_, atom):
        atom.charge = 1

    def radical(_, atom):
        atom.radical_electrons = 1

    def lone_pair(_, atom):
        atom.lone_pairs = 1

    def hydrogen(molecule, atom):
        hydrogen_atom = next(
            neighbor for neighbor in atom.bonds if neighbor.symbol == "H"
        )
        molecule.remove_atom(hydrogen_atom)

    def bond_order(_, atom):
        next(
            bond for neighbor, bond in atom.bonds.items() if neighbor.symbol == "C"
        ).order = 2.0

    def atom_stereo(molecule, atom):
        molecule.props["atom_stereo"] = {atom: "R"}

    def bond_stereo(molecule, atom):
        bond = next(
            bond for neighbor, bond in atom.bonds.items() if neighbor.symbol == "C"
        )
        molecule.props["bond_stereo"] = {bond: "E"}

    def multiplicity(molecule, _):
        molecule.multiplicity = 3

    baseline = keyed()
    for change in (
        element,
        isotope,
        charge,
        radical,
        lone_pair,
        hydrogen,
        bond_order,
        atom_stereo,
        bond_stereo,
        multiplicity,
    ):
        assert keyed(change) != baseline, change.__name__
    assert keyed(policy_change={"policy_id": "another-policy"}) != baseline
    assert (
        keyed(policy_change={"canonicalizer_version": "site-class-canonicalizer-v2"})
        != baseline
    )
    assert keyed(policy_change={"orientation": "participant-1"}) != baseline


def test_invalid_inputs_and_every_finite_policy_bound_fail_without_mutation():
    molecule = Molecule(smiles="CCC")
    root = heavy_atoms(molecule)[1]
    before = molecule.to_adjacency_list()

    cases = [
        ((), policy(), "missing_root"),
        ((("site", root),), policy(radius=5), "domain_exceeded"),
        ((("site", root),), policy(max_degree=1), "domain_exceeded"),
        ((("site", root),), policy(max_explicit_h=1), "domain_exceeded"),
        ((("site", root),), policy(max_heavy_atoms=1), "domain_exceeded"),
        ((("site", root),), policy(max_root_size=0), "domain_exceeded"),
        ((("site", root),), policy(allowed_elements=["H"]), "domain_exceeded"),
        ((("site", root),), policy(allowed_isotopes=[13]), "domain_exceeded"),
        ((("site", root),), policy(allowed_charges=[1]), "domain_exceeded"),
        ((("site", root),), policy(allowed_radicals=[1]), "domain_exceeded"),
        ((("site", root),), policy(allowed_lone_pairs=[1]), "domain_exceeded"),
        ((("site", root),), policy(allowed_bond_orders=[2.0]), "domain_exceeded"),
        ((("site", root),), policy(allowed_multiplicities=[2]), "domain_exceeded"),
        ((("site", root),), policy(instance_uuid="run-specific"), "invalid_policy"),
        (
            (("site", root),),
            policy(max_valence_by_element={}),
            "invalid_policy",
        ),
    ]
    for roots, selected_policy, code in cases:
        with pytest.raises(SiteClassError) as error:
            extract_site_class(molecule, roots, selected_policy)
        assert error.value.code == code
        assert molecule.to_adjacency_list() == before

    missing_spin = molecule.copy(deep=True)
    missing_spin.multiplicity = -187
    with pytest.raises(SiteClassError, match="missing_multiplicity"):
        extract_site_class(
            missing_spin, (("site", heavy_atoms(missing_spin)[1]),), policy()
        )

    unsupported = molecule.copy(deep=True)
    unsupported.props["unknown_chemistry"] = True
    with pytest.raises(SiteClassError, match="unsupported_field"):
        extract_site_class(
            unsupported, (("site", heavy_atoms(unsupported)[1]),), policy()
        )

    bad_stereo = molecule.copy(deep=True)
    bad_root = heavy_atoms(bad_stereo)[1]
    bad_stereo.props["atom_stereo"] = {bad_root: "tetrahedral-unknown"}
    with pytest.raises(SiteClassError, match="unsupported_stereo"):
        extract_site_class(bad_stereo, (("site", bad_root),), policy())

    bad_cap = molecule.copy(deep=True)
    cap = heavy_atoms(bad_cap)[0]
    cap.props["artificial_cap"] = True
    with pytest.raises(SiteClassError, match="malformed_port"):
        extract_site_class(
            bad_cap, (("site", heavy_atoms(bad_cap)[1]),), policy(radius=1)
        )

    orphan_port = molecule.copy(deep=True)
    orphan_port.atoms[0].props["boundary_port"] = {
        "element": "C",
        "bond_orders": [1.0],
    }
    with pytest.raises(SiteClassError, match="malformed_port"):
        extract_site_class(
            orphan_port, (("site", heavy_atoms(orphan_port)[1]),), policy()
        )

    non_hydrogen_cap = molecule.copy(deep=True)
    non_hydrogen_cap.atoms[0].props["artificial_cap"] = True
    non_hydrogen_cap.atoms[0].props["boundary_port"] = {
        "element": "C",
        "bond_orders": [1.0],
    }
    with pytest.raises(SiteClassError, match="malformed_port"):
        extract_site_class(
            non_hydrogen_cap,
            (("site", heavy_atoms(non_hydrogen_cap)[1]),),
            policy(),
        )

    nonnumeric_port = molecule.copy(deep=True)
    cap = next(atom for atom in nonnumeric_port.atoms if atom.symbol == "H")
    cap.props["artificial_cap"] = True
    cap.props["boundary_port"] = {"element": "C", "bond_orders": ["single"]}
    with pytest.raises(SiteClassError, match="malformed_port"):
        extract_site_class(
            nonnumeric_port,
            (("site", heavy_atoms(nonnumeric_port)[1]),),
            policy(),
        )


def test_linear_alkane_interior_has_one_class_across_chain_lengths():
    short = Molecule(smiles="CCCCCCCCC")
    long = Molecule(smiles="CCCCCCCCCCC")
    short_key = extract_site_class(
        short, (("interior", heavy_atoms(short)[4]),), policy(radius=2)
    )[0]
    long_key = extract_site_class(
        long, (("interior", heavy_atoms(long)[5]),), policy(radius=2)
    )[0]
    assert short_key == long_key


def methyl_branched_roots(molecule):
    carbons = heavy_atoms(molecule)
    heavy_degree = {
        atom: sum(neighbor.symbol != "H" for neighbor in atom.bonds) for atom in carbons
    }
    tertiary = [atom for atom in carbons if heavy_degree[atom] == 3]
    backbone_ch = tertiary[len(tertiary) // 2]
    pendant = next(
        neighbor for neighbor in backbone_ch.bonds if heavy_degree.get(neighbor) == 1
    )
    backbone_ch2 = next(
        neighbor
        for neighbor in backbone_ch.bonds
        if heavy_degree.get(neighbor) == 2
        and any(
            heavy_degree.get(other) == 3
            for other in neighbor.bonds
            if other is not backbone_ch
        )
    )
    return backbone_ch2, backbone_ch, pendant


def test_methyl_branched_alkane_site_classes_are_length_stable_and_distinct():
    short = Molecule(smiles="CC(C)CC(C)CC(C)CC(C)CC(C)C")
    long = Molecule(smiles="CC(C)CC(C)CC(C)CC(C)CC(C)CC(C)CC(C)C")
    short_roots = methyl_branched_roots(short)
    long_roots = methyl_branched_roots(long)

    short_keys = [
        extract_site_class(short, (("site", root),), policy(radius=2))[0]
        for root in short_roots
    ]
    long_keys = [
        extract_site_class(long, (("site", root),), policy(radius=2))[0]
        for root in long_roots
    ]
    assert short_keys == long_keys
    assert len(set(short_keys)) == 3

    for root_index in (0, 1):
        short_radical = short.copy(deep=True)
        long_radical = long.copy(deep=True)
        short_root = methyl_branched_roots(short_radical)[root_index]
        long_root = methyl_branched_roots(long_radical)[root_index]
        short_root.radical_electrons = 1
        long_root.radical_electrons = 1
        short_radical.multiplicity = 2
        long_radical.multiplicity = 2
        assert (
            extract_site_class(
                short_radical, (("site", short_root),), policy(radius=2)
            )[0]
            == extract_site_class(
                long_radical, (("site", long_root),), policy(radius=2)
            )[0]
        )


def aromatic_substituted_chain_roots(molecule):
    def aromatic(atom):
        return atom.props.get("inRing") is True

    backbone = [
        atom
        for atom in heavy_atoms(molecule)
        if atom.symbol == "C" and not aromatic(atom)
    ]
    backbone_set = set(backbone)
    benzylic = [
        atom for atom in backbone if any(aromatic(neighbor) for neighbor in atom.bonds)
    ]
    methylene = [
        atom
        for atom in backbone
        if atom not in benzylic
        and sum(neighbor in backbone_set for neighbor in atom.bonds) == 2
    ]
    return methylene[len(methylene) // 2], benzylic[len(benzylic) // 2]


def test_aromatic_substituted_oligomer_backbone_classes_are_length_stable():
    short = Molecule(
        smiles="CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)C"
    )
    long = Molecule(
        smiles="CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)C"
    )
    short_roots = aromatic_substituted_chain_roots(short)
    long_roots = aromatic_substituted_chain_roots(long)
    for short_root, long_root in zip(short_roots, long_roots):
        assert (
            extract_site_class(short, (("site", short_root),), policy(radius=2))[0]
            == extract_site_class(long, (("site", long_root),), policy(radius=2))[0]
        )


def test_bond_roots_and_boundary_signatures_are_canonical_pattern_data():
    molecule = Molecule(smiles="CCCCCCC")
    carbons = heavy_atoms(molecule)
    root_bond = carbons[3].bonds[carbons[4]]
    key, mapping, pattern = extract_site_class(
        molecule,
        (("atom-role", carbons[3]), ("bond-role", root_bond)),
        policy(radius=1),
    )
    assert len(key) == 64
    assert len(mapping) == 4
    assert pattern["ordered_root"][1]["kind"] == "bond"
    json.dumps(pattern, sort_keys=True)

    changed = molecule.copy(deep=True)
    changed_carbons = heavy_atoms(changed)
    changed_carbons[1].element = get_element("N")
    original_port_key = extract_site_class(
        molecule, (("site", carbons[3]),), policy(radius=1)
    )[0]
    changed_port_key = extract_site_class(
        changed, (("site", changed_carbons[3]),), policy(radius=1)
    )[0]
    assert original_port_key != changed_port_key
