"""Run: PYTHONPATH=$PWD python test/rmgpy/kmc/fixtures/i049_probe/direct.py

Verify the corrected catalogue/discovery coverage using independent graph checks.
"""

import json
from itertools import combinations

from common import benzylic_tail, describe
from rmgpy.kmc.compiler import linear_ps_smiles, ps_proxy_set, short_ps_molecule_catalogue
from rmgpy.molecule.molecule import Molecule


def main():
    catalogue = short_ps_molecule_catalogue(1)
    rows = []
    for units in (1, 2, 3):
        base = Molecule(smiles=linear_ps_smiles(units))
        backbone = [a for a in base.atoms if a.element.number == 6
                    and not base.is_atom_in_cycle(a)]
        actual = [m for m in catalogue["molecules"] if m["repeat_units"] == units]
        stored = [Molecule().from_adjacency_list(m["adjacency_list"]) for m in actual]
        exhaustive = []
        missing = []
        for count in range(3):
            for sites in combinations(range(2 * units), count):
                mol = base.copy(deep=True)
                atoms = [a for a in mol.atoms if a.element.number == 6
                         and not mol.is_atom_in_cycle(a)]
                for site in sites:
                    atom = atoms[site]
                    hydrogen = next(a for a in atom.edges if a.element.number == 1)
                    mol.remove_atom(hydrogen)
                    atom.radical_electrons += 1
                mol.update(sort_atoms=False)
                assert not any(mol.is_isomorphic(other) for other in exhaustive)
                exhaustive.append(mol)
                if not any(mol.is_isomorphic(other) for other in stored):
                    missing.append({"sites": sites,
                                    **describe(mol.to_adjacency_list(remove_h=False))})
        reflection_preserves = []
        for site, atom in enumerate(backbone):
            reflected = backbone[-1 - site]
            attached_ring = lambda a: any(base.is_atom_in_cycle(b) for b in a.edges
                                         if b.element.number == 6)
            reflection_preserves.append(attached_ring(atom) == attached_ring(reflected))
        assert not any(reflection_preserves)
        tail = Molecule(smiles=benzylic_tail(units))
        assert any(tail.is_isomorphic(m) for m in stored)
        rows.append({"units": units, "stored_count": len(stored),
                     "correct_count": len(exhaustive), "tail_present": tail.to_smiles(),
                     "missing": missing,
                     "stored": [{"sites": m["radical_sites"],
                                 **describe(m["adjacency_list"])} for m in actual]})
    assert [r["stored_count"] for r in rows] == [4, 11, 22]
    assert all(not r["missing"] for r in rows)
    assert [r["correct_count"] for r in rows] == [4, 11, 22]
    proxies = []
    for proxy in ps_proxy_set():
        participants = [describe(s.molecule[0].to_adjacency_list(remove_h=False))
                        for s in proxy.reactants]
        if "benzylic_end_radical" in proxy.site_type:
            assert any(r["terminal_benzylic"] for p in participants for r in p["radicals"])
        proxies.append({"site_type": proxy.site_type,
                        "participant_site_types": proxy.participant_site_types,
                        "smiles": [p["smiles"] for p in participants]})
    assert len(proxies) == 34
    print(json.dumps({"catalogue": rows, "proxies": proxies}, indent=2))
    print("PASS: n=1–3 catalogue 37/37; 34 declarations include both end orientations")


if __name__ == "__main__":
    main()
