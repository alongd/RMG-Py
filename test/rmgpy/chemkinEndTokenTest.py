#!/usr/bin/env python3
"""RMG's Chemkin reader ends the REACTIONS and THERM blocks only on a line
whose first token is Chemkin's END keyword, not on any line holding 'end'."""

from rmgpy.chemkin import load_chemkin_file
from rmgpy.molecule import Molecule

THERM_ENTRY = [
    "{0:<24}C   1H   4          G   200.000  6000.000 1000.00      1",
    " 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2",
    "-7.45375000E+02 3.35532000E+00 2.50000000E+00 0.00000000E+00 0.00000000E+00    3",
    " 0.00000000E+00 0.00000000E+00-7.45375000E+02 3.35532000E+00                   4",
]


def test_species_named_with_end_keep_their_reactions_and_thermo(tmp_path):
    names = {"blend(3)": "CC", "end_cap(4)": "C", "H(1)": "[H]"}
    therm = ["THERM ALL", "   300.000  1000.000  5000.000"]
    for name in names:
        therm += [line.format(name) for line in THERM_ENTRY]
    text = "\n".join(["ELEMENTS H C END", "SPECIES", *names, "END", *therm, "END",
                      "REACTIONS    KCAL/MOLE   MOLES",
                      "blend(3)=>end_cap(4)+end_cap(4)    1.0 0.0 0.0",
                      "end_cap(4)+H(1)=>blend(3)    2.0 0.0 0.0",
                      "H(1)+H(1)=>blend(3)    3.0 0.0 0.0",
                      "end", ""])
    (tmp_path / "chem.inp").write_text(text)
    (tmp_path / "dict.txt").write_text("\n".join(f"{n}\n{Molecule(smiles=s).to_adjacency_list()}"
                                                 for n, s in names.items()))
    species, reactions = load_chemkin_file(str(tmp_path / "chem.inp"), str(tmp_path / "dict.txt"),
                                           check_duplicates=False, generate_resonance_structures=False)
    assert [str(r) for r in reactions] == ["blend(3) => end_cap(4) + end_cap(4)", "end_cap(4) + H(1) => blend(3)",
                                           "H(1) + H(1) => blend(3)"]
    assert [r.kinetics.A.value_si for r in reactions][0] == 1.0
    assert all(s.thermo is not None for s in species), [s.label for s in species if s.thermo is None]
