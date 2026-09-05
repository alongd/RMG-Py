#!/usr/bin/env python
"""i213 probe: does one atom type for a bonded NEUTRAL argon make Ar2+ buildable?

Run from the repository root of this worktree, so that ./rmgrc is the file picked up:

    python docs/i213-bonded-argon-atomtype/probe.py

Prints, in order: a provenance block (which rmgpy is loaded, whether it is compiled, which
argon labels the LOADED binary actually holds, which database directory resolved), then each
measurement with its outcome.  Every failure is printed with its full traceback rather than
aborting the run, so one script produces the whole before/after record.
"""

import os
import sys
import traceback

SEP = "=" * 78


def section(title):
    print(f"\n{SEP}\n{title}\n{SEP}")


def report(name, fn):
    """Run fn; print its value, or the full traceback of whatever it raised."""
    try:
        value = fn()
    except Exception:
        print(f"[RAISED] {name}")
        traceback.print_exc(file=sys.stdout)
        return None
    print(f"[OK]     {name}: {value}")
    return value


# --------------------------------------------------------------------------------------
section("PROVENANCE")
# --------------------------------------------------------------------------------------
print(f"cwd                 = {os.getcwd()}")
print(f"python              = {sys.executable}")

import rmgpy
import rmgpy.molecule.atomtype as atomtype_module

print(f"rmgpy.__file__      = {rmgpy.__file__}")
print(f"atomtype.__file__   = {atomtype_module.__file__}")
print(f"atomtype compiled   = {atomtype_module.__file__.endswith('.so')}")

# A .so can predate the source edit.  Ask the loaded binary which labels it holds, rather
# than trusting the .py sitting next to it.
from rmgpy.molecule.atomtype import ATOMTYPES, AtomTypeError, get_atomtype

argon_labels = sorted(label for label in ATOMTYPES if label.startswith("Ar"))
print(f"argon labels LOADED = {argon_labels}")
for label in argon_labels:
    at = ATOMTYPES[label]
    print(
        f"    {label:6s} single={at.single} lone_pairs={at.lone_pairs} charge={at.charge} "
        f"specific={[t.label for t in at.specific]}"
    )

from rmgpy import settings

print(f"database.directory  = {settings['database.directory']}")
print(f"  exists            = {os.path.isdir(settings['database.directory'])}")
print(f"  realpath          = {os.path.realpath(settings['database.directory'])}")

from rmgpy.molecule import Molecule

# --------------------------------------------------------------------------------------
section("1. Ar2+ FROM A LOCALISED-CHARGE ADJACENCY LIST")
# --------------------------------------------------------------------------------------
# Ar2+ ground state is 2-Sigma-u+, a doublet cation.  Localising the hole on atom 1:
#   atom 1  Ar+ : 7 valence e- = 1 (bond) + 6 (p3)           -> u0 p3 c+1
#   atom 2  Ar  : 8 valence e- = 1 (bond) + 6 (p3) + 1 (u1)  -> u1 p3 c0
AR2_PLUS_ADJLIST = """1 Ar u0 p3 c+1 {2,S}
2 Ar u1 p3 c0 {1,S}
"""
print("adjacency list under test:")
print(AR2_PLUS_ADJLIST)

mol = report(
    "Molecule().from_adjacency_list(Ar2+)",
    lambda: Molecule().from_adjacency_list(AR2_PLUS_ADJLIST),
)

# The narrower control: the bonded NEUTRAL argon centre on its own, which is the atom the
# ticket says has no type.  Its cation partner is the positive control -- it must type today.
section("1b. SINGLE-ATOM CONTROLS")


# update_atomtypes stops at the FIRST atom it cannot resolve, so the two centres have to be
# asked about individually.  Build without the exception, then type each atom by hand.
def _unraised_pair():
    return Molecule().from_adjacency_list(AR2_PLUS_ADJLIST, raise_atomtype_exception=False)


def bonded_cation_argon():
    m = _unraised_pair()
    return get_atomtype(m.atoms[0], m.atoms[0].bonds).label


def bonded_neutral_argon():
    m = _unraised_pair()
    return get_atomtype(m.atoms[1], m.atoms[1].bonds).label


report("get_atomtype(bare Ar, u0 p4 c0)", lambda: Molecule().from_adjacency_list("1 Ar u0 p4 c0").atoms[0].atomtype.label)
report("get_atomtype(bare Ar+, u1 p3 c+1)", lambda: Molecule().from_adjacency_list("1 Ar u1 p3 c+1").atoms[0].atomtype.label)
report("get_atomtype(BONDED Ar+ centre, positive control)", bonded_cation_argon)
report("get_atomtype(BONDED neutral Ar centre, the gap)", bonded_neutral_argon)

if mol is not None:
    # ----------------------------------------------------------------------------------
    section("2. UPDATE, ROUND-TRIP, FORMULA, BALANCE")
    # ----------------------------------------------------------------------------------
    report("atom types", lambda: [a.atomtype.label for a in mol.atoms])
    report("mol.update()", lambda: (mol.update(), [a.atomtype.label for a in mol.atoms])[1])
    report("get_formula()", lambda: mol.get_formula())
    report("get_net_charge()", lambda: mol.get_net_charge())
    report("get_radical_count()", lambda: mol.get_radical_count())
    report("multiplicity", lambda: mol.multiplicity)
    report("is_linear()", lambda: mol.is_linear())
    report("get_molecular_weight()", lambda: mol.get_molecular_weight())

    def roundtrip():
        text = mol.to_adjacency_list()
        again = Molecule().from_adjacency_list(text)
        return f"\n---\n{text}--- isomorphic back: {mol.is_isomorphic(again)}"

    report("to_adjacency_list -> from_adjacency_list", roundtrip)

    def balance():
        from rmgpy.reaction import Reaction
        from rmgpy.species import Species

        ar = Species(label="Ar", molecule=[Molecule().from_adjacency_list("1 Ar u0 p4 c0")])
        ar_cation = Species(label="Ar+", molecule=[Molecule().from_adjacency_list("1 Ar u1 p3 c+1")])
        ar2_cation = Species(label="Ar2+", molecule=[mol.copy(deep=True)])
        electron = Species(label="e", molecule=[Molecule().from_adjacency_list("1 e u0 p0 c-1")])
        good = Reaction(reactants=[ar_cation, ar], products=[ar2_cation])
        bad = Reaction(reactants=[ar_cation], products=[ar2_cation, electron])  # negative control
        return (
            f"Ar+ + Ar <=> Ar2+ : {good.is_balanced()} "
            f"| negative control Ar+ <=> Ar2+ + e- : {bad.is_balanced()}"
        )

    report("is_balanced in a trial reaction", balance)

    # ----------------------------------------------------------------------------------
    section("3. SMILES / InChI ROUND-TRIP")
    # ----------------------------------------------------------------------------------
    def smiles_roundtrip():
        smi = mol.to_smiles()
        back = Molecule().from_smiles(smi)
        shape = [(a.symbol, a.charge, a.radical_electrons, a.lone_pairs) for a in back.atoms]
        return (
            f"to_smiles() = {smi!r} -> from_smiles gives {shape} "
            f"net charge {back.get_net_charge()} | isomorphic back: {mol.is_isomorphic(back)}"
        )

    report("SMILES round-trip (Ar2+)", smiles_roundtrip)

    def smiles_monatomic_control():
        """The separately-gated known defect: monatomic [Ar+] reads back as the dication."""
        cation = Molecule().from_adjacency_list("1 Ar u1 p3 c+1")
        smi = cation.to_smiles()
        back = Molecule().from_smiles(smi)
        return (
            f"to_smiles() = {smi!r} -> back net charge {back.get_net_charge()} (expected +1) "
            f"| isomorphic back: {cation.is_isomorphic(back)}"
        )

    report("SMILES round-trip control (monatomic Ar+)", smiles_monatomic_control)
    report("to_inchi()", lambda: mol.to_inchi())

    # ----------------------------------------------------------------------------------
    section("4. THERMOCHEMISTRY: DOES GROUP ADDITIVITY RETURN A SILENT NUMBER?")
    # ----------------------------------------------------------------------------------
    from rmgpy.data.rmg import RMGDatabase
    from rmgpy.species import Species

    # NOTE: RMGDatabase.load() only reaches load_thermo() when surface=True (rmgpy/data/rmg.py
    # line 112), so the thermo database is loaded here directly rather than through load().
    db = RMGDatabase()
    db.load_thermo(os.path.join(settings["database.directory"], "thermo"),
                   thermo_libraries=[], depository=False, surface=False)

    def _describe(thermo):
        return (
            "\n    !!! A NUMBER WAS RETURNED FOR Ar2+ !!!"
            f"\n    H298    = {thermo.get_enthalpy(298.0) / 4184.0:.3f} kcal/mol"
            f"\n    S298    = {thermo.get_entropy(298.0) / 4.184:.3f} cal/mol/K"
            f"\n    Cp(300) = {thermo.get_heat_capacity(300.0) / 4.184:.3f} cal/mol/K"
            f"\n    comment = {thermo.comment!r}"
        )

    def _ar2_species():
        spc = Species(label="Ar2+", molecule=[mol.copy(deep=True)])
        spc.generate_resonance_structures()
        return spc

    # the group-additivity estimator itself
    report("thermo.get_thermo_data_from_groups(Ar2+)",
           lambda: _describe(db.thermo.get_thermo_data_from_groups(_ar2_species())))
    # the full lookup RMG actually calls for a species (library -> QM -> GA)
    report("thermo.get_thermo_data(Ar2+)",
           lambda: _describe(db.thermo.get_thermo_data(_ar2_species())))
    # the HBI radical-saturation step is reached first, so also call the group summation
    # DIRECTLY on Ar2+, bypassing HBI, to see whether the group tree itself would yield a number
    report("thermo.compute_group_additivity_thermo(Ar2+) [HBI bypassed]",
           lambda: _describe(db.thermo.compute_group_additivity_thermo(mol.copy(deep=True))))
    # positive control: the same two calls on a species GA is known to handle
    control = Species(label="ethane", molecule=[Molecule().from_smiles("CC")])
    control.generate_resonance_structures()
    report("CONTROL thermo.get_thermo_data_from_groups(ethane)",
           lambda: _describe(db.thermo.get_thermo_data_from_groups(control)))
    # negative control: monatomic Ar, an element GA has no group for either
    argon = Species(label="Ar", molecule=[Molecule().from_adjacency_list("1 Ar u0 p4 c0")])
    report("CONTROL thermo.get_thermo_data_from_groups(monatomic Ar)",
           lambda: _describe(db.thermo.get_thermo_data_from_groups(argon)))

section("DONE")
