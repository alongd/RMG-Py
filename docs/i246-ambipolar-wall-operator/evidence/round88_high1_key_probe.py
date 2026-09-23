"""Round 88 HIGH 1 probe: InChI layer order, and the remove-/q-/p-layers key.

The premise (verified here, not assumed): InChI orders its layers
formula / c / h / q / p / stereo / i, so the isotope layer /i sits AFTER the charge
layer /q. Truncating at the first /q therefore drops /i from a CHARGED species while a
neutral keeps it -- the key becomes a different rule depending on charge, and a 13C
cation keys identically to an ordinary-carbon neutral. Removing only /q and /p (a
first-char filter on the /-split tokens) keeps /i on both charge states."""
import os, sys
sys.path.insert(0, os.environ['PYTHONPATH'])
from rmgpy.molecule import Molecule


def inchi(s):
    try:
        return Molecule().from_smiles(s).to_inchi()
    except Exception as e:
        return f"<ERR {e}>"


def old_key(ic):  # truncate at first /q or /p (the round-79/83 key)
    cut = len(ic)
    for sep in ('/q', '/p'):
        idx = ic.find(sep)
        if idx != -1 and idx < cut:
            cut = idx
    return ic[:cut]


def new_key(ic):  # remove the /q and /p layers, keep every other layer
    return '/'.join(p for p in ic.split('/') if p[:1] not in ('q', 'p'))


cases = [
    ("Ar",              "[Ar]"),
    ("neutral 13C-DME", "[13CH3]O[CH3]"),
    ("13C-DME cation",  "[13CH3][O+][CH3]"),
    ("12C-DME",         "COC"),
    ("ethanol",         "CCO"),
]
print(f"{'name':18s} {'inchi':40s} {'OLD key':30s} {'NEW key'}")
for name, s in cases:
    ic = inchi(s)
    print(f"{name:18s} {ic:40s} {old_key(ic):30s} {new_key(ic)}")


def keyed(fn, s):
    return fn(inchi(s))


print()
print("OLD: 13C-DME cation == 12C-DME neutral (the transmutation bug):",
      keyed(old_key, "[13CH3][O+][CH3]") == keyed(old_key, "COC"))
print("NEW: 13C-DME cation == 12C-DME neutral (must be False):",
      keyed(new_key, "[13CH3][O+][CH3]") == keyed(new_key, "COC"))
print("NEW: 13C-DME neutral == 13C-DME cation (charge independence kept):",
      keyed(new_key, "[13CH3]O[CH3]") == keyed(new_key, "[13CH3][O+][CH3]"))
print("NEW: Ar == Ar+ (electronic states still coincide):",
      keyed(new_key, "[Ar]") == keyed(new_key, "[Ar+]"))
