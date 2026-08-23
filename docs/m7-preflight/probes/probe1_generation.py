"""Generation-arrow probe: does a neutral organolithium generate a Cation_* reaction
producing a cation, and how is the electron carried (side + sign)?"""
import os, sys
from rmgpy import settings
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.species import Species

DB = "/home/alon/Code/RMG-database-plasma/input"
settings['database.directory'] = DB

CATION_FAMILIES = [
    'Cation_R_Recombination',
    'Cation_Li_Abstraction',
    'Cation_Addition_MultipleBond',
    'Cation_Addition_MultipleBond_Disprop',
    'Cation_NO_Substitution',
    'Cation_NO_Ring_Opening',
]

kdb = KineticsDatabase()
kdb.load(os.path.join(DB, 'kinetics'), families=CATION_FAMILIES, depositories=[])
print("loaded families:", sorted(kdb.families.keys()))

def spc(label, smiles):
    s = Species(label=label).from_smiles(smiles)
    s.generate_resonance_structures()
    return s

# Neutral organolithium seed (methyllithium) — analogous to training "CH3 + Li <=> CH3Li"
seeds = {
    'CH3Li': 'C[Li]',
    'C2H5Li': 'CC[Li]',
    'NH2Li': 'N[Li]',
    'LiH':   '[LiH]',
}

for label, smi in seeds.items():
    try:
        s = spc(label, smi)
    except Exception as e:
        print(f"\n### {label} ({smi}) FAILED to build species: {e!r}")
        continue
    print(f"\n### seed {label} ({smi})  charge={s.molecule[0].get_net_charge()}")
    try:
        rxns = kdb.generate_reactions_from_families([s], only_families=None)
    except Exception as e:
        print(f"  generate_reactions_from_families raised: {type(e).__name__}: {e}")
        continue
    if not rxns:
        print("  (no reactions generated)")
    for rxn in rxns:
        # Is any product a cation?
        prod_charges = [m.molecule[0].get_net_charge() for m in rxn.products]
        react_charges = [m.molecule[0].get_net_charge() for m in rxn.reactants]
        cation = any(c > 0 for c in prod_charges) or any(c > 0 for c in react_charges)
        e_meta = getattr(rxn, 'electrons', 0)
        fam = getattr(rxn, 'family', '?')
        fwd = getattr(rxn, 'is_forward', None)
        tag = "  *** CATION ***" if cation else ""
        print(f"  [{fam}] fwd={fwd} electrons_meta={e_meta} | "
              f"{' + '.join(str(r) for r in rxn.reactants)} => "
              f"{' + '.join(str(p) for p in rxn.products)}"
              f"  Rq={react_charges} Pq={prod_charges}{tag}")
