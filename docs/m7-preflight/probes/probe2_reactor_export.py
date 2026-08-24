"""Arrow-by-arrow probe on the ONE generated Cation reaction from neutral CH3Li:
   electron sign (item 2), reactor placement (item 3 break), export (item 4)."""
import os, traceback
from rmgpy import settings
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.species import Species

DB = "/home/alon/Code/RMG-database-plasma/input"
settings['database.directory'] = DB
CAT = ['Cation_R_Recombination']

kdb = KineticsDatabase()
kdb.load(os.path.join(DB, 'kinetics'), families=CAT, depositories=[])

s = Species(label='CH3Li').from_smiles('C[Li]'); s.generate_resonance_structures()
rxns = kdb.generate_reactions_from_families([s], only_families=None)
rxn = rxns[0]
print("REACTION:", rxn)
print("  family      :", rxn.family)
print("  is_forward  :", rxn.is_forward)
print("  electrons   :", getattr(rxn, 'electrons', 0))
print("  reactants   :", [(str(r), r.molecule[0].get_net_charge()) for r in rxn.reactants])
print("  products    :", [(str(p), p.molecule[0].get_net_charge()) for p in rxn.products])

# --- item 2: electron placement + sign via electron_balance.expand_electrons ---
print("\n== electron_balance.expand_electrons (export-side view) ==")
try:
    from rmgpy.electron_balance import expand_electrons, check_electron_balance
    # build a species_list that includes an electron pseudo-species
    e = Species(label='e-').from_adjacency_list("1 e u1 p0 c-1")
    slist = list(rxn.reactants) + list(rxn.products) + [e]
    R, P = expand_electrons(rxn, slist)
    print("  expand_electrons -> reactants:", [str(x) for x in R])
    print("  expand_electrons -> products :", [str(x) for x in P])
    equation = "%s <=> %s" % (' + '.join(str(x) for x in R), ' + '.join(str(x) for x in P))
    try:
        check_electron_balance(rxn, R, P, equation)
        print("  check_electron_balance: PASS (charge/electron balanced)")
    except Exception as ex:
        print("  check_electron_balance: FAIL ->", type(ex).__name__, "->", ex)
except Exception:
    print("  expand_electrons path raised:"); traceback.print_exc()

# --- item 3: reactor electron placement (resolve_electron_placement) ---
print("\n== electron_placement.resolve_electron_placement (reactor-side) ==")
try:
    from rmgpy import electron_placement
    e = Species(label='e-').from_adjacency_list("1 e u1 p0 c-1")
    slist = list(rxn.reactants) + list(rxn.products) + [e]
    view = electron_placement.resolve_electron_placement(rxn, slist)
    print("  RESOLVED OK ->", view)
except Exception as ex:
    print("  RAISED:", type(ex).__name__)
    print("  MESSAGE:", ex)

# --- item 4: chemkin export of the charged reaction ---
print("\n== chemkin.write_reaction_string (export-side) ==")
try:
    from rmgpy.chemkin import write_reaction_string
    e = Species(label='e-').from_adjacency_list("1 e u1 p0 c-1")
    slist = list(rxn.reactants) + list(rxn.products) + [e]
    txt = write_reaction_string(rxn, species_list=slist)
    print("  write_reaction_string ->", txt)
except Exception as ex:
    print("  RAISED:", type(ex).__name__, "->", ex)
