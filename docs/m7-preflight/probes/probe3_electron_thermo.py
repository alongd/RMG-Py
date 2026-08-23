"""Is the electron-thermo break a config omission or a u1/u0 representation mismatch?"""
import os
from rmgpy import settings
settings['database.directory'] = "/home/alon/Code/RMG-database-plasma/input"
from rmgpy.data.thermo import ThermoDatabase
from rmgpy.species import Species

DB = "/home/alon/Code/RMG-database-plasma/input"
tdb = ThermoDatabase()
tdb.load(os.path.join(DB, 'thermo'),
         libraries=['electrocatThermo', 'LithiumPrimaryThermo', 'LithiumAdditionalThermo', 'primaryThermoLibrary'])

for lbl, adj in [('e- (u1)', '1 e u1 p0 c-1'), ('e- (u0)', '1 e u0 p0 c-1')]:
    s = Species(label=lbl).from_adjacency_list(adj)
    print(f"\n=== {lbl}: is_electron={s.is_electron()} ===")
    # direct library hit?
    try:
        lib_hit = tdb.get_thermo_data_from_libraries(s)
        print("  library thermo:", "HIT" if lib_hit else "MISS", (lib_hit[0].comment if lib_hit else ""))
    except Exception as ex:
        print("  library lookup raised:", type(ex).__name__, ex)
    # full path (may fall to groups)
    try:
        t = tdb.get_thermo_data(s)
        print("  get_thermo_data OK  H298=%.1f kJ/mol  source=%s" % (t.get_enthalpy(298)/1000.0, getattr(t,'comment','')[:60]))
    except Exception as ex:
        print("  get_thermo_data RAISED:", type(ex).__name__, str(ex)[:120])
