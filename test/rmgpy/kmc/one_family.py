"""Test support (I-165): one kinetics family, the thermo database and the
prepared rate rules, loaded once per family name (module-cached).  A test
that needs the RMG database loads exactly one family through ``family_db``."""

from pathlib import Path

DB = str(Path("/home/alon/runs/phase2b/database"))
_DB = {}


def available() -> bool:
    return Path(DB).is_dir()


def family_db(name):
    """One kinetics family + thermo + prepared rules (module-cached)."""
    if name not in _DB:
        from rmgpy.data.rmg import RMGDatabase
        from rmgpy.kmc.compiler import prepare_rate_rules
        db = RMGDatabase()
        db.load_kinetics(DB + "/input/kinetics", reaction_libraries=[], seed_mechanisms=None,
                         kinetics_families=[name], kinetics_depositories=["training"])
        db.load_thermo(DB + "/input/thermo", thermo_libraries=["primaryThermoLibrary"],
                       depository=True)
        prepare_rate_rules(db.kinetics, db.thermo, verbose=False)
        _DB[name] = db
    return _DB[name]
