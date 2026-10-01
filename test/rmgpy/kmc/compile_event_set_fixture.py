"""Compile the real PS event set in an isolated hash-seeded interpreter."""

import argparse
from pathlib import Path

from rmgpy.data.rmg import RMGDatabase
from rmgpy.kmc.compiler import (
    EventSetCompiler,
    PS_FAMILY_CANDIDATES,
    PS_PROXY_UNITS,
    ps_proxy_set,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("database")
    parser.add_argument("output")
    args = parser.parse_args()
    database_path = Path(args.database)
    family_root = database_path / "input/kinetics/families"
    family_universe = sorted(
        path.name for path in family_root.iterdir() if (path / "groups.py").is_file()
    )

    database = RMGDatabase()
    database.load_kinetics(
        str(database_path / "input/kinetics"),
        reaction_libraries=[],
        seed_mechanisms=None,
        kinetics_families=list(PS_FAMILY_CANDIDATES),
        kinetics_depositories=["training"],
    )
    database.load_thermo(
        str(database_path / "input/thermo"),
        thermo_libraries=["primaryThermoLibrary"],
        depository=True,
    )
    proxies = ps_proxy_set(PS_PROXY_UNITS)
    active, excluded, reactions = EventSetCompiler.discover_family_reactions(
        database.kinetics,
        proxies,
        PS_FAMILY_CANDIDATES,
        family_universe=family_universe,
    )
    compiler = EventSetCompiler(
        database.kinetics,
        proxies,
        active,
        excluded_families=excluded,
        database_path=database_path,
        thermo_database=database.thermo,
        reaction_cache=reactions,
    )
    path, _ = compiler.write_artifact(args.output)
    print(path)


if __name__ == "__main__":
    main()
