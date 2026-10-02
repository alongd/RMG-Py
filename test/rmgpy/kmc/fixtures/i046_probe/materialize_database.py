"""Materialize the existing probe allowlist at the dispatch pin, read-only."""

import argparse
import importlib.util
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
DATABASE = Path("/home/alon/Code/RMG-database")
DATABASE_SHA = "4a12d36fcdc193ede82c8d1ab5c1653495d445bc"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scratch", type=Path)
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location("i034", HERE.parent / "i034_probe/run_probe.py")
    prior = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prior)
    snapshot = prior.snapshot_database(DATABASE, args.scratch / "database")
    paths = subprocess.check_output(
        ["git", "-C", str(DATABASE), "ls-tree", "-r", "--name-only", DATABASE_SHA,
         "--", "input/kinetics/families"], text=True,
    ).splitlines()
    universe = sorted({Path(path).parent.name for path in paths if path.endswith("/groups.py")})
    (args.scratch / "family-universe.json").write_text(json.dumps(universe) + "\n")
    (args.scratch / "snapshot.json").write_text(json.dumps(snapshot, indent=2) + "\n")
    print(json.dumps(snapshot, sort_keys=True))


if __name__ == "__main__":
    main()
