"""Compile the PS artifact cold or through the encode/decode provider and print its id.

    python generation_parity_driver.py DATABASE OUTPUT_DIR {cold,roundtrip}

Both modes use the one product cold-compile procedure with identical inputs;
only the generation provider differs.  Equal artifact ids prove the codec is
transparent to every consumer of generation output.
"""

import argparse
import hashlib
import logging
import time
from pathlib import Path

from rmgpy.kmc.compile import CompileRequest, compile_artifact
from rmgpy.kmc.compiler import PS_FAMILY_CANDIDATES, PS_PROXY_UNITS
from rmgpy.kmc.generation_codec import DecodeLimits, roundtrip_generation_provider_factory

# Padded @5 contexts yield single-call payloads of hundreds of MB (about 150 KB
# per reaction); the library default (256 MiB) is deliberately smaller.
LIMITS = DecodeLimits(
    max_bytes=8 << 30, max_tokens=10**10, max_atoms=10**9, max_bonds=10**9,
    max_molecules=10**8, max_species=10**8, max_reactions=10**8, max_references=10**11,
    max_kinetics=10**8, max_thermo=10**8,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("database")
    parser.add_argument("output")
    parser.add_argument("mode", choices=("cold", "roundtrip"))
    args = parser.parse_args()
    logger = logging.getLogger("rmgpy.kmc.compiler")
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)
    logger.propagate = False
    started = time.time()
    path, artifact_id = compile_artifact(
        CompileRequest(
            database_path=str(Path(args.database).resolve()),
            output_dir=str(Path(args.output).resolve()),
            families=tuple(PS_FAMILY_CANDIDATES),
            proxy_units=PS_PROXY_UNITS,
        ),
        generation_provider_factory=(
            roundtrip_generation_provider_factory(LIMITS) if args.mode == "roundtrip" else None
        ),
    )
    print(f"mode={args.mode} artifact={path}")
    print(f"artifact_id={artifact_id}")
    print(f"file_sha256={hashlib.sha256(Path(path).read_bytes()).hexdigest()}")
    print(f"wall_seconds={time.time() - started:.0f}")


if __name__ == "__main__":
    main()
