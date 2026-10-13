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
from rmgpy.kmc.generation_codec import roundtrip_generation_provider_factory


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
            roundtrip_generation_provider_factory() if args.mode == "roundtrip" else None
        ),
    )
    print(f"mode={args.mode} artifact={path}")
    print(f"artifact_id={artifact_id}")
    print(f"file_sha256={hashlib.sha256(Path(path).read_bytes()).hexdigest()}")
    print(f"wall_seconds={time.time() - started:.0f}")


if __name__ == "__main__":
    main()
