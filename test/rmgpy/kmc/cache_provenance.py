"""Reuse expensive public-RMG results only across unchanged generator code."""

import functools
import hashlib
import json
import os
from pathlib import Path
import subprocess


def supplied_artifact():
    """Use an explicit artifact or fail validation; never fall back to compiling."""
    from rmgpy.kmc.compiler import compiler_source_hash, validate_artifact

    supplied = os.environ.get("RMG_KMC_ARTIFACT")
    if not supplied:
        return None
    path = Path(supplied)
    payload = path.read_bytes()
    assert path.stem == hashlib.sha256(payload).hexdigest()
    artifact = json.loads(payload)
    assert artifact["provenance"]["compiler_sources_sha256"] == compiler_source_hash()
    validate_artifact(artifact)
    return path, artifact


@functools.lru_cache(maxsize=None)
def generator_code_unchanged(repository, origin, current):
    changes = subprocess.check_output(
        [
            "git",
            "-C",
            str(repository),
            "diff",
            "--name-only",
            origin,
            current,
            "--",
            "rmgpy",
        ],
        text=True,
    ).splitlines()
    changes.extend(
        subprocess.check_output(
            ["git", "-C", str(repository), "diff", "--name-only", "--", "rmgpy"],
            text=True,
        ).splitlines()
    )
    return not any(not path.startswith("rmgpy/kmc/") for path in changes)
