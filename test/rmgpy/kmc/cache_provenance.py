"""Reuse expensive public-RMG results only across unchanged generator code."""

import functools
import subprocess


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
