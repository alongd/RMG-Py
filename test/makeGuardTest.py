#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@northeastern.edu) and the RMG Team            #
# (rmg_dev@mit.edu)                                                           #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
Regression lock for the Makefile's shared-environment guard.

`pip install -e .` rewrites the *environment's* editable-install record rather than
installing into the checkout that ran it. Because many worktrees of this repository
share one conda environment, the old default goal (`make` -> `pip install -e .`)
silently repointed every other worktree's ``import rmgpy`` at whichever tree ran it.

Testing that by letting real ``pip`` run is exactly the thing that must never happen,
and asserting only "``make`` fails" is worthless: on a machine whose editable-install
record is already broken, ``make`` fails for an unrelated reason and the assertion is
green before the guard exists.

So these tests run ``make`` inside a sandbox whose ``PATH`` puts stand-in ``python``
and ``pip`` executables ahead of the real ones. Nothing real is installed, built, or
removed; instead every invocation is recorded and classified, so the tests can ask the
only question that matters: *did this target reach the editable-install operation?*

The harness is kept honest by ``test_harness_detects_hazard_on_preguard_snapshot``,
which points it at the committed pre-guard Makefile and requires it to report the
hazard. That test is the RED-before, preserved as an executable artifact: if it ever
stops finding the hazard in the old Makefile, the harness has gone blind and the
green results below mean nothing.
"""

import hashlib
import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_DIR = Path(__file__).resolve().parent / "makeGuard"
LIVE_MAKEFILE = REPO_ROOT / "Makefile"
PREGUARD_MAKEFILE = FIXTURE_DIR / "preguard_Makefile.snapshot"

# Events the stand-ins emit when a recipe reaches something that would have
# mutated the shared environment.
HAZARD_EVENTS = ("PIP_INSTALL_REACHED", "PIP_UNINSTALL_REACHED", "PIP_OTHER_REACHED")

# Every target that reached `pip install -e .` or `pip uninstall` before the guard,
# whether directly or through a prerequisite.
HAZARDOUS_TARGETS = [
    "",  # the default goal, i.e. bare `make`
    "all",
    "install",
    ".installed",
    "clean",
    "scoop",
] + [f"eg{i}" for i in range(11)]

MAKE = shutil.which("make")

pytestmark = pytest.mark.skipif(MAKE is None, reason="`make` is not available on PATH")


class MakeRun:
    """Result of one sandboxed `make` invocation."""

    def __init__(self, returncode, stdout, stderr, events, argv):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        self.events = events
        self.argv = argv

    @property
    def output(self):
        return self.stdout + self.stderr

    @property
    def hazards(self):
        return [line for line in self.events if line.startswith(HAZARD_EVENTS)]


def _run_make(sandbox, makefile, target="", variables=(), extra_env=None, cwd=None):
    """
    Run `make` on `makefile` with stand-in `python`/`pip` ahead of the real ones.

    `sandbox` is a fresh directory; `cwd` defaults to a copy of the Makefile inside it,
    so nothing outside the sandbox is written. Pass `cwd` explicitly to exercise the
    guard from a real checkout (the Makefile there is used as-is).
    """
    bin_dir = sandbox / "bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    for name in ("python", "pip"):
        stand_in = bin_dir / name
        shutil.copy(FIXTURE_DIR / f"fake_{name}", stand_in)
        stand_in.chmod(0o755)

    if cwd is None:
        cwd = sandbox / "work"
        cwd.mkdir(parents=True, exist_ok=True)
        shutil.copy(makefile, cwd / "Makefile")

    argv_log = sandbox / "argv.log"
    event_log = sandbox / "event.log"
    argv_log.write_text("")
    event_log.write_text("")

    env = {
        "PATH": f"{bin_dir}{os.pathsep}{os.defpath}",
        "HOME": str(sandbox),
        "MAKE_GUARD_ARGV_LOG": str(argv_log),
        "MAKE_GUARD_EVENT_LOG": str(event_log),
    }
    env.update(extra_env or {})

    command = [MAKE, "--no-print-directory"]
    if target:
        command.append(target)
    command.extend(variables)

    completed = subprocess.run(
        command,
        cwd=str(cwd),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=120,
    )
    return MakeRun(
        completed.returncode,
        completed.stdout,
        completed.stderr,
        event_log.read_text().splitlines(),
        argv_log.read_text().splitlines(),
    )


def _manifest(directory):
    """sha256 of every file under `directory`, keyed by relative path."""
    manifest = {}
    for path in sorted(Path(directory).rglob("*")):
        if path.is_file():
            manifest[str(path.relative_to(directory))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return manifest


def _write_editable_link_state(site_packages, state, tmp_path):
    """
    Build a stand-in `site-packages` in one of the four states the guard must be
    indifferent to. The records are the pre-PEP-660 setuptools-develop pair: pip
    reads the `.egg-link`, `import` follows `easy-install.pth`.
    """
    site_packages.mkdir(parents=True, exist_ok=True)
    if state == "absent":
        return
    if state == "current_checkout":
        target = REPO_ROOT
    elif state == "other_worktree":
        target = tmp_path / "another-worktree"
        target.mkdir(parents=True, exist_ok=True)
    elif state == "stale_broken":
        target = tmp_path / "worktree-deleted-weeks-ago"  # deliberately never created
    else:
        raise ValueError(state)
    (site_packages / "reactionmechanismgenerator.egg-link").write_text(f"{target}\n.\n")
    (site_packages / "easy-install.pth").write_text(f"{target}\n")


# ---------------------------------------------------------------------------
# The RED-before, kept executable.
# ---------------------------------------------------------------------------


def test_harness_detects_hazard_on_preguard_snapshot(tmp_path):
    """
    Bare `make` on the pre-guard Makefile reaches `pip install -e .`.

    This is the observation that makes every other test in this file meaningful.
    Note the exit status: zero. The pre-guard default goal did not merely *try* to
    install, it would have *succeeded*, repointing the shared environment. Any
    apparent safety on a live machine came from the environment's own broken
    editable-install record, not from this Makefile.
    """
    run = _run_make(tmp_path / "sandbox", PREGUARD_MAKEFILE)

    assert any(e.startswith("PIP_INSTALL_REACHED") for e in run.events), (
        "harness failed to detect the known hazard in the pre-guard Makefile; "
        f"events were {run.events}"
    )
    assert any("-e" in line and "pip install" in line for line in run.argv)
    assert run.returncode == 0, "pre-guard default goal was expected to run to completion"


# ---------------------------------------------------------------------------
# Proofs 1, 2, 9: the guard refuses, names the replacement, and mutates nothing.
# ---------------------------------------------------------------------------


def test_bare_make_refuses_before_pip(tmp_path):
    run = _run_make(tmp_path / "sandbox", LIVE_MAKEFILE)

    assert run.returncode != 0
    assert run.hazards == [], f"guarded Makefile still reached pip: {run.hazards}"
    assert "Refusing to modify the shared RMG environment." in run.output


def test_guard_message_names_the_safe_replacement(tmp_path):
    run = _run_make(tmp_path / "sandbox", LIVE_MAKEFILE)

    assert "make build" in run.output
    assert "unsafe-install-shared-env" in run.output


def test_guard_runs_before_any_recipe_that_could_mutate(tmp_path):
    """The refusal is the *first* thing that happens: nothing at all is executed."""
    run = _run_make(tmp_path / "sandbox", LIVE_MAKEFILE)

    assert run.argv == [], f"default goal executed commands before refusing: {run.argv}"


# ---------------------------------------------------------------------------
# Proof 3: the four editable-install-link states.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("state", ["absent", "current_checkout", "other_worktree", "stale_broken"])
def test_guard_is_indifferent_to_the_editable_link_state(tmp_path, state):
    """
    The guard must not depend on the environment's editable-install record — that
    record being broken is the accident the guard replaces, and a guard keyed on a
    path would reproduce it. Here the record is planted in all four states; the
    refusal must be byte-identical, and the record itself must be untouched.
    """
    site_packages = tmp_path / "site-packages"
    _write_editable_link_state(site_packages, state, tmp_path)
    before = _manifest(site_packages)

    run = _run_make(
        tmp_path / "sandbox",
        LIVE_MAKEFILE,
        extra_env={"PYTHONPATH": str(site_packages)},
    )

    assert run.returncode != 0
    assert run.hazards == []
    assert "Refusing to modify the shared RMG environment." in run.output
    assert "make build" in run.output
    assert _manifest(site_packages) == before, "the guard touched the editable-install record"


def test_guard_refusal_is_identical_across_all_link_states(tmp_path):
    """Same message, same exit path, whatever the environment's record says."""
    outputs = {}
    for state in ["absent", "current_checkout", "other_worktree", "stale_broken"]:
        site_packages = tmp_path / state / "site-packages"
        _write_editable_link_state(site_packages, state, tmp_path / state)
        run = _run_make(
            tmp_path / state / "sandbox",
            LIVE_MAKEFILE,
            extra_env={"PYTHONPATH": str(site_packages)},
        )
        outputs[state] = (run.returncode, run.output)

    distinct = set(outputs.values())
    assert len(distinct) == 1, f"guard behaviour varied with the link state: {outputs}"


# ---------------------------------------------------------------------------
# Proof 3 (coverage), 8: every hazardous target, and the scoped clean.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("target", HAZARDOUS_TARGETS)
def test_hazardous_target_no_longer_reaches_pip(tmp_path, target):
    """
    The guard covers the hazardous *path*, not the literal command `make`. Each of
    these reached `pip install -e .` or `pip uninstall` before the change, directly
    or through a prerequisite. None may reach it now. Targets whose later steps fail
    in the sandbox (the `eg*` examples cannot find `examples/`) are still valid here:
    the assertion is about pip, not about the target completing.
    """
    run = _run_make(tmp_path / "sandbox", LIVE_MAKEFILE, target=target)

    assert run.hazards == [], f"`make {target}` reached the shared environment: {run.hazards}"


def test_clean_is_scoped_and_does_not_uninstall(tmp_path):
    """
    Proof 8: the scoped-clean route stays usable, and no longer removes the package
    from the environment other workers are importing from.
    """
    run = _run_make(tmp_path / "sandbox", LIVE_MAKEFILE, target="clean")

    assert run.returncode == 0
    assert run.hazards == []
    assert any("utilities.py clean" in line for line in run.argv), run.argv


# ---------------------------------------------------------------------------
# Proofs 4, 5: `make build` stays available and stays out of the environment.
# ---------------------------------------------------------------------------


def test_build_compiles_in_place_without_touching_the_environment(tmp_path):
    run = _run_make(tmp_path / "sandbox", LIVE_MAKEFILE, target="build")

    assert run.returncode == 0
    assert run.hazards == []
    assert any("setup.py" in line and "build_ext" in line and "--inplace" in line for line in run.argv), run.argv


def test_build_leaves_a_planted_editable_link_untouched(tmp_path):
    site_packages = tmp_path / "site-packages"
    _write_editable_link_state(site_packages, "current_checkout", tmp_path)
    before = _manifest(site_packages)

    run = _run_make(
        tmp_path / "sandbox",
        LIVE_MAKEFILE,
        target="build",
        extra_env={"PYTHONPATH": str(site_packages)},
    )

    assert run.returncode == 0
    assert _manifest(site_packages) == before


# ---------------------------------------------------------------------------
# The retained editable-install targets: conspicuous name, explicit opt-in.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("target", ["unsafe-install-shared-env", "unsafe-uninstall-shared-env"])
def test_editable_install_targets_require_the_opt_in(tmp_path, target):
    run = _run_make(tmp_path / "sandbox", LIVE_MAKEFILE, target=target)

    assert run.returncode != 0
    assert run.hazards == []
    assert "Refusing to modify the shared RMG environment." in run.output


@pytest.mark.parametrize("wrong", ["1", "true", "YES", "y", ""])
def test_opt_in_is_not_satisfied_by_a_near_miss(tmp_path, wrong):
    run = _run_make(
        tmp_path / "sandbox",
        LIVE_MAKEFILE,
        target="unsafe-install-shared-env",
        variables=[f"CONFIRM_SHARED_ENV_MUTATION={wrong}"],
    )

    assert run.returncode != 0
    assert run.hazards == []


def test_opt_in_reaches_the_install_when_explicitly_confirmed(tmp_path):
    """
    The escape hatch has to actually work, or the maintenance procedure the guard
    points people at is a dead end. Real pip is never involved: the stand-in records
    the call and returns success.
    """
    run = _run_make(
        tmp_path / "sandbox",
        LIVE_MAKEFILE,
        target="unsafe-install-shared-env",
        variables=["CONFIRM_SHARED_ENV_MUTATION=yes"],
    )

    assert run.returncode == 0
    assert any(e.startswith("PIP_INSTALL_REACHED") for e in run.events), run.events


# ---------------------------------------------------------------------------
# Proof 7: from a linked git worktree, against the real checkout.
# ---------------------------------------------------------------------------


def test_guard_holds_when_invoked_from_the_real_checkout(tmp_path):
    """
    Runs `make` with no arguments in the actual repository directory, which on the
    development machines this guard exists for is a linked git worktree. The
    stand-in `pip` is still ahead of the real one on PATH — if the guard ever
    regressed, this test must report it rather than perform the install it is
    testing for.
    """
    run = _run_make(tmp_path / "sandbox", LIVE_MAKEFILE, cwd=REPO_ROOT)

    assert run.returncode != 0
    assert run.hazards == []
    assert run.argv == []
    assert "make build" in run.output


def test_repository_is_a_linked_worktree_or_a_plain_checkout():
    """
    Records which of the two shapes proof 7 was observed in, so a green run from a
    plain clone is not mistaken for a green run from a linked worktree.
    """
    git_marker = REPO_ROOT / ".git"
    assert git_marker.exists()
    shape = "linked worktree" if git_marker.is_file() else "plain checkout"
    print(f"guard proof 7 observed in a {shape}: {REPO_ROOT}")


# ---------------------------------------------------------------------------
# The fixture that keeps the RED-before honest.
# ---------------------------------------------------------------------------


def test_preguard_snapshot_still_contains_the_hazard_it_documents():
    """If someone 'tidies' the snapshot, the RED-before stops proving anything."""
    text = PREGUARD_MAKEFILE.read_text()
    assert "pip install --no-build-isolation -vv -e ." in text
    assert text.count("all: check") == 1
