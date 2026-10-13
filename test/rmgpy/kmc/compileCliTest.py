"""The one product compile procedure: hook, spawned child, CE parity."""

import copy
import json
import os
import signal
import subprocess
import sys
import textwrap
import time
from pathlib import Path

import pytest

from rmgpy.kmc.compiler import (
    EventSetCompiler,
    _production_padded_root,
    ps_proxy_set,
)

REPOSITORY = Path(__file__).resolve().parents[3]
DATABASE = Path(os.environ.get("RMG_DATABASE_PATH", "/nonexistent"))
needs_database = pytest.mark.skipif(
    not (DATABASE / "input/kinetics/families").is_dir(),
    reason="RMG_DATABASE_PATH is not a database checkout",
)
# A single family keeps the real compile small.
SUBSET = ["intra_H_migration"]


def _environment(**extra):
    environment = dict(os.environ)
    environment["PYTHONHASHSEED"] = "0"
    environment["PYTHONPATH"] = str(REPOSITORY) + os.pathsep + environment.get(
        "PYTHONPATH", ""
    )
    environment.pop("RMG_DATABASE_SHA", None)
    environment.update(extra)
    return environment


def _write_request_files(directory):
    families = directory / "families.json"
    families.write_text(json.dumps(SUBSET))
    return families


class _Boom(Exception):
    pass


class _ForbiddenDatabase:
    """A database whose own generator must never be reached."""

    families = {name: object() for name in SUBSET + ["fake"]}

    def __init__(self):
        self.direct_calls = 0

    def generate_reactions_from_families(self, *args, **kwargs):
        self.direct_calls += 1
        raise AssertionError("generation bypassed the explicit hook")


# --- every generation call site goes through the hook -----------------------


def test_discovery_routes_every_generation_call_through_the_hook():
    database = _ForbiddenDatabase()
    calls = []

    def provider(reactants, products=None, only_families=None, resonance=True):
        calls.append(only_families)
        return []

    EventSetCompiler.discover_family_reactions(
        database, ps_proxy_set()[:2], SUBSET, generation_provider=provider
    )
    assert calls
    assert database.direct_calls == 0


def test_per_proxy_generation_routes_through_the_hook():
    database = _ForbiddenDatabase()
    calls = []

    def provider(reactants, products=None, only_families=None, resonance=True):
        calls.append(only_families)
        return []

    proxy = ps_proxy_set()[0]
    compiler = EventSetCompiler(
        database, [proxy], SUBSET, generation_provider=provider
    )
    assert list(compiler._generate(compiler.proxies[0])) == []
    assert calls == [[family] for family in SUBSET]
    assert database.direct_calls == 0


def test_padded_root_regeneration_routes_through_the_hook():
    from rmgpy.data.kinetics.family import TemplateReaction
    from rmgpy.molecule.molecule import Molecule
    from rmgpy.species import Species

    reactant = Species(molecule=[Molecule(smiles="[CH3]")])
    product = Species(molecule=[Molecule(smiles="[CH3]")])
    expected = TemplateReaction(
        reactants=[reactant], products=[product], family="fake", degeneracy=1.0
    )
    derived = copy.deepcopy(expected)
    derived.degeneracy = 7.0
    seen = []

    def provider(reactants, products=None, only_families=None, resonance=True):
        seen.append((only_families, resonance))
        return [derived]

    database = _ForbiddenDatabase()
    matched, count = _production_padded_root(
        database, expected, generation_provider=provider
    )
    assert matched.degeneracy == 7.0 and count == 1
    assert seen == [(["fake"], False)]
    assert database.direct_calls == 0


# --- real database: the hook sees real calls and its errors abort ------------


@needs_database
def test_real_compile_hook_receives_calls_and_its_error_aborts(tmp_path):
    from rmgpy.kmc.compile import CompileRequest, compile_artifact

    request = CompileRequest(
        database_path=str(DATABASE),
        output_dir=str(tmp_path / "out"),
        families=tuple(SUBSET),
        family_universe=tuple(SUBSET),
    )
    count = []
    defaults = []

    def factory(default):
        defaults.append(default)

        def provider(reactants, products=None, only_families=None, resonance=True):
            count.append(only_families)
            if len(count) == 3:
                raise _Boom("hook says stop")
            return default(reactants, products, only_families, resonance)

        return provider

    with pytest.raises(_Boom, match="hook says stop"):
        compile_artifact(request, generation_provider_factory=factory)
    assert len(count) == 3
    assert len(defaults) == 1 and callable(defaults[0])
    assert not list((tmp_path / "out").glob("*"))


# --- CLI versus CE ------------------------------------------------------------


def _start_cli(output, families):
    return subprocess.Popen(
        [
            sys.executable, "-m", "rmgpy.kmc.compile", str(DATABASE), str(output),
            "--families", str(families), "--family-universe", str(families),
            "--database-sha", "parity-sha",
        ],
        cwd=REPOSITORY, env=_environment(), stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, text=True,
    )


def _start_ce(tmp_path, output, families, cache_root):
    script = tmp_path / "run_ce.py"
    script.write_text(textwrap.dedent(f"""
        import json, sys
        sys.path.insert(0, {str(REPOSITORY / "test/rmgpy/kmc")!r})
        import compile_event_set_fixture as driver
        driver.PS_FAMILY_CANDIDATES = tuple(json.load(open({str(families)!r})))
        sys.argv = ["ce", {str(DATABASE)!r}, {str(output)!r},
                    "--family-universe", {str(families)!r},
                    "--database-sha", "parity-sha"]
        driver.main()
    """))
    # The fixture derives its cache identity from the working repository.
    return subprocess.Popen(
        [sys.executable, str(script)], cwd=REPOSITORY,
        env=_environment(RMG_KMC_CACHE_ROOT=str(cache_root)),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )


@needs_database
@pytest.mark.skipif(
    os.environ.get("RMG_KMC_SLOW") != "1", reason="set RMG_KMC_SLOW=1"
)
def test_cli_and_ce_cold_compile_are_byte_identical(tmp_path):
    families = _write_request_files(tmp_path)
    cli_out, ce_out = tmp_path / "cli", tmp_path / "ce"
    cli = _start_cli(cli_out, families)
    ce = _start_ce(tmp_path, ce_out, families, tmp_path / "empty-cache")
    cli_stdout, cli_stderr = cli.communicate()
    ce_stdout, ce_stderr = ce.communicate()
    assert cli.returncode == 0, cli_stderr[-3000:]
    assert ce.returncode == 0, ce_stderr[-3000:]
    cli_artifact = next(cli_out.glob("*.json"))
    ce_artifact = next(ce_out.glob("*.json"))
    print(f"artifact ids: cli={cli_artifact.stem} ce={ce_artifact.stem}")
    assert cli_artifact.name == ce_artifact.name
    assert cli_artifact.read_bytes() == ce_artifact.read_bytes()
    assert cli_artifact.stem in cli_stdout


# --- process isolation --------------------------------------------------------


def _children_of(pid):
    """Direct children that are spawned compile children (not the tracker)."""
    children = []
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            stat = (entry / "stat").read_text()
        except OSError:
            continue
        fields = stat[stat.rindex(")") + 2:].split()
        try:
            command = (entry / "cmdline").read_bytes()
        except OSError:
            continue
        if int(fields[1]) == pid and b"spawn_main" in command:
            children.append(int(entry.name))
    return children


def _alive(pid):
    try:
        stat = Path(f"/proc/{pid}/stat").read_text()
    except OSError:
        return False
    return stat[stat.rindex(")") + 2:].split()[0] != "Z"


@needs_database
def test_parent_sigkill_reaps_the_compile_child(tmp_path):
    families = _write_request_files(tmp_path)
    parent = subprocess.Popen(
        [
            sys.executable, "-m", "rmgpy.kmc.compile", str(DATABASE),
            str(tmp_path / "out"), "--families", str(families),
        ],
        cwd=REPOSITORY, env=_environment(),
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    try:
        deadline = time.time() + 120
        children = []
        while time.time() < deadline and not children:
            children = _children_of(parent.pid)
            time.sleep(0.2)
        assert children, "the CLI never spawned a compile child"
        time.sleep(3)  # let the child get properly under way
        assert all(_alive(pid) for pid in children)
        parent.send_signal(signal.SIGKILL)
        parent.wait()
        deadline = time.time() + 10
        while time.time() < deadline and any(_alive(pid) for pid in children):
            time.sleep(0.2)
        assert not any(_alive(pid) for pid in children), "orphaned compile child"
    finally:
        parent.kill()
        for pid in children:
            try:
                os.kill(pid, signal.SIGKILL)
            except OSError:
                pass


def test_child_exception_gives_nonzero_exit_traceback_and_no_artifact(tmp_path):
    families = _write_request_files(tmp_path)
    output = tmp_path / "out"
    result = subprocess.run(
        [
            sys.executable, "-m", "rmgpy.kmc.compile",
            str(tmp_path / "no-such-database"), str(output),
            "--families", str(families),
        ],
        cwd=REPOSITORY, env=_environment(), capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "Traceback" in result.stderr
    assert "compile child failed" in result.stderr
    assert not output.exists() or not list(output.iterdir())


# --- importing has no side effects -------------------------------------------


def test_import_loads_no_database_and_mutates_no_global(tmp_path):
    script = tmp_path / "probe.py"
    script.write_text(textwrap.dedent("""
        import gc, logging, os, sys
        import rmgpy.data.rmg, rmgpy.kmc.compiler as compiler
        before_env = dict(os.environ)
        before_globals = dict(vars(compiler))
        logger = logging.getLogger("rmgpy.kmc.compiler")
        before_logger = (list(logger.handlers), logger.level, logger.propagate)
        before_root = (list(logging.getLogger().handlers), logging.getLogger().level)
        import rmgpy.kmc.compile
        assert dict(os.environ) == before_env
        assert dict(vars(compiler)) == before_globals
        assert (list(logger.handlers), logger.level, logger.propagate) == before_logger
        assert (list(logging.getLogger().handlers), logging.getLogger().level) == before_root
        assert not [o for o in gc.get_objects()
                    if type(o) is rmgpy.data.rmg.RMGDatabase]
        print("clean")
    """))
    result = subprocess.run(
        [sys.executable, str(script)], cwd=REPOSITORY, env=_environment(),
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr[-2000:]
    assert result.stdout.strip().endswith("clean")
