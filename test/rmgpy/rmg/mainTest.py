#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
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

import array
import builtins
import errno
import fcntl
import logging
import os
import resource
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import textwrap
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from rmgpy import get_path, settings
import rmgpy.util as util_module
from rmgpy.data.rmg import RMGDatabase
from rmgpy.exceptions import InputError
from rmgpy.rmg import input as rmg_input
from rmgpy.rmg.main import RMG, RMG_Memory, initialize_log, make_profile_graph
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
from rmgpy.species import Species

originalPath = get_path()


PROGRAMMATIC_INVALID_DEPOSITORIES = [["all"], ["default"], [""], [None], [[]]]
NUMPY_INVALID_DEPOSITORIES = [
    np.array(["all"]),
    np.array(["default"]),
    np.array([""]),
    np.array([None]),
    np.array([[]]),
]


def _minimal_programmatic_rmg(tmp_path, kinetics_depositories):
    rmg = RMG(output_directory=str(tmp_path))
    rmg.database_directory = str(tmp_path)
    rmg.thermo_libraries = []
    rmg.transport_libraries = []
    rmg.reaction_libraries = []
    rmg.statmech_libraries = []
    rmg.seed_mechanisms = []
    rmg.kinetics_families = []
    rmg.kinetics_depositories = kinetics_depositories
    rmg.adsorption_groups = "adsorptionPt111"
    rmg.initial_species = []
    rmg.reaction_systems = []
    rmg.reaction_model = type("Model", (), {"vibrational_manifolds": []})()
    rmg.simulator_settings_list = [SimulatorSettings()]
    rmg.model_settings_list = [ModelSettings(tol_move_to_core=0.1)]
    return rmg


@pytest.mark.parametrize("selection", PROGRAMMATIC_INVALID_DEPOSITORIES)
def test_load_database_validates_programmatic_depositories(monkeypatch, tmp_path, selection):
    """Directly configured depository lists obey the same contract as input decks."""
    class Database:
        def load(self, **kwargs):
            raise AssertionError("invalid depositories reached the database loader")

    import rmgpy.rmg.main as main
    monkeypatch.setattr(main, "RMGDatabase", Database)
    rmg = _minimal_programmatic_rmg(tmp_path, selection)

    with pytest.raises(InputError):
        rmg.load_database()


@pytest.mark.parametrize("selection", NUMPY_INVALID_DEPOSITORIES)
def test_depository_directives_reject_numpy_arrays(monkeypatch, tmp_path, selection):
    """Array equality cannot turn a non-string value into a directive."""
    class Database:
        def load(self, **kwargs):
            raise AssertionError("invalid depositories reached the database loader")

    import rmgpy.rmg.main as main
    monkeypatch.setattr(main, "RMGDatabase", Database)
    rmg = _minimal_programmatic_rmg(tmp_path, selection)

    with pytest.raises(InputError):
        rmg_input.normalize_kinetics_depositories(selection)
    with pytest.raises(InputError):
        rmg.load_database()

    output = tmp_path / "invalid-array.py"
    with pytest.raises(InputError):
        rmg.save_input(str(output))
    assert not output.exists()


@pytest.mark.parametrize("selection", [["NIST"], "all", None])
def test_load_database_accepts_programmatic_depositories(monkeypatch, tmp_path, selection):
    """Documented programmatic selections still reach the database unchanged."""
    class Database:
        def load(self, **kwargs):
            self.load_kwargs = kwargs

        kinetics = type("Kinetics", (), {"families": {}, "libraries": {}})()
        thermo = type("Thermo", (), {"libraries": {}, "adsorption_groups": None})()
        forbidden_structures = type("Forbidden", (), {"entries": {}})()

    import rmgpy.rmg.main as main
    monkeypatch.setattr(main, "RMGDatabase", Database)
    rmg = _minimal_programmatic_rmg(tmp_path, selection)

    rmg.load_database()

    assert rmg.database.load_kwargs["kinetics_depositories"] == selection
    assert rmg.kinetics_depositories == selection
    if isinstance(selection, list):
        assert rmg.database.load_kwargs["kinetics_depositories"] is not selection


@pytest.mark.parametrize("selection", PROGRAMMATIC_INVALID_DEPOSITORIES)
def test_uncertainty_reload_validates_programmatic_depositories(monkeypatch, tmp_path, selection):
    """The uncertainty family reload refuses malformed programmatic selections."""
    class Kinetics:
        families = {}

        def load_families(self, path, families, depositories):
            raise AssertionError("invalid depositories reached the family loader")

    class Uncertainty:
        def __init__(self, output_directory):
            pass

        def extract_sources_from_model(self):
            pass

    import rmgpy.rmg.main as main
    monkeypatch.setattr(main, "Uncertainty", Uncertainty)
    rmg = _minimal_programmatic_rmg(tmp_path, selection)
    rmg.database = type("Database", (), {"kinetics": Kinetics(), "thermo": object()})()
    rmg.reaction_model.core = type("Core", (), {"species": [], "reactions": []})()
    rmg.uncertainty = {"global": False, "local": True, "uncorrelated": False, "correlated": False}

    with pytest.raises(InputError):
        rmg.run_uncertainty_analysis()


@pytest.mark.parametrize("selection", [["NIST"], "all", None])
def test_uncertainty_reload_accepts_programmatic_depositories(monkeypatch, tmp_path, selection):
    """The uncertainty reload preserves valid programmatic selections and list ownership."""
    class Kinetics:
        families = {}
        forwarded = None

        def load_families(self, path, families, depositories):
            self.forwarded = depositories

    class Uncertainty:
        def __init__(self, output_directory):
            pass

        def extract_sources_from_model(self):
            pass

    import rmgpy.rmg.main as main
    monkeypatch.setattr(main, "Uncertainty", Uncertainty)
    rmg = _minimal_programmatic_rmg(tmp_path, selection)
    rmg.database = type("Database", (), {"kinetics": Kinetics(), "thermo": object()})()
    rmg.reaction_model.core = type("Core", (), {"species": [], "reactions": []})()
    rmg.uncertainty = {"global": False, "local": True, "uncorrelated": False, "correlated": False}

    rmg.run_uncertainty_analysis()

    assert rmg.database.kinetics.forwarded == selection
    assert rmg.kinetics_depositories == selection
    if isinstance(selection, list):
        assert rmg.database.kinetics.forwarded is not selection


@pytest.mark.parametrize("selection", PROGRAMMATIC_INVALID_DEPOSITORIES)
def test_save_input_validates_programmatic_depositories(tmp_path, selection):
    """A refused programmatic selection cannot damage an existing input deck."""
    rmg = _minimal_programmatic_rmg(tmp_path, ["NIST"])
    output = tmp_path / "existing.py"
    rmg.save_input(str(output))
    original = output.read_bytes()
    rmg.kinetics_depositories = selection

    with pytest.raises(InputError):
        rmg.save_input(str(output))

    assert output.read_bytes() == original


def test_save_input_encoding_refusal_preserves_existing_deck(tmp_path):
    """A UTF-8 encoding refusal cannot damage an existing input deck."""
    rmg = _minimal_programmatic_rmg(tmp_path, ["NIST"])
    species = Species(label="N2").from_smiles("N#N")
    rmg.initial_species = [species]
    output = tmp_path / "existing.py"
    rmg.save_input(str(output))
    original = output.read_bytes()
    species.label = "\ud835\udc41\u2082"

    with pytest.raises(UnicodeEncodeError):
        rmg.save_input(str(output))

    assert output.read_bytes() == original


@pytest.mark.parametrize(("species_count", "failure_phase"), [(1, "close"), (12, "write")])
def test_save_input_io_failure_preserves_destinations(tmp_path, species_count, failure_phase):
    """Real close and write failures cannot replace or create an input deck."""
    rmg = _minimal_programmatic_rmg(tmp_path, ["NIST"])
    rmg.initial_species = [
        Species(label=f"alkane_{index}").from_smiles("C" * index)
        for index in range(1, species_count + 1)
    ]
    existing = tmp_path / "existing.py"
    absent = tmp_path / "absent.py"
    rmg.save_input(str(existing))
    original = existing.read_bytes()
    assert len(original) < 8192 if failure_phase == "close" else len(original) > 8192
    original_files = set(tmp_path.iterdir())

    def fail_save(target):
        old_limit = resource.getrlimit(resource.RLIMIT_FSIZE)
        old_signal = signal.signal(signal.SIGXFSZ, signal.SIG_IGN)

        try:
            resource.setrlimit(resource.RLIMIT_FSIZE, (512, old_limit[1]))
            with pytest.raises(OSError) as raised:
                rmg.save_input(str(target))
            assert raised.value.errno == errno.EFBIG
        finally:
            resource.setrlimit(resource.RLIMIT_FSIZE, old_limit)
            signal.signal(signal.SIGXFSZ, old_signal)

    fail_save(existing)
    assert existing.read_bytes() == original

    fail_save(absent)
    assert not absent.exists()
    assert set(tmp_path.iterdir()) == original_files


def test_save_input_removes_failed_dangling_symlink_target(tmp_path):
    rmg = _minimal_programmatic_rmg(tmp_path, ["NIST"])
    link = tmp_path / "input.py"
    target = tmp_path / "target.py"
    link.symlink_to(target.name)
    old_limit = resource.getrlimit(resource.RLIMIT_FSIZE)
    old_signal = signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
    try:
        resource.setrlimit(resource.RLIMIT_FSIZE, (4, old_limit[1]))
        with pytest.raises(OSError) as raised:
            rmg.save_input(str(link))
        assert raised.value.errno == errno.EFBIG
    finally:
        resource.setrlimit(resource.RLIMIT_FSIZE, old_limit)
        signal.signal(signal.SIGXFSZ, old_signal)

    assert link.is_symlink()
    assert not target.exists()
    assert list(tmp_path.glob(".rmg-*")) == []


@pytest.mark.parametrize("link_kind", ["relative", "absolute", "dangling"])
def test_save_input_preserves_symlink_and_updates_target(tmp_path, link_kind):
    """Saving through a symlink updates its target without replacing the link."""
    rmg = _minimal_programmatic_rmg(tmp_path, ["NIST"])
    target = tmp_path / "target.py"
    link = tmp_path / "input.py"
    if link_kind != "dangling":
        target.write_bytes(b"old target\n")
    link_target = str(target) if link_kind == "absolute" else target.name
    os.symlink(link_target, link)

    rmg.save_input(str(link))

    assert link.is_symlink()
    assert target.exists()
    assert b"kineticsDepositories = ['NIST']" in target.read_bytes()


@pytest.mark.parametrize("mode", [0o644, 0o664, 0o600, 0o755])
def test_save_input_preserves_existing_mode_and_owner(tmp_path, mode):
    """In-place rewriting retains the destination's mode and owner."""
    rmg = _minimal_programmatic_rmg(tmp_path, ["NIST"])
    output = tmp_path / "input.py"
    output.write_bytes(b"old deck\n")
    os.chmod(output, mode)
    original_owner = (output.stat().st_uid, output.stat().st_gid)

    rmg.save_input(str(output))

    assert output.stat().st_mode & 0o777 == mode
    assert (output.stat().st_uid, output.stat().st_gid) == original_owner


def test_save_input_preserves_access_acl_and_extended_attributes(tmp_path):
    """In-place rewriting preserves access ACLs and extended attributes."""
    if not hasattr(os, "setxattr") or not shutil.which("setfacl") or not shutil.which("getfacl"):
        pytest.skip("filesystem xattrs or ACL command-line tools are unavailable")
    rmg = _minimal_programmatic_rmg(tmp_path, ["NIST"])
    output = tmp_path / "input.py"
    output.write_bytes(b"old deck\n")
    os.setxattr(output, "user.rmg-provenance", b"calibration-v1")
    subprocess.run(["setfacl", "-m", f"u:{os.getuid()}:rw", str(output)], check=True)
    attributes = {name: os.getxattr(output, name) for name in os.listxattr(output)}
    access_acl = subprocess.check_output(["getfacl", "-cp", str(output)])

    rmg.save_input(str(output))

    assert {name: os.getxattr(output, name) for name in os.listxattr(output)} == attributes
    assert subprocess.check_output(["getfacl", "-cp", str(output)]) == access_acl


@pytest.mark.parametrize("file_umask", [0o022, 0o002])
def test_save_input_new_file_obeys_umask(tmp_path, file_umask):
    """A newly saved deck gets the normal create mode, not forced 0600."""
    rmg = _minimal_programmatic_rmg(tmp_path, ["NIST"])
    output = tmp_path / "input.py"
    previous_umask = os.umask(file_umask)
    try:
        rmg.save_input(str(output))
    finally:
        os.umask(previous_umask)

    assert output.stat().st_mode & 0o777 == 0o666 & ~file_umask


def test_save_input_refuses_read_only_existing_target(tmp_path):
    """A writable directory cannot authorize replacing a read-only deck."""
    rmg = _minimal_programmatic_rmg(tmp_path, ["NIST"])
    output = tmp_path / "input.py"
    output.write_bytes(b"old deck\n")
    os.chmod(output, 0o444)
    original = output.read_bytes()

    with pytest.raises(PermissionError, match="not writable"):
        rmg.save_input(str(output))

    assert output.read_bytes() == original
    assert output.stat().st_mode & 0o777 == 0o444
    assert set(tmp_path.iterdir()) == {output}


def test_save_input_updates_every_hard_link(tmp_path):
    """Saving rewrites the existing inode shared by every hard link."""
    rmg = _minimal_programmatic_rmg(tmp_path, ["NIST"])
    output = tmp_path / "input.py"
    alias = tmp_path / "canonical.py"
    output.write_bytes(b"old deck\n")
    os.link(output, alias)
    rmg.save_input(str(output))

    assert os.path.samefile(output, alias)
    assert output.read_bytes() == alias.read_bytes()
    assert b"kineticsDepositories = ['NIST']" in alias.read_bytes()


def test_save_input_follows_symlink_parent_before_dotdot(tmp_path):
    """The kernel, not lexical normalization, determines symlink/.. traversal."""
    rmg = _minimal_programmatic_rmg(tmp_path, ["NIST"])
    outside = tmp_path / "outside"
    child = outside / "child"
    child.mkdir(parents=True)
    canonical = outside / "input.py"
    canonical.write_bytes(b"old canonical\n")
    sibling = tmp_path / "input.py"
    sibling.write_bytes(b"old sibling\n")
    link = tmp_path / "link"
    link.symlink_to(outside.name + "/child")

    rmg.save_input(str(link / ".." / "input.py"))

    assert b"kineticsDepositories = ['NIST']" in canonical.read_bytes()
    assert sibling.read_bytes() == b"old sibling\n"


@pytest.mark.parametrize("flag", ["d", "A", "S"])
def test_save_input_preserves_inode_flags_and_atime(tmp_path, flag):
    """In-place saves retain inode policy and the pre-save access time."""
    if shutil.which("chattr") is None:
        pytest.skip("chattr is unavailable on this platform")
    rmg = _minimal_programmatic_rmg(tmp_path, ["NIST"])
    output = tmp_path / "input.py"
    output.write_bytes(b"old deck\n")
    result = subprocess.run(["chattr", "+" + flag, str(output)], capture_output=True, text=True)
    if result.returncode:
        pytest.skip("filesystem refused chattr +{}: {}".format(flag, result.stderr.strip()))
    old_time = 946684800000000000
    os.utime(output, ns=(old_time, old_time))

    def inode_flags():
        value = array.array("L", [0])
        with output.open("rb") as stream:
            fcntl.ioctl(stream.fileno(), 0x80086601, value, True)
        return value[0]

    before_flags = inode_flags()
    before_atime = output.stat().st_atime_ns
    try:
        rmg.save_input(str(output))
        assert inode_flags() == before_flags
        assert output.stat().st_atime_ns == before_atime
    finally:
        subprocess.run(["chattr", "-" + flag, str(output)], check=True)


def test_save_input_restores_existing_bytes_after_midwrite_enospc(monkeypatch, tmp_path):
    """A partial write to the original is restored from its durable backup."""
    rmg = _minimal_programmatic_rmg(tmp_path, ["NIST"])
    output = tmp_path / "input.py"
    rmg.save_input(str(output))
    original = output.read_bytes()
    real_open = builtins.open
    failed = False

    class FullDiskWriter:
        def __init__(self, stream):
            self.stream = stream

        def __enter__(self):
            return self

        def __exit__(self, exception_type, exception, traceback):
            self.stream.close()

        def __getattr__(self, name):
            return getattr(self.stream, name)

        def write(self, data):
            nonlocal failed
            failed = True
            self.stream.write(data[: max(1, len(data) // 2)])
            self.stream.flush()
            raise OSError(errno.ENOSPC, "simulated full disk")

    def open_with_full_disk(path, mode="r", *args, **kwargs):
        nonlocal failed
        stream = real_open(path, mode, *args, **kwargs)
        if not failed and mode == "wb" and os.fspath(path) == str(output):
            return FullDiskWriter(stream)
        return stream

    monkeypatch.setattr(builtins, "open", open_with_full_disk)
    with pytest.raises(OSError) as raised:
        rmg.save_input(str(output))

    assert raised.value.errno == errno.ENOSPC
    assert failed
    assert output.read_bytes() == original
    assert set(tmp_path.iterdir()) == {output}


def test_save_input_retains_backup_when_restore_fails(monkeypatch, tmp_path):
    """A failed restoration reports and retains the durable backup."""
    rmg = _minimal_programmatic_rmg(tmp_path, ["NIST"])
    output = tmp_path / "input.py"
    rmg.save_input(str(output))
    original = output.read_bytes()
    real_open = builtins.open
    write_attempts = 0

    class FailingWriter:
        def __init__(self, stream):
            self.stream = stream

        def __enter__(self):
            return self

        def __exit__(self, exception_type, exception, traceback):
            self.stream.close()

        def __getattr__(self, name):
            return getattr(self.stream, name)

        def write(self, data):
            self.stream.write(data[:1])
            self.stream.flush()
            raise OSError(errno.ENOSPC, "simulated full disk")

    def open_with_failed_restore(path, mode="r", *args, **kwargs):
        nonlocal write_attempts
        if mode == "wb" and os.fspath(path) == str(output):
            write_attempts += 1
            if write_attempts == 2:
                raise OSError(errno.EIO, "simulated restore failure")
            return FailingWriter(real_open(path, mode, *args, **kwargs))
        return real_open(path, mode, *args, **kwargs)

    monkeypatch.setattr(builtins, "open", open_with_failed_restore)
    with pytest.raises(OSError, match="backup was retained") as raised:
        rmg.save_input(str(output))

    backups = list(tmp_path.glob(".rmg-backup-*"))
    assert write_attempts == 2
    assert len(backups) == 1
    assert str(backups[0]) in str(raised.value)
    assert backups[0].read_bytes() == original


@pytest.mark.parametrize("kind", ["directory", "fifo", "self-loop", "two-loop"])
def test_save_input_refuses_nonregular_and_unresolved_destinations(tmp_path, kind):
    """Saving cannot destroy special files or unresolved symbolic-link cycles."""
    rmg = _minimal_programmatic_rmg(tmp_path, ["NIST"])
    output = tmp_path / "input.py"
    if kind == "directory":
        output.mkdir()
    elif kind == "fifo":
        os.mkfifo(output)
    elif kind == "self-loop":
        output.symlink_to(output.name)
    else:
        other = tmp_path / "other.py"
        output.symlink_to(other.name)
        other.symlink_to(output.name)
    original_entries = {path.name: path.lstat().st_mode for path in tmp_path.iterdir()}

    with pytest.raises(OSError):
        rmg.save_input(str(output))

    assert {path.name: path.lstat().st_mode for path in tmp_path.iterdir()} == original_entries
    if kind == "fifo":
        assert stat.S_ISFIFO(output.lstat().st_mode)
    if "loop" in kind:
        assert output.is_symlink()


@pytest.mark.parametrize("kind", ["self-loop", "two-loop"])
def test_save_input_refuses_loop_hidden_behind_dotdot_path(tmp_path, kind):
    """A dangling spelling that reaches a loop is not mistaken for a new file."""
    rmg = _minimal_programmatic_rmg(tmp_path, ["NIST"])
    output = tmp_path / "input.py"
    if kind == "self-loop":
        output.symlink_to(output.name)
    else:
        other = tmp_path / "other.py"
        output.symlink_to(other.name)
        other.symlink_to(output.name)
    alias = tmp_path / "alias.py"
    alias.symlink_to("missing/../input.py")
    before = {path.name: os.readlink(path) for path in tmp_path.iterdir()}

    with pytest.raises(OSError):
        rmg.save_input(str(alias))

    assert {path.name: os.readlink(path) for path in tmp_path.iterdir()} == before


def test_save_input_refuses_trailing_separator_without_changing_target(tmp_path):
    """Canonicalizing a symlink target must not erase an invalid path suffix."""
    rmg = _minimal_programmatic_rmg(tmp_path, ["NIST"])
    output = tmp_path / "input.py"
    output.write_bytes(b"old deck\n")
    original = output.read_bytes()

    with pytest.raises(OSError):
        rmg.save_input(str(output) + os.sep)

    assert output.read_bytes() == original
    assert set(tmp_path.iterdir()) == {output}


def test_save_input_accepts_maximum_length_basename(tmp_path):
    """Temporary names must not extend an already maximum-length basename."""
    rmg = _minimal_programmatic_rmg(tmp_path, ["NIST"])
    output = tmp_path / ("x" * 252 + ".py")

    rmg.save_input(str(output))

    assert b"kineticsDepositories = ['NIST']" in output.read_bytes()


@pytest.mark.parametrize("signum", [signal.SIGINT, signal.SIGTERM])
def test_save_input_signal_sweep_preserves_complete_destination(tmp_path, signum):
    """Every save statement settles handles and backups before delivering signals."""
    rmg = _minimal_programmatic_rmg(tmp_path, ["NIST"])
    save_names = {
        "_write_input_file_atomically",
        "write_file_with_backup",
        "_write_files_with_backups",
        "_prepare_save",
        "_open_backup",
        "_rewrite_save",
        "_restore_save",
        "_remove_backup",
        "__enter__",
        "__exit__",
    }
    reached = set()

    def collect(frame, event, arg):
        if event == "line" and frame.f_code.co_name in save_names:
            reached.add((frame.f_code, frame.f_lineno))
        return collect

    seed = tmp_path / "seed.py"
    seed.write_bytes(b"old deck\n")
    try:
        sys.settrace(collect)
        rmg.save_input(str(seed))
    finally:
        sys.settrace(None)

    class SaveSignal(BaseException):
        pass

    for code, line in sorted(reached, key=lambda item: (item[0].co_name, item[1])):
        directory = tmp_path / (code.co_name + "-" + str(line))
        directory.mkdir()
        output = directory / "input.py"
        original = b"old deck\n"
        output.write_bytes(original)
        fired = False
        error = None

        def interrupt_at_line(frame, event, arg):
            nonlocal fired
            if not fired and frame.f_code is code and event == "line" and frame.f_lineno == line:
                fired = True
                os.kill(os.getpid(), signum)
            return interrupt_at_line

        previous_handler = signal.signal(signum, lambda caught, frame: (_ for _ in ()).throw(SaveSignal(caught)))
        try:
            sys.settrace(interrupt_at_line)
            rmg.save_input(str(output))
        except BaseException as caught:
            error = caught
        finally:
            sys.settrace(None)
            signal.signal(signum, previous_handler)

        assert fired
        assert isinstance(error, SaveSignal)
        assert set(directory.iterdir()) <= {output}
        assert output.read_bytes() == original or b"kineticsDepositories = ['NIST']" in output.read_bytes()
        descriptor_root = "/proc/self/fd"
        if os.path.isdir(descriptor_root):
            open_paths = []
            for descriptor in os.listdir(descriptor_root):
                try:
                    open_paths.append(os.readlink(os.path.join(descriptor_root, descriptor)))
                except OSError:
                    pass
            assert not any(str(directory) in path for path in open_paths)


def test_programmatic_default_depositories_save_and_reread(monkeypatch, tmp_path):
    """A validated programmatic directive saves a deck accepted by the input reader."""
    class Database:
        def load(self, **kwargs):
            assert kwargs["kinetics_depositories"] == ["training"]

        kinetics = type("Kinetics", (), {"families": {}, "libraries": {}})()
        thermo = type("Thermo", (), {"libraries": {}, "adsorption_groups": None})()
        forbidden_structures = type("Forbidden", (), {"entries": {}})()

    import rmgpy.rmg.main as main
    monkeypatch.setattr(main, "RMGDatabase", Database)
    rmg = _minimal_programmatic_rmg(tmp_path, "default")
    rmg.load_database()
    saved = tmp_path / "saved.py"
    rmg.save_input(str(saved))

    reread = RMG(input_file=str(saved), output_directory=str(tmp_path / "reread"))
    rmg_input.read_input_file(str(saved), reread)

    assert reread.kinetics_depositories == ["training"]


def test_load_database_does_not_allow_database_to_mutate_depositories(monkeypatch, tmp_path):
    """The database loader may append training, but never to the run's own list."""
    class Database:
        def load(self, **kwargs):
            kwargs["kinetics_depositories"].append("training")
        kinetics = type("Kinetics", (), {"families": {}, "libraries": {}})()
        thermo = type("Thermo", (), {"libraries": {}})()

    import rmgpy.rmg.main as main
    monkeypatch.setattr(main, "RMGDatabase", Database)
    rmg = RMG()
    rmg.database_directory = str(tmp_path)
    rmg.reaction_libraries = []
    rmg.initial_species = []
    rmg.reaction_systems = []
    rmg.reaction_model = type("Model", (), {"vibrational_manifolds": []})()
    rmg.simulator_settings_list = [SimulatorSettings()]
    rmg.model_settings_list = [ModelSettings(tol_move_to_core=0.1)]
    rmg.kinetics_depositories = ["NIST"]
    rmg.adsorption_groups = "adsorptionPt111"
    rmg.load_database()
    assert rmg.kinetics_depositories == ["NIST"]
    output = tmp_path / "input.py"
    rmg.save_input(str(output))
    assert "kineticsDepositories = ['NIST']" in output.read_text()


def test_load_database_treats_unset_depositories_as_training_enabled(monkeypatch, tmp_path):
    """An untouched RMG object keeps the historical training-rule default."""
    class Family:
        auto_generated = False
        forward_template = type("Template", (), {"reactants": [], "products": []})()
        training_calls = 0

        def add_rules_from_training(self, thermo_database):
            self.training_calls += 1

        def fill_rules_by_averaging_up(self, verbose=False):
            pass

    family = Family()

    class Database:
        def load(self, **kwargs):
            self.load_kwargs = kwargs

        kinetics = type("Kinetics", (), {"families": {"test": family}, "libraries": {}})()
        thermo = type("Thermo", (), {"libraries": {}, "adsorption_groups": None})()
        forbidden_structures = type("Forbidden", (), {"entries": {}})()

    import rmgpy.rmg.main as main
    monkeypatch.setattr(main, "RMGDatabase", Database)
    rmg = RMG()
    rmg.database_directory = str(tmp_path)
    rmg.reaction_libraries = []
    rmg.initial_species = []
    rmg.reaction_systems = []
    rmg.reaction_model = type("Model", (), {"vibrational_manifolds": []})()
    rmg.simulator_settings_list = [SimulatorSettings()]
    rmg.model_settings_list = [ModelSettings(tol_move_to_core=0.1)]
    rmg.adsorption_groups = "adsorptionPt111"

    assert rmg.kinetics_depositories is None
    rmg.load_database()

    assert rmg.database.load_kwargs["kinetics_depositories"] is None
    assert family.training_calls == 1


def test_uncertainty_reload_does_not_allow_database_to_mutate_depositories(monkeypatch, tmp_path):
    """The uncertainty reload receives a copy of the run's depository list."""
    class Kinetics:
        families = {}

        def load_families(self, path, families, depositories):
            depositories.append("training")

    class Uncertainty:
        def __init__(self, output_directory):
            pass

        def extract_sources_from_model(self):
            pass

    import rmgpy.rmg.main as main
    monkeypatch.setattr(main, "Uncertainty", Uncertainty)
    rmg = RMG(output_directory=str(tmp_path))
    rmg.database_directory = str(tmp_path)
    rmg.database = type("Database", (), {
        "kinetics": Kinetics(),
        "thermo": object(),
    })()
    rmg.reaction_model = type("Model", (), {
        "core": type("Core", (), {"species": [], "reactions": []})(),
        "vibrational_manifolds": [],
    })()
    rmg.initial_species = []
    rmg.reaction_systems = []
    rmg.simulator_settings_list = [SimulatorSettings()]
    rmg.model_settings_list = [ModelSettings(tol_move_to_core=0.1)]
    rmg.kinetics_families = ["H_Abstraction"]
    rmg.kinetics_depositories = ["NIST"]
    rmg.adsorption_groups = "adsorptionPt111"
    rmg.uncertainty = {
        "global": False,
        "local": True,
        "uncorrelated": False,
        "correlated": False,
    }

    rmg.run_uncertainty_analysis()

    assert rmg.kinetics_depositories == ["NIST"]
    output = tmp_path / "input.py"
    rmg.save_input(str(output))
    assert "kineticsDepositories = ['NIST']" in output.read_text()


@pytest.mark.functional
def test_external_library_seed_restart_round_trip(tmp_path, monkeypatch):
    """A seed restart reopens its external library by absolute source path."""
    _external_library_restart_round_trip(tmp_path, monkeypatch, three_reactions=False)


@pytest.mark.functional
def test_external_library_three_reaction_core_and_edge_restart(tmp_path, monkeypatch):
    """Core and restart_edge share a single load of a three-reaction source."""
    _external_library_restart_round_trip(tmp_path, monkeypatch, three_reactions=True)


def _external_library_restart_round_trip(tmp_path, monkeypatch, three_reactions):
    library = tmp_path / "external" / "tiny-external"
    library.mkdir(parents=True)
    (library / "reactions.py").write_text(textwrap.dedent("""
        name = "tiny-external"
        entry(index=1, label="H2 + O2 <=> H + HO2",
              kinetics=Arrhenius(A=(1.0e6, 'm^3/(mol*s)'), n=0, Ea=(0, 'kJ/mol')))
    """))
    (library / "dictionary.txt").write_text(textwrap.dedent("""
        H2
        1 H u0 p0 c0 {2,S}
        2 H u0 p0 c0 {1,S}

        O2
        1 O u1 p2 c0 {2,S}
        2 O u1 p2 c0 {1,S}

        H
        1 H u1 p0 c0

        HO2
        1 O u1 p2 c0 {2,S}
        2 O u0 p2 c0 {1,S} {3,S}
        3 H u0 p0 c0 {2,S}
    """))
    if three_reactions:
        with (library / 'reactions.py').open('a') as f:
            f.write("entry(index=2, label='H2 + O2 <=> OH + OH', "
                    "kinetics=Arrhenius(A=(1.0e6, 'm^3/(mol*s)'), n=0, Ea=(0, 'kJ/mol')))\n")
            f.write("entry(index=3, label='H2O2 + H2 <=> H2O + H2O', "
                    "kinetics=Arrhenius(A=(1.0e6, 'm^3/(mol*s)'), n=0, Ea=(0, 'kJ/mol')))\n")
        with (library / 'dictionary.txt').open('a') as f:
            f.write(textwrap.dedent("""
                OH
                1 O u1 p2 c0 {2,S}
                2 H u0 p0 c0 {1,S}

                H2O2
                1 O u0 p2 c0 {2,S} {3,S}
                2 O u0 p2 c0 {1,S} {4,S}
                3 H u0 p0 c0 {1,S}
                4 H u0 p0 c0 {2,S}

                H2O
                1 O u0 p2 c0 {2,S} {3,S}
                2 H u0 p0 c0 {1,S}
                3 H u0 p0 c0 {1,S}
            """))
    source = os.path.realpath(str(library))
    first = tmp_path / "first"
    first.mkdir()
    input_file = first / "input.py"
    input_file.write_text(textwrap.dedent("""
        database(thermoLibraries=['primaryThermoLibrary'], reactionLibraries=[{source!r}],
                 seedMechanisms=[], kineticsDepositories=['training'],
                 kineticsFamilies=['H_Abstraction'], kineticsEstimator='rate rules')
        species(label='H2', reactive=True, structure=SMILES('[H][H]'))
        species(label='O2', reactive=True, structure=SMILES('[O][O]'))
        simpleReactor(temperature=(1000,'K'), pressure=(1,'bar'),
                      initialMoleFractions={{'H2': .67, 'O2': .33}}, terminationTime=(1e-12,'s'))
        simulator(atol=1e-16, rtol=1e-8)
        model(toleranceKeepInEdge=0, toleranceMoveToCore=1e-3, toleranceInterruptSimulation=1e-3)
        options(name='external-restart', generateOutputHTML=False, generatePlots=False,
                saveSimulationProfiles=False, saveEdgeSpecies=False)
    """.format(source=source)))
    RMG(input_file=str(input_file), output_directory=str(first)).execute()
    assert source in (first / 'seed' / 'seed' / 'reactions.py').read_text()

    restart_dir = tmp_path / "restart"
    restart_dir.mkdir()
    restart_input = first / "restart_from_seed.py"
    monkeypatch.chdir(restart_dir)
    restart = RMG(input_file=str(restart_input), output_directory=str(restart_dir))
    restart.execute()
    reactions = restart.database.kinetics.libraries['restart'].get_library_reactions()
    assert restart.database.kinetics.external_library_labels[source] == 'tiny-external'
    assert all(rxn.library != source for rxn in reactions if hasattr(rxn, 'library'))
    assert any(rxn.library == 'tiny-external' for rxn in reactions)

    seed_only_input = tmp_path / 'seed_only_restart.py'
    seed_only_input.write_text(
        restart_input.read_text()
        .replace("restartFromSeed(path='seed')", "restartFromSeed(path={!r})".format(str(first / 'seed')))
        .replace("reactionLibraries=[{!r}]".format(source), "reactionLibraries=[]")
    )
    seed_only = tmp_path / 'seed-only'
    seed_only.mkdir()
    restarted_from_seed = RMG(input_file=str(seed_only_input), output_directory=str(seed_only))
    if three_reactions:
        from rmgpy.data.kinetics.library import KineticsLibrary
        loads = []
        original_load = KineticsLibrary.load

        def tracked_load(self, filename, *args, **kwargs):
            if os.path.realpath(filename) == os.path.join(source, 'reactions.py'):
                loads.append(filename)
            return original_load(self, filename, *args, **kwargs)

        with patch.object(KineticsLibrary, 'load', tracked_load):
            restarted_from_seed.execute()
            # Startup loads restart_edge but does not currently admit it. Exercise
            # its public admission path against the actual first job's seed output.
            edge = restarted_from_seed.database.kinetics.libraries['restart_edge'].get_library_reactions()
            assert any(getattr(rxn, 'library', None) == 'tiny-external' for rxn in edge)
            restarted_from_seed.reaction_model.add_reaction_library_to_edge('restart_edge')
        assert len(loads) == 1
        assert len(restarted_from_seed.database.kinetics.libraries['tiny-external'].entries) == 3
        assert sum(label == 'tiny-external' for label, kind in
                   restarted_from_seed.database.kinetics.library_order) == 1
    else:
        restarted_from_seed.execute()
    assert restarted_from_seed.database.kinetics.external_library_labels[source] == 'tiny-external'

    shutil.move(str(library), str(tmp_path / 'moved-external'))
    missing = tmp_path / 'missing'
    missing.mkdir()
    with pytest.raises(IOError, match='External kinetics library .* recorded by restart seed is missing'):
        RMG(input_file=str(seed_only_input), output_directory=str(missing)).execute()


@pytest.mark.functional
class TestMain:
    @classmethod
    def setup_class(cls):
        """A function that is run ONCE before all unit tests in this class."""
        cls.testDir = os.path.join(originalPath, "..", "test", "rmgpy", "test_data", "mainTest")
        cls.outputDir = "output"
        cls.databaseDirectory = settings["database.directory"]

        # Read database content through symlinks, but keep export destinations
        # local so saveSeedToDatabase exercises the real writer safely.
        cls.seedExport = tempfile.TemporaryDirectory(prefix='rmg-main-seeds-')
        database_root = cls.seedExport.name
        for name in os.listdir(cls.databaseDirectory):
            if name != 'kinetics':
                os.symlink(os.path.join(cls.databaseDirectory, name), os.path.join(database_root, name))
        kinetics_root = os.path.join(database_root, 'kinetics')
        os.mkdir(kinetics_root)
        for name in os.listdir(os.path.join(cls.databaseDirectory, 'kinetics')):
            if name != 'libraries':
                os.symlink(os.path.join(cls.databaseDirectory, 'kinetics', name), os.path.join(kinetics_root, name))
        libraries_root = os.path.join(kinetics_root, 'libraries')
        os.mkdir(libraries_root)
        for name in os.listdir(os.path.join(cls.databaseDirectory, 'kinetics', 'libraries')):
            if name not in ('testSeed', 'testSeed_edge'):
                os.symlink(os.path.join(cls.databaseDirectory, 'kinetics', 'libraries', name),
                           os.path.join(libraries_root, name))
        cls.seedKinetics = os.path.join(libraries_root, 'testSeed')
        cls.seedKineticsEdge = os.path.join(libraries_root, 'testSeed_edge')
        assert not os.path.exists(cls.seedKinetics)
        assert not os.path.exists(cls.seedKineticsEdge)

        output_path = os.path.join(cls.testDir, cls.outputDir)
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        os.mkdir(output_path)

        cls.rmg = RMG(input_file=os.path.join(cls.testDir, 'input.py'), output_directory=output_path)
        with patch.dict(settings, {'database.directory': database_root}):
            cls.rmg.execute()

    @classmethod
    def teardown_class(cls):
        """A function that is run ONCE after all unit tests in this class."""
        # Reset module level database
        import rmgpy.data.rmg

        rmgpy.data.rmg.database = None

        # Remove output directory
        shutil.rmtree(os.path.join(cls.testDir, cls.outputDir))

        cls.seedExport.cleanup()

    def test_rmg_execute(self):
        """Test that RMG.execute completed successfully."""
        assert isinstance(self.rmg.database, RMGDatabase)
        assert self.rmg.done

    def test_rmg_increases_reactions(self):
        """Test that RMG.execute increases reactions and species."""
        assert len(self.rmg.reaction_model.core.reactions) > 0
        assert len(self.rmg.reaction_model.core.species) > 1
        assert len(self.rmg.reaction_model.edge.reactions) > 0
        assert len(self.rmg.reaction_model.edge.species) > 0

    def test_rmg_seed_mechanism_creation(self):
        """Test that the expected seed mechanisms are created in output directory."""
        seed_dir = os.path.join(self.testDir, self.outputDir, "seed")
        assert os.path.exists

        assert os.path.exists(os.path.join(seed_dir, "seed"))  # kinetics library folder made

        assert os.path.exists(os.path.join(seed_dir, "seed", "dictionary.txt"))  # dictionary file made
        assert os.path.exists(os.path.join(seed_dir, "seed", "reactions.py"))  # reactions file made

    def test_rmg_seed_edge_mechanism_creation(self):
        """Test that the expected seed mechanisms are created in output directory."""
        seed_dir = os.path.join(self.testDir, self.outputDir, "seed")
        assert os.path.exists

        assert os.path.exists(os.path.join(seed_dir, "seed_edge"))  # kinetics library folder made

        assert os.path.exists(os.path.join(seed_dir, "seed_edge", "dictionary.txt"))  # dictionary file made
        assert os.path.exists(os.path.join(seed_dir, "seed_edge", "reactions.py"))  # reactions file made

    def test_rmg_seed_library_creation(self):
        """Test that seed mechanisms are created in the correct database locations."""
        assert self.rmg.save_seed_to_database
        assert os.path.isfile(os.path.join(self.seedKinetics, 'reactions.py'))
        assert os.path.isfile(os.path.join(self.seedKinetics, 'dictionary.txt'))

    def test_rmg_seed_edge_library_creation(self):
        """Test that edge seed mechanisms are created in the correct database locations."""
        assert self.rmg.save_seed_to_database
        assert os.path.isfile(os.path.join(self.seedKineticsEdge, 'reactions.py'))
        assert os.path.isfile(os.path.join(self.seedKineticsEdge, 'dictionary.txt'))

    def test_rmg_rms_mechanism_files_creation(self):
        """Test that rms mechanisms are created in the correct location."""
        assert os.path.exists(os.path.join(self.rmg.output_directory,"rms"))
        assert len(os.listdir(os.path.join(self.rmg.output_directory,"rms"))) != 0
        
    def test_rmg_seed_works(self):
        """Test that the created seed libraries work.

        Note: Since this test modifies the class level RMG instance,
        it can cause other tests to fail if run out of order."""
        # Load the seed libraries into the database
        self.rmg.database.load(
            path=self.databaseDirectory,
            thermo_libraries=[],
            reaction_libraries=[self.seedKinetics, self.seedKineticsEdge],
            seed_mechanisms=[self.seedKinetics, self.seedKineticsEdge],
            kinetics_families="default",
            kinetics_depositories=[],
            depository=False,
        )

        self.rmg.reaction_model = CoreEdgeReactionModel()
        self.rmg.reaction_model.add_reaction_library_to_edge("testSeed")  # try adding seed as library
        assert len(self.rmg.reaction_model.edge.species) > 0
        assert len(self.rmg.reaction_model.edge.reactions) > 0

        self.rmg.reaction_model = CoreEdgeReactionModel()
        self.rmg.reaction_model.add_seed_mechanism_to_core("testSeed")  # try adding seed as seed mech
        assert len(self.rmg.reaction_model.core.species) > 0
        assert len(self.rmg.reaction_model.core.reactions) > 0

        self.rmg.reaction_model = CoreEdgeReactionModel()
        self.rmg.reaction_model.add_reaction_library_to_edge("testSeed_edge")  # try adding seed as library
        assert len(self.rmg.reaction_model.edge.species) > 0
        assert len(self.rmg.reaction_model.edge.reactions) > 0

        self.rmg.reaction_model = CoreEdgeReactionModel()
        self.rmg.reaction_model.add_seed_mechanism_to_core("testSeed_edge")  # try adding seed as seed mech
        assert len(self.rmg.reaction_model.core.species) > 0
        assert len(self.rmg.reaction_model.core.reactions) > 0

    def test_rmg_memory(self):
        """
        test that RMG Memory objects function properly
        """
        for rxnsys in self.rmg.reaction_systems:
            Rmem = RMG_Memory(rxnsys, None)
            Rmem.generate_cond()
            Rmem.get_cond()
            Rmem.add_t_conv_N(1.0, 0.2, 2)
            Rmem.generate_cond()
            Rmem.get_cond()

    def test_make_cantera_input_file_from_ck(self):
        """
        This test ensures that a usable Cantera input file is created via the Chemkin to Cantera conversion.
        """
        import cantera as ct

        cantera_files = os.path.join(self.rmg.output_directory, "cantera_from_ck")
        files = os.listdir(cantera_files)
        for f in files:
            if ".yaml" in f:
                try:
                    ct.Solution(os.path.join(cantera_files, f))
                except:
                    assert False, "The output Cantera file is not loadable in Cantera."
    
    def test_make_cantera_input_file_directly_1(self):
        """
        This tests to ensure that a usable Cantera input file is created via direct yaml writer 1.
        """
        import cantera as ct

        cantera_files = os.path.join(self.rmg.output_directory, "cantera1")
        files = os.listdir(cantera_files)
        for f in files:
            if ".yaml" in f:
                try:
                    ct.Solution(os.path.join(cantera_files, f))
                except:
                    assert False, "The output Cantera file is not loadable in Cantera."

    def test_make_cantera_input_file_directly_2(self):
        """
        This tests to ensure that a usable Cantera input file is created via direct yaml writer 2.
        """
        import cantera as ct

        cantera_files = os.path.join(self.rmg.output_directory, "cantera2")
        files = os.listdir(cantera_files)
        for f in files:
            if ".yaml" in f:
                try:
                    ct.Solution(os.path.join(cantera_files, f))
                except:
                    assert False, "The output Cantera file is not loadable in Cantera."

    def test_cantera_input_files_match_chemkin_later(self):
        """
        Copy the Cantera YAML files (generated directly by RMG and converted from Chemkin)
        to the test data directory so that yaml_cantera1Test can compare them.
        """
        # Copy RMG-generated YAML 1 to test data directory
        cantera_dir = os.path.join(self.rmg.output_directory, "cantera1")
        rmg_yaml_path = os.path.join(cantera_dir, 'chem_annotated.yaml')
        assert os.path.exists(rmg_yaml_path), f"RMG-generated Cantera YAML file {rmg_yaml_path} not found"
        test_data_cantera_target = os.path.join(self.testDir, '..', 'yaml_writer_data', 'cantera1', 'from_main_test.yaml')
        os.makedirs(os.path.dirname(test_data_cantera_target), exist_ok=True)
        shutil.copy(rmg_yaml_path, test_data_cantera_target)

        # Copy RMG-generated YAML 2 to test data directory
        cantera_dir = os.path.join(self.rmg.output_directory, "cantera2")
        rmg_yaml_path = os.path.join(cantera_dir, 'chem_annotated.yaml')
        assert os.path.exists(rmg_yaml_path), f"RMG-generated Cantera YAML file {rmg_yaml_path} not found"
        test_data_cantera_target = os.path.join(self.testDir, '..', 'yaml_writer_data', 'cantera2', 'from_main_test.yaml')
        os.makedirs(os.path.dirname(test_data_cantera_target), exist_ok=True)
        shutil.copy(rmg_yaml_path, test_data_cantera_target)

        # Copy chemkin-converted YAML to test data directory
        cantera_from_ck_dir = os.path.join(
            self.rmg.output_directory, "cantera_from_ck"
        )
        ck_yaml_path = os.path.join(cantera_from_ck_dir, "chem_annotated.yaml")
        assert os.path.exists(ck_yaml_path), f"Chemkin-converted YAML file {ck_yaml_path} not found"
        test_data_chemkin_target = os.path.join(self.testDir, '..', 'yaml_writer_data', 'ck2yaml', 'from_main_test.yaml')
        os.makedirs(os.path.dirname(test_data_chemkin_target), exist_ok=True)
        shutil.copy(ck_yaml_path, test_data_chemkin_target)


@pytest.mark.functional
class TestRestartWithFilters:
    @classmethod
    def setup_class(cls):
        """A function that is run ONCE before all unit tests in this class."""
        cls.testDir = os.path.join(originalPath, "..", "test", "rmgpy", "test_data", "restartTest")
        cls.outputDir = os.path.join(cls.testDir, "output_w_filters")
        cls.databaseDirectory = settings["database.directory"]

        try:
            os.mkdir(cls.outputDir)
        except FileExistsError:
            pass  # output directory already exists
        initialize_log(logging.INFO, os.path.join(cls.outputDir, "RMG.log"))

        cls.rmg = RMG(
            input_file=os.path.join(cls.testDir, "restart_w_filters.py"),
            output_directory=os.path.join(cls.outputDir),
        )

    def test_restart_with_filters(self):
        """
        Test that the RMG restart job with filters included completed without problems
        """
        self.rmg.execute()
        with open(os.path.join(self.outputDir, "RMG.log"), "r") as f:
            assert "MODEL GENERATION COMPLETED" in f.read()

    @classmethod
    def teardown_class(cls):
        """A function that is run ONCE after all unit tests in this class."""
        # Reset module level database
        import rmgpy.data.rmg

        rmgpy.data.rmg.database = None

        # Remove output directory
        shutil.rmtree(cls.outputDir)


@pytest.mark.functional
class TestRestartNoFilters:
    @classmethod
    def setup_class(cls):
        """A function that is run ONCE before all unit tests in this class."""
        cls.testDir = os.path.join(originalPath, "..", "test", "rmgpy", "test_data", "restartTest")
        cls.outputDir = os.path.join(cls.testDir, "output_no_filters")
        cls.databaseDirectory = settings["database.directory"]

        os.mkdir(cls.outputDir)
        initialize_log(logging.INFO, os.path.join(cls.outputDir, "RMG.log"))

        cls.rmg = RMG(
            input_file=os.path.join(cls.testDir, "restart_no_filters.py"),
            output_directory=os.path.join(cls.outputDir),
        )

    def test_restart_no_filters(self):
        """
        Test that the RMG restart job with no filters included completed without problems
        """
        self.rmg.execute()
        with open(os.path.join(self.outputDir, "RMG.log"), "r") as f:
            assert "MODEL GENERATION COMPLETED" in f.read()

    @classmethod
    def teardown_class(cls):
        """A function that is run ONCE after all unit tests in this class."""
        # Reset module level database
        import rmgpy.data.rmg

        rmgpy.data.rmg.database = None

        # Remove output directory
        shutil.rmtree(cls.outputDir)


@pytest.mark.functional
class TestMainFunctions:
    @classmethod
    def setup_class(cls):
        """A function that is run ONCE before all unit tests in this class."""
        cls.testDir = os.path.join(originalPath, "..", "test", "rmgpy", "test_data", "mainTest")
        cls.outputDir = os.path.join(cls.testDir, "output")
        cls.databaseDirectory = settings["database.directory"]

        os.makedirs(os.path.join(cls.testDir, cls.outputDir), exist_ok=True)

        cls.max_iter = 10

        cls.rmg = RMG(
            input_file=os.path.join(cls.testDir, "superminimal_input.py"),
            output_directory=cls.outputDir,
        )

        cls.rmg.execute(max_iterations=cls.max_iter)

    def test_save_seed_modulus(self):
        """
        Test that saveSeedModulus argument from superminimal_input.py saved the correct number of seeds
        """
        path = os.path.join(self.outputDir, "previous_seeds")
        num_dir_actual = sum(os.path.isdir(os.path.join(path, i)) for i in os.listdir(path))
        num_dir_expected = self.max_iter // 2 + 1  # +1 is for saving iteration 0
        assert num_dir_actual == num_dir_expected

    def test_max_iter(self):
        """
        Test the command line argument of -i
        """
        df = pd.read_excel(os.path.join(self.outputDir, "statistics.xls"))
        num_rows = df.shape[0]

        num_iter_actual = num_rows
        num_iter_expected = self.max_iter + 2 # +2 is for saving iteration 0, and the final after the loop ends.
        assert num_iter_actual == num_iter_expected

    @classmethod
    def teardown_class(cls):
        """A function that is run ONCE after all unit tests in this class."""
        # Reset module level database
        import rmgpy.data.rmg

        rmgpy.data.rmg.database = None

        # Remove output directory
        shutil.rmtree(cls.outputDir)


class TestProfiling:
    @classmethod
    def setup_class(cls):
        """A function that is run ONCE before all unit tests in this class."""
        # Making the profile graph requires a display. See if one is available first
        cls.display_found = False

        try:
            cls.display_found = bool(os.environ["DISPLAY"])
        except KeyError:  # This means that no display was found
            pass
        cls.test_dir = os.path.join(originalPath, "..", "test", "rmgpy", "test_data", "mainTest")

    @patch("rmgpy.rmg.main.logging")
    def test_make_profile_graph(self, mock_logging):
        """
        Test that `make_profile_graph` function behaves properly given the current display state
        """
        profile_file = os.path.join(self.test_dir, "RMG.profile")
        make_profile_graph(profile_file)
        if self.display_found:  # Check that the profile graph was made
            assert os.path.exists(os.path.join(self.test_dir, "RMG.profile.dot.pdf"))
        else:  # We can't test making a profile graph on this system, but at least test that this was recognized
            mock_logging.warning.assert_called_with(
                "Could not find a display, which is required in order to generate "
                "the profile graph. This "
                "is likely due to this job being run on a remote server without performing X11 forwarding "
                "or running the job through a job manager like SLURM.\n\n The graph can be generated later "
                "by running with the postprocessing flag `rmg.py -P input.py` from any directory/computer "
                "where both the input file and RMG.profile file are located and a display is available.\n\n"
                "Note that if the postprocessing flag is specified, this will force the graph generation "
                "regardless of if a display was found, which could cause this program to crash or freeze."
            )

    @classmethod
    def teardown_class(cls):
        """A function that is run ONCE after all unit tests in this class."""

        if cls.display_found:  # Remove output PDF
            os.remove(os.path.join(cls.test_dir, "RMG.profile.dot.pdf"))
            os.remove(os.path.join(cls.test_dir, "RMG.profile.dot"))
            os.remove(os.path.join(cls.test_dir, "RMG.profile.dot.ps2"))


class TestCanteraOutputConversion:
    """
    Tests if we can convert Chemkin files to Cantera files without crashing.
    (Or raising an exception for bad files.)
    """
    def setup_class(self):
        self.chemkin_files = {
            """ELEMENTS
	H
	D /2.014/
	T /3.016/
	C
	CI /13.003/
	O
	OI /18.000/
	N

END

SPECIES
    ethane(1)       
    CH3(4)          
END

THERM ALL
   300.000  1000.000  5000.000

ethane(1)               H 6  C 2            G100.000   5000.000  954.52        1
 4.58987205E+00 1.41507042E-02-4.75958084E-06 8.60284590E-10-6.21708569E-14    2
-1.27217823E+04-3.61762003E+00 3.78032308E+00-3.24248354E-03 5.52375224E-05    3
-6.38573917E-08 2.28633835E-11-1.16203404E+04 5.21037799E+00                   4

CH3(4)                  H 3  C 1            G100.000   5000.000  1337.62       1
 3.54144859E+00 4.76788187E-03-1.82149144E-06 3.28878182E-10-2.22546856E-14    2
 1.62239622E+04 1.66040083E+00 3.91546822E+00 1.84153688E-03 3.48743616E-06    3
-3.32749553E-09 8.49963443E-13 1.62856393E+04 3.51739246E-01                   4

END



REACTIONS    KCAL/MOLE   MOLES

CH3(4)+CH3(4)=ethane(1)                             8.260e+17 -1.400    1.000    

END
""": True,
            """ELEMENTS
	CI /13.003/
	O
	OI /18.000/
	N

END

SPECIES
    ethane(1)       
    CH3(4)          
END

THERM ALL
   300.000  1000.000  5000.000

ethane(1)               H 6  C  2            G100.000   5000.000  954.52        1
 4.58987205E+00 1.41507042E-02-4.75958084E-06 8.60284590E-10-6.21708569E-14    2
-1.27217823E+04-3.61762003E+00 3.78032308E+00-3.24248354E-03 5.52375224E-05    3
-6.38573917E-08 2.28633835E-11-1.16203404E+04 5.21037799E+00                   4

CH3(4)                  H 3  C 1            G100.000   5000.000  1337.62       1
 3.54144859E+00 4.76788187E-03-1.82149144E-06 3.28878182E-10-2.22546856E-14    2
 1.62239622E+04 1.66040083E+00 3.91546822E+00 1.84153688E-03 3.48743616E-06    3
-3.32749553E-09 8.49963443E-13 1.62856393E+04 3.51739246E-01                   4

END



REACTIONS    KCAL/MOLE   MOLES

CH3(4)+CH3(4)=ethane(1)                             8.260e+17 -1.400    1.000    

END
""": False,
            """ELEMENTS
	H
	D /2.014/
	T /3.016/
	C
	CI /13.003/
	O
	OI /18.000/
	N

END

SPECIES
    ethane(1)       
    CH3(4)          
END

THERM ALL
   300.000  1000.000  5000.000

ethane(1)               H 6  C 2            G100.000   5000.000  954.52        1
 4.58987205E+00 1.41507042E-02-4.75958084E-06 8.60284590E-10-6.21708569E-14    2
-1.27217823E+04-3.61762003E+00 3.78032308E+00-3.24248354E-03 5.52375224E-05    3
-6.38573917E-08 2.28633835E-11-1.16203404E+04 5.21037799E+00                   4

END

REACTIONS    KCAL/MOLE   MOLES

CH3(4)+CH3(4)=ethane(1)                             8.260e+17 -1.400    1.000    

END
""": False,
        }
        self.rmg = RMG()
        self.dir_name = "temp_dir_for_testing"
        self.rmg.output_directory = os.path.join(originalPath, "..", "test", "rmgpy", "test_data", self.dir_name)

        self.tran_dat = """
! Species         Shape    LJ-depth  LJ-diam   DiplMom   Polzblty  RotRelaxNum Data     
! Name            Index    epsilon/k_B sigma     mu        alpha     Zrot      Source   
ethane(1)           2     252.301     4.302     0.000     0.000     1.500    ! GRI-Mech
CH3(4)              2     144.001     3.800     0.000     0.000     0.000    ! GRI-Mech
        """

    def teardown_class(self):
        os.chdir(originalPath)
        # try to remove the tree. If testChemkinToCanteraConversion properly
        # ran, the files should already be removed.
        try:
            shutil.rmtree(self.dir_name)
        except OSError:
            pass
        # go back to the main RMG-Py directory
        os.chdir("..")

    def test_chemkin_to_cantera_conversion(self):
        """
        Tests that good and bad chemkin files raise proper exceptions
        """

        from cantera.ck2yaml import InputError

        for ck_input, works in self.chemkin_files.items():
            os.chdir(originalPath)
            os.mkdir(self.dir_name)
            os.chdir(self.dir_name)

            f = open("chem001.inp", "w")
            f.write(ck_input)
            f.close()

            f = open("tran.dat", "w")
            f.write(self.tran_dat)
            f.close()

            if works:
                self.rmg.generate_cantera_files_from_chemkin(os.path.join(os.getcwd(), "chem001.inp"))
            else:
                with pytest.raises(InputError):
                    self.rmg.generate_cantera_files_from_chemkin(os.path.join(os.getcwd(), "chem001.inp"))

            # clean up
            os.chdir(originalPath)
            shutil.rmtree(self.dir_name)
