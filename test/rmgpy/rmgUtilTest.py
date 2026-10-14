#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
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
import collections
import errno
import fcntl
import inspect
import json
import os
import resource
import shutil
import signal
import stat
import subprocess
import sys
import textwrap

import pytest
import yaml

import rmgpy.util as util_module
from rmgpy.exceptions import SpeciesIdentityError
from rmgpy.rmg import input as rmg_input
from rmgpy.util import strip_yaml_notes, write_files_atomically, write_file_with_backup


class UtilTest:
    @pytest.mark.parametrize("mode", ["term-only", "ignore-int-then-term", "return-both"])
    def test_atomic_writer_delivers_every_held_signal(self, tmp_path, mode):
        target = tmp_path / "output"
        target.write_bytes(b"OLD")
        child = textwrap.dedent('''
            import json
            import os
            import signal
            import subprocess
            import sys
            from rmgpy import util

            mode, path = sys.argv[1:]
            delivered = []
            fired = False
            if mode == "return-both":
                def handler(signum, frame):
                    delivered.append(signum)
                signal.signal(signal.SIGINT, handler)
                signal.signal(signal.SIGTERM, handler)
            else:
                signal.signal(signal.SIGINT, signal.SIG_IGN)
                signal.signal(signal.SIGTERM, signal.SIG_DFL)
            expected = {s: signal.getsignal(s) for s in (signal.SIGINT, signal.SIGTERM)}

            def trace(frame, event, arg):
                global fired
                if frame.f_code is util._rewrite_save.__code__ and event == "call" and not fired:
                    fired = True
                    sender = "import os, signal, time; "
                    if mode != "term-only":
                        sender += "os.kill({}, signal.SIGINT); time.sleep(0.05); ".format(os.getpid())
                    sender += "os.kill({}, signal.SIGTERM)".format(os.getpid())
                    subprocess.run([sys.executable, "-c", sender], check=True, timeout=10)
                if frame.f_code is util._DeferredFileSignals.__exit__.__code__ and event == "call":
                    print("PENDING", json.dumps([s for s, f in frame.f_locals["self"].pending]), flush=True)
                return trace

            try:
                sys.settrace(trace)
                util.write_files_atomically([(path, "NEW")])
            finally:
                sys.settrace(None)
            print("RETURNED", json.dumps({
                "delivered": delivered,
                "handlers_restored": all(signal.getsignal(s) == h for s, h in expected.items()),
            }), flush=True)
        ''')
        result = subprocess.run(
            [sys.executable, "-c", child, mode, str(target)],
            capture_output=True, text=True, timeout=30,
        )
        assert result.stderr == ""
        pending = [signal.SIGTERM] if mode == "term-only" else [signal.SIGINT, signal.SIGTERM]
        assert result.stdout.splitlines()[0] == "PENDING " + json.dumps(pending)
        assert target.read_bytes() == b"NEW"
        assert not list(tmp_path.glob(".rmg-*"))
        if mode == "return-both":
            assert result.returncode == 0
            returned = json.loads(result.stdout.splitlines()[1].removeprefix("RETURNED "))
            assert returned == {"delivered": pending, "handlers_restored": True}
        else:
            assert result.returncode == -signal.SIGTERM

    @pytest.mark.parametrize("link_kind", ["relative", "absolute", "dangling"])
    def test_atomic_writer_preserves_symlink_target_semantics(self, tmp_path, link_kind):
        target = tmp_path / "target.txt"
        link = tmp_path / "output.txt"
        if link_kind != "dangling":
            target.write_text("old")
        os.symlink(str(target) if link_kind == "absolute" else target.name, link)

        write_files_atomically([(link, "new")])

        assert link.is_symlink()
        assert target.read_text() == "new"

    def test_atomic_writer_preserves_mode_and_refuses_read_only_target(self, tmp_path):
        output = tmp_path / "output.txt"
        output.write_text("old")
        os.chmod(output, 0o664)
        owner = (output.stat().st_uid, output.stat().st_gid)

        write_files_atomically([(output, "new")])
        assert output.stat().st_mode & 0o777 == 0o664
        assert (output.stat().st_uid, output.stat().st_gid) == owner

        os.chmod(output, 0o444)
        with pytest.raises(PermissionError, match="not writable"):
            write_files_atomically([(output, "refused")])
        assert output.read_text() == "new"
        assert output.stat().st_mode & 0o777 == 0o444

    def test_atomic_writer_updates_every_hard_link(self, tmp_path):
        output = tmp_path / "output.txt"
        alias = tmp_path / "canonical.txt"
        output.write_text("old")
        os.link(output, alias)

        write_files_atomically([(output, "new")])

        assert os.path.samefile(output, alias)
        assert output.read_text() == alias.read_text() == "new"

    def test_atomic_writer_follows_symlink_parent_before_dotdot(self, tmp_path):
        outside = tmp_path / "outside"
        child = outside / "child"
        child.mkdir(parents=True)
        canonical = outside / "output.txt"
        canonical.write_text("old canonical")
        sibling = tmp_path / "output.txt"
        sibling.write_text("old sibling")
        link = tmp_path / "link"
        link.symlink_to(outside.name + "/child")

        write_files_atomically([(link / ".." / "output.txt", "new")])

        assert canonical.read_text() == "new"
        assert sibling.read_text() == "old sibling"

    @pytest.mark.parametrize("flag", ["d", "A", "S"])
    def test_atomic_writer_preserves_inode_flags_and_atime(self, tmp_path, flag):
        if shutil.which("chattr") is None:
            pytest.skip("chattr is unavailable on this platform")
        output = tmp_path / "output.txt"
        output.write_text("old")
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
            write_files_atomically([(output, "new")])
            assert inode_flags() == before_flags
            assert output.stat().st_atime_ns == before_atime
        finally:
            subprocess.run(["chattr", "-" + flag, str(output)], check=True)

    def test_atomic_save_production_has_no_extended_attribute_calls(self):
        source = inspect.getsource(util_module)

        assert "listxattr" not in source
        assert "getxattr" not in source
        assert "setxattr" not in source
        assert "removexattr" not in source

    @pytest.mark.parametrize("kind", ["directory", "fifo", "self-loop", "two-loop"])
    def test_atomic_writer_refuses_nonregular_and_unresolved_destinations(self, tmp_path, kind):
        output = tmp_path / "output.txt"
        if kind == "directory":
            output.mkdir()
        elif kind == "fifo":
            os.mkfifo(output)
        elif kind == "self-loop":
            output.symlink_to(output.name)
        else:
            other = tmp_path / "other.txt"
            output.symlink_to(other.name)
            other.symlink_to(output.name)
        original_entries = {path.name: path.lstat().st_mode for path in tmp_path.iterdir()}

        with pytest.raises(OSError):
            write_files_atomically([(output, "new")])

        assert {path.name: path.lstat().st_mode for path in tmp_path.iterdir()} == original_entries
        if kind == "fifo":
            assert stat.S_ISFIFO(output.lstat().st_mode)
        if "loop" in kind:
            assert output.is_symlink()

    @pytest.mark.parametrize("kind", ["self-loop", "two-loop"])
    def test_atomic_writer_refuses_loop_hidden_behind_dotdot_path(self, tmp_path, kind):
        output = tmp_path / "output.txt"
        if kind == "self-loop":
            output.symlink_to(output.name)
        else:
            other = tmp_path / "other.txt"
            output.symlink_to(other.name)
            other.symlink_to(output.name)
        alias = tmp_path / "alias.txt"
        alias.symlink_to("missing/../output.txt")
        before = {path.name: os.readlink(path) for path in tmp_path.iterdir()}

        with pytest.raises(OSError):
            write_files_atomically([(alias, "new")])

        assert {path.name: os.readlink(path) for path in tmp_path.iterdir()} == before

    def test_atomic_writer_preserves_acl_and_extended_attributes(self, tmp_path):
        if not hasattr(os, "setxattr") or not shutil.which("setfacl") or not shutil.which("getfacl"):
            pytest.skip("filesystem xattrs or ACL command-line tools are unavailable")
        output = tmp_path / "output.txt"
        output.write_text("old")
        os.setxattr(output, "user.rmg-provenance", b"calibration-v1")
        subprocess.run(["setfacl", "-m", f"u:{os.getuid()}:rw", str(output)], check=True)
        attributes = {name: os.getxattr(output, name) for name in os.listxattr(output)}
        access_acl = subprocess.check_output(["getfacl", "-cp", str(output)])

        write_files_atomically([(output, "new")])

        assert {name: os.getxattr(output, name) for name in os.listxattr(output)} == attributes
        assert subprocess.check_output(["getfacl", "-cp", str(output)]) == access_acl

    def test_atomic_writer_does_not_inherit_directory_default_acl(self, tmp_path):
        if not hasattr(os, "setxattr") or not shutil.which("setfacl") or not shutil.which("getfacl"):
            pytest.skip("filesystem xattrs or ACL command-line tools are unavailable")
        subprocess.run(["setfacl", "-m", f"d:u:{os.getuid()}:rwx", str(tmp_path)], check=True)
        output = tmp_path / "output.txt"
        output.write_text("old")
        subprocess.run(["setfacl", "-b", str(output)], check=True)
        os.chmod(output, 0o660)
        attributes = {name: os.getxattr(output, name) for name in os.listxattr(output)}
        access_acl = subprocess.check_output(["getfacl", "-cp", str(output)])

        write_files_atomically([(output, "new")])

        assert {name: os.getxattr(output, name) for name in os.listxattr(output)} == attributes
        assert subprocess.check_output(["getfacl", "-cp", str(output)]) == access_acl

    @pytest.mark.parametrize("file_umask", [0o022, 0o002])
    def test_atomic_writer_new_file_obeys_umask(self, tmp_path, file_umask):
        output = tmp_path / "output.txt"
        previous_umask = os.umask(file_umask)
        try:
            write_files_atomically([(output, "new")])
        finally:
            os.umask(previous_umask)
        assert output.stat().st_mode & 0o777 == 0o666 & ~file_umask

    def test_atomic_writer_accepts_maximum_length_basename(self, tmp_path):
        output = tmp_path / ("x" * 255)

        write_files_atomically([(output, "new")])

        assert output.read_text() == "new"

    def test_atomic_writer_refuses_trailing_separator(self, tmp_path):
        output = tmp_path / "output.txt"
        output.write_text("old")

        with pytest.raises(OSError):
            write_files_atomically([(str(output) + os.sep, "new")])

        assert output.read_text() == "old"

    @pytest.mark.parametrize("route", ["single", "multi"])
    def test_atomic_writer_removes_failed_dangling_symlink_target(self, tmp_path, route):
        link = tmp_path / "output.txt"
        target = tmp_path / "target.txt"
        link.symlink_to(target.name)
        old_limit = resource.getrlimit(resource.RLIMIT_FSIZE)
        old_handler = signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
        try:
            resource.setrlimit(resource.RLIMIT_FSIZE, (4, old_limit[1]))
            with pytest.raises(OSError) as raised:
                if route == "single":
                    write_file_with_backup(link, b"0123456789")
                else:
                    write_files_atomically([(link, "0123456789")])
            assert raised.value.errno == errno.EFBIG
        finally:
            resource.setrlimit(resource.RLIMIT_FSIZE, old_limit)
            signal.signal(signal.SIGXFSZ, old_handler)

        assert link.is_symlink()
        assert not target.exists()
        assert list(tmp_path.glob(".rmg-*")) == []

    @pytest.mark.parametrize("alias_kind", ["lexical", "dangling-links"])
    def test_atomic_writer_refuses_aliases_of_absent_destination(self, tmp_path, alias_kind):
        target = tmp_path / "target.txt"
        if alias_kind == "lexical":
            first = tmp_path / "target.txt"
            second = os.fspath(tmp_path) + "/./target.txt"
            links = []
        else:
            first = tmp_path / "first.txt"
            second = tmp_path / "second.txt"
            first.symlink_to(target.name)
            second.symlink_to(target.name)
            links = [first, second]

        with pytest.raises(SpeciesIdentityError, match="alias the same file"):
            write_files_atomically([(first, "FIRST"), (second, "SECOND")])

        assert not target.exists()
        assert all(link.is_symlink() for link in links)
        assert list(tmp_path.glob(".rmg-*")) == []

    @pytest.mark.parametrize("outcome", ["success", "efbig"])
    @pytest.mark.parametrize("path_kind", ["direct", "dangling", "chain"])
    def test_atomic_writer_refuses_destination_created_after_preparation(
        self, tmp_path, path_kind, outcome
    ):
        target = tmp_path / "target.txt"
        destination = target
        links = []
        if path_kind != "direct":
            destination = tmp_path / "output.txt"
            destination.symlink_to("target.txt" if path_kind == "dangling" else "middle.txt")
            links.append(destination)
        if path_kind == "chain":
            middle = tmp_path / "middle.txt"
            middle.symlink_to(target.name)
            links.append(middle)

        creator = (
            "from pathlib import Path; import sys; "
            "stream = Path(sys.argv[1]).open('xb'); "
            "stream.write(b'CREATED BY OTHER TOOL'); stream.close()"
        )
        old_limit = resource.getrlimit(resource.RLIMIT_FSIZE)
        old_handler = signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
        creator_exit = None
        fired = False

        def create_between_preparation_and_open(frame, event, arg):
            nonlocal creator_exit, fired
            if frame.f_code is util_module._rewrite_save.__code__ and event == "call" and not fired:
                fired = True
                resource.setrlimit(resource.RLIMIT_FSIZE, old_limit)
                creator_exit = subprocess.run(
                    [sys.executable, "-c", creator, os.fspath(target)], check=True
                ).returncode
                if outcome == "efbig":
                    resource.setrlimit(resource.RLIMIT_FSIZE, (4, old_limit[1]))
            return create_between_preparation_and_open

        try:
            if outcome == "efbig":
                resource.setrlimit(resource.RLIMIT_FSIZE, (4, old_limit[1]))
            sys.settrace(create_between_preparation_and_open)
            with pytest.raises(FileExistsError) as raised:
                write_files_atomically([(destination, "0123456789")])
        finally:
            sys.settrace(None)
            resource.setrlimit(resource.RLIMIT_FSIZE, old_limit)
            signal.signal(signal.SIGXFSZ, old_handler)

        assert raised.value.errno == errno.EEXIST
        assert os.fspath(destination) in str(raised.value)
        assert fired and creator_exit == 0
        assert target.read_bytes() == b"CREATED BY OTHER TOOL"
        assert all(link.is_symlink() for link in links)
        assert list(tmp_path.glob(".rmg-*")) == []

    @pytest.mark.parametrize("signum", [signal.SIGINT, signal.SIGTERM])
    def test_atomic_writer_signal_sweep_preserves_complete_destination(self, tmp_path, signum):
        save_names = {
            "write_files_atomically",
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

        seed = tmp_path / "seed.txt"
        seed.write_text("old")
        try:
            sys.settrace(collect)
            write_files_atomically([(seed, "new")])
        finally:
            sys.settrace(None)

        class SaveSignal(BaseException):
            pass

        for code, line in sorted(reached, key=lambda item: (item[0].co_name, item[1])):
            directory = tmp_path / (code.co_name + "-" + str(line))
            directory.mkdir()
            output = directory / "output.txt"
            output.write_text("old")
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
                write_files_atomically([(output, "new")])
            except BaseException as caught:
                error = caught
            finally:
                sys.settrace(None)
                signal.signal(signum, previous_handler)

            assert fired
            assert isinstance(error, SaveSignal)
            assert set(directory.iterdir()) <= {output}
            assert output.read_text() in {"old", "new"}
            descriptor_root = "/proc/self/fd"
            if os.path.isdir(descriptor_root):
                open_paths = []
                for descriptor in os.listdir(descriptor_root):
                    try:
                        open_paths.append(os.readlink(os.path.join(descriptor_root, descriptor)))
                    except OSError:
                        pass
                assert not any(str(directory) in path for path in open_paths)

    def test_save_context_signal_sweep_restores_handlers_and_filesystem_errors(self, tmp_path):
        functions = [
            write_files_atomically,
            write_file_with_backup,
            rmg_input._write_input_file_atomically,
            util_module._write_files_with_backups,
            util_module._open_backup,
            util_module._prepare_save,
            util_module._rewrite_save,
            util_module._restore_save,
            util_module._remove_backup,
            util_module._DeferredFileSignals.__init__,
            util_module._DeferredFileSignals.__enter__,
            util_module._DeferredFileSignals.__exit__,
            util_module._DeferredFileSignals._restore_handlers,
            util_module._DeferredFileSignals._deliver_pending,
        ]
        codes = set()

        def add_code(code):
            codes.add(code)
            for constant in code.co_consts:
                if isinstance(constant, type(code)):
                    add_code(constant)

        for function in functions:
            add_code(function.__code__)

        class SaveSignal(BaseException):
            pass

        original_handlers = {
            signum: signal.getsignal(signum)
            for signum in (signal.SIGINT, signal.SIGTERM, signal.SIGXFSZ)
        }
        old_limit = resource.getrlimit(resource.RLIMIT_FSIZE)
        serial = 0
        totals = collections.Counter()
        divergences = []

        def has_efbig(error):
            seen = set()
            while error is not None and id(error) not in seen:
                seen.add(id(error))
                if isinstance(error, OSError) and error.errno == errno.EFBIG:
                    return True
                error = error.__cause__ or error.__context__
            return False

        try:
            for route in ("input", "util"):
                for scenario in ("success", "restore", "absent-failure", "multi-restore"):
                    if route == "input" and scenario == "multi-restore":
                        continue
                    for signum in (signal.SIGINT, signal.SIGTERM):
                        def prepare():
                            nonlocal serial
                            serial += 1
                            directory = tmp_path / str(serial)
                            directory.mkdir()
                            output = directory / "output.txt"
                            old = b"OLD" + b"q" * 1021
                            if scenario != "absent-failure":
                                output.write_bytes(old)
                            if scenario == "multi-restore":
                                (directory / "first.txt").write_bytes(b"first old")
                            return directory, output, old

                        def call(output):
                            if scenario != "success":
                                resource.setrlimit(resource.RLIMIT_FSIZE, (4096, old_limit[1]))
                            if route == "input":
                                rmg_input._write_input_file_atomically(os.fspath(output), b"N" * 20000)
                            elif scenario == "multi-restore":
                                write_files_atomically([
                                    (output.parent / "first.txt", "first new"),
                                    (output, "N" * 20000),
                                ])
                            else:
                                write_files_atomically([(output, "N" * 20000)])

                        _, output, _ = prepare()
                        events = []
                        counts = collections.Counter()

                        def discover(frame, event, arg):
                            if frame.f_code in codes and event in ("line", "return"):
                                key = (frame.f_code, event, frame.f_lineno)
                                counts[key] += 1
                                events.append((key, counts[key]))
                            return discover

                        signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
                        try:
                            sys.settrace(discover)
                            call(output)
                        except BaseException:
                            pass
                        finally:
                            sys.settrace(None)
                            resource.setrlimit(resource.RLIMIT_FSIZE, old_limit)

                        for target, occurrence in events:
                            directory, output, old = prepare()
                            seen = 0
                            fired = False
                            error = None

                            def interrupt(frame, event, arg):
                                nonlocal seen, fired
                                if (frame.f_code, event, frame.f_lineno) == target:
                                    seen += 1
                                    if seen == occurrence:
                                        fired = True
                                        os.kill(os.getpid(), signum)
                                return interrupt

                            def handler(caught, frame):
                                raise SaveSignal(caught)

                            expected_handlers = {signal.SIGINT: handler, signal.SIGTERM: handler}
                            for caught, expected in expected_handlers.items():
                                signal.signal(caught, expected)
                            try:
                                sys.settrace(interrupt)
                                call(output)
                            except BaseException as caught:
                                error = caught
                            finally:
                                sys.settrace(None)
                                resource.setrlimit(resource.RLIMIT_FSIZE, old_limit)
                                leaked_handlers = [
                                    caught
                                    for caught, expected in expected_handlers.items()
                                    if signal.getsignal(caught) is not expected
                                ]
                                for caught, expected in expected_handlers.items():
                                    signal.signal(caught, expected)

                            data = output.read_bytes() if output.exists() else None
                            backups = list(directory.glob(".rmg-*"))
                            open_files = []
                            if os.path.isdir("/proc/self/fd"):
                                for descriptor in os.listdir("/proc/self/fd"):
                                    try:
                                        destination = os.readlink(os.path.join("/proc/self/fd", descriptor))
                                    except OSError:
                                        continue
                                    if os.fspath(directory) + os.sep in destination:
                                        open_files.append(destination)
                            first_complete = (
                                scenario != "multi-restore"
                                or (directory / "first.txt").read_bytes() in (b"first old", b"first new")
                            )
                            filesystem_error_replaced = (
                                has_efbig(error)
                                and (
                                    target[0].co_name in {"__exit__", "_deliver_pending"}
                                    or (
                                        target[0].co_name == "_restore_handlers"
                                        and occurrence == 1
                                    )
                                )
                                and not (isinstance(error, OSError) and error.errno == errno.EFBIG)
                            )
                            if isinstance(error, OSError) and error.errno == errno.EFBIG:
                                totals["original_efbig"] += 1
                            totals["targets"] += 1
                            totals["fired"] += int(fired)
                            if (
                                not fired
                                or leaked_handlers
                                or backups
                                or open_files
                                or data not in (old, b"N" * 20000, None)
                                or not first_complete
                                or filesystem_error_replaced
                            ):
                                divergences.append({
                                    "route": route,
                                    "scenario": scenario,
                                    "signal": int(signum),
                                    "function": target[0].co_name,
                                    "event": target[1],
                                    "line": target[2],
                                    "occurrence": occurrence,
                                    "error": repr(error),
                                    "leaked_handlers": [int(item) for item in leaked_handlers],
                                    "backups": [os.fspath(item) for item in backups],
                                    "open_files": open_files,
                                    "length": None if data is None else len(data),
                                    "filesystem_error_replaced": filesystem_error_replaced,
                                })
        finally:
            sys.settrace(None)
            resource.setrlimit(resource.RLIMIT_FSIZE, old_limit)
            for signum, handler in original_handlers.items():
                signal.signal(signum, handler)

        assert totals["targets"] == totals["fired"]
        assert totals["original_efbig"] > 0
        print("SIGNAL_R11_SUMMARY", dict(totals), "divergences", len(divergences))
        assert not divergences, json.dumps(divergences[:20], indent=2)

    def run_strip(self, tmp_path, source_text):
        source_path = tmp_path / "chem_annotated.yaml"
        destination_path = tmp_path / "chem.yaml"
        yaml.safe_load(source_text)  # just to check it's valid YAML
        source_path.write_text(source_text)
        yaml.safe_load(source_path.read_text())
        strip_yaml_notes(source_path, destination_path)

        result = destination_path.read_text()
        yaml.safe_load(result)  # just to check it's valid YAML
        return result

    def test_strip_yaml_notes_removes_block_style_note(self, tmp_path):
        source = """species:
- name: CH4
  composition: {C: 1, H: 4}
  note: RMG-generated species note
- name: O2
  composition: {O: 2}
"""
        expected = """species:
- name: CH4
  composition: {C: 1, H: 4}
- name: O2
  composition: {O: 2}
"""

        assert self.run_strip(tmp_path, source) == expected

    def test_strip_yaml_notes_removes_block_style_multiline_note(self, tmp_path):
        source = """reactions:
- equation: CH4 + O2 <=> CH3 + HO2
  rate-constant: {A: 1.0e+06, b: 0.0, Ea: 10000.0}
  note: |
    Estimated by RMG.
    Includes database comments.
- equation: H + O2 <=> O + OH
  rate-constant: {A: 1.0e+14, b: 0.0, Ea: 15000.0}
"""
        expected = """reactions:
- equation: CH4 + O2 <=> CH3 + HO2
  rate-constant: {A: 1.0e+06, b: 0.0, Ea: 10000.0}
- equation: H + O2 <=> O + OH
  rate-constant: {A: 1.0e+14, b: 0.0, Ea: 15000.0}
"""

        assert self.run_strip(tmp_path, source) == expected

    def test_strip_yaml_notes_removes_single_line_flow_note(self, tmp_path):
        source = """species:
- name: Ar
  transport: {model: gas, note: RMG transport, geometry: atom, diameter: 3.33, well-depth: 136.5}
"""
        expected = """species:
- name: Ar
  transport: {model: gas, geometry: atom, diameter: 3.33, well-depth: 136.5}
"""

        assert self.run_strip(tmp_path, source) == expected

    def test_strip_yaml_notes_removes_single_line_flow_note_with_quoted_comma(self, tmp_path):
        source = """species:
- name: Ar
  transport: {model: gas, note: 'comma here, danger!', geometry: atom, diameter: 3.33, well-depth: 136.5}
"""
        expected = """species:
- name: Ar
  transport: {model: gas, geometry: atom, diameter: 3.33, well-depth: 136.5}
"""

        assert self.run_strip(tmp_path, source) == expected

    def test_strip_yaml_notes_removes_wrapped_flow_note(self, tmp_path):
        source = """species:
- name: CH4
  transport: {model: gas, geometry: nonlinear, diameter: 3.746,
    note: RMG transport note,
    well-depth: 141.4}
"""
        expected = """species:
- name: CH4
  transport: {model: gas, geometry: nonlinear, diameter: 3.746,
    well-depth: 141.4}
"""

        assert self.run_strip(tmp_path, source) == expected

    def test_strip_yaml_notes_removes_wrapped_flow_multiline_note(self, tmp_path):
        source = """species:
- name: CH4
  transport: {model: gas, geometry: nonlinear, diameter: 3.746,
    note: RMG transport note
      with wrapped detail,
    well-depth: 141.4}
"""
        expected = """species:
- name: CH4
  transport: {model: gas, geometry: nonlinear, diameter: 3.746,
    well-depth: 141.4}
"""

        assert self.run_strip(tmp_path, source) == expected
