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

import dis
import os
import signal
import stat
import subprocess
import sys

import pytest
import yaml

import rmgpy.util as util_module
from rmgpy.util import strip_yaml_notes, write_files_atomically


class UtilTest:
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

    def test_atomic_writer_refuses_multiply_linked_target(self, tmp_path):
        output = tmp_path / "output.txt"
        alias = tmp_path / "canonical.txt"
        output.write_text("old")
        os.link(output, alias)

        with pytest.raises(OSError, match="multiple hard links"):
            write_files_atomically([(output, "new")])

        assert os.path.samefile(output, alias)
        assert output.read_text() == alias.read_text() == "old"

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

    def test_atomic_writer_preserves_acl_and_extended_attributes(self, tmp_path):
        output = tmp_path / "output.txt"
        output.write_text("old")
        os.setxattr(output, "user.rmg-provenance", b"calibration-v1")
        subprocess.run(["setfacl", "-m", f"u:{os.getuid()}:rw", str(output)], check=True)
        attributes = {name: os.getxattr(output, name) for name in os.listxattr(output)}
        access_acl = subprocess.check_output(["getfacl", "-cp", str(output)])

        write_files_atomically([(output, "new")])

        assert {name: os.getxattr(output, name) for name in os.listxattr(output)} == attributes
        assert subprocess.check_output(["getfacl", "-cp", str(output)]) == access_acl

    def test_atomic_writer_removes_inherited_default_acl_from_replacement(self, tmp_path):
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

    def test_atomic_writer_refuses_unretained_extended_attribute(self, monkeypatch, tmp_path):
        output = tmp_path / "output.txt"
        output.write_text("old")
        os.setxattr(output, "user.rmg-provenance", b"calibration-v1")
        original = output.read_bytes()
        monkeypatch.setattr(os, "setxattr", lambda path, name, value: None)

        with pytest.raises(OSError, match="retain required filesystem metadata"):
            write_files_atomically([(output, "new")])

        assert output.read_bytes() == original
        assert os.getxattr(output, "user.rmg-provenance") == b"calibration-v1"
        assert set(tmp_path.iterdir()) == {output}

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

    @pytest.mark.parametrize("target_exists", [True, False])
    @pytest.mark.parametrize("traced_function", [write_files_atomically, util_module._write_files_atomically_owned])
    def test_atomic_writer_interrupt_sweep_preserves_complete_destination(self, tmp_path, target_exists, traced_function):
        lines = sorted({line for _, line in dis.findlinestarts(traced_function.__code__)})

        for line in lines:
            directory = tmp_path / str(line)
            directory.mkdir()
            output = directory / "output.txt"
            if target_exists:
                output.write_text("old")
            fired = False
            error = None

            def interrupt_at_line(frame, event, arg):
                nonlocal fired
                if not fired and frame.f_code is traced_function.__code__ and event == "line" and frame.f_lineno == line:
                    fired = True
                    os.kill(os.getpid(), signal.SIGINT)
                return interrupt_at_line

            try:
                sys.settrace(interrupt_at_line)
                write_files_atomically([(output, "new")])
            except BaseException as caught:
                error = caught
            finally:
                sys.settrace(None)

            if fired:
                assert isinstance(error, KeyboardInterrupt)
            assert set(directory.iterdir()) <= {output}
            if output.exists():
                assert output.read_text() in {"old", "new"}

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
