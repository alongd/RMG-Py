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

import argparse
import logging
import os.path
import re
import shutil
import time
import tempfile
from functools import wraps


class Subject(object):
    """Subject in Observer Pattern"""

    def __init__(self):
        self._observers = []

    """
    Call this method when your (self-implemented)
    observer class should start listening to the Subject
    class.

    e.g.:

    listener = YourOwnListener()
    subject.attach(listener)
    """

    def attach(self, observer):
        if not observer in self._observers:
            self._observers.append(observer)

    """
    Call this method when your (self-implemented)
    observer class should stop listening to the Subject
    class.

    e.g.:
    listener = YourOwnListener()
    subject.attach(listener)

    ...<do some work>...

    subject.detach(listener)

    """

    def detach(self, observer):
        try:
            self._observers.remove(observer)
        except ValueError:
            pass

    """
    Call this method in classes that implement
    Subject, when the data that your interested in,
    is available.

    e.g.:
    class YourClass(Subject):
        ...
        def simulate(...)
            <stuff is being done>

            self.notify()

            <continue doing other stuff>            


    Make sure that your listener class implements the update(subject)
    method!

    e.g.:

    class YourOwnListener(object):
        def __init__(self):
            self.data = []

        def update(self, subject):
            self.data.append(subject.data)

    """

    def notify(self, modifier=None):
        for observer in self._observers:
            if modifier != observer:
                observer.update(self)


def _strip_wrapped_flow_yaml_notes(text):
    """Strip wrapped flow-style YAML notes without a nested regex."""
    lines = text.splitlines(keepends=True)
    stripped_lines = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if (
            stripped_lines
            and stripped_lines[-1].rstrip().endswith(",")
            and line.lstrip(" \t").startswith("note:")
        ):
            end_index = i
            while end_index < len(lines):
                if "}" in lines[end_index]:
                    comma_index = stripped_lines[-1].rfind(",")
                    brace_index = lines[end_index].find("}")
                    stripped_lines[-1] = stripped_lines[-1][:comma_index] + lines[end_index][brace_index:]
                    i = end_index + 1
                    break
                if lines[end_index].rstrip().endswith(","):
                    i = end_index + 1
                    break
                end_index += 1
            else:
                stripped_lines.append(line)
                i += 1
                continue
            continue

        stripped_lines.append(line)
        i += 1

    return "".join(stripped_lines)


def _find_flow_yaml_delimiter(text, start):
    """Find the next comma or closing brace outside quoted text."""
    quote = None
    i = start
    while i < len(text):
        char = text[i]
        if quote:
            if char == quote:
                if quote == "'" and i + 1 < len(text) and text[i + 1] == "'":
                    i += 2
                    continue
                quote = None
            elif quote == '"' and char == "\\":
                i += 2
                continue
        elif char in ("'", '"'):
            quote = char
        elif char in ",}":
            return i
        i += 1
    return -1


def _strip_single_line_flow_yaml_notes(text):
    """Strip single-line flow-style YAML notes without splitting quoted commas."""
    lines = text.splitlines(keepends=True)
    stripped_lines = []
    note_pattern = re.compile(r'([,{])[ \t]*note:')
    for line in lines:
        search_start = 0
        while True:
            match = note_pattern.search(line, search_start)
            if not match:
                break
            value_start = match.end()
            value_end = _find_flow_yaml_delimiter(line, value_start)
            if value_end == -1:
                break
            delimiter = line[value_end]
            if delimiter == ",":
                line = line[:match.start()] + line[value_end:]
                search_start = match.start()
            elif match.group(1) == ",":
                line = line[:match.start()] + line[value_end:]
                search_start = match.start()
            else:
                line = line[:match.start() + 1] + line[value_end:]
                search_start = match.start() + 1
        stripped_lines.append(line)

    return "".join(stripped_lines)


def make_output_subdirectory(output_directory, folder):
    """
    Create a subdirectory `folder` in the output directory. If the folder
    already exists (e.g. from a previous job) its contents are deleted.
    """
    dirname = os.path.join(output_directory, folder)
    if os.path.exists(dirname):
        # The directory already exists, so delete it (and all its content!)
        shutil.rmtree(dirname)
    os.mkdir(dirname)


def strip_yaml_notes(src, dst):
    """Read a YAML file, strip ``note:`` fields, and write the
    result to *dst*. Preserves formatting (block style, key
    ordering, etc.) - important when the source is the carefully
    crafted ck2yaml output.

    Three patterns are handled:
      1. Block-style:        ``  note: ...`` on its own line,
         possibly followed by deeper-indented continuation lines
         (multi-line literal/folded scalars).
      2. Single-line flow:   ``{..., note: foo, ...}``  -> ``{..., ...}``
      3. Wrapped flow:       a flow mapping that wraps with the
         trailing ``,`` at the end of one line and
         ``    note: foo`` on the next  -> drop the note field.
    """
    if not os.path.exists(src):
        return
    with open(src) as f:
        text = f.read()
    # State-bearing adjacency notes are identity, not removable commentary.
    from rmgpy.export import refuse_resolved_species
    if re.search(r'(?im)\b(?:electronicstate|vibrationallevel)\s+', text):
        import yaml
        from rmgpy.molecule import Molecule
        data = yaml.safe_load(text)
        references = []
        for record in data.get('species', []):
            note = record.get('note', '')
            if re.search(r'(?im)^(?:electronicstate|vibrationallevel)\s+', note):
                references.append(Molecule().from_adjacency_list(note))
        refuse_resolved_species(references, 'strip_yaml_notes ' + src + ' -> ' + dst)
    # Wrapped flow style: a flow mapping that wraps after a trailing comma,
    # with ``note: value`` on the next line.
    text = _strip_wrapped_flow_yaml_notes(text)
    text = _strip_single_line_flow_yaml_notes(text)
    # Block style: ``  note: ...\n`` plus deeper-indented
    # continuation lines.
    text = re.sub(r'^( +)note:.*\n(?:\1 +[^\n]*\n)*', '', text, flags=re.MULTILINE)
    with open(dst, "w") as f:
        f.write(text)


def timefn(fn):
    @wraps(fn)
    def measure_time(*args, **kwargs):
        t1 = time.time()
        result = fn(*args, **kwargs)
        t2 = time.time()
        logging.info("@timefn: {} took {:.2f} seconds".format(fn.__name__, t2 - t1))
        return result

    return measure_time


def parse_command_line_arguments(command_line_args=None):
    """
    Parse the command-line arguments being passed to RMG Py. This uses the
    :mod:`argparse` module, which ensures that the command-line arguments are
    sensible, parses them, and returns them.
    """

    parser = argparse.ArgumentParser(description='Reaction Mechanism Generator (RMG) is an automatic chemical reaction '
                                                 'mechanism generator that constructs kinetic models composed of '
                                                 'elementary chemical reaction steps using a general understanding of '
                                                 'how molecules react.')

    parser.add_argument('file', metavar='FILE', type=str, nargs=1,
                        help='a file describing the job to execute')

    # Options for controlling the amount of information printed to the console
    # By default a moderate level of information is printed; you can either
    # ask for less (quiet), more (verbose), or much more (debug)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-q', '--quiet', action='store_true', help='only print warnings and errors')
    group.add_argument('-v', '--verbose', action='store_true', help='print more verbose output')
    group.add_argument('-d', '--debug', action='store_true', help='print debug information')

    # Add options for controlling what directories files are written to
    parser.add_argument('-o', '--output-directory', type=str, nargs=1, default='',
                        metavar='DIR', help='use DIR as output directory')

    # Add restart option
    parser.add_argument('-r', '--restart', type=str, nargs=1, metavar='path/to/seed/', help='restart RMG from a seed',
                        default='')

    parser.add_argument('-p', '--profile', action='store_true',
                        help='run under cProfile to gather profiling statistics, and postprocess them if job completes')
    parser.add_argument('-P', '--postprocess', action='store_true',
                        help='postprocess profiling statistics from previous [failed] run; does not run the simulation')

    parser.add_argument('-t', '--walltime', type=str, nargs=1, default=None,
                        metavar='DD:HH:MM:SS', help='set the maximum execution time')

    parser.add_argument('-i', '--maxiter', type=int, nargs=1, default=None,
                        help='set the maximum number of RMG iterations')

    # Add option to select max number of processes for reaction generation
    parser.add_argument('-n', '--maxproc', type=int, nargs=1, default=1,
                        help='max number of processes used during reaction generation')

    # Add option to output a folder that stores the details of each kinetic database entry source
    parser.add_argument('-k', '--kineticsdatastore', action='store_true',
                        help='output a folder, kinetics_database, that contains a .txt file for each reaction family '
                             'listing the source(s) for each entry')

    args = parser.parse_args(command_line_args)

    # Process args to set correct default values and format

    # For output and scratch directories, if they are empty strings, set them
    # to match the input file location
    args.file = args.file[0]

    # If walltime was specified, retrieve this string from the element 1 list
    if args.walltime:
        args.walltime = args.walltime[0]

    if args.restart:
        args.restart = args.restart[0]

    if args.maxiter:
        args.maxiter = args.maxiter[0]

    if args.maxproc != 1:
        args.maxproc = args.maxproc[0]

    # Set directories
    input_directory = os.path.abspath(os.path.dirname(args.file))

    if args.output_directory == '':
        args.output_directory = input_directory
    # If output directory was specified, retrieve this string from the element 1 list
    else:
        args.output_directory = args.output_directory[0]

    if args.postprocess:
        args.profile = True

    return args


def as_list(item, default=None):
    """
    Wrap the given item in a list if it is not None and not already a list.

    Args:
        item: the item to be put in a list
        default (optional): a default value to return if the item is None
    """
    if isinstance(item, list):
        return item
    elif item is None:
        return default
    else:
        return [item]


def get_reaction_collider(side):
    """Return the complete ``(+identifier)`` group, including nested parentheses.

    State tags and indices may both use parentheses inside an identifier.
    An unclosed group is rejected rather than interpreted as a shorter name.
    """
    start = side.find('(+')
    if start < 0:
        return None
    depth = 0
    for index in range(start, len(side)):
        if side[index] == '(':
            depth += 1
        elif side[index] == ')':
            depth -= 1
            if depth == 0:
                return side[start:index + 1]
    raise ValueError('Unclosed third body collider in reaction side "{0}".'.format(side))


def write_files_atomically(entries):
    """Stage complete files, preserve modes and roll back recoverable failures.

    Individual renames are atomic. Generation stamps make a crash or concurrent
    split detectable by readers; this function does not lock other writers.
    Only this invocation's temporary files and backups are removed on failure.
    """
    entries = list(entries)
    destinations = [os.path.realpath(os.path.abspath(path)) for path, _ in entries]
    from rmgpy.exceptions import SpeciesIdentityError
    if len(set(destinations)) != len(destinations):
        raise SpeciesIdentityError('Output destinations alias the same file: {0}.'.format(destinations))
    for index, path in enumerate(destinations):
        for other in destinations[:index]:
            if os.path.exists(path) and os.path.exists(other) and os.path.samefile(path, other):
                raise SpeciesIdentityError('Output destinations alias the same file: {0}, {1}.'.format(path, other))
    staged, backups, replaced = [], {}, []
    try:
        for path, content in entries:
            path = os.path.abspath(path)
            directory = os.path.dirname(path)
            os.makedirs(directory, exist_ok=True)
            fd, temporary = tempfile.mkstemp(prefix='.' + os.path.basename(path) + '.', suffix='.tmp', dir=directory)
            # Recreate our reserved random path using the kernel's normal umask.
            # O_EXCL prevents overwriting anything created in the intervening window.
            os.close(fd)
            os.unlink(temporary)
            fd = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o666)
            staged.append((temporary, path, content))
            with os.fdopen(fd, 'w', encoding='utf-8') as stream:
                stream.write(content)
            if os.path.exists(path):
                fd, backup = tempfile.mkstemp(prefix='.' + os.path.basename(path) + '.', suffix='.bak', dir=directory)
                os.close(fd)
                backups[path] = backup
                shutil.copy2(path, backup)
                shutil.copymode(path, temporary)
        for temporary, path, content in staged:
            os.replace(temporary, path)
            replaced.append((path, content))
    except BaseException:
        for path, content in reversed(replaced):
            # Do not overwrite another writer's completed generation during rollback.
            with open(path, encoding='utf-8') as stream:
                still_ours = stream.read() == content
            if not still_ours:
                continue
            if path in backups:
                os.replace(backups[path], path)
                del backups[path]
            else:
                os.unlink(path)
        raise
    finally:
        for temporary in [item[0] for item in staged] + list(backups.values()):
            if os.path.exists(temporary):
                os.unlink(temporary)


_GENERATION_MARKER = 'RMG-PAIR-GENERATION'
_GENERATION_PATTERN = re.compile(r'^(?:#|!|//) RMG-PAIR-GENERATION (v1-sha256:[0-9a-f]{64})$')


def stamp_file_generation(entries):
    """Stamp (path, content, comment-prefix) entries with one content digest.

    The digest covers ordered UTF-8 contents with length boundaries, excluding
    paths and the stamp itself, so identical output has identical generations.
    """
    import hashlib
    entries = list(entries)
    digest = hashlib.sha256()
    for _, content, _ in entries:
        data = content.encode('utf-8')
        digest.update(len(data).to_bytes(8, 'big'))
        digest.update(data)
    token = 'v1-sha256:' + digest.hexdigest()
    return [(path, ('{0} {1} {2}\n'.format(prefix, _GENERATION_MARKER, token) if prefix else '') + content)
            for path, content, prefix in entries]


def generation_token(content):
    """Read a generation marker, preserving evidence of malformed markers."""
    lines = [line for line in content.splitlines()
             if re.match(r'^(?:#|!|//) RMG-PAIR-GENERATION(?:\s|$)', line)]
    if not lines:
        return None
    if len(lines) != 1:
        return 'invalid:multiple:' + repr(lines)
    match = _GENERATION_PATTERN.fullmatch(lines[0])
    return match.group(1) if match else 'invalid:' + lines[0]


def strip_generation_marker(content):
    """Remove only the save-generation header, retaining other legacy bytes."""
    return ''.join(line for line in content.splitlines(keepends=True)
                   if not _GENERATION_PATTERN.fullmatch(line.rstrip('\r\n')))


def read_generation_files(kinetics_paths, dictionary_path, *, strip_markers=True):
    """Read matching snapshots; retain exact marker text for external parsers if requested."""
    from rmgpy.exceptions import GenerationMismatchError
    with open(dictionary_path, encoding='utf-8') as stream:
        dictionary = stream.read()
    second = generation_token(dictionary)
    snapshots = {}
    for path in kinetics_paths:
        with open(path, encoding='utf-8') as stream:
            kinetics = stream.read()
        first = generation_token(kinetics)
        if first != second or (first is not None and first.startswith('invalid:')):
            raise GenerationMismatchError(
                'Generation mismatch: kinetics "{0}" token={1!r}; dictionary "{2}" token={3!r}.'.format(
                    path, first, dictionary_path, second))
        snapshots[path] = strip_generation_marker(kinetics) if strip_markers else kinetics
    return snapshots, strip_generation_marker(dictionary) if strip_markers else dictionary


def read_generation_pair(kinetics_path, dictionary_path):
    """Validate the exact snapshots that will be parsed, without reopening paths.

    Token-free pairs retain legacy semantics. Partial, malformed or mismatched
    stamps refuse with both paths and both tokens in GenerationMismatchError.
    """
    snapshots, dictionary = read_generation_files([kinetics_path], dictionary_path)
    return snapshots[kinetics_path], dictionary


def parse_reaction_equation(equation):
    """Split <=>, => or = without consuming a participant's first character."""
    from rmgpy.exceptions import DatabaseError
    pieces = re.split(r'(<=>|=>)', equation)
    if len(pieces) == 1:
        depth = 0
        for index, character in enumerate(equation):
            depth += (character == '(') - (character == ')')
            if character == '=' and depth == 0:
                pieces = [equation[:index], '=', equation[index + 1:]]
                break
    if len(pieces) != 3 or not pieces[0].strip() or not pieces[2].strip():
        raise DatabaseError('Invalid reaction equation: {0!r}.'.format(equation))
    return pieces[0].strip(), pieces[2].strip(), pieces[1] != '=>'
