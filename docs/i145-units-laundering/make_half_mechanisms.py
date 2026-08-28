#!/usr/bin/env python
# encoding: utf-8

"""
Split the lithium deck's Chemkin output into two half-mechanisms, one reaction each.

    python docs/i145-units-laundering/make_half_mechanisms.py <deck-output-dir>

Writes ``<deck-output-dir>/merge-half0/`` and ``<deck-output-dir>/merge-half1/``, each a
complete Chemkin model (elements, species, thermo, one reaction) plus a copy of the
deck's ``species_dictionary.txt``. These are the two inputs ``scripts/mergeModels.py``
is then run on, reproducing the shape of the failure the audit reported: two ordinary
model halves, an ordinary merge, no plasma-specific invocation and no special flags.

The split is purely textual -- the REACTIONS block is partitioned by reaction entry and
everything above it is copied verbatim -- so neither half's numbers are touched here.
"""

import os
import shutil
import sys


def split_reactions_block(text):
    """Return ``(head, header_line, [entry, ...], tail)`` for a Chemkin file."""
    lines = text.splitlines(keepends=True)
    start = next(i for i, line in enumerate(lines)
                 if line.strip().startswith('REACTIONS'))
    end = next(i for i in range(start + 1, len(lines))
               if lines[i].strip().upper() == 'END')
    head = ''.join(lines[:start])
    header_line = lines[start]
    tail = ''.join(lines[end:])

    entries, current = [], []
    for line in lines[start + 1:end]:
        stripped = line.strip()
        if not stripped:
            if current:
                current.append(line)
            continue
        # A reaction entry begins on a line carrying '=' outside a comment.
        code = stripped.split('!')[0]
        if '=' in code and not stripped.startswith('!'):
            if current:
                entries.append(''.join(current))
                current = []
        current.append(line)
    if current:
        entries.append(''.join(current))
    return head, header_line, entries, tail


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    out_dir = sys.argv[1]
    chem_path = os.path.join(out_dir, 'chemkin', 'chem.inp')
    dict_path = os.path.join(out_dir, 'chemkin', 'species_dictionary.txt')

    with open(chem_path) as handle:
        text = handle.read()
    head, header_line, entries, tail = split_reactions_block(text)
    print('reaction entries found: {0}'.format(len(entries)))
    if len(entries) != 2:
        print('expected exactly 2; got {0}'.format(len(entries)))
        for i, entry in enumerate(entries):
            print('--- entry {0} ---\n{1}'.format(i, entry))
        return 1

    for index, entry in enumerate(entries):
        half_dir = os.path.join(out_dir, 'merge-half{0}'.format(index))
        if os.path.isdir(half_dir):
            shutil.rmtree(half_dir)
        os.makedirs(half_dir)
        with open(os.path.join(half_dir, 'chem.inp'), 'w') as handle:
            handle.write(head + header_line + '\n' + entry + tail)
        shutil.copy(dict_path, os.path.join(half_dir, 'species_dictionary.txt'))
        print('wrote {0}'.format(os.path.join(half_dir, 'chem.inp')))
        print(entry.rstrip())
    return 0


if __name__ == '__main__':
    sys.exit(main())
