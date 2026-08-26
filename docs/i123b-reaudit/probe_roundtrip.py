#!/usr/bin/env python
# encoding: utf-8

"""
Can the mechanism the I-123 deck writes be read back?

The first audit could not ask this question: the Chemkin writer refused every plasma
rate law, so there was no file. The writer now emits one. This probe reads the deck's
own artifacts back with the same readers RMG and Cantera use, and reports what each
one does -- it asserts nothing about what *should* happen, it measures what does.

Three round trips are attempted, each independently:

  A. RMG's own Chemkin reader on the emitted `chem.inp` + `species_dictionary.txt`.
  B. The same file with the `TDEP/.../` auxiliary line removed, to isolate whether the
     TDEP line is the whole of the obstacle or only part of it.
  C. Cantera's native reader on the emitted `cantera2/chem_annotated.yaml`.

Usage, from the RMG-Py worktree root:

    python docs/i123b-reaudit/probe_roundtrip.py <deck-output-directory>
"""

import os
import sys
import traceback

from rmgpy import settings
from rmgpy.chemkin import load_chemkin_file


def banner(text):
    print('\n' + '-' * 78)
    print(text)
    print('-' * 78)


def attempt(label, fn):
    print('\n{0}'.format(label))
    try:
        result = fn()
    except Exception as exc:  # noqa: BLE001 -- reporting the class and text is the point
        print('  OUTCOME : RAISED {0}'.format(type(exc).__name__))
        for line in str(exc).strip().splitlines():
            print('            {0}'.format(line))
        print('  (traceback tail)')
        for line in traceback.format_exc().strip().splitlines()[-4:]:
            print('            {0}'.format(line))
        return None
    print('  OUTCOME : OK')
    return result


def describe(species, reactions):
    print('  species  : {0}  {1}'.format(len(species), [str(s) for s in species]))
    print('  reactions: {0}'.format(len(reactions)))
    for rxn in reactions:
        print('     {0}   kinetics={1}  electrons={2}'.format(
            rxn, type(rxn.kinetics).__name__, getattr(rxn, 'electrons', 'n/a')))


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    out = os.path.abspath(sys.argv[1])

    print('=' * 78)
    print('database.directory = {0}'.format(settings['database.directory']))
    print('deck output        = {0}'.format(out))
    print('=' * 78)

    chem = os.path.join(out, 'chemkin', 'chem.inp')
    dictionary = os.path.join(out, 'chemkin', 'species_dictionary.txt')
    transport = os.path.join(out, 'chemkin', 'tran.dat')

    banner('A. RMG Chemkin reader on the file RMG just wrote')
    with open(chem) as f:
        original = f.read()
    print('  reaction block as written:')
    for line in original.split('REACTIONS')[-1].strip().splitlines():
        if line.strip() and not line.startswith('END'):
            print('     | {0}'.format(line.rstrip()))
    a = attempt('load_chemkin_file(chem.inp)',
                lambda: load_chemkin_file(chem, dictionary_path=dictionary,
                                          transport_path=transport))
    if a:
        describe(*a[:2])

    banner('B. the same file with the TDEP auxiliary line removed')
    stripped_path = os.path.join(out, 'chemkin', 'chem_no_tdep.inp')
    kept = [ln for ln in original.splitlines(keepends=True)
            if not ln.strip().upper().startswith('TDEP/')]
    removed = len(original.splitlines()) - len(kept)
    with open(stripped_path, 'w') as f:
        f.writelines(kept)
    print('  removed {0} TDEP line(s) -> {1}'.format(removed, stripped_path))
    b = attempt('load_chemkin_file(chem_no_tdep.inp)',
                lambda: load_chemkin_file(stripped_path, dictionary_path=dictionary,
                                          transport_path=transport))
    if b:
        describe(*b[:2])

    banner('C. Cantera native reader on the yaml RMG wrote alongside it')

    def read_cantera():
        import cantera as ct
        path = os.path.join(out, 'cantera2', 'chem_annotated.yaml')
        gas = ct.Solution(path)
        return gas

    c = attempt('cantera.Solution(cantera2/chem_annotated.yaml)', read_cantera)
    if c is not None:
        print('  species  : {0}  {1}'.format(c.n_species, c.species_names))
        print('  reactions: {0}'.format(c.n_reactions))
        for i in range(c.n_reactions):
            print('     {0}'.format(c.reaction(i).equation))

    banner('SUMMARY')
    print('  A  RMG Chemkin reader, file as written      : {0}'.format('READ' if a else 'REFUSED'))
    print('  B  RMG Chemkin reader, TDEP line removed    : {0}'.format('READ' if b else 'REFUSED'))
    print('  C  Cantera native reader on RMG''s own yaml  : {0}'.format(
        'READ' if c is not None else 'REFUSED'))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
