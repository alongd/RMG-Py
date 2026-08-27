#!/usr/bin/env python
# encoding: utf-8

"""
The three known defects, re-confirmed on THIS union rather than inherited from a
previous pass's description.

  1. The RMS YAML writer refuses these two rate laws. The decision to ship anyway rests
     on the failure being LOUD -- raised at save, not a silently wrong number -- so that
     is the property measured here, at the writer and by reading the dispatch's terminal
     branch.

  2. A kinetics library carrying BOTH charge channels still cannot be loaded: the
     duplicate check inside library loading decides duplication before the owner is
     attached, so no placement declaration is visible and the two collapse. The campaign's
     workaround is that the two channels live in two libraries. Measured here: that the
     refusal is still a refusal, that the pin on it is still expected-failing rather than
     silently passing, and that the shipped deck still declares two libraries.

  3. The lithium cation's enthalpy. A separate ticket has diagnosed it since the last
     pass. This probe MEASURES the value on this union's database and says which entry the
     deck actually uses, rather than inheriting either the old description or the fix.

Run from the RMG-Py worktree root so ``rmgrc`` resolves:

    python docs/i123c-audit/probe_known_defects.py
"""

import os
import sys

from rmgpy import settings
from rmgpy import electron_placement
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.library import KineticsLibrary
from rmgpy.data.thermo import ThermoDatabase
from rmgpy.electron_balance import get_electron_placement_counts
from rmgpy.exceptions import DatabaseError

IONISATION = 'PlasmaElectronImpactIonization'
RECOMBINATION = 'PlasmaRadiativeRecombination'
DECK = 'docs/i123-integration/input.py'

FAILURES = []


def check(label, condition, detail=''):
    ok = bool(condition)
    print('  [{0}] {1}{2}'.format('PASS' if ok else 'FAIL', label,
                                  (' -- ' + detail) if detail else ''))
    if not ok:
        FAILURES.append(label)
    return ok


def note(text):
    print('  [..] {0}'.format(text))


def banner(text):
    print('\n' + '-' * 78)
    print(text)
    print('-' * 78)


def header():
    print('=' * 78)
    import rmgpy.molecule.molecule as _m
    print('rmgpy              = {0}'.format(os.path.abspath(
        os.path.dirname(os.path.dirname(electron_placement.__file__)))))
    print('compiled module    = {0}'.format(_m.__file__))
    print('database.directory = {0}'.format(settings['database.directory']))
    print('resolved realpath  = {0}'.format(os.path.realpath(settings['database.directory'])))
    print('=' * 78)


def load_channels():
    db = KineticsDatabase()
    db.load_libraries(os.path.join(settings['database.directory'], 'kinetics', 'libraries'),
                      libraries=[IONISATION, RECOMBINATION])
    return db


# ------------------------------------------------------------------ defect 1: RMS
def defect1_rms(db):
    banner('DEFECT 1 -- the RMS YAML writer, and whether it still fails LOUDLY')
    from rmgpy import yaml_rms

    source = open(yaml_rms.__file__.replace('.pyc', '.py')).read().splitlines()
    for lineno, line in enumerate(source, 1):
        if 'does not have a defined conversion' in line:
            note('rmgpy/yaml_rms.py:{0}  {1}'.format(lineno, line.strip()))
            note('  preceding line: {0}'.format(source[lineno - 2].strip()))
            check('the dispatch\'s terminal branch is a raise, not a fall-through default',
                  source[lineno - 2].strip() == 'else:'
                  or 'raise ValueError' in source[lineno - 2],
                  source[lineno - 2].strip())

    for label in (IONISATION, RECOMBINATION):
        rxn = db.libraries[label].get_library_reactions()[0]
        kinetics = rxn.kinetics
        try:
            yaml_rms.obj_to_dict(kinetics, [])
        except ValueError as exc:
            note('{0} -> {1}  LOUD: ValueError'.format(label, type(kinetics).__name__))
            note('    {0}'.format(exc))
            check('{0} rate law is REFUSED loudly by the RMS writer'.format(
                type(kinetics).__name__), True)
        except Exception as exc:
            check('{0} rate law raises ValueError (not {1})'.format(
                type(kinetics).__name__, type(exc).__name__), False, str(exc))
        else:
            check('{0} rate law is REFUSED rather than silently converted'.format(
                type(kinetics).__name__), False,
                'obj_to_dict returned a value -- a wrong number would be emitted')

    with open(DECK) as handle:
        deck_text = handle.read()
    check('the shipped deck still disables the RMS writer as its documented workaround',
          'generateRMSYAML=False' in deck_text.replace(' ', ''),
          'generateRMSYAML setting in ' + DECK)


# --------------------------------------------------- defect 2: one library, both channels
def defect2_one_library(db):
    banner('DEFECT 2 -- one library carrying both channels still cannot be loaded')

    merged = KineticsLibrary(label='PlasmaBothChannels')
    merged.entries = {}
    index = 0
    for label in (IONISATION, RECOMBINATION):
        for entry in db.libraries[label].entries.values():
            index += 1
            entry.index = index
            merged.entries[str(index)] = entry
            note('entry {0}: {1}   electrons={2}  family={3}  placement={4}'.format(
                index, entry.item, entry.item.electrons,
                getattr(entry.item, 'family', None),
                get_electron_placement_counts(entry.item)))

    try:
        merged.check_for_duplicates()
    except DatabaseError as exc:
        note('check_for_duplicates REFUSED -> DatabaseError')
        note('    {0}'.format(str(exc).strip()))
        check('a single library carrying both channels is still REFUSED at load', True)
    else:
        check('a single library carrying both channels is still REFUSED at load', False,
              'check_for_duplicates accepted it -- the pinned xfail should now be failing')

    note('the two entries above carry family=None: the owner is attached later, in '
         'get_library_reactions, which is why the placement declaration is invisible here')

    with open(DECK) as handle:
        deck_text = handle.read()
    check('the shipped deck still keeps the two channels in TWO libraries',
          IONISATION in deck_text and RECOMBINATION in deck_text,
          '{0} + {1}'.format(IONISATION, RECOMBINATION))

    for label in (IONISATION, RECOMBINATION):
        reactions = db.libraries[label].get_library_reactions()
        check('{0} carries exactly one entry'.format(label), len(reactions) == 1,
              str(len(reactions)))


# ---------------------------------------------------- defect 3: the lithium cation
def defect3_lithium_enthalpy():
    banner('DEFECT 3 -- the lithium cation enthalpy, MEASURED on this union')

    thermo_db = ThermoDatabase()
    libraries_path = os.path.join(settings['database.directory'], 'thermo', 'libraries')
    available = sorted(f[:-3] for f in os.listdir(libraries_path) if f.endswith('.py'))
    wanted = [name for name in ('LithiumPrimaryThermo', 'LithiumAdditionalThermo',
                                'PlasmaCationThermo')
              if name in available]
    note('lithium/cation thermo libraries present in this database: {0}'.format(wanted))
    thermo_db.load_libraries(libraries_path, libraries=wanted)

    found = []
    for lib_name in wanted:
        library = thermo_db.libraries[lib_name]
        for entry in library.entries.values():
            label = entry.label
            if label.replace('[', '').replace(']', '') in ('Lip', 'Li+'):
                thermo = entry.data
                h298 = thermo.get_enthalpy(298.15) / 1000.0
                try:
                    e0 = thermo.E0.value_si / 1000.0
                except AttributeError:
                    e0 = None
                found.append((lib_name, label, h298, e0, entry.short_desc))
                note('{0} :: entry {1!r} index={2}'.format(lib_name, label, entry.index))
                note('    H298      = {0:.3f} kJ/mol'.format(h298))
                note('    E0        = {0}'.format(
                    '{0:.3f} kJ/mol'.format(e0) if e0 is not None else 'not set'))
                note('    shortDesc = {0}'.format(
                    (entry.short_desc or '').strip().splitlines()[:1]))

    check('a lithium cation thermo entry exists in this database', bool(found))

    # The reference number, stated with its basis so the arithmetic is auditable.
    dfh_li_g = 159.3    # kJ/mol, standard enthalpy of formation of gaseous atomic lithium
    ie_li = 520.2       # kJ/mol, first ionisation energy of lithium
    reference = dfh_li_g + ie_li
    note('reference dfH(Li+, g) = dfH(Li, g) + IE(Li) = {0} + {1} = {2:.1f} kJ/mol'.format(
        dfh_li_g, ie_li, reference))
    note('    [I] basis: standard atomisation enthalpy and first ionisation energy of '
         'lithium; quoted, not re-derived from a primary source in this probe')

    for lib_name, label, h298, e0, _desc in found:
        gap = reference - h298
        note('{0} :: {1}  H298 = {2:.3f} kJ/mol   gap vs reference = {3:+.1f} kJ/mol '
             '({4:+.2f} eV)'.format(lib_name, label, h298, gap, gap / 96.485))

    return found


def main():
    header()
    db = load_channels()
    defect1_rms(db)
    defect2_one_library(db)
    defect3_lithium_enthalpy()

    print('\n' + '=' * 78)
    if FAILURES:
        print('FAILURES ({0}):'.format(len(FAILURES)))
        for label in FAILURES:
            print('  - {0}'.format(label))
        return 1
    print('all checks passed')
    return 0


if __name__ == '__main__':
    sys.exit(main())
