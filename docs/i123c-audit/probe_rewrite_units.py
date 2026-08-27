#!/usr/bin/env python
# encoding: utf-8

"""
What happens when a plasma mechanism that has been READ BACK is written out again.

The audit asked what the newly-working Chemkin read-back means for a reloaded mechanism.
The representation answer is in ``probe_readback_identity.py``. This probe answers the
other half: whether the reloaded object survives a second trip through the writer.

It does not. The rate constant comes back 10**6 too small, silently, and the file that
carries it declares the units it is not in.

The mechanism, stated so it can be checked rather than believed:

  * ``rmgpy/chemkin.pyx::_plasma_arrhenius_for_chemkin`` reduces a
    :class:`TwoTemperaturePlasma` to the modified-Arrhenius form Chemkin can write, and
    builds that Arrhenius with ``A=(kinetics.A.value_si, kinetics.A.units)`` -- the value
    in SI, paired with whatever units string the object happens to declare.
  * Every ``TwoTemperaturePlasma`` RMG builds for itself declares SI units:
    ``to_two_temp_plasma`` on each of the plasma rate laws passes
    ``A=(arr.A.value_si, "m^3/(mol*s)")``. For those the pairing is accidentally right and
    the defect is invisible.
  * The Chemkin READER is the one producer of a ``TwoTemperaturePlasma`` whose units are
    NOT SI: it rebuilds the rate law from the file's own ``cm^3/(mol*s)`` numbers. Pairing
    an SI value with ``cm^3/(mol*s)`` understates a second-order rate by exactly 10**6.

So the writer defect is old and the reader is new, and it is the reader that makes the
defect reachable. This probe measures both halves separately so the attribution is not an
inference.

Run from the repo root so ``rmgrc`` resolves:

    python docs/i123c-audit/probe_rewrite_units.py <deck-output-dir>
"""

import os
import sys

import rmgpy.kinetics as _kinetics
from rmgpy import settings
from rmgpy import electron_placement
from rmgpy.chemkin import (load_chemkin_file, save_chemkin_file,
                           _plasma_arrhenius_for_chemkin)

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


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    out_dir = sys.argv[1]
    chem_path = os.path.join(out_dir, 'chemkin', 'chem.inp')
    dict_path = os.path.join(out_dir, 'chemkin', 'species_dictionary.txt')

    print('=' * 78)
    import rmgpy.molecule.molecule as _m
    print('rmgpy              = {0}'.format(os.path.abspath(
        os.path.dirname(os.path.dirname(electron_placement.__file__)))))
    print('compiled module    = {0}'.format(_m.__file__))
    print('database.directory = {0}'.format(settings['database.directory']))
    print('resolved realpath  = {0}'.format(os.path.realpath(settings['database.directory'])))
    print('chemkin file       = {0}'.format(chem_path))
    print('=' * 78)

    # ------------------------------------------------------------------ the unit pairing
    banner('the writer\'s unit pairing, driven directly on two TwoTemperaturePlasma objects')
    si = _kinetics.TwoTemperaturePlasma(
        A=(1.49e12, 'm^3/(mol*s)'), n=-0.267,
        Ea_g=(6.784e5, 'J/mol'), Ea_e=(6.784e5, 'J/mol'), T0=(1.0, 'K'))
    cgs = _kinetics.TwoTemperaturePlasma(
        A=(1.49e18, 'cm^3/(mol*s)'), n=-0.267,
        Ea_g=(6.784e5, 'J/mol'), Ea_e=(6.784e5, 'J/mol'), T0=(1.0, 'K'))
    note('the same rate law, declared two ways:')
    note('    SI  : A = {0}   value_si = {1:.6e}'.format(si.A, si.A.value_si))
    note('    cgs : A = {0}   value_si = {1:.6e}'.format(cgs.A, cgs.A.value_si))
    check('the two objects ARE the same rate law',
          abs(si.A.value_si - cgs.A.value_si) / si.A.value_si < 1e-9,
          '{0:.6e} vs {1:.6e} in SI'.format(si.A.value_si, cgs.A.value_si))

    arr_si, _note_si = _plasma_arrhenius_for_chemkin(si)
    arr_cgs, _note_cgs = _plasma_arrhenius_for_chemkin(cgs)
    note('reduced for Chemkin:')
    note('    from the SI-declared object  : A = {0}   value_si = {1:.6e}'.format(
        arr_si.A, arr_si.A.value_si))
    note('    from the cgs-declared object : A = {0}   value_si = {1:.6e}'.format(
        arr_cgs.A, arr_cgs.A.value_si))
    ratio = arr_cgs.A.value_si / arr_si.A.value_si
    note('    ratio of the two reductions  : {0:.6e}'.format(ratio))
    check('the reduction is UNIT-DEPENDENT -- the same rate law reduces to two different '
          'numbers depending on the units string it declares',
          abs(ratio - 1e-6) / 1e-6 < 1e-6,
          'ratio = {0:.6e}, i.e. 1e-6 for a second-order rate'.format(ratio))
    note('rmgpy/chemkin.pyx: A=(kinetics.A.value_si, kinetics.A.units) -- an SI VALUE '
         'paired with the object\'s DECLARED units. Right only when the object already '
         'declares SI.')

    # ---------------------------------------------------------------- the live round trip
    banner('the live round trip: RMG writes, RMG reads, RMG writes again')
    with open(chem_path) as handle:
        original = handle.read()
    print('  as first written:')
    print(original[original.index('REACTIONS'):])

    species, reactions = load_chemkin_file(chem_path, dict_path)
    print('  read back:')
    for rxn in reactions:
        kin = rxn.kinetics
        print('    {0}'.format(rxn))
        print('      kinetics : {0}'.format(type(kin).__name__))
        print('      A        : {0}   (value_si = {1:.6e})'.format(kin.A, kin.A.value_si))
        print('      k(11604.5 K, Te=11604.5 K) = {0:.6e}'.format(
            kin.get_rate_coefficient(11604.5, 11604.5)
            if _accepts_two(kin) else kin.get_rate_coefficient(11604.5)))

    rewritten = os.path.join(out_dir, 'chemkin', 'chem_rewritten.inp')
    save_chemkin_file(rewritten, species, reactions)
    with open(rewritten) as handle:
        second = handle.read()
    print('  as written the SECOND time:')
    print(second[second.index('REACTIONS'):])

    first_As = _arrhenius_columns(original)
    second_As = _arrhenius_columns(second)
    note('A column, first write : {0}'.format(first_As))
    note('A column, second write: {0}'.format(second_As))
    check('the two files carry the same A column', first_As == second_As,
          'first {0} vs second {1}'.format(first_As, second_As))
    if first_As and second_As and len(first_As) == len(second_As):
        for a1, a2 in zip(first_As, second_As):
            note('    {0:.6e} -> {1:.6e}   ratio {2:.6e}'.format(a1, a2, a2 / a1))

    check('both files declare the same units in their REACTIONS header',
          _units_header(original) == _units_header(second),
          '{0!r} vs {1!r}'.format(_units_header(original), _units_header(second)))

    # ------------------------------------------------------- the rate the third trip sees
    banner('the rate a downstream consumer of the SECOND file would use')
    species2, reactions2 = load_chemkin_file(rewritten, dict_path)
    for rxn1, rxn2 in zip(reactions, reactions2):
        k1 = rxn1.kinetics.get_rate_coefficient(11604.5, 11604.5)
        k2 = rxn2.kinetics.get_rate_coefficient(11604.5, 11604.5)
        note('{0}'.format(rxn1))
        note('    k after one round trip  = {0:.6e}'.format(k1))
        note('    k after two round trips = {0:.6e}   ratio {1:.6e}'.format(k2, k2 / k1))
        check('the rate survives a second round trip', abs(k2 / k1 - 1.0) < 1e-9,
              'ratio {0:.6e}'.format(k2 / k1))

    print('\n' + '=' * 78)
    if FAILURES:
        print('FAILURES ({0}):'.format(len(FAILURES)))
        for label in FAILURES:
            print('  - {0}'.format(label))
        return 1
    print('all checks passed')
    return 0


def _accepts_two(kinetics):
    return isinstance(kinetics, _kinetics.TwoTemperaturePlasma)


def _units_header(text):
    for line in text.splitlines():
        if line.strip().startswith('REACTIONS'):
            return line.strip()
    return None


def _arrhenius_columns(text):
    values = []
    body = text[text.index('REACTIONS'):]
    for line in body.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith('!') or stripped.startswith('REACTIONS'):
            continue
        if stripped.upper().startswith('TDEP/') or stripped.upper().startswith('END'):
            continue
        parts = stripped.split()
        if len(parts) >= 4:
            try:
                values.append(float(parts[-3]))
            except ValueError:
                continue
    return values


if __name__ == '__main__':
    sys.exit(main())
