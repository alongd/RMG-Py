#!/usr/bin/env python
# encoding: utf-8

"""
STAGE B of the i123e end-to-end measurement: THREE successive Chemkin round trips.

Starting from the deck's own chemkin/chem.inp + species_dictionary.txt (+ tran.dat),
perform load -> write -> load -> write -> load -> write. At EACH of the three
readbacks, for BOTH charge channels:

  * assert the deserialized rate-law CLASS explicitly (the rate-law trap: a Chemkin
    TDEP/ reaction reads back as TwoTemperaturePlasma, whose electron-temperature
    entry point is get_rate_coefficient_two_temp(T, Te); the library classes
    VoronovEIArrhenius/BadnellRRArrhenius use get_rate_coefficient_electron_temp(Te);
    probing the wrong spelling silently falls through to get_rate_coefficient(T) at
    the GAS temperature, ~33 orders off, with ratios that still pass);
  * evaluate k at (T=1000 K, Te=11604.5 K) with the class-correct method;
  * print class + absolute value.

Values are compared across all three trips and against the canonical
(pre-serialisation) library values from PlasmaElectronImpactIonization and
PlasmaRadiativeRecombination. Any per-trip drift (the historical blocker was every
rate 1e6x too small, compounding per trip) is flagged.

Usage: python docs/i123e-audit/evidence/probe_three_roundtrips.py <deck-output-dir>
Round-trip artifacts are written under <deck-output-dir>/../roundtrips/.
"""

import os
import sys

from rmgpy import settings
from rmgpy.chemkin import (load_chemkin_file, save_chemkin_file,
                           save_species_dictionary, save_transport_file)
from rmgpy.data.kinetics.database import KineticsDatabase

T = 1000.0
TE = 11604.5
IONISATION = 'PlasmaElectronImpactIonization'
RECOMBINATION = 'PlasmaRadiativeRecombination'

FAILURES = []


def check(label, condition, detail=''):
    ok = bool(condition)
    print('  [{0}] {1}{2}'.format('PASS' if ok else 'FAIL', label,
                                  (' -- ' + detail) if detail else ''))
    if not ok:
        FAILURES.append(label)
    return ok


def class_correct_k(kinetics):
    """Return (method_name, k) selecting the entry point by CLASS, loudly."""
    name = type(kinetics).__name__
    if name == 'TwoTemperaturePlasma':
        return 'get_rate_coefficient_two_temp(T, Te)', \
            kinetics.get_rate_coefficient_two_temp(T, TE)
    if name in ('VoronovEIArrhenius', 'BadnellRRArrhenius'):
        return 'get_rate_coefficient_electron_temp(Te)', \
            kinetics.get_rate_coefficient_electron_temp(TE)
    raise AssertionError('unexpected rate-law class at readback: {0}'.format(name))


def channel_of(rxn):
    """'ionisation' or 'recombination'.

    The canonical library form carries its electrons as a SCALAR (1 reactant, 1
    product on both channels), so participant counts cannot distinguish them there;
    the rate-law class can. The reloaded Chemkin form is TwoTemperaturePlasma on
    both channels, so there the participant counts (explicit electrons) decide.
    """
    cls = type(rxn.kinetics).__name__
    if cls == 'VoronovEIArrhenius':
        return 'ionisation'
    if cls == 'BadnellRRArrhenius':
        return 'recombination'
    if len(rxn.products) > len(rxn.reactants):
        return 'ionisation'
    return 'recombination'


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    deck = os.path.abspath(sys.argv[1])
    work = os.path.join(os.path.dirname(deck), 'roundtrips')

    print('=' * 78)
    print('database.directory = {0}'.format(settings['database.directory']))
    print('deck output        = {0}'.format(deck))
    print('roundtrip workdir  = {0}'.format(work))
    print('T = {0} K   Te = {1} K'.format(T, TE))
    print('=' * 78)

    # ---- canonical, pre-serialisation library values -------------------------
    print('\nCANONICAL LIBRARY VALUES (pre-serialisation)')
    db = KineticsDatabase()
    db.load_libraries(os.path.join(settings['database.directory'],
                                   'kinetics', 'libraries'),
                      libraries=[IONISATION, RECOMBINATION])
    canonical = {}
    for label in (IONISATION, RECOMBINATION):
        for rxn in db.libraries[label].get_library_reactions():
            cls = type(rxn.kinetics).__name__
            method, k = class_correct_k(rxn.kinetics)
            ch = channel_of(rxn)
            canonical[ch] = k
            print('  {0:<14} {1:<36} class={2:<22} {3}'.format(ch, str(rxn), cls, method))
            print('  {0:<14} k(T={1:g}, Te={2:g}) = {3:.10e}'.format('', T, TE, k))

    # ---- three successive trips ---------------------------------------------
    chem = os.path.join(deck, 'chemkin', 'chem.inp')
    dic = os.path.join(deck, 'chemkin', 'species_dictionary.txt')
    tran = os.path.join(deck, 'chemkin', 'tran.dat')

    table = []  # (trip, channel, class, k)
    for trip in (1, 2, 3):
        print('\n' + '-' * 78)
        print('TRIP {0}: load {1}'.format(trip, chem))
        species, reactions = load_chemkin_file(chem, dictionary_path=dic,
                                               transport_path=tran)
        print('  loaded {0} species, {1} reactions'.format(len(species), len(reactions)))
        check('trip {0}: exactly 2 reactions read back'.format(trip),
              len(reactions) == 2, '{0}'.format(len(reactions)))
        seen = set()
        for rxn in reactions:
            cls = type(rxn.kinetics).__name__
            check('trip {0}: rate-law class asserted explicitly'.format(trip),
                  cls in ('TwoTemperaturePlasma', 'VoronovEIArrhenius',
                          'BadnellRRArrhenius'), cls)
            method, k = class_correct_k(rxn.kinetics)
            ch = channel_of(rxn)
            seen.add(ch)
            table.append((trip, ch, cls, k))
            print('  {0:<14} {1:<36} class={2:<22} {3}'.format(ch, str(rxn), cls, method))
            print('  {0:<14} k = {1:.10e}   electrons={2}'.format(
                '', k, getattr(rxn, 'electrons', 'n/a')))
        check('trip {0}: both channels present'.format(trip),
              seen == {'ionisation', 'recombination'}, str(sorted(seen)))

        # write for the next trip
        outdir = os.path.join(work, 'trip{0}'.format(trip))
        os.makedirs(outdir, exist_ok=True)
        next_chem = os.path.join(outdir, 'chem.inp')
        next_dic = os.path.join(outdir, 'species_dictionary.txt')
        next_tran = os.path.join(outdir, 'tran.dat')
        save_chemkin_file(next_chem, species, reactions)
        save_species_dictionary(next_dic, species)
        save_transport_file(next_tran, species)
        print('  wrote {0}'.format(next_chem))
        chem, dic, tran = next_chem, next_dic, next_tran

    # ---- table and drift analysis -------------------------------------------
    print('\n' + '=' * 78)
    print('TABLE: trip | channel | class | k(T=1000, Te=11604.5)')
    for trip, ch, cls, k in table:
        print('  {0}    | {1:<14} | {2:<22} | {3:.10e}'.format(trip, ch, cls, k))

    print('\nDRIFT ANALYSIS (vs canonical library value)')
    for ch in ('ionisation', 'recombination'):
        ks = [k for (t, c, cls, k) in table if c == ch]
        ref = canonical[ch]
        print('  {0}: canonical = {1:.10e}'.format(ch, ref))
        for i, k in enumerate(ks, 1):
            ratio = k / ref
            print('    trip {0}: k = {1:.10e}   k/canonical = {2:.10g}'.format(i, k, ratio))
        check('{0}: no drift across trips (trip1 == trip2 == trip3)'.format(ch),
              ks[0] == ks[1] == ks[2],
              ' / '.join('{0:.10e}'.format(k) for k in ks))
        # The offset of the readback from the canonical library form is NOT drift:
        # it is the documented lossy export fit. chemkin.pyx's
        # _plasma_arrhenius_for_chemkin() reduces the pure-Te plasma forms via
        # to_arrhenius(), a modified-Arrhenius fit of k(Te) over the class's own
        # validity window, and writes that loss into the file as an inline note.
        # Report the one-time fit deviation, and require that it is incurred ONCE
        # (at export) rather than compounding per trip.
        rel = abs(ks[0] - ref) / ref
        print('    one-time export-fit deviation from canonical: {0:.3%} '
              '(documented lossy fit, incurred at export, not per trip)'.format(rel))
        check('{0}: fit loss incurred once, not compounded (trip3 == trip1)'.format(ch),
              ks[2] == ks[0])
        check('{0}: NOT the historical 1e6x-per-trip collapse'.format(ch),
              all(abs(k / ref) > 1e-3 for k in ks))

    print('\n' + '=' * 78)
    if FAILURES:
        print('STAGE B: {0} CHECK(S) FAILED'.format(len(FAILURES)))
        for f in FAILURES:
            print('  - {0}'.format(f))
        return 1
    print('STAGE B: three round trips stable; classes asserted; no drift')
    print('=' * 78)
    return 0


if __name__ == '__main__':
    sys.exit(main())
