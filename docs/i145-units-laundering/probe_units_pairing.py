#!/usr/bin/env python
# encoding: utf-8

"""
I-145 -- the laundered rate constant: does a mechanism survive being read and written?

Run from the repo root, so ``rmgrc`` resolves:

    python docs/i145-units-laundering/probe_units_pairing.py <deck-output-dir>

``<deck-output-dir>`` is the output directory of a run of
``docs/i123-integration/input.py`` -- it must contain ``chemkin/chem.inp`` and
``chemkin/species_dictionary.txt``.

Exit 0 when every check passes; 1 otherwise. Before the I-145 repair this probe exits 1
on sections C and D; after it, on nothing.

WHAT EACH SECTION ESTABLISHES

  A. THE READER IS CORRECT -- measured twice, independently, so that "the writer is the
     wrong half" is a measurement and not an assumption.

     A1 against the FILE. The file states its own convention in its header
        (``REACTIONS KCAL/MOLE MOLES``, Chemkin lengths in cm). Recompute
        ``k = A * Te**n * exp(-Ea/(R*Te))`` from the printed triple by hand under that
        convention, and compare with what the reloaded object returns. This checks the
        reader against the artifact, using no RMG machinery but the object under test.

        ``R`` is taken from :mod:`rmgpy.constants` rather than hardcoded. RMG carries the
        2010 CODATA value 8.314472 J/(mol*K); hardcoding the 2018 value
        8.31446261815324 makes this check fail by a ratio of 1.0000079, which is the two
        constants' disagreement and nothing about the reader. Measured: with Ea = 162.15
        kcal/mol at Te = 11604.5 K the two constants alone predict 1.0000079342 and the
        observed discrepancy was 1.000007934.

     A2 against the GENERATOR, to the printed digit. Load the two kinetics libraries the
        deck used, take the ORIGINAL ``VoronovEIArrhenius`` / ``BadnellRRArrhenius``
        entries, and push them through the very reduction the writer calls. Then assert,
        exactly rather than to a tolerance:

          - formatting the generator's reduction with the writer's own format string
            reproduces the file's printed triple character for character -- so the file
            faithfully carries the generated rate to the precision it prints;
          - the reloaded object's parameters, re-expressed in the file's convention,
            reproduce that same printed triple -- so the reader faithfully rebuilds it.

        Together those two make "the reader reproduces the generated rate to every
        printed digit" a measurement. The residual k(Te) discrepancy against the
        unrounded generator is reported and bounded, not asserted away: the file prints
        ``n`` to three decimals, so a rounding of 5e-4 in ``n`` moves k(Te) by
        ``exp(5e-4 * ln(11604.5))`` = 1.0047, i.e. up to 0.47%. The observed
        discrepancies are 0.23% and 0.38%, both inside that bound, so they are the
        artifact's precision and not the reader's error.

  B. THE MISPAIRING, DRIVEN DIRECTLY. The same rate law is built twice -- once declaring
     SI units, once declaring cgs -- and pushed through
     ``_plasma_arrhenius_for_chemkin``. A reduction that depends on which units string
     the source happens to declare is the defect stated in one measurement. The factor is
     10**6 because the rate coefficient carries m**3 in its dimensionality and
     1 m**3 = (100 cm)**3 = 10**6 cm**3; a third-order rate carries m**6 and would be off
     by 10**12.

  C. THREE SUCCESSIVE ROUND TRIPS. Write, read, write, read, write, read. The rate
     constant must be identical at every stage. Before the repair each trip multiplies by
     10**-6 and the third file is 10**-18 of the first, silently, at exit status 0.

  D. THE UNITS HEADER IS HONOURED. Every file in the chain declares the same header, so
     the numbers under it must also be the same. This is the check that forbids the
     wrong repair -- relabelling the header to match wrong numbers would satisfy C in
     isolation but not C and D together.

  E. NEGATIVE CONTROL. An ordinary non-plasma mechanism -- no electrons, no TDEP, plain
     Arrhenius -- must round-trip byte-identically, and the repair must not move a single
     byte of it. The probe writes it three times and compares; the "unchanged by the
     repair" half is the separate before/after comparison of the emitted files, which the
     ticket's evidence records.
"""

import math
import os
import shutil
import sys

import rmgpy.constants as constants
import rmgpy.kinetics as _kinetics
from rmgpy import electron_placement
from rmgpy import settings
from rmgpy.chemkin import (load_chemkin_file, save_chemkin_file,
                           save_species_dictionary, _plasma_arrhenius_for_chemkin)
from rmgpy.data.kinetics.database import KineticsDatabase

R = constants.R       # J/(mol*K) -- RMG's own value, not a hardcoded one; see A1 above
T_GAS = 11604.5       # K -- the deck's electron temperature; along T = Te the
TE = 11604.5          #      TwoTemperaturePlasma reduction is exact

# The exact format string rmgpy/chemkin.pyx::write_kinetics_entry uses for the
# modified-Arrhenius triple. Reused here so the comparison is against the writer's own
# precision rather than against a precision this probe invented.
TRIPLE_FORMAT = '{0:<9.3e} {1:<9.3f} {2:<9.3f}'

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


def note(text):
    print('  [..] {0}'.format(text))


def banner(text):
    print('\n' + '=' * 78)
    print(text)
    print('=' * 78)


def parse_reaction_lines(text):
    """
    Return ``[(equation, A, n, Ea, triple), ...]`` from a Chemkin REACTIONS block, in
    file order.

    ``A`` is as printed, i.e. in the file's own cm/mol/s convention; ``Ea`` as printed,
    i.e. kcal/mol under a ``KCAL/MOLE`` header. ``triple`` is the three printed fields
    re-joined in the writer's own format, for character-exact comparison.
    """
    entries = []
    body = text[text.index('REACTIONS'):]
    for line in body.splitlines():
        stripped = line.strip()
        if (not stripped or stripped.startswith('!')
                or stripped.startswith('REACTIONS')
                or stripped.upper().startswith('TDEP/')
                or stripped.upper().startswith('END')):
            continue
        parts = stripped.split()
        if len(parts) < 4:
            continue
        try:
            A, n, Ea = float(parts[-3]), float(parts[-2]), float(parts[-1])
        except ValueError:
            continue
        entries.append((' '.join(parts[:-3]), A, n, Ea,
                        TRIPLE_FORMAT.format(A, n, Ea)))
    return entries


def units_header(text):
    for line in text.splitlines():
        if line.strip().startswith('REACTIONS'):
            return line.strip()
    return None


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    out_dir = sys.argv[1]
    chem_path = os.path.join(out_dir, 'chemkin', 'chem.inp')
    dict_path = os.path.join(out_dir, 'chemkin', 'species_dictionary.txt')

    import rmgpy.molecule.molecule as _m
    print('=' * 78)
    print('rmgpy              = {0}'.format(os.path.abspath(
        os.path.dirname(os.path.dirname(electron_placement.__file__)))))
    print('compiled module    = {0}'.format(_m.__file__))
    print('database.directory = {0}'.format(settings['database.directory']))
    print('resolved realpath  = {0}'.format(os.path.realpath(settings['database.directory'])))
    print('chemkin file       = {0}'.format(chem_path))
    print('=' * 78)

    with open(chem_path) as handle:
        original = handle.read()
    printed = parse_reaction_lines(original)
    species, reactions = load_chemkin_file(chem_path, dict_path)

    # ------------------------------------------------------------------------ section A1
    banner('A1. THE READER, CHECKED AGAINST THE FILE\'S OWN STATED CONVENTION')
    note('header: {0!r}'.format(units_header(original)))
    note('so the printed triple means k = A * Te**n * exp(-Ea/(R*Te)) with A in')
    note('cm**3/(mol*s) for a second-order rate and Ea in kcal/mol.')
    note('R = {0!r} J/(mol*K), taken from rmgpy.constants'.format(R))
    for rxn, (equation, A_cm, n, Ea_kcal, _triple) in zip(reactions, printed):
        k_hand_cgs = A_cm * TE ** n * math.exp(-Ea_kcal * 4184.0 / (R * TE))
        k_hand_si = k_hand_cgs / 1.0e6  # cm**3 -> m**3 for a second-order rate
        k_obj = rxn.kinetics.get_rate_coefficient(T_GAS, TE)
        note('{0}'.format(equation))
        note('    printed        : A = {0:.6e} cm^3/(mol*s), n = {1}, Ea = {2} kcal/mol'
             .format(A_cm, n, Ea_kcal))
        note('    reloaded object: A = {0}   value_si = {1:.6e}'
             .format(rxn.kinetics.A, rxn.kinetics.A.value_si))
        note('    k by hand from the file  = {0:.6e} m^3/(mol*s)'.format(k_hand_si))
        note('    k from the reloaded obj  = {0:.6e} m^3/(mol*s)'.format(k_obj))
        check('reloaded k matches hand arithmetic on the printed triple: {0}'.format(equation),
              abs(k_obj / k_hand_si - 1.0) < 1e-9,
              'ratio {0:.9f}'.format(k_obj / k_hand_si))
        check('the reloaded A declares the file\'s units, not SI: {0}'.format(equation),
              rxn.kinetics.A.units == 'cm^3/(mol*s)',
              'units = {0!r}'.format(rxn.kinetics.A.units))

    # ------------------------------------------------------------------------ section A2
    banner('A2. THE READER, CHECKED AGAINST THE GENERATOR, TO THE PRINTED DIGIT')
    db = KineticsDatabase()
    db.load_libraries(os.path.join(settings['database.directory'], 'kinetics', 'libraries'),
                      libraries=[IONISATION, RECOMBINATION])
    generated = {}
    for lib_name in (IONISATION, RECOMBINATION):
        library = db.libraries[lib_name]
        for entry in library.entries.values():
            kin = entry.data
            if not isinstance(kin, (_kinetics.VoronovEIArrhenius,
                                    _kinetics.BadnellRRArrhenius)):
                continue
            arrhenius, _note = _plasma_arrhenius_for_chemkin(kin)
            # Exactly what write_kinetics_entry prints for this reduction.
            triple = TRIPLE_FORMAT.format(
                arrhenius.A.value_si / (arrhenius.T0.value_si ** arrhenius.n.value_si)
                * arrhenius.A.get_conversion_factor_from_si_to_cm_mol_s(),
                arrhenius.n.value_si,
                arrhenius.Ea.value_si / 4184.0)
            generated[type(kin).__name__] = (
                str(entry.item), arrhenius.get_rate_coefficient(TE), triple)
    note('the two library entries, reduced by the very function the writer calls,')
    note('then formatted with the writer\'s own format string:')
    for cls_name, (item, k_gen, triple) in sorted(generated.items()):
        note('    {0:<22s} {1:<20s} {2}   k(Te) = {3:.6e}'.format(
            cls_name, item, triple, k_gen))

    # The deck's first reaction is the ionisation (Voronov), the second the
    # recombination (Badnell); the file's own TDEP notes name them.
    expected_order = ['VoronovEIArrhenius', 'BadnellRRArrhenius']
    check('both library rate laws were found', len(generated) == 2,
          'found {0}'.format(sorted(generated)))
    if len(generated) == 2:
        for rxn, cls_name, entry in zip(reactions, expected_order, printed):
            _equation, _A, _n, _Ea, file_triple = entry
            _item, k_gen, gen_triple = generated[cls_name]
            k_read = rxn.kinetics.get_rate_coefficient(T_GAS, TE)

            note('{0}'.format(rxn))
            note('    triple the GENERATOR reduces to : {0!r}'.format(gen_triple))
            note('    triple the FILE carries         : {0!r}'.format(file_triple))
            check('the file carries the generated triple, character for character: '
                  '{0}'.format(rxn), gen_triple == file_triple,
                  '{0!r} vs {1!r}'.format(gen_triple, file_triple))

            # Re-express the RELOADED object in the file's own convention and format it
            # the same way. This is the reader's output measured on the writer's terms.
            kin = rxn.kinetics
            read_triple = TRIPLE_FORMAT.format(
                kin.A.value_si / (kin.T0.value_si ** kin.n.value_si)
                * kin.A.get_conversion_factor_from_si_to_cm_mol_s(),
                kin.n.value_si,
                kin.Ea_g.value_si / 4184.0)
            note('    triple the RELOADED object gives: {0!r}'.format(read_triple))
            check('the reader rebuilds the file\'s triple, character for character: '
                  '{0}'.format(rxn), read_triple == file_triple,
                  '{0!r} vs {1!r}'.format(read_triple, file_triple))

            # And the residual against the UNROUNDED generator, bounded rather than
            # asserted away: n is printed to three decimals, so a 5e-4 rounding in n
            # moves k(Te) by exp(5e-4*ln(Te)).
            bound = math.exp(5e-4 * math.log(TE)) - 1.0
            ratio = k_read / k_gen
            note('    k as GENERATED (unrounded library object) = {0:.6e}'.format(k_gen))
            note('    k as READ BACK (from the file)            = {0:.6e}'.format(k_read))
            note('    ratio = {0:.9f}; bound from printing n to 3 decimals = {1:.9f}'
                 .format(ratio, 1.0 + bound))
            check('the residual against the unrounded generator is inside the bound the '
                  'file\'s own printed precision sets: {0}'.format(rxn),
                  abs(ratio - 1.0) <= bound,
                  'residual {0:.3e} vs bound {1:.3e}'.format(abs(ratio - 1.0), bound))
    note('')
    note('THEREFORE: the read-back half is correct -- it rebuilds every digit the file')
    note('carries, and the file carries every digit the generator produced. Anything the')
    note('round trip loses is lost on the way OUT, which is what the rest of this probe')
    note('measures.')

    # ------------------------------------------------------------------------- section B
    banner('B. THE WRITER\'S UNIT PAIRING, DRIVEN DIRECTLY')
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

    arr_si, _ = _plasma_arrhenius_for_chemkin(si)
    arr_cgs, _ = _plasma_arrhenius_for_chemkin(cgs)
    note('reduced for Chemkin:')
    note('    from the SI-declared object  : A = {0}   value_si = {1:.6e}'.format(
        arr_si.A, arr_si.A.value_si))
    note('    from the cgs-declared object : A = {0}   value_si = {1:.6e}'.format(
        arr_cgs.A, arr_cgs.A.value_si))
    ratio = arr_cgs.A.value_si / arr_si.A.value_si
    note('    ratio of the two reductions  : {0:.6e}   (1 m**3 = 10**6 cm**3)'.format(ratio))
    check('the reduction is UNIT-INDEPENDENT: the same rate law reduces to the same '
          'number however its units are declared',
          abs(ratio - 1.0) < 1e-12,
          'ratio = {0:.9e}'.format(ratio))
    check('the reduction preserves the source rate law exactly (SI-declared source)',
          abs(arr_si.A.value_si / si.A.value_si - 1.0) < 1e-12,
          '{0:.6e} vs {1:.6e}'.format(arr_si.A.value_si, si.A.value_si))
    check('the reduction preserves the source rate law exactly (cgs-declared source)',
          abs(arr_cgs.A.value_si / cgs.A.value_si - 1.0) < 1e-12,
          '{0:.6e} vs {1:.6e}'.format(arr_cgs.A.value_si, cgs.A.value_si))

    # ---------------------------------------------------------------------- sections C+D
    banner('C+D. THREE SUCCESSIVE ROUND TRIPS')
    scratch = os.path.join(out_dir, 'roundtrip')
    if os.path.isdir(scratch):
        shutil.rmtree(scratch)
    os.makedirs(scratch)

    texts = [original]
    rate_sets = [[rxn.kinetics.get_rate_coefficient(T_GAS, TE) for rxn in reactions]]
    current_species, current_reactions = species, reactions
    for trip in range(1, 4):
        path = os.path.join(scratch, 'chem_trip{0}.inp'.format(trip))
        save_chemkin_file(path, current_species, current_reactions)
        with open(path) as handle:
            texts.append(handle.read())
        current_species, current_reactions = load_chemkin_file(path, dict_path)
        rate_sets.append([rxn.kinetics.get_rate_coefficient(T_GAS, TE)
                          for rxn in current_reactions])

    labels = ['as generated'] + ['after trip {0}'.format(i) for i in (1, 2, 3)]
    for label, text in zip(labels, texts):
        note('{0:<14s} header {1!r}'.format(label, units_header(text)))
        for equation, A, n, Ea, _triple in parse_reaction_lines(text):
            note('    {0:<40s} {1:>11.3e} {2:>8.3f} {3:>10.3f}'.format(equation, A, n, Ea))

    for i, label in enumerate(labels):
        note('{0:<14s} k = {1}'.format(
            label, ', '.join('{0:.6e}'.format(k) for k in rate_sets[i])))

    base_rates = rate_sets[0]
    for trip in (1, 2, 3):
        these = rate_sets[trip]
        ratios = [b / a for a, b in zip(base_rates, these)]
        check('the rate constant survives round trip {0}'.format(trip),
              all(abs(r - 1.0) < 1e-12 for r in ratios),
              'ratios ' + ', '.join('{0:.9e}'.format(r) for r in ratios))

    base_A = [entry[1] for entry in parse_reaction_lines(texts[0])]
    for trip in (1, 2, 3):
        these_A = [entry[1] for entry in parse_reaction_lines(texts[trip])]
        check('the printed A column survives round trip {0}'.format(trip),
              base_A == these_A,
              '{0} vs {1}'.format(base_A, these_A))

    headers = [units_header(text) for text in texts]
    check('every file in the chain declares the same units header',
          len(set(headers)) == 1, '{0}'.format(set(headers)))

    # ------------------------------------------------------------------------- section E
    banner('E. NEGATIVE CONTROL: AN ORDINARY NON-PLASMA MECHANISM')
    control_dir = os.path.join('test', 'rmgpy', 'test_data', 'chemkin',
                               'chemkin_py', 'minimal')
    control_chem = os.path.join(control_dir, 'chem.inp')
    control_dict = os.path.join(control_dir, 'species_dictionary.txt')
    note('mechanism: {0}'.format(os.path.abspath(control_chem)))
    c_species, c_reactions = load_chemkin_file(control_chem, control_dict)
    note('{0} species, {1} reactions; kinetics classes: {2}'.format(
        len(c_species), len(c_reactions),
        sorted({type(rxn.kinetics).__name__ for rxn in c_reactions})))
    check('the control mechanism carries no plasma rate law',
          not any(isinstance(rxn.kinetics, (_kinetics.TwoTemperaturePlasma,
                                            _kinetics.VoronovEIArrhenius,
                                            _kinetics.BadnellRRArrhenius,
                                            _kinetics.ElectronCollisionPlasma))
                  for rxn in c_reactions),
          'classes {0}'.format(sorted({type(r.kinetics).__name__ for r in c_reactions})))
    check('the control mechanism carries no electrons',
          all(getattr(rxn, 'electrons', 0) == 0 for rxn in c_reactions))

    control_scratch = os.path.join(out_dir, 'negative-control')
    if os.path.isdir(control_scratch):
        shutil.rmtree(control_scratch)
    os.makedirs(control_scratch)
    control_texts = []
    cur_s, cur_r = c_species, c_reactions
    for trip in range(1, 4):
        path = os.path.join(control_scratch, 'chem_trip{0}.inp'.format(trip))
        # The written file uses the writer's own species identifiers, so the dictionary
        # has to be re-emitted alongside it; reusing the input fixture's dictionary
        # leaves species unresolved on the way back in.
        dict_out = os.path.join(control_scratch, 'dict_trip{0}.txt'.format(trip))
        save_chemkin_file(path, cur_s, cur_r)
        save_species_dictionary(dict_out, cur_s)
        with open(path, 'rb') as handle:
            control_texts.append(handle.read())
        cur_s, cur_r = load_chemkin_file(path, dict_out)
    for trip in (2, 3):
        check('the control mechanism is byte-identical after round trip {0}'.format(trip),
              control_texts[trip - 1] == control_texts[0],
              '{0} bytes vs {1} bytes'.format(len(control_texts[trip - 1]),
                                              len(control_texts[0])))
    note('control file digest (for the before/after comparison): {0}'.format(
        __import__('hashlib').sha256(control_texts[0]).hexdigest()))

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
