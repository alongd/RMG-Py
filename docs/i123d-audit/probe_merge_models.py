#!/usr/bin/env python
# encoding: utf-8

"""
``scripts/mergeModels.py``, unmodified, on two mechanisms each carrying one half of the
lithium charge network.

    python docs/i123d-audit/probe_merge_models.py <deck-output-dir>

The deck output directory must contain ``chemkin/chem.inp`` and
``chemkin/species_dictionary.txt``. The two halves are produced by
``docs/i145-units-laundering/make_half_mechanisms.py``, which partitions the REACTIONS
block textually and copies everything above it verbatim, so no number is touched on the
way in.

What this establishes, and why both halves of it are needed:

  1. BOTH CHANNELS SURVIVE. A merged model that carries only one of the two reactions is
     the I-134 shape -- the sink discarded as a duplicate of its own source -- seen through
     the shipped merge script rather than through the generator. The count alone is the
     check that would have caught it.

  2. BOTH RATE CONSTANTS ARE UNCHANGED. A merged model that carries both reactions with
     the wrong numbers is the I-145 shape -- a read-then-write that relabels a rate law's
     pre-exponential with units it is not in. ``mergeModels`` reads two Chemkin files and
     writes one, so it is a full round trip with a shipped script in the middle, and it is
     the exact invocation through which that defect was reachable with no plasma-specific
     flag. k(Te) is compared against the deck's own file, reaction by reaction, at the
     electron temperature the deck runs at.

The script is invoked as a SUBPROCESS, by path, exactly as a user would invoke it -- not
imported and not monkeypatched -- so that what is measured is the shipped entry point.

The resolved database directory is printed at the head of the run rather than trusted
from configuration.
"""

import os
import subprocess
import sys

from rmgpy import settings
from rmgpy.chemkin import load_chemkin_file

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TE = 11604.5  # K -- the deck's electron temperature; 1 eV in kelvin

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


def rate_at(kinetics, T, Te):
    """k for a rate law that may or may not know about an electron temperature.

    The order of these branches matters and was got wrong once. A Chemkin ``TDEP/``
    reaction reads back as a ``TwoTemperaturePlasma``, whose electron-temperature entry
    point is ``get_rate_coefficient_two_temp(T, Te)`` -- it has no
    ``get_rate_coefficient_electron_temp``, which is the name the two *library* rate laws
    (Voronov, Badnell) use. Probing only for the library spelling silently falls through
    to ``get_rate_coefficient(T)``, i.e. the gas temperature, and for the ionisation
    channel that is k = 8.6e-25 at 1000 K against 1.08e+08 at Te -- 33 orders of magnitude
    away from the number the reactor uses. It is still a consistent comparison, so the
    ratios come out right either way, but the printed label would have been a lie.
    """
    two_temp = getattr(kinetics, 'get_rate_coefficient_two_temp', None)
    if two_temp is not None:
        return two_temp(T, Te)
    electron_temp = getattr(kinetics, 'get_rate_coefficient_electron_temp', None)
    if electron_temp is not None:
        return electron_temp(Te)
    return kinetics.get_rate_coefficient(T)


def signature(rxn):
    """A stable, order-insensitive label for a reaction, for pairing across files."""
    def side(species_list):
        return '+'.join(sorted(str(s) for s in species_list))
    return '{0}=>{1}'.format(side(rxn.reactants), side(rxn.products))


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    deck = os.path.abspath(sys.argv[1])

    print('=' * 78)
    print('database.directory = {0}'.format(settings['database.directory']))
    print('resolved           = {0}'.format(
        os.path.realpath(settings['database.directory'])))
    print('repo               = {0}'.format(REPO))
    print('deck output        = {0}'.format(deck))
    print('=' * 78)

    chem = os.path.join(deck, 'chemkin', 'chem.inp')
    dic = os.path.join(deck, 'chemkin', 'species_dictionary.txt')
    if not (os.path.isfile(chem) and os.path.isfile(dic)):
        print('deck output is missing chemkin/chem.inp or chemkin/species_dictionary.txt')
        return 2

    banner('THE DECK ITSELF -- the reference both halves are measured against')
    deck_species, deck_reactions = load_chemkin_file(chem, dic)
    reference = {}
    for rxn in deck_reactions:
        reference[signature(rxn)] = rate_at(rxn.kinetics, 1000.0, TE)
        print('    {0:<44} k(Te={1:g}) = {2:.10e}'.format(
            str(rxn), TE, reference[signature(rxn)]))
    check('the deck carries both charge channels', len(deck_reactions) == 2,
          '{0} reaction(s)'.format(len(deck_reactions)))

    banner('SPLIT -- one reaction per half, textually, numbers untouched')
    split = subprocess.run(
        [sys.executable,
         os.path.join(REPO, 'docs', 'i145-units-laundering', 'make_half_mechanisms.py'),
         deck],
        cwd=REPO, capture_output=True, text=True)
    sys.stdout.write(split.stdout)
    sys.stderr.write(split.stderr)
    check('make_half_mechanisms exited 0', split.returncode == 0,
          'exit {0}'.format(split.returncode))

    half0 = os.path.join(deck, 'merge-half0')
    half1 = os.path.join(deck, 'merge-half1')
    for half in (half0, half1):
        h_species, h_reactions = load_chemkin_file(
            os.path.join(half, 'chem.inp'), os.path.join(half, 'species_dictionary.txt'))
        check('{0} carries exactly one reaction'.format(os.path.basename(half)),
              len(h_reactions) == 1, '{0} reaction(s)'.format(len(h_reactions)))
        for rxn in h_reactions:
            k = rate_at(rxn.kinetics, 1000.0, TE)
            ref = reference.get(signature(rxn))
            note('{0:<40} k = {1:.10e}'.format(str(rxn), k))
            check('splitting did not move {0}'.format(signature(rxn)),
                  ref is not None and k == ref,
                  'half {0:.10e} vs deck {1:.10e}'.format(k, ref if ref else float('nan')))

    banner('MERGE -- scripts/mergeModels.py, unmodified, as a subprocess')
    workdir = os.path.join(deck, 'merge-run')
    if not os.path.isdir(workdir):
        os.makedirs(workdir)
    cmd = [sys.executable, os.path.join(REPO, 'scripts', 'mergeModels.py'),
           '--model1', os.path.join(half0, 'chem.inp'),
           os.path.join(half0, 'species_dictionary.txt'),
           '--model2', os.path.join(half1, 'chem.inp'),
           os.path.join(half1, 'species_dictionary.txt')]
    print('    $ ' + ' '.join(cmd))
    merged = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    sys.stdout.write(merged.stdout)
    sys.stderr.write(merged.stderr)
    check('mergeModels exited 0', merged.returncode == 0,
          'exit {0}'.format(merged.returncode))

    out_chem = os.path.join(workdir, 'chem.inp')
    out_dic = os.path.join(workdir, 'species_dictionary.txt')
    if not os.path.isfile(out_chem):
        for root, _dirs, files in os.walk(workdir):
            if 'chem.inp' in files:
                out_chem = os.path.join(root, 'chem.inp')
                out_dic = os.path.join(root, 'species_dictionary.txt')
                break
    check('mergeModels wrote a Chemkin file', os.path.isfile(out_chem), out_chem)
    if not os.path.isfile(out_chem):
        return 1

    banner('THE MERGED MODEL -- both channels, and both rate constants')
    m_species, m_reactions = load_chemkin_file(out_chem, out_dic)
    for rxn in m_reactions:
        print('    {0:<44} k(Te={1:g}) = {2:.10e}'.format(
            str(rxn), TE, rate_at(rxn.kinetics, 1000.0, TE)))
    check('the merged model carries BOTH charge channels', len(m_reactions) == 2,
          '{0} reaction(s)'.format(len(m_reactions)))

    seen = set()
    for rxn in m_reactions:
        sig = signature(rxn)
        seen.add(sig)
        k = rate_at(rxn.kinetics, 1000.0, TE)
        ref = reference.get(sig)
        if ref is None:
            check('merged reaction {0} came from the deck'.format(sig), False,
                  'no deck counterpart')
            continue
        ratio = k / ref if ref else float('nan')
        check('{0} rate constant unchanged by the merge'.format(sig), k == ref,
              'merged {0:.10e} vs deck {1:.10e}   ratio {2:.10g}'.format(k, ref, ratio))
    for sig in reference:
        check('{0} survived the merge'.format(sig), sig in seen)

    print('\n' + '=' * 78)
    if FAILURES:
        print('MERGE MODELS: {0} CHECK(S) FAILED'.format(len(FAILURES)))
        for f in FAILURES:
            print('  - {0}'.format(f))
        return 1
    print('MERGE MODELS: both channels survived and neither rate constant moved')
    print('=' * 78)
    return 0


if __name__ == '__main__':
    sys.exit(main())
