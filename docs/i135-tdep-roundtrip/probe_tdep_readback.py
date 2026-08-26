#!/usr/bin/env python3
"""
I-135 probe: can the Chemkin file RMG just wrote be read back?

Usage:  python docs/i135-tdep-roundtrip/probe_tdep_readback.py <deck-output-dir>

Run from the repo root so ``rmgrc`` resolves; the resolved database directory is
printed at the head of the run rather than trusted from configuration.

Three trips, the same three the I-123 re-audit measured, so the before/after
comparison is like for like:

  A. RMG's own ``load_chemkin_file`` on ``<out>/chemkin/chem.inp`` as written.
  B. the same file with the ``TDEP/.../`` auxiliary line removed -- the control
     that isolates the TDEP line as the whole of the obstacle.
  C. Cantera's ``ck2yaml`` on the same file, the ecosystem reader this project
     does not own.

For trip A the probe also reports what the reconstructed rate law *is*: whether
it is evaluated at the electron temperature, and how its k compares with the
modified-Arrhenius triple actually written on the reaction line. A trip that
"succeeds" by reading the reaction back as a plain gas-temperature Arrhenius is
a worse outcome than a refusal, so the probe measures that explicitly instead of
reporting a bare READ.
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


def read_with_rmg(chem_path, dict_path):
    try:
        species, reactions = load_chemkin_file(chem_path, dict_path)
    except Exception as exc:
        return None, None, exc
    return species, reactions, None


def describe(reactions):
    for rxn in reactions:
        kin = rxn.kinetics
        print('  reaction  : {0}'.format(rxn))
        print('  kinetics  : {0}'.format(type(kin).__name__))
        print('  uses Te   : {0}'.format(getattr(kin, 'uses_electron_temperature', False)))
        print('  electrons : {0}'.format(rxn.electrons))
        for T in (11604.5, 20000.0, 100000.0):
            try:
                print('    k({0:>9.1f} K) = {1:.6e}'.format(T, kin.get_rate_coefficient(T)))
            except Exception as exc:  # pragma: no cover - diagnostic only
                print('    k({0:>9.1f} K) raised {1}'.format(T, exc))


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    out_dir = sys.argv[1]
    chem_path = os.path.join(out_dir, 'chemkin', 'chem.inp')
    dict_path = os.path.join(out_dir, 'chemkin', 'species_dictionary.txt')

    print('database.directory (resolved) : {0}'.format(settings['database.directory']))
    print('rmgpy package                 : {0}'.format(
        os.path.dirname(os.path.dirname(os.path.abspath(load_chemkin_file.__code__.co_filename)))))
    print('chemkin file                  : {0}'.format(chem_path))

    banner('the reaction lines as written')
    with open(chem_path) as f:
        text = f.read()
    in_reactions = False
    for line in text.splitlines():
        if line.strip().startswith('REACTIONS'):
            in_reactions = True
        if in_reactions and line.strip():
            print('  ' + line.rstrip())

    banner('A. RMG load_chemkin_file on the file as written')
    species_a, reactions_a, err_a = read_with_rmg(chem_path, dict_path)
    if err_a is None:
        print('  READ: {0} species, {1} reactions'.format(len(species_a), len(reactions_a)))
        describe(reactions_a)
    else:
        print('  REFUSED: {0}: {1}'.format(type(err_a).__name__, err_a))
        traceback.print_exc(file=sys.stdout)

    banner('B. the same file with the TDEP auxiliary line removed')
    stripped = os.path.join(out_dir, 'chemkin', 'chem_no_tdep.inp')
    kept = [ln for ln in text.splitlines(True) if not ln.strip().upper().startswith('TDEP/')]
    with open(stripped, 'w') as f:
        f.writelines(kept)
    print('  removed {0} TDEP line(s)'.format(len(text.splitlines()) - len(kept)))
    species_b, reactions_b, err_b = read_with_rmg(stripped, dict_path)
    if err_b is None:
        print('  READ: {0} species, {1} reactions'.format(len(species_b), len(reactions_b)))
        describe(reactions_b)
    else:
        print('  REFUSED: {0}: {1}'.format(type(err_b).__name__, err_b))

    banner('C. Cantera ck2yaml on the file as written')
    try:
        from cantera import ck2yaml
        import cantera
        print('  cantera {0} from {1}'.format(cantera.__version__, cantera.__file__))
        out_yaml = os.path.join(out_dir, 'chemkin', 'probe_ck2yaml.yaml')
        ck2yaml.convert_mech(chem_path, out_name=out_yaml, quiet=True, permissive=True)
        print('  CONVERTED -> {0}'.format(out_yaml))
    except Exception as exc:
        print('  REFUSED: {0}: {1}'.format(type(exc).__name__, str(exc).splitlines()[-1]))

    banner('summary')
    print('  A  RMG reader, file as written           : {0}'.format('READ' if err_a is None else 'REFUSED'))
    print('  B  RMG reader, TDEP line removed         : {0}'.format('READ' if err_b is None else 'REFUSED'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
