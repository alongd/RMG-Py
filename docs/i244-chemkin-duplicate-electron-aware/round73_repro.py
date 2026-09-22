"""
Round 73's three reported regressions, reproduced from scratch on either arm.

Self-contained: it builds its own species and reactions, so it can be run against
any checkout of RMG-Py without a scratchpad helper. Which arm it is running on is
decided by a *value*, not by a path: the branch defines
``chemkin_duplicate_group_key``, the base does not, and the banner also names which
function currently holds the group-level authority.

Run it from the root of the checkout under test:

    python docs/i244-chemkin-duplicate-electron-aware/round73_repro.py

Defect 1 is asserted by loading the deck with ``cantera.Solution``, not by
comparing deck text and not by ``ck2yaml.convert_mech`` alone: conversion writes
the YAML, and ``Kinetics::checkDuplicates`` only runs when the mechanism is
loaded. The rendered text differs between the two entries either way
(``H(1)+OH(2)=>...`` against ``OH(2)+H(1)=>...``), so a text comparison sees
nothing at all.
"""
import os
import tempfile
from types import SimpleNamespace

import rmgpy.chemkin as cm
from rmgpy.chemkin import render_chemkin_file, save_chemkin
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial

ARM = 'BRANCH' if hasattr(cm, 'chemkin_duplicate_group_key') else 'BASE'
print('arm                    :', ARM)
print('chemkin module         :', cm.__file__)
# Which engine is loaded, by a value rather than a count. The group-level authority
# is `chemkin_duplicate_flags` once Defect 2 is repaired and `mark_duplicate_reactions`
# before that, so ask both and print what answered -- a probe that names only the
# older one reads as BASE on a repaired engine.
for _name in ('chemkin_duplicate_flags', 'mark_duplicate_reactions'):
    _doc = getattr(getattr(cm, _name, None), '__doc__', '') or ''
    if 'group-level' in _doc:
        print('group-level authority  :', _name)
        break
else:
    print('group-level authority  : none (pairwise engine)')


def spc(label, index, smiles):
    """An RMG Species carrying the thermo the Chemkin writer needs."""
    s = Species(label=label, molecule=[Molecule(smiles=smiles)])
    s.index = index
    coeffs = [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, -11.7246]
    s.thermo = NASA(
        polynomials=[NASAPolynomial(coeffs=coeffs, Tmin=(200, 'K'), Tmax=(1000, 'K')),
                     NASAPolynomial(coeffs=coeffs, Tmin=(1000, 'K'), Tmax=(6000, 'K'))],
        Tmin=(200, 'K'), Tmax=(6000, 'K'))
    return s


H = spc('H', 1, '[H]')
OH = spc('OH', 2, '[OH]')
H2 = spc('H2', 3, '[H][H]')
O = spc('O', 4, '[O]')
SPECIES = [H, OH, H2, O]


def rxn(reactants, A, dup=False):
    r = LibraryReaction(reactants=list(reactants), products=[H2, O],
                        kinetics=Arrhenius(A=(A, 'cm^3/(mol*s)'), n=0, Ea=(0, 'kcal/mol')),
                        library='SomeLibrary', reversible=False)
    r.duplicate = dup
    return r


def equations(text):
    return [line.split()[0] for line in text.splitlines()
            if '=>' in line and not line.lstrip().startswith('!')]


def cantera_verdict(text):
    """Convert the deck AND load it. Only the load runs ``checkDuplicates``."""
    import cantera as ct
    import cantera.ck2yaml as ck2yaml
    with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', '/tmp')) as d:
        inp = os.path.join(d, 'chem.inp')
        out = os.path.join(d, 'chem.yaml')
        with open(inp, 'w') as f:
            f.write(text)
        try:
            ck2yaml.convert_mech(input_file=inp, out_name=out, quiet=True)
        except Exception as exc:
            return 'ck2yaml.convert_mech REJECTED -- %s' % str(exc).strip().splitlines()[0][:110]
        converted = 'ck2yaml.convert_mech ACCEPTED'
        try:
            ct.Solution(out)
        except Exception as exc:
            line = ([l for l in str(exc).splitlines() if 'uplicate' in l]
                    or [str(exc)[:110]])[0].strip()
            return '%s; cantera.Solution REJECTED -- %s' % (converted, line)
        return '%s; cantera.Solution LOADED (deck is valid)' % converted


print('\n--- DEFECT 1: a pre-marked pair whose reactant LISTS differ only in order ---')
a = rxn([H, OH], 1e12, dup=True)
b = rxn([OH, H], 2e12, dup=True)
text = render_chemkin_file(SPECIES, [a, b], verbose=False)
eqs = equations(text)
print('flags after the writer :', (a.duplicate, b.duplicate))
print('rendered equations     :', eqs)
print('identical as text?     :', len(eqs) == 2 and eqs[0] == eqs[1])
print('DUPLICATE lines        :', text.count('DUPLICATE'))
print('VERDICT                :', cantera_verdict(text))

print('\n--- DEFECT 3: a generator argument ---')
try:
    from_list = render_chemkin_file(SPECIES, [rxn([H, OH], 1e12)], verbose=False)
    from_gen = render_chemkin_file(SPECIES, (x for x in [rxn([H, OH], 1e12)]), verbose=False)
    n_list, n_gen = len(equations(from_list)), len(equations(from_gen))
    print('from a list            :', n_list, 'reaction(s)')
    print('from a generator       :', n_gen, 'reaction(s)')
    print('VERDICT                :',
          'MECHANISM SILENTLY LOST' if n_gen < n_list else 'ok -- the generator is honoured')
except Exception as exc:
    print('RAISED                 :', type(exc).__name__, str(exc)[:80])
    print('VERDICT                : loud failure, no silent loss')

print('\n--- DEFECT 2: core save, then core+edge save, over shared reaction objects ---')
core = rxn([H, OH], 1e12)
edge = rxn([H, OH], 2e12)
model = SimpleNamespace(core=SimpleNamespace(species=SPECIES, reactions=[core]),
                        edge=SimpleNamespace(species=[], reactions=[edge]),
                        output_species_list=[], output_reaction_list=[],
                        surface_site_density=None)
with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', '/tmp')) as d:
    save_chemkin(model, os.path.join(d, 'core.inp'), os.path.join(d, 'core-v.inp'),
                 save_edge_species=False)
    after_core = (core.duplicate, edge.duplicate)
    save_chemkin(model, os.path.join(d, 'edge.inp'), os.path.join(d, 'edge-v.inp'),
                 save_edge_species=True)
    after_edge = (core.duplicate, edge.duplicate)
print('after the core-only save  :', after_core)
print('after the core+edge save  :', after_edge)
from rmgpy.yaml_cantera2 import reaction_to_dict_list
entry = reaction_to_dict_list(core, SPECIES)[0]
print('core-only Cantera entry   : duplicate=%s  %s'
      % (entry.get('duplicate', False), entry['equation']))
print('VERDICT                   :',
      'LONE duplicate:true SHIPPED TO CANTERA' if entry.get('duplicate') else 'ok')
