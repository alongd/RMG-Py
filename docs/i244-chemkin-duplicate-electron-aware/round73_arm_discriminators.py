"""
Four values that say which engine is actually loaded, printed from the built
extension rather than read off the source.

Run inside each RED arm. Exactly one line must differ from the repaired engine's,
and it must be the arm's own. A build that silently did not take, or a restore
that reverted more than it was meant to, shows up here as the wrong line moving --
which is the failure this file exists to catch, having already happened once.
"""
import logging

logging.getLogger().setLevel(logging.ERROR)

import rmgpy.chemkin as cm
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial


def spc(label, index, smiles):
    s = Species(label=label, molecule=[Molecule(smiles=smiles)])
    s.index = index
    coeffs = [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, -11.7246]
    s.thermo = NASA(
        polynomials=[NASAPolynomial(coeffs=coeffs, Tmin=(200, 'K'), Tmax=(1000, 'K')),
                     NASAPolynomial(coeffs=coeffs, Tmin=(1000, 'K'), Tmax=(6000, 'K'))],
        Tmin=(200, 'K'), Tmax=(6000, 'K'))
    return s


H, OH, H2, O = (spc('H', 1, '[H]'), spc('OH', 2, '[OH]'),
                spc('H2', 3, '[H][H]'), spc('O', 4, '[O]'))
SPECIES = [H, OH, H2, O]


def rxn(reactants, A, dup=False):
    r = LibraryReaction(reactants=list(reactants), products=[H2, O],
                        kinetics=Arrhenius(A=(A, 'cm^3/(mol*s)'), n=0, Ea=(0, 'kcal/mol')),
                        library='SomeLibrary', reversible=False)
    r.duplicate = dup
    return r


# D1 -- a permuted pair keys alike (repaired: True)
a, b = rxn([H, OH], 1e12), rxn([OH, H], 2e12)
print('D1 permuted pair keys alike        :',
      cm.chemkin_duplicate_group_key(a) == cm.chemkin_duplicate_group_key(b))

# D2 -- rendering leaves the flags alone (repaired: True)
c, d = rxn([H, OH], 1e12, dup=True), rxn([OH, H], 2e12, dup=True)
lonely = rxn([H2, O], 3e12, dup=True)   # two reactants, so the A-factor units match
cm.render_chemkin_file(SPECIES, [c, d, lonely], verbose=False)
print('D2 render leaves flags untouched   :',
      (c.duplicate, d.duplicate, lonely.duplicate) == (True, True, True))

# D3 -- a generator survives the render (repaired: True)
n_list = cm.render_chemkin_file(SPECIES, [rxn([H, OH], 1e12)], verbose=False).count('=>')
n_gen = cm.render_chemkin_file(SPECIES, (x for x in [rxn([H, OH], 1e12)]),
                               verbose=False).count('=>')
print('D3 generator survives the render   :', n_gen == n_list and n_list > 0)

# D4 -- the authority overrides the flags it is handed (repaired: True)
wrong = [rxn([H, OH], 1e12, dup=True)]            # lone, flagged: answer must be False
pair = [rxn([H, OH], 1e12), rxn([H, OH], 2e12)]   # unflagged pair: answer must be True
print('D4 authority overrides its input   :',
      cm.chemkin_duplicate_flags(wrong) == [False]
      and cm.chemkin_duplicate_flags(pair) == [True, True])
