#!/usr/bin/env python
"""I-148 repair stress: monotonicity (8.1) and guard-not-weakened (8.2).

Prints the resolved rmgpy/database path at the head, then:
  8.1 -- for a spread of reactions, compares get_placement_owner (I-148, family
         then library) against the PRE-I-148 behaviour (family only), asserting
         (a) every reaction that resolved before resolves IDENTICALLY, and
         (b) only reactions whose current attribution (family) names no
             declaration can NEWLY resolve.
  8.2 -- feeds a plasma-shaped reaction with its family stripped to the reactor
         resolver resolve_electron_placement and asserts it RAISES its named
         ElectronPlacementError rather than net-deriving.
Object-level; does not use the pytest scratch path.
"""
import os, sys
import rmgpy
from rmgpy import settings
print('[HEAD] rmgpy.__file__     =', rmgpy.__file__)
print('[HEAD] database.directory =', settings['database.directory'])

from rmgpy.electron_balance import get_placement_owner, get_placement_declaration
from rmgpy.electron_placement import (FAMILY_ELECTRON_PLACEMENT,
                                      resolve_electron_placement)
from rmgpy.exceptions import ElectronPlacementError
from rmgpy.reaction import Reaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.species import Species
from rmgpy.molecule import Molecule

print('\n[HEAD] FAMILY_ELECTRON_PLACEMENT =', dict(FAMILY_ELECTRON_PLACEMENT))


def pre_i148_owner(reaction):
    """The behaviour before I-148: only Reaction.family was consulted."""
    fam = getattr(reaction, 'family', None)
    return fam if (fam and fam in FAMILY_ELECTRON_PLACEMENT) else None


# ---- build a spread of reactions -------------------------------------------
Li = Species(label='Li', molecule=[Molecule(smiles='[Li]')])
Lip = Species(label='Li+', molecule=[Molecule(smiles='[Li+]')])
H = Species(label='H', molecule=[Molecule(smiles='[H]')])
H2 = Species(label='H2', molecule=[Molecule(smiles='[H][H]')])

cases = []
# (label, reaction, expected_pre, expected_post, note)
# ordinary non-plasma template reaction: family not declared
r_ord = TemplateReaction(reactants=[H, H2], products=[H2, H], family='H_Abstraction')
cases.append(('ordinary H_Abstraction (family-generated)', r_ord, None, None,
              'undeclared family; untouched by I-148'))
# plasma family-generated reaction (declared family)
r_famplasma = TemplateReaction(reactants=[Li], products=[Lip], family='PlasmaElectronImpactIonization')
cases.append(('plasma family PlasmaElectronImpactIonization', r_famplasma,
              'PlasmaElectronImpactIonization', 'PlasmaElectronImpactIonization',
              'declared family; family branch, identical before/after'))
# library reaction: LibraryReaction.__init__ sets family=library
r_lib = LibraryReaction(reactants=[Li], products=[Lip], library='PlasmaElectronImpactIonization')
cases.append(('library PlasmaElectronImpactIonization (fresh load)', r_lib,
              'PlasmaElectronImpactIonization', 'PlasmaElectronImpactIonization',
              'library sets family=library; resolves via family branch'))
# seed-RENAMED reaction: get_library_reactions overwrites family='seed', but
# preserves the original owner in .library. THIS is the I-148 case.
r_seed = LibraryReaction(reactants=[Li], products=[Lip], library='PlasmaElectronImpactIonization')
r_seed.family = 'seed'                     # container rename the seed reload performs
r_seed.library = 'PlasmaElectronImpactIonization'  # preserved provenance
cases.append(('seed-renamed (family="seed", library preserved)', r_seed,
              None, 'PlasmaElectronImpactIonization',
              'family names nothing -> newly resolves via library provenance (I-148 fix)'))
# seed-renamed with NO usable provenance at all
r_lost = LibraryReaction(reactants=[Li], products=[Lip], library='seed')
r_lost.family = 'seed'
r_lost.library = 'seed'
cases.append(('seed-renamed, no provenance (family=library="seed")', r_lost,
              None, None, 'no attribution names a declaration -> None, both'))

print('\n==== 8.1 MONOTONICITY ====')
mono_ok = True
newly = []
for label, rxn, exp_pre, exp_post, note in cases:
    pre = pre_i148_owner(rxn)
    post = get_placement_owner(rxn)
    # (a) resolved-before must resolve identically
    resolved_before_ok = (pre is None) or (post == pre)
    # (b) newly-resolving must have had an undeclared family
    newly_resolved = (pre is None and post is not None)
    if newly_resolved:
        newly.append(label)
        newly_ok = (pre_i148_owner(rxn) is None)  # its family named no declaration
    else:
        newly_ok = True
    ok = (pre == exp_pre and post == exp_post and resolved_before_ok and newly_ok)
    mono_ok = mono_ok and ok
    print(f'  [{"OK" if ok else "FAIL"}] {label}')
    print(f'         pre(family-only)={pre!r}  post(I-148)={post!r}  '
          f'(expected pre={exp_pre!r} post={exp_post!r})  {note}')
print(f'  newly-resolving reactions: {newly}')
print(f'  MONOTONICITY: {"HOLDS" if mono_ok else "VIOLATED"} '
      f'(every resolved-before resolves identically; only undeclared-family reactions newly resolve)')

print('\n==== 8.2 GUARD NOT WEAKENED ====')
# A plasma-shaped reaction with NO family attribution must RAISE, never net-derive.
# A plain Reaction has no `family` attribute at all: getattr(...,'family',None) is
# None, which is exactly the no-provenance condition the guard must refuse.
r_nofam = Reaction(reactants=[Li], products=[Lip])
r_nofam.electrons = 1            # net: one electron produced (ionisation-shaped)
assert getattr(r_nofam, 'family', None) is None, 'plain Reaction should carry no family'
elec = Species(label='e-', molecule=[Molecule().from_adjacency_list('1 e u1 p0 c-1')])
species_list = [Li, Lip, elec]
try:
    view = resolve_electron_placement(r_nofam, species_list)
    print(f'  [FAIL] resolver returned a view instead of raising: {view}')
    guard_ok = False
except ElectronPlacementError as e:
    msg = str(e)
    names_reaction = 'family' in msg and ('cannot be inferred' in msg or 'no family' in msg)
    guard_ok = True
    print('  [OK] resolve_electron_placement RAISED ElectronPlacementError, no net-derive fallback')
    print('       message:', msg.replace('\n', ' ')[:200])
    print(f'       names the failure (family-declared, not net-derived): {names_reaction}')

print('\n==== RESULT ====')
print('8.1 monotonicity:', 'HOLDS' if mono_ok else 'VIOLATED')
print('8.2 guard:', 'INTACT (raises, no net fallback)' if guard_ok else 'WEAKENED')
sys.exit(0 if (mono_ok and guard_ok) else 1)
