#!/usr/bin/env python
# encoding: utf-8

"""
The three consumers of the repaired identity predicate that were fixed by construction
and never exercised.

`docs/i134-duplicate-electrons.md` Sec 9 names them and marks all three UNKNOWN:

  1. `ReactionModel.merge`  (`rmgpy/rmg/model.py:129`, `either_direction=True`) -- reached
     by the repair through `Reaction.is_isomorphic`, but no `mergeModels.py` was ever run
     over two plasma models. If merge were still electron-blind, the same drop this whole
     line of work exists to fix would reappear there silently.
  2. `rmgpy/rmg/pdep.py:927,939` -- the same predicate, consumed when a pressure-dependent
     network reaction is checked against an existing core/edge library reaction.
  3. `Reaction.is_isomorphic(check_template_rxn_products=True)`
     (`rmgpy/reaction.py:707-727`) -- the products-only generation shortcut, which still
     compares the NET `electrons` scalar and negates it on `is_forward`. The question here
     is scoping, not repair: is that shortcut reachable in a plasma generation path?

Stages:

  STAGE 1  MERGE, in memory. Two `ReactionModel`s, one carrying the electron-impact
           ionisation and one the radiative recombination, merged with the real
           `ReactionModel.merge`. Both channels must survive. Two controls run beside it:
           a genuine duplicate must still collapse, and a genuine reverse pair must still
           collapse -- otherwise "both survived" would only mean merge had stopped
           deduplicating at all.

  STAGE 2  PDEP. Reachability first: whether any shipped plasma deck turns pressure
           dependence on, and whether anything in `rmgpy/rmg/pdep.py` ever gives a
           `PDepReaction` a non-zero electron count. Then the comparison pdep.py actually
           makes, driven directly, so the verdict is measured rather than assumed latent.

  STAGE 3  THE GENERATION SHORTCUT. Whether `check_template_rxn_products=True` is reached
           when a real plasma family generates reactions out of the database, and -- if it
           is -- what `is_forward` is on the reactions that reach it, since that flag is
           what the doubted negation keys on.

Run from the RMG-Py worktree root so `rmgrc` resolves:

    python docs/i123c-audit/probe_identity_consumers.py
"""

import copy
import os
import sys

import rmgpy.data.kinetics.family as family_module
from rmgpy import settings
from rmgpy import electron_placement
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.electron_balance import get_electron_placement_counts
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import ReactionModel
from rmgpy.rmg.pdep import PDepReaction
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial

IONISATION = 'PlasmaElectronImpactIonization'
RECOMBINATION = 'PlasmaRadiativeRecombination'
ATTACHMENT = 'Plasma_Electron_Attachment'

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


def flat_thermo():
    """Constant-Cp NASA polynomial. `ReactionModel.merge` compares the thermo of any two
    species it judges common, so every species here needs one; nothing below depends on
    these being the real numbers."""
    coeffs = [2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    return NASA(
        polynomials=[NASAPolynomial(coeffs=coeffs, Tmin=(100, 'K'), Tmax=(5000, 'K')),
                     NASAPolynomial(coeffs=coeffs, Tmin=(5000, 'K'), Tmax=(20000, 'K'))],
        Tmin=(100, 'K'), Tmax=(20000, 'K'))


def header():
    print('=' * 78)
    print('rmgpy              = {0}'.format(os.path.abspath(
        os.path.dirname(os.path.dirname(electron_placement.__file__)))))
    import rmgpy.molecule.molecule as _m
    print('compiled module    = {0}'.format(_m.__file__))
    print('database.directory = {0}'.format(settings['database.directory']))
    print('resolved realpath  = {0}'.format(os.path.realpath(settings['database.directory'])))
    print('exists             = {0}'.format(os.path.isdir(settings['database.directory'])))
    print('=' * 78)


def load_channels():
    db = KineticsDatabase()
    db.load_libraries(os.path.join(settings['database.directory'], 'kinetics', 'libraries'),
                      libraries=[IONISATION, RECOMBINATION])
    ionisation = db.libraries[IONISATION].get_library_reactions()[0]
    recombination = db.libraries[RECOMBINATION].get_library_reactions()[0]
    return ionisation, recombination


# ------------------------------------------------------------------ stage 1: merge
def stage1_merge(ionisation, recombination):
    print('\nSTAGE 1 -- ReactionModel.merge over two models, one channel each')

    for rxn in (ionisation, recombination):
        for spc in list(rxn.reactants) + list(rxn.products):
            if spc.thermo is None:
                spc.thermo = flat_thermo()

    note('ionisation    : {0}  electrons={1}  placement={2}  owner={3}'.format(
        ionisation, ionisation.electrons,
        get_electron_placement_counts(ionisation), ionisation.family))
    note('recombination : {0}  electrons={1}  placement={2}  owner={3}'.format(
        recombination, recombination.electrons,
        get_electron_placement_counts(recombination), recombination.family))

    source_model = ReactionModel()
    source_model.species = list(ionisation.reactants) + list(ionisation.products)
    source_model.reactions = [ionisation]

    sink_model = ReactionModel()
    sink_model.species = list(recombination.reactants) + list(recombination.products)
    sink_model.reactions = [recombination]

    merged = source_model.merge(sink_model)
    print('    merged model: {0} species, {1} reactions'.format(
        len(merged.species), len(merged.reactions)))
    for rxn in merged.reactions:
        print('      - {0}   electrons={1}  kinetics={2}'.format(
            rxn, rxn.electrons, type(rxn.kinetics).__name__))

    check('merge keeps BOTH channels', len(merged.reactions) == 2,
          '{0} reaction(s) in the merged model'.format(len(merged.reactions)))
    kinds = sorted(type(r.kinetics).__name__ for r in merged.reactions)
    check('the merged model carries both rate laws',
          kinds == ['BadnellRRArrhenius', 'VoronovEIArrhenius'], str(kinds))
    check('the cation is consumed somewhere in the merged model',
          any(any(s.is_isomorphic(ionisation.products[0]) for s in r.reactants)
              for r in merged.reactions),
          'a loss channel for the cation survived the merge')

    # ---- control A: a genuine duplicate must still collapse
    twin = copy.deepcopy(ionisation)
    for spc in list(twin.reactants) + list(twin.products):
        if spc.thermo is None:
            spc.thermo = flat_thermo()
    dup_a = ReactionModel()
    dup_a.species = list(ionisation.reactants) + list(ionisation.products)
    dup_a.reactions = [ionisation]
    dup_b = ReactionModel()
    dup_b.species = list(twin.reactants) + list(twin.products)
    dup_b.reactions = [twin]
    merged_dup = dup_a.merge(dup_b)
    check('CONTROL A: merging a model with itself still deduplicates',
          len(merged_dup.reactions) == 1,
          '{0} reaction(s)'.format(len(merged_dup.reactions)))

    # ---- control B: a genuine reverse pair must still collapse
    true_reverse = copy.deepcopy(ionisation)
    reactants, products = list(true_reverse.reactants), list(true_reverse.products)
    true_reverse.reactants = products
    true_reverse.products = reactants
    true_reverse.electrons = -ionisation.electrons
    for spc in list(true_reverse.reactants) + list(true_reverse.products):
        if spc.thermo is None:
            spc.thermo = flat_thermo()
    rev_a = ReactionModel()
    rev_a.species = list(ionisation.reactants) + list(ionisation.products)
    rev_a.reactions = [ionisation]
    rev_b = ReactionModel()
    rev_b.species = list(true_reverse.reactants) + list(true_reverse.products)
    rev_b.reactions = [true_reverse]
    merged_rev = rev_a.merge(rev_b)
    note('true reverse of the ionisation: {0}  electrons={1}  placement={2}'.format(
        true_reverse, true_reverse.electrons, get_electron_placement_counts(true_reverse)))
    check('CONTROL B: a genuine reverse pair still collapses in merge',
          len(merged_rev.reactions) == 1,
          '{0} reaction(s)'.format(len(merged_rev.reactions)))


# ------------------------------------------------------------------- stage 2: pdep
def stage2_pdep(ionisation, recombination):
    print('\nSTAGE 2 -- the pressure-dependence consumer')

    # -- reachability, part 1: does any shipped plasma deck ask for pressure dependence?
    root = os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
    decks = []
    for base, _dirs, files in os.walk(os.path.join(root, 'docs')):
        for fname in files:
            if fname.startswith('input') and fname.endswith('.py'):
                decks.append(os.path.join(base, fname))
    plasma_decks = []
    pdep_decks = []
    for deck in sorted(decks):
        with open(deck) as handle:
            text = handle.read()
        if 'plasmaReactor(' in text:
            plasma_decks.append(deck)
            if 'pressureDependence(' in text:
                pdep_decks.append(deck)
    note('plasma decks found: {0}'.format(len(plasma_decks)))
    for deck in plasma_decks:
        note('    {0}'.format(os.path.relpath(deck, root)))
    check('no shipped plasma deck turns pressure dependence on',
          not pdep_decks, str(pdep_decks))

    # -- reachability, part 2: can a PDepReaction ever carry an electron count?
    import rmgpy.rmg.pdep as pdep_module
    with open(pdep_module.__file__.replace('.pyc', '.py')) as handle:
        pdep_source = handle.read()
    assignments = [line.strip() for line in pdep_source.splitlines()
                   if 'electrons' in line]
    note('every mention of `electrons` in rmgpy/rmg/pdep.py:')
    for line in assignments:
        note('    {0}'.format(line))
    probe_pdep = PDepReaction(reactants=list(ionisation.reactants),
                             products=list(ionisation.products))
    check('a PDepReaction is born with electrons = 0',
          probe_pdep.electrons == 0, str(probe_pdep.electrons))
    check('and therefore with placement counts (0, 0)',
          get_electron_placement_counts(probe_pdep) == (0, 0),
          str(get_electron_placement_counts(probe_pdep)))

    # -- the comparison pdep.py actually makes, driven directly
    print('    the pdep.py:927/939 comparison, driven directly:')
    for label, other in (('ionisation (electrons=+1)', ionisation),
                         ('recombination (electrons=-1)', recombination)):
        verdict = other.is_isomorphic(probe_pdep, either_direction=True)
        note('    LibraryReaction {0}  vs  neutral PDepReaction over the same heavy '
             'species -> is_isomorphic(either_direction=True) = {1}'.format(label, verdict))
        check('the {0} does not match a neutral network reaction'.format(label),
              verdict is False, str(verdict))

    # -- and what a charged PDepReaction would do, if anything ever made one
    charged_pdep = PDepReaction(reactants=list(ionisation.reactants),
                                products=list(ionisation.products),
                                electrons=ionisation.electrons)
    note('a hand-built PDepReaction with electrons=+1 has placement counts {0} '
         '(no owner => the net-derived fallback)'.format(
             get_electron_placement_counts(charged_pdep)))
    note('    ionisation vs that charged PDepReaction -> {0}'.format(
        ionisation.is_isomorphic(charged_pdep, either_direction=True)))
    note('    the ionisation library DECLARES (1, 2); the ownerless PDepReaction derives '
         '(0, 1) from the net scalar, so they differ -- this is the declared/undeclared '
         'seam, not a defect of the repair')


# --------------------------------------------------------- stage 3: the shortcut
import rmgpy.data.kinetics.database as kinetics_db_module

_REAL_FIND_DEGENERATE = family_module.find_degenerate_reactions
SHORTCUT_CALLS = []


def _recording_find_degenerate(rxn_list, *args, **kwargs):
    record = {
        'n_offered': len(rxn_list),
        'charged': [r.electrons for r in rxn_list if getattr(r, 'electrons', 0)],
        'is_forward': [getattr(r, 'is_forward', None) for r in rxn_list],
        'comparisons_possible': len(rxn_list) > 1,
    }
    result = _REAL_FIND_DEGENERATE(rxn_list, *args, **kwargs)
    record['n_returned'] = len(result)
    SHORTCUT_CALLS.append(record)
    return result


def stage3_shortcut():
    print('\nSTAGE 3 -- is the products-only generation shortcut reachable in a plasma path?')

    print('    static reachability: who passes check_template_rxn_products=True')
    root = os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
    for rel in ('rmgpy/data/kinetics/common.py', 'rmgpy/data/kinetics/family.py'):
        with open(os.path.join(root, rel)) as handle:
            for lineno, line in enumerate(handle, 1):
                if 'check_template_rxn_products=True' in line:
                    note('    {0}:{1}  {2}'.format(rel, lineno, line.strip()))
    note('    both of those sit inside find_degenerate_reactions '
         '(rmgpy/data/kinetics/common.py) and its family-side callers, which is the only '
         'production route to the shortcut')

    database = RMGDatabase()
    database.load(
        path=settings['database.directory'],
        thermo_libraries=['primaryThermoLibrary'],
        reaction_libraries=[],
        seed_mechanisms=[],
        kinetics_families=[ATTACHMENT],
        kinetics_depositories=['training'],
        depository=False,
    )

    o2 = Species(molecule=[Molecule().from_adjacency_list(
        'multiplicity 3\n1 O u1 p2 c0 {2,S}\n2 O u1 p2 c0 {1,S}\n')])
    o2.generate_resonance_structures()

    # find_degenerate_reactions is imported by value into two module namespaces; the
    # production route RMG itself takes is the one in rmgpy/data/kinetics/database.py.
    family_module.find_degenerate_reactions = _recording_find_degenerate
    kinetics_db_module.find_degenerate_reactions = _recording_find_degenerate
    try:
        SHORTCUT_CALLS.clear()
        generated = database.kinetics.generate_reactions_from_families(
            [o2], only_families=[ATTACHMENT])
    finally:
        family_module.find_degenerate_reactions = _REAL_FIND_DEGENERATE
        kinetics_db_module.find_degenerate_reactions = _REAL_FIND_DEGENERATE

    print('    {0} reaction(s) generated by {1} from triplet O2, through '
          'generate_reactions_from_families -- the route RMG itself takes'.format(
              len(generated), ATTACHMENT))
    for rxn in generated:
        print('      - {0}   electrons={1}  is_forward={2}  degeneracy={3}  type={4}'.format(
            rxn, getattr(rxn, 'electrons', None), getattr(rxn, 'is_forward', None),
            getattr(rxn, 'degeneracy', None), type(rxn).__name__))

    check('the family generated at least one charged reaction',
          any(getattr(r, 'electrons', 0) for r in generated),
          str([getattr(r, 'electrons', 0) for r in generated]))

    print('    find_degenerate_reactions calls recorded during that generation:')
    for i, record in enumerate(SHORTCUT_CALLS):
        note('    call {0}: offered={1} returned={2} charged_electrons={3} '
             'is_forward={4} comparisons_possible={5}'.format(
                 i, record['n_offered'], record['n_returned'], record['charged'],
                 record['is_forward'], record['comparisons_possible']))

    reached = bool(SHORTCUT_CALLS)
    check('the degeneracy path -- the shortcut\'s only production caller -- IS reached '
          'during plasma family generation', reached, str(len(SHORTCUT_CALLS)))

    charged_reached = any(record['charged'] for record in SHORTCUT_CALLS)
    check('and it is reached carrying CHARGED reactions', charged_reached)

    compared = sum(1 for r in SHORTCUT_CALLS if r['comparisons_possible'])
    check('and with more than one reaction offered, so a comparison actually happens',
          compared > 0, '{0} of {1} calls offered >1 reaction'.format(
              compared, len(SHORTCUT_CALLS)))

    forwards = set()
    for record in SHORTCUT_CALLS:
        forwards.update(record['is_forward'])
    note('is_forward values seen on reactions entering the shortcut: {0}'.format(
        sorted(str(f) for f in forwards)))
    note('the doubted negation (`-self.electrons` when `is_forward` is falsy) is only '
         'taken when is_forward is falsy on a reaction reaching this path')
    check('every reaction reaching the shortcut on this path is is_forward=True, so the '
          'doubted negation is never taken here', forwards == {True}, str(forwards))

    # The shortcut's own behaviour, driven directly on a charged pair, so the scoping
    # answer is not left resting on whether the family happened to generate two.
    print('    the shortcut driven directly on a charged pair:')
    charged = [r for r in generated if getattr(r, 'electrons', 0)]
    if charged:
        first = charged[0]
        twin = copy.deepcopy(first)
        flipped = copy.deepcopy(first)
        flipped.electrons = -first.electrons
        note('    same electrons  -> {0}'.format(
            first.is_isomorphic(twin, check_identical=False, strict=False,
                                check_template_rxn_products=True)))
        note('    opposite electrons -> {0}'.format(
            first.is_isomorphic(flipped, check_identical=False, strict=False,
                                check_template_rxn_products=True)))
        check('the shortcut separates a charged reaction from its charge-reversed twin',
              not first.is_isomorphic(flipped, check_identical=False, strict=False,
                                      check_template_rxn_products=True))


def main():
    header()
    ionisation, recombination = load_channels()
    stage1_merge(ionisation, recombination)
    stage2_pdep(ionisation, recombination)
    stage3_shortcut()

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
