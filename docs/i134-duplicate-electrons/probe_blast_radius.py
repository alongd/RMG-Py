#!/usr/bin/env python
# encoding: utf-8

"""
How far does the electron-blind duplicate check reach?

`docs/i123b-reaudit.md` Sec 7.2 established that RMG discards the lithium radiative
recombination as a duplicate of the electron-impact ionisation, and that forcing
`Reaction.electrons` to any value changes nothing. What it did NOT establish is the
*extent*: whether that blindness is a property of the one code path that two kinetics
libraries happen to meet, or a property of RMG's core-model duplicate check as a whole --
in which case every charged reaction RMG generates is exposed, including the ones families
produce.

That distinction decides the size of the ticket, so this probe measured it before anything
was changed. The answer was GENERAL: the blindness lives in the one predicate all four of
`check_for_existing_reaction`'s return sites share, two of which are the family branches.

The checks below now assert the REPAIRED behaviour, so this file is a standing regression
artifact rather than a one-shot measurement. What it looked like before the repair is
preserved verbatim in `evidence/probe-blast-radius-BEFORE.log`, where every verdict here is
inverted -- that log is the record of the defect, and the two should be read side by side.
Line numbers quoted from `rmgpy/rmg/model.py` in that log are one lower than the current
ones, because the repair added an import above them.

  STAGE 1  Which return site in `check_for_existing_reaction` fires for the real library
           pair -- recorded at runtime from the caller's line number, not inferred by
           reading the function. Before the repair it was `model.py:521` and the sink was
           discarded; now no site fires and both channels survive.

  STAGE 2  Where the electron count enters the comparison. `generate_reaction_id` is
           label-only by design and stays that way; `are_identical_species_references` is
           the single predicate all four return sites share, and is where the electrons
           are now compared.

  STAGE 3  BLAST RADIUS A: the FAMILY branches (`model.py:471` and `model.py:485`), driven
           with a real `KineticsFamily` and a real family-generated reaction out of the
           database -- non-dissociative electron attachment, `O2 + e- => O2-`. Three
           second reactions are offered against it: its true reverse (which SHOULD
           collapse, and is the control), the same-direction charge reversal, and the
           neutral transformation over the same heavy species.

  STAGE 4  BLAST RADIUS B: what is actually being confused with what. The electron
           PLACEMENT declaration -- the `(reactant_count, product_count)` pair the owner
           states, which the plasma reactor and both export writers already read -- says
           the lithium pair is not a reverse pair at all: reversing the ionisation's
           `(1, 2)` gives `(2, 1)`, three-body recombination, order 3; the radiative
           channel is `(1, 0)`, order 2. Their NET electron counts are nevertheless exactly
           equal and opposite, which is why no comparison built on the net scalar can tell
           them apart.

  STAGE 5  RMG's OWN canonical identity predicate, `Reaction.is_isomorphic`, had the same
           hole. It had already been taught about electrons -- but only about the NET
           scalar, so it collapsed the lithium pair too. The gap was one notion of identity
           short, not one function short, which is why both predicates now read the same
           per-side placement helper.

Run from the RMG-Py worktree root so `rmgrc` resolves:

    python docs/i134-duplicate-electrons/probe_blast_radius.py
"""

import inspect
import logging
import os

import rmgpy.rmg.model as model_module
from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.molecule import Molecule
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.electron_balance import get_placement_declaration
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species

logging.basicConfig(level=logging.CRITICAL)

IONISATION = 'PlasmaElectronImpactIonization'
RECOMBINATION = 'PlasmaRadiativeRecombination'
ATTACHMENT = 'Plasma_Electron_Attachment'

FAILURES = []
NOTES = []


def check(label, condition, detail=''):
    mark = 'PASS' if condition else 'FAIL'
    if not condition:
        FAILURES.append(label)
    print('  [{0}] {1}{2}'.format(mark, label, ' -- {0}'.format(detail) if detail else ''))


def note(text):
    NOTES.append(text)
    print('  [..] {0}'.format(text))


# ----------------------------------------------------------------------------------
# Instrumentation. `are_identical_species_references` is called from exactly four places,
# all inside `check_for_existing_reaction`. Recording the caller's line number therefore
# names the branch that fired, without editing the function under test.
# ----------------------------------------------------------------------------------

_REAL_IDENTICAL = model_module.are_identical_species_references
CALL_SITES = []


def _recording_identical(rxn1, rxn2):
    result = _REAL_IDENTICAL(rxn1, rxn2)
    CALL_SITES.append((inspect.currentframe().f_back.f_lineno, result))
    return result


def instrument():
    CALL_SITES.clear()
    model_module.are_identical_species_references = _recording_identical


def uninstrument():
    model_module.are_identical_species_references = _REAL_IDENTICAL


def deciding_site():
    """The line number of the call whose True verdict ended the search, or None."""
    for lineno, result in CALL_SITES:
        if result:
            return lineno
    return None


# ----------------------------------------------------------------------------------
# Loading
# ----------------------------------------------------------------------------------


def load_database():
    db = RMGDatabase()
    db.load(
        path=settings['database.directory'],
        thermo_libraries=['primaryThermoLibrary'],
        kinetics_families=[ATTACHMENT],
        reaction_libraries=[IONISATION, RECOMBINATION],
        kinetics_depositories=['training'],
    )
    return db


def sole_reaction(db, name):
    reactions = db.kinetics.libraries[name].get_library_reactions()
    assert len(reactions) == 1, '{0} has {1} entries, expected 1'.format(name, len(reactions))
    return reactions[0]


def fresh_model(db):
    m = CoreEdgeReactionModel()
    m.kinetics_database = db.kinetics
    return m


def add(m, rxn):
    """Return (reaction, is_new) exactly as add_reaction_library_to_edge does."""
    return m.make_new_reaction(rxn, generate_thermo=False, generate_kinetics=False)


def name_of(s):
    label = getattr(s, 'label', None)
    if label:
        return label
    try:
        return s.to_smiles()
    except AttributeError:
        return str(s)


def eq(rxn):
    return '{0} => {1}   [electrons={2:+d}]'.format(
        ' + '.join(name_of(s) for s in rxn.reactants),
        ' + '.join(name_of(s) for s in rxn.products),
        rxn.electrons)


def build_reverse(rxn):
    """Build the physical reverse of `rxn` the way `family.py:4249-4251` does: the two
    species lists swap and `electrons` is negated. Still a TemplateReaction of the same
    family, because the reverse of a family reaction still belongs to that family."""
    rev = TemplateReaction(
        reactants=[s for s in rxn.products],
        products=[s for s in rxn.reactants],
        electrons=-rxn.electrons,
        family=rxn.family,
        template=getattr(rxn, 'template', None),
        degeneracy=rxn.degeneracy,
    )
    rev.is_forward = False
    return rev


def build_same_direction_flipped(rxn):
    """Same reactants, same products, opposite electron stoichiometry. Physically this is
    a different reaction (it moves charge the other way across the same heavy species);
    it is the shape that reaches `model.py:471`, the same-direction branch."""
    flipped = TemplateReaction(
        reactants=[s for s in rxn.reactants],
        products=[s for s in rxn.products],
        electrons=-rxn.electrons,
        family=rxn.family,
        template=getattr(rxn, 'template', None),
        degeneracy=rxn.degeneracy,
    )
    flipped.is_forward = rxn.is_forward
    return flipped


# ----------------------------------------------------------------------------------


def main():
    print('=' * 84)
    print('database.directory = {0}'.format(settings['database.directory']))
    print('rmgpy package      = {0}'.format(os.path.dirname(model_module.__file__)))
    print('=' * 84)

    db = load_database()

    # ---------------------------------------------------------------- stage 1
    print('\nSTAGE 1 -- does the lithium sink reach the model? (recorded, not read)')
    m = fresh_model(db)
    instrument()
    _, ion_new = add(m, sole_reaction(db, IONISATION))
    ion_sites = list(CALL_SITES)
    CALL_SITES.clear()
    _, rec_new = add(m, sole_reaction(db, RECOMBINATION))
    rec_sites = list(CALL_SITES)
    site = deciding_site()
    uninstrument()

    print('    ionisation offered first : is_new={0}, identity calls at lines {1}'.format(
        ion_new, [s[0] for s in ion_sites]))
    print('    recombination offered next: is_new={0}, identity calls at lines {1}'.format(
        rec_new, [(s[0], s[1]) for s in rec_sites]))
    check('the source enters the model', ion_new is True, 'is_new={0}'.format(ion_new))
    check('the SINK enters the model too', rec_new is True, 'is_new={0}'.format(rec_new))
    check('no duplicate verdict fired against the sink',
          site is None,
          'deciding site model.py:{0} (it was model.py:521 before the repair)'.format(site))

    # ---------------------------------------------------------------- stage 2
    print('\nSTAGE 2 -- where does the electron count enter the comparison?')
    id_body = inspect.getsource(model_module.generate_reaction_id)
    check('generate_reaction_id stays label-only, by design',
          'electrons' not in id_body, '{0} occurrences'.format(id_body.count('electrons')))
    identical_body = inspect.getsource(_REAL_IDENTICAL)
    check('are_identical_species_references compares electron PLACEMENT',
          'get_electron_placement_counts' in identical_body,
          'the single predicate all four return sites share')
    check('...and does not fall back to the net scalar there',
          'reaction.electrons' not in identical_body
          and '.electrons ==' not in identical_body,
          'per-side counts only')

    ion = sole_reaction(db, IONISATION)
    rec = sole_reaction(db, RECOMBINATION)
    id_ion = model_module.generate_reaction_id(ion)
    id_rec = model_module.generate_reaction_id(rec)
    print('    generate_reaction_id(ionisation)    = {0}'.format(id_ion))
    print('    generate_reaction_id(recombination) = {0}'.format(id_rec))
    check('the two reaction IDs are exact reverses of one another',
          list(id_ion) == list(id_rec)[::-1],
          'so any reverse-match branch matches them')

    # ---------------------------------------------------------------- stage 3
    print('\nSTAGE 3 -- BLAST RADIUS A: are FAMILY-generated reactions affected?')
    fam = db.kinetics.families[ATTACHMENT]
    check('a real charged family loaded', fam is not None, ATTACHMENT)
    check('the family declares a non-zero electron stoichiometry',
          fam.electrons != 0, 'electrons={0:+d}'.format(fam.electrons))

    def fresh_attachment():
        """Regenerate the reaction from the family every time. `make_new_reaction`
        replaces a reaction's reactant/product lists with the model's own Species
        objects, and `are_identical_species_references` compares object identity, so a
        reaction reused across two models would carry the first model's objects into the
        second and report a verdict that is an artefact of the probe."""
        o2 = Species(molecule=[Molecule().from_adjacency_list("""
multiplicity 3
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
""")])
        o2.generate_resonance_structures()
        generated = fam.generate_reactions([o2.molecule])
        assert generated, 'the family generated nothing from O2'
        return generated[0]

    attach = fresh_attachment()
    check('the family generates a real reaction from a real reactant',
          attach is not None, 'from O2 via {0}'.format(ATTACHMENT))
    print('    generated by {0}: {1}'.format(ATTACHMENT, eq(attach)))
    check('the generated reaction carries the family electron declaration',
          attach.electrons == fam.electrons,
          'Reaction.electrons={0:+d}'.format(attach.electrons))

    def offer_second(second):
        """Put a real family reaction in a fresh model, then offer `second`. Returns
        (is_new, deciding line number)."""
        m2 = fresh_model(db)
        instrument()
        _, first_new = add(m2, fresh_attachment())
        assert first_new is True, 'the attachment did not enter the model'
        CALL_SITES.clear()
        _, second_new = add(m2, second)
        site = deciding_site()
        uninstrument()
        return second_new, site

    # -- 3a: the CONTROL. The true reverse of the attachment is the detachment; a
    #    reaction and its reverse are the same reaction, and collapsing them is what the
    #    duplicate check is FOR. This is here to be kept, not to be fixed.
    detach = build_reverse(fresh_attachment())
    print('    3a CONTROL: true reverse (family.py:4249-4251 construction): {0}'.format(eq(detach)))
    detach_new, site_3a = offer_second(detach)
    check('CONTROL: a genuine reverse pair collapses (correct behaviour)',
          detach_new is False,
          'is_new={0}, deciding site model.py:{1}'.format(detach_new, site_3a))
    note('the family REVERSE branch is reached, at rmgpy/rmg/model.py:{0}'.format(site_3a))

    # -- 3b: same heavy species, same direction, opposite charge transfer. `A + e- => B`
    #    and `A => B + e-` are different reactions: different molecularity, different
    #    charge balance. Nothing about them is a reverse pair.
    flipped = build_same_direction_flipped(fresh_attachment())
    print('    3b same direction, opposite charge transfer: {0}'.format(eq(flipped)))
    flipped_new, site_3b = offer_second(flipped)
    check('a same-direction electron-flipped reaction now SURVIVES',
          flipped_new is True,
          'is_new={0}, deciding site model.py:{1} (it was dropped at model.py:471 '
          'before the repair)'.format(flipped_new, site_3b))

    # -- 3c: the fully general shape, needing no exotic data at all. A NEUTRAL reaction
    #    over the same heavy species as a charged one. `A => B` (electrons=0) and
    #    `A + e- => B` (electrons=-1) share every field the duplicate check reads.
    neutral = build_same_direction_flipped(fresh_attachment())
    neutral.electrons = 0
    print('    3c the same heavy transformation with NO electron: {0}'.format(eq(neutral)))
    neutral_new, site_3c = offer_second(neutral)
    check('a NEUTRAL reaction over the same heavy species now SURVIVES',
          neutral_new is True,
          'is_new={0}, deciding site model.py:{1} (it was dropped at model.py:471 '
          'before the repair)'.format(neutral_new, site_3c))
    note('the neutral/charged collision needed no plasma library at all, which is what '
         'made this a general defect rather than a plasma one')

    # ---------------------------------------------------------------- stage 4
    print('\nSTAGE 4 -- BLAST RADIUS B: what exactly is being confused with what?')
    print('    Electron PLACEMENT, from rmgpy.electron_balance.get_placement_declaration --')
    print('    the (reactant_count, product_count) pair the owner declares, which is what')
    print('    the plasma reactor and both export writers already read:')
    for label, rxn in (('ionisation   ', ion), ('recombination', rec),
                       ('attachment   ', attach)):
        print('        {0}  family={1:32s} placement={2}  net electrons={3:+d}'.format(
            label, str(getattr(rxn, 'family', None)), get_placement_declaration(rxn),
            rxn.electrons))

    ion_placement = get_placement_declaration(ion)
    rec_placement = get_placement_declaration(rec)
    check('the lithium pair is NOT a reverse pair: reversing the ionisation placement '
          'does not give the recombination placement',
          tuple(reversed(ion_placement)) != rec_placement,
          '{0} reversed is {1}, but recombination is {2} -- the reverse of the ionisation '
          'is three-body recombination (order 3), not the radiative channel (order 2)'.format(
              ion_placement, tuple(reversed(ion_placement)), rec_placement))
    check('...yet their NET electron counts are exactly equal and opposite',
          ion.electrons == -rec.electrons,
          '{0:+d} vs {1:+d} -- so no comparison built on the net scalar alone can '
          'separate them'.format(ion.electrons, rec.electrons))

    instances = [
        ('library / lithium', 'ionisation vs radiative recombination',
         rec_new is False, False),
        ('family / oxygen', 'attachment vs its true reverse', detach_new is False, True),
        ('family / oxygen', 'attachment vs same-direction flip', flipped_new is False, False),
        ('family / oxygen', 'attachment vs the neutral transformation',
         neutral_new is False, False),
    ]
    print('')
    for path, what, collapsed, should_collapse in instances:
        print('    {0:20s} {1:42s} collapsed={2!s:5s} correct={3}'.format(
            path, what, collapsed, should_collapse))
    check('exactly the genuine reverse pair collapses, and nothing else',
          all(collapsed is should for _, _, collapsed, should in instances),
          'the repair reaches the family branches, not just the two lithium libraries')

    # ---------------------------------------------------------------- stage 5
    print("\nSTAGE 5 -- RMG's own canonical identity predicate, repaired the same way")
    # Fresh, un-added objects: `make_new_reaction` swaps a reaction's lists for the
    # model's own Species, and comparing a mutated reaction against an unmutated one
    # would compare Species against Molecule.
    attach_iso = fresh_attachment()
    pairs = [
        ('ionisation vs recombination', sole_reaction(db, IONISATION),
         sole_reaction(db, RECOMBINATION), False),
        ('attachment vs its true reverse', attach_iso, build_reverse(fresh_attachment()), True),
        ('attachment vs same-direction flip', attach_iso,
         build_same_direction_flipped(fresh_attachment()), False),
    ]
    for label, a, b, should_be_same in pairs:
        iso = a.is_isomorphic(b, either_direction=True)
        print('    is_isomorphic({0:36s}) = {1!s:5s}  (correct answer: {2})'.format(
            label, iso, should_be_same))
        check('is_isomorphic is correct for: {0}'.format(label), iso is should_be_same)
    note('is_isomorphic used to compare the NET scalar, never the per-side placement, and '
         'so collapsed the lithium pair for the same reason model.py did -- see the same '
         'three lines in evidence/probe-blast-radius-BEFORE.log. Both predicates now read '
         'the same placement helper, which is what keeps ReactionModel.merge, '
         'KineticsLibrary.check_for_duplicates and the pressure-dependence code consistent '
         'with the model builder.')

    print('\n' + '=' * 84)
    for n in NOTES:
        print('NOTE: {0}'.format(n))
    if FAILURES:
        print('BLAST-RADIUS PROBE: {0} CHECK(S) FAILED'.format(len(FAILURES)))
        for f in FAILURES:
            print('   - {0}'.format(f))
    else:
        print('BLAST-RADIUS PROBE: all checks passed.')
    print('=' * 84)
    return 1 if FAILURES else 0


if __name__ == '__main__':
    raise SystemExit(main())
