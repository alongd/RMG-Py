#!/usr/bin/env python3
"""
Round 111 probe -- the repair is installed at the call site, not in the thing repaired.

    MPLCONFIGDIR=/tmp/claude-1000/mpl PYTHONPATH=$PWD \\
        ~/anaconda3/envs/rmg_env/bin/python \\
        docs/i221-saturation/probes/round111_second_entry_point_probe.py

Round 110 completed the five lossy reducers and installed them inside `copy()`, by giving
one `pickle.Pickler` subclass a ``dispatch_table``. A dispatch table belongs to the pickler
that carries it, so every *other* transport in RMG still gets the lossy reducers. The
sharpest of those is the one production uses on every parallel run: `Pool.map`
(``rmgpy/rmg/react.py:69``) serialises the generated reactions with `multiprocessing`'s own
`ForkingPickler`, which has never heard of round 110.

Beside it, two more repairs that stopped at the object in front of them:

* base `Reaction.copy()` builds one ``id(original) -> copy`` map over reactants and
  products together, so a `Species` that is on **both** sides -- which is every
  electron-impact reaction this campaign generates, ``Li + e- => Li+ + e- + e-`` -- has its
  reactant-side `pairs` member overwritten by the product copy;
* `_warn_unattributable()` enumerates **loaded** families only, so a library rate whose
  ``family:`` line is missing passes in silence when the quarantined family is on disk and
  not loaded -- which is the ordinary foreign-seed case, and the one the gate exists for.

And the two subclasses nobody carried: `DepositoryReaction` and `PDepReaction` inherit a
`copy()` that returns a base `Reaction`.

``REPRODUCED`` means the defect is present. Controls must hold at both tips. Exit 1 while
any finding is present, 2 if a control breaks.
"""

import logging
import os
import pickle
import shutil
import sys
import tempfile
from multiprocessing.reduction import ForkingPickler

logging.disable(logging.INFO)

from rmgpy import settings
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.depository import DepositoryReaction, KineticsDepository
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.data.kinetics.quarantine import (QUARANTINE_FILENAME, check_quarantine,
                                            load_family_quarantine)
from rmgpy.electron_placement import resolve_electron_placement
from rmgpy.molecule import Molecule
from rmgpy.molecule.fragment import CuttingLabel, Fragment
from rmgpy.reaction import Reaction
from rmgpy.rmg.pdep import PDepReaction
from rmgpy.species import Species

FAMILY = 'Cation_R_Recombination'
IONISATION = 'Plasma_Electron_Impact_Ionization'

results = []
controls = []
notes = []


def record(name, reproduced, detail):
    results.append((name, reproduced, detail))
    print("\n{0}  {1}".format("REPRODUCED " if reproduced else "absent     ", name))
    print("    " + detail.replace("\n", "\n    "))


def control(name, ok, detail):
    controls.append((name, ok, detail))
    print("\n{0}  control: {1}".format("holds      " if ok else "BROKEN     ", name))
    print("    " + detail.replace("\n", "\n    "))


def note(name, detail):
    notes.append((name, detail))
    print("\nnote        {0}".format(name))
    print("    " + detail.replace("\n", "\n    "))


# ---------------------------------------------------------------------------------------
# Fixtures. Everything that can come out of the database comes out of the database.
# ---------------------------------------------------------------------------------------

def loaded_families():
    families_path = os.path.join(settings['database.directory'], 'kinetics', 'families')
    database = KineticsDatabase()
    database.load_families(path=families_path, families=[FAMILY, IONISATION])
    return database


def generated(family):
    """
    A reaction the family generated, from molecules that have been through production's
    own `generate_resonance_structures()` -- which is what assigns atom ids.

    Called with **default** arguments, so the state under test is the state production
    hands its consumers. Round 110's own verification probe passed ``delete_labels=False``
    and reported green on a path production never takes.
    """
    reactants = []
    for smiles in ('[Li+]', '[CH3]'):
        species = Species(molecule=[Molecule(smiles=smiles)])
        species.generate_resonance_structures()
        reactants.append(species.molecule[0])
    return family.generate_reactions(reactants)[0]


def ionisation_view(family):
    """
    The campaign's own reaction, as the reactor sees it: ``Li + e- => Li+ + e- + e-``.

    `resolve_electron_placement` appends the *same* canonical electron object to both
    sides (``electron_placement.py:_place_declared_electrons``), which is the production
    source of a species on both sides. Nothing here is hand-assembled: the reaction is
    generated, its rate comes from the family, and the placement is the shipped resolver's.
    """
    reaction = family.generate_reactions([Molecule(smiles='[Li]')])[0]
    reaction.kinetics = family.get_kinetics(
        reaction, template_labels=reaction.template, degeneracy=reaction.degeneracy)[0][0]
    reaction.ensure_species()
    electron = Species(label='e', molecule=[Molecule(smiles='e')])
    view = resolve_electron_placement(
        reaction, [electron] + list(reaction.reactants) + list(reaction.products))
    view.generate_pairs()
    return view


def atoms_of(reaction):
    atoms = []
    for structure in list(reaction.reactants) + list(reaction.products):
        for molecule in (structure.molecule if isinstance(structure, Species)
                         else [structure]):
            atoms.extend(molecule.atoms)
    return atoms


def through_the_parallel_transport(value):
    """
    The value as `Pool.map` would hand it back: `multiprocessing`'s own pickler.

    ``ForkingPickler`` is what ``react.py:69`` serialises with, in both directions. Using
    it directly rather than starting a pool keeps the measurement about the transport and
    not about the pool; `the_parallel_run_loses_it` below does start one.
    """
    return ForkingPickler.loads(ForkingPickler.dumps(value))


def payload_of(atoms):
    return [(atom.id, dict(atom.props)) for atom in atoms]


# ---------------------------------------------------------------------------------------
# HIGH 1 -- the transport beside the one that was repaired
# ---------------------------------------------------------------------------------------

def the_parallel_transport_drops_the_payload(family):
    reaction = generated(family)
    before = payload_of(atoms_of(reaction))
    after = payload_of(atoms_of(through_the_parallel_transport(reaction)))
    lost = [(b, a) for b, a in zip(before, after) if b != a]
    record('the parallel transport drops the per-atom payload `copy()` now keeps',
           bool(lost),
           'of {0} atoms, {1} came back different. ids {2} -> {3}; props {4} -> {5}'.format(
               len(before), len(lost),
               [b[0] for b in before[:4]], [a[0] for a in after[:4]],
               [b[1] for b in before[:2]], [a[1] for a in after[:2]])
           if lost else
           'all {0} atoms came back with their ids and props intact'.format(len(before)))


def the_parallel_run_loses_it(database):
    """
    The same measurement through an actual `Pool`, so the finding is about production's
    path and not about a pickler chosen to represent it.
    """
    import rmgpy.data.rmg
    from rmgpy.rmg.react import react

    previous = getattr(rmgpy.data.rmg, 'database', None)

    class _Stub:
        pass

    stub = _Stub()
    stub.kinetics = database
    rmgpy.data.rmg.database = stub
    try:
        reactants = []
        for smiles in ('[Li+]', '[CH3]'):
            species = Species(molecule=[Molecule(smiles=smiles)])
            species.generate_resonance_structures()
            reactants.append(species)
        task = [((reactants[0], reactants[1]), [FAMILY])]
        serial = react(list(task), procnum=1)[0]
        parallel = react(list(task), procnum=2)[0]
    finally:
        rmgpy.data.rmg.database = previous

    if not serial or not parallel:
        record('a parallel run differs from the serial run it is supposed to reproduce',
               False, 'the generation produced nothing, so nothing was compared')
        return
    one = payload_of(atoms_of(serial[0]))
    two = payload_of(atoms_of(parallel[0]))
    record('a parallel run differs from the serial run it is supposed to reproduce',
           one != two,
           'serial ids {0} props {1}\nparallel ids {2} props {3}'.format(
               [x[0] for x in one[:4]], [x[1] for x in one[:2]],
               [x[0] for x in two[:4]], [x[1] for x in two[:2]]))


def the_parallel_transport_cannot_carry_a_fragment():
    fragment = Fragment().from_smiles_like_string('CCR')
    reaction = Reaction(reactants=[fragment],
                        products=[Fragment().from_smiles_like_string('[CH3]')])
    try:
        after = through_the_parallel_transport(reaction)
    except Exception as error:                                   # noqa: BLE001
        record('a fragment reaction cannot cross the parallel transport at all', True,
               'it raised {0}: {1}. A parallel run over fragment chemistry dies where '
               'the serial run succeeds.'.format(type(error).__name__, error))
        return
    kinds = [type(atom).__name__ for atom in after.reactants[0].atoms]
    wrong = (not isinstance(after.reactants[0], Fragment)
             or 'CuttingLabel' not in kinds)
    record('a fragment reaction cannot cross the parallel transport at all', wrong,
           'it came back a {0} holding {1}'.format(type(after.reactants[0]).__name__, kinds))


def the_parallel_transport_cannot_carry_a_surface_molecule():
    molecule = Molecule(smiles='C')
    molecule.metal = 'Pt'
    molecule.facet = '111'
    reaction = Reaction(reactants=[molecule], products=[Molecule(smiles='C')])
    try:
        after = through_the_parallel_transport(reaction)
    except Exception as error:                                   # noqa: BLE001
        record('a surface reaction cannot cross the parallel transport at all', True,
               'it raised {0}: {1}'.format(type(error).__name__, error))
        return
    back = after.reactants[0]
    record('a surface reaction cannot cross the parallel transport at all',
           (back.metal, back.facet) != ('Pt', '111'),
           'it came back with metal={0!r} facet={1!r}'.format(back.metal, back.facet))


# ---------------------------------------------------------------------------------------
# HIGH 2 -- one species on both sides
# ---------------------------------------------------------------------------------------

def the_copy_maps_a_reactant_onto_a_product(view):
    copied = view.copy()
    detached = []
    for index, pair in enumerate(copied.pairs):
        if not any(species is pair[0] for species in copied.reactants):
            detached.append((index, str(pair[0])))
    try:
        copied.reactants.index(copied.pairs[0][0])
        raised = None
    except ValueError as error:
        raised = error
    record('a species on both sides loses its reactant-side pair member',
           bool(detached),
           '{0!s}: pairs {1}; detached {2}; reactants.index(pair[0]) -> {3}'.format(
               view, [(str(a), str(b)) for a, b in copied.pairs], detached,
               'ValueError' if raised is not None else 'ok')
           if detached else
           'every pairs[0] member is one of the copy\'s own reactants')


def the_copy_aliases_the_collider(view):
    reaction = view.copy()
    reaction.specific_collider = Species(label='M', molecule=[Molecule(smiles='[He]')])
    copied = reaction.copy()
    record('the deep copy aliases specific_collider',
           copied.specific_collider is reaction.specific_collider,
           'copy().specific_collider is reaction.specific_collider -> {0}. The docstring '
           'says "deep copy"; an edit to the copy\'s collider reaches the '
           'original.'.format(copied.specific_collider is reaction.specific_collider))


# ---------------------------------------------------------------------------------------
# HIGH 3 -- the family is on disk and the enumeration cannot see it
# ---------------------------------------------------------------------------------------

MANIFEST = '''
name = "{label}/quarantine"
state = "QUARANTINED FOR TESTING"
appliesToKineticsClass = "Marcus"
reason = "a reason that must reach the warning"
shortDesc = "short"
longDesc = "long"
'''


def make_marcus():
    from rmgpy.kinetics import Marcus
    return Marcus(A=(1.73e06, 'm^3/(mol*s)'), n=2,
                  lmbd_i_coefs=[21824.5, -0.0341626, -0.0013254, 4.92966e-07],
                  beta=(1.2e10, '1/m'), wr=(0, 'kJ/mol'), wp=(0, 'kJ/mol'),
                  lmbd_o=(0, 'J/mol'), comment='no authorship recorded here')


class _Recorder(logging.Handler):
    def __init__(self):
        logging.Handler.__init__(self)
        self.messages = []

    def emit(self, record_):
        self.messages.append(record_.getMessage())


def _database_on_disk(root, label, manifest=True, unreadable=False):
    """A families directory carrying one family, with or without a readable manifest."""
    family_path = os.path.join(root, 'kinetics', 'families', label)
    os.makedirs(family_path)
    if manifest:
        target = os.path.join(root, 'elsewhere.py')
        with open(target, 'w') as handle:
            handle.write(MANIFEST.format(label=label))
        if unreadable:
            # A symlink is refused by the O_NOFOLLOW descent -- the manifest is there and
            # cannot be read, which is exactly the state that must not read as absence.
            os.symlink(target, os.path.join(family_path, QUARANTINE_FILENAME))
        else:
            shutil.copy(target, os.path.join(family_path, QUARANTINE_FILENAME))
    return family_path


def _with_database_directory(root, families):
    """Point `settings` at `root` and the loaded-family registry at `families`."""
    import rmgpy.data.rmg

    class _Kinetics:
        def __init__(self, f):
            self.families = f

    class _Database:
        def __init__(self, f):
            self.kinetics = _Kinetics(f)

    previous_directory = settings['database.directory']
    previous_database = getattr(rmgpy.data.rmg, 'database', None)
    settings['database.directory'] = root
    rmgpy.data.rmg.database = _Database(families)
    return previous_directory, previous_database


def _restore(previous):
    import rmgpy.data.rmg
    settings['database.directory'] = previous[0]
    rmgpy.data.rmg.database = previous[1]


def _clear_caches():
    from rmgpy.data.kinetics import quarantine as module
    module._DISK_QUARANTINE_CACHE.clear()
    module._DISK_ANY_QUARANTINE_CACHE.clear()
    module._UNATTRIBUTED_WARNED.clear()
    module._UNANSWERED_WARNED.clear()
    module._UNSAFE_MANIFESTS_WARNED.clear()


def an_unattributed_rate_passes_in_silence():
    root = tempfile.mkdtemp()
    _database_on_disk(root, 'A_Family_Quarantined_On_Disk')
    previous = _with_database_directory(root, {})       # nothing loaded, manifest on disk
    recorder = _Recorder()
    logging.getLogger().addHandler(recorder)
    logging.disable(logging.NOTSET)
    try:
        _clear_caches()
        reaction = LibraryReaction(
            reactants=[Species(label='Lip', molecule=[Molecule(smiles='[Li+]')])],
            products=[Species(label='CH3Li', molecule=[Molecule(smiles='C[Li]')])],
            library='a_seed_with_no_family_line',
            kinetics=make_marcus(),
            entry=Entry(index=1, label='no authorship', long_desc='nothing here'))
        check_quarantine(reaction, 'test')
    finally:
        logging.disable(logging.INFO)
        logging.getLogger().removeHandler(recorder)
        _restore(previous)
        shutil.rmtree(root, ignore_errors=True)
        _clear_caches()
    said = [m for m in recorder.messages if 'quarantin' in m.lower()]
    record('an unattributed rate passes in silence when the family is on disk, unloaded',
           not said,
           'check_quarantine() logged nothing about a quarantine. The manifest is on '
           'disk, its criterion is the rate\'s own kinetics class, and the entry records '
           'no author -- the exact state the gate exists to report.'
           if not said else
           'it warned: {0}'.format(said[0][:180]))


# ---------------------------------------------------------------------------------------
# MEDIUM 1 -- an unreadable manifest is an absent one
# ---------------------------------------------------------------------------------------

def an_unreadable_manifest_reads_as_an_absent_one():
    root = tempfile.mkdtemp()
    unreadable = _database_on_disk(root, 'Unreadable', manifest=True, unreadable=True)
    absent = _database_on_disk(root, 'Absent', manifest=False)
    previous = _with_database_directory(root, {})
    recorder = _Recorder()
    logging.getLogger().addHandler(recorder)
    logging.disable(logging.NOTSET)
    try:
        _clear_caches()
        refused = load_family_quarantine('Unreadable', unreadable)
        missing = load_family_quarantine('Absent', absent)

        class _Family:
            label = 'Unreadable'

            def __init__(self, q):
                self.quarantine = q

        reaction = LibraryReaction(
            reactants=[Species(label='Lip', molecule=[Molecule(smiles='[Li+]')])],
            products=[Species(label='CH3Li', molecule=[Molecule(smiles='C[Li]')])],
            library='some_library', kinetics=make_marcus())
        before = len(recorder.messages)
        check_quarantine(reaction, 'test', family=_Family(refused))
        answered_clean = not [m for m in recorder.messages[before:]
                              if 'Cannot tell' in m or 'cannot be consulted' in m]
    finally:
        logging.disable(logging.INFO)
        logging.getLogger().removeHandler(recorder)
        _restore(previous)
        shutil.rmtree(root, ignore_errors=True)
        _clear_caches()
    record('a manifest that could not be read is indistinguishable from one that is absent',
           refused is missing and answered_clean,
           'load_family_quarantine() returned {0!r} for the refused manifest and {1!r} '
           'for the family that has none, and check_quarantine(..., family=...) then '
           'answered without a word: an error became a clean bill of '
           'health.'.format(refused, missing))


# ---------------------------------------------------------------------------------------
# MEDIUM 2 -- the two subclasses the workaround never reached
# ---------------------------------------------------------------------------------------

def _depository_reaction():
    depository = KineticsDepository(label='Some_Family/training')
    entry = Entry(index=7, label='a training reaction')
    reaction = DepositoryReaction(
        reactants=[Species(label='Lip', molecule=[Molecule(smiles='[Li+]')])],
        products=[Species(label='CH3Li', molecule=[Molecule(smiles='C[Li]')])],
        depository=depository, family='Some_Family', entry=entry)
    reaction.allow_max_rate_violation = True
    reaction.rank = 3
    reaction.is_forward = True
    reaction.comment = 'a comment'
    reaction.label = 'a label'
    return reaction


def _pdep_reaction():
    reaction = PDepReaction(
        reactants=[Species(label='Lip', molecule=[Molecule(smiles='[Li+]')])],
        products=[Species(label='CH3Li', molecule=[Molecule(smiles='C[Li]')])],
        network='a network stand-in')
    reaction.allow_max_rate_violation = True
    reaction.rank = 3
    reaction.comment = 'a comment'
    return reaction


def _subclass_copy_loses_the_class(name, reaction, extra):
    copied = reaction.copy()
    wrong = type(copied) is not type(reaction)
    lost = [field for field in extra
            if getattr(copied, field, None) is not getattr(reaction, field, None)]
    record('{0}.copy() returns a base Reaction'.format(name), wrong or bool(lost),
           'copy() returned a {0}; {1} did not survive'.format(
               type(copied).__name__, lost or 'nothing else'))


def _subclass_reduce_drops_fields(name, reaction, fields):
    back = pickle.loads(pickle.dumps(reaction))
    lost = [field for field in fields
            if getattr(back, field, None) != getattr(reaction, field, None)]
    record('{0}.__reduce__ drops fields it was never told about'.format(name), bool(lost),
           '{0} came back changed: {1}'.format(
               lost, [(f, getattr(reaction, f, None), getattr(back, f, None))
                      for f in lost]))


# ---------------------------------------------------------------------------------------
# Controls -- these must hold at the base and at the tip
# ---------------------------------------------------------------------------------------

def control_the_serial_copy_still_carries_the_payload(family):
    reaction = generated(family)
    before = payload_of(atoms_of(reaction))
    after = payload_of(atoms_of(reaction.copy()))
    control("round 110's property: copy() keeps the payload",
            before == after,
            '{0} atoms, ids {1}, props {2}'.format(
                len(before), [b[0] for b in before[:4]], [b[1] for b in before[:2]]))


def control_distinct_species_were_never_the_problem():
    reactant = Species(label='A', molecule=[Molecule(smiles='[Li+]')])
    product = Species(label='B', molecule=[Molecule(smiles='C[Li]')])
    reaction = Reaction(reactants=[reactant], products=[product])
    reaction.generate_pairs()
    copied = reaction.copy()
    ok = all(any(species is pair[0] for species in copied.reactants)
             and any(species is pair[1] for species in copied.products)
             for pair in copied.pairs)
    control('a reaction whose species are all distinct copies correctly', ok,
            'every pair member is one of the copy\'s own species -- which is why the '
            'existing test could not see the defect above')


def control_the_fixture_really_shares_one_species(view):
    shared = [str(r) for r in view.reactants for p in view.products if r is p]
    control('the fixture really does hold one Species object on both sides', bool(shared),
            'shared: {0}; this is the resolver\'s own view object, not an '
            'assembled one'.format(shared))


def control_an_absent_manifest_is_still_an_answer():
    root = tempfile.mkdtemp()
    path = _database_on_disk(root, 'Ordinary_Family', manifest=False)
    previous = _with_database_directory(root, {})
    try:
        _clear_caches()
        from rmgpy.data.kinetics.quarantine import resolve_quarantine
        answer = resolve_quarantine('Ordinary_Family')
        loaded = load_family_quarantine('Ordinary_Family', path)
    finally:
        _restore(previous)
        shutil.rmtree(root, ignore_errors=True)
        _clear_caches()
    control('a family with no manifest is still answered, and answered clean',
            answer == (None, True) and loaded is None,
            'resolve_quarantine -> {0!r}, load_family_quarantine -> {1!r}'.format(
                answer, loaded))


def control_a_loaded_quarantine_still_refuses():
    from rmgpy.exceptions import QuarantinedKineticsError
    root = tempfile.mkdtemp()
    path = _database_on_disk(root, 'Loaded_Family')
    quarantine = load_family_quarantine('Loaded_Family', path)

    class _Family:
        label = 'Loaded_Family'

        def __init__(self, q):
            self.quarantine = q

    previous = _with_database_directory(root, {'Loaded_Family': _Family(quarantine)})
    try:
        _clear_caches()
        reaction = LibraryReaction(
            reactants=[Species(label='Lip', molecule=[Molecule(smiles='[Li+]')])],
            products=[Species(label='CH3Li', molecule=[Molecule(smiles='C[Li]')])],
            library='some_library', kinetics=make_marcus())
        try:
            check_quarantine(reaction, 'test', family=_Family(quarantine))
            refused = None
        except QuarantinedKineticsError as error:
            refused = error
    finally:
        _restore(previous)
        shutil.rmtree(root, ignore_errors=True)
        _clear_caches()
    control('the gate still refuses a rate from a loaded quarantined family',
            refused is not None,
            'check_quarantine raised {0}'.format(type(refused).__name__ if refused
                                                 else 'nothing'))


def control_the_parallel_transport_still_carries_an_ordinary_object():
    reaction = Reaction(reactants=[Species(label='A', molecule=[Molecule(smiles='C')])],
                        products=[Species(label='B', molecule=[Molecule(smiles='C')])],
                        degeneracy=3)
    back = through_the_parallel_transport(reaction)
    control('the parallel transport still round-trips an ordinary reaction',
            len(back.reactants) == 1 and back.degeneracy == 3,
            'reactants {0}, degeneracy {1}'.format(len(back.reactants), back.degeneracy))


# ---------------------------------------------------------------------------------------
# The census
# ---------------------------------------------------------------------------------------

def the_census():
    from rmgpy.data.kinetics.family import _LOSSY_REDUCERS
    note('the hand-enumerated half, and how an unlisted lossy class fails',
         'the table holds {0} classes: {1}. Membership is a measurement for `Atom` alone '
         '-- `fields_the_reducer_drops` probes it -- and a judgement for the other four. '
         'A class that starts losing a field tomorrow and is not in the table loses it '
         'SILENTLY: nothing raises, the field simply arrives as its default. The census '
         'test in quarantineTest.py is what turns that into a red test rather than a '
         'wrong number.'.format(
             len(_LOSSY_REDUCERS),
             ', '.join(sorted(cls.__name__ for cls in _LOSSY_REDUCERS))))


def main():
    database = loaded_families()
    family = database.families[FAMILY]
    view = ionisation_view(database.families[IONISATION])

    the_parallel_transport_drops_the_payload(family)
    the_parallel_run_loses_it(database)
    the_parallel_transport_cannot_carry_a_fragment()
    the_parallel_transport_cannot_carry_a_surface_molecule()
    the_copy_maps_a_reactant_onto_a_product(view)
    the_copy_aliases_the_collider(view)
    an_unattributed_rate_passes_in_silence()
    an_unreadable_manifest_reads_as_an_absent_one()
    _subclass_copy_loses_the_class('DepositoryReaction', _depository_reaction(),
                                   ('depository', 'entry', 'family'))
    _subclass_copy_loses_the_class('PDepReaction', _pdep_reaction(), ('network',))
    _subclass_reduce_drops_fields('DepositoryReaction', _depository_reaction(),
                                  ('allow_max_rate_violation', 'rank', 'comment', 'label',
                                   'elementary_high_p', 'allow_pdep_route',
                                   'network_kinetics'))
    _subclass_reduce_drops_fields('PDepReaction', _pdep_reaction(),
                                  ('allow_max_rate_violation', 'rank', 'comment',
                                   'elementary_high_p', 'allow_pdep_route'))

    control_the_serial_copy_still_carries_the_payload(family)
    control_distinct_species_were_never_the_problem()
    control_the_fixture_really_shares_one_species(view)
    control_an_absent_manifest_is_still_an_answer()
    control_a_loaded_quarantine_still_refuses()
    control_the_parallel_transport_still_carries_an_ordinary_object()

    the_census()

    reproduced = [name for name, flag, _ in results if flag]
    broken = [name for name, ok, _ in controls if not ok]
    print("\n" + "=" * 78)
    print("findings reproduced : {0} of {1}".format(len(reproduced), len(results)))
    print("controls holding    : {0} of {1}".format(
        len(controls) - len(broken), len(controls)))
    print("noted, not findings : {0}".format(len(notes)))
    for name in reproduced:
        print("  REPRODUCED  {0}".format(name))
    for name in broken:
        print("  CONTROL BROKEN  {0}".format(name))
    if broken:
        return 2
    return 1 if reproduced else 0


if __name__ == '__main__':
    sys.exit(main())
