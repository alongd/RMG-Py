#!/usr/bin/env python3
"""
Round 110 probe -- the round trip carries the graph and drops the payload.

    MPLCONFIGDIR=/tmp/claude-1000/mpl PYTHONPATH=$PWD \\
        ~/anaconda3/envs/rmg_env/bin/python \\
        docs/i221-saturation/probes/round110_transport_payload_probe.py

Round 108 made `copy()` reproduce the reaction through one `pickle` round trip so that two
fields referring to one object still refer to one object on the other side. That property
holds. The transport it chose, though, delegates to each class's own ``__reduce__``, and
those reducers are incomplete for one class, *wrong* for two, and misrouted for one field
pair:

* ``Atom.__reduce__`` never mentions ``id``, ``coords`` or ``props``. Atom ids drive
  resonance-structure correspondence and ``props`` carries ``'inRing'``, which feeds group
  matching, so this is a silent wrong answer rather than a crash.
* ``Fragment`` and ``CuttingLabel`` inherit reducers naming ``Molecule`` and ``Atom`` as
  the class to rebuild; a cutting label's symbol is then handed to ``get_element`` and
  raises ``KeyError: 'R'``.
* ``Molecule.__reduce__`` passes ``metal`` and ``facet`` into ``__init__``'s ``inchi`` and
  ``smiles`` parameters -- the fifth and sixth positional, not the seventh and eighth --
  so a surface molecule is rebuilt from ``metal='Pt'`` read as an InChI.
* ``Species.__reduce__`` omits ``aug_inchi``, ``creation_iteration``,
  ``explicitly_allowed`` and ``symmetry_number``.

Two further findings are this branch's own regressions rather than upstream's, and one
belongs to the partition rather than the transport. All of them are measured here against
the shipped database, on reactions the family itself generated.

``REPRODUCED`` means the defect is present. Controls must hold at both tips. Exit 1 while
any finding is present, 2 if a control breaks.
"""

import copy as copy_module
import logging
import os
import pickle
import sys
import tempfile

import numpy as np

logging.disable(logging.INFO)

from rmgpy import settings
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import KineticsLibrary, LibraryReaction
from rmgpy.data.kinetics.quarantine import QUARANTINE_FILENAME, _read_manifest
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.molecule.atomtype import ATOMTYPES
from rmgpy.molecule.fragment import CuttingLabel, Fragment
from rmgpy.molecule.molecule import Atom, Bond
from rmgpy.reaction import Reaction
from rmgpy.species import Species

FAMILY = 'Cation_R_Recombination'

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
# Fixtures. Everything below that can come out of the family comes out of the family.
# ---------------------------------------------------------------------------------------

def loaded_family():
    families_path = os.path.join(settings['database.directory'], 'kinetics', 'families')
    database = KineticsDatabase()
    database.load_families(path=families_path, families=[FAMILY])
    return database.families[FAMILY]


def generated(family, delete_labels=True):
    """
    A reaction the family generated, from molecules that have been through production's
    own `generate_resonance_structures()` -- which is what assigns atom ids. Nothing here
    is assigned by the fixture.
    """
    reactants = []
    for smiles in ('[Li+]', '[CH3]'):
        species = Species(molecule=[Molecule(smiles=smiles)])
        species.generate_resonance_structures()
        reactants.append(species.molecule[0])
    return family.generate_reactions(reactants, delete_labels=delete_labels)[0]


def atoms_of(reaction):
    atoms = []
    for structure in list(reaction.reactants) + list(reaction.products):
        for molecule in (structure.molecule if isinstance(structure, Species)
                         else [structure]):
            atoms.extend(molecule.atoms)
    return atoms


def copied(reaction):
    """`reaction.copy()`, or the exception it raised."""
    try:
        return reaction.copy(), None
    except Exception as error:                                   # noqa: BLE001
        return None, error


# ---------------------------------------------------------------------------------------
# HIGH 1 -- the per-atom state the old path preserved
# ---------------------------------------------------------------------------------------

def the_copy_erases_atom_state(family):
    reaction = generated(family, delete_labels=False)
    before = [(atom.id, dict(atom.props)) for atom in atoms_of(reaction)]
    after_reaction, error = copied(reaction)
    if error is not None:
        record('copy() erases per-atom state that the old path preserved', True,
               'copy() raised {0}: {1} -- the state cannot even be compared'.format(
                   type(error).__name__, error))
        return
    after = [(atom.id, dict(atom.props)) for atom in atoms_of(after_reaction)]
    lost_ids = [(b[0], a[0]) for b, a in zip(before, after) if b[0] != a[0]]
    lost_props = [(b[1], a[1]) for b, a in zip(before, after) if b[1] != a[1]]
    record('copy() erases per-atom state that the old path preserved',
           bool(lost_ids or lost_props),
           'atom ids changed on {0} of {1} atoms {2}; props changed on {3} {4}. Ids drive '
           'resonance-structure correspondence and props carries \'inRing\', which feeds '
           'group matching.'.format(
               len(lost_ids), len(before), lost_ids[:3], len(lost_props), lost_props[:2])
           if (lost_ids or lost_props) else
           'every atom of the copy carries the id and the props the original\'s did')


def control_the_state_is_productions_and_the_old_path_kept_it(family):
    reaction = generated(family, delete_labels=False)
    atoms = atoms_of(reaction)
    present = [atom for atom in atoms if atom.id != -1 or atom.props]
    old = []
    for structure in reaction.reactants:
        old.extend(structure.copy(deep=True).atoms)
    kept = all(a.id == b.id and a.props == b.props
               for a, b in zip(atoms_of(reaction)[:len(old)], old))
    control('the atom state is production\'s, and Molecule.copy(deep=True) keeps it',
            bool(present) and kept,
            '{0} of {1} atoms carry an assigned id or non-empty props before any copy, '
            'and the deep copy `copy()` used before round 108 reproduces both. So the '
            'loss above is the new transport\'s, not an absence in the fixture.'.format(
                len(present), len(atoms))
            if present and kept else
            'the fixture does not carry the state under test ({0} atoms with state, old '
            'path kept it: {1})'.format(len(present), kept))


# ---------------------------------------------------------------------------------------
# HIGH 2 -- a fragment reaction cannot be copied at all
# ---------------------------------------------------------------------------------------

def a_fragment_reaction_cannot_be_copied():
    fragment = Fragment().from_smiles_like_string('CCR')
    reaction = TemplateReaction(reactants=[fragment],
                                products=[Fragment().from_smiles_like_string('[CH3]')],
                                family='A_Family')
    reaction.labeled_atoms = {'reactants': {}, 'products': {}}
    after, error = copied(reaction)
    if error is not None:
        record('a fragment reaction cannot be copied at all', True,
               'copy() raised {0}: {1}. `CuttingLabel` inherits `Atom.__reduce__`, which '
               'names `Atom` as the class to rebuild and hands it the cutting label\'s '
               'symbol as an element.'.format(type(error).__name__, error))
        return
    kinds = [type(atom).__name__ for atom in after.reactants[0].atoms]
    wrong = (not isinstance(after.reactants[0], Fragment)
             or 'CuttingLabel' not in kinds)
    record('a fragment reaction cannot be copied at all', wrong,
           'the copy came back as {0} with atoms {1}'.format(
               type(after.reactants[0]).__name__, sorted(set(kinds)))
           if wrong else
           'the copy is a Fragment and keeps its CuttingLabel, and the cutting label '
           'keeps its name ({0})'.format(
               [a.name for a in after.reactants[0].atoms
                if isinstance(a, CuttingLabel)]))


def control_the_previous_path_could_copy_a_fragment():
    fragment = Fragment().from_smiles_like_string('CCR')
    deep = copy_module.deepcopy(fragment)
    ok = (isinstance(deep, Fragment)
          and any(isinstance(atom, CuttingLabel) for atom in deep.atoms))
    control('the deep copy this replaced could copy a fragment', ok,
            'deepcopy(fragment) is a Fragment whose cutting label survives, so fragment '
            'copying is a working feature and not something that never worked'
            if ok else 'deepcopy of a fragment does not preserve it either')


# ---------------------------------------------------------------------------------------
# The regressions this branch introduced, found while measuring the two above
# ---------------------------------------------------------------------------------------

def the_familys_own_output_cannot_be_copied(family):
    reaction = generated(family)                     # production defaults
    after, error = copied(reaction)
    record('a reaction the family generated cannot be copied', error is not None,
           'copy() raised {0}: {1}. `family.py:2651` deletes `labeled_atoms` once the '
           'labels have been read back, so every reaction `generate_reactions()` returns '
           'is missing it -- and the deepened set is applied from the table rather than '
           'from what the object holds.'.format(type(error).__name__, error)
           if error is not None else
           'copy() reproduces the family\'s own output, `labeled_atoms` absent and all')


def a_surface_molecule_cannot_be_copied():
    molecule = Molecule(smiles='CC')
    molecule.metal, molecule.facet = 'Pt', '111'
    reaction = TemplateReaction(reactants=[molecule],
                                products=[Molecule(smiles='[CH3]')], family='A_Family')
    reaction.labeled_atoms = {'reactants': {}, 'products': {}}
    after, error = copied(reaction)
    if error is not None:
        record('a reaction over a surface molecule cannot be copied', True,
               'copy() raised {0}: {1}. `Molecule.__reduce__` passes metal and facet into '
               '__init__(..., inchi, smiles), which are the fifth and sixth positional '
               'parameters rather than the seventh and eighth.'.format(
                   type(error).__name__, error))
        return
    lost = (after.reactants[0].metal, after.reactants[0].facet) != ('Pt', '111')
    record('a reaction over a surface molecule cannot be copied', lost,
           'the copy came back with metal={0!r} facet={1!r}'.format(
               after.reactants[0].metal, after.reactants[0].facet) if lost else
           'the copy keeps metal=\'Pt\' and facet=\'111\'')


def a_species_loses_its_own_state():
    reactant = Species(label='ethane', molecule=[Molecule(smiles='CC')])
    reactant.symmetry_number = 2
    reactant.aug_inchi = 'an augmented inchi'
    reactant.creation_iteration = 3
    reactant.explicitly_allowed = True
    reaction = TemplateReaction(
        reactants=[reactant],
        products=[Species(label='ethyl', molecule=[Molecule(smiles='C[CH2]')])],
        family='A_Family')
    reaction.labeled_atoms = {'reactants': {}, 'products': {}}
    after, error = copied(reaction)
    if error is not None:
        record('a copied Species loses state its own reducer omits', True,
               'copy() raised {0}: {1}'.format(type(error).__name__, error))
        return
    copy = after.reactants[0]
    lost = [name for name in ('symmetry_number', 'aug_inchi', 'creation_iteration',
                              'explicitly_allowed')
            if getattr(copy, name) != getattr(reactant, name)]
    record('a copied Species loses state its own reducer omits', bool(lost),
           'lost: {0}'.format(lost) if lost else
           'symmetry_number, aug_inchi, creation_iteration and explicitly_allowed all '
           'survive the copy')


# ---------------------------------------------------------------------------------------
# MEDIUM 1 -- the partition's default
# ---------------------------------------------------------------------------------------

def an_unclassified_field_is_aliased(family):
    reaction = generated(family, delete_labels=False)
    after, error = copied(reaction)
    if error is not None:
        record('a field in neither table is silently aliased', True,
               'copy() raised {0}: {1}'.format(type(error).__name__, error))
        return
    aliased = [name for name in ('template',)
               if getattr(after, name, None) is getattr(reaction, name, None)
               and getattr(reaction, name, None) is not None]
    record('a field in neither table is silently aliased', bool(aliased),
           '{0} is the *same list object* on the copy and the original ({1!r}), so an '
           'edit to one is an edit to the other'.format(aliased, reaction.template)
           if aliased else
           'template is a list of its own on the copy; the partition has no silent '
           'default left')


def an_unknown_field_is_accepted_without_a_word():
    reaction = TemplateReaction(
        reactants=[Species(label='a', molecule=[Molecule(smiles='CC')])],
        products=[Species(label='b', molecule=[Molecule(smiles='C[CH2]')])],
        family='A_Family')
    reaction.labeled_atoms = {'reactants': {}, 'products': {}}
    reaction.a_field_nobody_classified = ['a mutable one']
    after, error = copied(reaction)
    if error is None:
        shared = getattr(after, 'a_field_nobody_classified', None) is \
            reaction.a_field_nobody_classified
        record('a field nobody classified is carried by reference rather than refused',
               shared,
               'the copy shares the original\'s list; the default for an unclassified '
               'field is alias, and it is silent' if shared else
               'the copy has a list of its own, but nothing announced the classification')
    else:
        record('a field nobody classified is carried by reference rather than refused',
               False,
               'copy() refused it loudly: {0}: {1}'.format(
                   type(error).__name__, str(error).split('.')[0]))


# ---------------------------------------------------------------------------------------
# MEDIUM 2 -- the manifest's identity is sampled before its bytes are read
# ---------------------------------------------------------------------------------------

MANIFEST = ("state = 'QUARANTINED FOR QUANTITATIVE PLASMA USE'\n"
            "reason = '{0}'\n"
            "kinetics_class = 'Marcus'\n")


def the_identity_is_sampled_before_the_content():
    with tempfile.TemporaryDirectory() as root:
        family_path = os.path.join(root, 'kinetics', 'families', 'A_Family')
        os.makedirs(family_path)
        path = os.path.join(family_path, QUARANTINE_FILENAME)
        with open(path, 'w') as handle:
            handle.write(MANIFEST.format('the reason the first version gave'))

        rewritten = MANIFEST.format('a reason nobody has ever approved, written in place')
        real_fstat = os.fstat
        state = {'rewrites': 0}

        def rewriting_fstat(fd):
            info = real_fstat(fd)
            if state['rewrites'] == 0 and os.path.samestat(info, os.stat(path)):
                state['rewrites'] += 1
                # In place, between the identity call and the read: same inode, new bytes.
                with open(path, 'w') as handle:
                    handle.write(rewritten)
            return info

        os.fstat = rewriting_fstat
        try:
            content, identity = _read_manifest(family_path)
        finally:
            os.fstat = real_fstat

        if content is None:
            record('the manifest\'s identity is sampled before its bytes are read', False,
                   'the read was refused once the file changed under it, so nothing was '
                   'cached and the caller gets an unanswered question rather than an '
                   'answer keyed to the wrong version')
            return
        current = os.stat(path)
        stale = identity != (current.st_mtime_ns, current.st_ctime_ns,
                             current.st_size, current.st_ino)
        record('the manifest\'s identity is sampled before its bytes are read',
               stale and rewritten in content,
               'the returned content is the rewritten version and the returned identity '
               'is the version before it. `resolve_quarantine` caches the first under the '
               'second, so a later lookup is a HIT that returns a quarantine the file on '
               'disk does not state -- round 99\'s defect one layer down.'
               if stale else
               'content and identity describe the same version of the file')


# ---------------------------------------------------------------------------------------
# The addendum -- allow_max_rate_violation
# ---------------------------------------------------------------------------------------

def base_reaction_copy_drops_the_flag():
    reaction = Reaction(
        reactants=[Species(label='a', molecule=[Molecule(smiles='CC')])],
        products=[Species(label='b', molecule=[Molecule(smiles='C[CH2]')])])
    reaction.allow_max_rate_violation = True
    reaction.rank = 5
    reaction.is_forward = True
    reaction.pairs = [(reaction.reactants[0], reaction.products[0])]
    copy = reaction.copy()
    lost = [name for name in ('allow_max_rate_violation', 'rank', 'is_forward')
            if getattr(copy, name) != getattr(reaction, name)]
    record('Reaction.copy drops allow_max_rate_violation', bool(lost),
           'lost: {0}. A rate deliberately allowed past the collision limit comes back '
           'forbidden to exceed it, and the default is the value that hides the '
           'loss.'.format(lost) if lost else
           'allow_max_rate_violation, rank and is_forward all survive Reaction.copy')


def base_reaction_copy_severs_pairs():
    reaction = Reaction(
        reactants=[Species(label='a', molecule=[Molecule(smiles='CC')])],
        products=[Species(label='b', molecule=[Molecule(smiles='C[CH2]')])])
    reaction.pairs = [(reaction.reactants[0], reaction.products[0])]
    copy = reaction.copy()
    owned = list(copy.reactants) + list(copy.products)
    stray = [pair for pair in copy.pairs
             for member in pair if not any(member is s for s in owned)]
    record('Reaction.copy severs pairs from the species it copied', bool(stray),
           'a pairs member is a species the copy does not own; Species.__eq__ is '
           'identity, so reactants.index(pair[0]) raises' if stray else
           'every pairs member is one of the copy\'s own species')


def control_the_loader_already_carries_the_flag():
    """The addendum's first site, measured on *this* branch rather than assumed."""
    item = Reaction(reactants=[Species(label='a', molecule=[Molecule(smiles='CC')])],
                    products=[Species(label='b', molecule=[Molecule(smiles='C[CH2]')])])
    item.allow_max_rate_violation = True
    entry = Entry(index=1, label='a <=> b', item=item,
                  data=Arrhenius(A=(1.0, 's^-1'), n=0, Ea=(0, 'kJ/mol')))
    library = KineticsLibrary(label='L')
    library.entries = {'a <=> b': entry}
    library.auto_generated = False
    reaction = library.get_library_reactions()[0]
    control('get_library_reactions already carries allow_max_rate_violation here',
            reaction.allow_max_rate_violation is True,
            'the constructor call does omit the flag, and `_carry_entry_fields` -- round '
            '102\'s derived carry -- supplies it immediately afterwards. The review is '
            'right about the constructor and right for any branch without that carry; on '
            'this one the site is already closed, at both tips.'
            if reaction.allow_max_rate_violation is True else
            'the flag does NOT reach the reaction the loader builds')


# ---------------------------------------------------------------------------------------
# deepcopy, the other entry point
# ---------------------------------------------------------------------------------------

def deepcopy_severs_the_labelled_atoms(family):
    reaction = generated(family, delete_labels=False)
    deep = copy_module.deepcopy(reaction)
    owned = atoms_of(deep)
    stray = [label for group in deep.labeled_atoms.values()
             for label, atom in group.items()
             if not any(atom is a for a in owned)]
    record('deepcopy(reaction) severs the labelled atoms', bool(stray),
           'labels whose atom is inside no molecule the copy owns: {0}. `deepcopy` '
           'recurses into `Molecule.__deepcopy__`, which discards the memo. '
           '`family.py:3854` and `database.py:755` both take this path.'.format(stray)
           if stray else
           'every labelled atom is inside one of the copy\'s own molecules, and the copy '
           'is still a {0}'.format(type(deep).__name__))


# ---------------------------------------------------------------------------------------
# Controls on what must not have changed
# ---------------------------------------------------------------------------------------

def control_round_108s_properties_still_hold(family):
    reaction = generated(family, delete_labels=False)
    after, error = copied(reaction)
    if error is not None:
        control('round 108\'s identity properties still hold', False,
                'copy() raised {0}'.format(type(error).__name__))
        return
    owned = list(after.reactants) + list(after.products)
    atoms = atoms_of(after)
    pairs_ok = all(any(member is s for s in owned)
                   for pair in after.pairs for member in pair)
    labels_ok = all(any(atom is a for a in atoms)
                    for group in after.labeled_atoms.values()
                    for atom in group.values())
    fresh = not any(a is b for a in owned
                    for b in list(reaction.reactants) + list(reaction.products))
    control('round 108\'s identity properties still hold',
            pairs_ok and labels_ok and fresh,
            'every pairs member is one of the copy\'s own species, every labelled atom is '
            'inside one of its molecules, and nothing is the original\'s object'
            if pairs_ok and labels_ok and fresh else
            'pairs {0}, labelled atoms {1}, distinct from the original {2}'.format(
                pairs_ok, labels_ok, fresh))


def control_the_atom_types_stay_interned(family):
    reaction = generated(family, delete_labels=False)
    after, error = copied(reaction)
    if error is not None:
        control('the copied atoms keep the interned atom types', False,
                'copy() raised {0}'.format(type(error).__name__))
        return
    interned = all(a.atomtype is b.atomtype
                   for a, b in zip(atoms_of(reaction), atoms_of(after)))
    control('the copied atoms keep the interned atom types', interned,
            'every copied atom\'s atomtype is the same object as the original\'s, which '
            'is what `Atom.__reduce__` buys by storing the label and looking it up again'
            if interned else
            'an atom type came back as a different object')


def control_a_whole_object_reducer_would_de_intern_them():
    """
    The negative control on the choice above: why `Atom` keeps its own reducer and is
    completed rather than replaced.
    """
    interned = ATOMTYPES['Cs']
    by_value = pickle.loads(pickle.dumps(interned))
    broken = (by_value is not interned
              and not interned.is_specific_case_of(by_value))
    control('carrying an atom type by value would break group matching', broken,
            'a pickled `AtomType` is a different object, `AtomType` inherits identity '
            'equality, and `is_specific_case_of` over the de-interned copy is False. A '
            'reducer that carried `Atom` state wholesale would buy the three lost fields '
            'at the price of silent group mismatches, so `Atom` keeps its own reducer and '
            'has what it drops added back.'
            if broken else
            'a pickled AtomType still compares as the interned one, so the interning '
            'this design protects has stopped mattering')


def control_value_comparison_is_blind_to_all_of_this(family):
    reaction = generated(family, delete_labels=False)
    after, error = copied(reaction)
    if error is not None:
        control('value comparison sees none of this', True,
                'copy() raised, so there is nothing to compare -- the point stands '
                'vacuously here and is made on the round 108 test set in the findings')
        return
    by_value = ([str(m) for m in after.reactants] == [str(m) for m in reaction.reactants]
                and sorted(after.labeled_atoms) == sorted(reaction.labeled_atoms))
    control('value comparison sees none of this', by_value,
            'the reactants compare equal as strings and labeled_atoms by key, at this tip '
            'and at the base alike. Every round 108 assertion of that shape is listed in '
            'the findings with what it cannot catch.')


# ---------------------------------------------------------------------------------------
# The census
# ---------------------------------------------------------------------------------------

PLAIN = (type(None), bool, int, float, complex, str, bytes, np.ndarray)


def is_plain(value):
    if isinstance(value, PLAIN):
        return True
    if isinstance(value, (list, tuple, set, frozenset)):
        return all(is_plain(item) for item in value)
    if isinstance(value, dict):
        return all(is_plain(key) and is_plain(item) for key, item in value.items())
    return False


def writable(cls, _cache={}):
    if cls not in _cache:
        probe = cls.__new__(cls)
        names = set()
        for name in dir(cls):
            if name.startswith('_'):
                continue
            if type(getattr(cls, name, None)).__name__ != 'getset_descriptor':
                continue
            try:
                setattr(probe, name, None)
            except AttributeError:
                continue
            except Exception:                                    # noqa: BLE001
                pass
            names.add(name)
        _cache[cls] = frozenset(names)
    return _cache[cls]


def plain_state(obj):
    state = {}
    for name in sorted(writable(type(obj)) | set(getattr(obj, '__dict__', {}))):
        try:
            value = getattr(obj, name)
        except Exception:                                        # noqa: BLE001
            continue
        if is_plain(value):
            state[name] = value
    return state


def same(one, other):
    if isinstance(one, np.ndarray) or isinstance(other, np.ndarray):
        return np.array_equal(one, other)
    try:
        return one is other or bool(one == other)
    except Exception:                                            # noqa: BLE001
        return False


def census_of_one(obj):
    try:
        back = pickle.loads(pickle.dumps(obj, pickle.HIGHEST_PROTOCOL))
    except Exception as error:                                   # noqa: BLE001
        return 'RAISES {0}: {1}'.format(type(error).__name__, error)
    if type(back) is not type(obj):
        return 'rebuilt as {0}'.format(type(back).__name__)
    before = plain_state(obj)
    lost = [name for name, value in before.items()
            if not same(value, getattr(back, name, None))]
    return ', '.join(sorted(lost)) if lost else 'nothing'


def the_census(family):
    reaction = generated(family, delete_labels=False)
    atom = atoms_of(reaction)[1]
    bond = list(atom.edges.values())[0] if atom.edges else None
    species = Species(label='x', molecule=[Molecule(smiles='CC')])
    species.symmetry_number, species.aug_inchi = 2, 'AI'
    species.creation_iteration, species.explicitly_allowed = 3, True
    surface = Molecule(smiles='CC')
    surface.metal, surface.facet = 'Pt', '111'
    probe_atom = Atom(element='C')
    probe_atom.id, probe_atom.props = -424242, {'inRing': False}
    probe_atom.coords = np.array([1.5, 2.5, 3.5])
    probe_atom.terminal = probe_atom.ignore = True

    rows = [
        ('Atom (production)', census_of_one(atom)),
        ('Atom (every plain field set)', census_of_one(probe_atom)),
        ('Bond', census_of_one(bond) if bond is not None else 'no bond in the fixture'),
        ('Molecule (ordinary)', census_of_one(reaction.reactants[1])),
        ("Molecule (metal='Pt')", census_of_one(surface)),
        ('Species', census_of_one(species)),
        ('Fragment', census_of_one(Fragment().from_smiles_like_string('CCR'))),
        ('CuttingLabel', census_of_one(
            [a for a in Fragment().from_smiles_like_string('CCR').atoms
             if isinstance(a, CuttingLabel)][0])),
        ('TransitionState', census_of_one(reaction.transition_state)
         if reaction.transition_state is not None else 'None on this reaction'),
        ('Arrhenius', census_of_one(Arrhenius(A=(1.0, 's^-1'), n=0, Ea=(0, 'kJ/mol'),
                                              comment='a comment'))),
        ('TemplateReaction', census_of_one(reaction)),
        ('LibraryReaction', census_of_one(LibraryReaction(
            reactants=list(reaction.reactants), products=list(reaction.products),
            library='L'))),
    ]
    note('every class that crosses the transport, and what its own reducer drops',
         '\n'.join('  {0:<30} {1}'.format(name, lost) for name, lost in rows)
         + '\n\n'
         'Plain fields only: an object whose __eq__ is identity compares unequal to its '
         'own faithful copy, so reference fields are audited by round 108\'s `is` '
         'assertions instead of here.')


def note_the_residue():
    note('what this repair does NOT reach',
         'A plain `pickle.dumps(reaction)` still reaches the lossy reducers. `__reduce__` '
         'cannot dictate how the objects nested inside its state are pickled, so only a '
         'pickler that carries the dispatch table -- which is what `copy()` uses -- gets '
         'the complete ones. Closing that would take either a change to '
         '`rmgpy/molecule/`, which this round is forbidden, or a process-wide '
         '`copyreg.pickle()` registration, which would change the behaviour of code that '
         'never imported this module.')


def main():
    family = loaded_family()

    the_copy_erases_atom_state(family)
    a_fragment_reaction_cannot_be_copied()
    the_familys_own_output_cannot_be_copied(family)
    a_surface_molecule_cannot_be_copied()
    a_species_loses_its_own_state()
    an_unclassified_field_is_aliased(family)
    an_unknown_field_is_accepted_without_a_word()
    the_identity_is_sampled_before_the_content()
    base_reaction_copy_drops_the_flag()
    base_reaction_copy_severs_pairs()
    deepcopy_severs_the_labelled_atoms(family)

    control_the_state_is_productions_and_the_old_path_kept_it(family)
    control_the_previous_path_could_copy_a_fragment()
    control_the_loader_already_carries_the_flag()
    control_round_108s_properties_still_hold(family)
    control_the_atom_types_stay_interned(family)
    control_a_whole_object_reducer_would_de_intern_them()
    control_value_comparison_is_blind_to_all_of_this(family)

    the_census(family)
    note_the_residue()

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
