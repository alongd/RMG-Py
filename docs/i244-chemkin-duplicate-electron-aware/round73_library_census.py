"""
Census: which installed kinetics-library entries change duplicate group when the
Chemkin group key keys each side as a MULTISET of participants rather than as an
ordered list?

Unit of the census is one library, because that is the unit in which a permuted
pair can be *authored*: two entries of the same library, same species on each
side, written in a different order. Entries of different libraries are not pooled,
because no deck is assembled that way.

Surrogate for identity. ``chemkin_duplicate_group_key`` keys species by ``id()``,
which is right for a deck -- where one Species object stands for one species --
and useless across freshly loaded library entries, which each build their own
objects. Inside one library the species *label* is the identity namespace: it is
what the library's own species dictionary is keyed by, and it is what the Chemkin
identifier is built from. So this census substitutes the label for ``id()`` and
leaves every other term of the key untouched. Labels that resolve to more than one
molecular formula inside a library would break that substitution; the run reports
how many it found.

Output: for each library, the number of entries whose duplicate ANSWER differs
between the two key forms, and the number of groups that merge. The denominator is
every reaction loaded.
"""
import logging
import os
import sys
import traceback
from collections import defaultdict

logging.getLogger().setLevel(logging.ERROR)

from rmgpy import settings
from rmgpy.chemkin import chemkin_duplicate_group_key
from rmgpy.data.kinetics import KineticsDatabase
from rmgpy.electron_balance import get_electron_placement_counts

LIBRARY_ROOT = os.path.join(settings['database.directory'], 'kinetics', 'libraries')
print('database libraries root :', LIBRARY_ROOT)
print('key in use is multiset  :', 'multiset' in (chemkin_duplicate_group_key.__doc__ or ''))


def side_token(species_list):
    return tuple(spc.label for spc in species_list)


def key(reaction, multiset):
    """The shipped key with ``id(spc)`` replaced by the species label."""
    kinetics = reaction.kinetics
    pressure_dependent = bool(kinetics is not None and kinetics.is_pressure_dependent())
    reactant_electrons, product_electrons = get_electron_placement_counts(reaction)
    reactants = side_token(reaction.reactants)
    products = side_token(reaction.products)
    if multiset:
        reactants, products = tuple(sorted(reactants)), tuple(sorted(products))
    reactant_side = (reactants, reactant_electrons)
    product_side = (products, product_electrons)
    if reaction.reversible:
        sides = tuple(sorted([reactant_side, product_side]))
    else:
        sides = (reactant_side, product_side)
    collider = reaction.specific_collider.label if reaction.specific_collider is not None else None
    return (reaction.__class__, collider, pressure_dependent, bool(reaction.reversible), sides)


def answers(reactions, multiset):
    groups = defaultdict(list)
    for index, reaction in enumerate(reactions):
        groups[key(reaction, multiset)].append(index)
    flags = [False] * len(reactions)
    for members in groups.values():
        for index in members:
            flags[index] = len(members) > 1
    return flags, groups


database = KineticsDatabase()
names = []
for root, dirs, files in os.walk(LIBRARY_ROOT):
    if 'reactions.py' in files:
        names.append(os.path.relpath(root, LIBRARY_ROOT))
names.sort()
print('candidate libraries     :', len(names))

total_reactions = 0
total_changed = 0
total_merged_groups = 0
failed = []
label_collisions = 0
changed_detail = []

for name in names:
    library_name = name
    try:
        database.load_libraries(LIBRARY_ROOT, libraries=[library_name])
        library = database.libraries[library_name]
        reactions = library.get_library_reactions()
    except Exception:
        failed.append((library_name, traceback.format_exc().strip().splitlines()[-1]))
        continue

    formulas = defaultdict(set)
    for reaction in reactions:
        for spc in reaction.reactants + reaction.products:
            try:
                formulas[spc.label].add(spc.molecule[0].get_formula())
            except Exception:
                pass
    collisions = [label for label, fs in formulas.items() if len(fs) > 1]
    label_collisions += len(collisions)

    ordered_flags, ordered_groups = answers(reactions, multiset=False)
    multiset_flags, multiset_groups = answers(reactions, multiset=True)

    changed = [i for i in range(len(reactions)) if ordered_flags[i] != multiset_flags[i]]
    merged = len(ordered_groups) - len(multiset_groups)

    total_reactions += len(reactions)
    total_changed += len(changed)
    total_merged_groups += merged

    if changed or merged or collisions:
        print('  %-52s reactions=%-6d changed=%-4d groups_merged=%-4d label_collisions=%d'
              % (library_name, len(reactions), len(changed), merged, len(collisions)))
        for i in changed[:6]:
            changed_detail.append('    %s :: %s  (ordered=%s -> multiset=%s)'
                                  % (library_name, reactions[i],
                                     ordered_flags[i], multiset_flags[i]))

print('')
print('=== CENSUS ===')
print('libraries loaded              :', len(names) - len(failed), 'of', len(names))
print('reactions examined            :', total_reactions)
print('entries whose answer CHANGES  :', total_changed)
print('groups merged by the multiset :', total_merged_groups)
print('label/formula collisions      :', label_collisions)
for line in changed_detail:
    print(line)
if failed:
    print('libraries that failed to load :', len(failed))
    for library_name, why in failed:
        print('    %-40s %s' % (library_name, why[:100]))
sys.stdout.flush()
