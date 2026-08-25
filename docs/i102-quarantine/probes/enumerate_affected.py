#!/usr/bin/env python3
"""
Enumerate the quarantined kinetics entries **from the database**.

Nothing here knows which entries are affected. It loads every kinetics family in
``database.directory``, asks each one whether it carries a quarantine manifest, and
lets the manifest's own criterion pick the entries out of what is actually loaded.
Change the data and the answer changes with it; that is the point.

Usage (from a directory whose ``rmgrc`` points at the database under test)::

    python enumerate_affected.py

Exits 0 and prints the affected set. Exits 1 if no quarantine manifest is found at
all, since that means the marker is not where the gate will look for it.
"""

import os
import sys

from rmgpy import settings
from rmgpy.data.kinetics.database import KineticsDatabase


def main():
    families_path = os.path.join(settings['database.directory'], 'kinetics', 'families')
    print('database.directory : {0}'.format(settings['database.directory']))
    print('families           : {0}'.format(families_path))

    database = KineticsDatabase()
    database.load_families(path=families_path, families='all')

    quarantined = {label: family for label, family in database.families.items()
                   if family.quarantine is not None}

    print('families loaded    : {0}'.format(len(database.families)))
    print('families quarantined: {0}'.format(len(quarantined)))
    print()

    total = 0
    for label in sorted(quarantined):
        family = quarantined[label]
        q = family.quarantine
        affected = q.affected_entries(family)
        n_rules, n_training = len(affected['rules']), len(affected['training'])
        total += n_rules + n_training

        print('=' * 78)
        print('family                 : {0}'.format(label))
        print('state                  : {0}'.format(q.state))
        print('appliesToKineticsClass : {0} -> {1}'.format(q.kinetics_class_name, q.kinetics_class))
        print('reason                 : {0}'.format(q.reason))
        print('manifest               : {0}'.format(q.path))
        print('-' * 78)
        print('rate rules affected    : {0}'.format(n_rules))
        for entry in sorted(affected['rules'], key=lambda e: e.index):
            print('    [{0:>3}] {1:<50} {2}'.format(
                entry.index, entry.label, type(entry.data).__name__))
        print('training entries affected: {0}'.format(n_training))
        for entry in sorted(affected['training'], key=lambda e: e.index):
            print('    [{0:>3}] {1:<50} {2}'.format(
                entry.index, entry.label, type(entry.data).__name__))

        # The complement, so a reader can see the criterion is selective and not
        # simply matching everything in the family.
        all_rules = sum(len(entries) for entries in family.rules.entries.values())
        all_training = sum(len(d.entries) for d in family.depositories
                           if d.label.endswith('/training'))
        print('-' * 78)
        print('unaffected in this family: {0} rate rules, {1} training entries'.format(
            all_rules - n_rules, all_training - n_training))

    print('=' * 78)
    print('TOTAL AFFECTED ENTRIES : {0}'.format(total))

    if not quarantined:
        print('NO QUARANTINE MANIFEST FOUND -- the gate would find nothing to enforce.',
              file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
