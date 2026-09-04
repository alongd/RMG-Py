#!/usr/bin/env python
"""I-067 round 2: recompute the recipe-action census against the PRODUCTION family database.

Round 1 reported "8 actions" from a grep over families/*/groups.py. The review said the
executor supports NINE symbols. Both can be true -- supported vocabulary vs. vocabulary
actually used -- and this script computes each separately instead of asserting either.

It parses every recipe literally with ast, so an action symbol that a regex would miss
(different quoting, line-wrapped list, nested in a dict) still gets counted, and any symbol
outside the executor's known set is reported as UNKNOWN rather than silently dropped.
"""
import ast
import os
import re
import sys
from collections import Counter

DB = '/home/alon/Code/RMG-database/input/kinetics/families'

# The symbols ReactionRecipe._apply actually dispatches on (family.py:526,536,549).
EXECUTOR_SYMBOLS = {
    'CHANGE_BOND', 'FORM_BOND', 'BREAK_BOND',
    'GAIN_RADICAL', 'LOSE_RADICAL',
    'GAIN_CHARGE', 'LOSE_CHARGE',
    'GAIN_PAIR', 'LOSE_PAIR',
}

# Actions that would change the element composition of the product set. The executor has
# no branch for any of these; the census exists to prove none appears in the database.
ATOM_CHANGING = {'ADD_ATOM', 'REMOVE_ATOM', 'ADD_GROUP', 'REMOVE_GROUP'}


def recipe_actions(path):
    """Yield every action symbol in every ReactionRecipe literal in `path`."""
    src = open(path, encoding='utf-8', errors='replace').read()
    for m in re.finditer(r'\brecipe\s*\(\s*actions\s*=\s*(\[)', src):
        start = m.start(1)
        depth = 0
        for i in range(start, len(src)):
            if src[i] == '[':
                depth += 1
            elif src[i] == ']':
                depth -= 1
                if depth == 0:
                    end = i + 1
                    break
        else:
            print('UNTERMINATED recipe literal in %s' % path, file=sys.stderr)
            continue
        try:
            actions = ast.literal_eval(src[start:end])
        except Exception as exc:
            print('UNPARSED recipe in %s: %s' % (path, exc), file=sys.stderr)
            continue
        yield [a for a in actions if isinstance(a, (list, tuple)) and a]


PER_ATOM = {'GAIN_RADICAL', 'LOSE_RADICAL', 'GAIN_CHARGE', 'LOSE_CHARGE', 'GAIN_PAIR', 'LOSE_PAIR'}
BOND = {'CHANGE_BOND', 'FORM_BOND', 'BREAK_BOND'}


def nonunique_label_families(recipes):
    """Families where a per-atom action targets a label the recipe itself uses twice.

    ReactionRecipe._apply (family.py:531-536) loops `for atom in
    struct.get_labeled_atoms(label)` and applies the delta to EVERY match, so the net
    radical/charge/pair change is `change * multiplicity`. A recipe that names the same
    label as both centers of a bond action proves multiplicity >= 2 for that label,
    without needing the template group. This is the counter-example to "the net delta is
    a per-family constant."
    """
    out = {}
    for fam, actions in recipes.items():
        twice = set()
        for a in actions:
            if a[0] in BOND and len(a) >= 4 and a[1] == a[3]:
                twice.add(a[1])
        hits = [a for a in actions if a[0] in PER_ATOM and len(a) >= 2 and a[1] in twice]
        if hits:
            out[fam] = hits
    return out


def main():
    counts = Counter()
    recipes = {}
    files = 0
    for fam in sorted(os.listdir(DB)):
        gp = os.path.join(DB, fam, 'groups.py')
        if not os.path.isfile(gp):
            continue
        files += 1
        actions = [a for rec in recipe_actions(gp) for a in rec]
        if actions:
            recipes[fam] = actions
        counts.update(str(a[0]).strip() for a in actions)

    print('families scanned:            %d' % files)
    print('families with a recipe:      %d' % len(recipes))
    print('recipe actions parsed:       %d' % sum(counts.values()))
    if not sum(counts.values()):
        print('\nFAIL: parsed zero actions -- the extractor matched nothing, so every '
              'downstream count below is vacuous.', file=sys.stderr)
        return 2
    print('distinct symbols USED in DB: %d' % len(counts))
    print('symbols the EXECUTOR knows:  %d' % len(EXECUTOR_SYMBOLS))
    print()
    for sym, n in counts.most_common():
        tag = '' if sym in EXECUTOR_SYMBOLS else '   <-- UNKNOWN TO EXECUTOR'
        print('  %-14s %6d%s' % (sym, n, tag))

    unused = sorted(EXECUTOR_SYMBOLS - set(counts))
    print()
    print('executor symbols NOT used by any production family: %s' % (unused or 'none'))

    unknown = sorted(set(counts) - EXECUTOR_SYMBOLS)
    atomch = sorted(set(counts) & ATOM_CHANGING)
    print('symbols used but not dispatched by the executor:    %s' % (unknown or 'none'))
    print('ATOM-CHANGING actions found:                        %s' % (atomch or 'none'))

    nonuniq = nonunique_label_families(recipes)
    print()
    print('--- radical/charge/pair delta: is it a per-family CONSTANT? ---')
    print('families where a per-atom action targets a provably non-unique label: %d of %d'
          % (len(nonuniq), len(recipes)))
    for fam in sorted(nonuniq)[:12]:
        print('  %-38s %s' % (fam, nonuniq[fam]))
    if len(nonuniq) > 12:
        print('  ... and %d more' % (len(nonuniq) - 12))

    ok = not unknown and not atomch
    print()
    print('VERDICT (composition): %s' % ('CONSERVED by every production recipe -- no action '
                                         'adds or removes an atom' if ok else 'CLAIM BROKEN'))
    print('VERDICT (radical delta): %s'
          % ('NOT a per-family constant -- %d families apply a per-atom action to a label '
             'that matches multiple atoms, so the delta is change x multiplicity, known '
             'only at match time' % len(nonuniq) if nonuniq
             else 'no counter-example found in the database'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
