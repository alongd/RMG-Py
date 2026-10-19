"""Pair proxies over the family-tree classes of each participant.

Two-participant compile proxies over all single structures grow
quadratically.  The rate RMG assigns a bimolecular reaction is read from the
family tree's nodes for the pair, and those nodes depend on each participant
only through what the tree's groups match on it.  So participants that every
requested family's tree reads alike are interchangeable as pair partners:
one witness per class suffices, and the pairs are formed over the class
witnesses.  Nothing here names a polymer.

``tree_fingerprint(family, molecule)``: what the family reads of one
participant in its bimolecular directions; equal fingerprints give equal
template sets against any partner:

- **Per-reactant trees** (the forward template's reactants are tree tops,
  rate rules): for every slot and every labelling RMG makes of the
  participant (``_match_reactant_to_template`` over its reactive resonance
  forms), the node ``descend_tree`` reaches from that slot's top, as
  ``get_reaction_template`` descends one reactant at a time.
- **One merged Root** (an auto-generated tree, a bimolecular template of one
  entry): ``get_reaction_template`` descends the two reactants merged into
  one structure.  A node's group then matches the pair exactly when its
  connected components embed in the two participants.  So the fingerprint
  holds, for every slot labelling and every tree node, whether each
  component, and the whole group, embeds in the participant alone.  Labels
  that occur in more than one Root component are left unpinned (RMG's
  merged descent sees them as a list); every other label is pinned.
- **A bimolecular reverse direction** that RMG generates (a family that is
  not its own reverse): the rate is the forward tree's on the products,
  which no single participant determines. A participant that fills a slot
  of that template is its own class, so its pairs stay all pairs.

A family with no bimolecular direction reads nothing (empty fingerprint).

``class_pairs(families, molecules)``: the classes (equal fingerprints under
every family) and one pair per two classes whose members can fill two
different slots of one bimolecular direction.  A participant whose
fingerprint is empty under every family reacts in no bimolecular direction
and is paired with nothing.

**Which members stand for a class pair.** RMG drops a reaction whose
products are its reactants. A radical paired with its own chain capped at the
radical site trades one hydrogen between two copies of the same species,
so that pair loses the reaction at that site while every other member pair
keeps it. A class pair is therefore stood for by its first member pair in
order whose participants are not one another with a radical capped by a
hydrogen (``_hydrogen_partners``); only when every member pair is such a
pair is it stood for by the first.  Every other member pair then reaches
only templates the chosen one reaches.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Mapping, Sequence


def _slot_templates(template, reactant_num):
    """The per-slot template reactants of a direction (split, as RMG splits
    a merged one-entry template of a bimolecular direction), and whether
    the template is one merged Root."""
    entries = template.reactants
    if len(entries) == 1 and reactant_num == 2:
        return entries[0].item.split(), True
    return list(entries), False


def _labelled(molecule, mapping):
    copy = molecule.copy(deep=True)
    index = {atom: i for i, atom in enumerate(molecule.atoms)}
    for atom in copy.atoms:
        atom.label = ""
    for atom, group_atom in mapping.items():
        copy.atoms[index[atom]].label = group_atom.label
    return copy


def _labellings(family, molecule, slot_template) -> list:
    from rmgpy.kmc import rmg_moves as rm

    out = []
    for form in rm.reactant_forms(molecule, "all"):
        for mapping in family._match_reactant_to_template(form, slot_template):
            out.append(_labelled(form, mapping))
    return out


def _embeds(group, structure, pinned) -> bool:
    """Whether ``group`` embeds in ``structure`` with its ``pinned`` labels on
    the structure's atoms of the same label (as ``match_node_to_structure``,
    the structure's labelled atoms that the group does not name set aside)."""
    centers = group.get_all_labeled_atoms()
    atoms = structure.get_all_labeled_atoms()
    initial = {}
    for label, center in centers.items():
        if label in pinned:
            if label not in atoms or isinstance(atoms[label], list) or isinstance(center, list):
                return False
            if not atoms[label].is_specific_case_of(center):
                return False
            initial[atoms[label]] = center
    flagged = [atom for label, atom in atoms.items()
               if label not in centers and not isinstance(atom, list)]
    for atom in flagged:
        atom.ignore = True
    try:
        return bool(structure.is_subgraph_isomorphic(group, initial))
    finally:
        for atom in flagged:
            atom.ignore = False


def _merged_reading(family, labelled, pinned, nodes) -> tuple:
    from rmgpy.molecule.group import Group

    reading = []
    for label in nodes:
        group = family.groups.entries[label].item
        if not isinstance(group, Group):
            raise ValueError(f"{family.label}: merged-tree node {label} is not a group")
        reading.append((tuple(_embeds(c, labelled, pinned) for c in group.split()),
                        _embeds(group, labelled, pinned)))
    return tuple(reading)


def _directions(family):
    """(name, template, reactant number, generated reverse) per direction RMG
    generates (a family that declares no reactant or product number has its
    template's)."""
    forward = family.reactant_num or len(family.forward_template.reactants)
    out = [("forward", family.forward_template, forward, False)]
    if family.reverse_template is not None and not family.own_reverse:
        reverse = family.product_num or len(family.reverse_template.reactants)
        out.append(("reverse", family.reverse_template, reverse, True))
    return out


def tree_fingerprint(family, molecule, identity=None) -> frozenset:
    """What ``family`` reads of one participant in its bimolecular directions
    (module doc): a frozenset of ``(direction, slot, reading)``.  ``identity``
    names the participant in the reading of a generated bimolecular reverse
    direction (default: its own adjacency list)."""
    out = set()
    for direction, template, number, reverse in _directions(family):
        if number != 2:
            continue
        slots, merged = _slot_templates(template, number)
        if len(slots) != 2:
            continue
        if reverse:
            name = identity if identity is not None else molecule.to_adjacency_list()
            for slot, slot_template in enumerate(slots):
                if _labellings(family, molecule, slot_template):
                    out.add((direction, slot, ("identity", name)))
            continue
        if merged:
            counts = defaultdict(int)
            for component in slots:
                for label in component.get_all_labeled_atoms():
                    counts[label] += 1
            pinned = {label for label, n in counts.items() if n == 1}
            nodes = sorted(family.groups.entries)
            for slot, component in enumerate(slots):
                for labelled in _labellings(family, molecule, component):
                    out.add((direction, slot, _merged_reading(family, labelled, pinned, nodes)))
        else:
            for slot, entry in enumerate(slots):
                for labelled in _labellings(family, molecule, entry):
                    node = family.groups.descend_tree(labelled, labelled.get_all_labeled_atoms(),
                                                      root=entry, strict=True)
                    out.add((direction, slot, None if node is None else node.label))
    return frozenset(out)


def class_pairs(families: Mapping, molecules: Mapping) -> tuple:
    """``(pairs, classes)`` over ``molecules`` (name -> molecule, in order):
    ``classes`` the participants with equal fingerprints under every family
    of ``families`` (name -> loaded family), in order; ``pairs`` one member
    pair ``(a, b)`` (a before or equal to b in order) per two classes that
    fill two different slots of one bimolecular direction, chosen as the
    module doc says."""
    names = list(molecules)
    fingerprints = {name: tuple(tree_fingerprint(families[f], molecules[name], identity=name)
                                for f in sorted(families))
                    for name in names}
    classes = {}
    for name in names:
        classes.setdefault(fingerprints[name], []).append(name)
    groups = list(classes.values())
    order = {name: k for k, name in enumerate(names)}

    def slots(name):
        return [{(direction, slot) for direction, slot, _reading in fingerprint}
                for fingerprint in fingerprints[name]]

    pairs = []
    for i, first in enumerate(groups):
        for second in groups[i:]:
            if not any(_two_slots(a, c) for a, c in zip(slots(first[0]), slots(second[0]))):
                continue
            candidates = all_pairs(first) if first is second else \
                tuple((a, c) for a in first for c in second)
            chosen = next((pair for pair in candidates
                           if not _hydrogen_partners(molecules[pair[0]], molecules[pair[1]])),
                          candidates[0])
            pairs.append(tuple(sorted(chosen, key=order.__getitem__)))
    return tuple(pairs), tuple(tuple(members) for members in groups)


def _hydrogen_partners(first, second) -> bool:
    """Whether one molecule is the other with one radical site capped by a
    hydrogen."""
    from rmgpy.molecule import Atom, Bond

    for radical, other in ((first, second), (second, first)):
        if len(radical.atoms) + 1 != len(other.atoms):
            continue
        for index, atom in enumerate(radical.atoms):
            if not atom.radical_electrons:
                continue
            capped = radical.copy(deep=True)
            site = capped.atoms[index]
            site.decrement_radical()
            hydrogen = Atom(element="H")
            capped.add_atom(hydrogen)
            capped.add_bond(Bond(site, hydrogen, order=1))
            capped.update_multiplicity()
            if capped.is_isomorphic(other, strict=False):
                return True
    return False


def _two_slots(first: set, second: set) -> bool:
    return any(d1 == d2 and s1 != s2 for d1, s1 in first for d2, s2 in second)


def all_pairs(names: Sequence[str]) -> tuple:
    """Every unordered pair of ``names``, each with itself included."""
    names = list(names)
    return tuple((a, b) for i, a in enumerate(names) for b in names[i:])
