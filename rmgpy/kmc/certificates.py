"""Complete rate-selection certificates: what RMG read when it chose a rate.

A certificate is built per channel member and per RMG mapping (one directed
``rmg_moves`` move: template slot order, the resonance form of every
reactant, one subgraph isomorphism per template reactant).  For each mapping
it records every tree RMG descends (``KineticsGroups.get_reaction_template``
-> ``Database.descend_tree``, rmgpy/data/base.py): the ROOT check of every top
node tried and every child predicate tried, matched or failed, with

* the predicate's semantics (``nodes``: group atoms with their atom types,
  radical electrons, lone pairs, charges and ``props`` such as ``inRing``;
  group bonds with their orders; logic nodes with their components, in RMG's
  evaluation order and with its short-circuit);
* how RMG anchors it: ``anchored`` (every group label has one witness atom:
  RMG's initial map), ``partial``, ``unanchored`` (a label is a list on either
  side, e.g. R_Recombination's two ``*``: RMG builds no initial map for it and
  the subgraph search runs over the WHOLE structure) or ``labels_missing``
  (strict matching fails before any search);
* ``candidates``: for every free (unanchored) labelled group atom, every
  witness atom of the structure RMG could map it to (``is_specific_case_of``),
  and ``ignored``: labelled witness atoms RMG hides from the search;
* the **boundary signature**: the ball of the predicate's radius around the
  anchors and candidates (every embedding lies in it: each group atom is
  within the radius of a labelled group atom, and a subgraph isomorphism
  never lengthens a path) and, in the descent's atom table, every ball atom's
  element, atom type, electrons, label, ``inRing`` and FULL bond list,
  including the bonds that leave the ball;
* the whole-structure guards ``Molecule.is_subgraph_isomorphic`` reads before
  the search (group multiplicity, radical count, element counts).

Radius unit: ordinary bonds of the witness graph (``RADIUS_UNIT``); explicit
hydrogens are graph atoms (``H_HANDLING``).

Rate source.  ``KineticsFamily.get_kinetics`` (``return_all_kinetics=False``,
the compiler's call) searches every depository first, by WHOLE-reaction
isomorphism (``Reaction.is_isomorphic``: every reactant and product species
isomorphic, multiplicity included), and stops at the first depository with a
match; only then does it estimate from the tree.  Every search is recorded
(entries searched, matches, the selected entry by RMG's rank rule).  A
depository match makes the certificate ``whole_species``: its support is the
whole reactant and product species, never a window.  For a rule-estimated
channel ``depository_exclusion`` states, per record, whether any depository
entry could match a long-chain analogue (see its docstring); unknown is
refused.  The tree estimate (rate rules, averages, auto-generated walk) is
recorded as in ``channels.selection_certificate``.

Spin and domain: every species of the channel carries its spin verdict
(``domain.spin_verdict``) and its v1 domain verdict (``domain.check``); the
multiplicity convention of the bridge is a global guard.

``recertify`` regenerates the mapped reaction on an extended witness, re-runs
the selection and compares through the atom map (``selection_view``).
"""

import hashlib
import json
from collections import defaultdict
from typing import Dict, Iterable, Optional, Sequence, Tuple

from rmgpy.kmc import channels as chm
from rmgpy.kmc import domain as dom
from rmgpy.kmc import rmg_moves as rm

CERTIFICATE_VERSION = "kmc-rate-certificate/2"
REFERENCE_SCHEMA = "kmc-rate-certificate-ref/1"
RADIUS_UNIT = ("ordinary bonds of the witness graph (bond order ignored); the radius of a "
               "group is the largest bond distance from its labelled atoms to any group atom")
H_HANDLING = ("explicit: hydrogens are graph atoms, a bond to H is one bond, an H inside a "
              "ball is listed with its own row")
REPEAT_UNIT_FORMULA = (("C", 8), ("H", 8))
WHOLE_SPECIES = "whole_species"
TREE = "tree_predicates"


class CertificateError(ValueError):
    """The certificate cannot be built (bad input)."""


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def canonical_json(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def certificate_sha256(certificate: dict) -> str:
    return hashlib.sha256(canonical_json(certificate).encode("utf-8")).hexdigest()


def _digest(value) -> str:
    return hashlib.sha256(repr(value).encode("utf-8")).hexdigest()[:16]


def aromatic_root(molecule):
    """The input written in its most aromatic resonance form (most benzene
    bonds; ties: least integer-keyed signature), atom ids kept.  RMG can return
    a Kekule input as an extra REACTIVE form (I-134), so the forms enumerated
    from a Kekule input depend on how the witness was written; the forms
    enumerated from the most aromatic form do not.  Every certificate is built
    from this representation."""
    from rmgpy.molecule.resonance import generate_resonance_structures
    if len(molecule.atoms) < 3:
        return molecule
    base = molecule.copy(deep=True)
    for atom, original in zip(base.atoms, molecule.atoms):
        atom.id = original.id
    forms = generate_resonance_structures(base, keep_isomorphic=True, filter_structures=True)
    ids = sorted(a.id for a in molecule.atoms)
    keyed = []
    for f in forms:
        if sorted(a.id for a in f.atoms) != ids:
            raise CertificateError("resonance generation lost the atom ids")
        keyed.append((-sum(1 for b in f.get_all_edges() if b.is_benzene()),
                      rm._signature(rm.molecule_rows(f)), f))
    keyed.sort(key=lambda x: (x[0], x[1]))
    return rm.molecule_from_rows(rm.molecule_rows(keyed[0][2]))


def witness_molecules(graphs: Sequence[str], offset: int = 0):
    """RMG Molecules of adjacency lists; atom ids are the adjacency positions,
    consecutive across the participants (the certificate's id namespace for a
    record)."""
    from rmgpy.molecule import Molecule
    mols = []
    for graph in graphs:
        m = Molecule().from_adjacency_list(graph)
        for atom in m.atoms:
            atom.id = offset
            offset += 1
        mols.append(m)
    return mols


def _atom_row(atom):
    atomtype = atom.atomtype.label if atom.atomtype is not None else None
    in_ring = atom.props.get("inRing") if hasattr(atom, "props") else None
    bonds = sorted([b.id, float(bond.order)] for b, bond in atom.bonds.items())
    return [atom.id, atom.element.symbol, atomtype, int(atom.radical_electrons),
            int(atom.charge), int(atom.lone_pairs), atom.label or "", in_ring, bonds]


def _group_semantics(group):
    index = {a: k for k, a in enumerate(group.atoms)}
    atoms = []
    for k, a in enumerate(group.atoms):
        atoms.append({"index": k, "label": a.label or "",
                      "atomtypes": sorted(t.label for t in a.atomtype),
                      "radical_electrons": sorted(a.radical_electrons or []),
                      "lone_pairs": sorted(a.lone_pairs or []),
                      "charge": sorted(a.charge or []),
                      "props": {str(k2): v for k2, v in sorted((a.props or {}).items())
                                if k2 != "reg_dim_r" and isinstance(v, (bool, int, str, float))}})
    bonds = []
    for a in group.atoms:
        for b, bond in a.bonds.items():
            if index[a] < index[b]:
                bonds.append([index[a], index[b], sorted(float(o) for o in bond.order)])
    return {"kind": "group", "atoms": atoms, "bonds": sorted(bonds),
            "multiplicity": sorted(group.multiplicity or []),
            "radical_count": int(getattr(group, "radicalCount", 0) or 0),
            "element_count": {str(k): int(v) for k, v in
                              sorted((getattr(group, "elementCount", {}) or {}).items())}}


def _radius(group):
    return chm._group_radius(group)


def _ball(structure, centres, radius):
    return chm._ball(structure, centres, radius)


# ---------------------------------------------------------------------------
# predicates and descents
# ---------------------------------------------------------------------------
class _Recorder:
    """Records predicates of one structure (one descent): the atom table and
    the node semantics."""

    def __init__(self, database, structure, atoms, nodes):
        self.database = database
        self.structure = structure
        self.atoms = atoms           # RMG's label dict (get_all_labeled_atoms)
        self.nodes = nodes           # shared node-semantics dict
        self.table = {}              # atom id -> row

    def _note(self, ball):
        for a in ball:
            if a.id not in self.table:
                self.table[a.id] = _atom_row(a)

    def node_semantics(self, entry):
        from rmgpy.data.base import LogicNode
        if entry.label in self.nodes:
            return
        item = entry.item
        if isinstance(item, LogicNode):
            comps = [c if isinstance(c, str) else repr(c) for c in item.components]
            self.nodes[entry.label] = {"kind": "logic", "operator": type(item).__name__,
                                       "invert": bool(item.invert), "components": comps}
        else:
            self.nodes[entry.label] = _group_semantics(item)

    def evaluate(self, node):
        """(predicate, verdict) for a tree node, by RMG's own call; logic
        nodes are expanded in RMG's order and checked against RMG's verdict."""
        from rmgpy.data.base import LogicNode
        entry = self.database.entries[node] if isinstance(node, str) else node
        verdict = bool(self.database.match_node_to_structure(entry, self.structure, self.atoms,
                                                              strict=True))
        self.node_semantics(entry)
        if isinstance(entry.item, LogicNode):
            pred = self._logic(entry.label, entry.item)
            if pred["raw"] != verdict:
                raise chm.UnsupportedSelection(
                    f"logic node {entry.label}: recorded verdict {pred['raw']} != RMG's {verdict}")
            pred.pop("raw")
            pred["node"] = entry.label
            pred["result"] = "matched" if verdict else "failed"
            return pred, verdict
        pred = self._group(entry.label, entry.item, verdict)
        return pred, verdict

    def _logic(self, label, logic):
        from rmgpy.data.base import LogicNode
        comps = []
        raw = None
        is_or = type(logic).__name__ == "LogicOr"
        for comp in logic.components:
            if isinstance(comp, LogicNode):
                sub = self._logic(repr(comp), comp)
                value = sub.pop("raw")
                sub["node"] = repr(comp)
                sub["result"] = "matched" if value else "failed"
            else:
                sub, value = self.evaluate(comp)
            comps.append(sub)
            if is_or and value:
                raw = True
                break
            if not is_or and not value:
                raw = False
                break
        if raw is None:
            raw = not is_or
        raw = raw != bool(logic.invert)
        return {"logic": type(logic).__name__, "invert": bool(logic.invert),
                "components": comps, "evaluated": len(comps),
                "of": len(logic.components), "raw": raw}

    def _group(self, label, group, verdict):
        structure, atoms = self.structure, self.atoms
        centers = group.get_all_labeled_atoms()
        anchors, free, missing = {}, [], []
        for lab, center in sorted(centers.items()):
            if lab not in atoms:
                missing.append(lab)
                continue
            atom = atoms[lab]
            if isinstance(center, list) or isinstance(atom, list):
                free.append(lab)
            else:
                anchors[lab] = atom
        ignored = sorted(a.id for lab, a in _flat_labels(structure) if lab not in centers)
        guards = {
            "multiplicity": {"group": sorted(group.multiplicity or []),
                             "structure": int(structure.multiplicity),
                             "read": bool(group.multiplicity)},
            "radical_count": {"group_minimum": int(group.radicalCount or 0),
                              "structure": int(structure.get_radical_count()),
                              "read": bool(group.radicalCount and group.radicalCount > 0)},
            "element_count": {"group_minimum": {str(k): int(v) for k, v in
                                                sorted((group.elementCount or {}).items())},
                              "structure": {str(k): int(v) for k, v in sorted(
                                  structure.get_element_count().items())
                                  if k in (group.elementCount or {})},
                              "read": bool(group.elementCount)},
        }
        out = {"node": label, "result": "matched" if verdict else "failed", "logic": None,
               "guards": guards, "ignored": ignored}
        if missing:
            out.update({"anchoring": "labels_missing", "labels_missing": missing,
                        "anchors": {}, "free_labels": [], "candidates": {}, "radius": None,
                        "ball": []})
            if verdict:
                raise chm.UnsupportedSelection(f"node {label} matched with labels missing")
            return out
        radius = _radius(group)
        if radius is None:
            raise chm.UnsupportedSelection(f"node {label} has atoms not connected to a label")
        ignored_ids = set(ignored)
        candidates = {}
        for lab in free:
            gatoms = centers[lab] if isinstance(centers[lab], list) else [centers[lab]]
            hits = sorted({a.id for a in structure.atoms if a.id not in ignored_ids
                           and any(a.is_specific_case_of(g) for g in gatoms)})
            candidates[lab] = hits
        by_id = {a.id: a for a in structure.atoms}
        centres = [a for a in anchors.values()] + [by_id[i] for ids in candidates.values()
                                                    for i in ids]
        ball = _ball(structure, centres, radius)
        self._note(ball)
        anchoring = "anchored" if not free else ("unanchored" if not anchors else "partial")
        out.update({"anchoring": anchoring,
                    "scan": "ball" if not free else
                            "whole structure (candidate atoms of the free labels)",
                    "anchors": {lab: a.id for lab, a in sorted(anchors.items())},
                    "free_labels": free, "candidates": candidates, "radius": radius,
                    "ball": sorted([a.id, d] for a, d in ball.items())})
        if verdict:
            out["embeddings"] = self._embeddings(group, anchors, centers, label)
        return out

    def _embeddings(self, group, anchors, centers, label):
        """RMG's subgraph search with RMG's initial map and ignore flags."""
        structure = self.structure
        initial = {anchors[lab]: centers[lab] for lab in anchors}
        flagged = [a for lab, a in _flat_labels(structure) if lab not in centers]
        for a in flagged:
            a.ignore = True
        try:
            maps = structure.find_subgraph_isomorphisms(group, initial_map=initial)
        finally:
            for a in flagged:
                a.ignore = False
        if not maps:
            raise chm.UnsupportedSelection(f"node {label}: RMG matched but no embedding found")
        index = {g: k for k, g in enumerate(group.atoms)}
        keyed = sorted(sorted([index[g], s.id] for s, g in m.items()) for m in maps)
        return {"count": len(keyed), "least": keyed[0],
                "images": sorted({s.id for m in maps for s in m})}


def _flat_labels(structure):
    out = []
    for lab, a in structure.get_all_labeled_atoms().items():
        for atom in (a if isinstance(a, list) else [a]):
            out.append((lab, atom))
    return out


def _descend(recorder, root):
    """``Database.descend_tree`` from ``root`` (root check included)."""
    root_pred, ok = recorder.evaluate(root)
    out = {"top": root.label, "root": root_pred, "steps": [], "leaf": None}
    if not ok:
        return out
    node = root
    while True:
        children, hits = [], []
        for child in node.children:
            pred, hit = recorder.evaluate(child)
            children.append(pred)
            if hit:
                hits.append(child)
        step = {"node": node.label, "children": children}
        out["steps"].append(step)
        if len(hits) == 1:
            step["choice"], step["rule"] = hits[0].label, "single match"
            node = hits[0]
        elif len(hits) > 1:
            step["choice"], step["rule"] = hits[0].label, f"first of {len(hits)} overlapping matches"
            node = hits[0]
        elif node.children and node.children[-1].label.startswith("Others-"):
            step["choice"], step["rule"] = node.children[-1].label, "no match: Others- fallback"
            out["leaf"] = node.children[-1].label
            return out
        else:
            step["choice"], step["rule"] = None, "no child matched: stop"
            out["leaf"] = node.label
            return out


def template_descents(family, reactant_molecules, nodes):
    """Every descent ``groups.get_reaction_template`` performs on labelled
    reactants (template order), with root checks, as recorded descents."""
    from rmgpy.data.base import LogicNode
    database = family.groups
    tops = []
    for entry in database.top:
        if entry not in tops:
            tops.append(entry)
    out = []
    if len(tops) == 1 and len(reactant_molecules) > 1:
        merged = reactant_molecules[0]
        for other in reactant_molecules[1:]:
            merged = merged.merge(other)
        atoms = merged.get_all_labeled_atoms()
        rec = _Recorder(database, merged, atoms, nodes)
        d = _descend(rec, tops[0])
        d.update({"structure": "merged reactants", "labelled": chm._labelled_ids(atoms),
                  "atoms": sorted(rec.table.values())})
        out.append(d)
        return out
    for entry in tops:
        labels = [] if isinstance(entry.item, LogicNode) else \
            list(entry.item.get_all_labeled_atoms())
        for k, mol in enumerate(reactant_molecules):
            if not all(mol.contains_labeled_atom(lab) for lab in labels):
                out.append({"top": entry.label, "structure": f"template reactant {k}",
                            "skipped": "reactant lacks the top node's labels "
                                       f"{sorted(labels)} (RMG does not try it)",
                            "root": None, "steps": [], "leaf": None, "labelled": {},
                            "atoms": []})
                continue
            atoms = mol.get_all_labeled_atoms()
            rec = _Recorder(database, mol, atoms, nodes)
            d = _descend(rec, entry)
            d.update({"structure": f"template reactant {k}",
                      "labelled": chm._labelled_ids(atoms), "atoms": sorted(rec.table.values())})
            out.append(d)
    return out


def descent_leaves(descents):
    return [d["leaf"] for d in descents if d.get("leaf") is not None]


# ---------------------------------------------------------------------------
# depository search
# ---------------------------------------------------------------------------
def _species(mol):
    from rmgpy.species import Species
    spc = Species(molecule=[mol.copy(deep=True)])
    if len(mol.atoms) > 1:
        spc.generate_resonance_structures(keep_isomorphic=False)
    return spc


def depository_search(family, reactants, products):
    """RMG's depository search (``get_kinetics_from_depository`` for every
    depository, in order, stopping at the first with a match), recorded."""
    from rmgpy.reaction import Reaction
    rxn = Reaction(reactants=[_species(m) for m in reactants],
                   products=[_species(m) for m in products])
    searches = []
    selected = None
    for depository in family.depositories:
        entries = list(depository.entries.values())
        matches = []
        for entry in entries:
            if entry.item.is_isomorphic(rxn):
                matches.append({"index": int(entry.index), "label": entry.label,
                                "rank": entry.rank,
                                "is_forward": bool(entry.item.is_isomorphic(
                                    rxn, either_direction=False))})
        searches.append({"depository": depository.label, "entries_searched": len(entries),
                         "matches": matches})
        if matches:
            ranked = matches
            if any(m["rank"] == 0 for m in ranked) and not all(m["rank"] == 0 for m in ranked):
                ranked = [m for m in ranked if m["rank"] != 0]
            ranked = sorted(ranked, key=lambda m: (m["rank"] if m["rank"] is not None else 0,
                                                   m["index"]))
            selected = dict(ranked[0], depository=depository.label)
            break
    return {"searches": searches, "selected": selected,
            "identified_by": "Reaction.is_isomorphic (either direction): every reactant and "
                             "product species isomorphic (Species.is_isomorphic over resonance "
                             "structures; multiplicity, charges and electrons compared)",
            "reads": "whole reactant and product species",
            "stops_at_first_depository_with_a_match": True}


def max_depository_heavy_atoms(family) -> int:
    """Largest heavy-atom count of any species of any depository entry."""
    best = 0
    for depository in family.depositories:
        for entry in depository.entries.values():
            for spc in list(entry.item.reactants) + list(entry.item.products):
                mol = spc.molecule[0] if hasattr(spc, "molecule") else spc
                best = max(best, sum(1 for a in mol.atoms if not a.is_hydrogen()))
    return best


def _formula(mol):
    out = defaultdict(int)
    for a in mol.atoms:
        out[a.element.symbol] += 1
    return dict(out)


def _depository_formulas(family):
    rows = []
    for depository in family.depositories:
        for entry in depository.entries.values():
            sides = []
            for side in (entry.item.reactants, entry.item.products):
                sides.append(tuple(sorted(tuple(sorted(_formula(
                    s.molecule[0] if hasattr(s, "molecule") else s).items())) for s in side)))
            rows.append((depository.label, int(entry.index), sides[0], sides[1]))
    return rows


def depository_exclusion(family, reactants, products, *, cache=None,
                         repeat=REPEAT_UNIT_FORMULA) -> dict:
    """Can any depository entry match a long-chain analogue of this reaction?

    A long-chain analogue adds ``n_i >= 0`` repeat units (``repeat``, C8H8 for
    PS) to species ``i``, at least one unit in total, and never removes atoms.
    ``Reaction.is_isomorphic`` needs every species of the entry isomorphic to
    a species of the reaction, so (a) equal heavy-atom counts and (b) equal
    formulas.  Argument, in order:

    1. ``excluded:size`` if a species of the witness has more heavy atoms than
       every species of every depository entry (``H_max``): an analogue keeps
       that species at least as large, so no entry can match it.
    2. ``excluded:formula`` if no entry has, on either orientation, the
       species formula multisets of ANY growth vector (each species grown
       independently on both sides, bounded by ``H_max``; a superset of the
       real analogues, so the exclusion is sound).
    3. Otherwise ``refused:unknown`` with the candidate entries."""
    cache = {} if cache is None else cache
    key = ("hmax", family.label)
    if key not in cache:
        cache[key] = max_depository_heavy_atoms(family)
    hmax = cache[key]
    heavy = [sum(1 for a in m.atoms if not a.is_hydrogen()) for m in list(reactants) + list(products)]
    out = {"repeat_unit": dict(repeat), "depository_max_heavy_atoms": hmax,
           "witness_max_heavy_atoms": max(heavy) if heavy else 0}
    if heavy and max(heavy) > hmax:
        out.update({"verdict": "excluded:size",
                    "argument": f"a witness species has {max(heavy)} heavy atoms > {hmax}, the "
                                "largest species of any depository entry; long-chain analogues "
                                "only add atoms, and Reaction.is_isomorphic needs every species "
                                "isomorphic"})
        return out
    fkey = ("formulas", family.label)
    if fkey not in cache:
        cache[fkey] = _depository_formulas(family)
    rows = cache[fkey]
    rep = dict(repeat)
    rep_heavy = sum(v for k, v in rep.items() if k != "H")

    def growths(mols):
        options = []
        for m in mols:
            f = _formula(m)
            h = sum(v for k, v in f.items() if k != "H")
            per = []
            n = 0
            while h + n * rep_heavy <= hmax:
                g = dict(f)
                for k, v in rep.items():
                    g[k] = g.get(k, 0) + n * v
                per.append((n, tuple(sorted(g.items()))))
                n += 1
            options.append(per)
        return options

    def side_sets(mols):
        import itertools
        sets = set()
        for combo in itertools.product(*growths(mols)):
            sets.add((sum(n for n, _f in combo) > 0, tuple(sorted(f for _n, f in combo))))
        return sets

    rs, ps = side_sets(reactants), side_sets(products)
    candidates = []
    for label, index, er, ep in rows:
        for a, b in ((er, ep), (ep, er)):
            hit_r = [g for g, f in rs if f == a]
            hit_p = [g for g, f in ps if f == b]
            if hit_r and hit_p and (any(hit_r) or any(hit_p)):
                candidates.append({"depository": label, "index": index})
                break
    if not candidates:
        out.update({"verdict": "excluded:formula",
                    "argument": "no depository entry has the species formulas of any growth of "
                                "the witness species by whole repeat units (bounded by the "
                                "largest depository species); formula equality is necessary "
                                "for Reaction.is_isomorphic"})
        return out
    out.update({"verdict": "refused:unknown", "candidates": candidates,
                "argument": "a depository entry has the formulas of a grown witness; whether "
                            "it is isomorphic to an actual long-chain analogue is not decided"})
    return out


# ---------------------------------------------------------------------------
# certifier
# ---------------------------------------------------------------------------
class Certifier:
    """Builds certificates; caches depository searches (by whole-species
    identity), descents (by labelled structure) and resonance forms."""

    def __init__(self, families, resonance: str = "all", max_channel_sets: int = 64):
        self.families = families
        self.resonance = resonance
        self.max_channel_sets = max_channel_sets
        self._channel_keys = []
        self._records = {}
        self.cache = {}
        self._descents = {}
        self._depository = {}

    # -- forms and labelled structures --------------------------------------
    def channel_set(self, family, reactants):
        """Forward channels on the aromatic roots of ``reactants`` (cached by
        the witness's integer-keyed signature; bounded)."""
        roots = [aromatic_root(m) for m in reactants]
        key = ("channels", family.label, tuple(rm._signature(rm.molecule_rows(m))
                                                for m in roots))
        if key not in self.cache:
            if len(self._channel_keys) >= self.max_channel_sets:
                self.cache.pop(self._channel_keys.pop(0), None)
            self.cache[key] = chm.build_channels(family, roots, "forward",
                                                 resonance=self.resonance)
            self._channel_keys.append(key)
        return self.cache[key]

    def _forms(self, channel_set, k):
        key = ("forms", id(channel_set), k)
        if key not in self.cache:
            self.cache[key] = (channel_set, rm.reactant_forms(channel_set.reactants[k],
                                                               channel_set.resonance))
        return self.cache[key][1]

    def _labelled(self, channel_set, move):
        label_of = {a: lab for _s, lab, a in move.labels}
        out = []
        for k in range(len(channel_set.reactants)):
            src = self._forms(channel_set, k)[move.forms[k]]
            form = src.copy(deep=True)
            for atom, original in zip(form.atoms, src.atoms):
                atom.id = original.id
                atom.label = label_of.get(original.id, "")
            out.append(form)
        return out

    def _descent_key(self, family, molecules):
        return (family.label, tuple(sorted(
            (rm._signature(rm.molecule_rows(m)),
             tuple(sorted((a.id, a.label) for a in m.atoms if a.label))) for m in molecules)))

    def descents(self, family, channel_set, move, nodes, store):
        """Refs of the recorded descents of one mapping (template direction)."""
        mols = self._labelled(channel_set, move)
        key = self._descent_key(family, mols)
        if key not in self._descents:
            self._descents[key] = template_descents(family, mols, nodes)
        ref = _digest(key)
        store[ref] = self._descents[key]
        return ref

    # -- channel ------------------------------------------------------------
    def certify(self, family, reactants, *, pathway=None, products=None, member_labels=None,
                channel_set=None, species=None) -> dict:
        """Certificate of the forward (family-template direction) channel(s) on
        ``reactants`` selected by ``pathway`` + product ``products`` (species
        keys), or by ``member_labels`` (a set of frozensets of (label, atom id):
        every channel with a member one of whose mappings has those labels)."""
        species = species or chm.SpeciesKeys()
        if channel_set is None:
            channel_set = self.channel_set(family, reactants)
        chosen = []
        for ch in channel_set.channels:
            if pathway is not None and ch.pathway != tuple(pathway):
                continue
            if products is not None and ch.products != tuple(sorted(products)):
                continue
            if member_labels is not None:
                if not any(_label_set(m) in member_labels
                           for m in self._member_moves(channel_set, ch, all_members=True)):
                    continue
            chosen.append(ch)
        cert = {"version": CERTIFICATE_VERSION, "family": family.label,
                "radius_unit": RADIUS_UNIT, "h_handling": H_HANDLING,
                "multiplicity_convention": dom.MULTIPLICITY_CONVENTION,
                "witness": {"molecules": [sorted(a.id for a in m.atoms)
                                          for m in channel_set.reactants]},
                "channels": [], "descents": {}, "nodes": {}, "species": {}, "gaps": []}
        if not chosen:
            cert["gaps"].append("no channel on the witness matches the selector")
        for ch in chosen:
            cert["channels"].append(self._channel(family, channel_set, ch, cert, species))
        cert["global_guards"] = _global_guards(cert)
        cert["gaps"].extend(completeness_gaps(cert))
        cert["complete"] = not cert["gaps"]
        # completeness is about what RMG read; whether a long-chain analogue
        # could switch the rate to a depository entry is a separate verdict
        cert["depository_exclusion_refused"] = [
            ";".join(ch["pathway"]) for ch in cert["channels"]
            if not ch.get("whole_species")
            and (ch.get("depository_exclusion") or {}).get("verdict", "").startswith("refused")]
        return cert

    def certify_channel(self, family, channel_set, channel, species=None) -> dict:
        return self.certify(family, channel_set.reactants, pathway=channel.pathway,
                            products=channel.products, channel_set=channel_set, species=species)

    def _member_moves(self, channel_set, channel, all_members=False):
        keys = rm.product_key
        members = set(channel.members)
        return [m for m in channel_set.moves
                if m.status in chm.ALLOWED_STATUSES and m.template is not None
                and tuple(m.template) == channel.pathway and m.products is not None
                and keys(m) in members]

    def _channel(self, family, channel_set, ch, cert, species):
        keys = rm.product_key
        out = {"direction": ch.direction, "pathway": list(ch.pathway), "kind": ch.kind,
               "reactants": list(ch.reactants), "products": list(ch.products),
               "size": ch.size, "mappings": list(ch.mappings), "support": ch.support,
               "ambiguous_with": [list(p) for p in ch.ambiguous_with]}
        if ch.direction != "forward":
            out["gap"] = "certificates are built in the family template (forward) direction"
            return out
        by_outcome = defaultdict(list)
        for m in self._member_moves(channel_set, ch):
            by_outcome[keys(m)].append(m)
        others = channel_set.by_outcome()
        members = []
        for outcome in ch.members:
            moves = sorted(by_outcome.get(outcome, ()), key=rm.Move.sort_key)
            maps = []
            for m in moves:
                try:
                    ref = self.descents(family, channel_set, m, cert["nodes"], cert["descents"])
                except chm.UnsupportedSelection as exc:
                    maps.append({"order": list(m.order), "forms": list(m.forms),
                                 "labels": [list(x) for x in m.labels], "status": m.status,
                                 "descent": None, "gap": str(exc)})
                    continue
                maps.append({"order": list(m.order), "forms": list(m.forms),
                             "labels": [list(x) for x in m.labels], "status": m.status,
                             "descent": ref,
                             "leaves": descent_leaves(cert["descents"][ref])})
            amb = sorted({c.pathway for c in others.get(outcome, ()) if c.pathway != ch.pathway})
            members.append({"outcome": _digest(outcome), "mappings": maps,
                            "ambiguous_with": [list(p) for p in amb]})
        out["members"] = members
        # species verdicts
        rep = ch.representatives[0]
        mols = list(channel_set.reactants) + [rm.molecule_from_rows(p) for p in
                                              rm.split_rows(rep.product_rows())]
        for mol in mols:
            k = species(mol)
            if k not in cert["species"]:
                spin = dom.spin_verdict(mol)
                cert["species"][k] = {
                    "spin": {"status": spin.status, "reason": spin.reason,
                             "radical_electrons": spin.radical_electrons,
                             "convention": spin.convention},
                    "domain": dom.check(mol).as_dict()}
        # rate source
        reactants = list(channel_set.reactants)
        products = [rm.molecule_from_rows(p) for p in rm.split_rows(rep.product_rows())]
        dkey = (family.label, tuple(sorted(ch.reactants)), tuple(sorted(ch.products)))
        if dkey not in self._depository:
            self._depository[dkey] = depository_search(family, reactants, products)
        out["depository"] = self._depository[dkey]
        selected = out["depository"]["selected"]
        if selected is not None:
            out["source_kind"] = "training_exact_match"
            out["whole_species"] = True
            out["support_scope"] = WHOLE_SPECIES
            out["selection"] = {"training_entry": selected}
        else:
            out["whole_species"] = False
            out["support_scope"] = TREE
            out["depository_exclusion"] = depository_exclusion(family, reactants, products,
                                                                cache=self.cache)
            try:
                kind, sel = self._tree_selection(family, ch)
                out["source_kind"], out["selection"] = kind, sel
            except chm.UnsupportedSelection as exc:
                out["source_kind"], out["selection"] = "unsupported", {"reason": str(exc)}
        out["rmg_check"] = self._rmg_check(family, channel_set, ch, out)
        return out

    def _tree_selection(self, family, ch):
        template = family.retrieve_template(list(ch.pathway))
        if family.auto_generated:
            walk = chm.autogenerated_walk(family, ch.pathway)   # before RMG: RMG may not terminate
            kinetics, _entry = family.get_kinetics_for_template(template, degeneracy=1)
            sel = chm._autogenerated_selection(family, ch.pathway, kinetics, walk=walk)
            return "autogenerated_tree", sel
        kinetics, _entry = family.get_kinetics_for_template(template, degeneracy=1)
        kind, sel = chm._rate_rule_selection(family, ch.pathway, kinetics)
        sel["supplying_node"] = sel["supplying_templates"][0] \
            if len(sel["supplying_templates"]) == 1 else None
        return kind, sel

    def _rmg_check(self, family, channel_set, ch, out):
        """RMG's own get_kinetics on the representative must agree with the
        recorded depository result (run only when the walk terminates)."""
        if out.get("source_kind") == "unsupported":
            return {"skipped": "selection unsupported"}
        rxn, _r, _p = chm.template_reaction(family, channel_set, ch.representatives[0])
        rxn.reactants = [_species(s.molecule[0]) for s in rxn.reactants]
        rxn.products = [_species(s.molecule[0]) for s in rxn.products]
        _k, source, entry, _f = family.get_kinetics(rxn, template_labels=list(ch.pathway),
                                                     degeneracy=1, return_all_kinetics=False)
        rmg_dep = None if isinstance(source, str) else int(entry.index)
        ours = out["depository"]["selected"]
        agree = (rmg_dep is None and ours is None) or (ours is not None and rmg_dep == ours["index"])
        if not agree:
            out.setdefault("gaps", []).append(f"RMG selected depository entry {rmg_dep}, the "
                                              f"recorded search {ours}")
        return {"rmg_depository_entry": rmg_dep, "agrees": agree}

    # -- records ------------------------------------------------------------
    def certify_record(self, record) -> dict:
        record = _as_dict(record)
        rs = record.get("rate_source") or {}
        kind = rs.get("kind")
        head = {"event_id": record.get("event_id"), "family": record.get("family"),
                "template": record.get("template"), "rate_source_kind": kind,
                "orientation": record.get("orientation")}
        if kind not in ("RMG family estimate", "reference-thermo reverse"):
            return {**head, "certificate": None, "source_kind": "not_rmg_selection",
                    "complete": False, "gaps": [f"rate source {kind!r} is not an RMG family "
                                                "selection (library or archived rate)"]}
        family = self.families[record["family"]]
        side, mols, others = self.template_side(record)
        pathway = tuple(record["template"].split(";")) if record["template"] else None
        rkey = (record["family"], record["template"],
                tuple(record[side]), tuple(_graph_field(record, side, other=True)))
        if rkey not in self._records:
            species = chm.SpeciesKeys()
            products = tuple(sorted(species(m) for m in others))
            self._records[rkey] = self.certify(family, mols, pathway=pathway,
                                               products=products, species=species)
        cert = self._records[rkey]
        padded = (record.get("proxy_padding") or {}).get("status") == "padded"
        return {**head, "template_side": side, "padded": padded,
                "role": "rate_source" if kind == "RMG family estimate" else "reverse_of_rate_source",
                "certificate": cert, "complete": cert["complete"], "gaps": list(cert["gaps"])}

    def template_side(self, record):
        """(field name, family-template reactant molecules, other side): the
        as_generated record's reactants, the reversed record's products (the
        rate witness graphs when padded)."""
        record = _as_dict(record)
        reactant_field = "rate_witness_reactant_graphs" if record.get(
            "rate_witness_reactant_graphs") else "reactant_graphs"
        product_field = "rate_witness_product_graphs" if record.get(
            "rate_witness_product_graphs") else "product_graphs"
        if record.get("orientation", "as_generated") == "as_generated":
            a, b = reactant_field, product_field
        else:
            a, b = product_field, reactant_field
        mols = witness_molecules(record[a])
        others = witness_molecules(record[b], offset=sum(len(m.atoms) for m in mols))
        return a, mols, others

    def recertify(self, record, witness, atom_map, certificate=None) -> dict:
        """Re-run the selection of a record on an extended witness.

        ``witness``: molecules of the record's family-template reactants on
        the extension (for a reversed record: the extension of its products,
        the side the rate is selected on), with atom ids; ``atom_map``: record
        certificate witness id -> extended witness id.  See
        ``recertify_certificate``."""
        record = _as_dict(record)
        old = certificate or self.certify_record(record)
        old_cert = old.get("certificate") if "certificate" in old else old
        if old_cert is None:
            return {"certificate": None, "same_selection": False, "flips": [],
                    "reason": "the record has no RMG selection to recertify"}
        return self.recertify_certificate(self.families[record["family"]], old_cert,
                                          witness, atom_map)

    def recertify_certificate(self, family, old_cert, witness, atom_map) -> dict:
        """Regenerate the mapped reaction on ``witness``, re-run the selection
        and compare through the map.

        The new certificate covers every channel of the extension holding a
        mapping whose labelled atoms are the transported labelled atoms of a
        certified mapping.  ``same_selection``: for every certified mapping,
        the extension has a mapping with the transported labels and an equal
        ``mapping_view`` (pathway, channel kind and support, source kind,
        supplying node(s), contributors, depository entry, and the tree
        descents node by node with each verdict).  Atom ids and whole-species
        keys are never compared; members the extension adds (not images of
        certified mappings) are new events and are not compared."""
        transport = {int(k): int(v) for k, v in atom_map.items()}
        wanted = set()
        unmapped = set()
        for ch in old_cert["channels"]:
            for mem in ch.get("members", []):
                for mp in mem["mappings"]:
                    ids = [a for _s, _l, a in mp["labels"]]
                    unmapped.update(a for a in ids if a not in transport)
                    wanted.add(frozenset((lab, transport.get(a)) for _s, lab, a in mp["labels"]))
        if unmapped:
            return {"certificate": None, "same_selection": False, "flips": [],
                    "reason": f"reacting atoms {sorted(unmapped)} are not in the atom map"}
        new = self.certify(family, witness, member_labels=wanted)
        old_view = mapping_view(old_cert, transport)
        new_view = mapping_view(new)
        missing = sorted(str(sorted(k)) for k in old_view if k not in new_view)
        changed = sorted(str(sorted(k)) for k in old_view
                         if k in new_view and new_view[k] != old_view[k])
        return {"certificate": new, "same_selection": not missing and not changed,
                "missing_mappings": missing, "changed_mappings": changed,
                "flips": flipped_predicates(old_cert, new, transport)}

    def certify_reverse_channel(self, family, channel_set, channel):
        """Certificates of a reverse-direction channel: RMG selects its rate in
        the family template direction, on the move's products.  One
        certificate per member: the forward channel rebuilt on that member's
        products (atom ids kept) with the channel's pathway and the reverse
        channel's reactant species as products.  Returns [(certificate,
        template-side molecules)]."""
        out = []
        for move in channel.representatives:
            mols = [rm.molecule_from_rows(p) for p in rm.split_rows(move.product_rows())]
            species = chm.SpeciesKeys()
            cert = self.certify(family, mols, pathway=channel.pathway,
                                products=channel.reactants, species=species)
            out.append((cert, mols))
        return out


def _graph_field(record, side, other=False):
    names = {"reactant_graphs": "product_graphs", "product_graphs": "reactant_graphs",
             "rate_witness_reactant_graphs": "rate_witness_product_graphs",
             "rate_witness_product_graphs": "rate_witness_reactant_graphs"}
    return record[names[side]] if other else record[side]


def _label_set(move):
    return frozenset((lab, a) for _s, lab, a in move.labels)


def _as_dict(record):
    if isinstance(record, dict):
        return record
    return record.to_dict()


# ---------------------------------------------------------------------------
# completeness, guards, views
# ---------------------------------------------------------------------------
def _iter_predicates(pred):
    yield pred
    for c in pred.get("components", []) or []:
        yield from _iter_predicates(c)


def iter_descent_predicates(descent):
    if descent.get("root"):
        yield from _iter_predicates(descent["root"])
    for step in descent.get("steps", []):
        for c in step["children"]:
            yield from _iter_predicates(c)


def completeness_gaps(cert) -> list:
    """Every reason the certificate is incomplete (empty: complete)."""
    gaps = []
    for ch in cert["channels"]:
        name = ";".join(ch["pathway"])
        if ch.get("gap"):
            gaps.append(f"{name}: {ch['gap']}")
            continue
        if len(ch.get("members", [])) != ch["size"]:
            gaps.append(f"{name}: {len(ch.get('members', []))} members certified of {ch['size']}")
        for k, mem in enumerate(ch.get("members", [])):
            if not mem["mappings"]:
                gaps.append(f"{name}: member {k} has no mapping")
            if len(mem["mappings"]) != ch["mappings"][k]:
                gaps.append(f"{name}: member {k} has {len(mem['mappings'])} mappings certified "
                            f"of {ch['mappings'][k]}")
            for mp in mem["mappings"]:
                if mp.get("descent") is None:
                    gaps.append(f"{name}: member {k}: {mp.get('gap')}")
                    continue
                descents = cert["descents"][mp["descent"]]
                if mp["leaves"] != list(ch["pathway"]):
                    gaps.append(f"{name}: member {k}: recorded leaves {mp['leaves']}")
                for d in descents:
                    if d.get("skipped"):
                        continue
                    if d.get("root") is None:
                        gaps.append(f"{name}: descent from {d['top']} without root check")
                    table = {row[0] for row in d["atoms"]}
                    for p in iter_descent_predicates(d):
                        if p.get("logic"):
                            continue
                        no_candidates = p.get("anchoring") == "unanchored" and not any(
                            p.get("candidates", {}).values())
                        if p.get("anchoring") != "labels_missing" and not no_candidates and (
                                p.get("radius") is None or not p.get("ball")):
                            gaps.append(f"{name}: predicate {p['node']} without a ball")
                        if any(i not in table for i, _d in p.get("ball", [])):
                            gaps.append(f"{name}: predicate {p['node']} ball atom without a row")
        if ch.get("depository") is None:
            gaps.append(f"{name}: depository search not recorded")
        if ch.get("source_kind") in (None, "unsupported"):
            gaps.append(f"{name}: selection unsupported: "
                        f"{(ch.get('selection') or {}).get('reason')}")
        gaps.extend(f"{name}: {g}" for g in ch.get("gaps", []))
    return gaps


def _global_guards(cert):
    read = defaultdict(set)
    for ref, descents in cert["descents"].items():
        for d in descents:
            for p in iter_descent_predicates(d):
                for name, g in (p.get("guards") or {}).items():
                    if g.get("read"):
                        read[name].add(p["node"])
    guards = [{"guard": "multiplicity_convention", "value": dom.MULTIPLICITY_CONVENTION,
               "read_by": "the bridge (rmg_moves.molecule_from_rows) and every whole-species "
                          "comparison (depository search, product identity)"},
              {"guard": "depository_whole_species",
               "read_by": "KineticsFamily.get_kinetics before any tree estimate"},
              {"guard": "inRing", "read_by": "atom predicates with an inRing prop (ring "
                                             "perception over the whole molecule; the value is "
                                             "recorded per ball atom)"}]
    for name in ("multiplicity", "radical_count", "element_count"):
        guards.append({"guard": f"group_{name}",
                       "read_by": sorted(read.get(name, ())),
                       "where": "Molecule.is_subgraph_isomorphic, whole structure"})
    return guards


def _descent_signature(descents, transport):
    out = []
    for d in descents:
        if d.get("skipped"):
            out.append((d["top"], d["structure"], "skipped"))
            continue
        steps = tuple((s["node"], tuple((c["node"], c["result"]) for c in s["children"]),
                       s["choice"], s["rule"]) for s in d["steps"])
        root = d["root"]["result"] if d.get("root") else None
        out.append((d["top"], d["structure"], root, steps, d["leaf"]))
    return tuple(out)


def _channel_selection(ch):
    sel = ch.get("selection") or {}
    dep = ch.get("depository") or {}
    return (tuple(ch["pathway"]), ch["kind"], ch.get("source_kind"), ch.get("support"),
            canonical_json(sel.get("supplying_node")),
            canonical_json([c.get("template") for c in sel.get("contributors", [])]),
            canonical_json([w.get("moved") for w in sel.get("parent_walk", [])]),
            (dep.get("selected") or {}).get("index"))


def mapping_view(cert, transport: Optional[Dict[int, int]] = None):
    """{transported mapping labels: (channel selection, descent signature)}
    over every certified mapping (``transport``: old id -> new id)."""
    t = (lambda a: a) if transport is None else (lambda a: transport.get(a, ("unmapped", a)))
    out = defaultdict(set)
    for ch in cert["channels"]:
        sel = _channel_selection(ch)
        for mem in ch.get("members", []):
            for mp in mem["mappings"]:
                labels = frozenset((lab, t(a)) for _s, lab, a in mp["labels"])
                sig = _descent_signature(cert["descents"][mp["descent"]], transport) \
                    if mp.get("descent") else None
                out[labels].add((sel, sig))
    return {k: frozenset(v) for k, v in out.items()}


def selection_view(cert, transport: Optional[Dict[int, int]] = None):
    """The whole certificate's selection, ids transported (``transport``) or
    raw: the set of (mapping labels, channel selection, descent signature)."""
    return frozenset(mapping_view(cert, transport).items())


def flipped_predicates(old, new, transport):
    """(transported labels, node, old result, new result) for every certified
    mapping of ``old`` (labels transported) whose predicate verdict differs
    on the corresponding mapping of ``new``, or that was tried in only one."""
    def index(cert, tr):
        out = defaultdict(dict)
        for ch in cert["channels"]:
            for mem in ch.get("members", []):
                for mp in mem["mappings"]:
                    if not mp.get("descent"):
                        continue
                    labels = frozenset((lab, tr(a)) for _s, lab, a in mp["labels"])
                    for d in cert["descents"][mp["descent"]]:
                        for p in iter_descent_predicates(d):
                            out[labels][(d["top"], p["node"])] = p["result"]
        return out
    a = index(old, lambda x: transport.get(x, ("unmapped", x)))
    b = index(new, lambda x: x)
    flips = []
    for labels in sorted(a, key=lambda s: sorted(map(str, s))):
        pa, pb = a.get(labels, {}), b.get(labels, {})
        for key in sorted(set(pa) | set(pb)):
            ra, rb = pa.get(key, "not tried"), pb.get(key, "not tried")
            if ra != rb:
                flips.append({"labels": sorted([lab, str(x)] for lab, x in labels),
                              "top": key[0], "node": key[1], "old": ra, "new": rb})
    return flips


# ---------------------------------------------------------------------------
# record binding
# ---------------------------------------------------------------------------
def certificate_reference(record_certificate: dict) -> dict:
    """The serialized reference an EventRecord carries (``rate_certificate``):
    the sha256 of the record certificate (canonical JSON), where it lives,
    and a summary."""
    cert = record_certificate.get("certificate")
    kinds = sorted({ch.get("source_kind") for ch in (cert or {}).get("channels", [])} - {None})
    return {"schema": REFERENCE_SCHEMA,
            "sha256": certificate_sha256(record_certificate),
            "store": "artifact.rate_certificates",
            "template_side": record_certificate.get("template_side"),
            "padded": bool(record_certificate.get("padded")),
            "complete": bool(record_certificate.get("complete")),
            "source_kinds": kinds or [record_certificate.get("source_kind")],
            "whole_species": any(ch.get("whole_species") for ch in (cert or {}).get("channels", []))}


def attach_rate_certificates(records, families, certifier: Optional[Certifier] = None):
    """Certify every record from its actual rate source (the rate witness
    graphs when padded) and attach the reference.  Returns (records, store):
    EventRecords with ``rate_certificate`` set (event ids unchanged: the
    reference is outside the event hash) and ``{sha256: record certificate}``."""
    from dataclasses import replace
    certifier = certifier or Certifier(families)
    store = {}
    out = []
    for record in records:
        try:
            rc = certifier.certify_record(record)
        except Exception as exc:  # reported, never silently dropped
            rc = {"event_id": _as_dict(record).get("event_id"), "certificate": None,
                  "complete": False, "gaps": [f"{type(exc).__name__}: {exc}"],
                  "source_kind": "error"}
        ref = certificate_reference(rc)
        store[ref["sha256"]] = rc
        out.append(replace(record, rate_certificate=ref) if not isinstance(record, dict)
                   else dict(record, rate_certificate=ref))
    return out, store


# ---------------------------------------------------------------------------
# +1 repeat-unit extension of a witness (PS)
# ---------------------------------------------------------------------------
def extend_witness(molecules, protected: Iterable[int], *, shift: int = 100000):
    """Add one PS repeat unit (C8H8) to every chain molecule (one with a
    phenyl) at its saturated chain end farthest from ``protected`` atom ids:
    a CH3 head becomes CH3-CH(Ph)-CH2-, a PhCH2 tail becomes -CH(Ph)-CH2-CH2Ph
    (head-to-tail).  One H of the end carbon is removed.  Every kept atom is
    renumbered ``id + shift``; new atoms get fresh ids.  Returns (molecules,
    atom_map old id -> new id, info per molecule)."""
    protected = set(protected)
    next_id = shift + 1 + max(a.id for m in molecules for a in m.atoms)
    new_mols, atom_map, info = [], {}, []
    for mol in molecules:
        rows = rm.molecule_rows(mol)
        has_ring = any(mol.is_atom_in_cycle(a) for a in mol.atoms)
        new_rows = {}
        if not has_ring:
            for i, row in rows.items():
                new_rows[i + shift] = row[:5] + (tuple((j + shift, o) for j, o in row[5]),)
                atom_map[i] = i + shift
            new_mols.append(rm.molecule_from_rows(new_rows))
            info.append({"extended": False, "reason": "no phenyl (not a chain)"})
            continue
        by_id = {a.id: a for a in mol.atoms}
        ring = {a.id for a in mol.atoms if mol.is_atom_in_cycle(a)}
        ends = []
        for a in mol.atoms:
            if not a.is_carbon() or a.id in ring or a.radical_electrons:
                continue
            heavy = [b for b in a.bonds if not b.is_hydrogen()]
            hs = [b for b in a.bonds if b.is_hydrogen()]
            if any(bond.order != 1.0 for bond in a.bonds.values()) or not hs:
                continue
            backbone = [b for b in heavy if b.id not in ring]
            phenyl = [b for b in heavy if b.id in ring]
            if len(backbone) == 1 and not phenyl and len(hs) == 3:
                ends.append((a, "head"))
            elif len(backbone) == 1 and len(phenyl) == 1 and len(hs) == 2:
                ends.append((a, "tail"))
        if not ends:
            for i, row in rows.items():
                new_rows[i + shift] = row[:5] + (tuple((j + shift, o) for j, o in row[5]),)
                atom_map[i] = i + shift
            new_mols.append(rm.molecule_from_rows(new_rows))
            info.append({"extended": False, "reason": "no saturated chain end"})
            continue
        dist = chm._distances(mol, sorted(protected & set(by_id))) if protected & set(by_id) \
            else {}
        # the end farthest from the protected (reacting) atoms; prefer an end
        # with an H that is not protected (removing a reacting H unmaps it)
        end, kind = max(ends, key=lambda e: (
            any(b.is_hydrogen() and b.id not in protected for b in e[0].bonds),
            dist.get(e[0].id, 10 ** 6), -e[0].id))
        hs = sorted(b.id for b in end.bonds if b.is_hydrogen())
        free = [h for h in hs if h not in protected]
        h_removed = (free or hs)[0]
        edges = {i: dict(row[5]) for i, row in rows.items()}
        del edges[h_removed]
        del edges[end.id][h_removed]
        base = {i: rows[i][:5] for i in rows if i != h_removed}

        def add(element, nbrs):
            nonlocal next_id
            i = next_id
            next_id += 1
            base[i] = (element, None, 0, 0, 0)
            edges[i] = {}
            for j, o in nbrs:
                edges[i][j] = o
                edges[j][i] = o
            return i

        def phenyl_on(c):
            ring_ids = [add("C", [(c, 1.0)])]
            for _ in range(5):
                ring_ids.append(add("C", [(ring_ids[-1], 1.5)]))
            edges[ring_ids[0]][ring_ids[-1]] = 1.5
            edges[ring_ids[-1]][ring_ids[0]] = 1.5
            for r in ring_ids[1:]:
                add("H", [(r, 1.0)])

        if kind == "head":   # CH3-  ->  CH3-CH(Ph)-CH2-
            ca = add("C", [(end.id, 1.0)])
            add("H", [(ca, 1.0)])
            phenyl_on(ca)
            cb = add("C", [(ca, 1.0)])
            for _ in range(3):
                add("H", [(cb, 1.0)])
        else:                # PhCH2-  ->  -CH(Ph)-CH2-CH2Ph
            ca = add("C", [(end.id, 1.0)])
            for _ in range(2):
                add("H", [(ca, 1.0)])
            cb = add("C", [(ca, 1.0)])
            for _ in range(2):
                add("H", [(cb, 1.0)])
            phenyl_on(cb)
        final = {}
        for i, row in base.items():
            ni = i + shift if i in rows else i
            final[ni] = row + (tuple(sorted(((j + shift if j in rows else j), o)
                                            for j, o in edges[i].items())),)
            if i in rows:
                atom_map[i] = ni
        new_mols.append(rm.molecule_from_rows(final))
        info.append({"extended": True, "end": kind, "end_atom": end.id,
                     "distance_to_protected": dist.get(end.id), "removed_h": h_removed})
    return new_mols, atom_map, info


# module-level conveniences ---------------------------------------------------
def certify_channel(family, channel_set, channel, certifier: Optional[Certifier] = None) -> dict:
    certifier = certifier or Certifier({family.label: family}, channel_set.resonance)
    return certifier.certify_channel(family, channel_set, channel)


def certify_record(record, families, certifier: Optional[Certifier] = None) -> dict:
    return (certifier or Certifier(families)).certify_record(record)


def recertify(record, witness, atom_map, *, families, certificate=None,
              certifier: Optional[Certifier] = None) -> dict:
    return (certifier or Certifier(families)).recertify(record, witness, atom_map, certificate)
